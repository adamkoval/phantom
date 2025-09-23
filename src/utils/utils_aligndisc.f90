!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module misalignedutils
!
! Routines to manage misaligned and off-center discs
!  Can handle gas disc, gas disc + sinks and warped discs
!
! :References:
!
! :Owner: Adam Koval
!
! :Runtime parameters: None
!
! :Dependencies: centreofmass, part, physcon, vectorutils
!
 use centreofmass, only:reset_centreofmass, get_total_angular_momentum
 
 implicit none

 public :: rotate_align, rotate_coordinates

 private

contains

subroutine rotate_align(xyzh, vxyzu, npart, sphere_radius)
 use part,    only:nptmass, xyzmh_ptmass, vxyz_ptmass, isdead_or_accreted

 real,             intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: sphere_radius
 integer,          intent(in) :: npart

 integer :: i
 integer :: imax_mass
 real, dimension(3) :: disc_center

 real, dimension(3) :: L_tot, L_tot_local
 integer :: npartlocal, ilocal
 real, allocatable :: xyzh_local(:,:), vxyzu_local(:,:)
 real :: rotate_about_z, rotate_about_y

 ! Center disc and get angular momentum vector
 if (nptmass > 0) then
    ! Find most massive sink
    imax_mass = maxloc(xyzmh_ptmass(4,1:nptmass), 1)
    disc_center = xyzmh_ptmass(1:3,imax_mass)
    do i=1,npart
      xyzh(1:3,i) = xyzh(1:3,i) - disc_center
    enddo
    do i=1,nptmass
      xyzmh_ptmass(1:3,i) = xyzmh_ptmass(1:3,i) - disc_center
    enddo
    print*, 'Centering disc on most massive sink, ID ', imax_mass, ' at ', disc_center

    ! Find local L around centered sink
    npartlocal = 0
    do i=1, npart
       if (.not.isdead_or_accreted(xyzh(4,i))) then
          if (sqrt(xyzh(1,i)**2 + xyzh(2,i)**2 + xyzh(3,i)**2) <= sphere_radius) then
            npartlocal = npartlocal + 1
          endif
       endif
    enddo

    print*, 'Found', npartlocal, ' active particles within ', sphere_radius, 'AU of sink'

    ! Allocate local arrays
    allocate(xyzh_local(4,npartlocal))
    allocate(vxyzu_local(4,npartlocal))

    ! Populate local arrays
    ilocal=0
    do i=1, npart
      if (.not.isdead_or_accreted(xyzh(4,i))) then
          if (sqrt(xyzh(1,i)**2 + xyzh(2,i)**2 + xyzh(3,i)**2) <= sphere_radius) then
            ilocal = ilocal + 1
            xyzh_local(:,ilocal) = xyzh(:,i)
            vxyzu_local(:,ilocal) = vxyzu(:,i)
          endif
       endif
    enddo

    call get_total_angular_momentum(xyzh_local,vxyzu_local,npartlocal,L_tot_local,&
                                 xyzmh_ptmass,vxyz_ptmass,nptmass)
    L_tot = L_tot_local
   !  print*, 'Disc centered on sink'
 else
    call reset_centreofmass(npart,xyzh,vxyzu)
    call get_total_angular_momentum(xyzh,vxyzu,npart,L_tot)
    print*, 'Disc centered on gas COM'
 endif

 call rotate_coordinates(npart,xyzh,vxyzu,L_tot,rotate_about_z,rotate_about_y)

end subroutine rotate_align

!-------------------------------------------
!+
! Rotates particle coordinates to align disc plane with xy-plane
!+
!-------------------------------------------
subroutine rotate_coordinates(npart, xyzh, vxyzu, L_tot, rotate_about_z, rotate_about_y)
 use vectorutils, only:rotatevec
 use physcon, only:pi

 integer, intent(in) :: npart
 real, intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real, intent(in) :: L_tot(3)
 real, intent(out) :: rotate_about_z, rotate_about_y

 integer :: i
 real, dimension(3) :: temp, pos_vec, vel_vec
 real :: temp_mag, L_tot_mag, L_tot_rotated(3)

 ! Rotate so that L_tot is along z-axis
 temp = (/L_tot(1),L_tot(2),0./)
 temp_mag = sqrt(dot_product(temp,temp))

 if (temp_mag > tiny(temp_mag) .and. abs(temp(2)) > tiny(temp(2))) then
   rotate_about_z = -acos(dot_product((/1.,0.,0./),temp/temp_mag))*temp(2)/abs(temp(2))
 else
   rotate_about_z = 0.
 endif

 ! Now rotate about y-axis to get L_tot along z-axis
 L_tot_rotated = L_tot
 call rotatevec(L_tot_rotated,(/0.,0.,1./),rotate_about_z)
 L_tot_mag = sqrt(dot_product(L_tot_rotated,L_tot_rotated))
 if (L_tot_mag > tiny(L_tot_mag)) then
   rotate_about_y = -acos(dot_product((/0.,0.,1./),L_tot_rotated/L_tot_mag))
 else
   rotate_about_y = 0.
 endif

 print*, 'Original disc angular momentum vector L_tot:', L_tot
 print*, 'After z-rotation L_tot_rotated:', L_tot_rotated
 print*, 'Rotation angles - about z:', rotate_about_z*180./pi, ' about y:', rotate_about_y*180./pi

 ! Rotate all particle positions and velocities
 do i=1, npart
   ! Rotate positions
   pos_vec = xyzh(1:3,i)
   call rotatevec(pos_vec,(/0.,0.,1./),rotate_about_z)
   call rotatevec(pos_vec,(/0.,1.,0./),rotate_about_y)
   xyzh(1:3,i) = pos_vec

   ! Rotate velocities
   vel_vec = vxyzu(1:3,i)
   call rotatevec(vel_vec,(/0.,0.,1./),rotate_about_z)
   call rotatevec(vel_vec,(/0.,1.,0./),rotate_about_y)
   vxyzu(1:3,i) = vel_vec
 enddo

 ! Verify rotation
 call get_total_angular_momentum(xyzh,vxyzu,npart,L_tot_rotated)
 print*, 'Final L_tot after rotation:', L_tot_rotated
 print*, 'Should be close to (0,0,|L|):', sqrt(dot_product(L_tot_rotated,L_tot_rotated))

end subroutine rotate_coordinates

end module misalignedutils