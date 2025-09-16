!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module analysis
!
! Analysis routine for discs by DF, adapted from a routine by CJN.
! Edited for use with variable gammai and mui and more accurate alpha_ss calc by AKY
! Further edited to properly treat misaligned discs by AK
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - nbins : *Number of radial bins*
!   - rin   : *Inner Disc Radius*
!   - rout  : *Outer Disc Radius*
!
! :Dependencies: dim, eos, eos_stamatellos, getneighbours, infile_utils,
!   io, kernel, part, physcon, prompting, units, vectorutils, centreofmass
!
 use getneighbours,    only:generate_neighbour_lists, read_neighbours, write_neighbours, &
                           neighcount,neighb,neighmax
 use eos_stamatellos,  only:du_store
 use centreofmass,     only:reset_centreofmass,get_total_angular_momentum
 implicit none
 character(len=20), parameter, public :: analysistype = 'disc_stresses_misaligned'
 public :: do_analysis, radial_binning, calc_gravitational_forces

 ! Variables for radial binning
 integer :: nbins
 real    :: rin,rout,dr,hlim
 integer, allocatable,dimension(:)   :: ipartbin
 real,    allocatable,dimension(:)   :: rad,ninbin,sigma,csbin,vrbin,vphibin,omega
 real,    allocatable,dimension(:)   :: H,toomre_q,epicyc,part_scaleheight,tcool,h_smooth
 real,    allocatable,dimension(:)   :: alpha_reyn,alpha_grav,alpha_mag,alpha_art
 real,    allocatable,dimension(:)   :: rpart,phipart,vrpart,vphipart, gr,gphi,Br,Bphi
 real,    allocatable,dimension(:,:) :: gravxyz,zsetgas
 real,    allocatable,dimension(:)   :: mean_z,mean_vz,rms_z,rms_vz

 logical :: write_neighbour_list = .true.  ! Write the neighbour list to file, if true


 private

contains

subroutine do_analysis(dumpfile,numfile,xyzh,vxyzu,pmass,npart,time,iunit)
 use io,      only:fatal
 use part,    only:gravity,mhd,eos_vars,nptmass,xyzmh_ptmass,vxyz_ptmass,isdead_or_accreted
 use eos,     only:ieos
 use eos_stamatellos, only:eos_file,read_optab,optable
 use readwrite_dumps, only:write_smalldump

 character(len=*), intent(in) :: dumpfile
 real,             intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: pmass,time
 integer,          intent(in) :: npart,iunit,numfile

 character(len=9) :: output
 integer          :: ierr

 integer :: i
 integer :: imax_mass
 real, dimension(3) :: disc_center

 real, dimension(3) :: L_tot, L_tot_local
 integer :: npartlocal, ilocal
 real, allocatable :: xyzh_local(:,:), vxyzu_local(:,:)
 real, parameter :: sphere_radius = 100.0
 real :: rotate_about_z, rotate_about_y

 ! Code calculates the following alphas:
 ! Reynolds stress: 2*dvr*dvphi/3*cs^2
 ! Gravitational stress: gr*gphi/(4*pi*cs^2*Grho)
 ! Maxwell stress: 2*-Br*Bphi/3(4*pi*rho*cs^2)
 ! Numerical stress: 
! Print the analysis being done
 write(*,'("Performing analysis type ",A)') analysistype
 write(*,'("Input file name is ",A)') dumpfile

 write(output,"(a4,i5.5)") 'stresses_',numfile
 write(*,'("Output file name is ",A)') output

! Read analysis options
 call read_analysis_options
 print *, "ieos=", ieos
 if (ieos==24 .and. .not. allocated(optable)) call read_optab(eos_file,ierr)

 if (mhd) print*, 'This is an MHD dump: will calculate Maxwell Stress'

 ! Center disc and get angular momentum vector
 if (nptmass > 0) then
   !  call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)
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

 if (gravity) then
! Calculate gravitational forces for all particles (gradient of the potential)
    call calc_gravitational_forces(dumpfile,npart,xyzh,vxyzu)
 else
    allocate(gravxyz(3,npart))
    gravxyz(:,:) = 0.0
 endif

! Calculate the radial and tangential velocities, forces and fields
 call transform_to_cylindrical(npart,xyzh,vxyzu)

! Bin particles by radius
 call radial_binning(npart,xyzh,vxyzu,pmass,eos_vars)

! Calculate stresses
 call calc_stresses(npart,xyzh,vxyzu,pmass)

! Write out data to file
 call write_radial_data(iunit,output,time,L_tot,rotate_about_z,rotate_about_y)

! Write small dump
 print*, 'Writing small dump of rotated coordinates to ', trim(dumpfile)//'_rot'
 call write_smalldump(time,dumpfile//'_rot')

! End of analysis
 call deallocate_arrays

 print '(a,a)', 'Analysis complete for dump ',dumpfile

end subroutine do_analysis


!-------------------------------------------
!+
! Rotates particle coordinates to align disc plane with xy-plane
!+
!-------------------------------------------
subroutine rotate_coordinates(npart,xyzh,vxyzu,L_tot,rotate_about_z,rotate_about_y)
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


!-------------------------------------------
!+
! Read options for analysis from file
!+
!-------------------------------------------
subroutine read_analysis_options
 use prompting,    only:prompt
 use infile_utils, only:open_db_from_file,inopts,read_inopt,write_inopt,close_db
 use io,           only:fatal
 use eos,          only:ieos
 type(inopts), allocatable :: db(:)
 integer :: ierr,nerr,ieos_analysis
 logical :: inputexist
 character(len=21) :: inputfile
 integer, parameter :: iunit = 10

! Check for existence of input file
 inputfile = trim(analysistype)//'.options'
 inquire(file=inputfile, exist=inputexist)

 if (inputexist) then
    nerr = 0
    print '(a,a,a)', "Parameter file ",inputfile, " found: reading analysis options"
    call open_db_from_file(db,inputfile,iunit,ierr)
    call read_inopt(ieos_analysis,'ieos',db,errcount=nerr)
    call read_inopt(nbins,'nbins',db,errcount=nerr)
    call read_inopt(rin,'rin',db,errcount=nerr)
    call read_inopt(rout,'rout',db,errcount=nerr)
    call read_inopt(hlim,'hlim',db,errcount=nerr)
    call close_db(db)
    if (nerr > 0) then
       call fatal(trim(analysistype),'Error in reading '//trim(inputfile))
    endif

 else

    print '(a,a,a)', "Parameter file ",inputfile, " NOT found"
    nbins = 128; rin = 1.; rout = 100.
    call prompt('Enter ieos to use for analysis (e.g. 24): ', ieos_analysis)
    call prompt('Enter the number of radial bins: ', nbins)
    call prompt('Enter the disc inner radius: ', rin)
    call prompt('Enter the disc outer radius: ', rout)
    call prompt('Enter scale height limit for tcool averages: ', hlim) !from midplane

! Write choices to new inputfile

    open(unit=iunit,file=inputfile,status='new',form='formatted')
    write(iunit,"(a)") '# parameter options for analysis of '//trim(analysistype)
    call write_inopt(ieos_analysis,'ieos','ieos for analysis',iunit)
    call write_inopt(nbins,'nbins','Number of radial bins',iunit)
    call write_inopt(rin,'rin','Inner Disc Radius',iunit)
    call write_inopt(rout,'rout','Outer Disc Radius',iunit)
    call write_inopt(hlim,'hlim','Scale height limit for tcool averages',iunit)
    close(iunit)
 endif

 if (ieos /= ieos_analysis) then
    print *, "Changing ieos from", ieos, "to", ieos_analysis
    ieos = ieos_analysis
 endif
 print*, 'Inner Disc Radius (code units): ', rin
 print*, 'Outer Disc Radius (code units): ', rout
 print*, 'Number of bins: ', nbins
 print*, 'Using particles of h <', hlim 
 print*, 'Using ieos=',ieos
 
end subroutine read_analysis_options


!---------------------------------------------------
!+
! Computes gravitational force vector from potential gradients
! Requires neighbours to compute SPH derivatives
!+
!---------------------------------------------------
subroutine calc_gravitational_forces(dumpfile,npart,xyzh,vxyzu)

 use dim, only:gravity,maxp
 use part, only:poten,igas,iphase,maxphase,rhoh,massoftype,iamgas
 use kernel, only: get_kernel,get_kernel_grav1,cnormk

 character(len=*),intent(in) :: dumpfile
 real,intent(in) :: xyzh(:,:),vxyzu(:,:)
 integer,intent(in) :: npart


 integer :: j,k,igrav,ipart
 real,dimension(3) :: dr
 real :: rij,rij2, hj1,hj21,hj41,q2i,qi
 real :: rhoi, rhoj, wabi, grkerni, dphidhi, grpmrho1

 logical :: iamgasi,existneigh
 character(100) :: neighbourfile

 if (allocated(gravxyz))deallocate(gravxyz)
 allocate(gravxyz(3,npart))

 gravxyz(:,:) = 0.0

! Construct neighbour lists for derivative calculations

 neighbourfile = 'neigh_'//trim(dumpfile)
 inquire(file=neighbourfile,exist = existneigh)

 if (existneigh.eqv..true.) then
    print*, 'Neighbour file ', trim(neighbourfile), ' found'
    call read_neighbours(neighbourfile,npart)

 else

    ! If there is no neighbour file, generate the list
    print*, 'No neighbour file found: generating'

    call generate_neighbour_lists(xyzh,vxyzu,npart,dumpfile,write_neighbour_list)

 endif

 print*, 'Neighbour lists acquired'
 print*, 'Calculating gravitational force from potential derivative'

 do ipart=1,npart

    rhoi = rhoh(xyzh(4,ipart),massoftype(igas))

    ! Smoothing lengths
    hj1 = 1.0/xyzh(4,ipart)
    hj21 = hj1*hj1
    hj41 = hj21*hj21

    over_neighbours: do k = 1, neighcount(ipart)

       if (k>neighmax) exit over_neighbours
       j = neighb(ipart,k)

       if (maxphase==maxp) then
          iamgasi = iamgas(iphase(j))
       else
          iamgasi = .true.
       endif

       if (ipart==j) cycle
       if (.not.iamgasi) cycle

       ! Calculate gradient of SPH kernel
       ! Separation of particles

       rij2 = 0.0
       do igrav = 1,3
          dr(igrav) = xyzh(igrav,ipart) - xyzh(igrav,j)
          rij2 = rij2 + dr(igrav)*dr(igrav)
       enddo

       rij = sqrt(rij2)
       dr(:) = dr(:)/rij



       ! r/h
       q2i = rij2*hj21
       qi = rij*hj1

       !--kernel and gradient
       if (gravity) then
          call get_kernel_grav1(q2i,qi,wabi,grkerni,dphidhi)
       else
          call get_kernel(q2i,qi,wabi,grkerni)
       endif

       ! Prefactor to SPH sum = grad*mass/rho
       rhoj = rhoh(xyzh(4,j),massoftype(igas))

       grpmrho1 = grkerni*massoftype(igas)/rhoi

       ! Now calculate each gravitational force vector element
       ! F = - grad phi

       grav_force: do igrav = 1,3
          gravxyz(igrav,ipart) = gravxyz(igrav,ipart) - grpmrho1*hj41*dr(igrav)*(poten(j)-poten(ipart))
       enddo grav_force

       ! End of matrix loop

    enddo over_neighbours
    ! End of loop over nearest neighbours

    gravxyz(:,ipart) = gravxyz(:,ipart)/(massoftype(igas))

 enddo
! End loop over all particles

end subroutine calc_gravitational_forces


!---------------------------------------------------
!+
! Calculates cylindrical transformation for positions, forces and fields
!+
!---------------------------------------------------
subroutine transform_to_cylindrical(npart,xyzh,vxyzu)

 use part, only: mhd,gravity,Bxyz
 implicit none

 integer, intent(in) :: npart
 real,intent(in) ::xyzh(:,:),vxyzu(:,:)
 integer :: ipart

 allocate(rpart(npart))
 allocate(phipart(npart))
 allocate(vrpart(npart))
 allocate(vphipart(npart))
 allocate(gr(npart))
 allocate(gphi(npart))
 allocate(Br(npart))
 allocate(Bphi(npart))

 rpart(:) = 0.0
 phipart(:) = 0.0
 vrpart(:) = 0.0
 vphipart(:) = 0.0
 gr(:) = 0.0
 gphi(:) = 0.0
 Br(:) = 0.0
 Bphi(:) = 0.0

 do ipart=1,npart

    ! Calculate positions first
    rpart(ipart) = sqrt(xyzh(1,ipart)*xyzh(1,ipart) + &
        xyzh(2,ipart)*xyzh(2,ipart))

    phipart(ipart) = atan2(xyzh(2,ipart),xyzh(1,ipart))

    ! Now velocities
    vrpart(ipart) = vxyzu(1,ipart)*cos(phipart(ipart)) + &
        vxyzu(2,ipart)*sin(phipart(ipart))

    vphipart(ipart) = vxyzu(2,ipart)*cos(phipart(ipart)) - &
        vxyzu(1,ipart)*sin(phipart(ipart))

    ! Now gravitational forces
    if (gravity) then
       gr(ipart) = gravxyz(1,ipart)*cos(phipart(ipart)) + &
           gravxyz(2,ipart)*sin(phipart(ipart))

       gphi(ipart) = gravxyz(2,ipart)*cos(phipart(ipart)) - &
           gravxyz(1,ipart)*sin(phipart(ipart))
    endif

    ! Finally, B-Fields
    if (mhd) then
       Br(ipart) = Bxyz(1,ipart)*cos(phipart(ipart)) + &
           Bxyz(2,ipart)*sin(phipart(ipart))
       Bphi(ipart) = Bxyz(2,ipart)*cos(phipart(ipart)) - &
           Bxyz(1,ipart)*sin(phipart(ipart))
    endif

 enddo

end subroutine transform_to_cylindrical

!---------------------------------------------------------------
!+
! Bins particles by radius (and computes some binned quantities)
!+
!---------------------------------------------------------------

subroutine radial_binning(npart,xyzh,vxyzu,pmass,eos_vars)
 use physcon, only:pi
 use eos,     only:get_spsound,ieos
 use part,    only:rhoh,isdead_or_accreted
 use units,   only:utime

 integer,intent(in) :: npart
 real,intent(in) :: pmass
 real,intent(in) :: xyzh(:,:),vxyzu(:,:),eos_vars(:,:)

 integer :: ibin,ipart,nbinned,iallocerr
 real :: area,csi
 logical :: do_tcool

 print '(a,I4)', 'Carrying out radial binning, number of bins: ',nbins
 do_tcool = .false.
 if (allocated(du_store)) then
    do_tcool = .true.
    print *, "Calculating t_cool"
 endif
 if (do_tcool) print *, "size=", size(du_store)
 allocate(ipartbin(npart))
 allocate(rad(nbins))
 allocate(ninbin(nbins))
 allocate(sigma(nbins))
 allocate(csbin(nbins))
 allocate(omega(nbins))
 allocate(vrbin(nbins))
 allocate(vphibin(nbins))
 allocate(part_scaleheight(nbins))
 allocate(tcool(nbins))
 allocate(h_smooth(nbins))
 allocate(mean_z(nbins))
 allocate(rms_z(nbins))
 allocate(mean_vz(nbins))
 allocate(rms_vz(nbins))

 ipartbin(:) = 0
 ninbin(:) = 0.0
 sigma(:) = 0.0
 csbin(:) = 0.0
 omega(:) = 0.0
 vrbin(:) = 0.0
 vphibin(:) = 0.0
 part_scaleheight(:) = 0.0
 tcool(:) = 0.0
 h_smooth(:) = 0.
 mean_z(:) = 0.0
 rms_z(:) = 0.0
 mean_vz(:) = 0.0
 rms_vz(:) = 0.0

 allocate(zsetgas(npart,nbins),stat=iallocerr)
 ! If you don't have enough memory to allocate zsetgas, then calculate H the slow way with less memory.
 if (iallocerr/=0) then
    write(*,'(/,a)') ' WARNING: Could not allocate memory for array zsetgas!'
    write(*,'(a)')   '          (It possibly requires too much memory)'
    write(*,'(a,/)') '          Try calculate scaleheight the slow way.'
 endif

! Set up radial bins

 dr = (rout-rin)/real(nbins-1)
 do ibin=1,nbins
    rad(ibin)=rin + real(ibin-0.5)*dr
 enddo

! Now bin all particles

 nbinned = 0
 do ipart=1,npart

    ! i refers to particle, ii refers to bin
    if (.not.isdead_or_accreted(xyzh(4,ipart))) then ! IF ACTIVE

       ibin = int((rpart(ipart)-rad(1))/dr + 1)

       ! Check particle in the binning range - if not skip it
       if (ibin > nbins) ibin=0
       if (ibin < 1)  ibin=0

       if (ibin<=0) cycle

       nbinned = nbinned +1

       ninbin(ibin) = ninbin(ibin) +1
       ipartbin(ipart) = ibin
       
       csi = get_spsound(ieos,xyzh(1:3,ipart),rhoh(xyzh(4,ipart),pmass),vxyzu(:,ipart))
       csbin(ibin) = csbin(ibin) + csi
       h_smooth(ibin) = h_smooth(ibin) + xyzh(4,ipart) 
       area = pi*((rad(ibin)+0.5*dr)**2-(rad(ibin)- 0.5*dr)**2)
       sigma(ibin) = sigma(ibin) + pmass/area

       mean_z(ibin) = mean_z(ibin) + xyzh(3,ipart)
       mean_vz(ibin) = mean_vz(ibin) + vxyzu(3,ipart)

       vrbin(ibin) = vrbin(ibin) + vrpart(ipart)
       vphibin(ibin) = vphibin(ibin) + vphipart(ipart)
       omega(ibin) = omega(ibin) + vphipart(ipart)/rad(ibin)
       zsetgas(int(ninbin(ibin)),ibin) = xyzh(3,ipart)
       if (do_tcool) then
          if (abs(xyzh(3,ipart)) < hlim*csi/(vphipart(ipart)/rad(ibin)) .and. abs(du_store(ipart)) > 0.) then
             ! include only particles < hlim
             tcool(ibin) = tcool(ibin) + vxyzu(4,ipart)/du_store(ipart)
          endif
       endif
    endif

 enddo

 call calculate_H(nbins,part_scaleheight,zsetgas,int(ninbin))
 part_scaleheight(:) = part_scaleheight(:)
 print*, nbinned, ' particles have been binned'

 where(ninbin(:)/=0)
    csbin(:) = csbin(:)/ninbin(:)
    vrbin(:) = vrbin(:)/ninbin(:)
    vphibin(:) = vphibin(:)/ninbin(:)
    omega(:) = omega(:)/ninbin(:)
    h_smooth(:) = h_smooth(:)/ninbin(:)
    mean_z(:) = mean_z(:)/ninbin(:)
    mean_vz(:) = mean_vz(:)/ninbin(:)
 end where
 if (do_tcool) then
     where(ninbin(:)/=0)
       tcool(:) = tcool(:)/ninbin(:)
       tcool(:) = abs(tcool(:))*utime
    end where
 endif

 ! Calculate rms values
 do ipart=1,npart
   if (.not.isdead_or_accreted(xyzh(4,ipart))) then
      ibin = ipartbin(ipart)
      if (ibin>0 .and. ibin<=nbins) then
         rms_z(ibin) = rms_z(ibin) + (xyzh(3,ipart)-mean_z(ibin))**2
         rms_vz(ibin) = rms_vz(ibin) + (vxyzu(3,ipart)-mean_vz(ibin))**2
      endif
   endif
 enddo

 ! Normalise and bin rms values
 where(ninbin(:)/=0)
    rms_z(:) = sqrt(rms_z(:)/(ninbin(:)-1))
    rms_vz(:) = sqrt(rms_vz(:)/(ninbin(:)-1))
 elsewhere
    rms_z(:) = 0.0
    rms_vz(:) = 0.0
 end where

 print*, 'Binning Complete'

end subroutine radial_binning


!--------------------------------------------------------------
!+
! Calculates the disc stresses (alphas) and other radial quantities
!+
!--------------------------------------------------------------
subroutine calc_stresses(npart,xyzh,vxyzu,pmass)
 use physcon, only: pi,gg,kb_on_mh
 use units,   only: print_units, umass,udist,utime,unit_velocity,unit_density,unit_Bfield
 use dim,     only: gravity
 use part,    only: mhd,rhoh,alphaind,imu,itemp
! use eos,     only: ieos

 implicit none

 integer, intent(in) :: npart
 real,intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,intent(in) :: pmass
 
 integer :: ibin,ipart
 real :: cs2, dvr,dvphi,Keplog,rhopart,unit_force,beta_sph,Hbin
 beta_sph = 2.
 print*, 'Calculating Disc Stresses'

 allocate(alpha_reyn(nbins))
 allocate(alpha_grav(nbins))
 allocate(alpha_mag(nbins))
 allocate(alpha_art(nbins))
 allocate(toomre_q(nbins))
 allocate(epicyc(nbins))
 allocate(H(nbins))

 alpha_reyn(:) = 0.0
 alpha_grav(:) = 0.0
 alpha_mag(:) = 0.0
 alpha_art(:) = 0.0
 toomre_q(:) = 0.0
 epicyc(:) = 0.0

! Convert data into cgs
 print*, 'Converting data to cgs units'

 call print_units

 sigma(:) = sigma(:)*umass/(udist*udist)
 !if (ieos /= 24) then
    csbin(:) = csbin(:)*unit_velocity
 !endif

 omega(:) = omega(:)/utime
!  h_smooth(:) = h_smooth(:)*udist

 Keplog = 1.5
 unit_force = (unit_velocity/utime)

 if (gravity) then
    gr(:) = gr(:)*unit_force
    gphi(:) = gphi(:)*unit_force
 endif

 if (mhd) then
    Br(:) = Br(:)*unit_Bfield
    Bphi(:) = Bphi(:)*unit_Bfield
 endif

 do ipart=1,npart
    ibin = ipartbin(ipart)

    if (ibin<=0) cycle

    dvr = (vrpart(ipart) - vrbin(ibin))*unit_velocity
    dvphi = (vphipart(ipart) -vphibin(ibin))*unit_velocity
    rhopart = rhoh(xyzh(4,ipart),pmass)*unit_density

    alpha_reyn(ibin) = alpha_reyn(ibin) + dvr*dvphi

!  To do: Handle constant alpha_sph
    Hbin = part_scaleheight(ibin)*udist
    alpha_art(ibin) = alpha_art(ibin) + 0.06*alphaind(1,ipart)*h_smooth(ibin)*udist/Hbin &
         + (beta_sph*0.04*h_smooth(ibin)*h_smooth(ibin)*udist*udist/Hbin/Hbin)

    if (gravity) alpha_grav(ibin) = alpha_grav(ibin) + gr(ipart)*gphi(ipart)/rhopart
    if (mhd) alpha_mag(ibin) = alpha_mag(ibin) + Br(ipart)*Bphi(ipart)/rhopart

 enddo

! Normalise stresses (by number of particles in each bin)
! Also calculate other properties that rely on binned data
! (e.g. Toomre Q)

 do ibin=1,nbins
    cs2 = csbin(ibin)*csbin(ibin)

    ! calculate epicyclic frequency
    if (ibin>1) then
       epicyc(ibin) = (rad(ibin)*rad(ibin)*omega(ibin)-rad(ibin-1)*rad(ibin-1)*omega(ibin-1))/(rad(ibin)-rad(ibin-1))
    else
       epicyc(ibin) = rad(ibin)*omega(ibin)
    endif

    if (epicyc(ibin)*omega(ibin)>=0.0) then
       epicyc(ibin) = sqrt(2.0*epicyc(ibin)*omega(ibin)/rad(ibin))
    else
       epicyc(ibin) = -sqrt(2.0*abs(epicyc(ibin)*omega(ibin))/rad(ibin))
    endif

    if (sigma(ibin)>0.0) then
       toomre_q(ibin) = csbin(ibin)*epicyc(ibin)/(pi*gg*sigma(ibin))
    else
       toomre_q(ibin) = 1.0e30
    endif

    if (abs(omega(ibin))>0.0) then
       H(ibin) = csbin(ibin)/omega(ibin)
    else
       H(ibin) = 1.0e30
    endif

    if (ninbin(ibin) >0) then
       alpha_reyn(ibin) = alpha_reyn(ibin)/(Keplog*cs2*ninbin(ibin))
       alpha_art(ibin) = alpha_art(ibin)/ninbin(ibin)
       
       if (gravity) alpha_grav(ibin) = alpha_grav(ibin)/(Keplog*4.0*pi*gg*ninbin(ibin)*cs2)
       if (mhd) alpha_mag(ibin) = alpha_mag(ibin)/(Keplog*cs2*ninbin(ibin))
    else
       alpha_reyn(ibin) = 0.0
       alpha_art(ibin) = 0.0
       if (gravity) alpha_grav(ibin) = 0.0
       if (mhd) alpha_mag(ibin) = 0.0
    endif
 enddo

! Convert H back to code units

 H(:) = H(:)/udist

 print*, 'Stresses calculated'
end subroutine calc_stresses

!--------------------------------------------------------------
!+
! Writes radially binned data to file
!+
!--------------------------------------------------------------
subroutine write_radial_data(iunit,output,time,L_tot,rotate_about_z,rotate_about_y)
 use physcon, only: pi

 implicit none

 integer, intent(in) :: iunit
 real, intent(in) :: time
 real, intent(in) :: L_tot(3)
 real, intent(in) :: rotate_about_z, rotate_about_y
 character(len=*) :: output
 integer :: ibin

 print '(a,a)', 'Writing to file ',output
 open(iunit,file=output)
 write(iunit,'("# Disc Stress data at t = ",es20.12)') time
 write(iunit,'("# Disc angular momentum vector: ",3(es20.12,1x))') L_tot(1),L_tot(2),L_tot(3)
 write(iunit,'("# Rotation angles - about z (deg): ",es20.12)') rotate_about_z*180./pi
 write(iunit,'("# Rotation angles - about y (deg): ",es20.12)') rotate_about_y*180./pi
 write(iunit,"('#',18(1x,'[',i2.2,1x,a11,']',2x))") &
       1,'radius (AU)', &
       2,'sigma (cgs)', &
       3,'cs (cgs)', &
       4,'omega (cgs)', &
       5,'epicyc (cgs)',&
       6,'H (AU)',&
       7,'Toomre Q',&
       8,'alpha_reyn',&
       9,'alpha_grav',&
       10,'alpha_mag',&
       11,'alpha_art',&
       12,'particle H (au)',&
       13,'t_cool',&
       14,'<h>',&
       15,'mean_z',&
       16,'rms_z',&
       17,'mean_vz',&
       18,'rms_vz'

 do ibin=1,nbins
    write(iunit,'(18(es18.10,1X))') rad(ibin),sigma(ibin),csbin(ibin), &
            omega(ibin),epicyc(ibin),H(ibin), abs(toomre_q(ibin)),alpha_reyn(ibin), &
            alpha_grav(ibin),alpha_mag(ibin),alpha_art(ibin),part_scaleheight(ibin),&
            tcool(ibin),h_smooth(ibin),mean_z(ibin),rms_z(ibin),mean_vz(ibin),rms_vz(ibin)
 enddo

 close(iunit)

 print '(a)', 'File Write complete'

end subroutine write_radial_data

subroutine calculate_H(nbin,H,zsetgas,ninbin)
! copied from utils disc
 integer, intent(in)  :: nbin
 real,    intent(out) :: H(:)
 real,    intent(in)  :: zsetgas(:,:)
 integer, intent(in)  :: ninbin(:)
 integer :: ii
 real    :: meanzii

 do ii = 1,nbin
    if (ninbin(ii)==0) then
       meanzii = 0.
    else
       meanzii = sum(zsetgas(1:ninbin(ii),ii))/real(ninbin(ii))
    endif
    H(ii) = sqrt(sum(((zsetgas(1:ninbin(ii),ii)-meanzii)**2)/(real(ninbin(ii)-1))))
 enddo

end subroutine calculate_H

!--------------------------------------------------------
!+
! Deallocate arrays
!+
!-------------------------------------------------------
subroutine deallocate_arrays
 implicit none

 deallocate(gravxyz)
 deallocate(neighcount,neighb)
 deallocate(ipartbin,ninbin,rad)
 deallocate(rpart,phipart,vrpart,vphipart)
 deallocate(gr,gphi,Br,Bphi,vrbin,vphibin)
 deallocate(sigma,csbin,H,toomre_q,omega,epicyc)
 deallocate(alpha_reyn,alpha_grav,alpha_mag,alpha_art)
 deallocate(part_scaleheight,h_smooth)
 deallocate(mean_z,mean_vz,rms_z,rms_vz)
 if (allocated(tcool)) deallocate(tcool)


end subroutine deallocate_arrays
!-------------------------------------------------------
end module analysis
