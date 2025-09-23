!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! Rotate particle coordinates to align disc plane with xy-plane
! and center disc on most massive sink (or COM if no sinks)
!
! :References: None
!
! :Owner: Adam Koval
!
! :Runtime parameters: None
!
! :Dependencies: 
!
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

contains

subroutine modify_dump(npart, npartoftype, massoftype, xyzh, vxyzu)
 use readwrite_dumps, only:write_smalldump, write_fulldump
 use misalignedutils, only:rotate_align
 use prompting,       only:prompt

 integer, intent(in) :: npart, npartoftype(:)
 real, intent(in) :: massoftype(:)
 real, intent(inout) :: xyzh(:,:), vxyzu(:,:)

 real :: sphere_radius
 sphere_radius = 100.0 ! default value

 print*,' *** Aligning disc in dump file ***'
 call prompt('What radius (in au) should be used to select particles for angular momentum calculation? (default 100 au):',&
             sphere_radius)
 call rotate_align(xyzh, vxyzu, npart, sphere_radius)

end subroutine modify_dump

end module moddump

