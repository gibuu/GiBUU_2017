! test case for Ticket #72
! http://gibuu.physik.uni-giessen.de/GiBUU/ticket/72
!
! when compiled with FORT=gfortran ARCH=32, running "./test_ticket72.x < job" produces the following backtrace:
!
! Backtrace for this error:
!   + function __quadpack_MOD_qage (0x8CCAF16)
!     at line 689 of file quadpack.f90
!   + function __quadpack_MOD_qag (0x8CCB022)
!     at line 336 of file quadpack.f90
!   + function __mesonwidthvacuum_MOD_semistablefinalstate (0x85F27D4)
!     at line 349 of file mesonWidthVacuum.f90
!   + function __mesonwidthvacuum_MOD_vacuumwidth (0x85F403C)
!     at line 135 of file mesonWidthVacuum.f90
!   + function test (0x80482F9)
!     at line 18 of file test_ticket72.f90
!   + function __libc_start_main (0x8DEA780)
!     at line 258 of file libc-start.c

program test
  use IDTable, only: phi
  use particleProperties, only: InitParticleProperties, meson
  use decayChannels
  use mesonWidthVacuum, only: vacuumWidth
  implicit none

  real :: gammaTotal
  real, dimension(1:size(decays2body_meson)) :: ratio2
  real, dimension(1:size(decays3body_meson)) :: ratio3

  call InitParticleProperties

  call vacuumWidth (meson(phi)%minmass, phi, gammaTotal, ratio2, ratio3)

end program test 
