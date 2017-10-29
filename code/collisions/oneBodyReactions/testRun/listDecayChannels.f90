!*******************************************************************************
! some test program to perform MC decays of a given particle and list all decay
! channels.
!
! here done for the phi meson, but could be extended to all other unstable
! particles.
!
! needs jobcard 'listDecayChannels.job' for switching on 'varirate'
!*******************************************************************************

program main
  use inputGeneral, only: readInputGeneral
  use version, only: printVersion
  use particleProperties, only: initParticleProperties, hadron
  use IdTable, only: phi
  use particleDefinition
  use master_1body, only: decayParticle
  use output, only: writeParticle
  use twoBodyStatistics, only: rate
  use callstack, only: system_command

  implicit none

  type(particle) :: resonance
  type(particle),dimension(1:20) :: finalState
  logical :: pauliFlag, collisionFlag
  integer :: iMC
  integer, parameter :: nMC = 1000000



  call PrintVersion
  call readInputGeneral
  call initParticleProperties

  call setToDefault(resonance)

  resonance%ID = phi
  resonance%mass = hadron(phi)%mass
  resonance%momentum(0) = resonance%mass

  call setToDefault(finalState)

  do iMC = 1,nMC
     call decayParticle(resonance, finalState, collisionFlag, pauliFlag, .true., 99.9)

     if (.not.collisionFlag) then
        write(*,*) 'collisionFlag = .false.'
     end if

!     call writeParticle(6,0,finalState)

     call rate((/resonance/),finalState,99.9)

  end do

  ! this is tricky: a different time enforces producing output of the stats!
  call rate((/resonance/),finalState,999.9) ! do output

  ! now we print everything to screen by just printing the file:

  write(*,*) '# number events:',nMC
  call system_command("cat VariRate.rates.dat")

end program main
