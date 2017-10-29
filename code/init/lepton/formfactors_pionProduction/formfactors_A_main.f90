!******************************************************************************
!****m* /formfactors_A_main
! NAME
! module formfactors_A_main
!
! PURPOSE
! This module administrates the formfactors for gamma* N -> N pi.
!
! USES
! * formfactors_A_pi0_neutron
! * formfactors_A_pi0_protron
! * formfactors_A_piPlus_neutron
! * formfactors_A_piMinus_proton
!
! INPUTS
! None
!******************************************************************************
module formfactors_A_main

  implicit none
  private

  public :: getA, cleanUp

contains

  !****************************************************************************
  !****f* formfactors_A_main/getA
  ! NAME
  ! function getA(pionCharge_out,nucCharge_out,thetaIn,sIn,QSquaredIn)
  !
  ! PURPOSE
  ! This function returns the invariant Amplitudes A_1, ... A_6 of the MAID analysis.
  !
  ! INPUTS
  ! * real  ::  thetaIn    ! Theta scattering angle of the pion relative to q in the CM system of the hadronic vertex
  ! * real  ::  QSquaredIn ! Q^2=-q^mu q_mu  of the gamma
  ! * real  ::  sIn        ! Mandelstam s of the hadronic vertex: gamma* N-> pi N
  ! * integer :: pionCharge_out,nucCharge_out ! charges of outgoing pion and nucleon
  !
  ! OUTPUT
  ! complex, dimension(1:6) :: getA   ! The complex amplitudes A_1 to A_6 in units of:
  ! * [GeV**-2] for A(1)
  ! * [GeV**-4] for A(2) and A(5)
  ! * [GeV**-3] for A(3), A(4) and A(5)
  !
  !****************************************************************************
  function getA(pionCharge_out,nucCharge_out,thetaIn,sIn,QSquaredIn)
    use formfactors_A_piZero_neutron, only: getA_piZero_neutron
    use formfactors_A_piZero_proton, only: getA_piZero_proton
    use formfactors_A_piPlus_neutron, only: getA_piPlus_neutron
    use formfactors_A_piMinus_proton, only: getA_piMinus_proton

    complex, dimension(1:6) :: getA
    real  ::  thetaIn    ! Theta scattering angle of the pion relative to q in the CM system of the hadronic vertex
    real  ::  QSquaredIn,QSquared ! Q^2=-q^mu q_mu  of the gamma
    real  ::  sIn        ! Mandelstam s of the hadronic vertex: gamma* N-> pi N
    integer :: pionCharge_out,nucCharge_out ! charges of outgoing pion and nucleon

    if (QsquaredIn.le.0) then
       if (QsquaredIn.lt.-1E-4) then
          write(*,*) "Warning Q^2 less than zero in getA!", QsquaredIn
       end if
       QSquared=1E-12
    else
       QSquared=QSquaredIn
    end if

    select case (pionCharge_out)
    case (1)
       select case (nucCharge_out)
       case (0)
          getA=getA_piPlus_neutron(thetaIn,sIn,QSquared)
       case default
          getA=0.
          call error()
       end select
    case (-1)
       select case (nucCharge_out)
       case (1)
          getA=getA_piMinus_proton(thetaIn,sIn,QSquared)
       case default
          getA=0.
          call error()
       end select
    case (0)
       select case (nucCharge_out)
       case (0)
          getA=getA_piZero_neutron(thetaIn,sIn,QSquared)
       case (1)
          getA=getA_piZero_proton(thetaIn,sIn,QSquared)
       case default
          getA=0.
          call error()
       end select
    case default
       call error()
    end select

  contains

    subroutine error()
      write(*,'(A)') 'Invalid charges in formfactors_A_main/getA'
      write(*,'(A,I4)') 'Outgoing pion charge   =', pionCharge_out
      write(*,'(A,I4)') 'Outgoing nucleon charge=', nucCharge_out
      write(*,'(A)') 'Severe problem: Stopping'
      stop
    end subroutine error

  end function getA


  subroutine cleanUp
    use formfactors_A_piZero_neutron, only: cleanupPi0n => cleanup
    use formfactors_A_piZero_proton, only: cleanupPi0p => cleanup
    use formfactors_A_piPlus_neutron, only: cleanupPipn => cleanup
    use formfactors_A_piMinus_proton, only: cleanupPimp => cleanup
    call cleanupPi0n
    call cleanupPi0p
    call cleanupPipn
    call cleanupPimp
  end subroutine


end module formfactors_A_main
