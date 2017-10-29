!******************************************************************************
!****m* /NNPION_NN
! NAME
! module NNPION_NN
! PURPOSE
! Implements the Matrix elements for N N Pion -> N N.
!******************************************************************************
module NNPION_NN

  implicit none
  private

  public :: gamma_NNPion_NN, matrix_NN_NNPion

contains

  !****************************************************************************
  !****f* NNPION_NN/gamma_NNPion
  ! NAME
  ! function gamma_NNPion_NN(srts,Epion,Enucleon,rhoProton,rhoNeutron,chargePion,chargeNucleon,OutPutFlag) Result(gamma)
  ! PURPOSE
  ! Evaluates Gammas for pion absorption on a pair of nucleons. There exist the following possible channels:
  ! * 1: p p pi^0-> p p  ; which is equivalent to n n -> n n pi^0
  ! * 2: n n pi^+-> n p  ; which is equivalent to  p n -> p p pi^-
  ! * 3: p n pi^0-> p n  ;
  ! * 4: p n pi^+-> p p  ; which is equivalent to n n -> p n pi^-
  ! OUTPUT
  ! * real :: gamma ! absorption rate of pion on one specific pair on nucleons of the input
  ! INPUTS
  ! * real, intent(in) :: srts ! SQRT(s)
  ! * real, intent(in) :: Epion ! Energy of pion
  ! * real, intent(in),dimension(1:2) :: Enucleon !  Energies of nucleons
  ! * real, intent(in) :: rhoProton
  ! * real, intent(in) :: rhoNeutron
  ! * integer, intent(in) :: chargePion
  ! * integer, dimension(1:2), intent(in) :: chargeNucleon
  ! * real,intent(out) :: gamma
  ! * logical, optional,intent(in) :: outputFlag
  !****************************************************************************
  function gamma_NNPion_NN(srts,Epion,Enucleon,rhoProton,rhoNeutron,chargePion,chargeNucleon,OutPutFlag) Result(gamma)

    use constants, only: pi, mN, mPi, hbarc, GeVSquared_times_mb
    use twoBodyTools, only: pCM

    real, intent(in) :: srts, Epion, Enucleon(1:2), rhoProton, rhoNeutron
    integer, intent(in) :: chargePion, chargeNucleon(1:2)
    logical, optional, intent(in) :: outputFlag
    real :: gamma

    integer :: k
    real, dimension(1:4) :: matrixelements,rhoTwo
    logical, save,dimension(-1:1) :: openFlag=.true.
    character(20), dimension(-1:1), parameter :: fileName = (/ 'gamma_PiMinus.dat', &
                                                               'gamma_PiNull.dat ', &
                                                               'gamma_PiPlus.dat ' /)

    rhoTwo=0.
    gamma=0.

    select case (chargePion)
    case (-1)
      select case (sum(chargeNucleon))
      case (2)  ! pi- pp -> pn
        rhoTwo(2)=rhoProton**2/2.
        k=2
      case (1)  ! pi- pn -> nn
        rhoTwo(4)=rhoProton*rhoNeutron
        k=4
      case default
        return
      end select
    case (0)
      select case (sum(chargeNucleon))
      case (2)  ! pi0 pp -> pp
        rhoTwo(1)=(rhoProton**2)/2.
        k=1
      case (0)  ! pi0 nn -> nn
        rhoTwo(1)=(rhoNeutron**2)/2.
        k=1
      case (1)  ! pi0 pn -> pn
        rhoTwo(3)=rhoProton*rhoNeutron
        k=3
      end select
    case (1)
      select case (sum(chargeNucleon))
      case (0)  ! pi+ nn -> pn
        rhoTwo(2)=rhoNeutron**2/2.
        k=2
      case (1)  ! pi+ pn -> pp
        rhoTwo(4)=rhoProton*rhoNeutron
        k=4
      case default
        return
      end select
    end select

    matrixElements = matrix_NN_NNPion(srts)

    gamma = Min(100000., gamma+matrixElements(k)*pCM(srts,mN,mN)/ srts/4. /pi*rhoTwo(k)/8./ePion/eNucleon(1)/eNucleon(2) &
                         * hbarc**6 * GeVSquared_times_mb)

     ! Make output
     if (present(outputFlag)) then
        if (OutPutFlag) then
           if (openFlag(chargePion)) then
              open(11, File=fileName(chargePion))
              openFlag(chargePion)=.false.
           else
              open(11, File=fileName(chargePion),position='Append')
           end if
           write(11,'(5F15.3,2I5)') ePion-mPi,eNucleon,srts,gamma,chargeNucleon
           close(11)
        end if
     end if

  end function gamma_NNPion_NN


  !****************************************************************************
  !****f* NNPION_NN/gamma_NNPion_NN_output
  ! NAME
  ! function gamma_NNPion_NN_output(srts,Epion,rhoProton,rhoNeutron,chargePion,OutPutFlag) result (gamma)
  ! NOTES
  ! Evaluates Gamma for the channels:
  ! * 1: p p pi^0-> p p  ; which is equivalent to n n -> n n pi^0
  ! * 2: n n pi^+-> n p  ; which is equivalent to  p n -> p p pi^-
  ! * 3: p n pi^0-> p n  ;
  ! * 4: p n pi^+-> p p  ; which is equivalent to n n -> p n pi^-
  ! The nucleons are considered to rest.
  ! This routine is mainly suited for output.
  ! OUTPUT
  ! * real, dimension(1:4),intent(out) :: gamma
  ! INPUTS
  ! * real, intent(in) :: srts ! SQRT(s)
  ! * real, intent(in) :: Epion ! Energy of pion
  ! * real, intent(in) :: rhoProton
  ! * real, intent(in) :: rhoNeutron
  ! * integer, intent(in) :: chargePion
  ! * logical, optional,intent(in) :: outputFlag  -  if .true. then results are written to files "gamma_NNPiPlus_NN.dat","gamma_NNPiNull_NN.dat","gamma_NNPiMinus_NN.dat"
  !****************************************************************************
!   function gamma_NNPion_NN_output(srts,Epion,rhoProton,rhoNeutron,chargePion,OutPutFlag) result (gamma)
!
!     use constants, only: pi, mN, mPi, hbarc, GeVSquared_times_mb
!     use twoBodyTools, only: pCM
!
!     real, intent(in) :: srts, Epion, rhoProton, rhoNeutron
!     integer, intent(in) :: chargePion
!     logical, optional,intent(in) :: outputFlag
!     real, dimension(1:4) :: gamma
!
!     real :: Enucleon(1:2), matrixelements(1:4), rhoTwo(1:4)
!     logical, save, dimension(-1:1) :: openFlag=.true.
!     character(25), dimension(-1:1), parameter :: fileName = (/ 'gamma_NNPiMinus_NN.dat', &
!                                                                'gamma_NNPiNull_NN.dat ', &
!                                                                'gamma_NNPiPlus_NN.dat ' /)
!
!     matrixElements = matrix_NN_NNPion(srts)
!
!     eNucleon=mN  ! nucleons rest
!
!     Select Case(chargePion)
!     Case(-1)
!        rhoTwo(1)=0. ! channel not open
!        rhoTwo(2)=rhoProton**2/2.
!        rhoTwo(3)=0. ! channel not open
!        rhoTwo(4)=rhoProton*rhoNeutron
!     Case(0)
!        rhoTwo(1)=(rhoProton**2+rhoNeutron**2)/2.
!        rhoTwo(2)=0. ! channel not open
!        rhoTwo(3)=rhoProton*rhoNeutron
!        rhoTwo(4)=0. ! channel not open
!     Case(1)
!        rhoTwo(1)=0. ! channel not open
!        rhoTwo(2)=rhoNeutron**2/2.
!        rhoTwo(3)=0. ! channel not open
!        rhoTwo(4)=rhoProton*rhoNeutron
!     end select
!
!     gamma(1:4) = Min(100000., matrixElements(1:4)*pCM(srts,mN,mN)/ srts/4. /pi*rhoTwo(1:4)/8./ePion/eNucleon(1)/eNucleon(2) &
!                               * hbarc**6 * GeVSquared_times_mb)
!
!      ! Make output
!      IF(Present(outputFlag)) then
!         If(OutPutFlag) then
!            If(Openflag(chargePion)) then
!               Open(11, File=fileName(chargePion))
!               openFlag(chargePion)=.false.
!               write(11,'(A)') '# Ekin of pion,mass of nucleon, energy of Nucleon, srts,gamma(1:4), sum(gamma)'
!               write(11,'(A)') '# gamma(1): p p pi^0-> p p  +  n n pi^0 -> nn'
!               write(11,'(A)') '# gamma(2): n n pi^+-> n p  ; which is equivalent to  p p pi^- ->pn '
!               write(11,'(A)') '# gamma(3): p n pi^0-> p n   '
!               write(11,'(A)') '# gamma(4): p n pi^+-> p p  ; which is equivalent to p n pi^- -> nn'
!            else
!               Open(11, File=fileName(chargePion),position='Append')
!            end if
!            write(11,'(9F15.6)') ePion-mPi,eNucleon,srts,gamma, Sum(gamma)
!            close(11)
!         end if
!      end IF
!   end function gamma_NNPion_NN_output


  !****************************************************************************
  !****f* NNPION_NN/matrix_NN_NNPion
  ! NAME
  ! function matrix_NN_NNPion(srts,outputFlag) result (matrixElements)
  ! PURPOSE
  ! Implements the Matrix elements for N N Pion -> N N in units of  "GeV**2 * mB"  for the channels:
  ! * 1: p p -> p p pi^0 ; which is equivalent to n n -> n n pi^0
  ! * 2: p n -> n n pi^+ ; which is equivalent to  p n -> p p pi^-
  ! * 3: p n -> p n pi^0 ;
  ! * 4: p p -> p n pi^+ ; which is equivalent to n n -> p n pi^-
  ! INPUTS
  ! * real, intent(in) ::srts
  ! OUTPUT
  ! * real, dimension(1:4), intent(out) :: matrixElements
  !****************************************************************************
  function matrix_NN_NNPion(srts,outputFlag) result (matrixElements)

    use barBar_barbarMes, only: NN_NNpi_direct
    use threeBodyPhaseSpace, only: Integrate_3bodyPS
    use constants, only: pi, mN, mPi
    use twoBodyTools, only: pCM

    real, intent(in) ::srts
    logical, optional, intent(in) ::outputFlag
    real, dimension(1:4) :: matrixElements

    real, dimension (1:4) :: Sigma_NNPion
    real :: ps
    logical, save :: openFlag=.true.

    ! detailed balance: pi N N -> N N
    ! (see p.46 of Effenbergers diploma thesis)
    ! or page 56 of Buss' thesis
    ps = Integrate_3bodyPS (srts, mN, mN, mPi)
    ! ps=\int dm12 dm23 /s (!)

    Sigma_NNPion = NN_NNpi_direct (srts)
    ! Field Sigma_NNPion:
    ! 1: p p -> p p pi^0 ; which is equivalent to n n -> n n pi^0
    ! 2: p n -> n n pi^+ ; which is equivalent to  p n -> p p pi^-
    ! 3: p n -> p n pi^0 ;
    ! 4: p p -> p n pi^+ ; which is equivalent to n n -> p n pi^-

    if (ps>1E-12) then
      matrixElements(1:4) = 64.*(2.*pi)**3*pCM(srts,mN,mN)*srts*Sigma_NNPion(1:4)/ps &   ! units GeV**2 mB
                            * (/1.,2.,1.,0.5/)                                           ! factor for identical particles
    else
      matrixElements(1:4) = 0.
    end if

     ! Make output
     if (present(outputFlag)) then
        if (outputFlag) then
           if (Openflag) then
              open(11, File='matrixelements_NN_NNPION.dat')
              openFlag=.false.
           else
              open(11, File='matrixelements_NN_NNPION.dat',position='append')
           end if
           write(11,'(5F25.9)') srts, matrixElements
           close(11)
        end if
     end if

   end function matrix_NN_NNPion

end module NNPION_NN
