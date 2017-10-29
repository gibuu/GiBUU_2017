!******************************************************************************
!****m* /electronPionProduction_kine
! NAME
! module  electronPionProduction_kine
! PURPOSE
! Contains kinematical considerations for gamma* nucleon -> pion nucleon events.
!******************************************************************************
module  electronPionProduction_kine
  use minkowski, only: pair => SP

  implicit none

  private

  public :: getV_out
  public :: get_dV_Pi_dk
  public :: get_k_abs_improved
  public :: get_k_abs
  public :: getKinematics_eN

  logical, save :: debug=.false.
  integer, save :: numberTries=0

  !****************************************************************************
  !****g* electronPionProduction_kine/pionPot
  ! SOURCE
  !
  logical ,save :: pionPot=.true.
  !
  ! PURPOSE
  ! * Switch pion potential explicitly on and off in this module.
  ! * Only for debugging.
  !****************************************************************************
  public :: pionPot

contains

  !****************************************************************************
  !****s* electronPionProduction_kine/getKinematics_eN
  ! NAME
  ! subroutine getKinematics_eN(eN, pionCharge, nucleon_out_charge,
  ! phi_k, theta_k, k, pf, twoRoots, success, pionNucleonSystem, bothSolutions)
  !
  ! PURPOSE
  ! Evaluates the full kinematics for a
  ! electron nucleon -> electron pion nucleon reaction.
  !
  ! The angles of the outgoing pion are measured relative to the momentum
  ! transfer "q". All angles in degree.
  !
  ! INPUTS
  ! * type(electronNucleon_event) :: eN  -- underlying electron nucleon event
  ! * integer                     :: pionCharge --
  ! * integer                     :: nucleon_out_charge --
  ! * real                        :: phi_k,theta_k -- Pion angles [degrees]
  ! * integer           , OPTIONAL:: pionNucleonSystem --
  !   1 = pion angles phi_k and theta_k are defined in lab,
  !   2 = they are defined in the CM frame of the incoming nucleon and photon
  !
  ! RESULT
  ! * real, dimension(0:3) :: k    -- Momentum of outgoing pion
  ! * real, dimension(0:3) :: pf   -- Momentum of outgoing nucleon
  ! * logical              :: twoRoots --
  !   true if there were two possible kinematics
  ! * logical              :: success --
  !   Flag shows whether the kinematics could be established.
  ! * real, dimension(1:2,0:3), OPTIONAL :: bothSolutions --
  !   if two roots are found then this variable returns both solutions
  !
  ! NOTES
  ! Especially if the pionNucleonSystem is not chosen to be the CM-frame (=2),
  ! then two solutions for a given pair (phi_k,theta_k) may exist. If this is
  ! the case, then we choose randomly one of them and set twoRoots=.true.!
  !****************************************************************************
  subroutine getKinematics_eN(eN,pionCharge,nucleon_out_charge,phi_k,theta_k,&
       & k,pf,twoRoots,success,pionNucleonSystem,bothSolutions)
    use vector, only: absVec
    use particleDefinition
    use random, only: rn_trueFalse
    use degRad_conversion, only: degrees, radian
    use rotation, only: rotateTo
    use eN_eventDefinition, only: electronNucleon_event
    use lorentzTrafo, only: lorentz

    type(electronNucleon_event) , intent(in) :: eN
    integer                     , intent(in) :: pionCharge,nucleon_out_charge
    real                        , intent(in) :: phi_k,theta_k
    integer           , optional, intent(in) :: pionNucleonSystem
    real, dimension(0:3)        , intent(out):: k,pf
    logical                     , intent(out):: success
    logical                     , intent(out):: twoRoots
    real, dimension(1:2,0:3)    , intent(out), optional :: bothSolutions

    real :: k_abs
    real, dimension(1:3) :: k_unit
    real, dimension(0:3) :: q
    real, dimension(1:3) :: betaTOCM
    real, dimension (0:3) :: pcm, qcm
    type(particle) :: nucleonCM
    real, dimension (1:2,0:3) :: k2
    integer :: numRoots

    logical,parameter :: debug_kin=.false.
    logical           :: pionAngles_inCM

    if (present(pionNucleonSystem)) then
       select case (pionNucleonSystem)
       case (1)
          pionAngles_inCM=.false.
       case (2)
          pionAngles_inCM=.true.
       case default
          write(*,*) 'Wrong value for pionNucleonSystem in getKinematics_EN',pionNucleonSystem
          write(*,*) 'STOP!'
          stop
       end select
    else
       pionAngles_inCM=.false.
    end if

    success=.false.

    ! 1) Photon Momentum
    q=eN%boson%momentum

    if (debug_kin) then
       write(*,*)
       write(*,'(A)')'###############################################'
       write(*,'(A)') 'Virtual photon :'
       write(*,'(A,4F9.5)') 'q=',q
       write(*,*)
       write(*,'(A,F9.5)') 'Q^2=',-pair(q,q)
       write(*,*)
    end if

    ! Define unit vector in direction of outgoing pion:
    ! Here z-axis is chosen in the direction of q(1:3)
    betaTOCM(1:3)=(q(1:3)+eN%nucleon%momentum(1:3))/(q(0)+eN%nucleon%momentum(0))
    if (pionAngles_inCM) then
       ! Boost event to CM
       qCM=q
       call lorentz(betaTOCM, qCM, ' getKinematics_eN')
       pCM=eN%nucleon%momentum
       call lorentz(betaTOCM, pCM, ' getKinematics_eN')
       nucleonCM=eN%nucleon
       nucleonCM%momentum=pCM
       ! Construct pion unit vector in CM system
       k_unit(1:3)=(/sin(radian(theta_k))*cos(radian(phi_k)),&
            &        sin(radian(theta_k))*sin(radian(phi_k)),&
            &        cos(radian(theta_k))/)

       ! Rotate this pion vector to a system where the CM q is defining the z-axis
       k_unit = rotateTo (qCM(1:3), k_unit)
       if (debug) write(*,'(A,F9.5)') 'theta_pion=',degrees(acos(Dot_product(k_unit,q(1:3)) &
            & /sqrt(Dot_product(q(1:3),q(1:3)))))

       ! Solve energy momentum conservation to get absolute value of k(1:3):
       k_abs=get_K_abs(k_unit(1),k_unit(2),k_unit(3),qCM,nucleonCM,pionCharge,nucleon_out_charge,success,k,betaToCF_in=-betaToCM)

       ! Boost k back
       call lorentz(-betaTOCM, k, ' getKinematics_eN')
       twoRoots=.false.
    else
       k_unit(1:3)=(/sin(radian(theta_k))*cos(radian(phi_k)),&
            &        sin(radian(theta_k))*sin(radian(phi_k)),&
            &        cos(radian(theta_k))/)
       ! Rotate this pion vector to a system where q is defining the z-axis
       k_unit = rotateTo (q(1:3), k_unit)
       if (debug) write(*,'(A,F9.5)') 'theta_pion=',degrees(acos(Dot_product(k_unit,q(1:3)) &
            & /sqrt(Dot_product(q(1:3),q(1:3)))))
       ! Solve energy momentum conservation to get absolute value of k(1:3):
       ! Here we need to care, since in the calculation frame there may exist two valid |k|'s for one value of (theta_lab,phi_lab). This
       ! comes due to the boost of the CM-frame: given a large enough boost, even those events which travel backwards in CM frame of the pion
       ! nucleon system travel forwards in the calculation frame. So there are two |k|'s for one value of (theta_lab,phi_lab): one is the CM
       ! frame forward-solution and the other is the backward-solution.
       call get_k_abs_improved(numRoots,k_unit(1),k_unit(2),k_unit(3),q,eN%nucleon,pionCharge,  &
            &  nucleon_out_charge,success,k2,(/0.,0.,0./))
       if (numRoots.eq.1) then
          k=k2(1,:)
          twoRoots=.false.
       else if (numRoots.eq.2) then
          ! Monte Carlo Decision:
          ! We choose randomly one of the solutions.
          if (rn_trueFalse()) then
             k=k2(1,:)
          else
             k=k2(2,:)
          end if
          twoRoots=.true.
       else
          k=0.
          twoRoots=.false.
       end if
       K_abs=absVec(k(1:3))
    end if

    numberTries=numberTries+1
    if (success) then
       ! Check solutions
       if (k_abs.lt.0.0000001) then
          write(*,*) 'Warning: Outgoing pion momentum 0!'
          write(*,*) "k_abs=", k_abs
       end if
       pf=q+eN%nucleon%momentum-k
    else
       k=0
       pf=0
    end if
    if (debug_kin) write(*,*) "getKinematics_eN: k_abs=", k_abs,success

    if (present(bothSolutions)) bothSolutions=k2

  end subroutine getKinematics_eN

  !****************************************************************************
  !****s* electronPionProduction_kine/get_k_abs_improved
  ! NAME
  ! function get_k_abs_improved(numRoots, ak, bk, ck, q, initNuc, pionCharge,
  ! nucleon_out_charge, success, kout, betaToCF)
  !
  ! PURPOSE
  ! Evaluates for gamma nucleon -> pion nucleon the absolute value of the
  ! pion momentum, when given a unit vector (ak,bk,ck) in its direction.
  ! The value of gamma momentum "q", the intial nucleus "initNuc" and the
  ! charge of the pion "pionCharge" are input.
  !
  ! The input momenta are given in the calculation system, and also the
  ! output k is given in this system.
  !
  ! INPUTS
  ! * real,intent(in) :: ak,bk,ck -- unit-vector (ak,bk,ck)
  !   in direction of pion momentum
  ! * real, dimension(0:3) :: q --
  ! * type(particle)  :: initNuc --
  ! * integer         :: pionCharge --
  ! * integer         :: nucleon_out_charge --
  ! * real, dimension(1:3), OPTIONAL  :: betaToCF --
  !   Velocity of CF frame in the frame where initNuc, q and ak,bk,ck are
  !   defined
  !
  ! OUTPUT
  ! * integer :: numRoots  -- Number of roots
  ! * logical :: success   -- Flag if the kinematics could be established.
  ! * real, dimension(1:2,0:3) :: kout -- found roots
  !
  ! NOTES
  ! We use a Newton-Algorithm to solve the energy and momentum conservation
  ! condition. Searches also for two roots!!
  !
  ! Solves q(0)+pi(0)=pf(0)+k(0) at the position "position".
  !****************************************************************************
  subroutine get_k_abs_improved(numRoots,ak,bk,ck,q,initNuc,pionCharge,&
       & nucleon_out_charge,success,kout,betaToCF)
    use particleDefinition
    use random

    real               , intent(in) :: ak,bk,ck
    type(particle)     , intent(in) :: initNuc
    integer            , intent(in) :: pionCharge,nucleon_out_charge
    real, dimension(0:3),intent(in) :: q
    real, dimension(1:3),intent(in) :: betaToCF

    logical                           , intent(out) :: success
    real, dimension(1:2,0:3)          , intent(out) :: kout       ! found roots
    integer                           , intent(out) :: numRoots   ! Number of roots

    real, parameter :: eps=0.005
    real :: secondRoot, firstRoot, low, high

    firstRoot= get_k_abs(ak,bk,ck,q,initNuc,pionCharge,nucleon_out_charge,success,kout(1,:),betaToCF,0.,100000.)
    if (.not.success) then
       numRoots=0
       kout=0.
       success=.false.
       return
    else
       ! Search for second root:
       low=firstRoot+eps
       high=1000000000000000.
       secondRoot=get_k_abs(ak,bk,ck,q,initNuc,pionCharge,nucleon_out_charge,success,kout(2,:),betaToCF,low,high)
       if (success) then
          numRoots=2
          success=.true.
          !write(*,*) 'number of roots',numRoots, firstRoot,secondRoot
          return
       end if
       low=0.
       high=firstRoot-eps
       secondRoot=get_k_abs(ak,bk,ck,q,initNuc,pionCharge,nucleon_out_charge,success,kout(2,:),betaToCF,low,high)
       if (success) then
          numRoots=2
          success=.true.
          !write(*,*) 'number of roots',numRoots, firstRoot,secondRoot
          return
       else
          numRoots=1
          kout(2,:)=0.
          success=.true.
          !write(*,*) 'number of roots',numRoots, firstRoot,secondRoot
       end if
    end if


  end subroutine get_k_abs_improved


  !****************************************************************************
  !****s* electronPionProduction_kine/get_k_abs
  ! NAME
  ! real function get_k_abs(ak, bk, ck, q, initNuc, pionCharge,
  ! nucleon_out_charge, success, kout, betaToCF, low, high)
  !
  ! PURPOSE
  ! Evaluates for gamma nucleon -> pion nucleon the absolute value of the
  ! pion momentum, when given a unit vector (ak,bk,ck) in its direction.
  ! The value of gamma momentum "q", the intial nucleus "initNuc" and the
  ! charge of the pion "pionCharge" are input.
  !
  ! The input (and output) momenta are given in the calculation system.
  !
  ! INPUTS
  ! * real,intent(in) :: ak,bk,ck -- unit-vector (ak,bk,ck)
  !   in direction of pion momentum
  ! * real, dimension(0:3) :: q --
  ! * type(particle)  :: initNuc --
  ! * integer         :: pionCharge --
  ! * integer         :: nucleon_out_charge --
  ! * real, dimension(1:3), OPTIONAL  :: betaToCF --
  !   Velocity of CF frame in the frame where initNuc, q and ak,bk,ck are
  !   defined
  ! * real   ,OPTIONAL:: low,high --
  !
  ! OUTPUT
  ! * logical :: success   -- Flag if the kinematics could be established.
  ! * real, dimension(0:3), OPTIONAL  :: kout --
  !
  ! NOTES
  ! * We use a Newton-Algorithm to solve the energy and momentum conservation
  !   condition.
  ! * DOES NOT WORK IF THERE IS MORE THAN ONE SOLUTION, which happens if
  !   there are large boosts involved, i.e. if |q| and |q_0| are large.
  ! * Solves q(0)+pi(0)=pf(0)+k(0) at the position "position".
  !****************************************************************************
  real function get_k_abs(ak,bk,ck,q,initNuc,pionCharge,nucleon_out_charge,&
       & success,kout,betaToCF_in,low_in,high_in)
    use particleDefinition
    use random

    real,intent(in) ::  ak, bk,ck
    logical, intent(out) :: success
    type(particle) ,intent(in) :: initNuc
    integer , intent(in):: pionCharge,nucleon_out_charge

    real, dimension(0:3),intent(in) :: q

    real, dimension(0:3),intent(out),optional  :: kout
    real, dimension(1:3),intent(in) ,optional  :: betaToCF_in

    real, dimension(0:3) :: k
    real, dimension(1:3) :: f, k_abs
    real, dimension(1:3) :: betaToCF
    real, optional, intent(in) :: low_in,high_in

    real :: low,high
    real, parameter :: maxError=0.0005
    integer,parameter :: maxcounter=10000
    !integer,parameter :: maxcounter=5000
    logical, parameter :: debugFlag=.false.
    logical, parameter :: writeFlag=.false. ! Write the iteration point (kabs,errorEnergyConservation) to fort.12
    integer :: counter
    integer :: numRestarts

    low = 0.     !=0 eV as lower bound for abs(vec(k))
    if (present(low_in)) low=low_in

    high = 1000. !=1TeV upper bound
    if (present(high_in)) high=high_in

    betaToCF=0.
    if (present(betaToCF_in)) betaToCF=betaToCF_in

    ! Default return values: failure
    success = .false.
    get_k_abs = 1111111111.

    !Make some ansatz to evaluate starting values
    k_abs(1)=low+0.2
    if (.not.set_k(k_abs(1),k,initNuc,pionCharge,betaToCF)) then
!       write(*,*) 'k_abs Problem 1'
       return ! FAILURE !!!!
    end if
    f(1)=error_energyConservation(k,q,initNuc,pionCharge,nucleon_out_charge,betaToCF,k_abs(1),low,high)
    if (writeFlag) write(12,*) k_abs(1),f(1)

    if (abs(f(1)).lt.maxError.and.k_abs(1).gt.0) then !  Successful!!!
       if (DebugFlag) write(*,*) 'SuccesFull with f(1)',f(1),k_abs(1)
       get_k_abs=k_abs(1)
       success=.true.
       if (present(kout)) kout=k
       return
    end if

    k_abs(2)=low+0.5
    if (.not.set_k(k_abs(2),k,initNuc,pionCharge,betaToCF)) then
!       write(*,*) 'k_abs Problem 2'
       return ! FAILURE !!!!
    end if
    f(2)=error_energyConservation(k,q,initNuc,pionCharge,nucleon_out_charge,betaToCF,k_abs(2),low,high)
    if (writeFlag) write(12,*) k_abs(2),f(2)

    if (abs(f(2)).lt.maxError.and.k_abs(2).gt.0) then !  Successful!!!
       if (DebugFlag) write(*,*) 'SuccessFull with f(2)',f(2),k_abs(2)
       get_k_abs=k_abs(2)
       success=.true.
       if (present(kout)) kout=k
       return
    end if

    if (debugFlag) Print * , 'Begin Iteration'

    ! Begin Iteration

    counter=0

    numRestarts=0
    do
       if (abs(f(2)-f(1)).eq.0) then
          ! Make new start guesses: This case usually occurs if the solution for k_abs is negative.
          k_abs(1)=rn()*10.
          k_abs(2)=rn()*10.
          if (.not. set_k(k_abs(1),k,initNuc,pionCharge,betaToCF)) then
!             write(*,*) 'k_abs Problem x1: ',k_abs(1)
             return ! FAILURE !!!!
          end if
          f(1)=error_energyConservation(k,q,initNuc,pionCharge,nucleon_out_charge,betaToCF,k_abs(1),low,high)

          if (.not. set_k(k_abs(2),k,initNuc,pionCharge,betaToCF)) then
!             write(*,*) 'k_abs Problem x2: ',k_abs(2)
             return ! FAILURE !!!!
          end if
          f(2)=error_energyConservation(k,q,initNuc,pionCharge,nucleon_out_charge,betaToCF,k_abs(2),low,high)

          numRestarts=numRestarts+1
          if (numRestarts.eq.100) then
             if (present(kout)) kout=k
             return ! FAILURE !!!!
          end if
       end if

       k_abs(3)=k_abs(2)-f(2)*(k_abs(2)-k_abs(1))/(f(2)-f(1))
       if (.not. set_k(k_abs(3),k,initNuc,pionCharge,betaToCF)) then
!          write(*,*) 'k_abs Problem x3: ',k_abs(3)
          return ! FAILURE !!!!
       end if
       f(3)=error_energyConservation(k,q,initNuc,pionCharge,nucleon_out_charge,betaToCf,k_abs(3),low,high)
       if (writeFlag) write(12,*) k_abs(3),f(3)

       if (abs(f(3)).lt.maxError.and.k_abs(3).gt.0) then !Successful!!!
          get_k_abs=k_abs(3)
          success=.true.
          if (present(kout)) kout=k
          return
       end if

       if (counter.gt.maxCounter) then
          if (present(kout)) kout=k
          return
       end if
       counter=counter+1
       k_abs(1)=k_abs(2)
       k_abs(2)=k_abs(3)
       f(1)=f(2)
       f(2)=f(3)
    end do



  contains

!!$    subroutine scan()
!!$      real :: x,e
!!$      open(222,file="get_k_abs_scan.dat")
!!$      x=-1.
!!$      do
!!$         x=x+0.01
!!$         if (.not. set_k(x,k,initNuc,pionCharge,betaToCF)) then
!!$            write(*,*) 'k_abs Problem scan'
!!$         end if
!!$         e=error_energyConservation(k,q,initNuc,pionCharge,nucleon_out_charge,betaToCf,k_abs(3),low,high)
!!$
!!$         write(222,*) x,e
!!$         if(x.gt.10) exit
!!$      end do
!!$      close(222)
!!$    end subroutine scan

!!$    subroutine errorMessage
!!$      write(*,*)ak,bk,ck,q
!!$      write(*,*) 'k_abs(1:3):',k_abs
!!$      write(*,*)'k' , k
!!$    end subroutine errorMessage

    logical function set_k(absolut,k,initNuc,pionCharge,betaToCF)

      use energyCalc, only: energyDetermination
      use idTable, only: pion
      use particleDefinition
      use constants, only: mPi

      real               ,intent(in)  :: absolut
      integer            ,intent(in)  :: pionCharge
      type(particle)     ,intent(in)  :: initNuc
      real,dimension(0:3),intent(out) :: k
      real,dimension(1:3),intent(in)  :: betaToCF
      type(particle)                  :: outPion

      set_k = .true. ! default return value: success

      k(1:3) = (/ak,bk,ck/)*absolut

      if (abs(absolut) > 1000.0) then
         k(0)=sqrt(mPi**2+absolut**2)
         set_k = .false.
         return ! FAILURE !!!!
      end if

      if (pionpot) then
         ! Define outgoing pion
         outPion=initNuc
         outPion%ID=pion
         outPion%momentum(1:3)=k(1:3)
         outPion%charge=pionCharge
         outPion%mass=mPi
         call energyDetermination(outPion,betaToCF,warn=.false.)
         k(0)=outPion%momentum(0)
      else
         k(0)=sqrt(mPi**2+absolut**2)
      end if
    end function set_k

  end function get_k_abs




  !****************************************************************************
  !****s* electronPionProduction_kine/error_energyConservation
  ! NAME
  ! real function error_energyConservation(k, q, initial_nucleon, pionCharge,
  ! nucleon_out_charge, betaToCF, kabs, low_in, high_in)
  !
  ! PURPOSE
  ! Evaluates the errror in the energy conservation condition, after using
  ! momentum conservation :
  !   error= [E(pion)+E(final nucleon) ] - [ E(virtual photon)+E(initial nucleon) ]
  !
  ! INPUTS
  ! * real, dimension(0:3),intent(in) :: k
  ! * real, dimension(0:3),intent(in) :: q
  ! * type(particle),intent(in) :: initial_nucleon
  ! * integer, intent(in) :: pionCharge
  ! * integer, intent(in) :: nucleon_out_charge
  ! * real, dimension(1:3),intent(in) :: betaToCF
  ! * real, intent(in)   :: kabs
  ! * real, intent(in) :: low_in
  ! * real, intent(in) :: high_in
  !
  ! OUTPUT
  !
  !
  !****************************************************************************
  real function error_energyConservation(k, q, initial_nucleon, pionCharge,&
       & nucleon_out_charge,betaToCF,kabs,low_in,high_in)
    use particleDefinition
    use energyCalc, only: energyDetermination
    use idTable, only: nucleon
    use constants, only: mN

    type(particle),intent(in) :: initial_nucleon
    integer, intent(in) :: pionCharge,nucleon_out_charge
    real, intent(in)   :: kabs

    real, dimension(0:3),intent(in) :: k,q
    real, dimension(1:3),intent(in) :: betaToCF
    real :: penalty
    real, parameter :: penalty_constant=1.
    real, intent(in) :: low_in, high_in

    !real, dimension(0:3) :: pf
    type(particle) :: final_nucleon
    logical :: debugFlag=.false.

    if (debugFlag) write(*,*) 'in error_energyConservation'

    ! Construct the final state nucleon:
    final_nucleon%momentum(1:3)=q(1:3)+initial_nucleon%momentum(1:3)-k(1:3)
    final_nucleon%mass=mN
    final_nucleon%position(1:3)=initial_nucleon%position(1:3)
    final_nucleon%ID=nucleon
    final_nucleon%charge=nucleon_out_charge
    final_nucleon%perturbative=initial_nucleon%perturbative

!!$    If(AbsMom(final_nucleon).gt.1000) then
!!$       write(*,*) "Particle with large momentum"
!!$
!!$       write(*,*) 'q: ',q(1:3)
!!$       write(*,*) 'k: ',k(1:3)
!!$       write(*,*) 'N: ',initial_nucleon%momentum(1:3)
!!$    end If

    call energyDetermination(final_nucleon,betaToCF, warn=.false.)

    if (debugFlag) write(*,*) 'final_nucleon%momentum(0)=',final_nucleon%momentum(0)

    error_energyConservation=q(0)+initial_nucleon%momentum(0)-k(0)-final_nucleon%momentum(0)

    if (kabs.lt.low_in) then
       ! Make penalty since we are in a forbidden regime
       penalty=error_energyConservation/abs(error_energyConservation)*penalty_constant*abs(kabs-low_in)
       error_energyConservation=error_energyConservation+ penalty
    else if (kabs.gt.high_in) then
       penalty=error_energyConservation/abs(error_energyConservation)*penalty_constant*abs(kabs-high_in)
       error_energyConservation=error_energyConservation+ penalty
    end if

  end function error_energyConservation



  subroutine get_dV_Pi_dk(dV,k_abs,initNuc,pionCharge)
    ! returns the momentum derivative of the potential of the pion
    use potentialModule, only: potential_LRF
    use particleDefinition
    use idtable, only: pion

    real, intent(out) :: dV
    real, intent(in) :: k_abs
    type(particle),intent(in) :: initNuc
    integer, intent(in) :: pionCharge
    type(particle) :: outPion
    real, dimension(-2:2) :: V
    real,parameter :: dp=0.02
    integer :: i

    ! Define outgoing pion
    outPion=initNuc
    outPion%ID=pion
    outPion%momentum(1:3)=(/k_abs,0.,0./)
    outPion%charge=pionCharge
    ! Energy is not well defined since the potential should not depend on it:
    outPion%momentum(0)=0.

    do i=-2,2
       outPion%momentum(1:3)=(/k_abs+float(i)*dp,0.,0./)
       V(i)=potential_LRF(outPion)
    end do

    dV=(V(-2)-8.0*V(-1)+8.0*V(1)-V(2))  /(12.0*dp)
    if (debug) then
       write(*,*) 'dV_out (1st order)', (V(1)-V(-1))/(2.*dp)
       write(*,*) 'dV_out (2nd order)', dV
    end if

  end subroutine get_dV_Pi_dk



  subroutine getV_out(V_out,dV_out,pf,initNuc,pionCharge)
    ! returns the potential of the nucleon and its derivative at the kinematics of the outgoing nucleon
    use potentialModule, only:potential_LRF
    use particleDefinition

    real, intent(out) :: V_out, dV_out
    real, intent(in), dimension(0:3) :: pf
    type(particle),intent(in) :: initNuc
    integer, intent(in) :: pionCharge
    type(particle) :: outNuc
    real :: absMom2
    real, dimension(-2:2) :: V
    real,parameter :: dp=0.02
    integer :: i

    ! Define outgoing nucleon
    outNuc=initNuc
    outNuc%momentum=pf
    outNuc%charge=initNuc%charge-pionCharge
    absMom2=sqrt(Dot_product(pf(1:3),pf(1:3)))

    do i=-2,2
       outNuc%momentum(1:3)=(/0.,0.,float(i)*dp+absMom2/)
       V(i)=potential_LRF(outNuc)
    end do

    V_out=V(0)

    dV_out=(V(-2)-8.0*V(-1)+8.0*V(1)-V(2))/(12.0*dp)
    if (debug) then
       write(*,*) 'dV_out (1st order)', (V(1)-V(-1))/(2.*dp)
       write(*,*) 'dV_out (2nd order)', dV_out
    end if

  end subroutine getV_out


end module electronPionProduction_kine
