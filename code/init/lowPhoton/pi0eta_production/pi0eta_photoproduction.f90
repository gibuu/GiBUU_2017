!******************************************************************************
!****m* /pi0eta_photoproduction
! NAME
! module pi0eta_photoproduction
! PURPOSE
! Includes a parametrization for gamma p -> p pi0 eta cross sections.
!******************************************************************************
module pi0eta_photoproduction

  implicit none
  private

  logical, save                     :: initFlag=.true.
  real, save, dimension (1:100,1:2) :: sigma_field      ! Field to store input data
  integer, save                     :: max_index=0      ! Number of data points in input file
  real, save                        :: threshold=0.     ! Threshold for the process
  ! File name where the input data are stored
  character(100) ,parameter         :: partial_filename=trim("/photo_pi0eta/gamma_p_TO_p_pi0_eta.dat")

  public :: event_gamma_p_to_p_pi0_eta,  sigma_gamma_p_to_p_pi0_eta

contains


  !****************************************************************************
  !****s* pi0eta_photoproduction/event_gamma_p_to_p_pi0_eta
  ! NAME
  ! subroutine event_gamma_p_to_p_pi0_eta(nuc,gamma,sigTot,srts,srtFree,betaCM,mediumAtPosition,fState,flag,perturbative)
  ! PURPOSE
  ! Generates a pi^0 eta nucleon event.
  ! INPUTS
  ! * type(particle),intent(in):: nuc,gamma       ! incoming nucleon, photon
  ! * real,intent(in):: sigTot                    ! total cross section
  ! * real,intent(in):: srts                      ! sqrt(s)
  ! * real,intent(in):: srtFree                   ! sqrt(s) in vacuum
  ! * real,intent(in):: betaCM(3)                 ! beta for boost from CM to Lab frame
  ! * logical , intent(in ) :: perturbative       ! .true.= final state consists of perturbative particles
  ! * type(medium),intent(in) :: mediumAtPosition
  !
  ! OUTPUT
  ! * type(particle),intent(out):: fState(3)       ! final state particles
  ! * logical,intent(out):: flag                   ! true => success
  !****************************************************************************
  subroutine event_gamma_p_to_p_pi0_eta(nuc,gamma,sigTot,srts,srtFree,betaCM,mediumAtPosition,fState,flag,perturbative)
    use particleDefinition, only: particle
    use IdTable, only: nucleon,pion,eta
    use lorentzTrafo, only: lorentzCalcBeta
    use mediumDefinition, only: medium
    use dichteDefinition, only: dichte
    use densityModule, only: densityAt
    use master_2body, only: setKinematics
    use mediumModule, only: getMediumCutOff
    use RMF, only: getRMF_flag
    use constants, only: mN

    type(particle),intent(in):: nuc,gamma
    real,intent(in):: sigTot
    real,intent(in):: srts
    real,intent(in):: srtFree
    real,intent(in):: betaCM(3)
    logical , intent(in ) :: perturbative
    type(medium),intent(in) :: mediumAtPosition
    type(particle),intent(out):: fState(3)
    logical,intent(out):: flag

    type(particle):: iState(2)     ! initial state part. (gamma + N)
    real:: betaLRF(3)
    type(dichte) :: density

    flag=.false.
    ! set up initial state
    iState=(/gamma,nuc/)
    ! set up final state
    fState%perturbative=perturbative

    fState(1)%id=nucleon
    fState(2)%id=eta
    fState(3)%id=pion

    fState(1)  %charge=nuc%charge
    fState(2:3)%charge=0

    fState(1)%position=nuc%position
    fState(2)%position=nuc%position
    fState(3)%position=nuc%position

    fState(1)%mass=mN


    fState%perweight=sigTot
    ! Evaluate boost to LRF
    density=densityAt(nuc%position)
    if (density%baryon(0).gt.getMediumCutOff()/100. .and. .not.getRMF_flag() ) then
       betaLRF = lorentzCalcBeta (density%baryon, 'vecmesProduction')
    else
       betaLRF=0.
    end if
    ! assign the masses and kinematics
    call setKinematics(srts,srtFree,betaLRF,-betaCM,mediumAtPosition,iState,fState,flag)
    fState(2)%in_Formation=.true.
    fState(2)%formationTime=-999 ! use old formation time concept by default
  end subroutine event_gamma_p_to_p_pi0_eta



  !****************************************************************************
  !****f* pi0eta_photoproduction/sigma_gamma_p_to_p_pi0_eta
  ! NAME
  ! function sigma_gamma_p_to_p_pi0_eta(W) result(sigma)
  ! PURPOSE
  ! Returns a splined curve which goes through the data given in the file
  ! trim(path_to_input)//"/photo_pi0eta/gamma_p_TO_p_pi0_eta.dat"
  ! INPUTS
  ! real, intent(in) :: W -- Center of mass energy in vacuum
  !
  ! OUTPUT
  ! real :: sigma    -- Xsection in mu b
  !
  ! NOTES
  ! If the point W lies outside the input data, then we set sigma=0 for
  ! W<smallest energy and assume linear interpolation for W>W_max.
  !****************************************************************************
  function sigma_gamma_p_to_p_pi0_eta(W) result(sigma)
    use spline, only: Bsplint2
    use particleProperties, only: hadron
    use IDTABLE, only: eta
    use constants, only: mN, mPi

    real, intent(in) :: W
    real :: sigma

    if (initFlag) then
       call readData()
       initFlag=.false.
       threshold=hadron(eta)%mass+mPi+mN
    end if

    if (W.le.sigma_field(1,1).or.W.le.threshold) then
       sigma=0.
    else if (W.ge.sigma_field(max_Index,1)) then
       write(*,*) 'WARNING: sigma_gamma_p_to_p_pi0_eta out of bounds. W=',W
       sigma=max(0.,sigma_field(max_Index,2)+(sigma_field(max_Index,2)-sigma_field(max_Index-1,2))/&
            & (sigma_field(max_Index,1)-sigma_field(max_Index-1,1))*(W-sigma_field(max_Index,1)))
    else
       sigma = Bsplint2(sigma_field(1:max_Index,1),sigma_field(1:max_Index,2),W)
    end if
  end function sigma_gamma_p_to_p_pi0_eta


  !****************************************************************************
  !****s* pi0eta_photoproduction/readData
  ! NAME
  ! subroutine readData()
  ! PURPOSE
  ! Reads the data in the file
  ! trim(path_to_input)//"/photo_pi0eta/gamma_p_TO_p_pi0_eta.dat"
  ! to the field sigma_field. Max_index is set to the number of data points.
  !****************************************************************************
  subroutine readData()
    use inputGeneral, only: path_to_input

    integer :: ios
    integer :: i
    character(200) :: filename

    filename=trim(path_to_input)//trim(partial_filename)

    open(100,file=trim(filename),status='old',ioStat=ios)

    if (ios.ne.0) then
       write(*,*) 'ERROR in  sigma_gamma_p_to_p_pi0_eta'
       write(*,*) 'File', filename, " is not available"
       stop
    end if

    max_index=0
    do i=1,100
       Read(100,*,iostat=ios) sigma_field(i,1), sigma_field(i,2)
       if (ios.ne.0) then
          max_index=i-1
          exit
       end if
    end do
    if (max_index.eq.0) then
       write(*,*) 'ERROR in  sigma_gamma_p_to_p_pi0_eta, module pi0eta_photoproduction'
       write(*,*) 'Upper bound for sigma_field is too small!'
       write(*,*) 'You must increase this bound!!! STOP'
       stop
    end if

    write(*,*)
    write(*,*) 'Input for gamma p -> p pi^0 eta:'
    write(*,'(A,A14,A15)') '#', 'W[GeV]', 'sigma[mu b]'
    do i=1,max_index
       write(*,'(2F15.8)') sigma_Field(i,:)
    end do
    write(*,*) max_index
    write(*,*)

  end subroutine readData


end module pi0eta_photoproduction
