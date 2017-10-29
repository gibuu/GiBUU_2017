!******************************************************************************
!****m* /gamma2Pi_Xsections
! NAME
! module gamma2Pi_Xsections
! PURPOSE
! Includes the cross sections for gamma Nucleon -> Nucleon Pion Pion processes.
!******************************************************************************

module gamma2Pi_Xsections
  !****************************************************************************
  !****n* gamma2Pi_Xsections/gamma_2Pi_Xsections
  ! NAME
  ! NAMELIST gamma_2Pi_Xsections
  ! PURPOSE
  ! Includes:
  ! * experimentalXsections
  !****************************************************************************

  use inputGeneral
  use cl_splines

  implicit none
  private

  !****************************************************************************
  !****g* gamma2Pi_Xsections/experimentalXsections
  ! SOURCE
  !
  logical , save :: experimentalXsections=.true.
  !
  ! PURPOSE
  ! * If .true. then the Xsections are taken from the experiment
  ! * If .false. then the theoretical values are given
  !****************************************************************************

  logical , save :: initFlag=.true.

  public :: gamma2pi, cleanUp

  integer, save :: counter_CLerror=0

  type(tspline),save :: s1,s2,s3,s4,s5,s6

contains

  !****************************************************************************
  !****f* gamma2Pi_Xsections/init
  ! NAME
  ! subroutine init
  ! PURPOSE
  ! Read out job card to initialize module parameters.
  !****************************************************************************
  subroutine init

    use output

    integer :: ios

    NAMELIST /gamma_2Pi_Xsections/ experimentalXsections

    call Write_ReadingInput('gamma_2Pi_Xsections',0)

    rewind(5)
    read(5,nml=gamma_2Pi_xsections,IOSTAT=IOS)
    call Write_ReadingInput('gamma_2Pi_Xsections',0,IOS)

    write(*,*) 'Set experimentalXsections       :' ,experimentalXsections
    call Write_ReadingInput('gamma_2Pi_Xsections',1)

  end subroutine init


  subroutine cleanUp
    call cl_cleanupSpline(s1)
    call cl_cleanupSpline(s2)
    call cl_cleanupSpline(s3)
    call cl_cleanupSpline(s4)
    call cl_cleanupSpline(s5)
    call cl_cleanupSpline(s6)
  end subroutine


  !****************************************************************************
  !****s* gamma2Pi_Xsections/gamma2pi
  ! NAME
  ! subroutine gamma2pi(qnuk,srts,sig2pi,betaToLRF,mediumAtPosition,position)
  !
  ! PURPOSE
  ! Calculation of 2 pion cross sections
  ! INPUTS
  ! * integer :: qnuk      -- charge of nucleon
  ! * real    :: srts      -- SQRT(s)
  ! OUTPUT
  ! * sig2pi(0:3) Cross sections for gamma+nucleon->nucleon+2Pi production
  ! * sig2pi(1) -> nucleon piMinus piPlus
  ! * sig2pi(2) -> nucleon piPlus piNull or nucleon piMinus piNull
  ! * sig2pi(3) -> nucleon piNull piNull
  ! * sig2pi(0) -> Total Xsection into nucleon+2 Pions
  !
  ! NOTES
  ! * Threshold improvement by Pascal Muehlich
  ! * Units of cross sections : [10^-6 barn]
  !****************************************************************************

  subroutine gamma2pi(qnuk,srts,sig2pi,betaToLRF,mediumAtPosition,position)

    use gamma2Pi_Xsections_luis, only: minva
    use mediumDefinition
    use constants, only: mN, mPi

    integer, intent(in)              :: qnuk
    real, intent(in)                 :: srts
    type(medium) ,intent(in) :: mediumAtPosition
    real, dimension(1:3) ,intent(in) :: betaToLRF
    real, dimension(1:3) ,intent(in) :: position

    real, dimension(0:3),intent(out)  :: sig2pi


    real :: photonEnergy
    real :: pLab!,resmass
    real, dimension(1:3) :: sigma
    logical, dimension(1:3) :: flags
    real :: sigma_doublecharged,sigma_singlecharged,sigma_Uncharged
    logical , parameter :: improvedThreshold=.true.
!   real, external :: gP_ProtonPiPlusPiMinus , gP_NeutronPiPlusPiNull
!   real, external :: gP_ProtonPiNullPiNull
!   real, external :: gN_NeutronPiPlusPiMinus
!   real, external :: gN_ProtonPiMinusPiNull
!   real, external :: gN_NeutronPiNullPiNull


    if (initFlag) then
       call init
       initFlag=.false.
    end if

    sig2pi=0.

    if (srts.gt.mN+2.*mPi) then
       plab=(srts**2-mN**2)/2./mN
       ! Experimental 2 pion cross sections:
       ! Double charged : two charged pions
       ! Single charged : one charged pion
       ! uncharged      : only pi^0 's

       sigma_doubleCharged = 0.
       sigma_singlecharged = 0.
       sigma_unCharged     = 0.

       if (experimentalXsections) then
          if (qnuk.eq.1) then
             sigma_doubleCharged = gP_ProtonPiPlusPiMinus (plab)
             sigma_singlecharged = gP_NeutronPiPlusPiNull (plab)
             sigma_unCharged     = gP_ProtonPiNullPiNull  (plab)
          else if (qnuk.eq.0) then
             sigma_doublecharged= gN_NeutronPiPlusPiMinus (plab)
             sigma_singlecharged= gN_ProtonPiMinusPiNull  (plab)
             sigma_unCharged=     gN_NeutronPiNullPiNull  (plab)
          else
             write(*,*) 'Charge not valid in gamma2pi' , qnuk
          end if

          if (improvedThreshold) then
             !*     improved threshold behaviour  (take care of region where there are no data!)
             call thres2pi(srts,sigma,qnuk,flags)
             if (sigma(1).gt.0.and.flags(1)) sigma_DoubleCharged=sigma(1)
             if (sigma(2).gt.0.and.flags(2)) sigma_SingleCharged=sigma(2)
             if (sigma(3).gt.0.and.flags(3)) sigma_UnCharged    =sigma(3)
          end if
       else
          photonEnergy=(srts**2-mN**2)/2./mN
          ! Luis Xsextions by Amplitudes
          if (qnuk.eq.1) then
             !1: gamma p -> pi+ pi- p
             !2: gamma p -> pi+ pi0 n
             !3: gamma p -> pi0 pi0 p
             sigma_doubleCharged = minva(photonEnergy,1,betaToLRF,mediumAtPosition,position)
             sigma_singlecharged = minva(photonEnergy,2,betaToLRF,mediumAtPosition,position)
             sigma_unCharged     = minva(photonEnergy,3,betaToLRF,mediumAtPosition,position)
          else if (qnuk.eq.0) then
             !4: gamma n -> pi+ pi- n
             !5: gamma n -> pi- pi0 p
             !6: gamma n -> pi0 pi0 n
             sigma_doublecharged= minva(photonEnergy,4,betaToLRF,mediumAtPosition,position)
             sigma_singlecharged= minva(photonEnergy,5,betaToLRF,mediumAtPosition,position)
             sigma_unCharged=     minva(photonEnergy,6,betaToLRF,mediumAtPosition,position)
          else
             write(*,*) 'Charge not valid in gamma2pi' , qnuk
          end if

       end if

       sig2pi(1)=max(sigma_doublecharged,0.) ! piMinus piPlus
       sig2pi(2)=max(sigma_singlecharged,0.) ! piPlus piNull or piMinus piNull
       sig2pi(3)=max(sigma_Uncharged,0.)      ! piNull piNull
       sig2pi(0)=Sum(sig2pi(1:3))    ! Total Xsection into nucleon+2 Pions
    end if

    return

  end subroutine gamma2pi


  !****************************************************************************
  !*2 pion off the proton
  !****************************************************************************

  function gP_ProtonPiPlusPiMinus(plab,minPlab,minSig) Result (sigma)
    !****f* gamma2Pi_Xsections/gP_ProtonPiPlusPiMinus
    ! Subroutine for calculation of gamma p -> p pi+ pi- cross section
    ! NOTES
    ! Reads out data file and uses a spline with it.
    ! OUTPUT
    ! * real, optional     :: minPlab,minSig ! plab and sigma of lowest lying data point in the data file
    !**************************************************************************

    real, intent(in) :: plab
    real :: sigma
    integer, parameter :: Dateiende=-1

    real, dimension(1:200),save :: pLabField, sigmaField!, derivativeField

    logical, save :: initFlag=.true.
    integer :: ios
    real,save :: maximalMomentum,maximalCrossSection
    integer, save :: i
    logical, parameter :: debug=.false.
    logical            :: success
    integer            :: error
    real, optional,intent(out)     :: minPlab,minSig

    if (initFlag) then
       if (debug) write(*,*) 'Initializing gP_ProtonPiPlusPiMinus'
       open(100,file=trim(path_to_Input)//'/photo_twoPi/gamp-ppipm.dat',status='old')

       do i=lbound(plabField,dim=1),ubound(plabField,dim=1)
          read(100,*,IOSTAT=IOS) plabField(i),sigmaField(i)
          if (IOS.eq.DATEIENDE) exit
       end do
       i=i-1
       maximalMomentum=pLabField(i)
       maximalCrossSection=sigmaField(i)
       if (debug) write(*,*) 'Finished Initializing gP_ProtonPiPlusPiMinus',maximalMomentum, maximalCrossSection

       close(100)

       s1=cl_initSpline(plabField(lbound(plabField,dim=1):i),sigmaField(lBound(plabField,dim=1):i))
       initFlag=.false.
    end if


    sigma=cl_spline(s1,plab,success,error)
    if (.not.success.and.error.gt.0.and.counter_CLerror.lt.100) then
       call cl_error(error,' gP_ProtonPiPlusPiMinus',plab)
       counter_CLerror=counter_CLerror+1
       if (counter_clerror.eq.100) write(*,*) 'Now I STOP outputting this error!!!!'
    end if


    if (plab.gt.maximalMomentum) then
       !write(*,*) 'Error: Momentum out of bounds in gP_ProtonPiPlusPiMinus'
       sigma=maximalCrossSection
    end if
    ! Return lowest lying data point
    if (present(minPlab).and.present(minSig)) then
       minPlab=s1%x(lbound(s1%x,dim=1))
       minSig =s1%y(lbound(s1%y,dim=1),1)
    end if
  end function gP_ProtonPiPlusPiMinus

  !****************************************************************************

  function gP_NeutronPiPlusPiNull(plab,minPlab,minSig) Result (sigma)
    !****f* gamma2Pi_Xsections/gP_NeutronPiPlusPiNull
    ! Subroutine for calculation of gamma p -> n pi+ pi0  cross section
    ! NOTES
    ! Reads out data file and uses a spline with it.
    ! OUTPUT
    ! * real, optional     :: minPlab,minSig ! plab and sigma of lowest lying data point in the data file
    !**************************************************************************

    real, intent(in) :: plab
    real  :: sigma
    integer, parameter :: Dateiende=-1

    real, dimension(1:200),save :: pLabField, sigmaField!, derivativeField

    logical, save :: initFlag=.true.
    integer ::  ios
    real,save :: maximalMomentum,maximalCrossSection
    integer, save :: i
    logical            :: success
    integer            :: error
    real, optional     :: minPlab,minSig


    if (initFlag) then
       open(100,file=trim(path_to_Input)//'/photo_twoPi/gamp-npip0.dat',status='old')

       do i=lbound(plabField,dim=1),ubound(plabField,dim=1)
          read(100,*,IOSTAT=IOS) plabField(i),sigmaField(i)
          if (IOS.eq.DATEIENDE)      exit
       end do
       i=i-1
       maximalMomentum=pLabField(i)
       maximalCrossSection=sigmaField(i)

       close(100)
       s2=cl_initSpline(plabField(lbound(plabField,dim=1):i),sigmaField(lBound(plabField,dim=1):i))

       initFlag=.false.
    end if


    sigma=cl_spline(s2,plab,success,error)

    if (.not.success.and.error.gt.0.and.counter_CLerror.lt.100) then
       call cl_error(error,' gP_NeutronPiPlusPiNull',plab)
       counter_CLerror=counter_CLerror+1
       if (counter_clerror.eq.100) write(*,*) 'Now I STOP outputting this error!!!!'
    end if


    if (plab.gt.maximalMomentum) then
       !write(*,*) 'Error: Momentum out of bounds in gP_NeutronPiPlusPiNull'
       sigma=maximalCrossSection
    end if
    ! Return lowest lying data point
    if (present(minPlab).and.present(minSig)) then
       minPlab=s2%x(lbound(s2%x,dim=1))
       minSig =s2%y(lbound(s2%y,dim=1),1)
    end if
  end function gP_NeutronPiPlusPiNull

  !****************************************************************************

  function gP_ProtonPiNullPiNull(plab,minPlab,minSig) Result (sigma)
    !****f* gamma2Pi_Xsections/gP_ProtonPiNullPiNul
    ! Subroutine for calculation of gamma p -> p pi0 pi0  cross section
    ! NOTES
    ! Reads out data file and uses a spline with it.
    ! OUTPUT
    ! * real, optional     :: minPlab,minSig ! plab and sigma of lowest lying data point in the data file
    !**************************************************************************

    real, intent(in) :: plab
    real  :: sigma
    integer, parameter :: Dateiende=-1

    real, dimension(1:200),save :: pLabField, sigmaField!, derivativeField

    logical, save :: initFlag=.true.
    integer ::  ios
    real,save :: maximalMomentum,maximalCrossSection
    integer, save :: i
    logical            :: success
    integer            :: error
    real, optional     :: minPlab,minSig


    if (initFlag) then
       open(100,file=trim(path_to_Input)//'/photo_twoPi/gamp-ppi00.dat',status='old')


       do i=lbound(plabField,dim=1),ubound(plabField,dim=1)
          read(100,*,IOSTAT=IOS) plabField(i),sigmaField(i)
          if (IOS.eq.DATEIENDE)      exit
       end do
       i=i-1
       maximalMomentum=pLabField(i)
       maximalCrossSection=sigmaField(i)

       close(100)
       s3=cl_initSpline(plabField(lbound(plabField,dim=1):i),sigmaField(lBound(plabField,dim=1):i))

       initFlag=.false.
    end if


    sigma=cl_spline(s3,plab,success,error)

    if (.not.success.and.error.gt.0.and.counter_CLerror.lt.100) then
       call cl_error(error,' gP_ProtonPiNullPiNull',plab)
       counter_CLerror=counter_CLerror+1
       if (counter_clerror.eq.100) write(*,*) 'Now I STOP outputting this error!!!!'
    end if

    if (plab.gt.maximalMomentum) then
       !write(*,*) 'Error: Momentum out of bounds in gP_ProtonPiNullPiNull'
       sigma=maximalCrossSection
    end if
    ! Return lowest lying data point
    if (present(minPlab).and.present(minSig)) then
       minPlab=s3%x(lbound(s3%x,dim=1))
       minSig =s3%y(lbound(s3%y,dim=1),1)
    end if
  end function gP_ProtonPiNullPiNull

  !****************************************************************************
  !*2 pion off the neutron
  !****************************************************************************

  function gN_NeutronPiPlusPiMinus(plab,minPlab,minSig) Result (sigma)
    !****f* gamma2Pi_Xsections/gN_NeutronPiPlusPiMinus
    ! Subroutine for calculation of gamma n -> n pi+ pi-  cross section
    ! NOTES
    ! Reads out data file and uses a spline with it.
    ! OUTPUT
    ! * real, optional     :: minPlab,minSig ! plab and sigma of lowest lying data point in the data file
    !**************************************************************************

    real, intent(in) :: plab
    real  :: sigma
    integer, parameter :: Dateiende=-1

    real, dimension(1:200),save :: pLabField, sigmaField!, derivativeField

    logical, save :: initFlag=.true.
    integer :: ios
    real,save :: maximalMomentum,maximalCrossSection
    integer, save :: i
    logical            :: success
    integer            :: error
    real, optional     :: minPlab,minSig


    if (initFlag) then
       open(100,file=trim(path_to_Input)//'/photo_twoPi/gamn-npipm.dat',status='old')

       do i=lbound(plabField,dim=1),ubound(plabField,dim=1)
          read(100,*,IOSTAT=IOS) plabField(i),sigmaField(i)
          if (IOS.eq.DATEIENDE)   exit
       end do
       i=i-1
       maximalMomentum=pLabField(i)
       maximalCrossSection=sigmaField(i)

       close(100)
       s4=cl_initSpline(plabField(lbound(plabField,dim=1):i),sigmaField(lBound(plabField,dim=1):i))
       initFlag=.false.
    end if

    sigma=cl_spline(s4,plab,success,error)


    if (.not.success.and.error.gt.0.and.counter_CLerror.lt.100) then
       call cl_error(error,' gN_NeutronPiPlusPiMinus',plab)
       counter_CLerror=counter_CLerror+1
       if (counter_clerror.eq.100) write(*,*) 'Now I STOP outputting this error!!!!'
    end if

    if (plab.gt.maximalMomentum) then
       !write(*,*) 'Error: Momentum out of bounds in gN_NeutronPiPlusPiMinus'
       sigma=maximalCrossSection
    end if
    ! Return lowest lying data point
    if (present(minPlab).and.present(minSig)) then
       minPlab=s4%x(lbound(s4%x,dim=1))
       minSig =s4%y(lbound(s4%y,dim=1),1)
    end if
  end function gN_NeutronPiPlusPiMinus

  !****************************************************************************

  function gN_ProtonPiMinusPiNull(plab,minPlab,minSig) Result (sigma)
    !****f* gamma2Pi_Xsections/gN_ProtonPiMinusPiNull
    ! Subroutine for calculation of gamma n -> p pi- pi0  cross section
    ! NOTES
    ! Reads out data file and uses a spline with it.
    ! OUTPUT
    ! * real, optional     :: minPlab,minSig ! plab and sigma of lowest lying data point in the data file
    !**************************************************************************

    real, intent(in) :: plab
    real  :: sigma
    integer, parameter :: Dateiende=-1

    real, dimension(1:200),save :: pLabField, sigmaField!, derivativeField

    logical, save :: initFlag=.true.
    integer ::  ios
    real,save :: maximalMomentum,maximalCrossSection
    integer, save :: i
    logical            :: success
    integer            :: error
    real, optional     :: minPlab,minSig


    if (initFlag) then
       open(100,file=trim(path_to_Input)//'/photo_twoPi/gamn-ppim0.dat',status='old')

       do i=lbound(plabField,dim=1),ubound(plabField,dim=1)
          read(100,*,IOSTAT=IOS) plabField(i),sigmaField(i)
          if (IOS.eq.DATEIENDE)    exit
       end do
       i=i-1
       maximalMomentum=pLabField(i)
       maximalCrossSection=sigmaField(i)

       close(100)
       s5=cl_initSpline(plabField(lbound(plabField,dim=1):i),sigmaField(lBound(plabField,dim=1):i))
       initFlag=.false.
    end if


    sigma=cl_spline(s5,plab,success,error)

    if (.not.success.and.error.gt.0.and.counter_CLerror.lt.100) then
       call cl_error(error,' gN_ProtonPiMinusPiNull',plab)
       counter_CLerror=counter_CLerror+1
       if (counter_clerror.eq.100) write(*,*) 'Now I STOP outputting this error!!!!'
    end if

    if (plab.gt.maximalMomentum) then
       !write(*,*) 'Error: Momentum out of bounds in gN_ProtonPiMinusPiNull'
       sigma=maximalCrossSection
    end if
    ! Return lowest lying data point
    if (present(minPlab).and.present(minSig)) then
       minPlab=s5%x(lbound(s5%x,dim=1))
       minSig =s5%y(lbound(s5%y,dim=1),1)
    end if
  end function gN_ProtonPiMinusPiNull

  !****************************************************************************


  function gN_NeutronPiNullPiNull(plab,minPlab,minSig) Result (sigma)
    !****f* gamma2Pi_Xsections/gN_NeutronPiNullPiNull
    ! Subroutine for calculation of gamma n -> n pi0 pi0  cross section
    ! NOTES
    ! Reads out data file and uses a spline with it.
    ! OUTPUT
    ! * real, optional     :: minPlab,minSig ! plab and sigma of lowest lying data point in the data file
    !**************************************************************************

    real, intent(in) :: plab
    real  :: sigma
    integer, parameter :: Dateiende=-1

    real, dimension(1:200),save :: pLabField, sigmaField!, derivativeField

    logical, save :: initFlag=.true.
    integer :: ios
    real,save :: maximalMomentum,maximalCrossSection
    integer, save :: i
    logical            :: success
    integer            :: error
    real, optional     :: minPlab,minSig


    if (initFlag) then
       open(100,file=trim(path_to_Input)//'/photo_twoPi/gamn-npi00.dat',status='old')

       do i=lbound(plabField,dim=1),ubound(plabField,dim=1)
          read(100,*,IOSTAT=IOS) plabField(i),sigmaField(i)
          if (IOS.eq.DATEIENDE)   exit
       end do
       i=i-1
       maximalMomentum=pLabField(i)
       maximalCrossSection=sigmaField(i)

       close(100)
       s6=cl_initSpline(plabField(lbound(plabField,dim=1):i),sigmaField(lBound(plabField,dim=1):i))
       initFlag=.false.
    end if

    sigma=cl_spline(s6,plab,success,error)
    if (.not.success.and.error.gt.0.and.counter_CLerror.lt.100) then
       call cl_error(error,' gN_NeutronPiNullPiNull',plab)
       counter_CLerror=counter_CLerror+1
       if (counter_clerror.eq.100) write(*,*) 'Now I STOP outputting this error!!!!'
    end if

    if (plab.gt.maximalMomentum) then
       !write(*,*) 'Error: Momentum out of bounds in gN_ProtonPiMinusPiNull'
       sigma=maximalCrossSection
    end if
    ! Return lowest lying data point
    if (present(minPlab).and.present(minSig)) then
       minPlab=s6%x(lbound(s6%x,dim=1))
       minSig =s6%y(lbound(s6%y,dim=1),1)
    end if
  end function gN_NeutronPiNullPiNull


  !****************************************************************************
  !****s* gamma2Pi_Xsections/thres2pi
  ! NAME
  ! subroutine thres2pi(srts,sigma,qnuk,inThresholdRegion)
  ! NOTES
  ! * This subroutine provides a reasonable threshold behaviour for
  !   the two pion photoproduction cross sections. gamma nucleon -> nucleon+2pi
  !
  ! * Retrieves the lowest data point of all the parametrizations and fits a phase space curve to it.
  !
  ! AUTHOR
  ! * Oliver Buss
  !
  ! INPUTS
  ! * qnuK : charge of nucleon
  ! * srts : sqrt(s)
  !
  ! OUTPUT
  ! * real :: sigma(1:3)
  !   cross sections (1=double charged, 2= single charged, 3= pi^0 pi^0)
  !
  ! * logical, dimension(1:3), intent (out) :: inThresholdRegion
  !   .true. if sqrt(s) is smaller than lowest data point in the channel (1=double charged, 2= single charged, 3= pi^0 pi^0)
  !****************************************************************************
  subroutine thres2pi(srts,sigma,qnuk,inThresholdRegion)

    use constants, only: mPi, mN
    use threeBodyPhaseSpace, only: Integrate_3bodyPS
    use output, only: Write_InitStatus

    real   , intent(in) :: srts
    integer, intent(in) :: qnuk
    real, dimension(1:3)   , intent (out) :: sigma
    logical, dimension(1:3), intent (out) :: inThresholdRegion

    real :: flux,ener
    real, dimension(0:1,1:3), save :: matrixSquared=0.
    real, dimension(0:1,1:3), save :: min_egamma
    real, dimension(0:1,1:3), save :: min_sigma
    logical, save :: firstTime=.true.
    integer :: i

    if (firstTime) then
       call evaluate_msquared
       firsttime=.false.
    end if

    sigma=0.0
    ener=(srts**2-mn**2)/2./mn
    flux=4./(mn*ener)
    sigma = matrixSquared(qnuk,:) * flux * Integrate_3bodyPS(srts,mN,mPi,mPi)
    do i=1,3
       inThresholdRegion(i)=(ener.lt.min_egamma(qnuk,i))
    end do

  contains
    subroutine evaluate_MSquared
        ! Evaluates the (Matrix element squared) by fitting to the data and extracts the minima of the data points

        real :: dummy,plab_dummy
        integer :: i,j
        real :: srts
        real :: flux
        plab_dummy=0.5
        dummy = gP_ProtonPiPlusPiMinus (plab_dummy,min_egamma(1,1),min_sigma(1,1))
        dummy = gP_NeutronPiPlusPiNull (plab_dummy,min_egamma(1,2),min_sigma(1,2))
        dummy = gP_ProtonPiNullPiNull  (plab_dummy,min_egamma(1,3),min_sigma(1,3))
        dummy = gN_NeutronPiPlusPiMinus (plab_dummy,min_egamma(0,1),min_sigma(0,1))
        dummy = gN_ProtonPiMinusPiNull  (plab_dummy,min_egamma(0,2),min_sigma(0,2))
        dummy = gN_NeutronPiNullPiNull  (plab_dummy,min_egamma(0,3),min_sigma(0,3))

        do i=0,1
           do j=1,3
              flux=4./(mn*min_egamma(i,j))
              srts=sqrt(mn**2+2.*mn*min_egamma(i,j))
              matrixSquared(i,j) = min_sigma(i,j) / flux / Integrate_3bodyPS(srts,mN,mPi,mPi)
           end do
        end do

        call Write_InitStatus('2Pi threshold',0)
        write(*,*) '     Minimal data points:    Egamma          sigma          |M|^2'
        write(*,'(A,f10.4,1p,2E15.4)') ' gamma P -> P pi^- pi^+  :', min_egamma(1,1),min_sigma(1,1),matrixSquared(1,1)
        write(*,'(A,f10.4,1p,2E15.4)') ' gamma P -> N pi^+ pi^0  :', min_egamma(1,2),min_sigma(1,2),matrixSquared(1,2)
        write(*,'(A,f10.4,1p,2E15.4)') ' gamma P -> P pi^0 pi^0  :', min_egamma(1,3),min_sigma(1,3),matrixSquared(1,3)
        write(*,'(A,f10.4,1p,2E15.4)') ' gamma N -> N pi^+ pi^-  :', min_egamma(0,1),min_sigma(0,1),matrixSquared(0,1)
        write(*,'(A,f10.4,1p,2E15.4)') ' gamma N -> P pi^- pi^0  :', min_egamma(0,2),min_sigma(0,2),matrixSquared(0,2)
        write(*,'(A,f10.4,1p,2E15.4)') ' gamma N -> N pi^0 pi^0  :', min_egamma(0,3),min_sigma(0,3),matrixSquared(0,3)
        call Write_InitStatus('2Pi threshold',1)

      end subroutine evaluate_MSquared
  end subroutine thres2pi
end module gamma2Pi_Xsections
