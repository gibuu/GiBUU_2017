!******************************************************************************
!****m* /pionXsection
! NAME
! pionXsection
! PURPOSE
! Module with all routines necessary to make output for pion induced processes.
! Valid up to 180 MeV incoming pion energy (above that 2pi production becomes
! important).Therefore there could be a pipi finalstate, which could lead to
! double counting.
!******************************************************************************
module pionXsection

  implicit none
  private

  Public:: pionXsectionAnalysis

  !****************************************************************************
  !****g* pionXsection/debug
  ! SOURCE
  logical, parameter :: debug=.false.
  ! PURPOSE
  ! If .true. debug infos are produced.
  !****************************************************************************

  !****************************************************************************
  !****g* pionXsection/CMFrame
  ! SOURCE
  logical, save :: CMFrame = .false.
  ! PURPOSE
  ! If .true. Xsection is evaluated in CM-Frame of the incoming pion and a
  ! resting nucleon, else in calculation frame.
  !****************************************************************************

  logical, save :: initflag=.true.

  !****************************************************************************
  !****g* pionXsection/dsigma_dOmegadE_switch
  ! SOURCE
  logical, save :: dsigma_dOmegadE_switch=.false.
  ! PURPOSE
  ! If .true. then dsigma/dOmega and dSigma/dOmega/dE are evaluated.
  !****************************************************************************


  !****************************************************************************
  !****g* pionXsection/twoPi_switch
  ! SOURCE
  logical, save :: twoPi_switch=.false.
  ! PURPOSE
  ! If .true. then 2Pi output is evaluated.
  !****************************************************************************


contains


  !****************************************************************************
  !****s* pionXsection/readinput
  ! NAME
  ! subroutine readinput
  ! INPUTS
  ! NONE
  ! OUTPUT
  ! NONE
  ! PURPOSE
  ! This subroutine reads the namelist "pionAnalysis" out of the jobcard.
  ! Only called once to initialize the module.
  !****************************************************************************
  subroutine readinput
    use output, only: Write_ReadingInput
    integer :: IOS

    !**************************************************************************
    !****n* pionXsection/pionAnalysis
    ! NAME
    ! NAMELIST /pionAnalysis/
    ! PURPOSE
    ! Includes the input switches:
    ! * CMFrame
    ! * dsigma_dOmegadE_switch
    ! * twoPi_switch
    !**************************************************************************
    NAMELIST /pionAnalysis/ CMFrame, dsigma_dOmegadE_switch, twoPi_switch

    call Write_ReadingInput('pionAnalysis',0)
    rewind(5)
    read(5,nml=pionAnalysis,IOSTAT=IOS)
    call Write_ReadingInput('pionAnalysis',0,IOS)
    write(*,*) '  Evaluating dsigma/dOmega and dSigma/dOmega/dE ? ', dsigma_dOmegadE_switch
    write(*,*) '  Evaluating 2 Pion output                      ? ', twoPi_switch

    write(*,*) '  Switch for CM-Frame :',CMFrame
    call Write_ReadingInput('pionAnalysis',1)

  end subroutine readinput


  !****************************************************************************
  !****s* pionXsection/pionXsectionAnalysis
  ! NAME
  ! subroutine pionXsectionAnalysis(pertParticles,finalFlag)
  ! INPUTS
  ! type(particle), intent(in),dimension(:,:)  :: pertParticles        ! Particles which shall be analyzed
  ! logical, intent(in) :: finalFlag                                   ! if .true. than the final output for a series of calls will be done
  ! OUTPUT
  !  * File './pionInduced_xSections.dat'   with absorption, reaction and quasielastic cross sections
  !  * File 'pionInduced_dTheta'//RealToChar(getEkin()*1000.) //'MeV.dat' with quasielastic cross section : dsigma/dtheta in mB/radian
  !
  ! NOTES
  ! This routine is based on "pionout" in inpion.f by M.Effenberger
  ! PURPOSE
  ! * This subroutine produces output for pion-nucleus scattering, including impact
  ! Parameter integration.
  ! * Gives dSigma/dOmega for Pions in file "pionInduced_dTheta... .dat"
  ! * Gives Sigma_reac und Sigma_abs  for Pions in file "pionInduced_Xsections.dat"
  ! USES
  ! To generate 2pi Xsection output we use the routine:
  !  * use LowPhotonAnalysis, only : twoPi_Output
  !****************************************************************************
  subroutine pionXsectionAnalysis  (pertParticles,finalFlag)
    use LowPhotonAnalysis, only: twoPi_Output
    use particleDefinition
    use idTable
    use history, only: history_getGeneration
    use output
    use constants, only: pi, mN, mPi
    use initPion, only: getEkin, getTotalPerweight
    use lorentzTrafo, only: lorentz
    use collisionNumbering, only:pert_numbering
    use inputGeneral, only:localEnsemble

    type(particle), intent(in),dimension(:,:)  :: pertParticles
    logical, intent(in) :: finalFlag

    ! Local
    real, save,dimension(-1:1,0:10)  ::  quasiElastic_Xsection, quasiElastic_error,quasielastic_total_error
    real, dimension(-1:1) :: sigma_QE, sigma_QE_error
    integer :: generation
    real, save :: absorption_xSection, sigmaTotalSave
    integer, save :: numberRuns=0
    integer, save :: absEvents
    real :: totalPerweight
    real :: sigma_Absorption, sigmaTotal
    real,save :: sigma_Absorption_error, sigma_Total_error
    real, dimension (-1:1,0:10) :: sigma_quasi
    integer :: index, ensemble
    real :: sigTot!,reaction
    logical, save :: startFlag=.true.
    logical, save :: firstFlag=.true.

    real :: theta ! Angle to z-Axis
    real, parameter :: deltaTheta=1. ! Delta(Theta)
    real, dimension (-1:1,0:180),save :: sigma_dTheta  ! first index : charge, second index angle
    integer :: thetaIndex
    integer :: i
    real :: absMom_perpendicular
    character(60) :: fileName
    logical, save :: firstFlag_all=.true.
    real, dimension(0:3) :: mom
    real, dimension(1:3) :: beta
    !logical,save :: writeHeader=.true.

    if (initflag) then
       call readinput
       initflag=.false.
    end if

    totalPerweight=getTotalPerweight()

    write(*,*) '** In pionXsectionAnalysis'
    write(*,*) ' Total Perweight=' , totalPerweight

    ! **********initialization******************************

    if (startFlag) then
       write(*,*) '  Initializing pion analysis'
       absorption_xSection=0.
       sigmaTotalSave=0.
       numberRuns=0
       do i=-1,1
          quasiElastic_xSection(i,:)=0.
          quasiElastic_error(i,:)=0.
       end do
       absEvents=0
       sigma_dTheta(-1,:)=0.
       sigma_dTheta( 0,:)=0.
       sigma_dTheta( 1,:)=0.
       sigma_absorption_error=0.
       sigma_Total_error=0.
    end if


    numberRuns=numberRuns+1

    !******************absorption cross section*******************************

    !***absorption heisst: Es kommt kein Pion aus der Reaktion***********
    !** Check ob Pion in Kern gefangen wird: Ekin > Epot!!**************

    ! (1) Evaluate absorption cross section of specific run
    sigma_Absorption=totalPerweight ! black disc = no pion escaped
    sigmaTotal=totalPerweight ! black disc = no pion escaped

    ensembleLoop : do ensemble=1,size(pertParticles,dim=1)
       indexLoop : do index=1,size(pertParticles,dim=2)
          ! Count pions that excaped and reduce absorption cross section due to them
          if (pertParticles(ensemble,index)%ID.eq.pion) then
             sigma_Absorption=sigma_Absorption-pertParticles(ensemble,index)%perweight
             absEvents=absEvents+1
             ! If(pertParticles(ensemble,index)%lastCollisionTime.le.0.0001) then ! particle did not interact
             if (pertParticles(ensemble,index)%event(1).eq.pert_numbering()) then ! particle did not interact
                sigmaTotal=sigmaTotal-pertParticles(ensemble,index)%perweight
             end if
          end if
       end do indexLoop
    end do ensembleLoop



    ! Write output of every single run
    if (firstFlag_all) then
       open(345,File='./pionInduced_xSections_all.dat')
       firstFlag_all=.false.
       write(345,*) '# elab, Sigma  abs , Sigma_total'
    else
       open(345,File='./pionInduced_xSections_all.dat',position='append')
    end if
    write(345,'(8F10.4 ,2I12)') getEKin(),sigma_Absorption, sigmaTotal
    close(345)

    ! FOR ERROR ANALYSIS, used for standard deviations among single runs
    sigma_total_error=sigma_total_error+sigmaTotal**2
    sigma_absorption_error=sigma_absorption_error+sigma_Absorption**2


    ! (2) Add up cross sections of all runs
    absorption_xSection=absorption_XSection+sigma_Absorption
    sigmaTotalsave=sigmaTotalSave+sigmaTotal


    !**Quasielastic Crosssection********************************************

    ! (1) Evaluate cross section of some specific run
    do i=-1,1
       sigma_Quasi(i,:)=0
    end do
    ensembleSchleife : do ensemble=1,size(pertParticles,dim=1)
       indexSchleife : do index=1,size(pertParticles,dim=2)
          if (pertParticles(ensemble,index)%ID <= 0) cycle indexSchleife
          ! Count pions that escaped and look wether they 've been undergoing collisions
          !          if((pertParticles(ensemble,index)%ID.eq.pion).and.(pertParticles(ensemble,index)%lastCollisionTime.ge.0.0001)) then
          if ((pertParticles(ensemble,index)%ID.eq.pion).and.(pertParticles(ensemble,index)%event(1) &
               & .ne.pert_numbering())) then
             generation=min(ubound(sigma_Quasi,dim=2),history_getGeneration(pertParticles(ensemble,index)%history))
             sigma_Quasi(pertParticles(ensemble,index)%charge,generation)=&
                  & sigma_Quasi(pertParticles(ensemble,index)%charge,generation)&
                  & +pertParticles(ensemble,index)%perweight
          end if

          if (pertParticles(ensemble,index)%ID.eq.pion) then
             mom=pertParticles(ensemble,index)%momentum
             if (CMFrame) then
                !Boost momenta from lab frame to cm frame
                beta=(/0.,0., sqrt((mPi+getekin())**2-mPi**2)/&
                     & (getekin()+mPi+mN)/)
                call lorentz(beta,mom)
             end if

             ! Evaluating the angular distribution
             absMom_perpendicular=SQRT(mom(1)**2+mom(2)**2)
             if (abs(mom(3)).lt.0.000001) then
                theta=90.
             else if (mom(3).gt.0.) then
                theta=ATAN(absMom_perpendicular/mom(3))
                theta=theta/pi*180.
             else
                theta=pi-ATAN(absMom_perpendicular/Abs(mom(3)))
                theta=theta/pi*180.
             end if
             if (debug) write(*,*) mom(1:3) , 'theta=',theta
             thetaIndex=Max(Min(180,NINT(theta/deltaTheta)),0)

             sigma_dTheta(pertParticles(ensemble,index)%charge,thetaIndex)=&
                  & sigma_dTheta(pertParticles(ensemble,index)%charge,thetaIndex) &
                  &       +pertParticles(ensemble,index)%perweight/(deltaTheta*pi/180.)
          end if
       end do indexSchleife
    end do ensembleSchleife

    write(*,*) 'sigmaQuasi=', sigma_Quasi

    ! (2) Add up cross sections of all runs
    do i=-1,1
       quasiElastic_xSection(i,:)=quasiElastic_XSection(i,:)+sigma_Quasi(i,:)
       quasiElastic_error(i,:)=quasiElastic_error(i,:)+sigma_Quasi(i,:)**2
    end do

    ! (3) Divide by number of runs at the end
    if (finalFlag) then
       do i=-1,1
          quasiElastic_xSection(i,:)=quasiElastic_xSection(i,:)/float(numberRuns)
       end do
       sigma_dTheta(-1,:) =sigma_dTheta (-1,:)/float(numberRuns)
       sigma_dTheta( 0,:) =sigma_dTheta ( 0,:)/float(numberRuns)
       sigma_dTheta( 1,:) =sigma_dTheta ( 1,:)/float(numberRuns)
       absorption_xSection=absorption_xSection/float(numberRuns)
       sigmaTotalsave=sigmaTotalSave / float(numberruns)

       if (numberruns.gt.1) then
          ! Evaluate standard deviations
          sigma_total_error     =sqrt(abs(sigma_Total_error-sigmaTotalSave**2*float(numberruns))&
               & /float(numberruns-1)/float(numberruns))
          do i=-1,1
             quasielastic_total_error(i,:) =sqrt(abs(quasiElastic_error(i,:)&
                  & -quasiElastic_xSection(i,:)**2*float(numberruns))&
                  &/float(numberruns-1)/float(numberruns))
          end do
          sigma_absorption_error=sqrt(abs(sigma_absorption_error-absorption_xSection**2*float(numberruns))/&
               & float(numberruns-1)/float(numberruns))
       else
          sigma_total_error=9999.
          sigma_absorption_error=9999.
          do i=-1,1
             quasielastic_total_error(i,:) =9999.
          end do
       end if


       if (firstFlag) then
          open(140,File='./pionInduced_xSections.dat')
          open(240,File='./pionInduced_QE_generation.dat')
          write(140,'(A)') '# elab, Sigma piMinus, Sigma piNull,Sigma piPlus, Sigma_QElastic' // &
                        &  ', absorption_xSection, sigma Total, sigma Total(check), absorption Events ,number of runs'// &
                        &  ', error of quasiElastic(-1:1), error of absorption_xSection, error of sigma Total'
          write(240,'(A)') '# Quasi elastic cross section as function of "generation"'
          write(240,'(A)') '# elab, Sigma piMinus(0:10), Sigma piNull(0:10),Sigma piPlus(0:10)'
          write(240,'(A)') '# E.g. "Sigma piMinus(5)" includes only pi^- particles of 5th generation'
          firstFlag=.false.
       else
          open(140,File='./pionInduced_xSections.dat',position='append')
          open(240,File='./pionInduced_QE_generation.dat',position='append')
       end if

       do i=-1,1
          sigma_QE(i)=Sum(quasiElastic_xSection(i,:))
          sigma_QE_error(i)=Sum(quasiElastic_total_error(i,:))
       end do

       sigTot=sum(sigma_QE)+absorption_xSection

       write(140,'(20G12.4)') getEKin(),sigma_QE, sum(sigma_QE), absorption_xSection, &
            & sigTot, sigmaTotalsave,absEvents, numberRuns,sigma_QE_error,&
            & sigma_absorption_error, sigma_total_error
       close(140)

       write(240,'(33F10.4)') getEKin(),quasiElastic_xSection(-1,:),&
            & quasiElastic_xSection(0,:),quasiElastic_xSection(1,:)
       close(240)


       ! Angular distributions
       fileName='pionInduced_dTheta'//RealToChar(getEkin()*1000.) //'MeV.dat'
       !       Print *, filename
       open(170,File=filename)
       write(170,*) '# elab [GeV] =',getEkin()
       write(170,*) '# theta angle[°] , Sigma piMinus, Sigma piNull,Sigma piPlus,Sigma_QElastic'
       do i=lbound(sigma_dtheta,dim=2),ubound(sigma_dtheta,dim=2)
          write(170,'(5F12.3)') float(i)*deltaTheta, sigma_dTheta(-1:1,i), sum(sigma_dTheta(:,i))
       end do
       close(170)
    else
       open(140,file='pionInduced_prelim.dat',position='append')
       write(140,'(A)') '# elab, Sigma piMinus, Sigma piNull,Sigma piPlus, Sigma_QElastic' // &
            &  ', absorption_xSection, sigma Total, number of runs,'

       do i=-1,1
          sigma_QE(i)=Sum(quasiElastic_xSection(i,:))
       end do

       sigTot=sum(sigma_QE)+absorption_xSection

       write(140,'(7F10.4 ,1I12,5F10.4)') getEKin(),sigma_QE/float(numberRuns), &
            & Sum(sigma_QE)/float(numberRuns), absorption_xSection/float(numberRuns)&
            & ,sigmaTotalSave/float(numberruns), &
            & numberRuns
       close(140)


       open(240,File='./pionInduced_QE_generation_prelim.dat',position='append')
       write(240,'(A)') '# Quasi elastic cross section as function of "generation"'
       write(240,'(A)') '# elab, Sigma piMinus(0:10), Sigma piNull(0:10),Sigma piPlus(0:10)'
       write(240,'(A)') '# E.g. "Sigma piMinus(5)" includes only pi^- particles of 5th generation'
       write(240,'(33F10.4)') getEKin(),quasiElastic_xSection(-1,:)/float(numberruns) &
            & ,quasiElastic_xSection(0,:)/float(numberruns),quasiElastic_xSection(1,:)/float(numberruns)
       close(240)


       fileName='pionInduced_dTheta'//RealToChar(getEkin()*1000.) //'MeV.dat'
       open(170,File=filename)
       write(170,*) '# elab [GeV] =',getEkin()
       write(170,*) '# theta angle[°] , Sigma piMinus, Sigma piNull,Sigma piPlus,Sigma_QElastic'
       do i=lbound(sigma_dtheta,dim=2),ubound(sigma_dtheta,dim=2)
          write(170,'(5F12.3)') float(i)*deltaTheta, sigma_dTheta(-1:1,i)/float(numberRuns), &
               & sum(sigma_dTheta(:,i))/float(numberRuns)
       end do
       close(170)
    end if



    ! dSigma/dOmega/dE_pion *********************************************************
    if (dsigma_dOmegadE_switch) then

       call dSigmadOmegadE(pertParticles,finalFlag)
       if (.not.localEnsemble) then
          ! dSigma/dOmega/dE_pion out of NN collisions*********************************************************
          call dSigmadOmegadE_NN(pertParticles,finalFlag)

          ! dSigma/dOmega/dE_pion out of Delta decays*********************************************************
          call dSigmadOmegadE_Resonance(pertParticles,finalFlag)

          ! dSigma/dOmega/dE_pion out of N pion direct collisions*********************************************************
          call dSigmadOmegadE_pionNuk(pertParticles,finalFlag)
       end if
    end if

    ! dsigma(pi pi)/dm and sigma(pi pi)
    if (twoPi_switch) call twoPi_output(pertParticles,finalFlag,getEkin())
    !*****************************Write data to file***************************************

    startFlag=.false.
    if (finalFlag) then
       startFlag=.true.
    end if
    return
  end subroutine pionXsectionAnalysis


  !****************************************************************************
  !****s* pionXsection/dSigmadOmegadE
  ! NAME
  ! subroutine dSigmadOmegadE(pertParticles,flag)
  ! INPUTS
  ! type(particle), dimension(:,:)  :: pertParticles  ! The particle vector
  ! logical :: flag ! Whether it's the final run
  ! OUTPUT
  ! The files:
  ! * pionInduced_dOmegadE_piMinus*MeV.dat
  ! * pionInduced_dOmegadE_piNull*MeV.dat
  ! * pionInduced_dOmegadE_piPlus*MeV.dat
  !
  ! where * is the kinetic energy of the incoming pions.
  ! PURPOSE
  ! This subroutine produces output for pion-nucleus scattering, including
  ! impact parameter integration.
  ! Gives dSigma/dOmega/dE(pion) for Pions in file.
  !****************************************************************************
  subroutine dSigmadOmegadE(pertParticles,finalFlag)

    use idTable
    use output
    use constants, only: pi, mN, mPi
    use particleDefinition
    use initPion, only: getEkin
    use lorentzTrafo, only: lorentz

    type(particle), intent(in),dimension(:,:)  :: pertParticles
    logical, intent(in) :: finalFlag

    integer, parameter :: ekin_dim=50
    integer, parameter :: theta_dim=18

    !real, save :: delta_ekin=0.240/float(ekin_dim)
    !real, save :: delta_theta=180./float(theta_dim)

    real, dimension(0:ekin_Dim,0:Theta_dim,-1:1),save  :: dsigma,dsigma_total, dsigma_squared,error

    real, dimension(0:Theta_dim,-1:1),save  :: squared_dTheta, new_dTheta,old_dTheta


    real :: theta, absMom_perpendicular,theta_degree,ekin

    integer :: ekinIndex, thetaIndex
    !logical :: first

    integer :: i,j,k
    integer,save :: numberruns=0
    logical, save :: startFlag=.true.

    integer :: ensemble, index,file
    real, dimension(1:3) :: beta
    real, dimension(0:3) :: mom
    real :: error_dTheta_old,error_dTheta_new

    real, save :: delta_ekin !=0.240/float(ekin_dim)
    real, save :: delta_theta!=180./float(theta_dim)
    delta_ekin=0.240/float(ekin_dim)
    delta_theta=180./float(theta_dim)



    ! Index 1 : Ekin of final pion, Index 2 : theta of final pion, Index 3: charge of final pion

    if (startFlag) then
       do i=lbound(dsigma,dim=1),ubound(dsigma,dim=1)
          do j=lbound(dsigma,dim=2),ubound(dsigma,dim=2)
             do k=lbound(dsigma,dim=3),ubound(dsigma,dim=3)
                dsigma_total(i,j,k)=0.
                dsigma_squared(i,j,k)=0.
             end do
          end do
       end do
       do i=lbound(squared_dTheta,dim=1),ubound(squared_dTheta,dim=1)
          do j=lbound(squared_dTheta,dim=2),ubound(squared_dTheta,dim=2)
             squared_dTheta(i,j)=0.
             old_dTheta(i,j)=0.
             new_dTheta(i,j)=0.
          end do
       end do
       if (debug) then
          open(22,File='pionAnaTest.txt')
          write(22,*) 'dTheta=', delta_theta
          write(22,*) 'dEkin=', delta_ekin
          close(22)
       end if
       numberruns=0
       !first=.true.
       startFlag=.false.
    end if

    numberruns=numberruns+1

    do i=lbound(dsigma,dim=1),ubound(dsigma,dim=1)
       do j=lbound(dsigma,dim=2),ubound(dsigma,dim=2)
          do k=lbound(dsigma,dim=3),ubound(dsigma,dim=3)
             dsigma(i,j,k)=0.
          end do
       end do
    end do

    ensembleSchleife : do ensemble=1,size(pertParticles,dim=1)
       indexSchleife : do index=1,size(pertParticles,dim=2)

          if (pertParticles(ensemble,index)%ID <= 0) cycle indexSchleife

          !          if((pertParticles(ensemble,index)%ID.eq.pion).and.(pertParticles(ensemble,index)%lastCollisionTime.ge.0.0001)) then
          if (pertParticles(ensemble,index)%ID.eq.pion) then
             ! Evaluating the angular distribution
             mom=pertParticles(ensemble,index)%momentum(0:3)
             if (CMFrame) then
                !Boost momenta from lab frame to cm frame
                beta=(/0.,0., sqrt((mPi+getekin())**2-mPi**2)/&
                     & (getekin()+mPi+mN)/)
                call lorentz(beta,mom)
             end if

             absMom_perpendicular=SQRT(mom(1)**2+mom(2)**2)
             if (abs(mom(3)).lt.0.000001) then
                theta=pi/2.
                theta_degree=90.
             else if (mom(3).gt.0.) then
                theta=ATAN(absMom_perpendicular/mom(3))
                theta_degree=theta/pi*180.
             else
                theta=pi-ATAN(absMom_perpendicular/Abs(mom(3)))
                theta_degree=theta/pi*180.
             end if
             if (debug) write(*,*) pertParticles(ensemble,index)%momentum(1:3) , 'theta=',theta


             if (cmFrame) then
                ekin=sqrt(pertParticles(ensemble,index)%mass**2+mom(1)**2+mom(2)**2+mom(3)**2)&
                     & -pertParticles(ensemble,index)%mass
             else
                ekin=kineticEnergy(pertParticles(ensemble,index))
             end if

             thetaIndex=Max(Min(180,NINT(theta_degree/delta_Theta)),0)
             ekinIndex=Max(Min(ekin_dim,NINT(ekin/delta_ekin)),0)

             if (abs(sin(theta)).gt.0) then
                dsigma(ekinIndex,thetaIndex,pertParticles(ensemble,index)%charge)=&
                     & dsigma(ekinIndex,thetaIndex,pertParticles(ensemble,index)%charge) &
                     &       + pertParticles(ensemble,index)%perweight/ &
                     & (delta_Theta*pi/180.*sin(theta)*2.*pi)/delta_ekin
             end if
          end if
       end do indexSchleife
    end do ensembleSchleife


    ! FOR ERROR ANALYSIS, used for standard deviations among single runs
    do i=lbound(dsigma,dim=1),ubound(dsigma,dim=1)
       do j=lbound(dsigma,dim=2),ubound(dsigma,dim=2)
          do k=lbound(dsigma,dim=3),ubound(dsigma,dim=3)
             dsigma_squared(i,j,k)=dsigma_squared(i,j,k)+dsigma(i,j,k)**2
          end do
       end do
    end do


    ! (2) Add up cross sections of all runs
    do i=lbound(dsigma,dim=1),ubound(dsigma,dim=1)
       do j=lbound(dsigma,dim=2),ubound(dsigma,dim=2)
          dsigma_total(i,j,:)=dsigma_total(i,j,:)+dsigma(i,j,:)
       end do
    end do


    ! (3) Write output

    open(144,File='./pionInduced_dOmegadE_piMinus'//RealToChar(getEkin()*1000.) //'MeV.dat')
    open(145,File='./pionInduced_dOmegadE_piNull'//RealToChar(getEkin()*1000.) //'MeV.dat')
    open(146,File='./pionInduced_dOmegadE_piPlus'//RealToChar(getEkin()*1000.) //'MeV.dat')
    write(144,'(A,f7.4,A,I4)') '# Final Pi^- distribution. Initial energy:', getEKin(), &
         &  '. Number of runs:', numberruns
    write(145,'(A,f7.4,A,I4)') '# Final Pi^0 distribution. Initial energy:', getEKin(),  &
         &'. Number of runs:', numberruns
    write(146,'(A,f7.4,A,I4)') '# Final Pi^+ distribution. Initial energy:', getEKin(),  &
         &'. Number of runs:', numberruns
    write(144,'(2A)') '# ePion [GeV], theta [degree], dsigma/dOmega/dE [mB/GeV/sr],' &
         & ,' error of dsigma/dOmega/dE[mB/GeV/sr]'
    write(145,'(2A)') '# ePion [GeV], theta [degree], dsigma/dOmega/dE [mB/GeV/sr],' &
         & ,' error of dsigma/dOmega/dE[mB/GeV/sr]'
    write(146,'(2A)') '# ePion [GeV], theta [degree], dsigma/dOmega/dE [mB/GeV/sr],' &
         & ,' error of dsigma/dOmega/dE[mB/GeV/sr]'


     do i=lbound(dsigma_total,dim=1),ubound(dsigma_total,dim=1)
       do j=lbound(dsigma_total,dim=2),ubound(dsigma_total,dim=2)
          do k=lbound(dsigma_total,dim=3),ubound(dsigma_total,dim=3)
             if (numberruns.gt.1) then
                ! Evaluate standard deviations of mean value
                error(i,j,k) =sqrt(abs(dsigma_squared(i,j,k)-dsigma_total(i,j,k)**2/float(numberruns)) &
                     & /float(numberruns)/float(numberruns-1))
             else
                error(i,j,k)  =999.
             end if

             if (debug) then
                if ((j.eq.lbound(dsigma_total,dim=2)).or.(j.eq.ubound(dsigma_total,dim=2))) then
                   write(*,*) float(i)*delta_ekin, float(j)*delta_theta, dsigma_total(i,j,k)*2.,  error*2
                else
                   write(*,*) float(i)*delta_ekin, float(j)*delta_theta, dsigma_total(i,j,k),  error
                end if
             end if

             file=145+k
             if ((j.eq.lbound(dsigma_total,dim=2)).or.(j.eq.ubound(dsigma_total,dim=2))) then
                write(file,'(4F18.7)') float(i)*delta_ekin, float(j)*delta_theta, &
                     & dsigma_total(i,j,k)/float(numberRuns)*2., error(i,j,k)*2.
             else
                write(file,'(4F18.7)') float(i)*delta_ekin, float(j)*delta_theta, &
                     & dsigma_total(i,j,k)/float(numberRuns), error(i,j,k)
             end if


             write(file,*)
          end do
       end do
    end do
    close(144)
    close(145)
    close(146)

    ! Integrated over DE

    open(244,File='./pionInduced_dOmega_piMinus'//RealToChar(getEkin()*1000.) //'MeV.dat')
    open(245,File='./pionInduced_dOmega_piNull'//RealToChar(getEkin()*1000.) //'MeV.dat')
    open(246,File='./pionInduced_dOmega_piPlus'//RealToChar(getEkin()*1000.) //'MeV.dat')
    write(244,'(A,f7.4,A,I4)') '# Final Pi^- distribution. Initial energy:', &
& getEKin(), '. Number of runs:', numberruns
    write(245,'(A,f7.4,A,I4)') '# Final Pi^0 distribution. Initial energy:', &
&  getEKin(), '. Number of runs:', numberruns
    write(246,'(A,f7.4,A,I4)') '# Final Pi^+ distribution. Initial energy:',  &
& getEKin(), '. Number of runs:', numberruns
    write(244,'(A)') '# theta [degree], dsigma/dOmega [mB/sr], error of dsigma/dOmega[mB/sr]'
    write(245,'(A)') '# theta [degree], dsigma/dOmega [mB/sr], error of dsigma/dOmega[mB/sr]'
    write(246,'(A)') '# theta [degree], dsigma/dOmega [mB/sr], error of dsigma/dOmega[mB/sr]'



    ! new_dTheta    =Present  sigma_dTheta of all combined runs
    !               = 1/numberRuns Sum_{i=1,numberRuns} ( sigma of single run)
    ! old_dTheta    =Previous sigma_dTheta of all combined runs
    !               = 1/(numberRuns-1) Sum_{i=1,numberRuns-1} ( sigma of single run)
    ! squared_dTheta=Sum over all runs of  (sigma_dTheta of a single run)
    !               =Sum_Runs [ ( new_dTheta*(numberRuns) - old_dTheta*(numberRuns-1) )**2 ]

    do j=lbound(dsigma_total,dim=2),ubound(dsigma_total,dim=2)
       do k=lbound(dsigma_total,dim=3),ubound(dsigma_total,dim=3)
          if (debug) then
             write(*,*) float(j)*delta_theta, sum(dsigma_total(:,j,k)),  error
          end if
          file=245+k

          if ((j.eq.lbound(dsigma_total,dim=2)).or.(j.eq.ubound(dsigma_total,dim=2))) then
             new_dTheta(j,k)=sum(dsigma_total(:,j,k))/float(numberRuns)*delta_ekin*2.
             if (numberruns.gt.1) then
                squared_dTheta(j,k)=squared_dTheta(j,k)+(new_dTheta(j,k)*numberRuns-old_dTheta(j,k)*(numberRuns-1))**2
             else
                squared_dTheta(j,k)=new_dTheta(j,k)**2
             end if
             error_dTheta_old=sum(error(:,j,k))*delta_ekin*2.
          else
             new_dTheta(j,k)=sum(dsigma_total(:,j,k))/float(numberRuns)*delta_ekin
             if (numberruns.gt.1) then
                squared_dTheta(j,k)=squared_dTheta(j,k)+(new_dTheta(j,k)*numberRuns-old_dTheta(j,k)*(numberRuns-1))**2
             else
                squared_dTheta(j,k)=new_dTheta(j,k)**2
             end if
             error_dTheta_old=sum(error(:,j,k))*delta_ekin
          end if
          if (numberruns.gt.1) then
             error_dTheta_new=sqrt(abs(squared_dTheta(j,k)-new_dTheta(j,k)**2*numberRuns) &
                  & /float(numberruns)/float(numberruns-1))
          else
             error_dTheta_new=999.
          end if
          write(file,'(4F18.7)') float(j)*delta_theta, new_dTheta(j,k),error_dTheta_old,error_dTheta_new
       end do
    end do

    do j=lbound(old_dTheta,dim=1),ubound(old_dTheta,dim=1)
       do k=lbound(old_dTheta,dim=2),ubound(old_dTheta,dim=2)
          old_dTheta(j,k)=new_dTheta(j,k)
       end do
    end do

    close(244)
    close(245)
    close(246)


    if (finalflag) then
       startFlag=.true.
    end if


  end subroutine dSigmadOmegadE




  !****************************************************************************
  !****s* pionXsection/dSigmadOmegadENN
  ! NAME
  ! subroutine dSigmadOmegadE_NN(pertParticles,finalFlag)
  ! INPUTS
  ! type(particle), intent(in),dimension(:,:)  :: pertParticles  ! The particle vector
  ! logical,intent(in) :: flag ! Whether it's the final run
  ! OUTPUT
  ! The files:
  ! * pionInduced_dOmegadE_piMinus_NN*MeV.dat
  ! * pionInduced_dOmegadE_piNull_NN*MeV.dat
  ! * pionInduced_dOmegadE_piPlus_NN*MeV.dat
  !
  ! where * is the kinetic energy of the incoming pions.
  ! PURPOSE
  ! This subroutine produces output for pion-nucleus scattering,
  ! including impactp arameter integration
  !
  ! Gives dSigma/dOmega/dE(pion) for Pions in file.
  ! Only those pions are considered which are produced in NN collisions.
  ! NOTES
  ! Makes only sense if useStatistics=.true. in collisionTerm
  !****************************************************************************
  subroutine dSigmadOmegadE_NN(pertParticles,finalFlag)

    use idTable
    use output
    use constants, only: pi, mN, mPi
    use particleDefinition
    use initPion, only: getEkin
    use statistics, only: getInfo
    use collisionNumbering, only:pert_numbering
    use lorentzTrafo, only: lorentz

    type(particle), intent(in),dimension(:,:)  :: pertParticles
    logical, intent(in) :: finalFlag

    integer, parameter :: ekin_dim=50
    integer, parameter :: theta_dim=18

    !real, save :: delta_ekin=0.240/float(ekin_dim)
    !real, save :: delta_theta=180./float(theta_dim)

    real, dimension(0:ekin_Dim,0:Theta_dim,-1:1),save  :: dsigma,dsigma_total, dsigma_squared

    real :: theta, absMom_perpendicular,error,theta_degree,ekin

    integer :: ekinIndex, thetaIndex
    !logical :: first

    integer :: i,j,k
    integer,save :: numberruns=0
    logical, save :: startFlag=.true.

    integer :: ensemble, index,file

    logical :: getInfoFlag
    real, dimension(1:3) :: prodPlace
    integer, dimension(1:3) :: prodParticles

    real , dimension(1:3) :: beta
    real , dimension(0:3) :: mom


    real, save :: delta_ekin!=0.240/float(ekin_dim)
    real, save :: delta_theta!=180./float(theta_dim)
    delta_ekin=0.240/float(ekin_dim)
    delta_theta=180./float(theta_dim)


    ! Index 1 : Ekin of final pion, Index 2 : theta of final pion, Index 3: charge of final pion

    if (startFlag) then
       do i=lbound(dsigma,dim=1),ubound(dsigma,dim=1)
          do j=lbound(dsigma,dim=2),ubound(dsigma,dim=2)
             do k=lbound(dsigma,dim=3),ubound(dsigma,dim=3)
                dsigma_total(i,j,k)=0.
                dsigma_squared(i,j,k)=0.
             end do
          end do
       end do
       if (debug) then
          open(22,File='pionAnaTest.txt')
          write(22,*) 'dTheta=', delta_theta
          write(22,*) 'dEkin=', delta_ekin
          close(22)
       end if
       numberruns=0
       !first=.true.
       startFlag=.false.
    end if

    numberruns=numberruns+1

    do i=lbound(dsigma,dim=1),ubound(dsigma,dim=1)
       do j=lbound(dsigma,dim=2),ubound(dsigma,dim=2)
          do k=lbound(dsigma,dim=3),ubound(dsigma,dim=3)
             dsigma(i,j,k)=0.
          end do
       end do
    end do

    ensembleSchleife : do ensemble=1,size(pertParticles,dim=1)
       indexSchleife : do index=1,size(pertParticles,dim=2)

          if (pertParticles(ensemble,index)%ID <= 0) cycle indexSchleife

          if (pertParticles(ensemble,index)%event(1).ne.pert_numbering()) then
             !          if(pertParticles(ensemble,index)%lastCollisionTime.ge.0.0001) then
             if (pertParticles(ensemble,index)%ID.eq.pion) then
                call getInfo(pertParticles(ensemble,index)%number,ensemble,prodPlace,prodParticles,getInfoFlag)
                if (.not.getInfoFlag) then
                   !write(*,*) 'Problem in getting the information with getinfo'
                   cycle indexSchleife
                end if
                if (.not.((prodParticles(1).eq.nucleon).and.(prodParticles(2).eq.nucleon))) then
                   ! Pion not produced in NN collision
                   cycle indexSchleife
                end if

                ! Evaluating the angular distribution


                mom=pertParticles(ensemble,index)%momentum(0:3)
                if (CMFrame) then
                   !Boost momenta from lab frame to cm frame
                   beta=(/0.,0., sqrt((mPi+getekin())**2-mPi**2)/&
                        & (getekin()+mPi+mN)/)
                   call lorentz(beta,mom)
                end if

                absMom_perpendicular=SQRT(mom(1)**2+mom(2)**2)
                if (abs(mom(3)).lt.0.000001) then
                   theta=pi/2.
                   theta_degree=90.
                else if (mom(3).gt.0.) then
                   theta=ATAN(absMom_perpendicular/mom(3))
                   theta_degree=theta/pi*180.
                else
                   theta=pi-ATAN(absMom_perpendicular/Abs(mom(3)))
                   theta_degree=theta/pi*180.
                end if
                if (debug) write(*,*) pertParticles(ensemble,index)%momentum(1:3) , 'theta=',theta


                if (cmFrame) then
                   ekin=sqrt(pertParticles(ensemble,index)%mass**2+mom(1)**2+mom(2)**2+mom(3)**2)&
                        & -pertParticles(ensemble,index)%mass
                else
                   ekin=kineticEnergy(pertParticles(ensemble,index))
                end if

                thetaIndex=Max(Min(180,NINT(theta_degree/delta_Theta)),0)
                ekinIndex=Max(Min(ekin_dim,NINT(ekin/delta_ekin)),0)

                if (abs(sin(theta)).gt.0) then
                   dsigma(ekinIndex,thetaIndex,pertParticles(ensemble,index)%charge)=&
                        & dsigma(ekinIndex,thetaIndex,pertParticles(ensemble,index)%charge) &
                        &       + pertParticles(ensemble,index)%perweight/(delta_Theta*pi/180.*sin(theta)*2.*pi)/delta_ekin
                end if
             end if
          end if
       end do indexSchleife
    end do ensembleSchleife


    ! FOR ERROR ANALYSIS, used for standard deviations among single runs
    do i=lbound(dsigma,dim=1),ubound(dsigma,dim=1)
       do j=lbound(dsigma,dim=2),ubound(dsigma,dim=2)
          do k=lbound(dsigma,dim=3),ubound(dsigma,dim=3)
             dsigma_squared(i,j,k)=dsigma_squared(i,j,k)+dsigma(i,j,k)**2
          end do
       end do
    end do


    ! (2) Add up cross sections of all runs
    do i=lbound(dsigma,dim=1),ubound(dsigma,dim=1)
       do j=lbound(dsigma,dim=2),ubound(dsigma,dim=2)
          dsigma_total(i,j,:)=dsigma_total(i,j,:)+dsigma(i,j,:)
       end do
    end do


    ! (3) Write output

    open(144,File='./pionInduced_dOmegadE_piMinus_NN'//RealToChar(getEkin()*1000.) //'MeV.dat')
    open(145,File='./pionInduced_dOmegadE_piNull_NN'//RealToChar(getEkin()*1000.) //'MeV.dat')
    open(146,File='./pionInduced_dOmegadE_piPlus_NN'//RealToChar(getEkin()*1000.) //'MeV.dat')
    write(144,'(A,f7.4,A,I4)') '# Final Pi^- distribution. Initial energy:', getEKin(), '. Number of runs:', numberruns
    write(145,'(A,f7.4,A,I4)') '# Final Pi^0 distribution. Initial energy:', getEKin(), '. Number of runs:', numberruns
    write(146,'(A,f7.4,A,I4)') '# Final Pi^+ distribution. Initial energy:', getEKin(), '. Number of runs:', numberruns
    write(144,'(A)') '# ePion [GeV], theta [degree], dsigma/dOmega/dE [mB/GeV/sr], error of dsigma/dOmega/dE[mB/GeV/sr]'
    write(145,'(A)') '# ePion [GeV], theta [degree], dsigma/dOmega/dE [mB/GeV/sr], error of dsigma/dOmega/dE[mB/GeV/sr]'
    write(146,'(A)') '# ePion [GeV], theta [degree], dsigma/dOmega/dE [mB/GeV/sr], error of dsigma/dOmega/dE[mB/GeV/sr]'

    do i=lbound(dsigma_total,dim=1),ubound(dsigma_total,dim=1)
       do j=lbound(dsigma_total,dim=2),ubound(dsigma_total,dim=2)
          do k=lbound(dsigma_total,dim=3),ubound(dsigma_total,dim=3)
             if (numberruns.gt.1) then
                ! Evaluate standard deviations of mean value
                error =sqrt(abs(dsigma_squared(i,j,k)-dsigma_total(i,j,k)**2/float(numberruns))&
                     & /float(numberruns)/float(numberruns-1))
             else
                error  =999.
             end if

             if (debug) then
                write(*,*) float(i)*delta_ekin, float(j)*delta_theta, dsigma_total(i,j,k),  error
             end if

             file=145+k
             if ((j.eq.lbound(dsigma_total,dim=2)).or.(j.eq.ubound(dsigma_total,dim=2))) then
                write(file,'(4F18.7)') float(i)*delta_ekin, float(j)*delta_theta, &
                     & dsigma_total(i,j,k)/float(numberRuns)*2., error*2.
             else
                write(file,'(4F18.7)') float(i)*delta_ekin, float(j)*delta_theta, &
                     & dsigma_total(i,j,k)/float(numberRuns), error
             end if

             write(file,*)
          end do
       end do
    end do
    close(144)
    close(145)
    close(146)

    ! Integrated over DE

    open(244,File='./pionInduced_dOmega_piMinus_NN'//RealToChar(getEkin()*1000.) //'MeV.dat')
    open(245,File='./pionInduced_dOmega_piNull_NN'//RealToChar(getEkin()*1000.) //'MeV.dat')
    open(246,File='./pionInduced_dOmega_piPlus_NN'//RealToChar(getEkin()*1000.) //'MeV.dat')
    write(244,'(A,f7.4,A,I4)') '# Final Pi^- distribution. Initial energy:', getEKin(), '. Number of runs:', numberruns
    write(245,'(A,f7.4,A,I4)') '# Final Pi^0 distribution. Initial energy:', getEKin(), '. Number of runs:', numberruns
    write(246,'(A,f7.4,A,I4)') '# Final Pi^+ distribution. Initial energy:', getEKin(), '. Number of runs:', numberruns
    write(244,'(A)') '# theta [degree], dsigma/dOmega [mB/sr], error of dsigma/dOmega[mB/sr]'
    write(245,'(A)') '# theta [degree], dsigma/dOmega [mB/sr], error of dsigma/dOmega[mB/sr]'
    write(246,'(A)') '# theta [degree], dsigma/dOmega [mB/sr], error of dsigma/dOmega[mB/sr]'

    do j=lbound(dsigma_total,dim=2),ubound(dsigma_total,dim=2)
       do k=lbound(dsigma_total,dim=3),ubound(dsigma_total,dim=3)
          if (debug) then
             write(*,*) float(j)*delta_theta, sum(dsigma_total(:,j,k)),  error
          end if
          file=245+k
          if ((j.eq.lbound(dsigma_total,dim=2)).or.(j.eq.ubound(dsigma_total,dim=2))) then
             write(file,'(4F18.7)') float(j)*delta_theta, sum(dsigma_total(:,j,k))/float(numberRuns)*delta_ekin*2.
          else
             write(file,'(4F18.7)') float(j)*delta_theta, sum(dsigma_total(:,j,k))/float(numberRuns)*delta_ekin
          end if

       end do
    end do

    close(244)
    close(245)
    close(246)

    if (finalflag) then
       startFlag=.true.
    end if


  end subroutine dSigmadOmegadE_NN


  !****************************************************************************
  !****s* pionXsection/dSigmadOmegadE_Resonance
  ! NAME
  ! subroutine dSigmadOmegadE_Resonance(pertParticles,finalFlag)
  ! INPUTS
  ! type(particle), intent(in),dimension(:,:)  :: pertParticles  ! The particle vector
  ! logical,intent(in) :: flag ! Whether it's the final run
  ! OUTPUT
  ! The files:
  ! * pionInduced_dOmegadE_piMinus_R*MeV.dat
  ! * pionInduced_dOmegadE_piNull_R*MeV.dat
  ! * pionInduced_dOmegadE_piPlus_R*MeV.dat
  !
  ! where * is the kinetic energy of the incoming pions.
  ! PURPOSE
  ! This subroutine produces output for pion-nucleus scattering,
  ! including impact parameter integration
  !
  ! Gives dSigma/dOmega/dE(pion) for Pions in file.
  ! Only those pions are considered which are produced in Resonance decays.
  ! NOTES
  ! Makes only sense if useStatistics=.true. in collisionTerm
  !****************************************************************************
  subroutine dSigmadOmegadE_Resonance(pertParticles,finalFlag)

    use idTable
    use output
    use constants, only: pi, mN, mPi
    use particleDefinition
    use initPion, only: getEkin
    use statistics, only: getInfo
    use collisionNumbering, only:pert_numbering
    use lorentzTrafo, only: lorentz

    type(particle), intent(in),dimension(:,:)  :: pertParticles
    logical, intent(in) :: finalFlag

    integer, parameter :: ekin_dim=50
    integer, parameter :: theta_dim=18

    !real, save :: delta_ekin=0.240/float(ekin_dim)
    !real, save :: delta_theta=180./float(theta_dim)

    real, dimension(0:ekin_Dim,0:Theta_dim,-1:1,Delta:F37_1950),save  :: dsigma,dsigma_total, dsigma_squared
    real, dimension(Delta:F37_1950),save  :: error,summe

    real :: theta, absMom_perpendicular,theta_degree,ekin

    integer :: ekinIndex, thetaIndex
    !logical :: first

    integer :: i,j,k
    integer,save :: numberruns=0
    logical, save :: startFlag=.true.

    integer :: ensemble, index,file

    logical :: getInfoFlag
    real, dimension(1:3) :: prodPlace
    integer, dimension(1:3) :: prodParticles
    real , dimension(1:3) :: beta
    real , dimension(0:3) :: mom

    real, save :: delta_ekin !=0.240/float(ekin_dim)
    real, save :: delta_theta !=180./float(theta_dim)
    delta_ekin=0.240/float(ekin_dim)
    delta_theta=180./float(theta_dim)


    ! Index 1 : Ekin of final pion, Index 2 : theta of final pion, Index 3: charge of final pion

    if (startFlag) then
       do i=lbound(dsigma,dim=1),ubound(dsigma,dim=1)
          do j=lbound(dsigma,dim=2),ubound(dsigma,dim=2)
             do k=lbound(dsigma,dim=3),ubound(dsigma,dim=3)
                dsigma_total(i,j,k,:)=0.
                dsigma_squared(i,j,k,:)=0.
             end do
          end do
       end do
       if (debug) then
          open(22,File='pionAnaTest.txt')
          write(22,*) 'dTheta=', delta_theta
          write(22,*) 'dEkin=', delta_ekin
          close(22)
       end if
       numberruns=0
       !first=.true.
       startFlag=.false.
    end if

    numberruns=numberruns+1

    do i=lbound(dsigma,dim=1),ubound(dsigma,dim=1)
       do j=lbound(dsigma,dim=2),ubound(dsigma,dim=2)
          do k=lbound(dsigma,dim=3),ubound(dsigma,dim=3)
             dsigma(i,j,k,:)=0.
          end do
       end do
    end do

    ensembleSchleife : do ensemble=1,size(pertParticles,dim=1)
       indexSchleife : do index=1,size(pertParticles,dim=2)

          if (pertParticles(ensemble,index)%ID <= 0) cycle indexSchleife

          if (pertParticles(ensemble,index)%event(1).ne.pert_numbering()) then
             !          if(pertParticles(ensemble,index)%lastCollisionTime.ge.0.0001) then
             if (pertParticles(ensemble,index)%ID.eq.pion) then
                call getInfo(pertParticles(ensemble,index)%number,ensemble,prodPlace,prodParticles,getInfoFlag)
                if (.not.getInfoFlag) then
                   !write(*,*) 'Problem in getting the information with getinfo'
                   cycle indexSchleife
                end if
                if (  (.not.((prodParticles(1).ge.delta).and.(prodParticles(1).le.F37_1950).and.(prodParticles(2).eq.0))).and. &
                     (.not.((prodParticles(2).ge.delta).and.(prodParticles(2).le.F37_1950).and.(prodParticles(1).eq.0)))     ) then
                   ! Pion not produced in resonance decay
                   cycle indexSchleife
                end if

                ! Evaluating the angular distribution
                mom=pertParticles(ensemble,index)%momentum(0:3)
                if (CMFrame) then
                   !Boost momenta from lab frame to cm frame
                   beta=(/0.,0., sqrt((mPi+getekin())**2-mPi**2)/&
                        & (getekin()+mPi+mN)/)
                   call lorentz(beta,mom)
                end if

                absMom_perpendicular=SQRT(mom(1)**2+mom(2)**2)
                if (abs(mom(3)).lt.0.000001) then
                   theta=pi/2.
                   theta_degree=90.
                else if (mom(3).gt.0.) then
                   theta=ATAN(absMom_perpendicular/mom(3))
                   theta_degree=theta/pi*180.
                else
                   theta=pi-ATAN(absMom_perpendicular/Abs(mom(3)))
                   theta_degree=theta/pi*180.
                end if
                if (debug) write(*,*) pertParticles(ensemble,index)%momentum(1:3) , 'theta=',theta


                if (cmFrame) then
                   ekin=sqrt(pertParticles(ensemble,index)%mass**2+mom(1)**2+mom(2)**2+mom(3)**2)-pertParticles(ensemble,index)%mass
                else
                   ekin=kineticEnergy(pertParticles(ensemble,index))
                end if

                thetaIndex=Max(Min(180,NINT(theta_degree/delta_Theta)),0)
                ekinIndex=Max(Min(ekin_dim,NINT(ekin/delta_ekin)),0)



                if (abs(sin(theta)).gt.0) then
                   if (prodParticles(2).eq.0) then
                      dsigma(ekinIndex,thetaIndex,pertParticles(ensemble,index)%charge,prodParticles(1))=&
                           & dsigma(ekinIndex,thetaIndex,pertParticles(ensemble,index)%charge,prodParticles(1)) &
                           &       + pertParticles(ensemble,index)%perweight/(delta_Theta*pi/180.*sin(theta)*2.*pi)/delta_ekin
                   else if (prodParticles(1).eq.0) then
                      dsigma(ekinIndex,thetaIndex,pertParticles(ensemble,index)%charge,prodParticles(2))=&
                           & dsigma(ekinIndex,thetaIndex,pertParticles(ensemble,index)%charge,prodParticles(2)) &
                           &       + pertParticles(ensemble,index)%perweight/(delta_Theta*pi/180.*sin(theta)*2.*pi)/delta_ekin
                   end if
                end if
             end if
          end if
       end do indexSchleife
    end do ensembleSchleife


    ! FOR ERROR ANALYSIS, used for standard deviations among single runs
    do i=lbound(dsigma,dim=1),ubound(dsigma,dim=1)
       do j=lbound(dsigma,dim=2),ubound(dsigma,dim=2)
          do k=lbound(dsigma,dim=3),ubound(dsigma,dim=3)
             dsigma_squared(i,j,k,:)=dsigma_squared(i,j,k,:)+dsigma(i,j,k,:)**2
          end do
       end do
    end do


    ! (2) Add up cross sections of all runs
    do i=lbound(dsigma,dim=1),ubound(dsigma,dim=1)
       do j=lbound(dsigma,dim=2),ubound(dsigma,dim=2)
          do k=lbound(dsigma,dim=3),ubound(dsigma,dim=3)
             dsigma_total(i,j,k,:)=dsigma_total(i,j,k,:)+dsigma(i,j,k,:)
          end do
       end do
    end do


    ! (3) Write output

    open(144,File='./pionInduced_dOmegadE_piMinus_R'//RealToChar(getEkin()*1000.) //'MeV.dat')
    open(145,File='./pionInduced_dOmegadE_piNull_R'//RealToChar(getEkin()*1000.) //'MeV.dat')
    open(146,File='./pionInduced_dOmegadE_piPlus_R'//RealToChar(getEkin()*1000.) //'MeV.dat')
    write(144,'(A,f7.4,A,I4)') '# Final Pi^- distribution. Initial energy:', getEKin(), '. Number of runs:', numberruns
    write(145,'(A,f7.4,A,I4)') '# Final Pi^0 distribution. Initial energy:', getEKin(), '. Number of runs:', numberruns
    write(146,'(A,f7.4,A,I4)') '# Final Pi^+ distribution. Initial energy:', getEKin(), '. Number of runs:', numberruns
    write(144,'(A)') '# ePion [GeV], theta [degree], dsigma/dOmega/dE [mB/GeV/sr], error of dsigma/dOmega/dE[mB/GeV/sr]'
    write(145,'(A)') '# ePion [GeV], theta [degree], dsigma/dOmega/dE [mB/GeV/sr], error of dsigma/dOmega/dE[mB/GeV/sr]'
    write(146,'(A)') '# ePion [GeV], theta [degree], dsigma/dOmega/dE [mB/GeV/sr], error of dsigma/dOmega/dE[mB/GeV/sr]'

    do i=lbound(dsigma_total,dim=1),ubound(dsigma_total,dim=1)
       do j=lbound(dsigma_total,dim=2),ubound(dsigma_total,dim=2)
          do k=lbound(dsigma_total,dim=3),ubound(dsigma_total,dim=3)
             if (numberruns.gt.1) then
                ! Evaluate standard deviations of mean value
                error =sqrt(abs(dsigma_squared(i,j,k,:)-dsigma_total(i,j,k,:)**2/float(numberruns))/float(numberruns)&
                     & /float(numberruns-1))
             else
                error  =999.
             end if

             if (debug) then
                write(*,*) float(i)*delta_ekin, float(j)*delta_theta, dsigma_total(i,j,k,:),  error
             end if

             file=145+k
             if ((j.eq.lbound(dsigma_total,dim=2)).or.(j.eq.ubound(dsigma_total,dim=2))) then
                write(file,'(4F62.7)') float(i)*delta_ekin, float(j)*delta_theta, &
                     & dsigma_total(i,j,k,:)/float(numberRuns)*2., error(:)*2.
             else
                write(file,'(4F62.7)') float(i)*delta_ekin, float(j)*delta_theta, &
                     & dsigma_total(i,j,k,:)/float(numberRuns), error(:)
             end if
             write(file,*)
          end do
       end do
    end do
    close(144)
    close(145)
    close(146)

    ! Integrated over DE

    open(244,File='./pionInduced_dOmega_piMinus_Delta'//RealToChar(getEkin()*1000.) //'MeV.dat')
    open(245,File='./pionInduced_dOmega_piNull_Delta'//RealToChar(getEkin()*1000.) //'MeV.dat')
    open(246,File='./pionInduced_dOmega_piPlus_Delta'//RealToChar(getEkin()*1000.) //'MeV.dat')
    write(244,'(A,f7.4,A,I4)') '# Final Pi^- distribution. Initial energy:', getEKin(), '. Number of runs:', numberruns
    write(245,'(A,f7.4,A,I4)') '# Final Pi^0 distribution. Initial energy:', getEKin(), '. Number of runs:', numberruns
    write(246,'(A,f7.4,A,I4)') '# Final Pi^+ distribution. Initial energy:', getEKin(), '. Number of runs:', numberruns
    write(244,'(A)') '# theta [degree], dsigma/dOmega [mB/sr], error of dsigma/dOmega[mB/sr]'
    write(245,'(A)') '# theta [degree], dsigma/dOmega [mB/sr], error of dsigma/dOmega[mB/sr]'
    write(246,'(A)') '# theta [degree], dsigma/dOmega [mB/sr], error of dsigma/dOmega[mB/sr]'

    do j=lbound(dsigma_total,dim=2),ubound(dsigma_total,dim=2)
       do k=lbound(dsigma_total,dim=3),ubound(dsigma_total,dim=3)
          if (debug) then
             !             write(*,*) float(j)*delta_theta, sum(dsigma_total(:,j,k)),  error
          end if
          file=245+k
          summe=0.
          do i=lbound(dsigma_total,dim=1),ubound(dsigma_total,dim=1)
             summe=summe+dsigma_total(i,j,k,:)
          end do

          if ((j.eq.lbound(dsigma_total,dim=2)).or.(j.eq.ubound(dsigma_total,dim=2))) then
             write(file,'(62F18.7)') float(j)*delta_theta, summe/float(numberRuns)*delta_ekin*2.
          else
             write(file,'(62F18.7)') float(j)*delta_theta, summe/float(numberRuns)*delta_ekin
          end if


       end do
    end do

    close(244)
    close(245)
    close(246)





    if (finalflag) then
       startFlag=.true.
    end if


  end subroutine dSigmadOmegadE_Resonance



  !****************************************************************************
  !****s* pionXsection/dSigmadOmegadE_pionNuk
  ! NAME
  ! subroutine dSigmadOmegadE_pionNuk(pertParticles,finalFlag)
  ! INPUTS
  ! type(particle), intent(in),dimension(:,:)  :: pertParticles  ! The particle vector
  ! logical,intent(in) :: flag ! Whether it's the final run
  ! OUTPUT
  ! The files:
  ! * pionInduced_dOmegadE_piMinus_piN*MeV.dat
  ! * pionInduced_dOmegadE_piNull_piN*MeV.dat
  ! * pionInduced_dOmegadE_piPlus_piN*MeV.dat
  !
  ! where * is the kinetic energy of the incoming pions.
  ! PURPOSE
  ! This subroutine produces output for pion-nucleus scattering,
  ! including impact parameter integration
  !
  ! Gives dSigma/dOmega/dE(pion) for Pions in file.
  ! Only those pions are considered which are produced in pi N collisions.
  ! NOTES
  ! Makes only sense if useStatistics=.true. in collisionTerm
  !****************************************************************************
  subroutine dSigmadOmegadE_pionNuk(pertParticles,finalFlag)

    use idTable
    use output
    use constants, only: pi, mN, mPi
    use particleDefinition
    use initPion, only: getEkin
    use statistics, only: getInfo
    use collisionNumbering, only:pert_numbering
    use lorentzTrafo, only: lorentz

    type(particle), intent(in),dimension(:,:)  :: pertParticles
    logical, intent(in) :: finalFlag

    integer, parameter :: ekin_dim=50
    integer, parameter :: theta_dim=18

    !real, save :: delta_ekin=0.240/float(ekin_dim)
    !real, save :: delta_theta=180./float(theta_dim)

    real, dimension(0:ekin_Dim,0:Theta_dim,-1:1),save  :: dsigma,dsigma_total, dsigma_squared

    real :: theta, absMom_perpendicular,error,theta_degree,ekin

    integer :: ekinIndex, thetaIndex
    !logical :: first

    integer :: i,j,k
    integer,save :: numberruns=0
    logical, save :: startFlag=.true.

    integer :: ensemble, index,file

    logical :: getInfoFlag
    real, dimension(1:3) :: prodPlace
    integer, dimension(1:3) :: prodParticles

    real,dimension(1:3) :: beta
    real, dimension(0:3) :: mom

    real, save :: delta_ekin   !=0.240/float(ekin_dim)
    real, save :: delta_theta  !=180./float(theta_dim)
    delta_ekin=0.240/float(ekin_dim)
    delta_theta=180./float(theta_dim)


    ! Index 1 : Ekin of final pion, Index 2 : theta of final pion, Index 3: charge of final pion

    if (startFlag) then
       do i=lbound(dsigma,dim=1),ubound(dsigma,dim=1)
          do j=lbound(dsigma,dim=2),ubound(dsigma,dim=2)
             do k=lbound(dsigma,dim=3),ubound(dsigma,dim=3)
                dsigma_total(i,j,k)=0.
                dsigma_squared(i,j,k)=0.
             end do
          end do
       end do
       if (debug) then
          open(22,File='pionAnaTest.txt')
          write(22,*) 'dTheta=', delta_theta
          write(22,*) 'dEkin=', delta_ekin
          close(22)
       end if
       numberruns=0
       !first=.true.
       startFlag=.false.
    end if

    numberruns=numberruns+1

    do i=lbound(dsigma,dim=1),ubound(dsigma,dim=1)
       do j=lbound(dsigma,dim=2),ubound(dsigma,dim=2)
          do k=lbound(dsigma,dim=3),ubound(dsigma,dim=3)
             dsigma(i,j,k)=0.
          end do
       end do
    end do

    ensembleSchleife : do ensemble=1,size(pertParticles,dim=1)
       indexSchleife : do index=1,size(pertParticles,dim=2)

          if (pertParticles(ensemble,index)%ID <= 0) cycle indexSchleife

          if (pertParticles(ensemble,index)%event(1).ne.pert_numbering()) then
             !          if(pertParticles(ensemble,index)%lastCollisionTime.ge.0.0001) then
             if (pertParticles(ensemble,index)%ID.eq.pion) then
                call getInfo(pertParticles(ensemble,index)%number,ensemble,prodPlace,prodParticles,getInfoFlag)
                if (.not.getInfoFlag) then
                   !write(*,*) 'Problem in getting the information with getinfo'
                   cycle indexSchleife
                end if
                if (  (.not.((prodParticles(1).eq.nucleon).and.(prodParticles(2).eq.pion))).and. &
                     (.not.((prodParticles(2).eq.nucleon).and.(prodParticles(1).eq.pion)))     ) then
                   ! Pion not produced in N pion collision
                   cycle indexSchleife
                end if

                mom=pertParticles(ensemble,index)%momentum(0:3)
                ! Evaluating the angular distribution             mom=pertParticles(ensemble,index)%momentum(0:3)
                if (CMFrame) then
                   !Boost momenta from lab frame to cm frame
                   beta=(/0.,0., sqrt((mPi+getekin())**2-mPi**2)/&
                        & (getekin()+mPi+mN)/)
                   call lorentz(beta,mom)
                end if

                absMom_perpendicular=SQRT(mom(1)**2+mom(2)**2)
                if (abs(mom(3)).lt.0.000001) then
                   theta=pi/2.
                   theta_degree=90.
                else if (mom(3).gt.0.) then
                   theta=ATAN(absMom_perpendicular/mom(3))
                   theta_degree=theta/pi*180.
                else
                   theta=pi-ATAN(absMom_perpendicular/Abs(mom(3)))
                   theta_degree=theta/pi*180.
                end if
                if (debug) write(*,*) pertParticles(ensemble,index)%momentum(1:3) , 'theta=',theta


                if (cmFrame) then
                   ekin=sqrt(pertParticles(ensemble,index)%mass**2+mom(1)**2+mom(2)**2+mom(3)**2)&
                        & -pertParticles(ensemble,index)%mass
                else
                   ekin=kineticEnergy(pertParticles(ensemble,index))
                end if

                thetaIndex=Max(Min(180,NINT(theta_degree/delta_Theta)),0)
                ekinIndex=Max(Min(ekin_dim,NINT(ekin/delta_ekin)),0)

                if (abs(sin(theta)).gt.0) then
                   dsigma(ekinIndex,thetaIndex,pertParticles(ensemble,index)%charge)=&
                        & dsigma(ekinIndex,thetaIndex,pertParticles(ensemble,index)%charge) &
                        &       + pertParticles(ensemble,index)%perweight/(delta_Theta*pi/180.*sin(theta)*2.*pi)/delta_ekin
                end if
             end if
          end if
       end do indexSchleife
    end do ensembleSchleife


    ! FOR ERROR ANALYSIS, used for standard deviations among single runs
    do i=lbound(dsigma,dim=1),ubound(dsigma,dim=1)
       do j=lbound(dsigma,dim=2),ubound(dsigma,dim=2)
          do k=lbound(dsigma,dim=3),ubound(dsigma,dim=3)
             dsigma_squared(i,j,k)=dsigma_squared(i,j,k)+dsigma(i,j,k)**2
          end do
       end do
    end do


    ! (2) Add up cross sections of all runs
    do i=lbound(dsigma,dim=1),ubound(dsigma,dim=1)
       do j=lbound(dsigma,dim=2),ubound(dsigma,dim=2)
          dsigma_total(i,j,:)=dsigma_total(i,j,:)+dsigma(i,j,:)
       end do
    end do


    ! (3) Write output

    open(144,File='./pionInduced_dOmegadE_piMinus_piN'//RealToChar(getEkin()*1000.) //'MeV.dat')
    open(145,File='./pionInduced_dOmegadE_piNull_piN'//RealToChar(getEkin()*1000.) //'MeV.dat')
    open(146,File='./pionInduced_dOmegadE_piPlus_piN'//RealToChar(getEkin()*1000.) //'MeV.dat')
    write(144,'(A,f7.4,A,I4)') '# Final Pi^- distribution. Initial energy:', getEKin(), '. Number of runs:', numberruns
    write(145,'(A,f7.4,A,I4)') '# Final Pi^0 distribution. Initial energy:', getEKin(), '. Number of runs:', numberruns
    write(146,'(A,f7.4,A,I4)') '# Final Pi^+ distribution. Initial energy:', getEKin(), '. Number of runs:', numberruns
    write(144,'(A)') '# ePion [GeV], theta [degree], dsigma/dOmega/dE [mB/GeV/sr], error of dsigma/dOmega/dE[mB/GeV/sr]'
    write(145,'(A)') '# ePion [GeV], theta [degree], dsigma/dOmega/dE [mB/GeV/sr], error of dsigma/dOmega/dE[mB/GeV/sr]'
    write(146,'(A)') '# ePion [GeV], theta [degree], dsigma/dOmega/dE [mB/GeV/sr], error of dsigma/dOmega/dE[mB/GeV/sr]'

    do i=lbound(dsigma_total,dim=1),ubound(dsigma_total,dim=1)
       do j=lbound(dsigma_total,dim=2),ubound(dsigma_total,dim=2)
          do k=lbound(dsigma_total,dim=3),ubound(dsigma_total,dim=3)
             if (numberruns.gt.1) then
                ! Evaluate standard deviations of mean value
                error =sqrt(abs(dsigma_squared(i,j,k)-dsigma_total(i,j,k)**2/float(numberruns))/float(numberruns)&
                     & /float(numberruns-1))
             else
                error  =999.
             end if

             if (debug) then
                write(*,*) float(i)*delta_ekin, float(j)*delta_theta, dsigma_total(i,j,k),  error
             end if

             file=145+k
             if ((j.eq.lbound(dsigma_total,dim=2)).or.(j.eq.ubound(dsigma_total,dim=2))) then
                write(file,'(4F18.7)') float(i)*delta_ekin, float(j)*delta_theta, dsigma_total(i,j,k)&
                     & /float(numberRuns)*2., error*2
             else
                write(file,'(4F18.7)') float(i)*delta_ekin, float(j)*delta_theta,&
                     &  dsigma_total(i,j,k)/float(numberRuns), error
             end if
             write(file,*)
          end do
       end do
    end do
    close(144)
    close(145)
    close(146)

    ! Integrated over DE

    open(244,File='./pionInduced_dOmega_piMinus_piN'//RealToChar(getEkin()*1000.) //'MeV.dat')
    open(245,File='./pionInduced_dOmega_piNull_piN'//RealToChar(getEkin()*1000.) //'MeV.dat')
    open(246,File='./pionInduced_dOmega_piPlus_piN'//RealToChar(getEkin()*1000.) //'MeV.dat')
    write(244,'(A,f7.4,A,I4)') '# Final Pi^- distribution. Initial energy:', getEKin(), '. Number of runs:', numberruns
    write(245,'(A,f7.4,A,I4)') '# Final Pi^0 distribution. Initial energy:', getEKin(), '. Number of runs:', numberruns
    write(246,'(A,f7.4,A,I4)') '# Final Pi^+ distribution. Initial energy:', getEKin(), '. Number of runs:', numberruns
    write(244,'(A)') '# theta [degree], dsigma/dOmega [mB/sr], error of dsigma/dOmega[mB/sr]'
    write(245,'(A)') '# theta [degree], dsigma/dOmega [mB/sr], error of dsigma/dOmega[mB/sr]'
    write(246,'(A)') '# theta [degree], dsigma/dOmega [mB/sr], error of dsigma/dOmega[mB/sr]'

    do j=lbound(dsigma_total,dim=2),ubound(dsigma_total,dim=2)
       do k=lbound(dsigma_total,dim=3),ubound(dsigma_total,dim=3)
          if (debug) then
             write(*,*) float(j)*delta_theta, sum(dsigma_total(:,j,k)),  error
          end if
          file=245+k

          if ((j.eq.lbound(dsigma_total,dim=2)).or.(j.eq.ubound(dsigma_total,dim=2))) then
             write(file,'(4F18.7)') float(j)*delta_theta, sum(dsigma_total(:,j,k))/float(numberRuns)*delta_ekin*2.
          else
             write(file,'(4F18.7)') float(j)*delta_theta, sum(dsigma_total(:,j,k))/float(numberRuns)*delta_ekin
          end if



       end do
    end do

    close(244)
    close(245)
    close(246)


    if (finalflag) then
       startFlag=.true.
    end if


  end subroutine dSigmadOmegadE_pionNuk


end module pionXsection
