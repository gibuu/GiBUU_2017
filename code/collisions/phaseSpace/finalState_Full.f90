!******************************************************************************
!****m* /finalState_Full
! NAME
! module finalState_Full
!
! PURPOSE
! This module contains the routines "massAss" and "assMass" which are used to
! choose the momentum and masses in a final state prescription. Read chapter
! 4.7 of Lehr's Dr. thesis to understand the principles of this decision.
!******************************************************************************
module finalState_Full

  implicit none
  private

  public :: massAss_Full
  public :: assMass_Full
  public :: massAss_nBody

  !****************************************************************************
  !****g* finalState_Full/maxbwd_scalingFactor
  ! SOURCE
  !
  real, save :: maxbwd_scalingFactor=1.
  !
  ! PURPOSE
  ! * Rescales maxBWD
  !****************************************************************************


  !****************************************************************************
  !****g* finalState_Full/silentMode
  ! SOURCE
  !
  logical, save :: silentMode=.true.
  !
  ! PURPOSE
  ! * Switches error messages off in massAss. Errors can still be seen
  !   looking at massAssStatus.dat
  !****************************************************************************


  !****************************************************************************
  !****g* finalState_Full/NYK_isotropic
  ! SOURCE
  !
  logical, save :: NYK_isotropic = .false.
  !
  ! PURPOSE
  ! If .true., the angular distribution in Nucleon-Hyperon-Kaon production
  ! is assumed to be isotropic. If .false., a non-isotropic distribution
  ! is used, as described in Larionov/Mosel, Phys.Rev. C 72 (2005) 014901.
  ! See also momenta_in_3Body_BYK.
  !****************************************************************************


  ! Type definition for self-analysis:
  type channelType
     logical       :: init=.false.
     character(40) :: description=''
     real          :: max=0.
     real          :: max_duringRun=0.
     integer       :: numEvents         =0
     integer       :: numEvents_problem =0
     integer       :: numEvents_critical=0
     real          :: value=0.
     real          :: valueSQUARED=0.
     integer       :: numSteps = 0
  end type channelType


  logical, save :: debugFlag=.false.
  logical, save :: init_readinput=.true.

CONTAINS

  !****************************************************************************
  !****s* finalState_Full/readinput
  ! NAME
  ! subroutine readinput
  ! PURPOSE
  ! This subroutine reads in the namelist "FinalState_Full".
  ! Only called once to initialize the module.
  !****************************************************************************
  subroutine readinput
    use output, only: Write_ReadingInput

    integer :: IOS

    !**************************************************************************
    !****n* finalState_Full/FinalState_Full
    ! NAME
    ! NAMELIST /FinalState_Full/
    ! PURPOSE
    ! This namelist includes the following switches:
    ! * maxbwd_scalingFactor
    ! * silentMode
    ! * NYK_isotropic
    !**************************************************************************
    NAMELIST /finalState_Full/ maxbwd_scalingFactor, silentMode, NYK_isotropic

    call Write_ReadingInput('finalState_Full',0)
    rewind(5)
    read(5,nml=finalState_Full,IOSTAT=IOS)
    call Write_ReadingInput('finalState_Full',0,IOS)
    write(*,*) 'Scaling factor for maxBWD:', maxbwd_scalingFactor
    if (silentMode) write(*,*) 'WARNING: Silent mode. Errors are only reported to massAssStatus.dat'
    write(*,*) 'NYK_isotropic: ', NYK_isotropic
    call Write_ReadingInput('finalState_Full',1)

    init_readinput = .false.

  end subroutine readinput

  !****************************************************************************

  function initChannel(descr,m) result(c)
    type(channelType) :: c       ! channel
    character(*)      :: descr   ! description of channel
    real              :: m       ! maximal function value for this channel
    c%init=.true.
    c%description=descr
    c%max=m
    c%max_duringRun=0.
  end function initChannel

  !****************************************************************************

  subroutine protocolChannel(c,m,maxi,pairIn,pairOut,masses,steps)
    use particleDefinition
    type(channelType) :: c    ! channel
    real              :: m    ! evaluated function value for the channel
    real              :: maxi ! maximal value
    type(particle), dimension (1:2), intent(in) :: pairIn,pairOut
    real, intent(in),dimension(1:2)  :: masses
    integer, intent(in) :: steps
    c%max=max(maxi,c%max)
    c%max_duringRun=max(m,c%max_duringRun)
    c%numEvents=c%numEvents+1
    c%value=c%value+m
    c%valueSQUARED=c%valueSQUARED+m**2
    c%numSteps = c%numSteps + steps
    if (m.gt.maxi) then
       c%numEvents_critical=c%numEvents_critical+1
       if (.not. silentMode) then
         write(*,'(A,3G12.3)') "Warning in massAss (protocolChannel): ",m,maxi,m/maxi
         write(*,*) ' Incoming particles :', pairIn(1:2)%ID
         write(*,*) ' Outgoing particles :', pairOut(1:2)%ID
         write(*,*) ' Outgoing masses    :', masses
       end if
    else if (m.gt.maxi*0.8) then
       c%numEvents_problem=c%numEvents_problem+1
    end if
  end subroutine protocolChannel


  !****************************************************************************
  !****s* finalState_Full/massAss_Full
  ! NAME
  ! subroutine massAss_Full (srts, medium_AtCollision, pairIn, pairOut, betaToLRF, betaToCM, L, successFlag)
  !
  ! PURPOSE
  ! This routine selects the masses of a 2 particle final state
  ! according to phase space x spectral functions, simultaneously the
  ! scattering angle is choosen (necessary if spectral functions are
  ! momentum dependent). See J. Lehr Dr. Thesis.
  ! We regard the process a+b->c+d.
  !
  ! INPUTS
  ! * real                           :: srts       -- sqrt(s)
  ! * type(medium)                   :: medium_AtCollision -- medium information : density, temperature,...
  ! * type(particle), dimension(1:2) :: pairIn    -- incoming particles a and b
  ! * type(particle), dimension(1:2) :: pairOut   -- outgoing particles c and d
  !                                                  (only id's, charges and antiflags are used at input)
  ! * real, dimension(1:3)           :: betaToLRF -- beta for boost to LRF
  ! * real, dimension(1:3)           :: betaToCM  -- beta for boost to CM-Frame
  ! * integer                        :: L         -- angular momentum in final state
  !
  ! OUTPUT
  ! * logical :: successFlag -- .true. if mass assignment was successful
  ! * type(particle), dimension(1:2) :: pairOut --
  !   final state particles with full kinematics in the CM frame.
  !   ID, momentum(1:3),charge and mass are now defined.
  !****************************************************************************
  subroutine massAss_Full (srts, medium_AtCollision, pairIn, pairOut, betaToLRF, betaToCM, L, successFlag)

    use IDTable, only: rho,omegaMeson,phi,isMeson,isBaryon
    use mediumDefinition
    use particleDefinition
    use particleProperties, only: hadron
    use lorentzTrafo, only: lorentz
    use MesonWidthMedium, only: WidthMesonMedium, get_MediumSwitchMesons
    use BaryonWidthMedium, only: WidthBaryonMedium, get_MediumSwitch_coll
    use mesonPotentialModule, only: vecMes_massShift

    real, intent(in)                           :: srts
    type(medium), intent(in)                   :: medium_AtCollision
    type(particle), dimension(1:2), intent(in) :: pairIn
    real, dimension(1:3), intent(in)           :: betaToLRF, betaToCM
    integer, intent(in)                        :: L

    type(particle), dimension(1:2), intent(INout) :: pairOut
    logical, intent(out)                          :: successFlag

    integer, dimension(1:2)   :: idOUT,chargeOUT ! id's and charges of outgoing particles
    real, dimension(1:3,1:2)  :: momentum        ! CM momentum of outgoing particles
    real, dimension(1:2)      :: masses          ! masses of outgoing particles
    real, dimension(1:2)      :: spotOUT         ! scalar potential of produced particles
    type(channelType), dimension(1:18), save :: channels
    integer :: channelNumber,i,k
    integer,save :: numCalls=0
    logical :: massFailure_nucleon
    logical, dimension (1:2)  :: flagStable
    real, dimension(1:2)  :: gamma_pole, mass_pole, minMass, maxMass
    real, dimension(0:3) :: pLRF
    real, dimension(1:3) :: betaCMtoLAB
    real :: gamTot, spectral_max, intfac_max

    if (init_readinput) call readinput

    numCalls=numCalls+1

    idOUT(1:2)=pairOut(1:2)%Id
    chargeOUT(1:2)=pairOut(1:2)%charge

    ! set final state potential
    do i=1,2
      select case (idOUT(i))
      case (rho,omegaMeson,phi)
        ! set scalar potential for vector mesons (mass shift!)
        spotOut(i) = vecMes_massShift(idOUT(i),medium_AtCollision%density)
      case default
        ! in all other cases: neglect potential
        spotOut(i) = 0.
      end select
    end do

    call throwDice (successFlag, momentum, masses, massFailure_nucleon, .false.)  ! Monte Carlo decision for final state

    ! If the mass decision failed for a final state including a nucleon, then we assume the nucleon on-shell and try once more:
    if (massFailure_nucleon) call throwDice (successFlag, momentum, masses, massFailure_nucleon, .true.)

    if (successFlag) then
      do i=1,2
        pairOUT(i)%momentum(0) = sqrt((masses(i)+spotOut(i))**2+Dot_Product(momentum(:,i),momentum(:,i)))
        pairOUT(i)%momentum(1:3) = momentum(:,i)
        pairOUT(i)%mass = masses(i)
      end do
    end if

    !    call WriteParticle(6,2,1, pairout(1))
    !    call WriteParticle(6,2,2, pairout(2))
    !    write(*,*) '~~~~~~~~~~'

    ! Write Status of the module to file:
    ! * Which channels have been initialized, how often?,
    !   how often did it fail?, how often was it close to failure?
    if (mod(numCalls,1000)==0) then
       ! Print out status every 1000th time to save system time for opening and closing
       call writeChannels()
    end if

  contains

    !**************************************************************************
    !****s* massAss/writeInput
    ! NAME
    ! subroutine writeInput()
    ! PURPOSE
    ! This subroutine prints the input to screen for debugging!
    !**************************************************************************
    subroutine writeInput()
      use mediumDefinition, only: writeMedium
      use output, only: writeParticle_debug
      write(*,*) 'srts',srts
      call writeParticle_debug(pairIn(1))
      call writeParticle_debug(pairIn(2))
      call writeMedium(medium_AtCollision)
      write(*,*) 'IDs',idOUT
      write(*,*) 'charges',chargeOUT
      write(*,'(A,2E15.5)') 'spot=',spotOUT
      write(*,'(A,3E15.5)') 'betaToLRF=',betaToLRF
      write(*,'(A,3E15.5)') 'betaToCM=',betaToCM
      call writeParticle_debug(pairOUT(1))
      call writeParticle_debug(pairOUT(2))
    end subroutine writeInput

    !**************************************************************************
    !****s* massAss/getMaximum
    ! NAME
    ! subroutine getMaximum(maxBWD)
    ! PURPOSE
    ! This subroutine evaluates the maximal value of
    ! $\frac{p_{cd}}{p^vacuum_{cd}}  \A_c  \A_d   \frac{d\mu_c}{dy_c} \frac{d\mu_d}{dy_d}$
    ! in process a b -> c d. It's a purely empirical value.
    ! RESULT
    ! real :: maxBWD
    !**************************************************************************
    subroutine getMaximum(maxBWD)
      use IdTable, only: nucleon, Delta, nres, nsres, photon, pion, rho, omegaMeson, phi, D35_2350, D35_1930, isBaryonResonance
      use mesonWidthMedium, only: get_MediumSwitchMesons
      use baryonWidthMedium, only: get_MediumSwitch_coll

      real, intent(out)  :: maxBWD
      if (((pairIn(1)%ID.eq.photon).or.(pairIn(2)%ID.eq.photon)) &
           & .and.((idOut(1).eq.delta.and.idOut(2).ne.pion.and.isMeson(idOut(2))) &
           .or.(idOut(2).eq.delta.and.idOut(1).ne.pion.and.isMeson(idOut(1))))) then
         ! gamma N -> Delta meson(no pion)
         maxbwd=12.
         channelNumber=1
         if (.not.channels(channelNumber)%init) &
              & channels(channelNumber)=initChannel('gamma N -> Delta meson(no pion)',maxbwd*maxbwd_scalingFactor)
      else if (get_MediumSwitchMesons() .and. (pairOut(1)%Id==omegaMeson .or. pairOut(2)%Id==omegaMeson)) then
         ! omega with in-medium width
         maxbwd=4.
         channelNumber=2
         if (.not.channels(channelNumber)%init) &
              & channels(channelNumber)=initChannel('omega with in-medium width',maxbwd*maxbwd_scalingFactor)
      else if (get_MediumSwitchMesons() .and. (pairOut(1)%Id==phi .or. pairOut(2)%Id==phi)) then
         ! phi with in-medium width
         maxbwd=4.
         channelNumber=3
         if (.not.channels(channelNumber)%init) &
              & channels(channelNumber)=initChannel('phi with in-medium width',maxbwd*maxbwd_scalingFactor)
      else if (pairIn(1)%Id.ge.nres+4.and.pairIn(1)%Id.le.1+nres+nsres .and. pairIn(2)%Id.eq.0) then
         ! Y^* decay:
         maxbwd=12.
         channelNumber=4
         if (.not.channels(channelNumber)%init) &
              & channels(channelNumber)=initChannel('Y^* decay',maxbwd*maxbwd_scalingFactor)
      else if ((pairOut(1)%Id==rho .and. pairOut(2)%Id==Delta) &
                & .or.(pairOut(1)%Id==Delta .and. pairOut(2)%Id==rho)) then
         ! rho Delta
         maxbwd=20.
         channelNumber=5
         if (.not.channels(channelNumber)%init) &
              & channels(channelNumber)=initChannel('rho and Delta',maxbwd*maxbwd_scalingFactor)
      else if (pairIn(1)%Id==D35_2350 .and. pairIn(2)%Id==0) then
         maxbwd=20.
         channelNumber=6
         if (.not.channels(channelNumber)%init) &
              & channels(channelNumber)=initChannel('D35_2350 decay',maxbwd*maxbwd_scalingFactor)
      else if (pairIn(1)%Id==D35_1930 .and. pairIn(2)%Id==0) then
         maxbwd=20.
         channelNumber=7
         if (.not.channels(channelNumber)%init) &
              & channels(channelNumber)=initChannel('D35_1930 decay',maxbwd*maxbwd_scalingFactor)
      else if (get_MediumSwitch_coll().and. isBaryon(pairOut(1)%ID).and. isBaryon(pairOut(2)%ID).and. &
                  & (pairOut(1)%Id.eq.nucleon.or.pairOut(2)%Id.eq.nucleon)) then
         ! Collisional width: 2 Baryon final states with nucleon in outgoing channel
         if (pairOut(1)%Id+pairOut(2)%Id.gt.2) then
            ! N-Resonance or Resonance-Resonance scattering for offshell nucleons
            maxbwd=50
            channelNumber=8
            if (.not.channels(channelNumber)%init) &
                 & channels(channelNumber)=initChannel('coll width: NR or RR',maxbwd*maxbwd_scalingFactor)
         else
            ! NN scattering for offshell nucleons
            maxbwd=100.
            channelNumber=9
            if (.not.channels(channelNumber)%init) &
                 & channels(channelNumber)=initChannel('coll width: NN',maxbwd*maxbwd_scalingFactor)
         end if
      else if (get_MediumSwitch_coll().and. (pairOut(1)%Id==nucleon .or. pairOut(2)%Id==nucleon)) then
         ! Collisional width: Nucleon meson final states
         maxbwd=40.
         channelNumber=10
         if (.not.channels(channelNumber)%init) &
              & channels(channelNumber)=initChannel('coll width: N mes',maxbwd*maxbwd_scalingFactor)
      else if (pairIn(1)%Id>0 .and. pairIn(1)%Id<=nres+1 .and. pairIn(2)%Id==0) then
         ! decay of nucleon resonances
         maxbwd=6.
         channelNumber=11
         if (.not.channels(channelNumber)%init) &
              & channels(channelNumber)=initChannel('Res dec',maxbwd*maxbwd_scalingFactor)
      else if (.not.(medium_AtCollision%useMedium .and. get_MediumSwitchMesons()) .and. &
                ((flagStable(2) .and. (IdOut(1)==rho .or. IdOut(1)==omegaMeson .or. IdOut(1)==phi)) .or. &
                 (flagStable(1) .and. (IdOut(2)==rho .or. IdOut(2)==omegaMeson .or. IdOut(2)==phi)))) then
         ! vector meson (with vacuum width) + stable particle (e.g. omega-N, phi-N)
         if (flagstable(1)) then
           k=2
         else if (flagstable(2)) then
           k=1
         end if
         pLRF(1:3)=0.
         pLRF(0)=maxMass(k)+spotOut(k)
         betaCMtoLab=-betaToCM
         call lorentz(betaCMtoLab(1:3), plrf(0:3), 'finalState(1)')         ! boost from CM to Lab
         call lorentz(betaToLRF(1:3), plrf(0:3), 'finalState(2)')         ! boost from Lab to LRF
         ! calculate gamtot, spectral and intfac for mass(k)=maxMass(k)
         gamtot=WidthMesonMedium(IDOut(k),maxMass(k),plrf(0:3),medium_ATCollision)
         spectral_max=maxMass(k)**2*gamtot*gamma_pole(k)/ ((mass_pole(k)**2-maxMass(k)**2)**2 + gamtot**2*maxMass(k)**2)
         intfac_max=gamma_pole(k)**2/ ((maxMass(k)-mass_pole(k))**2+gamma_pole(k)**2/4.)/4.
         ! assumption: "spectral/intfac" is maximal for largest possible mass
         maxbwd=spectral_max/intfac_max*1.05
         channelNumber=12
         if (.not.channels(channelNumber)%init) &
              & channels(channelNumber)=initChannel('VM + stable in vac.',maxbwd*maxbwd_scalingFactor)
      else if (((isMeson(pairIn(1)%Id).and.isBaryon(pairIn(2)%Id)) .or. (isMeson(pairIn(2)%Id).and.isBaryon(pairIn(1)%Id))) &
               .and. (idOut(1)==phi.or.idOut(2)==phi)) then
         ! baryon meson -> baryon phi
         maxbwd = 24.
         channelNumber = 13
         if (.not.channels(channelNumber)%init) &
           channels(channelNumber)=initChannel('baryon meson -> baryon phi ',maxbwd*maxbwd_scalingFactor)
      else if (pairIN(1)%ID==nucleon .and. pairIN(2)%ID==nucleon .and. &
               pairOUT(1)%ID==nucleon .and. isBaryonResonance(pairOUT(2)%ID)) then
         ! NN -> N Res
         maxbwd = 6.
         channelNumber = 14
         if (.not.channels(channelNumber)%init) &
           channels(channelNumber) = initChannel('NN -> NR',maxbwd*maxbwd_scalingFactor)
      else if (pairIN(1)%ID==nucleon .and. pairIN(2)%ID==nucleon .and. &
               pairOUT(1)%ID==Delta .and. isBaryonResonance(pairOUT(2)%ID)) then
         ! NN -> Delta Res
         maxbwd = 36.
         channelNumber = 15
         if (.not.channels(channelNumber)%init) &
           channels(channelNumber) = initChannel('NN -> DR',maxbwd*maxbwd_scalingFactor)
         ! K K* -> rho omega
      else if (min(pairIN(1)%ID,pairIN(2)%ID)==110 .and. &
               max(pairIN(1)%ID,pairIN(2)%ID)==111 .and. &
               min(pairOUT(1)%ID,pairOUT(2)%ID)==103 .and. &
               max(pairOUT(1)%ID,pairOUT(2)%ID)==105 ) then
         maxbwd = 70.
         channelNumber = 16
         if (.not.channels(channelNumber)%init) &
           channels(channelNumber) = initChannel('KK* -> rho omega',maxbwd*maxbwd_scalingFactor)
         ! K K* -> omega omega
      else if (min(pairIN(1)%ID,pairIN(2)%ID)==110 .and. &
               max(pairIN(1)%ID,pairIN(2)%ID)==111 .and. &
               min(pairOUT(1)%ID,pairOUT(2)%ID)==105 .and. &
               max(pairOUT(1)%ID,pairOUT(2)%ID)==105 ) then
         maxbwd = 40.
         channelNumber = 17
         if (.not.channels(channelNumber)%init) &
           channels(channelNumber) = initChannel('KK* -> omega omega',maxbwd*maxbwd_scalingFactor)

      else
         ! default
         maxbwd = 4.
         channelNumber = 18
         if (.not.channels(channelNumber)%init) &
           channels(channelNumber)=initChannel('Default',maxbwd*maxbwd_scalingFactor)
      end if
      maxbwd=maxbwd*maxbwd_scalingFactor
    end subroutine getMaximum

    !**************************************************************************
    !****s* massAss/checkMaximum
    ! NAME
    ! subroutine checkMaximum(value, maxValue)
    ! PURPOSE
    ! This subroutine checks wether the maximal value is really the maximal
    ! value.
    ! If not, then error messages are written to standard out and
    ! 'massAssError.dat'. After a critical number of errors the program is
    ! terminated!!
    !**************************************************************************
    subroutine checkMaximum(value, maxValue)
      use baryonWidthMedium, only: get_MediumSwitch_coll
      use CallStack, only: Traceback

      real, intent(in) :: value, maxValue
      integer, save :: numberOfErrors = 0
      integer, parameter :: MaxNumberOfErrors = 100

      if (value>maxValue) then
         if (.not. silentMode) then
            write(*,*) 'Warning in massass: Value > maxValue.'
            write(*,*) 'This causes problem in Monte-Carlo decision. Modify Maxvalues!!!'
            write(*,*) ' Incoming particles :', pairIn(1:2)%ID
            write(*,*) ' Outgoing particles :', IDOut(1:2)

            open(243,File='massAssError.dat',position='Append',status='unknown')
            write(243,*) 'Warning in massass: Value > maxValue.'
            write(243,*) 'This causes problem in Monte-Carlo decision. Modify Maxvalues!!!'
            write(243,*) ' Incoming particles :', pairIn(1:2)%ID
            write(243,*) ' Outgoing particles :', IDOut(1:2)
            write(243,*) 'Value=',Value, 'MaxBwd=',maxValue
            write(243,*) 'for offshell particles: medium_AtCollision%density',medium_AtCollision%density
            write(243,*) '*****************************************************************************'
            close(243)
         end if
         numberOfErrors=numberOfErrors+1
         if ((.not. get_MediumSwitch_coll()) .and. (MaxNumberOfErrors<numberOfErrors)) then
            write(*,*) 'channel: ', channelNumber
            call Traceback('This problem occured now too often! Stop')
         end if
      end if
    end subroutine checkMaximum

    !**************************************************************************
    !****s* massAss/throwDice
    ! NAME
    ! subroutine throwDice (success, momentum, mass, massFailure_Nuc, treatNuc_onshell)
    ! PURPOSE
    ! We consider a + b -> c + d. The kinematics of c and d are determined.
    ! Here the masses and momenta of the particles are choosen by Monte-Carlo
    ! decision.
    ! Then p_cd/p_cd_vacuumMax*spectral(c)*spectral(d)/(intfac(c)*intfac(d)*maxValue)
    ! is evaluated which our probability to accept a choosen final state.
    ! We then decide to reject or accept the final state.
    ! If rejected new masses and momenta are choosen until we find a final
    ! state that is accepted.
    ! If it's not possible to find a final state, than successs=.false. is set
    ! and the routine is terminated.
    ! INPUTS
    ! * real :: maxValue --
    !   maximum of the function
    ! * logical :: treatNuc_onshell  --
    !   true: set nucleon masses to onshell value
    !
    ! OUTPUT
    ! * real, dimension(1:3,1:2) :: momentum -- of the final state
    !   (first index: momentum, second: particle)
    ! * real, dimension(1:2)     :: mass     -- of the final state
    ! * logical                  :: success  --
    !   .true. if mass and momentum could be set
    ! * logical                  :: massFailure_nuc --
    !   .true. if process failed when setting the mass of a final state
    !   involving a nucleon
    !**************************************************************************
    subroutine throwDice (success, momentum, mass, massFailure_Nuc, treatNuc_onshell)
      use IdTable
      use random, only: rn
      use winkelVerteilung, only: winkel
      use mesonWidthMedium, only: WidthMesonMedium
      use twoBodyTools, only: pCM, pCM_sqr
      use CALLSTACK, only: TRACEBACK
      use constants, only: rhoNull
      use distributions, only: BlattWeisskopf
      use baryonWidthVacuum, only: interactionRadius

      real, intent(out),dimension(1:2) :: mass
      real, dimension(1:3,1:2),intent(out) :: momentum
      logical, intent(out) :: success, massFailure_nuc
      logical, intent(in)  :: treatNuc_onshell

      real :: maxValue,probability,p_cd_vacuumMax,p_cd,helper,blw,blw_max
      real, dimension(1:2) :: yMax, yMin, y
      real, dimension(1:2) :: Spectral, intfac   ! spectral functions
      real, dimension(1:3) :: pscatt
      integer, parameter :: maxMassSteps=60000
      integer, parameter :: maxMomSteps=300000
      integer :: momSteps, massSteps
      integer, save :: err_count_mass=0, err_count_mom=0, err_count_mass_nucleon=0
      logical :: winkelflag

      success=.false.
      massFailure_nuc=.false.

      do k=1,2
         select case (idOut(k))
         case (nucleon)
            if (get_MediumSwitch_coll()) then
               if (treatNuc_onshell) then
                  gamma_Pole(k)=0.
               else
                  gamma_Pole(k)=0.035*medium_AtCollision%density/rhoNull  ! 35 MeV is empirical value!!!
               end if
            else
               gamma_Pole(k) = hadron(idOut(k))%width
            end if
            mass_Pole(k)  = hadron(idOut(k))%mass
            minMass(k)    = hadron(idOut(k))%minmass

         case (Delta:nbar)
            gamma_Pole(k) = hadron(idOut(k))%width
            mass_Pole(k)  = hadron(idOut(k))%mass
            minMass(k)    = hadron(idOut(k))%minmass

         case (pion:eta,sigmaMeson,etaPrime,etaC:dSStar_minus)
            mass_Pole(k)  = hadron(idOut(k))%mass
            gamma_Pole(k) = WidthMesonMedium(idOut(k),hadron(idOut(k))%mass,(/0.,0.,0.,0./),medium_AtCollision)
            minMass(k)    = hadron(idOut(k))%minmass

         case (rho,omegaMeson,phi)
            mass_Pole(k)  = hadron(idOut(k))%mass
            gamma_Pole(k) = WidthMesonMedium(idOut(k),hadron(idOut(k))%mass,(/0.,0.,0.,0./),medium_AtCollision)
            if (get_MediumSwitchMesons() .and. medium_AtCollision%useMedium) then
               minMass(k)    = - spotOut(k) + 0.01
            else
               minMass(k)    = hadron(idOut(k))%minmass
            end if

         case (photon)
            gamma_Pole(k) = 0.
            mass_Pole(k)  = 0.

         case default
            write(*,*) 'strange particle ID in massAss:',k,IDOut
            call TRACEBACK('strange particle ID in massAss')
         end select

         ! Check whether particles are regarded as stable
         if (gamma_pole(k) < 1e-03) then
            flagStable(k)=.true.
            minmass(k)=mass_Pole(k)
         else
            flagStable(k)=.false.
         end if
      end do

      if (debugFlag) write(*,*) flagStable, IdOut, minMass

      ! Check kinematics
      if (sum(MinMass) > srts) then
         success=.false.
         return
      end if

      ! Do variable transformation mass-> y
      do k=1,2
         if (.not. flagStable(k)) then

            if (IdOUT(k)==nucleon) then
              maxmass(k)=min(1.5,srts-minmass(-k+3)-spotOut(k))
            else
              maxmass(k)=srts-minmass(-k+3)-spotOut(k)
            end if

            if (maxmass(k) < minmass(k)) then
               write(*,*) 'problems in massass maxmass.lt.minmass',    srts,maxmass(k),minmass(k)
               write(*,*) idOUT
               write(*,*) spotOUT
               stop
            end if
            ymax(k) = 2.*atan((maxmass(k)-mass_pole(k))/gamma_pole(k)*2.)
            ymin(k) = 2.*atan((minmass(k)-mass_pole(k))/gamma_pole(k)*2.)
         else
            maxmass(k)=mass_pole(k)
         end if
      end do

      call getMaximum(maxValue)


      ! MONTE CARLO STARTS :::

      if (.not.(flagstable(1).and.flagstable(2))) then
            ! Evaluate maximal momentum of particles c and d in CM-frame (in the vacuum, but with final state potentials)
         p_cd_vacuumMax = pCM_sqr (srts**2, (minmass(1)+spotOut(1))**2, (minmass(2)+spotOut(2))**2)
         if (p_cd_vacuumMax<0.) then
            write(*,*) 'problems in massass s too low',srts,minmass
            stop
         end if
         p_cd_vacuumMax = sqrt(p_cd_vacuumMax)
         blw_max = BlattWeisskopf(p_cd_vacuumMax*interactionRadius, L)**2

         ! Start monte carlo to find momenta and masses
         momLoop : do momSteps=1,maxMomSteps
            ! Monte carlo to find the masses
            massLoop : do massSteps=1,maxMassSteps
               ! Loop is necessary because both masses have to be chosen simultaneously*
               do k=1,2
                  if (.not. flagStable(k)) then
                     y(k)=ymin(k)+rn()*(ymax(k)-ymin(k))
                     mass(k)=.5*tan(y(k)/2.)*gamma_pole(k)+ mass_pole(k)
                  else
                     mass(k)=mass_pole(k)
                  end if
               end do
               if (Sum(mass)+Sum(spotOut)<srts) exit massLoop
            end do massLoop
            if (massSteps >= maxMassSteps) then
               !Failure of the algorithm, no valid solution for the masses could be found
               success=.false.
               if ((IdOUT(1)==nucleon .and. .not.flagstable(1)) .or. (IdOUT(2)==nucleon .and. .not.flagstable(2))) then
                  write(*,*) 'Warning in massAss: Mass failure with broad nucleon in final state!', err_count_mass_nucleon
                  if (.not. silentMode) then
                     write(*,*) ' Minimal Masses', minMass
                     write(*,*) ' Maximal Masses', maxMass
                     call writeInput()
                  end if
                  err_count_mass_nucleon=err_count_mass_nucleon+1
                  massFailure_nuc=.true.
                  if (err_count_mass_nucleon > 10000) then
                     stop 'massAss, Mass iteration for nucleon final state: more than 10.000 errors!!'
                  end if
               else
                  err_count_mass=err_count_mass+1
                  write(*,*) 'Warning : Mass Iteration in massAss failed'
                  write(*,*) ' Minimal Masses', minMass
                  write(*,*) ' Maximal Masses', maxMass
                  call writeInput()
                  if (err_count_mass > 10000) then
                     stop 'massAss, Mass iteration'
                  end if
               end if
               return
            end if

            ! cm momentum of produced particles
            p_cd = pCM(srts,mass(1)+spotOut(1),mass(2)+spotOut(2))
            if (p_cd/p_cd_vacuumMax>=1.) then
               write(*,*) "Problem with p_cd's in massAss", p_cd,p_cd_vacuumMax,srts,idOut(1:2),mass(1:2),minMass(1:2),spotOut(1:2)
            end if

            ! determine scattering angle
            pairOut(1:2)%mass=mass(1:2)
            pscatt = winkel (pairIN, pairOut, srts, betaToCM, medium_AtCollision, winkelflag)
            if (.not. winkelflag) then
               if (get_MediumSwitch_coll() .or. get_MediumSwitchMesons()) then
                  ! kinematics doesn't allow this process
                  success=.false.
                  return
               else
                  write(*,*) 'massAss -> winkelflag=.false. but mediumSwitch=.false. -> should not happen -> stop'
                  write(*,*) pairIN(1),pairIN(2),get_MediumSwitch_coll(),get_MediumSwitchMesons()
                  stop
               end if
            end if

            ! Calculate spectral functions
            do k=1,2
               if (.not. flagStable(k)) then
                  ! Determine momenta in LRF for evaluation of the width
                  pLRF(1:3)=(-1)**(k+1)*pscatt(:)*p_cd
                  plrf(0)=sqrt((mass(k)+spotOut(k))**2+plrf(1)**2+ plrf(2)**2+plrf(3)**2)
                  betaCMtoLab=-betaToCM
                  call lorentz(betaCMtoLab(1:3), plrf(0:3), 'finalState(1)')         ! boost from CM to Lab
                  call lorentz(betaToLRF(1:3), plrf(0:3), 'finalState(2)')         ! boost from Lab to LRF

                  if (isMeson(idOut(k))) then
                     gamtot=WidthMesonMedium(IDOut(k),mass(k), plrf(0:3) ,medium_ATCollision)
                  else
                     gamtot=WidthBaryonMedium(IDOut(k),mass(k),plrf(0:3) ,medium_ATCollision)
                  end if
                  spectral(k)=mass(k)**2*gamtot*gamma_pole(k)/ ((mass_pole(k)**2-mass(k)**2)**2 + gamtot**2 *mass(k)**2)
                  intfac(k)=gamma_pole(k)**2/ ((mass(k)-mass_pole(k))**2+gamma_pole(k)**2/4.)/4.
               else
                  spectral(k)=1.
                  intfac(k)=1.
               end if
            end do

            blw = BlattWeisskopf(p_cd*interactionRadius, L)**2

            helper = spectral(1) * spectral(2) / intfac(1) / intfac(2)
            probability = p_cd/p_cd_vacuumMax * blw/blw_max * helper / maxValue
            if (probability>=1.) call checkMaximum(helper, maxValue)

            if (probability>rn()) then
               ! Success
               success=.true.
               exit momLoop
            end if
         end do momLoop

         ! Protocol the maximum
         call protocolChannel(channels(channelNumber),helper,maxValue,pairIn,pairOut,mass,momSteps)

         if (momsteps < maxMomSteps) then
            ! Success
            momentum(:,1)=pscatt*p_cd
            momentum(:,2)=-pscatt*p_cd
            success=.true.
            if (debugFlag) write(*,*) momSteps
         else
            err_count_mom=err_count_mom+1
            write(*,*) 'Momentum iteration in massass failed. Momsteps=', momsteps
            write(*,*) IDout, chargeOUT
            write(*,*) SPotOut, srts, minmass
            write(*,*) probability, helper, blw, spectral(1:2), intfac(1:2), p_cd, p_cd_vacuumMax
            write(*,*) 'MaxValue=',maxValue
            write(*,*) pairIn(1:2)%ID, L
            success=.false.
            if (err_count_mom > 10000) then
               stop 'In massass, momentum iteration'
            end if
            return
         end if

      else ! stable particles in final state

         mass(1)=mass_pole(1)
         mass(2)=mass_pole(2)
         p_cd = pCM_sqr (srts**2, (mass(1)+spotOut(1))**2, (mass(2)+spotOut(2))**2)
         if (p_cd<0.) then
            ! failure : particles can't be produced with this srts
            success=.false.
            return
         end if

         p_cd = sqrt (p_cd)
         success = .true.

         ! determine scattering angle
         pairOut(1:2)%mass=mass(1:2)
         pscatt = winkel (pairIN, pairOut, srts, betaToCM, medium_AtCollision)

         momentum(:,1)=pscatt*p_cd
         momentum(:,2)=-pscatt*p_cd
      end if

    end subroutine throwDice


    subroutine writeChannels()
      integer, parameter :: iFile = 98
      integer :: n,events_tot,steps_tot
      type(channelType) :: c
      open(iFile,file='massAssStatus.dat')
      write(iFile,'(A20,8A12)')    '____desription____', '___max___', 'max_in_Run','num_Events' &
           & ,'failures[%]','large[%]','___mean___','deviation', 'num_steps'
      events_tot = 0
      steps_tot = 0
      do n=lbound(channels,dim=1),ubound(channels,dim=1)
        c = channels(n)
        if (c%init .and. c%numEvents>0) then
          write(iFile,'(A20,8G12.5)')    c%description, c%max, c%max_duringRun, c%numEvents,&
               float(c%numEvents_critical)/float(c%numEvents)*100., float(c%numEvents_problem )/float(c%numEvents)*100., &
               c%value/float(c%numEvents), &
               sqrt(max(0.,1./float(c%numEvents)*(c%valueSquared-c%value**2/float(c%numEvents)))), &
               float(c%numSteps)/c%numEvents
          events_tot = events_tot + c%numEvents
          steps_tot = steps_tot + c%numSteps
        end if
      end do
      write(iFile,*)
      write(iFile,'(A35,G12.5)') 'Total number of events: ',events_tot
      if (events_tot.gt.0) write(iFile,'(A35,G12.5)') 'Overall average number of steps: ',float(steps_tot)/events_tot
      close(iFile)
    end subroutine writeChannels

  end subroutine massAss_Full

  !****************************************************************************
  !****s* finalState_Full/assMass_Full
  ! NAME
  ! subroutine assMass_Full (srts, medium_AtCollision, pairIn, tripleOut, spotOut, betaToLRF, betaToCM, flag)
  !
  ! PURPOSE
  ! This routine selects the masses of a 3 particle final state
  ! according to phase space x spectral functions, simultaneously the
  ! scattering angle is choosen (necessary if spectral functions are
  ! momentum dependent). See J. Lehr Dr. Thesis.
  ! We regard the process a+b->c+d+e.
  !
  ! INPUTS
  ! * real                           :: srts       ! sqrt(s)
  ! * type(medium)                   :: medium_AtCollision ! medium information : density, temperature,...
  ! * type(particle), dimension(1:2) :: pairIn     ! incoming particles a and b
  ! * type(particle), dimension(1:3) :: tripleOut  ! outgoing particles c,d and e
  ! *                                              ! (only id's, charges and antiflags are used at input)
  ! * real, dimension(1:3)           :: spotOut    ! scalar potential of produced particles

  ! * real, dimension(1:3)           :: betaToLRF  ! beta for boost to LRF
  ! * real, dimension(1:3)           :: betaToCM   ! beta for boost to CM-Frame
  !
  ! RESULT
  ! * logical :: flag         ! set to .true. is mass assignment was successful
  ! * type(particle), dimension (1:3) ,intent(out) :: tripleOut  !
  !   final state particles with almost full kinematics in CM frame.
  !   ID, momentum(1:3),charge and mass are defined.
  !****************************************************************************
  subroutine assMass_Full (srts, medium_AtCollision, pairIn, tripleOut, spotOut, betaToLRF, betaToCM, flag)
    use mediumDefinition
    use particleDefinition

    real, intent(in)                              :: srts
    type(medium), intent(in)                      :: medium_AtCollision
    type(particle), dimension(1:2), intent(in)    :: pairIn
    real, dimension(1:3), intent(in)              :: spotOut
    real, dimension(1:3), intent(in)              :: betaToLRF, betaToCM

    type(particle), dimension(1:3), intent(INout) :: tripleOut
    logical, intent(out)                          :: flag

    real :: maxValue
    integer, dimension(1:3) :: idOUT,chargeOUT ! Id's and charges of outgoing particles
    real, dimension(1:3,1:3) :: momentum  ! CM momentum of outgoing particles
    real, dimension(1:3) :: masses ! masses of outgoing particles
    integer :: i

    if (init_readinput) call readinput

    idOUT(1:3)     = tripleOut(1:3)%Id
    chargeOUT(1:3) = tripleOut(1:3)%charge

    call getMaximum (maxValue)

    call throwDice (maxValue, flag, momentum, masses)                ! Monte Carlo decision for final state

    if (flag) then
      do i=1,3
        tripleOUT(i)%momentum(0) = sqrt((masses(i)+spotOut(i))**2+Dot_Product(momentum(:,i),momentum(:,i)))
        tripleOUT(i)%momentum(1:3) = momentum(:,i)
        tripleOUT(i)%mass = masses(i)
      end do
    end if

  contains

    !**************************************************************************
    !****s* assMass/getMaximum
    ! NAME
    ! subroutine getMaximum(maxBWD)
    ! PURPOSE
    ! This subroutine evaluates the maximal value for the phase space decision
    ! variable.
    ! RESULT
    ! real :: maxBWD
    !**************************************************************************
    subroutine getMaximum(maxBWD)
      use IdTable, only: nucleon, photon, rho, omegaMeson, phi
      use baryonWidthMedium, only: get_MediumSwitch_coll
      real, intent(out)  :: maxBWD
      if ((pairIn(1)%ID==photon .or. pairIn(2)%ID==photon) .and. &
          (idOut(1)==rho .or. idOut(1)==omegaMeson .or. idOut(1)==phi .or. &
           idOut(2)==rho .or. idOut(2)==omegaMeson .or. idOut(2)==phi .or. &
           idOut(3)==rho .or. idOut(3)==omegaMeson .or. idOut(3)==phi)) then
         ! gamma N -> Vectormeson + X
         maxbwd = 2.
      else if (get_MediumSwitch_coll() .and. (idOut(1)==nucleon .or. idOut(2)==nucleon .or. idOut(3)==nucleon)) then
         ! offshell nucleons
         maxbwd = 100. * maxbwd_scalingFactor
      else
         maxbwd = 8.
      end if
    end subroutine getMaximum

    !**************************************************************************
    !****s* assMass/checkMaximum
    ! NAME
    ! subroutine checkMaximum(value, maxValue)
    ! PURPOSE
    ! This subroutine checks wether the maximal value is really the maximal
    ! value.
    ! If not, then error messages are written to standard out and
    ! 'massAssError.dat'. After a critical number of errors the program is
    ! terminated!!
    !**************************************************************************
    subroutine checkMaximum(value, maxValue)
      use baryonWidthMedium, only: get_MediumSwitch_coll

      real, intent(in) :: value, maxValue
      integer, save :: numberOfErrors=0
      integer, parameter :: MaxNumberOfErrors=100

      if (value>maxValue) then
         write(*,*) 'Warning in assMass: Value > maxValue.'
         write(*,*) 'This causes problem in Monte-Carlo decision. Modify Maxvalues!!!'
         write(*,*) ' Incoming particles :', pairIn(1:2)%ID
         write(*,*) ' Outgoing particles :', IDOut(1:3)

         open(243,File='assMassError.dat',position='Append',status='unknown')
         write(243,*) 'Warning in assMass: Value > maxValue.'
         write(243,*) 'This causes problem in Monte-Carlo decision. Modify Maxvalues!!!'
         write(243,*) ' Incoming particles :', pairIn(1:2)%ID
         write(243,*) ' Outgoing particles :', IDOut(1:3)
         write(243,*) 'Value=',Value, 'MaxBwd=',maxValue
         write(243,*) '*****************************************************************************'
         close(243)
         numberOfErrors=numberOfErrors+1
         if ((.not.get_MediumSwitch_coll()) .and. (MaxNumberOfErrors<numberOfErrors)) then
            write(*,*) 'assMass/checkMaximum: This problem occured now too often! Stop'
            stop
         end if
      end if
    end subroutine checkMaximum

    !**************************************************************************
    !****s* assMass/throwDice
    ! NAME
    ! subroutine throwDice (maxValue, success, momentum, mass)
    ! PURPOSE
    ! We consider a + b -> c + d + e. The kinematics of c, d & e are
    ! determined.
    ! Here the masses and momenta of the particles are choosen by Monte-Carlo
    ! decision.
    ! pcmAbs(1)/ pcmAbs(3)/ecm(1)/eCM(3) *spectral(c)*spectral(d)*spectral(d)/ intfac(c) /intfac(d)/intfac(e) /maxValue is
    ! evaluated which is our probability to accept a chosen final state.
    ! We then decide to reject or accept the final state. If rejected new
    ! masses and momenta are choosen until we find a final state that is
    ! accepted.
    ! If it is not possible to find a final state, than successs=.false. is set
    ! and the routine is terminated.
    !**************************************************************************
    subroutine throwDice (maxValue, success, momentum, mass)
      use IdTable, only: photon, nucleon, Delta, nbar, pion, nmes, Lambda, SigmaResonance, Kaon, isMeson
      use particleProperties, only: hadron
      use MesonWidthMedium, only: WidthMesonMedium
      use BaryonWidthMedium, only: WidthBaryonMedium, get_MediumSwitch_coll
      use nBodyPhaseSpace, only: momenta_in_3BodyPS, momenta_in_3Body_BYK
      use random, only: rn
      use lorentzTrafo, only: lorentz
      use constants, only: pi,rhoNull
      use rotation, only: rotateYZ
      use CALLSTACK, only: TRACEBACK

      real, intent(in) :: maxValue
      real, intent(out),dimension(1:3) :: mass
      real, dimension(1:3,1:3),intent(out) :: momentum  ! first index: momentum index, second: particle
      logical, intent(out) :: success

      real :: probability, gamTot, cosTheta_13, sinTheta_13 , phi, cosTheta
      real, dimension(1:3)  :: Spectral, gamma_pole, mass_pole, minMass, maxMass, intfac
      logical, dimension (1:3) :: flagStable
      real, dimension(1:3) :: yMax, yMin, y, pcmAbs, ECM
      real, dimension(0:3) :: pLRF, pcm
      real, dimension(1:3) :: betaCMtoLab
      integer :: k, massSteps, momSteps, outerSteps
      integer, parameter :: maxMassSteps=1000
      integer, parameter :: maxMomSteps=100000
      integer, parameter :: maxOuterSteps=60000

      do k=1,3
         select case (idOut(k))
         case (nucleon)
            if (get_MediumSwitch_coll()) then
               gamma_pole(k) =0.035*medium_AtCollision%density/rhoNull   ! 35 MeV is empirical value!!!
            else
               gamma_pole(k) = hadron(idOut(k))%width
            end if
            mass_pole(k)  = hadron(idOut(k))%mass
            minMass(k)    = hadron(idOut(k))%minmass
         case (Delta:nbar,pion:pion+nMes-1)
            gamma_pole(k) = hadron(idOut(k))%width
            mass_pole(k)  = hadron(idOut(k))%mass
            minMass(k)    = hadron(idOut(k))%minmass
         case (photon)
            gamma_pole(k) =0.
            mass_pole(k)  =0.

         case default
            write(*,*) k,IDOut
            call TRACEBACK('strange particle ID in massAss')
         end select

         ! Check whether particles are regarded as stable
         if (gamma_pole(k)<1e-03) then
            flagStable(k)=.true.
            minmass(k)=mass_Pole(k)
         else
            flagStable(k)=.false.
         end if
      end do

      if (flagstable(1).and.flagstable(2).and.flagstable(3)) then

         ! simplest case: all particles are stable
         mass(1:3) = mass_pole(1:3)
         if (sum(mass) > srts) then
            success = .false.  ! Not possible by kinematics
         else
            if ((.not. NYK_isotropic) .and. &
                IdOut(1)==nucleon .and. IdOut(3)==Kaon .and. &
                (IdOut(2)==Lambda  .or. IdOut(2)==SigmaResonance)) then
              pcm = pairIn(1)%momentum
              call lorentz(betaToCM, pcm) ! boost from Lab to CM
              momentum = momenta_in_3Body_BYK (srts, pcm(1:3), mass(1:3))
            else
              ! set momenta in phase space by random
              momentum = momenta_in_3BodyPS (srts, mass(1:3))
            end if
            success = .true.
         end if

      else

         ! Do variable transformation mass-> y
         do k=1,3
            if (.not.flagStable(k)) then
              maxmass(k)=srts-spotOut(k)-minmass(mod(k,3)+1)-minmass(mod(k+1,3)+1)
              if (idout(k)==nucleon) maxmass(k) = min(1.5,maxmass(k))

              if (maxmass(k)<minmass(k)) then
                write(*,*) 'problems in assMass maxmass.lt.minmass',    srts,maxmass(k),minmass(k)
                write(*,*) idOUT
                write(*,*) spotOUT
                success=.false.
                return
              end if
              ymax(k)=2.*atan((maxmass(k)-mass_pole(k))/gamma_pole(k)*2.)
              ymin(k)=2.*atan((minmass(k)-mass_pole(k))/gamma_pole(k)*2.)
            else
              maxmass(k)=mass_pole(k)
            end if
         end do

         ! Monte Carlo start :
         outerloop : do outerSteps=1,maxOuterSteps

            momLoop : do MomSteps=1,maxMomSteps
               ! Choose masses by choosing y according to a flat distribution
               massLoop : do MassSteps=1,maxMassSteps
                  do k=1,3
                     if (.not.flagStable(k)) then
                        y(k)=ymin(k)+rn()*(ymax(k)-ymin(k))
                        mass(k)=.5*tan(y(k)/2.)*gamma_pole(k)+ mass_pole(k)
                     else
                        mass(k)=mass_pole(k)
                     end if
                  end do
                  if (Sum(mass)+Sum(spotOut).lt.srts) exit massLoop
               end do massLoop

               if (massSteps>=maxMassSteps) then
                  !Failure of the algorithm, no valid solution for the masses could be found
                  write(*,*) 'Warning : Mass Iteration in massAss failed'
                  write(*,*) ' Minimal Masses', minMass
                  write(*,*) ' Maximal Masses', maxMass
                  success=.false.
                  write(*,*) 'ymin,ymax :' , ymin, ymax
                  stop 'assMass, Mass iteration'
                  return
               end if


               ! Choose absolute momenta by random, maximal value of momenta in CM-Frame is given by srts (assuming masses are zero)
               pcmAbs(1)=max(rn()*srts,1E-6)
               pcmAbs(3)=max(rn()*srts,1E-6)

               Ecm(1)=sqrt((mass(1)+spotOut(1))**2+pcmAbs(1)**2)
               Ecm(3)=sqrt((mass(3)+spotOut(3))**2+pcmAbs(3)**2)
               ! Due to energy conservation in CM-Frame :
               Ecm(2)=srts-Ecm(1)-Ecm(3)

               if (Ecm(2)<mass(2)+spotOut(2)) cycle momLoop        ! new try if mass > energy

               pcmAbs(2)=sqrt(Ecm(2)**2-(mass(2)+spotOut(2))**2)

               ! Check wether momentum conservation is possible :
               ! Evaluate Cosinus(theta) with theta between particles 1&3.  Using (pCM(1)+pCM(2)+pCM(3))**2=0 in vector notation.
               cosTheta_13 = (pcmAbs(2)**2-pcmAbs(1)**2-pcmAbs(3)**2) / ( 2.* pcmAbs(1) * pcmAbs(3) )
               if (abs(cosTheta_13)<=1.) exit momLoop
            end do momLoop

            if (momSteps >= maxMomSteps) then
               write(*,*) 'Momentum iteration in assMass failed. Momsteps=', momsteps
               write(*,*) IDOUT, SPotOut, chargeOUT,srts
               success=.false.
               stop
            end if

            !******************************************************************
            ! Set momenta
            ! (1) First choose arbitrarily z-axis in direction of momentum of third particle, later we'll randomly rotate coord-system
            momentum(:,3)=(/0.,0.,pcmAbs(3)/)

            phi=rn()*2.*pi
            sinTheta_13=SQRT(1.-cosTheta_13**2)
            momentum(:,1)=(/sinTheta_13*cos(phi),sinTheta_13*sin(phi),CosTheta_13/)*pCmAbs(1)

            ! Use momentum conservation in CM-Frame
            momentum(:,2)=-momentum(:,1)-momentum(:,3)

            ! Check consistency
            if (abs(pcmAbs(2)**2-Dot_Product(momentum(:,2),momentum(:,2)))>1E-5) then
               write(*,*) 'Consistency error in assMass calculations'
               write(*,*) momentum(:,2),Dot_Product(momentum(:,2),momentum(:,2)),".ne.", pcmAbs(2)
               stop 'assmass, consistency check'
            end if

            ! (2) Now we rotate the coordinate system by some random angle
            phi=rn()*2.*pi
            cosTheta=2.*(rn()-0.5)
            do k=1,3
               momentum(1:3,k) = rotateYZ (0.0, phi, momentum(1:3,k), cosTheta=cosTheta)
            end do
            ! Finished setting momenta
            !******************************************************************

            ! Calculate spectral functions
            do k=1,3
               if (.not.flagStable(k)) then
                  ! Determine momenta in LRF for evaluation of the width
                  pLRF(1:3)=momentum(:,k)
                  plrf(0)=sqrt((mass(k)+spotOut(k))**2+plrf(1)**2+plrf(2)**2+plrf(3)**2)
                  betaCMtoLab=-betaToCM
                  call lorentz(betaCMtoLab(1:3), plrf(0:3), 'finalState(3)')         ! boost from CM to Lab
                  call lorentz(betaToLRF(1:3), plrf(0:3), 'finalState(4)')             ! boost from Lab to LRF
                  if (isMeson(idOut(k))) then
                     gamtot=WidthMesonMedium(IDOut(k),mass(k), plrf(0:3) ,medium_ATCollision)
                  else
                     gamtot=WidthBaryonMedium(IDOut(k),mass(k),plrf(0:3) ,medium_ATCollision)
                  end if
                  spectral(k)=mass(k)**2*gamtot*gamma_pole(k)/ ((mass_pole(k)**2-mass(k)**2)**2+gamtot**2*  mass(k)**2)
                  intfac(k)=gamma_pole(k)**2/ ((mass(k)-mass_pole(k))**2+gamma_pole(k)**2/4.)/4.
               else
                  spectral(k)=1.
                  intfac(k)=1.
               end if
            end do

            probability = spectral(1)*spectral(2)*spectral(3) / (intfac(1)*intfac(2)*intfac(3)) &
                          * pcmAbs(1)*pcmAbs(3)/(ecm(1)*ecm(3)*maxValue)

            if (probability>=1) call checkMaximum(spectral(1)*spectral(2)*spectral(3)/ (intfac(1)*intfac(2)*intfac(3)) &
                                                  * pcmAbs(1)*pcmAbs(3)/(ecm(1)*ecm(3)), maxValue)

            if (probability>rn()) exit outerLoop ! success!!!

         end do OuterLoop

         if (outerSteps < maxOuterSteps) then
            success=.true.
         else
            success=.false.
            write(*,*) 'Problem in assMass. Not enough iterations to find good final state'
            write(*,*) ' outerSteps=', outerSteps
            write(*,*) ' MaxOuterSteps=', MaxouterSteps
            return
            !                stop ' assMass, outerloop'
         end if

      end if

    end subroutine throwDice

  end subroutine assMass_Full


  !****************************************************************************
  !****s* finalState_Full/massass_nBody
  ! NAME
  ! subroutine massass_nBody(srts,betaToLRF,betaToCM,mediumAtCollision,finalState,success)
  ! PURPOSE
  ! Determines the masses and momenta of outgoing particles in the c.m. frame
  ! of incoming particles.
  ! INPUTS
  ! * real :: srts                                ! Invariant energy
  ! * real, dimension(1:3) :: betaToLRF           ! Boost from computational frame to LRF
  ! * real, dimension(1:3) :: betaToCM            ! Boost from computational frame to CM frame
  ! * type(medium) :: mediumAtCollision           ! Medium info at collision point
  ! * type(particle), dimension(:) :: finalState  ! Id's and charges of the outgoing particles
  ! OUTPUT
  ! * logical                      :: success     ! .true. if masses and momenta are set
  ! * type(particle), dimension(:) :: finalState  ! Masses and four-momenta of the final state particles
  ! NOTES
  ! n=2,...,6 outgoing particles can be treated.
  ! NOTES
  ! This subroutine has a similar structure as subroutine massAss.
  !****************************************************************************
  subroutine massass_nBody(srts,betaToLRF,betaToCM,mediumAtColl,finalState,success)

  use particleDefinition, only: particle
  use particleProperties, only: hadron
  use mediumDefinition
  use random, only: rn
  use nBodyPhaseSpace, only: momenta_in_nBodyPS, integrate_nBodyPS
  use lorentzTrafo, only: lorentz
  use MesonWidthMedium, only: WidthMesonMedium
  use BaryonWidthMedium, only: WidthBaryonMedium
  use IdTable, only: rho, isMeson
  use constants, only: mPi

  real, intent(in) :: srts
  real, dimension(1:3), intent(in) :: betaToLRF
  real, dimension(1:3), intent(in) :: betaToCM
  type(medium), intent(in) :: mediumAtColl
  type(particle), dimension(:), intent(inout) :: finalState
  logical, intent(out)                        :: success

  logical, allocatable, dimension(:) :: flagStable
  real, allocatable, dimension(:) :: gamma_pole, mass_pole, minMass, maxMass, normMass,&
                                      &ymax, ymin, mass
  real, allocatable, dimension(:,:) :: pn

  real :: sumMinMass, phaseSpace_norm, phaseSpace, y, probability, spectral, gamtot, intfac
  real, dimension(0:3) :: plrf

  logical, parameter :: flagMomDep=.false.  ! if .true. --- take into account momentum dependence
                                            ! in spectral functions (takes more CPU time)

  real :: maxValue                        ! Adjustable parameter to have acceptance probability less than 1
  integer, parameter :: maxMassSteps=1000
  integer, parameter :: maxAttempts=20000
  integer :: i, massSteps
  integer :: k,n
  logical :: flagAllStable
  logical, save :: flagPlot=.true.


  if (allocated(flagStable)) deallocate(flagStable)
  if (allocated(gamma_pole)) deallocate(gamma_pole)
  if (allocated(mass_pole)) deallocate(mass_pole)
  if (allocated(minMass)) deallocate(minMass)
  if (allocated(maxMass)) deallocate(maxMass)
  if (allocated(normMass)) deallocate(normMass)
  if (allocated(ymin)) deallocate(ymin)
  if (allocated(ymax)) deallocate(ymax)
  if (allocated(mass)) deallocate(mass)
  if (allocated(pn)) deallocate(pn)

  n=size(finalState,dim=1)
  allocate(flagStable(1:n),gamma_pole(1:n),mass_pole(1:n),minMass(1:n),maxMass(1:n))
  allocate(normMass(1:n))
  allocate(ymin(1:n),ymax(1:n),mass(1:n))
  allocate(pn(1:3,1:n))

  ! Check whether there are unstable particles in finalState:
  flagAllStable=.true.
  do k=1,n

     gamma_Pole(k) = hadron(finalState(k)%Id)%width
     mass_Pole(k)  = hadron(finalState(k)%Id)%mass
     minMass(k)    = hadron(finalState(k)%Id)%minmass

     if (gamma_pole(k).lt.1.e-03) then
        flagStable(k)=.true.
        minMass(k)=mass_Pole(k)
     else
        flagStable(k)=.false.
        flagAllStable=.false.
     end if

  end do

  sumMinMass=sum(minMass(1:n))

  ! Check kinematics:
  if (sumMinMass.gt.srts) then
     success=.false.
     return
  end if

  if (flagAllStable) then  ! Only stable outgoing particles
     finalState(1:n)%mass=mass_Pole(1:n)
     ! Monte-Carlo sampling of the phase space:
     call momenta_in_nBodyPS(srts,finalState(1:n)%mass,pn(1:3,1:n))
     do k=1,n
        finalState(k)%momentum(1:3)=pn(1:3,k)
        finalState(k)%momentum(0)=sqrt(finalState(k)%mass**2+dot_product(pn(1:3,k),pn(1:3,k)))
     end do
     success=.true.
     return
  end if

  ! There are unstable outgoing particles

  ! Determine the normalization factor for the acceptance probability:
  maxValue=4.

  ! Do variable transformation mass-> y
  do k=1,n
     if (.not.flagStable(k)) then
        maxMass(k)= srts - sumMinMass + minMass(k)
        if (maxMass(k).lt.minMass(k)) then
           write(*,*) 'problems in setKinematics_annihilation maxmass.lt.minmass',&
                    & srts,maxmass(k),minmass(k)
           write(*,*) finalState(1:n)%Id
           stop
        end if
        ymax(k)= 2.*atan((maxmass(k)-mass_pole(k))/gamma_pole(k)*2.)
        ymin(k)= 2.*atan((minmass(k)-mass_pole(k))/gamma_pole(k)*2.)
     else
        maxMass(k)=mass_pole(k)
     end if
  end do

!  Testing:
  if (flagPlot) then
     call plot
     flagPlot=.false.
  end if

  ! Determine phase space volume for normalization of acceptance probability:
  do k=1,n
    if (finalState(k)%Id.eq.rho) then
       normMass(k)=mass_Pole(k)-gamma_Pole(k)
    else
       normMass(k)=max(mPi,mass_Pole(k)-3.*gamma_Pole(k))
    end if
  end do

  if (sum(normMass(1:n)).gt.srts) normMass(1:n)=0.5*(normMass(1:n)+minMass(1:n))
  if (sum(normMass(1:n)).gt.srts) normMass(1:n)=minMass(1:n)

  phaseSpace_norm = integrate_nBodyPS (srts, normMass(1:n))

  ! Start monte carlo to find masses and (optionally) momenta:
  attemptLoop : do i=1,maxAttempts

      ! Monte-Carlo choice of the outgoing masses:
      massLoop : do massSteps=1,maxMassSteps

         do k=1,n
            if (.not.flagStable(k)) then
                y=ymin(k)+rn()*(ymax(k)-ymin(k))
                mass(k)=.5*tan(y/2.)*gamma_pole(k)+ mass_pole(k)
            else
                mass(k)=mass_pole(k)
            end if
         end do

         if (Sum(mass(1:n)).lt.srts) exit massLoop

         if (massSteps.eq.maxMassSteps) then
            write(*,*) 'problems in massass_nBody: mass choice is impossible',&
            & srts, sum(mass_Pole(1:n)), finalState(1:n)%Id
            success=.false.
            return
         end if

      end do massLoop

      if (flagMomDep) then
         ! Monte-Carlo sampling of the phase space:
         call momenta_in_nBodyPS(srts,mass(1:n),pn(1:3,1:n))
      else
         pn=0.
      end if

      ! Determine the phase space factor:
      phaseSpace = integrate_nBodyPS (srts, mass(1:n))

      ! Calculate spectral functions:
      probability=1.
      do k=1,n
         if (.not.flagStable(k)) then
            ! Determine momenta in LRF for evaluation of the width
            plrf(1:3)=pn(1:3,k)
            plrf(0)=sqrt(mass(k)**2+plrf(1)**2+plrf(2)**2+plrf(3)**2)
            call lorentz(-betaToCM, plrf, 'massass_nBody(1)')       ! boost from CM to comput. frame
            call lorentz(betaToLRF, plrf, 'massass_nBody(2)')       ! boost from comput. frame to LRF
            if (isMeson(finalState(k)%Id)) then
               gamtot=WidthMesonMedium(finalState(k)%Id, mass(k), plrf(0:3), mediumAtColl)
            else
               gamtot=WidthBaryonMedium(finalState(k)%Id, mass(k), plrf(0:3), mediumAtColl)
            end if
            spectral=mass(k)**2*gamtot*gamma_pole(k) / ( (mass_pole(k)**2-mass(k)**2)**2 &
                          &                             +gamtot**2*mass(k)**2 )
            intfac=gamma_pole(k)**2/ ((mass(k)-mass_pole(k))**2+gamma_pole(k)**2/4.)/4.
         else
            spectral=1.
            intfac=1.
         end if
         probability=probability*spectral/intfac
      end do

      probability=probability*phaseSpace/phaseSpace_norm/maxValue

      ! write(*,*)' probability: ', probability

      if (probability.gt.1.) then
         write(*,*) ' problem in massass_nBody, probability.gt.1: ', probability
         write(*,*) srts, sum(mass_Pole(1:n)), finalState(1:n)%Id
         maxValue=maxValue*2.
         write(*,*) ' maxValue is readjusted to ', maxValue
         write(*,*) ' If this problem repeats too often, choose another initial maxValue'
         cycle attemptLoop
      end if

      if (rn().lt.probability) then
         success=.true.
         exit attemptLoop
      end if

      if (i.eq.maxAttempts) then
         write(*,*) 'problems in massass_nBody: too many attempts',&
                 & srts, sum(mass_Pole(1:n)), finalState(1:n)%Id
         success=.false.
         return
      end if

  end do attemptLoop

  if (.not.flagMomDep) call momenta_in_nBodyPS(srts,mass(1:n),pn(1:3,1:n))

  ! Set the outgoing particle vector:

  do k=1,n
     finalState(k)%mass=mass(k)
     finalState(k)%momentum(1:3)=pn(1:3,k)
     finalState(k)%momentum(0)=sqrt(mass(k)**2+dot_product(pn(1:3,k),pn(1:3,k)))
  end do


  contains

     subroutine plot

     use constants, only: mPi

     integer :: i
     real :: gamma_pole, mass_pole, mass

     open(1,file='massass_nBody.dat',status='unknown')
     ! Determine maximum possible phase space volume:
     phaseSpace_norm = integrate_nBodyPS (3.024, (/mPi,mPi,mPi,0.619,0.619,0.619/))
     do i = 1,2000
        mass=2*mPi+0.001*float(i)
        plrf(1:3)=0.
        plrf(0)=mass
        call lorentz(-betaToCM, plrf, 'massass_nBody(1)')       ! boost from CM to comput. frame
        call lorentz(betaToLRF, plrf, 'massass_nBody(2)')       ! boost from comput. frame to LRF
        gamtot=WidthMesonMedium(rho,mass,plrf(0:3),mediumAtColl)
        gamma_pole=hadron(rho)%width
        mass_pole=hadron(rho)%mass
        spectral=mass**2*gamtot*gamma_pole / ( (mass_pole**2-mass**2)**2 &
                          &                             +gamtot**2*mass**2 )
        intfac=gamma_pole**2/ ((mass-mass_pole)**2+gamma_pole**2/4.)/4.
        phaseSpace = integrate_nBodyPS (3.024, (/mPi,mPi,mPi,mass,mass,mass/))
        write(1,'(5(1x,e13.6))') mass,spectral,intfac,phaseSpace,phaseSpace_norm
     end do

     end subroutine plot


  end subroutine massass_nBody





end module finalState_Full
