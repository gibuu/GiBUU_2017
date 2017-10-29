!******************************************************************************
!****m* /ElementaryAnalysis
! NAME
! module ElementaryAnalysis
!
! PURPOSE
! Contains an output program for elementary particle collision.
!******************************************************************************
module ElementaryAnalysis

  implicit none
  private

  public :: DoElementaryAnalysis


  !****************************************************************************
  !****g* ElementaryAnalysis/DoOutChannels
  ! SOURCE
  !
  logical, save :: DoOutChannels = .false.
  ! PURPOSE
  ! switch on/off: reporting of all final state channels
  !****************************************************************************

  !****************************************************************************
  !****g* ElementaryAnalysis/DoH2d
  ! SOURCE
  !
  logical, save :: DoH2d=.false.
  ! PURPOSE
  ! if .true. than make output of 2-dimensional histograms
  ! (they could produce files of size 240 mb)
  !****************************************************************************

  !****************************************************************************
  !****g* ElementaryAnalysis/Do45ForAllEvents
  ! SOURCE
  !
  logical, save :: Do45ForAllEvents = .false.
  ! PURPOSE
  ! flag to decide,
  ! whether DoElementaryAnalysis4(5).dat is written for all events
  ! or just for events, where the output channel consist of pions
  !****************************************************************************

  !****************************************************************************
  !****g* ElementaryAnalysis/DodNNbar
  ! SOURCE
  !
  logical, save :: DodNNbar=.false.
  ! PURPOSE
  !****************************************************************************

  !****************************************************************************
  !****g* ElementaryAnalysis/DoPanda
  ! SOURCE
  !
  logical, save :: DoPanda=.false.
  ! PURPOSE
  ! if .true., elementary analysis for channels with S = -2 and -3 (Xi, Omega)
  !****************************************************************************

  !****************************************************************************
  !****g* ElementaryAnalysis/Dodsigdt
  ! SOURCE
  !
  logical, save :: Dodsigdt=.false.
  ! PURPOSE
  !****************************************************************************

  !****************************************************************************
  !****g* ElementaryAnalysis/Do2Part
  ! SOURCE
  !
  logical, save :: Do2Part=.false.
  ! PURPOSE
  !****************************************************************************


  logical, save :: FlagReadInput = .true.

contains

  !****************************************************************************
  !****s* ElementaryAnalysis/readInput
  ! NAME
  ! subroutine readInput
  ! PURPOSE
  ! Reads input in jobcard out of namelist "Elementary_Analysis".
  !****************************************************************************
  subroutine readInput
    use output, only: Write_ReadingInput

    integer :: ios

    !**************************************************************************
    !****n* ElementaryAnalysis/Elementary_Analysis
    ! NAME
    ! NAMELIST /Elementary_Analysis/
    ! PURPOSE
    ! Includes the switches:
    ! * DoOutChannels
    ! * DoH2d
    ! * Do45ForAllEvents
    ! * DodNNbar
    ! * DoPanda
    ! * Dodsigdt
    ! * Do2Part
    !**************************************************************************
    NAMELIST /Elementary_Analysis/ DoOutChannels, DoH2d, Do45ForAllEvents, &
                                   Dodsigdt, Do2Part, DoPanda

    call Write_ReadingInput('Elementary_Analysis',0)
    rewind(5)
    read(5,nml=Elementary_Analysis,IOSTAT=IOS)
    call Write_ReadingInput('Elementary_Analysis',0,IOS)

    write(*,*) 'DoOutChannels    : ',DoOutChannels
    write(*,*) 'DoH2d            : ',DoH2d
    write(*,*) 'Do45ForAllEvents : ',Do45ForAllEvents
    write(*,*) 'DodNNbar         : ',DodNNbar
    write(*,*) 'Dodsigdt         : ',Dodsigdt
    write(*,*) 'Do2Part          : ',Do2Part
    write(*,*) 'DoPanda          : ',DoPanda

    call Write_ReadingInput('Elementary_Analysis',1)
    FlagReadInput = .false.
  end subroutine readInput

  !****************************************************************************
  !****s* ElementaryAnalysis/DoElementaryAnalysis
  ! NAME
  ! subroutine DoElementaryAnalysis (realparticles, finalFlag)
  !
  ! PURPOSE
  ! Computes various cross sections for a collision of two elementary
  ! particles.
  !
  ! INPUTS
  ! * type(particle), dimension(:,:) :: realparticles  --  real particle vector
  ! * logical :: finalFlag -- .true. if it is the last call for one specific
  !   energy, therefore final output must be made.
  !
  ! RESULT
  ! * total, elastic and strangeness production cross sections (mbarn)
  ! * pion (pi^-,pi^0,pi^+),
  !   kaon (K^0,K^+), antikaon (K^-,Kbar^0),
  !   pion Lambda and pion Sigma (S^-,S^0,S^+)
  !   production cross sections (mbarn)
  !
  ! All those results are printed to the files:
  ! * 'DoElementaryAnalysis1.dat'
  ! * 'DoElementaryAnalysis2.dat'
  ! * 'DoElementaryAnalysis3.dat'
  ! * 'DoElementaryAnalysis4.dat'
  ! * 'DoElementaryAnalysis5.dat'
  ! * ....etc
  ! See the headers of these files for their contents.
  !
  ! NOTES
  ! * The files "OutChannels.nnn.dat" list the cross section for every
  !   ouput channel that occured in the run !
  !****************************************************************************
  subroutine DoElementaryAnalysis (realparticles, finalFlag)

    use IdTable
    use particleDefinition
    use twoBodyTools, only: MomentumTransfer2, IsElastic, IsChargeExchange
    use initElementary, only: impactParameter, srts, p_lab, particleID, target, projectile, siggeo, particleCharge, particleAnti, &
                              ekin_lab, particlemass
    use ParticleProperties, only: PartName,isStrange
    use XsectionRatios, only: getShift0
    use output, only: intToChar
    use AnaEventDefinition, only: tAnaEvent
    use AnaEvent, only: event_add, event_CLEAR
    use preEventDefinition
    use lorentzTrafo, only: lorentz
    use histf90
    use hist2Df90
    use particlePointerListDefinition
    use particlePointerList, only: ParticleList_INIT
    use PreEvListDefinition
    use PreEvList, only: CreateSortedPreEvent, PreEvList_CLEAR, PreEvList_INIT, PreEvList_INSERT, PreEvList_Print
    use constants, only: pi, mN, mPi

    type(particle), dimension(:,:), intent(in), target :: realparticles
    logical, intent (in)             :: finalFlag

    character(15),dimension(2) :: name
    real, save :: sigma_total,sigma_elastic,sigma_strangeness
    real, save :: sigma_CEX,sigma_PROD_pbarp,sigma_ANN,sigma_Y
    real, save :: sigma_Lambda_LambdaBar,sigma_Lambda_LambdaBar_X,sigma_Lambda_LambdaBar_Pions(1:10), &
                  sigma_Lambda_LambdaBar_KKbar_X,sigma_Lambda_LambdaBar_NNbar_X,sigma_Lambda_SigmaBar_cc
    real, save :: sigma_pion(-1:1),sigma_kaon(0:1,1:2),sigma_antikaon(-1:0,1:2),sigma_KKbar(1:19),sigma_pipi(1:4)
    real, save :: sigma_NucNuc,sigma_NucBarNuc,sigma_1pion_pbarp(0:1,-1:1),sigma_2pion_pbarp, &
                  sigma_3pion_pbarp(-1:1),sigma_4pion_pbarp,sigma_2pion,sigma_2pion_Eta,sigma_2pion_KKbar
    real, save :: sigma_antikaon_nucleon
    real, save :: sigma_pion_Lambda,sigma_pion_Sigma(-1:1)
    integer, parameter :: N_max=100             ! Maximum number of particles produced in an event
    real, save :: P_Npion(0:N_max)             ! Pion multiplicity distribution
    real, save :: sigma_prong(0:N_max)         ! topological Xsections
    integer, parameter :: Nmom=100             ! Number of momentum bins
    real, parameter :: dmom=0.02               ! Momentum bin (GeV/c)
    real, save, dimension(1:Nmom,0:10) :: dNpiondMom  ! Pion momentum distribution
    real, save :: fnorm, pion_events

    real, save :: sigma_LpToS0p,sigma_SmpToLn,sigma_SmpToS0n !Lambda Nuk <--> Sigma Nuk
    real, save :: sigma_I0LL,sigma_I1SL !Xi Nuk --> Lambda Lambda, Lambda Sigma
    real, save :: sigma_Xi_XiBar,sigma_Xi_XiBar_X !pBarp --> XiBarXi(+X)
    real, save :: sigma_Xi0_Xi0Bar, sigma_Xi0_Xi0Bar_X, sigma_JPsi
    real, save :: sigma_LS0, sigma_Xi0n, sigma_LSm

    real, save :: sigma_OmegaBaryon, sigma_OmegaToXi
    real, save :: sigma_OmegaP_to_Xi0Lambda,sigma_OmegaP_to_Xi0Sigma0,sigma_OmegaP_to_XimSigmap
    real, save :: sigma_OmegaN_to_XimLambda,sigma_OmegaN_to_XimSigma0,sigma_OmegaN_to_Xi0Sigmam
    real, save :: sigma_OmegaBaryonX, sigma_OmegaToXiX
    real, save :: sigma_OmegaP_to_Xi0LambdaX,sigma_OmegaP_to_Xi0Sigma0X,sigma_OmegaP_to_XimSigmapX
    real, save :: sigma_OmegaN_to_XimLambdaX,sigma_OmegaN_to_XimSigma0X,sigma_OmegaN_to_Xi0SigmamX

    real, save :: sigma_KbarNuk_to_XiK(1:5),sigma_KbarNuk_to_XiMinus_X !KBar Nuk --> Xi K

    integer :: i,j,k,l,numpart,totCharge,numPions,numEtas,numBaryons,numChPart,ibin
    integer :: numPi(-1:1), numN(0:1), numNbar(-1:0), numLambda, numLambdaBar
    integer :: numSigma(-1:1), numSigmaBar(-1:1), numXi(-1:0), numXiBar(0:1)
    integer :: numKaon(0:1), numKaonBar(-1:0)

    integer :: numOmega, numOmegaBar !strange Omega (anti)baryon (S=-3/+3)!

    integer :: numensembles
    integer, save :: ncall=0   ! counter of calls
    integer, save :: isu=0     ! counter of subsequent runs at fixed energy
    integer, save :: nwrite=1  ! counter of calls with writing
    logical :: flag3
    logical :: flag_not_only_pions, flag_inel
    real :: ptot(0:3), momentum(0:3), beta(1:3), s, momentumAbs, pinumAv
    type(tAnaEvent) :: event
    type(tPreEvListEntry),save :: preEv
    type(tPreEvList), save :: ListPreEv
    type(particle) :: first, second    ! First and second particles (with nonzero Id's) in the order as stored
                                       ! in realParticles array
    type(particle), POINTER :: pPart

    type(histogram2D), save :: h2dMomentPion(-1:1),h2dEnergyPion(-1:1)
    type(histogram), save   :: hEnergyPion(-1:1)
    type(histogram), save   :: dsigdt_elastic,dsigdt_CEX,dsigdt_Delta,dsigdt_DeltaBar,dsigdt_NucleonBar
    type(histogram), save   :: dsigdt_LambdaBarLambda,dsigdt_SigmaBarLambda
    type(histogram), save   :: dsigdm_Delta,dsigdm_DeltaBar
    type(histogram), save   :: dNpiondcosTheta    ! Pion polar angle distribution
    type(histogram), save   :: dNpiondPhi    ! Pion azimuthal angle distribution
    type(histogram), save   :: dNrhodM    ! Rho-meson mass distribution
    type(histogram), save   :: dNomegadM    ! Omega-meson mass distribution
    type(histogram), save   :: dNNbardEkin(1:5)  ! Antinucleon kin. energy distribution

    character(3), dimension(-1:1), parameter :: piName  = (/'pi-','pi0','pi+'/)
    real :: pT2,pL,pT2max,Emax,Q2,cosTheta,Phii,Ekin

    if (FlagReadInput) call readInput

    ncall=ncall+1
    isu=isu+1

    if (impactParameter >= 0.) then
       write(*,*) ' DoElementaryAnalysis can not be used for fixed impact parameter'
       write(*,*) ' You should set  impactParameter < 0 for integration mode !!!'
       stop
    end if

    if (DoOutChannels .and. ncall==1) then !===== Global Init =====
       call ParticleList_INIT(event%particleList)
       call PreEvList_INIT(ListPreEv)
    end if

    if (isu==1) then   !===== Init per Energy =====
       sigma_total=0.
       sigma_elastic=0.
       sigma_CEX=0.
       sigma_PROD_pbarp=0.
       sigma_ANN=0.
       sigma_Y=0.
       sigma_Lambda_LambdaBar=0.
       sigma_Lambda_LambdaBar_X=0.
       sigma_Lambda_LambdaBar_Pions=0.
       sigma_Lambda_LambdaBar_KKbar_X=0.
       sigma_Lambda_LambdaBar_NNbar_X=0.
       sigma_Lambda_SigmaBar_cc=0.
       sigma_JPsi=0.
       sigma_pion=0.
       sigma_NucNuc=0.
       sigma_NucBarNuc=0.
       sigma_1pion_pbarp=0.
       sigma_2pion_pbarp=0.
       sigma_3pion_pbarp=0.
       sigma_4pion_pbarp=0.
       sigma_kaon=0.
       sigma_antikaon=0.
       sigma_antikaon_nucleon=0.
       sigma_strangeness=0.
       sigma_pion_Lambda=0.
       sigma_pion_Sigma=0.
       P_Npion=0.
       sigma_prong=0.
       pion_events=0.
       dNpiondMom=0.
       sigma_2pion=0.
       sigma_2pion_Eta=0.
       sigma_2pion_KKbar=0.
       sigma_KKbar=0.
       sigma_pipi=0.

       PANDA0 : if (DoPanda) then

          !Kbar+N -> Xi + K
          sigma_KbarNuk_to_XiK=0.
          sigma_KbarNuk_to_XiMinus_X=0.

          !p+pBar -> Xi + XiBar (+X)
          sigma_Xi_XiBar_X=0.
          sigma_Xi_XiBar=0.
          sigma_Xi0_Xi0Bar=0.
          sigma_Xi0_Xi0Bar_X=0.

          !Hyperon+Nucleon scattering:
          sigma_LpToS0p = 0.
          sigma_SmpToLn = 0.
          sigma_SmpToS0n = 0.
          sigma_I0LL = 0.0
          sigma_I1SL = 0.0
          sigma_LS0 = 0.
          sigma_Xi0n = 0.
          sigma_LSm = 0.

          !channels relevant for Omega baryon (S=-3) production:
          sigma_OmegaBaryon = 0.0
          sigma_OmegaBaryonX = 0.0

          !\Omega N --> \Xi hyperon(S=-1)
          sigma_OmegaToXi = 0.0

          sigma_OmegaP_to_Xi0Lambda = 0.0
          sigma_OmegaP_to_Xi0Sigma0 = 0.0
          sigma_OmegaP_to_XimSigmap = 0.0
          sigma_OmegaN_to_XimLambda = 0.0
          sigma_OmegaN_to_XimSigma0 = 0.0
          sigma_OmegaN_to_Xi0Sigmam = 0.0

          sigma_OmegaToXiX = 0.0

          sigma_OmegaP_to_Xi0LambdaX = 0.0
          sigma_OmegaP_to_Xi0Sigma0X = 0.0
          sigma_OmegaP_to_XimSigmapX = 0.0
          sigma_OmegaN_to_XimLambdaX = 0.0
          sigma_OmegaN_to_XimSigma0X = 0.0
          sigma_OmegaN_to_Xi0SigmamX = 0.0

       end if PANDA0

       if (DoOutChannels) call PreEvList_CLEAR(ListPreEv)

       pT2max = srts/2
       do i=1,2
          if (isBaryon(particleId(i))) then
             pT2max = pT2max - mN
          else
             pT2max = pT2max - mPi
          end if
       end do
       pT2max = (pT2max)**2
       Emax = max(srts,p_lab+1.0)

       do i=-1,1
         call CreateHist(hEnergyPion(i), 'dsigma/dE, <pT2>(E) '//piName(i), &
             & 0.0, Emax, 0.02)

         if (DoH2d) then
            call CreateHist2D(h2dMomentPion(i), 'dsigma/dpL dpT2, '//piName(i), &
                & (/-1.0,0.0/), (/p_lab+1.0,pT2max/), (/0.02,0.02/) , .true.)
            call CreateHist2D(h2dEnergyPion(i), 'dsigma/dE dpT2, '//piName(i), &
                & (/0.0,0.0/), (/Emax,pT2max/), (/0.02,0.02/) , .true.)
         end if

       end do

       call CreateHist(dsigdt_elastic,'dsigma/dt (mb/GeV**2) vs. -t (GeV**2) elastic',0.,2.,0.01)
       call CreateHist(dsigdt_CEX,'dsigma/dt (mb/GeV**2) vs. -t (GeV**2) CEX',0.,1.,0.01)
       call CreateHist(dsigdt_NucleonBar,'dsigma/dt (mb/GeV**2) vs. -t (GeV**2) for NBar',0.,1.,0.001)
       call CreateHist(dsigdt_Delta,'dsigma/dt (mb/GeV**2) vs. -t (GeV**2) for Delta',0.,1.,0.01)
       call CreateHist(dsigdt_DeltaBar,'dsigma/dt (mb/GeV**2) vs. -t (GeV**2) for DeltaBar',0.,1.,0.01)
       call CreateHist(dsigdt_LambdaBarLambda,'dsigma/dt (mb/GeV**2) vs. -t (GeV**2) for LambdaBar',0.,1.5,0.01)
       call CreateHist(dsigdt_SigmaBarLambda,'dsigma/dt (mb/GeV**2) vs. -t (GeV**2) for SigmaBar',0.,1.5,0.01)
       call CreateHist(dsigdm_Delta,'dsigma/dm (mb/GeV) vs. m_delta (GeV) for Delta',1.,2.,0.01)
       call CreateHist(dsigdm_DeltaBar,'dsigma/dm (mb/GeV) vs. m_delta (GeV) for DeltaBar',1.,2.,0.01)
       call CreateHist(dNpiondcosTheta,'dNpion/dcosTheta vs. Theta for annihilation',-1.,1.,0.02)
       call CreateHist(dNpiondPhi,'dNpion/dPhi vs. Phi for annihilation',-pi,pi,2.*pi/100.)
       call CreateHist(dNrhodM,'dNrho/dM vs. M',0.2,1.5,0.01)
       call CreateHist(dNomegadM,'dNomega/dM vs. M',0.2,1.5,0.01)

       if (DodNNbar) then
          do i=1,5
             call CreateHist(dNNbardEkin(i),'dNNbar/dEkin vs. Ekin',0.,2.,0.01)
          end do
       end if

    end if

    ! C.M. velocity of the colliding particles:
    beta(1:3)=(/0.,0.,p_lab/)/sqrt(srts**2+p_lab**2)

    numensembles=size(realParticles,dim=1)

    Ensemble_loop : do i = 1,numensembles

       numpart=0
       totCharge=0
       numPions=0
       numEtas=0
       numBaryons=0
       numChPart=0
       numPi(-1:1)=0
       numN(0:1)=0
       numNbar(-1:0)=0
       numLambda=0
       numLambdaBar=0
       numSigma(-1:1)=0
       numSigmaBar(-1:1)=0
       numXi(-1:0)=0
       numXiBar(0:1)=0
       numKaon(0:1)=0
       numKaonBar(-1:0)=0

       numOmega = 0
       numOmegaBar = 0

       flag3=.false.
       flag_not_only_pions=.false.
       ptot=0.

       if (DoOutChannels) call event_CLEAR(event)

       Particle_loop1 : do j = 1,size(realParticles,dim=2)

          if (realParticles(i,j)%ID <= 0) cycle Particle_loop1

          numpart=numpart+1
          if (numpart.eq.1) first=realParticles(i,j)
          if (numpart.eq.2) second=realParticles(i,j)

          if (numpart > N_max) then
            write(*,*) ' Too many outgoing particles in DoElementaryAnalysis !!!'
            stop
          end if

          if (realParticles(i,j)%event(1)<=2) cycle Ensemble_loop

          pPart => realParticles(i,j)

          if (DoOutChannels) call event_add(event,pPart)

          ptot=ptot+pPart%momentum
          Q2=MomentumTransfer2(pPart,projectile)

          totCharge=totCharge+pPart%charge

          if (pPart%charge .ne. 0) numChPart=numChPart+1

          select case (pPart%ID)
          case (nucleon)
             if (.not.pPart%antiParticle) then
               numN(pPart%charge)=numN(pPart%charge)+1
             else
               numNbar(pPart%charge)=numNbar(pPart%charge)+1
               call AddHist(dsigdt_NucleonBar,Q2,1.)
             end if
          case (delta)
             if (.not.pPart%antiParticle) then
               call AddHist(dsigdt_Delta,Q2,1.)
               call AddHist(dsigdm_Delta,pPart%mass,1.)
             else
               call AddHist(dsigdt_DeltaBar,Q2,1.)
               call AddHist(dsigdm_DeltaBar,pPart%mass,1.)
             end if
          case (Lambda)
             if (.not.pPart%antiParticle) then
               numLambda=numLambda+1
             else
               numLambdaBar=numLambdaBar+1
             end if
          case (SigmaResonance)
             if (.not.pPart%antiParticle) then
               numSigma(pPart%charge)=numSigma(pPart%charge)+1
             else
               numSigmaBar(pPart%charge)=numSigmaBar(pPart%charge)+1
             end if
          case (Xi)
             if (.not.pPart%antiParticle) then
               numXi(pPart%charge)=numXi(pPart%charge)+1
             else
               numXiBar(pPart%charge)=numXiBar(pPart%charge)+1
             end if
          case (omegaResonance)
             if (.not.pPart%antiParticle) then
               numOmega=numOmega+1
             else
               numOmegaBar=numOmegaBar+1
             end if
          case (pion)
             sigma_pion(pPart%charge)=sigma_pion(pPart%charge)+1.
             numPions=numPions+1
             numPi(pPart%charge)=numPi(pPart%charge)+1

             pL = pPart%momentum(3)
             pT2= pPart%momentum(1)**2+pPart%momentum(2)**2

             if (DoH2d) then
                call AddHist2D(h2dMomentPion(pPart%charge), (/pL,pT2/), 1.0)
                call AddHist2D(h2dEnergyPion(pPart%charge), (/pPart%momentum(0),pT2/), 1.0)
             end if
             call AddHist(hEnergyPion(pPart%charge), pPart%momentum(0), 1.0, pT2)
          case (eta)
             numEtas=numEtas+1
          case (kaon)
             numKaon(pPart%charge)=numKaon(pPart%charge)+1
          case (kaonBar)
             numKaonBar(pPart%charge)=numKaonBar(pPart%charge)+1
          case (JPsi)
             sigma_JPsi=sigma_JPsi+1.
          end select

          if (isStrange(pPart%ID)) flag3 = .true. ! do not reset to .false. !

          if (pPart%ID.ne.pion) flag_not_only_pions=.true.

          if (isBaryon(pPart%ID) .and. .not.pPart%antiParticle) numBaryons=numBaryons+1

       end do Particle_loop1

       !if( totCharge.ne.particleCharge(1)+particleCharge(2) .and. numBaryons.ne.0 ) then
       !  write(*,*)' In DoElementaryAnalysis: charge not conserved !'
       !  write(*,*)' Summary charge of incoming particles:', particleCharge(1)+particleCharge(2)
       !  write(*,*)' Total charge of the produced event:', totCharge
       !  Particle_loop_print : do j = 1,size(realParticles,dim=2)
       !     if(realParticles(i,j)%ID <= 0) cycle Particle_loop_print
       !     write(*,*) realParticles(i,j)%ID, realParticles(i,j)%charge
       !  end do Particle_loop_print
       !  stop
       !end if

       sigma_total=sigma_total+1.

       s=ptot(0)**2-dot_product(ptot(1:3),ptot(1:3))

       flag_inel=.true.

       if (numpart==2) then

         if (IsElastic((/projectile,target/),(/first,second/)) .and. abs(sqrt(s)-srts)<1.e-08 ) then
              sigma_elastic=sigma_elastic+1.
              flag_inel=.false.
              if ( IsSamePart(first,projectile) ) then
                  Q2=MomentumTransfer2(first,projectile)
              else
                  Q2=MomentumTransfer2(second,projectile)
              end if
              call AddHist(dsigdt_elastic,Q2,1.)
         end if

         if (IsChargeExchange(projectile,target,first,second)) then
              sigma_CEX=sigma_CEX+1.
              flag_inel=.false.
              if (first%Id.eq.projectile%Id .and. (first%antiParticle.eqv.projectile%antiParticle)) then
                  Q2=MomentumTransfer2(first,projectile)
              else
                  Q2=MomentumTransfer2(second,projectile)
              end if
              call AddHist(dsigdt_CEX,Q2,1.)
         end if

         if (first%Id.eq.Lambda .and. second%Id.eq.Lambda) then
             if (first%antiParticle.eqv.projectile%antiParticle) then
                 Q2=MomentumTransfer2(first,projectile)
             else
                 Q2=MomentumTransfer2(second,projectile)
             end if
             call AddHist(dsigdt_LambdaBarLambda,Q2,1.)
         end if


         if (first%Id.eq.SigmaResonance .and. second%Id.eq.Lambda .or. &
           &first%Id.eq.Lambda .and. second%Id.eq.SigmaResonance  ) then
             if (first%antiParticle.eqv.projectile%antiParticle) then
                 Q2=MomentumTransfer2(first,projectile)
             else
                 Q2=MomentumTransfer2(second,projectile)
             end if
             call AddHist(dsigdt_SigmaBarLambda,Q2,1.)
         end if



       end if

       if (numChPart <= N_max) sigma_prong(numChPart)=sigma_prong(numChPart)+1.

       if (flag3) sigma_strangeness=sigma_strangeness+1.

       if (numpart==2.and.numPions==1) then
          if (numLambda.eq.1) then
            sigma_pion_Lambda=sigma_pion_Lambda+1.
          else if (numSigma(-1).eq.1) then
            sigma_pion_Sigma(-1)=sigma_pion_Sigma(-1)+1.
          else if (numSigma(0).eq.1) then
            sigma_pion_Sigma(0)=sigma_pion_Sigma(0)+1.
          else if (numSigma(1).eq.1) then
            sigma_pion_Sigma(1)=sigma_pion_Sigma(1)+1.
          end if
       end if


       if ( numLambda > 0 .or. numLambdaBar > 0 .or. sum(numSigma(:)) > 0 &
         & .or. sum(numSigmaBar(:)) > 0 ) sigma_Y=sigma_Y+1.

       if (numpart==2 .and. sum(numN(:))==1 .and. sum(numKaonBar(:))==1) &
         &  sigma_antikaon_nucleon=sigma_antikaon_nucleon+1.

       if (numpart==2 .and. sum(numN(:))==2) sigma_NucNuc=sigma_NucNuc+1.

       if (numpart==2 .and. sum(numN(:))==1 .and. sum(numNbar(:))==1) sigma_NucBarNuc=sigma_NucBarNuc+1.

       if (numpart==3 .and. sum(numNbar(:))==1) then
         do k=0,1
           do l=-1,1
             if (numN(k)==1 .and. numPi(l)==1) sigma_1pion_pbarp(k,l)=sigma_1pion_pbarp(k,l)+1.
           end do
         end do
       else if (numpart==4 .and. numNbar(-1)==1 .and. numN(1)==1&
              & .and. numPi(-1)==1 .and. numPi(1)==1) then ! pbar p pi- pi+
         sigma_2pion_pbarp=sigma_2pion_pbarp+1.
       else if (numpart==5 .and. sum(numNbar(:))==1 .and. sum(numN(:))==1&
                            & .and. sum(numPi(:))==3) then
         if (numNbar(-1)==1 .and. numN(1)==1 .and. numPi(0)==1) then ! pbar p pi+ pi- pi0
           sigma_3pion_pbarp(0)=sigma_3pion_pbarp(0)+1.
         else if (numNbar(-1)==1 .and. numN(0)==1 .and. numPi(0)==0) then ! pbar n pi+ pi+ pi-
           sigma_3pion_pbarp(1)=sigma_3pion_pbarp(1)+1.
         else if (numNbar(0)==1 .and. numN(1)==1 .and. numPi(0)==0) then ! nbar p pi- pi- pi+
           sigma_3pion_pbarp(-1)=sigma_3pion_pbarp(-1)+1.
         end if
       else if (numpart==6 .and.  numNbar(-1)==1 .and. numN(1)==1&
              & .and. numPi(-1)==2 .and. numPi(1)==2) then ! pbar p 2pi- 2pi+
         sigma_4pion_pbarp=sigma_4pion_pbarp+1.
       end if

       if (numpart>2 .and. sum(numNbar(:))==1 .and. sum(numN(:))==1) sigma_PROD_pbarp=sigma_PROD_pbarp+1.

       if (numBaryons==0) then
          sigma_ANN=sigma_ANN+1.
          sigma_kaon(:,1)=sigma_kaon(:,1)+numKaon(:)
          sigma_antikaon(:,1)=sigma_antikaon(:,1)+numKaonBar(:)
       else
          sigma_kaon(:,2)=sigma_kaon(:,2)+numKaon(:)
          sigma_antikaon(:,2)=sigma_antikaon(:,2)+numKaonBar(:)
       end if

       if (numLambda==1 .and. numLambdaBar==1) then
          sigma_Lambda_LambdaBar_X=sigma_Lambda_LambdaBar_X+1.
          if (numpart==2) then
            sigma_Lambda_LambdaBar=sigma_Lambda_LambdaBar+1.
          else if (numPions==numpart-2 .and. numPions <= 10) then
            sigma_Lambda_LambdaBar_Pions(numPions)=sigma_Lambda_LambdaBar_Pions(numPions)+1.
          else if (numpart >=4 .and. sum(numKaon(:))==1 .and.  sum(numKaonBar(:))==1) then
            sigma_Lambda_LambdaBar_KKbar_X=sigma_Lambda_LambdaBar_KKbar_X+1.
          else if (numpart >=4 .and. sum(numN(:))==1 .and.  sum(numNBar(:))==1) then
            sigma_Lambda_LambdaBar_NNbar_X=sigma_Lambda_LambdaBar_NNbar_X+1.
          end if
          !Particle_loop_Lambda_LambdaBar : do j = 1,size(realParticles,dim=2)
          !   if(realParticles(i,j)%ID <= 0) cycle Particle_loop_Lambda_LambdaBar
          !   write(*,*) 'Lambda_LambdaBar-event:',i,realParticles(i,j)%ID,realParticles(i,j)%charge
          !End Do Particle_loop_Lambda_LambdaBar
       end if

       PANDA1 : if (DoPanda) then

          if (numXi(-1)==1 .and. numXiBar(+1)==1) then
             sigma_Xi_XiBar_X=sigma_Xi_XiBar_X+1.
             if (numpart==2) then
                sigma_Xi_XiBar=sigma_Xi_XiBar+1.
             end if
          end if

          if (numXi(0)==1 .and. numXiBar(0)==1) then
             sigma_Xi0_Xi0Bar_X=sigma_Xi0_Xi0Bar_X+1.
             if (numpart==2) then
                sigma_Xi0_Xi0Bar=sigma_Xi0_Xi0Bar+1.
             end if
          end if

          !hyperon+nucleon scattering:
          if (numpart==2) then

             !Lambda(Sigma) + N channels:
             if (numSigma(0)==1) sigma_LpToS0p = sigma_LpToS0p + 1.
             if (numN(0)==1 .and. numLambda==1) sigma_SmpToLn = sigma_SmpToLn + 1.
             if (numN(0)==1 .and. numSigma(0)==1) sigma_SmpToS0n = sigma_SmpToS0n + 1.

             !Xi + N channels:

             !Xi^- p, Xi^0 n --> Lambda Lambda :
             if (numLambda==2) sigma_I0LL = sigma_I0LL + 1.

             !Xi^- p, Xi^0 n --> Lambda Sigma^0 :
             if (numLambda==1 .and. numSigma(0)==1) sigma_LS0 = sigma_LS0 + 1.

             !Xi^- p --> Xi^0 n :
             if (numXi(0)==1 .and. numN(0)==1) sigma_Xi0n = sigma_Xi0n + 1.

             !Xi^0 p --> Lambda Sigma^+ :
             if (numLambda==1 .and. numSigma(1)==1) sigma_I1SL = sigma_I1SL + 1.

             !Xi^- n --> Lambda Sigma^- :
             if (numLambda==1 .and. numSigma(-1)==1) sigma_LSm = sigma_LSm + 1.

          end if

          !K^- proton --> Xi+Kaon channels:
          if (numpart==2) then
             if (numXi(-1)==1 .and. numKaon(1)==1) sigma_KbarNuk_to_XiK(1) = sigma_KbarNuk_to_XiK(1) + 1.
             if (numXi(0)==1  .and. numKaon(0)==1) sigma_KbarNuk_to_XiK(2) = sigma_KbarNuk_to_XiK(2) + 1.
          end if
          !K^- proton --> Xi+Kaon+1pion channels:
          if (numpart==3) then
             if (numXi(-1)==1 .and. numKaon(0)==1 .and. numPi(1)==1) &
                  & sigma_KbarNuk_to_XiK(3) = sigma_KbarNuk_to_XiK(3) + 1.
             if (numXi(-1)==1 .and. numKaon(1)==1 .and. numPi(0)==1) &
                  & sigma_KbarNuk_to_XiK(4) = sigma_KbarNuk_to_XiK(4) + 1.
             if (numXi(0)==1 .and. numKaon(1)==1 .and. numPi(-1)==1) &
                  & sigma_KbarNuk_to_XiK(5) = sigma_KbarNuk_to_XiK(5) + 1.
          end if
          if (numXi(-1)==1) sigma_KbarNuk_to_XiMinus_X = sigma_KbarNuk_to_XiMinus_X + 1.


       end if PANDA1

       if ( (numLambda==1 .and. numSigmaBar(0)==1) .or. (numLambdaBar==1 .and. numSigma(0)==1) ) then
          sigma_Lambda_SigmaBar_cc = sigma_Lambda_SigmaBar_cc + 1.
       end if

       if (numPions==2) then
          sigma_2pion=sigma_2pion+1.
          if (numpart==2) then
             if (numPi(-1)==1.and.numPi(0)==1) then        ! pi^- pi^0
                sigma_pipi(1)=sigma_pipi(1)+1.
             else if (numPi(0)==2) then                    ! pi^0 pi^0
                sigma_pipi(2)=sigma_pipi(2)+1.
             else if (numPi(-1)==1.and.numPi(1)==1) then   ! pi^- pi^+
                sigma_pipi(3)=sigma_pipi(3)+1.
             else if (numPi(0)==1.and.numPi(1)==1) then    ! pi^0 pi^+
                sigma_pipi(4)=sigma_pipi(4)+1.
             end if
          end if
          if (numpart==3 .and. numEtas==1) sigma_2pion_Eta=sigma_2pion_Eta+1.
          if (numpart==4 .and. sum(numKaon(:))==1 .and. sum(numKaonBar(:))==1)&
             &sigma_2pion_KKbar=sigma_2pion_KKbar+1.
       end if

       !*** KKbar production in Nbar N annihilation: ************************************

       if ( sum(numKaon(:))==1 .and. sum(numKaonBar(:))==1 ) then

          if (numpart==2) then   ! K Kbar

             if (numKaon(0)==1 .and. numKaonBar(-1)==1) then
                sigma_KKbar(1)=sigma_KKbar(1)+1.           !  K^0 K^-
             else if (numKaon(1)==1 .and. numKaonBar(-1)==1) then
                sigma_KKbar(2)=sigma_KKbar(2)+1.           !  K^+ K^-
             else if (numKaon(0)==1 .and. numKaonBar(0)==1) then
                sigma_KKbar(3)=sigma_KKbar(3)+1.           !  K^0 Kbar^0
             else if (numKaon(1)==1 .and. numKaonBar(0)==1) then
                sigma_KKbar(4)=sigma_KKbar(4)+1.           !  K^+ Kbar^0
             end if

          else if (  numKaon(0)==1.and.numKaonBar(-1)==1  .or. &
                   &numKaon(1)==1.and.numKaonBar(0)==1 ) then

             if (numpart==3 .and. numPions==1) then
                sigma_KKbar(5)=sigma_KKbar(5)+0.5          ! pi^{+-} K^{-+} K^0_S   or   pi^0 K^- K^0_S
             else if (numpart==4 .and. numPions==2) then
                if (numPi(1)==1.or.numPi(-1)==1) then
                   sigma_KKbar(6)=sigma_KKbar(6)+0.5          ! pi^0 pi^{+-} K^{-+} K^0_S   or   pi^+ pi^- K^- K^0_S
                else if (numPi(-1)==2) then
                   sigma_KKbar(7)=sigma_KKbar(7)+0.5          ! pi^- pi^- K^+ K^0_S
                end if
             else if (numpart==5 .and. numPions==3) then
                if (numPi(0)<=1 .and. numPi(1)>=1) then
                   sigma_KKbar(8)=sigma_KKbar(8)+0.5          ! pi^{-+} pi^{+-} pi^{+-} K^{-+} K^0_S   or  pi^+ pi^0 pi^- K^- K^0_S
                else if (numPi(0)==1) then
                   sigma_KKbar(9)=sigma_KKbar(9)+0.5          ! pi^0 pi^- pi^- K^+ K^0_S
                end if
             else if (numpart==6 .and. numPions==4 .and. numPi(0)==1) then
                sigma_KKbar(10)=sigma_KKbar(10)+0.5          ! pi^{-+} pi^{+-} pi^{+-} pi^0 K^{-+} K^0_S
             end if

          else if (numKaon(0)==1 .and. numKaonBar(0)==1) then

             if (numpart==3 .and. numPions==1) then
                sigma_KKbar(11)=sigma_KKbar(11)+0.25                 ! pi^0 2K^0_S   or  pi^- 2K^0_S
             else if (numpart==4 .and. numPions==2 .and. numPi(-1)==1) then
                sigma_KKbar(12)=sigma_KKbar(12)+0.25                 ! pi^- pi^+ 2K^0_S   or  pi^0 pi^- 2K^0_S
             else if (numpart==5 .and. numPions==3 .and. numPi(0)<=1) then
                sigma_KKbar(13)=sigma_KKbar(13)+0.25               ! pi^- pi^0 pi^+ 2K^0_S   or   pi^+ pi^- pi^- 2K^0_S
             else if (numpart==6 .and. numPions==4) then
                if (numPi(0)==0) then
                   sigma_KKbar(14)=sigma_KKbar(14)+0.25               ! pi^- pi^- pi^+ pi^+ 2K^0_S
                else if (numPi(0)==2) then
                   sigma_KKbar(15)=sigma_KKbar(15)+0.25               ! pi^+ pi^0 pi^0 pi^- 2K^0_S
                end if
             end if

          else if (numKaon(1)==1 .and. numKaonBar(-1)==1) then

             if (numpart==3 .and. numPions==1) then
                sigma_KKbar(16)=sigma_KKbar(16)+1.                 ! pi^0 K^+ K^-   or   pi^- K^+ K^-
             else if (numpart==4 .and. numPions==2 .and. numPi(-1)==1) then
                sigma_KKbar(17)=sigma_KKbar(17)+1.                 ! pi^- pi^+ K^+ K^-   or  pi^- pi^0  K^+ K^-
             else if (numpart==5 .and. numPions==3 .and. numPi(0)<=1) then
                sigma_KKbar(18)=sigma_KKbar(18)+1.                 ! pi^- pi^0 pi^+ K^+ K^-  or  pi^+ pi^- pi^- K^+ K^-
             else if (numpart==6 .and. numPions==4 .and. numPi(0)==0) then
                sigma_KKbar(19)=sigma_KKbar(19)+1.                 ! pi^- pi^- pi^+ pi^+ K^+ K^-
             end if

          end if

       end if


       PANDA2 : if (DoPanda) then

          !Omega production channels
          if (numOmega==1) then
             !p + pbar --> Omega + OmegaBar (+X)
             !Kbar + Lambda --> Omega + K
             !Kbar + Sigma --> Omega + K
             !Kbar + Xi --> Omega + \pi (non-strange mesons)
             ! A+B --> \Omega^-
             if (numpart==2) sigma_OmegaBaryon = sigma_OmegaBaryon + 1.
             sigma_OmegaBaryonX = sigma_OmegaBaryonX + 1.

          end if

          !rescattering channels including Omega baryon (exclusive):
          if (numpart==2) then

             if ( (numXi(-1)==1 .or. numXi(0)==1) .and. &
                  & (numLambda==1 .or. numSigma(-1)==1 .or. numSigma(0)==1 .or. numSigma(1)==1) ) then
                sigma_OmegaToXi = sigma_OmegaToXi + 1.
             end if

             !Omega^- + proton channels:
             if (numXi(0)==1 .and. numLambda==1) sigma_OmegaP_to_Xi0Lambda = sigma_OmegaP_to_Xi0Lambda + 1.
             if (numXi(0)==1 .and. numSigma(0)==1) sigma_OmegaP_to_Xi0Sigma0 = sigma_OmegaP_to_Xi0Sigma0 + 1.
             if (numXi(-1)==1 .and. numSigma(1)==1) sigma_OmegaP_to_XimSigmap = sigma_OmegaP_to_XimSigmap + 1.

             !Omega^- + neutron channels:
             if (numXi(-1)==1 .and. numLambda==1) sigma_OmegaN_to_XimLambda = sigma_OmegaN_to_XimLambda + 1.
             if (numXi(-1)==1 .and. numSigma(0)==1) sigma_OmegaN_to_XimSigma0 = sigma_OmegaN_to_XimSigma0 + 1.
             if (numXi(0)==1 .and. numSigma(-1)==1) sigma_OmegaN_to_Xi0Sigmam = sigma_OmegaN_to_Xi0Sigmam + 1.

          else

             if ( (numXi(-1)==1 .or. numXi(0)==1) .and. &
                  & (numLambda==1 .or. numSigma(-1)==1 .or. numSigma(0)==1 .or. numSigma(1)==1) ) then
                sigma_OmegaToXiX = sigma_OmegaToXiX + 1.
             end if

             !Omega^- + proton channels:
             if (numXi(0)==1 .and. numLambda==1) sigma_OmegaP_to_Xi0LambdaX = sigma_OmegaP_to_Xi0LambdaX + 1.
             if (numXi(0)==1 .and. numSigma(0)==1) sigma_OmegaP_to_Xi0Sigma0X = sigma_OmegaP_to_Xi0Sigma0X + 1.
             if (numXi(-1)==1 .and. numSigma(1)==1) sigma_OmegaP_to_XimSigmapX = sigma_OmegaP_to_XimSigmapX + 1.

             !Omega^- + neutron channels:
             if (numXi(-1)==1 .and. numLambda==1) sigma_OmegaN_to_XimLambdaX = sigma_OmegaN_to_XimLambdaX + 1.
             if (numXi(-1)==1 .and. numSigma(0)==1) sigma_OmegaN_to_XimSigma0X = sigma_OmegaN_to_XimSigma0X + 1.
             if (numXi(0)==1 .and. numSigma(-1)==1) sigma_OmegaN_to_Xi0SigmamX = sigma_OmegaN_to_Xi0SigmamX + 1.

          end if

       end if PANDA2

       !#######################################################################
       if ( (.not.flag_not_only_pions .and. numBaryons.eq.0) .or. Do45ForAllEvents) then

          pion_events=pion_events+1.

          P_Npion(numPions)=P_Npion(numPions)+1.
!          P_Npion(numChPart)=P_Npion(numChPart)+1.

          Particle_loop2 : do j = 1,size(realParticles,dim=2)

             if (realParticles(i,j)%ID <= 0) cycle Particle_loop2

             if (realParticles(i,j)%ID .ne. pion) cycle Particle_loop2

             pPart => realParticles(i,j)

             momentum=pPart%momentum
             call lorentz(beta,momentum)

             momentumAbs=sqrt(dot_product(momentum(1:3),momentum(1:3)))
             ibin=nint((momentumAbs-dmom/2.)/dmom)+1
             if (ibin.ge.1 .and. ibin.le.Nmom .and. pPart%charge.ne.0) then
                dNpiondMom(ibin,0)=dNpiondMom(ibin,0)+1.
                if (numPions.ge.1 .and. numPions.le.10) &
                     & dNpiondMom(ibin,numPions)=dNpiondMom(ibin,numPions)+1.
             end if

             cosTheta=momentum(3)/momentumAbs
             call AddHist(dNpiondcosTheta,cosTheta,1.)

             Phii=atan2(momentum(2),momentum(1))
             call AddHist(dNpiondPhi,Phii,1.)

          end do Particle_loop2

       end if
       !#######################################################################

       !#######################################################################
       Particle_loop3 : do j = 1,size(realParticles,dim=2)
          if (realParticles(i,j)%ID <= 0) cycle Particle_loop3
          if (realParticles(i,j)%ID.eq.rho) call AddHist(dNrhodM,realParticles(i,j)%mass,1.)
          if (realParticles(i,j)%ID.eq.omegaMeson) call AddHist(dNomegadM,realParticles(i,j)%mass,1.)
       end do Particle_loop3
       !#######################################################################

       !#######################################################################

       if (DodNNbar) then
          if (numBaryons.gt.0) then
             do j = 1,size(realParticles,dim=2)
                if (realParticles(i,j)%ID <= 0) cycle
                if ( isBaryon(realParticles(i,j)%ID) .and. &
                     &realParticles(i,j)%antiParticle ) then
                   momentum=realParticles(i,j)%momentum
                   call lorentz(beta,momentum)
                   Ekin=momentum(0)-realParticles(i,j)%mass
                   call AddHist(dNNbardEkin(1),Ekin,1.)
                   if (flag_inel) then
                      call AddHist(dNNbardEkin(2),Ekin,1.)
                   else
                      call AddHist(dNNbardEkin(3),Ekin,1.)
                   end if
                   if (numpart.eq.3.and.numpions.eq.1) then
                      call AddHist(dNNbardEkin(4),Ekin,1.)
                   else if (numpart.eq.4.and.numpions.eq.2) then
                      call AddHist(dNNbardEkin(5),Ekin,1.)
                   end if
                end if
             end do
          end if
       end if

       !#######################################################################

       if (DoOutChannels) then
          if (CreateSortedPreEvent(event,PreEv%preE)) then
             PreEv%weight = 1.0
             call PreEvList_INSERT(ListPreEv,PreEv)
          end if
       end if

    end do Ensemble_loop


    fnorm=siggeo/float(numensembles)/float(isu)

    do i=-1,1

       if (DoH2d) then
          call WriteHist2D_Gnuplot(h2dMomentPion(i),141, add=1e-20,mul=fnorm,&
               &file='DoElementaryAnalysis6_p.'//piName(i)//'.'//intToChar(nwrite)//'.dat')
          call WriteHist2D_Gnuplot(h2dEnergyPion(i),141, add=1e-20,mul=fnorm,&
               &file='DoElementaryAnalysis6_E.'//piName(i)//'.'//intToChar(nwrite)//'.dat')
       end if

! T.G. --> this stuff produces 3000 output files!!!!! i commented out for the moment.
!       call WriteHist(hEnergyPion(i),141, add=1e-20,mul=fnorm,DoAve=.true.,&
!            &file='DoElementaryAnalysis7_E.'//piName(i)//'.'//intToChar(nwrite)//'.dat')
    end do

    if (finalFlag) then

       if (ncall==isu) then

          name(1) = PartName(particleID(1),particleCharge(1),particleAnti(1))
          name(2) = PartName(particleID(2),particleCharge(2),particleAnti(2))

          !********************************************************************
          open(31,file='DoElementaryAnalysis1.dat')
          write(31,'(A,I5,"[=>",A,"]")') '# Id of 1-st particle=             ', particleId(1),name(1)
          write(31,'(A,I5)')             '# charge of 1-st particle=         ', particleCharge(1)
          write(31,'(A,I5,"[=>",A,"]")') '# Id of 2-nd particle=             ', particleId(2),name(2)
          write(31,'(A,I5)')             '# charge of 2-nd particle=         ', particleCharge(2)
          write(31,'(A,E15.5)')          '# Geometrical cross section, mbarn:', siggeo
          write(31,'(A,I5)')             '# Number of ensembles:             ', numensembles
          write(31,'(A,I5)')             '# Number of runs at fixed energy:  ', isu
          write(31,'(6A)')&
               & '# Ek_lab:  p_lab:   srts:   tot:   elas:    CEX:     PROD_pbarp:   ANN:   Y:',&
               & '  pi-:   pi0:   pi+:   K+:    K-:    K0(ann,nonann):  Kbar0(ann,nonann):',&
               & '  LambdaLambdaBar:     LambdaLambdaBarX:',&
               & '  LambdaLambdaBarPi:  ... LambdaLambdaBar10Pi:',&
               & '  LambdaLambdaBarKKbarX:  LambdaLambdaBarNNbarX:',&
               & '  Str:  piLambda: piSigma^-:  piSigma^0:  piSigma^+:  sigma_prong:   sigma_Lambda_Sigma0Bar_cc'
          close(31)
          !********************************************************************
          open(31,file='DoElementaryAnalysis2.dat')
          write(31,'(A,I5,"[=>",A,"]")')  '# Id of 1-st particle=               ', particleId(1),name(1)
          write(31,'(A,I5)')              '# charge of 1-st particle=           ', particleCharge(1)
          write(31,'(A,I5,"[=>",A,"]")')  '# Id of 2-nd particle=               ', particleId(2),name(2)
          write(31,'(A,I5)')              '# charge of 2-nd particle=           ', particleCharge(2)
          write(31,'(A,(E15.5,1x))')      '# Geometrical cross section, mbarn:  ', siggeo
          write(31,'(A,I5)')              '# Number of ensembles:               ', numensembles
          write(31,'(A,I5)')              '# Number of runs at fixed energy:    ', isu
          write(31,'(3A)')                '# Elab:   p_lab:  srts:  KbarN:  NucNuc:  NucBarNuc:',&
                                         &' sigma_1pion_pbarp: sigma_2pion_pbarp:',&
                                         &' sigma_3pion_pbarp: sigma_4pion_pbarp:'
          close(31)
          !********************************************************************
          open(31,file='DoElementaryAnalysis3.dat')
          write(31,*)'# shift:', getShift0()
          close(31)
          !********************************************************************
          open(31,file='DoElementaryAnalysis4.dat')
          write(31,'(A,I5,"[=>",A,"]")')  '# Id of 1-st particle=               ', particleId(1),name(1)
          write(31,'(A,I5)')              '# charge of 1-st particle=           ', particleCharge(1)
          write(31,'(A,I5,"[=>",A,"]")')  '# Id of 2-nd particle=               ', particleId(2),name(2)
          write(31,'(A,I5)')              '# charge of 2-nd particle=           ', particleCharge(2)
          write(31,'(A,(E15.5,1x))')      '# Geometrical cross section, mbarn:  ', siggeo
          write(31,'(A,I5)')              '# Number of ensembles:               ', numensembles
          write(31,'(A,I5)')              '# Number of runs at fixed energy:    ', isu
          close(31)
          !********************************************************************
          open(31,file='DoElementaryAnalysis5.dat')
          write(31,'(A,I5,"[=>",A,"]")')  '# Id of 1-st particle=               ', particleId(1),name(1)
          write(31,'(A,I5)')              '# charge of 1-st particle=           ', particleCharge(1)
          write(31,'(A,I5,"[=>",A,"]")')  '# Id of 2-nd particle=               ', particleId(2),name(2)
          write(31,'(A,I5)')              '# charge of 2-nd particle=           ', particleCharge(2)
          write(31,'(A,(E15.5,1x))')      '# Geometrical cross section, mbarn:  ', siggeo
          write(31,'(A,I5)')              '# Number of ensembles:               ', numensembles
          write(31,'(A,I5)')              '# Number of runs at fixed energy:    ', isu
          close(31)
          !********************************************************************
          open(31,file='DoElementaryAnalysis6.dat')
          write(31,'(A,I5,"[=>",A,"]")')  '# Id of 1-st particle=               ', particleId(1),name(1)
          write(31,'(A,I5)')              '# charge of 1-st particle=           ', particleCharge(1)
          write(31,'(A,I5,"[=>",A,"]")')  '# Id of 2-nd particle=               ', particleId(2),name(2)
          write(31,'(A,I5)')              '# charge of 2-nd particle=           ', particleCharge(2)
          write(31,'(A,(E15.5,1x))')      '# Geometrical cross section, mbarn:  ', siggeo
          write(31,'(A,I5)')              '# Number of ensembles:               ', numensembles
          write(31,'(A,I5)')              '# Number of runs at fixed energy:    ', isu
          write(31,'(A)')                 '# Elab:   p_lab:  srts:  tot:   JPsi:'
          close(31)
          !********************************************************************

          if (Do2Part) then

             open(31,file='TwoPionChannels.dat')
             write(31,'(A,I5,"[=>",A,"]")') '# Id of 1-st particle=             ', particleId(1),name(1)
             write(31,'(A,I5)')             '# charge of 1-st particle=         ', particleCharge(1)
             write(31,'(A,I5,"[=>",A,"]")') '# Id of 2-nd particle=             ', particleId(2),name(2)
             write(31,'(A,I5)')             '# charge of 2-nd particle=         ', particleCharge(2)
             write(31,'(A,E15.5)')          '# Geometrical cross section, mbarn:', siggeo
             write(31,'(A,I5)')             '# Number of ensembles:             ', numensembles
             write(31,'(A,I5)')             '# Number of runs at fixed energy:  ', isu
             write(31,'(2A)')&
                  & '# Ek_lab:  p_lab:   srts:   sigma_2pion_incl:  sigma_2pion_eta:   sigma_2pion_KKbar:',&
                  & '  pi0 pi-:   pi0 pi0:    pi+ pi-:  pi+ pi0:'
             close(31)
             !*****************************************************************
             open(31,file='KKbarChannels.dat')
             write(31,'(A,I5,"[=>",A,"]")') '# Id of 1-st particle=             ', particleId(1),name(1)
             write(31,'(A,I5)')             '# charge of 1-st particle=         ', particleCharge(1)
             write(31,'(A,I5,"[=>",A,"]")') '# Id of 2-nd particle=             ', particleId(2),name(2)
             write(31,'(A,I5)')             '# charge of 2-nd particle=         ', particleCharge(2)
             write(31,'(A,E15.5)')          '# Geometrical cross section, mbarn:', siggeo
             write(31,'(A,I5)')             '# Number of ensembles:             ', numensembles
             write(31,'(A,I5)')             '# Number of runs at fixed energy:  ', isu
             write(31,'(A)')'# Column No.: Quantity:'
             write(31,'(A)')'# 1           Ek_lab'
             write(31,'(A)')'# 2           p_lab'
             write(31,'(A)')'# 3           srts'
             write(31,'(A)')'# 4           K^- K^0'
             write(31,'(A)')'# 5           K^- K^+'
             write(31,'(A)')'# 6           Kbar^0 K^0'
             write(31,'(A)')'# 7           Kbar^0 K^+'
             write(31,'(A)')'# 8           pi^{+-} K^{-+} K^0_S   or   pi^0 K^- K^0_S'
             write(31,'(A)')'# 9           pi^0 pi^{+-} K^{-+} K^0_S   or   pi^+ pi^- K^- K^0_S'
             write(31,'(A)')'# 10          pi^- pi^- K^+ K^0_S'
             write(31,'(A)')'# 11          pi^{-+} pi^{+-} pi^{+-} K^{-+} K^0_S   or  pi^+ pi^0 pi^- K^- K^0_S'
             write(31,'(A)')'# 12          pi^0 pi^- pi^- K^+ K^0_S'
             write(31,'(A)')'# 13          pi^{-+} pi^{+-} pi^{+-} pi^0 K^{-+} K^0_S'
             write(31,'(A)')'# 14          pi^0 2K^0_S   or  pi^- 2K^0_S'
             write(31,'(A)')'# 15          pi^- pi^+ 2K^0_S   or  pi^0 pi^- 2K^0_S'
             write(31,'(A)')'# 16          pi^- pi^0 pi^+ 2K^0_S   or   pi^+ pi^- pi^- 2K^0_S'
             write(31,'(A)')'# 17          pi^- pi^- pi^+ pi^+ 2K^0_S'
             write(31,'(A)')'# 18          pi^+ pi^0 pi^0 pi^- 2K^0_S'
             write(31,'(A)')'# 19          pi^0 K^+ K^-   or   pi^- K^+ K^-'
             write(31,'(A)')'# 20          pi^- pi^+ K^+ K^-   or  pi^- pi^0  K^+ K^-'
             write(31,'(A)')'# 21          pi^- pi^0 pi^+ K^+ K^-  or  pi^+ pi^- pi^- K^+ K^-'
             write(31,'(A)')'# 22          pi^- pi^- pi^+ pi^+ K^+ K^-'
             close(31)
             !*****************************************************************
          end if

          if (DoPanda) then

             open(31,file='Panda1.dat')
             write(31,'(A,I5,"[=>",A,"]")')  '# Id of 1-st particle=               ', particleId(1),name(1)
             write(31,'(A,I5)')              '# charge of 1-st particle=           ', particleCharge(1)
             write(31,'(A,I5,"[=>",A,"]")')  '# Id of 2-nd particle=               ', particleId(2),name(2)
             write(31,'(A,I5)')              '# charge of 2-nd particle=           ', particleCharge(2)
             write(31,'(A,(E15.5,1x))')      '# Geometrical cross section, mbarn:  ', siggeo
             write(31,'(A,I5)')              '# Number of ensembles:               ', numensembles
             write(31,'(A,I5)')              '# Number of runs at fixed energy:    ', isu
             write(31,'(A)')                '# Elab:   p_lab:  srts:  elastic:  inelastic:'
             close(31)
             !*****************************************************************
             open(31,file='Panda2.dat')
             write(31,'(A,I5,"[=>",A,"]")')  '# Id of 1-st particle=               ', particleId(1),name(1)
             write(31,'(A,I5)')              '# charge of 1-st particle=           ', particleCharge(1)
             write(31,'(A,I5,"[=>",A,"]")')  '# Id of 2-nd particle=               ', particleId(2),name(2)
             write(31,'(A,I5)')              '# charge of 2-nd particle=           ', particleCharge(2)
             write(31,'(A,(E15.5,1x))')      '# Geometrical cross section, mbarn:  ', siggeo
             write(31,'(A,I5)')              '# Number of ensembles:               ', numensembles
             write(31,'(A,I5)')              '# Number of runs at fixed energy:    ', isu
             write(31,'(A)')                '# Elab:   p_lab:  srts:  elastic:  inelastic:'
             close(31)

          end if

       end if

       !Normas:
       sigma_total=sigma_total*fnorm
       sigma_elastic=sigma_elastic*fnorm
       sigma_prong=sigma_prong*fnorm
       sigma_pion=sigma_pion*fnorm
       sigma_kaon=sigma_kaon*fnorm
       sigma_antikaon=sigma_antikaon*fnorm
       sigma_strangeness=sigma_strangeness*fnorm
       sigma_pion_Lambda=sigma_pion_Lambda*fnorm
       sigma_pion_Sigma=sigma_pion_Sigma*fnorm
       sigma_antikaon_nucleon=sigma_antikaon_nucleon*fnorm
       sigma_NucNuc=sigma_NucNuc*fnorm
       sigma_NucBarNuc=sigma_NucBarNuc*fnorm
       sigma_1pion_pbarp=sigma_1pion_pbarp*fnorm
       sigma_2pion_pbarp=sigma_2pion_pbarp*fnorm
       sigma_3pion_pbarp=sigma_3pion_pbarp*fnorm
       sigma_4pion_pbarp=sigma_4pion_pbarp*fnorm
       sigma_CEX=sigma_CEX*fnorm
       sigma_PROD_pbarp=sigma_PROD_pbarp*fnorm
       sigma_ANN=sigma_ANN*fnorm
       sigma_Y=sigma_Y*fnorm
       sigma_Lambda_LambdaBar=sigma_Lambda_LambdaBar*fnorm
       sigma_Lambda_LambdaBar_X=sigma_Lambda_LambdaBar_X*fnorm
       sigma_Lambda_LambdaBar_Pions=sigma_Lambda_LambdaBar_Pions*fnorm
       sigma_Lambda_LambdaBar_KKbar_X=sigma_Lambda_LambdaBar_KKbar_X*fnorm
       sigma_Lambda_LambdaBar_NNbar_X=sigma_Lambda_LambdaBar_NNbar_X*fnorm
       sigma_Lambda_SigmaBar_cc=sigma_Lambda_SigmaBar_cc*fnorm
       sigma_JPsi=sigma_JPsi*fnorm
       sigma_2pion=sigma_2pion*fnorm
       sigma_2pion_Eta=sigma_2pion_Eta*fnorm
       sigma_2pion_KKbar=sigma_2pion_KKbar*fnorm
       sigma_pipi=sigma_pipi*fnorm
       sigma_KKbar=sigma_KKbar*fnorm

       PANDA3 : if (DoPanda) then

          sigma_Xi_XiBar=sigma_Xi_XiBar*fnorm
          sigma_Xi_XiBar_X=sigma_Xi_XiBar_X*fnorm
          sigma_Xi0_Xi0Bar=sigma_Xi0_Xi0Bar*fnorm
          sigma_Xi0_Xi0Bar_X=sigma_Xi0_Xi0Bar_X*fnorm

          sigma_I0LL = fnorm*sigma_I0LL
          sigma_I1SL = fnorm*sigma_I1SL
          sigma_SmpToLn = sigma_SmpToLn*fnorm
          sigma_LpToS0p = sigma_LpToS0p*fnorm
          sigma_SmpToS0n = sigma_SmpToS0n*fnorm

          sigma_LS0 =  sigma_LS0*fnorm
          sigma_Xi0n = sigma_Xi0n*fnorm
          sigma_LSm = sigma_LSm*fnorm

          sigma_KbarNuk_to_XiK=sigma_KbarNuk_to_XiK*fnorm
          sigma_KbarNuk_to_XiMinus_X=sigma_KbarNuk_to_XiMinus_X*fnorm

          sigma_OmegaBaryon = sigma_OmegaBaryon*fnorm

          sigma_OmegaToXi = sigma_OmegaToXi*fnorm

          sigma_OmegaP_to_Xi0Lambda = sigma_OmegaP_to_Xi0Lambda*fnorm
          sigma_OmegaP_to_Xi0Sigma0 = sigma_OmegaP_to_Xi0Sigma0*fnorm
          sigma_OmegaP_to_XimSigmap = sigma_OmegaP_to_XimSigmap*fnorm
          sigma_OmegaN_to_XimLambda = sigma_OmegaN_to_XimLambda*fnorm
          sigma_OmegaN_to_XimSigma0 = sigma_OmegaN_to_XimSigma0*fnorm
          sigma_OmegaN_to_Xi0Sigmam = sigma_OmegaN_to_Xi0Sigmam*fnorm

          sigma_OmegaBaryonX = sigma_OmegaBaryonX*fnorm

          sigma_OmegaToXiX = sigma_OmegaToXiX*fnorm

          sigma_OmegaP_to_Xi0LambdaX = sigma_OmegaP_to_Xi0LambdaX*fnorm
          sigma_OmegaP_to_Xi0Sigma0X = sigma_OmegaP_to_Xi0Sigma0X*fnorm
          sigma_OmegaP_to_XimSigmapX = sigma_OmegaP_to_XimSigmapX*fnorm
          sigma_OmegaN_to_XimLambdaX = sigma_OmegaN_to_XimLambdaX*fnorm
          sigma_OmegaN_to_XimSigma0X = sigma_OmegaN_to_XimSigma0X*fnorm
          sigma_OmegaN_to_Xi0SigmamX = sigma_OmegaN_to_Xi0SigmamX*fnorm


       end if PANDA3


       if (pion_events.gt.0) then
          P_Npion(:)=P_Npion(:)/pion_events                  !   *fnorm/sigma_total
          dNpiondMom(:,:)=dNpiondMom(:,:)/pion_events/dmom   !   *fnorm/sigma_total/dmom
       end if

       open(31,file='DoElementaryAnalysis1.dat',position='append')
       write(31,5) ekin_lab, p_lab, srts, sigma_total, sigma_elastic, sigma_CEX, sigma_PROD_pbarp,&
            & sigma_ANN, sigma_Y,&
            & sigma_pion(-1:1), sum(sigma_kaon(1,:)), sum(sigma_antikaon(-1,:)),&
            & sigma_kaon(0,1),sigma_kaon(0,2),sigma_antikaon(0,1),sigma_antikaon(0,2),&
            & sigma_Lambda_LambdaBar,sigma_Lambda_LambdaBar_X,&
            & sigma_Lambda_LambdaBar_Pions(1:10),sigma_Lambda_LambdaBar_KKbar_X,&
            & sigma_Lambda_LambdaBar_NNbar_X,&
            & sigma_strangeness, sigma_pion_Lambda, sigma_pion_Sigma(-1:1),&
            & sigma_prong(0:N_max),sigma_Lambda_SigmaBar_cc
       close(31)

       PANDA4 : if (DoPanda) then

          !Lambda(Sigma)+nucleon rescattering:
          open(31,file='Panda1.dat',position='append')
          write(31,'(8e13.6)') ekin_lab, p_lab, srts, sigma_total, sigma_elastic, &
               & sigma_LpToS0p,sigma_SmpToLn,sigma_SmpToS0n
          close(31)

          !Xi+nucleon rescattering:
          open(31,file='Panda2.dat',position='append')
          write(31,'(10e13.6)') ekin_lab, p_lab, srts, sigma_total, sigma_elastic, &
               & sigma_I0LL,sigma_LS0,sigma_Xi0n, &  !Xi^- p (Xi^0 n) channels
               & sigma_LSm, & !Xi^- n
               & sigma_I1SL !Xi^0 p
          close(31)

          !Xi-production in ppBar collisions:
          open(31,file='Panda3.dat',position='append') !p pBar --> Xi XiBar(+X)
          write(31,'(9e13.6)') ekin_lab, p_lab, srts, sigma_total, sigma_elastic, &
               & sigma_Xi_XiBar_X,sigma_Xi_XiBar, &
               & sigma_Xi0_Xi0Bar_X,sigma_Xi0_Xi0Bar
          close(31)

          !Xi-production in Antikaon+Nucleon collisions:
          open(31,file='Panda4.dat',position='append')
          write(31,'(11e13.6)') ekin_lab, p_lab, srts, sigma_total, sigma_elastic, &
               sigma_KbarNuk_to_XiK,sigma_KbarNuk_to_XiMinus_X
          close(31)

          !Omega(S=-3)-production:
          open(31,file='Panda5.dat',position='append')
          write(31,'(11e13.6)') ekin_lab, p_lab, srts, sigma_total, sigma_elastic, &
               & sigma_OmegaBaryon, sigma_OmegaBaryonX
          close(31)

          !Omega(S=-3) rescattering:
          open(31,file='Panda6.dat',position='append')
          write(31,'(20e13.6)') ekin_lab, p_lab, srts, sigma_total, sigma_elastic, &
               & sigma_OmegaToXiX, sigma_OmegaToXi, &
               & sigma_OmegaP_to_Xi0Lambda, sigma_OmegaP_to_Xi0LambdaX, &
               & sigma_OmegaP_to_Xi0Sigma0, sigma_OmegaP_to_Xi0Sigma0X, &
               & sigma_OmegaP_to_XimSigmap, sigma_OmegaP_to_XimSigmapX, &
               & sigma_OmegaN_to_XimLambda, sigma_OmegaN_to_XimLambdaX, &
               & sigma_OmegaN_to_XimSigma0, sigma_OmegaN_to_XimSigma0X, &
               & sigma_OmegaN_to_Xi0Sigmam, sigma_OmegaN_to_Xi0SigmamX
          close(31)

       end if PANDA4

5      format(37(1x,f8.3),101(1x,e13.6),1x,f8.3)

       open(31,file='DoElementaryAnalysis2.dat',position='append')
       write(31,5) ekin_lab,p_lab,srts,sigma_antikaon_nucleon,sigma_NucNuc,sigma_NucBarNuc,&
                  &sigma_1pion_pbarp(0,-1:1),sigma_1pion_pbarp(1,-1:1),&
                  &sigma_2pion_pbarp,sigma_3pion_pbarp(-1:1),sigma_4pion_pbarp
       close(31)

       open(31,file='DoElementaryAnalysis3.dat',position='append')
       write(31,'(2(e13.6,1x))') srts - particleMass(1) - particleMass(2), sigma_total
       close(31)

       open(31,file='DoElementaryAnalysis4.dat',position='append')
       write(31,*)'# Elab: ', ekin_lab
       write(31,*)'# p_lab: ', p_lab
       write(31,*)'# srts: ', srts
       write(31,*)'# sigma_total: ', sigma_total
       write(31,*)'# sigma_elastic: ', sigma_elastic
       write(31,*)'# N:    P_Npion:'
       if (pion_events.gt.0) then
          pinumAv=0.
          do j=0,N_max
             pinumAv=pinumAv+float(j)*P_Npion(j)
             write(31,*) j, P_Npion(j)
          end do
          write(31,*)'# Norma: ', sum(P_Npion(:))
          write(31,*)'# Average pion number:', pinumAv
       else
          write(31,*)
          write(31,*) '>>>>>>>>>>>>>>>>>>>> ATTENTION: pion_events==0 !'
          write(31,*)
       end if
       close(31)

       open(31,file='DoElementaryAnalysis5.dat',position='append')
       write(31,*)'# Elab: ', ekin_lab
       write(31,*)'# p_lab: ', p_lab
       write(31,*)'# srts: ', srts
       write(31,*)'# sigma_total: ', sigma_total
       write(31,*)'# sigma_elastic: ', sigma_elastic
       write(31,*)'# momentum, GeV/c:    dNpiondMom, c/GeV:'
       if (pion_events.gt.0) then
          do j=1,Nmom
             write(31,'(f6.3,8(1x,e10.3))')  (float(j)-0.5)*dmom, dNpiondMom(j,0), &
                  & dNpiondMom(j,1:6), sum(dNpiondMom(j,7:10))
          end do
          write(31,*)'# Norma: ', sum(dNpiondMom(:,0))*dmom
       else
          write(31,*)
          write(31,*) '>>>>>>>>>>>>>>>>>>>> ATTENTION: pion_events==0 !'
          write(31,*)
       end if
       close(31)

       open(31,file='DoElementaryAnalysis6.dat',position='append')
       write(31,'(1x,3(f9.5,1x),e14.8,1x,e15.8)') ekin_lab, p_lab, srts, sigma_total, sigma_JPsi
       close(31)

       if (Dodsigdt) then

          call WriteHist(dsigdt_elastic,31,mul=fnorm,file='dsigdt_elastic.dat')
          call WriteHist(dsigdt_CEX,31,mul=fnorm,file='dsigdt_CEX.dat')
          call WriteHist(dsigdt_NucleonBar,31,mul=fnorm,file='dsigdt_NucleonBar.dat')
          call WriteHist(dsigdt_Delta,31,mul=fnorm,file='dsigdt_Delta.dat')
          call WriteHist(dsigdt_DeltaBar,31,mul=fnorm,file='dsigdt_DeltaBar.dat')
          call WriteHist(dsigdt_LambdaBarLambda,31,mul=fnorm,file='dsigdt_LambdaBarLambda.dat')
          call WriteHist(dsigdt_SigmaBarLambda,31,mul=fnorm,file='dsigdt_SigmaBarLambda.dat')
          call WriteHist(dsigdm_Delta,31,mul=fnorm,file='dsigdm_Delta.dat')
          call WriteHist(dsigdm_DeltaBar,31,mul=fnorm,file='dsigdm_DeltaBar.dat')
          call WriteHist(dNpiondcosTheta,31,file='dNpiondcosTheta.dat')
          call WriteHist(dNpiondPhi,31,file='dNpiondPhi.dat')
          call WriteHist(dNrhodM,31,file='dNrhodM.dat')
          call WriteHist(dNomegadM,31,file='dNomegadM.dat')
          call WriteHist(dNNbardEkin(1),31,file='dNNbardEkin_all.dat')
          call WriteHist(dNNbardEkin(2),31,file='dNNbardEkin_inel.dat')
          call WriteHist(dNNbardEkin(3),31,file='dNNbardEkin_el.dat')
          call WriteHist(dNNbardEkin(4),31,file='dNNbardEkin_1pi.dat')
          call WriteHist(dNNbardEkin(5),31,file='dNNbardEkin_2pi.dat')

       end if

       if (Do2Part) then

          open(31,file='TwoPionChannels.dat',position='append')
          write(31,5) ekin_lab, p_lab, srts, sigma_2pion, sigma_2pion_eta, sigma_2pion_KKbar,&
               &sigma_pipi(1:4)
          close(31)

          open(31,file='KKbarChannels.dat',position='append')
          write(31,'(22(1x,e13.6))') ekin_lab, p_lab, srts, sigma_KKbar(1:19)
          close(31)

       end if

       isu=0

       if (DoOutChannels) then
         open(141,file='OutChannels.'//intToChar(nwrite)//'.dat', status='unknown')
         rewind(141)
         call PreEvList_Print(141,ListPreEv,fnorm)
         close(141)
       end if

       nwrite=nwrite+1

    end if

  end subroutine DoElementaryAnalysis


end module ElementaryAnalysis
