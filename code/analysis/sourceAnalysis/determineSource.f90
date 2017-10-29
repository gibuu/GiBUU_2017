!******************************************************************************
!****m* /determineSource
! NAME
! module determineSource
!
! PURPOSE
! Administrates the determination of the source(s) which will
! undergo multifragmentation.
! NOTES
! * 3 Methods of source-selection added through the variable "SelectionMethod
!   controlled from NAMELIST (see sourceAnalysis_main.f90).
! * This step has been neccessary to performe the analysis in the same way
!   as in the experiment (ALADIN/INDRA experiment currently used).
! * This module works in both, relativistic RMF- and non-relativistic Skyrme-mode.
!
!******************************************************************************
module determineSource

  use sourceTypeDefinition

  private

  !****************************************************************************
  !****g* determineSource/TheSource
  ! SOURCE
  !
  type(quelle), allocatable, dimension(:,:), SAVE :: TheSource
  ! PURPOSE
  ! The variable "type(quelle) TheSource", in which the major properties
  ! are stored, e.g. mass,charge,energy,flow,position,velocity.
  !****************************************************************************

  !****************************************************************************
  !****g* determineSource/MaxNumSources
  ! SOURCE
  !
  integer,SAVE :: MaxNumSources=10
  ! PURPOSE
  ! Max. possible number of existing sources in actual timestep.
  !****************************************************************************

  !****************************************************************************
  !****g* determineSource/MaxSizeSources
  ! SOURCE
  !
  integer,SAVE :: MaxSizeSources=10
  ! PURPOSE
  ! Max. possible mass number of the existing sources in actual timestep.
  !****************************************************************************

  !****************************************************************************
  !****g* determineSource/divideFireball_Flag
  ! SOURCE
  !
  Logical,SAVE :: divideFireball_Flag=.false.
  ! PURPOSE
  ! If true, fireball source is divided into pieces using coalescence.
  !****************************************************************************

  !****************************************************************************
  !****g* determineSource/TeilchenPosX
  ! SOURCE
  !
  real,allocatable,dimension(:,:) ,SAVE :: TeilchenPosX
  ! PURPOSE
  ! Field to store x-positions of particles at time=0fm/c. Important only
  ! if SelectionMethod==2 (geometrical selection).
  !****************************************************************************



  public :: Get_FragmentingSource,deallocate_source,Get_InitialPosX
  public :: TheSource,MaxNumSources

contains

  !****************************************************************************
  !****s* determineSource/Get_FragmentingSource
  ! NAME
  ! subroutine Get_FragmentingSource(numEnsemples,numParticles,&
  !       &  SelectionMethod,betaChoice,hyperSource,&
  !       &  particleVector,Spectator_cutoff,rho_cutoff,A_cutoff,stossParameter, &
  !       &  Numsources,sourceType,&
  !       &  targetNuc,projectileNuc)
  ! PURPOSE
  ! * hadron-nucleus collisions: only one source (compound nucleus).
  ! * Heavy-Ion collisions:  three possibilities (fireball,spectators).
  !   Separation between different sources based on number of
  !   successfull collisions and on rapidity.
  ! INPUTS
  ! *   type(particle), dimension(:,:), intent(in) :: particleVector
  ! *   integer,                        intent(in) :: numEnsemples,numParticles
  ! *   integer,                        intent(in) :: SelectionMethod
  ! *   integer,                        intent(in) :: betaChoice
  ! *   logical,                        intent(in) :: hyperSource
  ! *   real,                           intent(in) :: Spectator_cutoff,rho_cutoff
  ! *   integer,                        intent(in) :: A_cutoff
  ! *   real,                           intent(in) :: stossParameter
  ! *   type(tNucleus),                 pointer    :: targetNuc
  ! *   type(tNucleus),optional,        pointer    :: projectileNuc
  ! OUTPUT
  ! *   integer, dimension(1:numEnsemples),                intent(out) :: NumSources
  ! *   integer, dimension(1:numEnsemples,1:numParticles), intent(out) :: sourceType
  !
  !****************************************************************************
  subroutine Get_FragmentingSource(numEnsemples,numParticles,&
       &  SelectionMethod,betaChoice,hyperSource,&
       &  particleVector,Spectator_cutoff,rho_cutoff,A_cutoff,stossParameter, &
       &  Numsources,sourceType,&
       &  targetNuc,projectileNuc)
    !**************************************************************************
    use constants, only: rhoNull
    use IdTable, only: nucleon,Lambda,SigmaResonance
    use nucleusDefinition
    use particleDefinition, only: particle
    use history, only: history_getGeneration
    use dichteDefinition, only: dichte
    use densitymodule, only: DensityAt
    use clus, only: DoClustering

    implicit none
    !---------------------------------------------------------------------
    !Input-Output variables
    !---------------------------------------------------------------------
    type(particle), dimension(:,:), intent(in) :: particleVector
    integer,                        intent(in) :: numEnsemples,numParticles
    integer,                        intent(in) :: SelectionMethod
    integer,                        intent(in) :: betaChoice
    logical,                        intent(in) :: hyperSource
    real,                           intent(in) :: Spectator_cutoff,rho_cutoff
    integer,                        intent(in) :: A_cutoff
    real,                           intent(in) :: stossParameter
    type(tNucleus),                 pointer    :: targetNuc
    type(tNucleus),optional,        pointer    :: projectileNuc

    integer, dimension(1:numEnsemples),                intent(out) :: NumSources
    integer, dimension(1:numEnsemples,1:numParticles), intent(out) :: sourceType
    !---------------------------------------------------------------------
    !Local variables
    !---------------------------------------------------------------------
    type(dichte)         :: localDensity
    real, dimension(1:3) :: place,velo
    integer              :: i,j,stossZahl
    real                 :: ypart
    real,    save        :: yproj=0.,ytarg=0.,rhoCut
    logical, save        :: startFlag=.true.

    !---------------------------------------------------------------------
    if (startFlag) then
       rhoCut = rhoNull/rho_cutoff
       if (present(projectileNuc)) then
          yproj      = 0.5*log((1.+ProjectileNuc%velocity(3))/(1.-ProjectileNuc%velocity(3)))
          ytarg      = 0.5*log((1.+TargetNuc%velocity(3))/(1.-TargetNuc%velocity(3)))
       end if
       startFlag = .false.
    end if


    !---------------------------------------------------------------------
    if (SelectionMethod <= 2) then

    EventClass : if (.not.present(projectileNuc) ) then !Hadron-Nucleus

       Loop_Ensemples1 : do i=1,numEnsemples
          Loop_Particles1 : do j=1,numParticles

             if (ParticleVector(i,j)%ID == 0) cycle Loop_Particles1
             if (ParticleVector(i,j)%ID < 0) exit Loop_Particles1

             sourceType(i,j) = 999 !indicates that particle does not belonge to source

             if (ParticleVector(i,j)%Antiparticle) cycle Loop_Particles1

             if (.not.hyperSource) then
                if (ParticleVector(i,j)%ID .ne. nucleon) cycle Loop_Particles1 !only nucleons!!!!
             else
                if (      ParticleVector(i,j)%ID.ne.nucleon &
                  & .and. ParticleVector(i,j)%ID.ne.Lambda &
                  & .and. ParticleVector(i,j)%ID.ne.SigmaResonance &
                  & .or.  ParticleVector(i,j)%ID.eq.SigmaResonance &
                  & .and. ParticleVector(i,j)%charge.ne.0     ) cycle !only nucleons,Lambda's, and Sigma0's
             end if

             place(1:3)   = ParticleVector(i,j)%position(1:3)
             localDensity = DensityAt(place)
             if ( localDensity%baryon(0) < rhoCut ) cycle

             call setParticlesToSources(1)

             !               sourceType(i,j) = 1 !target-source

          end do Loop_Particles1
       end do Loop_Ensemples1

       !in hadron-nucleus collisions:
       if (SelectionMethod==0) then !the source=one residual nucleus
          NumSources(:) = 1
          MaxNumSources = 1
       else !the source=two sources (residual target+moving source)
          NumSources(:) = 2
          MaxNumSources = 2
       end if
       MaxSizeSources= targetNuc%mass + 10

    else !Heavy-Ion

       Loop_Ensemples2 : do i=1,numEnsemples

          Loop_Particles2 : do j=1,numParticles

             if (ParticleVector(i,j)%ID == 0) cycle Loop_Particles2
             if (ParticleVector(i,j)%ID < 0) exit Loop_Particles2

             sourceType(i,j) = 999 !indicates that particle does not belonge to source(s)

             if (ParticleVector(i,j)%Antiparticle) cycle Loop_Particles2

             if (.not.hyperSource) then
                if (ParticleVector(i,j)%ID .ne. nucleon) cycle Loop_Particles2 !only nucleons!!!!
             else
                if (      ParticleVector(i,j)%ID.ne.nucleon &
                  & .and. ParticleVector(i,j)%ID.ne.Lambda &
                  & .and. ParticleVector(i,j)%ID.ne.SigmaResonance &
                  & .or.  ParticleVector(i,j)%ID.eq.SigmaResonance &
                  & .and. ParticleVector(i,j)%charge.ne.0     ) cycle !only nucleons,Lambda's, and Sigma0's
             end if

             place(1:3)   = ParticleVector(i,j)%position(1:3)
             localDensity = DensityAt(place)

             if ( localDensity%baryon(0) < rhoCut ) cycle

             StossZahl = history_getGeneration(particleVector(i,j)%history)

             velo(1:3) = particleVector(i,j)%velocity(1:3)
             ypart = 0.5*log((1+velo(3))/(1-velo(3))) !rapidity along z-axis

             call setParticlesToSources(2)

          end do Loop_Particles2

       end do Loop_Ensemples2

       !in HIC: (1) target source
       !        (2) projectile source
       !        (3) fireball source (optional: division into many smaller sources)
       call Fireball(numEnsemples,numParticles,particleVector,sourceType,NumSources)

       MaxNumSources = get_MaxNumber(NumEnsemples,NumSources)

    end if EventClass

    else if (SelectionMethod==3) then

       !spanning-tree criterion
       call DoClustering(numEnsemples,numParticles,particleVector,hyperSource,sourceType,NumSources)
       MaxNumSources = get_MaxNumber(NumEnsemples,NumSources)

    else

       write(*,*) ' wrong fragment-source selection method: ', SelectionMethod
       stop

    end if
    !---------------------------------------------------------------------

    allocate(TheSource(1:NumEnsemples,1:MaxNumSources))

    call extractSource

    if (present(projectileNuc) ) MaxSizeSources = get_MaxSize(NumEnsemples,NumSources)


  contains

    subroutine setParticlesToSources(iType)
      !*** selection of spectators(proj/targ) & participants
      !*** 3 Methods possible
      !*** For spectator fragmentation use only rapidity criterion!!!
      !    (cuts depend on experimental situation, here: ALADIN/INDRA-GSI)
      !*** by DEFAULT: SelectionMethod==0.
      implicit none

      integer, intent(in) :: iType
      real                :: veloz
      real                :: PosX,CutProj,CutTarg

      select case (iType)

      case (1) !Hadron-Nucleus

         if (SelectionMethod==0) then !One source in its rest frame
            sourceType(i,j) = 1
         else !more than one source (valid only at very high energies)
            veloz = ParticleVector(i,j)%momentum(3)/ParticleVector(i,j)%momentum(0)
            if (veloz > 0.5) then
               sourceType(i,j) = 2 !projectile-like source
            else
               sourceType(i,j) = 1 !target-source inb its rest frame
            end if
         end if

      case (2) !Heavy-Ion
         if (SelectionMethod==0) then
            !rapidity criterion only (ALADIN: y/yp > 0.75 for spectators)
            if ( ypart/yproj > 0.75 ) then
               sourceType(i,j) = 2 !projectile-source
            else if ( ypart/ytarg > 0.75 ) then
               sourceType(i,j) = 1 !target-source
            else
               sourceType(i,j) = 3 !fireball-source
            end if
         else if (SelectionMethod==1) then
            !collision criterion only
            if ( (float(stossZahl) < Spectator_cutoff) .and. &
                 &  (ypart > 0.) ) then
               sourceType(i,j) = 2 !projectile-source
            else if ( (float(stossZahl) < Spectator_cutoff) .and. &
                 &  (ypart < 0.) ) then
               sourceType(i,j) = 1 !target-source
            else if ( float(stossZahl) > Spectator_cutoff ) then
               sourceType(i,j) = 3 !fireball-source
            end if
         else if (SelectionMethod==2) then
            !geometrical selection
            CutProj = targetNuc%radius     + targetNuc%surface     - stossParameter/2.
            CutTarg = ProjectileNuc%radius + projectileNuc%surface - stossParameter/2.
            PosX = TeilchenPosX(i,j)
            if ( (abs(PosX) > CutProj) .and. (PosX > 0.) ) then
               sourceType(i,j) = 2 !projectile-source
            else if ( (abs(PosX) > CutTarg) .and.(PosX < 0.) ) then
               sourceType(i,j) = 1 !target-source
            else
               sourceType(i,j) = 3 !fireball-source
            end if
         else
            write(*,*) 'Module determineSource,routine setParticlesToSources:'
            write(*,*) 'Invalid input for SelectionMethod:',SelectionMethod
            write(*,*) 'Termination of the program'
            STOP
         end if

      end select

    end subroutine setParticlesToSources

    subroutine extractSource

      use lorentzTrafo, only: lorentz
      use IdTable, only: nucleon
      use potentialModule, only: trueEnergy !non-relativistic Skyrme
      use baryonPotentialModule, only: getPotentialEQSType
      use densitymodule, only: true4Momentum_RMF
      use coulomb, only: emfoca
      use RMF, only: getRMF_flag

      implicit none
      !-------------------------------------------------------------------
      ! Local variables
      !-------------------------------------------------------------------
      type(particle)                          :: teilchenLRF
      real,    dimension(1:3)                 :: Ort,impuls,beta
      real,    dimension(0:3)                 :: pin,trueMomentum
      real,    dimension(0:3,1:MaxNumSources) :: P_tot,Pstar_tot
      Real,    dimension(1:3,1:MaxNumSources) :: SourcePosition
      real,    dimension(1:MaxNumSources)     :: MassNorm,v_rad
      integer, dimension(1:MaxNumSources)     :: SMass,SCharge
      integer, dimension(1:MaxNumSources)     :: SLambda,SSigma0
      real                                    :: Ex,rsq,energy_rad,gamma_rad
      integer                                 :: i,m,l,isSource !,ii1,ii2,ii3
      logical                                 :: dummyFlag
      !-------------------------------------------------------------------
      ! determine average properties of fragmenting source (A,Z,Ex)
      !-------------------------------------------------------------------
      Loop_Ensemples3 : do m=1,numEnsemples

         SMass(1:NumSources(m))              = 0   !mass of the source(s)
         SCharge(1:NumSources(m))            = 0   !charge of the source(s)
         SLambda(1:NumSources(m))            = 0   !number of Lambda-hyperons in the source(s)
         SSigma0(1:NumSources(m))            = 0   !number of Sigma0-hyperons in the source(s)
         P_tot(0:3,1:NumSources(m))          = 0.0 !Ex. energy of the source(s)
         Pstar_tot(0:3,1:NumSources(m))      = 0.0 !total 4-momentum of the source(s)
         SourcePosition(1:3,1:NumSources(m)) = 0.0 !center-of-mass position of the source(s)
         MassNorm(1:NumSources(m))           = 0.0 !total mass (GeV) of the source(s)
         v_rad(1:NumSources(m))              = 0.0 !radial flow (in units of [c]) of the source(s)
         energy_rad = 0.0 !energy of radial flow (GeV/A) of the source(s)

         Loop_Particles3 : do l=1,numParticles

            if (ParticleVector(m,l)%ID == 0) cycle Loop_Particles3
            if (ParticleVector(m,l)%ID < 0) exit Loop_Particles3

            isSource = SourceType(m,l)

            if (isSource == 999) cycle

            Pstar_tot(0:3,isSource) = Pstar_tot(0:3,isSource) + &
                 &                    ParticleVector(m,l)%momentum(0:3)

            Ort(1:3)  = ParticleVector(m,l)%position(1:3)
            impuls(1:3) = ParticleVector(m,l)%momentum(1:3)

            if ( getRMF_Flag() ) then !RMF-mode (relativistic non-linear Walecka)
               call true4Momentum_RMF(ParticleVector(m,l),trueMomentum,dummyflag)
               P_tot(0,isSource)   = P_tot(0,isSource) + trueMomentum(0) &
                    + 0.5*emfoca(Ort,impuls,ParticleVector(m,l)%charge,ParticleVector(m,l)%ID)
               P_tot(1:3,isSource) = P_tot(1:3,isSource) + trueMomentum(1:3)
            end if

            SMass(isSource) = SMass(isSource) + 1
            if (ParticleVector(m,l)%ID==nucleon .and. ParticleVector(m,l)%charge==1) then
                  SCharge(isSource) = SCharge(isSource) + 1
            else if (ParticleVector(m,l)%ID==Lambda) then
                  SLambda(isSource) = SLambda(isSource) + 1
            else if (ParticleVector(m,l)%ID==SigmaResonance) then
                  SSigma0(isSource) = SSigma0(isSource) + 1
            end if


            if (isSource < 3) then
               SourcePosition(:,isSource) = SourcePosition(:,isSource) + &
                    & ParticleVector(m,l)%mass*Ort(:)
            else
               if (divideFireball_flag) then
                  SourcePosition(1,isSource) = SourcePosition(1,isSource) + &
                       & sqrt( Ort(1)**2 + Ort(2)**2 + Ort(3)**2 )*ParticleVector(m,l)%mass
                  SourcePosition(2:3,isSource) = 0.0
               else
                  SourcePosition(:,isSource) = SourcePosition(:,isSource) + &
                       & ParticleVector(m,l)%mass*Ort(:)
               end if
            end if

            MassNorm(isSource) = MassNorm(isSource) + ParticleVector(m,l)%mass

            if (present(projectileNuc)) then !make sence only for HeavyIon-runs.
               if (isSource .le. 2) then !spectators do not experience radial expansion.
                  teilchenLRF%momentum(1:3) = 0.0 !force it to 0
                  teilchenLRF%momentum(0) = ParticleVector(m,l)%mass !just to not divide below by 0
               else
                  teilchenLRF = ParticleVector(m,l)
               end if
            else
               teilchenLRF = ParticleVector(m,l)
            end if

            !average radial flow velocity:
            rsq     = sqrt( dot_product(Ort(1:3),Ort(1:3)) )
            v_rad(isSource) = v_rad(isSource) + &
                 &   dot_product(Ort(1:3),teilchenLRF%momentum(1:3))&
                 & / ( teilchenLRF%momentum(0)*rsq )

         end do Loop_Particles3
         !-------------------------------------------------------------------
         ! Non-Relativistic Skyrme mode:
         ! The total energy of the source is calculated after boosting all
         ! momenta into the source's rest frame
         !-------------------------------------------------------------------
         SkyrmeMode : if ( .not. getRMF_flag() ) then

            Loop_Particles4 : do l=1,numParticles

               if (ParticleVector(m,l)%ID == 0) cycle Loop_Particles4
               if (ParticleVector(m,l)%ID < 0) exit Loop_Particles4

               isSource = SourceType(m,l)

               if (isSource == 999) cycle

               !boost parameters
               beta(1:3) = Pstar_tot(1:3,isSource) / Pstar_tot(0,isSource)
               pin(0:3)  = ParticleVector(m,l)%momentum(0:3)
               call lorentz(beta,pin,'determineSource-1') !boost to source rest frame (SRF)
               teilchenLRF%ID       = ParticleVector(m,l)%ID
               teilchenLRF%charge   = ParticleVector(m,l)%charge
               teilchenLRF%position = ParticleVector(m,l)%position
               teilchenLRF%momentum = pin
               teilchenLRF%Mass     = ParticleVector(m,l)%mass

               if (getPotentialEQSType()==6) then !Birger's potential(no Coulomb)
                  P_tot(0,isSource)   = P_tot(0,isSource) + &
                       & getEnergyBirger(teilchenLRF)
               else !Skyrme+Coulomb
                  Ort(1:3)    = ParticleVector(m,l)%position(1:3)
                  impuls(1:3) = pin(1:3)
                  P_tot(0,isSource) = P_tot(0,isSource) + trueEnergy(teilchenLRF) &
                                      + 0.5*emfoca(Ort,impuls,ParticleVector(m,l)%charge,ParticleVector(m,l)%ID)
               end if
               P_tot(1:3,isSource) = 0.0

            end do Loop_Particles4

         end if SkyrmeMode

         !-------------------------------------------------------------------
         do i=1,NumSources(m)

            if (SMass(i).ge.A_cutoff) then

               v_rad(i) = v_rad(i) / float(SMass(i))

               gamma_rad = 1./sqrt(1.-v_rad(i)**2)
               energy_rad = MassNorm(i)*(gamma_rad-1.) / float(SMass(i))

               P_tot(:,i) = P_tot(:,i) / float(SMass(i))
               Pstar_tot(:,i) = Pstar_tot(:,i) / float(SMass(i))

               if ( getRMF_Flag() ) then
                  !boost parameters
                  if (betaChoice==0) then
                      beta(1:3) = Pstar_tot(1:3,i) / Pstar_tot(0,i)
                  else
                      beta(1:3) = P_tot(1:3,i) / P_tot(0,i)
                  end if
                  pin(0:3)  = P_tot(0:3,i)
                  call lorentz(beta,pin,'determineSource-2') !boost to rest frame
                  P_tot(:,i) = pin(:)
               end if

               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               !Note:You have to performe in addition a groundState run with
               !     the source parameters (A,Z) to obtain the corresponding
               !     (with respect to the considered mean-field model and
               !     with the same number of ensemples) binding energy and
               !     then extract the excitation energy.
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               Ex = P_tot(0,i) !Total energy of the source(s)(GeV/A)


               TheSource(m,i)%status      = .true.
               TheSource(m,i)%Size        = SMass(i)
               TheSource(m,i)%Charge      = SCharge(i)
               TheSource(m,i)%nLambda     = SLambda(i)
               TheSource(m,i)%nSigma0     = SSigma0(i)
               TheSource(m,i)%position(:) = (SourcePosition(:,i) / MassNorm(i))
               TheSource(m,i)%velocity(:) = beta(:)
               TheSource(m,i)%ExEnergy    = Ex
               TheSource(m,i)%radEnergy   = energy_rad

            else

               TheSource(m,i)%status      = .false.
               TheSource(m,i)%Size        = 0
               TheSource(m,i)%Charge      = 0
               TheSource(m,i)%nLambda     = 0
               TheSource(m,i)%nSigma0     = 0
               TheSource(m,i)%position(:) = 99999. !undefined
               TheSource(m,i)%velocity(:) = 99999. !undefined
               TheSource(m,i)%ExEnergy    = 99999. !undefined
               TheSource(m,i)%radEnergy   = 99999. !undefined

            end if

         end do

         ! Restore SourceType for particles which belong to the sources
         ! below cut-off mass
         do l=1,numParticles

            if (ParticleVector(m,l)%ID == 0) cycle
            if (ParticleVector(m,l)%ID < 0) exit

            isSource = SourceType(m,l)

            if (isSource == 999) cycle

            if (.not.TheSource(m,isSource)%status) SourceType(m,l)=999

          end do

      end do Loop_Ensemples3

    end subroutine ExtractSource

  !****************************************************************************
  end subroutine Get_FragmentingSource
  !****************************************************************************

  !****************************************************************************
  !****f* determineSource/getEnergyBirger
  ! NAME
  ! function getEnergyBirger
  !
  ! PURPOSE
  ! Calculates the total energy per nucleon in Birger's model.
  !
  !****************************************************************************
  Real Function getEnergyBirger(teilchen)

    use particleDefinition
    use nucDLDA, only: getEParticleLaplace
    use dichteDefinition
    use densitymodule, only: densityAt,gridSpacing
    use baryonPotentialModule, only: rhoLaplace

    implicit none

    type(particle), intent(in) :: Teilchen

    type(dichte) :: pos
    real, dimension(1:4) :: Etemp
    real, dimension(1:3) :: rvec
    real :: rhocent,rhoptemp,ptemp,rabstemp
    real, dimension(1:3), save :: stepsize=(/0.,0.,0./)
    logical, save :: local_init=.true.
    integer :: l

    if (local_init) then
       stepSize=gridSpacing
       local_init=.false.
    end if
    rvec=teilchen%position(1:3)
    pos=densityAt(rvec)
    rhocent=SQRT(pos%baryon(0)**2-Dot_Product(pos%baryon(1:3),pos%baryon(1:3)))
    rhoptemp=0.
    ptemp=0.
    rabstemp=0.
    do l=1,3,1
       ptemp=ptemp+teilchen%momentum(l)**2
       rabstemp=rabstemp+rvec(l)**2
    end do
    rabstemp=SQRT(rabstemp)
    rhoptemp=rhoLaplace(rvec,stepsize)
    call getEParticleLaplace(Etemp,rhocent,rhoptemp,ptemp)

    getEnergyBirger = Etemp(1) / 1000. !in units of GeV!!!

    !**************************************************************************
  end Function getEnergyBirger    !***************************************
  !****************************************************************************


  !****************************************************************************
  !****s* determineSource/Fireball
  ! NAME
  ! subroutine Fireball
  !
  ! PURPOSE
  ! divides fireball source into smaller clusters using coalescence. We choose
  ! the coalescence parameter in such way to avoid the appereance
  ! of free nucleons.
  ! NOTES
  ! The fireball-source/pieces should be as compact as possible, since
  ! the statistical multifragmentation model (SMM) simulates a statistical
  ! break-up (~coalescence), apart de-excitation.
  !
  !****************************************************************************
  subroutine Fireball(numEnsemples,numParticles,pv,sourceType,NumSources)
    use particleDefinition, only: particle
    use densitymodule, only: gridPoints,gridSpacing
    implicit none

    integer, intent(in) :: numEnsemples,numParticles
    type(particle), dimension(:,:), intent(in) :: pv
    integer, dimension(1:numEnsemples,1:numParticles), intent(inOut) :: sourceType
    integer, dimension(:), intent(out) :: NumSources

    integer :: i,j,i1,isSource
    integer, dimension(1:numEnsemples,1:numParticles) :: ifrm,fsource,Source_Save
    real, dimension(1:3) :: r1
    real :: rsq,rmin,rmax
    logical :: flag
    !---------------------------------------------------------------------
    ! Do nothing if divideFireball_Flag=.false.
    !---------------------------------------------------------------------
    if (.not.divideFireball_Flag) then
       NumSources(:) = 3
       return
    end if
    !---------------------------------------------------------------------
    ! initialization of local variables
    !---------------------------------------------------------------------
    ifrm(:,:)    = 0
    fsource(:,:) = 999
    Source_save(:,:) = SourceType(:,:)
    !---------------------------------------------------------------------
    ! Do coalescence. The fireball source is divided in radial cells of
    ! unique radial flow expansion energy.
    !---------------------------------------------------------------------
    Loop_Ensemples1 : do i=1,numEnsemples

       NumSources(i) = 3 !index 1-->target, index 2-->projectile

       Loop_over_space : do i1=1,gridPoints(1)

          rmin = float(i1-1)*gridSpacing(1)
          rmax = float(i1)*gridSpacing(1)

          flag = .false.

          Loop_Particles1 : do j=1,numParticles

             if (pv(i,j)%ID == 0) cycle Loop_Particles1
             if (pv(i,j)%ID < 0) cycle Loop_Particles1

             isSource = SourceType(i,j)
             if (isSource /= 3) cycle Loop_Particles1 !only fireball particles
             if (ifrm(i,j) == 1) cycle Loop_Particles1 !don't count twice

             flag = .false.

             r1(1:3) = pv(i,j)%position(1:3)
             rsq     = sqrt( r1(1)**2 + r1(2)**2 + r1(2)**2 )

             if (rsq > rmin .and. rsq <= rmax) then
                ifrm(i,j) = 1
                fsource(i,j) = NumSources(i)
                flag = .true.
             end if

          end do Loop_Particles1

          if (flag) NumSources(i) = NumSources(i) + 1

       end do Loop_over_space

    end do Loop_Ensemples1
    !---------------------------------------------------------------------
    ! update the variable "sourceType(:,:)".
    ! Reminder: in the variable "sourceType(:,:)" we store the iformation,
    ! in which source (target,projectile,fireball) each particle belonges.
    !---------------------------------------------------------------------
    Loop_Ensemples3 : do i=1,numEnsemples
       Loop_Particles3 : do j=1,numParticles

          if (pv(i,j)%ID == 0) cycle
          if (pv(i,j)%ID < 0) exit

          if (source_save(i,j) < 3) cycle !only fireball particles here!!!

          sourceType(i,j) = fsource(i,j)

       end do Loop_Particles3
    end do Loop_Ensemples3
    !---------------------------------------------------------------------

    !**************************************************************************
  end subroutine Fireball !***********************************************
  !****************************************************************************

  !****************************************************************************
  !****s* determineSource/get_MaxNumber
  ! NAME
  ! function get_MaxNumber
  !
  ! PURPOSE
  ! Exctract the possible max. values for number and size of sources.
  ! Needed for allocation of the "type(quelle) TheSource" and of
  ! the variables which characterizes the source (sourceProperties.f90 module).
  !
  !****************************************************************************
  function get_MaxNumber(NumEnsemples,NumSources) Result(MaxNumber)
    implicit none

    integer,                            intent(in) :: NumEnsemples
    integer, dimension(1:NumEnsemples), intent(in) :: NumSources
    integer :: i,MaxNumber

    !max. number of valid sources
    MaxNumber = 0
    do i=1,NumEnsemples
       if (MaxNumber < NumSources(i)) then
          MaxNumber = NumSources(i)
       end if
    end do

    !**************************************************************************
  end function get_MaxNumber !********************************************
  !****************************************************************************

  !****************************************************************************
  !****s* determineSource/get_MaxSize
  ! NAME
  ! function get_MaxSize
  !
  ! PURPOSE
  ! Exctract the possible max. values for number and size of sources.
  ! Needed for allocation of the "type(quelle) TheSource" and of
  ! the variables which characterizes the source (sourceProperties.f90 module).
  !
  !****************************************************************************
  function get_MaxSize(NumEnsemples,NumSources) Result(MaxSize)
    implicit none

    integer,                            intent(in) :: NumEnsemples
    integer, dimension(1:NumEnsemples), intent(in) :: NumSources
    integer :: i,j,MaxSize

    !max. size
    MaxSize = 0

    Loop_over_Ensemples : do i=1,NumEnsemples

       Loop_over_Sources : do j=1,NumSources(i)

          if (.not.TheSource(i,j)%status) cycle Loop_over_Sources

          if (MaxSize < TheSource(i,j)%Size) then
             MaxSize = TheSource(i,j)%Size
          end if

       end do Loop_over_Sources

    end do Loop_over_Ensemples

    !**************************************************************************
  end function get_MaxSize !**********************************************
  !****************************************************************************

  !****************************************************************************
  !****s* determineSource/deallocate_source
  ! NAME
  ! subroutine deallocate_source
  !
  ! PURPOSE
  ! deallocates the type(quelle) TheSource, if allocated.
  ! This is needed, since the size of the type(quelle) TheSource is
  ! variable for each time step.
  !
  !****************************************************************************
  subroutine deallocate_source
    implicit none

    if (allocated(TheSource)) deallocate(TheSource)

    !**************************************************************************
  end subroutine deallocate_source !**************************************
  !****************************************************************************

  !****************************************************************************
  !****s* determineSource/Get_InitialPosX
  ! NAME
  ! subroutine Get_InitialPosX
  !
  ! PURPOSE
  ! stores initial x-positions of particles into the field teilchenPosX.
  ! NOTES
  ! Important only if geometrical selection is used (SelectionMethod==2).
  !
  !****************************************************************************
  subroutine Get_InitialPosX(NumEns,NumPart,teilchen)
    use particleDefinition
    implicit none
    integer,                        intent(in) :: NumEns,NumPart
    type(particle), dimension(:,:), intent(in) :: teilchen

    allocate(teilchenPosX(1:numEns,1:numPart))
    teilchenPosX(:,:) = teilchen(:,:)%position(1)

    !**************************************************************************
  end subroutine Get_InitialPosX !****************************************
  !****************************************************************************


end module determineSource
