!*******************************************************************************
!****m* /smmModule
! NAME
! module smmModule
!
! PURPOSE
! This is the main module for the Statistical Multifragmentation Model (SMM).
!*******************************************************************************
module smmModule

  implicit none
  PRIVATE

  PUBLIC :: main_SMM

contains


  !*************************************************************************
  !****s* smmModule/main_SMM
  ! NAME
  ! Subroutine main_SMM
  !
  ! PURPOSE
  ! * This is the main subroutine for the SMM. 
  ! * Extracts important parameters from input needed for the SMM.
  ! * Calls the SMM-code.
  ! * Calls analysis tools for analysing the produced fragment vector.
  !*************************************************************************
  subroutine main_SMM

    use typeDefinitions, only : cluster,particle,quelle
    use inputGeneral,    only : SubEvents,NumEnsemples,& 
         &                      RealParticles,& 
         &                      ALADIN_Flag,CHARMS_Flag,Get_Hyp, Get_GiBUUvec
    use inputSMM
    use velocityFieldModule, only : Get_collectiveFlow,Get_RadialFlowProfile
    use fragmentation
    use fragmentationHyp
    use fragmentationAnalysis

    !***********************************************************************
    !****g* main_SMM/TheSource
    ! SOURCE
    !
    type(quelle),allocatable, dimension(:,:), SAVE :: TheSource
    !
    ! PURPOSE
    ! The source vector: target(first index), fireball(second index) and 
    ! (third index) spectator source properties are stored into the 
    ! type(quelle).
    !
    !***********************************************************************

    !***********************************************************************
    !****g* main_SMM/FragmentVector
    ! SOURCE
    !
    type(cluster), allocatable, dimension(:,:), SAVE :: FragmentVector
    !
    ! PURPOSE
    ! The Fragment vector: ground-state fragments and free nucleons.
    ! NOTES
    ! "Emitted" particles are not contained here 
    ! (see comments on ParticleVector just below)
    !
    !***********************************************************************

    !***********************************************************************
    !****g* main_SMM/ParticleVector
    ! SOURCE
    !
    type(particle), allocatable, dimension(:,:), SAVE :: ParticleVector
    !
    ! PURPOSE
    ! The particle vector: original GiBUU vector of real particles
    ! NOTES
    ! It is needed only for selecting the "emitted" particles, which are 
    ! not considered in the analysis of statistical multifragmentation.
    !
    !***********************************************************************

    !***********************************************************************
    !****g* main_SMM/NumSources
    ! SOURCE
    !
    integer, allocatable, dimension(:), SAVE :: NumSources
    !
    ! PURPOSE
    ! Number of sources for each MC-event.
    !
    !***********************************************************************

    integer :: i,j,m,iii
    integer :: Masse,Ladung,idummy,SourceType,TypeOfSource
    real    :: Exc,Erad,x,y,z,vx,vy,vz,bb,impactParameter
    real    :: Ebind,Etot,ProjectileVelo,rdummy

    real, dimension(1:2) :: A_res

    logical,save :: printFlag=.false.

    !-----------------------------------------------------------------------*
    call get_InputSMM
    !-----------------------------------------------------------------------*
    !Allocate major fields
    !-----------------------------------------------------------------------*
    allocate(FragmentVector(1:(SMM_Events),1:500))
    allocate(ParticleVector(1:NumEnsemples,1:RealParticles))
    allocate(NumSources(1:NumEnsemples))
    allocate(TheSource(1:NumEnsemples,1:MaxNumSources))
    !-----------------------------------------------------------------------*
    !Open-Read source-files
    !-----------------------------------------------------------------------*

    write(*,'(A)') ' --- Start reading SourceInfo: "'//trim(PathToSMMInput)//trim(SourceInfo)//'":'

    open(Unit=2, File=trim(PathToSMMInput)//trim(SourceInfo), & 
         & Status='old', Action='read')
    !-----------------------------------------------------------------------*
    ! Start of statistical fragmentation analysis for all subsequent runs
    !-----------------------------------------------------------------------*
    Loop_over_SubEvents : do m=1,SubEvents

       write(*,'(A,i5)') 'Subsequent run = ',m

       if (Get_GiBUUvec) then
          !set ParticleVector to default values:
          call setToDefault(1,FragmentVector,ParticleVector)
          !Get the particle vector from the GiBUU data files:
          call Get_ParticleVector(m,SubEvents,ParticleVector)
       endif

       if ( .not. ALADIN_Flag ) then

          !calculates the collective flow field
          !needed for boosting "fireball"-particles - produced by SMM-code - to
          !the frame, in which fireball expands with collective velocity field.
          if (EventType==1) then
             call Get_collectiveFlow(ParticleVector)
             call Get_RadialFlowProfile(m)
          endif

       endif

       if (.not.getExcBalance) then
          !get the binding energy of fragmenting source.
          !Hadron-Nucleus: Ebind from JobCard 
          !(extracted from additional BUU ground state run)
          if (EventType /= 1) then !proton-induced
             Ebind=E_bind_Input
          endif
       endif

       Loop_over_Ensemples : do i=1,NumEnsemples

          read(2,100) NumSources(i)

          Loop_over_Sources : do j=1,NumSources(i)

             !j=1,2: Spectators (target/projectile), otherwise: fireball
             read(2,*) rdummy,Masse,Ladung,x,y,z,vx,vy,vz,Etot,Erad,bb,idummy,idummy

             impactParameter = bb

             !get the binding energy of fragmenting source.
             !Heavy-Ion collisions: use of approximate values in 
             !function BindingEnergy
             if (EventType==1) then !HIC
                Ebind = BindingEnergy(Masse,Ladung)
             endif

             if (Erad > 0.010) then
                call ExcitationEnergy(Etot,Ebind,Exc,CorrectExc,Delta_Exc,Erad)
                SourceType = 3 !fireball
             else
                if (Masse > 2) then
                   if (.not. GetExcBalance) then 
                      call ExcitationEnergy(Etot,Ebind,Exc,CorrectExc,Delta_Exc)
                   else
                      if (j>1) then
                         write(*,*) 'Main_SMM: '
                         write(*,*) 'Do not use energy balance method if more...'
                         write(*,*) 'than one sources exist!!!'
                         write(*,*) 'STOP'
                         STOP
                      endif
                      A_res(1)=Masse-Ladung
                      A_res(2)=Ladung
                      call ExcFromBalance(i,beamID,beamEnergy,A_init,A_res, &
                                          vx,vy,vz,ParticleVector,Exc)
                   end if
                else
                   Exc = 0.0 !no SMM for deuterons!
                endif
                if (vz < 0.) then 
                   SourceType = 1 !target
                else
                   SourceType = 2 !projectile
                endif
             endif

!!$             Exc = -0.1
!!$             Erad = 0.0
!!$             if (abs(vz) < 0.3) then
!!$                SourceType = 3
!!$             else
!!$                if (vz < 0.) then 
!!$                   SourceType = 1 !target
!!$                else
!!$                   SourceType = 2 !projectile
!!$                endif
!!$             endif

             !Countes events as function of mass/charge number/Excitation energy
             !used only for Hadron-Nucleus and stability of ground-state
             if (EventType==2 .and. Hysto_Flag) & 
                  & call HistOfSource(m,i,SubEvents,NumEnsemples,Masse,Ladung,Exc)


!             write(*,'(A,4i4,2x,3f8.3,2x,3f8.3,2x,3f9.4)') & 
!                  & 'i,j,A,Z,X,V,Etot,Erad,Ex = ',& 
!                  & i,j,Masse,Ladung,x,y,z,vx,vy,vz,& 
!                  & (Etot-0.938)*1000.,Erad*1000.,Exc*1000.

             TheSource(i,j)%status      = .true.
             TheSource(i,j)%charge      = Ladung
             TheSource(i,j)%Size        = Masse
             TheSource(i,j)%ExEnergy    = Exc
             TheSource(i,j)%radEnergy   = Erad
             TheSource(i,j)%Type        = SourceType
             TheSource(i,j)%velocity(1) = vx
             TheSource(i,j)%velocity(2) = vy
             TheSource(i,j)%velocity(3) = vz
             TheSource(i,j)%position(1) = x
             TheSource(i,j)%position(2) = y
             TheSource(i,j)%position(3) = z

          end do Loop_over_Sources

       end do Loop_over_Ensemples


       Loop_over_Ensemples2 : do i=1,NumEnsemples

          write(*,'(A,i5)') 'NumEnsemple = ',i

          !set FragmentVector to default values
          call setToDefault(2,FragmentVector,ParticleVector)

          !call main routine for SMM code
          write(*,*) 'calling SMM code...'
          call MultiFragmentation(SMM_Seed,SMM_Events,m,i,NumSources,& 
               & EventType,SMM_Flag,ALADIN_Flag,TheSource,FragmentVector)

          !call routine for coalescence between SMM-fragments & hyperons
          if (Get_Hyp) &
               call MultiFragmentationHyp(i,TheSource, &
                                          FragmentVector,ParticleVector)

          if (m==SubEvents .and. i==NumEnsemples) printFlag=.true.

          !print FragmentVector
          if (printFV_Flag) call printFragmentVector(m,FragmentVector)

          !analyze FragmentVector
          write(*,*) 'calling analysis routines...'

          if (Get_Hyp) &
               call detailedHypAnalysis(m,i,SubEvents,NumEnsemples,SMM_Events, &
                                        FragmentVector,particleVector,printflag)

          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          !~~~~~~~~~~~~~~~~~ALADIN/INDRA/CHARMS data analyzed ~~~~~~~~~~~~~~
          !
          ALADINorFOPI : if (ALADIN_Flag) then !Spectator fragmentation
          !
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


             call spectatorFragmentation(SubEvents,NumEnsemples,SMM_Events, & 
                  & i,FragmentVector,ParticleVector,printflag,nint(A_init(2)))
             call VariousDistributions(SubEvents,numEnsemples,SMM_Events, & 
                  & i,FragmentVector,ParticleVector,XiTrigger,printFlag,impactParameter)

             if (CHARMS_Flag) then
                do iii=1,MaxNumSources
                   if ( .not. TheSource(i,iii)%status ) cycle
                   if (TheSource(i,iii)%Type /= 2 ) CYCLE
                   ProjectileVelo = TheSource(i,iii)%velocity(3)
                   TypeOfSource   = TheSource(i,iii)%Type
                end do
                call velocityDistr_Charms(SubEvents,NumEnsemples,SMM_Events,& 
                     & FragmentVector,printflag, & 
                     & ProjectileVelo,TypeOfSource)
                call CharmsAnalysis(FragmentVector,printflag)
             endif

          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          else !FOPI data analyzed (Central collisions)
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                call VariousDistributions(SubEvents,numEnsemples,SMM_Events, & 
                     & i,FragmentVector,ParticleVector,XiTrigger,printFlag,impactParameter)
                call rapidityDistribution(SubEvents,numEnsemples,SMM_Events, & 
                     & i,FragmentVector,ParticleVector,printFlag)

          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          endif ALADINorFOPI
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          !Distributions versus excitation energy "Ex" (MeV) of initial residual nucleus
          !(only for hadron-induced reactions):
          if (EventType==2) & 
               & call ExDistributions(SubEvents,numEnsemples,SMM_Events, & 
               & i,FragmentVector,ParticleVector,TheSource,printFlag,impactParameter)

          !for completeness: spectra & Co. of (emitted) mesons:
          if (Get_GiBUUvec) & 
                call producedParticles(SubEvents,numEnsemples, &
                                       i,ParticleVector, &
                                       printFlag,impactParameter)

       end do Loop_over_Ensemples2


    end do Loop_over_SubEvents
    !-----------------------------------------------------------------------*

    close(unit=2,status='keep')

    deallocate(FragmentVector)
    deallocate(ParticleVector)
    deallocate(NumSources)
    deallocate(TheSource)


100 format(25x,1x,i5)

  end subroutine main_SMM !*************************************************


  !*************************************************************************
  !****s* smmModule/BindingEnergy
  ! NAME
  ! Function BindingEnergy
  !
  ! PURPOSE
  ! determines the binding energy of light nuclei
  ! NOTES
  ! * For light nuclei we use exp. values, otherwise RMF values. This is 
  !   because RMF does not give reasonable results for binding energies 
  !   for light clusters! This problem has to be resolved in the future 
  !   using a better groundState determination in GiBUU runs.
  !*************************************************************************
  real function BindingEnergy(A,Z)
    integer, intent(in) :: A,Z
    real :: Ebind
    !-----------------------------------------------------------------------*
    ! RMF-Values also for very light nuclei
    ! Note: These values are not realistic, however, they must be used, in 
    !       order to consistently calculate the excitation energy.
    !-----------------------------------------------------------------------*
    EBind = 11.5
    !Fragments with A=2 (deuteron) is assumed to be always bound.
    if (A==3) Ebind = 5.5
    if (A==4) Ebind = 7.5
    if (A==5) Ebind = 7.
    if (A==6) Ebind = 7.4
    if (A==7) Ebind = 8.
    if (A==8) Ebind = 8.3
    if (A==9) Ebind = 8.6
    if (A==10) Ebind = 9.
    if (A==11) Ebind = 9.3
    if (A==12) Ebind = 9.4
    if (A==13) Ebind = 9.95
!!$    !-----------------------------------------------------------------------*
!!$    ! exp-Values for very light nuclei, otherwise mean RMF value (-11.5MeV/A)
!!$    !-----------------------------------------------------------------------*
!!$    EBind = 11.5
!!$    if (A==2 .and. Z==1) Ebind = 1.112
!!$    if (A==3) Ebind = 2.65
!!$    if (A==4) then
!!$       Ebind = 3.21
!!$       if (Z==1) Ebind = 1.4
!!$       if (Z==2) Ebind = 7.074
!!$       if (Z==3) Ebind = 1.154
!!$    endif
!!$    if (A==5) then
!!$       Ebind = 4.03
!!$       if (Z==1) Ebind = 1.336
!!$       if (Z==2) Ebind = 5.481
!!$       if (Z==3) Ebind = 5.266
!!$    endif
!!$    if (A==6) then
!!$       Ebind = 3.16
!!$       if (Z==1) Ebind = 0.964
!!$       if (Z==2) Ebind = 4.878
!!$       if (Z==3) Ebind = 5.332
!!$       if (Z==4) Ebind = 4.487
!!$       if (Z==5) Ebind = 0.152
!!$    endif
!!$    if (A==7) then
!!$       Ebind = 4.0
!!$       if (Z==1) Ebind = 0.940
!!$       if (Z==2) Ebind = 4.119
!!$       if (Z==3) Ebind = 5.606
!!$       if (Z==4) Ebind = 5.371
!!$       if (Z==5) Ebind = 3.531
!!$    endif
!!$    if (A==8) then
!!$       Ebind = 4.8
!!$       if (Z==2) Ebind = 3.926
!!$       if (Z==3) Ebind = 5.159
!!$       if (Z==4) Ebind = 7.062
!!$       if (Z==5) Ebind = 4.717
!!$       if (Z==6) Ebind = 3.098
!!$    endif
!!$    if (A==9) then
!!$       Ebind = 5.05
!!$       if (Z==2) Ebind = 3.348
!!$       if (Z==3) Ebind = 5.038
!!$       if (Z==4) Ebind = 6.463
!!$       if (Z==5) Ebind = 6.257
!!$       if (Z==6) Ebind = 4.337
!!$    endif    
!!$    if (A==10) then
!!$       Ebind = 5.0
!!$       if (Z==2) Ebind = 3.033
!!$       if (Z==3) Ebind = 4.531
!!$       if (Z==4) Ebind = 6.498
!!$       if (Z==5) Ebind = 6.475
!!$       if (Z==6) Ebind = 6.032
!!$       if (Z==7) Ebind = 3.644
!!$    endif
    !-----------------------------------------------------------------------*
    Ebind = Ebind / 1000.
    BindingEnergy = Ebind !binding energy [GeV/A]

  end function BindingEnergy !**********************************************


  !*************************************************************************
  !****s* smmModule/ExcFromBalance
  ! NAME
  ! Subroutine ExcFromBalance
  !
  ! PURPOSE
  ! Calculates the excitation energy of the source(s) using energy balance 
  ! event-by-event
  !*************************************************************************
  subroutine ExcFromBalance(inum,beamID,beamEnergy,A_init,A_res, &
                            vx,vy,vz,PV,Exc)
    use TypeDefinitions, only : particle
    !-----------------------------------------------------------------------
    ! input/output:
    !-----------------------------------------------------------------------
    type(particle), dimension(:,:), intent(in) :: PV
    real, dimension(1:2), intent(in) :: A_init,A_res
    integer, intent(in) :: inum,beamID
    real, intent(in) :: beamEnergy,vx,vy,vz
    real, intent(out) :: Exc
    !-----------------------------------------------------------------------
    ! local variables
    !-----------------------------------------------------------------------
    real, parameter :: M_nuc=0.938
    integer :: j
    real    :: ekin_beam, ekin_part, BE1, BE2, vsq, ekin_coll, M_beam

    !-----------------------------------------------------------------------
    Select Case(abs(beamID))
    Case(1)
       M_beam=0.938
    Case(53)
       M_beam=1.315
    Case default
       write(*,*) 'from MainSMM/ExcFromBalance:'   
       write(*,*) 'not defined option for beamID. beamID = ',beamID
       write(*,*) 'STOP'
       STOP
    end Select
    !-----------------------------------------------------------------------
    ! energy of the beam:
    ekin_beam = beamEnergy+M_beam

    ! energy of emitted particles: 
    ekin_part = 0.0
    do j=1,size(PV,dim=2)
       ekin_part = ekin_part + PV(inum,j)%momentum(0) 
    end do

    !get binding energy from Bethe-Weizsaecker:
    BE1 = BWformula(A_init)/Sum(A_init(:))
    BE2 = BWformula(A_res)/Sum(A_res(:))

    !get recoil: 
    vsq = vx**2+vy**2+vz**2
    ekin_coll = vsq*(M_nuc-BE2)*Sum(A_res(:))/2.

    !excitation energy (GeV):
    Exc = ekin_beam + (M_nuc-BE1)*Sum(A_init(:))-(M_nuc-BE2)*Sum(A_res(:)) & 
      & - ekin_part - ekin_coll

    !exc. energy per nucleon (GeV/A):
    Exc = Exc / Sum(A_res(:))

!    write(*,'(a,3f12.4)') '*** ',Sum(A_res(:)), Exc*1000.
!    write(*,*)


    if (Exc < 0.0) then
       Exc = 0.0
    endif

!    write(*,'(7f12.4)') Sum(A_init(:)), Sum(A_res(:)), ekin_beam*1000., ekin_part*1000., ekin_coll*1000., BE1*1000., BE2*1000.
!    write(150,'(7f12.4)') Sum(A_init(:)), Sum(A_res(:)), ekin_beam*1000., ekin_part*1000., ekin_coll*1000., BE1*1000., BE2*1000.

  contains

    real function BWformula(A)
      real, dimension(1:2), intent(in) :: A
      real, parameter :: av=15.67, a0=17.23, ac=0.714, as=23.3 !MeV
      real :: M,N,Z

      M = A(1)+A(2)
      N = A(1)
      Z = A(2)

      BWformula = av*M - a0*M**(2./3.) - ac*Z*(Z-1.)*M**(-0.33333333) - & 
           & as*(N-Z)**2/M

      BWformula = abs(BWformula/1000.)

    end function BWformula

  !*************************************************************************
  end subroutine ExcFromBalance !*******************************************
  !*************************************************************************


  !*************************************************************************
  !****s* smmModule/ExcitationEnergy
  ! NAME
  ! Subroutine ExcitationEnergy
  !
  ! PURPOSE
  ! Calculates the excitation energy of the source(s).
  ! NOTES
  ! * For spectator sources no radial flow by definition.
  ! * For fireball sources we have to subtract from the total energy the 
  !   average energy of radial expansion, which enters also as input 
  !   in the SMM code.
  ! * However, due to strong radial dependence of the radial flow this 
  !   prescription hase to be regarded as an approximation.
  ! * The upper limit in the excitation energy of fireballs is neccessary 
  !   since the SMM code terminates otherwise (A.Botvina).
  ! * For fireball the corrected radial flow energy is calculated as well.
  !*************************************************************************
  subroutine ExcitationEnergy(Etot,Ebind,Exc,CorrectExc,Delta_Exc,Erad)
    real,           intent(in)    :: Etot,Ebind,Delta_Exc
    logical,        intent(in)    :: CorrectExc
    real,           intent(out)   :: Exc
    real, optional, intent(inOut) :: Erad
    !-----------------------------------------------------------------------
    real :: Exc_new,Erad_new,Ediff
    !-----------------------------------------------------------------------
    Exc  = Etot - 0.938 + Ebind  !reference energy, Ebind < 0 !!!

    !Exc_max ~ 0.025 GeV/A !!!!at higher excitation energies problems with SMM
    !(Botvina: use a lower value for Exc and distribute the rest energy into
    ! radial flow, in any case at so high excitation values
    !the system always dissintegrates into free nucleons...
    if (present(Erad)) then
       Exc_new = Exc
       Erad_new= Erad
       Ediff = Etot - 0.938 + Ebind - Erad
       if (Ediff > 25./1000.) then
          Exc_new = Etot - 0.938 + Ebind - (Ediff-25./1000.)
          Erad_new= Erad + (Ediff-25./1000.)
          if (Exc_new < 0. .or. Erad_new < 0.) then
             write(*,*) 'smmModule, ExcitationEnergy: '
             write(*,*) 'wrong determination of Exc/Erad!!!',& 
                  & Ediff,Exc,Erad,Exc_new,Erad_new
          endif
       endif
       Exc = Exc_new
       Erad= Erad_new
    endif

    !Correction due to violation on total energy (small increase of 
    !binding energy of ground-state nucleus).
    !The value for DeltaExc is taken from JobCard (GeV per nucleon!!!)
    if (CorrectExc) Exc = Exc - Delta_Exc !in units of GeV per nucleon

    if (Exc < 0.0) then
       Exc = 0.0
    endif


  !*************************************************************************
  end subroutine ExcitationEnergy  !****************************************
  !*************************************************************************

  !*************************************************************************
  !****s* smmModule/setToDefault
  ! NAME
  ! Subroutine setToDefault
  !
  ! PURPOSE
  ! * Initializes the fragment and particle vector.
  !*************************************************************************
  subroutine setToDefault(switch,FV,PV)
    use TypeDefinitions, only : cluster,particle

    integer,                        intent(in)    :: switch
    type(cluster),  dimension(:,:), intent(inout) :: FV
    type(particle), dimension(:,:), intent(inout) :: PV

    integer :: i,j

    if (switch==1) then !reset ParticleVector

       Loop_over_ensembles2: do i = 1,size(PV,dim=1)
          Loop_over_particles2 : do j = 1,size(PV,dim=2)

             PV(i,j)%number       = 0
             PV(i,j)%bornTime     = -999.
             PV(i,j)%lastCollTime = -999.
             PV(i,j)%collHis      = -999
             PV(i,j)%position     = 0.
             PV(i,j)%momentum     = 0.
             PV(i,j)%mass         = 0.
             PV(i,j)%ID           = 999
             PV(i,j)%charge       = 0
             PV(i,j)%event        = 0
             PV(i,j)%ensemple     = 0

          end do Loop_over_particles2
       end do Loop_over_ensembles2

    else !reset FragmentVector

       Loop_over_ensembles1: do i = 1,size(FV,dim=1)
          Loop_over_particles1 : do j = 1,size(FV,dim=2)

             FV(i,j)%position     = 0.
             FV(i,j)%momentum     = 0.
             FV(i,j)%mass         = 0.
             FV(i,j)%ID           = 0
             FV(i,j)%MassNumber   = 0
             FV(i,j)%ChargeNumber = 0
             FV(i,j)%HypNumber    = 0
             FV(i,j)%HypType      = '-'
             FV(i,j)%FreeBound    = .true.
             FV(i,j)%HypContent   = -999

          end do Loop_over_particles1
       end do Loop_over_ensembles1

    endif

  end subroutine setToDefault !*********************************************


  !***************************************************************************
  !****s* smmModule/Get_ParticleVector
  ! NAME
  ! Subroutine Get_ParticleVector
  !
  ! PURPOSE
  ! * Reads the (real) GiBUU particle vector.
  !***************************************************************************
  subroutine Get_ParticleVector(isu,SubEvents,ParticleVector)
    use TypeDefinitions, only : particle
    use inputGeneral,    only : NumEnsemples,pathToBUUInput, BUU_DataFile

    type(particle), dimension(:,:), intent(inout) :: ParticleVector
    integer, intent(in) :: isu,SubEvents

    integer :: j,m
    integer :: ios,ior,ioe,idp,iso,ActualEvent,subev
    integer :: GlobalID, collHis
    real    :: mass,x,y,z,px,py,pz,p0,prodTime,lastCollTime

    logical, SAVE :: openFlag=.true.

    integer,save :: Read_nEns, Read_Format
    integer :: mmax

!    integer :: mtest

    Read_nEns = -1
!    Read_Format = 1 ! = DetailedHyperonOutput = .false.
    Read_Format = 2 ! = DetailedHyperonOutput = .true.
    mmax = 0

    if (openFlag) then
       open(Unit=200, File=trim(PathToBUUInput)// BUU_DataFile, & 
            & Status='old', Action='read',Iostat=ios)

       if (ios /= 0) then
          write(*,*) 'Get_ParticleVector: BUU_DataFile Open failed: ios = ',ios
          write(*,*) 'Get_ParticleVector: !!! Termination of program NOW !!!'
          STOP
       endif
       openFlag = .false.

       write(*,'(A,A,A)') 'BUU file "',trim(PathToBUUInput)//trim(BUU_DataFile),'" opened...'

       read(200,*,iostat=ior) Read_nEns, Read_Format
       if (ior.ne.0) then
          write(*,*) 'No Information about nEns or Format found in file...'
          backspace(Unit=200)
       else
          if (Read_nEns.ne.NumEnsemples) then
             write(*,*) 'mismatch for #Ensembles: ',Read_nEns,NumEnsemples
             stop
          end if
       end if

    endif

    if (isu==1) then
       read(200,*)
       read(200,*)
    endif
   
!    if (.not.openFlag) backspace(Unit=1) !move one line back!!!!

!    mtest = 0

   ENSEMPLES : do j=1,NumEnsemples
      m=0
      Particles : do 


         select case(Read_Format)
         case(1)

            read(200,*,iostat=ior) idp,iso,mass,x,y,z,px,py,pz,ActualEvent,subev
            GlobalID     = 0
            prodTime     = -999.
            lastCollTime = -999.
            collHis      = -999

         case(2)
            read(200,*,iostat=ior) idp
            backspace(Unit=200)

            if (idp < 1000) then
               read(200,*,iostat=ior) idp,iso,mass,x,y,z,px,py,pz,ActualEvent,subev
               GlobalID     = 0
               prodTime     = -999.
               lastCollTime = -999.
               collHis      = -999
            else
               read(200,*,iostat=ior) GlobalID,idp,prodTime,lastCollTime,collHis, & 
                    & iso,mass,x,y,z,px,py,pz,ActualEvent,subev
            end if
         end select

         if (ActualEvent /= j) then
            backspace(Unit=200) !move one line back!!!!
            exit Particles
         endif

         if (ior<0) exit Particles !end of file, exit in any case...

         if ( idp==2 ) then !not still decayed resonances!
            write(*,*) 
            write(*,*) 'smmModule/Get_ParticleVector:'
            write(*,*) 'Subsequent run = ',isu,'  NumEns = ',j
            write(*,*) 'warning: resonance as emitted particle: idp,iso =  ',idp,iso
         endif

         m  = m + 1
         p0 = sqrt(mass**2+px**2+py**2+pz**2)
         ParticleVector(j,m)%number      = globalID
         ParticleVector(j,m)%bornTime    = prodTime
         ParticleVector(j,m)%lastCollTime= lastCollTime
         ParticleVector(j,m)%collHis     = collHis
         ParticleVector(j,m)%id          = idp
         ParticleVector(j,m)%event       = subev
         ParticleVector(j,m)%ensemple    = ActualEvent
         ParticleVector(j,m)%Charge      = iso
         ParticleVector(j,m)%mass        = mass
         ParticleVector(j,m)%position(1) = x ![fm]
         ParticleVector(j,m)%position(2) = y ![fm]
         ParticleVector(j,m)%position(3) = z ![fm]
         ParticleVector(j,m)%momentum(0) = p0 ![GeV]
         ParticleVector(j,m)%momentum(1) = px ![GeV]
         ParticleVector(j,m)%momentum(2) = py ![GeV]
         ParticleVector(j,m)%momentum(3) = pz ![GeV]
      end do Particles

      if (m>mmax) mmax = m

   end do ENSEMPLES

   write(*,*) 'm_max = ',mmax

   if (isu==SubEvents) then
      close(unit=200,status='keep',iostat=ioe)
      if(ioe /= 0) then
         write(*,*) 'Get_ParticleVector: cannot close BUU_DataFile, IOSTAT = ',ioe
         write(*,*) 'Get_ParticleVector: !!! Termination of program NOW !!!'
         STOP
      endif
   endif

   end subroutine Get_ParticleVector !*************************************

  !*************************************************************************
  !****s* smmModule/HistOfSource
  ! NAME
  ! Subroutine HistOfSource
  !
  ! PURPOSE
  ! * Distributes events with respect to mass-/charge-number & excitation energy.
  !*************************************************************************
  subroutine HistOfSource(ActualSubEvent,ActualNumEnsemple,&
                          SubEvents,NumEnsemples,Masse,Ladung,Excitation)

    integer, intent(in) :: Masse,Ladung,ActualSubEvent,ActualNumEnsemple,& 
         &                 SubEvents,NumEnsemples
    real, intent(in) :: Excitation !GeV/A

    integer :: im,ibinm,ibinz,ibine

    integer, parameter :: mval=200
    real,    parameter :: dm=1.
    integer, parameter :: zval=100
    real,    parameter :: dz=1.
    integer, parameter :: eval=400
    real,    parameter :: de=0.05 !E_x - bin in MeV/A !!!

    real, dimension(1:mval),SAVE :: dndm
    real, dimension(1:zval),SAVE :: dndz
    real, dimension(1:eval),SAVE :: dnde,dndee,Exc_bin
    logical,                SAVE :: Hysto_Flag=.true.

    real :: E_x,Wert,Norm,MeanEx

    if (Hysto_Flag) then
       call Hysto_Init
       Hysto_Flag = .false.
    endif

    !Distributions in A
    ibinm = int( float(Masse)/dm)
    ibinm = min( mval, max(1,ibinm) )
    dndm(ibinm) = dndm(ibinm) + 1.

    !Distributions in Z
    ibinz = int( float(Ladung)/dz)
    ibinz = min( zval, max(1,ibinz) )
    dndz(ibinz) = dndz(ibinz) + 1.

    !Distributions in E
    E_x = Excitation*1000. !in MeV/A
    if (E_x > 0.) then
       ibine = int( E_x/de)
       ibine = min( eval, max(1,ibine) )
       dnde(ibine) = dnde(ibine) + 1.
       dndee(ibine) = dndee(ibine) + E_x
    endif

    if ( ActualSubEvent == SubEvents .and. & 
       & ActualNumEnsemple == NumEnsemples ) then

       do im=1,mval
          write(100,1000) float(im),dndm(im)/dm/float(SubEvents*NumEnsemples)
       end do

       do im=1,zval
!          Wert = Wert + dndz(im)/dz/float(SubEvents*NumEnsemples)
          write(101,1000) float(im),dndz(im)/dz/float(SubEvents*NumEnsemples)
       end do

       Wert = 0.0
       Norm = 0.0
       do im=1,eval
          Wert = Wert + dndee(im)
          Norm = Norm + dnde(im)
       end do
       MeanEx = Wert/Norm

       write(102,*) '# <E_x> (MeV/A) : ',MeanEx
       do im=1,eval
          write(102,1000) Exc_bin(im),dnde(im)/de/float(SubEvents*NumEnsemples)
       end do

    endif

1000 format(2ES15.4)

  contains

    subroutine Hysto_Init
      integer :: i

      dndm(:) = 0.0
      dndz(:) = 0.0
      do i=1,eval
         Exc_bin(i) = float(i)*de
         dndee(i) = 0.0
         dnde(i) = 0.0
      end do
      Exc_Bin(:) = Exc_bin(:) - de/2.

    end subroutine Hysto_Init


  end subroutine HistOfSource !*********************************************



end module smmModule
