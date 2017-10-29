!***************************************************************************
!****m* /fragmentation
! NAME
! module fragmemtation

! FUNCTION
! This module prepares the input for the 
! Statistical Multifragmentation Model (SMM), runs this code for all 
! MC-events and constructs the FragmentVector.
! NOTES
! * Please take care of the parameter "EnsemplesMax": It is important for 
!   intermediate variables between Fortran90 and Fortran77 routines. It fixes 
!   the dimension of variables exchanged between Fortran90 and Fortran77 
!   routines, since in old Fortran77 allocatable arrays cannot be declared!!!
!   The value of the parameter "EnsemplesMax" must coincide with the 
!   corresponding parameter in the MAINMA subroutine.
! * SMM-code uses M_nuc=0.940 GeV as the bare nucleon mass !!!
!   take care on this when calculating observables related to kin. energies 
!   of SMM fragments.
!***************************************************************************

module fragmentation

  PRIVATE

  PUBLIC :: MultiFragmentation

contains 

  subroutine MultiFragmentation(SMM_Seed,SMM_Events, & 
       & SubRuns,NumEns,NumSources, & 
       & EventType,SMM_Flag,ALADIN_Flag,TheSource,FragmentVector)

    use typeDefinitions,     only : cluster,quelle
    use VelocityFieldModule, only : Get_VelocityFieldAt

    implicit none
    !---------------------------------------------------------------------
    !Input-Variables
    !---------------------------------------------------------------------
    integer,                       intent(in)    :: NumEns
    integer,                       intent(in)    :: SubRuns
    integer,                       intent(in)    :: SMM_Events
    integer,                       intent(in)    :: SMM_Seed
    integer,                       intent(in)    :: EventType
    integer,                       intent(in)    :: SMM_Flag
    integer,       dimension(:),   intent(in)    :: NumSources
    type(quelle),  dimension(:,:), intent(in)    :: TheSource
    logical,                       intent(in)    :: ALADIN_Flag
    type(cluster), dimension(:,:), intent(inout) :: FragmentVector
    !---------------------------------------------------------------------
    !Local Variables
    !---------------------------------------------------------------------
    real,    dimension(1:3) :: Pos,beta,collectiveFlow
    real,    dimension(0:3) :: pin,pinLRF

    integer                       :: prodParticles,totalStates
    integer, dimension(1:500)     :: Afinal,Zfinal,WhichProcess
    real,    dimension(0:3,1:500) :: Momentum
    real,    dimension(1:3,1:500) :: Ort

    integer :: SourceSize,SourceCharge,iflow
    real    :: Ex,Erad,Eflow
    real    :: SourceMass, gamma
    integer :: j,m,l,ii,n

    integer, SAVE :: LocalSeed

    logical :: ExFlag

    !---------------------------------------------------------------------
    ! Local initializations
    !---------------------------------------------------------------------
    if (NumEns==1) then
       LocalSeed = SMM_Seed !local Seed for SMM-RG
    endif
    ii = 0 !counting total SMM runs (=SMM_Events)

    !---------------------------------------------------------------------
    ! call SMM code event-by-event (=SubEvents*NumEnsemples = GiBUU-events)
    ! * For each GiBUU-run SMM_Events are performed
    ! * SMM USES ITS OWN RANDOM-NUMBER GENERATOR !!!
    !---------------------------------------------------------------------
    
    Loop_over_SMMruns : do n=1,SMM_Events

       LocalSeed = LocalSeed + 1

       ii = ii + 1

       totalStates = 0

       Loop_over_Sources : Do j=1,NumSources(NumEns)

          if ( .not. TheSource(NumEns,j)%status ) cycle Loop_over_Sources

          if (ALADIN_Flag) then
!             if (TheSource(NumEns,j)%Type /= 2) CYCLE Loop_over_Sources
             if (TheSource(NumEns,j)%Type == 3) CYCLE Loop_over_Sources
          endif

          SourceSize   = TheSource(NumEns,j)%Size
          SourceCharge = TheSource(NumEns,j)%Charge
          Erad         = TheSource(NumEns,j)%radEnergy
          Ex           = TheSource(NumEns,j)%ExEnergy

          if (EventType .ne. 2) then
             beta(1)      = TheSource(NumEns,j)%velocity(1)
             beta(2)      = TheSource(NumEns,j)%velocity(2)
             beta(3)      = TheSource(NumEns,j)%velocity(3)
          else
             beta(:) = 0.0 ! for p+X E_exc is calculated in target rest frame!!!
                           ! Do not boost 4-momenta of SMM fragments for p+X !!!
          endif

          !          if (j > 2) then
          if (Erad > 0.010) then
             iflow = 1
             Eflow = Erad
             if (ALADIN_Flag) then
                write(*,*) 'module fragmentation, routine MultiFragmentation:'
                write(*,*) ' * Wrong input-type for fragmenting source'
                write(*,*) '   when analyzing ALADIN/INDRA data!!!'
                write(*,*) 'This option has to be used ONLY for spectator fragmentation'
                write(*,*) 'Termination of the program'
                STOP
             endif
          else
             iflow = 0
             Eflow = 0.0
          endif

          if (EX < 0.00001) then
             EX = 0.0
          endif

          !calling of SMM code:
          if (EX > 0.0) then
             call MAINMA(j,LocalSeed,SourceSize,SourceCharge,EX,& 
                  &      SubRuns,Eflow,iflow,SMM_Flag,& 
                  &      prodParticles,Afinal,Zfinal, & 
                  &      WhichProcess, & 
                  &      Ort,momentum)
             ExFlag = .true.
          else
             ProdParticles = 1
             ExFlag = .false.
          endif

          if (EX .le. 0.00001 .and. prodParticles /= 1) then
             write(*,*) 'module fragmentation, routine MultiFragmentation:'
             write(*,*) ' * Wrong output for SMM-clusters!!!'
             write(*,*) 'Termination of the program'
             STOP             
          endif

          call checks(prodParticles)

          !---------------------------------------------------------------------
          ! * Construct FragmentVector (includes also free nucleons)
          ! * boost momenta of the final states into the global CMS
          ! * transforme also the positions into the CF
          !---------------------------------------------------------------------
          Loop_over_SMMparticles : do m=1,prodParticles

             l = totalStates + m

             if (EX > 0.0) then
                call checks(l)
                pin(0:3) = Momentum(0:3,m) !p^mu in moving frame
                call checks(1,pin(0),'vor boost')
             endif

             !SMM distributes particles considering no radial 
             !dependence of radial expansion! Thus, one has 
             !to distribute the produced in coordinate space, and 
             !then boost them acc. their radial flow beta(r).
             !Distribution in coordinate space still under study...
             if (EventType==1 .and. iflow==1) then
                if (ALADIN_Flag) then
                   write(*,*) 'module fragmentation, routine MultiFragmentation:'
                   write(*,*) ' * Wrong input-type for fragmenting source'
                   write(*,*) '   when analyzing ALADIN/INDRA data!!!'
                   write(*,*) 'This option has to be used ONLY for spectator fragmentation'
                   write(*,*) 'Termination of the program'
                   STOP
                endif
                Pos(1:3) = Ort(:,m) + TheSource(NumEns,j)%position(:)
                call Get_VelocityFieldAt(Pos,collectiveFlow)
                beta(1:3) = collectiveFlow(1:3) !radial flow profile
                if ((1.-Dot_Product(beta,beta)) .le.0.0000001) then
                   write(*,*) 'fragmentation module, MultiFragmentation:'
                   write(*,*) 'wrong determination of radial flow!!!'
                   write(*,'(A,6f12.5)') 'Pos,beta = ',&
                        & Pos(:),beta(:)
                   STOP
                endif
             endif

             if (EX > 0.0) then
                pinLRF(:) = pin(:)
                if (EventType .ne. 2) then 
                   call lorentz(-beta,pin,'MultiFragmentation') !boost back from LRF to global CF
                   call checks(1,pin(0),'nach boost')
                endif

                FragmentVector(ii,l)%position(:)  = Ort(:,m) + TheSource(NumEns,j)%position(:)
                FragmentVector(ii,l)%momentumLRF(:) = pinLRF(:) !p^mu in moving frame
                FragmentVector(ii,l)%momentum(:)  = pin(:)   !p^mu in global CMS-frame
                FragmentVector(ii,l)%Mass         = Afinal(m)*0.940 !see notes !!!
                FragmentVector(ii,l)%MassNumber   = Afinal(m)
                FragmentVector(ii,l)%ChargeNumber = Zfinal(m)
                FragmentVector(ii,l)%HypNumber    = 0
                FragmentVector(ii,l)%ID           = 1
                FragmentVector(ii,l)%Mechanism    = WhichProcess(m)
                FragmentVector(ii,l)%stableFlag   = .true.
                if ( (Afinal(m)==1) .and. (Zfinal(m)==1.or.Zfinal(m)==0)) then
                   FragmentVector(ii,l)%FreeBound    = .false.
                else
                   FragmentVector(ii,l)%FreeBound    = .true.
                endif

             else

                SourceMass = float(SourceSize)*0.940
                gamma      = 1./sqrt(1.-dot_product(beta(1:3),beta(1:3)))
                pin(1:3)   = beta(1:3) *gamma * SourceMass
                pin(0)     = sqrt(SourceMass**2+dot_product(pin(1:3),pin(1:3)))

                FragmentVector(ii,l)%position(:)  = TheSource(NumEns,j)%position(:)
                FragmentVector(ii,l)%momentum(:)  = pin(:)   !p^mu in global CMS-frame
                FragmentVector(ii,l)%Mass         = SourceMass
                FragmentVector(ii,l)%MassNumber   = SourceSize
                FragmentVector(ii,l)%ChargeNumber = SourceCharge
                FragmentVector(ii,l)%HypNumber    = 0
                FragmentVector(ii,l)%ID           = 1
                FragmentVector(ii,l)%Mechanism    = 0
                FragmentVector(ii,l)%stableFlag   = .false.
                if ( (SourceSize==1) .and. (SourceCharge==1.or.SourceCharge==0)) then
                   FragmentVector(ii,l)%FreeBound    = .false.
                else
                   FragmentVector(ii,l)%FreeBound    = .true.
                endif

             end if

          end do Loop_over_SMMparticles

          totalStates = totalStates + prodParticles

       end do Loop_over_Sources

    end do Loop_over_SMMruns

  !***********************************************************************
  end subroutine MultiFragmentation !*************************************
  !***********************************************************************

  !***********************************************************************
  !****s* fragmentation/checks
  ! NAME
  ! subroutine checks
  !***********************************************************************
  subroutine checks(value1,value2,Name)
    implicit none
    integer,               intent(in) :: value1
    real,        optional, intent(in) :: value2
    character(*),optional, intent(in)  :: Name

    if (value1 > 500) then
       write(*,*) 'module smm_Main, routine checkDimensions:'
       write(*,*) 'dimension overflow in fragment vector! ',value1
       write(*,*) 'Termination of program'
       STOP
    endif

    if (present(value2)) then
       if (value2 < 0.940) then
          write(*,*) 'module smm_Main, routine MultiFragmentation:'
          if (present(Name)) write(*,*) 'energy < M ????',Name
          write(*,*) 'Termination of program'
          STOP
       endif
    endif

  !***********************************************************************
  end subroutine checks !*************************************************
  !***********************************************************************

  !Loretz trafo routine, copied from GiBUU

  !***********************************************************************
  !****s* fragmentation/lorentz
  ! NAME
  ! subroutine lorentz(beta,fourVector,CallName)
  ! FUNCTION
  ! performs Lorentz transformation of fourVector into a 
  ! system which is traveling with the velocity beta(1:3) 
  ! INPUTS
  !  real,dimension(0:3),intent(inout) :: fourVector
  !  real,dimension(1:3),intent(in)    :: beta
  !  character(*),intent(in),optional  :: CallName
  ! RESULT
  !  real,dimension(0:3),intent(inout) :: fourVector
  !***********************************************************************
  subroutine lorentz(beta, fourVector, CallName)

    implicit none
    real,dimension(0:3),intent(inout) :: fourVector
    real,dimension(1:3),intent(in)    :: beta
    character(*),optional,intent(in)  :: CallName

    real  :: gamma
    real  :: betaFour

    !Evaluate gamma
    gamma = 1.0 - Dot_Product(beta,beta)
    if(gamma .gt. 0.0) then
       gamma = 1.0/sqrt(gamma)
    else
       write(*,*)
       if (Present(CallName)) then
          write(*,*) 'lorentz called by ',CallName
       endif
       write(*,*)'(1-beta**2) in lorentz less or equal zero: ', gamma,beta
       write(*,*)'beta=',  beta
       write(*,*)'Stop program'
       stop
    end if
    !beta*fourVector(1:3)
    betaFour = Dot_Product(beta, fourVector(1:3))
    !do transformation
    fourVector(1:3)= fourVector(1:3) + & 
         & gamma*beta(1:3)*(gamma/(gamma+1.)*betaFour-fourvector(0))
    fourVector(0)= gamma*(fourVector(0)-betaFour)
    return
  !***********************************************************************
  end subroutine lorentz
  !***********************************************************************

end module fragmentation
