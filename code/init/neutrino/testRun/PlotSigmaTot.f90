program PlotSigmaTot

  use constants
  use inputGeneral
  use particleDefinition
  use Coll_gammaN
  use VMMassPythia
  use hadronFormation, only : forceInitFormation
  use Coll_Pythia
  use CollTools

  use eN_eventDefinition
  use eN_event
  use eventGenerator_eN_lowEnergy
  use eventGenerator_eN_HiEnergy
  use idTable,only :nres

  use ParamEP
  use output
  use minkowski, only: abs4,abs4Sq


  use collisionTerm, only: collideMain
  use preEventDefinition
  use PreEvListDefinition
  use PreEvList, only: CreateSortedPreEvent,PreEvList_PrintEntry
  use lowElectron_origin, only : origin_singlePi,origin_doublePi,origin_DIS
  use insertion, only: GarbageCollection

  use Coll_nuN

  use PILCollected

  use AnaEventDefinition
  use AnaEvent, only : event_init,event_clear,event_add
  use MultiplicityAnalysis

  implicit none

  Type(electronNucleon_event), save :: eNev0
  Type(electronNucleon_event), save :: eNev
  type(particle) :: TargetNuc
  type(particle), dimension(1,1) :: realPart
  type(particle),dimension(1,1:50),TARGET  :: finalState

  type(tAnaEvent)            :: tEv

  real :: Ebeam,cross
  logical :: flagOK
  integer :: iEV

  ! Parameters to play with:

!  integer :: nEV=1000
  integer:: nEV=10
  logical :: doCC = .true.
!  logical :: doCC = .false.
!  logical :: doMassless=.true.
  logical :: doMassless=.false.
  integer :: iGeneration = 2
  integer:: cTarget = 1
  real :: Wcut = 1.6

  logical :: doMultAna = .true.

  NAMELIST /datatable/ nEV,doCC,doMassless,iGeneration,cTarget,Wcut


  integer :: iEn
  real :: SumAna1,SumAna2

  integer :: iPi,nPi
  integer :: iC1,iC2, iiC, i

  type(preEvent), dimension(1:50) :: PreE
  
  type(tPreEvListEntry), dimension(1:6) :: PreELE
  real :: SumAll,Sum1pi(1:6)

  character*(10) :: Prefix_MultAna
  type(particle), POINTER :: pPart

  call readinputGeneral
  call init_database
  call forceInitFormation
!  call InitParticleProperties

  call SetSomeDefaults_PY

  !...Reading input:

  call Write_ReadingInput('datatable',0)
  rewind(5)
  read(5,nml=datatable)
  write(*,'(A,3i7)') '  nEv   = ',nEv
  write(*,'(A,3i7)') '  iGen  = ',iGeneration
  write(*,*)
  write(*,'(A,3i7)') '  chargeTarget  = ',cTarget
  write(*,*)
  write(*,*) '  doCC        =',doCC
  write(*,*) '  doMassless  =',doMassless

  call Write_ReadingInput('datatable',1)
  write(*,*)

  !...Set up the target nucleon:

  call setToDefault(TargetNuc)
  TargetNuc%ID = 1
  TargetNuc%charge = cTarget
  TargetNuc%mass = 0.938
  TargetNuc%momentum = (/0.938, 0.0, 0.0, 0.0 /)
  TargetNuc%Position = 9999999d0

  realPart(1,1)%ID = 0

  !...Prepare the PreEvents to compare with:

  finalstate%ID = 0
  finalstate%antiparticle = .false.
  finalstate(1,(/1:2/))%ID=(/1,101/)

  write(*,*) 'FinalStates:'

  iiC=0
  do iC1=1,0,-1
     do iC2=1,-1,-1
        finalstate(1,(/1:2/))%charge=(/iC1,iC2/)
        iiC=iiC+1
        allocate(PreELE(iiC)%preE(6))
        if(.not.CreateSortedPreEvent(finalState(1,:),preELE(iiC)%preE)) then
           write(*,*)  'Error 1:',iiC
           call WriteParticle(6,1,finalState(1,:))
           stop
        endif
        call PreEvList_PrintEntry(6,preELE(iiC),1.0)
!        call WriteParticle(6,iiC,finalstate(1,:))
     end do
  end do

  !...Initialize the tEvent type

  call event_INIT(tEv)

  !...Set up the initial kinematics:

  do iEn=50,-4,-1
     Ebeam = 10d0**(iEn*0.05d0)
!  do iEn=20,3,-1
!  do iEn=100,3,-1
!     Ebeam = 0.1d0*iEn

     call PILCollected_ZERO

     eNev0%electron_in = (/Ebeam,0d0, 0d0,Ebeam/)
     eNev0%electron_out= eNev0%electron_in ! !!!! DUMMY !!!!
     eNev0%photon      = eNev0%electron_in-eNev0%electron_out
     
     eNev0%nucleon = TargetNuc
     
     eNev0%QSquared = -abs4Sq(eNev0%photon)
     eNev0%W = abs4(eNev0%photon+eNev0%nucleon%momentum)
     
     eNev0%nucleon_free      = eNev0%nucleon
     
     eNev0%W_free = abs4(eNev0%photon+eNev0%nucleon_free%momentum)
     
     call write_electronNucleon_event(eNev0,.FALSE.,.TRUE.)
     
     if (doMultAna) then
        write(Prefix_MultAna,'(f8.4)') Ebeam
        call Multiplicity_Reset
     endif

     !...Do some Events:
     
     SumAll = 0d0
     Sum1pi = 0d0

     do iEV=1,nEV
        
        eNev = eNev0
        call eNev_init_Target(eNev,TargetNuc,flagOK)
        if (.not.flagOK) cycle
        
        call DoColl_nuN_Py(eNev,finalState(1,:),flagOK, doCC, iGeneration, doMassless, cross)
        if (.not.flagOK) cycle

        if (eNev%W_free.lt.Wcut) cycle

!        call write_electronNucleon_event(eNev,.FALSE.,.FALSE.)
!        call WriteParticle(6,1,finalState(1,:))
!        write(*,*) 'cross section: (in mb)',cross
        SumAll=SumAll+cross

        !...Force all particles to decay...
        call collideMain(finalstate,realPart, 0.0, .true.)
        call GarbageCollection(finalstate)

        finalstate(1,:)%perweight = 1.0

        !...Do Multiplicity analysis (if desired)
        if (doMultAna) then
           call event_CLEAR(tEv)
           do i=1,size(finalstate,dim=2)
              if (finalstate(1,i)%ID.lt.0) exit
              if (finalstate(1,i)%ID.eq.0) cycle
              pPart => finalstate(1,i)
              call event_ADD(tEv,pPart)
           end do
           call Multiplicity_AddEvent(tEv)
        end if


        !...Compare with stored preEvents:
        if (.not.CreateSortedPreEvent(finalState(1,:),preE)) then
           write(*,*) 'Error PreE'
           call WriteParticle(6,1,finalState(1,:))
           cycle
        end if
        do iiC=1,6
           flagOK=.true.
           do i=1,5
              if (preE(i)%mass.ne.preELE(iiC)%preE(i)%mass) flagOK=.false. ! attention: abuse of mass
           end do
           if(flagOK) then
              Sum1pi(iiC)= Sum1pi(iiC)+cross
              exit
           end if
        end do

     end do

     call Multiplicity_Write(Prefix_MultAna)

     SumAll=SumAll/nEV
     Sum1pi=Sum1pi/nEV

     write(121,'(f12.5,1P,20e13.4)') Ebeam,SumAll,Sum(Sum1pi),Sum1pi
     write(*,'(f12.5,1P,20e13.4)') Ebeam,SumAll,Sum(Sum1pi),Sum1pi

     call DoAnaEstimate

     write(122,*) Ebeam,SumAna1,SumAna2
     write(*,*) 'Ana:',Ebeam,SumAna1,SumAna2

  end do

  write(*,*) 'done.'

contains

    subroutine DoAnaEstimate

    use CollTools
    IMPLICIT NONE

    integer :: iX, iY
    real :: x,y, Q2, nu, Eprime,cost, sigma0,sigma1,sigma2, srts2,W
    real, dimension(-25:25) :: xq


    real, parameter :: GF=1.166E-5 ! in GeV-2
    real, parameter :: hc2 = 0.389 ! in GeV2 mb


    SumAna1 = 0d0
    sumAna2 = 0d0

    srts2 = 0.938**2+2*0.938*Ebeam

    sigma0 = GF**2*srts2*hc2/(2*3.14)  ! in mb !!!  

!    write(123,*) Ebeam,sqrt(srts2),sigma0

    call SetSomeDefaults_PY

    do iX=2,99,2
       x = iX*0.01
       do iY=2,99,2
          y = iY*0.01
          
          nu = y*Ebeam
          Q2 = 2*0.938*nu * x
          Eprime = Ebeam - nu
!          cost = 1.0 -2*y ! Halzen Martin, wrong !!!!!
          cost = 1- Q2/(2*Ebeam*Eprime)

          if (cost.lt.-1) then
!             write(122,'(2F13.5,4E13.5,5e13.5)') x,y,Q2,nu,Eprime,cost,&
!               & 1e-20,1e-20,1e-20,1e-20,1e-20
             cycle
          end if

          W = 0.938**2-Q2+2*0.938*nu
          if (W.lt.Wcut) cycle

!          call pypdfu(2212, x,Q2,xq)
          call pypdfl(2212, x,Q2,xq)

          if (TargetNuc%charge .gt. 0) then ! *** proton ***
             sigma1 = 2*sigma0 * (xq(1)+ (1-y)**2*xq(-2)) ! x(d + (1-y)**2ubar) 
             sigma2 = 2*sigma0 * (xq(-1)+ (1-y)**2*xq(2)) ! x(dbar + (1-y)**2u)
          else                                 ! *** neutron ***
             sigma1 = 2*sigma0 * (xq(2)+ (1-y)**2*xq(-1)) ! x(u + (1-y)**2dbar) 
             sigma2 = 2*sigma0 * (xq(-2)+ (1-y)**2*xq(1)) ! x(ubar + (1-y)**2d)
          endif

          SumAna1=SumAna1+sigma1
          SumAna2=SumAna2+sigma2

!          write(122,'(2F13.5,4E13.5,5e13.5)') x,y,Q2,nu,Eprime,cost,&
!               & sigma0, sigma1,sigma2, 2*0.938*nu*Ebeam, 0.938*nu/(Ebeam-nu)

       end do
!       write(122,*)

    end do

    SumAna1=SumAna1*(2*0.01)*(2*0.01)
    SumAna2=SumAna2*(2*0.01)*(2*0.01)

  end subroutine DoAnaEstimate

end program PlotSigmaTot
