program CalcJoos76

  use constants
  use inputGeneral
  use particleProperties, only: initParticleProperties
  use particleDefinition
  use Coll_gammaN
  use PythiaSpecFunc
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
  use minkowski, only: abs4

  use collisionTerm, only: collideMain
  use preEventDefinition
  use PreEvListDefinition
  use PreEvList, only: CreateSortedPreEvent,PreEvList_PrintEntry
  use Electron_origin, only : origin_singlePi,origin_doublePi,origin_DIS
  use insertion, only: GarbageCollection

  use PILCollected

  implicit none

  integer :: iW,iQ2,iN,iNdone,iN2,iHH,iiW

  integer :: iWmin,iWmax,iWdelta
  integer :: iQ2min,iQ2max,iQ2delta
  integer :: iNmax

  real    :: W, Q2
  real :: nu,Wfree,eps,fT,Ebeam
  Type(electronNucleon_event), save :: eNev_InitData
  Type(electronNucleon_event), save :: eNev
  logical :: doQE,doRes,do1Pi,do2Pi,doDIS

  logical,save,dimension(2:nres+1) :: useRes=.true.
  integer :: channel
  logical :: flagOK
  real :: s,x

  real :: XS
  real :: hhh
  
  type(particle) :: TargetNuc
  type(particle), dimension(1,1) :: realPart
  type(particle),dimension(1,1:25)  :: finalState

  type(preEvent), dimension(1:25) :: PreE
  type(tPreEvListEntry), dimension(1:6) :: Pre1Pi
  type(tPreEvListEntry), dimension(1:12) :: Pre2Pi
  integer :: iC1,iC2,iC3, iiC, i
  real :: SumALL(0:5), Sum1pi(0:6,0:5),Sum2pi(0:12,0:5)

  NAMELIST /datatable/ iWmin,iWmax,iWdelta,iQ2min,iQ2max,iQ2delta,iNmax,&
       & doRes,do1Pi,do2Pi,doDIS,Ebeam


  call readinputGeneral
  call initParticleProperties
  call forceInitFormation
!  call InitParticleProperties

  call SetSomeDefaults_PY

  call setToDefault(TargetNuc)
  TargetNuc%ID = 1
  TargetNuc%charge = 1
  TargetNuc%mass = 0.938
  TargetNuc%momentum = (/0.938, 0.0, 0.0, 0.0 /)
  TargetNuc%Position = 9999999d0

  Ebeam = 7.2 ! Joos 1976

  iWmin   = 1100
  iWmax   = 1990
  iWdelta =   20

  iQ2min   =    0
  iQ2max   = 2000
  iQ2delta = 1000

  iNmax = 10000

  doQE  = .false.
  doRes = .true.
  do1Pi = .true.
  do2Pi = .true.
  doDIS = .true.

  call Write_ReadingInput('datatable',0)
  rewind(5)
  read(5,nml=datatable)
  write(*,*) ' Ebeam= ',Ebeam
  write(*,'(A,3i7)') '  iW   = ',iWmin,iWmax,iWdelta
  write(*,'(A,3i7)') '  iQ2  = ',iQ2min,iQ2max,iQ2delta
  write(*,'(A,3i7)') '  iN   = ',1,iNmax,1
  write(*,*)
  write(*,*) '  doQE   =',doQE
  write(*,*) '  doRes  =',doRes
  write(*,*) '  do1Pi  =',do1Pi
  write(*,*) '  do2Pi  =',do2Pi
  write(*,*) '  doDIS  =',doDIS

  call Write_ReadingInput('datatable',1)
  write(*,*)

  realPart(1,1)%ID = 0

  !...Prepare the PreEvents to compare with:

  write(  *,*) '1--Pion FinalStates:'
  write(100,*) '1--Pion FinalStates:'

  finalstate%ID = 0
  finalstate%antiparticle = .false.
  finalstate(1,1:2)%ID=(/1,101/)

  iiC=0
  do iC1=1,0,-1
     do iC2=1,-1,-1
        finalstate(1,1:2)%charge=(/iC1,iC2/)
        iiC=iiC+1
        allocate(Pre1Pi(iiC)%preE(6))
        if(.not.CreateSortedPreEvent(finalState(1,:),pre1Pi(iiC)%preE)) then
           write(*,*)  'Error 1:',iiC
           call WriteParticle(6,1,finalState(1,:))
           stop
        endif
        call PreEvList_PrintEntry(  6,pre1Pi(iiC),1.0)
        call PreEvList_PrintEntry(100,pre1Pi(iiC),1.0)
!        call WriteParticle(6,iiC,finalstate(1,:))
     end do
  end do

  write(  *,*) '2--Pion FinalStates:'
  write(100,*) '2--Pion FinalStates:'

  finalstate%ID = 0
  finalstate%antiparticle = .false.
  finalstate(1,1:3)%ID=(/1,101,101/)

  iiC=0
  do iC1=1,0,-1
     do iC2=1,-1,-1
        do iC3=iC2,-1,-1
           finalstate(1,1:3)%charge=(/iC1,iC2,iC3/)
           iiC=iiC+1
           allocate(Pre2Pi(iiC)%preE(6))
           if(.not.CreateSortedPreEvent(finalState(1,:),pre2Pi(iiC)%preE)) then
              write(*,*)  'Error 2:',iiC
              call WriteParticle(6,1,finalState(1,:))
              stop
           endif
           call PreEvList_PrintEntry(  6,pre2Pi(iiC),1.0)
           call PreEvList_PrintEntry(100,pre2Pi(iiC),1.0)
        end do
     end do
  end do

  do iQ2=iQ2min,iQ2max,iQ2delta
     Q2 = iQ2*0.001
!     if (Q2<0.01) Q2=0.1
     if (Q2<0.01) Q2=0.01

     iiW = 0
     do iW=iWmin,iWmax,iWdelta
        W = iW*0.001
        iiW = iiW+1

        eNev_InitData = eNev_init_BWQ(Ebeam,W,Q2)
        eNev = eNev_init_Target(eNev_InitData,TargetNuc,flagOK)

        call write_electronNucleon_event(eNev,.FALSE.,.FALSE.) 
        call eNeV_GetKinV(eNev, nu,Q2,W,Wfree,eps,fT)
        s=abs4(eNev%nucleon%momentum+eNev%electron_in)
        x=eNeV_Get_LightX(eNev)
        write(*,*) 'nu :',nu
        write(*,*) 'xB :',Q2/(2*0.938*nu)
        write(*,*) 'sqrt(s) :',s
        write(*,*) 'x  :',x
        write(*,*) 'eps:',eps

!        stop

        call PILCollected_ZERO

        iNdone=0
        SumALL = 0.0
        Sum1pi = 0.0
        Sum2pi = 0.0



        hhh = 1./(137.*4*3.14*(Ebeam*0.938)**2) * W*(W**2-0.938**2)/((1-eps)*Q2)

        write(  *,'(2f7.3,1P,12e13.4)') W,Q2, nu,Q2/(2*0.938*nu),eps,fT,hhh
        write(111,'(2f7.3,1P,12e13.4)') W,Q2, nu,Q2/(2*0.938*nu),eps,fT,hhh

        do iN=1,iNmax
           if (DoPR(2)) write(*,*) '=======iN =',iN

           if (W.lt.2.0) then
              call eventGen_eN_lowEnergy(eNev,doQE,doRes,useRes,do1Pi,do2Pi,doDIS,finalState(1,:),channel,flagOK,XS)
           else
              call eventGen_eN_HiEnergy(eNev,iN,(/1.,1.,1.,1./),finalState(1,:),channel,flagOK,XS)
           endif

           if (.not.flagOK) cycle

           if (DoPR(2)) write(*,*) 'channel:',channel
!           call WriteParticle(6,1,finalState(1,:))

           call collideMain(finalstate,realPart, 0.0, .true.)
!           call WriteParticle(6,1,finalState(1,:))

           call GarbageCollection(finalstate)

           if (.not.CreateSortedPreEvent(finalState(1,:),preE)) then
              write(*,*) 'Error PreE'
              call WriteParticle(6,1,finalState(1,:))
              cycle
           end if

           if (W.lt.2.0) then
              iHH=2
              select case(channel)
              case (origin_singlePi)
                 iHH=3
              case (origin_doublePi)
                 iHH=4
              case (origin_DIS)
                 iHH=5
              end select
           else
              iHH=2
              select case(channel)
              case (5000+origin_singlePi)
                 iHH=3
              case (5000+origin_doublePi)
                 iHH=4
              case (2004,5000+origin_DIS)
                 iHH=5
              end select
           endif

           !... Total Cross Section:

           iNdone = iNdone+1
           SumALL(  0) = SumALL(  0) + XS
           SumALL(iHH) = SumALL(iHH) + XS

           !... 1Pion final state:
           do iiC=1,6
              flagOK=.true.
              do i=1,5
                 if (preE(i)%mass.ne.pre1Pi(iiC)%preE(i)%mass) flagOK=.false. ! attention: abuse of mass
              end do
              if(flagOK) then
                 Sum1pi(0,  0)  = Sum1pi(0,  0)  +XS
                 Sum1pi(0,  iHH)= Sum1pi(0,  iHH)+XS
                 Sum1pi(iiC,0)  = Sum1pi(iiC,0)  +XS
                 Sum1pi(iiC,iHH)= Sum1pi(iiC,iHH)+XS
                 exit
              end if
           end do

           !... 2Pion final state:
           do iiC=1,12
              flagOK=.true.
              do i=1,5
                 if (preE(i)%mass.ne.pre2Pi(iiC)%preE(i)%mass) flagOK=.false. ! attention: abuse of mass
              end do
              if(flagOK) then
                 Sum2pi(0,  0)  = Sum2pi(0,  0)  +XS
                 Sum2pi(0,  iHH)= Sum2pi(0,  iHH)+XS
                 Sum2pi(iiC,0)  = Sum2pi(iiC,0)  +XS
                 Sum2pi(iiC,iHH)= Sum2pi(iiC,iHH)+XS
                 exit
              end if
           end do

        end do

        call eNeV_GetKinV(eNev, nu,Q2,W,Wfree,eps,fT)
        if (W.lt.2.0) then
           fT = fT/ ( 1e3* pi/(eNev%electron_out(0)*eNev%electron_in(0)))
        else
           fT = 1.0
        end if

        SumALL = SumALL/(iNmax*fT)
        Sum1pi = Sum1pi/(iNmax*fT)
        Sum2pi = Sum2pi/(iNmax*fT)
 

        write(121,'(2f7.3,i9,1P,12e13.4)') W,Q2, iNdone,SumALL

        do iHH=0,5
           write( 130+iHH,'(2f7.3,i9,1P,20e13.4)') W,Q2, -99,Sum1pi(:,iHH)
           write(1300+iHH,'(2f7.3,i9,1P,20e13.4)') W,Q2, -99,Sum2pi(:,iHH)
        end do

     end do

     if (iiW>1) then ! write empty lines for gnuplot (cf. index)
        write(111,*)
        write(111,*)
        write(121,*)
        write(121,*)
        do iHH=0,5
           write( 130+iHH,*)
           write( 130+iHH,*)
           write(1300+iHH,*)
           write(1300+iHH,*)
        end do
     end if

  end do
  write(*,*) 'Done.'


end program CalcJoos76
