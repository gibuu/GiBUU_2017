! This program calculates the total cross section as function of W and Q2.
! It also does a splitting into different channels:
! * W<2GeV: ...
! * W>2GeV: ...
! The program does a Rosenbluth separation, i.e. it calculates the cross
! section for two distinct values of epsilon and prints it for the 
! extrapolation eps=0 (==sigma_T) and eps=1 (==sigma_T+sigma_L).
!
! We added the printout of multiplicities in order to check the KNO scaling.
! (This is only done for the calculations with the smaller epsilon.)
!
! This program is a combination of "PlotLowXS.f90" and "PlotHighXS.f90",
! which became now obsolete. The treatment at the border W=2GeV has to
! be refined.

program plotAllXS

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

  use AnaEventDefinition
  use AnaEvent, only : event_init,event_clear,event_add
  use MultiplicityAnalysis

  use collisionTerm, only: collideMain
  use insertion, only: GarbageCollection

  implicit none

  integer :: iW,iQ2,iN,iN1,iN2, i

  integer :: iWmin,iWmax,iWdelta
  integer :: iQ2min,iQ2max,iQ2delta
  integer :: iNmax
  logical :: verbose

  real    :: W, Q2
  real :: nu,Wfree,eps,fT
  Type(electronNucleon_event), save :: eNev1,eNev2
  logical :: doQE,doRes,do1Pi,do2Pi,doDIS
  type(particle), dimension(1,1) :: realPart
  type(particle),dimension(1,1:50),TARGET  :: finalState
  logical,save,dimension(2:nres+1) :: useRes=.true.
  integer :: channel
  logical :: flagOK
  real :: s1,s2,x1,x2, R1,R2

  real :: XS, XS_sum1, XS_sum2
  real,dimension(0:4) :: XS_Arr, XS_Arr_sum1, XS_Arr_sum2
  ! please note: for W<2 we have dimension(1:5) 

  real, parameter :: eps1=0.05, eps2=0.99
  
  type(particle) :: TargetNuc
  integer :: chargeTarget

  logical :: isHigh

  logical :: doMultAna = .true.
  character*(20) :: Prefix_MultAna
  type(particle), POINTER :: pPart
  type(tAnaEvent)            :: tEv

  NAMELIST /datatable/ iWmin,iWmax,iWdelta,iQ2min,iQ2max,iQ2delta,iNmax,&
       & doRes,do1Pi,do2Pi,doDIS, verbose, chargeTarget


  call readinputGeneral
  call initParticleProperties
  call forceInitFormation
!  call InitParticleProperties

  call SetSomeDefaults_PY
  call event_INIT(tEv)

  iWmin   = 110
  iWmax   = 199
  iWdelta =   2

  iQ2min   =   0
  iQ2max   = 200
  iQ2delta = 100

  iNmax = 10000

  doQE  = .false.
  doRes = .true.
  do1Pi = .true.
  do2Pi = .true.
  doDIS = .true.

  verbose = .false.

  chargeTarget = 1

  call Write_ReadingInput('datatable',0)
  rewind(5)
  read(5,nml=datatable)
  write(*,'(A,3i7)') '  iW   = ',iWmin,iWmax,iWdelta
  write(*,'(A,3i7)') '  iQ2  = ',iQ2min,iQ2max,iQ2delta
  write(*,'(A,3i7)') '  iN   = ',1,iNmax,1
  write(*,*)
  write(*,*) '  doQE   =',doQE
  write(*,*) '  doRes  =',doRes
  write(*,*) '  do1Pi  =',do1Pi
  write(*,*) '  do2Pi  =',do2Pi
  write(*,*) '  doDIS  =',doDIS
  write(*,*)
  write(*,*) '  verbose=',verbose
  write(*,*)
  write(*,*) '  charge of target: ',chargeTarget

  call Write_ReadingInput('datatable',1)
  write(*,*)
  write(*,*) 'eps:',eps1,eps2

  call setToDefault(TargetNuc)
  TargetNuc%ID = 1
  TargetNuc%charge = chargeTarget
  TargetNuc%mass = 0.938
  TargetNuc%momentum = (/0.938, 0.0, 0.0, 0.0 /)
  TargetNuc%Position = 9999999d0

  realPart(1,1)%ID = 0

  do iQ2=iQ2min,iQ2max,iQ2delta
     Q2 = iQ2*0.01
     if (Q2<0.1) Q2=0.1

     do iW=iWmin,iWmax,iWdelta
        W = iW*0.01

        call eNev_init_eWQ(eNev1,eps1,W,Q2,flagOK)
        call eNev_init_Target(eNev1,TargetNuc,flagOK)

        call eNev_init_eWQ(eNev2,eps2,W,Q2,flagOK)
        call eNev_init_Target(eNev2,TargetNuc,flagOK)

        call write_electronNucleon_event(eNev1,.FALSE.,.FALSE.) 
        call eNeV_GetKinV(eNev1, nu,Q2,W,Wfree,eps,fT)
        s1=abs4(eNev1%nucleon%momentum+eNev1%lepton_in%momentum)
        x1=eNeV_Get_LightX(eNev1)
        write(*,*) 'nu :',nu
        write(*,*) 'xB :',Q2/(2*0.938*nu)
        write(*,*) 'sqrt(s) :',s1
        write(*,*) 'x  :',x1

        call write_electronNucleon_event(eNev2,.FALSE.,.FALSE.) 
        call eNeV_GetKinV(eNev2, nu,Q2,W,Wfree,eps,fT)
        s2=abs4(eNev2%nucleon%momentum+eNev2%lepton_in%momentum)
        x2=eNeV_Get_LightX(eNev2)
        write(*,*) 'nu :',nu
        write(*,*) 'xB :',Q2/(2*0.938*nu)
        write(*,*) 'sqrt(s) :',s2
        write(*,*) 'x  :',x2

!        stop

        XS_sum1=0.0
        XS_sum2=0.0

        XS_Arr_sum1=0.0
        XS_Arr_sum2=0.0
        
        iN1=0
        iN2=0

        write(*,'(2f7.3,1P,12e13.4)') W,Q2


        write(111,'(2f7.3,1P,12e13.4)') W,Q2, nu,Q2/(2*0.938*nu),&
             &s1,x1,(1-x1)*(s1**2-2*0.938**2),&
             &s2,x2,(1-x2)*(s2**2-2*0.938**2)

        isHigh = (W.ge.2.0)

        if (doMultAna) then
           write(Prefix_MultAna,'(2f8.4)') W,Q2
           call Multiplicity_Reset
        endif

        do iN=1,iNmax
           if (verbose) write(*,*) '=======iN =',iN

           if (isHigh) then
              call eventGen_eN_HiEnergy(eNev1,iN,(/1.,1.,1.,1./),finalState(1,:),channel,flagOK,XS,XS_Arr)
              XS_Arr=XS_Arr*1000
           else
              call eventGen_eN_lowEnergy(eNev1,doQE,doRes,useRes,do1Pi,do2Pi,doDIS,finalState(1,:),channel,flagOK,XS,XS_Arr)
           endif

           if (flagOK) then
              iN1=iN1+1
              XS_sum1=XS_sum1+XS
              XS_Arr_sum1=XS_Arr_sum1+XS_Arr

              !...Force all particles to decay...
              call collideMain(finalstate,realPart, 0.)
              call GarbageCollection(finalstate)

              finalstate%perweight = 1.0
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
           end if

           if (isHigh) then
              call eventGen_eN_HiEnergy(eNev2,iN,(/1.,1.,1.,1./),finalState(1,:),channel,flagOK,XS,XS_Arr)
              XS_Arr=XS_Arr*1000
           else
              call eventGen_eN_lowEnergy(eNev2,doQE,doRes,useRes,do1Pi,do2Pi,doDIS,finalState(1,:),channel,flagOK,XS,XS_Arr)
           endif

           if (flagOK) then
              iN2=iN2+1
              XS_sum2=XS_sum2+XS
              XS_Arr_sum2=XS_Arr_sum2+XS_Arr
           end if

        end do

        if (iN1.gt.0) then
           call eNeV_GetKinV(eNev1, nu,Q2,W,Wfree,eps,fT)
           if (isHigh) then
              fT = 1.0
           else
              fT = fT/ ( 1e3* pi/(eNev1%lepton_out%momentum(0)*eNev1%lepton_in%momentum(0)))
           endif
           XS_sum1=XS_sum1/(iN1*fT)
           XS_Arr_sum1=XS_Arr_sum1/(iN1*fT)
        end if

        if (iN2.gt.0) then
           call eNeV_GetKinV(eNev2, nu,Q2,W,Wfree,eps,fT)
           if (isHigh) then
              fT = 1.0
           else
              fT = fT/ ( 1e3* pi/(eNev2%lepton_out%momentum(0)*eNev2%lepton_in%momentum(0)))
           endif
           XS_sum2=XS_sum2/(iN2*fT)
           XS_Arr_sum2=XS_Arr_sum2/(iN2*fT)
        end if

        write(121,'(2f7.3,1P,12e13.4)') W,Q2, XS_sum1,XS_Arr_sum1, XS_sum2,XS_Arr_sum2

        write(122,'(2f7.3,1P,12e13.4)') W,Q2, (XS_sum1-XS_sum2)/(eps1-eps2), &
             & (XS_Arr_sum1-XS_Arr_sum2)/(eps1-eps2),&
             & (eps2*XS_sum1-eps1*XS_sum2)/(eps2-eps1), &
             & (eps2*XS_Arr_sum1-eps1*XS_Arr_sum2)/(eps2-eps1)


        call CalcParamEP(W,Q2,eps1, XS_sum1)
        call CalcParamEP(W,Q2,eps2, XS_sum2)

        write(221,'(2f7.3,1P,12e13.4)') W,Q2, XS_sum1, XS_sum2
        write(222,'(2f7.3,1P,12e13.4)') W,Q2, (XS_sum1-XS_sum2)/(eps1-eps2),&
             & (eps2*XS_sum1-eps1*XS_sum2)/(eps2-eps1)

        if (isHigh) then
           call CalcParamEP_ALLM(W,Q2,XS_sum1)
           call CalcParamEP_ALLM97(W,Q2,XS_sum2)
           call CalcParamEP_R1990(W,Q2,R1)
        
           write(301,'(2f7.3,1P,12e13.4)') W,Q2, Q2/(2*0.938*nu), XS_sum1, XS_sum2, R1
        endif

        if (doMultAna) call Multiplicity_Write(Prefix_MultAna)

     end do
     
     write(111,*)
     write(111,*)

     write(121,*)
     write(121,*)
     write(122,*)
     write(122,*)

     write(221,*)
     write(221,*)
     write(222,*)
     write(222,*)

     write(301,*)
     write(301,*)

  end do

  write(*,*) 'Done.'


end program plotAllXS
