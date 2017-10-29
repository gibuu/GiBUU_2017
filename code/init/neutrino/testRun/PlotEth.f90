program PlotEth

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

  use hist2Df90

  implicit none

  Type(electronNucleon_event), save :: eNev0
  Type(electronNucleon_event), save :: eNev
  type(particle) :: TargetNuc
  type(particle), dimension(1,1) :: realPart
  type(particle),dimension(1,1:50)  :: finalState

  real :: cross
  logical :: flagOK
  integer :: iEV,nEV10, nEV100

  ! Parameters to play with:

  real :: Ebeam = 2.0
!  integer :: nEV=1000
  integer,save:: nEV=-10
  logical :: doCC = .true.
!  logical :: doCC = .false.
!  logical :: doMassless=.true.
  logical :: doMassless=.false.
  integer :: iGeneration = 2
  integer:: cTarget = 1
  real :: Wcut1 = 1.6, Wcut2 = 1.65

  NAMELIST /datatable/ nEV,doCC,doMassless,iGeneration,cTarget,Wcut1,Wcut2,Ebeam


  integer :: iEn
  real :: SumAna1,SumAna2

  integer :: iPi,nPi
  integer :: iC1,iC2, iiC, i, iBin

  type(preEvent), dimension(1:25) :: PreE
  
  type(tPreEvListEntry), dimension(1:6) :: PreELE
  real :: SumAll,Sum1pi(1:6)

  type(histogram2D) :: h2d, h2d1pi,h2dXY

  real :: cost,Eprime,X,Y,nu,Q2,W, fak

  real, dimension(0:101) :: xMin,xMax

  call readinputGeneral
  call init_database
  call forceInitFormation
!  call InitParticleProperties

  call SetSomeDefaults_PY



  !...Reading input:

  call Write_ReadingInput('datatable',0)
  rewind(5)
  read(5,nml=datatable)
  write(*,*) ' Ebeam = ',Ebeam
  write(*,'(A,3i7)') '  nEv   = ',nEv
  write(*,'(A,3i7)') '  iGen  = ',iGeneration
  write(*,*)
  write(*,'(A,3i7)') '  chargeTarget  = ',cTarget
  write(*,*)
  write(*,*) ' doCC        =',doCC
  write(*,*) ' doMassless  =',doMassless
  write(*,*)
  write(*,*) ' Wcut: ',Wcut1,Wcut2

  call Write_ReadingInput('datatable',1)
  write(*,*)

  nEv100 = nEv / 100
  nEv10 = nEv / 10

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

  call CreateHist2D(h2d, "Eprime vs cos(theta)", (/0.0,0.0/),(/1.0,5.0/),(/0.01,0.01/))
  call CreateHist2D(h2d1pi, "Eprime vs cos(theta) (1pi)", (/0.0,0.0/),(/1.0,5.0/),(/0.01,0.01/))
  call CreateHist2D(h2dXY, "x vs y", (/0.0,0.0/),(/1.0,1.0/),(/0.01,0.01/))

  xMin =  99.9
  xMax = -99.9

  !...Set up the initial kinematics:

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
     

  !...Do some Events:
     
  SumAll = 0d0
  Sum1pi = 0d0
     
  do iEV=1,nEV

     CALL WriteStatusMC(iEv,nEv,nEv100,nEv10)
     call PILCollected_ZERO

     eNev = eNev0
     call eNev_init_Target(eNev,TargetNuc,flagOK)
     if (.not.flagOK) cycle
     
     call DoColl_nuN_Py(eNev,finalState(1,:),flagOK, doCC, iGeneration, doMassless, cross)
     if (.not.flagOK) cycle

     call eNeV_GetKinV(eNev, nu,Q2,W)
     
     if (W.le.Wcut1) cycle
     if (W.le.Wcut2) cross = cross*(W-Wcut1)/(Wcut2-Wcut1) 

!        call write_electronNucleon_event(eNev,.FALSE.,.FALSE.)
!        call WriteParticle(6,1,finalState(1,:))
!        write(*,*) 'cross section: (in mb)',cross
     SumAll=SumAll+cross

     !...Force all particles to decay...
     call collideMain(finalstate,realPart, 0.0, .true.)
     call GarbageCollection(finalstate)

     !...Compare with stored preEvents:
     if (.not.CreateSortedPreEvent(finalState(1,:),preE)) then
        write(*,*) 'Error PreE'
        call WriteParticle(6,1,finalState(1,:))
        cycle
     end if

     Eprime = eNev%electron_out(0)
     cost = eNeV_Get_CostLepton(eNev)
!     call write_electronNucleon_event(eNev)
!     write(*,*) 'Eprime,cost:',Eprime,cost

     X = Q2/(2*0.938*nu)
     Y = nu/Ebeam

!     write(*,*) '====>',X,Y

     call AddHist2d(h2d, (/cost,Eprime/), cross)
     if (cost.gt.-99.0) then
        call AddHist2d(h2dXY, (/X,Y/), cross)
        
        iBin = Y*100
        if ((iBin.ge.0).and.(iBin.le.100)) then
           if (X.lt.xMin(iBin)) xMin(iBin) = X
           if (X.gt.xMax(iBin)) xMax(iBin) = X
        endif

     endif

     do iiC=1,6
        flagOK=.true.
        do i=1,5
           if (preE(i)%mass.ne.preELE(iiC)%preE(i)%mass) flagOK=.false. ! attention: abuse of mass
        end do
        if(flagOK) then
           Sum1pi(iiC)= Sum1pi(iiC)+cross
           call AddHist2d(h2d1pi, (/cost,Eprime/), cross)
           exit
        end if
     end do

  end do

  SumAll=SumAll/nEV
  Sum1pi=Sum1pi/nEV
  
  write(121,'(f12.5,1P,20e13.4)') Ebeam,SumAll,Sum(Sum1pi),Sum1pi
  write(*,'(f12.5,1P,20e13.4)') Ebeam,SumAll,Sum(Sum1pi),Sum1pi

  fak = 1.0/nEV
  call WriteHist2D_Gnuplot(h2d,201,mul=fak,add=1e-20,file="Neutrino.h2D.theta_Eprime.dat")
  call WriteHist2D_Gnuplot(h2d1pi,201,mul=fak,add=1e-20,file="Neutrino.h2D.theta_Eprime.1pi.dat")
  call WriteHist2D_Gnuplot(h2dXY,201,mul=fak,add=1e-20,file="Neutrino.h2D.X_Y.dat",dump=.true.)
  
  open(13,file="Neutrino.MinMax.dat",status="unknown")
  do iBin=0,100
     y=iBin*0.01
     write(13,'(1P,3e13.4)') y,xMin(iBin),xMax(iBin)
  enddo
  close(13)

  write(*,*) 'done.'

contains

  subroutine WriteStatusMC(iEV, nEV, nEV_100, nEv_10)
    IMPLICIT NONE
    integer iEV, nEV, nEV_100, nEv_10
    
    
    integer n_Perc
    
    if (nEV.lt.100) return  ! dont write status
    if (mod(iEV,nEV_100) .eq. 0) then
       n_Perc = iEV/nEV_100
       
         write(*,*) 'PlotEth : ',n_Perc,'%'

       open(12,file='PlotEth.run',status='unknown')
       write(12,*) 'PlotEth.run : ',n_Perc,'%'
       close(12)
        
    endif

  end subroutine WriteStatusMC

end program PlotEth
