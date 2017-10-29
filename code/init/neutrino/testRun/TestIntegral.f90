program TestIntegral

  use constants
  use inputGeneral
  use particleDefinition
  use particleProperties, only: initParticleProperties
  use hadronFormation, only : forceInitFormation
  use eN_eventDefinition
  use eN_event
  use CollTools
  use output

  use neutrinoXsection

  implicit none

  integer,parameter ::maxIP=33   !number of possible final states:
  !1: nucleon (QE)
  !2-31: non-strange baryon resonance (as in IdTable)
  !32: pi neutron-background  (e.g. nu + n -> mu + pi+ + n)
  !33: pi proton-background   (e.g. nu + n -> mu + pi0 + p)

  type(particle) :: TargetNuc
  Type(electronNucleon_event), save :: eNev0
  Type(electronNucleon_event), save :: eNev(maxIP)

  integer :: IP
  logical :: flagOK
  real :: sigma1(maxIP)
  real :: sigma2(maxIP), sigmaE2(maxIP)
  real :: sigma3(maxIP), sigmaE3(maxIP)
  integer :: nC1(maxIP),nC2(maxIP),nC3(maxIP)

  ! Exchange:

  Type(electronNucleon_event), save :: Global_eNev
  real, save :: Global_costheta
  integer, save :: Global_IP
  integer, save :: Global_nCalls


  ! Parameters to play with:

  real :: Ebeam = 2.0
  logical :: doCC = .true.
  integer :: iGeneration = 2
  integer :: cTarget = 1
  real :: cost=0.8
  integer :: charge_out,pion_charge_out

  NAMELIST /datatable/ doCC,cTarget,iGeneration,Ebeam,cost

  call readinputGeneral
  call initParticleProperties
  call forceInitFormation
!  call InitParticleProperties

  call SetSomeDefaults_PY



  !...Reading input:

  call Write_ReadingInput('datatable',0)
  rewind(5)
  read(5,nml=datatable)
  write(*,*) ' Ebeam = ',Ebeam
  write(*,*) ' cost  = ',cost
  write(*,'(A,3i7)') '  iGen  = ',iGeneration
  write(*,'(A,3i7)') '  chargeTarget  = ',cTarget
  write(*,*) ' doCC          = ',doCC
  write(*,*)
  call Write_ReadingInput('datatable',1)
  write(*,*)


  !...Set up the target nucleon:

  call setToDefault(TargetNuc)
  TargetNuc%ID = 1
  TargetNuc%charge = cTarget
  TargetNuc%mass = 0.938
  TargetNuc%momentum = (/0.938, 0.0, 0.0, 0.0 /)
  TargetNuc%Position = 9999999d0


  !...Set up the event type

  if (doCC) then
     call eNev_SetProcess(eNev0, 2,iGeneration)
  else
     call eNev_SetProcess(eNev0, 3,iGeneration)
  endif
  call eNev_init_nuStep1(eNev0,TargetNuc)
  
  sigma1 = 0
  sigma2 = 0
  sigmaE2 = 0
  sigma3 = 0
  nC2 = 0
  do IP = 1,maxIP
     
     write(*,*) 'IP=',IP,'...'

     eNev(IP) = eNev0
     call eNev_init_nuStep2(eNev(IP),Ebeam,IP,flagOK)
     if (.not.flagOK) cycle
     if(.not.SetHadronCharge(eNev(IP),IP,charge_out,pion_charge_out)) cycle

     call IntegrateSpanish(eNev(IP),IP,cost,sigma1(IP))
     nC1(IP)=Global_nCalls
!!$     call IntegrateNew(eNev(IP),IP,cost,sigma2(IP),sigmaE2(IP))
!!$     nC2(IP)=Global_nCalls
     call IntegrateNewAdap(eNev(IP),IP,cost,sigma3(IP),sigmaE3(IP))
     nC3(IP)=Global_nCalls

  end do
     
  write(*,*) '=========final result:==========='

  do IP = 1,maxIP
     write(*,'(i3,1P,1e13.5,i5,2e13.5,i5,2e13.5,i5)') IP,&
          & sigma1(IP),nC1(IP),&
          & sigma2(IP),sigmaE2(IP),nC2(IP),&
          & sigma3(IP),sigmaE3(IP),nC3(IP)
  end do

  write(*,*) 'okay.'
  

contains

  subroutine IntegrateSpanish(eNev,IP,cost,sig)
    use gauss_integration
    use idTable, only: nucleon
    implicit none

    type(electronNucleon_event), intent(inout) :: eNev
    integer, intent(in) :: IP
    real, intent(in) :: cost
    real, intent(out) :: sig

    real :: costheta
    real :: xmin,xmax
    real, dimension(:),allocatable :: yy,xx
    integer :: intprec,l,n2
    type(particle), dimension(20)   :: OutPart

    sig = 0
    Global_nCalls = 0
    write(*,*) '----- spanish:'
    
    intprec=3
    if(IP.eq.nucleon) intprec=300
    
    allocate (xx(20*intprec))
    allocate (yy(20*intprec))

    costheta = cost
    call getEleptonLimits(IP,eNev,xmin,xmax)
    call sg20r(xmin,xmax,intprec,xx,n2)

    yy = 0.0
    do l=1,n2
       call eNev_init_nuStep3a(eNev,xx(l),costheta,flagOK)
       if (.not.flagOK) cycle
       call XsecdCosthetadElepton(eNev,IP,OutPart, yy(l))
    end do
    call rg20r(xmin,xmax,intprec,yy,sig)

    deallocate(xx)
    deallocate(yy)

    Global_nCalls = 20*intprec
    
  end subroutine IntegrateSpanish

  subroutine IntegrateNew(eNev,IP,cost,sig,sigE)
    use idTable, only: nucleon
    implicit none

    type(electronNucleon_event), intent(inout) :: eNev
    integer, intent(in) :: IP
    real, intent(in) :: cost
    real, intent(out) :: sig,sigE

    real :: xmin,xmax, hh1,hh2, r1,r2

    sig = 0
    Global_nCalls = 0
    write(*,*) '----- new:'

    call getEleptonLimits(IP,eNev,xmin,xmax)

    hh1 = eNev%lepton_in%momentum(0)
    hh2 = 0.938*hh1/(0.938+(1-cost)*hh1)

    write(*,*) 'Eprime <=',hh2,hh1
    
    xmax = min(xmax,hh2)

    if(IP.eq.nucleon) then
       xmin = max(0.0, xmax - 0.3)
       xmax = min(hh1, xmax + 0.3)
    end if

    Global_eNev = eNev
    Global_costheta = cost
    Global_IP = IP

    call DGS56P(Integrand, xmin,xmax,r1,r2)

    write(*,*) r1,r2
    sig = r1
    sigE = r2
  end subroutine IntegrateNew

  subroutine IntegrateNewAdap(eNev,IP,cost,sig,sigE)
    use idTable, only: nucleon
    implicit none

    type(electronNucleon_event), intent(inout) :: eNev
    integer, intent(in) :: IP
    real, intent(in) :: cost
    real, intent(out) :: sig,sigE

    real :: xmin,xmax, hh1,hh2, r1,r2

    sig = 0
    Global_nCalls = 0
    write(*,*) '----- new Adap:'

    call getEleptonLimits(IP,eNev,xmin,xmax)

    hh1 = eNev%lepton_in%momentum(0)
    hh2 = 0.938*hh1/(0.938+(1-cost)*hh1)

    write(*,*) 'Eprime <=',hh2,hh1
    
    xmax = min(xmax,hh2)

    if(IP.eq.nucleon) then
       xmin = max(0.0, xmax - 0.3)
       xmax = min(hh1, xmax + 0.3)
    end if

    Global_eNev = eNev
    Global_costheta = cost
    Global_IP = IP

    call DADAPT(Integrand, xmin,xmax,1,0.1,0.0, r1,r2)

    write(*,*) r1,r2
    sig = r1
    sigE = r2

    write(121,*)
    write(121,*)

  end subroutine IntegrateNewAdap

  real function Integrand(x)
    implicit none
    real :: x,h
    logical :: flagOK
    type(particle), dimension(20)   :: OutPart
    
    Global_nCalls=Global_nCalls+1
    Integrand = 0.0
    call eNev_init_nuStep3a(Global_eNev,x,Global_costheta,flagOK)
    if (.not.flagOK) return
    call XsecdCosthetadElepton(Global_eNev,Global_IP,OutPart, h)
    Integrand = h

    write(121,*) x,h

    return
    
  end function Integrand
  

end program TestIntegral
