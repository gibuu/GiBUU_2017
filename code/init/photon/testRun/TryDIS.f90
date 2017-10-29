program TryDIS

  use inputGeneral
  use particleDefinition
  use particleProperties, only: initParticleProperties
  use output
  use Coll_gammaN
  use PythiaSpecFunc
  use hadronFormation, only : forceInitFormation
  use Coll_Pythia
  use CollTools
  use collisionTerm
  use PILCollected
  use CallStack


  implicit none

  type(particle),dimension(1,1)  :: inPart               ! incoming nucleon
  type(particle),dimension(1,10) :: outPart  ! outgoing particles
  logical                      :: flagOK
  real                         :: W
  real                         :: Q2
  real                         :: eps
  real, dimension(1:4)         :: rVMD
  real, dimension(0:3)         :: pcm
  real, dimension(1:3)         :: beta
  logical                      :: DoDifr
  real, dimension(0:4)         :: Cross
  integer                      :: EventClass

  integer :: iMC
  integer, parameter :: nMC = 1000

  
  integer :: i, nN, npi
  real, dimension(-1:3) :: XS
  real :: nu


  call readinputGeneral
  call initParticleProperties
  call forceInitFormation
!  call InitParticleProperties

  call SetSomeDefaults_PY

  inPart%ID    = 1
  inPart%Charge= 1
  inPart%mass  = 0.938
  inPart%momentum(0) = 0.938

  call WriteParticle(6,1,1,inpart(1,1))

!  call DoSingleEvents
  call DoGrid

contains

  subroutine DoGrid
    integer :: iQ,iW

    eps = 0.9

    pcm = 0
    beta = 0
    DoDifr = .false.

    do iQ=1,45
       !    do iQ=22,45

       Q2=iQ*0.1

       do iW=110,200,2
          !  do iW=190,200,2
          W = real(iW)*0.01
          if (W.ge.2.0) W = 1.99

          nu = (W**2+Q2-0.938**2)/(2*0.938)

          XS = 0.0

          if (W.le.1.5) then
             write(*,'(0P,4f12.4,1P,5e12.3)') Q2,nu,W,eps,XS/nMC
             write(127,'(0P,4f12.4,1P,5e12.3)') Q2,nu,W,eps,XS/nMC
             cycle
          endif

          call Init_VM_Mass(W)

          call PILCollected_ZERO

          do iMC=1,nMC
             call TRACEBACK("DoColl_gammaN_Py changed.")
!             call DoColl_gammaN_Py(inPart(1,1),outPart(1,:),flagOK, W,Q2,eps,rVMD, pcm,beta, DoDifr, Cross,EventClass,minW=1.5)
             if (.not.flagOK) cycle
             call collideMain(outPart,inPart,0.0)

             nN = 0
             npi = 0
             do i=1,10
                select case (outpart(1,i)%ID)
                case (1)
                   nN = nN+1
                case (101)
                   npi = npi+1
                end select
             enddo

             XS(-1) = XS(-1) + Cross(0)
             if (nN.eq.1) then
                if (npi.gt.3) npi=3
                XS(npi) = XS(npi)+Cross(0)
             end if

          end do
          write(*,'(0P,4f12.4,1P,5e12.3)') Q2,nu,W,eps,XS/nMC
          write(127,'(0P,4f12.4,1P,5e12.3)') Q2,nu,W,eps,XS/nMC
       end do
       write(*,*)
       write(127,*)
    end do
        

  end subroutine DoGrid

  subroutine DoSingleEvents

    W = 1.7
    Q2 = 1.3
    eps =0.9

    pcm = 0
    beta = 0
    DoDifr = .false.

    call Init_VM_Mass(W)
    XS = 0.0

    do iMC=1,nMC

       call DoColl_gammaN_Py(inPart(1,1),outPart(1,:),flagOK, W,Q2,eps,rVMD, pcm,beta, DoDifr, Cross,EventClass,minW=1.5)
       !     call DoColl_gammaN_Py(inPart,outPart,flagOK, W,Q2,eps,rVMD, pcm,beta, DoDifr, Cross,EventClass)

!       write(*,*) ">>>",EventClass,Cross(0)
!       call WriteParticle(6,1,outPart(1,:))

       call collideMain(outPart,inPart,0.0)
!       call WriteParticle(6,1,outPart(1,:))


       nN = 0
       npi = 0
       do i=1,10
          select case (outpart(1,i)%ID)
          case (1)
             nN = nN+1
          case (101)
             npi = npi+1
          end select
       enddo

       XS(-1) = XS(-1) + Cross(0)
       if (nN.eq.1) then
          if (npi.gt.3) npi=3
          XS(npi) = XS(npi)+Cross(0)
       end if


    end do

    write(*,'(1P,5e12.3)') XS/nMC

  end subroutine DoSingleEvents

end program TryDIS
