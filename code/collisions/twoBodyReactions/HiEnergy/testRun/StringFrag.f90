program main
  use inputGeneral
  use particleDefinition
!  use particleProperties
  use particleProperties, only: initParticleProperties
  use CollTools
  use PythiaSpecFunc
  use version
  use hadronFormation
  use output

  implicit none


  integer :: iEv,nEv,i,iS,nS,i2
  real :: sS

  integer :: nPart(0:3),nPartMax(0:3),nPartAve(0:3)
  integer :: npi(-1:4),nN,nOmega,nnOmega
  integer :: nPartRho(0:11)

  logical :: isRho0

  COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
  integer N,NPAD,K
  double precision P,V
  SAVE /PYJETS/

  COMMON/PYINT1/MINT(400),VINT(400)
  integer MINT
  double precision VINT
  SAVE /PYINT1/

  COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
  integer MSTU,MSTJ
  double precision PARU,PARJ
  SAVE /PYDAT1/


  call PrintVersion
!!$  call SetSomeDefaults_PY
!!$  call PYLIST(12)
!!$  stop

  call readinputGeneral
  call forceInitFormation
  call InitParticleProperties

!!$  call SetSomeDefaults_PY
!!$
!!$!  call readinputGeneral
!!$  call forceInitFormation
!!$  call InitParticleProperties


!  call readinputGeneral
!  call init_database
!  call forceInitFormation
!!  call InitParticleProperties


!  call SetSomeDefaults_PY
!  call PYLIST(12)
!  call LUGIVE('MSTU(15)=2')
!  call SetSomeDefaults_Fr
!  call LULIST(12)
!  stop

!  call SetSwitchPythiaHermes(.true.)

  call SetSomeDefaults_PY

  DoPrLevel(1) = .FALSE.
  DoPrLevel(2) = .FALSE.
  
!  call PYGIVE('MSTJ(21)=2')
  call PYGIVE('MSTJ(21)=0')

  nS = 200

!  do iS=3,nS
!     sS = iS*0.1

!  do iS=110,400
!     sS = iS*0.01

  do iS=140,2000,10
     sS = iS*0.01



     nPartMax = 0
     nPartAve = 0
     npi = 0
     nnomega = 0
     nPartRho = 0

     nEv = 100000
!     nEv = 10000
!     nEv = 1000
!     nEv = 10

     call Init_VM_Mass(sS)

     do iEv=1,nEv

        if (useJetSetVec) call GetJetsetVecINIT

!        call PY2ENT(0, 1,-1, sS) ! u+ubar string

!        call PY2ENT(0, 2101,2, sS) ! proton-like
        call PY2ENT(0, 2103,2, sS) ! proton-like
!        call PY2ENT(0, 2203,1, sS) ! proton-like
        
!        call PYLIST(2)
!        call PYEDIT(1)
!        call PYLIST(2)
!        stop
        
        if (MINT(51).eq.2) cycle

        nPart=0

        do i=1,N
           if (K(i,1).ge.10) cycle
           if (abs(K(i,2)).lt.100) cycle
           if (abs(K(i,2)).lt.1000) then
              nPart(1)=nPart(1)+1
           else if (abs(K(i,2)).lt.10000) then
              if (K(i,2)>0) then
                 nPart(2)=nPart(2)+1
              else
                 nPart(3)=nPart(3)+1
              endif
           endif
        enddo
        nPart(0)=sum(nPart(1:3))

        nPartAve = nPartAve+nPart
        do i=1,3
           if (nPart(i) > nPartMax(i)) nPartMax(i) = nPart(i)
        enddo

!        if (N.gt.2) then
!           call PYLIST(1)
!           stop
!        end if

        if (MSTJ(21).eq.2) then

! The following works only correctly, if we set MSTJ(21) = 2 in
! the namelist pythia, i.e. switch decays on

           nN = 0
           do i=1,N
              if (K(i,1).ge.10) cycle
              select case(K(i,2))
              case(2112,2212)
                 nN = nN+1
              end select
           end do
           
           nOmega=0
           if (nN.eq.1) then
              i2 = 0
              do i=1,N
                 if (K(i,1).ge.10) cycle
                 select case(K(i,2))
                 case(111,211)
                    i2 = i2+1
                 case(223)
                    nOmega=nOmega+1
                 end select
              end do
              if (i2.lt.4) then
                 npi(i2) = npi(i2)+1
              else
                 npi(4) = npi(4)+1
              end if
              if (nOmega.gt.0) nnOmega=nnOmega+1
           else
              npi(-1)=npi(-1)+1
           end if

        end if

!!$        if (i2.eq.0) then
!!$           call PYLIST(2)
!!$        end if

        if (MSTJ(21).ne.2) then
! counting how many particles in addition to a rho0:

           
           isRho0 = .false.
           do i=1,N
              if (K(i,1).ge.10) cycle
              if (K(i,2).ne.113) cycle
              isRho0=.true.
              exit
           end do
           if (isRho0) then
              if (nPart(0).gt.11) then
                 nPartRho(11)=nPartRho(11)+1
              else
                 nPartRho(nPart(0))=nPartRho(nPart(0))+1
              end if
           end if

        end if

     enddo

     nPartRho(0)=sum(nPartRho(1:11))

     write(71,'(20f13.5)') sS,nPartAve/(1.0*nEv)
     write(72,'(20f13.5)') sS,nPartMax*1.0
     write(73,'(20f13.5)') sS,npi/(1.0*nEv),nnOmega/(1.0*nEv)
     write(74,'(20f13.5)') sS,nPartRho(0)/(1.0*nEv),nPartRho(1:11)/(1.0*nPartRho(0))

     write(6,'(A,20f13.5)') 'percentage of piN: ',sS,npi/(1.0*nEv),nnOmega/(1.0*nEv)
  enddo

end program main
