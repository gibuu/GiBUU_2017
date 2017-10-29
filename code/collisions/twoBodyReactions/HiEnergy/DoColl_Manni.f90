!******************************************************************************
!****m* /Coll_Manni
! PURPOSE
! Implement meson + baryon -- annihilation processes.
!******************************************************************************
module Coll_Manni

  IMPLICIT NONE
  private

  integer,save :: parton(2)         ! outgoing quark/diquark
  integer,save :: KFBar,KFMes       ! ingoing baryon and meson

  public :: DoColl_ManniProb, DoColl_Manni, DoColl_ManniCheck

  logical, parameter :: flag_Geiss = .false.
  ! if .true., use the energy dependent
  ! strangeness suppression factor from
  ! J. Geiss et al., NPA 644, 107 (1998)

  !****************************************************************************
  !****g* Coll_Manni/angDistribution
  ! SOURCE
  !
  integer, save :: angDistribution = 2
  ! PURPOSE
  ! Switch to select the angular distribution:
  ! * 1: isotropic
  ! * 2: diquark/quark aligned like baryon/meson
  !****************************************************************************

  !****************************************************************************
  !****g* Coll_Manni/itry_max
  ! SOURCE
  integer, save :: itry_max = 10
  ! PURPOSE
  ! maximum number of tries
  !****************************************************************************


  logical, save:: initFlag=.true.

contains

  !****************************************************************************
  !****s* Coll_Manni/readInput
  ! NAME
  ! subroutine readInput
  !
  ! PURPOSE
  ! Reads input in jobcard out of namelist "coll_Manni"
  !****************************************************************************
  subroutine readInput
    use output, only: Write_ReadingInput

    !**************************************************************************
    !****n* Coll_Manni/coll_Manni
    ! NAME
    ! NAMELIST /coll_Manni/
    ! PURPOSE
    ! Includes the switches:
    ! * angDistribution
    ! * itry_max
    !**************************************************************************
    NAMELIST /coll_Manni/ angDistribution,itry_max

    integer :: ios

    if (.not.initFlag) return

    call Write_ReadingInput('coll_Manni',0)
    rewind(5)
    read(5,nml=coll_Manni,IOSTAT=ios)
    call Write_ReadingInput('coll_Manni',0,ios)

    write(*,*) 'angDistribution = ',angDistribution
    write(*,*) 'itry_max        = ',itry_max

    call Write_ReadingInput('coll_Manni',1)

    initFlag = .false.

  end subroutine readInput



  !****************************************************************************
  !****f* Coll_Manni/DoColl_ManniProb
  ! NAME
  ! real function DoColl_ManniProb (inPart, srtS)
  !
  ! PURPOSE
  ! Calculate the probability for the Manni process.
  ! Reference: M. Wagner et al, Phys. Rev. C71 (2005) 034910
  !****************************************************************************
  real function DoColl_ManniProb (inPart, srtS)
    use particleDefinition
    use IDtable, only: pion,rho,omegaMeson,sigmaMeson,eta,etaPrime,isBaryon
    use ID_translation, only: KFfromBUU, SplitBaryon, SplitMeson
    use random, only: rn

    type(particle),dimension(:),intent(in)   :: inPart   ! incoming particles
    real,                       intent(in)   :: srtS

    integer :: iiM,iiB, IDBar,IDMes
    integer :: iqm(2),iqb(3) ! quark contents
    integer :: aq,rq,j,mi,hq
    logical :: manniflag,mwflag
    real    :: x,fac
    integer :: mq,mdiq           ! outgoing quark/diquark

    if (initFlag) call readInput

    DoColl_ManniProb = 0

    ! decide which particle is the meson

    iiB = 2
    if (isBaryon(inPart(1)%ID)) iiB=1
    iiM = 3 - iiB ! iiB=1 -> iiM=2, iiB=2 -> iiM=1

    IDBar = inPart(iiB)%ID
    if (inPart(iiB)%antiparticle) IDBar = -IDBar
    IDMes = inPart(iiM)%ID
    if (inPart(iiM)%antiparticle) IDMes = -IDMes


    ! no charmed particles:

    if (abs(IDMes).ge.114) return
    if (abs(IDBar).gt.55)  return

    ! convert BUU to KF and split into quark contents:

    KFBar = KFfromBUU (IDBar,inPart(iiB)%charge)
    KFMes = KFfromBUU (IDMes,inPart(iiM)%charge)

    call SplitBaryon(KFBar, iqb(1), iqb(2), iqb(3))
    call SplitMeson(KFMes, iqm(1), iqm(2))

    ! SplitMeson liefert immer (/1,-1/) und nie höhere Werte

    ! correct quark contents:

    x = rn()

    if (inPart(iiM)%charge == 0) then
       select case (IDMes)
       case (pion,rho,omegaMeson)
          if (x>0.5) iqm = (/2,-2/)
       end select
    end if

    select case (IDMes)
    case (sigmaMeson)
       if (x>0.5) then
          iqm = (/2,-2/)
       else
          iqm = (/1,-1/)
       end if
    case (eta)
       if (x<1./6.) then
          iqm = (/1,-1/)
       else if (x<2./6.) then
          iqm = (/2,-2/)
       else
          iqm = (/3,-3/)
       end if
    case (etaPrime)
       if (x<1./3.) then
          iqm = (/1,-1/)
       else if (x<2./3.) then
          iqm = (/2,-2/)
       else
          iqm = (/3,-3/)
       end if
    end select

    ! find quark antiquark to annihilate if possible

    if (iqm(1)<0) then
       aq=-iqm(1)
       rq=iqm(2)
    else
       aq=-iqm(2)
       rq=iqm(1)
    end if

    manniflag=.false.
    mwflag=.true.
    j=1
    do while(mwflag.and.j.lt.4)
       if (iqb(j).eq.aq) then
          manniflag=.true.
          mwflag=.false.
          mi=j
          iqb(j)=rq
       end if
       j=j+1
    end do

    if (.not.manniflag) return  ! no annihilation possible

    if (mi.ne.3) then
       hq=iqb(3)
       iqb(3)=iqb(mi)
       iqb(mi)=hq
    end if
    if (iqb(2).gt.iqb(1)) then
       hq=iqb(1)
       iqb(1)=iqb(2)
       iqb(2)=hq
    end if

    ! now the remaining three quarks are stored in q with the baryon quarks first

    ! double the cross-section
    ! if annihilation would have been possible in different ways

    fac=1.
    if (manniflag) then
       do j=1,2
          if (iqb(j).eq.aq) then
             fac=2.
          end if
       end do
    end if

    ! determine diquark and quark for jetset

    mdiq=1000*iqb(1)+100*iqb(2)+3
    mq=iqb(3)

    parton(iiB) = mdiq
    parton(iiM) = mq

    ! manni-cross-section(fraction):

    if ( .not.flag_Geiss ) then
      DoColl_ManniProb = fac/2.*max(1.2-0.2*srtS,0.) ! strangeness suppression 0.3
    else
      DoColl_ManniProb = max(0.85-0.17*srtS,0.) ! Geiss strangeness suppr.
    end if

!    DoColl_ManniProb = DoColl_ManniProb * fac

  end function DoColl_ManniProb


  !****************************************************************************
  !****s* Coll_Manni/DoColl_Manni
  ! NAME
  ! subroutine DoColl_Manni(inPart,outPart,flagOK, srtS,pcm,beta)
  !
  ! PURPOSE
  ! Perform the Manni-Process
  !
  ! INPUTS
  ! cf. DoColl_Pythia
  !
  ! OUTPUT
  ! cf. DoColl_Pythia
  !
  !****************************************************************************
  subroutine DoColl_Manni(inPart,outPart,flagOK, srtS,pcm,beta)
    use particleDefinition
    use constants, only: twoPi
    use random, only: rn
    use CollTools, only: SetSomeDefaults_PY, SetVectorFromPYJETS
    use hadronFormation, only: useJetSetVec
    use CollGetLeading, only: GetLeading_PY
    use twoBodyTools, only: IsElastic

    type(particle),dimension(:),intent(in)   :: inPart   ! incoming particles
    type(particle),dimension(:),intent(inout):: outPart  ! outgoing particles
    real,                       intent(in)   :: srtS
    real, dimension(0:3),       intent(in)   :: pcm
    real, dimension(1:3),       intent(in)   :: beta
    logical,                    intent(out)  :: flagOK

    integer :: iTry!,i,j
    real :: theta, phi
    real :: a, b

    COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
    integer MSTU,MSTJ
    double precision PARU,PARJ
    SAVE /PYDAT1/

    COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
    integer KCHG
    double precision PMAS,PARF,VCKM
    SAVE /PYDAT2/

    COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
    integer MDCY,MDME,KFDP
    double precision BRAT
    SAVE /PYDAT3/

    COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
    integer MSTP,MSTI
    double precision PARP,PARI
    SAVE /PYPARS/

    COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
    integer MSEL,MSELPD,MSUB,KFIN
    double precision CKIN
    SAVE /PYSUBS/

    COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
    integer N,NPAD,K
    double precision P,V
    SAVE /PYJETS/

    COMMON/PYINT1/MINT(400),VINT(400)
    integer MINT
    double precision VINT
    SAVE /PYINT1/


!!$    write(*,*) 'srtrs:',srts
!!$    call WriteParticle(6,99,1,inPart(1))
!!$    call WriteParticle(6,99,2,inPart(2))


    !...set up PYTHIA/JETSET:

    call SetSomeDefaults_PY

    flagOK = .FALSE.
    iTry = 0

    ! Set the strangeness suppression factor P(s)/P(u):
    if ( .not.flag_Geiss ) then
       parj(2)=0.30
    else
       a = -0.006666
       b = 0.433
       parj(2) = min(max(a*srts+b,0.3),0.4)
    end if

    do
       outPart%ID = 0 ! reset outgoing particles

       iTry=iTry+1
       if (iTry.ge.itry_max) then
          write(*,'(A,i3,G12.5,6i5)') 'DoColl_Manni: failure, itry=', iTry, srts, inPart(1)%ID, inPart(2)%ID, &
                                      inPart(1)%charge, inPart(2)%charge, parton(1), parton(2)
          return
       end if

       !... Generate THE EVENT:

       if (useJetSetVec) call GetJetsetVecINIT
       MINT(51) = 0     ! reset error flag

       call PY2ENT(0, parton(1),parton(2), srtS)

       if (MINT(51).eq.2) cycle

       MSTI(4) = 0                ! start search line in GetLeading_PY

       if (useJetSetVec) then
          call GetJetsetVec(.TRUE.)
!          call PYLIST(2)
!          call GetJetSetVec_List(6,1,N)

          call GetJetsetVecCheckT(-1d-5)

          call GetJetsetVecPYEDIT
       end if

       call GetLeading_PY         ! find leading particles

!!$       call PYLIST(2)

       call PYEDIT(1)             ! clean up event list

       !...Boost and rotate(randomly)(lu2ent puts first parton o z-axis)
       !...the whole event to final system

       if (angDistribution.eq.1) then ! isotropic
          phi=twopi*rn()
          theta=acos(2.*(rn()-0.5))
       else
          phi = atan2(pcm(2),pcm(1))
          theta = atan2(sqrt(pcm(1)**2+pcm(2)**2),pcm(3))
       end if

       call PYROBO(1,N,theta, phi, beta(1),beta(2),beta(3))

!!$       call PYLIST(2)
!!$       call PYEDIT(1)             ! clean up event list

       if (useJetSetVec) then
          call GetJetsetVecPYROBO(theta,phi, beta(1),beta(2),beta(3))

!!$       call PYLIST(2)
!!$       call GetJetSetVec_List(6,1,N)
!!$       stop
       end if

       !...Copy Particles to ouput-vector

       call SetVectorFromPYJETS(outPart, 0.0) !!!! HardScaleQ2 not yet set !!!!

       !...Test for elastic event
       if ((N==2) .and. IsElastic(inPart(1:2),outPart(1:2))) cycle

       !...exit the loop
       exit

    end do

    flagOK = .TRUE.

!!$    do i=1,30
!!$       if (OutPart(i)%ID > 0) call WriteParticle(6,99,i,outPart(i))
!!$    enddo
!!$
!!$    N=2
!!$    do j=1,2
!!$       do i=1,3
!!$          P(j,i) = inPart(j)%momentum(i)
!!$          P(j,4) = inPart(j)%momentum(0)
!!$          P(j,5) = inPart(j)%mass
!!$       end do
!!$       K(j,1)=1
!!$       K(j,2)=92
!!$    end do
!!$    call PYLIST(2)
!!$    call PYROBO(1,N,0.0, 0.0, -beta(1),-beta(2),-beta(3))
!!$    call PYROBO(1,N,0.0, -phi,0d0,0d0,0d0)
!!$    call PYROBO(1,N,-theta,0d0,0d0,0d0,0d0)
!!$    call PYLIST(2)
!!$    stop

  end subroutine DoColl_Manni

  !****************************************************************************
  !****s* Coll_Manni/DoColl_ManniCheck
  ! NAME
  ! subroutine DoColl_ManniCheck(inPart,outPart,flagOK)
  !
  ! PURPOSE
  ! Perform Checks on the Manni-Process: charge conservation,
  ! energy conservation, etc
  !
  !****************************************************************************
  subroutine DoColl_ManniCheck(inPart,outPart,flagOK)
    use particleDefinition
    use output, only: WriteParticle

    type(particle),dimension(:),intent(in)   :: inPart   ! incoming particles
    type(particle),dimension(:),intent(in)   :: outPart  ! outgoing particles
    logical,                    intent(out)  :: flagOK

    integer :: i,cIn, cOut

    flagOK = .TRUE.

    !...Charge conservation:

    cIn = 0
    do i=1,2
       cIn = cIn+inPart(i)%charge
    end do
    cOut = 0
    do i=1,size(outPart)
       if (outPart(i)%ID.ne.0) cOut = cOut+outPart(i)%charge
    end do

    if (cIn.ne.cOut) then
       write(*,*) 'DoColl_ManniCheck: charge !!!', cIn,cOut

       call WriteParticle(6,1,inPart)
       call WriteParticle(6,2,outPart)

       call PYLIST(2)

       call PYGIVE('MINT(51)=')

       stop
       flagOK = .FALSE.
    end if

    !...Energy conservation

    ! to be implemented...


  end subroutine DoColl_ManniCheck

  !****************************************************************************

end module Coll_Manni
