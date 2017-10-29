!******************************************************************************
!****m* /Coll_BaB
! NAME
! module Coll_BaB
! PURPOSE
! Implement Baryon-AntiBaryon Annihilation (for HiEnergy).
!******************************************************************************
module Coll_BaB

  IMPLICIT NONE
  private

  !****************************************************************************
  !****g* Coll_BaB/iset
  ! SOURCE
  !
  integer, save :: iset = 1
  ! PURPOSE
  ! Switch to choose an initialization of jets:
  ! * 1: phase space distribution, also the charge is conserved (new prescription)
  ! * 2: first jet along inPart(1) momentum, 3-d jet opposite,
  !      others orthogonal, charge is not conserved (old prescription)
  !****************************************************************************

  logical, save :: initFlag=.true.
  logical, parameter :: debugFlag=.false.

  public :: DoColl_BaB

contains


  !****************************************************************************
  !****s* Coll_BaB/readInput
  ! NAME
  ! subroutine readInput
  !
  ! PURPOSE
  ! Reads input in jobcard out of namelist "coll_BaB"
  !****************************************************************************
  subroutine readInput
    use output, only: Write_ReadingInput

    !**************************************************************************
    !****n* Coll_BaB/coll_BaB
    ! NAME
    ! NAMELIST /coll_BaB/
    ! PURPOSE
    ! Includes the switches:
    ! * iset
    !**************************************************************************
    NAMELIST /coll_BaB/ iset
    integer :: ios

    if (.not.initFlag) return

    call Write_ReadingInput('coll_BaB',0)
    rewind(5)
    read(5,nml=coll_BaB,IOSTAT=ios)
    call Write_ReadingInput('coll_BaB',0,ios)

    write(*,*) 'iset = ',iset

    call Write_ReadingInput('coll_BaB',1)

    initFlag = .false.

  end subroutine readInput


  !****************************************************************************
  !****s* Coll_BaB/DoColl_BaB
  ! NAME
  ! subroutine DoColl_BaB (inPart, outPart, flagOK, sqrtS, pcm, beta)
  !
  ! PURPOSE
  ! Perform a collision of particles given in "inPart" with energy "sqrtS" and
  ! return outgoing particles in "outPart".
  !
  ! "pcm" and "beta" are vectors used for Boost and Rotation of the event.
  !
  ! if "flagOK" is false, no event happened, the output in "outPart" should
  ! be neglected!
  !
  ! INPUTS
  ! * type(particle),dimension(:) :: inPart   -- incoming particles
  ! * real                        :: sqrtS    -- energy of ollision
  ! * real, dimension(0:3)        :: pcm      -- boost-vector
  ! * real, dimension(1:3)        :: beta     -- boost vector
  !
  ! OUTPUT
  ! * type(particle),dimension(:) :: outPart  ! outgoing particles
  ! * logical                     :: flagOK   ! event okay ?
  !
  ! NOTES
  ! cf. DoColl_Pythia
  !
  ! in order to understand the meaning of "pcm" and "beta":
  ! The (Pythia-)event is done in the restframe of the two given particles.
  ! Then a call to PYROBO according
  !       phi = atan2(pcm(2),pcm(1))
  !       theta = atan2(sqrt(pcm(1)**2+pcm(2)**2),pcm(3))
  !       call PYROBO(1,N, theta,phi, beta(1),beta(2),beta(3))
  ! is performed in order to transform the system into the desired
  ! (Lab-) system.
  !****************************************************************************
  subroutine DoColl_BaB (inPart, outPart, flagOK, sqrtS, pcm, beta)
    use particleDefinition, only: particle
    use IdTable, only: nucleon, Delta
    use constants, only: twopi
    use random, only: rn
    use CollTools, only: setSomeDefaults_PY, CheckUndecayedString, SetVectorFromPYJETS
    use hadronFormation, only: useJetSetVec
    use CollGetLeading, only: GetLeading_PY
    use nBodyPhaseSpace, only: momenta_in_4BodyPS
    use twoBodyTools, only: IsElastic, IsChargeExchange
    use ID_translation, only: KFfromBUU, SplitBaryon

    type(particle),dimension(:),intent(in)   :: inPart   ! incoming particles
    type(particle),dimension(:),intent(inout):: outPart  ! outgoing particles
    real,                       intent(in)   :: sqrtS
    real, dimension(0:3),       intent(in)   :: pcm
    real, dimension(1:3),       intent(in)   :: beta
    logical,                    intent(out)  :: flagOK

    integer :: i,j,iTry,i0,j0
    integer :: ID1,ID2, IZ1,IZ2
    integer :: KF1, KF2          ! KF-codes of incoming particles

    integer :: totS,totC         ! sums for incoming hadrons:
                                 !     strangeness, charm
    integer :: totAbsS, totAbsC  ! sums of absolute values
    integer :: iQ(2,3)           ! quark contents of inc. hadrons
    integer :: oQ1,oQ2, oQ3,oQ4  ! kf-codes of outgoing quarks
    real, dimension(4) :: mass   ! masses of outgoing quarks
    real, dimension(3,4) :: impuls ! three-momenta of outgoing quarks
    real, dimension(4) :: e      ! energies of outgoing quarks
    real :: x1,x2,x4,x12,x14     ! energy fractions of outgoing quarks

    real :: phii, theta, cost
    real :: hhh(2)
    real, parameter :: DeltaMM = 0.500**2-0.330**2

    real :: PYMASS


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


    if (initFlag) call readInput

    flagOK = .FALSE.
    outPart%ID = 0 ! reset outgoing particles


    !...Set ID etc

    if (inPart(1)%antiparticle) then ! (Shift antiparticle to second position)
       ID1=inPart(2)%ID
       ID2=-inPart(1)%ID
       IZ1=inPart(2)%charge
       IZ2=inPart(1)%charge
    else
       ID1=inPart(1)%ID
       ID2=-inPart(2)%ID
       IZ1=inPart(1)%charge
       IZ2=inPart(2)%charge
    end if

    !... check whether mesons:

    if (abs(ID1).gt.100 .or. abs(ID2).gt.100) then
       write(*,1000) 'not for mesons:',ID1,IZ1,ID2,IZ2
       stop
    end if

    !... check for charm and double strange:

    if (abs(ID1).ge.53 .or. abs(ID2).ge.53) then
       write(*,1000) 'not for charm or double-strange:',ID1,IZ1,ID2,IZ2
       return
    end if


    !...set up PYTHIA/JETSET:

    call SetSomeDefaults_PY

    !... Get KF and quark content of incoming hadrons:

    KF1 = KFfromBUU (ID1,IZ1)
    KF2 = KFFromBUU (ID2,IZ2)

    call SplitBaryon(KF1, iQ(1,1), iQ(1,2), iQ(1,3))
    call SplitBaryon(KF2, iQ(2,1), iQ(2,2), iQ(2,3))

    if (debugFlag) then
       write(*,*) ' 1-st baryon: ', KF1, iQ(1,1), iQ(1,2), iQ(1,3)
       write(*,*) ' 2-nd baryon: ', KF2, iQ(2,1), iQ(2,2), iQ(2,3)
    end if

    !... Calculate total charge, strangeness, charm of incoming hadrons:

    totAbsS = 0
    totAbsC = 0
    totS = 0
    totC = 0
    do i=1,2
       do j=1,3
          if (abs(iQ(i,j)).eq.3) then
             totAbsS = totAbsS +1
             totS = totS + iQ(i,j)
          else if (abs(iQ(i,j)).eq.4) then
             totAbsC = totAbsC +1
             totC = totC + iQ(i,j)
          end if
       end do
    end do
    totS = totS/3
    totC = totC/4

    if (totAbsC.ne.0) then
       write(*,1000) 'not for charm:',ID1,IZ1,ID2,IZ2
       return
    end if

    !... Set outgoing quarks: (+kinematics)

    select case (iset)

    case (1)  ! Isotropic angular distribution and also modified flavor setting procedure
             ! (similar to DoColl_Manni.f90)

       i0=0
       i_loop : do i=1,3
         do j=1,3
           if (iq(2,j).eq.-iq(1,i)) then
              i0=i
              j0=j
              exit i_loop
           end if
         end do
       end do i_loop

       if (debugFlag) write(*,*) ' i0, j0: ', i0, j0

       if (i0.eq.0) return  ! Annihilation impossible

       oQ1=0
       oQ3=0
       do i=1,3
         if (i.ne.i0) then
           if (oQ1.eq.0) then
             oQ1=iq(1,i)
           else if (oQ3.eq.0) then
             oQ3=iq(1,i)
           end if
         end if
       end do

       oQ2=0
       oQ4=0
       do j=1,3
         if (j.ne.j0) then
           if (oQ2.eq.0) then
             oQ2=iq(2,j)
           else if (oQ4.eq.0) then
             oQ4=iq(2,j)
           end if
         end if
       end do

       if (debugFlag) write(*,*) 'oQ: ',oq1,oq2,oq3,oq4

!      Phase space redistribution of the outgoing quark momenta:

       mass(1)=PYMASS(abs(oQ1))
       mass(2)=PYMASS(abs(oQ2))
       mass(3)=PYMASS(abs(oQ3))
       mass(4)=PYMASS(abs(oQ4))

       impuls = momenta_in_4BodyPS (sqrts, mass)

       e(1)=sqrt(mass(1)**2+dot_product(impuls(1:3,1),impuls(1:3,1)))
       e(2)=sqrt(mass(2)**2+dot_product(impuls(1:3,2),impuls(1:3,2)))
       e(4)=sqrt(mass(4)**2+dot_product(impuls(1:3,4),impuls(1:3,4)))

       x1=2.*e(1)/sqrts
       x2=2.*e(2)/sqrts
       x4=2.*e(4)/sqrts

       x12=2.*(e(1)*e(2)-dot_product(impuls(1:3,1),impuls(1:3,2)))/sqrts**2
       x14=2.*(e(1)*e(4)-dot_product(impuls(1:3,1),impuls(1:3,4)))/sqrts**2

    case (2)  ! Orthognal jet structure
             !... (the constraint is, that all momenta are the same!)

       select case (totAbsS)
       case (0) ! ... no s-quark: p/n/Delta + p~/n~/Delta~
          hhh = 0.5
          oQ1 = int(rn()+1.5)  ! 50%: 1=d-quark, 50%: 2=u-quark
          oQ2 = -oQ1
          oQ3 = int(rn()+1.5)  ! 50%: 1=d-quark, 50%: 2=u-quark
          oQ4 = -oQ3
       case (1) ! ... one s-quark
          hhh(1) = (4*DeltaMM+5*sqrtS**2-3*sqrt(sqrtS**2*(sqrtS**2+8*DeltaMM)))/32
          hhh(2) = hhh(1)+DeltaMM
          hhh = 2*sqrt(hhh)/sqrtS
          oQ1 = 2                ! u -quark
          oQ2 = -2               ! u~-Quark
          if (totS.lt.0) then
             oQ3 = 1             ! d -quark
             oQ4 = -3            ! s~-quark
          else
             oQ3 = 3             ! s -quark
             oQ4 = -1            ! d~-quark
          end if
       case (2)  ! ... two s-quarks or ss~
          hhh(1) = 0.5 - 2*DeltaMM/sqrtS**2
          hhh(2) = 0.5 + 2*DeltaMM/sqrtS**2
          if (totS.eq.0) then
             oQ1 = int(rn()+1.5) ! 50%: 1=d-quark, 50%: 2=u-quark
             oQ2 = -oQ1
             oQ3 = 3             ! s -quark
             oQ4 = -3            ! s~-quark
          else
             write(*,1000) 'not for this:',ID1,IZ1,ID2,IZ2
             return
          end if
       end select

       ! outgoing quarks are fixed

!       write(*,*) 'oQ: ',oq1,oq2,oq3,oq4

       !... Set the momenta/energies of the outgoing quarks:

       ! the following leads to:
       !    Quark 1: +z direction, energy = sqrts/4 (approx.)
       !    Quark 2: -x direction,   -"-
       !    Quark 3: -z direction,   -"-
       !    Quark 4: +x direction,   -"-

       x1 = hhh(1)
       x2 = hhh(1)
       if (abs(oQ4).eq.3) then
          x4 = hhh(2)
       else
          x4 = hhh(1)
       end if

       x12 = x1*x2/2
       x14 = x1*x4/2

    end select


    !... start the event loop::

    iTRY = 0
    do
       iTry = iTry+1
       if (iTRY.eq.100) then
          write(*,*) 'Problem in COLL_BaB: iTRY==100! :',ID1,IZ1,ID2,IZ2,sqrts
          return
       end if

       !... Generate THE EVENT:

       if (useJetSetVec) call GetJetsetVecINIT
       MINT(51) = 0     ! reset error flag

       call PY4ENT(0, oQ1,oQ2,oQ3,oQ4, sqrts, x1,x2,x4,x12,x14)
       call CheckUndecayedString()

       if (MINT(51).eq.2) cycle

       MSTI(4) = 0 ! start search line in GetLeading_PY

       call GetLeading_PY         ! find leading particles
       call PYEDIT(1)             ! clean up event list
!      call PYLIST(2)

       ! Destroy azimuthal correlations due to 4-th jet momentum
       ! lying in xz-plane  (see JETSET manual, subroutine LU4ENT):
       phii=twopi*rn()
       call PYROBO(1,N,0.,phii,0.,0.,0.)

       !...Rotate and boost the whole event to final system:

       select case (iset)

       case (1)
          phii=twopi*rn()
          cost=1.-2.*rn()
          theta=acos(cost)
       case (2)
          phii = atan2(pcm(2),pcm(1))
          theta = atan2(sqrt(pcm(1)**2+pcm(2)**2),pcm(3))
       end select

       call PYROBO(1,N, theta, phii, beta(1), beta(2), beta(3))

!      call PYLIST(2)

       !...Copy Particles to ouput-vector

       call SetVectorFromPYJETS(outPart, 0.0) !!!! HardScaleQ2 not yet set !!!!

       !...Correct Charges
       ! wird hier nicht benutzt

       if (N==2) then

          ! Test for elastic event
          if (IsElastic(InPart(1:2),OutPart(1:2))) cycle

          ! Test for charge-exchange event in antinucleon-nucleon
          ! antinucleon-delta or nucleon-antidelta collision
          if ((InPart(1)%Id+InPart(2)%Id<=nucleon+delta) .and. &
              (InPart(1)%antiParticle.neqv.InPart(2)%antiParticle) .and. &
              IsChargeExchange(InPart(1),InPart(2),OutPart(1),OutPart(2))) cycle

       end if

       !...exit the loop
       exit

    end do


    flagOK = .TRUE.



1000 FORMAT (1X,'ERROR: COLL_BaB ',A,'!',1X,'ID1(IZ1),ID2(IZ2) = ',2(i4,'(',i2,') '))

  end subroutine DoColl_BaB

  !****************************************************************************



  !****************************************************************************
  !****s* Coll_BaB/GetBalacedKT
  ! NAME
  ! subroutine GetBalancedKT(kTsquared, v)
  ! PURPOSE
  ! Select 3 transverse momenta vec{kT}_i = (kT_x, kT_y)_i randomly
  ! so that for all three |vec{kT}_i| follows a Gauss distribution and
  ! additionally all three vectors sum up to zero.
  ! INPUTS
  ! * real :: kTsquared         -- width of the distribution
  ! OUTPUT
  ! * real, dimension(3,2) :: v -- the three 2D vectors
  !****************************************************************************
!   subroutine GetBalancedKT(kTsquared, v)
!     use random, only: rnGauss
!
!     real, intent(in) :: kTsquared
!     real, dimension(3,2), intent(out) :: v
!
!     real :: sigma, sigmaS, sigmaD, Sx, Sy, Dx, Dy
!
!     sigma  = sqrt(3.0/4.0 * kTsquared)
!     sigmaS = sqrt(2.0/3.0) * sigma
!     sigmaD = sqrt(2.0)     * sigma
!
!     Sx = rnGauss(sigmaS, 0.0)
!     Sy = rnGauss(sigmaS, 0.0)
!     Dx = rnGauss(sigmaD, 0.0)
!     Dy = rnGauss(sigmaD, 0.0)
!
!     v(1,1) = (Sx+Dx)/2
!     v(1,2) = (Sy+Dy)/2
!     v(2,1) = (Sx-Dx)/2
!     v(2,2) = (Sy-Dy)/2
!     v(3,1) = -(v(1,1)+v(2,1))
!     v(3,2) = -(v(1,2)+v(2,2))
!
!   end subroutine GetBalancedKT


end module Coll_BaB
