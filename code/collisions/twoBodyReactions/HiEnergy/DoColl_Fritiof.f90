!******************************************************************************
!****m* /Coll_Fritiof
! PURPOSE
! Implement a + b -> X processes/reactions done by FRITIOF
!******************************************************************************
module Coll_Fritiof

  IMPLICIT NONE
  private
  public :: DoColl_Fritiof

  logical, parameter :: flag_Geiss = .false.   ! if .true., use the energy dependent
                                               ! strangeness suppression factor from
                                               ! J. Geiss et al., NPA 644, 107 (1998)
contains

  !****************************************************************************
  !****s* Coll_Fritiof/DoColl_Fritiof
  ! NAME
  ! subroutine DoColl_Fritiof (inPart, outPart, flagOK, sqrtS, pcm, beta)
  !
  ! PURPOSE
  ! Perform a collision of particles given in "inPart" with energy "sqrtS" and
  ! return outgoing particles in "outPart".
  ! "pcm" and "beta" are vectors used for Boost and Rotation of the event.
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
  !
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
  subroutine DoColl_Fritiof (inPart, outPart, flagOK, sqrtS, pcm, beta)
    use particleDefinition, only: particle
    use output, only: DoPR
    use IdTable, only: nucleon, Delta
    use CollTools, only: ReduceCharge, ConvertInPartFritiof, SetSomeDefaults_FR, SetVectorFromPYJETS, CorrectChargeVector
    use twoBodyTools, only: IsElastic, IsChargeExchange
    use CollGetLeading, only: GetLeading_FR
    use hadronFormation, only: useJetSetVec
    use ID_translation, only: KFfromBUU

    type(particle),dimension(:),intent(in)   :: inPart   ! incoming particles
    type(particle),dimension(:),intent(inout):: outPart  ! outgoing particles
    real,                       intent(in)   :: sqrtS
    real, dimension(0:3),       intent(in)   :: pcm
    real, dimension(1:3),       intent(in)   :: beta
    logical,                    intent(out)  :: flagOK

    integer :: iTry
    real :: theta, phi
    integer :: ID1,ID2, IZ1,IZ2
    integer :: iz1c, iz2c        ! reduced charges
    integer :: id1c, id2c        ! reduced particle IDs
    integer :: kf1, kf2          ! KF-codes of incoming particles
    integer :: DeltaQ            ! units of charges to be added

    real :: a, b

    COMMON/FRCODES/IPT(2),PACD(27),NNUC(27),NPROT(27),KCD(27),RO1(27,2),EXMA(9,2)
    integer IPT, NNUC, NPROT, KCD
    real RO1, EXMA
    character*4 PACD
    SAVE /FRCODES/

    COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
    integer MSTU,MSTJ
    real PARU,PARJ
    SAVE /LUDAT1/

    COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
    integer N,NPAD,K
    double precision P,V
    SAVE /PYJETS/


    !...set default output

    outPart%ID = 0 ! reset outgoing particles
    flagOK = .FALSE.

    !...write incoming particles

!!$    call WriteParticle(6,1,1,inPart(1))
!!$    call WriteParticle(6,1,2,inPart(2))

    ! Set the strangeness suppression factor P(s)/P(u):
    if (.not. flag_Geiss) then
       parj(2)=0.30
    else
       a = -0.006666
       b = 0.433
       parj(2) = min(max(a*sqrts+b,0.3),0.4)
    end if


    if (sqrtS<2.4 .and. ((inPart(1)%ID<100 .and. inPart(1)%ID>52) .or. &
                         (inPart(2)%ID<100 .and. inPart(2)%ID>52))) then
      write(*,*) 'DoColl_Fritiof: not for Xi,Omega,... and sqrts<2.4!'
      return
    end if

    !...Set ID etc
    ID1 = inPart(1)%ID
    ID2 = inPart(2)%ID
    if (inPart(1)%antiparticle) ID1 = -ID1
    if (inPart(2)%antiparticle) ID2 = -ID2
    IZ1 = inPart(1)%charge
    IZ2 = inPart(2)%charge

    !...reduce charge:
    IZ1c = ReduceCharge(ID1,IZ1)
    IZ2c = ReduceCharge(ID2,IZ2)

    DeltaQ = (iz1+iz2)-(iz1c+iz2c)

    !...convert input particles:
    ID1c = ConvertInPartFritiof(ID1)
    ID2c = ConvertInPartFritiof(ID2)

    !...transpose BUU-code -> KF:
    KF1 = KFfromBUU (ID1c,IZ1c)
    KF2 = KFfromBUU (ID2c,IZ2c)

    !...set up PYTHIA/JETSET:
    call SetSomeDefaults_FR

    !...set up FRITIOF
    KCD(1) = kf1              ! -- Projectile:  --
    NPROT(1) = iz1c
    EXMA(1,1) = inPart(1)%mass
    EXMA(1,2) = inPart(1)%mass

    KCD(2) = kf2              ! -- Target:  --
    NPROT(2) = iz2c
    EXMA(2,1) = inPart(2)%mass
    EXMA(2,2) = inPart(2)%mass


    !... Start the Event Loop
    iTry = 0
    do
       outPart%ID = 0 ! reset outgoing particles

       iTry=iTry+1
       if (iTry.ge.100) then
          write(*,*) 'DoColl_Fritiof: itry=',iTry
          return
       end if

       if (useJetSetVec) call GetJetsetVecINIT

       !...Generate THE EVENT:


       !       write(*,*) 'In DoColl_Fritiof:', sqrts, inPart(1)%mass, inPart(2)%mass

       call FRHADHAD('CMS','NEW1','NEW2',sqrts, inPart(1)%mass, inPart(2)%mass)

       !       call PYLIST(2)

       if (MSTU(24).ne.0) then
          if (iTRY.eq.99) write(*,*) 'DoColl_Fritiof: MSTU(24).ne.0)',MSTU(24)
          cycle
       end if

       if (useJetSetVec) then
          call GetJetsetVec(.FALSE.)
          !          call PYLIST(2)
          !          call GetJetSetVec_List(6,1,N)

          call GetJetsetVecCheckT(-1d-5)

          call GetJetsetVecPYEDIT
       end if

       call GetLeading_FR ! find leading particles

       call LUEDIT(1)     ! clean up event list

      !...Rotate and boost the whole event to final system:
       phi = atan2(pcm(2),pcm(1))
       theta = atan2(sqrt(pcm(1)**2+pcm(2)**2),pcm(3))
       call LUROBO(theta, phi, beta(1), beta(2), beta(3))

       !      write(*,*) 'LULIST'
       !      call LULIST(2)
       !      write(*,*) 'PYLIST'
       !      call PYLIST(2)


       if (useJetSetVec) then
          call GetJetsetVecPYROBO(theta,phi, beta(1),beta(2),beta(3))

          !         call PYLIST(2)
          !         call GetJetSetVec_List(6,1,N)
          !         stop
       end if


       !...Copy Particles to ouput-vector

       call SetVectorFromPYJETS(outPart, 0.0) !!!! HardScaleQ2 not yet set !!!!

       !...Correct Charges

       if (DeltaQ.ne.0) call CorrectChargeVector(outPart,DeltaQ)
       if (DeltaQ.ne.0) then
          if (DoPR(2)) write(*,*) 'DoColl_Fritiof: Charge correction failed. ReDo Event!!'

!!$         call WriteParticle(6,1,1,inPart(1))
!!$         call WriteParticle(6,1,2,inPart(2))
!!$         write(*,*) 'IN       :',ID1,ID2,IZ1,IZ2
!!$         write(*,*) 'Converted:',ID1c,ID2c,IZ1c,IZ2c,(iz1+iz2)-(iz1c+iz2c),DeltaQ
!!$         call LULIST(2)

          cycle
       end if

       if (N==2) then

          ! Test for elastic event
          if (IsElastic(InPart(1:2),OutPart(1:2))) cycle

          ! Test for charge exchange event for antinucleon-nucleon,
          ! antinucleon-delta or nucleon-antidelta collision
          if ((InPart(1)%Id+InPart(2)%Id<=nucleon+delta) .and. &
              (InPart(1)%antiParticle.neqv.InPart(2)%antiParticle) .and. &
              IsChargeExchange(InPart(1),InPart(2),OutPart(1),OutPart(2))) cycle

       end if

       !...exit the loop
       exit

    end do

    flagOK = .TRUE.

  end subroutine DoColl_Fritiof

  !****************************************************************************
end module Coll_Fritiof
