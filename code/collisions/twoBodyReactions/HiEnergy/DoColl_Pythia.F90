!******************************************************************************
!****m* /Coll_Pythia
! NAME
! module Coll_Pythia
! FUNCTION
! Implement a + b -> X processes/reactions done by PYTHIA.
!******************************************************************************
module Coll_Pythia

  IMPLICIT NONE
  private

  public :: DoColl_Pythia

contains


  !****************************************************************************
  !****s* Coll_Pythia/DoColl_Pythia
  ! NAME
  ! subroutine DoColl_Pythia(inPart, outPart, flagOK, sqrtS, pcm, beta)
  !
  ! PURPOSE
  ! Perform a collision of particles given in "inPart" with energy "sqrtS" and
  ! return outgoing particles in "outPart".
  ! "pcm" and "beta" are vectors used for boost and rotation of the event.
  ! If "flagOK" is false, no event happened and the output in "outPart" should
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
  ! In order to understand the meaning of "pcm" and "beta":
  ! The (Pythia) event is done in the rest frame of the two given particles.
  ! Then a call to PYROBO according
  !       phi = atan2(pcm(2),pcm(1))
  !       theta = atan2(sqrt(pcm(1)**2+pcm(2)**2),pcm(3))
  !       call PYROBO(1,N, theta,phi, beta(1),beta(2),beta(3))
  ! is performed in order to transform the system into the desired
  ! (lab) frame.
  !****************************************************************************
  subroutine DoColl_Pythia(inPart, outPart, flagOK, sqrtS, pcm, beta)
    use particleDefinition, only: particle
    use IdTable, only: nucleon, Delta
    use CollTools
    use twoBodyTools, only: IsElastic, IsChargeExchange
    use hadronFormation, only: useJetSetVec
    use output, only: DoPR, WriteParticle
    use CollGetLeading, only: GetLeading_PY
    use ID_translation, only: KFfromBUU, SplitKFtoQ
    use ParticleProperties, only: hadron

    type(particle),dimension(:),intent(in)   :: inPart   ! incoming particles
    type(particle),dimension(:),intent(inout):: outPart  ! outgoing particles
    real,                       intent(in)   :: sqrtS
    real, dimension(0:3),       intent(in)   :: pcm
    real, dimension(1:3),       intent(in)   :: beta
    logical,                    intent(out)  :: flagOK

    integer :: iTry, i
    real :: theta, phi
    integer :: ID1,ID2, IZ1,IZ2
    integer :: iz1c, iz2c        ! reduced charges
    integer :: id1c, id2c        ! reduced particle IDs
    integer :: kf1, kf2          ! KF-codes of incoming particles
    integer :: DeltaQ!,DeltaQ0    ! units of charges to be added
    logical :: EventIsAnti       ! particles converted to its anti ?

    type(particle) :: part1, part2 ! incoming particles

    integer :: Q1,Q2,Q3
    integer :: strangeQuarkIn, strangeQuarkOut

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


    character*20 Buf,cBeam,cTarget
    real :: EA,EB,pA


!---------------------------------------------------------------------



    !...set default output

    outPart%ID = 0 ! reset outgoing particles
    flagOK = .FALSE.
    EventIsAnti = .FALSE.


    !...set incoming particles:
     if (inPart(1)%antiparticle) then ! (Shift antiparticle to second position)
        part1 = inPart(2)
        part2 = inPart(1)
     else
        part1 = inPart(1)
        part2 = inPart(2)
     end if

     if (part1%antiparticle) then
        if (DoPr(1)) write(*,*) 'DoColl_Pythia: converting anti+anti: '
        call WriteParticle(6,0,1,inPart(1))
        call WriteParticle(6,0,2,inPart(2))

        call ConvertToAnti(part1)
        call ConvertToAnti(part2)
        EventIsAnti = .TRUE.

        call WriteParticle(6,0,1,part1)
        call WriteParticle(6,0,2,part2)
        stop
     end if


    !...no collision for baryonic resonances Xi,Omega,Sigma_c,Lambda_c,... for sqrts<3:
!     if (sqrtS<3.) then
!        if (((part1%ID<100).and.(part1%ID>52)) .or.((part2%ID<100).and.(part2%ID>52))) then
!           if (DoPr(2)) write(*,*) 'no PYTHIA collision for baryonic resonances Xi,Omega,Sigma_c,Lambda_c... for sqrts<3! ', &
!                                   part1%ID, part2%ID, sqrts
!           return
!        endif
!     endif

!!$    call WriteParticle(6,0,1,part1)
!!$    call WriteParticle(6,0,2,part2)
!!$    write(*,*) 'Sqrts:',sqrts
!!$    write(*,*) 'PCM:  ',PCM
!!$    write(*,*) 'Beta: ',BETA


100 continue

    !...Set ID etc

    ID1 = part1%ID
    ID2 = part2%ID
    if (part1%antiparticle) ID1 = -ID1
    if (part2%antiparticle) ID2 = -ID2
    IZ1 = part1%charge
    IZ2 = part2%charge

    !...reduce charge:

    IZ1c = ReduceCharge(ID1,IZ1)
    IZ2c = ReduceCharge(ID2,IZ2)

    DeltaQ = (iz1+iz2)-(iz1c+iz2c)

    !...convert input particles:

    ID1c = ConvertInPartPythia(ID1)
    ID2c = ConvertInPartPythia(ID2)

    if (ID1c == 0) return
    if (ID2c == 0) return


    if (ID2c < -31) then
!!$       write(*,*) 'DoColl_Pythia: no "X+AntiBaryon_s" collisions !!!'
!!$       call WriteParticle(6,0,1,part1)
!!$       call WriteParticle(6,0,2,part2)
!!$       write(*,*) 'Trying to do as   "AntiX+Baryon_s"          0 !!!'
!!$       write(*,*)
       if (DoPr(1)) write(*,*) 'DoColl_Pythia: Retry "X+AntiBaryon_s" as ANTI!!!'

       if (EventIsAnti) then
          if (DoPr(2)) write(*,*) 'DoColl_Pythia: Event is already ANTI. QUIT!'
          return
       end if

!!$       if ((ID1c .eq. 110).or.(ID1c .eq. 112)) then
!!$          if (DoPr(1)) write(*,*) 'DoColl_Pythia: ANTI not possible for "K+AntiBaryon_s". QUIT'
!!$          return
!!$       endif


       if (inPart(1)%antiparticle) then
          call ConvertToAnti(inPart(1),part1) ! converted particle is anti,
          call ConvertToAnti(inPart(2),part2) ! should be stored at pos 2
       else
          call ConvertToAnti(inPart(2),part1)
          call ConvertToAnti(inPart(1),part2)
       end if
       EventIsAnti = .TRUE.
       goto 100! LabelReDo
    end if


    !...transpose BUU-code -> KF:

    KF1 = KFfromBUU(ID1c,IZ1c)
    KF2 = KFfromBUU(ID2c,IZ2c)

    ! quark content has different sign!
    StrangeQuarkIn = -(sign(1,ID1c)*hadron(abs(ID1c))%strangeness + &
                       sign(1,ID2c)*hadron(abs(ID2c))%strangeness)

    !...set up PYTHIA/JETSET:

    call SetSomeDefaults_PY

    CKIN(3) = 0d0          ! min pThat (allow low pT events)
    MSTP(142) = 0          ! don't use reweighted events

    !... Initialize PYTHIA

    call PYNAME(KF1,Buf)
    write(cBeam,'(A)') Buf(1:10)
    call PYNAME(KF2,Buf)
    write(cTarget,'(A)') Buf(1:10)


    ! K^0, K^0~ -> K_L^0:
    if (KF1== 311) cBeam   = 'KL0'
    if (KF1==-311) cBeam   = 'KL0'
    if (KF2== 311) cTarget = 'KL0'
    if (KF2==-311) cTarget = 'KL0'

!    write(*,*) '--',cBeam,cTarget, sqrts
!    call WriteParticle(6,0,1,inPart(1))
!    call WriteParticle(6,0,2,inPart(2))


    call CalcMomentum(sqrts,part1%mass,part2%mass, EA,EB,pA)

    P(1,1) = 0d0
    P(1,2) = 0d0
    P(1,3) = pA
    P(1,4) = EA
    P(1,5) = part1%mass

    P(2,1) = 0d0
    P(2,2) = 0d0
    P(2,3) = -pA
    P(2,4) = EB
    P(2,5) = part2%mass


!    call PYINIT('CMS', cBeam,cTarget, sqrts)
    call PYINIT('5MOM', cBeam,cTarget, sqrts)

    !... Start the Event Loop

    iTry = 0
    do
       outPart%ID = 0 ! reset outgoing particles

       iTry=iTry+1
       if (iTry.ge.100) then
          write(*,*) 'DoColl_Pythia: itry=',iTry
          return
       end if

       if (useJetSetVec) call GetJetsetVecINIT

       !...Generate THE EVENT:

       call PYEVNT
!       call PYLIST(2)
!       stop

       if (MINT(51)==2) then
          N = 0
          ! print *,"failure in DoColl_Pythia!"
          return ! -> FAILURE
       end if

       !...Check Strangeness conservation

       strangeQuarkOut = 0
       do i=1,N
          if (K(i,1) >= 10) cycle
          call SplitKFtoQ(K(i,2), Q1,Q2,Q3)
          ! count all strange quarks
          if (abs(Q1) == 3) strangeQuarkOut = strangeQuarkOut + sign(1,Q1)
          if (abs(Q2) == 3) strangeQuarkOut = strangeQuarkOut + sign(1,Q2)
          if (abs(Q3) == 3) strangeQuarkOut = strangeQuarkOut + sign(1,Q3)
       end do

       if (strangeQuarkOut /= strangeQuarkIn) then
!!$          call PYLIST(2)
!!$          write(*,*) 'strangeQuarkIn  = ' , strangeQuarkIn
!!$          write(*,*) 'strangeQuarkOut = ' , strangeQuarkOut
!!$          stop

!!$          if (DoPr(1)) write(*,*) 'DoColl_Pythia: redo event, strangness not conserved.'

          cycle ! just redo the event to ensure strangeness conservation
       end if


       if (useJetSetVec) then
          call GetJetsetVec(.TRUE.)
!          call PYLIST(2)
!          call GetJetSetVec_List(6,1,N)

          call GetJetsetVecCheckT(-1d-5)

          call GetJetsetVecPYEDIT
       end if

       call GetLeading_PY         ! find leading particles
       call PYEDIT(1)             ! clean up event list

       !...Rotate and boost the whole event to final system:

       phi = atan2(pcm(2),pcm(1))
       theta = atan2(sqrt(pcm(1)**2+pcm(2)**2),pcm(3))
       call PYROBO(1,N, theta,phi, beta(1),beta(2),beta(3))

       if (useJetSetVec) then
          call GetJetsetVecPYROBO(theta,phi, beta(1),beta(2),beta(3))

!          call PYLIST(2)
!          call GetJetSetVec_List(6,1,N)
!          stop
       end if


       !...Copy Particles to ouput-vector

       call SetVectorFromPYJETS(outPart, real(max(1d0,VINT(52))) )

       ! If MSTI(1)==95 it was a low-pT-scattering and VINT(52)=0
       ! Otherwise VINT(52) holds something like
       !   Q2 = pT_hat^2 + (m_3^2+m_4^2)/2
       ! See PYTHIA manulal for details.
       !
       ! 1d0 is set as a (dummy) minimal value !!!!

       !...Correct Charges

!       DeltaQ0 = DeltaQ
       if (DeltaQ.ne.0) call CorrectChargeVector(outPart,DeltaQ)
       if (DeltaQ.ne.0) then
          if (DoPr(2)) then
             write(*,*) 'DoColl_Pythia: Charge correction failed. ReDo Event!!'
!             write(*,*) 'DeltaQ : ',DeltaQ0,DeltaQ
!             call PYLIST(2)
          end if

          cycle
       end if

       !...correct antiparticles if necessary

       if (EventIsAnti) then
!      write(*,*) 'DoColl_Pythia: Event was done as ANTI'
          call ConvertToAnti(outPart)
       end if

       if (N==2) then

          ! Test for elastic event
          if (IsElastic(InPart(1:2),OutPart(1:2))) cycle

          ! Test for charge-exchange event in antinucleon-nucleon,
          ! antinucleon-delta or nucleon-antidelta collision
          if ((InPart(1)%Id+InPart(2)%Id<=nucleon+delta) .and. &
              (InPart(1)%antiParticle.neqv.InPart(2)%antiParticle) .and. &
              IsChargeExchange(InPart(1),InPart(2),OutPart(1),OutPart(2))) cycle

       end if

       !...exit the loop
       exit

    end do


   flagOK = .TRUE.





!   call PYGIVE('MSTI(1)=')

   !  11 f + f' -> f + f' (QCD)
   !  12 f + fbar -> f' + fbar'
   !  13 f + fbar -> g + g
   !  28 f + g -> f + g
   !  53 g + g -> f + fbar
   !  68 g + g -> g + g
   !  95 Low-pT scattering

!   call PYGIVE('MSTI(1)=')
!   call PYGIVE('VINT(52)=')
!   call PYLIST(2)
!   write(*,*) '===================='

  end subroutine DoColl_Pythia



  !****************************************************************************
  !****is* Coll_Pythia/CalcMomentum
  ! NAME
  ! subroutine CalcMomentum(sqrts,mA,mB, EA,EB,pA)
  !
  ! PURPOSE
  ! calculate the momentum in a two particle system (fundamental kinematics)
  !
  ! INPUTS
  ! * sqrts,mA,mB, EA,EB -- cm energy, masses, energies
  !
  ! OUTPUT
  ! * pA -- momentum of particle A
  !
  ! NOTES
  ! This is very elementary. Check code in order not to duplicate routines!
  !****************************************************************************

  subroutine CalcMomentum(sqrts,mA,mB, EA,EB,pA)
    real, intent(in)  :: sqrts,mA,mB
    real, intent(out) :: EA,EB,pA

    real :: mA2,mB2,s

    s = sqrts**2
    mA2 = mA**2
    mB2 = mB**2

    EA = (s+(mA2-mB2))/(2*sqrts)
    EB = (s-(mA2-mB2))/(2*sqrts)
    pA = sqrt((s-mA2-mB2)**2-4*mA2*mB2)/(2*sqrts)

  end subroutine CalcMomentum


 !*****************************************************************************

end module Coll_Pythia



