!******************************************************************************
!****m* /AddDecay
! NAME
! Module AddDecay
! PURPOSE
! Perform additional decays, e.g. of Kaons.
!******************************************************************************
Module AddDecay

  implicit none
  private

  public :: PerformAddDecay

contains


  !****************************************************************************
  !****s* AddDecay/PerformAddDecay
  ! NAME
  ! subroutine PerformAddDecay(PartVec,time)
  !
  ! PURPOSE
  ! Loop over the whole particle vector and perform some (additional)
  ! decays, as e.g. kaon to pions.
  ! This is intended to be used with some "reconstruction" routines
  ! as used by experimentalists.
  !
  ! Here we use the JETSET part of PYTHIA, which is totally independend
  ! of any fragmentation stuff etc. This means for particles decaying by
  ! this routine, branching ratios and angular distributions are as given
  ! by PYTHIA parameters, cf. "PYLIST12.txt".
  !
  ! please note: we use PYTHIA, i.e. PYTHIAv6.2
  !
  ! INPUTS
  ! * type(particle),dimension(:,:) :: PartVec -- particle vector
  ! * real                          :: time    -- actual time in the code
  !
  ! OUTPUT
  ! * type(particle),dimension(:,:) :: PartVec -- particle vector
  !
  ! NOTES
  ! At the moment we are just decaying kaon and kaonBar (ID=110,111).
  ! It can be easily expanded to all other particles, if one sets
  ! the corresponding "stability flag" in array MDCY correctly.
  !****************************************************************************
  subroutine PerformAddDecay(PartVec,time)
    use particleDefinition
    use insertion, only: particlePropagated, setIntoVector
    use IdTable, only: isHadron
    use particleProperties, only: hadron
    use CollTools, only: SetSomeDefaults_PY, SetVectorFromPYJETS
    use inputGeneral, only: fullEnsemble

    type(particle),intent(INOUT),dimension(:,:) :: PartVec
    real, intent(in) :: time

    integer :: iEns,iPart, ID
    type(particle),dimension(10) :: outPart
    logical :: setFlag, NumbersAlreadySet

!!$    write(*,*) '...PerformAddDecay'

    call SetSomeDefaults_PY()

    do iEns = 1,size(PartVec,dim=1)
       do iPart = 1,size(PartVec,dim=2)
          ID = PartVec(iEns,iPart)%ID
          if (ID < 0) exit
          if (.not. isHadron(ID)) cycle
          if (iand(hadron(ID)%stability,4) .eq. 0) cycle

          call DoDecay(PartVec(iEns,iPart))

          NumbersAlreadySet = AcceptGuessedNumbers()
          if (particlePropagated(outPart(1))) then
             PartVec(iEns,iPart) = outPart(1)
          else
             PartVec(iEns,iPart)%ID = 0
          end if
          if (fullensemble) then
             call setIntoVector(outPart(2:), PartVec, &
                  setFlag,NumbersAlreadySet)
          else
             call setIntoVector(outPart(2:), PartVec(iEns:iEns,:), &
                  setFlag,NumbersAlreadySet)
          end if

       end do
    end do


  contains


    subroutine DoDecay(Part)
      use hadronFormation, only: useJetSetVec
      use history, only: setHistory
      use ID_translation, only: KFfromBUU
!       use particleProperties, only: partName

      type(particle) :: Part

      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      integer MSTU,MSTJ
      double precision PARU,PARJ
      SAVE /PYDAT1/

      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      integer MDCY,MDME,KFDP
      double precision BRAT
      SAVE /PYDAT3/

      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      integer MSTP,MSTI
      double precision PARP,PARI
      SAVE /PYPARS/

      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      integer N,NPAD,K
      double precision P,V
      SAVE /PYJETS/

      COMMON/PYINT1/MINT(400),VINT(400)
      integer MINT
      double precision VINT
      SAVE /PYINT1/

      integer :: i
      logical :: dummy

      !...reset PYTHIA arrays:

      N = 1
      K(1,1:5) = 0
      P(1,1:5) = 0.0
      V(1,1:5) = 0.0

      MINT(51) = 0     ! reset error flag

      !...set particle into PYTHIA array:

      K(1,1)   = 1
      K(1,2)   = KFfromBUU (Part)

      P(1,5)   = Part%mass
      P(1,4)   = Part%momentum(0)
      P(1,1:3) = Part%momentum(1:3)

      !...set PYTHIA parameters for decays:

      MSTJ(21) = 2              ! particle decay on/off
!!$      MDCY(113,1)=1             ! decay: K0
!!$      MDCY(112,1)=1             ! decay: K_S0
!!$!      MDCY(105,1)=1             ! decay: K_L0
!!$      MDCY(116,1)=1             ! decay: K+

      !...perform PYTHIA decay:

      if (useJetSetVec) call GetJetsetVecINIT
      call PYEXEC

      !...clean up PYTHIA output:

      call PYEDIT(1)
      K(1:N,4) = 0 ! normally nLead
      K(1:N,5) = 0 ! normally number of string

      !...move PYTHIA output to a GiBUU particle vector:

      outPart%ID = 0

      dummy = useJetSetVec
      useJetSetVec = .false. ! in order to avoid possible problems;
                             ! now we can do this, but not before PYEXEC
      call SetVectorFromPYJETS(outPart,0.0)
      useJetSetVec = dummy

      do i=1,N
         outPart(i)%position = Part%position
         outPart(i)%event = Part%event
      end do

      outPart(1:N)%scaleCS=1.
      outPart(1:N)%in_formation=.false.
      outPart(1:N)%firstevent = Part%firstevent
      outPart(1:N)%perturbative = Part%perturbative
      outPart(1:N)%lastCollisionTime = time
      outPart(1:N)%perWeight = Part%perWeight
      call setHistory (Part, outPart(1:N))

!       write (*,'(20A)') "AddDecay ",partName(Part),(partName(outPart(i)),i=1,N)

!!$      write(*,*)
!!$      call  WriteParticle(6,0,0,Part)
!!$      do i=1,N
!!$         if (outPart(i)%ID.ne.0) call  WriteParticle(6,99,i,outPart(i))
!!$      end do

       if (MINT(51)==2) then
          N = 0
          ! print *,"failure in AddDecay!"
          return ! -> FAILURE
       end if

      !...reset PYTHIA parameters:

      MSTJ(21) = 0              ! particle decay on/off
!!$      MDCY(113,1)=0             ! decay: K0
!!$      MDCY(112,1)=0             ! decay: K_S0
!!$      MDCY(105,1)=0             ! decay: K_L0
!!$      MDCY(116,1)=0             ! decay: K+

    end subroutine DoDecay

  end subroutine PerformAddDecay


end Module AddDecay
