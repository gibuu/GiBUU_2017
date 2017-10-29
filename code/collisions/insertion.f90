!******************************************************************************
!****m* /Insertion
! NAME
! module Insertion
! PURPOSE
! This module collects routines for inserting particles into the particle vectors.
!******************************************************************************
module Insertion

  implicit none
  private

  !****************************************************************************
  !****g* Insertion/minimumEnergy
  ! SOURCE
  !
  real, save, public :: minimumEnergy=0.005
  ! PURPOSE
  ! Minimal kinetic energy in GeV for produced perturbative nucleons.
  ! If their energy is below this threshold, then they are not propagated,
  ! i.e. they are not inserted in the particle vector.
  !
  ! NOTES
  ! This value was formerly given in the namelist "collisionterm".
  !****************************************************************************


  !****************************************************************************
  !****g* Insertion/propagateNoPhoton
  ! SOURCE
  !
  logical, save :: propagateNoPhoton=.true.
  ! PURPOSE
  ! If .true. then we eliminate all photons, such that they are not propagated and
  ! do not show up in the particle vector.
  ! If .false. then photons are explicitly propagated.
  !****************************************************************************


  public :: particlePropagated,GarbageCollection,FindLastUsed,setIntoVector

  logical, save :: initFlag=.true.

contains

  !****************************************************************************
  !****s*  Insertion/readInput
  ! NAME
  ! subroutine Insertion_SetInitialized
  ! PURPOSE
  ! just set initflag to .true.
  !****************************************************************************
  subroutine readInput
    use output

    integer :: ios

    !**************************************************************************
    !****n* Insertion/insertion
    ! NAME
    ! NAMELIST /insertion/
    ! PURPOSE
    ! Namelist for module Insertion includes:
    ! * minimumEnergy
    ! * propagateNoPhoton
    !**************************************************************************
    NAMELIST /insertion/ minimumEnergy,propagateNoPhoton

    rewind(5)
    call Write_ReadingInput('insertion',0)
    read(5,nml=insertion,iostat=ios)
    call Write_ReadingInput('insertion',0,ios)

    write(*,*) 'Minimal Energy for perturbative nucleons =', minimumEnergy
    write(*,*) 'Propagate No Photons?                       =', propagateNoPhoton
    call Write_ReadingInput('insertion',1)

    initFlag = .false.

  end subroutine readInput


  !****************************************************************************
  !****s* Insertion/GarbageCollection
  ! NAME
  ! subroutine GarbageCollection(partVec,DoCollHist)
  !
  ! PURPOSE
  ! Rearrange particles in the vector (per ensemble) in such a way,
  ! that there are no holes (i.e. entries with the special ID 'NOP') inbetween.
  ! In addition, all particles after the last one get the special ID 'EOV' in
  ! order to indicate, that no non-empty entries will follow.
  !
  ! INPUTS
  ! * type(particle), dimension(:,:) :: partVec
  ! * logical, OPTIONAL :: DoCollHist -- Flag whether to do additional
  !   rearrangements
  !
  ! OUTPUT
  ! * partVec changed
  !****************************************************************************
  subroutine GarbageCollection(partVec,DoCollHist,iiEns)

    use IdTable, only: EOV, NOP
    use particleDefinition
    use CollHistory, only: CollHist_DoGBC, CollHist_SetSize

    type(particle), dimension(:,:), intent(inOUT) ::  partVec
    logical, OPTIONAL, intent(in) :: DoCollHist
    integer, OPTIONAL, intent(in) :: iiEns

    integer :: i1,i2,iEns,nPart,iEns1,iEns2
    logical :: DoCH

    nPart = size(partVec,dim=2)

    DoCH = .false.
    if (present(DoCollHist)) DoCH = DoCollHist

    if (DoCH) call CollHist_SetSize(ubound(partVec))


    if (present(iiEns)) then
       iEns1 = iiEns
       iEns2 = iiEns
    else
       iEns1 = 1
       iEns2 = size(partVec,dim=1)
    end if

    EnsLoop: do iEns=iEns1,iEns2

       i1 = 1
       i2 = nPart
       Loop0: do

          ! find 'last particle':
          Loop2: do
             if (i1>=i2) exit Loop0 ! finished !!!
             if (partVec(iEns,i2)%ID > 0) exit Loop2 ! particle found
             if (partVec(iEns,i2)%ID == NOP) partVec(iEns,i2)%ID = EOV
             i2 = i2-1
          end do Loop2

          ! find 'first hole':
          Loop1: do
             if (i1>=i2) exit Loop0 ! finished !!!
             if (partVec(iEns,i1)%ID < 0) exit Loop0 ! finished !!!
             if (partVec(iEns,i1)%ID == NOP) exit Loop1 ! hole found
             i1 = i1+1
          end do Loop1

          ! copy 'last particle' into 'first hole':
          partVec(iEns,i1) = partVec(iEns,i2)
          partVec(iEns,i2)%ID = EOV

          if (DoCH) call CollHist_DoGBC(iEns,i1,i2)

          i2 = i2-1

       end do Loop0
    end do EnsLoop
  end subroutine GarbageCollection


  !****************************************************************************
  !****s* Insertion/FindLastUsed
  ! NAME
  ! function FindLastUsed(particles)
  !
  ! PURPOSE
  ! Returns the last (used) entry from particle vector particles.
  ! The particle with the next index has the special ID 'EOV'.
  !
  ! INPUTS
  ! * type(particle), dimension(:) :: particles -- particle vector
  !
  ! OUTPUT
  ! * integer :: FindLastUsed -- number of last used entry
  !
  ! NOTES
  ! Uses a (fast) bisection method!
  ! Therefore it relies on that "GarbageCollection" has been called.
  !****************************************************************************
  function FindLastUsed(particles)

    use particleDefinition

    integer :: FindLastUsed
    type(particle), dimension(:), intent(in) :: particles

    integer :: i1,i2,im

    i1 = 1
    i2 = size(particles)

    if (particles(i2)%ID>=0) then
       FindLastUsed = i2
       return
    end if

    do while(i2-i1 > 1)
       im = (i1+i2)/2
       if (particles(im)%ID>=0) then
          i1 = im
       else
          i2 = im
       end if
    end do
    FindLastUsed = i1
    return

  end function FindLastUsed


  !****************************************************************************
  !****s*  Insertion/setIntoVector
  ! NAME
  ! subroutine setIntoVector(finalState, partVec, flagOK,numberIsSet,numbers,positions)
  !
  ! PURPOSE
  ! This subroutine tries to find empty spaces in "partVec" and
  ! sets the elements of "finalState" into these holes.
  !
  ! INPUTS
  ! * type(particle),dimension(:,:) :: partVec
  ! * type(particle),dimension(:)   :: finalState
  ! * logical, optional             :: numberIsSet --
  !   .true. if finalstate has already %number set,
  !   .false. if this still needs to be done
  !
  ! OUTPUT
  ! * logical                       :: flagOK --
  !   true, if all Insertions were successful
  ! * integer,dimension(:),optional :: numbers --
  !   Vector with value of %number assigned to each final state entry
  ! * integer,dimension(2,:),optional :: positions --
  !   vector with positions, where particles were inserted (iEns,iPart)
  ! * partVec changed
  !
  ! NOTES
  ! Only particles for which the function "particlePropagated" returns true
  ! are really considered.
  !
  ! Concerning the variables "lastEnsemble, lastIndex" :
  ! By saving the index of the last hole in the vector we try to save time
  ! when searching for the next hole.
  !****************************************************************************
  subroutine setIntoVector(finalState, partVec, flagOK, numberIsSet, numbers, positions)

    use particleDefinition

    type(particle), dimension(:,:), intent(inOUT) ::  partVec
    type(particle), dimension(:),   intent(in)    ::  finalState
    logical, intent(out) :: flagOK
    integer, dimension(:),   intent(out), optional::  numbers
    integer, dimension(:,:), intent(out), optional::  positions
    logical, optional, intent(in) :: numberIsSet

    integer, save :: lastEnsemble=0
    integer, save :: lastIndex=0
    integer :: i, j,k, iEns, nEns, nPart

    ! Initialize at first call
    if (initFlag) call ReadInput

    nEns  = size(partVec,dim=1)
    nPart = size(partVec,dim=2)

    if (present(positions)) positions = 0

    ! Loop over final state particles :
    finalState_Loop :do k=lbound(finalState,dim=1),ubound(finalState,dim=1)

       flagOK=.false.

       ! (1) Decision if particle is propagated
       if (.not.particlePropagated(finalState(k))) then
!          if (finalState(k)%ID > 0) write(*,*) 'skip',finalState(k)%ID
          flagOK=.true.
          cycle
       end if

       ! (2) Find empty space in particle vector

       ensemble_Loop : do iEns=1,nEns
          i=Mod(iEns+lastEnsemble-1,nEns)+1

          ! (2a) First try: get the next from the last call

          j = lastIndex+1
          if (j > 1 .and. j <= nPart) then
             if (partVec(i,j-1)%Id>=0.and.partVec(i,j)%Id<0) then ! all is fine
                partVec(i,j)=finalState(k)
                if (present(numberIsSet)) then
                   if (.not.numberIsSet) call setNumber(partVec(i,j))
                else
                   call setNumber(partVec(i,j))
                end if
                lastEnsemble=i
                lastIndex   =j
                flagOK=.true.
                if (present(numbers)) numbers(k)=partVec(i,j)%number
                if (present(positions)) positions(1:2,k) = (/i,j/)
                exit ensemble_Loop
             end if
          end if

          ! (2b) Second Try: get the next from the Search

          j = FindLastUsed(partVec(i,:))+1
          if (j > 1 .and. j <= nPart) then
             if (partVec(i,j-1)%Id>=0.and.partVec(i,j)%Id<0) then ! all is fine
                partVec(i,j)=finalState(k)
                if (present(numberIsSet)) then
                   if (.not.numberIsSet) call setNumber(partVec(i,j))
                else
                   call setNumber(partVec(i,j))
                end if
                lastEnsemble=i
                lastIndex   =j
                flagOK=.true.
                if (present(numbers)) numbers(k)=partVec(i,j)%number
                if (present(positions)) positions(1:2,k) = (/i,j/)
                exit ensemble_Loop
             end if
          end if


          ! (2c) Third Try: Go through the vector step by step...

!          write(*,*) 'setIntoVector, step 2c reached! ',k
          do j=1,nPart
             if (partVec(i,j)%Id<=0) then
!                write(*,*) '*** ID<=0, j=',j, partVec(i,j)%Id

                partVec(i,j)=finalState(k)

                if (present(numberIsSet)) then
                   if (.not.numberIsSet) call setNumber(partVec(i,j))
                else
                   call setNumber(partVec(i,j))
                end if
                lastEnsemble=i
                lastIndex   =j
                flagOK=.true.
                if (present(numbers)) numbers(k)=partVec(i,j)%number
                if (present(positions)) positions(1:2,k) = (/i,j/)
                exit ensemble_Loop

             end if
          end do

          ! (3) No Hole found!

          write(*,*) 'setIntoVector, step 3 reached: no Hole found! ',k,iEns

       end do ensemble_Loop

       if (.not.flagOK) exit finalState_Loop     ! no hole could be found
    end do finalState_Loop

  end subroutine setIntoVector


  !****************************************************************************
  !****f*  Insertion/particlePropagated
  ! NAME
  ! logical function particlePropagated (Part)
  !
  ! PURPOSE
  ! Return .true. if "Part" is a particle which shall be propagated in
  ! the code: We do not propagate photons and very low-energetic
  ! perturbative nucleons (cf. minimumEnergy).
  !
  ! INPUTS
  ! * type(particle) :: Part
  !
  ! OUTPUT
  ! * function value
  !****************************************************************************
  logical function particlePropagated (Part)

    use particleDefinition
    use IDTable, only: photon, nucleon, isLepton

    type(particle), intent(in) :: Part

    particlePropagated=.false.

    ! Don't propagate photons and very low energetic perturbative nucleons.

    if (Part%ID <= 0) return
    if (isLepton(part%ID)) return
    if (Part%ID == photon .and. propagateNoPhoton) return

    if (Part%ID == nucleon .and. Part%perturbative .and. ((freeEnergy(Part)-Part%mass) < minimumEnergy)) return

    particlePropagated=.true.

  end function particlePropagated


end module Insertion
