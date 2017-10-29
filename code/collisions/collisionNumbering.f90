!******************************************************************************
!****m* /collisionNumbering
! NAME
! module collisionNumbering
!
! PURPOSE
! * includes the numbering of the particles according to their collision
!   partners.
! * This numbering is stored in the "%event" variable.
! * Also includes the routine which based on %event decides whether two
!   particles did collide before.
! * Include routines which store and output the number of collisions.
! NOTES
! In GiBUU we label each particle by the event it stems from.
! This labeling is different for real and perturbative particles.
!
! REAL PARTICLES:
! * We assign a unique number to each event.
!   Event #1 gets %event=firstEventNumber,
!   Event #2 gets %event=firstEventNumber+1,...
! * To check whether to particles stem from the same event,
!   we just need to check whether their %event entries coincide.
!
! PERTURBATIVE PARTICLES:
! *  The perturbative particles produced in some event get
!    as %event the %number of the realParticle which took part in the event
! *  To check whether a real and perturbative particle stem from
!    the same event, we check whether pert%event=real%number
!******************************************************************************
module collisionNumbering

  implicit none
  private

  !****************************************************************************
  !****g* collisionNumbering/firstEventNumber
  ! PURPOSE
  ! First number to be given to particles in a collision
  ! SOURCE
  !
  integer, parameter :: firstEventNumber=1000000
  !****************************************************************************

  !****************************************************************************
  !****g* collisionNumbering/CountedEvents
  ! PURPOSE
  ! Global storage of the number of events
  ! * first index  --- 0: integrated number, 1: from a given time step only
  ! * second index --- 0: all, 1: 1-Body, 2: 2-Body, 3: 3-Body
  ! * third index  --- 1: real particles, 2: perturbative particles
  !
  ! SOURCE
  !
  integer, dimension(0:1,0:3,1:2), save :: CountedEvents = 0
  !****************************************************************************

  public :: real_numbering, pert_numbering, reportEventNumber
  public :: pert_firstnumbering, pert_firstnumbering12
  public :: check_justCollided
  public :: getCountedEvents, nullCountedEvents, writeCountedEvents
  public :: getn_participants, nulln_participants
  public :: real_firstnumbering

  integer, save :: n_participants_rr=0
  integer, save :: n_participants_rp=0

contains

  !****************************************************************************

  !****************************************************************************
  !****f* collisionNumbering/getCountedEvents
  ! NAME
  ! integer function getCountedEvents(iStat, iN, iType)
  ! PURPOSE
  ! Returns number of counted events
  !
  ! INPUTS
  ! * integer :: iStat --- 0: integrated number, 1: from a given time step only
  ! * integer :: iN --- 0: all, 1: 1-Body, 2: 2-Body, 3: 3-Body
  ! * integer :: iType  --- 1: real particles, 2: perturbative particles
  !
  ! OUTPUT
  ! * integer
  !
  ! NOTES
  ! no range checking
  !****************************************************************************
  integer function getCountedEvents(iStat, iN, iType)
    integer, intent(in) :: iStat, iN, iType
    getCountedEvents = CountedEvents(iStat, iN, iType)
  end function getCountedEvents


  !****************************************************************************
  !****s* collisionNumbering/nullCountedEvents
  ! NAME
  ! subroutine nullCountedEvents(iflag)
  !
  ! PURPOSE
  ! Sets countedEvents back to zero
  !
  ! INPUTS
  ! * integer:: iflag -- =0: set to zero all countedEvents,
  !   =1: set to zero only the previous time step counter
  !****************************************************************************
  subroutine nullCountedEvents(iflag)
    integer, intent(in) :: iflag
    if (iflag.eq.0) then
      CountedEvents(:,:,:)=0
    else if (iflag.eq.1) then
      CountedEvents(1,:,:)=0
    end if
  end subroutine nullCountedEvents

  !****************************************************************************
  !****s* collisionNumbering/writeCountedEvents
  ! NAME
  ! subroutine writeCountedEvents(iChoice)
  !
  ! PURPOSE
  ! Write some values of the stored CountedEvents to stdout
  !
  ! INPUTS
  ! * integer :: iChoice -- selects the values to be printed
  ! * real, optional :: time --- actual time
  !****************************************************************************
  subroutine writeCountedEvents(iChoice, time)
    integer, intent(in) :: iChoice
    real, intent(in), optional :: time

    character(200), parameter :: format4 = '(A," -- all, 1-body, 2-body, 3-body:",4i9)'
    character(200), parameter :: format5a= '("number of ",A," events at time=",f12.3,": ",4i9)'
    character(200), parameter :: format5b= '("number of ",A," events (1-,2-,3-Body)      : ",4i9)'

    select case (iChoice)
    case (0)
       write(*,format4) 'Real        ',CountedEvents(0,0:3,1)
       write(*,format4) 'Perturbative',CountedEvents(0,0:3,2)
    case (1)
       write(*,format5a) 'REAL        ',time,CountedEvents(1,1:3,1)
       write(*,format5a) 'PERTURBATIVE',time,CountedEvents(1,1:3,2)
    case (2)
       write(*,format5b) 'REAL        ',CountedEvents(0,1:3,1)
       write(*,format5b) 'PERTURBATIVE',CountedEvents(0,1:3,2)
    end select


  end subroutine writeCountedEvents

  !****************************************************************************
  !****f* collisionNumbering/getn_participants
  ! NAME
  ! function getn_participants()
  !
  ! PURPOSE
  ! Returns number of participants.
  !
  ! RESULT
  ! * integer, dimension(1:2) :: getn_participants
  !****************************************************************************
  function getn_participants()
    integer, dimension(1:2) :: getn_participants
    getn_participants(1)=n_participants_rr
    getn_participants(2)=n_participants_rp
  end function getn_participants


  !****************************************************************************
  !****s* collisionNumbering/nulln_participants
  ! NAME
  ! subroutine nulln_participants
  !
  ! PURPOSE
  ! Sets to zero the number of participants.
  !****************************************************************************
  subroutine nulln_participants
    n_participants_rr=0
    n_participants_rp=0
  end subroutine nulln_participants


  !****************************************************************************
  !****s* collisionNumbering/ReportEventNumber
  ! NAME
  ! subroutine ReportEventNumber(InPart,OutPart, nr,time,code1,code2,weight)
  ! PURPOSE
  ! Used for statistics. Called e.g. by collsionTerm.
  !
  ! INPUTS
  ! * type(particle),dimension(:)   :: InPart -- incoming particles
  ! * type(particle),dimension(:)   :: OutPart-- outgoing particles
  ! * integer,dimension(1:2)        :: nr     -- EventNumber
  ! * real                          :: time   -- time of collision
  ! * integer                       :: code1  -- code of collision
  ! * integer,             optional :: code2  -- subcode of collision
  ! * real,                optional :: weight -- event-weight
  !
  ! possible values for code1 are:
  ! * 1-body:    11= real particle decay, 12= pert.particle decay
  ! * 2-body:   211= real+real ,         212= real+pert
  ! * 3-body:  3111= real+real+real,    3112= real+real+pert
  !
  ! possible values for code2 (2-body) are:
  ! *  -3: BaB    (Baryon-Antibaryon-Annihilation)
  ! *  -2: Manni  (Meson-Baryon-Annihilation)
  ! *  -1: Elastic
  ! *   1: FRITIOF
  ! *   2: PYTHIA
  !
  ! OUTPUT
  ! ---
  !
  ! NOTES
  ! You can (ab)use this routine for debugging purposes:
  ! listing the event if a particle with a given ID collides or
  ! if a particle, which is produced in some given event collides again,
  ! many possibilities...
  !
  ! But be aware: you have no information about finalState !
  !****************************************************************************
  subroutine ReportEventNumber(InPart,OutPart, nr,time,code1,code2,weight)

    use particleDefinition
    use collisionReporter, only: cR_Add
!     use output, only: WriteParticle

    type(particle),dimension(:),  intent(in) :: InPart
    type(particle),dimension(:),  intent(in) :: OutPart
    integer,dimension(1:2),       intent(in) :: nr
    real,                         intent(in) :: time
    integer,                      intent(in) :: code1
    integer,             optional,intent(in) :: code2
    real,                optional,intent(in) :: weight


!...local variables:
    integer :: code11 ! , eventclass
    real :: w;
!     integer :: i
    real, save :: time0_rr=-100., time0_rp=-100.
    integer, save :: number0
!     integer :: nOut

!...set some defaults:
    code11=0
    if (present(code2)) code11=code2
    w = 1.0
    if (present(weight)) w = weight

!    eventclass = -1
!    if (code1.lt.100) then
!       eventclass = 1
!    else if(code1.lt.300) then
!       eventclass = 2
!    else if(code1.ge.3000) then
!       eventclass = 3
!    else
!       write(*,*) 'Ooops! eventclass of ',code1,' ???'
!    endif

!...

    select case (code1)
    case (211)
       if (time.ne.time0_rr) then
          number0=getNumber()
          n_participants_rr=0
          time0_rr=time
       end if
       ! Count both colliding real particles, if they were created
       ! before the given time step:
       if (InPart(1)%number.le.number0) n_participants_rr=n_participants_rr+1
       if (InPart(2)%number.le.number0) n_participants_rr=n_participants_rr+1
    case (212)
       if (time.ne.time0_rp) then
          number0=getNumber()
          n_participants_rp=0
          time0_rp=time
       end if
       ! Count only perturbative particles:
       if (InPart(2)%number.le.number0) n_participants_rp=n_participants_rp+1
    end select

    select case (code1)
    case (11)     ! real particle decay:
       CountedEvents(:,1,1) = CountedEvents(:,1,1)+1
       CountedEvents(:,0,1) = CountedEvents(:,0,1)+1
    case (12)      ! perturbative particle decay:
       CountedEvents(:,1,2) = CountedEvents(:,1,2)+1
       CountedEvents(:,0,2) = CountedEvents(:,0,2)+1
    case (211)      ! real+real collision:
       CountedEvents(:,2,1) = CountedEvents(:,2,1)+1
       CountedEvents(:,0,1) = CountedEvents(:,0,1)+1
    case (212)      ! real+perturbative collision:
       CountedEvents(:,2,2) = CountedEvents(:,2,2)+1
       CountedEvents(:,0,2) = CountedEvents(:,0,2)+1
    case (3111)      ! real+real+real collision:
       CountedEvents(:,3,1) = CountedEvents(:,3,1)+1
       CountedEvents(:,0,1) = CountedEvents(:,0,1)+1
    case (3112)      ! real+real+perturbative collision:
       CountedEvents(:,3,2) = CountedEvents(:,3,2)+1
       CountedEvents(:,0,2) = CountedEvents(:,0,2)+1
    case default
       write(*,*) ' In ReportEventNumber: wrong code1', code1
       stop
    end select

!    write(*,'(A,i5,i3,f9.3,f9.3,3i9)') 'code,N :',code1,code11,time,sqrtS(InPart),InPart%number

    if (code1==211 .or. code1==212) call cR_Add (sqrtS(InPart), time, code11, w)


!    call WriteEventHistory


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! You can (ab)use this routine for debugging purposes:
!!!!! listing the event if a particle with a given ID collides or
!!!!! if a particle, which is produced in some given event collides again,
!!!!! many possibilities...
!!!!! Be aware: you have no information about finalState !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!$    do i=1,size(InPart,dim=1)
!!$       if ((InPart(i)%firstevent.eq.628001).or.(InPart(i)%firstevent.eq.628002)) then
!!$          call WriteEventHistory
!!$          exit
!!$       end if
!!$    end do

!    if (time.eq.13.) then
!       write(*,'(A,i5,i3,f9.3,f9.3,3i9)') 'code,N :',code1,code11,time,sqrtS(InPart),InPart%number
!       do i=1,size(InPart,dim=1)
!          call WriteParticle(6,99,i,InPart(i))
!       end do
!!$
!!$
!!$       if (code11.eq.2) then
!!$          call QYLIST(2)
!!$       else if (code11.eq.1) then
!!$          call LULIST(2)
!!$       endif
!!$
!    endif


!!$    do i=1,size(InPart,dim=1)
!!$       if (InPart(i)%number .eq. 100049) then
!!$          write(*,*) 'part ',InPart(i)%number,' just collided!!!!'
!!$          write(*,*) nr,time,code1,code11
!!$
!!$          call QYLIST(2)
!!$          call GetJetSetVec_List(6,0,0)
!!$
!!$          stop
!!$       endif
!!$   enddo

!   contains
!
!     subroutine WriteEventHistory
!
!       use particleProperties, only : partName
!
!       integer, save :: iEvent = 0
!       integer, parameter :: iFile = 133
!       character*15 :: NameIn(3),NameOut(30)
!
!       iEvent = iEvent +1
!
!       do i=size(OutPart),1,-1
!          nOut = i
!          if (OutPart(i)%ID .gt. 0) exit
!       end do
!
!       if (code1.lt.100) then
!          write(iFile,1001) iEvent,time,nr,code1,sqrtS(InPart), &
!               & InPart(1)%number, OutPart(1:nOut)%number
!       else if(code1.lt.1000) then
!          write(iFile,1002) iEvent,time,nr,code1,code11,sqrtS(InPart), &
!               & InPart(1:2)%number, OutPart(1:nOut)%number
!       else if(code1.lt.10000) then
!          write(iFile,1003) iEvent,time,nr,code1,sqrtS(InPart), &
!               & InPart(1:3)%number, OutPart(1:nOut)%number
!       endif
!
!       do i=1,nOut
!          NameOut(i) = PartName(OutPart(i))
!       end do
!       if (code1.lt.100) then
!          do i=1,1
!             NameIn(i) = PartName(InPart(i))
!          enddo
!          write(iFile,2001) iEvent,NameIn(1),NameOut(1:nOut)
!       else if(code1.lt.1000) then
!          do i=1,2
!             NameIn(i) = PartName(InPart(i))
!          enddo
!          write(iFile,2002) iEvent,NameIn(1:2),NameOut(1:nOut)
!       else if(code1.lt.10000) then
!          do i=1,3
!             NameIn(i) = PartName(InPart(i))
!          enddo
!          write(iFile,2003) iEvent,NameIn(1:3),NameOut(1:nOut)
!       end if
!
! 1001  format(i9,f5.2,'[',2i9,']:',i5,'    |',f9.3,'|',i9,' --> ',30i9)
! 1002  format(i9,f5.2,'[',2i9,']:',i5,'(',i2,')|',f9.3,'|',2i9,' --> ',30i9)
! 1003  format(i9,f5.2,'[',2i9,']:',i5,'    |',f9.3,'|',3i9,' --> ',30i9)
!
! 2001  format(i9,'   ',A15,' --> ',30A15)
! 2002  format(i9,'   ',2A15,' --> ',30A15)
! 2003  format(i9,'   ',3A15,' --> ',30A15)
!
!     end subroutine WriteEventHistory


  end subroutine ReportEventNumber


  !****************************************************************************
  !****f* collisionNumbering/pert_numbering
  ! NAME
  ! integer function pert_numbering(realteilchen)
  !
  ! PURPOSE
  ! Calculate the Event number of a perturbative particle which collided with
  ! the real particle "realteilchen".
  ! If realTeilchen is not given than the return value is -999,
  ! which can never be associated to any particle.
  !
  ! INPUTS
  ! * type(particle), optional :: realTeilchen
  !
  ! OUTPUT
  ! event number
  !****************************************************************************
  integer function pert_numbering(realteilchen)
    use particleDefinition
    type(particle),intent(in),optional :: realTeilchen
    if (present(realTeilchen)) then
       pert_numbering=realTeilchen%number
    else
       ! particle shall not be associated with any real particle
       pert_numbering=-999
    end if
  end function pert_numbering


  !****************************************************************************
  !****f* collisionNumbering/real_numbering
  ! NAME
  ! integer function real_numbering()
  !
  ! PURPOSE
  ! Calculate the Event number of a real particle.
  ! The event number for the real particles are unique.
  ! A collision happened at later time will cause a higher event number.
  !
  ! INPUTS
  ! ---
  !
  ! OUTPUT
  ! event number
  !
  ! NOTES
  ! The minimal value of numbering is firstEventNumber.
  !****************************************************************************
  integer function real_numbering()

    !**************************************************************************
    !****g* real_numbering/eventNumber
    ! SOURCE
    !
    integer, save  :: eventNumber=firstEventNumber
    !
    ! PURPOSE
    ! Number to be given to every real particle stemming from a collision
    ! event:
    ! * This number should by unique !
    ! * Therefore we raise  "eventnumber" by one for every event.
    ! * All produced particles then get this "eventnumber" set
    !   into "teilchen(i,j)%event".
    ! If two real particles agree in "teilchen(i,j)%event",
    ! "teilchen(k,l)%event" then they are not allowed to scatter to prevent
    ! multiple collisions.
    ! of one pair without interaction with the outside world.
    !**************************************************************************
    real_numbering=eventNumber
    eventNumber=eventNumber+1
  end function real_numbering

  !****************************************************************************
  !****f* collisionNumbering/pert_firstnumbering
  ! NAME
  ! integer function pert_numbering(realTeilchen,pertTeilchen)
  !
  ! PURPOSE
  ! calculate the value of firstevent, which has to be given to the
  ! final particles, when realTeilchen and pertTeilchen collided
  !
  ! possible return values:
  ! * realTeilchen%number, if pertTeilchen not given
  ! * pertTeilchen%firstEvent, if pertTeilchen%firstEvent.ne.0
  ! * some exotic value, if you wish
  ! * realTeilchen%number as default
  !
  ! INPUTS
  ! * type(particle)          :: realTeilchen
  ! * type(particle),optional :: pertTeilchen
  !
  ! NOTES
  ! This routine should be called only, if for a incoming particle we
  ! had %firstEvent==0;
  ! If not, it does also the right stuff, but slower.
  !
  ! if only realTeilchen is given, it is as calling
  ! "pert_numbering(realTeilchen)": return value =  realTeilchen%number
  !
  !****************************************************************************
  integer function pert_firstnumbering(realTeilchen,pertTeilchen)
    use particleDefinition
    use inputGeneral, only: eventType

    type(particle),intent(in)          :: realTeilchen
    type(particle),intent(in),optional :: pertTeilchen

    if (.not.present(pertTeilchen)) then
       pert_firstnumbering = realTeilchen%number
       return
    end if

    if (pertTeilchen%firstEvent.ne.0) then
       pert_firstnumbering = pertTeilchen%firstEvent
       return
    end if

    select case (eventType)
    case (12) ! --- HiPion
       pert_firstnumbering = pert_firstnumbering12()

!       write(*,'(A,4i10)') 'pert_firstnumbering',realTeilchen%number,pertTeilchen%number,pert_firstnumbering


    case default
       pert_firstnumbering = realTeilchen%number

    end select

  end function pert_firstnumbering


  !****************************************************************************
  !****f* collisionNumbering/pert_firstnumbering12
  ! NAME
  ! integer function pert_firstnumbering12(reset,DoNotInc)
  !
  ! PURPOSE
  ! just return an increasing number every call
  !
  ! INPUTS
  ! * logical, optional :: reset -- if .true. the internal counter is
  !   reset to zero
  ! * logical, optional :: DoNotInc -- if .true. the internal counter is
  !   not increased, the returned number is the number given before
  !
  ! NOTES
  ! * called by pert_numbering
  ! * used for eventType == 12 (HiPion)
  !****************************************************************************
  integer function pert_firstnumbering12(reset,DoNotInc)
    logical, optional :: reset
    logical, optional :: DoNotInc

    integer, save :: iCount = 0

    iCount = iCount + 1
    if (present(DoNotInc)) then
       if (DoNotInc) iCount = iCount - 1
    end if

    if (present(reset)) then
       if (reset) iCount=0
    end if

    pert_firstnumbering12 = iCount
    return
  end function pert_firstnumbering12


  !****************************************************************************
  !****f* collisionNumbering/real_firstnumbering
  ! NAME
  ! integer function real_firstnumbering(reset)
  !
  ! PURPOSE
  ! just return an increasing number every call
  !
  ! INPUTS
  ! * logical, optional :: reset -- if .true. the internal counter is
  !   reset to zero
  !****************************************************************************
  integer function real_firstnumbering(reset)
    logical, optional :: reset

    integer, save :: iCount = 0

    iCount = iCount + 1

    if (present(reset)) then
       if (reset) iCount=0
    end if

    real_firstnumbering = iCount
  end function real_firstnumbering


  !****************************************************************************
  !****f* collisionNumbering/check_justCollided
  ! NAME
  ! logical function check_justCollided(a,b)
  ! PURPOSE
  ! Checks wether particles a and b stem from the same collision,
  ! and no other collision happened yet.
  ! RESULT
  ! * true = if they are just created before by the very same event
  ! * false= if not...
  !****************************************************************************
  logical function check_justCollided(a,b)
    use particleDefinition
    type(particle), intent(in) :: a,b

    if (a%perturbative.neqv.b%perturbative) then !real-pert

       if (a%perturbative) then ! a is perturbative, b is real
         check_justCollided = (a%event(1)==b%number) .or. (a%event(2)==b%number)
       else ! b is perturbative, a is real
         check_justCollided = (b%event(1)==a%number) .or. (b%event(2)==a%number)
       end if

    else  ! a and b are both real or both perturbative
       check_justCollided = (a%event(1)==b%event(1)) .or. (a%event(1)==b%event(2)) .or. &
                            (a%event(2)==b%event(1)) .or. (a%event(2)==b%event(2))
    end if

  end function check_justCollided


end module collisionNumbering
