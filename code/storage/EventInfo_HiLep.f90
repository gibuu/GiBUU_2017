!******************************************************************************
!****m* /EventInfo_HiLep
! NAME
! module EventInfo_HiLep
!
! PURPOSE
! This module provides everything to store informations concerning
! reactions (not particles, cf. "module PartInfoList_nLead" etc.
!
! Here we imply:
! The information about iEnsemble and iPart is impressed via
!    ...%firstEvent = iEns*1000+iPart
!
! NOTES
! here we have a somehow more special situation, since the mudule-global
! variable "beamEnergy" is used in all calls.
!
!******************************************************************************
module EventInfo_HiLep

  private

  !****************************************************************************
  !****t* EventInfo_HiLep/tEventInfo_HiLep
  ! SOURCE
  !
  type tEventInfo_HiLep
     logical :: flagOK=.False.
     real    :: Weight=0.
     real    :: nu=0.
     real    :: Q2=0.
     real    :: eps=0.
     integer :: EventType=0     ! cf. module Electron_origin
     real    :: phi_Lepton=0.   ! for lepton induced events only
     real    :: W               ! the value including fermi mom and potentials
     real, dimension(3) :: pos
  end type tEventInfo_HiLep
  !
  ! PURPOSE
  ! This holds all information we want to store connected to
  ! HiLepton and HiPhoton induced reactions.
  !****************************************************************************

  type(tEventInfo_HiLep),save,dimension(:,:),allocatable :: ARR_HiLep

  !****************************************************************************
  !****g* EventInfo_HiLep/beamEnergy
  ! SOURCE
  !
  real,save :: beamEnergy = 0.
  !
  ! PURPOSE
  ! This variable stores the electron beam energy
  !****************************************************************************

  public :: EventInfo_HiLep_Init,EventInfo_HiLep_Store,EventInfo_HiLep_Get
!   PUBLIC :: EventInfo_HiLep_GetBeamEnergy
  public :: EventInfo_HiLep_Dump

contains

  !****************************************************************************
  !****s* EventInfo_HiLep/EventInfo_HiLep_Init
  ! NAME
  ! subroutine EventInfo_HiLep_Init(NumEnsembles,MassNum)
  !
  ! PURPOSE
  ! if necessary allocate memory. Reset the corresponding arrays.
  !
  ! INPUTS
  ! * integer :: NumEnsembles -- number of ensembles
  ! * integer :: MassNum      -- mass number of target (="A")
  !
  ! OUTPUT
  ! ---
  !****************************************************************************
  subroutine EventInfo_HiLep_Init(NumEnsembles,MassNum)
    implicit none
    integer,intent(in) :: NumEnsembles,MassNum
    logical,save :: initFlag=.true.

    if (initFlag) then
       allocate(ARR_HiLep(1:NumEnsembles,-1:MassNum))
       initFlag=.false.
    end if
    ARR_HiLep%flagOK=.False.
    ARR_HiLep%Weight=0.
    ARR_HiLep%nu=0.
    ARR_HiLep%Q2=0.
    ARR_HiLep%eps=0.
    ARR_HiLep%EventType=0
    ARR_HiLep%phi_Lepton=0.
    ARR_HiLep%W=0.
    ARR_HiLep%pos(1)=0.
    ARR_HiLep%pos(2)=0.
    ARR_HiLep%pos(3)=0.

    beamEnergy = 0.

  end subroutine EventInfo_HiLep_Init


  !****************************************************************************
  !****s* EventInfo_HiLep/EventInfo_HiLep_Store
  ! NAME
  ! subroutine EventInfo_HiLep_Store(i,j,Weight,nu,Q2,eps,EventType,Ebeam,phi_Lepton,W,pos)
  !
  ! PURPOSE
  ! Store the event info connected with ensemble "i" and nucleon "j":
  !
  !
  ! INPUTS
  ! * integer :: i -- actual Ensemble
  ! * integer :: j -- actual Nucleon
  ! * real    :: Weight -- weigt of event
  ! * real    :: nu,Q2,eps -- photonic variables
  ! * integer :: EventType -- cf. module Electron_origin
  ! * real    :: Ebeam        [OPTIONAL] -- energy of lepton beam
  ! * real    :: phi_Lepton   [OPTIONAL] -- phi-angle of lepton
  ! * real    :: W            [OPTIONAL] -- the value including fermi
  !   mom and potentials
  ! * real,dimension(3):: pos [OPTIONAL] -- position of the interaction
  ! OUTPUT
  ! ---
  !****************************************************************************
  subroutine EventInfo_HiLep_Store(i,j,Weight,nu,Q2,eps,EventType,Ebeam,phi_Lepton,W,pos)
    implicit none
    integer,intent(in)          :: i,j,EventType
    real,   intent(in)          :: nu,Q2,eps,Weight
    real,   intent(in),optional :: Ebeam,phi_Lepton,W
    real,dimension(3),intent(in),optional :: pos

    ARR_HiLep(i,j)%flagOK      =.TRUE.
    ARR_HiLep(i,j)%Weight      = Weight
    ARR_HiLep(i,j)%nu          = nu
    ARR_HiLep(i,j)%Q2          = Q2
    ARR_HiLep(i,j)%eps         = eps
    ARR_HiLep(i,j)%EventType   = EventType

!    write(*,'(A,2i4.0,1P,2e12.5)') 'EventInfo_HiLep_Store: ',i,j,nu,Q2

    if (present(Ebeam)) beamEnergy=Ebeam
    if (present(phi_Lepton))   ARR_HiLep(i,j)%phi_Lepton   = phi_Lepton
    if (present(W))            ARR_HiLep(i,j)%W            = W
    if (present(pos))          ARR_HiLep(i,j)%pos          = pos

  end subroutine EventInfo_HiLep_Store


  !****************************************************************************
  !****f* EventInfo_HiLep/EventInfo_HiLep_Get
  ! NAME
  ! logical function EventInfo_HiLep_Get(i,j,Weight,nu,Q2,eps,EventType,Ebeam,phi_Lepton,W,pos)
  !
  ! PURPOSE
  ! Get the event info stored connected with ensemble "i" and nucleon "j".
  !
  ! if j>1000, the ensemble index is set by "i = j/1000"
  !
  ! INPUTS
  ! * integer :: i -- actual Ensemble
  ! * integer :: j -- actual Nucleon
  !
  ! OUTPUT
  ! * real    :: Weight -- weigt of event
  ! * real    :: nu,Q2,eps -- photonic variables
  ! * integer :: EventType -- cf. module Electron_origin
  ! * real    :: Ebeam        [OPTIONAL] -- energy of lepton beam
  ! * real    :: phi_Lepton   [OPTIONAL] -- phi-angle of lepton
  ! * real    :: W            [OPTIONAL] -- the value including fermi
  !   mom and potentials
  ! * real,dimension(3):: pos [OPTIONAL] -- position of the interaction
  !
  !****************************************************************************
  logical function EventInfo_HiLep_Get(i,j,Weight,nu,Q2,eps,EventType,Ebeam,phi_Lepton,W,pos)
    implicit none
    integer,intent(in)           :: i,j
    integer,intent(out)          :: EventType
    real,   intent(out)          :: nu,Q2,eps,Weight
    real,   intent(out),optional :: Ebeam,phi_Lepton,W
    real,dimension(3),intent(out),optional :: pos

    integer :: ii,jj


    EventInfo_HiLep_Get = .FALSE.

    if (.not.allocated(ARR_HiLep)) return

    if (j<1000) then
       ii = i
       jj = j
    else
       ii = j/1000       ! == number Ensemble
       jj = j-ii*1000    ! == number realParticle
    end if

    Weight      = ARR_HiLep(ii,jj)%Weight
    nu          = ARR_HiLep(ii,jj)%nu
    Q2          = ARR_HiLep(ii,jj)%Q2
    eps         = ARR_HiLep(ii,jj)%eps
    EventType   = ARR_HiLep(ii,jj)%EventType

    if (present(Ebeam)) Ebeam=beamEnergy
    if (present(phi_Lepton))   phi_Lepton   = ARR_HiLep(ii,jj)%phi_Lepton
    if (present(W))            W            = ARR_HiLep(ii,jj)%W
    if (present(pos))          pos          = ARR_HiLep(ii,jj)%pos

    EventInfo_HiLep_Get = ARR_HiLep(ii,jj)%flagOK

!    write(*,'(A,2i4.0,1P,2e12.5)') 'EventInfo_HiLep_Get:   ',ii,jj,nu,Q2


  end function EventInfo_HiLep_Get


  !****************************************************************************
  !****f* EventInfo_HiLep/EventInfo_HiLep_GetBeamEnergy
  ! NAME
  ! real function EventInfo_HiLep_GetBeamEnergy()
  !
  ! PURPOSE
  ! return the value of the variable "beamEnergy"
  !
  ! INPUTS
  ! ---
  !
  ! OUTPUT
  ! value of the variable "beamEnergy"
  !****************************************************************************
!   real function EventInfo_HiLep_GetBeamEnergy()
!     EventInfo_HiLep_GetBeamEnergy=beamEnergy
!     return
!   end function EventInfo_HiLep_GetBeamEnergy

  !****************************************************************************
  !****s* EventInfo_HiLep/EventInfo_HiLep_Dump
  ! NAME
  ! subroutine EventInfo_HiLep_Dump
  !
  ! PURPOSE
  ! Dump the information to a file
  !****************************************************************************
  subroutine EventInfo_HiLep_Dump
    integer :: i,j

    open(103,file='EventInfo_HiLep.dat',position='Append')

    do i=1,size(ARR_HiLep,dim=1)
       do j=-1,ubound(ARR_HiLep,dim=2)
          if (.not.ARR_HiLep(i,j)%flagOK) cycle

          write(103,'(2i7,1P,4e13.5,0P,i5,1P,5e13.5)') i,j, &
               & ARR_HiLep(i,j)%Weight, &
               & ARR_HiLep(i,j)%nu, ARR_HiLep(i,j)%Q2, &
               & ARR_HiLep(i,j)%eps, ARR_HiLep(i,j)%EventType, &
               & ARR_HiLep(i,j)%phi_Lepton, ARR_HiLep(i,j)%W, &
               & ARR_HiLep(i,j)%pos

       end do
    end do

    close(103)

  end subroutine EventInfo_HiLep_Dump

end module EventInfo_HiLep
