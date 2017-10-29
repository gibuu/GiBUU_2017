!******************************************************************************
!****m* /Electron_origin
! NAME
! module Electron_origin
!
! PURPOSE
! Stores & returns the type of event for electron nucleon scattering
! according to %firstEvent
!
! The numbering scheme is:
! * -3000 : BaB
! * -2000 : Manni
! * -1000 : Elastic
! *  0xxx : Low Energy
! *  1000 : Fritiof
! *  2xxx : Pythia
!
! The "Low Energy" Subtopics (0xxx) are:
! * 00xx : gamma N -> xx  with xx = resonance ID
! * 010V : gamma N -> V N        ! == Nucleon
! * 011V : gamma N -> V N pi
! * 020V : gamma N -> V Delta    ! == Delta
! * 021V : gamma N -> V Delta pi
! * 0310 : gamma N -> Lambda K   ! == Strangeness
! * 0320 : gamma N -> Sigma K
! * 0330 : gamma N -> K+ K- N
! * 0410 : gamma N -> pi0 eta ...! == explicit channels
! * ... (free entries)
! * 0510 : gamma N N -> N N      ! == 2p2h
! * 0520 : gamma N N -> N Delta
! * ... (free entries)
! * 0810 : gamma N -> pi N       ! == Background
! * 0820 : gamma N -> pi pi N
! * 082V : gamma N -> pi pi N   via  V N
! * 0829 : gamma N -> pi pi N   via  pi Delta
!
! The vector meson "V" is given by:
! * 1 : rho
! * 2 : omega
! * 3 : phi
! * 4 : J/psi
!
! The "Pythia" Subtopics (2xxx) are (cf. MSTI(9) ):
! * 2001 : VMD
! * 2002 : direct
! * 2003 : anomalous
! * 2004 : DIS
!******************************************************************************

module Electron_origin

  implicit none

  private
  !****************************************************************************
  !****g* Electron_origin/origin_XXX
  ! PURPOSE
  ! Constants defined for ease of coding
  !
  ! SOURCE
  !
  integer, parameter :: origin_singlePi    =  810 ! gamma N -> pi N
  integer, parameter :: origin_doublePi    =  820 ! gamma N -> pi pi N
  integer, parameter :: origin_DIS         = 2004 ! DIS
  integer, parameter :: origin_vecmes      =  100 ! gamma N -> V N
  integer, parameter :: origin_vecmes_rho  =  101 !
  integer, parameter :: origin_vecmes_omega=  102 !
  integer, parameter :: origin_vecmes_phi  =  103 !
  integer, parameter :: origin_vecmes_Delta=  200 ! gamma N -> V Delta
  integer, parameter :: origin_pi0eta      =  410 ! gamma N -> pi0 eta ...
  integer, parameter :: origin_2p2hQE      =  510 ! gamma N N -> N N
  integer, parameter :: origin_2p2hDelta   =  520 ! gamma N N -> N Delta
  integer, parameter :: origin_fritiof     = 1000 ! gamma N -> X (via FRITIOF)
  !****************************************************************************

  integer , dimension (:), allocatable, save :: table

  public :: le_whichOrigin
  public :: le_isOmegaEvent
  public :: le_isRhoEvent
  public :: le_isPhiEvent
  public :: le_isDisEvent
  public :: le_isResEvent
  public :: le_isPiPiEvent
  public :: le_isPi0EtaEvent
  public :: cleanup

  public :: origin_singlePi,origin_doublePi,origin_DIS
  public :: origin_vecmes
  public :: origin_vecmes_rho,origin_vecmes_omega,origin_vecmes_phi,origin_vecmes_Delta
  public :: origin_pi0eta, origin_fritiof
  public :: origin_2p2hQE, origin_2p2hDelta


contains

  !****************************************************************************
  !****s* Electron_origin/le_whichOrigin
  ! NAME
  ! subroutine le_whichOrigin(Switch,index,info)
  !
  ! PURPOSE
  ! Stores & returns the type of event according to %firstEvent
  !
  ! INPUTS
  ! * integer :: index
  ! * integer,optional :: info
  ! * integer :: switch
  !
  ! USAGE
  ! If switch= ... :
  ! * 0 -- initialize the subroutine
  ! * 1 -- store the information ("info") of the type of event for  firstEvent=index
  ! * 2 -- returne the information ("info") of the type of event for  firstEvent=index
  !****************************************************************************
  subroutine le_whichOrigin(Switch,index,info)

    integer , dimension (:), allocatable :: dummy
    integer :: index
    integer,optional :: info
    integer :: switch

    select case (switch)
    case (0)
       if (allocated(table)) deallocate(table)
       allocate(table(0:index))
       table=0
    case (1)
       if (present(info)) then
          allo: do
             if (index.gt.ubound(table,dim=1)) then
                allocate(dummy(lbound(table,dim=1):uBound(table,dim=1)))
                dummy=table
                deallocate(table)
                allocate(table(lbound(dummy,dim=1):2*uBound(dummy,dim=1)))
                table(lbound(dummy,dim=1):uBound(dummy,dim=1))=dummy
                deallocate(dummy)
             else
                exit allo
             end if
          end do allo
          table(index)=info
       else
          write(*,*) 'error in le_whichOrigin',index,present(info),switch
       end if
    case (2)
       if (present(info)) then
          if (allocated(table)) then
             info=table(index)
          else
             info=999999999 ! no info available
          end if
       else
          write(*,*) 'error in le_whichOrigin',index,present(info),switch
       end if

    case default
       write(*,*) 'error in le_whichOrigin',switch
    end select


  end subroutine le_whichOrigin


  subroutine cleanUp
    if (allocated(table)) deallocate(table)
  end subroutine


  !****************************************************************************
  !****s* Electron_origin/le_isResEvent
  ! NAME
  ! logical function le_isResEvent(firstEvent)
  !
  ! PURPOSE
  ! Returns .true. if the event with %firstEvent=firstEvent was a resonance event, else .false. .
  !
  ! INPUTS
  ! integer, intent(in) :: firstEvent
  !****************************************************************************
  logical function le_isResEvent(firstEvent)

    use IdTable, only: isBaryonResonance

    integer,intent(in) :: firstEvent
    integer :: info

    call le_whichOrigin(2,firstEvent,info)
    le_isResEvent = isBaryonResonance(info)

  end function le_isResEvent


  !****************************************************************************
  !****s* Electron_origin/le_isPi0EtaEvent
  ! NAME
  ! logical function le_isPi0EtaEvent(firstEvent)
  !
  ! PURPOSE
  ! Returns .true. if the event with %firstEvent=firstEvent was a pi^0 eta production event, else .false. .
  !
  ! INPUTS
  ! integer, intent(in) :: firstEvent
  !****************************************************************************
  logical function le_isPi0EtaEvent(firstEvent)

    integer,intent(in) :: firstEvent
    integer :: info

    call le_whichOrigin(2,firstEvent,info)
    le_isPi0EtaEvent=(info.eq.origin_pi0eta)

  end function le_isPi0EtaEvent


  !****************************************************************************
  !****s* Electron_origin/le_isPiPiEvent
  ! NAME
  ! logical function le_isPiPiEvent(firstEvent)
  !
  ! PURPOSE
  ! Returns .true. if the event with %firstEvent=firstEvent was a 2pi production event, else .false. .
  !
  ! INPUTS
  ! integer, intent(in) :: firstEvent
  !****************************************************************************
  logical function le_isPiPiEvent(firstEvent)

    integer,intent(in) :: firstEvent
    integer :: info

    call le_whichOrigin(2,firstEvent,info)
    le_isPiPiEvent=(info.eq.origin_doublePi)

  end function le_isPiPiEvent


  !****************************************************************************
  !****s* Electron_origin/le_isOmegaEvent
  ! NAME
  ! logical function le_isOmegaEvent(firstEvent)
  !
  ! PURPOSE
  ! Returns .true. if the event with %firstEvent=firstEvent was a omega production event, else .false. .
  !
  ! INPUTS
  ! integer, intent(in) :: firstEvent
  !****************************************************************************
  logical function le_isOmegaEvent(firstEvent)

    integer,intent(in) :: firstEvent
    integer :: info

    call le_whichOrigin(2,firstEvent,info)
    le_isOmegaEvent=(info.eq.origin_vecmes_omega)

  end function le_isOmegaEvent
  !****************************************************************************
  !****s* Electron_origin/le_isphiEvent
  ! NAME
  ! logical function le_isPhiEvent(firstEvent)
  !
  ! PURPOSE
  ! Returns .true. if the event with %firstEvent=firstEvent was a phi production event, else .false. .
  !
  ! INPUTS
  ! integer, intent(in) :: firstEvent
  !****************************************************************************
  logical function le_isPhiEvent(firstEvent)

    integer,intent(in) :: firstEvent
    integer :: info

    call le_whichOrigin(2,firstEvent,info)
    le_isPhiEvent=(info.eq.origin_vecmes_phi)

  end function le_isPhiEvent


  !****************************************************************************
  !****s* Electron_origin/le_isDISEvent
  ! NAME
  ! logical function le_isDISEvent(firstEvent)
  !
  ! PURPOSE
  ! Returns .true. if the event with %firstEvent=firstEvent was a DIS event, else .false. .
  !
  ! INPUTS
  ! integer, intent(in) :: firstEvent
  !****************************************************************************
  logical function le_isDISEvent(firstEvent)

    integer,intent(in) :: firstEvent
    integer :: info

    call le_whichOrigin(2,firstEvent,info)
    le_isDISEvent=(info.eq.origin_DIS)

  end function le_isDISEvent

  !****************************************************************************
  !****s* Electron_origin/le_isRhoEvent
  ! NAME
  ! logical function le_isRhoEvent(firstEvent)
  !
  ! PURPOSE
  ! Returns .true. if the event with %firstEvent=firstEvent was a rho production event, else .false. .
  !
  ! INPUTS
  ! integer, intent(in) :: firstEvent
  !****************************************************************************
  logical function le_isRhoEvent(firstEvent)

    integer,intent(in) :: firstEvent
    integer :: info

    call le_whichOrigin(2,firstEvent,info)
    le_isRhoEvent=(info.eq.origin_vecmes_rho)

  end function le_isRhoEvent


end module Electron_origin
