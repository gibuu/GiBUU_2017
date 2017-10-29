!******************************************************************************
!****m* /finalStateModule
! NAME
! module finalStateModule
!
! PURPOSE
! This module implements the routines "massass" and "assmass" as switches
! to the full (old,slow) routines for momentum dependent spectral functions
! and to improved, fast routines, which are only valid for momentum independent
! spectral functions.
!******************************************************************************
module finalStateModule
  implicit none
  private

  public :: massAss
  public :: assMass

contains

  !****************************************************************************
  !****s* finalStateModule/massAss
  ! NAME
  ! subroutine massAss(srts,mediumAtColl,PartIn,PartOut,betaToLRF,betaToCM,L,flagOK)
  !
  ! PURPOSE
  ! This routine selects the masses of a 2-particle final state
  ! according to phase space x spectral functions.
  !
  ! We regard the process a+b->c+d.
  !
  ! INPUTS
  ! * real                           :: srts         -- sqrt(s)
  ! * type(medium)                   :: mediumAtColl -- medium information : density, temperature,...
  ! * type(particle), dimension(1:2) :: PartIn       -- incoming particles a and b
  ! * type(particle), dimension(1:2) :: PartOut      -- outgoing particles c and d
  !                                                  (only id's, charges and antiflags are used at input)
  ! * real, dimension(1:3)           :: betaToLRF    -- beta for boost to LRF
  ! * real, dimension(1:3)           :: betaToCM     -- beta for boost to CM-Frame
  ! * integer                        :: L            -- angular momentum in final state
  !
  ! OUTPUT
  ! * logical :: flagOK -- .true. if mass assignment was successful
  ! * type(particle), dimension(1:2) :: PartOut --
  !   final state particles with full kinematics in the CM frame.
  !   ID, momentum(1:3),charge and mass are now defined.
  !****************************************************************************
  subroutine massAss(srts,mediumAtColl,PartIn,PartOut,betaToLRF,betaToCM,L,flagOK)

    use IdTable, only: nucleon,delta
    use mediumDefinition
    use particleDefinition
    use finalState_Full, only: massAss_Full

    real,                           intent(in)    :: srts
    type(medium),                   intent(in)    :: mediumAtColl
    type(particle), dimension(1:2), intent(in)    :: PartIn
    real, dimension(1:3),           intent(in)    :: betaToLRF
    real, dimension(1:3),           intent(in)    :: betaToCM
    integer,                        intent(in)    :: L

    type(particle), dimension(1:2), intent(inout) :: PartOut
    logical,                        intent(out)   :: flagOK

    ! special treatment of N N -> N Delta because here model for matrix element
    ! is available:
    if (PartIn(1)%ID==nucleon .and. PartIn(2)%ID==nucleon .and.   &
        ((PartOut(1)%ID==nucleon .and. PartOut(2)%ID==Delta) .or. &
         (PartOut(2)%ID==nucleon .and. PartOut(1)%ID==Delta))) then

      call massass_NN_NDelta (srts, mediumAtColl, PartIn, PartOut, betaToCM, flagOK)
    else
      call massAss_Full(srts,mediumAtColl,PartIn,PartOut,betaToLRF,betaToCM,L,flagOK)
    end if

  end subroutine massAss


  !****************************************************************************
  !****s* finalStateModule/assMass
  ! NAME
  ! subroutine assMass(srts,mediumAtColl,PartIn,PartOut,spotOut,betaToLRF,betaToCM,flag)
  !
  ! This routine selects the masses of a 3 particle final state
  ! according to phase space x spectral functions.
  !
  ! We regard the process a+b->c+d+e.
  !
  ! INPUTS
  ! * real                :: srts       -- sqrt(s)
  ! * type(medium)        :: mediumAtColl -- medium information : density, temperature,...
  ! * type(particle), dimension(1:2) :: PartIn     -- incoming particles a and b
  ! * type(particle), dimension(1:3) :: PartOut  -- outgoing particles c,d and e
  !   (only id's, charges and antiflags are used at input)
  ! * real, dimension(1:3)           :: spotOUT    -- scalar potential of produced particles
  ! * real, dimension(1:3)           :: betaToLRF  -- beta for boost to LRF
  ! * real, dimension(1:3)           :: betaToCM   -- beta for boost to CM-Frame
  !
  ! RESULT
  ! * logical :: flagOK   -- set to .true. is mass assignment was successful
  ! * type(particle), dimension (1:3) ,intent(out) :: PartOut  --
  !   final state particles with almost full kinematics in CM frame.
  !   ID, momentum(1:3),charge and mass are defined.
  !****************************************************************************
  subroutine assMass(srts,mediumAtColl,PartIn,PartOut,spotOut,betaToLRF,betaToCM,flagOK)
    use mediumDefinition
    use particleDefinition
    use finalState_Full, only: assMass_Full

    real,                               intent(in)    :: srts
    type(medium),                       intent(in)    :: mediumAtColl
    real, dimension(1:3),               intent(in)    :: spotOUT
    type(particle), dimension (1:2),    intent(in)    :: PartIn
    real, dimension(1:3),               intent(in)    :: betaToLRF
    real, dimension(1:3),               intent(in)    :: betaToCM
    type(particle), dimension (1:3),    intent(inout) :: PartOut
    logical,                            intent(out)   :: flagOK

    call assMass_Full(srts,mediumAtColl,PartIn,PartOut,spotOut,betaToLRF,betaToCM,flagOK)

  end subroutine assMass


  !****************************************************************************
  !****************************************************************************
  subroutine massAss_NN_NDelta (srts, mediumAtColl, PartIn, PartOut, betaToCM, flagOK)

    use mediumDefinition
    use particleDefinition
    use IdTable, only: nucleon,delta
    use particleProperties, only: hadron
    use CALLSTACK, only: TRACEBACK
    use dimi, only: massNN_NDelta
    use NbarN_to_NbarDelta, only: massNbarN_NbarDelta
    use twoBodyTools, only: pCM
    use winkelVerteilung, only: winkel
    use constants, only: mN

    real,                           intent(in)    :: srts
    type(medium),                   intent(in)    :: mediumAtColl
    type(particle), dimension(1:2), intent(in)    :: PartIn
    type(particle), dimension(1:2), intent(inout) :: PartOut
    real, dimension(1:3),           intent(in)    :: betaToCM
    logical,                        intent(out)   :: flagOK

    real :: massDelta, p_cd
    real, dimension(1:3) :: pscatt

    flagOK = .false.

    if (srts.lt.(mN+hadron(delta)%minmass)) return

    if (.not.partIn(1)%antiParticle .and. .not.partIn(2)%antiParticle) then
      massDelta = massNN_NDelta(srts)
    else if (partIn(1)%antiParticle .neqv. partIn(2)%antiParticle) then
      massDelta = massNbarN_NbarDelta(srts)
    else
       call TRACEBACK('Anti-Anti not implemented')
    end if

    if (PartOut(1)%ID==nucleon) then
       partOut(1:2)%mass=(/mN, massDelta/)
    else
       partOut(1:2)%mass=(/massDelta, mN/)
    end if

    pscatt = winkel (partIn, partOut, srts, betaToCM, mediumAtColl)

    p_cd = pCM(srts,mN,massDelta)
    partOut(1)%momentum(1:3) =  pscatt*p_cd
    partOut(2)%momentum(1:3) = -pscatt*p_cd

    partOut(1:2)%momentum(0) = sqrt(partOut(1:2)%mass**2 + p_cd**2)

    flagOK = .true.

  end subroutine massAss_NN_NDelta

end module finalStateModule
