!***************************************************************************
!****m* /TmunuDefinition
! NAME
! module TmunuDefinition
! PURPOSE
! Here type(tTmunuNmu) is defined.
! This module includes also functions for this type.
!***************************************************************************
module TmunuDefinition
  implicit none
  private

  !***************************************************************************
  !****t* BoxAnalysis/tTmunuNmu
  ! NAME
  ! type tTmunuNmu
  ! PURPOSE
  ! This stores Tmunu, Nmu, Jmu:
  ! * Tmunu = 1/(nEns*V) sum_i (pmu pnu)/p0
  ! * Nmu = 1/(nEns*V) sum_i pmu/p0
  ! * Jmu = 1/(nEns*V) sum_i pmu/p0 * q_i
  ! The normalization has to be fixed manually.
  !
  ! It is assumed, that the baryon and the strangeness current are
  ! proportional to Nmu. Thus only the zeorth component is stored
  !
  ! SOURCE
  !
  type tTmunuNmu
     real, dimension(0:9) :: Tmunu = 0. ! the Tmunu tensor
     real, dimension(0:3) :: Nmu   = 0. ! the particle flow
     real, dimension(0:3) :: Jmu   = 0. ! the electrical charge flow
     real :: B  = 0. ! the baryon number
     real :: S  = 0. ! the strangeness
  end type tTmunuNmu
  !***************************************************************************

  ! string constants may be broken over multiple continuation lines:
  character(len=*), parameter :: headTmunu = "# 1: timestep &
         &2: T00 3: T11 4: T22 5: T33 &
         &6: T01 7: T02 8: T03 &
         &9: T21 10: T31 11: T32 &
         &12: N0 13: N1 14: N2 15: N3 &
         &16: J0 17: J1 18: J2 19: J3 &
         &20: B0 21: S0"

  PUBLIC :: tTmunuNmu
  PUBLIC :: headTmunu
  PUBLIC :: fillTmunu


contains

  !***************************************************************************
  !****s* TmunuDefinition/fillTmunu
  ! NAME
  ! subroutine fillTmunu(TmunuNmu, part, weight)
  ! PURPOSE
  ! add the momentum of the given particle to the entries of Tmunu, Nmu, Jmu
  !***************************************************************************
  subroutine fillTmunu(TmunuNmu, part, weight)
    use particleDefinition
    use particleProperties
    use IdTable
!    use Hagedorn, only: Hagedorn_IDtoBSI

    type(tTmunuNmu), intent(inOut) :: TmunuNmu
    type(particle), intent(in) :: part
    real, intent(in), optional :: weight

    integer :: qB,qS,qC,qIx2

    real :: w, oneE

    w = 1.0
    if (present(weight)) w = weight

    oneE = w/part%momentum(0)

    TmunuNmu%Tmunu(0:3) = TmunuNmu%Tmunu(0:3) + part%momentum(0:3)**2*oneE ! T00,T11,T22,T33
    TmunuNmu%Tmunu(4:6) = TmunuNmu%Tmunu(4:6) + part%momentum(1:3)*w ! T01,T02,T03
    TmunuNmu%Tmunu(7)   = TmunuNmu%Tmunu(7)   + part%momentum(2)*part%momentum(1)*oneE ! T21
    TmunuNmu%Tmunu(8)   = TmunuNmu%Tmunu(8)   + part%momentum(3)*part%momentum(1)*oneE ! T31
    TmunuNmu%Tmunu(9)   = TmunuNmu%Tmunu(9)   + part%momentum(3)*part%momentum(2)*oneE ! T32

    TmunuNmu%Nmu(0:3) = TmunuNmu%Nmu(0:3) + part%momentum(0:3) * oneE

    if (part%charge .ne. 0) then
       TmunuNmu%Jmu(0:3) = TmunuNmu%Jmu(0:3) + part%momentum(0:3) * part%charge * oneE
    endif

    !    call Hagedorn_IDtoBSI(part, qB,qS,qC,qIx2)
    qB = 0
    if (isBaryon(part%ID)) qB = 1
    qS = hadron(part%ID)%strangeness
    if (part%antiparticle) then
       qB = -qB
       qS = -qS
    end if

    TmunuNmu%B = TmunuNmu%B + qB*w
    TmunuNmu%S = TmunuNmu%S + qS*w

  end subroutine fillTmunu


end module TmunuDefinition
