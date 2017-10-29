!******************************************************************************
!****m* /formfactors_A_input
! NAME
! module formfactors_A_input
!
! PURPOSE
! This module includes input to the formfactors.
! In addition it collects routines common to all formfactors_A_* modules
!
! INPUTS
! NAMELIST 'formfactors_pion'
!
! NOTES
! here is also the right place to define the routines 'get_i' and
! 'get_neighbors', which are multiple defined.
!******************************************************************************
module formfactors_A_input

  implicit none
  private

  !****************************************************************************
  !****g* formfactors_A_input/which_MaidVersion
  ! SOURCE
  !
  integer, save :: which_MaidVersion=2
  !
  ! PURPOSE
  ! choice of MAID version:  1=2003, 2=2007
  !****************************************************************************

  integer, parameter :: MAID2003=1
  integer, parameter :: MAID2007=2


  public :: get_Filename, ValidateVars

  logical, save :: firstTime = .true.

contains

  subroutine readInput
    use output, only: Write_ReadingInput

    !**************************************************************************
    !****n* formfactors_A_input/formfactors_pion
    ! NAME
    ! NAMELIST formfactors_pion
    ! PURPOSE
    ! This namelist includes the following switches:
    ! * which_MaidVersion
    !**************************************************************************
    NAMELIST /formfactors_pion/ which_MaidVersion

    integer :: ios

    call Write_ReadingInput("formfactors_pion",0)
    rewind(5)
    read(5,nml=formfactors_pion,IOSTAT=ios)
    call Write_ReadingInput("formfactors_pion",0,ios)
    select case (which_MaidVersion)
    case (MAID2003)
       write(*,*) 'MAID 2003 is being used for gamma* N -> N pi'
    case (MAID2007)
       write(*,*) 'MAID 2007 is being used for gamma* N -> N pi'
    case default
       write(*,*) 'Wrong MAID switch!! which_MAIDVersion', which_MaidVersion
       write(*,*) 'Must be :'
       write(*,*) MAID2003, 'for MAID 2003 or ', MAID2007,' for MAID2007'
       stop
    end select
    call Write_ReadingInput("formfactors_pion",1)

    firstTime = .false.

  end subroutine readInput


  function get_Filename (str)
    use inputGeneral, only: path_to_input
    character(len=4), intent(in) :: str
    character(len=200) :: get_Filename

    if (firstTime) call readinput

    select case (which_MAIDVersion)
    case (MAID2003)
      get_Filename = trim(path_to_Input) // '/electronNucleon/' // trim(str) // '_03_21_91_19.out.bz2'
    case (MAID2007)
      get_Filename = trim(path_to_Input) // '/electronNucleon/' // trim(str) // '_07_51_91_19.out.bz2'
    case default
      stop 'formfactors_A_input/get_Filename'
    end select

  end function get_Filename


  !****************************************************************************
  !****s* formfactors_A_input/ValidateVars
  ! NAME
  ! subroutine ValidateVars
  ! PURPOSE
  ! Check whether the the given parameters are in the valid range.
  ! OUTPUT
  ! The routine returns .false. if W is above the supported range (else .true.).
  ! It stops with an error if a violation of the theta or Q**2 bounds occurs.
  !****************************************************************************
  logical function ValidateVars (W, theta, Q2, Wmin, Wmax, thetaMin, thetaMax, Q2min, Q2max, ModName)
    use callstack, only: TraceBack
    use output, only: DoPR

    real, intent(inOut) :: W, Q2
    real, intent(in) :: theta, Wmin, Wmax, thetaMin, thetaMax, Q2min, Q2max
    character(*), intent(in) :: ModName

!     if ((W>Wmax) .and. (W<Wmax*1.1)) then
!        if (DoPR(2)) write(*,'(3A)') 'WARNING: ',trim(ModName),', W>W_max, but less than 10%.'
!        W = Wmax
!     end if

    if ((Q2>Q2max) .and. (Q2<Q2max*1.1)) then
       if (DoPR(2)) write(*,'(3A)') 'WARNING: ',trim(ModName),', Q2>Q2max, but less than 10%.'
       Q2 = Q2max
    end if

    if (W>Wmax) then
      ValidateVars = .false.
      return
    end if

    if ((theta>thetaMax) .or. (Q2>Q2max)) then
       write(*,'(/,3A)') 'ERROR: ',trim(ModName),', out of upper bounds! (in vs. bounds)'
       write(*,*) '...theta:', theta,thetaMax,(theta>thetaMax)
       write(*,*) '...W:    ', W, Wmax,(W>Wmax)
       write(*,*) '...Q2:   ', Q2,Q2max,(Q2>Q2max)
       write(*,*) 'Stop!'
       call TraceBack()
       stop
    end if

    if ((W<Wmin) .or. (theta<thetaMin) .or. (Q2<Q2min)) then
       write(*,'(/,3A)') 'ERROR: ',trim(ModName),', out of lower bounds! (in vs. bounds)'
       write(*,*) '...theta:', theta,thetaMin,(theta.lt.thetaMin)
       write(*,*) '...W:    ', W, Wmin,(W.lt.Wmin)
       write(*,*) '...Q2:   ', Q2,Q2min,(Q2.lt.Q2min)
       write(*,*) 'Stop!'
       call TraceBack()
       stop
    end if

    ValidateVars = .true.

  end function ValidateVars

end module formfactors_A_input
