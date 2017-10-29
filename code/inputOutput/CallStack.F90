!******************************************************************************
!****m* /CallStack
! NAME
! module CallStack
!
! PURPOSE
! This module provides a wrapper around compiler-specific routines.
!******************************************************************************
module CallStack

  implicit none
  private

  public :: traceback, system_command

contains

  !****************************************************************************
  !****s* CallStack/TRACEBACK
  ! NAME
  ! subroutine TRACEBACK(string,user_exit_code)
  !
  ! PURPOSE
  ! Write out the call stack of the program, depending on the compiler.
  ! One currently gets a backtrace with the following compilers:
  ! * ifort (via the routine TRACEBACKQQ)
  ! * gfortran 4.7 or higher (via ABORT)
  ! See also the documentation of these routines in the respective compiler
  ! manual.
  !
  ! INPUTS
  ! * character*(*), optional :: string -- header line to write
  ! * integer, optional :: user_exit_code -- code whether return or stop
  !   program
  ! * By specifying a user exit code of -1, control returns to the calling
  !   program. Specifying a user exit code with a positive value requests that
  !   specified value be returned to the operating system. The default value
  !   is 0, which causes the application to abort execution.
  !****************************************************************************
  subroutine TRACEBACK(string,user_exit_code)
#ifdef ifort
    use IFCORE
#endif

    character*(*), intent(in), optional:: string
    integer, intent(in), optional :: user_exit_code

#ifdef ifort
!#warning "compiling with ifort"

    if (present(string)) then
       if (present(user_exit_code)) then
          call TRACEBACKQQ(string=string,user_exit_code=user_exit_code)
       else
          call TRACEBACKQQ(string=string)
       end if
    else
       if (present(user_exit_code)) then
          call TRACEBACKQQ(user_exit_code=user_exit_code)
       else
          call TRACEBACKQQ()
       end if
    end if

#else
!#warning "not compiling with ifort"
    if (present(string)) then
       write(*,'(A)') string
       flush(6)
    end if
# ifdef __GFORTRAN__
    if (present(user_exit_code)) then
       if (user_exit_code.eq.-1) then
          write(*,*) '--- no call stack trace possible ---'
          return
       end if
    end if
    ! ABORT is a GNU extension and gives a backtrace with gfortran 4.7 and above
    call abort()
# else

    write(*,*) '--- no call stack trace possible ---'

    if (present(user_exit_code)) then
       if (user_exit_code.eq.-1) return
       stop 123
    end if
    stop
# endif
#endif

  end subroutine TRACEBACK


  !****************************************************************************
  !****s* CallStack/SYSTEM_COMMAND
  ! NAME
  ! subroutine SYSTEM_COMMAND(string)
  !
  ! PURPOSE
  ! Runs a shell command, as given by the argument 'string'.
  !
  ! INPUTS
  ! * character*(*) :: string -- shell command to execute
  !
  ! NOTES
  ! This routine is compiler-dependent. If compiled with ifort, the
  ! Intel-specific function SYSTEMQQ is called. Otherwise the F08 intrinsic
  ! procedure EXECUTE_COMMAND_LINE is used.
  !****************************************************************************
  subroutine SYSTEM_COMMAND(string)
#ifdef ifort
    use IFPORT
#endif

    character*(*), intent(in):: string

#ifdef ifort
    logical :: res
    res = SYSTEMQQ(string)
#elif __GFORTRAN__
    call EXECUTE_COMMAND_LINE(string)
#else
    call SYSTEM(string)
#endif
  end subroutine SYSTEM_COMMAND

end module CallStack
