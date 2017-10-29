!***************************************************************************
!****m* /templateModule
! NAME
! module templateModule
!
! PURPOSE
! This module defines ...
!
! INPUTS
! The Namelist "templateModule_nl" in the Jobcard.
!
! NOTES
! In order to get best results from RoboDoc, you should respect the spelling
! of modules, subroutines, functions etc troughout all doku lines: 
! usage of capital or small letters should be consistent all over the files! 
!***************************************************************************

!***************************************************************************
!****n* templateModule/templateModule_nl
! NAME
! NAMELIST /templateModule_nl/
! PURPOSE
! This Namelist for module "templateModule" includes:
! * var1 -- The variable Nr1
! * var2 -- The variable Nr2
!***************************************************************************
module templateModule

  PRIVATE

  !*************************************************************************
  !****g* templateModule/var1
  ! SOURCE
  !
  integer, save :: var1
  ! PURPOSE
  ! The variable Nr1
  !
  ! This variable has the function blablabla
  !*************************************************************************

  !*************************************************************************
  !****g* templateModule/var2
  ! SOURCE
  !
  integer, save :: var2
  ! PURPOSE
  ! The variable Nr2
  !
  ! This variable has the function blablabla
  !*************************************************************************

  PUBLIC :: Sub1, Fun1

contains
  
  !*************************************************************************
  !****s* templateModule/Sub1
  ! NAME
  ! subroutine Sub1(par1,par2,par3)
  !
  ! PURPOSE
  ! The first Subroutine.
  !
  ! This Subroutine calculates ...
  ! ...
  !
  ! INPUTS
  ! * integer :: par1 -- The parameter Nr1
  ! * real    :: par2 -- The parameter Nr2
  !
  ! OUTPUT
  ! * real    :: par2 -- The parameter Nr2
  !                      is modified according ...
  ! * integer :: par3 -- The parameter Nr3
  !
  ! * The global parameter Xwidget is set to .TRUE.
  !
  !*************************************************************************
  subroutine Sub1(par1,par2,par3)

    integer,     intent(in)    :: par1 ! please do not place any explanation
    real,        intent(inout) :: par2 ! here. do it in the header
    integer,     intent(out)   :: par3

    ! ...

  end subroutine Sub1

  !*************************************************************************
  !****is* templateModule/Sub2
  ! NAME
  ! subroutine Sub2
  !
  ! PURPOSE
  ! The second Subroutine.
  !
  ! This Subroutine calculates ...
  ! ...
  ! INPUTS
  ! (-- no input --)
  !
  ! OUTPUT
  ! (-- no output --)
  !
  ! SIDE EFFECTS
  ! The global variable xyz is changed to ...
  !
  !*************************************************************************
  subroutine Sub2
  end subroutine Sub2

  !*************************************************************************
  !****f* templateModule/Fun1
  ! NAME
  ! real function Fun1()
  !
  ! PURPOSE
  ! The first Function
  !
  ! This Function calculates ...
  ! ...
  ! INPUTS
  ! (-- no input --)
  !
  ! OUTPUT
  ! returns the value of sin(pi)
  !
  ! SIDE EFFECTS
  ! (-- no side effects --)
  !
  !*************************************************************************
  real function Fun1()
  end function Fun1




end module templateModule
