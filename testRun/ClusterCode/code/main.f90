!***************************************************************************
!****m* /Main
! NAME
! program Main
!
! PURPOSE
! This is the main program for fragment formation.
!***************************************************************************

program Main

  use InputGeneral, only: Get_inputGeneral, TheModel
  use ModelDefinitions, only: Coalescence,SMM
  use coalescenceModule, only: main_Coalescence
  use smmModule, only: main_SMM
  use version, only: printVersion

  implicit none

  !------------------------------------------------------------------------*
  write(*,"(/,78('='))")
  call PrintVersion
  write(*,"(78('='),/)")
  !-------------------------------------------------------------------------*
  call Get_inputGeneral
  !-------------------------------------------------------------------------*
  Select Case(TheModel)

     Case(Coalescence)
        call main_Coalescence

     Case(SMM)
        call main_SMM

     Case default
        write(*,*) 'main.f90: no valid fragmentation model',TheModel
        write(*,*) 'STOP'
        STOP
        
     end Select

end program Main
