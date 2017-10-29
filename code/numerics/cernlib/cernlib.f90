!******************************************************************************
!****m* /cern_lib
! NAME
! module CERN_Lib
! PURPOSE
! This module is a wrapper for various CERNLIB routines:
! * 'DZEROX'
! * 'DCAUCH'
!******************************************************************************
module cern_lib

  Interface

     !*************************************************************************
     !****f* cern_lib/DZEROX
     ! NAME
     ! double precision function DZEROX(A0,B0,EPS,MAXF,F,MODE)
     ! PURPOSE
     ! Finding a Zero of a Function
     ! INPUTS
     ! * double precision :: A0, B0 -- specify the end points of the search interval
     ! * double precision :: EPS    -- accuracy parameter
     ! * integer          :: MAXF   -- maximum permitted number of references to the function F within the iteration loop
     ! * double precision, external :: F -- user-supplied FUNCTION
     ! * integer          :: MODE   -- Mode=1 or Mode=2 defines the algorithm for finding x_0
     ! OUTPUT
     ! * function value represents x_0, the found zero
     ! NOTES
     ! cf. original CERNLIB documentation
     !*************************************************************************
     double precision function DZEROX(A0,B0,EPS,MAXF,F,MODE)
       real(8), intent(in) :: A0, B0
       real(8), intent(in) :: EPS
       integer, intent(in) :: MAXF
       double precision, external :: F
       integer, intent(in) :: mode
     end function DZEROX

  end Interface


  Interface
     FUNCTION DCAUCH(F,A,B,S,EPS)
       IMPLICIT REAL(8) (A-H,O-Z)
       REAL(8),EXTERNAL :: F
     end FUNCTION DCAUCH
  end interface

end module cern_lib
