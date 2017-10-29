!***************************************************************************
!****m* /random
! NAME
! module random
! FUNCTION
! Contains all information and routines, which are necessary for the
! random number generation.
!
! NOTES
! Formerly we used the random number generator RAN3 from Numerical Recipes.
! We switched to the generator "PYR" provided by Pythia, which is a rewrite
! of the RANMAR routine from the Cernlib.
!
! Contrary to the old RAN3 behaviour, the new rn() excludes the endpoints 
! 0 and 1 !
!
! Copied from the GiBUU code.
!***************************************************************************
module random

  
  Private :: initRan

  logical,save :: resetRandom = .false. ! reread random generator, used by setRandom, useful for debugging

  !*************************************************************************
  !****f* random/rnExp
  ! NAME
  ! real function rnExp(A)
  ! real function rnExp(A,x1,x2)
  ! FUNCTION
  ! Evaluates a random number x according to exp(-A*x) distribution.
  ! If given, x is restricted to lie between x1 and x2
  ! INPUTS
  ! * real :: A -- slope parameter
  ! * real :: x1 -- minimal x-value
  ! * real :: x2 -- maximal x-value
  ! RESULT
  ! random number
  !*************************************************************************
  interface rnExp
     module procedure rnExp1, rnExp2
  end interface

  logical,save,private :: first=.true.

contains


  !*************************************************************************
  !****f* random/InitRan
  ! NAME
  ! subroutine InitRan
  ! FUNCTION
  ! Reads random seed out of namelist 'initRandom' in jobcard and initializes
  ! the random number generator
  !*************************************************************************
  subroutine InitRan

    IMPLICIT NONE

    !***********************************************************************
    !****g* InitRan/SEED
    ! PURPOSE
    ! Random seed
    ! SOURCE
    !
    integer,save ::  SEED=970525
    !***********************************************************************

    integer :: ios ! Flag to check that namelist is available

    !***********************************************************************
    !****n* InitRan/initRandom
    ! NAME
    ! NAMELIST /initRandom/ 
    ! PURPOSE
    ! Includes the input variable:
    ! * SEED
    !***********************************************************************
    NAMELIST /initRandom/ Seed,resetRandom

    rewind(5)
    read(5,nml=initRandom,IOSTAT=IOS)

    write(*,*) 'Seed: ',Seed
    write(*,*) 'resetRandom: ',resetRandom

    call InitPYR(Seed)

  end subroutine InitRan


  !*************************************************************************
  !****s* random/SetRandom
  ! NAME
  ! subroutine SetRandom
  ! PURPOSE
  ! write out/read random number generators
  ! NOTES
  ! This routine provides a shortcut for MC debugging:
  ! if an error occurs in run 12345 of 99999 runs/energy after 3 days
  ! of CPU time, the idea is to reset the random generators by reading
  ! in some files, that already run 1 reproduces this error.
  !
  ! In order to avoid some re-ordering of the random number lists, only
  ! in run number 3 the previously written files are read.
  ! The code stops in run 4.
  !*************************************************************************
  subroutine setRandom
    implicit none

    logical ReadPYR
    integer,save :: nRead = 0

    common/DebugWRT/ DoWRT
    logical DoWRT
    save/DebugWRT/

    if (first) then
       ! Initializes the random number generator 
       ! when it's called for the first time
       call InitRan
       first=.false.
    end if

    if (resetRandom) then
       nread = nread+1
       DoWRT = .false.
       if (nRead.lt.3) return
       DoWRT = .true.

       write(*,*) '=== ATTENTION: reset of random generators!!! ==='

       if (nRead.gt.3) then
          write(*,*) 'nRead>3. stop.'
          stop
       endif

       if (.not.ReadPYR()) then
          write(*,*) 'ReadPYR failed. stop.'
          stop
       endif
    else
       call DumpPYR()  ! write the PYTHIA ran.gen. to file
    endif

    !    write(*,*) rn(0), pyr(0)

  end subroutine setRandom


  !*************************************************************************
  !****f* random/rn_openInterval
  ! NAME
  ! real function rn_openInterval()
  ! FUNCTION
  ! evaluates a random number in (0,1)
  ! USAGE 
  ! (real)=rn_openInterval()
  ! Notes
  ! Finds random number which is in [1E-8,1-1E-8], which is approximately (0,1).
  !*************************************************************************
  real function rn_openInterval()
    implicit none
    real,parameter :: eps=1E-8
    do
       ! Find random number which is not 1 or 0:
       rn_openInterval=rn()
       if(abs(rn_openInterval-1.).gt.eps.and.abs(rn_openInterval).gt.eps) return
    end do
  end function rn_openInterval


  !*************************************************************************
  !****f* random/rn
  ! NAME
  ! real function rn()
  ! FUNCTION
  ! evaluates a random number in (0,1)
  ! USAGE 
  ! (real)=rn()
  ! NOTES
  ! ckecks whether the random generator should be initialised, otherwise
  ! it just calls PYR. Endpoints are excluded.
  !*************************************************************************
  real function rn()
    implicit none
    real PYR ! prototype

    if (first) then
       ! Initializes the random number generator 
       ! when it's called for the first time
       call InitRan
       first=.false.
    end if

    rn=PYR(0)
    return
  end function rn

  !*************************************************************************
  !****f* random/rnGauss
  ! NAME
  ! real function rnGauss(StdDev, Mean)
  ! FUNCTION
  ! evaluates a random number in according a Gauss distribution
  ! INPUTS
  ! * real :: StdDev -- standard deviation
  ! * real :: Mean   -- mean value of distribution
  ! RESULT
  ! random number
  !*************************************************************************
  real FUNCTION rnGauss(StdDev, Mean)
    IMPLICIT NONE
    real StdDev, Mean
    
    common /rnGauss_Dat/ x1,x2,isCached
    real x1,x2
    integer isCached
    save /rnGauss_Dat/
    
    data isCached /0/

    real v1,v2, w, y

    if (isCached.ne.0) then
       isCached = 0
       rnGauss = x2*StdDev+Mean
       return
    endif

10  v1 = 2.0*rn() - 1.0 
    v2 = 2.0*rn() - 1.0
    w = v1**2 + v2**2
    if (w.gt.1.0) goto 10

    y = sqrt( (-2.0*log(w))/w )
    x1 = v1 * y
    x2 = v2 * y

    isCached = 1
    rnGauss = x1*StdDev+Mean
    return
  end FUNCTION rnGauss


  !*************************************************************************
  !****f* random/rnCos
  ! NAME
  ! real function rnCos()
  ! FUNCTION
  ! Evaluates a random number x according to cos(x) distribution. 
  ! So cos(x) is assumed to be isotropic in [-1,1] 
  ! INPUTS
  ! * NONE
  ! RESULT
  ! random number
  !*************************************************************************
  FUNCTION rnCos() Result(x)
    IMPLICIT NONE
    real :: x
    real :: cosX
    cosX=1.-2.*rn()
    x=acos(cosX)
  end FUNCTION rnCos

  !*************************************************************************
  ! cf. interface random/rnExp
  !*************************************************************************
  FUNCTION rnExp1(A) Result(x)
    IMPLICIT NONE
    real, intent(in) :: A
    real :: x
    x= -Log(rn())/A
  end FUNCTION rnExp1
  !-------------------------------------------------------------------------
  FUNCTION rnExp2(A,x1,x2) Result(x)
    IMPLICIT NONE
    real, intent(in) :: A,x1,x2
    real :: x,r1,r2
    r1 = exp(-A*x1)
    r2 = exp(-A*x2)
    x = -log(r2+rn()*(r1-r2))/A
  end FUNCTION rnExp2

  !*************************************************************************



end module random

