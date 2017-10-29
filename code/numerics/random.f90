!******************************************************************************
!****m* /random
! NAME
! module random
! PURPOSE
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
!******************************************************************************
module random
  implicit none
  private


  !****************************************************************************
  !****g* random/Seed
  ! PURPOSE
  ! Random Seed (used to initialize the random number generator),
  ! accessible through the namelist 'initRandom'.
  ! If Seed is zero (default), then it is set via "SYSTEM_CLOCK()".
  !
  ! SOURCE
  !
  integer,save ::  Seed = 0
  !****************************************************************************


  !****************************************************************************
  !****g* random/resetRandom
  ! PURPOSE
  ! Reread random generator, used by setRandom, useful for debugging.
  !
  ! SOURCE
  !
  logical,save :: resetRandom = .false.
  !****************************************************************************


  !****************************************************************************
  !****f* random/rnExp
  ! NAME
  ! real function rnExp (A)
  ! real function rnExp (A, x1, x2)
  ! PURPOSE
  ! Evaluates a random number x according to an exp(A*x) distribution.
  ! If given, x is restricted to lie between x1 and x2 (it doesn't matter which
  ! one of the two is larger).
  ! INPUTS
  ! * real :: A -- slope parameter (can be negative)
  ! * real :: x1 -- minimal x-value
  ! * real :: x2 -- maximal x-value
  ! RESULT
  ! random number
  !****************************************************************************
  interface rnExp
     module procedure rnExp1, rnExp2
  end interface

  public :: rn, rnFlat, rnExp, rnGauss, rnCos, rn_openInterval, rnPower
  public :: rnOmega, rnOmega_anis, rnOmega_angles
  public :: setRandom, rn_truefalse, ranCharge

  logical, save :: first = .true.

contains


  !****************************************************************************
  !****f* random/InitRan
  ! NAME
  ! subroutine InitRan
  ! PURPOSE
  ! Reads random seed out of namelist 'initRandom' in jobcard and initializes
  ! the random number generator.
  !****************************************************************************
  subroutine InitRan
    use output, only: Write_ReadingInput

    integer :: ios
    integer(8) :: seed8, rate

    !**************************************************************************
    !****n* random/initRandom
    ! NAME
    ! NAMELIST /initRandom/
    ! PURPOSE
    ! Includes the input variables:
    ! * Seed
    ! * resetRandom
    !**************************************************************************
    NAMELIST /initRandom/ Seed,resetRandom

    call Write_ReadingInput('initRandom',0)
    rewind(5)
    read(5,nml=initRandom,IOSTAT=ios)
    call Write_ReadingInput('initRandom',0,ios)

    if (Seed == 0) then
      call system_clock (seed8, rate)              ! get 8-byte clock value
      write(*,*) 'Resetting Seed via system clock: ', seed8
      write(*,*) '                      precision: ',rate,'counts/sec'
      seed = int(mod(seed8,int(huge(seed),8)),4)   ! convert to 4-byte seed
    end if

    write(*,*) 'Seed: ', Seed
    write(*,*) 'resetRandom: ', resetRandom

    call InitPYR (Seed)

    call Write_ReadingInput('initRandom',1)

    first = .false.

  end subroutine InitRan


  !****************************************************************************
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
  !****************************************************************************
  subroutine setRandom

    logical,external :: ReadPYR
    integer,save :: nRead = 0

    common/DebugWRT/ DoWRT
    logical DoWRT
    save/DebugWRT/

    if (first) call InitRan

    if (resetRandom) then
       nread = nread+1
       DoWRT = .false.
       if (nRead.lt.3) return
       DoWRT = .true.

       write(*,*) '=== ATTENTION: reset of random generators!!! ==='

       if (nRead.gt.3) then
          write(*,*) 'nRead>3. stop.'
          stop
       end if

       if (.not.ReadPYR()) then
          write(*,*) 'ReadPYR failed. stop.'
          stop
       end if
    else
       call DumpPYR()  ! write the PYTHIA ran.gen. to file
    end if

    !    write(*,*) rn(0), pyr(0)

  end subroutine setRandom


  !****************************************************************************
  !****f* random/rn_openInterval
  ! NAME
  ! real function rn_openInterval()
  ! PURPOSE
  ! Evaluates a random number in (0,1).
  ! USAGE
  ! (real)=rn_openInterval()
  ! Notes
  ! Finds random number which is in [1E-8,1-1E-8], which is approximately (0,1).
  !****************************************************************************
  real function rn_openInterval()
    real,parameter :: eps=1E-8
    do
       ! Find random number which is not 1 or 0:
       rn_openInterval=rn()
       if (abs(rn_openInterval-1.).gt.eps.and.abs(rn_openInterval).gt.eps) return
    end do
  end function rn_openInterval


  !****************************************************************************
  !****f* random/rn
  ! NAME
  ! real function rn()
  ! PURPOSE
  ! Evaluates a random number in (0,1).
  ! USAGE
  ! (real)=rn()
  ! NOTES
  ! Checks whether the random generator should be initialised, otherwise
  ! it just calls PYR. Endpoints are excluded.
  !****************************************************************************
  real function rn()
    real,external :: PYR ! prototype
    if (first) call InitRan
    rn=PYR(0)
  end function rn


  !****************************************************************************
  !****f* random/rnFlat
  ! NAME
  ! real function rnFlat(xMin,xMax)
  ! PURPOSE
  ! Draws a random number from a flat (uniform) distribution on (xMin,xMax).
  ! INPUTS
  ! * real :: xMin,xMax -- minimum and maximum value
  !****************************************************************************
  real function rnFlat (xMin,xMax)
    real, intent(in) :: xMin, xMax
    rnFlat = XMin + rn()*(xMax-xMin)
  end function rnFlat


  !****************************************************************************
  !****f* random/rnGauss
  ! NAME
  ! real function rnGauss(StdDev, Mean)
  ! PURPOSE
  ! evaluates a random number in according a Gauss distribution
  ! INPUTS
  ! * real :: StdDev -- standard deviation
  ! * real :: Mean   -- mean value of distribution
  ! RESULT
  ! Random number.
  !****************************************************************************
  real FUNCTION rnGauss(StdDev, Mean)
    real,intent(in) :: StdDev, Mean
    real,save :: x1,x2
    integer,save :: isCached = 0
    real v1,v2, w, y

    if (isCached.ne.0) then
       isCached = 0
       rnGauss = x2*StdDev+Mean
       return
    end if

    do
      v1 = 2.0*rn() - 1.0
      v2 = 2.0*rn() - 1.0
      w = v1**2 + v2**2
      if (w.le.1.0) exit
    end do

    y = sqrt( (-2.0*log(w))/w )
    x1 = v1 * y
    x2 = v2 * y

    isCached = 1
    rnGauss = x1*StdDev+Mean
  end FUNCTION rnGauss


  !****************************************************************************
  !****f* random/rnCos
  ! NAME
  ! real function rnCos()
  ! PURPOSE
  ! Evaluates a random number x according to cos(x) distribution.
  ! So cos(x) is assumed to be isotropic in [-1,1]
  ! INPUTS
  ! * NONE
  ! RESULT
  ! Random number.
  !****************************************************************************
  FUNCTION rnCos() Result(x)
    real :: x
    real :: cosX
    cosX=1.-2.*rn()
    x=acos(cosX)
  end FUNCTION rnCos


  !****************************************************************************
  ! cf. interface random/rnExp
  !****************************************************************************
  real FUNCTION rnExp1 (A)
    real, intent(in) :: A
    rnExp1 = Log(rn())/A
  end FUNCTION rnExp1
  !-------------------------------------------------------------------------
  real FUNCTION rnExp2 (A, x1, x2)
    real, intent(in) :: A, x1, x2
    real :: r1,r2
    r1 = exp(A*x1)
    r2 = exp(A*x2)
    rnExp2 = log(r2+rn()*(r1-r2))/A
  end FUNCTION rnExp2
  !****************************************************************************


  !****************************************************************************
  !****f* random/rnPower
  ! NAME
  ! real function rnPower (n, xMin, xMax)
  ! PURPOSE
  ! Draws a random number according to a power-law distribution ~ x^n.
  ! INPUTS
  ! * real :: n         -- exponent in power law (arbitrary real number)
  ! * real :: xMin,xMax -- minimum and maximum value
  ! RESULT
  ! Random number between xMin and xMax.
  !****************************************************************************
  real function rnPower (n, xMin, xMax)
    real, intent(in) :: n, xMin, xMax
    real :: n1

    n1 = n + 1.
    if (abs(n1)<1E-3) then
      rnPower = xMin * (xMax/xMin)**rn()
    else if (xMin>0. .and. xMax>0.) then
      rnPower = (xMin**n1 + rn()*(xMax**n1-xMin**n1)) ** (1./n1)
    else
      rnPower = sign(1.,xMin) * (abs(xMin)**n1 + rn()*(abs(xMax)**n1-abs(xMin)**n1)) ** (1./n1)
    end if

  end function


  !****************************************************************************
  !****f* random/rnOmega
  ! NAME
  ! function rnOmega()
  ! PURPOSE
  ! Generates a unit vector with random angle (isotropical distribution).
  ! INPUTS
  ! * NONE
  ! RESULT
  ! real :: rnOmega(3) - random 3-dim. unit-vector
  !****************************************************************************
  function rnOmega()
    use constants, only: twopi
    real:: rnOmega(3)
    real:: phi,cost,sint
    phi=twopi*rn()                      ! phi
    cost=2.*rn()-1.                     ! cos(theta)
    sint=sqrt(1.-cost**2)               ! sin(theta)
    rnOmega(1:3) = (/ sint*cos(phi), sint*sin(phi), cost /)
  end function rnOmega


  !****************************************************************************
  !****f* random/rnOmega_anis
  ! NAME
  ! function rnOmega_anis (B)
  ! PURPOSE
  ! Generates a unit vector with random angle (anisotropical distribution).
  ! cos(theta) is chosen according to a probability distribution ~ 1 + B*cos**2(theta).
  ! INPUTS
  ! * real, intent(in) :: B  ---  parameter of prob. distr.
  ! RESULT
  ! real :: rnOmega_anis(3) - random 3-dim. unit-vector
  !****************************************************************************
  function rnOmega_anis (B)
    use constants, only: twopi
    real :: rnOmega_anis(3)                   ! return value: 3-vector
    real, intent(in) :: B                     ! parameter of prob. distr.
    real :: phi,cost,sint,x

    phi=twopi*rn()

    ! Choose cos(theta) according to f(cost)~(1+B*cost**2)
    do
       cost = 1 - 2.*rn()
       x = rn()
       if (x < (1+B*cost**2)/(2*(1+B/3.))) exit
    end do

    sint = sqrt(1.-cost**2)   ! sin(theta)
    rnOmega_anis(1:3) = (/ sint*cos(phi), sint*sin(phi), cost /)

  end function rnOmega_anis


  !****************************************************************************
  !****s* random/rnOmega_angles
  ! NAME
  ! subroutine rnOmega(theta,phi)
  ! PURPOSE
  ! Generates random angles theta and phi (istropical distribution) in degree.
  ! INPUTS
  ! * NONE
  ! OUTPUT
  ! * real :: theta, phi
  !****************************************************************************
  subroutine rnOmega_angles(theta,phi)
    use constants, only: pi
    real,intent(out) :: phi,theta
    real :: cost
    phi=rn()   *360.                   ! phi
    cost=2.*rn()-1.                     ! cos(theta)
    theta=acos(cost)*180./pi
  end subroutine rnOmega_angles
 !*****************************************************************************


  !****************************************************************************
  !****f* random/rn_trueFalse
  ! NAME
  ! function rn_trueFalse
  ! PURPOSE
  ! Gives true or false with same probability.
  ! INPUTS
  ! * NONE
  ! RESULT
  ! Logical.
  !****************************************************************************
  function rn_truefalse()
    logical:: rn_trueFalse
    rn_trueFalse=(rn().gt.0.5)
  end function rn_truefalse


  !****************************************************************************
  !****s* random/ranCharge
  ! NAME
  ! subroutine ranCharge (izmin, izmax, iztot, izout, flag)
  ! PURPOSE
  ! This subroutine distributes the charges randomly
  ! to a given total charge iztot.
  ! INPUTS
  ! * integer, intent(in)   :: izmin(:) ! Vector of minimal charges
  ! * integer, intent(in)   :: izmax(:) ! Vector of maximal charges
  ! * integer, intent(in)   :: iztot    ! total charge
  ! OUTPUT
  ! * real, dimension(:)    :: izout    ! Random charge configuration
  ! * logical, intent(out)  :: flag     ! =.false. if procedure failed
  !****************************************************************************
  subroutine rancharge (izmin, izmax, iztot, izout, flag)
    integer, intent(in)  :: izmin(:), izmax(:), iztot
    integer, intent(out) :: izout(:)
    logical, intent(out) :: flag

    integer :: n,totmax,totmin,numc,numf,i,i2,j,d,iztv
    integer, parameter :: numfm = 100
    integer :: ifo(numfm)

    ! Check Input
    if ((size(izout,dim=1)/=size(izmax,dim=1)) .or. (size(izout,dim=1)/=size(izmin,dim=1))) then
       write(*,*) 'Critical error in ranCharge, dimensions do not fit'
    else
       n=size(izout,dim=1)
    end if

    totmax=0
    totmin=0
    numc=1

    do i=1,n
       totmax=totmax+izmax(i)
       totmin=totmin+izmin(i)
       numc=numc*(izmax(i)-izmin(i)+1)
    end do

    if (totmin>iztot .or. totmax<iztot) then
       flag=.false.
       return
    else
       flag=.true.
    end if

    numf=0
    do i=0,numc-1
       i2=i
       iztv=0
       d=numc
       do j=1,n
          d=d/(izmax(j)-izmin(j)+1)
          izout(j)=izmin(j)+i2/d
          i2=i2-(i2/d)*d
          iztv=iztv+izout(j)
       end do
       !   *         write(*,*)i,'charges:',izout,numf
       if (iztv==iztot) then
          numf=numf+1
          if (numf>numfm) then
             write(*,*) 'problems in rancharge numf.gt.numfm'
             numf=numfm
          end if
          ifo(numf)=i
       end if
    end do
    ! *monte-carlo decision
    if (numf==0) then
       write(*,*) 'problems in rancharge numf=0'
       flag=.false.
       return
    end if
    i=ifo(int(rn()*numf)+1)
    i2=i
    d=numc
    do j=1,n
       d=d/(izmax(j)-izmin(j)+1)
       izout(j)=izmin(j)+i2/d
       i2=i2-(i2/d)*d
    end do
  end subroutine rancharge


end module random
