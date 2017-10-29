!******************************************************************************
!****m* /clebschGordan
! NAME
! module clebschGordan
!
! PURPOSE
! Provides Clebsch-Gordan factors
!******************************************************************************
module clebschGordan

  implicit none
  private

  public :: clebschSquared, CG

  logical, parameter :: debug = .false.

contains


  !****************************************************************************
  !****f* clebschGordan/clebschSquared
  ! NAME
  ! real function clebschSquared(is1, is2, is3, iz1, iz2, iz3)
  !
  ! PURPOSE
  ! This function determines the isospin-factors for the processes  1 2 -> 3,
  ! namely squares of the Clebsch-Gordan factors
  ! <is1, iz1; is2, iz2 | is3, iz3 >**2.
  ! The argument "iz3" is optional.
  ! If it is omitted, it's assumed to be iz3=iz1+iz2.
  !
  ! INPUTS
  ! * real :: is1,is2,is3 -- total isospins
  ! * real :: iz1,iz2,iz3 -- z-component of the isospins
  ! OUTPUT
  ! * real  --  Square of the clebsch gordan factor
  !****************************************************************************
  real function clebschSquared(is1, is2, is3, iz1, iz2, iz3_in)

    real, intent(in) :: is1,is2,is3,iz1,iz2
    real, intent(in), optional :: iz3_in
    real :: iz3
    real, parameter :: epsilon=0.001

    if (present(iz3_in)) then
      iz3=iz3_in
    else
      iz3=iz1+iz2
    end if
    if (debug) then
      ! Check input. Input must be of the form : j or j+0.5 with j being an integer,
      ! and is1,is2,is3 must be positive
      if (.not.(isInteger(2*is1)).or.(is1.lt.-epsilon)) then
         write(*,*) 'is1 is bad input',is1
         stop
      end if
      if (.not.(isInteger(2*is2)).or.(is2.lt.-epsilon)) then
         write(*,*) 'is2 is bad input',is2
         stop
      end if
      if (.not.(isInteger(2*is3)).or.(is3.lt.-epsilon)) then
         write(*,*) 'is3 is bad input',is3
         stop
      end if
      if (.not.(isInteger(2*iz1))) then
         write(*,*) 'iz1 is bad input',iz1
         stop
      end if
      if (.not.(isInteger(2*iz2))) then
         write(*,*) 'iz2 is bad input',iz2
         stop
      end if
      if (.not.(isInteger(2*iz3))) then
         write(*,*) 'iz3 is bad input',iz3
         stop
      end if

    end if

    clebschSquared = CG(nint(is1*2),nint(is2*2),nint(is3*2),nint(iz1*2),nint(iz2*2),nint(iz3*2))**2

  contains

    logical function isInteger(x)
      real, intent(in) ::x
      isInteger=(abs(x-Int(x))<epsilon)
    end function isInteger

  end function clebschSquared



  !****************************************************************************
  !****f* clebschGordan/CG
  ! NAME
  ! real function CG(j1, j2, j3, m1, m2, m3)
  !
  ! PURPOSE
  ! A function to compute arbitrary Clebsch-Gordan coefficients.
  ! Works for isospin values of up to j=50.
  ! The argument "m3" is optional.
  ! If it is omitted, it's assumed to be m3=m1+m2.
  !
  ! Input values are 'isospin times 2' to have integer values
  !
  ! INPUTS
  ! * integer :: j1,j2,j3 -- isospin (times 2)
  ! * integer :: m1,m2,m3 -- z-component of the isospin (times 2)
  ! OUTPUT
  ! * real  --  Clebsch-Gordan coefficient <j1,m1;j2,m2|j3,m3>
  !
  ! NOTE
  ! This code is based on the original 'f3j' routine by Caswell & Maximon,
  ! cf. http://archive.org/details/fortranprogramsf409casw
  !****************************************************************************
  double precision function CG(j1, j2, j3, m1, m2, m3_in)
    integer, intent(in) :: j1,j2,j3,m1,m2
    integer, intent(in), optional :: m3_in
    integer :: m3
    if (present(m3_in)) then
       m3=m3_in
       if (m3/=m1+m2) then
          CG = 0.
          return
       end if
    else
       m3=m1+m2
    end if
    if (j1<0 .or. j2<0 .or. j3<0 &
         .or. abs(m1)>j1 .or. abs(m2)>j2 .or. abs(m3)>j3 &
         .or. j3>j1+j2 .or. j3<abs(j1-j2) &
         .or. mod(j1+j2+j3,2)/=0 .or. mod(j1+m1,2)/=0 &
         .or. mod(j2+m2,2)/=0 .or. mod(j3+m3,2)/=0) then
      CG = 0.
    else if (j1==0 .or. j2==0) then
      CG = 1.
    else
      CG = (-1.d0)**((j1-j2+m3)/2) * dsqrt(dble(j3+1)) * f3j(j1,j2,j3,m1,m2,-m3)
    end if

  contains

    double precision function f3j(j1, j2, j3, m1, m2, m3)

      integer, intent(in) :: j1,j2,j3,m1,m2,m3
      integer :: mtri(9),n,kmin,kmax,mini,min2,min3,min4,min5,ncut,k,num
      double precision :: fn,uk,s,delog,ulog,plog,sig,slog,p

      logical, save :: init = .true.
      integer, parameter :: fl_max = 4*50 + 2 ! originally: 322
      ! with the above array size, the routine will work for all spin values up to 50
      double precision, save :: fl(fl_max) = 0.0d0

      f3j=0.0d0

      ! initialization: calculate fl(n)
      if (init) then
        fn=1.0d0
        do n=3,fl_max
          fn=fn+1.0d0
          fl(n)=fl(n-1)+log(fn)
        end do
        init = .false.
      end if

      mtri(1)=(j1+j2-j3)/2
      mtri(2)=(j1-j2+j3)/2
      mtri(3)=(-j1+j2+j3)/2
      mtri(4)=(j1+m1)/2
      mtri(5)=(j1-m1)/2
      mtri(6)=(j2+m2)/2
      mtri(7)=(j2-m2)/2
      mtri(8)=(j3+m3)/2
      mtri(9)=(j3-m3)/2

      kmin=max(-j3+j2-m1,-j3+j1+m2,0)/2
      if (j2-j3+m1<0) then
        kmax=j1+j2-j3
      else
        kmax=j1-m1
      end if
      kmax=min(j2+m2,kmax)/2

      mini=mtri(1)-kmin+1
      min2=mtri(5)-kmin+1
      min3=mtri(6)-kmin+1
      min4=(j3-j2+m1)/2+kmin
      min5=(j3-j1-m2)/2+kmin

      ! sum series in double precision
      uk=1.0d-10
      s=1.0d-10
      ncut=0
      kmax=kmax-kmin
      if (kmax>0) then
        do k=1,kmax
          uk=-uk*dble((mini-k)*(min2-k)*(min3-k))/dble((kmin+k)*(min4+k)*(min5+k))
          if (abs(uk)>=1.0d30) then
            uk=1.0d-10*uk
            s=1.0d-10*s
            ncut=ncut+1
          end if
          if (abs(uk)<1.0d-20) exit
          s=s+uk
        end do
      end if

      ! calculate delta functions
      delog=0.0d0
      do n=1,9
        delog=delog+fl(mtri(n)+1)
      end do
      num=(j1+j2+j3)/2+2
      delog=0.5d0*(delog-fl(num))
      ulog=-fl(kmin+1)-fl(mini)-fl(min2)-fl(min3)-fl(min4+1)-fl(min5+1)
      plog=delog+ulog
      if (plog+80.0d0<0 .or. ncut>0) then
        sig=sign(1.0d0,s)
        s=abs(s)
        slog=log(s)+dble(ncut+1)*log(1.0d+10)
        f3j=sig*exp(slog+plog)
      else
        s=s*1.0d+10
        p=exp(plog)
        f3j=p*s
      end if
      num=kmin+(j1-j2-m3)/2
      if (mod(num,2)/=0) f3j=-f3j

    end function f3j

  end function CG



end module clebschGordan
