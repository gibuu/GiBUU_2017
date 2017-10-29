! This is a benchmark for the Clebsch-Gordan routine, which compares it to the
! GSL implementation (gsl_sf_coupling_3j). It needs to be linked against the GSL,
! i.e. adding "-lgsl -lgslcblas" in the Makefile is required.
program BenchClebsch

use ClebschGordan, only: CG

implicit none

integer, parameter ::  N = 1E4
real, parameter :: JJJ = 2.5

real :: j1,j2,j,m1,m2,c,jmax
integer :: i
real :: clock_start, clock_finish, t1, t2

jmax=0.
open (42, file="bench.dat")

do while (jmax<=JJJ)
  print '(A,f3.1)',"jmax=",jmax

  call cpu_time(clock_start)
  do i=1,N
    j1=0.
    do while (j1<=jmax)
      j2=0.
      do while (j2<=jmax)
         j=abs(j1-j2)
         do while (j<=j1+j2)
            m1=-j1
            do while (m1<=j1)
               m2=-j2
               do while (m2<=j2)
                  c = CG_gsl(nint(j1*2),nint(j2*2),nint(j*2),nint(2*m1),nint(2*m2))
                  m2=m2+1
               end do
               m1=m1+1
            end do
            j=j+1.
         end do
         j2=j2+0.5
      end do
      j1=j1+0.5
    end do
  end do
  call cpu_time(clock_finish)
  t1 = clock_finish-clock_start
  print *,"CG_gsl : ",time_format(t1)

  call cpu_time(clock_start)
  do i=1,N
    j1=0.
    do while (j1<=jmax)
      j2=0.
      do while (j2<=jmax)
         j=abs(j1-j2)
         do while (j<=j1+j2)
            m1=-j1
            do while (m1<=j1)
               m2=-j2
               do while (m2<=j2)
                  c = CG(nint(j1*2),nint(j2*2),nint(j*2),nint(2*m1),nint(2*m2))
                  m2=m2+1
               end do
               m1=m1+1
            end do
            j=j+1.
         end do
         j2=j2+0.5
      end do
      j1=j1+0.5
    end do
  end do
  call cpu_time(clock_finish)
  t2 = clock_finish-clock_start
  print *,"CG     : ",time_format(t2)

  print *
  write(42,*) jmax,t1,t2
  flush(42)
  jmax=jmax+0.5
end do

close(42)


contains


  double precision function CG_gsl (j1, j2, j3, m1, m2, m3_in)
    use iso_c_binding, only: c_int, c_double

    integer, intent(in) :: j1,j2,j3,m1,m2
    integer, intent(in), optional :: m3_in

    interface
      real(c_double) function gsl_sf_coupling_3j (two_ja, two_jb, two_jc, two_ma, two_mb, two_mc) bind(c)
        import
        integer(c_int), value :: two_ja, two_jb, two_jc, two_ma, two_mb, two_mc
      end function
    end interface

    integer :: m3
    if (Present(m3_in)) then
      m3=m3_in
    else
      m3=m1+m2
    end if
    if (j1<0 .or. j2<0 .or. j3<0 .or. abs(m1)>j1 .or. abs(m2)>j2 .or. abs(m3)>j3 .or. &
        j3>j1+j2 .or. j3<abs(j1-j2) .or. m3/=m1+m2 .or. &
        mod(j1+j2+j3,2)/=0 .or. mod(j1+m1,2)/=0 .or. mod(j2+m2,2)/=0 .or. mod(j3+m3,2)/=0) then
      CG_gsl = 0.
    else if (j1==0 .or. j2==0) then
      CG_gsl = 1.
    else
      CG_gsl = (-1.d0)**((j1-j2+m3)/2) * dsqrt(dble(j3+1)) * gsl_sf_coupling_3j(j1,j2,j3,m1,m2,-m3)
    end if

  end function


  function time_format(t)
    real,intent(in)::t
    character::time_format*13
    integer:: ms,s,m,h
    ms=mod(int(t*1000),1000)
    s=mod(int(t),60)
    m=mod(int(t/60),60)
    h=mod(int(t/60/60),24)
    write(time_format,"(i2.2,a,i2.2,a,i2.2,a,i3.3)") h,":",m,":",s,".",ms
  end function


end program
