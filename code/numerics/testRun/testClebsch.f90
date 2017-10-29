program TestClebsch

use ClebschGordan, only: CG

implicit none

real :: j1,j2,j,m1,m2,c1 !,c2

real, parameter :: jmax = 8.

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
           c1 = CG(Nint(2*j1),nint(2*j2),nint(2*j),nint(2*m1),nint(2*m2))
!            c2 = CG_gsl(nint(2*j1),nint(2*j2),nint(2*j),nint(2*m1),nint(2*m2))
           print '(A,2F8.4,A,2F8.4,A,2F8.4,A,2F8.4)','<', j1,m1,' ;',j2,m2,' |',j,m1+m2,' > =',c1 !,c2
!            if (abs(c1-c2)>1E-6) stop 'FAILURE'
           m2=m2+1
         end do
         m1=m1+1
       end do
       print *
       j=j+1.
     end do
    j2=j2+0.5
  end do
  j1=j1+0.5
end do

! contains
!
!   double precision function CG_gsl (j1, j2, j3, m1, m2, m3_in)
!     use iso_c_binding, only: c_int, c_double
!
!     integer, intent(in) :: j1,j2,j3,m1,m2
!     integer, intent(in), optional :: m3_in
!
!     interface
!       function gsl_sf_coupling_3j(two_ja, two_jb, two_jc, two_ma, two_mb, two_mc) bind(c)
!       import
!       integer(c_int), value :: two_ja, two_jb, two_jc, two_ma, two_mb, two_mc
!       real(c_double) :: gsl_sf_coupling_3j
!       end function gsl_sf_coupling_3j
!     end interface
!
!     integer :: m3
!     if (Present(m3_in)) then
!       m3=m3_in
!     else
!       m3=m1+m2
!     end if
!     if (j1<0 .or. j2<0 .or. j3<0 .or. abs(m1)>j1 .or. abs(m2)>j2 .or. abs(m3)>j3 .or. &
!         j3>j1+j2 .or. j3<abs(j1-j2) .or. m3/=m1+m2 .or. &
!         mod(j1+j2+j3,2)/=0 .or. mod(j1+m1,2)/=0 .or. mod(j2+m2,2)/=0 .or. mod(j3+m3,2)/=0) then
!       CG_gsl = 0.
!     else if (j1==0 .or. j2==0) then
!       CG_gsl = 1.
!     else
!       CG_gsl = (-1.d0)**((j1-j2+m3)/2) * dsqrt(dble(j3+1)) * gsl_sf_coupling_3j(j1,j2,j3,m1,m2,-m3)
!     end if
!
!   end function

end program
