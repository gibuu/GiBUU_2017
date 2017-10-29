
program testFF
  use particleProperties, only: initParticleProperties
  implicit none

  call initParticleProperties
  call test_formfactor
!   call test_gammaWidth

contains

subroutine test_formfactor()
use FF_Delta_production
use Formfactor_ResProd
use particleproperties
use idtable

real, parameter :: dq=0.005
real, parameter :: dw=0.1
integer :: i  !,j
real :: Qs  !,W,c3v,c4v,c5v,c6v,c3a,c4a,c5a,c6a
!logical :: ff_set



open(10,file='deltaFF.dat')
open(11,file='deltaFF_neu.dat')
open(12,file='D13_FF.dat')

!do j=0,10
   !W=dw*j+1.1
   do i=0,400
      Qs=dq*i
    !  call formfactors_Delta(Qs,W,0,c3v,c4v,c5v,c6v,c3a,c4a,c5a,c6a)
    !  write(10,'(10F12.5)') Qs,W,c3v,c4v,c5v,c6v,c3a,c4a,c5a,c6a
!      write(11,'(10F12.5)') Qs, baryon(delta)%mass, getFormfactor_Res(Qs,baryon(delta)%mass,delta,1,0,FF_set)
!      write(12,'(10F12.5)') Qs, baryon(D13_1520)%mass, getFormfactor_Res(Qs,baryon(D13_1520)%mass,D13_1520,1,0,FF_set)
   end do
!end do

write(*,*) 'end of test FF delta'

end subroutine test_formfactor

! subroutine test_gammaWidth()
!   use Formfactor_ResProd
!   use idtable, only : delta
! 
!   integer :: i
!   real :: W,k_w,k_r, width
! 
!   do i=0,1000
!      W=1.+float(i)*0.01
!      width=width_gammaNR(delta,W,k_w,k_r)
!      write(100,'(4E12.4)') W, width, k_W, k_r
!   end do
! 
! end subroutine test_gammaWidth


end program testFF
