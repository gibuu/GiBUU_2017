
program test
  use idTable
  use inputGeneral
  use particleProperties
  implicit none

  write(*,*) "**********************************" 
  write(*,*) "**********************************" 
  write(*,*)
  Write(*,*)      "Testing the routines which generate the width of the resonances "
  write(*,*)
  write(*,*) "**********************************" 
  write(*,*) "**********************************" 
  write(*,*)
  call initParticleProperties()
  call readinputGeneral
  call testDelta_ticket71


contains


  subroutine testDelta_ticket71
   use baryonWidth
   use deltaWidth
   implicit none
   integer :: i, id,j
   real :: mass,sum,dens,mom1(1:3),mom2(1:3)
   integer :: index
   mom1=(/0.,0.,1./)
   mom2=(/0.,0.,0./)
   dens=0.

   Open(999,file='DeltaWidth_vacuum.dat')
  
   massLoop: do i=0,200
     mass=0.8+float(i)*0.01
     write(999,'(6F15.9)') mass, &
          & delta_fullWidth(mass,mom2,dens),delta_fullWidth(mass,mom1,dens)&
          & , FullWidthBaryon(2,mass)
   end do massLoop
   close(999)

  end subroutine testDelta_ticket71
end program test
