program testFF
use inputGeneral
use  formfactors_A
implicit none


integer :: i,j,k
real :: dT=2
real :: ds=0.002
real :: dQ=0.01
complex :: A(1:6)

call readinputGeneral()

! Theta const
do i=0,100
   do j=0,200
      A=getA(30.,(1.1+ds*i)**2,dQ*j)
      write(11,'(14E15.4)') (1.1+ds*i)**2,dQ*j, A(1:6)
   end do
   write(11,*)
end do

! s const
do i=0,90
   do j=0,200
      A=getA(dt*i,1.7,dQ*j)
      write(12,'(14E15.4)') dT*i,dQ*j, A(1:6)
   end do
   write(12,*)
end do

! q^2 const
do i=0,90
   do j=0,200
      A=getA(dt*i,(1.1+ds*j)**2,0.6)
      write(13,'(14E15.4)') dT*i,(1.1+ds*j)**2, A(1:6)
   end do
   write(13,*)
end do

! q^2 const
! W const
write(14,*) '# W=1.4, q^2=0.6'
do i=0,90
   A=getA(dt*i,(1.4)**2,0.6)
   write(14,'(15E11.3)') dt*i, A(1:6)
end do




end program testFF
