program test

call sub1
call sub2
call sub3
call sub4

end program



subroutine sub1
  use tabulation
  implicit none

  type(table_1dim) :: t1
  integer :: i,check
  real ::x,y
  real :: min=-2.
  real :: max=12.
  real :: delta=0.1

  open(11,file="test1.dat")
  if(.not.read_1dim(t1,"t1.dat",min,max,delta)) then 
     t1=init_table_1dim(min,max,delta,f)
  end if
  call printTable(t1,'t1.dat')
  do i=-1000,1000
     x=float(i)*0.02
     y=getValue(x,t1,check)
     write(11,'(4E15.4,I3)') x,f(x)-y,f(x),y,check
  end do
  close(11)
       
contains 

  real function f(x)
    implicit none
    real :: x
    f=sin(x)
    

  end function f


end subroutine sub1


subroutine sub2
  use tabulation
  implicit none

  type(table_2dim) :: t2
  integer :: i,j,check
  real ::x,y,z
  real :: min=-2.
  real :: max=2.
  real :: delta=0.07

  open(11,file="test2.dat")
  if(.not.read_2dim(t2,"t2.dat",(/min,min+1./),(/max,max+1./),(/delta,2.*delta/))) then 
     t2=init_table_2dim((/min,min+1./),(/max,max+1./),(/delta,2.*delta/),f)
  end if
  call printTable(t2,'t2.dat')
  do i=-30,30
  do j=-30,30
     x=float(i)*0.1
     y=float(j)*0.1
     z=getValue((/x,y/),t2,check)
     write(11,'(5E15.4,I3)') x,y,f(x,y)-z,f(x,y),z,check
  end do
  write(11,*)
  end do
  close(11)

       
contains 

  real function f(x,y)
    implicit none
    real :: x,y

    f=sin(x)+y
    

  end function f


end subroutine sub2



subroutine sub3
  use tabulation
  implicit none

  type(table_3dim) :: t3
  integer :: i,j,check
  real ::x,y,z
  real :: min=-2.
  real :: max=2.
  real :: delta=0.07

  open(11,file="test3.dat")
  if(.not.read_3dim(t3,"t3.dat",(/min,min+1.,min+1./),(/max,max+1.,max+1./),(/delta,2.*delta,2.*delta/))) then 
     t3=init_table_3dim((/min,min+1.,min+1./),(/max,max+1.,max+1./),(/delta,2.*delta,2.*delta/),f)
  end if
  call printTable(t3,'t3.dat')
  do i=-30,30
  do j=-30,30
     x=float(i)*0.1
     y=float(j)*0.1
     z=getValue((/x,y,1.3/),t3,check)
     write(11,'(5E15.4,I3)') x,y,f(x,y,1.3)-z,f(x,y,1.3),z,check
  end do
  write(11,*)
  end do
  close(11)

       
contains 

  real function f(x,y,z)
    implicit none
    real :: x,y,z

    f=sin(x)+y*z
    

  end function f


end subroutine sub3


subroutine sub4
  use tabulation
  implicit none

  type(table_4dim) :: t4
  integer :: i,j,check
  real ::x,y,z
  real :: min=-2.
  real :: max=2.
  real :: delta=0.07

  open(11,file="test4_1.dat")
  open(12,file="test4_2.dat")
  open(13,file="test4_3.dat")
  if(.not.read_4dim(t4,"t4.dat",(/min,min+1.,min+1.,min+1./),(/max,max+1.,max+1.,max+1./),(/delta,2.*delta,2.*delta,2.*delta/))) then 
!     t4=init_table_4dim((/min,min+1.,min+1.,min+1./),(/max,max+1.,max+1.,max+1./),(/delta,2.*delta,2.*delta,2.*delta/),f)
     t4=init_table((/min,min+1.,min+1.,min+1./),(/max,max+1.,max+1.,max+1./),(/delta,2.*delta,2.*delta,2.*delta/),f)
  end if
  call printTable(t4,'t4.dat')
  do i=-30,30
  do j=-30,30
     x=float(i)*0.1
     y=float(j)*0.1
     z=getValue((/x,y,1.3,-1./),t4,check)
     write(11,'(5E15.4,I3)') x,y,f(x,y,1.3,-1.)-z,f(x,y,1.3,-1.),z,check
     z=getValue((/x,1.3,y,-1./),t4,check)
     write(12,'(5E15.4,I3)') x,y,f(x,1.3,y,-1.)-z,f(x,1.3,y,-1.),z,check
     z=getValue((/x,1.3,-1.,y/),t4,check)
     write(13,'(5E15.4,I3)') x,y,f(x,1.3,-1.,y)-z,f(x,1.3,-1.,y),z,check
  end do
  write(11,*)
  write(12,*)
  write(13,*)
  end do
  close(11)
  close(12)
  close(13)

       
contains 

  real function f(x,y,z,a)
    implicit none
    real :: x,y,z,a

    f=sin(x)+y*z+3*a
    

  end function f


end subroutine sub4
