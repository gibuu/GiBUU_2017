program test
  use distributions
  implicit none
  real :: pole=938., width=1., x
  integer :: i

  ! write out Breit-Wigner distribution
  do i=90000,100000
     x=float(i)*0.01
     write(14,*) x, bW(x,pole,width)
  end do

  !write(*,*) bW_norm(pole,width,pole-1.,pole+1.)
  !write(*,*) bW_norm(pole,width,pole-10.,pole+10.)
  !write(*,*) bW_norm(pole,width,pole-100.,pole+100.)

  ! write out Gauss and Novosibirsk distributions
  do i=1,1000
    x=float(i)*0.001
    write(15,*) x, gauss(x,0.782,0.025), novo(x,0.782,0.025,-0.09)
  end do

end program test
