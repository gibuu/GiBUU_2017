program TryFluxMC

  use random
  use histf90

  implicit none

  type(histogram) :: H1,H2
  

  real :: x,y,II,JJ

  integer :: i,n

  call CreateHist(H1,"f(x)", 0.0,1.0,0.01)
  call CreateHist(H2,"x according f", 0.0,1.0,0.01)
  
  n = 1000000

  II = 0.
  do i=1,n
     x = rn()
     y = f(x)
     
     call AddHist(H1,x,y)

     II = II + g(x)*y

  end do

  write(*,*) n,II/n
  call WriteHist(H1,141,mul=1./n, add=1e-20)

  II = 0.
  JJ = 0.
  do i=1,n
     x = selectAccordF()
     call AddHist(H2,x,1d0/f(x))

     JJ = JJ + 1d0/f(x)
     II = II + g(x)

  enddo
  write(*,*) n,II/n*5./3.,II/n
  call WriteHist(H2,142,mul=1./n, add=1e-20)
  write(*,*) JJ/n*5./3.
  write(*,*) n,II/JJ

contains

  real function f(x)
    implicit none
    real x
    
    f = 5 * x**2
  end function f

  real function fInv(y)
    implicit none
    real y
    
    fInv = sqrt(y/5.0)
  end function fInv

  real function g(x)
    implicit none
    real x
    
    g = sin(x*3.141592654)
  end function g

  real function selectAccordF()
    implicit none

    real :: x,ymax
    
    ymax = 5.0

    do
       x = rn()
       if (rn()*ymax.lt.f(x)) exit
    end do
    selectAccordF = x

  end function selectAccordF


end program TryFluxMC
