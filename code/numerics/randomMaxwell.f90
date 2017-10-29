!******************************************************************************
!****m* /randomMaxwell
! NAME
! module randomMaxwell
! PURPOSE
! generate momentum for massive particles according a relativistic
! Maxwell distribution.
!
! Algorithm taken from M. Swisdak, Phys.Plasmas 20 (2013) 062110,
! http://inspirehep.net/record/1235671
!******************************************************************************
module randomMaxwell
  implicit none
  private

  public:: initMaxwell
  public:: rnMaxwell

  real :: A, pM, pM2, pLo, pHi
  real :: lambdaLo, lambdaHi
  real :: qM, qLo, qHi
  real :: fpM,logfpM
  real :: mass_

contains

  !****************************************************************************
  !****if* randomMaxwell/funk
  ! NAME
  ! real function funk(p)
  ! PURPOSE
  ! distribution function to sample from
  !****************************************************************************
  real function funk(p)
    real, intent(in) :: p
    real :: p2
    p2 = p**2
    funk = p2*exp(-A*p2/(1+sqrt(1+p2)))
  end function funk

  !****************************************************************************
  !****if* randomMaxwell/funkZero
  ! NAME
  ! real function funkZero(p)
  ! PURPOSE
  ! internal function for finding zeros
  !****************************************************************************
  real function funkZero(p)
    real, intent(in) :: p
    real :: p2
    p2 = p**2
    funkZero = log(p2)-A*p2/(1+sqrt(1+p2)) - logfpm + 1
  end function funkZero

  !****************************************************************************
  !****if* randomMaxwell/findLo
  ! NAME
  ! real function findLo(p)
  ! PURPOSE
  ! find the left solution funk(pLo)*exp(1) = funk(pM)
  !****************************************************************************
  real function findLo(p)
    real, intent(in) :: p
    real, dimension(1:3) :: x,y

    x = (/ p/2, p , 0.0 /)
    y(1) = funkZero(x(1))
    y(2) = funkZero(x(2))
    y(3) = 0.0

    do
       if (y(1)*y(2) < 0) exit ! zero is bracketed
       x(2) = x(1)
       x(1) = x(1)/2
       y(2) = y(1)
       y(1) = funkZero(x(1))
    end do

    do
       if (x(2)-x(1) < 1e-4) exit ! accuracy reached
       x(3) = (x(1)+x(2))/2
       y(3) = funkZero(x(3))

       if (y(3)*y(2) < 0) then
          x(1) = x(3)
          y(1) = y(3)
       else
          x(2) = x(3)
          y(2) = y(3)
       end if
    end do

    findLo = x(1)

  end function findLo

  !****************************************************************************
  !****if* randomMaxwell/findHi
  ! NAME
  ! real function findHi(p)
  ! PURPOSE
  ! find the right solution funk(pHi)*exp(1) = funk(pM)
  !****************************************************************************
  real function findHi(p)
    real, intent(in) :: p
    real, dimension(1:3) :: x,y

    x = (/ p, 2*p , 0.0 /)
    y(1) = funkZero(x(1))
    y(2) = funkZero(x(2))
    y(3) = 0.0

    do
       if (y(1)*y(2) < 0) exit ! zero is bracketed
       x(1) = x(2)
       x(2) = x(2)+p
       y(1) = y(2)
       y(2) = funkZero(x(2))
    end do

    do
       if (x(2)-x(1) < 1e-4) exit ! accuracy reached
       x(3) = (x(1)+x(2))/2
       y(3) = funkZero(x(3))

       if (y(3)*y(2) < 0) then
          x(1) = x(3)
          y(1) = y(3)
       else
          x(2) = x(3)
          y(2) = y(3)
       end if
    end do

    findHi = x(1)

  end function findHi

!!$  !*************************************************************************
!!$  subroutine PrintDist
!!$    integer :: i
!!$    real :: p,y
!!$
!!$    do i=1,100
!!$       p = i*0.1
!!$       y = funk(p)
!!$
!!$       write(113,*) p,y
!!$    end do
!!$
!!$    write(113,*) "# pLo: ",pLo,funk(pLo),funk(pLo)*exp(1.0)
!!$    write(113,*) "# pM:  ",pM,funk(pM)
!!$    write(113,*) "# pHi: ",pHi,funk(pHi),funk(pHi)*exp(1.0)
!!$
!!$    write(113,*) "# qLo: ",qLo
!!$    write(113,*) "# qM:  ",qM
!!$    write(113,*) "# qHi: ",qHi
!!$
!!$    write(113,*) "# lambdaLo: ",lambdaLo
!!$    write(113,*) "# lambdaHi: ",lambdaHi
!!$
!!$  end subroutine PrintDist

  !****************************************************************************
  !****s* randomMaxwell/initMaxwell
  ! NAME
  ! subroutine initMaxwell(mass,temp)
  ! PURPOSE
  ! initialize the random number generator with given mass and temperature
  !****************************************************************************
  subroutine initMaxwell(mass,temp)

    real, intent(in) :: mass, temp

    mass_ = mass
    A = mass/temp
    pM2 = 2*(1+sqrt(1+A**2))/A**2
    pM = sqrt(pM2)
    logfpM = log(pM2)-A*pM2/(1+sqrt(1+pM2))
    fpM = exp(logfpM)

    pLo = findLo(pM)
    pHi = findHi(pM)

    ! lambda(p) = f(p)/f'(p):
    lambdaLo =  1/ ( 2/pLo - A * pLo/sqrt(pLo**2+1))
    lambdaHi = -1/ ( 2/pHi - A * pHi/sqrt(pHi**2+1))

    qHi = lambdaHi/(pHi-pLo)
    qLo = lambdaLo/(pHi-pLo)
    qM = 1 - (qHi+qLo)

!    call printDist

  end subroutine initMaxwell

  !****************************************************************************
  !****s* randomMaxwell/rnMaxwell
  ! NAME
  ! real function rnMaxwell()
  ! PURPOSE
  ! return a random number according the relativistic Maxwell distribution,
  ! which has been initialized with a mass and a temperature
  !****************************************************************************
  real function rnMaxwell()

    use random, only: rn

    real :: U,V,X,Y,EE

    do
       U = rn()
       V = rn()

       if (U <= qM) then
          Y = U/qM
          X = (1-Y)*(pLo+lambdaLo) + Y*(pHi-lambdaHi)
          if (V*fpM <= funk(X)) exit ! ==> done
       else if (U <= qM+qHi) then
          EE = -log((U-qM)/qHi)
          X = pHi - lambdaHi*(1-EE)
          if (V*fpM*exp(-EE) <= funk(X) ) exit ! ==> done
       else
          EE = -log((U-(qM+qHi))/qLo)
          X = pLo + lambdaLo*(1-EE)
          if (V*fpM*exp(-EE) <= funk(X) ) exit ! ==> done
       end if

    end do

    rnMaxwell = X * mass_

  end function rnMaxwell

end module randomMaxwell
