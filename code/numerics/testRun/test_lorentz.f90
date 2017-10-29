program test

  call test_lorentz
  call test_boostToCM
end program test


subroutine test_boostToCM
  use random
  use vector
  use lorentzTrafo
  implicit none
  real, dimension(0:3) :: a,b
  integer ::i,j
  real :: mass
  do i=1,100000
     do j=1,3
        a(j)=rn()*10.
        b(j)=rn()*10.
     end do
     mass=rn()*10.
     a(0)=sqrt(mass**2+Sum(a(1:3)**2))
     mass=rn()*10.
     b(0)=sqrt(mass**2+Sum(b(1:3)**2))
     call boostToCM(a,b)
     write(100,'(8E12.3)') a,b
     if(absVec(a(1:3)+b(1:3)).gt.1E-5) then
        write(*,*) 'ERROR! NOT CM!'
        write(*,*) a
        write(*,*) b
     end if
  end do
 

end subroutine test_boostToCM

subroutine test_lorentz
  use minkowski, only : abs4
  use lorentzTrafo
  real, dimension(1:3) :: beta, beta2
  real :: srts
  real, dimension(0:3) :: mom, mom2



  mom=(/1.35369180089968,      -0.164540101269438,   4.761940348776917E-002,  0.381479778916432 /)
  mom2=(/abs4(mom),0.,0.,0./)
  write(*,*) mom
  beta=mom(1:3)/mom(0)

  call lorentz(beta, mom)
  write(*,*) mom
  call lorentz(-beta, mom2)
  write(*,*) mom2
end subroutine test_lorentz



