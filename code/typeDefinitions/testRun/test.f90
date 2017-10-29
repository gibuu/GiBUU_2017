program test
  use particleDefinition
  
  implicit none
  type(particle) :: teilchen
  type(particle),dimension(3:7) :: teilchenVector
  type(particle),dimension(1:2,8:10)  :: teilchenMatrix

  write(*,*) "**********************************" 
  write(*,*) "**********************************" 
  write(*,*)
  Write(*,*)      "Testing particleDefinition "
  write(*,*)
  write(*,*) "**********************************" 
  write(*,*) "**********************************" 

  call setNumber(teilchen)
  write(*,*) teilchen%number

  call setNumber(teilchenVector)
  write(*,*) teilchenVector%number

  call setNumber(teilchenMatrix)
  write(*,*) teilchenMatrix%number

end program test
