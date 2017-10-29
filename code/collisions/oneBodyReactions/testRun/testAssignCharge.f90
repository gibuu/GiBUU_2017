program testAssignCharge
  use master_1body, only: assignCharge
  use IdTable
  use particleProperties, only: initParticleProperties
  use particleDefinition

  implicit none

  integer :: j, charge
  type(particle),dimension(1:2) :: pair


  call initParticleProperties

  Write(*,*) 'Testing AssignCharge'

  pair%ID=(/119,101/)


  call assignCharge(pair,121,-1)

  write(*,*) pair%charge

  pair%ID=(/101,119/)
  call assignCharge(pair,121,-1)
  write(*,*) pair%charge


  pair%ID=(/ds_plus,101/)


  call assignCharge(pair,dsStar_plus,1)

  write(*,*) pair%charge
  pair%ID=(/101,ds_plus/)


  call assignCharge(pair,dsStar_plus,1)
  write(*,*) pair%charge

  do charge=-1,1
     do j=1,10
        pair%ID=(/pion,pion/)
        call assignCharge(pair,rho,charge)
        write(*,'(A,I2,A,I2,A,I2,A)') 'rho(',charge,')=> pion(',pair(1)%charge,') pion(',pair(2)%charge,')'
     end do
     write(*,*) '*********************************'
  end do


  write(*,*) '*********************************'
  write(*,*) '*********************************'
  write(*,*) '*********************************'
  write(*,*) '*********************************'


  pair%ID=(/pion,dMeson/)
  call assignCharge(pair,dStar,0)
  write(*,*) pair%charge

  pair%ID=(/pion,dBar/)
  call assignCharge(pair,dStarBar,0)
  write(*,*) pair%charge
  write(*,*) '*********************************'


  do j=1,10
     pair%ID=(/pion,dMeson/)
     call assignCharge(pair,dStar,1)
     write(*,*) pair%charge
  end do
  write(*,*) '*********************************'

  do j=1,10

     pair%ID=(/pion,dBar/)
     call assignCharge(pair,dStarBar,-1)
     write(*,*) pair%charge
  end do





end program testAssignCharge
