program testMink

  use constants, only: ii
  use minkowski
  use matrix_Module, only: unit2, unit4, trace, MatrixMult, dagger

  implicit none

  integer :: i,j, mu, nu
  complex,dimension(0:3,0:3) :: gam
  real,dimension(0:3) :: q,p
  
  complex, dimension(1:2,1:2,1:3), parameter :: sigma = reshape((/sigma1,sigma2,sigma3/),(/2,2,3/))

!   call printMatrices()
!   call testSlash()
!   call testDagger()
!   call testSigma()
!   call testMisc()
  call testLevi()


contains


  subroutine testSlash
    ! Test Slashing:
    write(*,*)
    write(*,*) 'Test 1:'
    q=(/1.23,1.9,2.3,7.1/)
    write(*,'(A,4F15.4)') 'q=',q
    write(*,*) 'tr [slashed(q) slashed(q)] - 4 q^mu q^_mu=0'
    write(*,'(A, 4F15.4)')' = ', trace(MatMul(slashed(q),slashed(q)))-4*SP(q,q)

    write(*,*)
    write(*,*) 'Test 2:'
    p=(/21.23,17.9,7.3,7.121/)
    write(*,'(A, 4F15.4)') 'p=',q
    write(*,'(A, 4F15.4)') 'q=',q
    write(*,*) 'tr [slashed(q) slashed(p)] - 4 q^mu p^_mu=0'
    write(*,'(A,4F15.4)')' = ', trace(MatMul(slashed(q),slashed(p)))-4*SP(q,p)
  end subroutine


  subroutine testDagger
    ! Test daggering
    write(*,*)
    write(*,*) 'Test 3: gamma(i)^dagger-gamma_0 gamma(i) gamma_0 = 0'
    do i=0,3
      gam=dagger(gamma(:,:,i))-MatrixMult(gamma0,gamma(:,:,i),gamma0)
      write(*,*) 'i=', i
      do j=0,3
        write(*,'(4("(",2F5.1,")"))') gam(j,:)
      end do
      write(*,*)
    end do
  end subroutine


  subroutine testSigma
    ! Test sigma^mu nu
    write(*,*)
    write(*,*) 'Test 3: sigma^mu nu-i/2*(2 gamma^mu gamma^nu-2g^mu nu) = 0'
    do mu=0,3
    do nu=0,3
      write(*,*) '#mu,nu=',mu,nu
      gam=sigma4(mu,nu)-ii*(MatrixMult(gamma(:,:,mu),gamma(:,:,nu))-metricTensor(mu,nu)*unit4)
      call printMatrixF(gam)
      write(*,*)
    end do
    end do
  end subroutine

  
  subroutine testMisc

    write(*,*) 's_z(initial)=-1/2'
    write(*,*)
    write(*,*) '(sigma4(1,0)+ii*sigma4(2,0))*(/0.,1.,0.,0. /)'
    write(*,*) MatMul(MatMul(sigma4(1,0)+ii*sigma4(2,0),gamma5),(/0.,1.,0.,0. /))
    write(*,*)
    write(*,*) '(sigma4(1,3)+ii*sigma4(2,3))*(/0.,1.,0.,0. /)'
    write(*,*) MatMul(MatMul(sigma4(1,3)+ii*sigma4(2,3),gamma5),(/0.,1.,0.,0. /))
    write(*,*)
    write(*,*) '(gamma(1)+ii*gamma(2))*(/0.,1.,0.,0. /)'
    write(*,*) MatMul(MatMul(gamma1+ii*gamma2,gamma5),(/0.,1.,0.,0. /))

    write(*,*)
    write(*,*) 's_z(initial)=1/2'
    write(*,*)
    write(*,*) 'MatMul(MatMul(sigma4(0,3),gamma(5)),(/1.,0.,0.,0. /))'
    write(*,*) MatMul(MatMul(sigma4(0,3),gamma5),(/1.,0.,0.,0. /))
    write(*,*)
    write(*,*) 'MatMul(MatMul(gamma(0),gamma(5)),(/1.,0.,0.,0. /))'
    write(*,*) MatMul(MatMul(gamma0,gamma5),(/1.,0.,0.,0. /))
    write(*,*)
    write(*,*) 'MatMul(MatMul(gamma(3),gamma(5)),(/1.,0.,0.,0. /))'
    write(*,*) MatMul(MatMul(gamma3,gamma5),(/1.,0.,0.,0. /))
    write(*,*)
    write(*,*)
    write(*,*) 'gamma(1)+i gamma(2)'
!     call printMatrixF(gamma(1)+ii*gamma(2))

 !    write(*,*) '10+i 20'
  !   call printMatrixF(sigma4(1,0)+ii*sigma4(2,0))
  !   write(*,*) '13+i 23'
   !  call printMatrixF(sigma4(1,3)+ii*sigma4(2,3))

  end subroutine


  subroutine testLevi
    use output, only: timeMeasurement
    integer :: i,l,a(1:4)

    ! (1) test a few arbitrary input combinations of the Levi-Civita tensor
    ! should be 1:
    print *,levi_civita(0,1,2,3), levi_civita(3,2,1,0), levi_civita(1,2,3,4), levi_civita(4,3,2,1)
    ! should be -1:
    print *,levi_civita(0,3,2,1),levi_civita(1,0,2,3),levi_civita(0,2,1,3),levi_civita(1,4,3,2)
    print *,levi_civita(2,1,3,4),levi_civita(1,3,2,4)
    !should be 0:
    print *,levi_civita(0,0,0,0),levi_civita(3,2,4,3)
    
    ! (2) do a small benchmark with random numbers
    call timeMeasurement(.true.)
    do i=1,1000000000
      a = mod ((/ irand(), irand(), irand(), irand() /), 4)
      l = levi_civita(a(1), a(2), a(3), a(4))
!       print *,a(1:4),l
    end do
    call timeMeasurement()

  end subroutine


  subroutine printMatrixF(A)
    integer :: j
    complex, intent(in),  dimension(0:3,0:3) :: a 
    
    do j=0,3
       write(*,'(4("(",2F5.1,")"))') A(j,:)
    end do
  end subroutine printMatrixF


  subroutine printMatrices()
    integer :: i,j,k
    complex,dimension(1:2,1:2) :: sig
    complex,dimension(0:3,0:3) :: gam
    !character, dimension(20) :: formatGamma='(4("(",2F5.1,")"))'

    write(*,*) '*********************'
    write(*,*) 'Sigmas='
    write(*,*) '*********************'
    do i=1,3
       write(*,'(A,I4,A)') '# Sigma(',i,')'
       sig=sigma(:,:,i)
       write(*,'(2("(",2F5.1,")"))') sig(1,:)
       write(*,'(2("(",2F5.1,")"))') sig(2,:)
       write(*,*)
    end do
    write(*,*)
    write(*,*)
    write(*,*) '*********************'
    write(*,*) 'Gammas='
    write(*,*) '*********************'
    do i=0,3
       write(*,'(A,I4,A)') '# Gamma(',i,')'
       do j=0,3
          write(*,'(4("(",2F5.1,")"))') gamma(j,:,i)
       end do
       write(*,*)
    end do

    write(*,'(A)') '# Gamma(5)'
    gam=ii*MatMul(MatMul(MatMul(gamma0,gamma1),gamma2),gamma3)
    do j=0,3
       write(*,'(4("(",2F5.1,")"))') gam(j,:)
    end do
    write(*,*)

    do j=0,3
       write(*,'(4("(",2F5.1,")"))') gamma5(j,:)
    end do
    write(*,*)

    do i=6,11
       write(*,'(A,I4,A)') '# Gamma(',i,')'
       do j=0,3
          write(*,'(4("(",2F5.1,")"))') gamma(j,:,i)
       end do
       write(*,*)
    end do

    write(*,*) '*********************'
    write(*,*) 'Unit Matrices='
    write(*,*) '*********************'

    write(*,'(A)') '# 4-dim'
    do j=0,3
       write(*,'(4("(",2F5.1,")"))') unit4(j,:)
    end do
    write(*,*)

    write(*,'(A)') '# 2-dim'
    write(*,'(2("(",2F5.1,")"))') unit2(1,:)
    write(*,'(2("(",2F5.1,")"))') unit2(2,:)

    write(*,*) '*********************'
    write(*,*) 'Sigma4='
    write(*,*) '*********************'

    do j=0,3
       do k=0,3
          write(*,'(A,I4,A,I4,A)') '# Sigma4(',j,',',k,')'
          gam = sigma4par(:,:,j,k)
          do i=0,3
             write(*,'(4("(",2F5.1,")"))') gam(i,:)
          end do
          write(*,*)
       end do
    end do

  end subroutine printMatrices


end program testMink
