program PhaseSpace2p2h

  use random, only : rnOmega
  use twoBodyTools, only : pCM_sqr
  use minkowski, only : abs3,abs4
  use lorentzTrafo, only : lorentz
  use inputGeneral
  use particleProperties, only: initParticleProperties

  implicit none


  real :: nu, Q2, Ebeam
  real, parameter :: pF = 0.250, mN=0.938
  real :: srts, pcm, pp
  real, dimension(0:3) :: q,p1,p2, pp1,pp2, sum, sum1, k,kprime
  logical :: blocked

  integer :: i1,i2
  integer, parameter :: n1=1000, n2=10
  integer :: sum0,sum00,sum01,sum02,sumblocked
  integer :: iNu
  integer, parameter :: nNu = 70
  real,parameter :: nuDelta=0.01 
  real :: I0, Iblocked
  real :: M2,M2sum,sumPS,sumPSDelta
  logical :: flagOK 

  call readInputGeneral
  call initParticleProperties


  Ebeam = 1.0
!  Ebeam = 0.262

  Q2 = 0.205
!  Q2 = 0.1

  ! write some kinematical constraints:

  nu = Q2/(2*mN)
  write(31,'(2f13.5,/,2f13.5,/)') nu,0.0, nu,1.0 ! vacuum peak position
  
  nu = Q2/(2*mN)*(sqrt(1+(pF/mN)**2)-(pF/mN)*sqrt(1+(4*mN**2/Q2)))
  write(31,'(2f13.5,/,2f13.5,/)') nu,0.0, nu,1.0 ! minimal nu for QE

  nu = Q2/(2*mN)*(sqrt(1+(pF/mN)**2)+(pF/mN)*sqrt(1+(4*mN**2/Q2)))
  write(31,'(2f13.5,/,2f13.5,/)') nu,0.0, nu,1.0 ! maximal nu for QE

  nu = (1.232**2-mN**2+Q2)/(2*mN)
  write(31,'(2f13.5,/,2f13.5,/)') nu,0.0, nu,1.0 ! vacuum Delta peak

  nu = Ebeam-Q2/(4*Ebeam)
  write(31,'(2f13.5,/,2f13.5,/)') nu,0.0, nu,1.0 ! maximal possible nu


  call DoQE()
  call Do2p2h()

contains
  subroutine DoQE()
    use quasiElastic_electron, only : matrixElement


    real :: cosalpha

    do iNu = 1,nNu*10
       nu = iNu*nuDelta/10
       if (nu.gt.Ebeam) cycle

       q = (/nu, 0, 0, sqrt(nu**2+Q2) /)
       call SetK(Ebeam,nu,Q2,flagOK)
       if (.not. flagOK) cycle


       I0 = 0.
       Iblocked = 0.
       M2sum = 0.0

       do i1=0,n1
          pp = pF*real(i1)/n1
          cosalpha = (2*sqrt(mN**2+pp**2)*q(0) - Q2)/(2*pp*q(3))

          if (cosalpha < -1) cycle
          if (cosalpha >  1) cycle
          
          I0 = I0 + pp**2*pF/n1

          p1(1:3) = (/0.0, pp*sqrt(1-cosalpha**2), pp*cosalpha/)
          p1(0) = sqrt(mN**2+pF**2)

          pp1 = p1 + q

          M2 = matrixElement(p1,pp1,k,kprime,1)
             
!!$             write(*,'(A,": ",1P,4e13.5)') "q    ",q
!!$             write(*,'(A,": ",1P,4e13.5)') "k    ",k
!!$             write(*,'(A,": ",1P,4e13.5)') "k'   ",kprime
!!$             write(*,'(A,": ",1P,4e13.5)') "p1   ",p1
!!$             write(*,'(A,": ",1P,4e13.5)') "pp1  ",pp1
!!$             write(*,*) M2
             
          M2Sum = M2Sum + M2*pp**2*pF/n1 

          blocked = .false.
          if (abs3(pp1).lt.pF) blocked = .true.

          if (.not.blocked) then

             Iblocked = Iblocked + pp**2*pF/n1 
             
          end if

       end do
       if (i0>0) then
          write(*,'(3f13.5)') nu,Iblocked/I0
          write(21,'(3f13.5)') nu,Iblocked/I0
          write(22,'(3f13.5)') nu,M2Sum/I0
       end if
    end do

  end subroutine DoQE

  subroutine Do2p2h

    do iNu = 1,nNu*10
       nu = iNu*nuDelta/10
       if (nu.gt.Ebeam) cycle

       q = (/nu, 0, 0, sqrt(nu**2+Q2) /)
       call SetK(Ebeam,nu,Q2,flagOK)
       if (.not. flagOK) cycle

       sum0 = 0
       sum00 = 0
       sum01 = 0
       sum02 = 0
       sumblocked = 0
       sumPS = 0
       sumPSDelta = 0

       do i1=1,n1

          p1(1:3) = pF*rnOmega()
          p1(0) = sqrt(mN**2+pF**2)

          p2(1:3) = pF*rnOmega()
          p2(0) = sqrt(mN**2+pF**2)

          sum = q + p1 + p2
          srts = abs4(sum)

          sum00 = sum00 + 1

          if (srts.gt.mN+1.232) sum02=sum02+1

          pcm = pcm_sqr(srts**2,mN**2,1.232**2)
          if (pcm > 0) then
             sumPSDelta = sumPSDelta + 3.14*sqrt(pcm)/srts
          end if

          pcm = pcm_sqr(srts**2,mN**2,mN**2)
          if (pcm < 0) then
!             write(*,*) 'pcm**2<0: ',pcm
             cycle
          end if
          pcm = sqrt(pcm)

          sum01 = sum01 + 1
          sumPS = sumPS + 3.14*pcm/srts



          do i2=1,n2

             sum0 = sum0+1

             pp1(1:3) = pcm*rnOmega()
             pp1(0) = sqrt(mN**2 + pcm**2)
             pp2(1:3) = -pp1(1:3)
             pp2(0) = pp1(0)

             call lorentz(-sum(1:3)/sum(0), pp1)
             call lorentz(-sum(1:3)/sum(0), pp2)

             sum1 = pp1+pp2

!!$        write(*,'(A,": ",1P,4e13.5)') "q    ",q
!!$        write(*,'(A,": ",1P,4e13.5)') "p1   ",p1
!!$        write(*,'(A,": ",1P,4e13.5)') "p2   ",p2
!!$        write(*,'(A,": ",1P,4e13.5)') "sum  ",sum
!!$        write(*,'(A,": ",1P,4e13.5)') "pp1  ",pp1
!!$        write(*,'(A,": ",1P,4e13.5)') "pp2  ",pp2
!!$        write(*,'(A,": ",1P,4e13.5)') "sum1 ",sum1

             blocked = .false.
             if (abs3(pp1).lt.pF) blocked = .true.
             if (abs3(pp2).lt.pF) blocked = .true.
!!$           write(*,*) 'blocked:',blocked
             if (.not.blocked) sumblocked=sumblocked+1
          end do
       end do

!!$  write(*,*) 'sum:    ',sum0
!!$  write(*,*) 'blocked:',sumblocked

       write(*,'(2f13.5)') nu,real(sumblocked)/sum0
       write(11,'(2f13.5)') nu,real(sumblocked)/sum0
       write(12,'(5f13.5)') nu,sumPS/sum01,sumPS/sum00,sumPSDelta/sum00
       write(13,'(5f13.5)') nu,real(sum01)/sum00,real(sum02)/sum00

    end do

  end subroutine Do2p2h

  subroutine SetK(Ebeam,nu,Q2,flagOK)
    real, intent(in) :: Ebeam,nu,Q2
    logical, intent(out) :: flagOK
    real :: pZ,pX

    flagOK = .false.
    k = 0.0
    pZ = (2*Ebeam*nu + Q2) / (2*sqrt(nu**2+Q2))
    pX = ( Ebeam**2-pZ**2 )
    if (pX<0) then
       write(*,*) 'pX-failure'
       return ! --- failure
    end if
    pX = sqrt( pX )

    k = (/Ebeam,    pX, 0E0, pZ/)
    kprime = k-q
    flagOK = .true.
  end subroutine SetK

end program PhaseSpace2p2h
