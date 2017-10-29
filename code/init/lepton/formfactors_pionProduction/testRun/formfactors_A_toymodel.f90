
module formfactors_A_toymodel


  public :: getA

  logical, parameter :: debug=.false.

contains


  function getA(thetaIn,sIn,QSquaredIn)
    implicit none
    logical, save :: initFlag=.true.
    complex, save , dimension(:,:),Allocatable :: A
    integer , parameter :: N=1000
    real, dimension(1:N),save :: Theta ,s, t
    !integer,save :: maxvalue
    complex, dimension(1:6) :: getA
    real, intent(in) :: thetaIn, QSquaredIn,sIn
    integer :: i,k

    if(initFlag) then
       call init_datFile
       initFlag=.false.
    end if

    do i=lBound(theta,dim=1),ubound(theta,dim=1)
       if(thetaIn.gt.theta(i)) then
          cycle
       else
          getA=A(:,i)
          if(debug) then 
             write(*,'(A,F12.5)') 's=',sIn
             write(*,'(A,F12.5)') 'W=',sqrt(sIN)
             write(*,'(A,F12.5)') 'Q^2=',qSquaredIn
             write(*,'(A,F12.5)') 'theta=',theta(i)
             do k=1,6
                write(*,'(A,I6,A,2E19.12)') 'A_{',k,'}(theta)=',A(k,i)
                select case(k)
                   case(1)
                      write(*,'(A,2E19.12,A)') '=',A(k,i)*0.197**2, 'fm^2'
                   case(2,5)
                      write(*,'(A,2E19.12,A)') '=',A(k,i)*0.197**4, 'fm^4'
                   case Default
                      write(*,'(A,2E19.12,A)') '=',A(k,i)*0.197**3, 'fm^3'
                end select
             end do 
          end if
          exit
       end if
    end do


  contains
    subroutine init_datFile
      !    open(100,file=trim(path_to_Input)//'/photo_twoPi/gamp-ppipm.dat',status='old')
      implicit none
      integer, parameter :: Dateiende=-1
      real, dimension(1:N),save ::ReSpalte1, ImSpalte1, ReSpalte2,ImSpalte2
      !logical, save :: initFlag=.true.
      integer :: ios
      integer :: i,k

      open(100,file='maid_A.html',status='old')

      i=1
      do 
         read(100,*,IOSTAT=IOS)  Theta(i),s(i),t(i), ReSpalte1(i), ImSpalte1(i), ReSpalte2(i), ImSpalte2(i)
         if(ios.eq.0) i=i+1
         if (IOS.eq.DATEIENDE) exit
      end do
      i=i-1
      !maxvalue=i-1
      close(100)

      Allocate(A(1:6,1:i/3))
      do k=1,i/3
         A(1,k)=CMPLX(ReSpalte1(k),ImSpalte1(k))
         A(2,k)=CMPLX(ReSpalte2(k),ImSpalte2(k))
         A(3,k)=CMPLX(ReSpalte1(i/3+k),ImSpalte1(i/3+k))
         A(4,k)=CMPLX(ReSpalte2(i/3+k),ImSpalte2(i/3+k))
         A(5,k)=CMPLX(ReSpalte1(i*2/3+k),ImSpalte1(i*2/3+k))
         A(6,k)=CMPLX(ReSpalte2(i*2/3+k),ImSpalte2(i*2/3+k))
      end do
      open(101,file='FF_A1.dat')
      open(102,file='FF_A2.dat')
      open(103,file='FF_A3.dat')
      open(104,file='FF_A4.dat')
      open(105,file='FF_A5.dat')
      open(106,file='FF_A6.dat')
      open(107,file='theta.dat')
      do k=1,6
         write(100+k,'(A20,I3,A)') '# theta , A(',k, ')'
      end do

      do i= lBound(A,dim=2),uBound(A,dim=2)
         do k=1,6
            write(100+k,'(1F8.3,2E15.8)') theta(i), A(k,i)
         end do
      end do
      write(107,*) theta(:)

      close(101)
      close(102)
      close(103)
      close(104)
      close(105)
      close(106)
      close(107)
      ! Convert everything from fm to GeV^-1
      ! 197 MeV fm =1 => fm=1/(0.197 GeV)
      A(1,:)=A(1,:)/(0.197**2)
      A(2,:)=A(2,:)/(0.197**4)
      A(3,:)=A(3,:)/(0.197**3)
      A(4,:)=A(4,:)/(0.197**3)
      A(5,:)=A(5,:)/(0.197**4)
      A(6,:)=A(6,:)/(0.197**3)

    end subroutine init_datFile
  end function getA




end module formfactors_A_toymodel
