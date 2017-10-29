program test

  use argonnev18
  use inputGeneral, only: readInputGeneral

  implicit none

  call readInputGeneral

  call argonne_test
  call testDeuteriumPot
  
contains

  !*************************************************************************
  ! NAME
  ! subroutine argonne_test()
  ! PURPOSE
  ! * This function tests the deuteron wave functions. 
  ! * Output: Normalization, sqrt(<k^2>, sqrt(<r^2>))=2* RMS radius
  !*************************************************************************
  subroutine argonne_test()

    real :: r,sum,rms,wf,k
    real,parameter :: dr=0.001
    real,parameter :: dk=0.001

    write(*,*) '*****************************************'
    write(*,*) 'Testing position-space wave function:'
    write(*,*) '*****************************************'
    sum=0.
    rms=0.
    r=0.
    do while (r<=20.)
       wf=argonne_WF_rSpace(r)
       write(11,*) r, wf
       rms=rms+r**2*wf*dr
       sum=sum+wf*dr
       r=r+dr
    end do
    write(*,*)'sqrt(<r^2>)=',sqrt(rms),'fm => Radius=',sqrt(rms)/2.,'fm'
    write(*,*)'Normalization=',sum
    write(*,*)

    write(*,*) '*****************************************'
    write(*,*) 'Testing momentum-space wave function:'
    write(*,*) '*****************************************'
    sum=0.
    rms=0.
    k=0.
    do while(k<=8.)
       wf=argonne_WF_kSpace(k)
       write(12,*) k, wf
       rms=rms+k**4*wf*dk
       sum=sum+k**2*wf*dk
       k=k+dk
    end do

    write(*,*)'sqrt(<k^2>)=',sqrt(rms),'GeV'
    write(*,*)'Normalization=',sum
    write(*,*)

  end subroutine argonne_test


  subroutine testDeuteriumPot
    integer :: i
    real :: r
    open(10,file="argonnePot.dat")
    do i=0,400
       r = 0.05*i
       write(*,'(5G20.3)') r,  argonne_deuteriumPot(r)
       write(10,'(5G20.3)') r,  argonne_deuteriumPot(r)
    end do
    close(10)
  end subroutine testDeuteriumPot


end program test
