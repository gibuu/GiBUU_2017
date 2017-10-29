program tester


call test


contains


  subroutine test
    use rotation
    use constants,only :pi
    implicit none
    real, dimension(1:3) :: x,y
    real :: theta, phi , th,ph

    x=(/ 1.,0.,0. /)

    phi=0
    theta=pi/4.
    call rotate(phi,theta,x,y)
    call wr2(x,y,phi,theta)

    phi=0.
    theta=-pi/4.
    call rotate(phi,theta,y,x)
    call wr2(y,x,phi,theta)


    phi=pi
    theta=pi/4.
    call rotate(phi,theta,x,y)
    call wr2(x,y,phi,theta)


    phi=2*pi
    theta=pi/4.
    call rotate(phi,theta,x,y)
    call wr2(x,y,phi,theta)


    phi=pi
    theta=pi/2.
    call rotate(phi,theta,x,y)
    call wr2(x,y,phi,theta)




    x=(/ 1.,0.,0. /)
    call get_phi_theta(x,th,ph)
    call wr(x,th,ph)


    x=(/ -1.,0.,0. /)
    call get_phi_theta(x,th,ph)
    call wr(x,th,ph)

    x=(/ 0.,1.,0. /)
    call get_phi_theta(x,th,ph)
    call wr(x,th,ph)


    x=(/ 0.,-1.,0. /)
    call get_phi_theta(x,th,ph)
    call wr(x,th,ph)


    x=(/ 0.,0.,-1. /)
    call get_phi_theta(x,th,ph)
    call wr(x,th,ph)

    x=(/ 0.,0.,1. /)
    call get_phi_theta(x,th,ph)
    call wr(x,th,ph)


    x=(/ 0.,-1.,1. /)
    call get_phi_theta(x,th,ph)
    call wr(x,th,ph)

    x=(/ 0.,-1.,-1. /)
    call get_phi_theta(x,th,ph)
    call wr(x,th,ph)



    x=(/ 0.,1.,1. /)
    call get_phi_theta(x,th,ph)
    call wr(x,th,ph)

    x=(/ 1.,1.,1. /)
    call get_phi_theta(x,th,ph)
    call wr(x,th,ph)


  end subroutine test

      subroutine wr(x,th,ph)
        use constants, only : pi
        implicit none
        real :: x(1:3),th, ph
        write(*,'(A,3F6.2)') 'x=',x
        write(*,'(A,F6.2,A,F6.2,A)') 'theta=',th/pi,'*pi.  phi=',ph/pi,'*pi'
        write(*,*)
      end subroutine wr
      subroutine wr2(x,y,th,ph)
        use constants, only : pi
        implicit none
        real :: x(1:3),y(1:3),th, ph
        write(*,'(A,3F6.2,A,3F6.2)') 'x=',x, '->',y
        write(*,'(A,F6.2,A,F6.2,A)') 'theta=',th/pi,'*pi.  phi=',ph/pi,'*pi'
        write(*,*)
      end subroutine wr2


end program tester
