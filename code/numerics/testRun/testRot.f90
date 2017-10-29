program testRot
! Test rotations for electron init
use degRad_conversion
use vector
!use lowElectron, only : rotateEvent
implicit none

real, parameter :: eIn =2.
real, parameter :: eOut=1.
real, parameter :: theta=30.
real :: phi=129.
real :: theta_q,phi_q
real, dimension(0:3) :: In, Out, q
character(20) :: Format


write(*,*) 'Please insert phi-angle of the scattered electron'
read *, phi
write(*,*) 'Phi=',phi

Format='(A10,4F15.5,A,F15.5)'

In=(/eIn,0.,0.,eIn/)
Out(0)=eOut
Out(1:3)=sphericalVector(theta,phi,eOut)
q=In-Out


write(*,*) AbsVec(In(1:3))

write(*,format) 'In',In    , ' AbsVec=', AbsVec(In(1:3))
write(*,format) 'Out',Out  , ' AbsVec=', AbsVec(out(1:3))
write(*,format) 'q',q      , ' AbsVec=', AbsVec(q(1:3))

call rotateEvent(In,Out,q,theta_q,phi_q)

write(*,format) 'In',In    , ' AbsVec=', AbsVec(In(1:3))
write(*,format) 'Out',Out  , ' AbsVec=', AbsVec(out(1:3))
write(*,format) 'q',q      , ' AbsVec=', AbsVec(q(1:3))

call rotateEventBack(In,Out,q,theta_q,phi_q)

write(*,format) 'In',In    , ' AbsVec=', AbsVec(In(1:3))
write(*,format) 'Out',Out  , ' AbsVec=', AbsVec(out(1:3))
write(*,format) 'q',q      , ' AbsVec=', AbsVec(q(1:3))

contains


  subroutine rotateEvent(a,b,q,theta_q,phi_q)
    use rotation, only : get_phi_Theta,rotateZY
    real, intent(out) :: theta_q,phi_q
    real, intent(inout),dimension(0:3) :: a,b,q
    ! Rotate coordinate system such that q lies along z-axis, like in figure 6 of Drechsel
    ! and Tiator ( J. Phys. G: Nucl. Part. Phys 18(1992),449-497)
    call get_phi_Theta(q(1:3),theta_q,phi_q)
    write(*,*) theta_q,phi_q
    write(*,format) 'Check get_phi_theta:' ,sphericalVector(degrees(theta_q),degrees(phi_q),AbsVec(q(1:3)))-q(1:3)
    q(1:3) = rotateZY(theta_q,phi_q,q(1:3))
    a(1:3) = rotateZY(theta_q,phi_q,a(1:3))
    b(1:3) = rotateZY(theta_q,phi_q,b(1:3))
  end subroutine rotateEvent

  subroutine rotateEventBACK(a,b,q,theta_q,phi_q)
    use rotation, only : get_phi_Theta,rotateYZ
    real, intent(in) :: theta_q,phi_q
    real, intent(inout),dimension(0:3) :: a,b,q
    ! Rotate coordinate system such that q lies along z-axis, like in figure 6 of Drechsel
    ! and Tiator ( J. Phys. G: Nucl. Part. Phys 18(1992),449-497)
    q(1:3) = rotateYZ(theta_q,phi_q,q(1:3))
    a(1:3) = rotateYZ(theta_q,phi_q,a(1:3))
    b(1:3) = rotateYZ(theta_q,phi_q,b(1:3))
  end subroutine rotateEventBACK


end program testRot
