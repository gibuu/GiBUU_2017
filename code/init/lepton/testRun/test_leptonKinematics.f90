program test
use leptonKinematics
implicit none

real :: e_i, e_f, theta, QS

e_i=0.7
e_f=0.57
theta=32.

Qs= evaluate_QSquared(theta,e_i,e_f)
write(*,*) e_i,e_f, theta, Qs
theta= evaluate_theta_lf(Qs,e_i,e_f)
write(*,*) e_i,e_f, theta, Qs

end program test
