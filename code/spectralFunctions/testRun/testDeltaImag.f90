program test
  use inputGeneral



  call readinputGeneral
  call init_Database

  call tester

end program test



subroutine tester
  ! This routine produces output to compare our imaginary part to Oset et al., especiallly fig. 12 of NPA468(1987), 631-652
  use selfenergy_imagPart
  use idTable, only : delta, pion,nucleon
  use particleProperties
  use  baryonWidthMedium_tables
  use deltaWidth
  use constants, only : pi, hbarc,rhoNull
  use inMediumWidth, only : evaluateCollisionBroadening

  implicit none
  real, parameter :: rhoN=0.75*rhoNull/2.
  real, parameter :: rhoP=0.75*rhoNull/2.

  real :: p_pi, T_pi, E_pi, m_delta ! pion momentum and kinetic energy and energy
  real :: k_f ! Fermi momentum
  real :: p_delta, E_delta ! delta momentum and energy
  integer, parameter :: particleID=delta
  integer :: i
  real :: osetGamma,imsig2,imsig3,imsigq,buuwidth,buu_elastic

  ! Move along piN scattering trajectory in E-p plane (see Osets paper for details)
  open(10,file='delta_selfenergy_imag.oset.dat')
  write(10,'(A)') '#T_pi, p_Delta, E_delta, m_delta,osetGamma,2.*imsig2, 2*imsig3,2*imsigq,buu_full,buu_elastic'

  do i=0,40
     p_pi=float(i)*0.04

     E_pi=sqrt(p_pi**2+meson(pion)%mass**2)
     T_pi=E_pi-meson(pion)%mass
     k_f=(1.5*pi**2*(rhon+rhop))**(1./3.)*hbarc

     P_delta=sqrt(3./5.*k_f**2+p_pi**2)
     E_delta=3./5.*k_f**2/2./baryon(nucleon)%mass+E_pi+baryon(nucleon)%mass
     m_delta=sqrt(-p_delta**2+E_delta**2)
     write(*,*) 'm,p=', m_delta,p_delta
     call deloset(m_delta,p_delta,rhoN+rhoP,imsig2,imsig3,imsigq)
     buuWidth=evaluateCollisionBroadening(particleID,p_delta,m_delta,rhoN*0.197**3,rhoP*0.197**3,1000,buu_elastic)
     osetGamma=2.*(imsig2+imsig3+imsigq)
     !     write(10,'(11G10.3)') T_pi, p_Delta, E_delta, m_delta,&
     !          & get_inMediumWidth(particleID,P_delta,m_delta,rhoN,rhoP,1,.false.),osetGamma &
     !  & ,2.*imsig2, 2*imsig3,2*imsigq,buuWidth,buu_elastic
     write(10,'(11G10.3)') T_pi, p_Delta, E_delta, m_delta &
          & ,osetGamma  ,2.*imsig2, 2*imsig3,2*imsigq,buuWidth,buu_elastic
     !selfenergy_imag(particleID,P_delta,E_delta,rhoN,rhoP)
  end do
  close(10)
  write(*,*) '************************************************************'
  write(*,*) 'Finished TEST 1/3'
  write(*,*) '************************************************************'



  ! Now keep mass of incoming delta constant: m=m_pole
  open(10,file='delta_selfenergy_imag.pole.dat')
  write(10,'(A)') '#p_Delta, E_delta, m_delta,osetGamma,2.*imsig2, 2*imsig3,2*imsigq,buu_full,buu_elastic'
  m_delta=baryon(delta)%mass
  do i=0,25
     P_delta=float(i)*0.04
     E_delta=sqrt(m_delta**2+p_delta**2)

     call deloset(m_delta,p_delta,rhoN+rhoP,imsig2,imsig3,imsigq)
     buuWidth=evaluateCollisionBroadening(particleID,p_delta,m_delta,rhoN*0.197**3,rhoP*0.197**3,1000,buu_elastic)
     osetGamma=2.*(imsig2+imsig3+imsigq)
     write(10,'(11G10.3)') p_Delta, E_delta, m_delta &
          & ,osetGamma  ,2.*imsig2, 2*imsig3,2*imsigq,buuWidth,buu_elastic
  end do
  close(10)
  write(*,*) '************************************************************'
  write(*,*) 'Finished TEST 2/3'
  write(*,*) '************************************************************'



  ! Now keep mass of incoming delta constant: m=1.15 GEV
  open(10,file='delta_selfenergy_imag.1_15.dat')
  write(10,'(A)') '#p_Delta, E_delta, m_delta,osetGamma,2.*imsig2, 2*imsig3,2*imsigq,buu_full,buu_elastic'
  do i=0,25
     P_delta=float(i)*0.04
     m_delta=1.15
     E_delta=sqrt(m_delta**2+p_delta**2)

     call deloset(m_delta,p_delta,rhoN+rhoP,imsig2,imsig3,imsigq)
     buuWidth=evaluateCollisionBroadening(particleID,p_delta,m_delta,rhoN*0.197**3,rhoP*0.197**3,1000,buu_elastic)
     osetGamma=2.*(imsig2+imsig3+imsigq)
     write(10,'(11G10.3)') p_Delta, E_delta, m_delta &
          & ,osetGamma  ,2.*imsig2, 2*imsig3,2*imsigq,buuWidth,buu_elastic
     write(*,'(11G10.3)') p_Delta, E_delta, m_delta &
          & ,osetGamma  ,2.*imsig2, 2*imsig3,2*imsigq,buuWidth,buu_elastic
  end do
  close(10)
  write(*,*) '************************************************************'
  write(*,*) 'Finished TEST 3/3'
  write(*,*) '************************************************************'



end subroutine tester
