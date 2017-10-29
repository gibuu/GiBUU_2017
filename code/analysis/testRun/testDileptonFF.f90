program testDileptonFF

  use Dilepton_Analysis, only: FF_omega_effe, FF_omega_terschluesen, FF_VMD, &
                               FF_eta, FF_etaprime_genesis, FF_etaprime_terschluesen, &
                               Delta_FF_Dipole, Delta_FF_VMD, Delta_FF_MAID, Delta_FF_Iachello
  use inputGeneral, only: readInputGeneral
  use particleProperties, only: initParticleProperties

  implicit none

  real, parameter :: dm = 0.001

  integer :: i
  real :: m_ee, msqr

  call readInputGeneral
  call initParticleProperties

!   call test_omega
  call test_etaprime
!   call test_Delta

contains

  subroutine test_omega
    m_ee = dm
    open (911,file="omegaDalitzFF.dat")
    do i=1,1000
      write (911,'(4ES13.5)') m_ee, FF_omega_effe(m_ee), FF_VMD(m_ee), FF_omega_terschluesen(0.782, m_ee)
      m_ee = m_ee + dm
    end do
    close (911)
  end subroutine

  subroutine test_etaprime
    m_ee = dm
    open (912,file="etaprimeDalitzFF.dat")
    do i=1,1200
      write (912,'(5ES13.5)') m_ee, FF_eta(m_ee), FF_VMD(m_ee), FF_etaprime_genesis(m_ee), FF_etaprime_terschluesen(m_ee)
      m_ee = m_ee + dm
    end do
    close (912)
  end subroutine

  subroutine test_Delta
    use NDeltaFF_Ramalho, only: NDeltaTL
    use distributions, only: markusPostFormfactor
    use constants, only: melec, mN, mPi
    use particleProperties, only: hadron
    use IdTable, only: Delta

    real :: FF(1:5),W
    open (913,file="DeltaDalitzFF.dat")
    m_ee = 0.
    msqr = 0.
    do i=1,200
      FF(1) = NDeltaTL (msqr, 1.23)
      FF(2) = NDeltaTL (msqr, 1.43)
      FF(3) = NDeltaTL (msqr, 1.63)
      FF(4) = NDeltaTL (msqr, 1.83)
      FF(5) = NDeltaTL (msqr, 2.03)
      write (913,'(13G12.5)') m_ee, Delta_FF_VMD(msqr),Delta_FF_Dipole(msqr), &
                              Delta_FF_Iachello(1.23,msqr),Delta_FF_Iachello(1.43,msqr),Delta_FF_Iachello(1.63,msqr), &
                              Delta_FF_Iachello(1.83,msqr),Delta_FF_Iachello(2.03,msqr), FF(1:5)
      m_ee = m_ee + 0.01
      msqr = m_ee**2
    end do
    close(913)

    open(914,file="DeltaDalitzFF_Wdep.dat")
    msqr=2*melec
    do i=1,100
      W = mN+mPi+i*0.01
      FF(1) = NDeltaTL (msqr, W)
      write (914,'(4G12.5)') W, Delta_FF_Iachello(W,msqr), FF(1), markusPostFormfactor (W, hadron(Delta)%mass, mN+mPi, 1.071)
    end do
    close(914)
  end subroutine

end
