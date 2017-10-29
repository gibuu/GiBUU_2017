!******************************************************************************
!****m* /NucExRTF
! NAME
! module NucExRTF
! PURPOSE
! This module includes a routine to calculate properties of ground state
! nuclei (p- & n-radii, density profiles, etc.) within
! Relativistic Thomas-Fermi (RTF) program from Horst Lenske.
! NOTES
! * Valid only in the RMF-mode!
! * Coulomb field is always switched on in RTF!
!   use "coulombFlag=true" in jobCard for consistency
!******************************************************************************
module NucExRTF

  implicit none
  private

  public :: NucExRTF_Main

contains

  !****************************************************************************
  !****s* NucExRTF/NucExRTF_Main
  ! NAME
  ! subroutine NucExRTF_Main(nuc%mass,nuc%charge)
  ! PURPOSE
  ! See notes in the header of this module.
  ! NOTES
  ! RTF-code gives also selfenergies (coupling*field) in units of
  ! [MeV] in file "fields.d".
  !****************************************************************************
  subroutine NucExRTF_Main(A_nuc,Z_nuc,RTF_dens)

    use RMF, only: getRMF_flag, getRMF_parSet, &
                    g_sigma, m_sigma, g_omega, m_omega, &
                    g_rho, m_rho, g_2, g_3
    use spline, only: Bsplint2
    use constants, only: pi, mN

    !-----------------------------------------------------------------------
    ! Input-Output variables
    !-----------------------------------------------------------------------
    integer, intent(in) :: A_nuc, Z_nuc
    real, dimension(0:2000,1:2), intent(out) :: RTF_dens
    !-----------------------------------------------------------------------
    ! Local variables
    !-----------------------------------------------------------------------
    real, parameter :: fourpi=4.*pi
    real, parameter :: rmax = 20.
    real, parameter :: dr = 0.1
    real, parameter :: dx = 0.01
    real, parameter :: epst = 1.e-08
    integer :: j
    real :: gs_rtf, gv_rtf, gr_rtf, ms_rtf, mv_rtf, mr_rtf, g2_rtf, g3_rtf
    real :: xn,dummy,mstar, pf  ! ,Ef
    real, dimension(1:200)  :: xx, y1, y2, y3, y4, y5, y6
    real, dimension(0:2000) :: x_new, sigmaPot, omegaPot, rhoPot, coulombPot
    real :: Ef_p,Ef_n,pf_p,pf_n

    !-----------------------------------------------------------------------
    ! Checks
    !-----------------------------------------------------------------------

    if (.not.getRMF_flag()) then
       write(*,*) 'Module NucExRTF / routine NucExRTF_Main:'
       write(*,*) &
            & 'Wrong parameter set for Relativistic Thomas-Fermi'
       write(*,*) &
            & 'ExRTF works only for relativistic mean-field!'
       write(*,*) 'STOP'
       STOP
    end if

    if (getRMF_parSet() == 3 .or. getRMF_parSet() == 6) then ! ExRTF does not run for the NL-Parameter set
                              ! from Lang et al.!
       write(*,*) 'Module NucExRTF / routine NucExRTF_Main:'
       write(*,*) &
            & 'Do not use Lang et al. parameters for ExRTF!!!'
       write(*,*) 'STOP'
       STOP
    end if
    !-----------------------------------------------------------------------
    ! Set variables for ExRTF code
    ! Meson-Baryon coupling constants: g_i**2/(4\pi) (dimensionless units)
    ! Non-linearities of isoscalar, scalar sigma field:
    ! -g_2/(g_sigma**3*M) (dimensionless units)
    !  g_3/(g_sigma**4) (dimensionless units)
    !-----------------------------------------------------------------------
    gs_rtf = g_sigma**2/fourpi
    gv_rtf = g_omega**2/fourpi
    gr_rtf = g_rho**2/fourpi
    ms_rtf = m_sigma*1000. !in units of MeV for ExRTF
    mv_rtf = m_omega*1000. !in units of MeV for ExRTF
    if (g_rho /= 0.0) then
       mr_rtf = m_rho*1000. !in units of MeV for ExRTF
    else
       mr_rtf = 0.0 !in units of MeV for ExRTF
       gr_rtf = 0.0
    end if
    g2_rtf = -g_2 / (g_sigma**3*mN)
    g3_rtf = g_3 / (g_sigma**4)
    !-----------------------------------------------------------------------
    !Generate input-file "ertf.input" for ExRTF code
    !-----------------------------------------------------------------------
    open(1,file='ertf.input')
    write(1,'(2(f5.1,2x))') float(A_nuc),float(Z_nuc)
    write(1,'(2(f5.1,2x))') rmax,dr
    write(1,'(A)') '''ertf.output'''
    write(1,*) '0,1,0,0,0,0,0,0,0,0'
    write(1,'(e12.4)') epst

    !vector - isovector (omega meson field)
    write(1,'(f7.3,1x,i1,1x,i1,f8.4,1x,i1)') &
         & mv_rtf,1,0,gv_rtf,0

    write(1,'(f7.3,1x,i1,1x,i1,f8.4,1x,i1)') &
!                                !          & 980.000,0,1,2.49,0   !scalar - isovector
         & 980.000,0,1,0.00001,0   !vector - isovector !dummy line for ExRTF

    !scalar - isovector (sigma meson field)
    write(1,'(f7.3,1x,i1,1x,i1,f8.4,1x,i1)') &
         & ms_rtf,0,0,gs_rtf,2

    ! non-linearities (power 3)
    write(1,'(i1,1x,e12.5)') &
         & 3,g2_rtf

    ! non-linearities (power 4)
    write(1,'(i1,1x,e12.5)') &
         & 4,g3_rtf

    !vector - isovector (rho meson field)
    if (g_rho .ne. 0.0) then
       write(1,'(f7.3,1x,i1,1x,i1,f8.4,1x,i1)') &
            & mr_rtf,1,1,gr_rtf,0
    end if

    !Photon field
    write(1,'(f7.2,1x,i1,1x,i2,f7.2,1x,i1)') &
         & 0.0,1,-1,1.0,0
    write(1,'(f7.3,1x,i1,1x,i1,f7.2,1x,i1)') &
         & -1.00,0,0,0.00,0        !Terminates Bosons
    !-----------------------------------------------------------------------
    ! Call ExRTF code.
    ! Relevant output: file "gsd.d", which contains density profiles for
    ! protons, neutrons separately.
    !-----------------------------------------------------------------------
    close(1)
    call ExRTF
    !-----------------------------------------------------------------------
    !Prepare density profiles for BUU initialization
    !x(:) : distance (fm)
    !y1(:): proton density profile
    !y2(:): neutron density profile
    !-----------------------------------------------------------------------
    open(2,file='gsd.d')
    read(2,*)
    read(2,*)
    do j=1,200
       read(2,*) xx(j),y1(j),y2(j),dummy,dummy,dummy,dummy,dummy,dummy
    end do
    close(2)
    do j=0,2000
       xn = float(j)*dx
       x_new(j)  = xn
       RTF_dens(j,1) = Bsplint2(xx,y1,xn)
       RTF_dens(j,2) = Bsplint2(xx,y2,xn)
    end do
    !-----------------------------------------------------------------------
    !Check
    !-----------------------------------------------------------------------
    open(3,file='fields.d')
    read(3,*)
    read(3,*)
    do j=1,200
       read(3,*) dummy,y3(j),dummy,y4(j),y5(j),y6(j),dummy,dummy
    end do
    close(3)
    do j=0,2000
       xn = float(j)*dx
       x_new(j)      = xn
       sigmaPot(j)   = Bsplint2(xx,y3,xn)
       omegaPot(j)   = Bsplint2(xx,y4,xn)
       rhoPot(j)     = Bsplint2(xx,y5,xn)
       coulombPot(j) = Bsplint2(xx,y6,xn)
    end do
    !-----------------------------------------------------------------------
    ! Output1: RTF densities
    ! Output2: fields & Fermi-Energies for later check.
    !-----------------------------------------------------------------------
    open(10,file='Fields_RTF.dat')
    write(10,1001)
1001 format('# RTF-Densities and RTF-fields'/ &
          & '#', 5x,'r',5x,'Pf_p',5x,'Pf_n',5x,'rhov_p',5x,'rhov_n', &
          & 5x,'Sigma',5x,'Omega',5x,'Rho',5x,'Coulomb', &
          & 5x,'Ef_p',5x,'Ef_n')
    do j=0,2000
       pf    = (1.5*pi**2*(RTF_dens(j,1)+RTF_dens(j,2)))**0.333333*197.33
       pf_p  = (3.*pi**2*RTF_dens(j,1))**0.333333*197.33
       pf_n  = (3.*pi**2*RTF_dens(j,2))**0.333333*197.33
       mstar = 938.-sigmaPot(j)
!        Ef    = sqrt(pf**2+mstar**2)+omegaPot(j)+rhoPot(j)
       Ef_p  = sqrt(pf_p**2+mstar**2)+omegaPot(j)+rhoPot(j)
       Ef_n  = sqrt(pf_n**2+mstar**2)+omegaPot(j)-rhoPot(j)

       write(10,1000) x_new(j),pf_p,pf_n,RTF_dens(j,:), &
            & sigmaPot(j)/1000.,omegaPot(j)/1000., &
            & rhoPot(j)/1000.,coulombPot(j)/1000., &
            & (Ef_p+coulombPot(j)-938.)/1000., (Ef_n-938.)/1000.
    end do

1000 format(40e15.5)
    close(10)

  !****************************************************************************
  end subroutine NucExRTF_Main !********************************************
  !****************************************************************************


end module NucExRTF
