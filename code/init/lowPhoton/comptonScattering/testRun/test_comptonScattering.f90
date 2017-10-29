
!*************************************************************************
!****p* /test_comptonScattering
! NAME
!  program test_comptonScattering
!
! PURPOSE
! * Tests the module "comptonScattering".
! * Prints out files named 
!   compton_"theta_gamma"_deg_"nucleon charge"_charge_cutoff_"cut_off"_GEV_CMframe.dat, e.g.
!   "compton_090_deg_001charge_cutoff1296GEV_CMframe.dat"
!   which give the cross section dsigma/dOmega_CM in [nb/sr] as a function of E_gamma. 
!   Several files obtained with different cut-off values are generated. 
! * This program is suited to fix Lambda
! * It can be steered via the namelist test_comptonScattering_nl
!*************************************************************************
program test_comptonScattering
  use comptonScattering
  use inputGeneral, only: readinputGeneral
  use particleProperties, only: initParticleProperties
  use output
  use formFactor_ResProd,  only : ff_resProd_setCutoff
  implicit none
  integer, save :: charge=0             ! Nucleon charge
  real   , save :: theta_gamma=130.7    ! Photon scattering angle
  integer, save :: frame=1              ! 1=photon angle in lab frame, 2=photon angle in CM frame
  
  real :: delta_egamma,e_gamma
  real :: cut_off

  call initParticleProperties
  call readinputGeneral
  call readinput()

  delta_egamma=0.03


  ! Loop over cut-offs
  cut_off=0.5
  cutoff_loop: do

     cut_off=cut_off*1.1
     if(cut_off.gt.2) exit

     open(200,file='compton_'//realToChar(theta_gamma)//'_deg_'//realToChar(float(charge))//'charge_cutoff'//&
          & realToChar4(cut_off*1000.) // 'GEV_CMframe.dat')
     call ff_resProd_setCutoff(cut_off)
     e_gamma=0.
     ! Loop over photon energy in lab frame
     egamma_loop: do 
        write(*,*) 'E_gamma=',e_gamma
        e_gamma=e_gamma+delta_egamma
        if(e_gamma.gt.0.8) exit
        write(200,*) e_gamma, theta_gamma,compton_crossSection_dOmega(e_gamma,theta_gamma,charge,2)*1E6
     end do egamma_loop

  end do cutoff_loop

  close(200)

contains

  subroutine readInput
    use output

    implicit none
    integer :: ios

    !***************************************************************************
    !****n* test_comptonScattering/test_comptonScattering_nl
    ! NAME
    ! NAMELIST /test_comptonScattering_nl/
    ! PURPOSE
    ! This Namelist for program "test_comptonScattering" includes:
    ! * charge
    ! * theta_gamma
    !***************************************************************************
    NAMELIST /test_comptonScattering_nl/ charge, theta_gamma, frame

    call Write_ReadingInput('test_comptonScattering_nl',0)

    rewind(5)
    read(5,nml=test_comptonScattering_nl,IOSTAT=ios)
    call Write_ReadingInput('test_comptonScattering_nl',0,ios)

    write(*,'(A20,L4)')    '  Charge=            ',Charge
    write(*,'(A20,F15.4)') '  Theta (degrees)=   ',Theta_gamma
    write(*,'(A20,I4)')    '  FRAME=             ',frame

    call Write_ReadingInput('test_comptonScattering_nl',1)

  end subroutine readInput


end program test_comptonScattering
