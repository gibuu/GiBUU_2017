!******************************************************************************
!****m* /selfEnergy_mesons
! NAME
! module selfEnergy_Mesons
!
! PURPOSE
! Includes the self energy (real and imaginary part) of the mesons
!******************************************************************************
module selfEnergy_mesons
  use tabulation
  use mediumDefinition
  implicit none
  private


  !****************************************************************************
  !****g* selfEnergy_mesons/dispersion
  ! SOURCE
  logical, save :: dispersion = .false.
  !
  ! PURPOSE
  ! Use dispersive real parts of the self energy.
  !****************************************************************************


  integer, parameter :: ID_min = 103, ID_max = 107
  type(table_2dim),save,dimension(ID_min:ID_max) :: realPartTable

  integer, save :: particleID
  type(medium), save :: med
  real, save :: pole

  public :: get_imagPart, get_realPart, tabulate_realPart, getDispersion

  logical, save :: first = .true.


contains


  logical function getDispersion()
    if (first) call readInput
    getDispersion = dispersion
  end function


  subroutine readInput
    use output, only: Write_ReadingInput, intToChar
    use inputGeneral, only: path_to_input

    !**************************************************************************
    !****n* selfEnergy_mesons/selfEnergyMesons
    ! NAME
    ! NAMELIST /selfEnergyMesons/
    ! PURPOSE
    ! Includes the switches:
    ! * dispersion
    !**************************************************************************
    NAMELIST /selfEnergyMesons/ dispersion
    character(200) :: filename
    integer :: IOS

    ! read namelist
    call Write_ReadingInput('selfEnergyMesons',0)
    rewind(5)
    read(5,nml=selfEnergyMesons,IOSTAT=IOS)

    write(*,'(A,L8)') 'Dispersion relation             ?', dispersion

    call Write_ReadingInput('selfEnergyMesons',0,IOS)

    ! read tabulation files
    if (dispersion) then
      do particleID = ID_min, ID_max
         filename = trim(path_to_input)//'/inMediumWidth/RealPart_'//intToChar(particleID)//'.dat.bz2'
         if (.not. readTable (realPartTable(particleID), filename)) then
            write(*,*) "Error reading real part of self energy for ID = ", particleID
            stop
         end if
      end do
    end if

    first = .false.
  end subroutine


  real function get_imagPart (ID, mass, mediumAtPosition)
    use mesonWidthMedium, only: widthMesonMedium
    integer, intent(in) :: ID
    real, intent(in) :: mass
    type(medium), intent(in) :: mediumAtPosition
    real, parameter :: mom(0:3) = (/0.,0.,0.,0./) ! dummy momentum variable (momentum is unused)
    get_imagPart = mass * widthMesonMedium (ID, mass, mom, mediumAtPosition)
  end function


  real function get_realPart (ID, mass, mediumAtPosition)
    integer, intent(in) :: ID
    real, intent(in) :: mass
    type(medium), intent(in) :: mediumAtPosition
    real(4) :: x(1:2)
    if (first) call readInput
    if (dispersion) then
      x = (/mass, mediumAtPosition%density/)
      get_realPart = interpolate2 (realPartTable(ID),x)
    else
      get_realPart = 0.
    end if
  end function


  real function integrand(s)
    use mesonWidthMedium, only: widthMesonMedium
    real, intent(in) :: s ! = m**2
    real, parameter :: mom(0:3) = (/0.,0.,0.,0./) ! dummy momentum variable (momentum is unused)
    integrand = - sqrt(s) * widthMesonMedium(particleID, sqrt(s), mom, med)/(s-pole)
  end function


  real function calc_realPart(m,dens)
    use constants, only: pi
    use particleProperties, only: hadron
    use quadPack, only: qawc

    real, intent(in) :: m,dens
    real :: a,b,c,int_lo,int_hi,abserr
    integer :: neval,ierr
    real, save :: pole1,pole2
    real, parameter :: cutoff = 7. ! upper limit for cauchy integral, limited by width tabulation
    real, parameter :: abs_acc = -1.
    real, parameter :: rel_acc = 1E-5
    real, parameter :: eps = 1E-3

    med%density = dens
    med%densityProton = dens/2.
    med%densityNeutron = dens/2.
    med%useMedium = .true.
    pole1 = min(m**2,hadron(particleID)%mass**2)   ! lower pole
    pole2 = max(m**2,hadron(particleID)%mass**2)   ! upper pole

    ! integration boundaries
    a = 0. ! meson(ID)%minMass**2
    b = (pole1 + pole2) / 2.
    c = cutoff**2

    if (m>cutoff-eps .or. abs(m-hadron(particleID)%mass)<eps) then ! m < meson(ID)%minmass+eps .or.
      calc_realPart = 0.
      return
    end if

    ! integrate over lower pole
    pole = pole2
    call qawc ( integrand, a, b, pole1, abs_acc, rel_acc, int_lo, abserr, neval, ierr )
    if (ierr /= 0) then
      write(*,*) "selfEnergy_mesons, QAWC1: ",particleID,m,ierr
    end if

    ! intergrate over upper pole
    pole = pole1
    call qawc ( integrand, b, c, pole2, abs_acc, rel_acc, int_hi, abserr, neval, ierr )
    if (ierr /= 0) then
      write(*,*) "selfEnergy_mesons, QAWC2: ",particleID,m,ierr
    end if

    calc_realPart = (m**2-hadron(particleID)%mass**2) * (int_lo+int_hi) / pi

  end function


  subroutine tabulate_realPart()
    use output, only: intToChar

    real, dimension (1:2), save :: mini,maxi,delta
    real, parameter :: max_rho    = 0.2  ! in fm**-3
    real, parameter :: max_mass   = 3.0  ! in GeV
    real, parameter :: delta_rho  = 0.02
    real, parameter :: delta_mass = 0.005

    mini = (/0.,0./)
    maxi = (/max_mass, max_rho/)
    delta = (/delta_mass, delta_rho/)

    do particleID = ID_min,ID_max
      realPartTable(particleID) = init_table_2dim (mini,maxi,delta,calc_realPart)
      call printTable(realPartTable(particleID),'RealPart_'//intToChar(particleID)//'.dat.bz2')
    end do

  end subroutine


end module
