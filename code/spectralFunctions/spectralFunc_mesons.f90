!******************************************************************************
!****m* /spectralFuncMesons
! NAME
! module spectralFuncMesons
!
! PURPOSE
! Includes the spectral functions of the mesons.
!
! NOTES
! Public routines: specFuncMes.
!******************************************************************************
module spectralFuncMesons
  implicit none
  private

  !****************************************************************************
  !****g* spectralFuncMesons/relativistic
  ! SOURCE
  logical, save :: relativistic = .true.
  ! PURPOSE
  ! * Use either relativistic or non relativistic spectral functions.
  !****************************************************************************

  !****************************************************************************
  !****g* spectralFuncMesons/debugFlag
  ! SOURCE
  logical, parameter :: debugFlag=.false.
  ! PURPOSE
  ! Switch for debug information
  !****************************************************************************

  logical, save :: initflag=.true.

  public :: specFuncMes

contains


  subroutine readinput
    use output, only: Write_ReadingInput
    integer :: IOS

    !**************************************************************************
    !****n* spectralFuncMesons/spectralFunctionMesons
    ! NAME
    ! NAMELIST /spectralFunctionMesons/
    ! PURPOSE
    ! Includes the switches:
    ! * relativistic
    !**************************************************************************
    NAMELIST /spectralFunctionMesons/ relativistic

    call Write_ReadingInput('spectralFunctionMesons',0)
    rewind(5)
    read(5,nml=spectralFunctionMesons,IOSTAT=IOS)

    write(*,'(A,L8)') 'Relativistic spectral functions ?', relativistic

    call Write_ReadingInput('spectralFunctionMesons',0,IOS)

    initFlag=.false.
  end subroutine readinput


  !****************************************************************************
  !****f* spectralFuncMesons/specFuncMes
  ! NAME
  ! real function specFuncMes(part,bareMass_out)
  !
  ! PURPOSE
  ! Returns the spectral function of a meson. Normalized such that integrating over m**2 gives 1.
  !
  ! INPUTS
  ! * type(particle)          :: part
  ! OUTPUT
  ! * real, optional          :: bareMass_out
  !
  !****************************************************************************
  real function specFuncMes(part,bareMass_out)
    use constants, only: pi
    use particleDefinition
    use IdTable, only: isMeson
    use particleProperties, only: hadron
    use minkowski, only: sp
    use mediumDefinition
    use mediumModule, only: mediumAt
    use mesonWidthMedium, only: WidthMesonMedium
    use potentialModule, only: scapot
    use selfenergy_mesons, only: get_realPart

    type(particle), intent(in) :: part
    real, optional, intent(out) :: bareMass_out  ! bare mass of the particle
    real :: width, invMass, invMassSquared, baremass, scalarPotential, realPart
    logical :: success
    type(medium) :: mediumAtPosition

    if (initFlag) then
       call readInput
       initFlag=.false.
    end if

    ! Check input
    if (.not.isMeson(part%ID)) then
       write(*,*) 'The specFuncMes routine is only designed for meson spectral functions! Stop!'
       write(*,*) 'ID=',part%ID
       stop
    end if

    ! Invariant mass of the particle:
    invMassSquared=sp(part%momentum,part%momentum)
    if (invMassSquared.le.0.) then
       specFuncMes=0.
       if (present(baremass_out)) baremass_out=0.
       return
    else
       invMass=sqrt(invMassSquared)
    end if

    scalarPotential=scapot(part,baremass,success)
    mediumAtPosition=mediumAt(part%position)
    width=WidthMesonMedium(part%ID,bareMass,part%momentum,mediumATposition)

    realPart = get_realPart(part%ID, bareMass, mediumAtPosition)

    if (relativistic) then
       specfuncMes=width/pi*invMass/((invMass**2-(hadron(part%ID)%mass+scalarPotential)**2-realPart)**2 + invMass**2*width**2)
    else
       specfuncMes=1./(4.*pi)/(hadron(part%ID)%mass+scalarPotential)*width/ &
                   ((invMass-(hadron(part%ID)%mass+scalarPotential))**2+width**2/4.)
    end if

    if (.not.success) then
       write(*,*) 'SpecFuncMes: No succes in scapot',part%ID,part%charge,part%momentum,part%position,baremass
       write(*,*) 'Setting specFuncMes=0.!'
       specfuncMes=0.
    end if

    if (present(baremass_out)) baremass_out=baremass

  end function specFuncMes


end module spectralFuncMesons
