!******************************************************************************
!****m* /mediumModule
! NAME
! module mediumModule

! PURPOSE
! Administrates the setup of the nuclear medium.
!******************************************************************************
module mediumModule

  use mediumDefinition

  implicit none
  private

  !****************************************************************************
  !****g* mediumModule/mediumCutOff
  ! SOURCE
  !
  real :: mediumCutOff=1.E-8
  ! PURPOSE
  ! If the density is lower than this value, then we treat the medium
  ! like vacuum.
  !****************************************************************************

  logical, save :: initFlag=.true.

  !****************************************************************************
  !****f* mediumModule/mediumAt1
  ! NAME
  ! type(medium) function mediumAt1(r)
  ! type(medium) function mediumAt2(density,r)
  ! PURPOSE
  ! Evaluate medium at some space-point r.
  ! INPUTS
  ! * real,dimension(1:3) :: r -- position where medium should be calculated
  ! * type(dichte) :: density -- density at this point
  ! NOTES
  ! if density not given, it is calculated via densityAt(r)
  !****************************************************************************
  Interface mediumAt
     Module Procedure mediumAt1,mediumAt2
  End Interface

  public :: getMediumCutOff,mediumAt

contains

  !****************************************************************************
  !****s* mediumModule/readInput
  ! NAME
  ! subroutine readInput
  ! PURPOSE
  ! Reads input in jobcard out of namelist "MediumModule".
  !****************************************************************************
  subroutine readInput
    use output, only: Write_ReadingInput

    integer :: ios

    !**************************************************************************
    !****n* mediumModule/MediumModule
    ! NAME
    ! NAMELIST MediumModule
    ! PURPOSE
    ! Includes the input switches:
    ! * mediumCutOff
    !**************************************************************************
    NAMELIST /mediumModule/ mediumCutOff

    call Write_ReadingInput('mediumModule',0)
    rewind(5)
    read(5,nml=mediumModule,IOSTAT=ios)
    call Write_ReadingInput('mediumModule',0,ios)

    write(*,*) 'Medium Cut-Off =', mediumCutOFF
    call Write_ReadingInput('mediumModule',1)

    initFlag=.false.
  end subroutine readInput


  real function getMediumCutOff()
    if (initFlag) call ReadInput
    getMediumCutOff=mediumCutOff
  end function


  !****************************************************************************
  ! cf. interface mediumAt
  !****************************************************************************
  type(medium) function mediumAt1(r)
    use dichteDefinition
    use densityModule, only: densityAt

    real,dimension(1:3),intent(in) :: r         ! position where medium should be calculated
    type(dichte):: density

    density=densityAt(r)
    mediumAt1=mediumAt2(density,r)

  end function mediumAt1
  !-------------------------------------------------------------------------
  type(medium) function mediumAt2(density,r)
    use dichteDefinition
    use minkowski, only: abs4
    use thermoDynamics, only: temperatureAt
    use RMF, only: getRMF_flag

    type(dichte),intent(in):: density
    real,dimension(1:3),intent(in) :: r         ! position where medium should be calculated

    if (initFlag) call ReadInput

    mediumAt2%density=abs4(density%baryon)
    mediumAt2%densityProton=abs4(density%proton)
    mediumAt2%densityNeutron=abs4(density%neutron)

    if (mediumAt2%densityProton+mediumAt2%densityNeutron.gt.mediumCutOff .and. .not.getRMF_flag()) then
       mediumAt2%useMedium=.true.
       mediumAt2%temperature=temperatureAt(r)
    else
       mediumAt2%useMedium=.false.
    end if

  end function mediumAt2
  !****************************************************************************

end module mediumModule
