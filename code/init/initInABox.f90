!******************************************************************************
!****m* /initInABox
! NAME
! module initInABox
! PURPOSE
! Initializes nucleons for a calculation in a box of nuclear matter
! ("periodic boundary conditions").
!******************************************************************************
module initInABox

  implicit none

  private

  !****************************************************************************
  !****g* initInABox/neutron_Density
  ! SOURCE
  !
  real, save :: neutron_Density = 0.084
  ! PURPOSE
  ! * neutron Density [fm^-3]
  !****************************************************************************

  !****************************************************************************
  !****g* initInABox/proton_Density
  ! SOURCE
  !
  real, save :: proton_Density = 0.084
  ! PURPOSE
  ! * proton Density [fm^-3]
  !****************************************************************************

  !****************************************************************************
  !****g* initInABox/fermiMotion
  ! SOURCE
  !
  logical, save :: fermiMotion = .true.
  ! PURPOSE
  ! * true  = switch on  Fermi motion
  ! * false = switch off Fermi motion
  !****************************************************************************

  !****************************************************************************
  !****g* initInABox/temp
  ! SOURCE
  !
  real, save :: temp = 0.
  ! PURPOSE
  ! If fermiMotion is true, this switch determines the temperature (in GeV)
  ! used in the Fermi distribution.
  !****************************************************************************

  !****************************************************************************
  !****g* initInABox/energy_density
  ! SOURCE
  !
  real, save :: energy_density = 0.
  ! PURPOSE
  ! Energy density in GeV/fm^3. If a finite positive number is given, the box
  ! will be boosted to a frame with the given energy density.
  !****************************************************************************

  !****************************************************************************
  !****g* initInABox/standing_wave_number
  ! SOURCE
  !
  integer, save :: standing_wave_number = 0
  ! PURPOSE
  ! If this number is larger than zero, the initial density distribution will
  ! not be uniform, but is modulated with a standing wave in z direction.
  ! The given number determines the number of oscillations throughout the box.
  ! The amplitude of the oscillations is currently fixed to be 20% of the
  ! (average) nucleon density.
  !****************************************************************************

  public :: initializeInABox
  public :: BoostToEps

contains

  !****************************************************************************
  !****s* initInABox/initInput
  ! NAME
  ! subroutine initInput
  ! PURPOSE
  ! Reads input out of jobcard. Namelist 'initInABox'.
  !****************************************************************************
  subroutine initInput

    use output, only: Write_ReadingInput
    use dichteDefinition
    use densityModule, only: get_densitySwitch,set_densityInput,densityAt
    use callstack, only: traceBack

    type(dichte) :: density

    !**************************************************************************
    !****n* initInABox/initInABox
    ! NAME
    ! NAMELIST initInABox
    ! PURPOSE
    ! Includes the input parameters:
    ! * proton_Density
    ! * neutron_Density
    ! * fermiMotion
    ! * temp
    ! * energy_density
    ! * standing_wave_number
    !**************************************************************************
    NAMELIST /initInABox/ proton_Density, neutron_Density, fermiMotion, temp, energy_density, standing_wave_number

    call Write_ReadingInput('initInABox',0)
    rewind(5)
    read(5,nml=initInABox)
    write(*,*) '  proton  Density   [1/fm^3] =', proton_Density
    write(*,*) '  neutron Density   [1/fm^3] =', neutron_Density
    write(*,*) '  fermiMotion =', fermiMotion
    if (fermiMotion) &
      write(*,*) '  temperature =', temp
    if (energy_density > 0.) &
      write(*,*) '  Energy Density  [GeV/fm^3] =', energy_density
    if (standing_wave_number > 0) &
      write(*,*) '  wave number =', standing_wave_number

    ! Taking care of the density routine
    select case (get_densitySwitch())
    case (1)
       ! ok
    case (3)
       write(*,'(A)') 'WARNING : Set parameters densityInput_XXX in '// &
            'densityModule to the values which are given here!!'
       call set_densityInput(proton_density,neutron_density)
    case default
       write(*,*) 'Wrong densitySwitch for box-calculations:', get_densitySwitch()
       call TRACEBACK('Only 1 or 3 allowed!')
    end select

    density=densityAt((/0.,0.,0./))
    write(*,'(A,4F8.3)') '  proton  Density with densityAt=',density%proton
    write(*,'(A,4F8.3)') '  neutron Density with densityAt=',density%neutron

    call Write_ReadingInput('initInABox',1)

  end subroutine initInput

  !****************************************************************************
  !****s* initInABox/initializeInABox
  ! NAME
  ! subroutine initializeInABox(teilchen)
  ! PURPOSE
  ! Initialize nucleons in a box.
  !****************************************************************************
  subroutine initializeInABox(teilchen)
    use particleDefinition
    use IdTable, only: nucleon
    use random, only: rn, rnFlat, rnOmega
    use densityModule, only: gridsize
    use insertion, only: GarbageCollection
    use constants, only: mN, pi
    use output, only: Write_InitStatus

    type(particle), dimension(:,:), intent(inOut) :: teilchen

    integer :: numberNeutrons, numberProtons, i,k,index,offset, producedProtons
    real :: volume, z

    call Write_InitStatus('box of nucleons',0)
    call initInput

    ! volume: factor two for each dimension, because box goes from -gridsize to +gridsize in each direction
    volume = 8. * gridsize(1)*gridsize(2)*gridsize(3)
    ! Calculate the number of particles which are to be initialized
    numberNeutrons = NINT(volume*neutron_density)
    numberProtons  = NINT(volume*proton_density)
    write(*,*) ' Number of neutrons per ensemble=',numberNeutrons
    write(*,*) ' Number of protons  per ensemble=',numberProtons
    write(*,'(A,3F9.3)') '  Gridsize   =', gridsize(:)
    write(*,*) ' Volume of box= (8*Gridsize**3) = ', volume
    write(*,*) ' Number Ensembles             =',size(teilchen(:,1))
    write(*,*) ' Number Particles per Ensemble=',size(teilchen(1,:))

    do i=1,size(teilchen,dim=1)  ! Loop over all ensembles
       producedProtons=0
       offset=0
       do k=1,numberNeutrons+numberProtons   ! Loop over all particles in the box
          do ! Search for empty space in particle vector
             index=k+offset
             if (teilchen(i,index)%ID > 0) then
                offset=offset+1
             else if (index>size(teilchen,dim=2)) then
                write(*,*) 'Real particle vector too small. Stop in initializeInABox.'
                stop
             else
                exit
             end if
          end do
          call setToDefault(teilchen(i,index)) ! set teilchen to its default values
          call setNumber(teilchen(i,index))    ! give each particle a unique number
          ! each particle needs a different %event, otherwise no collisions can happen
          Teilchen(i,index)%event = -Teilchen(i,index)%number
          Teilchen(i,index)%ID=Nucleon
          Teilchen(i,index)%antiparticle=.false.
          Teilchen(i,index)%perturbative=.false.
          Teilchen(i,index)%productionTime=0.
          Teilchen(i,index)%mass=mN
          call chooseCharge()
          Teilchen(i,index)%position(1) = rnFlat(-1.,1.) * gridSize(1)
          Teilchen(i,index)%position(2) = rnFlat(-1.,1.) * gridSize(2)
          z = rnFlat(-1.,1.) * gridSize(3)
          if (standing_wave_number > 0) then
            ! density follows standing wave in z direction
            do while (rnFlat(0., 1.2) > 1. + 0.2*sin(pi*standing_wave_number*z/gridSize(3)))  ! rejection loop
              z = rnFlat(-1.,1.) * gridSize(3)  ! re-sample z coordinate
            end do
          end if
          Teilchen(i,index)%position(3) = z
          call chooseMomentum(z)
       end do
       if (producedProtons/=numberProtons) then
          write(*,*) 'Problem in initPhaseSpace', producedProtons, numberProtons
          stop
       end if
    end do

    call GarbageCollection(teilchen)

    call Write_InitStatus('box of nucleons',1)

  contains

    ! Choose randomly the charge of the nucleon
    subroutine chooseCharge
      real :: probabilityProton
      probabilityProton=float(numberProtons-producedProtons)/float(numberProtons+numberNeutrons-k+1)
      if (rn()<=probabilityProton) then
         Teilchen(i,index)%charge=1
         producedProtons=producedProtons+1
      else
         Teilchen(i,index)%charge=0
      end if
    end subroutine chooseCharge

    ! Choose momentum randomly, according to (local) Fermi sphere
    subroutine chooseMomentum(z)
      use constants, only: hbarc
      use distributions, only: Fermi
      real, intent(in) :: z
      real :: pFermi, p, dens, mu0, mu, E, E_max, p_max
      if (fermiMotion) then
        if (teilchen(i,index)%charge==1) then
          dens = proton_density
        else
          dens = neutron_density
        end if
        if (standing_wave_number > 0) then
          ! local density (according to standing wave in z direction)
          dens = dens * (1. + 0.2*sin(pi*standing_wave_number*z/gridSize(3)))
        end if
        pFermi = (3.*pi**2*dens)**(1./3.) * hbarc   ! Fermi momentum
        if (temp>0.) then
          ! positive temperature: sample from Fermi distribution (using rejection sampling)
          mu0 = sqrt(pFermi**2 + mN**2)             ! Fermi energy
          mu = mu0 * (1 - pi**2/12.*(temp/mu0)**2)  ! non-relativistic estimate
          E_max = mu + 10*temp                      ! arbitrary energy cutoff
          p_max = sqrt(E_max**2 - mN**2)            ! corresponding momentum
          do
            p = p_max * rn()**(1./3.)          ! choose momentum according to p**2 dp
            E = sqrt(p**2 + mN**2)             ! calculate energy
            if (rn() < Fermi(E,mu,temp)) exit  ! check Fermi distribution
          end do
        else
          ! zero temperature: random momentum in sphere of radius pFermi
          p = pFermi * rn()**(1./3.)
        end if
        Teilchen(i,index)%momentum(1:3) = p * rnOmega()
      else
        ! particles at rest
        p = 0.
        Teilchen(i,index)%momentum(1:3) = 0.
      end if
      Teilchen(i,index)%momentum(0) = sqrt(p**2+mN**2)
      ! Assume vacuum dispersion relation:
      Teilchen(i,index)%velocity(1:3)=Teilchen(i,index)%momentum(1:3)/Teilchen(i,index)%momentum(0)
    end subroutine chooseMomentum


  end subroutine initializeInABox


  !****************************************************************************
  !****s* initInABox/BoostToEps
  ! NAME
  ! subroutine BoosToEps(Teilchen)
  ! PURPOSE
  ! Boost the particles along z-axis such that the energy density (=E/V)
  ! corresponds to the input parameter.
  ! This is only done if a finite positive energy density is given.
  !****************************************************************************
  subroutine BoostToEps(Teilchen)
    use particleDefinition
    use random, only: rn
    use densityModule, only: gridsize
    use constants, only: mN
    use lorentzTrafo, only: lorentz
    use output, only: Write_InitStatus

    type(particle), dimension(:,:), intent(inOut) :: Teilchen

    integer :: i,j,nEns
    real :: betaCM, SumE, betaV(1:3)

    ! only do the boost if energy_density is positive
    if (energy_density <= 0.) return

    call Write_InitStatus('BoostToEps',0)

    betaCM = ( (proton_Density+neutron_Density)*mN / energy_density )**2
    if (betaCM >= 1.0) then
       write(*,*) 'betaCM invalid! stop!'
       stop
    end if
    betaCM = sqrt(1-betaCM)

    write(*,*) 'betaCM = ',abs(betaCM)
    write(*,*) ' ---> gamma = ', 1./sqrt(1-betaCM**2)

    SumE = 0.0
    nEns = size(Teilchen,dim=1)
    do i=1,nEns  ! Loop over all ensembles
       do j=1,size(Teilchen,dim=2)  ! Loop over all particles
          if (Teilchen(i,j)%ID <= 0) cycle
          if (rn() > 0.5) then
             betaV = (/0., 0.,  betaCM/)
          else
             betaV = (/0., 0., -betaCM/)
          end if
          call lorentz(betaV,Teilchen(i,j)%momentum, "BoostToEps")

          SumE = SumE + Teilchen(i,j)%momentum(0)

       end do
    end do

    write(*,*) 'E/V = ',SumE/(nEns* 8*gridsize(1)*gridsize(2)*gridsize(3))

    call Write_InitStatus('BoostToEps',1)

  end subroutine BoostToEps


end module initInABox
