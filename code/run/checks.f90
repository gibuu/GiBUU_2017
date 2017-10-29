!******************************************************************************
!****m* /checks
! NAME
! module check
! PURPOSE
! Here some routines are collected, which do various checks, e.g. calculate the
! binding energy from the particle vector etc.
!******************************************************************************
module checks

  implicit none
  private

  real,save :: minimalRadius=1000.
  real,save :: maximalRadius=-1000.

  logical :: initFlag=.true.

  public :: checkRadius
  public :: evaluateTimeStep
  public :: CheckGridSize
  public :: evaluateTotal4Momentum_RMF
  public :: ChecksSetDefaulSwitches
  public :: ChecksSwitchRealRun
  public :: ChecksCallAll
  public :: ChecksCallEnergy

  !****************************************************************************
  !****g* checks/Do_CheckDensity
  ! SOURCE
  logical, save :: Do_CheckDensity = .false.
  ! PURPOSE
  ! Flag to indicate whether the density check routine should be called.
  !****************************************************************************

  !****************************************************************************
  !****g* checks/Do_CheckCoulomb
  ! SOURCE
  logical, save :: Do_CheckCoulomb = .false.
  ! PURPOSE
  ! Flag to indicate whether the Coulomb check routine should be called.
  !****************************************************************************

  !****************************************************************************
  !****g* checks/Do_CheckFermiSurface
  ! SOURCE
  logical, save :: Do_CheckFermiSurface = .false.
  ! PURPOSE
  ! Flag to indicate whether the Fermi-surface check routine should be called.
  !****************************************************************************

  !****************************************************************************
  !****g* checks/Do_CheckRadius
  ! SOURCE
  logical, save :: Do_CheckRadius = .false.
  ! PURPOSE
  ! Flag to indicate whether the radius check routine should be called.
  !****************************************************************************

  !****************************************************************************
  !****g* checks/Do_CheckMomentumDensity
  ! SOURCE
  logical, save :: Do_CheckMomentumDensity = .false.
  ! PURPOSE
  ! Flag to indicate whether the momentum-density check routine should be called.
  !****************************************************************************

  !****************************************************************************
  !****g* checks/Do_CheckEnergyLDA
  ! SOURCE
  logical, save :: Do_CheckEnergyLDA = .false.
  ! PURPOSE
  ! Flag to indicate whether the local density approximation check routine
  ! should be called.
  !****************************************************************************

  !****************************************************************************
  !****g* checks/Do_OccupiedReal
  ! SOURCE
  logical, save :: Do_OccupiedReal = .false.
  ! PURPOSE
  ! Flag to indicate whether the occupation check routine for the real particle
  ! vector should be called.
  !****************************************************************************

  !****************************************************************************
  !****g* checks/Do_OccupiedPert
  ! SOURCE
  logical, save :: Do_OccupiedPert = .false.
  ! PURPOSE
  ! Flag to indicate whether the occupation check routine for the perturbative
  ! particle vector should be called.
  !****************************************************************************

  !****************************************************************************
  !****g* checks/Do_CheckEnergy
  ! SOURCE
  logical, save :: Do_CheckEnergy = .false.
  ! PURPOSE
  ! Flag to indicate whether the energy check routine should be called.
  !****************************************************************************

  !****************************************************************************
  !****g* checks/Do_TachyonsReal
  ! SOURCE
  logical, save :: Do_TachyonsReal = .false.
  ! PURPOSE
  ! Flag to indicate whether the tachyon check routine for the real
  ! particle vector should be called.
  !****************************************************************************

  !****************************************************************************
  !****g* checks/Do_TachyonsPert
  ! SOURCE
  logical, save :: Do_TachyonsPert = .false.
  ! PURPOSE
  ! Flag to indicate whether the tachyon check routine for the perturbative
  ! particle vector should be called.
  !****************************************************************************

  !****************************************************************************
  !****g* checks/TachyonIsBlocking
  ! SOURCE
  logical, save :: TachyonIsBlocking = .false.
  ! PURPOSE
  ! Select whether the occurrence of a 'tachyon' in the check routines will stop
  ! the code or not (error messages will hopefully occur later in the code).
  !****************************************************************************

  !****************************************************************************
  !****g* checks/Do_CheckPertFlag
  ! SOURCE
  logical, save :: Do_CheckPertFlag = .true.
  ! PURPOSE
  ! Flag to indicate whether the flag '%perturbative' is set correctly in the
  ! particle vectors
  !****************************************************************************

  !****************************************************************************
  !****g* checks/Do_CheckConservation
  ! SOURCE
  logical, save :: Do_CheckConservation = .false.
  ! PURPOSE
  ! Flag to indicate whether conservation of energy/momentum, baryon number and
  ! strangeness between time steps for the real particles should be checked
  !****************************************************************************


  logical :: checkGS_Flag = .false. ! true: only for time=0, false: for all times

contains


  !****************************************************************************
  !****s* checks/ChecksSetDefaulSwitches
  ! NAME
  ! subroutine ChecksSetDefaulSwitches(eventtype)
  ! PURPOSE
  ! Set the default values of the switches "Do_Check*" according to the
  ! given eventtype.
  ! NOTES
  ! We repeat the specific codes snippet from GiBUU.f90, since not all
  ! conditions have been transported. Delete the commented lines after
  ! fulfilment.
  !****************************************************************************
  subroutine ChecksSetDefaulSwitches(eventtype)
    use EventTypes
    use RMF, only: getRMF_flag

    integer, intent(in) :: eventtype

    select case (eventtype)

    case (elementary,HeavyIon,hadron)
       Do_OccupiedReal = .true.

    case (inAbox)
       Do_OccupiedReal = .true.
!       Do_CheckDensity = .true.

    case (inAbox_pion,inAbox_delta)
       Do_CheckDensity = .true.
       Do_CheckCoulomb = .true.
       Do_OccupiedPert = .true.

    case (Box)
       Do_OccupiedReal = .true.

    case (LoPion,RealPhoton,LoLepton,neutrino)
       Do_checkMomentumDensity = .true.
       checkGS_Flag = .true.
       ! if CheckGS:Flag then...
       ! Do_CheckDensity = .true.
       ! Do_CheckCoulomb = .true.
       ! Do_CheckFermiSurface = .true.
       Do_OccupiedPert = .true.
       Do_TachyonsPert = .true.
       Do_TachyonsReal = .true.

    case (HiPion)
       Do_OccupiedPert = .true.
       Do_TachyonsPert = .true.
       Do_TachyonsReal = .true.

    case (groundState)
       Do_CheckDensity = .true.
       ! Do_CheckCoulomb = .true.
       Do_CheckRadius = .true.
       Do_checkMomentumDensity = .true.
       Do_checkEnergyLDA = .true.
       if (.not.getRMF_flag()) then
          Do_CheckFermiSurface = .true.
          Do_CheckEnergy = .true.
       end if
       Do_OccupiedReal = .true.

    case (ExternalSource)
       Do_OccupiedReal = .true.

    case default
       Do_OccupiedPert = .true.

    end select

    if (getRMF_flag()) then
       if (eventType.ne.elementary) then
          Do_CheckEnergy = .true.
       end if
    end if

!!$       select case (eventtype)
!!$
!!$       case (elementary,HeavyIon)
!!$
!!$          call occupied(realParticles,'realParticles') ! check occupation of the vector
!!$
!!$       case (inAbox_pion,inAbox_delta)
!!$
!!$          fName='Density_'//intTochar(timestep)//'.dat'
!!$          call checkDensity(realparticles,fName)
!!$          fName='Coulomb_'//intTochar(timestep)//'.dat'
!!$          call checkCoulomb(fName)
!!$          call occupied(pertParticles,'pertParticles') ! check occupation of the vector
!!$
!!$
!!$       case (LoPion,RealPhoton,neutrino)
!!$          call occupied(pertParticles,'pertParticles') ! check occupation of the vector
!!$          if (checkGS_Flag) then !in frozen mode check only at beginning ground state
!!$             !fName='Momentum_'//intTochar(timestep)//'.dat'
!!$             !call checkMomentumDensity(realparticles,fName)
!!$             !call checkRadius(realparticles,time)
!!$             fName='Density_'//intTochar(timestep)//'.dat'
!!$             call checkDensity(realparticles,fName)
!!$             fName='Coulomb_'//intTochar(timestep)//'.dat'
!!$             call checkCoulomb(fName)
!!$             fName='FermiSurface_'//intTochar(timestep)//'.dat'
!!$             call evaluateFermiSurface(fName)
!!$             checkGS_Flag = .false.
!!$          end if
!!$
!!$       case (HiPion)
!!$
!!$          !          fName='Density_'//intTochar(timestep)//'.dat'
!!$          !          call checkDensity(realparticles,fName)
!!$
!!$          call occupied(pertParticles,'pertParticles')
!!$
!!$       case (groundState)
!!$          fName='Density_'//intTochar(timestep)//'.dat'
!!$          call checkDensity(realparticles,fName)
!!$!          fName='Coulomb_'//intTochar(timestep)//'.dat'
!!$!          call checkCoulomb(fName)
!!$          call occupied(pertParticles,'pertParticles') ! check occupation of the vector
!!$          call checkRadius(realparticles,time)
!!$          fName='Momentum_'//intTochar(timestep)//'.dat'
!!$          call checkMomentumDensity(realparticles,fName)
!!$
!!$          select case (getPotentialEQSType())
!!$          case (6)
!!$             call checkEnergyLDA(realparticles)
!!$          case (8)
!!$             call checkEnergyLDAWelke(realparticles)
!!$          end select
!!$
!!$
!!$!          fName='RelMom_Deuterium_'//intTochar4(timestep)//'.dat'
!!$!          open(100,file=trim(fname))
!!$!          fName='RelPos_Deuterium_'//intTochar4(timestep)//'.dat'
!!$!          open(101,file=trim(fname))
!!$!          call deuteriumPL_getSpectra(deuterium_pointerList,100,101,time)
!!$!          close(100)
!!$!          close(101)
!!$
!!$          if( .not.getRMF_flag() ) then
!!$             fName='FermiSurface_'//intTochar(timestep)//'.dat'
!!$             call evaluateFermiSurface(fName)
!!$          end if
!!$
!!$       case default
!!$
!!$          call occupied(pertParticles,'pertParticles') ! check occupation of the vector
!!$
!!$       end select

  end subroutine ChecksSetDefaulSwitches


!!$  !*************************************************************************
!!$  !****s* checks/ChecksSwitchOffFrozen
!!$  ! NAME
!!$  ! subroutine ChecksSwitchOffFrozen()
!!$  ! PURPOSE
!!$  ! After having set all defaults and after the routines have been called
!!$  ! once at the first time step, we switch off
!!$  ! all checks just oncerning the realParticles vector, if
!!$  ! we have selected the frozen approximation.
!!$  !*************************************************************************
!!$  subroutine ChecksSwitchOffFrozen()
!!$
!!$    Do_OccupiedReal = .false.
!!$    Do_CheckDensity = .false.
!!$    Do_CheckRadius = .false.
!!$    Do_checkMomentumDensity = .false.
!!$    Do_checkEnergyLDA = .false.
!!$
!!$  end subroutine ChecksSwitchOffFrozen

  !****************************************************************************
  !****s* checks/ChecksSwitchRealRun
  ! NAME
  ! subroutine ChecksSwitchRealRun(isReal)
  ! PURPOSE
  ! Set Do_OccupiedReal,Do_OccupiedPert
  !****************************************************************************
  subroutine ChecksSwitchRealRun(isReal)

    logical, intent(in) :: isReal

    if (isReal) then
       Do_OccupiedReal = .true.
       Do_OccupiedPert = .false.
    else
       Do_OccupiedReal = .false.
       Do_OccupiedPert = .true.
    end if

  end subroutine ChecksSwitchRealRun

  !****************************************************************************
  !****s* checks/ChecksSwitchGS
  ! NAME
  ! subroutine ChecksSwitchGS(isFrozen)
  ! PURPOSE
  ! Set the Flag, whether Ground State checks are only done at timestep 0
  ! or for all timesteps.
  !****************************************************************************
!   subroutine ChecksSwitchGS(isFrozen)
!     logical, intent(in) :: isFrozen
!
!     checkGS_Flag = isFrozen
!   end subroutine ChecksSwitchGS

  !****************************************************************************
  !****s* checks/Init
  ! NAME
  ! subroutine init()
  ! PURPOSE
  ! Read the namelist "Checks".
  !****************************************************************************
  subroutine Init
    use output

    integer :: ios

    !**************************************************************************
    !****n* checks/Checks
    ! NAME
    ! NAMElist /Checks/
    ! PURPOSE
    ! Includes the following switches:
    ! * Do_CheckDensity
    ! * Do_CheckCoulomb
    ! * Do_CheckFermiSurface
    ! * Do_CheckRadius
    ! * Do_CheckMomentumDensity
    ! * Do_CheckEnergyLDA
    ! * Do_OccupiedReal
    ! * Do_OccupiedPert
    ! * Do_CheckEnergy
    ! * Do_TachyonsReal
    ! * Do_TachyonsPert
    ! * TachyonIsBlocking
    ! * Do_CheckPertFlag
    ! * Do_CheckConservation
    !**************************************************************************
    NAMELIST /checks/ Do_CheckDensity, Do_CheckCoulomb, Do_CheckFermiSurface, &
         Do_CheckRadius, Do_CheckMomentumDensity, Do_CheckEnergyLDA, &
         Do_OccupiedReal, Do_OccupiedPert, Do_CheckEnergy, &
         Do_TachyonsReal, Do_TachyonsPert, TachyonIsBlocking, &
         Do_CheckPertFlag, Do_CheckConservation

    call Write_ReadingInput('checks',0)
    rewind(5)
    read(5,nml=checks,iostat=ios)
    call Write_ReadingInput('checks',0,ios)

    write(*,'(A,L1)') ' Do_CheckDensity         = ',Do_CheckDensity
    write(*,'(A,L1)') ' Do_CheckCoulomb         = ',Do_CheckCoulomb
    write(*,'(A,L1)') ' Do_CheckFermiSurface    = ',Do_CheckFermiSurface
    write(*,'(A,L1)') ' Do_CheckRadius          = ',Do_CheckRadius
    write(*,'(A,L1)') ' Do_CheckMomentumDensity = ',Do_CheckMomentumDensity
    write(*,'(A,L1)') ' Do_CheckEnergyLDA       = ',Do_CheckEnergyLDA
    write(*,'(A,L1)') ' Do_OccupiedReal         = ',Do_OccupiedReal
    write(*,'(A,L1)') ' Do_OccupiedPert         = ',Do_OccupiedPert
    write(*,'(A,L1)') ' Do_CheckEnergy          = ',Do_CheckEnergy
    write(*,'(A,L1)') ' Do_CheckPertFlag        = ',Do_CheckPertFlag
    write(*,'(A,L1)') ' Do_CheckConservation    = ',Do_CheckConservation
    write(*,'(A,L1)') ' Do_TachyonsReal         = ',Do_TachyonsReal
    write(*,'(A,L1)') ' Do_TachyonsPert         = ',Do_TachyonsPert
    write(*,*)
    write(*,'(A,L1)') ' TachyonIsBlocking       = ',TachyonIsBlocking

    call Write_ReadingInput('checks',1)
    initflag = .false.
  end subroutine Init


  !****************************************************************************
  !****s* checks/ChecksCallEnergy
  ! NAME
  ! subroutine ChecksCallEnergy(time, realParticles)
  ! PURPOSE
  ! This is one of the upper-level routines, called to check the energies.
  !
  ! Major switches: 'Do_CheckEnergy'
  !****************************************************************************
  subroutine ChecksCallEnergy(time, realParticles)
    use particleDefinition
    use RMF, only: getRMF_flag

    real, intent(in) :: time
    type(particle), dimension(:,:), intent(in) :: realParticles

    if (initflag) call init()

    if (.not. Do_CheckEnergy) return

    call evaluateEnergy(realParticles,time)
    if (getRMF_flag()) then
       call evaluateTotal4Momentum_RMF(realParticles,time)
    else
       call evaluateBindingEnergy_Teis(realParticles,time)
    end if

  end subroutine ChecksCallEnergy


  !****************************************************************************
  !****s* checks/ChecksCallOccupied
  ! NAME
  ! subroutine ChecksCallOccupied(realParticles,pertParticles,text)
  ! PURPOSE
  ! This is one of the upper level routines, called to check the occupansies
  ! of the real and perturbative particle vectors.
  !
  ! Major switches: 'Do_OccupiedReal', 'Do_OccupiedPert'
  !****************************************************************************
  subroutine ChecksCallOccupied(realParticles,pertParticles,text)
    use particleDefinition

    type(particle),dimension(:,:),intent(in)  :: realParticles
    type(particle),dimension(:,:),intent(in)  :: pertParticles
    character(*), intent(in) :: text

    ! e.g. text = 'At the end of the run: '

    if (initflag) call init()

    if (Do_OccupiedReal) then
       call occupied(realParticles,text//'realParticles')
    end if
    if (Do_OccupiedPert) then
       call occupied(pertParticles,text//'pertParticles')
    end if

  end subroutine ChecksCallOccupied

  !****************************************************************************
  !****s* checks/ChecksCallAll
  ! NAME
  ! subroutine ChecksCallAll(timestep,time,realParticles,pertParticles)
  ! PURPOSE
  ! This is one of the upper level routines, called to perform major checks
  ! as indicated by the flags 'Do_Check???'.
  ! NOTES
  ! * Some test are done by other routines, as e.g. 'ChecksCallOccupied'
  ! * This routine is called *before* 'GarbageCollection', thus you can not
  !   rely on the fact, that there are no holes in the particle vector.
  !****************************************************************************
  subroutine ChecksCallAll(timestep,time,realParticles,pertParticles)
    use particleDefinition
    use output, only: intTochar
    use baryonPotentialModule, only: getPotentialEQSType

    integer, intent(in) :: timestep
    real, intent(in) :: time
    type(particle),dimension(:,:),intent(in)  :: realParticles
    type(particle),dimension(:,:),intent(in)  :: pertParticles

    character(50) :: fName

    if (initflag) call init()

    if (Do_CheckDensity) then
       fName='Density_'//intTochar(timestep)//'.dat'
       call checkDensity(realParticles,fName)
    end if

    if (Do_CheckCoulomb) then
       fName='Coulomb_'//intTochar(timestep)//'.dat'
       call checkCoulomb(fName)
    end if

    if (Do_CheckFermiSurface) then
       fName='FermiSurface_'//intTochar(timestep)//'.dat'
       call evaluateFermiSurface(fName)
    end if

    if (Do_CheckRadius) then
       call checkRadius(realparticles,time)
    end if

    if (Do_checkMomentumDensity) then
       fName='Momentum_'//intTochar(timestep)//'.dat'
       call checkMomentumDensity(realparticles,timestep,fName)
    end if

    if (Do_checkEnergyLDA) then
       select case (getPotentialEQSType())
       case (6)
          call checkEnergyLDA(realparticles)
       case (8)
          call checkEnergyLDAWelke(realparticles)
       end select
    end if

    if (Do_CheckPertFlag) then
       call ChecksCallPertFlag(realParticles, .false.)
       call ChecksCallPertFlag(pertParticles, .true.)
    end if

    call ChecksCallOccupied(realParticles,pertParticles,"")

    if ((checkGS_Flag).and.(Do_TachyonsReal)) then
       call ChecksCallTachyon(realParticles,"realParticles:")
    end if
    if (Do_TachyonsPert) then
       call ChecksCallTachyon(pertParticles,"pertParticles:")
    end if

    if (Do_CheckConservation) then
       call CheckConservation(realParticles)
    end if

  end subroutine ChecksCallAll


  !****************************************************************************
  !****s* checks/evaluateFermiSurface
  ! NAME
  ! subroutine evaluateFermiSurface(fileName)
  ! PURPOSE
  !
  !****************************************************************************
  subroutine evaluateFermiSurface(fileName)
    use particleDefinition
    use dichteDefinition
    use densitymodule, only: get_densitySwitch, set_densitySwitch, densityAt, fermiMomentum_sym, fermiMomentum_noIsospin
    use coulomb, only: emfoca
    use minkowski, only: abs4
    use baryonPotentialModule
    use Idtable, only: nucleon
    use constants, only: pi, mN
    use mediumDefinition

    character(*), intent(in),optional :: fileName
    type(particle) :: p
    real :: rhoNuc,pot, pf, average_fs, average,densInt,dp,fsInt
    real :: Ef_p,Ef_n,rhoNuc_p,rhoNuc_n,pf_p,pf_n,cpot,Efree_p,Efree_n
    type(dichte) :: dens
    real, parameter :: dr=0.05
    integer :: i,j,densitySwitch_Save
    type(medium)         :: med

    p%ID=nucleon
    p%position=0.
    p%momentum=0.
    p%mass=mN

    med%temperature    =0.
    med%useMedium      =.true.

    if (present(fileName)) then
       open(100,File=filename)
    else
       open(100,File='FermiSurface_init.dat')
    end if

    write(100,'(A1,6A20)') &
         & "#","x", "E_f_p", "E_f_n", "E_f_p(rel)", "E_f_n(rel)","Density"

    if ( .not.present(fileName)) then

       densitySwitch_Save = get_densitySwitch()
       call set_densitySwitch(2)

       do i=0,1000
          p%position(1) = float(i)*dr
          dens          = densityAt(p%position(1))
          rhoNuc        = abs4(dens%proton+dens%neutron)
          rhoNuc_p      = abs4(dens%proton)
          rhoNuc_n      = abs4(dens%neutron)

          med%density        = rhoNuc
          med%densityProton  = rhoNuc_p
          med%densityNeutron = rhoNuc_n

          p%momentum(1) = fermiMomentum_sym(rhoNuc)
          pf_p          = fermiMomentum_noIsospin(rhoNuc_p)
          pf_n          = fermiMomentum_noIsospin(rhoNuc_n)
          pot           = BaryonPotential(p,med,.false.)
          cpot = emfoca(p%position,p%momentum,1, 1)
          Ef_p          = pf_p**2/(2*mN)+pot + cpot
          Ef_n          = pf_n**2/(2*mN)+pot
          p%momentum(1) = fermiMomentum_noIsospin(rhoNuc_p)
          Efree_p       = freeEnergy(p)
          p%momentum(1) = fermiMomentum_noIsospin(rhoNuc_n)
          Efree_n       = freeEnergy(p)
          write(100,'(36E20.9)') p%position(1), &
               & Ef_p,Ef_n,Efree_p+pot-mN+cpot,Efree_n+pot-mN, &
               & rhoNuc_p, rhoNuc_n,pot
       end do

       !return to original value:
       call set_densitySwitch(densitySwitch_Save)

       return ! Exit this routine
    end if

    do i=0,1000
       p%position(1) = float(i)*dr
       dens          = densityAt(p%position(1))
       rhoNuc        = abs4(dens%proton+dens%neutron)
       rhoNuc_p      = abs4(dens%proton)
       rhoNuc_n      = abs4(dens%neutron)
       med%density        = rhoNuc
       med%densityProton  = rhoNuc_p
       med%densityNeutron = rhoNuc_n
       p%momentum(1) = fermiMomentum_sym(rhoNuc)
       pf_p          = fermiMomentum_noIsospin(rhoNuc_p)
       pf_n          = fermiMomentum_noIsospin(rhoNuc_n)
       pot           = BaryonPotential(p,med,.false.)
       cpot = emfoca(p%position,p%momentum,1, 1)
       Ef_p          = pf_p**2/(2*mN)+pot + cpot
       Ef_n          = pf_n**2/(2*mN)+pot
       p%momentum(1) = fermiMomentum_noIsospin(rhoNuc_p)
       Efree_p       = freeEnergy(p)
       p%momentum(1) = fermiMomentum_noIsospin(rhoNuc_n)
       Efree_n       = freeEnergy(p)
       write(100,'(36E20.9)') p%position(1), &
            & Ef_p,Ef_n,Efree_p+pot-mN+cpot,Efree_n+pot-mN, &
            & rhoNuc_p, rhoNuc_n,pot
    end do

    average=0.
    densInt=0.
    do i=0,1000
       p%position(1)=float(i)*dr
       dens=densityAt(p%position)
       rhoNuc=abs4(dens%proton+dens%neutron)
       rhoNuc_p      = abs4(dens%proton)
       rhoNuc_n      = abs4(dens%neutron)
       med%density        = rhoNuc
       med%densityProton  = rhoNuc_p
       med%densityNeutron = rhoNuc_n
       p%momentum(1)=fermiMomentum_sym(rhoNuc)
       pot=BaryonPotential(p,med,.false.)
       average=average+(freeEnergy(p)+pot-mN)*4.*pi/3.*rhoNuc*p%position(1)**2/dr
       densInt=densInt+4.*pi/3.*rhoNuc*p%position(1)**2/dr
    end do

    if (densInt.gt.1E-10) &
         write(100,*) '# Average surface energy:', average/densInt


    average=0.
    densInt=0.
    do i=0,1000
       p%position(1)=float(i)*dr
       dens=densityAt(p%position)
       rhoNuc=abs4(dens%proton+dens%neutron)
       rhoNuc_p      = abs4(dens%proton)
       rhoNuc_n      = abs4(dens%neutron)
       med%density        = rhoNuc
       med%densityProton  = rhoNuc_p
       med%densityNeutron = rhoNuc_n
       pf=fermiMomentum_sym(rhoNuc)
       dp=pf/1000.
       average_fs=0.
       fsInt=0.
       if (pf.gt.0.001) then
          do j=0,1000
             p%momentum(1)=float(j)*dp
             pot=BaryonPotential(p,med,.false.)
             average_fs=average_fs+pot*4.*pi/3.*p%momentum(1)**2/dp
             fsInt     =fsInt     +    4.*pi/3.*p%momentum(1)**2/dp
          end do
          average=average+average_fs*4.*pi/3.*p%position(1)**2/dr
          densInt=densInt+     fsInt*4.*pi/3.*p%position(1)**2/dr
       end if

    end do


    if (densInt.gt.1E-10) &
         write(100,*) '# Average single-partilce binding energy:', &
         average/densInt


  end subroutine evaluateFermiSurface


  !****************************************************************************
  !****s* checks/particleTrack
  ! NAME
  ! subroutine particleTrack(p,name)
  ! PURPOSE
  ! Prints the track of a particle
  !****************************************************************************
!   subroutine particleTrack(p,name)
!     use particleDefinition
!     use output, only : intToChar
!     type(particle) :: p
!     character(*),optional :: name
!
!     if(present(name)) then
!        Open(100,file='particle_track_'//trim(name)//'.dat', position='Append')
!     else
!        Open(100,file='particle_track.dat', position='Append')
!     end if
!     write(100,'(I8,7F20.5)') p%ID, p%position , p%momentum
!   end subroutine particleTrack


  !****************************************************************************
  !****s* checks/particleTrack_byNumber
  ! NAME
  ! subroutine particleTrack_byNumber(parts,number,time,DoNew)
  ! PURPOSE
  ! Prints the track of a particle with given number
  !****************************************************************************
!   subroutine particleTrack_byNumber(parts,number,time,DoNew)
!     use particleDefinition
!     use output, only : intToChar4
!
!     type(particle),dimension(:,:),intent(in)  :: parts
!     integer, intent(in) :: number
!     real, OPTIONAL, intent(in) :: time
!     logical, OPTIONAL, intent(in) :: DoNew
!
!
!     integer :: i,j
! !    logical :: found
!     logical,save :: firstTime=.true.
!     real :: t
!
!     t = 0
!     if (present(time)) t=time
!
!     if (present(DoNew)) firstTime = DoNew
!
!     if(firstTime) then
!        Open(100,file='particle_track_'//trim(intToChar4(number))//'.dat')
!        write(100,*) '#time, ID, position(1:3), momentum(0:3),mass,velocity(1:3)'
!        firstTime=.false.
!     else
!        Open(100,file='particle_track_'//trim(intToChar4(number))//'.dat', position='Append')
!     end if
!
! !    found=.false.
!     do j=lbound(parts,dim=2),ubound(parts,dim=2)
!        do i=lbound(parts,dim=1),ubound(parts,dim=1)
!           if(parts(i,j)%number.eq.number.and.parts(i,j)%ID.gt.0) then
!              write(100,'(f15.5,I8,11F13.5)') t, parts(i,j)%ID, parts(i,j)%position , &
!                   & parts(i,j)%momentum,parts(i,j)%mass,parts(i,j)%velocity
! !             found=.true.
! !             write(*,*) 'Particle found!'
! !             write(*,'(I8,11F13.5)') parts(i,j)%ID, parts(i,j)%position , &
! !                  & parts(i,j)%momentum,parts(i,j)%mass,parts(i,j)%velocity
!           else
! !             parts(i,j)%ID=0
!           end if
!        end do
!     end do
! !    if(.not.found) stop'particleTrack_byNumber'
!     close(100)
!   end subroutine particleTrack_byNumber


  !****************************************************************************
  !****s* checks/occupied
  ! NAME
  ! subroutine occupied(particleVector,name)
  ! PURPOSE
  ! Prints information, how many entries and of every ensemble is actually
  ! occupied by particles and the percentage of maximum occupation
  !****************************************************************************

  !****************************************************************************
  !****o* checks/Occupation.dat
  ! NAME
  ! file Occupation.dat
  ! PURPOSE
  ! Show how many entries of every ensemble is actually occupied by particles.
  ! If it is close to the maximal possible (lengthPert  for perturbative runs), one should
  ! increase parameter length_Perturbative in the namelist "input".
  ! Columns:
  ! * 1: ensemble number
  ! * 2: occupation
  !****************************************************************************

  subroutine occupied(particleVector,name)
    use particleDefinition

    type(particle),intent(in),dimension(:,:) :: particleVector
    integer :: i,j,countTotal,count, over80count
    character(*), intent(in) :: name

    countTotal=0  ! Counts all ocuppied spaces
    over80count=0 ! Counts ensembles which are more than 80% occupied

!     Open(178,File='Occupation.dat')

    !  write(*,*)
!     Write(178,*) '# Maximal possible number of entries is', size(particleVector,dim=2)
!     Write(178,*) '# 1: ensemble-number  2: entries'
    do i=lbound(particleVector,dim=1),ubound(particleVector,dim=1)
       count=0
       do j=lbound(particleVector,dim=2),ubound(particleVector,dim=2)
         if (particleVector(i,j)%ID > 0) count=count+1
       end do
!        Write(178,*) i,count
       countTotal=countTotal+count
       if (float(count)/float(size(particleVector,dim=2))>0.8) over80count=over80count+1
!        if (float(count)/float(size(particleVector,dim=2))>0.95) write(*,*) '#### Warning: Ensemble #',i,' is over 95% occupied'
    end do
!     Close(178)
    write(*,'(A,A,A,F5.1,A,F5.1,A)') '#### ', name,' is occupied by ', &
                                     float(countTotal*100)/float(size(particleVector,dim=1)*size(particleVector,dim=2)), &
                                     '% with ', float(over80count*100)/float(size(particleVector,dim=1)), '% ensembles over 80%.'
  end subroutine occupied


  !****************************************************************************
  !****s* checks/checkPauli
  ! NAME
  ! subroutine checkPauli(realparticles)
  ! PURPOSE
  !
  !****************************************************************************
!   subroutine checkPauli(realparticles)
!     use pauliBlockingModule, only: pauliBlocking
!     use particleDefinition
!     use constants, only: mN
!
!     type(particle),intent(in),dimension(:,:) :: realparticles
!
!     integer :: i,k
!     real :: x,p,z!,y
!     real,dimension(0:3) :: momentum
!     logical :: pauliFlag
!     Write(*,*) 'Check pauli-blocking'
!     !  call timeMeasurement(.true.)
!     Print *,'Print CheckPauliXP.dat'
!     Open(11,File='CheckPauliXP.dat')
!     Do k=-25,25
!        p=k*0.04
!        momentum=(/SQRT(p**2+mN**2),p,0.,0./)
!        Do i=-25,25
!           x=i*0.5
!           pauliflag=pauliBlocking(momentum,(/x,0.,0./),0,realparticles)
!           If (pauliflag) then
!              write(11,*) x,p,1
!           else
!              write(11,*) x,p,0
!           end if
!        end do
!     end do
!     Close(11)
!     !  call timeMeasurement
!     !  call timeMeasurement(.true.)
!     Print *,'Print CheckPauliXZ.dat'
!     Open(11,File='CheckPauliXZ.dat')
!     Do k=-25,25
!        z=k*1.0
!        p=0.01
!        momentum=(/SQRT(p**2+mN**2),p,0.,0./)
!        Do i=-25,25
!           x=i*0.5
!           pauliflag=pauliBlocking(momentum,(/x,0.,z/),0,realparticles)
!           If (pauliflag) then
!              write(11,*) x,z,1
!           else
!              write(11,*) x,z,0
!           end if
!        end do
!     end do
!     Close(11)
!     !  call timeMeasurement
!     !  call timeMeasurement(.true.)
!     Do k=-25,25
!        z=k*1.0
!        p=0.01
!        momentum=(/SQRT(p**2+mN**2),p,0.,0./)
!        Do i=-25,25
!           x=i*0.5
!           pauliflag=pauliBlocking(momentum,(/x,0.,z/),0,realparticles)
!        end do
!     end do
!     !  call timeMeasurement
!
!   end subroutine checkPauli


  !****************************************************************************
  !****s* checks/checkMomentumDensity
  ! NAME
  ! subroutine checkMomentumDensity(realParticles,filename)
  ! PURPOSE
  !
  !****************************************************************************
  subroutine checkMomentumDensity(realParticles,timestep,filename)
    use particleDefinition
    use histf90
    use hist2Df90
    use dichteDefinition
    use densitymodule, only: densityAt
    use minkowski, only: abs4
    use idtable, only: nucleon

    type(particle), dimension(:,:) :: realParticles
    integer, intent(in) :: timestep
    character(20),intent(in),optional :: filename

    type(histogram), save   :: momentumHist
    type(histogram2D), save :: momDens
    real, save :: integral,mean_momentum
    logical, save :: DoInit = .true.

    real :: dens_nuc, p
    real, parameter :: dp=0.005 ! grid spacing for histogram
    type(dichte) :: dens
    integer :: i,j

    if ((checkGS_Flag).and.(timestep > 0)) return

    write(*,*) 'Check momentum density'

    if (DoInit) then
       integral=0.
       mean_momentum=0.

       call CreateHist2D(momDens,'momentum vs. density', &
            (/0.,0./),(/0.4,0.18/),(/dp,0.005/))
       call CreateHist(momentumHist, 'momentum density',0.,0.5,dp)
       DoInit = .false.
    end if

    do i= lbound(realParticles,dim=1),ubound(realParticles,dim=1)
       do j= lbound(realParticles,dim=2),ubound(realParticles,dim=2)
          if (realParticles(i,j)%ID.ne.nucleon) cycle
          p = absMom(realParticles(i,j))
          call AddHist(momentumHist,p ,1.)
          dens=densityAt(realParticles(i,j)%position)
          dens_nuc=abs4(dens%proton+dens%neutron)
          call AddHist2d(momDens,(/p,dens_nuc/) ,1.)

          integral      = integral     +1.
          mean_momentum = mean_momentum+p
       end do
    end do

    if (present(fileName)) then
       open(121,File=filename)
    else
       open(121,File='CheckMomentumDensity.dat')
    end if
    write(121,*) '#### Mean Momentum=', mean_momentum/integral
    call WriteHist(momentumHist,121,mul=1./integral)
    close(121)


    open(121,File='Mom_vs_Dens.dat')
    call WriteHist2D_Gnuplot(momDens,121,mul=1./integral)
    close(121)

    if (.not. checkGS_Flag) then
       ! prepare everything for next call:
       call removeHist(momentumHist)
       call removeHist2D(momDens)
       DoInit = .true.
    end if

  end subroutine checkMomentumDensity



  !****************************************************************************
  !****s* checks/checkDensity
  ! NAME
  ! subroutine checkDensity(realParticles,filename)
  ! PURPOSE
  !
  !****************************************************************************
  subroutine checkDensity(realParticles,filename)
    use densityModule, only: densityAt,get_densitySwitch,set_densitySwitch,&
         updateDensity,gridPoints,gridSpacing
    use dichtedefinition
    use particleDefinition

    type(particle), dimension(:,:) :: realParticles
    character(20),intent(in),optional :: filename
    integer :: i,j,k
    real :: x,y,z
    real , dimension(1:3) :: place
    real :: meanvalue
    type(dichte) :: density

    integer :: save_densitySwitch
    write(*,*) 'Check density'

    ! set density to dynamic
    save_densitySwitch=get_densitySwitch()
    call set_densitySwitch(1)

    if (save_densitySwitch/=1) call updateDensity(realParticles)

    if (present(fileName)) then
       open(121,File=filename)
    else
       open(121,File='CheckDensity.dat')
    end if
    do j=-30,30
       z=j*0.25
       do i=-30,30
          x=i*0.25
          place=(/x,0.,z/)
          density=densityAt(place)
          write(121,fmt='(5f8.3)') x,z,&
               density%baryon(0),density%proton(0),density%neutron(0)
       end do
       write(121,*)
    end do
    close(121)

    meanValue=0.
    do j=-gridPoints(3),gridPoints(3)-1
       z=float(j)*gridSpacing(3)
       do i=-gridPoints(1),gridPoints(1)-1
          x=float(i)*gridSpacing(1)
          do k=-gridPoints(2),gridPoints(2)-1
             y=float(k)*gridSpacing(2)
             place=(/x,y,z/)
             density=densityAt(place)
             meanvalue=meanvalue+density%baryon(0)
          end do
       end do
    end do
    write(*,*) 'meanvalue of baryon density=', &
         meanValue/float(2*gridPoints(1)*2*gridPoints(2)*2*gridPoints(3))

    if (present(fileName)) then
       open(121,File='ZAxis'//filename)
    else
       open(121,File='ZAxisCheckDensity.dat')
    end if
    do j=-80,80
       z=j*0.15
       place=(/0.,0.,z/)
       density=densityAt(place)
       write(121,fmt='(5f8.3)') z,&
            density%baryon(0),density%proton(0),density%neutron(0)
    end do
    close(121)
    ! Set densitySwitch back to original value

    call set_densitySwitch(save_densitySwitch)

  end subroutine checkDensity


  !****************************************************************************
  !****s* checks/checkCoulomb
  ! NAME
  ! subroutine CheckCoulomb(filename)
  ! PURPOSE
  !
  !****************************************************************************
  subroutine CheckCoulomb(filename)

    use coulomb, only: emfoca

    character(20),intent(in),optional :: filename

    real :: x,y,z
    real,dimension(1:3) :: r,p
    integer :: i,j

    if (present(fileName)) then
       open(11,File=filename)
       open(12,File=trim(filename)//'.Zaxis.dat')
    else
       open(11,File='CheckCoulomb.dat')
       open(12,File='CheckCoulomb_Zaxis.dat')
    end if

    y=0.
    do i=-25,25
       x=i*1.
       do j=-25,25
          z=j*1.
          r=(/x,y,z/)
          p=(/0.,0.,0./)
          write(11,fmt='(4f15.7)') x,z,emfoca(r,p,-1,1),evalCoulombPot((/x,y,z/))
       end do
       write(11,*)
    end do
    close(11)
    y=0.
    x=0.
    do j=-100,100
       z=j*0.3
       r=(/x,y,z/)
       p=(/0.,0.,0./)
       write(12,fmt='(4f15.7)') z,emfoca(r,p,-1,1),evalCoulombPot((/x,y,z/))
    end do

  contains
    ! Analytical solution for Coulomb potential, assuming that rho=rho(|r|)
    real function evalCoulombPot(r_vec)
      use densityModule
      use dichteDefinition
      use constants, only: pi
      use vector, only: absVec

      real, intent(in),dimension(1:3) :: r_vec
      type(dichte) :: dens
      real, dimension(1:3) :: rp_vec
      real :: r,rp,phi, rhoP
      real, parameter :: elmcon=0.001439
      real, parameter :: dr=0.01 ! Step size in radius for integration

      evalCoulombPot=0.
      phi=0.
      rp=0.
      r=absVec(r_vec)
      do
         rp=rp+dr
         rp_vec=(/rp,0.,0./)
         dens=densityAt(rp_vec)
         rhoP=dens%proton(0)
         if (rp.ge.r) then
            phi=phi+rhoP*rp
         else
            phi=phi+rhoP*rp**2/r
         end if
         if (rp.gt.20) exit
      end do

      evalCoulombPot=4*pi*phi* elmcon*dr

    end function evalCoulombPot

  end subroutine checkCoulomb

  !****************************************************************************
  !****s* checks/checkPot
  ! NAME
  ! subroutine CheckPot
  ! PURPOSE
  !
  !****************************************************************************
!   subroutine CheckPot
!     use particleDefinition
!     use IdTable, only: nucleon, Delta, pion, phi
!     use yukawa, only: getYukawaAlpha, yukPot
!     use densitymodule, only: densityAt
!     use dichtedefinition
!     use potentialModule, only: potential_LRF
!     use constants, only: mN, rhoNull
!     use energyCalc, only: energyDetermination
!
!     integer i,j
!     real ::  x,y,p,dens
!     real,dimension(1:3) :: r,place
!     logical,parameter :: checkEnergy=.false.
!
!     type(particle) :: testParticle
!     type(dichte) :: density
!
!     Print *, "Yukawa-alpha:",getYukawaAlpha()
!
!     Write(*,*) 'Check yukawa'
!     Open(11,File='CheckYukawa.dat')
!     Do i=-25,25
!        x=i*0.5
!        Do j=-25,25
!           y=j*0.5
!           r=(/x,0.,y/)
!           write(11,fmt='(3f6.3)') x,y,yukPot(r)
!        end do
!     end do
!     Close(11)
!
!
!     testparticle%Id=nucleon
!     Open(11,File='./output/nucleonPot.dat')
!     Do i=0,6
!        dens=i*0.25*rhoNull
!        Do j=0,10
!           p=j*0.04
!           testparticle%position(1:3)=(/x,0.,0./)
!           testparticle%momentum(1:3)=(/p,0.,0./)
!           write(11,fmt='(3f6.3)') dens/rhoNull,p,potential_LRF(testParticle)
!        end do
!     end do
!     close(11)
!
!     testparticle%Id=phi
!     Open(11,File='./output/phiPot.dat')
!     Do i=0,25
!        x=i*0.5
!        Do j=0,25
!           p=j*0.04
!           testparticle%position(1:3)=(/x,0.,0./)
!           testparticle%momentum(1:3)=(/p,0.,0./)
!           place=testparticle%position(1:3)
!           density=densityAt(place)
!           write(11,fmt='(3f6.3)') density%baryon(0)/rhoNull,p,potential_LRF(testParticle)
!        end do
!     end do
!     Close(11)
!
!     testparticle%Id=pion
!     Open(11,File='./output/pionPot.dat')
!     Do i=0,6
!        dens=i*0.25*rhoNull
!        Do j=0,10
!           p=j*0.04
!           testparticle%position(1:3)=(/x,0.,0./)
!           testparticle%momentum(1:3)=(/p,0.,0./)
!           write(11,fmt='(3f6.3)') dens/rhoNull,p,potential_LRF(testParticle)
!        end do
!     end do
!     close(11)
!
!     testparticle%Id=delta
!     Open(11,File='./output/deltaPot.dat')
!     Do i=0,25
!        x=i*0.5
!        Do j=0,25
!           p=j*0.04
!           testparticle%position(1:3)=(/x,0.,0./)
!           testparticle%momentum(1:3)=(/p,0.,0./)
!           place=testparticle%position(1:3)
!           density=densityAt(place)
!           write(11,fmt='(3f6.3)') density%baryon(0),p,potential_LRF(testParticle)
!        end do
!     end do
!     Close(11)
!
!     If (checkEnergy) then
!        Write(*,*) 'Check EnergyDetermination'
!        testparticle%Id=nucleon
!        testparticle%charge=1
!        testparticle%mass=mN
!        testparticle%position(1:3)=(/0.,0.,0./)
!        p=0.
!        testparticle%momentum(0:3)=(/SQRT(mN**2+p**2),p,0.,0./)
!        place=testparticle%position(1:3)
!        density=densityAt(place)
!        write(*,*) density%baryon(0),testparticle%momentum(0:3)
!        call energyDetermination(testparticle)
!        write(*,*) density%baryon(0),testparticle%momentum(0:3)
!        Write(*,*) potential_LRF(testParticle)+mN
!        Write(*,*)
!
!        testparticle%position(1:3)=(/0.,0.,0./)
!        p=0.1
!        testparticle%momentum(0:3)=(/SQRT(mN**2+p**2),p,0.,0./)
!        place=testparticle%position(1:3)
!        density=densityAt(place)
!        write(*,*) density%baryon(0),testparticle%momentum(0:3)
!        call energyDetermination(testparticle)
!        write(*,*) density%baryon(0),testparticle%momentum(0:3)
!        Write(*,*) potential_LRF(testParticle,density%baryon(0))+SQRT(mN**2+p**2)
!        Write(*,*)
!
!        testparticle%position(1:3)=(/0.,0.,0./)
!        p=10.
!        testparticle%momentum(0:3)=(/SQRT(mN**2+p**2),p,0.,0./)
!        place=testparticle%position(1:3)
!        density=densityAt(place)
!        write(*,*) density%baryon(0),testparticle%momentum(0:3)
!        call energyDetermination(testparticle)
!        write(*,*) density%baryon(0),testparticle%momentum(0:3)
!        Write(*,*) potential_LRF(testParticle,density%baryon(0))+SQRT(mN**2+p**2)
!        Write(*,*)
!
!        testparticle%position(1:3)=(/1120.,0.,0./)
!        p=10.
!        testparticle%momentum(0:3)=(/SQRT(mN**2+p**2),p,0.,0./)
!        place=testparticle%position(1:3)
!        density=densityAt(place)
!        write(*,*) density%baryon(0),testparticle%momentum(0:3)
!        call energyDetermination(testparticle)
!        write(*,*) density%baryon(0),testparticle%momentum(0:3)
!        Write(*,*) potential_LRF(testParticle)+SQRT(mN**2+p**2)
!        Write(*,*)
!     end if
!   end subroutine CheckPot


  !****************************************************************************
  !****s* checks/checkRadius
  ! NAME
  ! subroutine checkRadius(teilchen, time)
  ! PURPOSE
  ! Average radius over all particle vector
  !****************************************************************************
  subroutine checkRadius(teilchen, time)

    use particleDefinition
    use IDTable

    type(particle), dimension(:,:),intent(in) :: teilchen
    real, intent(in), optional :: time
    integer :: i,j
    real :: Radius, Radius2, r
    integer :: numberNucleons
    logical,save :: firstTime=.true.

    if (firsttime) then
       open(20,file='MeanRadius.dat')
       firsttime=.false.
    else
       open(20,file='MeanRadius.dat',position='Append')
    end if

    numberNucleons=0
    Radius=0.
    Radius2=0.
    do j=1,size(teilchen,dim=2)
       do i=1,size(teilchen,dim=1)
          if (teilchen(i,j)%Id.eq.nucleon) then
             r = absPos(teilchen(i,j))
             Radius  = Radius  + r
             Radius2 = Radius2 + r**2
             numberNucleons=numberNucleons+1
          end if
       end do
    end do
    Radius  = Radius/float(numberNucleons)
    Radius2 = Radius2/float(numberNucleons)
    if (sqrt(Radius2).lt.minimalRadius) minimalradius=sqrt(Radius2)
    if (sqrt(Radius2).gt.maximalRadius) maximalradius=sqrt(Radius2)

    if (present(time)) then
       write(20,fmt='(3F14.5,I12,2F14.5)') time,Radius,sqrt(Radius2),&
            numberNucleons,minimalRadius,maximalRadius
    else
       write(20,fmt='(2F14.5,I6)') Radius,sqrt(Radius2),numberNucleons
    end if
    close(20)
  end subroutine checkRadius

  !****************************************************************************
  !****f* checks/checkRadius_getAmplitude
  ! NAME
  ! subroutine checkRadius_getAmplitude()
  ! PURPOSE
  ! return the value (MaximalRadius-MinimalRadius)
  !****************************************************************************
!   real function checkRadius_getAmplitude()
!     checkRadius_getAmplitude=abs(MaximalRadius-MinimalRadius)
!   end function checkRadius_getAmplitude

  !****************************************************************************
  !****s* checks/checkRadius_resetMinMax
  ! NAME
  ! subroutine checkRadius_resetMinMax()
  ! PURPOSE
  ! reset the values of MaximalRadius and MinimalRadius
  !****************************************************************************
!   subroutine checkRadius_resetMinMax()
!     MinimalRadius=100000.
!     MaximalRadius=-100000.
!   end subroutine checkRadius_resetMinMax


  !****************************************************************************
  !****f* checks/saveEpot
  ! NAME
  ! real function saveEpot(E,rho,i)
  ! PURPOSE
  ! Sums the potential Enrgy of all nucleons
  ! INPUTS
  ! * real :: E -- energy
  ! * real :: rho -- density
  ! * logical :: i -- .true.: store; .false.: retrieve and delete
  ! OUTPUT
  ! Energy in units of MeV
  !****************************************************************************
!   real function saveEpot(E,rho,i)
!     real, intent(in)::E,rho
!     logical, intent(in)::i
!     real, save::Energy=0.,sumdensity=0.
!
!     if (i) then
!        Energy=Energy+E
!        sumdensity=sumdensity+rho
!        saveEpot=0.
!     else
!        saveEpot=Energy/sumdensity
!        Energy=0.
!        sumdensity=0.
!     end if
!
!   end function saveEpot

  !****************************************************************************
  !****s* checks/checkEnergyLDA
  ! NAME
  ! subroutine checkEnergyLDA(teilchen)
  ! PURPOSE
  !
  !****************************************************************************
  subroutine checkEnergyLDA(teilchen)
    use particleDefinition
    use IDTable
    use nucleusDefinition
    use nucDLDA, only: getEParticleLaplace
    use nucleus, only: getTarget
    use dichteDefinition
    use densitymodule, only: densityAt,gridSpacing
    use baryonPotentialModule, only: rhoLaplace

    type(particle), dimension(:,:),intent(in) :: teilchen
    real,dimension(1:3)  :: rvec
    real, dimension(1:3), save :: stepsize=0.
    real, dimension(1:4) :: Etemp, EnergySum
    real::rhocent,rhoptemp,ptemp,rabstemp
    integer,save :: schritt=0
    integer :: i,j,l,numberParticles,notBound,outside,outside1,outside2
    logical,save :: firstTime=.true.
    type(dichte) :: pos
    type(tnucleus), pointer, save :: targetNucleus
    real, save :: kernRadius=0.
    EnergySum=0.
    numberParticles=0
    notBound=0
    outside=0
    outside1=0
    outside2=0
    schritt=schritt+1
    if (firsttime) then
       open(20,file='MeanEnergy.dat')
       targetNucleus => getTarget()
       kernRadius=targetNucleus%MaxDist
       stepsize=gridSpacing
       firsttime=.false.
    else
       open(20,file='MeanEnergy.dat',position='Append')
    end if


    do j=1,size(teilchen,dim=2)
       do i=1,size(teilchen,dim=1)
          if (teilchen(i,j)%Id > 0) then
             numberParticles=numberParticles+1
             rvec=teilchen(i,j)%position(1:3)
             pos=densityAt(rvec)
             rhocent=SQRT(pos%baryon(0)**2 &
                  -Dot_Product(pos%baryon(1:3),pos%baryon(1:3)))
             rhoptemp=0.
             ptemp=0.
             rabstemp=0.
             do l=1,3,1
                ptemp=ptemp+teilchen(i,j)%momentum(l)**2
                rabstemp=rabstemp+rvec(l)**2
             end do
             rabstemp=SQRT(rabstemp)
             rhoptemp=rhoLaplace(rvec,stepsize)
             call getEParticleLaplace(Etemp,rhocent,rhoptemp,ptemp)
             EnergySum=EnergySum+Etemp
             if (Etemp(1).ge.0.) then
                notBound=notBound+1
             end if
             if (rabstemp.gt.kernRadius) then
                outside=outside+1
                if (rabstemp.gt.(kernRadius+1.)) then
                   outside1=outside+1
                   if (rabstemp.gt.(kernRadius+2.)) then
                      outside2=outside+1
                   end if
                end if
             end if
          end if
       end do
    end do


    EnergySum=EnergySum/float(numberParticles)


    ! Ausgabe von Zeitschritt, Gesamtenergie, kinetischer Energie,
    ! Potentialenergie, Gradientenenergie, Anzahl der Teilchen,
    ! Anzahl der ungebundenen Teilchen,
    ! Teilchen die über den maximalen Radius der berechneten Dichteverteilung
    ! hinaus sind, dasselbe +1fm, +2fm
    write(20,fmt='(I5,4(1X,F10.4),5(1X,I7))') &
          schritt,EnergySum(1),EnergySum(2),EnergySum(3),EnergySum(4),&
          numberParticles,notBound,outside,outside1,outside2
    close(20)

  end subroutine checkEnergyLDA

  !****************************************************************************
  !****s* checks/checkEnergyLDAWelke
  ! NAME
  ! subroutine checkEnergyLDAWelke(teilchen)
  ! PURPOSE
  !
  !****************************************************************************
  subroutine checkEnergyLDAWelke(teilchen)
    use particleDefinition
    use IDTable
    use nucleusDefinition
    use nucDLDA, only: getEParticleLaplaceWelke
    use nucleus, only: getTarget
    use dichteDefinition
    use densitymodule, only: densityAt,gridSpacing
    use baryonPotentialModule, only: rhoLaplace

    type(particle), dimension(:,:),intent(in) :: teilchen
    real,dimension(1:3)  :: rvec
    real, dimension(1:3), save :: stepsize=0.
    real, dimension(1:5) :: Etemp, EnergySum
    real::rhocent,rhoptemp,ptemp,rabstemp,Energymomentum
    integer,save :: schritt=0
    integer :: i,j,l,numberParticles,notBound,outside,outside1,outside2
    logical,save :: firstTime=.true.
    type(dichte) :: pos
    type(tnucleus), pointer, save :: targetNucleus
    real, save :: kernRadius=0.
    EnergySum=0.
    numberParticles=0
    notBound=0
    outside=0
    outside1=0
    outside2=0
    schritt=schritt+1
    if (firsttime) then
       open(20,file='MeanEnergy.dat')
       targetNucleus => getTarget()
       kernRadius=targetNucleus%MaxDist
       stepsize=gridSpacing
       firsttime=.false.
    else
       open(20,file='MeanEnergy.dat',position='Append')
    end if

    Energymomentum=0.

    do j=1,size(teilchen,dim=2)
       do i=1,size(teilchen,dim=1)
          if (teilchen(i,j)%Id > 0) then
             Energymomentum=Energymomentum+teilchen(i,j)%momentum(0)
             numberParticles=numberParticles+1
             rvec=teilchen(i,j)%position(1:3)
             pos=densityAt(rvec)
             rhocent=SQRT(pos%baryon(0)**2 &
                  -Dot_Product(pos%baryon(1:3),pos%baryon(1:3)))
             rhoptemp=0.
             ptemp=0.
             rabstemp=0.
             do l=1,3,1
                ptemp=ptemp+teilchen(i,j)%momentum(l)**2
                rabstemp=rabstemp+rvec(l)**2
             end do
             rabstemp=SQRT(rabstemp)
             rhoptemp=rhoLaplace(rvec,stepsize)
             call getEParticleLaplaceWelke(Etemp,rhocent,rhoptemp,ptemp)
             EnergySum=EnergySum+Etemp
             if (Etemp(1).ge.0.) then
                notBound=notBound+1
             end if
             if (rabstemp.gt.kernRadius) then
                outside=outside+1
                if (rabstemp.gt.(kernRadius+1.)) then
                   outside1=outside+1
                   if (rabstemp.gt.(kernRadius+2.)) then
                      outside2=outside+1
                   end if
                end if
             end if
          end if
       end do
    end do


    EnergySum=EnergySum/float(numberParticles)
    Energymomentum=Energymomentum/float(numberParticles)*1000.-938.

    ! Ausgabe von Zeitschritt, Gesamtenergie, kinetischer Energie,
    ! Potentialenergie, Gradientenenergie, Anzahl der Teilchen,
    ! Anzahl der ungebundenen Teilchen, Teilchen die über den maximalen Radius
    ! der berechneten Dichteverteilung hinaus sind,
    ! dasselbe +1fm, +2fm, Energie durch Impulsdifferenz
    write(20,fmt='(I5,4(1X,F10.4),5(1X,I7),2(1X,F10.4))') &
         schritt,EnergySum(1),EnergySum(2),EnergySum(3),EnergySum(4), &
         numberParticles,notBound,outside,&
         outside1,outside2,EnergySum(5),Energymomentum
    close(20)

  end subroutine checkEnergyLDAWelke

  !****************************************************************************
  !****s* checks/evaluateEnergy
  ! NAME
  ! subroutine evaluateEnergy(teilchen,time)
  ! PURPOSE
  ! print the averaged kinetic energy sqrt(p^2+m^2) and also p0 per time step
  !****************************************************************************
  subroutine evaluateEnergy(teilchen,time)
    use particleDefinition
    use idTable, only: nucleon

    type(particle), dimension(:,:),intent(in) :: teilchen
    real, intent(in) :: time
    integer :: i,j
    real :: eKin, EKin0
    integer :: numberParticles
    logical,save :: firstTime=.true.

    if (firstTime) then
       open(20,file='Energy.dat')
       firstTime=.false.
    else
       open(20,file='Energy.dat',position='Append')
    end if

    Ekin  = 0.
    Ekin0 = 0.
    numberParticles=0
    ! Sum the kinetic energy over all ensembles
    do i=1,size(teilchen,dim=1)
       do j=1,size(teilchen,dim=2)
          if (teilchen(i,j)%ID.eq.nucleon) then
             Ekin=Ekin+KineticEnergy(teilchen(i,j))
             Ekin0=Ekin0+teilchen(i,j)%momentum(0)-teilchen(i,j)%mass
             numberParticles=numberParticles+1
          end if
       end do
    end do
    ! Divide the kinetic energy by the number of particles
    Ekin=Ekin/real(numberParticles)
    Ekin0=Ekin0/real(numberParticles)
    ! Printout the kinetic Energy
    write(20,fmt='(f9.4,3F13.6)') time,Ekin,Ekin0
    close(20)
  end subroutine evaluateEnergy


  !****************************************************************************
  !****s* checks/evaluateBindingEnergy_Teis
  ! NAME
  ! subroutine evaluateBindingEnergy_Teis(teilchen, time)
  ! PURPOSE
  ! Print binding energy of all nucleons (no Coulomb effects included).
  !****************************************************************************
  subroutine evaluateBindingEnergy_Teis(teilchen, time)
    use particleDefinition
    use potentialModule, only: trueEnergy
    use IdTable, only: isBaryon, nucleon
    use output, only: writeFileDocu
    use histf90
    use constants, only: mN

    type(particle), dimension(:,:), intent(in) :: teilchen
    real, intent(in) :: time

    real :: EE
    real, dimension(3) :: E_all, E_bar, E_nuc, E_nucB, EE_nucB, dE
    integer :: i, j, nEns, nPart
    type(histogram) :: hBA

    logical,save :: firstTime=.true.

    if (firstTime) then
       call writeFileDocu('BindingEnergy.dat', &
            'Includes binding energy of all nucleons, no Coulomb effects included')
       open(22,file='BindingEnergy.dat')
       write(22,*) '# 1: time'
       write(22,*) '# 2: total energy / number of all particles'
       write(22,*) '# 3: total baryon energy / net number of baryons'
       write(22,*) '# 4: total nucleon energy / number of nucleons'
       write(22,*) '# 5: total bound nucleon energy / number of bound nucleons'
       write(22,*) '# 6: number of all particles per ensemble'
       write(22,*) '# 7: net number of baryons per ensemble'
       write(22,*) '# 8: number of nucleons per ensemble'
       write(22,*) '# 9: number of bound nucleons per ensemble'

       call writeFileDocu('BindingEnergyArr.dat', &
            'The distribution connected with the average in "BindingEnergy.dat"')
       open(21,file='BindingEnergyArr.dat')

       firstTime=.false.
    else
       open(22,file='BindingEnergy.dat',position='Append')
       open(21,file='BindingEnergyArr.dat',position='Append')
       write(21,*)
       write(21,*)
    end if

    call CreateHist(hBA, "B/A for bound nucleons", -0.1,0.1,0.001)

    E_all=0.
    E_bar=0.
    E_nuc=0.
    E_nucB=0.

    nEns=size(teilchen,dim=1)
    nPart=size(teilchen,dim=2)

    ! Sum the full energy over all ensembles:
    do i=1,nEns
       EE_nucB=0.
       do j=1,nPart
          ! Sum up energies of all particles (not only nucleons):
          if (teilchen(i,j)%ID > 0) then
             EE = trueEnergy(teilchen(i,j))
!             EE = teilchen(i,j)%momentum(0)

!             dE = teilchen(i,j)%perWeight*(/1.,EE,0./)+(/0.,0.,1./)
             dE = (/1.,EE,1./)
             E_all=E_all+dE
             ! Count only baryons, since baryon number must be conserved:
             if (isBaryon(teilchen(i,j)%ID)) then
                if (.not.teilchen(i,j)%antiparticle) then
                   E_bar=E_bar+dE
                   if (teilchen(i,j)%ID==nucleon) then
                      E_nuc=E_nuc+dE
                      if (teilchen(i,j)%momentum(0)<=mN) then
                         E_nucB=E_nucB+dE
                         EE_nucB=EE_nucB+dE
                      end if
                   end if
                else
                   E_bar=E_bar+dE
                end if
             end if
          else if (teilchen(i,j)%ID < 0) then
             exit
          end if
       end do

       if (EE_nucB(3)>0.) &
            call AddHist(hBA,EE_nucB(2)/EE_nucB(1)-mN,EE_nucB(1)/EE_nucB(3))

    end do

    ! Printout the Energy per baryon:
    write(22,fmt='(F8.2,4ES14.6,4F10.3)') time, &
         E_all(2)/E_all(1), E_bar(2)/E_bar(1), &
         E_nuc(2)/E_nuc(1), E_nucB(2)/E_nucB(1), &
         E_all(3)/nEns,     E_bar(3)/nEns, &
         E_nuc(3)/nEns,     E_nucB(3)/nEns
    close(22)

    write(21,*) '# t=',time
    call WriteHist(hBA,21,mul=1./nEns)
    close(21)

  end subroutine evaluateBindingEnergy_Teis


  !****************************************************************************
  !****s* checks/checkTrajectories
  ! NAME
  ! subroutine checkTrajectories(teilchen)
  ! PURPOSE
  !
  !****************************************************************************
!   subroutine checkTrajectories(teilchen)
!     use particleDefinition
!     use IDTable
!
!     logical,save :: initFlag=.true.
!     type(particle), dimension(:,:),intent(in) :: teilchen
!     If (initFlag) then
!        Open(20,file='Trajectory_1.dat')
!        Open(21,file='Trajectory_2.dat')
!        Open(22,file='Trajectory_3.dat')
!        Open(23,file='Trajectory_3.dat')
!        initFlag=.false.
!     else
!        Open(20,file='Trajectory_1.dat',position='Append')
!        Open(21,file='Trajectory_2.dat',position='Append')
!        Open(22,file='Trajectory_3.dat',position='Append')
!        Open(23,file='Trajectory_3.dat',position='Append')
!     end if
!
!     Write(20,fmt='(3f8.2)') teilchen(1,1)%position(1:3)
!     If(size(teilchen,dim=2).ge.4) then
!        Write(21,fmt='(3f8.2)') teilchen(1,2)%position(1:3)
!        Write(22,fmt='(3f8.2)') teilchen(1,3)%position(1:3)
!        Write(23,fmt='(3f8.2)') teilchen(1,4)%position(1:3)
!     end if
!     close(20)
!     close(21)
!     close(22)
!     close(23)
!
!   end subroutine checkTrajectories


  !****************************************************************************
  !****s* checks/evaluateTimeStep
  ! subroutine evaluateTimeStep(iflag,coll_num,delta_T_max,time,delta_T_new)
  ! PURPOSE
  ! Compute time step using the frequency of two-body collisions
  ! INPUTS
  ! * integer :: iflag  -- 1: use real-real counted collisions;
  !   2: use real-perturbative counted collisions
  ! * real :: coll_num  -- allow .leq.coll_num collisions/particle
  ! * real :: delta_T_max -- maximal allowed time step (fm/c)
  ! * real :: time -- current time (fm/c)
  ! OUTPUT
  ! * real :: delta_T_new -- calculated time step
  ! NOTES
  ! If gridSpacing(3) <  delta_T_max, the upper limit on the time step
  ! is chosen as gridSpacing(3). This is important for high energy HIC
  ! due to Lorentz contraction along z-axis.
  !****************************************************************************
  subroutine evaluateTimeStep(iflag,coll_num,delta_T_max,time,delta_T_new)

    use collisionNumbering, only: GetCountedEvents,getn_participants
    use densitymodule, only: get_realGridSpacing

    real, intent(in) :: coll_num
    integer, intent(in) :: iflag
    real, intent(in) :: delta_T_max
    real, intent(in) :: time
    real, intent(out) :: delta_T_new

    real, parameter ::  kappa=1.   ! parameter which determines the slope of
                                   ! increasing time step when collision
                                   ! frequency drops

    integer, dimension(1:2) :: n_participants
    real, dimension(1:3) :: GridSpacing
    integer :: num_2body,num_part!,i
    real :: frequency,freq,slope,delta_T_max1
    real, save :: delta_T, freq_old
    real, save :: time_old=1000.
    logical, save :: flag

    n_participants=getn_participants()
    GridSpacing = get_realGridSpacing()

    delta_T_max1 = min(delta_T_max,GridSpacing(3))  ! Time step restriction due to grid spacing

    if ( time < time_old ) then
       delta_T = delta_T_max1
       freq_old=0.
       flag = .true.
    end if

    if (iflag.eq.1) then
       num_2body=2*getCountedEvents(1,2,1)
       num_part=n_participants(1)
    else if (iflag.eq.2) then
       num_2body=getCountedEvents(1,2,2)
       num_part=n_participants(2)
    else
       write(*,*) ' In checkTimeStep: wrong iflag', iflag
       stop
    end if

    if (delta_T.le.0.) then
       write(*,*) ' In checkTimeStep: delta_T=', delta_T
       stop
    end if

    if (num_part.gt.0) then
       frequency=float(num_2body)/float(num_part)/delta_T
    else
       frequency=0.
    end if

    freq=float(num_2body)/delta_T

    if ( flag .and. num_2body > 0 ) flag = .false.

    if ( flag ) then

       delta_T_new = delta_T

    else

       if (freq.ge.freq_old) then
          if (frequency.gt.0.) then
             delta_T_new=min(delta_T_max1,coll_num/frequency)
          else
             delta_T_new=delta_T_max1
          end if
       else
          if (freq.gt.0.) then
             slope=1.+(freq_old/freq-1.)*kappa
             delta_T_new=min(delta_T_max1,slope*delta_T)
          else
             delta_T_new=delta_T_max1
          end if
       end if

    end if

    delta_T=delta_T_new
    freq_old=freq
    time_old = time

    open(1,file='timestep.dat',status='unknown',position='append')
    write(1,'(f6.3,1x,I6,1x,I6,1x,e13.6,1x,f6.3)') time,num_2body,num_part,&
         freq,delta_T_new
    close(1)

  end subroutine evaluateTimeStep


  !****************************************************************************
  !****s* checks/evaluateTotal4Momentum_RMF
  ! subroutine evaluateTotal4Momentum_RMF(teilchen,time)
  ! PURPOSE
  ! Compute the total energy and momentum of the system within
  ! the relativistic mean field model.
  ! INPUTS
  ! * type(particle), dimension(:,:), intent(in) :: teilchen  ! particle array
  ! * real, intent(in) :: time ! current time (fm/c)
  ! OUTPUT
  ! Prints out the total 4-momentum vs time.
  !****************************************************************************
  subroutine evaluateTotal4Momentum_RMF(teilchen,time)

    use particleDefinition, only: particle
    use IdTable, only: isBaryon
    use densitymodule, only: true4Momentum_RMF
    use coulomb, only: emfoca

    type(particle), dimension(:,:), intent(in) :: teilchen
    real, intent(in) :: time

    integer :: i,j
    integer :: nEns, nPart
    real :: baryonNumber, charge, antibaryonNumber, baryons_inside_grid, meff, meff_aver
    real, dimension(0:3) :: Pstar_tot, P_tot ! Total kinetic and canonical 4-momenta of particles.
    real, dimension(0:3) :: momentum

    real :: P_tc
    real, dimension(1:3) :: place,impuls

    logical :: flag

    nEns = size(teilchen,dim=1)
    nPart = size(teilchen,dim=2)

    open(20,file='Total4Momentum_RMF.dat',position='Append')

    Pstar_tot = 0.
    P_tot = 0.
    P_tc = 0.
    baryonNumber = 0.
    charge = 0.
    antibaryonNumber = 0.
    baryons_inside_grid = 0.
    meff_aver = 0.

    ! Sum single-particle kinetic 4-momenta:

    Loop_over_ensembles : do i=1,nEns
       Loop_over_particles : do j=1,nPart

          if ( teilchen(i,j)%id == 0 ) then
             cycle Loop_over_particles
          else if ( teilchen(i,j)%id < 0 ) then
             exit Loop_over_particles
          end if

          Pstar_tot(0:3) = Pstar_tot(0:3) + teilchen(i,j)%momentum(0:3)

          call true4Momentum_RMF(teilchen(i,j),momentum,flag)

          charge = charge + float(teilchen(i,j)%charge)

          if ( isBaryon(teilchen(i,j)%ID) ) then
             meff = sqrt( teilchen(i,j)%momentum(0)**2 - dot_product( teilchen(i,j)%momentum(1:3), &
                  &teilchen(i,j)%momentum(1:3) ) )
             meff_aver = meff_aver + meff
             if ( .not.teilchen(i,j)%antiparticle ) then
                baryonNumber = baryonNumber + 1.
             else
                antibaryonNumber = antibaryonNumber + 1.
             end if
             if (flag) baryons_inside_grid = baryons_inside_grid + 1.
          end if

          P_tot(0:3) = P_tot(0:3) + momentum(0:3)

          ! coulomb contribution
          ! emfoca gives the electrostatic scalar potential from the solution
          ! of the Poisson-equation. By considering only the 0-component of
          ! the vector potential A^{\mu}=(\phi,\vec{A}) in T^{00}(x), one
          ! obtains V_c(x) = 1/2*e*j^{0}_{p}*A^{0} as the coulomb contribution
          ! to the energy density T^{00}.
          ! The total coulomb energy is derived from the space integration over
          ! V_c(x) leading to a simple sumation over
          ! the test particles (j(x) is represented in terms of delta functions
          ! within the test-particle formalism).
          place(1:3)  = teilchen(i,j)%position(1:3)
          impuls(1:3) = teilchen(i,j)%momentum(1:3)
          P_tc = P_tc + momentum(0) &
               + 0.5*emfoca(place,impuls,teilchen(i,j)%charge,teilchen(i,j)%ID)

       end do Loop_over_particles
    end do Loop_over_ensembles

    Pstar_tot = Pstar_tot / float(nEns)
    P_tot = P_tot / float(nEns)
    P_tc = P_tc / float(nEns)
    baryonNumber = baryonNumber / float(nEns)
    charge = charge / float(nEns)
    antibaryonNumber = antibaryonNumber / float(nEns)
    baryons_inside_grid = baryons_inside_grid / float(nEns)
    meff_aver = meff_aver / float(nEns) / baryonNumber

    write(20,'(11(1x,f13.6))')   time, &
         & baryonNumber, charge, antibaryonNumber, baryons_inside_grid, &
         & P_tc/(baryonNumber-antibaryonNumber),&
         & P_tot(0:3)/(baryonNumber-antibaryonNumber),meff_aver


    close(20)

  end subroutine evaluateTotal4Momentum_RMF


  !****************************************************************************
  !****s* checks/CheckGridSize
  ! subroutine CheckGridSize(teilchen,time,time_max,TheEventType,gridSpacing,
  ! gridSize)
  ! PURPOSE
  ! Checks if real particles start to escape out of grid.
  ! INPUTS
  ! * type(particle), dimension(:,:) :: teilchen        -- particle array
  ! * real :: time            -- current time (fm/c)
  ! * real :: time_max -- max. time of simulation (fm/c)
  ! * integer :: TheEventType    -- Eventclass
  ! * real :: gridSpacing     -- size of cells
  ! * real :: gridSize        -- Size of the box
  !****************************************************************************
  subroutine CheckGridSize(teilchen,time,time_max,TheEventType,gridSpacing, &
    gridSize)

    use particleDefinition
    use IdTable, only: isMeson
    use nucleusDefinition
    use nucleus, only: getTarget, getProjectile

    !--------------------------------------------------------------------------
    ! Input-Output Variables
    !--------------------------------------------------------------------------
    type(particle), dimension(:,:), intent(in) :: teilchen
    real,           dimension(1:3), intent(in) :: gridSpacing,gridSize
    real,   intent(in)   :: time,time_max
    integer,intent(in)   :: TheEventType
    !--------------------------------------------------------------------------
    !Local variables
    !--------------------------------------------------------------------------
    real, save,allocatable,dimension(:) :: OutOfGrid
    real,    save  :: SystemSize
    integer, save  :: affected_events
    logical, save  :: initGridFlag=.true.
    logical, save  :: WarningFlag=.false.
    real, dimension(1:3) :: PosOrig
    real           :: MeanValue
    integer        :: i,j,nEns,nPart
    type(tnucleus),pointer :: TheTarget,TheProjectile

    if (initGridFlag) then
       TheTarget     => getTarget()
       TheProjectile => getProjectile()
       SystemSize    = float(TheTarget%Mass)
       if (TheEventType==1) &
            SystemSize = float(TheTarget%Mass + TheProjectile%Mass)
       if (TheEventType==300) &
            SystemSize = float(TheTarget%Mass+1)
       allocate(OutOfGrid(1:Size(teilchen,dim=1)))
       OutOfGrid(:)    = 0.
       affected_events = 0
       initGridFlag    = .false.
    end if


    if (.not.WarningFlag) then !do the check ONLY at once per subsequent run
       nEns = size(teilchen,dim=1)
       nPart = size(teilchen,dim=2)

       do i=1,nEns
          do j=1,nPart
             if ( teilchen(i,j)%id == 0 ) then
                cycle
             else if ( teilchen(i,j)%id < 0 ) then
                exit
             end if
             if (isMeson(teilchen(i,j)%ID)) cycle ! only baryons
             posOrig(:)=abs( NINT(Teilchen(i,j)%position(:)/gridSpacing(:)) )
             if ( (posOrig(1) > gridSize(1)) .or. &
                  & (posOrig(2) > gridSize(2)) .or. &
                  & (posOrig(3) > gridSize(3)) )  &
                  OutOfGrid(i) = OutOfGrid(i) + 1.
          end do
          !more than 10% of baryons outside of grid=>event affected
          if ( (OutOfGrid(i)/SystemSize)*100. > 10. ) then
             affected_events = affected_events + 1
          end if
       end do


       if ( affected_events == nEns ) then
          MeanValue = Sum(OutOfGrid)/float(nEns)
          write(*,*)
          write(*,*) '   Warning in CheckGridSize, module checks.f90'
          write(*,*) '   more than 10% of particles/ensemple out of grid!'
          write(*,'(A,f7.2,A,f5.2,A,f5.2,A)') '   <Number of escaped baryons> = ', &
               & MeanValue,' corresponding to ',(MeanValue/SystemSize)*100., &
               & '% at time => ',time,' fm/c'
          write(*,*)
          WarningFlag     = .true.
          !initialize again OutOgGrid(:) and affected_events after the check:
          OutOfGrid(:)    = 0.
          affected_events = 0
       end if
    end if


    !--------------------------------------------------------------------------
    ! reset to init-values warningFlag,OutOfGrid(:) and affected_events
    ! for the next subsequent run (just to be sure...)
    !--------------------------------------------------------------------------
    if (abs(time-time_max) < 0.0001) then
       WarningFlag     = .false.
       OutOfGrid(:)    = 0.
       affected_events = 0
    end if

  end subroutine CheckGridSize

  !****************************************************************************
  !****is* checks/ChecksCallTachyon
  ! NAME
  ! subroutine ChecksCallTachyon(Parts,Text)
  ! PURPOSE
  ! This routine checks, whether particles with m_eff^2 = abs4Sq(momentum) < 0
  ! exist. These are called 'tachyons' in this context.
  !
  ! If the parameter 'TachyonIsBlocking' is set, this routines stops when
  ! tachyons occured.
  !****************************************************************************
  subroutine ChecksCallTachyon(Parts,Text)
    use particleDefinition
    use minkowski, only: abs4Sq
    use output, only: WriteParticle
    use callstack, only: traceback
!!$    use neutrinoProdInfo, only: neutrinoProdInfo_Dump

    type(particle),intent(in),dimension(:,:) :: Parts
    character(*), intent(in) :: Text

    integer :: i,j
    logical :: found = .false.

    do i=lbound(Parts,dim=1),ubound(Parts,dim=1)
       do j=lbound(Parts,dim=2),ubound(Parts,dim=2)
          if (Parts(i,j)%ID<=0) cycle

          if (abs4Sq(Parts(i,j)%momentum) < -1E-6) then
             found = .true.
             write(*,'(A,A)') Text,'   Tachyon found !!!'
             write(*,*) '  M^2 = ', abs4Sq(Parts(i,j)%momentum)
             call WriteParticle(6,i,j,Parts(i,j))

             ! since the problem with the tachyons occured within the neutrino
             ! business, I added the following line for debugging purposes.
             ! In principle, you can remove this line (and the 'use ...' above)
!!$             call neutrinoProdInfo_Dump(Parts(i,j)%firstevent)

          end if

       end do
    end do

    if ((found).and.(TachyonIsBlocking)) call TRACEBACK()

  end subroutine ChecksCallTachyon

  !****************************************************************************
  !****is* checks/ChecksCallPertFlag
  ! NAME
  ! subroutine ChecksCallPertFlag(particles, flag)
  ! PURPOSE
  ! This routine checks, whether the particle property '%perturbative'
  ! is set to the given value (.true./.false.( for all particles
  !****************************************************************************
  subroutine ChecksCallPertFlag(particles, flag)
    use particleDefinition
    use callstack, only: traceback
    use output, only: WriteParticle

    type(particle), dimension(:,:), intent(in) :: particles
    logical, intent(in) :: flag

    integer :: iEns, iPart

    do iEns=lbound(particles,dim=1),ubound(particles,dim=1)
       do iPart=lbound(particles,dim=2),ubound(particles,dim=2)
          if (particles(iEns,iPart)%ID <= 0) cycle
          if (particles(iEns,iPart)%perturbative .neqv. flag) then
             call WriteParticle(6,iEns,iPart,particles(iEns,iPart))
             write(*,*) 'checking for flag=',flag
             call traceback('flag %perturbative has wrong value!')
          end if
       end do
    end do


  end subroutine ChecksCallPertFlag

  !****************************************************************************
  !****is* checks/CheckConservation
  ! NAME
  ! subroutine CheckConservation(particles)
  ! PURPOSE
  ! This routine checks the conservation of energy/momentum, baryon number
  ! and strangeness between time steps.
  ! NOTES
  ! * This routine only does its job when called with 'realParticles'
  ! * Decays including photons violate energy/momentum conservation.
  !   You may switch these off by setting
  !     FileNameDecayChannels = "DecayChannels.noPhotonicMeson.dat"
  !   in the namelist 'initDatabase'
  !****************************************************************************
  subroutine CheckConservation(particles)
    use particleDefinition
    use callstack, only: traceback
    use IdTable, only: isBaryon, isMeson
    use particleProperties, only: hadron

    type(particle), dimension(:,:), intent(in), target :: particles

    integer :: iEns, iPart
    type(particle), pointer :: pPart
    real, save, dimension(0:3) :: MomTot = 0.
    real, dimension(0:3) :: MomTot_in = 0.
    integer, save :: qBTot = 0, qSTot = 0
    integer :: qBTot_in = 0, qSTot_in = 0
    logical, save :: first = .true.
    integer :: qB,qS,qC,qIx2
    logical :: flagOK

    write(*,*) 'Check Conservation...'

    MomTot_in = 0.
    qBTot_in = 0
    qSTot_in = 0

    do iEns=lbound(particles,dim=1),ubound(particles,dim=1)
       do iPart=lbound(particles,dim=2),ubound(particles,dim=2)
          if (particles(iEns,iPart)%ID <= 0) cycle
          pPart => particles(iEns,iPart)
          MomTot_in = MomTot_in + pPart%momentum

          if (isBaryon(pPart%Id)) then
             if (.not.pPart%antiparticle) then
                qBTot_in = qBTot_in + 1
                qSTot_in = qSTot_in + hadron(pPart%Id)%strangeness
             else
                qBTot_in = qBTot_in - 1
                qSTot_in = qSTot_in - hadron(pPart%Id)%strangeness
             end if
          else if (isMeson(pPart%Id)) then
             if (.not.pPart%antiparticle) then
                qSTot_in = qSTot_in + hadron(pPart%Id)%strangeness
             else
                qSTot_in = qSTot_in - hadron(pPart%Id)%strangeness
             end if
          end if

       end do
    end do

    if (.not.first) then
       flagOK = .true.
       if (any((MomTot - MomTot_in)>1e-3)) then
          write(*,*)
          write(*,*) 'MomTot   :',MomTot
          write(*,*) 'MomTot_in:',MomTot_in
          write(*,*) 'Diff     :',MomTot-MomTot_in
          flagOK = .false.
       end if
       if (qBTot_in /= qBTot) then
          write(*,*)
          write(*,*) 'qBTot   :',qBTot
          write(*,*) 'qBTot_in:',qBTot_in
          flagOK = .false.
       end if
       if (qSTot_in /= qSTot) then
          write(*,*)
          write(*,*) 'qSTot   :',qSTot
          write(*,*) 'qSTot_in:',qSTot_in
          flagOK = .false.
       end if

       if (.not.flagOK) call TraceBack("Conservation violated!")

    end if

    first = .false.
    MomTot = MomTot_in
    qBTot = qBTot_in
    qSTot = qSTot_in

  end subroutine CheckConservation

end module checks
