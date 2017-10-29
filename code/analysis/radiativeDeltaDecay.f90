!******************************************************************************
!****m* /radiativeDeltaStorage
! NAME
! module radiativeDeltaStorage
!
! PURPOSE
! This module stores information about the initial Delta needed for the
! analysis of the radiative Delta decay.
!******************************************************************************
module radiativeDeltaStorage
  implicit none
  private

  !****************************************************************************
  !****t* radiativeDeltaStorage/tradiativeDeltaStorage
  ! SOURCE
  !
  type tradiativeDeltaStorage
     logical :: flagOK=.False.
     real    :: radius=0.
     real    :: absthreemomentum=0.
     real    :: mass=0.
     integer :: charge=0
  end type tradiativeDeltaStorage
  !
  ! PURPOSE
  ! This holds all information we want to store about the intial Delta.
  !****************************************************************************

  type(tradiativeDeltaStorage),save,dimension(:),allocatable :: DeltaTrackingInfo

  public :: radiativeDeltaStorage_Init,radiativeDeltaStorage_Clear, &
            radiativeDeltaStorage_Store,radiativeDeltaStorage_Get

contains


  !****************************************************************************
  !****s* radiativeDeltaStorage/radiativeDeltaStorage_Init
  ! NAME
  ! subroutine radiativeDeltaStorage_Init(NumInitialEvents)
  !
  ! PURPOSE
  ! Allocate memory and reset the corresponding arrays.
  !
  ! INPUTS
  ! * integer :: NumInitialEvents -- number of possible initial events
  !****************************************************************************
  subroutine radiativeDeltaStorage_Init(NumInitialEvents)

    integer,intent(in) :: NumInitialEvents

    if (allocated(DeltaTrackingInfo)) then
       deallocate(DeltaTrackingInfo)
    end if
    allocate(DeltaTrackingInfo(1:NumInitialEvents))
    DeltaTrackingInfo%flagOK=.false.
    DeltaTrackingInfo%radius=0.
    DeltaTrackingInfo%absthreemomentum=0.
    DeltaTrackingInfo%mass=0.
    DeltaTrackingInfo%charge=0

  end subroutine RadiativeDeltaStorage_Init


  !****************************************************************************
  !****s* radiativeDeltaStorage/radiativeDeltaStorage_clear
  ! NAME
  ! subroutine radiativeDeltaStorage_clear
  !
  ! PURPOSE
  ! If necessary clear allocated memory.
  !****************************************************************************
  subroutine radiativeDeltaStorage_clear

    !logical,save :: initFlag=.true.

    if (allocated(DeltaTrackingInfo)) then
       deallocate(DeltaTrackingInfo)
    end if

  end subroutine RadiativeDeltaStorage_clear


  !****************************************************************************
  !****s* radiativeDeltaStorage/radiativeDeltaStorage_Store
  ! NAME
  ! subroutine radiativeDeltaStorage_Store(i,value,value_rec)
  !
  ! PURPOSE
  ! Store the event info connected with number "i".
  !****************************************************************************
  subroutine radiativeDeltaStorage_Store(i,radius,absthreemomentum,mass,charge)

    integer,intent(in)          :: i
    real,   intent(in)          :: radius,absthreemomentum,mass
    integer,intent(in)          :: charge

    DeltaTrackingInfo(i)%flagOK=.true.
    DeltaTrackingInfo(i)%radius=radius
    DeltaTrackingInfo(i)%absthreemomentum=absthreemomentum
    DeltaTrackingInfo(i)%mass=mass
    DeltaTrackingInfo(i)%charge=charge

  end subroutine RadiativeDeltaStorage_Store


  !****************************************************************************
  !****f* radiativeDeltaStorage/radiativeDeltaStorage_Get
  ! NAME
  ! logical function radiativeDeltaStorage_Get(i,value,value_rec)
  !
  ! PURPOSE
  ! Get the event info stored connected with number "i".
  !
  ! INPUTS
  ! * integer :: i -- actual number of event
  !
  ! OUTPUT
  ! stored values
  !****************************************************************************
  logical function radiativeDeltaStorage_Get(i,radius,absthreemomentum,mass,charge)

    integer,intent(in)           :: i
    real,   intent(out)          :: radius,absthreemomentum,mass
    integer,intent(out)          :: charge

    radiativeDeltaStorage_Get = .FALSE.

    if (.not.allocated(DeltaTrackingInfo)) return

    radius=DeltaTrackingInfo(i)%radius
    absthreemomentum=DeltaTrackingInfo(i)%absthreemomentum
    mass=DeltaTrackingInfo(i)%mass
    charge=DeltaTrackingInfo(i)%charge

    radiativeDeltaStorage_Get = DeltaTrackingInfo(i)%flagOK

  end function RadiativeDeltaStorage_Get

end module RadiativeDeltaStorage




!******************************************************************************
!****m* /radiativeDeltaDecay
! NAME
! module radiativeDeltaDecay
! PURPOSE
! Calculate cross sections for radiative Delta decay.
!******************************************************************************
module radiativeDeltaDecay
  use histMC, only: histogramMC
  use inputGeneral, only: num_Runs_SameEnergy
  implicit none
  private

  !****************************************************************************
  !****g* radiativeDeltaDecay/Enable
  ! NAME
  ! Enable
  ! PURPOSE
  ! If .true. the dilepton analysis will be performed, otherwise not.
  ! SOURCE
  !
  logical, save :: Enable=.false.
  !****************************************************************************

  type(histogramMC),save:: msigma,psigma,rsigma

  real, save :: tsigmaphoton0=0.
  real, save :: tsigmapion0=0.
  real, save :: tsigmaphoton1=0.
  real, save :: tsigmapion1=0.


  public::DoRadiativeDeltaDecay,radiativeDeltaDecay_write_CS


contains


  !****************************************************************************
  !****s* radiativeDeltaDecay/readInput
  ! NAME
  ! subroutine readInput
  ! PURPOSE
  ! Reads input in jobcard out of namelist "radiative_DeltaDecay". Possible Input:
  ! * Enable
  !****************************************************************************
  subroutine readInput
    use output, only: Write_ReadingInput

    integer :: ios
    NAMELIST /radiative_DeltaDecay/ Enable

    call Write_ReadingInput('radiative_DeltaDecay',0)
    rewind(5)
    read(5,nml=radiative_DeltaDecay,IOSTAT=IOS)
    call Write_ReadingInput('radiative_DeltaDecay',0,IOS)

    write(*,*) 'Enable radiative_DeltaDecay?', Enable

    call Write_ReadingInput('radiative_DeltaDecay',1)
  end subroutine readInput


  !****************************************************************************
  !****f* radiativeDeltaDecay/DoRadiativeDeltaDecay
  ! NAME
  ! subroutine DoRadiativeDeltaDecay(nt,pertPart)
  ! PURPOSE
  !
  ! INPUTS
  ! * integer nt - time step number
  ! * type(particle) pertPart(:,:) - list of perturbative particles
  !****************************************************************************
  subroutine DoRadiativeDeltaDecay(nt,pertPart)
    use particleDefinition, only: particle
    use nucleusDefinition, only: tNucleus
    use inputGeneral, only: delta_T,numTimeSteps
    use constants, only: mN !, pi, alphaQED
    use IdTable, only: delta
    use nucleus, only: getTarget
    use baryonWidth, only: FullWidthBaryon
    use histMC, only: CreateHistMC,copyDesc
    use radiativeDeltaStorage

    ! INPUT PARAMETERS:
    integer,intent(in)::nt
    type(particle),intent(in):: pertPart(:,:)
    ! CONSTANTS:
    !    real,parameter::g=5.44 ! Delta decay
    integer,parameter::nevent=10 ! number of events for dilepton simulation (to enhance stat.)
    ! VARIABLES:
    integer::i,j,k,id,charge
    real::gamma,mass,perw,gam,gtot
    !real::mt,kaell,f,q0
    !logical::flags
    real::weight
    logical,save::first=.true.
    type(tNucleus),pointer::targNuc
    real :: extNucRadius
    real :: radius_ini,absthreemomentum_ini,mass_ini
    integer :: charge_ini

    if (first) then
       call readInput
       first=enable
    end if

    if (.not. Enable) return

    if (first) then
       first=.false.
       targNuc => getTarget()
       extNucRadius=targNuc%radius*1.2  !give 20% extra

       ! create histograms
       call CreateHistMC(msigma,'radiativeDelta Cross Section dSigma/dM [mubarn/GeV]',1.,2.0,0.01,4)
       call CreateHistMC(psigma,'radiativeDelta Cross Section dSigma/dP_lab [mubarn/GeV]',0.,2.,0.01,4)
       call CreateHistMC(rsigma,'radiativeDelta Cross Section dSigma/dradius [mubarn]',0.,extNucRadius,0.05,4)
       msigma%yDesc(1:4) = (/ 'photons from Delta0', 'photons from Delta+', 'pi0 from Delta0    ', 'pi0 from Delta+    ' /)
       call CopyDesc(psigma,msigma)
       call CopyDesc(rsigma,msigma)
       msigma%xDesc='Mass of initial Delta [GeV]'
       psigma%xDesc='Momentum of initial Delta [GeV]'
       rsigma%xDesc='Radius of initial Delta [fm]'

    end if


!!! loop for calculating cross sections
    do k=1,nevent
       do i=lbound(pertPart,1),ubound(pertPart,1) ! loop over ensembles
          do j=lbound(pertPart,2),ubound(pertPart,2) ! loop over all perturbative particles
             id=pertPart(i,j)%ID
             if (id==-1) exit ! end of current ensemble reached, jump to next one
             if (id<-1 .or. id>1000) then
                ! error: should not occur
                write(*,*) "found bad ID in radiativeDeltaDecay: ",id
                stop
             end if

             if (id.ne.delta) cycle

             charge=pertPart(i,j)%charge
             mass=pertPart(i,j)%mass
             gamma=pertPart(i,j)%momentum(0)/mass
             perw=pertPart(i,j)%perWeight

             if (.not.(mass>0.)) then
                write(*,*) "Warning: Bad mass in radiativeDeltaDecay!",id,mass,pertPart(i,j)%mass
                stop
             end if

             if (mass.le.mN) then
                write(*,*) 'mass.le.mN', mass
                write(*,*) 'decay not possible'
                write(*,*) 'offshell nucleons not yet implemented'
             end if

             if (charge/=1.and.charge/=0) cycle


!!$             ! G. Wolf formula
!!$             q0=(mass**2-mN**2)/2./mass
!!$             f=-1.5*(mass+mN)/mN/((mN+mass)**2)
!!$             mt=alphaQED*4.*pi*f**2*g**2*mass**2/9./mN*q0**2*(5.*mass-3.*(q0+mN))
!!$             kaell=max(0.,mass**4+mN**4-2.*mass**2*mN**2)
!!$             gam=sqrt(kaell)/16./pi/mass**2*mN*2*mt

             !PDG value -> VERY ROUGH ESTIMATE!!!!!
             gtot=max(FullWidthBaryon(id,mass),0.001)
             gam=0.0056*gtot

             if (nt/=numTimeSteps) then
                weight=perw*gam/gamma*delta_T/0.197/float(nevent)
             else
                gtot=max(FullWidthBaryon(id,mass),0.001)
                weight=perw*gam/gtot/float(nevent)
             end if

             if (.not.radiativeDeltaStorage_Get(pertPart(i,j)%firstEvent,radius_ini,absthreemomentum_ini,mass_ini,charge_ini)) then
                write(*,*) 'error in radiativeDeltaStorage_Get -> FALSE, stop'
                stop
             end if

             call CS(radius_ini,absthreemomentum_ini,mass_ini,charge_ini,weight)

          end do
       end do
    end do

  end subroutine DoRadiativeDeltaDecay


  subroutine CS(radius,mom,mass,charge,pw)
    ! calculate cross sections (in microbarn/GeV)
    use inputGeneral, only: num_energies
    use histMC, only: AddHistMC

    ! input variables
    real,intent(in)::radius,mom,mass
    real,intent(in)::pw !perturbative weight
    integer, intent(in) :: charge


    if (num_energies.gt.1) then
       write(*,*) 'warning radiative Delta:num_energies.gt.1 -> not implemented yet -> stop'
       stop
    end if

    !cross sections (already averaged over num_runs_same_energy)

    call AddHistMC(msigma,mass,charge+1,pw/float(num_Runs_SameEnergy))
    call AddHistMC(psigma,mom,charge+1,pw/float(num_Runs_SameEnergy))
    call AddHistMC(rsigma,radius,charge+1,pw/float(num_Runs_SameEnergy))

    if (charge.eq.1) tsigmaphoton1=tsigmaphoton1+pw/float(num_Runs_SameEnergy)
    if (charge.eq.0) tsigmaphoton0=tsigmaphoton0+pw/float(num_Runs_SameEnergy)

  end subroutine CS


  subroutine radiativeDeltaDecay_write_CS(Particles)
    ! writes out the cross section (called at the end of the simulation)
    use histMC, only: AddHistMC,WriteHistMC
    use particleDefinition
    use Idtable, only: pion
    use radiativeDeltaStorage

    type(particle), intent(in), dimension(:,:) :: Particles
    integer :: i,j
    real :: radius,mom,mass
    integer :: charge
    real :: pw


    if (.not.enable) return

    do i=lbound(Particles,dim=1),ubound(Particles,dim=1)
       do j=lbound(Particles,dim=2),ubound(Particles,dim=2)
          if (Particles(i,j)%ID.ne.pion) cycle
          if (Particles(i,j)%charge.ne.0) cycle

          if (.not.radiativeDeltaStorage_Get(Particles(i,j)%firstEvent,radius,mom,mass,charge)) then
             write(*,*) 'error in radiativeDeltaStorage_Get -> FALSE, stop'
             stop
          end if

          pw=Particles(i,j)%perweight/float(num_Runs_SameEnergy)

          call AddHistMC(msigma,mass,charge+1+2,pw)
          call AddHistMC(psigma,mom,charge+1+2,pw)
          call AddHistMC(rsigma,radius,charge+1+2,pw)

          if (charge.eq.1) tsigmapion1=tsigmapion1+pw
          if (charge.eq.0) tsigmapion0=tsigmapion0+pw

       end do
    end do

    call WriteHistMC(msigma,'radiativeDeltaMass.dat')
    call WriteHistMC(psigma,'radiativeDeltaMom.dat')
    call WriteHistMC(rsigma,'radiativeDeltaRadius.dat')

    open(10,File='radiativeDeltaTotal.dat')
    write(10,*) '# photons from Delta0, photons from Delta+, pi0 from Delta0, pi0 from Delta+,vac-rate pi0/Delta'
    write(10,'(10g16.8)') tsigmaphoton0,tsigmaphoton1,tsigmapion0,tsigmapion1,(2./3.)/0.0056
    close(10)

    !note: vacuum estimate: Delta+->pi0 p clebsch-gordon sqrt(2/3) -> 2/3 * gammatot
    !                       Delta0->pi0 n                sqrt(2/3) -> 2/3 * gammatot
    !PDG value for radiative decay: 0.0056 * gammatot
    !therefore rate pi0/Delta=(2/3)/0.0056


  end subroutine RadiativeDeltaDecay_write_CS


end module radiativeDeltaDecay
