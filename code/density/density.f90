!******************************************************************************
!****m* /densitymodule
! NAME
! module densitymodule
!
! PURPOSE
! Administrates the calculation of the density.
!******************************************************************************
module densitymodule
  use dichteDefinition
  use constants, only: singlePrecision

  implicit none

  private

  !****************************************************************************
  !****g* densitymodule/densitySwitch
  ! SOURCE
  !
  integer, save :: densitySwitch=1
  ! PURPOSE
  ! This switch decides whether the density is static or dynamic during
  ! the run. ("Static" makes sense only for fixed target scenarios!)
  !
  ! One can use a static density if the nucleus stays roughly in its
  ! ground state during the collision.
  !
  ! possible values:
  ! * 0: Density is set to 0.
  ! * 1: Dynamic density according to test-particle distribution.
  ! * 2: Static density (not for heavy-ion collisions).
  ! * 3: Resting matter: Density is given by the two input parameters
  !   "densityInput_neutron" and "densityInput_proton".
  !****************************************************************************


  !****************************************************************************
  !****g* densitymodule/linearInterpolation
  ! SOURCE
  !
  logical, save :: linearInterpolation = .true.
  ! PURPOSE
  ! If this switch is 'true', then the dynamic-density mode uses linear
  ! interpolation to determine the density in between the gridpoints.
  !****************************************************************************


  !Input parameters for the grid where density is calculated upon :

  !****************************************************************************
  !****g* densitymodule/gridSize
  ! SOURCE
  !
  real, dimension(1:3), save, public :: gridSize = (/12.,12.,12./)
  ! PURPOSE
  ! Size of density grid in fm.
  !****************************************************************************

  !****************************************************************************
  !****g* densitymodule/gridPoints
  ! SOURCE
  !
  integer, dimension(1:3), save, public :: gridPoints = (/30,30,30/)
  ! PURPOSE
  ! Number of gridpoints in each space direction.
  !****************************************************************************


  !****************************************************************************
  !****g* densitymodule/gridSpacing
  ! SOURCE
  !
  real, dimension(1:3), save, public :: gridSpacing = 1.
  ! PURPOSE
  ! Spacing of the grid in each dimension, determined by gridSize/gridPoints.
  !****************************************************************************


  !****************************************************************************
  !****g* densitymodule/densityInput_proton
  ! SOURCE
  !
  real, save :: densityInput_proton=0.084
  ! PURPOSE
  ! Assumed proton density if densitySwitch=3
  !****************************************************************************

  !****************************************************************************
  !****g* densitymodule/densityInput_neutron
  ! SOURCE
  !
  real, save  :: densityInput_neutron=0.084
  ! PURPOSE
  ! Assumed neutron density if densitySwitch=3
  !****************************************************************************

  !****************************************************************************
  !****g* densitymodule/setnewsmearing
  ! SOURCE
  !
  logical, save :: setnewsmearing=.false.
  ! PURPOSE
  ! Readjust the smearing to a different width if .true.
  !****************************************************************************

  !****************************************************************************
  !****g* densitymodule/newsmearing
  ! SOURCE
  !
  real, save  :: newsmearing=1.
  ! PURPOSE
  ! Use a smearing width as in a grid wtih newsmearing times the gridspacing
  !****************************************************************************


  type(dichte),          save, Allocatable, dimension(:,:,:), public :: densityField      ! Field where density is saved in
  real(singlePrecision), save, Allocatable, dimension(:,:,:), public :: totalDensity      ! density of baryons plus density of antibaryons (fm^-3)

  ! Fields needed when the RMF is used:
  real,                  save, Allocatable, dimension(:,:,:),   public :: sigmaField      ! sigma-meson field (GeV)
  real(singlePrecision), save, Allocatable, dimension(:,:,:,:), public :: omegaField      ! omega-meson field (GeV)
  real(singlePrecision), save, Allocatable, dimension(:,:,:,:), public :: rhoField        ! rho-meson field (GeV)

  real(singlePrecision), save, Allocatable, dimension(:,:,:,:) :: baryonCurrent_old  ! baryon density    from previous time step (fm^-3)
  real(singlePrecision), save, Allocatable, dimension(:,:,:,:) :: omegaField_old     ! omega-meson field from previous time step (GeV)
  real(singlePrecision), save, Allocatable, dimension(:,:,:,:) :: rhoField_old       ! rho-meson field   from previous time step (GeV)

  real(singlePrecision), save, Allocatable, dimension(:,:,:) :: scalarDensity   ! scalar density (fm^-3)
  real(singlePrecision), save, Allocatable, dimension(:,:,:) :: d_scalarDensity_dsigma    ! partial derivative of the scalar density over sigma field (fm^-3/GeV)
  real,                  save, Allocatable, dimension(:,:,:) :: K_rmf, Diag_rmf, U_rmf    ! auxiliary arrays needed when RMF gradients are included

  real(singlePrecision), save, Allocatable, dimension(:,:,:,:), public :: fourMomentumDensity     ! four-momentum density field (not used in propagation)


  !****************************************************************************
  !****g* densitymodule/numberLargePoints
  ! SOURCE
  !
  integer, save, public :: numberLargePoints=2
  ! PURPOSE
  ! Number of points which are considered to the left and right to smear
  ! density on
  !****************************************************************************

  integer, parameter,public :: smallerGridPoints=2 !number of Gridpoints in small grid

  !normalized weights for smearing
  real, save,dimension(:,:),allocatable,public :: smearingWeights

  logical, save  :: initFlag=.true.   ! Checks whether input was already read in.

  type(dichte), parameter :: densZero = dichte(0.,0.,0.)

!!$  PUBLIC :: hadronicPotEnergy
  public :: densityAt
  public :: getBaryonDensity
  public :: updateDensity
  public :: updateRMF
  public :: storeFields
  public :: energyDeterminationRMF
  public :: Particle4Momentum_RMF, true4Momentum_RMF
  public :: SelfEnergy_scalar, SelfEnergy_vector, SelfEnergy_vector_old, DiracMass
  public :: boostToLRF
  public :: acceptGrid
  public :: getGridSpacing, getGridPoints, get_realGridSpacing
  public :: fermiMomentum_sym, fermiMomentum_noIsospin
  public :: get_densitySwitch, set_densitySwitch, set_densityInput
  public :: cleanup
  public :: FermiMomAt
  public :: GetGridIndex
  public :: writeDensityPlane

contains

!!$  !*************************************************************************
!!$  !****f* densitymodule/hadronicPotEnergy
!!$  ! NAME
!!$  ! real function hadronicPotEnergy(alpha,beta,tau)
!!$  ! PURPOSE
!!$  ! Evaluate the energy wich is stored in the hadronic potential for a potential which is
!!$  ! not momentum dependend.
!!$  ! INPUTS
!!$  ! * real alpha,beta,tau --  Potential parameters
!!$  ! RESULT
!!$  ! * Potential Energy in GEV
!!$  !*************************************************************************
!!$  function hadronicPotEnergy(alpha,beta,tau)
!!$    use constants, only : rhoNull
!!$
!!$    real :: hadronicPotEnergy
!!$    real, intent(in) :: alpha,beta,tau ! parameters of hadronic potential
!!$
!!$    ! Integrate the energy density
!!$    integer :: index_X,index_Y,index_Z
!!$    Do index_X=-gridPoints(1),gridPoints(1)
!!$       Do index_Y=-gridPoints(2),gridPoints(2)
!!$          Do index_Z=-gridPoints(3),gridPoints(3)
!!$             hadronicPotEnergy=hadronicPotEnergy &
!!$                  & +(alpha/2.*densityField(index_X,index_Y,index_Z)%baryon(0)**2/rhoNull &
!!$                  & +beta/(tau+1.)*(densityField(index_X,index_Y,index_Z)%baryon(0)**(tau+1.)) &
!!$                  & /(rhoNull**tau))*gridSpacing(1)*gridSpacing(2)*gridSpacing(3)
!!$          End do
!!$       End do
!!$    End do
!!$    hadronicPotEnergy=hadronicPotEnergy/1000. !in units of GeV
!!$  end function hadronicPotEnergy

  !****************************************************************************
  !****s* densitymodule/init
  ! NAME
  ! subroutine init
  !
  ! PURPOSE
  ! read the namelist, initalize (parts of) the module
  !****************************************************************************
  subroutine init
    use output, only: Write_ReadingInput
    use RMF, only: getRMF_flag, grad_flag, fourMomDen_flag
    use nucleusDefinition
    use nucleus, only: getTarget
    use inputGeneral, only: eventType
    use eventtypes, only: elementary, Box

    integer :: ios
    type(tNucleus), pointer :: targetNuc

    !**************************************************************************
    !****n* densitymodule/initDensity
    ! NAME
    ! NAMELIST /initDensity/
    ! PURPOSE
    ! Includes the input switches and variables:
    ! * densitySwitch
    ! * linearInterpolation
    ! * densityInput_proton
    ! * densityInput_neutron
    ! * gridSize
    ! * gridPoints
    ! * setnewsmearing
    ! * newsmearing
    ! * numberLargePoints
    !**************************************************************************
    NAMELIST /initDensity/ densitySwitch, linearInterpolation, densityInput_proton, densityInput_neutron, &
                           gridSize, gridPoints, setnewsmearing, newsmearing, numberLargePoints

    call Write_ReadingInput('initDensity',0)
    rewind(5)
    read(5,nml=initDensity,iostat=ios)
    call Write_ReadingInput('initDensity',0,ios)

    select case (eventType)
    case (elementary, Box)
      densitySwitch = 0
      write(*,*) 'densitySwitch is set to 0 for elementary target/Box'
    case default
      targetNuc => getTarget()
      if (targetNuc%mass==1) then
         densitySwitch = 0
         write(*,*) 'densitySwitch is set to 0 for nucleon target'
      end if
    end select

    write(*,'(A,I5,A)') ' Set densitySwitch to ',densitySwitch,'.'

    if (densitySwitch==3) then
       write(*,'(A)')        '   => We assume resting matter.'
       write(*,'(A,F5.3,A)') '      proton density:  ',densityInput_proton, '/fm^3.'
       write(*,'(A,F5.3,A)') '      neutron density: ',densityInput_neutron,'/fm^3.'
    end if
    write(*,'(A,L6,A)') ' Set linearInterpolation to  ',linearInterpolation,'.'
    if (setnewsmearing) then
       write(*,*) 'Smearingwidth set to ', newsmearing
    end if

    ! Set grid spacing and allocate vectors
    gridSpacing=gridSize/float(gridPoints)

    write(*,'(A,2(F6.2,","),F6.2,A)') ' The gridsize of the density grid is        = (', gridsize, ') fm'
    write(*,'(A,2(I6,","),I6,A)')     ' The number of gridpoints per dimension are = (', gridPoints, ') '
    write(*,'(A,2(F6.2,","),F6.2,A)') ' The grid spacing is                        = (', gridSpacing, ') fm'

    allocate(densityField(-gridPoints(1):gridPoints(1),-gridPoints(2):gridPoints(2),-gridPoints(3):gridPoints(3)))
    allocate(totalDensity(-gridPoints(1):gridPoints(1),-gridPoints(2):gridPoints(2),-gridPoints(3):gridPoints(3)))

    if ( getRMF_flag() ) then ! RMF used
        allocate(sigmaField(-gridPoints(1):gridPoints(1),-gridPoints(2):gridPoints(2),&
                         &  -gridPoints(3):gridPoints(3)))
        allocate(omegaField(-gridPoints(1):gridPoints(1),-gridPoints(2):gridPoints(2),&
                         &  -gridPoints(3):gridPoints(3),0:3))
        allocate(rhoField(-gridPoints(1):gridPoints(1),-gridPoints(2):gridPoints(2),&
                         &  -gridPoints(3):gridPoints(3),0:3))
        allocate(scalarDensity(-gridPoints(1):gridPoints(1),-gridPoints(2):gridPoints(2),&
                         &     -gridPoints(3):gridPoints(3)))
        allocate(d_scalarDensity_dsigma(-gridPoints(1):gridPoints(1),&
                         &              -gridPoints(2):gridPoints(2),&
                         &              -gridPoints(3):gridPoints(3)))
        allocate(omegaField_old(-gridPoints(1):gridPoints(1),&
                         &         -gridPoints(2):gridPoints(2),&
                         &         -gridPoints(3):gridPoints(3),0:3))
        allocate(rhoField_old(-gridPoints(1):gridPoints(1),&
                         &         -gridPoints(2):gridPoints(2),&
                         &         -gridPoints(3):gridPoints(3),0:3))
        allocate(baryonCurrent_old(-gridPoints(1):gridPoints(1),&
                         &         -gridPoints(2):gridPoints(2),&
                         &         -gridPoints(3):gridPoints(3),1:3))
        if (grad_flag) then
           allocate(K_rmf(-gridPoints(1):gridPoints(1),-gridPoints(2):gridPoints(2),&
                       &  -gridPoints(3):gridPoints(3)))
           allocate(Diag_rmf(-gridPoints(1):gridPoints(1),-gridPoints(2):gridPoints(2),&
                          &  -gridPoints(3):gridPoints(3)))
           allocate(U_rmf(-gridPoints(1):gridPoints(1),-gridPoints(2):gridPoints(2),&
                          &  -gridPoints(3):gridPoints(3)))
        end if
        if (fourMomDen_flag) then
          allocate(fourMomentumDensity(-gridPoints(1):gridPoints(1),&
                           &           -gridPoints(2):gridPoints(2),&
                           &           -gridPoints(3):gridPoints(3),0:3))
        end if
    end if

    call Write_ReadingInput('initDensity',1)
    initFlag=.false.

  end subroutine init


  subroutine cleanUp
    use RMF, only: getRMF_flag, grad_flag, fourMomDen_flag

    if (allocated(densityField)) deallocate(densityField, totalDensity)

    if ( getRMF_flag() ) then ! RMF used
        deallocate(sigmaField,omegaField,rhoField,scalarDensity,d_scalarDensity_dsigma,omegaField_old, &
             & rhoField_old, baryonCurrent_old)
        if (grad_flag) then
           deallocate(K_rmf,Diag_rmf,U_rmf)
        end if
        if (fourMomDen_flag) then
          deallocate(fourMomentumDensity)
        end if
    end if
    if (allocated(SmearingWeights)) deallocate(SmearingWeights)
  end subroutine


  !****************************************************************************
  !****s* densitymodule/acceptGrid
  ! NAME
  ! subroutine  acceptGrid(GridSpacing_x,GridSpacing_y,GridSpacing_z)
  !
  ! PURPOSE
  ! Accepts the grid spacings and recalculates
  ! the numbers of points in each direction. Needed for relativistic HIC,
  ! when initial nuclei are Lorentz-contracted. Called from the module
  ! initHeavyIon.
  !****************************************************************************
  subroutine acceptGrid(GridSpacing_x,GridSpacing_y,GridSpacing_z)
  use RMF, only: getRMF_flag, grad_flag, fourMomDen_flag
  use nucleusDefinition
  use nucleus, only: getTarget, getProjectile
  use inputGeneral, only: time_max

  type(tNucleus), pointer :: targetNuc, projectileNuc
  real, dimension(1:3) :: PosT, PosP, checkGrid

  real, intent(in) :: GridSpacing_x,GridSpacing_y,GridSpacing_z

  if (initFlag) call init

  if ( GridSpacing_x <= 0.) then
    write(*,*) ' acceptGrid: Wrong grid spacing along x-axis:', GridSpacing_x
    stop
  end if

  if ( GridSpacing_y <= 0.) then
    write(*,*) ' acceptGrid: Wrong grid spacing along y-axis:', GridSpacing_y
    stop
  end if

  if ( GridSpacing_z <= 0.) then
    write(*,*) ' acceptGrid: Wrong grid spacing along z-axis:', GridSpacing_z
    stop
  end if

  GridSpacing(1:3) = (/ GridSpacing_x, GridSpacing_y, GridSpacing_z /)

  ! Adjust gridSize for HIC;
  ! x- and y-directions: nuclei-radius + surface + 5fm
  ! z-direction in addition: (HIC) grid extension up to Pos_z = V_beam*time_max
  targetNuc => getTarget()
  projectileNuc => getProjectile()
  PosT(1:3) = targetNuc%position(1:3) + targetNuc%radius+targetNuc%surface+5.
  PosP(1:3) = projectileNuc%position(1:3) + projectileNuc%radius+projectileNuc%surface+5.
  PosT(3) = PosT(3) + targetNuc%velocity(3)*time_max
  PosP(3) = PosP(3) + projectileNuc%velocity(3)*time_max
  checkGrid(1:3) = aint( max(PosT(1:3),PosP(1:3)) )

  !check that nuclei will stay inside the spatial grid
  if ( gridSize(1) < checkGrid(1) .or. gridSize(2) < checkGrid(2) .or. &
     & gridSize(3) < checkGrid(3) ) then

     write(*,*) '  nuclei partially outside of spatial grid!'
     write(*,'(A,2(F6.2,","),F6.2,A)') '  checkGrid        = (', checkGrid, ') fm'
     write(*,*) '  readjust gridSize too'

     gridSize(1:3) = checkGrid(1:3)

  end if

  gridPoints = nint(gridSize/gridSpacing)

  write(*,*) ' New grid parameters:'
  write(*,'(A,2(F6.2,","),F6.2,A)') '  The gridsize of the density grid is        = (', gridsize, ') fm'
  write(*,'(A,2(I6,","),I6,A)')     '  The number of gridpoints per dimension are = (', gridPoints, ') '
  write(*,'(A,2(F6.2,","),F6.2,A)') '  The grid spacing is                        = (', gridSpacing, ') fm'

  deallocate(densityField,totalDensity)
  allocate(densityField(-gridPoints(1):gridPoints(1),-gridPoints(2):gridPoints(2),&
         &  -gridPoints(3):gridPoints(3)))
  allocate(totalDensity(-gridPoints(1):gridPoints(1),&
                       &-gridPoints(2):gridPoints(2),&
                       &-gridPoints(3):gridPoints(3)))

  if ( getRMF_flag() ) then ! RMF used

      deallocate( sigmaField, omegaField, rhoField, &
                & scalarDensity, d_scalarDensity_dsigma, &
                & omegaField_old, rhoField_old, baryonCurrent_old)

      allocate(sigmaField(-gridPoints(1):gridPoints(1),-gridPoints(2):gridPoints(2),&
                       &  -gridPoints(3):gridPoints(3)))
      allocate(omegaField(-gridPoints(1):gridPoints(1),-gridPoints(2):gridPoints(2),&
                       &  -gridPoints(3):gridPoints(3),0:3))
      allocate(rhoField(-gridPoints(1):gridPoints(1),-gridPoints(2):gridPoints(2),&
                       &  -gridPoints(3):gridPoints(3),0:3))
      allocate(scalarDensity(-gridPoints(1):gridPoints(1),-gridPoints(2):gridPoints(2),&
                       &     -gridPoints(3):gridPoints(3)))
      allocate(d_scalarDensity_dsigma(-gridPoints(1):gridPoints(1),&
                       &              -gridPoints(2):gridPoints(2),&
                       &              -gridPoints(3):gridPoints(3)))
      allocate(omegaField_old(-gridPoints(1):gridPoints(1),&
                       &      -gridPoints(2):gridPoints(2),&
                       &      -gridPoints(3):gridPoints(3),0:3))
      allocate(rhoField_old(-gridPoints(1):gridPoints(1),&
                       &      -gridPoints(2):gridPoints(2),&
                       &      -gridPoints(3):gridPoints(3),0:3))
        allocate(baryonCurrent_old(-gridPoints(1):gridPoints(1),&
                         &         -gridPoints(2):gridPoints(2),&
                         &         -gridPoints(3):gridPoints(3),1:3))
      if (grad_flag) then
           deallocate(K_rmf, Diag_rmf, U_rmf)
           allocate(K_rmf(-gridPoints(1):gridPoints(1),-gridPoints(2):gridPoints(2),&
                       &  -gridPoints(3):gridPoints(3)))
           allocate(Diag_rmf(-gridPoints(1):gridPoints(1),-gridPoints(2):gridPoints(2),&
                          &  -gridPoints(3):gridPoints(3)))
           allocate(U_rmf(-gridPoints(1):gridPoints(1),-gridPoints(2):gridPoints(2),&
                          &  -gridPoints(3):gridPoints(3)))
      end if
      if (fourMomDen_flag) then
          deallocate(fourMomentumDensity)
          allocate(fourMomentumDensity(-gridPoints(1):gridPoints(1),&
                           &           -gridPoints(2):gridPoints(2),&
                           &           -gridPoints(3):gridPoints(3),0:3))
      end if
  end if

  end subroutine acceptGrid


  !****************************************************************************
  !****f* densitymodule/getBaryonDensity
  ! NAME
  ! real function getBaryonDensity (pos)
  ! PURPOSE
  ! Evaluate the baryon density at a certain position.
  !****************************************************************************
  real function getBaryonDensity (pos)
    use idTable, only: nucleon
    use RMF, only: getRMF_flag, ModificationFactor
    real, dimension(1:3),intent(in) :: pos

    integer :: i(1:3)
    real :: factor

    i(:) = nint(pos(:)/gridSpacing(:))

    if (all(abs(i(:))<=gridPoints(:))) then

      if (getRMF_flag()) then
        ! Modification factor of the antinucleon coupling constants
        factor = ModificationFactor(nucleon,.true.)
      else
        ! Remember: w/o RMF totalDensity is (baryon-antibaryon) density
        factor = 1.
      end if

      ! density of the baryons (excluding antibaryons)
      getBaryonDensity = (densityField(i(1),i(2),i(3))%baryon(0) + factor*totalDensity(i(1),i(2),i(3))) / (1.+factor)

    else

      getBaryonDensity = 0.

    end if

  end function getBaryonDensity


  !****************************************************************************
  !****f* densitymodule/densityAt
  ! NAME
  ! type(dichte) function densityAt(r)
  ! PURPOSE
  ! Evaluate density at some space-point r.
  ! INPUTS
  !  * real, dimension(1:3), intent(in) :: r   ---   position where density should be calculated
  !****************************************************************************
  type(dichte) function densityAt(r)
    use densityStatic, only: staticDensity
    use nucleus, only: getTarget
    use inputGeneral, only: continousBoundaries

    real, dimension(1:3), intent(in) :: r

    real,dimension(1:3) :: rNew
    integer :: indizes(1:3), i

    if (initFlag) call init

    select case (densitySwitch)
    case (0) ! Assume no density

       densityAt = densZero

    case (1) ! dynamic density according to test-particle distribution

       if (GetGridIndex(r,indizes,0)) then
          ! inside grid
          if (linearInterpolation) then
            densityAt = interpolate(r)
          else ! no interpolation
            densityAt = densityField(indizes(1),indizes(2),indizes(3))
          end if
       else
          ! outside grid
          densityAt = densZero

          if (continousBoundaries) then
             ! Implement continous boundary conditions
             rNew=r
             do i=1,3
                if (r(i)>gridSize(i)) then
                   rNew(i)=r(i)-2.*gridsize(i)
                else if (r(i)<-gridSize(i)) then
                   rNew(i)=r(i)+2.*gridsize(i)
                end if
             end do
             if (GetGridIndex(rNew,indizes,0)) then
                ! inside grid
                if (linearInterpolation) then
                  densityAt = interpolate(rNew)
                else ! no interpolation
                  densityAt = densityField(indizes(1),indizes(2),indizes(3))
                end if
             end if
          end if
       end if

    case (2) ! static density of a resting target

       densityAt=staticDensity(r,getTarget())

    case (3) ! Use input density

       densityAt%baryon(0)= densityInput_proton+densityInput_neutron
       densityAt%neutron(0)=densityInput_neutron
       densityAt%proton(0)= densityInput_proton
       densityAt%baryon(1:3)= 0.
       densityAt%neutron(1:3)=0.
       densityAt%proton(1:3)= 0.

    case default
       write(*,*) 'This density switch is not valid:',densitySwitch
       write(*,*) 'Stop'
       stop

    end select



  contains

    !**************************************************************************
    ! For linear interpolation between the grid points.
    !**************************************************************************
    type(dichte) function interpolate (r)
      use inputGeneral, only: continousBoundaries

      real, dimension(1:3), intent(in) :: r  ! position where density should be calculated
      integer :: i,j,k,lowIndex(1:3)
      integer,dimension(0:7,1:3) :: grid     ! field with all corners of 3D-box where r is situated in
      real,dimension(1:3) :: gridPos         ! position of gridpoints
      real :: factor

      ! The point "r" where the density should be calculated
      ! is sitting in a 3D-box with lowest corner "LowIndex" on the
      ! density grid. We first construct this point. Then all other corners of the 3D-box
      ! are constructed.
      ! In the end a simple linear interpolation is used to make the density smooth
      ! inside the box.

      ! (1.) Construct Lowest lying point (most negative point!)
      LowIndex=Int(r/gridSpacing)  ! always chooses integers which are closer to the origin
      do i=1,3                     ! convert them to integers which are more negative!
         if (r(i)/gridsize(i).lt.0) then
            LowIndex(i)=LowIndex(i)-1
         end if
      end do

      ! (2.) Define GridPoints on the corners of the 3D box in which the point "r" is situated
      do i=0,1
         do j=0,1
            do k=0,1
               grid(i+2*j+4*k,1:3)=LowIndex+(/i,j,k/)
            end do
         end do
      end do

      ! (3.) Do linear interpolation
      interpolate = densZero

      do i=0,7
         gridPos(1:3)=grid(i,1:3)*gridSpacing(1:3) !position of grid point
         factor=1.  ! evaluate weight for linear interpolation for each grid point
         do j=1,3
            factor=factor*Abs(gridSpacing(j)-abs(r(j)-gridPos(j)))/gridspacing(j)
         end do
         if (continousBoundaries) then
            ! Implement continuous boundary conditions
            do j=1,3
               if (grid(i,j).lt.-gridPoints(j)) then
                  grid(i,j)=grid(i,j)+2*gridPoints(j)
               else if (grid(i,j).gt.gridPoints(j)) then
                  grid(i,j)=grid(i,j)-2*gridPoints(j)
               end if
            end do
         end if
         if ((abs(grid(i,1)).le.gridPoints(1)).and.(Abs(grid(i,2)).le.gridPoints(2))&
              & .and.(abs(grid(i,3)).le.gridPoints(3))) then
            interpolate=interpolate+(densityField(grid(i,1),grid(i,2),grid(i,3))*factor)
         end if
      end do
    end function interpolate

  end function densityAt



  !****************************************************************************
  !  Dynamic DENSITY
  !****************************************************************************

  !****************************************************************************
  !****f* densitymodule/updateDensity
  ! NAME
  ! subroutine updateDensity(parts)
  ! PURPOSE
  ! Updates the vector densityField which is used by densityAt and stores the
  ! density of the testparticles.
  ! INPUTS
  ! * type(particle) parts(:,:)
  !
  ! NOTE
  ! This routine does not gets faster when using OpenMP.
  !****************************************************************************
  subroutine updateDensity(Parts)

    use particledefinition
    use inputGeneral, only: continousboundaries
    use output, only: DoPr
    use idtable, only: isBaryon

    type(particle),dimension(:,:),intent(in),target :: Parts

    !Flag that controls whether density weights are already calculated
    logical, save :: FirstTime=.true.


    integer :: small,large
    integer :: index1, index2, index3
    integer,dimension(1:3) :: posOrig
    integer,dimension(1:3) :: indexSmall, indexC

    integer :: bEdge ! binary edge code
    type(particle),pointer :: pPart
    integer :: iEns,iPart, nEns,nPart
    real :: SW

    if (initFlag) call init

    if (FirstTime) then
       call initDensityWeights
       call FillStaticDensity
       FirstTime=.false.
    end if

    if (densitySwitch.ne.1) return

    if (DoPr(2)) write(*,*) 'Updating density field'

    densityField = densZero
    totalDensity = 0.


    nEns  = size(Parts,dim=1)
    nPart = size(Parts,dim=2)

    do iPart=1,nPart
       do iEns=1,nEns

          pPart => Parts(iEns,iPart)

          if (pPart%Id <  0) exit
          if (pPart%id <= 0) cycle
          if (.not.isBaryon(pPart%id)) cycle

          ! position in large grid:
          posOrig=NINT(pPart%position(1:3)/gridSpacing)
          ! position in small grid:
          indexSmall=NINT((pPart%position(1:3)/gridSpacing-float(posOrig)) &
               * (2.*float(smallerGridPoints)+0.999999))

          ! Test for errors:
          if (any((abs(indexSmall).gt.smallerGridPoints))) then
             write(*,*) 'Problem in updateDensity, module density.f90'
             write(*,*) IndexSmall, 'too big, this must not happen !!!'
             write(*,*) pPart%id
             write(*,*) pPart%position
             write(*,*) posOrig
             stop
          end if

          small=1+(smallerGridPoints+indexSmall(3))                               &
               & +(smallerGridPoints+indexSmall(2))*(2*smallerGridPoints+1)       &
               & +(smallerGridPoints+indexSmall(1))*(2*smallerGridPoints+1)**2
          large=0

          ! Smearing particle over points in Neighborhood:
          do index1=posOrig(1)-numberLargePoints,posOrig(1)+numberlargePoints
             do index2=posOrig(2)-numberLargePoints,posOrig(2)+numberlargePoints
                do index3=posOrig(3)-numberLargePoints,posOrig(3)+numberlargePoints
                   large=large+1
                   indexC = (/index1,index2,index3/)

                   if (continousBoundaries) then
                      ! Implement continous boundary conditions,
                      ! so every point outside the grid is attributed to a point inside.

                      where (indexC >  gridPoints) indexC = indexC - 2*gridPoints
                      where (indexC < -gridPoints) indexC = indexC + 2*gridPoints

                   end if

                   SW = smearingWeights(small,large)

                   if ( (all(abs(indexC) <= gridPoints)) &   ! Point is inside Grid
                        &.and.(SW.gt.0.0) ) then

                      bEdge = 0

                      if (continousBoundaries) then
                         ! When having continous boundaries, then every contribution
                         ! to a point at the edge of the grid must also be attributed
                         ! to its opposing point on the other side.
                         if (abs(indexC(1)) == gridPoints(1)) bEdge = bEdge+1
                         if (abs(indexC(2)) == gridPoints(2)) bEdge = bEdge+2
                         if (abs(indexC(3)) == gridPoints(3)) bEdge = bEdge+4
                      end if

                      call addTodensityfield( indexC )

                      ! the following may be optimized further!!!!
                      select case (bEdge)
                      case (1)
                         call addTodensityfield( indexC * (/ -1, 1, 1 /) )
                      case (2)
                         call addTodensityfield( indexC * (/  1,-1, 1 /) )
                      case (3)
                         call addTodensityfield( indexC * (/ -1, 1, 1 /) )
                         call addTodensityfield( indexC * (/  1,-1, 1 /) )
                         call addTodensityfield( indexC * (/ -1,-1, 1 /) )
                      case (4)
                         call addTodensityfield( indexC * (/  1, 1,-1 /) )
                      case (5)
                         call addTodensityfield( indexC * (/ -1, 1, 1 /) )
                         call addTodensityfield( indexC * (/  1, 1,-1 /) )
                         call addTodensityfield( indexC * (/ -1, 1,-1 /) )
                      case (6)
                         call addTodensityfield( indexC * (/  1,-1, 1 /) )
                         call addTodensityfield( indexC * (/  1, 1,-1 /) )
                         call addTodensityfield( indexC * (/  1,-1,-1 /) )
                      case (7)
                         call addTodensityfield( indexC * (/ -1, 1, 1 /) )
                         call addTodensityfield( indexC * (/  1,-1, 1 /) )
                         call addTodensityfield( indexC * (/ -1,-1, 1 /) )
                         call addTodensityfield( indexC * (/  1, 1,-1 /) )
                         call addTodensityfield( indexC * (/ -1, 1,-1 /) )
                         call addTodensityfield( indexC * (/  1,-1,-1 /) )
                         call addTodensityfield( indexC * (/ -1,-1,-1 /) )
                      end select

                   end if

                end do ! loop Index1
             end do ! loop Index2
          end do ! loop Index3

       end do ! loop particle
    end do ! loop ensembles

  contains

    !**************************************************************************


    subroutine addTodensityfield(ind)

      use idtable, only: nucleon, isBaryon
      use RMF, only: getRMF_flag, ModificationFactor

      integer, dimension(3), intent(in) :: ind
      real :: factor
      integer :: I1,I2,I3
      real,dimension(0:3) :: vec

      ! abbrev:
      I1 = ind(1)
      I2 = ind(2)
      I3 = ind(3)

      vec(0) = 1.
      vec(1:3) = pPart%velocity(1:3)
      vec = vec*SW

      if (getRMF_flag()) then

         ! Modification factor for the coupling constants (always positive !!!):
         factor = ModificationFactor(pPart%Id,pPart%antiparticle)
         if ( pPart%antiparticle ) factor = -factor

         totalDensity(I1,I2,I3)=totalDensity(I1,I2,I3)+SW

      else

         factor = 1.

         if ( .not.pPart%antiparticle ) then
            totalDensity(I1,I2,I3)=totalDensity(I1,I2,I3)+SW
         else
            totalDensity(I1,I2,I3)=totalDensity(I1,I2,I3)-SW
         end if

      end if

      densityField(I1,I2,I3)%baryon = densityField(I1,I2,I3)%baryon &
           + vec*factor

      if ( pPart%ID.eq.nucleon ) then
         if ( pPart%charge.eq.0 ) then
            ! Evaluate neutron density
            densityField(I1,I2,I3)%neutron = densityField(I1,I2,I3)%neutron &
                 + vec*factor
         else  ! if ( abs(pPart%charge).eq.1 ) then
            ! Evaluate proton density
            densityField(I1,I2,I3)%proton = densityField(I1,I2,I3)%proton &
                 & + vec*factor
         end if
      end if

    end subroutine addTodensityfield


    !**************************************************************************
    !****f* updateDensity/initDensityWeights
    ! NAME
    ! subroutine initDensityWeights
    ! PURPOSE
    ! * Initializes weights which are used to evaluate the densities at the gridpoints.
    ! * Each particle has some coordinate off the gridpoints. Therefore it is needed to
    !   define weights which tell how much a particle at position r contributes to the
    !   ith grid point. The particle are smeared with a gaussian distribution.
    ! * The width is chosen such that it is equal to the maximum of the gridspacings.
    ! INPUTS
    ! * type(particle) Parts(:,:)
    ! USES
    ! * idTable, particleDefinition
    !**************************************************************************
    subroutine initDensityWeights
      use output
      use Callstack, only: TRACEBACK

      !grid spacing in small grid:
      real, dimension(1:3) :: smallGridSpacing


      real, save, dimension(1:3) :: smearingWidth=0.5  ! Gaussian width, used to smear the density
      real, save, dimension(1:3) :: smearingCutOff=5.  ! (cutoff for smearing)**2

      real,dimension(1:3) :: rLarge        ! Point on large grid
      real,dimension(1:3) :: rSmall        ! Point on small grid
      integer :: indexSmall1,indexSmall2,indexSmall3
      integer :: indexLarge1,indexLarge2,indexLarge3
      integer :: small,large
      real,dimension(1:3) :: rSquare        ! (distance of points on large and small grids)**2
      real :: norm             ! Normalization
      integer :: i
      logical,save :: firstTime=.true.

      call Write_InitStatus("Smearing weights for density",0)

      if (.not.firstTime) call TRACEBACK("WARNING: Allocating array again !!!",-1)
      allocate(smearingWeights(1:(2*smallerGridPoints+1)**3,1:(2*numberLargepoints+1)**3))

      ! Set smearing width in each direction equal to the grid spacing in this direction:

      if (setnewsmearing) then
         smearingWidth=gridSpacing*newsmearing
         smallGridSpacing=GridSpacing*newsmearing/(2.*float(smallerGridPoints)+1.)
      else
         smearingWidth=gridSpacing
         smallGridSpacing=GridSpacing/(2.*float(smallerGridPoints)+1.)
      end if


      ! Check that cut off is small than the maximal distance of a point on which I am smearing:
      do i = 1,3
         do
            if ( (float(numberlargePoints)+0.5)*gridSpacing(i).lt.sqrt(smearingCutOff(i)) ) then
!!$            If(firstTime) then
!!$               write (*,'("     ************************************************************")')
!!$               write (*,'("     Error! Smearing cut off too large!!! Enlarge numberlargePoints or decrease cut off!")')
!!$               write (*,'("     Setting smearing cut off to better value... ")')
!!$               write (*,'("     ************************************************************")')
!!$               firsttime=.false.
!!$            end if
               smearingCutOff(i)=smearingCutOff(i)*0.9
            else
               exit
            end if
         end do
      end do

      write(*,'(" Gridspacing in XYZ: ",3F8.5)') gridspacing(1:3)
      write(*,'(" Cut offs for smearing in XYZ: ",3F8.5," fm")') &
           sqrt(smearingCutOff)
      write(*,'(" Smearing with the function e^(-x^2/2/",F4.2,"-y^2/2/",F4.2,"-z^2/2/",F4.2,")")') &
           smearingWidth**2
      write(*,'(" Minimal distance of not-considered gridpoint: ",3F8.5)') &
           (float(numberlargePoints)+0.5)*gridSpacing
      write(*,*) "points in real  grid: ",numberLargePoints
      write(*,*) "points in small grid: ",smallerGridPoints


      smearingWeights = 0.0

      !Loop over points on small grid
      do indexSmall1=-smallerGridPoints,smallerGridPoints
         do indexSmall2=-smallerGridPoints,smallerGridPoints
            do indexSmall3=-smallerGridPoints,smallerGridPoints
               small=1+(smallerGridPoints+indexSmall3)& !translate 3dim coordinate in 3 letter number
                    & +(smallerGridPoints+indexSmall2)*(2*smallerGridPoints+1) &
                    & +(smallerGridPoints+indexSmall1)*(2*smallerGridPoints+1)**2
               rSmall=(/indexSmall1,indexSmall2,indexSmall3/)*smallGridSpacing

               Norm=0.
               large=0

               !Loop over points in the real grid:
               do indexLarge1=-numberLargePoints,numberLargePoints
                  do indexLarge2=-numberLargePoints,numberLargePoints
                     do indexLarge3=-numberLargePoints,numberLargePoints
                        large=large+1
                        rLarge=(/indexLarge1,indexLarge2,indexLarge3/)*GridSpacing

                        ! (Distance between points on small and real grid)**2:
                        rSquare(:) = (rLarge(:)-rSmall(:))**2

                        if (all(rSquare < SmearingCutoff)) then
                           smearingWeights(small,large)= &
                                exp(-sum(rSquare(:)/(2.*smearingWidth(:)**2)))
                        end if
                        Norm=Norm+smearingWeights(small,large)    !Normalization
                     end do
                  end do
               end do

               !Normalizing the weigths:
               smearingWeights(small,:) = smearingWeights(small,:)/Norm

            end do
         end do
      end do

      ! overall normalization factor:
      smearingWeights = smearingWeights &
           /(gridSpacing(1)*gridSpacing(2)*gridSpacing(3)) &
           /size(Parts(:,1))

      call Write_InitStatus("Smearing weights for density",1)
      firstTime = .false.

    end subroutine initDensityWeights

    !**************************************************************************
    !****f* updateDensity/FillStaticDensity
    ! NAME
    ! subroutine FillStaticDensity
    ! PURPOSE
    ! In the case of static density, this fills the fields densityField and
    ! totalDensity with values of the density parametrization. Needed in
    ! the case of RMF calculations and also for Coulomb potential.
    !**************************************************************************
    subroutine FillStaticDensity
      use densityStatic
      use nucleus, only: getTarget
      use output, only: Write_InitStatus

      real,dimension(1:3)    :: r         !position where density should be calculated

      integer :: Index1,Index2,Index3

      if (densitySwitch.ne.2) return
      call Write_InitStatus("FillStaticDensity",0)
      do Index1 = -gridPoints(1),gridPoints(1)
         r(1) = Index1*gridSpacing(1)
         do Index2 = -gridPoints(2),gridPoints(2)
             r(2) = Index2*gridSpacing(2)
            do Index3 = -gridPoints(3),gridPoints(3)
               r(3) = Index3*gridSpacing(3)
               densityField(Index1,Index2,Index3) = staticDensity(r,getTarget())
               totalDensity(Index1,Index2,Index3) = densityField(Index1,Index2,Index3)%Baryon(0)
            end do
         end do
      end do
      call Write_InitStatus("FillStaticDensity",1)

    end subroutine FillStaticDensity

  end subroutine updateDensity


  !****************************************************************************
  !****s* densitymodule/storeFields
  ! NAME
  ! subroutine storeFields
  ! PURPOSE
  ! Store the old values of the omega, rho and density fields.
  !****************************************************************************
  subroutine storeFields
    use output, only: DoPr
    integer :: k

    if (DoPr(2)) write(*,*) 'Store the old values of the omega, rho and density fields'

    omegaField_old = omegaField
    rhoField_old   = rhoField
    do k=1,3
      baryonCurrent_old(:,:,:,k) = densityField(:,:,:)%baryon(k)
    end do

  end subroutine storeFields


  !****************************************************************************
  !****s* densitymodule/updateRMF
  ! NAME
  ! subroutine updateRMF(Parts)
  ! PURPOSE
  ! Updates the sigmaField which is needed when the propagation
  ! with the RMF is done.
  ! Updates also the baryon velocities
  ! which are needed in the subsequent updating of the baryon 4-current
  ! by the subroutine updateDensity.
  ! INPUTS
  ! * type(particle), dimension(:,:) :: Parts
  ! USES
  ! idTable, particleDefinition,RMF,ADI,constants
  ! NOTES
  ! The scalarDensity is computed as well.
  !****************************************************************************
  subroutine updateRMF(Parts)

  use idTable
  use particleDefinition
  use RMF
  use ADI
  use constants, only: pi,hbarc,mN
  use output, only: DoPr

  type(particle), dimension(:,:) :: Parts

  integer :: Index1,Index2,Index3
  integer :: i,j,k,k_max
  integer :: small,large
  integer, dimension(1:3) :: posOrig
  integer, dimension(1:3) :: indexSmall
  integer :: niter

  real, dimension(1:3) :: rpos
  real, dimension(0:3) :: momentum
  real :: rho_lrf, pstar2, mstar, estar, pf_n, pf_p
  real :: funMax, fun, derfun, factor
  real :: eta_x, eta_y, tmp, lhs, error

  logical :: flagit
  logical, parameter :: debug = .false.
  logical, save :: flagini=.true.

  if (.not.flagini .and. densitySwitch.ne.1) then

     ! In the case of static density only
     ! update particle energies (E^*'s) and velocities:
     do j=1,Size(Parts,dim=1)
        do i=1,Size(Parts,dim=2)

          if ( Parts(j,i)%id == 0 ) then
             cycle
          else if ( Parts(j,i)%id < 0 ) then
             exit
          end if

          call energyDeterminationRMF( Parts(j,i) )

          Parts(j,i)%velocity(1:3) = Parts(j,i)%momentum(1:3) &
                                    & / Parts(j,i)%momentum(0)

        end do
     end do

     return

  end if

  if (DoPr(2)) write(*,*) 'Updating Relativistic Mean Fields'

  if (flagini) then

    ! Initialize sigmaField:

    do Index1 = -gridPoints(1),gridPoints(1)
      do Index2 = -gridPoints(2),gridPoints(2)
        do Index3 = -gridPoints(3),gridPoints(3)
          rho_lrf = densityField(Index1,Index2,Index3)%baryon(0)**2 &
                 &- dot_product(densityField(Index1,Index2,Index3)%baryon(1:3),&
                               &densityField(Index1,Index2,Index3)%baryon(1:3))
          rho_lrf = sqrt(max(0.,rho_lrf))
          sigmaField(Index1,Index2,Index3) = -fshift(rho_lrf)/g_sigma
        end do
      end do
    end do

    write(*,*) ' sigmaField initialized'

    do k = 0,3
      omegaField(:,:,:,k) = a_6/g_omega*densityField(:,:,:)%baryon(k)
    end do

    write(*,*) ' omegaField initialized'

    if (g_rho /= 0.0) then
       do k = 0,3
          rhoField(:,:,:,k) = &
               & a_7/g_rho*(densityField(:,:,:)%proton(k)-densityField(:,:,:)%neutron(k))
       end do
       write(*,*) ' rhoField initialized'
    else
       rhoField(:,:,:,:) = 0.0
    end if

    flagini=.false.

  end if


  if ((debug).and.(DoPR(2))) write(*,*) ' Iterations to find the sigma field:'

  flagit = .true.
  niter = 0
  Loop_over_iterations : do while(flagit)

       niter = niter + 1

       if (densitySwitch.eq.1) then

           scalarDensity(:,:,:) = 0.
           d_scalarDensity_dsigma(:,:,:) = 0.

           Loop_over_ensembles_1 : do j=1,Size(Parts,dim=1)
             Loop_over_particles_1 : do i=1,Size(Parts,dim=2)

               if ( Parts(j,i)%id == 0 ) then
                  cycle Loop_over_particles_1
               else if ( Parts(j,i)%id < 0 ) then
                  exit Loop_over_particles_1
               end if

               if ( Parts(j,i)%id >= pion) cycle  Loop_over_particles_1 ! Only baryons are accounted for presently

               ! Modification factor for the coupling constants:
               factor = ModificationFactor(Parts(j,i)%Id,Parts(j,i)%antiparticle)

               ! position in large grid:
               rpos = Parts(j,i)%position(1:3)/gridSpacing
               posOrig=NINT(rpos)

               if (      abs(posOrig(1)).gt.gridPoints(1) &
                  &.or. abs(posOrig(2)).gt.gridPoints(2) &
                  &.or. abs(posOrig(3)).gt.gridPoints(3) ) cycle Loop_over_particles_1

               pstar2 = dot_product(Parts(j,i)%momentum(1:3),Parts(j,i)%momentum(1:3))

               ! position in small grid:
               indexSmall=NINT((rpos-float(posOrig))*(2.*float(smallerGridPoints)+0.999999))

               ! Test for errors:
               if ((abs(indexSmall(1)).gt.smallerGridPoints)&
                    &.or.(abs(indexSmall(2)).gt.smallerGridPoints)&
                    &.or.(abs(indexSmall(3)).gt.smallerGridPoints)) then
                  write(*,*) 'Problem in updateRMF, module density.f90'
                  write(*,*) IndexSmall, 'too big, this must happen !!!'
                  stop
               end if

               small=1+(smallerGridPoints+indexSmall(3))                               &
                    & +(smallerGridPoints+indexSmall(2))*(2*smallerGridPoints+1)       &
                    & +(smallerGridPoints+indexSmall(1))*(2*smallerGridPoints+1)**2

               large=0

               ! Smearing particle over points in Neighborhood:
               do Index1=posOrig(1)-numberLargePoints,posOrig(1)+numberlargePoints
                  do Index2=posOrig(2)-numberLargePoints,posOrig(2)+numberlargePoints
                     do Index3=posOrig(3)-numberLargePoints,posOrig(3)+numberlargePoints

                        large=large+1

                        if (       abs(Index1).le.gridPoints(1) &
                          & .and. abs(Index2).le.gridPoints(2) &
                          & .and. abs(Index3).le.gridPoints(3) &
                          & .and. smearingWeights(small,large).gt.0. ) then

                           mstar = Parts(j,i)%mass + factor*g_sigma*sigmaField(Index1,Index2,Index3)
                           estar = sqrt( mstar**2 + pstar2 )

                           scalarDensity(Index1,Index2,Index3) = scalarDensity(Index1,Index2,Index3) &
                             & + mstar/estar*smearingWeights(small,large) * factor
                           d_scalarDensity_dsigma(Index1,Index2,Index3) &
                             & = d_scalarDensity_dsigma(Index1,Index2,Index3) &
                             & + g_sigma*pstar2/estar**3*smearingWeights(small,large) * factor**2

                        end if

                     end do
                  end do
               end do

             end do Loop_over_particles_1
           end do Loop_over_ensembles_1

       else ! Static density profiles are used:

           do Index1 = -gridPoints(1),gridPoints(1)
             do Index2 = -gridPoints(2),gridPoints(2)
               do Index3 = -gridPoints(3),gridPoints(3)

                 mstar = mN + g_sigma*sigmaField(Index1,Index2,Index3)
                 pf_n=(3.*pi**2*densityField(Index1,Index2,Index3)%neutron(0))**0.333333*hbarc
                 pf_p=(3.*pi**2*densityField(Index1,Index2,Index3)%proton(0))**0.333333*hbarc

                 scalarDensity(Index1,Index2,Index3) &
                       & = densityField(Index1,Index2,Index3)%neutron(0)*f(mstar/pf_n) &
                       & + densityField(Index1,Index2,Index3)%proton(0)*f(mstar/pf_p)

                 d_scalarDensity_dsigma(Index1,Index2,Index3) &
                       & = densityField(Index1,Index2,Index3)%neutron(0)*g_sigma/pf_n*fprime(mstar/pf_n) &
                       & + densityField(Index1,Index2,Index3)%proton(0)*g_sigma/pf_p*fprime(mstar/pf_p)


               end do
             end do
           end do

       end if

       funMax = 0.

       if (grad_flag) then
          tmp = -(gridSpacing(3)*m_sigma/hbarc)**2/g_sigma
          eta_x = (gridSpacing(3)/gridSpacing(1))**2
          eta_y = (gridSpacing(3)/gridSpacing(2))**2
       end if

       do Index1 = -gridPoints(1),gridPoints(1)
         do Index2 = -gridPoints(2),gridPoints(2)
           do Index3 = -gridPoints(3),gridPoints(3)

             fun = g_sigma*sigmaField(Index1,Index2,Index3) &
                 & + a_1*scalarDensity(Index1,Index2,Index3) &
                 & + a_2*sigmaField(Index1,Index2,Index3)**2 &
                 & + a_3*sigmaField(Index1,Index2,Index3)**3   ! function we want to have = 0.

             derfun = g_sigma*( 1. + a_4*sigmaField(Index1,Index2,Index3) &
                             &     + a_5*sigmaField(Index1,Index2,Index3)**2 ) &
                      & +a_1*d_scalarDensity_dsigma(Index1,Index2,Index3) ! d fun / d sigma

             if (.not.grad_flag) then ! local fields:

                if ( abs(fun) > funMax ) funMax = abs(fun)

                if ( derfun /= 0. ) &
                  & sigmaField(Index1,Index2,Index3) = sigmaField(Index1,Index2,Index3) &
                                                & - fun/derfun

             else ! gradients included:

                K_rmf(Index1,Index2,Index3) = tmp*( fun - derfun* sigmaField(Index1,Index2,Index3) )

                Diag_rmf(Index1,Index2,Index3) = -tmp*derfun/3.

                if ( abs(Index1).lt.gridPoints(1) .and. &
                   &abs(Index2).lt.gridPoints(2) .and. &
                   &abs(Index3).lt.gridPoints(3) ) then

                   lhs = ( -eta_x * (  sigmaField(Index1-1,Index2,Index3) &
                                   &+ sigmaField(Index1+1,Index2,Index3) ) &
                         & -eta_y * ( sigmaField(Index1,Index2-1,Index3) &
                                   &+ sigmaField(Index1,Index2+1,Index3) ) &
                         & -sigmaField(Index1,Index2,Index3-1) &
                         & -sigmaField(Index1,Index2,Index3+1) &
                         & + 2.*(eta_x+eta_y+1.) * sigmaField(Index1,Index2,Index3) ) / tmp

                   error = abs(fun - lhs)

                   if ( error .gt. funMax ) funMax = error

                end if

             end if

           end do
         end do
       end do

       if ((debug).and.(DoPR(2))) write(*,*) ' In updateRMF, niter, funMax: ', niter, funMax

       if (grad_flag) &
         & call ADI_solve_Douglas(sigmaField, K_rmf, Diag_rmf, eta_x, eta_y)

       if ( funMax <= 1.e-04 ) then
         flagit = .false.
       else if ( niter == 50 ) then
         write(*,*) ' In updateRMF: bad convergence after 50 iterations:',&
                   & funMax
         stop
       end if

  end do Loop_over_iterations

  ! Update particle energies (E^*'s) and velocities:
  Loop_over_ensembles_2 : do j=1,Size(Parts,dim=1)
    Loop_over_particles_2 : do i=1,Size(Parts,dim=2)

       if ( Parts(j,i)%id == 0 ) then
          cycle Loop_over_particles_2
       else if ( Parts(j,i)%id < 0 ) then
          exit Loop_over_particles_2
       end if

       call energyDeterminationRMF( Parts(j,i) )

       Parts(j,i)%velocity(1:3) = Parts(j,i)%momentum(1:3) &
                                 & / Parts(j,i)%momentum(0)

    end do Loop_over_particles_2
  end do Loop_over_ensembles_2

  call updateDensity(Parts)

  if (lorentz_flag) then
    k_max=3    ! With space components of the omega- and rho- fields
  else
    k_max=0    ! W/o space components of the omega- and rho- fields
  end if

  !-------------------------------------------------------------------------
  if ((debug).and.(DoPR(2))) write(*,*) ' Iterations to find the omega field:'
  !-------------------------------------------------------------------------
  do k = 0,k_max

     if (.not.grad_flag) then ! local fields:

        omegaField(:,:,:,k) = a_6/g_omega*densityField(:,:,:)%baryon(k)

     else ! gradients included:

        K_rmf(:,:,:) = gridSpacing(3)**2*g_omega*hbarc*densityField(:,:,:)%baryon(k)

        if (k.eq.0) Diag_rmf(:,:,:) = (gridSpacing(3)*m_omega/hbarc)**2 / 3.

        U_rmf(:,:,:) = omegaField(:,:,:,k)

        call ADI_solve_Douglas(U_rmf, K_rmf, Diag_rmf, eta_x, eta_y)

        omegaField(:,:,:,k) = U_rmf(:,:,:)

     end if

  end do

  if (g_rho /= 0.0) then
     !----------------------------------------------------------------------
     if ((debug).and.(DoPR(2))) write(*,*) ' Iterations to find the rho field:'
     !----------------------------------------------------------------------

     do k = 0,k_max

        if (.not.grad_flag) then ! local fields:

           rhoField(:,:,:,k) = &
                & a_7/g_rho*(densityField(:,:,:)%proton(k)-densityField(:,:,:)%neutron(k))

        else ! gradients included:

           K_rmf(:,:,:) = gridSpacing(3)**2*g_rho*hbarc*( densityField(:,:,:)%proton(k)&
                                                            & -densityField(:,:,:)%neutron(k) )

           if (k.eq.0) Diag_rmf(:,:,:) = (gridSpacing(3)*m_rho/hbarc)**2 / 3.

           U_rmf(:,:,:) = rhoField(:,:,:,k)

           call ADI_solve_Douglas(U_rmf, K_rmf, Diag_rmf, eta_x, eta_y)

           rhoField(:,:,:,k) = U_rmf(:,:,:)

        end if

     end do

  else
     do k = 0,k_max
        rhoField(:,:,:,k) = 0.0
     end do
  end if


  if ( fourMomDen_flag ) then ! Compute four-momentum density field:

       fourMomentumDensity(:,:,:,:)=0.

       Loop_over_ensembles_3 : do j=1,Size(Parts,dim=1)
         Loop_over_particles_3 : do i=1,Size(Parts,dim=2)

           if ( Parts(j,i)%id == 0 ) then
              cycle Loop_over_particles_3
           else if ( Parts(j,i)%id < 0 ) then
              exit Loop_over_particles_3
           end if

           ! position in large grid:
           rpos = Parts(j,i)%position(1:3)/gridSpacing
           posOrig=NINT(rpos)

           if (any(abs(posOrig).gt.gridPoints)) cycle Loop_over_particles_3

           ! position in small grid:
           indexSmall=NINT((rpos-float(posOrig))*(2.*float(smallerGridPoints)+0.999999))


           small=1+(smallerGridPoints+indexSmall(3))                               &
                & +(smallerGridPoints+indexSmall(2))*(2*smallerGridPoints+1)       &
                & +(smallerGridPoints+indexSmall(1))*(2*smallerGridPoints+1)**2

           large=0

           ! Smearing particle over points in Neighborhood:
           do Index1=posOrig(1)-numberLargePoints,posOrig(1)+numberlargePoints
              do Index2=posOrig(2)-numberLargePoints,posOrig(2)+numberlargePoints
                 do Index3=posOrig(3)-numberLargePoints,posOrig(3)+numberlargePoints

                    large=large+1

                    if (       abs(Index1).le.gridPoints(1) &
                      & .and. abs(Index2).le.gridPoints(2) &
                      & .and. abs(Index3).le.gridPoints(3) &
                      & .and. smearingWeights(small,large).gt.0. ) then

                       call true4Momentum_RMF(Parts(j,i),momentum)
                       fourMomentumDensity(Index1,Index2,Index3,0:3)&
                          &=fourMomentumDensity(Index1,Index2,Index3,0:3)&
                          &+momentum(0:3)*smearingWeights(small,large)

                    end if

                 end do
              end do
           end do

         end do Loop_over_particles_3
       end do Loop_over_ensembles_3

    end if

  end subroutine updateRMF



  !****************************************************************************
  !****s* densitymodule/energyDeterminationRMF
  ! NAME
  ! subroutine energyDeterminationRMF(Parts)
  ! PURPOSE
  ! This subroutine determines the one-particle energy E^*,
  ! which is the zeroth component of the kinetic four-momentum
  ! in the frame, where the space components of the kinetic four-momentum
  ! are given.
  ! INPUTS
  ! * type(particle),intent(inOut) :: Parts              ! Particle whose energy should be calculated.
  ! NOTES
  ! Should be used in RMF-mode. Please, notice, that not the full single-particle energy
  ! is computed here. The vector field contribution is missed in E^*.
  !****************************************************************************
  subroutine energyDeterminationRMF(Parts)

    use particleDefinition
    use idTable, only: kaon, kaonBar, isBaryon
    use RMF, only: ModificationFactor

    type(particle),intent(inOut) :: Parts ! particle whose energy should be calculated

    real, dimension(1:3) :: rpos
    real    :: pstar2, mstar, factor
    integer :: Index1, Index2, Index3


    ! Check input:
    if ( Parts%Id <= 0 ) then
      write(*,*) ' In energyDeterminationRMF: wrong input particle', Parts%Id
      stop
    end if

    factor = ModificationFactor(Parts%Id,Parts%antiparticle)

    pstar2 = dot_product(Parts%momentum(1:3),Parts%momentum(1:3))

    !for mesons mean-field optionally only for kaons and antikaons

    if ( (isBaryon(Parts%Id) .or. Parts%Id.eq.Kaon .or. Parts%Id.eq.kaonBar) &
       & .and. factor.gt.0. ) then

       ! position in large grid:
       rpos = Parts%position(1:3)/gridSpacing
       Index1=NINT(rpos(1))
       Index2=NINT(rpos(2))
       Index3=NINT(rpos(3))

       if (       abs(Index1).le.gridPoints(1) &
         & .and. abs(Index2).le.gridPoints(2) &
         & .and. abs(Index3).le.gridPoints(3)  ) then

         mstar = DiracMass(Index1,Index2,Index3,Parts%mass,Parts%id,Parts%charge,Parts%antiparticle)

       else

         mstar = Parts%mass

       end if

       Parts%momentum(0) = sqrt( mstar**2 + pstar2 )

    else

       Parts%momentum(0) = sqrt( Parts%mass**2 + pstar2 )

    end if

  end subroutine energyDeterminationRMF


  !****************************************************************************
  !****s* densitymodule/Particle4Momentum_RMF
  ! NAME
  ! subroutine Particle4Momentum_RMF(teilchen,momentum)
  ! PURPOSE
  ! This subroutine determines the canonical four-momentum
  ! in computational frame (i.e. where mesonic mean fields are given).
  ! INPUTS
  ! * type(particle),intent(in) :: teilchen              ! Particle
  ! OUTPUT
  ! * real, dimension(0:3), intent(out) :: momentum ! canonical 4-momentum of the particle
  ! NOTES
  ! Should be used in RMF-mode. It is supposed, that the kinetic four-momentum of the particle
  ! is already determined by the subroutine energyDeterminationRMF earlier.
  ! Electromagnetic part is not included.
  !****************************************************************************
  subroutine Particle4Momentum_RMF(teilchen,momentum)

    use particleDefinition
    use idTable, only: kaon,kaonBar,isBaryon
    use RMF, only: ModificationFactor

    type(particle), intent(in) :: teilchen        ! Particle
    real, dimension(0:3), intent(out) :: momentum ! canonical 4-momentum of the particle

    integer :: I1, I2, I3, k
    real :: fact
    real, dimension(0:3) :: V ! canonical 4-momentum of the particle

    ! Check input:
    if ( teilchen%Id <= 0 ) then
      write(*,*) ' In Particle4Momentum_RMF: wrong input particle', teilchen%Id
      stop
    end if

    fact = ModificationFactor(teilchen%Id,teilchen%antiparticle)

    momentum = teilchen%momentum

    if ( (isBaryon(teilchen%Id) .or. teilchen%id.eq.Kaon .or. teilchen%id.eq.kaonBar) &
       & .and. fact > 0. ) then

       ! position in large grid:
       I1 = NINT( teilchen%position(1) / gridSpacing(1) )
       I2 = NINT( teilchen%position(2) / gridSpacing(2) )
       I3 = NINT( teilchen%position(3) / gridSpacing(3) )

       if (       abs(I1).le.gridPoints(1) &
         & .and. abs(I2).le.gridPoints(2) &
         & .and. abs(I3).le.gridPoints(3)  ) then

          do k=0,3
             V(k) = SelfEnergy_vector(i1,i2,i3,k,teilchen%id,teilchen%charge,teilchen%antiparticle)
          end do

          momentum(0:3) = momentum(0:3) + V(0:3)

       end if

    end if

  end subroutine Particle4Momentum_RMF

  !****************************************************************************
  !****s* densitymodule/true4Momentum_RMF
  ! NAME
  ! subroutine true4Momentum_RMF(teilchen,momentum,inside_grid_flag)
  ! PURPOSE
  ! This subroutine determines the "true" sigle-particle 4-momentum,
  ! i.e. the 4-momentum which is additive to produce the total
  ! 4-momentum of the system in the computational frame.
  ! INPUTS
  ! * type(particle),intent(in) :: teilchen         ! particle whose "true" 4-momentum should be calculated.
  ! OUTPUT
  ! * real, dimension(0:3), intent(out) :: momentum ! "true" 4-momentum of the particle
  ! * logical, optional, intent(out) :: inside_grid_flag ! .true. if the particle is inside grid
  !                                                    ! .false. otherwise
  ! NOTES
  ! * Should be used in RMF-mode. The input particle kinetic 4-momentum must be given in
  !   the computational frame (where the density field is defined).
  ! * Formula for the fieldenergy updated; gradient terms replaced the corresponding sources using
  !   the meson-field equations and partial integrations.
  !****************************************************************************
  subroutine true4Momentum_RMF(teilchen,momentum,inside_grid_flag)

    use constants, only: hbarc
    use particleDefinition
    use idTable, only: kaon,kaonBar,isMeson
    use RMF, only: ModificationFactor, g_2, g_3, g_omega, g_rho, g_sigma, lorentz_flag

    type(particle),       intent(in) :: teilchen
    real, dimension(0:3), intent(out) :: momentum
    logical, optional,    intent(out) :: inside_grid_flag

    integer :: I1, I2, I3
    real    :: fieldEnergy, fact, bar_plus_antibar_density, factor  ! , isofact

    real, dimension(0:3) :: BaryonCurrent, IsospinCurrent

    !------------------------------------------------------------------------------------------------------

    if (present(inside_grid_flag)) inside_grid_flag = .false.

    factor = ModificationFactor(teilchen%Id,teilchen%antiparticle)

    if ( teilchen%ID <= 0 ) then

       write(*,*) 'In true4Momentum_RMF: wrong input particle Id', &
                 & teilchen%ID
       stop

    else if ( isMeson(teilchen%Id) .and. &
      &      ( teilchen%Id.ne.Kaon.and.teilchen%Id.ne.kaonBar .or. &
      &        (teilchen%id.eq.Kaon.or.teilchen%Id.eq.kaonBar).and.factor.eq.0. ) ) then

       momentum(0:3) = teilchen%momentum(0:3)
       return

    end if

    !------------------------------------------------------------------------------------------------------

    I1 = NINT( teilchen%position(1) / gridSpacing(1) )
    I2 = NINT( teilchen%position(2) / gridSpacing(2) )
    I3 = NINT( teilchen%position(3) / gridSpacing(3) )

    InsideGrid : if (  abs(I1) < gridPoints(1) .and. &
       & abs(I2) < gridPoints(2) .and. &
       & abs(I3) < gridPoints(3) ) then

       if (present(inside_grid_flag)) inside_grid_flag = .true.

       ! Energy density of the fields:
       ! valid in RMF with and without gradient terms

       BaryonCurrent(:)  = densityField(I1,I2,I3)%baryon(:)
       IsospinCurrent(:) = densityField(I1,I2,I3)%proton(:)-densityField(I1,I2,I3)%neutron(:)

       if (lorentz_flag) then  ! With space components of the omega- & rho-fields:

          fieldEnergy = &
              & 0.5*g_omega*( &
              &                 dot_product( omegaField(I1,I2,I3,1:3),BaryonCurrent(1:3) )&
              &               - BaryonCurrent(0)*omegaField(I1,I2,I3,0) ) &
              & + 0.5*g_rho*( &
              &                 dot_product(IsospinCurrent(1:3),rhoField(I1,I2,I3,1:3))&
              &               - IsospinCurrent(0)*rhoField(I1,I2,I3,0) ) &
              & - 0.5*g_sigma*scalarDensity(I1,I2,I3)*sigmaField(I1,I2,I3) &
              & - ( g_2*sigmaField(I1,I2,I3)**3/6. &
              &    +g_3*sigmaField(I1,I2,I3)**4/4. ) / hbarc**3

       else ! W/o space components of the omega- & rho-fields:

         fieldEnergy = - 0.5*g_omega*BaryonCurrent(0)*omegaField(I1,I2,I3,0) &
              &        - 0.5*g_rho*IsospinCurrent(0)*rhoField(I1,I2,I3,0) &
              &        - 0.5*g_sigma*scalarDensity(I1,I2,I3)*sigmaField(I1,I2,I3) &
              & - ( g_2*sigmaField(I1,I2,I3)**3/6. &
              &    +g_3*sigmaField(I1,I2,I3)**4/4. ) / hbarc**3

       end if

       fact = ModificationFactor(teilchen%Id,teilchen%antiparticle)

!        if (teilchen%ID==nucleon) then
!           if (teilchen%charge==0) then
!              isofact=-1.
!           else
!              isofact=1.
!           endif
!        else
!           isofact=0. !isospin sector presently only only for protons and neutrons
!        endif

       if ( teilchen%antiparticle .or. teilchen%ID==kaonBar ) fact=-fact

       bar_plus_antibar_density = totalDensity(I1,I2,I3)

       if ( bar_plus_antibar_density > 1.e-06 ) then

!          momentum(0) = teilchen%momentum(0) + g_omega*omegaField(I1,I2,I3,0)*fact &
!                       &+ isofact*fact*g_rho*rhoField(I1,I2,I3,0) &
!                       &+ fieldEnergy / bar_plus_antibar_density
          momentum(0) = teilchen%momentum(0) + SelfEnergy_vector(i1,i2,i3,0,teilchen%Id,teilchen%charge,teilchen%antiparticle) &
               &      + fieldEnergy / bar_plus_antibar_density
       else

          momentum(0) = teilchen%momentum(0)

       end if

       if (lorentz_flag) then
!         momentum(1:3) = teilchen%momentum(1:3) + g_omega*omegaField(I1,I2,I3,1:3)*fact &
!                       &+ isofact*fact*g_rho*rhoField(I1,I2,I3,1:3)
          momentum(1:3) = teilchen%momentum(1:3) &
                        + SelfEnergy_vector(i1,i2,i3,(/1,2,3/),teilchen%Id,teilchen%charge,teilchen%antiparticle)
       else
         momentum(1:3) = teilchen%momentum(1:3)
       end if

    else

       momentum(0:3) = teilchen%momentum(0:3)

    end if InsideGrid

  end subroutine true4Momentum_RMF


  !****************************************************************************
  ! calculates the effective (Dirac) mass of a particle
  !
  real function DiracMass(i1,i2,i3,barMass,id,charge,antiFlag)
    use IdTable, only: kaon,kaonBar
    use RMF, only: kaonpot_flag

    integer, intent(in) :: i1,i2,i3,id,charge
    logical, intent(in) :: antiFlag
    real,    intent(in) :: barMass

    integer :: i
    real :: Masse, S
    real, dimension(0:3) :: V

    S = SelfEnergy_scalar(i1,i2,i3,id,antiFlag)

    if ( kaonpot_flag .and. (id==kaon .or. id==kaonBar) ) then
       do i=0,3
          V(i) = Selfenergy_vector(i1,i2,i3,i,id,charge,antiFlag)
       end do
       Masse = sqrt(barMass**2-S+dot_product(V,V))
    else
       Masse = barMass + S
    end if

    DiracMass = Masse

  end function DiracMass
  !****************************************************************************


  !****************************************************************************
  ! calculates the vector component of the RMF selfenergy of a particle
  ! NOTE: gv_kaon = 3/(8 * f_pi*^2)
  !
  elemental real function SelfEnergy_vector(i1,i2,i3,k,id,charge,antiFlag)
    use constants, only: hbarc
    use IdTable, only: nucleon,kaon,kaonBar
    use RMF, only: ModificationFactor, g_omega, g_rho, gv_kaon, kaonpot_flag

    integer, intent(in) :: i1,i2,i3,k,id,charge
    logical, intent(in) :: antiFlag

    real :: field, fac, isofac

    if ( kaonpot_flag .and.  (id == kaon .or. id==kaonBar) ) then
       fac = 1.
       if (id==kaonBar) fac = -1.
       field = densityField(i1,i2,i3)%baryon(k) * gv_kaon * fac * hbarc**3 !GeV
    else
       fac = ModificationFactor(id,antiFlag)
       if (antiFlag) fac = -fac
       if (ID==nucleon) then
          if (charge==0) then
             isofac = -1.
          else
             isofac = 1.
          end if
       else
          isofac = 0.
       end if
       field = omegaField(i1,i2,i3,k)*g_omega*fac + rhoField(i1,i2,i3,k)*g_rho*isofac
    end if

    SelfEnergy_vector = field


  end function SelfEnergy_vector
  !****************************************************************************



  !****************************************************************************
  ! calculates the scalar component of the RMF selfenergy of a particle
  ! NOTE: for kaons selfenergy_scalar includes the Sigma_KN term only!
  !       gs_kaon = Sigma_KN / f_pi^2
  !
  real function SelfEnergy_scalar(i1,i2,i3,id,antiFlag)
    use constants, only: hbarc
    use IdTable, only: kaon,kaonBar
    use RMF, only: ModificationFactor, g_sigma, gs_kaon, kaonpot_flag

    integer, intent(in) :: i1,i2,i3,id
    logical, intent(in) :: antiFlag

    real :: field, fac

    if ( kaonpot_flag .and. (id == kaon .or. id==kaonBar) ) then
       field = scalarDensity(i1,i2,i3) * gs_kaon * hbarc**3 !GeV**2 !!
    else
       fac = ModificationFactor(id,antiFlag)
       field = sigmaField(i1,i2,i3)*g_sigma*fac !GeV
    end if

    SelfEnergy_scalar = field

  end function SelfEnergy_scalar
  !****************************************************************************

  !****************************************************************************
  ! calculates the vector component of the RMF selfenergy of a particle
  ! from the previous time step; needed for time derivatives in propagation_RMF
  ! NOTE: gv_kaon = 3/(8 * f_pi*^2)
  !
  real function SelfEnergy_vector_old(i1,i2,i3,k,id,charge,antiFlag)
    use constants, only: hbarc
    use IdTable, only: nucleon,kaon,kaonBar
    use RMF, only: ModificationFactor, g_omega, g_rho, gv_kaon, kaonpot_flag

    integer, intent(in) :: i1,i2,i3,k,id,charge
    logical, intent(in) :: antiFlag

    real :: field, fac, isofac

    if ( kaonpot_flag .and. (id == kaon .or. id==kaonBar) ) then
       fac = 1.
       if (id==kaonBar) fac = -1.
       field = baryonCurrent_old(i1,i2,i3,k) * gv_kaon * fac * hbarc**3 !GeV**2 !!!
    else
       fac = ModificationFactor(id,antiFlag)
       if (antiFlag) fac = -fac
       if (ID==nucleon) then
          if (charge==0) then
             isofac = -1.
          else
             isofac = 1.
          end if
       else
          isofac = 0.
       end if
       field = omegaField_old(i1,i2,i3,k)*g_omega*fac + rhoField_old(i1,i2,i3,k)*g_rho*isofac
    end if

    SelfEnergy_vector_old = field


  end function SelfEnergy_vector_old
  !****************************************************************************



  !****************************************************************************
  !****f* densitymodule/get_densitySwitch
  ! NAME
  ! function get_densitySwitch()
  ! PURPOSE
  ! If not yet done, reads input from jobcard and then returns densitySwitch.
  !****************************************************************************
  integer function get_densitySwitch()
    if (initFlag) call init
    get_densitySwitch=densitySwitch
  end function get_densitySwitch


  subroutine set_densitySwitch(d)
    integer,intent(in) :: d
    densitySwitch = d
  end subroutine set_densitySwitch


  subroutine set_densityInput(pro,neu)
    real, intent(in) :: pro, neu
    densityinput_proton = pro
    densityinput_neutron = neu
  end subroutine


  !****************************************************************************
  !****f* densitymodule/get_realGridSpacing
  ! NAME
  ! function get_realGridSpacing()
  ! PURPOSE
  ! returns GridSize
  ! RESULT
  ! real, dimension(1:3) gridSize in units of fermi
  !****************************************************************************
  function get_realGridSpacing()
    real, dimension(1:3) :: get_realGridSpacing
    if (initFlag) call init
    get_realGridSpacing(1:3)=gridSpacing(1:3)
  end function get_realGridSpacing


  !****************************************************************************
  !****f* densitymodule/getGridSpacing
  ! NAME
  ! function getGridSpacing()
  ! PURPOSE
  ! Evaluates GridSize for Derivatives.
  ! RESULT
  ! real, dimension(1:3) gridSize in units of fermi
  !****************************************************************************
  function getGridSpacing()
    real, dimension(1:3) :: getGridSpacing

    if (initFlag) call init
    if (linearInterpolation) then
       getGridSpacing(1:3)=gridSpacing(1:3)/4.
    else
       getGridSpacing(1:3)=gridSpacing(1:3)
    end if
  end function getGridSpacing


  !****************************************************************************
  !****f* densitymodule/getGridPoints
  ! NAME
  ! function getGridPoints()
  ! PURPOSE
  ! returns  GridPoints
  ! RESULT
  ! integer, dimension(1:3) gridPoints in units of fermi
  !****************************************************************************
  function getGridPoints()
    integer, dimension(1:3) :: getGridPoints

    if (initFlag) call init
    getGridPoints(1:3)=gridPoints(1:3)
  end function getGridPoints


  !****************************************************************************
  !****f* densitymodule/GetGridIndex
  ! NAME
  ! logical function GetGridIndex(r,ind,add)
  ! PURPOSE
  ! Returns .false. if the point is outside the grid (then 'ind' is invalid).
  ! Otherwise 'ind' is the index of the corresponding coordinate.
  ! NOTES
  ! The first tests according the real grid size is necessary, since it
  ! may give integer overfloats, if one only checks the integer values.
  !****************************************************************************
  logical function GetGridIndex(r,ind,add)
    real, dimension(1:3),intent(in) :: r
    integer, dimension(1:3),intent(out) :: ind
    integer,intent(in) :: add

    GetGridIndex = .false.

    if (any(abs(r).gt.gridSize+gridSpacing)) return ! --> failure
    ind = nint(r/gridSpacing)
    if (any(abs(ind).gt.gridpoints+add)) return ! --> failure

    GetGridIndex = .true.
  end function GetGridIndex


  !****************************************************************************
  !****s* densitymodule/boostToLRF
  ! NAME
  ! subroutine boostToLRF (teilchen, switch, density)
  ! PURPOSE
  ! Boosts particle between "Local Rest Frame (LRF)" and "calculation frame (CF)".
  ! INPUTS
  ! * type(particle), intent(inout) :: teilchen  ! Particle which is to be boosted
  ! * integer, intent(in) :: switch              ! 1= boost from CF to LRF
  !                                              ! 2= boost from LRF to CF
  ! NOTES
  ! The LRF is the frame in which the baryon current vanishes.
  ! If the density is very small, then no boost takes place.
  !****************************************************************************
  subroutine boostToLRF (teilchen, switch, density_in)
    use ParticleDefinition
    use lorentzTrafo, only: lorentz

    type(particle), intent(inout) :: teilchen
    integer, intent(in) :: switch
    type(dichte), intent(in), optional :: density_in

    real, dimension(1:3) :: beta
    type(dichte) :: density

    if (present(density_in)) then
      density=density_in
    else
      density=densityAt(teilchen%position(1:3))
    end if

    if (density%baryon(0)<1E-8) return  ! do nothing if density is too small

    beta(1:3) = density%baryon(1:3)/density%baryon(0)  ! beta of LRF in Calculation frame

    ! Boost particle's momentum
    select case (switch)
    case (1) ! The particle is boosted to LRF out of CF
       call lorentz( beta,teilchen%momentum(0:3), 'density(1)')
    case (2) ! The particle is boosted to calculation frame out of LRF.
       call lorentz(-beta,teilchen%momentum(0:3), 'density(2)')
    case default
       write(*,*) 'Wrong value for switch in boostToLRF:', switch
       Stop
    end select

  end subroutine boostToLRF


  !****************************************************************************
  !****f* densitymodule/fermiMomentum_sym
  ! NAME
  ! real function fermiMomentum_sym(rho)
  ! PURPOSE
  ! Evaluate the fermi momentum for symmetric nuclear matter. Assumption: rho_p=rho_n !!!
  ! INPUTS
  ! * real, intent(in) :: rho ! density in fm^-3
  ! RESULT
  ! * Fermi momentum in GeV
  !****************************************************************************
  real function fermiMomentum_sym(rho)
    use constants, only: hbarc, pi

    real, intent(in) :: rho ! Density in fm^-3
    ! Take care of improper input by checking rho:
    if (rho.lt.-1E-20) then
       write(*,*) 'ERROR in fermiMomentum_sym: Rho less than zero!!!'
       stop 'density.f90: fermiMomentum_sym'
    else if (rho.lt.1E-20) then
       fermiMomentum_sym=0.
    else
       fermiMomentum_sym=(3./2.*pi**2*rho*hbarc**3)**(1./3.)
    end if
  end function fermiMomentum_sym
  !****************************************************************************
  !****f* densitymodule/fermiMomentum_noIsospin
  ! NAME
  ! real function fermiMomentum_noIsospin(rho)
  ! PURPOSE
  ! Evaluate the fermi momentum for a gas of fermions.
  ! INPUTS
  ! * real, intent(in) :: rho ! density in fm^-3 (ATTENTION: Do not use rho(nucleon) but rho(proton) or rho(neutron)
  ! RESULT
  ! * Fermi momentum in GeV
  !****************************************************************************
  real function fermiMomentum_noIsospin(rho)
    use constants, only: hbarc, pi

    real, intent(in) :: rho ! Density in fm^-3
    ! Take care of improper input by checking rho:
    if (rho.lt.-1E-20) then
       write(*,*) 'ERROR in fermiMomentum_sym: Rho less than zero!!!'
       stop 'density.f90: fermiMomentum_sym'
    else if (rho.lt.1E-20) then
       fermiMomentum_noIsospin=0.
    else
       fermiMomentum_noIsospin=(3.*pi**2*rho*hbarc**3)**(1./3.)
    end if
  end function fermiMomentum_noIsospin

  !****************************************************************************
  !****f* densitymodule/FermiMomAt
  ! NAME
  ! real function FermiMomAt(pos,charge)
  !
  ! PURPOSE
  ! calculate the fermi momentum at some position. If no charge is given
  ! it uses the averaged proton/neutron densities, otherwise it uses
  ! the corresponding isospin channel
  !
  ! INPUTS
  ! * real, dimension(1:3) :: pos -- the space coordinates
  ! * integer, OPTIONAL :: charge -- the isospin channel
  !****************************************************************************
  real function FermiMomAt(pos,charge)
    use nucleusDefinition, only: tNucleus
    use minkowski, only: abs4
    use nucleus, only: getTarget
    use constants, only: pi,hbarc
    use densityStatic, only: staticDensity

    real,dimension(1:3),intent(in)    :: pos   !position where fermi mom should be calculated
    integer, intent(in), optional :: charge

    real :: rho
    type(tNucleus),pointer :: nuc
    type(dichte) :: dens

    if (initFlag) call init
    select case (densitySwitch)
    case (2) !Static density of a resting target
       nuc => getTarget()
       dens = staticDensity(pos,nuc)
       if (nuc%ReAdjustForConstBinding) then
          dens%proton = dens%proton/nuc%facP
          dens%neutron = dens%neutron/nuc%facN
          dens%baryon = dens%proton + dens%neutron
       end if
    case default
       dens = densityAt(pos)
    end select

    if (present(charge)) then
       select case (charge)
       case (0)
          rho = abs4(dens%neutron)
       case (1)
          rho = abs4(dens%proton)
       end select
    else
       rho = abs4(dens%baryon)/2
    end if
    FermiMomAt = (3.*pi**2*rho)**(1./3.)*hbarc

  end function FermiMomAt
  !****************************************************************************

  !****************************************************************************
  !****s* densitymodule/writeDensityPlane
  ! NAME
  !
  !
  ! PURPOSE
  ! Write the density to a file as cuts according some planes
  !
  ! INPUTS
  ! * character*(*) :: fName  -- name of file to write to
  ! * integer :: iPlane -- selection of plane to write (see below)
  !
  ! possible values for iPlane:
  ! * 1: yz-plane (through origin)
  ! * 2: xz-plane (through origin)
  ! * 3: xy-plane (through origin)
  !
  !****************************************************************************
  subroutine writeDensityPlane(fName,iPlane)

    use Callstack, only: TRACEBACK

    character*(*),intent(in) :: fName
    integer, intent(in) :: iPlane

    integer, parameter :: iFile = 645
    integer, dimension(1:3) :: indexOrigin
    integer :: i1, i2
    real, dimension(1:3) :: r

    if (.not.getGridIndex( (/ 0.,0.,0. /), indexOrigin, 0)) then
       call TraceBack('failure: Origin not in Grid!')
    end if


    open(iFile, file=fName, status='unknown')
    rewind(iFile)

    select case (iPlane)
    case (1)
       do i1 = -gridPoints(2),gridPoints(2)
          do i2 = -gridPoints(3),gridPoints(3)
             r = (/indexOrigin(1),i1,i2/) * gridSpacing
             write(iFile, '(3f12.4,1P,12e13.5)') r,densityField(indexOrigin(1),i1,i2)
          end do
          write(iFile,*)
       end do

    case (2)
       do i1 = -gridPoints(1),gridPoints(1)
          do i2 = -gridPoints(3),gridPoints(3)
             r = (/i1,indexOrigin(2),i2/) * gridSpacing
             write(iFile, '(3f12.4,1P,12e13.5)') r,densityField(i1,indexOrigin(2),i2)
          end do
          write(iFile,*)
       end do

    case (3)
       do i1 = -gridPoints(1),gridPoints(1)
          do i2 = -gridPoints(2),gridPoints(2)
             r = (/i1,i2,indexOrigin(3)/) * gridSpacing
             write(iFile, '(3f12.4,1P,12e13.5)') r,densityField(i1,i2,indexOrigin(2))
          end do
          write(iFile,*)
       end do

    end select



    close(iFile)

  end subroutine writeDensityPlane

end module densitymodule
