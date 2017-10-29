!******************************************************************************
!****m* /mesonWidthMedium_tables
! NAME
! module mesonWidthMedium_tables
! PURPOSE
! Implements the routines for the medium width of the mesons
!******************************************************************************
module mesonWidthMedium_tables

  use IDTable

  implicit none
  private

  public:: get_inMediumWidth,minimalWidth,numSteps_rhoN,numSteps_rhoP,numSteps_absP,numSteps_mass,max_rhoP,max_rhoN
  public:: getMaxMes,getMinMes,getMaxAbsP,cleanup


  !****************************************************************************
  !****g* mesonWidthMedium_tables/maxMes
  ! SOURCE
  !
  integer,save  :: maxMes=1000
  !
  ! PURPOSE
  ! Read the data table up to a maximum meson ID. ONLY FOR TESTING!!!
  !****************************************************************************

  !****************************************************************************
  !****g* mesonWidthMedium_tables/minMes
  ! SOURCE
  !
  integer,save  :: minMes=-1000
  !
  ! PURPOSE
  ! Read the data table starting at this minimal meson ID. ONLY FOR TESTING!!!
  !****************************************************************************


  !****************************************************************************
  !****g* mesonWidthMedium_tables/minimalWidth
  ! SOURCE
  !
  real,parameter  :: minimalWidth = 1E-6
  !
  ! PURPOSE
  ! * The minimal width a particle can get. Must be finite for numerics.
  !****************************************************************************

  ! Grid parameters:
  integer, parameter:: numSteps_rhoN=4
  integer, parameter:: numSteps_rhoP=4
  integer, parameter:: numSteps_absP=75
  integer, parameter:: numSteps_mass=400
  real, parameter   :: max_rhoP=0.1*0.197**3 ! In GEV**3
  real, parameter   :: max_rhoN=0.1*0.197**3 ! In GEV**3
  real, save, dimension(pion:pion+nmes-1) :: delta_mass,min_mass,max_mass,delta_absP,max_absP
  real, save ::  delta_rhoN,delta_rhoP

  ! the TABLE !
  real(4), save, dimension(:,:,:,:,:) , Allocatable:: widthTable

  logical,save  :: readInput_flag=.true.,readTable_flag=.true.

contains

  !****************************************************************************
  !****s* mesonWidthMedium_tables/readInput
  ! NAME
  ! subroutine readInput
  ! PURPOSE
  ! Reads input in jobcard out of namelist "mesonWidthMedium_tables".
  !****************************************************************************
  subroutine readInput

    use output

    integer :: ios ! checks file behavior

    NAMELIST /mesonWidthMedium_tables/ minMes,maxMes

    call Write_ReadingInput('mesonWidthMedium_tables',0)
    rewind(5)
    read(5,nml=mesonWidthMedium_tables,IOSTAT=IOS)
    call Write_ReadingInput('mesonWidthMedium_tables',0,IOS)

    write(*,*) 'maxMes                     ?    ', maxMes
    write(*,*) 'minMes                     ?    ', minMes

    call Write_ReadingInput('mesonWidthMedium_tables',1)

  end subroutine readInput


  subroutine cleanUp
    if (allocated(widthTable)) deallocate(widthTable)
  end subroutine


  integer function getMaxMes()
    if (readInput_flag) then
      call readInput
      readInput_flag=.false.
    end if
    getMaxMes=maxMes
  end function getMaxMes


  integer function getMinMes()
    if (readInput_flag) then
      call readInput
      readInput_flag=.false.
    end if
    getMinMes=minMes
  end function getMinMes


  real function getMaxAbsP(ID)
    integer,intent(in):: ID
    if (readTable_flag) then
      call readTable
      readTable_flag=.false.
    end if
    getMaxAbsP = max_absP(ID)
  end function getMaxAbsP


  !****************************************************************************
  !****f* mesonWidthMedium_tables/get_inMediumWidth
  ! NAME
  ! function get_inMediumWidth(particleID,absP,mass,rhoN,rhoP) result(width)
  !
  ! PURPOSE
  ! * Returns the in-medium width of mesons according to  Gamma=Gamma_free+sigma*rho*v
  ! * An average over the Fermi sea is performed
  ! * Sigma is given by the actual full collision term, with the incoming particles in the vacuum
  ! * Simplification: An average over the charge of the meson is performed.
  !
  ! INPUTS
  ! * integer, intent(in)              :: particleID    -- ID of meson
  ! * real  , intent(in)               :: absP          -- absolute Momentum (in LRF)
  ! * real ,intent(in)                 :: mass          -- Mass of meson !  p_mu p^mu = mass of the meson (offshell)
  ! * real, intent (in)                :: rhoN,rhoP     -- proton and neutron density in fm^-3
  ! * logical, optional, intent(in)    :: doMassCorrect_in ! only for debugging: switch on/off mass correction
  !
  ! OUTPUT
  ! * real :: collisional width
  !
  !****************************************************************************
  function get_inMediumWidth(particleID,absP,massIN,med) result(width)   ! doMassCorrect_in
    use idTable, only: pion,nmes
    use mediumDefinition
!    use particleDefinition, only : particle
!     use potentialModule, only : potential_LRF
    use tabulation, only: interpolate4

    integer, intent(in)              :: particleID
    real  , intent(in)               :: absP
    real ,intent(in)                 :: massIN           ! p_mu p^mu = mass of the resonance (offshell)
    !real, intent (in)                :: rhoN_in,rhoP_in  ! In fm^-3
    type(medium), intent(in)         :: med
!     logical, optional, intent(in)    :: doMassCorrect_in ! only for debugging
    real :: width

    real :: mass !,p_null,pot
!     type(particle) :: part
!     logical :: doMassCorrect

    real(4), dimension(1:4) :: X
    logical,save :: first=.true.
    real(4), dimension(1:4),save :: minX, maxX, deltaX
    real :: rhoN,rhoP

!     doMassCorrect=.true.
!     if(present(doMassCorrect_in)) doMassCorrect=domassCorrect_in

    mass=massIN
!     if(doMassCorrect) then
!        ! Shift the mass according to the potential:
!        ! (1) Define p_0
!        p_null  =   sqrt(absP**2+mass**2)
!        ! (2) Get the potential
!        !    define a particle:
!        part%ID=particleID
!        part%momentum(0)=0   ! potential shall not depend on the momentum!!!
!        part%momentum(1:3)= (/absP,0.,0./)
!        part%mass=mass
!        part%perturbative=.true.
!        !    get the potential for that particle
!        pot= potential_LRF(part,med)
!        !    remove the potential
!        p_null  =   p_null-pot
!        mass    =   sqrt(max(1E-4,p_null**2-absP**2))
!     end if

    if ((particleID.lt.pion).or.(particleID.gt.pion+nmes-1)) then
       write(*,*) 'ERROR: particleID out of bounds in get_inMediumWidth',particleID, pion,pion+nmes-1
    end if

    ! convert densities from fm^(-3) to GeV^3
    rhoN=med%densityNeutron*0.197**3
    rhoP=med%densityProton*0.197**3

    if (readInput_flag) then
       call readinput
       readInput_flag=.false.
    end if

    if (readTable_flag) then
       ! READ IN THE DATA FILES:
       call readTable
       readTable_flag=.false.
    end if

    if (first) then
       ! Initialize arrays containing grid information
       minX = (/ 0., min_mass(particleID), 0., 0. /)
       maxX = (/ max_absP(particleID), max_mass(particleID), max_rhoN, max_rhoP /)
       deltaX = (/ delta_absP(particleID), delta_Mass(particleID), delta_rhoN, delta_rhoP /)
       first = .false.
    else
       ! Only update absP & mass parameters of the grid information
       minX(2)  = min_mass(particleID)
       maxX(1:2) = (/ max_absP(particleID), max_mass(particleID) /)
       deltaX(1:2) = (/ delta_absP(particleID), delta_Mass(particleID) /)
    end if

    x(1:4) = (/ absP,mass,rhoN,rhoP /)
    width = interpolate4(X,minX,maxX,deltaX,widthTable(particleId,:,:,:,:))

    width=max(minimalWidth,width)

  end function get_inMediumWidth


  !****************************************************************************
  !****s* mesonWidthMedium_tables/readTable
  ! NAME
  ! subroutine readTable
  !
  ! PURPOSE
  ! * Reads tabulated arrays from files buuinput/inMediumWidth/InMediumWidth."particleID".dat.bz2
  ! OUTPUT
  ! * Initialized arrays widhtTable, min_mass, max_mass, delta_mass, max_absP, delta_absP
  !****************************************************************************
  subroutine readTable
    use inputGeneral, only: path_To_Input
    use output, only: intToChar,Write_ReadingInput
    use bzip

!     character(20) :: format
    character(5) ::raute
    character(200) ::fileName
    integer :: particleID,index_absP,index_mass,index_rhoN,ios
    logical :: fileFailure
    integer :: numSteps_absP_in,numSteps_mass_in,numSteps_rhoN_in,numSteps_rhoP_in,num_MonteCarlo_Points_in,&
               min_charge_loop_in,max_charge_loop_in
    real :: max_absP_in,max_rhoP_in,max_rhoN_in,minMass_in,maxMass_in
    type(bzFile) :: f
    character(len=100) :: buf
    integer :: ll

    if (readInput_flag) then
       call readinput
       readInput_flag=.false.
    end if

    allocate(widthTable(max(minMes,pion):min(maxMes,pion+nmes-1),0:numSteps_absP,0:numSteps_mass, &
             0:numSteps_rhoN,0:numSteps_rhoP))
    delta_rhoN=max_rhoN/float(numSteps_rhoN)
    delta_rhoP=max_rhoP/float(numSteps_rhoP)

!     format='(' // intToChar(numSteps_rhoP+1) // 'E15.4)'

    call Write_ReadingInput("InMediumWidth (meson)",0)

    fileFailure=.false.
       idLoop: do particleId=max(minMes,pion),min(maxMes,pion+nmes-1)
          fileName=trim(path_to_Input)//'/inMediumWidth/InMediumWidth.'//trim(intToChar(particleID))//'.dat.bz2'
          f = bzOpenR(trim(fileName))
          ll = 0
          call bzReadLine(f,buf,ll)
          read(buf(1:ll),*,iostat=ios) raute, numSteps_absP_in, numSteps_mass_in , numSteps_rhoN_in,&
               & numSteps_rhoP_in,num_MonteCarlo_Points_in,min_charge_loop_in,max_charge_loop_in
          if (ios.ne.0) then
             write(*,*) 'Error in opening input file (2): ',trim(fileName)
             fileFailure=.true.
             exit idLoop
          end if
          ll = 0
          call bzReadLine(f,buf,ll)
          read(buf(1:ll),*,iostat=ios) raute, max_absP_in, max_rhoP_in , max_rhoN_in,minMass_in, maxMass_in
          if (ios.ne.0) then
             write(*,*) 'Error in opening input file (3): ',trim(fileName)
             fileFailure=.true.
             exit idLoop
          end if
          ! set mass boundaries
          min_mass(particleID)=minMass_in
          max_mass(particleID)=maxMass_in
          delta_mass(particleID)=(max_mass(particleID)-min_mass(particleID))/float(numSteps_mass)
          max_absP(particleID)=max_absP_in
          delta_absP(particleID)=max_absP(particleID)/float(numSteps_absP)
          write(*,'(A,I3,A,F8.6,A,F8.6,A,F8.6,A,F8.6,A,F8.6,A,I2,A,I2,A,I4)') 'ID=',particleID, &
               '   minM=',minMass_in,' maxM=',maxMass_in,' dM=',delta_mass(particleID), &
               '   maxP=',max_absP_in,' dP=',delta_absP(particleID), &
               '   minCh=',min_charge_loop_in,' maxCh=',max_charge_loop_in, &
               '   MC_points=',num_MonteCarlo_Points_in
          ! check boundaries
          if ((numSteps_absP.ne.numSteps_absP_in) &
               & .or.(numSteps_rhoP.ne.numSteps_rhoP_in) &
               & .or.(numSteps_rhoN.ne.numSteps_rhoN_in) &
               & .or.(numSteps_mass.ne.numSteps_mass_in) &
               & .or.(abs(max_rhoP_in-max_rhoP).gt.0.00001) &
               & .or.(abs(max_rhoN_in-max_rhoN).gt.0.00001) &
               & ) then
             write(*,*) 'Error in opening input file: ',trim(fileName)
             write(*,*) 'Grid points differ!'
             if (numSteps_absP.ne.numSteps_absP_in)  write(*,*) 'absP',numSteps_absP,numSteps_absP_in
             if (numSteps_rhoP.ne.numSteps_rhoP_in)  write(*,*) 'rhoP',numSteps_rhoP,numSteps_rhoP_in
             if (numSteps_rhoN.ne.numSteps_rhoN_in)  write(*,*) 'rhoN',numSteps_rhoN,numSteps_rhoN_in
             if (numSteps_mass.ne.numSteps_mass_in) write(*,*) 'mass',numSteps_mass,numSteps_mass_in
             if (abs(max_rhoP_in-max_rhoP).gt.0.00001) write(*,*) 'max_rhop', max_rhoP,max_rhoP_in
             if (abs(max_rhoN_in-max_rhoN).gt.0.00001) write(*,*) 'max_rhon', max_rhoN,max_rhoN_in
             fileFailure=.true.
             exit idLoop
          end if

          ll = 0
          do index_absP=0,numSteps_absP
             do index_mass=0,numSteps_mass
                do index_rhoN=0,numSteps_rhoN
                   call bzReadLine(f,buf,ll)
                   read(buf(1:ll),*,iostat=ios) widthTable(particleID,index_absP,index_mass,index_rhoN,0:numSteps_rhoP)
                   if (IOS.ne.0) then
                      write(*,*) 'Error in opening input file: ',fileName
                      write(*,*) 'Error', particleID,ios,index_absP,index_mass,index_rhoN
                      fileFailure=.true.
                      exit idLoop
                   end if
                end do
             end do
          end do
       end do idLoop
    call bzCloseR(f)
    if (fileFailure) then
       write(*,*) 'You first need to tabulate the IN-MEDIUM WIDTH table. STOP!!'
       stop
    end if
    call Write_ReadingInput("InMediumWidth (meson)",1)

  end subroutine readTable


end module mesonWidthMedium_tables
