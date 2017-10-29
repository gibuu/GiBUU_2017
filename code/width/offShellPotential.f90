!******************************************************************************
!****m* /offShellPotential
! NAME
! module offShellPotential
! PURPOSE
! This module calculates the offshell potential.
!******************************************************************************
module offShellPotential

  implicit none
  private

  !****************************************************************************
  !****g* offShellPotential/OffShell_debug
  ! PURPOSE
  ! Switch on or off whether the debug information shall be given.
  ! SOURCE
  logical, parameter :: OffShell_debug=.false.
  !****************************************************************************

  !****************************************************************************
  !****g* offShellPotential/useOffShellPotentialBaryons
  ! PURPOSE
  ! Switch on or off whether the offshellness should be used for baryons.
  ! SOURCE
  logical, save :: useOffShellPotentialBaryons=.false.
  !
  ! NOTES
  ! * must be set to "TRUE" if mediumSwitch_coll
  !   (see module BaryonWidthMedium) is .true.
  ! * if .true. then delta_T (see module inputGeneral) must be <=0.05
  !   AND delta_P (see module propagation)  must be <=0.002;
  !   AND delta_E (see module propagation)  must be <=0.002;
  !   slows down propagation by a factor of 10
  !****************************************************************************

  !****************************************************************************
  !****g* offShellPotential/useOffShellPotentialMesons
  ! PURPOSE
  ! Switch on or off whether the offshellness should be used for mesons.
  ! SOURCE
  logical, save :: useOffShellPotentialMesons=.false.
  !****************************************************************************

  !****************************************************************************
  !****g* offShellPotential/offshell_cutoff
  ! PURPOSE
  ! If offshellParameter is less than this value, then we treat
  ! it as zero -> getoffshellMass returns the pole mass!
  ! SOURCE
  real   , parameter :: offshell_cutoff = 1E-5
  !****************************************************************************

  !****************************************************************************
  !****g* offShellPotential/max_offshellparameter
  ! PURPOSE
  ! The maximal value for the offshell parameter. Note: empirical value!
  ! This only applies to baryons. For mesons we have no restrictions
  ! on the offshell parameter.
  ! SOURCE
  real, save :: max_offshellparameter=5.
  !****************************************************************************

  !****************************************************************************
  !****g* offShellPotential/extrapolateBaryonWidth
  ! PURPOSE
  ! Whether to extrapolate the baryon width below minimal mass or not.
  ! SOURCE
  logical, save :: extrapolateBaryonWidth=.true.
  !****************************************************************************


  !****************************************************************************
  !****g* offShellPotential/relativistic
  ! PURPOSE
  ! * false: Use non-rel. off-shell parameter x=Delta m/Gamma,
  !   which obeys Stefan Leupold's non-rel. EOM.
  ! * true: Use rel. off-shell parameter x=Delta m^2/Gamma,
  !   which obeys Cassing's rel. EOM.
  ! SOURCE
  logical, save :: relativistic = .false.
  !****************************************************************************


  logical, save :: initFlag    =.true.

  public :: get_useOffShellPotentialBaryons, get_useOffShellPotentialMesons, treatParticleOffShell, &
            getOffShellParameter, HamiltonFunc_offshell, get_offshell_debug, setOffShellParameter



  !****************************************************************************
  !****s* offShellPotential/setOffShellParameter
  ! NAME
  ! Interface setOffShellParameter
  ! PURPOSE
  ! This is an interface which can be called like a subroutine using
  ! "call setOffShellParameter(p,success)".
  !
  ! It calculates the offshell parameter for a single particle or a set of particles
  ! If it fails for one of them, then it returns .false., otherwise .true. . The particles should
  ! be properly initialized, only the offshellparameter should be left-over to initialize.
  !
  ! INPUTS
  ! * type(particle) , intent(inout), dimension (:) :: p
  ! or:
  ! * type(particle) , intent(inout) :: p
  !
  ! OUTPUT
  ! * logical, intent(out):: success -- .true. if the offshell parameter could be evaluated for all particles
  !****************************************************************************
  Interface setOffShellParameter
     Module Procedure setOffShellParameter_1dim, setOffShellParameter_single
  end Interface



contains

  !****************************************************************************
  !****s* offShellPotential/readInput
  ! NAME
  ! subroutine readInput
  ! PURPOSE
  ! Reads input in jobcard out of namelist "OffShellPotential".
  !****************************************************************************
  subroutine readInput
    use output
    use baryonWidthMedium, only: get_MediumSwitch_coll
    use mesonWidthMedium, only: get_MediumSwitchMesons
    use inputGeneral
    use eventtypes

    integer :: ios

    !**************************************************************************
    !****n* offShellPotential/OffShellPotential
    ! NAME
    ! NAMELIST /OffShellPotential/
    ! PURPOSE
    ! Includes the switches:
    ! * useOffShellPotentialBaryons
    ! * useOffShellPotentialMesons
    ! * extrapolateBaryonWidth
    ! * max_offshellparameter
    ! * relativistic
    !**************************************************************************
    NAMELIST /offShellPotential/ useOffShellPotentialBaryons, useOffShellPotentialMesons, &
                                 extrapolateBaryonWidth, max_offshellparameter, relativistic

    call Write_ReadingInput('offShellPotential',0)
    rewind(5)
    read(5,nml=offShellPotential,IOSTAT=IOS)
    call Write_ReadingInput('offShellPotential',0,IOS)

    if (.not.get_MediumSwitch_coll() .and. useOffShellPotentialBaryons) then
       write(*,*) 'MediumSwitch_coll of baryons is switched off, therefore we now switch off useOffShellPotentialBaryons'
       useOffShellPotentialBaryons = .false.
    else if (get_MediumSwitch_coll() .and. .not.useOffShellPotentialBaryons) then
       write(*,*) 'MediumSwitch_coll of baryons is switched on, therefore useOffShellPotentialBaryons must be set to TRUE!'
       stop
    end if

    if (.not.get_MediumSwitchMesons() .and. useOffShellPotentialMesons) then
       write(*,*) 'medium width of mesons is switched off, therefore we now switch off useOffShellPotentialMesons'
       useOffShellPotentialMesons = .false.
    else if (get_MediumSwitchMesons() .and. .not.useOffShellPotentialMesons) then
       write(*,*) 'medium width of mesons is switched on, therefore useOffShellPotentialMesons must be set!'
       useOffShellPotentialMesons = .false.
    end if

    if ((useOffShellPotentialBaryons.or.useOffShellPotentialMesons) .and. .not.LRF_equals_CALC_frame) then
       write(*,*) 'STOP: Off-shell potential works only in the "LRF=Calculation frame" assumption'
       write(*,*) 'The more general case is not yet implemented!!'
       stop
    end if

    write(*,*) 'use offShell Potential for the baryons?', useOffShellPotentialBaryons

    write(*,*) 'use offShell Potential for the mesons?', useOffShellPotentialMesons

    write(*,*) 'extrapolate baryon width below minimal mass?', extrapolateBaryonWidth

    write(*,*) 'max. offshell parameter for baryons: ', max_offshellparameter

    write(*,*) 'use relativistic off-shell parameter?', relativistic

    call Write_ReadingInput('offShellPotential',1)

  end subroutine readInput


  !****************************************************************************
  !****f* offShellPotential/get_useOffShellPotentialBaryons()
  ! NAME
  ! logical function get_useOffShellPotentialBaryons()
  ! PURPOSE
  ! returns the value of useOffShellPotentialBaryons
  !****************************************************************************
  logical function get_useOffShellPotentialBaryons()
    if (initFlag) then
       call readInput
       initFlag=.false.
    end if
    get_useOffShellPotentialBaryons=useOffShellPotentialBaryons
  end function get_useOffShellPotentialBaryons


  !****************************************************************************
  !****f* offShellPotential/get_useOffShellPotentialMesons()
  ! NAME
  ! logical function get_useOffShellPotentialMesons()
  ! PURPOSE
  ! returns the value of useOffShellPotentialMesons
  !****************************************************************************
  logical function get_useOffShellPotentialMesons()
    if (initFlag) then
       call readInput
       initFlag=.false.
    end if
    get_useOffShellPotentialMesons=useOffShellPotentialMesons
  end function get_useOffShellPotentialMesons


  !****************************************************************************
  !****f* offShellPotential/get_OffShell_debug()
  ! NAME
  ! logical function get_OffShell_debug()
  ! PURPOSE
  ! returns the value of OffShell_debug
  !****************************************************************************
  logical function get_OffShell_debug()
    if (initFlag) then
       call readInput
       initFlag=.false.
    end if
    get_OffShell_debug=OffShell_debug
  end function get_OffShell_debug


  !****************************************************************************
  !****f* offShellPotential/treatParticleOffShell(partID,partOffShellParameter)
  ! NAME
  ! logical function treatParticleOffShell(partID,partOffShellParameter)
  ! PURPOSE
  ! returns for a given particle, whether it is treated as
  ! offShell particle (e.g. in propagation).
  !****************************************************************************
  logical function treatParticleOffShell(partID,partOffShellParameter)
    use IdTable, only: isBaryon, isMeson
    integer, intent(in) :: partID
    real, intent(in) :: partoffShellparameter
    if (initFlag) then
       call readInput
       initFlag=.false.
    end if

    if (useOffShellPotentialBaryons.and.isBaryon(partId).and.abs(partOffshellparameter).gt.offshell_cutoff) then
       treatParticleOffShell=.true.
    else if (useOffShellPotentialMesons.and.isMeson(partId).and.abs(partOffshellparameter).gt.offshell_cutoff) then
       treatParticleOffShell=.true.
    else
       treatParticleOffShell=.false.
    end if

  end function treatParticleOffShell


  !****************************************************************************
  !****f* offShellPotential/HamiltonFunc_offshell
  ! NAME
  ! real function HamiltonFunc_offshell(part)
  ! PURPOSE
  ! determines Hamilton function for the offshell potential prescription
  !****************************************************************************
  real function HamiltonFunc_offshell(part_in,outOfBounds,massIsDetermined_in,full_offShell_in)
    use particleDefinition
    use potentialModule, only: potential_LRF,massdetermination
    use minkowski, only: abs4
    use output, only: writeParticle_debug
    use IdTable, only: isMeson
    use selfenergy_mesons, only: get_realPart
    use mediumDefinition
    use mediumModule, only: mediumAt
    !use callStack, only: traceBack

    type(particle),intent(in) :: part_in
    type(particle) :: part
    real :: mass,rp
    logical :: success
    logical :: outOFBounds ! .true. if the width table is out of bounds
    logical, optional :: massIsDetermined_in, full_offShell_in
    logical :: massIsDetermined,full_offShell
    type(medium) :: med

    if (present(massIsDetermined_in)) then
       massIsDetermined=massIsDetermined_in
    else
       massIsDetermined=.false.
    end if
    if (present(full_offshell_in)) then
       full_offshell=full_offshell_in
    else
       full_offshell=.true.
    end if

    if (initFlag) then
       call readInput
       initFlag=.false.
    end if

    part=part_in

    if (treatParticleOffShell(part%ID,part%OffShellParameter)) then
       if (full_offshell) then
          !Offshell mass includes off shell potential!
          if (.not.massIsDetermined) then
             call massDetermination(part,success=success)
             if (.not. success) then
                hamiltonFunc_offshell=99999.
                !call traceBack('problem: No success in massDet in HamiltonFunc',-1)
                call writeParticle_debug(part)
                return
             end if
          end if
          med = mediumAt (part%position)
          med%useMedium = .true. ! to avoid threshold effects
          if (isMeson(part%ID)) then
            rp = get_realPart (part%ID, abs4(part%momentum), med)
          else
            rp = 0.
          end if
          mass=getOffShellMass(part%ID,part%offshellparameter,part%momentum,part%mass,med,outOfBounds)
          hamiltonFunc_offshell= sqrt(Dot_Product(part%momentum(1:3),part%momentum(1:3))+mass**2 + rp) + potential_LRF(part)
       else
          outofBounds=.false.
          ! no offshell potential
          hamiltonFunc_offshell= sqrt(Dot_Product(part%momentum(1:3),part%momentum(1:3))+part%mass**2) + potential_LRF(part)
       end if
    else
       write(*,*) 'HAMILTONFUNC_offshell should not be called -> STOP'
       call offShellErrorMessage(part)
       STOP
    end if

  end function HamiltonFunc_offshell



  !****************************************************************************
  !****f* offShellPotential/setOffShellParameter_single
  ! NAME
  ! subroutine setOffShellParameter_single(p,success)
  !
  ! PURPOSE
  ! This subroutine calculates the offshell parameter for one given particle. If it fails,
  ! then it returns .false. otherwise .true.
  !
  ! INPUTS
  ! * type(particle),intent(inout)  :: p
  !
  ! OUTPUT
  ! * logical, intent(out):: success -- .true. if the offshell parameter could be evaluated.
  !****************************************************************************
  subroutine setOffShellParameter_single(p,success)
    use particleDefinition

    type(particle),intent(inout)  :: p
    logical, intent(out):: success

    p%offshellParameter=getOffShellParameter(p%ID,p%mass,p%momentum,p%position,success)

  end subroutine setOffShellParameter_single



  !****************************************************************************
  !****f* offShellPotential/setOffShellParameter_1dim
  ! NAME
  ! subroutine setOffShellParameter_1dim(p,success)
  !
  ! PURPOSE
  ! This subroutine calculates the offshell parameter for a given set of particles. If it fails for one
  ! of them, then it returns .false. otherwise .true.
  !
  ! INPUTS
  ! * type(particle) , intent(inout), dimension (:) :: p
  !
  ! OUTPUT
  ! * logical, intent(out):: success -- .true. if the offshell parameter could be evaluated.
  !****************************************************************************
  subroutine setOffShellParameter_1dim(p,success)
    use particleDefinition

    type(particle) , intent(inout), dimension (:) :: p
    logical, intent(out):: success
    logical :: flagOK
    integer :: i
    success=.true.
    do i= lbound(p,dim=1),ubound(p,dim=1)
       if (p(i)%ID.le.0) cycle
       call setOffShellParameter_single(p(i),flagOk)
       if (.not.flagOk) then
          success=.false.
       end if
    end do
  end subroutine setOffShellParameter_1dim

  !****************************************************************************
  !****f* offShellPotential/getOffShellParameter
  ! NAME
  ! real function getOffShellParameter
  !
  ! PURPOSE
  ! This function calculates the offshell parameter to be set into particle%offshellparameter.
  ! When calling this routine, make sure, that momentum(0) is set correctly!!!!!!!!
  !
  ! INPUTS
  ! * integer, intent(in) :: partID   -- ID of particle
  ! * real, intent(in) :: bareMass     -- bare mass= %mass
  ! * real, dimension(0:3), intent(in) :: momentum  -- 4-momentum at production point
  ! * real, dimension(1:3), intent(in) :: position  -- production point
  !
  ! OUTPUT
  ! * logical, intent(out):: success -- .true. if the offshell parameter could be evaluated.
  ! * real ::                getOffShellParameter
  !****************************************************************************
  real function getOffShellParameter(partID,bareMass,momentum,position,success)
    use IDTable, only: rho, omegaMeson, phi, isBaryon
    use particleProperties, only: hadron
    use mediumDefinition
    use hist2Df90
    use mediumModule, only: mediumAt
    use minkowski, only: abs4
    use mesonWidthMedium, only: WidthMesonMedium

    integer, intent(in) :: partID
    real, intent(in) :: bareMass
    real, dimension(0:3), intent(in) :: momentum
    real, dimension(1:3), intent(in) :: position
    logical, intent(out):: success

    real :: width, poleMass

    type(medium) :: mediumAtPosition

    !for debugging
    integer,save :: numcalls=0
    type(histogram2D),save :: histo_nucl
    type(histogram2D),save :: histo_delta,histo_delta_rho
    logical, save :: initHistFirst=.true.
    logical :: outofBounds
    integer, save :: counter_too_large=0

    if (initFlag) then
       call readInput
       initFlag=.false.
    end if

    numcalls=numcalls+1

    !set standard output
    getOffShellParameter=0.
    success=.false.

    if (partID.le.0) return

    mediumAtPosition = mediumAt(position)
    mediumAtPosition%useMedium = .true. ! to avoid threshold effects

    !BARYONS --------------------------------------------------
    if (useOffShellPotentialBaryons.and.isBaryon(partId)) then
       width=getBaryonWidth(partID,bareMass,momentum,mediumAtPosition,outofBounds)
    !MESONS --------------------------------------------------
    else if (useOffShellPotentialMesons.and.(partID==rho .or. partID==omegaMeson .or. partID==phi)) then
       width=WidthMesonMedium(partID,baremass,momentum,mediumAtposition)
    else
       success=.true.
       return
    end if

    poleMass=hadron(partID)%mass

    if (width.gt.0) then
       if (relativistic) then
         getOffShellParameter=(bareMass**2-poleMass**2)/(2.*abs4(momentum)*width)
       else
         getOffShellParameter=(bareMass-poleMass)/width
       end if
    else
       getOffShellParameter=0.
       success=.true.
       return
    end if


    !cutoff: we do not allow particles which are too far offshell
    if (isBaryon(partId) .and. abs(getoffshellparameter)>max_offshellparameter) then
       ! Declare failure: particle is too far off-shell
       if (counter_too_large.lt.100) then
          write(*,'(A,2G13.5)') 'WARNING: off-shell parameter too large!', getoffshellparameter, max_offshellparameter
          write(*,'(12g13.5)') partID,bareMass,momentum,poleMass,width,mediumAtPosition%density
          counter_too_large=counter_too_large+1
          if (counter_too_large==100) write(*,*) 'This happended now 100 times. Stopping output!'
       end if
       getoffshellparameter=0.
       success=.false.
    else
       success=.true.
    end if



    !for debugging:
    if (initHistFirst.and.offshell_debug) then
       call CreateHist2D(histo_nucl, 'off shell parameter',(/-20.,0./),(/20.,0.7/),(/0.05,0.01/))
       call CreateHist2D(histo_delta, 'off shell parameter',(/-20.,0./),(/20.,0.7/),(/0.05,0.01/))
       call CreateHist2D(histo_delta_rho, 'off shell parameter',(/-20.,0./),(/20.,0.2/),(/0.05,0.005/))
       initHistFirst=.false.
    end if

    if (offshell_debug.and.partID.eq.1) call AddHist2D(histo_nucl, (/getoffshellparameter,width/),1.)
    if (offshell_debug.and.partID.eq.2) call AddHist2D(histo_delta, (/getoffshellparameter,width/),1.)
    if (offshell_debug.and.partID.eq.2) call AddHist2D(histo_delta_rho, (/getoffshellparameter,mediumAtPosition%density/),1.)

    if (mod(numCalls,1000).eq.0.and.offshell_debug) then
       open(199,file='offshell_params_nucl.dat')
       call WriteHist2D_Gnuplot(histo_nucl,199,mul=histo_nucl%xBin(1)*histo_nucl%xBin(2))
       close(199)
       open(199,file='offshell_params_delta.dat')
       call WriteHist2D_Gnuplot(histo_delta,199,mul=histo_delta%xBin(1)*histo_delta%xBin(2))
       close(199)
       open(199,file='offshell_params_delta_rho.dat')
       call WriteHist2D_Gnuplot(histo_delta_rho,199,mul=histo_delta_rho%xBin(1)*histo_delta_rho%xBin(2))
       close(199)
    end if


  end function getOffShellParameter


  !****************************************************************************
  !****f* offShellPotential/getOffShellMass
  !NAME
  !real function getOffShellMass
  !
  !FUNCTION
  !This function calculates the offshellmass of a particle depending on its offshellparameter
  !and the current kinematics -> full four-momentum needed!!!
  !
  ! INPUTS
  ! * integer, intent(in) :: partID   -- ID of particle
  ! * real, intent(in) :: offshellparameter
  ! * real, intent(in) :: bareMass     -- bare mass= %mass
  ! * real, dimension(0:3), intent(in) :: momentum  -- 4-momentum
  ! * type(medium),intent(in) :: mediumAtPosition   -- medium information
  !
  ! OUTPUT
  ! * real ::                getOffShellMass
  ! * real, optional :: getOffShellness
  !****************************************************************************
  real function getOffShellMass(partID,offshellparameter,momentum,baremass,mediumAtPosition,outOfBounds)
    use particleDefinition
    use IdTable, only: isBaryon, isMeson
    use particleProperties, only: hadron
    use baryonWidthMedium_tables, only: get_minMass
    use mediumDefinition
    use dichteDefinition
    use minkowski, only: abs4, SP
    use callstack, only: traceback
    use mesonWidthMedium, only: WidthMesonMedium

    integer, intent(in) :: partID
    real, intent(in) :: offshellparameter
    real, dimension(0:3),intent(in) :: momentum
    real, intent(in) :: baremass
    type(medium),intent(in) :: mediumAtPosition

    logical, intent(out) :: outOfBounds
    real :: width

    if (initFlag) then
       call readInput
       initFlag=.false.
    end if

    outOfBounds=.false.
    if (treatParticleOffShell(partID,OffShellParameter)) then

       ! Baryons
       if (isBaryon(partId)) then

          if (SP(momentum,momentum).le.0.) then
             write(*,*) 'strange abs4 in offshell -> STOP in offshellPotential/getOffshellMass'
             call errorMessage_SPleZero()
             call traceback()
             stop
          end if
          width = getBaryonWidth(partID,bareMass,momentum,mediumAtPosition,outOfBounds)
          if (relativistic) then
            getOffShellMass=sqrt(max(hadron(partID)%mass**2+offshellparameter*2.*abs4(momentum)*width,0.))
          else
            getOffShellMass=hadron(partID)%mass+offshellparameter*width
          end if
          if (.not.extrapolateBaryonWidth) then
             getOffShellMass=max(getOffShellMass,get_minMass(partID))
             if (outofbounds) getOffShellMass=get_minMass(partID)
          end if

       ! Mesons
       else if (isMeson(partId)) then

          if (SP(momentum,momentum).le.0.) then
             write(*,*) 'Warning: strange abs4 in offshellPotential/getOffshellMass'
             call errorMessage_SPleZero()
             outOfBounds = .true.
             getOffShellMass = 0.
             return
          end if
          width = WidthMesonMedium(partID,baremass,momentum,mediumAtposition)
          if (relativistic) then
            getOffShellMass=sqrt(max(hadron(partID)%mass**2+offshellparameter*2.*abs4(momentum)*width,0.))
          else
            getOffShellMass=max(hadron(partID)%mass+offshellparameter*width,0.)
          end if

       else
          write(*,*) 'getOffShellMass: no meson/baryon! CRITICAL ERROR: SHOULD NOT BE CALLED: Strange ID!! -> STOP', partID
          stop
       end if

    else if (partID/=0) then
       write(*,*) 'getOffShellMass: treatParticleOffShell=.false.! CRITICAL ERROR: SHOULD NOT BE CALLED!! -> STOP',  &
            & partID,OffShellParameter
       stop
    else
       write(*,*) 'WARNING: getOffShellMass: ID=0'
       getOffShellMass=0
    end if

  contains

    subroutine errorMessage_SPleZero()
      write(*,*) 'p^mu p_mu is less or equal zero in getOffShellMass:', SP(momentum,momentum)
      write(*,*) 'partID,bareMass:'
      write(*,*) partID,bareMass
      write(*,*) 'momentum:'
      write(*,*) momentum
      write(*,*) 'offshellparameter:'
      write(*,*) offshellparameter
      write(*,*) 'density:'
      write(*,*) mediumAtPosition%density
    end subroutine errorMessage_SPleZero


  end function getOffShellMass


  !****************************************************************************
  !****f* offShellPotential/getBaryonWidth
  ! NAME
  ! real function getBaryonWidth
  !
  ! PURPOSE
  ! returns the baryonWidth. In case, extrapolateBaryonWidth is set to .true.
  ! the baryon width is extrapolated to masses below minimalmass.
  ! IMPORTANT: *** Baryon Width is calculated in the Lab Frame! ***
  !
  !
  ! INPUTS
  ! * integer, intent(in) :: partID   -- ID of particle
  ! * real, intent(in) :: bareMass     -- bare mass= %mass
  ! * real, dimension(0:3), intent(in) :: momentum  -- 4-momentum
  ! * type(medium),intent(in) :: mediumAtPosition   -- medium information
  !
  ! OUTPUT
  ! * real ::                getBaryonWidth
  !****************************************************************************
  real function getBaryonWidth(partID,bareMass,momentum,mediumAtPosition,outofBounds)
    use IdTable, only: isBaryon
    use particleProperties, only: hadron
    use baryonwidthmedium, only: WidthBaryonMedium
    use mediumDefinition
    use minkowski, only: abs4
    use twoBodyTools, only: pCM
    use constants, only: mN

    integer, intent(in) :: partID
    real, intent(in) :: bareMass
    real, dimension(0:3), intent(in) :: momentum
    type(medium),intent(in) :: mediumAtPosition
    logical,intent(out) :: outofBounds
    real :: width
    real :: x1,x2,y1,y2,m,b,s_cut,s_real
    integer :: Method=3
    real :: p_squared
    if (initFlag) then
       call readInput
       initFlag=.false.
    end if

    !Baryons
    if (isBaryon(partId)) then

       if (extrapolateBaryonWidth.and.baremass.lt.(hadron(partid)%minmass+0.015)) then
          x1=hadron(partid)%minmass+0.015
          y1=WidthBaryonMedium(partID,x1,momentum,mediumATposition,outofBounds)
          select case (method)
          case (1)
             ! tina's method
             x2=hadron(partid)%minmass+0.014
             y2=WidthBaryonMedium(partID,x2,momentum,mediumATposition,outofBounds)
             m=(y1-y2)/(x1-x2)
             b=y1-m*x1
             width=m*baremass+b
             if (offshell_debug) write(*,*) 'fakewidth=',Width
             width=max(0.001,width)  !minimalwidth as set in baryonWidthMedium_tables
             outofBounds=.false.
          case (2)
             ! Constant width
             width=y1
             getBaryonWidth=width
             outofBounds=.true.
          case (3)
             ! Assume that width is due to NN absorption below threshold
             p_squared=Dot_Product(momentum(1:3),momentum(1:3))
             s_cut =sqrt((mN+sqrt(x1**2      +p_squared) )**2-p_Squared )
             s_real=sqrt((mN+sqrt(baremass**2+p_squared) )**2-p_squared )
             width=y1*pcm(s_real,mN,mN)/pcm(s_cut,mN,mN)
             outofBounds=.false.
          end select
       else
          width=WidthBaryonMedium(partID,baremass,momentum,mediumATposition,outofBounds)
          if (offshell_debug) write(*,*) 'width=',Width,partID,baremass,momentum,mediumATposition,outofBounds
       end if

       ! Gamma_lab = Gamma_restframe / gamma
       getBaryonWidth = width * abs4(momentum)/momentum(0)

    else

       write(*,*) 'getBaryonWidth: no baryon! CRITICAL ERROR: SHOULD NOT BE CALLED: Strange ID!! -> STOP', partID
       stop
    end if

  end function getBaryonWidth



  !****************************************************************************
  !****s* offShellPotential/offShellErrorMessage(partID,bareMass,momentum,position)
  ! NAME
  ! subroutine  offShellErrorMessage(partID,bareMass,momentum,position)
  ! PURPOSE
  ! Routine returns particle information.
  ! INPUTS
  ! * type(particle), intent(in) :: part
  !****************************************************************************
  subroutine offShellErrorMessage(part)
    use particleDefinition
    use minkowski, only: SP

    type(particle), intent(in) :: part
    character(20) :: form
    form='(A,4G12.5)'

    write(*,'(A)')'offShellErrorMessage: printing particle properties...:'
    write(*,form)'Particle ID:                ', part%ID
    write(*,form)'Particle mass:              ', part%mass
    write(*,form)'Particle offshellparameter: ', part%offshellparameter
    write(*,form)'Particle number:            ', part%number
    write(*,form)'Particle firstevent:        ', part%firstevent
    write(*,form)'Particle history:           ', part%history
    write(*,form)'Particle momentum:          ', part%momentum
    write(*,form)'Particle perturbative?:     ', part%perturbative
    write(*,form)'Particle charge:            ', part%Charge
    write(*,form)'Particle position:          ', part%position
    write(*,form)'SP(momentum,momentum):      ', SP(part%momentum,part%momentum)
    write(*,form)'SQ(abs(SP(momentum,momentum):',sqrt(abs(SP(part%momentum,part%momentum)))
    write(*,*)

  end subroutine offShellErrorMessage



end module offShellPotential
