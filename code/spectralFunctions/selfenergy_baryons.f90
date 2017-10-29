!******************************************************************************
!****m* /selfenergy_baryons
! NAME
! module selfenergy_baryons
! PURPOSE
! * Implements the routines for the real part of the self energy for baryons
! * Based on dispersion relations
!******************************************************************************
module selfenergy_baryons
  use idTable, only: nres
  use tabulation, only: table_4dim
  use mediumDefinition
  implicit none
  private

  !  Public :: calc_RealPart,integrand
  public :: get_RealPart, selfenergy_imag, cleanup

  !****************************************************************************
  !****g* selfenergy_baryons/writeLocal
  ! SOURCE
  !
  logical,save  :: writeLocal=.false.
  !
  ! PURPOSE
  ! * Tables are outputted to local directory, not to buuinput
  !****************************************************************************

  !****************************************************************************
  !****g* selfenergy_baryons/debugFlag
  ! SOURCE
  !
  logical, parameter :: debugFlag=.false.
  !
  ! PURPOSE
  ! * Switch for debug information
  !****************************************************************************


  !****************************************************************************
  !****g* selfenergy_baryons/rel_accuracy
  ! SOURCE
  !
  real, save        :: rel_accuracy=0.05
  !
  ! PURPOSE
  ! Relative accuracy for resonance self energy
  !****************************************************************************

  !****************************************************************************
  !****g* selfenergy_baryons/rel_accuracy_nuc
  ! SOURCE
  !
  real, parameter   :: rel_accuracy_nuc=0.05
  !
  ! PURPOSE
  ! Relative accuracy for nucleon self energy
  !****************************************************************************

  !****************************************************************************
  !****g* selfenergy_baryons/intSolver
  ! SOURCE
  !
  integer, save :: intSolver=1
  !
  ! PURPOSE
  ! Decide on the numerical package to be used for the Cauchy integral:
  ! * 1=quadpack routine
  ! * 2=cernlib routine
  !****************************************************************************

  integer, parameter :: quad=1
  integer, parameter :: cernlib=2

  !****************************************************************************
  !****g* selfenergy_baryons/makeTable
  ! SOURCE
  !
  logical, save :: makeTable=.true.
  !
  ! PURPOSE
  ! Switch on/off the usage of an input tabulation
  !****************************************************************************


  !****************************************************************************
  !****g* selfenergy_baryons/noDispersion
  ! SOURCE
  !
  logical, save :: noDispersion=.false.
  !
  ! PURPOSE
  ! Switch on/off the usage dispersion relations
  !****************************************************************************

  !****************************************************************************
  !****g* selfenergy_baryons/extrapolateAbsP
  ! SOURCE
  !
  logical,save  :: extrapolateAbsP=.false.
  !
  ! PURPOSE
  ! if(true) then set absP to maxAbsP if absP is larger
  !****************************************************************************



  logical,save  :: readInput_flag=.true.


  ! Globals to set up the integrand:
  real, save      :: absP_g
  real, save      :: E_g
  real, save      :: E_pole_g
  integer, save   :: particleID_g
  integer,save    :: switch_g
  type(medium),save :: med_g


  logical, save :: tellErrors=.false.

  !****************************************************************************
  !****g* selfenergy_baryons/maxRes
  ! SOURCE
  integer,save :: maxRes=100
  ! PURPOSE
  !
  !****************************************************************************

  !****************************************************************************
  !****g* selfenergy_baryons/minRes
  ! SOURCE
  integer,save :: minRes=-100
  ! PURPOSE
  !
  !****************************************************************************


  integer,save :: funcID

  type(table_4dim),save,dimension(1:nres+1):: realPartTable

contains

  !****************************************************************************
  !****s* selfenergy_baryons/readInput
  ! NAME
  ! subroutine readInput
  ! PURPOSE
  ! Reads input in jobcard out of namelist "selfenergy_realPart".
  !****************************************************************************
  subroutine readInput
    use output, only: Write_ReadingInput

    integer :: ios ! checks file behavior

    !**************************************************************************
    !****n* selfenergy_baryons/selfenergy_realPart
    ! NAME
    ! NAMELIST /selfenergy_realPart/
    ! PURPOSE
    ! Includes the switches:
    ! * rel_accuracy
    ! * intSolver
    ! * makeTable
    ! * noDispersion
    ! * maxRes
    ! * minRes
    ! * extrapolateAbsP
    ! * writeLocal
    !**************************************************************************
    NAMELIST /selfenergy_realPart/ rel_accuracy, intSolver, makeTable, noDispersion, maxRes, minRes, &
                                   extrapolateAbsP, writeLocal

    call Write_ReadingInput('selfenergy_realPart',0)
    rewind(5)
    read(5,nml=selfenergy_realPart,IOSTAT=IOS)
    call Write_ReadingInput('selfenergy_realPart',0,IOS)

    if (noDispersion) then
       write(*,*) 'No Dispersion relation !!!!!!!! noDispersion=',noDispersion
    else
       write(*,*) 'Relative accuracy            ', rel_accuracy
       write(*,*) 'Which cauchy integral solver?', intSolver
       if (intSolver.eq.quad) then
          write(*,*) '  --> QUADPACK'
       else if (intSolver.eq.cernlib) then
          write(*,*) '  --> CERNLIB'
       else
          write(*,*) 'wrong input!'
          stop
       end if
       write(*,*) 'Make table                  ?', makeTable
       write(*,*) 'Maxres                      ?', maxRes
       write(*,*) 'Minres                      ?', minRes
    end if

    if (writeLocal)write(*,*)  ' WARNING : OUTPUT to local directory, and not to buuinput!'

    if (debugflag) write(*,*) 'debugging mode switched on in realPart'

    if (extrapolateAbsP) write(*,*) 'absP is set to maximum tabulated value'

    call Write_ReadingInput('selfenergy_realPart',1)

  end subroutine readInput

  !****************************************************************************
  !****************************************************************************
  subroutine cleanUp
    use tabulation, only: cleanupTable
    integer :: ID
    if (noDispersion) return
    do Id=max(minRes,1),min(maxRes,nres+1)
      call cleanupTable(realPartTable(Id))
    end do
  end subroutine cleanUp


  !****************************************************************************
  !****f* selfenergy_baryons/selfenergy_Imag
  ! NAME
  ! function selfenergy_Imag(particleID,absP,E,med,pos) result(imagSelf)
  !
  ! PURPOSE
  ! * Returns the imaginary part of the self energy of baryons.
  !
  ! INPUTS
  ! * integer      :: particleID    -- ID of baryon
  ! * real         :: absP          -- absolute Momentum
  ! * real         :: E             -- Energy of baryon
  ! * type(medium) :: med           -- medium information
  ! * real, dimension(3), OPTIONAL:: pos -- position of particle
  !
  ! OUTPUT
  ! * real :: imagSelf
  !****************************************************************************
  function selfenergy_imag(particleID,absP,E,med,pos) result(imagSelf)
    use idTable, only: nres
    use mediumDefinition

    integer, intent(in)              :: particleID
    real, intent(in)                 :: absP
    real, intent(in)                 :: E
    type(medium), intent(in)         :: med
    real, dimension(1:3),intent(in),OPTIONAL :: pos

    real :: imagSelf
    real :: mass

    imagSelf=0.

    if ((particleID.lt.1).or.(particleID.gt.nres+1)) then
       write(*,*) 'ERROR: particleID out of bounds in getFull_inMediumWidth',particleID, nres+1
    end if

    if (E.lt.0) return

    if (E**2-absP**2.lt.0) then
       !write(*,*) 'WARNING E**2-absP**2 less than zero', E, absP
       return
    else if ((E**2-absP**2).gt.7**2) then
       return
    else if ((E**2-absP**2).gt.5**2) then
       mass=sqrt(E**2-absP**2)
 !      imagSelf=(7.-mass)/2.*(-1.)*get_inMediumWidth(particleID,absP,5.,med,3)*mass
       imagSelf=(7.-mass)/2.*(-1.)*getWidth()*mass
       return
    else
       mass=sqrt(E**2-absP**2)
!       imagSelf=-get_inMediumWidth(particleID,absP,mass,med,3)*mass
       imagSelf=-getWidth()*mass
    end if

    contains
      real function getWidth()
        ! Assumptions: Resting nuclear matter, isotropic matter
        use baryonWidthMedium
        use potentialModule, only: massDetermination

        real,dimension(0:3) :: momLRF        ! momentum of resonance in LRF
        real :: baremass
        logical :: success

        momLRF=(/E,absp,0.,0./)

        baremass=-1234. !strange value set here on purpose!! is of no relevance, but might help
                        !debugging for nucleons in WidthBaryonMedium
        ! for nucleons: baremass is set to hadron(nucleon)%mass in WidthBaryonMedium_tables
        !               since we apply the assumption that the nucleon width is not mass dependent
        ! => no need to call massDetermination for nucleons
        if (particleID.ne.1) then
           if (present(pos)) then
              call massDetermination(particleID,momLRF,med,baremass,verbose=.false.,success=success,pos=pos)
           else
              call massDetermination(particleID,momLRF,med,baremass,verbose=.false.,success=success)
           end if
        else
           success=.true.
        end if
        if (.not.success) then
           getWidth=0.
        else
           getWidth=WidthBaryonMedium(particleID,baremass,momLRF,med)
        end if
      end function getWidth

  end function selfenergy_imag


  !****************************************************************************
  !****f* selfenergy_baryons/get_RealPart
  ! NAME
  ! function get_RealPart(particleID,absP,mass,med,pos) result(realPart)
  !
  ! PURPOSE
  ! * Returns the real part of the baryon self energy according to
  !   dispersion relations
  !
  ! INPUTS
  ! * integer  :: particleID    -- ID of baryon
  ! * real     :: absP          -- absolute Momentum
  ! * real     :: mass          -- Mass of baryon
  ! * type(medium) :: med       -- medium information
  ! * real, dimension(3), OPTIONAL:: pos -- position of particle
  !
  ! OUTPUT
  ! * real :: realPart
  !****************************************************************************
  function get_RealPart(particleID,absP,mass,med,pos) result(realPart)
    use inputGeneral, only: path_to_input
    use output, only: intToChar
    use idTable, only: nres, delta_idTable =>delta
    use mediumDefinition
    use particleProperties, only: hadron
    use idTable, only: nucleon
    use baryonWidthMedium, only: get_MediumSwitch,get_MediumSwitch_coll, &
         & get_MediumSwitch_Delta,get_MediumSwitch_proton_neutron
    use baryonWidthMedium_tables, only: get_deltaOset
    use mediumDefinition
    use tabulation, only: init_table_4dim, readTable, getValue, printTable

    integer, intent(in)              :: particleID
    real, intent(in)                 :: absP
    real, intent(in)                 :: mass
    type(medium), intent(in)         :: med
    real, dimension(3), intent(in), optional :: pos
    real :: realPart

    real, parameter               :: max_absP=2.9
    real, parameter               :: max_rhoP=0.1 ! In fm**-3
    real, parameter               :: max_rhoN=0.1 ! In fm**-3

    real, save ::  delta_rhoN=0.025
    real, save ::  delta_rhoP=0.025
    real, save ::  delta_absP=0.04
    real, dimension(1:nres+1) :: delta_mass
    real, dimension(1:nres+1) :: min_mass
    real, dimension(1:nres+1) :: max_mass
    integer, save  :: numSteps_mass=400

    logical,save  :: firstTime=.true.
    integer :: outOfBounds,ID,intSolver_in,ios
    real(4), dimension (1:4) :: x
    real(4), dimension (1:4),save :: delta, minimum,maximum
    character(100) :: filename
    real :: energy,rel_accuracy_nuc_in, rel_accuracy_in,rhoN_new,rhoP_new,dummy
    logical :: noDispersion_in,forceNewCalculation
    logical :: in_MediumSwitch,in_MediumSwitch_coll,in_MediumSwitch_Delta,in_MediumSwitch_proton_neutron
    integer, save :: outofbounds_counter=0

    if (readInput_flag) then
       call readinput
       readInput_flag=.false.
    end if

    if ((particleID.lt.1).or.(particleID.gt.nres+1)) then
       write(*,*) 'ERROR: particleID out of bounds in getFull_inMediumWidth',particleID, nres+1
       realPart=0.
       return
    end if

    if (makeTable.and.(.not.nodispersion).and.(.not.(particleID.eq.delta_idTable.and.get_deltaOset()))) then
       ! Dispersion must we switched off for the Oset delta width since the width is not meaningful for high energy!!!!
       if (firstTime) then
          ! Generate lookup table
          ! (a) Check whether there is a table preconfigured which was using the same switches
          filename=trim(path_to_input)//'/inMediumWidth/RealPart_switches.dat'
          open(77,file=filename,iostat=ios,status='old')
          if (ios.ne.0) then
             forceNewCalculation=.true.
             close(77)
             write(*,*) 'Failed to read real part switches from file (1) -> Redo the parametrizations'
             !call writeParameters()
          else
             read(77,*,iostat=ios) rel_accuracy_in,rel_accuracy_nuc_in,intSolver_in,noDispersion_in
             if (ios.ne.0) then
                forceNewCalculation=.true.
                close(77)
                write(*,*) 'Failed to read real part switches from file (2) -> Redo the parametrizations'
                !call writeParameters()
             else
                if (abs(rel_accuracy-rel_accuracy_in).gt.1E-4.or.abs(rel_accuracy_nuc-rel_accuracy_nuc_in).gt.1E-4 &
                     & .or.intSolver.ne.intSolver_in.or.(nodispersion.neqv.noDispersion_in)) then
                   forceNewCalculation=.true.
                   close(77)
                   write(*,*) 'Real part switches changed -> Redo the parametrizations'
                   !call writeParameters()
                else
                   forceNewCalculation=.false.
                   write(*,*) 'Same real part switches as the parametrization'
                end if
             end if
             if (.not.forceNewCalculation) then
                read(77,*,iostat=ios) in_MediumSwitch,in_MediumSwitch_coll,&
                     & in_MediumSwitch_Delta,in_MediumSwitch_proton_neutron
                if (ios.ne.0) then
                   forceNewCalculation=.true.
                   close(77)
                   write(*,*) 'Failed to read real part switches from file (2) -> Redo the parametrizations'
                   !call writeParameters()
                else
                   if (     (in_MediumSwitch.neqv.get_MediumSwitch() ).or. &
                        &  (in_MediumSwitch.and.( &
                        &       (in_MediumSwitch_coll.neqv.get_MediumSwitch_coll()   ).or. &
                        &       ((in_MediumSwitch_delta.neqv.get_MediumSwitch_delta()).and.(particleID.eq.delta_idTable )).or. &
                        &       (in_MediumSwitch_proton_neutron.neqv.get_MediumSwitch_proton_neutron() )&
                        &  )) &
                        &) then

                      forceNewCalculation=.true.
                      close(77)
                      write(*,*) 'Width switches changed -> Redo the parametrizations'
                      write(*,*) 'in_MediumSwitch', in_MediumSwitch,get_MediumSwitch()
                      write(*,*) 'in_MediumSwitch_coll', in_MediumSwitch_coll,get_MediumSwitch_coll()
                      write(*,*) 'in_MediumSwitch_Delta',in_MediumSwitch_Delta,get_MediumSwitch_delta()
                      write(*,*) 'in_MediumSwitch_proton_neutron',in_MediumSwitch_proton_neutron,get_MediumSwitch_proton_neutron()
                      !call writeParameters()
                   else
                      forceNewCalculation=.false.
                      write(*,*) 'Same width switches as the parametrization'
                   end if
                end if
             end if
          end if
          ! (b) read out or generate table
          if (forceNewCalculation) call writeParameters(.false.)
          do Id=max(minRes,1),min(maxRes,nres+1)
             write(*,*) "ID = ",ID
             min_mass(ID)=hadron(ID)%minmass+0.01
             if (ID.eq.nucleon) min_mass(ID)=0.2
             max_mass(ID)=2.5
             delta_mass(ID)=(max_mass(ID)-min_mass(ID))/float(numSteps_mass)

             minimum=(/0.,min_Mass(ID),0.,0./)
             maximum=(/max_absP,max_Mass(ID),max_rhoN,max_rhoP/)
             delta=(/delta_absP,delta_Mass(ID),delta_rhoN,delta_rhoP/)
             if (writeLocal) then
                filename='RealPart_'//intToChar(ID)//'.dat.bz2'
             else
                filename=trim(path_to_input)//'/inMediumWidth/RealPart_'//intToChar(ID)//'.dat.bz2'
             end if
             if (forceNewCalculation) then
                funcID=ID
                realPartTable(ID)=init_table_4dim(minimum,maximum,delta,func,.true.)
                call printTable(realPartTable(ID),trim(filename))
             else if (.not.readTable(realPartTable(ID),trim(filename),minimum,delta)) then
                funcID=ID
                realPartTable(ID)=init_table_4dim(minimum,maximum,delta,func,.true.)
                call printTable(realPartTable(ID),trim(filename))
             end if
          end do
          if (forceNewCalculation) call writeParameters(.true.)

          firstTime=.false.
          write(*,*) 'real part table finished'
       end if

       rhon_new=max(med%densityNeutron,1.E-8)
       rhop_new=max(med%densityProton,1.E-8)

       x=(/absP,mass,rhoN_new,rhoP_new/)

       if (extrapolateAbsP) then
          if (absP.gt.max_absP)  x=(/max_absP,mass,rhoN_new,rhoP_new/)
       end if

       realPart=getValue(x,realPartTable(particleId),outOfBounds)

       if (outOfBounds.ne.0) then
          if (outofbounds_counter.lt.100) then
             write(*,'(A,4F9.4)') 'WARNING: Out of bounds in get_inMediumRealPart',&
                  & absP,mass,med%densityNeutron,med%densityProton
             write(*,'(A,4F9.4)') 'Min=                                          ',realPartTable(particleId)%min
             write(*,'(A,4F9.4)') 'Max=                                          ',realPartTable(particleId)%max
             outofbounds_counter=outofbounds_counter+1
             if (outofbounds_counter.eq.100) write(*,*) 'THIS HAPPENED NOW 100 times: STOPPPING OUTPUT NOW !!!'
          end if

          if (mass.lt.min_mass(particleID)) then
             realPart=0.  !the width is anyway zero below min_mass(particleID)=hadron(ID)%minmass+0.01 - avoids time-consuming calculation of realPart
          else
             energy=sqrt(mass**2+absP**2)
             if (present(pos)) then
                realPart=calc_RealPart(particleID,absP,energy,med,pos=pos)
             else
                realPart=calc_RealPart(particleID,absP,energy,med)
             end if
          end if
       else
          ! The real part table include the EQS=5 contribution at the pole. The following
          ! line removes that contribution and adds the contribution according to the present EQS type.

          if (present(pos)) then
             realPart=realPart+get_realSelf_AtPole(particleID,absP,med,dummy,pos=pos)   &
                  &           -get_realSelf_AtPole(particleID,absP,med,dummy,5,pos=pos)
          else
             realPart=realPart+get_realSelf_AtPole(particleID,absP,med,dummy)   &
                  &           -get_realSelf_AtPole(particleID,absP,med,dummy,5)
          end if
       end if
    else

       energy=sqrt(mass**2+absP**2)
       if (present(pos)) then
          realPart=calc_RealPart(particleID,absP,energy,med,pos=pos)
       else
          realPart=calc_RealPart(particleID,absP,energy,med)
       end if
    end if

  contains

    subroutine writeParameters(flag)
      logical, intent(in) :: flag
      if (writeLocal) then
         filename='RealPart_switches.dat'
      else
         filename=trim(path_to_input)//'/inMediumWidth/RealPart_switches.dat'
      end if
      open(107,file=filename,iostat=ios)
      if (flag) then
         write(107,'(2E20.5,I4,L8)') rel_accuracy,rel_accuracy_nuc,intSolver,noDispersion
         write(107,'(4L8)') get_MediumSwitch(),get_MediumSwitch_coll(), &
              & get_MediumSwitch_Delta(),get_MediumSwitch_proton_neutron()
      else
         write(107,*) 'Inconsistent real parts!!'
      end if
      close(107)
    end subroutine writeParameters

  end function get_RealPart


  real function func(a,b,c,d)
    use mediumDefinition

    real,intent(in) :: a,b,c,d
    type(medium)       :: med

    med%density = c+d
    med%densityProton = d
    med%densityNeutron = c
    med%useMedium = .true.

    func=calc_RealPart(funcID,a,sqrt(a**2+b**2),med)
  end function func

  !****************************************************************************
  !****f* selfenergy_baryons/calc_RealPart
  ! NAME
  ! function RealPart(particleID,absP,E,med,pos) result(realSelf)
  !
  !
  ! PURPOSE
  ! * Returns the real part of the self energy of baryons according to dispersion relations.
  ! * Units: [GeV**2]
  !
  ! INPUTS
  ! * integer      :: particleID    -- ID of baryon
  ! * real         :: absP          -- absolute Momentum
  ! * real         :: E             -- Energy of baryon
  ! * type(medium) :: med           -- medium information
  ! * real, dimension(3), OPTIONAL:: pos -- position of particle
  !
  ! OUTPUT
  ! * real :: realSelf
  !
  !****************************************************************************
  function calc_RealPart(particleID,absP,E,med,pos) result(realSelf)
    use idTable, only: nres,delta
    use particleProperties, only: hadron
    use particleDefinition
    use potentialModule
    use baryonWidthMedium_tables, only: get_DeltaOset
    use constants, only: pi
    use mediumDefinition

    integer, intent(in)              :: particleID
    real, intent(in)                 :: absP
    real, intent(in)                 :: E
    type(medium), intent(in)       :: med
    real, dimension(3), intent(in), optional :: pos
    real :: realSelf
    real :: realSelf_atPole,pot_atPole,E_pole

    if ((particleID.lt.1).or.(particleID.gt.nres+1)) then
       write(*,*) 'ERROR: particleID out of bounds in inMediumRealPart',particleID, nres+1
       realSelf=0.
       return
    end if

    if (readInput_flag) then
       call readinput
       readInput_flag=.false.
    end if

    ! Get the potential at the pole position
    if (present(pos)) then
       realSelf_atPole=get_realSelf_AtPole(particleID,absP,med,pot_AtPole,pos=pos)
    else
       realSelf_atPole=get_realSelf_AtPole(particleID,absP,med,pot_AtPole)
    end if
    E_pole=sqrt(absP**2+hadron(particleID)%mass**2)+pot_atPole

    if (noDispersion.or.(particleID.eq.delta.and.get_deltaOset())) then
       ! Dispersion must we switched off for the Oset width since the width is not meaningful for high energy!!!!
       realSelf=realSelf_atPole
    else
       realSelf=realSelf_atPole+(E-E_pole)/pi*principalValue(particleID,absP,E,E_pole,med)
    end if
  end function calc_RealPart



  !****************************************************************************
  !****f* selfenergy_baryons/get_realSelf_AtPole
  ! NAME
  ! real function get_realSelf_AtPole(particleID,absP,med,pot_AtPole,EQS,pos)
  !
  ! PURPOSE
  ! * Returns the real part of the self energy of baryons at the pole energy
  ! * Units: [GeV**2]
  !
  ! INPUTS
  ! * integer           :: particleID    -- ID of baryon
  ! * real              :: absP          -- absolute momentum
  ! * type(medium)      :: med           -- medium information
  ! * integer,optional  :: EQS           -- if present then we use this EQS
  !   type to evaluate the potential, else the default is used
  ! * real, dimension(3), OPTIONAL:: pos -- position of particle
  !
  ! OUTPUT
  ! * real :: realSelf at the pole position with mass=vacuum mass
  ! * real :: pot_AtPole  -- potential at pole position
  !
  !****************************************************************************
  real function get_realSelf_AtPole(particleID,absP,med,pot_AtPole,EQS,pos)
    use idTable, only: nres
    use particleProperties, only: hadron
    use particleDefinition
    use baryonPotentialModule, only: baryonPotential
    use mediumDefinition

    integer, intent(in)            :: particleID
    real, intent(in)               :: absP
    type(medium), intent(in)       :: med
    integer,intent(in),optional    :: EQS
    real, dimension(3), intent(in), optional :: pos
    real, intent(out)              :: pot_AtPole

    real :: potScalar_atPole
    type(particle) :: part

    if ((particleID.lt.1).or.(particleID.gt.nres+1)) then
       write(*,*) 'ERROR: particleID out of bounds in inMediumRealPart',particleID, nres+1
       get_realSelf_AtPole=0.
       return
    end if

    if (readInput_flag) then
       call readinput
       readInput_flag=.false.
    end if

    ! Define a particle at the pole:
    part%ID=particleID
    part%momentum(1:3)= (/absP,0.,0./)
    part%momentum(0)=-9999999999. ! Dummy, just to make somebody wonder if he wants to use momentum(0) in baryonPotential
    part%mass=hadron(particleID)%mass
    part%perturbative=.false.
    if (present(pos)) part%position = pos

    ! Get the potential at the pole position
    if (present(EQS)) then
       pot_atPole= BaryonPotential(part,med,.not.(present(pos)),EQS)
    else
       pot_atPole= BaryonPotential(part,med,.not.(present(pos)))
    end if
    potScalar_atPole=sqrt((sqrt(absP**2+hadron(particleID)%mass**2)+pot_atPole)**2-absP**2)-hadron(particleID)%mass

    get_realSelf_atPole=2.*hadron(particleID)%mass*potScalar_atPole+potScalar_atPole**2

  end function get_realSelf_AtPole


  !****************************************************************************
  !****f* selfenergy_baryons/principalValue
  ! NAME
  ! function principalValue(particleID,absP,E,E_pole,med) result(princi)
  !
  ! PURPOSE
  ! * Returns the principal value of the dispersion relation integral
  !
  ! INPUTS
  ! * integer      :: particleID    -- ID of baryon
  ! * real         :: absP          -- absolute Momentum
  ! * real         :: E             -- energy
  ! * real         :: E_pole        -- Energy at pole = subtraction point
  ! * type(medium) :: med -- medium information
  !
  ! OUTPUT
  ! * real :: princi
  !****************************************************************************
  real function principalValue(particleID,absP,E,E_pole,med) result(princi)
    use quadPack
    use particleProperties, only: hadron
    use cern_lib, only: dcauch
    use mediumDefinition

    integer, intent(in)              :: particleID
    real  , intent(in)               :: absP
    real ,intent(in)                 :: E_pole, E
!    real, intent (in)                :: rhoN,rhoP ! In fm^-3
    type(medium), intent(in)       :: med

    !real :: dummy
    real :: result_low,result_up, abserr
    integer :: neval,ier
    real, parameter :: cutOff=20.
    real  :: lowerCutOff
    real :: a,b

    real :: abs_acc, rel_acc

    abs_acc=rel_accuracy_nuc*(E**2-absP**2-hadron(particleID)%mass**2)
    rel_acc=rel_accuracy

    lowerCutOff=hadron(particleID)%minmass+0.002
    ! initialise integrand

    if (E_pole.gt.E) then
       ! Treat the pole at "E"
       if (debugFlag) then
          write(*,*) 'E<E_pole', E,E_pole
       end if
       select case (IntSolver)
       case (quad)
          a=lowerCutoff
          b=E+(E_pole-E)/2.
          if (a.gt.b) then
             write(*,*) 'cutoff too small: a>b',a,b
             a=b-0.01
             !stop
          end if
          call init_integrand(absP,med,E,E_pole,particleID,2)
          call qawc ( integrand, a, b, E, abs_acc, rel_acc, result_low, abserr, neval, ier )
          if (debugFlag.or.ier.ne.0)  call errorMessage_qawc(neval,ier,absErr,result_low,'low')
          ! Treat the pole at "E_pole"
          a=E+(E_pole-E)/2.
          b=cutoff
          if (a.gt.b) then
             write(*,*) 'cutoff too small: a>b',a,b
             stop
          end if
          call init_integrand(absP,med,E,E_pole,particleID,3)
          call qawc ( integrand, a, b, E_pole,abs_acc, rel_acc, result_up, abserr, neval, ier )
          if (debugFlag.or.ier.ne.0)  call errorMessage_qawc(neval,ier,absErr,result_up,'up')
       case (cernlib)
          !          write(*,*) 'cernlib'
          call init_integrand(absP,med,E,E_pole,particleID,1)
          a=lowerCutOff
          b=E+(E_pole-E)/2.
          if (a.gt.b) then
             write(*,*) 'cutoff too small: a>b',a,b
             stop
          end if
          result_low=dCAUCH(integrandD,dble(A),dble(B),dble(E),dble(abs_acc))
          a=E+(E_pole-E)/2.
          b=cutoff
          if (a.gt.b) then
             write(*,*) 'cutoff too small: a>b',a,b
             stop
          end if
          result_up=dCAUCH(integrandD,dble(A),dble(B),dble(E_pole),dble(abs_acc))
       end select

    else if (abs(E_pole-E).lt.0.00001) then
       write(*,*) 'WARNING: E_pole=E'
       princi=0.
       return
    else
       ! Treat the pole at "E_pole"
       select case (IntSolver)
       case (quad)
          a=lowerCutOff
          b=E_pole+(E-E_pole)/2.
          if (a.gt.b) then
             write(*,*) 'cutoff too small: a>b',a,b
             stop
          end if
          if (debugFlag) then
             write(*,*) 'E>E_pole', E,E_pole
             write(*,*) 'a,b', a,b
          end if
          call init_integrand(absP,med,E,E_pole,particleID,3)
          call qawc ( integrand, a, b, E_pole,abs_acc, rel_acc, result_low, abserr, neval, ier )
          if (debugFlag.or.ier.ne.0)  call errorMessage_qawc(neval,ier,absErr,result_low,'low 2')
          ! Treat the pole at "E"
          a=E_pole+(E-E_pole)/2.
          b=cutoff
          if (a.gt.b) then
             write(*,*) 'cutoff too small: a>b',a,b
             stop
          end if
          if (debugFlag)  write(*,*) 'a,b', a,b
          call init_integrand(absP,med,E,E_pole,particleID,2)
          call qawc ( integrand, a, b, E, abs_acc, rel_acc, result_up, abserr, neval, ier )
          if (debugFlag.or.ier.ne.0)  call errorMessage_qawc(neval,ier,absErr,result_up,'up 2')

       case (cernlib)
          call init_integrand(absP,med,E,E_pole,particleID,1)
          a=lowerCutOff
          b=E_pole+(E-E_pole)/2.
          if (a.gt.b) then
             write(*,*) 'cutoff too small: a>b',a,b,'cernlib low'
             stop
          end if
          !          write(*,*) '1',a,b
          result_low=dCAUCH(integrandD,dble(A),dble(B),dble(E_pole),dble(abs_acc))
          a=E_pole+(E-E_pole)/2.
          b=cutoff
          if (a.gt.b) then
             write(*,*) 'cutoff too small: a>b',a,b,'cernlib up'
             stop
          end if
          !          write(*,*) '2',a,b
          result_up=dCAUCH(integrandD,dble(A),dble(B),dble(E),dble(abs_acc))
       end select


    end if
    princi=result_low+result_up

    contains
      subroutine errorMessage_qawc(neval,ier,absErr,resu,who)
        real :: absErr, resu
        integer :: ier, neval
        character(*) :: who
        real :: relError

        if (abs(resu).gt.0) then
           relError=abs(absErr/resu)
        else
           relError=1.
        end if

        if (debugFlag.or.relError.gt.0.1) then
           if (relError.gt.0.3) write(*,*) 'CRITICAL ERROR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
           write(*,'(3A,2I6,E15.4,2(A,E15.4))')  'Problem in qawc ',who,':', neval, ier,absErr, &
                & 'result=',resu, ' rel. Error=', relError
        end if

      end subroutine errorMessage_qawc


  end function principalValue

  !****************************************************************************
  !****s* selfenergy_baryons/init_integrand
  ! NAME
  ! subroutine init_integrand(absP_in,med_in,E_in,E_pole_in,particleID_in,switch_in)
  !
  ! PURPOSE
  ! * Sets global variables which are used by the integrand of the dispersion integral
  !****************************************************************************
  subroutine init_integrand(absP_in,med_in,E_in,E_pole_in,particleID_in,switch_in)
    use mediumDefinition

    real, intent(in)    :: absP_in,E_pole_in,E_in
    integer, intent(in) :: particleID_in, switch_in
    type(medium), intent(in) :: med_in

    ! Reset the constants:
    absP_g=absP_in
    med_g=med_in
    E_pole_g=E_pole_in
    E_g=E_in
    particleID_g=particleID_IN
    switch_g=switch_in
    return
  end subroutine init_integrand

  !****************************************************************************
  !****f* selfenergy_baryons/integrand
  ! NAME
  ! real function integrand(Eprime)
  !
  ! PURPOSE
  ! * Integrand of the dispersion relation
  !****************************************************************************
  real function integrand(Eprime)

    real :: Eprime

    select case (switch_g)
    case (1)
       ! Full integrand
       if (abs(Eprime-E_g).lt.0.0000000001) then
          if (tellErrors) write(*,*) 'Error eprime=E', E_g, Eprime,E_pole_g
          integrand=0.
          return
       else if (abs(Eprime-E_pole_g).lt.0.0000000001) then
          if (tellErrors) write(*,*) 'Error eprime=E_pole', E_g, Eprime,E_pole_g
          integrand=0.
          return
       end if
       integrand=1./(Eprime-E_g)/(Eprime-E_pole_g)  *selfEnergy_Imag(particleID_g,absP_g,EPrime,med_g)
    case (2)
       integrand=1.             /(Eprime-E_pole_g)  *selfEnergy_Imag(particleID_g,absP_g,EPrime,med_g)

    case (3)
       integrand=1./(Eprime-E_g)                    *selfEnergy_Imag(particleID_g,absP_g,EPrime,med_g)
    end select

  end function integrand

  double precision function integrandD(Eprime)
    real :: Eprime
    integrandD = dble(integrand(Eprime))
  end function integrandD

end module selfenergy_baryons
