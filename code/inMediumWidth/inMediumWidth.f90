!******************************************************************************
!****m* /inMediumWidth
! NAME
! module inMediumWidth
! PURPOSE
! * Implements the routines for the medium width of the baryons
! * Generates tabulated input files
! * see also code/inMediumWidth/testRun/tabulateImagPart.f90 which
!   creates an executable using this routine to tabulate the width
!******************************************************************************
module inMediumWidth

  private

  !****************************************************************************
  !****g* inMediumWidth/debugFlag
  ! SOURCE
  !
  logical,save  :: debugFlag=.false.
  !
  ! PURPOSE
  ! * Switch for debug information
  !****************************************************************************

  !****************************************************************************
  !****g* inMediumWidth/writeLocal
  ! SOURCE
  !
  logical,save  :: writeLocal=.false.
  !
  ! PURPOSE
  ! * If .true. then output files are written to ./"filename", else to buuinput.
  !****************************************************************************

  !****************************************************************************
  !****g* inMediumWidth/maxRes
  ! SOURCE
  !
  integer,save  :: maxRes=1000
  !
  ! PURPOSE
  ! Read the data table up to a maximum resonance ID. ONLY FOR TESTING!!!
  !****************************************************************************

  !****************************************************************************
  !****g* inMediumWidth/minRes
  ! SOURCE
  !
  integer,save  :: minRes=-1000
  !
  ! PURPOSE
  ! Read the data table starting at this minimal resonance ID. ONLY FOR TESTING!!!
  !****************************************************************************


  ! flags for mesons
  integer, save :: maxMes = 1000      ! max. meson ID to tabulate
  integer, save :: minMes = -1000     ! min. meson ID to tabulate
  real,    save :: max_absP_Mes = 3.0 ! max. momentum for tabulation
  integer, save :: num_MonteCarlo_Points_mesons = 250 ! number of monte-carlo points
  integer, save :: cmin = 0, cmax = 0 ! bounds of charge loop

  logical,save  :: readInput_flag=.true.



  public :: tabulate_inMediumWidth_baryons,tabulate_inMediumWidth_mesons
contains


  !****************************************************************************
  !****s* inMediumWidth/readInput
  ! NAME
  ! subroutine readInput
  ! PURPOSE
  ! Reads input in jobcard out of namelist "inMediumWidth".
  !****************************************************************************
  subroutine readInput

    use output

    implicit none
    integer :: ios ! checks file behavior

    NAMELIST /inMediumWidth/ debugFlag,maxRes,minRes,writeLocal, &
                               maxMes,minMes,max_absP_Mes,num_MonteCarlo_Points_mesons,cmin,cmax

    call Write_ReadingInput('inMediumWidth',0)
    rewind(5)
    read(5,nml=inMediumWidth,IOSTAT=IOS)
    call Write_ReadingInput('inMediumWidth',0,IOS)

    write(*,*) 'Debugging information        ?   ', debugFlag
    write(*,*) 'maxRes                       ?   ', maxRes
    write(*,*) 'minRes                       ?   ', minRes
    write(*,*) 'write locally                ?   ', writeLocal

    write(*,*) 'maxMes                       ?   ', maxMes
    write(*,*) 'minMes                       ?   ', minMes
    write(*,*) 'max_absP_Mes                 ?   ', max_absP_Mes
    write(*,*) 'num_MonteCarlo_Points_mesons ?   ', num_MonteCarlo_Points_mesons
    write(*,*) 'cmin                         ?   ', cmin
    write(*,*) 'cmax                         ?   ', cmax

    call Write_ReadingInput('inMediumWidth',1)
  end subroutine readInput


  !****************************************************************************
  !****s* inMediumWidth/tabulate_inMediumWidth_baryons
  ! NAME
  ! subroutine tabulate_inMediumWidth_baryons()
  !
  ! PURPOSE
  ! * Tabulates the in-medium width of baryons according to
  !   Gamma=Gamma_free*Pauli_Blocking+sigma*rho*v *Pauli_Blocking
  ! * An average over the Fermi sea is performed
  ! * Sigma is given by the actual full collision term, with the incoming particles in the vacuum
  ! * Simplification: An average over the charge of the baryon is performed.
  !
  ! INPUTS
  ! * NONE
  !
  ! OUTPUT
  ! * Files "inMediumWidth/InMediumWidth.particleID.1.dat.bz2" and
  !   "inMediumWidth/InMediumWidth.particleID.2.dat.bz2"
  !****************************************************************************
  subroutine tabulate_inMediumWidth_baryons()
    use idTable, only: nres,nucleon
    use mediumDefinition
    use particleProperties, only: hadron
    use inputGeneral, only: path_To_Input
    use output, only: intToChar,line
    use baryonWidthMedium_tables, only: numSteps_rhoN,numSteps_rhoP,numSteps_absP,numsteps_mass,max_absP,max_rhoP,max_rhoN, &
         & num_MonteCarlo_Points, get_min_charge_loop, get_max_charge_loop
    use bzip

    implicit none

    real, save ::  delta_rhoN,delta_rhoP,delta_absP
    real, save, dimension(1:nres+1) :: delta_mass,min_mass,max_mass
    real, save, dimension(:,:,:,:,:) , Allocatable:: widthTable
    ! (:,:,:,:,:,1) = collisional broadening
    ! (:,:,:,:,:,2) = pauli blocked vacuum width

    integer :: particleID,index_absP,index_mass, index_rhoN,index_rhoP,dummy,i
    integer, parameter :: showColl=1,showPauli=2
    character(20) :: format
    character(200) ::fileName
    real :: rhoN, rhoP,absP, mass
    type(bzFile) :: f
    character(len=100) :: buf

    if (readInput_flag) then
       call readinput
       readInput_flag=.false.

       write(*,*) 'tabulation parameters:'
       write(*,*) 'numSteps_rhoN= ',numSteps_rhoN
       write(*,*) 'numSteps_rhoP= ',numSteps_rhoP
       write(*,*) 'numSteps_absP= ',numSteps_absP
       write(*,*) 'numSteps_mass= ',numSteps_mass
       write(*,*) 'max_absP= ',max_absP
       write(*,*) 'max_rhoP= ',max_rhoP
       write(*,*) 'max_rhoN= ',max_rhoN
       write(*,*) 'num_MonteCarlo_Points= ',num_MonteCarlo_Points
       ! The following construction is necessary to avoid a recursive I/O operation:
       dummy=get_min_charge_loop()
       write(*,*) 'min_charge_loop= ',dummy
       dummy=get_max_charge_loop()
       write(*,*) 'max_charge_loop= ',dummy
    end if

    allocate(widthTable(0:numSteps_absP,0:numSteps_mass,0:numSteps_rhoN,0:numSteps_rhoP,1:2))
    delta_rhoN=max_rhoN/float(numSteps_rhoN)
    delta_rhoP=max_rhoP/float(numSteps_rhoP)
    delta_absP=max_absP/float(numSteps_absP)

    format='(' // intToChar(numSteps_rhoP+1) // 'E12.4)'

    ! Set the bounds for the parametrization of the mass
    do particleId=1,nres+1
       !         min_mass(particleID)=max(minimalMass(particleID)+0.02,hadron(particleID)%mass-max(0.1,hadron(particleID)%width*3.))
       !         max_mass(particleID)=hadron(particleID)%mass+max(0.1,hadron(particleID)%width*3.)
       if (particleID.eq.nucleon) then
          min_mass(particleID)=0.80
       else
          min_mass(particleID)=hadron(particleID)%minmass+0.01
       end if
       max_mass(particleID)=7.
       delta_mass(particleID)=(max_mass(particleID)-min_mass(particleID))/float(numSteps_mass)
    end do

    write(*,*) 'Tabulating IN-MEDIUM WIDTH table for BARYONS'
    do particleId=max(minRes,1),min(maxRes,nres+1)
       do index_absP=0,numSteps_absP
          absP=float(index_absP)*delta_absP
          write(*,*) (particleID-1)*(numSteps_absP+1)+index_absP,'/',(nres+1)*(numSteps_absP+1)
          write(*,'(A,F15.4,2(A,I8))') 'absP=', absP,' steps:',index_absP,'/',numsteps_absP
          do index_mass=0,numSteps_mass
             mass=min_mass(particleID)+float(index_mass)*delta_mass(particleID)
             write(*,'(A,F15.4,2(A,I8))') 'mass=', mass,' steps:',index_Mass,'/',numsteps_mass
             do index_rhoN=0,numSteps_rhoN
                rhoN=float(index_rhoN)*delta_rhoN
                do index_rhoP=0,numSteps_rhoP
                   if (mass.lt.hadron(particleID)%minmass) then
                      widthTable(index_absP,index_mass,index_rhoN,index_rhoP,showPauli)= 0.
                      widthTable(index_absP,index_mass,index_rhoN,index_rhoP,showColl)=  0.
                   else
                      rhoP=float(index_rhoP)*delta_rhoP
                      widthTable(index_absP,index_mass,index_rhoN,index_rhoP,showPauli)=&
                           & get_pauliBlockedDecayWidth(particleID,absP,mass,rhoN,rhoP)
                      widthTable(index_absP,index_mass,index_rhoN,index_rhoP,showColl)=&
                           & evaluateCollisionBroadening_baryons(particleID,absP,mass,rhoN,rhoP)
                   end if
                end do
             end do
          end do
       end do

       ! Write the table to files:
       do i=1,2
          if (writelocal) then
             fileName='./InMediumWidth.'//trim(intToChar(particleID))&
                  & //'_'//trim(intToChar(i))//'.dat.bz2'
          else
             fileName=trim(path_to_Input)//'/inMediumWidth/InMediumWidth.'//trim(intToChar(particleID))&
                  & //'_'//trim(intToChar(i))//'.dat.bz2'
          end if
          f = bzOpenW(trim(fileName))
          write(buf,'(A,7I6)')'#', numSteps_absP, numSteps_mass , numSteps_rhoN, numSteps_rhoP,num_MonteCarlo_Points,&
               & get_min_charge_loop(),get_max_charge_loop()
          call bzWriteLine(f,buf)
          write(buf,'(A,5E15.4)')'#', max_absP, max_rhoP , max_rhoN,min_mass(particleID),max_mass(particleID)
          call bzWriteLine(f,buf)
          do index_absP=0,numSteps_absP
             do index_mass=0,numSteps_mass
                do index_rhoN=0,numSteps_rhoN
                   write(buf,format) widthTable(index_absP,index_mass,index_rhoN,0:numSteps_rhoP,i)
                   call bzWriteLine(f,buf)
                end do
             end do
          end do
          call bzCloseW(f)
       end do
       write(*,*) '... finished particleID=',particleID
    end do
    write(*,*) '... finished tabulating IN-MEDIUM WIDTH table for BARYONS'
    write(*,*) line
    write(*,*)

  end subroutine tabulate_inMediumWidth_baryons



  !****************************************************************************
  !****s* inMediumWidth/tabulate_inMediumWidth_mesons
  ! NAME
  ! subroutine tabulate_inMediumWidth_mesons()
  !
  ! PURPOSE
  ! * Tabulates the in-medium width of mesons according to Gamma=Gamma_free+sigma*rho*v
  ! * An average over the Fermi sea is performed
  ! * Sigma is given by the actual full collision term, with the incoming particles in the vacuum
  ! * Simplification: An average over the charge of the meson is performed.
  !
  ! INPUTS
  ! * NONE
  !
  ! OUTPUT
  ! * Files "inMediumWidth/InMediumWidth.particleID.dat.bz2"
  !****************************************************************************
  subroutine tabulate_inMediumWidth_mesons()
    use idTable, only: nmes,pion
    use inputGeneral, only: path_To_Input
    use output, only: intToChar,line
    use mesonWidthMedium_tables, only: numSteps_rhoN,numSteps_rhoP,numSteps_absP,numsteps_mass,max_rhoP,max_rhoN
    use bzip

    implicit none

    real, parameter :: min_mass=0.01
    real, parameter :: max_mass=3.0
    real ::  delta_rhoN,delta_rhoP,delta_absP,delta_mass
    real, dimension(:,:,:,:), allocatable :: widthTable
    integer :: index_absP,index_mass,index_rhoN,index_rhoP,particleID
    character(20) :: format
    character(200) :: fileName
    real :: rhoN, rhoP, absP, mass
    type(bzFile) :: f
    character(len=100) :: buf

    if (readInput_flag) then
       call readinput
       readInput_flag=.false.

       write(*,*) 'tabulation parameters:'
       write(*,*) 'numSteps_rhoN= ',numSteps_rhoN
       write(*,*) 'numSteps_rhoP= ',numSteps_rhoP
       write(*,*) 'numSteps_absP= ',numSteps_absP
       write(*,*) 'numSteps_mass= ',numSteps_mass
       write(*,*) 'max_rhoP= ',max_rhoP
       write(*,*) 'max_rhoN= ',max_rhoN
    end if

    allocate(widthTable(0:numSteps_absP,0:numSteps_mass,0:numSteps_rhoN,0:numSteps_rhoP))
    delta_rhoN=max_rhoN/float(numSteps_rhoN)
    delta_rhoP=max_rhoP/float(numSteps_rhoP)
    delta_absP=max_absP_Mes/float(numSteps_absP)
    delta_mass=(max_mass-min_mass)/float(numSteps_mass)

    format='(' // intToChar(numSteps_rhoP+1) // 'E12.4)'

    write(*,*) 'Tabulating IN-MEDIUM WIDTH table for MESONS'
    do particleId=max(pion,minMes),min(pion+nmes-1,maxMes)
       write(*,*) 'particleID=',particleID
       do index_absP=0,numSteps_absP
          absP=float(index_absP)*delta_absP
          write(*,'(A,F15.4,2(A,I8))') 'absP=', absP,' steps:',index_absP,'/',numsteps_absP
          do index_mass=0,numSteps_mass
             mass=min_mass+float(index_mass)*delta_mass
             write(*,'(A,F15.4,2(A,I8))') 'mass=', mass,' steps:',index_Mass,'/',numsteps_mass
             do index_rhoN=0,numSteps_rhoN
                rhoN=float(index_rhoN)*delta_rhoN
                do index_rhoP=0,numSteps_rhoP
                   rhoP=float(index_rhoP)*delta_rhoP
                   widthTable(index_absP,index_mass,index_rhoN,index_rhoP)=&
                            & evaluateCollisionBroadening_mesons(particleID,absP,mass,rhoN,rhoP)
                end do
             end do
             write(*,*) '     w: ',widthTable(index_absP,index_mass,1,1),'...', &
                        widthTable(index_absP,index_mass,numSteps_rhoN,numSteps_rhoP)
          end do
       end do

       ! Write the table to file:
       if (writelocal) then
          fileName='./InMediumWidth.'//trim(intToChar(particleID))//'.dat.bz2'
       else
          fileName=trim(path_to_Input)//'/inMediumWidth/InMediumWidth.'//trim(intToChar(particleID))//'.dat.bz2'
       end if
       f = bzOpenW(trim(fileName))
       write(buf,'(A,7I6)')'#', numSteps_absP,numSteps_mass,numSteps_rhoN,numSteps_rhoP,num_MonteCarlo_Points_mesons,cmin,cmax
       call bzWriteLine(f,buf)
       write(buf,'(A,5E15.4)')'#', max_absP_Mes, max_rhoP , max_rhoN,min_mass,max_mass
       call bzWriteLine(f,buf)
       do index_absP=0,numSteps_absP
          do index_mass=0,numSteps_mass
             do index_rhoN=0,numSteps_rhoN
                write(buf,format) widthTable(index_absP,index_mass,index_rhoN,0:numSteps_rhoP)
                call bzWriteLine(f,buf)
             end do
          end do
       end do
       call bzCloseW(f)

       write(*,*) '... finished particleID=',particleID
    end do
    write(*,*) '... finished tabulating IN-MEDIUM WIDTH table for MESONS'
    write(*,*) line
    write(*,*)

    deallocate(widthTable)

  end subroutine tabulate_inMediumWidth_mesons



  !****************************************************************************
  !****f* inMediumWidth/get_pauliBlockedDecayWidth
  ! NAME
  ! function get_pauliBlockedDecayWidth(particleID,momentumLRF,mass,rhoN,rhoP) result(gammaDecay)
  !
  ! PURPOSE
  ! * Returns the in-medium width of baryons according to  Gamma=Gamma_free*Pauli_Blocking
  ! * An average over the Fermi sea is performed
  ! * Simplification: An average over the charge of the baryon is performed.
  !
  ! INPUTS
  ! * integer, intent(in)              :: particleID    -- ID of baryon
  ! * real  , intent(in)               :: absP          -- absolute Momentum
  ! * real ,intent(in)                 :: mass          -- Mass of baryon
  ! * real, intent (in)                :: rhoN,rhoP     -- proton and neutron density in fm^-3
  !
  ! OUTPUT
  ! * real :: gammaDecay
  !****************************************************************************
  function get_pauliBlockedDecayWidth(particleID,momLRF,mass,rhoN,rhoP) result(gammaDecay)
    !****f*
    ! NAME
    ! PURPOSE
    ! NOTES
    !**************************************************************************
    use random
    use master_1Body, only: decayParticle
    use particleProperties, only: validcharge_id
    use idTable, only: nucleon
    use particleDefinition
    use baryonWidthMedium_tables, only: num_MonteCarlo_Points,get_min_charge_loop,get_max_charge_loop

    implicit none


    integer, intent(in)              :: particleID
    real, intent(in)                 :: momLRF
    real ,intent(in)                 :: mass

    real :: gammaDecay

    type(particle)                   :: part
    type(particle),dimension(1:10)   :: finalState


    logical:: pauliBlocking

    logical :: collisionFlag, pauliFlag, finalFlag
    real :: rhoN,rhoP ! neutron and proton density
    real :: gamma, fermi_n, fermi_p

    integer :: i,j

    call setToDefault(part)
    part%ID=particleID
    part%mass=mass
    part%momentum(1:2)= 0
    part%momentum(3)=momLRF
    part%momentum(0) = freeEnergy(part)
    part%velocity=part%momentum(1:3)/part%momentum(0)
    part%position=999.   ! outside nucleus

    fermi_p=fermiMom(rhoP)
    fermi_n=fermiMom(rhoN)
    finalFlag=.true.

    gammaDecay=0.
    monteCarloLoop :  do i=1,num_MonteCarlo_Points

       !***********************************************************************
       ! PRELIMINARY : AVERAGING OVER CHARGE
       do
          part%charge=NINT(float(get_min_charge_loop())+float(get_max_charge_loop()-get_min_charge_loop())*rn())
          if (validCharge_ID(part%ID,part%charge)) exit
       end do
       !***********************************************************************

       finalState%ID=0

       call decayParticle(part,finalState,collisionFlag,pauliFlag,finalFlag,0.,gamma)
       if (.not.collisionFlag) cycle monteCarloLoop

       ! Check Pauli-Blocking
       if (.not.pauliFlag) then
          pauliBlocking=.false.
          pauliLoop : do j=lbound(finalState,dim=1),ubound(finalState,dim=1)
             if (finalState(j)%ID.eq.nucleon.and..not.finalstate(j)%antiparticle) then
                if (finalState(j)%charge.eq.0) then
                   if (AbsMom(finalState(j)).lt.fermi_n) then
                      pauliBlocking=.true.
                      exit pauliLoop
                   end if
                else if (finalState(j)%charge.eq.1) then
                   if (AbsMom(finalState(j)).lt.fermi_p) then
                      pauliBlocking=.true.
                      exit pauliLoop
                   end if
                else
                   write(*,*) 'error in inMediumWidth.charge.', finalstate(j)
                   stop 'error in inMediumWidth.charge.'
                end if
             end if
          end do pauliLoop

          if (pauliBlocking) cycle monteCarloLoop
       end if
       gammaDecay=gammaDecay+gamma
    end do monteCarloLoop
    gammaDecay=gammaDecay/float(num_MonteCarlo_Points)
  end function get_pauliBlockedDecayWidth




  !****************************************************************************
  !****f* inMediumWidth/evaluateCollisionBroadening_baryons
  ! NAME
  ! function evaluateCollisionBroadening_baryons(particleID,momentumLRF,mass,rhoN,rhoP) RESULT(gcoll)
  !
  ! PURPOSE
  ! * Returns the in-medium width of baryons according to Gamma=sigma*rho*v *Pauli_Blocking
  ! * An average over the Fermi sea is performed
  ! * Sigma is given by the actual full collision term, with the incoming particles in the vacuum
  ! * Simplification: An average over the charge of the baryon is performed.
  !
  ! INPUTS
  ! * integer, intent(in)              :: particleID    -- ID of baryon
  ! * real  , intent(in)               :: absP          -- absolute Momentum
  ! * real ,intent(in)                 :: mass          -- Mass of baryon
  ! * real, intent (in)                :: rhoN,rhoP     -- proton and neutron density in GeV^3
  ! * integer, optional, intent (in)   :: monte_in     -- if this input is given, then the number of Monte Carlo points is chosen
  !   according to this input; otherwise it is equal to "num_MonteCarlo_Points".
  !
  ! OUTPUT
  ! * real :: gcoll
  ! * real , optional,intent(out)      :: gcoll_elastic_out -- elastic collisional width
  !****************************************************************************
  function evaluateCollisionBroadening_baryons(particleID,momLRF,mass,rhoN,rhoP,monte_in,gcoll_elastic_out) RESULT(gcoll)
    use particleDefinition
    use random
    use constants, only: GeVSquared_times_mb, mN
    use particleProperties, only: validCharge_ID
    use master_2Body, only: generateFinalState
    use IDTable, only: nucleon
    use baryonWidthMedium_tables, only: get_min_charge_loop,get_max_charge_loop

    implicit none

    integer, intent(in)              :: particleID
    real,    intent(in)              :: momLRF
    real ,intent(in)                 :: mass
    real, intent(in)                 :: rhoN,rhoP
    integer, optional ,intent(in)    :: monte_in ! Input for Monte carlo points
    real , optional,intent(out)      :: gcoll_elastic_out ! elastic collisional width
    real :: gColl

    type(particle)                   :: part,nuc

    integer :: i,j,charge,num_MonteCarlo_Points_local

    real    :: stringFactor
    integer :: numEnsembles
    type(particle), dimension(1:2)  :: pair
    type(particle), dimension(1:100) :: finalState
    real    :: time
    logical :: collisionFlag
    logical :: HiEnergyFlag  ! .true. if fritiof was used
    integer :: HiEnergyType  ! 0:LowEnergy, 1:Fritiof, 2:Pythia
    real    :: sigmaTot

    real :: fermiMomentum, pnuc,cost,vrel, fermi_p, fermi_n,rho

    real :: gammaLorentz

    logical :: pauliBlocking
    logical :: usePauli
    logical,parameter :: vac_check=.false. ! only for debugging!!!

    !real ::   imsig2,imsig3,imsigq,absP
    real :: gcoll_elastic
    logical :: pauliIsUsedforXsection!,elastic
    real :: srts
    integer :: pauliblocked

    if (readInput_flag) then
       call readinput
       readInput_flag=.false.
    end if

    if (present(monte_in)) then
       num_MonteCarlo_Points_local=monte_in
    else
       num_MonteCarlo_Points_local= 250 !for speed up, was "num_MonteCarlo_Points"
    end if

    stringFactor=1.
    numEnsembles=1
    time=999.

    call setToDefault(part)
    part%ID=particleID

    if (particleID.eq.nucleon) then
       ! We don't know the off-shell cross sections. Therefore we
       ! assume that the nucleon width is independent of mass:
       part%mass=mN
    else
       part%mass=mass
    end if
    part%momentum(1:2)= 0
    part%momentum(3)=momLRF
    part%momentum(0) = freeEnergy(part)
    part%velocity=part%momentum(1:3)/part%momentum(0)
    part%position=999.   ! outside nucleus


    call setToDefault(nuc)
    nuc%ID=nucleon
    nuc%mass=mN
    nuc%position=part%position

    gColl=0.
    gColl_elastic=0.

    fermi_p=fermiMom(rhoP)
    fermi_n=fermiMom(rhoN)



    nucleonChargeLoop: do charge=0,1
       if (vac_check.and.charge.eq.0) cycle
       if (charge.eq.1) then
          fermiMomentum=fermi_p
          rho=rhoP
       else if (charge.eq.0) then
          fermiMomentum=fermi_n
          rho=rhoN
       end if
       if (rho.lt.1E-10) cycle

       nuc%charge=charge


       usePauli=.true.
       pauliblocked=0

       !check if srts>4 (much above smooth transition region) and might run into Pythia
       !if yes, reduce number of Monte-Carlo-points to 1 and neglect Pauli blocking
       !(only cross section needed)
       !take minimal srts: both momenta parallel

       nuc%momentum(1) = 0.
       nuc%momentum(2) = 0.
       nuc%momentum(3) = fermiMomentum
       nuc%momentum(0) = freeEnergy(nuc)
       srts=sqrtS(nuc,part)
       if (debugFlag) write(*,*) 'min srts', srts
       if (srts.gt.4.) then
          num_MonteCarlo_Points_local=1
          usePauli=.false.
          if (debugFlag) write(*,*) 'speedup for srts', srts
       end if


       monteCarloLoop :  do i=1,num_MonteCarlo_Points_local
          ! Setting up the incoming nucleon
          !  * momentum of incoming nucleon:

          !********************************************************************
          ! PRELIMINARY : AVERAGING OVER CHARGE
          do
             if (vac_check) then
                part%charge=1 ! Set resonance charge to 1 for debugging
             else
                part%charge=NINT(float(get_min_charge_loop())+float(get_max_charge_loop()-get_min_charge_loop())*rn())
             end if
             if (validCharge_ID(part%ID,part%charge)) exit
          end do
          !********************************************************************

          pnuc=(rn())**(1./3.)*fermiMomentum
          cost=(rn()-0.5)*2.
          nuc%momentum(1) = pnuc*sqrt(max(1.-cost**2,0.))
          nuc%momentum(2) = 0.
          nuc%momentum(3) = pnuc*cost

          if (vac_check) nuc%momentum(1:3) = 0.!set nuc momentum to 0 for debugging

          nuc%momentum(0) = freeEnergy(nuc)
          nuc%velocity    = nuc%momentum(1:3)/nuc%momentum(0)

          pair(1)=part
          pair(2)=nuc

          call setToDefault(finalState)
          finalState%ID=0

          call generateFinalState (pair, finalState, stringFactor, numEnsembles, time, collisionFlag, HiEnergyFlag, HiEnergyType, &
                                   sigTot_out=sigmaTot, pauliIncluded_out=pauliIsUsedforXsection)
          if (.not. collisionFlag) cycle monteCarloLoop

          ! Check Pauli-Blocking
          pauliBlocking=.false.
          if (usePauli.and.(.not.pauliIsUsedforXsection)) then
             pauliLoop : do j=lbound(finalState,dim=1),ubound(finalState,dim=1)

                if (finalState(j)%ID.eq.nucleon.and..not.finalstate(j)%antiparticle) then
                   if (finalState(j)%charge.eq.0) then
                      if (AbsMom(finalState(j)).lt.fermi_n) then
                         pauliBlocking=.true.
                         exit pauliLoop
                      end if
                   else if (finalState(j)%charge.eq.1) then
                      if (AbsMom(finalState(j)).lt.fermi_p) then
                         pauliBlocking=.true.
                         exit pauliLoop
                      end if
                   else
                      write(*,*) 'error in inMediumWidth.charge.',finalState(j)
                      stop 'error in inMediumWidth.charge.'
                   end if
                end if
             end do pauliLoop
             if (pauliBlocking) then
                pauliblocked=pauliblocked+1
                if (debugFlag) write(*,*) 'Pauli blocked'
                cycle monteCarloLoop
             end if
          end if
          ! Evaluate relative velocity:
          vrel=sqrt(Dot_Product(nuc%velocity-part%velocity,nuc%velocity-part%velocity))


          if (debugFlag) then
             write(*,*) 'In:',pair%ID
             write(*,*) 'OUT:',finalState(1:3)%ID
             write(*,*) 'sigma:',sigmaTot
             write(*,*) 'vrel:',vrel
             write(*,*) 'rho:',rho
             write(*,*) 'flags:',collisionFlag,HiEnergyFlag,HiEnergyType
             write(*,*) 'Gamma:',vrel*rho*sigmaTot* GeVSquared_times_mb
             write(*,'(4G18.3)') absmom(nuc), absMom(part),momLRF
             write(*,'(4G18.3)') part%mass, nuc%mass
             write(*,'(4G18.3)') finalstate(1)%mass, finalstate(2)%mass
             write(*,'(4G18.3)') AbsMom(finalState(1)),AbsMom(finalState(2)),fermi_n,fermi_p
          end if

          ! Evaluate width in GeV:
          gColl=gColl+ vrel*rho*sigmaTot* GeVSquared_times_mb

          if (particleID.ne.nucleon) then
             if ((finalState(1)%ID.eq.nucleon.or.finalState(2)%ID.eq.nucleon)&
                  & .and.(finalState(1)%ID.eq.particleID.or.finalState(2)%ID.eq.particleID) &
                  & .and.(finalState(3)%ID.le.0) ) then
                ! Elastic event
                gColl_elastic=gColl_elastic+ vrel*rho*sigmaTot* GeVSquared_times_mb
             end if
          else
             if (finalState(1)%ID.eq.nucleon.and.finalState(2)%ID.eq.nucleon.and.finalState(3)%ID.le.0) then
                ! Elastic event
                gColl_elastic=gColl_elastic+ vrel*rho*sigmaTot* GeVSquared_times_mb
             end if
          end if

       end do monteCarloLoop

       if (debugFlag) write(*,*) 'Pauli blocked events', pauliblocked

    end do nucleonChargeLoop
    gColl=gColl/float(num_MonteCarlo_Points_local)
    gColl_elastic=gColl_elastic/float(num_MonteCarlo_Points_local)
    ! Transform the widht into the particles rest-frame
    gammaLorentz = 1./sqrt( 1. - dot_product(part%velocity(1:3),part%velocity(1:3)) )
    if (vac_check) gammaLorentz=1.
    gColl=gColl*gammaLorentz
    if (present(gcoll_elastic_out)) gcoll_elastic_out=gcoll_elastic*gammaLorentz

  end function evaluateCollisionBroadening_baryons



!******************************************************************************
  !****f* inMediumWidth/evaluateCollisionBroadening_mesons
  ! NAME
  ! function evaluateCollisionBroadening_mesons(particleID,momentumLRF,mass,rhoN,rhoP) RESULT(gcoll)
  !
  ! PURPOSE
  ! * Returns the in-medium width of mesons according to Gamma=sigma*rho*v
  ! * An average over the Fermi sea is performed
  ! * Sigma is given by the actual full collision term, with the incoming particles in the vacuum
  ! * Simplification: An average over the charge of the meson is performed.
  !
  ! INPUTS
  ! * integer, intent(in)              :: particleID    -- ID of particle
  ! * real  , intent(in)               :: absP          -- absolute Momentum
  ! * real ,intent(in)                 :: mass          -- Mass of meson
  ! * real, intent (in)                :: rhoN,rhoP     -- proton and neutron density in GeV^3
  !
  ! OUTPUT
  ! * real :: gcoll
  !****************************************************************************
  function evaluateCollisionBroadening_mesons(particleID,momLRF,mass,rhoN,rhoP) RESULT(gcoll)
    use particleDefinition
    use random
    use constants, only: GeVSquared_times_mb, mN
    use particleProperties, only: validCharge_ID
    use master_2Body, only: generateFinalState
    use IDTable, only: nucleon

    implicit none

    integer, intent(in)              :: particleID
    real,    intent(in)              :: momLRF
    real ,intent(in)                 :: mass
    real, intent(in)                 :: rhoN,rhoP
    real :: gColl

    type(particle)                   :: part,nuc

    integer :: i,j,charge,num_MonteCarlo_Points_local

    real    :: stringFactor=1.
    integer :: numEnsembles=1
    type(particle), dimension(1:2)  :: pair
    type(particle), dimension(1:100) :: finalState
    real    :: time=999.
    logical :: collisionFlag
    logical :: HiEnergyFlag  ! .true. if fritiof was used
    integer :: HiEnergyType  ! 0:LowEnergy, 1:Fritiof, 2:Pythia
    real    :: sigmaTot

    real :: fermiMomentum, pnuc,cost,vrel, fermi_p, fermi_n,rho

    real :: gammaLorentz

    logical :: pauliBlocking
    logical :: usePauli

    logical :: pauliIsUsedforXsection
    real :: srts

    if (readInput_flag) then
       call readinput
       readInput_flag=.false.
    end if

    num_MonteCarlo_Points_local = num_MonteCarlo_Points_mesons

    call setToDefault(part)
    part%ID=particleID

    part%mass=mass
    part%momentum(1:2)= 0
    part%momentum(3)=momLRF
    part%momentum(0) = freeEnergy(part)
    part%velocity=part%momentum(1:3)/part%momentum(0)
    part%position=999.   ! outside nucleus

    call setToDefault(nuc)
    nuc%ID=nucleon
    nuc%mass=mN
    nuc%position=part%position

    gColl=0.

    fermi_p=fermiMom(rhoP)
    fermi_n=fermiMom(rhoN)

    nucleonChargeLoop: do charge=0,1
       if (charge.eq.1) then
          fermiMomentum=fermi_p
          rho=rhoP
       else if (charge.eq.0) then
          fermiMomentum=fermi_n
          rho=rhoN
       end if
       if (rho.lt.1E-10) cycle

       nuc%charge=charge

       usePauli=.true.

       !check if srts>4 (much above smooth transition region) and might run into Pythia
       !if yes, reduce number of Monte-Carlo-points to 1 and neglect Pauli blocking
       !(only cross section needed)
       !take minimal srts: both momenta parallel

       nuc%momentum(1) = 0.
       nuc%momentum(2) = 0.
       nuc%momentum(3) = fermiMomentum
       nuc%momentum(0) = freeEnergy(nuc)
       srts=sqrtS(nuc,part)
       if (debugFlag) write(*,*) 'min srts', srts
       if (srts>4.) then
          num_MonteCarlo_Points_local=1
          usePauli=.false.
          if (debugFlag) write(*,*) 'speedup for srts', srts
       end if


       monteCarloLoop :  do i=1,num_MonteCarlo_Points_local

          !********************************************************************
          ! PRELIMINARY : AVERAGING OVER CHARGE
          do
             part%charge=NINT(float(cmin)+float(cmax-cmin)*rn())
             if (validCharge_ID(part%ID,part%charge)) exit
          end do
          !********************************************************************

          ! Setting up the incoming nucleon
          !  * momentum of incoming nucleon:
          pnuc=(rn())**(1./3.)*fermiMomentum
          cost=(rn()-0.5)*2.
          nuc%momentum(1) = pnuc*sqrt(max(1.-cost**2,0.))
          nuc%momentum(2) = 0.
          nuc%momentum(3) = pnuc*cost
          nuc%momentum(0) = freeEnergy(nuc)
          nuc%velocity    = nuc%momentum(1:3)/nuc%momentum(0)

          pair(1)=part
          pair(2)=nuc
          call setToDefault(finalState)
          finalState%ID=0
          call generateFinalState (pair, finalState, stringFactor, numEnsembles, time, collisionFlag, HiEnergyFlag, HiEnergyType, &
                                   sigTot_out=sigmaTot, pauliIncluded_out=pauliIsUsedforXsection)
          if (.not.collisionFlag) cycle monteCarloLoop

          ! Check Pauli-Blocking
          pauliBlocking=.false.
          if (usePauli.and.(.not.pauliIsUsedforXsection)) then
             pauliLoop : do j=lbound(finalState,dim=1),ubound(finalState,dim=1)
                if (finalState(j)%ID.eq.nucleon.and..not.finalstate(j)%antiparticle) then
                   if (finalState(j)%charge==0 .and. AbsMom(finalState(j))<fermi_n) then
                      pauliBlocking=.true.
                      exit pauliLoop
                   else if (finalState(j)%charge==1 .and. AbsMom(finalState(j))<fermi_p) then
                      pauliBlocking=.true.
                      exit pauliLoop
                   end if
                end if
             end do pauliLoop
             if (pauliBlocking) then
                cycle monteCarloLoop
             end if
          end if
          ! Evaluate relative velocity:
          vrel=sqrt(Dot_Product(nuc%velocity-part%velocity,nuc%velocity-part%velocity))

          if (debugFlag) then
             write(*,*) 'In:',pair%ID
             write(*,*) 'OUT:',finalState(1:3)%ID
             write(*,*) 'sigma:',sigmaTot
             write(*,*) 'vrel:',vrel
             write(*,*) 'rho:',rho
             write(*,*) 'flags:',collisionFlag,HiEnergyFlag,HiEnergyType
             write(*,*) 'Gamma:',vrel*rho*sigmaTot* GeVSquared_times_mb
             write(*,'(4G18.3)') absmom(nuc), absMom(part),momLRF
             write(*,'(4G18.3)') part%mass, nuc%mass
             write(*,'(4G18.3)') finalstate(1)%mass, finalstate(2)%mass
             write(*,'(4G18.3)') AbsMom(finalState(1)),AbsMom(finalState(2)),fermi_n,fermi_p
          end if

          ! Evaluate width in GeV:
          gColl=gColl+ vrel*rho*sigmaTot* GeVSquared_times_mb

       end do monteCarloLoop

    end do nucleonChargeLoop
    gammaLorentz = 1./sqrt( 1. - dot_product(part%velocity(1:3),part%velocity(1:3)) )
    gColl=gColl*gammaLorentz/float(num_MonteCarlo_Points_local)

  end function evaluateCollisionBroadening_mesons



  !****************************************************************************
  !****f* inMediumWidth/fermiMom
  ! NAME
  ! real function fermiMom(rho)
  !
  ! PURPOSE
  ! * Returns the fermi momentum.
  !
  ! INPUTS
  ! * real, intent (in)                :: rho     -- proton or neutron density in GeV^-3
  !
  ! OUTPUT
  ! * Fermi momentum in GeV
  !****************************************************************************
  real function fermiMom(rho)
    use constants, only: pi
    implicit none
    real :: rho
    fermiMom=(3.*pi**2*rho)**(1./3.)
  end function fermiMom

end module inMediumWidth
