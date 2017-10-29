!***************************************************************************
!****m* /InputSMM
! NAME
! module InputSMM

! FUNCTION
! Defines the major input-parameters  and reads them from a JobCard.
! NOTES
! ALL INPUT PARAMETERS FOR THE SOURCES ARE EXTRACTED FROM DYNAMICAL 
! GiBUU-CALCULATIONS.
!***************************************************************************

Module InputSMM

  !*************************************************************************
  !****n* inputSMM/inputSMM
  ! NAME
  ! NAMELIST /inputSMM/
  ! PURPOSE
  ! Includes the input switches:
  ! * SMM_Seed
  ! * EventType
  ! * SMM_Flag
  ! * printFV_Flag
  ! * CorrectExc
  ! * Delta_Exc
  ! * E_Bind_Input
  ! * Hysto_Flag
  ! * HypCoalaMethod
  ! * RadiusHypC, RadiusHypP
  ! * GetExcBalance
  ! * beamID
  ! * beamEnergy
  ! * A_init
  ! * MaxNumSources
  ! * SourceInfo
  ! * PathToSMMInput
  !*************************************************************************


  !*************************************************************************
  !****g* inputSMM/SMM_Seed
  ! SOURCE
  !
  Integer,save :: SMM_Seed = 12345
  !
  ! PURPOSE
  ! Initialization for Random-Number Generator used in SMM-code.
  !
  !*************************************************************************  

  !*************************************************************************
  !****g* inputSMM/SMM_Events
  ! SOURCE
  !
  Integer,save :: SMM_Events = 50
  !
  ! PURPOSE
  ! MC-Events for the statistical multifragmentation code.
  !
  !*************************************************************************  

  !*************************************************************************
  !****g* inputSMM/SMM_Flag
  ! SOURCE
  !
  Integer,save :: SMM_Flag = 1
  !
  ! PURPOSE
  ! * Two options for SMM:
  ! * (1) SMM_Flag = 0 : Only evaporation
  ! * (2) Otherwise    : Fermi-break up + evaporation 
  ! * Case (2) is the most realistic one. 
  !
  !*************************************************************************  

  !*************************************************************************
  !****g* inputSMM/EventType
  ! SOURCE
  !
  Integer,save :: EventType = 1
  !
  ! PURPOSE
  ! =1 : Heavy-Ion Collision / =2: Hadron-Nucleus collision.
  !
  ! NOTES
  ! * In Case of EventType==2 there is now the possibility to test the 
  !   stability of an initialized nucleus in its BUU ground state at 
  !   different times. 
  ! * In the ideal case of a perfect initialization the SMM-code would 
  !   give only one "fragment", i.e. the initialized nucleus. However, 
  !   due to spurious fluctuations the initialized nucleus is excited 
  !   and thus the SMM-code let it dissociate! 
  !*************************************************************************

  !*************************************************************************
  !****g* inputSMM/MaxNumSources
  ! SOURCE
  !
  Integer,save :: MaxNumSources = 10
  !
  ! PURPOSE
  ! Max. number of sources, needed to allocate the "type(quelle) TheSource".
  !
  !*************************************************************************

  !*************************************************************************
  !****g* inputSMM/SourceInfo
  ! SOURCE
  !
  character(50), save :: SourceInfo='SourceFile.dat'
  !
  ! PURPOSE
  ! Properties of target/Fireball/Spectator sources are stored into these file.
  !
  !*************************************************************************

  !*************************************************************************
  !****g* inputSMM/PathToSMMInput
  ! SOURCE
  !
  character(100), save :: PathToSMMInput='~/'
  !
  ! PURPOSE
  ! Link to the input-directory for the source(s) files "SourceInfo"
  !
  !*************************************************************************

  !*************************************************************************
  !****g* inputSMM/printFV_Flag
  ! SOURCE
  !
  Logical, save :: printFV_Flag=.false.
  !
  ! PURPOSE
  ! If true, print the fragmentVector.
  ! NOTES
  ! Use this option only if neccessary, this output file occupies a lot of 
  ! memory (several 100 of MBs!!!)
  !
  !*************************************************************************

  !*************************************************************************
  !****g* inputSMM/CorrectExc
  ! SOURCE
  !
  Logical, save :: CorrectExc=.false.
  !
  ! PURPOSE
  ! If true, Correction factor for excitation energy used.
  ! NOTES
  ! This is due to the fact, that energy is not perfectly conserved. For 
  ! collisions with low excitation small changes in total energy (due to 
  ! non-conservation) may have significant influence on 
  ! statistical fragmentation (see also the notes on the 
  ! variable "EventType").
  !
  !*************************************************************************

  !*************************************************************************
  !****g* inputSMM/Delta_Exc
  ! SOURCE
  !
  Real, save :: Delta_Exc=0.0
  !
  ! PURPOSE
  ! Correction value for calculating excitation energy (in MeV per Nucleon).
  !
  !*************************************************************************

  !*************************************************************************
  !****g* inputSMM/E_Bind
  ! SOURCE
  !
  Real, save :: E_Bind_Input=1000000.
  !
  ! PURPOSE
  ! Binding energy per nucleon in units of MeV.
  ! NOTES
  ! * Used only for Hadron-Nucleus EventType.
  ! * The precise value of the binding energy is important for 
  !   proton-induced reactions with low excitation. In case 
  !   of heavy-ion collisions, we use approximate values given 
  !   in function BindingEnergy (Main_SMM.f90-module).
  ! * Please, be carefull which value you are using for E_Bind 
  !   (see notes below).
  ! * We need the value of the binding energy of the fragmenting nucleus, in 
  !   order to calculate its excitation energy (as the difference between the 
  !   total energy and the binding energy). Since the"BUU-ground-state" 
  !   properties of a source (given by A and Z) depend on technical parameters, 
  !   such as number of ensemples, gridsize of box etc, we have to explicitly 
  !   give here the "BUU" binding energy, as extracted from a ground-state 
  !   calculation.
  !
  !*************************************************************************

  !*************************************************************************
  !****g* inputSMM/Hysto_Flag
  ! SOURCE
  !
  Logical, save :: Hysto_Flag=.false.
  !
  ! PURPOSE
  ! If true, Z- and A-hystogramms are created (usefull only for checking 
  ! groundstate-calculations. 
  !
  !*************************************************************************

  !*************************************************************************
  !****g* inputSMM/RadiusHypC
  ! SOURCE
  !
  Real, save :: RadiusHypC=3.0  !estimation!!!
  !
  ! PURPOSE
  ! Coalescence radius in coordinate space [fm] (roughly estimated) 
  ! for coalescence between SMM fragments and strange baryons.
  !
  !*************************************************************************

  !*************************************************************************
  !****g* inputSMM/RadiusHypP
  ! SOURCE
  !
  Real, save :: RadiusHypP=1.4 !~Fermi momentum
  !
  ! PURPOSE
  ! Coalescence radius in momentum space [1/fm] (roughly estimated) 
  ! for coalescence between SMM fragments and strange baryons.
  !
  !*************************************************************************

  !*************************************************************************
  !****g* inputSMM/HypCoalaMethod
  ! SOURCE
  !
  Integer, save :: HypCoalaMethod=1
  !
  ! PURPOSE
  ! Method for coalescence between SMM fragments and strange baryons:
  ! HypCoalaMethod=1 ==> coalescence in coordinate & momentum space (default)
  ! HypCoalaMethod=2 ==> coalescence only in coordinate space
  !
  !*************************************************************************

  !*************************************************************************
  !****g* inputSMM/GetExcBalance
  ! SOURCE
  !
  Logical, save :: GetExcBalance=.true.
  !
  ! PURPOSE
  ! Switch for calculating excitation energy via energy balance (event-by-event).
  !
  !*************************************************************************

  !*************************************************************************
  !****g* inputSMM/beamID
  ! SOURCE
  !
  Integer, save :: beamID=1
  !
  ! PURPOSE
  ! ID of the beam particle. So far: (anti)proton, Xi (relevant for PANDA).
  !
  !*************************************************************************

  !*************************************************************************
  !****g* inputSMM/beamEnergy
  ! SOURCE
  !
  Real, save :: beamEnergy=0.0
  !
  ! PURPOSE
  ! kinetic energy of beam particle (LAB-frame), in units of GeV.
  !
  !*************************************************************************

  !*************************************************************************
  !****g* inputSMM/A_init
  ! SOURCE
  !
  Real, dimension(1:2), save :: A_init=0.0
  !
  ! PURPOSE
  ! Neutron (1) and proton (2) number of target nucleus.
  !
  !*************************************************************************

  !*************************************************************************
  !****g* inputSMM/XiTrigger
  ! SOURCE
  !
  Logical, save :: XiTrigger=.false.
  !
  ! PURPOSE
  ! If true, only events with one Xi are selected. Some observables (not everything)
  ! are then calculated for these "exclusive" events.
  !
  !*************************************************************************

!==========================================================================*
CONTAINS
!==========================================================================*

  !*************************************************************************
  !Read input-parameters ***************************************************
  !*************************************************************************
  subroutine Get_InputSMM
    use InputGeneral, only : SubEvents, NumEnsemples,Get_Hyp,ALADIN_Flag, & 
         & pathToBUUInput,BUU_DataFile
    implicit none
    integer :: ios
    !-----------------------------------------------------------------------
    NAMELIST /inputSMM/SMM_Seed,SMM_Events,EventType,SMM_Flag,printFV_Flag, &
         & CorrectExc,Delta_Exc,E_Bind_Input,Hysto_Flag, & 
         & HypCoalaMethod, RadiusHypC, RadiusHypP, & 
         & GetExcBalance, beamID, beamEnergy, A_init, XiTrigger, & 
         & MaxNumSources,SourceInfo,PathToSMMInput

    write(*,*) '----- Reading Namelist "inputSMM": Start -----'
    rewind(5)
    read(5,nml=inputSMM,IOSTAT=ios)
    !-----------------------------------------------------------------------
    if(ios.ne.0) then
       write(*,*) 'Error in  namelist "inputSMM" : This namelist is crucial. STOP!'
       stop
    end if
    RadiusHypP = RadiusHypP*0.19733 !Conversion 1/fm-->GeV
    Delta_Exc = Delta_Exc/1000. !Conversion MeV-->GeV
    E_Bind_Input = E_Bind_Input/1000.
    !-----------------------------------------------------------------------

999 format(//,78('='),/,20('='),'  START LISTING OF CONTROL PARAMETERS ',20('='),/,78('='))

997 format('  * SMM: Full Statistical (calls ZVEC1  in SMM-code)')
996 format('  * SMM: Only evaporation (calls EVANUC in SMM-code)')
998 format('  * Number of GiBUU Events to be analyzed = ',i5)
995 format('  * Number of SMM Events                  = ',i5)

1001 format('  * Hadron-Nucleus reactions: compound nucleus to be analyzed  :',1x,a)
1002 format('  * Heavy ion reactions: source file to be analyzed       :',1x,a)
1005 format('  * GiBUU file to be analyzed (only emitted particles)    :',1x,a)
1006 format('  * Source file to be analyzed                            :',1x,a)

1007 format('  * Value for Binding energy [GeV/A]: E_Bind_Input        :',1x,f8.4)
1008 format('  * Correction value for Excitation Energy [GeV/A]        :',1x,f8.4)

1009 format('  * Method for coalescence between fragments & Hyperons   :',1x,i2)
1010 format('  * Coalescence radius in coordinate space (fm)           :',1x,f8.4)
1011 format('  * Coalescence radius in momentum space (GeV)            :',1x,f8.4)
1012 format('  * Hypercluster formation is switched off')

9999 format(78('='),/,20('='),'   END LISTING OF CONTROL PARAMETERS  ',20('='),/,78('='),//)


    write(*,999)
    write(*,998) SubEvents*NumEnsemples
    write(*,995) SMM_Events

    if (SMM_Flag==0) then 
       write(*,996)
    else
       write(*,997)
    endif

    if (GetExcBalance .and. EventType==1) then
       write(*,*) 'Module inputSMM, routine Get_InputSMM:'
       write(*,*) 'Do not use energy balance in HIC-runs!!!'
       write(*,*) 'STOP'
       STOP
    endif

    write(*,*) 
    if (EventType==2 .and. .not.GetExcBalance) then 
       if (E_Bind_Input==1000.) then
          write(*,*) 'Module inputSMM, routine Get_InputSMM:'
          write(*,*) ' * For proton-induced reactions you must give'
          write(*,*) '   the precise value of the binding energy per nucleon'
          write(*,*) '   in your JobCard!!! '
          write(*,*) ' * Please run a BUU ground state with exact the same parameters'
          write(*,*) '   as for the proton-induced events analyzed here,'
          write(*,*) '   and extract the binding energy at t=0fm/c from the appropriate file'
          write(*,*) '!!! Termination of the ClusterCode !!!'
          STOP
       else
          write(*,1007) E_Bind_input
       endif
    endif
    write(*,*) ' * Correction for Excitation energy : ',CorrectExc
    if (CorrectExc) write(*,1008) Delta_Exc


    if (EventType==2) then
       write(*,1001) trim(PathToSMMInput)// SourceInfo
    else if (EventType==1) then
       write(*,1002) trim(PathToSMMInput)// SourceInfo
    else
       write(*,*) 'Module inputSMM, routine Get_InputSMM:'
       write(*,*) 'Wrong input for EventType: ',EventType
       write(*,*) '!!! Termination of the program !!!'
       STOP
    endif

    if ( .not. Get_Hyp ) then
       write(*,1012) 
    else
       write(*,1009) HypCoalaMethod
       write(*,1010) RadiusHypC
       write(*,1011) RadiusHypP
    endif

    if ( (EventType==1) .and. (.not.ALADIN_Flag) ) then
       write(*,*)
       write(*,*) 'Module inputSMM, routine Get_InputSMM:'
       write(*,*) ' * Fragmentation analysis of fireballs is still under construction...'
       write(*,*) '   STOP'
       write(*,*)
       STOP
    end if


    write(*,1005) trim(PathToBUUInput)// trim(BUU_DataFile)
    write(*,1006) trim(PathToSMMInput)// trim(SourceInfo)

    write(*,9999)

    write(*,*) '----- Reading Namelist "inputSMM": END   -----'
    write(*,*)

  !*************************************************************************
  end subroutine Get_InputSMM !*********************************************
  !*************************************************************************


end module InputSMM
