!***************************************************************************
!****m* /InputCoalescence
! NAME
! module InputCoalescence
!
! PURPOSE
! This is the main module for reading the parameters/switches concerning 
! the phase-space coalescence.
!***************************************************************************
Module InputCoalescence

  !*************************************************************************
  !****g* InputCoalescence/Get_Analysis
  ! SOURCE
  !
  Logical,save :: Get_Analysis = .false.
  !
  ! PURPOSE
  ! Switch for Analysis of Coalescence-Events
  !
  !*************************************************************************

  !*************************************************************************
  !****g* InputCoalescence/Get_Flow
  ! SOURCE
  !
  Logical,save :: Get_Flow = .false.
  !
  ! PURPOSE
  ! Switch for calculate Flow-Observables
  !
  !*************************************************************************

  !*************************************************************************
  !****g* InputCoalescence/Get_Zdist
  ! SOURCE
  !
  Logical,save :: Get_Zdist = .false.
  !
  ! PURPOSE
  ! Switch for calculate ChargeDistributions
  !
  !*************************************************************************

  !*************************************************************************
  !****g* InputCoalescence/Get_Model
  ! SOURCE
  !
  Integer,save :: Get_Model = 0
  !
  ! PURPOSE
  ! Switch for Mean-Field Model (1=Skyrme-MDI, 2=RMF)
  ! NOTES 
  ! Precently only RMF has been fully implemented. This is only important 
  ! in calculating the binding energies of fragments.
  !
  !*************************************************************************

  !*************************************************************************
  !****g* InputCoalescence/Get_Asy
  ! SOURCE
  !
  Logical,save :: Get_Asy = .true.
  !
  ! PURPOSE
  ! Switch for checking N/Z-Ratio of produced clusters
  ! NOTES
  ! Should always be true!!!
  !
  !*************************************************************************

  !*************************************************************************
  !****g* InputCoalescence/Get_GEM
  ! SOURCE
  !
  Logical,save :: Get_GEM = .false.
  !
  ! PURPOSE
  ! Switch for statistical de-excitation
  ! NOTES
  ! Generalized Evaporation Model (S. Furihata). 
  ! Fermi break-up not included!!
  !
  !*************************************************************************

  !*************************************************************************
  !****g* InputCoalescence/MaxGEM
  ! SOURCE
  !
  Integer,save :: MaxGEM = 50
  !
  ! PURPOSE
  ! Max. steps for Evaporation
  ! NOTES
  ! * Evaporation is repeated max. MaxGEM times, until all fragments 
  !   are stable. 
  ! * MaxGEM=50 is a good estimation. Usualy only 10-20 max. steps are needed.
  !
  !*************************************************************************

  !*************************************************************************
  !****g* InputCoalescence/Get_Asy
  ! SOURCE
  !
  Logical,save :: Get_Output = .false.
  !
  ! PURPOSE
  ! Switch for write down the FragmentVector into an output file
  ! NOTES
  ! needed only for HypHI- and FOPI-Collaborations.
  !
  !*************************************************************************

  !*************************************************************************
  !****g* InputCoalescence/R_c
  ! SOURCE
  !
  Real,save :: R_c = 3.5
  !
  ! PURPOSE
  ! Coalescence Parameter in Coordinate Space (fm)
  !
  !*************************************************************************

  !*************************************************************************
  !****g* InputCoalescence/P_c
  ! SOURCE
  !
  Real,save :: P_c = 1.3
  !
  ! PURPOSE
  ! Coalescence Parameter in Momentum Space (1/fm)
  !
  !*************************************************************************

  !*************************************************************************
  !****g* InputCoalescence/Out1
  ! SOURCE
  !
  Character(60),save :: Out1='Flows'
  !
  ! PURPOSE
  ! Output File (contains flow observables)
  !
  !*************************************************************************

  !*************************************************************************
  !****g* InputCoalescence/Out2
  ! SOURCE
  !
  Character(60),save :: Out2='Pt-Spectra'
  !
  ! PURPOSE
  ! Output File (contains flow observables)
  !
  !*************************************************************************
  
  !*************************************************************************
  !****g* InputCoalescence/Out3
  ! SOURCE
  !
  Character(60),save :: Out3='ChargeDist'
  !
  ! PURPOSE
  ! Output File (contains Charge Distributions)
  !
  !*************************************************************************

  !*************************************************************************
  !****g* InputCoalescence/Out4
  ! SOURCE
  !
  Character(60),save :: Out4='TeilchenVector'
  !
  ! PURPOSE
  ! Output File (contains the fragment vector)
  !
  !*************************************************************************

  !*************************************************************************
  !****g* InputCoalescence/Out5
  ! SOURCE
  !
  Character(60),save :: Out5='dNdY'
  !
  ! PURPOSE
  ! Output File (Rapidity spectra for various particles)
  !
  !*************************************************************************

  !*************************************************************************
  !****g* InputCoalescence/Out6
  ! SOURCE
  !
  Character(60),save :: Out6='dNdY'
  !
  ! PURPOSE
  ! Output File (Rapidity spectra for various particles)
  !
  !*************************************************************************

  !*************************************************************************
  !****g* InputCoalescence/Out7
  ! SOURCE
  !
  Character(60),save :: Out7='dNdY'
  !
  ! PURPOSE
  ! Output File (Rapidity spectra for various hyperclusters)
  !
  !*************************************************************************

  !*************************************************************************
  !****g* InputCoalescence/PathToCLusterInput
  ! SOURCE
  !
  character(100), save :: PathToCLusterInput
  !
  ! PURPOSE
  ! Link to the input-directory for the evaporation part (see notes above).
  !
  !*************************************************************************


!==========================================================================*
CONTAINS
!==========================================================================*

  !*************************************************************************
  !Read input-parameters ***************************************************
  !*************************************************************************
  subroutine Get_InputCoalescence
    use InputGeneral, only : RealParticles,SubEvents,NumEnsemples,& 
         &                   pathToBUUInput,BUU_DataFile,& 
         &                   Get_Hyp
    implicit none
    integer :: ios
    !-----------------------------------------------------------------------
    NAMELIST /inputCoalescence/Get_Analysis,Get_Flow,Get_Zdist,& 
         & Get_Model, & 
         & Get_Asy,Get_GEM,Get_Output,R_c,P_c,  & 
         & MaxGEM,BUU_DataFile,Out1,Out2,Out3,Out4,Out5,Out6,Out7,&
         & PathToClusterInput

    rewind(5)
    read(5,nml=inputCoalescence,IOSTAT=ios)
    !-----------------------------------------------------------------------
    P_c = P_c*0.19733 !conversion from fm**-1 --> GeV !!!!!!!!!!
    !-----------------------------------------------------------------------
    if(ios.ne.0) then
       write(*,*) 'Error in  namelist "inputCoalescence" : This namelist is crucial. STOP!'
       stop
    end if
    !-----------------------------------------------------------------------

999  format(//,78('='),/,20('='),'  START LISTING OF CONTROL PARAMETERS ',20('='),/,78('='))
1000 format('  * Number of subsequent runs              = ',i5)
1001 format('  * Number of ensembles                    = ',i5)
1002 format('  * Length of ParticleVector               = ',i5)
1003 format('  * Coalescence-Parameter                  = ',f8.4,' [fm]  in coordinate space')
1004 format('  * Coalescence-Parameter                  = ',f8.4,' [GeV] in  momentum  space')
1007 format('  * Switch for Asymmetry Check             = ',l3)
1011 format('  * Switch for statistical GEM-Model decay = ',l3)
1009 format('  * Switch for Hyperon-Cluster Coalescence = ',l3)
1010 format('  * Switch for Analysis(Flows, spectra...) = ',l3)
1012 format('  * GiBUU file to be analyzed              : ',/, & 
          & 4x,a)
9999 format(78('='),/,20('='),'   END LISTING OF CONTROL PARAMETERS  ',20('='),/,78('='),//)


    write(*,999)
    write(*,1000) SubEvents
    write(*,1001) NumEnsemples
    write(*,1002) RealParticles
    write(*,1003) R_c
    write(*,1004) P_c
    if (.not.Get_Asy) then
       write(*,1007) Get_Asy
       write(*,*) '  ==> No check on N/Z-asymmetry!'
       write(*,*) '      !!!Be sure what you are doing by setting Get_Asy = .false.'
       write(*,*) '      !!!You will obtain in this case always strange clusters!!!'
    else
       write(*,1007) Get_Asy
       write(*,*) '    ==> You have chosen check for Asymmetry'
       write(*,*) '    ==> For clusters with unrealistic asymmetry'
       write(*,*) '        the nearest N/Z-stable element in the'
       write(*,*) '        NNDC-Table is found'
    endif
    if (.not.Get_GEM) then
       write(*,1011) Get_GEM
       write(*,*) '    ==> No use of Statistical Evaporation Model (S. Furihata)'
    else
       write(*,1011) Get_GEM
       write(*,*) '    ==> Statistical Evaporation (GEM) of hot clusters'
    endif
    if (.not.Get_Hyp) then
       write(*,1009) Get_Hyp
       write(*,*) '    ==> No Hyperon-Cluster Coalescence'
    else
       write(*,1009) Get_Hyp
       write(*,*) '    ==> With Hyperon-Cluster Coalescence'
    endif
    if (.not.Get_Analysis) then
       write(*,1010) Get_Analysis
       write(*,*) '    ==>No Analysis'
    else
       write(*,1010) Get_Analysis
       write(*,*) '    ==> Switch for Collective Flows  = ',Get_Flow
       write(*,*) '    ==> Switch for A- and Z-Distr.   = ',Get_ZDist
       write(*,*) '    ==> The boost parameters are initialized in ClusterAnalysis module'
    endif
    if (Get_Model==1) then
       write(*,*) ' * Non-Relativistic Mean-Field, Skyrme-soft'
    else if (Get_Model==2) then
       write(*,*) ' * Relativistic Mean-Field, Non-Linear Walecka Model'
    else 
       write(*,*) '   * Get_Model = ',Get_Model
       write(*,*) '   * This option for Mean-field is not defined!!!'
       write(*,*) '     Program will terminate NOW!!!'
       STOP
    endif

    write(*,1012) trim(PathToBUUInput)// trim(BUU_DataFile)

    write(*,9999)

  !*************************************************************************
  end subroutine Get_InputCoalescence !*************************************
  !*************************************************************************


end module InputCoalescence
