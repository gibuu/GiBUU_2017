!***************************************************************************
!****m* /fragmentationAnalysis
! NAME
! module fragmentationAnalysis

! FUNCTION
! This module contains some usefull analysis routines.
! NOTES
! This module uses all variables of the module InputForAnalysis.
!***************************************************************************

module fragmentationAnalysis

  use InputForAnalysis
  implicit none

  PRIVATE

  PUBLIC :: printFragmentVector, VariousDistributions, & 
       &    rapidityDistribution,spectatorFragmentation, & 
       &    charmsAnalysis, VelocityDistr_Charms, & 
       &    detailedHypAnalysis, producedParticles, & 
       &    ExDistributions

  real, parameter :: pi=3.141592654

  Logical, SAVE :: initParameters = .true.

contains

  !*************************************************************************
  !****s* FragmentationAnalysis/printFragmentVector
  ! FUNCTION
  ! Prints out the fragment vector.
  ! 
  !*************************************************************************
  subroutine printFragmentVector(isu,FragmentVector)
    use typeDefinitions, only : cluster

    integer, intent(in) :: isu
    type(cluster), dimension(:,:), intent(in) :: FragmentVector

    integer :: i,j,k

    open(103,file='auaufra.all',position='Append')
    Loop_over_ensembles: do i = 1,size(FragmentVector,dim=1)
       Loop_over_particles : do j = 1,size(FragmentVector,dim=2)
          If ( FragmentVector(i,j)%id == 0 ) then
             cycle Loop_over_particles
          end If
          write(103,50) FragmentVector(i,j)%ID,FragmentVector(i,j)%chargeNumber,&
               FragmentVector(i,j)%massNumber,FragmentVector(i,j)%HypNumber,&
               (FragmentVector(i,j)%position(k), k=1,3),&
               (FragmentVector(i,j)%momentum(k)/FragmentVector(i,j)%momentum(0), k=1,3),i,isu
       enddo Loop_over_particles
    enddo Loop_over_ensembles
    close(103)

50  format(i4,1x,i4,1x,i4,1x,i4,1x,3(1x,f8.3),3(1x,f8.3),1x,i4,1x,i4)

  !*************************************************************************
  end subroutine printFragmentVector !**************************************
  !*************************************************************************


  !*************************************************************************
  !****s* FragmentationAnalysis/detailedHypAnalysis
  ! FUNCTION
  ! Prints out detailed information on bound hyperons of all hyperfragments.
  !*************************************************************************
  subroutine detailedHypAnalysis(is,in,SubEvents,NumEnsemples,SMM_Events, &
                                 FragmentVector,particleVector,printflag)
    use typeDefinitions, only : cluster,particle

    integer, intent(in) :: is,in
    integer, intent(in) :: SubEvents,NumEnsemples,SMM_Events
    logical, intent(in) :: printflag
    type(cluster),  dimension(:,:), intent(in) :: FragmentVector
    type(particle), dimension(:,:), intent(in) :: particleVector

    integer :: i,j,HypIndex1
    integer, dimension(1:2) :: HypIndex2
    integer, save :: HyperCount=0, HyperCountSingle=0, HyperCountDouble=0, HyperCountAll=0

    open(103,file='HyperInfo.dat',position='Append')

    Loop_over_ensembles: do i = 1,size(FragmentVector,dim=1)
       Loop_over_particles : do j = 1,size(FragmentVector,dim=2)

          If ( FragmentVector(i,j)%id == 0 ) cycle Loop_over_particles
          If ( FragmentVector(i,j)%HypNumber == 0 ) cycle Loop_over_particles

          HyperCount = HyperCount + 1

          write(103,'(a,1x,i8)') 'subseqent run         = ',is
          write(103,'(a,1x,i8)') 'number of ens.        = ',in
          write(103,'(a,1x,i8)') 'number of SMM-Event   = ',i
          write(103,*)
          write(103,'(a,1x,i8)') 'Index of hypercluster = ',j
          write(103,*)  'stableFlag            = ',FragmentVector(i,j)%stableFlag
          write(103,'(a,1x,i8)') 'HypNumber             = ',FragmentVector(i,j)%HypNumber
          write(103,'(a,1x,a)')  'HypType               = ',FragmentVector(i,j)%HypType
          write(103,'(a,1x,i8)') 'Mass number           = ',FragmentVector(i,j)%MassNumber
          write(103,'(a,1x,i8)') 'Charge number         = ',FragmentVector(i,j)%ChargeNumber
          write(103,'(a,1x,3f12.4)') 'pos.                  = ',FragmentVector(i,j)%position
          write(103,'(a,1x,3f12.4)') 'Vel.                  = ',& 
               & FragmentVector(i,j)%momentum(1:3)/FragmentVector(i,j)%momentum(0)
          write(103,*)
          write(103,*) 'Counting index of Hyp.= ',HyperCount
          write(103,*)

          if (FragmentVector(i,j)%HypNumber==1) then

             HyperCountSingle = HyperCountSingle + 1

             HypIndex1 = FragmentVector(i,j)%HypContent(1)
             write(103,'(a,1x,i8)') 'Index of bound hyperon = ',HypIndex1
             write(103,'(a)') 'Information on bound Hyperon:'
             write(103,'(a,1x,i8)') 'ID       = ',particleVector(in,HypIndex1)%number
             write(103,'(a,1x,f8.4)') 'born     = ',particleVector(in,HypIndex1)%bornTime
             write(103,'(a,1x,f8.4)') 'lastColl = ',particleVector(in,HypIndex1)%lastCollTime
             write(103,'(a,1x,i10)')  'collHis  = ',particleVector(in,HypIndex1)%collHis
             write(103,'(a,1x,3f12.4)') 'Pos.   = ',particleVector(in,HypIndex1)%position
             write(103,'(a,1x,3f12.4)') 'Vel.   = ', & 
                  & particleVector(in,HypIndex1)%momentum(1:3)/particleVector(in,HypIndex1)%momentum(0)
          else

             HyperCountDouble = HyperCountDouble + 1

             HypIndex2(1:2) = FragmentVector(i,j)%HypContent(1:2)
             write(103,'(a,1x,i8)') 'Index of bound hyperon #1 = ',HypIndex2(1)
             write(103,'(a)') 'Information on bound Hyperons:'
             write(103,'(a,1x,i8)') 'ID       = ',particleVector(in,HypIndex2(1))%number
             write(103,'(a,1x,f8.4)') 'born     = ',particleVector(in,HypIndex2(1))%bornTime
             write(103,'(a,1x,f8.4)') 'lastColl = ',particleVector(in,HypIndex2(1))%lastCollTime
             write(103,'(a,1x,i10)') 'collHis  = ',particleVector(in,HypIndex2(1))%collHis
             write(103,'(a,1x,3f12.4)') 'Pos.   = ',particleVector(in,HypIndex2(1))%position
             write(103,'(a,1x,3f12.4)') 'Vel.   = ', & 
                  & particleVector(in,HypIndex2(1))%momentum(1:3)/particleVector(in,HypIndex2(1))%momentum(0)
             write(103,*)
             write(103,'(a,1x,i8)') 'Index of bound hyperon #2 = ',HypIndex2(2)
             write(103,*) 'Information on bound Hyperons:'
             write(103,'(a,1x,i8)') 'ID       = ',particleVector(in,HypIndex2(2))%number
             write(103,'(a,1x,f8.4)') 'born     = ',particleVector(in,HypIndex2(2))%bornTime
             write(103,'(a,1x,f8.4)') 'lastColl = ',particleVector(in,HypIndex2(2))%lastCollTime
             write(103,'(a,1x,i10)') 'collHis  = ',particleVector(in,HypIndex2(2))%collHis
             write(103,'(a,1x,3f12.4)') 'Pos.   = ',particleVector(in,HypIndex2(2))%position
             write(103,'(a,1x,3f12.4)') 'Vel.   = ', & 
                  & particleVector(in,HypIndex2(2))%momentum(1:3)/particleVector(in,HypIndex2(2))%momentum(0)
          endif
          write(103,*)
          write(103,'(a)') '##################################################################'
          write(103,*)

       enddo Loop_over_particles
    enddo Loop_over_ensembles
    close(103)

    i = in
    Particle_Loop2 : do j=1,size(ParticleVector,dim=2)

       if (ParticleVector(i,j)%ID .eq. 999) cycle !empty 

       if (ParticleVector(i,j)%ID==32) HyperCountAll = HyperCountAll + 1

    end do Particle_Loop2

    if (printFlag) then
       open(104,file='HyperInfoMul.dat',position='Append')

       write(104,'(a)') 'Total number of hyperons: '
       write(104,'(e12.4)') float(HyperCountAll)/float(SubEvents*NumEnsemples*SMM_Events)
       write(104,'(a)') 'Total number of hyperons bound in fragments: '
       write(104,'(e12.4)') float(HyperCount)/float(SubEvents*NumEnsemples*SMM_Events)
       write(104,'(a)') 'Number of hyperons bound in single-Lambda fragments: '
       write(104,'(e12.4)') float(HyperCountSingle)/float(SubEvents*NumEnsemples*SMM_Events)
       write(104,'(a)') 'Number of hyperons bound in double-Lambda fragments: '
       write(104,'(e12.4)') float(HyperCountDouble)/float(SubEvents*NumEnsemples*SMM_Events)

       close(104)
    endif

  !*************************************************************************
  end subroutine detailedHypAnalysis !**************************************
  !*************************************************************************

  !*************************************************************************
  !****s* FragmentationAnalysis/producedParticles
  ! FUNCTION
  ! Calculates various observables of produced mesons (pion,kaon,antikaon,...) 
  ! & produced baryons (hyperons, cascade-particles,...)
  ! 
  !*************************************************************************
  subroutine producedParticles(SubEvents,NumEnsemples, &
                               irun,ParticleVector,printflag,impactParameter)

    use typeDefinitions,  only : particle
    !-----------------------------------------------------------------------
    integer,                        intent(in) :: SubEvents,NumEnsemples
    integer,                        intent(in) :: irun
    type(particle), dimension(:,:), intent(in) :: ParticleVector
    logical,                        intent(in) :: printflag
    Real,                           intent(in) :: impactParameter
    !-----------------------------------------------------------------------
    integer :: i,j,ibine,ibiny,pID
    real    :: kinEn,k0f,k1f,k2f,k3f,yb,zz

    real, SAVE :: Norm=0.

    integer, parameter :: ipart=6 !!!!!!!!!!!! it was 5 !!!!!!!!!!!!!!!

    !kin. energy spectra:
    integer, parameter                    :: eVal = 250
    real,    parameter                    :: dekin = 0.01
    real, dimension(1:ipart,0:eVal), SAVE :: dNdEkin
    real, dimension(0:eVal),         SAVE :: Ebin

    !rapidity spectra:
    real,    parameter :: yMax  = -6.
    integer, parameter :: nrap  = 60
    real,    parameter :: dy    = 0.1

    real, dimension(1:ipart,-nrap:nrap), SAVE :: dndy
    real, dimension(-nrap:nrap),         SAVE :: ybin  

    real, dimension(1:3) :: xp
    real :: xsq

    logical, SAVE :: initParameters = .true.
    logical, SAVE :: initFlag4 = .true.

    integer :: numLL
    logical :: eventLL
    !-----------------------------------------------------------------------

    if (initParameters) then
       call GetAnalysisParameters
       initParameters = .false.
    endif

    if (initFlag4) then
       call SFInit
       initFlag4 = .false.
    endif

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! loop over BUU particle vector (produced mesons & baryons)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    i = irun

    numLL = 0
    eventLL = .false.

    !select events with 2 Lambdas:
    do j=1,size(ParticleVector,dim=2)
       if (ParticleVector(i,j)%ID .eq. 999) cycle !empty 
       pID = ParticleVector(i,j)%ID
       if (pID > 100) then
          eventLL = .false.
          exit
       endif
       if (pID==32) then
          numLL = numLL + 1
       endif
    end do
    if (numLL==2) eventLL=.true.

    Particle_Loop2 : do j=1,size(ParticleVector,dim=2)

       if (ParticleVector(i,j)%ID .eq. 999) cycle !empty 

       k0f= sqrt( & 
            & ParticleVector(i,j)%mass**2 + & 
            & Dot_Product(ParticleVector(i,j)%momentum(1:3),& 
            & ParticleVector(i,j)%momentum(1:3)))
       k3f= ParticleVector(i,j)%momentum(3)
       k1f= ParticleVector(i,j)%momentum(1)
       k2f= ParticleVector(i,j)%momentum(2)
       zz = float(ParticleVector(i,j)%Charge)
       pID = ParticleVector(i,j)%ID
       yb = 0.5*log( (k0f+k3f)/(k0f-k3f) ) / yproj
       kinEn = k0f - ParticleVector(i,j)%mass

       xp(1:3) = ParticleVector(i,j)%position(1:3)
       xsq = sqrt( dot_product(xp(:),xp(:)) )

       if (yb < 0.) then
          ibiny = int( yb/dy - 1.)
       else
          ibiny = int( yb/dy)
       endif
       if (abs(ibiny) > nrap) then
          write(*,*) 'upss... ',yb,ibiny
          STOP
       endif

       ibine = int( KinEn/dekin)
       ibine = min( eVal, max(0,ibine) )

       Select Case(pID)

       Case(101) !all pions
          call spectra(1)

       Case(32) !Lambdas
          call spectra(2)

          if (eventLL) then !events with 2 Lambdas !!!!!!!!!!!!!!
             call spectra(6)
          endif

       Case(53) !Xi's
          call spectra(3)

       Case(-32) !AntiLambdas
          call spectra(4)

       Case(-53) !AntiXi's
          call spectra(5)

       end Select
    end do Particle_Loop2

    !-----------------------------------------------------------------------
    FinalOutput : if (printFlag) then !-------------------------------------
    !-----------------------------------------------------------------------

       open(10,file='dNdY_producedParticles.dat')
       write(10,*) '# b,db:',ImpactParameter,impactParameter_bin
       write(10,*) '# y-bin; pions; Lambdas; Xi; Lambdabar; Xibar'
       do i=-nrap,nrap
          write(10,100) ybin(i),dndy(:,i)/dy
       end do
       close(10)

       open(10,file='dNdE_producedParticles.dat')
       write(10,*) '# b,db:',ImpactParameter,impactParameter_bin
       write(10,*) '# Ekin-bin; pions; Lambdas; Xi; Lambdabar; Xibar'
       do i=0,eVal
          write(10,100) Ebin(i),dndekin(:,i)/dekin
       end do
       close(10)

100    format(20e12.4)

    !-----------------------------------------------------------------------
    endif FinalOutput !-----------------------------------------------------
    !-----------------------------------------------------------------------

  contains 

    subroutine SFInit

      integer :: i
      real    :: ee, ystep

      Norm = 1./float(SubEvents*NumEnsemples)

      dNdEkin(:,:) = 0.0
      dndy(:,:)    = 0.0

      ee = 0.0
      do i=0,Eval
         Ebin(i)      = ee
         ee           = ee + dekin
      end do
      do i=0,Eval
         Ebin(i) = Ebin(i) + dekin/2.
      end do

      ystep = yMax
      do i=-nrap,nrap
         ybin(i) = ystep
         ystep = ystep + dy
      end do
      do i=-nrap,nrap
         ybin(i) = ybin(i) + dy/2.
      end do


    end subroutine SFInit


    subroutine Spectra(icount)

      integer, intent(in) :: icount

      dndy(icount,ibiny) = dndy(icount,ibiny) +  Norm
      dndekin(icount,ibine) = dndekin(icount,ibine) +  Norm

    end subroutine Spectra


  !*************************************************************************
  end subroutine ProducedParticles !****************************************
  !*************************************************************************

  !*************************************************************************
  !****s* FragmentationAnalysis/VariousDistributions
  ! FUNCTION
  ! Calculates various observables as function of the atomic (Z), 
  ! mass (A) and neutron (N) numbers and kinetic energy spectra of 
  ! various particles (protons, neutrons and mesons). 
  ! NOTES
  ! * Useful only for proton-induced reactions (Rejmund et al.)
  ! * (1) Do not normalize hier <E_kin>, since we performe the SMM analysis 
  !       for different impact parameters separately now. The normalization 
  !       is done after we have selected ALL the runs from b=b_min up to b=b_max !!!
  ! * Update(09.12.07): All observables contain information on 
  !   different mechanisms of fragmentation, saved in the last index:
  ! * [0,2]-->0:1: only fission-like, 2: only evaporation-like, 0: other proc.
  !   processes
  ! * Update (03.04.08) : Impact parameter weighted integration is performed 
  !   afterwards, see (1).
  ! 
  !*************************************************************************
  subroutine VariousDistributions(SubEvents,NumEnsemples,SMM_Events, & 
       & irun,FragmentVector,ParticleVector,XiTrigger,printflag,impactParameter)

    use typeDefinitions,  only : cluster,particle

    integer,                        intent(in) :: SubEvents,NumEnsemples
    integer,                        intent(in) :: irun,SMM_Events
    type(cluster),  dimension(:,:), intent(in) :: FragmentVector
    type(particle), dimension(:,:), intent(in) :: ParticleVector
    logical,                        intent(in) :: printflag, XiTrigger
    Real,                           intent(in) :: impactParameter

    integer :: i,j,l,ibin,izz,ibine,ibint,ibiny,model,pID
    real    :: kinEn,k0f,k1f,k2f,k3f,yb,zz,aa,nn,izz2

    integer :: HypFra
    character(2) :: HypType

    real, SAVE :: NormBUU=0., NormSMM=0., Norm2=0.

    real, dimension(1:9), SAVE :: XStot !total production yields of different heavy elements

    integer, parameter                   :: zval = 200 !common for Z-, N- & A-distributions
    real,    parameter                   :: dz = 1.
    real, dimension(zval),          SAVE :: zbin
    !distributions in Z (1. index: bin, 2. index: see variable "model"):
    real, dimension(zval,0:2),      SAVE :: Zdistr,Zdistr2, Zdistr_Xi
    !PANDA: distributions in Z & A including fragments with one (L) and two (LL) Lambdas:
    real, dimension(zval),          SAVE :: ZdistrL, ZdistrLL, AdistrL, AdistrLL
    real, dimension(zval),          SAVE :: ZdistrL_Xi, ZdistrLL_Xi, AdistrL_Xi, AdistrLL_Xi
    !For index explanation see Zdistr
    real, dimension(zval,0:2),      SAVE :: EkinZ,EkinNorm !distributions in Z
    real, dimension(zval,0:2),      SAVE :: Adistr, Adistr_Xi    !distributions in A
    real, dimension(1:13,zval,0:2), SAVE :: Adistr2              !distributions in A for different isotopes

    logical, SAVE :: initFlag1=.true.

    integer, parameter                           :: zmax=15,nval=140
    real,    parameter                           :: dn=1.
    real, dimension(1:nval),                SAVE :: nbin
    real, dimension(1:(zmax+2),1:nval,0:2), SAVE :: Ndistr2

    integer, parameter                        :: eVal = 250
    real,    parameter                        :: dekin = 0.005
    real, dimension(1:2,1:5,1:eVal,0:2), SAVE :: dNdEkin
    real, dimension(1:eVal),             SAVE :: ekinBin
    real, dimension(1:eVal,0:2),         SAVE :: Edistr,EdistrN
    real, dimension(1:eVal,0:20),        SAVE :: dndekinZ

    !Neutron,Proton & pion double differential spectra
    !separate origin of nucleon spectra (statistical-->SMM; pre-equilibrium-->BUU)
    integer, parameter                       :: tVal=16
!    real,    parameter                       :: ThetaBin=10.
    real,    parameter                       :: DEGRAD=180./pi
    real, dimension(1:eVal,1:tVal,1:2), SAVE :: dNdEkin2 
    real, dimension(1:eVal),            SAVE :: dNdEkin2_Pions
    real, dimension(1:eVal,1:tVal,1:2), SAVE :: dNdEkin2_SMM,dNdEkin2_BUU
    real, save :: NormSpectra=0.
    real :: theta

    !Variables for PANDA-Analysis:

    real,    parameter :: yMax  = -3.
    integer, parameter :: nrap  = 30
    real,    parameter :: dy    = 0.1

    real, dimension(-nrap:nrap,0:20), SAVE :: dndyZ
    real, dimension(-nrap:nrap),      SAVE :: ybin  

    !various distributions:
    integer, parameter                 :: mulMax=50
    real,    parameter                 :: dmul=1.

    !1.Index--> mult. bin, 2.Index--> 0:total, 1:evaporated (SMM), 2:emitted (BUU)
    real, dimension(0:mulMax,0:2),SAVE :: Panda_Ndistr !emitted neutrons
    real, dimension(0:mulMax,0:2),SAVE :: Panda_Zdistr !emitted protons
    real, dimension(0:mulMax,0:2),SAVE :: Panda_CP     !charged particles
    real, dimension(0:mulMax,0:2),SAVE :: Panda_Z1     !particles with Z=1

    real, dimension(0:mulMax),SAVE :: Panda_Z2     !particles with Z=2
    real, dimension(0:mulMax),SAVE :: Panda_IMF    !IMF's (5 < A < 25)

    real, dimension(0:2) :: freeNeutrons_RUN, freeProtons_RUN
    real, dimension(0:2) :: CP_RUN, Z1_RUN
    real ::Z2_RUN, IMF_RUN

    logical :: eventXi
    integer :: numXi

!    integer :: ibinRUN, ibinBUU
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (initParameters) then
       call GetAnalysisParameters
       initParameters = .false.
    endif

    if (initFlag1) then
       call ZdistInit
       call GetParametersSpectra
       initFlag1 = .false.
    endif

    !-----------------------------------------------------------------------
    ! Trigger: events with one Xi-Particle with E_beam < 0.08 GeV
    !-----------------------------------------------------------------------
    if (XiTrigger) then
       eventXi = .false.
       numXi = 0
       i = irun
       do j=1,size(ParticleVector,dim=2)
          if (ParticleVector(i,j)%ID .eq. 999) cycle !empty 
          k0f= sqrt( & 
               & ParticleVector(i,j)%mass**2 + & 
               & Dot_Product(ParticleVector(i,j)%momentum(1:3),& 
               & ParticleVector(i,j)%momentum(1:3)))
          pID = ParticleVector(i,j)%ID
          !       if (pID==53 .and. (k0f-ParticleVector(i,j)%mass) < 0.0801 ) then
          if (pID==53 ) then
             numXi = numXi + 1
             write(*,*) '############# ',i,pID,numXi
          endif
       end do

       if (numXi > 0) then
          eventXi=.true.
       endif
       if (numXi > 1) then
          write(*,*) 'more than one low energy Xi in event with number ',irun
          write(*,*) 'numXi = ',numXi
       endif
    else
       eventXi = .true.
       numXi   = 1
    endif


    !-----------------------------------------------------------------------
    ! init. of local variables needed for PANDA-Analysis:
    ! (pBar + X)
    !-----------------------------------------------------------------------
    freeNeutrons_RUN = 0.0 !emitted neutrons
    freeProtons_RUN  = 0.0 !emitted protons
    CP_RUN           = 0.0 !charged particles
    Z1_RUN           = 0.0 !particles with Z=1
    Z2_RUN           = 0.0 !particles with Z=2
    IMF_RUN          = 0.0 !particles with 5 < A < 25

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! LOOP OVER SMM EVENTS &  PRODUCED SMM PARTICLES 
    ! NORMALIZATION: NormSMM
    ! (NormSMM: Normalization for Loop over additional SMM-runs per BUU event!)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Ensemples_Loop1 : do i=1,size(FragmentVector,dim=1) 

       Particle_Loop1 : do j=1,size(FragmentVector,dim=2)

          if (FragmentVector(i,j)%ID==0) cycle

          k0f= sqrt( & 
               & FragmentVector(i,j)%mass**2 + & 
               & Dot_Product(FragmentVector(i,j)%momentum(1:3),& 
               & FragmentVector(i,j)%momentum(1:3)))
          k3f= FragmentVector(i,j)%momentum(3)
          k1f= FragmentVector(i,j)%momentum(1)
          k2f= FragmentVector(i,j)%momentum(2)
          zz = float(FragmentVector(i,j)%ChargeNumber) !charge number
          aa = float(FragmentVector(i,j)%MassNumber)   !mass number
          nn = aa - zz !neutron number

          HypFra  = FragmentVector(i,j)%HypNumber
          HypType = FragmentVector(i,j)%HypType

          yb = 0.5*log( (k0f+k3f)/(k0f-k3f) )/yproj

          !Separate info on fragmentation process:
          !1-> only fission / 2-> only evaporation / 0->others (fragmentation/break-up/deexcitation)
          model = FragmentVector(i,j)%Mechanism
          !check
          if (model < 0 .or. model > 2) then
             write(*,*) 'module fragmentationAnalysis / routine VariousDistributions:'
             write(*,*) 'Invalid fragmentation model: model =  ',model
             write(*,*) 'Termination of the programme'
             STOP
          endif

          !rapidity distributions (in Z, no free protons):
          if (yb < 0.) then
             ibiny = int( yb/dy - 1.)
          else
             ibiny = int( yb/dy)
          endif
          if (abs(ibiny) > nrap) then
             write(*,*) 'upss... ',yb,ibiny
             STOP
          endif
          if (FragmentVector(i,j)%ChargeNumber < 21) & 
               & dndyZ(ibiny,FragmentVector(i,j)%ChargeNumber) =  & 
               & dndyZ(ibiny,FragmentVector(i,j)%ChargeNumber) + NormSMM


          !inclusive production cross sections of different elements
          !for comparison with data from paper: Yu.E.Titarenko et al.
          !considered elements: Co,Mn,Sc,K,Be
          if (zz==27.) then
             if (aa==57.) XStot(1) = XStot(1) + NormSMM
             if (aa==56.) XStot(2) = XStot(2) + NormSMM
             if (aa==55.) XStot(3) = XStot(3) + NormSMM
          endif
          if (zz==25. .and. aa==54.) then
             XStot(4) = XStot(4) + NormSMM
          endif
          if (zz==21.) then
             if (aa==48.) XStot(5) = XStot(5) + NormSMM
             if (aa==47.) XStot(6) = XStot(6) + NormSMM
             if (aa==44.) XStot(7) = XStot(7) + NormSMM
          endif
          if (zz==19. .and. aa==42.) then
             XStot(8) = XStot(8) + NormSMM
          endif
          if (zz==4. .and. aa==7.) then
             XStot(9) = XStot(9) + NormSMM
          endif


          !For test: kinEn1 coincides now with KinEn. 
          !kinEn1 = FragmentVector(i,j)%momentum(0) - FragmentVector(i,j)%mass
          KinEn = k0f - FragmentVector(i,j)%mass

          !Distributions in A
          ibin = int( aa/dz)
          ibin = min( zval, max(1,ibin) )
          Adistr(ibin,model) = Adistr(ibin,model) + NormSMM

          !(PANDA) Mass distribution for single-Lambda hypernuclei:
          if(HypFra==1 .and. trim(HypType)=='L') AdistrL(ibin) = AdistrL(ibin) + NormSMM

          !(PANDA) Mass distribution for double-Lambda hypernuclei:
          if(HypFra==2 .and. trim(HypType)=='LL') AdistrLL(ibin) = AdistrLL(ibin) + NormSMM

          !(PANDA) trigger events with one Xi:
          if (eventXi) then
             Adistr_Xi(ibin,model) = Adistr_Xi(ibin,model) + NormSMM
             if(HypFra==1 .and. trim(HypType)=='L') AdistrL_Xi(ibin) = AdistrL_Xi(ibin) + NormSMM
             if(HypFra==2 .and. trim(HypType)=='LL') AdistrLL_Xi(ibin) = AdistrLL_Xi(ibin) + NormSMM
          endif

          if (zz==20.) Adistr2(1,ibin,model) = Adistr2(1,ibin,model)   + NormSMM
          if (zz==25.) Adistr2(2,ibin,model) = Adistr2(2,ibin,model)   + NormSMM
          if (zz==30.) Adistr2(3,ibin,model) = Adistr2(3,ibin,model)   + NormSMM
          if (zz==35.) Adistr2(4,ibin,model) = Adistr2(4,ibin,model)   + NormSMM
          if (zz==40.) Adistr2(5,ibin,model) = Adistr2(5,ibin,model)   + NormSMM
          if (zz==45.) Adistr2(6,ibin,model) = Adistr2(6,ibin,model)   + NormSMM
          if (zz==50.) Adistr2(7,ibin,model) = Adistr2(7,ibin,model)   + NormSMM
          if (zz==55.) Adistr2(8,ibin,model) = Adistr2(8,ibin,model)   + NormSMM
          if (zz==60.) Adistr2(9,ibin,model) = Adistr2(9,ibin,model)   + NormSMM
          if (zz==65.) Adistr2(10,ibin,model) = Adistr2(10,ibin,model) + NormSMM
          if (zz==70.) Adistr2(11,ibin,model) = Adistr2(11,ibin,model) + NormSMM
          if (zz==75.) Adistr2(12,ibin,model) = Adistr2(12,ibin,model) + NormSMM
          if (zz==80.) Adistr2(13,ibin,model) = Adistr2(13,ibin,model) + NormSMM

          !Distributions in N

          ibin = int( nn/dn)
          ibin = min( nval, max(1,ibin) )

          if (nn > 1) then !no elements below this neutron number
             izz2 = 5.
             do izz=1,zmax !loop over fixed charge number (Z=5-75)
                if (zz==izz2) Ndistr2(izz,ibin,model) = Ndistr2(izz,ibin,model)+NormSMM
                izz2 = izz2 + 5.
             end do
             if (zz==80.) Ndistr2(zmax+1,ibin,model) = Ndistr2(zmax+1,ibin,model) + NormSMM
             if (zz==97.) Ndistr2(zmax+2,ibin,model) = Ndistr2(zmax+2,ibin,model) + NormSMM
          endif

          !Kin. Energy spectra of free protons & neutrons
          !here low-energy part (from statistical decay of residual nucleus)
          ibine = int( KinEn/dekin)
          ibine = min( eVal, max(1,ibine) )
          if (ibine < 1 .or. ibine > eVal) then
             write(*,*) '1-something is wrong in analysis-module (VariousDistributions)!!!', & 
                  & kinEn,ibine
             STOP
          endif

          !kin. energy distributions (in Z, no free protons):
          if (FragmentVector(i,j)%ChargeNumber < 21) & 
               & dndekinZ(ibine,FragmentVector(i,j)%ChargeNumber) =  & 
               & dndekinZ(ibine,FragmentVector(i,j)%ChargeNumber) + NormSMM


          if (aa==1) then
             theta=atan2(sqrt(k1f**2+k2f**2),k3f)*degrad
             !          ibint = int( theta/thetabin)
             !          ibint = min( tVal, max(1,ibint) )
             call PolarWinkelBin(theta,ibint)
             if (ibint /= 0) then
                if (zz==1.) then
                   dNdEkin2(ibine,ibint,1)     = dNdEkin2(ibine,ibint,1) + NormSMM
                   dNdEkin2_SMM(ibine,ibint,1) = dNdEkin2_SMM(ibine,ibint,1) + NormSMM
                endif
                if (zz==0.) then
                   dNdEkin2(ibine,ibint,2)     = dNdEkin2(ibine,ibint,2) + NormSMM
                   dNdEkin2_SMM(ibine,ibint,2) = dNdEkin2_SMM(ibine,ibint,2) + NormSMM
                endif
             end if
          end if

          !Distributions in Z

          !PANDA-Analysis (pBar+X):
          if (aa==1. .and. zz==0.) then
             freeNeutrons_RUN(0) = freeNeutrons_RUN(0) + 1. 
             freeNeutrons_RUN(1) = freeNeutrons_RUN(1) + 1. 
          endif
          if (aa==1. .and. zz==1.) then
             freeProtons_RUN(0) = freeProtons_RUN(0) + 1. 
             freeProtons_RUN(1) = freeProtons_RUN(1) + 1. 
          endif
          if (zz.ge.1) then
             CP_RUN(0) = CP_RUN(0) + 1.
             CP_RUN(1) = CP_RUN(1) + 1.
          endif
          if (zz.eq.1) then
             Z1_RUN(0) = Z1_RUN(0) + 1.
             Z1_RUN(1) = Z1_RUN(1) + 1.
          endif
          if (zz.eq.2) Z2_RUN = Z2_RUN + 1.
          if (aa.gt.5 .and. aa.lt.25) IMF_RUN = IMF_RUN + 1.

          if (zz.eq.0.) cycle

          ibin = int( zz/dz)
          ibin = min( zval, max(1,ibin) )

          !Charge particle nultiplicity distribution:
          Zdistr(ibin,model) = Zdistr(ibin,model) + NormSMM

          !(PANDA) Charge particle nultiplicity distribution for single-Lambda hypernuclei:
          if(HypFra==1 .and. trim(HypType)=='L') ZdistrL(ibin) = ZdistrL(ibin) + NormSMM

          !(PANDA) Charge particle nultiplicity distribution for double-Lambda hypernuclei:
          if(HypFra==2 .and. trim(HypType)=='LL') ZdistrLL(ibin) = ZdistrLL(ibin) + NormSMM

          !(PANDA) trigger events with one Xi-particle:
          if (eventXi) then
             Zdistr_Xi(ibin,model) = Zdistr_Xi(ibin,model) + NormSMM
             if(HypFra==1 .and. trim(HypType)=='L') ZdistrL_Xi(ibin) = ZdistrL_Xi(ibin) + NormSMM
             if(HypFra==2 .and. trim(HypType)=='LL') ZdistrLL_Xi(ibin) = ZdistrLL_Xi(ibin) + NormSMM
          endif

          EkinZ(ibin,model) = EkinZ(ibin,model) + KinEn    !Sum of E_kin at each dZ-bin
          EkinNorm(ibin,model) = EkinNorm(ibin,model) + 1. !Normalization (SEE NOTES!!!)
 
          !Kin. Energy spectra of different isotopes
          Edistr(ibine,model) = Edistr(ibine,model) + KinEn !total
          EdistrN(ibine,model) = EdistrN(ibine,model) + 1.

          !dNdEkin for different Z-regions
          if (zz > 5. .and. zz < 15.) then
             dNdEkin(1,1,ibine,model) = dNdEkin(1,1,ibine,model) + KinEn 
             dNdEkin(2,1,ibine,model) = dNdEkin(2,1,ibine,model) + 1.    
          endif
          if (zz>15. .and. zz<25.) then
             dNdEkin(1,2,ibine,model) = dNdEkin(1,2,ibine,model) + KinEn 
             dNdEkin(2,2,ibine,model) = dNdEkin(2,2,ibine,model) + 1.    
          endif
          if (zz>25. .and. zz<35.) then
             dNdEkin(1,3,ibine,model) = dNdEkin(1,3,ibine,model) + KinEn 
             dNdEkin(2,3,ibine,model) = dNdEkin(2,3,ibine,model) + 1.    
          endif
          if (zz>35. .and. zz<45.) then
             dNdEkin(1,4,ibine,model) = dNdEkin(1,4,ibine,model) + KinEn 
             dNdEkin(2,4,ibine,model) = dNdEkin(2,4,ibine,model) + 1.    
          endif
          if (zz>45. .and. zz<55.) then
             dNdEkin(1,5,ibine,model) = dNdEkin(1,5,ibine,model) + KinEn 
             dNdEkin(2,5,ibine,model) = dNdEkin(2,5,ibine,model) + 1.    
          endif

          if (abs(yb) .le. 0.5) Zdistr2(ibin,model) = Zdistr2(ibin,model) + NormSMM

       end do Particle_Loop1

    end do Ensemples_Loop1

    !PANDA:
    ! normalize here, if extra SMM events are requested!
    freeNeutrons_RUN(0:2) = freeNeutrons_RUN(0:2)/float(SMM_Events) !emitted neutrons
    freeProtons_RUN(0:2)  = freeProtons_RUN(0:2)/float(SMM_Events) !emitted neutrons
    CP_RUN(0:2)           = CP_RUN(0:2)/float(SMM_Events) !charged particles
    Z1_RUN(0:2)           = Z1_RUN(0:2)/float(SMM_Events) !particles with Z=1
    Z2_RUN           = Z2_RUN/float(SMM_Events) !particles with Z=2
    IMF_RUN          = IMF_RUN/float(SMM_Events) !particles with 5 < A < 25


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! LOOP OVER BUU PARTICLES (only emitted ones)
    ! NORMALIZATION: NormBUU
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    i = irun
    Particle_Loop2 : do j=1,size(ParticleVector,dim=2)

       if (ParticleVector(i,j)%ID .eq. 999) cycle !empty 
!       if (ParticleVector(i,j)%ID .ne. 1) cycle !only nucleons
!!!       if (ParticleVector(i,j)%SD .ne. 999) cycle !only emitted nucleons

       k0f= sqrt( & 
            & ParticleVector(i,j)%mass**2 + & 
            & Dot_Product(ParticleVector(i,j)%momentum(1:3),& 
            & ParticleVector(i,j)%momentum(1:3)))
       k3f= ParticleVector(i,j)%momentum(3)
       k1f= ParticleVector(i,j)%momentum(1)
       k2f= ParticleVector(i,j)%momentum(2)
       zz = float(ParticleVector(i,j)%Charge)
       pID = ParticleVector(i,j)%ID

!       kinEn1 = ParticleVector(i,j)%momentum(0) - ParticleVector(i,j)%Mass
       KinEn = k0f - ParticleVector(i,j)%Mass

       !Kin. Energy spectra of emitted protons & neutrons
       !here high-energy part (from pre-equilibrium emission)
       ibine = int( KinEn/dekin)
       ibine = min( eVal, max(1,ibine) )
       if (ibine < 1 .or. ibine > eVal) then
          write(*,*) '1-something is wrong in analysis-module (VariousDistributions)!!!', & 
               & kinEn,ibine
          STOP
       endif

       theta=atan2(sqrt(k1f**2+k2f**2),k3f)*degrad
       !          ibint = int( theta/thetabin)
       !          ibint = min( tVal, max(1,ibint) )
       call PolarWinkelBin(theta,ibint)

       Select Case(pID)
       Case(1) !nucleons
          if (zz==1. .and. ibint /= 0) then
             dNdEkin2(ibine,ibint,1)     = dNdEkin2(ibine,ibint,1) + NormBUU
             dNdEkin2_BUU(ibine,ibint,1) = dNdEkin2_BUU(ibine,ibint,1) + NormBUU
          endif
          if (zz==0. .and. ibint /= 0) then
             dNdEkin2(ibine,ibint,2)     = dNdEkin2(ibine,ibint,2) + NormBUU
             dNdEkin2_BUU(ibine,ibint,2) = dNdEkin2_BUU(ibine,ibint,2) + NormBUU
          endif

       Case(101) !pions
          dNdEkin2_pions(ibine) = dNdEkin2_pions(ibine) + KinEn * NormBUU
       end Select

       !--------------------------------------------------------------------
       if (pID > 1) cycle !below only free nucleons (protons, neutrons)
       !--------------------------------------------------------------------

       ibin = int( zz/dz)
       ibin = min( zval, max(1,ibin) )

       !just to be sure...:
       if (ibin > 1) then
          write(*,*) 'do not count resonances!',pID,ibin,zz
       endif

       !PANDA-Analysis (pBar+X):
       if (zz.ge.1) then
          CP_RUN(0) = CP_RUN(0) + 1.
          CP_RUN(2) = CP_RUN(2) + 1.
       endif
       if (zz.eq.1) then
          Z1_RUN(0) = Z1_RUN(0) + 1.
          Z1_RUN(2) = Z1_RUN(2) + 1.
          freeProtons_RUN(0) = freeProtons_RUN(0) + 1. !total
          freeProtons_RUN(1) = freeProtons_RUN(1) + 1. !only pre-equilibrium
       endif

       !protons:
       if (zz .ne. 0.) then

          yb = 0.5*log( (k0f+k3f)/(k0f-k3f) )/yproj

          Zdistr(ibin,0) = Zdistr(ibin,0) + NormBUU

          if (abs(yb) .le. 0.5) Zdistr2(ibin,0) = Zdistr2(ibin,0) + NormBUU

       !neutrons: 
       else 
          freeNeutrons_RUN(0) = freeNeutrons_RUN(0) +1. !total
          freeNeutrons_RUN(2) = freeNeutrons_RUN(2) +1. !only pre-equilibrium
       endif

    end do Particle_Loop2

    !-----------------------------------------------------------------------
    !PANDA-Analysis (pBar+X):
    !-----------------------------------------------------------------------

    do i=0,2

       !Neutron multiplicity distributions:
       ibin = int( freeNeutrons_RUN(i)/dmul)
       ibin = min( mulMax, max(0,ibin) )
       Panda_Ndistr(ibin,i) = Panda_Ndistr(ibin,i) + 1.

       !Proton multiplicity distributions:
       ibin = int( freeProtons_RUN(i)/dmul)
       ibin = min( mulMax, max(0,ibin) )
       Panda_Zdistr(ibin,i) = Panda_Zdistr(ibin,i) + 1.

       !Charged-particle  multiplicity distributions:
       ibin = int( CP_RUN(i)/dmul)
       ibin = min( mulMax, max(0,ibin) )
       Panda_CP(ibin,i) = Panda_CP(ibin,i) + 1.
       ibin = int( Z1_RUN(i)/dmul)
       ibin = min( mulMax, max(0,ibin) )
       Panda_Z1(ibin,i) = Panda_Z1(ibin,i) + 1.

    end do

    !Z=2-particle  multiplicity distributions:
    ibin = int( Z2_RUN/dmul)
    ibin = min( mulMax, max(0,ibin) )
    Panda_Z2(ibin) = Panda_Z2(ibin) + 1.
    !IMF-particle  multiplicity distributions:
    ibin = int( IMF_RUN/dmul)
    ibin = min( mulMax, max(0,ibin) )
    Panda_IMF(ibin) = Panda_IMF(ibin) + 1.

    !-----------------------------------------------------------------------
    FinalOutput : if (printFlag) then !-------------------------------------
    !-----------------------------------------------------------------------

       open(40,file='LadungsVerteilung_allProcesses.dat')
       write(40,*) '# b,db:',ImpactParameter,impactParameter_bin
       do i=1,zval
          if (i>100) exit !for charge distributions
          write(40,100) zbin(i),Sum(Zdistr(i,0:2))/dz,& 
               & Sum(Zdistr2(i,0:2))/dz, & 
               & Sum(Adistr(i,0:2))/dz, & 
               & Sum(EkinZ(i,0:2)),Sum(EkinNorm(i,0:2))
       end do
       open(41,file='LadungsVerteilung_otherProcesses.dat')
       write(41,*) '# b,db:',ImpactParameter,impactParameter_bin
       do i=1,zval
          if (i>100) exit !for charge distributions
          write(41,100) zbin(i),Zdistr(i,0)/dz,& 
               & Zdistr2(i,0)/dz, & 
               & Adistr(i,0)/dz,EkinZ(i,0),EkinNorm(i,0)
       end do
       open(42,file='LadungsVerteilung_OnlyFission.dat')
       write(42,*) '# b,db:',ImpactParameter,impactParameter_bin
       do i=1,zval
          if (i>100) exit !for charge distributions
          write(42,100) zbin(i),Zdistr(i,1)/dz,& 
               & Zdistr2(i,1)/dz, & 
               & Adistr(i,1)/dz,EkinZ(i,1),EkinNorm(i,1)
       end do
       open(43,file='LadungsVerteilung_OnlyEvaporation.dat')
       write(43,*) '# b,db:',ImpactParameter,impactParameter_bin
       do i=1,zval
          if (i>100) exit !for charge distributions
          write(43,100) zbin(i),Zdistr(i,2)/dz,& 
               & Zdistr2(i,2)/dz, & 
               & Adistr(i,2)/dz,EkinZ(i,2),EkinNorm(i,2)
       end do
100    format(6e12.4)

       close(40)
       close(41)
       close(42)
       close(43)

       open(50,file='MassVerteilung_allProcesses.dat')
       write(50,*) '# b,db:',ImpactParameter,impactParameter_bin
       do i=1,zval
          write(50,101) zbin(i),& 
               & Sum(Adistr(i,0:2))/dz, & 
               & (Sum(Adistr2(l,i,0:2))/dz,l=1,13)
       end do
       open(51,file='MassVerteilung_otherProcesses.dat')
       write(51,*) '# b,db:',ImpactParameter,impactParameter_bin
       do i=1,zval
          write(51,101) zbin(i),& 
               & Adistr(i,0)/dz, & 
               & (Adistr2(l,i,0)/dz,l=1,13)
       end do
       open(52,file='MassVerteilung_OnlyFission.dat')
       write(52,*) '# b,db:',ImpactParameter,impactParameter_bin
       do i=1,zval
          write(52,101) zbin(i),& 
               & Adistr(i,1)/dz, & 
               & (Adistr2(l,i,1)/dz,l=1,13)
       end do
       open(53,file='MassVerteilung_OnlyEvaporation.dat')
       write(53,*) '# b,db:',ImpactParameter,impactParameter_bin
       do i=1,zval
          write(53,101) zbin(i),& 
               & Adistr(i,2)/dz, & 
               & (Adistr2(l,i,2)/dz,l=1,13)
       end do
101    format(15e12.4)

       close(50)
       close(51)
       close(52)
       close(53)

       open(60,file='NeutronVerteilung_allProcesses.dat')
       write(60,*) '# b,db:',ImpactParameter,impactParameter_bin
       do i=1,nval
          write(60,102) nbin(i),& 
               & (Sum(Ndistr2(l,i,0:2))/dn,l=1,zmax+2)
       end do
       open(61,file='NeutronVerteilung_otherProcesses.dat')
       write(61,*) '# b,db:',ImpactParameter,impactParameter_bin
       do i=1,nval
          write(61,102) nbin(i),& 
               & (Ndistr2(l,i,0)/dn,l=1,zmax+2)
       end do
       open(62,file='NeutronVerteilung_OnlyFission.dat')
       write(62,*) '# b,db:',ImpactParameter,impactParameter_bin
       do i=1,nval
          write(62,102) nbin(i),& 
               & (Ndistr2(l,i,1)/dn,l=1,zmax+2)
       end do
       open(63,file='NeutronVerteilung_OnlyEvaporation.dat')
       write(63,*) '# b,db:',ImpactParameter,impactParameter_bin
       do i=1,nval
          write(63,102) nbin(i),& 
               & (Ndistr2(l,i,2)/dn,l=1,zmax+2)
       end do

       !PANDA-Analysis (pBar+X):
       open(64,file='MultiplicityDistributions_PANDA.dat')
       write(64,*) '# b,db:',ImpactParameter,impactParameter_bin
       do i=0,MulMax
          write(64,102) float(i)*dmul, & 
               & (Panda_Zdistr(i,l)*Norm2/dmul, l=0,2) , & 
               & (Panda_Ndistr(i,l)*Norm2/dmul, l=0,2) , & 
               & (Panda_CP(i,l)*Norm2/dmul, l=0,2) ,  & 
               & (Panda_Z1(i,l)*Norm2/dmul, l=0,2) ,  & 
               & Panda_Z2(i)*Norm2/dmul, & 
               & Panda_IMF(i)*Norm2/dmul
       end do
102    format(25e12.4)

       open(65,file='LadungsVerteilung_PANDA.dat')
       write(65,*) '# b,db:',ImpactParameter,impactParameter_bin
       do i=1,zval
          if (i>100) exit !for charge distributions
          write(65,'(7e12.4)') zbin(i), & 
               & Sum(Zdistr(i,0:2))/dz, ZdistrL(i)/dz, ZdistrLL(i)/dz, & 
               & Sum(Zdistr_Xi(i,0:2))/dz, ZdistrL_Xi(i)/dz, ZdistrLL_Xi(i)/dz 
       end do

       open(66,file='MassVerteilung_PANDA.dat')
       write(66,*) '# b,db:',ImpactParameter,impactParameter_bin
       do i=1,zval
          write(66,'(7e12.4)') zbin(i), & 
               & Sum(Adistr(i,0:2))/dz , AdistrL(i)/dz, AdistrLL(i)/dz, & 
               & Sum(Adistr_Xi(i,0:2))/dz , AdistrL_Xi(i)/dz, AdistrLL_Xi(i)/dz
       end do


       close(60)
       close(61)
       close(62)
       close(63)
       close(64)
       close(65)
       close(66)

       open(70,file='EnergieVerteilung_allProcesses.dat')
       write(70,*) '# b,db:',ImpactParameter,impactParameter_bin
       do i=1,eVal
          write(70,103) ekinBin(i),Sum(Edistr(i,0:2)),Sum(EdistrN(i,0:2)), & 
               & (Sum(dNdEkin(1,j,i,0:2)),Sum(dNdEkin(2,j,i,0:2)),j=1,5)
       end do
       open(71,file='EnergieVerteilung_otherProcesses.dat')
       write(71,*) '# b,db:',ImpactParameter,impactParameter_bin
       do i=1,eVal
          write(71,103) ekinBin(i),Edistr(i,0),EdistrN(i,0), & 
               & (dNdEkin(1,j,i,0),dNdEkin(2,j,i,0),j=1,5)
       end do
       open(72,file='EnergieVerteilung_OnlyFission.dat')
       write(72,*) '# b,db:',ImpactParameter,impactParameter_bin
       do i=1,eVal
          write(72,103) ekinBin(i),Edistr(i,1),EdistrN(i,1), & 
               & (dNdEkin(1,j,i,1),dNdEkin(2,j,i,1),j=1,5)
       end do
       open(73,file='EnergieVerteilung_OnlyEvaporation.dat')
       write(73,*) '# b,db:',ImpactParameter,impactParameter_bin
       do i=1,eVal
          write(73,103) ekinBin(i),Edistr(i,2),EdistrN(i,2), & 
               & (dNdEkin(1,j,i,2),dNdEkin(2,j,i,2),j=1,5)
       end do
103    format(32e12.4)

       close(70)
       close(71)
       close(72)
       close(73)

       open(70,file='KineticEnergySpectraZ.dat')
       write(70,*) '# b,db:',ImpactParameter,impactParameter_bin
       do i=1,eVal
          write(70,103) ekinBin(i),dndekinZ(i,:)/dekin
       end do
       close(70)

       open(70,file='dNdY_Z.dat')
       write(70,*) '# b,db:',ImpactParameter,impactParameter_bin
       do i=-nrap,nrap
          write(70,103) ybin(i),dndyZ(i,:)/dy
       end do
       close(70)

       !Double differential particle spectra
       open(80,file='KineticEnergySpectra_Protons.dat')
       open(81,file='KineticEnergySpectra_Neutrons.dat')
       open(82,file='KineticEnergySpectra_Pions.dat')
       open(83,file='KineticEnergySpectra_ProtonsSMM.dat')
       open(84,file='KineticEnergySpectra_ProtonsBUU.dat')
       open(85,file='KineticEnergySpectra_NeutronsSMM.dat')
       open(86,file='KineticEnergySpectra_NeutronsBUU.dat')
       do i=1,7
          write(79+i,*) '# b,db:',ImpactParameter,impactParameter_bin
       end do
       do i=1,eVal
          write(80,103) ekinBin(i), & 
               & (dNdEkin2(i,j,1)*NormSpectra,j=1,tVal)
          write(81,103) ekinBin(i), & 
               & (dNdEkin2(i,j,2)*NormSpectra,j=1,tVal)
          write(82,103) ekinBin(i), & 
               & dNdEkin2_Pions(i)*NormSpectra
          write(83,103) ekinBin(i), & 
               & (dNdEkin2_SMM(i,j,1)*NormSpectra,j=1,tVal)
          write(84,103) ekinBin(i), & 
               & (dNdEkin2_BUU(i,j,1)*NormSpectra,j=1,tVal)
          write(85,103) ekinBin(i), & 
               & (dNdEkin2_SMM(i,j,2)*NormSpectra,j=1,tVal)
          write(86,103) ekinBin(i), & 
               & (dNdEkin2_BUU(i,j,2)*NormSpectra,j=1,tVal)
       end do

       do i=1,7
          close(79+i)
       end do

       open(90,file='TotalYields.dat')
       write(90,*) '# b,db:',ImpactParameter,impactParameter_bin
       write(90,103) XStot(1:9)
       
       

    !-----------------------------------------------------------------------
    endif FinalOutput !-----------------------------------------------------
    !-----------------------------------------------------------------------

  contains 

    subroutine ZdistInit

      integer :: i
      real :: wert,thetaBin_sr,ystep

      NormBUU = 1./float(SubEvents*NumEnsemples)
      NormSMM = 1./float(SubEvents*NumEnsemples*SMM_Events)
      Norm2   = 1./float(SubEvents*NumEnsemples)

      !initialize XStot (total yield of several elements)
      XStot(1:9) = 0.0

      !initialize distributions in Z- and A-number
      Zdistr(1:zval,0:2)       = 0.
      ZdistrL(1:zval)          = 0.
      ZdistrLL(1:zval)         = 0.

      Zdistr_Xi(1:zval,0:2)       = 0.
      ZdistrL_Xi(1:zval)          = 0.
      ZdistrLL_Xi(1:zval)         = 0.

      Zdistr2(1:zval,0:2)      = 0.
      Adistr(1:zval,0:2)       = 0.
      AdistrL(1:zval)          = 0.
      AdistrLL(1:zval)         = 0.

      Adistr_Xi(1:zval,0:2)       = 0.
      AdistrL_Xi(1:zval)          = 0.
      AdistrLL_Xi(1:zval)         = 0.

      EkinZ(1:zval,0:2)        = 0.
      EkinNorm(1:zval,0:2)     = 0.
      Adistr2(1:13,1:zval,0:2) = 0.
      do i=1,zval
         zbin(i) = float(i)*dz
      end do
      
      !initialize distributions in neutron number
      Ndistr2(1:(zmax+2),1:nval,0:2) = 0.
      do i=1,nval
         nbin(i) = float(i)*dn
      end do
      Panda_Ndistr(0:mulMax,0:2) = 0.0
      Panda_Zdistr(0:mulMax,0:2) = 0.0
      Panda_CP(0:mulMax,0:2)     = 0.0
      Panda_Z1(0:mulMax,0:2)     = 0.0
      Panda_Z2(0:mulMax)     = 0.0
      Panda_IMF(0:mulMax)    = 0.0

      !initialize distributions in kinetic energy
      dNdEkin(1:2,1:5,1:eVal,0:2)  = 0.0
      dndekinZ(1:eVal,0:20)        = 0.0
      Edistr(1:eVal,0:2)           = 0.0
      EdistrN(1:eVal,0:2)          = 0.0
      dNdEkin2(1:eval,1:tval,1:2)  = 0.0
      dNdEkin2_pions(1:eval)       = 0.0
      dNdEkin2_SMM(1:eval,1:tval,1:2) = 0.0
      dNdEkin2_BUU(1:eval,1:tval,1:2) = 0.0
      wert = 0.0
      do i=1,eVal
         ekinBin(i) = wert
         wert = wert + dekin
      end do
      ekinBin(:) = ekinBin(:) + dekin/2.

      !normalization of double differential spectra:
!      thetaBin_sr = 2.*pi*(1.-cos(pi*thetabin/180.)) 
      thetaBin_sr = 2.*pi*(1.-cos(pi*deltaTheta*2./180.)) 
      NormSpectra = 1./thetaBin_sr/dekin

      dndyZ(:,:) = 0.0
      ystep = yMax
      do i=-nrap,nrap
         ybin(i) = ystep
         ystep = ystep + dy
      end do
      do i=-nrap,nrap
         ybin(i) = ybin(i) + dy/2.
      end do


    end subroutine ZdistInit



    subroutine PolarWinkelBin(theta,ibint)

      integer :: i
      real,    intent(in)  :: theta
      integer, intent(out) :: ibint

      ibint = 0

      do i=1,ThetaMax
         if (Angles(i)==1000.) cycle !not readed from jobCard
         if (i==1) then
            if (theta > Angles(i) .and. & 
                 & theta < (Angles(i)+2.*deltaTheta) ) then
               ibint = i
               exit
            end if
         else
            if (theta > (Angles(i)-deltaTheta) .and. & 
                 & theta < (Angles(i)+deltaTheta) ) then
               ibint = i
               exit
            endif
         end if
      end do

    end subroutine PolarWinkelBin



  !*************************************************************************
  end subroutine VariousDistributions !*************************************
  !*************************************************************************



  !*************************************************************************
  !****s* FragmentationAnalysis/ExDistributions
  ! FUNCTION
  ! Calculates various observables as function of excitation energy.
  ! 
  !*************************************************************************
  subroutine ExDistributions(SubEvents,NumEnsemples, SMM_Events,  & 
       & irun,FragmentVector,ParticleVector,TheSource,printflag,impactParameter)

    use typeDefinitions

    integer,                        intent(in) :: SubEvents,NumEnsemples,SMM_Events
    integer,                        intent(in) :: irun
    type(cluster),  dimension(:,:), intent(in) :: FragmentVector
    type(particle), dimension(:,:), intent(in) :: ParticleVector
    type(quelle),   dimension(:,:), intent(in) :: TheSource
    logical,                        intent(in) :: printflag
    Real,                           intent(in) :: impactParameter

    integer :: i,j,ibin,pID
    real    :: zz,aa,nn,wert_CP,wert_IMF,Ex

    integer, parameter             :: ExMax = 150 
    real,    parameter             :: Ex_step = 10.
    real, dimension(1:ExMax), SAVE :: Ex_bin
    real, dimension(1:ExMax), SAVE :: dNdEx,dNdEx_CP,dNdEx_IMF !,Norm_CP,Norm_IMF

    real, save :: Norm=0.0
    logical, SAVE :: initFlag1=.true.

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (initFlag1) then
       call ExInit
       initFlag1 = .false.
    endif

    if ( .not. TheSource(iRun,1)%status ) then
       write(*,*) 'ExDistributions: '
       write(*,*) 'Actual event does not contain any source.'
       return
    endif
    
    !total excitation energy of the residual source (MeV):
    Ex = ( TheSource(iRun,1)%ExEnergy * TheSource(iRun,1)%Size )*1000.

    ibin = int(Ex/Ex_step)
    ibin = min(ExMax, max(1,ibin) )

    !Excitation energy distribution:
    dNdEx(ibin) = dNdEx(ibin) + 1.

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! LOOP OVER SMM EVENTS &  PRODUCED SMM PARTICLES 
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Ensemples_Loop1 : do i=1,size(FragmentVector,dim=1) 

       Particle_Loop1 : do j=1,size(FragmentVector,dim=2)

          if (FragmentVector(i,j)%ID==0) cycle

          aa = float(FragmentVector(i,j)%MassNumber)   !mass number
          zz = float(FragmentVector(i,j)%ChargeNumber)   !charge number
          nn = aa - zz !neutron number

          !Charged Particles (CP):
          if (zz >= 1.) then
             dNdEx_CP(ibin) = dNdEx_CP(ibin) + 1./float(SMM_Events) !zz
!             Norm_CP(ibin)  = Norm_CP(ibin) + 1.
          endif

          !Intermediate Mass Fragments (IMF):
          if ( (zz >= 5.).and.(zz <= 25.) ) then
             dNdEx_IMF(ibin) = dNdEx_IMF(ibin) + 1./float(SMM_Events) !zz
!             Norm_IMF(ibin)  = Norm_IMF(ibin) + 1.
          endif

       end do Particle_Loop1

    end do Ensemples_Loop1

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! LOOP OVER BUU PARTICLES (only emitted ones)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    i = irun
    Particle_Loop2 : do j=1,size(ParticleVector,dim=2)

       if (ParticleVector(i,j)%ID .eq. 999) cycle !empty 

       zz = float(ParticleVector(i,j)%Charge)
       pID = ParticleVector(i,j)%ID

       !--------------------------------------------------------------------
       if (pID > 1) cycle ! only free nucleons
       !--------------------------------------------------------------------

       if (zz == 1. ) then
             dNdEx_CP(ibin) = dNdEx_CP(ibin) + 1. !zz
!             Norm_CP(ibin)  = Norm_CP(ibin) + 1.
       endif

    end do Particle_Loop2

    !-----------------------------------------------------------------------
    FinalOutput : if (printFlag) then !-------------------------------------
    !-----------------------------------------------------------------------

       open(40,file='ExDistributions.dat')
       write(40,*) '# b,db:',ImpactParameter,impactParameter_bin
       do i=1,ExMax
!!$          if (Norm_CP(i) /= 0.0) then
!!$             wert_CP = dNdEx_CP(i) / Norm_CP(i)
!!$          else
!!$             wert_CP = 0.0
!!$          endif
!!$          if (Norm_IMF(i) /= 0.0) then
!!$             wert_IMF = dNdEx_IMF(i) / Norm_IMF(i)
!!$          else
!!$             wert_IMF = 0.0
!!$          endif

          write(40,'(4e15.5)') Ex_bin(i), & 
               & dNdEx(i)/Norm/Ex_step, & 
               & dNdEx_CP(i)/Norm/Ex_step, & 
               & dNdEx_IMF(i)/Norm/Ex_step
       end do
       close(40)
    !-----------------------------------------------------------------------
    endif FinalOutput !-----------------------------------------------------
    !-----------------------------------------------------------------------

  contains 

    subroutine ExInit

      integer :: i

      Norm = float(SubEvents*NumEnsemples)

      dNdEx(:) = 0.0
      dNdEx_CP(:) = 0.0
      dNdEx_IMF(:) = 0.0

!      Norm_CP(:) = 0.0
!      Norm_IMF(:) = 0.0

      do i=1,ExMax
         Ex_bin(i) = float(i)*Ex_step
      end do
      
    end subroutine ExInit

  !*************************************************************************
  end subroutine ExDistributions !******************************************
  !*************************************************************************




  !*************************************************************************
  !****s* FragmentationAnalysis/rapidityDistribution
  ! FUNCTION
  ! Calculates various observables as function of the rapidity.
  ! 
  !*************************************************************************
  subroutine rapidityDistribution(SubEvents,NumEnsemples,SMM_Events, &
                                  irun,FragmentVector,ParticleVector,printflag)

    use typeDefinitions,  only : cluster,particle

    integer,                        intent(in) :: SubEvents,NumEnsemples
    integer,                        intent(in) :: irun,SMM_Events
    type(cluster),  dimension(:,:), intent(in) :: FragmentVector
    type(particle), dimension(:,:), intent(in) :: ParticleVector
    logical,                        intent(in) :: printflag


    integer, parameter :: nrap = 41, ip = 6
    real,    parameter :: dy = 0.1

    real, save,dimension(nrap,0:ip,1:3) :: dndyfra  !rapidity (x,y,z)
    real, save,dimension(nrap)        :: ybin  !rapidity bin

    logical, SAVE :: initFlag2=.true.

    real, SAVE :: NormBUU=0., NormSMM=0.

    integer, dimension(1:3) :: ibin
    real,    dimension(1:3) :: yb
    integer :: i,j,ir,charge,mass
    real    :: k0f,k1f,k2f,k3f
    logical :: PartType
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (initParameters) then
       call GetAnalysisParameters
       initParameters = .false.
    endif

    if (initFlag2) then
       call rapInit
       initFlag2 = .false.
    endif

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! LOOP OVER SMM EVENTS &  PRODUCED SMM PARTICLES 
    ! NORMALIZATION: NormSMM
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Ensemples_Loop1 : do i=1,size(FragmentVector,dim=1)
       Particle_Loop1 : do j=1,size(FragmentVector,dim=2)

          if (FragmentVector(i,j)%ID==0) cycle

          mass     = FragmentVector(i,j)%MassNumber
          charge   = FragmentVector(i,j)%ChargeNumber
          k0f= sqrt( & 
               & FragmentVector(i,j)%mass**2 + & 
               & Dot_Product(FragmentVector(i,j)%momentum(1:3),& 
               & FragmentVector(i,j)%momentum(1:3)))
          k1f      = FragmentVector(i,j)%momentum(1)
          k2f      = FragmentVector(i,j)%momentum(2)
          k3f      = FragmentVector(i,j)%momentum(3)
          PartType = FragmentVector(i,j)%FreeBound


          yb(3) = 0.5*log( (k0f+k3f)/(k0f-k3f) )/yproj
          yb(1) = 0.5*log( (k0f+k1f)/(k0f-k1f) )/yproj
          yb(2) = 0.5*log( (k0f+k2f)/(k0f-k2f) )/yproj

          loop_over_y1 : do ir=1,3 !loop over y_x,y_y,y_z

             if (abs(yb(ir)) .ge. 2.0) cycle

             ibin(ir) = int( yb(ir)/dy + 21. )
             ibin(ir) = min( nrap, max(1,ibin(ir)) )

             !------------------------------------------ charged particles
             if (charge.ge.1) then
                dndyfra(ibin(ir),0,ir) = & 
                     & dndyfra(ibin(ir),0,ir) + charge * NormSMM
             endif
             !------------------------------------------ free protons
             if (charge.eq.1 .and. mass.eq.1) then
                dndyfra(ibin(ir),1,ir) = & 
                     & dndyfra(ibin(ir),1,ir) + NormSMM
             endif
             !------------------------------------------ deuterons
             if (charge.eq.1 .and. mass.eq.2) then   
                dndyfra(ibin(ir),2,ir) = & 
                     & dndyfra(ibin(ir),2,ir) + NormSMM
             endif
             !------------------------------------------ He-isotopes
             if (charge.eq.2) then
                if (mass.eq.3) then
                   dndyfra(ibin(ir),3,ir) = & 
                        & dndyfra(ibin(ir),3,ir) + NormSMM !He3
                endif
                if (mass.eq.4) then
                   dndyfra(ibin(ir),4,ir) = & 
                        & dndyfra(ibin(ir),4,ir) + NormSMM !He4
                endif
                if (mass.eq.6) then
                   dndyfra(ibin(ir),5,ir) = & 
                        & dndyfra(ibin(ir),5,ir) + NormSMM !He6
                endif
             endif
             !------------------------------------------ Li-isotopes
             if (charge.eq.3) then
                dndyfra(ibin(ir),6,ir) = & 
                     & dndyfra(ibin(ir),6,ir) + NormSMM !Li
             endif
          end do loop_over_y1

       end do Particle_Loop1
    end do Ensemples_Loop1

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! LOOP OVER BUU PARTICLES (only emitted ones)
    ! NORMALIZATION: NormBUU
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    i = irun
    Particle_Loop2 : do j=1,size(ParticleVector,dim=2)

       if (ParticleVector(i,j)%ID .eq. 999) cycle !empty
       if (ParticleVector(i,j)%ID .ne. 1)   cycle !only nucleons
!!!       if (ParticleVector(i,j)%SD .ne. 999) cycle !only emitted nucleons --> not used any more

       charge   = ParticleVector(i,j)%Charge
       k0f= sqrt( & 
            & ParticleVector(i,j)%mass**2 + & 
            & Dot_Product(ParticleVector(i,j)%momentum(1:3),& 
            & ParticleVector(i,j)%momentum(1:3)))
       k1f      = ParticleVector(i,j)%momentum(1)
       k2f      = ParticleVector(i,j)%momentum(2)
       k3f      = ParticleVector(i,j)%momentum(3)

       yb(3) = 0.5*log( (k0f+k3f)/(k0f-k3f) )/yproj
       yb(1) = 0.5*log( (k0f+k1f)/(k0f-k1f) )/yproj
       yb(2) = 0.5*log( (k0f+k2f)/(k0f-k2f) )/yproj

       loop_over_y2 : do ir=1,3 !loop over y_x,y_y,y_z

          if (abs(yb(ir)) .ge. 2.0) cycle

          ibin(ir) = int( yb(ir)/dy + 21. )
          ibin(ir) = min( nrap, max(1,ibin(ir)) )

          !------------------------------------------ free protons
          if (charge.eq.1) then
             dndyfra(ibin(ir),0,ir) = & 
                  & dndyfra(ibin(ir),0,ir) + NormBUU
             dndyfra(ibin(ir),1,ir) = & 
                  & dndyfra(ibin(ir),1,ir) + NormBUU
          endif
       end do loop_over_y2

    end do Particle_Loop2


    if (printFlag) then
       open(4,file='RapiditaetsVerteilungen.dat')
       do i=1,nrap
          do j=0,ip
             dndyFra(i,j,1:3) = dndyFra(i,j,1:3)/dy
          end do

          write(4,100)  ybin(i),& 
               &        (dndyFra(i,j,3),j=0,ip),&
               &        (dndyFra(i,j,1),j=0,ip),& 
               &        (dndyFra(i,j,2),j=0,ip)
       end do
100    format(50e12.4)

       close(4)
    endif

  contains 

    !-----------------------------------------------------------------------
    subroutine rapInit

      integer :: iloc1,iloc2
      real :: ystep

      NormBUU = 1./float(SubEvents*NumEnsemples)
      NormSMM = 1./float(SubEvents*NumEnsemples*SMM_Events)

      do iloc1 = 1,nrap
         do iloc2 = 0,ip
            dndyfra(iloc1,iloc2,1:3) = 0.0
         end do
      end do
      ystep = -2.0
      do iloc1=1,nrap
         ybin(iloc1) = ystep
         ystep = ystep + dy
      end do
      do iloc1=1,nrap
         ybin(iloc1) = ybin(iloc1) + dy/2.
      end do

    end subroutine RapInit

    !***********************************************************************
  end subroutine RapidityDistribution !*************************************
  !*************************************************************************


  !*************************************************************************
  !****s* FragmentationAnalysis/spectatorFragmentation
  ! FUNCTION
  ! This routine calculates observables related to fragmentation 
  ! (1) from projectile spectators in low-energy HIC
  !     (ALADIN/INDRA experiments; Spokesmann: W. Trautmann)
  ! (2) from projectile spectators in high-energy HIC
  !     (HypHI experiments; Spokesman: T. Saito)
  ! (3) from decay of residual nuclei in proton-induced reactions
  !     (JPARC-Experiment; Spokesman: T. Saito)
  ! NOTES
  ! Following observables are considered:
  ! * ALADIN/INDRA:
  !   Zbound (charged-weighted sum for SF with Z.ge.2)
  !   IMF (Intermediate Mass Fragments defined as SF in the range: 3.le.Z.le.30)
  !   Energy spectra for light clusters emitted from SF
  !   Isotope- & Isotone-Ratios (d/t,3He/4He,6Li/7Li,...to be continued...)
  ! * HypHI & JPARC:
  !   Rapidity & kinetic energy distributions of fragments with and 
  !   without strangeness degree of freedom are calculated as well.
  !   
  !*************************************************************************
  subroutine spectatorFragmentation(SubEvents,NumEnsemples,SMM_Events,irun,&
                                    FragmentVector,ParticleVector,printflag,Ztarget)

    use typeDefinitions,  only : cluster,particle

    integer,                        intent(in) :: SubEvents,NumEnsemples,irun
    integer,                        intent(in) :: SMM_Events,Ztarget
    type(cluster),  dimension(:,:), intent(in) :: FragmentVector
    type(particle), dimension(:,:), intent(in) :: ParticleVector
    logical,                        intent(in) :: printflag

    integer, parameter :: ipart = 8
    integer, parameter :: ipartZ = 100
    integer, parameter :: Eval  = 80
    real,    parameter :: Estep = 0.0025
    real,    parameter :: yMax  = -6.
    integer, parameter :: nrap  = 60
    real,    parameter :: dy    = 0.1

    !for PANDA:
    real,    parameter :: yMax2  = -2.
    integer, parameter :: nrap2  = 80
    real,    parameter :: dy2    = 0.025

    real,                                  SAVE :: Zbound
    real,                                  SAVE :: IMF
    real, dimension(1:ipart),              SAVE :: Yield
    real, dimension(1:Eval,1:ipart),       SAVE :: dNdEkin
    real, dimension(1:Eval),               SAVE :: Ebin
    real, dimension(-nrap:nrap,1:ipart),   SAVE :: dndysp
    real, dimension(-nrap:nrap,1:2),       SAVE :: dndy_BUU
    real, dimension(-nrap:nrap),           SAVE :: ybin
    !for PANDA:
    real, dimension(-nrap2:nrap2,1:ipartZ),  SAVE :: dndyZ,dndyZ_L,dndyZ_LL 
    real, dimension(-nrap2:nrap2),           SAVE :: dNdY_Fall, dNdY_FL, dNdY_FLL
    real, dimension(-nrap2:nrap2),           SAVE :: ybin2
    logical,                               SAVE :: initFlag3=.true.

!    real, dimension(1:ipart),            SAVE :: YieldL,YieldLL,YieldP
    real, dimension(1:ipart),            SAVE :: YieldL,YieldP
    real, dimension(1:Eval,1:ipart),     SAVE :: dNdEkinL,dNdEkinLL
    real, dimension(1:Eval),             SAVE :: dNdE_Fall, dNdE_FL, dNdE_FLL
!    real, dimension(-nrap:nrap,1:ipart), SAVE :: dndyspL,dndyspLL,dndyspP
    real, dimension(-nrap:nrap,1:ipart), SAVE :: dndyspL,dndyspP

    integer :: i,j,ibin,ibiny,aa,zz,pID,ibiny2,ibinFra
    real    :: k0f,k1f,k2f,k0flrf,k3f,yb,Ekin, EkinFra
    integer :: HypFra
    character(2) :: HypType

    real, SAVE :: NormSMM=0.,NormBUU=0.

!    real, dimension(0:3) :: pin
!    real, dimension(1:3) :: beta_beam
!    real :: veloPF

    real :: PionProb
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (initParameters) then
       call GetAnalysisParameters
       initParameters = .false.
    endif

    if (initFlag3) then
       call SFInit
       if (Ztarget > ipartZ) then
          write(*,*) 'fragmentationAnalysis/VariousDistr.:'
          write(*,*) 'Ztarget > ipartZ !!!'
          write(*,*) 'Ztarget, ipartZ = ', Ztarget, ipartZ
          write(*,*) 'STOP'
          STOP
       endif
       initFlag3 = .false.
    endif

    Ensemples_Loop1 : do i=1,size(FragmentVector,dim=1)
       Particle_Loop1 : do j=1,size(FragmentVector,dim=2)

          if (FragmentVector(i,j)%ID==0) cycle

          k0f= sqrt( & 
               & FragmentVector(i,j)%mass**2 + & 
               & Dot_Product(FragmentVector(i,j)%momentum(1:3),& 
               & FragmentVector(i,j)%momentum(1:3)))
          k3f = FragmentVector(i,j)%momentum(3)

          EkinFra = k0f - FragmentVector(i,j)%mass

          yb  = 0.5*log( (k0f+k3f)/(k0f-k3f) ) 

          zz = FragmentVector(i,j)%ChargeNumber
          aa = FragmentVector(i,j)%MassNumber

          if (zz > ipartZ) then
             write(*,*) 'fragmentationAnalysis / spect.fragmentation:'
             write(*,*) 'zz > ipartZ! zz,ipartZ = ',zz, ipartZ
             write(*,*) 'STOP'
             STOP
          endif

          HypFra  = FragmentVector(i,j)%HypNumber
          HypType = FragmentVector(i,j)%HypType

          !probability of hypernuclei from pion+cluster-->LambdaCluster + Kaon
          PionProb = FragmentVector(i,j)%pionic 

          !spectra in moving frame of projectile (ALADIN experiment)
          k0flrf  = sqrt( & 
               & FragmentVector(i,j)%mass**2 + & 
               & Dot_Product(FragmentVector(i,j)%momentumLRF(1:3),& 
               & FragmentVector(i,j)%momentumLRF(1:3)))
          Ekin    = k0flrf - FragmentVector(i,j)%mass



!          write(*,*) i,j,aa,ekinFra,yb


          if (zz.eq.0) cycle

          !Zbound: charge-weighted sum for Z ge 2
          if (zz.ge.2) Zbound = Zbound + float(zz) * NormSMM

          !IMF (Intermediate Mass Fragments): 3 le Z le 30
          if (zz.ge.3 .and. zz.le.30) IMF = IMF + NormSMM

          ibin = int( Ekin/Estep + 1. )
          ibin = min( Eval, max(1,ibin) )

          ibinFra = int( EkinFra/Estep + 1. )
          ibinFra = min( Eval, max(1,ibin) )

          if (yb < 0.) then
             ibiny = int( yb/dy - 1.)
          else
             ibiny = int( yb/dy)
          endif

          if (yb < 0.) then
             ibiny2 = int( yb/dy2 - 1.)
          else
             ibiny2 = int( yb/dy2)
          endif

          if (zz==1) then !H-like fragments & hyperfragments
             if (aa==2) call momentumSpectra(1)
             if (aa==3) call momentumSpectra(2)
             if (aa==4) call momentumSpectra(3)
             if (aa==5) call momentumSpectra(4)
          end if
          if (zz==2) then !He-like fragments & hyperfragments (only Lambda!)
             if (aa==2) call momentumSpectra(5)
             if (aa==3) call momentumSpectra(6)
             if (aa==4) call momentumSpectra(7)
             if (aa==5) call momentumSpectra(8)
          end if

          !sorting fragments and hyperfragments acc. Z:
          call momentumSpectraZ(zz)

          !spectra of all fragments, single- & double-Lambda fragments:
          call momentumSpectraALL

          if (PionProb > 0.0) then
             if (zz==1) then !H-like hyperfragments from direct pion scattering
                if (aa==2) call momentumSpectra2(1)
                if (aa==3) call momentumSpectra2(2)
                if (aa==4) call momentumSpectra2(3)
                if (aa==5) call momentumSpectra2(4)
             end if
             if (zz==2) then !He-like hyperfragments from direct pion scattering
                if (aa==2) call momentumSpectra2(5)
                if (aa==3) call momentumSpectra2(6)
                if (aa==4) call momentumSpectra2(7)
                if (aa==5) call momentumSpectra2(8)
             end if
          end if

       end do Particle_Loop1
    end do Ensemples_Loop1

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! LOOP OVER BUU PARTICLES (only emitted ones)
    ! NORMALIZATION: NormBUU
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    i = irun
    Particle_Loop2 : do j=1,size(ParticleVector,dim=2)

       if (ParticleVector(i,j)%ID .eq. 999) cycle !empty 

       k0f= sqrt( & 
            & ParticleVector(i,j)%mass**2 + & 
            & Dot_Product(ParticleVector(i,j)%momentum(1:3),& 
            & ParticleVector(i,j)%momentum(1:3)))
       k3f= ParticleVector(i,j)%momentum(3)
       k1f= ParticleVector(i,j)%momentum(1)
       k2f= ParticleVector(i,j)%momentum(2)
       zz = ParticleVector(i,j)%Charge
       pID = ParticleVector(i,j)%ID
       yb = 0.5*log( (k0f+k3f)/(k0f-k3f) )/yproj 

       if (yb < 0.) then
          ibiny = int( yb/dy - 1.)
       else
          ibiny = int( yb/dy)
       endif

       Select Case(pID)

       Case(1) !emitted protons
          if (zz==1) call momentumSpectra_BUU(1)

       Case(32) !Lambdas
          call momentumSpectra_BUU(2)

       end Select


    end do Particle_Loop2


    PrintResults : if (printFlag) then
       open(4,file='SpectatorFragmentation.dat')
       write(4,101) Zbound
       write(4,102) IMF
       if (yield(4).ne.0.) write(4,103) yield(1) / yield(4)
       if (yield(3).ne.0.) write(4,104) yield(2) / yield(3)
       write(4,105) yield(5) / yield(6)
       do i=1,Eval
          write(4,100) Ebin(i),(dNdEkin(i,j)/estep,j=1,ipart)
       end do
       close(4)

       open(14,file='SpectatorFragmentation_dNdY.dat')
       do i=-nrap,nrap
          write(14,100) ybin(i),(dndysp(i,j)/dy,j=1,ipart)
       end do
       close(14)

       open(16,file='SpectatorFragmentationHyp.dat')
       open(15,file='SpectatorFragmentationHyp_dNdY.dat')
       write(15,106) yieldL(:)
!       write(15,107) yieldLL(:)
       write(15,107) yieldP(:)
       do i=-nrap,nrap
!          write(15,100) ybin(i), & 
!               & (dndyspL(i,j)/dy,j=1,ipart), & 
!               & (dndyspLL(i,j)/dy,j=1,ipart)
          write(15,100) ybin(i), & 
               & (dndyspL(i,j)/dy,j=1,ipart), & 
               & (dndyspP(i,j)/dy,j=1,ipart)
       end do
       do i=1,Eval
          write(16,100) Ebin(i), & 
               & (dNdEkinL(i,j)/estep,j=1,ipart), & 
               & (dNdEkinLL(i,j)/estep,j=1,ipart)
       end do
       close(15)

       open(17,file='ParticlesBUU_dNdY.dat')
       do i=-nrap,nrap
          write(17,100) ybin(i), (dndy_BUU(i,j)/dy,j=1,2)
       end do
       
       close(17)

       !PANDA Analysis:
       open(17, file='dNdYZ_PANDA.dat')
       do i=-nrap2,nrap2
          write(17,'(100e15.5)') ybin2(i),(dndyZ(i,j)/dy2,j=1,Ztarget)
       end do
       close(17)
       open(17, file='dNdYZ_LHyp_PANDA.dat')
       do i=-nrap2,nrap2
          write(17,'(100e15.5)') ybin2(i),(dndyZ_L(i,j)/dy2,j=1,Ztarget)
       end do
       close(17)
       open(17, file='dNdYZ_LLHyp_PANDA.dat')
       do i=-nrap2,nrap2
          write(17,'(100e15.5)') ybin2(i),(dndyZ_LL(i,j)/dy2,j=1,Ztarget)
       end do
       close(17)

       open(17, file='dNdY_ALL_PANDA.dat')
       do i=-nrap2,nrap2
          write(17,'(4e15.5)') ybin2(i), dNdY_Fall(i)/dy2, dNdY_FL(i)/dy2, dNdY_FLL(i)/dy2
       end do
       close(17)

       open(17, file='dNdEkin_ALL_PANDA.dat')
       do i=1,Eval
          write(17,'(4e15.5)') Ebin(i), dNdE_Fall(i)/estep, dNdE_FL(i)/estep, dNdE_FLL(i)/estep
       end do
       close(17)


       !End of PANDA Analysis

100    format(50e12.4)
101    format('# <Zbound>       = ',f15.6)
102    format('# <IMF>          = ',f15.6)
103    format('# <d>/<t>        = ',f15.6)
104    format('# <3He>/<4He>    = ',f15.6)
105    format('# <6Li>/<7Li>    = ',f15.6)
106    format('# (3-6)H & (3-6)HeL  = ',8e15.6)
107    format('# (3-6)H & (3-6)HeLL = ',8e15.6)

    endif PrintResults

  contains 

    subroutine SFInit

      integer :: i,j
      real    :: ee, ystep

      NormSMM = 1./float(SubEvents*NumEnsemples*SMM_Events)
      NormBUU = 1./float(SubEvents*NumEnsemples)

      Zbound   = 0.0
      IMF      = 0.0
      yield(:) = 0.0
      yieldL(:) = 0.0
!      yieldLL(:) = 0.0
      yieldP(:) = 0.0
      ee       = 0.0
      do i=1,Eval
         dNdEkin(i,:) = 0.0
         dNdEkinL(i,:)= 0.0
         dNdEkinLL(i,:)= 0.0
         dNdE_Fall(i) = 0.0
         dNdE_FL(i) = 0.0
         dNdE_FLL(i) = 0.0
         Ebin(i)      = ee
         ee           = ee + estep
      end do
      do i=1,Eval
         Ebin(i) = Ebin(i) + estep/2.
      end do

      do i = -nrap,nrap
         do j = 1,ipart
            dndysp(i,j) = 0.0
         end do
         dNdY_Fall(i) = 0.0
         dNdY_FL(i) = 0.0
         dNdY_FLL(i) = 0.0
         dndyspL(i,:) = 0.0
!         dndyspLL(i,:) = 0.0
         dndyspP(i,:) = 0.0
         dndy_BUU(i,:) = 0.0
      end do
      ystep = yMax
      do i=-nrap,nrap
         ybin(i) = ystep
         ystep = ystep + dy
      end do
      do i=-nrap,nrap
         ybin(i) = ybin(i) + dy/2.
      end do

      do i = -nrap2,nrap2
         dNdY_Fall(i) = 0.0
         dNdY_FL(i) = 0.0
         dNdY_FLL(i) = 0.0
         dndyZ(i,:) = 0.0
         dndyZ_L(i,:) = 0.0
         dndyZ_LL(i,:) = 0.0
      end do
      ystep = yMax2
      do i=-nrap2,nrap2
         ybin2(i) = ystep
         ystep = ystep + dy2
      end do
      do i=-nrap2,nrap2
         ybin2(i) = ybin2(i) + dy2/2.
      end do

    end subroutine SFInit

    subroutine momentumSpectraZ(icount)

      integer, intent(in) :: icount

      if (abs(ibiny2) > nrap2) return

      dndyZ(ibiny2,icount) = dndyZ(ibiny2,icount) +  NormSMM

      if(HypFra==1 .and. trim(HypType)=='L') then !single-Lambda hypernuclei
         dndyZ_L(ibiny2,icount) = dndyZ_L(ibiny2,icount) +  NormSMM                
      endif

      if(HypFra==2 .and. trim(HypType)=='LL') then !double-Lambda hypernuclei
         dndyZ_LL(ibiny2,icount) = dndyZ_LL(ibiny2,icount) +  NormSMM                
      endif

    end subroutine momentumSpectraZ

    subroutine momentumSpectra(icount)

      integer, intent(in) :: icount

      if (abs(ibiny) < nrap) dndysp(ibiny,icount) = dndysp(ibiny,icount) +  NormSMM
      dNdEkin(ibin,icount) = dNdEkin(ibin,icount) +  NormSMM

      if(HypFra==1 .and. trim(HypType)=='L') then !single-Lambda hypernuclei
         yieldL(icount)        = yieldL(icount)        +  NormSMM
         dNdEkinL(ibin,icount) = dNdEkinL(ibin,icount) +  NormSMM
         if (abs(ibiny) < nrap) dndyspL(ibiny,icount) = dndyspL(ibiny,icount) +  NormSMM                
      endif

      if(HypFra==2 .and. trim(HypType)=='LL') then !double-Lambda hypernuclei
!         yieldLL(icount)        = yieldLL(icount)        +  NormSMM
         dNdEkinLL(ibin,icount) = dNdEkinLL(ibin,icount) +  NormSMM
!         dndyspLL(ibiny2,icount) = dndyspLL(ibiny,icount) +  NormSMM                
      endif

    end subroutine momentumSpectra

    subroutine momentumSpectraALL

      !all fragments:
      dNdE_Fall(ibinFra) = dNdE_Fall(ibinFra) +  NormSMM

      if(HypFra==1 .and. trim(HypType)=='L') then !single-Lambda hypernuclei
         dNdE_FL(ibinFra) = dNdE_FL(ibinFra) +  NormSMM
      endif

      if(HypFra==2 .and. trim(HypType)=='LL') then !double-Lambda hypernuclei
         dNdE_FLL(ibinFra) = dNdE_FLL(ibinFra) +  NormSMM
      endif

      if (abs(ibiny2) > nrap2) return

      dNdY_Fall(ibiny2) = dNdY_Fall(ibiny2) +  NormSMM

      if(HypFra==1 .and. trim(HypType)=='L') then !single-Lambda hypernuclei
         dNdY_FL(ibiny2) = dNdY_FL(ibiny2) +  NormSMM
      endif

      if(HypFra==2 .and. trim(HypType)=='LL') then !double-Lambda hypernuclei
         dNdY_FLL(ibiny2) = dNdY_FLL(ibiny2) +  NormSMM
      endif

    end subroutine momentumSpectraALL

    subroutine momentumSpectra_BUU(icount)

      integer, intent(in) :: icount

      dndy_BUU(ibiny,icount) = dndy_BUU(ibiny,icount) +  NormBUU

    end subroutine momentumSpectra_BUU


    subroutine momentumSpectra2(icount)

      integer, intent(in) :: icount

      yieldP(icount) = yieldP(icount) +  NormSMM*PionProb
      dndyspP(ibiny,icount) = dndyspP(ibiny,icount) +  NormSMM*PionProb                

    end subroutine momentumSpectra2



  !*************************************************************************
  end subroutine SpectatorFragmentation !***********************************
  !*************************************************************************

  !/////////////////////////////////////////////////////////////////////////
  !/////////////////// ROUTINES RELATED TO CHARMS EXPERIMENTS //////////////
  !/////////////////////////////////////////////////////////////////////////

  !*************************************************************************
  !****s* FragmentationAnalysis/velocityDistr_Charms
  ! FUNCTION
  ! This routine calculates observables related to fragmentation 
  ! (1) from projectile spectators in low-energy HIC
  !     (CHARMS Experiment)
  !   
  !*************************************************************************
  subroutine velocityDistr_Charms(SubEvents,NumEnsemples,SMM_Events,& 
       & FragmentVector,printflag, ProjectileVelo,TypeOfSource)

    use typeDefinitions,  only : cluster

    integer,                        intent(in) :: SubEvents,NumEnsemples
    integer,                        intent(in) :: SMM_Events
    type(cluster),  dimension(:,:), intent(in) :: FragmentVector
    logical,                        intent(in) :: printflag
    real,                           intent(in) :: ProjectileVelo
    integer,                        intent(in) :: TypeOfSource

    real,    parameter :: yMaxs = -5.
    integer, parameter :: nraps = 50
    real,    parameter :: ds    = 0.1

    real, dimension(-nraps:nraps),    SAVE :: ybins  
    real, dimension(-nraps:nraps,1:9),SAVE :: dndyC

    logical,SAVE :: initFlag4=.true.

    integer :: i,j,ibins,aa,zz
    real    :: k0f,k3f,yb,ySRF,betaSRF,ySpect

    real, SAVE :: NormSMM=0.

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (initParameters) then
       call GetAnalysisParameters
       initParameters = .false.
    endif

    if (initFlag4) then
       call CharmsInit
       initFlag4 = .false.
    endif

    ySpect = 0.5*log( (1. + ProjectileVelo)/(1. - ProjectileVelo) ) 
    !test (just to be sure):
    if (ySpect < 0.0 .or. abs(ySpect) < 0.05) then
       write(*,*) 'Module fragmentationAnalysis / routine velocityDistr_Charms:'
       write(*,*) '!!! wrong type of source is analyzed !!!'
       write(*,*) 'Type of the source = ',TypeOfSource
       write(*,*) '!!! STOP !!!'
       STOP
    end if

    Ensemples_Loop1 : do i=1,size(FragmentVector,dim=1)
       Particle_Loop1 : do j=1,size(FragmentVector,dim=2)

          if (FragmentVector(i,j)%ID==0) cycle

          k0f= sqrt( & 
               & FragmentVector(i,j)%mass**2 + & 
               & Dot_Product(FragmentVector(i,j)%momentum(1:3),& 
               & FragmentVector(i,j)%momentum(1:3)))
          k3f = FragmentVector(i,j)%momentum(3)

          yb  = 0.5*log( (k0f+k3f)/(k0f-k3f) ) 

          zz = FragmentVector(i,j)%ChargeNumber
          aa = FragmentVector(i,j)%MassNumber

          !velocity [cm/ns] of projectile fragments in projectile's rest frame (SRF)
          if (yb > 0.0) then !projectile-like fragments only !!!
             ySRF    = yb - ySpect
             betaSRF = (exp(2.*ySRF)-1.)/(exp(2.*ySRF)+1.)*30.

             if (betaSRF < 0.) then
                ibins = int( betaSRF/ds - 1.)
             else
                ibins = int( betaSRF/ds)
             endif

             !          write(*,*) 'ibins, betaSRF',ibins, betaSRF

             if (abs(ibins) .le. nraps) then
                if (aa==6  .and. zz==3)  dndyC(ibins,1) = dndyC(ibins,1) + NormSMM
                if (aa==10 .and. zz==5)  dndyC(ibins,2) = dndyC(ibins,2) + NormSMM
                if (aa==12 .and. zz==6)  dndyC(ibins,3) = dndyC(ibins,3) + NormSMM
                if (aa==16 .and. zz==8)  dndyC(ibins,4) = dndyC(ibins,4) + NormSMM
                if (aa==20 .and. zz==10)  dndyC(ibins,5) = dndyC(ibins,5) + NormSMM
                if (aa==24 .and. zz==12) dndyC(ibins,6) = dndyC(ibins,6) + NormSMM
                if (aa==30 .and. zz==14) dndyC(ibins,7) = dndyC(ibins,7) + NormSMM
                if (aa==41 .and. zz==20) dndyC(ibins,8) = dndyC(ibins,8) + NormSMM
                if (aa==54 .and. zz==25) dndyC(ibins,9) = dndyC(ibins,9) + NormSMM
             end if

          end if

       end do Particle_Loop1
    end do Ensemples_Loop1

    PrintResults : if (printFlag) then

       open(180,file='Rapidity_of_Isotopes.dat')
       do i=-nraps,nraps
          write(180,110) ybins(i),dndyC(i,:)/ds
       end do
       close(180)
110    format(50e12.4)

    endif PrintResults
    
    
  contains 

    subroutine CharmsInit

      integer :: i
      real    :: ysteps

      NormSMM = 1./float(SubEvents*NumEnsemples*SMM_Events)

      ysteps = yMaxs
      do i=-nraps,nraps
         ybins(i) = ysteps
         ysteps = ysteps + ds
      end do
      do i=-nraps,nraps
         ybins(i) = ybins(i) + ds/2.
      end do
      dndyC(:,:)   = 0.0


    end subroutine CharmsInit

  !*************************************************************************
  end subroutine VelocityDistr_Charms !*************************************
  !*************************************************************************

  !*************************************************************************
  !****s* FragmentationAnalysis/CharmsAnalysis
  ! FUNCTION
  ! This routine calculates observables related to spectator fragmentation 
  ! only. It is usefull for a comparison with Experiments of the 
  ! Charms-Collaboration (http://www-win.gsi.de/charms/theses.htm).
  !*************************************************************************
  subroutine CharmsAnalysis(FragmentVector,printflag)

    use typeDefinitions,  only : cluster

    type(cluster),  dimension(:,:), intent(in) :: FragmentVector
    logical, intent(in) :: printflag

    logical, save :: init_Flag=.true.

    integer, save :: iCount
    real,    save, dimension(1:3) :: MeanPos
    real,    save, dimension(0:3) :: MeanVelo
    integer, save :: MeanMass, MeanCharge

    integer :: i,j,aa,zz,MaxMas
    integer, dimension(1:Size(FragmentVector,dim=1)) :: MaxPos
    logical, dimension(1:Size(FragmentVector,dim=1)) :: goodEvent
    real, dimension(0:3) :: kf

    real, dimension(1:3) :: beta_beam
    real, dimension(0:3) :: pin
    !-----------------------------------------------------------------------

    if (initParameters) then
       call GetAnalysisParameters
       initParameters = .false.
    endif

    if (init_Flag) then
       call initRoutine
       init_Flag = .false.
    endif

    !-----------------------------------------------------------------------
    Ensemples_Loop1 : do i=1,size(FragmentVector,dim=1)
       goodEvent(i) = .false.
       MaxMas    = 1000
       MaxPos(i) = 0
       Particle_Loop1 : do j=1,size(FragmentVector,dim=2)

          if (FragmentVector(i,j)%ID==0) cycle

          aa = FragmentVector(i,j)%MassNumber

          if (aa < 35) cycle !only events with A > 35

          goodEvent(i) = .true.

          !find heaviest fragment
          if (MaxMas > aa) then
             MaxMas = aa
             MaxPos(i) = j
          endif
          
       end do Particle_Loop1
    end do Ensemples_Loop1

    !-----------------------------------------------------------------------

    Ensemples_Loop2 : do i=1,size(FragmentVector,dim=1)

       if (.not.goodEvent(i)) cycle

       j = MaxPos(i)

       if (FragmentVector(i,j)%ID==0) then
          write(*,*) 'fragmentationAnalysis/charms-part:'
          write(*,*) 'wrong selection of heaviest fragment...'
          write(*,*) 'Position of Amax= ',j
          write(*,*) '!!! Termination of the programe !!!'
          STOP
       endif

       kf(0)= sqrt( & 
            & FragmentVector(i,j)%mass**2 + & 
            & Dot_Product(FragmentVector(i,j)%momentum(1:3),& 
            & FragmentVector(i,j)%momentum(1:3)))
       kf(1:3) = FragmentVector(i,j)%momentum(1:3)

       zz = FragmentVector(i,j)%ChargeNumber
       aa = FragmentVector(i,j)%MassNumber

       Icount        = iCount        + 1
       MeanMass      = MeanMass      + aa
       MeanCharge    = MeanCharge    + zz
       MeanPos(:)    = meanPos(:)    + FragmentVector(i,j)%position(:)
       MeanVelo(:) = MeanVelo(:) + kf(:)

    end do Ensemples_Loop2
    !-----------------------------------------------------------------------

    PrintResults : if (printFlag) then

       open(100, file='CharmsAnalysis.dat')

       write(100,'(A)') '#################################################'
       write(100,'(A)') '#  Output of analysis for CHARMS-collaboration  #'
       write(100,'(A)') '#  Print mean properties of heaviest fragments  #'
       write(100,'(A)') '#################################################'

       write(100,'(A,1x,f12.6)')  '<A> = ',float( MeanMass ) / float(iCount)
       write(100,'(A,1x,f12.6)')  '<Z> = ',float(MeanCharge) / float(iCount)
       write(100,'(A,1x,3f12.6)') '<X> = ',MeanPos / float(iCount)
       write(100,'(A,1x,3f12.6)') '<V> = ',MeanVelo(1:3)/MeanVelo(0)

       pin(:)         = MeanVelo(:)/float(iCount)
       beta_beam(3)   = uz_Proj / u0_Proj
       beta_beam(1:2) = 0.0
       call boost(beta_beam,pin)

       write(100,'(A,1x,2f12.6)') '*** <V_z>(PF) = ', & 
            & (pin(3) / pin(0)),(pin(3) / pin(0))*30.

    end if PrintResults

    contains

      subroutine initRoutine

        iCount=0
        MeanMass=0
        MeanCharge=0
        MeanVelo=0.
        MeanPos=0.

      end subroutine initRoutine

  !*************************************************************************
  end subroutine CharmsAnalysis !*******************************************
  !*************************************************************************

  !*************************************************************************
  subroutine boost(beta_beam,pin)
  !*************************************************************************
    real,dimension(0:3),intent(inout) :: pin
    real,dimension(1:3),intent(in) :: beta_beam
    real  :: gamma
    real  :: betaFour

    !Evaluate gamma
    gamma = 1.0 - Dot_Product(beta_beam,beta_beam)
    if(gamma .gt. 0.0) then
       gamma = 1.0/sqrt(gamma)
    else
       write(*,*)
       write(*,*)'(1-beta**2) in lorentz less or equal zero: ', gamma
       write(*,*)'beta=',  beta_beam
       write(*,*)'Stop program'
       stop
    end if
    !beta*pin(1:3)
    betaFour = Dot_Product(beta_beam, pin(1:3))
    !do transformation
    pin(1:3)= pin(1:3) + & 
         & gamma*beta_beam(1:3)*(gamma/(gamma+1.)*betaFour-pin(0))
    pin(0)= gamma*(pin(0)-betaFour)

  !*************************************************************************
  end subroutine boost !****************************************************
  !*************************************************************************




end module fragmentationAnalysis
