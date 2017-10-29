!******************************************************************************
!****m* /HeavyIonAnalysis
! NAME
! module HeavyIonAnalysis
!
! PURPOSE
! Contains output routines for heavy ion collisions.
!******************************************************************************
module HeavyIonAnalysis

  implicit none
  private

  !****************************************************************************
  !****g* HeavyIonAnalysis/flag_outputReal
  ! PURPOSE
  ! If .true., then the output of the real particle vector
  ! will be written to the file 'DoHIA.dat'.
  ! SOURCE
  !
  logical, save :: flag_outputReal = .false.
  !****************************************************************************


  !****************************************************************************
  !****g* HeavyIonAnalysis/flag_outputPert
  ! PURPOSE
  ! If .false., then the output of the perturbative particle vector
  ! will be written to the file 'DoHIA_pert.dat'.
  ! SOURCE
  !
  logical, save :: flag_outputPert = .false.
  !****************************************************************************


  !****************************************************************************
  !****g* HeavyIonAnalysis/flag_outputDetailed
  ! PURPOSE
  ! Print out more detailed information at each time step
  ! from subroutine HeavyIon_evol:
  ! * rhorad_*.dat
  ! * rhoz_*.dat
  ! * rhozx_*.dat
  ! * Fields_*.dat
  ! * pauli_*.dat
  ! * dens_max.dat
  ! SOURCE
  !
  logical, save :: flag_outputDetailed = .false.
  !****************************************************************************


  !****************************************************************************
  !****g* HeavyIonAnalysis/pionAnalysis
  ! PURPOSE
  ! This flag generates various pion spectra (p_T, m_T, y, etc).
  ! The analysis operates under the assumption of a fixed target,
  ! and expects the collision to be performed in the CMS system
  ! (cf. cmsFlag in namelist /heavyIon/).
  ! The analysis matches the one applied to the HADES data in
  ! Agakishiev et al., Eur.Phys.J. A40 (2009) 45-49.
  ! SOURCE
  !
  logical, save :: pionAnalysis = .false.
  !****************************************************************************


  !****************************************************************************
  !****g* HeavyIonAnalysis/rapBinning
  ! PURPOSE
  ! Rapidity binning for the pion analysis (only used if pionAnalysis = .true.).
  ! The numbers represent the binning borders in y0. For each of the seven
  ! y0 bins, a separate mT spectrum will be generated.
  ! SOURCE
  !
  real, dimension(0:7), save :: rapBinning = (/ -0.75, -0.45, -0.15, 0.15, 0.45, 0.75, 1.05, 1.35 /)
  !****************************************************************************


  !****************************************************************************
  !****g* HeavyIonAnalysis/KaonAnalysis
  ! PURPOSE
  ! This flag generates various Kaon spectra and Kaon-related analyses.
  ! SOURCE
  !
  logical, save :: KaonAnalysis = .false.
  !****************************************************************************

  !****************************************************************************
  !****g* HeavyIonAnalysis/DensityPlot
  ! PURPOSE
  ! This flag select printing the density for several time steps
  ! SOURCE
  !
  logical, save :: DensityPlot = .false.
  !****************************************************************************

  !****************************************************************************
  !****g* HeavyIonAnalysis/NucleonMassPlot
  ! PURPOSE
  ! This flag select printing the (invariant) mass of the nucleons for several
  ! time steps
  ! SOURCE
  !
  logical, save :: NucleonMassPlot = .false.
  !****************************************************************************

  !***************************************************************************
  !****g* HeavyIonAnalysis/do_Tmunu
  ! SOURCE
  logical,save :: do_Tmunu=.false.
  ! PURPOSE
  ! Switch for Tmunu output.
  !***************************************************************************

  !***************************************************************************
  !****g* HeavyIonAnalysis/rotateZ_Tmunu
  ! SOURCE
  logical,save :: rotateZ_Tmunu=.false.
  ! PURPOSE
  ! select, whether the particles are first rotated to be aligned to the
  ! z-axis
  !***************************************************************************

  !***************************************************************************
  !****g* HeavyIonAnalysis/correctPot_Tmunu
  ! SOURCE
  integer,save :: correctPot_Tmunu = 0
  ! PURPOSE
  ! select, whether the energy is corrected for the potential or not:
  ! * 0: no correction
  ! * 1: full potential added to p0
  ! * 2: only U_b/2+U_r added to p0
  !***************************************************************************

  logical, save :: initFlag=.true.


  real, dimension(:), allocatable, save :: &
       arrMultSet2, arrMultSet4, arrMultSet5

  ! string constants may be broken over multiple continuation lines:
  character(*), parameter :: Form5 = &
       "(i4,1x,i2,1x,i8,2(1x,i4),2x,f6.3,1x,i1,3(1x,f8.3),3(1x,f8.3)&
       &1x,i5,1x,i4,1x,f5.2)"
  character(*), parameter :: Form6 = &
       "(i4,1x,i2,1x,i10,2(1x,i4),2x,f8.5,1x,i1,4(1x,f10.5)&
       &1x,e13.6,1x,i5,1x,i4,1x,f5.2)"

  public :: DoHeavyIonAnalysis, DoHeavyIonAnalysisTime, HeavyIon_evol

contains


  !****************************************************************************
  !****s* HeavyIonAnalysis/init
  ! NAME
  ! subroutine init
  ! PURPOSE
  ! Reads in the namelist "HICanalysis_Input"
  ! INPUTS
  ! * (none)
  ! OUTPUT
  ! * Initializes global module variables
  !****************************************************************************
  subroutine init
    use output, only: Write_ReadingInput

    integer :: ios

    !**************************************************************************
    !****n* HeavyIonAnalysis/HICanalysis_Input
    ! NAME
    ! NAMELIST /HICanalysis_Input/
    ! PURPOSE
    ! Includes the switches:
    ! * flag_outputReal
    ! * flag_outputPert
    ! * flag_outputDetailed
    ! * pionAnalysis
    ! * rapBinning
    ! * KaonAnalysis
    ! * DensityPlot
    ! * NucleonMassPlot
    ! * do_Tmunu
    ! * rotateZ_Tmunu
    ! * correctPot_Tmunu
    !**************************************************************************
    NAMELIST /HICanalysis_Input/ &
         flag_outputReal, flag_outputPert, flag_outputDetailed, &
         pionAnalysis, rapBinning, KaonAnalysis, &
         DensityPlot, NucleonMassPlot, &
         do_Tmunu, rotateZ_Tmunu, correctPot_Tmunu

    call Write_ReadingInput('HICanalysis_Input',0)
    rewind(5)
    read(5,nml=HICanalysis_Input,iostat=ios)
    call Write_ReadingInput('HICanalysis_Input',0,ios)

    write(*,*) 'flag_outputReal    : ', flag_outputReal
    write(*,*) 'flag_outputPert    : ', flag_outputPert
    write(*,*) 'flag_outputDetailed: ', flag_outputDetailed
    write(*,*) 'pionAnalysis       : ', pionAnalysis
    if (pionAnalysis) write(*,'(A,8F7.3)') ' rapBinning         : ', rapBinning
    write(*,*) 'KaonAnalysis       : ', KaonAnalysis
    write(*,*) 'DensityPlot        : ', DensityPlot
    write(*,*) 'NucleonMassPlot    : ', NucleonMassPlot
    write(*,*) 'do Tmunu           : ', do_Tmunu
    if (do_Tmunu) then
       write(*,*) '  Tmunu: rotateZ   : ', rotateZ_Tmunu
       write(*,*) '  Tmunu: correctPot: ', correctPot_Tmunu
    end if

    call Write_ReadingInput('HICanalysis_Input',1)

    initFlag=.false.

  end subroutine init


  !****************************************************************************
  !****s* HeavyIonAnalysis/DoHeavyIonAnalysis
  ! NAME
  ! subroutine DoHeavyIonAnalysis(realParts,pertParts,finalFlag)
  !
  ! PURPOSE
  ! Makes the output of the test particle id's, charges, masses, positions and
  ! momenta on disk.
  ! Also does some statistical analysis (see subroutine histo1).
  !
  ! INPUTS
  ! * type(particle), dimension(:,:), :: realParts -- real particle vector
  ! * type(particle), dimension(:,:), :: pertParts -- perturb. particle vector
  ! * logical, intent (in)     :: finalFlag -- if .true., it is the last call
  !   for one specific energy, therefore final output must be made.
  !
  ! RESULT
  ! The test particle infos are printed into file 'DoHIA.dat'.
  ! The statistical results are printed to the files:
  ! * 'DoHIA1.dat'
  ! * 'DoHIA2.dat'
  !****************************************************************************
  subroutine DoHeavyIonAnalysis(realParts, pertParts, finalFlag)

    use IdTable, only: nucleon, pion, EOV, NOP, isMeson
    use particleDefinition
    use initHeavyIon, only: b_HI => b
    use initHadron, only: b_had => b, particleId, antiparticle, perturbative
    use initElementary, only: b_ele => impactParameter
    use twoBodyStatistics, only: sqrts_distribution, rate
    use inputGeneral, only: eventtype
    use eventtypes, only: HeavyIon, hadron, elementary
    use history, only: history_getParents

    type(particle), dimension(:,:), intent(in), target :: realParts
    type(particle), dimension(:,:), intent(in), target :: pertParts
    logical,                        intent(in)         :: finalFlag

    integer :: i,j,k,nEns,nPart,indFree,parents(1:3)
    real :: factor,stossParameter
    integer, save :: isu=0     ! counter of subsequent runs
    type(particle), dimension(1:2) :: dummy

    type(particle), POINTER :: pPart

    if (initFlag) call init

    isu = isu + 1

    select case (eventtype)
    case (HeavyIon)
       stossParameter = b_HI
    case (hadron)
       stossParameter = b_had
    case (elementary)
       stossParameter = b_ele
    case default
       write(*,*) ' Problem in DoHeavyIonAnalysis: ', &
            'impact parameter is not defined'
       stossParameter=0.
    end select

    nEns = size(realParts,dim=1)
    nPart = size(realParts,dim=2)


    if (flag_outputReal) then

       !***********************************************************************
       !****o* HeavyIonAnalysis/DoHIA.dat
       ! NAME
       ! file DoHIA.dat
       ! PURPOSE
       ! Contains the full dump of the real particle vector at the end of the
       ! simulation, including particle ID, charge, position, momentum, etc.
       !
       ! Please note, that this file is in some old-fashioned output style.
       !
       ! DEPRECATED. may be deleted sooner or later.
       !***********************************************************************
       open(30,file='DoHIA.dat',position='Append')

       do i = 1,nEns
          do j = 1,nPart

             pPart => realParts(i,j)
             if (pPart%ID <= 0) cycle

             factor = merge(-1., 1.,pPart%antiparticle)
             indFree = merge(1, 0, isFree(realParts,i,j,pPart))

             parents = history_getParents(pPart%history)
             do k=1,2
                if (parents(k)>200) parents(k)=200-parents(k)
             end do

             write(30,Form5) &
                  int(factor)*pPart%ID,pPart%charge,&
                  pPart%event(1),parents(1:2),pPart%mass,&
                  indFree,pPart%position(1:3),&
                  pPart%momentum(1:3),&
                  i,isu,stossParameter


          end do

       end do

       close(30)
    end if


    if (flag_outputPert) then

       !***********************************************************************
       !****o* HeavyIonAnalysis/DoHIA_pert.dat
       ! NAME
       ! file DoHIA_pert.dat
       ! PURPOSE
       ! Contains the full dump of the perturbative particle vector at the end
       ! of the simulation, including particle ID, charge, position, momentum,
       ! etc.
       !***********************************************************************
       open(30,file='DoHIA_pert.dat',position='Append')
       do i = 1,nEns
          do j = 1,size(pertParts,dim=2) ! nPart for pert vector

             pPart => pertParts(i,j)
             if (pPart%ID == EOV) exit
             if (pPart%ID == NOP) cycle
             if (pPart%ID == nucleon) cycle

             factor = merge(-1., 1.,pPart%antiparticle)
             indFree = merge(1, 0, isFree(realParts,i,j,pPart))

             parents = history_getParents(pPart%history)
             do k=1,2
                if (parents(k)>200) parents(k)=200-parents(k)
             end do

             write(30,Form6) &
                  int(factor)*pPart%ID,pPart%charge,&
                  pPart%event(1),parents(1:2),pPart%mass,&
                  indFree, (pPart%momentum(k), k=0,3),&
                  pPart%perWeight,i,isu,stossParameter

          end do

       end do

       close(30)

    end if

    if (finalFlag) then
       call sqrts_distribution(dummy,0,.true.)
       call rate(dummy,dummy,0.,.true.)
       call Spectra
       call countParts(realParts, -99.9, "final_")
    end if

    if (eventtype==hadron &
         .and. particleId==nucleon &
         .and. antiparticle &
         .and. .not.perturbative ) then
       call histo1
    end if

    if (pionAnalysis) call analyze_Pions
    if (KaonAnalysis) call analyze_Kaons


  contains

    !**************************************************************************
    !***is* DoHeavyIonAnalysis/histo1
    !**************************************************************************
    subroutine histo1

      integer, parameter :: N_max=200  ! Max number of pions produced in event
      real, save :: P_Npion(0:N_max)   ! Pion multiplicity distribution
      integer, parameter :: Nmom=100   ! Number of momentum bins
      real, parameter :: dmom=0.02     ! Momentum bin (GeV/c)
      real, save, dimension(1:Nmom,0:10) :: dNpiondMom  ! Pion momentum distr.
      real, save :: pion_events

      integer :: numPions,ibin
      logical :: flag_not_only_pions
      real :: fnorm,pinumAv,momentumAbs
      type(particle), POINTER :: pPart


      if (isu.eq.1) then
         pion_events=0.
         P_Npion=0.
         dNpiondMom=0.
         open(36,file='FewPionEvents.dat')
      end if

      Ensemble_loop : do i = 1,nEns

         numPions=0
         flag_not_only_pions=.false.

         Particle_loop1 : do j = 1,size(realParts,dim=2)

            if (realParts(i,j)%ID <= 0) cycle Particle_loop1

            if (realParts(i,j)%ID.eq.pion) then
               numPions=numPions+1
            else if (isMeson(realParts(i,j)%ID)) then
               flag_not_only_pions=.true.
            end if

         end do Particle_loop1


         if (.not.flag_not_only_pions) then

            pion_events=pion_events+1.
            if (numPions.le.N_max)  P_Npion(numPions)=P_Npion(numPions)+1.

            Particle_loop2 : do j = 1,size(realParts,dim=2)

               pPart => realParts(i,j)
               if (pPart%ID <= 0) cycle Particle_loop2

               if (numPions.le.1) then
                  factor = merge( -1., 1., pPart%antiparticle)

                  write(36,5) int(factor)*pPart%ID,pPart%charge,&
                       pPart%mass,&
                       pPart%position(1:3),&
                       pPart%momentum(1:3),&
                       i,isu
5                 format(i4,1x,i2,1x,f6.3,3(1x,f8.3),3(1x,f8.3),1x,i5,1x,i4)
               end if

               if (pPart%ID .ne. pion) cycle Particle_loop2


               momentumAbs=sqrt(dot_product(pPart%momentum(1:3),pPart%momentum(1:3)))
               ibin=nint((momentumAbs-dmom/2.)/dmom)+1
               if (ibin.ge.1 .and. ibin.le.Nmom .and. pPart%charge.ne.0) then
                  dNpiondMom(ibin,0)=dNpiondMom(ibin,0)+1.
                  if (numPions.ge.1 .and. numPions.le.10) &
                       & dNpiondMom(ibin,numPions)=dNpiondMom(ibin,numPions)+1.
               end if

            end do Particle_loop2

         end if


      end do Ensemble_loop

      if (finalFlag) then

         ! Normas:
         !fnorm=1./float(nEns)/float(isu)
         if (pion_events.gt.0.) then
            fnorm=1./pion_events
         else
            fnorm=0.
         end if

         dNpiondMom(:,:)=dNpiondMom(:,:)*fnorm/dmom
         P_Npion(:)=P_Npion(:)*fnorm


         open(36,file='DoHIA1.dat')
         write(36,*)'# Distribution of events in pion number'
         write(36,*)'# Number of events: ', nEns*isu
         write(36,*)'# Npion:    P_Npion:'
         pinumAv=0.
         do j=0,N_max
            pinumAv=pinumAv+float(j)*P_Npion(j)
            write(36,*) j, P_Npion(j)
         end do
         write(36,*)'# Norma: ', sum(P_Npion(:))
         write(36,*)'# Average pion number:', pinumAv
         close(36)

         open(36,file='DoHIA2.dat')
         write(36,*)'# Charged pion momentum distribution'
         write(36,*)'# Number of events: ', nEns*isu
         write(36,*)'# momentum, GeV/c:    dNpiondMom, c/GeV:'
         do j=1,Nmom
            write(36,'(f6.3,8(1x,e10.3))')  (float(j)-0.5)*dmom, dNpiondMom(j,0), &
                 & dNpiondMom(j,1:6), sum(dNpiondMom(j,7:10))
         end do
         write(36,*)'# Norma: ', sum(dNpiondMom(:,0))*dmom
         close(36)

      end if

    end subroutine histo1

    !**************************************************************************
    !***is* DoHeavyIonAnalysis/analyze_Pions
    !**************************************************************************
    subroutine analyze_Pions
      use histMC
      use IDTable, only: EOV
      use minkowski, only: abs4, abs3
      use inputGeneral, only: num_Runs_SameEnergy, num_energies
      use initHeavyIon, only: ekin_lab_Projectile
      use constants, only: mN
      use nucleusDefinition
      use nucleus, only: getProjectile
      use histMC_avg

      type(histogramMC), save :: hist_pT, hist_y, hist_y0, hist_cost
      type(histogramMC), dimension(0:7), save :: hist_mT       ! 0=total; 1:7=different rap. bins
      type(histogramMC_avg), save :: hist_v1_y, hist_v2_y, hist_v1_u, hist_v2_u
      logical, save :: first = .true.
      real, save :: w, y_cms, u_proj
      integer :: i, j, ch, rb
      real :: mom(0:3), pT, mT, m, y, y0, cost, E, p, pabs, v1, v2, ut0, beta_proj, beta(1:3), gamma
      character(40) :: str
      type(tnucleus), pointer :: proj
      integer, save :: n = 0
      real :: fac

      if (first) then
         call CreateHistMC(hist_pT, 'Pion transverse momentum spectra: dN/dpT',&
              0.,1.,0.01,3)
         call CreateHistMC(hist_y,  'Pion rapidity spectra: dN/dy', &
              -3.,3.,0.1 ,3)
         call CreateHistMC(hist_y0, 'Pion rapidity spectra: dN/dy0', &
              -3.,3.,0.1 ,3)
         call CreateHistMC(hist_mT, 'Pion transverse mass spectra: mT^(-2)*dN/dmT', &
              0.,1.,0.01,3)
         call CreateHistMC(hist_cost, 'Pion polar angle spectra: dN/dcos(theta) with 0.2<p<0.8', &
              -1.,1.,0.05,3)
         call CreateHistMC_avg(hist_v1_y, 'Directed flow of pions: v1(y0)', &
              -2.,2.,0.2 ,3)
         call CreateHistMC_avg(hist_v2_y, 'Elliptic flow of pions: v2(y0)', &
              -2.,2.,0.2 ,3)
         call CreateHistMC_avg(hist_v1_u, 'Directed flow of pions: v1(ut0)', &
              0.,4.,0.2 ,3)
         call CreateHistMC_avg(hist_v2_u, 'Elliptic flow of pions: v2(ur0)', &
              0.,4.,0.2 ,3)
         do i=1,7
            write(str,'(A,i1,A,f5.2,A,f5.2,A)') &
                 " (rap. bin #",i,": y0 = ",rapBinning(i-1),&
                 " ... ",rapBinning(i)," GeV)"
            hist_mT(i)%name = trim(hist_mT(i)%name) // str
         end do
         hist_pT%xDesc = 'p_T [GeV]'
         hist_pT%yDesc = (/ "pi-", "pi0", "pi+" /)
         call CopyDesc(hist_y, hist_pT)
         call CopyDesc(hist_y0, hist_pT)
         call CopyDesc(hist_cost, hist_pT)
         call CopyDesc(hist_mT, hist_pT)
         hist_y%xDesc    = 'y'
         hist_y0%xDesc   = 'y0'
         hist_cost%xDesc = 'cos(theta_cm)'
         hist_mT%xDesc   = 'm_T - m [GeV]'
         hist_v1_y%xDesc   = 'y0'
         hist_v1_y%yDesc = (/ "pi-", "pi0", "pi+" /)
         call CopyDesc_avg(hist_v2_y, hist_v1_y)
         call CopyDesc_avg(hist_v1_u, hist_v1_y)
         call CopyDesc_avg(hist_v2_u, hist_v1_y)
         hist_v1_u%xDesc   = 'ut0'
         hist_v2_u%xDesc   = 'ut0'
         w = 1. / (nEns*num_Runs_SameEnergy*num_energies)   ! weight factor
         E = 2*mN + ekin_lab_Projectile
         p = sqrt(ekin_lab_Projectile**2 + 2.*mN*ekin_lab_Projectile)  ! assumption: fixed target
         y_cms = 0.5*log((E+p)/(E-p))
         write(*,*) "analyze_Pions: y_cms = ", y_cms
         proj => getProjectile()
         beta_proj = sqrt(sum(proj%velocity**2))
         u_proj = beta_proj / sqrt(1. - beta_proj**2)
         write(*,*) "analyze_Pions: u_proj = ", u_proj
         first = .false.
      end if

      do i=1,nEns
         do j=1,nPart
            if (realParts(i,j)%ID == EOV) exit
            if (realParts(i,j)%ID /= pion) cycle

            ch = realParts(i,j)%charge+2
            mom = realParts(i,j)%momentum
            beta = mom(1:3) / mom(0)
            gamma = 1./ sqrt( 1. - dot_product(beta,beta) )
            m = abs4(mom)                                 ! mass
            pabs = abs3(mom)                              ! absolute momentum
            y = 0.5*log((mom(0)+mom(3))/(mom(0)-mom(3)))  ! rapidity
            y0 = y/y_cms                                  ! 'normalized' rapidity
            pt = sqrt(mom(1)**2+mom(2)**2)                ! transverse momentum
            mt = sqrt(m**2+pt**2)                         ! transverse mass
            cost = mom(3)/pabs
            v1 = mom(1)/pt                       ! v1 = < px/pT >
            v2 = (mom(1)**2 - mom(2)**2)/pt**2   ! v2 = < (px**2-py**2)/pT**2 >
            ut0 = gamma * sqrt(beta(1)**2+beta(2)**2) / u_proj  ! transverse comp. of scaled four-velocity

            call AddHistMC(hist_pT,    pt,   ch, w)
            call AddHistMC(hist_y,     y,    ch, w)
            call AddHistMC(hist_y0,    y0,   ch, w)
            call AddHistMC(hist_mT(0), mt-m, ch, w/mt**2)
            if (pabs>0.2 .and. pabs<0.8) &
                 call AddHistMC(hist_cost,  cost, ch, w)
            if (ut0>1.0 .and. ut0<4.2) then
              call AddHistMC_avg(hist_v1_y, y0, ch, v1)
              call AddHistMC_avg(hist_v2_y, y0, ch, v2)
            end if
            if (y0>-1.8 .and. y0<0.) then
              call AddHistMC_avg(hist_v1_u, ut0, ch, v1)
              call AddHistMC_avg(hist_v2_u, ut0, ch, v2)
            end if
            ! determine rap. bin
            rb = 0
            do while (rb < 8)
               if (y0 < rapBinning(rb)) exit
               rb = rb + 1
            end do
            if (rb>0 .and. rb<8) &
                 call AddHistMC(hist_mT(rb), mt-m, ch, w/mt**2)
         end do
      end do

      n = n + 1                                        ! count runs
      fac = num_Runs_SameEnergy*num_energies/float(n)  ! multiplication factor for writing histograms

      call WriteHistMC(hist_pT,  'PionPt.dat',   mul=fac)
      call WriteHistMC(hist_y,   'PionY.dat',    mul=fac)
      call WriteHistMC(hist_y0,  'PionY0.dat',   mul=fac)
      call WriteHistMC(hist_cost,'PionCost.dat', mul=fac)
      do i=0,7
         str = ""
         if (i>0) write(str,'(A,i1)') "_rapBin",i
         call WriteHistMC(hist_mT(i), 'PionMt'//trim(str)//'.dat', mul=fac)
      end do
      call WriteHistMC_avg(hist_v1_y,'PionV1_y.dat')
      call WriteHistMC_avg(hist_v2_y,'PionV2_y.dat')
      call WriteHistMC_avg(hist_v1_u,'PionV1_u.dat')
      call WriteHistMC_avg(hist_v2_u,'PionV2_u.dat')
    end subroutine analyze_Pions

    !**************************************************************************
    !***is* DoHeavyIonAnalysis/analyze_Kaons
    !**************************************************************************
    subroutine analyze_Kaons
      use histMC
      use initHeavyIon, only: ekin_lab_Projectile
      use idTable, only: Kaon, Kaonbar, isBaryon
      use constants, only: mN
      use minkowski, only: abs4, abs3
      use inputGeneral, only: num_Runs_SameEnergy, num_energies

      type(histogramMC), save :: hist_p, hist_pT, hist_y, hist_y0, &
           hist_cost, hist_mt, hist_parents
      logical, save :: first = .true.
      real, save :: w, y_cms
      integer :: i, j, ch, parents(1:3), channel
      real :: mom(0:3), pT, mT, m, y, y0, cost, E, p, pabs

      if (first) then
         call CreateHistMC(hist_p,  'Kaon momentum spectra: dN/dp', &
              0.,1.5,0.05,4)
         call CreateHistMC(hist_pT, 'Kaon transverse momentum spectra: dN/dpT',&
              0.,1.,0.01,4)
         call CreateHistMC(hist_y,  'Kaon rapidity spectra: dN/dy', &
              -3.,3.,0.1 ,4)
         call CreateHistMC(hist_y0, 'Kaon rapidity spectra: dN/dy0', &
              -3.,3.,0.1 ,4)
         call CreateHistMC(hist_mT, 'Kaon transverse mass spectra: mT^(-2)*dN/dmT', &
              0.,1.,0.01,4)
         call CreateHistMC(hist_cost, 'Kaon polar angle spectra: dN/dcos(theta)', &
              -1.,1.,0.05,4)
         call CreateHistMC(hist_parents, 'Kaon parent sources', &
              0.5,5.5,1.,4)
         hist_p%xDesc = 'p [GeV]'
         hist_p%yDesc = (/ "K0   ", "K+   ", "K-   ", "K0bar" /)
         call CopyDesc(hist_pT,      hist_p)
         call CopyDesc(hist_y,       hist_p)
         call CopyDesc(hist_y0,      hist_p)
         call CopyDesc(hist_cost,    hist_p)
         call CopyDesc(hist_mT,      hist_p)
         call CopyDesc(hist_parents, hist_p)
         hist_pT%xDesc   = 'p_T [GeV]'
         hist_y%xDesc    = 'y'
         hist_y0%xDesc   = 'y0'
         hist_cost%xDesc = 'cos(theta_cm)'
         hist_mT%xDesc   = 'm_T - m [GeV]'
         hist_parents%xDesc = "source channel"
         w = 1. / (nEns*num_Runs_SameEnergy*num_energies)   ! weight factor
         E = 2*mN + ekin_lab_Projectile
         p = sqrt(ekin_lab_Projectile**2 + 2.*mN*ekin_lab_Projectile)  ! assumption: fixed target
         y_cms = 0.5*log((E+p)/(E-p))
         write(*,*) "analyze_Kaons: y_cms = ",y_cms
         first = .false.
      end if

      do i=1,nEns
         do j=1,nPart
            if (realParts(i,j)%ID == EOV) exit
            if (realParts(i,j)%ID == Kaon) then
               ch = realParts(i,j)%charge+1
            else if (realParts(i,j)%ID == Kaonbar) then
               ch = realParts(i,j)%charge+4
            else
               cycle
            end if

            mom = realParts(i,j)%momentum
            m = abs4(mom)                                 ! mass
            pabs = abs3(mom)                              ! absolute momentum
            y = 0.5*log((mom(0)+mom(3))/(mom(0)-mom(3)))  ! rapidity
            y0 = y/y_cms                                  ! 'normalized' rapidity
            pt = sqrt(mom(1)**2+mom(2)**2)                ! transverse momentum
            mt = sqrt(m**2+pt**2)                         ! transverse mass
            cost = mom(3)/pabs
            call AddHistMC(hist_p,    pabs, ch, w)
            call AddHistMC(hist_pT,     pt, ch, w)
            call AddHistMC(hist_y,       y, ch, w)
            call AddHistMC(hist_y0,     y0, ch, w)
            call AddHistMC(hist_mT,   mt-m, ch, w/mt**2)
            call AddHistMC(hist_cost, cost, ch, w)
            ! determine where it came from
            parents = history_getParents (realParts(i,j)%history)
            if (parents(2) > 0) then
              ! 2-body collisions
              if (isBaryon(parents(1)) .and. isBaryon(parents(2))) then
                channel = 1  ! channel 1: BB
              else if (isMeson(parents(1)) .and. isMeson(parents(2))) then
                channel = 2  ! channel 2: mm
              else
                channel = 3  ! channel 3: mB
              end if
!               print *,"analyze_Kaons (2-body):", ch, parents(1:3)
            else if (isMeson(parents(1))) then
              channel = 4   ! channel 4: meson decay
!               print *,"analyze_Kaons (meson dec.):", ch, parents(1:3)
            else if (isBaryon(parents(1))) then
              channel = 5   ! channel 5: baryon decay
!               print *,"analyze_Kaons (baryon dec.):", ch, parents(1:3)
            else
              write(*,*) "analyze_Kaons: unknown parents!", ch, parents(1:3), realParts(i,j)%history
              stop
            end if
            call AddHistMC(hist_parents, float(channel), ch, w)
         end do
      end do
      call WriteHistMC(hist_p,       'KaonPlab.dat')
      call WriteHistMC(hist_pT,      'KaonPt.dat')
      call WriteHistMC(hist_y,       'KaonY.dat')
      call WriteHistMC(hist_y0,      'KaonY0.dat')
      call WriteHistMC(hist_mT,      'KaonMt.dat')
      call WriteHistMC(hist_cost,    'KaonCost.dat')
      call WriteHistMC(hist_parents, 'KaonParents.dat')
    end subroutine analyze_Kaons


    !**************************************************************************
    !***is* DoHeavyIonAnalysis/Spectra
    !**************************************************************************
    subroutine Spectra

      use histMPf90
      use IDtable, only: isHadron

      type(histogramMP), save :: hMP_ESet2, hMP_ESet4, hMP_ESet5
      type(histogramMP), save :: hMP_mTSet2, hMP_mTSet4, hMP_mTSet5

      integer :: i,j, nEns,nPart
      real :: mulFak, mom, mom0, mT, y
      type(particle), POINTER :: pPart
      real, parameter :: dy = 0.05

      call CreateHistMP(hMP_ESet2, "dN/pE dE", 0.0, 2.5, 0.02, 2)
      call CreateHistMP(hMP_ESet4, "dN/pE dE", 0.0, 2.5, 0.02, 4)
      call CreateHistMP(hMP_ESet5, "dN/pE dE", 0.0, 2.5, 0.02, 5)

      call CreateHistMP(hMP_mTSet2, "dN/mT2 dmT dy", 0.0, 2.5, 0.02, 2)
      call CreateHistMP(hMP_mTSet4, "dN/mT2 dmT dy", 0.0, 2.5, 0.02, 4)
      call CreateHistMP(hMP_mTSet5, "dN/mT2 dmT dy", 0.0, 2.5, 0.02, 5)

      nEns  = size(realParts,dim=1)
      nPart = size(realParts,dim=2)
      mulfak = 1.0/(nEns)

      do i=1,nEns
         do j=1,nPart
            pPart => realParts(i,j)
            if (pPart%Id <  0) exit
            if (pPart%Id <= 0) cycle

            mom = absMom(pPart)
            mom0 = pPart%momentum(0)
            mT = sqrt(max(1e-15,pPart%momentum(0)**2-pPart%momentum(3)**2))
            y = rapidity(pPart)

            call AddHistMP(hMP_ESet2, pPart, mom0, 1.0/(mom0*mom), 1.0)
            call AddHistMP(hMP_ESet4, pPart, mom0, 1.0/(mom0*mom), 1.0)
            call AddHistMP(hMP_ESet5, pPart, mom0, 1.0/(mom0*mom), 1.0)

            if (abs(y) < dy) then
               call AddHistMP(hMP_mTSet2, pPart, mT, 1.0/(2*dy*mT**2), 1.0)
               call AddHistMP(hMP_mTSet4, pPart, mT, 1.0/(2*dy*mT**2), 1.0)
               call AddHistMP(hMP_mTSet5, pPart, mT, 1.0/(2*dy*mT**2), 1.0)
            end if


         end do
      end do

      call WriteHistMP(hMP_ESet2, file='E_Set2_final.dat', &
           add=1e-20, mul=mulfak, iColumn=1)
      call WriteHistMP(hMP_ESet4, file='E_Set4_final.dat', &
           add=1e-20, mul=mulfak, iColumn=1)
      call WriteHistMP(hMP_ESet5, file='E_Set5_final.dat', &
           add=1e-20, mul=mulfak, iColumn=1)

      call WriteHistMP(hMP_mTSet2, file='mT_Set2_final.dat', &
           add=1e-20, mul=mulfak, iColumn=1)
      call WriteHistMP(hMP_mTSet4, file='mT_Set4_final.dat', &
           add=1e-20, mul=mulfak, iColumn=1)
      call WriteHistMP(hMP_mTSet5, file='mT_Set5_final.dat', &
           add=1e-20, mul=mulfak, iColumn=1)


    end subroutine Spectra


  end subroutine DoHeavyIonAnalysis


  !****************************************************************************
  !****s* HeavyIonAnalysis/DoHeavyIonAnalysisTime
  ! NAME
  ! subroutine DoHeavyIonAnalysisTime(realParts, time)
  !
  ! PURPOSE
  ! Do some analysis for every time step
  !
  ! INPUTS
  ! * type(particle), dimension(:,:) :: realParts --- real particle vector
  ! * real :: time --- actual time
  !****************************************************************************
  subroutine DoHeavyIonAnalysisTime(realParts, time)
    use particleDefinition
    use inputGeneral, only: time_max, delta_T, timeForOutput, timeSequence
    use output, only: intTochar
    use IdTable, only: EOV, NOP

    type(particle), dimension(:,:), intent(in) :: realParts
    real, intent(in) :: time

    integer, save    :: isut=0      ! number of subsequent runs
    integer, save    :: itime=0

    if (time==0.) itime=0

    if (time>= timeForOutput) then

      if (mod(itime*delta_T,timeSequence)<1E-4) then
        call writeParticleVectorHI
        call RapiditySpectra
      end if

      itime = itime + 1

    end if

    if (abs(time-time_max) < 1E-4) isut = isut + 1

  contains

    !**************************************************************************
    !****is* DoHeavyIonAnalysisTime/writeParticleVectorHI
    ! NAME
    ! subroutine writeParticleVectorHI()
    !
    ! PURPOSE
    ! write out the particle vectors in the format used for HI analysis
    !**************************************************************************
    subroutine writeParticleVectorHI
      use minkowski, only: abs4
      use inputGeneral, only: eventtype
      use eventtypes, only: HeavyIon, hadron
      use initHeavyIon, only: b_HI => b
      use initHadron, only: b_had => b

      integer :: i, j, fact, indFree
      real    :: stossParameter

      !************************************************************************
      !****o* HeavyIonAnalysis/DoHIATime___.dat
      ! NAME
      ! file DoHIATime___.dat
      ! PURPOSE
      ! Contains the full dump of the real particle vector at a certain time
      ! step.
      ! The time is given in the file name and in the first line of the file.
      ! The columns have the following meaning:
      ! * 1    = particle ID
      ! * 2    = charge
      ! * 3    = vac. mass in GeV
      ! * 4    = eff. mass in GeV
      ! * 5    = isFree (particle is bound or free?)
      ! * 6-8  = position (x,y,z) in fm
      ! * 9-11 = momentum (px,py,pz) in GeV
      ! * 12   = ensemble no.
      ! * 13   = run no.
      ! * 14   = impact parameter in fm
      !************************************************************************
      open(103,file='DoHIATime'//intToChar(nint(time/delta_T))//'.dat',position='Append')
      if (isut==0) then
         write(103,*) '# time = ',time,' fm/c'
         write(103,*)
      end if

      if (eventtype==HeavyIon) then
         stossParameter = b_HI
      else if (eventtype==hadron) then
         stossParameter = b_had
      else
         write(*,*) ' Problem in DoHeavyIonAnalysisTime: impact parameter is not defined'
         stossParameter=0.
      end if

      Loop_over_ensembles: do i = 1,size(realParts,dim=1)
         Loop_over_particles : do j = 1,size(realParts,dim=2)
            if (realParts(i,j)%id == NOP) cycle Loop_over_particles
            if (realParts(i,j)%id == EOV) exit Loop_over_particles
            if (realParts(i,j)%ID > 0) then
               if (realParts(i,j)%antiparticle) then
                  fact = -1
               else
                  fact = 1
               end if
               if (isFree(realParts,i,j,realParts(i,j))) then
                  indFree=1
               else
                  indFree=0
               end if
               write(103,50) fact*realParts(i,j)%ID, realParts(i,j)%charge, &
                             realParts(i,j)%mass, abs4(realParts(i,j)%momentum)**2, indFree, &
                             realParts(i,j)%position(1:3), realParts(i,j)%momentum(1:3), &
                             i, isut+1, stossParameter
            end if
         end do Loop_over_particles
      end do Loop_over_ensembles
      close(103)

50    format(2i4,2f9.4,i2,6f9.3,2i6,f6.2)

    end subroutine writeParticleVectorHI


    !**************************************************************************
    !****is* DoHeavyIonAnalysisTime/RapiditySpectra
    ! NAME
    ! subroutine RapiditySpectra()
    !
    ! PURPOSE
    ! generate and write out some rapidity spectra
    !**************************************************************************
    subroutine RapiditySpectra
      use idTable, only: nucleon
      use histf90

      integer, parameter :: nrap = 60
      real, parameter :: ystart=-6.0, dy=0.1

      type(histogram), allocatable, save :: dNdy(:)
      logical, save :: RapidityInitFLAG=.true.

      integer :: i, j, k, N, nEns
      real :: yb, p0, pz
      character(len=100) :: title

      if (RapidityInitFLAG) then
        N = int((time_max - timeForOutput)/timeSequence)
        ! print *,"RapiditySpectra: allocating histograms: ", N, time_max, timeForOutput, timeSequence
        allocate(dNdy(0:N))
        do k=0,N
          write(title,'(A,f7.2,A)') "nucleon rapidity distribution dN/dy at time = ", timeForOutput + k*timeSequence, " fm/c"
          call createHist (dNdy(k),  title, ystart, ystart+2*nrap*dy, dy)
        end do
        RapidityInitFLAG = .false.
      end if

      k = nint((time - timeForOutput) / timeSequence)

      nEns = size(realParts,dim=1)

      Loop_over_ensembles: do i = 1,nEns
         Loop_over_particles : do j = 1,size(realParts,dim=2)
            if (realParts(i,j)%id == NOP) cycle Loop_over_particles
            if (realParts(i,j)%id == EOV) exit Loop_over_particles

            if (realParts(i,j)%ID == nucleon) then

               p0 = realParts(i,j)%momentum(0)
               pz = realParts(i,j)%momentum(3)
               yb = 0.5*log( (p0+pz)/(p0-pz) )

               call addHist(dNdy(k), yb, 1.0)
            end if
         end do Loop_over_particles
      end do Loop_over_ensembles

      call writeHist(dNdy(k),  file='RapidityDistributions'//intTochar(nint(time/delta_T))//'.dat', mul = 1./nEns)

    end subroutine RapiditySpectra


  end subroutine DoHeavyIonAnalysisTime

  !****************************************************************************
  !****s* HeavyIonAnalysis/HeavyIon_evol
  ! NAME
  ! subroutine HeavyIon_evol(realParts, time, timestep)
  !
  ! PURPOSE
  ! time evolution of some global observables
  !
  ! INPUTS
  ! * type(particle), dimension(:,:) :: realParts --- real particle vector
  ! * real :: time --- actual time
  ! * integer :: timestep --- number of time step
  !
  ! NOTES
  ! * This routine is probably abused and its meaning is mixed with
  !   'DoHeavyIonAnalysisTime'. Thus it should be merged into this one and be
  !   deleted
  !****************************************************************************
  subroutine HeavyIon_evol(realParts, time, timestep)

    use IdTable
    use particleDefinition
    use densitymodule
    use RMF, only: getRMF_flag, fourMomDen_flag, g_omega, g_rho, g_sigma, &
         ModificationFactor
    use particleProperties, only: hadron
    use output, only: realTochar,intTochar
    use inputGeneral, only: delta_T, eventtype
    use constants, only: pi, hbarc, rhoNull, mN
    use coulomb, only: emfoca
    use PauliBlockingModule, only: pauliBlocking
    use initHadron, only: b,z,p_lab,E_bind
    use collisionNumbering, only: GetCountedEvents
    use eventtypes, only: ET_hadron => hadron
    use thermoDynamics, only: temperatureAt, muAt
    use histf90, only: histogram, createHist, addHist, writeHist

    type(particle), dimension(:,:), intent(in), target  :: realParts
    real, intent(in) :: time
    integer, intent(in) :: timestep

    type(particle), POINTER :: pPart

    real, dimension(-121:121,-2:2) :: parnum, parnum_free   ! particle numbers
    integer :: nEns,i,j,k,id,charge,Index1,Index2,Index3,npart
    integer, save :: icall=0     ! counter of calls
    integer, save :: isut=0      ! number of subsequent run
    real :: rhobar_points,rhobar_gauss,rhorad,rhorad_n,rhorad_p,r1,r2,r, &
         velrad1,velrad2, &
         rhoz,rhoz_bar,rhoz_antibar,rho_bar_max,rho_antibar_max,rholrf, &
         endens,rhobar_local,m_inv
    real :: rhoz_Bar_points,rhoz_AntiBar_points,rhoz_BoundBar_points, &
         rhoz_TargetBar_points,&
         &energy,fnorm
    real :: mstar, pf, sigma_rad, factor !,Ef
    real, dimension(1:3) :: p2_aver

    ! Check conservation laws:
    real, dimension(0:3) :: p ! total 4-momentum of all particles
    real :: baryon_number,charge_number,strangeness

    ! Properties of the bound system:
    real, parameter :: rho_min_bound=0.01*rhoNull ! minimal density
    real :: B_bound, N_bound, Z_bound ! baryon, neutron and charge numbers
    real :: rho_bound                ! mass averaged density
    real :: E_kinColl_bound          ! collective kinetic energy
    real, dimension(0:3) :: P_bound  ! total 4-momentum (w/o Coulomb contribution in P_bound(0))
    real :: E_Coul_bound             ! Coulomb energy

    real :: cpot
    real, dimension(1:3) :: place,impuls,velColl
    real, dimension(0:3) :: momentum

    logical :: flag_output, blockFlag

    real :: rhoz_n, rhoz_p, Ef_p, Ef_n, pf_p, pf_n, temp, mub
    real :: Ecoul   ! EcoulNorm
    real :: pion_ratio

    real, dimension(1:3)     :: number_of_baryons, Qzz_aver, r2_aver, r_aver_sq
    real, dimension(1:3)     :: number_of_neutrons,number_of_protons,np_ratio
    real, dimension(1:3,1:3) :: r_aver
    real, dimension(1:6) :: probability

    real :: mulfak

    type(histogram) :: histMnucleon
    integer :: nNucleons

    if (initFlag) call init

    if (icall==0) then
       open(32,file='evol.dat')
       write(32,*)'# Column No., quantity:'
       write(32,*)'# 1  time'
       write(32,*)'# 2  nucleon multiplicity'
       write(32,*)'# 3  delta'
       write(32,*)'# 4  Higher Nonstrange Baryonic Resonances'
       write(32,*)'# 5  pi'
       write(32,*)'# 6  K'
       write(32,*)'# 7  Kbar'
       write(32,*)'# 8  Lambda'
       write(32,*)'# 9  Lambda free'
       write(32,*)'# 10 Sigma'
       write(32,*)'# 11 Sigma free'
       write(32,*)'# 12 Xi'
       write(32,*)'# 13 Xi free'
       write(32,*)'# 14 XiStar'
       write(32,*)'# 15 XiStar free'
       write(32,*)'# 16 K^*'
       write(32,*)'# 17 Kbar^*'
       write(32,*)'# 18 Y^* (S=-1 only)'
       write(32,*)'# 19 LambdaBar'
       write(32,*)'# 20 SigmaBar'
       write(32,*)'# 21 Ybar^* (S=1 only)'
       write(32,*)'# 22 XiBar'
       write(32,*)'# 23 XiBar^*'
       write(32,*)'# 24 J/Psi'
       write(32,*)'# 25 D'
       write(32,*)'# 26 Dbar'
       write(32,*)'# 27 D^*'
       write(32,*)'# 28 Dbar^*'

       open(33,file='dens.dat',status='unknown')
       write(33,'(A)') '# time rhocen_baryons rhocen_antibaryons ptot '// &
            'baryon_number charge strangeness '// &
            'temperature(central) mu_B(central)'

       open(40,file='conserv.dat')
       write(40,*)'# time: energy:  px:  py:  pz:  baryon number:  '// &
            'charge: strangeness:'

       if (getRMF_flag() .and. fourMomDen_flag ) then
          open(43,file='BoundSystem.dat')
          write(43,*)'# Properties of the bound system:'
          write(43,*)'# time:   B:   N:   Z:  P(0:3), AGeV:   '// &
               'E_kinColl, AGeV:  E_Coul, AGeV:  <rho>, fm ^-3:'
       end if

       open(44,file='collisions.dat')
       write(44,*)'# time N_2body(integrated) N_2body(timestep)'

       open(50,file='isospinRatios.dat')

    end if

    icall = icall + 1
    if ( abs(time-delta_T) < 1.e-06 ) isut=isut+1
    flag_output=.false.

    parnum= 0.
    parnum_free= 0.
    rhobar_points= 0.
    endens= 0.
    p2_aver= 0.
    r2_aver=0.
    r_aver=0.0
    Qzz_aver=0.
    p(:)=0.
    baryon_number=0.
    charge_number=0.
    strangeness=0.
    number_of_baryons=0.
    number_of_neutrons=0.
    number_of_protons=0.
    E_Coul_bound=0.

    call createHist(histMnucleon, 'dN/dM nucleon', 0.8,0.95,0.002)
    nNucleons = 0

    nEns = size(realParts,dim=1)
    ensemble_loop : do i = 1,nEns
       particle_loop : do j = 1,size(realParts,dim=2)

          pPart => realParts(i,j)

          if (pPart%ID == EOV) exit particle_loop
          if (pPart%ID == NOP) cycle particle_loop

          if (NucleonMassPlot) then
             if (pPart%ID == 1 .and. .not.pPart%antiparticle) then
                call addHist(histMnucleon, sqrtS(pPart), 1.0)
                nNucleons = nNucleons+1
             end if
          end if

          factor = merge(-1., 1., pPart%antiparticle)
          id = pPart%ID*int(factor)
          charge = pPart%charge
          if (abs(id).le.121) then
             parnum(id,charge) = parnum(id,charge) + 1.
             if (id.eq.Lambda .or. id.eq.SigmaResonance .or. id.eq.Xi .or. id.eq.XiStar) then
                if (IsFree(realParts,i,j,pPart)) parnum_free(id,charge) = parnum_free(id,charge) + 1.
             end if
          end if
          p(:)=p(:)+pPart%momentum(:)
          if (isBaryon(pPart%ID)) then
             baryon_number=baryon_number+factor
             strangeness=strangeness+real(hadron(pPart%ID)%strangeness)*factor
          else if (isMeson(pPart%ID)) then
             strangeness=strangeness+real(hadron(pPart%ID)%strangeness)
          end if
          charge_number=charge_number+real(pPart%charge)

          if (all(abs(pPart%position) < 0.5)) then
             !       Central baryon density:
             if (isBaryon(pPart%ID)) then
                rhobar_points=  rhobar_points + factor
             end if
             !       Central energy density:
             endens= endens + pPart%momentum(0)
          end if

          if (isBaryon(pPart%ID)) then

             !       Average squares of momentum components :
             p2_aver(1:3)= p2_aver(1:3) + pPart%momentum(1:3)**2

             if (.not.pPart%antiparticle) then

                !         Indexes in large grid:
                Index1=NINT(pPart%position(1)/gridSpacing(1))
                Index2=NINT(pPart%position(2)/gridSpacing(2))
                Index3=NINT(pPart%position(3)/gridSpacing(3))

                if (        abs(Index1).le.gridPoints(1) &
                     & .and. abs(Index2).le.gridPoints(2) &
                     & .and. abs(Index3).le.gridPoints(3)  ) then

                   if (getRMF_flag()) then
                      factor=ModificationFactor(nucleon,.true.)    ! Modification factor of the antinucleon coupling constants
                      rhobar_local=(   densityField(Index1,Index2,Index3)%baryon(0)  &
                           & + factor*totalDensity(Index1,Index2,Index3) )/(1.+factor)  ! density of the baryons
                   else
                      ! w/o RMF "totalDensity" is actually (baryon-antibaryon) density and
                      ! densityField(0,0,0)%baryon(0) is (baryon+antibaryon) density (see density.f90: addTodensityfield)
                      rhobar_local=(   densityField(Index1,Index2,Index3)%baryon(0)  &
                           & + totalDensity(Index1,Index2,Index3) )/2.                 ! density of the baryons
                   end if

                   !cut on density (0.1*rho_sat)
                   if (rhobar_local.gt.0.1*rhoNull) then
                      number_of_baryons(1)=number_of_baryons(1)+1.
                      r2_aver(1)= r2_aver(1) &
                           & + dot_product(pPart%position,pPart%position)
                      r_aver(1,:)= r_aver(1,:) + pPart%position(:)
                      Qzz_aver(1)=Qzz_aver(1) + 2.*pPart%position(3)**2 &
                           &- pPart%position(1)**2 - pPart%position(2)**2
                      if (pPart%ID==1 .and. pPart%charge==0) &
                           & number_of_neutrons(1)=number_of_neutrons(1)+1.
                      if (pPart%ID==1 .and. pPart%charge==1) &
                           & number_of_protons(1)=number_of_protons(1)+1.
                   end if

                   !cut on density (0.01*rho_sat)
                   if (rhobar_local.gt.0.01*rhoNull) then
                      number_of_baryons(2)=number_of_baryons(2)+1.
                      r2_aver(2) = r2_aver(2) &
                           & + dot_product(pPart%position,pPart%position)
                      r_aver(2,:)= r_aver(2,:) + pPart%position(:)
                      Qzz_aver(2)=Qzz_aver(2) + 2.*pPart%position(3)**2 &
                           &- pPart%position(1)**2 - pPart%position(2)**2
                      if (pPart%ID==1 .and. pPart%charge==0) &
                           & number_of_neutrons(2)=number_of_neutrons(2)+1.
                      if (pPart%ID==1 .and. pPart%charge==1) &
                           & number_of_protons(2)=number_of_protons(2)+1.
                   end if

                   !cut on binding energy (E<0)
                   if (getRMF_flag()) then
                      call Particle4Momentum_RMF(pPart,momentum)
                      energy=momentum(0)
                   else
                      energy=pPart%momentum(0)
                   end if
                   place(1:3)= pPart%position(1:3)
                   impuls(1:3)= pPart%momentum(1:3)
                   energy = energy + emfoca(place,impuls,pPart%charge,pPart%ID)
                   if ( energy - pPart%mass < 0.) then
                      number_of_baryons(3)=number_of_baryons(3)+1.
                      r2_aver(3) = r2_aver(3) &
                           & + dot_product(pPart%position,pPart%position)
                      r_aver(3,:)= r_aver(3,:) + pPart%position(:)
                      if (pPart%ID==1 .and. pPart%charge==0) &
                           & number_of_neutrons(3)=number_of_neutrons(3)+1.
                      if (pPart%ID==1 .and. pPart%charge==1) &
                           & number_of_protons(3)=number_of_protons(3)+1.
                   end if

                   if (rhobar_local.gt.rho_min_bound) then
                      place(1:3)= pPart%position(1:3)
                      impuls(1:3)= pPart%momentum(1:3)
                      E_Coul_bound = E_Coul_bound + 0.5*emfoca(place,impuls,pPart%charge,pPart%ID) !see notes in evaluateTotal4Momentum_RMF
                   end if

                end if

             end if

          end if

       end do particle_loop
    end do ensemble_loop

    mulfak = 1.0/float(nEns)

    if (DensityPlot) then
       if (mod(timestep,5)==0) then
          call writeDensityPlane('density_YZ_'//intTochar(timestep)//'.dat',1)
          call writeDensityPlane('density_XZ_'//intTochar(timestep)//'.dat',2)
          call writeDensityPlane('density_XY_'//intTochar(timestep)//'.dat',3)
       end if
    end if

    if (NucleonMassPlot) then
       if (mod(timestep,5)==0) then
          if (nNucleons>0) then
             call WriteHist(histMnucleon, &
                  file='Mnucleon_'//intTochar(timestep)//'.dat', &
                  mul=1.0/nNucleons)
          end if
       end if
    end if

    if (do_Tmunu) then
       if (mod(timestep,5)==0) then
          call doTmunu(realParts,timestep)
       end if
    end if

    call countParts(realParts, time, "")


    parnum= parnum*mulfak
    parnum_free= parnum_free*mulfak
    rhobar_points= rhobar_points*mulfak/1.**3
    endens= endens*mulfak/1.**3
    p2_aver(:)= p2_aver(:)*mulfak
    r2_aver(:)= r2_aver(:)/number_of_baryons(:)
    do i=1,3
       r_aver(i,:)  = r_aver(i,:)/number_of_baryons(i)
       r_aver_sq(i) = sqrt(dot_product(r_aver(i,:),r_aver(i,:)))
    end do
    Qzz_aver(:)=Qzz_aver(:)/number_of_baryons(:)
    p(:)=p(:)*mulfak
    baryon_number=baryon_number*mulfak
    charge_number=charge_number*mulfak
    strangeness=strangeness*mulfak
    number_of_baryons(:)=number_of_baryons(:)*mulfak
    number_of_protons(:)=number_of_protons(:)*mulfak
    number_of_neutrons(:)=number_of_neutrons(:)*mulfak
    E_Coul_bound=E_Coul_bound*mulfak

    !n/p-ratio:
    do i=1,3
       if ( number_of_protons(i) .ne. 0.0 ) then
          np_ratio(i) = number_of_neutrons(i) / number_of_protons(i)
       else
          np_ratio(i) = 0.0
       end if
    end do

    !pion-ratio pi^{-}/pi^{+}:
    if ( parnum(pion,+1).ne.0.0 ) then
       pion_ratio = parnum(pion,-1)/parnum(pion,+1)
    else
       pion_ratio = 0.0
    end if

    write(32,5) time,sum(parnum(nucleon,:)),sum(parnum(delta,:)),&
         sum(parnum(P11_1440:F37_1950,:)),&
         sum(parnum(pion,:)),sum(parnum(kaon,:)),sum(parnum(kaonBar,:)),&
         sum(parnum(Lambda,:)),sum(parnum_free(Lambda,:)),&
         sum(parnum(SigmaResonance,:)),sum(parnum_free(SigmaResonance,:)),&
         sum(parnum(Xi,:)),sum(parnum_free(Xi,:)),&
         sum(parnum(XiStar,:)),sum(parnum_free(XiStar,:)),&
         sum(parnum(kaonStar,:)),sum(parnum(kaonStarBar,:)),&
         sum(parnum(Sigma_1385:Sigma_1915,:)),&
         sum(parnum(-Lambda,:)),sum(parnum(-SigmaResonance,:)),&
         sum(parnum(-Sigma_1915:-Sigma_1385,:)),&
         sum(parnum(-Xi,:)),sum(parnum(-XiStar,:)),&
         sum(parnum(JPsi,:)),sum(parnum(dMeson,:)),sum(parnum(dBar,:)),&
         sum(parnum(dStar,:)),sum(parnum(dStarBar,:))

    write(50,5) time,(np_ratio(i),i=1,3),pion_ratio

5   format(1P,50(1x,e13.6))

    rhobar_gauss= 0.
    do k = -2,2
       do j = -2,2
          do i = -2,2
             rhobar_gauss= rhobar_gauss + densityField(i,j,k)%baryon(0)
          end do
       end do
    end do
    rhobar_gauss= rhobar_gauss/5.**3

    if (getRMF_flag()) then
       factor=ModificationFactor(nucleon,.true.)    ! Modification factor of the antinucleon coupling constants
       rhoz_bar=(densityField(0,0,0)%baryon(0)+factor*totalDensity(0,0,0))/(1.+factor)  ! density of the baryons
       rhoz_antibar=(totalDensity(0,0,0)-densityField(0,0,0)%baryon(0))/(1.+factor)  ! density of the antibaryons
    else
       ! w/o RMF "totalDensity" is actually (baryon-antibaryon) density and
       ! densityField(0,0,0)%baryon(0) is (baryon+antibaryon) density (see density.f90: addTodensityfield)
       rhoz_bar=(densityField(0,0,0)%baryon(0)+totalDensity(0,0,0))/2.
       rhoz_antibar=(densityField(0,0,0)%baryon(0)-totalDensity(0,0,0))/2.
    end if

    temp = temperatureAt ((/0.,0.,0./))
    mub = muAt (rhoz_bar, temp)
    write(33,5) time, rhoz_bar, rhoz_antibar, p(:), baryon_number, charge_number, strangeness, temp, mub

    write(40,5) time,(r_aver_sq(i), sqrt(r2_aver(i)), number_of_baryons(i), i=1,3)

    if (getRMF_flag() .and. fourMomDen_flag ) then
       P_bound(:)=0.
       E_kinColl_bound=0.
       B_bound=0.
       N_bound=0.
       Z_bound=0.
       rho_bound=0.
       do Index1=-gridpoints(1),gridpoints(1)
          do Index2=-gridPoints(2),gridPoints(2)
             do Index3=-gridPoints(3),gridPoints(3)
                rhobar_local=densityField(Index1,Index2,Index3)%baryon(0)
                if (rhobar_local.gt.rho_min_bound) then
                   P_bound(0:3)=P_bound(0:3)+fourMomentumDensity(Index1,Index2,Index3,0:3)
                   m_inv=fourMomentumDensity(Index1,Index2,Index3,0)**2 &
                        &-dot_product(fourMomentumDensity(Index1,Index2,Index3,1:3),&
                        &fourMomentumDensity(Index1,Index2,Index3,1:3))
                   m_inv=sqrt(max(0.,m_inv))
                   E_kinColl_bound=E_kinColl_bound+fourMomentumDensity(Index1,Index2,Index3,0)&
                        &-m_inv
                   B_bound=B_bound+rhobar_local
                   N_bound=N_bound+densityField(Index1,Index2,Index3)%neutron(0)
                   Z_bound=Z_bound+densityField(Index1,Index2,Index3)%proton(0)
                   rho_bound=rho_bound+rhobar_local**2
                end if
             end do
          end do
       end do
       P_bound(0:3)=P_bound(0:3)*gridSpacing(1)*gridSpacing(2)*gridSpacing(3)
       E_kinColl_bound=E_kinColl_bound*gridSpacing(1)*gridSpacing(2)*gridSpacing(3)
       B_bound=B_bound*gridSpacing(1)*gridSpacing(2)*gridSpacing(3)
       N_bound=N_bound*gridSpacing(1)*gridSpacing(2)*gridSpacing(3)
       Z_bound=Z_bound*gridSpacing(1)*gridSpacing(2)*gridSpacing(3)
       rho_bound=rho_bound*gridSpacing(1)*gridSpacing(2)*gridSpacing(3)/B_bound

       write(43,5) time, B_bound, N_bound, Z_bound, P_bound(:)/B_bound, E_kinColl_bound/B_bound,&
            &E_Coul_bound/B_bound, rho_bound
    end if

    write(44,5) time, float(getCountedEvents(0,2,1))/float(nEns), float(getCountedEvents(1,2,1))/float(nEns)

    ! flush all files, so that output is actually written to disk
    flush(32)
    flush(33)
    flush(40)
    flush(44)
    flush(50)

    !-----------------------------------------------------------------------------------------------
    ! more detailed output, if flag_outputDetailed=.true.
    !-----------------------------------------------------------------------------------------------
    if ( .not. flag_outputDetailed ) return
    !-----------------------------------------------------------------------------------------------

    if ( abs(time-delta_T) < 1.e-06 .or. abs(time-nint(time)) < 0.5*delta_T ) then
       flag_output=.true.
       open(34,file='rhorad_'//realTochar(time+1.e-06)//'_'//intTochar(isut)//'.dat',position='append')
       write(34,'(A,3(f10.4,2x))')'# gridSpacing, fm:', gridSpacing
       write(34,'(A,f10.4)')'# time, fm/c: ', time
       write(34,'(A,i3)')'# run number: ', isut
       open(35,file='rhoz_'//realTochar(time+1.e-06)//'_'//intTochar(isut)//'.dat',position='append')
       write(35,'(A,3(f10.4,2x))')'# gridSpacing, fm:', gridSpacing
       write(35,'(A,f10.4)')'# time, fm/c: ', time
       write(35,'(A,i3)')'# run number: ', isut
       open(38,file='rhozx_'//realTochar(time+1.e-06)//'_'//intTochar(isut)//'.dat',position='append')
       write(38,'(A,3(f10.4,2x))')'# gridSpacing, fm:', gridSpacing
       write(38,'(A,f10.4)')'# time, fm/c: ', time
       write(38,'(A,i3)')'# run number: ', isut
       if (getRMF_flag()) then
          open(36,file='Fields_'//realTochar(time+1.e-06)//'_'//intTochar(isut)//'.dat',position='append')
          write(36,'(A,3(f10.4,2x))')'# gridSpacing, fm:', gridSpacing
          write(36,'(A,f10.4)')'# time, fm/c: ', time
          write(36,'(A,i3)')'# run number: ', isut
       end if
       open(37,file='pauli_'//realTochar(time+1.e-06)//'_'//intTochar(isut)//'.dat',position='append')
       write(37,'(A,f10.4)')'# time, fm/c: ', time
       write(37,'(A,i3)')'# run number: ', isut
       write(37,*)'# Momentum, GeV/c:    PauliBlock:'
    end if

    if (icall == 1) then
       open(39,file='dens_max.dat',status='unknown')
       if (eventtype == ET_hadron) &
            write(39,'(A,4(e15.8,1x))')" # hadron's b, z, p_lab, E_bind: ", b, z, p_lab, E_bind
       write(39,*)'# time  rho_bar_max  rho_antibar_max'
    end if


    Final_Output : if (flag_output) then

       write(34,*) '# r, fm:   rhorad:  rhorad_n:   rhorad_p:  velrad1/c:   velrad2/c:  PauliBlock(1-6):'
       do k= 0,min(gridpoints(1),gridpoints(2),gridpoints(3))

          r1=max(0.,(float(k)-0.5)*gridSpacing(1))
          r2=(float(k)+0.5)*gridSpacing(1)

          npart=0
          rhorad=0.
          rhorad_n=0.
          rhorad_p=0.
          velrad1=0.
          velrad2=0.
          sigma_rad=0.

          do i = 1,nEns
             particle_loop1 : do j = 1,size(realParts,dim=2)

                pPart => realParts(i,j)

                if (pPart%ID == EOV) exit particle_loop1
                if (pPart%ID == NOP) cycle particle_loop1

                place(1:3)= pPart%position(1:3)
                r=sqrt(dot_product(place(1:3),place(1:3)))

                if ( r.gt.r1 .and. r.le.r2 ) then

                   npart=npart+1

                   !             Indexes in large grid:
                   Index1=NINT(place(1)/gridSpacing(1))
                   Index2=NINT(place(2)/gridSpacing(2))
                   Index3=NINT(place(3)/gridSpacing(3))

                   if (        abs(Index1).le.gridPoints(1) &
                        & .and. abs(Index2).le.gridPoints(2) &
                        & .and. abs(Index3).le.gridPoints(3)  ) then

                      rhobar_local=densityField(Index1,Index2,Index3)%baryon(0)
                      rhorad=rhorad+rhobar_local

                      rhorad_n=rhorad_n+densityField(Index1,Index2,Index3)%neutron(0)
                      rhorad_p=rhorad_p+densityField(Index1,Index2,Index3)%proton(0)

                      if (allocated(sigmaField)) &
                           & sigma_rad=sigma_rad+sigmaField(Index1,Index2,Index3)

                      if (rhobar_local.gt.1.e-06) then
                         velColl(1:3)=densityField(Index1,Index2,Index3)%baryon(1:3) &
                              &/rhobar_local
                         if (r.gt.1.e-06) velrad1=velrad1+dot_product(velColl(1:3),place(1:3))/r
                      end if

                      if (allocated(fourMomentumDensity)) then
                         if (fourMomentumDensity(Index1,Index2,Index3,0).gt.1.e-06) then
                            velColl(1:3)=fourMomentumDensity(Index1,Index2,Index3,1:3) &
                                 &/fourMomentumDensity(Index1,Index2,Index3,0)
                            if (r.gt.1.e-06) velrad2=velrad2+dot_product(velColl(1:3),place(1:3))/r
                         end if
                      end if

                   end if

                end if


             end do particle_loop1
          end do

          if (npart.gt.0) then
             rhorad=rhorad/float(npart)
             rhorad_n=rhorad_n/float(npart)
             rhorad_p=rhorad_p/float(npart)
             sigma_rad=sigma_rad/float(npart)
             velrad1=velrad1/float(npart)
             velrad2=velrad2/float(npart)
          end if

          mstar = mN + g_sigma*sigma_rad

          !********** Pauli blocking check, radial dependence: ****************************
          do i=1,6
             momentum(1:3)=0.
             momentum(0)=sqrt(mN**2+dot_product(momentum(1:3),momentum(1:3)))
             place(1:3)=0.
             if (i.eq.1) place(1)=float(k)*gridSpacing(1)
             if (i.eq.2) place(1)=-float(k)*gridSpacing(1)
             if (i.eq.3) place(2)=float(k)*gridSpacing(2)
             if (i.eq.4) place(2)=-float(k)*gridSpacing(2)
             if (i.eq.5) place(3)=float(k)*gridSpacing(3)
             if (i.eq.6) place(3)=-float(k)*gridSpacing(3)
             blockFlag=pauliBlocking(momentum,place,1,realParts,probability(i))
          end do
          !********************************************************************

          write(34,5) float(k)*gridSpacing(1), rhorad, rhorad_n, rhorad_p, velrad1, velrad2,&
               &probability   !, mstar

       end do


       !********** Pauli blocking check, momentum dependence: *****************************
       place(1:3)=0.
       do k=1,201
          do i=1,6
             momentum(1:3)=0.
             if (i.eq.1) momentum(1)=0.00240*float(k-1)
             if (i.eq.2) momentum(1)=-0.00240*float(k-1)
             if (i.eq.3) momentum(2)=0.00240*float(k-1)
             if (i.eq.4) momentum(2)=-0.00240*float(k-1)
             if (i.eq.5) momentum(3)=0.00240*float(k-1)
             if (i.eq.6) momentum(3)=-0.00240*float(k-1)
             momentum(0)=sqrt(mN**2+dot_product(momentum(1:3),momentum(1:3)))
             blockFlag=pauliBlocking(momentum,place,1,realParts,probability(i))
          end do
          write(37,5) 0.00240*float(k-1),probability
       end do
       !***********************************************************************

       if (getRMF_flag()) then
          write(35,'(2A)') &
               & '# z, fm:     rhoz_bar:   rhoz_antibar:',&
               &'rhoz_Bar_points:  rhoz_AntiBar_points:  rhoz_BoundBar_points: rhoz_TargetBar_points:'
       else
          write(35,'(2A)') &
               & '# z, fm:     rhoz_bar:   rhoz_antibar:   rholrf:',&
               &'rhoz_Bar_points:  rhoz_AntiBar_points:  rhoz_BoundBar_points: rhoz_TargetBar_points:'
       end if

       Loop_over_zGrid1 : do k= -gridpoints(3),gridpoints(3)

          rhoz= densityField(0,0,k)%baryon(0)

          !proton & neutron densities
          rhoz_n=densityField(0,0,k)%neutron(0)
          rhoz_p=densityField(0,0,k)%proton(0)

          if (getRMF_flag()) then

             factor=ModificationFactor(nucleon,.true.)    ! Modification factor of the antinucleon coupling constants

             rhoz_bar=(rhoz+factor*totalDensity(0,0,k))/(1.+factor)  ! density of the baryons

             rhoz_antibar=(totalDensity(0,0,k)-rhoz)/(1.+factor)     ! density of the antibaryons

             pf=(1.5*pi**2*max(0.,rhoz_bar))**0.333333*hbarc

             pf_p=(3.*pi**2*max(0.,rhoz_p))**0.333333*hbarc
             pf_n=(3.*pi**2*max(0.,rhoz_n))**0.333333*hbarc


             mstar=mN+g_sigma*sigmaField(0,0,k)

             !           Ef=sqrt(pf**2+mstar**2)+g_omega*omegaField(0,0,k,0)

             Ef_p=sqrt(pf_p**2+mstar**2)+g_omega*omegaField(0,0,k,0)
             Ef_n=sqrt(pf_n**2+mstar**2)+g_omega*omegaField(0,0,k,0)

          else

             ! Remember: w/o RMF totalDensity is (baryon-antibaryon) density
             rhoz_bar=(rhoz+totalDensity(0,0,k))/2.         ! density of the baryons
             rhoz_antibar=(rhoz-totalDensity(0,0,k))/2.     ! density of the antibaryons

             rholrf= sqrt(max( densityField(0,0,k)%baryon(0)**2 &
                  -densityField(0,0,k)%baryon(1)**2 &
                  -densityField(0,0,k)%baryon(2)**2 &
                  -densityField(0,0,k)%baryon(3)**2, 0. ))

          end if

          rhoz_Bar_points=0.
          rhoz_AntiBar_points=0.
          rhoz_BoundBar_points=0.
          rhoz_TargetBar_points=0.
          !        ECoul = 0.0
          !        EcoulNorm = 0.0

          do i = 1,nEns
             particle_loop2 : do j = 1,size(realParts,dim=2)

                pPart => realParts(i,j)

                if (pPart%ID == EOV) exit particle_loop2
                if (pPart%ID == NOP) cycle particle_loop2

                if (pPart%ID >= pion) cycle particle_loop2   ! Count only (anti)baryons

                place(1:3)= pPart%position(1:3)

                if ( abs(place(1)) <= 0.5*gridSpacing(1) .and. &
                     &abs(place(2)) <= 0.5*gridSpacing(2) .and. &
                     &abs(place(3)-float(k)*gridSpacing(3)) <= 0.5*gridSpacing(3) ) then

                   if ( .not.pPart%antiparticle ) then
                      rhoz_Bar_points=rhoz_bar_points+1.
                      if (getRMF_flag()) then
                         call Particle4Momentum_RMF(pPart,momentum)
                         energy=momentum(0)
                      else
                         energy=pPart%momentum(0)
                      end if
                      impuls(1:3)= pPart%momentum(1:3)
                      cpot = emfoca(place,impuls,pPart%charge,pPart%ID)
                      energy=energy+cpot
                      !for testing E_F=E_F(r)
                      !                   if (pPart%charge==1) then
                      !                      Ecoul = Ecoul + cpot
                      !                      EcoulNorm = EcoulNorm + 1.
                      !                   end if
                      ! write(*,*)'Id, energy:', pPart%ID, energy
                      if (energy-pPart%mass < 0) rhoz_BoundBar_points=rhoz_BoundBar_points+1.
                      if (pPart%event(1)==1) &
                           & rhoz_TargetBar_points=rhoz_TargetBar_points+1.
                   else
                      rhoz_antibar_points=rhoz_antibar_points+1.
                   end if

                end if

             end do particle_loop2
          end do

          if (getRMF_flag()) then
             !Coulomb contribution to the Fermi-Energy:
             place=(/0.,0.,float(k)*gridSpacing(3)/)
             impuls=(/0.,0.,pf/)
             Ecoul = emfoca(place,impuls,1,1)
          end if

          fnorm=1./(float(nEns)*gridSpacing(1)*gridSpacing(2)*gridSpacing(3))
          rhoz_Bar_points=rhoz_Bar_points*fnorm
          rhoz_AntiBar_points=rhoz_AntiBar_points*fnorm
          rhoz_BoundBar_points=rhoz_BoundBar_points*fnorm
          rhoz_TargetBar_points=rhoz_TargetBar_points*fnorm

          if (getRMF_flag()) then

             write(35,5)  float(k)*gridSpacing(3), &
                  & rhoz_bar, rhoz_antibar, &
                  & rhoz_Bar_points,rhoz_AntiBar_points,&
                  & rhoz_BoundBar_points,rhoz_TargetBar_points

             write(36,5) float(k)*gridSpacing(3), &
                  & pf_p, pf_n, rhoz_p, rhoz_n, &
                  & -g_sigma*sigmaField(0,0,k), &
                  & g_omega*omegaField(0,0,k,0), &
                  & g_rho*rhoField(0,0,k,0),ECoul, &
                  & (Ef_p+Ecoul+g_rho*rhoField(0,0,k,0)-mN), &
                  & (Ef_n-g_rho*rhoField(0,0,k,0)-mN)
          else

             write(35,5)  float(k)*gridSpacing(3),rhoz_bar,rhoz_antibar,rholrf,rhoz_Bar_points,&
                  &rhoz_AntiBar_points,rhoz_BoundBar_points,rhoz_TargetBar_points

          end if

       end do Loop_over_zGrid1

    end if Final_Output

    if (flag_output) write(38,'(A)') '# z, fm:  x, fm:     rhoz_bar:   rhoz_antibar:'
    if (getRMF_flag()) factor=ModificationFactor(nucleon,.true.)
    rho_bar_max=-0.1
    rho_antibar_max=-0.1

    Loop_over_yGrid2 : do j= -gridpoints(2),gridpoints(2)
       Loop_over_zGrid2 : do k= -gridpoints(3),gridpoints(3)
          Loop_over_xGrid2 : do i= -gridpoints(1),gridpoints(1)

             rhoz= densityField(i,j,k)%baryon(0)

             if (getRMF_flag()) then
                rhoz_bar=(rhoz+factor*totalDensity(i,j,k))/(1.+factor)  ! density of the baryons
                rhoz_antibar=(totalDensity(i,j,k)-rhoz)/(1.+factor)     ! density of the antibaryons
             else
                rhoz_bar=(rhoz+totalDensity(i,j,k))/2.         ! density of the baryons
                rhoz_antibar=(rhoz-totalDensity(i,j,k))/2.     ! density of the antibaryons
             end if

             if (rhoz_bar.gt.rho_bar_max) rho_bar_max=rhoz_bar
             if (rhoz_antibar.gt.rho_antibar_max) rho_antibar_max=rhoz_antibar

             if (flag_output .and. j.eq.0) then
                write(38,5)  float(k)*gridSpacing(3),float(i)*gridSpacing(1),&
                     & rhoz_bar, rhoz_antibar
                if (i.eq.gridpoints(1))  write(38,*)
             end if

          end do Loop_over_xGrid2
       end do Loop_over_zGrid2
    end do Loop_over_yGrid2

    write(39,5) time, rho_bar_max, rho_antibar_max

  end subroutine HeavyIon_evol


  !*** Determine whether j-th particle of i-th parallel ensemble is free or not
  logical function IsFree(realParts,i,j,teilchen)

    use particleDefinition
    use densitymodule, only: Particle4Momentum_RMF
    use RMF, only: getRMF_flag
    use coulomb, only: emfoca

    type(particle), dimension(:,:), intent(in)  :: realParts
    integer, intent (in) :: i, j
    type(particle), intent(in)  ::  teilchen  ! particle of interest

    integer, parameter :: imode=1   ! 1 - criterion according to dstmin
    ! 2 - criterion according to binding energy
    real, parameter :: dstmin = 3.  ! minimum distance for free particle (fm)
    real :: tmp,dist2,energy
    integer :: j1
    real, dimension(0:3) :: momentum
    real, dimension(1:3) :: place

    place=teilchen%position

    if (imode.eq.1) then

       tmp=100.
       do j1=1,size(realParts,dim=2)
          if (j1.ne.j .or. teilchen%perturbative) then
             dist2=(realParts(i,j1)%position(1)-place(1))**2 &
                  &+(realParts(i,j1)%position(2)-place(2))**2 &
                  &+(realParts(i,j1)%position(3)-place(3))**2
             if (dist2.lt.tmp) tmp=dist2
          end if
       end do
       if (tmp.gt.dstmin**2) then
          IsFree=.true.
       else
          IsFree=.false.
       end if

    else if (imode.eq.2) then

       if (getRMF_flag()) then
          call Particle4Momentum_RMF(teilchen,momentum)
          energy=momentum(0)
       else
          energy=teilchen%momentum(0)
       end if
       energy=energy+emfoca(place,(/0.,0.,0./),realParts(i,j)%charge,realParts(i,j)%ID)
       if (energy.gt.teilchen%mass) then
          IsFree=.true.
       else
          IsFree=.false.
       end if

    end if

  end function IsFree


  subroutine countParts(Parts,time, prefix)

    use histMPf90, only: Map2HistMP, Map2HistMP_getN, WriteHistMP_Names
    use particleDefinition, only: particle

    type(particle),dimension(:,:),intent(in), target :: Parts
    real, intent(in) :: time
    character*(*),intent(in) :: prefix


    integer :: nHist, nEns, nPart, i, j, iID
    type(particle), POINTER :: pPart
    real :: mulfak

    if (.not.allocated(arrMultSet2)) then

       nHist = Map2HistMP_getN(2)
       allocate( arrMultSet2(0:nHist) )

       open(123,file="Mult_Set2.dat", status="unknown")
       call WriteHistMP_Names(2,123)
       close(123)
       open(123,file="final_Mult_Set2.dat", status="unknown")
       call WriteHistMP_Names(2,123)
       close(123)

       nHist = Map2HistMP_getN(4)
       allocate( arrMultSet4(0:nHist) )

       open(123,file="Mult_Set4.dat", status="unknown")
       call WriteHistMP_Names(4,123)
       close(123)
       open(123,file="final_Mult_Set4.dat", status="unknown")
       call WriteHistMP_Names(4,123)
       close(123)

       nHist = Map2HistMP_getN(5)
       allocate( arrMultSet5(0:nHist) )

       open(123,file="Mult_Set5.dat", status="unknown")
       call WriteHistMP_Names(5,123)
       close(123)
       open(123,file="final_Mult_Set5.dat", status="unknown")
       call WriteHistMP_Names(5,123)
       close(123)

    end if

    arrMultSet2 = 0.0
    arrMultSet4 = 0.0
    arrMultSet5 = 0.0

    nEns  = size(Parts,dim=1)
    nPart = size(Parts,dim=2)
    mulfak = 1.0/nEns

    do i=1,nEns
       do j=1,nPart
          pPart => Parts(i,j)
          if(pPart%Id <  0) exit
          if(pPart%Id == 0) cycle

          arrMultSet2(0) = arrMultSet2(0) + 1.0
          arrMultSet4(0) = arrMultSet4(0) + 1.0
          arrMultSet5(0) = arrMultSet5(0) + 1.0

          iID = Map2HistMP(pPart, 2)
          if (iID>0) then
             arrMultSet2(iID) = arrMultSet2(iID) + 1.0
          endif
          iID = Map2HistMP(pPart, 4)
          if (iID>0) then
             arrMultSet4(iID) = arrMultSet4(iID) + 1.0
          endif
          iID = Map2HistMP(pPart, 5)
          if (iID>0) then
             arrMultSet5(iID) = arrMultSet5(iID) + 1.0
          endif

       end do
    end do


    open(123,file=prefix//"Mult_Set2.dat",status="old",position='append')
    write(123,'(f11.5,1P,100E12.4,0P)') time, &
         & arrMultSet2(1:)*mulfak,arrMultSet2(0)*mulfak
    close(123)

    open(123,file=prefix//"Mult_Set4.dat",status="old",position='append')
    write(123,'(f11.5,1P,100E12.4,0P)') time, &
         & arrMultSet4(1:)*mulfak,arrMultSet4(0)*mulfak
    close(123)

    open(123,file=prefix//"Mult_Set5.dat",status="old",position='append')
    write(123,'(f11.5,1P,100E12.4,0P)') time, &
         & arrMultSet5(1:)*mulfak,arrMultSet5(0)*mulfak
    close(123)


  end subroutine countParts

  subroutine doTmunu(realParts,timestep)

    use IdTable
    use particleDefinition
    use TmunuDefinition
    use output, only: intTochar
    use rotation, only: rotateTo, rotateFrom
    use constants, only: pi
    use densityModule, only: boostToLRF
    use potentialModule, only: potential_LRF, trueEnergy

    type(particle), dimension(:,:), intent(in), target  :: realParts
    integer, intent(in) :: timestep

    logical, save :: initFlagTmunu=.true.

    integer :: iEns,iPart
    type(particle) :: part
    type(tTmunuNmu), dimension(:), allocatable, save :: ArrX,ArrY,ArrZ


    real, parameter :: dX = 0.2
    integer :: iBin
    integer, parameter :: nBin = 25
    real, save :: mulFak = 1.0

    type(tTmunuNmu), save :: tTmunuNmu0 ! used to reset the array !
    real :: w


    if (initFlagTmunu) then

       allocate(ArrX(nBin))
       allocate(ArrY(nBin))
       allocate(ArrZ(nBin))

       mulFak = 1.0/(size(realParts,dim=1)*dX**3) ! = 1/(nEns*V)
       initFlagTmunu = .false.
    end if

    ArrX = tTmunuNmu0
    ArrY = tTmunuNmu0
    ArrZ = tTmunuNmu0


    do iEns=1,size(realParts,dim=1)
       do iPart=1,size(realParts,dim=2)

          if (realParts(iEns,iPart)%ID == EOV) exit
          if (realParts(iEns,iPart)%ID == NOP) cycle

          part = realParts(iEns,iPart) ! create local copy
          w = 1.0

          select case (correctPot_Tmunu)

          case (1)
             call boostToLRF(part,1)  ! boost from calculation frame to LRF
             part%momentum(0) = part%momentum(0) - potential_LRF(part)
             call boostToLRF(part,2)  ! boost from LRF to calculation frame

          case (2)

             part%momentum(0) = trueEnergy(part, .true.)

          end select

          if (rotateZ_Tmunu) then

             part%momentum(1:3) = rotateFrom( part%position(1:3), &
                  part%momentum(1:3) )
             part%position(1:3) = rotateFrom( part%position(1:3), &
                  part%position(1:3) )

             w = w * dX**2/(4*pi* part%position(3)**2)
          end if

          if ( part%position(1)>0 &
               .and. abs(part%position(2))<dX/2 &
               .and. abs(part%position(3))<dX/2 ) then

             iBin = int( part%position(1) / dX )+1
             if (iBin <= nBin) call fillTmunu(ArrX(iBin), part)
          end if

          if ( part%position(2)>0 &
               .and. abs(part%position(1))<dX/2 &
               .and. abs(part%position(3))<dX/2 ) then

             iBin = int( part%position(2) / dX )+1
             if (iBin <= nBin) call fillTmunu(ArrY(iBin), part)
          end if

          if ( part%position(3)>0 &
               .and. abs(part%position(1))<dX/2 &
               .and. abs(part%position(2))<dX/2 ) then

             iBin = int( part%position(3) / dX )+1
             if (iBin <= nBin) call fillTmunu(ArrZ(iBin), part, w)
          end if

       end do
    end do

    open(123,file='Tmunu_'//intTochar(timestep)//'.dat', status="unknown")
    write(123,'(A)') headTmunu

    do iBin=1,nBin
       write(123,'(f11.4,1P,100E14.6,0P)') iBin*dX - dX/2, &
            & ArrX(iBin)%Tmunu(:)*mulfak, &
            & ArrX(iBin)%Nmu(:)*mulfak, &
            & ArrX(iBin)%Jmu(:)*mulfak, &
            & ArrX(iBin)%B*mulfak, ArrX(iBin)%S*mulfak
    end do
    write(123,*)
    write(123,*)

    do iBin=1,nBin
       write(123,'(f11.4,1P,100E14.6,0P)') iBin*dX - dX/2, &
            & ArrY(iBin)%Tmunu(:)*mulfak, &
            & ArrY(iBin)%Nmu(:)*mulfak, &
            & ArrY(iBin)%Jmu(:)*mulfak, &
            & ArrY(iBin)%B*mulfak, ArrY(iBin)%S*mulfak
    end do
    write(123,*)
    write(123,*)

    do iBin=1,nBin
       write(123,'(f11.4,1P,100E14.6,0P)') iBin*dX - dX/2, &
            & ArrZ(iBin)%Tmunu(:)*mulfak, &
            & ArrZ(iBin)%Nmu(:)*mulfak, &
            & ArrZ(iBin)%Jmu(:)*mulfak, &
            & ArrZ(iBin)%B*mulfak, ArrZ(iBin)%S*mulfak
    end do
    write(123,*)
    write(123,*)


    close(123)

  end subroutine doTmunu


end module HeavyIonAnalysis
