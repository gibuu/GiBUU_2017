!******************************************************************************
!****m* /AntibaryonWidth
! NAME
! module AntibaryonWidth
! PURPOSE
! Includes the routines administrating the simultaneous annihilation
! and calculation of annihilation and total widths of an antibaryon.
!******************************************************************************
module AntibaryonWidth

  implicit none
  private

  public :: DoAnnihilation

contains

  !****************************************************************************
  !****s* AntibaryonWidth/decideOnAnnihilation
  ! NAME
  ! subroutine decideOnAnnihilation(parts,time,flag,widthAverageOut,widthTotalAverageOut,
  ! baryonDensityAverageOut,sigmaAnniAverageOut,vRelAverageOut,
  ! srtS_vacuumAverageOut)
  ! PURPOSE
  ! Decide whether to annihilate an antiproton at a given time step.
  ! INPUTS
  ! * type(particle), dimension(:,:) :: parts --- real particle vector
  ! * real                           :: time  --- actual time step
  ! OUTPUT
  ! *  logical:: flag -- if .true. --- do annihilation
  ! *  real, optional :: widthAverageOut --- annihilation width of an antiproton
  !    averaged over all parallel ensembles
  ! *  real, optional :: widthTotalAverageOut --- total width of an antiproton
  !    averaged over all parallel ensembles
  ! *  real, optional :: baryonDensityAverageOut --- average density of baryons
  ! *  real, optional :: sigmaAnniAverageOut --- average annihilation cross section
  ! *  real, optional :: vRelAverageOut --- average relative velocity
  ! *  real, optional :: srtS_vacuumAverageOut --- average vacuum c.m. energy
  ! NOTES
  ! Presently suitable for the real particles and for the parallel ensemble
  ! mode.
  !****************************************************************************
  subroutine decideOnAnnihilation(parts,time,flag, &
       widthAverageOut,widthTotalAverageOut, &
       baryonDensityAverageOut,sigmaAnniAverageOut,vRelAverageOut, &
       srtS_vacuumAverageOut)

    use particleDefinition, only: particle
    use IdTable, only: isMeson
    use inputGeneral, only: delta_T
    use random, only: rn

    type(particle), intent(in), dimension(:,:) :: parts
    real, intent(in)                           :: time
    logical, intent(OUT):: flag
    real, optional, intent(OUT) :: widthAverageOut
    real, optional, intent(OUT) :: widthTotalAverageOut
    real, optional, intent(OUT) :: baryonDensityAverageOut
    real, optional, intent(OUT) :: sigmaAnniAverageOut
    real, optional, intent(OUT) :: vRelAverageOut
    real, optional, intent(OUT) :: srtS_vacuumAverageOut

    integer :: numberEnsembles, numberParticles
    integer :: ensemble, index, numAntiBaryons
    real :: widthAverage,widthTotalAverage,baryonDensityAverage,sigmaAnniAverage
    real :: vRelAverage,srtS_vacuumAverage
    real :: width,widthTotal,baryonDensity,sigmaAnni,vRel,srtS_vacuum,x,Pann

    numberEnsembles=size(parts,dim=1)
    numberParticles=size(parts,dim=2)

    widthAverage=0.
    widthTotalAverage=0.
    baryonDensityAverage=0.
    sigmaAnniAverage=0.
    vRelAverage=0.
    srtS_vacuumAverage=0.

    ensemble_loop : do ensemble=1,numberEnsembles

       numAntiBaryons=0

       index_loop :   do index=1,numberParticles

          if (parts(ensemble,index)%ID < 0) exit index_loop !!! ID==-1
          if (parts(ensemble,index)%ID <= 0) cycle index_loop

          if (.not.parts(ensemble,index)%antiparticle) cycle index_loop

          if (isMeson(parts(ensemble,index)%Id)) then
             write(*,*) ' Meson with antiparticle=.true. in decideOnAnnihilation !!!!'
             stop
          end if

          ! Antibaryon is selected, determine its width:
          numAntiBaryons=numAntiBaryons+1
          call  gamma(parts(ensemble,index),time,width,widthTotal,baryonDensity,sigmaAnni,&
               &vRel,srtS_vacuum)
          widthAverage= widthAverage + width
          widthTotalAverage= widthTotalAverage + widthTotal
          baryonDensityAverage= baryonDensityAverage + baryonDensity
          sigmaAnniAverage= sigmaAnniAverage + baryonDensity*sigmaAnni
          vRelAverage= vRelAverage + baryonDensity*vRel
          srtS_vacuumAverage= srtS_vacuumAverage + baryonDensity*srtS_vacuum

       end do index_loop

       if (numAntiBaryons.ne.1) then
          write(*,*) ' In decideOnAnnihilation: wrong number of antibaryons: ', numAntiBaryons
          write(*,*) ' in ensemble #: ', ensemble
          stop
       end if

    end do ensemble_loop

    widthAverage=widthAverage/float(numberEnsembles)
    widthTotalAverage=widthTotalAverage/float(numberEnsembles)
    baryonDensityAverage=baryonDensityAverage/float(numberEnsembles)
    if (baryonDensityAverage.gt.0.) then
       sigmaAnniAverage= sigmaAnniAverage/float(numberEnsembles)/baryonDensityAverage
       vRelAverage= vRelAverage/float(numberEnsembles)/baryonDensityAverage
       srtS_vacuumAverage=srtS_vacuumAverage/float(numberEnsembles)/baryonDensityAverage
    else
       sigmaAnniAverage=0.
       vRelAverage=0.
       srtS_vacuumAverage=0.
    end if

    ! Probability of annihilation:
    Pann= 1. - exp(-widthAverage*delta_T)

    x=rn()

    if ( x.lt.Pann ) then
       flag=.true.
    else
       flag=.false.
    end if

    if (present(widthAverageOut)) widthAverageOut=widthAverage
    if (present(widthTotalAverageOut)) widthTotalAverageOut=widthTotalAverage
    if (present(baryonDensityAverageOut)) baryonDensityAverageOut=baryonDensityAverage
    if (present(sigmaAnniAverageOut)) sigmaAnniAverageOut=sigmaAnniAverage
    if (present(vRelAverageOut)) vRelAverageOut=vRelAverage
    if (present(srtS_vacuumAverageOut)) srtS_vacuumAverageOut=srtS_vacuumAverage

  contains

    !**************************************************************************
    !****s* decideOnAnnihilation/gamma
    ! NAME
    ! subroutine gamma(antiBaryon,time,width,widthTotal,baryonDensityOut,sigmaAnniOut,vRelOut,srtS_vacuumOut)
    ! PURPOSE
    ! Determine the annihilation width of an antibaryon.
    ! INPUTS
    ! * type(particle) :: antiBaryon     ! antibaryon particle
    ! * real           :: time           ! actual time step
    ! OUTPUT
    ! * real :: width                   ! width w/r to annihilation (c/fm)
    ! * real :: widthTotal              ! total width (c/fm)
    ! * real, optional :: baryonDensityOut ! density of baryons (fm**-3)
    ! * real, optional :: sigmaAnniOut     ! annihilation cross section (mb)
    ! * real, optional :: vRelOut          ! relative velocity (c)
    ! * real, optional :: srtS_vacuumOut   ! c.m. energy (GeV) used in
    !   calculation of annihilation cross section
    !**************************************************************************
    subroutine gamma(antiBaryon,time,width,widthTotal,baryonDensityOut,sigmaAnniOut,vRelOut,srtS_vacuumOut)

      use constants, only: pi
      use XsectionRatios, only: accept_event
      use twoBodyTools, only: sqrtS_free, pCM, get_PInitial
      use particleDefinition, only: sqrtS, setToDefault
      use RMF, only: getRMF_flag
      use master_2Body, only: XsectionMaster, setKinematics, setKinematicsHiEnergy
      use densitymodule, only: getBaryonDensity
      use IdTable, only: pion
      use preEventDefinition
      use mediumDefinition
      use dichteDefinition
      use nucleusDefinition
      use nucleus, only: getTarget
      use densityStatic, only: staticDensity
      use Annihilation, only: annihilate
      use lorentzTrafo, only: lorentzCalcBeta
      use pauliBlockingModule, only: checkPauli

      type(particle), intent(in) :: antiBaryon
      real, intent(in)           :: time
      real, intent(out) :: width
      real, intent(out) :: widthTotal
      real, optional, intent(out) :: baryonDensityOut
      real, optional, intent(out) :: sigmaAnniOut
      real, optional, intent(out) :: vRelOut
      real, optional, intent(out) :: srtS_vacuumOut

      ! 1 --- use point particle density,
      ! 2 --- use smeared by gaussians density
      ! 3 --- use static density
      integer, parameter :: imode=2

      ! Radius (fm) of a sphere around antibaryon
      ! to look for annihilation partners.
      real, parameter ::    R=1.
      real, parameter ::    V=4./3.*pi*R**3  ! Volume of the sphere.
      integer, parameter :: baryonNumber_ens_max=10 ! Maximum number of the baryons
      ! found in the sphere in an ensemble
      integer, parameter :: nTrials=1   ! Number of the tried collisions to determine
      ! in-medium reduced cross sections.
      integer, parameter :: nFinal=20                     ! Size of the final state array
      type(particle), dimension(1:nFinal) :: finalState   ! Array of the final state particles

      type(medium) :: mediumAtColl
      type(particle) :: Baryon
      type(particle), dimension(:), allocatable :: BaryonsFound
      type(preEvent), dimension(1:4) :: chosenEvent
      type(dichte) :: density
      integer :: ensemble, index, nSuccess, n, nloop, HiEnergyType
      integer :: baryonNumber, maxBaryonNumber, i, i_max, maxId
      real :: baryonDensity, vRel
      real :: d2, srtS_vacuum, mstar_antibar, mstar_bar, srtS_star, srtS
      real, dimension(0:7) :: sigs
      real :: sigmaNonAnni
      real :: sigmaTotalRed, sigmaAnniRed, sigmaNonAnniRed
      real, dimension(0:3) :: momLRF
      real, dimension(1:3) :: betaToLRF,betaToCM
      logical :: collisionFlag
      logical :: HiEnergyFlag

      maxBaryonNumber=numberEnsembles*baryonNumber_ens_max
      !          maxBaryonNumber=baryonNumber_ens_max
      allocate(BaryonsFound(1:maxBaryonNumber))

      baryonNumber=0
      ensemble_loop : do ensemble=1,numberEnsembles

         index_loop :   do index=1,numberParticles

            if (parts(ensemble,index)%ID < 0) exit index_loop !!! ID==-1
            if (parts(ensemble,index)%ID <= 0) cycle index_loop

            if ( isMeson(parts(ensemble,index)%Id) .or.&
                 &parts(ensemble,index)%antiparticle ) cycle index_loop

            d2=dot_product( antiBaryon%position(1:3)-parts(ensemble,index)%position(1:3), &
                 & antiBaryon%position(1:3)-parts(ensemble,index)%position(1:3) )

            if (d2 .gt. R**2) cycle index_loop

            baryonNumber=baryonNumber+1

            if (baryonNumber.le.maxBaryonNumber) &
                 & BaryonsFound(baryonNumber)=parts(ensemble,index)


         end do index_loop

      end do ensemble_loop

      if (baryonNumber.eq.0) then
         deallocate(BaryonsFound)
         if (present(baryonDensityOut)) baryonDensityOut=0.
         if (present(VrelOut)) vRelOut=0.
         if (present(sigmaAnniOut)) sigmaAnniOut=0.
         if (present(srtS_vacuumOut)) srtS_vacuumOut=0.
         width=0.
         widthTotal=0.
         return
      end if

      select case (imode)
      case (1)
         baryonDensity=float(baryonNumber)/float(numberEnsembles)/V
         !baryonDensity=float(baryonNumber)/V
      case (2)
         baryonDensity = getBaryonDensity(antiBaryon%position)
      case (3)
         density=staticDensity(antiBaryon%position,getTarget())
         baryonDensity=density%baryon(0)
      case default
         write(*,*) ' In gamma (width of the antiproton): wrong imode '
         stop
      end select

      ! Random choice of the baryon from the set of found baryons:
      i_max=min(baryonNumber,maxBaryonNumber)
      x=rn()*float(i_max)
      i=nint(x)
      if (i.eq.0) i=i_max

      if (i.lt.1 .or. i.gt.maxBaryonNumber) then
         write(*,*) 'In width: index i out of bounds, i= ', i
         stop
      end if

      Baryon=BaryonsFound(i)  ! Baryon is selected

      betaToCM = lorentzCalcBeta (antiBaryon%momentum + Baryon%momentum, 'width')
      srtS=sqrtS((/antiBaryon,Baryon/),"width, srtS")

      if (.not.getRMF_flag() ) then
         srtS_vacuum=sqrtS_free((/antiBaryon,Baryon/))
      else   ! RMF mode:
         mstar_antibar = sqrtS(antiBaryon,'width, mstar_antibar')
         mstar_bar = sqrtS(Baryon,'width, mstar_bar')
         srtS_star = srtS
         srtS_vacuum = srtS_star - mstar_antibar - mstar_bar  + antiBaryon%mass + Baryon%mass
      end if

      mediumAtColl%useMedium=.true.
      mediumAtColl%density=baryonDensity

      ! Dummy settings:
      momLRF=0.
      betaToLRF=0.

      ! Simulate annnihilation:

      call XsectionMaster(srtS_vacuum,(/antiBaryon,Baryon/),mediumAtColl, &
           momLRF, chosenEvent, sigs, HiEnergyFlag)
      if (sigs(0) .lt. sigs(3)) then
         write(6,*)'srtS_vacuum, Total,Anni:',srtS_vacuum,sigs(0),sigs(3)
      end if

      sigmaNonAnni = sigs(0)-sigs(3)

      nSuccess=0
      loop_over_trials_anni : do n=1,nTrials
         call setToDefault(finalState)
         call annihilate(antiBaryon,Baryon,time,finalState,collisionFlag,HiEnergyType)
         if (.not.collisionFlag) cycle loop_over_trials_anni
         if ( .not.accept_event((/antiBaryon,Baryon/),finalState) ) then
            write(*,*) 'In width: event not accepted, annihilation'
            cycle loop_over_trials_anni
         end if
         nSuccess=nSuccess+1
      end do loop_over_trials_anni
      ! In-medium reduction:
      sigmaAnniRed=sigmaAnni*float(nSuccess)/float(nTrials)
      !write(*,*)'nTrials, nSuccess, sigmaAnniRed: ', nTrials, nSuccess, sigmaAnniRed

      ! Simulate scattering, CEX or inelastic production:
      nSuccess=0
      n=0
      loop_over_trials : do nloop=1,10 ! Upper limit must be set more than nTrials

         call setToDefault(finalState)

         call XsectionMaster(srtS_vacuum,(/antiBaryon,Baryon/),mediumAtColl, &
              momLRF, chosenEvent, sigs, HiEnergyFlag)
         if (sigs(0) .lt. sigs(3)) then
            write(6,*)'srtS_vacuum, Total,Anni:',srtS_vacuum,sigs(0),sigs(3)
         end if

         finalState(1:4)%ID=chosenEvent(1:4)%ID
         finalState(1:4)%charge=chosenEvent(1:4)%charge
         finalState(1:4)%antiParticle=chosenEvent(1:4)%antiParticle
         finalState(1:4)%mass=chosenEvent(1:4)%mass

         if (.not.HiEnergyFlag) then  !**** NOT HIGH ENERGY ****

            if ( finalState(1)%Id.eq.pion .and. finalState(2)%Id.eq.pion .and.&
                 & finalState(3)%Id.eq.pion ) cycle loop_over_trials ! Exclude annihilation

            n=n+1

            call ResetPosition((/antiBaryon,Baryon/),finalState,maxId)

            HiEnergyType=0

            call setKinematics(srtS,srtS_vacuum, betaToLRF,betaToCM,&
                 mediumAtColl,(/antiBaryon,Baryon/),&
                 finalState(1:maxID),collisionFlag)

         else  !**** HIGH ENERGY (presently w/o energy conservation) ****

            call setKinematicsHiEnergy(srtS,srtS_vacuum, sigs, betaToCM,&
                 & (/antiBaryon,Baryon/),&
                 & time, finalState, collisionFlag, HiEnergyType)

            if (HiEnergyType.eq.-3) cycle loop_over_trials ! Exclude annihilation

            n=n+1

            call ResetPosition((/antiBaryon,Baryon/),finalState,maxId)

         end if

         if (.not.collisionFlag) cycle loop_over_trials

         if ( .not.accept_event((/antiBaryon,Baryon/),finalState) ) then
            write(*,*) 'In width: event not accepted, nonannihilation'
            cycle loop_over_trials
         end if

         flag = checkPauli(finalState,parts)
         if (.not.flag) cycle loop_over_trials

         nSuccess=nSuccess+1

         if (n.eq.nTrials) exit loop_over_trials

      end do loop_over_trials

      ! In-medium reduction:
      if (n.eq.nTrials) then
         sigmaNonAnniRed=sigmaNonAnni*float(nSuccess)/float(nTrials)
      else
         sigmaNonAnniRed=0.
      end if
      !write(*,*)'n, nSuccess, sigmaNonAnniRed: ', n, nSuccess, sigmaNonAnniRed

      sigmaTotalRed=sigmaAnniRed+sigmaNonAnniRed

      ! Determine relative velocity:
      if ( .not.getRMF_flag() ) then
         vRel= get_pInitial((/antiBaryon,Baryon/),0) * sqrtS(antiBaryon,Baryon) &
              &/antiBaryon%momentum(0)/Baryon%momentum(0)
      else
         vRel= pcm(srtS_star,mstar_antibar,mstar_bar) * srtS_star &
              &/antiBaryon%momentum(0)/Baryon%momentum(0)
      end if

      width=baryonDensity*0.1*sigmaAnniRed*vRel
      widthTotal=baryonDensity*0.1*sigmaTotalRed*vRel

      deallocate(BaryonsFound)

      if (present(baryonDensityOut)) baryonDensityOut=baryonDensity
      if (present(VrelOut)) vRelOut=vRel
      if (present(sigmaAnniOut)) sigmaAnniOut=sigmaAnniRed
      if (present(srtS_vacuumOut)) srtS_vacuumOut=srtS_vacuum

    end subroutine gamma


    !**************************************************************************
    !****s* decideOnAnnihilation/ResetPosition
    ! NAME
    ! subroutine ResetPosition
    ! PURPOSE
    ! * calculate maxID, i.e. the number of finalState-particles
    ! * set position of finalState-particles:
    !   It loops over the FS particles vector and sets for the first
    !   FS-particle, which is baryon/meson as the first IS-particle,
    !   the spatial position as the first IS particle.
    !   Looping continues and the same is done for the next FS particle,
    !   which is of the same type as the secon IS particle.
    !   All other particles get as position the position of the interaction.
    !**************************************************************************
    subroutine ResetPosition(pair,finalState,maxId)

      use IdTable, only: isBaryon

      type(particle), dimension(1:2), intent(in) :: pair
      type(particle), dimension(:) :: finalState
      integer, intent(out) :: maxId

      integer :: k
      real, dimension(1:3) :: position
      logical :: flag1,flag2

      maxId= ubound(finalState,dim=1)
      do k=lbound(finalState,dim=1),ubound(finalState,dim=1)
         if (finalState(k)%Id.eq.0) then
            maxID=k-1
            exit
         end if
      end do

      position=0.5*(pair(1)%position+pair(2)%position)

      flag1= .false.
      flag2= .false.
      do k=lbound(finalState,dim=1),maxID
         if ( .not.flag1 .and. &
              isBaryon(finalState(k)%ID) .and. &
              (pair(1)%antiparticle.eqv.finalState(k)%antiparticle) ) then
            finalState(k)%position = pair(1)%position
            flag1= .true.
         else if ( .not.flag2 .and. &
              isBaryon(finalState(k)%ID) .and. &
              (pair(2)%antiparticle.eqv.finalState(k)%antiparticle) ) then
            finalState(k)%position = pair(2)%position
            flag2 = .true.
         else
            finalState(k)%position = position
         end if
      end do

    end subroutine ResetPosition


    !**************************************************************************
    !****s* decideOnAnnihilation/gamma_exact
    ! NAME
    ! subroutine gamma_exact(antiBaryon,time,width)
    ! PURPOSE
    ! Determine the annihilation width of an antibaryon.
    ! INPUTS
    ! * type(particle), intent(in) :: antiBaryon     ! antibaryon particle
    ! * real, intent(in)           :: time           ! actual time step
    ! OUTPUT
    ! * real, intent(out) :: width                   ! width w/r to annihilation (c/fm)
    !**************************************************************************
!     subroutine gamma_exact(antiBaryon,time,width)
!
!       use constants, only: pi
!       use XsectionRatios, only: accept_event
!       use twoBodyTools, only: sqrtS_free, pCM, get_PInitial
!       use particleDefinition, only: sqrtS
!       use RMF, only: getRMF_flag
!       use barAntiBar, only: sigmaBarAntiBar
!       use Annihilation
!       use mediumDefinition
!
!       type(particle), intent(in) :: antiBaryon     ! antibaryon particle
!       real, intent(in)           :: time           ! actual time step
!       real, intent(out) :: width                   ! width w/r to annihilation (c/fm)
!
!       real, parameter ::    R=1.        ! Radius (fm) of a sphere around antibaryon
!       ! to look for annihilation partners.
!       real, parameter ::    V=4./3.*pi*R**3  ! Volume of the sphere.
!       integer, parameter :: nTrials=1   ! Number of the tried annihilations to determine
!       ! in-medium reduction of the annihilation
!       ! cross section.
!       integer, parameter :: nFinal=20                     ! Size of the final state array
!
!       type(particle) :: Baryon
!       type(particle), dimension(1:nFinal) :: finalState   ! Array of the final state particles
!       type(medium) :: mediumAtColl
!       integer :: ensemble, index, nSuccess, n, HiEnergyType
!       real :: sum, d2, srtS_vacuum, mstar_antibar, mstar_bar, srtS_star, vRel
!       real :: sigmaTotal, sigmaElast, sigmaAnni
!       logical :: collisionFlag
!
!       sum=0.
!       ensemble_loop : do ensemble=1,numberEnsembles
!
!          index_loop :   do index=1,numberParticles
!
!             if (parts(ensemble,index)%ID < 0) exit index_loop !!! ID==-1
!             if (parts(ensemble,index)%ID <= 0) cycle index_loop
!
!             if ( isMeson(parts(ensemble,index)%Id) .or.&
!                  &parts(ensemble,index)%antiparticle ) cycle index_loop
!
!             d2=dot_product( antiBaryon%position(1:3)-parts(ensemble,index)%position(1:3), &
!                  & antiBaryon%position(1:3)-parts(ensemble,index)%position(1:3) )
!
!             if (d2 .gt. R**2) cycle index_loop
!
!             Baryon=parts(ensemble,index)   ! Baryon is selected
!
!             if (.not.getRMF_flag() ) then
!                srtS_vacuum=sqrtS_free((/antiBaryon,Baryon/))
!             else   ! RMF mode:
!                mstar_antibar = sqrtS(antiBaryon,'width, mstar_antibar')
!                mstar_bar = sqrtS(Baryon,'width, mstar_bar')
!                srtS_star = sqrtS((/antiBaryon,Baryon/),"width, srtS_star")
!                srtS_vacuum = srtS_star - mstar_antibar - mstar_bar  + antiBaryon%mass + Baryon%mass
!             end if
!
!             mediumAtColl%useMedium=.false.  ! Dummy setting
!
!             ! Determine annihilation cross section in vacuum:
!             call sigmaBarAntiBar(srtS_vacuum,(/antiBaryon,Baryon/),mediumAtColl,sigTotal=sigmaTotal,&
!                  &sigElastic=sigmaElast,sigAnnihilation=sigmaAnni)
!
!             nSuccess=0
!             loop_over_trials : do n=1,nTrials
!
!                ! Annihilate the antibaryon and baryon into some final state particles:
!                call annihilate(antiBaryon,Baryon,time,finalState,collisionFlag,HiEnergyType)
!
!                if (.not.collisionFlag) cycle loop_over_trials
!
!                if ( .not.accept_event((/antiBaryon,Baryon/),finalState) ) then
!                   write(*,*) 'In width: event not accepted'
!                   cycle loop_over_trials
!                end if
!
!                nSuccess=nSuccess+1
!
!             end do loop_over_trials
!
!             ! In-medium reduction of the annihilation cross section:
!             sigmaAnni=sigmaAnni*float(nSuccess)/float(nTrials)
!
!             if (sigmaAnni.gt.0.) then
!                ! Determine relative velocity:
!                if ( .not.getRMF_flag() ) then
!                   vRel= get_pInitial((/antiBaryon,Baryon/),0) * sqrtS(antiBaryon,Baryon) &
!                        &/antiBaryon%momentum(0)/Baryon%momentum(0)
!                else
!                   vRel= pcm(srtS_star,mstar_antibar,mstar_bar) * srtS_star &
!                        &/antiBaryon%momentum(0)/Baryon%momentum(0)
!                end if
!                sum= sum + vRel*sigmaAnni
!             end if
!
!          end do index_loop
!
!       end do ensemble_loop
!
!       width= 0.1*sum/float(numberEnsembles)/V
!
!     end subroutine gamma_exact


  end subroutine decideOnAnnihilation



  !****************************************************************************
  !****s* AntibaryonWidth/DoAnnihilation
  ! NAME
  ! subroutine DoAnnihilation(parts,time)
  ! PURPOSE
  ! Looks for an antibaryon and annihilates it with the closest baryon.
  ! INPUTS
  ! * type(particle), dimension(:,:) :: parts --- real particle vector
  ! * real :: time --- actual time step
  ! * real :: annihilationTime --- time when annihilation is enforced
  ! OUTPUT
  ! * type(particle), dimension(:,:) :: parts --- real particle vector changed
  ! NOTES
  ! Presently suitable for the real particles and for the parallel ensemble
  ! mode.
  !****************************************************************************
  subroutine DoAnnihilation(parts,time,annihilationTime)

    use particleDefinition, only: particle, AcceptGuessedNumbers, sqrtS
    use IdTable, only: isMeson, isBaryon
    use collisionNumbering
    use output, only: WriteParticleVector, DoPr
    use RMF, only: getRMF_flag
    use densitymodule, only: updateRMF, storeFields
    use Insertion, only: particlePropagated, setIntoVector
    use XsectionRatios, only: accept_event
    use Annihilation
    use inputGeneral, only: delta_T

    type(particle), intent(inOUT), target, dimension(:,:) :: parts
    real, intent(in)                                      :: time
    real, intent(in)                                      :: annihilationTime

    type(particle), pointer :: antibar, bar

    integer, parameter :: nFinal=20                     ! Size of the final state array
    type(particle), dimension(1:nFinal) :: finalState   ! Array of the final state particles

    integer :: numberEnsembles, numberParticles
    integer :: ensemble, index, index1, index_min

    real :: d2_min, d2, width, widthTotal, baryonDensity, sigmaAnni, vRel, srtS_vacuum
    integer :: n, i, number
    integer :: HiEnergyType
    logical :: numbersAlreadySet, setFlag, collisionFlag, flag1

    logical, parameter :: flagStochastic=.false.       ! if .true. --- decide by Monte Carlo,
    ! whether annihilation will happen according
    ! to the antibaryon annihilation width
    logical, save :: flagAnni=.false.
    real, save :: timePrevious=1000.
    real, save :: widthIntegral=0.
    real :: Psurv

    write(*,*) 'In DoAnnihilation, time: ', time

    if (.not.flagAnni .or. time-timePrevious.lt.1.e-06) then
       call decideOnAnnihilation(parts,time,flagAnni,width,widthTotal,baryonDensity,&
            &sigmaAnni,vRel,srtS_vacuum)
       widthIntegral=widthIntegral+width*delta_T
       Psurv=exp(-widthIntegral)
       if (Psurv.le.1.e-99) Psurv=0.
       if (time-timePrevious.lt.1.e-06) then
          open(1,file='antiBaryonWidth.dat',status='unknown')
       else
          open(1,file='antiBaryonWidth.dat',position='Append')
       end if
       write(1,'(1x,8(e13.6,1x))') time, Psurv, width, widthTotal, &
            & baryonDensity, sigmaAnni, vRel, srtS_vacuum
       timePrevious=time
    end if

    if (.not.flagStochastic) then
       if (time.gt.annihilationTime-1.e-06) then
          flagAnni=.true.
       else
          flagAnni=.false.
       end if
    end if

    if (.not.flagAnni) return

    flag1=.false.

    numberEnsembles=size(parts,dim=1)
    numberParticles=size(parts,dim=2)

    ensemble_loop : do ensemble=1,numberEnsembles

       index_loop :   do index=1,numberParticles

          if (parts(ensemble,index)%ID < 0) exit index_loop !!! ID==-1
          if (parts(ensemble,index)%ID <= 0) cycle index_loop

          if (.not.parts(ensemble,index)%antiparticle) cycle index_loop

          if (.not.isBaryon(parts(ensemble,index)%Id)) then
             write(*,*) ' Meson with antiparticle=.true. in DoAnnihilation !!!!'
             stop
          end if

          ! Antibaryon is selected

          antibar=>parts(ensemble,index)

          d2_min=10000.
          index_min=0

          index1_loop :   do index1=1,numberParticles

             if (parts(ensemble,index1)%ID < 0) exit index1_loop !!! ID==-1
             if (parts(ensemble,index1)%ID <= 0) cycle index1_loop

             if (      parts(ensemble,index1)%antiparticle &
                  & .or. isMeson(parts(ensemble,index1)%Id) ) cycle index1_loop

             d2=dot_product( antibar%position(1:3)-parts(ensemble,index1)%position(1:3), &
                  & antibar%position(1:3)-parts(ensemble,index1)%position(1:3) )

             if (d2 <= d2_min) then
                d2_min=d2
                index_min=index1
             end if

          end do index1_loop

          if (index_min.ne.0) then
             ! Baryon is selected
             bar=>parts(ensemble,index_min)
          else
             write(*,*) ' No baryon found in in DoAnnihilation !!!!'
             stop
          end if

          n=0
          Loop_over_annihilation_events : do

             n=n+1

             if (n.gt.100) then
                write(*,*) ' In DoAnnihilation: no annihilation after ',n-1,' trials'
                write(*,*) ' antibar%Id, antibar%charge: ', antibar%Id, antibar%charge
                write(*,*) ' bar%Id, bar%charge: ', bar%Id, bar%charge
                write(*,*) ' srtS_star: ', sqrtS((/antibar,bar/),"DoAnnihilation, srtS_star")
                exit Loop_over_annihilation_events
             end if

             ! Annihilate the antibaryon and baryon into some final state particles:
             call annihilate(antibar,bar,time,finalState,collisionFlag,HiEnergyType)

             if (.not.collisionFlag) cycle Loop_over_annihilation_events

             if ( .not.accept_event((/antibar,bar/),finalState) ) then
                write(*,*) ' Event not accepted'
                collisionFlag=.false.
                cycle Loop_over_annihilation_events
             end if

             exit Loop_over_annihilation_events

          end do Loop_over_annihilation_events


          if (collisionFlag) then
             flag1=.true.
          else
             cycle index_loop
          end if

          ! Label event by eventNumber:
          number=real_numbering()
          finalState%event(1)=number
          finalState%event(2)=number

          NumbersAlreadySet = AcceptGuessedNumbers()

          do i=1,nFinal
             if (finalState(i)%id <= 0) exit
             if (particlePropagated(finalState(i))) then
                ! Find empty space in the particle vector:
                call setIntoVector( finalState(i:i), &
                     & parts(ensemble:ensemble,:), setFlag, NumbersAlreadySet )
                ! Check that setting into real particle vector worked out :
                if (.not.setFlag) then
                   write(*,*) 'Real particle vector too small!'
                   write(*,*) size(finalState), lbound(parts), ubound(parts)
                   write(*,*) 'Dumping real particles to files "RealParticles_stop_*!"'
                   call WriteParticleVector("RealParticles_stop",parts)
                   stop ' annihilation'
                end if
             end if
          end do

          ! Destroy the baryon and antibaryon test particles:
          antibar%Id=0
          bar%Id=0

       end do index_loop

       if ( mod(ensemble,100).eq.0 .and. flag1 ) then
          if ( getRMF_flag() ) then
             if (DoPr(2)) write(*,*) 'Updating Relativistic Mean Fields after annihilation'
             call updateRMF(parts)
          end if
          flag1=.false.
       end if

    end do ensemble_loop

    if (getRMF_flag()) call storeFields

    call analyse

  contains

    subroutine analyse

      use IdTable

      integer, parameter :: Npion_max=20             ! Maximum number of pions produced in an event
      real, save :: sigma_Npion(0:Npion_max)

      integer :: Npion
      real :: Npion_aver

      integer, save :: ncalls=0

      ncalls=ncalls+1

      if (ncalls.eq.1) then
         open(2,file='annihilation_analyse.dat',status='unknown')
      end if

      sigma_Npion=0.
      Npion_aver=0.

      do ensemble=1,numberEnsembles

         Npion=0

         index_loop :   do index=1,numberParticles

            if (parts(ensemble,index)%ID < 0) exit index_loop !!! ID==-1
            if (parts(ensemble,index)%ID <= 0) cycle index_loop

            if ( parts(ensemble,index)%ID.eq.pion ) Npion=Npion+1

         end do index_loop

         Npion_aver=Npion_aver+float(Npion)
         if ( Npion.le.Npion_max ) sigma_Npion(Npion)=sigma_Npion(Npion)+1.

      end do

      write(2,*)'# time: ', time
      Npion_aver=Npion_aver/float(numberEnsembles)
      write(2,*)'# Average number of pions: ',  Npion_aver
      write(2,*)'# Number of ensembles: ', numberEnsembles
      write(2,*)'# N_pion:    N_events:'
      do Npion=0,Npion_max
         write(2,*) Npion, sigma_Npion(Npion)
      end do

    end subroutine analyse

  end subroutine DoAnnihilation


end module AntibaryonWidth
