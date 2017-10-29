!******************************************************************************
!****m* /resonanceCrossSections
! NAME
! module resonanceCrossSections
! PURPOSE
! Includes all routines to evaluate "baryon meson -> R -> X" cross sections
! with resonances R in the intermediate state.
!******************************************************************************
module resonanceCrossSections
  implicit none
  private

  public :: barMes2resonance, barMes_R_barMes, sigma_npi_n2pi_resonances

  !****************************************************************************
  !****g* resonanceCrossSections/fullPropagator
  ! PURPOSE
  ! Includes also the real parts in the resonance propagator. In former works
  ! (i.e. in the old Efffenberger code) this has been neglected.
  ! It should be set to .true. only if mediumSwitch_coll=.true. in the
  ! namelist width_Baryon.
  !
  ! SOURCE
  !
  logical, save :: fullPropagator=.false.
  !****************************************************************************

  logical, save :: initInput_Flag = .true.

  logical, parameter :: debug = .false.

contains


  subroutine readInput
    use output, only:Write_ReadingInput

    integer :: ios

    !**************************************************************************
    !****n* resonanceCrossSections/ResonanceCrossSections
    ! NAME
    ! NAMELIST ResonanceCrossSections
    ! PURPOSE
    ! Includes parameters:
    ! * fullPropagator
    !**************************************************************************
    NAMELIST /ResonanceCrossSections/ fullPropagator

    call Write_ReadingInput("ResonanceCrossSections",0)
    rewind(5)
    read(5,nml=ResonanceCrossSections,IOSTAT=ios)
    call Write_ReadingInput("ResonanceCrossSections",0,ios)
    write(*,*) '  Use the full self energy in the propagator of the resonances?',fullPropagator
    write(*,'(A)') '   (It should be TRUE only if mediumSwitch_coll=.true. in the namelist width_Baryon)'

    call Write_ReadingInput("ResonanceCrossSections",1)

    initInput_Flag = .false.

  end subroutine readInput


  !****************************************************************************
  !****f* resonanceCrossSections/resonanceMass
  ! NAME
  ! real function resonanceMass (ID, charge, media, FourMomentum, position, perturbative)
  ! PURPOSE
  ! Determines the bare (vacuum) mass of a resonance,
  ! if its 4-momentum (effective mass) is known.
  ! INPUTS
  ! * integer, intent(in) :: ID               ! ID of resonance
  ! * integer, intent(in) :: charge           ! charge of resonance
  ! * type(medium) :: media                   ! Medium at the point of the resonance
  ! * real, intent(in) :: FourMomentum(0:3)   ! 4-momentum of resonance in the LRF
  ! * real, intent(in) :: position(1:3)       ! position of resonance in calculation frame
  ! RESULT
  ! * resonance mass in GeV
  ! NOTES
  ! This was known as "detmass" in the old code.
  !****************************************************************************
  real function resonanceMass (ID, charge, media, FourMomentum, position, perturbative)
    use particleDefinition
    use mediumDefinition
    use baryonPotentialModule, only: baryonPotential
    use coulomb, only: emfoca

    integer, intent(in) :: ID, charge
    type(medium) :: media
    real, intent(in) :: FourMomentum(0:3), position(1:3)
    logical, intent(in) :: perturbative

    type(particle) :: resonance
    real :: potLRF

    ! determine hadronic potential
    if (media%useMedium) then
      resonance%position = position
      resonance%perturbative = perturbative
      resonance%charge = charge
      resonance%Id = ID
      resonance%momentum(0:3) = FOURmomentum(0:3)
      potLRF = BaryonPotential (resonance, media, .false.)
    else
      potLRF = 0.
    end if

    ! add Coulomb potential
    potLRF = potLRF + emfoca (position, FourMomentum(1:3), charge, ID)

    ! evaluate the 'bare' mass of the resonance in the LRF
    resonanceMass = sqrt( (FourMomentum(0)-potLRF)**2 - sum(FourMomentum(1:3)**2) )

  end function resonanceMass


  !****************************************************************************
  !****f* resonanceCrossSections/barMes2resonance
  ! NAME
  ! function barMes2resonance (idMeson_ini, idBaryon_ini, chargeMeson_ini, chargeBaryon_ini, propagated,
  ! mediumAtCollision, momentumLRF, masses, mesonMass, baryonMass, position, perturbative, srts) result (sigma)
  ! PURPOSE
  ! Evaluates baryon meson -> Resonance cross sections by summing over the resonances final state.
  ! One gets baryon meson -> X by summing over the (baryon meson -> R-> X) channels for all resonances.
  ! See eq. (2.52) in Effenberger Phd.
  ! INPUTS
  ! * integer, intent(in) :: idMeson_ini,idBaryon_ini
  ! * integer, intent(in) :: chargeMeson_ini,chargeBaryon_ini
  ! * logical, intent(in) :: propagated
  !   if .true., then onsider only propagated particles in final state, otherwise
  !   consider all particles in final state
  ! * type(medium), intent(in) :: mediumAtCollision      ! Medium information
  ! * real, intent(in),dimension(0:3) :: momentumLRF     ! Total  Momentum in LRF=Momentum of resonance in LRF
  ! * real, intent(in)     ::  mesonMass, baryonMass     ! masses of colliding particles.
  ! * real, intent(in), dimension(1:3) :: position       ! Position where resonance is produced
  ! * logical, intent(in)              :: perturbative   ! flag for the perturbative nature of the resonance
  ! * real, intent(in), optional          :: srts        ! sqrt(s) of the collision -- needed in RMF mode
  ! RESULT
  ! * real, intent(out),  dimension(Delta:nbar) :: sigma   --- Cross section for each intermediate baryon resonance in units of mB. The index denotes the resonance.
  ! * real, intent(out),  dimension(Delta:nbar) :: masses  --- Masses of the produced resonances; index denotes the resonance.
  !****************************************************************************
  function barMes2resonance (idMeson_ini, idBaryon_ini, chargeMeson_ini, chargeBaryon_ini, propagated, mediumAtColl, &
                             momLRF, masses, mesonMass, baryonMass, position, perturbative, srts) result (sigma)

    use mediumDefinition, only: medium, vacuum
    use baryonWidthMedium, only: WidthBaryonMedium, partialWidthBaryonMedium
    use particleProperties, only: hadron, isNonExotic
    use ClebschGordan, only: CG
    use twoBodyTools, only: pCM_sqr
    use constants, only: pi, GeVSquared_times_mb
    use IdTable, only: Delta, nbar, isMeson, isBaryon
    use RMF, only: getRMF_flag
    use selfenergy_baryons, only: get_realPart,selfEnergy_imag
    use minkowski, only: abs3, abs4

    integer, intent(in)              :: idMeson_ini, idBaryon_ini, chargeMeson_ini, chargeBaryon_ini
    logical, intent(in)              :: propagated
    type(medium), intent(in)         :: mediumAtColl
    real, intent(in), dimension(0:3) :: momLRF
    real, dimension(Delta:nbar), intent(out) :: masses
    real, intent(in)                 :: mesonMass, baryonMass, position(1:3)
    logical, intent(in)              :: perturbative
    real, intent(in), optional       :: srts
    real, dimension(Delta:nbar)      :: sigma

    integer :: resID, chargeRes, i1, i2, i3, iz1, iz2, iz3
    real :: momCM2, gamma_In, gammaTot, spinFactor, isoFactor, absMom_LRF, fullMass, PI_REAL, PI_IMAG

    if (initInput_Flag) call readInput()

    ! Check Input
    if (.not.(isMeson(idMeson_ini).and.isBaryon(idBaryon_ini))) then
       write(*,*) 'Error in in input of "resonanceXsection" '
       write(*,*) idMeson_ini,idBaryon_ini
       stop
    end if

    sigma=0.
    masses=0.

    absMom_LRF = abs3(momLRF)
    if (fullPropagator) fullmass = abs4(momLRF)

    ! Loop over intermediate resonances R : m B->  R -> X
    do resId=Delta,nbar
       if (.not.hadron(resId)%usedForXsections) cycle ! Exclude resonances
       if (propagated.and.(.not.(hadron(resId)%propagated))) cycle

       if (abs(hadron(idMeson_ini)%strangeness + hadron(idBaryon_ini)%strangeness - hadron(resID)%strangeness) > 0) cycle   ! strangeness conservation
       if (abs(hadron(idMeson_ini)%charm       + hadron(idBaryon_ini)%charm       - hadron(resID)%charm)       > 0) cycle   ! charm conservation

       ! Convert charges to isospins in SU(4) & evaluate isospin factors for initial state
       ! Q=Y/2+I_3 with Y=S+C+B  =>  I_3=Q-(S+B+C)/2.
       i1  = hadron(idMeson_ini)%isoSpinTimes2
       iz1 = nint(2.*chargeMeson_ini) - hadron(idmeson_ini)%strangeness - hadron(idmeson_ini)%charm
       i2  = hadron(idBaryon_ini)%isoSpinTimes2
       iz2 = nint(2.*chargeBaryon_ini-1.) - hadron(idbaryon_ini)%strangeness - hadron(idbaryon_ini)%charm
       i3  = hadron(resID)%isoSpinTimes2
       iz3 = iz1 + iz2

       if (debug) write(*,*) i1,i2,i3,iz1,iz2
       isoFactor = CG(i1,i2,i3,iz1,iz2,iz3)**2
       if (isoFactor<1E-3) cycle   ! discard isospin-forbidden states

       ! Evaluate the spin factors
       spinFactor=(2.*hadron(resId)%spin+1.)/(2.*hadron(idBaryon_ini)%spin+1.)/(2.*hadron(idMeson_ini)%spin+1.)

       ! Evaluate mass of the resonance and store it
       chargeRes=chargeMeson_ini+chargeBaryon_ini
       if (.not. getRMF_flag()) then
         masses(resID) = resonanceMass (resID, chargeRes, mediumAtColl, momLRF, position, perturbative)
         if (.not. masses(resID)>0.) cycle   ! this can happen if the potential gets too strong
       else if (present(srts)) then
         masses(resID)=srts
       else
         write(*,*) ' In barMes2resonance: srts must present in RMF mode !!!'
         stop
       end if

       if (debug) write(*,*) masses(resID), resID

       momCM2 = pCM_sqr (masses(resID)**2, mesonMass**2, baryonMass**2)    ! Evaluate CM-Momentum
       if (momCM2 <= 0.) cycle

       if (debug) write(*,*) 'momCM2=',momCM2,'mass=',masses(ResID),'mesmass=', mesonMass,'barmass=', baryonMass

       ! use inWidth since there are both masses of the incoming particles given
       ! And don't use medium since only the out-width shall be dressed due to the final states
       gamma_In = partialWidthBaryonMedium (resID,masses(resID),.true.,IDmeson_ini, &
                                            IDbaryon_ini,momLRF,vacuum,baryonMass,mesonMass)
       gammaTot = WidthBaryonMedium (resID,masses(resID),momLRF,mediumAtColl)

       if (abs(gamma_In)<1E-6 .or. abs(gammaTot)<1E-6) cycle

       if (resID==Delta.and.debug) write(300,'(6F8.4)')  momCM2, gamma_In, gammaTot, masses(resID)

       if (fullPropagator.and.isNonExotic(resID)) then
          PI_imag=selfenergy_imag(resID,absMom_LRF,momLRF(0),mediumAtColl)
          PI_real=   get_RealPart(resID,absMom_LRF,fullmass      ,mediumAtColl)
          sigma(resID) = isoFactor*spinFactor*4.*pi/momCM2*masses(ResID)*gamma_In*(-PI_imag) &
                         / ((fullMass**2-hadron(resID)%mass**2- PI_real)**2+PI_imag**2) / GeVSquared_times_mb
       else
          ! neglect real part: old Effenberger ansatz
          sigma(resID) = isoFactor*spinFactor*4.*pi/momCM2*masses(ResID)**2*gamma_In*gammaTot &
                         / ((masses(resID)**2-hadron(resID)%mass**2)**2+gammaTot**2*masses(resID)**2) / GeVSquared_times_mb
       end if

    end do
  end function barMes2resonance



  !****************************************************************************
  !****f* resonanceCrossSections/barMes_R_barMes
  ! NAME
  ! real function barMes_R_barMes (idMes_in, idBar_in, idMes_out, idBar_out, chargeMes_in, chargeBar_in, chargeMes_out, chargeBar_out,
  ! background, propagated, MediumAtCollision, momentumLRF, mesonMass_in, baryonMass_in,
  ! position, perturbative, srts)
  !
  ! PURPOSE
  ! Evaluates contribution of resonances to baryon meson -> baryon meson according to the Breit-Wigner formula,
  ! cf. eq. (2.52) in Effenberger Phd.
  ! One gets B m -> B m by summing over the (b M -> R->b M) channels for all resonances.
  !
  ! INPUTS
  ! * integer, intent(in) :: idMes_in, idBar_in, idMes_out, idBar_out                 --- IDs of particles in initial and final state.
  ! * integer, intent(in) :: chargeMes_in, chargeBar_in, chargeMes_out, chargeBar_out --- Charges of particles in initial and final state.
  ! * logical, intent(in) :: background ---
  !   .true. = Compute background cross section : Therefore consider only non-propagated particles in intermediate state;
  !   .false. = Full Xsection : Consider all particles in intermediate state
  ! * logical, intent(in) :: propagated ---
  !   .true. = Consider only propagated particles in intermediate state;
  !   .false. = Consider all particles in intermediate state
  !
  ! Information for In-Medium modifications :
  ! * type(medium), intent(in) :: mediumAtCollision      ! Medium information at the point of collision
  ! * real, intent(in),dimension(0:3) :: momentumLRF     ! Total  Momentum in LRF=Momentum of resonance in LRF
  ! * real, intent(in) :: mesonMass_in, baryonMass_in    ! masses of colliding particles
  ! * real, intent(in),dimension(1:3) :: position        ! Position where resonance is produced
  ! * logical, intent(in)              :: perturbative   ! flag for the perturbative nature of the resonance
  ! * real, intent(in), optional                :: srts  ! sqrt(s) of the collision -- needed in RMF mode
  !
  ! NOTES
  ! Note that it is not useful to set both background and propagated to .true. : Result=0.
  ! RESULT
  ! * resonance contribution to baryon meson -> baryon meson (in mb)
  !****************************************************************************
  function barMes_R_barMes (idMes_in, idBar_in, idMes_out, idBar_out, chargeMes_in, chargeBar_in, chargeMes_out, chargeBar_out, &
                            background, propagated, mediumAtColl, momLRF, mesonMass_in, baryonMass_in, &
                            position, perturbative, srts) result (sigma)

    use baryonWidthMedium, only: WidthBaryonMedium, partialWidthBaryonMedium
    use mediumDefinition, only: medium, vacuum
    use twoBodyTools, only: pCM_sqr
    use particleProperties, only: hadron, isNonExotic
    use ClebschGordan, only: CG
    use constants, only: pi, GeVSquared_times_mb
    use RMF, only: getRMF_flag
    use selfenergy_baryons, only: get_realPart, selfEnergy_imag
    use minkowski, only: abs3, abs4
    use IdTable, only: Delta, nbar, isBaryon, isMeson

    integer, intent(in)        :: idMes_in, idBar_in, idMes_out, idBar_out, &
                                  chargeMes_in, chargeBar_in, chargeMes_out, chargeBar_out
    logical, intent(in)        :: background, propagated
    type(medium), intent(in)   :: mediumAtColl
    real, intent(in)           :: momLRF(0:3), mesonMass_in, baryonMass_in, position(1:3)
    logical, intent(in)        :: perturbative
    real, intent(in), optional :: srts
    real :: sigma

    integer :: resID, chargeRes, i1, i2, i3, iz1, iz2, iz3
    real :: momCM2, Gamma_In, Gamma_Out, Gamma_Tot, spinFactor, isoFactor, massRes, absMom_LRF, fullMass, PI_REAL, PI_IMAG

    if (initInput_Flag) call readInput()

    sigma = 0.

    ! Check Input
    if (.not.(isMeson(idMes_in) .and. isMeson(idMes_out) .and. isBaryon(idBar_in) .and. isBaryon(idBar_out))) then
       write(*,*) 'Error in in input of "resonanceXsection" '
       write(*,*) idMes_in, idBar_in, idMes_out, idBar_out
       stop
    end if

    if (abs(chargeMes_in+chargeBar_in-chargeMes_out-chargeBar_out) > 0) return        ! charge conservation

    if (abs(hadron(idMes_in)%strangeness +hadron(idBar_in)%strangeness &
           -hadron(idMes_out)%strangeness-hadron(idBar_out)%strangeness) > 0) return  ! strangeness conservation

    if (abs(hadron(idMes_in)%charm+hadron(idBar_in)%charm &
           -hadron(idMes_out)%charm-hadron(idBar_out)%charm) > 0) return              ! charm conservation

    chargeRes = chargeMes_in + chargeBar_in
    absMom_LRF = abs3(momLRF)
    if (fullPropagator) fullmass = abs4(momLRF)

    ! Loop over intermediate resonances R : m B -> R -> m' B'
    do resID = Delta, nbar

       if (.not. hadron(resID)%usedForXsections) cycle             ! Exclude resonances
       if (background .and. hadron(resID)%propagated) cycle        ! Exclude propagated resonances for background
       if (propagated .and. .not. hadron(resID)%propagated) cycle  ! Exclude non-propagated resonances, use only the propagated ones

       ! check conservation laws
       if (abs(hadron(idMes_in)%strangeness + hadron(idBar_in)%strangeness - hadron(resID)%strangeness) > 0) cycle   ! strangeness conservation
       if (abs(hadron(idMes_in)%charm       + hadron(idBar_in)%charm       - hadron(resID)%charm)       > 0) cycle   ! charm conservation

       ! Convert charges to isospins in SU(4) & evaluate isospin factors
       ! Q=Y/2+I_3 with Y=S+C+B  =>  I_3=Q-(S+C+B)/2.

       ! (a) Initial State
       i1  = hadron(idMes_in)%isoSpinTimes2
       iz1 = nint(2.*chargeMes_in) - hadron(idMes_in)%strangeness - hadron(idMes_in)%charm
       i2  = hadron(idBar_in)%isoSpinTimes2
       iz2 = nint(2.*chargeBar_in-1.) - hadron(idBar_in)%strangeness - hadron(idBar_in)%charm
       i3  = hadron(resID)%isoSpinTimes2
       iz3 = iz1 + iz2
       if (debug) then
         write(*,*) 'In resonanceCrossSections'
         write(*,*) i1,i2,i3,iz1,iz2
       end if
       isoFactor = CG(i1,i2,i3,iz1,iz2,iz3)**2
       if (isoFactor<1E-3) cycle   ! discard isospin-forbidden states

       ! (b) Final State
       i1  = hadron(idMes_out)%isoSpinTimes2
       iz1 = nint(2.*chargeMes_out) - hadron(idMes_out)%strangeness - hadron(idMes_out)%charm
       i2  = hadron(IDBar_out)%isoSpinTimes2
       iz2 = nint(2.*chargeBar_out-1.) - hadron(idBar_out)%strangeness - hadron(idBar_out)%charm
       i3  = hadron(resID)%isoSpinTimes2
       iz3 = iz1 + iz2
       isoFactor = isoFactor * CG(i1,i2,i3,iz1,iz2,iz3)**2
       if (isoFactor<1E-3) cycle   ! discard isospin-forbidden states

       ! determine resonance mass
       if (.not. getRMF_flag()) then
         massRes = resonanceMass (resID, chargeRes, mediumAtColl, momLRF, position, perturbative)
         if (.not. massRes>0.) cycle   ! this can happen if the potential gets too strong
       else if (present(srts)) then
         massres = srts
       else
         write(*,*) ' In barMes_R_barMes: srts must present in RMF mode !!!'
         stop
       end if
       if (debug) write(*,*) massres

       momCM2 = pCM_sqr (massRes**2, mesonMass_in**2, baryonMass_in**2)    ! Evaluate CM-Momentum
       if (momCM2 <= 0.) cycle

       ! Evaluate the spin factors
       spinFactor = (2.*hadron(resID)%spin+1.) / ((2.*hadron(idBar_in)%spin+1.)*(2.*hadron(idMes_in)%spin+1.))

       ! Evaluate the partial decay width for incoming and outgoing channel.
       ! Use inWidth since there are both masses of the incoming particles given.
       ! And don't use medium since only the out-width shall be dressed due to the final states.

       Gamma_In = partialWidthBaryonMedium (resID, massRes, .true., idMes_in, idBar_in, &
                                            momLRF, vacuum, baryonMass_in, mesonMass_in)

       Gamma_Out = partialwidthBaryonMedium (resID, massRes, .false., idMes_out, idBar_out, &
                                             momLRF, mediumAtColl)

       Gamma_Tot = WidthBaryonMedium (resID, massRes, momLRF, mediumAtColl)

       ! Evaluate cross section
       if (fullPropagator.and.isNonExotic(resID)) then
          PI_imag=selfenergy_imag(resID,absMom_LRF,momLRF(0),mediumAtColl)
          PI_real=   get_RealPart(resID,absMom_LRF,fullmass      ,mediumAtColl)
          sigma = sigma + isoFactor*spinFactor*4.*pi/momCM2*massRes**2*Gamma_In*Gamma_Out &
                          / ((fullMass**2-hadron(resID)%mass**2-PI_real)**2+PI_imag**2) / GeVSquared_times_mb
       else
          ! neglect real part: old Effenberger ansatz
          sigma = sigma + isoFactor*spinFactor*4.*pi/momCM2*massRes**2*Gamma_In*Gamma_Out &
                          / ((massRes**2-hadron(resID)%mass**2)**2+Gamma_Tot**2*massRes**2) / GeVSquared_times_mb
       end if

    end do

  end function barMes_R_barMes


  !****************************************************************************
  !****s* resonanceCrossSections/sigma_npi_n2pi_resonances
  ! NAME
  ! subroutine sigma_npi_n2pi_resonances(srts,charge_iniPion,background,sigmaTotal)
  ! PURPOSE
  ! Evaluates contribution of resonances to proton Pion -> Nucleon Pion Pion in the VACUUUM assumption [Gamma's are not the medium modified ones]
  ! One gets pion p -> pi pi N by summing over the (pion p -> R->pion Delta), (pion p -> R-> N rho),
  ! (pion p -> R-> N sigma) and (pion p -> R-> pion P11_1440) channels.
  ! All resonances but the P_11(1440) decay fully into N pi. Therefore the final states are in the end
  ! N Pi Pi states. For the P_11(1440) the ratio of N pi pi final states is given by its decay ratio into N pi.
  !
  ! INPUTS
  ! * logical, intent(in) :: background ---
  !   .true. = Compute background cross section : Therefore consider only non-propagated particles in intermediate state;
  !   .false. = Full Xsection : Consider all particles in intermediate state
  ! * real, intent(in) :: srts ---
  !   sqrt(s) in the process
  ! * integer, intent(in) :: charge_iniPion ---
  !   charge of incoming pion
  !
  ! RESULT
  ! * real, intent(out),dimension(-2:2) :: sigmaTotal ---
  !   cross sction in mB
  !
  ! Meaning of index in sigmaTotal :
  ! * -2 : pi+ pi- in final state
  ! * -1 : pi0 pi- in final state
  ! *  0 : pi0 pi0 in final state
  ! *  1 : pi+ pi0 in final state
  ! *  2 : pi+ pi+ in final state
  !****************************************************************************
  subroutine sigma_npi_n2pi_resonances (srts, charge_iniPion, background, sigmaTotal)
      use baryonWidth, only: partialwidthBaryon, FullWidthBaryon
      use particleProperties, only: hadron, isNonExotic
      use IdTable, only: nucleon, Delta, P11_1440, nbar, pion, rho, sigmaMeson
      use constants, only: pi, mN, mPi, GeVSquared_times_mb
      use selfenergy_baryons, only: get_realPart, selfEnergy_imag
      use mediumDefinition

      real, intent(in) :: srts
      integer, intent(in)  :: charge_iniPion
      logical, intent(in) :: background
      real, intent(out),dimension(-2:2) :: sigmaTotal

      real, dimension(1:4) :: gamma_Out,sigma
      real :: gammaTot, momCM, gamma_In, spinFactor, absMom_LRF, fullMass, PI_REAL, PI_IMAG, energy
      integer :: resID,decID
      type(medium) :: med

      if (initInput_Flag) call readInput()

      absMom_LRF=sqrt(((srts**2-mN**2-mPi**2)/(2*mN))**2-mPi**2)
      if (fullPropagator) then
         fullmass  =srts
         energy=sqrt(srts**2+absMom_LRF**2)
         med%temperature    =0.
         med%useMedium      =.true.
         med%density        = 0
         med%densityProton  = 0
         med%densityNeutron = 0
      end if

      !check input
      if (.not.(abs(charge_iniPion).le.1)) then
         write(*,*) 'Problems with input in npi_n2pi_resonances', charge_iniPion
      end if

      !Evaluate CM-Momentum
      momCM=Sqrt(max((srts**2-(mPi+mN)**2) &
           &                      *(srts**2-(mPi-mN)**2)  /4./srts**2,1E-8))


      ! Evaluate cross section for all channels
      sigmaTotal(-2:2)=0.
      ! Loop over intermediate resonances R : m B->  R -> m' B'
      do resId=1,nbar
         if (.not.hadron(resId)%usedForXsections) cycle   ! Exclude resonances
         if (background.and.hadron(resId)%propagated) cycle ! Exclude propagated resonances for background contribution

         ! Evaluate the spin factors
         spinFactor=(2.*hadron(resId)%spin+1.)/(2.*hadron(nucleon)%spin+1.)/(2.*hadron(pion)%spin+1.)
         ! Evaluate the partial decay width for incoming and outgoing channel
         gamma_In=partialwidthBaryon(resID,srts,.true.,pion,nucleon)
         ! Get pion N -> pi pi N by summing just over the ...
         ! pion N -> pion Delta, pion N -> N rho, pion N -> N sigma, pion N -> pion P11_1440
         ! ... channels. All resonances but the P_11(1440) decay fully into N pi.
         ! We correct for the P_11(1440) later.
         gamma_Out(1)=partialwidthBaryon(resID,srts,.false.,pion,Delta)
         gamma_Out(2)=partialwidthBaryon(resID,srts,.false.,rho,nucleon)
         gamma_Out(3)=partialwidthBaryon(resID,srts,.false.,sigmaMeson,nucleon)
         gamma_Out(4)=partialwidthBaryon(resID,srts,.false.,pion,P11_1440)
         ! Formula(2.52) Effenberger Phd.for each channel

         gammaTot=FullWidthBaryon(resID,srts)

         do decId=1,4
            if (fullPropagator.and.isNonExotic(resID)) then
               PI_imag=selfenergy_imag(resID,absMom_LRF,energy        ,med)
               PI_real=   get_RealPart(resID,absMom_LRF,fullmass      ,med)
               sigma(decID) = spinFactor*4.*pi/momCm**2*srts**2*Gamma_In*Gamma_Out(decID) &
                              / ((fullMass**2-hadron(resID)%mass**2-PI_real)**2+PI_imag**2) / GeVSquared_times_mb
            else
               ! neglect real part: old Effenberger ansatz
               sigma(decID) = spinFactor*4.*pi/momCm**2*srts**2*Gamma_In*Gamma_Out(decID) &
                              / ((srts**2-hadron(resID)%mass**2)**2+gammaTot**2*srts**2) / GeVSquared_times_mb
            end if
         end do
         ! correction for N(1440)
         ! multiply the cross section by the N_1440 decay ratio into Nucleon Pion
         ! P11(1440) is the only resonance which has different decay channels than n pi
         ! other channels of P11(1440) might lead to 3 pion nucleon final states.
         sigma(4)=sigma(4)*hadron(P11_1440)%decays(1)

         ! differentiate different isospinChannels and make isospin calculus for the different channnels
         select case (hadron(resID)%isoSpinTimes2)
         case (1) ! Isospin=1/2
            ! not possible to scatter pi^{+} proton via a I=1/2 channel, only contributions of pi^{-} and pi^{0}
            if (charge_iniPion.eq.0) then
               sigmaTotal(0)=sigmaTotal(0)  + 1./3.* (2./9.*sigma(1)+1./3.*sigma(3)+ 1./9.*sigma(4))
               sigmaTotal(-2)=sigmaTotal(-2)+ 1./3.* (5./9.*sigma(1)+1./3.*sigma(2)+ 2./3.*sigma(3)+4./9.*sigma(4))
               sigmaTotal(1)=sigmaTotal(1)  + 1./3.* (2./9.*sigma(1)+2./3.*sigma(2)+ 4./9.*sigma(4))
            else if (charge_iniPion.eq.-1) then
               sigmaTotal(0)=sigmaTotal(0)  + 2./3.* (2./9.*sigma(1)+1./3.*sigma(3)+ 1./9.*sigma(4))
               sigmaTotal(-2)=sigmaTotal(-2)+ 2./3.* (5./9.*sigma(1)+1./3.*sigma(2)+ 2./3.*sigma(3)+4./9.*sigma(4))
               sigmaTotal(-1)=sigmaTotal(-1) + 2./3.* (2./9.*sigma(1)+2./3.*sigma(2)+ 4./9.*sigma(4))
            end if
         case (3) ! Isospin=3/2
            if (charge_iniPion.eq.0) then
                sigmaTotal(0)=sigmaTotal(0)  +  2./3.* (2./45.*sigma(1)+2./9.*sigma(4))
                sigmaTotal(-2)=sigmaTotal(-2)+  2./3.* (26./45.*sigma(1)+2./3.*sigma(2)+2./9.*sigma(4))
                sigmaTotal(1)=sigmaTotal(1)  +  2./3.* (17./45.*sigma(1)+1./3.*sigma(2)+5./9.*sigma(4))
            else if (charge_iniPion.eq.-1) then
                sigmaTotal(0)=sigmaTotal(0)+   1./3.*(2./45.*sigma(1)+2./9.*sigma(4))
                sigmaTotal(-2)=sigmaTotal(-2)+ 1./3.*(26./45.*sigma(1)+2./3.*sigma(2)+2./9.*sigma(4))
                sigmaTotal(-1)=sigmaTotal(-1)+  1./3.*(17./45.*sigma(1)+1./3.*sigma(2)+5./9.*sigma(4))
            else if (charge_iniPion.eq.1) then
                sigmaTotal(1)=sigmaTotal(1)+    13./15.*sigma(1)+sigma(2)+1./3.*sigma(4)
                sigmaTotal(2)=sigmaTotal(2)+     2./15.*sigma(1)+2./3.*sigma(4)
            end if
         end select
      end do
   end subroutine sigma_npi_n2pi_resonances

  end module resonanceCrossSections
