!******************************************************************************
!****m* /baryonWidthVacuum
! NAME
! module baryonWidthVacuum
! NOTES
! Module which calculates the partial and full widths of the baryon resonances
! in dependence of their mass.
! Their pole mass is given by their mass in 'particleProperties' and the
! widths at this pole mass as well.
! Everything corresponds to the vacuum situation. The resulting width is
! therefore always the width in vacuum!
! As mass of the resonance we use the four-vector definition: p_mu p^mu= mass**2
! Prescription according to Manley et al. Phys. Rev. D45 (1992) 4002.
!******************************************************************************
module baryonWidthVacuum

  use constants, only: hbarc

  implicit none
  private

  !****************************************************************************
  !****g* baryonWidthVacuum/interactionRadius
  ! SOURCE
  !
  real, parameter, public :: interactionRadius = 1./hbarc
  ! PURPOSE
  ! interaction radius for Blatt-Weisskopf functions: r = 1fm
  !****************************************************************************

  !****************************************************************************
  !****g* baryonWidthVacuum/srts_srt_switch
  ! SOURCE
  !
  logical, parameter, public :: srts_srt_switch=.false.
  ! PURPOSE
  ! Modifies the width according to S. Leupold's definition of the width,
  ! one especially has to exchange s against sqrt(s) in the denominator of
  ! Formula 2.76 of Effenbergers Phd
  !****************************************************************************

  !****************************************************************************
  !****g* baryonWidthVacuum/deltaRho_cutoff
  ! SOURCE
  !
  real, save :: deltaRho_cutoff=0.85
  ! PURPOSE
  ! * Cut off parameter for the decay of a resonance into delta rho.
  ! * Units of GeV
  !****************************************************************************

  !****************************************************************************
  !****g* baryonWidthVacuum/meson_cutoff
  ! SOURCE
  !
  real, save :: meson_cutoff=1.6
  ! PURPOSE
  ! * Cut off parameter for the decay of a resonance into a baryon and an
  !   unstable meson.
  ! * Units of GeV
  !****************************************************************************

  !****************************************************************************
  !****g* baryonWidthVacuum/baryon_cutoff
  ! SOURCE
  !
  real, save :: baryon_cutoff=2.0
  ! PURPOSE
  ! * Cut off parameter for the decay of a resonance into an unstable baryon
  !   and a meson.
  ! * Units of GeV
  !****************************************************************************

  !****************************************************************************
  !****g* baryonWidthVacuum/use_cutoff
  ! SOURCE
  !
  logical, save :: use_cutoff = .true.
  ! PURPOSE
  ! * Switch on and off the use of cut off parameters.
  ! * These cut-offs are necessary when working with dispersion relations to
  !   deduce the real part.
  !****************************************************************************


  !****************************************************************************
  !****g* baryonWidthVacuum/Delta_width
  ! SOURCE
  !
  integer, save :: Delta_width = 1
  ! PURPOSE
  ! Select a parametrization for the Delta width:
  ! * 1 = Manley   (GiBUU default, cf. Manley/Saleski, Phys. Rev. D 45, 1992)
  ! * 2 = Dmitriev (Dmitriev/Sushkov/Gaarde, Nucl. Phys. A 459, 1986)
  ! * 3 = Moniz    (Koch/Moniz/Ohtsuka, Ann. of Phys. 154, 1984)
  ! * 4 = Verwest  (Phys. Lett. B 83, 1979)
  ! * 5 = UrQMD    (Bass et al., Prog. Part. Nucl. Phys. 41, 1998)
  !****************************************************************************


  logical, save :: initFlag=.true.


  public :: vacuumWidth
  public :: getParameters


contains


  !****************************************************************************
  !****s* baryonWidthVacuum/getParameters
  ! NAME
  ! subroutine getParameters(a,b,c,flag)
  ! PURPOSE
  ! return the values of ...
  !****************************************************************************
  subroutine getParameters(a,b,c,flag,delta)
    real,intent(out) :: a,b,c
    logical,intent(out) :: flag
    integer,intent(out) :: delta

    if (initFlag) call readinput

    a = deltaRho_cutoff
    b = meson_cutoff
    c = baryon_cutoff
    flag = use_cutoff
    delta = Delta_width
  end subroutine getParameters


  !****************************************************************************
  !****s* baryonWidthVacuum/readInput
  ! NAME
  ! subroutine readInput
  ! PURPOSE
  ! Reads input in jobcard out of namelist "BaryonWidthVacuum".
  !****************************************************************************
  subroutine readInput
    use output, only: Write_ReadingInput

    integer :: ios ! checks file behavior

    !**************************************************************************
    !****n*  baryonWidthVacuum/BaryonWidthVacuum
    ! NAME
    ! NAMELIST /BaryonWidthVacuum/
    ! PURPOSE
    ! Includes the input switches:
    ! * use_cutoff
    ! * deltaRho_cutoff
    ! * baryon_cutoff
    ! * meson_cutoff
    ! * Delta_width
    !**************************************************************************
    NAMELIST /baryonWidthVacuum/ deltaRho_cutoff, baryon_cutoff, meson_cutoff, use_cutoff, Delta_width

    call Write_ReadingInput('baryonWidthVacuum',0)
    rewind(5)
    read(5,nml=baryonWidthVacuum,IOSTAT=IOS)
    call Write_ReadingInput('baryonWidthVacuum',0,IOS)

    if (use_cutOff) then
       write(*,*) 'delta-Rho cutoff =         ', deltaRho_Cutoff
       write(*,*) 'meson cutoff =             ', meson_Cutoff
       write(*,*) 'baryoncutoff =             ', baryon_Cutoff
    else
       write(*,*) 'No additional cutoffs'
    end if
    write(*,*) 'Delta width parametrization: ', Delta_width

    call Write_ReadingInput('baryonWidthVacuum',1)
    initFlag = .false.

  end subroutine readInput


  !****************************************************************************
  !****if* baryonWidthVacuum/monopoleFormfactor
  ! NAME
  ! real function monopoleFormfactor(mass,pole,lambda)
  ! PURPOSE
  ! ...
  !****************************************************************************
!   real function monopoleFormfactor (mass, pole, lambda)
!     real, intent(in) :: mass, pole, lambda
!
!     if (use_cutOff) then
!        monopoleFormfactor=lambda**4/(lambda**4+(pole**2-mass**2)**2)
!     else
!        monopoleFormfactor=1.
!     end if
!   end function monopoleFormfactor



  !****************************************************************************
  !****s* baryonWidthVacuum/vacuumWidth
  ! NAME
  ! function vacuumWidth (mass, partID, ratio, rho_AB_atPole) result (gammaTotal)
  ! PURPOSE
  ! This routine calculates the total vacuum decay width of a baryonic
  ! resonance, i.e. the sum of all partial decay widths, and the branching
  ! ratios for each channel as a function of the offshell mass (m^2 = p_mu p^mu).
  ! Parameters taken from Manley et al. Phys. Rev. D45 (1992) 4002 and PDG.
  ! The mass dependence of the resonances is treated according to Manley.
  ! INPUTS
  ! * real :: mass        -- offshell mass of Resonance in GeV
  ! * integer :: partID   -- ID of resonance
  ! OUTPUT
  ! * real, dimension(1:nDecays) :: ratio -- branching ratio = partial width / full width
  !   for all 2 body decay channels of the baryons
  ! * real, dimension(1:nDecays) :: rho_AB_atPole --
  !   rho_ab(poleMass) of specific decay channel according
  !   Effenberger Dr. thesis, eq 2.76
  ! * real :: gammaTotal -- full width in GeV
  !****************************************************************************
  function vacuumWidth (mass, partID, ratio, rho_AB_atPole) result (gammaTotal)
    use DecayChannels, only: Decay2BodyBaryon
    use particleProperties, only: hadron, nDecays, getAngularMomentum_baryon
    use IdTable, only: Delta

    real, intent(in) :: mass
    integer, intent(in) :: partID
    real, intent(out), dimension(nDecays) :: ratio, rho_AB_atPole
    real :: gammaTotal

    real, parameter :: widthCutOff = 1E-10
    real, dimension(nDecays) :: gamma
    real :: partialWidth, poleMass, mass1,mass2
    integer :: L, idStable, idUnStable, k, dID
    logical, parameter :: debug=.false.

    if (initFlag) call readinput

    gamma=0.
    do k=1,nDecays
       dID = hadron(partID)%decaysID(k)
       if (dID<=0) cycle ! 2Body: >0, 3Body: <0

       partialWidth = hadron(partID)%decays(k) * hadron(partID)%width  !in GeV
       if (partialWidth < widthCutOff) cycle

       poleMass=hadron(partID)%mass
       mass1=hadron(Decay2BodyBaryon(dID)%ID(1))%mass  !mass of produced meson
       mass2=hadron(Decay2BodyBaryon(dID)%ID(2))%mass !mass of produced baryon
       ! Set Angular Momentum of decaying particles
       L = getAngularMomentum_baryon(dID,partID)

       if (debug) write(*,*) "Masses", mass1,mass2, polemass, "angular=",L

       if (Decay2BodyBaryon(dID)%stable(1).and.Decay2BodyBaryon(dID)%stable(2)) then
          !  decay into 2 stable particles

          if (partID==Delta .and. dID==1) then
            ! Delta -> pi N: choose between different parametrizations
            select case (Delta_width)
            case (1) ! default: Manley
              call stableFinalState (mass, poleMass, mass1, mass2, L, partialWidth, gamma(k), rho_AB_atPole(k))
            case (2) ! Dmitriev
              call Delta_Dmitriev (mass, partialWidth, gamma(k), rho_AB_atPole(k))
            case (3) ! Moniz
              call Delta_Moniz (mass, partialWidth, gamma(k), rho_AB_atPole(k))
            case (4) ! Verwest
              call Delta_Verwest (mass, partialWidth, gamma(k), rho_AB_atPole(k))
            case (5) ! UrQMD
              call Delta_Bass (mass, partialWidth, gamma(k), rho_AB_atPole(k))
            case default
              write(*,*) "bad value of Delta_width: ", Delta_width
              stop
            end select
          else
            ! Manley parametrization
            call stableFinalState (mass, poleMass, mass1, mass2, L, partialWidth, gamma(k), rho_AB_atPole(k))
          end if

          if (debug) write(*,*) "Stable, decayChannel:",dID,"Gamma=",gamma(k)

       else if (Decay2BodyBaryon(dID)%stable(1).neqv.Decay2BodyBaryon(dId)%stable(2)) then
          !  decay into 1 stable particle + 1 unstable particle
          if (Decay2BodyBaryon(dID)%stable(1)) then
             ! particle 1 is stable
             IdStable=Decay2BodyBaryon(dID)%ID(1)
             IdUnStable=Decay2BodyBaryon(dID)%ID(2)
          else
             ! particle 2 is stable
             IdStable=Decay2BodyBaryon(dID)%ID(2)
             IdUnStable=Decay2BodyBaryon(dID)%ID(1)
          end if
          call semiStableFinalState(mass,polemass,IdStable,IdUnstable,L,partialWidth,gamma(k),rho_AB_atPole(k))

          if (debug) write(*,*) "Stable=", idStable, "Unstable=", idunStable,"Gamma=", gamma(k)

       else
          !   decay into 2 unstable particles
          if (dID<15 .or. dID>18) then
             ! no delta rho channel!!!
             write(*,*) 'critical error in manley'
             write(*,*) 'two unstable particles in the final state, but they are not delta and rho'
             write(*,*) 'DecayChannel:', k
             write(*,*) 'Stopping BUU'
             stop
          else
             call rhoDeltaFinalState(mass,poleMass,L,partialWidth,gamma(k),rho_AB_atPole(k))
          end if
          if (debug) write(*,*) "DeltaRho, Gamma=", gamma(k)
       end if
    end do

    gammaTotal=Sum(gamma(:))

    if (debug) write(*,*) gammaTotal

    if (gammaTotal /= 0) then
       ratio(:)=gamma(:)/gammaTotal
    else
       ratio(:)=0.
    end if

    if (debug) write(*,*)

  end function vacuumWidth


  !****************************************************************************
  !****is* baryonWidthVacuum/stableFinalState
  ! NAME
  ! subroutine stableFinalState (mass, poleMass, mass1, mass2, L, partialWidth_pole, partialWidth_mass, rho_AB_Pole)
  ! PURPOSE
  ! calculate the partial width for the decay into two stable particles.
  ! Manley et al. Phys. Rev. D45 (1992) 4002.
  ! NOTES
  ! decay only allowed if:
  ! * (1) decay width>widthCutOff and ...
  ! * (2) mass of resonance > Sum off masses of decay products
  ! INPUTS
  ! * integer :: L              -- Angular Momentum
  ! * real :: polemass          -- pole mass of mother resonance
  ! * real :: mass              -- mass of mother resonance
  ! * real :: mass1             -- mass of first decay product
  ! * real :: mass2             -- mass of second decay product
  ! * real :: partialWidth_pole -- partial width of decayChannel at pole mass
  ! OUTPUT
  ! * real :: partialWidth_mass -- partial width at real mass
  ! * real :: rho_AB_Pole       -- Equation 2.76 of Effe Dr. evaluated at pole
  !****************************************************************************
  subroutine stableFinalState (mass, poleMass, mass1, mass2, L, partialWidth_pole, partialWidth_mass, rho_AB_Pole)
    use distributions, only: BlattWeisskopf
    use twoBodyTools, only: pCM

    real, intent(in) :: mass, poleMass, mass1, mass2
    integer, intent(in) :: L
    real, intent(in) :: partialWidth_pole
    real, intent(out) :: partialWidth_mass, rho_AB_Pole

    real :: p_ab_mass, p_ab_pole, rho_ab_mass

    if (initFlag) call readinput

    ! Determine momentum of outgoing particles in Restframe of Resonance
    p_ab_mass = pCM(mass,mass1,mass2)
    p_ab_pole = pCM(polemass,mass1,mass2)

    ! Evaluate rho_ab according to equation 2.76 in Effenbergers Dr. thesis
    ! rho_ab(mu)=p_ab/mu * BlattWeisskopf(pab*interactionRadius,L)
    if (srts_srt_switch) then
       rho_ab_mass=p_ab_mass/mass**2*(BlattWeisskopf(p_ab_mass*interactionRadius,L))**2
       rho_ab_pole=p_ab_pole/polemass**2*(BlattWeisskopf(p_ab_pole*interactionRadius,L))**2
    else
       rho_ab_mass=p_ab_mass/mass*(BlattWeisskopf(p_ab_mass*interactionRadius,L))**2
       rho_ab_pole=p_ab_pole/polemass*(BlattWeisskopf(p_ab_pole*interactionRadius,L))**2
    end if

    if (mass <= mass1 + mass2) then
       partialWidth_Mass=0.
    else
       partialWidth_Mass=partialWidth_pole*rho_ab_mass/rho_ab_pole
    end if
  end subroutine stableFinalState


  !****************************************************************************
  !****is* baryonWidthVacuum/Delta_Dmitriev
  ! NAME
  ! subroutine Delta_Dmitriev (mass, partialWidth_pole, partialWidth_mass, rho_AB_pole)
  ! PURPOSE
  ! Calculate the Delta width according to:
  ! Dmitriev/Sushkov/Gaarde, Nucl. Phys. A 459 (1986) 503-524
  ! INPUTS
  ! * real :: mass              -- mass of  mother resonance
  ! * real :: partialWidth_pole -- partial width of decayChannel at pole mass
  ! OUTPUT
  ! * real :: partialWidth_mass -- partial width at real mass
  ! * real :: rho_AB_Pole       -- rho_ab as defined in eq. (2.75) of Effenberger PhD evaluated at pole
  !****************************************************************************
  subroutine Delta_Dmitriev (mass, partialWidth_pole, partialWidth_mass, rho_AB_pole)
    use particleProperties, only: hadron
    use IdTable, only: Delta
    use constants, only: mN, mPi
    use twoBodyTools, only: pCM

    real, intent(in) :: mass, partialWidth_pole
    real, intent(out) :: partialWidth_mass, rho_AB_pole

    real :: k,k0,rho_AB_mass
    real, parameter :: kappa = 0.2

    k  = pCM(mass,mN,mPi)
    k0 = pCM(hadron(Delta)%mass,mN,mPi)

    rho_AB_mass = k**3  / (k**2  + kappa**2)
    rho_AB_pole = k0**3 / (k0**2 + kappa**2)

    if (mass < mN+mPi) then
      partialWidth_mass = 0.
    else
      partialWidth_mass = partialWidth_pole * rho_AB_mass / rho_AB_pole
    end if

  end subroutine Delta_Dmitriev


  !****************************************************************************
  !****is* baryonWidthVacuum/Delta_Moniz
  ! NAME
  ! subroutine Delta_Moniz (mass, partialWidth_pole, partialWidth_mass, rho_AB_pole)
  ! PURPOSE
  ! Calculate the Delta width according to:
  ! Koch/Moniz/Ohtsuka, Ann. of Phys. 154 (1984) 99
  ! INPUTS
  ! * real :: mass              -- mass of  mother resonance
  ! * real :: partialWidth_pole -- partial width of decayChannel at pole mass
  ! OUTPUT
  ! * real :: partialWidth_mass -- partial width at real mass
  ! * real :: rho_AB_Pole       -- rho_ab as defined in eq. (2.75) of Effenberger PhD evaluated at pole
  !****************************************************************************
  subroutine Delta_Moniz (mass, partialWidth_pole, partialWidth_mass, rho_AB_pole)
    use particleProperties, only: hadron
    use IdTable, only: Delta
    use constants, only: mN, mPi
    use twoBodyTools, only: pCM

    real, intent(in) :: mass, partialWidth_pole
    real, intent(out) :: partialWidth_mass, rho_AB_pole

    real :: k,k0,rho_AB_mass
    real, parameter :: d = 0.09

    k  = pCM(mass,mN,mPi)
    k0 = pCM(hadron(Delta)%mass,mN,mPi)

    rho_AB_mass = k**3  / mass               / (k**2  + d)**2
    rho_AB_pole = k0**3 / hadron(Delta)%mass / (k0**2 + d)**2

    if (mass < mN+mPi) then
      partialWidth_mass = 0.
    else
      partialWidth_mass = partialWidth_pole * rho_AB_mass / rho_AB_pole
    end if

  end subroutine Delta_Moniz


  !****************************************************************************
  !****is* baryonWidthVacuum/Delta_Verwest
  ! NAME
  ! subroutine Delta_Verwest (mass, partialWidth_pole, partialWidth_mass, rho_AB_pole)
  ! PURPOSE
  ! Calculate the Delta width according to:
  ! B.J. Verwest, Phys. Lett. B 83 (1979) 161.
  ! INPUTS
  ! * real :: mass              -- mass of  mother resonance
  ! * real :: partialWidth_pole -- partial width of decayChannel at pole mass
  ! OUTPUT
  ! * real :: partialWidth_mass -- partial width at real mass
  ! * real :: rho_AB_Pole       -- rho_ab as defined in eq. (2.75) of Effenberger PhD evaluated at pole
  !****************************************************************************
  subroutine Delta_Verwest (mass, partialWidth_pole, partialWidth_mass, rho_AB_pole)
    use particleProperties, only: hadron
    use IdTable, only: Delta
    use constants, only: mN, mPi
    use twoBodyTools, only: pCM

    real, intent(in) :: mass, partialWidth_pole
    real, intent(out) :: partialWidth_mass, rho_AB_pole

    real :: k,k0,rho_AB_mass

    k  = pCM(mass,mN,mPi)
    k0 = pCM(hadron(Delta)%mass,mN,mPi)

    rho_AB_mass = k**3  * (sqrt(k**2  + mPi**2)+mN)
    rho_AB_pole = k0**3 * (sqrt(k0**2 + mPi**2)+mN)

    if (mass < mN+mPi) then
      partialWidth_mass = 0.
    else
      partialWidth_mass = partialWidth_pole * rho_AB_mass / rho_AB_pole
    end if

  end subroutine Delta_Verwest


  !****************************************************************************
  !****is* baryonWidthVacuum/Delta_Bass
  ! NAME
  ! subroutine Delta_Bass (mass, partialWidth_pole, partialWidth_mass, rho_AB_pole)
  ! PURPOSE
  ! Calculate the Delta width as used in UrQMD, cf.:
  ! Bass et al., Prog. Part. Nucl. Phys. 41 (1998) 255-369.
  ! INPUTS
  ! * real :: mass              -- mass of  mother resonance
  ! * real :: partialWidth_pole -- partial width of decayChannel at pole mass
  ! OUTPUT
  ! * real :: partialWidth_mass -- partial width at real mass
  ! * real :: rho_AB_Pole       -- rho_ab as defined in eq. (2.75) of Effenberger PhD evaluated at pole
  !****************************************************************************
  subroutine Delta_Bass (mass, partialWidth_pole, partialWidth_mass, rho_AB_pole)
    use particleProperties, only: hadron
    use IdTable, only: Delta
    use constants, only: mN, mPi
    use twoBodyTools, only: pCM

    real, intent(in) :: mass, partialWidth_pole
    real, intent(out) :: partialWidth_mass, rho_AB_pole

    real :: k,k0,rho_AB_mass

    k  = pCM(mass,mN,mPi)
    k0 = pCM(hadron(Delta)%mass,mN,mPi)

    rho_AB_mass = k**3  / mass               / (1+0.2*(k/k0)**2)
    rho_AB_pole = k0**3 / hadron(Delta)%mass / 1.2

    if (mass < mN+mPi) then
      partialWidth_mass = 0.
    else
      partialWidth_mass = partialWidth_pole * rho_AB_mass / rho_AB_pole
    end if
  end subroutine Delta_Bass


  !****************************************************************************
  !****s* baryonWidthVacuum/semistableFinalState
  ! NAME
  ! subroutine semistableFinalState (mass, polemass, idStable, idUnstable, L, partialWidth_pole, partialWidth_mass, rho_AB_Pole)
  ! PURPOSE
  ! Calculate the partial width for the decay into
  ! one stable and one unstable particle.
  !
  ! Calculates the partial width dependend on mass of a resonance decaying into
  ! one stable and one unstable decay product.
  ! According to Manley and Effenberger' Dr. Thesis equation 2.76.
  ! INPUTS
  ! * real     :: mass -- mass of mother resonance
  ! * real     :: polemass -- polemass of mother resonance
  ! * integer  :: idUnstable -- ID of unstable decay product
  ! * integer  :: idStable -- ID of stable decay product
  ! * integer  :: L -- angular momentum of final state
  ! * real     :: partialWidth_pole -- partial width at pole of resonance
  ! OUTPUT
  ! * real :: partialWidth_mass -- partial width at real mass
  ! * real :: rho_AB_Pole -- Equation 2.76 of Effe Dr. evaluated at pole
  !****************************************************************************
  subroutine semistableFinalState (mass, polemass, idStable, idUnstable, L, partialWidth_pole, partialWidth_mass, rho_AB_Pole)
    use IdTable, only: isBaryon
    use ParticleProperties, only: hadron
    use rhoABIntegrandVariables, only: idUnstable_copy, L_copy, massStable, massUnStable, gammaUnStable, srts
    use quadpack, only: qag
    use distributions, only: markusPostFormFactor

    real, intent(in) :: mass, polemass
    integer, intent(in) :: idStable, idUnstable, L
    real, intent(in) :: partialWidth_pole
    real, intent(out) :: partialWidth_mass, rho_AB_Pole

    real :: lowerBound, upperBound    ! lower and upper bound of integrals
    real :: rho_ab_mass               ! Equation 2.76 Effenberger Dr.

    real :: cutOff
    logical,parameter :: debug=.false.

    real,   parameter :: qag_absError=0.
    real,   parameter :: qag_relError=0.00001
    integer,parameter :: qag_key=6
    real              :: qag_error_estimation
    integer           :: qag_neval,qag_error_code

    if (initFlag) call readinput

    ! setting for rho_ab_integrand
    L_copy=L
    idUnstable_copy=idUnstable

    massStable = hadron(idStable)%mass

    massUnStable=hadron(idUnStable)%mass
    gammaUnStable=hadron(idUnStable)%width

    if (isBaryon(idUnStable)) then
      cutOff = baryon_cutOff
    else
      cutOff = meson_cutOff
    end if

    if (debug) write(*,'(A,3F8.4)') "In semiStable:", massStable, massUnstable, gammaUnstable
    !    If(debug) write(*,'(A,2I4)'),"In semiStable:", idStable, idUnstable

    ! Use kinematics at the actual mass of the resonance :
    ! upperBound is the maximal mass of the unstable particle in the decay product.
    ! lowerBound is the minimal mass of the unstable particle (i.e. the rest masses of its stable decay products)
    srts=mass
    lowerBound=hadron(IdUnstable)%minmass
    upperBound=mass-massStable
    !    Print *, "lower&upperbound=", upperbound, lowerBound

    if (upperBound - lowerBound < 1E-12) then
      rho_ab_mass = 0.
    else
      call qag (rho_ab_Integrand, lowerBound, upperBound, qag_absError, qag_relError, qag_key, rho_ab_mass, &
                qag_error_estimation, qag_neval, qag_error_code)

      if (qag_error_code /= 0) then
        write(*,*) 'Error with QAG in baryonWidthVacuum (mass):', qag_error_code, qag_error_estimation, qag_neval
        write(*,*) mass, poleMass, IdStable, IdUnStable, L, lowerBound, upperBound, rho_ab_mass
      end if
    end if


    ! Use kinematics at the pole mass of the resonance :
    ! upperBound is the maximal mass of the unstable particle in the decay product.
    srts=polemass
    lowerBound=hadron(IdUnstable)%minmass
    upperBound=poleMass-massStable
    !    Print *, "lower&upperbound=", upperbound, lowerBound
    call qag (rho_AB_Integrand, lowerBound, upperBound, qag_absError, qag_relError, qag_key, rho_ab_pole, &
              qag_error_estimation, qag_neval, qag_error_code)

    if (qag_error_code /= 0) then
      write(*,*) 'Error with QAG in baryonWidthVacuum (pole):', qag_error_code, qag_error_estimation, qag_neval
      write(*,*) poleMass, IdStable, IdUnStable, L, lowerBound, upperBound, rho_ab_pole
    end if

    if (debug) then
       write(*,*) "Rho_AB_pole=", rho_ab_pole
       write(*,*) "Rho_AB_mass=", rho_ab_mass
    end if

    if ((mass <= massStable + lowerBound) .or. (rho_ab_pole==0.)) then
       partialWidth_Mass=0.
    else
       if (debug) write(*,*) partialWidth_pole,rho_ab_mass,rho_ab_pole
       partialWidth_mass = partialWidth_pole * rho_ab_mass / rho_ab_pole
       if (use_cutOff) partialWidth_mass = partialWidth_mass &
                       * markusPostFormFactor(mass,poleMass,lowerBound+massStable,cutOff)**2
    end if

  end subroutine semistableFinalState


  !****************************************************************************
  !****if* baryonWidthVacuum/rho_AB_Integrand
  ! NAME
  ! real function rho_AB_Integrand(mu)
  ! NOTES
  ! This is the integrand of rho_ab(mu) in Effenbergers Dr.-Thesis, page 27,
  ! formula 2.76.
  ! Only in the case of only one unstable particle!
  !****************************************************************************
  real function rho_AB_Integrand (mu)
    use IdTable
    use particleProperties, only: hadron
    use rhoABIntegrandVariables, only: idUnstable_copy, L_copy, massStable, massUnStable, gammaUnStable, srts
    use distributions, only: BlattWeisskopf,RelBW
    use constants, only: mN, mPi
    use twoBodyTools, only: pcm
    use mesonWidthVacuum, only: vacuumWidthMeson => vacuumWidth

    real, intent(in) :: mu ! Mass of unstable particle
    real :: p_AB, corrector, qcmUnstable_mu, qcmUnstable_pole,gammaCorrect
    !    real, parameter :: radiusRho=5.3
    logical, parameter :: debug=.false.

    ! "srtS" is the mass of mother Resonance. Set in the mother subroutine.

    ! It decays into a stable particle with mass "massstable" and
    ! an unstable particle with mass "mu".
    ! p_AB is the momentum of the two decay products in the restframe of the mother resonance.

    ! massStable, massUnstable and gammaUnstable are all defined in the mother subroutine.

    ! If p_ab imaginary then this kinematical situation is forbidden
    ! Therefore p_ab is then set to 0 which corresponds later to
    ! set rho_AB_integrand to zero.

    p_AB=pcm(srtS,mu,massStable)

    ! Redefine the gamma of the unstable particle, since the unstable particle is supposed off-shell:

    select case (idUnstable_copy)
    case (Delta) ! Delta -> pion + nuclon
       ! this delta has mass mu, so we need gammaUnstable(mu)
       ! The qcm momentum of pion and nucleon in the delta restframe are at the real mass "mu" of the delta:
       qcmUnstable_mu=pcm(mu,mN,mPi)
       ! And at the delta's polemass:
       qcmUnstable_pole=pcm(massUnstable,mN,mPi)

       corrector = qcmUnstable_mu/qcmUnstable_pole * massUnstable/mu &
                   * (BlattWeisskopf(qcmUnstable_mu*InteractionRadius,1)/BlattWeisskopf(qcmUnstable_pole*InteractionRadius,1))**2
       if (debug) write(*,*) "Correcting for delta width:",  gammaUnstable,corrector
       gammaCorrect=gammaUnstable*corrector

    case (rho,sigmaMeson,kaonStarBar)
       ! use vacuum width of meson
       gammaCorrect = vacuumWidthMeson (mu, idUnstable_copy)

    case (P11_1440,lambda_1520,lambda_1820)
       ! No Correction
       gammaCorrect=gammaUnstable

    case (Sigma_1385) ! Sigma(1385) -> pion + Lambda
       ! This Sigma* has mass mu, so we need gammaUnstable(mu).
       ! The qcm momentum of pion and Lambda in the restframe of Sigma* are at the real mass "mu" of the Sigma*:
       qcmUnstable_mu=pcm(mu,hadron(lambda)%mass,mPi)

       ! And at the Sigma* polemass:
       qcmUnstable_pole=pcm(massUnstable,hadron(lambda)%mass,mPi)

       corrector = qcmUnstable_mu/qcmUnstable_pole * massUnstable/mu &
                   * (BlattWeisskopf(qcmUnstable_mu*InteractionRadius,1)/BlattWeisskopf(qcmUnstable_pole*InteractionRadius,1))**2
       if (debug) write(*,*) "Correcting for delta width:",  gammaUnstable,corrector
       gammaCorrect=gammaUnstable*corrector

    case (lambda_1405) ! Lambda(1405) -> pion + Sigma
        qcmUnstable_mu=pcm(mu,mPi,hadron(SigmaResonance)%mass)
        qcmUnstable_pole=pcm(massUnstable,mPi,hadron(SigmaResonance)%mass)
        corrector=qcmUnstable_mu/qcmUnstable_pole * massUnstable/mu
        if (debug) write(*,*) "Correcting for  Lambda(1405) width:",  gammaUnstable,corrector
        gammaCorrect=gammaUnstable*corrector

    case default

       write(*,*) "Wrong Resonance as unstable Particle in rho_ab_integrand :",idUnstable_copy
       write(*,*) "This resonance is not yet implemented in subroutine 'semistableFinalState' :",idUnstable_copy

    end select

    ! Define gamma as function of mu
    if (srts_srt_switch) then
       rho_AB_Integrand = RelBW(mu,massUnStable,gammaCorrect) * p_AB/srtS**2 * BlattWeisskopf(p_AB*interactionRadius,L_copy)**2
    else
       rho_AB_Integrand = RelBW(mu,massUnStable,gammaCorrect) * p_AB/srtS * BlattWeisskopf(p_AB*interactionRadius,L_copy)**2
    end if

    if (debug) then
      write(*,*) RelBW(mu,massUnStable,gammaCorrect),p_AB,srtS,BlattWeisskopf(p_AB*interactionRadius,L_copy)**2
      write(*,*) "rho_AB_Integrand=",rho_AB_Integrand
    end if

  end function rho_AB_Integrand


  !****************************************************************************
  !****s* baryonWidthVacuum/rhoDeltaFinalState
  ! NAME
  ! subroutine rhoDeltaFinalState(mass,poleMass,L,partialWidth_pole,partialWidth_mass,rho_AB_atPole)
  ! PURPOSE
  ! Calculate Resonance decays into delta and rho.
  ! INPUTS
  ! * integer :: L -- Angular Momentum
  ! * real :: mass  -- mass of resonance
  ! * real :: partialWidth_pole -- partial widht of decayChannel at pole mass
  ! * real :: poleMass -- pole mass of resonance
  ! OUTPUT
  ! * real :: partialWidth_mass -- partial width at real mass
  ! * real :: rho_AB_atPole -- Equation 2.76 of Effe Dr. thesis evaluated at pole mass
  !****************************************************************************
  subroutine rhoDeltaFinalState (mass, poleMass, L, partialWidth_pole, partialWidth_mass, rho_AB_atPole)
    use constants, only: mN, mPi
    use distributions, only: markusPostFormFactor

    real, intent(in) :: mass, poleMass
    integer, intent(in) :: L
    real, intent(in) :: partialWidth_pole
    real, intent(out) :: partialWidth_mass, rho_AB_atPole

    real :: rho_ab_mass

    if (initFlag) call readinput

    rho_AB_atPole = deltarho(poleMass,L,.true.)

    if (mass <= mN + 3*mPi) then
       partialWidth_mass=0.
    else
       rho_AB_mass = deltarho(mass,L,.true.)
       !partialWidth_mass=partialWidth_pole*rho_AB_mass/rho_AB_atPole*monopoleFormFactor(mass,poleMass,0.5)
       partialWidth_mass = partialWidth_pole * rho_AB_mass / rho_AB_atPole
       if (use_cutOff) partialWidth_mass = partialWidth_mass &
                       * (markusPostFormFactor(mass,poleMass,mN+3*mPi,deltaRho_cutoff))**2
    end if

  contains
    !**************************************************************************
    !****if* rhoDeltaFinalState/deltaRho
    ! NAME
    ! function deltaRho (srts, ang, flag) result (integral)
    ! NOTES
    ! This routine calculates the integral of the Delta rho final state
    ! It calculates rho_AB as quoted in equation 2.76 in Effenbergers Dr. thesis
    !**************************************************************************
    function deltaRho (srts, ang, flag) result (integral)
      use IdTable, only: rho, sigmaMeson, Delta
      use particleProperties, only: hadron
      use distributions, only: BlattWeisskopf
      use constants, only: pi, twopi, mN, mPi, hbarc
      use decayChannels, only: get_rhoDelta_is_sigmaDelta

      real, intent(in) :: srts
      integer, intent(in) :: ang ! angularMomentum
      logical, intent(in) :: flag ! flag to switch on the relativistic spectral function
      real :: Integral

      real, parameter :: rint = 1./hbarc
      real, parameter :: dy = pi/30.
      real, parameter :: r = 5.3

      real :: gamres0  ! Width of the rho
      real :: gamres1  ! Width of the delta
      real :: mres0    ! Mass of the rho
      real :: mres1    ! Mass of the delta
      integer :: numd,numr,i,j
      real :: delmin,delmax,rhomin,rhomass,pfinal2,pfinal,q2,q,gamrho,gamdel,qr,sigmamr,sigmamd,rhomax,delmass, &
              ymax1,ymin1,ymax2,ymin2,y1,y2,dy1,dy2,intfac0,intfac1

      if (get_rhoDelta_is_sigmaDelta()) then
        mres0 = hadron(sigmaMeson)%mass
        gamres0 = hadron(sigmaMeson)%width
      else
        mres0 = hadron(rho)%mass
        gamres0 = hadron(rho)%width
      end if

      mres1 = hadron(delta)%mass
      gamres1 = hadron(delta)%width

      delmin=mN+mPi
      delmax=srts-2.*mPi
      ymax1=2.*atan((delmax-mres1)/gamres1*2.)
      ymin1=2.*atan((delmin-mres1)/gamres1*2.)

      numd=max(int((ymax1-ymin1)/dy),1)
      dy1=(ymax1-ymin1)/float(numd)
      integral=0
      do i=1,numd
         y1=ymin1+(float(i)-0.5)*dy1
         delmass=0.5*tan(y1/2.)*gamres1+mres1
         delmass=min(max(delmass,delmin),delmax)

         rhomin=2.*mPi
         rhomax=srts-delmass

         ymax2=2.*atan((rhomax-mres0)/gamres0*2.)
         ymin2=2.*atan((rhomin-mres0)/gamres0*2.)

         numr=max(int((ymax2-ymin2)/dy),1)
         dy2=(ymax2-ymin2)/float(numr)
         do j=1,numr

            y2=ymin2+(float(j)-0.5)*dy2
            rhomass=0.5*tan(y2/2.)*gamres0+mres0
            rhomass=min(max(rhomass,rhomin),rhomax)

            pfinal2=(srts**2-(delmass+rhomass)**2)*(srts**2-(delmass-rhomass)**2)/(4.*srts**2)

            if (pfinal2<0) then
               write(*,*) 'problems in deltarho pfinal2',pfinal2
               stop
            end if

            pfinal=sqrt(pfinal2)

            if (flag) then
               q2=rhomass**2/4.-mPi**2
               q=sqrt(max(q2,0.))
               qr=sqrt(mres0**2/4.-mPi**2)
               gamrho=gamres0*mres0/rhomass*(q/qr)**3*(1.+r**2*qr**2)/(1.+r**2*q**2)

               q2=(delmass**2-(mN+mPi)**2)*(delmass**2-(mN-mPi)**2)/(4.*delmass**2)
               q=sqrt(max(q2,0.))
               qr=sqrt((mres1**2-(mN+mPi)**2)*(mres1**2-(mN-mPi)**2)/(4.*mres1**2))
               gamdel=gamres1*mres1/delmass*q/qr*(BlattWeisskopf(q*rint,1)/BlattWeisskopf(qr*rint,1))**2

               intfac0=gamres0/((rhomass-mres0)**2+gamres0**2/4.)
               intfac1=gamres1/((delmass-mres1)**2+gamres1**2/4.)
               sigmamr=2./pi*rhomass**2*gamrho/((rhomass**2-mres0**2)**2+rhomass**2*gamrho**2)/intfac0
               sigmamd=2./pi*delmass**2*gamdel/((delmass**2-mres1**2)**2+delmass**2*gamdel**2)/intfac1
            else
               sigmamr=1./twopi
               sigmamd=1./twopi
            end if
            integral=integral+pfinal*(BlattWeisskopf(pfinal*rint,ang))**2*sigmamr*sigmamd*dy1*dy2
         end do
      end do

      if (srts_srt_switch) then
         integral=integral/srts**2
      else
         integral=integral/srts
      end if

    end function deltaRho

  end subroutine rhoDeltaFinalState


end module baryonWidthVacuum
