!******************************************************************************
!****m* /mesonWidthVacuum
! NAME
! module mesonWidthVacuum
! NOTES
! Module which calculates the partial and full widths of the meson resonances
! in dependence of their mass.
! Their pole mass is given by their mass in 'particleProperties' and the
! widths at this pole mass as well.
! Everything corresponds to the vacuum situation. The resulting width is
! therefore always the width in vacuum!
! As mass of the resonance we use the four-vector definition: p_mu p^mu= mass**2
! Prescription according to Manley et al. Phys. Rev. D45 (1992) 4002.
! The In-Width is assumed to be the outwidth.
!******************************************************************************
module mesonWidthVacuum

  use constants, only: hbarc
  implicit none
  private

  logical, parameter :: debug=.false.
  real, parameter :: interactionRadius = 1./hbarc       ! Interaction radius for Blatt Weisskopf functions (1 fm)
  logical, parameter :: srts_srt_switch = .false.       ! Modifies the width according to S. Leupold's definition of the width.
                                                        ! One especially has to exchange s against sqrt(s)
                                                        ! in the denominator of Formula 2.76 of Effenbergers Phd.


  !****************************************************************************
  !****g* mesonWidthVacuum/omega_width
  ! SOURCE
  !
  integer, save :: omega_width = 1
  ! PURPOSE
  ! Select a parametrization for the omega vacuum width:
  ! * 1 = GiBUU default (a la Manley)
  ! * 2 = Muehlich
  !****************************************************************************


  logical, save :: initFlag = .true.


  public :: vacuumWidth, dileptonWidth


contains


  !****************************************************************************
  !****s* mesonWidthVacuum/readInput
  ! NAME
  ! subroutine readInput
  ! PURPOSE
  ! Reads input in jobcard out of namelist "MesonWidthVacuum".
  !****************************************************************************
  subroutine readInput
    use output, only: Write_ReadingInput

    !**************************************************************************
    !****n*  mesonWidthVacuum/MesonWidthVacuum
    ! NAME
    ! NAMELIST /MesonWidthVacuum/
    ! PURPOSE
    ! Includes the input switches:
    ! * omega_width
    !**************************************************************************
    NAMELIST /mesonWidthVacuum/ omega_width

    integer :: ios

    call Write_ReadingInput('mesonWidthVacuum',0)
    rewind(5)
    read(5,nml=mesonWidthVacuum,IOSTAT=ios)
    call Write_ReadingInput('mesonWidthVacuum',0,ios)

    write(*,*) 'omega width parametrization: ', omega_width

    call Write_ReadingInput('mesonWidthVacuum',1)
    initFlag = .false.

  end subroutine readInput


  !****************************************************************************
  !****s* mesonWidthVacuum/vacuumWidth
  ! NAME
  ! function vacuumWidth (mass, partID, ratio) result (gammaTotal)
  ! NOTES
  ! This routine calculates the total vacuum decay width of a mesonic
  ! resonance, i.e. the sum of all partial decay widths, and the branching
  ! ratios for each channel as a function of the offshell mass (m^2 = p_mu p^mu).
  ! Parameters taken from Manley et al. Phys. Rev. D45 (1992) 4002 and PDG.
  ! The mass dependence of the resonances is treated according to Manley.
  ! INPUTS
  ! * real,    intent(in) :: mass     -- offshell mass of Resonance in GeV
  ! * integer, intent(in) :: partID   -- ID of resonance
  ! OUTPUT
  ! * real, dimension(1:nDecays), optional :: ratio -- branching ratio,
  !   i.e. partial width / full width for all decay channels of the mesons
  ! * real :: gammaTotal -- full width (mass)  in GeV at the mass of "mass"
  !****************************************************************************
  function vacuumWidth (mass, partID, ratio) result (gammaTotal)
    use DecayChannels, only: Decay2BodyMeson, Decay3BodyMeson
    use particleProperties, only: hadron, nDecays, get_rho_dilep, getAngularMomentum_meson
    use IdTable, only: photon, dstar, dStarbar, rho, omegaMeson

    real, intent(in) :: mass
    integer, intent(in) :: partID
    real, intent(out), dimension(nDecays), optional :: ratio
    real :: gammaTotal

    real, parameter :: widthCutOff = 1E-11
    real, dimension(nDecays) :: gamma
    real :: partialWidth, poleMass
    real, dimension(1:3) :: massOut
    integer :: L, k, dID

    if (initFlag) call readInput

    ! Initialize
    gammaTotal=0.
    gamma=0.
    if (present(ratio)) ratio=0.

    if ((partID==dstar) .or. (partID==dStarbar)) return

    poleMass=hadron(partID)%mass

    do k=1,nDecays
       dID = hadron(partID)%decaysID(k)
       select case (dID)
       case (0)
          cycle

       case (1:) ! ==== 2Body: >0 ====
          partialWidth=hadron(partID)%decays(k)*hadron(partID)%width !in GeV
          if (partialWidth<widthCutOff) cycle

          massOut = 0.
          if (Decay2BodyMeson(dID)%ID(1) /= photon) massOut(1)=hadron(Decay2BodyMeson(dID)%ID(1))%mass
          if (Decay2BodyMeson(dID)%ID(2) /= photon) massOut(2)=hadron(Decay2BodyMeson(dID)%ID(2))%mass

          select case (dID)
          case (1,3,4,6,12)
             ! in the scalar-scalar channels use the mass dependent width
             ! Set Angular Momentum of decaying particles
             L = getAngularMomentum_meson(dID,partID)
             if (debug) write(*,*) "Masses", massOut(1:2), polemass, "angular=",L
             if (partID==omegaMeson .and. omega_width == 2) then
                if (dID==1) then
                   gamma(k) = omegaVacuum (mass**2, 1)   ! omega -> 2pi
                else if (dID==6) then
                   gamma(k) = omegaVacuum (mass**2, 2)   ! omega -> pi0 gamma
                end if
             else
                ! decay into 2 stable particles
                gamma(k) = stableFinalState (mass, poleMass, massOut(1),massOut(2), L, partialWidth)
             end if

          case (2)
             ! -> rho pi
             L = getAngularMomentum_meson(dID,partID)
             if (debug) write(*,*) "Masses", massOut(1:2), polemass, "angular=",L
             gamma(k) = semiStableFinalState(mass,poleMass,Decay2BodyMeson(dID)%ID(1:2),L,partialWidth)

          case default
             ! in other channels assume width to be constant in mass
             if (mass>sum(massOut(1:2))) gamma(k)=partialWidth
          end select

          if (debug) write(*,*) "Stable, decayChannel:",dID,"Gamma=",gamma(k)

       case (:-1) ! ==== 3Body: <0 ====
          partialWidth=hadron(partID)%decays(k)*hadron(partID)%width !in GeV
          if (partialWidth<widthCutOff) cycle

          massOut = 0.
          if (Decay3BodyMeson(-dID)%ID(1) /= photon) massOut(1)=hadron(Decay3BodyMeson(-dID)%ID(1))%mass
          if (Decay3BodyMeson(-dID)%ID(2) /= photon) massOut(2)=hadron(Decay3BodyMeson(-dID)%ID(2))%mass
          if (Decay3BodyMeson(-dID)%ID(3) /= photon) massOut(3)=hadron(Decay3BodyMeson(-dID)%ID(3))%mass

          if (mass < sum(massOut(1:3))) cycle

          select case (-dID)
          case (2,3)
            if (partID==omegaMeson .and. omega_width==2 .and. dID==-2) then
               gamma(k) = omegaVacuum (mass**2, 3)  ! omega -> 3pi
            else
               gamma(k) = threePi (partID, partialWidth, mass)
            end if
          case default
            gamma(k)=partialWidth
          end select

       end select
    end do

    gammaTotal=Sum(gamma(:))

    if (partID == rho .and. get_rho_dilep()) gammaTotal = gammaTotal + dileptonWidth (partID, mass)

    if (debug) write(*,*) gammaTotal

    if (present(ratio)) then
      if (gammaTotal > 0) then
        ratio(:)=gamma(:)/gammaTotal
      else
        ratio(:)=0.
      end if
    end if

    if (debug) write(*,*)

  end function vacuumWidth


  !****************************************************************************
  !****f* mesonWidthVacuum/stableFinalState
  ! NAME
  ! real function stableFinalState(mass, poleMass, mass1, mass2, L, partialWidth_pole)
  ! NOTES
  ! Resonance decays into stable particles.
  ! Decay only allowed if:
  ! * (1) decay width > widthCutOff and ...
  ! * (2) mass of resonance > Sum of masses of decay products
  !
  ! Manley et al. Phys. Rev. D45 (1992) 4002.
  ! OUTPUT
  ! * returns the partial width for the specified channel at the given mass
  !****************************************************************************
  real function stableFinalState(mass, poleMass, mass1, mass2, L, partialWidth_pole)
    use distributions, only: BlattWeisskopf
    use twoBodyTools, only: pCM

    integer, intent(in) :: L !Angular Momentum
    real, intent(in) :: polemass  ! pole mass of mother resonance
    real, intent(in) :: mass  !mass of  mother resonance
    real, intent(in) :: mass1 !mass of first decay product
    real, intent(in) :: mass2 !mass of second decay product
    real, intent(in) :: partialWidth_pole ! partial widht of decayChannel at pole mass
    real :: p_ab_mass,p_ab_pole,rho_ab_mass,rho_ab_pole

    if (mass <= mass1+mass2) then
      ! decay not allowed
      stableFinalState = 0.
      return
    end if

    ! Determine momentum of outgoing particles in Restframe of Resonance
    p_ab_mass=pCM(mass,mass1,mass2)
    p_ab_pole=pCM(polemass,mass1,mass2)

    ! Evaluate rho_ab according to equation 2.76 in Effenbergers Dr. thesis
    ! rho_ab(mu)=p_ab/mu * BlattWeisskopf(pab*interactionRadius,L)
    if (srts_srt_switch) then
       rho_ab_mass=p_ab_mass/mass**2*(BlattWeisskopf(p_ab_mass*interactionRadius,L))**2
       rho_ab_pole=p_ab_pole/polemass**2*(BlattWeisskopf(p_ab_pole*interactionRadius,L))**2
    else
       rho_ab_mass=p_ab_mass/mass*(BlattWeisskopf(p_ab_mass*interactionRadius,L))**2
       rho_ab_pole=p_ab_pole/polemass*(BlattWeisskopf(p_ab_pole*interactionRadius,L))**2
    end if

    stableFinalState = partialWidth_pole * rho_ab_mass/rho_ab_pole
  end function stableFinalState


  !****************************************************************************
  !****f* mesonWidthVacuum/semistableFinalState
  ! NAME
  ! real function semistableFinalState(mass,polemass,ID,L,partialWidth_pole)
  ! NOTES
  ! Resonance decays into one stable and one unstable particle.
  ! Calculates the partial width dependend on mass of a meson resonance
  ! decaying into one stable and one unstable decay product.
  ! According to Manley and Effenberger' Dr. Thesis equation 2.76.
  ! OUTPUT
  ! * returns the partial width for the specified channel at the given mass
  !****************************************************************************
  real function semistableFinalState(mass,polemass,ID,L,partialWidth_pole)

    use rhoABIntegrandVariables, only: L_copy, idUnstable_copy, massStable, massUnStable, gammaUnStable, srts
    use quadpack, only: qag
    use particleProperties, only: hadron

    real, intent(in) :: mass       ! mass of mother resonance
    real, intent(in) :: polemass   ! polemass of mother resonance
    integer, intent(in) :: ID(1:2) ! outgoing IDs: 1=stable, 2=unstable
    Integer, intent(in) :: L       ! angular momentum of final state
    real, intent(in) :: partialWidth_pole ! partial width at pole of mother resonance

    real :: LowerBound, upperBound ! lower and upper bound of integrals
    real :: rho_AB_pole, rho_AB_mass ! Equation 2.76 Effenberger Dr.

    logical,parameter :: debug=.false.
    real,   parameter :: qag_absError=0., qag_relError=0.00001
    integer,parameter :: qag_key=6
    real              :: qag_error_estimation
    integer           :: qag_neval,qag_error_code

    ! setting for rho_ab_integrand
    L_copy=L
    idUnstable_copy=ID(2)

    massStable=hadron(ID(1))%mass
    massUnStable=hadron(ID(2))%mass
    gammaUnStable=hadron(ID(2))%width

    if (debug) write(*,'(A,3F8.4)') "In semiStable:", massStable, massUnstable, gammaUnstable

    lowerbound = hadron(ID(2))%minmass
    if (mass <= lowerBound + massStable) then
      ! decay not allowed
      semistableFinalState = 0.
      return
    end if

    ! (1) Use kinematics at the actual mass of the resonance :
    ! Upperbound is the maximal mass of the unstable particle in the decay product.
    ! lowerbound is the minimal mass of the unstable particle : therefore the restmasses of its stable decay products
    srts=mass
    upperBound=mass-massStable
    call qag (rho_ab_Integrand, lowerBound, upperBound, qag_absError, qag_relError, qag_key, rho_ab_mass, &
              qag_error_estimation, qag_neval, qag_error_code)

    if (qag_error_code /= 0) then
      write(*,*) 'Error with QAG in mesonWidthVacuum (mass):', qag_error_code, qag_error_estimation, qag_neval
    end if

    ! (2) Use kinematics at the pole mass of the resonance :
    ! Upperbound is the maximal mass of the unstable particle in the decay product.
    srts=polemass
    upperBound=poleMass-massStable
    call qag (rho_ab_Integrand, lowerBound, upperBound, qag_absError, qag_relError, qag_key, rho_ab_pole, &
              qag_error_estimation, qag_neval, qag_error_code)

    if (qag_error_code /= 0) then
      write(*,*) 'Error with QAG in mesonWidthVacuum (pole):',qag_error_code, qag_error_estimation, qag_neval
    end if

    ! (3) Calculate partial width
    if (debug) then
       write(*,*) "Rho_AB_pole=", rho_ab_pole
       write(*,*) "Rho_AB_mass=", rho_ab_mass
    end if

    if (rho_ab_pole == 0.) then
       semistableFinalState = 0.
    else
       if (debug) write(*,*) partialWidth_pole,rho_ab_mass,rho_ab_pole
       semistableFinalState = partialWidth_pole * rho_ab_mass/rho_ab_pole
    end if

  end function semistableFinalState


  !****************************************************************************
  !****f* mesonWidthVacuum/rho_AB_Integrand
  ! NAME
  ! real function rho_AB_Integrand (mu)
  ! NOTES
  ! This is the integrand of rho_ab(mu) in Effenbergers Dr.-Thesis, page 27, formula 2.76.
  ! Only in the case of only one unstable particle, currently only used for decays into "rho pi"!
  !****************************************************************************
  real function rho_AB_Integrand (mu)

    use IdTable, only: rho
    use rhoABIntegrandVariables, only: idUnstable_copy, L_copy, massStable, massUnStable, gammaUnStable, srts
    use twoBodyTools, only: pCM_sqr
    use distributions, only: BlattWeisskopf,RelBW
    use constants, only: mPi

    real, intent(in) :: mu ! Mass of unstable particle
    real :: p_AB, gamma

    ! "srtS" is the mass of mother Resonance. Set in the mother subroutine.

    ! It decays into a stable particle with mass "massStable" and
    ! an unstable particle with mass "mu".
    ! p_AB is the momentum of the two decay products in the restframe of the mother resonance.

    ! massStable, massUnstable and gammaUnstable are all defined in the mother subroutine.

    p_AB = pCM_sqr(srtS**2, mu**2, massStable**2)
    if (p_AB < 0.) then
      ! this situation is kinematically forbidden
      rho_AB_Integrand = 0.
      return
    end if
    p_AB = sqrt (p_AB)

    ! Redefine the gamma of the unstable particle, since the unstable particle is supposed off-shell:

    if (idUnstable_copy==rho) then
       ! then the rho decays into two pions
       ! this rho has mass mu, so we need gammaUnstable(mu)
       if (mu>2*mPi) then
         gamma = stableFinalState (mu, massUnstable, mPi, mPi, 1, gammaUnstable)
       else
         gamma = 0.
       end if
    else
       write(*,*) "Wrong Resonance as unstable Particle in rho_ab_integrand :",idUnstable_copy
       write(*,*) "This resonance is not yet implemented in subroutine 'semistableFinalState' :",idUnstable_copy
       Stop
    end if

    ! Define gamma as function of mu
    if (srts_srt_switch) then
       rho_AB_Integrand = RelBW(mu,massUnStable,gamma) * p_AB/srtS**2 * BlattWeisskopf(p_AB*interactionRadius,L_copy)**2
    else
       rho_AB_Integrand = RelBW(mu,massUnStable,gamma) * p_AB/srtS    * BlattWeisskopf(p_AB*interactionRadius,L_copy)**2
    end if

    if (debug) then
      write(*,*) RelBW(mu,massUnStable,gamma),p_AB,srtS,BlattWeisskopf(p_AB*interactionRadius,L_copy)**2
      write(*,*) "rho_AB_Integrand=",rho_AB_Integrand
    end if
  end function rho_AB_Integrand


  !****************************************************************************
  ! 3 pi decay width
  ! according to PDG 2006, equ. (38.19), assuming constant matrix element
  ! and isotropic angles
  !****************************************************************************
  real function threePi(ID,ppwidth,srts)

    use particleProperties, only: hadron
    use gauss_integration, only: sg20r, rg20r
    use twoBodyTools, only: pCM
    use constants, only: mPi

    integer,intent(in):: ID       ! ID of decaying particle
    real,intent(in)   :: ppwidth  ! "partial pole width"
    real,intent(in)   :: srts     ! sqrt(s) = mass of decaying particle
    integer::i,ns
    integer,parameter::n=3
    real,dimension(1:20*n)::absi,orde
    real :: m,M0,resu
    real, save :: norm
    real, parameter :: L = 0.5

    M0 = hadron(ID)%mass

    if (srts>3.*mPi) then
       ! on-shell (has to be determined only once)
       call sg20r(2.*mPi,M0-mPi,n,absi,ns)
       do i=1,ns
          m=absi(i)
          orde(i)=2.*m*pCM(M0,mPi,m)*pCM(m,mPi,mPi)/M0**2
       end do
       call rg20r(2.*mPi,M0-mPi,n,orde,norm)
       ! at actual mass sqrt(s)
       call sg20r(2.*mPi,srts-mPi,n,absi,ns)
       do i=1,ns
          m=absi(i)
          orde(i)=2.*m*pCM(srts,mPi,m)*pCM(m,mPi,mPi)/srts**2
       end do
       call rg20r(2.*mPi,srts-mPi,n,orde,resu)
       ! multiply with partial on-shell width and cut-off factor
       threePi = ppwidth * resu/norm * (L**2+M0**2)/(L**2+srts**2)
    else
       threePi=0.
    end if
  end function threePi


  !****************************************************************************
  !****f* mesonWidthVacuum/omegaVacuum
  ! NAME
  ! real function omegaVacuum (s, r)
  ! INPUTS
  ! real, intent(in)    :: s     ! m^2 of the omega meson
  ! integer, intent(in) :: r     ! select channel
  ! NOTES
  ! Returns the partial omega widths in the vacuum:
  ! * r=0 : total width
  ! * r=1 : omega -> pi pi
  ! * r=2 : omega -> pi0 gamma
  ! * r=3 : omega -> rho pi (as an approximation for omega -> 3pi)
  ! See: Pascal Muehlich, PhD thesis, chapter 8.2.1.
  !****************************************************************************
  real function omegaVacuum (s, r)
    use constants, only: pi, mPi
    use IDTable, only: rho, omegaMeson
    use particleProperties, only: hadron
    use twoBodyTools, only: pCM

    real,    intent(in) :: s
    integer, intent(in) :: r
    ! branching ratios:
    real,parameter::rpipi=0.017
    real,parameter::rpigam=0.087
    real,parameter::rrhopi=0.896
    ! parameters for rho pi decay:
    real,parameter::frho=6.14
    !    real,parameter::gomropi=2.165

    real :: mome

    mome = hadron(omegaMeson)%mass

    select case (r)
    case (1) ! pi pi
       omegaVacuum=gpipi(s)
    case (2) ! pi0 gamma
       omegaVacuum=gpigam(s)
    case (3) ! rho pi
       omegaVacuum=grhopi(s)   ! min(grhopi(s),hadron(omegaMeson)%width)
    case (0) ! total width
       omegaVacuum=gpipi(s)+gpigam(s)+grhopi(s)
    end select

  contains

    !
    ! omega -> pi pi width
    !
    real function gpipi(s)
      real,intent(in)::s
      if (sqrt(s)>2.*mPi) then
         gpipi=hadron(omegaMeson)%width*rpipi*(pCM(sqrt(s),mPi,mPi)/pCM(mome,mPi,mPi))**3
      else
         gpipi=0.
      end if
    end function gpipi

    !
    ! omega -> rho pi width
    !
    real function grhopi(s)
      use gauss_integration, only: sg20r, rg20r
      real,intent(in)::s
      integer::i,ns
      integer,parameter::n=3
      real,dimension(1:20*n)::absi,orde
      real::m,resu,grpi
      logical,save::first=.true.
      real,save::norm
      if (sqrt(s)>3.*mPi) then
         if (first) then
            ! on-shell (has to be determined only once)
            call sg20r(2.*mPi,mome-mPi,n,absi,ns)
            do i=1,ns
               m=absi(i)
               grpi=(pCM(mome,mPi,m))**3
               orde(i)=2.*m*arho(m)*grpi
            end do
            call rg20r(2.*mPi,mome-mPi,n,orde,norm)
            first=.false.
         end if
         ! at actual omega mass sqrt(s)
         call sg20r(2.*mPi,sqrt(s)-mPi,n,absi,ns)
         do i=1,ns
            m=absi(i)
            grpi=(pCM(sqrt(s),mPi,m))**3
            orde(i)=2.*m*arho(m)*grpi
         end do
         call rg20r(2.*mPi,sqrt(s)-mPi,n,orde,resu)
         ! multiply with partial on-shell width
         grhopi=hadron(omegaMeson)%width*rrhopi*resu/norm
      else
         grhopi=0.
      end if
    end function grhopi

    !
    ! omega -> pi0 gamma width
    !
    real function gpigam(s)
      real,intent(in)::s
      if (sqrt(s)>mPi) then
         gpigam=hadron(omegaMeson)%width*rpigam*((s-mPi**2)/sqrt(s))**3/((mome**2-mPi**2)/mome)**3
      else
         gpigam=0.
      end if
    end function gpigam

    !
    ! rho -> pi pi width
    !
    real function gamrho(s)
      real,intent(in)::s
      if (s>(2.*mPi)**2) then
         gamrho=1./(6.*pi)*frho**2*pCM(sqrt(s),mPi,mPi)**3/s
      else
         gamrho=0.
      end if
    end function gamrho

    !
    ! rho vacuum spectral function
    !
    real function arho(m)
      real,intent(in)::m
      real::mrho
      mrho = hadron(rho)%mass
      arho=1./pi*mrho*gamrho(m**2)/((m**2-mrho**2)**2+(mrho*gamrho(m**2))**2)
    end function arho

  end function omegaVacuum


  !****************************************************************************
  !****s* mesonWidthVacuum/dileptonWidth
  ! NAME
  ! real function dileptonWidth (ID, mass)
  ! NOTES
  ! Calculate the dilepton decay width of the vector mesons rho, omega and phi (V -> e+e-)
  ! according to strict vector-meson dominance (VMD).
  ! References:
  !  * http://arxiv.org/abs/1203.3557v2, eq. (11) and table 3.
  !  * Effenberger PhD, chapter 2.7.1
  ! INPUTS
  ! * integer, intent(in) :: ID   --- vector meson ID (should be 103, 105 or 107)
  ! * real, intent(in) :: mass    --- vector meson mass in GeV
  ! OUTPUT
  ! * width in GeV
  !****************************************************************************
  real function dileptonWidth (ID, mass)
    use particleProperties, only: hadron
    integer, intent(in) :: ID
    real, intent(in) :: mass
    real, parameter :: CV(3) = (/ 9.078E-06, 7.666E-07, 1.234E-06 /)
    dileptonWidth = CV((ID-101)/2)*hadron(id)%mass**4/mass**3
  end function dileptonWidth


end module mesonWidthVacuum
