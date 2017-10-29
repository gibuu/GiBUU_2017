!******************************************************************************
!****m* /parBarMes_HighEnergy
! NAME
! module parBarMes_HighEnergy
! PURPOSE
! Includes all routines which are parametrizations of
! "baryon meson -> X" data for the high-energy region.
!******************************************************************************
module parBarMes_HighEnergy

  implicit none
  private

  public :: paramBarMesHE
  public :: paramBarMesHE_v
  public :: paramBarMesHE_pion
!  PUBLIC :: paramBarMesHE_kaon

contains


  !****************************************************************************
  !****s*  parBarMes_HighEnergy/paramBarMesHE
  ! NAME
  ! subroutine paramBarMesHE (srts, idM, idB, izM, izB, media, sigma, sigmael)
  !
  ! PURPOSE
  ! Calculates elastic and total high-energy meson-baryon cross sections.
  !
  ! INPUTS
  ! * real :: srts -- SQRT(s) (free)
  ! * integer :: idM,idB -- ID of meson and baryon
  ! * integer :: izM,izB -- charge of meson and baryon
  ! * type(medium) :: media  -- medium at position
  ! OUTPUT
  ! * real :: sigma -- sigma(m B -> X) (in mb)
  ! * real :: sigel -- sigma(m B -> m B) (in mb)
  !
  ! NOTES
  ! * Expects that anti-baryons have (-) sign in the ID (e.g. antinucleon=-1).
  !   Expects that mesons are no anti-mesons!!
  ! * The elastic cross section of V N -> V N is derived from the
  !   gamma N -> V N cross section.
  ! * Here we call corresponding single standing routines and perform
  !   some modifications in order to guarantee a smoother transition
  !   from low to high energy.
  ! * This routine was formerly known as "fritziCS".
  !****************************************************************************
  subroutine paramBarMesHE (srts, idM, idB, izM, izB, media, sigma, sigmael)
    use IdTable
    use mediumDefinition
    use photonXSections, only: calcXS_gammaN2VN
    use particleProperties, only: hadron
    use constants, only: alphaQED, mN, mK
    use output, only: DoPR
    use parametrizationsBarMes, only: sigelkp

    real,         intent(in)  :: srts
    integer,      intent(in)  :: idM,idB,izM,izB
    type(medium), intent(in)  :: media
    real,         intent(out) :: sigma,sigmael

    integer :: iidM,ivec, izB2,izM2

    integer, dimension(103:109), parameter :: iivec = (/1,0,2,0,3,0,4/) ! rho,omega,phi or J/psi?

    real, dimension(4) :: sigVM1
    real, parameter, dimension(1:4) :: gv = (/2.2, 23.6,18.4,11.5/) !  VMD coupling constants

    real, dimension(-2:2) :: sigK1, sigK2
    real, dimension(-1:1) :: sigP1, sigP2

    real :: xs,xc ! scaling of XS according strangeness/charm content
    real :: E,plab

    !==== reset all values:

    sigmael=0.
    sigma=0.

    ! select whether rho is rho0 or not...
    iidM = idM
    if (iidM==rho .and. izM/=0) iidM = 0


    select case (iidM)

    case (rho,omegaMeson,phi,JPsi) !==== Vector mesons ====

       if (abs(idB)==nucleon) then
          ivec=iivec(iidM) ! rho0,omega,phi or J/psi?

          !---- Elastic cross section via photon cross section : ----

          call calcXS_gammaN2VN(srts,media,sigVM1)
          sigmael = sigVM1(ivec)/alphaQED*gv(ivec)/1000.

          ! connection of Hi and Lo cross sections / "Flunschen":

          select case (ivec)
          case (1) ! == rho0 ==
             if (srts<=5.0) sigmael = min(3.15,sigmael)
          case (2) ! == omega ==
             if (srts<=4.1) sigmael = 25.*exp(-srts+0.8)+3.1
          end select

          !---- Total cross section: ----

          call paramBarMesHE_v(srts,sigVM1)
          sigma = sigVM1(ivec)
       end if

    case (kaon:kaonStarBar) !==== Kaons ====

       call paramBarMesHE_kaon(srts,sigK1,sigK2)
       izB2 = izB * sign(1,idB)
       izM2 = 0

       select case (idM)
       case (kaon,kaonStar)
          izM2 = -2
          if ((izB2<=0 .and. izM==1) .or. (izB2>0 .and. izM==0)) izM2 = -1

          if (abs(idB)==nucleon .and. idM==kaon) then
             if (srts<=2.614 .and. idB>0) then
                E=(srts**2-mK**2-mN**2)/(2.*mN)
                plab=sqrt(E**2-mK**2)
                if (izM+izB==1) then
                   sigmael=sigelkp(plab,2)
                else
                   sigmael=sigelkp(plab,1)
                end if
             else if (srts<=4.8) then
                sigmael=1/(srts-1.98)+2.8
             else
                sigmael=sigK2(sign(2,-idB)) ! -2 für baryon, +2 für anti-baryon
             end if
          end if

       case (kaonBar,kaonStarBar)
          izM2 = 2
          if ((izB2<=0 .and. izM==-1) .or. (izB2>0 .and. izM==0)) izM2 = 1

          if (abs(idB)==nucleon .and. idM==kaonBar) then
             if (izM2*sign(1,idB)==2 .and. srts<=2.85) then
                sigmael=20*exp(-srts+1.32)
             else if (izM2*sign(1,idB)==1 .and. srts<=4.469) then
                ! sigmael= 20*exp(-srts)+2.5
                sigmael=3.31  ! this works better for the high-momentum asymtotic
                              ! of the data on K- n --> K- n cross section.
                              ! There is still a jump between this value
                              ! and low-energy parameterization. To be improved.
             else
                sigmael=sigK2(sign(2,idB))
             end if
          end if
       end select

       izM2=izM2*sign(1,idB)

!       if (idM.eq.kaonBar.and.izM2.eq.1.and.srts.le.3.) then
!          sigma=20*exp(-srts+2)+14.5

       if (idM==kaonStarBar .and. izM2==2 .and. srts<=4.4) then
          sigma = max (0., 20.*(1.-exp(-2*(srts-2.06)))+3.)
       else if (idM==kaonStarBar .and. izM2==1 .and. srts<=5.95) then
          sigma = max (0., 20.*(1.-exp(-2*(srts-2.18))))
       else
          sigma=sigK1(izM2)
       end if

    case default !==== all other mesons/ nonstrange baryons====

       call paramBarMesHE_pion(srts,sigP1,sigP2)

       izM2 = min(max((2*izB-1)*izM,-1),1)

       xs=0.
       xc=0.

       if (abs(idB)==nucleon .and. idM==pion) then
          sigmael=sigP2(izM2) ! elastic pi N scattering
       else
          select case (idM)
          case (etaPrime)
             xs = 1.0
          case (etaC)
             xc = 1.0
          case (dMeson:dStarBar)
             xc = 0.5
          case (dS_plus:dSStar_Minus)
             xs = 0.5
             xc = 0.5
          end select
       end if

       sigma = sigP1(izM2) * (1.0-0.4*xs-0.5*xc)
       ! Cross sections are scaled according to Falter Phd page 72

    end select


    ! First Check:
    if (sigma<0) then
      if (DoPR(2)) write(*,'(A,f8.4,4i5,2f12.5,2i5)') 'PROBLEMS in paramBarMesHE: sigma<0: ', &
                   srts,idM,izM,idB,izB,sigma,sigmael,izM2,izB2
      sigma=0
   end if

    ! Now set elastic cross section, if not yet set:
    if (sigmael<1e-8) then
       sigmael = 0.039*sigma**(3./2.)

       select case (idM)
       case (eta)
          if (srts<=3.8) sigmael = 50*exp(-2.5*srts+4.6)+5.0

       case (rho)
          izM2=min(max((2*abs(izB)-1)*izM,-1),1)*sign(1,idB)
          select case (izM2)
          case (1)
             if (srts<=3.3) sigmael=min(5.37,sigmael)
          case (-1)
             if (srts<=6.06) sigmael=1.0-exp(-(srts-2.5))+4.0
          end select
       end select

    end if

    ! Scale cross section according strangeness/charm content of baryon:
    if (abs(idB)<Lambda) then
       xs=0.
       xc=0.
    else if (abs(idB)<Xi) then
       xs=1./3.
       xc=0.
    else if (abs(idB)<omegaResonance) then
       xs=2./3.
       xc=0.
    else if (abs(idB)<Lambda_CPlus) then
       xs=1.
       xc=0.
    else
       xc=1./3.
       xs=abs(hadron(abs(idB))%strangeness)/3.
    end if

    sigma   = sigma  *(1.-0.4*xs-0.5*xc)
    sigmael = sigmael*(1.-0.4*xs-0.5*xc)**(3./2.)

    sigma = max (sigma, sigmael)

  end subroutine paramBarMesHE


  !****************************************************************************
  !****s* parBarMes_HighEnergy/paramBarMesHE_v
  ! NAME
  ! subroutine paramBarMesHE_v(srts,sigma)
  ! PURPOSE
  ! Pythia parametrizations of the total vector meson nucleon cross sections
  ! INPUTS
  ! * real :: srts -- SQRT(s) (free)
  ! OUTPUT
  ! * real, dimension(4) :: sigma -- sigma(V N -> X) (in mb)
  ! NOTES
  ! These are pomeron/reggeon parametrizations for V N -> X cross section
  ! according Donnachie and Landshoff for V = rho0, omega, phi, J/Psi.
  ! The parameters are taken from PYTHIA.
  !
  ! Thresholds are not considered.
  !
  ! This was formerly known as "vectot"
  !****************************************************************************
  subroutine paramBarMesHE_v(srts,sigma)

    real, intent(in) :: srts
    real, dimension(4), intent(out) :: sigma

    real :: s

    real, parameter :: eps=0.0808, eta=-0.4525
    real, dimension(4), parameter :: X = (/13.63, 13.63, 10.01, 1.001/)
    real, dimension(4), parameter :: Y = (/31.79, 31.79, -1.52, -0.152/)

    s = srts**2
    sigma = X * s**eps + Y * s**eta

  end subroutine paramBarMesHE_v



  !****************************************************************************
  !****s* parBarMes_HighEnergy/paramBarMesHE_pion
  ! NAME
  ! subroutine paramBarMesHE_pion(srts,sigma,sigel)
  ! INPUTS
  ! * real :: srts -- SQRT(s) (free)
  ! OUTPUT
  ! * real, dimension(-1:1) :: sigma -- sigma(pi N -> X) (in mb)
  ! * real, dimension(-1:1) :: sigel -- sigma(pi N -> pi N) (in mb)
  ! The indices correspond to the charge states:
  ! * -1 -> pi-  p
  ! * 0  -> ( pi+ p + pi- p)/2
  ! * +1 -> pi p
  ! NOTES
  ! This are more or less brute force fits to experimental data.
  !
  ! Reference ???
  !
  ! This was formerly known as "pitot"
  !****************************************************************************
  subroutine paramBarMesHE_pion(srts,sigma,sigel)
    use constants, only: mN, mPi

    real, intent(in) :: srts
    real, dimension(-1:1), intent(out) :: sigma,sigel


    integer i,j,ich,nch
    parameter(nch=2)
    real a(nch,5),plab,a2(nch,5)

    data((a(i,j),j=1,5),i=1,nch)     /16.4,  19.3, -0.42, 0.19, 0.,    33.0,  14.0, -1.36, 0.456, -4.03/
    data((a2(i,j),j=1,5),i=1,nch)      /0.,    11.4, -0.4,  0.079, 0.,    1.76,   11.2, -0.64, 0.043, 0./

    sigma = 0
    sigel = 0

    plab=sqrt((srts**2-mN**2-mPi**2)**2/(4.*mN**2)-mPi**2)

    !==== pi+ p ====
    ich=1
    sigma(1)=min( a(ich,1)+ a(ich,2)*plab**a(ich,3)+         a(ich,4)*(log(plab))**2+         a(ich,5)*log(plab),30.)
    sigel(1)=min(a2(ich,1)+a2(ich,2)*plab**a2(ich,3)+       a2(ich,4)*(log(plab))**2+        a2(ich,5)*log(plab),sigma(1))

    ! somehow follow the structure, but do not overshoot:
    if (plab < 1.45) then
       sigma(1)=0.0 ! 15.0 + (plab-0.75)*(35.0-15.0)/(1.45-0.75)
    else if (plab < 1.85) then
       sigma(1)=35.0 + (plab-1.45)*(30.0-35.0)/(1.85-1.45)
    end if
    if (plab < 1.45) then
       sigel(1)=0.0 ! 8.5 + (plab-0.85)*(16.0- 8.5)/(1.45-0.85)
    else if (plab < 2.00) then
       sigel(1)=16.0 + (plab-1.45)*( 8.5-16.0)/(2.00-1.45)
    end if

    !==== pi- p ====
    ich=2
    sigma(-1)=min( a(ich,1)+ a(ich,2)*plab**a(ich,3)+        a(ich,4)*(log(plab))**2+      a(ich,5)*log(plab),35.)
    sigel(-1)=min(a2(ich,1)+a2(ich,2)*plab**a2(ich,3)+      a2(ich,4)*(log(plab))**2+     a2(ich,5)*log(plab),sigma(-1))

    sigma(0)=(sigma(1)+sigma(-1))/2.
    sigel(0)=(sigel(1)+sigel(-1))/2.

  end subroutine paramBarMesHE_pion


  !****************************************************************************
  !****s* parBarMes_HighEnergy/paramBarMesHE_kaon
  ! NAME
  ! subroutine paramBarMesHE_kaon (srts, sigma, sigel)
  ! INPUTS
  ! * real :: srts -- SQRT(s) (free)
  ! OUTPUT
  ! * real, dimension(-2:2) :: sigma -- sigma(K N -> X) (in mb)
  ! * real, dimension(-2:2) :: sigel -- sigma(K N -> K N) (in mb)
  ! The indices correspond to the charge states:
  ! * -2 -> K+ p, K0  n, K- p~, K0~ n~
  ! * -1 -> K+ n, K0  p, K- n~, K0~ p~
  ! * +1 -> K- n, K0~ p, K+ n~, K0  p~
  ! * +2 -> K- p, K0~ n, K+ p~, K0  n~
  !
  ! NOTES
  ! These are more or less brute force fits to experimental data.
  ! Taken from PDG review, L. Montanet et al., PRD 50, 1173 (1994).
  ! This was formerly known as "kaontot".
  !****************************************************************************
  subroutine paramBarMesHE_kaon (srts, sigma, sigel)
    use constants, only: mN, mK

    real, intent(in) :: srts
    real, dimension(-2:2), intent(out) :: sigma, sigel

    integer :: i,j
    real :: plab
    real, dimension(4,5) :: a1,a2

    ! Parameters from Falter PhD Table 4.4, with Typo for K+ n, corrected from -1.3 to -0.89

    ! total
    data((a1(i,j),j=1,5),i=1,4) /18.1,  0., 0., 0.26, -1.,    &
                                 18.7,  0., 0., 0.21, -0.89,  &
                                 25.2,  0., 0., 0.38, -2.9,   &
                                 32.1,  0., 0., 0.66, -5.6/
    ! elastic
    data((a2(i,j),j=1,5),i=1,4) / 5.,    8.1, -1.8,  0.16, -1.3,  &
                                  5.,    8.1, -1.8,  0.16, -1.3,  &
                                  7.3,    0.,   0.,  0.29, -2.4,  &
                                  7.3,    0.,   0.,  0.29, -2.4 /

    sigma = 0
    sigel = 0

    plab=sqrt((srts**2-mN**2-mK**2)**2/(4.*mN**2)-   mK**2)

    !==== K+ p
    sigma(-2) = a1(1,1) + a1(1,2)*plab**a1(1,3) + a1(1,4)*(log(plab))**2 + a1(1,5)*log(plab)
    sigel(-2) = a2(1,1) + a2(1,2)*plab**a2(1,3) + a2(1,4)*(log(plab))**2 + a2(1,5)*log(plab)

    !==== K+ n
    sigma(-1) = a1(2,1) + a1(2,2)*plab**a1(2,3) + a1(2,4)*(log(plab))**2 + a1(2,5)*log(plab)
    sigel(-1) = a2(2,1) + a2(2,2)*plab**a2(2,3) + a2(2,4)*(log(plab))**2 + a2(2,5)*log(plab)

    !==== K- n
    sigma(1) = a1(3,1) + a1(3,2)*plab**a1(3,3) + a1(3,4)*(log(plab))**2 + a1(3,5)*log(plab)
    sigel(1) = a2(3,1) + a2(3,2)*plab**a2(3,3) + a2(3,4)*(log(plab))**2 + a2(3,5)*log(plab)

    !==== K- p
    sigma(2) = a1(4,1) + a1(4,2)*plab**a1(4,3) + a1(4,4)*(log(plab))**2 + a1(4,5)*log(plab)
    sigel(2) = a2(4,1) + a2(4,2)*plab**a2(4,3) + a2(4,4)*(log(plab))**2 + a2(4,5)*log(plab)

    sigma(:) = min (sigma(:), 30.)
    sigel(:) = min (sigel(:), sigma(:))

  end subroutine paramBarMesHE_kaon


end module parBarMes_HighEnergy
