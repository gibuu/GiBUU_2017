!******************************************************************************
!****m*  /threeBodyPhaseSpace
! NAME
! module threeBodyPhaseSpace
! NOTES
! Includes routines which are necessary for the three-body phase-space.
! The subroutine momenta_in_3BodyPS is put into the  module nBodyPhaseSpace.
!******************************************************************************
module threeBodyPhaseSpace

  implicit none
  private

  public :: Integrate_3BodyPS,Integrate_3bodyPS_Resonance

contains


  !****************************************************************************
  !****f*  threeBodyPhaseSpace/Integrate_3bodyPS
  ! NAME
  ! function Integrate_3bodyPS (srts, mass1, mass2, mass3) result (ps)
  ! PURPOSE
  ! Evaluates Integral over three body phase space in vacuum
  !    ps=d\Phi_3 *16*(2*pi)**7       (d\Phi_3 from PDG)
  ! INPUTS
  ! * real :: srts -- sqrt(s) in the problem
  ! * real :: mass1,mass2,mass3 -- masses of the three particles
  ! OUTPUT
  ! * real :: ps = Integral of phase space
  ! NOTES
  ! Formerly known as "bops3".Is now faster than bops3.
  ! The code is so hard to read, since I wanted to make it faster.
  ! This routine is called often, therefore the optimization became important!
  !****************************************************************************
  function Integrate_3bodyPS (srts, mass1, mass2, mass3) result (ps)

    real, intent(in) :: srts, mass1, mass2, mass3
    real :: ps  ! Integral of the 3body-phase-space

    real m12,e3star,e1star !,m13min,m13max
    integer :: i
    integer, parameter :: nm12 = 100
    real :: sqrtM12,sqrtE1,sqrtE3,m1m2Squared,factor

    ps=0.

    !*Int. von  von m12 [(m1+m2)**2,(srts-m3)**2]
    if (srts.le.mass1+mass2+mass3) return

    !  ps=Integral[ 1/srts**2 dm12 dm13]
    !  ps=Integral[ (m13max-m13min)/srts**2 dm12]
    !  ps=SUM_i [( 1./srts**2* (m13max(i)-m13min(i))*((srts-mass3)**2-   m1m2Squared)/float(nm12))] with i=1,...,nm12
    !  ps=SUM_i [(m13max(i)-m13min(i))]        /srts**2* *((srts-mass3)**2-   m1m2Squared)/float(nm12))   with i=1,...,nm12
    !  First we evaluate the sum and then the factor is multiplied

    ! (1) Evaluate Integral( m13max-m13min)
    m1m2Squared=(mass1+mass2)**2
    factor=1./float(nm12)*((srts-mass3)**2- m1m2Squared)
    do i=1,nm12
       !  m12=(float(i)-0.5)/float(nm12)*((srts-mass3)**2- (mass1+mass2)**2)+(mass1+mass2)**2
       m12=(float(i)-0.5)*factor+m1m2Squared

       m12=max(m12,0.)
       sqrtM12=sqrt(m12)
       e3star=(srts**2-m12-mass3**2)/2./sqrtM12
       e1star=(m12+mass1**2-mass2**2)/2./sqrtM12
       sqrtE1=sqrt(max(e1star**2-mass1**2,0.))
       sqrtE3=sqrt(max(e3star**2-mass3**2,0.))

!       m13min=(e1star+e3star)**2-(sqrtE1+sqrtE3)**2
!       m13max=(e1star+e3star)**2-(sqrtE1-sqrtE3)**2
!       m13min=(e1star+e3star)**2
!       m13max=m13min-(sqrtE1-sqrtE3)**2
!       m13min=m13min-(sqrtE1+sqrtE3)**2

       ! ps=ps+ m13max-m13min=ps+4*sqrtE1*sqrtE3
       ps=ps+4*sqrtE1*sqrtE3

    end do
    ! (2) Multiply with factor which does not depend on the loop index :
    ps=ps* 1./srts**2*  ((srts-mass3)**2-   m1m2Squared)/float(nm12)

  end function Integrate_3bodyPS



  !****************************************************************************
  !****f*  threeBodyPhaseSpace/Integrate_3bodyPS_Resonance
  ! NAME
  ! function Integrate_3bodyPS_Resonance (srts, mass1, mass2, resonanceID, scalarPotential) result (ps)
  ! PURPOSE
  ! Evaluates Integral over three body phase space in vacuum with a resonance
  !  among the three particles. Therefore one has to Integrate over the
  ! mass of the resonance as well.
  ! ps=Integral d(massResonance) d\Phi_3 *16*(2*pi)**7 (d\Phi_3 from PDG)
  ! INPUTS
  ! * real :: mass1,mass2 = masses of the three particles
  ! * integer :: idRes = Id of the resonance
  ! * real,optional :: scalarPotential = scalarPotential of the resonance
  ! * srts : sqrt(s) in the problem
  ! OUTPUT
  ! * real,dimension(1:2) :: ps = Integral of phase space
  ! * ps(1) : Full width in the nominator of the spectral function
  ! * ps(2) : Only (nucleon kaonBar) width in nominator of spectral function
  ! NOTES
  ! Formerly known as "massInt"
  !
  ! Be careful : since there is a mass cut off on the masses, there is no
  ! normalization off the spectral functions any more
  !****************************************************************************
  function Integrate_3bodyPS_Resonance (srts, mass1, mass2, resonanceID, scalarPotential) result (ps)
    use idTable, only: nucleon, kaonBar, isMeson
    use particleProperties, only: hadron
    use baryonWidthMedium, only: WidthBaryonMedium, partialWidthBaryonMedium
    use MesonWidthMedium, only: WidthMesonMedium
    use mediumDefinition, only: vacuum
    use constants, only: pi

    real,    intent(in) :: srts, mass1, mass2
    integer, intent(in) :: resonanceID
    real,    intent(in), optional :: scalarPotential
    real, dimension(1:2) :: ps

    real :: minmass,maxmass,psa,mass, gamtot,spectral(1:2) !,dma,dm
    real :: intfac,y,ymax,ymin,mres0,gamres0,dya
    integer :: i,nm
    real,dimension(0:3)  :: momLRF
    logical :: mesonFlag ! If resonance is a meson
    real,  parameter :: dy=2*pi/100.

    momLRF=0.

    ! Standard output:
    ps=0.

    ! Check Input, check wether resonance is meson or baryon
    mesonFlag = isMeson(resonanceID)
    minmass=hadron(resonanceID)%minmass

    ! pole-mass and pole-width
    mres0=hadron(resonanceID)%mass
    gamres0=hadron(resonanceID)%width

    maxmass=srts-mass1-mass2-scalarPotential

    if (gamres0<1e-02) then
       ! Assume resonance to be stable
       ps(1) = Integrate_3bodyPS (srts, mass1, mass2, mres0+scalarPotential)
       ps(2) = 0.
    else if (maxmass<=minmass) then
       ps=0.
    else
       ! integrate over mass distribution
       ps=0.
       ymax=2.*atan((maxmass-mres0) /gamres0*2.)
       ymin=2.*atan((minmass-mres0) /gamres0*2.)
       nm=max(int((ymax-ymin)/dy),1)
       dya=(ymax-ymin)/float(nm)
       do i=1,nm
          y=ymin+(float(i)-0.5)*dya
          mass=.5*tan(y/2.)*gamres0+mres0
          mass=min(max(mass,minmass),maxmass)

          ! Evaluate Gamma(mass+scalarPotential)
          if (mesonFlag) then
             gamtot=WidthMesonMedium(resonanceID,mass,momLRF,vacuum)
          else
             gamtot=WidthBaryonMedium(resonanceID,mass,momLRF,vacuum)
          end if

          spectral(1)=2./pi*mass**2*gamtot/((mass**2- mres0**2)**2+gamtot**2*mass**2)
          if (mesonFlag) then
             spectral(2)=0
          else
             ! spectral function * Gamma( nucleon kaonBar )/Gamma_tot
             spectral(2)=spectral(1)*partialWidthBaryonMedium(resonanceID,mass,.false.,kaonBar,nucleon,momLRF,vacuum)/gamtot
           end if

          psa = Integrate_3bodyPS (srts, mass1, mass2, mass+scalarPotential)

          intfac=gamres0/((mass-mres0)**2+gamres0**2/4.)

          ps(1)=ps(1)+psa*spectral(1)*dya/intfac
          ps(2)=ps(2)+psa*spectral(2)*dya/intfac

       end do
    end if
  end function Integrate_3bodyPS_Resonance

end module threeBodyPhaseSpace
