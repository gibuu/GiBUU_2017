!******************************************************************************
!****m* /barAntiBar
! PURPOSE
! Implements baryon+antibaryon -> X cross section
! NOTES
!******************************************************************************
module barAntiBar

  implicit none
  private

  !****************************************************************************
  !****g* barAntiBar/fact_LambdaBar
  ! SOURCE
  !
  real, save :: fact_LambdaBar=1.
  ! PURPOSE
  ! Enhancement factor of pbar p -> Lambda LambdaBar cross section
  ! (for larger statistics)
  !****************************************************************************

  !****************************************************************************
  !****g* barAntiBar/fact_JPsi
  ! SOURCE
  !
  real, save :: fact_JPsi=1.
  ! PURPOSE
  ! Enhancement factor of pbar p -> J/Psi cross section (for larger statistics)
  !****************************************************************************

  !****************************************************************************
  !****g* barAntiBar/fact_JPsi_width
  ! SOURCE
  !
  real, save :: fact_JPsi_width=1.
  ! PURPOSE
  ! Enhancement factor of the J/Psi total width (for larger statistics)
  !****************************************************************************

  !****************************************************************************
  !****g* barAntiBar/useAnni
  ! SOURCE
  !
  logical,save :: useAnni = .true.
  ! PURPOSE
  ! Flag whether to perform Baryon-Antibarion annihilation or not at all
  !****************************************************************************


  logical, save :: initFlag=.true.

  public :: sigmaBarAntiBar

contains

  !****************************************************************************
  !****s* barAntiBar/init
  ! NAME
  ! subroutine init
  ! PURPOSE
  ! Reads in namelist "barAntiBar_input"
  ! INPUTS
  ! * (none)
  ! OUTPUT
  ! * Initializes global module variables
  !****************************************************************************
  subroutine init
    use output

    integer :: ios

    !**************************************************************************
    !****n* barAntiBar/barAntiBar_input
    ! NAME
    ! NAMELIST barAntiBar_input
    ! PURPOSE
    ! Namelist which includes the input variables:
    ! * fact_LambdaBar
    ! * fact_JPsi
    ! * fact_JPsi_width
    ! * useAnni
    !**************************************************************************
    NAMELIST /barAntiBar_input/ fact_LambdaBar,fact_JPsi,fact_JPsi_width, &
         useAnni

    call Write_ReadingInput('barAntiBar_input',0)
    rewind(5)
    read(5,nml=barAntiBar_input,iostat=ios)
    call Write_ReadingInput('barAntiBar_input',0,ios)

    if (fact_LambdaBar.ne.1.) then
       write(*,*) ' ATTENTION: pbar p -> LambdaBar Lambda cross section'
       write(*,*) ' is rescaled by a factor of ', fact_LambdaBar
    end if

    if (fact_JPsi.ne.1.) then
       write(*,*) ' ATTENTION: pbar p -> J/Psi cross section'
       write(*,*) ' is rescaled by a factor of ', fact_JPsi
    end if

    if (fact_JPsi_width.ne.1.) then
       write(*,*) ' ATTENTION: J/Psi total width'
       write(*,*) ' is rescaled by a factor of ', fact_JPsi_width
    end if

    write(*,*) 'use Anni: ',useAnni

    call Write_ReadingInput('barAntiBar_input',1)

    initFlag = .false.
  end subroutine init


  !****************************************************************************
  !****s* barAntiBar/sigmaBarAntiBar
  ! NAME
  ! subroutine sigmaBarAntiBar(srts,teilchenIN,mediumATcollision,
  ! sigTotal,sigElast,sigCEX,sigAnni,
  ! sigLambdaBar,sigSigmaBar,sigXiBar,sigOmegaBar,
  ! sigJPsi,sigProduction,sigHyperon )
  !
  ! PURPOSE
  ! Computes total, elastic and other cross sections for baryon+antibaryon
  ! collisions.
  ! INPUTS
  ! * real :: srts --- sqrt(s) of collision (GeV)
  ! * type(particle),dimension(1:2) :: teilchenIn --- colliding particles
  ! * type(medium) :: mediumATcollision  --- medium infos at collision point
  ! OUTPUT
  ! * real :: sigTotal     --- total cross section (mb)
  ! * real :: sigElastic   --- elastic cross section (mb)
  ! * real :: sigCEX      ---
  !   charge exchange cross section (mb)
  ! * real :: sigAnni ---
  !   annihilation into mesons cross section (mb)
  ! * real :: sigLambdaBar ---
  !   exclusive (anti)hyperon production: Bbar+B -> Lambda+Lambdabar
  ! * real :: sigSigmaBar ---
  !   exclusive (anti)hyperon production :
  !   Bbar+B -> Lambda+Sigma0Bar, LambdaBar+Sigma0 (mb)
  ! * real :: sigXiBar ---
  !   exclusive (anti)cascade production: Bbar+B -> Xi+XiBar
  ! * real :: sigOmegaBar ---
  !   exclusive (anti)Omega (S=-3) production: Bbar+B -> Omega+OmegaBar
  ! * real :: sigJPsi ---
  !   J/Psi production cross section (mb): Bbar+B -> J/Psi
  ! * real, optional :: sigProduction   ---
  !   production Bbar+B -> B+Bbar+mesons cross section (mb)
  ! * real, optional :: sigHyperon      ---
  !   (anti)hyperon production Bbar+B -> Y+Ybar+mesons,
  !   B+Ybar+Kbar, Bbar+Y+K cross section (mb)
  !****************************************************************************
  subroutine sigmaBarAntiBar(srts,partIn,mediumAtColl, &
                           & sigTotal,sigElast,sigCEX,sigAnni, &
                           & sigLambdaBar,sigSigmaBar,sigXiBar,sigOmegaBar, &
                           & sigJPsi,sigProduction,sigHyperon )

    ! called with:
    ! * Tot,Elast,CEX,Anni,Lbar,Sbar,Xibar,Omegabar,Jpsi ! master_2body.f90
    ! * Tot,Elast,CEX,Anni,Lbar,Sbar,Xibar,Omegabar,Jpsi,Prod,Hyperon ! antiBarBar_main.f90


    use mediumDefinition
    use particleDefinition
    use constants, only: rhoNull, mN
    use twoBodyTools, only: pcm

    real, intent(in)  ::           srts
    type(particle),dimension(1:2), intent(in)    :: partIn
    type(medium),                  intent(in)    :: mediumAtColl
    real, intent(out) :: sigTotal
    real, intent(out) :: sigElast
    real, intent(out) :: sigCEX
    real, intent(out) :: sigAnni
    real, intent(out) :: sigLambdaBar
    real, intent(out) :: sigSigmaBar
    real, intent(out) :: sigXiBar
    real, intent(out) :: sigOmegaBar
    real, intent(out) :: sigJPsi
    real, optional, intent(out)  :: sigProduction
    real, optional, intent(out)  :: sigHyperon


    real :: m1,m2,s,p,u,p12
    real :: xsPROD, xsY
    integer :: totCharge, totId, nDelta
    integer :: I3_Delta

    ! momentum above which the PDG parameterization is working (GeV/c):
    real, parameter :: pCut=1.85
    real :: E,v_rel

    ! If .true. then annihilation cross section is medium-modified according to
    ! E. Hernandez and E. Oset, Z. Phys. A 341, 201 (1992)
    logical, parameter :: densityDependence_Oset=.false.


    if (initFlag) call init

    xsProd = 0.
    xsY = 0.

    s=srts**2
    m1=partIn(1)%mass
    m2=partIn(2)%mass
    E=(s-m1**2-m2**2)/(2.*m2)       ! Energy of 1-st particle in lab frame
    p=sqrt(E**2-m1**2)              ! Momentum of 1-st particle in lab frame
    v_rel=p/E
    p=v_rel*mN/sqrt(1-v_rel**2)  ! Momentum of antiproton for the same relative velocity

    p12=pcm(srts,m1,m2)   ! c.m. momentum of colliding particles

    !s0=4.*0.938**2
    !p=sqrt(max(s**2/s0-s,1e-06))           ! lab. momentum

    totCharge= sum(partIn(1:2)%charge)
    totId= sum(partIn(1:2)%ID)

    if (totId.eq.3) then
       if (partIn(1)%ID.eq.2) then
          nDelta=1
       else
          nDelta=2
       end if
       if (partIn(nDelta)%antiparticle) then
          I3_Delta=2*partIn(nDelta)%charge+1  ! twice I3 projection by Gell-Mann - Nishidjima formula
       else
          I3_Delta=2*partIn(nDelta)%charge-1
       end if
    end if

    sigElast = paramElast(p)
    sigCEX = paramCEX(totID,totCharge, p)
    sigAnni = paramAnni(totID,totCharge, p)

    if (densityDependence_Oset.and.mediumAtColl%useMedium) then
       u = mediumAtColl%density/rhoNull
       sigAnni = sigAnni * ( 1. + 4.59*u + 10.6*u**2 + 12.8*u**3 )
    end if

    if (.not.useAnni) sigAnni = 0.

    if (present(sigProduction) .or. p <= pCut) then
       xsProd = paramProd(totID, p)
    end if

    if (present(sigHyperon) .or. p <= pCut) then
       xsY = paramHyperon(totID, p)
    end if

    sigLambdaBar = paramLambdaBar(totID,totCharge, s)
    sigSigmaBar = paramSigmaBar(totID,totCharge, s, I3_Delta)
    sigXiBar = paramXiBar(totID,totCharge, s, I3_Delta)
    sigOmegaBar = paramOmegaBar(totID,totCharge, s)
    sigJPsi = paramJPsi(totID,totCharge, s, p12)

    if ( p <= pCut ) then
       sigTotal= sigElast + sigCEX + sigAnni + xsPROD + xsY
    else
       sigTotal = paramTotal(p)
    end if

    if (present(sigProduction)) sigProduction = xsPROD
    if (present(sigHyperon)) sigHyperon = xsY

    ! Enhancement factors (to get more statistics):
    sigLambdaBar=sigLambdaBar*fact_LambdaBar   ! ??? after sigTotal ???
    sigJPsi=sigJPsi*fact_JPsi                  ! ??????????????????????

  end subroutine sigmaBarAntiBar

  ! for the doku:
  !
  !    integer, intent(in) :: totID     ! sum(ID)
  !    integer, intent(in) :: totCharge ! sum(Charge)
  !    real, intent(in)    :: p         ! momentum
  !    real, intent(in)    :: s         ! sqrt^2


  real function paramTotal(p)
    real, intent(in)    :: p         ! momentum

    ! Parameterization of total cross sections from PDG,
    ! L. Montanet et al., PRD 50, 1173 (1994) (see p. 1335);
    ! see also T. Falter et al., PRC 70, 054609 (2004):

    paramTotal = 38.4+77.6*p**(-0.64)+0.26*log(p)**2-1.2*log(p)

  end function paramTotal

  real function paramElast(p)
    real, intent(in)    :: p         ! momentum

    paramElast = 0.

    if ( p < 2.03 ) then
       ! Parameterization from J. Cugnon and J. Vandermeulen,
       ! Annales de Physique (France) 14, 49 (1989);
       ! see also C.B. Dover et al., Prog. Part. Nucl. Phys. 29, 87 (1992):
       ! paramElast = 42.3/p**0.54 + 4.3*exp(-(p-1.5)**2)

       paramElast = 40.0/p**0.56 + 5.8*exp(-(p-1.85)**2)   ! A.L.

    else
       ! Use PDG parameterization:
       paramElast = 10.2+52.7*p**(-1.16)+0.125*log(p)**2-1.28*log(p)
    end if

  end function paramElast

  real function paramCEX(totID,totCharge, p)
    integer, intent(in) :: totID     ! sum(ID)
    integer, intent(in) :: totCharge ! sum(Charge)
    real, intent(in)    :: p         ! momentum

    paramCEX = 0.

    if ( (totId.eq.2 .and. totCharge.eq.0) &
         .or. &
         (totId.eq.3 .and. abs(totCharge).ne.2) ) then

       ! Parameterizations for pbar+p -> nbar+n from J. Cugnon and J. Vandermeulen,
       ! Annales de Physique (France) 14, 49 (1989);
       ! see also C.B. Dover et al., Prog. Part. Nucl. Phys. 29, 87 (1992):

       if ( p < 0.5 ) then
          if ( p <= 0.1 ) then
             !               paramCEX= 0.
          else
             paramCEX= 10.9*(p-0.1)/p**1.6
          end if
       else
          paramCEX= 7.1/p**0.9       ! This is actually not good at p > 3 GeV/c
          ! (above experiment)
       end if
    end if

  end function paramCEX

  real function paramAnni(totID,totCharge, p)
    integer, intent(in) :: totID     ! sum(ID)
    integer, intent(in) :: totCharge ! sum(Charge)
    real, intent(in)    :: p         ! momentum

    paramAnni = 0.

    if (p.lt.0.51) then
       if ( totId.eq.2 .and. totCharge.ne.0 &
            &.and. p.lt.0.382 ) then
          ! Parameterization for low energy nbar+p annihilation cross section
          ! from T. Armstrong et al., PRD 36, 659 (1987):
          paramAnni = 41.4 + 29./p
       else
          ! pbar+p data parameterization by A.L.:
          paramAnni = 51.52/p**0.85 + 0.034/p**2.94
       end if
    else if (p.lt.6.34) then
       paramAnni = 88.8/p**0.4 - 24.2             ! A.L.
    else
       paramAnni = 38./p**0.5 + 24./p**1.1        ! Cugnon
    end if

  end function paramAnni

  real function paramProd(totID, p)
    integer, intent(in) :: totID     ! sum(ID)
    real, intent(in)    :: p         ! momentum

    paramProd = 0.

    if (totId <= 3) then
       if ( p <= 0.793 ) then
          !           paramProd = 0.
       else
          paramProd = 30.*(p-0.793)**1.5/(2.+(p-0.793)**1.5)
       end if
    end if

  end function paramProd

  real function paramHyperon(totID, p)
    integer, intent(in) :: totID     ! sum(ID)
    real, intent(in)    :: p         ! momentum

    paramHyperon = 0.

    if (totId <= 3) then
       if ( p <= 1.435 ) then
          !              paramHyperon = 0.
       else
          paramHyperon = 3.*(p-1.435)/(10.+(p-1.435))
       end if
    end if

  end function paramHyperon

  real function paramLambdaBar(totID,totCharge, s)
    integer, intent(in) :: totID     ! sum(ID)
    integer, intent(in) :: totCharge ! sum(Charge)
    real, intent(in)    :: s         ! sqrt^2

    real :: exc ! excess energy

    paramLambdaBar = 0.

    if (totId.eq.2 .and. totCharge.eq.0) then

       if (s <= 4.982) then ! Lambda+LambdaBar Threshold [GeV**2]
          !            paramLambdaBar = 0.0
       else if (s < 5.071) then
          exc = sqrt(s) - 2.232
          paramLambdaBar = 0.0357345*sqrt(exc)+9.00539*exc**1.5
       else
          paramLambdaBar = max(0.010,0.725872*(s/4.982-1.)**0.774485*(4.982/s)**3.35044)
       end if

    end if
  end function paramLambdaBar

  real function paramSigmaBar(totID,totCharge, s, I3_Delta)
    integer, intent(in) :: totID     ! sum(ID)
    integer, intent(in) :: totCharge ! sum(Charge)
    real, intent(in)    :: s         ! sqrt^2
    integer, intent(in) :: I3_Delta  ! ...

    paramSigmaBar = 0.

    if ( totId <= 3 .and. abs(totCharge) <= 1 ) then

       if (s <= 5.327) then ! Lambda+Sigma0 Threshold [GeV**2]
          !             paramSigmaBar = 0.0
       else
          ! XS is fitted to Lambda+Sigma0Bar+c.c.
          !paramSigmaBar = 2.43292*(s/5.327-1.)**1.2983*(5.327/s)**11.1481
          paramSigmaBar = 0.183665*(s/5.327-1.)**0.436713*(5.327/s)**1.85006
       end if

       ! Multiplicative factors from isospin relations:
       if (totCharge /= 0) then
          if (totId .eq. 2) then
             paramSigmaBar=2.*paramSigmaBar
          else if (totId .eq. 3) then
             if (abs(I3_Delta).eq.3) then
                paramSigmaBar=1.5*paramSigmaBar
             else
                paramSigmaBar=0.5*paramSigmaBar
             end if
          end if
       end if
    end if
  end function paramSigmaBar

  real function paramXiBar(totID,totCharge, s, I3_Delta)
    integer, intent(in) :: totID     ! sum(ID)
    integer, intent(in) :: totCharge ! sum(Charge)
    real, intent(in)    :: s         ! sqrt^2
    integer, intent(in) :: I3_Delta  ! ...

    paramXiBar = 0.

    if ( totId <= 3 .and. abs(totCharge) <= 1 ) then
       if (s <= 6.9169) then        ! Xi+XiBar threshold [GeV^2]
          !             paramXiBar=0.
       else
          paramXiBar=0.004  !mb
       end if


       ! DeltaBar+N or Delta+Nbar collisions -->
       ! multiplicative factors from isospin relations:
       if (totId .eq. 3) then
          if (totCharge == 0) then
             paramXiBar=0.5*paramXiBar
          else if (abs(I3_Delta).eq.1) then
             paramXiBar=0.25*paramXiBar
          else
             paramXiBar=0.75*paramXiBar
          end if
       end if
    end if

  end function paramXiBar

  real function paramOmegaBar(totID,totCharge, s)
    integer, intent(in) :: totID     ! sum(ID)
    integer, intent(in) :: totCharge ! sum(Charge)
    real, intent(in)    :: s         ! sqrt^2

    paramOmegaBar = 0.

    if (totId==2 .and. totCharge == 0) then
       if (s <= 11.823) then        ! Omega+OmegaBar threshold [GeV^2]
!             paramOmegaBar=^0.0
       else
          paramOmegaBar=7.6e-05*(s/11.1823-1.)**1.7212*(11.1823/s)**6.5033
       end if
    end if

  end function paramOmegaBar

  real function paramJPsi(totID,totCharge, s, p12)

    use ParticleProperties, only: hadron
    use IdTable, only: JPsi
    use constants, only: pi

    integer, intent(in) :: totID     ! sum(ID)
    integer, intent(in) :: totCharge ! sum(Charge)
    real, intent(in)    :: s         ! sqrt^2
    real, intent(in)    :: p12       ! = pcm

    ! J/Psi total width (GeV) and branching ratio to pbar p (nbar n):
    real, parameter :: Gam0 = 92.9e-06
    real, parameter :: Br = 2.2e-03

    ! J/Psi partial width to pbar p (nbar n):
    real, parameter :: Gam_pbarp = Gam0 * Br

    real :: Gam

    paramJPsi = 0.

    if ( totId.eq.2 .and. totCharge.eq.0 ) then
       Gam = Gam0 * fact_JPsi_width
       paramJPsi = 3.*pi*0.389/p12**2*s*Gam_pbarp*Gam &
            &/((s-hadron(JPsi)%mass**2)**2+s*Gam**2)
    end if

  end function paramJPsi



end module barAntiBar
