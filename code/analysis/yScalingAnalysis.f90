!******************************************************************************
!****m* /yScalingAnalysis
! NAME
! module yScalingAnalysis
!
! PURPOSE
! * This module concludes the superscaling analysis of lepton induced processes
! * For details on superscaling see Ivan Lappo-Danilevski diploma thesis at
! * theorie.uni-giessen.de
!******************************************************************************

module yScalingAnalysis
  use initNeutrino, only: max_finalstate_ID
  implicit none
  private

  Public:: yScaling_Analyze

  !****************************************************************************
  !****g* yScalingAnalysis/debug
  ! SOURCE
  !
  logical, parameter :: debug=.false.
  ! PURPOSE
  ! switches on and off debug infos: If .true. debug infos are produced,
  ! if .false. not
  !****************************************************************************

  !****************************************************************************
  !****g* yScalingAnalysis/analyze
  ! SOURCE
  !
  logical, save :: analyze=.false.
  ! PURPOSE
  ! Determines wether the y-scaling analysis is performed
  !****************************************************************************

  !****************************************************************************
  !****g* yScalingAnalysis/optionalOutput
  ! SOURCE
  !
  logical, save :: optionalOutput=.false.
  ! PURPOSE
  ! Determines wether in addition to the standard 'scaling_analysis.dat' other histograms
  ! will be generated. E.g.
  ! * 'single_nucleon.dat' - a table for comparing nucleon-knockout with fully inclusive
  !   x sections
  ! * 'scaling_info.dat' - general parameters of the analysis, to be used for quick
  !   analyis
  ! * 'scaling_delta.dat' - output to be used for analyis of scaling function in
  !   resonance excitation region
  !****************************************************************************

  !****************************************************************************
  !****g* yScalingAnalysis/variable
  ! SOURCE
  !
  integer, save :: variable=1
  ! PURPOSE
  ! determines which kind of scaling variable will be used
  ! (cf. Donnelly, Sick 1999):
  ! * 1) RFG full variable Psi
  ! * 2) RFG approximation Psi
  ! * 3) PWIA full Upsilon (y/kf)
  ! * 9) evaluation will be done for all variables, output written to
  !   seperate files
  !****************************************************************************

  !****************************************************************************
  !****g* yScalingAnalysis/kFermi
  ! SOURCE
  !
  real, save :: kFermi=0.2251
  ! PURPOSE
  ! Nucleon Fermi momentum in nucleus. If none specified 0.2251 will be used, except
  ! if densitySwitch_static is set to 5, then fermiMomentum_input is used.
  ! The 0.225_1_ aims at preventing confusion whith delibaretely set differences between
  ! kFermi and fermiMomentum_input
  !****************************************************************************

  !****************************************************************************
  !****g* yScalingAnalysis/E_shift
  ! SOURCE
  !
  real, save :: E_shift=0.020
  ! PURPOSE
  ! Energy correction to account for binding effects, otherwise neglected in RFG model
  !****************************************************************************

  !****************************************************************************
  !****g* yScalingAnalysis/mA
  ! SOURCE
  !
  real, save :: mA=-1
  ! PURPOSE
  ! Tabulated mass of target nucleus
  !****************************************************************************

  !****************************************************************************
  !****g* yScalingAnalysis/mA_1
  ! SOURCE
  !
  real, save :: mA_1=-1
  ! PURPOSE
  ! Tabulated mass of rest nucleus (target nucleus with knocked out nucleon)
  !****************************************************************************

  logical, save :: initflag=.true.


contains

  !****************************************************************************
  !****s* yScalingAnalysis/readinput
  ! NAME
  ! subroutine readinput
  ! INPUTS
  ! NONE
  ! OUTPUT
  ! NONE
  ! PURPOSE
  ! This subroutine reads out the namelist "YScalingAnalysis".
  ! Only called once to initialize the module
  !****************************************************************************
  subroutine readinput
    use output
    use nucleusDefinition
    use nucleus, only: GetTarget

    integer :: IOS
    type(tnucleus),pointer :: NucInfo

    !**************************************************************************
    !****n* yScalingAnalysis/YScalingAnalysis
    ! NAME
    ! NAMELIST /YScalingAnalysis/
    ! PURPOSE
    ! This namelist includes:
    ! * analyze
    ! * optionalOutput
    ! * variable
    ! * kFermi
    ! * E_shift
    !**************************************************************************
    NAMELIST /yScalingAnalysis/ analyze, optionalOutput, variable, kFermi, E_shift

    call Write_ReadingInput('yScalingAnalysis',0)
    rewind(5) !5=jobcard
    read(5,nml=yScalingAnalysis,IOSTAT=IOS)
    call Write_ReadingInput('yScalingAnalysis',0,IOS) !message about success of reading

    NucInfo => GetTarget()
    if (NucInfo%densitySwitch_static == 5) then
      if (kFermi == 0.2251) then
        kFermi = NucInfo%fermiMomentum_input
      end if
    end if

    write(*,*) 'Performing yScaling Analysis: ', analyze
    if (analyze) then
       write(*,*) '   optionalOutput: ',optionalOutput
       write(*,*) '   variable:       ',variable
       write(*,*) '   kFermi:         ',kFermi
       write(*,*) '   E_shift:        ',E_shift
    end if

    call Write_ReadingInput('yScalingAnalysis',1)
  end subroutine readinput


  !****************************************************************************
  !****s* yScalingAnalysis/debug_prompt
  ! NAME
  ! subroutine debug_prompt
  ! INPUTS
  ! NONE
  ! OUTPUT
  ! NONE
  ! PURPOSE
  ! This subroutine prints out basic debug informations for the y-scaling analysis
  ! (Not in use at the moment)
  !****************************************************************************
!   subroutine debug_prompt()
!     write(*,*) '################### SCALING ANALYSIS DEBUG #######################'
!   end subroutine debug_prompt


  !****************************************************************************
  !****f* yScalingAnalysis/calculateY
  ! NAME
  ! function calculateY
  ! INPUTS
  ! * integer, intent(in) :: var  -- type of scaling variable (cf. variable)
  ! * real, intent(in) :: Li -- initial lepton energy
  ! * real, intent(in) :: Lf -- final lepton energy
  ! * real, intent(in) :: Qsquared -- four-momentum transfer
  ! OUTPUT
  ! * real :: calculateY -- resulting scaling variable
  ! PURPOSE
  ! This function calculates the scaling variable
  !****************************************************************************

  function calculateY(var,Li,Lf,Qsquared)
    use constants, only: mN


    integer, intent(in) :: var
    real, intent(in) :: Li !initial-state lepton energy
    real, intent(in) :: Lf !final-state lepton energy
    real, intent(in) :: Qsquared
    real :: calculateY
    real :: eta, lam, tau, kap, ws,q, psi0, omega, W

    q = sqrt((Li-Lf)**2 + Qsquared)
    ! meaning of sign: q = sqrt( energy transfer - four-momentum transfer^2 )

    kap = q/(2.*mN)                ! dimensionless momentum transfer
    eta = kFermi/mN                ! dimensionless Fermi momentum

    ws = Li - Lf - E_shift                           ! shifted energy transfer
    lam = ws/(2.*mN)               ! dimensionless energy transfer
    tau = kap**2 - lam**2                            ! dimensionless 4-momentum transfer

    if (debug) then
        write(*,*) 'yScaling : Kinematical variables q, omega_shift, lambda, tau, &
        & eta, kappa', q, ws, lam, tau, eta, kap
    end if

    select case (var)
    case (1) ! full RFG Psi
      calculateY = 1./sqrt(sqrt(1.+eta**2)-1.)*(lam-tau)/sqrt((1.+lam)*tau+kap*sqrt(tau*(1.+tau)))
    case (2) ! approximated RFG Psi
      psi0 = 2./eta*(sqrt(lam*(lam+1.))-kap)
      calculateY = psi0 *(1+sqrt( 1.+1/(4.*kap**2) )*1./2.*eta*psi0)
    case (3) ! full PWIA Upsilon
      omega = Li - Lf
      if ((mA < 0) .or. (mA_1 < 0)) then
          write(*,*) 'Wrong masses in yScalingAnalysis.f90. Target mass, rest mass: ', mA,mA_1
          write(*,*) 'STOP!'
          stop
      end if
      W = sqrt((mA+omega)**2-q**2)
      if (debug) write(*,*) 'Calculating PWIA scaling variable for mA, mA_1,W',mA,mA_1,W
      calculateY = 1./(2.*W**2)* ( (mA+omega) * sqrt(W**2 - (mA_1+mN)**2 )&
      & *sqrt(W**2 - (mA_1 -mN)**2)-q*(W**2+mA_1**2 - mN**2) ) / kFermi
    case default
      write(*,*) 'Wrong variable type in yScalingAnalysis.f90', variable
      write(*,*) 'STOP!'
      stop
    end select

  end function calculateY

  !****************************************************************************
  !****f* yScalingAnalysis/calculateYRes
  ! NAME
  ! function calculateYRes
  ! INPUTS
  ! * real, intent(in)              :: Li -- initial lepton energy
  ! * real, intent(in)              :: Lf -- final lepton energy
  ! * real, intent(in)              :: Qsquared -- four-momentum transfer
  ! * real, intent(in)              :: mass -- resonance mass
  ! OUTPUT
  ! * real                          :: calculateYRes -- resulting scaling variable
  ! PURPOSE
  ! This function calculates the scaling variable for the resonance region
  !****************************************************************************

  function calculateYRes(Li,Lf,Qsquared,mass)
    use constants, only: mN

    real, intent(in) :: Li !initial lepton Energy
    real, intent(in) :: Lf !final ----
    real, intent(in) :: Qsquared
    real, intent(in) :: mass
    real :: calculateYRes
    real :: eta, lam, tau, kap, ws,q, betha, rho !, psi0, omega, W

    q = sqrt((Li-Lf)**2 + Qsquared)
    ! meaning of sign: q = sqrt( energy transfer - four-momentum transfer^2 )

    kap = q/(2.*mN)                ! dimensionless momentum transfer
    eta = kFermi/mN                ! dimensionless Fermi momentum

    ws = Li - Lf - E_shift                           ! shifted energy transfer
    lam = ws/(2.*mN)               ! dimensionless energy transfer
    tau = kap**2 - lam**2                            ! dimensionless 4-momentum transfer

    ! cf. C. Maieron Physical Review C, Vol. 65, 025502
    betha = (mass**2 - mN**2)/(4.*mN**2)
    rho = 1. + betha/tau

    if (debug) then
        write(*,*) 'yScalingRes : Kinematical variables q, omega_shift, lambda, tau, eta, kappa', q, ws, lam, tau, eta, kap
    end if


    calculateYRes = 1./sqrt(sqrt(1.+eta**2)-1.)*(lam-tau*rho)/sqrt((1.+lam*rho)*tau+kap*sqrt(tau*(tau*rho**2+1.)))

  end function calculateYRes


  !****************************************************************************
  !****f* yScalingAnalysis/mottXsection
  ! NAME
  ! function mottXsection
  ! INPUTS
  ! * real, intent(in)                 :: Ei -- incoming lepton energy
  ! * real, intent(in)                 :: theta -- scattering angle
  ! * integer, intent(in)              :: Z -- target charge
  ! OUTPUT
  ! * real                             :: mottXsection
  ! PURPOSE
  ! This function calculates the Mott cross section for an
  ! ultrarelativistic spin 1/2 lepton. Result dimension is 1/Gev^2
  !****************************************************************************

  function mottXsection(Ei,theta,Z)

    use degRad_conversion, only: radian
    use constants, only: alphaQED

    real, intent(in)         :: Ei
    real, intent(in)         :: theta
    integer, intent(in), optional    :: Z
    real :: mottXsection
    real :: m
    real :: Zdummy ! Fortran: mandatory dummy variables, "mottXsection(Ei,theta,Z=1)" not allowed

    if (present (Z)) then
      Zdummy=Z
    else
      Zdummy=1
    end if

    if (debug) write(*,*) 'Calculating mott Xsection for Ei,theta,Z',Ei,theta,Zdummy

    m = ( alphaQED*Zdummy*cos(radian(theta/2.))/(2.*Ei*sin(radian(theta/2.))**2) )**2
    if (debug) write(*,*) m

    mottXsection = m
  end function mottXsection


  !****************************************************************************
  !****s* yScalingAnalysis/cc1Xsection
  ! NAME
  ! function cc1Xsection
  ! INPUTS
  ! * real, intent(in)                :: q -- 3-momentum transfer [GeV]
  ! * real, intent(in)                :: p -- momentum of struck nucleon [GeV]
  ! * real, intent(in)                :: Eps -- excitation energy [GeV]
  ! * real, intent(in)                :: mAi -- rest mass of target nucleus [GeV]
  ! * real, intent(in)                :: mA_1i -- rest mass of rest nucleus [GeV]
  ! OUTPUT
  ! * real, intent(out)               :: wl -- longitudinal response
  ! * real, intent(out)               :: wt -- transverse response
  ! PURPOSE
  ! This subroutine calculates the off-shell single-nucleon responses according to de Forests cc1 perscription
  ! (1983, Nucl. Phys. A  392, 232-248). But following the conventions of Donnelly and Sick
  ! (in Phys. Rev. C 60, 065502).
  !****************************************************************************
  subroutine cc1Xsection(q,p,Eps,m0A,m0A_1,wl,wt)

    use nucleusDefinition
    use nucleus, only: GetTarget
    use FF_QE_nucleonScattering, only: formfactors_QE
    use constants, only: mN

    real, intent(in)  :: q, p, Eps, m0A, m0A_1
    real, intent(out) :: wl, wt

    type(tnucleus),pointer :: NucInfo
    integer :: Z,N
    real :: E,En,Eb,kap,tau,taub,eta,del2,GEp,GMp,GEn,GMn,Ge2,Gm2,delG,W2,delW1,delW2, dummy1, dummy2, Qsquared !,W1

    NucInfo => GetTarget()

    Z = NucInfo%charge
    N = NucInfo%mass-NucInfo%charge

    E = m0A - sqrt(m0A_1**2+p**2) - Eps
    En = sqrt( (q+p)**2 + mN**2 )
    Eb = sqrt( mN**2 + p**2 )

    kap = q/(2.*mN)
    tau = kap**2 - ( (En - E)/(2.*mN) )**2
    taub = kap**2 - ( (En - Eb)/(2.*mN) )**2
    eta = p/mN

    del2 = taub/kap**2 * ( (En+Eb)/(2.* mN) )**2 - (1.+taub)

    ! Get form factors from most recent parametrization
    QSquared = tau * (4.*mN**2)
    call formfactors_QE(QSquared,1,0, dummy1, dummy2, GE=GEn,GM= GMn)
    call formfactors_QE(QSquared,1,1, dummy1, dummy2, GE=GEp, GM=GMp)

    Ge2 = Z*GEp**2 + N*GEn**2
    Gm2 = Z*GMp**2 + N*GMn**2
    delG = Z*GEp*GMp+N*GEn*GMn

!     W1 = tau*Gm2
    W2 = 1./(1.+tau)*(Ge2+tau*Gm2)
    delW1 = (taub-tau)/(1.+tau)**2*(Ge2+Gm2-2.*delG)
    delW2 = (taub-tau)/(1.+tau)**2*(Ge2-Gm2)

    wl = 1./(2.*kap*sqrt(1.+eta**2))* (kap**2/taub)* ( Ge2+del2*(W2+delW1)+ (1.+taub)*delW1+(1.+tau)*delW2 )
    wt = 1./(2.*kap*sqrt(1.+eta**2))* ( 2.*taub*Gm2+del2*(W2+delW1) )

  end subroutine cc1Xsection


  !****************************************************************************
  !****f* yScalingAnalysis/scalingFunction
  ! NAME
  ! function scalingFunction
  ! INPUTS
  ! * integer                          :: var -- type of scaling variable (cf. variable)
  ! * real                             :: Li -- initial lepton energy
  ! * real                             :: Lf -- final lepton energy
  ! * real                             :: Qsquared -- four-momentum transfer
  ! * integer                          :: process_ID -- 1,2,3 for EM,CC,NC see initNeutrino.f90
  ! * integer                          :: flavor_ID -- 1,2,3 for e,mu,tau
  ! OUTPUT
  ! * real                             :: scalingFunction
  ! PURPOSE
  ! This function calculates the scaling function f for different models
  !****************************************************************************
  real function scalingFunction(var,Xsection,Li,Lf,Qsquared,theta,process_ID,flavor_ID, GL, GT)

    use nucleus, only: GetTarget
    use FF_QE_nucleonScattering, only: formfactors_QE
    use nucleusDefinition
    use degRad_conversion, only: radian
    use constants, only: GeVSquared_times_mb,melec, mmuon,mtau, GF,coscab, pi, mN

    integer, intent(in) :: var
    real, intent(in) :: xsection, Li, Lf, Qsquared, theta
    integer, intent(in) :: process_ID, flavor_ID
    real, optional , intent(out) :: GL, GT

    real :: GEp, GMp, GEn, GMn, mott, q, kap, vt, vl, tau, delta, psi, &
            ksi, W2, wl,wt, y, dummy1, dummy2, ml, sig0 , tanthetaTh2, kp , FA,GA,GAp,FP,GP,GE,GM, &
            VCC,VCL,VLL,VTp, nu, rho, rhoP, XVVL,XAACL,XT,XTp, F1,F2, &
            RVVL, RVVT, RAACC, RAACL, RAALL, RAAT, RVATp, chi
    integer :: Z,N
    type(tnucleus),pointer :: NucInfo

    NucInfo => GetTarget()

    Z = NucInfo%charge
    N = NucInfo%mass-NucInfo%charge

    select case (process_ID)
    case (1)
      select case (flavor_ID)
      case (1)
        select case (var)
        case (1) ! full RFG Psi
          q = sqrt((Li-Lf)**2 + Qsquared)
          mott = mottXsection(Li,theta)

          kap = q/(2.*mN)
          tau = Qsquared/((2.*mN)**2)

          vl = (Qsquared/(q**2))**2
          vt = Qsquared/(2.*q**2)+(tan(radian(theta/2.)))**2

          ksi = sqrt( 1. + (kFermi/mN)**2 )-1.
          psi = calculateY(var,Li,Lf,Qsquared)
          delta = ksi*(1.-psi**2)*( sqrt(tau*(1.+tau))/kap +1./3.*ksi*(1.-psi**2)*tau/(kap**2) )

        ! Get form factors from most recent parametrization
        call formfactors_QE(QSquared,1,0, dummy1, dummy2, GE=GEn, GM=GMn)
        call formfactors_QE(QSquared,1,1, dummy1, dummy2, GE=GEp, GM=GMp)

          W2 = 1./(1.+tau)*( Z*GEp**2+N*GEn**2+tau*(Z* GMp**2 + N*GMn**2) )
          GL = (kap**2/tau)*(Z*GEp**2+N*GEn**2 + W2*delta)/( 2.*kap*(1.+ksi*(1.+psi**2)/2.) )
          GT = (2.*tau*(Z* GMp**2 + N*GMn**2) + W2*delta)/( 2.*kap*(1.+ksi*(1.+psi**2)/2.) )

          if (debug) write(*,*) 'Calculating scaling Function for q',q,'kap',kap,'tau',tau,'vl',vl,'vt',vt, &
            & 'GEp',GEp,'GMp',GMp,'GEn',GEn,'Gmn',Gmn, 'GL',GL,'GT',GT, &
            'W2',W2,'ksi',ksi,'psi',psi,'delta',delta,'Z',Z, 'N',N

          scalingFunction = kFermi * Xsection / ( mott*(vl*GL+vt*GT) )* &
           & GeVSquared_times_mb ! additional conversion since  Xsection is given in mbarn = fm^2 * 0.1

        case (2) ! approximate RFG Psi
          q = sqrt((Li-Lf)**2 + Qsquared)
          mott = mottXsection(Li,theta)

          kap = q/(2.*mN)
          tau = Qsquared/((2.*mN)**2)

          vl = (Qsquared/(q**2))**2
          vt = Qsquared/(2.*q**2)+(tan(radian(theta/2.)))**2

        ! Get form factors from most recent parametrization
        call formfactors_QE(QSquared,1,0, dummy1, dummy2, GE=GEn, GM=GMn)
        call formfactors_QE(QSquared,1,1, dummy1, dummy2, GE=GEp, GM=GMp)

          if (debug) write(*,*) 'Calculating scaling Function for q',q,'kap',kap,'tau',tau,'vl',vl,'vt',vt, &
            & 'GEp',GEp,'GMp',GMp,'GEn',GEn,'Gmn',Gmn, 'Z',Z, 'N',N

          scalingFunction = kFermi * Xsection / &
             & ( mott*( kap/(2.*tau)*vl*(Z*GEp**2+N*GEn**2) + tau/kap*vt*(Z* GMp**2 + N*GMn**2) ) ) * &
             & GeVSquared_times_mb ! additional conversion since  Xsection is given in mbarn = fm^2 * 0.1

        case (3) ! full PWIA Upsilon
          q = sqrt((Li-Lf)**2 + Qsquared)
          mott = mottXsection(Li,theta)

          vl = (Qsquared/(q**2))**2
          vt = Qsquared/(2.*q**2)+(tan(radian(theta/2.)))**2

          y = calculateY(var,Li,Lf,Qsquared)*kFermi

          call cc1Xsection(q,-y,0.,mA,mA_1,wl,wt)
          !evaluate off-shell prescription at minimum nucleon momentum and minimum excitation energy

          if (debug) write(*,*) 'Calculating scaling Function for q,vl,vt,wl,wt:', &
           &  q,vl,vt,wl,wt

          scalingFunction = kFermi * Xsection / ( mott*(vl*wl+vt*wt) )* &
          & GeVSquared_times_mb ! additional conversion since  Xsection is given in mbarn = fm^2 * 0.1

        case default
          write(*,*) 'Wrong variable type in yScalingAnalysis.f90', variable
          write(*,*) 'STOP!'
          stop
        end select
      case default
        write(*,*) 'Unsupported particle flavor',flavor_ID,'for process ID',process_ID,'in yScalingAnalysis.f90'
        write(*,*) 'STOP!'
        stop
      end select
    case (2)
    select case (flavor_ID)
      case (1)
        ml = melec
      case (2)
        ml = mmuon
      case (3)
        ml = mtau
    end select
    select case (var)
      case (1)
        chi=1.

        q = sqrt((Li-Lf)**2 + Qsquared)
        tau = Qsquared/((2.*mN)**2)

        if (4.*Li*Lf == Qsquared) write(*,*) 'WARNING : backward scattering yields infinite scaling function'
        tanthetaTh2= Qsquared/(4*Li*Lf-Qsquared) ! tan(theta~/2)**2, while theta~ approx theta

        kp = sqrt(Lf**2-ml**2)
        nu=(Li-Lf)/q
        rho = 1. - nu**2
        rhoP=q/(Li+Lf)
        delta = ml/sqrt(Qsquared)

        VCC=1.-tanthetaTh2*delta**2
        VCL=nu + 1./rhoP*tanthetaTh2*delta**2
        VLL=nu**2+tanthetaTh2*(1.+2.*nu/rhoP + rho*delta**2)*delta**2
        VT=(rho/2.+tanthetaTh2)-1./rhoP*tanthetaTh2*(nu+rho*rhoP/2.*delta**2)*delta**2
        VTp=tanthetaTh2/rhoP*(1.-nu*rhoP*delta**2)
        VL=VCC-2.*nu*VCL+nu**2*VLL

        sig0 = (GF *coscab)**2 / (2.*pi**2)* ( kp * cos(atan2(sqrt(Qsquared),sqrt(4.*Li*Lf-Qsquared))) )**2
        call formfactors_QE(QSquared,2,0, F1, F2,FA=FA,FP=FP,GE=GE, GM=GM)
        GA=FA; GP= 2.*FP ! conversion due to difference in definition between Leitner diss. 4.1
                                   ! and Amaro C71, 015501 (2005)
        GAp=GA-tau*GP

        RVVL = 1./rho * GE**2
        RVVT = 2.*tau*GM**2
        RAACC = nu**2/rho * GAp**2
        RAACL = -nu/rho*GAp**2
        RAALL=1./rho*GAp**2
        RAAT=2.*(1.+tau)*GA**2
        RVATp=2.*sqrt(tau*(1.+tau))*GM*GA

        XVVL=VL*RVVL
        XAACL=VCC*RAACC+2.*VCL*RAACL+VLL*RAALL
!        XAACL = tanthetaTh2*GAp**2*(1.+delta**2)*delta**2
        XT=VT*(RVVT + RAAT)
        XTp=2.*VTp*RVATp

        if (debug) write(*,*) 'Calculating scaling function for neutrino case:', &
         & 'ml',ml,'q', q,'tau', tau,'tanthetaTh2',tanthetaTh2,'kp', kp,'nu', nu ,'rho',rho,'rhoP',rhoP,'delta',delta, &
         & 'VCC', VCC, 'VCL', VCL, 'VLL', VLL,'VT',VT,'VTp',Vtp,'VL', &
         & VL,'sig0',sig0,'GE',GE,'GM',GM,'GA', GA, 'GAp',GAp,'XVVL',XVVL, &
         & 'XAACL', XAACL, 'XT', XT, 'XTp',XTp,'F1',F1,'F2',F2

        scalingFunction = 1/2.* kFermi*q/(mN*N) * Xsection / (sig0 * (XVVL+XAACL+XT+chi*XTp) ) * &
        ! additional factor 1/2
         & GeVSquared_times_mb ! additional conversion since  Xsection is given in mbarn = fm^2 * 0.1
      case (2)
        chi=1.

        q = sqrt((Li-Lf)**2 + Qsquared)
        tau = Qsquared/((2.*mN)**2)

        if (4.*Li*Lf == Qsquared) write(*,*) 'WARNING : backward scattering yields infinite scaling function'
        tanthetaTh2= tan(radian(theta)/2)**2

        kp = sqrt(Lf**2-ml**2)
        nu=(Li-Lf)/q
        rho = 1. - nu**2
        rhoP=q/(Li+Lf)
        delta = ml/sqrt(Qsquared)

        VCC=1.
        VCL=nu
        VLL=nu**2
        VT=(rho/2.+tanthetaTh2)
        VTp=sqrt(tanthetaTh2*(rho+tanthetaTh2))
        VL=rho**2

        sig0 = (GF *coscab)**2 / (2.*pi**2)* ( kp * cos(atan2(sqrt(Qsquared),sqrt(4.*Li*Lf-Qsquared))) )**2
        call formfactors_QE(QSquared,2,0, F1, F2,FA=FA,FP=FP,GE=GE, GM=GM)
        GA=FA; GP= 2.*FP ! conversion due to difference in definition between Leitner diss. 4.1
                                   ! and Amaro C71, 015501 (2005)
        GAp=GA-tau*GP

        RVVL = 1./rho * GE**2
        RVVT = 2.*tau*GM**2
        RAACC = nu**2/rho * GAp**2
        RAACL = -nu/rho*GAp**2
        RAALL=1./rho*GAp**2
        RAAT=2.*(1.+tau)*GA**2
        RVATp=2.*sqrt(tau*(1.+tau))*GM*GA

        XVVL=VL*RVVL
        XAACL=0
!        XAACL=VCC*RAACC+2.*VCL*RAACL+VLL*RAAL
!        XAACL = tanthetaTh2*GAp**2*(1.+delta**2)*delta**2
        XT=VT*(RVVT + RAAT)
        XTp=2.*VTp*RVATp

        if (debug) write(*,*) 'Calculating scaling function for neutrino case:', &
         & 'ml',ml,'q', q,'tau', tau,'tanthetaTh2',tanthetaTh2,'kp', kp,'nu', nu ,'rho',rho,'rhoP',rhoP,'delta',delta, &
         & 'VCC', VCC, 'VCL', VCL, 'VLL', VLL,'VT',VT,'VTp',Vtp,  &
         & 'VL',VL,'sig0',sig0,'GE',GE,'GM',GM,'GA', GA, 'GAp',GAp,'XVVL',XVVL, &
         & 'XAACL', XAACL, 'XT', XT, 'XTp',XTp,'F1',F1,'F2',F2

        scalingFunction = 1/2.* kFermi*q/(mN*N) * Xsection / (sig0 * (XVVL+XAACL+XT+chi*XTp) ) * &
         ! additional factor 1/2
         & GeVSquared_times_mb ! additional conversion since  Xsection is given in mbarn = fm^2 * 0.1
      case (3)
        write(*,*) 'Ignoring request to perform scaling-function computation for PWIA scenario.'
        write(*,*) 'Setting scaling function to 0'
        scalingFunction = 0
      case default
        write(*,*) 'Wrong variable type in yScalingAnalysis.f90', variable
        write(*,*) 'STOP!'
        stop
      end select
    end select

  end function scalingFunction

  !****************************************************************************
  !****s* yScalingAnalysis/scaling_Analyze
  ! NAME
  ! subroutine scaling_Analyze(Particles,finalFlag,num_runs_sameEnergy,eventType)
  ! INPUTS
  ! * particles
  ! * finalFlag
  ! * num_runs_sameEnergy
  ! * eventType
  ! This subroutine produces output for lepton-nucleus y-scaling analysis
  !****************************************************************************
  subroutine yScaling_Analyze(Particles,finalFlag,eventType)

    use lowElectron, only: le_get_energy_li, &
                             & le_get_energy_lf, le_get_theta_lf, le_get_QSquared, low_Ele_incXsection
    use leptonKinematics, only: evaluate_Qsquared
    use initNeutrino, only: get_init_namelist, get_runtime_vars
    use neutrinoXsection, only: get_xsection_namelist
    use particleDefinition
    use particleProperties, only: hadron
    use IDTable, only: nucleon,delta
    use eventtypes, only: Neutrino,LoLepton
    use nucleus, only: GetTarget
    use nucleusDefinition
    use tabNuclearMass, only: read_mass
    use constants, only: GeVSquared_times_mb, pi
    use degRad_conversion, only: degrees

    type(particle), intent(in), dimension(:,:) :: Particles
    logical, intent(in) :: finalFlag
    integer, intent(in) :: eventType

    ! Local variables:
    integer :: i,j
    integer,save :: numberOfCalls=0
    integer,save :: numberOfFinals=0
    integer, save :: Z,N, process_ID, flavor_ID, nuXsectionMode !target and collision properties
    real, dimension(0:1) ::            tsigmanucleon=0.
    real,dimension(0:1),save ::        sum_tsigmanucleon=0.
    real,dimension(0:max_finalstate_ID) :: sigabsArr

    real :: totX = 0.

    real :: lep_en_f,lep_en_i,lep_QSquared,lep_theta_f,y,scafu, mott, add_xsection, GL, GT
    type(tnucleus),pointer :: NucInfo

    if (initflag) then

      call readinput
      initflag=.false.
      if (.not.analyze) return

      if (eventType == Neutrino) then
        call get_init_namelist(Gprocess_ID=process_ID,Gflavor_ID= flavor_ID,GnuXsectionMode=nuXsectionMode)
        if (nuXsectionMode /= 1) then
              write(*,*) 'Unsupported nuXsectionMode in yScalingAnalysis.f90: ', nuXsectionMode
              write(*,*) 'STOP!'
              stop
        end if
        if (process_ID == 1) then     !! this check is also performed in scalingFunction
            if (flavor_ID /= 1) then
              write(*,*) 'Unsupported flavor_ID in yScalingAnalysis.f90: ', flavor_ID
              write(*,*) 'STOP!'
              stop
            end if
         end if
      end if

      if (eventType == LoLepton) then
        process_ID = 1
        flavor_ID = 1


      end if

      NucInfo => GetTarget()
      N = NucInfo%mass-NucInfo%charge
      Z = NucInfo%charge

      if ((variable==3) .or. (variable==9)) then !get tabulated mass data in case of PWIA anylsis
        mA = read_mass(N,Z)
        mA_1 = read_mass(N,Z-1)
        write(*,*) 'Initial mass',mA,'filal mass',mA_1
      end if

      numberOfCalls=0
      numberOfFinals=0

      open(10,File='scaling_analysis.dat')
      write(10,*)'# Doing scaling analysis (cf. Donnelly, Sick 1999) for variable', variable
      write(10,*)'# 1) for RFG full variable Psi'
      write(10,*)'# 2) for approximated RFG variable Psi'
      write(10,*)'# 3) for full PWIA variable Upsilon'
      write(10,*)'# 9) for all analyses t0 be perormed simultaneously and the output to be stored'
      write(10,*)'#    in separate files scaling_analysis_<variable>.dat'
      write(10,*)'#'
      write(10,*)'# order of entries: 1) scaling variable; 2) scaling function;'
      write(10,*)'#                   3) li;  4) lf; 5) Q^2; 6) theta;'
      write(10,*)'#                   7) inclusive double diff cross section; 8) Mott x section;'
      write(10,*)'#                   9) lonitudinal formfactor GL; 10) transversal formfactor GT'
      write(10,*)'# Units are GeV(^2) or none, except for [dsigma/(dOmega*dE)]: mb/Gev'
      close(10)

      if (optionalOutput) then
        open(10,File='single_nucleon.dat')
        write(10,*)'# Comparing single-nucleon-knockout with fully inclusive cross section'
        write(10,*)'# NOTE: should only differ for timesteps > 0, due to formation time of res. and FSI'
        write(10,*)'# order of entries: incoming lepton energy, outgoing lepton energy,'
        write(10,*)'#       sigma single-neutron, sigma single-proton, sigma single-nucleon,'
        write(10,*)'#       sigma inclusive, ratio sig_sn/sig_inc'
        write(10,*)'# Units are [Gev] and [mb/Gev]'
        close(10)

        open(10,File='scaling_delta.dat')
        write(10,*)'# Doing full RFG scaling analysis in resonance region'
        write(10,*)'# (cf. Maieron, Donnelly, Sick 2002)'
        write(10,*)'# Scaling function is extracted using analyisis scripts'
        write(10,*)'# order of entries: li, lf, Q^2, theta, inclusive double diff cross section, '
        write(10,*)'#                   Mott x section, scaling variable'
        write(10,*)'# Units are GeV(^2) or none, except for [dsigma/(dOmega*dE*A)]: mb/Gev'
        close(10)

        !**********************************************************************
        !****o* yScalingAnalysis/scaling_info.dat
        ! NAME
        ! file scaling_info.dat
        ! PURPOSE
        ! * Collection of information concerning the scaling analysis
        ! * order of entries: kFermi, E_shift, N ,Z
        ! NOTE
        !**********************************************************************

        open(10,File='scaling_info.dat')
        write(10,*)'# Collection of information concerning the scaling analysis'
        write(10,*)'# order of entries: kFermi, N, Z'
        write(10,*)'# Units are GeV(^2) or none, except for [dsigma/(dOmega*dE*A)]: mb/Gev'
        write(10,'(10g13.5)') kFermi, E_shift, N, Z
        close(10)

      end if

    end if

    if (.not.analyze) return ! necessary since analysis called every run and not only on init

    numberOfCalls=numberOfCalls+1

    write(*,*) '################### SCALING ANALYSIS #######################'
    write(*,*) ' number of calls: ', numberofCalls,' number of finals: ', numberofFinals


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !nullify cross section
    tsigmanucleon=0.

    !calculate nucleon knockout cross-section
    do i=lbound(Particles,dim=1),ubound(Particles,dim=1)
      do j=lbound(Particles,dim=2),ubound(Particles,dim=2)
        if (Particles(i,j)%ID.le.0) cycle

        if (Particles(i,j)%ID.eq.nucleon) then

          !total cross section
          tsigmanucleon(Particles(i,j)%charge)=tsigmanucleon(Particles(i,j)%charge)+Particles(i,j)%perweight

        end if

      end do
    end do

    sum_tsigmanucleon=sum_tsigmanucleon+tsigmanucleon

    select case (eventType)
    case (LoLepton)
      add_xsection = low_Ele_incXsection
    case (Neutrino)
      call get_runtime_vars(GsigabsArr=sigabsArr)
      select case (process_ID)
      case (1)
        ! conversion since neutrinoXsection counts dsigma/dOmega in nanobarn, not dsigma/dOmega in milibarn
        ! for EM processes; see "end subroutine XsecdCosthetadElepton" in neutrinoXsection
        add_xsection = sigabsArr(0)/(2*pi*1E6)
      case (2)
        ! conversion since neutrinoXsection counts dsigma/dOmega in 10^-38 cm^2,
        ! not dsigma/dOmega in milibarn (10^-27 cm^2)
        ! for CC/NC processes; see "end subroutine XsecdCosthetadElepton" in neutrinoXsection
        add_xsection = sigabsArr(0)/(2*pi*1E11)
      end select
    end select

    totX = totX + add_xsection

    ! At the end of each set of same energy runs the analyis is performed and written
    if (finalFlag) then

      ! average the summed up cross section, multiply by A, since cross sections
      ! are given as mb/Gev per nucleon in LowElectron
      totX=totX/real(numberofcalls)*(N+Z)

    select case (eventType)
    case (LoLepton)
      lep_en_i = le_get_energy_li()
      lep_en_f = le_get_energy_lf()
      lep_QSquared = le_get_QSquared()
      lep_theta_f = le_get_theta_lf()
    case (Neutrino)
      call get_xsection_namelist(Genu=lep_en_i,Gelepton=lep_en_f,Gcostheta=lep_theta_f)
      lep_theta_f=degrees(acos(lep_theta_f))
      lep_Qsquared=evaluate_Qsquared(lep_theta_f,lep_en_i,lep_en_f)
    end select

      mott=mottXsection(lep_en_i,lep_theta_f)/GeVSquared_times_mb

      select case (variable)
        case (1:3)
          !********************************************************************
          !****o* yScalingAnalysis/scaling_analysis.dat
          ! NAME
          ! file scaling_analysis.dat
          ! PURPOSE
          !   Evaluate scaling variable and scaling function (as described by Donnelly
          !   and Sick in Phys. Rev. C 60, 065502).
          ! PURPOSE
          ! Columns:
          ! * 1 : scaling variable
          ! * 2 : scaling function
          ! * 3 : initial lepton energy [Gev]
          ! * 4 : final lepton energy [Gev]
          ! * 5 : momentum transfer Q^2 [Gev^2]
          ! * 6 : angle of outgoing lepton [°]
          ! * 7 : inclusive cross section [mb/GeV]
          ! * 8 : mott cross section [mb/GeV]
          ! * 9 : Longitudinal single-nucleon response GL
          ! * 10: Transversal single-nucleon response GT
          !********************************************************************
          open(11,File='scaling_analysis.dat',position='append')
          y = calculateY(variable,lep_en_i, lep_en_f, lep_QSquared)
          scafu = scalingFunction(variable, totX, lep_en_i, lep_en_f, &
          & lep_QSquared, lep_theta_f,process_ID,flavor_ID, GL=GL,GT=GT)
          write(11,'(10g13.5)') y, scafu,lep_en_i, lep_en_f, lep_QSquared, lep_theta_f, totX, mott,GL,GT
          close(11)
        case (9)
          !********************************************************************
          !****o* yScalingAnalysis/scaling_analysis_fullRFG.dat
          ! NAME
          ! file scaling_analysis_fullRFG.dat
          ! PURPOSE
          !   Specific analysis with full RFG scaling variable computation
          !   See : yScalingAnalysis/scaling_analysis.dat
          !********************************************************************
          open(11,File='scaling_analysis_fullRFG.dat',position='append')
          y = calculateY(1,lep_en_i, lep_en_f, lep_QSquared)
          scafu = scalingFunction(1, totX, lep_en_i, lep_en_f, lep_QSquared, lep_theta_f, &
          & process_ID,flavor_ID, GL=GL,GT=GT)
          write(11,'(10g13.5)') y, scafu,lep_en_i, lep_en_f, lep_QSquared, lep_theta_f, totX, mott,GL,GT
          close(11)
          !********************************************************************
          !****o* yScalingAnalysis/scaling_analysis_approxRFG.dat
          ! NAME
          ! file scaling_analysis_approxRFG.dat
          ! PURPOSE
          !   Specific analysis with approximate RFG scaling variable computation
          !   See : yScalingAnalysis/scaling_analysis.dat
          !********************************************************************
          open(11,File='scaling_analysis_approxRFG.dat',position='append')
          y = calculateY(2,lep_en_i, lep_en_f, lep_QSquared)
          scafu = scalingFunction(2, totX, lep_en_i, lep_en_f, lep_QSquared, lep_theta_f, &
          & process_ID,flavor_ID, GL=GL,GT=GT)
          write(11,'(10g13.5)') y, scafu,lep_en_i, lep_en_f, lep_QSquared, lep_theta_f, totX, mott,GL,GT
          close(11)
          !********************************************************************
          !****o* yScalingAnalysis/scaling_analysis_fullPWIA.dat
          ! NAME
          ! file scaling_analysis_fullPWIA.dat
          ! PURPOSE
          !   Specific analysis with full PWIA scaling variable computation
          !   See : yScalingAnalysis/scaling_analysis.dat
          !********************************************************************
          open(11,File='scaling_analysis_fullPWIA.dat',position='append')
          y = calculateY(3,lep_en_i, lep_en_f, lep_QSquared)
          scafu = scalingFunction(3, totX, lep_en_i, lep_en_f, lep_QSquared, lep_theta_f, &
          & process_ID,flavor_ID, GL=GL,GT=GT)
          write(11,'(10g13.5)') y, scafu,lep_en_i, lep_en_f, lep_QSquared, lep_theta_f, totX, mott,GL,GT
          close(11)
        case default
          write(*,*) 'Wrong variable type in yScalingAnalysis.f90', variable
          write(*,*) 'STOP!'
          stop
      end select

      if (optionalOutput) then
        !**********************************************************************
        !****o* yScalingAnalysis/single_nucleon.dat
        ! NAME
        ! file single_nucleon.dat
        ! PURPOSE
        ! * A table for comparing nucleon-knockout with fully inclusive x sections
        ! * order of entries: incoming lepton energy, outgoing lepton energy,
        !      sigma single-neutron, sigma single-proton, sigma single-nucleon
        !      sigma inclusive, ratio sig_sn/sig_inc
        ! * Units are [mb/Gev]
        ! NOTE
        ! should only differ significantly for timesteps > 0
        !**********************************************************************
        sum_tsigmanucleon=sum_tsigmanucleon*(N+Z)/real(numberofcalls)
        open(11,File='single_nucleon.dat',position='append')
        write(11,'(10g13.5)') lep_en_i,lep_en_f, sum_tsigmanucleon(0)/real(numberofcalls), &
        & sum_tsigmanucleon(1), &
        & (sum_tsigmanucleon(0)+ sum_tsigmanucleon(1)), totX, &
        & (sum_tsigmanucleon(0)+ sum_tsigmanucleon(1))/totX
        close(11)

          !********************************************************************
          !****o* yScalingAnalysis/scaling_delta.dat
          ! NAME
          ! file scaling_delta.dat
          ! PURPOSE
          !   Evaluate scaling variable for the resonance region and particulary
          !   the delta (as described by Maieron, Donnelly
          !   and Sick in Phys. Rev. C 65, 025502).
          ! PURPOSE
          ! Columns:
          ! * 1 : initial lepton energy [Gev]
          ! * 2 : final lepton energy [Gev]
          ! * 3 : momentum transfer Q^2 [Gev^2]
          ! * 4 : angle of outgoing lepton [°]
          ! * 5 : inclusive cross section [mb/GeV]
          ! * 6 : mott cross section [mb/GeV]
          ! * 7 : scaling variable
          ! * [8] : scaling function (to be filled in by 'delta_scaling.py'
          ! * [9] : cross section as estimated by simple RFG model (to be filled in by 'delta_scaling.py')
          !********************************************************************
        sum_tsigmanucleon=sum_tsigmanucleon*(N+Z)/real(numberofcalls)
        open(11,File='scaling_delta.dat',position='append')
        y = calculateYRes(lep_en_i, lep_en_f, lep_QSquared,hadron(delta)%mass)
        write(11,'(10g13.5)') lep_en_i, lep_en_f, lep_QSquared, lep_theta_f, totX, mott, y
        close(11)
      end if
      sum_tsigmanucleon=0.
      totx=0.
      numberOfCalls=0
      numberOfFinals=numberOfFinals+1
    end if

    write(*,*) '################### SCALING ANALYSIS FINISHED #######################'

  end subroutine yScaling_Analyze

end module yScalingAnalysis
