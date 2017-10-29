!******************************************************************************
!****m* /FF_QE_nucleonScattering
! NAME
! module FF_QE_nucleonScattering
! PURPOSE
! Provides the electro-weak form factors for the process
! lepton nucleon -> lepton' nucleon'
! and the Sachs form factors.
!******************************************************************************
module FF_QE_nucleonScattering

  implicit none
  private

  !****************************************************************************
  !****m* FF_QE_nucleonScattering/parametrization
  ! SOURCE
  !
  integer,save :: parametrization=3
  !
  ! PURPOSE
  ! * 0 = dipole approximation
  ! * 1 = BBA03 parametrization
  ! * 2 = BBBA05 parametrization
  ! * 3 = BBBA07 parametrization
  !****************************************************************************


  !****************************************************************************
  !****m* FF_QE_nucleonScattering/useNonStandardMA
  ! SOURCE
  !
  logical,save::useNonStandardMA=.false.
  !
  ! PURPOSE
  ! if one wants to use a specific axial mass, set this to true and choose value
  ! for MA_in
  !****************************************************************************


  !****************************************************************************
  !****m* FF_QE_nucleonScattering/MA_in
  ! SOURCE
  !
  real,save :: MA_in=1.0
  !
  ! PURPOSE
  ! axial mass (only if useNonStandardMA=.true.)
  !****************************************************************************

  !****************************************************************************
  !****m* FF_QE_nucleonScattering/MV2
  ! SOURCE
  !
  real,save :: MV2=0.71
  !
  ! PURPOSE
  ! vector mass squared in the dipole parametrization of the vector form
  ! factors
  !****************************************************************************



  !****************************************************************************
  !****m* FF_QE_nucleonScattering/deltas
  ! SOURCE
  !
  real,save :: deltas=-0.15
  !
  ! PURPOSE
  ! strange contribution to the axial ff.
  !****************************************************************************

  !****************************************************************************
  !****m* FF_QE_nucleonScattering/axialMonopole
  ! SOURCE
  !
  logical,save :: axialMonopole=.false.
  !
  ! PURPOSE
  ! use axial ff. of Gari, Kaulfuss PLB 138 (1984)
  !****************************************************************************

  ! constants: magnetic moments of proton and neutron
  real, parameter :: mup = 2.793
  real, parameter :: mun = -1.913


  logical,save :: initFlag=.true.

  public :: formfactors_QE

contains

  subroutine initInput
    use output

    integer :: ios

    !**************************************************************************
    !****n* FF_QE_nucleonScattering/ff_QE
    ! NAME
    ! NAMELIST ff_QE
    ! PURPOSE
    ! Includes parameters for neutrino matrix elements:
    ! * parametrization
    ! * MV2
    ! * MA_in
    ! * useNonStandardMA
    ! * deltas
    ! * axialMonopole
    !**************************************************************************
    NAMELIST /ff_QE/ parametrization,MV2,MA_in,&
         useNonStandardMA,deltas,axialMonopole

    call Write_ReadingInput('ff_QE',0)

    rewind(5)
    read(5,nml=ff_QE,IOSTAT=ios)
    call Write_ReadingInput("ff_QE",0,ios)

    write(*,*) 'Formfactor Parametrization for QE:'
    select case (parametrization)
    case (0)
       write(*,*) '     => dipole approximation'
    case (1)
       write(*,*) '     => BBA03 parametrization'
    case (2)
       write(*,*) '     => BBBA05 parametrization'
    case (3)
       write(*,*) '     => BBBA07 parametrization'
    case default
       write(*,*) ' strange value for QE FF parametrization -> STOP', parametrization
       stop
    end select

    if (useNonStandardMA) write(*,*) 'useNonStandardMA: nucleon axial mass: ', MA_in

    write(*,*) 'deltas= ', deltas

    if (axialMonopole) write(*,*) 'use axial ff. of Gari, Kaulfuss PLB 138 (1984)'

    call Write_ReadingInput('ff_QE',1)

  end subroutine initInput




  !****************************************************************************
  !****s* FF_QE_nucleonScattering/formfactors_QE
  ! NAME
  ! subroutine formfactors_QE(QSquared,processID,initialState_charge,F1,F2,FA,FP)
  !
  ! PURPOSE
  ! Provides the electro-weak form factors for the process
  ! lepton nucleon -> lepton' nucleon'.
  !
  ! INPUTS
  ! * real, intent(in) ::    QSQuared ! =Q^2 : virtuality of gauge boson in units of GeV**2
  !                                             =(-1)*Mandelstam-t
  ! * integer, intent(in) :: initialState_charge ! Specifies the charge of the incoming nucleon (either 1 or 0)
  ! * integer, intent(in) :: processID
  !
  ! "processID" specifies the reaction type: EM, CC, NC
  !
  ! OUTPUT
  ! * real, intent(out) :: F1,F2            ! Vector Form factors
  ! * real, intent(out),optional :: FA,FP   ! Axial Form factors (not necessary as input if processID=3)
  !****************************************************************************
  subroutine formfactors_QE(QSquared,processID,initialState_charge,F1,F2,FA,FP,GE,GM)
    use constants, only: sinsthweinbg, mN, MPi
    use leptonicID

    real, intent(in) :: QSQuared ! =Q^2 : virtuality of gauge boson
    integer, intent(in) :: processID, initialState_charge ! Specifies the reaction type

    real, intent(out) :: F1,F2   ! Vector Form factors
    real, intent(out),optional :: FA,FP   ! Axial Form factors


    ! Sachs form factors:
    real :: GEp,GMp,GEn,GMn
    real,optional ::GE,GM

    real :: tau
    real :: MAS, F10S, F20S,MVS
    real :: FQs, F1S, F2S, FAS
    real :: F1p, F1n, F2p, F2n
    real :: tau3
    real :: MA


    ! axial parameters:
    real, parameter :: gA=-1.2695 !PDG value
    real, parameter :: MA03=1.0   !Bodek fit BBA03 hep-ex/0308005
    real, parameter :: MA07=0.999 !Kuzmin fit based on BBBA07: arXiv:0712.4384 [hep-ph]

    real, parameter :: lambda1=0.85
    real, parameter :: lambda2=1.38


    if (initFlag) then
       call initInput
       initFlag=.false.
    end if

    !*** Check input:
    if (.not.((initialState_charge.eq.0).or.(initialState_charge.eq.1))) then
       write(*,*) 'Error in formfactors_QE! Strange initialState_charge',initialState_charge
       stop
    end if

    if (QSquared.lt.0) then
       write(*,*) 'QSquared less than zero in formfactors_QE. STOP!', Qsquared
       stop
    end if

    !*** Set Sachs form factors:

    tau = QSquared/(4.*mN**2)

    select case (parametrization)
    case (0)
       call DipoleFF(QSquared, GEp, GMp, GEn, GMn)
       MA=MA03

    case (1)
       ! hep-ex/0308005 - BBA2003 formfactors     valid up to Q**2= 6 GeV**2
       call BBA2003(QSquared, GEp, GMp, GEn, GMn)
       MA=MA03

    case (2)
       ! hep-ex/0602017 v3 (2006) -- BBBA2005 formfactors
       call BBBA2005(tau, GEp, GMp, GEn, GMn)
       MA=MA03  !no refitting has been done by Bodek et al compared to the BBA03 fit

    case (3)
       ! arXiv:0708.1827 [hep-ex] and arXiv:0708.1946 [hep-ex] -- BBA2007 formfactors
       call BBBA2007(tau, GEp, GMp, GEn, GMn)
       MA=MA07

    case default
       write(*,*) 'Wrong parametrization in  formfactors_QE',parametrization,'STOP!!!'
       stop
    end select

    if (useNonStandardMA) MA=MA_in  !overwrite MA with the value given in the jobcard

    if (initialState_charge==1) then
       if (present(GE)) GE=GEp
       if (present(GM)) GM=GMp
    else if (initialState_charge==0) then
       if (present(GE)) GE=GEn
       if (present(GM)) GM=GMn
    end if

    !*** Get Pauli and Dirac form factors:
    select case (processID)

    case (EM,antiEM)
       if (present(fa)) fa=0.
       if (present(fp)) fp=0.

       if (initialState_charge.eq.proton) then
          F1 = ((GEp)+tau*(GMp))/(1.+tau)
          F2 = ((GMp) - (GEp))/(1.+tau)
       else if (initialState_charge.eq.neutron) then
          F1 = ((GEn)+tau*(GMn))/(1.+tau)
          F2 = ((GMn) - (GEn))/(1.+tau)
       else
          write(*,*) 'Error in formfactors_QE! Strange initialState_charge',initialState_charge
       end if


    case (CC,antiCC)
       if (.not.(present(fa).and.present(fp))) then
          write(*,*) 'Error in formfactors_QE! Fa or fp are not defined',present(fa),present(fp),processID
          stop
       end if

       if (initialState_charge.eq.proton) then
          ! can be nonzero for proton in u-channel
          F1 = ((GEp-GEn)+tau*(GMp-GMn))/(1.+tau) !0.
          F2 = ((GMp-GMn) - (GEp-GEn))/(1.+tau) !0.
          FA = gA/(1.+QSquared/MA**2)**2 !0.
          FP=  ((2.*mN**2)/(mPi**2. + QSquared))*FA !0.
          return
       else if (initialState_charge.eq.neutron) then
          F1 = ((GEp-GEn)+tau*(GMp-GMn))/(1.+tau)
          F2 = ((GMp-GMn) - (GEp-GEn))/(1.+tau)
          FA = gA/(1.+QSquared/MA**2)**2
          if (axialMonopole) FA = gA/(1.+QSquared/MA**2)*(lambda1**2/(lambda1**2+QSquared))*lambda2**4/(lambda2**4+QSquared**2)
          FP= ((2.*mN**2)/(mPi**2. + QSquared))*FA
       else
          write(*,*) 'Error in formfactors_QE! Strange initialState_charge',initialState_charge
       end if


    case (NC,antiNC)
       if (.not.(present(fa).and.present(fp))) then
          write(*,*) 'Error in formfactors_QE! Fa or fp are not defined',present(fa),present(fp),processID
          stop
       end if

       F10S=0.
       F20S=0.

       !fit I
       !MAS=1.012
       !F10S=0.53
       !F20S=-0.4
       !deltas=-0.21

       !fit II
       !MAS=1.049
       !F10S=0.
       !F20S=0.
       !deltas=-0.15

       !fit III
       !   MAS=1.00
       !   F10S=0.
       !   F20S=0.
       !   deltas=0.

       MVS=0.843
       MAS=MA

       FQS=1./((1.+tau)*(1.+QSquared/MVS**2)**2)
       F1S=F10S*QSquared*FQS
       F2S=F20S*FQS
       FAS=deltas/(1. + QSquared/MAS**2)**2

       F1p=(tau*GMp+GEp)/(1.+tau)
       F1n=(tau*GMn+GEn)/(1.+tau)
       F2p=(GMp-GEp)/(1.+tau)
       F2n=(GMn-GEn)/(1.+tau)

       if (initialState_charge.eq.proton) then
          tau3=1.  !proton
          F1=(0.5-2.*sinsthweinbg)*F1p-0.5*F1n - F1S/2.
          F2=(0.5-2.*sinsthweinbg)*F2p-0.5*F2n - F2S/2.
       else if (initialState_charge.eq.neutron) then
          tau3=-1.  !neutron
          F1=(0.5-2.*sinsthweinbg)*F1n-0.5*F1p - F1S/2.
          F2=(0.5-2.*sinsthweinbg)*F2n-0.5*F2p - F2S/2.
       else
          write(*,*) 'Error in formfactors_QE! Strange initialState_charge',initialState_charge
       end if
       FA = gA*tau3/(2.*(1.+QSquared/MA**2.)**2.)+FAS/2.
       if (axialMonopole)  &
            & FA = gA*tau3/2./(1.+QSquared/MA**2)*(lambda1**2/(lambda1**2+QSquared))*lambda2**4/(lambda2**4+QSquared**2)+FAS/2.
       FP=((2.*mN**2)/(mPi**2.+QSquared))*FA


    case default
       write(*,*) 'Error in formfactors_QE! Invalid process ID:', processID
    end select


  end subroutine formfactors_QE


  !****************************************************************************
  !****s* FF_QE_nucleonScattering/DipoleFF
  ! NAME
  ! subroutine DipoleFF(QSquared, GEp, GMp, GEn, GMn)
  !
  ! PURPOSE
  ! * Provides the Sachs form factors in dipole approximation.
  !
  ! INPUTS
  ! *  real, intent(in) :: QSquared     ! = Q^2 : virtuality of gauge boson in units of GeV**2
  !
  ! OUTPUT
  ! *  real, intent(out) :: GEp, GMp, GEn, GMn            ! Sachs Form factors
  !
  !****************************************************************************
  subroutine DipoleFF(QSquared, GEp, GMp, GEn, GMn)
    real, intent(in) :: QSquared
    real, intent(out) :: GEp, GMp, GEn, GMn ! Sachs form factors
    GEp = 1./(1.+QSquared/MV2)**2           ! dipole
    GMp = GEp * mup
    GEn = 0.
    GMn = GEp * mun
  end subroutine DipoleFF


  !****************************************************************************
  !****s* FF_QE_nucleonScattering/BBA2003
  ! NAME
  ! subroutine BBA2003(tau, GEp, GMp, GEn, GMn)
  !
  ! PURPOSE
  ! * Provides the BBA2003 formfactors according to hep-ex/0308005,
  ! valid up to Q**2= 6 GeV**2
  !
  ! INPUTS
  ! *  real, intent(in) ::    QSquared   ! = Q^2 : virtuality of gauge boson in units of GeV**2
  !
  ! OUTPUT
  ! *  real, intent(out) :: GEp, GMp, GEn, GMn            ! Sachs Form factors
  !
  !****************************************************************************
  subroutine BBA2003(QSquared, GEp, GMp, GEn, GMn)
    use constants, only: mN
    real, intent(in)  :: QSquared
    real, intent(out) :: GEp, GMp, GEn, GMn
    real, parameter :: mvecsq=0.843**2
    real :: tau
    tau = QSquared/(4.*mN**2)
    GEp=1./(1.+3.253*QSquared+1.422*QSquared**2+0.08582*QSquared**3+ &
         0.3318*QSquared**4-0.09371*QSquared**5+0.01076*QSquared**6)
    GMp=mup/(1.+3.104*QSquared+1.428*QSquared**2+0.1112*QSquared**3- &
         0.006981*QSquared**4+0.0003705*QSquared**5-7.063e-6*QSquared**6)
    GEn=(-mun*0.942*tau)/(1.+4.61*tau)*(1.+QSquared/mvecsq)**(-2)       !parametrization of Krutov (hep-ph/0202183)
    GMn=mun/(1.+3.043*QSquared+0.8548*QSquared**2+0.6806*QSquared**3- &
         0.1287*QSquared**4+0.008912*QSquared**5)
  end subroutine BBA2003


  !****************************************************************************
  !****s* FF_QE_nucleonScattering/BBBA2005
  ! NAME
  ! subroutine BBBA2005(tau, GEp, GMp, GEn, GMn)
  !
  ! PURPOSE
  ! * Provides the BBBA2005 formfactors according to hep-ex/0602017 v3 (2006)
  !
  ! INPUTS
  ! *  real, intent(in) ::    tau
  !
  ! OUTPUT
  ! *  real, intent(out) :: GEp, GMp, GEn, GMn            ! Sachs Form factors
  !
  !****************************************************************************
  subroutine BBBA2005(tau, GEp, GMp, GEn, GMn)
    real, intent(in)  :: tau
    real, intent(out) :: GEp, GMp, GEn, GMn

    real, dimension(0:2,1:4) :: a
    real, dimension(1:4,1:4) :: b

    real, dimension(0:4) :: tauField
    integer :: i

    do i=0,4
       tauField(i)=tau**i
    end do

    a(:,1)=(/1.,-0.0578 ,0.  /)
    a(:,2)=(/1., 0.15   ,0.  /)
    a(:,3)=(/0., 1.25   ,1.3 /)
    a(:,4)=(/1., 1.81   ,0.  /)

    b(:,1)=(/11.1 ,   13.6,   33.  ,   0. /)
    b(:,2)=(/11.1 ,   19.6,    7.54,   0. /)
    b(:,3)=(/-9.86,  305. , -758.  , 802. /)
    b(:,4)=(/14.1 ,   20.7,   68.7 ,   0. /)

    GEP=Dot_Product(a(:,1),tauField(0:2))/(1.+Dot_Product(b(:,1),tauField(1:4)))
    GMP=Dot_Product(a(:,2),tauField(0:2))/(1.+Dot_Product(b(:,2),tauField(1:4)))
    GEN=Dot_Product(a(:,3),tauField(0:2))/(1.+Dot_Product(b(:,3),tauField(1:4)))
    GMN=Dot_Product(a(:,4),tauField(0:2))/(1.+Dot_Product(b(:,4),tauField(1:4)))

    GMP=GMP*mup
    GMn=GMn*mun

  end subroutine BBBA2005


  !****************************************************************************
  !****s* FF_QE_nucleonScattering/BBBA2007
  ! NAME
  ! subroutine BBBA2007(tau, GEp, GMp, GEn, GMn)
  !
  ! PURPOSE
  ! * Provides the BBBA2007 formfactors according
  !   to arXiv:0708.1827 [hep-ex] and arXiv:0708.1946 [hep-ex]
  ! * c++ file download from www.pas.rochester.edu/~bodek/FF/
  !
  ! INPUTS
  ! *  real, intent(in) ::    tau
  !
  ! OUTPUT
  ! *  real, intent(out) :: GEp, GMp, GEn, GMn            ! Sachs Form factors
  !
  !****************************************************************************
  subroutine BBBA2007(tau, GEp, GMp, GEn, GMn)
    real, intent(in)  :: tau
    real, intent(out) :: GEp, GMp, GEn, GMn

    real,dimension(2,4) :: kellyPara !first dimension: Gep, Gmp
    real,dimension(0:6) :: nodes
    real,dimension(2) :: GKelly
    real,dimension(4) :: GLagrange,part
    real, dimension(4,0:6) :: lagrangeParameter !first dimension: Gep,Gmp,Gmn,Gen, we use their fit with d/u=0.0
    integer :: i,j
    real :: xi

    !Kelly parametrization
    kellyPara(1,:)=(/ -0.24, 10.98, 12.82, 21.97 /)
    kellyPara(2,:)=(/ 0.1717, 11.26, 19.32, 8.33 /)
    GKelly=(1+kellyPara(:,1)*tau)/(1+kellyPara(:,2)*tau+kellyPara(:,3)*tau**2+kellyPara(:,4)*tau**3)


    !Nachtman variable xi
    xi=0.
    if (tau.gt.0.) xi=2./(1.+sqrt(1.+1/tau))


    !Lagrange parameterization
    lagrangeParameter(1,:)=(/ 1., .9927, .9898, .9975, .9812, .9340, 1. /)
    lagrangeParameter(2,:)=(/ 1., 1.0011, .9992, .9974, 1.0010, 1.0003,1. /)
    lagrangeParameter(3,:)=(/ 1., .9958, .9877, 1.0193, 1.0350, .9164, .7300 /)
    lagrangeParameter(4,:)=(/1., 1.1011, 1.1392, 1.0203, 1.1093, 1.5429, 0.9706/)

    GLagrange=0.
    nodes= (/ 0., 1./6., 2./6., 3./6., 4./6., 5./6., 1./)
    do i=0,6
       part = lagrangeParameter(:,i)
       do j=0,6
          if (i.ne.j) then
             part = part*(xi - nodes(j))/(nodes(i)-nodes(j))
          end if
       end do
       GLagrange=GLagrange+part
    end do


    !now: BBBA2007 form factors
    Gep=GLagrange(1)*GKelly(1)
    Gmp=GLagrange(2)*GKelly(2)
    Gmn=GLagrange(3)*Gmp
    Gen=GLagrange(4)*Gep*(1.7*tau)/(1+3.3*tau)  !last factor: Galster factor see Nucl Phys B32 221 (1971)

    GMP=GMP*mup
    GMn=GMn*mun

  end subroutine BBBA2007

end module FF_QE_nucleonScattering
