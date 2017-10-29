!******************************************************************************
!****m* /hypNuc_hypNuc
! NOTES
! This module includes the calculation of
! (a) Lambda(Sigma) nucleon --> Lambda(Sigma) nucleon (S=-1 channels)
! (b) Xi nucleon --> Lambda(Sigma) Lambda(Sigma) (S=-2 channels)
! cross sections including isospin exchange.
! The properties of all final channels are calculated.
!******************************************************************************
module hypNuc_hypNuc

  use preEventDefinition
  implicit none

  private

  public :: hypNuc_hypNuc_Main
  public :: get_Channels_YN

  ! ID's & charges of final states.
  ! 1-index: Number of isospin channels
  ! 2-index: number of particle in 2-body final state (Baryon,Hyperon)
  Integer, dimension(1:4,1:2), SAVE :: IdsOut,ChargesOut
  ! YN->YN Xsection (in units of mb)

  !Model switch for Xi nucleon xsections:
  Integer, Save :: XiN_Set = 2

contains

  !****************************************************************************
  !****f* hypNuc_hypNuc/hypNuc_hypNuc_Main
  ! NAME
  ! function hypNuc_hypNuc_Main (srts, teilchenIN) result (sigma_yn)
  ! PURPOSE
  ! * This routine prepares the final state for each isospin channel.
  ! INPUTS
  ! * type(particle), dimension (1:2), intent(in) :: teilchenIN  -- incoming particles
  ! * real,                            intent(in) :: srts        -- sqrt(s) [GeV]
  ! OUTPUT
  ! * see global variables of this module.
  ! NOTES
  ! Included reactions:
  ! * Lambda Nucleon -> Lambda Nucleon
  ! * Lambda Nucleon -> Sigma  Nucleon
  ! * Sigma  Nucleon -> Sigma  Nucleon
  ! * Sigma  Nucleon -> Lambda Nucleon
  ! * Xi     Nucleon -> Xi     Nucleon
  ! * Xi     Nucleon -> Lambda Lambda
  ! * Xi     Nucleon -> Lambda Sigma
  !****************************************************************************
  function hypNuc_hypNuc_Main (srts, partIn) result (sigma_yn)
    use particleDefinition
    use particleProperties, only: hadron
    use IdTable, only: nucleon,Lambda,SigmaResonance,Xi
    use random, only: rn
    !-----------------------------------------------------------------------------
    real,           intent(in)                  :: srts
    type(particle), intent(in), dimension (1:2) :: partIn
    real, dimension(1:4) :: sigma_yn
    !-----------------------------------------------------------------------------
    integer :: inuc, ihyp
    real    :: M_Nuc,M_Hyp

    logical, save :: hypNucInit = .true.

    !-----------------------------------------------------------------------------
    ! initialization:
    !-----------------------------------------------------------------------------
    sigma_yn(:) = 0.
    IdsOut(:,:) = 0
    ChargesOut(:,:) = 9999

    if (hypNucInit) then
       call get_XiN_Set
       hypNucInit=.false.
    end if
    !-----------------------------------------------------------------------------
    ! Check input: only Lambda/Sigma/Xi + Nucleon scattering!
    !-----------------------------------------------------------------------------
    if ( hadron(partIn(1)%ID)%strangeness==0) then
       inuc = 1
    else
       inuc = 2
    end if
    if (partIn(inuc)%ID /= nucleon) return  ! No nucleon resonances!
    ihyp = partIn(3-inuc)%ID
    if ( (ihyp > SigmaResonance) .and. (ihyp /= Xi) ) return
    !-----------------------------------------------------------------------------
    M_Nuc = partIn(inuc)%mass
    M_Hyp = partIn(3-inuc)%mass
    !-----------------------------------------------------------------------------
    ! elastic channel:
    IdsOut    (1,1:2) = partIn(1:2)%ID
    ChargesOut(1,1:2) = partIn(1:2)%charge
    !-----------------------------------------------------------------------------
    ! Prepare all isospin channels and get the PreEvent-vector:
    !-----------------------------------------------------------------------------
    select case (sum(partIn(:)%ID))

       case ((Lambda+nucleon))

          ! Lambda proton:
          if (partIn(inuc)%charge==1) then

             ! Lambda proton-->Lambda proton:
             sigma_yn(1) = xsectionYN(srts,M_Nuc,M_Hyp,1)

             ! Lambda proton-->Sigma^+ neutron:
             sigma_yn(2) = xsectionYN(srts,M_Nuc,M_Hyp,2)
             IdsOut(2,inuc) = nucleon
             IdsOut(2,3-inuc) = SigmaResonance
             ChargesOut(2,inuc) = 0
             ChargesOut(2,3-inuc) = 1

             ! Lambda proton-->Sigma^0 proton:
             sigma_yn(3) = xsectionYN(srts,M_Nuc,M_Hyp,3)
             IdsOut(3,inuc) = nucleon
             IdsOut(3,3-inuc) = SigmaResonance
             ChargesOut(3,inuc) = 1
             ChargesOut(3,3-inuc) = 0

          ! Lambda neutron:
          else if (partIn(inuc)%charge==0) then

             ! Lambda neutron-->Lambda neutron:
             sigma_yn(1) = xsectionYN(srts,M_Nuc,M_Hyp,4)

             ! Lambda neutron-->Sigma^- proton:
             sigma_yn(2) = xsectionYN(srts,M_Nuc,M_Hyp,5)
             IdsOut(2,inuc) = nucleon
             IdsOut(2,3-inuc) = SigmaResonance
             ChargesOut(2,inuc) = 1
             ChargesOut(2,3-inuc) = -1

             ! Lambda neutron-->Sigma^0 neutron:
             sigma_yn(3) = xsectionYN(srts,M_Nuc,M_Hyp,6)
             IdsOut(3,inuc) = nucleon
             IdsOut(3,3-inuc) = SigmaResonance
             ChargesOut(3,inuc) = 0
             ChargesOut(3,3-inuc) = 0

          end if

       case ((SigmaResonance+nucleon))

          ! Sigma^0 proton:
          if ( partIn(3-inuc)%charge==0 .and. partIn(inuc)%charge==1 ) then
             ! Sigma^0 proton-->Sigma^0 proton:
             sigma_yn(1) = xsectionYN(srts,M_Nuc,M_Hyp,7)

             ! Sigma^0 proton-->Lambda proton:
             sigma_yn(2) = xsectionYN(srts,M_Nuc,M_Hyp,8)
             IdsOut(2,inuc) = nucleon
             IdsOut(2,3-inuc) = Lambda
             ChargesOut(2,inuc) = 1
             ChargesOut(2,3-inuc) = 0

             ! Sigma^0 proton-->Sigma^+ neutron:
             sigma_yn(3) = xsectionYN(srts,M_Nuc,M_Hyp,9)
             IdsOut(3,inuc) = nucleon
             IdsOut(3,3-inuc) = SigmaResonance
             ChargesOut(3,inuc) = 0
             ChargesOut(3,3-inuc) = 1

          ! Sigma^0 neutron:
          else if ( partIn(3-inuc)%charge==0 .and. partIn(inuc)%charge==0 ) then
             ! Sigma^0 neutron-->Sigma^0 neutron:
             sigma_yn(1) = xsectionYN(srts,M_Nuc,M_Hyp,10)

             ! Sigma^0 neutron-->Sigma^- proton:
             sigma_yn(2) = xsectionYN(srts,M_Nuc,M_Hyp,11)
             IdsOut(2,inuc) = nucleon
             IdsOut(2,3-inuc) = SigmaResonance
             ChargesOut(2,inuc) = 1
             ChargesOut(2,3-inuc) = -1

             ! Sigma^0 neutron-->Lambda neutron:
             sigma_yn(3) = xsectionYN(srts,M_Nuc,M_Hyp,12)
             IdsOut(3,inuc) = nucleon
             IdsOut(3,3-inuc) = Lambda
             ChargesOut(3,inuc) = 0
             ChargesOut(3,3-inuc) = 0

          ! Sigma^- proton:
          else if ( partIn(3-inuc)%charge==-1 .and. partIn(inuc)%charge==1 ) then
             ! Sigma^- proton-->Sigma^- proton:
             sigma_yn(1) = xsectionYN(srts,M_Nuc,M_Hyp,13)

             ! Sigma^- proton-->Lambda neutron:
             sigma_yn(2) = xsectionYN(srts,M_Nuc,M_Hyp,14)
             IdsOut(2,inuc) = nucleon
             IdsOut(2,3-inuc) = Lambda
             ChargesOut(2,inuc) = 0
             ChargesOut(2,3-inuc) = 0

             ! Sigma^- proton-->Sigma^0 neutron:
             sigma_yn(3) = xsectionYN(srts,M_Nuc,M_Hyp,15)
             IdsOut(3,inuc) = nucleon
             IdsOut(3,3-inuc) = SigmaResonance
             ChargesOut(3,inuc) = 0
             ChargesOut(3,3-inuc) = 0

          ! Sigma^- neutron:
          else if ( partIn(3-inuc)%charge==-1 .and. partIn(inuc)%charge==0 ) then
             ! Sigma^- neutron-->Sigma^- neutron:
             sigma_yn(1) = xsectionYN(srts,M_Nuc,M_Hyp,16)

          ! Sigma^+ proton:
          else if ( partIn(3-inuc)%charge==1 .and. partIn(inuc)%charge==1 ) then
             ! Sigma^+ proton-->Sigma^+ proton:
             sigma_yn(1) = xsectionYN(srts,M_Nuc,M_Hyp,17)

          ! Sigma^+ neutron:
          else if ( partIn(3-inuc)%charge==1 .and. partIn(inuc)%charge==0 ) then
             ! Sigma^+ neutron-->Sigma^+ neutron:
             sigma_yn(1) = xsectionYN(srts,M_Nuc,M_Hyp,18)

             ! Sigma^+ neutron-->Lambda proton:
             sigma_yn(2) = xsectionYN(srts,M_Nuc,M_Hyp,19)
             IdsOut(2,inuc) = nucleon
             IdsOut(2,3-inuc) = Lambda
             ChargesOut(2,inuc) = 1
             ChargesOut(2,3-inuc) = 0

             ! Sigma^+ neutron-->Sigma^0 proton:
             sigma_yn(3) = xsectionYN(srts,M_Nuc,M_Hyp,20)
             IdsOut(3,inuc) = nucleon
             IdsOut(3,3-inuc) = SigmaResonance
             ChargesOut(3,inuc) = 1
             ChargesOut(3,3-inuc) = 0

          end if

       case ((Xi+nucleon))
          ! Xi^- p & Xi^0 n (I=0) --> elastic & Lambda Lambda
          if ( ( partIn(3-inuc)%charge==-1 .and. partIn(inuc)%charge==1 ) .or. &
               ( partIn(3-inuc)%charge==0  .and. partIn(inuc)%charge==0 ) ) then

             if (XiN_Set == 1) then !Rijken/Yamamoto, nucl-th/0608074

                ! Xi^- p & Xi^0 n (I=0) --> elastic
                sigma_yn(1) = xsectionYN(srts,M_Nuc,M_Hyp,21)

                ! Xi^- p & Xi^0 n (I=0) --> Lambda Lambda
                sigma_yn(2) = xsectionYN(srts,M_Nuc,M_Hyp,22)
                IdsOut(2,inuc) = Lambda
                IdsOut(2,3-inuc) = Lambda
                ChargesOut(2,inuc) = 0
                ChargesOut(2,3-inuc) = 0

             else if (XiN_Set == 2) then !Fujiwara et al., PRC64, 054001

                ! Xi^- p & Xi^0 n (I=0) --> elastic
                sigma_yn(1) = xsectionYN(srts,M_Nuc,M_Hyp,25)

                ! Xi^- p & Xi^0 n (I=0) --> Lambda Lambda
                sigma_yn(2) = xsectionYN(srts,M_Nuc,M_Hyp,26)
                IdsOut(2,inuc) = Lambda
                IdsOut(2,3-inuc) = Lambda
                ChargesOut(2,inuc) = 0
                ChargesOut(2,3-inuc) = 0

                ! Xi^- p & Xi^0 n (I=0) --> Lambda Sigma^0
                sigma_yn(3) = xsectionYN(srts,M_Nuc,M_Hyp,27)
                IdsOut(3,inuc) = Lambda
                IdsOut(3,3-inuc) = SigmaResonance
                ChargesOut(3,inuc) = 0
                ChargesOut(3,3-inuc) = 0

                ! Xi^- p (I=0) --> Xi^0 n
                if ( partIn(3-inuc)%charge==-1 .and. partIn(inuc)%charge==1 ) then
                   sigma_yn(4) = xsectionYN(srts,M_Nuc,M_Hyp,28)
                   IdsOut(4,inuc) = Xi
                   IdsOut(4,3-inuc) = nucleon
                   ChargesOut(4,inuc) = 0
                   ChargesOut(4,3-inuc) = 0
                end if

                ! Xi^0 n (I=0) --> Xi^- p (detailed balance)
                if ( partIn(3-inuc)%charge==0 .and. partIn(inuc)%charge==0 ) then
                   sigma_yn(4) = xsectionYN(srts,M_Nuc,M_Hyp,29)
                   IdsOut(4,inuc) = Xi
                   IdsOut(4,3-inuc) = nucleon
                   ChargesOut(4,inuc) = -1
                   ChargesOut(4,3-inuc) = 1
                end if

             else

                write(*,*) 'module HypBar_HypBar.f90/function hypNuc_hypNuc_Main: '
                write(*,*) 'wrong input for variable XiN_Set! XiN_Set = ',XiN_Set
                write(*,*) 'STOPPING execution'
                STOP

             end if

             ! Xi^- + n --> Xi^- n & Lambda Sigma^-
             ! (only for YN_Set=2, Fujiwara et al., PRC64, 054001)
             else if ( (XiN_Set==2) .and. &
                  & (partIn(3-inuc)%charge==-1 .and. partIn(inuc)%charge==0) ) then

                ! Xi^- + n --> Xi^- n
                sigma_yn(1) = xsectionYN(srts,M_Nuc,M_Hyp,30)

                ! Xi^- + n --> Lambda Sigma^-
                sigma_yn(2) = xsectionYN(srts,M_Nuc,M_Hyp,31)
                IdsOut(2,inuc) = Lambda
                IdsOut(2,3-inuc) = SigmaResonance
                ChargesOut(2,inuc) = 0
                ChargesOut(2,3-inuc) = -1

          ! Xi + Nucleon (I=1) --> Xi N & Lambda Sigma
          ! (for Xi^0 p channel only one model exists (Rijken/Yamamoto, nucl-th/0608074))
          else if ( partIn(3-inuc)%charge==0 .and. partIn(inuc)%charge==1 ) then
             ! Xi + Nucleon (I=1) --> Xi + Nucleon
             sigma_yn(1) = xsectionYN(srts,M_Nuc,M_Hyp,23)

             ! Xi + Nucleon (I=1) --> Lambda Sigma
             sigma_yn(2) = xsectionYN(srts,M_Nuc,M_Hyp,24)
             if (rn() < 0.5) then
                IdsOut(2,inuc) = Lambda
                IdsOut(2,3-inuc) = SigmaResonance
                ChargesOut(2,inuc) = 0
                ChargesOut(2,3-inuc) = 1
             else
                IdsOut(2,3-inuc) = Lambda
                IdsOut(2,inuc) = SigmaResonance
                ChargesOut(2,3-inuc) = 0
                ChargesOut(2,inuc) = 1
             end if

          end if

       case default

          write(*,*) 'HypBar_HypBar_Main: undefined channel!', partIn%ID
          STOP

       end select

  end function hypNuc_hypNuc_Main


  !****************************************************************************
  !****s* hypNuc_hypNuc/get_XiN_Set
  ! NAME
  ! returns the value of XiN_Set
  ! PURPOSE
  ! Reads from Namelist 'XiN_input' the value for the variable XiN_Set, which
  ! determines the model usage for the Xi nucleon cross sections.
  !****************************************************************************
  subroutine get_XiN_Set
    use output
    integer :: ios

    NAMELIST /XiN_input/XiN_Set

       call Write_ReadingInput('XiN_input',0)
       rewind(5)
       read(5,nml=XiN_input,iostat=ios)
       call Write_ReadingInput('XiN_input',0,ios)
       write(*,*) ' Set XiN_Set to', XiN_Set,'.'
       call Write_ReadingInput('SMM_input',1)

  end subroutine get_XiN_Set



  !****************************************************************************
  !****f* hypNuc_hypNuc/get_Channels_YN
  ! NAME
  ! function get_Channels_YN (InChan) result (ch)
  ! PURPOSE
  ! Returns channel infos (IDs and charges) as a PreEvent vector.
  ! INPUTS
  ! integer :: InChan -- Number of final channel
  !****************************************************************************
  function get_Channels_YN (InChan) result (ch)
    integer, intent(in) :: InChan
    type(preEvent), dimension(1:2) :: ch

    ch(1:2)%ID = IdsOut(InChan,1:2)
    ch(1:2)%charge = ChargesOut(InChan,1:2)

  end function get_Channels_YN


  !****************************************************************************
  !****f* hypNuc_hypNuc/xsectionYN
  ! NAME
  ! real function xsectionYN (srts, M_Nuc, M_Hyp, ichan)
  ! PURPOSE
  ! * This function calculates the hyperon-baryon -> hyperon baryon cross sections
  !   for a particular isospin channel "ichan".
  ! * The cross sections are given in mb.
  ! INPUTS
  ! * integer, intent(in) :: ichan  -- initial channel
  ! * real,    intent(in) :: M_Nuc,M_Hyp -- masses of incoming particles
  ! * real,    intent(in) :: srts -- sqrt(s) [GeV]
  ! NOTES
  ! * Fits to experimentally known cross sections.
  ! * Not for all channels the XS were experimentally known. The XS's for the unknown
  !   channels were calculated using detailed balance or charge symmetry.
  ! * For the important channels, where a Lambda is converted into a Sigma
  !   (or a Sigma to a Lambda), the Xsections are known.
  ! * References for exp. data on hyperon+nucleon scattering:
  !   Landolt-Boernstein, New Series I/12b, p. 323: "Hyperon induced reactions"
  !   M.M. Nagels, T.A. Rijken, J.J. De Swart, Annals of Physics 79 (1973) 338.
  !   T.A. Rijken, Y. Yamamoto, Phys. Rev. C73 (2006) 044008.
  ! * Xi+N scattering cross sections: no exp. data available here; fits to theoretical
  !   calculations: ESC04a model,
  !   T.A. Rijken, Y. Yamamoto, nucl-th/0608074 (30.08.2006), Tables IX,X,XI & XII.
  ! * Xi+N scattering cross sections UPDATE: theoretical calculations from
  !   Fujiwara et al., PRC64, 054001 including exclusive more channels.
  !****************************************************************************
  real function xsectionYN (srts, M_Nuc, M_Hyp, ichan)

    use particleProperties, only: hadron
    use IdTable, only: Lambda, SigmaResonance, Xi
    use twoBodyTools, only: pCM
    use constants, only: mN

    real,    intent(in) :: srts, M_Nuc, M_Hyp
    integer, intent(in) :: ichan

    real :: plab,p_ab,p_cd,balFac,pmev, x

    if (ichan<1 .or. ichan>31) then
       write(*,*) 'hypNuc_hypNuc/xsection: ichan not defined! ichan = ',ichan
       STOP
    end if

    ! hyperon as beam-particle
    plab = sqrt( (srts**2-(m_hyp-m_nuc)**2)*(srts**2-(m_hyp+m_nuc)**2) )/(2.*m_nuc)
    pmev = plab*1000.
    x    = pmev

    p_ab = pCM (srts, M_Nuc, M_Hyp)  ! needed for detailed balance

    xsectionYN = 0.

    select case (ichan)

    ! Lambda/Sigma + Nucleon channels:

    case (1,4)

       ! Lambda + N --> Lambda + N   (parametrization by Oliver Arnold, following Rijken et al, PRC 59, 21-40, 1999)
       if (plab < 0.4) then
         xsectionYN = 203.56*exp(-14.47*plab**2) + 253.88*exp(-76.19*plab**2)
       else
          xsectionYN = 14.4*plab**(-0.12)
       end if

    case (2,3,5,6,8,12)

       ! Lambda + N --> Sigma0 + N  (ichan==3,6)
       !       if (plab < 0.436) then
       if (plab < 0.8) then
          xsectionYN = 30.*plab**(4.9)
       else if (plab > 0.8 .and. plab < 1.206) then
          xsectionYN = 5.7*plab**(-2.5) !orig.
       else if (plab > 1.206) then
          xsectionYN = min (19., 9.9456-10.254*plab+4.1135*plab**2)
       end if

       if (ichan==2 .or. ichan==5) then
         ! Lambda + Proton  --> Sigma+ + Neutron  (ichan==2)
         ! Lambda + Neutron --> Sigma- + Proton   (ichan==5)
         xsectionYN = 0.5*xsectionYN
       else if (ichan==8 .or. ichan==12) then
         ! Sigma0 + N --> Lambda + N
         p_cd = pCM(srts,hadron(lambda)%mass,mN)
         if (p_ab<1E-10) then
           write(*,*) 'WARNING: pInitial is zero in xsectionYN', p_ab
           balFac= 0.
         else
           balFac= (p_cd/p_ab)**2
         end if
         xsectionYN = xsectionYN*balFac ! detailed balance
       end if

    case (7,10,13,18)

       ! Sigma0 + N --> Sigma0 + N  (ichan==7,10)
       ! Sigma- + p --> Sigma- + p  (ichan==13)
       ! Sigma+ + n --> Sigma+ + n  (ichan==18)
       xsectionYN =  13.5*plab**(-1.25)
       if (plab>1.0) xsectionYN =  13.5

    case (9,11,15,20)

       ! Sigma- + Proton  --> Sigma0 + Neutron (ichan==15)
       ! Sigma+ + Neutron --> Sigma0 + Proton  (ichan==20)
       xsectionYN = 13.5*plab**(-1.25)
       if (plab > 1.0) xsectionYN = 13.5

       if (ichan==9 .or. ichan==11) then
         ! Sigma0 + Proton  --> Sigma+ + Neutron (ichan==9)
         ! Sigma0 + Neutron --> Sigma- + Proton  (ichan==11)
         p_cd = pCM(srts,hadron(SigmaResonance)%mass,mN)
         if (p_ab<1E-10) then
           write(*,*) 'WARNING: pInitial is zero in xsectionYN', p_ab
           balFac= 0.
         else
           balFac= (p_cd/p_ab)**2
         end if
         xsectionYN = xsectionYN*balFac ! detailed balance
       end if

    case (14,19)

       ! Sigma- + Proton  --> Lambda + Neutron (ichan==14)
       ! Sigma+ + Neutron --> Lambda + Proton  (ichan==19)
       xsectionYN = 13.2*plab**(-1.18)
       if (plab > 1.0) xsectionYN = 13.2

    case (16,17)

       ! Sigma- + Neutron --> Sigma- + Neutron (ichan==16)
       ! Sigma+ + Proton  --> Sigma+ + Proton  (ichan==17)
       xsectionYN = 38.*plab**(-0.62)

    ! Xi Nucleon channels according
    ! T.A. Rijken, Y. Yamamoto, nucl-th/0608074 (30.08.2006), Tables IX,X,XI & XII.
    case (21)
       ! Xi N --> Xi N (I=0)
       xsectionYN = 17.3886*exp(-0.01*pmev)+2572.95*pmev**(-1.18381)

    case (22)
       ! Xi N --> Lambda Lambda (I=0)
       xsectionYN = 417.996*exp(-0.00567813*pmev) - 340.722*pmev**(-0.359183) + 31.*pmev**0.02

    case (23)
       ! Xi N --> Xi N (I=1)
       xsectionYN = 328.083*exp(-0.00418763*pmev)+0.00602*pmev

    case (24)
       ! Xi N --> Sigma Lambda (I=1)
       if (pmev > 590.) xsectionYN = 4.626*(pmev/590.+1.)**(-0.604771)

    ! Xi Nucleon channels according Fujiwara et al., PRC64, 054001:

    !Xi^- p --> Xi^- p, Xi^0 n --> Xi^0 n :
    case (25)
       if (x .lt. 208.23) then
          xsectionYN = 43179.9/(x**1.4509+293.093)
       else if (x.gt.400. .and. x .le. 570.) then
          xsectionYN = -24.9623*x**0.165388 + 58.*x**0.0515513
       else if (x .gt. 570. .and. x .le. 600.5) then
          xsectionYN = 50./((x-589.5)**2+25.) + 11.
       else
          xsectionYN = 865.604*x**(-0.780425) + 0.46935*x**0.357044
       end if
       if (x.gt.1000.) xsectionYN = 9.4737488

       !write(500,*) x,xsectionYN

    ! Xi^- p & Xi^0 n (I=0) --> Lambda Lambda
    case (26)
       if (x .le. 200.) then
          xsectionYN = 555307./(x**2.25814+6206.04)
       else if (x.gt.200. .and. x.lt.1000.) then
          xsectionYN = 6.06073*x**0.69767 - 5.74077*x**0.705059
       else
          xsectionYN = 2.39538574
       end if

    ! Xi^- p & Xi^0 n (I=0) --> Lambda Sigma^0
    case (27)
       if (x .le.586.2) then
          xsectionYN=0.0
       else if (x.gt.586.2 .and. x.lt.1000.) then
          xsectionYN = 18.1576*(x/586.2-1.)**0.296965*(586.2/x)**2.31161
       else
          xsectionYN = 4.76379824
       end if

    ! Xi^- p (I=0) --> Xi^0 n + detailed balance
    case (28, 29)
       if (x.le.233.333) then
          xsectionYN = 302.336/( 0.0615418*x + 0.177857 )
       else if (x.gt.233.333 .and. x.le.343.137) then
          xsectionYN = 186597./( 9.72237*x + 6711.97 )
       else if (x.gt.343.137 .and. x.le.600.) then
          xsectionYN = 59.8056/( 0.00847445*x + 0.265243 )
       else if (x.gt.600. .and. x.lt.1000.) then
          xsectionYN = 1136.84/( 0.120506*x + 23.3462 )
       else
          xsectionYN = 7.90283298
       end if
       if (ichan==29) then
          p_cd = pCM(srts,hadron(Xi)%mass,mN)
          if (p_ab<1E-10) then
             write(*,*) 'WARNING: pInitial is zero in xsectionYN', p_ab
             balFac= 0.
          else
             balFac= (p_cd/p_ab)**2
          end if
          xsectionYN = xsectionYN*balFac ! detailed balance
       end if

    ! Xi^- + n --> Xi^- n
    case (30)
       if (x.le.571.) then
          xsectionYN = 6.373005+0.873676*x**(1.65697)- &
           &  0.875998*x**(1.65655)
       else if (x.gt.571. .and. x.le.618.627) then
          xsectionYN = 15.5/((x-592.)**2-1.) + 12.0975
       else if (x.gt.618.627 .and. x.lt.1000.) then
          xsectionYN = -34.545 + 7.2751 * log(x)
       else
          xsectionYN = 15.7096138
       end if

    ! Xi^- + n --> Lambda Sigma^-
    case (31)
       if (x.le.587.93) then
          xsectionYN = 1.0e-10
       else if (x.gt.587.93 .and. x.lt.1000.) then
          xsectionYN = 25.*(x/587.93-1.)**(0.18) * (587.93/x)**1.5
       else
          xsectionYN = 10.5716934
       end if

    case default

       write(*,*) 'HypBar_HypBar/xsectionsYN: undefined channel! ichan = ', ichan
       STOP

    end select

  end function xsectionYN


end module hypNuc_hypNuc
