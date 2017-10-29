!******************************************************************************
!****m* /nu_MatrixElementRS
! NAME
! module nu_MatrixElementRS
!
! PURPOSE
! This module calculates the squared spin summed and averaged Rein-Sehgal
! matrix element for lepton induced resonance excitation.
!
!******************************************************************************


module nu_MatrixElementRS
  implicit none

  private

  logical, save :: initFlag=.true.

  real,save :: MVRS = 0.84
  real,save :: MARS = 1.032
  logical, save :: vector_FF_switch=.true.
  logical, save :: axial_FF_switch=.true.


  public :: nu_MaEl_RS


contains

  subroutine readInput_nu_MatrixElementRS
    use output

    integer :: ios

    NAMELIST /neutrino_matrixelement_RS/  MVRS,MARS,vector_FF_switch,axial_FF_switch
    call Write_ReadingInput('neutrino_matrixelement_RS',0)
    rewind(5)
    read(5,nml=neutrino_matrixelement_RS,IOSTAT=ios)
    call Write_ReadingInput("neutrino_matrixelement_RS",0,ios)
    call Write_ReadingInput('neutrino_matrixelement_RS',1)

    write(*,*) 'values for Rein-Sehgal-FF: MV=',MVRS,' MA=',MARS

    if (.not.vector_FF_switch) write(*,*) 'WARNING: vector FF are switched off!!!'
    if (.not.axial_FF_switch) write(*,*) 'WARNING: axial FF are switched off!!!'

  end subroutine readInput_nu_MatrixElementRS


  real function nu_MaEl_RS(process_ID,finalstate_ID,charge_in,k_in,k_out,p_in,p_out)
    use idtable, only: nucleon,delta, &
                        P11_1440,S11_1535,S11_1650,D13_1520,D13_1700,D15_1675,P11_1710,P13_1720,F15_1680,F17_1990,S31_1620, &
                        D33_1700,P31_1910,P33_1600,P33_1920,F35_1905,F37_1950
    use constants, only: alphaQED,coscab,GF,mN,pi,sinsthweinbg
    use ParticleProperties, only: hadron
    use minkowski, only: SP,abs4
    use neutrino_IDTable
    use leptonicID

    real, dimension(0:3), intent(in) :: k_in, k_out, p_in, p_out
    real:: s_in,t_in,Mi,Mf!,ml_in,ml_out
    integer, intent(in) :: finalstate_ID, charge_in
    integer, intent(in) :: process_ID

    real :: coupling=0.
    real :: Qs, M,epsilon,costheta
    integer :: n
    real :: enu, nu, elepton, qq, qsq, omega
    real :: root_half_omega, m_target,W, U, V
    real :: z, xi, MV, MA, ffcorr,b,c
    real :: ga,gv,lambda,tv,t,rv,r
    real :: s,ta,ra,rminus,rplus,tminus,tplus
    real :: fminus1,fminus3,fplus1,fplus3,f0plus,f0minus

    real :: at,bt

    real :: sigm,sigp,sign

    if (initFlag) then
       call readInput_nu_MatrixElementRS
       initFlag=.false.
    end if

    if (finalstate_ID.eq.nucleon) then
       write(*,*) 'finalstate_ID.eq.nucleon -> this should not happen -> stop'
       stop
    end if

    !note: better don't use Rein and Sehgal for tau neutrinos!!!

    nu_MaEl_RS=0.

    t_in=SP(k_in-k_out,k_in-k_out)
    s_in=SP(k_in+p_in,k_in+p_in)
    Mi=abs4(p_in)
    Mf=abs4(p_out)
    M=mN

    Qs=-t_in


    enu=(s_in-Mi**2)/(2.*Mi)
    nu=0.5*(Mf**2-Mi**2+QS)/Mi
    elepton=enu-nu
    qq=sqrt(nu**2+Qs)
    qsq=-Qs

    costheta=1.-Qs/2./enu/elepton
    epsilon=1./(1.+2.*(1.+nu**2/Qs)*(1.-costheta)/(1.+costheta))

    omega=1.05
    root_half_omega=sqrt(omega/2.)

    m_target=Mi  !MNucleon?
    w=Mf


    z=0.75

    xi = sinsthweinbg

    MV=MVRS
    MA=MARS


    n=2
    if (finalstate_ID.eq.delta) n=0
    if (finalstate_ID.eq.D13_1520.or.finalstate_ID.eq.S11_1535.or.finalstate_ID.eq.S31_1620.or.  &
         & finalstate_ID.eq.S11_1650.or.finalstate_ID.eq.D13_1700.or.finalstate_ID.eq.D15_1675  &
         & .or.finalstate_ID.eq.D33_1700) n=1


    !PCAC value of 1.2 in GA ??
    ffcorr=sqrt(1.+Qs/(4.*M**2))
    gv=ffcorr**(1-2*n)*(1./(1.+Qs/MV**2))**2
    ga=ffcorr**(1-2*n)*(1./(1.+Qs/MA**2))**2


    if (.not.vector_FF_switch) gv=0.
    if (.not.axial_FF_switch) ga=0.


    lambda=(m_target*qq)/w/root_half_omega
    tv=root_half_omega*gv/(3.*w)
    T=tv
    rv=sqrt(2.)*(m_target/w)*(((w+m_target)*qq)/((w+m_target)**2-qsq))*gv
    R=rv
    s=(-qsq/qq**2)*((3.*w*m_target+qsq-m_target**2)/(6.*m_target**2))*gv
    ta=(2./3.)*z*root_half_omega*(m_target/w)*(qq/((w+m_target)**2-qsq))*ga
    ra=(z*sqrt(2.)/(6.*w))*(w+m_target+(2*n*omega*w)/((w+m_target)**2-qsq) )*ga
    b=(z/(3.*w))*root_half_omega*(1.+(abs(w**2-m_target**2+qsq)/((w+m_target)**2-qsq)))*ga
    c=(z/(6.*m_target*qq))*(w**2-m_target**2+n*omega*(abs(w**2-m_target**2+qsq)/((w+m_target)**2-qsq)))*ga

    rminus=-(rv-ra)
    rplus=-(rv+ra)
    tminus=-(tv-ta)
    tplus=-(tv+ta)


    !note the different sign!!! in U and V in contrast to Rein Sehgal,
    !but in agreement with Nuance and Neugen, data and my former calculation
    if (process_ID.gt.0) then !neutrino
       U=(enu+elepton-qq)/(2.*enu)
       V=(enu+elepton+qq)/(2.*enu)
    else   !antineutrino
       V=(enu+elepton-qq)/(2.*enu)
       U=(enu+elepton+qq)/(2.*enu)
    end if


    select case (process_ID)
    case (NC,antiNC)
       coupling=GF**2/4.

       if (charge_in.eq.proton) nu_MaEl_RS=nu_MaEl_NCprot()
       if (charge_in.eq.neutron) nu_MaEl_RS=nu_MaEl_NCneut()

    case (CC)
       coupling=GF**2/4.*coscab**2

       if (process_ID.gt.0) then
          if (charge_in.eq.proton) then
             if (hadron(finalstate_ID)%isoSpinTimes2.eq.3) then
                nu_MaEl_RS=3.*nu_MaEl_CC()
             else if (hadron(finalstate_ID)%isoSpinTimes2.eq.1) then
                nu_MaEl_RS=0.
             end if
          end if
          if (charge_in.eq.neutron) nu_MaEl_RS=nu_MaEl_CC()
       else !antineutrino
          if (charge_in.eq.neutron) then
             if (hadron(finalstate_ID)%isoSpinTimes2.eq.3) then
                nu_MaEl_RS=3.*nu_MaEl_CC()
             else if (hadron(finalstate_ID)%isoSpinTimes2.eq.1) then
                nu_MaEl_RS=0.
             end if
          end if
          if (charge_in.eq.proton) nu_MaEl_RS=nu_MaEl_CC()
       end if


    case (EM)
     !  coupling=(4.*pi*alphaQED)**2/Qs**2*1./4.

       coupling=64.*pi**2*alphaQED**2*Mf**2/Qs/(1.-epsilon)

       if (charge_in.eq.proton) nu_MaEl_RS=nu_MaEl_EMprot()
       if (charge_in.eq.neutron) nu_MaEl_RS=nu_MaEl_EMneut()

    end select

    nu_MaEl_RS=coupling*nu_MaEl_RS


  contains

    real function nu_MaEl_NCprot()

      select case (finalstate_ID)

      case (delta)

         fMinus1 =  -sqrt(2.) * (rminus + 2 * xi * R)
         fPlus1  =   sqrt(2.) * (rplus  + 2 * xi * R)
         fMinus3 =  -Sqrt(6.) * (rminus + 2 * xi * R)
         fPlus3  =   Sqrt(6.) * (rplus  + 2 * xi * R)
         f0Minus =  2.*sqrt(2.) * c
         f0Plus  =   f0Minus

      case (S11_1535)

         at       = sqrt(3./2.) * (1.-2.*xi) * Lambda * S
         bt       = Sqrt(2./3.) * (Lambda * c - 3*b)

         fMinus1 =  Sqrt(3.) * (Tminus + 2*xi*T) + Sqrt(2./3.) * (Lambda * (Rminus + 3*xi*R))
         fPlus1  = -1.*Sqrt(3.) * (Tplus  + 2*xi*T) - Sqrt(2./3.) * (Lambda * (Rplus  + 3*xi*R))
         fMinus3 =  0.
         fPlus3  =  0.
         f0Minus = -at + bt
         f0Plus  =  at + bt


      case (D13_1520)

         at       = Sqrt(3.) * (1.-2*xi) * Lambda * S
         bt       = (2./Sqrt(3.)) * Lambda * c

         fMinus1 = sqrt(3./2.) * (Tminus + 2*xi*T) - sqrt(4./3.) * Lambda * (Rminus + 3*xi*R)
         fPlus1  = sqrt(3./2.) * (Tplus  + 2*xi*T) - sqrt(4./3.) * Lambda * (Rplus  + 3*xi*R)
         fMinus3 = 3./Sqrt(2.) * (Tminus + 2*xi*T)
         fPlus3  = 3./Sqrt(2.) * (Tplus  + 2*xi*T)
         f0Minus = -at + bt
         f0Plus  = -at - bt


      case (S11_1650)

         fMinus1 =  1./Sqrt(24.) * Lambda * Rminus
         fPlus1  = -1./Sqrt(24.) * Lambda * Rplus
         fMinus3 =  0.
         fPlus3  =  0.
         f0Minus = -1./Sqrt(6.) * (Lambda * c - 3*b)
         f0Plus  =  f0Minus


      case (D13_1700)

         fMinus1 =  1./Sqrt(120.) * Lambda * Rminus
         fPlus1  =  1./Sqrt(120.) * Lambda * Rplus
         fMinus3 =  3./Sqrt(40.)  * Lambda * Rminus
         fPlus3  =  3./Sqrt(40.)  * Lambda * Rplus
         f0Minus =  1./Sqrt(30.)  * Lambda * c
         f0Plus  =  -1.* f0Minus


      case (D15_1675)

         fMinus1 = -Sqrt(3./40.) * Lambda * Rminus
         fPlus1  =  Sqrt(3./40.) * Lambda * Rplus
         fMinus3 = -Sqrt(3./20.) * Lambda * Rminus
         fPlus3  =  Sqrt(3./20.) * Lambda * Rplus
         f0Minus =  Sqrt(3./10.) * Lambda * c
         f0Plus  =  f0Minus

      case (S31_1620)

         at       = sqrt(3./2.) * (1.-2.*xi) * Lambda * S
         bt       = 1./Sqrt(6.) * (Lambda * c - 3*b)

         fMinus1 =  Sqrt(3.) * (Tminus + 2*xi*T) - 1./Sqrt(6.) * Lambda * (Rminus + 2*xi*R)
         fPlus1  = -Sqrt(3.) * (Tplus  + 2*xi*T) + 1./Sqrt(6.) * Lambda * (Rplus  + 2*xi*R)
         fMinus3 =  0.
         fPlus3  =  0.
         f0Minus = -at-bt
         f0Plus  =  at-bt


      case (D33_1700)

         at       = Sqrt(3.) * (1.-2.*xi) * Lambda * S
         bt       = 1./Sqrt(3.) * Lambda * c

         fMinus1 = sqrt(3./2.) * (Tminus + 2*xi*T) + 1./Sqrt(3.) * Lambda * (Rminus + 2*xi*R)
         fPlus1  = sqrt(3./2.) * (Tplus  + 2*xi*T) + 1./Sqrt(3.) * Lambda * (Rplus  + 2*xi*R)
         fMinus3 = 3./Sqrt(2.) * (Tminus + 2*xi*T)
         fPlus3  = 3./Sqrt(2.) * (Tplus  + 2*xi*T)
         f0Minus = -at-bt
         f0Plus  = -at+bt


      case (P11_1440)

         at       = 0.25 * Sqrt(3.) * (1-4*xi) * Lambda**2 * S
         bt       = (5./12.)*Sqrt(3.) * (Lambda**2 * c - 2 * Lambda * b)

         fMinus1 = -(5./12.)*Sqrt(3.) * Lambda**2 * (Rminus + (12./5.)*xi*R)
         fPlus1  = -(5./12.)*Sqrt(3.) * Lambda**2 * (Rplus  + (12./5.)*xi*R)
         fMinus3 = 0.
         fPlus3  = 0.
         f0Minus = -at+bt
         f0Plus  = -at-bt


      case (P33_1600)

         fMinus1 =  1./Sqrt(6.) * Lambda**2 * (Rminus + 2*xi*R)
         fPlus1  = -1./Sqrt(6.) * Lambda**2 * (Rplus  + 2*xi*R)
         fMinus3 =  1./Sqrt(2.) * Lambda**2 * (Rminus + 2*xi*R)
         fPlus3  = -1./Sqrt(2.) * Lambda**2 * (Rplus  + 2*xi*R)
         f0Minus = -Sqrt(2./3.) * (Lambda**2 * c - 2 * Lambda * b)
         f0Plus  =  f0Minus


      case (P13_1720)

         at       = Sqrt(3./20.) * (1.-4*xi) * Lambda**2 * S
         bt       = Sqrt(5./12.) * (Lambda**2 * c - 5 * Lambda * b)

         fMinus1 = -Sqrt(27./40.) * Lambda * (Tminus + 4*xi*T) - Sqrt(5./12.) * Lambda**2 * (Rminus + (12./5.)*xi*R)
         fPlus1  =  Sqrt(27./40.) * Lambda * (Tplus  + 4*xi*T) + Sqrt(5./12.) * Lambda**2 * (Rplus  + (12./5.)*xi*R)
         fMinus3 =  3./Sqrt(40.)  * Lambda * (Tminus + 4*xi*T)
         fPlus3  = -3./Sqrt(40.)  * Lambda * (Tplus  + 4*xi*T)
         f0Minus =  at-bt
         f0Plus  = -at-bt


      case (F15_1680)

         at       = 3./Sqrt(40.) * (1.-4*xi)* Lambda**2 * S
         bt       = Sqrt(5./8.) * Lambda**2 * c

         fMinus1 = -3./Sqrt(20.) * Lambda * (Tminus + 4. * xi *T) + Sqrt(5./8.) * Lambda**2 * (Rminus + (12./5.) * xi * R)
         fPlus1  = -3./Sqrt(20.) * Lambda * (Tplus  + 4. * xi *T) + Sqrt(5./8.) * Lambda**2 * (Rplus  + (12./5.) * xi * R)
         fMinus3 = -Sqrt(18./20.) * Lambda * (Tminus + 4. * xi *T)
         fPlus3  = -Sqrt(18./20.) * Lambda * (Tplus  + 4. * xi *T)
         f0Minus =  at - bt
         f0Plus  =  at + bt


      case (P31_1910)

         fMinus1 = -1./Sqrt(15.) * Lambda**2 * (Rminus + 2*xi*R)
         fPlus1  = -1./Sqrt(15.) * Lambda**2 * (Rplus  + 2*xi*R)
         fMinus3 =  0.
         fPlus3  =  0.
         f0Minus = -Sqrt(4./15.) * (Lambda**2 * c - 5 * Lambda * b)
         f0Plus  = -1.* f0Minus


      case (P33_1920)

         fMinus1 =  1./Sqrt(15.) * Lambda**2 * (Rminus + 2*xi*R)
         fPlus1  = -1./Sqrt(15.) * Lambda**2 * (Rplus  + 2*xi*R)
         fMinus3 = -1./Sqrt(5.)  * Lambda**2 * (Rminus + 2*xi*R)
         fPlus3  =  1./Sqrt(5.)  * Lambda**2 * (Rplus  + 2*xi*R)
         f0Minus = -(2./Sqrt(15.)) * (Lambda**2 * c - 5 * Lambda * b)
         f0Plus  =  f0Minus


      case (F35_1905)

         fMinus1 =  1./Sqrt(35.)  * Lambda**2 * (Rminus + 2*xi*R)
         fPlus1  =  1./Sqrt(35.)  * Lambda**2 * (Rplus  + 2*xi*R)
         fMinus3 =  Sqrt(18./35.) * Lambda**2 * (Rminus + 2*xi*R)
         fPlus3  =  Sqrt(18./35.) * Lambda**2 * (Rplus  + 2*xi*R)
         f0Minus =  2./Sqrt(35.)  * Lambda**2 * c
         f0Plus  =  -1. * f0Minus


      case (F37_1950)

         fMinus1 =  -Sqrt(6./35.) * Lambda**2 * (Rminus + 2*xi*R)
         fPlus1  =   Sqrt(6./35.) * Lambda**2 * (Rplus  + 2*xi*R)
         fMinus3 =  -Sqrt(2./7.)  * Lambda**2 * (Rminus + 2*xi*R)
         fPlus3  =   Sqrt(2./7.)  * Lambda**2 * (Rplus  + 2*xi*R)
         f0Minus = 2*Sqrt(6./35.) * Lambda**2 * c
         f0Plus  =  f0Minus


      case (P11_1710)

         at       = Sqrt(3./8.) * (1.-2*xi) *  Lambda**2 * S
         bt       = 1./Sqrt(6.) * (Lambda**2 * c - 2 * Lambda * b)

         fMinus1 =  1./Sqrt(6.) *  Lambda**2 * (Rminus + 3*xi*R)
         fPlus1  =  1./Sqrt(6.) *  Lambda**2 * (Rplus  + 3*xi*R)
         fMinus3 =  0.
         fPlus3  =  0.
         f0Minus =  at-bt
         f0Plus  =  at+bt


      case (F17_1990)

         fMinus1 = 0.
         fPlus1  = 0.
         fMinus3 = 0.
         fPlus3  = 0.
         f0Minus = 0.
         f0Plus  = 0.

      case default

         !write(*,*) 'not implemented -> STOP'
         nu_MaEl_NCprot=0.
         return
      end select


      sigm=U**2*Qs/qq**2*(fMinus3**2+fMinus1**2)
      sigp=V**2*Qs/qq**2*(fPlus3**2+fPlus1**2)
      sign=2.*U*V*Mi**2/Mf**2*(f0Plus**2+f0Minus**2)

      nu_MaEl_NCprot=32.*Mf**2*enu**2*(sigm+sigp+sign)


    end function nu_MaEl_NCprot


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    real function nu_MaEl_NCneut()

      select case (finalstate_ID)

      case (delta)

         fMinus1 =  -sqrt(2.) * (rminus + 2 * xi * R)
         fPlus1  =   sqrt(2.) * (rplus  + 2 * xi * R)
         fMinus3 =  -Sqrt(6.) * (rminus + 2 * xi * R)
         fPlus3  =   Sqrt(6.) * (rplus  + 2 * xi * R)
         f0Minus =  2.*sqrt(2.) * c
         f0Plus  =   f0Minus


      case (S11_1535)

         at       = sqrt(3./2.) * (1.-2*xi) * Lambda * S
         bt       = Sqrt(2./3.) * (Lambda * c - 3*b)

         fMinus1 = -1.*Sqrt(3.) * (Tminus + 2*xi*T) - Sqrt(2./3.) * Lambda * (Rminus + xi*R)
         fPlus1  =     Sqrt(3.) * (Tplus  + 2*xi*T) + Sqrt(2./3.) * Lambda * (Rplus  + xi*R)
         fMinus3 =  0.
         fPlus3  =  0.
         f0Minus =  at-bt
         f0Plus  = -at-bt


      case (D13_1520)

         at       = Sqrt(3.) * (1.-2*xi) * Lambda * S
         bt       = sqrt(4./3.) * Lambda * c

         fMinus1 = -sqrt(3./2.) * (Tminus + 2* xi * T) + sqrt(4./3.) * Lambda  * (Rminus + xi * R)
         fPlus1  = -sqrt(3./2.) * (Tplus  + 2* xi * T) + sqrt(4./3.) * Lambda  * (Rplus  + xi * R)
         fMinus3 = -3./Sqrt(2.) * (Tminus + 2* xi * T)
         fPlus3  = -3./Sqrt(2.) * (Tplus  + 2* xi * T)
         f0Minus =  at - bt
         f0Plus  =  at + bt


      case (S11_1650)

         fMinus1 = -1./Sqrt(24.) * Lambda * (Rminus + 4* xi * R)
         fPlus1  =  1./Sqrt(24.) * Lambda * (Rplus  + 4* xi * R)
         fMinus3 =  0.
         fPlus3  =  0.
         f0Minus =  1./Sqrt(6.) * (Lambda * c - 3*b)
         f0Plus  =  f0Minus


      case (D13_1700)

         fMinus1 = -1./Sqrt(120.) * Lambda * (Rminus + 4* xi * R)
         fPlus1  = -1./Sqrt(120.) * Lambda * (Rplus  + 4* xi * R)
         fMinus3 = -3./Sqrt(40.)  * Lambda * (Rminus + 4* xi * R)
         fPlus3  = -3./Sqrt(40.)  * Lambda * (Rplus  + 4* xi * R)
         f0Minus = -1./Sqrt(30.)  * Lambda * c
         f0Plus  =  -1.* f0Minus


      case (D15_1675)

         fMinus1 =  Sqrt(3./40.) * Lambda * (Rminus + 4* xi * R)
         fPlus1  = -Sqrt(3./40.) * Lambda * (Rplus  + 4* xi * R)
         fMinus3 =  Sqrt(3./20.) * Lambda * (Rminus + 4* xi * R)
         fPlus3  = -Sqrt(3./20.) * Lambda * (Rplus  + 4* xi * R)
         f0Minus = -Sqrt(3./10.) * Lambda * c
         f0Plus  =  f0Minus


      case (S31_1620)

         at       = sqrt(3./2.) * (1.-2.*xi) * Lambda * S
         bt       = 1./Sqrt(6.) * (Lambda * c - 3*b)

         fMinus1 =  Sqrt(3.) * (Tminus + 2*xi*T) - 1./Sqrt(6.) * Lambda * (Rminus + 2*xi*R)
         fPlus1  = -Sqrt(3.) * (Tplus  + 2*xi*T) + 1./Sqrt(6.) * Lambda * (Rplus  + 2*xi*R)
         fMinus3 =  0.
         fPlus3  =  0.
         f0Minus = -at-bt
         f0Plus  =  at-bt


      case (D33_1700)

         at       = Sqrt(3.) * (1.-2.*xi) * Lambda * S
         bt       = 1./Sqrt(3.) * Lambda * c

         fMinus1 = sqrt(3./2.) * (Tminus + 2*xi*T) + 1./Sqrt(3.) * Lambda * (Rminus + 2*xi*R)
         fPlus1  = sqrt(3./2.) * (Tplus  + 2*xi*T) + 1./Sqrt(3.) * Lambda * (Rplus  + 2*xi*R)
         fMinus3 = 3./Sqrt(2.) * (Tminus + 2*xi*T)
         fPlus3  = 3./Sqrt(2.) * (Tplus  + 2*xi*T)
         f0Minus = -at-bt
         f0Plus  = -at+bt


      case (P11_1440)

         at       = 0.25*Sqrt(3.) * Lambda**2 * S
         bt       = (5./12.)*Sqrt(3.) * (Lambda**2 * c - 2 * Lambda * b)

         fMinus1 = (5./12.)*Sqrt(3.) * Lambda**2 * (Rminus + (8./5.)* xi * R)
         fPlus1  = (5./12.)*Sqrt(3.) * Lambda**2 * (Rplus  + (8./5.)* xi * R)
         fMinus3 = 0.
         fPlus3  = 0.
         f0Minus = at - bt
         f0Plus  = at + bt


      case (P33_1600)

         fMinus1 =  1./Sqrt(6.) * Lambda**2 * (Rminus + 2*xi*R)
         fPlus1  = -1./Sqrt(6.) * Lambda**2 * (Rplus  + 2*xi*R)
         fMinus3 =  1./Sqrt(2.) * Lambda**2 * (Rminus + 2*xi*R)
         fPlus3  = -1./Sqrt(2.) * Lambda**2 * (Rplus  + 2*xi*R)
         f0Minus = -Sqrt(2./3.) * (Lambda**2 * c - 2 * Lambda * b)
         f0Plus  =  f0Minus


      case (P13_1720)

         at       = Sqrt(3./20.) * Lambda**2 * S
         bt       = Sqrt(5./12.) * (Lambda**2 * c - 5 * Lambda * b)

         fMinus1 =  Sqrt(27./40.) * Lambda * Tminus + Sqrt(5./12.) * Lambda**2 * (Rminus + (8./5.)*xi*R)
         fPlus1  = -Sqrt(27./40.) * Lambda * Tplus - Sqrt(5./12.) * Lambda**2 * (Rplus  + (8./5.)*xi*R)
         fMinus3 = -Sqrt(9./40.) * Lambda * Tminus
         fPlus3  =  Sqrt(9./40.) * Lambda * Tplus
         f0Minus = -at+bt
         f0Plus  =  at+bt


      case (F15_1680)

         at       = 3./Sqrt(40.) * Lambda**2 * S
         bt       = Sqrt(5./8.)  * Lambda**2 * c

         fMinus1 =  3./Sqrt(20.) * Lambda * Tminus - Sqrt(5./8.) * Lambda**2 * (Rminus + (8./5.)*xi*R)
         fPlus1  =  3./Sqrt(20.) * Lambda * Tplus - Sqrt(5./8.) * Lambda**2 * (Rplus  + (8./5.)*xi*R)
         fMinus3 =  Sqrt(18./20.) * Lambda * Tminus
         fPlus3  =  Sqrt(18./20.) * Lambda * Tplus
         f0Minus =  -at+bt
         f0Plus  =  -at-bt


      case (P31_1910)

         fMinus1 = -1./Sqrt(15.) * Lambda**2 * (Rminus + 2*xi*R)
         fPlus1  = -1./Sqrt(15.) * Lambda**2 * (Rplus  + 2*xi*R)
         fMinus3 =  0.
         fPlus3  =  0.
         f0Minus = -Sqrt(4./15.) * (Lambda**2 * c - 5 * Lambda * b)
         f0Plus  = -1.* f0Minus


      case (P33_1920)

         fMinus1 =  1./Sqrt(15.) * Lambda**2 * (Rminus + 2*xi*R)
         fPlus1  = -1./Sqrt(15.) * Lambda**2 * (Rplus  + 2*xi*R)
         fMinus3 = -1./Sqrt(5.)  * Lambda**2 * (Rminus + 2*xi*R)
         fPlus3  =  1./Sqrt(5.)  * Lambda**2 * (Rplus  + 2*xi*R)
         f0Minus = -(2./Sqrt(15.)) * (Lambda**2 * c - 5 * Lambda * b)
         f0Plus  =  f0Minus


      case (F35_1905)

         fMinus1 =  1./Sqrt(35.)  * Lambda**2 * (Rminus + 2*xi*R)
         fPlus1  =  1./Sqrt(35.)  * Lambda**2 * (Rplus  + 2*xi*R)
         fMinus3 =  Sqrt(18./35.) * Lambda**2 * (Rminus + 2*xi*R)
         fPlus3  =  Sqrt(18./35.) * Lambda**2 * (Rplus  + 2*xi*R)
         f0Minus =  2./Sqrt(35.)  * Lambda**2 * c
         f0Plus  =  -1. * f0Minus


      case (F37_1950)

         fMinus1 =  -Sqrt(6./35.) * Lambda**2 * (Rminus + 2*xi*R)
         fPlus1  =   Sqrt(6./35.) * Lambda**2 * (Rplus  + 2*xi*R)
         fMinus3 =  -Sqrt(2./7.)  * Lambda**2 * (Rminus + 2*xi*R)
         fPlus3  =   Sqrt(2./7.)  * Lambda**2 * (Rplus  + 2*xi*R)
         f0Minus = 2*Sqrt(6./35.) * Lambda**2 * c
         f0Plus  =  f0Minus


      case (P11_1710)

         at       = Sqrt(3./8.) * (1.-2*xi) * Lambda**2 * S
         bt       = 1./Sqrt(6.) * (Lambda**2 * c - 2 * Lambda * b)

         fMinus1 = -1./Sqrt(6.) * Lambda**2 * (Rminus + xi*R)
         fPlus1  = -1./Sqrt(6.) * Lambda**2 * (Rplus  + xi*R)
         fMinus3 = 0.
         fPlus3  = 0.
         f0Minus = -at+bt
         f0Plus  = -at-bt


      case (F17_1990)

         fMinus1 = 0.
         fPlus1  = 0.
         fMinus3 = 0.
         fPlus3  = 0.
         f0Minus = 0.
         f0Plus  = 0.


      case default

         !write(*,*) 'not implemented -> STOP'
         nu_MaEl_NCneut=0.
         return

      end select


      sigm=U**2*Qs/qq**2*(fMinus3**2+fMinus1**2)
      sigp=V**2*Qs/qq**2*(fPlus3**2+fPlus1**2)
      sign=2.*U*V*Mi**2/Mf**2*(f0Plus**2+f0Minus**2)

      nu_MaEl_NCneut=32.*Mf**2*enu**2*(sigm+sigp+sign)


    end function nu_MaEl_NCneut


    real function nu_MaEl_CC()

      select case (finalstate_ID)

      case (delta)

         fMinus1 = sqrt(2.) * Rminus
         fPlus1  = -sqrt(2.) * Rplus
         fMinus3 = Sqrt(6.) * Rminus
         fPlus3  = -Sqrt(6.) * Rplus
         f0Minus = -2*sqrt(2.) * c
         f0Plus  = f0Minus


      case (S11_1535)

         fMinus1 =  2.*Sqrt(3.) * Tminus + 4./Sqrt(6.) * Lambda * Rminus
         fPlus1  = -2.*Sqrt(3.) * Tplus  - 4./Sqrt(6.) * Lambda * Rplus
         fMinus3 =  0.
         fPlus3  =  0.
         f0Minus = -Sqrt(6.) * Lambda * S+2 * Sqrt(2./3.) * (Lambda * c - 3.* b)
         f0Plus  =  Sqrt(6.) * Lambda * S+2 * Sqrt(2./3.) * (Lambda * c - 3.* b)


      case (D13_1520)

         at = 2.* Sqrt(3.) * Lambda * S
         bt = (4./Sqrt(3.))* Lambda * c

         fMinus1 =  Sqrt(6.) * Tminus - 4./Sqrt(3.) * Lambda * Rminus
         fPlus1  =  Sqrt(6.) * Tplus  - 4./Sqrt(3.) * Lambda * Rplus
         fMinus3 =  6./sqrt(2.) * Tminus
         fPlus3  =  6./sqrt(2.) * Tplus
         f0Minus =  -at+bt
         f0Plus  =  -at-bt


      case (S11_1650)

         fMinus1 =  1./Sqrt(6.) * Lambda * Rminus
         fPlus1  = -1./Sqrt(6.) * Lambda * Rplus
         fMinus3 =  0.
         fPlus3  =  0.
         f0Minus = -Sqrt(2./3.) * (Lambda * c - 3.* b)
         f0Plus  =  f0Minus


      case (D13_1700)

         fMinus1 =  1./Sqrt(30.) * Lambda * Rminus
         fPlus1  =  1./Sqrt(30.) * Lambda * Rplus
         fMinus3 =  3./Sqrt(10.) * Lambda * Rminus
         fPlus3  =  3./Sqrt(10.) * Lambda * Rplus
         f0Minus =  Sqrt(2./15.) * Lambda * c
         f0Plus  =  -1. * f0Minus


      case (D15_1675)

         fMinus1 = -Sqrt(3./10.) * Lambda * Rminus
         fPlus1  =  Sqrt(3./10.) * Lambda * Rplus
         fMinus3 = -Sqrt(3./5.)  * Lambda * Rminus
         fPlus3  =  Sqrt(3./5.)  * Lambda * Rplus
         f0Minus =  Sqrt(6./5.)  * Lambda * c
         f0Plus  =  f0Minus


      case (S31_1620)

         at = sqrt(3./2.) * Lambda * S
         bt = 1./Sqrt(6.) * (Lambda * c - 3.* b)

         fMinus1 = -Sqrt(3.) * Tminus + 1./Sqrt(6.) * Lambda * Rminus
         fPlus1  =  Sqrt(3.) * Tplus  - 1./Sqrt(6.) * Lambda * Rplus
         fMinus3 =  0.
         fPlus3  =  0.
         f0Minus =  at+bt
         f0Plus  = -at+bt


      case (D33_1700)

         at = Sqrt(3.)   * Lambda * S
         bt = 1./Sqrt(3.) * Lambda * c

         fMinus1 = -sqrt(3./2.) * Tminus - 1./Sqrt(3.) * Lambda * Rminus
         fPlus1  = -sqrt(3./2.) * Tplus  - 1./Sqrt(3.) * Lambda * Rplus
         fMinus3 = -3./Sqrt(2.) * Tminus
         fPlus3  = -3./Sqrt(2.) * Tplus
         f0Minus =  at + bt
         f0Plus  =  at - bt


      case (P11_1440)

         at  = Sqrt(3./4.) * Lambda**2 * S
         bt  = 5.*Sqrt(3.)/6. * (Lambda**2 * c - 2 * Lambda * b)

         fMinus1 =  -5.*Sqrt(3.)/6. * Lambda**2 * Rminus
         fPlus1  =  -5.*Sqrt(3.)/6. * Lambda**2 * Rplus
         fMinus3 =   0.
         fPlus3  =   0.
         f0Minus =  -at+bt
         f0Plus  =  -at-bt


      case (P33_1600)

         fMinus1 = -1./Sqrt(6.) * Lambda**2 * Rminus
         fPlus1  =  1./Sqrt(6.) * Lambda**2 * Rplus
         fMinus3 = -1./Sqrt(2.) * Lambda**2 * Rminus
         fPlus3  =  1./Sqrt(2.) * Lambda**2 * Rplus
         f0Minus =  Sqrt(2./3.) * (Lambda**2 * c - 2 * Lambda * b)
         f0Plus  =  f0Minus


      case (P13_1720)

         at       = Sqrt(3./5.) * Lambda**2 * S
         bt       = Sqrt(5./3.) * (Lambda**2 * c - 5 * Lambda * b)

         fMinus1 =  -Sqrt(27./10.) * Lambda * Tminus - Sqrt(5./3.) * Lambda**2 * Rminus
         fPlus1  =   Sqrt(27./10.) * Lambda * Tplus + Sqrt(5./3.) * Lambda**2 * Rplus
         fMinus3 =   3./Sqrt(10.) * Lambda * Tminus
         fPlus3  =  -3./Sqrt(10.) * Lambda * Tplus
         f0Minus =   at-bt
         f0Plus  =  -at-bt


      case (F15_1680)

         at   = Sqrt(9./10.) * Lambda**2 * S
         bt   = Sqrt(5./2.)  * Lambda**2 * c

         fMinus1 = -3./Sqrt(5.)  * Lambda * Tminus + Sqrt(5./2.) * Lambda**2 * Rminus
         fPlus1  = -3./Sqrt(5.)  * Lambda * Tplus + Sqrt(5./2.) * Lambda**2 * Rplus
         fMinus3 = -Sqrt(18./5.) * Lambda * Tminus
         fPlus3  = -Sqrt(18./5.) * Lambda * Tplus
         f0Minus =  at - bt
         f0Plus  =  at + bt


      case (P31_1910)

         fMinus1 =  1./Sqrt(15.) * Lambda**2 * Rminus
         fPlus1  =  1./Sqrt(15.) * Lambda**2 * Rplus
         fMinus3 =  0.
         fPlus3  =  0.
         f0Minus =  2./Sqrt(15.) * (Lambda**2 * c - 5 * Lambda * b)
         f0Plus  = -1.* f0Minus


      case (P33_1920)

         fMinus1 = -1./Sqrt(15.) * Lambda**2 * Rminus
         fPlus1  =  1./Sqrt(15.) * Lambda**2 * Rplus
         fMinus3 =  1./Sqrt(5.)  * Lambda**2 * Rminus
         fPlus3  = -1./Sqrt(5.)  * Lambda**2 * Rplus
         f0Minus =  2./Sqrt(15.) * (Lambda**2 * c - 5 * Lambda * b)
         f0Plus  =  f0Minus


      case (F35_1905)

         fMinus1 =  -1./Sqrt(35.)  * Lambda**2 * Rminus
         fPlus1  =  -1./Sqrt(35.)  * Lambda**2 * Rplus
         fMinus3 =  -Sqrt(18./35.) * Lambda**2 * Rminus
         fPlus3  =  -Sqrt(18./35.) * Lambda**2 * Rplus
         f0Minus =  -2./Sqrt(35.)  * Lambda**2 * c
         f0Plus  =  -1.* f0Minus


      case (F37_1950)

         fMinus1 =  Sqrt(6./35.)  * Lambda**2 * Rminus
         fPlus1  = -Sqrt(6./35.)  * Lambda**2 * Rplus
         fMinus3 =  Sqrt(2./7.)   * Lambda**2 * Rminus
         fPlus3  = -Sqrt(2./7.)   * Lambda**2 * Rplus
         f0Minus = -Sqrt(24./35.) * Lambda**2 * c
         f0Plus  =  f0Minus


      case (P11_1710)

         at  = sqrt(3./2.) * Lambda**2 * S
         bt  = Sqrt(2./3.) * (Lambda**2 * c - 2 * Lambda * b)

         fMinus1 = Sqrt(2./3.) * Lambda**2 * Rminus
         fPlus1  = Sqrt(2./3.) * Lambda**2 * Rplus
         fMinus3 = 0.
         fPlus3  = 0.
         f0Minus = at - bt
         f0Plus  = at + bt


      case (F17_1990)

         fMinus1 =  -Sqrt(3./35.) * Lambda**2 * Rminus
         fPlus1  =   Sqrt(3./35.) * Lambda**2 * Rplus
         fMinus3 =  -1./Sqrt(7.)  * Lambda**2 * Rminus
         fPlus3  =   1./Sqrt(7.)  * Lambda**2 * Rplus
         f0Minus =   Sqrt(6./35.) * Lambda**2 * c
         f0Plus  =   f0Minus


      case default


         !write(*,*) 'not implemented -> STOP'
         nu_MaEl_CC=0.
         return

      end select

      sigm=U**2*Qs/qq**2*(fMinus3**2+fMinus1**2)
      sigp=V**2*Qs/qq**2*(fPlus3**2+fPlus1**2)
      sign=2.*U*V*Mi**2/Mf**2*(f0Plus**2+f0Minus**2)

      nu_MaEl_CC=32.*Mf**2*enu**2*(sigm+sigp+sign)


    end function nu_MaEl_CC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    real function nu_MaEl_EMprot()

      select case (finalstate_ID)

      case (delta)

         fPlus1  =  Sqrt(2.) * R
         fPlus3  =  Sqrt(6.) * R
         fMinus1 = -1 * fPlus1
         fMinus3 = -1 * fPlus3
         f0Minus =  0.
         f0Plus  =  0.


      case (S11_1535)

         fMinus1 =  Sqrt(3.) * T + sqrt(3./2.) * Lambda * R
         f0Minus = -sqrt(3./2.) * Lambda * S
         fPlus1  = -1. * fMinus1
         f0Plus  = -1. * f0Minus
         fMinus3 =  0.
         fPlus3  =  0.


      case (D13_1520)

         fMinus1 =  sqrt(3./2.) * T - Sqrt(3.) * Lambda * R
         fMinus3 =  3./Sqrt(2.) * T
         f0Minus = -Sqrt(3.) * Lambda * S
         fPlus1  =  fMinus1
         fPlus3  =  fMinus3
         f0Plus  =  f0Minus


      case (S11_1650)

         fMinus1 = 0.
         fPlus1  = 0.
         fMinus3 = 0.
         fPlus3  = 0.
         f0Minus = 0.
         f0Plus  = 0.


      case (D13_1700)

         fMinus1 = 0.
         fPlus1  = 0.
         fMinus3 = 0.
         fPlus3  = 0.
         f0Minus = 0.
         f0Plus  = 0.


      case (D15_1675)

         fMinus1 = 0.
         fPlus1  = 0.
         fMinus3 = 0.
         fPlus3  = 0.
         f0Minus = 0.
         f0Plus  = 0.


      case (S31_1620)

         fMinus1 =  Sqrt(3.) * T - 1./Sqrt(6.) * Lambda * R
         f0Minus = -sqrt(3./2.) * Lambda * S
         fPlus1  = -1. * fMinus1
         f0Plus  = -1. * f0Minus
         fMinus3 = 0.
         fPlus3  = 0.


      case (D33_1700)

         fMinus1 =  sqrt(3./2.) * T + 1./Sqrt(3.) * Lambda * R
         fMinus3 =  3./Sqrt(2.) * T
         f0Minus = -Sqrt(3.) * Lambda * S
         fPlus1  =  fMinus1
         fPlus3  =  fMinus3
         f0Plus  =  f0Minus


      case (P11_1440)

         fMinus1 = -0.5*Sqrt(3.) * Lambda**2 * R
         fPlus1  =  fMinus1
         fMinus3 =  0.
         fPlus3  =  0.
         f0Minus = -0.5*Sqrt(3.) * Lambda**2 * S
         f0Plus  =  f0Minus


      case (P33_1600)

         fMinus1 = 1./Sqrt(6.) * Lambda**2 * R
         fMinus3 = 1./Sqrt(2.) * Lambda**2 * R
         fPlus1  = -1. * fMinus1
         fPlus3  = -1. * fMinus3
         f0Minus = 0.
         f0Plus  = 0.


      case (P13_1720)

         fMinus1 = -Sqrt(27./10.) * Lambda * T - Sqrt(3./5.) * Lambda**2 * R
         fMinus3 =  3./Sqrt(10.) * Lambda * T
         f0Minus =  Sqrt(3./5.)  * Lambda**2 * S
         fPlus1  = -1. * fMinus1
         fPlus3  = -1. * fMinus3
         f0Plus  = -1. * f0Minus


      case (F15_1680)

         fMinus1 =  -3./Sqrt(5.)  * Lambda * T + 3./Sqrt(10.) * Lambda**2 * R
         fMinus3 =  -Sqrt(18./5.) * Lambda * T
         f0Minus =   3./Sqrt(10.) * Lambda**2 * S
         fPlus1  =  fMinus1
         fPlus3  =  fMinus3
         f0Plus  =  f0Minus


      case (P31_1910)

         fMinus1 = -1./Sqrt(15.) * Lambda**2 * R
         fPlus1  = fMinus1
         fMinus3 = 0.
         fPlus3  = 0.
         f0Minus = 0.
         f0Plus  = 0.


      case (P33_1920)

         fMinus1 =  1./Sqrt(15.) * Lambda**2 * R
         fMinus3 = -1./Sqrt(5.)  * Lambda**2 * R
         fPlus1  = -1.* fMinus1
         fPlus3  = -1.* fMinus3
         f0Minus =  0.
         f0Plus  =  0.


      case (F35_1905)

         fMinus1 = 1./Sqrt(35.)  * Lambda**2 * R
         fMinus3 = Sqrt(18./35.) * Lambda**2 * R
         fPlus1  = fMinus1
         fPlus3  = fMinus3
         f0Minus = 0.
         f0Plus  = 0.


      case (F37_1950)

         fMinus1 = -Sqrt(6./35.) * Lambda**2 * R
         fMinus3 = -Sqrt(2./7.)  * Lambda**2 * R
         fPlus1  = -1. * fMinus1
         fPlus3  = -1. * fMinus3
         f0Minus = 0.
         f0Plus  = 0.


      case (P11_1710)

         fMinus1 = Sqrt(3./8.) * Lambda**2 * R
         f0Minus = Sqrt(3./8.) * Lambda**2 * S
         fPlus1  = fMinus1
         f0Plus  = f0Minus
         fMinus3 = 0.
         fPlus3  = 0.


      case (F17_1990)

         fMinus1 = 0.
         fPlus1  = 0.
         fMinus3 = 0.
         fPlus3  = 0.
         f0Minus = 0.
         f0Plus  = 0.


      case default

         !write(*,*) 'not implemented -> STOP'
         nu_MaEl_EMprot=0.
         return

      end select

!      sigm=U**2*Qs/qq**2*(fMinus3**2+fMinus1**2)
!      sigp=V**2*Qs/qq**2*(fPlus3**2+fPlus1**2)
!      sign=2.*U*V*Mi**2/Mf**2*(f0Plus**2+f0Minus**2)
!      nu_MaEl_EMprot=32.*Mf**2*enu**2*(sigm+sigp+sign)

      sigm=0.5*(fMinus3**2+fMinus1**2)
      sigp=0.5*(fPlus3**2+fPlus1**2)
      sign=epsilon*Qs/qq**2*Mf**2/Mi**2*(f0Plus**2+f0Minus**2)
      nu_MaEl_EMprot=sigm+sigp+sign



    end function nu_MaEl_EMprot


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    real function nu_MaEl_EMneut()

      select case (finalstate_ID)

      case (delta)

         fPlus1  =  sqrt(2.) * R
         fPlus3  =  sqrt(6.) * R
         fMinus1 = -1. * fPlus1
         fMinus3 = -1. * fPlus3
         f0Minus =  0.
         f0Plus  =  0.


      case (S11_1535)

         fPlus1  =  sqrt(3.) * T + 1./sqrt(6.) * Lambda * R
         f0Minus =  sqrt(3./2.) * Lambda * S
         fMinus1 = -1. * fPlus1
         f0Plus  = -1. * f0Minus
         fMinus3 =  0.
         fPlus3  =  0.



      case (D13_1520)

         fMinus1 = -sqrt(3./2.) * T + sqrt(1./3.) * Lambda * R
         fMinus3 = -sqrt(9./2.) * T
         f0Minus =  sqrt(3.) * Lambda * S
         fPlus1  =  fMinus1
         fPlus3  =  fMinus3
         f0Plus  =  f0Minus


      case (S11_1650)

         fPlus1  =  sqrt(1./6.) * Lambda * R
         fMinus1 = -1 * fPlus1
         fMinus3 =  0.
         fPlus3  =  0.
         f0Minus =  0.
         f0Plus  =  0.


      case (D13_1700)

         fMinus1 = -(1./Sqrt(30.)) * Lambda * R
         fMinus3 = -(3./Sqrt(10.)) * Lambda * R
         fPlus1  =  fMinus1
         fPlus3  =  fMinus3
         f0Minus =  0.
         f0Plus  =  0.


      case (D15_1675)

         fMinus1 = sqrt(3./10.) * Lambda * R
         fMinus3 = sqrt(3./5.)  * Lambda * R
         fPlus1  = -1 * fMinus1
         fPlus3  = -1 * fMinus3
         f0Minus =  0.
         f0Plus  =  0.


      case (S31_1620)

         fMinus1 =  sqrt(3.) * T - sqrt(1./6.) * Lambda * R
         f0Minus = -sqrt(3./2.) * Lambda * S
         fPlus1  = -1. * fMinus1
         f0Plus  = -1. * f0Minus
         fMinus3 = 0.
         fPlus3  = 0.


      case (D33_1700)

         fMinus1 =  sqrt(3./2.) * T + sqrt(1./3.) * Lambda * R
         fMinus3 =  sqrt(9./2.) * T
         f0Minus = -sqrt(3.) * Lambda * S
         fPlus1  =  fMinus1
         fPlus3  =  fMinus3
         f0Plus  =  f0Minus


      case (P11_1440)

         fMinus1 = sqrt(1./3.) * Lambda**2 * R
         fPlus1  = fMinus1
         fMinus3 =  0.
         fPlus3  =  0.
         f0Minus =  0.
         f0Plus  =  0.


      case (P33_1600)

         fMinus1 = sqrt(1./6.) * Lambda**2 * R
         fMinus3 = sqrt(1./2.) * Lambda**2 * R
         fPlus1  = -1. * fMinus1
         fPlus3  = -1. * fMinus3
         f0Minus = 0.
         f0Plus  = 0.


      case (P13_1720)

         fMinus1 = sqrt(4./15.) * Lambda**2 * R
         fPlus1  = -1 * fMinus1
         fMinus3 =  0.
         fPlus3  =  0.
         f0Minus =  0.
         f0Plus  =  0.


      case (F15_1680)

         fMinus1 =  -sqrt(2./5.) * Lambda**2 * R
         fPlus1  =  fMinus1
         fMinus3 =  0.
         fPlus3  =  0.
         f0Minus =  0.
         f0Plus  =  0.


      case (P31_1910)

         fMinus1 =  -sqrt(1./15.) * Lambda**2 * R
         fPlus1  =  fMinus1
         fMinus3 =  0.
         fPlus3  =  0.
         f0Minus =  0.
         f0Plus  =  0.


      case (P33_1920)

         fMinus1 =  sqrt(1./15.) * Lambda**2 * R
         fMinus3 = -sqrt(1./5.)  * Lambda**2 * R
         fPlus1  = -1.* fMinus1
         fPlus3  = -1.* fMinus3
         f0Minus =  0.
         f0Plus  =  0.


      case (F35_1905)

         fMinus1 = sqrt(1./35.)  * Lambda**2 * R
         fMinus3 = sqrt(18./35.) * Lambda**2 * R
         fPlus1  = fMinus1
         fPlus3  = fMinus3
         f0Minus = 0.
         f0Plus  = 0.


      case (F37_1950)

         fMinus1 = -Sqrt(6./35.) * Lambda**2 * R
         fMinus3 = -Sqrt(2./7.) * Lambda**2 * R
         fPlus1  = -1. * fMinus1
         fPlus3  = -1. * fMinus3
         f0Minus = 0.
         f0Plus  = 0.


      case (P11_1710)

         fMinus1 = -Sqrt(1./24.) * Lambda**2 * R
         f0Minus = -Sqrt(3./8.)  * Lambda**2 * S
         fPlus1  = fMinus1
         f0Plus  = f0Minus
         fMinus3 = 0.
         fPlus3  = 0.



      case (F17_1990)

         fMinus1 = Sqrt(3./35.) *  Lambda**2 * R
         fPlus1  = -1 * fMinus1
         fMinus3 = 1./Sqrt(7.)  *  Lambda**2 * R
         fPlus3  = -1 * fMinus3
         f0Minus =  0.
         f0Plus  =  0.


      case default

         !write(*,*) 'not implemented -> STOP'
         nu_MaEl_EMneut=0.
         return

      end select


!      sigm=U**2*Qs/qq**2*(fMinus3**2+fMinus1**2)
!      sigp=V**2*Qs/qq**2*(fPlus3**2+fPlus1**2)
!      sign=2.*U*V*Mi**2/Mf**2*(f0Plus**2+f0Minus**2)
!      nu_MaEl_EMneut=32.*Mf**2*enu**2*(sigm+sigp+sign)

      sigm=0.5*(fMinus3**2+fMinus1**2)
      sigp=0.5*(fPlus3**2+fPlus1**2)
      sign=epsilon*Qs/qq**2*Mf**2/Mi**2*(f0Plus**2+f0Minus**2)
      nu_MaEl_EMneut=sigm+sigp+sign



    end function nu_MaEl_EMneut






  end function nu_MaEl_RS

end module nu_MatrixElementRS
