!******************************************************************************
!****m* /matrixElementQE
! NAME
! module matrixElementQE
!
! PURPOSE
! calculates the spin summed and averaged matrix element for EM, CC and NC
! QE scattering: l + N -> l' + N'
! Details are given below.
!******************************************************************************
module matrixElementQE
  implicit none
  private

  public:: matrixElementforQE

  !****************************************************************************
  !****g* neutrinoXsection/useQEextraterm
  ! SOURCE
  logical, save :: useQEextraterm=.true.
  ! PURPOSE
  ! switch on/off an extra term appearing in the current
  ! due to different masses of in- and outgoing nucleons
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoXsection/useCorrelations
  ! SOURCE
  logical, save :: useCorrelations=.false.
  ! PURPOSE
  ! switch on/off RPA correlations according to
  ! Nieves, Amaro, Valverde, PRC70, 055503 (2004)
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoXsection/nievesCorr_para
  ! SOURCE
  integer, save :: nievesCorr_para=2
  ! PURPOSE
  ! if RPA correlations are switched on, this parameter decides which set
  ! of varibles to use:
  ! * 1: modified Nieves et al., PRC70, 055503 (2004)
  ! * 2: original Nieves et al., PRC70, 055503 (2004)
  ! * 3: Tselyaev, Speth et al., PRC75, 014315 (2007)
  !****************************************************************************


  !****************************************************************************
  !****g* neutrinoXsection/gp
  ! SOURCE
  real, save :: gp=0.63
  ! PURPOSE
  ! vary gp if RPA correlations are switched on
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoXsection/withScalarInt
  ! SOURCE
  logical, save :: withScalarInt=.true.
  ! PURPOSE
  ! switch on/off scalar interactions
  !****************************************************************************

  logical, save :: initFlag=.true.


contains


  !****************************************************************************
  !****s* matrixElementQE/readInput_MatrixElementQE
  ! NAME
  ! subroutine readInput_MatrixElementQE
  ! PURPOSE
  ! This subroutine reads namelist 'MatrixElementQE'.
  !****************************************************************************
  subroutine readInput_MatrixElementQE
    use output, only: Write_ReadingInput

    integer :: ios

    !**************************************************************************
    !****n* matrixElementQE/MatrixElementQE
    ! NAME
    ! NAMELIST /MatrixElementQE/
    ! PURPOSE
    ! This Namelist includes:
    ! * useQEextraterm
    ! * useCorrelations
    ! * nievesCorr_para
    ! * gp
    ! * withScalarInt
    !**************************************************************************
    NAMELIST /matrixElementQE/  useQEextraterm,useCorrelations,nievesCorr_para,gp,withScalarInt

    call Write_ReadingInput('matrixElementQE',0)
    rewind(5)
    read(5,nml=matrixElementQE,IOSTAT=ios)
    call Write_ReadingInput("matrixElementQE",0,ios)
    call Write_ReadingInput('matrixElementQE',1)

    write(*,*) 'use extraterm in QE',useQEextraterm

    if (useCorrelations) then
       write(*,*) 'Nieves correlations are used -> extraterm is automatically included'
       if (nievesCorr_para.eq.1) then
          write(*,*) 'L. Alvarez-Ruso parameters'
       else if (nievesCorr_para.eq.2) then
          write(*,*) 'Nieves original parameters'
       else if (nievesCorr_para.eq.3) then
          write(*,*) 'updated Speth parameters'
       else
          write(*,*) 'not valid input for nievesCorr_para -> stop', nievesCorr_para
       end if
       write(*,*) 'gp=',gp
       write(*,*) 'withScalarInt', withScalarInt
    else
       write(*,*) 'no correlations'
    end if

  end subroutine readInput_MatrixElementQE


  !****************************************************************************
  !****f* matrixElementQE/MatrixElementforQE
  ! NAME
  ! real function MatrixElementforQE(process_ID,charge_in,k_in,k_out,p_in,p_out,position)
  !
  ! PURPOSE
  ! This function returns the spin summed and averaged matrix element for
  ! CC, NC, EM induced quasielastic scattering.
  !
  ! The notation is as follows:
  ! * MatrixElementforQE=coupling * L_munu H^munu
  !   with coupling as in line 173
  ! * L_munu=Tr((kslash + ml)Atilde(kprslash+mlpr)A) [*1/2 for EM spin average]
  !   with A_mu=gamma_mu(1-a gamma_5) where a=0 for EM, 1 for NC,CC neutrino, -1 for NC,CC antineutrino
  ! * H^munu=1/2 Tr(...) same as above with hadronic current instead
  !
  ! note: Peskin notation for spin sums (no 1/2M !!)
  !
  ! more details in Tinas notes!
  !
  !
  ! INPUTS
  ! * integer              :: process_ID   -- CC, NC or EM
  ! * integer              :: charge_in    -- charge of incoming nucleon
  ! * real, dimension(0:3) :: p_in         -- momentum of incoming nucleon
  ! * real, dimension(0:3) :: k_in         -- momentum of incoming lepton
  ! * real, dimension(0:3) :: p_out        -- momentum of outgoing nucleon
  ! * real, dimension(0:3) :: k_out        -- momentum of outgoing lepton
  ! * real, dimension(1:3), optional :: position     -- position of incoming nucleon
  !****************************************************************************
  real function MatrixElementforQE(process_ID,charge_in,k_in,k_out,p_in,p_out,position)
    use FF_QE_nucleonScattering
    use leptonicID
    use Minkowski, only:SP, abs4
    use constants, only: pi, alphaQED, mN, GF, coscab

    real, dimension(0:3), intent(in) :: k_in, k_out, p_in, p_out
    integer, intent(in) :: charge_in
    integer, intent(in) :: process_ID
    real, dimension(1:3), intent(in),optional :: position

    real :: F1=0., F2=0., FA=0., FP=0.
    real :: Mi2=0.,Mf2=0.,M2=0.
    real :: extraterm
    real :: s,t,Mi,Mf,ml_in,ml_out,M,Qs
    real :: a
    real :: coupling=0.
    real :: LH,dLH

    if (initFlag) then
       call readInput_MatrixElementQE
       initFlag=.false.
    end if

    MatrixElementforQE=0.

    !reaction possible??
    if (charge_in.eq.proton.and.process_ID.eq.CC) return
    if (charge_in.eq.neutron.and.process_ID.eq.antiCC) return


    if (useQEextraterm) then  !use the extra current conservating term in the vector current
       extraterm=1.
    else
       extraterm=0.
    end if

    select case (process_ID)
    case (NC,CC)
       a=-1.
    case (antiNC,antiCC)
       a=1.
    case (EM,antiEM)
       a=0.
    end select

    t=SP(k_in-k_out,k_in-k_out)
    Qs=-t
    s=SP(k_in+p_in,k_in+p_in)
    ml_in=sqrt(max(SP(k_in,k_in),0.))
    ml_out=sqrt(max(SP(k_out,k_out),0.))
    Mi=abs4(p_in)
    Mf=abs4(p_out)

    if (Qs.lt.10e-8) return

    select case (process_ID)
    case (NC,antiNC)
       coupling=GF**2/2.
    case (CC,antiCC)
       coupling=GF**2/2.*coscab**2
    case (EM,antiEM)
       coupling=(4.*pi*alphaQED)**2/Qs**2
    end select

    if (.not.present(position).and.useCorrelations) then
       write(*,*) 'QE with correlations, position is not present',present(position)
       stop
    end if

    if (useCorrelations.and.process_ID.lt.0.) then
       write(*,*) 'correlations for anti-leptons not yet implemented -> STOP'
       stop
    end if

    call formfactors_QE(Qs,process_ID,charge_in,F1,F2,FA,FP)

    M=mN
    M2=M**2
    Mi2=Mi**2
    Mf2=Mf**2

    if (abs(process_ID).eq.EM) coupling=coupling/2./2.
    !factor 1./2. from leptonic spin averaging
    !additional 1/2. required - see Mathematica notebook
    !leptonic tensor was calculated including axial part and thus had a
    !factor of 2 in there, has to be divided out for electrons


    LH = (1/(M2*t**2))*(2*(2*extraterm*F1*M*(Mf - Mi)*(Mf2 - Mi2 + t)*                       &
         (2*F1*M*(Mf + Mi) + F2*t) - t*(8*F1**2*Mf*(Mi - Mf)*M2 + 2*F1*F2*((3*Mf + Mi)*t -   &
         (Mf - Mi)**2*(Mf + Mi))*M + t*(((3*Mf - Mi)*(Mf + Mi) + t)*F2**2 + 16*Fa*Fp*M*Mf +  &
         4*Fp**2*(t - (Mf - Mi)**2))))*(ml_in)**4 + (4*(2*extraterm*F1*M*(Mf + Mi)*(2*F1*M*  &
         (Mf + Mi) + F2*t)*(Mf - Mi)**2 + t*(4*F1**2*((Mf - Mi)**2 - 2*t)*M2 + 2*F1*F2*      &
         (Mf + Mi)*((Mf - Mi)**2 - 2*t)*M - t*(((Mf + Mi)**2 - t)*F2**2 + 8*Fa**2*M2 +       &
         8*Fa*Fp*M*(Mf + Mi) + 4*Fp**2*(t - (Mf - Mi)**2))))* ml_out**2 + t*(8*F1**2*        &
         (-2*t**2 - ((Mf - Mi)*(2*Mi + extraterm*(Mf + Mi)) + 4*s)*t + (Mf - Mi)*(Mf + Mi)*  &
         ((extraterm + 2)*Mf2 + extraterm*Mi2 - 2*(extraterm + 1)*s))*M2 + 4*F1*t*((Mf - Mi) &
         *(8*a*Fa*M*(Mf + Mi) + F2*((extraterm + 1)*Mf2 + (extraterm - 3)*Mi2 - 2*(extraterm &
         - 1)*s)) + (8*a*Fa*M - (extraterm + 1)*F2*Mf + (extraterm - 3)*F2* Mi)*t)*M +       &
         2*t*((t**2 + ((Mf - 3*Mi)*(Mf + Mi) + 4*s)* t - 2*(Mf - Mi)*(Mf + Mi)*(Mf2 +        &
         Mi2 - 2*s))* F2**2 + 8*Fa**2*M2*((Mf + Mi)**2 - 2*s - t) + 4*Fp**2*t*(t - (Mf - Mi) &
         **2) + 8*Fa*M* (-2*Fp*(Mf - Mi)*(Mf2 - s) + 2*Fp*Mf*t + a*F2*(Mf + Mi)*(Mf2 - Mi2   &
         + t)))))*(ml_in)**2 +  2*(2*extraterm*F1*M*(Mf - Mi)*(Mf2 - Mi2 - t)* (2*F1*M*      &
         (Mf + Mi) + F2*t) - t*(8*F1**2*(Mf - Mi)*Mi*M2 + 2*F1*F2*((Mf + 3*Mi)*t - (Mf - Mi) &
         **2*(Mf + Mi))*M + t*((t - (Mf - 3*Mi)*(Mf + Mi))*F2**2 + 4*Fp*(4*Fa*M*Mi + Fp*     &
         (t - (Mf - Mi)**2)))))*(ml_out)**4 + 2*t*((F2**2 + 4*Fp**2)*t**3 + ((-((3*Mf - Mi)* &
         (Mf + Mi) - 4*s))*F2**2 - 8*F1**2*M2 - 8*Fa**2*M2 - 4*Fp**2*(Mf - Mi)**2 +          &
         8*Fa*M*(2*Fp*Mi + a*F2*(Mf + Mi)) + F1*M*(16*a*Fa*M - 2*F2*(3*Mf + Mi)))*t**2 +     &
         2*((Mf - Mi)*(Mf + Mi)*(Mf2 + Mi2 - 2*s)*F2**2 + 4*F1**2*M**2*(Mf2 - Mi*Mf - 2*s)   &
         + 4*Fa**2*M**2*((Mf + Mi)**2 - 2*s) - 4*Fa*M*(Mf - Mi)* (a*F2*(Mf + Mi)**2 + 2*Fp*  &
         (s - Mi2)) - F1*M*(Mf - Mi)*(8*a*Fa*M*(Mf + Mi) + F2*(-3*Mf2 + Mi2 + 2*s)))*t -     &
         8*F1**2*M**2*(Mf - Mi)*(Mf + Mi)*(Mi2 - s) - 2*extraterm*F1*M*(Mf - Mi)* (Mf2 + Mi2 &
         - 2*s - t)*(2*F1*M*(Mf + Mi) + F2*t))* ml_out**2 + 2*t**2*(2*t*(-Mf**4 + (2*s + t)  &
         *Mf2 + 2*Mi*t*Mf - Mi**4 - 2*s*(s + t) + Mi2*(2*s + t))*F2**2 +  8*a*Fa*M*(Mf + Mi) &
         *(Mf2 + Mi2 - 2*s - t)*t*F2 +  8*F1**2*M2*(t**2 - ((Mf - Mi)**2 - 2*s)*t + 2*(Mf2 - &
         s)*(Mi2 - s)) + 8*Fa**2*M2* (t**2 - ((Mf + Mi)**2 - 2*s)*t + 2*(Mf2 - s)*(Mi2 - s)) &
         +  4*F1*M*t*(4*a*Fa*M*(Mf2 + Mi2 - 2*s - t) + 2*F2*(Mf + Mi)*(t - (Mf - Mi)**2))))


    if (useCorrelations) then
       ! with correlations a la Nieves, Valverde
       call correlationsNieves(process_ID,charge_in,p_in,k_in,k_out,ml_out,s,t,position,dLH)
       LH=LH+dLH
    end if

    MatrixElementforQE = 1./2.*coupling*LH


  end function MatrixElementforQE

  !****************************************************************************
  !****s* matrixElementQE/correlationsNieves
  ! NAME
  ! subroutine correlationsNieves(process_ID,charge_in,p_in,k_in,k_out,ml_out,s_in,t_in,position,dLH)
  !
  ! PURPOSE
  ! This subroutines calculates the RPA correlations for CC, NC and EM processes
  ! which then have to be added to the "bare" matrix element.
  !
  ! Calculation from Nieves, Amaro and Valverde, PRC 70, 055503 (2004)
  ! code basically written by L. Alvarez-Ruso
  !
  !
  ! INPUTS
  ! * integer              :: process_ID   -- CC, NC or EM
  ! * integer              :: charge_in    -- charge of incoming nucleon
  ! * real, dimension(0:3) :: p_in         -- momentum of incoming nucleon
  ! * real, dimension(0:3) :: k_in         -- momentum of incoming lepton
  ! * real, dimension(0:3) :: k_out        -- momentum of outgoing lepton
  ! * real, dimension(1:3) :: position     -- position of incoming nucleon
  ! * real                 :: ml_out       -- mass of outgoing lepton
  ! * real                 :: s_in         -- Mandelstam s
  ! * real                 :: t_in         -- Mandelstam t
  !
  ! OUTPUT
  ! * real                 :: dLH          -- correlation part of the matrix element
  !****************************************************************************
  subroutine correlationsNieves(process_ID,charge_in,p_in,k_in,k_out,ml_out,s_in,t_in,position,dLH)
    !Calculates the correlations for QE processes according to Nieves, Valverde
    !code basically written by Luis Alvarez Ruso
    use constants, only: pi, hbarc, mPi, mN
    use idTable, only: rho
    use ParticleProperties, only: hadron
    use densityModule, only: densityAt
    use dichteDefinition
    use FF_QE_nucleonScattering
    use leptonicID

    integer, intent(in) :: process_ID          ! process: EM, NC, CC
    integer, intent(in) :: charge_in           ! charge of initial nucleon
    real,dimension(0:3),intent(in) :: p_in     ! energy an momentum of the initial nucleon
    real,dimension(0:3),intent(in) :: k_in     ! energy an momentum of the initial lepton
    real,dimension(0:3),intent(in) :: k_out    ! energy an momentum of the final lepton
    real,intent(in) :: s_in,t_in               ! Mandelstam variables
    real,intent(in) :: ml_out                  ! final lepton mass
    real, dimension(1:3), intent(in) :: position
    real,intent(out) :: dLH                    ! correlation part of the matrix element

    real :: F1=0.,F2=0.,FA=0.,FP=0.            ! vector and axial nucleon form factors
    real :: F1_p,F2_p,Fa_p,Fp_p,F1_n,F2_n,Fa_n,Fp_n  ! form factors on protons and neutrons: needed for NC and EM
    real :: F1_s,F2_s,Fa_s,F1_a,F2_a,Fa_a
    complex,dimension(0:3,0:3) :: dH           ! H part with correlations
    real :: pFave                              ! averaged fermi momentum: (pFp+pFn)/2
    real,dimension(0:3) :: k,p,kf              ! various 4 momenta
    real,dimension(0:3) :: q                   ! q=k-kf 4 momentum transfer
    real :: qm                                 ! absolute value of 3-momentum transfer
    real :: cosa,sina
    real,dimension(0:3) :: pr                  ! rotated nucleon momentum
    real,dimension(0:3) :: vp                  ! vector product (k x p)
    real :: t
    real :: rho0ave,rhoave                     ! averaged densities at the origin and at r
    real, dimension(1:3),parameter :: origin=(/0.,0.,0./)
    type(dichte) :: density,density0

    real, parameter :: mn_fm = mN/hbarc
    real, parameter :: mpi_fm= mPi/hbarc

    real,parameter :: hb=hbarc*1000.      ! (MeV x fm)

    real :: f,Lpi,Lrho,crho,c0,fpin,fpex,f0in,f0ex,g0in,g0ex

    real :: mrho_fm

    f=0.95
    Lpi=1200./hb
    !gp=0.63

    if (nievesCorr_para.eq.1) then
       !Luis set of parameters
       Lrho=1400./hb
       crho=3.94
       c0=380./hb
       fpin=0.33
       fpex=0.45
       f0in=0.07
       f0ex=-2.15
       g0in=0.575
       g0ex=0.575
    else if (nievesCorr_para.eq.2) then
       ! parameters in Amaro et al.
       Lrho=2500./hb
       crho=2.
       c0=380./hb
       fpin=0.33
       fpex=0.45
       f0in=0.07
       f0ex=-2.15
       g0in=0.575
       g0ex=0.575
    else if (nievesCorr_para.eq.3) then
       !New values from Tselyaev, Speth, etc nucl-th/06120642
       Lrho=2500./hb
       crho=2.
       c0=300./hb
       fpin=0.76
       fpex=2.30
       f0in=-0.002
       f0ex=-1.35
       g0in=0.05
       g0ex=0.05
    else
       write(*,*) 'no valid input for nievesCorr_para -> stop', nievesCorr_para
       stop
    end if

!!$    gp=0.
!!$    Lpi=0./hb
!!$    Lrho=0./hb
!!$    crho=0.
!!$    c0=380./hb
!!$    fpin=0.
!!$    fpex=0.
!!$    f0in=0.
!!$    f0ex=0.
!!$    g0in=0.
!!$    g0ex=0.


    !convert everything to inverse Fermi
    k=k_in/hbarc
    p=p_in/hbarc
    kf=k_out/hbarc
    t=t_in/hbarc**2

    mrho_fm = hadron(rho)%mass/hbarc

    !need rho0ave,rhoave which are the averaged densities at the origin and at r
    density=densityAt(position)
    density0=densityAt(origin)
    rhoave=(density%proton(0)+density%neutron(0))/2.
    rho0ave=(density0%proton(0)+density0%neutron(0))/2.

    !Rotation to the Lab frame with q paralel to the z axis as in the N-V paper
    q=k-kf
    qm=sqrt(sp(q,q))

    call vector_prod(k,q,vp)
    cosa=sp(k,q)/k(0)/qm
    sina=vp(2)/k(0)/qm

    call rotate(cosa,sina,p,pr)



    select case (process_ID)

    case (CC)

       call formfactors_QE(-t_in,process_ID,charge_in,F1,F2,FA,FP)
       Fp=Fp/mn_fm !Luis uses a different notation for Fp

       call Hadronic_Correl(q(0),qm,t,p(0),pr(3),dH)

       dLH=8*k(0)*(2*k(0) - q(0) - cosa*qm)*dH(0,0) + 8*k(0)*(-4*cosa*k(0) + 2*cosa*q(0) + 2*qm)*dH(0,3) + &
            & 8*k(0)*(2*k(0) - 2*cosa**2*k(0) - q(0) + cosa*qm)*dH(1,1) + 8*k(0)*((0,2)*cosa*q(0) - (0,2)*qm)*dH(1,2) + &
            & 8*k(0)*(2*cosa**2*k(0) - q(0) - cosa*qm)*dH(3,3)


    case (NC)

       call formfactors_QE(-t_in,process_ID,proton,F1_p,F2_p,FA_p,FP_p)
       Fp_p=Fp_p/mn_fm !Luis uses a different notation for Fp

       call formfactors_QE(-t_in,process_ID,neutron,F1_n,F2_n,FA_n,FP_n)
       Fp_n=Fp_n/mn_fm !Luis uses a different notation for Fp

       F1_s=F1_p+F1_n
       F1_a=F1_p-F1_n

       F2_s=F2_p+F2_n
       F2_a=F2_p-F2_n

       Fa_s=Fa_p+Fa_n
       Fa_a=Fa_p-Fa_n

       call Hadronic_Correl(q(0),qm,t,p(0),pr(3),dH)

       dLH=8*k(0)*(2*k(0) - q(0) - cosa*qm)*dH(0,0) + 8*k(0)*(-4*cosa*k(0) + 2*cosa*q(0) + 2*qm)*dH(0,3) + &
            & 8*k(0)*(2*k(0) - 2*cosa**2*k(0) - q(0) + cosa*qm)*dH(1,1) + 8*k(0)*((0,2)*cosa*q(0) - (0,2)*qm)*dH(1,2) + &
            & 8*k(0)*(2*cosa**2*k(0) - q(0) - cosa*qm)*dH(3,3)


    case (EM)

       call formfactors_QE(-t_in,process_ID,proton,F1_p,F2_p)
       call formfactors_QE(-t_in,process_ID,neutron,F1_n,F2_n)

       F1_s=F1_p+F1_n
       F1_a=F1_p-F1_n

       F2_s=F2_p+F2_n
       F2_a=F2_p-F2_n

       call Hadronic_Correl(q(0),qm,t,p(0),pr(3),dH)

       dLH=4*k(0)*(2*k(0) - q(0) - cosa*qm)*dH(0,0) + 4*k(0)*(-4*cosa*k(0) + 2*cosa*q(0) + 2*qm)*dH(0,3) + &
            & 4*k(0)*(2*k(0) - 2*cosa**2*k(0) - q(0) + cosa*qm)*dH(1,1) + 4*k(0)*(2*cosa**2*k(0) - q(0) - cosa*qm)*dH(3,3)


    case default
       write(*,*) 'correlations only implemented for EM, NC, CC, no anti-leptons'
       write(*,*) 'problem with process_ID in correlations -> STOP', process_ID
       stop

    end select

    !now transform back into GeV^4
    dLH=dLH*hbarc**4


  contains

    subroutine rotate(cosa,sina,p,pr)
      real,intent(in)::cosa,sina
      real,dimension(0:3),intent(in)::p
      real,dimension(0:3),intent(out)::pr
      pr(0)=p(0)
      pr(1)=cosa*p(1)-sina*p(3)
      pr(2)=p(2)
      pr(3)=sina*p(1)+cosa*p(3)
    end subroutine rotate

    subroutine Hadronic_Correl(q0,qm,q2,E,pz,dH)
      real,intent(in)::q0,qm,q2,E,pz
      complex,dimension(0:3,0:3),intent(out)::dH
      complex::uN,uD,u ! nucleon and Delta Lindhard functions and the sum
      complex,parameter::ui=(0.,1.)
      real::reuN,imuN  ! re and im parts of the nucleon Lindhard function with a gap
      real::vL,vT      ! longitudinal and transverse potentials
      real::CN,CL,CT,EN,DN   ! corrections due to correlations
      real::fpr,f0,g0
      real::D

      !initialize:
      reuN=0.
      imuN=0.
      vL=0.
      vT=0.
      uN=(0.,0.)
      uD=(0.,0.)

      pFave=(3.*pi**2*rhoave)**(0.33333)

      !Correlations (without gap)
      !call ulind(q0/mpi_fm,qm/mpi_fm,pFave/mpi_fm,uN,uD)
      !u=mpi_fm**2*(uN+uD)

      !Correlations (with gap)
      call ulind(q0/mpi_fm,qm/mpi_fm,pFave/mpi_fm,uN,uD)
      call drelin(q0,qm,pFave,5./hb,reuN)
      call dmalin(q0,qm,pFave,5./hb,imuN)
      uN=reuN+ui*imuN
      u=mpi_fm**2*uD+uN

      fpr=(rhoave/rho0ave)*fpin+(1.-rhoave/rho0ave)*fpex
      call potential(q0,qm,q2,vL,vT)
      CN=abs(1.-c0*fpr*uN)**(-2)
      CL=abs(1.-u*VL)**(-2)
      CT=abs(1.-u*VT)**(-2)

      dH=0.

      if (.not.withScalarInt) CN=1.


      select case (process_ID)

      case (CC)

         D=1./(mpi_fm**2-q2)
         dH(0,0)=4.*F1*(2.*F1*E**2-F2*qm**2)*(CN-1.)-8.*Fa**2*(q0)**2*D*(q2*D+2.)*mn_fm**2*(CL-1.)
         dH(0,3)=4.*E*(2.*pz+qm)*(F1**2*(CN-1.)+Fa**2*(CL-1.))-8.*Fa**2*q0*qm*D*(q2*D+2.)*mn_fm**2*(CL-1.)
         dH(3,3)=8.*Fa**2*(1.-qm**2*D*(q2*D+2))*mn_fm**2*(CL-1.)
         dH(1,1)=-2.*F2*(F2-2.*F1)*q2*(CT-1)+8.*Fa**2*mn_fm**2*(CT-1.)
         dH(1,2)=-8.*ui*Fa*(F1+F2)*qm*E*(CT-1.)


      case (NC)

         f0=(rhoave/rho0ave)*f0in+(1.-rhoave/rho0ave)*f0ex
         g0=(rhoave/rho0ave)*g0in+(1.-rhoave/rho0ave)*g0ex

         DN=abs(1.-c0*f0*uN)**(-2)
         EN=abs(1.-c0*g0*uN)**(-2)

         if (.not.withScalarInt) then
            DN=1.
            EN=1.
         end if

         dH(0,0)=4.*(F1_a**2*(CN-1.)+F1_s**2*(DN-1.))*E**2-2.*(F1_a*F2_a*(CN-1.)+F1_s*F2_s*(DN-1.))*qm**2
         dH(0,3)=(2.*F1_a**2*(CN-1.)+2.*F1_s**2*(DN-1.)+2.*Fa_a**2*(CL-1.)+2.*Fa_s**2*(EN-1.))*E*(2.*pz+qm)
         dH(3,3)=4.*(Fa_a**2*(CL-1.)+Fa_s**2*(EN-1.))*mn_fm**2
         dH(1,1)=-(F2_a**2*(CT-1.)+F2_s**2*(EN-1.))*q2-2.*(F1_a*F2_a*(CT-1.)+F1_s*F2_s*(EN-1.))*q2+ &
              & 4.*(Fa_a**2*(CT-1.)+Fa_s**2*(EN-1.))*mn_fm**2
         dH(1,2)=-4.*ui*(Fa_a*(F1_a+F2_a)*(CT-1.)+Fa_s*(F1_s+F2_s)*(EN-1.))*qm*E


      case (EM)

         f0=(rhoave/rho0ave)*f0in+(1.-rhoave/rho0ave)*f0ex
         g0=(rhoave/rho0ave)*g0in+(1.-rhoave/rho0ave)*g0ex

         DN=abs(1.-c0*f0*uN)**(-2)
         EN=abs(1.-c0*g0*uN)**(-2)

         if (.not.withScalarInt) then
            DN=1.
            EN=1.
         end if

         dH(0,0)=4.*(F1_a**2*(CN-1.)+F1_s**2*(DN-1.))*E**2-2.*(F1_a*F2_a*(CN-1.)+F1_s*F2_s*(DN-1.))*qm**2
         dH(0,3)=(2.*F1_a**2*(CN-1.)+2.*F1_s**2*(DN-1.))*E*(2.*pz+qm)
         dH(3,3)=0.
         dH(1,1)=-(F2_a**2*(CT-1.)+F2_s**2*(EN-1.))*q2-2.*(F1_a*F2_a*(CT-1.)+F1_s*F2_s*(EN-1.))*q2
         dH(1,2)=0.


      case default
         write(*,*) 'correlations only implemented for EM, NC, CC, no anti-leptons'
         write(*,*) 'problem with process_ID in correlations -> STOP', process_ID
         stop

      end select

    end subroutine Hadronic_Correl

    subroutine potential(Q0,Q,t,vL,vT)
      real,intent(in)::Q0,Q,t
      real,intent(out)::vL,vT
      real::ffpi,ffrho
      ffpi=(Lpi**2-mpi_fm**2)/(Lpi**2-t)
      ffrho=(Lrho**2-mrho_fm**2)/(Lrho**2-t)
      vL=(f/mpi_fm)**2*(Q**2/(t-mpi_fm**2)*ffpi**2+gp)
      vT=(f/mpi_fm)**2*(Q**2/(t-mrho_fm**2)*crho*ffrho**2+gp)
    end subroutine potential

    function sp(p1,p2)
      !Scalar product of two 3-momenta
      real, dimension(0:3), intent(in):: p1,p2
      real sp
      sp=p1(1)*p2(1)+p1(2)*p2(2)+p1(3)*p2(3)
    end function sp

    subroutine vector_prod(p1,p2,vp)
      !Vector product of two 3-momenta
      real, dimension(0:3), intent(in):: p1,p2
      real, dimension(0:3)::vp
      vp(1)=p1(2)*p2(3)-p1(3)*p2(2)
      vp(2)=p1(3)*p2(1)-p1(1)*p2(3)
      vp(3)=p1(1)*p2(2)-p1(2)*p2(1)
    end subroutine vector_prod

  end subroutine correlationsNieves


  !****************************************************************************
  !*      LINHARD FUNCTIONS
  !****************************************************************************
  SUBROUTINE ULIND(QZR,Q,XKF,CUFUN,cudel)
    use constants, only: pii => pi

    real, intent(in) :: QZR,Q,XKF
    complex, intent(out) :: CUFUN,cudel

    real :: XMN,XMD,WRES,FNS,FDS,QA,QZ,QZA,RUN,AM,AP,TERF,TERS,DCOCI,YMN
    real :: B,RO,FAC,A,TDIR,TCROS,SQS,QFR,GAMH
    complex :: CYI,CCROS,CDIR,CAC,CAPC
    complex :: cunuc


    XMN=6.7348   ! nucleon mass (in units of mpi)
    XMD=8.8315   ! delta mass        " "
    WRES=XMD-XMN
    FNS=1.
    FDS=4.52
    CYI=(0.,1.)
    !PII=3.141592
    QA=Q/XKF
    QZ=abs(QZR)
    QZA=abs(QZ*XMN/XKF**2)


    if (Q.EQ.0.) then
       RUN=FNS*XMN*XKF/PII**2*2./3.*(QA/QZA)**2
    else
       AM=QZA/QA-QA/2.
       AP=QZA/QA+QA/2.
       if (abs(1.-AM).LT.0.00001) AM=1.00001
       if (abs(1.-AP).LT.0.00001) AP=1.00001
       TERF=abs((1.+AM)/(1.-AM))+1.e-15
       TERS=abs((1.+AP)/(1.-AP))+1.e-15

       if (QZA-0.00001.le.0.) then
          RUN=FNS*XMN*XKF/PII**2*(-1.+(1.-AM**2)/(2. &
               & *QA)*ALOG(TERF)-(1.-AP**2)/(2.*QA)*ALOG(TERS))
       else
          DCOCI=abs(QA/QZA)

          if (DCOCI-0.1.le.0.) then
             RUN=FNS*XMN*XKF/PII**2*2./3.*(QA/QZA)**2
          else
             RUN=FNS*XMN*XKF/PII**2*(-1.+(1.-AM**2)/(2. &
                  & *QA)*ALOG(TERF)-(1.-AP**2)/(2.*QA)*ALOG(TERS))
          end if
       end if
    end if

    YMN=0.

    if ((QA**2/2.+QA).GE.QZA.AND.QZA.GE.(QA**2/2.-QA).AND.QA.GE.2.) &
         & YMN=-2.*XMN*XKF*FNS/(4.*PII*QA)*(1.-AM**2)

    if (QA.LT.2..AND.(QA+QA**2/2.).GE.QZA.AND.QZA.GE.(QA-QA** &
         & 2/2.))YMN=-2.*XMN*XKF*FNS/(4.*PII*QA)*(1.-AM**2)

    if (QA.LT.2..AND.0..LE.QZA.AND.QZA.LE.(QA-QA**2/2.)) &
         & YMN=-2.*XMN*XKF*FNS/(4.*PII*QA)*2.*QZA

    CUNUC=CMPLX(RUN,YMN)

    B=Q/XMD
    RO=XKF**3*2./(3.*PII**2)

    if (Q.EQ.0.) then
       FAC=0.
    else
       FAC=(4./3.*XKF/(2.*PII))**2/B**3*FDS
    end if

    if (abs(QZ)-1..le.0.) then
       A=(QZ-WRES-Q**2/(2.*XMD))/XKF
       AP=(-QZ-WRES-Q**2/(2.*XMD))/XKF

       if (abs(B/A)-0.1.le.0.) then
          TDIR=4./9.*FDS*RO/(A*XKF)
       else
          TERF=abs((A+B)/(A-B))+1.e-15
          TDIR=FAC*(B*A+(B**2-A**2)/2.*ALOG(TERF))
       end if

       if (abs(B/AP)-0.1.le.0.) then
          TCROS=4./9.*FDS*RO/(AP*XKF)
       else
          TERS=abs((AP-B)/(AP+B))+1.e-15
          TCROS=FAC*(B*AP-(B**2-AP**2)/2.*ALOG(TERS))
       end if

       CUDEL=TDIR+TCROS
    else
       SQS=XMN+abs(QZ)
       QFR=SQRT((QZ**2-1.)/(1.+QZ/XMN))
       GAMH=1./(3.*4.*PII)*FDS*XMN/SQS*QFR**3

       if (QZ.GT.0.) CAC=(QZ-WRES+CYI*GAMH-Q**2/(2.*XMD))/XKF
       if (QZ.GT.0.) CAPC=(-QZ-WRES-Q**2/(2.*XMD))/XKF
       if (QZ.LE.0.) CAC=(QZ-WRES-Q**2/(2.*XMD))/XKF
       if (QZ.LE.0.) CAPC=(-QZ-WRES+CYI*GAMH-Q**2/(2.*XMD))/XKF

       if (CABS(B/CAC)-0.1.le.0) then
          CDIR=4./9.*FDS*RO/(CAC*XKF)
       else
          CDIR=FAC*(B*CAC+(B**2-CAC**2)/2.*CLOG((CAC+B)/(CAC-B)))
       end if

       if (CABS(B/CAPC)-0.1.le.0) then
          CCROS=4./9.*FDS*RO/(CAPC*XKF)
       else
          CCROS=FAC*(B*CAPC-(B**2-CAPC**2)/2.*CLOG((CAPC-B)/(CAPC+B)))
       end if

       CUDEL=CDIR+CCROS
    end if
    !       CUFUN=CUNUC+CUDEL
    cufun=cunuc
  END SUBROUTINE ULIND

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE Dmalin(QU0,DQ,DKF,DGAP,DLINIM)
    use constants, only: dpi => pi

    real, intent(in) :: QU0,DQ,DKF,DGAP
    real, intent(out) :: DLINIM

    real :: DMPI,DMNU,DKF2,DZ,DZPRI,DC1,DC2,DFACTOR,Q0 ! ,DELTA, DMNU2,DNU, DX
    DMPI=0.70729
    DMNU=6.735*DMPI
    !DPI=3.14159
    DKF2=DKF*DKF
    !DMNU2=DMNU*DMNU
    !DNU=2.*DMNU*QU0/DKF2
    !DX=DQ/DKF
    !DELTA=2.*DMNU*DGAP/DKF2

    if (QU0.GT.DGAP) THEN
       Q0=QU0-DGAP
    else if ((QU0.GT.0.).AND.(QU0.LT.DGAP)) then
       DLINIM=0.
       return
    else
       Q0=-QU0
    end if

    if (DQ.LT.1.e-7) then
       DLINIM=0.
       return
    end if

    DZ=abs((DMNU/(DQ*DKF))*(Q0-DQ*DQ/(2.*DMNU)))
    DZPRI=abs((DMNU/(DQ*DKF))*(-Q0-DQ*DQ/(2.*DMNU)))
    if (DZ.LT.1.) then
       DC1=(1.-DZ**2)
    else
       DC1=0.
    end if
    if (DZPRI.LT.1.) then
       DC2=(1.-DZPRI**2)
    else
       DC2=0.
    end if

    DFACTOR=-(1./2.)*DKF2*DMNU*Q0/(DPI*DQ*ABS(Q0))

    DLINIM=DFACTOR*(DC1-DC2)

  END SUBROUTINE Dmalin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE DRELIN(Q0,DQ,DKF,DGAP,DREALIN)
    use constants, only: dpi => pi

    real, intent(in) :: Q0,DQ,DKF,DGAP
    real, intent(out) :: DREALIN

    real :: DMPI,DMNU,DKF2,DNU,DX,DELTA ! ,DMNU2
    real :: DC1,DC2,DC3,DC4,DC5,DC6
    real :: DFACTOR

    DMPI=0.70729
    DMNU=6.735*DMPI
    !DPI=3.14159

    if (DQ.LT.1.e-5) then
       DREALIN=0.
       return
    end if

    DKF2=DKF*DKF
    !DMNU2=DMNU*DMNU
    DNU=2.*DMNU*Q0/DKF2
    DX=DQ/DKF
    DELTA=2.*DMNU*DGAP/DKF2

    DC1=abs((DNU-DELTA+DX*DX-2.*DX)/(DNU-DELTA))
    DC2=abs((DNU+DELTA-DX*DX+2.*DX)/(DNU+DELTA))
    DC3=abs((DNU-DELTA-DX*DX-2.*DX)/(DNU-DELTA+DX*DX-2.*DX))
    DC4=abs((DNU+DELTA+DX*DX+2.*DX)/(DNU+DELTA-DX*DX+2.*DX))
    DC5=abs((DNU-DELTA-DX*DX-2.*DX)/(DNU-DELTA-DX*DX+2.*DX))
    DC6=abs((DNU+DELTA+DX*DX+2.*DX)/(DNU+DELTA+DX*DX-2.*DX))

    DFACTOR=-2.*DMNU*DKF/(DPI*DPI*2.*DX)
    if (DX.LE.2.) then
       DREALIN=DFACTOR*(DX+0.5*DELTA+0.5*(DNU-DELTA)*LOG(DC1)-0.5*(DNU+DELTA)*   &
            &  LOG(DC2)+0.5*(1.-0.25*(((DNU-DELTA)/DX)-DX)**2)*LOG(DC3)+0.5*     &
            &  (1.-0.25*(((DNU+DELTA)/DX)+DX)**2)*LOG(DC4))

    else if (DX.GE.2.) then
       DREALIN=DFACTOR*(DX+(DELTA/DX)+0.5*(1.-0.25*(((DNU-DELTA)/DX)-DX)**2)*LOG(DC5)+ &
            &  0.5*(1.-0.25*(((DNU+DELTA)/DX)+DX)**2)*LOG(DC6))
    else
       DREALIN=0.
    end if

  END SUBROUTINE DRELIN


!!$
!!$
!!$  !*********************************************************************
!!$  !*      LINHARD FUNCTIONS
!!$  !**********************************************************************
!!$  SUBROUTINE ULIND(QZR,Q,XKF,CUFUN,cudel)
!!$    !       SUBROUTINE ULIND(QZR,Q,XKF,CUFUN)
!!$    implicit real(a-b,d-h,o-z)
!!$    IMPLICIT COMPLEX (C)
!!$    XMN=6.7348   ! nucleon mass (in units of mpi)
!!$    XMD=8.8315   ! delta mass        " "
!!$    WRES=XMD-XMN
!!$    FNS=1.
!!$    FDS=4.52
!!$    CYI=(0.,1.)
!!$    PII=3.141592
!!$    QA=Q/XKF
!!$    QZ=ABS(QZR)
!!$    QZA=ABS(QZ*XMN/XKF**2)
!!$    IF(Q.EQ.0.) GO TO 15
!!$    AM=QZA/QA-QA/2.
!!$    AP=QZA/QA+QA/2.
!!$    IF(ABS(1.-AM).LT.0.00001) AM=1.00001
!!$    IF(ABS(1.-AP).LT.0.00001) AP=1.00001
!!$    TERF=ABS((1.+AM)/(1.-AM))+1.e-15
!!$    TERS=ABS((1.+AP)/(1.-AP))+1.e-15
!!$    IF(QZA-0.00001) 16,16,18
!!$18  DCOCI=ABS(QA/QZA)
!!$    IF(DCOCI-0.1) 15,15,16
!!$15  RUN=FNS*XMN*XKF/PII**2*2./3.*(QA/QZA)**2
!!$    GO TO 17
!!$16  RUN=FNS*XMN*XKF/PII**2*(-1.+(1.-AM**2)/(2. &
!!$         & *QA)*ALOG(TERF)-(1.-AP**2)/(2.*QA)*ALOG(TERS))
!!$17  YMN=0.
!!$    IF((QA**2/2.+QA).GE.QZA.AND.QZA.GE.(QA**2/2.-QA).AND.QA.GE.2.) &
!!$         & YMN=-2.*XMN*XKF*FNS/(4.*PII*QA)*(1.-AM**2)
!!$    IF(QA.LT.2..AND.(QA+QA**2/2.).GE.QZA.AND.QZA.GE.(QA-QA** &
!!$         & 2/2.))YMN=-2.*XMN*XKF*FNS/(4.*PII*QA)*(1.-AM**2)
!!$    IF(QA.LT.2..AND.0..LE.QZA.AND.QZA.LE.(QA-QA**2/2.)) &
!!$         & YMN=-2.*XMN*XKF*FNS/(4.*PII*QA)*2.*QZA
!!$    CUNUC=CMPLX(RUN,YMN)
!!$    B=Q/XMD
!!$    RO=XKF**3*2./(3.*PII**2)
!!$    IF(Q.EQ.0.) GOTO 250
!!$    FAC=(4./3.*XKF/(2.*PII))**2/B**3*FDS
!!$250 IF(ABS(QZ)-1.) 2,2,3
!!$2   A=(QZ-WRES-Q**2/(2.*XMD))/XKF
!!$    AP=(-QZ-WRES-Q**2/(2.*XMD))/XKF
!!$    IF(ABS(B/A)-0.1) 25,25,26
!!$25  TDIR=4./9.*FDS*RO/(A*XKF)
!!$    GO TO 27
!!$26  TERF=ABS((A+B)/(A-B))+1.e-15
!!$    TDIR=FAC*(B*A+(B**2-A**2)/2.*ALOG(TERF))
!!$27  IF(ABS(B/AP)-0.1) 28,28,29
!!$28  TCROS=4./9.*FDS*RO/(AP*XKF)
!!$    GO TO 30
!!$29  TERS=ABS((AP-B)/(AP+B))+1.e-15
!!$    TCROS=FAC*(B*AP-(B**2-AP**2)/2.*ALOG(TERS))
!!$30  CUDEL=TDIR+TCROS
!!$    GO TO 10
!!$3   SQS=XMN+ABS(QZ)
!!$    QFR=SQRT((QZ**2-1.)/(1.+QZ/XMN))
!!$    GAMH=1./(3.*4.*PII)*FDS*XMN/SQS*QFR**3
!!$    IF(QZ.GT.0.) CAC=(QZ-WRES+CYI*GAMH-Q**2/(2.*XMD))/XKF
!!$    IF(QZ.GT.0.)CAPC=(-QZ-WRES-Q**2/(2.*XMD))/XKF
!!$    IF(QZ.LE.0.) CAC=(QZ-WRES-Q**2/(2.*XMD))/XKF
!!$    IF(QZ.LE.0.) CAPC=(-QZ-WRES+CYI*GAMH-Q**2/(2.*XMD))/XKF
!!$    IF(CABS(B/CAC)-0.1) 35,35,36
!!$35  CDIR=4./9.*FDS*RO/(CAC*XKF)
!!$    GO TO 37
!!$36  CDIR=FAC*(B*CAC+(B**2-CAC**2)/2.*CLOG((CAC+B)/(CAC-B)))
!!$37  IF(CABS(B/CAPC)-0.1) 38,38,39
!!$38  CCROS=4./9.*FDS*RO/(CAPC*XKF)
!!$    GO TO 40
!!$39  CCROS=FAC*(B*CAPC-(B**2-CAPC**2)/2.*CLOG((CAPC-B)/(CAPC+B)))
!!$40  CUDEL=CDIR+CCROS
!!$10  CONTINUE
!!$    !       CUFUN=CUNUC+CUDEL
!!$    cufun=cunuc
!!$    RETURN
!!$  END SUBROUTINE ULIND



end module matrixElementQE
