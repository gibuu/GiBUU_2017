module electronNuc_to_Delta

  implicit none

  Public :: elepi

  private
  logical :: debug=.false.

  real,parameter::pi=3.1416
  complex,parameter::ui=(0.,1.)
  real,parameter::hc=1000.
  real,parameter::hb=197.33    ! (MeV x fm)
  real,parameter::hbg=hb/hc    ! (GeV x fm)

  real,parameter::alpha=1./137.

  !       nucleon
  real,parameter::mn=939./hc

  !       pion
  real,parameter::mpi=139.6/hc

  !       delta
  real,parameter::mdel=1232./hc  ! mass
  real,parameter::fs=2.13        ! DeltaNpi coupling

  real,parameter::Lam=1.         ! cut off mass for the s-channel form factor


  integer, parameter::nfs =0    ! swich for this form factor: 0(1) <-> without(with) fs

contains

  subroutine elepi(n,e,ef,thf,thpi,phipi,cs)
    use particleDefinition
    use electronPionProduction_kine, only :getKinematics
    use degRad_conversion, only : degrees
    use idTable, only : nucleon
    use constants, only: mN

    integer,intent(in)::n   ! defines the isospin factor 
    ! n=1: e + p -> e + p pi0
    ! n=2  e + p -> e + n pi+
    ! n=3  e + n -> e + n pi0
    ! n=4  e + n -> e + p pi-


    real,intent(in)::e      ! incoming electron energy [GeV]
    real,intent(in)::ef     ! outgoing electron energy [GeV]
    real,intent(in)::thf    ! Theta angle of the outgoing electron (its phi angle is set to zero, ie. zx is the scattering plane). [Radians]
    real,intent(in)::thpi   ! Theta angle of the outgoing pion [Radians]
    real,intent(in)::phipi  ! Phi angle of the outgoing pion [Radians]

    real,intent(out)::cs    ! 5-d cross section dsigma/d_ef/dOmega_ef/dOmega_pi [nb/(GeV str^2)]

    !        Other variables
    real::ti                ! isospin factor
    real,dimension(0:3)::k  ! incoming electron 4-momentum
    real,dimension(0:3)::kf ! outgoing electron 4-momentum
    real,dimension(0:3)::p  ! initial nucleon 4-momentum
    real,dimension(0:3)::pf ! final nucleon 4-momentum
    real,dimension(0:3)::ppi! pion 4-momentum 
    real,dimension(0:3)::q  ! transfered 4-momentum   
    
    real,dimension(1:2,0:3)::pf_arr,ppi_arr 
  

    logical :: success
    logical,parameter :: original=.true.

    type(particle) :: initNuc
    integer :: pionCharge,numsolutions

    select case (n)
    case (1)
       ti=sqrt(2./3.)
    case (2)
       ti=-sqrt(1./3.)
    case (3)
       ti=sqrt(2./3.)
    case (4)               
       ti=sqrt(1./3.)
    end select

    if(original) then
       call getKinematics_local(e, ef,thf,phipi,thpi,k,p,kf,pf,ppi,q,success)
    else
       select case (n)
       case (1)
          pionCharge=0
          initNuc%charge=1
       case (2)
          pionCharge=1
          initNuc%charge=1
       case (3)
          pionCharge=0
          initNuc%charge=0
       case (4)
          pionCharge=-1
          initNuc%charge=0
       end select
       initNuc%ID=nucleon
       initNuc%momentum=(/mN,0.,0.,0./)
       call getKinematics(pionCharge,initNuc,e,ef,degrees(thf),degrees(phipi),degrees(thpi),  &
            &  k,kf,ppi_arr,pf_arr,success,numSolutions)
       p=initNuc%momentum
    end if
    
    pf=pf_arr(1,:)
    ppi=ppi_arr(1,:) 

    if(success) then
       call crosec(cs,k,p,kf,pf,ppi,q,ti)
    else 
       cs=0.
    end if
  end subroutine elepi


  subroutine crosec(cs_s,k,p,kf,pf,ppi,q,ti)

    !           real,intent(in)::ppim   ! modulus of the pion 3-momentum
    real,intent(out)::cs_s
    real,intent(in),dimension(0:3)::k  ! incoming electron 4-momentum
    real,intent(in),dimension(0:3)::kf ! outgoing electron 4-momentum
    real,intent(in),dimension(0:3)::p  ! initial nucleon 4-momentum
    real,intent(in),dimension(0:3)::pf ! final nucleon 4-momentum
    real,intent(in),dimension(0:3)::ppi! pion 4-momentum 
    real,intent(in),dimension(0:3)::q  ! transfered 4-momentum   
    real, intent(in) :: ti


    real :: cosa 
    complex::gd    ! Delta propagator
    real::a,t
    real::ampl2    ! amplitude squared, summed and averaged over spins
    real :: ppim, q_abs
    real,dimension(0:3)::pdel  ! delta 4-momentum   
    real::qm                ! modulus of the 3-momentum tranfer

    real::W                 ! Delta invariant mass

    real::s!,t               ! Mandelstam variables for the 2-body part e + N -> e + Delta


    qm=sqrt(sp(q,q))

    s=sp4(p+k,p+k)
    if (s.lt.(mn+mpi)**2) then 
       write(*,*)'below pi production threshold' 
       stop
    end if

    ppim=sqrt(Dot_product(ppi(1:3),ppi(1:3)))
    q_abs=sqrt(Dot_product(q(1:3),q(1:3)))


    if(debug) then
       write(*,'(A3,4F15.6)') 'k',k  
       write(*,'(A3,4F15.6)') 'p',p
       write(*,'(A3,4E15.6)') 'kf',kf
       write(*,'(A3,4F15.6)') 'pf',pf
       write(*,'(A3,4F15.6)') 'ppi',ppi
       write(*,'(A3,4F15.6)') 'q',ppi
    end if
    !        angle between q and ppi

    if(q_abs.gt.0.and.ppim.gt.0) then
       cosa=Dot_Product(q(1:3),ppi(1:3))/ppim/q_abs
    else
       cosa=0
    end if
    pdel=pf+ppi

    W=sqrt(sp4(pdel,pdel))

    if (W.lt.(mn+mpi)) then 
       write(*,*)'Delta invariant mass too small W=',W
       stop
    end if

    gd=1./(W**2-mdel**2+ui*W*gamd(W))


    t=sp4(q,q)
    a=(ti*4.*pi*alpha)**2*(fs/mpi)**2*(abs(gd))**2/t**2*ffs(W)**2
    ampl2=1./4.*a*LH(ppi,pdel,p,k,kf,pf,W,t,s)
    cs_s=1./(2.*pi)**5*kf(0)*ppim**3/(4.*8.*k(0)*mn)/abs(ppim**2*pf(0)+ppi(0)*(ppim**2-qm*ppim*cosa))*ampl2    ! in GeV(-3)
    cs_s=cs_s*(hbg)**2*1.e7     ! in nb/GeV/str**2

  end subroutine crosec


  real function LH(ppi,pdel,p,k,kf,pf,W,t,s)
    use FF_Delta_production

    real,dimension(0:3),intent(in)::ppi,pdel,p,k,pf,kf
    real,intent(in)::W, t,s
    real::spr1,spr2,spr3,spr4
    real::C3v,C4v,C5v,c6v

    call formfactors_Delta(-t,W,1,c3v,c4v,c5v,c6v)

    spr1=sp4(pdel,ppi)
    spr2=sp4(p,ppi)
    spr3=sp4(k,ppi)
    spr4=sp4(kf,ppi)

    LH=(2*(W*(mn + W) - spr1)*(2*C3v**2*mn**2*(24*mn**3*s*W**3*spr2**2 - 24*mn*s**2*W**3*spr2**2 - 24*mn*s*t*W**3*spr2**2 - &
         & 12*mn**2*t*W**4*spr2**2 + 12*t**2*W**4*spr2**2 - 24*mn**3*W**5*spr2**2 + 24*mn*s*W**5*spr2**2 + &
         & 12*t*W**6*spr2**2 + 12*mn**4*t*W**2*spr2*spr1 - 24*mn**2*s*t*W**2*spr2*spr1 + 24*s**2*t*W**2*spr2*spr1 - &
         & 24*mn**2*t**2*W**2*spr2*spr1 + 24*s*t**2*W**2*spr2*spr1 + 12*t**3*W**2*spr2*spr1 - 48*mn**3*s*W**3*spr2*spr1 + &
         & 48*mn*s**2*W**3*spr2*spr1 + 12*mn**3*t*W**3*spr2*spr1 + 48*mn*s*t*W**3*spr2*spr1 - 12*mn*t**2*W**3*spr2*spr1 + &
         & 12*mn**2*t*W**4*spr2*spr1 - 24*s*t*W**4*spr2*spr1 - 12*t**2*W**4*spr2*spr1 + 48*mn**3*W**5*spr2*spr1 - &
         & 48*mn*s*W**5*spr2*spr1 - 12*mn*t*W**5*spr2*spr1 - 4*mn**6*t*spr1**2 + 8*mn**4*s*t*spr1**2 - 8*mn**2*s**2*t*spr1**2 + &
         & 12*mn**4*t**2*spr1**2 - 16*mn**2*s*t**2*spr1**2 + 8*s**2*t**2*spr1**2 - 12*mn**2*t**3*spr1**2 + 8*s*t**3*spr1**2 + &
         & 4*t**4*spr1**2 - 12*mn**4*t*W**2*spr1**2 + 32*mn**2*s*t*W**2*spr1**2 - 24*s**2*t*W**2*spr1**2 +&
         &  16*mn**2*t**2*W**2*spr1**2 - &
         & 32*s*t**2*W**2*spr1**2 - 4*t**3*W**2*spr1**2 + 24*mn**3*s*W**3*spr1**2 - 24*mn*s**2*W**3*spr1**2 -&
         &  12*mn**3*t*W**3*spr1**2 - &
         & 24*mn*s*t*W**3*spr1**2 + 20*mn*t**2*W**3*spr1**2 - 4*mn**2*t*W**4*spr1**2 + 24*s*t*W**4*spr1**2 + &
         & 12*t**2*W**4*spr1**2 - &
         & 24*mn**3*W**5*spr1**2 + 24*mn*s*W**5*spr1**2 + 12*mn*t*W**5*spr1**2 - 12*t*W**6*spr1**2 + 12*W**3*(mn + W)*spr4* &
         & (-2*(mn**4 - mn**2*(s + 2*t) - mn**3*W + mn*(s + t)*W + t*(s + t - W**2))*spr2 + &
         & (2*mn**4 - mn**2*(2*s + 3*t) - 2*mn**3*W + 2*mn*(s + t)*W + t*(t - W**2))*spr1) + 12*W**3*(mn + W)*spr3* &
         & (-2*(-t + (mn - W)**2)*(-t + mn*(mn + W))*spr4 - 2*(s*t + mn*W*(s - W**2) + mn**2*(-s + W**2))*spr2 + &
         & (2*mn*W*(s - W**2) + t*(-t + W**2) + mn**2*(-2*s + t + 2*W**2))*spr1) + &
         & t*W**2*((mn**2 - t)*(mn**4 + 2*s**2 + 2*s*t + t**2 - 2*mn**2*(s + t)) - 12*mn*s*(mn**2 - s - t)*W + &
         & (9*mn**4 + 18*s**2 + 20*s*t + t**2 - 10*mn**2*(2*s + t))*W**2 + 4*mn*(3*mn**2 - 3*s - 5*t)*W**3 + &
         & (mn**2 - 9*(2*s + t))*W**4 + 9*W**6)*mpi**2) - 2*C3v*C4v*mn*W*(12*mn**4*s*W**2*spr2**2 - &
         & 12*mn**2*s**2*W**2*spr2**2 - 24*mn**2*s*t*W**2*spr2**2 + 12*s**2*t*W**2*spr2**2 + 12*s*t**2*W**2*spr2**2 - &
         & 24*mn**3*s*W**3*spr2**2 + 24*mn*s**2*W**3*spr2**2 + 24*mn*s*t*W**3*spr2**2 - 12*mn**4*W**4*spr2**2 + & 
         & 24*mn**2*s*W**4*spr2**2 - &
         & 12*s**2*W**4*spr2**2 + 24*mn**2*t*W**4*spr2**2 - 24*s*t*W**4*spr2**2 - 12*t**2*W**4*spr2**2 + 24*mn**3*W**5*spr2**2 - &
         & 24*mn*s*W**5*spr2**2 - 12*mn**2*W**6*spr2**2 + 12*s*W**6*spr2**2 - 12*t*W**6*spr2**2 - 24*mn**4*s*W**2*spr2*spr1 + & 
         & 24*mn**2*s**2*W**2*spr2*spr1 - &
         & 6*mn**4*t*W**2*spr2*spr1 + 72*mn**2*s*t*W**2*spr2*spr1 - 48*s**2*t*W**2*spr2*spr1 + 12*mn**2*t**2*W**2*spr2*spr1 - &
         & 48*s*t**2*W**2*spr2*spr1 - 6*t**3*W**2*spr2*spr1 + 48*mn**3*s*W**3*spr2*spr1 - 48*mn*s**2*W**3*spr2*spr1 - & 
         & 24*mn**3*t*W**3*spr2*spr1 - &
         & 48*mn*s*t*W**3*spr2*spr1 + 24*mn*t**2*W**3*spr2*spr1 + 24*mn**4*W**4*spr2*spr1 - 48*mn**2*s*W**4*spr2*spr1 + & 
         
         & 24*s**2*W**4*spr2*spr1 - &
         & 36*mn**2*t*W**4*spr2*spr1 + 72*s*t*W**4*spr2*spr1 + 12*t**2*W**4*spr2*spr1 - 48*mn**3*W**5*spr2*spr1 + &
         &  48*mn*s*W**5*spr2*spr1 + &
         & 24*mn*t*W**5*spr2*spr1 + 24*mn**2*W**6*spr2*spr1 - 24*s*W**6*spr2*spr1 - 6*t*W**6*spr2*spr1 + 4*mn**6*t*spr1**2 - & 
         & 8*mn**4*s*t*spr1**2 + &
         & 8*mn**2*s**2*t*spr1**2 - 12*mn**4*t**2*spr1**2 + 16*mn**2*s*t**2*spr1**2 - 8*s**2*t**2*spr1**2 +& 
         &  12*mn**2*t**3*spr1**2 - 8*s*t**3*spr1**2 - &
         & 4*t**4*spr1**2 - 4*mn**5*t*W*spr1**2 + 8*mn**3*s*t*W*spr1**2 - 8*mn*s**2*t*W*spr1**2 + & 
         & 8*mn**3*t**2*W*spr1**2 - 8*mn*s*t**2*W*spr1**2 - &
         & 4*mn*t**3*W*spr1**2 + 12*mn**4*s*W**2*spr1**2 - 12*mn**2*s**2*W**2*spr1**2 + 6*mn**4*t*W**2*spr1**2 - &
         & 56*mn**2*s*t*W**2*spr1**2 + &
         & 36*s**2*t*W**2*spr1**2 + 44*s*t**2*W**2*spr1**2 - 6*t**3*W**2*spr1**2 - 24*mn**3*s*W**3*spr1**2 +&
         &  24*mn*s**2*W**3*spr1**2 + &
         & 24*mn**3*t*W**3*spr1**2 + 32*mn*s*t*W**3*spr1**2 - 32*mn*t**2*W**3*spr1**2 - 12*mn**4*W**4*spr1**2 + &
         & 24*mn**2*s*W**4*spr1**2 - 12*s**2*W**4*spr1**2 + &
         & 16*mn**2*t*W**4*spr1**2 - 48*s*t*W**4*spr1**2 - 8*t**2*W**4*spr1**2 + 24*mn**3*W**5*spr1**2 -  &
         & 24*mn*s*W**5*spr1**2 - 28*mn*t*W**5*spr1**2 - &
         & 12*mn**2*W**6*spr1**2 + 12*s*W**6*spr1**2 + 18*t*W**6*spr1**2 - 6*W**2*spr3* &
         & (2*(-t + (mn - W)**2)*(-mn**2 + t + W**2)**2*spr4 - &
         & 2*(s*(mn**2 - t)**2 + 2*mn*s*(-mn**2 + t)*W + (-mn**4 + mn**2*t + 2*s*t)*W**2 + 2*mn*(mn**2 + s - t)*W**3 - &
         & (s + t)*W**4 - 2*mn*W**5 + W**6)*spr2 + &
         & (-t**3 + 2*s*t*W**2 - (2*s + t)*W**4 + 2*W**6 + mn**4*(2*s - t - 2*W**2) + 4*mn*W*(s - W**2)*(t + W**2) + &
         & 2*mn**2*t*(-s + t + W**2) + mn**3*(-4*s*W + 4*W**3))*spr1) - 6*W**2*spr4* &
         & (2*((mn**2 - t)**2*(mn**2 - s - t) - 2*mn*(mn**2 - t)*(mn**2 - s - t)*W + (mn**2 - 2*s - t)*t*W**2 + &
         & 2*mn*(mn**2 - s - t)*W**3 - (mn**2 - s - 2*t)*W**4)*spr2 + &
         & (-2*mn**6 + mn**4*(2*s + 5*t) + 4*mn**5*W + 4*mn*(s + t)*W*(t + W**2) - 4*mn**3*W*(s + 2*t + W**2) + &
         & (t - W**2)*(t**2 + (2*s + 3*t)*W**2) - 2*mn**2*(t*(s + 2*t) + t*W**2 - W**4))*spr1) - &
         & t*W**2*((mn**2 - t)*(mn**4 - 4*s**2 + mn**2*(4*s - 2*t) - 4*s*t + t**2) - 2*mn*(5*mn**4 + 4*s**2 +   &
         & 4*s*t + 5*t**2 - 2*mn**2*(2*s + 5*t))*W + &
         & (3*mn**4 + 12*s**2 + 8*s*t - 9*t**2 + mn**2*(-8*s + 6*t))*W**2 + 4*mn*(3*mn**2 + 2*s - 5*t)*W**3 - &
         & (5*mn**2 + 12*s - t)*W**4 - 10*mn*W**5 + 9*W**6)*mpi**2) + 2*C4v*C5v*(-t +  &
         &  (mn - W)**2)*W**2*(-12*mn**2*s*W**2*spr2**2 + &
         & 12*s**2*W**2*spr2**2 - 6*mn**2*t*W**2*spr2**2 + 12*s*t*W**2*spr2**2 + 6*t**2*W**2*spr2**2 +   &
         & 12*mn**2*W**4*spr2**2 - 12*s*W**4*spr2**2 + &
         & 6*t*W**4*spr2**2 + 3*mn**4*t*spr2*spr1 - 12*mn**2*s*t*spr2*spr1 + 12*s**2*t*spr2*spr1 - 6*mn**2*t**2*spr2*spr1 + &
         & 12*s*t**2*spr2*spr1 + 3*t**3*spr2*spr1 + 24*mn**2*s*W**2*spr2*spr1 - 24*s**2*W**2*spr2*spr1 +   &
         & 6*mn**2*t*W**2*spr2*spr1 - &
         & 36*s*t*W**2*spr2*spr1 - 6*t**2*W**2*spr2*spr1 - 24*mn**2*W**4*spr2*spr1 + 24*s*W**4*spr2*spr1 + 3*t*W**4*spr2*spr1 - &
         & 5*mn**4*t*spr1**2 + 16*mn**2*s*t*spr1**2 - 16*s**2*t*spr1**2 - 16*s*t**2*spr1**2 + 5*t**3*spr1**2 -   &
         & 12*mn**2*s*W**2*spr1**2 + &
         & 12*s**2*W**2*spr1**2 + 28*s*t*W**2*spr1**2 + 6*t**2*W**2*spr1**2 + 12*mn**2*W**4*spr1**2 -   &
         & 12*s*W**4*spr1**2 - 11*t*W**4*spr1**2 + &
         & 3*spr3*(2*W**2*(2*(-t**2 + (mn - W)**2*(mn + W)**2)*spr4 + &
         & (t**2 + (2*s + t)*W**2 - 2*W**4 - mn**2*(2*s + t - 2*W**2))*spr2) + &
         & ((mn**2 - t)*(mn**2 - 2*s - t)*t + 2*((s - t)*t + mn**2*(2*s + t))*W**2 - (4*(mn**2 + s) + 3*t)*W**4 + 4*W**6)*spr1) + &
         & 3*spr4*(-2*W**2*(-2*mn**4 + t**2 - (2*s + t)*W**2 + mn**2*(2*s + t + 2*W**2))*spr2 + &
         & (mn**4*(t - 4*W**2) - 2*mn**2*(t - 2*W**2)*(s + t + W**2) + (t - W**2)*(t*(2*s + t) + (4*s + 3*t)*W**2))*spr1) + &
         & t*W**2*(5*mn**4 + 4*s**2 + 4*s*t - 5*t**2 - 4*s*W**2 + 5*W**4 - 2*mn**2*(2*s + 3*W**2))*mpi**2) + &
         & C4v**2*(-t + (mn - W)**2)*W**2*(-12*mn**2*s*W**2*spr2**2 + 12*s**2*W**2*spr2**2 + 12*s*t*W**2*spr2**2 + &
         & 12*mn**2*W**4*spr2**2 - 12*s*W**4*spr2**2 + 12*W**2*(-mn**2 + t + W**2)*spr3*((-mn**2 + t + W**2)*spr4 + &
         & (s - W**2)*(spr2 - spr1)) + 12*(mn**2 - s - t)*W**2*(mn**2 - t - W**2)*spr4*(spr2 - spr1) + &
         & 24*mn**2*s*W**2*spr2*spr1 - 24*s**2*W**2*spr2*spr1 - 12*mn**2*t*W**2*spr2*spr1 - 24*s*t*W**2*spr2*spr1 + &
         & 12*t**2*W**2*spr2*spr1 - 24*mn**2*W**4*spr2*spr1 + 24*s*W**4*spr2*spr1 + 12*t*W**4*spr2*spr1 - 2*mn**4*t*spr1**2 + &
         & 4*mn**2*s*t*spr1**2 - 4*s**2*t*spr1**2 + 4*mn**2*t**2*spr1**2 - 4*s*t**2*spr1**2 - 2*t**3*spr1**2 - &
         & 12*mn**2*s*W**2*spr1**2 + 12*s**2*W**2*spr1**2 + 12*mn**2*t*W**2*spr1**2 + 16*s*t*W**2*spr1**2 - 16*t**2*W**2*spr1**2 + &
         & 12*mn**2*W**4*spr1**2 - 12*s*W**4*spr1**2 - 14*t*W**4*spr1**2 + t*W**2*(5*mn**4 + 4*s**2 + 4*s*(t - W**2) + &
         & 5*(t + W**2)**2 - 2*mn**2*(2*s + 5*t + 3*W**2))*mpi**2) + &
         & C5v**2*(-t + (mn - W)**2)*(-12*mn**2*s*W**4*spr2**2 + 12*s**2*W**4*spr2**2 - 12*mn**2*t*W**4*spr2**2 + &
         & 12*s*t*W**4*spr2**2 - 12*t**2*W**4*spr2**2 + 12*mn**2*W**6*spr2**2 - 12*s*W**6*spr2**2 + 12*t*W**6*spr2**2 + &
         & 6*mn**4*t*W**2*spr2*spr1 - 24*mn**2*s*t*W**2*spr2*spr1 + 24*s**2*t*W**2*spr2*spr1 + 24*s*t**2*W**2*spr2*spr1 - &
         & 6*t**3*W**2*spr2*spr1 + 24*mn**2*s*W**4*spr2*spr1 - 24*s**2*W**4*spr2*spr1 + 24*mn**2*t*W**4*spr2*spr1 - &
         & 48*s*t*W**4*spr2*spr1 + 12*t**2*W**4*spr2*spr1 - 24*mn**2*W**6*spr2*spr1 + 24*s*W**6*spr2*spr1 - 6*t*W**6*spr2*spr1 - &
         & 16*mn**2*s*t**2*spr1**2 + 16*s**2*t**2*spr1**2 + 16*s*t**3*spr1**2 - 8*mn**4*t*W**2*spr1**2 +   &
         & 28*mn**2*s*t*W**2*spr1**2 - &
         & 28*s**2*t*W**2*spr1**2 - 44*s*t**2*W**2*spr1**2 - 8*t**3*W**2*spr1**2 - 12*mn**2*s*W**4*spr1**2 + &
         & 12*s**2*W**4*spr1**2 - 12*mn**2*t*W**4*spr1**2 + 40*s*t*W**4*spr1**2 + 16*t**2*W**4*spr1**2 +   &
         & 12*mn**2*W**6*spr1**2 - 12*s*W**6*spr1**2 - &
         & 8*t*W**6*spr1**2 + 6*W**2*(mn**2 + t - W**2)*spr4*(2*(mn**2 - s)*W**2*spr2 + &
         & (mn**2*(t - 2*W**2) - (2*s + t)*(t - W**2))*spr1) + 6*W**2*(mn**2 + t - W**2)*spr3* &
         & (2*W**2*((mn**2 + t - W**2)*spr4 - (s + t - W**2)*spr2) + (mn**2*t - (2*s + t - 2*W**2)*(t - W**2))*spr1) + &
         & t*W**2*(5*mn**4*W**2 + 2*mn**2*(t - W**2)*(2*s + 3*W**2) + (t - W**2)*(-4*s*(s + t) +   &
         & (4*s + 5*t)*W**2 - 5*W**4))*mpi**2) - &
         & 2*C3v*C5v*mn*W*(12*mn**4*s*W**2*spr2**2 - 12*mn**2*s**2*W**2*spr2**2 - 12*s**2*t*W**2*spr2**2 - &
         & 12*s*t**2*W**2*spr2**2 - 24*mn**3*s*W**3*spr2**2 + 24*mn*s**2*W**3*spr2**2 - 12*mn**3*t*W**3*spr2**2 + &
         & 24*mn*s*t*W**3*spr2**2 + 12*mn*t**2*W**3*spr2**2 - 12*mn**4*W**4*spr2**2 + 24*mn**2*s*W**4*spr2**2 - &
         & 12*s**2*W**4*spr2**2 + 12*mn**2*t*W**4*spr2**2 + 24*t**2*W**4*spr2**2 + 24*mn**3*W**5*spr2**2 -   &
         & 24*mn*s*W**5*spr2**2 + 12*mn*t*W**5*spr2**2 - &
         & 12*mn**2*W**6*spr2**2 + 12*s*W**6*spr2**2 - 24*t*W**6*spr2**2 + 6*mn**5*t*W*spr2*spr1 - 24*mn**3*s*t*W*spr2*spr1 + &
         & 24*mn*s**2*t*W*spr2*spr1 - 12*mn**3*t**2*W*spr2*spr1 + 24*mn*s*t**2*W*spr2*spr1 + 6*mn*t**3*W*spr2*spr1 -   &
         & 24*mn**4*s*W**2*spr2*spr1 + &
         & 24*mn**2*s**2*W**2*spr2*spr1 - 12*mn**4*t*W**2*spr2*spr1 + 48*mn**2*s*t*W**2*spr2*spr1 - 24*s**2*t*W**2*spr2*spr1 - &
         & 24*s*t**2*W**2*spr2*spr1 + 12*t**3*W**2*spr2*spr1 + 48*mn**3*s*W**3*spr2*spr1 - 48*mn*s**2*W**3*spr2*spr1 +   &
         & 12*mn**3*t*W**3*spr2*spr1 - &
         & 72*mn*s*t*W**3*spr2*spr1 - 12*mn*t**2*W**3*spr2*spr1 + 24*mn**4*W**4*spr2*spr1 - 48*mn**2*s*W**4*spr2*spr1 +   &
         & 24*s**2*W**4*spr2*spr1 - &
         & 24*mn**2*t*W**4*spr2*spr1 + 48*s*t*W**4*spr2*spr1 - 24*t**2*W**4*spr2*spr1 - 48*mn**3*W**5*spr2*spr1 +   &
         & 48*mn*s*W**5*spr2*spr1 + &
         & 6*mn*t*W**5*spr2*spr1 + 24*mn**2*W**6*spr2*spr1 - 24*s*W**6*spr2*spr1 + 12*t*W**6*spr2*spr1 +   &
         & 4*mn**6*t*spr1**2 - 8*mn**4*s*t*spr1**2 + &
         & 8*mn**2*s**2*t*spr1**2 - 4*mn**4*t**2*spr1**2 + 32*mn**2*s*t**2*spr1**2 - 24*s**2*t**2*spr1**2 - &
         & 4*mn**2*t**3*spr1**2 - 24*s*t**3*spr1**2 + 4*t**4*spr1**2 - 10*mn**5*t*W*spr1**2 + &
         & 32*mn**3*s*t*W*spr1**2 - 32*mn*s**2*t*W*spr1**2 - 32*mn*s*t**2*W*spr1**2 + 10*mn*t**3*W*spr1**2 +   &
         & 12*mn**4*s*W**2*spr1**2 - &
         & 12*mn**2*s**2*W**2*spr1**2 + 12*mn**4*t*W**2*spr1**2 - 56*mn**2*s*t*W**2*spr1**2 +   &
         & 36*s**2*t*W**2*spr1**2 + 60*s*t**2*W**2*spr1**2 + &
         & 4*t**3*W**2*spr1**2 - 24*mn**3*s*W**3*spr1**2 + 24*mn*s**2*W**3*spr1**2 + 56*mn*s*t*W**3*spr1**2 +  &
         &  12*mn*t**2*W**3*spr1**2 - 12*mn**4*W**4*spr1**2 + &
         & 24*mn**2*s*W**4*spr1**2 - 12*s**2*W**4*spr1**2 + 16*mn**2*t*W**4*spr1**2 - 48*s*t*W**4*spr1**2 -  &
         &  20*t**2*W**4*spr1**2 + 24*mn**3*W**5*spr1**2 - &
         & 24*mn*s*W**5*spr1**2 - 22*mn*t*W**5*spr1**2 - 12*mn**2*W**6*spr1**2 + 12*s*W**6*spr1**2 +   &
         & 12*t*W**6*spr1**2 - 6*W*spr4* &
         & (2*W*(mn**6 - mn**4*(s + t) - 2*mn**5*W - mn**2*(t - W**2)**2 + (s + t)*(t - W**2)**2 + &
         & mn**3*W*(2*s + t + 2*W**2) + mn*W*(t**2 - (2*s + t)*W**2))*spr2 - &
         & (2*mn**6*W - 2*mn**4*(s + t)*W + mn**5*(t - 4*W**2) - 2*mn**2*W*(t - W**2)**2 + 2*(s + t)*W*(t - W**2)**2 - &
         & 2*mn**3*(t - 2*W**2)*(s + t + W**2) + mn*(t - W**2)*(t*(2*s + t) + (4*s + 3*t)*W**2))*spr1) - 6*W*spr3* &
         & (2*(-t + (mn - W)**2)*W*(mn**2 - t - W**2)*(mn**2 + t - W**2)*spr4 - &
         & 2*W*(-(mn**3*W*(2*s + t - 2*W**2)) + mn**4*(s - W**2) - (s - W**2)*(t - W**2)**2 + &
         & mn*W*(t**2 + (2*s + t)*W**2 - 2*W**4))*spr2 - &
         & (mn**5*t + 2*W*(s - W**2)*(t - W**2)**2 + mn**4*(-2*s*W + 2*W**3) + &
         & mn*(t - W**2)*(t*(2*s + t) + (4*s - t)*W**2 - 4*W**4) - 2*mn**3*(t*(s + t) - (2*s + t)*W**2 + 2*W**4))*spr1) + &
         & t*W**2*(-mn**6 + 10*mn**5*W + mn**4*(-4*s + t - 3*W**2) - 4*mn**3*(2*s*W + 3*W**3) - &
         & (t - W**2)*(-12*s**2 - 12*s*t + t**2 + 4*(3*s + 2*t)*W**2 - 9*W**4) + &
         & 2*mn*W*(4*s**2 + 4*s*t - 5*t**2 - 4*s*W**2 + 5*W**4) + &
         & mn**2*(4*s**2 - 8*s*t + t**2 + 8*s*W**2 - 6*t*W**2 + 5*W**4))*mpi**2)))/(9.*mn**4*W**4)


  end function LH



  function gamd(w)
    !        Delta -> N pi

    real::gamd
    real,intent(in)::w 
    real::sdel,qcma

    sdel=w**2
    qcma=qcm(sdel)
    if (sdel.le.(mn+mpi)**2) then
       gamd=0.
    else
       !             
       gamd=1./6./pi*(fs/mpi)**2*(sqrt(mn**2+qcma**2)+mn)/2./w*qcma**3*ffs(w)**2

       !              gamd=1./(6.*pi)*(fs/mpi)**2*mn*qcm(sdel)**3/w     ! non relativistic
    end if
  end function gamd

  real function ffs(w)    ! s-channel form factor from M. Post

    real, intent(in)::w

    select case (nfs)
    case (0)
       ffs=1.
    case (1)
       ffs=Lam**4/(Lam**4+(w**2-mdel**2)**2)
    case default
       stop 'ffs'
    end select

  end function ffs

  function qcm(sdel)
    !          Returns the 3-momentum of the pion formed after the decay of a
    !          resonance (R->N pi) of inv. mass s in the rest frame of the resonance

    real::qcm
    real,intent(in)::sdel

    qcm=sqrt((sdel-mpi**2-mn**2)**2-4.*mpi**2*mn**2)/2./sqrt(sdel)

  end function qcm




  subroutine getKinematics_local(e, ef,thf,phipi,thpi,k,p,kf,pf,ppi,q,success)

    real::cosa,sina         ! cos and sin of the angle between q and ppi         
    real::disc              ! discriminant of the eq. for the pion momentum
    real::ppi_pl,ppi_mn     ! the two solutions of the pion momentum 


    logical, intent(out) :: success
    real, intent(out), dimension(0:3) :: k,p,kf,pf,ppi,q
    real, intent(in) :: e, ef,thf,phipi,thpi
    real :: qm,t

    success=.false.

    k(0)=e
    k(1)=0.
    k(2)=0.
    k(3)=e

    p(0)=mn
    p(1)=0.
    p(2)=0.
    p(3)=0.

    kf(0)=ef
    kf(1)=ef*sin(thf)
    kf(2)=0.
    kf(3)=ef*cos(thf)


    q=k-kf
    t=sp4(q,q)
    qm=sqrt(sp(q,q))

    !        angle between q and ppi
    cosa=(q(1)*sin(thpi)*cos(phipi)+q(2)*sin(thpi)*sin(phipi)+q(3)*cos(thpi))/qm
    sina=sqrt(1.-cosa**2)    ! always appears squared, so its sign is not relevant

    !        Determination of the pion momentum
    !        In general, for a given set of kinematical variables,  there are two solutions that have to be added

    disc=(mpi**2-t-2.*mn*q(0))**2-4.*mpi**2*(mn**2+(qm*sina)**2)

    if (disc.lt.0.) then 
       !       write(*,*)'disc<0 -> no solutions'
       return
    endif

    ppi_pl=((2.*mn*q(0)+t+mpi**2)*qm*cosa+(mn+q(0))*sqrt(disc))/2./((mn+q(0))**2-(qm*cosa)**2)
    ppi_mn=((2.*mn*q(0)+t+mpi**2)*qm*cosa-(mn+q(0))*sqrt(disc))/2./((mn+q(0))**2-(qm*cosa)**2)


    if (ppi_pl.ge.0.) then
       ppi(0)=sqrt(mpi**2+ppi_pl**2)
       ppi(1)=ppi_pl*sin(thpi)*cos(phipi)
       ppi(2)=ppi_pl*sin(thpi)*sin(phipi)
       ppi(3)=ppi_pl*cos(thpi)
    else if (ppi_mn.ge.0.) then
       ppi(0)=sqrt(mpi**2+ppi_mn**2)
       ppi(1)=ppi_mn*sin(thpi)*cos(phipi)
       ppi(2)=ppi_mn*sin(thpi)*sin(phipi)
       ppi(3)=ppi_mn*cos(thpi)
    end if
    pf=p+q-ppi

    success=.true.
  end subroutine getKinematics_local


  function sp4(p1,p2)

    real, dimension(0:3), intent(in):: p1,p2
    real sp4
    sp4=p1(0)*p2(0)-p1(1)*p2(1)-p1(2)*p2(2)-p1(3)*p2(3)
  end function sp4


  function sp(p1,p2)
    !        Scalar product of two 3-momenta

    real, dimension(0:3), intent(in):: p1,p2
    real sp

    sp=p1(1)*p2(1)+p1(2)*p2(2)+p1(3)*p2(3)

  end function sp

end module electronNuc_to_Delta
