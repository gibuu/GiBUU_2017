!******************************************************************************
!****m* /VecMesWinkel_Muehlich
! NAME
! module VecMesWinkel_Muehlich
! PURPOSE
! This module provides angular distributions for gamma N -> V N,
! according to Muehlich PhD, Appendix E.
!******************************************************************************
module VecMesWinkel_Muehlich

  use constants, only: mN, mPi
  implicit none

  private

  ! coupling constants
  real,parameter::eqed=0.303

  public :: SelectCost

contains


  real function SelectCost(srts,idv,idb,fn,massv,massb) result(cost)
    ! selects scattering angle cos(t) for gamma N -> V N / V Delta (low energy)
    ! calls: photorho,photomega,photophi,photodel
    use random, only: rn
    use twoBodyTools, only: pcm
    real,    intent(in) :: srts         ! sqrt(s)
    integer, intent(in) :: idv,idb      ! IDs of outgoing vector meson & baryon
    integer, intent(in) :: fn           ! charge of outgoing baryon
    real,    intent(in) :: massv,massb  ! masses of outgoing vector meson & baryon

    real::s,tmin,tmax,t,pini,pout,m2max,m2
    integer::ivec,niter
    logical::flag

    ! settings
    ivec=(idv-101)/2
    if (ivec/=1.and.ivec/=2.and.ivec/=3) then
        write(*,*) 'problem SelectCost ivec',ivec,idv,idb
        stop
    end if
    s=srts**2
    pini=pcm(srts,mn,0.)
    pout=pcm(srts,massv,massb)
    if (pini<=0..or.pout<=0.) then
        cost = 1.-2.*rn()
        return
    end if
    tmin=massv**2-2.*(pini*sqrt(massv**2+pout**2)-pini*pout)
    tmax=massv**2-2.*(pini*sqrt(massv**2+pout**2)+pini*pout)

    ! determine maximum
    if (idb==1) then
!         if(ivec==1)then
!           call photorho(s,tmin,massv,fn,m2max)
        if (ivec==2) then
          m2max = photomega(s,tmin,massv,fn)
        else if (ivec==3) then
          m2max = photophi(s,tmin,massv,fn)
        end if
!     else if(idb==2)then
!         call photodel(s,tmin,massv,massb,ivec,m2max)
    end if

    ! monte carlo
    flag=.true.
    niter=0
    do while(flag)
        niter=niter+1
        t=tmin+rn()*(tmax-tmin)
        if (idb==1) then
!           if(ivec==1)then
!               call photorho(s,t,massv,fn,m2)
          if (ivec==2) then
            m2 = photomega(s,t,massv,fn)
          else if (ivec==3) then
            m2 = photophi(s,t,massv,fn)
          end if
!         else if(idb==2)then
!           call photodel(s,t,massv,massb,ivec,m2)
        end if
        if (m2>=rn()*m2max)flag=.false.
        if (niter>10000.and.flag) then
          flag=.false.
          t=tmin+rn()*(tmax-tmin)
          write(*,*) 'problem SelectCost niter',niter,massv,srts
        end if
    end do

    ! result
    cost=max(min((t-massv**2+2.*pini*sqrt(massv**2+pout**2))/(2.*pini*pout),1.),-1.)

  end function SelectCost


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!   subroutine photorho(s,t,mv,inuc,msq)
!     real,intent(in)::s,t,mv
!     integer,intent(in)::inuc
!     real,intent(out)::msq
!     real::pini,pout,dsdt,dndt,sigm,srt,cost
!     real::errabs=0.,errel=0.01,errest
!     external::dqdag !qdag
!     real(kind=8)::norm
!     real(kind=8),external::dsrho
!     real,parameter::matrix=0.16
!     integer,save::first=0
!     integer,parameter::n=100
!     real,dimension(1:n),save::xn,yn,y2n
!     real,parameter::mrho=0.77
!     ! variables for dsrho
!     real,save::srts
!
!     pini=qcm(sqrt(s),mn,0.)
!     pout=qcm(sqrt(s),mn,mv)
!     if(pout==0.)then
!        msq=0.
!        return
!     end if
!     cost=(t-2*mn**2+2*sqrt(mn**2+pini**2)*sqrt(mn**2+pout**2))/2./pini/pout
!     if(first==0)then
!        call ininorm
!        first=1
!     end if
!     sigm=1./(pini*s)*matrix*pout !mb
!     srts=sqrt(s)              !important: after ininorm!!!
!     if(srts>(mn+mrho+0.02))then
!        if(srts<=xn(1))then
!           norm=yn(1)
!        else if(srts<xn(n))then
!           call splint_nr(xn,yn,y2n,n,srts,norm)
!        else
!           norm=xn(n)
!        end if
!        if(norm/=0.)then
!           dndt=dsrho(cost)/norm*pi/(pini*pout)
!           msq=sigm*dndt*64*pi*s*pini**2*4.*GeVSquared_times_mb
!        else
!           msq=0.
!        end if
!     else
!        dndt=1./(4.*pi)*pi/(pini*pout) !below threshold: isotropic
!        msq=sigm*dndt*64*pi*s*pini**2*4.*GeVSquared_times_mb
!     end if
!
!   contains
!
!     subroutine ininorm
!       ! calculates integral I[ds/do,cost=-1:1]*2pi
!       real,parameter::srtmin=mn+mrho,srtmax=4.
!       integer::i
!       write(*,*)'init of normalization for photorho'
!       do i=1,n
!          srts=srtmin+float(i-1)/float(n-1)*(srtmax-srtmin) !comrho
!          norm=0.
!          if(pout>0.)then
!             call erset(0,-1,0)
!             call dqdag(dsrho,-1.,1.,errabs,errel,1,norm,errest)
!          end if
!          xn(i)=srts
!          yn(i)=norm*2.*pi
!       end do
!       call spline_nr(xn,yn,n,5.e+30,5.e+30,y2n)
!     end subroutine ininorm
!
!   end subroutine photorho
!
!
!   function dsrho(cost)
!     ! dsigma/dOmega [microbarn/sr] on-shell
!     real(kind=8)::dsrho
!     real,intent(in)::cost
!     real::dsdt,pini,pout,t
!     pini=qcm(srts,mn,0.)
!     pout=qcm(srts,mn,mrho)
!     if(pout>0.)then
!        t=2.*mn**2-2.*(sqrt(mn**2+pini**2)*sqrt(mn**2+pout**2)-pini*pout*cost)
!        call dsig(srts,t,dsdt,1)
!        dsrho=dsdt*pini*pout/pi
!     else
!        dsrho=0.
!     end if
!   end function dsrho

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real function photomega(s,t,mv,inuc) result(msq)
    ! calculation of squared matrix element for omega photoproduction
    use constants, only: ii
    real,intent(in)::s,t,mv
    integer,intent(in)::inuc

    complex::m1,m2,m3,m4,m5,m6,m7,m8,m9
    complex::Ds,Du,Dsd,Dud
    real::fft,ffb,ff12p,u,mr,gamma,fn,kn

    ! couplings
    real, parameter :: gpnn=12.85
    real, parameter :: fvpg=1.78*eqed
    real, parameter :: kv=-0.94
    real, parameter :: kneut=-1.91,kprot=1.79
    real, parameter :: gv=3.94
    ! cutoffs
    real, parameter :: lamb=0.96
    real, parameter :: as=0.6
    real, parameter :: lpnn=0.73
    real, parameter :: lvpg=1.
    ! resonance parameters
    real, parameter :: gv2=-1.655
    real, parameter :: g1=1.2
    real, parameter :: mres=1.749
    real, parameter :: gamres=0.445
    real, parameter :: bratio=0.134
    real, parameter :: ag=0.5,av=0.5
    real, parameter :: lg=1.69,lv=2.

    u=2.*mn**2+mv**2-s-t
    mr=mres
    fn=float(inuc)
    if (inuc==0) then
       kn=kneut
    else
       kn=kprot
    end if
    ! width and propagators
    gamma=reswidth(s,mv)
    Ds=1./(s-mr**2+ii*gamma/2.)
    Du=1./(u-mr**2+ii*gamma/2.)
    Dsd=conjg(Ds)
    Dud=conjg(Du)
    ! form factors
    fft=form_tchannel(t)
    ffb=form_born(s,u)
    ff12p=(ag*resform(s,lg)+(1.-ag)*resform(u,lg))*(av*resform(s,lv)+(1.-av)*resform(u,lv))
    ! pion pole diagram
    m1=-((fft**2*fvpg**2*gpnn**2*(mv**2 - t)**2*t)/(mv**2*(mpi**2 - t)**2))
    ! born diagrams
    m2=(eqed**2*ffb**2*gv**2*(-8*fn*kn*(1 + kv)*mn**2*(mn**2 - s)*(mn**2   &
        + mv**2 - s - t)*(-2*mn**2*(mv**2 - t)**2 + kv*(mn**4*(mv**2 +     &
        2*t) + s*(-mv**2 + s + t)*(mv**2 + 2*t) + mn**2*(2*mv**4 -         &
        mv**2*(2*s + t) - t*(4*s + t)))) + kn**2*(mn**2 - s)*(mn**2        &
        + mv**2 - s - t)*(8*mn**2*(mn**2 - s)*t*(-mn**2 - mv**2 + s + t)   &
        + 8*kv*mn**2*(-(mn**4*(mv**2 + 2*t)) - s*(-mv**2 + s + t)*(mv**2   &
        + 2*t) + mn**2*(-2*mv**4 + mv**2*(2*s + t) + t*(4*s + t))) +       &
        kv**2*(2*mn**8 - 4*mn**6*(mv**2 + 2*(s + t)) + s*(-mv**2 + s + t)  &
        *(2*s*(s + t) - mv**2*(2*s + t)) + mn**2*(-mv**6 - 8*s*(s + t)**2  &
        + mv**4*(4*s + t) + 2*mv**2*s*(2*s + 3*t)) + mn**4*(-10*mv**4 +    &
        mv**2*(4*s + 7*t) + 2*(6*s**2 + 10*s*t + t**2)))) + 2*fn**2*mn**2* &
        (16*mn**10 + 4*mn**8*((8 + kv**2)*mv**2 - 16*s - kv**2*t) +        &
        8*mn**6*((7 + kv + kv**2)*mv**4 + 12*s**2 + 2*(2 + kv**2)*s*t +    &
        (3 + kv + kv**2)*t**2 - 2*mv**2*((6 + kv**2)*s - 2*(-2 + kv)*t))   &
        -kv**2*s*(mv**2 - s - t)*(mv**6 - 4*mv**4*s - 4*s*t*(s + t) +      &
        mv**2*(4*s**2 + 8*s*t + t**2)) + mn**4*((24 + 7*kv*(8 + 3*kv))*    &
        mv**6 - 4*mv**4*(2*(6 + kv*(2 + 3*kv))*s + (10 + kv*(18 + 7*kv))*  &
        t) + mv**2*(24*(4 + kv**2)*s**2 + 16*(2 + (-4 + kv)*kv)*s*t +      &
        (24 + kv*(24 + 13*kv))*t**2) - 4*(16*s**3 + 2*(8 + 3*kv**2)*s**2*  &
        t + 2*(2 + kv*(2 + 3*kv))*s*t**2 + (2 + kv*(2 + kv))*t**3)) +      &
        mn**2*(8*kv*s*(-mv**2 + s + t)*(mv**4 + 4*mv**2*t + t**2) +        &
        8*s*(-mv**2 + s + t)*(mv**4 - 2*mv**2*s + 2*s**2 + 2*s*t + t**2)   &
        + kv**2*(3*mv**8 - 5*mv**6*(2*s + t) + 8*s*t*(s + t)*(2*s + t) +   &
        mv**4*(24*s**2 + 8*s*t + 3*t**2) - mv**2*(16*s**3 + 32*s**2*t +    &
        10*s*t**2 + t**3))))))/(2.*mn**4*(mn**2 - s)**2*(mn**2 + mv**2 -   &
        s - t)**2)
    ! resonance contribution
    m3=(eqed**2*ff12p**2*g1**2*gv2**2*(Du*(Dud*(mn**2 + mv**2 - s - t)* &
        (2*mn**6 - 6*mn**3*mr*mv**2 - 2*mv**4*s + 4*mv**2*s**2 - 2*s**3 &
        - mv**4*t + 5*mv**2*s*t - 4*s**2*t + mv**2*t**2 - 2*s*t**2 +    &
        6*mn*mr*mv**2*(-mv**2 + s + t) + mn**4*(mv**2 - 2*(3*s + t)) +  &
        mn**2*(mr**2*(mv**2 - 2*t) + 6*s*(s + t) - mv**2*(5*s + t)) -   &
        mr**2*(-2*t*(s + t) + mv**2*(s + 2*t))) + Dsd*(2*mn**5*mr*mv**2 &
        - mn**6*(mv**2 - 2*t) + 2*mn**3*mr*mv**2*(mv**2 - 2*s - t) -    &
        mr**2*s*(mv**2 - 2*t)*(-mv**2 + s + t) + mn**4*(mv**4 - mr**2*  &
        (mv**2 - 2*t) + mv**2*(2*s - t) - 4*s*t) + mn*mr*mv**2*(mv**4 + &
        2*s**2 + 2*s*t + t**2 - 2*mv**2*(s + t)) + mn**2*(-(s*(mv**2 -  &
        2*t)*(-mv**2 + s + t)) + mr**2*(mv**4 + mv**2*(2*s - t) -       &
        4*s*t))))+ Ds*(Dsd*(mn**2 - s)*(2*mn**6 - 6*mn**3*mr*mv**2 +    &
        6*mn*mr*mv**2*s - mn**4*(mv**2 + 6*s) + mr**2*(mv**4 + 2*s*t -  &
        mv**2*(s + t)) + s*(-2*s*(s + t) + mv**2*(2*s + t)) + mn**2*    &
        (-mv**4 - mv**2*s + mr**2*(mv**2 - 2*t) + 2*s*(3*s + t))) +     &
        Dud*(2*mn**5*mr*mv**2 - mn**6*(mv**2 - 2*t) + 2*mn**3*mr*mv**2* &
        (mv**2 - 2*s - t) - mr**2*s*(mv**2 - 2*t)*(-mv**2 + s + t) +    &
        mn**4*(mv**4 - mr**2*(mv**2 - 2*t) + mv**2*(2*s - t) - 4*s*t) + &
        mn*mr*mv**2*(mv**4 + 2*s**2 + 2*s*t + t**2 - 2*mv**2*(s + t)) + &
        mn**2*(-(s*(mv**2 - 2*t)*(-mv**2 + s + t)) + mr**2*(mv**4 +     &
        mv**2*(2*s - t) - 4*s*t))))))/(4.*mn**4)
    ! born-pion interference
    m4=(eqed*ffb*fft*fvpg*gpnn*gv*(mv**2 - t)*(kn*(mn**4*t + kv*(mn**2 - s)* &
        (mn**2 + mv**2 - s - t)*(mv**2 + t) + s*t*(-mv**2 + s + t) -         &
        mn**2*(mv**4 - 3*mv**2*t + 2*t*(s + t))) + fn*(kv*mn**4*t +          &
        kv*s*t*(-mv**2 + s + t) - mn**2*((2 + kv)*mv**4 - (4 + 3*kv)*        &
        mv**2*t + 2*t*(t + kv*(s + t))))))/(mn*mv*(mn**2 - s)*(mpi**2 -      &
        t)*(mn**2 + mv**2 - s - t))
    ! pion-born interference
    m5=(eqed*ffb*fft*fvpg*gpnn*gv*(mv**2 - t)*(kn*(mn**4*t + kv*(mn**2 - s)* &
        (mn**2 + mv**2 - s - t)*(mv**2 + t) + s*t*(-mv**2 + s + t) -         &
        mn**2*(mv**4 - 3*mv**2*t + 2*t*(s + t))) + fn*(kv*mn**4*t +          &
        kv*s*t*(-mv**2 + s + t) - mn**2*((2 + kv)*mv**4 - (4 + 3*kv)*        &
        mv**2*t + 2*t*(t + kv*(s + t))))))/(mn*mv*(mn**2 - s)*(mpi**2 -      &
        t)*(mn**2 + mv**2 - s - t))
    ! resonance-pion interference
    m6=(eqed*ff12p*fft*fvpg*g1*gpnn*gv2*(Ds*(-mn**2 +s) + Du*(mn**2 + mv**2  &
        - s - t))*(mv**2 - t)*(mn*mv**2 + mr*t))/(2.*mn**2*mv*(mpi**2 - t))
    ! pion-resonance interference
    m7=(eqed*ff12p*fft*fvpg*g1*gpnn*gv2*(Dsd*(-mn**2 + s) + Dud*(mn**2 +     &
        mv**2-s-t))*(mv**2-t)*(mn*mv**2+mr*t))/(2.*mn**2*mv*(mpi**2 - t))
    ! born-resonance interference
    m8=(eqed**2*ff12p*ffb*g1*gv*gv2*(Dud*(2*fn*mn*(4*mn**9 + 2*mn**6*mr*((-2 &
        + kv)*mv**2 - kv*t) - 2*mn**7*(2*(-1 + kv)*mv**2 + 8*s + (2 + kv)*   &
        t) + mn**4*mr*((-8 + 3*kv)*mv**4 + mv**2*((8 - 6*kv)*s + 4*(3 -      &
        2*kv)*t) + 2*t*(3*kv*s - 2*t + 2*kv*t)) + kv*mr*s*(-mv**6 + 2*t*     &
        (s + t)**2 + mv**4*(3*s + 5*t) - 2*mv**2*(s**2 + 4*s*t + 3*t**2))    &
        - mn**2*mr*(mv**4*((-4 + 6*kv)*s + 3*kv*t) + 2*kv*t*(3*s**2 +        &
        4*s*t + t**2) + mv**2*((4 - 6*kv)*s**2 + 4*(1 - 4*kv)*s*t -          &
        5*kv*t**2)) - mn*(mv**2 - s - t)*(4*s*(s + t)*(-mv**2 + s + t) +     &
        kv*(mv**6 + 2*s*t*(s + t) + mv**2*(2*s + t)**2 - mv**4*(5*s +        &
        2*t))) - mn**3*((4 + 7*kv)*mv**6 + 2*s*(s + t)*(8*s + 3*(2 +         &
        kv)*t) + 4*mv**2*((-5 + 3*kv)*s**2 + (-3 + 4*kv)*s*t + (1 + 2*kv)*   &
        t**2) - mv**4*(8*t + 5*kv*(4*s + 3*t))) + mn**5*(4*(6*s**2 +         &
        5*s*t + t**2 - mv**2*(4*s + t)) + kv*(-11*mv**4 + 2*t*(3*s + t)      &
        + 2*mv**2*(6*s + 5*t)))) + kn*(mn**2 + mv**2 - s - t)*(-4*mn*        &
        (mn**4*mr*t + mn**2*mr*(mv**2 - 2*s - t)*t + mn**5*(mv**2 + t) +     &
        mr*s*t*(-mv**2 + s + t) + 2*mn**3*(mv**4 - s*t - mv**2*(s + t)) +    &
        mn*s*(-mv**4 + mv**2*s + t*(s + t))) + kv*(2*mn**8 - 2*mn**5*mr*     &
        (mv**2 + 2*t) - 2*mn*mr*s*(-mv**2 + s + t)*(mv**2 + 2*t) -           &
        2*mn**6*(mv**2 + 4*s + 2*t) + mn**4*(-6*mv**4 + 5*mv**2*t + 12*s*    &
        (s+t)) + 2*mn**3*mr*(-2*mv**4 + mv**2*(2*s + t) + t*(4*s + t)) +     &
        mn**2*(-mv**6 + mv**4*(2*s + t) + 2*mv**2*s*(3*s + 2*t) - 4*s*       &
        (2*s**2 + 3*s*t + t**2)) + s*(2*s*(s + t)**2 + mv**4*(2*s + t) -     &
        mv**2*(4*s**2 + 5*s*t + t**2))))) + Dsd*(-(kn*(mn**2 - s)*(-4*mn*    &
        (mn**4*mr*t + mn**2*mr*(mv**2 - 2*s - t)*t + mn**5*(mv**2 + t) +     &
        mr*s*t*(-mv**2 + s + t) + 2*mn**3*(mv**4 - s*t - mv**2*(s + t)) +    &
        mn*s*(-mv**4 + mv**2*s + t*(s + t))) + kv*(2*mn**8 - 2*mn**5*mr*     &
        (mv**2 + 2*t) - 2*mn*mr*s*(-mv**2 + s + t)*(mv**2 + 2*t) -           &
        2*mn**6*(mv**2 + 4*s + 2*t) + mn**4*(-6*mv**4 + 5*mv**2*t +          &
        12*s*(s + t)) + 2*mn**3*mr*(-2*mv**4 + mv**2*(2*s + t) + t*(4*s      &
        + t)) + mn**2*(-mv**6 + mv**4*(2*s + t) + 2*mv**2*s*(3*s + 2*t) -    &
        4*s*(2*s**2 + 3*s*t + t**2)) + s*(2*s*(s + t)**2 + mv**4*(2*s +      &
        t) - mv**2*(4*s**2 + 5*s*t + t**2))))) + 2*fn*mn*(4*mn**9 -          &
        2*mn**6*mr*((2 + kv)*mv**2 - kv*t) + 2*mn**7*(2*(1 + kv)*mv**2 -     &
        8*s + (2 + kv)*t) + kv*mr*s*(mv**6 - 2*s*t*(s + t) - mv**4*(3*s +    &
        t) + 2*mv**2*s*(s + 2*t)) + mn*s*(2*s*(s + t)*(2*s - kv*t) -         &
        mv**2*(2*s + t)*(2*(2 + kv)*s - kv*t) + mv**4*((4 + 3*kv)*s -        &
        kv*t)) + mn**5*((8 + 5*kv)*mv**4 - 4*mv**2*((4 + 3*kv)*s + (2 +      &
        kv)*t) + 2*s*(12*s - (2 + 3*kv)*t)) + mn**4*mr*(-((8 + 3*kv)*        &
        mv**4) + 2*mv**2*((4 + 3*kv)*s + 2*(3 + kv)*t) - 2*t*(2*t + kv*      &
        (3*s + t))) + mn**3*(4*s*(-2*mv**4 - 4*s**2 - s*t + t**2 + mv**2*    &
        (5*s + t)) + kv*(mv**6 + 4*mv**2*s*(3*s + t) + 2*s*t*(3*s + t) -     &
        mv**4*(8*s + t))) - mn**2*mr*(4*mv**2*s*(-mv**2 + s + t) + kv*       &
        (2*mv**6 - 3*mv**4*(2*s + t) - 2*s*t*(3*s + 2*t) + mv**2*(6*s**2     &
        + 8*s*t + t**2)))))))/(4.*mn**4*(mn**2-s)*(mn**2 + mv**2 - s - t))
    ! resonance-born interference
    m9=(eqed**2*ff12p*ffb*g1*gv*gv2*(Du*(2*fn*mn*(4*mn**9 + 2*mn**6*mr*((-2  &
        + kv)*mv**2 - kv*t) - 2*mn**7*(2*(-1 + kv)*mv**2 + 8*s + (2 + kv)*   &
        t) + mn**4*mr*((-8 + 3*kv)*mv**4 + mv**2*((8 - 6*kv)*s + 4*(3 -      &
        2*kv)*t) + 2*t*(3*kv*s - 2*t + 2*kv*t)) + kv*mr*s*(-mv**6 + 2*t*     &
        (s + t)**2 + mv**4*(3*s + 5*t) - 2*mv**2*(s**2 + 4*s*t + 3*t**2))    &
        - mn**2*mr*(mv**4*((-4 + 6*kv)*s + 3*kv*t) + 2*kv*t*(3*s**2 +        &
        4*s*t + t**2) + mv**2*((4 - 6*kv)*s**2 + 4*(1 - 4*kv)*s*t -          &
        5*kv*t**2)) - mn*(mv**2 - s - t)*(4*s*(s + t)*(-mv**2 + s + t) +     &
        kv*(mv**6 + 2*s*t*(s + t) + mv**2*(2*s + t)**2 - mv**4*(5*s +        &
        2*t))) - mn**3*((4 + 7*kv)*mv**6 + 2*s*(s + t)*(8*s + 3*(2 + kv)*    &
        t) + 4*mv**2*((-5 + 3*kv)*s**2 + (-3 + 4*kv)*s*t + (1 + 2*kv)*       &
        t**2) - mv**4*(8*t + 5*kv*(4*s + 3*t))) + mn**5*(4*(6*s**2 +         &
        5*s*t + t**2 - mv**2*(4*s + t)) + kv*(-11*mv**4 + 2*t*(3*s + t)      &
        + 2*mv**2*(6*s + 5*t)))) + kn*(mn**2 + mv**2 - s - t)*(-4*mn*        &
        (mn**4*mr*t + mn**2*mr*(mv**2 - 2*s - t)*t + mn**5*(mv**2 + t) +     &
        mr*s*t*(-mv**2 + s + t) + 2*mn**3*(mv**4 - s*t - mv**2*(s + t)) +    &
        mn*s*(-mv**4 + mv**2*s + t*(s + t))) + kv*(2*mn**8 - 2*mn**5*mr*     &
        (mv**2 + 2*t) - 2*mn*mr*s*(-mv**2 + s + t)*(mv**2 + 2*t) -           &
        2*mn**6*(mv**2 + 4*s + 2*t) + mn**4*(-6*mv**4 + 5*mv**2*t +          &
        12*s*(s + t)) + 2*mn**3*mr*(-2*mv**4 + mv**2*(2*s + t) + t*(4*s      &
        + t)) + mn**2*(-mv**6 + mv**4*(2*s + t) + 2*mv**2*s*(3*s + 2*t) -    &
        4*s*(2*s**2 + 3*s*t + t**2)) + s*(2*s*(s + t)**2 + mv**4*(2*s +      &
        t) - mv**2*(4*s**2 + 5*s*t + t**2))))) + Ds*(-(kn*(mn**2 - s)*       &
        (-4*mn*(mn**4*mr*t + mn**2*mr*(mv**2 - 2*s - t)*t + mn**5*(mv**2     &
        + t) + mr*s*t*(-mv**2 + s + t) + 2*mn**3*(mv**4 - s*t - mv**2*(s     &
        + t)) + mn*s*(-mv**4 + mv**2*s + t*(s + t))) + kv*(2*mn**8 -         &
        2*mn**5*mr*(mv**2 + 2*t) - 2*mn*mr*s*(-mv**2 + s + t)*(mv**2 +       &
        2*t) - 2*mn**6*(mv**2 + 4*s + 2*t) + mn**4*(-6*mv**4 + 5*mv**2*t     &
        + 12*s*(s + t)) + 2*mn**3*mr*(-2*mv**4 + mv**2*(2*s + t) + t*        &
        (4*s + t)) + mn**2*(-mv**6 + mv**4*(2*s + t) + 2*mv**2*s*(3*s +      &
        2*t) - 4*s*(2*s**2 + 3*s*t + t**2)) + s*(2*s*(s + t)**2 + mv**4*     &
        (2*s + t) - mv**2*(4*s**2 + 5*s*t + t**2))))) + 2*fn*mn*(4*mn**9     &
        - 2*mn**6*mr*((2 + kv)*mv**2 - kv*t) + 2*mn**7*(2*(1 + kv)*mv**2     &
        - 8*s + (2 + kv)*t) + kv*mr*s*(mv**6 - 2*s*t*(s + t) - mv**4*        &
        (3*s + t) + 2*mv**2*s*(s + 2*t)) + mn*s*(2*s*(s + t)*(2*s - kv*t)    &
        - mv**2*(2*s + t)*(2*(2 + kv)*s - kv*t) + mv**4*((4 + 3*kv)*s -      &
        kv*t)) + mn**5*((8 + 5*kv)*mv**4 - 4*mv**2*((4 + 3*kv)*s + (2 +      &
        kv)*t) + 2*s*(12*s - (2 + 3*kv)*t)) + mn**4*mr*(-((8 + 3*kv)*        &
        mv**4) + 2*mv**2*((4 + 3*kv)*s + 2*(3 + kv)*t) - 2*t*(2*t + kv*      &
        (3*s + t))) + mn**3*(4*s*(-2*mv**4 - 4*s**2 - s*t + t**2 + mv**2*    &
        (5*s + t)) + kv*(mv**6 + 4*mv**2*s*(3*s + t) + 2*s*t*(3*s + t) -     &
        mv**4*(8*s + t))) - mn**2*mr*(4*mv**2*s*(-mv**2 + s + t) + kv*       &
        (2*mv**6 - 3*mv**4*(2*s + t) - 2*s*t*(3*s + 2*t) + mv**2*(6*s**2     &
        + 8*s*t+t**2)))))))/(4.*mn**4*(mn**2 - s)*(mn**2 + mv**2 - s - t))

    ! coherent sum
    msq=real(m1+m2+m3+m4+m5+m6+m7+m8+m9)

  contains

    real function reswidth(s,mass)
      ! total width of resonance
      use constants, only: pi
      use twoBodyTools, only: pcm
      real,intent(in)::s,mass
      real::trace,b2,gam1,gam2
      b2=gv2/(2.*mn)
      trace=-2*(mass**2-(mn-Sqrt(s))**2)*(b2**2*(mass**2+2*(mn+Sqrt(s))**2))
      gam1=1./(8.*pi)*pcm(sqrt(s),mn,mass)/s*trace/2.*(resform(s,lv))**2
      gam2=(1.-bratio)*gamres*(pcm(sqrt(s),mn,mpi)/pcm(mres,mn,mpi))**3
      !& *(resform(s,1.)/resform(mres**2,1.))**2
      reswidth=gam1+gam2
    end function reswidth

    real function resform(x,cutoff)
      ! resonance form factor (RNG and RNV vertices)
      real,intent(in)::x,cutoff
      resform=cutoff**4/(cutoff**4+(x-mres**2)**2)
    end function resform

    real function form_tchannel(t)
      ! t-channel form factors
      real,intent(in)::t
      form_tchannel=(lpnn**2-mpi**2)/(lpnn**2-t)*(lvpg**2-mpi**2)/(lvpg**2-t)
    end function form_tchannel

    real function form_born(s,u)
      ! s- and u-channel form factors
      real,intent(in)::s,u
      form_born=as*lamb**4/(lamb**4+(s-mn**2)**2)+(1.-as)*lamb**4/(lamb**4+(u-mn**2)**2)
    end function form_born

  end function photomega

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real function photophi(scms,t,mp,fn) result(msq)
    ! A.I.TITOV, T.S.H.LEE, H.TOKI, O.STRELTSOVA
    ! PHYS.REV.C 60, 035205 (1999)
    ! "STRCTURE OF THE PHI PHOTOPRODUCTION AMPLITUDE AT A FEW GEV"
    use constants, only: pi,ii,GeVSquared_times_mb
    use twoBodyTools, only: pcm

    real,    intent(in) :: scms, t, mp
    integer, intent(in) :: fn

    real :: pin,pout,tmin,dsig,u,alph1,msqr,gpinn,f,f1,kp,kq,pq,fs,fu !,kcm,sth
    complex :: x1,x1d,x3,x4,x5,x3d,x5d,sum,y1,y2,y3,y4,kapan,sum6 !,x4d
    real, parameter :: alph10=1.08,alphb=0.25
    real, parameter :: gphigapi=-0.141,meta=0.547
    real, parameter :: getann=-3.527,gphigaeta=-0.707
    real, parameter :: c1=2.09,s1=1.6,bslope=1.
    real, parameter :: gphinn=-0.24,kapaphi=0.2
    real, parameter :: a1=0.9,a2=0.1
    real, parameter :: lambda=1.87

    pin=pcm(sqrt(scms),mn,0.)
    pout=pcm(sqrt(scms),mn,mp)
    tmin=mp**2-2.*(pin*sqrt(mp**2+pout**2)-pin*pout)
    if (abs(fn)==1) then
      kapan=1.79
      gpinn=-13.26
    else if (fn==0) then
      kapan=-1.91
      gpinn=13.26
    end if
!    sth=(mp+mn)**2
!    kcm=(scms-mn**2)/(2.*sqrt(scms))
    u=2.*mn**2+mp**2-scms-t
    kp=(scms-mn**2)/2.
    kq=(mp**2-t)/2.
    pq=(mn**2+mp**2-u)/2.
    alph1=alph10+alphb*t
    f=(4.*mn**2-2.8*t)/((4.*mn**2-t)*(1.-(t/0.7))**2)*exp(bslope*(t-tmin)/2.)
    f1=1./sqrt(1./(2.*mn**2*mp**2)*(kp*(kp*mp**2+kq**2)+2.*kp*kq*(pq-mp**2)-kq**2*(pq+mn**2))*4.*mn**2)
    x1=c1*f*f1*cexp(-ii*pi/2.*alph1)*((scms-s1)*alphb)**(alph1)
    x1d=conjg(x1)
    x3=-ii/(t-mpi**2)*eqed/mp*gpinn*gphigapi*(0.7**2-mpi**2)/(0.7**2-t)*(0.77**2-mpi**2)/(0.77**2-t)
    x3d=conjg(x3)
    x4=-ii/(t-meta**2)*eqed/mp*getann*gphigaeta*(1.**2-meta**2)/(1.**2-t)*(0.9**2-meta**2)/(0.9**2-t)
!    x4d=conjg(x4)
    x5=x3+x4
    x5d=conjg(x5)
    fs=lambda**4/(lambda**4+(scms-mn**2)**2)
    fu=lambda**4/(lambda**4+(u-mn**2)**2)
    y1=ii*kapan/2./mn
    y2=ii*kapaphi/2./mn
    y3=eqed*gphinn/(scms-mn**2)*(a1*fs+a2*fu)
    y4=eqed*gphinn/(u-mn**2)*(a1*fs+a2*fu)
    x3=x5
    x3d=x5d
    sum6=(2*mn**4*(2*mp**2 - t)*x1*x1d - t*(2*scms**2 + 2*scms*t +     &
         t**2)*x1*x1d +                                                &
         2*mn**2*(mp**4 + 2*scms*t - mp**2*(4*scms + t))*x1*x1d -      &
         mp**6*t*x3*x3d-mp**4*(4*scms*x1*x1d+t*x1*x1d-2*t**2*x3*x3d) + &
         mp**2*(4*scms**2*x1*x1d + 6*scms*t*x1*x1d + 2*t**2*x1*x1d -   &
         t**3*x3*x3d))/mp**2
    sum=sum6+sum2(scms,t,y1,y2,y3,y4,fn)+sum1(scms,t,y1,y2,y3,y4,x1,x5,fn)
    msqr=real(sum)
    dsig=1./(64.*pi*(scms-mn**2)**2)*msqr/GeVSquared_times_mb
    msq=dsig*pin**2*64.*pi*scms*4.*GeVSquared_times_mb

  end function photophi


  complex function sum1(s,t,y1,y2,y3,y4,x1,x3,fn)
      real, intent(in) :: s,t  !,u
      complex, intent(in) :: y1,y2,y3,y4,x1,x3
      integer, intent(in) :: fn

      complex :: x1d,x3d,part1,part2,part3,part4
      real, parameter :: mp = 1.02

      x1d=conjg(x1)
      x3d=conjg(x3)

      part1=(-2*x1*(fn*(mn**6*(y3 + y4) +                               &
            (0,3)*mn**3*mp**2*(mp**2 - t)*y2*(y3 + y4) +                &
            mp**4*(2*s*y3 - 3*s*y4 - t*y4) +                            &
            mp**2*(-(s*t*(y3 - 5*y4)) - s**2*(y3 - 4*y4) + 2*t**2*y4) - &
            (s**2 + s*t + t**2)*(t*y4 + s*(y3 + y4)) -                  &
            (0,1)*mn*mp**2*(mp**2 - t)*y2*                              &
            (mp**2*(y3 - 4*y4) - t*(y3 - 4*y4) + 3*s*(y3 + y4)) +       &
            mn**4*(t*y3 - 3*s*(y3 + y4) + mp**2*(-3*y3 + 2*y4)) +       &
            mn**2*(t**2*y3 + 2*s*t*y4 + mp**4*(-2*y3 + y4) +            &
            3*s**2*(y3 + y4) + mp**2*(4*s*y3 + t*y3 - 6*s*y4 - t*y4))   &
            ) + y1*((0,1)*mn**3*(3*mp**4 - 2*mp**2*t - t**2)*           &
            (y3 + y4) - 2*mn**6*mp**2*y2*(y3 + y4) +                    &
            (0,1)*mn*(3*mp**4 - 2*mp**2*t - t**2)*                      &
            ((mp**2 - t)*y4 - s*(y3 + y4)) +                            &
            mn**4*mp**2*y2*                                             &
            (2*mp**2*y3 + t*(-3*y3 + y4) + 6*s*(y3 + y4)) -             &
            mn**2*mp**2*y2*                                             &
            (mp**4*(y3 - 3*y4) - 2*s*t*(y3 - 3*y4) + t**2*(y3 - y4) +   &
            6*s**2*(y3 + y4) - 2*mp**2*(t*(y3 - 2*y4) + 2*s*y4)) +      &
            mp**2*y2*(mp**4*(2*s + t)*y4 +                              &
            (2*s + t)*(2*s*t*y4 + t**2*y4 + s**2*(y3 + y4)) -           &
            2*mp**2*(3*s*t*y4 + t**2*y4 + s**2*(y3 + 2*y4))))))/mp**2

      part2=2*x3*(y1*((0,1)*mn**3*(mp**4 - t**2)*y2*(y3 - y4) +         &
            mn**4*t*(y3 + y4) -                                         &
            (0,1)*mn*(mp**4 - t**2)*y2*(s*(y3 - y4) + (mp**2 - t)*y4) + &
            t*((mp**2 - t)**2*y4 + 2*s*(-mp**2 + t)*y4 +                &
            s**2*(y3 + y4)) +                                           &
            mn**2*(-2*mp**2*t*y3 + mp**4*(y3 + y4) +                    &
            t*(t*(y3 - y4) - 2*s*(y3 + y4)))) +                         &
            fn*((0,1)*mn*(mp**2-t)**2*(y3+y4)+mn**4*t*y2*(y3 + y4) +    &
            t*y2*((mp**2 - t)**2*y4 + 2*s*(-mp**2 + t)*y4 +             &
            s**2*(y3 + y4)) +                                           &
            mn**2*y2*(-2*mp**2*t*y3 + mp**4*(y3 + y4) +                 &
            t*(t*(y3 - y4) - 2*s*(y3 + y4)))))

      part3=(-2*x1d*(fn*(mn**6*(y3 + y4) +                              &
            (0,3)*mn**3*mp**2*(mp**2 - t)*y2*(y3 + y4) +                &
            mp**4*(2*s*y3 - 3*s*y4 - t*y4) +                            &
            mp**2*(-(s*t*(y3 - 5*y4)) - s**2*(y3 - 4*y4) + 2*t**2*y4) - &
            (s**2 + s*t + t**2)*(t*y4 + s*(y3 + y4)) -                  &
            (0,1)*mn*mp**2*(mp**2 - t)*y2*                              &
            (mp**2*(y3 - 4*y4) - t*(y3 - 4*y4) + 3*s*(y3 + y4)) +       &
            mn**4*(t*y3 - 3*s*(y3 + y4) + mp**2*(-3*y3 + 2*y4)) +       &
            mn**2*(t**2*y3 + 2*s*t*y4 + mp**4*(-2*y3 + y4) +            &
            3*s**2*(y3 + y4) + mp**2*(4*s*y3 + t*y3 - 6*s*y4 - t*y4))   &
            ) + y1*((0,1)*mn**3*(3*mp**4 - 2*mp**2*t - t**2)*           &
            (y3 + y4) - 2*mn**6*mp**2*y2*(y3 + y4) +                    &
            (0,1)*mn*(3*mp**4 - 2*mp**2*t - t**2)*                      &
            ((mp**2 - t)*y4 - s*(y3 + y4)) +                            &
            mn**4*mp**2*y2*                                             &
            (2*mp**2*y3 + t*(-3*y3 + y4) + 6*s*(y3 + y4)) -             &
            mn**2*mp**2*y2*                                             &
            (mp**4*(y3 - 3*y4) - 2*s*t*(y3 - 3*y4) + t**2*(y3 - y4) +   &
            6*s**2*(y3 + y4) - 2*mp**2*(t*(y3 - 2*y4) + 2*s*y4)) +      &
            mp**2*y2*(mp**4*(2*s + t)*y4 +                              &
            (2*s + t)*(2*s*t*y4 + t**2*y4 + s**2*(y3 + y4)) -           &
            2*mp**2*(3*s*t*y4 + t**2*y4 + s**2*(y3 + 2*y4))))))/mp**2

      part4=-2*x3d*(y1*((0,1)*mn**3*(mp**4 - t**2)*y2*(y3 - y4) +       &
            mn**4*t*(y3 + y4) -                                         &
            (0,1)*mn*(mp**4 - t**2)*y2*(s*(y3 - y4) + (mp**2 - t)*y4) + &
            t*((mp**2 - t)**2*y4 + 2*s*(-mp**2 + t)*y4 +                &
            s**2*(y3 + y4)) +                                           &
            mn**2*(-2*mp**2*t*y3 + mp**4*(y3 + y4) +                    &
            t*(t*(y3 - y4) - 2*s*(y3 + y4)))) +                         &
            fn*((0,1)*mn*(mp**2-t)**2*(y3+y4)+mn**4*t*y2*(y3+y4) +      &
            t*y2*((mp**2 - t)**2*y4 + 2*s*(-mp**2 + t)*y4 +             &
            s**2*(y3 + y4)) +                                           &
            mn**2*y2*(-2*mp**2*t*y3 + mp**4*(y3 + y4) +                 &
            t*(t*(y3 - y4) - 2*s*(y3 + y4)))))

      sum1=part1+part2+part3+part4

  end function


  complex function sum2(s,t,y1,y2,y3,y4,fn)
      real, intent(in):: s,t  !,u
      complex, intent(in) :: y1,y3,y2,y4
      integer, intent(in) :: fn

      real, parameter :: mp = 1.02

      sum2=4*(fn*y1*(3*(mn**2 - s)*                                &
           ((0,2)*mn*(mn**2 - mp**2 - s + t) +                     &
           3*(mn**4 + s*(-mp**2 + s) -                             &
           mn**2*(mp**2 + 2*s - 2*t))*y2 -                         &
           (0,1)*mn*(-mp**4 - 2*mp**2*s + mp**2*t - 2*s*t +        &
           2*mn**2*(mp**2 + t))*y2**2)*y3**2 +                     &
           ((0,2)*mn*(6*mn**4 - mp**4 - 6*mp**2*s + 6*s**2 +       &
           6*mn**2*(mp**2 - 2*s - t) + 2*mp**2*t + 6*s*t -         &
           t**2) + (mn**4*(5*mp**2 + t) +                          &
           s*(5*mp**2 + t)*(-mp**2 + s + t) -                      &
           mn**2*(9*mp**4 + 2*mp**2*(5*s - 12*t) +                 &
           t*(2*s + 15*t)))*y2 -                                   &
           (0,1)*mn*(mp**2 - t)*                                   &
           (4*mn**4 - 3*mp**4 - 4*mn**2*(mp**2 + 2*s - t) +        &
           4*s*(s + t) + mp**2*(-4*s + 3*t))*y2**2)*y3*y4 -        &
           3*(mn**2 + mp**2 - s - t)*                              &
           ((0,-2)*mn*(mn**2 + 2*mp**2 - s - 2*t) +                &
           3*(mn**4 - mn**2*(mp**2 + 2*s) -                        &
           (mp**2 - s - t)*(s + t))*y2 +                           &
           (0,1)*mn*(3*mp**4 + 2*mn**2*(mp**2 + t) -               &
           2*t*(s + t) - mp**2*(2*s + t))*y2**2)*y4**2) +          &
           y1**2*((mn**2 - s)**2*                                  &
           (4*mn**2 + 2*t - (0,6)*mn*(mp**2 + t)*y2 +              &
           (2*mn**4 + 2*s*(s + t) - mp**2*(2*s + t) -              &
           2*mn**2*(3*mp**2 + 2*s + t))*y2**2)*y3**2 +             &
           2*mn*(4*mn*(mn**2 - s)*(mn**2 + mp**2 - s - t) -        &
           (0,2)*(mp**2 - t)*                                      &
           (mn**4 + mn**2*(-mp**2 - 2*s + t) +                     &
           s*(-mp**2 + s + t))*y2 +                                &
           mn*(mp**6 + 4*mn**4*t - 2*mp**4*t +                     &
           mp**2*t*(-4*s + t) + 4*s*t*(s + t) +                    &
           4*mn**2*(mp**4 - (mp**2 + 2*s)*t))*y2**2)*y3*y4 +       &
           (mn**2 + mp**2 - s - t)**2*                             &
           (4*mn**2 + 2*t - (0,6)*mn*(mp**2 + t)*y2 +              &
           (2*mn**4 + 2*s*(s + t) - mp**2*(2*s + t) -              &
           2*mn**2*(3*mp**2 + 2*s + t))*y2**2)*y4**2) +            &
           fn**2*(s*(2*(-mp**2 + s + t) +                          &
           (mp**4 + 2*s*t - mp**2*(s + t))*y2**2)*y3**2 +          &
           2*mp**2*(-2*t + (mp**2*(-s + t) + s*(s + t))*y2**2)*y3* &
           y4 + (mp**2 - s - t)*                                   &
           (-2*s + (-2*t*(s + t) + mp**2*(s + 2*t))*y2**2)*y4**2   &
           + (0,6)*mn**5*y2*(y3 + y4)**2 +                         &
           4*mn**6*y2**2*(y3 + y4)**2 +                            &
           mn**4*((6 + (-9*mp**2 + 2*(-4*s + t))*y2**2)*y3**2 -    &
           2*(-8 + (11*mp**2 + 4*(2*s + t))*y2**2)*y3*y4 -         &
           (-14 + (17*mp**2 + 8*s + 6*t)*y2**2)*y4**2) +           &
           mn**2*(-((-6*mp**2 + 2*t +                              &
           (3*mp**4 + mp**2*(6*s - t) + 4*s*(-s + t))*y2**2        &
           )*y3**2) +                                              &
           2*(2*mp**2 - 2*t +                                      &
           (3*mp**4 + 4*s*(s + t) + mp**2*(-6*s + 3*t))*y2**2      &
           )*y3*y4 -                                               &
           (-10*mp**2 + 8*s + 6*t +                                &
           (7*mp**4 - 4*s*(s + t) - mp**2*(2*s + 5*t))*y2**2)      &
           *y4**2) -                                               &
           (0,6)*mn**3*y2*(y3 + y4)*                               &
           (-(t*(y3 - 3*y4)) + 2*s*(y3 + y4) +                     &
           2*mp**2*(2*y3 + y4)) +                                  &
           (0,2)*mn*y2*(3*s**2*(y3 + y4)**2 -                      &
           3*s*(y3 + y4)*(t*(y3 - 3*y4) + 2*mp**2*y4) +            &
           y4*(mp**2*t*(11*y3 - 9*y4) - t**2*(y3 - 6*y4) +         &
           mp**4*(2*y3 + 3*y4)))))

  end function


end module VecMesWinkel_Muehlich
