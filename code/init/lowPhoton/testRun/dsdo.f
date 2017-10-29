      subroutine dsdo(srts,x,dsigdo,dsigdot,finmass,pinitial)
*this subroutine calulates dsigma/dOmega for one-pion photoproduction
*input: srts
*       theta
*       ka: 1  p pi0
*           2  n pi+
*           3  p pi-
*           4  n pi0
*output: dsigdo(ka,0:nres2) = partial cross sections without interferences
*        dsigdot(ka) = total cross sections
      implicit none
      include "common"
      include "cominput"
      integer nsrt,nres2,lmax
      parameter (nsrt=67,nres2=4,lmax=4)
*variables for spline
      real srts2(nsrt),d1(4),d2(4)
      real zwer(nsrt),zwei(nsrt),zwmr(nsrt),zwmi(nsrt)
      real dumer(nsrt),dumei(nsrt),dummr(nsrt),dummi(nsrt)
      real zspler(4,0:2*lmax+1,nsrt),zsplei(4,0:2*lmax+1,nsrt),
     &     zsplmr(4,0:2*lmax+1,nsrt),zsplmi(4,0:2*lmax+1,nsrt)
      complex ample(0:2*lmax+1),amplm(0:2*lmax+1)
      complex emul(4,0:2*lmax+1,nsrt),mmul(4,0:2*lmax+1,nsrt)
      real emulr(nsrt),emuli(nsrt),mmulr(nsrt),mmuli(nsrt)
      real imzw,rezw,fac2
*resonances
      real br1440
      parameter (br1440=0.69)
      real a12(nres2,0:1),a32(nres2,0:1),ruhe(nres2),breite(nres2),
     &     br(nres2),x2(nres2),q0r(nres2),k0r(nres2),jang,alpha,
     &     cleb,aa(nres2,4),ba(nres2,4)
      real ratio(numchr)
      real reself,pci2,pcf2
      integer ang1(nres2),ang2(nres2),iso(nres2),signu(nres2),qn,qpi,
     &     ires(nres2)
      real finmass(2),finmass2(2),resmass,pinitial,g1pi
*amplitudes
      complex ampl(0:lmax,-1:1,4),bmpl(0:lmax,-1:1,4)
      complex hn(4),hsp(4),hsa(4),hd(4)
      complex hnr,hspr,hsar,hdr
      complex amplr,bmplr
      real diff1,diff2

      real gg,gtot

      real dsigdo(4,0:nres2),dsigdot(4)
*
      integer ka,l2,i2,k5,l3,i,kmin,kmax,qnuk,qnuk2,qn2,j,k2
      real qpion,k,srts,srt2,theta,x
      real mn,mpi

      save qpion,k,kmin,kmax
*linux
      save a12,a32,ruhe,breite,br,x2,q0r,k0r,aa,ba,ang1,ang2,iso,
     &     signu,ires,ampl,bmpl,srts2,zspler,zsplei,zsplmr,zsplmi,
     &     emul,mmul
*2pion variables
      integer n2pip1,n2pip2,n2pip3,n2pin1,n2pin2,n2pin3,n2pit
!      parameter(n2pip1=40,n2pip2=57,n2pip3=60,n2pin1=15,n2pin2=20,
!     &     n2pin3=24)

*PM+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
      parameter(n2pip1=68,n2pip2=49,n2pip3=79,n2pin1=68,n2pin2=49,!
     &     n2pin3=24)                                             !
*PM+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
      !änderung

      real plab1(n2pip1),plab2(n2pip2),plab3(n2pip3)
      real plab1n(n2pin1),plab2n(n2pin2),plab3n(n2pin3)
      real s2pipm(n2pip1),s2pip0(n2pip2),s2pi00(n2pip3)
      real s2pipmn(n2pin1),s2pim0n(n2pin2),s2pi00n(n2pin3)
      real t2pipm(n2pip1),t2pip0(n2pip2),t2pi00(n2pip3)
      real t2pipmn(n2pin1),t2pim0n(n2pin2),t2pi00n(n2pin3)
      real deltasig,deltap
      real sig2pi00,sig2pipm,sig2pip0,sigrt,sigres(nres2),sig2pi(0:3),
     &     plab,r2pi1,r2pi2,r2pi3
      real sigvec(3),sigveci(3)
      real sigma(3)
*linux
      save plab1,plab2,plab3,plab1n,plab2n,plab3n,s2pipm,s2pip0,
     &     s2pi00,s2pipmn,s2pim0n,s2pi00n,t2pipm,t2pip0,t2pi00,
     &     t2pipmn,t2pim0n,t2pi00n

      logical :: onlyTwoPi


*********************************************************************


      theta=acos(max(min(x,1.),-1.))
      do ka=1,4
         dsigdot(ka)=0.
         do j=0,nres2
            dsigdo(ka,j)=0.
         end do
      end do

      do ka=kmin,kmax
         hn(ka)=cmplx(0.,0.)
         hsp(ka)=cmplx(0.,0.)
         hsa(ka)=cmplx(0.,0.)
         hd(ka)=cmplx(0.,0.)
         do l2=0,lmax-1
            hn(ka)=hn(ka)+cos(theta/2.)*sqrt(2.)*(ampl(l2,1,ka)-
     &           ampl(l2+1,-1,ka))*
     &           (diff1(x,l2)-diff1(x,l2+1))
            hsa(ka)=hsa(ka)+sin(theta/2.)*sqrt(2.)*
     &           (ampl(l2,1,ka)+ampl(l2+1,-1,ka))*
     &           (diff1(x,l2)+diff1(x,l2+1))
            if(l2.ne.0) then
               hsp(ka)=hsp(ka)+cos(theta/2.)*sin(theta)/
     &              sqrt(2.)*(bmpl(l2,1,ka)-bmpl(l2+1,-1,ka))*
     &              (diff2(x,l2)-diff2(x,l2+1))
               hd(ka)=hd(ka)+sin(theta/2.)*sin(theta)/sqrt(2.)*
     &              (bmpl(l2,1,ka)+bmpl(l2+1,-1,ka))*
     &              (diff2(x,l2)+diff2(x,l2+1))
            end if
         end do
*background
         dsigdo(ka,0)=qpion/2./pinitial*(abs(hn(ka))**2+
     &        abs(hd(ka))**2+
     &        abs(hsp(ka))**2+abs(hsa(ka))**2)

         do i2=1,nres2
            if(iso(i2).eq.1) then
               resmass=finmass(1)
            else
               resmass=finmass(2)
            end if
*            if((i2.eq.1).and.(idelpot.ge.1)) then
*               reself=2.*et*(vdel-vnuk)
*            else
               reself=0.
*            end if
            gg=breite(i2)*(k/k0r(i2))**(2*ang2(i2))*((k0r(i2)**2+
     &           x2(i2)**
     &           2)/(k**2+x2(i2)**2))**ang2(i2)
            if(i2.eq.3) gg=breite(i2)*k/k0r(i2)
            if((i2.ne.1.or.iwsc.eq.0).and.(i2.ne.2.or.imed1520.eq.0))
     &           then
               call manley(srts,ires(i2),0,0,ratio,gtot,0)
*               call manley(resmass,ires(i2),0,0,ratio,gtot,0)
               g1pi=gtot*ratio(1)
            else
               if(i2.eq.1) then
                  call manley(srts,2,1,0,ratio,g1pi,1)
                  call manley(srts,2,0,0,ratio,gtot,1)
*                  call manley(resmass,2,1,0,ratio,g1pi,1)
*                  call manley(resmass,2,0,0,ratio,gtot,1)
               else if(i2.eq.2) then
                  call manley(srts,7,1,0,ratio,g1pi,0)
                  call manley(srts,7,0,0,ratio,gtot,1)
*                  call manley(resmass,7,1,0,ratio,g1pi,0)
*                  call manley(resmass,7,0,0,ratio,gtot,1)
               end if
*               write(401,*)srts,gtot
            end if
            amplr=aa(i2,ka)*(k0r(i2)*q0r(i2)/(k*qpion))**(1./2.)*
     &           resmass*(g1pi*gg)**
     &           (1./2.)*(cmplx(ruhe(i2)**2-resmass**2,resmass*
     &           gtot))/
     &           ((ruhe(i2)**2-resmass**2)**2+
     &           resmass**2*gtot**2)
            bmplr=ba(i2,ka)*(k0r(i2)*q0r(i2)/(k*qpion))**(1./2.)*
     &           resmass*(g1pi*gg)**
     &           (1./2.)*(cmplx(ruhe(i2)**2-resmass**2,resmass*
     &           gtot))/
     &           ((ruhe(i2)**2-resmass**2)**2+
     &           resmass**2*gtot**2)
            l2=ang1(i2)
            hspr=cmplx(0.,0.)
            hdr=cmplx(0.0)
            if(signu(i2).eq.1) then
               hnr=cos(theta/2.)*sqrt(2.)*amplr*
     &              (diff1(x,l2)-diff1(x,l2+1))
               hsar=sin(theta/2.)*sqrt(2.)*amplr*
     &              (diff1(x,l2)+diff1(x,l2+1))
               if(l2.ne.0) then
                  hspr=cos(theta/2.)*sin(theta)/sqrt(2.)*bmplr*
     &                 (diff2(x,l2)-diff2(x,l2+1))
                  hdr=sin(theta/2.)*sin(theta)/sqrt(2.)*bmplr*
     &                 (diff2(x,l2)+diff2(x,l2+1))
               end if
            else if(signu(i2).eq.-1) then
               hnr=-cos(theta/2.)*sqrt(2.)*amplr*
     &              (diff1(x,l2-1)-diff1(x,l2))
               hsar=sin(theta/2.)*sqrt(2.)*amplr*
     &              (diff1(x,l2-1)+diff1(x,l2))
               hspr=-cos(theta/2.)*sin(theta)/sqrt(2.)*bmplr*
     &              (diff2(x,l2-1)-diff2(x,l2))
               hdr=sin(theta/2.)*sin(theta)/sqrt(2.)*bmplr*
     &              (diff2(x,l2-1)+diff2(x,l2))
            else
               write (*,*)signu(i2),i2
               stop
            end if
            dsigdo(ka,i2)=qpion/2./pinitial*(abs(hnr)**2+abs(hsar)**2+
     &           abs(hspr)**2+abs(hdr)**2)

            hn(ka)=hn(ka)+hnr
            hsa(ka)=hsa(ka)+hsar
            hsp(ka)=hsp(ka)+hspr
            hd(ka)=hd(ka)+hdr
*Ende des Loops ueber Resonanzen
         end do
*Ende des Loops ueber Kanaele
      end do
*Aufsummation des 1 pi Anteils
      do ka=kmin,kmax
*dsigdot ist der kohaerent aufaddierte Wirkungsquerschnitt
         dsigdot(ka)=qpion/2./pinitial*(abs(hn(ka))**2+
     &        abs(hsa(ka))**2+
     &        abs(hsp(ka))**2+abs(hd(ka))**2)
      end do
      return
*************************************************************************
      entry gamma2pi(qnuk,srts,sigres,sig2pi,finmass2,onlyTwoPi)
*entry for calculation of 2 pion cross sections
*input: qnuk,srts
*output: sigres(nres2),sig2pi(0:3)
*      write(*,*)'in gamma2pi'
      sig2pi00=0.
      sig2pipm=0.
      sig2pip0=0.
      if(srts.gt.rmass+2.*pmass) then
         resmass=srts
*subtract first vacuum resonance contributions from experimental data
         do i2=1,nres2
            call manley(srts,ires(i2),0,0,ratio,gtot,0)

            gg=(k/k0r(i2))**(2*ang2(i2))*((k0r(i2)**2+x2(i2)**2)/
     &           (k**2+x2(i2)**2))**(ang2(i2))
            if(i2.eq.3) gg=k/k0r(i2)
            sigrt=resmass**2*gtot*gg/((resmass**2-
     &           ruhe(i2)
     &           **2)**2+resmass**2*gtot**2)*(k0r(i2)/k)**2*
     &           2*mn/ruhe(i2)*(a32(i2,qnuk)**2+a12(i2,qnuk)**2)*
     &           389.
*     p pi0 pi0
            sig2pi00=sigrt*(2./9.*ratio(numchrs+1)+
     &           1./3.*ratio(numchrs+3)+1./9.*ratio(numchrs+4)*br1440)
     &           +sig2pi00
*     p pi+ pi-
            sig2pipm=sigrt*(5./9.*ratio(numchrs+1)+
     &           1./3.*ratio(numchrs+2)+4./9.*ratio(numchrs+4)*br1440+
     &           2./3.*ratio(numchrs+3))
     &           +sig2pipm
*     n pi+ pi0
            sig2pip0=sigrt*(2./9.*ratio(numchrs+1)+
     &           2./3.*ratio(numchrs+2)+4./9.*ratio(numchrs+4)*br1440)
     &           +sig2pip0

         end do
         plab=(srts**2-rmass**2)/2./rmass
*     experimental 2 pion cross sections
         if(qnuk.eq.1) then
            if ((plab.ge.plab1(1)).and.(plab.le.plab1(n2pip1)))
     &           then
               call splint2(plab1,s2pipm,t2pipm,n2pip1,plab,r2pi1)
            else if(plab.lt.plab1(1)) then
               r2pi1=0.
            else
               r2pi1=s2pipm(n2pip1)
            end if
            if ((plab.ge.plab2(1)).and.(plab.le.plab2(n2pip2)))
     &           then
               call splint2(plab2,s2pip0,t2pip0,n2pip2,plab,r2pi2)
            else if(plab.lt.plab2(1)) then
               r2pi2=0.
            else
               r2pi2=s2pip0(n2pip2)
            end if
            if ((plab.ge.plab3(1)).and.(plab.le.plab3(n2pip3)))
     &           then
               call splint2(plab3,s2pi00,t2pi00,n2pip3,plab,r2pi3)
            else if(plab.lt.plab3(1)) then
               r2pi3=0.
            else
               r2pi3=s2pi00(n2pip3)
            end if
         else
            if ((plab.ge.plab1n(1)).and.(plab.le.
     &           plab1n(n2pin1))) then
               call splint2(plab1n,s2pipmn,t2pipmn,n2pin1,plab,
     &              r2pi1)
            else if(plab.lt.plab1n(1)) then
               r2pi1=0.
            else
               r2pi1=s2pipmn(n2pin1)
            end if
            if ((plab.ge.plab2n(1)).and.(plab.le.
     &           plab2n(n2pin2))) then
               call splint2(plab2n,s2pim0n,t2pim0n,n2pin2,plab,
     &              r2pi2)
            else if(plab.lt.plab2n(1)) then
               r2pi2=0.
            else
               r2pi2=s2pim0n(n2pin2)
            end if
            if ((plab.ge.plab3n(1)).and.(plab.le.
     &           plab3n(n2pin3))) then
               call splint2(plab3n,s2pi00n,t2pi00n,n2pin3,plab,
     &              r2pi3)
            else if(plab.lt.plab3n(1)) then
               r2pi3=0.
            else
               r2pi3=s2pi00n(n2pin3)
            end if
         end if
* improved threshold behaviour (25.04.03)

         if(plab<=0.45)then
            call thres2pi(srts,sigma,qnuk)
            if(sigma(1)/=0.)r2pi3=sigma(1)
            if(sigma(2)/=0.)r2pi1=sigma(2)
            if(sigma(3)/=0.)r2pi2=sigma(3)
         end if


         If (.not.onlyTwoPi) then
*     subtract explicit contribution of rho-meson
            call vecmes(srts,sigvec,sigveci,0,0.)
            if(iphoch.ne.3) then
               r2pi1=r2pi1-sigvec(1)
               if(sigvec(1)/r2pi1.gt.0.01) Print *,r2pi1,sigvec(1)
            end if

            if(irrespar/=3)then
               r2pi1=max(r2pi1-sig2pipm,0.)
               r2pi2=max(r2pi2-sig2pip0,0.)
               r2pi3=max(r2pi3-sig2pi00,0.)
               if(sig2pipm/r2pi1.gt.0.01) Print *,r2pi1,sig2pipm
               if(sig2pip0/r2pi2.gt.0.01) Print *,r2pi2,sig2pip0
               if(sig2pi00/r2pi3.gt.0.01) Print *,r2pi3,sig2pi00
            end if
         end if
      else
         r2pi1=0.
         r2pi2=0.
         r2pi3=0.
      end if

      sig2pi(0)=r2pi1+r2pi2+r2pi3
      sig2pi(1)=r2pi1
      sig2pi(2)=r2pi2 !piPlus piNull
      sig2pi(3)=r2pi3 !piNull piNull

      if(irrespar==3)then
         sigres=0.
         return
      end if

*calculate (in-medium) resonance contributions
      do i2=1,nres2
         if(iso(i2).eq.1) then
            resmass=finmass2(1)
         else
            resmass=finmass2(2)
         end if
         if((i2.ne.1.or.iwsc.ne.1).and.(i2.ne.2.or.imed1520.eq.0))
     &        then
            call manley(srts,ires(i2),0,0,ratio,gtot,0)
*            call manley(resmass,ires(i2),0,0,ratio,gtot,0)
            g1pi=gtot*ratio(1)
         else
            if(i2.eq.1) then
               call manley(srts,2,1,0,ratio,g1pi,1)
               call manley(srts,2,0,0,ratio,gtot,1)
*               call manley(resmass,2,1,0,ratio,g1pi,1)
*               call manley(resmass,2,0,0,ratio,gtot,1)
            else if(i2.eq.2) then
               call manley(srts,7,1,0,ratio,g1pi,0)
               call manley(srts,7,0,0,ratio,gtot,1)
*               call manley(resmass,7,1,0,ratio,g1pi,0)
*               call manley(resmass,7,0,0,ratio,gtot,1)
            end if
         end if
         gg=(k/k0r(i2))**(2*ang2(i2))*((k0r(i2)**2+x2(i2)**2)/
     &        (k**2+x2(i2)**2))**(ang2(i2))
         if(i2.eq.3) gg=k/k0r(i2)
         sigres(i2)=resmass**2*(gtot-g1pi)*gg/((resmass**2-
     &        ruhe(i2)
     &        **2)**2+resmass**2*gtot**2)*(k0r(i2)/k)**2*
     &        2*mn/ruhe(i2)*(a32(i2,qnuk)**2+a12(i2,qnuk)**2)*
     &        389.

      end do
      return
**********************************************************************
      entry inidsdo(qnuk2,srt2)
*Da Multipole quadratisch in cross-section eingehen, hat sich
*lineare Interpolation nicht bewaehrt
*      write(*,*)'in inidsdo',qnuk2,srt2
      pci2=(srt2**2-rmass**2)**2/(4.*srt2**2)
      k=sqrt(pci2)
      pcf2=((srt2**2-rmass**2-pmass**2)**2-4.*rmass**2*pmass**2)/
     &     (4.*srt2**2)
      qpion=sqrt(max(pcf2,0.))
      kmin=-2*qnuk2+3
      kmax=kmin+1
      do ka=kmin,kmax
         do i2=0,2*lmax+1
            do k5=1,nsrt
               dumer(k5)=real(emul(ka,i2,k5))
               dumei(k5)=aimag(emul(ka,i2,k5))
               dummr(k5)=real(mmul(ka,i2,k5))
               dummi(k5)=aimag(mmul(ka,i2,k5))

               zwer(k5)=zspler(ka,i2,k5)
               zwei(k5)=zsplei(ka,i2,k5)
               zwmr(k5)=zsplmr(ka,i2,k5)
               zwmi(k5)=zsplmi(ka,i2,k5)
            end do

*            write(*,*)'vor splint2'
            if(srt2.lt.srts2(nsrt)) then
               call splint2(srts2,dumer,zwer,nsrt,srt2,rezw)
               call splint2(srts2,dumei,zwei,nsrt,srt2,imzw)
               ample(i2)=cmplx(rezw,imzw)

               call splint2(srts2,dummr,zwmr,nsrt,srt2,rezw)
               call splint2(srts2,dummi,zwmi,nsrt,srt2,imzw)
               amplm(i2)=cmplx(rezw,imzw)
            else
               ample(i2)=cmplx(dumer(nsrt),dumei(nsrt))
               amplm(i2)=cmplx(dummr(nsrt),dummi(nsrt))
            end if
*            write(*,*)'nach splint2'
         end do
*     build amplitudes A, B
         do l2=0,lmax-1
            l3=2*l2+1
            ampl(l2,1,ka)=0.5*(float(l2+2)*ample(l3)+
     &           float(l2)*amplm(l3))
            bmpl(l2,1,ka)=ample(l3)-amplm(l3)
            l3=2*(l2+1)
            ampl(l2+1,-1,ka)=0.5*(float(l2+2)*amplm(l3)-
     &           float(l2)*ample(l3))
            bmpl(l2+1,-1,ka)=ample(l3)+amplm(l3)
         end do
         fac2=1.
         if(ka.eq.2) fac2=-1.
         do l2=0,lmax-1
            ampl(l2,1,ka)=ampl(l2,1,ka)/(0.197*1000.)*fac2
            ampl(l2+1,-1,ka)=ampl(l2+1,-1,ka)/(0.197*1000.)*fac2
            bmpl(l2,1,ka)=bmpl(l2,1,ka)/(0.197*1000.)*fac2
            bmpl(l2+1,-1,ka)=bmpl(l2+1,-1,ka)/(0.197*1000.)*fac2
         end do

*         write(*,*)'vor subtraction of resonance contr'
*     subtract resonance contributions
         do i2=1,nres2

            call manley(srt2,ires(i2),0,0,ratio,gtot,0)

            gg=breite(i2)*(k/k0r(i2))**(2*ang2(i2))*((k0r(i2)**2+
     &           x2(i2)**2)/(k**2+x2(i2)**2))**ang2(i2)
            if(i2.eq.3) gg=breite(i2)*k/k0r(i2)
            ampl(ang1(i2),signu(i2),ka)=ampl(ang1(i2),signu(i2),ka)
     &           -aa(i2,ka)*(k0r(i2)*q0r(i2)/(k*qpion))**(1./2.)*
     &           srt2*(gtot*gg*ratio(1))**
     &           (1./2.)*(cmplx(ruhe(i2)**2-srt2**2,srt2*gtot))/
     &           ((ruhe(i2)**2-srt2**2)**2+srt2**2*gtot**2)

            bmpl(ang1(i2),signu(i2),ka)=bmpl(ang1(i2),signu(i2),ka)
     &           -ba(i2,ka)*(k0r(i2)*q0r(i2)/(k*qpion))**(1./2.)*
     &           srt2*(gtot*ratio(1)*gg)**
     &           (1./2.)*(cmplx(ruhe(i2)**2-srt2**2,srt2*gtot))/
     &           ((ruhe(i2)**2-srt2**2)**2+srt2**2*gtot**2)
         end do
      end do
      return
      entry inidsdo2


      print *, 'In Inidsdo2'
*Einlesen der Helizitaetsamplituden und Resonanzparameter
      a32(1,1)=-0.26
      a12(1,1)=-0.141
      a32(2,1)=0.163
      a12(2,1)=-0.022
      a32(3,1)=0.

*      a12(3,1)=0.068
*      a12(3,1)=0.125

      a12(3,1)=5.82/1000.*(151./.43)**0.5


      a32(4,1)=0.135
      a12(4,1)=-0.014
      a32(1,0)=-0.26
      a12(1,0)=-0.141
      a32(2,0)=-0.137
      a12(2,0)=-0.062
      a32(3,0)=0.

*      a12(3,0)=-0.059
*      a12(3,0)=-0.1

      a12(3,0)=-sqrt(2./3.)*a12(3,1)


      a32(4,0)=-0.035
      a12(4,0)=0.027
      data (ruhe(i),i=1,nres2) / 1.232, 1.524, 1.534, 1.684 /
      data (ang1(i),i=1,nres2) / 1, 2, 0, 3 /
      data (ang2(i),i=1,nres2) / 1, 1, 1, 2 /
      data (breite(i),i=1,nres2) / 0.118, 0.124, 0.151, 0.139 /
      data (signu(i),i=1,nres2) / 1, -1, 1, -1 /
      data (iso(i),i=1,nres2) / 3, 1, 1, 1 /
      data (br(i),i=1,nres2) / 1., 0.59, 0.51, 0.7 /
      data (ires(i),i=1,nres2) / 2, 7, 4, 16 /
      mn=rmass
      mpi=pmass
      do i=1,nres2
        x2(i)=0.3
        q0r(i)=sqrt(((ruhe(i)**2-mn**2-mpi**2)**2-4.*mn**2*mpi**2)
     &       /(4.*ruhe(i)**2))
        k0r(i)=(ruhe(i)**2-mn**2)/(2.*ruhe(i))
*        write(*,*)i,k0r(i)
        jang=float(ang1(i))+float(signu(i))/2.
        alpha=(1./pi*k0r(i)/q0r(i)/(2.*jang
     &   +1.)*mn/ruhe(i)/breite(i))**(1./2.)
        do ka=1,4
          if(ka.eq.1) then
            qn=1
            qpi=0
          else if(ka.eq.2) then
            qn=1
            qpi=1
          else if(ka.eq.3) then
            qn=0
            qpi=-1
          else if(ka.eq.4) then
            qn=0
            qpi=0
          end if
          qn2=qn-qpi
          if(iso(i).eq.3) then
            if(qpi.eq.0) then
              cleb=sqrt(2./3.)
            else
              cleb=sqrt(1./3.)
            end if
          else
            if(qpi.eq.0) then
              if(qn2.eq.0) then
                cleb=sqrt(1./3.)
              else
                cleb=-sqrt(1./3.)
              end if
            else
              cleb=float(qpi)*sqrt(2./3.)
            end if
          end if

          aa(i,ka)=-float(signu(i))*alpha*cleb*a12(i,qn)
          if(jang.eq.0.5) then
             ba(i,ka)=0.
          else
             ba(i,ka)=float(signu(i))*4*alpha*((2.*jang-1)*
     &            (2.*jang+3.))**(-.5)*
     &            cleb*a32(i,qn)
          end if
       end do
      end do

      !write(*,*)'vor open'
      open(13,file='../buuinput/photo_mult.dat',
     &     status='unknown')
      !write(*,*)'nach open'
      do i=1,nsrt
         read(13,*)srts2(i)
         do j=0,2*lmax+1
            read(13,*)(d1(k2),k2=1,4)
            read(13,*)(d2(k2),k2=1,4)
            do k2=1,4
               emul(k2,j,i)=cmplx(d1(k2),d2(k2))
            end do
            read(13,*)(d1(k2),k2=1,4)
            read(13,*)(d2(k2),k2=1,4)
            do k2=1,4
               mmul(k2,j,i)=cmplx(d1(k2),d2(k2))
            end do
         end do
      end do
      !write(*,*)'nach do open'
      close(13)
*Spline-Init
      do j=0,2*lmax+1
         do k2=1,4
            do i=1,nsrt
               emulr(i)=real(emul(k2,j,i))
               emuli(i)=aimag(emul(k2,j,i))
               mmulr(i)=real(mmul(k2,j,i))
               mmuli(i)=aimag(mmul(k2,j,i))
            end do
            call spline2(srts2,emulr,nsrt,zwer)
            call spline2(srts2,emuli,nsrt,zwei)
            call spline2(srts2,mmulr,nsrt,zwmr)
            call spline2(srts2,mmuli,nsrt,zwmi)
            do i=1,nsrt
               zspler(k2,j,i)=zwer(i)
               zsplei(k2,j,i)=zwei(i)
               zsplmr(k2,j,i)=zwmr(i)
               zsplmi(k2,j,i)=zwmi(i)
            end do
         end do
      end do
      !print *, 'nach spline'
*2 pi off the proton
c      open(13,file='../buuinput/ppipm.dat',status='unknown')
!      open(13,file=pg_ppmrr.dat',status='unknown')
      open(13,file='../buuinput/twoPi/gamp-ppipm.dat',status='unknown')
      read(13,*)n2pit
       !print *, 'hallo welt '
      if(n2pit.ne.n2pip1) then
         write(*,*)'Fehler in gamini',n2pip1,n2pit
         stop
      end if
       !print *, 'hallo welt 1'

      do i=1,n2pip1
!        read(13,*)plab1(i),deltap,s2pipm(i),deltasig
        read(13,*)plab1(i),s2pipm(i)
      end do
      close(13)
       !print *, 'hallo welt 2'
c      open(13,file='../buuinput/ppip0.dat',status='unknown')
!      open(13,file='../buuinput/ppip0_neu.dat',status='unknown')
      open(13,file='../buuinput/twoPi/gamp-npip0.dat',status='unknown')
      read(13,*)n2pit
      !print *, 'hallo welt 3'

      if(n2pit.ne.n2pip2) then
         write(*,*)'Fehler in gamini',n2pip2,n2pit
         stop
      end if

      do i=1,n2pip2
c        read(13,*)plab2(i),deltap,s2pip0(i),deltasig
         read(13,*)plab2(i),s2pip0(i)
      end do
      close(13)
      !Print *, 'Hallo 1'
c      open(13,file='../buuinput/ppi00.dat',status='unknown')
c      open(13,file='../buuinput/ppi00_neu.dat',status='unknown')
!      open(13,file='../buuinput/kottu-pi0pi0-xsec.dat',status='unknown')
      open(13,file='../buuinput/twoPi/gamp-ppi00.dat',status='unknown')
      read(13,*)n2pit
      if(n2pit.ne.n2pip3) then
         write(*,*)'Fehler in gamini',n2pip3,n2pit
         stop
      end if
      do i=1,n2pip3
c        read(13,*)plab3(i),s2pi00(i),deltasig
         read(13,*)plab3(i),s2pi00(i)
      end do
      close(13)
      call spline3(plab1,s2pipm,n2pip1,t2pipm)
      call spline3(plab2,s2pip0,n2pip2,t2pip0)
      call spline3(plab3,s2pi00,n2pip3,t2pi00)

*2 pion off the neutron
!      open(13,file='../buuinput/npipm.dat',status='unknown')
      open(13,file='../buuinput/twoPi/gamn-npipm.dat',status='unknown')
      read(13,*)n2pit
      if(n2pit.ne.n2pin1) then
         write(*,*)'Fehler in gamini',n2pin1,n2pit
         stop
      end if
      do i=1,n2pin1
!        read(13,*)plab1n(i),s2pipmn(i),deltasig
        read(13,*)plab1n(i),s2pipmn(i)
      end do
      close(13)
!      open(13,file='../buuinput/npim0.dat',status='unknown')
      open(13,file='../buuinput/twoPi/gamn-ppim0.dat',status='unknown')
      read(13,*)n2pit
      if(n2pit.ne.n2pin2) then
         write(*,*)'Fehler in gamini',n2pin2,n2pit
         stop
      end if
      !Print *, 'Hallo 2'
      do i=1,n2pin2
!        read(13,*)plab2n(i),s2pim0n(i),deltasig
        read(13,*)plab2n(i),s2pim0n(i)
      end do
      close(13)
c      open(13,file='../buuinput/npi00.dat',status='unknown')
!      open(13,file='../buuinput/npi00_neu.dat',status='unknown')
c      open(13,file='../buuinput/npi00_neu_min.dat',status='unknown')
c      open(13,file='../buuinput/npi00_neu_max.dat',status='unknown')
      open(13,file='../buuinput/twoPi/gamn-npi00.dat',status='unknown')
      read(13,*)n2pit
      if(n2pit.ne.n2pin3) then
         write(*,*)'Fehler in gamini',n2pin3,n2pit
         stop
      end if
      !Print *, 'Hallo 3'
      do i=1,n2pin3
        read(13,*)plab3n(i),s2pi00n(i)
      end do
      close(13)
      call spline3(plab1n,s2pipmn,n2pin1,t2pipmn)
      call spline3(plab2n,s2pim0n,n2pin2,t2pim0n)
      call spline3(plab3n,s2pi00n,n2pin3,t2pi00n)
      print *, 'Am Ende von inidsdo2'
      return
      end
*********************************************************************
      subroutine thres2pi(srts,sigma,qnuk)

*     this subroutine provides a reasonable threshold behaviour for
*     the two pion photoproduction cross sections

      implicit none
      real srts,sigma(3),ener,mn,mp,flux,ps,matrix2(6),thres(6)
      parameter(mn=0.938,mp=0.138)
      integer ich,i,qnuk
      data matrix2 /34.5, 800.0, 150.0, 16.0, 2550.0, 298.0/
      data thres /0.37, 0.42, 0.407, 0.415, 0.373, 0.436/

      sigma=0.0
      ener=(srts**2-mn**2)/2./mn
      flux=2.*mn*ener
      call bops3(ps,srts,mn,mp,mp)
      do i=1,3
         ich=i+3*(1-qnuk)
         if(ener<=thres(ich))then
            sigma(i)=ps/flux*matrix2(ich)
         end if
      end do

      return
      end
*********************************************************************
      real function pleg(x,l)
      implicit none
      integer l
      real x
      if(l.lt.0) then
        pleg=0.
      else if(l.eq.0) then
        pleg=1.
      else if(l.eq.1) then
        pleg=x
      else if(l.eq.2) then
        pleg=(3.*x**2-1.)/2.
      else if(l.eq.3) then
        pleg=(5.*x**3-3.*x)/2.
      else if(l.eq.4) then
        pleg=(35.*x**4-30.*x**2+3.)/8.
      else if(l.eq.5) then
        pleg=(63.*x**5-70.*x**3+15.*x)/8.
      else if(l.ge.6) then
        write(*,*)'error l too large',l
        stop
      end if
      return
      end
**************************************************
      real function diff1(x,l)
      integer l
      real x
      if(l.lt.1) then
        diff1=0.
      else if(l.eq.1) then
        diff1=1.
      else if(l.eq.2) then
        diff1=3.*x
      else if(l.eq.3) then
        diff1=(15.*x**2-3.)/2.
      else if(l.eq.4) then
        diff1=(140.*x**3-60.*x)/8.
      else if(l.eq.5) then
        diff1=(315.*x**4-210.*x**2+15.)/8.
      else if(l.ge.6) then
        write(*,*)'error l too large',l
        stop
      end if
      return
      end
*********************************************
      real function diff2(x,l)
      integer l
      real x
      if(l.lt.2) then
        diff2=0.
      else if(l.eq.2) then
        diff2=3.
      else if(l.eq.3) then
        diff2=15.*x
      else if(l.eq.4) then
        diff2=(420.*x**2-60.)/8.
      else if(l.eq.5) then
        diff2=(1260.*x**3-420.*x)/8.
      else if(l.ge.6) then
        write(*,*)'error l too large',l
        stop
      end if
      return
      end
*********************************************************************
      subroutine vecmes(srts,sig,sigi,imed,rho)
      implicit none
      include "common"
      include "cominput"
      real srts,sig(3),matrix(3),pinitial,ps(4),pfinal2,maxcontr,
     &     sigi(3),matrix2(3),ps2(2)
      integer i,ires,imed
      real rho,spot(3)

      if(imed.eq.1) then
         if(ipotvec.ge.1) then
            call vecspot(rho,0.,spot(1),3)
            call vecspot(rho,0.,spot(2),5)
            call vecspot(rho,0.,spot(3),7)
         end if
      else
         do i=1,3
            spot(i)=0.
         end do
      end if
      matrix(1)=0.16
      pfinal2=(srts**2-(rmass+mesprop1(5,1)+spot(2))**2)*(srts**2-
     &     (rmass-mesprop1(5,1)+spot(2))**2)/4./srts**2
      if(pfinal2.lt.0) pfinal2=0
      matrix(2)=0.08*pfinal2/(2.*(srts-1.73)**2+pfinal2)

      matrix(3)=0.004

      matrix2(1)=0.5
      matrix2(2)=0.5
      matrix2(3)=0.

      pinitial=(srts**2-rmass**2)/2./srts

      do i=1,3
         if(i.eq.1) then
            ires=3+idbmax
         else if(i.eq.2) then
            ires=5+idbmax
         else if(i.eq.3) then
            ires=7+idbmax
         end if
         call massint2(ps,ires,srts,rmass,maxcontr,spot(i))
         sig(i)=1000.*matrix(i)*ps(1)/pinitial/srts**2
         call massint(ps2,ires,srts,rmass,pmass,spot(i))
         sigi(i)=1000.*matrix2(i)*ps2(1)/pinitial/srts

      end do


      return
      end
***************************************************************
      subroutine vecmesa(srts,idv,massv,cost)
      implicit none
      include "common"
      include "cominput"
      integer mbi
      parameter (mbi=5)
      real bi(mbi,2),pini,pfinal2,tmax,tmin,c1,x0,x,t,rn,
     &     pfinal,srts,massv,cost,egamma
      logical flag
      integer i,j,idv
      data ((bi(i,j),j=1,2),i=1,mbi)
*       e_gamma  slope
     &     /1.8, 5.75,
     &     2.5, 5.43,
     &     3.5, 6.92,
     &     4.5, 8.1,
     &     5.8, 7.9/

      egamma=(srts**2-rmass**2)/2./rmass
      flag=.true.
      i=0
      do while(flag)
         i=i+1
         if(i.eq.mbi.or.egamma.lt.bi(i,1)) flag=.false.
      end do

      pini=(srts**2-rmass**2)/2./srts
      pfinal2=(srts**2-(rmass+massv)**2)*(srts**2-(rmass-massv)**2)/
     &     4./srts**2
      if(pfinal2.le.0) then
         write(*,*)'problems in vecmesa',srts,rmass,massv
         stop
      end if
      pfinal=sqrt(pfinal2)

      tmax=2.*rmass**2-2.*(sqrt(rmass**2+pini**2)*sqrt(rmass**2+
     &     pfinal2)-pini*pfinal)
      tmin=2.*rmass**2-2.*(sqrt(rmass**2+pini**2)*sqrt(rmass**2+
     &     pfinal2)+pini*pfinal)


      if(tmax.le.tmin) then
         t=tmin
      else
         c1=bi(i,2)/(exp(bi(i,2)*tmax)-exp(bi(i,2)*tmin))
         x0=-c1/bi(i,2)*exp(bi(i,2)*tmin)

         x=rn(iseed)

         t=log(bi(i,2)/c1*(x-x0))/bi(i,2)
      end if

      cost=(t-2*rmass**2+2*sqrt(rmass**2+pini**2)*sqrt(rmass**2+
     &     pfinal2))/2./pini/pfinal
      if(abs(cost).gt.1) then
         write(*,*)'problems in vecmesa cost',cost
         cost=sign(1.,cost)
      end if
      return
      end
*******************************************************
      subroutine vecmesm(srts,minmass,mass1,mass2,id1,mass,rhor,pr,
     &     spot)
      implicit none
      include "common"
      include "cominput"
      include "comwidth"
      real srts,minmass,mass1,mass2,mass,massmax2,psam,x,xx,gamtot,rn,
     &     spectral,psa,ratio(nmesch2+nmesch3),bwmes,rhor,pr,spot
      integer id1,idm
      logical flag

      idm=id1-idbmax
      if(mesprop1(idm,2).lt.1e-03) then
*
         mass=mesprop1(idm,1)
      else
         massmax2=srts-mass1-mass2-spot
         if(massmax2.le.minmass) then
            write(*,*)'problems in vecmesm masses',
     &           minmass,massmax2
            stop
         end if

         call bops3(psam,srts,mass1,mass2,minmass)

         flag=.true.
         do while(flag)
            x=rn(iseed)
            xx=rn(iseed)
            mass=minmass+(massmax2-minmass)*x

            if((iwsc2.ge.1.and.idm.eq.3).or.
     &           (iwsc2.ge.2.and.idm.eq.5)) then
               rhores=rhor
               pres=pr
               gamtot=bwmes((mass+spot)**2,idm,0,0,0,ratio,1)
            else
               gamtot=bwmes((mass+spot)**2,idm,0,0,0,ratio,0)
            end if
            spectral=mass**2*gamtot/((mass**2-
     &           mesprop1(idm,1)**2)**2+gamtot**2*mass**2)

            call bops3(psa,srts,mass1,mass2,mass+spot)

            if(spectral*psa.gt.xx*psam/mesprop1(idm,2)) then
               flag=.false.
            end if
         end do
      end if
      return
      end
