*
*     Relativistic Nuclear Mean Field (qhd-2)
*     (Relativistic Thomas-Fermi Approximation)
*     Meson Fields (in the order as read in):
*     Vector: Omega Meson (Isoscalar)
*     Scalar: Sigma Meson (Isoscalar, with 3rd and 4th order
*                          Self Interaction)
*     Vector: Rho   Meson (Isovector)
*     EleMag: Photon      (Static Coulomb Interaction)
*     H. Lenske, Munich, January, 1989 (Version 1.0)
*     Revised: Dec.  9, 1993 (Version 3.0)   > Coulomb Interaction included              <
*     Revised: Dec. 20, 1993 (Version 3.1)   > Energy density, B.E. corrected            <
*     Revised: Jun. 17, 1997 (Version 3.2)   > Adjusted to use on PC                     <
*     Revised: Jul. 04, 2001 (Version 4.0)   > Infinite ASYMMETRIC Matter Eos added      <
*                                            > Order of Fields changed to scalar,vector  <
*                                            > Scalar-Isovector delta Meson included     <
*     Revised: March 02, 2002 (ERTF1_0)      > 2nd Order DGL for number density          <
*                                            >   -->  Extended RTF <--                   <
*     Revised: March 04, 2002 (ERTF1_1)      > optional 1st deriv. term added to DGL     <
*                                                                                        <
*     Revised: June  03, 2003 (ERTF1_2)      > KOUT removed from the input list          <
*
ctg      Program ExRTF
      subroutine ExRTF
c      implicit double precision (a-h,o-z)
      Parameter (KRX=400,KBX=6)
cTG      character*255 wkdir
      character*125 prnout,dataf,MTYPE(KBX)
      CHARACTER*4 NAT1,NAT
      common/cntrl/kout(10),epst
      common/pi4pi/pi,fpi,pisq
      common/boson/bmass(KBX),Jbo(KBX),Ibo(KBX),ISV(4)
      common/couple/G4pi(KBX),bcpl(KBX),rnge(KBX),IBX
      common/bnonl/NLINT(KBX),NLPOW(4,KBX),ANL(4,KBX),XNL(4,KBX)
      common/const/Xmss(4),Xmass,alfa,hbc,hbcm,hbcm3
      common/radius/rmax,dr,rst(KRX),rwgt(KRX),nrx
      common/SVF/bfld(KRX,KBX)
      common/VectorD/rho_v(KRX,4),xkf(KRX,2)
      common/ScalarD/rho_s(KRX,4),xstar(KRX,2)
      common/Nucl/ Amass,AN,AZ,b(2),rd
      dimension chem(2)
c
      DIMENSION NAT1(0:3),NAT(111)
      Data Xmss(1),Xmss(2)/938.27231d0,939.56563d0/
      DATA NAT1,NAT/'  N ','  P ','  D ','  T ',
     *        ' H','He','Li','Be',' B',' C',' N',' O',' F','Ne','Na',
     *        'Mg','Al','Si',' P',' S','Cl','Ar',' K','Ca','Sc','Ti',
     *        ' V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As',
     *        'Se','Br','Kr','Rb','Sr',' Y','Zr','Nb','Mo','Tc','Ru',
     *        'Rh','Pd','Ag','Cd','In','Sn','Sb','Te',' I','Xe','Cs',
     *        'Ba','La','Ce','Pr','Nd','Pr','Sm','Eu','Gd','Tb','Dy',
     *        'Ho','Er','Tm','Yb','Lu','Hf','Ta',' W','Re','Os','Ir',
     *        'Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra',
     *        'Ac','Th','Pa',' U','Np','Pu','Am','Cm','Bk','Cf','Es',
     *        'Fm',
     *        'Md','Nb','Lr','Rf','Ha','Sg','Ns','Hs','Mt','X ','X '/
      nat1(0)=nat1(0)
c
      pi =4*atan(1.d0)
      fpi=4*pi
      pisq=pi*pi
      Xmass=0.5*(Xmss(1)+Xmss(2))
c
      write(*,1001)
 1001 format(10x,40('*')/10x,'*',3x,'Relativistic Nuclear Mean-Field',
     1 4x,'*'/10x,'*',8x,'Thomas-Fermi-Calculation',6x,'*'/
     2 10x,40('*'))
!TG      it5=5
      it5=80
      write(*,102)
  102 format(5x,'filename for input data:')
ctg      read(*,*,end=1010)dataf
ctg 1010 open (5,file=dataf,status='OLD')
cTG      open (5,file='ertf.input',status='OLD')
      open (80,file='ertf.input')
      rewind(80)
      dataf = 'ertf.input'

c
cTG      tstart=timer()
c
cTG      call SHOWTIME
c
c Collect results on E, <r^2> etc from consequetive runs:
c
      open(29,file='earm.d',status='UNKNOWN')
c
      write(29,2901)
 2901 format('#',3x,'Mass',7x,'Z',10x,'E(A)',8x,'E(A)/A',
     & 6x,'ChemP(p)',6x,'ChemP(n)',
     & 6x,'<r^2>(p)',6x,'<r^2>(n)',6x,'<r^2>(A)')
c
cTG    1 call INPUT(wkdir,prnout,it5)
    1 call INPUT(prnout,it5)
c
c      tinit=timer()
c
      NN=AN
      NZ=AZ
      NA=Amass
      IZ=min(111,NZ)
      write(8,600)NA,NAT(IZ),NZ,NN,dataf,prnout
  600 format(//,1x,70('*')/
     1 ' *',
     & 10x,'Extended Relativistic Thomas-Fermi Calculation',12x,'*'/
     & ' *',15x,'( Program ERTF, rel.1.2, Jume  2003 )',16x,'*'/
     & ' *',15x,'(       H. Lenske, U. Giessen       )',16x,'*'/
     2 ' *',15x,'(     EoS of Asymmetric Matter      )',16x,'*'/
     2 ' *',15x,'( sigma,delta,omega and rho Mesons  )',16x,'*'/
     2 ' *',15x,'(   Static Coulomb Interaction      )',16x,'*'/
     2 ' *',15x,'(   RTF Densities and Mean-Field    )',16x,'*'/
     2 ' *',15x,'( 2nd order DGL for Number Density  )',16x,'*'/
     2 ' *',15x,'                                     ',16x,'*'/
     & ' *',22X,i3,'-',a2,'(Z=',i3,',N=',i3,')        ',19x,'*'/
     2 ' *',15x,'                                     ',16x,'*'/
cTG     3 ' *', 9x,' Working Directory : ',a38,'*'/
     4 ' *', 9x,' Input  Dataset    : ',a38,'*'/
     3 ' *', 9x,' Output Dataset    : ',a38,'*'/
     4 1x,70('*'))
      write(8,619)
  619 format(////1x,30('*'),'  Mesons  ',30('*'))

      close(it5) !TG close unit 80!!!!!!!!!!!!
C
C CONVENTION: 1. SCALAR-ISOSCLAR  (sigma)  /2. SCALAR-ISOVECTOR (delta; optional)
c             3. ISOSCALAR-VECTOR (omega)  /4. ISOVECTOR-VECTOR (rho  )
C             5. COULOMB (photon)
c
c          ===> Ordering according to this scheme is done in sub INPUT <===
C
      Mtype(1)='Scalar-IsoScalar'
      Mtype(2)='Scalar-IsoVector'
      Mtype(3)='Vector-IsoScalar'
      Mtype(4)='Vector-IsoVector'
      Mtype(5)='Vector-Photon   '
      write(8,6000)
 6000 format(/
     & ' No.',5x,'Mass',4x,'Range','  g^2/4pi',6x,'C^2',2x,'J',2x,'I',
     & 8x,'Type'/
     & '    ',4x,'[MeV]',5x,'[fm]')
      do I=1,IBX
      Rnge(i)=0.d0
      If(Bmass(i).gt.0.d0)then
      bcpl(i)=fpi*G4pi(i)*(Xmass/bmass(i))**2
      Rnge(i)=hbc/Bmass(i)
      endif
      If(JBo(i).eq.0.and.IBo(i).eq. 0)KB=1
      If(JBo(i).eq.0.and.IBo(i).eq. 1)KB=2
      If(JBo(i).eq.1.and.IBo(i).eq. 0)KB=3
      If(JBo(i).eq.1.and.IBo(i).eq. 1)KB=4
      If(JBo(i).eq.1.and.IBo(i).eq.-1)KB=5
      write(8,6001)I,Bmass(i),Rnge(i),G4pi(i),Bcpl(i),JBo(i),IBo(i),
     & Mtype(KB)
      enddo
 6001 format(i4,4f9.4,2i3,3x,a16)
c
      do I=1,IBX
      If(JBo(i).eq.0.and.IBo(i).eq. 0)KB=1
      If(JBo(i).eq.0.and.IBo(i).eq. 1)KB=2
      If(JBo(i).eq.1.and.IBo(i).eq. 0)KB=3
      If(JBo(i).eq.1.and.IBo(i).eq. 1)KB=4
      If(JBo(i).eq.1.and.IBo(i).eq.-1)KB=5
      If(NLINT(i).ne.0)then
      write(8,6010)Mtype(kb)
      do n=1,NLINT(i)
      XNL(n,i)=ANL(n,i)*Xmass/NLPOW(n,i)
      write(8,6011)n,NLPOW(n,i),ANL(n,i),XNL(n,i)
      enddo
      endif
      enddo
 6010 format(//'  Non-Linear Self-Interactions in ',a16,' Channel:'/
     & 10x,' Term',' Power','  Coefficient','  Transformed'/
     & 10x,'     ','      ','             ','     [MeV]   ')
 6011 format(10x,i5,i6,5e13.5)
      write(8,6013)
 6013 format(72('*'))
c
c EoS symmetric Matter:
c
      Xsi=0.5
      open(20,file='eos-50.d',status='UNKNOWN')
      open(21,file='mst-50.d',status='UNKNOWN')
      open(22,file='kin-50.d',status='UNKNOWN')
      open(23,file='mfd-50.d',status='UNKNOWN')
c
      call EoS(Xsi)
c
c EoS asymmetric Matter (Pb-Case: Z/A=82/208):
c
      Xsi=82.d0/208.d0
      open(20,file='eos-39.d',status='UNKNOWN')
      open(21,file='mst-39.d',status='UNKNOWN')
      open(22,file='kin-39.d',status='UNKNOWN')
      open(23,file='mfd-39.d',status='UNKNOWN')
c
      call EoS(Xsi)
c
c EoS asymmetric Matter (Z/A=0.25):
c
      Xsi=0.25
      open(20,file='eos-25.d',status='UNKNOWN')
      open(21,file='mst-25.d',status='UNKNOWN')
      open(22,file='kin-25.d',status='UNKNOWN')
      open(23,file='mfd-25.d',status='UNKNOWN')
c
      call EoS(Xsi)
c
c EoS asymmetric Matter (Z/A=0.10):
c
      Xsi=0.10
      open(20,file='eos-10.d',status='UNKNOWN')
      open(21,file='mst-10.d',status='UNKNOWN')
      open(22,file='kin-10.d',status='UNKNOWN')
      open(23,file='mfd-10.d',status='UNKNOWN')
c
      call EoS(Xsi)
c
c EoS Pure Neutron Matter:
c
      Xsi=0.0
      open(20,file='eos-00.d',status='UNKNOWN')
      open(21,file='mst-00.d',status='UNKNOWN')
      open(22,file='kin-00.d',status='UNKNOWN')
      open(23,file='mfd-00.d',status='UNKNOWN')
c
      call EoS(Xsi)
c
      CALL init
c
      CALL greenf
c
      CALL selfcon(itera,chem)
c
      call Energy(Itera,Chem)
c
cTG      TimeA=timer()-tinit
cTG      TimeT=timer()-tstart
cTG      write(8,699)TimeA,TimeT
cTG  699 format(/' ==> Time for this run  :',f8.2,' sec'/
cTG     #        ' ==> Time ellapsed      :',f8.2,' sec')
c
cTG      close(6)
c
cTG      call SHOWTIME
c
cTG      goto 1   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
      close(8)
      close(20)
      close(21)
      close(22)
      close(23)
      close(29)
c
cTG      stop
c
      return !TG
      end
c---------------------------------------------------------------------
      subroutine DFS(Amass,Z)
c
c Proton/Neutron Density Parameters from Skyrme Systematics
c Fermi Shape are used with half-density radius and diffuseness:
c
c R_{\tau}=R_{0\tau}*A^(1/3)+R_{1\tau}+R_{2\tau}*AS
c a_{\tau}=a_{0\tau}        +a_{1\tau}*AS
c
c A=Mass Number, AS=(N-Z)/A Asymmetry, \tau=1,2 for p,n
c
c Convention: Protons = 1, Neutrons = 2, Average = 3
c rhod contains the central densities (from normalization to Z/N numbers)
c RMS    contains the mass   rms-radii (obtained analytically)
c RCHRG  contains the charge rms-radii (obtained analytically)
c
c      implicit double precision (a-h,o-z)
      COMMON/DensP/rd(3),ad(3),rhod(3),RMS(3),RCHRG(3)
      dimension rnp(3,2),anp(2,2),xmsq(2),rmsd(3),rmsc(3),zpn(2)
      data rnp/1.2490,-0.5401,-0.9582,
     &         1.2131,-0.4415, 0.8931/
      data anp/0.4899,-0.1236, 0.4686,0.0741/
      data xmsq/0.7429d0,-0.113d0/
c
      pi=atan(1.d0)*4.d0
      fpi=4.d0*pi
      A1=Amass
      an=a1-z
      MM=A1
      NZ=Z
      NN=MM-NZ
      rd(3)=0.d0
      ad(3)=0.d0
      anz=z
      dnz=(an-z)/a1
      a3=a1**(1./3.)
* Proton/Neutron Parameters
      do 10 i=1,2
      rd(i) =a3*rnp(1,i)+rnp(2,i)+rnp(3,i)*dnz
      ad(i) =anp(1,i)+anp(2,i)*dnz
      x=(pi*ad(i)/rd(i))**2
c Vol.Integral of a Fermi distrib. (analytical)
      rhod(i)=0.75*anz/(rd(I)**3*pi*(1.+x))
c <r^2>        of a Fermi distrib. (analytical)
      rms(i)=0.2d0*rd(i)**5*(1.d0+x/(0.3d0)*(1.+0.7d0*x))*
     &       fpi*rhod(i)/anz
c
c charge radius (xmsq are the nucleonic <r^2> of p/n)
c
      if(i.eq.1)then
      rch=rms(i)+xmsq(i)
      else
      rch=xmsq(i)
      endif
      RCHRG(i)=rch
c
c Prepare averaged R_0 and diffuseness:
c
      rd(3)=rd(3)+anz*rd(i)
      ad(3)=ad(3)+anz*ad(i)
      RCHRG(3)=RCHRG(3)+anz*rch
      anz=an
  10  continue
c
c Averaged R_0 and diffuseness (do averaging):
c
      RCHRG(3)=RCHRG(3)/Z
      rd(3)=rd(3)/a1
      ad(3)=ad(3)/a1
      rhod(3)=0.75*a1/(rd(3)**3*(1.+(pi*ad(3)/rd(3))**2)*pi)
      i=3
      x=(pi*ad(i)/rd(i))**2
      rms(i)=0.2d0*rd(i)**5*(1.d0+x/(0.3d0)*(1.+0.7d0*x))*
     &       fpi*rhod(i)/a1
      rms(3)=(Z*rms(i)+AN*rms(2))/Amass
      drms  =rms(2)-rms(1)
      crms  =RCHRG(2)-RCHRG(1)
c
c CM - CORRECTION FOR <r**2>
c
      cma=1.d0-1.d0/a1
      zpn(2)=-        z/a1**2
      zpn(1)= (Amass-z)/a1**2
      do i=1,2
      rmsd(i)=cma*rms(i)  +zpn(i)*drms
      rmsc(i)=cma*RCHRG(i)+zpn(i)*crms
      enddo
c
      rmsd(3)=(Z*rmsd(1)+(Amass-Z)*rmsd(2))/Amass
      rmsc(3)=(Z*rmsc(1)+(Amass-Z)*rmsc(2))/Z
c
      do 50 i=1,3
      If(RCHRG(i).gt.0.d0)RCHRG(i)=sqrt(RCHRG(i))
      If(Rmsc (i).gt.0.d0)Rmsc (i)=sqrt( Rmsc(i))
      rmsd(i)=sqrt(rmsd(i))
   50 rms(i)=sqrt(rms(i))
c
      WRITE(8,6005)NZ,NN,MM,
     & (rd(I),i=1,3),(rd(i)/a3,i=1,3),(ad(I),i=1,3),
     $ (RMS(i),i=1,3),(RCHRG(i),i=1,3),
     $ (RMSD(i),i=1,3),(RMSC(i),i=1,3),
     & (rhod(I),I=1,3)
 6005 FORMAT(/28X,'Density Parameters from Systematics:'/
     & 28x,'   Protons','  Neutrons',3x,'Average'/
     &10x,'   Number        :',3I10/
     &10x,'   R             :',3f10.4,'  [fm]'/
     &10x,'   R(reduced)    :',3f10.4,'  [fm]'/
     &10x,'   a             :',3f10.4,'  [fm]'/
     &10x,'   <r^2>         :',3f10.4,'  [fm]'/
     &10x,'   <r^2>(charge) :',3f10.4,'  [fm] <== (Neutrons: [fm^2])'/
     &10x,'   <r^2>         :',3f10.4,'  [fm] <== (recoil corrected)'/
     &10x,'   <r^2>(charge) :',3f10.4,'  [fm] <== (recoil corrected)'/
     &10x,'   rho(central)  :',3f10.4,'  [A/fm^3]')
      return
      end
c----------------------------------------------------------------------
      SUBROUTINE Energy(Itera,Chem)
c      implicit double precision (a-h,o-z)
      Parameter (KBX=6,KRX=400)
c
c  Energies and Energy Density
C  Note: Total Divergences of Meson Fields are integrated out for Etot, Edens
c
      common/cntrl/kout(10),epst
      common/pi4pi/pi,fpi,pisq
      common/boson/bmass(KBX),Jbo(KBX),Ibo(KBX),ISV(4)
      common/couple/G4pi(KBX),bcpl(KBX),rnge(KBX),IBX
      common/bnonl/NLINT(KBX),NLPOW(4,KBX),ANL(4,KBX),XNL(4,KBX)
      common/const/Xmss(4),Xmass,alfa,hbc,hbcm,hbcm3
      common/radius/rmax,dr,rst(KRX),rwgt(KRX),nrx
      common/SVF/bfld(KRX,KBX)
      common/VectorD/rho_v(KRX,4),xkf(KRX,2)
      common/ScalarD/rho_s(KRX,4),xstar(KRX,2)
      common/Nucl/ Amass,AN,AZ,b(2),rd
      dimension Edens(KRX),Eplot(KRX,2),drv(6),chem(2),Efield(KBX)
      dimension Bmss(KRX,KBX),XV(10),RMS_v(3)
c
      I=1
c
      open (1,file='fields.d',status='UNKNOWN')
c
      Do 40 k=1,KBX
       Efield(k)=0.d0
   40 continue
c
      do n=1,KRX
      do k=1,KBX
       Bmss(n,k)=1.d0
      enddo
      enddo
c
      Etot=0.d0
      TEbary=0.d0
      TEmes=0.d0
      TEint=0.0
      do 100 n=1,nrx
      Emes=0.d0
      do 60 k=1,IBX
       if(bcpl(k).eq.0.d0)goto 60
       if(bcpl(k).ne.0.d0)bnorm=1.d0/bcpl(k)
       if(rnge(k).le.0.d0)rinv=0.0
       if(rnge(k).gt.0.d0)rinv=1.d0/rnge(k)
       if(n.eq.1.or.n.eq.nrx)then
        if(n.eq.1  )drv(k)=(bfld(n+1,k)-bfld(n  ,k))/dr
        if(n.eq.nrx)drv(k)=(bfld(n  ,k)-bfld(n-1,k))/dr
       else
        drv(k)=0.5d0*(bfld(n+1,k)-bfld(n-1,k))/dr
       endif
       if(k.le.2)then
        Emes=Emes+(drv(k)**2+(rinv*bfld(n,k))**2)*rnge(k)**2*bnorm
c
c  nonlinear scalar self-interaction
c   (Factor >2< - see Normal. of Emes where this FACTOR IS DIVIDED!)
c   (NOT TO BE INCLUDED IN W'(phi)!)
c
        If(NLINT(k).ne.0)then
        w=0.d0
       DO 10 j=1,NLINT(k)
        npow=NLPOW(j,k)
        Bmss(n,k)=Bmss(n,k)+bcpl(k)*ANL(j,k)*bfld(n,k)**(npow-2)
   10   w=w+ANL(j,k)*bfld(n,k)**npow/npow
        Emes=Emes-w*2.d0
        endif
       else
        Emes=Emes-(drv(k)**2+(rinv*bfld(n,k))**2)*bnorm*rnge(k)**2
       endif
        If(k.le.2)Efield(k)=Efield(k)+rwgt(n)*bfld(n,k)*rho_s(n,2+k)
        If(k.gt.2)then
         If(k.lt.IBX)Efield(k)=Efield(k)+rwgt(n)*bfld(n,k)*rho_v(n,k)
         If(k.eq.IBX)Efield(k)=Efield(k)+rwgt(n)*bfld(n,k)*rho_v(n,1)
        endif
   60 continue
c
      Emes=Emes*xmass*0.5d0/hbcm3
c
      Spn=xstar(n,1)*frhos(xkf(n,1),xstar(n,1))+
     & xstar(n,2)*frhos(xkf(n,2),xstar(n,2))
      Etp=xkf(n,1)**3*sqrt(xkf(n,1)**2+xstar(n,1)**2)
      Etn=xkf(n,2)**3*sqrt(xkf(n,2)**2+xstar(n,2)**2)
      Ebary=0.25d0*(Spn+(Etp+Etn)/pisq)*xmass/hbcm3
      Eint=bfld(n,1)*rho_s(n,3)+bfld(n,2)*rho_s(n,4)+
     & bfld(n,3)*rho_v(n,3)+bfld(n,4)*rho_v(n,4)+
     & bfld(n,5)*rho_v(n,1)
c
c Add sigma Self-Interaction:
c
      Eint=Eint-Ephi(Bfld(n,1),1)
c
      Eint=Eint*xmass/hbcm3
      Edens(n)=Ebary+Eint*0.5d0
c
      Etot=Etot+rwgt(n)*Edens(n)
      TEbary=TEbary+rwgt(n)*Ebary
      TEmes=TEmes+rwgt(n)*Emes
      TEint=TEint+rwgt(n)*Eint
c
      Eplot(n,1)=0.d0
      Eplot(n,2)=0.d0
      if(rho_v(n,3).gt.1.d-06)then
c
c Nuclear Energy Density per particle as a function of radius:
c
      Eplot(n,1)=(Edens(n)*hbcm3-Xmass*rho_v(n,3))/rho_v(n,3)
c
c Subtract Coulomb --> Nuclear EoS
c
      Eplot(n,2)=(Edens(n)*hbcm3
     & -xmass*rho_v(n,3)-0.5*bfld(n,IBX)*rho_v(n,1)*xmass)/rho_v(n,3)
      else
      Eplot(n,1)=0.d0
      Eplot(n,2)=0.d0
      endif
c
  100 continue
c
      Etot=Etot*fpi
c
      Do 110 k=1,IBX
       Efield(k)=Efield(k)*fpi*Xmass/(Amass*hbcm3)
  110 continue
      TEbary=TEbary*fpi/Amass-Xmass
      TEmes=TEmes*fpi/Amass
      TEint=TEint*fpi/Amass
      chem(1)=(-1.d0+chem(1))*Xmass
      chem(2)=(-1.d0+chem(2))*Xmass
      Ebnd   =Etot-Xmass*Amass
      EoA    =Ebnd/Amass
c
      write(8,7011)EoA,Ebnd,TEbary
 7011 format(//'  >>>> Binding Energy, Total Energy',
     & ' and Chemical Potentials:  <<<<'/
     & 7x,'Binding Energy (per Baryon)  -  E/A:',f10.5,' [MeV]'/
     & 7x,'Total Binding Energy               :',f10.2,' [MeV]'/
     & 7x,'Baryon  Energy (per Baryon)        :',f10.3,' [MeV]')
C
      ChemA=(B(1)*Chem(1)+B(2)*Chem(2))/Amass
      write(8,7012)chem,ChemA
 7012 format(
     & 7x,'Chemical Potential(Protons )       :',f10.5,' [MeV]'/
     & 7x,'Chemical Potential(Neutrons)       :',f10.5,' [MeV]'/
     & 7x,'Chemical Potential(average )       :',f10.5,' [MeV]'/)
c
c <r^2>:
c
      do i=1,3
      RMS_v(i)=0.d0
      enddo
c
      do n=1,NRX
      rsq=rst(n)*rst(n)
      do i=1,3
      RMS_v(i)=RMS_v(i)+rsq*Rwgt(n)*rho_v(n,i)
      enddo
      enddo
c
      RMS_v(1)=fpi*RMS_v(1)/AZ/hbcm3
      RMS_v(2)=fpi*RMS_v(2)/AN/hbcm3
      RMS_v(3)=fpi*RMS_v(3)/AMass/hbcm3
c
      do i=1,3
      RMS_v(i)=sqrt(RMS_v(i))
      enddo
C
      write(29,2900)AMass,AZ,Ebnd,EoA,Chem,RMS_v
 2900 format(2f8.2,10e14.6)
C
       write(8,7013)(Efield(k),k=1,IBX),TEint
 7013 format('  >>>> Baryon-Meson Interaction (per Baryon) <<<<'/
     & 7x,'sigma  :',f10.5,' [MeV]'/
     & 7x,'delta  :',f10.5,' [MeV]'/
     & 7x,'omega  :',f10.5,' [MeV]'/
     & 7x,'rho    :',f10.5,' [MeV]'/
     & 7x,'Coulomb:',f10.5,' [MeV]'/
     & 7x,24('-')/
     & 7x,'Eint   :',f10.5,' [MeV]'/)
C
 1100 format(f10.2,',',i3)
      Ifo=0
      If(Ifo.gt.0)write(8,700)Amass,itera
  700 format(///10x,'Meson Fields  -  Mass=',f8.1,4x,
     1 'No. of Iterations:',i4/
     2 7x,'r',7x,'sigma',7x,'delta',7x,'omega',8x,'Rho ',5x,'Coulomb',
     & 5x,'M*(p)/M',5x,'M*(n)/M')
c
      write(1,1700)Amass,itera
 1700 format('#',9x,'Meson Fields  -  Mass=',f8.1,4x,
     1 'No. of Iterations:',i4/'#',
     2 6x,'r',9x,'sigma',9x,'delta',9x,'omega',10x,'Rho ',7x,'Coulomb',
     & 7x,'M*(p)/M',7x,'M*(n)/M')
c
      open(2,file='umf-nr.d',status='unknown')
c
      write(2,2010)
 2010 format('#',9x,'Schroedinger Potentials:'/
     & '#',6x,'r',11x,'U0',11x,'U1',6x,'r*Vso_0',6x,'r*Vso_1')
c
      DO 1200 n=1,nrx
      m=n-1
      a1=bfld(n,1)*xmass
      a2=bfld(n,2)*xmass
      a3=bfld(n,3)*xmass
      a4=bfld(n,4)*xmass
      a5=bfld(n,5)*xmass
c Schroedinger Potential (central ;   isoscalar and isovector)
      Usc=a3-a1
      Uvc=a4-a2
c Schroedinger Potential (r*U(s.o.) ; isoscalar and isovector)
      mm=max(1,n-1)
      mp=min(nrx,mm+2)
      mm=max(1,mp-2)
      Um =bfld(mm,1)+bfld(mm,3)
      Up =bfld(mp,1)+bfld(mp,3)
      Usos=hbc*0.5*(Up-Um)/(2.*dr)
      Um =bfld(mm,2)+bfld(mm,4)
      Up =bfld(mp,2)+bfld(mp,4)
      Usov=hbc*0.5*(Up-Um)/(2.*dr)
c
      write(1,1101)rst(n),a1,a2,a3,a4,a5,xstar(n,1),xstar(n,2)
      write(2,1101)rst(n),Usc,Uvc,Usos,Usov
c
      If(Ifo.le.0)goto 1200
      if(n.eq.1.or.m.eq.5*(m/5))then
      write(8,701)rst(n),a1,a2,a3,a4,a5,xstar(n,1),xstar(n,2)
      endif
c
 1200 CONTINUE
      close(1)
      close(2)
 1101 format(f8.3,10e14.6)
  701 format(f8.2,10e12.4)
C
      Ido=0
      if(Ido.ne.0)then
      write(8,750)Amass
  750 format(///10x,'DENSITIES and ENERGY DENSITY - Mass=',f8.1//
     2 7x,'R',7x,'RHOB0',5x,'Edens/A',6x,'PROTON',5x,
     & 'NEUTRON',7x,'RHOB1',8x,'RHOS')
      endif
      Nzero=2*nrx
c
* Densities: rho_s,rho_B,rho_p,rho_n
c
      open(2,file='gsd.d',status='unknown')
c
      do m=1,4
      XV(m)=-0.5d0*rho_v(nrx,m)*rst(nrx)**2
      enddo
c
      write(2,20000)
20000 format('# RTF Densities'/
     a    '#', 9x,'r',9x,'rhov_p',9x,'rhov_n',
     a     9x,'rh0v_0',9x,'rh0v_1',9x,
     a     'rhos_p',
     a     9x,'rhos_n',9x,'rhos_0',9x,'rhos_1')
c
      DO 1250 n=1,nrx
      m=n-1
      dvp=rho_v(n,1)/hbcm3
      dvn=rho_v(n,2)/hbcm3
      dv0=rho_v(n,3)/hbcm3
      dv1=rho_v(n,4)/hbcm3
      dsp=rho_s(n,1)/hbcm3
      dsn=rho_s(n,2)/hbcm3
      ds0=rho_s(n,3)/hbcm3
      ds1=(rho_s(n,2)-rho_s(n,1))/hbcm3
      XMp=xstar(n,1)
      XMn=xstar(n,2)
c
      write(2,3000)rst(n),dvp,dvn,dv0,dv1,dsp,dsn,ds0,ds1
c
      rsq=rst(n)*rst(n)
      do j=1,4
      XV(j)=XV(j)+rsq*rho_v(n,j)
      enddo
      If(dv0.eq.0.d0.and.Nzero.gt.nrx)Nzero=n
      If(n.gt.Nzero+5)goto 1250
c
      If(Ido.le.0)goto 1250
      if(n.eq.1.or.m.eq.5*(m/5))then
      write(8,701)rst(n),dv0,Eplot(n,1),dvp,dvn,dv1,ds0
      endif
c
 1250 continue
 3000 format(f10.5,10e15.7)
      close(2)
      do m=1,5
      XV(m)=XV(m)*dr*fpi/hbcm3
      enddo
      write(8,6900)(XV(m),m=1,4)
 6900 format(/,7x,'Density Volume Integrals:'/
     &         7x,'Vector-Proton   :',f8.2/
     &         7x,'Vector-Neutron  :',f8.2/
     &         7x,'Vector-Isoscalar:',f8.2/
     &         7x,'Vector-Isovector:',f8.2)
c
* Energy density
c
      open(2,file='Edens.d',status='unknown')
      write(2,1255)
 1255 format('#',8x,'r',11x,'rho',11x,'E/A',4x,'w.o. Coul.')
      DO 1260 n=1,nrx
      a1=rho_v(n,3)/hbcm3
      a2=Eplot(n,1)
      a3=Eplot(n,2)
      if(a1.gt.0.d0)
     & write(2,1101)rst(n),a1,a2,a3
 1260 continue
      close(2)
* ms*(r)
      open(2,file='bmass.d',status='unknown')
      write(2,1265)
 1265 format('#',8x,'r',10x,'dens',9x,'sigma',9x,'delta',11x,'omega',
     & 11x,'rho  ')
      do n=1,nrx
      a1=rho_v(n,3)/hbcm3
      write(2,1101)rst(n),a1,(sqrt(Bmss(n,k)),k=1,4)
      enddo
      close(2)
c
      return
      end
c-----------------------------------------------------------------------
      REAL FUNCTION frhos(xkf,xstar)
c      implicit double precision (a-h,o-z)
      common/pi4pi/pi,fpi,pisq
      xkf2=xkf*xkf
      xm2 =xstar*xstar
      ef=sqrt(xkf2+xm2)
      If(xkf.lt.0.d0.or.xstar.le.0.d0)then
       write(8,629)xkf,xstar,ef
  629 format(' xkf,xstar:',2e13.5,' ef:',e13.5)
      stop
      endif
      frhos=0.5*xstar*(xkf*ef-xm2*log((xkf+ef)/xstar))/pisq
      return
      end
c---------------------------------------------------------------------
      SUBROUTINE SelfCon(iter,chem)
c      implicit double precision (a-h,o-z)
      Parameter (Ipr=0,KBX=6,KRX=400,ITERX=100)
      common/cntrl/kout(10),eps
      common/pi4pi/pi,fpi,pisq
      common/SVF/bfld(KRX,KBX)
      common/VectorD/rho_v(KRX,4),xkf(KRX,2)
      common/ScalarD/rho_s(KRX,4),xstar(KRX,2)
      common/boson/bmass(KBX),Jbo(KBX),Ibo(KBX),ISV(4)
      common/couple/G4pi(KBX),bcpl(KBX),rnge(KBX),ibx
      common/bnonl/NLINT(KBX),NLPOW(4,KBX),ANL(4,KBX),XNL(4,KBX)
      common/Nucl/ Amass,AN,AZ,b(2),rd
      common/const/Xmss(4),Xmass,alfa,hbc,hbcm,hbcm3
      common/radius/rmax,dr,rst(KRX),wst(KRX),nrx
      common/MRatio/RM(2)
      dimension Dnew(KRX,4),Snew(KRX,4),Bnew(KRX,KBX),Vnew(KRX)
      Dimension tts(ITERX),perfb(ITERX),perfs(ITERX)
      Dimension rms(3),xms(3),chem(2)
      dimension rms_pn(ITERX,2),Vlam(ITERX,2),Xold(KRX,2),Bold(KRX,KBX)
      dimension Upot(KRX,2)
c
c MVch =0 --> conventinal RTF:
c MVch =1 --> extended    RTF:
c
      MVch=0
c
      RM(1)=Xmss(1)/Xmass
      RM(2)=Xmss(2)/Xmass
c
      write(8,199)
  199 format(///)
c
c  start of the iteration part
c
cTG      tinit=timer()
      Istor=0
      open(7,file='selfc.d',status='unknown')
      iter=0
    1 iter=iter+1
      If(Istor.le.0)rewind 7
      if(iter.gt.iterx)goto 9000
cTG      tt1=timer()
      DO 10 n=1,nrx
      Bnew(n,5)=Bfld(n,5)
      do 10 k=1,4
      Bnew(n,k)=Bfld(n,k)
      Dnew(n,k)=rho_v(n,k)
      Snew(n,k)=rho_s(n,k)
   10 continue
c
c  Self-consistent calculation of the chemical potential,
c  the baryon densities and Fermi-momenta
c  AW and BW are iteration weights
c
      AW=0.4
      BW=1.d0-AW
      PWR=1.d0/3
c
      If(MVch.eq.0)then
      phas=-1.
      do 15 k=1,2
      phas=-phas
      cfac=0.5*(1.+phas)
      do 12 n=1,nrx
   12 Vnew(n)=Bnew(n,3)+phas*Bnew(n,4)+cfac*Bnew(n,IBX)
C
      CALL VCTF(nrx,wst,Vnew,xstar(1,k),Dnew(1,k),B(k),xkf(1,k),chem(k))
C
      Vlam(iter,k)=chem(k)
      do n=1,nrx
      Dnew(n,k)=AW*rho_v(n,k)+BW*Dnew(n,k)
      Xkf(n,k) =(3*pisq*Dnew(n,k))**pwr
      enddo
   15 continue
      endif
c
      If(MVch.eq.1)then
c
c      call VchD(Dnew,XKf,Bnew,Chem,Upot)
c
      do k=1,2
      Vlam(iter,k)=chem(k)
      do n=1,nrx
      Dnew(n,k)=AW*rho_v(n,k)+BW*Dnew(n,k)
      Xkf(n,k) =(3*pisq*Dnew(n,k))**pwr
      enddo
      enddo
      endif
c
c Update IsoScalar and IsoVector Vector Densities:
c
      do 18 n=1,nrx
      Dnew(n,3)=Dnew(n,1)+Dnew(n,2)
   18 Dnew(n,4)=Dnew(n,1)-Dnew(n,2)
c
c Update vector fields:
c
      do 20 k=3,4
c
      CALL Mfield(Dnew(1,k),Bnew(1,k),k,nrx)
c
   20 continue
c
c Update Coulomb Field:
c
      call Mfield(Dnew(1,1),bnew(1,IBX),IBX,nrx)
c
c Updating IsoScalar and Isovector Scalar densities to new kf(r;q):
c Solve coupled Scalar Field equations in INFINITE MATTER approximation:
c
      Mode=1
      If(Mode.eq.1)then
      do n=1,NRX
      Xold(n,1)=Xstar(n,1)
      Xold(n,2)=Xstar(n,2)
      do i=1,IBX
      Bold(n,i)=Bfld(n,i)
      enddo
      Zfp=Xkf(n,1)
      Zfn=Xkf(n,2)
      XMp=Xstar(n,1)
      XMn=Xstar(n,2)
      dvp=Dnew(n,1)
      dvn=Dnew(n,2)
c
      CALL SCD(Zfp,Zfn,dvp,dvn,XMp,XMn,dsp,dsn)
c
      Xstar(n,1)=XMp
      Xstar(n,2)=XMn
      Snew(n,1)=dsp
      Snew(n,2)=dsn
c
      s00=-0.5*(XMp+XMn-RM(1)-RM(2))
      s01=-0.5*(XMp-XMn-RM(1)+RM(2))
      Bnew(n,1)=s00
      Bnew(n,2)=s01
c
      enddo
c
c Update IsoScalar and IsoVector Scalar Densities:
c
      do n=1,nrx
      Snew(n,3)=Snew(n,1)+Snew(n,2)
      Snew(n,4)=Snew(n,1)-Snew(n,2)
      enddo
c
      endif
c
c Update IsoScalar and IsoVector Scalar Fields:
c
      CALL Sfield(xkf,Snew(1,3),Bnew,xstar,1,nrx)
      CALL Sfield(xkf,Snew(1,4),Bnew,xstar,2,nrx)
c
c Update proton and neutron Scalar Densities:
c
      do n=1,nrx
      Snew(n,1)=0.5d0*(Snew(n,3)+Snew(n,4))
      Snew(n,2)=0.5d0*(Snew(n,3)-Snew(n,4))
      enddo
c
c
c  check convergence of vector (rho_v) and scalar (rho_s)
c  densities :   <(rho(old)-rho(new))**2>/<rho(old)**2>  <  epst
c
      sumV=0.
      sumS=0.
      Vint=0.
      Sint=0.
      rms(1)=0.d0
      rms(2)=0.d0
      xms(1)=0.d0
      xms(2)=0.d0
      DO 40 n=1,nrx
      WR=wst(n)
c      WR=1.d0
      sumV=sumV+WR*(Dnew(n,3)-rho_v(n,3))**2
      sumS=sumS+WR*(Snew(n,3)-rho_s(n,3))**2
      Vint=Vint+WR*rho_v(n,3)**2
      Sint=Sint+WR*rho_s(n,3)**2
      rsq=rst(n)*rst(n)
* <r^2> for p/n :
      rms(1)=rms(1)+rsq*wst(n)*Dnew(n,1)
      rms(2)=rms(2)+rsq*wst(n)*Dnew(n,2)
      xms(1)=xms(1)+    wst(n)*Dnew(n,1)
      xms(2)=xms(2)+    wst(n)*Dnew(n,2)
   40 CONTINUE
c
      do n=1,NRX
c
c UPDATE densities:
c
      do k=1,4
      rho_v(n,k)=Dnew(n,k)
      rho_s(n,k)=Snew(n,k)
      enddo
c
c UPDATE Fields:
c
      do k=1,IBX
      Bfld (n,k)=Bnew(n,k)
      enddo
c
      enddo
c
c Check Performance:
c
      sumV=sumV/Vint
      sumS=sumS/Sint
c
      xms(3)=xms(1)+xms(2)
      rms(3)=(rms(1)+rms(2))/xms(3)
      rms(1)=rms(1)/xms(1)
      rms(2)=rms(2)/xms(2)
c
      rms_pn(Iter,1)=rms(1)
      rms_pn(Iter,2)=rms(2)
c
      perfb(iter)=sumV
      perfs(iter)=sumS
c
c Accumulated results of iteration sequence:
c Swichted off for KOUT(1) > 0 (Input option!)
c
      IF(KOUT(1).eq.0)then
c
      write(7,701)sumV,sumS
      write(7,659)iter
  701 format('#',3x,'drhob=',e12.4,3x,'drhos=',e12.4)
cTG  659 format('#'10x,'RESULTS FOR ITERATION No.',i4/
  659 format('#',10x,'RESULTS FOR ITERATION No.',i4/ !replaced by TG: missing comma added
     1 '#',8x,'r',6x,'rhov_0',6x,'rhos_0',7x,'omega',7x,'sigma',
     & 6x,'Umf(nr)')
c
      DO 45 n=1,nrx,2
      Umf=Bnew(n,3)-Bnew(n,1)
      Vp=Bnew(n,3)+Bnew(n,4)+Bnew(n,IBX)
      Vn=Bnew(n,3)-Bnew(n,4)
      Ep=Bold(n,3)+Bold(n,4)+Bold(n,IBX)+sqrt(Xkf(n,1)**2+Xold(n,1)**2)
      En=Bold(n,3)-Bold(n,4)            +sqrt(Xkf(n,2)**2+Xold(n,2)**2)
      write(7,661)rst(n),Dnew(n,3),Snew(n,3),Bnew(n,3),Bnew(n,1),Ep,En
   45 CONTINUE
      write(7,7023)
 7023 format()
      EndIf
c
  661 format(f9.2,10e12.4)
      If(Ipr.gt.0)write(8,7001)(sqrt(rms(k)),k=1,3)
 7001 format(/'  <r^2> :',3f10.5,' [fm]')
c
      Dlam=10.
      If(Iter.gt.1)then
      DLamP=(Vlam(Iter-1,1)-Vlam(iter,1))/Vlam(iter-1,1)
      DLamN=(Vlam(Iter-1,2)-Vlam(iter,2))/Vlam(iter-1,2)
      Dlam=0.5d0*(abs(DlamP)+abs(DlamN))
      endif
c
cTG      tts(iter)=timer()-tt1
c
      If(Dlam.lt.eps)goto 9000
      if(sumV.gt.eps.or.sumS.gt.eps)goto 1
c
 9000 close(7)
      If(MVch.eq.0)write(8,716)
      If(MVch.eq.1)write(8,718)
  716 format(//10x,'          Relativistic Thomas-Fermi Calculation')
  718 format(//10x,'EXTENDED Relativistic Thomas-Fermi Calculation')
      write(8,720)
  720 format(/5x,'PERFORMANCE of the Iteration Sequence - RMS-DEV. of',
     & ' Vector and Scalar Densities:'//
     1 5x,'iter.',7x,'drhob',7x,'drhos',4x,'<r^2>(p)',4x,'<r^2>(n)',
     & 3x,'Lambda(p)',3x,'Lambda(n)',
     & 2x,'Time (sec)')
      DO 300 k=1,min(iter,ITERX)
  300 write(8,721)k,perfb(k),perfs(k)
     & ,sqrt(rms_pn(k,1)),sqrt(rms_pn(k,2)),Vlam(k,1),Vlam(k,2),tts(k)
  721 format(5x,i5,2e12.4,2f12.4,2e12.4,f12.4)
c
cTG      ttotal=timer()-tinit
cTG      write(8,722)ttotal
cTG  722 format(/9x,'TIME for Iterations:',12x,f12.3,' secs')
      If(Iter.le.ITERX)return
 9005 write(8,9600)iter-1
 9600 format(9x,'NO SELFCONSISTENCY OBTAINED IN',i4,' ITERATIONS')
      return
      end
c---------------------------------------------------------------------
      SUBROUTINE VchD(rho_v,Xkf,Bfld,Vch,Upot)
c      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      Parameter (KRX=400,KBX=6,ITERX=50,ITERDX=200)
C
c Extended RTF: Using
c    M*^2 Psibar gamma_0 Psi = alpha*{(pslash Psibar) gamma_0 (pslash Psi)}+beta*{p^2(Psibar gamma_0 Psi)}
c  where alpha+beta=1 and with the infinite matter relation
c    M*^2  = (lambda-V)^2-k_F^2
c  a (NON-LINEAR) 2nd order DEQ is derived
c    d''-2b(kf(rho)^2+M*^2-(lambda-V0)^2)d=0
c  for the proton/neutron vector densities, where
c    d(r):=r*rho_v(r),
c    b   :=1/(1-beta)=1/alpha.
c  For beta=1/2, kappa^2=M^2-lambda^2 and r -> infinity:
c    d(r) -> exp(-2*kappa*r) .
c kf,M*, lambda and V0 are in units of M, i.e. dimensionless;
c rho_v(r) is in units of (hbc/M)^3, i.e. dimensionless;
c --> leads to a factor of (M/hbc)^2 for the potential term!
C
      common/pi4pi/pi,fpi,pisq
      common/ScalarD/rho_s(KRX,4),xstar(KRX,2)
      common/boson/bmass(KBX),Jbo(KBX),Ibo(KBX),ISV(4)
      common/couple/G4pi(KBX),bcpl(KBX),rnge(KBX),ibx
      common/Nucl/ Amass,AN,AZ,b(2),rdA
      common/const/Xmss(4),Xmass,alfa,hbc,hbcm,hbcm3
      common/radius/rmax,dr,rst(KRX),rwgt(KRX),nrx
      Dimension Vch(2),F(KRX),U(KRX),EFSQ(KRX)
      Dimension W(2*ITERX+1),Vst(2*ITERX+1)
      Dimension Vfld(KRX),DelD(ITERDX)
      Dimension rho_v(KRX,4),Xkf(KRX,2),Bfld(KRX,KBX),Upot(KRX,2)
      Data betaL,eps,epsD/0.5d0,1.d-08,1.d-08/
      Data Xid/0.15d0/
c
      hbcm2=hbcm*hbcm
      Fmsq=1.d0/hbcm2
      Yid=1.-Xid
      H=Rst(2)-Rst(1)
      Dfac=3*pisq
      pwr=1.d0/3
      Match=rdA/h+1.1
      bL=1.d0/(1.d0-betaL)
      Alfa0=0.
      MFail=0
      MPRID=-1
c
c Isospin Loop:
c
      do 1000 I3=1,2
c
c Total Vector Fields:
c
      If(I3.eq.1)then
      do n=1,NRX
      Vfld(n)=Bfld(n,3)+Bfld(n,4)+Bfld(n,5)
      enddo
      endif
c
      If(I3.eq.2)then
      do n=1,NRX
      Vfld(n)=Bfld(n,3)-Bfld(n,4)
      enddo
      endif
c
c Initialize:
c
      Vch(i3)=0.d0
      do n=1,NRX
      Ef=sqrt(Xkf(n,i3)**2+Xstar(n,i3)**2)
      Epn=Vfld(n)+Ef
c      Vch(i3)=Vch(i3)+rwgt(n)*rho_v(n,I3)*Epn
      enddo
c      Vch(i3)=fpi*Vch(i3)/(B(i3)*hbcm3)
c
c Iteration for non-linearity due to kf(rho):
c
      Id=0
   10 Id=Id+1
c
      M=0
c
      do n=1,NRX
      EFSQ(n)=Xkf(n,i3)**2+Xstar(n,i3)**2
      enddo
c
c Search for Chemical Potential for given potential u:
c
c      Vch(i3)=1.d0-30./XMass
      Vch(i3)=-20.d0
      dV=abs(Vch(i3)-1.d0)/ITERX
      Vpn=Vch(i3)-(Iterx-10)*dV
      do I=1,2*ITERX+1
      Vpn=Vpn+dV*1.5
      Vst(i)=Vpn
      V2=Vpn-1.d0
c
      do n=1,NRX
      U(n)=-2*bL*(EFSQ(n)-(Vpn-Vfld(n))**2)*Fmsq
      enddo
      write(*,*)'  V2:',V2
      Wdet=0.d0
c
c      call H2EQ(NRX,U,F,WDet,V2,NODE,H,Alfa0,MATCH)
c
      W(i)=WDet
c
      If(i.gt.1.and.M.le.0)then
       If(W(i)*W(i-1).lt.0.d0)M=i
      endif
      enddo
c
      If(M.gt.0) goto 100
c
      write(*,*)'   Initialization for Chem. Pot. Lambda failed! '
      write(8,*)'            -->    sub VchD:  <--'
      write(8,*)'   Initialization for Chem. Pot. Lambda failed! '
      write(8,*)' Wronski-Det for Matching on data set -> wronski.d'
      open(15,file='wronski.d',status='UNKNOWN')
      write(15,1501)
      do i=1,ITERX
      write(15,1500)i,Vst(i),W(i)
      enddo
 1500 format(i5,10e14.6)
 1501 format('# No.',6x,'lambda/M',6x,'WronskiD')
      close(15)
      MFail=1
      stop
c
  100 continue
c
      If(MPRID.gt.0)then
      If(i3.eq.1)open(15,file='wdet_p.d',status='UNKNOWN')
      If(i3.eq.2)open(15,file='wdet_n.d',status='UNKNOWN')
      write(15,1501)
      do i=1,2*ITERX+1
      write(15,1500)i,Vst(i),W(i)
      enddo
      close(15)
      endif
c
      V1=Vst(m-1)
      V2=Vst(m  )
      W1=W(m-1)
      W2=W(m  )
      dW=(W2-W1)/(V2-V1)
c
      If(MPRID.gt.0)then
      write(*,628)i3,id
      write(8,628)i3,id
  628 format(/' Lambda Iteration - i3:',i3,' for density step:',i3/
     & '  No.',7x,'lambda',6x,'dlambda',9x,'W/dW')
      endif
c
      I=0
  110 I=I+1
c
      If(MPRID.gt.0)then
      write(*,629)I,V2,dV,W2/dW
      write(8,629)I,V2,dV,W2/dW
  629 format(i5,10e13.5)
      endif
c
      If(I.gt.ITERX)then
       write(*,*)' Search for Chemical Potential failed!'
       stop
      endif
      dW=(W2-W1)/(V2-V1)
      dV=-W2/dW
      If(abs(dV/V2).lt.eps.and.abs(dV).lt.eps)goto 120
      V0=V2+dV
c
      do n=1,NRX
      U(n)=-2*bL*(EFSQ(n)-(V0-Vfld(n))**2)*Fmsq
      enddo
      V0P=V0-1.d0
c
c      call H2EQ(NRX,U,F,W0,V0P,NODE,H,Alfa0,MATCH)
c
c      If(W1*W0.lt.0.d0)then
c       W2=W0
c       V2=V0
c      endif
c      If(W2*W0.le.0.d0)then
c       W1=W0
c       V1=V0
c      endif
c
      If(V1.lt.Vst(m-1))then
        Write(*,*)' WARNING: V1 < V(m-1)!!'
        write(*,*)' V1      ,V2    :',V1,V2
        write(*,*)' Vst(m-1),Vst(m):',Vst(m-1),Vst(m)
        write(*,*)' W1      ,W2    :',W1,W2
        write(*,*)' W(m-1)  ,W(m)  :',W(m-1),W(m)
        stop
      endif
c
      If(V2.gt.Vst(m))then
        Write(*,*)' WARNING: V2 > V(m)!!'
        write(*,*)' V1      ,V2    :',V1,V2
        write(*,*)' Vst(m-1),Vst(m):',Vst(m-1),Vst(m)
        write(*,*)' W1      ,W2    :',W1,W2
        write(*,*)' W(m-1)  ,W(m)  :',W(m-1),W(m)
        stop
      endif
c
      If(Node.ne.1)  then
      write(*,*)' Node:',node
      write(8,*)' Node:',node
      write(8,*)' sub VchD: rho_v with higher nodes!'
      write(*,*)' m :',m
      write(*,*)' V1,V2:',Vst(m-1),Vst(m)
      write(*,*)' W1,W2:',W(m-1),W(m)
      write(8,*)' m :',m
      write(8,*)' V1,V2:',Vst(m-1),Vst(m)
      write(8,*)' W1,W2:',W(m-1),W(m)
      write(8,*)' ChemP,W:',V2,W2
      write(8,*)'     check data set >>dpn.d<< '
      open(15,file='wronski.d',status='UNKNOWN')
      write(15,1501)
      do i=1,ITERX
      write(15,1500)i,Vst(i),W(i)
      enddo
      do n=1,NRX
      F(n)=F(n)/rst(n)**(1+Alfa0)
      Dnorm=Dnorm+Rwgt(n)*F(n)
      enddo
      Dnorm=B(I3)*hbcm3/(fpi*Dnorm)
      If(F(2)*Dnorm.lt.0.d0)Dnorm=-Dnorm
      open(23,file='dpn.d',status='UNKNOWN')
      DO n=1,NRX
      F(n)=Dnorm*F(n)
      rho_v(n,i3)=Xid*F(n)+Yid*rho_v(n,i3)
      write(23,2300)rst(n),F(n),rho_v(n,i3),Upot(n,i3)
      enddo
      stop
      endif
c
      goto 110
c
  120 Vch(i3)=V2
      Dnorm=0.d0
c
c
      do n=1,NRX
      F(n)=F(n)/rst(n)**(1+Alfa0)
      Dnorm=Dnorm+Rwgt(n)*F(n)
      enddo
      Dnorm=B(I3)*hbcm3/(fpi*Dnorm)
      If(F(2)*Dnorm.lt.0.d0)Dnorm=-Dnorm
      DelD(id)=0.d0
c
      do n=1,NRX
      F(n)=Dnorm*F(n)
      DelD(id)=DelD(id)+(F(n)-rho_v(n,I3))**2
      rho_v(n,i3)=Xid*F(n)+Yid*rho_v(n,i3)
      Upot (n,i3)=U(n)
      Xkf  (n,i3)=abs(Dfac*rho_v(n,i3))**pwr
      enddo
c
      If(MFail.ne.0)stop
c
      DelD(id)=DelD(id)*fpi
      If(MPRID.gt.0)then
      write(*,629)Id,DelD(id),Vch(i3)
      write(8,629)Id,DelD(id),Vch(i3)
      endif
      If(DelD(id).lt.epsD)goto 1000
      If(ID.lt.IterDx)goto 10
      write(8,*)' sub VchD: rho_v did not converge!'
      write(8,*)'     check data set >>dpn.d<< '
c
      open(23,file='dpn.d',status='UNKNOWN')
c
      DO n=1,NRX
      write(23,2300)rst(n),F(n),rho_v(n,i3),Upot(n,i3)
      enddo
 2300 format(f8.2,10e13.5)
      close(23)
      stop
c
 1000 continue
c
      return
      end
c---------------------------------------------------------------------
      SUBROUTINE H2EQ(NDIM,U,F,WDet,E,NODE,H,Alfa0,MATCH)
C
C  SOLUTION OF A HOMOGENOUS 2.ND-ORDER DIFF.EQ. WITH
C  THE MODIFIED NUMEROV METHOD:
c         f''+u(r,E)f=0
C
c      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      Parameter (KRX=400)
      common/radius/rmax,dr,rst(KRX),rwgt(KRX),nrx
      DIMENSION U(NDIM),F(NDIM),D(2),S(KRX),W(KRX)
      DATA EPS/1.D-02/
c
      H2=H*H
      Q=H2/12.D0
      z=1.+Rst(2)/Rst(1)
      AA=Alfa0*(Alfa0+1)
c
      do n=1,ndim
       rsq=rst(n)*rst(n)
       W(n)=U(n)-AA/rsq
       s(n)=(1.d0-Q*W(n))
      enddo
C
C  START AT THE ORIGIN
C
      G2=0.0
      F(1)=(abs(W(1))*Rst(1))**(1+Alfa0)
      G1=F(1)*(2-S(1))
      G2=F(1)*(1.-(1+Alfa0)*z+0.5*AA*z**2)*(2-S(2))
      DO 100 N=2,MATCH+1
      G=2.*G1-G2-H2*F(N-1)*W(N-1)
      F(N)=G*S(n)
      G2=G1
  100 G1=G
C
C FIX THE PHASE TO BE >0 AT THE ORIGIN
C
      IF(F(2).LT.0.D0)THEN
      DO N=1,MATCH+1
       F(N)=-F(N)
      enddo
      ENDIF
C
C Return if for Positive E the full asymptotic part is required :
c
      IF(MATCH.ge.NDIM-10)goto 399
C
C  LOGARITHMIC DERIVATIVE (INTERNAL SOLUTION)
C  (EQUI-DISTANT 3-POINT FORMULA)
C
C     IF(F(MATCH).EQ.0.0)MATCH=MATCH-1
      INTX=MATCH
      FINT=F(INTX)
      D(1)=(F(INTX+1)-F(INTX-1))*0.5/(H*Fint)
      F2=F(INTX-1)
c      FMA(1)=FINT
c      DMA(1)=FINT*D(1)
C
C  EXTERNAL:
C
      NEXTR=NDIM-INTX+1
      if(E.ge.0.d0)then
       F(NDIM  )=0.0
       F(NDIM-1)=EPS
      else
       zn=sqrt(abs(W(NDIM  )))
       zm=sqrt(abs(W(NDIM-1)))
       F(NDIM  )=eps
       F(NDIM-1)=eps*exp((NDIM*(zn-zm)+zm)*h)
      endif
      G2=F(NDIM  )*(2.d0-S(ndim  ))
      G1=F(NDIM-1)*(2.d0-S(ndim-1))
      DO 200 N=2,NEXTR
      K=NDIM-N
      G=2.*G1-G2-H2*F(K+1)*W(K+1)
      F(K)=G*S(K)
      G2=G1
  200 G1=G
C
C  LOGARITHMIC DERIVATIVE (EXTERNAL SOLUTION)
C
      Fext=F(intx)
      D(2)=(F(INTX+1)-F(INTX-1))*0.5/(H*Fext)
      WDet=D(1)-D(2)
      F(INTX-1)=F2
c
      A=1.d0/Fint
      do n=1,INTX-1
      F(n)=A*F(n)
      enddo
c
      A=1.d0/Fext
      DO N=INTX,NDIM
      F(N)=A*F(N)
      enddo
* Normalize:
      Wnorm=0.d0
      do n=1,NDIM
      WNorm=Wnorm+F(n)**2
      enddo
      Wnorm=(Wnorm-0.5*(F(NDIM)**2+F(1)**2))*H
      If(Wnorm.le.0.d0)then
        write(*,*)' Error in Wnorm:',Wnorm
        stop
      endif
      Xnorm=1.d0/sqrt(Wnorm)
c
      do n=1,NDIM
      F(n)=F(n)*Xnorm
      enddo
c
      WDet=F(intx)*WDet*Xnorm
c
  399 NODE=1
      DO 400 N=2,NDIM-1
      IF(F(N)*F(N-1).GT.0.D0)GOTO 400
      NODE=NODE+1
  400 CONTINUE
c
      RETURN
      END
c---------------------------------------------------------------------
      SUBROUTINE VCTF(nrx,wst,Vpn,xstar,rho,AM,qf,x)
c
c Chemical Potential in Thomas-Fermi Approximation for q=p,n:
c  lambda_q = V0_q + sqrt(kf^2(q)+M*^2_q) ==> kf^2(q) = (lambda_q-V0_q)^2-M*^2_q
c  rho_q(lambda)=kf^3(q)/3pi^2
c  and giving rho=rho(r;lambda)
c  where lambda is fixed by <rho>=B  (B:=AM in input list!)
c ==>  M*(q) and kf(q) in units of M !
c
c      implicit double precision (a-h,o-z)
      Parameter (Ipr=0,KRX=400)
      common/const/Xmss(4),Xmass,alfa,hbc,hbcm,hbcm3
      common/pi4pi/pi,fpi,pisq
      common/radius/rmax,dr,rst(KRX),rwgt(KRX),nrxx
      dimension wst(nrx),Vpn(KRX),xstar(KRX),rho(KRX),qf(KRX)
      dimension s(3)
      DATA eps,iterx/1.e-06,50/
c
      f(xx,bb,cc)=xx*(xx-2.*bb)-cc
      estar(zz,qq)=sqrt(zz**2+qq**2)
c
      If(ipr.gt.0)write(8, 698)
  698 format(/8x,'ITERATION FOR CHEMICAL POTENTIAL')
c
      amm=AM*hbcm3
      A=4./(3.*pi*amm)
c
      xl=0.0
      DO 5 n=1,nrx
    5 xl=max(estar(xstar(n),qf(n))+Vpn(n),xl)
      x=min(0.9,0.85*xl)
c
      dx=1.2*x/iterx
      nsuch=0
      iter=0
    1 iter=iter+1
      if(iter.gt.iterx)then
       write(*,*)' sub VCTF: Initialization for lambda failed!'
       write(8,*)' sub VCTF: Initialization for lambda failed!'
       stop
      endif
c
      DO 3 i=1,3
    3 s(i)=0.
c
c Integrals S_i:
c
      DO n=1,nrx
      sqm=xstar(n)**2
      Vsq=  Vpn(n)**2
      qs=(x-Vpn(n))**2-sqm
      qs=Max(0.,qs)
      q=A*sqrt(qs)
c
c S_1 and S_2:
c
      p=1.
      DO i=1,2
      s(i)=s(i)+wst(n)*q*p
      p=p*Vpn(n)
      enddo
c
      s(3)=s(3)+wst(n)*q*(sqm-Vsq)
c
      enddo
c
      if(s(1).eq.0.0)s(1)=eps
      b=s(2)/s(1)
      c=(1.+s(3))/s(1)
      f1=f(x,b,c)
      if(nsuch.eq.1)goto 30
      if(iter.eq.1)goto 28
      nsuch=1
      If(IPR.gt.0.and.s(1).ne.eps)then
       write(8,629)iter,x,f1
  629  format(i4,8e13.5)
      endif
      if(f0*f1.le.0.0)goto 30
   28 nsuch=0
      x0=x
      f0=f1
      x=x+dx
      goto 1
   30 df=(f1-f0)/(x-x0)
      if(df.eq.0.)then
      write(*,*)' ==> df == 0 !'
      goto 40
      endif
      dx=-f0/df
      z=x0+dx
      x0=x
      f0=f1
      x=z
      If(ipr.gt.0)write(8, 699)iter,x0,f0,dx
  699 format(3x,'i:',i4,'  x0=',e12.4,'  f0=',e12.4,'  dx=',e12.4)
      if(abs(dx).gt.eps.or.abs(f1).gt.eps)goto 1
   40 f1=f(x,b,c)
c
c RTF densities:
c
      Dfac=1.d0/(3*pisq)
      cmass=0.d0
      DO 100 n=1,nrx
      qs=(x-Vpn(n))**2-xstar(n)**2
      qs=Max(0.,qs)
      qf(n)=sqrt(qs)
      rho(n)=Dfac*qs*qf(n)
c      write(8,697)rst(n),rho(n),qf(n),qs,Vpn(n),Xstar(n)
  697 format(f8.2,10e13.5)
      cmass=cmass+wst(n)*rho(n)
  100 CONTINUE
      cmass=cmass*fpi/hbcm3
  600 format('  iterat.:',i4,3x,'mass=',f8.2,3x,'lambda=',e12.4,3x,
     & 'f(lam)=',e12.4)
      If(ipr.gt.0)write(8,600)iter,cmass,x,f1
      return
 9000 write(8,9600)iter,x0,f0,x,f1,dx
 9600 format(3x,'NO CONVERGENCE FOR LAMBDA IN',i4,' ITERATIONS'/
     1 3x,'last values:',
     1 'x0,f0=',2e12.4,3x,'x1,f1,=',2e12.4,3x,'dx=',e12.4)
      return
      end
c---------------------------------------------------------------------
      SUBROUTINE Mfield(rho,phi,mboson,nrx)
c
c Compute Mean-Fields for given source term by Green Function Method
c
c      implicit double precision (a-h,o-z)
      Parameter (KBX=6,KRX=400)
      common/const/Xmss(4),Xmass,alfa,hbc,hbcm,hbcm3
      common/boson/bmass(KBX),Jbo(KBX),Ibo(KBX),ISV(4)
      common/couple/G4pi(KBX),bcpl(KBX),rnge(KBX),ibx
      common/bnonl/NLINT(KBX),NLPOW(4,KBX),ANL(4,KBX),XNL(4,KBX)
      common/radius/rmax,dr,rst(KRX),rwgt(KRX),nrxx
      common/green/b(KRX,KBX),h(KRX,KBX),gnn(KRX,KBX)
      dimension s(KRX),s1(KRX),s2(KRX),t(KRX)
      dimension rho(nrx),phi(nrx)
c
      hsq=dr*dr/12.
      drb=dr
      if(rnge(mboson).gt.0.d0)then
      drb=drb/rnge(mboson)**3
      hsq=hsq/rnge(mboson)**2
      else
      drb=drb/hbcm**2
      hsq=hsq/hbcm**2
      endif
c
c  compute the source term :
c
      DO 15 n=1,nrx
      s(n)=rho(n)*bcpl(mboson)
   15 t(n)=rst(n)**2*s(n)
c
c  compute basic integrals :
c (Progators include (M/hbc)^3 such that fields are dimensionless!)
c
      s1(1)=0.
      s2(nrx)=0.
      s2total=0.
      DO 30 n=1,nrx-1
      k=nrx-n
      s2total=s2total+t(n)*h(n,mboson)*drb
      s1(n+1)=s1(n)+t(n+1)*b(n+1,mboson)*drb
   30 s2(k)=s2(k+1)+t(k+1)*h(k+1,mboson)*drb
c
c  take care of the special procedure for r=0.
c
      s(1)=-s(1)*b(1,mboson)
      aa=s2total/s2(1)-1.d0
      if(abs(aa).gt.1.e-04)write(8, 665)mboson,s2(1),s2total,aa
  665 format(3x,'numerical instability - Boson:',i3,
     & '  s2(1),s2total=',2e13.5,3x,'ratio=',e13.5/)
      s2(1)=s2total
c
c
      DO 35 n=1,nrx
      phi(n)=h(n,mboson)*s1(n)+b(n,mboson)*s2(n)-hsq*s(n)
   35 CONTINUE
      return
      end
c----------------------------------------------------------------------
      SUBROUTINE greenf
c      implicit double precision (a-h,o-z)
      Parameter (KBX=6,KRX=400)
      common/const/Xmss(4),Xmass,alfa,hbc,hbcm,hbcm3
      common/boson/bmass(KBX),Jbo(KBX),Ibo(KBX),ISV(4)
      common/couple/G4pi(KBX),bcpl(KBX),rnge(KBX),ibx
      common/bnonl/NLINT(KBX),NLPOW(4,KBX),ANL(4,KBX),XNL(4,KBX)
      common/radius/rmax,dr,rst(KRX),rwgt(KRX),nrx
      common/green/b(KRX,KBX),h(KRX,KBX),gnn(KRX,KBX)
c
c  compute the propagators in r-space
c
      rnge(IBX)=0.d0
      DO 100 i=1,ibx
      if(rnge(i).eq.0.d0)goto 50
      r=1./rnge(i)
      dm=r*dr
      DO 10 n=1,nrx
      x=r*rst(n)
c
c  Gnn is the Diagonal Part of the Integral-Equation Kernel
c
      gnn(n,i)=(exp(-x)*sinh(x)-dm)*dm*bcpl(i)
      h(n,i)=exp(-x)/x
      if(abs(x).lt.1.e-04)goto 5
      b(n,i)=sinh(x)/x
      goto 10
    5 b(n,i)=1.+(x**2/6)*(1.+x**2/20)
      h(n,i)=0.0
   10 CONTINUE
      goto 100
c
c Coulomb:
c
   50 gnn(n,i)=0.d0
      DO 20 n=1,nrx
      b(n,i)=1.d0
   20 h(n,i)=1.d0/rst(n)
c
  100 CONTINUE
c
      return
      end
c-----------------------------------------------------------------------
cTG      SUBROUTINE input(wkdir,prnout,in)
      SUBROUTINE input(prnout,in)
c      implicit double precision (a-h,o-z)
      Parameter (KBX=6,KRX=400)
cTG      CHARACTER*255           wkdir
      character*125 prnout
      common/cntrl/kout(10),epst
      common/pi4pi/pi,fpi,pisq
      common/boson/bmass(KBX),Jbo(KBX),Ibo(KBX),ISV(4)
      common/couple/G4pi(KBX),bcpl(KBX),rnge(KBX),ibx
      common/bnonl/NLINT(KBX),NLPOW(4,KBX),ANL(4,KBX),XNL(4,KBX)
      common/const/Xmss(4),Xmass,alfa,hbc,hbcm,hbcm3
      common/radius/rmax,dr,rst(KRX),rwgt(KRX),nrx
      common/Nucl/ Amass,AN,AZ,b(2),rd
      Dimension INDB(20)
c
      Data XMp,XMn/938.27231d0,939.56563d0/

c$$$      write(*,102)
c$$$  102 format(5x,'filename for input data:')
c$$$ctg      read(*,*,end=1010)dataf
c$$$ctg 1010 open (5,file=dataf,status='OLD')
c$$$
c$$$cTG      open (5,file='ertf.input',status='OLD')
c$$$      rewind(5)
c$$$      open (5,file='ertf.input')
c$$$      dataf = 'ertf.input'

c
      Xmss(1)=XMp
      XMss(2)=XMn
      Xmss(3)=0.5d0*(Xmss(1)+Xmss(2))
      Xmss(4)=0.5d0*(Xmss(1)-Xmss(2))
      Xmass=Xmss(3)
c
      pi=4*atan(1.d0)
      fpi=4*pi
      pisq=pi*pi
c
      HBC =197.327053d0
      hbcm=hbc/Xmass
      hbcm3=hbcm**3
      ALFA=137.035989d0
c
      do I=1,KBX
      NLINT(i)=0
      Bcpl(i)=0.d0
      do K=1,4
      NLPOW(k,i)=0
      anl  (k,i)=0.d0
      enddo
      enddo
c
      do n=1,10
      Kout(n)=0
      enddo
c
      DO 6 i=1,1
      read(in,*)Amass,AZ
    6 AN=Amass-AZ

      if(amass.le.0.d0)stop
c Rmax, Dr
      read(in,*)rmax,dr
      Nmax=Rmax/dr
      Nrx=2*(Nmax/2)+1
c file name (store results for print-output):
      read(in,*)prnout
      open (8,file=prnout,status='UNKNOWN')
c Working Directory
cTG      read(in,*)wkdir
c
cTG      call MKD(wkdir)
c
c empty line (previously occupied by KOUT):
c
      read(in,*) text
c
c  accuracy:
c
      read(in,*)epst
c
c Mass,Spin,Isospin,g^2/4pi,non-linearity
c NOTE:
c Jbo:= 0  for Scalar Mesons
c Jbo:= 1  for Vector Mesons
c Ibo:= 0  for Iso-Scalar Mesons
c Ibo:= 1  for Iso-Vector Mesons
c Ibo:=-1  for the photon (and Bmass:=0)
c
      I=0
   20 i=i+1
      read(in,*)bmass(i),Jbo(i),Ibo(i),G4pi(i),NLINT(i)
      If(Bmass(i).lt.0.d0)goto 30
      INDB(i)=i
C
c Photon:
c
      If(Bmass(i).eq.0.d0)then
      G4pi(i)=1.d0/ALFA
      BCPL(I)=fpi*G4pi(i)
      Jbo (I)= 1
      Ibo (I)=-1
      endif

c
c non-linear couplings:
c
      if(NLINT(i).ne.0)then
      DO 25 j=1,NLINT(i)
      read(in,*)NLPOW(j,i),ANL(j,i)
   25 CONTINUE
c
      endif
c
      goto 20
c
   30 IBX=I-1
c
c "Canonical" Order of Bosons:
c --> Scalar-IsoScalar,Scalar-IsoVector,Vector-IsoScalar,Vector-IsoVector,Photon
c
      KB=1
      I1=0
      do I=1,IBX
c
c Scalar-IsoScalar (KB=1):
c
      If(Jbo(i).eq.0.and.Ibo(i).eq.0)then
      I1=I1+1
c
c Count the No. of Bosons of the same type:
c
      ISV(KB)=I1
c
      BI=BMass(i)
      GI=G4pi(i)
      JI=Jbo(i)
      II=Ibo(i)
      NI=NLINT(i)
c
      BP=BMass(i1)
      GP=G4pi(i1)
      JP=Jbo(i1)
      IP=Ibo(i1)
      NP=NLINT(i1)
c
      BMass(i1)=BI
      G4pi (i1)=GI
      JBo  (i1)=JI
      IBo  (i1)=II
      NLINT(i1)=NI
      BMass(i )=BP
      G4pi (i )=GP
      JBo  (i )=JP
      IBo  (i )=Ip
      NLINT(i)=NP
      do j=1,4
      np =NLPOW(j,i1)
      ap =ANL(j,i1)
      anl   (j,i1)=anl   (j,i)
      NLPOW(j,i1)=NLPOW(j,i)
      anl   (j,i )=ap
      NLPOW(j,i )=np
      enddo
      endif
c
      enddo
c
      KB=2
      I1=ISV(KB-1)
      do I=ISV(KB-1)+1,IBX
c
c Scalar-IsoVector (KB=2):
c
      If(Jbo(i).eq.0.and.Ibo(i).eq.1)then
      I1=I1+1
c
c Count the No. of Bosons of the same type:
c
      ISV(KB)=I1-ISV(KB-1)
c
      BI=BMass(i)
      GI=G4pi(i)
      JI=Jbo(i)
      II=Ibo(i)
      NI=NLINT(i)
c
      BP=BMass(i1)
      GP=G4pi(i1)
      JP=Jbo(i1)
      IP=Ibo(i1)
      NP=NLINT(i1)
c
      BMass(i1)=BI
      G4pi (i1)=GI
      JBo  (i1)=JI
      IBo  (i1)=II
      NLINT(i1)=NI
      BMass(i )=BP
      G4pi (i )=GP
      JBo  (i )=JP
      IBo  (i )=Ip
      NLINT(i)=NP
      do j=1,4
      np =NLPOW(j,i1)
      ap =ANL(j,i1)
      anl   (j,i1)=anl   (j,i)
      NLPOW(j,i1)=NLPOW(j,i)
      anl   (j,i )=ap
      NLPOW(j,i )=np
      enddo
      endif
c
      enddo
c
      KB=3
      I1=ISV(KB-1)
      do I=ISV(KB-1)+1,IBX
c
c Vector-IsoScalar (KB=3):
c
      If(Jbo(i).eq.0.and.Ibo(i).eq.1)then
      I1=I1+1
c
c Count the No. of Bosons of the same type:
c
      ISV(KB)=I1-ISV(KB-1)
c
      BI=BMass(i)
      GI=G4pi(i)
      JI=Jbo(i)
      II=Ibo(i)
      NI=NLINT(i)
c
      BP=BMass(i1)
      GP=G4pi(i1)
      JP=Jbo(i1)
      IP=Ibo(i1)
      NP=NLINT(i1)
c
      BMass(i1)=BI
      G4pi (i1)=GI
      JBo  (i1)=JI
      IBo  (i1)=II
      NLINT(i1)=NI
      BMass(i )=BP
      G4pi (i )=GP
      JBo  (i )=JP
      IBo  (i )=Ip
      NLINT(i)=NP
      do j=1,4
      np =NLPOW(j,i1)
      ap =ANL(j,i1)
      anl   (j,i1)=anl   (j,i)
      NLPOW(j,i1)=NLPOW(j,i)
      anl   (j,i )=ap
      NLPOW(j,i )=np
      enddo
      endif
c
      enddo
c
      KB=4
      I1=ISV(KB-1)
      do I=ISV(KB-1)+1,IBX
c
c Vector-IsoVector (KB=2):
c
      If(Jbo(i).eq.0.and.Ibo(i).eq.1)then
      I1=I1+1
c
c Count the No. of Bosons of the same type:
c
      ISV(KB)=I1-ISV(KB-1)
c
      BI=BMass(i)
      GI=G4pi(i)
      JI=Jbo(i)
      II=Ibo(i)
      NI=NLINT(i)
c
      BP=BMass(i1)
      GP=G4pi(i1)
      JP=Jbo(i1)
      IP=Ibo(i1)
      NP=NLINT(i1)
c
      BMass(i1)=BI
      G4pi (i1)=GI
      JBo  (i1)=JI
      IBo  (i1)=II
      NLINT(i1)=NI
      BMass(i )=BP
      G4pi (i )=GP
      JBo  (i )=JP
      IBo  (i )=Ip
      NLINT(i)=NP
      do j=1,4
      np =NLPOW(j,i1)
      ap =ANL(j,i1)
      anl   (j,i1)=anl   (j,i)
      NLPOW(j,i1)=NLPOW(j,i)
      anl   (j,i )=ap
      NLPOW(j,i )=np
      enddo
      endif
c
      enddo
      return
      end
c---------------------------------------------------------------------
      SUBROUTINE init
c      implicit double precision (a-h,o-z)
      Parameter (KBX=6,KRX=400)
      common/cntrl/kout(10),epst
      common/pi4pi/pi,fpi,pisq
      COMMON/DensP/rd(3),ad(3),rhod(3),RMS(3),RCHRG(3)
      common/boson/bmass(KBX),Jbo(KBX),Ibo(KBX),ISV(4)
      common/couple/G4pi(KBX),bcpl(KBX),rnge(KBX),ibx
      common/bnonl/NLINT(KBX),NLPOW(4,KBX),ANL(4,KBX),XNL(4,KBX)
      common/const/Xmss(4),Xmass,alfa,hbc,hbcm,hbcm3
      common/radius/rmax,dr,rst(KRX),rwgt(KRX),nrx
      common/SVF/bfld(KRX,KBX)
      common/VectorD/rho_v(KRX,4),xkf(KRX,2)
      common/ScalarD/rho_s(KRX,4),xstar(KRX,2)
      common/Nucl/ amass,an,AZ,b(2),rda
      dimension f(2),df(2)
      Dimension Y0(2),Ef(2),Vch(2)
      Dimension Upot(KRX,2)
c
      call DFS(Amass,AZ)
c
      fkf=3.*pisq
      Y0(1)=Xmss(1)/Xmass
      Y0(2)=Xmss(2)/Xmass
      Xsi=AZ/Amass
c
      open(23,file='mfd-a.d',status='UNKNOWN')
c
      write(23,2300)Amass,Xsi
 2300 format('#     Mean-Fields  for A :',f5.1,'  Z/A:',f5.1/
     &       '#',6x,'r',8x,'rho_v',8x,'rho_s',9x,'sigma',9x,'delta',9x,
     &        'omega',11x,'rho',9x,'M*(p)',9x,'M*(n)')
C
c Coulomb:
c
      RC    =Rchrg(3)/sqrt(0.6d0)
      UCEXT =AZ*(1./ALFA)*(HBC/XMASS)
      UCINT =UCEXT/(2.*RC)
C
      smass=1./hbcm3
      rdA=rd(3)
      b(1)=AZ
      b(2)=AN
      do 5 j=1,2
      f(j)=exp(-(dr+rd(j))/ad(j))
    5 df(j)=exp(dr/ad(j))
c
      Ef(1)=0.d0
      Ef(2)=0.d0
c
      Ipr=0
      If(Ipr.gt.0)then
      write(8,610)Amass,AZ,AN,rd(3)
  610 format(//15x,'Selfconsistent Calculation for A:',
     1 f6.1,' (Z:',f6.1,', N:',f6.1,' )'/
     2 28x,'Initial Densities and Fields   (rd=',f8.4,' (fm))'//
     1 10x,'r',8x,'rhob',8x,'rhos',7x,'sigma',8x,
     2 'delta',8x,'omega',8x,' rho ',8x,'M*(p)',8x,'M*(n)')
      else
      write(8,6110)Amass,AZ,AN,Xsi,rd(3)
 6110 format(//3x,'Selfconsistent Calculation for A:',
     1 f6.1,' (Z:',f6.1,', N:',f6.1,' )'/
     & 3x,'Asymmetry Z/A:',f6.3/
     2 3x,'Initial Densities and Fields (rd:',f8.4,' [fm]) on MFD-A.d')
      endif
c
      r=-dr+1.e-06
c      r=0.d0
      w=4.
      DO 20 n=1,nrx
      r=r+dr
      rst(n)=r
      w=6.-w
      ww=w
c
c Simpson Weights, including r^2:
c
      if(n.eq.1.or.n.eq.nrx)ww=1.
      rwgt(n)=r**2*ww*dr/3.
      do j=1,2
      f(j)=f(j)*df(j)
      rho_v(n,j)=rhoD(j)/(1.+f(j))
* Kf in Units of M:
      Xkf(n,j)=hbc*(fkf*rho_v(n,j))**(1./3.)/Xmass
      enddo
c
c Iso-scalar and Iso-vector Vector-densities:
c
      rho_v(n,3)=rho_v(n,1)+rho_v(n,2)
      rho_v(n,4)=rho_v(n,1)-rho_v(n,2)
c
      do m=1,IBX
      If(Jbo(m).eq.0.and.Ibo(m).eq. 0)B00=Bcpl(m)
      If(Jbo(m).eq.0.and.Ibo(m).eq. 1)B01=Bcpl(m)
      If(Jbo(m).eq.1.and.Ibo(m).eq. 0)B10=Bcpl(m)
      If(Jbo(m).eq.1.and.Ibo(m).eq. 1)B11=Bcpl(m)
      If(Jbo(m).eq.1.and.Ibo(m).eq.-1)C1M=Bcpl(m)
      enddo
c
      Zfp=Xkf(n,1)
      Zfn=Xkf(n,2)
      dvp=rho_v(n,1)*hbcm3
      dvn=rho_v(n,2)*hbcm3
*
*  nuclear matter effective mass
*  solution of the coupled transcendental equations for M*(p) and M*(n)
*
      CALL SCD(Zfp,Zfn,dvp,dvn,XMp,XMn,dsp,dsn)
*
*  Effective Masses and IsoScalar/IsoVector Scalar densities:
*
      Xstar(n,1)=XMp
      Xstar(n,2)=XMn
      rho_s(n,1)=dsp/hbcm3
      rho_s(n,2)=dsn/hbcm3
      rho_s(n,3)=rho_s(n,1)+rho_s(n,2)
      rho_s(n,4)=rho_s(n,1)-rho_s(n,2)
*
* Scalar Fields g*Phi/M (pure numbers!):
*
      s00=-0.5*(XMp+XMn-Y0(1)-Y0(2))
      s01=-0.5*(XMp-XMn-Y0(1)+Y0(2))
      Bfld(n,1)=s00
      Bfld(n,2)=s01
*
* Vector Fields g*V/M (pure numbers!):
*
      Bfld(n,3)=B10*hbcm3*rho_v(n,3)
      Bfld(n,4)=B11*hbcm3*rho_v(n,4)
C
C COULOMB (in units of M, pure numers!)
C
      IF(R.LT.RC)THEN
       UC=UCINT*(3.-(R/RC)**2)
      ELSE
       UC=UCEXT/R
      ENDIF
      bfld(n,5)=UC
c
c Fermi-energies (averaged over nuclear volume):
c
      Ep=Bfld(n,3)+Bfld(n,4)+Bfld(n,5)+sqrt(Xkf(n,1)**2+XMp**2)
      En=Bfld(n,3)-Bfld(n,4)          +sqrt(Xkf(n,2)**2+XMn**2)
      Ef(1)=Ef(1)+Rwgt(n)*rho_v(n,1)*Ep
      Ef(2)=Ef(2)+Rwgt(n)*rho_v(n,2)*En
c
      If(IPR.gt.0)then
      write(8,611)r,rho_v(n,3),rho_s(n,3),(bfld(n,k),k=1,4),
     & (xstar(n,k),k=1,2)
  611 format(3x,f8.3,8e13.5)
      endif
c
      write(23,2000)r,rho_v(n,3),rho_s(n,3),(Bfld(n,k),k=1,4),
     & (Xstar(n,k),k=1,2)
 2000 format(f8.2,10e14.6)
c
c Densities in units of 1/M^3, i.e. pure numbers!:
c
      do 20 j=1,4
      rho_v(n,j)=rho_v(n,j)*hbcm3
      rho_s(n,j)=rho_s(n,j)*hbcm3
   20 CONTINUE
c
      close(23)
c
      Ef(1)=fpi*Ef(1)*Xmass/AZ
      Ef(2)=fpi*Ef(2)*Xmass/AN
      Sp   =Ef(1)-Xmss(1)
      Sn   =Ef(2)-Xmss(2)
      Efm  =(AZ*Ef(1)+AN*Ef(2))/Amass
      Sm   =(AZ*Sp   +AN*Sn   )/Amass
      write(8,6120)Ef(1),Sp,Ef(2),Sn,Efm,Sm,Ef(1)/Xmass,Ef(2)/Xmass
 6120 format(3x,'Initial Fermi Energies and Separation Energies:'/
     &       3x,'Ef(p)  :',f10.4,'  S(p)   :',f10.4,' [MeV]'/
     &       3x,'Ef(n)  :',f10.4,'  S(n)   :',f10.4,' [MeV]'/
     &       3x,'Ef(A)  :',f10.4,'  S(A)   :',f10.4,' [MeV]'/
     &       3x,'Ef(p)  :',f10.6,'  Ef(n)  :',f10.6,' [1/M]')
c
c TF densities etc.:
c
      open(23,file='dtf.d',status='UNKNOWN')
c
      Xp=0.d0
      Xn=0.d0
      Dfac=1.d0/(3*pisq*hbcm**3)
      do n=1,nrx
      Vp=bfld(n,3)+bfld(n,4)+bfld(n,5)
      Vn=bfld(n,3)-bfld(n,4)
      Ep=sqrt(Xkf(n,1)**2+Xstar(n,1)**2)
      En=sqrt(Xkf(n,2)**2+Xstar(n,2)**2)
      P1=(Ef(1)/Xmass-Vp)**2+Xstar(n,1)**2
     &   -2*Xstar(n,1)*(Ef(1)/Xmass-Vp)*rho_s(n,1)/rho_v(n,1)
      P2=(Ef(2)/Xmass-Vn)**2+Xstar(n,2)**2
     &   -2*Xstar(n,2)*(Ef(2)/Xmass-Vn)*rho_s(n,2)/rho_v(n,2)
      Q1=(Ef(1)/Xmass-Vp)**2-Xstar(n,1)**2
      Q2=(Ef(2)/Xmass-Vn)**2-Xstar(n,2)**2
      Wp=Xkf(n,1)**2-Q1
      Wn=Xkf(n,2)**2-Q2
      Q1=max(0.,Q1)
      Q2=max(0.,Q2)
      dtfp=Q1*sqrt(Q1)*Dfac
      dtfn=Q2*sqrt(Q2)*Dfac
      Xp=Xp+Rwgt(n)*dtfp
      Xn=Xn+Rwgt(n)*dtfn
      write(23,2000)rst(n),dtfp,dtfn,Vp,Vn,Wp,Wn,P1/hbcm**2,P2/hbcm**2
c      write(23,2000)rst(n),dtfp,dtfn,Vp,Vn,Xstar(n,1),Xstar(n,2)
      enddo
      Xp=fpi*Xp
      Xn=fpi*Xn
      write(*,*)' Xp:',Xp,' Xn:',Xn
      write(8,6121)Xp,Xn
 6121 format(3x,'ATF(p) :',f10.4,'  ATF(n) :',f10.4,' [Num]')
c
      close(23)
c
c      call VchD(rho_v,Xkf,Bfld,Vch,Upot)
c
      write(8,6122)(Vch(i)*Xmass,i=1,2),((Vch(i)-1)*Xmass,i=1,2)
 6122 format(/'   Chemical Potentials and Separation Energies:'/
     &        '   lambda(p):',f12.4,'   lambda(n):',f12.4,' [MeV]'/
     &        '   S(p)     :',f12.4,'   S(n)     :',f12.4,' [MeV]')
      open(23,file='rhoD.d',status='UNKNOWN')
      do n=1,NRX
      write(23,2000)rst(n),rho_v(n,1)/hbcm3,rho_v(n,2)/hbcm3
     & ,Upot(n,1),Upot(n,2)
      enddo
      close(23)
c
      return
      end
c---------------------------------------------------------------------
      REAL FUNCTION scalar(xkfp,xkfn,Xmp,Xmn,KB)
*
* Source term for the scalar meson field
* M* and kf in units of M !
*
c      implicit double precision (a-h,o-z)
      Parameter (KBX=6)
      common/pi4pi/pi,fpi,pisq
      common/const/Xmss(4),Xmass,alfa,hbc,hbcm,hbcm3
      common/MRatio/RM(2)
      common/boson/bmass(KBX),Jbo(KBX),Ibo(KBX),ISV(4)
      common/couple/G4pi(KBX),bcpl(KBX),rnge(KBX),ibx
      common/bnonl/NLINT(KBX),NLPOW(4,KBX),ANL(4,KBX),XNL(4,KBX)
c
      f(a,b)=sqrt(a**2+b**2)
      g(a,b)=a*(f(a,b)*b-a**2*log((b+f(a,b))/a))/(2*pisq)
c
      If(kb.eq.1)Iso= 1
      If(kb.eq.2)Iso=-1
c
      w=0.
      if(NLINT(KB).eq.0)goto 15
      Yp=RM(1)
      Yn=RM(2)
      phi=-0.5*(Xmp-Yp+Iso*(Xmn-Yn))
c
c  nonlinear scalar self-interaction
c
      DO 10 i=1,NLINT(KB)
      npow=NLPOW(i,KB)-1
   10 w=w+ANL(i,KB)*phi**npow
c
   15 scalar=g(Xmp,xkfp)+Iso*g(Xmn,xkfn)-w
c
      return
      end
c-----------------------------------------------------------------------
      SUBROUTINE Sfield(xkf,rho,phi,xstar,mboson,nrx)
c
c Solve Scalar Field Equation by Green Function Method
c M* and kf in units of M !
c
c      implicit double precision (a-h,o-z)
      Parameter (Ipr=0,KBX=6,KRX=400)
      common/boson/bmass(KBX),Jbo(KBX),Ibo(KBX),ISV(4)
      common/couple/G4pi(KBX),bcpl(KBX),rnge(KBX),ibx
      common/bnonl/NLINT(KBX),NLPOW(4,KBX),ANL(4,KBX),XNL(4,KBX)
      common/radius/rmax,dr,rst(KRX),rwgt(KRX),nrxx
      common/green/b(KRX,6),h(KRX,6),gnn(KRX,6)
      common/MRatio/RM(2)
      dimension s(KRX),t(KRX)
      dimension rho(nrx),phi(KRX,KBX),Xkf(KRX,2),Xstar(KRX,2)
      DATA nonx,eps/50,1.e-06/
c
      mb=mboson
c
      If(Bcpl(mb).eq.0.d0)return
c
      nonl=NLINT(mboson)
      If(ipr.gt.0)write(8,600)
  600 format(//,5x,'ITERATION OF THE NONLINEAR SIGMA FIELD EQ.'//
     1 4x,'iter.',9x,'dvol',9x,'drms',10x,'vol',10x,'rms')
c
      If(MB.eq.1)Iso= 1
      If(MB.eq.2)Iso=-1
c
      non=0
   10 non=non+1
      if(non.gt.nonx)goto 9000
c
c  compute the source term for non-linear coupling
c
      DO 15 n=1,nrx
      s(n)=scalar(xkf(n,1),xkf(n,2),xstar(n,1),xstar(n,2),mb)
   15 CONTINUE
c
      CALL Mfield(s,t,mboson,nrx)
c
      v1 =0.
      v2 =0.
      u10=0.
      u20=0.
      DO 35 n=1,nrx
      Xm1=Xstar(n,1)
      Xm2=Xstar(n,2)
      Xf1=Xkf(n,1)
      Xf2=Xkf(n,2)
      snn=gnn(n,mb)*drhos(Xf1,Xf2,Xm1,Xm2,phi(n,mb),mb)
      d=phi(n,mb)-t(n)
      v1=v1+(rst(n)*d)**2
      v2=v2+(rst(n)**2*d)**2
      u10=u10+(rst(n)*phi(n,mb))**2
      u20=u20+(rst(n)**2*phi(n,mb))**2
      phi(n,mb)=phi(n,mb)-d/(1.-snn)
      xstar(n,1)=RM(1)-phi(n,1)-phi(n,2)
      xstar(n,2)=RM(2)-phi(n,1)+phi(n,2)
c      write(*,629)rst(n),xkf(n,1),xkf(n,2),xstar(n,1),xstar(n,2),
c     # phi(n,mb),phi(n,3-mb)
 629  format(f8.3,8e12.4)
c      If(n.eq.20*(N/20))pause
      rho(n)=frhos(xkf(n,1),xstar(n,1))+Iso*frhos(xkf(n,2),xstar(n,2))
   35 CONTINUE
      v1=v1*dr
      v2=v2*dr
      u10=u10*dr
      u20=u20*dr
      d1=abs(v1/u10)
      d2=abs(v2/u20)
      If(ipr.gt.0)write(8,601)non,d1,d2,u10,u20
  601 format(5x,i4,4e13.5)
      if(d1.lt.eps.and.d2.lt.eps)goto 100
      goto 10
  100 return
 9000 write(8,690)non
  690 format(2x,'+++++++++++++  no convergence in',i4,' iterations')
      return
      end
c---------------------------------------------------------------------
      REAL FUNCTION drhos(zk1,zk2,Zm1,Zm2,Phi,mb)
c
c Derivative of isoscalar/isovector rho_s w.r.t. the scalar fields Phi_i
c NOTE: M*(p)=M*(M(p)/M-(x0+x1)) ; M*(n)=M*(M(n)/M-(x0-x1))
c       rho_s(0):=rho_s(p)+rho_s(n)    ;  rho_s(1):=rho_s(p)-rho_s(n)
c       M := average nucleon mass
c       x0:=g_0 Phi_0/M (sigma-field, mb=1)  ;  x1:=g_1 Phi_1/M (delta-field, mb=2)
c ALSO: M*(q) and kf(q) in units of M !
c
c      implicit double precision (a-h,o-z)
      Parameter (KBX=6)
      common/pi4pi/pi,fpi,pisq
      common/MRatio/RM(2)
      common/bnonl/NLINT(KBX),NLPOW(4,KBX),ANL(4,KBX),XNL(4,KBX)
      dimension xk(2),Xm(2)
c
      xk(1)=Zk1
      xk(2)=Zk2
      Xm(1)=Zm1
      Xm(2)=Zm2
c
      drhos=0.
      do 10 i=1,2
      xk2=Xk(i)*Xk(i)
      xm2=Xm(i)*Xm(i)
      ef=sqrt(xk2+xm2)
      drhos=drhos
     1-(xk(i)*ef+xm2*(2*(xk(i)/ef)-3*log((xk(i)+ef)/Xm(i))))
   10 continue
      drhos=drhos*0.5/pisq
      if(NLINT(MB).eq.0)return
c
c Non-Linear Self-Interaction:
c
      DO 20 n=1,NLINT(MB)
      npw=NLPOW(n,MB)-1
      drhos=drhos-ANL(n,MB)*npw*phi**(npw-1)
   20 CONTINUE
      return
c
      end
c----------------------------------------------------------------------
c$$$     commented out by TG
c$$$      double precision function Timer()
c$$$      implicit double precision (a-h,o-z)
c$$$c
c$$$c Timer Function for WIN only!
c$$$c
c$$$      Integer*2 Ihour,Imin,Isec,Ihund
c$$$c
c$$$      call gettim(Ihour,Imin,Isec,Ihund)
c$$$c
c$$$      Timer=3600.d0*Ihour+60.d0*Imin+DFloat(Isec)+DFloat(Ihund)*0.01d0
c$$$c
c$$$      return
c$$$      end
c----------------------------------------------------------------------
c$$$      SUBROUTINE SHOWTIME
c$$$      INTEGER*2 hour, minute, second, hund, thour
c$$$c
c$$$c SystemTime for WIN only!
c$$$c
c$$$      CHARACTER*1 ap
c$$$c
c$$$      call gettim(hour,minute,second,hund)
c$$$c
c$$$
c$$$      IF (hour .GT. 12) THEN
c$$$         ap = 'p'
c$$$         thour = hour - 12
c$$$      ELSE
c$$$         thour = hour
c$$$         ap = 'a'
c$$$      ENDIF
c$$$
c$$$      WRITE (*, 9002) hour, minute, second, hund
c$$$
c$$$9002  FORMAT (/' --> Time ',I2, ':', I2.2, ':', I2.2, ':', I2.2,
c$$$     #  ' ', A, 'm'/)
c$$$
c$$$      END
c----------------------------------------------------------------------
ctg      subroutine MKD(wkdir)
ctg      CHARACTER*255           wkdir
ctg      LOGICAL*4               result,MAKEDIRQQ,CHANGEDIRQQ
ctgc
ctgc Change/Make directory
ctgc
ctg      result = CHANGEDIRQQ(wkdir)
ctg      if(result)then
ctg         WRITE (*,*) 'Successfully changed to directory'
ctg         write (*,6000)wkdir
ctg      ELSE
ctgC
ctgC     Make a new directory
ctgC
ctg         result = MAKEDIRQQ (wkdir)
ctg         IF (result) THEN
ctg            WRITE (*,*) 'Successfully created new subdirectory'
ctg            WRITE (*,6100)wkdir
ctg         endif
ctg 6000    format(' ',a122)
ctg 6100    format(' --->',a122)
ctgC
ctgC        Change to new directory
ctgC
ctg         result = CHANGEDIRQQ(wkdir)
ctg         IF (result) THEN
ctg            WRITE (*,*) 'Successfully changed to new directory'
ctg         else
ctg            WRITE(*,*) ' Could not change to new directory'
ctg            stop
ctg         endif
ctgC
ctg      endif
ctgc
ctg      return
ctg      end
c---------------------------------------------------------------------
      SUBROUTINE EoS(Xsi)
c      implicit double precision (a-h,o-z)
      Parameter (KDX=500,KBX=6)
c
c EoS for infinite matter at asymmetry Xsi=Z/A
c
      common/cntrl/kout(10),epst
      common/pi4pi/pi,fpi,pisq
      common/boson/bmass(KBX),Jbo(KBX),Ibo(KBX),ISV(4)
      common/couple/G4pi(KBX),bcpl(KBX),rnge(KBX),ibx
      common/bnonl/NLINT(KBX),NLPOW(4,KBX),ANL(4,KBX),XNL(4,KBX)
      common/const/Xmss(4),Xmass,alfa,hbc,hbcm,hbcm3
c
      dimension EoA(KDX),Eint(KDX,4)
      Dimension rho_v(KDX,4),Xkf(KDX,2)
      Dimension rho_s(KDX,4),Phi_s(KDX,2),Xstar(KDX,2)
      Dimension Umf(KDX,4),Dens(4),Y0(2)
      Dimension Press(KDX),Compr(KDX)
c      Dimension Phi_v(KDX,2),Umf(KDX,4),Dens(4)
c
      Data Dnm/0.159/
c
      write(20,2000)Xsi
 2000 format('#     EoS for Z/A :',f4.2/
     &       '#',8x,'rho_v',9x,'E/A-M',9x,'ESelf',11x,'E00',
     &       11x,'E01',11x,'E10',11x,'E11')
      write(21,2100)Xsi
 2100 format('#     rho_s and M* for Z/A :',f4.2/
     &       '#',8x,'rho_v',9x,'M*(p)',9x,'M*(n)',5x,'ms*(Pole)',
     &        5x,'ms*(mech)',6x,'rho_s(p)',6x,'rho_s(n)')
      write(22,2200)Xsi
 2200 format('#     Tk  and Eint for Z/A :',f4.2/
     &       '#',8x,'rho_v',9x,'Tk(p)',9x,'Tk(n)',9x,'   Tk',
     &        8x,'EintBM',9x,'DMass',4x,'ENL(sigma)')
      write(23,2300)Xsi
 2300 format('#     Mean-Fields  for Z/A :',f4.2/
     &       '#',8x,'rho_v',9x,'sigma',9x,'delta',9x,'omega',
     &        11x,'rho')
c
      Dfac=3*pisq
      pwr=1.d0/3
      Y0(1)=Xmss(1)/Xmass
      Y0(2)=Xmss(2)/Xmass
      Dmax=6*Dnm
      delD=0.005
      NDX=Dmax/delD
      NDX=min(KDX,NDX+1)
      delD=Dmax/(NDX-1)
c
c density:
c
      rnm=-delD+1.d-06
      do 1000 nd=1,NDX
      rnm=rnm+delD
      rho_v(nd,3)=rnm
      rho_v(nd,4)=rnm*(2*Xsi-1)
      rho_v(nd,1)=rnm*xsi
      rho_v(nd,2)=rnm*(1-xsi)
      If(rho_v(nd,1).gt.0.d0)then
      Xkf(nd,1)=(Dfac*rho_v(nd,1))**pwr
      Xkf(nd,2)=(Dfac*rho_v(nd,2))**pwr
      else
      Xkf(nd,1)=0.d0
      Xkf(nd,2)=(Dfac*rho_v(nd,2))**pwr
      endif
c
      Pkfp=Xkf(nd,1)*hbc
      Pkfn=Xkf(nd,2)*hbc
      Zfp=Pkfp/XMass
      Zfn=Pkfn/XMass
      dnmp=rho_v(nd,1)*hbcm3
      dnmn=rho_v(nd,2)*hbcm3
      dnms=dnmp+dnmn
      dnmv=dnmp-dnmn
c
c omega and rho Vector Fields in [MeV]:
c
      Umf(nd,3)=bcpl(3)*dnms*Xmass
      Umf(nd,4)=bcpl(4)*dnmv*Xmass
c
c
c Scalar Densities and M*:
c
      call SCD(Zfp,Zfn,dnmp,dnmn,XMp,XMn,dscp,dscn)
c
      Xstar(nd,1)=XMp
      Xstar(nd,2)=XMn
      rho_s(nd,1)=dscp
      rho_s(nd,2)=dscn
      rho_s(nd,3)=dscp+dscn
      rho_s(nd,4)=dscp-dscn
c
c Scalar Fields in [MeV]:
c s00 := 1-(M*(p)+M*(n))/2M := g*sigma/M !
c s01 :=  -(M*(p)-M*(n)+dM0(pn))/2M := g*delta/M !
c These procedures account properly the Non-Linearity terms!
c
      yp=Xstar(nd,1)
      yn=Xstar(nd,2)
      If(Xsi.gt.0.d0)then
c      s00=1-0.5d0*(yp+yn)
      s00=-0.5d0*(yp+yn-Y0(1)-Y0(2))
      s01=-0.5d0*(yp-yn-Y0(1)+Y0(2))
      else
      s01=-Bcpl(2)*rho_s(nd,2)
      s00=Y0(2)-(yn+s01)
      endif
      Phi_s(nd,1)=s00*Xmass
      Phi_s(nd,2)=s01*Xmass
      Umf(nd,1)=Phi_s(nd,1)
      Umf(nd,2)=Phi_s(nd,2)
c
c re-scale to 1/fm^3 :
c
      do k=1,4
      rho_s(nd,k)=rho_s(nd,k)/hbcm3
      enddo
c
      Dens(1)=rho_s(nd,3)
      Dens(2)=rho_s(nd,4)
      Dens(3)=rho_v(nd,3)
      Dens(4)=rho_v(nd,4)
c
c Interaction Energies:
c
      EBM=0.d0
      do 100 kb=1,4
      Eint(nd,kb)=Umf(nd,kb)*Dens(kb)
      EBM=EBM+Eint(nd,kb)
  100 continue
c
      do kb=1,4
      Eint(nd,kb)=Eint(nd,kb)/rho_v(nd,3)
      enddo
      EBM=EBM/rho_v(nd,3)
c
c Kinetic Energy per particle:
c
      Xfp=Zfp/Xstar(nd,1)
      Xfn=Zfn/Xstar(nd,2)
      dMS=(Xstar(nd,1)-1)*Xmss(1)*Xsi+(Xstar(nd,2)-1)*Xmss(2)*(1-Xsi)
      ESelf=0.5*EBM+dMS
      Tkp=Fk(Xfp)*Xstar(nd,1)*Xmss(1)*Xsi
      Tkn=Fk(Xfn)*Xstar(nd,2)*Xmss(2)*(1-Xsi)
      Tkin=Tkp+Tkn
c
c Energy density from Non-Linear Scalar Self-Interactions:
c
      theta=Umf(nd,1)/Xmass
      ENL=Xmass*Ephi(theta,1)/(hbcm3*rho_v(nd,3))
c
c Effective sigma "Pole" Mass and mechanical Mass:
c
      Xpsc=1.d0
      Xmsc=1.d0
      If(theta.ne.0.d0)then
      X2p=1+Bcpl(1)*d1Wphi(theta,1)/theta
      X2m=1+Bcpl(1)*d2Wphi(theta,1)
      Xpsc=sqrt(X2p)
      If(X2m.ge.0.d0)Xmsc= sqrt(X2m)
      If(X2m.lt.0.d0)Xmsc=-sqrt(abs(X2m))
      endif
      ESelf=ESelf-0.5*ENL
c
      EoA(nd)=Tkin+dMS+0.5d0*(EBM-ENL)
c
      write(20,2001)rho_v(nd,3),EoA(nd),2*ESelf,(Eint(nd,kb),kb=1,4)
      write(21,2001)rho_v(nd,3),(Xstar(nd,i),i=1,2),Xpsc,Xmsc,
     & (rho_s(nd,k),k=1,2)
      write(22,2001)rho_v(nd,3),Tkp,Tkn,Tkin,EBM,dMS,ENL
      write(23,2001)rho_v(nd,3),(Umf(nd,k),K=1,4)
c
 1000 continue
c
 2001 format(10e14.6)
      close(20)
      close(21)
      close(22)
      close(23)
c
c Equilibrium point:
c
      do nd=2,NDX-1
      d1EoA=0.5*(EoA(nd+1)-EoA(nd-1))/delD
      d2EoA=    (EoA(nd+1)+EoA(nd-1)-2*EoA(nd))/delD**2
      Press(nd)=rho_v(nd,3)**2*d1EoA
      Compr(nd)=9*rho_v(nd,3)**2*d2EoA
      enddo
c
      n0=0
      do nd=2,NDX-2
      If(Press(nd)*Press(nd+1).lt.0.d0)n0=nd
      enddo
c
      If(n0.gt.2)then
       If(Compr(n0).gt.0.d0)then
      d1Pr=0.5*(Press(n0+1)-Press(n0-1))/delD
      d2Pr=(Press(n0+1)+Press(n0-1)-2*Press(n0))/delD**2
      rho_eq=rho_v(n0,3)-Press(n0)/d1Pr
      Xf_eq=(1.5d0*pisq*rho_eq)**pwr
      dComp=0.5*(Compr(n0+1)-Compr(n0-1))/delD
      C_eq=Compr(n0)+(rho_eq-rho_v(n0,3))*dComp
      P_eq=0.5*(rho_eq-rho_v(n0,3))**2*d2Pr
      d1E=0.5*(EoA(n0+1)-EoA(n0-1))/delD
      d2E=(EoA(n0+1)+EoA(n0-1)-2*EoA(n0))/delD**2
      x=(rho_eq-rho_v(n0,3))
      E_eq=EoA(n0)+x*d1E+0.5*x*x*d2E
      write(8,6040)Xsi,rho_eq,Xf_eq,E_eq,C_eq,P_eq
       else
      write(8,6041)Xsi
       endif
      else
      write(8,6041)Xsi
      endif
c
 6040 format(/' Equilibrium Properties for Xsi:',f5.2/
     &       ' rho     :',f9.5,' [1/fm^3  ]'/
     &       ' kF      :',f9.5,' [1/fm    ]'/
     &       ' E/A     :',f9.3,' [MeV     ]'/
     &       ' Compress:',f9.2,' [MeV     ]'/
     &       ' Pressure:',f9.4,' [MeV/fm^3]')
 6041 format(/' Equilibrium not found for  Xsi:',f5.2)
c
      return
      end
c---------------------------------------------------------------------
      SUBROUTINE SCD(Zfp,Zfn,dvp,dvn,XMp,XMn,dsp,dsn)
c      implicit double precision (a-h,o-z)
      Parameter (KBX=6)
c
c Scalar densities for asymmetric infinite matter :
c dvp ; dvn are the p/n vector densities
c Zfp:=Kfp/M ; Zfn:=Kfn/M ; XMp:=M*(p)/M ; XMn:=M*(n)/M
c -> Zp:=Zfp/XMp=kfp/M*(p), Zn:=Zfn/XMn=kfn/M*(n)
c dsp and dsn are the scalar p/n densities (x (hbc/M)^3 -> pure numbers!)
c -> implies for scalar fields:= Phi/M
c
      common/pi4pi/pi,fpi,pisq
      common/boson/bmass(KBX),Jbo(KBX),Ibo(KBX),ISV(4)
      common/couple/G4pi(KBX),bcpl(KBX),rnge(KBX),ibx
      common/bnonl/NLINT(KBX),NLPOW(4,KBX),ANL(4,KBX),XNL(4,KBX)
      common/const/Xmss(4),Xmass,alfa,hbc,hbcm,hbcm3
      Dimension DPN(2),d1DPN(2,2),ds(-1:1,2)
      Dimension X0(2),XPN(-1:1,2)
      Data ITERX,dm,eps/100,1.d-03,1.d-08/
c
      cp=bcpl(1)+bcpl(2)
      cm=bcpl(1)-bcpl(2)
      Xmass=Xmss(3)
      X0(1)=Xmss(1)/Xmass
      X0(2)=Xmss(2)/Xmass
      dm2=2*dm
c
c Starting Values (for scalar densities with M^* = M and then improving iteratively):
c ==> delta-Meson is allowed!
c
      dsp=Fs(zfp)*dvp
      dsn=Fs(zfn)*dvn
c
c Scalar mean-fields (in units of M -> pure numbers!):
c
      phi_s=bcpl(1)*(dsp+dsn)
      phi_v=bcpl(2)*(dsp-dsn)
      XM =1.d0-phi_s
      XM=max(0.1,XM)
      XPN(0,1)=XM-1.d0+X0(1)
      XPN(0,2)=XM-1.d0+X0(2)
c
      Iter=0
    1 Iter=Iter+1
      do k=-1,1
      ZM=XPN(0,1)+k*dm
      XPN(k,1)=ZM
      zp=Zfp/ZM
      ds(k,1)=dvp*Fs(zp)
      enddo
c
      do k=-1,1
      ZM=XPN(0,2)+k*dm
      XPN(k,2)=ZM
      zn=Zfn/ZM
      ds(k,2)=dvn*Fs(zn)
      enddo
c
c Non-Linear Self-Interactions:
c
      yp=XPN(0,1)
      yn=XPN(0,2)
      ypn=-0.5*(yp+yn-X0(1)-X0(2))
      WNL =      d1Wphi(ypn,1)*Bcpl(1)
      dWNL=0.5d0*d2Wphi(ypn,1)*Bcpl(1)
c
      do i=1,2
      DPN(i)=XPN(0,i)-X0(i)+ds(0,i)*cp+ds(0,3-i)*cm-WNL
      d1DPN(i,i  )=1.d0+cp*(ds(1,  i)-ds(-1,  i))/dm2+dWNL
      d1DPN(i,3-i)=     cm*(ds(1,3-i)-ds(-1,3-i))/dm2+dWNL
      enddo
c
c      write(36,629)iter,dMP,dMN,XPN(0,1),XPN(0,2),ds(0,1),ds(0,2)
c
      det=d1DPN(1,1)*d1DPN(2,2)-d1DPN(1,2)*d1DPN(2,1)
      dMP=-(DPN(1)*d1DPN(2,2)-DPN(2)*d1DPN(1,2))/det
      dMN=-(DPN(2)*d1DPN(1,1)-DPN(1)*d1DPN(2,1))/det
      XPN(0,1)=XPN(0,1)+dMP
      XPN(0,2)=XPN(0,2)+dMN
  629 format(i4,2e13.5,2x,2e13.5,2x,2e13.5)
      If((abs(dMP/XPN(0,1)).gt.eps.or.abs(dMN/XPN(0,2)).gt.eps).
     & and.Iter.le.ITERX)goto 1
      If(Iter.gt.ITERX)then
       write(8,*)' sub SCD: Iteration for M* failed!'
       stop
      endif
c
      XMp=XPN(0,1)
      XMn=XPN(0,2)
      dsp=ds (0,1)
      dsn=ds (0,2)
c
      return
      end
c---------------------------------------------------------------------
cTG      double precision Function Fs(z)
      Real Function Fs(z)
c      implicit double precision (a-h,o-z)
c
c Generic Scalar Density Profile Function, rho(vector) extracted!
c -> rho_s=rho_v*Fs(z) , z=kf/M*, for Fs(z) see below:
c
      z2=z*z
      z3=z2*z
      z4=z3*z
      z6=z4*z2
c
c Series expansions:
c
      If(z.lt.0.01d0)then
c
c Fs(z) for z -> 0:
c
c                         2         4         6      7
c              1 - 3/10 z  + 9/56 z  - 5/48 z  + O(z )
c
      Fs=(1-0.3d0*z2+9*z4/56-5*z6/48)
c
      goto 100
c
      endif
c
      If(z.gt.100.)then
c
c Fs(z) for z -> infinity:
c
c
c           3/4 - 3/2 ln(2) - 3/2 ln(z)         1     15  1     35   1        1
c 3/2 1/z + --------------------------- - 9/16 ---- + -- ---- - --- ---- + O(---)
c                        3                       5    64   7    256   9       11
c                       z                       z         z          z       z
c
      Fs=3.*(1.+(0.5-log(2.*z))/z2-3./(8.*z4)+5./(32.*z6))/(2.*z)
c
      goto 100
c
      endif
c
c Full expression for other z-values:
c
      Fs=3.*(z*sqrt(1.+z2)-log(z+sqrt(1.+z2)))/(2.*z3)
c
  100 continue
c
      return
      end
c---------------------------------------------------------------------
cTG      double precision Function Fk(z)
      real Function Fk(z)
c      implicit double precision (a-h,o-z)
c
c Generic Kinetc Energy Density Profile Function, rho(vector) extracted!
c -> <E-M*>=rho_v*Fk(z)*M^* , z=kf/M*, for full Fk(z) see below:
c
      z2=z*z
      z3=z2*z
      z4=z3*z
      z6=z4*z2
c
c Series expansions:
c
      If(z.lt.0.01)then
c
c Fk(z) for z -> 0:
c
c                       2         4         6      7
c                 3/10 z  - 3/56 z  + 1/48 z  + O(z )
c
cTG      Fk=0.3d0*z2-3*z4/56+z6/48
      Fk=0.3*z2-3.*z4/56.+z6/48.
c
      return
      endif
c
      If(z.gt.100.)then
c
c Fk(z) for z -> infinity:
c
c
c                       3/32 - 3/8 ln(2) - 3/8 ln(z)         1     15   1            1        1
c 3/4 z - 1 + 3/4 1/z + ---------------------------- - 3/32 ---- + --- ---- - 7/512 ---- + O(---)
c                                     3                       5    512   7            9       11
c                                    z                       z          z            z       z
c
      Fk=3.*(1.+1./z2+(0.125-log(2.*z)/2)/z4-1./(8.*z6)+5./(128.*z6*z2))
     & *z/4-1.
c
      return
      endif
c
c Full expression for other z-values:
c
      sqz=sqrt(1.+z2)
      Fk=3.*(z*sqz*(2.*z2+1.)-log(z+sqz))/(8.*z3)-1.
c
      return
      end
c---------------------------------------------------------------------
cTG      double precision Function Wphi(z,NB)
      Real Function Wphi(z,NB)
c      implicit double precision (a-h,o-z)
      Parameter (KBX=6)
c
c Non-Linear scalar Self-Interaction Potential:
c the TRANSFORMED coefficients XNL are used -> Wphi in MeV!
c -> W(z)=x_{n} z^n = x_3 z^3 + x_4 z^4 + .. ; z=g*Phi_s/M (pure number!)
c[ corresponds to theta/f_pi in Chiral Theory where M=g*f_pi ]
c
      common/boson/bmass(KBX),Jbo(KBX),Ibo(KBX),ISV(4)
      common/couple/G4pi(KBX),bcpl(KBX),rnge(KBX),ibx
      common/bnonl/NLINT(KBX),NLPOW(4,KBX),ANL(4,KBX),XNL(4,KBX)
c
      Wphi=0.
c
      If(NLINT(NB).gt.0)then
      do n=1,NLINT(NB)
      npw =NLPOW(n,NB)
      xn  =XNL(n,NB)
      Wphi=Wphi+xn*z**npw
      enddo
      endif
c
      return
      end
c---------------------------------------------------------------------
cTG      double precision Function d1Wphi(z,NB)
      Real Function d1Wphi(z,NB)
c      implicit double precision (a-h,o-z)
      Parameter (KBX=6)
c
c 1st Derivative of Non-Linear scalar Self-Interaction Potential:
c the UN-TRANSFORMED coefficients ANL are used!
c -> dW(z)=a_n z^(n-1 = a_3 z^2 + a_4 z^3 + .. ; z=g*Phi_s/M (pure number!)
c[ corresponds to theta/f_pi in Chiral Theory where M=g*f_pi ]
c
      common/boson/bmass(KBX),Jbo(KBX),Ibo(KBX),ISV(4)
      common/couple/G4pi(KBX),bcpl(KBX),rnge(KBX),ibx
      common/bnonl/NLINT(KBX),NLPOW(4,KBX),ANL(4,KBX),XNL(4,KBX)
c
      d1Wphi=0.
c
      If(NLINT(NB).gt.0)then
      do n=1,NLINT(NB)
      npw  =NLPOW(n,NB)-1
      xn   =ANL(n,NB)
      d1Wphi=d1Wphi+xn*z**npw
      enddo
      endif
c
      return
      end
c---------------------------------------------------------------------
cTG      double precision Function d2Wphi(z,NB)
      Real Function d2Wphi(z,NB)
c      implicit double precision (a-h,o-z)
      Parameter (KBX=6)
c
c 2nd Derivative of Non-Linear scalar Self-Interaction Potential:
c the UN-TRANSFORMED coefficients ANL are used!
c -> dW(z)=(n-1)a_n z^(n-2) = 2 a_3 z + 3 a_4 z^2 + .. ; z=g*Phi_s/M (pure number!)
c[ corresponds to theta/f_pi in Chiral Theory where M=g*f_pi ]
c NOTE: z=1-(xp+xn)/2  ==> additional factor -1/2 when evaluated w.r. to xp or xn in sub SCD
c x_tau=M*(tau)/M
c
      common/boson/bmass(KBX),Jbo(KBX),Ibo(KBX),ISV(4)
      common/couple/G4pi(KBX),bcpl(KBX),rnge(KBX),ibx
      common/bnonl/NLINT(KBX),NLPOW(4,KBX),ANL(4,KBX),XNL(4,KBX)
c
      d2Wphi=0.
c
      If(NLINT(NB).gt.0)then
      do n=1,NLINT(NB)
      npw  =NLPOW(n,NB)-2
      xn   =ANL(n,NB)*(npw+1)
      d2Wphi=d2Wphi+xn*z**npw
      enddo
      endif
c
      return
      end
c---------------------------------------------------------------------
cTG      double precision Function Ephi(z,NB)
      Real Function Ephi(z,NB)
c      implicit double precision (a-h,o-z)
      Parameter (KBX=6)
c
c Contribution of Non-Linear scalar Self-Interaction to Energy Density:
c
c Ephi(z)=z*d1Wphi(z)-2*Wphi(z):
c -> Ephi(z)=(1/3) a_3 z^3 + (1/2) a_4 z^4 + .. ; weights: 1-2/NLPOW
c
c the UN-TRANSFORMED coefficients ANL are used!
c -> z=g*Phi_s/M is a pure number!
c[ corresponds to theta/f_pi in Chiral Theory where M=g*f_pi ]
c
      common/boson/bmass(KBX),Jbo(KBX),Ibo(KBX),ISV(4)
      common/couple/G4pi(KBX),bcpl(KBX),rnge(KBX),ibx
      common/bnonl/NLINT(KBX),NLPOW(4,KBX),ANL(4,KBX),XNL(4,KBX)
c
      Ephi=0.
c
      If(NLINT(NB).gt.0)then
      do n=1,NLINT(NB)
      npw =NLPOW(n,NB)
      xn  =ANL(n,NB)*(1-2.d0/npw)
      Ephi=Ephi+xn*z**npw
      enddo
      endif
c
      return
      end
