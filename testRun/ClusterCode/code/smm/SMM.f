CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     THIS CODE WAS SENT BY ALEXANDRE BOTVINA WITH THE FOLLOWING MESSAGE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

cref. code: giv_smmgan_gsi_021.for
c
cDear Theo,
c
c
cI have attached the SMM code (in fortran), which you can use for secondary
cdeexcitation also, since it includes (automatically) evaporation, fission
cof heavy nuclei, fermi-break of light nuclei. This code is rather fully
cdescribed in J.P.Bondorf et al., Phys. Rep. , v.257, 133 (1995).
cYou can use also A.S.Botvina et al., Nucl. Phys. A475: 663-686 (1987).
c
cI recommend to keep the same parameters which are now in the code. You should
cgive in the beginning only:
c
c# Mass number and charge of the excited nucleus : IA0,IZ0
c# ETOT is excitation energy in GeV/nucleon
c# TPR0 is initial kinetic energy of nucleus (z-axis is assumed as flight
c  direction). May be easy for you to take TPR0=0, i.e., calculate the
c  deexcitation in the center of mass of the nucleus, and then recalculate
c  all momenta into the system you like
c# number of Monte-Carlo deexcitation events : ITNUM
c
cIf you want change other parameters - it is better to contact me.
c
c
c
c In the MAINMA subroutine:  ZVEZD1 is the main subroutine for Monte-Carlo
c simulation of disintegration of an excited nucleus.
c
c Its input: EE - excitation energy (GeV), AA - mass number, ZZ - charge,
c PX,PY,PZ - momenta of the nucleus (GeV/c), KSTART - start number of
c SPTU array.
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c Its output:  (KSTART-1) is total number of produced particles, they
c contained in SPT(10,500) array.
c SPT(4-6,n) - cos(THETA), sin(FI) and cos(FI) of
c the angles THETA and FI giving direction of particles' movement.
c SPT(7-9,n) - particles' kinetic energy (GeV), charge and mass (GeV).
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c If you want to use evaporation only, then you can uncomment EVANUC
c which stays just before ZVEZD1 in MAINMA (the input/output in EVANUC is
c the same) and comment ZVEZD1.
c
c
cFrom this output (SPT array) you can extract info about all particles
cproduced as a result of the deexcitation, event by event. Please,
clook how it is done in the code in the MAINMA. The code should work
cimmediately without extra subroutines.
c
cBest regards, Alexander

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     END OF THE ORIGINAL MESSAGE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     STARTING OF THE OTIGINAL CODE
C     CHANGES ARE LABELED AS "T.G" or "GAITANOS"
c     PLEASE, DO NOT REMOVE THE ADDITIONAL COMMENTS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c --- Feb.2005 GSI: change ZVEZD1- includ.small number of neutrons (<=2) in
c propagation with charged part., i.e. IP0=IP. Slightly influence spectra!
c --- Apr.2004 GSI: given to Lee Sobotka. To V.Uzhinsky - Dubna also?
c -------------- adapted for LINUX in GSI. Sept 2002. --------------
c *** 24 Jan 2000. GIV_SMMGAN  made from SMM_GAN_5 for Saclay. ***
c SMM_GAN.for made from smm_nat.for, version of July 98
c this version of SMM prepared for NATO project. Feb.98, Bologna.
c Made from the version of SMMAL7.for (last correct. Jan.98, NBI)
c Without preequilibrium parts but with flow and rotation.
c (The version ;9 of 5 March was used for finding connection between proton
c and neutron energies.)
c July 98, implemented in GANIL for checking energy ditribution depndence
c of c (Moretto) parameter and IMF correlation functions.
c -- 31 July, 98. GANIL: Change in PLACE to make disk consistent with
c rotation in Y-Z plane (compression along X not Z axis). See CELX.
c -- March 99, GANIL: new geometry - rotated sigar, see CELY,CEL
c -- March 99, new coordinate generation (PLACM) depending on targ./proj. --
c -- April 99 (Bo), introducing non-uniform flow (FLZ-param. along Z-axis --
c -- + add mass asymetry of Z-coord. of max fragment (POS0). --
c  ****** SMM_GAN_5 (made from smm_gan_4) ******
c --December 99 (Ganil). Modernizations: 1) DISNE2 + DISN02 (with FLZ) instead
cof DISNET + DISNE0; 2) new SELECA (SELECA_P - poissons); 3) (IP-1) instead of
c IP in FINDT.




CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     STATISTICAL MULTIFRAGMENTATION CODE (called event-by-event)
C     - changes are labeled by T.G. or GAITANOS (DO NOT REMOVE THESE COMMENTS)
C
C     INPUT:
C     ia0,iz0 (mass, charge of fragmenting source)
C     etot (excitation energy of fragmenting source in GeV/ia0)
C     OUTPUT:
C     Particles (number of final states (produced particles) per ensemple)
C     TheirCarges,TheirMasses,Momentum (charge,mass and 4-momentum of final states)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      SUBROUTINE MAINMA(ISOURCE,iseedi,IA0,IZ0,ETOT,
     &     SubSequentRuns,eflow_input,iflow_input, SMM_Flag,
     &     Particles,TheirMasses,TheirCharges,TheirMechanism,
     &     Ort,momentum)


C MULTIFRAGMENTION of moving source (NBI-95), with preequil.particles
c -- + ZVEZD2 -- FOR ALADIN experiment. (see JPRE and also E11(A)).

      COMMON /THEOKSTART/KSTART

c      common /randomnumber/iseedi

      COMMON /BLOKC/SPT(10,500) /BLFINL/A0F,Z0F,E0F,TF,FL,SRF,WW3F
      COMMON /BLMULF/IMULF /BLOKFR/FRA(500),SPSR,S3SR,Z3SR,MFRAG
      COMMON /BLOKF1/EXFRA(500),EKFRA(500)
      COMMON /BLOKF4/SAAM,SEKAM,SEXAM,SSEKX,SSEX
      COMMON /BLPARA/EPSIL0,FKAPPA,FKACOL
      COMMON /BLFLOW/IFLOW,AAA,EFLOW /BLAMOM/IAMOM,ANMOM,EMOM
      COMMON /BLTMIC/TMICR,TMICR2,NMICR,TMICT
      COMMON /BDMR0/IAAA(100),IZZZ(100),NFRAGG
     *,PXXX(100),PYYY(100),PZZZ(100)
      COMMON /BLJJJ/JJJ
      DIMENSION HISTA(400),HISTZ(200),PNUCL(3),PP(3)
      DIMENSION MCOR(20),MCI3(20)
      DIMENSION SFZFR(200),MZFR(200),ZMZFR(200)
      DIMENSION SFZB3(200),ZMZB3(200),MZB3(200)
      DIMENSION APAT(500),ZPAT(500),EKIP(500),TETP(500),PHIP(500)
      DIMENSION EKZ(500),ETRZ(500),NUMZ(500),RPHI(26)
      DIMENSION TZB(200),TZB3(200),MZB(200),ANPH0(500)
      DIMENSION VPZT(46),VPZ20(46),VPZ30(46),VPZ40(46),VPZ50(46)
      REAL NUMBZ


      COMMON /BLCEN/ XC(500),YC(500),ZC(500)

c      COMMON /BLCEN_THEO/ XC_THEO(500),YC_THEO(500),ZC_THEO(500) !their positions

      common / acceptEvent/Error_5


      real momentum(0:3,500),momsq
      real Ort(1:3,500)

      logical Error_5

      integer Particles
      integer TheirCharges(500)
      integer TheirMasses(500)
      integer TheirMechanism(500) !1->fission- / 2->evaporation-like processes
      INTEGER SubSequentRuns
      integer SMM_Flag

C----------------------------------------------------------------------------

c   INPUT for calc. with flow and rotation:
c   existing or not add. flow and rotation IFLOW, IAMOM (1 or 0)
cT.G.      IFLOW=0
         iflow = iflow_input
      IAMOM=0
c     energy of flow and rotation: EFLOW, EMOM in GeV/N
cT.G.      EFLOW=0.000
      eflow = eflow_input
      EMOM=0.000
c  for nonuniform flow along Z-axis: FLZ times larger than to X,Y
      FLZ=1.
c   main INPUT:
c   kinetic energy of the source in GeV/N : TPR0
      TPR0=0.0
c   mass number and charge of the source : IA0,IZ0
c      IA0=207
c      IZ0=82
c   total excit. energy (GeV/N) per nucleon : ETOT
c      ETOT=0.004
c      ETOT=0.0005*IITT
c -- see changing or not of coulomb volume : V_c=V_0+V_f  --
c -- in ZVEZD1: FRAMIK,FRAGTE (,FIXT)
      FKACOL=2.
c parametrization of E11(i): for central collisions- IPERI=0 ,
c for peripherial collisions- IPERI=1
cTG      IPERI=1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (iflow.eq.0) then
         IPERI = 1
      else
         IPERI = 0
      endif


      ITNUM = 1 !NumEnsemples

      NEVENT=1
      DO 500 IT=1,NEVENT
      ANMOMS=0.
      MFRAG=0
      SPSR=0.
      S3SR=0.
      Z3SR=0.
      DO 600 IJIJ=1,500
      EXFRA(IJIJ)=0.
      EKFRA(IJIJ)=0.
 600  FRA(IJIJ)=0.
      SAAM=0.
      SEKAM=0.
      SEXAM=0.
      SSEKX=0.
      SSEX=0.
      CALL DMR00
      CALL BPOCH1
      CALL BSPEC1
      INUC=0
      NUMB=0
      NUMBFR=0
      ITWRI=0
      A00=IA0
      Z00=IZ0
      IMULF=1
      A0F=1.
      Z0F=1.
      E0F=1.
      TPR=TPR0*IA0
      ETOT0=ETOT*A00
      PPR=SQRT(TPR*(TPR+2.*(IA0*0.940+ETOT0)))
      V00=PPR/SQRT(PPR**2+(0.94*A00+ETOT0)**2)
c  **********************************
      !WRITE 100 131,IMULF,FKACOL,A00,Z00
 131  FORMAT(2X,' multifrag.: IMULF=',I2,
     *'  FKACOL=',F5.1,3X,'   A00,Z00 =',2F6.0)
      IZMAXF=20
      !WRITE 100 132,TPR,ETOT0,PPR,V00
 132  FORMAT(1X,'TPR=',F8.4,' tot.ex.- ETOT0=',
     *F7.3,3X,'tot.moment- PPR=',F8.3,' cms veloc.V00=',F6.4)
      !WRITE 100 133,IZMAXF
 133  FORMAT(2X,'definition of IMF: Z=3-',I2)
      CALL THELI0
      CALL BFRENI
CCCCCCCCCCCCCCCCCCCCTHEO GAITANOS
c      iseedi = abs(time())
c      write(*,*) ' §§§§§§§§§§§§§§§§ iseedi = ',iseedi
      CALL RNDM1(iseedi)
CCCCCCCCCCCCCCCCCCCCTHEO GAITANOS
c      CALL RNDM1(0) !<-- Original
      DO 120 I=1,400
 120  HISTA(I)=0.
      DO 121 I=1,200
      IF(I.LE.20) MCOR(I)=0
      IF(I.LE.20) MCI3(I)=0
      MZFR(I)=0.
      ZMZFR(I)=0.
      SFZFR(I)=0.
      MZB3(I)=0.
      ZMZB3(I)=0.
      SFZB3(I)=0.
      MZB(I)=0.
      TZB(I)=0.
      TZB3(I)=0.
 121  HISTZ(I)=0.
      ICFI=0
      AFIS=0.
      SZBOU=0.
      SZBOU3=0.
      SZMAX=0.
      SMIMF=0.
      SMIMF2=0.
      RMZ=0.
      SNC=0.
      SAMAX=0.
      EKMAX=0.
      SNMAX=0.
      EKINZ=0.
      EKINZ2=0.
      NUMBZ=0
c **** generating the thermal source (AN - mass number, ZN - charge,
c **** ENEXT - excitation energy in GeV)
      AN=A00
      ZN=Z00
      ENEXT=ETOT0
      EN0=ENEXT/AN
c      PZE=V00*(ENEXT+0.94*AN)/SQRT(1.-V00**2)
      PZE=PPR
      CALL DELAM(A00,Z00,DLM00,DSH00,BAR00)
      !WRITE 100 607,DLM00
 607  FORMAT(2X,'  DLM00=',F7.2)
c -----------------------------------------
      INUC=INUC+1
c      IF(INUC.EQ.1.OR.MOD(INUC,500).EQ.0.)
      !WRITE 100 400,ITNUM,AN,ZN,EN0
 400  FORMAT(2X,'ITNUM=',I5,
     *'  AN,ZN,EN0=',3F9.4)
!      IF(IFLOW.EQ.1.OR.IAMOM.EQ.1) !WRITE 100 401,IFLOW,IAMOM,EFLOW,EMOM
 401  FORMAT(2X,'IFLOW,IAMOM=',2I2,'  including in EN0 - EFLOW,EMOM=',
     *2F9.4)
      DO 139 I=1,500
      IF(I.LE.26) RPHI(I)=0.
      IF(I.LE.46) VPZ20(I)=0.
      IF(I.LE.46) VPZ30(I)=0.
      IF(I.LE.46) VPZ40(I)=0.
      IF(I.LE.46) VPZ50(I)=0.
      IF(I.LE.46) VPZT(I)=0.
      EKZ(I)=0.
      ETRZ(I)=0.
      NUMZ(I)=0
 139  CONTINUE
c *******  open file to write - VAX system ********
      JWRITE=0
      IF(JWRITE.EQ.0) GO TO 181
c      OPEN(51,NAME='FILE51',STATUS='UNKNOWN',FORM='UNFORMATTED')
!orig      OPEN(51,NAME='FILE95')
      OPEN(51,file='FILE95')   !TG
      !WRITE 100 1830
 1830 FORMAT(1X,'   open file for writing')
 181  CONTINUE
c **************************************
      DO 410 JJJ=1,ITNUM
c --- for checking mom. and energ. balance --------
      PTX=0.
      PTY=0.
      PTZ=0.
      ETM=0.
      ETK=0.
c -------------------------------------------------
      NUMB=NUMB+1
      IBIG=0
      NZ330=0
      IZBOU=0
      IZBOU3=0
      NPAT=0
      JMAX=0
      ZMAX=0.
      IMF=0
      DO 200 I=1,10
      DO 200 J=1,500
 200  SPT(I,J)=0.
      KSTART=1
      AA=AN
      IF(AA.LT.1.) GO TO 410
      ZZ=ZN
      EE=ENEXT
c *** reduction of energy by flow and rotation ***
      IF(IFLOW.EQ.1) EE=EE-EFLOW*AA
      IF(IAMOM.EQ.1) EE=EE-EMOM*AA
c ********************************
      PNUCL(1)=0.
      PNUCL(2)=0.
      PNUCL(3)=PZE
      PX=PNUCL(1)
      PY=PNUCL(2)
      PZ=PNUCL(3)
      IF(JJJ.GT.2.OR.INUC.GT.2) GO TO 157
      !WRITE 100 154,AA,ZZ,EE
 154  FORMAT(2X,'MULTIFRAG. AA,ZZ, EE(GeV)=',2F6.1,1X,F9.4)
!      IF(IFLOW.EQ.1.OR.IAMOM.EQ.1) !WRITE 100 159,IFLOW,IAMOM,EFLOW,EMOM
 159  FORMAT(2X,'IFLOW,IAMOM=',2I2,'  add.energ.(GeV/N) EFLOW,EMOM =',
     *2F9.4)
 157  CONTINUE
ccc      !WRITE 100 502,JJJ,EE,AA,ZZ,KSTART
ccc 502  FORMAT(2X,'ZVEZD1: JJJ=',I4,' EE,AA,ZZ,KSTART=',3F9.3,I5)

      if (SMM_Flag.eq.0) then
         CALL EVANUC(EE,AA,ZZ,PX,PY,PZ,KSTART) !only evaporation
      else
         CALL ZVEZD1(EE,AA,ZZ,PX,PY,PZ,KSTART) !full SMM code (includes everything!!!)
      endif
c ---------------------------------------
      ANMOMS=ANMOMS+ANMOM
      KST1=KSTART-1
      IF(KST1.LT.1) GO TO 410
!      IF(JJJ.LE.2.AND.INUC.LE.2) !WRITE 100 8
 8    FORMAT(1X,'    BREAK UP.')
      DO 105 J=1,KST1
      IF(IT.GT.2) GO TO 170
c      IF(JJJ.LE.2.AND.INUC.LE.2) PRINT 7,ISOURCE,JJJ,
c     &     (SPT(JPER,J),JPER=1,10)
c      IF(INUC.LE.2) PRINT 7,ISOURCE,JJJ,
c     &     (SPT(JPER,J),JPER=1,10)
 7    FORMAT(1X,'SPT=',2I5,1x,10F12.5)
 170  CONTINUE
      IZ=INT(SPT(8,J)+0.5)
      IA=INT(SPT(9,J)/0.94+0.5)
      IF(SPT(4,J).GT.1.) SPT(4,J)=1.
      IF(SPT(4,J).LT.-1.) SPT(4,J)=-1.
      IF(SPT(5,J).GT.1.) SPT(5,J)=1.
      IF(SPT(5,J).LT.-1.) SPT(5,J)=-1.
      IF(SPT(6,J).GT.1.) SPT(6,J)=1.
      IF(SPT(6,J).LT.-1.) SPT(6,J)=-1.
      COT=SPT(4,J)
      A=IA
      Z=IZ
      IFIS=INT(SPT(3,J)+0.1)
      IF(IFIS.EQ.1) AFIS=AFIS+A
      IF(IFIS.EQ.1) ICFI=ICFI+1
      IF(SPT(8,J).GT.ZMAX) JMAX=J
      IF(SPT(8,J).GT.ZMAX) ZMAX=SPT(8,J)
      IF(IA.LT.1.OR.IA.GT.400) GO TO 304
      HISTA(IA)=HISTA(IA)+1.
 304  CONTINUE
      IF(IZ.LT.0.OR.IZ.GT.199) GO TO 305
      HISTZ(IZ+1)=HISTZ(IZ+1)+1.
 305  CONTINUE
      IF(IZ.GE.3) IBIG=IBIG+1
      IF(IZ.GE.3.AND.IZ.LE.IZMAXF) NZ330=NZ330+1
      IF(IZ.GE.1) NPAT=NPAT+1
      IF(IZ.GE.3) IZBOU3=IZBOU3+IZ
      IF(IZ.GE.2) IZBOU=IZBOU+IZ
      IF(IZ.GE.3.AND.IZ.LE.20) IMF=IMF+1
 105  CONTINUE
      SNMAX=SNMAX+1.
      SAMAX=SAMAX+INT(SPT(9,JMAX)/0.94+0.5)
      EKMAX=EKMAX+1000.*SPT(7,JMAX)
      CALL THELI1(KSTART)
      NUMBFR=NUMBFR+1
      SZBOU3=SZBOU3+IZBOU3
      SZBOU=SZBOU+IZBOU
      SZMAX=SZMAX+ZMAX
      SMIMF=SMIMF+NZ330
      SMIMF2=SMIMF2+NZ330*NZ330
      IF(IZBOU.GT.0.) RMZ=RMZ+FLOAT(NZ330)/FLOAT(IZBOU)
      SNC=SNC+NPAT
      IF(IBIG.GT.0.AND.IBIG.LE.19) MCOR(IBIG+1)=MCOR(IBIG+1)+1
      IF(NZ330.GE.0.AND.NZ330.LE.19) MCI3(NZ330+1)=MCI3(NZ330+1)+1
      IF(IZBOU.GE.1.AND.IZBOU.LE.200)
     *MZB(IZBOU)=MZB(IZBOU)+1
      IF(IZBOU.GE.1.AND.IZBOU.LE.200)
     *TZB(IZBOU)=TZB(IZBOU)+TMICT
      IF(IZBOU3.GE.1.AND.IZBOU3.LE.200)
     *TZB3(IZBOU3)=TZB3(IZBOU3)+TMICT
      IF(IZBOU3.GE.1.AND.IZBOU3.LE.200)
     *MZB3(IZBOU3)=MZB3(IZBOU3)+1
      IF(IZBOU3.GE.1.AND.IZBOU3.LE.200)
     *ZMZB3(IZBOU3)=ZMZB3(IZBOU3)+ZMAX
      IF(IZBOU3.GE.1.AND.IZBOU3.LE.200)
     *SFZB3(IZBOU3)=SFZB3(IZBOU3)+NZ330
      IF(NPAT.GE.1.AND.NPAT.LE.200)
     *MZFR(NPAT)=MZFR(NPAT)+1
      IF(NPAT.GE.1.AND.NPAT.LE.200)
     *ZMZFR(NPAT)=ZMZFR(NPAT)+ZMAX
      IF(NPAT.GE.1.AND.NPAT.LE.200)
     *SFZFR(NPAT)=SFZFR(NPAT)+NZ330
      DO 135 J=1,KST1
      IZ=INT(SPT(8,J)+0.5)
      IA=INT(SPT(9,J)/0.94+0.5)
      IF(IA.LE.0) GO TO 135
      AIA=IA
      ZIZ=IZ
      CALL DELAM(AIA,ZIZ,DLIAZ,DSIAZ,BRIAZ)
      ETM=ETM+0.001*DLIAZ
      COT=SPT(4,J)
      SIT=SQRT(1.-COT*COT)
      CALL TINP(PP,COT,SIT,SPT(6,J),SPT(5,J),SPT(7,J),SPT(9,J))
      PTX=PTX+PP(1)
      PTY=PTY+PP(2)
      PTZ=PTZ+PP(3)
      ETK=ETK+SPT(7,J)
c ------ checking isotropical fly-out for max fragm. -------
      IF(J.NE.JMAX) GO TO 1350
      VPZ=30.*(PP(3)/SPT(9,J))
      CALL HIST1(VPZ,-5.,5.,0.25,VPZT,46,1.)
      IF(ZMAX.GE.15.AND.ZMAX.LT.25)
     *CALL HIST1(VPZ,-5.,5.,0.25,VPZ20,46,1.)
      IF(ZMAX.GE.25.AND.ZMAX.LT.35)
     *CALL HIST1(VPZ,-5.,5.,0.25,VPZ30,46,1.)
      IF(ZMAX.GE.35.AND.ZMAX.LT.45)
     *CALL HIST1(VPZ,-5.,5.,0.25,VPZ40,46,1.)
      IF(ZMAX.GE.45.)
     *CALL HIST1(VPZ,-5.,5.,0.25,VPZ50,46,1.)
 1350 CONTINUE
c ----------------------------------------------------------
 135  CONTINUE
!      IF(JJJ.LE.2) !WRITE 100 608,PTX,PTY,PTZ,ETM,ETK
 608  FORMAT(2X,'sum. balance:      PTX, PTY, PTZ=',3F9.3/
     *2X,' ETM , ETK=',2F10.4)
c -- writing the event --
c      !WRITE 100 914,IZ0,IA0,ETOT,KST1,NEUTR,IMF
c 914  FORMAT(2X,'IZ0,IA0=',2I5,' ETOT=',F9.3,
c     *' kst1=',I3,' NEUTR=',I3,' IMF=',I2)
      IF(JWRITE.EQ.0) GO TO 182
c ----  writing only events with large Z_max ---
c      IF(ZMAX.LE.47.) GO TO 182
      ITWRI=ITWRI+1
c ----------------------------------------------
      NUMPAT=KST1
ccc      WRITE(51,136)NUMPAT
ccc 136  FORMAT(I3)
      ENEXW=ENEXT*1000.
c      WRITE(51,*)(NUMPAT,AN,ZN,ENEXW,FKACOL)
      NEUTR=KST1-NPAT
      WRITE(51,*)IZ0,IA0,ETOT,KST1,NEUTR,IMF
      PHIRAN=360.0*RNDM(-1)
      DO 134 J=1,KST1
      IZP=INT(SPT(8,J)+0.5)
      IF(IZP.LE.0.) GO TO 134
      IAP=INT(SPT(9,J)/0.94+0.5)
      COT=SPT(4,J)
      SIT=SQRT(1.-COT*COT)
      COTE=SPT(4,J)
      COFI=SPT(6,J)
      ANTET=(180./3.14159)*ACOS(COTE)
      ANPHI=(180./3.14159)*ACOS(COFI)
      IF(SPT(5,J).LT.0.) ANPHI=360.-ANPHI
      IF(ANPHI.GE.360.) ANPHI=ANPHI-360.
      ANPH0(J)=ANPHI
c --- add. randomization on PHI (since rotation is only along X) -------
      IF(IAMOM.EQ.0) GO TO 912
      IF(IAMOM.GT.0) ANPHI=ANPHI+PHIRAN
      IF(IAMOM.GT.0.AND.ANPHI.GE.360.) ANPHI=ANPHI-360.
      AOFI=(3.14159/180.)*ANPHI
      SPT(6,J)=COS(AOFI)
      SPT(5,J)=SIN(AOFI)
 912  CONTINUE
c ----------------------------------------------------------------------
      CALL TINP(PP,COT,SIT,SPT(6,J),SPT(5,J),SPT(7,J),SPT(9,J))
      EKIN=1000.*SPT(7,J)
c -------------------------------------
ccc      WRITE(51,137)IAP,IZP,EKIN,ANTET,ANPHI
ccc 137  FORMAT(I4,I5,F8.1,F7.2,F7.2)
c      WRITE(51,*)(IAP,IZP,EKIN,ANTET,ANPHI)
      WRITE(51,*)IAP,IZP,PP(1),PP(2),PP(3)
 134  CONTINUE
cc      TYPE*,'   writing one event in file    NUMPAT=',NUMPAT
c      TYPE*,'   writing one event in file    NPAT=',NPAT
 182  CONTINUE
c -----------------------
      EKINZ0=0.
      NUMPAT=KST1
c      PHIRAN=360.0*RNDM(-1)
      NFRAGG=0
      DO 138 J=1,KST1
      IZP=INT(SPT(8,J)+0.5)
      IAP=INT(SPT(9,J)/0.94+0.5)
      IF(IZP.LE.0) GO TO 910
      NFRAGG=NFRAGG+1
      IAAA(NFRAGG)=IAP
      IZZZ(NFRAGG)=IZP
      COT=SPT(4,J)
      SIT=SQRT(1.-COT*COT)
      CALL TINP(PP,COT,SIT,SPT(6,J),SPT(5,J),SPT(7,J),SPT(9,J))
      PXXX(NFRAGG)=1000.*PP(1)
      PYYY(NFRAGG)=1000.*PP(2)
      PZZZ(NFRAGG)=1000.*PP(3)
 910  CONTINUE
      EKIN=1000.*SPT(7,J)
      COTE=SPT(4,J)
      SITE=1.-COTE*COTE
      COFI=SPT(6,J)
      ANTET=(180./3.14159)*ACOS(COTE)
      ANPHI=(180./3.14159)*ACOS(COFI)
      IF(SPT(5,J).LT.0.) ANPHI=360.-ANPHI
c      IF(JJJ.LE.2) !WRITE 100 911,IAP,IZP,SPT(4,J),SPT(5,J),SPT(6,J),
c     *ANPH0(J),ANPHI
c 911  FORMAT(2X,'IAP=',I3,' IZP=',I2,'  rot: spt(4-5-6)=',3F7.3,
c     *' FIold=',F6.1,' FInew=',F6.1)
c --- add. randomization on PHI (since rotation is only along X) -------
cc      IF(IAMOM.GT.0) ANPHI=ANPHI+PHIRAN
cc      IF(IAMOM.GT.0.AND.ANPHI.GE.360.) ANPHI=ANPHI-360.
c -------------------------------------
c new part. characteristics: IAP,IZP,EKIN,ANTET,ANPHI
      APAT(J)=IAP
      ZPAT(J)=IZP
      EKIP(J)=EKIN
      TETP(J)=ANTET
      PHIP(J)=ANPHI
      IF(IZP.GE.1) EKINZ0=EKINZ0+EKIN*SITE
      IF(IZP.GE.1) EKINZ=EKINZ+EKIN*SITE
      IF(IZP.GE.1) NUMBZ=NUMBZ+1
      IF(IZP.GE.1) EKZ(IZP)=EKZ(IZP)+EKIN
      IF(IZP.GE.1) ETRZ(IZP)=ETRZ(IZP)+EKIN*(1.-COTE*COTE)
      IF(IZP.GE.1) NUMZ(IZP)=NUMZ(IZP)+1
 138  CONTINUE
      EKINZ2=EKINZ2+EKINZ0*EKINZ0
      CALL DMR
      CALL CPOCH1(NUMPAT)
      CALL CSPEC1(NUMPAT)
      NUMPA1=NUMPAT-1
      DO 141 J1=1,NUMPA1
      J11=J1+1
      IAP1=APAT(J1)
      IZP1=ZPAT(J1)
      DO 142 J2=J11,NUMPAT
      IAP2=APAT(J2)
      IZP2=ZPAT(J2)
      IF(IAP1.EQ.4.AND.IZP1.EQ.2.AND.IAP2.EQ.4.AND.IZP2.EQ.2) GO TO 143
      GO TO 142
 143  CONTINUE
      DPHI=ABS(PHIP(J1)-PHIP(J2))
      IF(DPHI.GT.180.) DPHI=360.-DPHI
      CALL HIST1(DPHI,0.,180.,9.,RPHI,26,1.)
 142  CONTINUE
 141  CONTINUE


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      write(*,*) '//////////////////////////////'
c      write(*,*) jjj,kstart-1

      do ikk=1,kstart-1
         momentum(0,ikk) = SPT(7,ikk) + SPT(9,ikk) !total energy: E_kin + Mass
         momsq = sqrt(momentum(0,ikk)**2-SPT(9,ikk)**2) !modul of 3-momentum
         SinTheta = sqrt(1.-SPT(4,ikk)**2)              !sin(Theta)
         momentum(1,ikk) = momsq * SPT(6,ikk)* SinTheta
         momentum(2,ikk) = momsq * SPT(5,ikk)* SinTheta
         momentum(3,ikk) = momsq * SPT(4,ikk)
         momentum(0,ikk) = sqrt(SPT(9,ikk)**2+
     &        momentum(1,ikk)**2+
     &        momentum(2,ikk)**2+momentum(3,ikk)**2)

         TheirCharges(ikk) = int(SPT(8,ikk))
         TheirMasses(ikk)  = int(SPT(9,ikk)/0.940)

         if (SPT(1,ikk).eq.1.) then
            if (SPT(3,ikk).eq.1.) TheirMechanism(ikk)=1 !fission-like fragments
            if (SPT(3,ikk).eq.0.) TheirMechanism(ikk)=2 !evaporation-like fragments
         else
            TheirMechanism(ikk) = 0
         endif

      end do

!      call PLACE3_TG(iflow,ia0,(kstart-1),TheirMasses)

      Particles = kstart-1
      if (Particles.gt.500) then
         write(*,*) 'from SMM.f : dimension overflow!!!!',
     &        kstart,particles
         write(*,*) '!!! TERMINATION OF PROGRAM !!!'
         STOP
      endif

c      Source_Radius = 1.5*float(ia0)**0.333333333

      ! * store positions of produced particles
      ! * check that momenta of produced particles
      !   are determined correctly
      !   (Massive SPT(7,4-6,...) <--> momentum(0:3,...)
      do ikk=1,kstart-1
c         Ort(1,ikk) = 0.0 !XC_THEO(ikk)
c         Ort(2,ikk) = 0.0 !YC_THEO(ikk)
c         Ort(3,ikk) = 0.0 !ZC_THEO(ikk)
         Ort(1,ikk) = XC(ikk)
         Ort(2,ikk) = YC(ikk)
         Ort(3,ikk) = ZC(ikk)

         COT=SPT(4,ikk)
         SIT=SQRT(1.-COT*COT)
         CALL TINP(PP,COT,SIT,SPT(6,ikk),SPT(5,ikk),
     &             SPT(7,ikk),SPT(9,ikk))
         PXXX_TEST = PP(1)
         PYYY_TEST = PP(2)
         PZZZ_TEST = PP(3)
         EKIN_TEST = SPT(7,ikk)
         ETOT_TEST = SPT(7,ikk) + SPT(9,ikk)
         if ( (abs(ETOT_TEST-momentum(0,ikk)).gt.0.000001) .or.
     &        (abs(PXXX_TEST-momentum(1,ikk)).gt.0.000001) .or.
     &        (abs(PYYY_TEST-momentum(2,ikk)).gt.0.000001) .or.
     &        (abs(PZZZ_TEST-momentum(3,ikk)).gt.0.000001) ) then
            write(*,*) 'SMM-Code,routine MAINMA:'
            write(*,*) 'Wrong determination of 4-momenta'
            write(*,'(A,i5,2x,4f12.6)') '*** Index,momentum(0:3) = ',
     &           ikk,momentum(0,ikk),momentum(1,ikk),
     &           momentum(2,ikk),momentum(3,ikk)
            write(*,'(A,i5,2x,4f12.6)') '*** Index,P_TEST        = ',
     &           ikk,ETOT_TEST,PXXX_TEST,PYYY_TEST,PZZZ_TEST
            write(*,*) 'Termination of the program'
            STOP
         endif

c         write(*,1234) isource,TheirMasses(ikk),TheirCharges(ikk),
c     &        (momentum(ijj,ikk),ijj=0,3),
c     &        EKIN_TEST+SPT(9,ikk),PXXX_TEST,PYYY_TEST,PZZZ_TEST

c            write(*,1233) SubSequentRuns,isource,Particles,
c     &           SPT(1,ikk),
c     &           TheirMasses(ikk),
c     &           TheirCharges(ikk),
c     &           (momentum(ijj,ikk),ijj=0,3)

      end do

c 1233 format(3i5,5x,f8.4,5x,2i5,5x,4f8.4)
c 1234 format('SPT(1,...)=0: ',5i5,5x,4f8.4)
c 1235 format('SPT(1,...)=1: ',5i5,5x,4f8.4)
c 1236 format('SPT(1,...)=2: ',5i5,5x,4f8.4)
c 1237 format('SPT(1,...)=3: ',5i5,5x,4f8.4)
c 1234 FORMAT(1X,3i5,2x,'FRAGMENTS=',2i5,25F8.3)


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c
 410  CONTINUE
c ******* end of writing in file ******
      IF(JWRITE.EQ.0) GO TO 183
      CLOSE(51)
      !WRITE 100 1831,ITNUM,ITWRI
 1831 FORMAT(1X,' WRITE in file finished  ITNUM=',I7,' ITWRI=',I7)
 183  CONTINUE
c ***************************************
      !WRITE 100 5,INUC,NUMB,NUMBFR
 5    FORMAT(2X,' NUMBER OF NUCL. - INUC=',I6,' NUMB=',I6,
     *' NUMBFR=',I5)
c
      IF(NUMB.GT.0) EKINZ=EKINZ/FLOAT(NUMB)
      IF(NUMB.GT.0) EKINZ2=EKINZ2/FLOAT(NUMB)
      IF(NUMB.GT.0) NUMBZ=NUMBZ/FLOAT(NUMB)
      DDD=EKINZ2-EKINZ*EKINZ
      IF(DDD.LE.0.) EKINZ2=0.
      IF(DDD.GT.0.) EKINZ2=SQRT(DDD)
      !WRITE 100 750,NUMBZ,EKINZ,EKINZ2
 750  FORMAT(1X,' numb.of charg.part. in event NUMBZ=',F6.2/2X,
     *'their aver.transv.energ. EKINZ=',F7.3,'  and variance =',F7.3)
      IF(IT.GT.2) GO TO 171
      DO 602 IJIJ=1,IA0
      IF(FRA(IJIJ).GT.0.) EXFRA(IJIJ)=EXFRA(IJIJ)/FRA(IJIJ)
      IF(FRA(IJIJ).GT.0.) EKFRA(IJIJ)=EKFRA(IJIJ)/(IJIJ*FRA(IJIJ))
      IF(MFRAG.GT.0) FRA(IJIJ)=FRA(IJIJ)/FLOAT(MFRAG)
 602  CONTINUE
      IF(MFRAG.GT.0) SAAM=SAAM/FLOAT(MFRAG)
      IF(MFRAG.GT.0) SEKAM=SEKAM/(SAAM*FLOAT(MFRAG))
      IF(MFRAG.GT.0) SEXAM=SEXAM/FLOAT(MFRAG)
      IF(MFRAG.GT.0) SSEKX=SSEKX/FLOAT(IA0*MFRAG)
      IF(MFRAG.GT.0) SSEX=SSEX/FLOAT(IA0*MFRAG)
      IF(MFRAG.GT.0) SPSR=SPSR/FLOAT(MFRAG)
      IF(MFRAG.GT.0) S3SR=S3SR/FLOAT(MFRAG)
      IF(MFRAG.GT.0) Z3SR=Z3SR/FLOAT(MFRAG)
      !WRITE 100 601,MFRAG,SPSR,S3SR,Z3SR
 601  FORMAT(2x,'numb. of zvezd1-events: MFRAG=',I6/
     *2X,'avr. before second.de-excit.: SPSR=',F6.2,
     *' S3SR=',F6.2,' Z3SR=',F6.2)
      !WRITE 100 631,SAAM,SEXAM,SEKAM,SSEKX,SSEX
 631  FORMAT(2X,'hot max fragm: SAAM=',F6.2/
     *1X,' its exit(MeV/n) SEXAM=',F6.2,' ekin(MeV/n) SEKAM=',F6.2/
     *2X,'sum: kin.(coul)+exit.energ.of all prim.frag.(MeV/n) SSEKX=',
     *F6.2/2X,'sum of exit.energ.all prim.frag.(MeV/n) SSEX=',F6.2)
      !WRITE 100 623
 623  FORMAT(12X,'mass distribution of hot fragments')
      !WRITE 100 603,(FRA(IJIJ),IJIJ=1,IA0)
 603  FORMAT(3x,'FRA=',10F7.3)
      !WRITE 100 621
 621  FORMAT(12X,'excitation energies (MeV/N) of hot fragments')
      !WRITE 100 611,(EXFRA(IJIJ),IJIJ=1,IA0)
 611  FORMAT(1x,'EXFRA=',10F7.3)
      CALL BFRENE(IA0,MFRAG)
      !WRITE 100 622
 622  FORMAT(12X,'kinetic energies (MeV/N) of hot fragments')
      !WRITE 100 612,(EKFRA(IJIJ),IJIJ=1,IA0)
 612  FORMAT(1x,'EKFRA=',10F7.2)
 171  CONTINUE
c
      IF(ICFI.GT.0) AFFS=AFIS/FLOAT(ICFI)
      !WRITE 100 140,ICFI,AFFS
 140  FORMAT(2X,'NUMBER OF FIS.FRAGM.  ICFI=',I6,5X,
     *'MEAN FIS.FRAGM.  AFFS=',F8.3)
      !WRITE 100 312
 312  FORMAT(5X,'MASS  YIELD - A (FROM 1 TO IA0)')
      !WRITE 100 310,(HISTA(I),I=1,IA0)
c      CALL SORTNA(IA0,HISTA)
      !WRITE 100 311
 311  FORMAT(2X,'CHARGE  YIELD - Z (FROM 0 TO IZ0)')
      !WRITE 100 310,(HISTZ(I),I=1,IZ0)
 310  FORMAT(1X,10F8.0)
      !WRITE 100 314
 314  FORMAT(2X,'EVENTS WITH EMISSION OF  0,1,2,... HEAVY (Z.GE.3)
     *FRAGMENTS')
      !WRITE 100 315,(MCOR(I),I=1,20)
 315  FORMAT(1X,10I8)
      IF(SAMAX.GT.0.) EKMAX=EKMAX/SAMAX
      IF(SNMAX.GT.0.) SAMAX=SAMAX/SNMAX
      !WRITE 100 1020, SNMAX,SAMAX,EKMAX
 1020 FORMAT(2X,'numb.max.Z-fragm. SNMAX=',F7.0,' mass SAMAX=',F6.2,
     *' kin.energ.(MeV/u) EKMAX=',F7.2)
      !WRITE 100 397,IZMAXF
 397  FORMAT(2X,'DISTR.ON NUMB. OF FRAG. 3.LE.Z.LE.',I2,' (0-19)')
      !WRITE 100 315,(MCI3(I),I=1,20)
      IF(NUMB.GT.0) ANMOMS=ANMOMS/NUMB
      IF(NUMBFR.GT.0) SZBOU=SZBOU/NUMBFR
      IF(NUMBFR.GT.0) SZBOU3=SZBOU3/NUMBFR
      IF(NUMBFR.GT.0) SZMAX=SZMAX/NUMBFR
      IF(NUMBFR.GT.0) SMIMF=SMIMF/NUMBFR
      IF(NUMBFR.GT.0) SMIMF2=SMIMF2/NUMBFR
      IF(NUMBFR.GT.0) RMZ=RMZ/NUMBFR
      IF(NUMBFR.GT.0) SNC=SNC/NUMBFR
      DDD=SMIMF2-SMIMF*SMIMF
      IF(DDD.LE.0.) SMIMF2=0.
      IF(DDD.GT.0.) SMIMF2=SQRT(DDD)
      !WRITE 100 1111,SMIMF,SMIMF2
 1111 FORMAT(2X,'mean IMF: SMIMF=',F7.4,'  variance: SMIMF2=',F7.4)
      !WRITE 100 111,SZBOU,SZBOU3,SZMAX,SMIMF,SNC,RMZ
 111  FORMAT(1X,' IN EVENTS : SZBOU=',F7.3,' SZBOU3=',F7.3,
     *' SZMAX=',F7.3,' SMIMF=',F7.3,' SNC=',F7.3,' RMZ=',F7.3)
      !WRITE 100 112,ANMOMS
 112  FORMAT(2X,' avr.ang.mom. ANMOMS (h-bar) =',F8.2)
      !WRITE 100 901
 901  FORMAT(2X,'DISTRIB.OF EVENTS IN BOUND CHARGE OF FRAGM. ',
     *'WITH Z>2 : ZB3  -  MZB3(1-IZ0)')
      !WRITE 100 315,(MZB3(I),I=1,IZ0)
      !WRITE 100 821
 821  FORMAT(2X,'DISTR.OF EVENTS WITH FRAGMENT. IN NPAT',
     *'  : MZFR(1-IZ0)')
      !WRITE 100 315,(MZFR(I),I=1,IZ0)
      DO 820 I=1,IZ0
      IF(MZFR(I).GT.0) ZMZFR(I)=ZMZFR(I)/FLOAT(MZFR(I))
      IF(MZFR(I).GT.0) SFZFR(I)=SFZFR(I)/FLOAT(MZFR(I))
      IF(MZB3(I).GT.0) ZMZB3(I)=ZMZB3(I)/FLOAT(MZB3(I))
      IF(MZB3(I).GT.0) SFZB3(I)=SFZB3(I)/FLOAT(MZB3(I))
      IF(MZB(I).GT.0) TZB(I)=TZB(I)/FLOAT(MZB(I))
      IF(MZB3(I).GT.0) TZB3(I)=TZB3(I)/FLOAT(MZB3(I))
 820  CONTINUE
      !WRITE 100 900
 900  FORMAT(2X,'*****************************************')
 703  FORMAT(2X,10F7.2)
      !WRITE 100 702
 702  FORMAT(2X,'AVERAGE M-can TEMP. IN ZBOU - TZB(1-IZ0)')
      !WRITE 100 703,(TZB(I),I=1,IZ0)
      !WRITE 100 704
 704  FORMAT(2X,'AVERAGE M-can TEMP. IN ZB3 - TZB3(1-IZ0)')
      !WRITE 100 703,(TZB3(I),I=1,IZ0)
      !WRITE 100 902
 902  FORMAT(2X,'AVERAGE MAX CHARGE IN ZB3 - ZMZB3(1-IZ0)')
      !WRITE 100 805,(ZMZB3(I),I=1,IZ0)
      !WRITE 100 903
 903  FORMAT(2X,'AVER.MULT. OF M_IMF IN ZB3 - SFZB3(1-IZ0)')
      !WRITE 100 805,(SFZB3(I),I=1,IZ0)
      !WRITE 100 900
      !WRITE 100 823
 823  FORMAT(2X,'DISTR.OF AVER.MAX CHARGE IN EVENTS WITH FRAGM.',
     *' IN NPAT  : ZMZFR(1-IZ0)')
      !WRITE 100 805,(ZMZFR(I),I=1,IZ0)
      !WRITE 100 825
 825  FORMAT(2X,'AVR.MULT. of M_imf ',
     *' IN NPAT  : SFZFR(1-IZ0)')
      !WRITE 100 805,(SFZFR(I),I=1,IZ0)
 805  FORMAT(2X,10F7.2)
      !WRITE 100 144
 144  FORMAT(2X,'total charge yield NUMZ(1-IZ0):')
      !WRITE 100 145,(NUMZ(I),I=1,IZ0)
 145  FORMAT(2X,10I7)
      DO 146 I=1,IZ0
      IF(NUMZ(I).GT.0) EKZ(I)=EKZ(I)/(2.*I*FLOAT(NUMZ(I)))
      IF(NUMZ(I).GT.0) ETRZ(I)=ETRZ(I)/(2.*I*FLOAT(NUMZ(I)))
 146  CONTINUE
      !WRITE 100 1147
 1147 FORMAT(2X,'avr.kinetic energy (MeV/2Z)   EKZ(1-IZ0):')
      !WRITE 100 805,(EKZ(I),I=1,IZ0)
      !WRITE 100 147
 147  FORMAT(2X,'avr.transvers.energy (MeV/2Z)   ETRZ(1-IZ0):')
      !WRITE 100 805,(ETRZ(I),I=1,IZ0)
      !WRITE 100 148
 148  FORMAT(2X,'azimuth.angl.correl.of 4-He (all events)',2X,
     *' RPHI (0-180 deg., step-9):')
      !WRITE 100 149,(RPHI(I),I=1,20)
      !WRITE 100 160,(RPHI(I),I=21,26)
 149  FORMAT(1X,10F8.0)
 160  FORMAT(1X,6F13.3)
c ----------------------------------
      !WRITE 100 1351
 1351 FORMAT(2X,'max fragm.distrib.in Z-veloc. All Z_max- VPZT ',
     *'(-5 +5cm/ns, 0.25):')
      !WRITE 100 149,(VPZT(I),I=1,40)
      !WRITE 100 160,(VPZT(I),I=41,46)
      !WRITE 100 1352
 1352 FORMAT(2X,'max fragm.in Z-veloc. Z_max=15-25 - VPZ20 ',
     *'(-5 +5cm/ns, step 0.25):')
      !WRITE 100 149,(VPZ20(I),I=1,40)
      !WRITE 100 160,(VPZ20(I),I=41,46)
      !WRITE 100 1353
 1353 FORMAT(2X,'max fragm.in Z-veloc. Z_max=25-35 - VPZ30 ',
     *'(-5 +5cm/ns, step 0.25):')
      !WRITE 100 149,(VPZ30(I),I=1,40)
      !WRITE 100 160,(VPZ30(I),I=41,46)
      !WRITE 100 1354
 1354 FORMAT(2X,'max fragm.in Z-veloc. Z_max=35-45 - VPZ40 ',
     *'(-5 +5cm/ns, step 0.25):')
      !WRITE 100 149,(VPZ40(I),I=1,40)
      !WRITE 100 160,(VPZ40(I),I=41,46)
      !WRITE 100 1355
 1355 FORMAT(2X,'max fragm.in Z-veloc. Z_max=45-.. - VPZ50 ',
     *'(-5 +5cm/ns, step 0.25):')
      !WRITE 100 149,(VPZ50(I),I=1,40)
      !WRITE 100 160,(VPZ50(I),I=41,46)
c ----------------------------------
      CALL THELIP
c      !WRITE 100 180
c 180  FORMAT(2X,'******* DMR ******** for fragments Z=3-20')
c      CALL DMRPR
c      !WRITE 100 165
c 165  FORMAT(2X,'******* POCH1 ********')
c      CALL EPOCH1
      !WRITE 100 166
 166  FORMAT(2X,'******* SPEC1 ********')
      CALL ESPEC1
 500  CONTINUE
      RETURN
      END


      SUBROUTINE HIST1(X,A,B,H,RX,N,W)
      DIMENSION RX(N)
      LE=(B-A)/H+6
      IF(LE-N) 11,11,12
 12   continue
!WRITE 100 13,LE,N
 13   FORMAT(5X,'ERROR IN DIMENSION IN HIST1  LE=',I4,'  N=',I4)
      RETURN
 11   CONTINUE
      RX(N)=RX(N)+X
      IF(X-A) 1,2,2
 1    RX(N-4)=RX(N-4)+W
      RETURN
 2    IF(X-B) 4,3,3
 3    RX(N-2)=RX(N-2)+W
      RETURN
 4    L=(X-A)/H
      RX(L+1)=RX(L+1)+W
      RX(N-1)=RX(N-1)+X
      RX(N-3)=RX(N-3)+W
      RETURN
      END


      SUBROUTINE SORTNA(IA0,RN)
      DIMENSION RN(400),R(400)
      R(1)=RN(1)
      R(2)=(RN(2)+RN(3)+RN(4)+RN(5))/4.
      IA1=IA0-1
      K=(IA1-5)/5
      IA4=K*5+5+1
      DO 1 I=1,K
      RR=0.
      DO 2 J=1,5
 2    RR=RR+RN(I*5+J)
      RR=RR/5.
      R(2+I)=RR
 1    CONTINUE
      IS=2+K+1
      R(IS)=0.
      IF(IA4.GT.IA1) GO TO 4
      DO 3 I=IA4,IA1
 3    R(IS)=R(IS)+RN(I)
      R(IS)=R(IS)/FLOAT(IA1-IA4+1)
      IS=IS+1
 4    CONTINUE
      R(IS)=RN(IA0)
      !WRITE 100 5
 5    FORMAT(2X,'AVERAGE MASS YIELD (IN  A=1,2-5,5-10,...,IA0)
     *ON CALCULATED EVENTS')
      !WRITE 100 6,(R(I),I=1,IS)
 6    FORMAT(5X,10F10.3)
      RETURN
      END

      FUNCTION COLHOT(L,RADNCL)
      COMMON /BLHT06/ZJ(68)/BLHT10/ZFJ(68)/BLHT09/AFJ(68)
     */BLHT05/AJ(68)/BLHT03/U,A,Z
       IF(L-1)1,1,2
1      COLHOT=0.
       RETURN
2      TEMP1=1.44/RADNCL
       COLHOT=TEMP1*((ZJ(L)*ZFJ(L))/(AJ(L)**.33333+AFJ(L)**.33333))
       IF(COLHOT) 3,3,4
3       COLHOT=0.
4       CONTINUE
       RETURN
       END


      FUNCTION GAMHOT(J,A,U,AM,RADNCL)
       COMMON /BLHT09/AFJ(68)/BLHT15/RJ(68)/BLHT14/GAN(68)
      COMMON /BLHT11/VJ(68) /BLHT05/AJ(68)/BLHT06/ZJ(68)
      COMMON /BLHT10/ZFJ(68) /BLHTEX/EXTF(68)
      COMMON /BLHT20/FREP(68)
      IF(RJ(J).LE.0.OR.U.LE.0.) GO TO 12
      PER=2.*SQRT(AM*A*U)
      AMFR=0.5*AM*(1.+(40.**0.333333)/(AJ(J)**0.333333))
      IF(J.GT.6) PER=PER-2.*SQRT(AMFR*AJ(J)*EXTF(J))
      RN=1.5
      CC=0.2
      IF(J.GT.2.AND.J.LE.6)  CC=0.1
      IF(J.GT.6) CC=    (AJ(J)/AFJ(J))**0.6666667
       IF(J-1)1,2,1
2      ALFA=.76+2.2/AFJ(1)**.33333
       BETA=(2.12/AFJ(1)**.66667-.05)/ALFA
      GO TO 3
1      ALFA=1.+CC
       BETA=0.
       GO TO 3
3      Q1=AM*AFJ(J)
       Q2=AM*AFJ(J)*RJ(J)
       Q3=(GAN(J)*AFJ(J)**.66667)*(ALFA/Q1**2)*(AFJ(J)/(AFJ(J)+AJ(J)))
      Q3=Q3*AJ(J)*(3.1416*RN**2)/(2.*41.5*3.1416**2)
      Q3=Q3*FREP(J)
       Q4=(2.*BETA*Q1-3.)/2.+Q2
       Q5=(2.*BETA*Q1-3.)*(SQRT(Q2)-.5)+2.*Q2
      IF(PER-160.) 20,20,21
 20   PEX1=Q4*EXP(-PER)
      GO TO 22
 21   PEX1=0
 22   PP2=PER-2.*SQRT(Q2)
      IF(PP2-160.) 23,23,24
 23   PEX2=Q5*EXP(-PP2)
      GO TO 25
 24   PEX2=0.
 25   GAMHOT=Q3*(PEX1+PEX2)
      IF(J.LE.2) RETURN
      TFORM=AJ(J)*AJ(J)*10.
      GFORM=0.21*940./TFORM
      IF(GAMHOT.GT.GFORM) GAMHOT=GFORM*FREP(J)
      RETURN
 12   GAMHOT=0.
      RETURN
      END


       FUNCTION TKIHOT(L,AM)
       COMMON /BLHT09/AFJ(68)/BLHT15/RJ(68)/BLHT11/VJ(68)
       COMMON /BLHT05/AJ(68)
       RJL=RJ(L)
       RB=4.*AM*AFJ(L)*RJL
       PP1=SQRT(RB)
 5     B1=RNDM(-1)
       IF(B1.LE.0.OR.B1.GT.1) B1=RNDM(-1)
       IF(PP1-160.) 21,21,22
 21    PEX1=EXP(-PP1)
       GO TO 23
 22    PEX1=0.
 23    RK=1.+(1./PP1)*ALOG(B1+(1.-B1)*PEX1)
       IF(L-1) 1,2,1
2      BETA=(2.12/AFJ(1)**0.66667-0.05)/(0.76+2.2/AFJ(1)**0.33333)
       Q1=1.+BETA/RJL
       Q2=Q1*SQRT(Q1)
       FRK=(((3.*SQRT(3.))/2.)/Q2)*(Q1*RK-RK**3)
       GO TO 3
1     FRK=((3.*SQRT(3.))/2.)*(RK-RK**3)
       GO TO 3
 3     B2=RNDM(-1)
       IF(B2-FRK) 4,4,5
4      TKIHOT=  RJL*(1.-RK**2)+VJ(L)
       RETURN
      END


       SUBROUTINE EVAHOT(ENEXT,ATWGHT,CHARGE,PNX,PNY,PNZ,KHOT)
C   EVAPORATION OF HOT FRAGMENTS
       COMMON /BLHT05/AJ(68)
     */BLHT06/ZJ(68)/BLHT14/GAN(68)/BLHT11/VJ(68)/BLHT15/RJ(68)
      COMMON /BLHT03/U,A,Z/BLHT09/AFJ(68)/BLHT10/ZFJ(68)
      COMMON /BLANGL/ANGL(4) /BLHTEX/EXTF(68)
      COMMON /SMPHOT/ SMPA(100),SMPZ(100),SMPE(100),SMPP(100,3)
      COMMON /BLSECON/LSECON
      DIMENSION GJ(68),BJ(68),GG(68),VN(3),PN(3),PNL(3),PP(3),PPL(3)
      DIMENSION GGJ(68)
      U=ENEXT*1000.
      A=ATWGHT
      Z=CHARGE
      REMN=940.*A
      PNL(1)=PNX*1000.
      PNL(2)=PNY*1000.
      PNL(3)=PNZ*1000.
      ENL=SQRT(PNL(1)**2+PNL(2)**2+PNL(3)**2+REMN**2)
      VN(1)=PNL(1)/ENL
      VN(2)=PNL(2)/ENL
      VN(3)=PNL(3)/ENL
      AM=0.125
      RADNCL=1.5
      DO 1 III=1,100
      SMPA(III)=0.
      SMPZ(III)=0.
      SMPE(III)=0.
      SMPP(III,1)=0.
      SMPP(III,2)=0.
      SMPP(III,3)=0.
 1    CONTINUE
      LSECON=0
      IN=68
      DO 20 K=1,100
      UA=U/A
      IF(UA.LE.1.7) GO TO 11
      IF(A.LE.16.) GO TO 11
      DO 4 I=1,IN
      VJ(I)=0.
 4    RJ(I)=-1000.
      CALL DELAM(A,Z,DL1,DSHEL1,BAR1)
      DO 6 I=1,IN
      AFJ(I)=A-AJ(I)
      ZFJ(I)=Z-ZJ(I)
        IF(AFJ(I).LT.ZFJ(I)) GO TO 6
        IF(AFJ(I).LT.AJ(I).OR.ZFJ(I).LT.ZJ(I)) GO TO 6
      VJ(I)=COLHOT(I,RADNCL)
      CALL DELAM(AFJ(I),ZFJ(I),DL2,DSHEL2,BAR2)
      CALL DELAM(AJ(I),ZJ(I),DL3,DSHEL3,BAR3)
      BJ(I)=DL2+DL3-DL1
      RR=U-(BJ(I)+VJ(I))
      EXTF(I)=0.
      IF(RR.GT.0.AND.I.GT.6) EXTF(I)=RR*AJ(I)/A
      RJ(I)=RR-EXTF(I)
 6    CONTINUE
       DO 7 I=1,IN
      GJ(I)=GAMHOT(I,A,U,AM,RADNCL)
7     CONTINUE
      G=0.
      DO 10 I=1,IN
10    G=G+GJ(I)
C ***************
c     IF(G.LE.0.) GO TO 11
c     !WRITE 100 33,U,A,Z
c33   FORMAT(2X,' U=',F7.2,' A=',F5.1,' Z=',F5.1)
c     !WRITE 100 30
c30   FORMAT(10X,'absolute widths:')
c     !WRITE 100 31,(GJ(I),I=1,IN)
c31   FORMAT(2X,' GJ=',10E10.3)
c     DO 32 I=1,IN
c     GGJ(I)=GJ(I)/G
c32   CONTINUE
c     !WRITE 100 35
c35   FORMAT(10X,'relative widths:')
c     !WRITE 100 34,(GGJ(I),I=1,IN)
c34   FORMAT(2X,'GGJ=',10F7.4)
C ***************
        IF(G) 11,11,12
 11   PNX=PNL(1)*0.001
      PNY=PNL(2)*0.001
      PNZ=PNL(3)*0.001
      ENEXT=U*0.001
      ATWGHT=A
      CHARGE=Z
      SMPA(K)=A
      SMPZ(K)=Z
      SMPE(K)=ENEXT
      SMPP(K,1)=PNX
      SMPP(K,2)=PNY
      SMPP(K,3)=PNZ
      IF(K.GT.1) LSECON=1
      KHOT=K
      RETURN
 12   CONTINUE
      DO 13 J=2,IN
13    GJ(J)=GJ(J-1)+GJ(J)
      BB=RNDM(-1)
      B=BB*G
      DO 14 J=1,IN
        IF(B-GJ(J)) 15,14,14
15    LM=J
      GO TO 16
 14   CONTINUE
 16   continue
      IF(LM.LT.1.OR.LM.GT.IN) WRITE(*,300) LM
 300  FORMAT(2X,' ERROR IN EPAHOT - LM=',I5)
      IF(LM.LT.1.OR.LM.GT.IN) then
         write(*,300) LM
         STOP !TG
         GO TO 11
      endif

      EP1=TKIHOT(LM,AM)
      EP2=ZJ(LM)
      EP3=940.*AJ(LM)
c ********
c     !WRITE 100 36,LM,AJ(LM),ZJ(LM),EP1,EXTF(LM)
c36   FORMAT(2X,'LM=',I3,' AJ,ZJ=',2F5.1,' EP1=',F8.3,' EXTF=',
c    *F8.3)
c ********
      U=U-BJ(LM)-EP1-EXTF(LM)
      A=AFJ(LM)
      Z=ZFJ(LM)
C  VPM - RELATIVE VELOCITY OF FRAGMENT
      VPM=SQRT((2.*EP1)/(EP3*AFJ(LM)/(AFJ(LM)+AJ(LM))))
      CALL ISANGL
C  IN CMS
      PP(1)=VPM*ANGL(4)*ANGL(3)*EP3/(1.+AJ(LM)/A)
      PP(2)=VPM*ANGL(4)*ANGL(2)*EP3/(1.+AJ(LM)/A)
      PP(3)=VPM*ANGL(1)*EP3/(1.+AJ(LM)/A)
      PN(1)=-PP(1)
      PN(2)=-PP(2)
      PN(3)=-PP(3)
      EP=SQRT(PP(1)**2+PP(2)**2+PP(3)**2+EP3**2)
      EN=SQRT(PN(1)**2+PN(2)**2+PN(3)**2+(940.*A)**2)
      CALL CLPV(PP,VN,PPL,EP)
      CALL CLPV(PN,VN,PNL,EN)
      SMPA(K)=AJ(LM)
      SMPZ(K)=ZJ(LM)
      SMPE(K)=EXTF(LM)*0.001
      SMPP(K,1)=PPL(1)*0.001
      SMPP(K,2)=PPL(2)*0.001
      SMPP(K,3)=PPL(3)*0.001
      ENL=SQRT(PNL(1)**2+PNL(2)**2+PNL(3)**2+(940.*A)**2)
      VN(1)=PNL(1)/ENL
      VN(2)=PNL(2)/ENL
      VN(3)=PNL(3)/ENL
20    CONTINUE
      !WRITE 100 21,U,A,Z
21    FORMAT(35X,37HMASSIVS SMP EXCEEDED AFTER EVAHOT    /40X,2HU=,
     *F10.5,4H  A=,F5.1,4H  Z=,F4.1)
       RETURN
       END
      BLOCK DATA B1
      COMMON /BLHT05/AJ(68) /BLHT06/ZJ(68) /BLHT14/GAN(68)
      COMMON /BLHT20/FREP(68)
      DATA AJ/
     * 1, 1, 2, 3, 3, 4, 6, 6, 7, 7, 8, 9, 9,10,11,11,12,13,13,14,
     *15,15,16,17,18,18,19,20,20,21,22,22,23,24,24,25,26,26,27,28,
     *28,29,30,30,31,32,32,33,33,34,35,35,36,36,36,37,37,38,39,39,
     *40,40,40,41,41,41,42,42/
      DATA ZJ/
     * 0, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7,
     * 7, 8, 8, 8, 8, 9, 9, 9,10,10,10,11,11,11,12,12,12,13,13,13,
     *14,14,14,15,15,15,16,15,16,16,16,17,16,17,18,17,18,18,18,19,
     *18,19,20,18,19,20,19,20/
      DATA FREP/6*1.,4*0.5,1.,2*0.5,1.,2*0.5,1.,2*0.5,1.,
     * 2*0.5,1.,1.,2*0.5,1.,2*0.5,1.,2*0.5,1.,2*0.5,1.,2*0.5,1.,0.5,
     * 0.5,1.,2*0.5,1.,4*0.5,1.,2*0.5,3*0.333,2*0.5,1.,2*0.5,
     * 6*0.333,2*0.5/
      DATA GAN/2.,2.,3.,2.,2.,1.,62*1./
      END

c         change ZVEZD1, GAMMA1, E11, EPSIL0 and TC
C  PROGRAM OF MULTIFRAGMENTION (WITH RETURN). SMM1 for Soviet VAX
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
      SUBROUTINE ZVEZD1(EE,AA,ZZ,P1,P2,P3,KST)
C MAIN PROGRAM OF MULTIFRAGMENTATION. INPUT: RES.NUCL. - AA,ZZ ; ITS
C EXCITATION - EE (GEV); MOMENTA IN LAB.SYS. - P1,P2,P3 (GEV/C);
C KST - FREE NUMBER IN SPT-MASSIVE (KST=1 - IN THE BEGIN.);
C OUTPUT: ALL PART. AND FRAGM. IN MASSIVE SPT(9,500).
C   *** CHANGED SWITCH.MICR.TO MACR.: INSTEAD SR00 - WW3  ***
c ****** implementation of radial flow (nonrelativistic) *****
c ** and rotation. EFLOW, EMOM in GeV/N. ANMOM in h-bar. *****
c ******  (with relativistic correction - RELCOR)       *****
      REAL N,MU,K,LT,L,KC
      COMMON /BLOKZV/SRN1(500),SRN2(500),SR1,SR2 /ADDNUC/IZPN,IAPN
      COMMON /BLFINL/A0F,Z0F,E0F,TF,FL,SRF,WW3F /BLIMA/IMA
      COMMON /BLOKN/N(500) /BLOKFR/FRA(500),SPSR,S3SR,Z3SR,MFRAG
      COMMON /BLOKF1/EXFRA(500),EKFRA(500)
      COMMON /BLOKF4/SAAM,SEKAM,SEXAM,SSEKX,SSEX
      COMMON /BLPARA/ EPSIL0,FKAPPA,FKACOL /BLMOM/AMOM,AM(500),IAMO
      COMMON /BL503/WCOMP,WW2,WW3,WW4,WW5,WEVAP,WW4G
      COMMON /BL502/M2,M3,M4,M5,MC2,MC3,MC4,MC5
      COMMON /BLMULF/IMULF /BENTRO/SSR /BLHKST/IHFKST(500,2)
      COMMON /BLERR/ERRMI,ERRMA,ERRFR,ERRRE
      COMMON /SMPHOT/ SMPA(100),SMPZ(100),SMPE(100),SMPP(100,3)
      COMMON /BLFLOW/IFLOW,A00,EFLOW /BLAMOM/IAMOM,ANMOM,EMOM
ccc      COMMON /BLJJJ/JJJ
      DIMENSION PA(500),PAZ(500),AMP(500),ZMP(500),EXPA(500)
      DIMENSION PN0(3,500),PNC(3),PNL(3),PN(3),VN(3)
      IAMO=0
      AMOM=0.
      A0=AA
      Z0=ZZ
      IA0=INT(A0+0.5)
C      DO 166 I7=1,IA0
C      AI7=I7
C      IF(I7.LE.4) AM(I7)=0.
C      IF(I7.GT.4) AM(I7)=0.4*2.69*(AI7**1.66667)
C 166  CONTINUE
      SR00=2.60
      IMA4=0
      IF(IMA.EQ.4.AND.A0.LE.110.) IMA4=1
      IF(IMA4.EQ.1) SR00=3.30
      IF(A0.LT.0.9) RETURN
      IREP=0
      E00=(EE/A0)*1000.
      IF(AA-16.) 107,107,108
 107  A=AA
      Z=ZZ
      ERRMI=0
      ERRMA=0
      ERRFR=0
      IF(Z.GT.10.OR.(A-Z).GT.11) GO TO 105
      U=EE
      PN(1)=P1
      PN(2)=P2
      PN(3)=P3
      FRA(IA0)=FRA(IA0)+1.
      EXFRA(IA0)=EXFRA(IA0)+EE*1000./FLOAT(IA0)
      EKFRA(IA0)=EKFRA(IA0)+1000.*(P1**2+P2**2+P3**2)/(2.*0.94*IA0)
      MFRAG=MFRAG+1
      IF(Z.GE.1.) SPSR=SPSR+1.
      IF(Z.GE.3.) S3SR=S3SR+1.
      IF(Z.GE.3.) Z3SR=Z3SR+Z
      SAAM=SAAM+IA0
      SEKAM=SEKAM+1000.*(P1**2+P2**2+P3**2)/(2.*0.94*IA0)
      SEXAM=SEXAM+EE*1000./FLOAT(IA0)
      SSEKX=SSEKX+1000.*(EE+(P1**2+P2**2+P3**2)/(2.*0.94*IA0))
      SSEX=SSEX+1000.*EE
      GO TO 157
 108  CONTINUE
      IF(IMULF.LE.0) GO TO 163
      IF(E00)163,163,165
 163  ERRMI=0
      ERRMA=0
      ERRFR=0
      GO TO 105
 165  CONTINUE
ccc      REMN=0.94*A0
      EDIN=0.
      IF(IFLOW.GE.1) EDIN=EDIN+EFLOW*A0
      IF(IAMOM.GE.1) EDIN=EDIN+EMOM*A0
      REMN=0.940*A0+EE+EDIN
      EWHOLE=SQRT(P1*P1+P2*P2+P3*P3+REMN*REMN)
      VN(1)=P1/EWHOLE
      VN(2)=P2/EWHOLE
      VN(3)=P3/EWHOLE
      K=FKAPPA
      KC=FKACOL
      IF(IMA.EQ.0) WW3=1.
      IF(IMA.EQ.0) GO TO 102
      IF(A0.EQ.A0F.AND.Z0.EQ.Z0F.AND.E00.EQ.E0F) GO TO 60
      ERRMI=0
      ERRMA=0
      ERRFR=0
      E0I=E00
C      IF(IAMO.EQ.1) E0I=E0I+41.5*(AMOM**2)/(2.*AM(IA0)*A0)
CCC      CALL FRATE3(A0,Z0,E0I,K,KC,EPSIL0,TK,SR)
      CALL FRAMIK(A0,Z0,E0I,K,KC,EPSIL0,TK,SR)
      SR1=SR
      DO 201 IJ=1,IA0
 201  SRN1(IJ)=N(IJ)
C *******
c      !WRITE 100 170,SR,TK,SSR
c 170  FORMAT(2X,'AFTER FRAMIK    SR=',F9.4,' TSR=',F9.4,' MEV',
c     *' ENTROP. SSR=',F9.4)
c      !WRITE 100 121,(N(I),I=1,IA0)
c      !WRITE 100 169,WCOMP,WEVAP,WW2,WW3,WW4,WW4G
c 169  FORMAT(2X,' WCOMP=',F8.5,' WEVAP=',F8.5,' WW2=',F8.5,' WW3=',
c     *F8.5,' WW4=',F8.5,' WW4G=',F8.5)
c      !WRITE 100 168,M2,M3,M4,MC4
c 168  FORMAT(5X,' M2=',I10,' M3=',I10,' M4=',I10,' MC4=',I10)
C *******
      GO TO 62
 60   TK=TF
      L=FL
      SR=SRF
      WW3=WW3F
 62   CONTINUE
      IF(IMA4.EQ.1) GO TO 204
      IF(WW3-0.55) 101,101,102
 204  IF(SR-SR00) 101,101,102
 101  CONTINUE
      CALL CHOSAZ(A0,Z0,TK,IP,PA,PAZ)
      GO TO 103
 102  CONTINUE
      IF(A0.EQ.A0F.AND.Z0.EQ.Z0F.AND.E00.EQ.E0F) GO TO 61
      CALL FRAGTE(A0,Z0,E00,K,KC,EPSIL0,TK,MU,L,SR)
      SR2=SR
      DO 202 IJ=1,IA0
 202  SRN2(IJ)=N(IJ)
C *******
c      !WRITE 100 171,SR,TK,SSR
c 171  FORMAT(2X,'AFTER FRAGTE    SR=',F9.4,' TK=',F8.3,
c     *' MEV    ENTROP. SSR=',F9.4)
c      !WRITE 100 121,(N(I),I=1,IA0)
C *******
 61   CONTINUE
      IF(IMA4.EQ.1) GO TO 207
      WWW=4.*(WW3-0.55)
      IF(WWW.LT.0.) WWW=0.
      IF(WWW.GT.1.) WWW=1.
      DO 203 IJ=1,IA0
 203  N(IJ)=(1.-WWW)*SRN1(IJ)+WWW*SRN2(IJ)
      SR=(1.-WWW)*SR1+WWW*SR2
      CALL SELECT(A0,Z0,N,TK,L,IP,PA,PAZ)
      GO TO 103
 207  CONTINUE
      CALL SELECT(A0,Z0,SRN2,TK,L,IP,PA,PAZ)
 103  CONTINUE
      IPROB=1
      DO 500 I=1,IP
      IA=INT(PA(I)+0.5)
      IZ=INT(PAZ(I)+0.5)
      IF(IA.GT.1.AND.IZ.GE.IA) IPROB=0
 500  IF(IA.GT.1.AND.IZ.LE.0) IPROB=0
      IF(IPROB.EQ.1) GO TO 501
      IF(IMA4.EQ.1) GO TO 205
      IF(WW3-0.55) 101,101,61
 205  IF(SR-SR00) 101,101,61
 501  CONTINUE
      A0F=A0
      Z0F=Z0
      E0F=E00
      TF=TK
      FL=L
      SRF=SR
      WW3F=WW3
C *******
ccc      IF(JJJ.LE.71) GO TO 504
c      !WRITE 100 120,SR,TK
c 120  FORMAT(5X,'SR=',F7.3,' TK=',F7.3,' MEV',15X,'N(A)')
c      !WRITE 100 121,(N(I),I=1,IA0)
c 121  FORMAT(2X,10F9.4)
c      !WRITE 100 122,IP
c 122  FORMAT(5X,' IP=',I3,' FRAGMENTS - PA,PAZ:')
c      !WRITE 100 123,(PA(I),PAZ(I),I=1,IP)
c 123  FORMAT(2X,20F6.1)
ccc 504  CONTINUE
C *******
      IF(IP-1) 105,105,106
 105  A=AA
      Z=ZZ
      U=EE
      PX=P1
      PY=P2
      PZ=P3
      FRA(IA0)=FRA(IA0)+1.
      EXFRA(IA0)=EXFRA(IA0)+EE*1000./FLOAT(IA0)
      EKFRA(IA0)=EKFRA(IA0)+1000.*(P1**2+P2**2+P3**2)/(2.*0.94*IA0)
      MFRAG=MFRAG+1
      IF(Z.GE.1.) SPSR=SPSR+1.
      IF(Z.GE.3.) S3SR=S3SR+1.
      IF(Z.GE.3.) Z3SR=Z3SR+Z
      SAAM=SAAM+IA0
      SEKAM=SEKAM+1000.*(P1**2+P2**2+P3**2)/(2.*0.94*IA0)
      SEXAM=SEXAM+EE*1000./FLOAT(IA0)
      SSEKX=SSEKX+1000.*(EE+(P1**2+P2**2+P3**2)/(2.*0.94*IA0))
      SSEX=SSEX+1000.*EE
      GO TO 156
 106  CONTINUE
      T=TK
      CALL FINDT(E00,A0,Z0,PA,PAZ,IP,T,ETP,ECOLB,ECP,EEX,EB)
      IF(ERRFR.LT.0.5) GO TO 111
C *******
ccc      IF(JJJ.LE.71) GO TO 505
c      !WRITE 100 800,IREP,TK,T,EB,ECOLB
c      !WRITE 100 801,(PA(I),I=1,IP)
c 801  FORMAT(2X,29F4.0)
c 800  FORMAT(2X,'IREP=',I2,' TK,T=',2F9.5,' EB,ECOLB=',2F10.4,'   PA:')
ccc 505  CONTINUE
C *******
      IREP=IREP+1
      IF(IREP.GT.10) GO TO 160
      IF(IMA4.EQ.1) GO TO 206
      IF(WW3-0.55) 101,101,61
 206  IF(SR-SR00) 101,101,61
 111  CONTINUE
c **** including implementation additional number of protons.
cc     IZPN=INT(0.5*(79.-ZZ))
cc     IAPN=INT(0.5*(197.-AA))
      IF(IZPN.LT.1.OR.IAPN.LT.1) GO TO 411
      IPOLD1=IP+1
      IPNEW=IP+IZPN
      DO 410 J=IPOLD1,IPNEW
      PA(J)=1.
      PAZ(J)=1.
 410  CONTINUE
      IP=IPNEW
      ZZ0=IZPN+ZZ
      AA0=IAPN+AA
      ECOEFC=(0.6*1.44/(1.17*(1.+KC)**0.333333))
      ECOLN=ECOEFC*ZZ0*ZZ0/(AA0**0.333333)
      DO 412 J4=1,IP
      ECOLN=ECOLN-ECOEFC*PAZ(J4)*PAZ(J4)/(PA(J4)**0.333333)
 412  CONTINUE
      ECOLB=ECOLN
 411  CONTINUE
c  *******
      I=0
      INEUTR=0
      DO 400 J=1,IP
      IF(PAZ(J).LE.0.) GO TO 400
      I=I+1
      AMP(I)=0.94*PA(J)
      ZMP(I)=PAZ(J)
 400  CONTINUE
      IP0=I
      DO 401 J=1,IP
      IF(PAZ(J).GT.0.) GO TO 401
      INEUTR=INEUTR+1
      AMP(IP0+INEUTR)=0.94*PA(J)
      ZMP(IP0+INEUTR)=0.
 401  CONTINUE
C *******
c      !WRITE 100 124,T,IP,INEUTR,IP0,ECOLB,ECP,EEX,EB
c 124  FORMAT(2X,'T=',F7.3,' MEV  IP=',I3,' INEUTR=',I3,' IP0=',I3,
c     *' ECOLB=',F9.4,' MEV  ECP=',F9.4,'  EEX=',F9.4,'  EB=',F9.4)
C *******
      SPSR=SPSR+IP0
      MFRAG=MFRAG+1
      IF(IP0.LE.1) IP0=IP
c including small number of neutrons in general propagation algorithm
      IF(INEUTR.LE.2) IP0=IP
CCC     TP=0.001*(ECOLB+1.5*IP0*T)
CCC     CALL DISIMP(IP0,AMP,PN0,TP)
      ECOLB=0.001*ECOLB
cc      CALL CULIMP(IP0,A0,Z0,AMP,ZMP,PN0,T,ECOLB)
cc      IF(IP0.GE.IP) GO TO 104
cc      TN=0.001*1.5*INEUTR*T
cc      CALL DISNET(INEUTR,IP0,AMP,PN0,T,TN)
c **** implementation of flow and rotation *************
      CALL POS0(IP,IP0,A0,Z0,AMP,ZMP)
      CALL CULIM0(IP0,A0,Z0,AMP,ZMP,PN0,T,ECOLB)
      IF(IP0.GE.IP) GO TO 104
      TN=0.001*1.5*INEUTR*T
c      CALL DISNE0(INEUTR,IP0,AMP,PN0,T,TN)
      CALL DISN02(INEUTR,IP0,AMP,PN0,T,TN)
c --- statistic on energy of hot fragments before coulomb ---
      CALL BFREN(INEUTR,IP0,AMP,PN0)
c **********************************************************
 104  CONTINUE
c --- finding new surface parameters
      CALL PARAM2(T,EPSIL0)
c --- statistics on microcanonical temperature
      CALL FIXT(T,K,IP,AMP,ZMP)
c --- statistics on hot fragments
      CALL HOTFR(IP0,A0,Z0,AMP,ZMP,PN0,T,ECOLB)
c -------------------------------
c **** relativistic correction ****
c      CALL RELCOR(IP,AMP,ZMP,PN0)
c ****
ccc      IF(JJJ.GE.72) !WRITE 100 506
ccc 506  FORMAT(2X,'print after DISNE0')
      KSTI=1
      IAAM=0
      DO 35 I=1,IP
      DO 36 J=1,3
 36   PNC(J)=PN0(J,I)
C *******
c      EPAT=(PNC(1)**2+PNC(2)**2+PNC(3)**2)/(2.*AMP(I))
c      !WRITE 100 125,IP,AMP(I),(PNC(J),J=1,3),EPAT
c 125  FORMAT(2X,'PARTICLE IP=',I3,' AMP=',F7.3,'  PNC(1-3) (GEV/C)=',
c     *3F7.3,' EPAT (GEV)=',F7.3)
C *******
      EN=SQRT(PNC(1)**2+PNC(2)**2+PNC(3)**2+AMP(I)**2)
      CALL CLPV(PNC,VN,PNL,EN)
C *******
c      EPAT=(PNL(1)**2+PNL(2)**2+PNL(3)**2)/(2.*AMP(I))
c      !WRITE 100 126,IP,AMP(I),(PNL(J),J=1,3),EPAT
c 126  FORMAT(2X,'B  [.C.K. IP=',I3,' AMP=',F7.3,'  PNL(1-3) (GEV/C)=',
c     *3F7.3,' EPAT (GEV)=',F7.3)
C *******
      A=AMP(I)/0.94
      Z=ZMP(I)
      IA=INT(A+0.5)
      IZ=INT(Z+0.5)
      IF(IZ.GE.3) S3SR=S3SR+1.
      IF(IZ.GE.3) Z3SR=Z3SR+Z
      IF(IA.GT.0.AND.IA.LE.500) FRA(IA)=FRA(IA)+1.
CCC     GO TO 35
      IF(IA-3) 51,51,52
 51   U=0.
      GO TO 53
 52   U=UEVA(A,T)
 53   CONTINUE
      IF(IA.GT.0.AND.IA.LE.500.AND.IA.GT.IAAM) EXAM=U*1000./FLOAT(IA)
      IF(IA.GT.0.AND.IA.LE.500.AND.IA.GT.IAAM)
     *EKAM=1000.*(PNL(1)**2+PNL(2)**2+PNL(3)**2)/(2.*0.94*IA)
      IF(IA.GT.0.AND.IA.LE.500.AND.IA.GT.IAAM) IAAM=IA
      IF(IA.GT.0.AND.IA.LE.500)
     *EKX=1000.*(U+(PNL(1)**2+PNL(2)**2+PNL(3)**2)/(2.*0.94*IA))
      SSEKX=SSEKX+EKX
      IF(IA.GT.0.AND.IA.LE.500) SSEX=SSEX+1000.*U
      IF(IA.GT.0.AND.IA.LE.500)
     *EXFRA(IA)=EXFRA(IA)+U*1000./FLOAT(IA)
      IF(IA.GT.0.AND.IA.LE.500.AND.I.LE.500)
     *EXPA(I)=EXPA(I)+U*1000./FLOAT(IA)
      IF(IA.GT.0.AND.IA.LE.500) EKFRA(IA)=EKFRA(IA)
     *+1000.*(PNL(1)**2+PNL(2)**2+PNL(3)**2)/(2.*0.94*IA)
      PX=PNL(1)
      PY=PNL(2)
      PZ=PNL(3)
      CALL EVAHOT(U,A,Z,PX,PY,PZ,KHT)
      DO 703 I3=1,KHT
      U=SMPE(I3)
      A=SMPA(I3)
      Z=SMPZ(I3)
      PX=SMPP(I3,1)
      PY=SMPP(I3,2)
      PZ=SMPP(I3,3)
      IA3=INT(A+0.5)
      IZ3=INT(Z+0.5)
c *******************
c      !WRITE 100 706,A,Z,U,PX,PY,PZ
c 706  FORMAT(2X,'AFTER EVAHOT: A,Z=',2F5.1,' U=',F8.4,
c     *' PX,PY,PZ=',3F8.4)
c *******************
      IF(IZ3.GT.10.OR.(IA3-IZ3).GT.11) GO TO 705
      IF(IA3-16) 704,704,705
 704  PN(1)=PX
      PN(2)=PY
      PN(3)=PZ
c -------- changing secondary de-excit. ----------
      CALL RAZVAL(U,A,Z,PN,KST)
c      CALL EVANUC(U,A,Z,PX,PY,PZ,KST)
c ------------------------------------------------
      GO TO 703
 705  CONTINUE
      CALL EVANUC(U,A,Z,PX,PY,PZ,KST)
 703  CONTINUE
c ---------------
      KSTF=KST-1
      IHFKST(I,1)=KSTI
      IHFKST(I,2)=KSTF
      KSTI=KST
c ---------------
 35   CONTINUE
      SAAM=SAAM+IAAM
      SEKAM=SEKAM+EKAM
      SEXAM=SEXAM+EXAM
      KST1=KST-1
      CALL HOTFR6(IP,KST1,A0,Z0,AMP,ZMP,PN0,T,ECOLB,EXPA)
      GO TO 160
c -------- changing secondary de-excit. ----------
 157  CONTINUE
      CALL RAZVAL(U,A,Z,PN,KST)
c      PX=PN(1)
c      PY=PN(2)
c      PZ=PN(3)
c      CALL EVANUC(U,A,Z,PX,PY,PZ,KST)
c ------------------------------------------------
      GO TO 160
 156  CALL EVAHOT(U,A,Z,PX,PY,PZ,KHT)
      DO 700 I2=1,KHT
      U=SMPE(I2)
      A=SMPA(I2)
      Z=SMPZ(I2)
      PX=SMPP(I2,1)
      PY=SMPP(I2,2)
      PZ=SMPP(I2,3)
      IA2=INT(A+0.5)
      IZ2=INT(Z+0.5)
c *******************
c      !WRITE 100 706,A,Z,U,PX,PY,PZ
c *******************
      IF(IZ2.GT.10.OR.(IA2-IZ2).GT.11) GO TO 701
      IF(IA2-16) 702,702,701
 702  PN(1)=PX
      PN(2)=PY
      PN(3)=PZ
c -------- changing secondary de-excit. ----------
      CALL RAZVAL(U,A,Z,PN,KST)
c      CALL EVANUC(U,A,Z,PX,PY,PZ,KST)
c ------------------------------------------------
      GO TO 700
 701  CALL EVANUC(U,A,Z,PX,PY,PZ,KST)
 700  CONTINUE
      GO TO 160
 160  continue
!IF(ERRMI.GT.0.5.OR.ERRMA.GT.0.5) !WRITE 100 161,AA,ZZ,EE,ERRMI,ERRMA
 161  FORMAT(2X,'SUSPECTING EVENT: AA,ZZ,EE=',3F9.4,' ERRMI,ERRMA=',
     *2F9.1)
!      IF(ERRFR.GT.0.5) !WRITE 100 162,AA,ZZ,EE,ERRFR
 162  FORMAT(5X,'EVENT EXCLUDED: AA,ZZ,EE=',3F9.4,' ERRFR=',F9.1)
      RETURN
      END

      BLOCK DATA B2
      COMMON /BLMULF/IMULF /BLIMA/IMA /BLPARA/EPSIL0,FKAPPA,FKACOL
      COMMON /RESNU0/EMIN,EMAX,YMIN,YMAX /BPLACE/RFI /ADDNUC/IZPN,IAPN
      COMMON /BLJPRE/JPRE /CELIPS/CEL,IPOSF1
      COMMON /BPERI/IPERI
      DATA CEL/1./,IPOSF1/0/
      DATA IMULF/1/,IMA/3/,EPSIL0/16./,FKAPPA/1./,FKACOL/2./
      DATA IZPN/0/,IAPN/0/,RFI/0./,EMIN/0./,EMAX/0./,YMIN/0./,YMAX/0./
      DATA JPRE/0/,IPERI/0/
      END

      SUBROUTINE FRAMIK(AA,ZZ,E00,K,KC,EPSIL0,TSR,MSR)
C HAXO[[EH[E [O [HEP[[[ E00 (MEV/N) M[[[T[[PA[MEHTHO[O COCTO[H[[
C B O[[ACT[ MA[[X [HEP[[[ BO[[[[[EH[[
C   ***********    KBA[[-M[KPOKAHOH[KA    *********
      REAL  K,KC,LT,MSR
      COMMON /BL500/W,W2(250),W3(25000),W4(10000) /BLIMA/IMA
      COMMON /BL501/IAF2(250),IAF3(25000),IAF4(10000)
      COMMON /BL502/M2,M3,M4,M5,MC2,MC3,MC4,MC5
      COMMON /BL503/WCOMP,WW2,WW3,WW4,WW5,WEVAP,WW4G
      COMMON /BLDIS/ J(20),IEND,JZ(20),IENDZ
      COMMON /BLOK1/W0,TC,RADNCL,G0 /BLOK0/A0,Z0,V0,XZ0,IA0,E0,E0Q
      COMMON /BLOK3/E11(500),XZ(500),GA(500) /BLOKT/TCON
      COMMON /BLOK2/B0,BT,DBT,LT,CP,PKP013 /BLOKN/SRW(500)
      COMMON /BLERR/ERRMI,ERRMA,ERRFR,ERRRE /BENTRO/SSR
      COMMON /BPERI/IPERI
      ERRMI=0.
      A0=AA
      Z0=ZZ
      IA0=INT(A0+0.5)
      E11(1)=0.
      DO 100 I=2,IA0
      A=I
      IF(IPERI.EQ.0) E11(I)=EPSIL0*(1.+3./(A-1))
C      IF(I.LE.50) E11(I)=EPSIL0*(1.+0.12*(A**0.33333-5.**0.33333))
c 100  IF(I.GT.50) E11(I)=EPSIL0*(1.+0.12*(50.**0.33333-5.**0.33333))
      IF(IPERI.EQ.1) E11(I)=EPSIL0/(1.+0.0020*(A/25.)**2)
 100  CONTINUE
      TCON=SQRT(E00/0.125)
      T=0.
      CALL PARAM(T,K,KC,EPSIL0)
      XZ0=Z0/A0
      A013=A0**0.3333333
      V0=(4.*3.1416/3.)*A0*RADNCL**3
c --- NOT changing of coulomb volume : V_c=V_0+V_f -----
      PKP013=1./(1.+KC)**0.3333333
      CP=0.6*1.44*(1.-PKP013)/RADNCL
c ------------------------------------
      E0Q=A0*(-W0+G0*(1.-2*XZ0)**2)+B0*A013*A013
     *+0.6*1.44*Z0*Z0/(RADNCL*A013)
      DO 1 I=1,IA0
      A=I
      GA(I)=1.
      XZ(I)=XZ0
      SRW(I)=0.
 1    CONTINUE
      XZ(2)=0.5
      XZ(3)=0.5
      XZ(4)=0.5
      GA(2)=3.
      GA(3)=4.
      GA(4)=1.
      GA(1)=4.
      W=0.
      MSR=0.
      TSR=0.
      SSR=0.
      WW2=0.
      WW3=0.
      WW4=0.
      WW5=0.
      WEVAP=0.
      WW4G=0.
      CALL SCOMPA(E00,K,KC,EPSIL0,SCOMP)
      WCOMP=EXP(SCOMP-SCOMP)
      W=W+WCOMP
      WEVAP=WEVAP+WCOMP
      MSR=MSR+WCOMP*1.
      TSR=TSR+TCON*WCOMP
      SSR=SSR+WCOMP*SCOMP
      SRW(IA0)=SRW(IA0)+WCOMP
      M2=0
      M3=0
      M4=0
      M5=0
      MMAX=3
      IF(IMA.EQ.2) MMAX=2
      IF(IMA.EQ.4.AND.IA0.LE.110) MMAX=4
      DO 10 M=2,MMAX
      MC=0
      J(M)=IA0
      M1=M-1
      DO 2 I=1,M1
 2    J(I)=0
 4    CALL DISA(M)
      IF(IEND.EQ.1) GO TO 10
      MC=MC+1
c --- changing of coulomb volume : V_c=V_0+V_f -----
c      KC=((1.+1.4*(M**0.333333-1)/(RADNCL*A013))**3-1.)
c      PKP013=1./(1.+KC)**0.3333333
c      CP=0.6*1.44*(1.-PKP013)/RADNCL
c ------------------------------------
      CALL CALCWE(M,J,SCOMP,E00,K,KC,EPSIL0,PR,ECOLB)
C ***********
C      IF(MC.LE.20) !WRITE 100 101,M,TCON,PR,ECOLB,(J(I),I=1,M)
C 101  FORMAT(2X,'M=',I3,' TCON=',F8.3,' PR=',F8.5,' ECOLB=',F8.3,' J=',
C     *5I4)
C ***********
      W=W+PR
      DO 6 I=1,M
      IA=J(I)
 6    SRW(IA)=SRW(IA)+PR
      MSR=MSR+M*PR
      TSR=TSR+TCON*PR
      IF(PR.GT.0.) SSR=SSR+PR*(SCOMP+ALOG(PR))
      IF(M.EQ.2) GO TO 22
      IF(M.EQ.3) GO TO 23
      IF(M.EQ.4) GO TO 24
      IF(M.EQ.5) GO TO 25
 22   MC2=MC
      WW2=WW2+PR
      IF(J(1).LE.4) WEVAP=WEVAP+PR
      M2=M2+1
      W2(M2)=PR
      IAF2(M2)=J(2)+1000*J(1)
      GO TO 4
 23   MC3=MC
      WW3=WW3+PR
      M3=M3+1
      W3(M3)=PR
      IAF3(M3)=J(3)+1000*J(2)+1000000*J(1)
      GO TO 4
 24   MC4=MC
      WW4=WW4+PR
      M4=M4+1
      IF(J(3).LE.4) WW4G=WW4G+PR
      IF(M4.GT.10000) GO TO 4
      W4(M4)=PR
      IAF4(M4)=J(4)+1000*J(3)+100000*J(2)+10000000*J(1)
      GO TO 4
 25   MC5=MC
      WW5=WW5+PR
      M5=M5+1
      GO TO 4
 10   CONTINUE
      IF(M4.GT.10000) M4=10000
      DO 12 I=1,IA0
 12   SRW(I)=SRW(I)/W
      DO 63 I8=1,M2
 63   W2(I8)=W2(I8)/W
      DO 64 I8=1,M3
 64   W3(I8)=W3(I8)/W
      DO 65 I8=1,M4
 65   W4(I8)=W4(I8)/W
      MSR=MSR/W
      TSR=TSR/W
      SSR=SSR/W
      WCOMP=WCOMP/W
      WW2=WW2/W
      WW3=WW3/W
      WW4=WW4/W
      WW4G=WW4G/W
      WW5=WW5/W
      WEVAP=WEVAP/W
      RETURN
      END


      SUBROUTINE CALCWE(M,J,SCOMP,E00,K,KC,EPSIL0,PR,ECOLB)
C  O[PE[E[EH[E BEPO[THOCT[ KOH[[[[PA[[[ PR=EXP(S(J(1-M),E0)
C       [ K[[OHOBCKO[O [AP[EPA ECOLB (MEV).
C  BXO[:M [PA[MEHTOB B MACC[BE J(1-M),[HEP[.BO[[.-E00 (MEV/N),
C  [APAMETP[ MO[E[[ K,KC,EPSIL0 [ [HTPO[[[ KOM[.[[PA SCOMP ([[[ [EPEC-
C  [ETA PR).
      REAL LT,K,KC
      COMMON /BLOK1/W0,TC,RADNCL,G0 /BLOK0/A0,Z0,V0,XZ0,IA0,E0,E0Q
      COMMON /BLOK3/E11(500),XZ(500),GA(500) /BLOKT/TCON
      COMMON /BLOK2/B0,BT,DBT,LT,CP,PKP013 /BLMOM/AMOM,AM(500),IAMO
      COMMON /BLERR/ERRMI,ERRMA,ERRFR,ERRRE
      DIMENSION EA(500),ECOLA(500),SAIN(500),J(20)
      E0=E00*A0
      A013=A0**0.333333
c --- keep constant or vary with multiplic. translational volume ---
c      VF=K*V0
      VF=V0*((1.+1.4*(M**0.333333-1)/(RADNCL*A013))**3-1.)
c ---
      M1=M-1
      FAKT=1.
      DO 32 I=1,M1
      FFF=1.
      IM=I+1
      DO 33 II=IM,M
      IF(J(I).EQ.J(II)) FFF=FFF+1.
 33   CONTINUE
      FAKT=FAKT*FFF
 32   CONTINUE
      PRFAKT=1.
      PRFAKT=PRFAKT*FAKT
      ECON=0.
      ECOLB=0.
      PRGA=1.
      PRA32=1.
      PAM32=1.
C      SAM=0.
      DO 34 I=1,M
      IA=J(I)
C      SAM=SAM+AM(IA)
C      IF(AM(IA).GT.0.) PAM32=PAM32*(AM(IA)*SQRT(AM(IA)))
      A=IA
      PRGA=PRGA*GA(IA)
      PRA32=PRA32*(A*SQRT(A))
      IF(IA-1) 41,41,42
 41   ECOLA(IA)=CP*XZ(IA)*XZ(IA)
      EA(IA)=ECOLA(IA)
      GO TO 50
 42   IF(IA-4) 43,43,44
 43   ECOLA(IA)=CP*XZ(IA)*XZ(IA)*(A**(5./3.))
      IF(IA.EQ.2) EA(IA)=-2.796+ECOLA(IA)
      IF(IA.EQ.3) EA(IA)=-9.224+ECOLA(IA)
      IF(IA.EQ.4) EA(IA)=-30.11+ECOLA(IA)
      GO TO 50
 44   ECOLA(IA)=CP*XZ(IA)*XZ(IA)*(A**(5./3.))
      EA(IA)=(-W0+G0*(1.-2.*XZ(IA))**2)*A+
     &B0*A**0.666667+ECOLA(IA)
 50   CONTINUE
      ECON=ECON+EA(IA)
      ECOLB=ECOLB+ECOLA(IA)
 34   CONTINUE
C      EMOM=41.5*(AMOM**2)/(2*SAM)
      ECOLB=ECOLB+0.6*1.44*Z0*Z0*PKP013/(RADNCL*A013)
      ECON=ECON+0.6*1.44*Z0*Z0*PKP013/(RADNCL*A013)
C      IF(IAMO.EQ.1) ECON=ECON+EMOM
      DO 51 I=1,M
      IA=J(I)
      A=IA
      ECOLB=ECOLB-0.6*1.44*XZ(IA)*XZ(IA)*(A**(5./3.))/RADNCL
 51   CONTINUE
      IF((E0+E0Q-ECON).LT.0.003) GO TO 56
      T=SQRT(E00/0.125)
      IF(T.LT.0.0012) T=0.0012
      HT=0.5
      K1=0
      K2=0
      ID=0
      ICNT=0
 24   CONTINUE
      IF(ICNT.GT.120) GO TO 71
      CALL PARAM(T,K,KC,EPSIL0)
      ECON=0.
      DO 61 I=1,M
      IA=J(I)
      A=IA
      IF(IA-1) 55,55,52
 55   EA(IA)=ECOLA(IA)
      GO TO 60
 52   IF(IA-4) 53,53,54
 53   IF(IA.EQ.2) EA(IA)=-2.796+ECOLA(IA)
      IF(IA.EQ.3) EA(IA)=-9.224+ECOLA(IA)
      IF(IA.EQ.4) EA(IA)=-30.11+ECOLA(IA)+4.*T*T/E11(4)
      GO TO 60
 54   EA(IA)=(-W0+T*T/E11(IA)+G0*(1.-2.*XZ(IA))**2)*A+
     &(BT-T*DBT)*A**0.666667+ECOLA(IA)
 60   CONTINUE
      ECON=ECON+EA(IA)
 61   CONTINUE
      ECON=ECON+0.6*1.44*Z0*Z0*PKP013/(RADNCL*A013)+1.5*T*(M-1)
C      IF(IAMO.EQ.1) ECON=ECON+EMOM+1.5*T*(M-1)
      D=(E0+E0Q-ECON)/E0
      IF(ABS(D).LT.0.003) GO TO 29
      ICNT=ICNT+1
      H=SIGN(HT,D)
      IF(D) 21,21,22
 21   K1=1
      IF(K1.EQ.1.AND.K2.EQ.1) ID=ID+1
      IF(ID.GT.29) GO TO 71
      T=T+H/2**ID
      IF(T.GE.0.001) GO TO 24
      K2=1
      H=HT
      GO TO 21
 22   K2=1
      IF(K1.EQ.1.AND.K2.EQ.1) ID=ID+1
      IF(ID.GT.29) GO TO 71
      T=T+H/2**ID
      IF(T.GE.0.001) GO TO 24
      K1=1
      H=HT
      GO TO 22
 71   ERRMI=ERRMI+1
C      !WRITE 100 72,ICNT,ID,D,T
C 72   FORMAT(5X,'ERROR IN CALCWE  ICNT=',I3,' ID=',I3,' D=',E12.5,
C     *'  T=',F8.3)
 56   PR=0.
      RETURN
 29   TCON=T
      SCON=0.
      DO 62 I=1,M
      IA=J(I)
      A=IA
      IF(IA-4) 63,63,64
 63   SAIN(IA)=0.
      IF(IA.EQ.4) SAIN(IA)=2*T*4./E11(IA)
      GO TO 65
 64   SAIN(IA)=2*T*A/E11(IA)-DBT*A**0.666667
 65   CONTINUE
      SCON=SCON+SAIN(IA)
 62   CONTINUE
      LT=16.15/SQRT(T)
      STRAN=ALOG(PRA32/PRFAKT)+(M-1)*ALOG(VF/LT**3)+1.5*(M-1)-
     *ALOG(A0*SQRT(A0))
      IF(STRAN.LT.0.) STRAN=0.
      SCON=SCON+ALOG(PRGA)+STRAN
C      SAMOM=1.5*(M-1)-3*(M-1)*ALOG(LT)+ALOG(PAM32)-
C     *ALOG(SAM*SQRT(SAM))
C      IF(SAMOM.LT.0.) SAMOM=0.
C      IF(IAMO.EQ.1) SCON=SCON+SAMOM
      PR=EXP(SCON-SCOMP)
C     IF(M.LE.2) !WRITE 100 100,M,(J(I),I=1,M),ECOLB,T,STRAN,SAMOM,SCON,PR
C100  FORMAT(2X,'M=',I3,' J=',2I3,' ECOLB=',F7.3,' T=',F7.3,' STRAN=',
C    *F7.3,' SAMOM=',F7.3,' SCON=',F7.3,' PR=',E10.3)
      RETURN
      END
      SUBROUTINE DISA(K)
C B[[OP O[HO[O PAC[PE[E[EH[[ A0 H[K[OHOB [O K [PA[MEHTAM C [[ETOM HEPA[-
C [[[[MOCT[. KOH[[[[PA[[[: J(1)...J(K)-CO[EP[[T MACC[ [PA[MEHTOB.
C HA[. [C[OB[E: J(K)=A0,J(K-1)=...J(1)=O -[A[AETC[ [EPE[ DISA.
C [P[[HAK KOH[A IEND=1 O[HA[AET,[TO BCE KOH[[[[PA[[[ [[E [C[EP[AH[.
C [P[[HAK BO[MO[HO[ KOH[[[[PA[[[: IEND=0 .
      COMMON /BLDIS/ J(20),IEND,JZ(20),IENDZ
      L=0
 1    L=L+1
      IF(L-K) 3,4,4
 4    IEND=1
      RETURN
 3    JL=J(L)
      JK=J(K)
      J(L)=J(L)+1
      J(K)=J(K)-1
      IF(J(L).GT.J(L+1).OR.J(K-1).GT.J(K)) GO TO 2
      IEND=0
      RETURN
 2    J(L)=1
      J(K)=JK+JL-1
      GO TO 1
      END
      SUBROUTINE CHOSAZ(A0,Z0,T,IP,PA,PAZ)
C B[[OP KAHA[A PA[BA[A [[ BCEX KOH[[[[PA[[[ [O M=4 .
C [A[O[HEH[E MACC[BA PA(1-IP) - B [OP[[KE [[[BAH[[ A.
      DIMENSION PA(500),PAZ(500)
      COMMON /BL500/W,W2(250),W3(25000),W4(10000)
      COMMON /BL501/IAF2(250),IAF3(25000),IAF4(10000)
      COMMON /BL502/M2,M3,M4,M5,MC2,MC3,MC4,MC5
      COMMON /BL503/WCOMP,WW2,WW3,WW4,WW5,WEVAP,WW4G
      WCOWW2=WCOMP+WW2
      WCOWW3=WCOWW2+WW3
      BR=RNDM(-1)
      IP=1
      IF(BR.GT.WCOMP.AND.BR.LE.WCOWW2) IP=2
      IF(BR.GT.WCOWW2.AND.BR.LE.WCOWW3) IP=3
      IF(BR.GT.WCOWW3.AND.WW4.GT.0.) IP=4
      GO TO (1,2,3,4),IP
 1    PA(1)=A0
      PAZ(1)=Z0
      RETURN
 2    WP=WCOMP
      JJ=0
      DO 10 I=1,M2
      WP=WP+W2(I)
      JJ=JJ+1
      IF(WP.GT.BR) GO TO 11
 10   CONTINUE
 11   MCH=JJ
      PA(2)=IAF2(MCH)/1000
      PA(1)=IAF2(MCH)-INT(PA(2))*1000
      GO TO 5
 3    WP=WCOMP+WW2
      JJ=0
      DO 12 I=1,M3
      WP=WP+W3(I)
      JJ=JJ+1
      IF(WP.GT.BR) GO TO 13
 12   CONTINUE
 13   MCH=JJ
      PA(3)=IAF3(MCH)/1000000
      PA(2)=(IAF3(MCH)-INT(PA(3))*1000000)/1000
      PA(1)=IAF3(MCH)-INT(PA(3))*1000000-INT(PA(2))*1000
      GO TO 5
 4    WP=WCOMP+WW2+WW3
      JJ=0
      DO 14 I=1,M4
      JJ=JJ+1
      WP=WP+W4(I)
      IF(WP.GT.BR) GO TO 15
 14   CONTINUE
 15   MCH=JJ
      PA(4)=IAF4(MCH)/10000000
      PA(3)=(IAF4(MCH)-INT(PA(4))*10000000)/100000
      PA(2)=(IAF4(MCH)-INT(PA(4))*10000000-INT(PA(3))*100000)/1000
      PA(1)=IAF4(MCH)-INT(PA(4))*10000000-INT(PA(3))*100000
     *-INT(PA(2))*1000
 5    CONTINUE
      CALL CHOSZ(A0,Z0,T,IP,PA,PAZ)
      RETURN
      END
      SUBROUTINE CHOSZ(A0,Z0,T,IP,PA,PAZ)
      DIMENSION PA(500),PAZ(500)
C O[PE[E[EH[E [AP[[OB IP [PA[MEHTOB B MACC[BE PAZ(1-IP)
      IZ0=INT(Z0+0.5)
      G0=25.
 5    ISZ=0
      DO 1 I=1,IP
      IA=INT(PA(I)+0.5)
 3    CALL RANNOR(BR1,BR2)
      CC=8.*G0
      ZM=PA(I)*Z0/A0
c      IF(PA(I).GT.1.5.AND.PA(I).LT.4.5) ZM=0.5*PA(I)
      IF(IA.EQ.2) ZM=0.5*PA(I)
      DZ=SQRT(PA(I)*T/CC)
      Z=BR1*DZ+ZM
      IZ=INT(Z+0.5)
      IF(IZ.GE.0.AND.IZ.LE.IA) GO TO 4
      Z=BR2*DZ+ZM
      IZ=INT(Z+0.5)
      IF(IZ.GE.0.AND.IZ.LE.IA) GO TO 4
      GO TO 3
 4    PAZ(I)=IZ
      ISZ=ISZ+INT(PAZ(I)+0.5)
 1    CONTINUE
      DELZ=IZ0-ISZ
      IF(ABS(DELZ).GT.1.1) GO TO 5
      PAZ(1)=PAZ(1)+DELZ
      RETURN
      END
      SUBROUTINE PARAM(T,K,KC,EPSIL0)
      REAL K,KC
      COMMON /BLOK1/W0,TC,RADNCL,G0 /BLOK0/A0,Z0,V0,XZ0,IA0,E0,E0Q
      COMMON /BLOK3/E11(500),XZ(500),GA(500)
      COMMON /BLOK2/B0,BT,DBT,LT,CP,PKP013
      W0=16.
      B0=18.
      TC=18.
c      TC=1000000.
      G0=25.
      RADNCL=1.17
      IF(T.LT.TC) GO TO 81
      BT=0.
      DBT=0.
      GO TO 82
 81   CONTINUE
      BT=B0*((TC*TC-T*T)/(TC*TC+T*T))**1.25
      DBT=B0*(-5*T*TC*TC/(TC*TC+T*T)**2)*((TC*TC-T*T)/(TC*TC+T*T))**.25
 82   CONTINUE
      IF(EPSIL0.GT.9999.) BT=B0
      IF(EPSIL0.GT.9999.) DBT=0.
      RETURN
      END
      SUBROUTINE PARAM2(T,EPSIL0)
      COMMON /BLOK1/W0,TC,RADNCL,G0
      COMMON /BLOK2/B0,BT,DBT,LT,CP,PKP013
      COMMON /BLOK20/BT2,DBT2
c for FINDT : finding surface energy of fragments
      IF(T.LT.TC) GO TO 81
      BT2=0.
      DBT2=0.
      GO TO 82
 81   CONTINUE
      BT2=B0*((TC*TC-T*T)/(TC*TC+T*T))**1.25
      DBT2=B0*(-5*T*TC*TC/(TC*TC+T*T)**2)*
     &((TC*TC-T*T)/(TC*TC+T*T))**.25
 82   CONTINUE
      IF(EPSIL0.GT.9999.) BT2=B0
      IF(EPSIL0.GT.9999.) DBT2=0.
      RETURN
      END
      SUBROUTINE SCOMPA(E00,K,KC,EPSIL0,SCOMP)
C  O[PE[E[EH[E [HTPO[[[ KOM[A[H[ COCTO[H[[
      COMMON /BLOK1/W0,TC,RADNCL,G0 /BLOK0/A0,Z0,V0,XZ0,IA0,E0,E0Q
      COMMON /BLOK3/E11(500),XZ(500),GA(500) /BLOKT/TCON
      COMMON /BLOK2/B0,BT,DBT,LT,CP,PKP013 /BLMOM/AMOM,AM(500),IAMO
      COMMON /BLERR/ERRMI,ERRMA,ERRFR,ERRRE
      REAL K,KC,LT
      E0=E00*A0
C      EMOM=41.5*(AMOM**2)/(2*AM(IA0))
C      IF((E0-EMOM).LT.0.003.AND.IAMO.EQ.1) GO TO 73
      A013=A0**0.333333
      T=SQRT(E00/0.125)
      IF(T.LT.0.0012) T=0.0012
      HT=0.5
      K1=0
      K2=0
      ID=0
      ICNT=0
 24   CONTINUE
      IF(ICNT.GT.120) GO TO 71
      CALL PARAM(T,K,KC,EPSIL0)
      ECON=(-W0+T*T/E11(IA0)+G0*(1.-2.*XZ(IA0))**2)*A0+
     *(BT-T*DBT)*A013*A013+0.6*1.44*Z0*Z0/(RADNCL*A013)
C      IF(IAMO.EQ.1) ECON=ECON+EMOM
      D=(E0+E0Q-ECON)/E0
      IF(ABS(D).LT.0.003) GO TO 29
      ICNT=ICNT+1
      H=SIGN(HT,D)
      IF(D) 21,21,22
 21   K1=1
      IF(K1.EQ.1.AND.K2.EQ.1) ID=ID+1
      IF(ID.GT.29) GO TO 71
      T=T+H/2**ID
      IF(T.GE.0.001) GO TO 24
      K2=1
      H=HT
      GO TO 21
 22   K2=1
      IF(K1.EQ.1.AND.K2.EQ.1) ID=ID+1
      IF(ID.GT.29) GO TO 71
      T=T+H/2**ID
      IF(T.GE.0.001) GO TO 24
      K1=1
      H=HT
      GO TO 22
 71   ERRMI=ERRMI+1
C      !WRITE 100 72,ICNT,ID,D,T
C 72   FORMAT(5X,'ERROR IN SCOMPA  ICNT=',I3,' ID=',I3,' D=',E12.5,
C     *'  T=',F8.3)
 73   SCOMP=0.
      RETURN
 29   TCON=T
      SCOMP=2*T*A0/E11(IA0)-DBT*A013*A013
      RETURN
      END
      SUBROUTINE FRAGTE(AA,ZZ,E00,K,KC,EPSIL0,TK,MU,L,SR)
C HAXO[[EH[E TEM[EPAT[P[ TK [ BCE[O M[[[T[[PA[MEHTHO[O COCTO[H[[ N(A),XZ
C [O [HEP[[[ BO[[[[[EH[[ (HA H[K[OH) E00 [[PA A0,Z0 [ [APAM.CB.O[[EMA K.
      REAL N,MU,K,LT,L,KC
      COMMON /BLOK1/W0,TC,RADNCL,G0 /BLOK0/A0,Z0,V0,XZ0,IA0,E0,E0Q
      COMMON /BLOK2/B0,BT,DBT,LT,CP,PKP013
      COMMON /BLOK3/E11(500),XZ(500),GA(500)
      COMMON /BLERR/ERRMI,ERRMA,ERRFR,ERRRE /BLJPRE/JPRE
      COMMON /BPERI/IPERI
      ERRMA=0.
      A0=AA
      Z0=ZZ
      IA0=INT(A0+0.5)
      SM=(1.+2.31*(E00-3.5))*A0/100.
      IF(SM.LT.2.) SM=2.
c --- keep constant or vary with multiplic. translational volume ---
      K=(1.+1.44*(SM**0.333333-1.)/(1.17*A0**0.333333))**3-1.
c --- changing of coulomb volume : V_c=V_0+V_f -----
c      KC=K
c ---
      E11(1)=0.
      DO 100 I=2,IA0
      A=I
      IF(IPERI.EQ.0) E11(I)=EPSIL0*(1.+3./(A-1))
C      IF(I.LE.50) E11(I)=EPSIL0*(1.+0.12*(A**0.33333-5.**0.33333))
C 100  IF(I.GT.50) E11(I)=EPSIL0*(1.+0.12*(50.**0.33333-5.**0.33333))
      IF(IPERI.EQ.1) E11(I)=EPSIL0/(1.+0.0020*(A/25.)**2)
 100  CONTINUE
      T=SQRT(E00/0.12)
      IF(JPRE.EQ.1) T=2.*E00/3.
      IF(T.LT.0.0012) T=0.0012
      HT=1.0
      IF(JPRE.EQ.1) HT=2.0
      ICNT=0
      K1=0
      K2=0
      ID=0
 24   CONTINUE
      IF(ICNT.GT.200) GO TO 71
      CALL FRAG(A0,Z0,T,K,KC,EPSIL0,MU,L,E0,SR)
C *********
c      !WRITE 100 110,E00,E0,T
c 110  FORMAT(2X,'**FRAGTE** E00=',F11.5,' E0=',F11.5,' T=',F11.5)
C *********
      D=(E00-E0)/E00
      IF(ABS(D).LT.0.003) GO TO 29
      ICNT=ICNT+1
      H=SIGN(HT,D)
      IF(D) 21,21,22
 21   K1=1
      IF(K1.EQ.1.AND.K2.EQ.1) ID=ID+1
      IF(ID.GT.30) GO TO 71
      T=T+H/2**ID
      IF(T.GE.0.001) GO TO 24
      K2=1
      H=HT
      GO TO 21
 22   K2=1
      IF(K1.EQ.1.AND.K2.EQ.1) ID=ID+1
      IF(ID.GT.30) GO TO 71
      T=T+H/2**ID
      IF(T.GE.0.001) GO TO 24
      K1=1
      H=HT
      GO TO 22
 29   TK=T
      RETURN
 71   ERRMA=ERRMA+1
C      !WRITE 100 72,ICNT,D,T,ID
C 72   FORMAT(5X,'ERROR IN FRAGTE   ICNT=',I3,' D=',E12.5,' T=',F8.3,
C     *'             ID=',I2)
      TK=T
      RETURN
      END
      SUBROUTINE FINDT(E00,A0,Z0,PA,PAZ,IP,TK,ETP,ECOL,ECP,EEX,EB)
C HAXO[[EH[E [O IP B[[PAHH[M [PA[MEHTAM [X TEM[EPAT[P[ T [ [HEP[[[,
C C[[TA[,[TO [O[HA[ [HEP[[[ BO[[[[[EH[[ C[CTEM[ E0*A0 - [[KC[POBAHA.
      REAL N,LT
      COMMON /BLOK1/W0,TC,RADNCL,G0
      COMMON /BLOK2/B0,BT,DBT,LT,CP,PKP013
      COMMON /BLOK3/E11(500),XZ(500),GA(500)
      COMMON /BLOK20/BT2,DBT2 /BLPARA/EPSIL0,FKAPPA,FKACOL
      COMMON /BLERR/ERRMI,ERRMA,ERRFR,ERRRE /BLJPRE/JPRE
      DIMENSION PA(500),PAZ(500)
      ERRFR=0.
      A013=A0**0.333333
      E0=E00
      HT=0.5
      IF(JPRE.EQ.1) HT=1.0
      ICNT=0
      K1=0
      K2=0
      ID=0
      T=TK
      IF(T.LT.0.0012) T=0.0012
      CALL DELAM(A0,Z0,DLM0,DSH0,BAR)
 24   CONTINUE
      IF(ICNT.GT.200) GO TO 71
      EEX=0.
      EB=-DLM0
      ECP=0.
c      ETP=1.5*T*IP
      ETP=1.5*T*(IP-1)
      CALL PARAM2(T,EPSIL0)
      DO 107 I=1,IP
      CALL DELAM(PA(I),PAZ(I),DLM,DSH,BAR)
      EB=EB+DLM
      PA03=PA(I)**0.333333
      ECP=ECP+(0.6*1.44/RADNCL)*PAZ(I)*PAZ(I)/PA03
      IA=INT(PA(I)+0.1)
      IF(IA.LE.3) GO TO 107
      IF(JPRE.EQ.1) GO TO 107
c ---       finding surface energy
c before June 96:     ESUF=PA(I)*T*T*2.5*B0/(TC*TC*PA03)
c                     IF(DBT.EQ.0.) ESUF=0.
      ESUF=(BT2-T*DBT2-B0)*PA03*PA03
c ---
      EEX=EEX+PA(I)*T*T/E11(IA)+ESUF
      IF(IA.EQ.4) EEX=EEX-ESUF
 107  CONTINUE
      ECOL=(0.6*1.44/RADNCL)*Z0*Z0*PKP013/A013-ECP*PKP013
      E0TP=ETP+EEX+EB+ECOL
      E0TPA0=E0TP/A0
C ********
c      !WRITE 100 110,E0,E0TPA0,T,ETP,EEX,EB,ECOL
c 110  FORMAT(2X,'*FINDT* E0=',F11.5,' E0TP/A0=',F11.5,' T=',F11.5,3X,
c     *'ETP,EEX,EB,ECOL=',4F11.5)
C ********
      D=(E0-E0TP/A0)/E0
      IF(ABS(D).LT.0.003) GO TO 29
      ICNT=ICNT+1
      H=SIGN(HT,D)
      IF(D) 21,21,22
 21   K1=1
      IF(K1.EQ.1.AND.K2.EQ.1) ID=ID+1
      IF(ID.GT.30) GO TO 71
      T=T+H/2**ID
      IF(T.GE.0.001) GO TO 24
      K2=1
      H=HT
      GO TO 21
 22   K2=1
      IF(K1.EQ.1.AND.K2.EQ.1) ID=ID+1
      IF(ID.GT.30) GO TO 71
      T=T+H/2**ID
      IF(T.GE.0.001) GO TO 24
      K1=1
      H=HT
      GO TO 22
 29   TK=T
      RETURN
 71   ERRFR=ERRFR+1
C      !WRITE 100 72,ICNT,ID,D,T,IP
C 72   FORMAT(5X,'ERROR IN FINDT    ICNT=',I3,' ID=',I3,' D=',E12.5,
C     *' T=',F9.3,'  IP=',I3)
      RETURN
      END
      FUNCTION UEVA(A,T)
C O[PE[E[EH[E [HEP[[[ BO[[[[[EH[[ [PA[MEHTOB [EPE[ [X [C[APEH[EM [[[
C [EPM[-PA[BA[OM (UEVA B [[B, TEM[EPAT[PA T B M[B).
      COMMON /BLOK2/B0,BT,DBT,LT,CP,PKP013 /BLOK1/W0,TC,RADNCL,G0
      COMMON /BLOK3/E11(500),XZ(500),GA(500)
      COMMON /BLOK20/BT2,DBT2
      IA=INT(A+0.5)
c ---       finding surface energy
c before June 96:     ESUF=A*T*T*2.5*B0/(TC*TC*A**0.333333))
c                     IF(DBT.EQ.0.) ESUF=0.
      ESUF=(BT2-T*DBT2-B0)*A**0.666667
c ---
      UEVA=A*T*T/E11(IA)+ESUF
      IF(IA.EQ.4) UEVA=A*T*T/E11(IA)
      UEVA=0.001*UEVA
      RETURN
      END
      SUBROUTINE FRAG(AA,ZZ,T,K,KC,EPSIL0,MU,L,EE,SR)
C [PA[MEHTA[[[ [[PA A0,Z0 [P[ [AHHO[ TEM[EPAT[PE T [ [APAM.CB. O[[EMA K.
C [APAMETP[ E11(A) [A[A[TC[ PAHEE B /BLOK2/. HA B[XO[E: X[M.[OT. MU [ L,
C [HEP[[[ E0(M[B/N),[HTPO[[[ S (1/N),CPE[H.[[C[O [PA[MEHTOB SR; MHO[ECTB
C HOCT[ N(A) [ [AP[[OBOE PAC[PE[E[EH[E XZ(A) - HAXO[[TC[ B /BLOK2/.
      REAL N,MU,K,LT,L,KC
      COMMON /BLOK2/B0,BT,DBT,LT,CP,PKP013 /BLOKN/N(500) /BENTRO/SSR
      COMMON /BLOK1/W0,TC,RADNCL,G0 /BLOK0/A0,Z0,V0,XZ0,IA0,E0,E0Q
      COMMON /BLOK3/E11(500),XZ(500),GA(500) /BLOKE/E(500)
      A0=AA
      Z0=ZZ
      TK=T
      CALL PARAM(T,K,KC,EPSIL0)
      A013=A0**0.33333333
      A023=A013*A013
      R0=RADNCL*A013
      V0=(4.*3.1415926/3.)*R0*R0*R0
      P0=A0/V0
      R=R0*(1.+KC)**0.33333333
      PKP013=R0/R
      XZ0=Z0/A0
      E0Q=A0*(-W0+G0*(1.-2.*XZ0)**2)+B0*A023+0.6*1.44*Z0*Z0/R0
      P=P0/(1.+KC)
      CP=(0.6*1.44/RADNCL)*(1.-PKP013)
      LT=16.15/SQRT(TK)
      CALL XZA(TK,K,MU,L)
      CALL ENERGA(TK,K)
      CALL ENTROP(TK,K,S)
      SSR=S
      E0T=0.
      DO 10 IT=1,IA0
      E0T=E0T+N(IT)*E(IT)
   10 CONTINUE
      E0T=E0T+0.6*1.44*Z0*Z0/R
      E0=(E0T-E0Q)/A0
      EE=E0
      SR=0.
      DO 31 J=1,IA0
      SR= SR+N(J)
 31   CONTINUE
      RETURN
      END
      SUBROUTINE CALCMU(T,K,MU,L)
C B[[[C[EH[E MU=MU(T) [C[O[[[[[ [PABHEH[E  N(1)*1+...+N(A0)*A0=A0
      REAL N,MU,K,LT,L
      COMMON /BLOK1/W0,TC,RADNCL,G0 /BLOK0/A0,Z0,V0,XZ0,IA0,E0,E0Q
      COMMON /BLOK2/B0,BT,DBT,LT,CP,PKP013    /BLOKN/N(500)
      COMMON /BLOK3/E11(500),XZ(500),GA(500)
      COMMON /BLERR/ERRMI,ERRMA,ERRFR,ERRRE /BLJPRE/JPRE
      HM=4.
      IF(JPRE.EQ.1) HM=16.
      MU=-W0-T*T/E11(5)-L*XZ(5)+G0*(1.-2*XZ(5))**2
     *+(2./3.)*BT/1.71+5.*CP*XZ(5)*XZ(5)*2.92/3.-1.5*T/5.
      IF(JPRE.EQ.1) MU=-A0/2.
      IF(ABS(MU).LE.HM) MU=-HM-0.01
      KK=0
      K1=0
      K2=0
      ID=0
 24   CONTINUE
      CALL MULTIA(T,K,MU,L)
      A0T=0.
      DO 10 IT=1,IA0
      A0T=A0T+IT*N(IT)
 10   CONTINUE
C *********
c      !WRITE 100 110,A0,A0T,MU
c 110  FORMAT(2X,'**CALCMU** A0=',F7.2,' A0T=',F7.2,' MU=',F8.3)
C *********
      KK=KK+1
      D=(A0-A0T)/A0
      IF(ABS(D).LT.0.001) GO TO 25
      IF(KK.GT.400.OR.ID.GT.30) GO TO 61
      H=SIGN(HM,D)
      IF(D) 21,21,22
 21   K1=1
      IF(K1.EQ.1.AND.K2.EQ.1) ID=ID+1
      MU=MU+H/2**ID
      GO TO 24
 22   K2=1
      IF(K1.EQ.1.AND.K2.EQ.1) ID=ID+1
      MU=MU+H/2**ID
      GO TO 24
 25   CONTINUE
      RETURN
   61 ERRMA=ERRMA+1
C      !WRITE 100 62,KK,ID,D,MU,T
C   62 FORMAT(5X,'ERROR IN CALCMU  KK=',I3,' ID=',I2,' D=',E12.5,' MU=',
C     *F8.3,' T=',F8.3)
      RETURN
       END
c ----- subroutines: Multia, Energa, Entrop, Xza, Selecz ------
      SUBROUTINE MULTIA(T,K,MU,L)
c --- corrected 3-H/3-He ratio , May 97 (GSI) ---
C B[[[C[EH[E MHO[ECTBEHHOCTE[ [PA[MEHTOB - N(A)
      REAL N,MU,K,LT,L
      COMMON /BLOK1/W0,TC,RADNCL,G0 /BLOK0/A0,Z0,V0,XZ0,IA0,E0,E0Q
      COMMON /BLOK2/B0,BT,DBT,LT,CP,PKP013   /BLOKN/N(500)
      COMMON /BLOK3/E11(500),XZ(500),GA(500)
      COMMON /BLYPRE/YNPR,YPPR,YNEQ,YPEQ
      COMMON /BLOKY/YN,YP,YD,YT,YH,YA /BLOKXZ/XZN,XZP,XZD,XZT,XZH,XZAL
      COMMON /BLJPRE/JPRE
      YN=(2*K*V0/(LT*LT*LT))*EXP(MU/T)
      YP=(2*K*V0/(LT*LT*LT))*EXP((MU+L-CP)/T)
      IF(JPRE.EQ.0) YNEQ=YN
      IF(JPRE.EQ.0) YPEQ=YP
      IF(JPRE.EQ.1) YNPR=YN
      IF(JPRE.EQ.1) YPPR=YP
      N(1)=YN+YP
      YD=3*2*SQRT(2.)*K*V0/(LT**3)
      YD=YD*EXP((2.796+2*(MU+L*XZD)-CP*XZD*XZD*(2.**(5./3.)))/T)
cc old: Y3=4*3*SQRT(3.)*K*V0/(LT**3)
cc      Y3=Y3*EXP((9.224+3*(MU+L*XZ(3))-CP*XZ(3)*XZ(3)*(3.**(5./3.)))/T)
      YT=2*3*SQRT(3.)*K*V0/(LT**3)
      YT=YT*EXP((8.48+0.51+3*(MU+L*XZT)-CP*XZT*XZT*(3.**(5./3.)))/T)
      YH=2*3*SQRT(3.)*K*V0/(LT**3)
      YH=YH*EXP((7.72+2.04+3*(MU+L*XZH)-CP*XZH*XZH*(3.**(5./3.)))/T)
      YA=8*K*V0/(LT**3)
      IF(JPRE.EQ.0)
     *YA=YA*EXP((30.11+4*(MU+L*XZAL+T*T/E11(4))
     *-CP*XZAL*XZAL*(4.**(5./3.)))/T)
      IF(JPRE.EQ.1)
     *YA=YA*EXP((30.11+4*(MU+L*XZAL)-CP*XZAL*XZAL*(4.**(5./3.)))/T)
      N(2)=YD
      N(3)=YT+YH
      N(4)=YA
      DO 2 I=5,IA0
      N(I)=0.
      IF(JPRE.EQ.1) GO TO 2
      A=I
      A23=A**0.66666667
      VNE=(MU+L*XZ(I)+W0+T*T/E11(I)-G0*(1.-2*XZ(I))**2)*A
     *-BT*A23-CP*XZ(I)*XZ(I)*A*A23
      VNE=VNE/T
      IF(VNE.GT.30.0) GO TO 1
      VN=EXP(VNE)
      N(I)=(K*V0*SQRT(A)*A/(LT*LT*LT))*VN
      IF(N(I).LT.1.0E-30) N(I)=0.
      GO TO 2
 1    N(I)=999.
 2    CONTINUE
      RETURN
      END
      SUBROUTINE ENERGA(T,K)
c --- corrected 3-H/3-He ratio , May 97 (GSI) ---
C B[[[C[EH[E [HEP[[[ [PA[MEHTOB  E(A)
      REAL N,MU,K,LT,L
      COMMON /BLOK1/W0,TC,RADNCL,G0 /BLOK0/A0,Z0,V0,XZ0,IA0,E0,E0Q
      COMMON /BLOK2/B0,BT,DBT,LT,CP,PKP013    /BLOKE/E(500)
      COMMON /BLOK3/E11(500),XZ(500),GA(500)
      COMMON /BLJPRE/JPRE
      COMMON /BLOKY/YN,YP,YD,YT,YH,YA /BLOKXZ/XZN,XZP,XZD,XZT,XZH,XZAL
      E(1)=1.5*T+CP*XZ(1)
      E(2)=-2.796+CP*XZ(2)*XZ(2)*(2.**(5./3.))+1.5*T
cc old:  E(3)=-9.224+CP*XZ(3)*XZ(3)*(3.**(5./3.))+1.5*T
      ET=-8.48-0.51+CP*XZT*XZT*(3.**(5./3.))+1.5*T
      EH=-7.72-2.04+CP*XZH*XZH*(3.**(5./3.))+1.5*T
      E(3)=(YT*ET+YH*EH)/(YT+YH)
      E(4)=-30.11+CP*XZ(4)*XZ(4)*(4.**(5./3.))+1.5*T+4*T*T/E11(4)
      IF(JPRE.EQ.1) E(4)=E(4)-4*T*T/E11(4)
      DO 2 I=5,IA0
      E(I)=0.
      IF(JPRE.EQ.1) GO TO 2
      A=I
      A23=A**0.66666667
      EV=A*(T*T/E11(I)-W0+G0*(1.-2*XZ(I))**2)
      ES=(BT-T*DBT)*A23
      EC=CP*A23*A*XZ(I)*XZ(I)
      ET=1.5*T
      E(I)=EV+ES+EC+ET
 2    CONTINUE
      RETURN
      END
      SUBROUTINE ENTROP(T,K,S)
c --- corrected 3-H/3-He ratio , May 97 (GSI) ---
C B[[[C[EH[E [HTPO[[[ [[PA [OC[E [PA[MEHTA[[[ - S
      REAL N,MU,K,LT,L
      COMMON /BLOK1/W0,TC,RADNCL,G0 /BLOK0/A0,Z0,V0,XZ0,IA0,E0,E0Q
      COMMON /BLOK2/B0,BT,DBT,LT,CP,PKP013    /BLOKN/N(500)
      COMMON /BLOK3/E11(500),XZ(500),GA(500)
      COMMON /BLJPRE/JPRE
      COMMON /BLOKY/YN,YP,YD,YT,YH,YA
      S=0.
      IF(YN.GT.0.) S=YN*(2.5+ALOG(2*K*V0/(LT*LT*LT*YN)))
      IF(YP.GT.0.) S=S+YP*(2.5+ALOG(2*K*V0/(LT*LT*LT*YP)))
      IF(YD.GT.0.) S=S+YD*(2.5+ALOG(3*K*V0*2*SQRT(2.)/(YD*LT**3)))
cc old:  IF(Y3.GT.0.) S=S+Y3*(2.5+ALOG(4*K*V0*3*SQRT(3.)/(Y3*LT**3)))
      IF(YT.GT.0.) S=S+YT*(2.5+ALOG(2*K*V0*3*SQRT(3.)/(YT*LT**3)))
      IF(YH.GT.0.) S=S+YH*(2.5+ALOG(2*K*V0*3*SQRT(3.)/(YH*LT**3)))
      IF(YA.GT.0.AND.JPRE.EQ.0)
     &S=S+YA*(2.5+ALOG(8*K*V0/(YA*LT**3))+8*T/E11(4))
      IF(YA.GT.0.AND.JPRE.EQ.1)
     &S=S+YA*(2.5+ALOG(8*K*V0/(YA*LT**3)))
      IF(JPRE.EQ.1) RETURN
      DO 1 I=5,IA0
      A=I
      IF(N(I).LE.0.) GO TO 1
      SV=A*2*T/E11(I)
      SS=-DBT*A**0.6666667
      ST=2.5+ALOG(K*V0*SQRT(A)*A/(LT*LT*LT*N(I)))
      S=S+(SV+SS+ST)*N(I)
 1    CONTINUE
      RETURN
      END
      SUBROUTINE XZA(T,K,MU,L)
c --- corrected 3-H/3-He ratio , May 97 (GSI) ---
C B[[[C[EH[E XZ(A) - [O[[ Z B [PA[MEHTE A,[C[O[[[[[ [PABHEH[E
C                 XZ(1)*1*N(1)+...+XZ(A0)*A0*N(A0)=Z0  ([P[[[[[EHHO)
      REAL N,MU,K,LT,L
      COMMON /BLOK1/W0,TC,RADNCL,G0 /BLOK0/A0,Z0,V0,XZ0,IA0,E0,E0Q
      COMMON /BLOK2/B0,BT,DBT,LT,CP,PKP013     /BLOKN/N(500)
      COMMON /BLOK3/E11(500),XZ(500),GA(500)
      COMMON /BLERR/ERRMI,ERRMA,ERRFR,ERRRE  /BLJPRE/JPRE
      COMMON /BLOKY/YN,YP,YD,YT,YH,YA /BLOKXZ/XZN,XZP,XZD,XZT,XZH,XZAL
      DIMENSION C1(500),C2(500)
      K1=0
      K2=0
      KK=0
      ID=0
      HM=1.
      IF(JPRE.EQ.1) HM=2.
      DO 1 I=5,IA0
      A=I
      CC=8.*G0+2.*CP*A**0.66666667
      C1(I)=4.*G0/CC
      C2(I)=1./CC
 1    CONTINUE
      L=XZ0*(8.*G0+2.*CP*A0**0.6666667)-4.*G0
 6    CONTINUE
      XZN=0.
      XZP=1.
      XZD=0.5
      XZT=0.333333
      XZH=0.666667
      XZAL=0.5
      DO 2 I=5,IA0
 2    XZ(I)=C1(I)+L*C2(I)
      CALL CALCMU(T,K,MU,L)
      XZ(1)=YP/(YP+YN)
      XZ(2)=XZD
      XZ(3)=(YT+2.*YH)/(3.*YT+3*YH)
      XZ(4)=XZAL
      Z0T=0.
      DO 3 I=1,IA0
 3    Z0T=Z0T+XZ(I)*I*N(I)
C *********
c      !WRITE 100 110,Z0,Z0T,L
c 110  FORMAT(2X,'**XZA** Z0=',F7.2,' Z0T=',F7.2,' L=',F8.3)
C *********
      KK=KK+1
      D=(Z0-Z0T)/Z0
      IF(ABS(D).LT.0.002) GO TO 9
      IF(KK.GT.200.OR.ID.GT.30) GO TO 7
      H=SIGN(HM,D)
      IF(D) 4,4,5
 4    K1=1
      IF(K1.EQ.1.AND.K2.EQ.1) ID=ID+1
      L=L+H/2**ID
      GO TO 6
 5    K2=1
      IF(K1.EQ.1.AND.K2.EQ.1) ID=ID+1
      L=L+H/2**ID
      GO TO 6
 9    CONTINUE
      RETURN
 7    ERRMA=ERRMA+1
C      !WRITE 100 8,KK,ID,D,L,MU,T
C 8    FORMAT(5X,'ERROR IN XZA KK=',I2,' ID=',I2,' D=',E12.5,' L=',F8.3,
C     *' MU=',F8.3,' T=',F8.3)
      RETURN
      END


      SUBROUTINE SELECZ(T,L,IZ0,IP,PA,PAZ)
c --- corrected 3-H/3-He ratio , May 97 (GSI) ---
      COMMON /BLOK2/B0,BT,DBT,LT,CP,PKP013
      COMMON /BLYPRE/YNPR,YPPR,YNEQ,YPEQ
      COMMON /BLJPRE/JPRE
      COMMON /BLOKY/YN,YP,YD,YT,YH,YA
      REAL PA(500),PAZ(500),L,LT

      common / acceptEvent/Error_5

      logical Error_5

C O[PE[E[EH[E [AP[[OB IP [PA[MEHTOB: B MACC[BE PAZ(1-IP). MACC[ [PA[MEHT
C [AH[ B MACC[BE PA(1-IP),B [OP[[KE [[[BAH[[ A. [AP[[[ [PA[[C[ [[ HOPMA[
C PAC[PE[E[EH[[ C [EHTPOM ZM [ [[C[EPC[E[ DZ. OCTATOK PO[[[P[[A ( OT[[[[
C C[MM[ PAZ OT Z0 HE [O[EE 0,+1,-1 ) [O[AB[[[C[ B PA(1).

      Loop_3 = 0
      Loop_5 = 0

      Error_5 = .false.

      G0=25.
 5    ISZ=0
      Loop_5 = Loop_5 + 1

      if (loop_5 .gt. 1000) then
         write(*,'(A,2x,i6)') '****** Selecz: Loop_5 : ',Loop_5
         Error_5 = .true.
         return
      endif

c      if (Error_5) then
c         write(*,'(A,i5)') '##### IP = ', IP
c         write(*,'(A,i5)') '##### IZ0 = ', IZ0
c         write(*,'(A,f12.4)') '##### T = ', T
c         write(*,'(A,f12.4)') '##### L = ', L
c      endif

      DO 1 I=1,IP
      IA=INT(PA(I)+0.5)
      IF(IA-1) 2,2,6
 2    BR=RNDM(-1)
      PAZ(I)=0.
cc      IF(BR.GT.(YN/(YN+YP))) PAZ(I)=1.
      IF(JPRE.EQ.0.AND.BR.GT.(YNEQ/(YNEQ+YPEQ))) PAZ(I)=1.
      IF(JPRE.EQ.1.AND.BR.GT.(YNPR/(YNPR+YPPR))) PAZ(I)=1.
      ISZ=ISZ+INT(PAZ(I)+0.5)
      GO TO 1
 6    IF(IA-2) 7,7,8
 7    PAZ(I)=1.
      ISZ=ISZ+INT(PAZ(I)+0.5)
      GO TO 1
 8    IF(IA-3) 9,9,10
 9    BR=RNDM(-1)
      PAZ(I)=1.
      IF(BR.GT.(YT/(YT+YH))) PAZ(I)=2.
cc old: IF(BR.GT.0.5) PAZ(I)=2.
      ISZ=ISZ+INT(PAZ(I)+0.5)
      GO TO 1
 10   IF(IA-4) 11,11,3
 11   PAZ(I)=2.
      ISZ=ISZ+INT(PAZ(I)+0.5)
      GO TO 1
 3    CALL RANNOR(BR1,BR2)
      Loop_3 = Lopp_3 + 1

      if (loop_3 .gt. 1000) then
         write(*,'(A,2x,i6)') 'Selecz: Loop_3 : ',Loop_3
         return
c         STOP
      endif

      CC=8.*G0+2.*CP*PA(I)**0.666667
      ZM=PA(I)*(4.*G0+L)/CC
cc      IF(PA(I).GT.1.5.AND.PA(I).LT.4.5) ZM=0.5*PA(I)
      DZ=SQRT(PA(I)*T/CC)
      Z=BR1*DZ+ZM
      IZ=INT(Z+0.5)
      IF(IZ.GE.0.AND.IZ.LE.IA) GO TO 4
      Z=BR2*DZ+ZM
      IZ=INT(Z+0.5)
      IF(IZ.GE.0.AND.IZ.LE.IA) GO TO 4
      GO TO 3
 4    PAZ(I)=IZ
      ISZ=ISZ+INT(PAZ(I)+0.5)
 1    CONTINUE


c      if (Error_5) then
c         do I=1,IP
c            write(*,'(A,i5,2x,2f12.4)')
c     &           '##### I,PA(I),PAZ(I) = ',I,PA(I), PAZ(I)
c         end do
c      endif
!      if (Error_5) then
!         write(*,*) 'STOP in SELECZ'
!         STOP
!      endif

      DELZ=IZ0-ISZ
      IF(ABS(DELZ).LT.0.1) RETURN
      IF(ABS(DELZ).GT.1.1) GO TO 5
      IF(PA(1).LT.4.5) GO TO 5
      PAZ(1)=PAZ(1)+DELZ




      RETURN
      END

c -----------------------------------------

      SUBROUTINE RANNOR(A,B)
      Y=RNDM(-1)
      Z=RNDM(-1)
      X=6.283185*Z
      A1=SQRT(-2.0*ALOG(Y))
      A=A1*SIN(X)
      B=A1*COS(X)
      RETURN
      END

      SUBROUTINE SELECT(A0,Z0,N,T,L,IP,PA,PAZ)
C C[[[A[HOE PA[[[EH[E [[PA (A0,Z0) HA [PA[MEHT[,B COOTBETCTB[[ C [HK[[[[
C PAC[PE[E[EH[EM N(A) [ HOPMA[[H[M PAC[PE[E[EH[EM [PA[MEHTOB [O Z.
C HA B[XO[E: [[C[O [PA[M.-IP,[X MACC[ -PA(1-IP) /B [OP[[KE [[[BAH[[ A/,
C [ [AP[[[ -PAZ(1-IP).
      INTEGER NA(500)
      REAL N(500),L
      DIMENSION PA(500),PAZ(500)
      IA0=INT(A0+0.5)
      IZ0=INT(Z0+0.5)
 14   CONTINUE
      CALL SELECA(N,IA0,NA,IP)
C  [A[O[HEH[E MACC[BA PA ( [[ MACC[BA NA )
      IS=0
      DO 305 I=1,IP
      PA(I)=0.
 305  PAZ(I)=0.
      DO 306 IA=1,IA0
      IF(NA(IA).EQ.0) GO TO 306
      JIA=NA(IA)
      DO 307 J=1,JIA
 307  PA(IS+J)=IA
      IS=IS+JIA
 306  CONTINUE
C [EPEH[MEPA[[[ MACC[BA PA B [OP[[KE [[[BAH[[ A
      DO 10 J=1,IP
      PAMAX=0.
      DO 11 I=J,IP
      IF(PA(I)-PAMAX) 11,11,12
 12   IM=I
      PAMAX=PA(IM)
 11   CONTINUE
      PA(IM)=PA(J)
      PA(J)=PAMAX
 10   CONTINUE
      IMINZ=0
      IMAXZ=0
      DO 13 IJ=1,IP
      IMAXZ=IMAXZ+INT((PA(IJ)-0.1)*0.5)+1
 13   IMINZ=IMINZ+INT((PA(IJ)+0.1)*0.5)
      IF(PA(1).LT.5.AND.IMINZ.GE.IZ0) GO TO 14
      IF(PA(1).LT.5.AND.IMAXZ.LE.IZ0) GO TO 14
      CALL SELECZ(T,L,IZ0,IP,PA,PAZ)
      RETURN
      END

      SUBROUTINE SELECA(N,IA0,NA,IP)
c -- newest, taking into account that multiplicity distribution should
c be a Poisson in grand canonics. 24 Nov. 99, Bo.
C C[[[A[H[[ B[[OP [[C[A [PA[MEHTOB IP [ [X MHO[ECTBEHHOCTE[ [O A -NA(A)
C [[ [[C[A H[K[OHOB IA0 [ [HK[[[[BH[X MHO[ECTBEHHOCTE[ -N(A).
      REAL WPOIS(500)
      REAL GG(500),N(500)
      INTEGER NA(500)
c      IMAXMA=0
      IMAXMA=1
      IMAXMI=0
      GH=0.
      DO 2 I=1,IA0
 2    GH=GH+N(I)
c      SQGH=SQRT(GH)+0.5
cc      SQGH=SQRT(GH)+1.0
c      SQGH=2.*SQRT(GH)+1.0
c ------------------------------------
      WPOIS(1)=GH*EXP(-GH)
      DO 11 I=2,IA0
      WPOIS(I)=WPOIS(I-1)*GH/FLOAT(I)
 11   CONTINUE
      MMAX=INT(GH)
      WPMAX=WPOIS(MMAX)
c ------------------------------------
      GG(1)=N(1)
      DO 3 I=2,IA0
 3    GG(I)=GG(I-1)+N(I)
 4    ISA=0
      DO 1 I=1,IA0
 1    NA(I)=0
      IP=0
 5    BRND=RNDM(-1)*GH
      DO 6 I=1,IA0
      IF(BRND-GG(I)) 7,6,6
 7    IA=I
      GO TO 8
 6    CONTINUE
 8    IP=IP+1
      NA(IA)=NA(IA)+1
      ISA=ISA+IA
      IAK=IA0-ISA
      IF(IAK.EQ.0) GO TO 10
      IF(IAK.GT.IMAXMA) GO TO 5
      IF(IAK.LT.IMAXMI) GO TO 4
      IP=IP+1
      NA(IAK)=NA(IAK)+1
 10   CONTINUE
c      IF(ABS(GH-IP).GT.SQGH) GO TO 4
c ------------------------------------
      BRNP=RNDM(-1)*WPMAX
      IF(BRNP.GT.WPOIS(IP)) GO TO 4
c ------------------------------------
c      IF(IP.EQ.1) GO TO 4
      IF(IP.LE.2) GO TO 4
      RETURN
      END

      SUBROUTINE DISNET(INET,IP0,AMP,PN0,T,TN)
C HAXO[[EH[E [HEP[[[ [ [M[[[[COB INET HE[TPOHOB [O [X TEM[EPAT[PE T(MEV)
C [O[HO[ K[H.[HEP[[[ TN(GEV),C[[TA[ PAC[PE[E[EH[E MAKCBE[OBCK[M. AMP(IP0
C -[X MACC[ (B GEV);[M[[[[C[ (GEV/C),B C[MME [A[[[E 0,[O[AB[[[TC[ B MACC
C PN0(3,500) C (IP0+1).
      DIMENSION AMP(500),PN0(3,500),ANL(3),PX(500),PY(500),PZ(500)
      DIMENSION PI(3),PJ(3),PN(3),AR(3),BR(3)
      IF(INET.LE.0) RETURN
      IP01=IP0+1
      IP0M=IP0+INET
      IF(INET-1) 1,1,2
 1    P=SQRT(2.*AMP(IP01)*TN)
      CALL ISOTR(ANL)
      PN0(1,IP01)=P*ANL(1)
      PN0(2,IP01)=P*ANL(2)
      PN0(3,IP01)=P*ANL(3)
      RETURN
 2    IF(INET-2) 3,3,4
 3    P=SQRT(2.*(AMP(IP01)*AMP(IP0M)/(AMP(IP01)+AMP(IP0M)))*TN)
      CALL ISOTR(ANL)
      PN0(1,IP01)=P*ANL(1)
      PN0(2,IP01)=P*ANL(2)
      PN0(3,IP01)=P*ANL(3)
      PN0(1,IP0M)=-PN0(1,IP01)
      PN0(2,IP0M)=-PN0(2,IP01)
      PN0(3,IP0M)=-PN0(3,IP01)
      RETURN
 4    IP0M2=IP0M-2
      ES=0.
      PSX=0.
      PSY=0.
      PSZ=0.
      FEMT=SQRT(0.5*T)*EXP(-0.5)
      DO 5 I=IP01,IP0M2
 7    E=RNDM(-1)*9.*T
      FE=SQRT(E)*EXP(-E/T)
      FERAND=RNDM(-1)*FEMT
      IF(FERAND-FE) 6,6,7
 6    P=SQRT(2.*E*0.001*AMP(I))
      CALL ISOTR(ANL)
      PX(I)=P*ANL(1)
      PY(I)=P*ANL(2)
      PZ(I)=P*ANL(3)
      ES=ES+E*0.001
      PSX=PSX+PX(I)
      PSY=PSY+PY(I)
      PSZ=PSZ+PZ(I)
 5    CONTINUE
      I1=IP0M2+1
      I2=IP0M2+2
      PPX=-PSX
      PPY=-PSY
      PPZ=-PSZ
      P=SQRT(PPX*PPX+PPY*PPY+PPZ*PPZ)
      EE=TN-ES
      EM=P*P/(2.*(AMP(I1)+AMP(I2)))
      IF(EE.LE.EM) GO TO 4
      H=1.+AMP(I2)/AMP(I1)
      CTM12=H*(1.-2.*AMP(I2)*EE/(P*P))
 11   CT1=1.-2.*RNDM(-1)
      IF(CT1*CT1-CTM12) 11,11,12
 12   IF(CTM12) 13,17,17
 13   IZN=1
      GO TO 15
 17   IF(CT1) 11,14,14
 14   IF(RNDM(-1)-0.5) 16,16,13
 16   IZN=-1
 15   CONTINUE
      P1=(P*CT1+IZN*SQRT(P*P*CT1*CT1-P*P*CTM12))/H
      P2=SQRT(P1*P1+P*P-2.*P1*P*CT1)
      PHI=6.28318*RNDM(-1)
      ST1=SQRT(1.-CT1*CT1)
      CPHI1=COS(PHI)
      SPHI1=SIN(PHI)
      CPHI2=-CPHI1
      SPHI2=-SPHI1
      CT2=(P*P+P2*P2-P1*P1)/(2.*P*P2)
      IF(CT2.GT.-1.AND.CT2.LT.1) GO TO 20
      ST2=0.
      GO TO 21
 20   ST2=SQRT(1.-CT2*CT2)
 21   CONTINUE
      PI(1)=P1*ST1*CPHI1
      PI(2)=P1*ST1*SPHI1
      PI(3)=P1*CT1
      PJ(1)=P2*ST2*CPHI2
      PJ(2)=P2*ST2*SPHI2
      PJ(3)=P2*CT2
      AR(1)=PPX
      AR(2)=PPY
      AR(3)=PPZ
      BR(1)=1.
      BR(2)=0.
      BR(3)=0.
      CALL ROTOR(AR,BR,PI,PN)
      PX(I1)=PN(1)
      PY(I1)=PN(2)
      PZ(I1)=PN(3)
      CALL ROTOR(AR,BR,PJ,PN)
      PX(I2)=PN(1)
      PY(I2)=PN(2)
      PZ(I2)=PN(3)
      PSX=PSX+PX(I1)+PX(I2)
      PSY=PSY+PY(I1)+PY(I2)
      PSZ=PSZ+PZ(I1)+PZ(I2)
      ES=ES+(PX(I1)**2+PY(I1)**2+PZ(I1)**2)/(2.*AMP(I1))
     *+(PX(I2)**2+PY(I2)**2+PZ(I2)**2)/(2.*AMP(I2))
      DO 8 I=IP01,IP0M
      PN0(1,I)=PX(I)
      PN0(2,I)=PY(I)
      PN0(3,I)=PZ(I)
 8    CONTINUE
      RETURN
      END

      SUBROUTINE ROTOR (AR,BR,PSTAR,PR)
C--------[[OK [OBOPOTA BEKTOPA PSTAR [[ PA[O[E[ C.[.[. B [CXO[H[[ (PR)
C--------AR, BR-BEKTOPA B [CXO[HO[ C.[.[.
C    BLOCK OF ROTATION.
      DIMENSION AR(3),BR(3),PSTAR(3),PR(3),AN(3)
      SP = 0.
      DO 31 IR=1,3
      SP = SP+AR(IR)*BR(IR)
   31 CONTINUE
      AMOD = SQRT(AR(1)**2+AR(2)**2+AR(3)**2)
      ALPHA1 = SP/AMOD
      BMOD2 = BR(1)**2+BR(2)**2+BR(3)**2
      ALPHA2 = SQRT(BMOD2-ALPHA1**2)
      AN(1) = AR(2)*BR(3)-AR(3)*BR(2)
      AN(2) = AR(3)*BR(1)-AR(1)*BR(3)
      AN(3) = AR(1)*BR(2)-AR(2)*BR(1)
      PR(1)=PSTAR(1)*BR(1)/ALPHA2+(PSTAR(3)-ALPHA1*PSTAR(1)/ALPHA2)
     1*AR(1)/AMOD+(PSTAR(2)*AN(1))/(ALPHA2*AMOD)
      PR(2)=PSTAR(1)*BR(2)/ALPHA2+(PSTAR(3)-ALPHA1*PSTAR(1)/ALPHA2)
     1*AR(2)/AMOD+(PSTAR(2)*AN(2))/(ALPHA2*AMOD)
      PR(3)=PSTAR(1)*BR(3)/ALPHA2+(PSTAR(3)-ALPHA1*PSTAR(1)/ALPHA2)
     1*AR(3)/AMOD+(PSTAR(2)*AN(3))/(ALPHA2*AMOD)
      RETURN
      END

      SUBROUTINE CULIMP(IP0,A0,ZC0,AM,CH,PN0,TEMP,ECOLB)
      DIMENSION X0(500),Y0(500),Z0(500),VX0(500),VY0(500),VZ0(500),
     *VX(500),VY(500),VZ(500),AM(500),CH(500),VN(500),RS(500),
     *PN0(3,500),PNT(3,500)
      COMMON /BLPAT/PA(500),PAZ(500),IP  /BPLACE/RFI
      COMMON /BLCEN2/XC(500),YC(500),ZC(500)  /BLDT/DT,IT
c     DT=10.
      DT=2.
      RN=1.17
      RN0=RN
      IF(RFI.GT.0.001) RN0=RFI
      RSYS=2.0*RN*(A0**0.333333)
      IP=IP0
      DO 50 I=1,IP
      PA(I)=AM(I)/0.94
 50   PAZ(I)=CH(I)
      CALL PLACE2(RSYS,RN0)
      DO 1 I=1,IP
      X0(I)=XC(I)
      Y0(I)=YC(I)
      Z0(I)=ZC(I)
 1    CONTINUE
      ET=1.5*IP*TEMP*0.001
c      CALL DISNET(IP,0,AM,PNT,TEMP,ET)
      CALL DISNE2(IP,0,AM,PNT,TEMP,ET)
      DO 40 I=1,IP
      VX0(I)=PNT(1,I)/AM(I)
      VY0(I)=PNT(2,I)/AM(I)
      VZ0(I)=PNT(3,I)/AM(I)
 40   CONTINUE
      CALL CULON(X0,Y0,Z0,VX0,VY0,VZ0,VX,VY,VZ,AM,CH,VN,RS,IP,ECOLB,ET)
      DO 7 I=1,IP
      PN0(1,I)=AM(I)*VX(I)
      PN0(2,I)=AM(I)*VY(I)
      PN0(3,I)=AM(I)*VZ(I)
 7    CONTINUE
      RETURN
      END

      SUBROUTINE CULON(X0,Y0,Z0,VX0,VY0,VZ0,VX,VY,VZ,AM,CH,VN,RS,M,EC,
     *ET)
      DIMENSION X0(500),Y0(500),Z0(500),VX0(500),VY0(500),VZ0(500),
     *X(500),Y(500),Z(500),VX(500),VY(500),VZ(500),AM(500),CH(500),
     *R(500,500),F(500,500),FX(500,500),FY(500,500),FZ(500,500),
     *FSX(500),FSY(500),FSZ(500),AX(500),AY(500),AZ(500),VS(500),
     *VN(500),RS(500),VSX(500),VSY(500),VSZ(500)
      COMMON /BLDT/DT,IT
c ---- following Coulomb acceleration in time for fragm. ----
      COMMON /BLACCL/ACL(50,500),ITN(50),JJN,NSIGN,ZSIGN(500)
      JJN=0
      JTOLD=0
      NSIGN=0
c -----------------
      IT=0
      TN=0.0
      TS=0.0
      DO 9 I=1,M
      AX(I)=0.0
      AY(I)=0.0
      AZ(I)=0.0
      VS(I)=SQRT(VX0(I)*VX0(I)+VY0(I)*VY0(I)+VZ0(I)*VZ0(I))
      VX(I)=VX0(I)
      VY(I)=VY0(I)
      VZ(I)=VZ0(I)
      X(I)=X0(I)
      Y(I)=Y0(I)
      Z(I)=Z0(I)
      DO 9 J=1,M
      F(I,J)=0.0
      R(I,J)=0.0
      FX(I,J)=0.0
      FY(I,J)=0.0
      FZ(I,J)=0.0
 9    CONTINUE
      IF(EC.LE.0.) RETURN
      NSIGN=M
 10   CONTINUE
      DO 1 I=1,M
      DO 1 J=1,M
      IF(I.EQ.J) GO TO 1
      R(I,J)=SQRT((X(I)-X(J))*(X(I)-X(J))+(Y(I)-Y(J))*(Y(I)-Y(J))+
     *(Z(I)-Z(J))*(Z(I)-Z(J)))
 1    CONTINUE
      ECULT=0.0
      DO 2 I=1,M
      DO 2 J=1,M
      IF(I.EQ.J) GO TO 2
      F(I,J)=1.44*(CH(I)*CH(J))/(R(I,J)*R(I,J))
      ECULT=ECULT+(1.44*0.5)*(CH(I)*CH(J))/R(I,J)
 2    CONTINUE
      DO 3 I=1,M
      DO 3 J=1,M
      IF(I.EQ.J) GO TO 3
      FX(I,J)=F(I,J)*(X(I)-X(J))/R(I,J)
      FY(I,J)=F(I,J)*(Y(I)-Y(J))/R(I,J)
      FZ(I,J)=F(I,J)*(Z(I)-Z(J))/R(I,J)
 3    CONTINUE
      DO 5 I=1,M
      FSX(I)=0.0
      FSY(I)=0.0
 5    FSZ(I)=0.0
      DO 4 I=1,M
      DO 4 J=1,M
      IF(I.EQ.J) GO TO 4
      FSX(I)=FSX(I)+FX(I,J)
      FSY(I)=FSY(I)+FY(I,J)
      FSZ(I)=FSZ(I)+FZ(I,J)
 4    CONTINUE
C     !WRITE 100 21,((FSX(I),FSY(I),FSZ(I)),I=1,M)
C21   FORMAT(2X,'FXYZ=',9(E11.4,1X))
c ----------------
      JJO=JJN
      JTNEW=TN
      IF((JTOLD+10).GT.JTNEW) GO TO 40
      JTOLD=JTNEW
      JJN=JJN+1
      IF(JJN.LE.50) ITN(JJN)=JTNEW
 40   CONTINUE
c ----------------
      DO 6 I=1,M
      AX(I)=FSX(I)/(AM(I)*1000.0)
      AY(I)=FSY(I)/(AM(I)*1000.0)
      AZ(I)=FSZ(I)/(AM(I)*1000.0)
c ----------------
      IF(JJN.LE.JJO.OR.JJN.GT.50) GO TO 41
      ZSIGN(I)=CH(I)
      ACL(JJN,I)=SQRT(AX(I)**2+AY(I)**2+AZ(I)**2)
 41   CONTINUE
c ----------------
 6    CONTINUE
      TN=TS+DT
      IT=IT+1
      DO 7 I=1,M
      VSX(I)=VX(I)
      VSY(I)=VY(I)
 7    VSZ(I)=VZ(I)
      DO 8 I=1,M
      VX(I)=VSX(I)+AX(I)*(TN-TS)
      VY(I)=VSY(I)+AY(I)*(TN-TS)
      VZ(I)=VSZ(I)+AZ(I)*(TN-TS)
      VS(I)=SQRT(VSX(I)*VSX(I)+VSY(I)*VSY(I)+VSZ(I)*VSZ(I))
      VN(I)=SQRT(VX(I)*VX(I)+VY(I)*VY(I)+VZ(I)*VZ(I))
      X(I)=X(I)+(VSX(I)+VX(I))*(TN-TS)*0.5
      Y(I)=Y(I)+(VSY(I)+VY(I))*(TN-TS)*0.5
      Z(I)=Z(I)+(VSZ(I)+VZ(I))*(TN-TS)*0.5
 8    RS(I)=SQRT(X(I)*X(I)+Y(I)*Y(I)+Z(I)*Z(I))
C     !WRITE 100 20,(VN(I),I=1,M),ECULT,(RS(I),I=1,M)
C     !WRITE 100 23,IT
c     IF(IT.GT.50) GO TO 11
      IF(IT.GT.50) DT=4.
      IF(IT.GT.75) DT=10.
      IF(IT.GT.125) DT=15.
      IF(IT.GT.150) DT=20
      IF(IT.GT.200) GO TO 11
      TS=TN
      GO TO 10
 11   CONTINUE
C20   FORMAT(1X,' VN=',3(1X,E11.4),3X,'ECULT=',E11.4,
C    *' RS=',3(1X,E11.4))
C23   FORMAT(2X,'IT=',I5)
      EE=0.
      DO 30 I=1,M
 30   EE=EE+(AM(I)*0.5)*VN(I)*VN(I)
      ETA=SQRT((EC+ET)/EE)
      DO 31 I=1,M
      VX(I)=ETA*VX(I)
      VY(I)=ETA*VY(I)
      VZ(I)=ETA*VZ(I)
 31   CONTINUE
      RETURN
      END

      SUBROUTINE PLACE(RSYS,RN)
      COMMON /BLPAT/ PA(500),PAZ(500),IP
      COMMON /BLCEN/ XC(500),YC(500),ZC(500)
      COMMON /CELIPS/CEL,IPOSF1
      JJJJ=0
C FOR PROLONGED ELLIPSOID in X-axis CELX>1
      CEL=1
c       CELY=2.
c      CELX=0.5
c      IF(IPOSF1.EQ.1) RSYS=1.26*RSYS
C FINDING FRAGMENT POSITIONS IN THE SYSTEM
 1    I=1
      R=(RSYS-RN*PA(I)**0.3333333)*(RNDM(-1))**0.3333333
      IF(IPOSF1.EQ.1) R=0.
      CT=1.-2.*RNDM(-1)
      PHI=6.28318*RNDM(-1)
      ST=SQRT(1-CT**2)
      CPHI=COS(PHI)
      SPHI=SIN(PHI)
      XC(1)=R*ST*CPHI
      YC(1)=R*ST*SPHI
C      ZC(1)=R*CT
      ZC(1)=R*CT*CEL
c ------------
c      XC(1)=R*CT*CELX
c      YC(1)=R*ST*CPHI
c      ZC(1)=R*ST*SPHI
c      XC(1)=R*ST*SPHI
c      YC(1)=R*CT*CELY
c      ZC(1)=R*ST*CPHI
c ------------
      IF(IP.LE.1) RETURN
 3    K=I
      KK=0
      I=I+1
 4    CONTINUE
      KK=KK+1
      IF(KK.GT.1000) GO TO 2
      R=(RSYS-RN*PA(I)**0.3333333)*(RNDM(-1))**0.3333333
      CT=1.-2.*RNDM(-1)
      PHI=6.28318*RNDM(-1)
      ST=SQRT(1-CT**2)
      CPHI=COS(PHI)
      SPHI=SIN(PHI)
      XC(I)=R*ST*CPHI
      YC(I)=R*ST*SPHI
C      ZC(I)=R*CT
      ZC(I)=R*CT*CEL
c      XC(I)=R*CT*CELX
c      YC(I)=R*ST*CPHI
c      ZC(I)=R*ST*SPHI
c      XC(I)=R*ST*SPHI
c      YC(I)=R*CT*CELY
c      ZC(I)=R*ST*CPHI
      DO 5 J=1,K
      RR=(XC(I)-XC(J))**2+(YC(I)-YC(J))**2+(ZC(I)-ZC(J))**2
      RMIN=RN*PA(I)**0.3333333+RN*PA(J)**0.3333333
      IF(RR-RMIN*RMIN) 4,5,5
 5    CONTINUE
      IF(I.EQ.IP) RETURN
      GO TO 3
 2    CONTINUE
      JJJJ=JJJJ+1
      IF(JJJJ.GE.1000) then
         WRITE(*,10) KK,I,PA(I),IP
         STOP
      endif
 10   FORMAT(3X,'ERROR IN PLACE - KK=',I4,' I=',I3,' PA(I)=',F5.1,
     &     ' IP=',I3)
      IF(IPOSF1.EQ.1) RSYS=1.1*RSYS
      GO TO 1
      END

      SUBROUTINE PLACE2(RSYS,RN)
      COMMON /BLPAT/ PA(500),PAZ(500),IP
      COMMON /BLCEN2/ XC(500),YC(500),ZC(500)
C FOR PROLONGED ELLIPSOID CEL>1
      CEL=1.
C FINDING FRAGMENT POSITIONS IN THE SYSTEM
 1    I=1
      R=(RSYS-RN*PA(I)**0.3333333)*(RNDM(-1))**0.3333333
      CT=1.-2.*RNDM(-1)
      PHI=6.28318*RNDM(-1)
      ST=SQRT(1-CT**2)
      CPHI=COS(PHI)
      SPHI=SIN(PHI)
      XC(1)=R*ST*CPHI
      YC(1)=R*ST*SPHI
C      ZC(1)=R*CT
      ZC(1)=R*CT*CEL
 3    K=I
      KK=0
      I=I+1
 4    CONTINUE
      KK=KK+1
      IF(KK.GT.1000) GO TO 2
      R=(RSYS-RN*PA(I)**0.3333333)*(RNDM(-1))**0.3333333
      CT=1.-2.*RNDM(-1)
      PHI=6.28318*RNDM(-1)
      ST=SQRT(1-CT**2)
      CPHI=COS(PHI)
      SPHI=SIN(PHI)
      XC(I)=R*ST*CPHI
      YC(I)=R*ST*SPHI
C      ZC(I)=R*CT
      ZC(I)=R*CT*CEL
      DO 5 J=1,K
      RR=(XC(I)-XC(J))**2+(YC(I)-YC(J))**2+(ZC(I)-ZC(J))**2
      RMIN=RN*PA(I)**0.3333333+RN*PA(J)**0.3333333
      IF(RR-RMIN*RMIN) 4,5,5
 5    CONTINUE
      IF(I.EQ.IP) RETURN
      GO TO 3
C 2    CONTINUE
 2    WRITE(*,10) KK,I,PA(I),IP
      STOP
 10   FORMAT(3X,'ERROR IN PLACE2 - KK=',I4,' I=',I3,' PA(I)=',F5.1,
     *     ' IP=',I3)
      GO TO 1
      END
      SUBROUTINE RAZVAL(UP,AP,ZP,PN,KST)
C [EPM[ PA[BA[ [E[KO[O [[PA C A=AP [ Z=ZP. [HEP[[[ BO[[[[[EH[[-UP ([[B).
C HA[. [M[[[[C-PN(3) ([[B/C). KST-HOMEP [A[O[HEH[[ MACC[BA SPT.
C ******** HOBA[ BEPC[[ - [ACT[[[ [A[A[TC[ B TA[[[[E *********
C          ( [[[T[BA[TC[ [ACT[[[ [O A=16 )
      COMMON /BLOK79/SFR(6,99) /ILEVRA/ILEVRA
      COMMON /SPDM/ MS(12,11),DM(12,11) /FKAP/FKAP /BLOKC/SPTU(10,500)
      COMMON /WWW/WK(20),WA(20),WN(100) /BLENN/ENN(100) /SUMN/NK(20),NF
      COMMON /BLJ/J(20),IEND /BLIAN/IAN(100) /BLIZN/IZN(100)
      COMMON /BLSECON/LSECON
      DIMENSION PN(3),VN(3),PNC(3),PNL(3),PN0(3,500),AMK(500)
      DIMENSION IPA(20),IPZ(20),WCON(2000),JCON(13000)
      KSTOLD=KST
      ILEVRA=2
      FKAP=1.
      IA0=INT(AP+0.5)
      IZ0=INT(ZP+0.5)
      NF=0
      DO 69 I=1,100
 69   WN(I)=0.
      DO 70 I=1,20
      NK(I)=0
      WA(I)=0.
 70   WK(I)=0.
      IF(AP.LT.1.9) GO TO 37
      PQ=PN(1)**2+PN(2)**2+PN(3)**2
      EAP=SQRT(PQ+(0.94*AP)**2)
      DO 40 I=1,3
 40   VN(I)=PN(I)/EAP
      CALL HELPRA
      DAZ=UP+GMASS(IA0,IZ0)
      SRQ=0.
      K2=0
      K00=IA0
c      IF(K00.GT.6) K00=6
      DO 71 K=2,K00
      J(1)=0
      DO 72 I=2,K
 72   J(I)=1
 78   CONTINUE
      CALL DISNUM(IA0,IZ0,K)
c      IF(IEND.EQ.1.AND.IA0.LE.4) !WRITE 100 7,(J(I),I=1,K)
c 7    FORMAT(2X,'**** IEND=1 ****  J=',12I5)
      IF(IEND.EQ.1) GO TO 71
      IF(NF.GE.2000.OR.K2.GT.12980) GO TO 71
      NF=NF+1
      NK(K)=NK(K)+1
      TQ=BCUL(K)
      RMQ=WEIGHT(IA0,K,DAZ,TQ)
c      !WRITE 100 300,NF,IA0,IZ0,DAZ,TQ,RMQ
c 300  FORMAT(1X,'NF,IA0,IZ0=',3I4,' DAZ,TQ,RMQ=',2F10.4,E12.4)
c      !WRITE 100 3,(J(I),I=1,K)
c 3    FORMAT(1X,'selec.config.J=',15I4)
      SRQ=SRQ+RMQ
      WCON(NF)=RMQ
      DO 73 I=1,K
      JI=J(I)
      JCON(K2+I)=JI
      IA=IAN(JI)
      IZ=IZN(JI)
      WN(JI)=WN(JI)+RMQ
 73   WA(IA)=WA(IA)+RMQ
      K2=K2+K
      WK(K)=WK(K)+RMQ
      GO TO 78
 71   CONTINUE
      IF(SRQ.LE.0.) GO TO 37
      DO 33 I=1,IA0
      WK(I)=WK(I)/SRQ
 33   WA(I)=WA(I)/SRQ
      DO 34 I=1,100
 34   WN(I)=WN(I)/SRQ
      DO 32 I=1,NF
 32   WCON(I)= WCON(I)/SRQ
 210                 BR=RNDM(-1)
C     !WRITE 100 106,K2,BR
C106  FORMAT(2X,'[[C[O BCEX [ACT[[ B MACC[BE K2=',I5,' C[[[. BR=',F8.5)
      SR=0.
      DO 13 N=1,NF
      NB=N
      SR=SR+WCON(N)
      IF(BR.LT.SR) GO TO 14
 13   CONTINUE
   14 CONTINUE
      NNP=0
      NN=0
      DO 61 I=2,K00
      NNP=NNP+NK(I)*I
      NN=NN+NK(I)
      IF(NN.GE.NB) GO TO 62
 61   CONTINUE
 62   K=I
      NPK1=(NB-(NN-NK(K))-1)*K
      NPK0=(NNP-NK(K)*K)+NPK1
c      !WRITE 100 101,NN,NB,NPK1,NPK0,K
c101   FORMAT(2X,'NN=',I5,' NB=',I5,' NPK1=',I5,' NPK0=',I5,' K=',I5)
      TN=DAZ
      DO 63 I=1,K
      JI=JCON(NPK0+I)
c      !WRITE 100 105,JI
c 105  FORMAT(5X,'JI=',I5)
      IPA(I)=IAN(JI)
      IPZ(I)=IZN(JI)
      AMK(I)=GMASS(IPA(I),IPZ(I))+0.001*ENN(JI)
      TN=TN-AMK(I)
 63   CONTINUE
      IF(TN.LT.0.) GO TO 210
c      !WRITE 100 102,(IPA(I),IPZ(I),I=1,K)
c 102  FORMAT(2X,'IPA,IPZ=',6(2I5,5X))
c      !WRITE 100 103,UP,DAZ,TN,(AMK(I),I=1,K)
c 103  FORMAT(2X,'UP,DAZ,TN=',3F8.4,' AMK=',8F8.4)
      CALL DISIMP(K,AMK,PN0,TN)
      DO 35 I=1,K
      DO 36 II=1,3
 36   PNC(II)=PN0(II,I)
      EN=SQRT(PNC(1)**2+PNC(2)**2+PNC(3)**2+AMK(I)**2)
      CALL CLPV(PNC,VN,PNL,EN)
      CALL PINT(PNL,CT,ST,CF,SF,TK,AMK(I))
C     !WRITE 100 104,AMK(I),TK,(PNL(II),II=1,3)
C104  FORMAT(2X,'AMK=',F8.4,' TK=',F8.4,' PNL=',3F8.4)
      SPTU(4,KST)=CT
      SPTU(5,KST)=SF
      SPTU(6,KST)=CF
      SPTU(7,KST)=TK
      SPTU(8,KST)=IPZ(I)
      SPTU(9,KST)=0.94*INT(AMK(I)/0.94+0.5)
      KST=KST+1
 35   CONTINUE
      GO TO 200
C      RETURN
 37   CONTINUE
      IF(AP.LT.0.99) RETURN
      WK(1)=1.
      WA(IA0)=1.
      NF=1
      NK(1)=1
      DO 38 II=1,100
      IA1=IAN(II)
      IZ1=IZN(II)
      IF(IA0.EQ.IA1.AND.IZ0.EQ.IZ1) WN(II)=WN(II)+1.
      IF(IA0.EQ.IA1.AND.IZ0.EQ.IZ1) GO TO 39
 38   CONTINUE
 39   CONTINUE
      CMN=0.94*AP
      CALL PINT(PN,CT,ST,CF,SF,TMN,CMN)
      SPTU(4,KST)=CT
      SPTU(5,KST)=SF
      SPTU(6,KST)=CF
      SPTU(7,KST)=TMN
      SPTU(8,KST)=ZP
      SPTU(9,KST)=CMN
      IF(LSECON.GT.0) SPTU(2,KST)=1.
      KST=KST+1
c      GO TO 200
      RETURN
 200  KST1=KST-1
C      !WRITE 100 203,KSTOLD,KST1
C 203  FORMAT(5X,'KSTOLD=',I3,' KST1=',I3,' ***** SPTU (BEFORE) *****')
C      !WRITE 100 204,((SPTU(I1,I2),I1=4,9),I2=KSTOLD,KST1)
C 204  FORMAT(2X,6F10.4)
      CALL POVRAZ(KSTOLD,KST1,KST2)
      KSS=KSTOLD-1
      DO 201 IPOV=1,KST2
      KSS=KSS+1
      DO 201 JPOV=1,6
      JPOVN=JPOV+3
      SPTU(JPOVN,KSS)=SFR(JPOV,IPOV)
      SPTU(2,KSS)=1.
 201  CONTINUE
      KST=KSS+1
C      !WRITE 100 205,KST2
C 205  FORMAT(5X,' KST2=',I3,' **** SFR ****')
C      !WRITE 100 204,((SFR(I1,I2),I1=1,6),I2=1,KST2)
C      !WRITE 100 206,KSTOLD,KSS
C 206  FORMAT(5X,' KSTOLD=',I3,' KSS=',I3,' **** SPTU (AFTER) ****')
C      !WRITE 100 204,((SPTU(I1,I2),I1=4,9),I2=KSTOLD,KSS)
      RETURN
      END
      FUNCTION WEIGHT(MA,K,DAZ,TQ)
C O[PE[E[EH[E [P[BE[EHHO[ BEPO[THOCT[ PAC[A[A MA HA K [PA[MEHTOB. HOMEP
C [PA[M. [AH[ B MACC[BE J(I). DAZ-[HEP[[[ MA,BK[[[A[ [HEP[[[ [OKO[. TQ-
C K[[OHOBCK[[ [AP[EP. AKV-[P[BE[EHH[[ O[[EM [[[ O[HO[O [PA[MEHTA. AKM-[P
C BE[EH[E C[[HOB[X [ MACCOB[X [AKTOPOB. RPM-[AKTOP[A[ TO[[ECTBEHH[X [[E-
C HOB. GAM-O[PATHA[ [AMMA-[[HK[[[. TN-K[HET[[ECKA[ [HEP[[[ OCKO[KOB.
      COMMON /SPDM/ MS(12,11),DM(12,11) /VAK/VAK /GAF/GAF(20)
      COMMON /BLJ/J(20),IEND /BLIAN/IAN(100) /BLIZN/IZN(100)
      COMMON /BLISP/ISP(100) /BLENN/ENN(100) /FKAP/FKAP
      BNQ=DAZ
      RMQ=0.
      SPM=1.
      AKM=1.
      AKV=MA*VAK
c --- correction on multiple-dependent freez-out (translat.) volume ---
      AMA03=(FLOAT(MA))**0.333333
      AKV=(AKV/FKAP)*((1.+1.4*(K**0.333333-1)/(1.17*AMA03))**3-1.)
c ---
      DO 40 I=1,K
      JI=J(I)
      IA=IAN(JI)
      IZ=IZN(JI)
      AKM=AKM*IA
      SPM=SPM*ISP(JI)
      BNQ=BNQ-GMASS(IA,IZ)-0.001*ENN(JI)
 40   CONTINUE
      TN=BNQ
      IF(TN.LE.TQ) GO TO 20
      TN=TN-TQ
      AKM=AKM/MA
      AKM=AKM*SQRT(AKM)*SPM
      IF(K.GT.2) GO TO 29
      RMQ=1.1283792*AKV*AKM*SQRT(TN)
      IF(J(1).EQ.J(2)) RMQ=0.5*RMQ
      GO TO 20
 29   CONTINUE
      KI=K-1
      RPM=1.
      VMK=1.
      TEQ=2.71828183*TN/(1.5*K-2.5)
      VTK=AKV*TEQ*SQRT(TEQ)
      DO 21 I=1,KI
      MRS=1
      IK=I+1
      VMK=VMK*VTK
      DO 22 I1=IK,K
      IF(J(I).EQ.J(I1)) MRS=MRS+1
 22   CONTINUE
      RPM=RPM*MRS
 21   CONTINUE
      GAM=GAF(K)
      RMQ=VMK*GAM*AKM/(TEQ*RPM)
 20   WEIGHT=RMQ
      RETURN
      END
      FUNCTION GMASS(MA,MZ)
C O[PE[E[EH[E MACC[ [E[K[X [[EP B [[B.([C[O[[[[[TC[ [[[EPO[H[E E[[H[[[)
C ( OCHOBHOE COCTO[H[E)
      COMMON /SPDM/ MS(12,11),DM(12,11)
      I2=MZ
      I=MA-I2
      IF(I.LT.0.OR.I.GT.11) GO TO 10
      IF(I2.LT.0.OR.I2.GT.10) GO TO 10
      DMP=DM(I+1,I2+1)
      IF(DMP.GT.90) GO TO 10
      GMASS=0.931437*(I+I2)+0.001*DMP
      RETURN
c 10   GMASS=0.939*MA
c ---- botvina change ----
 10   CONTINUE
      AA=MA
      ZZ=MZ
      CALL DELAM(AA,ZZ,DLMPP,DSH,BAR)
      GMASS=0.931437*AA+0.001*DLMPP
c ----
      RETURN
      END
      FUNCTION BCUL(K)
C  O[PE[E[EH[E K[[OHOBCKO[O [AP[EPA [[[ PA[BA[A HA KOH[[[[PA[[[ J(1-K)
C  [[ K [ACT[[.
      COMMON /FKAP/FKAP /BLJ/J(20),IEND /BLIAN/IAN(100) /BLIZN/IZN(100)
      DIMENSION IA(20),IZ(20)
      A=0.
      Z=0.
      COEF=(3.*1.44/(5.*1.3))*(1./(1.+FKAP)**0.333333)
      DO 1 I=1,K
      JI=J(I)
      IA(I)=IAN(JI)
      A=A+IA(I)
      IZ(I)=IZN(JI)
      Z=Z+IZ(I)
 1    CONTINUE
      EC=0.
      DO 4 I=1,K
      AP=IA(I)
      ZP=IZ(I)
      EC=EC+ZP*ZP/(AP**0.333333)
 4    CONTINUE
      EC=Z*Z/(A**0.333333)-EC
      BCUL=0.001*COEF*EC
      IF(BCUL.LT.0.) BCUL=0.
      RETURN
      END
      SUBROUTINE HELPRA
C MS(I,J)-MACC[B C[[HOB[X [AKTOPOB [E[K[X [[EP =(2*S+1),[[E I=A-Z+1,
C J=Z+1. DM(I,J)-MACC[B [E[EKTOB MACC B [[[EPO[H[X E[[H[[AX. [[[T[BA[TC[
C [[PA C 0.LE.Z.LE.10 [ Z.LE.A.LE.(11+Z); MO[HO C[[TAT[ PA[BA[ [E[K[X
C [[EP C Z0.LE.10 [ (A0-Z0).LE.11.
C ([[[T[BA[TC[ C[[H[             OCHOBH[X      COCTO[H[[,   [EPETC[ [HEP
C [[[ OCHOBHO[O COCTO[H[[ - [[[ [[EP HA[[HA[ C [[T[[).
C KPOME TO[O,[POBO[[TC[ BC[OMO[ATE[[H[E B[[[C[EH[[ [[[ O[PE[E[EH[[ [AMMA
C [[HK[[[ ([O CT[P[[H[[) [ [P[BE[EHHO[O O[[EMA.
      COMMON /VAK/VAK /GAF/GAF(20) /SPDM/MS(12,11),DM(12,11)
      COMMON /FKAP/FKAP /ILEVRA/ILEVRA
      COMMON /BLIAN1/IAN1(100) /BLIZN1/IZN1(100) /BLENN1/ENN1(100)
      COMMON /BLIAN2/IAN2(100) /BLIZN2/IZN2(100) /BLENN2/ENN2(100)
      COMMON /BLIAN/IAN(100) /BLIZN/IZN(100) /BLENN/ENN(100)
      COMMON /BLISP1/ISP1(100) /BLISP2/ISP2(100) /BLISP/ISP(100)
      DIMENSION IPA(12),JPZ(12),DMP(12),MSP(12),AMP(12),MD(12,11),
     *MS1(12,11)
      DATA MS1/
     * 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     * 2, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     * 0, 2, 1, 4, 1, 0, 1, 0, 0, 0, 0, 0,
     * 0, 0, 4, 3, 4, 5, 4, 0, 0, 0, 0, 0,
     * 0, 0, 1, 4, 1, 4, 1, 2, 1, 0, 0, 0,
     * 0, 0, 0, 5, 4, 7, 4, 3, 4, 0, 0, 0,
     * 0, 0, 0, 4, 1, 4, 1, 2, 1, 2, 1, 0,
     * 0, 0, 0, 0, 0, 3, 2, 3, 2, 5, 2, 0,
     * 0, 0, 0, 0, 0, 4, 1, 2, 1, 6, 1, 6,
     * 0, 0, 0, 0, 0, 0, 0, 1, 6, 3, 2, 5,
     * 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 1, 4/
      DATA MD/
     *99000,
     * 8071,99999,99999,99999,99999,99999,99999,99999,99999,99999,99999,
     * 7289,
     *13136,14950,25920,33790,99999,99999,99999,99999,99999,99999,99999,
     *99000,
     *14931, 2425,11390,17597,26110,31603,42030,50130,65000,75240,89260,
     *99000,
     *25130,11680,14087,14907,20947,24955,33250,40940,52940,61570,72280,
     *99000,
     *31700,18375,15769, 4942,11348,12607,20179,25072,35720,40720,51210,
     *99000,
     *99999,27940,22920,12416,12052, 8668,13370,16560,24230,29410,37960,
     *99000,
     *99999,35093,28912,15702,10650,00000, 3125, 3020, 9873,13693,17560,
     *99000,99999,99999,99999,25450,17338, 5346, 2864,  102, 5683, 7871,
     *13274,99000,99999,99999,42700,32070,23110, 8010, 2860,-4737, -810,
     * -783, 3332,99000,99999,99999,99999,39700,33380,17610,10693, 1952,
     *  873,-1486,  -16,99000,99999,99999,99999,49400,36400,24110,16470,
     * 5320, 1750,-7040,-5730/
      DO 22 I=1,12
      DO 22 I2=1,11
      MS(I,I2)=MS1(I,I2)
 22   DM(I,I2)=0.001*MD(I,I2)
      VL=1.3/(0.21*SQRT(0.94))
      VOL=FKAP*VL*VL*VL
      VAK=VOL*SQRT(2./3.14159265)*0.333333
      GAF(1)=0.
      GAF(2)=1./SQRT(3.1416)
      DO 10 K=3,20
      QK=1./(1.5*K-2.5)
      GQ=1.+QK*(1./12.+QK*(1./288.-QK*(139./51840.)))
 10   GAF(K)=SQRT(QK*0.1591549)/GQ
      IF(ILEVRA-1) 31,30,31
 30   DO 32 I=1,100
      IAN(I)=IAN1(I)
      IZN(I)=IZN1(I)
      ISP(I)=ISP1(I)
 32   ENN(I)=ENN1(I)
      RETURN
 31   DO 33 I=1,100
      IAN(I)=IAN2(I)
      IZN(I)=IZN2(I)
      ISP(I)=ISP2(I)
 33   ENN(I)=ENN2(I)
      RETURN
      END
      SUBROUTINE DISNUM(IA0,IZ0,K)
C  HAXO[[EH[E K [ACT[[ B KOH[[[[PA[[[ B B[[E [X HOMEPOB J(1-K) B MAC-
C  C[BAX IAN (B [OP[[.BO[POC.A) [ IZN (B [OP[[.BO[POC.Z [P[ [AHHOM A)
C  BO[MO[HA[ KOH[[[[PA[[[ IEND=0, KOHE[ IEND=1.
C  HA[.KOH[[[.J(2-K)=1, J(1)=0 [A[AETC[ [EPE[ DISNUM.
      COMMON /BLIAN/IAN(100) /BLIZN/IZN(100) /BLJ/J(20),IEND
 5    L=0
 1    L=L+1
      IF(L.LT.K) GO TO 10
      IF(L.GE.K) GO TO 11
 10   J(L)=J(L)+1
      IF(J(L).GT.J(L+1)) GO TO 2
      ISA=0
      DO 12 I=1,K
      JI=J(I)
 12   ISA=ISA+IAN(JI)
      IF(ISA.GT.IA0) GO TO 2
      IF(ISA.LT.IA0) GO TO 5
 4    ISZ=0
      DO 13 I=1,K
      JI=J(I)
 13   ISZ=ISZ+IZN(JI)
      IF(ISZ.NE.IZ0) GO TO 5
      IEND=0
      RETURN
 2    J(L)=1
      GO TO 1
 11   J(L)=J(L)+1
      IF(J(L).GT.100) GO TO 15
      ISA=0
      DO 14 I=1,K
      JI=J(I)
 14   ISA=ISA+IAN(JI)
      IF(ISA.LT.IA0) GO TO 5
      IF(ISA.EQ.IA0) GO TO 4
 15   IEND=1
      RETURN
      END
      BLOCK DATA B3
C  TA[[[[[,[A[A[[[E MACC[ (IAN),[AP[[[ (IZN),C[[H.MHO[.(ISP),[HEP[.OCH.
C  COCT.(ENN) [ACT[[ .
      COMMON /BLIAN1/IAN1(100) /BLIZN1/IZN1(100) /BLENN1/ENN1(100)
      COMMON /BLIAN2/IAN2(100) /BLIZN2/IZN2(100) /BLENN2/ENN2(100)
      COMMON /BLISP1/ISP1(100) /BLISP2/ISP2(100)
C ******   OCHOBH.+ [CTO[[[B[E BO[[[[[. COCTO[H[[  ******
      DATA IAN2/
     * 1, 1, 2, 3, 3, 4, 5, 5, 6, 6,   6, 7, 7, 7, 7, 8, 8, 8, 9, 9,
     *10,10,10,10,10,10,10,10,10,10,  10,10,
     *                                11,11,11,11,11,11,11,11,11,11,
     *11,11,11,11,11,11,11,11,11,12,  12,12,12,12,12,13,13,13,13,13,
     *14,14,14,14,14,14,14,14,14,14,  14,14,14,14,14,15,15,15,15,15,
     *15,15,15,15,15,15,15,15,15,15,  16,16,16,16,16,16,16,16/
      DATA IZN2/
     * 0, 1, 1, 1, 2, 2, 2, 3, 2, 3,   3, 3, 3, 4, 4, 3, 3, 4, 4, 5,
     * 4, 4, 4, 4, 4, 5, 5, 5, 5, 5,   6, 6,
     *                                 5, 5, 5, 5, 5, 5, 5, 5, 6, 6,
     * 6, 6, 6, 6, 6, 6, 6, 6, 6, 5,   5, 5, 5, 6, 6, 6, 6, 6, 6, 7,
     * 6, 6, 6, 6, 6, 7, 7, 7, 7, 7,   7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
     * 7, 7, 7, 7, 7, 8, 8, 8, 8, 8,   7, 7, 7, 7, 8, 8, 8, 8/
      DATA ISP2/
     * 2, 2, 3, 2, 2, 1, 4, 4, 1, 3,   1, 4, 2, 4, 2, 5, 3, 1, 4, 4,
     * 1, 5, 8, 1, 5, 7, 3, 1, 3, 5,   3, 5,
     *                                 4, 2, 6, 4,10, 6, 4, 6, 4, 2,
     * 6, 4, 2, 8, 6, 4, 4, 6, 8, 3,   5, 5, 4, 1, 5, 2, 2, 4, 6, 2,
     * 1, 3, 8, 6, 5, 3, 1, 3, 1, 5,   3, 7, 3, 7, 5, 2, 8, 4,10, 8,
     * 2, 4,14,14, 8, 2, 8, 4,10, 8,   5, 1, 7, 3, 1, 8, 5, 3/
      DATA ENN2/
     * 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,16.76,16.66, 0.00, 0.00,
     * 3.56, 0.00, 0.48, 0.00, 0.43, 0.00, 0.98, 0.00, 0.00, 0.00,
     * 0.00, 3.37, 5.96, 6.18, 6.26, 0.00, 0.72, 1.74, 2.15, 3.59,
     * 0.00, 3.35,
     * 0.00, 2.13, 4.44, 5.02, 6.76, 7.29, 7.98, 8.56, 0.00, 2.00,
     * 4.32, 4.80, 6.34, 6.48, 6.90, 7.50, 8.10, 8.42, 8.66, 0.00,
     * 0.95, 1.67, 2.65, 0.00, 4.44, 0.00, 3.09, 3.68, 3.85, 0.00,
     * 0.00, 6.09, 6.69, 6.96, 7.34, 0.00, 2.31, 3.95, 4.92, 5.11,
     * 5.69, 5.83, 6.20, 6.44, 7.03, 0.00, 5.28, 6.32, 7.22, 7.57,
     * 8.31, 8.57, 9.15, 9.79,10.00, 0.00, 5.22, 6.18, 6.83, 7.28,
     * 0.00, 0.12, 0.30, 0.40, 0.00, 6.10, 6.92, 7.12/
C  ******** OCHOBH. COCTO[H[[ *******
      DATA IAN1/
     * 1, 1, 2, 3, 3, 4, 5, 5, 6, 6,   7, 7, 8, 8, 9, 9,10,10,10,11,11,
     *12,12,13,13,14,14,15,15,16,16,  17,18,67*0/
      DATA IZN1/
     * 0, 1, 1, 1, 2, 2, 2, 3, 2, 3,   3, 4, 3, 4, 4, 5, 4, 5, 6, 5, 6,
     * 5, 6, 6, 7, 6, 7, 7, 8, 7, 8,   8, 8,67*0/
      DATA ISP1/
     * 2, 2, 3, 2, 2, 1, 4, 4, 1, 3,   4, 4, 5, 1, 4, 4, 1, 7, 3, 4, 4,
     * 3, 1, 2, 2, 1, 3, 2, 2, 5, 1,   6, 1,67*0/
      DATA ENN1/6*0.,16.76,16.66,92*0./
      END
      SUBROUTINE PINT(P,CT,ST,CF,SF,T,CM)
C O[PE[E[EH[E [O MACCE CM [ [M[[[[C[ P(3) [ACT[[[,EE K[H. [HEP[[[ T
C [ [[[OB - CT,ST,CF,SF ,O[PE[E[[[[[X HA[PAB[EH[E [O[ETA.
      DIMENSION P(3)
      PZ=P(3)**2
      PQ=P(1)**2+P(2)**2+PZ
      IF(PQ.LE.0.) GO TO 11
      CTQ=PZ/PQ
      T=SQRT(PQ+CM**2)-CM
      IF(CTQ.GE.1.) GO TO 10
      CT=SQRT(CTQ)
      IF(P(3).LE.0.) CT=-CT
      ST=SQRT(1.-CTQ)
      PMH=ST*SQRT(PQ)
      CF=P(1)/PMH
      SF=P(2)/PMH
      RETURN
 10   ST=0.
      CT=1.
      SF=0.
      CF=1.
      RETURN
 11   ST=0.
      CT=1.
      SF=0.
      CF=1.
      T=0.
      RETURN
      ENTRY TINP(P,CT,ST,CF,SF,T,CM)
C O[PE[E[EH[E [O MACCE CM,K[H.[HEP[[[ T [ [[[AM [ACT[[[,EE [M[[[[CA P(3)
      PM=SQRT(T*(T+2.*CM))
      P(3)=PM*CT
      PM=PM*ST
      P(1)=PM*CF
      P(2)=PM*SF
      RETURN
      END
      SUBROUTINE CLPV(P1,V,P2,E1)
C O[PE[E[EH[E [M[[[[CA [ACT[[[ P2(3) B [.C.K ,C[[TA[,[TO CTAPA[ C.K.([DE
C [M[[[[C [ACT[[[ P1(3)) [B[[ETC[ CO CKOP.V(3).E1-[O[H.[HEP[.[ACT[[[ B C
      DIMENSION P1(3),V(3),P2(3)
      SPV=0.
      V2=0.
      DO 10 I=1,3
      V2=V2+V(I)**2
      SPV=SPV+P1(I)*V(I)
 10   CONTINUE
      GAM=1./SQRT(1.-V2)
      TEMP=GAM*(SPV*GAM/(GAM+1.)+E1)
      DO 20 K=1,3
 20   P2(K)=P1(K)+V(K)*TEMP
      RETURN
      END
      SUBROUTINE ISOTR(ANL)
      DIMENSION ANL(3)
      CT=1.-2.*RNDM(-1)
      ST=SQRT(ABS(1.-CT*CT))
      ANL(3)=CT
      FI=6.2831843*RNDM(-1)
      ANL(1)=COS(FI)*ST
      ANL(2)=SIN(FI)*ST
      RETURN
      END
      SUBROUTINE DISIMP(K,AMK,PN0,TN)
C O[PE[E[EH[E [M[[[[COB K [PA[MEHTOB- PN0(3,1-K) ([[B/C) B C[CTEME [X
C [.[.  AMK(1-K) ([[B) -[X MACC[. TN ([[B) - [X [O[HA[ K[H. [HEP[[[.
C (AHA[O[[[HO [.[.  KO[[[OB[-1970, C.229-239)
      DIMENSION ANL(3),PNC(3),VRS(3),PN0(3,500),AMK(500)
      SMK=0.
      DO 18 I=1,K
 18   SMK=SMK+AMK(I)
      DO 11 I=1,3
 11   VRS(I)=0.
      TKN=TN
      KL=K-1
      DO 10 L=1,KL
      LK=K-L+1
      AMP=AMK(LK)
      SMK=SMK-AMP
      AMR=SMK
      TPR=TKN
      IF(LK.LT.3) GO TO 20
      TKM=TKN*RNKSI(LK-1)
      TPR=TKN-TKM
 20   PMC=SQRT(2.*((AMP*AMR)/(AMP+AMR))*TPR)
      CALL ISOTR(ANL)
      DO 12 I=1,3
      PNC(I)=PMC*ANL(I)
      PN0(I,LK)=VRS(I)*AMP+PNC(I)
      VRS(I)=VRS(I)-PNC(I)/AMR
 12   CONTINUE
      TKN=TKM
 10   CONTINUE
      DO 21 I=1,3
 21   PN0(I,1)=VRS(I)*AMR
      RETURN
      END
      FUNCTION RNKSI(K)
C HAXO[[EH[E RNKSI(K) O[PE[E[[[[EE K[H.[HEP[[[ K [ACT[[ -T(K) B C[CTEME
C [X [.[. [EPE[ K[H.[HEP[[[ K+1 [ACT[[ -T(K+1) B C[CTEMME [X [.[.
C           (   T(K)=RNKSI(K)*T(K+1)    )
      CSIM=(3.*K-5.)/(3.*K-4.)
      PEX=1.5*K-2.5
      FCSIM=SQRT(1.-CSIM)*CSIM**PEX
 20   CSI=RNDM(-1)
      FCSI =SQRT(1.-CSI)*CSI**PEX
      RF=FCSIM*RNDM(-1)
      IF(RF.GT.FCSI) GO TO 20
      RNKSI=CSI
      RETURN
      END
      SUBROUTINE POVRAZ(KSTOLD,KST1,KST2)
C  [OBTOPH[[ PA[BA[ HECTA[[[[H[X [PA[MEHTOB. BXO[:KST1 [ACT[[ B MACC[BE
C  SPTU(10,500), B[XO[: KST2 [ACT[[ B MACC[BE SFR(6,99)
      COMMON /BLOKC/SPTU(10,500) /BLOK79/SFR(6,99)
      DIMENSION P9(3,500),A9(500),PC9(3),PL9(3),V9(3),Z9(500)
      DO 1 I=1,6
      DO 1 J=1,99
 1    SFR(I,J)=0.
      KST2=0
      DO 2 I=KSTOLD,KST1
      CMN=SPTU(9,I)
      ZP=SPTU(8,I)
      TMN=SPTU(7,I)
      CF=SPTU(6,I)
      SF=SPTU(5,I)
      CT=SPTU(4,I)
      IF(CT.LT.1.) ST=SQRT(1.-CT*CT)
      IF(CT.GE.1.) ST=0.
      IA=INT(CMN/0.94+0.5)
      IZ=ZP
      GAM=1.+TMN/CMN
      V=SQRT((GAM**2-1.)/GAM**2)
      V9(3)=V*CT
      V9(2)=V*ST*SF
      V9(1)=V*ST*CF
      IF(IA.EQ.5.AND.IZ.EQ.2) GO TO 21
      IF(IA.EQ.5.AND.IZ.EQ.3) GO TO 22
      IF(IA.EQ.8.AND.IZ.EQ.4) GO TO 23
      IF(IA.EQ.9.AND.IZ.EQ.5) GO TO 24
      KST2=KST2+1
      DO 3 J=1,6
      J1=J+3
 3    SFR(J,KST2)=SPTU(J1,I)
      GO TO 2
 21   K9=2
      A9(1)=3.728173
      Z9(1)=2.
      A9(2)=0.939508
      Z9(2)=0.
      T9=0.001*(16.76+11.39-8.071-2.425)
      GO TO 100
 22   K9=2
      A9(1)=3.728173
      Z9(1)=2.
      A9(2)=0.938726
      Z9(2)=1.
      T9=0.001*(16.66+11.68-7.289-2.425)
      GO TO 100
 23   K9=2
      A9(1)=3.728173
      Z9(1)=2.
      A9(2)=3.728173
      Z9(2)=2.
      T9=0.001*(4.942-2*2.425)
      GO TO 100
 24   K9=3
      A9(1)=3.728173
      Z9(1)=2.
      A9(2)=3.728173
      Z9(2)=2.
      A9(3)=0.938726
      Z9(3)=1.
      T9=0.001*(12.416-7.289-2*2.425)
 100  CONTINUE
      CALL DISIMP(K9,A9,P9,T9)
      DO 101 J=1,K9
      DO 102 II=1,3
 102  PC9(II)=P9(II,J)
      E9=SQRT(PC9(1)**2+PC9(2)**2+PC9(3)**2+A9(J)**2)
      CALL CLPV(PC9,V9,PL9,E9)
      CALL PINT(PL9,C1,S1,C2,S2,T,A9(J))
      KST2=KST2+1
      SFR(1,KST2)=C1
      SFR(2,KST2)=S2
      SFR(3,KST2)=C2
      SFR(4,KST2)=T
      SFR(5,KST2)=Z9(J)
      SFR(6,KST2)=A9(J)
      IF(A9(J).GT.0.5) SFR(6,KST2)=0.94*INT(A9(J)/0.94+0.5)
 101  CONTINUE
 2    CONTINUE
C     !WRITE 100 11,KST1,KST2
C11   FORMAT(25X,'[OC[E POVRAZ               KST1=',I3,' KST2=',I3)
C     DO 12 I=1,KST2
C     !WRITE 100 10,SFR(6,I),SFR(5,I),SFR(4,I),SFR(1,I),SFR(3,I),SFR(2,I)
C10   FORMAT(2X,' A=',F7.3,'  Z=',F7.3,'  T=',F9.5,'  C1=',F9.5,
C    *'       C2,S2=',2F9.5)
C12   CONTINUE
      RETURN
      END
       SUBROUTINE EVANUC(ENEXT,ATWGHT,CHARGE,PNX,PNY,PNZ,KSTART)
       COMMON /BL0999/RNCL/BL1000/AM,AMF,AMCB(35)/AMNUCL/AMNEV(3)
       COMMON /BL1001/T1Y(130)/BL1002/T2XY(200)/BL1005/AJ(35)
     */BL1006/ZJ(35)/BL1014/GAM(35)/BL1011/VJ(35)/BL1015/RJ(35)
      COMMON /BL1003/U,A,Z/BL1009/AFJ(35)/BL1010/ZFJ(35)
      COMMON/BLOKC/SPTU(10,500)/BLANGL/ANGL(4)/BL1100/AMS,AMFS
      COMMON/BLFIS/FIS,EF,SNT(6,6)/BLCAN/ICAN
      COMMON /EXILEV/ ENLEV(20,35),SPLEV(20,35)
      COMMON /BLSECON/LSECON
       DIMENSION GJ(35),BJ(35),GG(35),VN(3),PN(3),PNL(3),PP(3),PPL(3)
      DIMENSION DSH2(35)
C [C[APEH[E [E[K[X [PA[MEHTOB [[ T[[.[[PA [O BA[CKO[[[. [[[T[BAETC[
C COXPAHEH[E [HEP[[[,[M[[[[CA [ MOMEHTA [M[[[[CA(K[ACC[[.). [E[EH[E
C [[[T[BAETC[ KAK BO[MO[H[[ KAHA[ [PEKPA[A[[[[ [C[AP[TE[[H[[ KACKA[.
C PA[MEPHOCT[: [HEP[[[-[[B,[M[[[[C-[[B/C,CKOPOCT[-1/C,MACCA-AT.HOMEP.
C       ******************************
C  ICAN=0 -[[[T[BAETC[ 6 KAHA[OB [ [E[EH[E; ICAN=1 - 32 KAHA[A [ [E[EH[E
C  ICAN=2 - 6 KAHA[OB, [E[EH[E, 34 [ 35 KAHA[; ICAN=3 ([[[ [P.) -BCE KAH
C       ******************************
      ILEVEL=3
      CALL DATEVP(ILEVEL)
      FIS=0.
      SNT(1,1)=ENEXT
      SNT(2,1)=ATWGHT
      SNT(3,1)=CHARGE
      SNT(4,1)=PNX
      SNT(5,1)=PNY
      SNT(6,1)=PNZ
      U=ENEXT*1000.
      A=ATWGHT
      Z=CHARGE
      AMNEV(1)=0.
      AMNEV(2)=0.
      AMNEV(3)=0.
      REMN=940.*A
      PNL(1)=PNX*1000.
      PNL(2)=PNY*1000.
      PNL(3)=PNZ*1000.
      ENL=SQRT(PNL(1)**2+PNL(2)**2+PNL(3)**2+REMN**2)
      VN(1)=PNL(1)/ENL
      VN(2)=PNL(2)/ENL
      VN(3)=PNL(3)/ENL
      DO 20 K=KSTART,500
      IN=35
      DO 4 I=1,IN
      VJ(I)=0.
      AMCB(I)=0.
C      DSH2(I)=0.
 4    RJ(I)=-1000.
      CALL DELAM(A,Z,DL1,DSHEL1,BARIER)
      AM=AMS
C      IF(U.GT.0.) AM=AMS*(1.+(1.-EXP(-0.06*U))*DSHEL1/U)
      IF(Z.GE.89.) AMF=1.04*AMS
      IF(Z.GT.85.AND.Z.LT.89) AMF=AMS*(1.04+0.01*(89.-Z))
      IF(Z.LE.85.) AMF=1.08*AMS
      DO 6 I=1,IN
      IF(I.EQ.33) GO TO 6
      AFJ(I)=A-AJ(I)
      ZFJ(I)=Z-ZJ(I)
        IF(ICAN.EQ.0.AND.I.GT.6) GO TO 6
        IF(ICAN.EQ.1.AND.I.GT.32) GO TO 6
        IF(ICAN.EQ.2.AND.I.GT.6.AND.I.LE.32) GO TO 6
        IF(AFJ(I).LT.ZFJ(I)) GO TO 6
        IF(AFJ(I).LT.AJ(I).OR.ZFJ(I).LT.ZJ(I)) GO TO 6
      JAF=INT(AFJ(I)+0.5)
      JZF=INT(ZFJ(I)+0.5)
      JNF=JAF-JZF
      ODD=11.*(2+2*(JNF/2)-JNF+2*(JZF/2)-JZF)/SQRT(AFJ(I))
      IF(RNCL.GE.0.1) RADNCL=RNCL
      IF(RNCL.LT.0.1)
     *RADNCL=2.173*(1+0.006103*ZJ(I)*ZFJ(I))/(1+0.009443*ZJ(I)*ZFJ(I))
      VJ(I)=COLOMB(I,RADNCL)
      CALL DELAM(AFJ(I),ZFJ(I),DL2,DSHEL2,BAR2)
C      DSH2(I)=DSHEL2
      CALL DELAM(AJ(I),ZJ(I),DL3,DSHEL3,BAR3)
      BJ(I)=DL2+DL3-DL1
c ******** in accordance with DELA8
      IF(A.LE.55.) RJ(I)=U-(BJ(I)+VJ(I))-ODD
      TCOEF=1.-(A-55.)/10.
      IF(A.GT.55.AND.A.LT.65.) RJ(I)=U-(BJ(I)+VJ(I))-ODD*TCOEF
      IF(A.GE.65.) RJ(I)=U-(BJ(I)+VJ(I))
c ********
C     RJ(I)=U-(BJ(I)+VJ(I))-ODD
C     RJ(I)=U-(BJ(I)+VJ(I))
C      IF(RJ(I).GT.0.)
C     *AMCB(I)=AMS*(1.+(1.-EXP(-0.06*RJ(I)))*DSHEL2/RJ(I))
 6    CONTINUE
      EF=BARIER
C      AMCB(33)=AMF
      IF(A.GT.100.) RJ(33)=U-BARIER
       DO 7 I=1,IN
      IF(RNCL.GE.0.1) RADNCL=RNCL
      IF(RNCL.LT.0.1)
     *RADNCL=2.173*(1+0.006103*ZJ(I)*ZFJ(I))/(1+0.009443*ZJ(I)*ZFJ(I))
      AMC=AM
C      AMC=AMCB(I)
      GJ(I)=GAMMAC(I,A,U,AM,AMC,AMF,RADNCL)
      RJI=RJ(I)
      GAMI=GAM(I)
      DO 71 IL=1,20
      IF(SPLEV(IL,I).LT.0.1) GO TO 70
      GAM(I)=SPLEV(IL,I)*AJ(I)
      RJ(I)=RJI-ENLEV(IL,I)
C      IF(RJ(I).GT.0.)
C     *AMC=AMS*(1.+(1.-EXP(-0.06*RJ(I)))*DSH2(I)/RJ(I))
      GJ(I)=GJ(I)+GAMMAC(I,A,U,AM,AMC,AMF,RADNCL)
 71   CONTINUE
 70   RJ(I)=RJI
      GAM(I)=GAMI
7     CONTINUE
      G=0.
      DO 10 I=1,IN
10    G=G+GJ(I)
        IF(G) 11,11,12
 11   PNX=PNL(1)*0.001
      PNY=PNL(2)*0.001
      PNZ=PNL(3)*0.001
      ENEXT=U*0.001
      ATWGHT=A
      CHARGE=Z
      REMN=940.*A
      SPTU(9,K)=REMN*0.001
      SPTU(8,K)=Z
      CALL PINT(PNL,CT,ST,CF,SF,T,REMN)
      SPTU(4,K)=CT
      SPTU(5,K)=SF
      SPTU(6,K)=CF
      SPTU(7,K)=T*0.001
      SPTU(1,K)=1.
      SPTU(3,K)=FIS
      IF(K.GT.KSTART) SPTU(2,K)=2.
      IF(LSECON.GT.0) SPTU(2,K)=2.
      KSTART=K+1
      RETURN
 12   CONTINUE
      DO 13 J=2,IN
13    GJ(J)=GJ(J-1)+GJ(J)
      BB=RNDM(-1)
      B=BB*G
      DO 14 J=1,IN
        IF(B-GJ(J)) 15,14,14
15    LM=J
      GO TO 16
 14   CONTINUE
 16     IF(LM-33) 18,17,18
 17   FIS=1.
      KSTART=K
      CALL FISION(U,A,Z,PNL(1),PNL(2),PNL(3),KSTART)
      RETURN
 18   AMC=AM
C18    AMC=AMCB(LM)
      EP1=TKIN(LM,AMC,AMF)
      EP2=ZJ(LM)
      EP3=940.*AJ(LM)
C     ALMAX=0.219327*1.2*(AFJ(LM)**0.333333+AJ(LM)**0.333333)
C    +*SQRT((AJ(LM)*AFJ(LM)/(AJ(LM)+AFJ(LM)))*(EP1-VJ(LM)))
C     BR=RNDM(-1)
C     AL=ALMAX*SQRT(BR)
C     BR1=1.-2.*RNDM(-1)
C     BR2=6.283185*RNDM(-1)
C     AMNEV(1)=AMNEV(1)-AL*SQRT(1.-BR1**2)*COS(BR2)
C     AMNEV(2)=AMNEV(2)-AL*SQRT(1.-BR1**2)*SIN(BR2)
C     AMNEV(3)=AMNEV(3)-AL*BR1
      U=U-BJ(LM)-EP1
      A=AFJ(LM)
      Z=ZFJ(LM)
C  VPM-OTHOC[TE[[HA[ CKOPOCT[ [PA[MEHTA
      VPM=SQRT((2.*EP1)/(EP3*AFJ(LM)/(AFJ(LM)+AJ(LM))))
      CALL ISANGL
C  [M[[[[C [PA[MEHTA B C.[.[.
      PP(1)=VPM*ANGL(4)*ANGL(3)*EP3/(1.+AJ(LM)/A)
      PP(2)=VPM*ANGL(4)*ANGL(2)*EP3/(1.+AJ(LM)/A)
      PP(3)=VPM*ANGL(1)*EP3/(1.+AJ(LM)/A)
      PN(1)=-PP(1)
      PN(2)=-PP(2)
      PN(3)=-PP(3)
      EP=SQRT(PP(1)**2+PP(2)**2+PP(3)**2+EP3**2)
      EN=SQRT(PN(1)**2+PN(2)**2+PN(3)**2+(940.*A)**2)
      CALL CLPV(PP,VN,PPL,EP)
      CALL CLPV(PN,VN,PNL,EN)
      ENL=SQRT(PNL(1)**2+PNL(2)**2+PNL(3)**2+(940.*A)**2)
      VN(1)=PNL(1)/ENL
      VN(2)=PNL(2)/ENL
      VN(3)=PNL(3)/ENL
      CALL PINT(PPL,CT,ST,CF,SF,T,EP3)
      SPTU(4,K)=CT
      SPTU(5,K)=SF
      SPTU(6,K)=CF
      SPTU(7,K)=T*0.001
      SPTU(8,K)=EP2
      SPTU(9,K)=EP3*0.001
      SPTU(2,K)=2.
20    CONTINUE
      !WRITE 100 21,U,A,Z
21    FORMAT(35X,37HMASSIV SPT EXCEEDED AFTER EVAPORATION/40X,2HU=,
     *F10.5,4H  A=,F5.1,4H  Z=,F4.1)
       RETURN
       END
      SUBROUTINE DATEVP(I)
C I=1 (I<2) -[[ET OCHOBH[X COCTO[H[[ [C[AP[[[E[C[ [ACT[[[,I=2-[[ET
C [EPB[X BO[[[[[EHH[X COCTO[H[[ ([O [POBHE[,KOTOP[E [C[[CKA[T H[K[OH[),
C I=3 (I>2) -[CPE[HEHH[[ [[ET H[[H[X [POBHE[ [P[ B[COK[X [HEP[[[X BO[[.
      COMMON /BLGAM/ GAM1(35),GAM2(35) /BL1014/ GAM(35)
      IF(I-2) 1,1,2
 1    DO 3 J=1,35
 3    GAM(J)=GAM1(J)
      IF(I.EQ.2) CALL FILLEV
      GO TO 5
 2    DO 4 J=1,35
 4    GAM(J)=GAM2(J)
 5    RETURN
      END
       FUNCTION COLOMB(L,RADNCL)
       COMMON /BL1006/ZJ(35)/BL1010/ZFJ(35)/BL1009/AFJ(35)
     */BL1005/AJ(35)/BL1003/U,A,Z
C B[[[C[EH[E K[[OHOBCKO[O [AP[EPA [C[[CKAH[[ [PA[MEHTA COPTA L (M[B)
       IF(L-1)1,1,2
1      COLOMB=0.
       RETURN
2      TEMP1=1.44/RADNCL
       COLOMB=TEMP1*((ZJ(L)*ZFJ(L))/(AJ(L)**.33333+AFJ(L)**.33333))
       IF(COLOMB) 3,3,4
3       COLOMB=0.
4       CONTINUE
       RETURN
       END
      SUBROUTINE ISANGL
C [[OTPO[H[[ PO[[[P[[ [[[OB B[[ETA [PA[MEHTA B C[CTEME OCTAT.[[PA.
      COMMON /BLANGL/ ANGL(4)
      ANGL(1)=1.-2.*RNDM(-1)
      ANGL(4)=SQRT(1.-ANGL(1)**2)
      F=2.*3.1415926536*RNDM(-1)
      ANGL(2)=SIN(F)
      ANGL(3)=COS(F)
      RETURN
      END
      FUNCTION GAMMAC(J,A,U,AM,AMC,AMF,RADNCL)
      COMMON /BLSEC/ISEC
C ISEC=2 - COOTBETCTB[ET [[ET[ [O[[AP[EPHO[O [C[[CKAH[[ [PA[MEHTOB C A.G
      IF(ISEC.NE.2.OR.J.LE.6) GAMMAC=
     *GAMMA1(J,A,U,AM,AMC,AMF,RADNCL)
      IF(ISEC.EQ.2.AND.J.GT.6) GAMMAC=
     *GAMMA2(J,A,U,AM,AMC,AMF,RADNCL)
      RETURN
      END
      FUNCTION TKIN(L,AMC,AMF)
      COMMON /BLSEC/ISEC
C ISEC=2 -COOTBETCTB[ET [[ET[ [O[[AP[EPHO[O [C[[CKAH[[ [PA[MEHTOB C A.GT
      IF(ISEC.NE.2.OR.L.LE.6) TKIN=TKIN1(L,AMC,AMF)
      IF(ISEC.EQ.2.AND.L.GT.6) TKIN=TKIN2(L,AMC,AMF)
      RETURN
      END
      FUNCTION GAMMA1(J,A,U,AM,AMC,AMF,RADNCL)
C  B[[[C[EH[E A[CO[[THO[ [[P[H[ [C[[CKAH[[ [PA[MEHTA COPTA J
C  ( B E[. (1./[P[BE[.[OCT.[[AHKA) M[B )
C  [[P[HA [E[EH[[ B TEX [E E[[H[[AX
       COMMON /BL1009/AFJ(35)/BL1015/RJ(35)/BL1014/GAM(35)
      COMMON /BL1011/VJ(35) /BL1005/AJ(35)/BL1006/ZJ(35)
      COMMON /BL1010/ZFJ(35)
      IF(RJ(J).LE.0.OR.U.LE.0.) GO TO 12
      PER=2.*SQRT(AM*A*U)
      IF(J.EQ.33) GO TO 10
      RN=1.5
      CC=0.2
      IF(J.GT.2.AND.J.LE.6)  CC=0.1
      IF(J.GT.6) CC=    (AJ(J)/AFJ(J))**0.6666667
       IF(J-1)1,2,1
2      ALFA=.76+2.2/AFJ(1)**.33333
       BETA=(2.12/AFJ(1)**.66667-.05)/ALFA
      GO TO 3
1      ALFA=1.+CC
       BETA=0.
       GO TO 3
3      Q1=AMC*AFJ(J)
       Q2=AMC*AFJ(J)*RJ(J)
       Q3=(GAM(J)*AFJ(J)**.66667)*(ALFA/Q1**2)*(AFJ(J)/(AFJ(J)+AJ(J)))
      Q3=Q3*(3.1416*RN**2)/(2.*41.5*3.1416**2)
       Q4=(2.*BETA*Q1-3.)/2.+Q2
       Q5=(2.*BETA*Q1-3.)*(SQRT(Q2)-.5)+2.*Q2
      IF(PER-160.) 20,20,21
 20   PEX1=Q4*EXP(-PER)
      GO TO 22
 21   PEX1=0
 22   PP2=PER-2.*SQRT(Q2)
      IF(PP2-160.) 23,23,24
 23   PEX2=Q5*EXP(-PP2)
      GO TO 25
 24   PEX2=0.
 25   GAMMA1=Q3*(PEX1+PEX2)
      RETURN
 10   Q1=2.*SQRT(AMF*A*RJ(J))
      Q2=1./(4.*3.1416)
      IF(PER-160.) 30,30,31
 30   PEX1=EXP(-PER)
      GO TO 32
 31   PEX1=0.
 32   PP2=PER-Q1
      IF(PP2-160.) 33,33,34
 33   PEX2=(Q1-1.)*EXP(-PP2)
      GO TO 35
 34   PEX2=0.
 35   GAMMA1=(Q2/(AMF*A))*(PEX1+PEX2)
      TFIS=10000.
      GFIS=0.21*940./TFIS
      IF(GAMMA1.GT.GFIS) GAMMA1=GFIS
      RETURN
 12   GAMMA1=0.
      RETURN
       END
       FUNCTION TKIN1(L,AMC,AMF)
C  PO[[[P[[ K[H.[HEP[[[ B[[ETEB[E[O [PA[MEHTA COPTA L B C.[.[. (M[B)
       COMMON /BL1009/AFJ(35)/BL1015/RJ(35)/BL1011/VJ(35)
       COMMON /BL1005/AJ(35)
       RJL=RJ(L)
       RB=4.*AMC*AFJ(L)*RJL
       PP1=SQRT(RB)
 5     B1=RNDM(-1)
       IF(B1.LE.0.OR.B1.GT.1) B1=RNDM(-1)
       IF(PP1-160.) 21,21,22
 21    PEX1=EXP(-PP1)
       GO TO 23
 22    PEX1=0.
 23    RK=1.+(1./PP1)*ALOG(B1+(1.-B1)*PEX1)
       IF(L-1) 1,2,1
2      BETA=(2.12/AFJ(1)**0.66667-0.05)/(0.76+2.2/AFJ(1)**0.33333)
       Q1=1.+BETA/RJL
       Q2=Q1*SQRT(Q1)
       FRK=(((3.*SQRT(3.))/2.)/Q2)*(Q1*RK-RK**3)
       GO TO 3
1     FRK=((3.*SQRT(3.))/2.)*(RK-RK**3)
       GO TO 3
 3     B2=RNDM(-1)
       IF(B2-FRK) 4,4,5
4      TKIN1=   RJL*(1.-RK**2)+VJ(L)
       RETURN
      END
       FUNCTION TKIN2(L,AMC,AMF)
C  PO[[[P[[ K[H.[HEP[[[ B[[ETEB[E[O [PA[MEHTA COPTA L B C.[.[. (M[B)
       COMMON /BL1009/AFJ(35)/BL1015/RJ(35)/BL1011/VJ(35)
       COMMON /BL1005/AJ(35) /BL1010/ZFJ(35)
C  [[[T[BAETC[ BO[MO[HOCT[ [O[[AP[EPHO[O [C[[CKAH[[ [PA[MEHTOB.
      DCOL=1.
      DALF=0.869+9.91/ZFJ(L)
      RJM=RJ(L)+VJ(L)
      EX=VJ(L)+DCOL
      Q1=AMC*AFJ(L)
      EY=VJ(L)/5.
      PER=2.*SQRT(Q1*RJM)
      IF(RJ(L).LE.DCOL) GO TO 1
      EM=VJ(L)+(SQRT(1.+4.*Q1*RJ(L))-1.)/(2.*Q1)
      IF(EM-EX) 1,1,2
 2    FM=(EM-VJ(L))*EXP(2.*SQRT(Q1*(RJM-EM))-PER)
      GO TO 3
 1    EM=RJM-Q1/(DALF*DALF)
      IF(EM.GT.EX) EM=EX
      IF(EM.LE.EY) EM=EY
      FM=DCOL*EXP(2.*SQRT(Q1*(RJM-EM))+DALF*(EM-EX)-PER)
 3    ER=EY+(RJM-EY)*RNDM(-1)
      FR=FM*RNDM(-1)
      IF(ER.GT.EX) FTR=(ER-VJ(L))*EXP(2.*SQRT(Q1*(RJM-ER))-PER)
      IF(ER.LE.EX) FTR=DCOL*EXP(2.*SQRT(Q1*(RJM-ER))+DALF*(ER-EX)-PER)
      IF(FR-FTR) 4,4,3
 4    TKIN2=ER
       RETURN
      END
      FUNCTION GAMMA2(J,A,U,AM,AMC,AMF,RADNCL)
C  B[[[C[EH[E A[CO[[THO[ [[P[H[ [C[[CKAH[[ [PA[MEHTA COPTA J
C  ( B E[. (1./[P[BE[.[OCT.[[AHKA) M[B )
C  [[P[HA [E[EH[[ B TEX [E E[[H[[AX
       COMMON /BL1009/AFJ(35)/BL1015/RJ(35)/BL1014/GAM(35)
      COMMON /BL1011/VJ(35) /BL1005/AJ(35)/BL1006/ZJ(35)
      COMMON /BL1010/ZFJ(35)
C  [[[T[BAETC[ BO[MO[HOCT[ [O[[AP[EPHO[O [C[[CKAH[[ [PA[MEHTOB.
      IF(U.LE.0.) GO TO 12
      PER=2.*SQRT(AM*A*U)
      IF(J.EQ.33) GO TO 10
      DCOL=1.
      DALF=0.869+9.91/ZFJ(J)
      RJM=RJ(J)+VJ(J)
      EY=VJ(J)/5.
      IF(RJM.LE.EY.OR.J.EQ.1) GO TO 12
       CC=    (AJ(J)/AFJ(J))**0.6666667
       ALFA=1.+CC
       Q1=AMC*AFJ(J)
      BETA=DCOL
       Q3=(GAM(J)*AFJ(J)**.66667)*(ALFA/Q1**2)*(AFJ(J)/(AFJ(J)+AJ(J)))
      Q3=Q3*(3.1416*RADNCL**2)/(2.*41.5*3.1416**2)
      RJ(J)=RJ(J)-DCOL
      GAMMA3=0.
      IF(RJ(J).LE.0.) GO TO 8
       Q2=AMC*AFJ(J)*RJ(J)
       Q4=(2.*BETA*Q1-3.)/2.+Q2
       Q5=(2.*BETA*Q1-3.)*(SQRT(Q2)-.5)+2.*Q2
      GAMMA3=Q3*(Q4*EXP(-PER)+Q5*EXP(2.*SQRT(Q2)-PER))
 8    RJ(J)=RJ(J)+DCOL
      EX=VJ(J)+DCOL
      EM=RJM-Q1/(DALF*DALF)
      IF(EM.GE.EX) GO TO 13
      IF(EM.LT.EX.AND.EM.GT.EY) GO TO 16
      IF(EM.LE.EY) GO TO 17
 13   SQ=SQRT(Q1/(RJM-EX))
      F1X=DALF-SQ
      F2X=SQ/(2.*(RJM-EX))
      IF(F1X.GE.(0.5*F2X)) GO TO 14
      CSI=0.48/SQRT(0.5*F2X)
      F1CSI=DALF-SQRT(Q1/(RJM-EX+CSI))
      GAMMA4=Q3*Q1*Q1*DCOL*SQRT(6.2832/F2X)*
     &EXP(2.*SQRT(Q1*(RJM-EX))-PER-F1CSI*CSI)
      GO TO 15
 14   CSI=0.693/F1X
      SQCSI=SQRT(Q1/(RJM-EX+CSI))
      F2CSI=SQCSI/(2.*(RJM-EX+CSI))
      GAMMA4=Q3*2*Q1*Q1*(DCOL/F1X)*
     &EXP(2.*SQRT(Q1*(RJM-EX))-PER-F2CSI*CSI*CSI/2.)
      GO TO 15
 16   SQ=SQRT(Q1/(RJM-EM))
      F2M=SQ/(2.*(RJM-EM))
      CSI=0.48/SQRT(0.5*F2M)
      F1CSI=DALF-SQRT(Q1/(RJM-EM+CSI))
      GAMMA4=Q3*2*Q1*Q1*DCOL*SQRT(6.2832/F2M)*
     &EXP(DALF*(EM-EX)+2.*SQRT(Q1*(RJM-EM))-PER-F1CSI*CSI)
      GO TO 15
 17   GAMMA4=0.
      GO TO 15
 15   GAMMA2=GAMMA3+GAMMA4
      RETURN
 10   IF(RJ(J).LE.0.) GO TO 12
      Q1=2.*SQRT(AMF*A*RJ(J))
      Q2=1./(4.*3.1416)
      GAMMA2=(Q2/(AMF*A))*((Q1-1.)*EXP(Q1-PER)+EXP(-PER))
      RETURN
 12   GAMMA2=0.
      RETURN
       END
      SUBROUTINE FILLEV
C  [A[O[HEH[E BO[[[[[EHH[X [POBHE[ [E[K[X [PA[MEHTOB
      COMMON /EXILEV/ENLEV(20,35),SPLEV(20,35)
      ENLEV(1,10)=3.56
      SPLEV(1,10)=1.
      ENLEV(1,11)=0.48
      SPLEV(1,11)=2.
      ENLEV(1,12)=0.98
      SPLEV(1,12)=3.
      ENLEV(1,13)=0.43
      SPLEV(1,13)=2.
      ENLEV(1,16)=3.37
      SPLEV(1,16)=5.
      ENLEV(2,16)=5.96
      SPLEV(2,16)=8.
      ENLEV(3,16)=6.18
      SPLEV(3,16)=1.
      ENLEV(4,16)=6.26
      SPLEV(4,16)=5.
      ENLEV(1,18)=0.72
      SPLEV(1,18)=3.
      ENLEV(2,18)=1.74
      SPLEV(2,18)=1.
      ENLEV(3,18)=2.15
      SPLEV(3,18)=3.
      ENLEV(4,18)=3.59
      SPLEV(4,18)=5.
      ENLEV(1,19)=2.13
      SPLEV(1,19)=2.
      ENLEV(2,19)=4.44
      SPLEV(2,19)=6.
      ENLEV(3,19)=5.02
      SPLEV(3,19)=4.
      ENLEV(4,19)=6.76
      SPLEV(4,19)=10.
      ENLEV(5,19)=7.29
      SPLEV(5,19)=6.
      ENLEV(6,19)=7.98
      SPLEV(6,19)=4.
      ENLEV(7,19)=8.56
      SPLEV(7,19)=6.
      ENLEV(1,20)=0.95
      SPLEV(1,20)=5.
      ENLEV(2,20)=1.67
      SPLEV(2,20)=5.
      ENLEV(3,20)=2.65
      SPLEV(3,20)=4.
      ENLEV(1,21)=2.00
      SPLEV(1,21)=2.
      ENLEV(2,21)=4.32
      SPLEV(2,21)=6.
      ENLEV(3,21)=4.80
      SPLEV(3,21)=4.
      ENLEV(4,21)=6.34
      SPLEV(4,21)=2.
      ENLEV(5,21)=6.48
      SPLEV(5,21)=8.
      ENLEV(6,21)=6.90
      SPLEV(6,21)=6.
      ENLEV(7,21)=7.50
      SPLEV(7,21)=4.
      ENLEV(8,21)=8.10
      SPLEV(8,21)=4.
      ENLEV(9,21)=8.42
      SPLEV(9,21)=6.
      ENLEV(10,21)=8.66
      SPLEV(10,21)=8.
      ENLEV(1,22)=4.44
      SPLEV(1,22)=5.
      ENLEV(1,23)=3.09
      SPLEV(1,23)=2.
      ENLEV(2,23)=3.68
      SPLEV(2,23)=4.
      ENLEV(3,23)=3.85
      SPLEV(3,23)=6.
      ENLEV(1,24)=6.09
      SPLEV(1,24)=3.
      ENLEV(2,24)=6.69
      SPLEV(2,24)=8.
      ENLEV(3,24)=6.96
      SPLEV(3,24)=6.
      ENLEV(4,24)=7.34
      SPLEV(4,24)=5.
      ENLEV(1,26)=2.31
      SPLEV(1,26)=1.
      ENLEV(2,26)=3.95
      SPLEV(2,26)=3.
      ENLEV(3,26)=4.92
      SPLEV(3,26)=1.
      ENLEV(4,26)=5.11
      SPLEV(4,26)=5.
      ENLEV(5,26)=5.69
      SPLEV(5,26)=3.
      ENLEV(6,26)=5.83
      SPLEV(6,26)=7.
      ENLEV(7,26)=6.20
      SPLEV(7,26)=3.
      ENLEV(8,26)=6.44
      SPLEV(8,26)=7.
      ENLEV(9,26)=7.03
      SPLEV(9,26)=5.
      ENLEV(1,27)=5.28
      SPLEV(1,27)=8.
      ENLEV(2,27)=6.32
      SPLEV(2,27)=4.
      ENLEV(3,27)=7.22
      SPLEV(3,27)=10.
      ENLEV(4,27)=7.57
      SPLEV(4,27)=8.
      ENLEV(5,27)=8.31
      SPLEV(5,27)=2.
      ENLEV(6,27)=8.57
      SPLEV(6,27)=4.
      ENLEV(7,27)=9.15
      SPLEV(7,27)=14.
      ENLEV(8,27)=9.79
      SPLEV(8,27)=14.
      ENLEV(9,27)=10.00
      SPLEV(9,27)=8.
      ENLEV(1,28)=0.12
      SPLEV(1,28)=1.
      ENLEV(2,28)=0.30
      SPLEV(2,28)=7.
      ENLEV(3,28)=0.40
      SPLEV(3,28)=3.
      ENLEV(1,29)=5.22
      SPLEV(1,29)=8.
      ENLEV(2,29)=6.18
      SPLEV(2,29)=4.
      ENLEV(3,29)=6.83
      SPLEV(3,29)=10.
      ENLEV(4,29)=7.28
      SPLEV(4,29)=8.
      ENLEV(1,30)=6.10
      SPLEV(1,30)=8.
      ENLEV(2,30)=6.92
      SPLEV(2,30)=5.
      ENLEV(3,30)=7.12
      SPLEV(3,30)=3.
      ENLEV(1,31)=0.87
      SPLEV(1,31)=2.
      ENLEV(2,31)=3.06
      SPLEV(2,31)=2.
      ENLEV(3,31)=3.84
      SPLEV(3,31)=6.
      ENLEV(1,32)=1.98
      SPLEV(1,32)=5.
      ENLEV(2,32)=3.57
      SPLEV(2,32)=10.
      ENLEV(3,32)=3.92
      SPLEV(3,32)=5.
      ENLEV(4,32)=4.46
      SPLEV(4,32)=3.
      ENLEV(5,32)=5.10
      SPLEV(5,32)=7.
      ENLEV(6,32)=5.33
      SPLEV(6,32)=13.
      ENLEV(7,32)=5.53
      SPLEV(7,32)=5.
      ENLEV(8,32)=6.20
      SPLEV(8,32)=3.
      ENLEV(9,32)=6.38
      SPLEV(9,32)=12.
      ENLEV(10,32)=6.88
      SPLEV(10,32)=1.
      RETURN
      END
      SUBROUTINE FISION(ENEXT,ATWGHT,CHARGE,PFX,PFY,PFZ,KSTART)
C  FISSION  (IN - MEV, OUT - GEV)
      COMMON/BLOKC/SPTU(10,500)/BLANGL/ANGL(4)
      COMMON/BLFIS/FIS,EF,SNT(6,6)
      UT=ENEXT
      AT=ATWGHT
      ZT=CHARGE
      VFX=PFX/(940.*AT)
      VFY=PFY/(940.*AT)
      VFZ=PFZ/(940.*AT)
      SNT(1,2)=ENEXT*0.001
      SNT(2,2)=ATWGHT
      SNT(3,2)=CHARGE
      SNT(4,2)=PFX*0.001
      SNT(5,2)=PFY*0.001
      SNT(6,2)=PFZ*0.001
      JJ=0
 1    CONTINUE
      JJ=JJ+1
      IF(JJ.GT.100) GO TO 10
      AR=AFIS(UT,AT,ZT,EF)
      IA1=INT(AR+0.5)
      ZR=ZFIS(AT,ZT,AR)
      IZ1=INT(ZR+0.5)
      AF1=IA1
      ZF1=IZ1
      AF2=AT-AF1
      ZF2=ZT-ZF1
      IA2=INT(AF2+0.5)
      IZ2=INT(ZF2+0.5)
      EKF=FISKIN(UT,AT,ZT,AF1,AF2,ZF1,ZF2)
      AMAS1=940.*AF1
      AMAS2=940.*AF2
      EK1=AF2*EKF/AT
      EK2=AF1*EKF/AT
      CALL DELAM(AT,ZT,DL0,DSH0,BAR0)
      CALL DELAM(AF1,ZF1,DL1,DSH1,BAR1)
      CALL DELAM(AF2,ZF2,DL2,DSH2,BAR2)
      ET=UT+DL0-DL1-DL2
C      !WRITE 100 2,AT,ZT,AF1,ZF1,AF2,ZF2,UT,DL0,DL1,DL2,ET,EKF
C 2    FORMAT(2X,'AT,ZT=',2F6.1,' AF1,ZF1=',2F6.1,' AF2,ZF2=',2F6.1,
C     *' UT,DL0=',2F6.1,' DL1,DL2=',2F6.1,' ET,EKF=',2F6.1)
      IF(ET.LT.EKF) GO TO 1
      U1=0.001*(ET-EKF)*AF1/AT
      U2=0.001*(ET-EKF)*AF2/AT
      CALL ISANGL
      VM=SQRT(2.*EK1/AMAS1)
      VX1=VM*ANGL(4)*ANGL(3)
      VY1=VM*ANGL(4)*ANGL(2)
      VZ1=VM*ANGL(1)
      PX1=AMAS1*(VX1+VFX)*0.001
      PY1=AMAS1*(VY1+VFY)*0.001
      PZ1=AMAS1*(VZ1+VFZ)*0.001
      TEMP=-1.*(AMAS1/AMAS2)
      VX2=TEMP*VX1
      VY2=TEMP*VY1
      VZ2=TEMP*VZ1
      PX2=AMAS2*(VX2+VFX)*0.001
      PY2=AMAS2*(VY2+VFY)*0.001
      PZ2=AMAS2*(VZ2+VFZ)*0.001
C  BEFORE EVAPORATION
      SNT(1,3)=U1
      SNT(2,3)=AF1
      SNT(3,3)=ZF1
      SNT(4,3)=PX1
      SNT(5,3)=PY1
      SNT(6,3)=PZ1
      SNT(1,4)=U2
      SNT(2,4)=AF2
      SNT(3,4)=ZF2
      SNT(4,4)=PX2
      SNT(5,4)=PY2
      SNT(6,4)=PZ2
C  AFTER EVAPORATION
C      !WRITE 100 3,KSTART
C 3    FORMAT(2X,' BEFORE EVAP.FRAGM. KSTART=',I4)
      CALL EVAPIN(U1,AF1,ZF1,PX1,PY1,PZ1,KSTART)
      SNT(1,5)=U1
      SNT(2,5)=AF1
      SNT(3,5)=ZF1
      SNT(4,5)=PX1
      SNT(5,5)=PY1
      SNT(6,5)=PZ1
C      !WRITE 100 4,KSTART
C 4    FORMAT(2X,'AFT.EVAP.FIRST, BEF.EVAP.SEC.FRAGM. KSTART=',I4)
      CALL EVAPIN(U2,AF2,ZF2,PX2,PY2,PZ2,KSTART)
      SNT(1,6)=U2
      SNT(2,6)=AF2
      SNT(3,6)=ZF2
      SNT(4,6)=PX2
      SNT(5,6)=PY2
      SNT(6,6)=PZ2
C      !WRITE 100 5,KSTART
C 5    FORMAT(2X,' AFTER EVAP.FRAGM. KSTART=',I4)
      RETURN
 10   continue !WRITE 100 11,ET,EKF,JJ,UT,AT,ZT
 11   FORMAT(2X,'ERROR IN FISION: ET=',F6.1,'.LT.EKF=',F6.1,
     *' JJ=',I4,' UT,AT,ZT=',3F6.1)
      RETURN
      END
      FUNCTION YIELDF(A0,A,F,AS,A1,A2,SIGS,SIG1,SIG2)
      XSIM=EXP(-0.5*(A-AS)**2/SIGS**2)
      XASIM=EXP(-0.5*(A-A2)**2/SIG2**2)
     *     +EXP(-0.5*(A-(A0-A2))**2/SIG2**2)
     * +0.5*EXP(-0.5*(A-A1)**2/SIG1**2)
     * +0.5*EXP(-0.5*(A-(A0-A1))**2/SIG1**2)
      IF(F.GT.1000.) YIELDF=XSIM
      IF(F.LT.0.001) YIELDF=XASIM
      IF(F.GE.0.001.AND.F.LE.1000) YIELDF=F*XSIM+XASIM
      RETURN
      END
      FUNCTION ZFIS(A,Z,AF)
      COMMON /BLFIS0/F,FH,AS,A1,A2,SIGS,SIG1,SIG2,ZP,SIGZ,EAVR,EMAX
      IF(AF.GE.134.) DZ=-0.45
      IF(AF.LE.(A-134.)) DZ=0.45
      IF(AF.GT.(A-134.).AND.AF.LT.134.) DZ=-0.45*(AF-0.5*A)/(134-0.5*A)
      ZP=AF*Z/A+DZ
      SIGZ=0.6
 3    CALL RANNOR(X,Y)
      XZ=ZP+SIGZ*X
      IF(XZ.LT.1.OR.XZ.GT.(Z-1).OR.XZ.GT.AF) XZ=ZP+SIGZ*Y
      IF(XZ.LT.1.OR.XZ.GT.(Z-1).OR.XZ.GT.AF) GO TO 3
      ZFIS=INT(XZ+0.5)
      RETURN
      END
       SUBROUTINE EVAPIN(ENEXT,ATWGHT,CHARGE,PNX,PNY,PNZ,KSTART)
C   EVAPORATION WITHOUT FISSION FOR SUBR.-FISION
       COMMON /BL0999/RNCL/BL1000/AM,AMF,AMCB(35)/AMNUCL/AMNEV(3)
       COMMON /BL1001/T1Y(130)/BL1002/T2XY(200)/BL1005/AJ(35)
     */BL1006/ZJ(35)/BL1014/GAM(35)/BL1011/VJ(35)/BL1015/RJ(35)
      COMMON /BL1003/U,A,Z/BL1009/AFJ(35)/BL1010/ZFJ(35)
      COMMON/BLOKC/SPTU(10,500)/BLANGL/ANGL(4)/BL1100/AMS,AMFS
      COMMON/BLFIS/FIS,EF,SNT(6,6) /BLCAN/ICAN
      COMMON /EXILEV/ ENLEV(20,35),SPLEV(20,35)
       DIMENSION GJ(35),BJ(35),GG(35),VN(3),PN(3),PNL(3),PP(3),PPL(3)
      DIMENSION DSH2(35)
      FIS=1.
      ILEVEL=3
      CALL DATEVP(ILEVEL)
      U=ENEXT*1000.
      A=ATWGHT
      Z=CHARGE
      AMNEV(1)=0.
      AMNEV(2)=0.
      AMNEV(3)=0.
      REMN=940.*A
      PNL(1)=PNX*1000.
      PNL(2)=PNY*1000.
      PNL(3)=PNZ*1000.
      ENL=SQRT(PNL(1)**2+PNL(2)**2+PNL(3)**2+REMN**2)
      VN(1)=PNL(1)/ENL
      VN(2)=PNL(2)/ENL
      VN(3)=PNL(3)/ENL
      DO 20 K=KSTART,500
      IN=35
      DO 4 I=1,IN
      VJ(I)=0.
      AMCB(I)=0.
C      DSH2(I)=0.
 4    RJ(I)=-1000.
      CALL DELAM(A,Z,DL1,DSHEL1,BAR1)
      AM=AMS
C      IF(U.GT.0.) AM=AMS*(1.+(1.-EXP(-0.06*U))*DSHEL1/U)
      DO 6 I=1,IN
      IF(I.EQ.33) GO TO 6
      AFJ(I)=A-AJ(I)
      ZFJ(I)=Z-ZJ(I)
        IF(ICAN.EQ.0.AND.I.GT.6) GO TO 6
        IF(ICAN.EQ.1.AND.I.GT.32) GO TO 6
        IF(ICAN.EQ.2.AND.I.GT.6.AND.I.LE.32) GO TO 6
        IF(AFJ(I).LT.ZFJ(I)) GO TO 6
        IF(AFJ(I).LT.AJ(I).OR.ZFJ(I).LT.ZJ(I)) GO TO 6
      JAF=INT(AFJ(I)+0.5)
      JZF=INT(ZFJ(I)+0.5)
      JNF=JAF-JZF
      ODD=11.*(2+2*(JNF/2)-JNF+2*(JZF/2)-JZF)/SQRT(AFJ(I))
      IF(RNCL.GE.0.1) RADNCL=RNCL
      IF(RNCL.LT.0.1)
     *RADNCL=2.173*(1+0.006103*ZJ(I)*ZFJ(I))/(1+0.009443*ZJ(I)*ZFJ(I))
      VJ(I)=COLOMB(I,RADNCL)
      CALL DELAM(AFJ(I),ZFJ(I),DL2,DSHEL2,BAR2)
C      DSH2(I)=DSHEL2
      CALL DELAM(AJ(I),ZJ(I),DL3,DSHEL3,BAR3)
      BJ(I)=DL2+DL3-DL1
c ******** in accordance with DELA8
      IF(A.LE.55.) RJ(I)=U-(BJ(I)+VJ(I))-ODD
      TCOEF=1.-(A-55.)/10.
      IF(A.GT.55.AND.A.LT.65.) RJ(I)=U-(BJ(I)+VJ(I))-ODD*TCOEF
      IF(A.GE.65.) RJ(I)=U-(BJ(I)+VJ(I))
c ********
C     RJ(I)=U-(BJ(I)+VJ(I))-ODD
C     RJ(I)=U-(BJ(I)+VJ(I))
C      IF(RJ(I).GT.0.)
C     *AMCB(I)=AMS*(1.+(1.-EXP(-0.06*RJ(I)))*DSHEL2/RJ(I))
 6    CONTINUE
      RJ(33)=-1001.
       DO 7 I=1,IN
      IF(RNCL.GE.0.1) RADNCL=RNCL
      IF(RNCL.LT.0.1)
     *RADNCL=2.173*(1+0.006103*ZJ(I)*ZFJ(I))/(1+0.009443*ZJ(I)*ZFJ(I))
      AMC=AM
C      AMC=AMCB(I)
      GJ(I)=GAMMAC(I,A,U,AM,AMC,AMF,RADNCL)
      RJI=RJ(I)
      GAMI=GAM(I)
      DO 71 IL=1,20
      IF(SPLEV(IL,I).LT.0.1) GO TO 70
      GAM(I)=SPLEV(IL,I)*AJ(I)
      RJ(I)=RJI-ENLEV(IL,I)
C      IF(RJ(I).GT.0.)
C     *AMC=AMS*(1.+(1.-EXP(-0.06*RJ(I)))*DSH2(I)/RJ(I))
      GJ(I)=GJ(I)+GAMMAC(I,A,U,AM,AMC,AMF,RADNCL)
 71   CONTINUE
 70   RJ(I)=RJI
      GAM(I)=GAMI
7     CONTINUE
      G=0.
      DO 10 I=1,IN
10    G=G+GJ(I)
        IF(G) 11,11,12
 11   PNX=PNL(1)*0.001
      PNY=PNL(2)*0.001
      PNZ=PNL(3)*0.001
      ENEXT=U*0.001
      ATWGHT=A
      CHARGE=Z
      REMN=940.*A
      SPTU(9,K)=REMN*0.001
      SPTU(8,K)=Z
      CALL PINT(PNL,CT,ST,CF,SF,T,REMN)
      SPTU(4,K)=CT
      SPTU(5,K)=SF
      SPTU(6,K)=CF
      SPTU(7,K)=T*0.001
      SPTU(1,K)=1.
      SPTU(3,K)=FIS
      SPTU(2,K)=2.
      KSTART=K+1
      RETURN
 12   CONTINUE
      DO 13 J=2,IN
13    GJ(J)=GJ(J-1)+GJ(J)
      BB=RNDM(-1)
      B=BB*G
      DO 14 J=1,IN
        IF(B-GJ(J)) 15,14,14
15    LM=J
      GO TO 16
 14   CONTINUE
 16     IF(LM-33) 18,17,18
 17   continue !WRITE 100 300
 300  FORMAT(2X,' ERROR IN EPAPIN - FALSE FISSION')
      GO TO 11
 18   AMC=AM
C18    AMC=AMCB(LM)
      EP1=TKIN(LM,AMC,AMF)
      EP2=ZJ(LM)
      EP3=940.*AJ(LM)
      U=U-BJ(LM)-EP1
      A=AFJ(LM)
      Z=ZFJ(LM)
C  VPM-OTHOC[TE[[HA[ CKOPOCT[ [PA[MEHTA
      VPM=SQRT((2.*EP1)/(EP3*AFJ(LM)/(AFJ(LM)+AJ(LM))))
      CALL ISANGL
C  [M[[[[C [PA[MEHTA B C.[.[.
      PP(1)=VPM*ANGL(4)*ANGL(3)*EP3/(1.+AJ(LM)/A)
      PP(2)=VPM*ANGL(4)*ANGL(2)*EP3/(1.+AJ(LM)/A)
      PP(3)=VPM*ANGL(1)*EP3/(1.+AJ(LM)/A)
      PN(1)=-PP(1)
      PN(2)=-PP(2)
      PN(3)=-PP(3)
      EP=SQRT(PP(1)**2+PP(2)**2+PP(3)**2+EP3**2)
      EN=SQRT(PN(1)**2+PN(2)**2+PN(3)**2+(940.*A)**2)
      CALL CLPV(PP,VN,PPL,EP)
      CALL CLPV(PN,VN,PNL,EN)
      ENL=SQRT(PNL(1)**2+PNL(2)**2+PNL(3)**2+(940.*A)**2)
      VN(1)=PNL(1)/ENL
      VN(2)=PNL(2)/ENL
      VN(3)=PNL(3)/ENL
      CALL PINT(PPL,CT,ST,CF,SF,T,EP3)
      SPTU(4,K)=CT
      SPTU(5,K)=SF
      SPTU(6,K)=CF
      SPTU(7,K)=T*0.001
      SPTU(8,K)=EP2
      SPTU(9,K)=EP3*0.001
      SPTU(2,K)=2.
20    CONTINUE
      !WRITE 100 21,U,A,Z
21    FORMAT(35X,37HMASSIV SPT EXCEEDED AFTER EVAPORATION/40X,2HU=,
     *F10.5,4H  A=,F5.1,4H  Z=,F4.1)
       RETURN
       END
      FUNCTION AFIS(E,A,Z,EF)
C IMPROVED VERSION TAKING INTO ACCOUNT DEEXCIT. BEFORE FISSION
      COMMON /BLFIS0/F,FH,AS,A1,A2,SIGS,SIG1,SIG2,ZP,SIGZ,EAVR,EMAX
      A1=134.
      A2=141.
      AS=A/2.
      IZ=INT(Z+0.5)
      IF(A.GE.235.) SIG2=5.6
      IF(A.LT.235.) SIG2=5.6+0.096*(A-235.)
      SIG1=0.5*SIG2
C     IF(E.LE.15.) SIGS=EXP((0.2412*E+0.8252)/2.)
      SIGS=EXP((0.01106*E+4.2772)/2.)
      IF(SIGS.GT.20.) SIGS=20.
      YA=2.*EXP(-0.5*(A2-AS)**2/SIG2**2)+EXP(-0.5*(A1-AS)**2/SIG1**2)
      YS=EXP(-0.5*(AS-(A1+A2)/2.)**2/SIGS**2)
      IF(IZ.GE.90) GO TO 1
      IF(IZ.EQ.89) GO TO 2
      IF(IZ.LE.88.AND.IZ.GE.82) GO TO 3
      IF(IZ.LT.82) GO TO 4
 1    IF(E.LE.16.25) F=EXP(0.5385*E-9.9564)
      IF(E.GT.16.25) F=EXP(0.09197*E-2.7003)
      GO TO 60
 2    F=EXP(0.09197*E-1.0808)
      GO TO 60
 3    X=EF-7.5
      IF(X.LT.0.) X=0.
      F=EXP(0.09197*(E-X)-1.0808)
      T1=1.03*F-YA
      T2=1.-YS*F
      IF(T2.LE.0.0001) T2=0.0001
      IF(T1.LE.0.0001) T1=0.0001
      FH=T1/T2
      IF(A.LT.227.) FH=FH*EXP(0.3*(227-A))
      GO TO 5
 4    FH=1001.
      GO TO 5
 60   T1=1.03*F-YA
      T2=1.-YS*F
      IF(T2.LE.0.0001) T2=0.0001
      IF(T1.LE.0.0001) T1=0.0001
      FH=T1/T2
 5    CONTINUE
      C2A=A2+3.72*SIG2
      C2S=AS+3.72*SIGS
      C2=AMAX1(C2A,C2S)
      IF(FH.GT.1000.) C2=C2S
      IF(FH.LT.0.001) C2=C2A
      C1=A-C2
      IF(C1.LT.30.) C2=A-30.
      IF(C1.LT.30.) C1=30.
      AM1=(AS+A1)/2.
      AM2=(A1+A2)/2.
      XM1=YIELDF(A,AS,FH,AS,A1,A2,SIGS,SIG1,SIG2)
      XM2=YIELDF(A,AM1,FH,AS,A1,A2,SIGS,SIG1,SIG2)
      XM3=YIELDF(A,A1,FH,AS,A1,A2,SIGS,SIG1,SIG2)
      XM4=YIELDF(A,AM2,FH,AS,A1,A2,SIGS,SIG1,SIG2)
      XM5=YIELDF(A,A2,FH,AS,A1,A2,SIGS,SIG1,SIG2)
      PMAX=AMAX1(XM1,XM2,XM3,XM4,XM5)
 6    X=C1+RNDM(-1)*(C2-C1)
      PX=YIELDF(A,X,FH,AS,A1,A2,SIGS,SIG1,SIG2)
      Y=RNDM(-1)
      IF(Y-PX/PMAX)7,7,6
 7    AFIS=INT(X+0.5)
      RETURN
      END
      FUNCTION FISKIN(E,A,Z,AF1,AF2,ZF1,ZF2)
C IMPROVED VERSION TAKING INTO ACCOUNT THE CHANNELS SEPARATELY
      COMMON /BLFIS0/F,FH,AS,A1,A2,SIGS,SIG1,SIG2,ZP,SIGZ,EAVR,EMAX
      AFM=AF1
      IF(AF2.GT.AF1) AFM=AF2
      IF(AFM.LT.(A/2.)) AFM=A-AFM
      EAVR=0.1071*Z*Z/A**0.33333+22.2
      P1=0.5*EXP(-0.5*(AFM-A1)**2/SIG1**2)
      P2=EXP(-0.5*(AFM-A2)**2/SIG2**2)
      PS=FH*EXP(-0.5*(AFM-AS)**2/SIGS**2)
      PAS=P1+P2
      IF(FH.GT.1000.) PAS=0.
      IF(FH.LT.0.001) PS=0.
      PSY=PS/(PAS+PS)
      ISY=0
      BRND=RNDM(-1)
      IF(BRND.LE.PSY) ISY=1
      IF(ISY.EQ.0) SIGE=10.
      IF(ISY.GE.1) SIGE=8.
      PPAS=2.*0.5*SIG1+2.*SIG2
      PPSY=FH*SIGS
      PPP=PPAS+PPSY
      XAS=PPAS/PPP
      XSY=PPSY/PPP
      EAVRSY=EAVR-12.5*XAS
      EAVRAS=EAVR+12.5*XSY
      IF(ISY.EQ.0) GO TO 10
      IF(ISY.GE.1) GO TO 11
 10   EAV0=EAVRAS
      A11=A1-0.7979*SIG1
      A12=A1+0.7979*SIG1
      A21=A2-0.7979*SIG2
      A22=A2+0.7979*SIG2
      EEF11=EEFIS(A,A11,ISY)
      EEF12=EEFIS(A,A12,ISY)
      EEF21=EEFIS(A,A21,ISY)
      EEF22=EEFIS(A,A22,ISY)
      TEMPL=0.5*SIG1*EEF11+0.5*SIG1*EEF12+SIG2*EEF21+SIG2*EEF22
      EMAX=EAV0*PPAS/TEMPL
      GO TO 12
 11   EAV0=EAVRSY
      AS0=AS+0.7979*SIGS
      EEF00=EEFIS(A,AS0,ISY)
      TEMPL=FH*SIGS*EEF00
      EMAX=EAV0*PPSY/TEMPL
 12   EKAFM=EMAX*EEFIS(A,AFM,ISY)
      JJJJ=0
 3    CALL RANNOR(X,Y)
      JJJJ=JJJJ+1
      IF(JJJJ.GT.100) EKIN=EAVR
      IF(JJJJ.GT.100) GO TO 650
      EKIN=EKAFM+SIGE*X
      IF(EKIN.LT.(EAVR-3.72*SIGE).OR.EKIN.GT.(EAVR+3.72*SIGE))
     *EKIN=EKAFM+SIGE*Y
      IF(EKIN.LT.(EAVR-3.72*SIGE).OR.EKIN.GT.(EAVR+3.72*SIGE))
     *GO TO 3
 650  FISKIN=EKIN
      RETURN
      END
      FUNCTION EEFIS(A,AFM,ISY)
C IMPROVED VERSION TAKING INTO ACCOUNT THE CHANNELS SEPARATELY
      IF(ISY.GE.1) GO TO 1
      A00=134.
      B1=23.5
      GO TO 2
 1    B1=5.32
      A00=A/2.
 2    TEMP1=10./A
      A0010=A00+10.
      AS=A/2.
      IF(AFM.GE.AS.AND.AFM.LE.A0010) EE=1.-B1*((AFM-A00)/A)**2
      IF(AFM.GT.A0010) EE=1.-B1*TEMP1**2-2.*TEMP1*B1*(AFM-A0010)/A
      EEFIS=EE
      RETURN
      END
      SUBROUTINE DELAM(ANUCL,ZNUCL,DM,SH,BARIER)
      COMMON /FUSR/ BARR,SMASS,SHLL
      IA=INT(ANUCL+0.5)
      IZ=INT(ZNUCL+0.5)
      IN=IA-IZ
      IF(IZ.LT.0.OR.IN.LT.0) GO TO 8
      ODD=-11.*(1+2*(IN/2)-IN+2*(IZ/2)-IZ)/SQRT(ANUCL)
      IF(IA.LT.65) GO TO 7
      IF(IN.GT.240.OR.IZ.GT.240) GO TO 7
      CM=0.
      SHLL=0.
      CBAR=0.
      NBAR=0
      CALL LYMASM(IZ,IA,CM,CBAR,NBAR)
      IF(NBAR.EQ.1) CBAR=BARR
      DM=CM
      SH=SHLL+ODD
      BARIER=CBAR
c ************* inclusion DELA8
      CALL DELA8(ANUCL,ZNUCL,DM8,SH8)
      DM=DM8
c *************
      RETURN
 7    CONTINUE
      CALL DELA1(ANUCL,ZNUCL,DM1,SH1)
      BAR1=1001.
      IF(IA.GE.65) BAR1=BARIE1(ANUCL,ZNUCL)
      BARIER=BAR1
      DM=DM1
      SH=SH1
c ************* inclusion DELA8
      IF(ANUCL.GT.55.) CALL DELA8(ANUCL,ZNUCL,DM8,SH8)
      TCOEF=1.-(ANUCL-55.)/10.
      IF(ANUCL.GT.55.AND.ANUCL.LT.65.) DM=DM1*TCOEF+DM8*(1.-TCOEF)
c *************
      RETURN
 8    continue !WRITE 100 9,IZ,IN
 9    FORMAT(2X,'ERROR IN DELAM: IZ=',I5,' IN=',I5)
      RETURN
      END
      SUBROUTINE DELA8(X,Y,DELTA8,DSHEL8)
      DSHEL8=0.
      I=INT(X+0.5)
      J=INT(Y+0.5)
      L=I-J
      DELTA8=1001.
      IF(I.LE.1.OR.J.LE.0.OR.L.LE.0) RETURN
      X03=X**0.3333333
      X06=X03*X03
      X13=X06*X06
      Y13=Y**1.3333333
      D=(X-2.*Y)/X
      DR=1-0.62025/X06
      ES=(25.8357-44.2355*D*D)*DR*DR*X06
      EC=0.779*(Y*(Y-1.)/X03)*(1.-1.5849/X06+1.2273/X+1.5772/X13)
      EEX=-0.4323*(Y13/X03)*(1.-0.57811/X03-0.14518/X06+0.49597/X)
      EV=X*(-17.0354+31.4506*D*D)
      DELTA8=8.071*X-0.783*Y+EV+ES+EC+EEX
      RETURN
      END
      SUBROUTINE DELA1(X,Y,DELTAM,DSHEL)
C  B[[[C[EH[E (M-A) CO[[ACTHO KAMEPOH[-57 (B [[[EPO[H[X E[[H[[AX)
C XOPO[EE COOTBETCTB[E [P[ Z.LE.8 [ CO[[ACOBAH[E [P[ [O[[[[X Z
C      ([[[ [E[E[ [C[APEH[[ [[ [E[K[X [[EP [ [X CAM[X)
      COMMON /BL1001/T1Y(130) /BL1002/T2XY(200)
      I=INT(X+0.5)
      J=INT(Y+0.5)
      L=I-J
      DELTAM=1001.
      IF(J.EQ.0) DELTAM=8.071*I
      IF(L.EQ.0) DELTAM=7.289*J
      IF(J.GT.10)  GO TO 11
      IF(I.EQ.1.AND.J.EQ.0) DELTAM=8.071
      IF(J.EQ.0)   GO TO 11
      GO TO (1,2,3,4,5,6,7,8,9,10),J
 1    IF(I.EQ.1.AND.J.EQ.1) DELTAM=7.289
      IF(I.EQ.2.AND.J.EQ.1) DELTAM=13.136
      IF(I.EQ.3.AND.J.EQ.1) DELTAM=14.950
      IF(I.EQ.4.AND.J.EQ.1) DELTAM=25.920
      IF(I.EQ.5.AND.J.EQ.1) DELTAM=33.790
      GO TO 11
 2    IF(I.EQ.3.AND.J.EQ.2) DELTAM=14.931
      IF(I.EQ.4.AND.J.EQ.2) DELTAM=2.425
      IF(I.EQ.5.AND.J.EQ.2) DELTAM=28.150
      IF(I.EQ.6.AND.J.EQ.2) DELTAM=17.593
      IF(I.EQ.7.AND.J.EQ.2) DELTAM=26.111
      IF(I.EQ.8.AND.J.EQ.2) DELTAM=31.599
      GO TO 11
 3    IF(I.EQ.5.AND.J.EQ.3) DELTAM=28.340
      IF(I.EQ.6.AND.J.EQ.3) DELTAM=14.087
      IF(I.EQ.7.AND.J.EQ.3) DELTAM=14.908
      IF(I.EQ.8.AND.J.EQ.3) DELTAM=20.947
      GO TO 11
 4    IF(I.EQ.6.AND.J.EQ.4) DELTAM=18.375
      IF(I.EQ.7.AND.J.EQ.4) DELTAM=15.769
      IF(I.EQ.8.AND.J.EQ.4) DELTAM=4.942
      IF(I.EQ.9.AND.J.EQ.4) DELTAM=11.348
      IF(I.EQ.10.AND.J.EQ.4) DELTAM=12.607
      IF(I.EQ.11.AND.J.EQ.4) DELTAM=20.176
      IF(I.EQ.12.AND.J.EQ.4) DELTAM=25.072
      GO TO 11
 5    IF(I.EQ.9.AND.J.EQ.5) DELTAM=12.415
      IF(I.EQ.10.AND.J.EQ.5) DELTAM=12.052
      IF(I.EQ.11.AND.J.EQ.5) DELTAM=8.668
      IF(I.EQ.12.AND.J.EQ.5) DELTAM=13.370
      IF(I.EQ.13.AND.J.EQ.5) DELTAM=16.562
      IF(I.EQ.14.AND.J.EQ.5) DELTAM=24.230
      GO TO 11
 6    IF(I.EQ.10.AND.J.EQ.6) DELTAM=15.702
      IF(I.EQ.11.AND.J.EQ.6) DELTAM=10.650
      IF(I.EQ.12.AND.J.EQ.6) DELTAM=0.
      IF(I.EQ.13.AND.J.EQ.6) DELTAM=3.125
      IF(I.EQ.14.AND.J.EQ.6) DELTAM=3.020
      GO TO 11
 7    IF(I.EQ.13.AND.J.EQ.7) DELTAM=5.346
      IF(I.EQ.14.AND.J.EQ.7) DELTAM=2.864
      IF(I.EQ.15.AND.J.EQ.7) DELTAM=0.102
      IF(I.EQ.16.AND.J.EQ.7) DELTAM=5.683
      IF(I.EQ.17.AND.J.EQ.7) DELTAM=7.871
      IF(I.EQ.18.AND.J.EQ.7) DELTAM=13.274
      IF(I.EQ.19.AND.J.EQ.7) DELTAM=15.790
      GO TO 11
 8    IF(I.EQ.14.AND.J.EQ.8) DELTAM=8.007
      IF(I.EQ.15.AND.J.EQ.8) DELTAM=2.863
      IF(I.EQ.16.AND.J.EQ.8) DELTAM=-4.737
      IF(I.EQ.17.AND.J.EQ.8) DELTAM=-0.809
      IF(I.EQ.18.AND.J.EQ.8) DELTAM=-0.783
      IF(I.EQ.19.AND.J.EQ.8) DELTAM=3.332
      IF(I.EQ.20.AND.J.EQ.8) DELTAM=3.799
      GO TO 11
 9    IF(I.EQ.17.AND.J.EQ.9) DELTAM=1.952
      IF(I.EQ.18.AND.J.EQ.9) DELTAM=0.873
      IF(I.EQ.19.AND.J.EQ.9) DELTAM=-1.486
      IF(I.EQ.20.AND.J.EQ.9) DELTAM=-0.016
      IF(I.EQ.21.AND.J.EQ.9) DELTAM=-0.046
      GO TO 11
 10   IF(I.EQ.19.AND.J.EQ.10) DELTAM=1.752
      IF(I.EQ.20.AND.J.EQ.10) DELTAM=-7.041
 11   DSHEL=0.
      IF(I.LE.1.OR.J.LE.0.OR.L.LE.0) RETURN
      X03=X**0.3333333
      X06=X03*X03
      X13=X06*X06
      Y13=Y**1.3333333
      D=(X-2.*Y)/X
      DR=1-0.62025/X06
      ES=(25.8357-44.2355*D*D)*DR*DR*X06
      EC=0.779*(Y*(Y-1.)/X03)*(1.-1.5849/X06+1.2273/X+1.5772/X13)
      EEX=-0.4323*(Y13/X03)*(1.-0.57811/X03-0.14518/X06+0.49597/X)
      EV=X*(-17.0354+31.4506*D*D)
      IF(J.LE.0.OR.J.GT.130.OR.L.LE.0.OR.L.GT.200) GO TO 22
      T1=T1Y(J)
      T2=T2XY(L)
      GO TO 23
 22   T1=0.
      T2=0.
 23   CONTINUE
      DELTA0=8.071*X-0.783*Y+EV+ES+EC+EEX+T1+T2
      IF(DELTAM.LT.1000.) GO TO 26
      DSHEL=T1+T2
      DELTAM=DELTA0
      RETURN
 26   DSHEL=DELTAM-DELTA0+T1+T2
      RETURN
      END
      BLOCK DATA B4
C [A[A[TC[ [APAMETP[ OCTAT.[[PA,[C[AP[[[[XC[ [PA[MEHTOB,O[PATHO[O
C CE[EH[[([OCTPOBCK[[-59),O[O[O[.[O[PABK[(KAMEPOH-57)
      COMMON /BL1001/T1Y(130)/BL1002/T2XY(200)/BL1014/GAM(35)
      COMMON /BL0999/RNCL  /BL1100/AMS,AMFS /BLGAM/GAM1(35),GAM2(35)
      COMMON /BL1005/AJ(35) /BL1006/ZJ(35) /BLCAN/ICAN /BLSEC/ISEC
      COMMON /EXILEV/ ENLEV(20,35),SPLEV(20,35)
      DIMENSION TT11(65),TT12(65)
      DIMENSION TT21(50),TT22(50),TT23(50),TT24(50)
      EQUIVALENCE (T1Y(1),TT11(1)),(T1Y(66),TT12(1))
      EQUIVALENCE (T2XY(1),TT21(1)),(T2XY(51),TT22(1)),(T2XY(101),
     *TT23(1)),(T2XY(151),TT24(1))
      DATA GAM2/2.,2.,6.,6.,6.,4.,20.,30.,20.,54.,73.,101.,73.,8.,
     *146.,100.,100.,343.,174.,393.,186.,61.,202.,113.,213.,233.,180.,
     *696.,194.,120.,458.,590.,0.,2790.,216./
      DATA GAM1/2.,2.,6.,6.,6.,4.,20., 6.,20.,18.,28., 40.,28.,8.,
     * 36., 10., 36., 70., 44., 36., 44.,12., 26., 14., 26., 42., 30.,
     * 80., 30., 16.,102., 18.,0.,54.,216./
      DATA GAM/2.,2.,6.,6.,6.,4.,20., 6.,20.,18.,28., 40.,28.,8.,
     * 36., 10., 36., 70., 44., 36., 44.,12., 26., 14., 26., 42., 30.,
     * 80., 30., 16.,102., 18.,0.,54.,216./
      DATA AJ/1.,1.,2.,3.,3.,4.,5.,6.,5.,6.,7.,8.,7.,8.,9.,10.,9.,
     *10.,11.,12.,11.,12.,13.,14.,13.,14.,15.,16.,15.,16.,17.,18.,0.,
     *18.,24./
      DATA ZJ/0.,1.,1.,1.,2.,2.,2.,2.,3.,3.,3.,3.,4.,4.,4.,4.,5.,5.,
     *5.,5.,6.,6.,6.,6.,7.,7.,7.,7.,8.,8.,8.,8.,0.,9.,11./
      DATA RNCL/0.0/,AMS/0.125/,AMFS/0.125/,ICAN/1/,ISEC/1/
      DATA ENLEV/700*0./,SPLEV/700*0./
      DATA TT11/
     *  20.80,  15.80,  21.00,  16.80,  19.80,
     *  16.50,  18.80,  16.50,  18.50,  17.20,
     *  18.26,  15.05,  16.01,  12.04,  13.27,
     *  11.09,  12.17,  10.26,  11.04,   8.41,
     *   9.79,   7.36,   8.15,   5.63,   5.88,
     *   3.17,   3.32,   0.82,   1.83,   0.97,
     *   2.33,   1.27,   2.92,   1.61,   2.91,
     *   1.35,   2.40,   0.89,   1.74,   0.36,
     *   0.95,  -0.65,  -0.04,  -1.73,  -0.96,
     *  -2.87,  -2.05,  -4.05,  -3.40,  -5.72,
     *  -3.75,  -4.13,  -2.42,  -2.85,  -1.01,
     *  -1.33,   0.54,  -0.02,   1.74,   0.75,
     *   2.24,   1.00,   1.98,   0.79,   1.54/
      DATA TT12/
     *   0.39,   1.08,   0.00,   0.78,  -0.35,
     *   0.58,  -0.55,   0.59,  -0.61,   0.59,
     *  -0.35,   0.32,  -0.96,  -0.52,  -2.08,
     *  -2.46,  -3.64,  -1.55,  -0.96,   0.97,
     *   0.88,   2.37,   1.75,   2.72,   1.90,
     *   2.55,   1.46,   1.93,   0.86,   1.17,
     *   0.08,   0.39,  -0.76,  -0.39,  -1.51,
     *  -1.17,  -2.36,  -1.95,  -3.06,  -2.62,
     *  -3.55,  -2.95,  -3.75,  -3.07,  -3.79,
     *  -3.06,  -3.77,  -3.05,  -3.78,  -3.12,
     *  -3.90,  -3.35,  -4.24,  -3.86,  -4.92,
     *  -5.06,  -6.77,  -7.41,  -9.18, -10.16,
     * -11.12,  -9.76,  -9.23,  -7.96,  -7.65/
      DATA TT21/
     *  -8.40, -12.90,  -8.00, -11.90,  -9.20,
     * -12.50, -10.80, -13.60, -11.20, -12.20,
     * -12.81, -15.40, -13.07, -15.80, -13.81,
     * -14.98, -12.63, -13.76, -11.37, -12.38,
     *  -9.23,  -9.65,  -7.64,  -9.17,  -8.05,
     *  -9.72,  -8.87, -10.76,  -8.64,  -8.89,
     *  -6.60,  -7.13,  -4.77,  -5.33,  -3.06,
     *  -3.79,  -1.72,  -2.79,  -0.93,  -2.19,
     *  -0.52,  -1.90,  -0.45,  -2.20,  -1.22,
     *  -3.07,  -2.42,  -4.37,  -3.94,  -6.08/
      DATA TT22/
     *  -4.49,  -4.50,  -3.14,  -2.93,  -1.04,
     *  -1.36,   0.69,   0.21,   2.11,   1.33,
     *   3.29,   2.46,   4.30,   3.32,   4.79,
     *   3.62,   4.97,   3.64,   4.63,   3.07,
     *   4.06,   2.49,   3.30,   1.46,   2.06,
     *   0.51,   0.74,  -1.18,  -1.26,  -3.54,
     *  -3.97,  -5.26,  -4.18,  -3.71,  -2.10,
     *  -1.70,  -0.08,  -0.18,   0.94,   0.27,
     *   1.13,   0.08,   0.91,  -0.31,   0.49,
     *  -0.78,   0.08,  -1.15,  -0.23,  -1.41/
      DATA TT23/
     *  -0.42,  -1.55,  -0.55,  -1.66,  -0.66,
     *  -1.73,  -0.75,  -1.74,  -0.78,  -1.69,
     *  -0.78,  -1.60,  -0.75,  -1.46,  -0.67,
     *  -1.26,  -0.51,  -1.04,  -0.53,  -1.84,
     *  -2.42,  -4.52,  -4.76,  -6.33,  -6.76,
     *  -7.81,  -5.80,  -5.37,  -3.63,  -3.35,
     *  -1.75,  -1.88,  -0.61,  -0.90,   0.09,
     *  -0.32,   0.55,  -0.13,   0.70,  -0.06,
     *   0.49,  -0.20,   0.40,  -0.22,   0.36,
     *  -0.09,   0.58,   0.12,   0.75,   0.15/
      DATA TT24/
     *   0.70,   0.17,   1.11,   0.89,   1.85,
     *   1.62,   2.54,   2.29,   3.20,   2.91,
     *   3.84,   3.53,   4.48,   4.15,   5.12,
     *   4.78,   5.75,   5.39,   6.31,   5.91,
     *   6.87,   6.33,   7.13,   6.61,   7.30,
     *   6.31,   6.27,   4.83,   4.49,   2.85,
     *   2.32,   0.58,  -0.11,  -0.98,   0.81,
     *   1.77,   3.37,   4.13,   5.60,   6.15,
     *   7.29,   7.35,   7.95,   7.67,   8.16,
     *   7.83,   8.31,   8.01,   8.53,   8.27/
       END
      FUNCTION BARIE1(A,Z)
C  B[[[C[EH[E [AP[EPA [E[EH[[ [[PA A,Z (M[B)
C  ( BARASHENKOV,ILJINOV,TONEEV -1972)
      COMMON /BL1001/T1Y(130) /BL1002/T2XY(200)
      I=INT(A+0.1)
      J=INT(Z+0.1)
      L=I-J
      X=Z*Z/A
      IF(X.LE.33.5) BF0=12.5+4.7*(33.5-X)**0.75
      IF(X.GT.33.5) BF0=12.5-2.7*(X-33.5)**0.666667
      IF((2*(J/2)).LT.J) D=0.
      IF((2*(J/2)).EQ.J) D=-0.5
      IF((2*(L/2)).EQ.L) D=D
      IF((2*(L/2)).LT.L) D=D+1.
      IF(J.LE.0.OR.J.GT.130.OR.L.LE.0.OR.L.GT.200) GO TO 22
      T1=T1Y(J)
      T2=T2XY(L)
      GO TO 23
 22   T1=0.
      T2=0.
 23   CONTINUE
      BARIE1=BF0+D-T1-T2
      IF(BARIE1.LE.0.) BARIE1=0.
      RETURN
      END
      SUBROUTINE LYMASM(IZ,IA,CMASS,CBARR,NOBARR)
C
C     WILLIAM D. MYERS - 6 JULY 1970
C
      COMMON /BLOCX/ IPARQ
      COMMON /FFSS/  ENEX
      COMMON /FUSR/   BARR,SMASS,SHLL
C     COMMON /BLOCXX/ ARQ,F1Q,F1MQ,SUFNUC
C     COMMON /CCC/    WOTNUC,VOLNUC,COULMB,
C    *                A,Z,UN,A1,A2,A3,GGMMA,A3RT2,A3RT,ZSQ,
C    *                ODDEV,SYM,PARMAS,ACOR
      DIMENSION EM(10),EMP(10),XK(10),Y(2),F(2)
      DATA ZVT/1.6666666666/, ZT/.3333333333/,
     1 ZTT/.6666666667/, SR5/2.2360679775/
      DATA EM/0.00,2.00,8.00,14.00,28.00,50.00,82.00,126.00,
     1        184.00,258.00/
      DATA CAY1/0.0/, CAY2/0.0/, CAY3/2.0/, CAY4/11.0/, CAY5/8.07144/,
     * CAY6/7.28899/,
     *            D/.444/, C/5.8/, SMALC/.325/
C
C     DMASS = REMAINDER AFTER CM - NO SHELL EFFECTS SUBTRACTED
C     SHLL = CALCULATED SHELL EFFECT
C     DIFMAS= DMASS - SHLL
C
C------------------------------
      A1 = 15.4941
      IPARQ=0
C..... IPARQ=0  PARAMETRS MYERS-SWIATECKI
      IF(IPARQ.NE.0) GO TO 121
      A2 = 17.9439
      A3 = 0.7053
      GAMMA = 1.7826
      GO TO 126
 121  CONTINUE
C..... IPARQ=1  PARAMETRS KRAPPE-NIX
      IF(IPARQ.NE.1) GO TO 123
      A2 = 24.70
      A3 = 0.74476032
      GAMMA = 4.0
      GO TO 126
 123  CONTINUE
C..... IPARQ=2  PARAMETRS PAULI-LEDERGERBER
      IF(IPARQ.NE.2) GO TO 124
      A2 = 19.008
      A3 = 0.720
      GAMMA = 2.840
      GO TO 126
 124  CONTINUE
C..... IPARQ=3  BF(T)   PARAMETRS MYERS-SWIATECKI
      ALEVEL=0.1
      AMPAR=ALEVEL*FLOAT(IA)
      TSQ=ENEX/AMPAR
      A2 = 17.9439*(1.-0.0063157*TSQ)
      A3 = 0.7053*(1.-0.001*TSQ)
      GAMMA = 1.7826
 126  CONTINUE
C------------------------------
      IF(IZ.NE.0) GO TO 15
      CMASS=0.0
      RETURN
   15 NOBARR=0
      DO 1 I=1,10
      EMP(I)=EM(I)**ZVT
    1 CONTINUE
      DO 2 I=1,9
      XK(I)=.600*(EMP(I+1)-EMP(I)) /(EM(I+1)-EM(I))
    2 CONTINUE
C
C     FOR DEFINITIONS OF CAY1 AND RZ,SEE UCRL-11980
C
      CAY1=3.28637900*A3**3
      RZ=.86398700/A3
    5 Z= FLOAT(IZ)
      ZSQ=Z**2
      N=IA-IZ
      UN= FLOAT(N)
      A= FLOAT(IA)
      A3RT=A**ZT
      A3RT2=A3RT**2
      A2RT= SQRT(A)
      SYM=((UN-Z)/A)**2
      ACOR=1.00-GAMMA*SYM
      PARMAS=CAY5*UN+CAY6*Z
      VOLNUC=-A1*ACOR*A
      SUFNUC=A2*ACOR*A3RT2
      COULMB=A3*ZSQ/A3RT
      FUZSUR=-CAY1*ZSQ/A
      ODDEV=-(1.00+2.00*(N/2)-UN+2.00*(IZ/2)-Z)/A2RT*CAY4
      WTERM=-CAY2*A3RT2* EXP(-CAY3*SYM)
      WOTNUC=PARMAS+COULMB+FUZSUR+ODDEV+WTERM
      SMASS=WOTNUC+VOLNUC+SUFNUC
      SPW=SUFNUC+WTERM
      C2=SPW/A3RT2
      X=.5*COULMB/SPW
      IF(X.GE.1.00) GO TO 4
C------------------------------
      ARQ=X
      IF(IPARQ.EQ.0) BARR=SUFNUC*XI(X)
      IF(IPARQ.EQ.2) BARR=SUFNUC*XI(X)
      IF(IPARQ.EQ.3) BARR=SUFNUC*XI(X)
      IF(IPARQ.EQ.1) BARR=SUFNUC*XIMOD(X)
      IF(IPARQ.EQ.0) F1Q=XI(X)
      IF(IPARQ.EQ.1) F1MQ=XIMOD(X)
C------------------------------
      GO TO 6
    4 BARR=0.0
    6 Y(1)=UN
      Y(2)=Z
      DO 31 J=1,2
      DO 32 I=1,9
      IF(Y(J)-EM(I+1)) 3,3,32
   32 CONTINUE
      !WRITE 100 332,J
  332 FORMAT('1FAILURE IN LYMASS - Y(',I1,') EXCEEDS LAST MAGIC NO.')
      STOP
    3 F(J)=XK(I)*(Y(J)-EM(I))-.600*(Y(J)**ZVT-EMP(I))
   31 CONTINUE
      S=(2.00/A)**ZTT*(F(1)+F(2))-SMALC*A3RT
      C2D2=C2*D**2
      EE=(C2D2+C2D2)*(1.00-X)
      FF=.425917710*C2D2*D*(1.00+X+X)/A3RT
      SSHELL=C*S
      V=SSHELL/EE
      EPS=1.500*FF/EE
      IF(EE*(1.00-3.00*V).LE.0.00) GO TO 51
      QCALC=0.00
      THETA=0.00
      SHLL=SSHELL
      GO TO 52
C
C       ESTIMATE THETA
C
   51 TO=1.00
C
C       ITERATE TO FIND EQUILIBRIUM THETA
C
  101 DO 725 IPQ=1,10
      TO2=TO**2
C----------------------------------------
C     IF (TO2.GT.170.) !WRITE 100 500, IZ,IA
C 500 FORMAT(1X,'LYMASM',2X,2I5)
C----------------------------------------
C     EXMT2= EXP(-TO2)
      EXMT2= 1.E-20
      IF((ABS(TO2)).LT.30.) EXMT2= EXP(-TO2)
C
      T=TO-(1.00-EPS*TO-V*(3.00-TO2-TO2)*EXMT2) /
     1(-EPS+V*TO*(10.00-4.00*TO2)*EXMT2)
      IF(T.LE.0.00) GO TO 728
      IF( ABS(T-TO) .LT.1.E-4) GO TO 732
      TO=T
  725 CONTINUE
      GO TO 729
  732 T2=T**2
C     EXT2= EXP(-T2)
      EXT2= 1.E-20
      IF((ABS(T2)).LT.30.) EXT2= EXP(-T2)
C
      TEST=EE*(1.00-EPS*(T+T)-V*((4.00*T2-12.00)*T2+3.00)* EXT2)
      IF(TEST.GT.0.00) GO TO 81
  728 TO=.100
      DO 100 I=1,20
      TO2=TO**2
      GL=EE*(1.00-EPS*TO-V*(3.00-TO2-TO2)* EXP(-TO2))
      IF(GL.GT.0.00) GO TO 101
  100 CONTINUE
  729 CMASS=SMASS
      CBARR=0.00
      NOBARR=1
      RETURN
   81 THETA=T
      ALPHA0=D*SR5/A3RT
      ALPHA=ALPHA0*THETA
      SIGMA=ALPHA*(1.00+ALPHA/14.00)
      EXS= EXP(SIGMA+SIGMA)- EXP(-SIGMA)
      QCALC=4.E-3*Z*(RZ*A3RT)**2*EXS
      T2=T**2
      SHLL=T2*(EE-FF*T) + SSHELL*(1.00-T2-T2)* EXP(-T2)
   52 CMASS=SMASS+SHLL
      CBARR=BARR-SHLL
      RETURN
      END
      FUNCTION XI(Z)
C
C     6-POINT LAGRANGE INTERPOLATION
C
      DIMENSION Y(51)
      DATA Y/.25900,.255200,.250700,.245100,.2400,.23400,.228500,
     1      .22200,.21600,.2100,.20300,.196800,.1900,.18300,.175800,
     2      .1692400,.1620300,.1547800,.147500,.1401900,.1328400,
     3      .1254500,.1180100,.1105200,.1029600,.0953500,.0876800,
     4      .0799900,.0722900,.064600,.0569500,.0493700,.0419300,
     5      .0347600,.0281100,.0223600,.0176200,.0137300,.0105600,
     6      .0079800,.0059100,.0042500,.0029600,.0019700,.0012300,
     7      7.1E-4,3.6E-4,1.5E-4,4.E-5,1.E-5,0.00/
C
C     THE X VALUES ARE EVENLY SPACED - X = 0(.02)1
C
      ZBH=Z*50.00
      M=IFIX(ZBH)
      DEL=ZBH- FLOAT(M)
      M=M+1
      IF(M.LE.51) GO TO 105
      M=51
  100 XI=Y(M)
      RETURN
  105 IF (DEL.LT.1.E-4) GO TO 100
      IF(M.GE.3) GO TO 110
      DEL=DEL- FLOAT(3-M)
      M=3
      GO TO 115
  110 IF(M.LE.48) GO TO 115
      DEL=DEL+ FLOAT(M-48)
      M=48
  115 DM3=DEL-3.00
      PROD=DM3*DEL
      W6=1.00/(1.2E2*DM3)
      DM2=DM3+1.00
      PROD=DM2*PROD
       W5=-1.00/(24.00*DM2)
      DM1=DM2+1.00
      PROD=DM1*PROD
      W4=1.00/(12.00*DM1)
      DP1=DM1+2.00
      PROD=DP1*PROD
      W2=1.00/(24.00*DP1)
      DP2=DP1+1.00
      PROD=DP2*PROD
      W1=-1.00/(1.2E2*DP2)
      W3=-1.00/(12.00*DEL)
      XI=PROD*(W1*Y(M-2)+W2*Y(M-1)+W3*Y(M)+W4*Y(M+1)+W5*Y(M+2)
     1 +W6*Y(M+3))
      RETURN
      END
      FUNCTION XIMOD(Z)
C
C     6-POINT LAGRANGE INTERPOLATION
C      IN MODIFIED LIQUID-DROP FORMULA
C         ( KRAPPE [ NIX --  IAEA-SM-174/12 )
C
      DIMENSION Y(51)
      DATA  Y/
     1    0.12200, 0.12100, 0.11980, 0.11830, 0.11690, 0.11520, 0.1133,
     2    0.11130, 0.10900, 0.10670, 0.10420, 0.10150, 0.09850, 0.09540,
     3    0.09180, 0.08780, 0.08350, 0.07900, 0.07460, 0.06960, 0.06470,
     4    0.05960, 0.05420, 0.04880, 0.04350, 0.03880, 0.03400, 0.02920,
     5    0.02460, 0.02020, 0.01580, 0.01220, 0.00900, 0.00660, 0.00490,
     6    0.00360, 0.00280, 0.00220, 0.00180, 0.00140, 0.00100, 0.00090,
     7    0.00060, 0.00040, 0.00020, 0.00010, 0.00000, 0.00000, 0.00000,
     8    0.00000, 0.00000/
C
C     THE X VALUES ARE EVENLY SPACED - X = 0(.02)1
C
      ZBH=Z*50.00
      M=IFIX(ZBH)
      DEL=ZBH- FLOAT(M)
      M=M+1
      IF(M.LE.51) GO TO 105
      M=51
  100 XIMOD=Y(M)
      RETURN
  105 IF (DEL.LT.1.E-4) GO TO 100
      IF(M.GE.3) GO TO 110
      DEL=DEL- FLOAT(3-M)
      M=3
      GO TO 115
  110 IF(M.LE.48) GO TO 115
      DEL=DEL+ FLOAT(M-48)
      M=48
  115 DM3=DEL-3.00
      PROD=DM3*DEL
      W6=1.00/(1.2E2*DM3)
      DM2=DM3+1.00
      PROD=DM2*PROD
       W5=-1.00/(24.00*DM2)
      DM1=DM2+1.00
      PROD=DM1*PROD
      W4=1.00/(12.00*DM1)
      DP1=DM1+2.00
      PROD=DP1*PROD
      W2=1.00/(24.00*DP1)
      DP2=DP1+1.00
      PROD=DP2*PROD
      W1=-1.00/(1.2E2*DP2)
      W3=-1.00/(12.00*DEL)
      XIMOD=PROD*(W1*Y(M-2)+W2*Y(M-1)+W3*Y(M)+W4*Y(M+1)+W5*Y(M+2)
     1 +W6*Y(M+3))
      RETURN
      END

c **** random number generator for kiae-computer *****

C***************** GGUBFS     *******************************
C   IMSL ROUTINE NAME   - GGUBFS
C
C-----------------------------------------------------------------------
C
C   COMPUTER            - DG7/SINGLE
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - BASIC UNIFORM (0,1) RANDOM NUMBER GENERATOR
C                           FUNCTION FORM OF GGUBS
C
C   USAGE               - FUNCTION GGUBFS()
C
C   ARGUMENTS    GGUBFS - RESULTANT DEVIATE.
C                DSEED  - INPUT/OUTPUT DOUBLE PRECISION VARIABLE
C                           ASSIGNED AN INTEGER VALUE IN THE
C                           EXCLUSIVE RANGE (1.D0, 2147483647.D0).
C                           DSEED IS REPLACED BY A NEW VALUE TO BE
C                           USED IN A SUBSEQUENT CALL.
C
C   PRECISION/HARDWARE  - SINGLE/ALL
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
C
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.
C
C-----------------------------------------------------------------------
C
      REAL FUNCTION GGUBFS()
C                                  SPECIFICATIONS FOR ARGUMENTS
      DOUBLE PRECISION   DSEED
      COMMON /EE/   DSEED,ISEED,RANCL
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION   D2P31M,D2P31
C                                  D2P31M=(2**31) - 1
C                                  D2P31 =(2**31)(OR AN ADJUSTED VALUE)
      DATA               D2P31M/2147483647.D0/
      DATA               D2P31 /2147483711.D0/
      SAVE D2P31M,D2P31
C                                  FIRST EXECUTABLE STATEMENT
      DSEED = DMOD(16807.D0*DSEED,D2P31M)
      GGUBFS = DSEED / D2P31
      RANCL = RANCL + 1.
      RETURN
      END

      FUNCTION RNDM(IX)
      DOUBLE PRECISION   DSEED
      COMMON /EE/   DSEED,ISEED,RANCL
      RNDM=GGUBFS()
      RETURN
      END

      SUBROUTINE RNDM1(I)
      DOUBLE PRECISION   DSEED
      COMMON /EE/   DSEED,ISEED,RANCL
      ISEED=I
      IF(ISEED.EQ.0) ISEED=711
      DSEED=DBLE(ISEED)
      RETURN
      END

      SUBROUTINE DISNE0(INET,IP0,AMP,PN0,T,TN)
c *** implem.of flow.and rotat. EFLOW, EMOM in GeV/N. ANMOM in h-bar ***
C HAXO[[EH[E [HEP[[[ [ [M[[[[COB INET HE[TPOHOB [O [X TEM[EPAT[PE T(MEV)
C [O[HO[ K[H.[HEP[[[ TN(GEV),C[[TA[ PAC[PE[E[EH[E MAKCBE[OBCK[M. AMP(IP0
C -[X MACC[ (B GEV);[M[[[[C[ (GEV/C),B C[MME [A[[[E 0,[O[AB[[[TC[ B MACC
C PN0(3,500) C (IP0+1).
      COMMON /BLFLOW/IFLOW,A00,EFLOW /BLFLCM/ALF0,FLOWC,XCM,YCM,ZCM
      COMMON /BLAMOM/IAMOM,ANMOM,EMOM /BLAMCM/BMOM,WR,EMOMC
      COMMON /BLPAT/PA(500),PAZ(500),IPT
      COMMON /BLCEN/XC(500),YC(500),ZC(500)
      COMMON /CELIPS/CEL,IPOSF1
      COMMON /BLFLZ/FLZ
      DIMENSION AMP(500),PN0(3,500),ANL(3),PX(500),PY(500),PZ(500)
      DIMENSION PI(3),PJ(3),PN(3),AR(3),BR(3),VV(3)
      IF(INET.LE.0) RETURN
      IP01=IP0+1
      IP0M=IP0+INET
      IF(INET-1) 1,1,2
 1    P=SQRT(2.*AMP(IP01)*TN)
      CALL ISOTR(ANL)
      VV(1)=0.
      VV(2)=0.
      VV(3)=0.
      IF(IFLOW.GE.1) VV(1)=VV(1)+ALF0*(XC(IP01)-XCM)
      IF(IFLOW.GE.1) VV(2)=VV(2)+ALF0*(YC(IP01)-YCM)
c      IF(IFLOW.GE.1) VV(3)=VV(3)+ALF0*(ZC(IP01)-ZCM)
c assuming that at CEL>1 the flow everywhere on the border is the same
c      IF(IFLOW.GE.1) VV(3)=VV(3)+ALF0*(ZC(IP01)/CEL-ZCM)
c -- for nonuniform flow: along Z-axis
      IF(IFLOW.GE.1) VV(3)=VV(3)+FLZ*ALF0*(ZC(IP01)-ZCM)
c ------------------------------------
      IF(IAMOM.GE.1) VV(3)=VV(3)-WR*(YC(IP01)-YCM)
      IF(IAMOM.GE.1) VV(2)=VV(2)+WR*(ZC(IP01)-ZCM)
      PN0(1,IP01)=P*ANL(1)+AMP(IP01)*VV(1)
      PN0(2,IP01)=P*ANL(2)+AMP(IP01)*VV(2)
      PN0(3,IP01)=P*ANL(3)+AMP(IP01)*VV(3)
      RETURN
 2    IF(INET-2) 3,3,4
 3    P=SQRT(2.*(AMP(IP01)*AMP(IP0M)/(AMP(IP01)+AMP(IP0M)))*TN)
      CALL ISOTR(ANL)
      VV(1)=0.
      VV(2)=0.
      VV(3)=0.
      IF(IFLOW.GE.1) VV(1)=VV(1)+ALF0*(XC(IP01)-XCM)
      IF(IFLOW.GE.1) VV(2)=VV(2)+ALF0*(YC(IP01)-YCM)
c      IF(IFLOW.GE.1) VV(3)=VV(3)+ALF0*(ZC(IP01)-ZCM)
c assuming that at CEL>1 the flow everywhere on the border is the same
c      IF(IFLOW.GE.1) VV(3)=VV(3)+ALF0*(ZC(IP01)/CEL-ZCM)
c -- for nonuniform flow: along Z-axis
      IF(IFLOW.GE.1) VV(3)=VV(3)+FLZ*ALF0*(ZC(IP01)-ZCM)
c ------------------------------------
      IF(IAMOM.GE.1) VV(3)=VV(3)-WR*(YC(IP01)-YCM)
      IF(IAMOM.GE.1) VV(2)=VV(2)+WR*(ZC(IP01)-ZCM)
      PN0(1,IP01)=P*ANL(1)+AMP(IP01)*VV(1)
      PN0(2,IP01)=P*ANL(2)+AMP(IP01)*VV(2)
      PN0(3,IP01)=P*ANL(3)+AMP(IP01)*VV(3)
      VV(1)=0.
      VV(2)=0.
      VV(3)=0.
      IF(IFLOW.GE.1) VV(1)=VV(1)+ALF0*(XC(IP0M)-XCM)
      IF(IFLOW.GE.1) VV(2)=VV(2)+ALF0*(YC(IP0M)-YCM)
c      IF(IFLOW.GE.1) VV(3)=VV(3)+ALF0*(ZC(IP0M)-ZCM)
c assuming that at CEL>1 the flow everywhere on the border is the same
c      IF(IFLOW.GE.1) VV(3)=VV(3)+ALF0*(ZC(IP0M)/CEL-ZCM)
c -- for nonuniform flow: along Z-axis
      IF(IFLOW.GE.1) VV(3)=VV(3)+FLZ*ALF0*(ZC(IP0M)-ZCM)
c ------------------------------------
      IF(IAMOM.GE.1) VV(3)=VV(3)-WR*(YC(IP0M)-YCM)
      IF(IAMOM.GE.1) VV(2)=VV(2)+WR*(ZC(IP0M)-ZCM)
      PN0(1,IP0M)=-P*ANL(1)+AMP(IP0M)*VV(1)
      PN0(2,IP0M)=-P*ANL(2)+AMP(IP0M)*VV(2)
      PN0(3,IP0M)=-P*ANL(3)+AMP(IP0M)*VV(3)
      RETURN
 4    IP0M2=IP0M-2
      ES=0.
      PSX=0.
      PSY=0.
      PSZ=0.
      FEMT=SQRT(0.5*T)*EXP(-0.5)
      DO 5 I=IP01,IP0M2
 7    E=RNDM(-1)*9.*T
      FE=SQRT(E)*EXP(-E/T)
      FERAND=RNDM(-1)*FEMT
      IF(FERAND-FE) 6,6,7
 6    P=SQRT(2.*E*0.001*AMP(I))
      CALL ISOTR(ANL)
      PX(I)=P*ANL(1)
      PY(I)=P*ANL(2)
      PZ(I)=P*ANL(3)
      ES=ES+E*0.001
      PSX=PSX+PX(I)
      PSY=PSY+PY(I)
      PSZ=PSZ+PZ(I)
 5    CONTINUE
      I1=IP0M2+1
      I2=IP0M2+2
      PPX=-PSX
      PPY=-PSY
      PPZ=-PSZ
      P=SQRT(PPX*PPX+PPY*PPY+PPZ*PPZ)
      EE=TN-ES
      EM=P*P/(2.*(AMP(I1)+AMP(I2)))
      IF(EE.LE.EM) GO TO 4
      H=1.+AMP(I2)/AMP(I1)
      CTM12=H*(1.-2.*AMP(I2)*EE/(P*P))
 11   CT1=1.-2.*RNDM(-1)
      IF(CT1*CT1-CTM12) 11,11,12
 12   IF(CTM12) 13,17,17
 13   IZN=1
      GO TO 15
 17   IF(CT1) 11,14,14
 14   IF(RNDM(-1)-0.5) 16,16,13
 16   IZN=-1
 15   CONTINUE
      P1=(P*CT1+IZN*SQRT(P*P*CT1*CT1-P*P*CTM12))/H
      P2=SQRT(P1*P1+P*P-2.*P1*P*CT1)
      PHI=6.28318*RNDM(-1)
      ST1=SQRT(1.-CT1*CT1)
      CPHI1=COS(PHI)
      SPHI1=SIN(PHI)
      CPHI2=-CPHI1
      SPHI2=-SPHI1
      CT2=(P*P+P2*P2-P1*P1)/(2.*P*P2)
      IF(CT2.GT.-1.AND.CT2.LT.1) GO TO 20
      ST2=0.
      GO TO 21
 20   ST2=SQRT(1.-CT2*CT2)
 21   CONTINUE
      PI(1)=P1*ST1*CPHI1
      PI(2)=P1*ST1*SPHI1
      PI(3)=P1*CT1
      PJ(1)=P2*ST2*CPHI2
      PJ(2)=P2*ST2*SPHI2
      PJ(3)=P2*CT2
      AR(1)=PPX
      AR(2)=PPY
      AR(3)=PPZ
      BR(1)=1.
      BR(2)=0.
      BR(3)=0.
      CALL ROTOR(AR,BR,PI,PN)
      PX(I1)=PN(1)
      PY(I1)=PN(2)
      PZ(I1)=PN(3)
      CALL ROTOR(AR,BR,PJ,PN)
      PX(I2)=PN(1)
      PY(I2)=PN(2)
      PZ(I2)=PN(3)
      PSX=PSX+PX(I1)+PX(I2)
      PSY=PSY+PY(I1)+PY(I2)
      PSZ=PSZ+PZ(I1)+PZ(I2)
      ES=ES+(PX(I1)**2+PY(I1)**2+PZ(I1)**2)/(2.*AMP(I1))
     *+(PX(I2)**2+PY(I2)**2+PZ(I2)**2)/(2.*AMP(I2))
      DO 8 I=IP01,IP0M
      VV(1)=0.
      VV(2)=0.
      VV(3)=0.
      IF(IFLOW.GE.1) VV(1)=VV(1)+ALF0*(XC(I)-XCM)
      IF(IFLOW.GE.1) VV(2)=VV(2)+ALF0*(YC(I)-YCM)
c      IF(IFLOW.GE.1) VV(3)=VV(3)+ALF0*(ZC(I)-ZCM)
c assuming that at CEL>1 the flow everywhere on the border is the same
c      IF(IFLOW.GE.1) VV(3)=VV(3)+ALF0*(ZC(I)/CEL-ZCM)
c -- for nonuniform flow: along Z-axis
      IF(IFLOW.GE.1) VV(3)=VV(3)+FLZ*ALF0*(ZC(I)-ZCM)
c ------------------------------------
      IF(IAMOM.GE.1) VV(3)=VV(3)-WR*(YC(I)-YCM)
      IF(IAMOM.GE.1) VV(2)=VV(2)+WR*(ZC(I)-ZCM)
      PN0(1,I)=PX(I)+AMP(I)*VV(1)
      PN0(2,I)=PY(I)+AMP(I)*VV(2)
      PN0(3,I)=PZ(I)+AMP(I)*VV(3)
 8    CONTINUE
      RETURN
      END
      SUBROUTINE CULIM0(IP0,A0,ZC0,AM,CH,PN0,TEMP,ECOLB)
      DIMENSION X0(500),Y0(500),Z0(500),VX0(500),VY0(500),VZ0(500),
     *VX(500),VY(500),VZ(500),AM(500),CH(500),VN(500),RS(500),
     *PN0(3,500),PNT(3,500)
      COMMON /BLPAT/PA(500),PAZ(500),IPT
      COMMON /BLCEN/XC(500),YC(500),ZC(500)  /BLDT/DT,IT
      COMMON /BLFLOW/IFLOW,A00,EFLOW /BLFLCM/ALF0,FLOWC,XCM,YCM,ZCM
      COMMON /BLAMOM/IAMOM,ANMOM,EMOM /BLAMCM/BMOM,WR,EMOMC
      DT=2.
      IP=IP0
      DO 1 I=1,IP
      X0(I)=XC(I)
      Y0(I)=YC(I)
      Z0(I)=ZC(I)
 1    CONTINUE
      ET=1.5*IP*TEMP*0.001
c      CALL DISNE0(IP,0,AM,PNT,TEMP,ET)
      CALL DISN02(IP,0,AM,PNT,TEMP,ET)
c --- statistic on energy of hot fragments before coulomb ---
      CALL BFREN(IP,0,AM,PNT)
c -----------------------------------------------------------
      DO 40 I=1,IP
      VX0(I)=PNT(1,I)/AM(I)
      VY0(I)=PNT(2,I)/AM(I)
      VZ0(I)=PNT(3,I)/AM(I)
 40   CONTINUE
      ETT=ET
      IF(IFLOW.GE.1) ETT=ETT+FLOWC
      IF(IAMOM.GE.1) ETT=ETT+EMOMC
      CALL CULON(X0,Y0,Z0,VX0,VY0,VZ0,VX,VY,VZ,AM,CH,VN,RS,IP,ECOLB,
     *ETT)
      DO 7 I=1,IP
      PN0(1,I)=AM(I)*VX(I)
      PN0(2,I)=AM(I)*VY(I)
      PN0(3,I)=AM(I)*VZ(I)
 7    CONTINUE
      RETURN
      END
      SUBROUTINE POS0(IPT,IP0,A0,ZC0,AM,CH)
c *** implem.of flow and rotat. EFLOW, EMOM in GeV/N. ANMOM in h-bar ***
      DIMENSION AM(500),CH(500)
      COMMON /BLPAT/PA(500),PAZ(500),IP
      COMMON /BLCEN/XC(500),YC(500),ZC(500)
      COMMON /BLFLOW/IFLOW,A00,EFLOW /BLFLCM/ALF0,FLOWC,XCM,YCM,ZCM
      COMMON /BLAMOM/IAMOM,ANMOM,EMOM /BLAMCM/BMOM,WR,EMOMC
      COMMON /CELIPS/CEL,IPOSF1
      COMMON /BLFLZ/FLZ
      A00=A0
      RN=1.17
      RSYS=1.82*RN*(A0**0.333333)
c for possible extension at CEL=2
c      RSYS=1.44*RN*(A0**0.333333)
c for celx=0.5 in place
c      RSYS=2.*1.44*RN*(A0**0.333333)
c --------------------
      IP=IPT
      DO 50 I=1,IP
      PA(I)=AM(I)/0.94
 50   PAZ(I)=CH(I)
c      !WRITE 100 101,IP,(PA(I),I=1,10)
c 101  FORMAT(1x,'ip=',I3,' pa=',10F7.1)
c      CALL PLACM(RSYS,RN)
      CALL PLACE(RSYS,RN)
c --------------------------------------------
cc "projectile" max is placed in the right, "target" max - in the left along
cc the Z-axis: mass asymetry (new version instead PLACM - mirrow reflection)
c      JXX=0
c      XTP0=119./(119.+129.)
c      XXX=RNDM(-1)
c      IF(XXX.GT.XTP0.AND.ZC(1).LT.0.) JXX=1
c      IF(XXX.LE.XTP0.AND.ZC(1).GE.0.) JXX=1
c      IF(JXX.EQ.0) GO TO 7
c      DO 6 IXX=1,IP
c      ZC(IXX)=-ZC(IXX)
c 6    CONTINUE
c 7    CONTINUE
c --------------------------------------------
      XCM=0.
      YCM=0.
      ZCM=0.
      SPA=0.
      DO 1 I=1,IP
      XCM=XCM+PA(I)*XC(I)
      YCM=YCM+PA(I)*YC(I)
      ZCM=ZCM+PA(I)*ZC(I)
      SPA=SPA+PA(I)
 1    CONTINUE
      XCM=XCM/SPA
      YCM=YCM/SPA
      ZCM=ZCM/SPA
      SMR2=0.
      BMOM=0.
      DO 2 I=1,IP
c      RAD2=(XC(I)-XCM)**2+(YC(I)-YCM)**2+(ZC(I)-ZCM)**2
c assuming that at CEL>1 the flow everywhere on the border is the same
c      RAD2=(XC(I)-XCM)**2+(YC(I)-YCM)**2+(ZC(I)/CEL-ZCM)**2
c -- for nonuniform flow: along Z-axis
      RAD2=(XC(I)-XCM)**2+(YC(I)-YCM)**2+(FLZ*FLZ)*(ZC(I)-ZCM)**2
c ------------------------------------
      RZY2=(ZC(I)-ZCM)**2+(YC(I)-YCM)**2
      SMR2=SMR2+AM(I)*RAD2
      BMOM=BMOM+AM(I)*RZY2
 2    CONTINUE
      ALF0=SQRT(2.*EFLOW*A00/SMR2)
      WR=SQRT(2.*EMOM*A00/BMOM)
c      WR=(ANMOM/5.06)/BMOM
c      EMOM=0.5*BMOM*WR*WR
      ANMOM=5.06*BMOM*WR
      EMOMC=0.
      FLOWC=0.
      DO 3 I=1,IP0
c      RAD2=(XC(I)-XCM)**2+(YC(I)-YCM)**2+(ZC(I)-ZCM)**2
c assuming that at CEL>1 the flow everywhere on the border is the same
c      RAD2=(XC(I)-XCM)**2+(YC(I)-YCM)**2+(ZC(I)/CEL-ZCM)**2
c -- for nonuniform flow: along Z-axis
      RAD2=(XC(I)-XCM)**2+(YC(I)-YCM)**2+(FLZ*FLZ)*(ZC(I)-ZCM)**2
c ------------------------------------
      RZY2=(ZC(I)-ZCM)**2+(YC(I)-YCM)**2
      FLOWC=FLOWC+0.5*AM(I)*RAD2*ALF0*ALF0
      EMOMC=EMOMC+0.5*AM(I)*RZY2*WR*WR
3     CONTINUE
      RETURN
      END

c ***
      SUBROUTINE THELI0
      COMMON /BHELI1/Y6LI(46),Y7LI(46),Y3HE(46),Y4HE(46),
     *Y6LI0,Y7LI0,Y3HE0,Y4HE0,NUMB
      COMMON /BHELI2/A(6),B(6),IA1(6),IZ1(6),IA2(6),IZ2(6)
      COMMON /BHELI3/Y1(6),Y11(6),Y2(6),Y21(6)
      COMMON /BHELI4/Y1F(6),Y11F(6),Y2F(6),Y21F(6)
      COMMON /BHELI5/Y1E(6),Y11E(6),Y2E(6),Y21E(6)
      COMMON /BHELI6/Y1M(6),Y11M(6),Y2M(6),Y21M(6)
      NUMB=0
      Y6LI0=0.
      Y7LI0=0.
      Y3HE0=0.
      Y4HE0=0.
      DO 121 I=1,46
      Y6LI(I)=0.
      Y7LI(I)=0.
      Y3HE(I)=0.
      Y4HE(I)=0.
 121  CONTINUE
c --- ratio: R12=(Y(A1,Z1)/Y(A1+1,Z1))/(Y(A2,Z2)/Y(A2+1,Z2))
c --- assumption: R12=(1/a)*exp(B/T) , a - spin and mass factor,
c --- B=BE(A1,Z1)-BE(A1+1,Z1)-BE(A2,Z2)+BE(A2+1,Z2) (Huang,Xi,Lynch,Tsang,...)
c   6,7-Li/3,4-He   : a=2.18 , B=13.32 (MeV)
c   9,10-Be/3,4-He  : a=0.38 , B=13.76
c   2,3-H/3,4-He    : a=1.59 , B=14.29
c   11,12-B/3,4-He  : a=1.11 , B=17.20  corrected
c   10,11-B/3,4-He  : a=0.86 , B= 9.12  corrected
c   7,8-Li/3,4-He   : a=1.98 , B=18.54
      DO 122 I=1,6
      Y1(I)=0.
      Y11(I)=0.
      Y2(I)=0.
      Y21(I)=0.
      Y1F(I)=0.
      Y11F(I)=0.
      Y2F(I)=0.
      Y21F(I)=0.
      Y1E(I)=0.
      Y11E(I)=0.
      Y2E(I)=0.
      Y21E(I)=0.
      Y1M(I)=0.
      Y11M(I)=0.
      Y2M(I)=0.
      Y21M(I)=0.
 122  CONTINUE
      RETURN
      END
      BLOCK DATA B5
      COMMON /BHELI2/A(6),B(6),IA1(6),IZ1(6),IA2(6),IZ2(6)
      DATA A/ 2.18, 0.38, 1.59, 1.11, 0.86, 1.98/
      DATA B/13.32,13.76,14.29,17.20, 9.12,18.54/
      DATA IA1/ 6, 9, 2,11,10, 7/
      DATA IA2/ 3, 3, 3, 3, 3, 3/
      DATA IZ1/ 3, 4, 1, 5, 5, 3/
      DATA IZ2/ 2, 2, 2, 2, 2, 2/
      END
      SUBROUTINE THELIP
      COMMON /BLFINL/A0F,Z0F,E0F,TF,FL,SRF,WW3F
      COMMON /BHELI1/Y6LI(46),Y7LI(46),Y3HE(46),Y4HE(46),
     *Y6LI0,Y7LI0,Y3HE0,Y4HE0,NUMB
      COMMON /BHELI3/Y1(6),Y11(6),Y2(6),Y21(6)
      COMMON /BHELI4/Y1F(6),Y11F(6),Y2F(6),Y21F(6)
      COMMON /BHELI5/Y1E(6),Y11E(6),Y2E(6),Y21E(6)
      COMMON /BHELI6/Y1M(6),Y11M(6),Y2M(6),Y21M(6)
      DIMENSION TX(40)
      !WRITE 100 13,A0F,Z0F,E0F,TF
 13   FORMAT(2X,' A0F,Z0F=',2F6.1,'   E0F, TF =',2F9.3)
      !WRITE 100 1,NUMB
  1   FORMAT(1X,'in THELI :',
     *'  numb.of anal.ev. NUMB=',I7)
c      !WRITE 100 8
c 8    FORMAT(1X,'yield of 3-He in excit.energy (MeV/N) ',
c     *' (from 0 to 20, step-0.5)   Y3HE:')
c      !WRITE 100 5,(Y3HE(I),I=1,40)
c      !WRITE 100 9,(Y3HE(I),I=41,46)
c      !WRITE 100 10
c 10   FORMAT(1X,'yield of 4-He in excit.energy (MeV/N) ',
c     *' (from 0 to 20, step-0.5)   Y4HE:')
c      !WRITE 100 5,(Y4HE(I),I=1,40)
c      !WRITE 100 9,(Y4HE(I),I=41,46)
c      !WRITE 100 11
c 11   FORMAT(1X,'yield of 6-Li in excit.energy (MeV/N) ',
c     *' (from 0 to 20, step-0.5)   Y6LI:')
c      !WRITE 100 5,(Y6LI(I),I=1,40)
c      !WRITE 100 9,(Y6LI(I),I=41,46)
c      !WRITE 100 12
c 12   FORMAT(1X,'yield of 7-Li in excit.energy (MeV/N) ',
c     *' (from 0 to 20, step-0.5)   Y7LI:')
c      !WRITE 100 5,(Y7LI(I),I=1,40)
c      !WRITE 100 9,(Y7LI(I),I=41,46)
      DO 2 I=1,40
      TX(I)=0.
      XX=0.
      IF(Y6LI(I).GT.0.AND.Y7LI(I).GT.0.AND.Y3HE(I).GT.0.AND.
     *Y4HE(I).GT.0) XX=(Y6LI(I)*Y4HE(I))/(Y7LI(I)*Y3HE(I))
      IF(XX.GT.0.) TX(I)=16./(ALOG(2.18*XX))
 2    CONTINUE
c      !WRITE 100 3
c 3    FORMAT(2X,'HE-LI Temperature(MeV) in excit.energy(MeV/N).'
c     *' (from 0 to 20, step-0.5)  TX:')
c      !WRITE 100 4,(TX(I),I=1,40)
 4    FORMAT(1X,10F8.2)
 5    FORMAT(1X,10F8.0)
 9    FORMAT(1X,6F13.3)
      !WRITE 100 6,Y3HE0,Y4HE0,Y6LI0,Y7LI0
 6    FORMAT(2X,'total yields: Y3HE0=',F7.0,'  Y4HE0=',F7.0,
     *'  Y6LI0=',F6.0,'  Y7LI0=',F6.0)
      TEMHL0=0.
      XX0=0.
      IF(Y6LI0.GT.0.AND.Y7LI0.GT.0.AND.Y3HE0.GT.0.AND.
     *Y4HE0.GT.0) XX0=(Y6LI0*Y4HE0)/(Y7LI0*Y3HE0)
      IF(XX0.GT.0.) TEMHL0=16./(ALOG(2.18*XX0))
      !WRITE 100 7,TEMHL0
 7    FORMAT(1X,'avr. He-Li temperat.(*1.2): TEMHL0=',F7.2,
     *3X,'summary yields & temperatures:')
      CALL TISOT(Y1,Y11,Y2,Y21)
      !WRITE 100 30
 30   FORMAT(1X,'yields & temperatures at freezout (multifragm.):')
      CALL TISOT(Y1M,Y11M,Y2M,Y21M)
      !WRITE 100 31
 31   FORMAT(1X,'yields & temperatures after fermi-break-up :')
      CALL TISOT(Y1F,Y11F,Y2F,Y21F)
      !WRITE 100 32
 32   FORMAT(1X,'yields & temperatures after evaporation :')
      CALL TISOT(Y1E,Y11E,Y2E,Y21E)
      RETURN
      END
      SUBROUTINE TISOT(Y1,Y11,Y2,Y21)
      COMMON /BHELI2/A(6),B(6),IA1(6),IZ1(6),IA2(6),IZ2(6)
      DIMENSION Y1(6),Y11(6),Y2(6),Y21(6),RISO(6),TISO(6)
 15   FORMAT(11X,6(F7.0,1X))
      !WRITE 100 14
 14   FORMAT(1X,'yields of ','  6-Li',3X,' 9-Be',3X,'  2-H',3X,
     *' 11-B',3X,' 10-B',3X,' 7-Li')
      !WRITE 100 15,(Y1(L),L=1,6)
      !WRITE 100 16
 16   FORMAT(1X,'yields of ','  7-Li',3X,'10-Be',3X,'  3-H',3X,
     *' 12-B',3X,' 11-B',3X,' 8-Li')
      !WRITE 100 15,(Y11(L),L=1,6)
      !WRITE 100 17
 17   FORMAT(1X,'yields of ','  3-He',3X,' 3-He',3X,' 3-He',3X,
     *' 3-He',3X,' 3-He',3X,' 3-He')
      !WRITE 100 15,(Y2(L),L=1,6)
      !WRITE 100 18
 18   FORMAT(1X,'yields of ','  4-He',3X,' 4-He',3X,' 4-He',3X,
     *' 4-He',3X,' 4-He',3X,' 4-He')
      !WRITE 100 15,(Y21(L),L=1,6)
      DO 19 L=1,6
      TISO(L)=0.
      RISO(L)=0.
      IF(Y1(L).EQ.0.OR.Y11(L).EQ.0.OR.Y2(L).EQ.0.OR.Y21(L).EQ.0)
     *GO TO 19
      RISO(L)=(Y1(L)/Y11(L))/(Y2(L)/Y21(L))
      TISO(L)=B(L)/ALOG(A(L)*RISO(L))
 19   CONTINUE
      !WRITE 100 20,(RISO(L),L=1,6)
 20   FORMAT(1X,'isotop.ratio RISO(1-6)=',6F8.3)
      !WRITE 100 21,(TISO(L),L=1,6)
 21   FORMAT(1X,'isot.temper. TISO(MeV)=',6F8.3)
      RETURN
      END
      SUBROUTINE THELI1(KSTART)
      COMMON /BLOKC/SPT(10,500) /BLFINL/A0F,Z0F,E0F,TF,FL,SRF,WW3F
      COMMON /BENTRO/SSR
      COMMON /BHELI1/Y6LI(46),Y7LI(46),Y3HE(46),Y4HE(46),
     *Y6LI0,Y7LI0,Y3HE0,Y4HE0,NUMB
      COMMON /BHELI2/A(6),B(6),IA1(6),IZ1(6),IA2(6),IZ2(6)
      COMMON /BHELI3/Y1(6),Y11(6),Y2(6),Y21(6)
      COMMON /BHELI4/Y1F(6),Y11F(6),Y2F(6),Y21F(6)
      COMMON /BHELI5/Y1E(6),Y11E(6),Y2E(6),Y21E(6)
      COMMON /BHELI6/Y1M(6),Y11M(6),Y2M(6),Y21M(6)
      KST1=KSTART-1
      IF(KST1.LT.1) GO TO 410
      IF(A0F.LT.0.) GO TO 410
c      !WRITE 100 8,A0F,E0F
c 8    FORMAT(1X,' in THELI,  BREAK UP.  A0F, E0F =',F6.1,F9.3)
      UA=E0F
      DO 105 J=1,KST1
c      !WRITE 100 7,(SPT(JPER,J),JPER=1,10)
c 7    FORMAT(1X,'SPT=',10F8.3)
      IZ=INT(SPT(8,J)+0.5)
      IA=INT(SPT(9,J)/0.94+0.5)
      IF(IA.EQ.4.AND.IZ.EQ.2) Y4HE0=Y4HE0+1.
      IF(IA.EQ.3.AND.IZ.EQ.2) Y3HE0=Y3HE0+1.
      IF(IA.EQ.6.AND.IZ.EQ.3) Y6LI0=Y6LI0+1.
      IF(IA.EQ.7.AND.IZ.EQ.3) Y7LI0=Y7LI0+1.
      IF(IA.EQ.4.AND.IZ.EQ.2)
     *CALL HIST1(UA,0.,20.,0.5,Y4HE,46,1.)
      IF(IA.EQ.3.AND.IZ.EQ.2)
     *CALL HIST1(UA,0.,20.,0.5,Y3HE,46,1.)
      IF(IA.EQ.6.AND.IZ.EQ.3)
     *CALL HIST1(UA,0.,20.,0.5,Y6LI,46,1.)
      IF(IA.EQ.7.AND.IZ.EQ.3)
     *CALL HIST1(UA,0.,20.,0.5,Y7LI,46,1.)
 105  CONTINUE
      NUMB=NUMB+1
      DO 106 J=1,KST1
      IZ=INT(SPT(8,J)+0.5)
      IA=INT(SPT(9,J)/0.94+0.5)
      ILAB=INT(SPT(2,J)+0.1)
      DO 107 L=1,6
      IF(IA.EQ.IA1(L).AND.IZ.EQ.IZ1(L)) Y1(L)=Y1(L)+1.
      IF(IA.EQ.(IA1(L)+1).AND.IZ.EQ.IZ1(L)) Y11(L)=Y11(L)+1.
      IF(IA.EQ.IA2(L).AND.IZ.EQ.IZ2(L)) Y2(L)=Y2(L)+1.
      IF(IA.EQ.(IA2(L)+1).AND.IZ.EQ.IZ2(L)) Y21(L)=Y21(L)+1.
      IF(ILAB) 108,108,109
 108  IF(IA.EQ.IA1(L).AND.IZ.EQ.IZ1(L)) Y1M(L)=Y1M(L)+1.
      IF(IA.EQ.(IA1(L)+1).AND.IZ.EQ.IZ1(L)) Y11M(L)=Y11M(L)+1.
      IF(IA.EQ.IA2(L).AND.IZ.EQ.IZ2(L)) Y2M(L)=Y2M(L)+1.
      IF(IA.EQ.(IA2(L)+1).AND.IZ.EQ.IZ2(L)) Y21M(L)=Y21M(L)+1.
      GO TO 107
 109  IF(ILAB-1) 110,110,111
 110  IF(IA.EQ.IA1(L).AND.IZ.EQ.IZ1(L)) Y1F(L)=Y1F(L)+1.
      IF(IA.EQ.(IA1(L)+1).AND.IZ.EQ.IZ1(L)) Y11F(L)=Y11F(L)+1.
      IF(IA.EQ.IA2(L).AND.IZ.EQ.IZ2(L)) Y2F(L)=Y2F(L)+1.
      IF(IA.EQ.(IA2(L)+1).AND.IZ.EQ.IZ2(L)) Y21F(L)=Y21F(L)+1.
      GO TO 107
 111  IF(ILAB-2) 112,112,107
 112  IF(IA.EQ.IA1(L).AND.IZ.EQ.IZ1(L)) Y1E(L)=Y1E(L)+1.
      IF(IA.EQ.(IA1(L)+1).AND.IZ.EQ.IZ1(L)) Y11E(L)=Y11E(L)+1.
      IF(IA.EQ.IA2(L).AND.IZ.EQ.IZ2(L)) Y2E(L)=Y2E(L)+1.
      IF(IA.EQ.(IA2(L)+1).AND.IZ.EQ.IZ2(L)) Y21E(L)=Y21E(L)+1.
 107  CONTINUE
 106  CONTINUE
 410  CONTINUE
      RETURN
      END
      SUBROUTINE RELCOR(IP,AMP,ZMP,PN0)
c recalculation of nonrelativistic momenta into relativistic
c ones using kinetic energy of fragments found in nonrelat.way
      COMMON /BLERR/ERRMI,ERRMA,ERRFR,ERRRE
      DIMENSION AMP(500),ZMP(500),PN0(3,500),P(3)
      TKS=0.
      CMS=0.
      TK2=0.
      P11=0.
      P12=0.
      P13=0.
      P21=0.
      P22=0.
      P23=0.
      V10=0.
      V20=0.
      DO 1 I=1,IP
      P(1)=PN0(1,I)
      P(2)=PN0(2,I)
      P(3)=PN0(3,I)
      P11=P11+P(1)
      P12=P12+P(2)
      P13=P13+P(3)
      CM=AMP(I)
      CMS=CMS+CM
      PQ=P(1)**2+P(2)**2+P(3)**2
      V10=V10+SQRT(PQ)/CM
      T=PQ/(2.*CM)
      TKS=TKS+T
 1    CONTINUE
c *** calc. coef. ALF for new momenta ***
      E0S=CMS+TKS
      HALF=0.1
      ALF=1.
      KK=0
      K1=0
      K2=0
      ID=0
 24   CONTINUE
      E0T=0.
      DO 10 IT=1,IP
      P00=PN0(1,IT)**2+PN0(2,IT)**2+PN0(3,IT)**2
      E0T=E0T+SQRT(AMP(IT)**2+(ALF**2)*P00)
 10   CONTINUE
      KK=KK+1
      D=(E0S-E0T)/E0S
c      !WRITE 100 60,KK,K1,K2,ID,D,IP,ALF,E0S,E0T
c   60 FORMAT(2X,' *** KK,K1,K2=',3I3,' ID=',I2,' D=',E12.5,' IP=',
c     *I3,' ALF=',F5.2,' E0S,E0T=',2F8.3)
      IF(ABS(D).LT.0.001) GO TO 25
      IF(KK.GT.60.OR.ID.GT.30) GO TO 61
      H=SIGN(HALF,D)
      IF(D) 21,21,22
 21   K1=1
      IF(K1.EQ.1.AND.K2.EQ.1) ID=ID+1
      ALF=ALF+H/2**ID
      GO TO 24
 22   K2=1
      IF(K1.EQ.1.AND.K2.EQ.1) ID=ID+1
      ALF=ALF+H/2**ID
      GO TO 24
   61 ERRRE=ERRRE+1
c      !WRITE 100 62,KK,ID,D,IP,ALF
c   62 FORMAT(5X,'ERROR IN RELCOR  KK=',I2,' ID=',I2,' D=',E12.5,' IP=',
c     *I3,' ALF=',F8.3)
      ALF=1.
 25   CONTINUE
c *************************
      DO 4 I=1,IP
      PN0(1,I)=ALF*PN0(1,I)
      PN0(2,I)=ALF*PN0(2,I)
      PN0(3,I)=ALF*PN0(3,I)
      P(1)=PN0(1,I)
      P(2)=PN0(2,I)
      P(3)=PN0(3,I)
      CM=AMP(I)
      P21=P21+P(1)
      P22=P22+P(2)
      P23=P23+P(3)
      PQ=P(1)**2+P(2)**2+P(3)**2
      TK2=TK2+PQ/(2.*CM)
      CM0=SQRT(PQ+CM**2)
      V20=V20+SQRT(PQ)/CM0
 4    CONTINUE
c      !WRITE 100 2, IP,ALF,CMS,TKS,TK2,P11,P12,P13,P21,P22,P23
c 2    FORMAT(2X,'IP=',I3,' ALF=',F5.2,' CMS=',F8.3,
c     *' relat.transf.: TKS,TK2=',2F8.3,
c     *' P11,P12,P13=',3F8.3,' P21,P22,P23=',3F8.3)
      V10=V10/FLOAT(IP)
      V20=V20/FLOAT(IP)
c      !WRITE 100 3, V10,V20
c 3    FORMAT(2X,' nonrelat.: V10=',F8.3,' relat.: V20=',F8.3)
      RETURN
      END
      SUBROUTINE DMR00
c ... Y0(1-Z0) - yield of fragm. at M_imf=0
c ... YF(1-10,1-Z0) - yield of fragm. at M_imf=1,2,3,4,5,6,7,8,10
c ... R23(1-Z0) = YF(2,1-Z0)/YF(3,1-Z0)
      COMMON /BDMR1/Y0(100),YF(10,100),YFS(10),MULT(20)
      COMMON /BDMR2/R12(100),R23(100),R34(100),R45(100)
     *,R56(100),R67(100),R78(100),R89(100)
      COMMON /BDMR3/Y0T(100),YFT(10,100),YFST(10)
      COMMON /BDMR4/R12T(100),R23T(100),R34T(100),R45T(100)
     *,R56T(100),R67T(100),R78T(100),R89T(100)
      COMMON /BDMR5/ZMAX,NUMB,ZMA(10),ZIMF0,NIMF0,ZIMF(10),NIMF(10)
     *,ZMAI(10),ZMA2(10),EZ6(10),YZ6(10)
      ZMAX=0.
      ZIMF0=0.
      NUMB=0
      NIMF0=0
      DO 1 I=1,100
      IF(I.LE.20) MULT(I)=0
      IF(I.LE.10) YFS(I)=0.
      IF(I.LE.10) YFST(I)=0.
      IF(I.LE.10) ZMA(I)=0.
      IF(I.LE.10) NIMF(I)=0
      IF(I.LE.10) ZIMF(I)=0.
      IF(I.LE.10) ZMAI(I)=0.
      IF(I.LE.10) ZMA2(I)=0.
      IF(I.LE.10) EZ6(I)=0.
      IF(I.LE.10) YZ6(I)=0.
      Y0(I)=0.
      Y0T(I)=0.
      R12(I)=0.
      R23(I)=0.
      R34(I)=0.
      R45(I)=0.
      R56(I)=0.
      R67(I)=0.
      R78(I)=0.
      R89(I)=0.
      R12T(I)=0.
      R23T(I)=0.
      R34T(I)=0.
      R45T(I)=0.
      R56T(I)=0.
      R67T(I)=0.
      R78T(I)=0.
      R89T(I)=0.
      DO 1 J=1,10
      YF(J,I)=0.
      YFT(J,I)=0.
 1    CONTINUE
      RETURN
      END
      SUBROUTINE DMR
      COMMON /BDMR0/IAA(100),IZZ(100),NFRAG1
     *,PXX(100),PYY(100),PZZ(100)
      COMMON /BDMR1/Y0(100),YF(10,100),YFS(10),MULT(20)
      COMMON /BDMR2/R12(100),R23(100),R34(100),R45(100)
     *,R56(100),R67(100),R78(100),R89(100)
      COMMON /BDMR3/Y0T(100),YFT(10,100),YFST(10)
      COMMON /BDMR4/R12T(100),R23T(100),R34T(100),R45T(100)
     *,R56T(100),R67T(100),R78T(100),R89T(100)
      COMMON /BDMR5/ZMAX,NUMB,ZMA(10),ZIMF0,NIMF0,ZIMF(10),NIMF(10)
     *,ZMAI(10),ZMA2(10),EZ6(10),YZ6(10)
      MIMF=0
      ZX=0.
      JZX=0
      ZF=0.
      EK6=0.
      NZ6=0
      DO 1 IR=1,NFRAG1
      IA=IAA(IR)
      IZ=IZZ(IR)
      IF(IZ.GT.ZX) JZX=IR
      IF(IZ.GT.ZX) ZX=IZ
      IF(IZ.LT.3.OR.IZ.GT.20) GO TO 1
      IF(IZ.GE.3.AND.IZ.LE.20) MIMF=MIMF+1
      IF(IZ.GE.3.AND.IZ.LE.20) ZF=ZF+IZ
      IF(IZ.EQ.6) P00=PXX(IR)**2+PYY(IR)**2+PZZ(IR)**2
      IF(IZ.EQ.6) NZ6=NZ6+1
      IF(IZ.EQ.6) EK6=EK6+P00/(2.*938.*IA)
 1    CONTINUE
      ZX2=0.
      ZFX=0.
      DO 5 IR=1,NFRAG1
      IA=IAA(IR)
      IZ=IZZ(IR)
      IF(IR.EQ.JZX) GO TO 7
      IF(IZ.GT.ZX2) ZX2=IZ
 7    CONTINUE
      IF(IZ.LT.3.OR.IZ.GT.20) GO TO 5
      IF(IZ.GT.ZFX) ZFX=IZ
 5    CONTINUE
      ZMAX=ZMAX+ZX
      NUMB=NUMB+1
      IF(MIMF.LE.0) GO TO 6
      ZF=ZF/FLOAT(MIMF)
      ZIMF0=ZIMF0+ZF
      NIMF0=NIMF0+1
      IF(MIMF.GT.10) GO TO 6
      ZIMF(MIMF)=ZIMF(MIMF)+ZF
      ZMAI(MIMF)=ZMAI(MIMF)+ZFX
      NIMF(MIMF)=NIMF(MIMF)+1
      ZMA(MIMF)=ZMA(MIMF)+ZX
      ZMA2(MIMF)=ZMA2(MIMF)+ZX2
      EZ6(MIMF)=EZ6(MIMF)+EK6
      YZ6(MIMF)=YZ6(MIMF)+NZ6
 6    CONTINUE
      DO 2 IR=1,NFRAG1
      IA=IAA(IR)
      IZ=IZZ(IR)
      IF(IZ.LE.0.OR.IZ.GT.100) GO TO 4
      IF(MIMF.EQ.0) Y0T(IZ)=Y0T(IZ)+1.
      IF(MIMF.LE.0.OR.MIMF.GT.10) GO TO 4
      YFT(MIMF,IZ)=YFT(MIMF,IZ)+1.
      YFST(MIMF)=YFST(MIMF)+1.
 4    CONTINUE
      IF(IZ.LT.3.OR.IZ.GT.20) GO TO 3
      IF(MIMF.EQ.0) Y0(IZ)=Y0(IZ)+1.
      IF(MIMF.LE.0.OR.MIMF.GT.10) GO TO 3
      YF(MIMF,IZ)=YF(MIMF,IZ)+1.
      YFS(MIMF)=YFS(MIMF)+1.
 3    CONTINUE
 2    CONTINUE
      IF(MIMF.GE.0.AND.MIMF.LE.19)
     *MULT(MIMF+1)=MULT(MIMF+1)+1
      RETURN
      END
      SUBROUTINE DMRPR
      COMMON /BDMR1/Y0(100),YF(10,100),YFS(10),MULT(20)
      COMMON /BDMR2/R12(100),R23(100),R34(100),R45(100)
     *,R56(100),R67(100),R78(100),R89(100)
      COMMON /BDMR3/Y0T(100),YFT(10,100),YFST(10)
      COMMON /BDMR4/R12T(100),R23T(100),R34T(100),R45T(100)
     *,R56T(100),R67T(100),R78T(100),R89T(100)
      COMMON /BDMR5/ZMAX,NUMB,ZMA(10),ZIMF0,NIMF0,ZIMF(10),NIMF(10)
     *,ZMAI(10),ZMA2(10),EZ6(10),YZ6(10)
 2    FORMAT(2X,10F7.0)
 3    FORMAT(2X,10F7.4)
      !WRITE 100 5,(MULT(I),I=1,20)
 5    FORMAT(2X,'imf MULT(0-19)=',10I6)
      !WRITE 100 6,(YFS(I),I=1,10)
 6    FORMAT(2X,'YFS(1-10)=',10F7.0)
      !WRITE 100 116,(YFST(I),I=1,10)
 116  FORMAT(2X,'YFST(1-10)=',10F7.0)
      DO 7 J=1,10
      DO 8 II=1,100
      IF(YFS(J).GT.0) YF(J,II)=YF(J,II)/YFS(J)
      IF(YFST(J).GT.0) YFT(J,II)=YFT(J,II)/YFST(J)
 8    CONTINUE
 7    CONTINUE
      IZMT=79
      IZM=20
      !WRITE 100 10
 10   FORMAT(5X,'(at M_imf= 0)  charge yield:  Y0(1-ZM)  ')
      !WRITE 100 2,(Y0(I),I=1,IZM)
      DO 11 J=1,10
      !WRITE 100 12
 12   FORMAT(2X,'-------------------------------------- ')
      !WRITE 100 13,J,J
 13   FORMAT(5X,'(at M_imf=',I2,')  chagre yiled:  YF(',I2,
     *',1-ZM)  ')
      !WRITE 100 3,(YF(J,I),I=1,IZM)
 11   CONTINUE
      !WRITE 100 12
      !WRITE 100 100
 100  FORMAT(5X,'(at M_imf= 0)  charge yield:  Y0T(1-ZMT)  ')
      !WRITE 100 2,(Y0T(I),I=1,IZMT)
      DO 101 J=1,10
      !WRITE 100 12
      !WRITE 100 103,J,J
 103  FORMAT(5X,'(at M_imf=',I2,')  chagre yiled:  YFT(',I2,
     *',1-ZMT)  ')
      !WRITE 100 3,(YFT(J,I),I=1,IZMT)
 101  CONTINUE
      !WRITE 100 12
      DO 1 III=1,100
      IF(YF(2,III).GT.0.) R12(III)=YF(1,III)/YF(2,III)
      IF(YF(3,III).GT.0.) R23(III)=YF(2,III)/YF(3,III)
      IF(YF(4,III).GT.0.) R34(III)=YF(3,III)/YF(4,III)
      IF(YF(5,III).GT.0.) R45(III)=YF(4,III)/YF(5,III)
      IF(YF(6,III).GT.0.) R56(III)=YF(5,III)/YF(6,III)
      IF(YF(7,III).GT.0.) R67(III)=YF(6,III)/YF(7,III)
      IF(YF(8,III).GT.0.) R78(III)=YF(7,III)/YF(8,III)
      IF(YF(9,III).GT.0.) R89(III)=YF(8,III)/YF(9,III)
 1    CONTINUE
      DO 104 III=1,100
      IF(YFT(2,III).GT.0.) R12T(III)=YFT(1,III)/YFT(2,III)
      IF(YFT(3,III).GT.0.) R23T(III)=YFT(2,III)/YFT(3,III)
      IF(YFT(4,III).GT.0.) R34T(III)=YFT(3,III)/YFT(4,III)
      IF(YFT(5,III).GT.0.) R45T(III)=YFT(4,III)/YFT(5,III)
      IF(YFT(6,III).GT.0.) R56T(III)=YFT(5,III)/YFT(6,III)
      IF(YFT(7,III).GT.0.) R67T(III)=YFT(6,III)/YFT(7,III)
      IF(YFT(8,III).GT.0.) R78T(III)=YFT(7,III)/YFT(8,III)
      IF(YFT(9,III).GT.0.) R89T(III)=YFT(8,III)/YFT(9,III)
 104  CONTINUE
      !WRITE 100 15
 15   FORMAT(5X,'R12(1-ZM)=YF(1,1-ZM)/YF(2,1-ZM) : ')
      !WRITE 100 3,(R12(I),I=1,IZM)
      !WRITE 100 12
      !WRITE 100 4
 4    FORMAT(5X,'R23(1-ZM)=YF(2,1-ZM)/YF(3,1-ZM) : ')
      !WRITE 100 3,(R23(I),I=1,IZM)
      !WRITE 100 12
      !WRITE 100 9
 9    FORMAT(5X,'R34(1-ZM)=YF(3,1-ZM)/YF(4,1-ZM) : ')
      !WRITE 100 3,(R34(I),I=1,IZM)
      !WRITE 100 12
      !WRITE 100 14
 14   FORMAT(5X,'R45(1-ZM)=YF(4,1-ZM)/YF(5,1-ZM) : ')
      !WRITE 100 3,(R45(I),I=1,IZM)
      !WRITE 100 12
      !WRITE 100 16
 16   FORMAT(5X,'R56(1-ZM)=YF(5,1-ZM)/YF(6,1-ZM) : ')
      !WRITE 100 3,(R56(I),I=1,IZM)
      !WRITE 100 12
      !WRITE 100 17
 17   FORMAT(5X,'R67(1-ZM)=YF(6,1-ZM)/YF(7,1-ZM) : ')
      !WRITE 100 3,(R67(I),I=1,IZM)
      !WRITE 100 12
      !WRITE 100 18
 18   FORMAT(5X,'R78(1-ZM)=YF(7,1-ZM)/YF(8,1-ZM) : ')
      !WRITE 100 3,(R78(I),I=1,IZM)
      !WRITE 100 12
      !WRITE 100 19
 19   FORMAT(5X,'R89(1-ZM)=YF(8,1-ZM)/YF(9,1-ZM) : ')
      !WRITE 100 3,(R89(I),I=1,IZM)
      !WRITE 100 12
      !WRITE 100 105
 105  FORMAT(5X,'R12T(1-ZMT) : ')
      !WRITE 100 3,(R12T(I),I=1,IZMT)
      !WRITE 100 12
      !WRITE 100 106
 106  FORMAT(5X,'R23T(1-ZMT) : ')
      !WRITE 100 3,(R23T(I),I=1,IZMT)
      !WRITE 100 12
      !WRITE 100 107
 107  FORMAT(5X,'R34T(1-ZMT) : ')
      !WRITE 100 3,(R34T(I),I=1,IZMT)
      !WRITE 100 12
      !WRITE 100 108
 108  FORMAT(5X,'R45T(1-ZMT) : ')
      !WRITE 100 3,(R45T(I),I=1,IZMT)
      !WRITE 100 12
      !WRITE 100 109
 109  FORMAT(5X,'R56T(1-ZMT) : ')
      !WRITE 100 3,(R56T(I),I=1,IZMT)
      !WRITE 100 12
      !WRITE 100 110
 110  FORMAT(5X,'R67T(1-ZMT) : ')
      !WRITE 100 3,(R67T(I),I=1,IZMT)
      !WRITE 100 12
      !WRITE 100 111
 111  FORMAT(5X,'R78T(1-ZMT) : ')
      !WRITE 100 3,(R78T(I),I=1,IZMT)
      !WRITE 100 12
      !WRITE 100 118
 118  FORMAT(5X,'R89T(1-ZMT) : ')
      !WRITE 100 3,(R89T(I),I=1,IZMT)
      IF(NIMF0.GT.0) ZIMF0=ZIMF0/FLOAT(NIMF0)
      IF(NUMB.GT.0) ZMAX=ZMAX/FLOAT(NUMB)
      !WRITE 100 114,ZMAX,ZIMF0
 114  FORMAT(2X,'max.charge ZMAX=',F6.2,
     *' avr. IMF charge ZIMF0=',F6.2)
      !WRITE 100 113,(NIMF(I),I=1,10)
 113  FORMAT(2X,'NIMF(1-10)=',10I7)
      DO 112 II=1,10
      IF(NIMF(II).GT.0) ZIMF(II)=ZIMF(II)/FLOAT(NIMF(II))
      IF(NIMF(II).GT.0) ZMA(II)=ZMA(II)/FLOAT(NIMF(II))
      IF(NIMF(II).GT.0) ZMA2(II)=ZMA2(II)/FLOAT(NIMF(II))
      IF(NIMF(II).GT.0) ZMAI(II)=ZMAI(II)/FLOAT(NIMF(II))
      IF(YZ6(II).GT.0) EZ6(II)=EZ6(II)/YZ6(II)
 112  CONTINUE
      !WRITE 100 117,(ZMA(II),II=1,10)
 117  FORMAT(2X,' ZMA(1-10)=',10F7.2)
      !WRITE 100 120,(ZMA2(II),II=1,10)
 120  FORMAT(2X,'ZMA2(1-10)=',10F7.2)
      !WRITE 100 115,(ZIMF(II),II=1,10)
 115  FORMAT(2X,'ZIMF(1-10)=',10F7.2)
      !WRITE 100 119,(ZMAI(II),II=1,10)
 119  FORMAT(2X,'ZMAI(1-10)=',10F7.2)
      !WRITE 100 131,(YZ6(II),II=1,10)
 131  FORMAT(2X,' YZ6(1-10)=',10F7.0)
      !WRITE 100 130,(EZ6(II),II=1,10)
 130  FORMAT(2X,' EZ6(1-10)=',10F7.2)
      RETURN
      END
      SUBROUTINE CPOCH1(NUMPAT)
      COMMON /BLOKC/SPT(10,500)
      ZMAX=0.
      ZMA2=0.
      ZMA3=0.
      JMAX=0
      JMA2=0
      JMA3=0
      NZ330=0
      DO 2 IT=1,NUMPAT
      IZ=INT(SPT(8,IT)+0.5)
      ZM=IZ
      IF(IZ.GE.3.AND.IZ.LE.30) NZ330=NZ330+1
      IF(ZM.GT.ZMAX) JMAX=IT
      IF(ZM.GT.ZMAX) ZMAX=ZM
 2    CONTINUE
      DO 3 IT=1,NUMPAT
      IF(IT.EQ.JMAX) GO TO 3
      ZM=INT(SPT(8,IT)+0.5)
      IF(ZM.GT.ZMA2) JMA2=IT
      IF(ZM.GT.ZMA2) ZMA2=ZM
 3    CONTINUE
      DO 4 IT=1,NUMPAT
      IF(IT.EQ.JMAX) GO TO 4
      IF(IT.EQ.JMA2) GO TO 4
      ZM=INT(SPT(8,IT)+0.5)
      IF(ZM.GT.ZMA3) JMA3=IT
      IF(ZM.GT.ZMA3) ZMA3=ZM
 4    CONTINUE
      CALL POCH1(NUMPAT,NZ330,ZMAX,ZMA2,ZMA3,JMAX,JMA2,JMA3)
      RETURN
      END
      SUBROUTINE BPOCH1
      COMMON /BLPOC1/IPOC0,IPOC(5) /BLPOC3/PDOM(26)
      COMMON /BLPOC2/EKTS(5),EKT2S(5),EKCS(5),EKC2S(5)
      IPOC0=0
      DO 1 I=1,5
      IPOC(I)=0
      EKTS(I)=0.
      EKT2S(I)=0.
      EKCS(I)=0.
      EKC2S(I)=0.
 1    CONTINUE
      DO 2 I=1,26
      PDOM(I)=0.
 2    CONTINUE
c --- for writing in Pochodzalla format
c      OPEN(60,NAME='FIL123')
      RETURN
      END
      SUBROUTINE EPOCH1
      COMMON /BLPOC1/IPOC0,IPOC(5) /BLPOC3/PDOM(26)
      COMMON /BLPOC2/EKTS(5),EKT2S(5),EKCS(5),EKC2S(5)
      DIMENSION SIGET(5),SIGEC(5)
      !WRITE 100 1,IPOC0
 1    FORMAT(2X,'number of events with N_imf.GE.3 and Z_ma3.GE.8 -',
     *' IPOC0=',I5)
      !WRITE 100 2,(IPOC(I),I=1,5)
 2    FORMAT(2X,'IPOC(1-5)=',5I8)
      DO 3 I=1,5
      SIGET(I)=0.
      IF(IPOC(I).GT.0) EKTS(I)=EKTS(I)/FLOAT(IPOC(I))
      IF(IPOC(I).GT.0) EKT2S(I)=EKT2S(I)/FLOAT(IPOC(I))
      TT=EKT2S(I)-EKTS(I)*EKTS(I)
      IF(TT.GT.0.) SIGET(I)=SQRT(TT)
      SIGEC(I)=0.
      IF(IPOC(I).GT.0) EKCS(I)=EKCS(I)/FLOAT(IPOC(I))
      IF(IPOC(I).GT.0) EKC2S(I)=EKC2S(I)/FLOAT(IPOC(I))
      TC=EKC2S(I)-EKCS(I)*EKCS(I)
      IF(TC.GT.0.) SIGEC(I)=SQRT(TC)
 3    CONTINUE
      !WRITE 100 4,(EKTS(I),I=1,5)
 4    FORMAT(2X,'EKTS(1-5)=',5F8.3)
      !WRITE 100 5,(SIGET(I),I=1,5)
 5    FORMAT(1X,'SIGET(1-5)=',5F8.3)
      !WRITE 100 6,(EKCS(I),I=1,5)
 6    FORMAT(2X,'EKCS(1-5)=',5F8.3)
      !WRITE 100 7,(SIGEC(I),I=1,5)
 7    FORMAT(1X,'SIGEC(1-5)=',5F8.3)
      !WRITE 100 9
 9    FORMAT(2X,'distribution in DOM (0-1, step-0.05), PDOM(1-26):')
      !WRITE 100 8,(PDOM(I),I=1,26)
 8    FORMAT(1X,10F8.3)
      RETURN
      END
      SUBROUTINE POCH1(KST1,NZ330,ZMAX,ZMA2,ZMA3,JMAX,JMA2,JMA3)
      COMMON /BLOKC/SPT(10,500)
      COMMON /BLPOC1/IPOC0,IPOC(5) /BLPOC3/PDOM(26)
      COMMON /BLPOC2/EKTS(5),EKT2S(5),EKCS(5),EKC2S(5)
      DIMENSION P1(3),P2(3),P3(3),P1C(3),P2C(3),P3C(3),V0(3)
      DIMENSION PSO(3),PSN(3),V1C(3),V2C(3),V3C(3)
c ---- for checking momentum balance and relat. momenta ----
      DO 50 III=1,3
      PSO(III)=0.
      PSN(III)=0.
      V1C(III)=0.
      V2C(III)=0.
      V3C(III)=0.
 50   CONTINUE
c ----
      IF(NZ330.GE.3) GO TO 1
      RETURN
 1    IF(ZMA3.GE.7.5) GO TO 2
      RETURN
 2    IPOC0=IPOC0+1
      AMAX=SPT(9,JMAX)/0.94
      AMA2=SPT(9,JMA2)/0.94
      AMA3=SPT(9,JMA3)/0.94
      EC=ZMAX*ZMA2/(AMAX**0.33333+AMA2**0.33333)+
     *   ZMAX*ZMA3/(AMAX**0.33333+AMA3**0.33333)+
     *   ZMA2*ZMA3/(AMA2**0.33333+AMA3**0.33333)
      EC=1.44*EC/1.4
      IC=0
      IF(EC.GE.30.AND.EC.LT.60.) IC=1
      IF(EC.GE.60.AND.EC.LT.90.) IC=2
      IF(EC.GE.90.AND.EC.LT.120.) IC=3
      IF(EC.GE.120.AND.EC.LT.150.) IC=4
      IF(EC.GE.150.AND.EC.LT.180.) IC=5
      IF(IC.EQ.0) RETURN
      IPOC(IC)=IPOC(IC)+1
      EKT=SPT(7,JMAX)+SPT(7,JMA2)+SPT(7,JMA3)
      EKT2=EKT*EKT
      EKTS(IC)=EKTS(IC)+EKT
      EKT2S(IC)=EKT2S(IC)+EKT2
      CM1=AMAX*0.94
      CM2=AMA2*0.94
      CM3=AMA3*0.94
      CT1=SPT(4,JMAX)
      ST1=SQRT(1.-CT1*CT1)
      CF1=SPT(6,JMAX)
      SF1=SPT(5,JMAX)
      CT2=SPT(4,JMA2)
      ST2=SQRT(1.-CT2*CT2)
      CF2=SPT(6,JMA2)
      SF2=SPT(5,JMA2)
      CT3=SPT(4,JMA3)
      ST3=SQRT(1.-CT3*CT3)
      CF3=SPT(6,JMA3)
      SF3=SPT(5,JMA3)
      CALL TINP(P1,CT1,ST1,CF1,SF1,SPT(7,JMAX),CM1)
      CALL TINP(P2,CT2,ST2,CF2,SF2,SPT(7,JMA2),CM2)
      CALL TINP(P3,CT3,ST3,CF3,SF3,SPT(7,JMA3),CM3)
      P01=P1(1)+P2(1)+P3(1)
      P02=P1(2)+P2(2)+P3(2)
      P03=P1(3)+P2(3)+P3(3)
c ----
      PR12O=SQRT((P1(1)-P2(1))**2+(P1(2)-P2(2))**2+
     *(P1(3)-P2(3))**2)
      PR23O=SQRT((P2(1)-P3(1))**2+(P2(2)-P3(2))**2+
     *(P2(3)-P3(3))**2)
      PR31O=SQRT((P3(1)-P1(1))**2+(P3(2)-P1(2))**2+
     *(P3(3)-P1(3))**2)
      PSO(1)=P01
      PSO(2)=P02
      PSO(3)=P03
c ----
      EM0=SQRT(P01*P01+P02*P02+P03*P03+(CM1+CM2+CM3)**2)
      V0(1)=-P01/EM0
      V0(2)=-P02/EM0
      V0(3)=-P03/EM0
      EM1=SQRT(P1(1)**2+P1(2)**2+P1(3)**2+CM1**2)
      EM2=SQRT(P2(1)**2+P2(2)**2+P2(3)**2+CM2**2)
      EM3=SQRT(P3(1)**2+P3(2)**2+P3(3)**2+CM3**2)
      CALL CLPV(P1,V0,P1C,EM1)
      CALL CLPV(P2,V0,P2C,EM2)
      CALL CLPV(P3,V0,P3C,EM3)
      EM1C=SQRT(P1C(1)**2+P1C(2)**2+P1C(3)**2+CM1**2)
      EM2C=SQRT(P2C(1)**2+P2C(2)**2+P2C(3)**2+CM2**2)
      EM3C=SQRT(P3C(1)**2+P3C(2)**2+P3C(3)**2+CM3**2)
      DO 57 JJ=1,3
      V1C(JJ)=P1C(JJ)/EM1C
      V2C(JJ)=P2C(JJ)/EM2C
      V3C(JJ)=P3C(JJ)/EM3C
 57   CONTINUE
      V12C=SQRT((V1C(1)-V2C(1))**2+(V1C(2)-V2C(2))**2+
     *(V1C(3)-V2C(3))**2)
      V23C=SQRT((V2C(1)-V3C(1))**2+(V2C(2)-V3C(2))**2+
     *(V2C(3)-V3C(3))**2)
      V31C=SQRT((V3C(1)-V1C(1))**2+(V3C(2)-V1C(2))**2+
     *(V3C(3)-V1C(3))**2)
      SSS=0.5*(V12C+V23C+V31C)
      VEQ=(V12C+V23C+V31C)/3.
      FEQ=(SQRT(3.)/4.)*VEQ*VEQ
      FSS=SQRT(SSS*(SSS-V12C)*(SSS-V23C)*(SSS-V31C))
      DOM=FSS/FEQ
      WW=1.
      CALL HIST1(DOM,0.,1.,0.05,PDOM,26,WW)
c ----
      PR12N=SQRT((P1C(1)-P2C(1))**2+(P1C(2)-P2C(2))**2+
     *(P1C(3)-P2C(3))**2)
      PR23N=SQRT((P2C(1)-P3C(1))**2+(P2C(2)-P3C(2))**2+
     *(P2C(3)-P3C(3))**2)
      PR31N=SQRT((P3C(1)-P1C(1))**2+(P3C(2)-P1C(2))**2+
     *(P3C(3)-P1C(3))**2)
      PSN(1)=P1C(1)+P2C(1)+P3C(1)
      PSN(2)=P1C(2)+P2C(2)+P3C(2)
      PSN(3)=P1C(3)+P2C(3)+P3C(3)
c      !WRITE 100 51,PSO(1),PSO(2),PSO(3),PSN(1),PSN(2),PSN(3)
c 51   FORMAT(2X,' PSO(1-3)=',3F9.4,'  PSN(1-3)=',3F9.4)
c      !WRITE 100 52,PR12O,PR23O,PR31O,PR12N,PR23N,PR31N
c 52   FORMAT(2X,' PR12O,PR23O,PR31O=',3F9.4,
c     *'  PR12N,PR23N,PR31N=',3F9.4)
c ----
      CALL PINT(P1C,CT1C,ST1C,CF1C,SF1C,T1C,CM1)
      CALL PINT(P2C,CT2C,ST2C,CF2C,SF2C,T2C,CM2)
      CALL PINT(P3C,CT3C,ST3C,CF3C,SF3C,T3C,CM3)
      EKC=T1C+T2C+T3C
      EKC2=EKC*EKC
      EKCS(IC)=EKCS(IC)+EKC
      EKC2S(IC)=EKC2S(IC)+EKC2
c --- writing for Pochodzalla evaluat.programm
      T1C=1000.*T1C
      T2C=1000.*T2C
      T3C=1000.*T3C
      EKC=1000.*EKC
c --- writing for Pochodzalla format
c      WRITE(60,*)(ZMAX,AMAX,T1C)
c      WRITE(60,*)(ZMA2,AMA2,T2C)
c      WRITE(60,*)(ZMA3,AMA3,T3C,EKC)
c --------------------------------------------
c      !WRITE 100 53,IC,EKT,EKC
c 53   FORMAT(2X,' IC=',I2,'  EKT=',F8.4,' EKC=',F8.4)
c      !WRITE 100 54,AMAX,ZMAX,SPT(7,JMAX),T1C
c 54   FORMAT(2X,' AMAX,ZMAX=',2F6.1,' T1O,T1C=',2F8.4)
c      !WRITE 100 55,AMA2,ZMA2,SPT(7,JMA2),T2C
c 55   FORMAT(2X,' AMA2,ZMA2=',2F6.1,' T2O,T2C=',2F8.4)
c      !WRITE 100 56,AMA3,ZMA3,SPT(7,JMA3),T3C
c 56   FORMAT(2X,' AMA3,ZMA3=',2F6.1,' T3O,T3C=',2F8.4)
      RETURN
      END
      SUBROUTINE BSPEC1
      COMMON /BLSPE1/EN(56),ENT(56),EP(56),EPT(56),EA(56),EAT(56)
      COMMON /BLSPE2/EDE(56),ETR(56),EHE3(56)
      COMMON /BLSPE3/SEN,SEP,SED,SET,SEH,SEA,NN,NP,ND,NT,NH,NA
      NN=0
      NP=0
      ND=0
      NT=0
      NH=0
      NA=0
      SEN=0.
      SEP=0.
      SED=0.
      SET=0.
      SEH=0.
      SEA=0.
      DO 1 I=1,56
      EN(I)=0.
      ENT(I)=0.
      EP(I)=0.
      EPT(I)=0.
      EA(I)=0.
      EAT(I)=0.
      EDE(I)=0.
      ETR(I)=0.
      EHE3(I)=0.
 1    CONTINUE
      RETURN
      END
      SUBROUTINE ESPEC1
      COMMON /BLSPE1/EN(56),ENT(56),EP(56),EPT(56),EA(56),EAT(56)
      COMMON /BLSPE2/EDE(56),ETR(56),EHE3(56)
      COMMON /BLSPE3/SEN,SEP,SED,SET,SEH,SEA,NN,NP,ND,NT,NH,NA
      IF(NN.GT.0) SEN=SEN/FLOAT(NN)
      IF(NP.GT.0) SEP=SEP/FLOAT(NP)
      IF(ND.GT.0) SED=SED/FLOAT(ND)
      IF(NT.GT.0) SET=SET/FLOAT(NT)
      IF(NH.GT.0) SEH=SEH/FLOAT(NH)
      IF(NA.GT.0) SEA=SEA/FLOAT(NA)
      !WRITE 100 13, NN,NP,ND,NT,NH,NA
 13   FORMAT(2X,'numb.of part. N,P,D,T,H,A =',6I9)
      !WRITE 100 14, SEN,SEP,SED,SET,SEH,SEA
 14   FORMAT(2X,'their avr.kin.energ.(MeV) =',6F9.2)
      TEMN=0.
      TEMP=0.
      TEMD=0.
      TEMT=0.
      TEMH=0.
      TEMA=0.
      TTT30=0.
      TTT60=0.
      TTT30=(EN(15)+EN(16)+EN(17))/3.
      TTT60=(EN(30)+EN(31)+EN(32))/3.
      IF(TTT60.GT.0.AND.TTT30.GT.0) TTT=ALOG(TTT30/TTT60)
      IF(TTT.GT.0.) TEMN=30./TTT
      TTT30=0.
      TTT60=0.
      TTT30=(EP(15)+EP(16)+EP(17))/3.
      TTT60=(EP(30)+EP(31)+EP(32))/3.
      IF(TTT60.GT.0.AND.TTT30.GT.0) TTT=ALOG(TTT30/TTT60)
      IF(TTT.GT.0.) TEMP=30./TTT
      TTT30=0.
      TTT60=0.
      TTT30=(EDE(15)+EDE(16)+EDE(17))/3.
      TTT60=(EDE(30)+EDE(31)+EDE(32))/3.
      IF(TTT60.GT.0.AND.TTT30.GT.0) TTT=ALOG(TTT30/TTT60)
      IF(TTT.GT.0.) TEMD=30./TTT
      TTT30=0.
      TTT60=0.
      TTT30=(ETR(15)+ETR(16)+ETR(17))/3.
      TTT60=(ETR(30)+ETR(31)+ETR(32))/3.
      IF(TTT60.GT.0.AND.TTT30.GT.0) TTT=ALOG(TTT30/TTT60)
      IF(TTT.GT.0.) TEMT=30./TTT
      TTT30=0.
      TTT60=0.
      TTT30=(EHE3(15)+EHE3(16)+EHE3(17))/3.
      TTT60=(EHE3(30)+EHE3(31)+EHE3(32))/3.
      IF(TTT60.GT.0.AND.TTT30.GT.0) TTT=ALOG(TTT30/TTT60)
      IF(TTT.GT.0.) TEMH=30./TTT
      TTT30=0.
      TTT60=0.
      TTT30=(EA(15)+EA(16)+EA(17))/3.
      TTT60=(EA(30)+EA(31)+EA(32))/3.
      IF(TTT60.GT.0.AND.TTT30.GT.0) TTT=ALOG(TTT30/TTT60)
      IF(TTT.GT.0.) TEMA=30./TTT
      !WRITE 100 12, TEMN,TEMP,TEMD,TEMT,TEMH,TEMA
 12   FORMAT(2X,'inver.slop.temper.(30/60) =',6F9.2)
 10   FORMAT(1X,10F8.0)
 11   FORMAT(1X,6F12.4)
      !WRITE 100 1
 1    FORMAT(2X,'energy spect. of neutrons, total: EN',
     *' (0-25MeV, step-0.5MeV)')
      !WRITE 100 10,(EN(I),I=1,50)
      !WRITE 100 11,(EN(I),I=51,56)
      !WRITE 100 2
 2    FORMAT(2X,'energy spect. of neutrons, transv.: ENT',
     *' (0-25MeV, step-0.5MeV)')
      !WRITE 100 10,(ENT(I),I=1,50)
      !WRITE 100 11,(ENT(I),I=51,56)
      !WRITE 100 3
 3    FORMAT(2X,'energy spect. of protons, total: EP',
     *' (0-25MeV, step-0.5MeV)')
      !WRITE 100 10,(EP(I),I=1,50)
      !WRITE 100 11,(EP(I),I=51,56)
      !WRITE 100 4
 4    FORMAT(2X,'energy spect. of protons, transv.: EPT',
     *' (0-25MeV, step-0.5MeV)')
      !WRITE 100 10,(EPT(I),I=1,50)
      !WRITE 100 11,(EPT(I),I=51,56)
      !WRITE 100 7
 7    FORMAT(2X,'energy spect. of deutrons, total: EDE',
     *' (0-25MeV, step-0.5MeV)')
      !WRITE 100 10,(EDE(I),I=1,50)
      !WRITE 100 11,(EDE(I),I=51,56)
      !WRITE 100 8
 8    FORMAT(2X,'energy spect. of tritons, total: ETR',
     *' (0-25MeV, step-0.5MeV)')
      !WRITE 100 10,(ETR(I),I=1,50)
      !WRITE 100 11,(ETR(I),I=51,56)
      !WRITE 100 9
 9    FORMAT(2X,'energy spect. of HE-3, total: EHE3',
     *' (0-100MeV, step-2MeV)')
      !WRITE 100 10,(EHE3(I),I=1,50)
      !WRITE 100 11,(EHE3(I),I=51,56)
      !WRITE 100 5
 5    FORMAT(2X,'energy spect. of alphas, total: EA',
     *' (0-100MeV, step-2MeV)')
      !WRITE 100 10,(EA(I),I=1,50)
      !WRITE 100 11,(EA(I),I=51,56)
      !WRITE 100 6
 6    FORMAT(2X,'energy spect. of alphas, transv.: EAT',
     *' (0-100MeV, step-2MeV)')
      !WRITE 100 10,(EAT(I),I=1,50)
      !WRITE 100 11,(EAT(I),I=51,56)
      RETURN
      END
      SUBROUTINE CSPEC1(KST1)
      COMMON /BLOKC/SPT(10,500)
      COMMON /BLSPE1/EN(56),ENT(56),EP(56),EPT(56),EA(56),EAT(56)
      COMMON /BLSPE2/EDE(56),ETR(56),EHE3(56)
      COMMON /BLSPE3/SEN,SEP,SED,SET,SEH,SEA,NN,NP,ND,NT,NH,NA
      DO 1 I=1,KST1
      IA=INT(SPT(9,I)/0.94+0.5)
      IZ=INT(SPT(8,I)+0.5)
      EK=1000.*SPT(7,I)
      COT=SPT(4,I)
      SIT2=1.-COT*COT
      ET=EK*SIT2
      IF(IA.EQ.1.AND.IZ.EQ.0) NN=NN+1
      IF(IA.EQ.1.AND.IZ.EQ.0) SEN=SEN+EK
      IF(IA.EQ.1.AND.IZ.EQ.1) NP=NP+1
      IF(IA.EQ.1.AND.IZ.EQ.1) SEP=SEP+EK
      IF(IA.EQ.2.AND.IZ.EQ.1) ND=ND+1
      IF(IA.EQ.2.AND.IZ.EQ.1) SED=SED+EK
      IF(IA.EQ.3.AND.IZ.EQ.1) NT=NT+1
      IF(IA.EQ.3.AND.IZ.EQ.1) SET=SET+EK
      IF(IA.EQ.3.AND.IZ.EQ.2) NH=NH+1
      IF(IA.EQ.3.AND.IZ.EQ.2) SEH=SEH+EK
      IF(IA.EQ.4.AND.IZ.EQ.2) NA=NA+1
      IF(IA.EQ.4.AND.IZ.EQ.2) SEA=SEA+EK
      IF(IA.EQ.1.AND.IZ.EQ.0)
     *CALL HIST1(EK,0.,25.,0.5,EN,56,1.)
      IF(IA.EQ.1.AND.IZ.EQ.0)
     *CALL HIST1(ET,0.,25.,0.5,ENT,56,1.)
      IF(IA.EQ.1.AND.IZ.EQ.1)
     *CALL HIST1(EK,0.,25.,0.5,EP,56,1.)
      IF(IA.EQ.1.AND.IZ.EQ.1)
     *CALL HIST1(ET,0.,25.,0.5,EPT,56,1.)
      IF(IA.EQ.4.AND.IZ.EQ.2)
     *CALL HIST1(EK,0.,100.,2.,EA,56,1.)
      IF(IA.EQ.4.AND.IZ.EQ.2)
     *CALL HIST1(ET,0.,100.,2.,EAT,56,1.)
      IF(IA.EQ.3.AND.IZ.EQ.2)
     *CALL HIST1(EK,0.,100.,2.,EHE3,56,1.)
      IF(IA.EQ.2.AND.IZ.EQ.1)
     *CALL HIST1(EK,0.,25.,0.5,EDE,56,1.)
      IF(IA.EQ.3.AND.IZ.EQ.1)
     *CALL HIST1(EK,0.,25.,0.5,ETR,56,1.)
 1    CONTINUE
      RETURN
      END
      SUBROUTINE HOTFR(IP0,A0,Z0,AM,CH,PN0,TEMP,ECOLB)
c --- for hot fragments in the break-up vo;ume: their number,
c --- positions and the number of accompanied light part.
      DIMENSION XX(500),YY(500),ZZ(500),AM(500),CH(500),PN0(3,500)
      COMMON /BLPAT/PA(500),PAZ(500),IPT
      COMMON /BLCEN/XC(500),YC(500),ZC(500)
      COMMON /BLHFR1/YIMF(20),YZ710(20),DYZ710(20)
      COMMON /BLHFR2/YZT,RZT,YZF(20),RZF(20),DRZF(20)
      COMMON /BLHFR3/YLF(20),DYLF(20),YHF(20),DYHF(20),YSF(20)
      COMMON /BLHFR4/ZLF(20),ZHF(20),ZSF(20)
      COMMON /BFHFR5/CULC(20),TEMC(20),EF(20),ELF(20),E71(20),EHF(20)
      COMMON /BLHFR6/ALF(20),A71(20),AHF(20),ELFN(20),E71N(20),EHFN(20)
      DO 1 I=1,IPT
      XX(I)=XC(I)
      YY(I)=YC(I)
      ZZ(I)=ZC(I)
 1    CONTINUE
      ILF=0
      IHF=0
      ISF=0
      IMF=0
      IZ710=0
      EF0=0.
      EL0=0.
      E70=0.
      EH0=0.
      AL0=0.
      A70=0.
      AH0=0.
      DO 40 I=1,IPT
      IZ=INT(PAZ(I)+0.5)
      IA=INT(PA(I)+0.5)
      IF(IZ.GE.1.AND.IZ.LE.6) ILF=ILF+1
      IF(IZ.GE.1.AND.IZ.LE.6) AL0=AL0+IA
      IF(IZ.GE.11) IHF=IHF+1
      IF(IZ.GE.11) AH0=AH0+IA
      IF(IZ.GE.7.AND.IZ.LE.10) IZ710=IZ710+1
      IF(IZ.GE.7.AND.IZ.LE.10) A70=A70+IA
      IF(IZ.GE.5.AND.IZ.LE.22) IMF=IMF+1
      IF(IZ.GT.22) ISF=ISF+1
      EN0=(PN0(1,I)**2+PN0(2,I)**2+PN0(3,I)**2)/(2.*AM(I))
      IF(IZ.GE.1) EF0=EF0+EN0
      IF(IZ.GE.1.AND.IZ.LE.6) EL0=EL0+EN0
      IF(IZ.GE.7.AND.IZ.LE.10) E70=E70+EN0
      IF(IZ.GE.11) EH0=EH0+EN0
 40   CONTINUE
      IF(IZ710.EQ.0) GO TO 10
      YSF(IMF+1)=YSF(IMF+1)+ISF
      YIMF(IMF+1)=YIMF(IMF+1)+1.
      YLF(IMF+1)=YLF(IMF+1)+ILF
      DYLF(IMF+1)=DYLF(IMF+1)+ILF*ILF
      YZ710(IMF+1)=YZ710(IMF+1)+IZ710
      DYZ710(IMF+1)=DYZ710(IMF+1)+IZ710*IZ710
      YHF(IMF+1)=YHF(IMF+1)+IHF
      DYHF(IMF+1)=DYHF(IMF+1)+IHF*IHF
      CULC(IMF+1)=CULC(IMF+1)+ECOLB
      TEMC(IMF+1)=TEMC(IMF+1)+TEMP
      EF(IMF+1)=EF(IMF+1)+EF0
      ELF(IMF+1)=ELF(IMF+1)+EL0
      E71(IMF+1)=E71(IMF+1)+E70
      EHF(IMF+1)=EHF(IMF+1)+EH0
      ALF(IMF+1)=ALF(IMF+1)+AL0
      A71(IMF+1)=A71(IMF+1)+A70
      AHF(IMF+1)=AHF(IMF+1)+AH0
      DO 41 I=1,IPT
      RR=SQRT(XX(I)**2+YY(I)**2+ZZ(I)**2)
      IZ=INT(PAZ(I)+0.5)
      IA=INT(PA(I)+0.5)
      IF(IZ.GT.22) ZSF(IMF+1)=ZSF(IMF+1)+IZ
      IF(IZ.GE.1.AND.IZ.LE.6) ZLF(IMF+1)=ZLF(IMF+1)+IZ
      IF(IZ.GT.10) ZHF(IMF+1)=ZHF(IMF+1)+IZ
      IF(IZ.GE.7.AND.IZ.LE.10) YZF(IMF+1)=YZF(IMF+1)+1.
      IF(IZ.GE.7.AND.IZ.LE.10) RZF(IMF+1)=RZF(IMF+1)+RR
      IF(IZ.GE.7.AND.IZ.LE.10) DRZF(IMF+1)=DRZF(IMF+1)+RR*RR
      IF(IZ.GE.7.AND.IZ.LE.10) RZT=RZT+RR
      IF(IZ.GE.7.AND.IZ.LE.10) YZT=YZT+1.
 41   CONTINUE
 10   CONTINUE
      RETURN
      END
      SUBROUTINE HOTFRB
      COMMON /BLHFR1/YIMF(20),YZ710(20),DYZ710(20)
      COMMON /BLHFR2/YZT,RZT,YZF(20),RZF(20),DRZF(20)
      COMMON /BLHFR3/YLF(20),DYLF(20),YHF(20),DYHF(20),YSF(20)
      COMMON /BLHFR4/ZLF(20),ZHF(20),ZSF(20)
      COMMON /BFHFR5/CULC(20),TEMC(20),EF(20),ELF(20),E71(20),EHF(20)
      COMMON /BLHFR6/ALF(20),A71(20),AHF(20),ELFN(20),E71N(20),EHFN(20)
      YZT=0.
      RZT=0.
      DO 1 I=1,20
      ELFN(I)=0.
      E71N(I)=0.
      EHFN(I)=0.
      ALF(I)=0.
      A71(I)=0.
      AHF(I)=0.
      CULC(I)=0.
      TEMC(I)=0.
      EF(I)=0.
      ELF(I)=0.
      E71(I)=0.
      EHF(I)=0.
      YIMF(I)=0.
      YHF(I)=0.
      DYHF(I)=0.
      YLF(I)=0.
      DYLF(I)=0.
      YZ710(I)=0.
      DYZ710(I)=0.
      YSF(I)=0.
      YZF(I)=0.
      RZF(I)=0.
      DRZF(I)=0.
      ZSF(I)=0.
      ZHF(I)=0.
      ZLF(I)=0.
 1    CONTINUE
      RETURN
      END
      SUBROUTINE HOTFRE
      COMMON /BLHFR1/YIMF(20),YZ710(20),DYZ710(20)
      COMMON /BLHFR2/YZT,RZT,YZF(20),RZF(20),DRZF(20)
      COMMON /BLHFR3/YLF(20),DYLF(20),YHF(20),DYHF(20),YSF(20)
      COMMON /BLHFR4/ZLF(20),ZHF(20),ZSF(20)
      COMMON /BFHFR5/CULC(20),TEMC(20),EF(20),ELF(20),E71(20),EHF(20)
      COMMON /BLHFR6/ALF(20),A71(20),AHF(20),ELFN(20),E71N(20),EHFN(20)
      RSYS=1.82*1.17*(160.**0.3333333)
      IF(YZT.GT.0.) RZT=RZT/(YZT*RSYS)
      DO 1 I=1,20
      IF(ALF(I).GT.0.) ELFN(I)=ELF(I)/ALF(I)
      IF(AHF(I).GT.0.) EHFN(I)=EHF(I)/AHF(I)
      IF(A71(I).GT.0.) E71N(I)=E71(I)/A71(I)
      IF(YLF(I).GT.0.) ELF(I)=ELF(I)/YLF(I)
      IF(YHF(I).GT.0.) EHF(I)=EHF(I)/YHF(I)
      IF(YZ710(I).GT.0.) E71(I)=E71(I)/YZ710(I)
      IF(YSF(I).GT.0.) ZSF(I)=ZSF(I)/YSF(I)
      IF(YLF(I).GT.0.) ZLF(I)=ZLF(I)/YLF(I)
      IF(YHF(I).GT.0.) ZHF(I)=ZHF(I)/YHF(I)
      IF(YIMF(I).GT.0.) YLF(I)=YLF(I)/YIMF(I)
      IF(YIMF(I).GT.0.) DYLF(I)=DYLF(I)/YIMF(I)
      DDD=DYLF(I)-YLF(I)**2
      IF(DDD.GT.0.) DYLF(I)=SQRT(DDD)
      IF(DDD.LE.0.) DYLF(I)=0.
      IF(YIMF(I).GT.0.) YHF(I)=YHF(I)/YIMF(I)
      IF(YIMF(I).GT.0.) DYHF(I)=DYHF(I)/YIMF(I)
      DDD=DYHF(I)-YHF(I)**2
      IF(DDD.GT.0.) DYHF(I)=SQRT(DDD)
      IF(DDD.LE.0.) DYHF(I)=0.
      IF(YZF(I).GT.0.) RZF(I)=RZF(I)/(YZF(I)*RSYS)
      IF(YZF(I).GT.0.) DRZF(I)=DRZF(I)/(YZF(I)*RSYS*RSYS)
      DDD=DRZF(I)-RZF(I)**2
      IF(DDD.GT.0.) DRZF(I)=SQRT(DDD)
      IF(DDD.LE.0.) DRZF(I)=0.
      IF(YIMF(I).GT.0.) YZ710(I)=YZ710(I)/YIMF(I)
      IF(YIMF(I).GT.0.) DYZ710(I)=DYZ710(I)/YIMF(I)
      DDD=DYZ710(I)-YZ710(I)**2
      IF(DDD.GT.0.) DYZ710(I)=SQRT(DDD)
      IF(DDD.LE.0.) DYZ710(I)=0.
      IF(YIMF(I).GT.0.) YSF(I)=YSF(I)/YIMF(I)
      IF(YIMF(I).GT.0.) CULC(I)=CULC(I)/YIMF(I)
      IF(YIMF(I).GT.0.) TEMC(I)=TEMC(I)/YIMF(I)
      IF(YIMF(I).GT.0.) EF(I)=EF(I)/YIMF(I)
      IF(YIMF(I).GT.0.) ALF(I)=ALF(I)/YIMF(I)
      IF(YIMF(I).GT.0.) AHF(I)=AHF(I)/YIMF(I)
      IF(YIMF(I).GT.0.) A71(I)=A71(I)/YIMF(I)
 1    CONTINUE
      DO 15 I=1,20
      IF(EF(I).GT.0.) ELF(I)=ELF(I)/EF(I)
      IF(EF(I).GT.0.) E71(I)=E71(I)/EF(I)
      IF(EF(I).GT.0.) EHF(I)=EHF(I)/EF(I)
      IF(EF(I).GT.0.) ELFN(I)=ELFN(I)/EF(I)
      IF(EF(I).GT.0.) E71N(I)=E71N(I)/EF(I)
      IF(EF(I).GT.0.) EHFN(I)=EHFN(I)/EF(I)
 15   CONTINUE
      !WRITE 100 2
 2    FORMAT(2X,' **** for hot fragments ****')
      !WRITE 100 4,YZT,RZT
 4    FORMAT(2X,'total number of "Z=7-10" YZT=',F9.0,
     *'  their mean radius RZT/RSYS=',F9.4)
      !WRITE 100 6
 6    FORMAT(1X,'#IMF',3X,'YIMF',4X,'YLF',3X,'DYLF',
     *3X,'YZ710',1X,'DYZ710',5X,'YZF',4X,'RZF',3X,'DRZF',
     *5X,'YHF',3X,'DYHF')
      DO 5 I=1,20
      IMF0=I-1
      !WRITE 100 7,IMF0,YIMF(I),YLF(I),DYLF(I),YZ710(I),DYZ710(I),
!     *YZF(I),RZF(I),DRZF(I),YHF(I),DYHF(I)
 7    FORMAT(1X,I4,1X,F6.0,1X,F6.3,1X,F6.3,2X,F6.3,1X,F6.3,
     *2X,F6.0,1X,F6.3,1X,F6.3,2X,F6.3,1X,F6.3)
 5    CONTINUE
      !WRITE 100 8
 8    FORMAT(2X,' *****************************')
      !WRITE 100 9
 9    FORMAT(1X,'#IMF',3X,'TEMC',3X,'CULC',5X,'EF',
     *5X,'ELF',4X,'E71',4X,'EHF',4X,'ZLF',3X,'ZHF',3X,'ZSF',
     *3X,'YSF')
      DO 10 I=1,20
      IMF0=I-1
      !WRITE 100 11,IMF0,TEMC(I),CULC(I),EF(I),ELF(I),E71(I),EHF(I),
!     *ZLF(I),ZHF(I),ZSF(I),YSF(I)
 11   FORMAT(1X,I4,1X,F6.2,1X,F6.3,1X,F6.3,2X,F6.3,1X,F6.3,
     *1X,F6.3,2X,F5.2,1X,F5.2,1X,F5.2,1X,F5.3)
 10   CONTINUE
      !WRITE 100 8
      !WRITE 100 14
 14   FORMAT(1X,'#IMF',5X,'ALF',4X,'A71',4X,'AHF',
     *4X,'ELFN',3X,'E71N',3X,'EHFN')
      DO 12 I=1,20
      IMF0=I-1
      !WRITE 100 13,IMF0,ALF(I),A71(I),AHF(I),ELFN(I),E71N(I),EHFN(I)
 13   FORMAT(1X,I4,2X,F6.1,1X,F6.1,1X,F6.1,2X,F6.4,1X,F6.4,
     *1X,F6.4)
 12   CONTINUE
      !WRITE 100 8
      RETURN
      END
      SUBROUTINE HOTF6B
c --- for hot fragments (final as Z=6) in the break-up volume
c  (and for accompany fragments: Z=1-5, Z>6, Z>20 in case of IMFs: Z=3-20)
      COMMON /BLHF61/YIMF(20),YZ6(20),DYZ6(20)
      COMMON /BLHF62/Y6T,R6T,Y6F(20),R6F(20),DR6F(20)
      COMMON /BLHF63/YLF(20),DYLF(20),YHF(20),DYHF(20),YSF(20)
      COMMON /BLHF64/ZLF(20),ZHF(20),ZSF(20),Z6F(20)
      COMMON /BLHF65/CULC(20),TEMC(20),EF(20),ELF(20),E6F(20),EHF(20)
      COMMON /BLHF66/ALF(20),A6F(20),AHF(20),ELFN(20),E6FN(20),EHFN(20)
      COMMON /BLHF67/YCLF(20),DYCLF(20),YCHF(20),DYCHF(20),YCSF(20)
      COMMON /BLACC1/NACLF(20,50),NAC6F(20,50),NACHF(20,50)
      COMMON /BLACC0/ACLF(20,50),AC6F(20,50),ACHF(20,50)
      Y6T=0.
      R6T=0.
      DO 1 I=1,20
      ELFN(I)=0.
      E6FN(I)=0.
      EHFN(I)=0.
      ALF(I)=0.
      A6F(I)=0.
      AHF(I)=0.
      CULC(I)=0.
      TEMC(I)=0.
      EF(I)=0.
      ELF(I)=0.
      E6F(I)=0.
      EHF(I)=0.
      YIMF(I)=0.
      YHF(I)=0.
      DYHF(I)=0.
      YCHF(I)=0.
      DYCHF(I)=0.
      YLF(I)=0.
      DYLF(I)=0.
      YCLF(I)=0.
      DYCLF(I)=0.
      YZ6(I)=0.
      DYZ6(I)=0.
      YSF(I)=0.
      YCSF(I)=0.
      Y6F(I)=0.
      R6F(I)=0.
      DR6F(I)=0.
      ZSF(I)=0.
      ZHF(I)=0.
      ZLF(I)=0.
      Z6F(I)=0.
      DO 2 J=1,50
      ACLF(I,J)=0.
      AC6F(I,J)=0.
      ACHF(I,J)=0.
      NACLF(I,J)=0
      NAC6F(I,J)=0
      NACHF(I,J)=0
 2    CONTINUE
 1    CONTINUE
      RETURN
      END
      SUBROUTINE HOTF6E
c --- for hot fragments (final as Z=6) in the break-up volume
c  (and for accompany fragments: Z=1-5, Z>6, Z>20 in case of IMFs: Z=3-20)
      COMMON /BLHF61/YIMF(20),YZ6(20),DYZ6(20)
      COMMON /BLHF62/Y6T,R6T,Y6F(20),R6F(20),DR6F(20)
      COMMON /BLHF63/YLF(20),DYLF(20),YHF(20),DYHF(20),YSF(20)
      COMMON /BLHF64/ZLF(20),ZHF(20),ZSF(20),Z6F(20)
      COMMON /BLHF65/CULC(20),TEMC(20),EF(20),ELF(20),E6F(20),EHF(20)
      COMMON /BLHF66/ALF(20),A6F(20),AHF(20),ELFN(20),E6FN(20),EHFN(20)
      COMMON /BLHF67/YCLF(20),DYCLF(20),YCHF(20),DYCHF(20),YCSF(20)
      COMMON /BLACC1/NACLF(20,50),NAC6F(20,50),NACHF(20,50)
      COMMON /BLACC0/ACLF(20,50),AC6F(20,50),ACHF(20,50)
      COMMON /BLACCL/ACL(50,500),ITN(50),JJN,NSIGN,ZSIGN(500)
      RSYS=1.82*1.17*(160.**0.3333333)
      IF(Y6T.GT.0.) R6T=R6T/(Y6T*RSYS)
      DO 1 I=1,20
      IF(ALF(I).GT.0.) ELFN(I)=ELF(I)/ALF(I)
      IF(AHF(I).GT.0.) EHFN(I)=EHF(I)/AHF(I)
      IF(A6F(I).GT.0.) E6FN(I)=E6F(I)/A6F(I)
      IF(YLF(I).GT.0.) ELF(I)=ELF(I)/YLF(I)
      IF(YHF(I).GT.0.) EHF(I)=EHF(I)/YHF(I)
      IF(Y6F(I).GT.0.) E6F(I)=E6F(I)/Y6F(I)
      IF(YSF(I).GT.0.) ZSF(I)=ZSF(I)/YSF(I)
      IF(YLF(I).GT.0.) ZLF(I)=ZLF(I)/YLF(I)
      IF(Y6F(I).GT.0.) Z6F(I)=Z6F(I)/Y6F(I)
      IF(YHF(I).GT.0.) ZHF(I)=ZHF(I)/YHF(I)
      IF(YIMF(I).GT.0.) YLF(I)=YLF(I)/YIMF(I)
      IF(YIMF(I).GT.0.) YCLF(I)=YCLF(I)/YIMF(I)
      IF(YIMF(I).GT.0.) DYCLF(I)=DYCLF(I)/YIMF(I)
      DDD=DYCLF(I)-YCLF(I)**2
      IF(DDD.GT.0.) DYCLF(I)=SQRT(DDD)
      IF(DDD.LE.0.) DYCLF(I)=0.
      IF(YIMF(I).GT.0.) YHF(I)=YHF(I)/YIMF(I)
      IF(YIMF(I).GT.0.) YCHF(I)=YCHF(I)/YIMF(I)
      IF(YIMF(I).GT.0.) DYCHF(I)=DYCHF(I)/YIMF(I)
      DDD=DYCHF(I)-YCHF(I)**2
      IF(DDD.GT.0.) DYCHF(I)=SQRT(DDD)
      IF(DDD.LE.0.) DYCHF(I)=0.
      IF(Y6F(I).GT.0.) R6F(I)=R6F(I)/(Y6F(I)*RSYS)
      IF(Y6F(I).GT.0.) DR6F(I)=DR6F(I)/(Y6F(I)*RSYS*RSYS)
      DDD=DR6F(I)-R6F(I)**2
      IF(DDD.GT.0.) DR6F(I)=SQRT(DDD)
      IF(DDD.LE.0.) DR6F(I)=0.
      IF(YIMF(I).GT.0.) YZ6(I)=YZ6(I)/YIMF(I)
      IF(YIMF(I).GT.0.) DYZ6(I)=DYZ6(I)/YIMF(I)
      DDD=DYZ6(I)-YZ6(I)**2
      IF(DDD.GT.0.) DYZ6(I)=SQRT(DDD)
      IF(DDD.LE.0.) DYZ6(I)=0.
      IF(YIMF(I).GT.0.) Y6F(I)=Y6F(I)/YIMF(I)
      IF(YIMF(I).GT.0.) YSF(I)=YSF(I)/YIMF(I)
      IF(YIMF(I).GT.0.) YCSF(I)=YCSF(I)/YIMF(I)
      IF(YIMF(I).GT.0.) CULC(I)=CULC(I)/YIMF(I)
      IF(YIMF(I).GT.0.) TEMC(I)=TEMC(I)/YIMF(I)
      IF(YIMF(I).GT.0.) EF(I)=EF(I)/YIMF(I)
      IF(YIMF(I).GT.0.) ALF(I)=ALF(I)/YIMF(I)
      IF(YIMF(I).GT.0.) AHF(I)=AHF(I)/YIMF(I)
      IF(YIMF(I).GT.0.) A6F(I)=A6F(I)/YIMF(I)
 1    CONTINUE
      DO 15 I=1,20
      IF(EF(I).GT.0.) ELF(I)=ELF(I)/EF(I)
      IF(EF(I).GT.0.) E6F(I)=E6F(I)/EF(I)
      IF(EF(I).GT.0.) EHF(I)=EHF(I)/EF(I)
      IF(EF(I).GT.0.) ELFN(I)=ELFN(I)/EF(I)
      IF(EF(I).GT.0.) E6FN(I)=E6FN(I)/EF(I)
      IF(EF(I).GT.0.) EHFN(I)=EHFN(I)/EF(I)
 15   CONTINUE
      !WRITE 100 2
 2    FORMAT(2X,' **** for hot fragments /new - HOTFR6/ ****')
      !WRITE 100 4,Y6T,R6T
 4    FORMAT(2X,'tot.numb.of final "Z=6" Y6T=',F9.0,
     *'  their mean radius R6T/RSYS=',F9.4)
      !WRITE 100 6
 6    FORMAT(1X,'#IMF',3X,'YIMF',4X,'YLF',3X,'YCLF',
     *5X,'YZ6',3X,'DYZ6',5X,'Y6F',4X,'R6F',3X,'DR6F',
     *5X,'YHF',3X,'YCHF')
      DO 5 I=1,20
      IMF0=I-1
      !WRITE 100 7,IMF0,YIMF(I),YLF(I),YCLF(I),YZ6(I),DYZ6(I),
!     *Y6F(I),R6F(I),DR6F(I),YHF(I),YCHF(I)
 7    FORMAT(1X,I4,1X,F6.0,1X,F6.3,1X,F6.3,2X,F6.3,1X,F6.3,
     *2X,F6.2,1X,F6.3,1X,F6.3,2X,F6.3,1X,F6.3)
 5    CONTINUE
      !WRITE 100 8
 8    FORMAT(2X,' *****************************')
      !WRITE 100 9
 9    FORMAT(1X,'#IMF',3X,'TEMC',3X,'CULC',5X,'EF',
     *5X,'ELF',4X,'E6F',4X,'EHF',4X,'ZLF',3X,'Z6F',3X,'ZHF',
     *3X,'ZSF',3X,'YSF')
      DO 10 I=1,20
      IMF0=I-1
      !WRITE 100 11,IMF0,TEMC(I),CULC(I),EF(I),ELF(I),E6F(I),EHF(I),
!     *ZLF(I),Z6F(I),ZHF(I),ZSF(I),YSF(I)
 11   FORMAT(1X,I4,1X,F6.2,1X,F6.3,1X,F6.3,2X,F6.3,1X,F6.3,
     *1X,F6.3,2X,F5.2,1X,F5.2,1X,F5.2,1X,F5.2,1X,F5.3)
 10   CONTINUE
      !WRITE 100 8
      !WRITE 100 14
 14   FORMAT(1X,'#IMF',5X,'ALF',4X,'A6F',4X,'AHF',
     *4X,'ELFN',3X,'E6FN',3X,'EHFN')
      DO 12 I=1,20
      IMF0=I-1
      !WRITE 100 13,IMF0,ALF(I),A6F(I),AHF(I),ELFN(I),E6FN(I),EHFN(I)
 13   FORMAT(1X,I4,2X,F6.1,1X,F6.1,1X,F6.1,2X,F6.4,1X,F6.4,
     *1X,F6.4)
 12   CONTINUE
      !WRITE 100 16
 16   FORMAT(5X,'******* for accelration ********')
      !WRITE 100 21,(ITN(I),I=1,50)
 21   FORMAT(1X,'time ITN(1-50)=',10I6)
      DO 17 I=1,20
      DO 17 J=1,50
      IF(NACLF(I,J).GT.0) ACLF(I,J)=100*ACLF(I,J)/FLOAT(NACLF(I,J))
      IF(NAC6F(I,J).GT.0) AC6F(I,J)=100*AC6F(I,J)/FLOAT(NAC6F(I,J))
      IF(NACHF(I,J).GT.0) ACHF(I,J)=100*ACHF(I,J)/FLOAT(NACHF(I,J))
 17   CONTINUE
 19   FORMAT(1X,6(I4,1X,F7.5,1X))
      !WRITE 100 20
 20   FORMAT(2X,'NACLF and avr.accel. ACLF*100 for dif.times',
     *' at IMF=1,...,6')
      DO 18 J=1,50
      !WRITE 100 19,(NACLF(I,J),ACLF(I,J),I=2,7)
 18   CONTINUE
      !WRITE 100 22
 22   FORMAT(2X,'NAC6F and avr.accel. AC6F*100 for dif.times',
     *' at IMF=1,...,6')
      DO 23 J=1,50
      !WRITE 100 19,(NAC6F(I,J),AC6F(I,J),I=2,7)
 23   CONTINUE
      !WRITE 100 24
 24   FORMAT(2X,'NACHF and avr.accel. ACHF*100 for dif.times',
     *' at IMF=1,...,6')
      DO 25 J=1,50
      !WRITE 100 19,(NACHF(I,J),ACHF(I,J),I=2,7)
 25   CONTINUE
      !WRITE 100 8
      RETURN
      END
      SUBROUTINE HOTFR6(IPT,KST1,A0,Z0,AM,CH,PN0,TEMP,ECOLB,EXPA)
c --- for hot fragments (final as Z=6) in the break-up volume
c  (and for accompany fragments: Z=1-5, Z>6, Z>20 in case of IMFs: Z=3-20)
      DIMENSION XX(500),YY(500),ZZ(500),AM(500),CH(500),PN0(3,500)
      DIMENSION PA(500),PAZ(500),EXPA(500)
      COMMON /BLCEN/XC(500),YC(500),ZC(500) /BLHKST/IHFKST(500,2)
      COMMON /BLOKC/SPTU(10,500)
      COMMON /BLHF61/YIMF(20),YZ6(20),DYZ6(20)
      COMMON /BLHF62/Y6T,R6T,Y6F(20),R6F(20),DR6F(20)
      COMMON /BLHF63/YLF(20),DYLF(20),YHF(20),DYHF(20),YSF(20)
      COMMON /BLHF64/ZLF(20),ZHF(20),ZSF(20),Z6F(20)
      COMMON /BLHF65/CULC(20),TEMC(20),EF(20),ELF(20),E6F(20),EHF(20)
      COMMON /BLHF66/ALF(20),A6F(20),AHF(20),ELFN(20),E6FN(20),EHFN(20)
      COMMON /BLHF67/YCLF(20),DYCLF(20),YCHF(20),DYCHF(20),YCSF(20)
      COMMON /BLACC0/ACLF(20,50),AC6F(20,50),ACHF(20,50)
      COMMON /BLACC1/NACLF(20,50),NAC6F(20,50),NACHF(20,50)
      COMMON /BLACCL/ACL(50,500),ITN(50),JJN,NSIGN,ZSIGN(500)
c -----------------
c      !WRITE 100 26,JJN,NSIGN
c 26   FORMAT(2X,'hotfr6; c.accel.:  JJN=',I4,' NSIGN=',I4)
c      !WRITE 100 27,(ITN(I),I=1,50)
c 27   FORMAT(1X,'times:',15I4)
c      !WRITE 100 28,(ZSIGN(I),I=1,NSIGN)
c 28   FORMAT(1X,'ch.part.:',15F4.0)
c -----------------
c      IAT=0
c      IZT=0
c      IZHT=0
c      EKINT=0.
c      EKINH=0.
c -----------------
      IAHT=0
      RCMX=0.
      RCMY=0.
      RCMZ=0.
      DO 1 I=1,IPT
      XX(I)=XC(I)
      YY(I)=YC(I)
      ZZ(I)=ZC(I)
      IAH=INT(AM(I)/0.94+0.5)
      IZH=INT(CH(I)+0.5)
      PA(I)=IAH
      PAZ(I)=IZH
      IAHT=IAHT+IAH
      RCMX=RCMX+IAH*XX(I)
      RCMY=RCMY+IAH*YY(I)
      RCMZ=RCMZ+IAH*ZZ(I)
c ------------------
c      EN0=(PN0(1,I)**2+PN0(2,I)**2+PN0(3,I)**2)/(2.*AM(I))
c      EKINH=EKINH+EN0
c      IZHT=IZHT+IZH
c      !WRITE 100 22,IAH,IZH
c 22   FORMAT(2X,'hotfr6:  IAH,IZH=',2I4)
c ------------------
 1    CONTINUE
      RCMX=RCMX/FLOAT(IAHT)
      RCMY=RCMY/FLOAT(IAHT)
      RCMZ=RCMZ/FLOAT(IAHT)
c ------------------
c      RCM=SQRT(RCMX**2+RCMY**2+RCMZ**2)
c      !WRITE 100 21,IAHT,IZHT,EKINH,RCM
c 21   FORMAT(2X,'hotfr6:  IAHT,IZHT=',2I4,' EKINH(gev)=',F7.4,
c     *' RCM(fm)=',F6.2)
c ------------------
      IMF=0
      IZ6=0
      DO 2 I=1,KST1
      IZ=INT(SPTU(8,I)+0.5)
      IA=INT(SPTU(9,I)/0.94+0.5)
c ----------------------
c      EK=SPTU(7,I)
c      EKINT=EKINT+EK
c      IAT=IAT+IA
c      IZT=IZT+IZ
c      !WRITE 100 20,IA,IZ,EK
c 20   FORMAT(2X,'hotfr6:  IA,IZ=',2I4,' EK=',F7.4)
c ----------------------
      IF(IZ.GE.3.AND.IZ.LE.20) IMF=IMF+1
      IF(IZ.EQ.6) IZ6=IZ6+1
 2    CONTINUE
c ----------------------
c      !WRITE 100 23,IAT,IZT,IMF,IZ6,EKINT
c 23   FORMAT(2X,'hotfr6:  IAT,IZT=',2I4,' IMF=',I3,' IZ6=',I3,
c     *' EKINT(gev)=',F8.4)
c ----------------------
      IF(IZ6.EQ.0.OR.JJN.EQ.0) GO TO 4
      NILF=0
      NI6F=0
      NIHF=0
      NISF=0
      EF0=0.
      EL0=0.
      E60=0.
      EH0=0.
      AL0=0.
      A60=0.
      AH0=0.
      DO 5 I=1,IPT
      IAH=INT(PA(I)+0.5)
      IZH=INT(PAZ(I)+0.5)
      EX0=EXPA(I)
      EN0=(PN0(1,I)**2+PN0(2,I)**2+PN0(3,I)**2)/(2.*AM(I))
      RR=SQRT((XX(I)-RCMX)**2+(YY(I)-RCMY)**2+(ZZ(I)-RCMZ)**2)
      J1=IHFKST(I,1)
      J2=IHFKST(I,2)
c -------------------
c      !WRITE 100 25,IAH,IZH,EN0,RR,J1,J2
c 25   FORMAT(7X,'primar.fragm.:  IAH,IZH=',2I4,' EN0=',F7.4,
c     *' RR(fm)=',F5.2,' J1,J2=',2I4)
c -------------------
      ILF=0
      I6F=0
      IHF=0
      ISF=0
      DO 6 J=J1,J2
      IZ=INT(SPTU(8,J)+0.5)
      IA=INT(SPTU(9,J)/0.94+0.5)
      EK=SPTU(7,J)
c -------------------
c      !WRITE 100 24,IA,IZ,EK
c 24   FORMAT(20X,'second.fragm.:  IA,IZ=',2I4,' EK=',F7.4)
c -------------------
      IF(IZ.GE.1.AND.IZ.LE.5) ILF=ILF+1
      IF(IZ.EQ.6) I6F=I6F+1
      IF(IZ.GE.7) IHF=IHF+1
      IF(IZ.GT.20) ISF=ISF+1
 6    CONTINUE
      NILF=NILF+ILF
      NI6F=NI6F+I6F
      NIHF=NIHF+IHF
      NISF=NISF+ISF
ccc      IF(ILF.GE.1.AND.I6F.EQ.0.AND.IHF.EQ.0) AL0=AL0+IAH
ccc      IF(I6F.GE.1) A60=A60+IAH
ccc      IF(I6F.EQ.0.AND.IHF.GE.1) AH0=AH0+IAH
      IF(ILF.GE.1.AND.I6F.EQ.0.AND.IHF.EQ.0) AL0=AL0+IAH
      IF(ILF.GE.1.AND.I6F.EQ.0.AND.IHF.EQ.0)
     *YLF(IMF+1)=YLF(IMF+1)+1.
      IF(I6F.GE.1.AND.IHF.EQ.0) A60=A60+IAH
      IF(I6F.GE.1.AND.IHF.EQ.0) Y6F(IMF+1)=Y6F(IMF+1)+1.
      IF(IHF.GE.1) AH0=AH0+IAH
      IF(IHF.GE.1) YHF(IMF+1)=YHF(IMF+1)+1.
      IF(ISF.GE.1) YSF(IMF+1)=YSF(IMF+1)+1.
      IF(ILF.GE.1.OR.I6F.GE.1.OR.IHF.GE.1) EF0=EF0+EN0
      IF(ILF.GE.1.AND.I6F.EQ.0.AND.IHF.EQ.0) EL0=EL0+EN0
      IF(I6F.GE.1.AND.IHF.EQ.0) E60=E60+EN0
      IF(IHF.GE.1) EH0=EH0+EN0
      IF(I6F.GE.1.AND.IHF.EQ.0) R6F(IMF+1)=R6F(IMF+1)+RR
      IF(I6F.GE.1.AND.IHF.EQ.0) DR6F(IMF+1)=DR6F(IMF+1)+RR*RR
      IF(I6F.GE.1.AND.IHF.EQ.0) R6T=R6T+RR
      IF(I6F.GE.1.AND.IHF.EQ.0) Y6T=Y6T+1.
      IF(ILF.GE.1.AND.I6F.EQ.0.AND.IHF.EQ.0)
     *ZLF(IMF+1)=ZLF(IMF+1)+IZH
      IF(I6F.GE.1.AND.IHF.EQ.0) Z6F(IMF+1)=Z6F(IMF+1)+IZH
      IF(IHF.GE.1) ZHF(IMF+1)=ZHF(IMF+1)+IZH
      IF(ISF.GE.1) ZSF(IMF+1)=ZSF(IMF+1)+IZH
c --- acceleration of fragments at different moments ---
      IF(JJN.EQ.0.OR.I.GT.NSIGN) GO TO 3
      IZH2=INT(ZSIGN(I)+0.5)
c ----------------
c      !WRITE 100 29,I,IZH,IZH2,ILF,I6F,IHF,ISF
c  29  FORMAT(2X,'before ACLF. I=',I2,' IZH,IZH2=',2I3,
c     *' ILF,I6F,IHF,ISF=',4I4)
c      !WRITE 100 30,(ACL(J,I),J=1,50)
c 30   FORMAT(2X,'ACL=',5E11.3)
c ----------------
      IF(IZH.NE.IZH2) GO TO 3
      DO 11 JT=1,50
      IF(ILF.GE.1.AND.I6F.EQ.0.AND.IHF.EQ.0)
     *ACLF(IMF+1,JT)=ACLF(IMF+1,JT)+ACL(JT,I)
      IF(ILF.GE.1.AND.I6F.EQ.0.AND.IHF.EQ.0)
     *NACLF(IMF+1,JT)=NACLF(IMF+1,JT)+1
      IF(I6F.GE.1.AND.IHF.EQ.0)
     *AC6F(IMF+1,JT)=AC6F(IMF+1,JT)+ACL(JT,I)
      IF(I6F.GE.1.AND.IHF.EQ.0)
     *NAC6F(IMF+1,JT)=NAC6F(IMF+1,JT)+1
      IF(IHF.GE.1) ACHF(IMF+1,JT)=ACHF(IMF+1,JT)+ACL(JT,I)
      IF(IHF.GE.1) NACHF(IMF+1,JT)=NACHF(IMF+1,JT)+1
 11   CONTINUE
 3    CONTINUE
c -----------------
 5    CONTINUE
      YIMF(IMF+1)=YIMF(IMF+1)+1.
      YCSF(IMF+1)=YCSF(IMF+1)+NISF
      YCLF(IMF+1)=YCLF(IMF+1)+NILF
      DYCLF(IMF+1)=DYCLF(IMF+1)+NILF*NILF
      YZ6(IMF+1)=YZ6(IMF+1)+NI6F
      DYZ6(IMF+1)=DYZ6(IMF+1)+NI6F*NI6F
      YCHF(IMF+1)=YCHF(IMF+1)+NIHF
      DYCHF(IMF+1)=DYCHF(IMF+1)+NIHF*NIHF
      CULC(IMF+1)=CULC(IMF+1)+ECOLB
      TEMC(IMF+1)=TEMC(IMF+1)+TEMP
      EF(IMF+1)=EF(IMF+1)+EF0
      ELF(IMF+1)=ELF(IMF+1)+EL0
      E6F(IMF+1)=E6F(IMF+1)+E60
      EHF(IMF+1)=EHF(IMF+1)+EH0
      ALF(IMF+1)=ALF(IMF+1)+AL0
      A6F(IMF+1)=A6F(IMF+1)+A60
      AHF(IMF+1)=AHF(IMF+1)+AH0
 4    CONTINUE
      RETURN
      END
      SUBROUTINE FIXT(T,K,M,AMP,ZMP)
      REAL K,LT
      COMMON /BLTMIC/TMICR,TMICR2,NMICR,TMICT
      COMMON /BLOK1/W0,TC,RADNCL,G0 /BLOK0/A0,Z0,V0,XZ0,IA0,E0,E0Q
      COMMON /BLOK3/E11(500),XZ(500),GA(500)
      COMMON /BLOK20/BT2,DBT2 /BLSMIC/SMIC
      DIMENSION AMP(500),ZMP(500),IAM(500),IZM(500)
      DO 1 I=1,M
      IZM(I)=INT(ZMP(I)+0.5)
      IAM(I)=INT(AMP(I)/0.94+0.5)
 1    CONTINUE
      A013=A0**0.333333
c --- keep constant or vary with multiplic. translational volume ---
c      VF=K*V0
      VF=V0*((1.+1.4*(M**0.333333-1)/(RADNCL*A013))**3-1.)
c ---
      M1=M-1
c      FAKT=1.
      FAKTL=0.
      IF(M1.LE.0) GO TO 35
      DO 32 I=1,M1
      FFF=1.
      IM=I+1
      DO 33 II=IM,M
      IF(IAM(I).EQ.IAM(II).AND.IZM(I).EQ.IZM(II)) FFF=FFF+1.
 33   CONTINUE
c      FAKT=FAKT*FFF
      FAKTL=FAKTL+ALOG(FFF)
 32   CONTINUE
 35   CONTINUE
c      PRFAKT=1.
c      PRFAKT=PRFAKT*FAKT
      PFAKTL=FAKTL
      PRGA=1.
      PRA32=1.
      DO 34 I=1,M
      IA=IAM(I)
      A=IA
      PRGA=PRGA*GA(IA)
      PRA32=PRA32*(A*SQRT(A))
 34   CONTINUE
      SCON=0.
      DO 62 I=1,M
      IA=IAM(I)
      A=IA
      IF(IA-4) 63,63,64
 63   SAIN=0.
      IF(IA.EQ.4) SAIN=2*T*4./E11(IA)
      GO TO 65
 64   SAIN=2*T*A/E11(IA)-DBT2*A**0.666667
 65   CONTINUE
      SCON=SCON+SAIN
 62   CONTINUE
      LT=16.15/SQRT(T)
c      STRAN=ALOG(PRA32/PRFAKT)+(M-1)*ALOG(VF/LT**3)+1.5*(M-1)-
      STRAN=ALOG(PRA32)-PFAKTL+(M-1)*ALOG(VF/LT**3)+1.5*(M-1)-
     *ALOG(A0*SQRT(A0))
      IF(STRAN.LT.0.) STRAN=0.
      SCON=SCON+ALOG(PRGA)+STRAN
      SMIC=SMIC+SCON
      TMICT=T
      TMICR=TMICR+T
      NMICR=NMICR+1
      TMICR2=TMICR2+T*T
      RETURN
      END
      SUBROUTINE FIXTB
      COMMON /BLTMIC/TMICR,TMICR2,NMICR,TMICT /BLSMIC/SMIC
      COMMON /BLICOR/ YL5_0,YL5_16,YL6_2,YL6_4,YB8_3,YB8_17
      SMIC=0.
      YL5_0=0.
      YL5_16=0.
      YL6_2=0.
      YL6_4=0.
      YB8_3=0.
      YB8_17=0.
      TMICR=0.
      NMICR=0.
      TMICR2=0.
      RETURN
      END
      SUBROUTINE FIXTE
      COMMON /BLTMIC/TMICR,TMICR2,NMICR,TMICT /BLSMIC/SMIC
      COMMON /BLICOR/ YL5_0,YL5_16,YL6_2,YL6_4,YB8_3,YB8_17
      COMMON /BENTRO/SSR
      IF(NMICR.GT.0) TMICR=TMICR/FLOAT(NMICR)
      IF(NMICR.GT.0) SMIC=SMIC/FLOAT(NMICR)
      IF(NMICR.GT.0) TMICR2=TMICR2/FLOAT(NMICR)
      DDD=TMICR2-TMICR*TMICR
      IF(DDD.LE.0.) TMICR2=0.
      IF(DDD.GT.0.) TMICR2=SQRT(DDD)
      !WRITE 100 1,NMICR,TMICR,TMICR2
 1    FORMAT(2X,'over NMICR=',I5,' events: microcan.temp. TMICR=',
     *F6.2,' MeV, var. TMICR2=',F6.2)
      !WRITE 100 2,SSR,SMIC
 2    FORMAT(2X,'entropy: SSR=',F8.3,'  microcan. SMIC=',F8.3)
      !WRITE 100 112, YL5_0,YL5_16,YL6_2,YL6_4,YB8_3,YB8_17
 112  FORMAT(1X,'excit.states: YL5_0,YL5_16=',F6.0,F5.0,
     *' YL6_2,YL6_4=',F6.0,F5.0,' YB8_3,YB8_17=',F6.0,F5.0)
      T_Li5=0.
      T_Li6=0.
      T_Be8=0.
      IF(YL5_0.LE.0.OR.YL5_16.LE.0.) GO TO 120
      T_Li5=16.66/ALOG(YL5_0/YL5_16)
 120  CONTINUE
      IF(YL6_2.LE.0.OR.YL6_4.LE.0.) GO TO 121
      T_Li6=(4.31-2.17)/ALOG((2.*YL6_2)/(3.*YL6_4))
 121  CONTINUE
      IF(YB8_3.LE.0.OR.YB8_17.LE.0.) GO TO 122
      T_Be8=(17.64-3.04)/ALOG((3.*YB8_3)/(5.*YB8_17))
 122  CONTINUE
      !WRITE 100 111, T_Li5, T_Li6, T_Be8
 111  FORMAT(2X,'temperature unstable states: T_Li5, T_Li6, T_Be8 =',
     *3F7.2,' (MeV)')
      RETURN
      END
c
      SUBROUTINE BFRENI
      COMMON /BFRE/ EHFR(500),SHFR(500),SEKHF,EHMAX,SHMAX
      DO 1 I=1,500
      EHFR(I)=0.
      SHFR(I)=0.
 1    CONTINUE
      SEKHF=0.
      EHMAX=0.
      SHMAX=0.
      RETURN
      END
      SUBROUTINE BFREN(INE,IP0,AM,PNT)
      COMMON /BFRE/ EHFR(500),SHFR(500),SEKHF,EHMAX,SHMAX
      DIMENSION AM(500),PNT(3,500)
      I1=IP0+1
      I2=INE+IP0
      IAMAX=0
      EEMAX=0
      DO 1 I=I1,I2
      IA=INT(AM(I)/0.94+0.1)
      P2=PNT(1,I)**2+PNT(2,I)**2+PNT(3,I)**2
      EHFR(IA)=EHFR(IA)+1000.*P2/(2.*AM(I))
      SEKHF=SEKHF+1000.*P2/(2.*AM(I))
      SHFR(IA)=SHFR(IA)+1.
      IF(IA.GT.IAMAX) EAMAX=1000.*P2/(2.*AM(I))
      IF(IA.GT.IAMAX) IAMAX=IA
 1    CONTINUE
      IF(IAMAX.GT.4) EEMAX=EAMAX
      IF(IAMAX.GT.4) EHMAX=EHMAX+EEMAX
      IF(IAMAX.GT.4) SHMAX=SHMAX+IAMAX
      RETURN
      END
      SUBROUTINE BFRENE(IA0,NFRAG)
      COMMON /BFRE/ EHFR(500),SHFR(500),SEKHF,EHMAX,SHMAX
      DO 1 I=1,IA0
      IF(SHFR(I).GT.0.) EHFR(I)=EHFR(I)/(I*SHFR(I))
 1    CONTINUE
      IF(NFRAG.GT.0) SEKHF=SEKHF/(FLOAT(IA0*NFRAG))
      !WRITE 100 4,SEKHF
 4    FORMAT(1X,'sum kin.energ.hot fragm.before Coulomb (MeV/n)',
     *' SEKHF=',F6.2)
      IF(SHMAX.GT.0.) EHMAX=EHMAX/SHMAX
      !WRITE 100 2,EHMAX
 2    FORMAT(1X,' EHMAX(MeV/n)=',F6.2,
     *' kin.energies (MeV/N) of hot fragm. before Coulomb')
      !WRITE 100 3,(EHFR(I),I=1,IA0)
 3    FORMAT(2X,'EHFR=',10F7.2)
      RETURN
      END
c ***
      SUBROUTINE DISNE2(INET,IP0,AMP,PN0,T,TN)
c modernization: instead of generation on Maxwell-Boltzmann distribution
c (DISNET) the generation on dyrect phase space is used (DISIMP). Apr.98, Bo.
c (to avoid the problem of the last particle at Boltzmann generation).
      DIMENSION AMP(500),PN0(3,500),PN1(3,500),AMP1(500)
      IF(INET.LE.0) RETURN
      IF(IP0.LT.0.OR.(IP0+INET).GT.500) RETURN
      DO 1 I1=1,INET
      I0=I1+IP0
      AMP1(I1)=AMP(I0)
      PN1(1,I1)=0.
      PN1(2,I1)=0.
      PN1(3,I1)=0.
 1    CONTINUE
      CALL DISIMP(INET,AMP1,PN1,TN)
      DO 2 I1=1,INET
      I0=I1+IP0
      PN0(1,I0)=PN1(1,I1)
      PN0(2,I0)=PN1(2,I1)
      PN0(3,I0)=PN1(3,I1)
 2    CONTINUE
      RETURN
      END
      SUBROUTINE DISN02(INET,IP0,AMP,PN0,T,TN)
c *** implem.of flow.and rotat. EFLOW, EMOM in GeV/N. ANMOM in h-bar ***
c modernization: instead of generation on Maxwell-Boltzmann distribution
c (DISNET) the generation on dyrect phase space is used (DISIMP). Apr.98, Bo.
c (to avoid the problem of the last particle at Boltzmann generation).
C HAXO[[EH[E [HEP[[[ [ [M[[[[COB INET HE[TPOHOB [O [X TEM[EPAT[PE T(MEV)
C [O[HO[ K[H.[HEP[[[ TN(GEV),C[[TA[ PAC[PE[E[EH[E MAKCBE[OBCK[M. AMP(IP0
C -[X MACC[ (B GEV);[M[[[[C[ (GEV/C),B C[MME [A[[[E 0,[O[AB[[[TC[ B MACC
C PN0(3,500) C (IP0+1).
c -- Dec.99 Ganil, add FLZ - assymetric flow along Z-axis --
      COMMON /BLFLOW/IFLOW,A00,EFLOW /BLFLCM/ALF0,FLOWC,XCM,YCM,ZCM
      COMMON /BLAMOM/IAMOM,ANMOM,EMOM /BLAMCM/BMOM,WR,EMOMC
      COMMON /BLPAT/PA(500),PAZ(500),IPT
      COMMON /BLCEN/XC(500),YC(500),ZC(500)
      COMMON /BLFLZ/FLZ
      DIMENSION AMP(500),PN0(3,500),ANL(3),PX(500),PY(500),PZ(500)
      DIMENSION PN1(3,500),AMP1(500),VV(3)
      IF(INET.LE.0) RETURN
      IF(IP0.LT.0.OR.(IP0+INET).GT.500) RETURN
      IP01=IP0+1
      IP0M=IP0+INET
      IF(INET-1) 1,1,2
 1    P=SQRT(2.*AMP(IP01)*TN)
      CALL ISOTR(ANL)
      VV(1)=0.
      VV(2)=0.
      VV(3)=0.
      IF(IFLOW.GE.1) VV(1)=VV(1)+ALF0*(XC(IP01)-XCM)
      IF(IFLOW.GE.1) VV(2)=VV(2)+ALF0*(YC(IP01)-YCM)
c      IF(IFLOW.GE.1) VV(3)=VV(3)+ALF0*(ZC(IP01)-ZCM)
c -- for nonuniform flow: along Z-axis
      IF(IFLOW.GE.1) VV(3)=VV(3)+FLZ*ALF0*(ZC(IP01)-ZCM)
c -----------------------
      IF(IAMOM.GE.1) VV(3)=VV(3)-WR*(YC(IP01)-YCM)
      IF(IAMOM.GE.1) VV(2)=VV(2)+WR*(ZC(IP01)-ZCM)
      PN0(1,IP01)=P*ANL(1)+AMP(IP01)*VV(1)
      PN0(2,IP01)=P*ANL(2)+AMP(IP01)*VV(2)
      PN0(3,IP01)=P*ANL(3)+AMP(IP01)*VV(3)
      RETURN
 2    IF(INET-2) 3,3,4
 3    P=SQRT(2.*(AMP(IP01)*AMP(IP0M)/(AMP(IP01)+AMP(IP0M)))*TN)
      CALL ISOTR(ANL)
      VV(1)=0.
      VV(2)=0.
      VV(3)=0.
      IF(IFLOW.GE.1) VV(1)=VV(1)+ALF0*(XC(IP01)-XCM)
      IF(IFLOW.GE.1) VV(2)=VV(2)+ALF0*(YC(IP01)-YCM)
c      IF(IFLOW.GE.1) VV(3)=VV(3)+ALF0*(ZC(IP01)-ZCM)
c -- for nonuniform flow: along Z-axis
      IF(IFLOW.GE.1) VV(3)=VV(3)+FLZ*ALF0*(ZC(IP01)-ZCM)
c -----------------------
      IF(IAMOM.GE.1) VV(3)=VV(3)-WR*(YC(IP01)-YCM)
      IF(IAMOM.GE.1) VV(2)=VV(2)+WR*(ZC(IP01)-ZCM)
      PN0(1,IP01)=P*ANL(1)+AMP(IP01)*VV(1)
      PN0(2,IP01)=P*ANL(2)+AMP(IP01)*VV(2)
      PN0(3,IP01)=P*ANL(3)+AMP(IP01)*VV(3)
      VV(1)=0.
      VV(2)=0.
      VV(3)=0.
      IF(IFLOW.GE.1) VV(1)=VV(1)+ALF0*(XC(IP0M)-XCM)
      IF(IFLOW.GE.1) VV(2)=VV(2)+ALF0*(YC(IP0M)-YCM)
c      IF(IFLOW.GE.1) VV(3)=VV(3)+ALF0*(ZC(IP0M)-ZCM)
c -- for nonuniform flow: along Z-axis
      IF(IFLOW.GE.1) VV(3)=VV(3)+FLZ*ALF0*(ZC(IP0M)-ZCM)
c -----------------------
      IF(IAMOM.GE.1) VV(3)=VV(3)-WR*(YC(IP0M)-YCM)
      IF(IAMOM.GE.1) VV(2)=VV(2)+WR*(ZC(IP0M)-ZCM)
      PN0(1,IP0M)=-P*ANL(1)+AMP(IP0M)*VV(1)
      PN0(2,IP0M)=-P*ANL(2)+AMP(IP0M)*VV(2)
      PN0(3,IP0M)=-P*ANL(3)+AMP(IP0M)*VV(3)
      RETURN
 4    CONTINUE
      DO 5 I1=1,INET
      I0=I1+IP0
      AMP1(I1)=AMP(I0)
      PN1(1,I1)=0.
      PN1(2,I1)=0.
      PN1(3,I1)=0.
 5    CONTINUE
      CALL DISIMP(INET,AMP1,PN1,TN)
      DO 6 I1=1,INET
      I0=I1+IP0
      PX(I0)=PN1(1,I1)
      PY(I0)=PN1(2,I1)
      PZ(I0)=PN1(3,I1)
 6    CONTINUE
      DO 8 I=IP01,IP0M
      VV(1)=0.
      VV(2)=0.
      VV(3)=0.
      IF(IFLOW.GE.1) VV(1)=VV(1)+ALF0*(XC(I)-XCM)
      IF(IFLOW.GE.1) VV(2)=VV(2)+ALF0*(YC(I)-YCM)
c      IF(IFLOW.GE.1) VV(3)=VV(3)+ALF0*(ZC(I)-ZCM)
c -- for nonuniform flow: along Z-axis
      IF(IFLOW.GE.1) VV(3)=VV(3)+FLZ*ALF0*(ZC(I)-ZCM)
c -----------------------
      IF(IAMOM.GE.1) VV(3)=VV(3)-WR*(YC(I)-YCM)
      IF(IAMOM.GE.1) VV(2)=VV(2)+WR*(ZC(I)-ZCM)
      PN0(1,I)=PX(I)+AMP(I)*VV(1)
      PN0(2,I)=PY(I)+AMP(I)*VV(2)
      PN0(3,I)=PZ(I)+AMP(I)*VV(3)
 8    CONTINUE
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     FIREBALL MATTER:
C     In the following some routines are added which distribute the
C     produced particles (clusters/free nucleons) in x-space.
C     Different options have been checked/studied with the aim
C     to include the missing dynamics of the BUU runs (in terms
C     radial dependence of expansion) in this analysis
C     as precise as possible. Unfortunately, this still does not work
C     well, since from BUU we do not have information on radial distributions
C     of fragments, but only of single nucleons.
C
C     SPECTATOR MATTER: Here we have no problems, since spectators does not
C     experience violent dynamics.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE PLACE_TG(iflow,IRSYS,ipart,IAMass)


!      COMMON /BLCEN/ XC(500),YC(500),ZC(500) !their positions

      COMMON /BLCEN_THEO/ XC_THEO(500),YC_THEO(500),ZC_THEO(500) !their positions

      integer RSYS
      integer IAMass(500)
      dimension PAA(500)

      dimension ifrm(500),inew(500)

      JJJJ=0
      CEL=1.

!!$c       CELY=2.
!!$c      CELX=0.5
!!$c      IF(IPOSF1.EQ.1) RSYS=1.26*RSYS
!!$C FINDING FRAGMENT POSITIONS IN THE SYSTEM

      IP = ipart

      RN = 1.5
      Factor = 1.5
      if (iflow.eq.0) Factor = 2.0

c      RSYS = 2.0*1.7*(float(IRSYS))**0.33333333
      RSYS = Factor*RN*(float(IRSYS))**0.33333333

      GOTO 9876 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !try to distribute first the heaviest fragments, and then
      !the free nucleons. This seems crucial for a correct
      !sampling of produced particles in phase space.
      !still under study, influences the rapidity spectra in
      !very central HIC, but not absolute yields.

      do ii=1,IP
         ifrm(ii) = 0
      end do
      newlabel = 0

      do 400 ii=1,IP
         if (IAMass(ii).le.1) goto 400
         MaxMass = 0
         itemp   = 0
         do 401 jj=1,IP
            if (IAMass(jj).le.1) goto 401
            if (MaxMass.gt.IAMass(jj) .and. ifrm(jj).eq.0) then
               MaxMass = IAMass(jj)
               itemp   = jj
            endif
 401     continue
         if (itemp.ne.0) then
            newlabel = newlabel + 1
            inew(newlabel) = jj
            ifrm(jj) = 1
         endif
 400  continue

      do ii=1,IP
         if (IAMass(ii).le.1) then
            newlabel = newlabel + 1
            inew(newlabel) = ii
            ifrm(ii) = 1
         endif
      end do

      if (newlabel .ne. IP) then
         write(*,*) 'SMM.f/PLACE_THEO:
     &        wrong re-sorting of fragments!',newlabel
      endif

      do ii=1,IP

         write(*,*) iflow,ii,IAMass(ii),inew(ii),IAMass(inew(ii))

      end do

 9876 continue !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      do ii=1,IP
         PAA(ii) = float( IAMass(ii) )
      end do

      IPOSF = 1

 1    I=1

      xxx = RNDM(-1)

      if(PAA(I).eq.1) then
         RN=0.8
      else
         RN=1.5
      endif
      R=(RSYS-RN*PAA(I)**0.3333333)*xxx**0.3333333

      IF(IPOSF1.EQ.1) R=0.
      CT=1.-2.*RNDM(-1)
      PHI=6.28318*RNDM(-1)
      ST=SQRT(1-CT**2)
      CPHI=COS(PHI)
      SPHI=SIN(PHI)
      XC_THEO(1)=R*ST*CPHI * CEL
      YC_THEO(1)=R*ST*SPHI * CEL
!C      ZC(1)=R*CT
      ZC_THEO(1)=R*CT  !*CEL
      IF(IP.LE.1) RETURN
 3    K=I
      KK=0
      I=I+1
 4    CONTINUE
      KK=KK+1
      IF(KK.GT.1000) GO TO 2
      if(PAA(I).eq.1) then
         RN=0.8
      else
         RN=1.5
      endif
      R=(RSYS-RN*PAA(I)**0.3333333)*(RNDM(-1))**0.3333333
      CT=1.-2.*RNDM(-1)
      PHI=6.28318*RNDM(-1)
      ST=SQRT(1-CT**2)
      CPHI=COS(PHI)
      SPHI=SIN(PHI)
      XC_THEO(I)=R*ST*CPHI * CEL
      YC_THEO(I)=R*ST*SPHI * CEL
!C      ZC(I)=R*CT
      ZC_THEO(I)=R*CT  !*CEL
      DO 5 J=1,K
      RR=(XC_THEO(I)-XC_THEO(J))**2+(YC_THEO(I)-YC_THEO(J))**2+
     &      (ZC_THEO(I)-ZC_THEO(J))**2
      if(PAA(J).eq.1) then
         RN=0.8
      else
         RN=1.5
      endif
      RMIN=RN*PAA(I)**0.3333333+RN*PAA(J)**0.3333333
      if(PAA(J).ne.1) then
         IF(RR-RMIN*RMIN) 4,5,5
      else
         goto 5
      endif
 5    CONTINUE
      IF(I.EQ.IP) RETURN
      GO TO 3
 2    CONTINUE
      JJJJ=JJJJ+1
      IF(JJJJ.GE.1000) PRINT 10 ,KK,I,PAA(I),IP,isource
 10   FORMAT(3X,'ERROR IN PLACE - KK=',I4,' I=',I3,' PAA(I)=',F5.1,
     *' IP=',I3,'ISOURCE=',i3)
      IF(IPOSF1.EQ.1) RSYS=1.5*RSYS
      GO TO 1
      END


c***************************************************************************
      SUBROUTINE PLACE3_TG(iflow,IRSYS,ipart,IAMass)
c***************************************************************************


!      COMMON /BLCEN/ XC(500),YC(500),ZC(500) !their positions

      COMMON /BLCEN_THEO/ XC_THEO(500),YC_THEO(500),ZC_THEO(500) !their positions

      integer RSYS
      integer IAMass(500)
      dimension PAA(500)

      dimension ifrm(500),inew(500)   !!!,iplace(500)

      logical place_flag

      CEL=1.

!!$C FINDING FRAGMENT POSITIONS IN THE SYSTEM

      IP = ipart

      RN = 1.7
      Factor = 1.5
      if (iflow.eq.0) Factor = 2.0

c      RSYS = 2.0*1.7*(float(IRSYS))**0.33333333
      RSYS = Factor*RN*(float(IRSYS))**0.33333333

      !try to distribute first the heaviest fragments, and then
      !the free nucleons. This seems crucial for a correct
      !sampling of produced particles in phase space.

c      GOTO 1234

      IFRAGM = 0
      do ii=1,IP
         ifrm(ii) = 0
         if (IAMass(ii).gt.1) IFRAGM = IFRAGM + 1
      end do

      if (IFRAGM.eq.0) goto 403

      newlabel = 0
 400  continue
      if (newlabel.eq.IFRAGM) goto 402
      MaxMass = 0
      itemp   = 0
      do 401 jj=1,IP
         if (IAMass(jj).le.1) goto 401
         if (IAMass(jj).gt.MaxMass .and. ifrm(jj).eq.0) then
            MaxMass = IAMass(jj)
            itemp   = jj
         endif
 401  continue
      if (itemp.ne.0) then
         newlabel = newlabel + 1
         inew(newlabel) = itemp
         ifrm(itemp) = 1
      endif
      if (newlabel.lt.IFRAGM) goto 400
 402  continue


      do ii=1,IP
         if (IAMass(ii).le.1) then
            newlabel = newlabel + 1
            inew(newlabel) = ii
            ifrm(ii) = 1
         endif
      end do

      if (newlabel .ne. IP) then
         write(*,*) 'SMM.f/PLACE_THEO:
     &        wrong re-sorting of fragments!',newlabel
      endif

 403  continue

c 1234 continue

      do iii=1,IP
         ii = inew(iii)
c         iplace(ii) = 0
         PAA(ii) = float( IAMass(ii) )
      end do

c---------------------------------------------------------------------------

      do 1 ii=1,IP
         i=inew(ii)
c         if (iplace(i).eq.1) goto 1

c         write(*,*) isource,i,IAMass(i),icount

         if (IAMass(i).eq.1) then
            RN = 0.5
         else
            RN = 1.5
         endif

         icount = 0
 4       continue
         icount = icount + 1
         if (icount.gt.1000) then
            write(*,*) 'PLACE_THEO: wrong sampling in x-space!'
            write(*,*) iflow,i,IAMass(i),icount
            write(*,*) XC_THEO(I),YC_THEO(I),ZC_THEO(I),
     &           sqrt(RR),RMIN,dummy
            STOP
         endif
         R=(RSYS-RN*PAA(I)**0.3333333)*RNDM(-1)**0.3333333
c         R=RSYS*RNDM(-1)**0.3333333
         IF(i.EQ.1) R=0.
         CT=1.-2.*RNDM(-1)
         PHI=6.28318*RNDM(-1)
         ST=SQRT(1-CT**2)
         CPHI=COS(PHI)
         SPHI=SIN(PHI)
         XC_THEO(I)=R*ST*CPHI
         YC_THEO(I)=R*ST*SPHI
C     ZC(1)=R*CT
         ZC_THEO(I)=R*CT
c         iplace(i) = 1
         IF(IP.LE.1) goto 1
         if (ii.eq.1) goto 1
         place_flag = .false.
         jj = i
c         do 3 kk=1,jj-1
c            j = inew(kk)
c            if (IAMass(j).eq.1) then
c               RN = 0.5
c            else
c               RN = 1.5
c            endif
c            RR=(XC_THEO(I)-XC_THEO(j))**2+(YC_THEO(I)-YC_THEO(j))**2+
c     &           (ZC_THEO(I)-ZC_THEO(j))**2
c            RMIN=RN*PAA(I)**0.3333333+RN*PAA(j)**0.3333333
c            dummy = RR-RMIN*RMIN
c            IF(dummy.lt.0.001) place_flag = .true.
c 3       continue
c         if (place_flag .eq. .true.) goto 4
 1    continue
c---------------------------------------------------------------------------

c      do ii=1,IP
c         i=inew(ii)
c         RR = sqrt(XC_THEO(i)**2+YC_THEO(i)**2+ZC_THEO(i)**2)
c         write(*,1111) isource,IAMass(i),
c     &        XC_THEO(i),YC_THEO(i),ZC_THEO(i),RR
c      end do

c---------------------------------------------------------------------------

 1111 format('nach sampling: ',2i5,4f9.4)

c***************************************************************************
      END !PLACE3_THEO *****************************************************
c***************************************************************************



      SUBROUTINE PLACE2_TG(isource,IRSYS,ipart,IAMass)


!      COMMON /BLCEN/ XC(500),YC(500),ZC(500) !their positions

      COMMON /BLCEN_THEO/ XC_THEO(500),YC_THEO(500),ZC_THEO(500) !their positions

      integer RSYS
      integer IAMass(500)
      dimension PAA(500)

      dimension ifrm(500),inew(500)

      JJJJ=0
      CEL=1.

!!$c       CELY=2.
!!$c      CELX=0.5
!!$c      IF(IPOSF1.EQ.1) RSYS=1.26*RSYS
!!$C FINDING FRAGMENT POSITIONS IN THE SYSTEM

      IP = ipart

      RN = 1.5
      Factor = 1.
      if (isource.le.2) Factor = 2.0

c      RSYS = 2.0*1.7*(float(IRSYS))**0.33333333
      RSYS = Factor*RN*(float(IRSYS))**0.33333333

      !try to distribute first the heaviest fragments, and then
      !the free nucleons. This seems crucial for a correct
      !sampling of produced particles in phase space.
      !still under study, influences the rapidity spectra in
      !very central HIC, but not absolute yields.

      IFRAGM = 0
      do ii=1,IP
         ifrm(ii) = 0
         if (IAMass(ii).gt.1) IFRAGM = IFRAGM + 1
      end do

      if (IFRAGM.eq.0) goto 403

      newlabel = 0
 400  continue
      if (newlabel.eq.IFRAGM) goto 402
      MaxMass = 1000000
      itemp   = 0
      do 401 jj=1,IP
         if (IAMass(jj).le.1) goto 401
         if (MaxMass.gt.IAMass(jj) .and. ifrm(jj).eq.0) then
            MaxMass = IAMass(jj)
            itemp   = jj
         endif
 401  continue
      if (itemp.ne.0) then
         newlabel = newlabel + 1
         inew(newlabel) = itemp
         ifrm(itemp) = 1
      endif
      if (newlabel.lt.IFRAGM) goto 400
 402  continue


      do ii=1,IP
         if (IAMass(ii).le.1) then
            newlabel = newlabel + 1
            inew(newlabel) = ii
            ifrm(ii) = 1
         endif
      end do

      if (newlabel .ne. IP) then
         write(*,*) 'SMM.f/PLACE_THEO:
     &        wrong re-sorting of fragments!',newlabel
      endif

 403  continue

      do iii=1,IP
         ii = inew(iii)
         PAA(ii) = float( IAMass(ii) )
      end do

      IPOSF = 1

 1    I=1
      IN = inew(I)

      xxx = RNDM(-1)

c      if(PAA(IN).eq.1) then
c         RN=0.8
c      else
c         RN=1.5
c      endif
      R=(RSYS-RN*PAA(IN)**0.3333333)*xxx**0.3333333

      IF(IPOSF1.EQ.1) R=0.
      CT=1.-2.*RNDM(-1)
      PHI=6.28318*RNDM(-1)
      ST=SQRT(1-CT**2)
      CPHI=COS(PHI)
      SPHI=SIN(PHI)
      XC_THEO(IN)=R*ST*CPHI * CEL
      YC_THEO(IN)=R*ST*SPHI * CEL
!C      ZC(1)=R*CT
      ZC_THEO(IN)=R*CT  !*CEL
      IF(IP.LE.1) RETURN
 3    K=I
      KK=0
      I=I+1
      IN = inew(I)
 4    CONTINUE
      KK=KK+1
      IF(KK.GT.1000) GO TO 2
c      if(PAA(IN).eq.1) then
c         RN=0.8
c      else
c         RN=1.5
c      endif
      R=(RSYS-RN*PAA(IN)**0.3333333)*(RNDM(-1))**0.3333333
      CT=1.-2.*RNDM(-1)
      PHI=6.28318*RNDM(-1)
      ST=SQRT(1-CT**2)
      CPHI=COS(PHI)
      SPHI=SIN(PHI)
      XC_THEO(IN)=R*ST*CPHI * CEL
      YC_THEO(IN)=R*ST*SPHI * CEL
!C      ZC(I)=R*CT
      ZC_THEO(IN)=R*CT  !*CEL
      DO 5 J=1,K
         INN = inew(j)
      RR=(XC_THEO(IN)-XC_THEO(INN))**2+(YC_THEO(IN)-YC_THEO(INN))**2+
     &      (ZC_THEO(IN)-ZC_THEO(INN))**2
c      if(PAA(INN).eq.1) then
c         RN=0.8
c      else
c         RN=1.5
c      endif
      RMIN=RN*PAA(IN)**0.3333333+RN*PAA(INN)**0.3333333
c      if(PAA(INN).ne.1) then
      IF(RR-RMIN*RMIN) 4,5,5
c      else
c         goto 5
c      endif
 5    CONTINUE
      IF(I.EQ.IP) RETURN
      GO TO 3
 2    CONTINUE
      JJJJ=JJJJ+1
      IF(JJJJ.GE.1000) PRINT 10 ,KK,I,PAA(I),IP,isource
 10   FORMAT(3X,'ERROR IN PLACE - KK=',I4,' I=',I3,' PAA(I)=',F5.1,
     *' IP=',I3,'ISOURCE=',i3)
      IF(IPOSF1.EQ.1) RSYS=1.5*RSYS
      GO TO 1
      END
