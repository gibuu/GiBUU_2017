C***********************************************************************
C in this version, all occurences of PY... are replaced by RY...
C***********************************************************************

C*********************************************************************
C This is am modified version of the "original" PYTHIA v5.5
C Kai Gallmeister:
C * Since this code has to be compiled with the automatic conversion
C   real -> double, all occurences of /LUJETS/ are replaced by /PYJETS/,
C   the corresponding common block of PYTHIAv6.2
C   Due to this, many routines do not need to know any more, whether
C   FRITIOF or PYTHIA generated the output
C
C*********************************************************************
CCPH  This file has enlarged event record LUJETS size=4000
C*********************************************************************
C*********************************************************************
C*                                                                  **
C*                                                  January 1991    **
C*                                                                  **
C*           The Lund Monte Carlo for Hadronic Processes            **
C*                                                                  **
C*                        PYTHIA version 5.5                        **
C*                                                                  **
C*                        Hans-Uno Bengtsson                        **
C*     Department of Theoretical Physics, University of Uppsala     **
C*            Thunbergsvagen 3, S-752 38 Uppsala, Sweden            **
C*               INTERNET address HANSUNO@THEP.LU.SE                **
C*                       Tel. +18 - 18 32 44                        **
C*                                                                  **
C*                        Torbjorn Sjostrand                        **
C*                    CERN/TH, CH-1211 Geneva 23                    **
C*                BITNET/EARN address TORSJO@CERNVM                 **
C*                       Tel. +22 - 767 28 20                       **
C*                                                                  **
C*       Copyright Hans-Uno Bengtsson and Torbjorn Sjostrand        **
C*                                                                  **
C*********************************************************************
C*********************************************************************
C                                                                    *
C  List of subprograms in order of appearance, with main purpose     *
C  (S = subroutine, F = function, B = block data)                    *
C                                                                    *
C  S   PYINIT   to administer the initialization procedure           *
C  S   PYEVNT   to administer the generation of an event             *
C  S   PYSTAT   to print cross-section and other information         *
C  S   PYINKI   to initialize kinematics of incoming particles       *
C  S   PYINRE   to initialize treatment of resonances                *
C  S   PYXTOT   to give total, elastic and diffractive cross-sect.   *
C  S   PYMAXI   to find differential cross-section maxima            *
C  S   PYPILE   to select multiplicity of pileup events              *
C  S   PYRAND   to select subprocess and kinematics for event        *
C  S   PYSCAT   to set up kinematics and colour flow of event        *
C  S   PYSSPA   to simulate initial state spacelike showers          *
C  S   PYRESD   to perform resonance decays                          *
C  S   PYMULT   to generate multiple interactions                    *
C  S   PYREMN   to add on target remnants                            *
C  S   PYDIFF   to set up kinematics for diffractive events          *
C  S   PYDOCU   to compute cross-sections and handle documentation   *
C  S   PYFRAM   to perform boosts between different frames           *
C  S   PYWIDT   to calculate full and partial widths of resonances   *
C  S   PYOFSH   to calculate partial width into off-shell channels   *
C  S   PYKLIM   to calculate borders of allowed kinematical region   *
C  S   PYKMAP   to construct value of kinematical variable           *
C  S   PYSIGH   to calculate differential cross-sections             *
C  S   PYSTFU   to evaluate structure functions                      *
C  S   PYSTGA   to evaluate photon structure function                *
C  S   PYSPLI   to find flavours left in hadron when one removed     *
C  F   PYGAMM   to evaluate ordinary Gamma function Gamma(x)         *
C  S   PYWAUX   to evaluate auxiliary functions W1(s) and W2(s)      *
C  S   PYI3AU   to evaluate auxiliary function I3(s,t,u,v)           *
C  F   PYSPEN   to evaluate Spence (dilogarithm) function Sp(x)      *
C  S   PYTEST   to test the proper functioning of the package        *
C  B   PYDATA   to contain all default values                        *
C  S   PYKCUT   to provide dummy routine for user kinematical cuts   *
C  S   PYSTFE   to provide interface to Tung or user structure func. *
C                                                                    *
C*********************************************************************

      SUBROUTINE RYINIT(FRAME,BEAM,TARGET,WIN)

C...Initializes the generation procedure; finds maxima of the
C...differential cross-sections to be used for weighting.
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/LUDAT3/MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5)
      COMMON/LUDAT4/CHAF(500)
      CHARACTER CHAF*8
      COMMON/RYSUBS/MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200)
      COMMON/RYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/RYINT1/MINT(400),VINT(400)
      COMMON/RYINT2/ISET(200),KFPR(200,2),COEF(200,20),ICOL(40,4,2)
      COMMON/RYINT5/NGEN(0:200,3),XSEC(0:200,3)
      SAVE /LUDAT1/,/LUDAT2/,/LUDAT3/,/LUDAT4/
      SAVE /RYSUBS/,/RYPARS/,/RYINT1/,/RYINT2/,/RYINT5/
      DIMENSION ALAMIN(20),NFIN(20)
      CHARACTER*(*) FRAME,BEAM,TARGET
      CHARACTER CHFRAM*8,CHBEAM*8,CHTARG*8,CHMO(12)*3,CHLH(2)*6
      DATA ALAMIN/0.20,0.29,0.20,0.40,0.187,0.212,0.191,0.155,
     &0.22,0.16,0.16,0.26,0.36,7*0.2/,NFIN/20*4/
      DATA CHMO/'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep',
     &'Oct','Nov','Dec'/, CHLH/'lepton','hadron'/

C...Reset MINT and VINT arrays. Write headers.
      DO 100 J=1,400
      MINT(J)=0
  100 VINT(J)=0.
      IF(MSTP(127).GE.1) WRITE(MSTU(11),5000) MSTP(181),MSTP(182),
     &MSTP(185),CHMO(MSTP(184)),MSTP(183)
      MSTP(127)=0
      IF(MSTU(12).GE.1) CALL LULIST(0)
      IF(MSTP(122).GE.1) WRITE(MSTU(11),5100)

C...Identify beam and target particles and initialize kinematics.
      CHFRAM=FRAME//' '
      CHBEAM=BEAM//' '
      CHTARG=TARGET//' '
      CALL RYINKI(CHFRAM,CHBEAM,CHTARG,WIN)
      IF(MINT(65).EQ.1) GOTO 160

C...Select partonic subprocesses to be included in the simulation.
      IF(MSEL.NE.0) THEN
        DO 110 I=1,200
  110   MSUB(I)=0
      ENDIF
      IF(MINT(43).EQ.1.AND.(MSEL.EQ.1.OR.MSEL.EQ.2)) THEN
C...Lepton+lepton -> gamma/Z0 or W.
        IF(MINT(11)+MINT(12).EQ.0) MSUB(1)=1
        IF(MINT(11)+MINT(12).NE.0) MSUB(2)=1
      ELSEIF(MINT(43).LE.3.AND.(MSEL.EQ.1.OR.MSEL.EQ.2)) THEN
C...Lepton+hadron: deep inelastic scattering.
        MSUB(11)=1
      ELSEIF(MSEL.EQ.1) THEN
C...High-pT QCD processes:
        MSUB(11)=1
        MSUB(12)=1
        MSUB(13)=1
        MSUB(28)=1
        MSUB(53)=1
        MSUB(68)=1
        IF(MSTP(82).LE.1.AND.CKIN(3).LT.PARP(81)) MSUB(95)=1
        IF(MSTP(82).GE.2.AND.CKIN(3).LT.PARP(82)) MSUB(95)=1
      ELSEIF(MSEL.EQ.2) THEN
C...All QCD processes:
        MSUB(11)=1
        MSUB(12)=1
        MSUB(13)=1
        MSUB(28)=1
        MSUB(53)=1
        MSUB(68)=1
        MSUB(91)=1
        MSUB(92)=1
        MSUB(93)=1
        MSUB(95)=1
      ELSEIF(MSEL.GE.4.AND.MSEL.LE.8) THEN
C...Heavy quark production.
        MSUB(81)=1
        MSUB(82)=1
        MSUB(84)=1
        DO 120 J=1,MIN(8,MDCY(21,3))
  120   MDME(MDCY(21,2)+J-1,1)=0
        MDME(MDCY(21,2)+MSEL-1,1)=1
        MSUB(85)=1
        DO 130 J=1,MIN(8,MDCY(22,3))
  130   MDME(MDCY(22,2)+J-1,1)=0
        MDME(MDCY(22,2)+MSEL-1,1)=1
      ELSEIF(MSEL.EQ.10) THEN
C...Prompt photon production:
        MSUB(14)=1
        MSUB(18)=1
        MSUB(29)=1
      ELSEIF(MSEL.EQ.11) THEN
C...Z0/gamma* production:
        MSUB(1)=1
      ELSEIF(MSEL.EQ.12) THEN
C...W+/- production:
        MSUB(2)=1
      ELSEIF(MSEL.EQ.13) THEN
C...Z0 + jet:
        MSUB(15)=1
        MSUB(30)=1
      ELSEIF(MSEL.EQ.14) THEN
C...W+/- + jet:
        MSUB(16)=1
        MSUB(31)=1
      ELSEIF(MSEL.EQ.15) THEN
C...Z0 & W+/- pair production:
        MSUB(19)=1
        MSUB(20)=1
        MSUB(22)=1
        MSUB(23)=1
        MSUB(25)=1
      ELSEIF(MSEL.EQ.16) THEN
C...H0 production:
        MSUB(3)=1
        MSUB(102)=1
        MSUB(103)=1
        MSUB(123)=1
        MSUB(124)=1
      ELSEIF(MSEL.EQ.17) THEN
C...H0 & Z0 or W+/- pair production:
        MSUB(24)=1
        MSUB(26)=1
      ELSEIF(MSEL.EQ.18) THEN
C...H0 production; interesting processes in e+e-.
        MSUB(24)=1
        MSUB(103)=1
        MSUB(123)=1
        MSUB(124)=1
      ELSEIF(MSEL.EQ.19) THEN
C...H0, H'0 and A0 production; interesting processes in e+e-.
        MSUB(24)=1
        MSUB(103)=1
        MSUB(123)=1
        MSUB(124)=1
        MSUB(153)=1
        MSUB(171)=1
        MSUB(173)=1
        MSUB(174)=1
        MSUB(158)=1
        MSUB(176)=1
        MSUB(178)=1
        MSUB(179)=1
      ELSEIF(MSEL.EQ.21) THEN
C...Z'0 production:
        MSUB(141)=1
      ELSEIF(MSEL.EQ.22) THEN
C...W'+/- production:
        MSUB(142)=1
      ELSEIF(MSEL.EQ.23) THEN
C...H+/- production:
        MSUB(143)=1
      ELSEIF(MSEL.EQ.24) THEN
C...R production:
        MSUB(144)=1
      ELSEIF(MSEL.EQ.25) THEN
C...LQ (leptoquark) production.
        MSUB(145)=1
        MSUB(162)=1
        MSUB(163)=1
        MSUB(164)=1
      ELSEIF(MSEL.GE.35.AND.MSEL.LE.38) THEN
C...Production of one heavy quark (W exchange):
        MSUB(83)=1
        DO 140 J=1,MIN(8,MDCY(21,3))
  140   MDME(MDCY(21,2)+J-1,1)=0
        MDME(MDCY(21,2)+MSEL-31,1)=1
      ENDIF

C...Count number of subprocesses on.
      MINT(48)=0
      DO 150 ISUB=1,200
      IF(MINT(44).LT.4.AND.ISUB.GE.91.AND.ISUB.LE.96.AND.
     &MSUB(ISUB).EQ.1) THEN
        WRITE(MSTU(11),5200) ISUB,CHLH(MINT(41)),CHLH(MINT(42))
        STOP
      ELSEIF(MSUB(ISUB).EQ.1.AND.ISET(ISUB).EQ.-1) THEN
        WRITE(MSTU(11),5300) ISUB
        STOP
      ELSEIF(MSUB(ISUB).EQ.1.AND.ISET(ISUB).LE.-2) THEN
        WRITE(MSTU(11),5400) ISUB
        STOP
      ELSEIF(MSUB(ISUB).EQ.1) THEN
        MINT(48)=MINT(48)+1
      ENDIF
  150 CONTINUE
      IF(MINT(48).EQ.0) THEN
        WRITE(MSTU(11),5500)
        STOP
      ENDIF
      MINT(49)=MINT(48)-MSUB(91)-MSUB(92)-MSUB(93)-MSUB(94)

C...Maximum 4 generations; set maximum number of allowed flavours.
  160 MSTP(1)=MIN(4,MSTP(1))
      MSTU(114)=MIN(MSTU(114),2*MSTP(1))
      MSTP(54)=MIN(MSTP(54),2*MSTP(1))

C...Sum up Cabibbo-Kobayashi-Maskawa factors for each quark/lepton.
      DO 180 I=-20,20
      VINT(180+I)=0.
      IA=IABS(I)
      IF(IA.GE.1.AND.IA.LE.2*MSTP(1)) THEN
        DO 170 J=1,MSTP(1)
        IB=2*J-1+MOD(IA,2)
        IPM=(5-ISIGN(1,I))/2
        IDC=J+MDCY(IA,2)+2
  170   IF(MDME(IDC,1).EQ.1.OR.MDME(IDC,1).EQ.IPM) VINT(180+I)=
     &  VINT(180+I)+VCKM((IA+1)/2,(IB+1)/2)
      ELSEIF(IA.GE.11.AND.IA.LE.10+2*MSTP(1)) THEN
        VINT(180+I)=1.
      ENDIF
  180 CONTINUE

C...Choose Lambda value to use in alpha-strong.
      MSTU(111)=MSTP(2)
      IF(MSTP(3).GE.1) THEN
        ALAM=PARP(1)
        IF(MSTP(51).GE.1.AND.MSTP(51).LE.13) ALAM=ALAMIN(MSTP(51))
        PARP(1)=ALAM
        PARP(61)=ALAM
        PARU(112)=ALAM
        PARJ(81)=ALAM
        IF(MSTP(51).GE.1.AND.MSTP(51).LE.13) MSTU(112)=NFIN(MSTP(51))
      ENDIF

C...Initialize widths and partial widths for resonances.
      CALL RYINRE
      IF(MINT(65).EQ.1) GOTO 200

C...Reset variables for cross-section calculation.
      DO 190 I=0,200
      DO 190 J=1,3
      NGEN(I,J)=0
  190 XSEC(I,J)=0.

C...Find parametrized total cross-sections.
      IF(MINT(44).EQ.4) CALL RYXTOT

C...Maxima of differential cross-sections.
      IF(MSTP(121).LE.1) CALL RYMAXI

C...Initialize possibility of pileup events.
      IF(MSTP(131).NE.0) CALL RYPILE(1)

C...Initialize multiple interactions with variable impact parameter.
      IF(MINT(44).EQ.4.AND.(MINT(49).NE.0.OR.MSTP(131).NE.0).AND.
     &MSTP(82).GE.2) CALL RYMULT(1)
  200 IF(MSTP(122).GE.1) WRITE(MSTU(11),5600)

C...Formats for initialization information.
 5000 FORMAT(///20X,'The Lund Monte Carlo - RYTHIA version ',I1,'.',I1/
     &20X,'**  Last date of change:  ',I2,1X,A3,1X,I4,'  **'/)
 5100 FORMAT('1',18('*'),1X,'RYINIT: initialization of RYTHIA ',
     &'routines',1X,17('*'))
 5200 FORMAT(1X,'Error: process number ',I3,' not meaningful for ',A6,
     &'-',A6,' interactions.'/1X,'Execution stopped!')
 5300 FORMAT(1X,'Error: requested subprocess',I4,' not implemented.'/
     &1X,'Execution stopped!')
 5400 FORMAT(1X,'Error: requested subprocess',I4,' not existing.'/
     &1X,'Execution stopped!')
 5500 FORMAT(1X,'Error: no subprocess switched on.'/
     &1X,'Execution stopped.')
 5600 FORMAT(/1X,22('*'),1X,'RYINIT: initialization completed',1X,
     &22('*'))

      RETURN
      END

C*********************************************************************
      SUBROUTINE RYEVNT
      PARAMETER (KSZJ=4000)

C...Administers the generation of a high-pT event via calls to
C...a number of subroutines.
      
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      integer N,NPAD,K
      double precision P,V
      SAVE /PYJETS/


      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/RYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/RYINT1/MINT(400),VINT(400)
      COMMON/RYINT2/ISET(200),KFPR(200,2),COEF(200,20),ICOL(40,4,2)
      COMMON/RYINT5/NGEN(0:200,3),XSEC(0:200,3)
      SAVE /LUDAT1/,/LUDAT2/
      SAVE /RYPARS/,/RYINT1/,/RYINT2/,/RYINT5/

C...Initial values for some counters.
      MINT(5)=MINT(5)+1
      MINT(7)=0
      MINT(8)=0
      MINT(83)=0
      MINT(84)=MSTP(126)
      MSTU(24)=0
      MSTU70=0
      MSTJ14=MSTJ(14)

C...Loop over number of pileup events; check space left.
      IF(MSTP(131).LE.0) THEN
        NPILE=1
      ELSE
        CALL RYPILE(2)
        NPILE=MINT(81)
      ENDIF

      DO 170 IPILE=1,NPILE
      IF(MINT(84)+100.GE.MSTU(4)) THEN
        CALL LUERRM(11,
     &  '(RYEVNT:) no more space in LUJETS for pileup events')
        IF(MSTU(21).GE.1) GOTO 180
      ENDIF
      MINT(82)=IPILE

C...Generate variables of hard scattering.
  100 CONTINUE
c......jochen
      MINT(31)=0
      MINT(51)=0

      CALL RYRAND

      ISUB=MINT(1)
      IF(MSTP(111).EQ.-1) GOTO 160

      IF(ISUB.LE.90.OR.ISUB.GE.95) THEN
C...Hard scattering (including low-pT):
C...reconstruct kinematics and colour flow of hard scattering.

        CALL RYSCAT

        IF(MINT(51).EQ.1) GOTO 100

        IPU1=MINT(84)+1
        IPU2=MINT(84)+2
        IF(ISUB.EQ.95) GOTO 110

C...Showering of initial state partons (optional).
        IF(MSTP(61).GE.1.AND.MINT(44).GE.2) CALL RYSSPA(IPU1,IPU2)

C...Showering of final state partons (optional).
        IF(MSTP(71).GE.1.AND.ISET(ISUB).GE.2) THEN
          IPU3=MINT(84)+3
          IPU4=MINT(84)+4
          IF(ISET(ISUB).EQ.6) IPU4=-3
          QMAX=SQRT(PARP(71)*VINT(52))
          IF(ISET(ISUB).GE.3.AND.ISET(ISUB).LE.5) QMAX=PMAS(23,1)
          IF(ISUB.EQ.8.OR.ISUB.EQ.76.OR.ISUB.EQ.77.OR.ISUB.EQ.124.OR.
     &    ISUB.EQ.174.OR.ISUB.EQ.179) QMAX=PMAS(24,1)
          CALL LUSHOW(IPU3,IPU4,QMAX)
        ENDIF

C...Decay of final state resonances.
        IF(MSTP(41).GE.1) CALL RYRESD
        IF(MINT(51).EQ.1) GOTO 100
        MINT(52)=N

C...Multiple interactions.
        IF(MSTP(81).GE.1.AND.MINT(44).EQ.4) CALL RYMULT(6)
        MINT(53)=N

C...Hadron remnants and primordial kT.
  110   CALL RYREMN(IPU1,IPU2)
        IF(MINT(51).EQ.1) GOTO 100

      ELSE
C...Diffractive and elastic scattering.
        CALL RYDIFF
      ENDIF

C...Recalculate energies from momenta and masses (if desired).
      IF(MSTP(113).GE.1) THEN
        DO 120 I=MINT(83)+1,N
  120   IF(K(I,1).GT.0.AND.K(I,1).LE.10) P(I,4)=SQRT(P(I,1)**2+
     &  P(I,2)**2+P(I,3)**2+P(I,5)**2)
      ENDIF

C...Rearrange partons along strings, check invariant mass cuts.
      MSTU(28)=0
      IF(MSTP(111).LE.0) MSTJ(14)=-1
      CALL LUPREP(MINT(84)+1)
      MSTJ(14)=MSTJ14
      IF(MSTP(112).EQ.1.AND.MSTU(28).EQ.3) GOTO 100
      IF(MSTP(125).EQ.0.OR.MSTP(125).EQ.1) THEN
        DO 130 I=MINT(84)+1,N
        IF(K(I,2).NE.94) GOTO 130
        K(I+1,3)=MOD(K(I+1,4)/MSTU(5),MSTU(5))
        K(I+2,3)=MOD(K(I+2,4)/MSTU(5),MSTU(5))
  130   CONTINUE
        CALL LUEDIT(12)
        CALL LUEDIT(14)
        IF(MSTP(125).EQ.0) CALL LUEDIT(15)
        IF(MSTP(125).EQ.0) MINT(4)=0
        DO 150 I=MINT(83)+1,N
        IF(K(I,1).EQ.11.AND.K(I,4).EQ.0.AND.K(I,5).EQ.0) THEN
          DO 140 I1=I+1,N
          IF(K(I1,3).EQ.I.AND.K(I,4).EQ.0) K(I,4)=I1
  140     IF(K(I1,3).EQ.I) K(I,5)=I1
        ENDIF
  150   CONTINUE
      ENDIF

C...Introduce separators between sections in LULIST event listing.
      IF(IPILE.EQ.1.AND.MSTP(125).LE.0) THEN
        MSTU70=1
        MSTU(71)=N
      ELSEIF(IPILE.EQ.1) THEN
        MSTU70=3
        MSTU(71)=2
        MSTU(72)=MINT(4)
        MSTU(73)=N
      ENDIF

C...Perform hadronization (if desired).
      IF(MSTP(111).GE.1) THEN
        CALL LUEXEC
        IF(MSTU(24).NE.0) GOTO 100
      ENDIF
      IF(MSTP(125).EQ.0.OR.MSTP(125).EQ.1) CALL LUEDIT(14)

C...Store event information and calculate Monte Carlo estimates of
C...subprocess cross-sections.
  160 IF(IPILE.EQ.1) CALL RYDOCU

C...Set counters for current pileup event and loop to next one.
      MSTI(41)=IPILE
      IF(IPILE.GE.2.AND.IPILE.LE.10) MSTI(40+IPILE)=ISUB
      IF(MSTU70.LT.10) THEN
        MSTU70=MSTU70+1
        MSTU(70+MSTU70)=N
      ENDIF
      MINT(83)=N
      MINT(84)=N+MSTP(126)
  170 CONTINUE

C...Generic information on pileup events.
      IF(MSTP(131).EQ.1.AND.MSTP(133).GE.1) THEN
        PARI(91)=VINT(132)
        PARI(92)=VINT(133)
        PARI(93)=VINT(134)
        IF(MSTP(133).GE.2) PARI(93)=PARI(93)*XSEC(0,3)/VINT(131)
      ENDIF

C...Transform to the desired coordinate frame.
  180 CALL RYFRAM(MSTP(124))
      MSTU(70)=MSTU70
      PARU(21)=VINT(1)

      RETURN
      END

C***********************************************************************

      SUBROUTINE RYSTAT(MSTAT)

C...Prints out information about cross-sections, decay widths, branching
C...ratios, kinematical limits, status codes and parameter values.
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/LUDAT3/MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5)
      COMMON/RYSUBS/MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200)
      COMMON/RYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/RYINT1/MINT(400),VINT(400)
      COMMON/RYINT4/WIDP(21:40,0:40),WIDE(21:40,0:40),WIDS(21:40,3)
      COMMON/RYINT5/NGEN(0:200,3),XSEC(0:200,3)
      COMMON/RYINT6/PROC(0:200)
      CHARACTER PROC*28
      SAVE /LUDAT1/,/LUDAT2/,/LUDAT3/
      SAVE /RYSUBS/,/RYPARS/,/RYINT1/,/RYINT4/,/RYINT5/,/RYINT6/
      CHARACTER CHAU*16,CHPA(-40:40)*9,CHIN(2)*12,
     &STATE(-1:5)*4,CHKIN(21)*18
      DATA STATE/'----','off ','on  ','on/+','on/-','on/1','on/2'/,
     &CHKIN/' m_hard (GeV/c^2) ',' p_T_hard (GeV/c) ',
     &'m_finite (GeV/c^2)','   y*_subsystem   ','     y*_large     ',
     &'     y*_small     ','    eta*_large    ','    eta*_small    ',
     &'cos(theta*)_large ','cos(theta*)_small ','       x_1        ',
     &'       x_2        ','       x_F        ',' cos(theta_hard)  ',
     &'m''_hard (GeV/c^2) ','       tau        ','        y*        ',
     &'cos(theta_hard^-) ','cos(theta_hard^+) ','      x_T^2       ',
     &'       tau''       '/

C...Cross-sections.
      IF(MSTAT.LE.1) THEN
        WRITE(MSTU(11),5000)
        WRITE(MSTU(11),5100)
        WRITE(MSTU(11),5200) 0,PROC(0),NGEN(0,3),NGEN(0,1),XSEC(0,3)
        DO 100 I=1,200
        IF(MSUB(I).NE.1) GOTO 100
        WRITE(MSTU(11),5200) I,PROC(I),NGEN(I,3),NGEN(I,1),XSEC(I,3)
  100   CONTINUE
        WRITE(MSTU(11),5300) 1.-FLOAT(NGEN(0,3))/
     &  MAX(1.,FLOAT(NGEN(0,2)))

C...Decay widths and branching ratios.
      ELSEIF(MSTAT.EQ.2) THEN
        DO 110 KF=-40,40
        CALL LUNAME(KF,CHAU)
  110   CHPA(KF)=CHAU(1:9)
        WRITE(MSTU(11),5400)
        WRITE(MSTU(11),5500)
C...Off-shell branchings.
        DO 130 I=1,17
        KC=I
        IF(I.GE.9) KC=I+2
        IF(I.EQ.17) KC=21
        NGP=0
        IF(KC.LE.20) NGP=(MOD(KC,10)+1)/2
        IF(NGP.LE.MSTP(1)) WRITE(MSTU(11),5600) KC,CHPA(KC),PMAS(KC,1),
     &  0.,0.,STATE(MDCY(KC,1)),0.
        DO 120 J=1,MDCY(KC,3)
        IDC=J+MDCY(KC,2)-1
        NGP1=0
        IF(IABS(KFDP(IDC,1)).LE.20) NGP1=
     &  (MOD(IABS(KFDP(IDC,1)),10)+1)/2
        NGP2=0
        IF(IABS(KFDP(IDC,2)).LE.20) NGP2=
     &  (MOD(IABS(KFDP(IDC,2)),10)+1)/2
  120   IF(MDME(IDC,2).EQ.102.AND.NGP1.LE.MSTP(1).AND.NGP2.LE.MSTP(1))
     &  WRITE(MSTU(11),5700) IDC,CHPA(KFDP(IDC,1)),CHPA(KFDP(IDC,2)),
     &  0.,0.,STATE(MDME(IDC,1)),0.
  130   CONTINUE
C...On-shell decays.
        DO 150 I=1,10
        KC=I+22
        IF(I.EQ.4) KC=32
        IF(I.GE.5.AND.I.LE.8) KC=I+29
        IF(I.GE.9) KC=I+30
        IF(WIDE(KC,0).GT.0.) THEN
          WRITE(MSTU(11),5600) KC,CHPA(KC),PMAS(KC,1),WIDP(KC,0),1.,
     &    STATE(MDCY(KC,1)),1.
          DO 140 J=1,MDCY(KC,3)
          IDC=J+MDCY(KC,2)-1
          NGP1=0
          IF(IABS(KFDP(IDC,1)).LE.20) NGP1=
     &    (MOD(IABS(KFDP(IDC,1)),10)+1)/2
          NGP2=0
          IF(IABS(KFDP(IDC,2)).LE.20) NGP2=
     &    (MOD(IABS(KFDP(IDC,2)),10)+1)/2
  140     IF(NGP1.LE.MSTP(1).AND.NGP2.LE.MSTP(1)) WRITE(MSTU(11),5700)
     &    IDC,CHPA(KFDP(IDC,1)),CHPA(KFDP(IDC,2)),WIDP(KC,J),
     &    WIDP(KC,J)/WIDP(KC,0),STATE(MDME(IDC,1)),
     &    WIDE(KC,J)/WIDE(KC,0)
        ELSE
          WRITE(MSTU(11),5600) KC,CHPA(KC),PMAS(KC,1),WIDP(KC,0),1.,
     &    STATE(MDCY(KC,1)),0.
        ENDIF
  150   CONTINUE
        WRITE(MSTU(11),5800)

C...Allowed incoming partons/particles at hard interaction.
      ELSEIF(MSTAT.EQ.3) THEN
        WRITE(MSTU(11),5900)
        CALL LUNAME(MINT(11),CHAU)
        CHIN(1)=CHAU(1:12)
        CALL LUNAME(MINT(12),CHAU)
        CHIN(2)=CHAU(1:12)
        WRITE(MSTU(11),6000) CHIN(1),CHIN(2)
        DO 160 KF=-40,40
        CALL LUNAME(KF,CHAU)
  160   CHPA(KF)=CHAU(1:9)
        DO 170 I=-20,22
        IF(I.EQ.0) GOTO 170
        IA=IABS(I)
        IF(IA.GT.MSTP(54).AND.IA.LE.10) GOTO 170
        IF(IA.GT.10+2*MSTP(1).AND.IA.LE.20) GOTO 170
        WRITE(MSTU(11),6100) CHPA(I),STATE(KFIN(1,I)),CHPA(I),
     &  STATE(KFIN(2,I))
  170   CONTINUE
        WRITE(MSTU(11),6200)

C...User-defined limits on kinematical variables.
      ELSEIF(MSTAT.EQ.4) THEN
        WRITE(MSTU(11),6300)
        WRITE(MSTU(11),6400)
        SHRMAX=CKIN(2)
        IF(SHRMAX.LT.0.) SHRMAX=VINT(1)
        WRITE(MSTU(11),6500) CKIN(1),CHKIN(1),SHRMAX
        PTHMIN=MAX(CKIN(3),CKIN(5))
        PTHMAX=CKIN(4)
        IF(PTHMAX.LT.0.) PTHMAX=0.5*SHRMAX
        WRITE(MSTU(11),6600) CKIN(3),PTHMIN,CHKIN(2),PTHMAX
        WRITE(MSTU(11),6700) CHKIN(3),CKIN(6)
        DO 180 I=4,14
  180   WRITE(MSTU(11),6500) CKIN(2*I-1),CHKIN(I),CKIN(2*I)
        SPRMAX=CKIN(32)
        IF(SPRMAX.LT.0.) SPRMAX=VINT(1)
        WRITE(MSTU(11),6500) CKIN(31),CHKIN(15),SPRMAX
        WRITE(MSTU(11),6800)

C...Status codes and parameter values.
      ELSEIF(MSTAT.EQ.5) THEN
        WRITE(MSTU(11),6900)
        WRITE(MSTU(11),7000)
        DO 190 I=1,100
  190   WRITE(MSTU(11),7100) I,MSTP(I),PARP(I),100+I,MSTP(100+I),
     &  PARP(100+I)
      ENDIF

C...Formats for printouts.
 5000 FORMAT('1',9('*'),1X,'RYSTAT:  Statistics on Number of ',
     &'Events and Cross-sections',1X,9('*'))
 5100 FORMAT(/1X,78('=')/1X,'I',34X,'I',28X,'I',12X,'I'/1X,'I',12X,
     &'Subprocess',12X,'I',6X,'Number of points',6X,'I',4X,'Sigma',3X,
     &'I'/1X,'I',34X,'I',28X,'I',12X,'I'/1X,'I',34('-'),'I',28('-'),
     &'I',4X,'(mb)',4X,'I'/1X,'I',34X,'I',28X,'I',12X,'I'/1X,'I',1X,
     &'N:o',1X,'Type',25X,'I',4X,'Generated',9X,'Tried',1X,'I',12X,
     &'I'/1X,'I',34X,'I',28X,'I',12X,'I'/1X,78('=')/1X,'I',34X,'I',28X,
     &'I',12X,'I')
 5200 FORMAT(1X,'I',1X,I3,1X,A28,1X,'I',1X,I12,1X,I13,1X,'I',1X,1P,
     &E10.3,1X,'I')
 5300 FORMAT(1X,'I',34X,'I',28X,'I',12X,'I'/1X,78('=')//
     &1X,'********* Fraction of events that fail fragmentation ',
     &'cuts =',1X,F8.5,' *********'/)
 5400 FORMAT('1',17('*'),1X,'RYSTAT:  Decay Widths and Branching ',
     &'Ratios',1X,17('*'))
 5500 FORMAT(/1X,78('=')/1X,'I',29X,'I',13X,'I',12X,'I',6X,'I',12X,'I'/
     &1X,'I',1X,'Branching/Decay Channel',5X,'I',1X,'Width (GeV)',1X,
     &'I',7X,'B.R.',1X,'I',1X,'Stat',1X,'I',2X,'Eff. B.R.',1X,'I'/1X,
     &'I',29X,'I',13X,'I',12X,'I',6X,'I',12X,'I'/1X,78('='))
 5600 FORMAT(1X,'I',29X,'I',13X,'I',12X,'I',6X,'I',12X,'I'/1X,'I',1X,
     &I4,1X,A9,'(',1P,E8.2,0P,')',1X,'->',1X,'I',2X,1P,E10.3,0P,1X,
     &'I',1X,1P,E10.3,0P,1X,'I',1X,A4,1X,'I',1X,1P,E10.3,0P,1X,'I')
 5700 FORMAT(1X,'I',1X,I4,1X,A9,1X,'+',1X,A9,2X,'I',2X,1P,E10.3,0P,
     &1X,'I',1X,1P,E10.3,0P,1X,'I',1X,A4,1X,'I',1X,1P,E10.3,0P,1X,'I')
 5800 FORMAT(1X,'I',29X,'I',13X,'I',12X,'I',6X,'I',12X,'I'/1X,78('='))
 5900 FORMAT('1',7('*'),1X,'RYSTAT: Allowed Incoming Partons/',
     &'Particles at Hard Interaction',1X,7('*'))
 6000 FORMAT(/1X,78('=')/1X,'I',38X,'I',37X,'I'/1X,'I',1X,
     &'Beam particle:',1X,A12,10X,'I',1X,'Target particle:',1X,A12,7X,
     &'I'/1X,'I',38X,'I',37X,'I'/1X,'I',1X,'Content',6X,'State',19X,
     &'I',1X,'Content',6X,'State',18X,'I'/1X,'I',38X,'I',37X,'I'/1X,
     &78('=')/1X,'I',38X,'I',37X,'I')
 6100 FORMAT(1X,'I',1X,A9,5X,A4,19X,'I',1X,A9,5X,A4,18X,'I')
 6200 FORMAT(1X,'I',38X,'I',37X,'I'/1X,78('='))
 6300 FORMAT('1',12('*'),1X,'RYSTAT: User-Defined Limits on ',
     &'Kinematical Variables',1X,12('*'))
 6400 FORMAT(/1X,78('=')/1X,'I',76X,'I')
 6500 FORMAT(1X,'I',16X,1P,E10.3,0P,1X,'<',1X,A,1X,'<',1X,1P,E10.3,0P,
     &16X,'I')
 6600 FORMAT(1X,'I',3X,1P,E10.3,0P,1X,'(',1P,E10.3,0P,')',1X,'<',1X,A,
     &1X,'<',1X,1P,E10.3,0P,16X,'I')
 6700 FORMAT(1X,'I',29X,A,1X,'=',1X,1P,E10.3,0P,16X,'I')
 6800 FORMAT(1X,'I',76X,'I'/1X,78('='))
 6900 FORMAT('1',12('*'),1X,'RYSTAT: Summary of Status Codes and ',
     &'Parameter Values',1X,12('*'))
 7000 FORMAT(/3X,'I',4X,'MSTP(I)',9X,'PARP(I)',20X,'I',4X,'MSTP(I)',9X,
     &'PARP(I)'/)
 7100 FORMAT(1X,I3,5X,I6,6X,1P,E10.3,0P,18X,I3,5X,I6,6X,1P,E10.3)

      RETURN
      END

C*********************************************************************

      SUBROUTINE RYINKI(CHFRAM,CHBEAM,CHTARG,WIN)
      PARAMETER (KSZJ=4000)

C...Identifies the two incoming particles and sets up kinematics,
C...including rotations and boosts to/from CM frame.
      
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      integer N,NPAD,K
      double precision P,V
      SAVE /PYJETS/


      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/RYSUBS/MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200)
      COMMON/RYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/RYINT1/MINT(400),VINT(400)
      SAVE /LUDAT1/,/LUDAT2/
      SAVE /RYSUBS/,/RYPARS/,/RYINT1/
      CHARACTER CHFRAM*8,CHBEAM*8,CHTARG*8,CHCOM(3)*8,CHALP(2)*26,
     &CHIDNT(3)*8,CHTEMP*8,CHCDE(19)*8,CHINIT*76
      DIMENSION LEN(3),KCDE(19)
      DATA CHALP/'abcdefghijklmnopqrstuvwxyz',
     &'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
      DATA CHCDE/'e-      ','e+      ','nu_e    ','nu_e~   ',
     &'mu-     ','mu+     ','nu_mu   ','nu_mu~  ','tau-    ',
     &'tau+    ','nu_tau  ','nu_tau~ ','pi+     ','pi-     ',
     &'n0      ','n~0     ','p+      ','p~-     ','gamma   '/
      DATA KCDE/11,-11,12,-12,13,-13,14,-14,15,-15,16,-16,
     &211,-211,2112,-2112,2212,-2212,22/

C...Convert character variables to lowercase and find their length.
      CHCOM(1)=CHFRAM
      CHCOM(2)=CHBEAM
      CHCOM(3)=CHTARG
      DO 120 I=1,3
      LEN(I)=8
      DO 100 LL=8,1,-1
      IF(LEN(I).EQ.LL.AND.CHCOM(I)(LL:LL).EQ.' ') LEN(I)=LL-1
      DO 100 LA=1,26
  100 IF(CHCOM(I)(LL:LL).EQ.CHALP(2)(LA:LA)) CHCOM(I)(LL:LL)=
     &CHALP(1)(LA:LA)
      CHIDNT(I)=CHCOM(I)

C...Fix up bar, underscore and charge in particle name (if needed).
      DO 110 LL=1,6
      IF(CHIDNT(I)(LL:LL+2).EQ.'bar') THEN
        CHTEMP=CHIDNT(I)
        CHIDNT(I)=CHTEMP(1:LL-1)//'~'//CHTEMP(LL+3:8)//'  '
      ENDIF
  110 CONTINUE
      IF(CHIDNT(I)(1:2).EQ.'nu'.AND.CHIDNT(I)(3:3).NE.'_') THEN
        CHTEMP=CHIDNT(I)
        CHIDNT(I)='nu_'//CHTEMP(3:7)
      ELSEIF(CHIDNT(I)(1:2).EQ.'n ') THEN
        CHIDNT(I)(1:3)='n0 '
      ELSEIF(CHIDNT(I)(1:2).EQ.'n~') THEN
        CHIDNT(I)(1:3)='n~0'
      ELSEIF(CHIDNT(I)(1:2).EQ.'p ') THEN
        CHIDNT(I)(1:3)='p+ '
      ELSEIF(CHIDNT(I)(1:2).EQ.'p~'.OR.CHIDNT(I)(1:2).EQ.'p-') THEN
        CHIDNT(I)(1:3)='p~-'
      ENDIF
  120 CONTINUE

C...Identify free initialization.
      IF(CHCOM(1)(1:2).EQ.'no') THEN
        MINT(65)=1
        RETURN
      ENDIF

C...Set initial state. Error for unknown codes. Reset variables.
      N=2
      DO 140 I=1,2
      K(I,1)=1
      K(I,2)=0
      DO 130 J=1,19
  130 IF(CHIDNT(I+1).EQ.CHCDE(J)) K(I,2)=KCDE(J)
cCG      P(I,5)=ULMASS(K(I,2))
      P(I,5)=ULGAMASS(I,0)
      MINT(40+I)=1
      IF(IABS(K(I,2)).GT.100.OR.K(I,2).EQ.22) MINT(40+I)=2
      MINT(44+I)=MINT(40+I)
      IF(MSTP(11).GE.1.AND.IABS(K(I,2)).EQ.11) MINT(44+I)=3
      DO 140 J=1,5
  140 V(I,J)=0.
      IF(K(1,2).EQ.0) WRITE(MSTU(11),5000) CHBEAM(1:LEN(2))
      IF(K(2,2).EQ.0) WRITE(MSTU(11),5100) CHTARG(1:LEN(3))
      IF(K(1,2).EQ.0.OR.K(2,2).EQ.0) STOP
      DO 150 J=6,10
  150 VINT(J)=0.
      CHINIT=' '

C...Set up kinematics for events defined in CM frame.
      IF(CHCOM(1)(1:2).EQ.'cm') THEN
        IF(CHCOM(2)(1:1).NE.'e') THEN
          LOFFS=(31-(LEN(2)+LEN(3)))/2
          CHINIT(LOFFS+1:76)='RYTHIA will be initialized for a '//
     &    CHCOM(2)(1:LEN(2))//' on '//CHCOM(3)(1:LEN(3))//' collider'//
     &    ' '
        ELSE
          LOFFS=(30-(LEN(2)+LEN(3)))/2
          CHINIT(LOFFS+1:76)='RYTHIA will be initialized for an '//
     &    CHCOM(2)(1:LEN(2))//' on '//CHCOM(3)(1:LEN(3))//' collider'//
     &    ' '
        ENDIF
        IF(MSTP(122).GE.1) WRITE(MSTU(11),5200) CHINIT
        IF(MSTP(122).GE.1) WRITE(MSTU(11),5300) WIN
        S=WIN**2
        P(1,1)=0.
        P(1,2)=0.
        P(2,1)=0.
        P(2,2)=0.
        P(1,3)=SQRT(((S-P(1,5)**2-P(2,5)**2)**2-(2.*P(1,5)*P(2,5))**2)/
     &  (4.*S))
        P(2,3)=-P(1,3)
        P(1,4)=SQRT(P(1,3)**2+P(1,5)**2)
        P(2,4)=SQRT(P(2,3)**2+P(2,5)**2)

C...Set up kinematics for fixed target events.
      ELSEIF(CHCOM(1)(1:3).EQ.'fix') THEN
        LOFFS=(29-(LEN(2)+LEN(3)))/2
        CHINIT(LOFFS+1:76)='RYTHIA will be initialized for '//
     &  CHCOM(2)(1:LEN(2))//' on '//CHCOM(3)(1:LEN(3))//
     &  ' fixed target'//' '
        IF(MSTP(122).GE.1) WRITE(MSTU(11),5200) CHINIT
        IF(MSTP(122).GE.1) WRITE(MSTU(11),5400) WIN
        P(1,1)=0.
        P(1,2)=0.
        P(2,1)=0.
        P(2,2)=0.
        P(1,3)=WIN
        P(1,4)=SQRT(P(1,3)**2+P(1,5)**2)
        P(2,3)=0.
        P(2,4)=P(2,5)
        S=P(1,5)**2+P(2,5)**2+2.*P(2,4)*P(1,4)
        VINT(10)=P(1,3)/(P(1,4)+P(2,4))
        CALL LUROBO(0.,0.,0.,0.,-VINT(10))
        IF(MSTP(122).GE.1) WRITE(MSTU(11),5500) SQRT(S)

C...Set up kinematics for events in user-defined frame.
      ELSEIF(CHCOM(1)(1:3).EQ.'use') THEN
        LOFFS=(12-(LEN(2)+LEN(3)))/2
        CHINIT(LOFFS+1:76)='RYTHIA will be initialized for '//
     &  CHCOM(2)(1:LEN(2))//' on '//CHCOM(3)(1:LEN(3))//
     &  ' user-specified configuration'//' '
        IF(MSTP(122).GE.1) WRITE(MSTU(11),5200) CHINIT
        IF(MSTP(122).GE.1) WRITE(MSTU(11),5600)
        IF(MSTP(122).GE.1) WRITE(MSTU(11),5700) CHCOM(2),P(1,1),
     &  P(1,2),P(1,3)
        IF(MSTP(122).GE.1) WRITE(MSTU(11),5700) CHCOM(3),P(2,1),
     &  P(2,2),P(2,3)
        P(1,4)=SQRT(P(1,1)**2+P(1,2)**2+P(1,3)**2+P(1,5)**2)
        P(2,4)=SQRT(P(2,1)**2+P(2,2)**2+P(2,3)**2+P(2,5)**2)
        DO 160 J=1,3
  160   VINT(7+J)=(DBLE(P(1,J))+DBLE(P(2,J)))/DBLE(P(1,4)+P(2,4))
        CALL LUROBO(0.,0.,-VINT(8),-VINT(9),-VINT(10))
        VINT(7)=ULANGL(P(1,1),P(1,2))
        CALL LUROBO(0.,-VINT(7),0.,0.,0.)
        VINT(6)=ULANGL(P(1,3),P(1,1))
        CALL LUROBO(-VINT(6),0.,0.,0.,0.)
        S=P(1,5)**2+P(2,5)**2+2.*(P(1,4)*P(2,4)-P(1,3)*P(2,3))
        IF(MSTP(122).GE.1) WRITE(MSTU(11),5500) SQRT(S)

C...Unknown frame. Error for too low CM energy.
      ELSE
        WRITE(MSTU(11),5800) CHFRAM(1:LEN(1))
        STOP
      ENDIF
      IF(S.LT.PARP(2)**2) THEN
        WRITE(MSTU(11),5900) SQRT(S)
        STOP
      ENDIF

C...Save information on incoming particles.
      MINT(11)=K(1,2)
      MINT(12)=K(2,2)
      MINT(43)=2*MINT(41)+MINT(42)-2
      MINT(44)=MINT(43)
      IF(MINT(11).EQ.22) MINT(44)=MINT(44)-2
      IF(MINT(12).EQ.22) MINT(44)=MINT(44)-1
      MINT(47)=2*MIN(2,MINT(45))+MIN(2,MINT(46))-2
      IF(MIN(MINT(45),MINT(46)).EQ.3) MINT(47)=5
      VINT(1)=SQRT(S)
      VINT(2)=S
      VINT(3)=P(1,5)
      VINT(4)=P(2,5)
      VINT(5)=P(1,3)

C...Store constants to be used in generation.
      IF(MSTP(82).LE.1) VINT(149)=4.*PARP(81)**2/S
      IF(MSTP(82).GE.2) VINT(149)=4.*PARP(82)**2/S

C...Formats for initialization and error information.
 5000 FORMAT(1X,'Error: unrecognized beam particle ''',A,'''.'/
     &1X,'Execution stopped!')
 5100 FORMAT(1X,'Error: unrecognized target particle ''',A,'''.'/
     &1X,'Execution stopped!')
 5200 FORMAT(/1X,78('=')/1X,'I',76X,'I'/1X,'I',A76,'I')
 5300 FORMAT(1X,'I',18X,'at',1X,F10.3,1X,'GeV center-of-mass energy',
     &19X,'I'/1X,'I',76X,'I'/1X,78('='))
 5400 FORMAT(1X,'I',22X,'at',1X,F10.3,1X,'GeV/c lab-momentum',22X,'I')
 5500 FORMAT(1X,'I',76X,'I'/1X,'I',11X,'corresponding to',1X,F10.3,1X,
     &'GeV center-of-mass energy',12X,'I'/1X,'I',76X,'I'/1X,78('='))
 5600 FORMAT(1X,'I',76X,'I'/1X,'I',25X,'px (GeV/c)',3X,'py (GeV/c)',3X,
     &'pz (GeV/c)',15X,'I')
 5700 FORMAT(1X,'I',15X,A8,3(2X,F10.3,1X),14X,'I')
 5800 FORMAT(1X,'Error: unrecognized coordinate frame ''',A,'''.'/
     &1X,'Execution stopped!')
 5900 FORMAT(1X,'Error: too low CM energy,',F8.3,' GeV for event ',
     &'generation.'/1X,'Execution stopped!')

      RETURN
      END

C*********************************************************************

      SUBROUTINE RYINRE

C...Calculates full and effective widths of gauge bosons, stores
C...masses and widths, rescales coefficients to be used for
C...resonance production generation.
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/LUDAT3/MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5)
      COMMON/LUDAT4/CHAF(500)
      CHARACTER CHAF*8
      COMMON/RYSUBS/MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200)
      COMMON/RYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/RYINT1/MINT(400),VINT(400)
      COMMON/RYINT2/ISET(200),KFPR(200,2),COEF(200,20),ICOL(40,4,2)
      COMMON/RYINT4/WIDP(21:40,0:40),WIDE(21:40,0:40),WIDS(21:40,3)
      COMMON/RYINT6/PROC(0:200)
      CHARACTER PROC*28
      SAVE /LUDAT1/,/LUDAT2/,/LUDAT3/,/LUDAT4/
      SAVE /RYSUBS/,/RYPARS/,/RYINT1/,/RYINT2/,/RYINT4/,/RYINT6/
      DIMENSION WDTP(0:40),WDTE(0:40,0:5)

C...Born level couplings in MSSM Higgs doublet sector.
      IF(MSTP(4).EQ.2) THEN
        TANBE=PARU(141)
        RATBE=((1.-TANBE**2)/(1.+TANBE**2))**2
        SQMZ=PMAS(23,1)**2
        SQMW=PMAS(24,1)**2
        SQMH=PMAS(25,1)**2
        SQMA=SQMH*(SQMZ-SQMH)/(SQMZ*RATBE-SQMH)
        SQMHP=0.5*(SQMA+SQMZ+SQRT((SQMA+SQMZ)**2-4.*SQMA*SQMZ*RATBE))
        SQMHC=SQMA+SQMW
        IF(MIN(SQMA,SQMHP,SQMHC).LE.0.) THEN
          WRITE(MSTU(11),5000)
          STOP
        ENDIF
        PMAS(35,1)=SQRT(SQMHP)
        PMAS(36,1)=SQRT(SQMA)
        PMAS(37,1)=SQRT(SQMHC)
        ALSU=0.5*ATAN(2.*TANBE*(SQMA+SQMZ)/((1.-TANBE**2)*
     &  (SQMA-SQMZ)))
        BESU=ATAN(TANBE)
        PARU(161)=-SIN(ALSU)/COS(BESU)
        PARU(162)=COS(ALSU)/SIN(BESU)
        PARU(163)=PARU(161)
        PARU(164)=SIN(BESU-ALSU)
        PARU(165)=PARU(164)
        PARU(171)=COS(ALSU)/COS(BESU)
        PARU(172)=SIN(ALSU)/SIN(BESU)
        PARU(173)=PARU(171)
        PARU(174)=COS(BESU-ALSU)
        PARU(175)=PARU(174)
        PARU(176)=COS(2.*ALSU)*COS(BESU+ALSU)-2.*SIN(2.*ALSU)*
     &  SIN(BESU+ALSU)
        PARU(177)=COS(2.*BESU)*COS(BESU+ALSU)
        PARU(181)=-TANBE
        PARU(182)=-1./TANBE
        PARU(183)=PARU(181)
        PARU(184)=0.
        PARU(185)=PARU(184)
        PARU(186)=COS(BESU-ALSU)
        PARU(187)=SIN(BESU-ALSU)
        PARU(195)=COS(BESU-ALSU)
      ENDIF

C...Reset full and effective widths of gauge bosons.
      XW=PARU(102)
      DO 100 I=21,40
      DO 100 J=0,40
      WIDP(I,J)=0.
  100 WIDE(I,J)=0.

C...Loop over possible resonances.
      DO 130 I=1,10
      IF(I.LE.3) KC=I+22
      IF(I.EQ.4) KC=37
      IF(I.EQ.5) KC=36
      IF(I.EQ.6) KC=35
      IF(I.EQ.7) KC=32
      IF(I.EQ.8) KC=34
      IF(I.GE.9) KC=I+30

C...Check that no fourth generation channels on by mistake.
      IF(MSTP(1).LE.3) THEN
        DO 110 J=1,MDCY(KC,3)
        IDC=J+MDCY(KC,2)-1
        KFA1=IABS(KFDP(IDC,1))
        KFA2=IABS(KFDP(IDC,2))
  110   IF(KFA1.EQ.7.OR.KFA1.EQ.8.OR.KFA1.EQ.17.OR.KFA1.EQ.18.OR.KFA2.
     &  EQ.7.OR.KFA2.EQ.8.OR.KFA2.EQ.17.OR.KFA2.EQ.18) MDME(IDC,1)=-1
      ENDIF

C...Find mass and evaluate width.
      PMR=PMAS(KC,1)
      IF(KC.EQ.25.OR.KC.EQ.35.OR.KC.EQ.36) MINT(62)=1
      CALL RYWIDT(KC,PMR**2,WDTP,WDTE)

C...Evaluate suppression factors due to non-simulated channels.
      IF(KCHG(KC,3).EQ.0) THEN
        WIDS(KC,1)=((WDTE(0,1)+WDTE(0,2))**2+
     &  2.*(WDTE(0,1)+WDTE(0,2))*(WDTE(0,4)+WDTE(0,5))+
     &  2.*WDTE(0,4)*WDTE(0,5))/WDTP(0)**2
        WIDS(KC,2)=(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))/WDTP(0)
        WIDS(KC,3)=0.
      ELSE
        WIDS(KC,1)=((WDTE(0,1)+WDTE(0,2))*(WDTE(0,1)+WDTE(0,3))+
     &  (WDTE(0,1)+WDTE(0,2)+WDTE(0,1)+WDTE(0,3))*(WDTE(0,4)+WDTE(0,5))+
     &  2.*WDTE(0,4)*WDTE(0,5))/WDTP(0)**2
        WIDS(KC,2)=(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))/WDTP(0)
        WIDS(KC,3)=(WDTE(0,1)+WDTE(0,3)+WDTE(0,4))/WDTP(0)
        IF(KC.EQ.24) THEN
          VINT(91)=((WDTE(0,1)+WDTE(0,2))**2+2.*(WDTE(0,1)+WDTE(0,2))*
     &    (WDTE(0,4)+WDTE(0,5))+2.*WDTE(0,4)*WDTE(0,5))/WDTP(0)**2
          VINT(92)=((WDTE(0,1)+WDTE(0,3))**2+2.*(WDTE(0,1)+WDTE(0,3))*
     &    (WDTE(0,4)+WDTE(0,5))+2.*WDTE(0,4)*WDTE(0,5))/WDTP(0)**2
        ENDIF
      ENDIF

C...Find factors to give widths in GeV.
      AEM=ULALEM(PMR**2)
      IF(KC.EQ.23) THEN
        FAC=AEM/(48.*XW*(1.-XW))*PMR
      ELSEIF(KC.EQ.24) THEN
        FAC=AEM/(24.*XW)*PMR
      ELSEIF(KC.EQ.25.OR.KC.EQ.35.OR.KC.EQ.36) THEN
        FAC=AEM/(8.*XW)*(PMR/PMAS(24,1))**2*PMR
      ELSEIF(KC.EQ.32) THEN
        FAC=AEM/(48.*XW*(1.-XW))*PMR
      ELSEIF(KC.EQ.34) THEN
        FAC=AEM/(24.*XW)*PMR
      ELSEIF(KC.EQ.37) THEN
        FAC=AEM/(8.*XW)*(PMR/PMAS(24,1))**2*PMR
      ELSEIF(KC.EQ.39) THEN
        FAC=AEM/4.*PMR
      ELSEIF(KC.EQ.40) THEN
        FAC=AEM/(12.*XW)*PMR
      ENDIF

C...Translate widths into GeV and save them.
      DO 120 J=0,40
      WIDP(KC,J)=FAC*WDTP(J)
  120 WIDE(KC,J)=FAC*WDTE(J,0)

C...Set resonance widths and branching ratios in JETSET;
C...also on/off switch for decays in RYTHIA/JETSET.
      PMAS(KC,2)=WIDP(KC,0)
      PMAS(KC,3)=MIN(0.9*PMAS(KC,1),10.*PMAS(KC,2))
      MDCY(KC,1)=MSTP(41)
      DO 130 J=1,MDCY(KC,3)
      IDC=J+MDCY(KC,2)-1
      BRAT(IDC)=WIDE(KC,J)/WIDE(KC,0)
  130 CONTINUE

C...Find heaviest new quark flavour allowed in processes 81-84.
      KFLQM=1
      DO 140 I=1,MIN(8,MDCY(21,3))
      IDC=I+MDCY(21,2)-1
      IF(MDME(IDC,1).LE.0) GOTO 140
      KFLQM=I
  140 CONTINUE
      MINT(55)=KFLQM
      KFPR(81,1)=KFLQM
      KFPR(81,2)=KFLQM
      KFPR(82,1)=KFLQM
      KFPR(82,2)=KFLQM
      KFPR(83,1)=KFLQM
      KFPR(84,1)=KFLQM
      KFPR(84,2)=KFLQM

C...Find heaviest new fermion flavour allowed in process 85.
      KFLFM=1
      DO 150 I=1,MIN(12,MDCY(22,3))
      IDC=I+MDCY(22,2)-1
      IF(MDME(IDC,1).LE.0) GOTO 150
      KFLFM=KFDP(IDC,1)
  150 CONTINUE
      MINT(56)=KFLFM
      KFPR(85,1)=KFLFM
      KFPR(85,2)=KFLFM

C...Flavours of leptoquark: redefine charge and name.
      KFLQQ=KFDP(MDCY(39,2),1)
      KFLQL=KFDP(MDCY(39,2),2)
      KCHG(39,1)=KCHG(IABS(KFLQQ),1)*ISIGN(1,KFLQQ)+
     &KCHG(IABS(KFLQL),1)*ISIGN(1,KFLQL)
      CHAF(39)(4:4)=CHAF(IABS(KFLQQ))(1:1)
      CHAF(39)(5:7)=CHAF(IABS(KFLQL))(1:3)

C...Special cases in treatment of gamma*/Z0: redefine process name.
      IF(MSTP(43).EQ.1) THEN
        PROC(1)='f + f~ -> gamma*'
      ELSEIF(MSTP(43).EQ.2) THEN
        PROC(1)='f + f~ -> Z0'
      ELSEIF(MSTP(43).EQ.3) THEN
        PROC(1)='f + f~ -> gamma*/Z0'
      ENDIF

C...Special cases in treatment of gamma*/Z0/Z'0: redefine process name.
      IF(MSTP(44).EQ.1) THEN
        PROC(141)='f + f~ -> gamma*'
      ELSEIF(MSTP(44).EQ.2) THEN
        PROC(141)='f + f~ -> Z0'
      ELSEIF(MSTP(44).EQ.3) THEN
        PROC(141)='f + f~ -> Z''0'
      ELSEIF(MSTP(44).EQ.4) THEN
        PROC(141)='f + f~ -> gamma*/Z0'
      ELSEIF(MSTP(44).EQ.5) THEN
        PROC(141)='f + f~ -> gamma*/Z''0'
      ELSEIF(MSTP(44).EQ.6) THEN
        PROC(141)='f + f~ -> Z0/Z''0'
      ELSEIF(MSTP(44).EQ.7) THEN
        PROC(141)='f + f~ -> gamma*/Z0/Z''0'
      ENDIF

C...Special cases in treatment of WW -> WW: redefine process name.
      IF(MSTP(45).EQ.1) THEN
        PROC(77)='W+ + W+ -> W+ + W+'
      ELSEIF(MSTP(45).EQ.2) THEN
        PROC(77)='W+ + W- -> W+ + W-'
      ELSEIF(MSTP(45).EQ.3) THEN
        PROC(77)='W+/- + W+/- -> W+/- + W+/-'
      ENDIF

C...Format for error information.
 5000 FORMAT(1X,'Error: unphysical input tan^2(beta) and m_H ',
     &'combination'/1X,'Execution stopped!')

      RETURN
      END

C*********************************************************************

      SUBROUTINE RYXTOT

C...Parametrizes total, double diffractive, single diffractive and
C...elastic cross-sections for different energies and beams.
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/RYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/RYINT1/MINT(400),VINT(400)
      COMMON/RYINT5/NGEN(0:200,3),XSEC(0:200,3)
      SAVE /LUDAT1/
      SAVE /RYPARS/,/RYINT1/,/RYINT5/
      DIMENSION BCS(5,8),BCB(2,5),BCC(3)

C...The following data lines are coefficients needed in the
C...Block, Cahn parametrization of total cross-section and nuclear
C...slope parameter; see below.
      DATA ((BCS(I,J),J=1,8),I=1,5)/
     1 41.74, 0.66, 0.0000, 337.,  0.0, 0.0, -39.3, 0.48,
     2 41.66, 0.60, 0.0000, 306.,  0.0, 0.0, -34.6, 0.51,
     3 41.36, 0.63, 0.0000, 299.,  7.3, 0.5, -40.4, 0.47,
     4 41.68, 0.63, 0.0083, 330.,  0.0, 0.0, -39.0, 0.48,
     5 41.13, 0.59, 0.0074, 278., 10.5, 0.5, -41.2, 0.46/
      DATA ((BCB(I,J),J=1,5),I=1,2)/
     1 10.79, -0.049, 0.040, 21.5, 1.23,
     2  9.92, -0.027, 0.013, 18.9, 1.07/
      DATA BCC/2.0164346,-0.5590311,0.0376279/

C...Total cross-section and nuclear slope parameter for pp and p-pbar
      NFIT=MIN(5,MAX(1,MSTP(31)))
      SIGP=BCS(NFIT,1)+BCS(NFIT,2)*(-0.25*PARU(1)**2*
     &(1.-0.25*BCS(NFIT,3)*PARU(1)**2)+(1.+0.5*BCS(NFIT,3)*PARU(1)**2)*
     &(LOG(VINT(2)/BCS(NFIT,4)))**2+BCS(NFIT,3)*
     &(LOG(VINT(2)/BCS(NFIT,4)))**4)/
     &((1.-0.25*BCS(NFIT,3)*PARU(1)**2)**2+2.*BCS(NFIT,3)*
     &(1.+0.25*BCS(NFIT,3)*PARU(1)**2)*(LOG(VINT(2)/BCS(NFIT,4)))**2+
     &BCS(NFIT,3)**2*(LOG(VINT(2)/BCS(NFIT,4)))**4)+BCS(NFIT,5)*
     &VINT(2)**(BCS(NFIT,6)-1.)*SIN(0.5*PARU(1)*BCS(NFIT,6))
      SIGM=-BCS(NFIT,7)*VINT(2)**(BCS(NFIT,8)-1.)*
     &COS(0.5*PARU(1)*BCS(NFIT,8))
      REFP=BCS(NFIT,2)*PARU(1)*LOG(VINT(2)/BCS(NFIT,4))/
     &((1.-0.25*BCS(NFIT,3)*PARU(1)**2)**2+2.*BCS(NFIT,3)*
     &(1.+0.25*BCS(NFIT,3)*PARU(1)**2)+(LOG(VINT(2)/BCS(NFIT,4)))**2+
     &BCS(NFIT,3)**2*(LOG(VINT(2)/BCS(NFIT,4)))**4)-BCS(NFIT,5)*
     &VINT(2)**(BCS(NFIT,6)-1.)*COS(0.5*PARU(1)*BCS(NFIT,6))
      REFM=-BCS(NFIT,7)*VINT(2)**(BCS(NFIT,8)-1.)*
     &SIN(0.5*PARU(1)*BCS(NFIT,8))
      SIGMA=SIGP-ISIGN(1,MINT(11)*MINT(12))*SIGM
      RHO=(REFP-ISIGN(1,MINT(11)*MINT(12))*REFM)/SIGMA

C...Nuclear slope parameter B, curvature C:
      NFIT=1
      IF(MSTP(31).GE.4) NFIT=2
      BP=BCB(NFIT,1)+BCB(NFIT,2)*LOG(VINT(2))+
     &BCB(NFIT,3)*(LOG(VINT(2)))**2
      BM=BCB(NFIT,4)+BCB(NFIT,5)*LOG(VINT(2))
      B=BP-ISIGN(1,MINT(11)*MINT(12))*SIGM/SIGP*(BM-BP)
      VINT(121)=B
      C=-0.5*BCC(2)/BCC(3)*(1.-SQRT(MAX(0.,1.+4.*BCC(3)/BCC(2)**2*
     &(1.E-03*VINT(1)-BCC(1)))))
      VINT(122)=C

C...Elastic scattering cross-section (fixed by sigma-tot, rho and B).
      SIGEL=SIGMA**2*(1.+RHO**2)/(16.*PARU(1)*PARU(5)*B)

C...Single diffractive scattering cross-section from Goulianos:
      SIGSD=2.*0.68*(1.+36./VINT(2))*LOG(0.6+0.1*VINT(2))

C...Double diffractive scattering cross-section (essentially fixed by
C...sigma-sd and sigma-el).
      SIGDD=SIGSD**2/(3.*SIGEL)

C...Total non-elastic, non-diffractive cross-section.
      SIGND=SIGMA-SIGDD-SIGSD-SIGEL

C...Rescale for pions.
      IF(IABS(MINT(11)).EQ.211.AND.IABS(MINT(12)).EQ.211) THEN
        SIGMA=4./9.*SIGMA
        SIGDD=4./9.*SIGDD
        SIGSD=4./9.*SIGSD
        SIGEL=4./9.*SIGEL
        SIGND=4./9.*SIGND
      ELSEIF(IABS(MINT(11)).EQ.211.OR.IABS(MINT(12)).EQ.211) THEN
        SIGMA=2./3.*SIGMA
        SIGDD=2./3.*SIGDD
        SIGSD=2./3.*SIGSD
        SIGEL=2./3.*SIGEL
        SIGND=2./3.*SIGND
      ENDIF

C...Save cross-sections in common block RYPARA.
      VINT(101)=SIGMA
      VINT(102)=SIGEL
      VINT(103)=SIGSD
      VINT(104)=SIGDD
      VINT(106)=SIGND
      XSEC(95,1)=SIGND

      RETURN
      END

C*********************************************************************

      SUBROUTINE RYMAXI

C...Finds optimal set of coefficients for kinematical variable selection
C...and the maximum of the part of the differential cross-section used
C...in the event weighting.
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/RYSUBS/MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200)
      COMMON/RYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/RYINT1/MINT(400),VINT(400)
      COMMON/RYINT2/ISET(200),KFPR(200,2),COEF(200,20),ICOL(40,4,2)
      COMMON/RYINT3/XSFX(2,-40:40),ISIG(1000,3),SIGH(1000)
      COMMON/RYINT4/WIDP(21:40,0:40),WIDE(21:40,0:40),WIDS(21:40,3)
      COMMON/RYINT5/NGEN(0:200,3),XSEC(0:200,3)
      COMMON/RYINT6/PROC(0:200)
      CHARACTER PROC*28
      SAVE /LUDAT1/,/LUDAT2/
      SAVE /RYSUBS/,/RYPARS/,/RYINT1/,/RYINT2/,/RYINT3/,/RYINT4/,
     &/RYINT5/,/RYINT6/
      CHARACTER CVAR(4)*4
      DIMENSION NPTS(4),MVARPT(500,4),VINTPT(500,30),SIGSPT(500),
     &NAREL(7),WTREL(7),WTMAT(7,7),WTRELN(7),COEFU(7),COEFO(7),
     &IACCMX(4),SIGSMX(4),SIGSSM(3)
      DATA CVAR/'tau ','tau''','y*  ','cth '/

C...Select subprocess to study: skip cases not applicable.
      NPOSI=0
      VINT(143)=1.
      VINT(144)=1.
      XSEC(0,1)=0.
      DO 410 ISUB=1,200
      IF(ISUB.GE.91.AND.ISUB.LE.95) THEN
        XSEC(ISUB,1)=VINT(ISUB+11)
        IF(MSUB(ISUB).NE.1) GOTO 410
        NPOSI=NPOSI+1
        GOTO 400
      ELSEIF(ISUB.EQ.96) THEN
        IF(MINT(44).NE.4) GOTO 410
        IF(MSUB(95).NE.1.AND.MSTP(81).LE.0.AND.MSTP(131).LE.0) GOTO 410
      ELSEIF(ISUB.EQ.11.OR.ISUB.EQ.12.OR.ISUB.EQ.13.OR.ISUB.EQ.28.OR.
     &ISUB.EQ.53.OR.ISUB.EQ.68) THEN
        IF(MSUB(ISUB).NE.1.OR.MSUB(95).EQ.1) GOTO 410
      ELSE
        IF(MSUB(ISUB).NE.1) GOTO 410
      ENDIF
      MINT(1)=ISUB
      ISTSB=ISET(ISUB)
      IF(ISUB.EQ.96) ISTSB=2
      IF(MSTP(122).GE.2) WRITE(MSTU(11),5000) ISUB

C...Find resonances (explicit or implicit in cross-section).
      MINT(72)=0
      KFR1=0
      IF(ISTSB.EQ.1.OR.ISTSB.EQ.3.OR.ISTSB.EQ.5) THEN
        KFR1=KFPR(ISUB,1)
      ELSEIF(ISUB.EQ.24.OR.ISUB.EQ.25.OR.ISUB.EQ.171.OR.ISUB.EQ.176)
     &THEN
        KFR1=23
      ELSEIF(ISUB.EQ.23.OR.ISUB.EQ.26.OR.ISUB.EQ.172.OR.ISUB.EQ.177)
     &THEN
        KFR1=24
      ELSEIF(ISUB.GE.71.AND.ISUB.LE.77) THEN
        KFR1=25
        IF(MSTP(46).EQ.5) THEN
          KFR1=30
          PMAS(30,1)=PARP(45)
          PMAS(30,2)=PARP(45)**3/(96.*PARU(1)*246.**2)
        ENDIF
      ENDIF
      CKMX=CKIN(2)
      IF(CKMX.LE.0.) CKMX=VINT(1)
      IF(KFR1.NE.0) THEN
        IF(CKIN(1).GT.PMAS(KFR1,1)+20.*PMAS(KFR1,2).OR.
     &  CKMX.LT.PMAS(KFR1,1)-20.*PMAS(KFR1,2)) KFR1=0
      ENDIF
      IF(KFR1.NE.0) THEN
        TAUR1=PMAS(KFR1,1)**2/VINT(2)
        GAMR1=PMAS(KFR1,1)*PMAS(KFR1,2)/VINT(2)
        MINT(72)=1
        MINT(73)=KFR1
        VINT(73)=TAUR1
        VINT(74)=GAMR1
      ENDIF
      IF(ISUB.EQ.141) THEN
        KFR2=23
        TAUR2=PMAS(KFR2,1)**2/VINT(2)
        GAMR2=PMAS(KFR2,1)*PMAS(KFR2,2)/VINT(2)
        IF(CKIN(1).GT.PMAS(KFR2,1)+20.*PMAS(KFR2,2).OR.
     &  CKMX.LT.PMAS(KFR2,1)-20.*PMAS(KFR2,2)) KFR2=0
        IF(KFR2.NE.0.AND.KFR1.NE.0) THEN
          MINT(72)=2
          MINT(74)=KFR2
          VINT(75)=TAUR2
          VINT(76)=GAMR2
        ELSEIF(KFR2.NE.0) THEN
          KFR1=KFR2
          TAUR1=TAUR2
          GAMR1=GAMR2
          MINT(72)=1
          MINT(73)=KFR1
          VINT(73)=TAUR1
          VINT(74)=GAMR1
        ENDIF
      ENDIF

C...Find product masses and minimum pT of process.
      SQM3=0.
      SQM4=0.
      MINT(71)=0
      VINT(71)=CKIN(3)
      VINT(80)=1.
      IF(ISTSB.EQ.2.OR.ISTSB.EQ.4) THEN
        NBW=0
        DO 100 I=1,2
        IF(KFPR(ISUB,I).EQ.0) THEN
        ELSEIF(MSTP(42).LE.0.OR.PMAS(KFPR(ISUB,I),2).LT.PARP(41)) THEN
          IF(I.EQ.1) SQM3=PMAS(KFPR(ISUB,I),1)**2
          IF(I.EQ.2) SQM4=PMAS(KFPR(ISUB,I),1)**2
        ELSE
          NBW=NBW+1
        ENDIF
  100   CONTINUE
        IF(NBW.GE.1) THEN
          CALL RYOFSH(3,0,KFPR(ISUB,1),KFPR(ISUB,2),0.,PQM3,PQM4)
          IF(MINT(51).EQ.1) THEN
            WRITE(MSTU(11),5100) ISUB
            MSUB(ISUB)=0
            GOTO 410
          ENDIF
          SQM3=PQM3**2
          SQM4=PQM4**2
        ENDIF
        IF(MIN(SQM3,SQM4).LT.CKIN(6)**2) MINT(71)=1
        IF(MINT(71).EQ.1) VINT(71)=MAX(CKIN(3),CKIN(5))
        IF(ISUB.EQ.96.AND.MSTP(82).LE.1) VINT(71)=PARP(81)
        IF(ISUB.EQ.96.AND.MSTP(82).GE.2) VINT(71)=0.08*PARP(82)
      ELSEIF(ISTSB.EQ.6) THEN
        CALL RYOFSH(5,0,KFPR(ISUB,1),KFPR(ISUB,2),0.,PQM3,PQM4)
        IF(MINT(51).EQ.1) THEN
          WRITE(MSTU(11),5100) ISUB
          MSUB(ISUB)=0
          GOTO 410
        ENDIF
        SQM3=PQM3**2
        SQM4=PQM4**2
      ENDIF
      VINT(63)=SQM3
      VINT(64)=SQM4

C...Prepare for additional variable choices in 2 -> 3.
      IF(ISTSB.EQ.5) THEN
        VINT(201)=0.
        VINT(206)=0.
        VINT(204)=PMAS(23,1)
        IF(ISUB.EQ.124) VINT(204)=PMAS(24,1)
        VINT(209)=VINT(204)
      ENDIF

C...Number of points for each variable: tau, tau', y*, cos(theta-hat).
      NPTS(1)=2+2*MINT(72)
      IF(MINT(47).EQ.1) THEN
        IF(ISTSB.EQ.1.OR.ISTSB.EQ.2.OR.ISTSB.EQ.6) NPTS(1)=1
      ELSEIF(MINT(47).EQ.5) THEN
        IF(ISTSB.LE.2.OR.ISTSB.GE.6) NPTS(1)=NPTS(1)+1
      ENDIF
      NPTS(2)=1
      IF(ISTSB.GE.3.AND.ISTSB.LE.5) THEN
        IF(MINT(47).GE.2) NPTS(2)=2
        IF(MINT(47).EQ.5) NPTS(2)=3
      ENDIF
      NPTS(3)=1
      IF(MINT(47).GE.4) NPTS(3)=3
      IF(MINT(45).EQ.3) NPTS(3)=NPTS(3)+1
      IF(MINT(46).EQ.3) NPTS(3)=NPTS(3)+1
      NPTS(4)=1
      IF(ISTSB.EQ.2.OR.ISTSB.EQ.4.OR.ISTSB.EQ.6) NPTS(4)=5
      NTRY=NPTS(1)*NPTS(2)*NPTS(3)*NPTS(4)

C...Reset coefficients of cross-section weighting.
      DO 110 J=1,20
  110 COEF(ISUB,J)=0.
      COEF(ISUB,1)=1.
      COEF(ISUB,8)=0.5
      COEF(ISUB,9)=0.5
      COEF(ISUB,13)=1.
      COEF(ISUB,18)=1.
      MCTH=0
      MTAUP=0
      METAUP=0
      VINT(23)=0.
      VINT(26)=0.
      SIGSAM=0.

C...Find limits and select tau, y*, cos(theta-hat) and tau' values,
C...in grid of phase space points.
      CALL RYKLIM(1)
      METAU=MINT(51)
      NACC=0
      DO 140 ITRY=1,NTRY
      IF(METAU.EQ.1) GOTO 140
      IF(MOD(ITRY-1,NPTS(2)*NPTS(3)*NPTS(4)).EQ.0) THEN
        MTAU=1+(ITRY-1)/(NPTS(2)*NPTS(3)*NPTS(4))
        IF(MTAU.GT.2+2*MINT(72)) MTAU=7
        CALL RYKMAP(1,MTAU,0.5)
        IF(ISTSB.GE.3.AND.ISTSB.LE.5) CALL RYKLIM(4)
        METAUP=MINT(51)
      ENDIF
      IF(METAUP.EQ.1) GOTO 140
      IF(ISTSB.GE.3.AND.ISTSB.LE.5.AND.MOD(ITRY-1,NPTS(3)*NPTS(4)).
     &EQ.0) THEN
        MTAUP=1+MOD((ITRY-1)/(NPTS(3)*NPTS(4)),NPTS(2))
        CALL RYKMAP(4,MTAUP,0.5)
      ENDIF
      IF(MOD(ITRY-1,NPTS(3)*NPTS(4)).EQ.0) THEN
        CALL RYKLIM(2)
        MEYST=MINT(51)
      ENDIF
      IF(MEYST.EQ.1) GOTO 140
      IF(MOD(ITRY-1,NPTS(4)).EQ.0) THEN
        MYST=1+MOD((ITRY-1)/NPTS(4),NPTS(3))
        IF(MYST.EQ.4.AND.MINT(45).NE.3) MYST=5
        CALL RYKMAP(2,MYST,0.5)
        CALL RYKLIM(3)
        MECTH=MINT(51)
      ENDIF
      IF(MECTH.EQ.1) GOTO 140
      IF(ISTSB.EQ.2.OR.ISTSB.EQ.4.OR.ISTSB.EQ.6) THEN
        MCTH=1+MOD(ITRY-1,NPTS(4))
        CALL RYKMAP(3,MCTH,0.5)
      ENDIF
      IF(ISUB.EQ.96) VINT(25)=VINT(21)*(1.-VINT(23)**2)

C...Store position and limits.
      MINT(51)=0
      CALL RYKLIM(0)
      IF(MINT(51).EQ.1) GOTO 140
      NACC=NACC+1
      MVARPT(NACC,1)=MTAU
      MVARPT(NACC,2)=MTAUP
      MVARPT(NACC,3)=MYST
      MVARPT(NACC,4)=MCTH
      DO 120 J=1,30
  120 VINTPT(NACC,J)=VINT(10+J)

C...Normal case: calculate cross-section.
      IF(ISTSB.NE.5) THEN
        CALL RYSIGH(NCHN,SIGS)

C..2 -> 3: find highest value out of a number of tries.
      ELSE
        SIGS=0.
        DO 130 IKIN3=1,MSTP(129)
        CALL RYKMAP(5,0,0.)
        IF(MINT(51).EQ.1) GOTO 130
        CALL RYSIGH(NCHN,SIGTMP)
        IF(SIGTMP.GT.SIGS) SIGS=SIGTMP
  130   CONTINUE
      ENDIF

C...Store cross-section.
      SIGSPT(NACC)=SIGS
      IF(SIGS.GT.SIGSAM) SIGSAM=SIGS
      IF(MSTP(122).GE.2) WRITE(MSTU(11),5200) MTAU,MYST,MCTH,MTAUP,
     &VINT(21),VINT(22),VINT(23),VINT(26),SIGS
  140 CONTINUE
      IF(NACC.EQ.0) THEN
        WRITE(MSTU(11),5100) ISUB
        MSUB(ISUB)=0
        GOTO 410
      ELSEIF(SIGSAM.EQ.0.) THEN
        WRITE(MSTU(11),5300) ISUB
        MSUB(ISUB)=0
        GOTO 410
      ENDIF
      IF(ISUB.NE.96) NPOSI=NPOSI+1

C...Calculate integrals in tau over maximal phase space limits.
      TAUMIN=VINT(11)
      TAUMAX=VINT(31)
      ATAU1=LOG(TAUMAX/TAUMIN)
      IF(NPTS(1).GE.2) THEN
        ATAU2=(TAUMAX-TAUMIN)/(TAUMAX*TAUMIN)
      ENDIF
      IF(NPTS(1).GE.4) THEN
        ATAU3=LOG(TAUMAX/TAUMIN*(TAUMIN+TAUR1)/(TAUMAX+TAUR1))/TAUR1
        ATAU4=(ATAN((TAUMAX-TAUR1)/GAMR1)-ATAN((TAUMIN-TAUR1)/GAMR1))/
     &  GAMR1
      ENDIF
      IF(NPTS(1).GE.6) THEN
        ATAU5=LOG(TAUMAX/TAUMIN*(TAUMIN+TAUR2)/(TAUMAX+TAUR2))/TAUR2
        ATAU6=(ATAN((TAUMAX-TAUR2)/GAMR2)-ATAN((TAUMIN-TAUR2)/GAMR2))/
     &  GAMR2
      ENDIF
      IF(NPTS(1).GT.2+2*MINT(72)) THEN
        ATAU7=LOG(MAX(2E-6,1.-TAUMIN)/MAX(2E-6,1.-TAUMAX))
      ENDIF

C...Reset. Sum up cross-sections in points calculated.
      DO 270 IVAR=1,4
      IF(NPTS(IVAR).EQ.1) GOTO 270
      IF(ISUB.EQ.96.AND.IVAR.EQ.4) GOTO 270
      NBIN=NPTS(IVAR)
      DO 150 J1=1,NBIN
      NAREL(J1)=0
      WTREL(J1)=0.
      COEFU(J1)=0.
      DO 150 J2=1,NBIN
  150 WTMAT(J1,J2)=0.
      DO 160 IACC=1,NACC
      IBIN=MVARPT(IACC,IVAR)
      IF(IVAR.EQ.1.AND.IBIN.EQ.7) IBIN=3+2*MINT(72)
      IF(IVAR.EQ.3.AND.IBIN.EQ.5.AND.MINT(45).NE.3) IBIN=4
      NAREL(IBIN)=NAREL(IBIN)+1
      WTREL(IBIN)=WTREL(IBIN)+SIGSPT(IACC)

C...Sum up tau cross-section pieces in points used.
      IF(IVAR.EQ.1) THEN
        TAU=VINTPT(IACC,11)
        WTMAT(IBIN,1)=WTMAT(IBIN,1)+1.
        WTMAT(IBIN,2)=WTMAT(IBIN,2)+(ATAU1/ATAU2)/TAU
        IF(NBIN.GE.4) THEN
          WTMAT(IBIN,3)=WTMAT(IBIN,3)+(ATAU1/ATAU3)/(TAU+TAUR1)
          WTMAT(IBIN,4)=WTMAT(IBIN,4)+(ATAU1/ATAU4)*TAU/
     &    ((TAU-TAUR1)**2+GAMR1**2)
        ENDIF
        IF(NBIN.GE.6) THEN
          WTMAT(IBIN,5)=WTMAT(IBIN,5)+(ATAU1/ATAU5)/(TAU+TAUR2)
          WTMAT(IBIN,6)=WTMAT(IBIN,6)+(ATAU1/ATAU6)*TAU/
     &    ((TAU-TAUR2)**2+GAMR2**2)
        ENDIF
        IF(NBIN.GT.2+2*MINT(72)) THEN
          WTMAT(IBIN,NBIN)=WTMAT(IBIN,NBIN)+(ATAU1/ATAU7)*
     &    TAU/MAX(2E-6,1.-TAU)
        ENDIF

C...Sum up tau' cross-section pieces in points used.
      ELSEIF(IVAR.EQ.2) THEN
        TAU=VINTPT(IACC,11)
        TAUP=VINTPT(IACC,16)
        TAUPMN=VINTPT(IACC,6)
        TAUPMX=VINTPT(IACC,26)
        ATAUP1=LOG(TAUPMX/TAUPMN)
        ATAUP2=((1.-TAU/TAUPMX)**4-(1.-TAU/TAUPMN)**4)/(4.*TAU)
        WTMAT(IBIN,1)=WTMAT(IBIN,1)+1.
        WTMAT(IBIN,2)=WTMAT(IBIN,2)+(ATAUP1/ATAUP2)*(1.-TAU/TAUP)**3/
     &  TAUP
        IF(NBIN.GE.3) THEN
          ATAUP3=LOG(MAX(2E-6,1.-TAUPMN)/MAX(2E-6,1.-TAUPMX))
          WTMAT(IBIN,3)=WTMAT(IBIN,3)+(ATAUP1/ATAUP3)*
     &    TAUP/MAX(2E-6,1.-TAUP)
        ENDIF

C...Sum up y* cross-section pieces in points used.
      ELSEIF(IVAR.EQ.3) THEN
        YST=VINTPT(IACC,12)
        YSTMIN=VINTPT(IACC,2)
        YSTMAX=VINTPT(IACC,22)
        AYST0=YSTMAX-YSTMIN
        AYST1=0.5*(YSTMAX-YSTMIN)**2
        AYST2=AYST1
        AYST3=2.*(ATAN(EXP(YSTMAX))-ATAN(EXP(YSTMIN)))
        WTMAT(IBIN,1)=WTMAT(IBIN,1)+(AYST0/AYST1)*(YST-YSTMIN)
        WTMAT(IBIN,2)=WTMAT(IBIN,2)+(AYST0/AYST2)*(YSTMAX-YST)
        WTMAT(IBIN,3)=WTMAT(IBIN,3)+(AYST0/AYST3)/COSH(YST)
        IF(MINT(45).EQ.3) THEN
          TAUE=VINTPT(IACC,11)
          IF(ISTSB.GE.3.AND.ISTSB.LE.5) TAUE=VINTPT(IACC,16)
          YST0=-0.5*LOG(TAUE)
          AYST4=LOG(MAX(1E-6,EXP(YST0-YSTMIN)-1.)/
     &    MAX(1E-6,EXP(YST0-YSTMAX)-1.))
          WTMAT(IBIN,4)=WTMAT(IBIN,4)+(AYST0/AYST4)/
     &    MAX(1E-6,1.-EXP(YST-YST0))
        ENDIF
        IF(MINT(46).EQ.3) THEN
          TAUE=VINTPT(IACC,11)
          IF(ISTSB.GE.3.AND.ISTSB.LE.5) TAUE=VINTPT(IACC,16)
          YST0=-0.5*LOG(TAUE)
          AYST5=LOG(MAX(1E-6,EXP(YST0+YSTMAX)-1.)/
     &    MAX(1E-6,EXP(YST0+YSTMIN)-1.))
          WTMAT(IBIN,NBIN)=WTMAT(IBIN,NBIN)+(AYST0/AYST5)/
     &    MAX(1E-6,1.-EXP(-YST-YST0))
        ENDIF

C...Sum up cos(theta-hat) cross-section pieces in points used.
      ELSE
        RM34=2.*SQM3*SQM4/(VINTPT(IACC,11)*VINT(2))**2
        RSQM=1.+RM34
        CTHMAX=SQRT(1.-4.*VINT(71)**2/(TAUMAX*VINT(2)))
        CTHMIN=-CTHMAX
        IF(CTHMAX.GT.0.9999) RM34=MAX(RM34,2.*VINT(71)**2/
     &  (TAUMAX*VINT(2)))
        ACTH1=CTHMAX-CTHMIN
        ACTH2=LOG(MAX(RM34,RSQM-CTHMIN)/MAX(RM34,RSQM-CTHMAX))
        ACTH3=LOG(MAX(RM34,RSQM+CTHMAX)/MAX(RM34,RSQM+CTHMIN))
        ACTH4=1./MAX(RM34,RSQM-CTHMAX)-1./MAX(RM34,RSQM-CTHMIN)
        ACTH5=1./MAX(RM34,RSQM+CTHMIN)-1./MAX(RM34,RSQM+CTHMAX)
        CTH=VINTPT(IACC,13)
        WTMAT(IBIN,1)=WTMAT(IBIN,1)+1.
        WTMAT(IBIN,2)=WTMAT(IBIN,2)+(ACTH1/ACTH2)/MAX(RM34,RSQM-CTH)
        WTMAT(IBIN,3)=WTMAT(IBIN,3)+(ACTH1/ACTH3)/MAX(RM34,RSQM+CTH)
        WTMAT(IBIN,4)=WTMAT(IBIN,4)+(ACTH1/ACTH4)/MAX(RM34,RSQM-CTH)**2
        WTMAT(IBIN,5)=WTMAT(IBIN,5)+(ACTH1/ACTH5)/MAX(RM34,RSQM+CTH)**2
      ENDIF
  160 CONTINUE

C...Check that equation system solvable; else trivial way out.
      IF(MSTP(122).GE.2) WRITE(MSTU(11),5400) CVAR(IVAR)
      MSOLV=1
      WTRELS=0.
      DO 170 IBIN=1,NBIN
      IF(MSTP(122).GE.2) WRITE(MSTU(11),5500) (WTMAT(IBIN,IRED),
     &IRED=1,NBIN),WTREL(IBIN)
      IF(NAREL(IBIN).EQ.0) MSOLV=0
  170 WTRELS=WTRELS+WTREL(IBIN)
      IF(MSOLV.EQ.0) THEN
        DO 180 IBIN=1,NBIN
        COEFU(IBIN)=1.
        WTRELN(IBIN)=0.1
  180   IF(WTRELS.GT.0.) WTRELN(IBIN)=MAX(0.1,WTREL(IBIN)/WTRELS)

C...Solve to find relative importance of cross-section pieces.
      ELSE
        DO 190 IBIN=1,NBIN
  190   WTRELN(IBIN)=MAX(0.1,WTREL(IBIN)/WTRELS)
        DO 200 IRED=1,NBIN-1
        DO 200 IBIN=IRED+1,NBIN
        RQT=WTMAT(IBIN,IRED)/WTMAT(IRED,IRED)
        WTREL(IBIN)=WTREL(IBIN)-RQT*WTREL(IRED)
        DO 200 ICOE=IRED,NBIN
  200   WTMAT(IBIN,ICOE)=WTMAT(IBIN,ICOE)-RQT*WTMAT(IRED,ICOE)
        DO 220 IRED=NBIN,1,-1
        DO 210 ICOE=IRED+1,NBIN
  210   WTREL(IRED)=WTREL(IRED)-WTMAT(IRED,ICOE)*COEFU(ICOE)
  220   COEFU(IRED)=WTREL(IRED)/WTMAT(IRED,IRED)
      ENDIF

C...Normalize coefficients, with piece shared democratically.
      COEFSU=0.
      WTRELS=0.
      DO 230 IBIN=1,NBIN
      COEFU(IBIN)=MAX(0.,COEFU(IBIN))
      COEFSU=COEFSU+COEFU(IBIN)
  230 WTRELS=WTRELS+WTRELN(IBIN)
      IF(COEFSU.GT.0.) THEN
        DO 240 IBIN=1,NBIN
  240   COEFO(IBIN)=PARP(122)/NBIN+(1.-PARP(122))*0.5*
     &  (COEFU(IBIN)/COEFSU+WTRELN(IBIN)/WTRELS)
      ELSE
        DO 250 IBIN=1,NBIN
  250   COEFO(IBIN)=1./NBIN
      ENDIF
      IF(IVAR.EQ.1) IOFF=0
      IF(IVAR.EQ.2) IOFF=17
      IF(IVAR.EQ.3) IOFF=7
      IF(IVAR.EQ.4) IOFF=12
      DO 260 IBIN=1,NBIN
      ICOF=IOFF+IBIN
      IF(IVAR.EQ.1.AND.IBIN.GT.2+2*MINT(72)) ICOF=7
      IF(IVAR.EQ.3.AND.IBIN.EQ.4.AND.MINT(45).NE.3) ICOF=ICOF+1
  260 COEF(ISUB,ICOF)=COEFO(IBIN)
      IF(MSTP(122).GE.2) WRITE(MSTU(11),5600) CVAR(IVAR),
     &(COEFO(IBIN),IBIN=1,NBIN)
  270 CONTINUE

C...Find two most promising maxima among points previously determined.
      DO 280 J=1,4
      IACCMX(J)=0
  280 SIGSMX(J)=0.
      NMAX=0
      DO 340 IACC=1,NACC
      DO 290 J=1,30
  290 VINT(10+J)=VINTPT(IACC,J)
      IF(ISTSB.NE.5) THEN
        CALL RYSIGH(NCHN,SIGS)
      ELSE
        SIGS=0.
        DO 300 IKIN3=1,MSTP(129)
        CALL RYKMAP(5,0,0.)
        IF(MINT(51).EQ.1) GOTO 300
        CALL RYSIGH(NCHN,SIGTMP)
        IF(SIGTMP.GT.SIGS) SIGS=SIGTMP
  300   CONTINUE
      ENDIF
      IEQ=0
      DO 310 IMV=1,NMAX
  310 IF(ABS(SIGS-SIGSMX(IMV)).LT.1E-4*(SIGS+SIGSMX(IMV))) IEQ=IMV
      IF(IEQ.EQ.0) THEN
        DO 320 IMV=NMAX,1,-1
        IIN=IMV+1
        IF(SIGS.LE.SIGSMX(IMV)) GOTO 330
        IACCMX(IMV+1)=IACCMX(IMV)
  320   SIGSMX(IMV+1)=SIGSMX(IMV)
        IIN=1
  330   IACCMX(IIN)=IACC
        SIGSMX(IIN)=SIGS
        IF(NMAX.LE.1) NMAX=NMAX+1
      ENDIF
  340 CONTINUE

C...Read out starting position for search.
      IF(MSTP(122).GE.2) WRITE(MSTU(11),5700)
      SIGSAM=SIGSMX(1)
      DO 390 IMAX=1,NMAX
      IACC=IACCMX(IMAX)
      MTAU=MVARPT(IACC,1)
      MTAUP=MVARPT(IACC,2)
      MYST=MVARPT(IACC,3)
      MCTH=MVARPT(IACC,4)
      VTAU=0.5
      VYST=0.5
      VCTH=0.5
      VTAUP=0.5

C...Starting point and step size in parameter space.
      DO 380 IRPT=1,2
      DO 370 IVAR=1,4
      IF(NPTS(IVAR).EQ.1) GOTO 370
      IF(IVAR.EQ.1) VVAR=VTAU
      IF(IVAR.EQ.2) VVAR=VTAUP
      IF(IVAR.EQ.3) VVAR=VYST
      IF(IVAR.EQ.4) VVAR=VCTH
      IF(IVAR.EQ.1) MVAR=MTAU
      IF(IVAR.EQ.2) MVAR=MTAUP
      IF(IVAR.EQ.3) MVAR=MYST
      IF(IVAR.EQ.4) MVAR=MCTH
      IF(IRPT.EQ.1) VDEL=0.1
      IF(IRPT.EQ.2) VDEL=MAX(0.01,MIN(0.05,VVAR-0.02,0.98-VVAR))
      IF(IRPT.EQ.1) VMAR=0.02
      IF(IRPT.EQ.2) VMAR=0.002
      IMOV0=1
      IF(IRPT.EQ.1.AND.IVAR.EQ.1) IMOV0=0
      DO 360 IMOV=IMOV0,8

C...Define new point in parameter space.
      IF(IMOV.EQ.0) THEN
        INEW=2
        VNEW=VVAR
      ELSEIF(IMOV.EQ.1) THEN
        INEW=3
        VNEW=VVAR+VDEL
      ELSEIF(IMOV.EQ.2) THEN
        INEW=1
        VNEW=VVAR-VDEL
      ELSEIF(SIGSSM(3).GE.MAX(SIGSSM(1),SIGSSM(2)).AND.
     &VVAR+2.*VDEL.LT.1.-VMAR) THEN
        VVAR=VVAR+VDEL
        SIGSSM(1)=SIGSSM(2)
        SIGSSM(2)=SIGSSM(3)
        INEW=3
        VNEW=VVAR+VDEL
      ELSEIF(SIGSSM(1).GE.MAX(SIGSSM(2),SIGSSM(3)).AND.
     &VVAR-2.*VDEL.GT.VMAR) THEN
        VVAR=VVAR-VDEL
        SIGSSM(3)=SIGSSM(2)
        SIGSSM(2)=SIGSSM(1)
        INEW=1
        VNEW=VVAR-VDEL
      ELSEIF(SIGSSM(3).GE.SIGSSM(1)) THEN
        VDEL=0.5*VDEL
        VVAR=VVAR+VDEL
        SIGSSM(1)=SIGSSM(2)
        INEW=2
        VNEW=VVAR
      ELSE
        VDEL=0.5*VDEL
        VVAR=VVAR-VDEL
        SIGSSM(3)=SIGSSM(2)
        INEW=2
        VNEW=VVAR
      ENDIF

C...Convert to relevant variables and find derived new limits.
      IF(IVAR.EQ.1) THEN
        VTAU=VNEW
        CALL RYKMAP(1,MTAU,VTAU)
        IF(ISTSB.GE.3.AND.ISTSB.LE.5) CALL RYKLIM(4)
      ENDIF
      IF(IVAR.LE.2.AND.ISTSB.GE.3.AND.ISTSB.LE.5) THEN
        IF(IVAR.EQ.2) VTAUP=VNEW
        CALL RYKMAP(4,MTAUP,VTAUP)
      ENDIF
      IF(IVAR.LE.2) CALL RYKLIM(2)
      IF(IVAR.LE.3) THEN
        IF(IVAR.EQ.3) VYST=VNEW
        CALL RYKMAP(2,MYST,VYST)
        CALL RYKLIM(3)
      ENDIF
      IF(ISTSB.EQ.2.OR.ISTSB.EQ.4.OR.ISTSB.EQ.6) THEN
        IF(IVAR.EQ.4) VCTH=VNEW
        CALL RYKMAP(3,MCTH,VCTH)
      ENDIF
      IF(ISUB.EQ.96) VINT(25)=VINT(21)*(1.-VINT(23)**2)

C...Evaluate cross-section. Save new maximum. Final maximum.
      IF(ISTSB.NE.5) THEN
        CALL RYSIGH(NCHN,SIGS)
      ELSE
        SIGS=0.
        DO 350 IKIN3=1,MSTP(129)
        CALL RYKMAP(5,0,0.)
        IF(MINT(51).EQ.1) GOTO 350
        CALL RYSIGH(NCHN,SIGTMP)
        IF(SIGTMP.GT.SIGS) SIGS=SIGTMP
  350   CONTINUE
      ENDIF
      SIGSSM(INEW)=SIGS
      IF(SIGS.GT.SIGSAM) SIGSAM=SIGS
      IF(MSTP(122).GE.2) WRITE(MSTU(11),5800) IMAX,IVAR,MVAR,IMOV,
     &VNEW,VINT(21),VINT(22),VINT(23),VINT(26),SIGS
  360 CONTINUE
  370 CONTINUE
  380 CONTINUE
  390 CONTINUE
      IF(MSTP(121).EQ.1) SIGSAM=PARP(121)*SIGSAM
      XSEC(ISUB,1)=1.05*SIGSAM
  400 IF(ISUB.NE.96) XSEC(0,1)=XSEC(0,1)+XSEC(ISUB,1)
  410 CONTINUE

C...Print summary table.
      IF(NPOSI.EQ.0) THEN
        WRITE(MSTU(11),5900)
        STOP
      ENDIF
      IF(MSTP(122).GE.1) THEN
        WRITE(MSTU(11),6000)
        WRITE(MSTU(11),6100)
        DO 420 ISUB=1,200
        IF(MSUB(ISUB).NE.1.AND.ISUB.NE.96) GOTO 420
        IF(ISUB.EQ.96.AND.MINT(44).NE.4) GOTO 420
        IF(ISUB.EQ.96.AND.MSUB(95).NE.1.AND.MSTP(81).LE.0) GOTO 420
        IF(MSUB(95).EQ.1.AND.(ISUB.EQ.11.OR.ISUB.EQ.12.OR.ISUB.EQ.13.OR.
     &  ISUB.EQ.28.OR.ISUB.EQ.53.OR.ISUB.EQ.68)) GOTO 420
        WRITE(MSTU(11),6200) ISUB,PROC(ISUB),XSEC(ISUB,1)
  420   CONTINUE
        WRITE(MSTU(11),6300)
      ENDIF

C...Format statements for maximization results.
 5000 FORMAT(/1X,'Coefficient optimization and maximum search for ',
     &'subprocess no',I4/1X,'Coefficient modes     tau',10X,'y*',9X,
     &'cth',9X,'tau''',7X,'sigma')
 5100 FORMAT(1X,'Warning: requested subprocess ',I3,' has no allowed ',
     &'phase space.'/1X,'Process switched off!')
 5200 FORMAT(1X,4I4,F12.8,F12.6,F12.7,F12.8,1P,E12.4)
 5300 FORMAT(1X,'Warning: requested subprocess ',I3,' has vanishing ',
     &'cross-section.'/1X,'Process switched off!')
 5400 FORMAT(1X,'Coefficients of equation system to be solved for ',A4)
 5500 FORMAT(1X,1P,8E11.3)
 5600 FORMAT(1X,'Result for ',A4,':',7F9.4)
 5700 FORMAT(1X,'Maximum search for given coefficients'/2X,'MAX VAR ',
     &'MOD MOV   VNEW',7X,'tau',7X,'y*',8X,'cth',7X,'tau''',7X,'sigma')
 5800 FORMAT(1X,4I4,F8.4,F11.7,F9.3,F11.6,F11.7,1P,E12.4)
 5900 FORMAT(1X,'Error: no requested process has non-vanishing ',
     &'cross-section.'/1X,'Execution stopped!')
 6000 FORMAT(/1X,8('*'),1X,'RYMAXI: summary of differential ',
     &'cross-section maximum search',1X,8('*'))
 6100 FORMAT(/11X,58('=')/11X,'I',38X,'I',17X,'I'/11X,'I  ISUB  ',
     &'Subprocess name',15X,'I  Maximum value  I'/11X,'I',38X,'I',
     &17X,'I'/11X,58('=')/11X,'I',38X,'I',17X,'I')
 6200 FORMAT(11X,'I',2X,I3,3X,A28,2X,'I',2X,1P,E12.4,3X,'I')
 6300 FORMAT(11X,'I',38X,'I',17X,'I'/11X,58('='))

      RETURN
      END

C*********************************************************************

      SUBROUTINE RYPILE(MPILE)

C...Initializes multiplicity distribution and selects mutliplicity
C...of pileup events, i.e. several events occuring at the same
C...beam crossing.
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/RYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/RYINT1/MINT(400),VINT(400)
      SAVE /LUDAT1/
      SAVE /RYPARS/,/RYINT1/
      DIMENSION WTI(0:200)
      SAVE IMIN,IMAX,WTI,WTS

C...Sum of allowed cross-sections for pileup events.
      IF(MPILE.EQ.1) THEN
        VINT(131)=VINT(106)
        IF(MSTP(132).GE.2) VINT(131)=VINT(131)+VINT(104)
        IF(MSTP(132).GE.3) VINT(131)=VINT(131)+VINT(103)
        IF(MSTP(132).GE.4) VINT(131)=VINT(131)+VINT(102)
        IF(MSTP(133).LE.0) RETURN

C...Initialize multiplicity distribution at maximum.
        XNAVE=VINT(131)*PARP(131)
        IF(XNAVE.GT.120.) WRITE(MSTU(11),5000) XNAVE
        INAVE=MAX(1,MIN(200,NINT(XNAVE)))
        WTI(INAVE)=1.
        WTS=WTI(INAVE)
        WTN=WTI(INAVE)*INAVE

C...Find shape of multiplicity distribution below maximum.
        DO 100 I=INAVE-1,1,-1
        IF(MSTP(133).EQ.1) WTI(I)=WTI(I+1)*(I+1)/XNAVE
        IF(MSTP(133).GE.2) WTI(I)=WTI(I+1)*I/XNAVE
        IF(WTI(I).LT.1E-6) GOTO 110
        WTS=WTS+WTI(I)
        WTN=WTN+WTI(I)*I
  100   IMIN=I

C...Find shape of multiplicity distribution above maximum.
  110   DO 120 I=INAVE+1,200
        IF(MSTP(133).EQ.1) WTI(I)=WTI(I-1)*XNAVE/I
        IF(MSTP(133).GE.2) WTI(I)=WTI(I-1)*XNAVE/(I-1)
        IF(WTI(I).LT.1E-6) GOTO 130
        WTS=WTS+WTI(I)
        WTN=WTN+WTI(I)*I
  120   IMAX=I
  130   VINT(132)=XNAVE
        VINT(133)=WTN/WTS
        IF(MSTP(133).EQ.1.AND.IMIN.EQ.1) VINT(134)=
     &  WTS/(WTS+WTI(1)/XNAVE)
        IF(MSTP(133).EQ.1.AND.IMIN.GT.1) VINT(134)=1.
        IF(MSTP(133).GE.2) VINT(134)=XNAVE

C...Pick multiplicity of pileup events.
      ELSE
        IF(MSTP(133).LE.0) THEN
          MINT(81)=MAX(1,MSTP(134))
        ELSE
          WTR=WTS*PYR(0)
          DO 140 I=IMIN,IMAX
          MINT(81)=I
          WTR=WTR-WTI(I)
          IF(WTR.LE.0.) GOTO 150
  140     CONTINUE
  150     CONTINUE
        ENDIF
      ENDIF

C...Format statement for error message.
 5000 FORMAT(1X,'Warning: requested average number of events per bunch',
     &'crossing too large, ',1P,E12.4)

      RETURN
      END

C*********************************************************************

      SUBROUTINE RYRAND

C...Generates quantities characterizing the high-pT scattering at the
C...parton level according to the matrix elements. Chooses incoming,
C...reacting partons, their momentum fractions and one of the possible
C...subprocesses.
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/RYSUBS/MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200)
      COMMON/RYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/RYINT1/MINT(400),VINT(400)
      COMMON/RYINT2/ISET(200),KFPR(200,2),COEF(200,20),ICOL(40,4,2)
      COMMON/RYINT3/XSFX(2,-40:40),ISIG(1000,3),SIGH(1000)
      COMMON/RYINT4/WIDP(21:40,0:40),WIDE(21:40,0:40),WIDS(21:40,3)
      COMMON/RYINT5/NGEN(0:200,3),XSEC(0:200,3)
      SAVE /LUDAT1/,/LUDAT2/
      SAVE /RYSUBS/,/RYPARS/,/RYINT1/,/RYINT2/,/RYINT3/,/RYINT4/,
     &/RYINT5/
      DOUBLE PRECISION SH,SQM1,SQM2,SQM3,SQM4,SQLA12,SQLA34,
     &THTER1,THTER2,THL,THU,THM,SQMMIN,SQMMAX,SQMI,SQMJ,SQMF,
     &SQUA,QUAR,B,C,EXPTH,THARG,RATLOG,TH

C...Initial values, specifically for (first) semihard interaction.
      MINT(17)=0
      MINT(18)=0
      VINT(143)=1.
      VINT(144)=1.
      IF(MINT(82).EQ.1) NGEN(0,2)=NGEN(0,2)+1
      IF(MSUB(95).EQ.1.OR.MINT(82).GE.2) CALL RYMULT(2)
      ISUB=0
  100 MINT(51)=0

C...Choice of process type - first event of pileup.
      IF(MINT(82).EQ.1.AND.(ISUB.LE.90.OR.ISUB.GT.96)) THEN
        RSUB=XSEC(0,1)*PYR(0)
        DO 110 I=1,200
        IF(MSUB(I).NE.1) GOTO 110
        ISUB=I
        RSUB=RSUB-XSEC(I,1)
        IF(RSUB.LE.0.) GOTO 120
  110   CONTINUE
  120   IF(ISUB.EQ.95) ISUB=96

C...Choice of inclusive process type - pileup events.
      ELSEIF(MINT(82).GE.2.AND.ISUB.EQ.0) THEN
        RSUB=VINT(131)*PYR(0)
        ISUB=96
        IF(RSUB.GT.VINT(106)) ISUB=93
        IF(RSUB.GT.VINT(106)+VINT(104)) ISUB=92
        IF(RSUB.GT.VINT(106)+VINT(104)+VINT(103)) ISUB=91
      ENDIF
      IF(MINT(82).EQ.1) NGEN(0,1)=NGEN(0,1)+1
      IF(MINT(82).EQ.1) NGEN(ISUB,1)=NGEN(ISUB,1)+1
      MINT(1)=ISUB
      ISTSB=ISET(ISUB)

C...Find resonances (explicit or implicit in cross-section).
      MINT(72)=0
      KFR1=0
      IF(ISTSB.EQ.1.OR.ISTSB.EQ.3.OR.ISTSB.EQ.5) THEN
        KFR1=KFPR(ISUB,1)
      ELSEIF(ISUB.EQ.24.OR.ISUB.EQ.25.OR.ISUB.EQ.171.OR.ISUB.EQ.176)
     &THEN
        KFR1=23
      ELSEIF(ISUB.EQ.23.OR.ISUB.EQ.26.OR.ISUB.EQ.172.OR.ISUB.EQ.177)
     &THEN
        KFR1=24
      ELSEIF(ISUB.GE.71.AND.ISUB.LE.77) THEN
        KFR1=25
        IF(MSTP(46).EQ.5) THEN
          KFR1=30
          PMAS(30,1)=PARP(45)
          PMAS(30,2)=PARP(45)**3/(96.*PARU(1)*246.**2)
        ENDIF
      ENDIF
      CKMX=CKIN(2)
      IF(CKMX.LE.0.) CKMX=VINT(1)
      IF(KFR1.NE.0) THEN
        IF(CKIN(1).GT.PMAS(KFR1,1)+20.*PMAS(KFR1,2).OR.
     &  CKMX.LT.PMAS(KFR1,1)-20.*PMAS(KFR1,2)) KFR1=0
      ENDIF
      IF(KFR1.NE.0) THEN
        TAUR1=PMAS(KFR1,1)**2/VINT(2)
        GAMR1=PMAS(KFR1,1)*PMAS(KFR1,2)/VINT(2)
        MINT(72)=1
        MINT(73)=KFR1
        VINT(73)=TAUR1
        VINT(74)=GAMR1
      ENDIF
      IF(ISUB.EQ.141) THEN
        KFR2=23
        TAUR2=PMAS(KFR2,1)**2/VINT(2)
        GAMR2=PMAS(KFR2,1)*PMAS(KFR2,2)/VINT(2)
        IF(CKIN(1).GT.PMAS(KFR2,1)+20.*PMAS(KFR2,2).OR.
     &  CKMX.LT.PMAS(KFR2,1)-20.*PMAS(KFR2,2)) KFR2=0
        IF(KFR2.NE.0.AND.KFR1.NE.0) THEN
          MINT(72)=2
          MINT(74)=KFR2
          VINT(75)=TAUR2
          VINT(76)=GAMR2
        ELSEIF(KFR2.NE.0) THEN
          KFR1=KFR2
          TAUR1=TAUR2
          GAMR1=GAMR2
          MINT(72)=1
          MINT(73)=KFR1
          VINT(73)=TAUR1
          VINT(74)=GAMR1
        ENDIF
      ENDIF

C...Find product masses and minimum pT of process,
C...optionally with broadening according to a truncated Breit-Wigner.
      VINT(63)=0.
      VINT(64)=0.
      MINT(71)=0
      VINT(71)=CKIN(3)
      IF(MINT(82).GE.2) VINT(71)=0.
      VINT(80)=1.
      IF(ISTSB.EQ.2.OR.ISTSB.EQ.4) THEN
        NBW=0
        DO 130 I=1,2
        IF(KFPR(ISUB,I).EQ.0) THEN
        ELSEIF(MSTP(42).LE.0.OR.PMAS(KFPR(ISUB,I),2).LT.PARP(41)) THEN
          VINT(62+I)=PMAS(KFPR(ISUB,I),1)**2
        ELSE
          NBW=NBW+1
        ENDIF
  130   CONTINUE
        IF(NBW.GE.1) THEN
          CALL RYOFSH(4,0,KFPR(ISUB,1),KFPR(ISUB,2),0.,PQM3,PQM4)
          IF(MINT(51).EQ.1) GOTO 100
          VINT(63)=PQM3**2
          VINT(64)=PQM4**2
        ENDIF
        IF(MIN(VINT(63),VINT(64)).LT.CKIN(6)**2) MINT(71)=1
        IF(MINT(71).EQ.1) VINT(71)=MAX(CKIN(3),CKIN(5))
      ELSEIF(ISTSB.EQ.6) THEN
        CALL RYOFSH(6,0,KFPR(ISUB,1),KFPR(ISUB,2),0.,PQM3,PQM4)
        IF(MINT(51).EQ.1) GOTO 100
        VINT(63)=PQM3**2
        VINT(64)=PQM4**2
      ENDIF

C...Prepare for additional variable choices in 2 -> 3.
      IF(ISTSB.EQ.5) THEN
        VINT(201)=0.
        VINT(206)=0.
        VINT(204)=PMAS(23,1)
        IF(ISUB.EQ.124) VINT(204)=PMAS(24,1)
        VINT(209)=VINT(204)
      ENDIF

      IF(ISTSB.EQ.0) THEN
C...Double or single diffractive, or elastic scattering:
C...choose m^2 according to 1/m^2 (diffractive), constant (elastic)
        IS=INT(1.5+PYR(0))
        VINT(63)=VINT(3)**2
        VINT(64)=VINT(4)**2
        IF(ISUB.EQ.92.OR.ISUB.EQ.93) VINT(62+IS)=PARP(111)**2
        IF(ISUB.EQ.93) VINT(65-IS)=PARP(111)**2
        SH=VINT(2)
        SQM1=VINT(3)**2
        SQM2=VINT(4)**2
        SQM3=VINT(63)
        SQM4=VINT(64)
        SQLA12=(SH-SQM1-SQM2)**2-4D0*SQM1*SQM2
        SQLA34=(SH-SQM3-SQM4)**2-4D0*SQM3*SQM4
        THTER1=SQM1+SQM2+SQM3+SQM4-(SQM1-SQM2)*(SQM3-SQM4)/SH-SH
        THTER2=SQRT(MAX(0D0,SQLA12))*SQRT(MAX(0D0,SQLA34))/SH
        THL=0.5D0*(THTER1-THTER2)
        THU=0.5D0*(THTER1+THTER2)
        THM=MIN(MAX(THL,DBLE(PARP(101))),THU)
        JTMAX=0
        IF(ISUB.EQ.92.OR.ISUB.EQ.93) JTMAX=ISUB-91
        DO 140 JT=1,JTMAX
        MINT(13+3*JT-IS*(2*JT-3))=1
        SQMMIN=VINT(59+3*JT-IS*(2*JT-3))
        SQMI=VINT(8-3*JT+IS*(2*JT-3))**2
        SQMJ=VINT(3*JT-1-IS*(2*JT-3))**2
        SQMF=VINT(68-3*JT+IS*(2*JT-3))
        SQUA=0.5D0*SH/SQMI*((1D0+(SQMI-SQMJ)/SH)*THM+SQMI-SQMF-
     &  SQMJ**2/SH+(SQMI+SQMJ)*SQMF/SH+(SQMI-SQMJ)**2/SH**2*SQMF)
        QUAR=SH/SQMI*(THM*(THM+SH-SQMI-SQMJ-SQMF*(1D0-(SQMI-SQMJ)/SH))+
     &  SQMI*SQMJ-SQMJ*SQMF*(1D0+(SQMI-SQMJ-SQMF)/SH))
        SQMMAX=SQUA+SQRT(MAX(0D0,SQUA**2-QUAR))
        IF(ABS(QUAR/SQUA**2).LT.1.D-06) SQMMAX=0.5D0*QUAR/SQUA
        SQMMAX=MIN(SQMMAX,(DBLE(VINT(1))-SQRT(SQMF))**2)
        VINT(59+3*JT-IS*(2*JT-3))=SQMMIN*(SQMMAX/SQMMIN)**PYR(0)
  140   CONTINUE
C...Choose t-hat according to exp(B*t-hat+C*t-hat^2).
        SQM3=VINT(63)
        SQM4=VINT(64)
        SQLA34=(SH-SQM3-SQM4)**2-4D0*SQM3*SQM4
        THTER1=SQM1+SQM2+SQM3+SQM4-(SQM1-SQM2)*(SQM3-SQM4)/SH-SH
        THTER2=SQRT(MAX(0D0,SQLA12))*SQRT(MAX(0D0,SQLA34))/SH
        THL=0.5D0*(THTER1-THTER2)
        THU=0.5D0*(THTER1+THTER2)
        B=VINT(121)
        C=VINT(122)
        IF(ISUB.EQ.92.OR.ISUB.EQ.93) THEN
          B=0.5D0*B
          C=0.5D0*C
        ENDIF
        THM=MIN(MAX(THL,DBLE(PARP(101))),THU)
        EXPTH=0D0
        THARG=B*(THM-THU)
        IF(THARG.GT.-20D0) EXPTH=EXP(THARG)
  150   TH=THU+LOG(EXPTH+(1D0-EXPTH)*DBLE(PYR(0)))/B
        TH=MAX(THM,MIN(THU,TH))
        RATLOG=MIN((B+C*(TH+THM))*(TH-THM),(B+C*(TH+THU))*(TH-THU))
        IF(RATLOG.LT.LOG(PYR(0))) GOTO 150
        VINT(21)=1.
        VINT(22)=0.
        VINT(23)=MIN(1D0,MAX(-1D0,(2D0*TH-THTER1)/THTER2))

C...Note: in the following, by In is meant the integral over the
C...quantity multiplying coefficient cn.
C...Choose tau according to h1(tau)/tau, where
C...h1(tau) = c1 + I1/I2*c2*1/tau + I1/I3*c3*1/(tau+tau_R) +
C...I1/I4*c4*tau/((s*tau-m^2)^2+(m*Gamma)^2) +
C...I1/I5*c5*1/(tau+tau_R') +
C...I1/I6*c6*tau/((s*tau-m'^2)^2+(m'*Gamma')^2) +
C...I1/I7*c7*tau/(1.-tau), and
C...c1 + c2 + c3 + c4 + c5 + c6 + c7 = 1.
      ELSEIF(ISTSB.GE.1.AND.ISTSB.LE.6) THEN
        CALL RYKLIM(1)
        IF(MINT(51).NE.0) GOTO 100
        RTAU=PYR(0)
        MTAU=1
        IF(RTAU.GT.COEF(ISUB,1)) MTAU=2
        IF(RTAU.GT.COEF(ISUB,1)+COEF(ISUB,2)) MTAU=3
        IF(RTAU.GT.COEF(ISUB,1)+COEF(ISUB,2)+COEF(ISUB,3)) MTAU=4
        IF(RTAU.GT.COEF(ISUB,1)+COEF(ISUB,2)+COEF(ISUB,3)+COEF(ISUB,4))
     &  MTAU=5
        IF(RTAU.GT.COEF(ISUB,1)+COEF(ISUB,2)+COEF(ISUB,3)+COEF(ISUB,4)+
     &  COEF(ISUB,5)) MTAU=6
        IF(RTAU.GT.COEF(ISUB,1)+COEF(ISUB,2)+COEF(ISUB,3)+COEF(ISUB,4)+
     &  COEF(ISUB,5)+COEF(ISUB,6)) MTAU=7
        CALL RYKMAP(1,MTAU,PYR(0))

C...2 -> 3, 4 processes:
C...Choose tau' according to h4(tau,tau')/tau', where
C...h4(tau,tau') = c1 + I1/I2*c2*(1 - tau/tau')^3/tau' +
C...I1/I3*c3*1/(1 - tau'), and c1 + c2 + c3 = 1.
        IF(ISTSB.GE.3.AND.ISTSB.LE.5) THEN
          CALL RYKLIM(4)
          IF(MINT(51).NE.0) GOTO 100
          RTAUP=PYR(0)
          MTAUP=1
          IF(RTAUP.GT.COEF(ISUB,18)) MTAUP=2
          IF(RTAUP.GT.COEF(ISUB,18)+COEF(ISUB,19)) MTAUP=3
          CALL RYKMAP(4,MTAUP,PYR(0))
        ENDIF

C...Choose y* according to h2(y*), where
C...h2(y*) = I0/I1*c1*(y*-y*min) + I0/I2*c2*(y*max-y*) +
C...I0/I3*c3*1/cosh(y*) + I0/I4*c4*1/(1-exp(y*-y*max)) +
C...I0/I5*c5*1/(1-exp(-y*-y*min)), I0 = y*max-y*min,
C...and c1 + c2 + c3 + c4 + c5 = 1.
        CALL RYKLIM(2)
        IF(MINT(51).NE.0) GOTO 100
        RYST=PYR(0)
        MYST=1
        IF(RYST.GT.COEF(ISUB,8)) MYST=2
        IF(RYST.GT.COEF(ISUB,8)+COEF(ISUB,9)) MYST=3
        IF(RYST.GT.COEF(ISUB,8)+COEF(ISUB,9)+COEF(ISUB,10)) MYST=4
        IF(RYST.GT.COEF(ISUB,8)+COEF(ISUB,9)+COEF(ISUB,10)+
     &  COEF(ISUB,11)) MYST=5
        CALL RYKMAP(2,MYST,PYR(0))

C...2 -> 2 processes:
C...Choose cos(theta-hat) (cth) according to h3(cth), where
C...h3(cth) = c0 + I0/I1*c1*1/(A - cth) + I0/I2*c2*1/(A + cth) +
C...I0/I3*c3*1/(A - cth)^2 + I0/I4*c4*1/(A + cth)^2,
C...A = 1 + 2*(m3*m4/sh)^2 (= 1 for massless products),
C...and c0 + c1 + c2 + c3 + c4 = 1.
        CALL RYKLIM(3)
        IF(MINT(51).NE.0) GOTO 100
        IF(ISTSB.EQ.2.OR.ISTSB.EQ.4.OR.ISTSB.EQ.6) THEN
          RCTH=PYR(0)
          MCTH=1
          IF(RCTH.GT.COEF(ISUB,13)) MCTH=2
          IF(RCTH.GT.COEF(ISUB,13)+COEF(ISUB,14)) MCTH=3
          IF(RCTH.GT.COEF(ISUB,13)+COEF(ISUB,14)+COEF(ISUB,15)) MCTH=4
          IF(RCTH.GT.COEF(ISUB,13)+COEF(ISUB,14)+COEF(ISUB,15)+
     &    COEF(ISUB,16)) MCTH=5
          CALL RYKMAP(3,MCTH,PYR(0))
        ENDIF

C...2 -> 3 : select pT1, phi1, pT2, phi2, y3 for 3 outgoing.
        IF(ISTSB.EQ.5) THEN
          CALL RYKMAP(5,0,0.)
          IF(MINT(51).NE.0) GOTO 100
        ENDIF

C...Low-pT or multiple interactions (first semihard interaction).
      ELSEIF(ISTSB.EQ.9) THEN
        CALL RYMULT(3)
        ISUB=MINT(1)
      ENDIF

C...Choose azimuthal angle.
      VINT(24)=PARU(2)*PYR(0)

C...Check against user cuts on kinematics at parton level.
      MINT(51)=0
      IF(ISUB.LE.90.OR.ISUB.GT.100) CALL RYKLIM(0)
      IF(MINT(51).NE.0) GOTO 100
      IF(MINT(82).EQ.1.AND.MSTP(141).GE.1) THEN
        MCUT=0
        IF(MSUB(91)+MSUB(92)+MSUB(93)+MSUB(94)+MSUB(95).EQ.0)
     &  CALL RYKCUT(MCUT)
        IF(MCUT.NE.0) GOTO 100
      ENDIF

C...Calculate differential cross-section for different subprocesses.
      CALL RYSIGH(NCHN,SIGS)

C...Calculations for Monte Carlo estimate of all cross-sections.
      IF(MINT(82).EQ.1.AND.ISUB.LE.90.OR.ISUB.GE.96) THEN
        XSEC(ISUB,2)=XSEC(ISUB,2)+SIGS
      ELSEIF(MINT(82).EQ.1) THEN
        XSEC(ISUB,2)=XSEC(ISUB,2)+XSEC(ISUB,1)
      ENDIF

C...Multiple interactions: store results of cross-section calculation.
      IF(MINT(44).EQ.4.AND.MSTP(82).GE.3) THEN
        VINT(153)=SIGS
        CALL RYMULT(4)
      ENDIF

C...Weighting using estimate of maximum of differential cross-section.
C...Check that weight not negative!
      VIOL=SIGS/XSEC(ISUB,1)
      IF(MSTP(123).LE.0) THEN
        IF(VIOL.LT.-1E-3) THEN
          WRITE(MSTU(11),5000) VIOL,NGEN(0,3)+1
          WRITE(MSTU(11),5100) ISUB,VINT(21),VINT(22),VINT(23),VINT(26)
          STOP
        ENDIF
      ELSE
        IF(VIOL.LT.MIN(-1E-3,VINT(109))) THEN
          VINT(109)=VIOL
          WRITE(MSTU(11),5200) VIOL,NGEN(0,3)+1
          WRITE(MSTU(11),5100) ISUB,VINT(21),VINT(22),VINT(23),VINT(26)
        ENDIF
      ENDIF
      IF(VIOL.LT.PYR(0)) GOTO 100

C...Check for possible violation of estimated maximum of differential
C...cross-section used in weighting.
      IF(MSTP(123).LE.0) THEN
        IF(VIOL.GT.1.) THEN
          WRITE(MSTU(11),5300) VIOL,NGEN(0,3)+1
          WRITE(MSTU(11),5100) ISUB,VINT(21),VINT(22),VINT(23),VINT(26)
          STOP
        ENDIF
      ELSEIF(MSTP(123).EQ.1) THEN
        IF(VIOL.GT.VINT(108)) THEN
          VINT(108)=VIOL
          IF(VIOL.GT.1.) THEN
            WRITE(MSTU(11),5400) VIOL,NGEN(0,3)+1
            WRITE(MSTU(11),5100) ISUB,VINT(21),VINT(22),VINT(23),
     &      VINT(26)
          ENDIF
        ENDIF
      ELSEIF(VIOL.GT.VINT(108)) THEN
        VINT(108)=VIOL
        IF(VIOL.GT.1.) THEN
          XDIF=XSEC(ISUB,1)*(VIOL-1.)
          XSEC(ISUB,1)=XSEC(ISUB,1)+XDIF
          IF(MSUB(ISUB).EQ.1.AND.(ISUB.LE.90.OR.ISUB.GT.96))
     &    XSEC(0,1)=XSEC(0,1)+XDIF
          WRITE(MSTU(11),5400) VIOL,NGEN(0,3)+1
          WRITE(MSTU(11),5100) ISUB,VINT(21),VINT(22),VINT(23),VINT(26)
          IF(ISUB.LE.9) THEN
            WRITE(MSTU(11),5500) ISUB,XSEC(ISUB,1)
          ELSEIF(ISUB.LE.99) THEN
            WRITE(MSTU(11),5600) ISUB,XSEC(ISUB,1)
          ELSE
            WRITE(MSTU(11),5700) ISUB,XSEC(ISUB,1)
          ENDIF
          VINT(108)=1.
        ENDIF
      ENDIF

C...Multiple interactions: choose impact parameter.
      VINT(148)=1.
      IF(MINT(44).EQ.4.AND.(ISUB.LE.90.OR.ISUB.GE.96).AND.MSTP(82).GE.3)
     &THEN
        CALL RYMULT(5)
        IF(VINT(150).LT.PYR(0)) GOTO 100
      ENDIF
      IF(MINT(82).EQ.1.AND.MSUB(95).EQ.1) THEN
        IF(ISUB.LE.90.OR.ISUB.GE.95) NGEN(95,1)=NGEN(95,1)+1
        IF(ISUB.LE.90.OR.ISUB.GE.96) NGEN(96,2)=NGEN(96,2)+1
      ENDIF
      IF(ISUB.LE.90.OR.ISUB.GE.96) MINT(31)=MINT(31)+1

C...Choose flavour of reacting partons (and subprocess).
      RSIGS=SIGS*PYR(0)
      QT2=VINT(48)
      RQQBAR=PARP(87)*(1.-(QT2/(QT2+(PARP(88)*PARP(82))**2))**2)
      IF(ISUB.NE.95.AND.(ISUB.NE.96.OR.MSTP(82).LE.1.OR.
     &PYR(0).GT.RQQBAR)) THEN
        DO 160 ICHN=1,NCHN
        KFL1=ISIG(ICHN,1)
        KFL2=ISIG(ICHN,2)
        MINT(2)=ISIG(ICHN,3)
        RSIGS=RSIGS-SIGH(ICHN)
        IF(RSIGS.LE.0.) GOTO 170
  160   CONTINUE

C...Multiple interactions: choose qq~ preferentially at small pT.
      ELSEIF(ISUB.EQ.96) THEN
        CALL RYSPLI(MINT(11),21,KFL1,KFLDUM)
        CALL RYSPLI(MINT(12),21,KFL2,KFLDUM)
        MINT(1)=11
        MINT(2)=1
        IF(KFL1.EQ.KFL2.AND.PYR(0).LT.0.5) MINT(2)=2

C...Low-pT: choose string drawing configuration.
      ELSE
        KFL1=21
        KFL2=21
        RSIGS=6.*PYR(0)
        MINT(2)=1
        IF(RSIGS.GT.1.) MINT(2)=2
        IF(RSIGS.GT.2.) MINT(2)=3
      ENDIF

C...Reassign QCD process. Partons before initial state radiation.
  170 IF(MINT(2).GT.10) THEN
        MINT(1)=MINT(2)/10
        MINT(2)=MOD(MINT(2),10)
      ENDIF
      IF(MINT(82).EQ.1.AND.MSTP(111).GE.0) NGEN(MINT(1),2)=
     &NGEN(MINT(1),2)+1
      MINT(15)=KFL1
      MINT(16)=KFL2
      MINT(13)=MINT(15)
      MINT(14)=MINT(16)
      VINT(141)=VINT(41)
      VINT(142)=VINT(42)
      VINT(151)=0.
      VINT(152)=0.

C...Format statements for differential cross-section maximum violations.
 5000 FORMAT(1X,'Error: negative cross-section fraction',1P,E11.3,1X,
     &'in event',1X,I7,'.'/1X,'Execution stopped!')
 5100 FORMAT(1X,'ISUB = ',I3,'; Point of violation:'/1X,'tau =',1P,
     &E11.3,', y* =',E11.3,', cthe = ',0P,F11.7,', tau'' =',1P,E11.3)
 5200 FORMAT(1X,'Warning: negative cross-section fraction',1P,E11.3,1X,
     &'in event',1X,I7)
 5300 FORMAT(1X,'Error: maximum violated by',1P,E11.3,1X,
     &'in event',1X,I7,'.'/1X,'Execution stopped!')
 5400 FORMAT(1X,'Warning: maximum violated by',1P,E11.3,1X,
     &'in event',1X,I7)
 5500 FORMAT(1X,'XSEC(',I1,',1) increased to',1P,E11.3)
 5600 FORMAT(1X,'XSEC(',I2,',1) increased to',1P,E11.3)
 5700 FORMAT(1X,'XSEC(',I3,',1) increased to',1P,E11.3)

      RETURN
      END

C*********************************************************************

      SUBROUTINE RYSCAT
      PARAMETER (KSZJ=4000)

C...Finds outgoing flavours and event type; sets up the kinematics
C...and colour flow of the hard scattering.
      
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      integer N,NPAD,K
      double precision P,V
      SAVE /PYJETS/


      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/LUDAT3/MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5)
      COMMON/RYSUBS/MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200)
      COMMON/RYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/RYINT1/MINT(400),VINT(400)
      COMMON/RYINT2/ISET(200),KFPR(200,2),COEF(200,20),ICOL(40,4,2)
      COMMON/RYINT3/XSFX(2,-40:40),ISIG(1000,3),SIGH(1000)
      COMMON/RYINT4/WIDP(21:40,0:40),WIDE(21:40,0:40),WIDS(21:40,3)
      COMMON/RYINT5/NGEN(0:200,3),XSEC(0:200,3)
      SAVE /LUDAT1/,/LUDAT2/,/LUDAT3/
      SAVE /RYSUBS/,/RYPARS/,/RYINT1/,/RYINT2/,/RYINT3/,/RYINT4/,
     &/RYINT5/
      DIMENSION WDTP(0:40),WDTE(0:40,0:5),PMQ(2),Z(2),CTHE(2),PHI(2)

C...Convert H' or A process into equivalent H one.
      ISUB=MINT(1)
      IHIGG=1
      KFHIGG=25
      IF((ISUB.GE.151.AND.ISUB.LE.160).OR.(ISUB.GE.171.AND.
     &ISUB.LE.180)) THEN
        IHIGG=2
        IF(MOD(ISUB-1,10).GE.5) IHIGG=3
        KFHIGG=33+IHIGG
        IF(ISUB.EQ.151.OR.ISUB.EQ.156) ISUB=3
        IF(ISUB.EQ.152.OR.ISUB.EQ.157) ISUB=102
        IF(ISUB.EQ.153.OR.ISUB.EQ.158) ISUB=103
        IF(ISUB.EQ.171.OR.ISUB.EQ.176) ISUB=24
        IF(ISUB.EQ.172.OR.ISUB.EQ.177) ISUB=26
        IF(ISUB.EQ.173.OR.ISUB.EQ.178) ISUB=123
        IF(ISUB.EQ.174.OR.ISUB.EQ.179) ISUB=124
      ENDIF

C...Choice of subprocess, number of documentation lines.
      IDOC=6+ISET(ISUB)
      IF(ISUB.EQ.95) IDOC=8
      IF(ISET(ISUB).EQ.5.OR.ISET(ISUB).EQ.6) IDOC=9
      MINT(3)=IDOC-6
      IF(IDOC.GE.9.AND.ISET(ISUB).LE.4) IDOC=IDOC+2
      MINT(4)=IDOC
      IPU1=MINT(84)+1
      IPU2=MINT(84)+2
      IPU3=MINT(84)+3
      IPU4=MINT(84)+4
      IPU5=MINT(84)+5
      IPU6=MINT(84)+6

C...Reset K, P and V vectors. Store incoming particles.
      DO 100 JT=1,MSTP(126)+10
      I=MINT(83)+JT
      DO 100 J=1,5
      K(I,J)=0
      P(I,J)=0.
  100 V(I,J)=0.
      DO 110 JT=1,2
      I=MINT(83)+JT
      K(I,1)=21
      K(I,2)=MINT(10+JT)
      P(I,1)=0.
      P(I,2)=0.
      P(I,5)=VINT(2+JT)
      P(I,3)=VINT(5)*(-1)**(JT+1)
  110 P(I,4)=SQRT(P(I,3)**2+P(I,5)**2)
      MINT(6)=2
      KFRES=0

C...Store incoming partons in their CM-frame.
      SH=VINT(44)
      SHR=SQRT(SH)
      SHP=VINT(26)*VINT(2)
      SHPR=SQRT(SHP)
      SHUSER=SHR
      IF(ISET(ISUB).GE.3.AND.ISET(ISUB).LE.5) SHUSER=SHPR
      DO 120 JT=1,2
      I=MINT(84)+JT
      K(I,1)=14
      K(I,2)=MINT(14+JT)
      K(I,3)=MINT(83)+2+JT
      P(I,3)=0.5*SHUSER*(-1.)**(JT-1)
  120 P(I,4)=0.5*SHUSER

C...Copy incoming partons to documentation lines.
      DO 130 JT=1,2
      I1=MINT(83)+4+JT
      I2=MINT(84)+JT
      K(I1,1)=21
      K(I1,2)=K(I2,2)
      K(I1,3)=I1-2
      DO 130 J=1,5
  130 P(I1,J)=P(I2,J)

C...Choose new quark/lepton flavour for relevant annihilation graphs.
      IF(ISUB.EQ.12.OR.ISUB.EQ.53.OR.ISUB.EQ.54.OR.ISUB.EQ.58) THEN
        IGLGA=21
        IF(ISUB.EQ.58) IGLGA=22
        CALL RYWIDT(IGLGA,SH,WDTP,WDTE)
  140   RKFL=(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))*PYR(0)
        DO 150 I=1,MDCY(IGLGA,3)
        KFLF=KFDP(I+MDCY(IGLGA,2)-1,1)
        RKFL=RKFL-(WDTE(I,1)+WDTE(I,2)+WDTE(I,4))
        IF(RKFL.LE.0.) GOTO 160
  150   CONTINUE
  160   CONTINUE
        IF(ISUB.EQ.54) THEN
          IF((KCHG(IABS(KFLF),1)/2.)**2.LT.PYR(0)) GOTO 140
        ELSEIF(ISUB.EQ.58) THEN
          IF((KCHG(IABS(KFLF),1)/3.)**2.LT.PYR(0)) GOTO 140
        ENDIF
      ENDIF

C...Final state flavours and colour flow: default values.
      JS=1
      MINT(21)=MINT(15)
      MINT(22)=MINT(16)
      MINT(23)=0
      MINT(24)=0
      KCC=20
      KCS=ISIGN(1,MINT(15))

      IF(ISUB.LE.10) THEN
      IF(ISUB.EQ.1) THEN
C...f + f~ -> gamma*/Z0.
        KFRES=23

      ELSEIF(ISUB.EQ.2) THEN
C...f + f~' -> W+/- .
        KCH1=KCHG(IABS(MINT(15)),1)*ISIGN(1,MINT(15))
        KCH2=KCHG(IABS(MINT(16)),1)*ISIGN(1,MINT(16))
        KFRES=ISIGN(24,KCH1+KCH2)

      ELSEIF(ISUB.EQ.3) THEN
C...f + f~ -> H0 (or H'0, or A0).
        KFRES=KFHIGG

      ELSEIF(ISUB.EQ.4) THEN
C...gamma + W+/- -> W+/-.

      ELSEIF(ISUB.EQ.5) THEN
C...Z0 + Z0 -> H0.
        XH=SH/SHP
        MINT(21)=MINT(15)
        MINT(22)=MINT(16)
        PMQ(1)=ULMASS(MINT(21))
        PMQ(2)=ULMASS(MINT(22))
  170   JT=INT(1.5+PYR(0))
        ZMIN=2.*PMQ(JT)/SHPR
        ZMAX=1.-PMQ(3-JT)/SHPR-(SH-PMQ(JT)**2)/(SHPR*(SHPR-PMQ(3-JT)))
        ZMAX=MIN(1.-XH,ZMAX)
        Z(JT)=ZMIN+(ZMAX-ZMIN)*PYR(0)
        IF(-1.+(1.+XH)/(1.-Z(JT))-XH/(1.-Z(JT))**2.LT.
     &  (1.-XH)**2/(4.*XH)*PYR(0)) GOTO 170
        SQC1=1.-4.*PMQ(JT)**2/(Z(JT)**2*SHP)
        IF(SQC1.LT.1.E-8) GOTO 170
        C1=SQRT(SQC1)
        C2=1.+2.*(PMAS(23,1)**2-PMQ(JT)**2)/(Z(JT)*SHP)
        CTHE(JT)=(C2-(C2**2-C1**2)/(C2+(2.*PYR(0)-1.)*C1))/C1
        CTHE(JT)=MIN(1.,MAX(-1.,CTHE(JT)))
        Z(3-JT)=1.-XH/(1.-Z(JT))
        SQC1=1.-4.*PMQ(3-JT)**2/(Z(3-JT)**2*SHP)
        IF(SQC1.LT.1.E-8) GOTO 170
        C1=SQRT(SQC1)
        C2=1.+2.*(PMAS(23,1)**2-PMQ(3-JT)**2)/(Z(3-JT)*SHP)
        CTHE(3-JT)=(C2-(C2**2-C1**2)/(C2+(2.*PYR(0)-1.)*C1))/C1
        CTHE(3-JT)=MIN(1.,MAX(-1.,CTHE(3-JT)))
        PHIR=PARU(2)*PYR(0)
        CPHI=COS(PHIR)
        ANG=CTHE(1)*CTHE(2)-SQRT(1.-CTHE(1)**2)*SQRT(1.-CTHE(2)**2)*CPHI
        Z1=2.-Z(JT)
        Z2=ANG*SQRT(Z(JT)**2-4.*PMQ(JT)**2/SHP)
        Z3=1.-Z(JT)-XH+(PMQ(1)**2+PMQ(2)**2)/SHP
        Z(3-JT)=2./(Z1**2-Z2**2)*(Z1*Z3+Z2*SQRT(Z3**2-(Z1**2-Z2**2)*
     &  PMQ(3-JT)**2/SHP))
        ZMIN=2.*PMQ(3-JT)/SHPR
        ZMAX=1.-PMQ(JT)/SHPR-(SH-PMQ(3-JT)**2)/(SHPR*(SHPR-PMQ(JT)))
        ZMAX=MIN(1.-XH,ZMAX)
        IF(Z(3-JT).LT.ZMIN.OR.Z(3-JT).GT.ZMAX) GOTO 170
        KCC=22
        KFRES=25

      ELSEIF(ISUB.EQ.6) THEN
C...Z0 + W+/- -> W+/-.

      ELSEIF(ISUB.EQ.7) THEN
C...W+ + W- -> Z0.

      ELSEIF(ISUB.EQ.8) THEN
C...W+ + W- -> H0.
        XH=SH/SHP
  180   DO 200 JT=1,2
        I=MINT(14+JT)
        IA=IABS(I)
        IF(IA.LE.10) THEN
          RVCKM=VINT(180+I)*PYR(0)
          DO 190 J=1,MSTP(1)
          IB=2*J-1+MOD(IA,2)
          IPM=(5-ISIGN(1,I))/2
          IDC=J+MDCY(IA,2)+2
          IF(MDME(IDC,1).NE.1.AND.MDME(IDC,1).NE.IPM) GOTO 190
          MINT(20+JT)=ISIGN(IB,I)
          RVCKM=RVCKM-VCKM((IA+1)/2,(IB+1)/2)
          IF(RVCKM.LE.0.) GOTO 200
  190     CONTINUE
        ELSE
          IB=2*((IA+1)/2)-1+MOD(IA,2)
          MINT(20+JT)=ISIGN(IB,I)
        ENDIF
  200   PMQ(JT)=ULMASS(MINT(20+JT))
        JT=INT(1.5+PYR(0))
        ZMIN=2.*PMQ(JT)/SHPR
        ZMAX=1.-PMQ(3-JT)/SHPR-(SH-PMQ(JT)**2)/(SHPR*(SHPR-PMQ(3-JT)))
        ZMAX=MIN(1.-XH,ZMAX)
        IF(ZMIN.GE.ZMAX) GOTO 180
        Z(JT)=ZMIN+(ZMAX-ZMIN)*PYR(0)
        IF(-1.+(1.+XH)/(1.-Z(JT))-XH/(1.-Z(JT))**2.LT.
     &  (1.-XH)**2/(4.*XH)*PYR(0)) GOTO 180
        SQC1=1.-4.*PMQ(JT)**2/(Z(JT)**2*SHP)
        IF(SQC1.LT.1.E-8) GOTO 180
        C1=SQRT(SQC1)
        C2=1.+2.*(PMAS(24,1)**2-PMQ(JT)**2)/(Z(JT)*SHP)
        CTHE(JT)=(C2-(C2**2-C1**2)/(C2+(2.*PYR(0)-1.)*C1))/C1
        CTHE(JT)=MIN(1.,MAX(-1.,CTHE(JT)))
        Z(3-JT)=1.-XH/(1.-Z(JT))
        SQC1=1.-4.*PMQ(3-JT)**2/(Z(3-JT)**2*SHP)
        IF(SQC1.LT.1.E-8) GOTO 180
        C1=SQRT(SQC1)
        C2=1.+2.*(PMAS(24,1)**2-PMQ(3-JT)**2)/(Z(3-JT)*SHP)
        CTHE(3-JT)=(C2-(C2**2-C1**2)/(C2+(2.*PYR(0)-1.)*C1))/C1
        CTHE(3-JT)=MIN(1.,MAX(-1.,CTHE(3-JT)))
        PHIR=PARU(2)*PYR(0)
        CPHI=COS(PHIR)
        ANG=CTHE(1)*CTHE(2)-SQRT(1.-CTHE(1)**2)*SQRT(1.-CTHE(2)**2)*CPHI
        Z1=2.-Z(JT)
        Z2=ANG*SQRT(Z(JT)**2-4.*PMQ(JT)**2/SHP)
        Z3=1.-Z(JT)-XH+(PMQ(1)**2+PMQ(2)**2)/SHP
        Z(3-JT)=2./(Z1**2-Z2**2)*(Z1*Z3+Z2*SQRT(Z3**2-(Z1**2-Z2**2)*
     &  PMQ(3-JT)**2/SHP))
        ZMIN=2.*PMQ(3-JT)/SHPR
        ZMAX=1.-PMQ(JT)/SHPR-(SH-PMQ(3-JT)**2)/(SHPR*(SHPR-PMQ(JT)))
        ZMAX=MIN(1.-XH,ZMAX)
        IF(Z(3-JT).LT.ZMIN.OR.Z(3-JT).GT.ZMAX) GOTO 180
        KCC=22
        KFRES=25
      ENDIF

      ELSEIF(ISUB.LE.20) THEN
      IF(ISUB.EQ.11) THEN
C...f + f' -> f + f'; th = (p(f)-p(f))**2.
        IF(MINT(2).LE.4) THEN
          KCC=MINT(2)
          IF(MINT(15)*MINT(16).LT.0) KCC=KCC+2
        ELSEIF(MINT(2).EQ.5) THEN
          KCC=22
        ELSE
C...W exchange: need to mix flavours according to CKM matrix.
          DO 220 JT=1,2
          I=MINT(14+JT)
          IA=IABS(I)
          IF(IA.LE.10) THEN
            RVCKM=VINT(180+I)*PYR(0)
            DO 210 J=1,MSTP(1)
            IB=2*J-1+MOD(IA,2)
            IPM=(5-ISIGN(1,I))/2
            IDC=J+MDCY(IA,2)+2
            IF(MDME(IDC,1).NE.1.AND.MDME(IDC,1).NE.IPM) GOTO 210
            MINT(20+JT)=ISIGN(IB,I)
            RVCKM=RVCKM-VCKM((IA+1)/2,(IB+1)/2)
            IF(RVCKM.LE.0.) GOTO 220
  210       CONTINUE
          ELSE
            IB=2*((IA+1)/2)-1+MOD(IA,2)
            MINT(20+JT)=ISIGN(IB,I)
          ENDIF
  220     CONTINUE
          KCC=22
        ENDIF

      ELSEIF(ISUB.EQ.12) THEN
C...f + f~ -> f' + f~'; th = (p(f)-p(f'))**2.
        MINT(21)=ISIGN(KFLF,MINT(15))
        MINT(22)=-MINT(21)
        KCC=4

      ELSEIF(ISUB.EQ.13) THEN
C...f + f~ -> g + g; th arbitrary.
        MINT(21)=21
        MINT(22)=21
        KCC=MINT(2)+4

      ELSEIF(ISUB.EQ.14) THEN
C...f + f~ -> g + gamma; th arbitrary.
        IF(PYR(0).GT.0.5) JS=2
        MINT(20+JS)=21
        MINT(23-JS)=22
        KCC=17+JS

      ELSEIF(ISUB.EQ.15) THEN
C...f + f~ -> g + Z0; th arbitrary.
        IF(PYR(0).GT.0.5) JS=2
        MINT(20+JS)=21
        MINT(23-JS)=23
        KCC=17+JS

      ELSEIF(ISUB.EQ.16) THEN
C...f + f~' -> g + W+/-; th = (p(f)-p(W-))**2 or (p(f~')-p(W+))**2.
        KCH1=KCHG(IABS(MINT(15)),1)*ISIGN(1,MINT(15))
        KCH2=KCHG(IABS(MINT(16)),1)*ISIGN(1,MINT(16))
        IF(MINT(15)*(KCH1+KCH2).LT.0) JS=2
        MINT(20+JS)=21
        MINT(23-JS)=ISIGN(24,KCH1+KCH2)
        KCC=17+JS

      ELSEIF(ISUB.EQ.17) THEN
C...f + f~ -> g + H0; th arbitrary.
        IF(PYR(0).GT.0.5) JS=2
        MINT(20+JS)=21
        MINT(23-JS)=25
        KCC=17+JS

      ELSEIF(ISUB.EQ.18) THEN
C...f + f~ -> gamma + gamma; th arbitrary.
        MINT(21)=22
        MINT(22)=22

      ELSEIF(ISUB.EQ.19) THEN
C...f + f~ -> gamma + Z0; th arbitrary.
        IF(PYR(0).GT.0.5) JS=2
        MINT(20+JS)=22
        MINT(23-JS)=23

      ELSEIF(ISUB.EQ.20) THEN
C...f + f~' -> gamma + W+/-; th = (p(f)-p(W-))**2 or (p(f~')-p(W+))**2.
        KCH1=KCHG(IABS(MINT(15)),1)*ISIGN(1,MINT(15))
        KCH2=KCHG(IABS(MINT(16)),1)*ISIGN(1,MINT(16))
        IF(MINT(15)*(KCH1+KCH2).LT.0) JS=2
        MINT(20+JS)=22
        MINT(23-JS)=ISIGN(24,KCH1+KCH2)
      ENDIF

      ELSEIF(ISUB.LE.30) THEN
      IF(ISUB.EQ.21) THEN
C...f + f~ -> gamma + H0; th arbitrary.
        IF(PYR(0).GT.0.5) JS=2
        MINT(20+JS)=22
        MINT(23-JS)=25

      ELSEIF(ISUB.EQ.22) THEN
C...f + f~ -> Z0 + Z0; th arbitrary.
        MINT(21)=23
        MINT(22)=23

      ELSEIF(ISUB.EQ.23) THEN
C...f + f~' -> Z0 + W+/-; th = (p(f)-p(W-))**2 or (p(f~')-p(W+))**2.
        KCH1=KCHG(IABS(MINT(15)),1)*ISIGN(1,MINT(15))
        KCH2=KCHG(IABS(MINT(16)),1)*ISIGN(1,MINT(16))
        IF(MINT(15)*(KCH1+KCH2).LT.0) JS=2
        MINT(20+JS)=23
        MINT(23-JS)=ISIGN(24,KCH1+KCH2)

      ELSEIF(ISUB.EQ.24) THEN
C...f + f~ -> Z0 + H0 (or H'0, or A0); th arbitrary.
        IF(PYR(0).GT.0.5) JS=2
        MINT(20+JS)=23
        MINT(23-JS)=KFHIGG

      ELSEIF(ISUB.EQ.25) THEN
C...f + f~ -> W+ + W-; th = (p(f)-p(W-))**2.
        MINT(21)=-ISIGN(24,MINT(15))
        MINT(22)=-MINT(21)

      ELSEIF(ISUB.EQ.26) THEN
C...f + f~' -> W+/- + H0 (or H'0, or A0);
C...th = (p(f)-p(W-))**2 or (p(f~')-p(W+))**2.
        KCH1=KCHG(IABS(MINT(15)),1)*ISIGN(1,MINT(15))
        KCH2=KCHG(IABS(MINT(16)),1)*ISIGN(1,MINT(16))
        IF(MINT(15)*(KCH1+KCH2).GT.0) JS=2
        MINT(20+JS)=ISIGN(24,KCH1+KCH2)
        MINT(23-JS)=KFHIGG

      ELSEIF(ISUB.EQ.27) THEN
C...f + f~ -> H0 + H0.

      ELSEIF(ISUB.EQ.28) THEN
C...f + g -> f + g; th = (p(f)-p(f))**2.
        KCC=MINT(2)+6
        IF(MINT(15).EQ.21) KCC=KCC+2
        IF(MINT(15).NE.21) KCS=ISIGN(1,MINT(15))
        IF(MINT(16).NE.21) KCS=ISIGN(1,MINT(16))

      ELSEIF(ISUB.EQ.29) THEN
C...f + g -> f + gamma; th = (p(f)-p(f))**2.
        IF(MINT(15).EQ.21) JS=2
        MINT(23-JS)=22
        KCC=15+JS
        KCS=ISIGN(1,MINT(14+JS))

      ELSEIF(ISUB.EQ.30) THEN
C...f + g -> f + Z0; th = (p(f)-p(f))**2.
        IF(MINT(15).EQ.21) JS=2
        MINT(23-JS)=23
        KCC=15+JS
        KCS=ISIGN(1,MINT(14+JS))
      ENDIF

      ELSEIF(ISUB.LE.40) THEN
      IF(ISUB.EQ.31) THEN
C...f + g -> f' + W+/-; th = (p(f)-p(f'))**2; choose flavour f'.
        IF(MINT(15).EQ.21) JS=2
        I=MINT(14+JS)
        IA=IABS(I)
        MINT(23-JS)=ISIGN(24,KCHG(IA,1)*I)
        RVCKM=VINT(180+I)*PYR(0)
        DO 230 J=1,MSTP(1)
        IB=2*J-1+MOD(IA,2)
        IPM=(5-ISIGN(1,I))/2
        IDC=J+MDCY(IA,2)+2
        IF(MDME(IDC,1).NE.1.AND.MDME(IDC,1).NE.IPM) GOTO 230
        MINT(20+JS)=ISIGN(IB,I)
        RVCKM=RVCKM-VCKM((IA+1)/2,(IB+1)/2)
        IF(RVCKM.LE.0.) GOTO 240
  230   CONTINUE
  240   KCC=15+JS
        KCS=ISIGN(1,MINT(14+JS))

      ELSEIF(ISUB.EQ.32) THEN
C...f + g -> f + H0; th = (p(f)-p(f))**2.
        IF(MINT(15).EQ.21) JS=2
        MINT(23-JS)=25
        KCC=15+JS
        KCS=ISIGN(1,MINT(14+JS))

      ELSEIF(ISUB.EQ.33) THEN
C...f + gamma -> f + g; th=(p(f)-p(f))**2.
        IF(MINT(15).EQ.22) JS=2
        MINT(23-JS)=21
        KCC=24+JS
        KCS=ISIGN(1,MINT(14+JS))

      ELSEIF(ISUB.EQ.34) THEN
C...f + gamma -> f + gamma; th=(p(f)-p(f))**2.
        IF(MINT(15).EQ.22) JS=2
        KCC=22
        KCS=ISIGN(1,MINT(14+JS))

      ELSEIF(ISUB.EQ.35) THEN
C...f + gamma -> f + Z0; th=(p(f)-p(f))**2.
        IF(MINT(15).EQ.22) JS=2
        MINT(23-JS)=23
        KCC=22

      ELSEIF(ISUB.EQ.36) THEN
C...f + gamma -> f' + W+/-; th=(p(f)-p(f'))**2.
        IF(MINT(15).EQ.22) JS=2
        I=MINT(14+JS)
        IA=IABS(I)
        MINT(23-JS)=ISIGN(24,KCHG(IA,1)*I)
        IF(IA.LE.10) THEN
          RVCKM=VINT(180+I)*PYR(0)
          DO 250 J=1,MSTP(1)
          IB=2*J-1+MOD(IA,2)
          IPM=(5-ISIGN(1,I))/2
          IDC=J+MDCY(IA,2)+2
          IF(MDME(IDC,1).NE.1.AND.MDME(IDC,1).NE.IPM) GOTO 250
          MINT(20+JS)=ISIGN(IB,I)
          RVCKM=RVCKM-VCKM((IA+1)/2,(IB+1)/2)
          IF(RVCKM.LE.0.) GOTO 260
  250     CONTINUE
        ELSE
          IB=2*((IA+1)/2)-1+MOD(IA,2)
          MINT(20+JS)=ISIGN(IB,I)
        ENDIF
  260   KCC=22

      ELSEIF(ISUB.EQ.37) THEN
C...f + gamma -> f + H0.

      ELSEIF(ISUB.EQ.38) THEN
C...f + Z0 -> f + g.

      ELSEIF(ISUB.EQ.39) THEN
C...f + Z0 -> f + gamma.

      ELSEIF(ISUB.EQ.40) THEN
C...f + Z0 -> f + Z0.
      ENDIF

      ELSEIF(ISUB.LE.50) THEN
      IF(ISUB.EQ.41) THEN
C...f + Z0 -> f' + W+/-.

      ELSEIF(ISUB.EQ.42) THEN
C...f + Z0 -> f + H0.

      ELSEIF(ISUB.EQ.43) THEN
C...f + W+/- -> f' + g.

      ELSEIF(ISUB.EQ.44) THEN
C...f + W+/- -> f' + gamma.

      ELSEIF(ISUB.EQ.45) THEN
C...f + W+/- -> f' + Z0.

      ELSEIF(ISUB.EQ.46) THEN
C...f + W+/- -> f' + W+/-.

      ELSEIF(ISUB.EQ.47) THEN
C...f + W+/- -> f' + H0.

      ELSEIF(ISUB.EQ.48) THEN
C...f + H0 -> f + g.

      ELSEIF(ISUB.EQ.49) THEN
C...f + H0 -> f + gamma.

      ELSEIF(ISUB.EQ.50) THEN
C...f + H0 -> f + Z0.
      ENDIF

      ELSEIF(ISUB.LE.60) THEN
      IF(ISUB.EQ.51) THEN
C...f + H0 -> f' + W+/-.

      ELSEIF(ISUB.EQ.52) THEN
C...f + H0 -> f + H0.

      ELSEIF(ISUB.EQ.53) THEN
C...g + g -> f + f~; th arbitrary.
        KCS=(-1)**INT(1.5+PYR(0))
        MINT(21)=ISIGN(KFLF,KCS)
        MINT(22)=-MINT(21)
        KCC=MINT(2)+10

      ELSEIF(ISUB.EQ.54) THEN
C...g + gamma -> f + f~; th arbitrary.
        KCS=(-1)**INT(1.5+PYR(0))
        MINT(21)=ISIGN(KFLF,KCS)
        MINT(22)=-MINT(21)
        KCC=27
        IF(MINT(16).EQ.21) KCC=28

      ELSEIF(ISUB.EQ.55) THEN
C...g + Z0 -> f + f~.

      ELSEIF(ISUB.EQ.56) THEN
C...g + W+/- -> f + f~'.

      ELSEIF(ISUB.EQ.57) THEN
C...g + H0 -> f + f~.

      ELSEIF(ISUB.EQ.58) THEN
C...gamma + gamma -> f + f~; th arbitrary.
        KCS=(-1)**INT(1.5+PYR(0))
        MINT(21)=ISIGN(KFLF,KCS)
        MINT(22)=-MINT(21)
        KCC=21

      ELSEIF(ISUB.EQ.59) THEN
C...gamma + Z0 -> f + f~.

      ELSEIF(ISUB.EQ.60) THEN
C...gamma + W+/- -> f + f~'.
      ENDIF

      ELSEIF(ISUB.LE.70) THEN
      IF(ISUB.EQ.61) THEN
C...gamma + H0 -> f + f~.

      ELSEIF(ISUB.EQ.62) THEN
C...Z0 + Z0 -> f + f~.

      ELSEIF(ISUB.EQ.63) THEN
C...Z0 + W+/- -> f + f~'.

      ELSEIF(ISUB.EQ.64) THEN
C...Z0 + H0 -> f + f~.

      ELSEIF(ISUB.EQ.65) THEN
C...W+ + W- -> f + f~.

      ELSEIF(ISUB.EQ.66) THEN
C...W+/- + H0 -> f + f~'.

      ELSEIF(ISUB.EQ.67) THEN
C...H0 + H0 -> f + f~.

      ELSEIF(ISUB.EQ.68) THEN
C...g + g -> g + g; th arbitrary.
        KCC=MINT(2)+12
        KCS=(-1)**INT(1.5+PYR(0))

      ELSEIF(ISUB.EQ.69) THEN
C...gamma + gamma -> W+ + W-; th arbitrary.
        MINT(21)=24
        MINT(22)=-24
        KCC=21

      ELSEIF(ISUB.EQ.70) THEN
C...gamma + W+/- -> Z0 + W+/-; th=(p(W)-p(W))**2.
        IF(MINT(15).EQ.22) MINT(21)=23
        IF(MINT(16).EQ.22) MINT(22)=23
        KCC=21
      ENDIF

      ELSEIF(ISUB.LE.80) THEN
      IF(ISUB.EQ.71.OR.ISUB.EQ.72) THEN
C...Z0 + Z0 -> Z0 + Z0; Z0 + Z0 -> W+ + W-.
        XH=SH/SHP
        MINT(21)=MINT(15)
        MINT(22)=MINT(16)
        PMQ(1)=ULMASS(MINT(21))
        PMQ(2)=ULMASS(MINT(22))
  270   JT=INT(1.5+PYR(0))
        ZMIN=2.*PMQ(JT)/SHPR
        ZMAX=1.-PMQ(3-JT)/SHPR-(SH-PMQ(JT)**2)/(SHPR*(SHPR-PMQ(3-JT)))
        ZMAX=MIN(1.-XH,ZMAX)
        Z(JT)=ZMIN+(ZMAX-ZMIN)*PYR(0)
        IF(-1.+(1.+XH)/(1.-Z(JT))-XH/(1.-Z(JT))**2.LT.
     &  (1.-XH)**2/(4.*XH)*PYR(0)) GOTO 270
        SQC1=1.-4.*PMQ(JT)**2/(Z(JT)**2*SHP)
        IF(SQC1.LT.1.E-8) GOTO 270
        C1=SQRT(SQC1)
        C2=1.+2.*(PMAS(23,1)**2-PMQ(JT)**2)/(Z(JT)*SHP)
        CTHE(JT)=(C2-(C2**2-C1**2)/(C2+(2.*PYR(0)-1.)*C1))/C1
        CTHE(JT)=MIN(1.,MAX(-1.,CTHE(JT)))
        Z(3-JT)=1.-XH/(1.-Z(JT))
        SQC1=1.-4.*PMQ(3-JT)**2/(Z(3-JT)**2*SHP)
        IF(SQC1.LT.1.E-8) GOTO 270
        C1=SQRT(SQC1)
        C2=1.+2.*(PMAS(23,1)**2-PMQ(3-JT)**2)/(Z(3-JT)*SHP)
        CTHE(3-JT)=(C2-(C2**2-C1**2)/(C2+(2.*PYR(0)-1.)*C1))/C1
        CTHE(3-JT)=MIN(1.,MAX(-1.,CTHE(3-JT)))
        PHIR=PARU(2)*PYR(0)
        CPHI=COS(PHIR)
        ANG=CTHE(1)*CTHE(2)-SQRT(1.-CTHE(1)**2)*SQRT(1.-CTHE(2)**2)*CPHI
        Z1=2.-Z(JT)
        Z2=ANG*SQRT(Z(JT)**2-4.*PMQ(JT)**2/SHP)
        Z3=1.-Z(JT)-XH+(PMQ(1)**2+PMQ(2)**2)/SHP
        Z(3-JT)=2./(Z1**2-Z2**2)*(Z1*Z3+Z2*SQRT(Z3**2-(Z1**2-Z2**2)*
     &  PMQ(3-JT)**2/SHP))
        ZMIN=2.*PMQ(3-JT)/SHPR
        ZMAX=1.-PMQ(JT)/SHPR-(SH-PMQ(3-JT)**2)/(SHPR*(SHPR-PMQ(JT)))
        ZMAX=MIN(1.-XH,ZMAX)
        IF(Z(3-JT).LT.ZMIN.OR.Z(3-JT).GT.ZMAX) GOTO 270
        KCC=22

      ELSEIF(ISUB.EQ.73) THEN
C...Z0 + W+/- -> Z0 + W+/-.
        XH=SH/SHP
  280   JT=INT(1.5+PYR(0))
        I=MINT(14+JT)
        IA=IABS(I)
        IF(IA.LE.10) THEN
          RVCKM=VINT(180+I)*PYR(0)
          DO 290 J=1,MSTP(1)
          IB=2*J-1+MOD(IA,2)
          IPM=(5-ISIGN(1,I))/2
          IDC=J+MDCY(IA,2)+2
          IF(MDME(IDC,1).NE.1.AND.MDME(IDC,1).NE.IPM) GOTO 290
          MINT(20+JT)=ISIGN(IB,I)
          RVCKM=RVCKM-VCKM((IA+1)/2,(IB+1)/2)
          IF(RVCKM.LE.0.) GOTO 300
  290     CONTINUE
        ELSE
          IB=2*((IA+1)/2)-1+MOD(IA,2)
          MINT(20+JT)=ISIGN(IB,I)
        ENDIF
  300   PMQ(JT)=ULMASS(MINT(20+JT))
        MINT(23-JT)=MINT(17-JT)
        PMQ(3-JT)=ULMASS(MINT(23-JT))
        JT=INT(1.5+PYR(0))
        ZMIN=2.*PMQ(JT)/SHPR
        ZMAX=1.-PMQ(3-JT)/SHPR-(SH-PMQ(JT)**2)/(SHPR*(SHPR-PMQ(3-JT)))
        ZMAX=MIN(1.-XH,ZMAX)
        IF(ZMIN.GE.ZMAX) GOTO 280
        Z(JT)=ZMIN+(ZMAX-ZMIN)*PYR(0)
        IF(-1.+(1.+XH)/(1.-Z(JT))-XH/(1.-Z(JT))**2.LT.
     &  (1.-XH)**2/(4.*XH)*PYR(0)) GOTO 280
        SQC1=1.-4.*PMQ(JT)**2/(Z(JT)**2*SHP)
        IF(SQC1.LT.1.E-8) GOTO 280
        C1=SQRT(SQC1)
        C2=1.+2.*(PMAS(23,1)**2-PMQ(JT)**2)/(Z(JT)*SHP)
        CTHE(JT)=(C2-(C2**2-C1**2)/(C2+(2.*PYR(0)-1.)*C1))/C1
        CTHE(JT)=MIN(1.,MAX(-1.,CTHE(JT)))
        Z(3-JT)=1.-XH/(1.-Z(JT))
        SQC1=1.-4.*PMQ(3-JT)**2/(Z(3-JT)**2*SHP)
        IF(SQC1.LT.1.E-8) GOTO 280
        C1=SQRT(SQC1)
        C2=1.+2.*(PMAS(23,1)**2-PMQ(3-JT)**2)/(Z(3-JT)*SHP)
        CTHE(3-JT)=(C2-(C2**2-C1**2)/(C2+(2.*PYR(0)-1.)*C1))/C1
        CTHE(3-JT)=MIN(1.,MAX(-1.,CTHE(3-JT)))
        PHIR=PARU(2)*PYR(0)
        CPHI=COS(PHIR)
        ANG=CTHE(1)*CTHE(2)-SQRT(1.-CTHE(1)**2)*SQRT(1.-CTHE(2)**2)*CPHI
        Z1=2.-Z(JT)
        Z2=ANG*SQRT(Z(JT)**2-4.*PMQ(JT)**2/SHP)
        Z3=1.-Z(JT)-XH+(PMQ(1)**2+PMQ(2)**2)/SHP
        Z(3-JT)=2./(Z1**2-Z2**2)*(Z1*Z3+Z2*SQRT(Z3**2-(Z1**2-Z2**2)*
     &  PMQ(3-JT)**2/SHP))
        ZMIN=2.*PMQ(3-JT)/SHPR
        ZMAX=1.-PMQ(JT)/SHPR-(SH-PMQ(3-JT)**2)/(SHPR*(SHPR-PMQ(JT)))
        ZMAX=MIN(1.-XH,ZMAX)
        IF(Z(3-JT).LT.ZMIN.OR.Z(3-JT).GT.ZMAX) GOTO 280
        KCC=22

      ELSEIF(ISUB.EQ.74) THEN
C...Z0 + H0 -> Z0 + H0.

      ELSEIF(ISUB.EQ.75) THEN
C...W+ + W- -> gamma + gamma.

      ELSEIF(ISUB.EQ.76.OR.ISUB.EQ.77) THEN
C...W+ + W- -> Z0 + Z0; W+ + W- -> W+ + W-.
        XH=SH/SHP
  310   DO 330 JT=1,2
        I=MINT(14+JT)
        IA=IABS(I)
        IF(IA.LE.10) THEN
          RVCKM=VINT(180+I)*PYR(0)
          DO 320 J=1,MSTP(1)
          IB=2*J-1+MOD(IA,2)
          IPM=(5-ISIGN(1,I))/2
          IDC=J+MDCY(IA,2)+2
          IF(MDME(IDC,1).NE.1.AND.MDME(IDC,1).NE.IPM) GOTO 320
          MINT(20+JT)=ISIGN(IB,I)
          RVCKM=RVCKM-VCKM((IA+1)/2,(IB+1)/2)
          IF(RVCKM.LE.0.) GOTO 330
  320     CONTINUE
        ELSE
          IB=2*((IA+1)/2)-1+MOD(IA,2)
          MINT(20+JT)=ISIGN(IB,I)
        ENDIF
  330   PMQ(JT)=ULMASS(MINT(20+JT))
        JT=INT(1.5+PYR(0))
        ZMIN=2.*PMQ(JT)/SHPR
        ZMAX=1.-PMQ(3-JT)/SHPR-(SH-PMQ(JT)**2)/(SHPR*(SHPR-PMQ(3-JT)))
        ZMAX=MIN(1.-XH,ZMAX)
        IF(ZMIN.GE.ZMAX) GOTO 310
        Z(JT)=ZMIN+(ZMAX-ZMIN)*PYR(0)
        IF(-1.+(1.+XH)/(1.-Z(JT))-XH/(1.-Z(JT))**2.LT.
     &  (1.-XH)**2/(4.*XH)*PYR(0)) GOTO 310
        SQC1=1.-4.*PMQ(JT)**2/(Z(JT)**2*SHP)
        IF(SQC1.LT.1.E-8) GOTO 310
        C1=SQRT(SQC1)
        C2=1.+2.*(PMAS(24,1)**2-PMQ(JT)**2)/(Z(JT)*SHP)
        CTHE(JT)=(C2-(C2**2-C1**2)/(C2+(2.*PYR(0)-1.)*C1))/C1
        CTHE(JT)=MIN(1.,MAX(-1.,CTHE(JT)))
        Z(3-JT)=1.-XH/(1.-Z(JT))
        SQC1=1.-4.*PMQ(3-JT)**2/(Z(3-JT)**2*SHP)
        IF(SQC1.LT.1.E-8) GOTO 310
        C1=SQRT(SQC1)
        C2=1.+2.*(PMAS(24,1)**2-PMQ(3-JT)**2)/(Z(3-JT)*SHP)
        CTHE(3-JT)=(C2-(C2**2-C1**2)/(C2+(2.*PYR(0)-1.)*C1))/C1
        CTHE(3-JT)=MIN(1.,MAX(-1.,CTHE(3-JT)))
        PHIR=PARU(2)*PYR(0)
        CPHI=COS(PHIR)
        ANG=CTHE(1)*CTHE(2)-SQRT(1.-CTHE(1)**2)*SQRT(1.-CTHE(2)**2)*CPHI
        Z1=2.-Z(JT)
        Z2=ANG*SQRT(Z(JT)**2-4.*PMQ(JT)**2/SHP)
        Z3=1.-Z(JT)-XH+(PMQ(1)**2+PMQ(2)**2)/SHP
        Z(3-JT)=2./(Z1**2-Z2**2)*(Z1*Z3+Z2*SQRT(Z3**2-(Z1**2-Z2**2)*
     &  PMQ(3-JT)**2/SHP))
        ZMIN=2.*PMQ(3-JT)/SHPR
        ZMAX=1.-PMQ(JT)/SHPR-(SH-PMQ(3-JT)**2)/(SHPR*(SHPR-PMQ(JT)))
        ZMAX=MIN(1.-XH,ZMAX)
        IF(Z(3-JT).LT.ZMIN.OR.Z(3-JT).GT.ZMAX) GOTO 310
        KCC=22

      ELSEIF(ISUB.EQ.78) THEN
C...W+/- + H0 -> W+/- + H0.

      ELSEIF(ISUB.EQ.79) THEN
C...H0 + H0 -> H0 + H0.
      ENDIF

      ELSEIF(ISUB.LE.90) THEN
      IF(ISUB.EQ.81) THEN
C...q + q~ -> Q + Q~; th = (p(q)-p(Q))**2.
        MINT(21)=ISIGN(MINT(55),MINT(15))
        MINT(22)=-MINT(21)
        KCC=4

      ELSEIF(ISUB.EQ.82) THEN
C...g + g -> Q + Q~; th arbitrary.
        KCS=(-1)**INT(1.5+PYR(0))
        MINT(21)=ISIGN(MINT(55),KCS)
        MINT(22)=-MINT(21)
        KCC=MINT(2)+10

      ELSEIF(ISUB.EQ.83) THEN
C...f + q -> f' + Q; th = (p(f) - p(f'))**2.
        KFOLD=MINT(16)
        IF(MINT(2).EQ.2) KFOLD=MINT(15)
        KFAOLD=IABS(KFOLD)
        IF(KFAOLD.GT.10) THEN
          KFANEW=KFAOLD+2*MOD(KFAOLD,2)-1
        ELSE
          RCKM=VINT(180+KFOLD)*PYR(0)
          IPM=(5-ISIGN(1,KFOLD))/2
          KFANEW=-MOD(KFAOLD+1,2)
  340     KFANEW=KFANEW+2
          IDC=MDCY(KFAOLD,2)+(KFANEW+1)/2+2
          IF(MDME(IDC,1).EQ.1.OR.MDME(IDC,1).EQ.IPM) THEN
            IF(MOD(KFAOLD,2).EQ.0) RCKM=RCKM-VCKM(KFAOLD/2,(KFANEW+1)/2)
            IF(MOD(KFAOLD,2).EQ.1) RCKM=RCKM-VCKM(KFANEW/2,(KFAOLD+1)/2)
          ENDIF
          IF(KFANEW.LE.6.AND.RCKM.GT.0.) GOTO 340
        ENDIF
        IF(MINT(2).EQ.1) THEN
          MINT(21)=ISIGN(MINT(55),MINT(15))
          MINT(22)=ISIGN(KFANEW,MINT(16))
        ELSE
          MINT(21)=ISIGN(KFANEW,MINT(15))
          MINT(22)=ISIGN(MINT(55),MINT(16))
        ENDIF
        KCC=22

      ELSEIF(ISUB.EQ.84) THEN
C...g + gamma -> Q + Q~; th arbitary.
        KCS=(-1)**INT(1.5+PYR(0))
        MINT(21)=ISIGN(MINT(55),KCS)
        MINT(22)=-MINT(21)
        KCC=27
        IF(MINT(16).EQ.21) KCC=28

      ELSEIF(ISUB.EQ.85) THEN
C...gamma + gamma -> F + F~; th arbitary.
        KCS=(-1)**INT(1.5+PYR(0))
        MINT(21)=ISIGN(MINT(56),KCS)
        MINT(22)=-MINT(21)
        KCC=21
      ENDIF

      ELSEIF(ISUB.LE.100) THEN
      IF(ISUB.EQ.95) THEN
C...Low-pT ( = energyless g + g -> g + g).
        KCC=MINT(2)+12
        KCS=(-1)**INT(1.5+PYR(0))

      ELSEIF(ISUB.EQ.96) THEN
C...Multiple interactions (should be reassigned to QCD process).
      ENDIF

      ELSEIF(ISUB.LE.110) THEN
      IF(ISUB.EQ.101) THEN
C...g + g -> gamma*/Z0.
        KCC=21
        KFRES=22

      ELSEIF(ISUB.EQ.102) THEN
C...g + g -> H0 (or H'0, or A0).
        KCC=21
        KFRES=KFHIGG

      ELSEIF(ISUB.EQ.103) THEN
C...gamma + gamma -> H0 (or H'0, or A0).
        KCC=21
        KFRES=KFHIGG
      ENDIF

      ELSEIF(ISUB.LE.120) THEN
      IF(ISUB.EQ.111) THEN
C...f + f~ -> g + H0; th arbitrary.
        IF(PYR(0).GT.0.5) JS=2
        MINT(20+JS)=21
        MINT(23-JS)=25
        KCC=17+JS

      ELSEIF(ISUB.EQ.112) THEN
C...f + g -> f + H0; th = (p(f) - p(f))**2.
        IF(MINT(15).EQ.21) JS=2
        MINT(23-JS)=25
        KCC=15+JS
        KCS=ISIGN(1,MINT(14+JS))

      ELSEIF(ISUB.EQ.113) THEN
C...g + g -> g + H0; th arbitrary.
        IF(PYR(0).GT.0.5) JS=2
        MINT(23-JS)=25
        KCC=22+JS
        KCS=(-1)**INT(1.5+PYR(0))

      ELSEIF(ISUB.EQ.114) THEN
C...g + g -> gamma + gamma; th arbitrary.
        IF(PYR(0).GT.0.5) JS=2
        MINT(21)=22
        MINT(22)=22
        KCC=21

      ELSEIF(ISUB.EQ.115) THEN
C...g + g -> g + gamma; th arbitrary.
        IF(PYR(0).GT.0.5) JS=2
        MINT(23-JS)=22
        KCC=22+JS
        KCS=(-1)**INT(1.5+PYR(0))

      ELSEIF(ISUB.EQ.116) THEN
C...g + g -> gamma + Z0.

      ELSEIF(ISUB.EQ.117) THEN
C...g + g -> Z0 + Z0.

      ELSEIF(ISUB.EQ.118) THEN
C...g + g -> W+ + W-.
      ENDIF

      ELSEIF(ISUB.LE.140) THEN
      IF(ISUB.EQ.121) THEN
C...g + g -> f + f~ + H0 (f + f~ -> H0 as inner process).

      ELSEIF(ISUB.EQ.122) THEN
C...gamma + gamma -> f + f' + H0 (f + f~ -> H0 as inner process).

      ELSEIF(ISUB.EQ.123) THEN
C...f + f' -> f + f' + H0 (or H'0, or A0) (Z0 + Z0 -> H0 as
C...inner process).
        KCC=22
        KFRES=KFHIGG

      ELSEIF(ISUB.EQ.124) THEN
C...f + f' -> f" + f"' + H0 (or H'0, or A) (W+ + W- -> H0 as
C...inner process).
        DO 360 JT=1,2
        I=MINT(14+JT)
        IA=IABS(I)
        IF(IA.LE.10) THEN
          RVCKM=VINT(180+I)*PYR(0)
          DO 350 J=1,MSTP(1)
          IB=2*J-1+MOD(IA,2)
          IPM=(5-ISIGN(1,I))/2
          IDC=J+MDCY(IA,2)+2
          IF(MDME(IDC,1).NE.1.AND.MDME(IDC,1).NE.IPM) GOTO 350
          MINT(20+JT)=ISIGN(IB,I)
          RVCKM=RVCKM-VCKM((IA+1)/2,(IB+1)/2)
          IF(RVCKM.LE.0.) GOTO 360
  350     CONTINUE
        ELSE
          IB=2*((IA+1)/2)-1+MOD(IA,2)
          MINT(20+JT)=ISIGN(IB,I)
        ENDIF
  360   CONTINUE
        KCC=22
        KFRES=KFHIGG

      ELSEIF(ISUB.EQ.131) THEN
C...g + g -> Z0 + q + q~.
        MINT(21)=KFPR(131,1)
        MINT(22)=KFPR(131,2)
        MINT(23)=-MINT(22)
        KCC=MINT(2)+10
        KCS=1
      ENDIF

      ELSEIF(ISUB.LE.160) THEN
      IF(ISUB.EQ.141) THEN
C...f + f~ -> gamma*/Z0/Z'0.
        KFRES=32

      ELSEIF(ISUB.EQ.142) THEN
C...f + f~' -> W'+/- .
        KCH1=KCHG(IABS(MINT(15)),1)*ISIGN(1,MINT(15))
        KCH2=KCHG(IABS(MINT(16)),1)*ISIGN(1,MINT(16))
        KFRES=ISIGN(34,KCH1+KCH2)

      ELSEIF(ISUB.EQ.143) THEN
C...f + f~' -> H+/-.
        KCH1=KCHG(IABS(MINT(15)),1)*ISIGN(1,MINT(15))
        KCH2=KCHG(IABS(MINT(16)),1)*ISIGN(1,MINT(16))
        KFRES=ISIGN(37,KCH1+KCH2)

      ELSEIF(ISUB.EQ.144) THEN
C...f + f~' -> R.
        KFRES=ISIGN(40,MINT(15)+MINT(16))

      ELSEIF(ISUB.EQ.145) THEN
C...q + l -> LQ (leptoquark).
        KFRES=ISIGN(39,MINT(15))
        IF(IABS(MINT(16)).LE.8) KFRES=ISIGN(39,MINT(16))
      ENDIF

      ELSE
      IF(ISUB.EQ.161) THEN
C...f + g -> f' + H+/-; th = (p(f)-p(f'))**2.
        IF(MINT(15).EQ.21) JS=2
        I=MINT(14+JS)
        IA=IABS(I)
        MINT(23-JS)=ISIGN(37,KCHG(IA,1)*I)
        IB=IA+MOD(IA,2)-MOD(IA+1,2)
        MINT(20+JS)=ISIGN(IB,I)
        KCC=15+JS
        KCS=ISIGN(1,MINT(14+JS))

      ELSEIF(ISUB.EQ.162) THEN
C...q + g -> LQ + l~; LQ=leptoquark; th=(p(q)-p(LQ))^2.
        IF(MINT(15).EQ.21) JS=2
        MINT(20+JS)=ISIGN(39,MINT(14+JS))
        KFLQL=KFDP(MDCY(39,2),2)
        MINT(23-JS)=-ISIGN(KFLQL,MINT(14+JS))
        KCC=15+JS
        KCS=ISIGN(1,MINT(14+JS))

      ELSEIF(ISUB.EQ.163) THEN
C...g + g -> LQ + LQ~; LQ=leptoquark; th arbitrary.
        KCS=(-1)**INT(1.5+PYR(0))
        MINT(21)=ISIGN(39,KCS)
        MINT(22)=-MINT(21)
        KCC=MINT(2)+10

      ELSEIF(ISUB.EQ.164) THEN
C...q + q~ -> LQ + LQ~; LQ=leptoquark; th=(p(q)-p(LQ))**2.
        MINT(21)=ISIGN(39,MINT(15))
        MINT(22)=-MINT(21)
        KCC=4
      ENDIF
      ENDIF

      IF(IDOC.EQ.7) THEN
C...Resonance not decaying: store colour connection indices.
        I=MINT(83)+7
        K(IPU3,1)=1
        K(IPU3,2)=KFRES
        K(IPU3,3)=I
        P(IPU3,4)=SHUSER
        P(IPU3,5)=SHUSER
        K(IPU1,4)=IPU2
        K(IPU1,5)=IPU2
        K(IPU2,4)=IPU1
        K(IPU2,5)=IPU1
        K(I,1)=21
        K(I,2)=KFRES
        P(I,4)=SHUSER
        P(I,5)=SHUSER
        N=IPU3
        MINT(21)=KFRES
        MINT(22)=0

C...Special case: colour flow in q + g -> LQ.
        IF(IABS(KFRES).EQ.39) THEN
          IF(MINT(15).GT.0.AND.MINT(15).LE.8) THEN
            K(IPU1,4)=IPU3
            K(IPU3,4)=MSTU(5)*IPU1
          ELSEIF(MINT(15).LT.0.AND.MINT(15).GE.-8) THEN
            K(IPU1,5)=IPU3
            K(IPU3,5)=MSTU(5)*IPU1
          ELSEIF(MINT(16).GT.0) THEN
            K(IPU2,4)=IPU3
            K(IPU3,4)=MSTU(5)*IPU2
          ELSE
            K(IPU2,5)=IPU3
            K(IPU3,5)=MSTU(5)*IPU2
          ENDIF
        ENDIF

      ELSEIF(IDOC.EQ.8) THEN
C...2 -> 2 processes: store outgoing partons in their CM-frame.
        DO 370 JT=1,2
        I=MINT(84)+2+JT
        K(I,1)=1
        IF(IABS(MINT(20+JT)).LE.100) THEN
          IF(KCHG(IABS(MINT(20+JT)),2).NE.0) K(I,1)=3
        ENDIF
        K(I,2)=MINT(20+JT)
        K(I,3)=MINT(83)+IDOC+JT-2
        IF(IABS(K(I,2)).LE.22) THEN
          P(I,5)=ULMASS(K(I,2))
        ELSE
          P(I,5)=SQRT(VINT(63+MOD(JS+JT,2)))
        ENDIF
  370   CONTINUE
        IF(P(IPU3,5)+P(IPU4,5).GE.SHR) THEN
          KFA1=IABS(MINT(21))
          KFA2=IABS(MINT(22))
          IF((KFA1.GT.3.AND.KFA1.NE.21).OR.(KFA2.GT.3.AND.KFA2.NE.21))
     &    THEN
            MINT(51)=1
            RETURN
          ENDIF
          P(IPU3,5)=0.
          P(IPU4,5)=0.
        ENDIF
        P(IPU3,4)=0.5*(SHR+(P(IPU3,5)**2-P(IPU4,5)**2)/SHR)
        P(IPU3,3)=SQRT(MAX(0.,P(IPU3,4)**2-P(IPU3,5)**2))
        P(IPU4,4)=SHR-P(IPU3,4)
        P(IPU4,3)=-P(IPU3,3)
        N=IPU4
        MINT(7)=MINT(83)+7
        MINT(8)=MINT(83)+8

C...Rotate outgoing partons using cos(theta)=(th-uh)/lam(sh,sqm3,sqm4).
        CALL LUDBRB(IPU3,IPU4,ACOS(VINT(23)),VINT(24),0D0,0D0,0D0)

      ELSEIF(IDOC.EQ.9.AND.ISET(ISUB).EQ.5) THEN
C...2 -> 3 processes (alt 1): store outgoing partons in their CM frame.
        DO 380 JT=1,2
        I=MINT(84)+2+JT
        K(I,1)=1
        IF(IABS(MINT(20+JT)).LE.100) THEN
          IF(KCHG(IABS(MINT(20+JT)),2).NE.0) K(I,1)=3
        ENDIF
        K(I,2)=MINT(20+JT)
        K(I,3)=MINT(83)+IDOC+JT-3
        IF(IABS(K(I,2)).LE.22) THEN
          P(I,5)=ULMASS(K(I,2))
        ELSE
          P(I,5)=SQRT(VINT(63+MOD(JS+JT,2)))
        ENDIF
        PT=SQRT(MAX(0.,VINT(197+5*JT)-P(I,5)**2+VINT(196+5*JT)**2))
        P(I,1)=PT*COS(VINT(198+5*JT))
        P(I,2)=PT*SIN(VINT(198+5*JT))
  380   CONTINUE
        K(IPU5,1)=1
        K(IPU5,2)=KFRES
        K(IPU5,3)=MINT(83)+IDOC
        P(IPU5,5)=SHR
        P(IPU5,1)=-P(IPU3,1)-P(IPU4,1)
        P(IPU5,2)=-P(IPU3,2)-P(IPU4,2)
        PMS1=P(IPU3,5)**2+P(IPU3,1)**2+P(IPU3,2)**2
        PMS2=P(IPU4,5)**2+P(IPU4,1)**2+P(IPU4,2)**2
        PMS3=P(IPU5,5)**2+P(IPU5,1)**2+P(IPU5,2)**2
        PMT3=SQRT(PMS3)
        P(IPU5,3)=PMT3*SINH(VINT(211))
        P(IPU5,4)=PMT3*COSH(VINT(211))
        PMS12=(SHPR-P(IPU5,4))**2-P(IPU5,3)**2
        SQL12=(PMS12-PMS1-PMS2)**2-4.*PMS1*PMS2
        IF(SQL12.LE.0.) THEN
          MINT(51)=1
          RETURN
        ENDIF
        P(IPU3,3)=(-P(IPU5,3)*(PMS12+PMS1-PMS2)+
     &  VINT(213)*(SHPR-P(IPU5,4))*SQRT(SQL12))/(2.*PMS12)
        P(IPU4,3)=-P(IPU3,3)-P(IPU5,3)
        P(IPU3,4)=SQRT(PMS1+P(IPU3,3)**2)
        P(IPU4,4)=SQRT(PMS2+P(IPU4,3)**2)
        N=IPU5
        MINT(7)=MINT(83)+7
        MINT(8)=MINT(83)+8

      ELSEIF(IDOC.EQ.9) THEN
C...2 -> 3 processes: store outgoing partons in their CM frame.
        DO 390 JT=1,3
        I=MINT(84)+2+JT
        K(I,1)=1
        IF(IABS(MINT(20+JT)).LE.10.OR.MINT(20+JT).EQ.21) K(I,1)=3
        K(I,2)=MINT(20+JT)
        K(I,3)=MINT(83)+IDOC+JT-3
        IF(JT.EQ.1) THEN
          P(I,5)=SQRT(VINT(63))
        ELSE
          P(I,5)=PMAS(KFPR(ISUB,2),1)
        ENDIF
  390   CONTINUE
        P(IPU3,4)=0.5*(SHR+(VINT(63)-VINT(64))/SHR)
        P(IPU3,3)=SQRT(MAX(0.,P(IPU3,4)**2-P(IPU3,5)**2))
        P(IPU4,4)=0.5*SQRT(VINT(64))
        P(IPU4,3)=SQRT(MAX(0.,P(IPU4,4)**2-P(IPU4,5)**2))
        P(IPU5,4)=P(IPU4,4)
        P(IPU5,3)=-P(IPU4,3)
        N=IPU5
        MINT(7)=MINT(83)+7
        MINT(8)=MINT(83)+9

C...Rotate and boost outgoing partons.
        CALL LUDBRB(IPU4,IPU5,ACOS(VINT(83)),VINT(84),0D0,0D0,0D0)
        CALL LUDBRB(IPU4,IPU5,0.,0.,0D0,0D0,
     &  -DBLE(P(IPU3,3)/(SHR-P(IPU3,4))))
        CALL LUDBRB(IPU3,IPU5,ACOS(VINT(23)),VINT(24),0D0,0D0,0D0)

      ELSEIF(IDOC.EQ.11) THEN
C...Z0 + Z0 -> H0, W+ + W- -> H0: store Higgs and outgoing partons.
        PHI(1)=PARU(2)*PYR(0)
        PHI(2)=PHI(1)-PHIR
        DO 400 JT=1,2
        I=MINT(84)+2+JT
        K(I,1)=1
        IF(IABS(MINT(20+JT)).LE.10.OR.MINT(20+JT).EQ.21) K(I,1)=3
        K(I,2)=MINT(20+JT)
        K(I,3)=MINT(83)+IDOC+JT-2
        P(I,5)=ULMASS(K(I,2))
        IF(0.5*SHPR*Z(JT).LE.P(I,5)) P(I,5)=0.
        PABS=SQRT(MAX(0.,(0.5*SHPR*Z(JT))**2-P(I,5)**2))
        PTABS=PABS*SQRT(MAX(0.,1.-CTHE(JT)**2))
        P(I,1)=PTABS*COS(PHI(JT))
        P(I,2)=PTABS*SIN(PHI(JT))
        P(I,3)=PABS*CTHE(JT)*(-1)**(JT+1)
        P(I,4)=0.5*SHPR*Z(JT)
        IZW=MINT(83)+6+JT
        K(IZW,1)=21
        K(IZW,2)=23
        IF(ISUB.EQ.8) K(IZW,2)=ISIGN(24,LUCHGE(MINT(14+JT)))
        K(IZW,3)=IZW-2
        P(IZW,1)=-P(I,1)
        P(IZW,2)=-P(I,2)
        P(IZW,3)=(0.5*SHPR-PABS*CTHE(JT))*(-1)**(JT+1)
        P(IZW,4)=0.5*SHPR*(1.-Z(JT))
  400   P(IZW,5)=-SQRT(MAX(0.,P(IZW,3)**2+PTABS**2-P(IZW,4)**2))
        I=MINT(83)+9
        K(IPU5,1)=1
        K(IPU5,2)=KFRES
        K(IPU5,3)=I
        P(IPU5,5)=SHR
        P(IPU5,1)=-P(IPU3,1)-P(IPU4,1)
        P(IPU5,2)=-P(IPU3,2)-P(IPU4,2)
        P(IPU5,3)=-P(IPU3,3)-P(IPU4,3)
        P(IPU5,4)=SHPR-P(IPU3,4)-P(IPU4,4)
        K(I,1)=21
        K(I,2)=KFRES
        DO 410 J=1,5
  410   P(I,J)=P(IPU5,J)
        N=IPU5
        MINT(23)=KFRES

      ELSEIF(IDOC.EQ.12) THEN
C...Z0 and W+/- scattering: store bosons and outgoing partons.
        PHI(1)=PARU(2)*PYR(0)
        PHI(2)=PHI(1)-PHIR
        JTRAN=INT(1.5+PYR(0))
        DO 420 JT=1,2
        I=MINT(84)+2+JT
        K(I,1)=1
        IF(IABS(MINT(20+JT)).LE.10.OR.MINT(20+JT).EQ.21) K(I,1)=3
        K(I,2)=MINT(20+JT)
        K(I,3)=MINT(83)+IDOC+JT-2
        P(I,5)=ULMASS(K(I,2))
        IF(0.5*SHPR*Z(JT).LE.P(I,5)) P(I,5)=0.
        PABS=SQRT(MAX(0.,(0.5*SHPR*Z(JT))**2-P(I,5)**2))
        PTABS=PABS*SQRT(MAX(0.,1.-CTHE(JT)**2))
        P(I,1)=PTABS*COS(PHI(JT))
        P(I,2)=PTABS*SIN(PHI(JT))
        P(I,3)=PABS*CTHE(JT)*(-1)**(JT+1)
        P(I,4)=0.5*SHPR*Z(JT)
        IZW=MINT(83)+6+JT
        K(IZW,1)=21
        IF(MINT(14+JT).EQ.MINT(20+JT)) THEN
          K(IZW,2)=23
        ELSE
          K(IZW,2)=ISIGN(24,LUCHGE(MINT(14+JT))-LUCHGE(MINT(20+JT)))
        ENDIF
        K(IZW,3)=IZW-2
        P(IZW,1)=-P(I,1)
        P(IZW,2)=-P(I,2)
        P(IZW,3)=(0.5*SHPR-PABS*CTHE(JT))*(-1)**(JT+1)
        P(IZW,4)=0.5*SHPR*(1.-Z(JT))
        P(IZW,5)=-SQRT(MAX(0.,P(IZW,3)**2+PTABS**2-P(IZW,4)**2))
        IPU=MINT(84)+4+JT
        K(IPU,1)=3
        K(IPU,2)=KFPR(ISUB,JT)
        IF(ISUB.EQ.72.AND.JT.EQ.JTRAN) K(IPU,2)=-K(IPU,2)
        IF(ISUB.EQ.73.OR.ISUB.EQ.77) K(IPU,2)=K(IZW,2)
        K(IPU,3)=MINT(83)+8+JT
        IF(IABS(K(IPU,2)).LE.10.OR.K(IPU,2).EQ.21) THEN
          P(IPU,5)=ULMASS(K(IPU,2))
        ELSE
          P(IPU,5)=SQRT(VINT(63+MOD(JS+JT,2)))
        ENDIF
        MINT(22+JT)=K(IPU,2)
  420   CONTINUE
C...Find rotation and boost for hard scattering subsystem.
        I1=MINT(83)+7
        I2=MINT(83)+8
        BEXCM=(P(I1,1)+P(I2,1))/(P(I1,4)+P(I2,4))
        BEYCM=(P(I1,2)+P(I2,2))/(P(I1,4)+P(I2,4))
        BEZCM=(P(I1,3)+P(I2,3))/(P(I1,4)+P(I2,4))
        GAMCM=(P(I1,4)+P(I2,4))/SHR
        BEPCM=BEXCM*P(I1,1)+BEYCM*P(I1,2)+BEZCM*P(I1,3)
        PX=P(I1,1)+GAMCM*(GAMCM/(1.+GAMCM)*BEPCM-P(I1,4))*BEXCM
        PY=P(I1,2)+GAMCM*(GAMCM/(1.+GAMCM)*BEPCM-P(I1,4))*BEYCM
        PZ=P(I1,3)+GAMCM*(GAMCM/(1.+GAMCM)*BEPCM-P(I1,4))*BEZCM
        THECM=ULANGL(PZ,SQRT(PX**2+PY**2))
        PHICM=ULANGL(PX,PY)
C...Store hard scattering subsystem. Rotate and boost it.
        SQLAM=(SH-P(IPU5,5)**2-P(IPU6,5)**2)**2-4.*P(IPU5,5)**2*
     &  P(IPU6,5)**2
        PABS=SQRT(MAX(0.,SQLAM/(4.*SH)))
        CTHWZ=VINT(23)
        STHWZ=SQRT(MAX(0.,1.-CTHWZ**2))
        PHIWZ=VINT(24)-PHICM
        P(IPU5,1)=PABS*STHWZ*COS(PHIWZ)
        P(IPU5,2)=PABS*STHWZ*SIN(PHIWZ)
        P(IPU5,3)=PABS*CTHWZ
        P(IPU5,4)=SQRT(PABS**2+P(IPU5,5)**2)
        P(IPU6,1)=-P(IPU5,1)
        P(IPU6,2)=-P(IPU5,2)
        P(IPU6,3)=-P(IPU5,3)
        P(IPU6,4)=SQRT(PABS**2+P(IPU6,5)**2)
        CALL LUDBRB(IPU5,IPU6,THECM,PHICM,DBLE(BEXCM),DBLE(BEYCM),
     &  DBLE(BEZCM))
        DO 430 JT=1,2
        I1=MINT(83)+8+JT
        I2=MINT(84)+4+JT
        K(I1,1)=21
        K(I1,2)=K(I2,2)
        DO 430 J=1,5
  430   P(I1,J)=P(I2,J)
        N=IPU6
        MINT(7)=MINT(83)+9
        MINT(8)=MINT(83)+10
      ENDIF

      IF(IDOC.GE.8.AND.ISET(ISUB).NE.6) THEN
C...Store colour connection indices.
        DO 440 J=1,2
        JC=J
        IF(KCS.EQ.-1) JC=3-J
        IF(ICOL(KCC,1,JC).NE.0.AND.K(IPU1,1).EQ.14) K(IPU1,J+3)=
     &  K(IPU1,J+3)+MINT(84)+ICOL(KCC,1,JC)
        IF(ICOL(KCC,2,JC).NE.0.AND.K(IPU2,1).EQ.14) K(IPU2,J+3)=
     &  K(IPU2,J+3)+MINT(84)+ICOL(KCC,2,JC)
        IF(ICOL(KCC,3,JC).NE.0.AND.K(IPU3,1).EQ.3) K(IPU3,J+3)=
     &  MSTU(5)*(MINT(84)+ICOL(KCC,3,JC))
  440   IF(ICOL(KCC,4,JC).NE.0.AND.K(IPU4,1).EQ.3) K(IPU4,J+3)=
     &  MSTU(5)*(MINT(84)+ICOL(KCC,4,JC))

C...Copy outgoing partons to documentation lines.
        IMAX=2
        IF(IDOC.EQ.9) IMAX=3
        DO 450 I=1,IMAX
        I1=MINT(83)+IDOC-IMAX+I
        I2=MINT(84)+2+I
        K(I1,1)=21
        K(I1,2)=K(I2,2)
        IF(IDOC.LE.9) K(I1,3)=0
        IF(IDOC.GE.11) K(I1,3)=MINT(83)+2+I
        DO 450 J=1,5
  450   P(I1,J)=P(I2,J)

      ELSEIF(IDOC.EQ.9) THEN
C...Store colour connection indices.
        DO 460 J=1,2
        JC=J
        IF(KCS.EQ.-1) JC=3-J
        IF(ICOL(KCC,1,JC).NE.0.AND.K(IPU1,1).EQ.14) K(IPU1,J+3)=
     &  K(IPU1,J+3)+MINT(84)+ICOL(KCC,1,JC)+
     &  MAX(0,MIN(1,ICOL(KCC,1,JC)-2))
        IF(ICOL(KCC,2,JC).NE.0.AND.K(IPU2,1).EQ.14) K(IPU2,J+3)=
     &  K(IPU2,J+3)+MINT(84)+ICOL(KCC,2,JC)+
     &  MAX(0,MIN(1,ICOL(KCC,2,JC)-2))
        IF(ICOL(KCC,3,JC).NE.0.AND.K(IPU4,1).EQ.3) K(IPU4,J+3)=
     &  MSTU(5)*(MINT(84)+ICOL(KCC,3,JC))
  460   IF(ICOL(KCC,4,JC).NE.0.AND.K(IPU5,1).EQ.3) K(IPU5,J+3)=
     &  MSTU(5)*(MINT(84)+ICOL(KCC,4,JC))

C...Copy outgoing partons to documentation lines.
        DO 470 I=1,3
        I1=MINT(83)+IDOC-3+I
        I2=MINT(84)+2+I
        K(I1,1)=21
        K(I1,2)=K(I2,2)
        K(I1,3)=0
        DO 470 J=1,5
  470   P(I1,J)=P(I2,J)
      ENDIF

C...Low-pT events: remove gluons used for string drawing purposes.
      IF(ISUB.EQ.95) THEN
        K(IPU3,1)=K(IPU3,1)+10
        K(IPU4,1)=K(IPU4,1)+10
        DO 480 J=41,66
  480   VINT(J)=0.
        DO 490 I=MINT(83)+5,MINT(83)+8
        DO 490 J=1,5
  490   P(I,J)=0.
      ENDIF

      RETURN
      END

C*********************************************************************

      SUBROUTINE RYSSPA(IPU1,IPU2)
C...Generates spacelike parton showers.
      IMPLICIT DOUBLE PRECISION(D)
      PARAMETER (KSZJ=4000)

      
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      integer N,NPAD,K
      double precision P,V
      SAVE /PYJETS/


      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/RYSUBS/MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200)
      COMMON/RYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/RYINT1/MINT(400),VINT(400)
      COMMON/RYINT2/ISET(200),KFPR(200,2),COEF(200,20),ICOL(40,4,2)
      COMMON/RYINT3/XSFX(2,-40:40),ISIG(1000,3),SIGH(1000)
      SAVE /LUDAT1/,/LUDAT2/
      SAVE /RYSUBS/,/RYPARS/,/RYINT1/,/RYINT2/,/RYINT3/
      DIMENSION KFLS(4),IS(2),XS(2),ZS(2),Q2S(2),TEVS(2),ROBO(5),
     &XFS(2,-25:25),XFA(-25:25),XFB(-25:25),XFN(-25:25),WTAP(-25:25),
     &WTSF(-25:25),THE2(2),ALAM(2),DQ2(3),DPC(3),DPD(4),DPB(4)

C...Calculate maximum virtuality and check that evolution possible.
      IPUS1=IPU1
      IPUS2=IPU2
      ISUB=MINT(1)
      Q2E=VINT(52)
      IF(ISET(ISUB).EQ.1) THEN
        Q2E=Q2E/PARP(67)
      ELSEIF(ISET(ISUB).GE.3.AND.ISET(ISUB).LE.5) THEN
        Q2E=PMAS(23,1)**2
        IF(ISUB.EQ.8.OR.ISUB.EQ.76.OR.ISUB.EQ.77.OR.ISUB.EQ.124.OR.
     &  ISUB.EQ.174.OR.ISUB.EQ.179) Q2E=PMAS(24,1)**2
      ENDIF
      TMAX=LOG(PARP(67)*PARP(63)*Q2E/PARP(61)**2)
      IF(PARP(67)*Q2E.LT.MAX(PARP(62)**2,2.*PARP(61)**2).OR.
     &TMAX.LT.0.2) RETURN

C...Common constants and initial values. Save normal Lambda value.
      XE0=2.*PARP(65)/VINT(1)
      ALAMS=PARU(111)
      PARU(111)=PARP(61)
      NS=N
  100 N=NS
      DO 110 JT=1,2
      KFLS(JT)=MINT(14+JT)
      KFLS(JT+2)=KFLS(JT)
      XS(JT)=VINT(40+JT)
      ZS(JT)=1.
      Q2S(JT)=PARP(67)*Q2E
      TEVS(JT)=TMAX
      ALAM(JT)=PARP(61)
      THE2(JT)=100.
      DO 110 KFL=-25,25
  110 XFS(JT,KFL)=XSFX(JT,KFL)
      DSH=VINT(44)
      IF(ISET(ISUB).GE.3.AND.ISET(ISUB).LE.5) DSH=VINT(26)*VINT(2)

C...Pick up leg with highest virtuality.
  120 N=N+1
      JT=1
      IF(N.GT.NS+1.AND.Q2S(2).GT.Q2S(1)) JT=2
      KFLB=KFLS(JT)
      XB=XS(JT)
      DO 130 KFL=-25,25
  130 XFB(KFL)=XFS(JT,KFL)
      DSHR=2D0*SQRT(DSH)
      DSHZ=DSH/DBLE(ZS(JT))
      XE=MAX(XE0,XB*(1./(1.-PARP(66))-1.))
      IF(MINT(40+JT).EQ.1.OR.XB+XE.GE.0.999) THEN
        Q2B=0.
        GOTO 220
      ENDIF

C...Maximum Q2 without or with Q2 ordering. Effective Lambda and n_f.
      IF(MSTP(62).LE.1) THEN
        Q2B=0.5*(1./ZS(JT)+1.)*Q2S(JT)+0.5*(1./ZS(JT)-1.)*(Q2S(3-JT)-
     &  SNGL(DSH)+SQRT((SNGL(DSH)+Q2S(1)+Q2S(2))**2+8.*Q2S(1)*Q2S(2)*
     &  ZS(JT)/(1.-ZS(JT))))
        TEVB=LOG(PARP(63)*Q2B/ALAM(JT)**2)
      ELSE
        Q2B=Q2S(JT)
        TEVB=TEVS(JT)
      ENDIF
      ALSDUM=ULALPS(PARP(63)*Q2B)
      TEVB=TEVB+2.*LOG(ALAM(JT)/PARU(117))
      TEVBSV=TEVB
      ALAM(JT)=PARU(117)
      B0=(33.-2.*MSTU(118))/6.

C...Calculate Altarelli-Parisi and structure function weights.
      DO 140 KFL=-25,25
      WTAP(KFL)=0.
  140 WTSF(KFL)=0.
      IF(KFLB.EQ.21) THEN
        WTAPQ=16.*(1.-SQRT(XB+XE))/(3.*SQRT(XB))
        DO 150 KFL=-MSTP(54),MSTP(54)
        IF(KFL.EQ.0) WTAP(KFL)=6.*LOG((1.-XB)/XE)
  150   IF(KFL.NE.0) WTAP(KFL)=WTAPQ
      ELSE
        WTAP(0)=0.5*XB*(1./(XB+XE)-1.)
        WTAP(KFLB)=8.*LOG((1.-XB)*(XB+XE)/XE)/3.
      ENDIF
  160 WTSUM=0.
      IF(KFLB.NE.21) XFBO=XFB(KFLB)
      IF(KFLB.EQ.21) XFBO=XFB(0)
      DO 170 KFL=-MSTP(54),MSTP(54)
      WTSF(KFL)=XFB(KFL)/XFBO
  170 WTSUM=WTSUM+WTAP(KFL)*WTSF(KFL)
      WTSUM=MAX(0.0001,WTSUM)

C...Choose new t: fix alpha_s, alpha_s(Q^2), alpha_s(k_T^2).
  180 IF(MSTP(64).LE.0) THEN
        TEVB=TEVB+LOG(PYR(0))*PARU(2)/(PARU(111)*WTSUM)
      ELSEIF(MSTP(64).EQ.1) THEN
        TEVB=TEVB*EXP(MAX(-70.,LOG(PYR(0))*B0/WTSUM))
      ELSE
        TEVB=TEVB*EXP(MAX(-70.,LOG(PYR(0))*B0/(5.*WTSUM)))
      ENDIF
  190 Q2REF=ALAM(JT)**2*EXP(TEVB)
      Q2B=Q2REF/PARP(63)

C...Evolution ended or select flavour for branching parton.
      IF(Q2B.LT.PARP(62)**2) THEN
        Q2B=0.
      ELSE
        WTRAN=PYR(0)*WTSUM
        KFLA=-MSTP(54)-1
  200   KFLA=KFLA+1
        WTRAN=WTRAN-WTAP(KFLA)*WTSF(KFLA)
        IF(KFLA.LT.MSTP(54).AND.WTRAN.GT.0.) GOTO 200
        IF(KFLA.EQ.0) KFLA=21

C...Choose z value and corrective weight.
        IF(KFLB.EQ.21.AND.KFLA.EQ.21) THEN
          Z=1./(1.+((1.-XB)/XB)*(XE/(1.-XB))**PYR(0))
          WTZ=(1.-Z*(1.-Z))**2
        ELSEIF(KFLB.EQ.21) THEN
          Z=XB/(1.-PYR(0)*(1.-SQRT(XB+XE)))**2
          WTZ=0.5*(1.+(1.-Z)**2)*SQRT(Z)
        ELSEIF(KFLA.EQ.21) THEN
          Z=XB*(1.+PYR(0)*(1./(XB+XE)-1.))
          WTZ=1.-2.*Z*(1.-Z)
        ELSE
          Z=1.-(1.-XB)*(XE/((XB+XE)*(1.-XB)))**PYR(0)
          WTZ=0.5*(1.+Z**2)
        ENDIF

C...Option with resummation of soft gluon emission as effective z shift.
        IF(MSTP(65).GE.1) THEN
          RSOFT=6.
          IF(KFLB.NE.21) RSOFT=8./3.
          Z=Z*(TEVB/TEVS(JT))**(RSOFT*XE/((XB+XE)*B0))
          IF(Z.LE.XB) GOTO 180
        ENDIF

C...Option with alpha_s(k_T^2): demand k_T^2 > cutoff, reweight.
        IF(MSTP(64).GE.2) THEN
          IF((1.-Z)*Q2B.LT.PARP(62)**2) GOTO 180
          ALPRAT=TEVB/(TEVB+LOG(1.-Z))
          IF(ALPRAT.LT.5.*PYR(0)) GOTO 180
          IF(ALPRAT.GT.5.) WTZ=WTZ*ALPRAT/5.
        ENDIF

C...Option with angular ordering requirement.
        IF(MSTP(62).GE.3) THEN
          THE2T=(4.*Z**2*Q2B)/(VINT(2)*(1.-Z)*XB**2)
          IF(THE2T.GT.THE2(JT)) GOTO 180
        ENDIF

C...Weighting with new structure functions.
        CALL RYSTFU(MINT(10+JT),XB,Q2REF,XFN)
        IF(KFLB.NE.21) XFBN=XFN(KFLB)
        IF(KFLB.EQ.21) XFBN=XFN(0)
        IF(XFBN.LT.1E-20) THEN
          IF(KFLA.EQ.KFLB) THEN
            TEVB=TEVBSV
            WTAP(KFLB)=0.
            GOTO 160
          ELSEIF(TEVBSV-TEVB.GT.0.2) THEN
            TEVB=0.5*(TEVBSV+TEVB)
            GOTO 190
          ELSE
            XFBN=1E-10
            IF(KFLB.NE.21) XFN(KFLB)=XFBN
            IF(KFLB.EQ.21) XFN(0)=XFBN
          ENDIF
        ENDIF
        DO 210 KFL=-MSTP(54),MSTP(54)
  210   XFB(KFL)=XFN(KFL)
        XA=XB/Z
        CALL RYSTFU(MINT(10+JT),XA,Q2REF,XFA)
        IF(KFLA.NE.21) XFAN=XFA(KFLA)
        IF(KFLA.EQ.21) XFAN=XFA(0)
        IF(XFAN.LT.1E-20) GOTO 160
        IF(KFLA.NE.21) WTSFA=WTSF(KFLA)
        IF(KFLA.EQ.21) WTSFA=WTSF(0)
        IF(WTZ*XFAN/XFBN.LT.PYR(0)*WTSFA) GOTO 160
      ENDIF

C...Define two hard scatterers in their CM-frame.
  220 IF(N.EQ.NS+2) THEN
        DQ2(JT)=Q2B
        DPLCM=SQRT((DSH+DQ2(1)+DQ2(2))**2-4D0*DQ2(1)*DQ2(2))/DSHR
        DO 240 JR=1,2
        I=NS+JR
        IF(JR.EQ.1) IPO=IPUS1
        IF(JR.EQ.2) IPO=IPUS2
        DO 230 J=1,5
        K(I,J)=0
        P(I,J)=0.
  230   V(I,J)=0.
        K(I,1)=14
        K(I,2)=KFLS(JR+2)
        K(I,4)=IPO
        K(I,5)=IPO
        P(I,3)=DPLCM*(-1)**(JR+1)
        P(I,4)=(DSH+DQ2(3-JR)-DQ2(JR))/DSHR
        P(I,5)=-SQRT(SNGL(DQ2(JR)))
        K(IPO,1)=14
        K(IPO,3)=I
        K(IPO,4)=MOD(K(IPO,4),MSTU(5))+MSTU(5)*I
  240   K(IPO,5)=MOD(K(IPO,5),MSTU(5))+MSTU(5)*I

C...Find maximum allowed mass of timelike parton.
      ELSEIF(N.GT.NS+2) THEN
        JR=3-JT
        DQ2(3)=Q2B
        DPC(1)=P(IS(1),4)
        DPC(2)=P(IS(2),4)
        DPC(3)=0.5*(ABS(P(IS(1),3))+ABS(P(IS(2),3)))
        DPD(1)=DSH+DQ2(JR)+DQ2(JT)
        DPD(2)=DSHZ+DQ2(JR)+DQ2(3)
        DPD(3)=SQRT(DPD(1)**2-4D0*DQ2(JR)*DQ2(JT))
        DPD(4)=SQRT(DPD(2)**2-4D0*DQ2(JR)*DQ2(3))
        IKIN=0
        IF(Q2S(JR).GE.(0.5*PARP(62))**2.AND.DPD(1)-DPD(3).GE.
     &  1D-10*DPD(1)) IKIN=1
        IF(IKIN.EQ.0) DMSMA=(DQ2(JT)/DBLE(ZS(JT))-DQ2(3))*(DSH/
     &  (DSH+DQ2(JT))-DSH/(DSHZ+DQ2(3)))
        IF(IKIN.EQ.1) DMSMA=(DPD(1)*DPD(2)-DPD(3)*DPD(4))/(2.*
     &  DQ2(JR))-DQ2(JT)-DQ2(3)

C...Generate timelike parton shower (if required).
        IT=N
        DO 250 J=1,5
        K(IT,J)=0
        P(IT,J)=0.
  250   V(IT,J)=0.
        K(IT,1)=3
        K(IT,2)=21
        IF(KFLB.EQ.21.AND.KFLS(JT+2).NE.21) K(IT,2)=-KFLS(JT+2)
        IF(KFLB.NE.21.AND.KFLS(JT+2).EQ.21) K(IT,2)=KFLB
        P(IT,5)=ULMASS(K(IT,2))
        IF(SNGL(DMSMA).LE.P(IT,5)**2) GOTO 100
        IF(MSTP(63).GE.1) THEN
          P(IT,4)=(DSHZ-DSH-P(IT,5)**2)/DSHR
          P(IT,3)=SQRT(P(IT,4)**2-P(IT,5)**2)
          IF(MSTP(63).EQ.1) THEN
            Q2TIM=DMSMA
          ELSEIF(MSTP(63).EQ.2) THEN
            Q2TIM=MIN(SNGL(DMSMA),PARP(71)*Q2S(JT))
          ELSE
C'''Here remains to introduce angular ordering in first branching.
            Q2TIM=DMSMA
          ENDIF
          CALL LUSHOW(IT,0,SQRT(Q2TIM))
          IF(N.GE.IT+1) P(IT,5)=P(IT+1,5)
        ENDIF

C...Reconstruct kinematics of branching: timelike parton shower.
        DMS=P(IT,5)**2
        IF(IKIN.EQ.0) DPT2=(DMSMA-DMS)*(DSHZ+DQ2(3))/(DSH+DQ2(JT))
        IF(IKIN.EQ.1) DPT2=(DMSMA-DMS)*(0.5*DPD(1)*DPD(2)+0.5*DPD(3)*
     &  DPD(4)-DQ2(JR)*(DQ2(JT)+DQ2(3)+DMS))/(4.*DSH*DPC(3)**2)
        IF(DPT2.LT.0.) GOTO 100
        DPB(1)=(0.5*DPD(2)-DPC(JR)*(DSHZ+DQ2(JR)-DQ2(JT)-DMS)/
     &  DSHR)/DPC(3)-DPC(3)
        P(IT,1)=SQRT(SNGL(DPT2))
        P(IT,3)=DPB(1)*(-1)**(JT+1)
        P(IT,4)=(DSHZ-DSH-DMS)/DSHR
        IF(N.GE.IT+1) THEN
          DPB(1)=SQRT(DPB(1)**2+DPT2)
          DPB(2)=SQRT(DPB(1)**2+DMS)
          DPB(3)=P(IT+1,3)
          DPB(4)=SQRT(DPB(3)**2+DMS)
          DBEZ=(DPB(4)*DPB(1)-DPB(3)*DPB(2))/(DPB(4)*DPB(2)-DPB(3)*
     &    DPB(1))
          CALL LUDBRB(IT+1,N,0.,0.,0D0,0D0,DBEZ)
          THE=ULANGL(P(IT,3),P(IT,1))
          CALL LUDBRB(IT+1,N,THE,0.,0D0,0D0,0D0)
        ENDIF

C...Reconstruct kinematics of branching: spacelike parton.
        DO 260 J=1,5
        K(N+1,J)=0
        P(N+1,J)=0.
  260   V(N+1,J)=0.
        K(N+1,1)=14
        K(N+1,2)=KFLB
        P(N+1,1)=P(IT,1)
        P(N+1,3)=P(IT,3)+P(IS(JT),3)
        P(N+1,4)=P(IT,4)+P(IS(JT),4)
        P(N+1,5)=-SQRT(SNGL(DQ2(3)))

C...Define colour flow of branching.
        K(IS(JT),3)=N+1
        K(IT,3)=N+1
        ID1=IT
        IF((K(N+1,2).GT.0.AND.K(N+1,2).NE.21.AND.K(ID1,2).GT.0.AND.
     &  K(ID1,2).NE.21).OR.(K(N+1,2).LT.0.AND.K(ID1,2).EQ.21).OR.
     &  (K(N+1,2).EQ.21.AND.K(ID1,2).EQ.21.AND.PYR(0).GT.0.5).OR.
     &  (K(N+1,2).EQ.21.AND.K(ID1,2).LT.0)) ID1=IS(JT)
        ID2=IT+IS(JT)-ID1
        K(N+1,4)=K(N+1,4)+ID1
        K(N+1,5)=K(N+1,5)+ID2
        K(ID1,4)=K(ID1,4)+MSTU(5)*(N+1)
        K(ID1,5)=K(ID1,5)+MSTU(5)*ID2
        K(ID2,4)=K(ID2,4)+MSTU(5)*ID1
        K(ID2,5)=K(ID2,5)+MSTU(5)*(N+1)
        N=N+1

C...Boost to new CM-frame.
        CALL LUDBRB(NS+1,N,0.,0.,-DBLE((P(N,1)+P(IS(JR),1))/(P(N,4)+
     &  P(IS(JR),4))),0D0,-DBLE((P(N,3)+P(IS(JR),3))/(P(N,4)+
     &  P(IS(JR),4))))
        IR=N+(JT-1)*(IS(1)-N)
        CALL LUDBRB(NS+1,N,-ULANGL(P(IR,3),P(IR,1)),PARU(2)*PYR(0),
     &  0D0,0D0,0D0)
      ENDIF

C...Save quantities, loop back.
      IS(JT)=N
      Q2S(JT)=Q2B
      DQ2(JT)=Q2B
      IF(MSTP(62).GE.3) THE2(JT)=THE2T
      DSH=DSHZ
      IF(Q2B.GE.(0.5*PARP(62))**2) THEN
        KFLS(JT+2)=KFLS(JT)
        KFLS(JT)=KFLA
        XS(JT)=XA
        ZS(JT)=Z
        DO 270 KFL=-25,25
  270   XFS(JT,KFL)=XFA(KFL)
        TEVS(JT)=TEVB
      ELSE
        IF(JT.EQ.1) IPU1=N
        IF(JT.EQ.2) IPU2=N
      ENDIF
      IF(N.GT.MSTU(4)-MSTU(32)-10) THEN
        CALL LUERRM(11,'(RYSSPA:) no more memory left in LUJETS')
        IF(MSTU(21).GE.1) N=NS
        IF(MSTU(21).GE.1) RETURN
      ENDIF
      IF(MAX(Q2S(1),Q2S(2)).GE.(0.5*PARP(62))**2.OR.N.LE.NS+1) GOTO 120

C...Boost hard scattering partons to frame of shower initiators.
      DO 280 J=1,3
  280 ROBO(J+2)=(P(NS+1,J)+P(NS+2,J))/(P(NS+1,4)+P(NS+2,4))
      K(N+2,1)=1
      DO 290 J=1,5
  290 P(N+2,J)=P(NS+1,J)
      ROBOT=ROBO(3)**2+ROBO(4)**2+ROBO(5)**2
      IF(ROBOT.GE.0.999999) THEN
        ROBOT=1.00001*SQRT(ROBOT)
        ROBO(3)=ROBO(3)/ROBOT
        ROBO(4)=ROBO(4)/ROBOT
        ROBO(5)=ROBO(5)/ROBOT
      ENDIF
      CALL LUDBRB(N+2,N+2,0.,0.,-DBLE(ROBO(3)),-DBLE(ROBO(4)),
     &-DBLE(ROBO(5)))
      ROBO(2)=ULANGL(P(N+2,1),P(N+2,2))
      ROBO(1)=ULANGL(P(N+2,3),SQRT(P(N+2,1)**2+P(N+2,2)**2))
      CALL LUDBRB(MINT(83)+5,NS,ROBO(1),ROBO(2),DBLE(ROBO(3)),
     &DBLE(ROBO(4)),DBLE(ROBO(5)))

C...Store user information. Reset Lambda value.
      K(IPU1,3)=MINT(83)+3
      K(IPU2,3)=MINT(83)+4
      DO 300 JT=1,2
      MINT(12+JT)=KFLS(JT)
  300 VINT(140+JT)=XS(JT)
      PARU(111)=ALAMS

      RETURN
      END

C*********************************************************************

      SUBROUTINE RYRESD

C...Allows resonances to decay (including parton showers for hadronic
C...channels).
      IMPLICIT DOUBLE PRECISION(D)
      PARAMETER (KSZJ=4000)

      
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      integer N,NPAD,K
      double precision P,V
      SAVE /PYJETS/


      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/LUDAT3/MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5)
      COMMON/RYSUBS/MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200)
      COMMON/RYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/RYINT1/MINT(400),VINT(400)
      COMMON/RYINT2/ISET(200),KFPR(200,2),COEF(200,20),ICOL(40,4,2)
      COMMON/RYINT4/WIDP(21:40,0:40),WIDE(21:40,0:40),WIDS(21:40,3)
      SAVE /LUDAT1/,/LUDAT2/,/LUDAT3/
      SAVE /RYSUBS/,/RYPARS/,/RYINT1/,/RYINT2/,/RYINT4/
      DIMENSION IREF(10,6),KDCY(2),KFL1(2),KFL2(2),KEQL(2),NSD(2),
     &ILIN(6),HGZ(2,3),COUP(6,4),CORL(2,2,2),PK(6,4),PKK(6,6),
     &CTHE(2),PHI(2),WDTP(0:40),WDTE(0:40,0:5),DBEZQQ(3)
      COMPLEX FGK,HA(6,6),HC(6,6)

C...The F, Xi and Xj functions of Gunion and Kunszt
C...(Phys. Rev. D33, 665, plus errata from the authors).
      FGK(I1,I2,I3,I4,I5,I6)=4.*HA(I1,I3)*HC(I2,I6)*(HA(I1,I5)*
     &HC(I1,I4)+HA(I3,I5)*HC(I3,I4))
      DIGK(DT,DU)=-4.*D34*D56+DT*(3.*DT+4.*DU)+DT**2*(DT*DU/(D34*D56)-
     &2.*(1./D34+1./D56)*(DT+DU)+2.*(D34/D56+D56/D34))
      DJGK(DT,DU)=8.*(D34+D56)**2-8.*(D34+D56)*(DT+DU)-6.*DT*DU-
     &2.*DT*DU*(DT*DU/(D34*D56)-2.*(1./D34+1./D56)*(DT+DU)+
     &2.*(D34/D56+D56/D34))

C...Some general constants.
      XW=PARU(102)
      SQMZ=PMAS(23,1)**2
      SQMW=PMAS(24,1)**2
      SH=VINT(44)

C...Define initial two objects.
      ISUB=MINT(1)
      IF(ISET(ISUB).EQ.1.OR.ISET(ISUB).EQ.3.OR.ISET(ISUB).EQ.5) THEN
        IREF(1,1)=MINT(84)+2+MIN(3,ISET(ISUB))
        IREF(1,2)=0
        IREF(1,3)=MINT(83)+6+MIN(3,ISET(ISUB))
        IREF(1,4)=0
      ELSEIF(ISET(ISUB).EQ.2.OR.ISET(ISUB).EQ.4) THEN
        IREF(1,1)=MINT(84)+1+ISET(ISUB)
        IREF(1,2)=MINT(84)+2+ISET(ISUB)
        IREF(1,3)=MINT(83)+5+ISET(ISUB)
        IREF(1,4)=MINT(83)+6+ISET(ISUB)
      ELSEIF(ISET(ISUB).EQ.6) THEN
        IREF(1,1)=MINT(84)+3
        IREF(1,2)=0
        IREF(1,3)=MINT(83)+7
        IREF(1,4)=0
      ENDIF

C...Check if initial resonance has been moved (in resonance + jet).
      IF(K(IREF(1,1),1).GT.10) THEN
        KFA=IABS(K(IREF(1,1),2))
        KDA=MOD(K(IREF(1,1),4),MSTU(4))
        IF(KFA.GE.23.AND.KFA.LE.40.AND.KDA.GT.1) IREF(1,1)=KDA
      ENDIF
      IF(IREF(1,2).GT.0) THEN
        IF(K(IREF(1,2),1).GT.10) THEN
          KFA=IABS(K(IREF(1,2),2))
          KDA=MOD(K(IREF(1,2),4),MSTU(4))
          IF(KFA.GE.23.AND.KFA.LE.40.AND.KDA.GT.1) IREF(1,2)=KDA
        ENDIF
      ENDIF
      IREF(1,5)=0
      IREF(1,6)=0

C...Loop over decay history.
      NP=1
      IP=0
  100 IP=IP+1
      NINH=0
      JTMAX=2
      IF(IP.EQ.1.AND.IREF(1,2).EQ.0) JTMAX=1
      NSAV=N

C...Start treatment of one or two resonances in parallel.
  110 N=NSAV
      DO 140 JT=1,JTMAX
      ID=IREF(IP,JT)
      KDCY(JT)=0
      KFL1(JT)=0
      KFL2(JT)=0
      KEQL(JT)=0
      NSD(JT)=ID
      IF(ID.EQ.0) GOTO 140
      KFA=IABS(K(ID,2))
      IF(KFA.LT.23.OR.KFA.GT.40) GOTO 140
      IF(K(ID,1).GT.10.OR.MDCY(KFA,1).EQ.0) GOTO 140

C...Select decay channel.
      IF(ISET(ISUB).NE.6) THEN
        IF(ISUB.EQ.1.OR.ISUB.EQ.22.OR.ISUB.EQ.141) MINT(61)=1
        CALL RYWIDT(KFA,P(ID,5)**2,WDTP,WDTE)
        IF(KCHG(KFA,3).EQ.0) THEN
          IPM=2
        ELSE
          IPM=(5+ISIGN(1,K(ID,2)))/2
        ENDIF
        KFB=0
        IF(JTMAX.EQ.2) KFB=IABS(K(IREF(IP,3-JT),2))
        WDTE0S=WDTE(0,1)+WDTE(0,IPM)+WDTE(0,4)
        IF(KFB.EQ.KFA) WDTE0S=WDTE0S+WDTE(0,5)
        IF(WDTE0S.LE.0.) GOTO 140
        RKFL=WDTE0S*PYR(0)
        IDL=0
  120   IDL=IDL+1
        IDC=IDL+MDCY(KFA,2)-1
        RKFL=RKFL-(WDTE(IDL,1)+WDTE(IDL,IPM)+WDTE(IDL,4))
        IF(KFB.EQ.KFA) RKFL=RKFL-WDTE(IDL,5)
        IF(IDL.LT.MDCY(KFA,3).AND.RKFL.GT.0.) GOTO 120
      ELSE
        IDC=MINT(35)
      ENDIF

C...Read out and classify decay channel chosen.
      KFL1(JT)=KFDP(IDC,1)*ISIGN(1,K(ID,2))
      KFC1A=IABS(KFL1(JT))
      IF(KFC1A.GT.100) KFC1A=LUCOMP(KFC1A)
      IF(KCHG(KFC1A,3).EQ.0) KFL1(JT)=IABS(KFL1(JT))
      KFL2(JT)=KFDP(IDC,2)*ISIGN(1,K(ID,2))
      KFC2A=IABS(KFL2(JT))
      IF(KFC2A.GT.100) KFC2A=LUCOMP(KFC2A)
      IF(KCHG(KFC2A,3).EQ.0) KFL2(JT)=IABS(KFL2(JT))
      KDCY(JT)=2
      IF(IABS(KFL1(JT)).LE.10.OR.KFL1(JT).EQ.21) KDCY(JT)=1
      IF(IABS(KFL1(JT)).GE.23.AND.IABS(KFL1(JT)).LE.40) KDCY(JT)=3
      KEQL(JT)=MDME(IDC,1)
      NSD(JT)=N
      HGZ(JT,1)=VINT(111)
      HGZ(JT,2)=VINT(112)
      HGZ(JT,3)=VINT(114)

C...Select masses and check that mass sum not too large.
      IF(MSTP(42).LE.0.OR.(PMAS(KFC1A,2).LT.PARP(41).AND.
     &PMAS(KFC2A,2).LT.PARP(41))) THEN
        P(N+1,5)=PMAS(KFC1A,1)
        P(N+2,5)=PMAS(KFC2A,1)
        IF(P(N+1,5)+P(N+2,5)+PARJ(64).GT.P(ID,5)) THEN
          CALL LUERRM(13,'(RYRESD:) daughter masses too large')
          MINT(51)=1
          RETURN
        ENDIF
      ELSEIF(IP.EQ.1) THEN
        CALL RYOFSH(2,KFA,KFL1(JT),KFL2(JT),P(ID,5),P(N+1,5),P(N+2,5))
        IF(MINT(51).EQ.1) RETURN
      ELSE
  130   P(N+1,5)=ULMASS(KFL1(JT))
        P(N+2,5)=ULMASS(KFL2(JT))
        IF(P(N+1,5)+P(N+2,5)+PARJ(64).GT.P(ID,5)) GOTO 130
      ENDIF

C...Fill decay products, prepared for parton showers for quarks.
C...Special case, done by hand, for leptoquark.
      MSTU(10)=1
      IF(KFA.EQ.39) THEN
        MSTU(19)=1
        CALL LU2ENT(N+1,KFL1(JT),KFL2(JT),P(ID,5))
        K(N-1,1)=3
        K(N,1)=1
        ISID=4
        IF(K(ID,2).LT.0) ISID=5
        K(ID,ISID)=K(ID,ISID)+(N-1)
        K(N-1,ISID)=MSTU(5)*ID
      ELSEIF(KDCY(JT).EQ.1) THEN
        CALL LU2ENT(-(N+1),KFL1(JT),KFL2(JT),P(ID,5))
      ELSE
        CALL LU2ENT(N+1,KFL1(JT),KFL2(JT),P(ID,5))
      ENDIF
      MSTU(10)=2
  140 IF(KFA.GE.23.AND.KFA.LE.40.AND.KFL1(JT).EQ.0) NINH=NINH+1

C...Check for allowed combinations. Skip if no decays.
      IF(ISET(ISUB).EQ.6) GOTO 290
      IF(JTMAX.EQ.2) THEN
        IF(KEQL(1).EQ.4.AND.KEQL(2).EQ.4) GOTO 110
        IF(KEQL(1).EQ.5.AND.KEQL(2).EQ.5) GOTO 110
      ENDIF
      IF(JTMAX.EQ.1.AND.KDCY(1).EQ.0) GOTO 350
      IF(JTMAX.EQ.2.AND.KDCY(1).EQ.0.AND.KDCY(2).EQ.0) GOTO 350

C...Order incoming partons and outgoing resonances.
      IF(JTMAX.EQ.2.AND.MSTP(47).GE.1.AND.NINH.EQ.0) THEN
        ILIN(1)=MINT(84)+1
        IF(K(MINT(84)+1,2).GT.0) ILIN(1)=MINT(84)+2
        IF(K(ILIN(1),2).EQ.21) ILIN(1)=2*MINT(84)+3-ILIN(1)
        ILIN(2)=2*MINT(84)+3-ILIN(1)
        IMIN=1
        IF(IREF(IP,5).EQ.25.OR.IREF(IP,5).EQ.35.OR.IREF(IP,5).
     &  EQ.36) IMIN=3
        IMAX=2
        IORD=1
        IF(K(IREF(IP,1),2).EQ.23) IORD=2
        IF(K(IREF(IP,1),2).EQ.24.AND.K(IREF(IP,2),2).EQ.-24) IORD=2
        IAKIPD=IABS(K(IREF(IP,IORD),2))
        IF(IAKIPD.EQ.25.OR.IAKIPD.EQ.35.OR.IAKIPD.EQ.36) IORD=3-IORD
        IF(KDCY(IORD).EQ.0) IORD=3-IORD

C...Order decay products of resonances.
        DO 150 JT=IORD,3-IORD,3-2*IORD
        IF(KDCY(JT).EQ.0) THEN
          ILIN(IMAX+1)=NSD(JT)
          IMAX=IMAX+1
        ELSEIF(K(NSD(JT)+1,2).GT.0) THEN
          ILIN(IMAX+1)=N+2*JT-1
          ILIN(IMAX+2)=N+2*JT
          IMAX=IMAX+2
          K(N+2*JT-1,2)=K(NSD(JT)+1,2)
          K(N+2*JT,2)=K(NSD(JT)+2,2)
        ELSE
          ILIN(IMAX+1)=N+2*JT
          ILIN(IMAX+2)=N+2*JT-1
          IMAX=IMAX+2
          K(N+2*JT-1,2)=K(NSD(JT)+1,2)
          K(N+2*JT,2)=K(NSD(JT)+2,2)
        ENDIF
  150   CONTINUE

C...Find charge, isospin, left- and righthanded couplings.
        DO 170 I=IMIN,IMAX
        DO 160 J=1,4
  160   COUP(I,J)=0.
        KFA=IABS(K(ILIN(I),2))
        IF(KFA.EQ.0.OR.KFA.GT.20) GOTO 170
        COUP(I,1)=KCHG(KFA,1)/3.
        COUP(I,2)=(-1)**MOD(KFA,2)
        COUP(I,4)=-2.*COUP(I,1)*XW
        COUP(I,3)=COUP(I,2)+COUP(I,4)
  170   CONTINUE

C...Full propagator dependence and flavour correlations for 2 gamma*/Z.
        IF(ISUB.EQ.22) THEN
          DO 180 I=3,5,2
          I1=IORD
          IF(I.EQ.5) I1=3-IORD
          DO 180 J1=1,2
          DO 180 J2=1,2
  180     CORL(I/2,J1,J2)=COUP(1,1)**2*HGZ(I1,1)*COUP(I,1)**2/16.+
     &    COUP(1,1)*COUP(1,J1+2)*HGZ(I1,2)*COUP(I,1)*COUP(I,J2+2)/4.+
     &    COUP(1,J1+2)**2*HGZ(I1,3)*COUP(I,J2+2)**2
          COWT12=(CORL(1,1,1)+CORL(1,1,2))*(CORL(2,1,1)+CORL(2,1,2))+
     &    (CORL(1,2,1)+CORL(1,2,2))*(CORL(2,2,1)+CORL(2,2,2))
          COMX12=(CORL(1,1,1)+CORL(1,1,2)+CORL(1,2,1)+CORL(1,2,2))*
     &    (CORL(2,1,1)+CORL(2,1,2)+CORL(2,2,1)+CORL(2,2,2))
          IF(COWT12.LT.PYR(0)*COMX12) GOTO 110
        ENDIF
      ENDIF

C...Select angular orientation type - Z'/W' only.
      MZPWP=0
      IF(ISUB.EQ.141) THEN
        IF(PYR(0).LT.PARU(130)) MZPWP=1
        IF(IP.EQ.2.AND.IABS(K(IREF(2,1),2)).EQ.37) MZPWP=2
        IAKIR=IABS(K(IREF(2,2),2))
        IF(IP.EQ.2.AND.(IAKIR.EQ.25.OR.IAKIR.EQ.35.OR.IAKIR.EQ.36))
     &  MZPWP=2
        IF(IP.GE.3) MZPWP=2
      ELSEIF(ISUB.EQ.142) THEN
        IF(PYR(0).LT.PARU(136)) MZPWP=1
        IAKIR=IABS(K(IREF(2,2),2))
        IF(IP.EQ.2.AND.(IAKIR.EQ.25.OR.IAKIR.EQ.35.OR.IAKIR.EQ.36))
     &  MZPWP=2
        IF(IP.GE.3) MZPWP=2
      ENDIF

C...Select random angles (begin of weighting procedure).
  190 DO 200 JT=1,JTMAX
      IF(KDCY(JT).EQ.0) GOTO 200
      IF(JTMAX.EQ.1) THEN
        CTHE(JT)=VINT(13)+(VINT(33)-VINT(13)+VINT(34)-VINT(14))*PYR(0)
        IF(CTHE(JT).GT.VINT(33)) CTHE(JT)=CTHE(JT)+VINT(14)-VINT(33)
        PHI(JT)=VINT(24)
      ELSE
        CTHE(JT)=2.*PYR(0)-1.
        PHI(JT)=PARU(2)*PYR(0)
      ENDIF
  200 CONTINUE

      IF(JTMAX.EQ.2.AND.MSTP(47).GE.1.AND.NINH.EQ.0) THEN
C...Construct massless four-vectors.
        DO 210 I=N+1,N+4
        K(I,1)=1
        DO 210 J=1,5
  210   P(I,J)=0.
        DO 220 JT=1,JTMAX
        IF(KDCY(JT).EQ.0) GOTO 220
        ID=IREF(IP,JT)
        P(N+2*JT-1,3)=0.5*P(ID,5)
        P(N+2*JT-1,4)=0.5*P(ID,5)
        P(N+2*JT,3)=-0.5*P(ID,5)
        P(N+2*JT,4)=0.5*P(ID,5)
        CALL LUDBRB(N+2*JT-1,N+2*JT,ACOS(CTHE(JT)),PHI(JT),DBLE(P(ID,1)/
     &  P(ID,4)),DBLE(P(ID,2)/P(ID,4)),DBLE(P(ID,3)/P(ID,4)))
  220   CONTINUE

C...Store incoming and outgoing momenta, with random rotation to
C...avoid accidental zeroes in HA expressions.
        DO 230 I=1,IMAX
        K(N+4+I,1)=1
        P(N+4+I,4)=SQRT(P(ILIN(I),1)**2+P(ILIN(I),2)**2+P(ILIN(I),3)**2+
     &  P(ILIN(I),5)**2)
        P(N+4+I,5)=P(ILIN(I),5)
        DO 230 J=1,3
  230   P(N+4+I,J)=P(ILIN(I),J)
  240   THERR=ACOS(2.*PYR(0)-1.)
        PHIRR=PARU(2)*PYR(0)
        CALL LUDBRB(N+5,N+4+IMAX,THERR,PHIRR,0D0,0D0,0D0)
        DO 250 I=1,IMAX
        IF(P(N+4+I,1)**2+P(N+4+I,2)**2.LT.1E-4*P(N+4+I,4)**2) GOTO 240
        DO 250 J=1,4
  250   PK(I,J)=P(N+4+I,J)

C...Calculate internal products.
        IF(ISUB.EQ.22.OR.ISUB.EQ.23.OR.ISUB.EQ.25.OR.ISUB.EQ.141.OR.
     &  ISUB.EQ.142) THEN
          DO 260 I1=IMIN,IMAX-1
          DO 260 I2=I1+1,IMAX
          HA(I1,I2)=SQRT((PK(I1,4)-PK(I1,3))*(PK(I2,4)+PK(I2,3))/
     &    (1E-20+PK(I1,1)**2+PK(I1,2)**2))*CMPLX(PK(I1,1),PK(I1,2))-
     &    SQRT((PK(I1,4)+PK(I1,3))*(PK(I2,4)-PK(I2,3))/
     &    (1E-20+PK(I2,1)**2+PK(I2,2)**2))*CMPLX(PK(I2,1),PK(I2,2))
          HC(I1,I2)=CONJG(HA(I1,I2))
          IF(I1.LE.2) HA(I1,I2)=CMPLX(0.,1.)*HA(I1,I2)
          IF(I1.LE.2) HC(I1,I2)=CMPLX(0.,1.)*HC(I1,I2)
          HA(I2,I1)=-HA(I1,I2)
  260     HC(I2,I1)=-HC(I1,I2)
        ENDIF
        DO 270 I=1,2
        DO 270 J=1,4
  270   PK(I,J)=-PK(I,J)
        DO 280 I1=IMIN,IMAX-1
        DO 280 I2=I1+1,IMAX
        PKK(I1,I2)=2.*(PK(I1,4)*PK(I2,4)-PK(I1,1)*PK(I2,1)-
     &  PK(I1,2)*PK(I2,2)-PK(I1,3)*PK(I2,3))
  280   PKK(I2,I1)=PKK(I1,I2)
      ENDIF

      IF(MSTP(47).LE.0.OR.NINH.NE.0) THEN
C...Isotropic decay selected by user.
        WT=1.
        WTMAX=1.

      ELSEIF(IREF(IP,5).EQ.25.OR.IREF(IP,5).EQ.35.OR.
     &IREF(IP,5).EQ.36) THEN
C...Angular weight for H0 -> Z0 + Z0 or W+ + W- -> 4 quarks/leptons.
        WT=16.*PKK(3,5)*PKK(4,6)
        IF(IP.EQ.1) WTMAX=SH**2
        IF(IP.GE.2) WTMAX=P(IREF(IP,6),5)**4

      ELSEIF(ISUB.EQ.1) THEN
C...Angular weight for gamma*/Z0 -> 2 quarks/leptons.
        EI=KCHG(IABS(MINT(15)),1)/3.
        AI=SIGN(1.,EI+0.1)
        VI=AI-4.*EI*XW
        EF=KCHG(IABS(KFL1(1)),1)/3.
        AF=SIGN(1.,EF+0.1)
        VF=AF-4.*EF*XW
        ASYM=2.*(EI*AI*VINT(112)*EF*AF+4.*VI*AI*VINT(114)*VF*AF)/
     &  (EI**2*VINT(111)*EF**2+EI*VI*VINT(112)*EF*VF+
     &  (VI**2+AI**2)*VINT(114)*(VF**2+AF**2))
        WT=1.+ASYM*CTHE(1)*ISIGN(1,MINT(15)*KFL1(1))+CTHE(1)**2
        WTMAX=2.+ABS(ASYM)

      ELSEIF(ISUB.EQ.2) THEN
C...Angular weight for W+/- -> 2 quarks/leptons.
        WT=(1.+CTHE(1)*ISIGN(1,MINT(15)*KFL1(1)))**2
        WTMAX=4.

      ELSEIF(ISUB.EQ.15.OR.ISUB.EQ.19) THEN
C...Angular weight for f + f~ -> gluon/gamma + Z0 ->
C...-> gluon/gamma + 2 quarks/leptons.
        WT=((COUP(1,3)*COUP(3,3))**2+(COUP(1,4)*COUP(3,4))**2)*
     &  (PKK(1,3)**2+PKK(2,4)**2)+((COUP(1,3)*COUP(3,4))**2+
     &  (COUP(1,4)*COUP(3,3))**2)*(PKK(1,4)**2+PKK(2,3)**2)
        WTMAX=(COUP(1,3)**2+COUP(1,4)**2)*(COUP(3,3)**2+COUP(3,4)**2)*
     &  ((PKK(1,3)+PKK(1,4))**2+(PKK(2,3)+PKK(2,4))**2)

      ELSEIF(ISUB.EQ.16.OR.ISUB.EQ.20) THEN
C...Angular weight for f + f~' -> gluon/gamma + W+/- ->
C...-> gluon/gamma + 2 quarks/leptons.
        WT=PKK(1,3)**2+PKK(2,4)**2
        WTMAX=(PKK(1,3)+PKK(1,4))**2+(PKK(2,3)+PKK(2,4))**2

      ELSEIF(ISUB.EQ.22) THEN
C...Angular weight for f + f~ -> Z0 + Z0 -> 4 quarks/leptons.
        S34=P(IREF(IP,IORD),5)**2
        S56=P(IREF(IP,3-IORD),5)**2
        TI=PKK(1,3)+PKK(1,4)+S34
        UI=PKK(1,5)+PKK(1,6)+S56
        FGK135=ABS(FGK(1,2,3,4,5,6)/TI+FGK(1,2,5,6,3,4)/UI)**2
        FGK145=ABS(FGK(1,2,4,3,5,6)/TI+FGK(1,2,5,6,4,3)/UI)**2
        FGK136=ABS(FGK(1,2,3,4,6,5)/TI+FGK(1,2,6,5,3,4)/UI)**2
        FGK146=ABS(FGK(1,2,4,3,6,5)/TI+FGK(1,2,6,5,4,3)/UI)**2
        FGK253=ABS(FGK(2,1,5,6,3,4)/TI+FGK(2,1,3,4,5,6)/UI)**2
        FGK263=ABS(FGK(2,1,6,5,3,4)/TI+FGK(2,1,3,4,6,5)/UI)**2
        FGK254=ABS(FGK(2,1,5,6,4,3)/TI+FGK(2,1,4,3,5,6)/UI)**2
        FGK264=ABS(FGK(2,1,6,5,4,3)/TI+FGK(2,1,4,3,6,5)/UI)**2
        WT=
     &  CORL(1,1,1)*CORL(2,1,1)*FGK135+CORL(1,1,2)*CORL(2,1,1)*FGK145+
     &  CORL(1,1,1)*CORL(2,1,2)*FGK136+CORL(1,1,2)*CORL(2,1,2)*FGK146+
     &  CORL(1,2,1)*CORL(2,2,1)*FGK253+CORL(1,2,2)*CORL(2,2,1)*FGK263+
     &  CORL(1,2,1)*CORL(2,2,2)*FGK254+CORL(1,2,2)*CORL(2,2,2)*FGK264
        WTMAX=16.*((CORL(1,1,1)+CORL(1,1,2))*(CORL(2,1,1)+CORL(2,1,2))+
     &  (CORL(1,2,1)+CORL(1,2,2))*(CORL(2,2,1)+CORL(2,2,2)))*S34*S56*
     &  ((TI**2+UI**2+2.*SH*(S34+S56))/(TI*UI)-S34*S56*(1./TI**2+
     &  1./UI**2))

      ELSEIF(ISUB.EQ.23) THEN
C...Angular weight for f + f~' -> Z0 + W+/- -> 4 quarks/leptons.
        D34=P(IREF(IP,IORD),5)**2
        D56=P(IREF(IP,3-IORD),5)**2
        DT=PKK(1,3)+PKK(1,4)+D34
        DU=PKK(1,5)+PKK(1,6)+D56
        CAWZ=COUP(2,3)/SNGL(DT)-2.*(1.-XW)*COUP(1,2)/(SH-SQMW)
        CBWZ=COUP(1,3)/SNGL(DU)+2.*(1.-XW)*COUP(1,2)/(SH-SQMW)
        FGK135=ABS(CAWZ*FGK(1,2,3,4,5,6)+CBWZ*FGK(1,2,5,6,3,4))
        FGK136=ABS(CAWZ*FGK(1,2,3,4,6,5)+CBWZ*FGK(1,2,6,5,3,4))
        WT=(COUP(5,3)*FGK135)**2+(COUP(5,4)*FGK136)**2
        WTMAX=4.*D34*D56*(COUP(5,3)**2+COUP(5,4)**2)*(CAWZ**2*
     &  DIGK(DT,DU)+CBWZ**2*DIGK(DU,DT)+CAWZ*CBWZ*DJGK(DT,DU))

      ELSEIF(ISUB.EQ.24.OR.ISUB.EQ.171.OR.ISUB.EQ.176) THEN
C...Angular weight for f + f~ -> Z0 + H0 -> 2 quarks/leptons + H0
C...(or H'0, or A0).
        WT=((COUP(1,3)*COUP(3,3))**2+(COUP(1,4)*COUP(3,4))**2)*
     &  PKK(1,3)*PKK(2,4)+((COUP(1,3)*COUP(3,4))**2+(COUP(1,4)*
     &  COUP(3,3))**2)*PKK(1,4)*PKK(2,3)
        WTMAX=(COUP(1,3)**2+COUP(1,4)**2)*(COUP(3,3)**2+COUP(3,4)**2)*
     &  (PKK(1,3)+PKK(1,4))*(PKK(2,3)+PKK(2,4))

      ELSEIF(ISUB.EQ.25) THEN
C...Angular weight for f + f~ -> W+ + W- -> 4 quarks/leptons.
        D34=P(IREF(IP,IORD),5)**2
        D56=P(IREF(IP,3-IORD),5)**2
        DT=PKK(1,3)+PKK(1,4)+D34
        DU=PKK(1,5)+PKK(1,6)+D56
        CDWW=(COUP(1,3)*SQMZ/(SH-SQMZ)+COUP(1,2))/SH
        CAWW=CDWW+0.5*(COUP(1,2)+1.)/SNGL(DT)
        CBWW=CDWW+0.5*(COUP(1,2)-1.)/SNGL(DU)
        CCWW=COUP(1,4)*SQMZ/(SH-SQMZ)/SH
        FGK135=ABS(CAWW*FGK(1,2,3,4,5,6)-CBWW*FGK(1,2,5,6,3,4))
        FGK253=ABS(FGK(2,1,5,6,3,4)-FGK(2,1,3,4,5,6))
        WT=FGK135**2+(CCWW*FGK253)**2
        WTMAX=4.*D34*D56*(CAWW**2*DIGK(DT,DU)+CBWW**2*DIGK(DU,DT)-CAWW*
     &  CBWW*DJGK(DT,DU)+CCWW**2*(DIGK(DT,DU)+DIGK(DU,DT)-DJGK(DT,DU)))

      ELSEIF(ISUB.EQ.26.OR.ISUB.EQ.172.OR.ISUB.EQ.177) THEN
C...Angular weight for f + f~' -> W+/- + H0 -> 2 quarks/leptons + H0
C...(or H'0, or A0).
        WT=PKK(1,3)*PKK(2,4)
        WTMAX=(PKK(1,3)+PKK(1,4))*(PKK(2,3)+PKK(2,4))

      ELSEIF(ISUB.EQ.30.OR.ISUB.EQ.35) THEN
C...Angular weight for f + g -> f + Z0 -> f + 2 quarks/leptons
C...and for f + gamma -> f + Z0 -> f + 2 quarks/leptons.
        IF(K(ILIN(1),2).GT.0) WT=((COUP(1,3)*COUP(3,3))**2+
     &  (COUP(1,4)*COUP(3,4))**2)*(PKK(1,4)**2+PKK(3,5)**2)+
     &  ((COUP(1,3)*COUP(3,4))**2+(COUP(1,4)*COUP(3,3))**2)*
     &  (PKK(1,3)**2+PKK(4,5)**2)
        IF(K(ILIN(1),2).LT.0) WT=((COUP(1,3)*COUP(3,3))**2+
     &  (COUP(1,4)*COUP(3,4))**2)*(PKK(1,3)**2+PKK(4,5)**2)+
     &  ((COUP(1,3)*COUP(3,4))**2+(COUP(1,4)*COUP(3,3))**2)*
     &  (PKK(1,4)**2+PKK(3,5)**2)
        WTMAX=(COUP(1,3)**2+COUP(1,4)**2)*(COUP(3,3)**2+COUP(3,4)**2)*
     &  ((PKK(1,3)+PKK(1,4))**2+(PKK(3,5)+PKK(4,5))**2)

      ELSEIF(ISUB.EQ.31) THEN
C...Angular weight for f + g -> f' + W+/- -> f' + 2 quarks/leptons.
        IF(K(ILIN(1),2).GT.0) WT=PKK(1,4)**2+PKK(3,5)**2
        IF(K(ILIN(1),2).LT.0) WT=PKK(1,3)**2+PKK(4,5)**2
        WTMAX=(PKK(1,3)+PKK(1,4))**2+(PKK(3,5)+PKK(4,5))**2

      ELSEIF(ISUB.EQ.71.OR.ISUB.EQ.72.OR.ISUB.EQ.73.OR.ISUB.EQ.76.OR.
     &ISUB.EQ.77) THEN
C...Angular weight for V_L1 + V_L2 -> V_L3 + V_L4 (V = Z/W).
        WT=16.*PKK(3,5)*PKK(4,6)
        WTMAX=SH**2

      ELSEIF(ISUB.EQ.141) THEN
        IF(IP.EQ.1.AND.(KDCY(1).EQ.1.OR.KDCY(1).EQ.2)) THEN
C...Angular weight for f + f~ -> gamma*/Z0/Z'0 -> 2 quarks/leptons.
C...Couplings of incoming flavour.
          KFAI=IABS(MINT(15))
          EI=KCHG(KFAI,1)/3.
          AI=SIGN(1.,EI+0.1)
          VI=AI-4.*EI*XW
          KFAIC=1
          IF(KFAI.LE.10.AND.MOD(KFAI,2).EQ.0) KFAIC=2
          IF(KFAI.GT.10.AND.MOD(KFAI,2).NE.0) KFAIC=3
          IF(KFAI.GT.10.AND.MOD(KFAI,2).EQ.0) KFAIC=4
          VPI=PARU(119+2*KFAIC)
          API=PARU(120+2*KFAIC)
C...Couplings of final flavour.
          KFAF=IABS(KFL1(1))
          EF=KCHG(KFAF,1)/3.
          AF=SIGN(1.,EF+0.1)
          VF=AF-4.*EF*XW
          KFAFC=1
          IF(KFAF.LE.10.AND.MOD(KFAF,2).EQ.0) KFAFC=2
          IF(KFAF.GT.10.AND.MOD(KFAF,2).NE.0) KFAFC=3
          IF(KFAF.GT.10.AND.MOD(KFAF,2).EQ.0) KFAFC=4
          VPF=PARU(119+2*KFAFC)
          APF=PARU(120+2*KFAFC)
C...Asymmetry and weight.
          ASYM=2.*(EI*AI*VINT(112)*EF*AF+EI*API*VINT(113)*EF*APF+
     &    4.*VI*AI*VINT(114)*VF*AF+(VI*API+VPI*AI)*VINT(115)*
     &    (VF*APF+VPF*AF)+4.*VPI*API*VINT(116)*VPF*APF)/
     &    (EI**2*VINT(111)*EF**2+EI*VI*VINT(112)*EF*VF+
     &    EI*VPI*VINT(113)*EF*VPF+(VI**2+AI**2)*VINT(114)*
     &    (VF**2+AF**2)+(VI*VPI+AI*API)*VINT(115)*(VF*VPF+AF*APF)+
     &    (VPI**2+API**2)*VINT(116)*(VPF**2+APF**2))
          WT=1.+ASYM*CTHE(1)*ISIGN(1,MINT(15)*KFL1(1))+CTHE(1)**2
          WTMAX=2.+ABS(ASYM)
        ELSEIF(IP.EQ.1.AND.IABS(KFL1(1)).EQ.24) THEN
C...Angular weight for f + f~ -> Z' -> W+ + W-.
          RM1=P(NSD(1)+1,5)**2/SH
          RM2=P(NSD(1)+2,5)**2/SH
          CCOS2=-(1./16.)*((1.-RM1-RM2)**2-4.*RM1*RM2)*
     &    (1.-2.*RM1-2.*RM2+RM1**2+RM2**2+10.*RM1*RM2)
          CFLAT=-CCOS2+0.5*(RM1+RM2)*(1.-2.*RM1-2.*RM2+(RM2-RM1)**2)
          WT=CFLAT+CCOS2*CTHE(1)**2
          WTMAX=CFLAT+MAX(0.,CCOS2)
        ELSEIF(IP.EQ.1.AND.(KFL1(1).EQ.25.OR.KFL1(1).EQ.35.OR.
     &  IABS(KFL1(1)).EQ.37)) THEN
C...Angular weight for f + f~ -> Z' -> H0 + A0, H'0 + A0, H+ + H-.
          WT=1.-CTHE(1)**2
          WTMAX=1.
        ELSEIF(IP.EQ.1.AND.KDCY(1).EQ.3) THEN
C...Angular weight for f + f~ -> Z' -> Z0 + H0.
          RM1=P(NSD(1)+1,5)**2/SH
          RM2=P(NSD(1)+2,5)**2/SH
          FLAM2=MAX(0.,(1.-RM1-RM2)**2-4.*RM1*RM2)
          WT=1.+FLAM2*(1.-CTHE(1)**2)/(8.*RM1)
          WTMAX=1.+FLAM2/(8.*RM1)
        ELSEIF(MZPWP.EQ.0) THEN
C...Angular weight for f + f~ -> Z' -> W+ + W- -> 4 quarks/leptons
C...(W:s like if intermediate Z).
          D34=P(IREF(IP,IORD),5)**2
          D56=P(IREF(IP,3-IORD),5)**2
          DT=PKK(1,3)+PKK(1,4)+D34
          DU=PKK(1,5)+PKK(1,6)+D56
          FGK135=ABS(FGK(1,2,3,4,5,6)-FGK(1,2,5,6,3,4))
          FGK253=ABS(FGK(2,1,5,6,3,4)-FGK(2,1,3,4,5,6))
          WT=(COUP(1,3)*FGK135)**2+(COUP(1,4)*FGK253)**2
          WTMAX=4.*D34*D56*(COUP(1,3)**2+COUP(1,4)**2)*
     &    (DIGK(DT,DU)+DIGK(DU,DT)-DJGK(DT,DU))
        ELSEIF(MZPWP.EQ.1) THEN
C...Angular weight for f + f~ -> Z' -> W+ + W- -> 4 quarks/leptons
C...(W:s approximately longitudinal, like if intermediate H).
          WT=16.*PKK(3,5)*PKK(4,6)
          WTMAX=SH**2
        ELSE
C...Angular weight for f + f~ -> Z' -> H+ + H-, Z0 + H0, H0 + A0,
C...H'0 + A0 -> 4 quarks/leptons.
          WT=1.
          WTMAX=1.
        ENDIF

      ELSEIF(ISUB.EQ.142) THEN
        IF(IP.EQ.1.AND.(KDCY(1).EQ.1.OR.KDCY(1).EQ.2)) THEN
C...Angular weight for f + f~' -> W'+/- -> 2 quarks/leptons.
          KFAI=IABS(MINT(15))
          KFAIC=1
          IF(KFAI.GT.10) KFAIC=2
          VI=PARU(129+2*KFAIC)
          AI=PARU(130+2*KFAIC)
          KFAF=IABS(KFL1(1))
          KFAFC=1
          IF(KFAF.GT.10) KFAFC=2
          VF=PARU(129+2*KFAFC)
          AF=PARU(130+2*KFAFC)
          ASYM=8.*VI*AI*VF*AF/((VI**2+AI**2)*(VF**2+AF**2))
          WT=1.+ASYM*CTHE(1)*ISIGN(1,MINT(15)*KFL1(1))+CTHE(1)**2
          WTMAX=2.+ABS(ASYM)
        ELSEIF(IP.EQ.1.AND.IABS(KFL2(1)).EQ.23) THEN
C...Angular weight for f + f~' -> W'+/- -> W+/- + Z0.
          RM1=P(NSD(1)+1,5)**2/SH
          RM2=P(NSD(1)+2,5)**2/SH
          CCOS2=-(1./16.)*((1.-RM1-RM2)**2-4.*RM1*RM2)*
     &    (1.-2.*RM1-2.*RM2+RM1**2+RM2**2+10.*RM1*RM2)
          CFLAT=-CCOS2+0.5*(RM1+RM2)*(1.-2.*RM1-2.*RM2+(RM2-RM1)**2)
          WT=CFLAT+CCOS2*CTHE(1)**2
          WTMAX=CFLAT+MAX(0.,CCOS2)
        ELSEIF(IP.EQ.1.AND.KDCY(1).EQ.3) THEN
C...Angular weight for f + f~ -> W'+/- -> W+/- + H0.
          RM1=P(NSD(1)+1,5)**2/SH
          RM2=P(NSD(1)+2,5)**2/SH
          FLAM2=MAX(0.,(1.-RM1-RM2)**2-4.*RM1*RM2)
          WT=1.+FLAM2*(1.-CTHE(1)**2)/(8.*RM1)
          WTMAX=1.+FLAM2/(8.*RM1)
        ELSEIF(MZPWP.EQ.0) THEN
C...Angular weight for f + f~' -> W' -> W + Z0 -> 4 quarks/leptons
C...(W/Z like if intermediate W).
          D34=P(IREF(IP,IORD),5)**2
          D56=P(IREF(IP,3-IORD),5)**2
          DT=PKK(1,3)+PKK(1,4)+D34
          DU=PKK(1,5)+PKK(1,6)+D56
          FGK135=ABS(FGK(1,2,3,4,5,6)-FGK(1,2,5,6,3,4))
          FGK136=ABS(FGK(1,2,3,4,6,5)-FGK(1,2,6,5,3,4))
          WT=(COUP(5,3)*FGK135)**2+(COUP(5,4)*FGK136)**2
          WTMAX=4.*D34*D56*(COUP(5,3)**2+COUP(5,4)**2)*
     &    (DIGK(DT,DU)+DIGK(DU,DT)-DJGK(DT,DU))
        ELSEIF(MZPWP.EQ.1) THEN
C...Angular weight for f + f~' -> W' -> W + Z0 -> 4 quarks/leptons
C...(W/Z approximately longitudinal, like if intermediate H).
          WT=16.*PKK(3,5)*PKK(4,6)
          WTMAX=SH**2
        ELSE
C...Angular weight for f + f~ -> W' -> W + H0 -> whatever.
          WT=1.
          WTMAX=1.
        ENDIF

      ELSEIF(ISUB.EQ.145.OR.ISUB.EQ.162.OR.ISUB.EQ.163.OR.ISUB.EQ.164)
     &THEN
C...Isotropic decay of leptoquarks (assumed spin 0).
        WT=1.
        WTMAX=1.

C...Obtain correct angular distribution by rejection techniques.
      ELSE
        WT=1.
        WTMAX=1.
      ENDIF
      IF(WT.LT.PYR(0)*WTMAX) GOTO 190

C...Construct massive four-vectors using angles chosen. Mark decayed
C...resonances, add documentation lines. Shower evolution.
  290 DO 340 JT=1,JTMAX
      IF(KDCY(JT).EQ.0) GOTO 340
      ID=IREF(IP,JT)
      IF(ISET(ISUB).NE.6) THEN
        CALL LUDBRB(NSD(JT)+1,NSD(JT)+2,ACOS(CTHE(JT)),PHI(JT),
     &  DBLE(P(ID,1)/P(ID,4)),DBLE(P(ID,2)/P(ID,4)),
     &  DBLE(P(ID,3)/P(ID,4)))
      ELSE
C...Z + q + q~ : angles already known, in rest frame of system.
        DO 300 J=1,3
  300   DBEZQQ(J)=(P(ID,J)+P(ID+1,J)+P(ID+2,J))/(P(ID,4)+P(ID+1,4)+
     &  P(ID+2,4))
        K(N+1,1)=1
        DO 310 J=1,5
  310   P(N+1,J)=P(ID,J)
        CALL LUDBRB(N+1,N+1,0.,0.,-DBEZQQ(1),-DBEZQQ(2),-DBEZQQ(3))
        PHIZQQ=ULANGL(P(N+1,1),P(N+1,2))
        THEZQQ=ULANGL(P(N+1,3),SQRT(P(N+1,1)**2+P(N+1,2)**2))
        CALL LUDBRB(NSD(JT)+1,NSD(JT)+2,ACOS(VINT(81)),VINT(82),
     &  0D0,0D0,DBLE(SQRT(P(N+1,1)**2+P(N+1,2)**2+P(N+1,3)**2)/
     &  P(N+1,4)))
        CALL LUDBRB(NSD(JT)+1,NSD(JT)+2,THEZQQ,PHIZQQ,DBEZQQ(1),
     &  DBEZQQ(2),DBEZQQ(3))
      ENDIF
      K(ID,1)=K(ID,1)+10
C...Do not kill colour flow through leptoquark!
      IF(IABS(K(ID,2)).NE.39) THEN
        K(ID,4)=NSD(JT)+1
        K(ID,5)=NSD(JT)+2
      ENDIF
      IDOC=MINT(83)+MINT(4)
      DO 330 I=NSD(JT)+1,NSD(JT)+2
      I1=MINT(83)+MINT(4)+1
      K(I,3)=I1
      IF(MSTP(128).GE.1) K(I,3)=ID
      IF(MSTP(128).LE.1) THEN
        MINT(4)=MINT(4)+1
        K(I1,1)=21
        K(I1,2)=K(I,2)
        K(I1,3)=IREF(IP,JT+2)
        DO 320 J=1,5
  320   P(I1,J)=P(I,J)
      ENDIF
  330 CONTINUE
      IF(MSTP(71).GE.1.AND.KDCY(JT).EQ.1) CALL LUSHOW(NSD(JT)+1,
     &NSD(JT)+2,P(ID,5))

C...Check if new resonances were produced.
      IF(KDCY(JT).NE.3) GOTO 340
      NP=NP+1
      IREF(NP,1)=NSD(JT)+1
      IREF(NP,2)=NSD(JT)+2
      IREF(NP,3)=IDOC+1
      IREF(NP,4)=IDOC+2
      IREF(NP,5)=K(IREF(IP,JT),2)
      IREF(NP,6)=IREF(IP,JT)
  340 CONTINUE

C...Fill information for 2 -> 1 -> 2. Loop back if needed.
      IF(JTMAX.EQ.1.AND.KDCY(1).NE.0) THEN
        MINT(7)=MINT(83)+6+2*ISET(ISUB)
        MINT(8)=MINT(83)+7+2*ISET(ISUB)
        MINT(25)=KFL1(1)
        MINT(26)=KFL2(1)
        VINT(23)=CTHE(1)
        RM3=P(N-1,5)**2/SH
        RM4=P(N,5)**2/SH
        BE34=SQRT(MAX(0.,(1.-RM3-RM4)**2-4.*RM3*RM4))
        VINT(45)=-0.5*SH*(1.-RM3-RM4-BE34*CTHE(1))
        VINT(46)=-0.5*SH*(1.-RM3-RM4+BE34*CTHE(1))
        VINT(48)=0.25*SH*BE34**2*MAX(0.,1.-CTHE(1)**2)
        VINT(47)=SQRT(VINT(48))
      ENDIF
  350 IF(IP.LT.NP) GOTO 100

      RETURN
      END

C*********************************************************************

      SUBROUTINE RYMULT(MMUL)
      PARAMETER (KSZJ=4000)

C...Initializes treatment of multiple interactions, selects kinematics
C...of hardest interaction if low-pT physics included in run, and
C...generates all non-hardest interactions.
      
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      integer N,NPAD,K
      double precision P,V
      SAVE /PYJETS/


      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/RYSUBS/MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200)
      COMMON/RYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/RYINT1/MINT(400),VINT(400)
      COMMON/RYINT2/ISET(200),KFPR(200,2),COEF(200,20),ICOL(40,4,2)
      COMMON/RYINT3/XSFX(2,-40:40),ISIG(1000,3),SIGH(1000)
      COMMON/RYINT5/NGEN(0:200,3),XSEC(0:200,3)
      SAVE /LUDAT1/,/LUDAT2/
      SAVE /RYSUBS/,/RYPARS/,/RYINT1/,/RYINT2/,/RYINT3/,/RYINT5/
      DIMENSION NMUL(20),SIGM(20),KSTR(500,2),VINTSV(80)
      SAVE XT2,XT2FAC,XC2,XTS,IRBIN,RBIN,NMUL,SIGM

C...Initialization of multiple interaction treatment.
      IF(MMUL.EQ.1) THEN
        IF(MSTP(122).GE.1) WRITE(MSTU(11),5000) MSTP(82)
        ISUB=96
        MINT(1)=96
        VINT(63)=0.
        VINT(64)=0.
        VINT(143)=1.
        VINT(144)=1.

C...Loop over phase space points: xT2 choice in 20 bins.
  100   SIGSUM=0.
        DO 120 IXT2=1,20
        NMUL(IXT2)=MSTP(83)
        SIGM(IXT2)=0.
        DO 110 ITRY=1,MSTP(83)
        RSCA=0.05*((21-IXT2)-PYR(0))
        XT2=VINT(149)*(1.+VINT(149))/(VINT(149)+RSCA)-VINT(149)
        XT2=MAX(0.01*VINT(149),XT2)
        VINT(25)=XT2

C...Choose tau and y*. Calculate cos(theta-hat).
        IF(PYR(0).LE.COEF(ISUB,1)) THEN
          TAUT=(2.*(1.+SQRT(1.-XT2))/XT2-1.)**PYR(0)
          TAU=XT2*(1.+TAUT)**2/(4.*TAUT)
        ELSE
          TAU=XT2*(1.+TAN(PYR(0)*ATAN(SQRT(1./XT2-1.)))**2)
        ENDIF
        VINT(21)=TAU
        CALL RYKLIM(2)
        RYST=PYR(0)
        MYST=1
        IF(RYST.GT.COEF(ISUB,8)) MYST=2
        IF(RYST.GT.COEF(ISUB,8)+COEF(ISUB,9)) MYST=3
        CALL RYKMAP(2,MYST,PYR(0))
        VINT(23)=SQRT(MAX(0.,1.-XT2/TAU))*(-1)**INT(1.5+PYR(0))

C...Calculate differential cross-section.
        VINT(71)=0.5*VINT(1)*SQRT(XT2)
        CALL RYSIGH(NCHN,SIGS)
  110   SIGM(IXT2)=SIGM(IXT2)+SIGS
  120   SIGSUM=SIGSUM+SIGM(IXT2)
        SIGSUM=SIGSUM/(20.*MSTP(83))

C...Reject result if sigma(parton-parton) is smaller than hadronic one.
        IF(SIGSUM.LT.1.1*VINT(106)) THEN
          IF(MSTP(122).GE.1) WRITE(MSTU(11),5100) PARP(82),SIGSUM
          PARP(82)=0.9*PARP(82)
          VINT(149)=4.*PARP(82)**2/VINT(2)
          GOTO 100
        ENDIF
        IF(MSTP(122).GE.1) WRITE(MSTU(11),5200) PARP(82), SIGSUM

C...Start iteration to find k factor.
        YKE=SIGSUM/VINT(106)
        SO=0.5
        XI=0.
        YI=0.
        XK=0.5
        IIT=0
  130   IF(IIT.EQ.0) THEN
          XK=2.*XK
        ELSEIF(IIT.EQ.1) THEN
          XK=0.5*XK
        ELSE
          XK=XI+(YKE-YI)*(XF-XI)/(YF-YI)
        ENDIF

C...Evaluate overlap integrals.
        IF(MSTP(82).EQ.2) THEN
          SP=0.5*PARU(1)*(1.-EXP(-XK))
          SOP=SP/PARU(1)
        ELSE
          IF(MSTP(82).EQ.3) DELTAB=0.02
          IF(MSTP(82).EQ.4) DELTAB=MIN(0.01,0.05*PARP(84))
          SP=0.
          SOP=0.
          B=-0.5*DELTAB
  140     B=B+DELTAB
          IF(MSTP(82).EQ.3) THEN
            OV=EXP(-B**2)/PARU(2)
          ELSE
            CQ2=PARP(84)**2
            OV=((1.-PARP(83))**2*EXP(-MIN(87.,B**2))+2.*PARP(83)*
     &      (1.-PARP(83))*2./(1.+CQ2)*EXP(-MIN(87.,B**2*2./(1.+CQ2)))+
     &      PARP(83)**2/CQ2*EXP(-MIN(87.,B**2/CQ2)))/PARU(2)
          ENDIF
          PACC=1.-EXP(-MIN(87.,PARU(1)*XK*OV))
          SP=SP+PARU(2)*B*DELTAB*PACC
          SOP=SOP+PARU(2)*B*DELTAB*OV*PACC
          IF(B.LT.1..OR.B*PACC.GT.1E-6) GOTO 140
        ENDIF
        YK=PARU(1)*XK*SO/SP

C...Continue iteration until convergence.
        IF(YK.LT.YKE) THEN
          XI=XK
          YI=YK
          IF(IIT.EQ.1) IIT=2
        ELSE
          XF=XK
          YF=YK
          IF(IIT.EQ.0) IIT=1
        ENDIF
        IF(ABS(YK-YKE).GE.1E-5*YKE) GOTO 130

C...Store some results for subsequent use.
        VINT(145)=SIGSUM
        VINT(146)=SOP/SO
        VINT(147)=SOP/SP

C...Initialize iteration in xT2 for hardest interaction.
      ELSEIF(MMUL.EQ.2) THEN
        IF(MSTP(82).LE.0) THEN
        ELSEIF(MSTP(82).EQ.1) THEN
          XT2=1.
          XT2FAC=XSEC(96,1)/VINT(106)*VINT(149)/(1.-VINT(149))
        ELSEIF(MSTP(82).EQ.2) THEN
          XT2=1.
          XT2FAC=VINT(146)*XSEC(96,1)/VINT(106)*VINT(149)*(1.+VINT(149))
        ELSE
          XC2=4.*CKIN(3)**2/VINT(2)
          IF(CKIN(3).LE.CKIN(5).OR.MINT(82).GE.2) XC2=0.
        ENDIF

      ELSEIF(MMUL.EQ.3) THEN
C...Low-pT or multiple interactions (first semihard interaction):
C...choose xT2 according to dpT2/pT2**2*exp(-(sigma above pT2)/norm)
C...or (MSTP(82)>=2) dpT2/(pT2+pT0**2)**2*exp(-....).
        ISUB=MINT(1)
        IF(MSTP(82).LE.0) THEN
          XT2=0.
        ELSEIF(MSTP(82).EQ.1) THEN
          XT2=XT2FAC*XT2/(XT2FAC-XT2*LOG(PYR(0)))
        ELSEIF(MSTP(82).EQ.2) THEN
          IF(XT2.LT.1..AND.EXP(-XT2FAC*XT2/(VINT(149)*(XT2+
     &    VINT(149)))).GT.PYR(0)) XT2=1.
          IF(XT2.GE.1.) THEN
            XT2=(1.+VINT(149))*XT2FAC/(XT2FAC-(1.+VINT(149))*LOG(1.-
     &      PYR(0)*(1.-EXP(-XT2FAC/(VINT(149)*(1.+VINT(149)))))))-
     &      VINT(149)
          ELSE
            XT2=-XT2FAC/LOG(EXP(-XT2FAC/(XT2+VINT(149)))+PYR(0)*
     &      (EXP(-XT2FAC/VINT(149))-EXP(-XT2FAC/(XT2+VINT(149)))))-
     &      VINT(149)
          ENDIF
          XT2=MAX(0.01*VINT(149),XT2)
        ELSE
          XT2=(XC2+VINT(149))*(1.+VINT(149))/(1.+VINT(149)-
     &    PYR(0)*(1.-XC2))-VINT(149)
          XT2=MAX(0.01*VINT(149),XT2)
        ENDIF
        VINT(25)=XT2

C...Low-pT: choose xT2, tau, y* and cos(theta-hat) fixed.
        IF(MSTP(82).LE.1.AND.XT2.LT.VINT(149)) THEN
          IF(MINT(82).EQ.1) NGEN(0,1)=NGEN(0,1)-1
          IF(MINT(82).EQ.1) NGEN(ISUB,1)=NGEN(ISUB,1)-1
          ISUB=95
          MINT(1)=ISUB
          VINT(21)=0.01*VINT(149)
          VINT(22)=0.
          VINT(23)=0.
          VINT(25)=0.01*VINT(149)

        ELSE
C...Multiple interactions (first semihard interaction).
C...Choose tau and y*. Calculate cos(theta-hat).
          IF(PYR(0).LE.COEF(ISUB,1)) THEN
            TAUT=(2.*(1.+SQRT(1.-XT2))/XT2-1.)**PYR(0)
            TAU=XT2*(1.+TAUT)**2/(4.*TAUT)
          ELSE
            TAU=XT2*(1.+TAN(PYR(0)*ATAN(SQRT(1./XT2-1.)))**2)
          ENDIF
          VINT(21)=TAU
          CALL RYKLIM(2)
          RYST=PYR(0)
          MYST=1
          IF(RYST.GT.COEF(ISUB,8)) MYST=2
          IF(RYST.GT.COEF(ISUB,8)+COEF(ISUB,9)) MYST=3
          CALL RYKMAP(2,MYST,PYR(0))
          VINT(23)=SQRT(MAX(0.,1.-XT2/TAU))*(-1)**INT(1.5+PYR(0))
        ENDIF
        VINT(71)=0.5*VINT(1)*SQRT(VINT(25))

C...Store results of cross-section calculation.
      ELSEIF(MMUL.EQ.4) THEN
        ISUB=MINT(1)
        XTS=VINT(25)
        IF(ISET(ISUB).EQ.1) XTS=VINT(21)
        IF(ISET(ISUB).EQ.2.OR.ISET(ISUB).EQ.6)
     &  XTS=(4.*VINT(48)+2.*VINT(63)+2.*VINT(64))/VINT(2)
        IF(ISET(ISUB).GE.3.AND.ISET(ISUB).LE.5) XTS=VINT(26)
        RBIN=MAX(0.000001,MIN(0.999999,XTS*(1.+VINT(149))/
     &  (XTS+VINT(149))))
        IRBIN=INT(1.+20.*RBIN)
        IF(ISUB.EQ.96) NMUL(IRBIN)=NMUL(IRBIN)+1
        IF(ISUB.EQ.96) SIGM(IRBIN)=SIGM(IRBIN)+VINT(153)

C...Choose impact parameter.
      ELSEIF(MMUL.EQ.5) THEN
        IF(MSTP(82).EQ.3) THEN
          VINT(148)=PYR(0)/(PARU(2)*VINT(147))
        ELSE
          RTYPE=PYR(0)
          CQ2=PARP(84)**2
          IF(RTYPE.LT.(1.-PARP(83))**2) THEN
            B2=-LOG(PYR(0))
          ELSEIF(RTYPE.LT.1.-PARP(83)**2) THEN
            B2=-0.5*(1.+CQ2)*LOG(PYR(0))
          ELSE
            B2=-CQ2*LOG(PYR(0))
          ENDIF
          VINT(148)=((1.-PARP(83))**2*EXP(-MIN(87.,B2))+2.*PARP(83)*
     &    (1.-PARP(83))*2./(1.+CQ2)*EXP(-MIN(87.,B2*2./(1.+CQ2)))+
     &    PARP(83)**2/CQ2*EXP(-MIN(87.,B2/CQ2)))/(PARU(2)*VINT(147))
        ENDIF

C...Multiple interactions (variable impact parameter) : reject with
C...probability exp(-overlap*cross-section above pT/normalization).
        RNCOR=(IRBIN-20.*RBIN)*NMUL(IRBIN)
        SIGCOR=(IRBIN-20.*RBIN)*SIGM(IRBIN)
        DO 150 IBIN=IRBIN+1,20
        RNCOR=RNCOR+NMUL(IBIN)
  150   SIGCOR=SIGCOR+SIGM(IBIN)
        SIGABV=(SIGCOR/RNCOR)*VINT(149)*(1.-XTS)/(XTS+VINT(149))
        VINT(150)=EXP(-MIN(87.,VINT(146)*VINT(148)*SIGABV/VINT(106)))

C...Generate additional multiple semihard interactions.
      ELSEIF(MMUL.EQ.6) THEN
        ISUBSV=MINT(1)
        DO 160 J=11,80
  160   VINTSV(J)=VINT(J)
        ISUB=96
        MINT(1)=96

C...Reconstruct strings in hard scattering.
        NMAX=MINT(84)+4
        IF(ISET(ISUBSV).EQ.1) NMAX=MINT(84)+2
        NSTR=0
        DO 180 I=MINT(84)+1,NMAX
        KCS=KCHG(LUCOMP(K(I,2)),2)*ISIGN(1,K(I,2))
        IF(KCS.EQ.0) GOTO 180
        DO 170 J=1,4
        IF(KCS.EQ.1.AND.(J.EQ.2.OR.J.EQ.4)) GOTO 170
        IF(KCS.EQ.-1.AND.(J.EQ.1.OR.J.EQ.3)) GOTO 170
        IF(J.LE.2) THEN
          IST=MOD(K(I,J+3)/MSTU(5),MSTU(5))
        ELSE
          IST=MOD(K(I,J+1),MSTU(5))
        ENDIF
        IF(IST.LT.MINT(84).OR.IST.GT.I) GOTO 170
        IF(KCHG(LUCOMP(K(IST,2)),2).EQ.0) GOTO 170
        NSTR=NSTR+1
        IF(J.EQ.1.OR.J.EQ.4) THEN
          KSTR(NSTR,1)=I
          KSTR(NSTR,2)=IST
        ELSE
          KSTR(NSTR,1)=IST
          KSTR(NSTR,2)=I
        ENDIF
  170   CONTINUE
  180   CONTINUE

C...Set up starting values for iteration in xT2.
        XT2=VINT(25)
        IF(ISET(ISUBSV).EQ.1) XT2=VINT(21)
        IF(ISET(ISUBSV).EQ.2) XT2=(4.*VINT(48)+2.*VINT(63)+2.*VINT(64))/
     &  VINT(2)
        IF(ISET(ISUBSV).EQ.3.OR.ISET(ISUBSV).EQ.4) XT2=VINT(26)
        IF(MSTP(82).LE.1) THEN
          XT2FAC=XSEC(ISUB,1)*VINT(149)/((1.-VINT(149))*VINT(106))
        ELSE
          XT2FAC=VINT(146)*VINT(148)*XSEC(ISUB,1)/VINT(106)*
     &    VINT(149)*(1.+VINT(149))
        ENDIF
        VINT(63)=0.
        VINT(64)=0.
        VINT(143)=1.-VINT(141)
        VINT(144)=1.-VINT(142)

C...Iterate downwards in xT2.
  190   IF(MSTP(82).LE.1) THEN
          XT2=XT2FAC*XT2/(XT2FAC-XT2*LOG(PYR(0)))
          IF(XT2.LT.VINT(149)) GOTO 230
        ELSE
          IF(XT2.LE.0.01*VINT(149)) GOTO 230
          XT2=XT2FAC*(XT2+VINT(149))/(XT2FAC-(XT2+VINT(149))*
     &    LOG(PYR(0)))-VINT(149)
          IF(XT2.LE.0.) GOTO 230
          XT2=MAX(0.01*VINT(149),XT2)
        ENDIF
        VINT(25)=XT2

C...Choose tau and y*. Calculate cos(theta-hat).
        IF(PYR(0).LE.COEF(ISUB,1)) THEN
          TAUT=(2.*(1.+SQRT(1.-XT2))/XT2-1.)**PYR(0)
          TAU=XT2*(1.+TAUT)**2/(4.*TAUT)
        ELSE
          TAU=XT2*(1.+TAN(PYR(0)*ATAN(SQRT(1./XT2-1.)))**2)
        ENDIF
        VINT(21)=TAU
        CALL RYKLIM(2)
        RYST=PYR(0)
        MYST=1
        IF(RYST.GT.COEF(ISUB,8)) MYST=2
        IF(RYST.GT.COEF(ISUB,8)+COEF(ISUB,9)) MYST=3
        CALL RYKMAP(2,MYST,PYR(0))
        VINT(23)=SQRT(MAX(0.,1.-XT2/TAU))*(-1)**INT(1.5+PYR(0))

C...Check that x not used up. Accept or reject kinematical variables.
        X1M=SQRT(TAU)*EXP(VINT(22))
        X2M=SQRT(TAU)*EXP(-VINT(22))
        IF(VINT(143)-X1M.LT.0.01.OR.VINT(144)-X2M.LT.0.01) GOTO 190
        VINT(71)=0.5*VINT(1)*SQRT(XT2)
        CALL RYSIGH(NCHN,SIGS)
        IF(SIGS.LT.XSEC(ISUB,1)*PYR(0)) GOTO 190

C...Reset K, P and V vectors. Select some variables.
        DO 200 I=N+1,N+2
        DO 200 J=1,5
        K(I,J)=0
        P(I,J)=0.
  200   V(I,J)=0.
        RFLAV=PYR(0)
        PT=0.5*VINT(1)*SQRT(XT2)
        PHI=PARU(2)*PYR(0)
        CTH=VINT(23)

C...Add first parton to event record.
        K(N+1,1)=3
        K(N+1,2)=21
        IF(RFLAV.GE.MAX(PARP(85),PARP(86))) K(N+1,2)=
     &  1+INT((2.+PARJ(2))*PYR(0))
        P(N+1,1)=PT*COS(PHI)
        P(N+1,2)=PT*SIN(PHI)
        P(N+1,3)=0.25*VINT(1)*(VINT(41)*(1.+CTH)-VINT(42)*(1.-CTH))
        P(N+1,4)=0.25*VINT(1)*(VINT(41)*(1.+CTH)+VINT(42)*(1.-CTH))
        P(N+1,5)=0.

C...Add second parton to event record.
        K(N+2,1)=3
        K(N+2,2)=21
        IF(K(N+1,2).NE.21) K(N+2,2)=-K(N+1,2)
        P(N+2,1)=-P(N+1,1)
        P(N+2,2)=-P(N+1,2)
        P(N+2,3)=0.25*VINT(1)*(VINT(41)*(1.-CTH)-VINT(42)*(1.+CTH))
        P(N+2,4)=0.25*VINT(1)*(VINT(41)*(1.-CTH)+VINT(42)*(1.+CTH))
        P(N+2,5)=0.

        IF(RFLAV.LT.PARP(85).AND.NSTR.GE.1) THEN
C....Choose relevant string pieces to place gluons on.
          DO 220 I=N+1,N+2
          DMIN=1E8
          DO 210 ISTR=1,NSTR
          I1=KSTR(ISTR,1)
          I2=KSTR(ISTR,2)
          DIST=(P(I,4)*P(I1,4)-P(I,1)*P(I1,1)-P(I,2)*P(I1,2)-
     &    P(I,3)*P(I1,3))*(P(I,4)*P(I2,4)-P(I,1)*P(I2,1)-
     &    P(I,2)*P(I2,2)-P(I,3)*P(I2,3))/MAX(1.,P(I1,4)*P(I2,4)-
     &    P(I1,1)*P(I2,1)-P(I1,2)*P(I2,2)-P(I1,3)*P(I2,3))
          IF(ISTR.EQ.1.OR.DIST.LT.DMIN) THEN
            DMIN=DIST
            IST1=I1
            IST2=I2
            ISTM=ISTR
          ENDIF
  210     CONTINUE

C....Colour flow adjustments, new string pieces.
          IF(K(IST1,4)/MSTU(5).EQ.IST2) K(IST1,4)=MSTU(5)*I+
     &    MOD(K(IST1,4),MSTU(5))
          IF(MOD(K(IST1,5),MSTU(5)).EQ.IST2) K(IST1,5)=
     &    MSTU(5)*(K(IST1,5)/MSTU(5))+I
          K(I,5)=MSTU(5)*IST1
          K(I,4)=MSTU(5)*IST2
          IF(K(IST2,5)/MSTU(5).EQ.IST1) K(IST2,5)=MSTU(5)*I+
     &    MOD(K(IST2,5),MSTU(5))
          IF(MOD(K(IST2,4),MSTU(5)).EQ.IST1) K(IST2,4)=
     &    MSTU(5)*(K(IST2,4)/MSTU(5))+I
          KSTR(ISTM,2)=I
          KSTR(NSTR+1,1)=I
          KSTR(NSTR+1,2)=IST2
  220     NSTR=NSTR+1

C...String drawing and colour flow for gluon loop.
        ELSEIF(K(N+1,2).EQ.21) THEN
          K(N+1,4)=MSTU(5)*(N+2)
          K(N+1,5)=MSTU(5)*(N+2)
          K(N+2,4)=MSTU(5)*(N+1)
          K(N+2,5)=MSTU(5)*(N+1)
          KSTR(NSTR+1,1)=N+1
          KSTR(NSTR+1,2)=N+2
          KSTR(NSTR+2,1)=N+2
          KSTR(NSTR+2,2)=N+1
          NSTR=NSTR+2

C...String drawing and colour flow for qq~ pair.
        ELSE
          K(N+1,4)=MSTU(5)*(N+2)
          K(N+2,5)=MSTU(5)*(N+1)
          KSTR(NSTR+1,1)=N+1
          KSTR(NSTR+1,2)=N+2
          NSTR=NSTR+1
        ENDIF

C...Update remaining energy; iterate.
        N=N+2
        IF(N.GT.MSTU(4)-MSTU(32)-10) THEN
          CALL LUERRM(11,'(RYMULT:) no more memory left in LUJETS')
          IF(MSTU(21).GE.1) RETURN
        ENDIF
        MINT(31)=MINT(31)+1
        VINT(151)=VINT(151)+VINT(41)
        VINT(152)=VINT(152)+VINT(42)
        VINT(143)=VINT(143)-VINT(41)
        VINT(144)=VINT(144)-VINT(42)
        IF(MINT(31).LT.240) GOTO 190
  230   CONTINUE
        MINT(1)=ISUBSV
        DO 240 J=11,80
  240   VINT(J)=VINTSV(J)
      ENDIF

C...Format statements for printout.
 5000 FORMAT(/1X,'****** RYMULT: initialization of multiple inter',
     &'actions for MSTP(82) =',I2,' ******')
 5100 FORMAT(8X,'pT0 =',F5.2,' GeV gives sigma(parton-parton) =',1P,
     &E9.2,' mb: rejected')
 5200 FORMAT(8X,'pT0 =',F5.2,' GeV gives sigma(parton-parton) =',1P,
     &E9.2,' mb: accepted')

      RETURN
      END

C*********************************************************************

      SUBROUTINE RYREMN(IPU1,IPU2)
      PARAMETER (KSZJ=4000)

C...Adds on target remnants (one or two from each side) and
C...includes primordial kT for hadron beams.
      
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      integer N,NPAD,K
      double precision P,V
      SAVE /PYJETS/


      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/RYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/RYINT1/MINT(400),VINT(400)
      SAVE /LUDAT1/,/LUDAT2/
      SAVE /RYPARS/,/RYINT1/
      DIMENSION KFLCH(2),KFLSP(2),CHI(2),PMS(0:6),IS(2),ISN(2),ROBO(5),
     &PSYS(0:2,5),PMIN(0:2)

C...Find event type and remaining energy.
      ISUB=MINT(1)
      NS=N
      IF(MINT(44).NE.4.OR.MSTP(81).LE.0) THEN
        VINT(143)=1.-VINT(141)
        VINT(144)=1.-VINT(142)
      ENDIF

C...Define initial partons.
  100 DO 130 JT=1,2
      I=MINT(83)+JT+2
      IF(JT.EQ.1) IPU=IPU1
      IF(JT.EQ.2) IPU=IPU2
      K(I,1)=21
      K(I,2)=K(IPU,2)
      K(I,3)=I-2
      PMS(JT)=0.
      IF(MINT(47).EQ.1) THEN
        DO 110 J=1,5
  110   P(I,J)=P(I-2,J)
      ELSEIF(ISUB.EQ.95) THEN
        K(I,2)=21
      ELSE
        P(I,5)=P(IPU,5)

C...No primordial kT, or chosen according to truncated Gaussian or
C...exponential.
  120   IF(MINT(40+JT).EQ.1.OR.MSTP(91).LE.0) THEN
          PT=0.
        ELSEIF(MSTP(91).EQ.1) THEN
          PT=PARP(91)*SQRT(-LOG(PYR(0)))
        ELSE
          RPT1=PYR(0)
          RPT2=PYR(0)
          PT=-PARP(92)*LOG(RPT1*RPT2)
        ENDIF
        IF(PT.GT.PARP(93)) GOTO 120
        PHI=PARU(2)*PYR(0)
        P(I,1)=PT*COS(PHI)
        P(I,2)=PT*SIN(PHI)
        PMS(JT)=P(I,5)**2+P(I,1)**2+P(I,2)**2
      ENDIF
  130 CONTINUE
      IF(MINT(47).EQ.1) RETURN

C...Kinematics construction for initial partons.
      I1=MINT(83)+3
      I2=MINT(83)+4
      IF(ISUB.EQ.95) THEN
        SHS=0.
        SHR=0.
      ELSE
        SHS=VINT(141)*VINT(142)*VINT(2)+(P(I1,1)+P(I2,1))**2+
     &  (P(I1,2)+P(I2,2))**2
        SHR=SQRT(MAX(0.,SHS))
        IF((SHS-PMS(1)-PMS(2))**2-4.*PMS(1)*PMS(2).LE.0.) GOTO 100
        P(I1,4)=0.5*(SHR+(PMS(1)-PMS(2))/SHR)
        P(I1,3)=SQRT(MAX(0.,P(I1,4)**2-PMS(1)))
        P(I2,4)=SHR-P(I1,4)
        P(I2,3)=-P(I1,3)

C...Transform partons to overall CM-frame.
        ROBO(3)=(P(I1,1)+P(I2,1))/SHR
        ROBO(4)=(P(I1,2)+P(I2,2))/SHR
        CALL LUDBRB(I1,I2,0.,0.,-DBLE(ROBO(3)),-DBLE(ROBO(4)),0D0)
        ROBO(2)=ULANGL(P(I1,1),P(I1,2))
        CALL LUDBRB(I1,I2,0.,-ROBO(2),0D0,0D0,0D0)
        ROBO(1)=ULANGL(P(I1,3),P(I1,1))
        CALL LUDBRB(I1,I2,-ROBO(1),0.,0D0,0D0,0D0)
        CALL LUDBRB(I1,MINT(52),ROBO(1),ROBO(2),DBLE(ROBO(3)),
     &  DBLE(ROBO(4)),0D0)
        ROBO(5)=MAX(-0.999999,MIN(0.999999,(VINT(141)-VINT(142))/
     &  (VINT(141)+VINT(142))))
        CALL LUDBRB(I1,MINT(52),0.,0.,0D0,0D0,DBLE(ROBO(5)))
      ENDIF

C...Check minimum invariant mass of remnant system(s).
      PSYS(0,4)=P(I1,4)+P(I2,4)+0.5*VINT(1)*(VINT(151)+VINT(152))
      PSYS(0,3)=P(I1,3)+P(I2,3)+0.5*VINT(1)*(VINT(151)-VINT(152))
      PMS(0)=MAX(0.,PSYS(0,4)**2-PSYS(0,3)**2)
      PMIN(0)=SQRT(PMS(0))
      DO 140 JT=1,2
      PSYS(JT,4)=0.5*VINT(1)*VINT(142+JT)
      PMIN(JT)=0.
      IF(MINT(40+JT).EQ.1.AND.MSTP(11).EQ.0) GOTO 140
      CALL RYSPLI(MINT(10+JT),MINT(12+JT),KFLCH(JT),KFLSP(JT))
      IF(KFLCH(JT).NE.0) PMIN(JT)=PMIN(JT)+ULMASS(KFLCH(JT))
      IF(KFLSP(JT).NE.0) PMIN(JT)=PMIN(JT)+ULMASS(KFLSP(JT))
      IF(KFLCH(JT)*KFLSP(JT).NE.0) PMIN(JT)=PMIN(JT)+0.5*PARP(111)
      PMIN(JT)=SQRT(PMIN(JT)**2+P(MINT(83)+JT+2,1)**2+
     &P(MINT(83)+JT+2,2)**2)
  140 CONTINUE
      IF(PMIN(0)+PMIN(1)+PMIN(2).GT.VINT(1).OR.PMIN(1).GT.PSYS(1,4).OR.
     &PMIN(2).GT.PSYS(2,4)) THEN
        MINT(51)=1
        RETURN
      ENDIF

C...Loop over two remnants; skip if none there.
      I=NS
      DO 210 JT=1,2
      ISN(JT)=0
      IF(MINT(40+JT).EQ.1.AND.MSTP(11).EQ.0) GOTO 210
      IF(JT.EQ.1) IPU=IPU1
      IF(JT.EQ.2) IPU=IPU2

C...Store first remnant parton.
      I=I+1
      IS(JT)=I
      ISN(JT)=1
      DO 150 J=1,5
      K(I,J)=0
      P(I,J)=0.
  150 V(I,J)=0.
      K(I,1)=1
      K(I,2)=KFLSP(JT)
      K(I,3)=MINT(83)+JT
      P(I,5)=ULMASS(K(I,2))

C...First parton colour connections and kinematics.
      KCOL=KCHG(LUCOMP(KFLSP(JT)),2)
      IF(KCOL.NE.0) THEN
        K(I,1)=3
        KFLS=(3-KCOL*ISIGN(1,KFLSP(JT)))/2
        K(I,KFLS+3)=IPU
        K(IPU,6-KFLS)=MOD(K(IPU,6-KFLS),MSTU(5))+MSTU(5)*I
      ENDIF
      IF(KFLCH(JT).EQ.0) THEN
        P(I,1)=-P(MINT(83)+JT+2,1)
        P(I,2)=-P(MINT(83)+JT+2,2)
        PMS(JT)=P(I,5)**2+P(I,1)**2+P(I,2)**2
        PSYS(JT,3)=SQRT(MAX(0.,PSYS(JT,4)**2-PMS(JT)))*(-1)**(JT-1)
        P(I,3)=PSYS(JT,3)
        P(I,4)=PSYS(JT,4)

C...When extra remnant parton or hadron: store extra remnant.
      ELSE
        I=I+1
        ISN(JT)=2
        DO 160 J=1,5
        K(I,J)=0
        P(I,J)=0.
  160   V(I,J)=0.
        K(I,1)=1
        K(I,2)=KFLCH(JT)
        K(I,3)=MINT(83)+JT
        P(I,5)=ULMASS(K(I,2))

C...Find parton colour connections of extra remnant.
        KCOL=KCHG(LUCOMP(KFLCH(JT)),2)
        IF(KCOL.EQ.2) THEN
          K(I,1)=3
          K(I,4)=MSTU(5)*IPU+IPU
          K(I,5)=MSTU(5)*IPU+IPU
          K(IPU,4)=MOD(K(IPU,4),MSTU(5))+MSTU(5)*I
          K(IPU,5)=MOD(K(IPU,5),MSTU(5))+MSTU(5)*I
        ELSEIF(KCOL.NE.0) THEN
          K(I,1)=3
          KFLS=(3-KCOL*ISIGN(1,KFLCH(JT)))/2
          K(I,KFLS+3)=IPU
          K(IPU,6-KFLS)=MOD(K(IPU,6-KFLS),MSTU(5))+MSTU(5)*I
        ENDIF

C...Relative transverse momentum when two remnants.
  170   CALL LUPTDI(1,P(I-1,1),P(I-1,2))
        IF(IABS(MINT(10+JT)).LT.20) THEN
          P(I-1,1)=0.
          P(I-1,2)=0.
        ENDIF
        PMS(JT+2)=P(I-1,5)**2+P(I-1,1)**2+P(I-1,2)**2
        P(I,1)=-P(MINT(83)+JT+2,1)-P(I-1,1)
        P(I,2)=-P(MINT(83)+JT+2,2)-P(I-1,2)
        PMS(JT+4)=P(I,5)**2+P(I,1)**2+P(I,2)**2

C...Relative distribution of electron energy into electron plus parton.
        IMB=1
        IF(MOD(MINT(10+JT)/1000,10).NE.0) IMB=2
        IF(IABS(MINT(10+JT)).LT.20) THEN
          XHRD=VINT(140+JT)
          T=LOG(MIN(1E4,MAX(1.,VINT(52)))/0.16)
          NFE=1
          IF(VINT(52).GT.25.) NFE=2
          IF(VINT(52).GT.300.) NFE=3
          CALL RYSTGA(NFE,XHRD,T,XPGL1,XPQU1,XPQD1)
          CALL RYSTGA(NFE,0.999999,T,XPGL2,XPQU2,XPQD2)
          IF(KFLCH(JT).EQ.21) THEN
            WTMX=2.*MAX(XPGL1,XPGL2)
          ELSEIF(MOD(IABS(KFLCH(JT)),2).EQ.0) THEN
            WTMX=2.*MAX(XPQU1,XPQU2)
          ELSE
            WTMX=2.*MAX(XPQD1,XPQD2)
          ENDIF
  180     XE=XHRD**PYR(0)
          XG=MIN(0.999999,XHRD/XE)
          XPGA=1.+(1.-XE)**2
          CALL RYSTGA(NFE,XG,T,XPGL,XPQU,XPQD)
          IF(KFLCH(JT).EQ.21) THEN
            WT=XPGA*XPGL
          ELSEIF(MOD(IABS(KFLCH(JT)),2).EQ.0) THEN
            WT=XPGA*XPQU
          ELSE
            WT=XPGA*XPQD
          ENDIF
          IF(WT.LT.PYR(0)*WTMX) GOTO 180
          CHI(JT)=(XE-XHRD)/(1.-XHRD)

C...Relative distribution of energy for particle into two jets.
        ELSEIF(IABS(KFLCH(JT)).LE.10.OR.KFLCH(JT).EQ.21) THEN
          CHIK=PARP(92+2*IMB)
          IF(MSTP(92).LE.1) THEN
            IF(IMB.EQ.1) CHI(JT)=PYR(0)
            IF(IMB.EQ.2) CHI(JT)=1.-SQRT(PYR(0))
          ELSEIF(MSTP(92).EQ.2) THEN
            CHI(JT)=1.-PYR(0)**(1./(1.+CHIK))
          ELSEIF(MSTP(92).EQ.3) THEN
            CUT=2.*0.3/VINT(1)
  190       CHI(JT)=PYR(0)**2
            IF((CHI(JT)**2/(CHI(JT)**2+CUT**2))**0.25*(1.-CHI(JT))**CHIK
     &      .LT.PYR(0)) GOTO 190
          ELSE
            CUT=2.*0.3/VINT(1)
            CUTR=(1.+SQRT(1.+CUT**2))/CUT
  200       CHIR=CUT*CUTR**PYR(0)
            CHI(JT)=(CHIR**2-CUT**2)/(2.*CHIR)
            IF((1.-CHI(JT))**CHIK.LT.PYR(0)) GOTO 200
          ENDIF

C...Relative distribution of energy for particle into jet plus particle.
        ELSE
          IF(MSTP(92).LE.1) THEN
            IF(IMB.EQ.1) CHI(JT)=PYR(0)
            IF(IMB.EQ.2) CHI(JT)=1.-SQRT(PYR(0))
          ELSE
            CHI(JT)=1.-PYR(0)**(1./(1.+PARP(93+2*IMB)))
          ENDIF
          IF(MOD(KFLCH(JT)/1000,10).NE.0) CHI(JT)=1.-CHI(JT)
        ENDIF

C...Construct total transverse mass; reject if too large.
        PMS(JT)=PMS(JT+4)/CHI(JT)+PMS(JT+2)/(1.-CHI(JT))
        IF(PMS(JT).GT.PSYS(JT,4)**2) GOTO 170
        PSYS(JT,3)=SQRT(PSYS(JT,4)**2-PMS(JT))*(-1)**(JT-1)

C...Subdivide longitudinal momentum according to value selected above.
        PW1=CHI(JT)*(PSYS(JT,4)+ABS(PSYS(JT,3)))
        P(IS(JT)+1,4)=0.5*(PW1+PMS(JT+4)/PW1)
        P(IS(JT)+1,3)=0.5*(PW1-PMS(JT+4)/PW1)*(-1)**(JT-1)
        P(IS(JT),4)=PSYS(JT,4)-P(IS(JT)+1,4)
        P(IS(JT),3)=PSYS(JT,3)-P(IS(JT)+1,3)
      ENDIF
  210 CONTINUE
      N=I

C...Check if longitudinal boosts needed - if so pick two systems.
      PDEV=ABS(PSYS(0,4)+PSYS(1,4)+PSYS(2,4)-VINT(1))+
     &ABS(PSYS(0,3)+PSYS(1,3)+PSYS(2,3))
      IF(PDEV.LE.1E-6*VINT(1)) RETURN
      IF(VINT(143).GT.0.2.AND.VINT(144).GT.0.2) THEN
        IR=1
        IL=2
      ELSEIF(VINT(143).GT.0.2) THEN
        IR=1
        IL=0
      ELSEIF(VINT(144).GT.0.2) THEN
        IR=0
        IL=2
      ELSEIF(PMS(1)/PSYS(1,4)**2.GT.PMS(2)/PSYS(2,4)**2) THEN
        IR=1
        IL=0
      ELSE
        IR=0
        IL=2
      ENDIF
      IG=3-IR-IL

C...Construct longitudinal boosts.
      IF((IG.EQ.1.AND.ISN(1).EQ.0).OR.(IG.EQ.2.AND.ISN(2).EQ.0)) THEN
        PPB=VINT(1)
        PNB=VINT(1)
      ELSE
        PPB=VINT(1)-PSYS(IG,4)-PSYS(IG,3)
        PNB=VINT(1)-PSYS(IG,4)+PSYS(IG,3)
      ENDIF
      PMTB=PPB*PNB
      PMTR=PMS(IR)
      PMTL=PMS(IL)
      SQLAM=SQRT(MAX(0.,(PMTB-PMTR-PMTL)**2-4.*PMTR*PMTL))
      SQSGN=SIGN(1.,PSYS(IR,3)*PSYS(IL,4)-PSYS(IL,3)*PSYS(IR,4))
      RKR=(PMTB+PMTR-PMTL+SQLAM*SQSGN)/(2.*(PSYS(IR,4)+PSYS(IR,3))*PNB)
      RKL=(PMTB+PMTL-PMTR+SQLAM*SQSGN)/(2.*(PSYS(IL,4)-PSYS(IL,3))*PPB)
      BER=(RKR**2-1.)/(RKR**2+1.)
      BEL=-(RKL**2-1.)/(RKL**2+1.)

C...Perform longitudinal boosts.
      IF(IR.EQ.1) THEN
        CALL LUDBRB(IS(1),IS(1)+ISN(1)-1,0.,0.,0D0,0D0,DBLE(BER))
      ELSE
        CALL LUDBRB(I1,NS,0.,0.,0D0,0D0,DBLE(BER))
      ENDIF
      IF(IL.EQ.2) THEN
        CALL LUDBRB(IS(2),IS(2)+ISN(2)-1,0.,0.,0D0,0D0,DBLE(BEL))
      ELSE
        CALL LUDBRB(I1,NS,0.,0.,0D0,0D0,DBLE(BEL))
      ENDIF

      RETURN
      END

C*********************************************************************

      SUBROUTINE RYDIFF
      PARAMETER (KSZJ=4000)

C...Handles diffractive and elastic scattering.
      
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      integer N,NPAD,K
      double precision P,V
      SAVE /PYJETS/


      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/RYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/RYINT1/MINT(400),VINT(400)
      SAVE /LUDAT1/
      SAVE /RYPARS/,/RYINT1/

C...Reset K, P and V vectors. Store incoming particles.
      DO 100 JT=1,MSTP(126)+10
      I=MINT(83)+JT
      DO 100 J=1,5
      K(I,J)=0
      P(I,J)=0.
  100 V(I,J)=0.
      N=MINT(84)
      MINT(3)=0
      MINT(21)=0
      MINT(22)=0
      MINT(23)=0
      MINT(24)=0
      MINT(4)=4
      DO 110 JT=1,2
      I=MINT(83)+JT
      K(I,1)=21
      K(I,2)=MINT(10+JT)
      P(I,5)=VINT(2+JT)
      P(I,3)=VINT(5)*(-1)**(JT+1)
  110 P(I,4)=SQRT(P(I,3)**2+P(I,5)**2)
      MINT(6)=2

C...Subprocess; kinematics.
      SQLAM=(VINT(2)-VINT(63)-VINT(64))**2-4.*VINT(63)*VINT(64)
      PZ=SQRT(SQLAM)/(2.*VINT(1))
      DO 150 JT=1,2
      I=MINT(83)+JT
      PE=(VINT(2)+VINT(62+JT)-VINT(65-JT))/(2.*VINT(1))

C...Elastically scattered particle.
      IF(MINT(16+JT).LE.0) THEN
        N=N+1
        K(N,1)=1
        K(N,2)=K(I,2)
        K(N,3)=I+2
        P(N,3)=PZ*(-1)**(JT+1)
        P(N,4)=PE
        P(N,5)=P(I,5)

C...Diffracted particle: valence quark kicked out.
      ELSEIF(MSTP(101).EQ.1) THEN
        N=N+2
        K(N-1,1)=2
        K(N,1)=1
        K(N-1,3)=I+2
        K(N,3)=I+2
        CALL RYSPLI(K(I,2),21,K(N,2),K(N-1,2))
        P(N-1,5)=ULMASS(K(N-1,2))
        P(N,5)=ULMASS(K(N,2))
        SQLAM=(VINT(62+JT)-P(N-1,5)**2-P(N,5)**2)**2-
     &  4.*P(N-1,5)**2*P(N,5)**2
        P(N-1,3)=(PE*SQRT(SQLAM)+PZ*(VINT(62+JT)+P(N-1,5)**2-
     &  P(N,5)**2))/(2.*VINT(62+JT))*(-1)**(JT+1)
        P(N-1,4)=SQRT(P(N-1,3)**2+P(N-1,5)**2)
        P(N,3)=PZ*(-1)**(JT+1)-P(N-1,3)
        P(N,4)=SQRT(P(N,3)**2+P(N,5)**2)

C...Diffracted particle: gluon kicked out.
      ELSE
        N=N+3
        K(N-2,1)=2
        K(N-1,1)=2
        K(N,1)=1
        K(N-2,3)=I+2
        K(N-1,3)=I+2
        K(N,3)=I+2
        CALL RYSPLI(K(I,2),21,K(N,2),K(N-2,2))
        K(N-1,2)=21
        P(N-2,5)=ULMASS(K(N-2,2))
        P(N-1,5)=0.
        P(N,5)=ULMASS(K(N,2))
C...Energy distribution for particle into two jets.
  120   IMB=1
        IF(MOD(K(I,2)/1000,10).NE.0) IMB=2
        CHIK=PARP(92+2*IMB)
        IF(MSTP(92).LE.1) THEN
          IF(IMB.EQ.1) CHI=PYR(0)
          IF(IMB.EQ.2) CHI=1.-SQRT(PYR(0))
        ELSEIF(MSTP(92).EQ.2) THEN
          CHI=1.-PYR(0)**(1./(1.+CHIK))
        ELSEIF(MSTP(92).EQ.3) THEN
          CUT=2.*0.3/VINT(1)
  130     CHI=PYR(0)**2
          IF((CHI**2/(CHI**2+CUT**2))**0.25*(1.-CHI)**CHIK.LT.
     &    PYR(0)) GOTO 130
        ELSE
          CUT=2.*0.3/VINT(1)
          CUTR=(1.+SQRT(1.+CUT**2))/CUT
  140     CHIR=CUT*CUTR**PYR(0)
          CHI=(CHIR**2-CUT**2)/(2.*CHIR)
          IF((1.-CHI)**CHIK.LT.PYR(0)) GOTO 140
        ENDIF
        IF(CHI.LT.P(N,5)**2/VINT(62+JT).OR.CHI.GT.1.-P(N-2,5)**2/
     &  VINT(62+JT)) GOTO 120
        SQM=P(N-2,5)**2/(1.-CHI)+P(N,5)**2/CHI
        IF((SQRT(SQM)+PARJ(32))**2.GE.VINT(62+JT)) GOTO 120
        PZI=(PE*(VINT(62+JT)-SQM)+PZ*(VINT(62+JT)+SQM))/
     &  (2.*VINT(62+JT))
        PEI=SQRT(PZI**2+SQM)
        PQQP=(1.-CHI)*(PEI+PZI)
        P(N-2,3)=0.5*(PQQP-P(N-2,5)**2/PQQP)*(-1)**(JT+1)
        P(N-2,4)=SQRT(P(N-2,3)**2+P(N-2,5)**2)
        P(N-1,3)=(PZ-PZI)*(-1)**(JT+1)
        P(N-1,4)=ABS(P(N-1,3))
        P(N,3)=PZI*(-1)**(JT+1)-P(N-2,3)
        P(N,4)=SQRT(P(N,3)**2+P(N,5)**2)
      ENDIF

C...Documentation lines.
      K(I+2,1)=21
      IF(MINT(16+JT).EQ.0) K(I+2,2)=MINT(10+JT)
      IF(MINT(16+JT).NE.0) K(I+2,2)=10*(MINT(10+JT)/10)
      K(I+2,3)=I
      P(I+2,3)=PZ*(-1)**(JT+1)
      P(I+2,4)=PE
      P(I+2,5)=SQRT(VINT(62+JT))
  150 CONTINUE

C...Rotate outgoing partons/particles using cos(theta).
      CALL LUDBRB(MINT(83)+3,N,ACOS(VINT(23)),VINT(24),0D0,0D0,0D0)

      RETURN
      END

C*********************************************************************

      SUBROUTINE RYDOCU
      PARAMETER (KSZJ=4000)

C...Handles the decumentation of the process in MSTI and PARI,
C...and also computes cross-sections based on accumulated statistics.
      
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      integer N,NPAD,K
      double precision P,V
      SAVE /PYJETS/


      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/RYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/RYSUBS/MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200)
      COMMON/RYINT1/MINT(400),VINT(400)
      COMMON/RYINT2/ISET(200),KFPR(200,2),COEF(200,20),ICOL(40,4,2)
      COMMON/RYINT5/NGEN(0:200,3),XSEC(0:200,3)
      SAVE /LUDAT1/
      SAVE /RYSUBS/,/RYPARS/,/RYINT1/,/RYINT2/,/RYINT5/

C...Calculate Monte Carlo estimates of cross-sections.
      ISUB=MINT(1)
      IF(MSTP(111).NE.-1) NGEN(ISUB,3)=NGEN(ISUB,3)+1
      NGEN(0,3)=NGEN(0,3)+1
      XSEC(0,3)=0.
      DO 100 I=1,200
      IF(I.EQ.96) THEN
        XSEC(I,3)=0.
      ELSEIF(MSUB(95).EQ.1.AND.(I.EQ.11.OR.I.EQ.12.OR.I.EQ.13.OR.
     &I.EQ.28.OR.I.EQ.53.OR.I.EQ.68)) THEN
        XSEC(I,3)=XSEC(96,2)*NGEN(I,3)/MAX(1.,FLOAT(NGEN(96,1))*
     &  FLOAT(NGEN(96,2)))
      ELSEIF(NGEN(I,1).EQ.0) THEN
        XSEC(I,3)=0.
      ELSEIF(NGEN(I,2).EQ.0) THEN
        XSEC(I,3)=XSEC(I,2)*NGEN(0,3)/(FLOAT(NGEN(I,1))*
     &  FLOAT(NGEN(0,2)))
      ELSE
        XSEC(I,3)=XSEC(I,2)*NGEN(I,3)/(FLOAT(NGEN(I,1))*
     &  FLOAT(NGEN(I,2)))
      ENDIF
  100 XSEC(0,3)=XSEC(0,3)+XSEC(I,3)

C...Rescale to known low-pT cross-section for standard QCD processes.
      IF(MSUB(95).EQ.1) THEN
        NGENS=NGEN(91,3)+NGEN(92,3)+NGEN(93,3)+NGEN(94,3)+NGEN(95,3)
        XSECS=XSEC(91,3)+XSEC(92,3)+XSEC(93,3)+XSEC(94,3)+XSEC(95,3)
        XMAXS=XSEC(95,1)
        IF(MSUB(91).EQ.1) XMAXS=XMAXS+XSEC(91,1)
        IF(MSUB(92).EQ.1) XMAXS=XMAXS+XSEC(92,1)
        IF(MSUB(93).EQ.1) XMAXS=XMAXS+XSEC(93,1)
        IF(MSUB(94).EQ.1) XMAXS=XMAXS+XSEC(94,1)
        FAC=1.
        IF(NGENS.LT.NGEN(0,3)) FAC=(XMAXS-XSECS)/(XSEC(0,3)-XSECS)
        XSEC(11,3)=FAC*XSEC(11,3)
        XSEC(12,3)=FAC*XSEC(12,3)
        XSEC(13,3)=FAC*XSEC(13,3)
        XSEC(28,3)=FAC*XSEC(28,3)
        XSEC(53,3)=FAC*XSEC(53,3)
        XSEC(68,3)=FAC*XSEC(68,3)
        XSEC(0,3)=XSEC(91,3)+XSEC(92,3)+XSEC(93,3)+XSEC(94,3)+
     &  XSEC(95,1)
      ENDIF

C...Reset information on hard interaction.
      DO 110 J=1,200
      MSTI(J)=0
  110 PARI(J)=0.

C...Copy integer valued information from MINT into MSTI.
      DO 120 J=1,31
  120 MSTI(J)=MINT(J)

C...Store cross-section and kinematics variables in PARI.
      PARI(1)=XSEC(0,3)
      PARI(2)=XSEC(0,3)/MINT(5)
      PARI(11)=VINT(1)
      PARI(12)=VINT(2)
      IF(ISUB.NE.95) THEN
        DO 130 J=13,22
  130   PARI(J)=VINT(30+J)
        PARI(31)=VINT(141)
        PARI(32)=VINT(142)
        PARI(33)=VINT(41)
        PARI(34)=VINT(42)
        PARI(35)=PARI(33)-PARI(34)
        PARI(36)=VINT(21)
        PARI(37)=VINT(22)
        PARI(38)=VINT(26)
        PARI(41)=VINT(23)
        PARI(42)=2.*VINT(47)/VINT(1)
      ENDIF

C...Store information on scattered partons in PARI.
      IF(ISUB.NE.95.AND.MINT(7)*MINT(8).NE.0) THEN
        DO 140 IS=7,8
        I=MINT(IS)
        PARI(36+IS)=P(I,3)/VINT(1)
        PARI(38+IS)=P(I,4)/VINT(1)
        PR=MAX(1E-20,P(I,5)**2+P(I,1)**2+P(I,2)**2)
        PARI(40+IS)=SIGN(LOG(MIN((SQRT(PR+P(I,3)**2)+ABS(P(I,3)))/
     &  SQRT(PR),1E20)),P(I,3))
        PR=MAX(1E-20,P(I,1)**2+P(I,2)**2)
        PARI(42+IS)=SIGN(LOG(MIN((SQRT(PR+P(I,3)**2)+ABS(P(I,3)))/
     &  SQRT(PR),1E20)),P(I,3))
        PARI(44+IS)=P(I,3)/SQRT(1E-20+P(I,1)**2+P(I,2)**2+P(I,3)**2)
        PARI(46+IS)=ULANGL(P(I,3),SQRT(P(I,1)**2+P(I,2)**2))
        PARI(48+IS)=ULANGL(P(I,1),P(I,2))
  140   CONTINUE
      ENDIF

C...Store sum up transverse and longitudinal momenta.
      PARI(65)=2.*PARI(17)
      IF(ISUB.LE.90.OR.ISUB.GE.95) THEN
        DO 150 I=MSTP(126)+1,N
        IF(K(I,1).LE.0.OR.K(I,1).GT.10) GOTO 150
        PT=SQRT(P(I,1)**2+P(I,2)**2)
        PARI(69)=PARI(69)+PT
        IF(I.LE.MINT(52)) PARI(66)=PARI(66)+PT
        IF(I.GT.MINT(52).AND.I.LE.MINT(53)) PARI(68)=PARI(68)+PT
  150   CONTINUE
        PARI(67)=PARI(68)
        PARI(71)=VINT(151)
        PARI(72)=VINT(152)
        PARI(73)=VINT(151)
        PARI(74)=VINT(152)
      ELSE
        PARI(66)=PARI(65)
        PARI(69)=PARI(65)
      ENDIF

C...Store various other pieces of information into PARI.
      PARI(61)=VINT(148)
      PARI(81)=VINT(138)

C...Set information for LUTABU.
      IF(ISET(ISUB).EQ.1.OR.ISET(ISUB).EQ.3) THEN
        MSTU(161)=MINT(21)
        MSTU(162)=0
      ELSEIF(ISET(ISUB).EQ.5) THEN
        MSTU(161)=MINT(23)
        MSTU(162)=0
      ELSE
        MSTU(161)=MINT(21)
        MSTU(162)=MINT(22)
      ENDIF

      RETURN
      END

C*********************************************************************

      SUBROUTINE RYFRAM(IFRAME)

C...Performs transformations between different coordinate frames.
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/RYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/RYINT1/MINT(400),VINT(400)
      SAVE /LUDAT1/
      SAVE /RYPARS/,/RYINT1/

      IF(IFRAME.LT.1.OR.IFRAME.GT.2) THEN
        WRITE(MSTU(11),5000) IFRAME,MINT(6)
        RETURN
      ENDIF
      IF(IFRAME.EQ.MINT(6)) RETURN

      IF(MINT(6).EQ.1) THEN
C...Transform from fixed target or user specified frame to
C...CM-frame of incoming particles.
        CALL LUROBO(0.,0.,-VINT(8),-VINT(9),-VINT(10))
        CALL LUROBO(0.,-VINT(7),0.,0.,0.)
        CALL LUROBO(-VINT(6),0.,0.,0.,0.)
        MINT(6)=2

      ELSE
C...Transform from particle CM-frame to fixed target or user specified
C...frame.
        CALL LUROBO(VINT(6),VINT(7),VINT(8),VINT(9),VINT(10))
        MINT(6)=1
      ENDIF
      MSTI(6)=MINT(6)

 5000 FORMAT(1X,'Error: illegal values in subroutine RYFRAM.',1X,
     &'No transformation performed.'/1X,'IFRAME =',1X,I5,'; MINT(6) =',
     &1X,I5)

      RETURN
      END

C*********************************************************************

      SUBROUTINE RYWIDT(KFLR,SH,WDTP,WDTE)

C...Calculates full and partial widths of resonances.
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/LUDAT3/MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5)
      COMMON/RYSUBS/MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200)
      COMMON/RYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/RYINT1/MINT(400),VINT(400)
      COMMON/RYINT4/WIDP(21:40,0:40),WIDE(21:40,0:40),WIDS(21:40,3)
      SAVE /LUDAT1/,/LUDAT2/,/LUDAT3/
      SAVE /RYSUBS/,/RYPARS/,/RYINT1/,/RYINT4/
      DIMENSION WDTP(0:40),WDTE(0:40,0:5),MOFSV(3,2),WIDWSV(3,2),
     &WID2SV(3,2)
      SAVE MOFSV,WIDWSV,WID2SV
      DATA MOFSV/6*0/,WIDWSV/6*0./,WID2SV/6*0./

C...Some common constants.
      KFLA=IABS(KFLR)
      KFHIGG=25
      IHIGG=1
      IF(KFLA.EQ.35.OR.KFLA.EQ.36) THEN
        KFHIGG=KFLA
        IHIGG=KFLA-33
      ENDIF
      AEM=ULALEM(SH)
      XW=PARU(102)
      AS=ULALPS(SH)
      RADC=1.+AS/PARU(1)

C...Reset width information.
      DO 100 I=0,40
      WDTP(I)=0.
      DO 100 J=0,5
  100 WDTE(I,J)=0.

      IF(KFLA.EQ.21) THEN
C...QCD:
        DO 110 I=1,MDCY(21,3)
        IDC=I+MDCY(21,2)-1
        IF(MDME(IDC,1).LT.0) GOTO 110
        RM1=PMAS(IABS(KFDP(IDC,1)),1)**2/SH
        RM2=PMAS(IABS(KFDP(IDC,2)),1)**2/SH
        IF(SQRT(RM1)+SQRT(RM2).GT.1.) GOTO 110
        IF(I.LE.8) THEN
C...QCD -> q + q~
          WDTP(I)=(1.+2.*RM1)*SQRT(MAX(0.,1.-4.*RM1))
          WID2=1.
        ENDIF
        WDTP(0)=WDTP(0)+WDTP(I)
        IF(MDME(IDC,1).GT.0) THEN
          WDTE(I,MDME(IDC,1))=WDTP(I)*WID2
          WDTE(0,MDME(IDC,1))=WDTE(0,MDME(IDC,1))+WDTE(I,MDME(IDC,1))
          WDTE(I,0)=WDTE(I,MDME(IDC,1))
          WDTE(0,0)=WDTE(0,0)+WDTE(I,0)
        ENDIF
  110   CONTINUE

      ELSEIF(KFLA.EQ.22) THEN
C...QED photon.
        DO 120 I=1,MDCY(22,3)
        IDC=I+MDCY(22,2)-1
        IF(MDME(IDC,1).LT.0) GOTO 120
        RM1=PMAS(IABS(KFDP(IDC,1)),1)**2/SH
        RM2=PMAS(IABS(KFDP(IDC,2)),1)**2/SH
        IF(SQRT(RM1)+SQRT(RM2).GT.1.) GOTO 120
        IF(I.LE.8) THEN
C...QED -> q + q~.
          EF=KCHG(I,1)/3.
          FCOF=3.*RADC
          WDTP(I)=FCOF*EF**2*(1.+2.*RM1)*SQRT(MAX(0.,1.-4.*RM1))
          WID2=1.
        ELSEIF(I.LE.12) THEN
C...QED -> l+ + l-.
          EF=KCHG(9+2*(I-8),1)/3.
          WDTP(I)=EF**2*(1.+2.*RM1)*SQRT(MAX(0.,1.-4.*RM1))
          WID2=1.
        ENDIF
        WDTP(0)=WDTP(0)+WDTP(I)
        IF(MDME(IDC,1).GT.0) THEN
          WDTE(I,MDME(IDC,1))=WDTP(I)*WID2
          WDTE(0,MDME(IDC,1))=WDTE(0,MDME(IDC,1))+WDTE(I,MDME(IDC,1))
          WDTE(I,0)=WDTE(I,MDME(IDC,1))
          WDTE(0,0)=WDTE(0,0)+WDTE(I,0)
        ENDIF
  120   CONTINUE

      ELSEIF(KFLA.EQ.23) THEN
C...Z0:
        ICASE=1
        XWC=1./(16.*XW*(1.-XW))
        FACH=AEM/3.*XWC*SH
  130   CONTINUE
        IF(MINT(61).GE.1.AND.ICASE.EQ.2) THEN
          VINT(111)=0.
          VINT(112)=0.
          VINT(114)=0.
        ENDIF
        IF(MINT(61).EQ.1.AND.ICASE.EQ.2) THEN
          EI=KCHG(IABS(MINT(15)),1)/3.
          AI=SIGN(1.,EI)
          VI=AI-4.*EI*XW
          SQMZ=PMAS(23,1)**2
          HZ=FACH*WDTP(0)
          IF(MSTP(43).EQ.1.OR.MSTP(43).EQ.3) VINT(111)=1.
          IF(MSTP(43).EQ.3) VINT(112)=
     &    2.*XWC*SH*(SH-SQMZ)/((SH-SQMZ)**2+HZ**2)
          IF(MSTP(43).EQ.2.OR.MSTP(43).EQ.3) VINT(114)=
     &    XWC**2*SH**2/((SH-SQMZ)**2+HZ**2)
        ENDIF
        DO 140 I=1,MDCY(23,3)
        IDC=I+MDCY(23,2)-1
        IF(MDME(IDC,1).LT.0) GOTO 140
        RM1=PMAS(IABS(KFDP(IDC,1)),1)**2/SH
        RM2=PMAS(IABS(KFDP(IDC,2)),1)**2/SH
        IF(SQRT(RM1)+SQRT(RM2).GT.1.) GOTO 140
        IF(I.LE.8) THEN
C...Z0 -> q + q~
          EF=KCHG(I,1)/3.
          AF=SIGN(1.,EF+0.1)
          VF=AF-4.*EF*XW
          FCOF=3.*RADC
        ELSEIF(I.LE.16) THEN
C...Z0 -> l+ + l-, nu + nu~
          EF=KCHG(I+2,1)/3.
          AF=SIGN(1.,EF+0.1)
          VF=AF-4.*EF*XW
          FCOF=1.
        ENDIF
        BE34=SQRT(MAX(0.,1.-4.*RM1))
        IF(ICASE.EQ.1) THEN
          WDTP(I)=FCOF*(VF**2*(1.+2.*RM1)+AF**2*(1.-4.*RM1))*BE34
        ELSEIF(MINT(61).EQ.1.AND.ICASE.EQ.2) THEN
          WDTP(I)=FCOF*((EI**2*VINT(111)*EF**2+EI*VI*VINT(112)*
     &    EF*VF+(VI**2+AI**2)*VINT(114)*VF**2)*(1.+2.*RM1)+
     &    (VI**2+AI**2)*VINT(114)*AF**2*(1.-4.*RM1))*BE34
        ELSEIF(MINT(61).EQ.2.AND.ICASE.EQ.2) THEN
          FGGF=FCOF*EF**2*(1.+2.*RM1)*BE34
          FGZF=FCOF*EF*VF*(1.+2.*RM1)*BE34
          FZZF=FCOF*(VF**2*(1.+2.*RM1)+AF**2*(1.-4.*RM1))*BE34
        ENDIF
        WID2=1.
        IF(ICASE.EQ.1) WDTP(0)=WDTP(0)+WDTP(I)
        IF(MDME(IDC,1).GT.0) THEN
          IF((ICASE.EQ.1.AND.MINT(61).NE.1).OR.
     &    (ICASE.EQ.2.AND.MINT(61).EQ.1)) THEN
            WDTE(I,MDME(IDC,1))=WDTP(I)*WID2
            WDTE(0,MDME(IDC,1))=WDTE(0,MDME(IDC,1))+WDTE(I,MDME(IDC,1))
            WDTE(I,0)=WDTE(I,MDME(IDC,1))
            WDTE(0,0)=WDTE(0,0)+WDTE(I,0)
          ENDIF
          IF(MINT(61).EQ.2.AND.ICASE.EQ.2) THEN
            IF(MSTP(43).EQ.1.OR.MSTP(43).EQ.3) VINT(111)=
     &      VINT(111)+FGGF*WID2
            IF(MSTP(43).EQ.3) VINT(112)=VINT(112)+FGZF*WID2
            IF(MSTP(43).EQ.2.OR.MSTP(43).EQ.3) VINT(114)=
     &      VINT(114)+FZZF*WID2
          ENDIF
        ENDIF
  140   CONTINUE
        IF(MINT(61).GE.1) ICASE=3-ICASE
        IF(ICASE.EQ.2) GOTO 130

      ELSEIF(KFLA.EQ.24) THEN
C...W+/-:
        DO 150 I=1,MDCY(24,3)
        IDC=I+MDCY(24,2)-1
        IF(MDME(IDC,1).LT.0) GOTO 150
        RM1=PMAS(IABS(KFDP(IDC,1)),1)**2/SH
        RM2=PMAS(IABS(KFDP(IDC,2)),1)**2/SH
        IF(SQRT(RM1)+SQRT(RM2).GT.1.) GOTO 150
        IF(I.LE.16) THEN
C...W+/- -> q + q~'
          FCOF=3.*RADC*VCKM((I-1)/4+1,MOD(I-1,4)+1)
        ELSEIF(I.LE.20) THEN
C...W+/- -> l+/- + nu
          FCOF=1.
        ENDIF
        WDTP(I)=FCOF*(2.-RM1-RM2-(RM1-RM2)**2)*
     &  SQRT(MAX(0.,(1.-RM1-RM2)**2-4.*RM1*RM2))
        WID2=1.
        WDTP(0)=WDTP(0)+WDTP(I)
        IF(MDME(IDC,1).GT.0) THEN
          WDTE(I,MDME(IDC,1))=WDTP(I)*WID2
          WDTE(0,MDME(IDC,1))=WDTE(0,MDME(IDC,1))+WDTE(I,MDME(IDC,1))
          WDTE(I,0)=WDTE(I,MDME(IDC,1))
          WDTE(0,0)=WDTE(0,0)+WDTE(I,0)
        ENDIF
  150   CONTINUE

      ELSEIF(KFLA.EQ.25.OR.KFLA.EQ.35.OR.KFLA.EQ.36) THEN
C...H0 (or H'0, or A0):
        DO 190 I=1,MDCY(KFHIGG,3)
        IDC=I+MDCY(KFHIGG,2)-1
        IF(MDME(IDC,1).LT.0) GOTO 190
        RM1=PMAS(IABS(KFDP(IDC,1)),1)**2/SH
        RM2=PMAS(IABS(KFDP(IDC,2)),1)**2/SH
        IF(I.NE.16.AND.I.NE.17.AND.SQRT(RM1)+SQRT(RM2).GT.1.) GOTO 190
        IF(I.LE.8) THEN
C...H0 -> q + q~
          WDTP(I)=3.*RM1*(1.-4.*RM1)*SQRT(MAX(0.,1.-4.*RM1))*RADC
          IF(MSTP(37).EQ.1) WDTP(I)=WDTP(I)*
     &    (LOG(MAX(4.,PARP(37)**2*RM1*SH/PARU(117)**2))/
     &    LOG(MAX(4.,SH/PARU(117)**2)))**(24./(33.-2.*MSTU(118)))
          IF(MSTP(4).GE.1.OR.IHIGG.GE.2) THEN
            IF(MOD(I,2).EQ.1) WDTP(I)=WDTP(I)*PARU(151+10*IHIGG)**2
            IF(MOD(I,2).EQ.0) WDTP(I)=WDTP(I)*PARU(152+10*IHIGG)**2
          ENDIF
          WID2=1.
        ELSEIF(I.LE.12) THEN
C...H0 -> l+ + l-
          WDTP(I)=RM1*(1.-4.*RM1)*SQRT(MAX(0.,1.-4.*RM1))
          IF(MSTP(4).GE.1.OR.IHIGG.GE.2) WDTP(I)=WDTP(I)*
     &    PARU(153+10*IHIGG)**2
          WID2=1.
        ELSEIF(I.EQ.13) THEN
C...H0 -> g + g; quark loop contribution only
          ETARE=0.
          ETAIM=0.
          DO 160 J=1,2*MSTP(1)
          EPS=(2.*PMAS(J,1))**2/SH
          IF(EPS.LE.1.) THEN
            IF(EPS.GT.1.E-4) THEN
              ROOT=SQRT(1.-EPS)
              RLN=LOG((1.+ROOT)/(1.-ROOT))
            ELSE
              RLN=LOG(4./EPS-2.)
            ENDIF
            PHIRE=0.25*(RLN**2-PARU(1)**2)
            PHIIM=0.5*PARU(1)*RLN
          ELSE
            PHIRE=-(ASIN(1./SQRT(EPS)))**2
            PHIIM=0.
          ENDIF
          ETAREJ=0.5*EPS*(1.+(EPS-1.)*PHIRE)
          ETAIMJ=0.5*EPS*(EPS-1.)*PHIIM
          IF(MSTP(4).GE.1.OR.IHIGG.GE.2.AND.MOD(J,2).EQ.1) THEN
            ETAREJ=ETAREJ*PARU(151+10*IHIGG)
            ETAIMJ=ETAIMJ*PARU(151+10*IHIGG)
          ELSEIF(MSTP(4).GE.1.OR.IHIGG.GE.2) THEN
            ETAREJ=ETAREJ*PARU(152+10*IHIGG)
            ETAIMJ=ETAIMJ*PARU(152+10*IHIGG)
          ENDIF
          ETARE=ETARE+ETAREJ
          ETAIM=ETAIM+ETAIMJ
  160     CONTINUE
          ETA2=ETARE**2+ETAIM**2
          WDTP(I)=(AS/PARU(1))**2*ETA2
          WID2=1.
        ELSEIF(I.EQ.14) THEN
C...H0 -> gamma + gamma; quark, charged lepton and W loop contributions
          ETARE=0.
          ETAIM=0.
          DO 170 J=1,3*MSTP(1)+1
          IF(J.LE.2*MSTP(1)) THEN
            EJ=KCHG(J,1)/3.
            EPS=(2.*PMAS(J,1))**2/SH
          ELSEIF(J.LE.3*MSTP(1)) THEN
            JL=2*(J-2*MSTP(1))-1
            EJ=KCHG(10+JL,1)/3.
            EPS=(2.*PMAS(10+JL,1))**2/SH
          ELSE
            EPS=(2.*PMAS(24,1))**2/SH
          ENDIF
          IF(EPS.LE.1.) THEN
            IF(EPS.GT.1.E-4) THEN
              ROOT=SQRT(1.-EPS)
              RLN=LOG((1.+ROOT)/(1.-ROOT))
            ELSE
              RLN=LOG(4./EPS-2.)
            ENDIF
            PHIRE=0.25*(RLN**2-PARU(1)**2)
            PHIIM=0.5*PARU(1)*RLN
          ELSE
            PHIRE=-(ASIN(1./SQRT(EPS)))**2
            PHIIM=0.
          ENDIF
          IF(J.LE.2*MSTP(1)) THEN
            ETAREJ=0.5*3.*EJ**2*EPS*(1.+(EPS-1.)*PHIRE)
            ETAIMJ=0.5*3.*EJ**2*EPS*(EPS-1.)*PHIIM
            IF(MSTP(4).GE.1.OR.IHIGG.GE.2.AND.MOD(J,2).EQ.1) THEN
              ETAREJ=ETAREJ*PARU(151+10*IHIGG)
              ETAIMJ=ETAIMJ*PARU(151+10*IHIGG)
            ELSEIF(MSTP(4).GE.1.OR.IHIGG.GE.2) THEN
              ETAREJ=ETAREJ*PARU(152+10*IHIGG)
              ETAIMJ=ETAIMJ*PARU(152+10*IHIGG)
            ENDIF
          ELSEIF(J.LE.3*MSTP(1)) THEN
            ETAREJ=0.5*EJ**2*EPS*(1.+(EPS-1.)*PHIRE)
            ETAIMJ=0.5*EJ**2*EPS*(EPS-1.)*PHIIM
            IF(MSTP(4).GE.1.OR.IHIGG.GE.2) THEN
              ETAREJ=ETAREJ*PARU(153+10*IHIGG)
              ETAIMJ=ETAIMJ*PARU(153+10*IHIGG)
            ENDIF
          ELSE
            ETAREJ=-0.5-0.75*EPS*(1.+(EPS-2.)*PHIRE)
            ETAIMJ=0.75*EPS*(EPS-2.)*PHIIM
            IF(MSTP(4).GE.1.OR.IHIGG.GE.2) THEN
              ETAREJ=ETAREJ*PARU(155+10*IHIGG)
              ETAIMJ=ETAIMJ*PARU(155+10*IHIGG)
            ENDIF
          ENDIF
          ETARE=ETARE+ETAREJ
          ETAIM=ETAIM+ETAIMJ
  170     CONTINUE
          ETA2=ETARE**2+ETAIM**2
          WDTP(I)=(AEM/PARU(1))**2*0.5*ETA2
          WID2=1.
        ELSEIF(I.EQ.15) THEN
C...H0 -> gamma + Z0; quark, charged lepton and W loop contributions
          ETARE=0.
          ETAIM=0.
          DO 180 J=1,3*MSTP(1)+1
          IF(J.LE.2*MSTP(1)) THEN
            EJ=KCHG(J,1)/3.
            AJ=SIGN(1.,EJ+0.1)
            VJ=AJ-4.*EJ*XW
            EPS=(2.*PMAS(J,1))**2/SH
            EPSP=(2.*PMAS(J,1)/PMAS(23,1))**2
          ELSEIF(J.LE.3*MSTP(1)) THEN
            JL=2*(J-2*MSTP(1))-1
            EJ=KCHG(10+JL,1)/3.
            AJ=SIGN(1.,EJ+0.1)
            VJ=AJ-4.*EJ*XW
            EPS=(2.*PMAS(10+JL,1))**2/SH
            EPSP=(2.*PMAS(10+JL,1)/PMAS(23,1))**2
          ELSE
            EPS=(2.*PMAS(24,1))**2/SH
            EPSP=(2.*PMAS(24,1)/PMAS(23,1))**2
          ENDIF
          IF(EPS.LE.1.) THEN
            ROOT=SQRT(1.-EPS)
            IF(EPS.GT.1.E-4) THEN
              RLN=LOG((1.+ROOT)/(1.-ROOT))
            ELSE
              RLN=LOG(4./EPS-2.)
            ENDIF
            PHIRE=0.25*(RLN**2-PARU(1)**2)
            PHIIM=0.5*PARU(1)*RLN
            PSIRE=-(1.+0.5*ROOT*RLN)
            PSIIM=0.5*PARU(1)*ROOT
          ELSE
            PHIRE=-(ASIN(1./SQRT(EPS)))**2
            PHIIM=0.
            PSIRE=-(1.+SQRT(EPS-1.)*ASIN(1./SQRT(EPS)))
            PSIIM=0.
          ENDIF
          IF(EPSP.LE.1.) THEN
            ROOT=SQRT(1.-EPSP)
            IF(EPSP.GT.1.E-4) THEN
              RLN=LOG((1.+ROOT)/(1.-ROOT))
            ELSE
              RLN=LOG(4./EPSP-2.)
            ENDIF
            PHIREP=0.25*(RLN**2-PARU(1)**2)
            PHIIMP=0.5*PARU(1)*RLN
            PSIREP=-(1.+0.5*ROOT*RLN)
            PSIIMP=0.5*PARU(1)*ROOT
          ELSE
            PHIREP=-(ASIN(1./SQRT(EPSP)))**2
            PHIIMP=0.
            PSIREP=-(1.+SQRT(EPSP-1.)*ASIN(1./SQRT(EPSP)))
            PSIIMP=0.
          ENDIF
          FXYRE=EPS*EPSP/(8.*(EPS-EPSP))*(1.-EPS*EPSP/(EPS-EPSP)*(PHIRE-
     &    PHIREP)-2.*EPS/(EPS-EPSP)*(PSIRE-PSIREP))
          FXYIM=EPS*EPSP/(8.*(EPS-EPSP))*(-EPS*EPSP/(EPS-EPSP)*(PHIIM-
     &    PHIIMP)-2.*EPS/(EPS-EPSP)*(PSIIM-PSIIMP))
          F1RE=EPS*EPSP/(2.*(EPS-EPSP))*(PHIRE-PHIREP)
          F1IM=EPS*EPSP/(2.*(EPS-EPSP))*(PHIIM-PHIIMP)
          IF(J.LE.2*MSTP(1)) THEN
            ETAREJ=-3.*EJ*VJ*(FXYRE-0.25*F1RE)
            ETAIMJ=-3.*EJ*VJ*(FXYIM-0.25*F1IM)
            IF(MSTP(4).GE.1.OR.IHIGG.GE.2.AND.MOD(J,2).EQ.1) THEN
              ETAREJ=ETAREJ*PARU(151+10*IHIGG)
              ETAIMJ=ETAIMJ*PARU(151+10*IHIGG)
            ELSEIF(MSTP(4).GE.1.OR.IHIGG.GE.2) THEN
              ETAREJ=ETAREJ*PARU(152+10*IHIGG)
              ETAIMJ=ETAIMJ*PARU(152+10*IHIGG)
            ENDIF
          ELSEIF(J.LE.3*MSTP(1)) THEN
            ETAREJ=-EJ*VJ*(FXYRE-0.25*F1RE)
            ETAIMJ=-EJ*VJ*(FXYIM-0.25*F1IM)
            IF(MSTP(4).GE.1.OR.IHIGG.GE.2) THEN
              ETAREJ=ETAREJ*PARU(153+10*IHIGG)
              ETAIMJ=ETAIMJ*PARU(153+10*IHIGG)
            ENDIF
          ELSE
            ETAREJ=-(1.-XW)*(((1.+2./EPS)*XW/(1.-XW)-
     &      (5.+2./EPS))*FXYRE+(3.-XW/(1.-XW))*F1RE)
            ETAIMJ=-(1.-XW)*(((1.+2./EPS)*XW/(1.-XW)-
     &      (5.+2./EPS))*FXYIM+(3.-XW/(1.-XW))*F1IM)
            IF(MSTP(4).GE.1.OR.IHIGG.GE.2) THEN
              ETAREJ=ETAREJ*PARU(155+10*IHIGG)
              ETAIMJ=ETAIMJ*PARU(155+10*IHIGG)
            ENDIF
          ENDIF
          ETARE=ETARE+ETAREJ
          ETAIM=ETAIM+ETAIMJ
  180     CONTINUE
          ETA2=(ETARE**2+ETAIM**2)/(1.-XW)
          WDTP(I)=(AEM/PARU(1))**2*(1.-PMAS(23,1)**2/SH)**3/XW*ETA2
          WID2=WIDS(23,2)
        ELSEIF(I.LE.17) THEN
C...H0 -> Z0 + Z0, W+ + W-
          PM1=PMAS(IABS(KFDP(IDC,1)),1)
          PG1=PMAS(IABS(KFDP(IDC,1)),2)
          IF(MINT(62).GE.1) THEN
            IF(MSTP(42).EQ.0.OR.(4.*(PM1+10.*PG1)**2.LT.SH.AND.
     &      CKIN(46).LT.CKIN(45).AND.CKIN(48).LT.CKIN(47).AND.
     &      MAX(CKIN(45),CKIN(47)).LT.PM1-10.*PG1)) THEN
              MOFSV(IHIGG,I-15)=0
              WIDW=(1.-4.*RM1+12.*RM1**2)*SQRT(MAX(0.,1.-4.*RM1))
              WID2=1.
            ELSE
              MOFSV(IHIGG,I-15)=1
              RMAS=SQRT(MAX(0.,SH))
              CALL RYOFSH(1,KFLA,KFDP(IDC,1),KFDP(IDC,2),RMAS,WIDW,WID2)
              WIDWSV(IHIGG,I-15)=WIDW
              WID2SV(IHIGG,I-15)=WID2
            ENDIF
          ELSE
            IF(MOFSV(IHIGG,I-15).EQ.0) THEN
              WIDW=(1.-4.*RM1+12.*RM1**2)*SQRT(MAX(0.,1.-4.*RM1))
              WID2=1.
            ELSE
              WIDW=WIDWSV(IHIGG,I-15)
              WID2=WID2SV(IHIGG,I-15)
            ENDIF
          ENDIF
          WDTP(I)=WIDW/(2.*(18-I))
          IF(MSTP(4).GE.1.OR.IHIGG.GE.2) WDTP(I)=WDTP(I)*
     &    PARU(138+I+10*IHIGG)**2
          WID2=WID2*WIDS(7+I,1)
        ELSEIF(I.EQ.18.AND.KFLA.EQ.35) THEN
C...H'0 -> H0 + H0.
          WDTP(I)=PARU(176)**2*0.25*PMAS(23,1)**4/SH**2*
     &    SQRT(MAX(0.,1.-4.*RM1))
          WID2=WIDS(25,2)**2
        ELSEIF(I.EQ.19.AND.KFLA.EQ.35) THEN
C...H'0 -> A0 + A0.
          WDTP(I)=PARU(177)**2*0.25*PMAS(23,1)**4/SH**2*
     &    SQRT(MAX(0.,1.-4.*RM1))
          WID2=WIDS(36,2)**2
        ELSEIF(I.EQ.18.AND.KFLA.EQ.36) THEN
C...A0 -> Z0 + H0.
          WDTP(I)=PARU(186)**2*0.5*SQRT(MAX(0.,(1.-RM1-RM2)**2-
     &    4.*RM1*RM2))**3
          WID2=WIDS(23,2)*WIDS(25,2)
        ENDIF
        WDTP(0)=WDTP(0)+WDTP(I)
        IF(MDME(IDC,1).GT.0) THEN
          WDTE(I,MDME(IDC,1))=WDTP(I)*WID2
          WDTE(0,MDME(IDC,1))=WDTE(0,MDME(IDC,1))+WDTE(I,MDME(IDC,1))
          WDTE(I,0)=WDTE(I,MDME(IDC,1))
          WDTE(0,0)=WDTE(0,0)+WDTE(I,0)
        ENDIF
  190   CONTINUE

      ELSEIF(KFLA.EQ.32) THEN
C...Z'0:
        ICASE=1
        XWC=1./(16.*XW*(1.-XW))
        FACH=AEM/3.*XWC*SH
        VINT(117)=0.
  200   CONTINUE
        IF(MINT(61).GE.1.AND.ICASE.EQ.2) THEN
          VINT(111)=0.
          VINT(112)=0.
          VINT(113)=0.
          VINT(114)=0.
          VINT(115)=0.
          VINT(116)=0.
        ENDIF
        IF(MINT(61).EQ.1.AND.ICASE.EQ.2) THEN
          KFAI=IABS(MINT(15))
          EI=KCHG(KFAI,1)/3.
          AI=SIGN(1.,EI+0.1)
          VI=AI-4.*EI*XW
          KFAIC=1
          IF(KFAI.LE.10.AND.MOD(KFAI,2).EQ.0) KFAIC=2
          IF(KFAI.GT.10.AND.MOD(KFAI,2).NE.0) KFAIC=3
          IF(KFAI.GT.10.AND.MOD(KFAI,2).EQ.0) KFAIC=4
          VPI=PARU(119+2*KFAIC)
          API=PARU(120+2*KFAIC)
          SQMZ=PMAS(23,1)**2
          HZ=FACH*VINT(117)
          SQMZP=PMAS(32,1)**2
          HZP=FACH*WDTP(0)
          IF(MSTP(44).EQ.1.OR.MSTP(44).EQ.4.OR.MSTP(44).EQ.5.OR.
     &    MSTP(44).EQ.7) VINT(111)=1.
          IF(MSTP(44).EQ.4.OR.MSTP(44).EQ.7) VINT(112)=
     &    2.*XWC*SH*(SH-SQMZ)/((SH-SQMZ)**2+HZ**2)
          IF(MSTP(44).EQ.5.OR.MSTP(44).EQ.7) VINT(113)=
     &    2.*XWC*SH*(SH-SQMZP)/((SH-SQMZP)**2+HZP**2)
          IF(MSTP(44).EQ.2.OR.MSTP(44).EQ.4.OR.MSTP(44).EQ.6.OR.
     &    MSTP(44).EQ.7) VINT(114)=XWC**2*SH**2/((SH-SQMZ)**2+HZ**2)
          IF(MSTP(44).EQ.6.OR.MSTP(44).EQ.7) VINT(115)=
     &    2.*XWC**2*SH**2*((SH-SQMZ)*(SH-SQMZP)+HZ*HZP)/
     &    (((SH-SQMZ)**2+HZ**2)*((SH-SQMZP)**2+HZP**2))
          IF(MSTP(44).EQ.3.OR.MSTP(44).EQ.5.OR.MSTP(44).EQ.6.OR.
     &    MSTP(44).EQ.7) VINT(116)=XWC**2*SH**2/((SH-SQMZP)**2+HZP**2)
        ENDIF
        DO 210 I=1,MDCY(32,3)
        IDC=I+MDCY(32,2)-1
        IF(MDME(IDC,1).LT.0) GOTO 210
        RM1=PMAS(IABS(KFDP(IDC,1)),1)**2/SH
        RM2=PMAS(IABS(KFDP(IDC,2)),1)**2/SH
        IF(SQRT(RM1)+SQRT(RM2).GT.1..OR.MDME(IDC,1).LT.0) GOTO 210
        IF(I.LE.16) THEN
          IF(I.LE.8) THEN
C...Z'0 -> q + q~
            EF=KCHG(I,1)/3.
            AF=SIGN(1.,EF+0.1)
            VF=AF-4.*EF*XW
            VPF=PARU(123-2*MOD(I,2))
            APF=PARU(124-2*MOD(I,2))
            FCOF=3.*RADC
          ELSEIF(I.LE.16) THEN
C...Z'0 -> l+ + l-, nu + nu~
            EF=KCHG(I+2,1)/3.
            AF=SIGN(1.,EF+0.1)
            VF=AF-4.*EF*XW
            VPF=PARU(127-2*MOD(I,2))
            APF=PARU(128-2*MOD(I,2))
            FCOF=1.
          ENDIF
          BE34=SQRT(MAX(0.,1.-4.*RM1))
          IF(ICASE.EQ.1) THEN
            WDTPZ=FCOF*(VF**2*(1.+2.*RM1)+AF**2*(1.-4.*RM1))*BE34
            WDTP(I)=FCOF*(VPF**2*(1.+2.*RM1)+APF**2*(1.-4.*RM1))*BE34
          ELSEIF(MINT(61).EQ.1.AND.ICASE.EQ.2) THEN
            WDTP(I)=FCOF*((EI**2*VINT(111)*EF**2+EI*VI*VINT(112)*
     &      EF*VF+EI*VPI*VINT(113)*EF*VPF+(VI**2+AI**2)*VINT(114)*
     &      VF**2+(VI*VPI+AI*API)*VINT(115)*VF*VPF+(VPI**2+API**2)*
     &      VINT(116)*VPF**2)*(1.+2.*RM1)+((VI**2+AI**2)*VINT(114)*
     &      AF**2+(VI*VPI+AI*API)*VINT(115)*AF*APF+(VPI**2+API**2)*
     &      VINT(116)*APF**2)*(1.-4.*RM1))*BE34
          ELSEIF(MINT(61).EQ.2) THEN
            FGGF=FCOF*EF**2*(1.+2.*RM1)*BE34
            FGZF=FCOF*EF*VF*(1.+2.*RM1)*BE34
            FGZPF=FCOF*EF*VPF*(1.+2.*RM1)*BE34
            FZZF=FCOF*(VF**2*(1.+2.*RM1)+AF**2*(1.-4.*RM1))*BE34
            FZZPF=FCOF*(VF*VPF*(1.+2.*RM1)+AF*APF*(1.-4.*RM1))*BE34
            FZPZPF=FCOF*(VPF**2*(1.+2.*RM1)+APF**2*(1.-4.*RM1))*BE34
          ENDIF
          WID2=1.
        ELSEIF(I.EQ.17) THEN
C...Z'0 -> W+ + W-
          WDTPZP=PARU(129)**2*(1.-XW)**2*
     &    SQRT(MAX(0.,(1.-RM1-RM2)**2-4.*RM1*RM2))**3*
     &    (1.+10.*RM1+10.*RM2+RM1**2+RM2**2+10.*RM1*RM2)
          IF(ICASE.EQ.1) THEN
            WDTPZ=0.
            WDTP(I)=WDTPZP
          ELSEIF(MINT(61).EQ.1.AND.ICASE.EQ.2) THEN
            WDTP(I)=(VPI**2+API**2)*VINT(116)*WDTPZP
          ELSEIF(MINT(61).EQ.2) THEN
            FGGF=0.
            FGZF=0.
            FGZPF=0.
            FZZF=0.
            FZZPF=0.
            FZPZPF=WDTPZP
          ENDIF
          WID2=WIDS(24,1)
        ELSEIF(I.EQ.18) THEN
C...Z'0 -> H+ + H-
          CZC=2.*(1.-2.*XW)
          BE34C=(1.-4.*RM1)*SQRT(MAX(0.,1.-4.*RM1))
          IF(ICASE.EQ.1) THEN
            WDTPZ=0.25*PARU(142)**2*CZC**2*BE34C
            WDTP(I)=0.25*PARU(143)**2*CZC**2*BE34C
          ELSEIF(MINT(61).EQ.1.AND.ICASE.EQ.2) THEN
            WDTP(I)=0.25*(EI**2*VINT(111)+PARU(142)*EI*VI*VINT(112)*
     &      CZC+PARU(143)*EI*VPI*VINT(113)*CZC+PARU(142)**2*
     &      (VI**2+AI**2)*VINT(114)*CZC**2+PARU(142)*PARU(143)*
     &      (VI*VPI+AI*API)*VINT(115)*CZC**2+PARU(143)**2*
     &      (VPI**2+API**2)*VINT(116)*CZC**2)*BE34C
          ELSEIF(MINT(61).EQ.2) THEN
            FGGF=0.25*BE34C
            FGZF=0.25*PARU(142)*CZC*BE34C
            FGZPF=0.25*PARU(143)*CZC*BE34C
            FZZF=0.25*PARU(142)**2*CZC**2*BE34C
            FZZPF=0.25*PARU(142)*PARU(143)*CZC**2*BE34C
            FZPZPF=0.25*PARU(143)**2*CZC**2*BE34C
          ENDIF
          WID2=WIDS(37,1)
        ELSEIF(I.EQ.19) THEN
C...Z'0 -> Z0 + gamma.
        ELSEIF(I.EQ.20) THEN
C...Z'0 -> Z0 + H0
          FLAM=SQRT(MAX(0.,(1.-RM1-RM2)**2-4.*RM1*RM2))
          WDTPZP=PARU(145)**2*4.*ABS(1.-2.*XW)*(3.*RM1+0.25*FLAM**2)*
     &    FLAM
          IF(ICASE.EQ.1) THEN
            WDTPZ=0.
            WDTP(I)=WDTPZP
          ELSEIF(MINT(61).EQ.1.AND.ICASE.EQ.2) THEN
            WDTP(I)=(VPI**2+API**2)*VINT(116)*WDTPZP
          ELSEIF(MINT(61).EQ.2) THEN
            FGGF=0.
            FGZF=0.
            FGZPF=0.
            FZZF=0.
            FZZPF=0.
            FZPZPF=WDTPZP
          ENDIF
          WID2=WIDS(23,2)*WIDS(25,2)
        ELSEIF(I.EQ.21.OR.I.EQ.22) THEN
C...Z' -> H0 + A0 or H'0 + A0.
          BE34C=SQRT(MAX(0.,(1.-RM1-RM2)**2-4.*RM1*RM2))**3
          IF(I.EQ.21) THEN
            CZAH=PARU(186)
            CZPAH=PARU(188)
          ELSE
            CZAH=PARU(187)
            CZPAH=PARU(189)
          ENDIF
          IF(ICASE.EQ.1) THEN
            WDTPZ=CZAH**2*BE34C
            WDTP(I)=CZPAH**2*BE34C
          ELSEIF(MINT(61).EQ.1.AND.ICASE.EQ.2) THEN
            WDTP(I)=(CZAH**2*(VI**2+AI**2)*VINT(114)+CZAH*CZPAH*
     &      (VI*VPI+AI*API)*VINT(115)+CZPAH**2*(VPI**2+API**2)*
     &      VINT(116))*BE34C
          ELSEIF(MINT(61).EQ.2) THEN
            FGGF=0.
            FGZF=0.
            FGZPF=0.
            FZZF=CZAH**2*BE34C
            FZZPF=CZAH*CZPAH*BE34C
            FZPZPF=CZPAH**2*BE34C
          ENDIF
          IF(I.EQ.21) WID2=WIDS(25,2)*WIDS(36,2)
          IF(I.EQ.22) WID2=WIDS(35,2)*WIDS(36,2)
        ENDIF
        IF(ICASE.EQ.1) THEN
          VINT(117)=VINT(117)+WDTPZ
          WDTP(0)=WDTP(0)+WDTP(I)
        ENDIF
        IF(MDME(IDC,1).GT.0) THEN
          IF((ICASE.EQ.1.AND.MINT(61).NE.1).OR.
     &    (ICASE.EQ.2.AND.MINT(61).EQ.1)) THEN
            WDTE(I,MDME(IDC,1))=WDTP(I)*WID2
            WDTE(0,MDME(IDC,1))=WDTE(0,MDME(IDC,1))+WDTE(I,MDME(IDC,1))
            WDTE(I,0)=WDTE(I,MDME(IDC,1))
            WDTE(0,0)=WDTE(0,0)+WDTE(I,0)
          ENDIF
          IF(MINT(61).EQ.2.AND.ICASE.EQ.2) THEN
            IF(MSTP(44).EQ.1.OR.MSTP(44).EQ.4.OR.MSTP(44).EQ.5.OR.
     &      MSTP(44).EQ.7) VINT(111)=VINT(111)+FGGF*WID2
            IF(MSTP(44).EQ.4.OR.MSTP(44).EQ.7) VINT(112)=VINT(112)+
     &      FGZF*WID2
            IF(MSTP(44).EQ.5.OR.MSTP(44).EQ.7) VINT(113)=VINT(113)+
     &      FGZPF*WID2
            IF(MSTP(44).EQ.2.OR.MSTP(44).EQ.4.OR.MSTP(44).EQ.6.OR.
     &      MSTP(44).EQ.7) VINT(114)=VINT(114)+FZZF*WID2
            IF(MSTP(44).EQ.6.OR.MSTP(44).EQ.7) VINT(115)=VINT(115)+
     &      FZZPF*WID2
            IF(MSTP(44).EQ.3.OR.MSTP(44).EQ.5.OR.MSTP(44).EQ.6.OR.
     &      MSTP(44).EQ.7) VINT(116)=VINT(116)+FZPZPF*WID2
          ENDIF
        ENDIF
  210   CONTINUE
        IF(MINT(61).GE.1) ICASE=3-ICASE
        IF(ICASE.EQ.2) GOTO 200

      ELSEIF(KFLA.EQ.34) THEN
C...W'+/-:
        DO 220 I=1,MDCY(34,3)
        IDC=I+MDCY(34,2)-1
        IF(MDME(IDC,1).LT.0) GOTO 220
        RM1=PMAS(IABS(KFDP(IDC,1)),1)**2/SH
        RM2=PMAS(IABS(KFDP(IDC,2)),1)**2/SH
        IF(SQRT(RM1)+SQRT(RM2).GT.1.) GOTO 220
        IF(I.LE.20) THEN
          IF(I.LE.16) THEN
C...W'+/- -> q + q~'
            FCOF=3.*RADC*(PARU(131)**2+PARU(132)**2)*
     &      VCKM((I-1)/4+1,MOD(I-1,4)+1)
          ELSEIF(I.LE.20) THEN
C...W'+/- -> l+/- + nu
            FCOF=PARU(133)**2+PARU(134)**2
          ENDIF
          WDTP(I)=FCOF*0.5*(2.-RM1-RM2-(RM1-RM2)**2)*
     &    SQRT(MAX(0.,(1.-RM1-RM2)**2-4.*RM1*RM2))
          WID2=1.
        ELSEIF(I.EQ.21) THEN
C...W'+/- -> W+/- + Z0
          WDTP(I)=PARU(135)**2*0.5*(1.-XW)*(RM1/RM2)*
     &    SQRT(MAX(0.,(1.-RM1-RM2)**2-4.*RM1*RM2))**3*
     &    (1.+10.*RM1+10.*RM2+RM1**2+RM2**2+10.*RM1*RM2)
          IF(KFLR.GT.0) WID2=WIDS(24,2)*WIDS(23,2)
          IF(KFLR.LT.0) WID2=WIDS(24,3)*WIDS(23,2)
        ELSEIF(I.EQ.23) THEN
C...W'+/- -> W+/- + H0
          FLAM=SQRT(MAX(0.,(1.-RM1-RM2)**2-4.*RM1*RM2))
          WDTP(I)=PARU(146)**2*2.*(3.*RM1+0.25*FLAM**2)*FLAM
          IF(KFLR.GT.0) WID2=WIDS(24,2)*WIDS(25,2)
          IF(KFLR.LT.0) WID2=WIDS(24,3)*WIDS(25,2)
        ENDIF
        WDTP(0)=WDTP(0)+WDTP(I)
        IF(MDME(IDC,1).GT.0) THEN
          WDTE(I,MDME(IDC,1))=WDTP(I)*WID2
          WDTE(0,MDME(IDC,1))=WDTE(0,MDME(IDC,1))+WDTE(I,MDME(IDC,1))
          WDTE(I,0)=WDTE(I,MDME(IDC,1))
          WDTE(0,0)=WDTE(0,0)+WDTE(I,0)
        ENDIF
  220   CONTINUE

      ELSEIF(KFLA.EQ.37) THEN
C...H+/-:
        DO 230 I=1,MDCY(37,3)
        IDC=I+MDCY(37,2)-1
        IF(MDME(IDC,1).LT.0) GOTO 230
        RM1=PMAS(IABS(KFDP(IDC,1)),1)**2/SH
        RM2=PMAS(IABS(KFDP(IDC,2)),1)**2/SH
        IF(SQRT(RM1)+SQRT(RM2).GT.1.) GOTO 230
        IF(I.LE.4) THEN
C...H+/- -> q + q~'
          WDTP(I)=3.*RADC*((RM1*PARU(141)**2+RM2/PARU(141)**2)*
     &    (1.-RM1-RM2)-4.*RM1*RM2)*
     &    SQRT(MAX(0.,(1.-RM1-RM2)**2-4.*RM1*RM2))
          WID2=1.
        ELSEIF(I.LE.8) THEN
C...H+/- -> l+/- + nu
          WDTP(I)=((RM1*PARU(141)**2+RM2/PARU(141)**2)*(1.-RM1-RM2)-
     &    4.*RM1*RM2)*SQRT(MAX(0.,(1.-RM1-RM2)**2-4.*RM1*RM2))
          WID2=1.
        ELSEIF(I.EQ.9) THEN
C...H+/- -> W+/- + H0.
          WDTP(I)=PARU(195)**2*0.5*SQRT(MAX(0.,(1.-RM1-RM2)**2-
     &    4.*RM1*RM2))**3
          IF(KFLR.GT.0) WID2=WIDS(24,2)*WIDS(25,2)
          IF(KFLR.LT.0) WID2=WIDS(24,3)*WIDS(25,2)
        ENDIF
        WDTP(0)=WDTP(0)+WDTP(I)
        IF(MDME(IDC,1).GT.0) THEN
          WDTE(I,MDME(IDC,1))=WDTP(I)*WID2
          WDTE(0,MDME(IDC,1))=WDTE(0,MDME(IDC,1))+WDTE(I,MDME(IDC,1))
          WDTE(I,0)=WDTE(I,MDME(IDC,1))
          WDTE(0,0)=WDTE(0,0)+WDTE(I,0)
        ENDIF
  230   CONTINUE

      ELSEIF(KFLA.EQ.39) THEN
C...LQ (leptoquark).
        DO 240 I=1,MDCY(39,3)
        IDC=I+MDCY(39,2)-1
        IF(MDME(IDC,1).LT.0) GOTO 240
        RM1=PMAS(IABS(KFDP(IDC,1)),1)**2/SH
        RM2=PMAS(IABS(KFDP(IDC,2)),1)**2/SH
        IF(SQRT(RM1)+SQRT(RM2).GT.1.) GOTO 240
        WDTP(I)=PARU(151)*SQRT(MAX(0.,(1.-RM1-RM2)**2-4.*RM1*RM2))**3
        WID2=1.
        WDTP(0)=WDTP(0)+WDTP(I)
        IF(MDME(IDC,1).GT.0) THEN
          WDTE(I,MDME(IDC,1))=WDTP(I)*WID2
          WDTE(0,MDME(IDC,1))=WDTE(0,MDME(IDC,1))+WDTE(I,MDME(IDC,1))
          WDTE(I,0)=WDTE(I,MDME(IDC,1))
          WDTE(0,0)=WDTE(0,0)+WDTE(I,0)
        ENDIF
  240   CONTINUE

      ELSEIF(KFLA.EQ.40) THEN
C...R:
        DO 250 I=1,MDCY(40,3)
        IDC=I+MDCY(40,2)-1
        IF(MDME(IDC,1).LT.0) GOTO 250
        RM1=PMAS(IABS(KFDP(IDC,1)),1)**2/SH
        RM2=PMAS(IABS(KFDP(IDC,2)),1)**2/SH
        IF(SQRT(RM1)+SQRT(RM2).GT.1.) GOTO 250
        IF(I.LE.6) THEN
C...R -> q + q~'
          FCOF=3.*RADC
        ELSEIF(I.LE.9) THEN
C...R -> l+ + l'-
          FCOF=1.
        ENDIF
        WDTP(I)=FCOF*(2.-RM1-RM2-(RM1-RM2)**2)*
     &  SQRT(MAX(0.,(1.-RM1-RM2)**2-4.*RM1*RM2))
        WID2=1.
        WDTP(0)=WDTP(0)+WDTP(I)
        IF(MDME(IDC,1).GT.0) THEN
          WDTE(I,MDME(IDC,1))=WDTP(I)*WID2
          WDTE(0,MDME(IDC,1))=WDTE(0,MDME(IDC,1))+WDTE(I,MDME(IDC,1))
          WDTE(I,0)=WDTE(I,MDME(IDC,1))
          WDTE(0,0)=WDTE(0,0)+WDTE(I,0)
        ENDIF
  250   CONTINUE

      ENDIF
      MINT(61)=0
      MINT(62)=0

      RETURN
      END

C***********************************************************************

      SUBROUTINE RYOFSH(MOFSH,KFMO,KFD1,KFD2,PMMO,RET1,RET2)

C...Calculates partial width and differential cross-section maxima
C...of channels/processes not allowed on mass-shell, and selects
C...masses in such channels/processes.
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/LUDAT3/MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5)
      COMMON/RYSUBS/MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200)
      COMMON/RYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/RYINT1/MINT(400),VINT(400)
      COMMON/RYINT2/ISET(200),KFPR(200,2),COEF(200,20),ICOL(40,4,2)
      COMMON/RYINT5/NGEN(0:200,3),XSEC(0:200,3)
      SAVE /LUDAT1/,/LUDAT2/,/LUDAT3/
      SAVE /RYSUBS/,/RYPARS/,/RYINT1/,/RYINT2/,/RYINT5/
      DIMENSION KFD(2),MBW(2),PMD(2),PGD(2),PMG(2),PML(2),PMU(2),
     &PMH(2),ATL(2),ATU(2),ATH(2),RMG(2),INX1(100),XPT1(100),
     &FPT1(100),INX2(100),XPT2(100),FPT2(100),WDTP(0:40),
     &WDTE(0:40,0:5)

C...Find if particles equal, maximum mass, matrix elements, etc.
      ISUB=MINT(1)
      KFD(1)=IABS(KFD1)
      KFD(2)=IABS(KFD2)
      MEQL=0
      IF(KFD(1).EQ.KFD(2)) MEQL=1
      MLM=0
      IF(MOFSH.GE.2.AND.MEQL.EQ.1) MLM=INT(1.5+PYR(0))
      IF(MOFSH.LE.2) THEN
        NOFF=44
        PMMX=PMMO
      ELSE
        NOFF=40
        PMMX=VINT(1)
        IF(CKIN(2).GT.CKIN(1)) PMMX=MIN(CKIN(2),VINT(1))
      ENDIF
      MMED=0
      IF((KFMO.EQ.25.OR.KFMO.EQ.35.OR.KFMO.EQ.36).AND.MEQL.EQ.1.AND.
     &(KFD(1).EQ.23.OR.KFD(1).EQ.24)) MMED=1
      IF((KFMO.EQ.32.OR.IABS(KFMO).EQ.34).AND.(KFD(1).EQ.23.OR.
     &KFD(1).EQ.24).AND.(KFD(2).EQ.23.OR.KFD(2).EQ.24)) MMED=2
      IF((KFMO.EQ.32.OR.IABS(KFMO).EQ.34).AND.(KFD(2).EQ.25.OR.
     &KFD(2).EQ.35.OR.KFD(2).EQ.36)) MMED=3
      LOOP=1

C...Find where Breit-Wigners are required, else select discrete masses.
  100 DO 110 I=1,2
      KFCA=KFD(I)
      IF(KFCA.GT.100) KFCA=LUCOMP(KFCA)
      PMD(I)=PMAS(KFCA,1)
      PGD(I)=PMAS(KFCA,2)
      IF(MSTP(42).LE.0.OR.PGD(I).LT.PARP(41)) THEN
        MBW(I)=0
        PMG(I)=PMD(I)
      ELSE
        MBW(I)=1
      ENDIF
  110 CONTINUE

C...Find allowed mass range and Breit-Wigner parameters.
      DO 120 I=1,2
      IF(MOFSH.EQ.1.AND.LOOP.EQ.1.AND.MBW(I).EQ.1) THEN
        PML(I)=PARP(42)
        PMU(I)=PMMX-PARP(42)
        IF(MBW(3-I).EQ.0) PMU(I)=MIN(PMU(I),PMMX-PMD(3-I))
        IF(PMU(I).LT.PML(I)+PARJ(64)) MBW(I)=-1
      ELSEIF(MBW(I).EQ.1.OR.MOFSH.GE.5) THEN
        ILM=I
        IF(MLM.EQ.2) ILM=3-I
        PML(I)=MAX(CKIN(NOFF+2*ILM-1),PARP(42))
        IF(MOFSH.GE.5.AND.I.EQ.2) PML(I)=MAX(PML(I),2.*PMAS(KFD2,1))
        PMU(I)=PMMX-MAX(CKIN(NOFF+5-2*ILM),PARP(42))
        IF(MOFSH.GE.5.AND.I.EQ.1) PMU(I)=MIN(PMU(I),PMMX-2.*
     &  PMAS(KFD2,1))
        IF(CKIN(NOFF+2*ILM).GT.CKIN(NOFF+2*ILM-1)) PMU(I)=MIN(PMU(I),
     &  CKIN(NOFF+2*ILM))
        IF(MBW(3-I).EQ.0) PMU(I)=MIN(PMU(I),PMMX-PMD(3-I))
        IF(I.EQ.MLM) PMU(I)=MIN(PMU(I),0.5*PMMX)
        IF(MEQL.EQ.0) PMH(I)=MIN(PMU(I),0.5*PMMX)
        IF(PMU(I).LT.PML(I)+PARJ(64)) MBW(I)=-1
        IF(MBW(I).EQ.1) THEN
          ATL(I)=ATAN((PML(I)**2-PMD(I)**2)/(PMD(I)*PGD(I)))
          ATU(I)=ATAN((PMU(I)**2-PMD(I)**2)/(PMD(I)*PGD(I)))
          IF(MEQL.EQ.0) ATH(I)=ATAN((PMH(I)**2-PMD(I)**2)/(PMD(I)*
     &    PGD(I)))
        ENDIF
      ENDIF
  120 CONTINUE
      IF(MBW(1).LT.0.OR.MBW(2).LT.0.OR.(MBW(1).EQ.0.AND.MBW(2).EQ.0))
     &THEN
        CALL LUERRM(13,'(RYOFSH:) no allowed decay product masses')
        MINT(51)=1
        RETURN
      ENDIF

C...Calculation of partial width of resonance.
      IF(MOFSH.EQ.1) THEN

C..If only one integration, pick that to be the inner.
        IF(MBW(1).EQ.0) THEN
          PM2=PMD(1)
          PMD(1)=PMD(2)
          PGD(1)=PGD(2)
          PML(1)=PML(2)
          PMU(1)=PMU(2)
        ELSEIF(MBW(2).EQ.0) THEN
          PM2=PMD(2)
        ENDIF

C...Start outer loop of integration.
        IF(MBW(1).EQ.1.AND.MBW(2).EQ.1) THEN
          ATL2=ATAN((PML(2)**2-PMD(2)**2)/(PMD(2)*PGD(2)))
          ATU2=ATAN((PMU(2)**2-PMD(2)**2)/(PMD(2)*PGD(2)))
          NPT2=1
          XPT2(1)=1.
          INX2(1)=0
          FMAX2=0.
        ENDIF
  130   IF(MBW(1).EQ.1.AND.MBW(2).EQ.1) THEN
          PM2S=PMD(2)**2+PMD(2)*PGD(2)*TAN(ATL2+XPT2(NPT2)*(ATU2-ATL2))
          PM2=MIN(PMU(2),MAX(PML(2),SQRT(MAX(0.,PM2S))))
        ENDIF
        RM2=(PM2/PMMX)**2

C...Start inner loop of integration.
        PML1=PML(1)
        PMU1=MIN(PMU(1),PMMX-PM2)
        IF(MEQL.EQ.1) PMU1=MIN(PMU1,PM2)
        ATL1=ATAN((PML1**2-PMD(1)**2)/(PMD(1)*PGD(1)))
        ATU1=ATAN((PMU1**2-PMD(1)**2)/(PMD(1)*PGD(1)))
        IF(PML1+PARJ(64).GE.PMU1.OR.ATL1+1E-7.GE.ATU1) THEN
          FUNC2=0.
          GOTO 180
        ENDIF
        NPT1=1
        XPT1(1)=1.
        INX1(1)=0
        FMAX1=0.
  140   PM1S=PMD(1)**2+PMD(1)*PGD(1)*TAN(ATL1+XPT1(NPT1)*(ATU1-ATL1))
        PM1=MIN(PMU1,MAX(PML1,SQRT(MAX(0.,PM1S))))
        RM1=(PM1/PMMX)**2

C...Evaluate function value - inner loop.
        FUNC1=SQRT(MAX(0.,(1.-RM1-RM2)**2-4.*RM1*RM2))
        IF(MMED.EQ.1) FUNC1=FUNC1*((1.-RM1-RM2)**2+8.*RM1*RM2)
        IF(MMED.EQ.2) FUNC1=FUNC1**3*(1.+10.*RM1+10.*RM2+RM1**2+
     &  RM2**2+10.*RM1*RM2)
        IF(FUNC1.GT.FMAX1) FMAX1=FUNC1
        FPT1(NPT1)=FUNC1

C...Go to next position in inner loop.
        IF(NPT1.EQ.1) THEN
          NPT1=NPT1+1
          XPT1(NPT1)=0.
          INX1(NPT1)=1
          GOTO 140
        ELSEIF(NPT1.LE.8) THEN
          NPT1=NPT1+1
          IF(NPT1.LE.4.OR.NPT1.EQ.6) ISH1=1
          ISH1=ISH1+1
          XPT1(NPT1)=0.5*(XPT1(ISH1)+XPT1(INX1(ISH1)))
          INX1(NPT1)=INX1(ISH1)
          INX1(ISH1)=NPT1
          GOTO 140
        ELSEIF(NPT1.LT.100) THEN
          ISN1=ISH1
  150     ISH1=ISH1+1
          IF(ISH1.GT.NPT1) ISH1=2
          IF(ISH1.EQ.ISN1) GOTO 160
          DFPT1=ABS(FPT1(ISH1)-FPT1(INX1(ISH1)))
          IF(DFPT1.LT.PARP(43)*FMAX1) GOTO 150
          NPT1=NPT1+1
          XPT1(NPT1)=0.5*(XPT1(ISH1)+XPT1(INX1(ISH1)))
          INX1(NPT1)=INX1(ISH1)
          INX1(ISH1)=NPT1
          GOTO 140
        ENDIF

C...Calculate integral over inner loop.
  160   FSUM1=0.
        DO 170 IPT1=2,NPT1
  170   FSUM1=FSUM1+0.5*(FPT1(IPT1)+FPT1(INX1(IPT1)))*
     &  (XPT1(INX1(IPT1))-XPT1(IPT1))
        FUNC2=FSUM1*(ATU1-ATL1)/PARU(1)
  180   IF(MBW(1).EQ.1.AND.MBW(2).EQ.1) THEN
          IF(FUNC2.GT.FMAX2) FMAX2=FUNC2
          FPT2(NPT2)=FUNC2

C...Go to next position in outer loop.
          IF(NPT2.EQ.1) THEN
            NPT2=NPT2+1
            XPT2(NPT2)=0.
            INX2(NPT2)=1
            GOTO 130
          ELSEIF(NPT2.LE.8) THEN
            NPT2=NPT2+1
            IF(NPT2.LE.4.OR.NPT2.EQ.6) ISH2=1
            ISH2=ISH2+1
            XPT2(NPT2)=0.5*(XPT2(ISH2)+XPT2(INX2(ISH2)))
            INX2(NPT2)=INX2(ISH2)
            INX2(ISH2)=NPT2
            GOTO 130
          ELSEIF(NPT2.LT.100) THEN
            ISN2=ISH2
  190       ISH2=ISH2+1
            IF(ISH2.GT.NPT2) ISH2=2
            IF(ISH2.EQ.ISN2) GOTO 200
            DFPT2=ABS(FPT2(ISH2)-FPT2(INX2(ISH2)))
            IF(DFPT2.LT.PARP(43)*FMAX2) GOTO 190
            NPT2=NPT2+1
            XPT2(NPT2)=0.5*(XPT2(ISH2)+XPT2(INX2(ISH2)))
            INX2(NPT2)=INX2(ISH2)
            INX2(ISH2)=NPT2
            GOTO 130
          ENDIF

C...Calculate integral over outer loop.
  200     FSUM2=0.
          DO 210 IPT2=2,NPT2
  210     FSUM2=FSUM2+0.5*(FPT2(IPT2)+FPT2(INX2(IPT2)))*
     &    (XPT2(INX2(IPT2))-XPT2(IPT2))
          FSUM2=FSUM2*(ATU2-ATL2)/PARU(1)
          IF(MEQL.EQ.1) FSUM2=2.*FSUM2
        ELSE
          FSUM2=FUNC2
        ENDIF

C...Save result; second integration for user-selected mass range.
        IF(LOOP.EQ.1) WIDW=FSUM2
        WID2=FSUM2
        IF(LOOP.EQ.1.AND.(CKIN(46).GE.CKIN(45).OR.CKIN(48).GE.CKIN(47).
     &  OR.MAX(CKIN(45),CKIN(47)).GE.1.01*PARP(42))) THEN
          LOOP=2
          GOTO 100
        ENDIF
        RET1=WIDW
        RET2=WID2/WIDW

C...Select two decay product masses of a resonance.
      ELSEIF(MOFSH.EQ.2) THEN
  220   DO 230 I=1,2
        IF(MBW(I).EQ.0) GOTO 230
        PMBW=PMD(I)**2+PMD(I)*PGD(I)*TAN(ATL(I)+PYR(0)*(ATU(I)-ATL(I)))
        PMG(I)=MIN(PMU(I),MAX(PML(I),SQRT(MAX(0.,PMBW))))
        RMG(I)=(PMG(I)/PMMX)**2
  230   CONTINUE
        IF((MEQL.EQ.1.AND.PMG(MAX(1,MLM)).GT.PMG(MIN(2,3-MLM))).OR.
     &  PMG(1)+PMG(2)+PARJ(64).GT.PMMX) GOTO 220

C...Weight with matrix element (if none known, use beta factor).
        FLAM=SQRT(MAX(0.,(1.-RMG(1)-RMG(2))**2-4.*RMG(1)*RMG(2)))
        IF(MMED.EQ.1) THEN
          WTBE=FLAM*((1.-RMG(1)-RMG(2))**2+8.*RMG(1)*RMG(2))
        ELSEIF(MMED.EQ.2) THEN
          WTBE=FLAM**3*(1.+10.*RMG(1)+10.*RMG(2)+RMG(1)**2+
     &    RMG(2)**2+10.*RMG(1)*RMG(2))
        ELSEIF(MMED.EQ.3) THEN
          WTBE=FLAM*(RMG(1)+FLAM**2/12.)
        ELSE
          WTBE=FLAM
        ENDIF
        IF(WTBE.LT.PYR(0)) GOTO 220
        RET1=PMG(1)
        RET2=PMG(2)

C...Find suitable set of masses for initialization of 2 -> 2 processes.
      ELSEIF(MOFSH.EQ.3) THEN
        IF(MBW(1).NE.0.AND.MBW(2).EQ.0) THEN
          PMG(1)=MIN(PMD(1),0.5*(PML(1)+PMU(1)))
          PMG(2)=PMD(2)
        ELSEIF(MBW(2).NE.0.AND.MBW(1).EQ.0) THEN
          PMG(1)=PMD(1)
          PMG(2)=MIN(PMD(2),0.5*(PML(2)+PMU(2)))
        ELSE
          IDIV=-1
  240     IDIV=IDIV+1
          PMG(1)=MIN(PMD(1),0.1*(IDIV*PML(1)+(10-IDIV)*PMU(1)))
          PMG(2)=MIN(PMD(2),0.1*(IDIV*PML(2)+(10-IDIV)*PMU(2)))
          IF(IDIV.LE.9.AND.PMG(1)+PMG(2).GT.0.9*PMMX) GOTO 240
        ENDIF
        RET1=PMG(1)
        RET2=PMG(2)

C...Evaluate importance of excluded tails of Breit-Wigners.
        IF(MEQL.EQ.0.AND.MBW(1).EQ.1.AND.MBW(2).EQ.1.AND.PMD(1)+PMD(2).
     &  GT.PMMX.AND.PMH(1).GT.PML(1).AND.PMH(2).GT.PML(2)) MEQL=2
        IF(MEQL.LE.1) THEN
          VINT(80)=1.
          DO 250 I=1,2
  250     IF(MBW(I).NE.0) VINT(80)=VINT(80)*1.25*(ATU(I)-ATL(I))/PARU(1)
        ELSE
          VINT(80)=(1.25/PARU(1))**2*MAX((ATU(1)-ATL(1))*
     &    (ATH(2)-ATL(2)),(ATH(1)-ATL(1))*(ATU(2)-ATL(2)))
        ENDIF
        IF(ISUB.EQ.22.AND.MSTP(43).NE.2) VINT(80)=4.*VINT(80)
        IF(MEQL.GE.1) VINT(80)=2.*VINT(80)

C...Pick one particle to be the lighter (if improves efficiency).
      ELSEIF(MOFSH.EQ.4) THEN
        IF(MEQL.EQ.0.AND.MBW(1).EQ.1.AND.MBW(2).EQ.1.AND.PMD(1)+PMD(2).
     &  GT.PMMX.AND.PMH(1).GT.PML(1).AND.PMH(2).GT.PML(2)) MEQL=2
  260   IF(MEQL.EQ.2) MLM=INT(1.5+PYR(0))

C...Select two masses according to Breit-Wigner + flat in s + 1/s.
        DO 270 I=1,2
        IF(MBW(I).EQ.0) GOTO 270
        PMV=PMU(I)
        IF(MEQL.EQ.2.AND.I.EQ.MLM) PMV=PMH(I)
        ATV=ATU(I)
        IF(MEQL.EQ.2.AND.I.EQ.MLM) ATV=ATH(I)
        RBR=PYR(0)
        IF(ISUB.EQ.22.AND.MSTP(43).NE.2) RBR=2.*RBR
        IF(RBR.LT.0.8) THEN
          PMSR=PMD(I)**2+PMD(I)*PGD(I)*TAN(ATL(I)+PYR(0)*(ATV-ATL(I)))
          PMG(I)=MIN(PMV,MAX(PML(I),SQRT(MAX(0.,PMSR))))
        ELSEIF(RBR.LT.0.9) THEN
          PMG(I)=SQRT(MAX(0.,PML(I)**2+PYR(0)*(PMV**2-PML(I)**2)))
        ELSEIF(RBR.LT.1.5) THEN
          PMG(I)=PML(I)*(PMV/PML(I))**PYR(0)
        ELSE
          PMG(I)=SQRT(MAX(0.,PML(I)**2*PMV**2/(PML(I)**2+PYR(0)*
     &    (PMV**2-PML(I)**2))))
        ENDIF
  270   CONTINUE
        IF((MEQL.GE.1.AND.PMG(MAX(1,MLM)).GT.PMG(MIN(2,3-MLM))).OR.
     &  PMG(1)+PMG(2)+PARJ(64).GT.PMMX) THEN
          NGEN(0,1)=NGEN(0,1)+1
          NGEN(MINT(1),1)=NGEN(MINT(1),1)+1
          GOTO 260
        ENDIF
        RET1=PMG(1)
        RET2=PMG(2)

C...Give weight for selected mass distribution.
        VINT(80)=1.
        DO 280 I=1,2
        IF(MBW(I).EQ.0) GOTO 280
        PMV=PMU(I)
        IF(MEQL.EQ.2.AND.I.EQ.MLM) PMV=PMH(I)
        ATV=ATU(I)
        IF(MEQL.EQ.2.AND.I.EQ.MLM) ATV=ATH(I)
        F0=PMD(I)*PGD(I)/((PMG(I)**2-PMD(I)**2)**2+
     &  (PMD(I)*PGD(I))**2)/PARU(1)
        F1=1.
        F2=1./PMG(I)**2
        F3=1./PMG(I)**4
        FI0=(ATV-ATL(I))/PARU(1)
        FI1=PMV**2-PML(I)**2
        FI2=2.*LOG(PMV/PML(I))
        FI3=1./PML(I)**2-1./PMV**2
        IF(ISUB.EQ.22.AND.MSTP(43).NE.2) THEN
          VINT(80)=VINT(80)*20./(8.+(FI0/F0)*(F1/FI1+6.*F2/FI2+
     &    5.*F3/FI3))
        ELSE
          VINT(80)=VINT(80)*10./(8.+(FI0/F0)*(F1/FI1+F2/FI2))
        ENDIF
        VINT(80)=VINT(80)*FI0
  280   CONTINUE
        IF(MEQL.GE.1) VINT(80)=2.*VINT(80)

      ELSEIF(MOFSH.EQ.5) THEN
C...Find suitable set of masses for initialization of 2 -> 3 process.
        IDIV=6
  290   IDIV=IDIV-1
        IF(MBW(1).EQ.0) THEN
          PMG(1)=PMD(1)
        ELSE
          PMSR=PMD(1)**2+PMD(1)*PGD(1)*TAN(ATL(1)+0.1*IDIV*(ATU(1)-
     &    ATL(1)))
          PMG(1)=MIN(PMU(1),MAX(PML(1),SQRT(MAX(0.,PMSR))))
        ENDIF
        PMG(2)=PML(2)*(PMU(2)/PML(2))**(0.1*IDIV)
        IF(IDIV.GE.1.AND.PMG(1)+PMG(2).GT.0.9*PMMX) GOTO 290
        RET1=PMG(1)
        RET2=PMG(2)

C...Evaluate size of selected phase space volume.
        VINT(80)=2.*LOG(PMU(2)/PML(2))
        IF(MBW(1).NE.0) VINT(80)=VINT(80)*1.25*(ATU(1)-ATL(1))/PARU(1)

C...Pick decay angles.
        VINT(81)=0.
        VINT(82)=0.5*PARU(1)
        VINT(83)=1.
        VINT(84)=0.

C...Select flavour of resonance decays.
        KFA=KFPR(ISUB,1)
        CALL RYWIDT(KFA,PMG(1)**2,WDTP,WDTE)
        IF(KCHG(KFA,3).EQ.0) THEN
          IPM=2
        ELSE
          IPM=(5+ISIGN(1,KFA))/2
        ENDIF
        WDTE0S=WDTE(0,1)+WDTE(0,IPM)+WDTE(0,4)
        IF(WDTE0S.LE.0.) THEN
          CALL LUERRM(12,'(RYOFSH:) no allowed resonace decay channel')
          MINT(51)=1
          RETURN
        ENDIF
        WDTEC=0.
        DO 300 IDL=1,MDCY(KFA,3)
        WDTEK=WDTE(IDL,1)+WDTE(IDL,IPM)+WDTE(IDL,4)
        IF(WDTEK.GT.WDTEC) THEN
          IDC=IDL+MDCY(KFA,2)-1
          WDTEC=WDTEK
        ENDIF
  300   CONTINUE
        MINT(35)=IDC

C...Compensating factor for all flavours.
        KFL=IABS(KFDP(IDC,1))
        QFL=KCHG(KFL,1)/3.
        AFL=SIGN(1.,QFL+0.1)
        VFL=AFL-4.*PARU(102)*QFL
        WDTEK=VFL**2+AFL**2
        VINT(80)=VINT(80)*WDTE0S/WDTEK

      ELSEIF(MOFSH.EQ.6) THEN
C...Select two masses, one basically Breit-Wigner, other dm^2/m^2.
        IF(MBW(1).NE.0) THEN
          RBR=PYR(0)
          IF(RBR.LT.0.8) THEN
            PMSR=PMD(1)**2+PMD(1)*PGD(1)*TAN(ATL(1)+PYR(0)*
     &      (ATU(1)-ATL(1)))
            PMG(1)=MIN(PMU(1),MAX(PML(1),SQRT(MAX(0.,PMSR))))
          ELSEIF(RBR.LT.0.9) THEN
            PMG(1)=SQRT(MAX(0.,PML(1)**2+PYR(0)*(PMU(1)**2-PML(1)**2)))
          ELSE
            PMG(1)=PML(1)*(PMU(1)/PML(1))**PYR(0)
          ENDIF
        ENDIF
        PMG(2)=PML(2)*(PMU(2)/PML(2))**PYR(0)
        IF(SQRT(MAX(0.,1.-(PML(2)/PMG(2))**2)).LT.PYR(0).OR.
     &  PMG(1)+PMG(2)+PARJ(64).GT.PMMX) THEN
          MINT(51)=1
          RETURN
        ENDIF
        RET1=PMG(1)
        RET2=PMG(2)

C...Give weight for selected mass distribution.
        VINT(80)=2.*LOG(PMU(2)/PML(2))
        IF(MBW(1).NE.0) THEN
          F0=PMD(1)*PGD(1)/((PMG(1)**2-PMD(1)**2)**2+
     &    (PMD(1)*PGD(1))**2)/PARU(1)
          F1=1.
          F2=1./PMG(1)**2
          FI0=(ATU(1)-ATL(1))/PARU(1)
          FI1=PMU(1)**2-PML(1)**2
          FI2=2.*LOG(PMU(1)/PML(1))
          VINT(80)=VINT(80)*10.*FI0/(8.+(FI0/F0)*(F1/FI1+F2/FI2))
        ENDIF

C...Select decay angles.
        VINT(81)=2.*PYR(0)-1.
        VINT(82)=PARU(2)*PYR(0)
        VINT(83)=2.*PYR(0)-1.
        VINT(84)=PARU(2)*PYR(0)

C...Select flavour of resonance decays.
        KFA=KFPR(ISUB,1)
        CALL RYWIDT(KFA,PMG(1)**2,WDTP,WDTE)
        IF(KCHG(KFA,3).EQ.0) THEN
          IPM=2
        ELSE
          IPM=(5+ISIGN(1,KFA))/2
        ENDIF
        WDTE0S=WDTE(0,1)+WDTE(0,IPM)+WDTE(0,4)
        IF(WDTE0S.LE.0.) THEN
          CALL LUERRM(12,'(RYOFSH:) no allowed resonace decay channel')
          MINT(51)=1
          RETURN
        ENDIF
        RKFL=WDTE0S*PYR(0)
        IDL=0
  320   IDL=IDL+1
        IDC=IDL+MDCY(KFA,2)-1
        RKFL=RKFL-(WDTE(IDL,1)+WDTE(IDL,IPM)+WDTE(IDL,4))
        IF(IDL.LT.MDCY(KFA,3).AND.RKFL.GT.0.) GOTO 320
        MINT(35)=IDC

C...Compensating factor for all flavours.
        KFL=IABS(KFDP(IDC,1))
        QFL=KCHG(KFL,1)/3.
        AFL=SIGN(1.,QFL+0.1)
        VFL=AFL-4.*PARU(102)*QFL
        WDTEK=VFL**2+AFL**2
        VINT(80)=VINT(80)*WDTE0S/WDTEK
      ENDIF

      RETURN
      END

C***********************************************************************

      SUBROUTINE RYKLIM(ILIM)
      PARAMETER (KSZJ=4000)

C...Checks generated variables against pre-set kinematical limits;
C...also calculates limits on variables used in generation.
      
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      integer N,NPAD,K
      double precision P,V
      SAVE /PYJETS/


      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/LUDAT3/MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5)
      COMMON/RYSUBS/MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200)
      COMMON/RYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/RYINT1/MINT(400),VINT(400)
      COMMON/RYINT2/ISET(200),KFPR(200,2),COEF(200,20),ICOL(40,4,2)
      SAVE /LUDAT1/,/LUDAT2/,/LUDAT3/
      SAVE /RYSUBS/,/RYPARS/,/RYINT1/,/RYINT2/

C...Common kinematical expressions.
      MINT(51)=0
      ISUB=MINT(1)
      ISTSB=ISET(ISUB)
      IF(ISUB.EQ.96) GOTO 120
      SQM3=VINT(63)
      SQM4=VINT(64)
      IF(ILIM.NE.1) THEN
        TAU=VINT(21)
        RM3=SQM3/(TAU*VINT(2))
        RM4=SQM4/(TAU*VINT(2))
        BE34=SQRT(MAX(1E-20,(1.-RM3-RM4)**2-4.*RM3*RM4))
      ENDIF
      PTHMIN=CKIN(3)
      IF(MIN(SQM3,SQM4).LT.CKIN(6)**2) PTHMIN=MAX(CKIN(3),CKIN(5))

      IF(ILIM.EQ.0) THEN
C...Check generated values of tau, y*, cos(theta-hat), and tau' against
C...pre-set kinematical limits.
        YST=VINT(22)
        CTH=VINT(23)
        TAUP=VINT(26)
        TAUE=TAU
        IF(ISTSB.GE.3.AND.ISTSB.LE.5) TAUE=TAUP
        X1=SQRT(TAUE)*EXP(YST)
        X2=SQRT(TAUE)*EXP(-YST)
        XF=X1-X2
        IF(TAU*VINT(2).LT.CKIN(1)**2) MINT(51)=1
        IF(CKIN(2).GE.0..AND.TAU*VINT(2).GT.CKIN(2)**2) MINT(51)=1
        IF(X1.LT.CKIN(21).OR.X1.GT.CKIN(22)) MINT(51)=1
        IF(X2.LT.CKIN(23).OR.X2.GT.CKIN(24)) MINT(51)=1
        IF(XF.LT.CKIN(25).OR.XF.GT.CKIN(26)) MINT(51)=1
        IF(YST.LT.CKIN(7).OR.YST.GT.CKIN(8)) MINT(51)=1
        IF(ISTSB.EQ.2.OR.ISTSB.EQ.4.OR.ISTSB.EQ.6) THEN
          PTH=0.5*BE34*SQRT(TAU*VINT(2)*MAX(0.,1.-CTH**2))
          EXPY3=MAX(1.E-10,(1.+RM3-RM4+BE34*CTH)/
     &    MAX(1.E-10,(1.+RM3-RM4-BE34*CTH)))
          EXPY4=MAX(1.E-10,(1.-RM3+RM4-BE34*CTH)/
     &    MAX(1.E-10,(1.-RM3+RM4+BE34*CTH)))
          Y3=YST+0.5*LOG(EXPY3)
          Y4=YST+0.5*LOG(EXPY4)
          YLARGE=MAX(Y3,Y4)
          YSMALL=MIN(Y3,Y4)
          ETALAR=10.
          ETASMA=-10.
          STH=SQRT(MAX(0.,1.-CTH**2))
          EXSQ3=SQRT(MAX(1E-20,((1.+RM3-RM4)*COSH(YST)+BE34*SINH(YST)*
     &    CTH)**2-4.*RM3))
          EXSQ4=SQRT(MAX(1E-20,((1.-RM3+RM4)*COSH(YST)-BE34*SINH(YST)*
     &    CTH)**2-4.*RM4))
          IF(STH.LT.1.E-6) GOTO 100
          EXPET3=((1.+RM3-RM4)*SINH(YST)+BE34*COSH(YST)*CTH+EXSQ3)/
     &    (BE34*STH)
          EXPET4=((1.-RM3+RM4)*SINH(YST)-BE34*COSH(YST)*CTH+EXSQ4)/
     &    (BE34*STH)
          ETA3=LOG(MIN(1.E10,MAX(1.E-10,EXPET3)))
          ETA4=LOG(MIN(1.E10,MAX(1.E-10,EXPET4)))
          ETALAR=MAX(ETA3,ETA4)
          ETASMA=MIN(ETA3,ETA4)
  100     CTS3=((1.+RM3-RM4)*SINH(YST)+BE34*COSH(YST)*CTH)/EXSQ3
          CTS4=((1.-RM3+RM4)*SINH(YST)-BE34*COSH(YST)*CTH)/EXSQ4
          CTSLAR=MIN(1.,MAX(CTS3,CTS4))
          CTSSMA=MAX(-1.,MIN(CTS3,CTS4))
          IF(PTH.LT.PTHMIN) MINT(51)=1
          IF(CKIN(4).GE.0..AND.PTH.GT.CKIN(4)) MINT(51)=1
          IF(YLARGE.LT.CKIN(9).OR.YLARGE.GT.CKIN(10)) MINT(51)=1
          IF(YSMALL.LT.CKIN(11).OR.YSMALL.GT.CKIN(12)) MINT(51)=1
          IF(ETALAR.LT.CKIN(13).OR.ETALAR.GT.CKIN(14)) MINT(51)=1
          IF(ETASMA.LT.CKIN(15).OR.ETASMA.GT.CKIN(16)) MINT(51)=1
          IF(CTSLAR.LT.CKIN(17).OR.CTSLAR.GT.CKIN(18)) MINT(51)=1
          IF(CTSSMA.LT.CKIN(19).OR.CTSSMA.GT.CKIN(20)) MINT(51)=1
          IF(CTH.LT.CKIN(27).OR.CTH.GT.CKIN(28)) MINT(51)=1
        ENDIF
        IF(ISTSB.GE.3.AND.ISTSB.LE.5) THEN
          IF(TAUP*VINT(2).LT.CKIN(31)**2) MINT(51)=1
          IF(CKIN(32).GE.0..AND.TAUP*VINT(2).GT.CKIN(32)**2) MINT(51)=1
        ENDIF

C...Additional p_T cuts on 2 -> 3 process.
        IF(ISTSB.EQ.6) THEN
          KFQ=KFPR(131,2)
          PMQQ=SQRT(VINT(64))
          PMQ=PMAS(KFQ,1)
          PZQ=SQRT(MAX(0.,(0.5*PMQQ)**2-PMQ**2))
          DO 110 I=MINT(84)+1,MINT(84)+2
          K(I,1)=1
          P(I,1)=0.
          P(I,2)=0.
          P(I,3)=PZQ*(-1.)**(I-1)
          P(I,4)=0.5*PMQQ
  110     P(I,5)=PMQ
          PEQQ=0.5*SQRT(TAU*VINT(2))*(1.+(VINT(64)-VINT(63))/
     &    (TAU*VINT(2)))
          PZQQ=SQRT(MAX(0.,PEQQ**2-VINT(64)))
          CALL LUDBRB(MINT(84)+1,MINT(84)+2,ACOS(VINT(83)),VINT(84),
     &    0D0,0D0,-DBLE(PZQQ/PEQQ))
          CALL LUDBRB(MINT(84)+1,MINT(84)+2,ACOS(VINT(23)),VINT(24),
     &    0D0,0D0,0D0)
          PTQ2=SQRT(P(MINT(84)+1,1)**2+P(MINT(84)+1,2)**2)
          PTQ3=SQRT(P(MINT(84)+2,1)**2+P(MINT(84)+2,2)**2)
          PTMNQ=MIN(PTQ2,PTQ3)
          PTMXQ=MAX(PTQ2,PTQ3)
          IF(PTMNQ.LT.CKIN(51)) MINT(51)=1
          IF(CKIN(52).GE.0..AND.PTMNQ.GT.CKIN(52)) MINT(51)=1
          IF(PTMXQ.LT.CKIN(53)) MINT(51)=1
          IF(CKIN(54).GE.0..AND.PTMXQ.GT.CKIN(54)) MINT(51)=1
          VINT(85)=PTMNQ
          VINT(86)=PTMXQ
        ENDIF

      ELSEIF(ILIM.EQ.1) THEN
C...Calculate limits on tau
C...0) due to definition
        TAUMN0=0.
        TAUMX0=1.
C...1) due to limits on subsystem mass
        TAUMN1=CKIN(1)**2/VINT(2)
        TAUMX1=1.
        IF(CKIN(2).GE.0.) TAUMX1=CKIN(2)**2/VINT(2)
C...2) due to limits on pT-hat (and non-overlapping rapidity intervals)
        TM3=SQRT(SQM3+PTHMIN**2)
        TM4=SQRT(SQM4+PTHMIN**2)
        YDCOSH=1.
        IF(CKIN(9).GT.CKIN(12)) YDCOSH=COSH(CKIN(9)-CKIN(12))
        TAUMN2=(TM3**2+2.*TM3*TM4*YDCOSH+TM4**2)/VINT(2)
        TAUMX2=1.
C...3) due to limits on pT-hat and cos(theta-hat)
        CTH2MN=MIN(CKIN(27)**2,CKIN(28)**2)
        CTH2MX=MAX(CKIN(27)**2,CKIN(28)**2)
        TAUMN3=0.
        IF(CKIN(27)*CKIN(28).GT.0.) TAUMN3=
     &  (SQRT(SQM3+PTHMIN**2/(1.-CTH2MN))+
     &  SQRT(SQM4+PTHMIN**2/(1.-CTH2MN)))**2/VINT(2)
        TAUMX3=1.
        IF(CKIN(4).GE.0..AND.CTH2MX.LT.1.) TAUMX3=
     &  (SQRT(SQM3+CKIN(4)**2/(1.-CTH2MX))+
     &  SQRT(SQM4+CKIN(4)**2/(1.-CTH2MX)))**2/VINT(2)
C...4) due to limits on x1 and x2
        TAUMN4=CKIN(21)*CKIN(23)
        TAUMX4=CKIN(22)*CKIN(24)
C...5) due to limits on xF
        TAUMN5=0.
        TAUMX5=MAX(1.-CKIN(25),1.+CKIN(26))
        VINT(11)=MAX(TAUMN0,TAUMN1,TAUMN2,TAUMN3,TAUMN4,TAUMN5)
        VINT(31)=MIN(TAUMX0,TAUMX1,TAUMX2,TAUMX3,TAUMX4,TAUMX5)
        IF(MINT(47).EQ.1.AND.(ISTSB.EQ.1.OR.ISTSB.EQ.2.OR.ISTSB.EQ.6))
     &  THEN
          VINT(11)=0.99999
          VINT(31)=1.00001
        ENDIF
        IF(VINT(31).LE.VINT(11)) MINT(51)=1

      ELSEIF(ILIM.EQ.2) THEN
C...Calculate limits on y*
        TAUE=TAU
        IF(ISTSB.GE.3.AND.ISTSB.LE.5) TAUE=VINT(26)
        TAURT=SQRT(TAUE)
C...0) due to kinematics
        YSTMN0=LOG(TAURT)
        YSTMX0=-YSTMN0
C...1) due to explicit limits
        YSTMN1=CKIN(7)
        YSTMX1=CKIN(8)
C...2) due to limits on x1
        YSTMN2=LOG(MAX(TAUE,CKIN(21))/TAURT)
        YSTMX2=LOG(MAX(TAUE,CKIN(22))/TAURT)
C...3) due to limits on x2
        YSTMN3=-LOG(MAX(TAUE,CKIN(24))/TAURT)
        YSTMX3=-LOG(MAX(TAUE,CKIN(23))/TAURT)
C...4) due to limits on xF
        YEPMN4=0.5*ABS(CKIN(25))/TAURT
        YSTMN4=SIGN(LOG(MAX(1E-20,SQRT(1.+YEPMN4**2)+YEPMN4)),CKIN(25))
        YEPMX4=0.5*ABS(CKIN(26))/TAURT
        YSTMX4=SIGN(LOG(MAX(1E-20,SQRT(1.+YEPMX4**2)+YEPMX4)),CKIN(26))
C...5) due to simultaneous limits on y-large and y-small
        YEPSMN=(RM3-RM4)*SINH(CKIN(9)-CKIN(11))
        YEPSMX=(RM3-RM4)*SINH(CKIN(10)-CKIN(12))
        YDIFMN=ABS(LOG(MAX(1E-20,SQRT(1.+YEPSMN**2)-YEPSMN)))
        YDIFMX=ABS(LOG(MAX(1E-20,SQRT(1.+YEPSMX**2)-YEPSMX)))
        YSTMN5=0.5*(CKIN(9)+CKIN(11)-YDIFMN)
        YSTMX5=0.5*(CKIN(10)+CKIN(12)+YDIFMX)
C...6) due to simultaneous limits on cos(theta-hat) and y-large or
C...   y-small
        CTHLIM=SQRT(MAX(0.,1.-4.*PTHMIN**2/(BE34**2*TAUE*VINT(2))))
        RZMN=BE34*MAX(CKIN(27),-CTHLIM)
        RZMX=BE34*MIN(CKIN(28),CTHLIM)
        YEX3MX=(1.+RM3-RM4+RZMX)/MAX(1E-10,1.+RM3-RM4-RZMX)
        YEX4MX=(1.+RM4-RM3-RZMN)/MAX(1E-10,1.+RM4-RM3+RZMN)
        YEX3MN=MAX(1E-10,1.+RM3-RM4+RZMN)/(1.+RM3-RM4-RZMN)
        YEX4MN=MAX(1E-10,1.+RM4-RM3-RZMX)/(1.+RM4-RM3+RZMX)
        YSTMN6=CKIN(9)-0.5*LOG(MAX(YEX3MX,YEX4MX))
        YSTMX6=CKIN(12)-0.5*LOG(MIN(YEX3MN,YEX4MN))
        VINT(12)=MAX(YSTMN0,YSTMN1,YSTMN2,YSTMN3,YSTMN4,YSTMN5,YSTMN6)
        VINT(32)=MIN(YSTMX0,YSTMX1,YSTMX2,YSTMX3,YSTMX4,YSTMX5,YSTMX6)
        IF(MINT(47).EQ.1) THEN
          VINT(12)=-0.00001
          VINT(32)=0.00001
        ELSEIF(MINT(47).EQ.2) THEN
          VINT(12)=0.99999*YSTMX0
          VINT(32)=1.00001*YSTMX0
        ELSEIF(MINT(47).EQ.3) THEN
          VINT(12)=-1.00001*YSTMX0
          VINT(32)=-0.99999*YSTMX0
        ENDIF
        IF(VINT(32).LE.VINT(12)) MINT(51)=1

      ELSEIF(ILIM.EQ.3) THEN
C...Calculate limits on cos(theta-hat)
        YST=VINT(22)
C...0) due to definition
        CTNMN0=-1.
        CTNMX0=0.
        CTPMN0=0.
        CTPMX0=1.
C...1) due to explicit limits
        CTNMN1=MIN(0.,CKIN(27))
        CTNMX1=MIN(0.,CKIN(28))
        CTPMN1=MAX(0.,CKIN(27))
        CTPMX1=MAX(0.,CKIN(28))
C...2) due to limits on pT-hat
        CTNMN2=-SQRT(MAX(0.,1.-4.*PTHMIN**2/(BE34**2*TAU*VINT(2))))
        CTPMX2=-CTNMN2
        CTNMX2=0.
        CTPMN2=0.
        IF(CKIN(4).GE.0.) THEN
          CTNMX2=-SQRT(MAX(0.,1.-4.*CKIN(4)**2/(BE34**2*TAU*VINT(2))))
          CTPMN2=-CTNMX2
        ENDIF
C...3) due to limits on y-large and y-small
        CTNMN3=MIN(0.,MAX((1.+RM3-RM4)/BE34*TANH(CKIN(11)-YST),
     &  -(1.-RM3+RM4)/BE34*TANH(CKIN(10)-YST)))
        CTNMX3=MIN(0.,(1.+RM3-RM4)/BE34*TANH(CKIN(12)-YST),
     &  -(1.-RM3+RM4)/BE34*TANH(CKIN(9)-YST))
        CTPMN3=MAX(0.,(1.+RM3-RM4)/BE34*TANH(CKIN(9)-YST),
     &  -(1.-RM3+RM4)/BE34*TANH(CKIN(12)-YST))
        CTPMX3=MAX(0.,MIN((1.+RM3-RM4)/BE34*TANH(CKIN(10)-YST),
     &  -(1.-RM3+RM4)/BE34*TANH(CKIN(11)-YST)))
        VINT(13)=MAX(CTNMN0,CTNMN1,CTNMN2,CTNMN3)
        VINT(33)=MIN(CTNMX0,CTNMX1,CTNMX2,CTNMX3)
        VINT(14)=MAX(CTPMN0,CTPMN1,CTPMN2,CTPMN3)
        VINT(34)=MIN(CTPMX0,CTPMX1,CTPMX2,CTPMX3)
        IF(VINT(33).LE.VINT(13).AND.VINT(34).LE.VINT(14)) MINT(51)=1

      ELSEIF(ILIM.EQ.4) THEN
C...Calculate limits on tau'
C...0) due to kinematics
        TAPMN0=TAU
        TAPMX0=1.
C...1) due to explicit limits
        TAPMN1=CKIN(31)**2/VINT(2)
        TAPMX1=1.
        IF(CKIN(32).GE.0.) TAPMX1=CKIN(32)**2/VINT(2)
        VINT(16)=MAX(TAPMN0,TAPMN1)
        VINT(36)=MIN(TAPMX0,TAPMX1)
        IF(MINT(47).EQ.1) THEN
          VINT(16)=0.99999
          VINT(36)=1.00001
        ENDIF
        IF(VINT(36).LE.VINT(16)) MINT(51)=1

      ENDIF
      RETURN

C...Special case for low-pT and multiple interactions:
C...effective kinematical limits for tau, y*, cos(theta-hat).
  120 IF(ILIM.EQ.0) THEN
      ELSEIF(ILIM.EQ.1) THEN
        IF(MSTP(82).LE.1) VINT(11)=4.*PARP(81)**2/VINT(2)
        IF(MSTP(82).GE.2) VINT(11)=PARP(82)**2/VINT(2)
        VINT(31)=1.
      ELSEIF(ILIM.EQ.2) THEN
        VINT(12)=0.5*LOG(VINT(21))
        VINT(32)=-VINT(12)
      ELSEIF(ILIM.EQ.3) THEN
        IF(MSTP(82).LE.1) ST2EFF=4.*PARP(81)**2/(VINT(21)*VINT(2))
        IF(MSTP(82).GE.2) ST2EFF=0.01*PARP(82)**2/(VINT(21)*VINT(2))
        VINT(13)=-SQRT(MAX(0.,1.-ST2EFF))
        VINT(33)=0.
        VINT(14)=0.
        VINT(34)=-VINT(13)
      ENDIF

      RETURN
      END

C*********************************************************************

      SUBROUTINE RYKMAP(IVAR,MVAR,VVAR)

C...Maps a uniform distribution into a distribution of a kinematical
C...variable according to one of the possibilities allowed. It is
C...assumed that kinematical limits have been set by a RYKLIM call.
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/RYSUBS/MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200)
      COMMON/RYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/RYINT1/MINT(400),VINT(400)
      COMMON/RYINT2/ISET(200),KFPR(200,2),COEF(200,20),ICOL(40,4,2)
      SAVE /LUDAT1/,/LUDAT2/
      SAVE /RYSUBS/,/RYPARS/,/RYINT1/,/RYINT2/

C...Convert VVAR to tau variable.
      ISUB=MINT(1)
      ISTSB=ISET(ISUB)
      IF(IVAR.EQ.1) THEN
        TAUMIN=VINT(11)
        TAUMAX=VINT(31)
        IF(MVAR.EQ.3.OR.MVAR.EQ.4) THEN
          TAURE=VINT(73)
          GAMRE=VINT(74)
        ELSEIF(MVAR.EQ.5.OR.MVAR.EQ.6) THEN
          TAURE=VINT(75)
          GAMRE=VINT(76)
        ENDIF
        IF(MINT(47).EQ.1.AND.(ISTSB.EQ.1.OR.ISTSB.EQ.2.OR.ISTSB.EQ.6))
     &  THEN
          TAU=1.
        ELSEIF(MVAR.EQ.1) THEN
          TAU=TAUMIN*(TAUMAX/TAUMIN)**VVAR
        ELSEIF(MVAR.EQ.2) THEN
          TAU=TAUMAX*TAUMIN/(TAUMIN+(TAUMAX-TAUMIN)*VVAR)
        ELSEIF(MVAR.EQ.3.OR.MVAR.EQ.5) THEN
          RATGEN=(TAURE+TAUMAX)/(TAURE+TAUMIN)*TAUMIN/TAUMAX
          TAU=TAURE*TAUMIN/((TAURE+TAUMIN)*RATGEN**VVAR-TAUMIN)
        ELSEIF(MVAR.EQ.4.OR.MVAR.EQ.6) THEN
          AUPP=ATAN((TAUMAX-TAURE)/GAMRE)
          ALOW=ATAN((TAUMIN-TAURE)/GAMRE)
          TAU=TAURE+GAMRE*TAN(ALOW+(AUPP-ALOW)*VVAR)
        ELSE
          AUPP=LOG(MAX(2E-6,1.-TAUMAX))
          ALOW=LOG(MAX(2E-6,1.-TAUMIN))
          TAU=1.-EXP(AUPP+VVAR*(ALOW-AUPP))
        ENDIF
        VINT(21)=MIN(TAUMAX,MAX(TAUMIN,TAU))

C...Convert VVAR to y* variable.
      ELSEIF(IVAR.EQ.2) THEN
        YSTMIN=VINT(12)
        YSTMAX=VINT(32)
        TAUE=VINT(21)
        IF(ISTSB.GE.3.AND.ISTSB.LE.5) TAUE=VINT(26)
        IF(MINT(47).EQ.1) THEN
          YST=0.
        ELSEIF(MINT(47).EQ.2) THEN
          YST=-0.5*LOG(TAUE)
        ELSEIF(MINT(47).EQ.3) THEN
          YST=0.5*LOG(TAUE)
        ELSEIF(MVAR.EQ.1) THEN
          YST=YSTMIN+(YSTMAX-YSTMIN)*SQRT(VVAR)
        ELSEIF(MVAR.EQ.2) THEN
          YST=YSTMAX-(YSTMAX-YSTMIN)*SQRT(1.-VVAR)
        ELSEIF(MVAR.EQ.3) THEN
          AUPP=ATAN(EXP(YSTMAX))
          ALOW=ATAN(EXP(YSTMIN))
          YST=LOG(TAN(ALOW+(AUPP-ALOW)*VVAR))
        ELSEIF(MVAR.EQ.4) THEN
          YST0=-0.5*LOG(TAUE)
          AUPP=LOG(MAX(1E-6,EXP(YST0-YSTMIN)-1.))
          ALOW=LOG(MAX(1E-6,EXP(YST0-YSTMAX)-1.))
          YST=YST0-LOG(1.+EXP(ALOW+VVAR*(AUPP-ALOW)))
        ELSE
          YST0=-0.5*LOG(TAUE)
          AUPP=LOG(MAX(1E-6,EXP(YST0+YSTMIN)-1.))
          ALOW=LOG(MAX(1E-6,EXP(YST0+YSTMAX)-1.))
          YST=LOG(1.+EXP(AUPP+VVAR*(ALOW-AUPP)))-YST0
        ENDIF
        VINT(22)=MIN(YSTMAX,MAX(YSTMIN,YST))

C...Convert VVAR to cos(theta-hat) variable.
      ELSEIF(IVAR.EQ.3) THEN
        RM34=2.*VINT(63)*VINT(64)/(VINT(21)*VINT(2))**2
        IF(ISUB.EQ.83) RM34=MAX(1E-20,RM34)
        RSQM=1.+RM34
        IF(2.*VINT(71)**2/(VINT(21)*VINT(2)).LT.0.0001) RM34=MAX(RM34,
     &  2.*VINT(71)**2/(VINT(21)*VINT(2)))
        CTNMIN=VINT(13)
        CTNMAX=VINT(33)
        CTPMIN=VINT(14)
        CTPMAX=VINT(34)
        IF(MVAR.EQ.1) THEN
          ANEG=CTNMAX-CTNMIN
          APOS=CTPMAX-CTPMIN
          IF(ANEG.GT.0..AND.VVAR*(ANEG+APOS).LE.ANEG) THEN
            VCTN=VVAR*(ANEG+APOS)/ANEG
            CTH=CTNMIN+(CTNMAX-CTNMIN)*VCTN
          ELSE
            VCTP=(VVAR*(ANEG+APOS)-ANEG)/APOS
            CTH=CTPMIN+(CTPMAX-CTPMIN)*VCTP
          ENDIF
        ELSEIF(MVAR.EQ.2) THEN
          RMNMIN=MAX(RM34,RSQM-CTNMIN)
          RMNMAX=MAX(RM34,RSQM-CTNMAX)
          RMPMIN=MAX(RM34,RSQM-CTPMIN)
          RMPMAX=MAX(RM34,RSQM-CTPMAX)
          ANEG=LOG(RMNMIN/RMNMAX)
          APOS=LOG(RMPMIN/RMPMAX)
          IF(ANEG.GT.0..AND.VVAR*(ANEG+APOS).LE.ANEG) THEN
            VCTN=VVAR*(ANEG+APOS)/ANEG
            CTH=RSQM-RMNMIN*(RMNMAX/RMNMIN)**VCTN
          ELSE
            VCTP=(VVAR*(ANEG+APOS)-ANEG)/APOS
            CTH=RSQM-RMPMIN*(RMPMAX/RMPMIN)**VCTP
          ENDIF
        ELSEIF(MVAR.EQ.3) THEN
          RMNMIN=MAX(RM34,RSQM+CTNMIN)
          RMNMAX=MAX(RM34,RSQM+CTNMAX)
          RMPMIN=MAX(RM34,RSQM+CTPMIN)
          RMPMAX=MAX(RM34,RSQM+CTPMAX)
          ANEG=LOG(RMNMAX/RMNMIN)
          APOS=LOG(RMPMAX/RMPMIN)
          IF(ANEG.GT.0..AND.VVAR*(ANEG+APOS).LE.ANEG) THEN
            VCTN=VVAR*(ANEG+APOS)/ANEG
            CTH=RMNMIN*(RMNMAX/RMNMIN)**VCTN-RSQM
          ELSE
            VCTP=(VVAR*(ANEG+APOS)-ANEG)/APOS
            CTH=RMPMIN*(RMPMAX/RMPMIN)**VCTP-RSQM
          ENDIF
        ELSEIF(MVAR.EQ.4) THEN
          RMNMIN=MAX(RM34,RSQM-CTNMIN)
          RMNMAX=MAX(RM34,RSQM-CTNMAX)
          RMPMIN=MAX(RM34,RSQM-CTPMIN)
          RMPMAX=MAX(RM34,RSQM-CTPMAX)
          ANEG=1./RMNMAX-1./RMNMIN
          APOS=1./RMPMAX-1./RMPMIN
          IF(ANEG.GT.0..AND.VVAR*(ANEG+APOS).LE.ANEG) THEN
            VCTN=VVAR*(ANEG+APOS)/ANEG
            CTH=RSQM-1./(1./RMNMIN+ANEG*VCTN)
          ELSE
            VCTP=(VVAR*(ANEG+APOS)-ANEG)/APOS
            CTH=RSQM-1./(1./RMPMIN+APOS*VCTP)
          ENDIF
        ELSEIF(MVAR.EQ.5) THEN
          RMNMIN=MAX(RM34,RSQM+CTNMIN)
          RMNMAX=MAX(RM34,RSQM+CTNMAX)
          RMPMIN=MAX(RM34,RSQM+CTPMIN)
          RMPMAX=MAX(RM34,RSQM+CTPMAX)
          ANEG=1./RMNMIN-1./RMNMAX
          APOS=1./RMPMIN-1./RMPMAX
          IF(ANEG.GT.0..AND.VVAR*(ANEG+APOS).LE.ANEG) THEN
            VCTN=VVAR*(ANEG+APOS)/ANEG
            CTH=1./(1./RMNMIN-ANEG*VCTN)-RSQM
          ELSE
            VCTP=(VVAR*(ANEG+APOS)-ANEG)/APOS
            CTH=1./(1./RMPMIN-APOS*VCTP)-RSQM
          ENDIF
        ENDIF
        IF(CTH.LT.0.) CTH=MIN(CTNMAX,MAX(CTNMIN,CTH))
        IF(CTH.GT.0.) CTH=MIN(CTPMAX,MAX(CTPMIN,CTH))
        VINT(23)=CTH

C...Convert VVAR to tau' variable.
      ELSEIF(IVAR.EQ.4) THEN
        TAU=VINT(21)
        TAUPMN=VINT(16)
        TAUPMX=VINT(36)
        IF(MINT(47).EQ.1) THEN
          TAUP=1.
        ELSEIF(MVAR.EQ.1) THEN
          TAUP=TAUPMN*(TAUPMX/TAUPMN)**VVAR
        ELSEIF(MVAR.EQ.2) THEN
          AUPP=(1.-TAU/TAUPMX)**4
          ALOW=(1.-TAU/TAUPMN)**4
          TAUP=TAU/MAX(1E-7,1.-(ALOW+(AUPP-ALOW)*VVAR)**0.25)
        ELSE
          AUPP=LOG(MAX(2E-6,1.-TAUPMX))
          ALOW=LOG(MAX(2E-6,1.-TAUPMN))
          TAUP=1.-EXP(AUPP+VVAR*(ALOW-AUPP))
        ENDIF
        VINT(26)=MIN(TAUPMX,MAX(TAUPMN,TAUP))

C...Selection of extra variables needed in 2 -> 3 process:
C...pT1, pT2, phi1, phi2, y3 for three outgoing particles.
C...Since no options are available, the functions of RYKLIM
C...and RYKMAP are joint for these choices.
      ELSEIF(IVAR.EQ.5) THEN

C...Read out total energy and particle masses.
        MINT(51)=0
        SHP=VINT(26)*VINT(2)
        SHPR=SQRT(SHP)
        PM1=VINT(201)
        PM2=VINT(206)
        PM3=SQRT(VINT(21))*VINT(1)
        IF(PM1+PM2+PM3.GT.0.9999*SHPR) THEN
          MINT(51)=1
          RETURN
        ENDIF
        PMRS1=VINT(204)**2
        PMRS2=VINT(209)**2

C...Select transverse momenta according to dpT^2/(pT^2+M^2)^2.
        PTSMX1=((SHP-PM1**2-(PM2+PM3)**2)**2-(2.*PM1*(PM2+PM3))**2)/
     &  (4.*SHP)
        IF(CKIN(52).GT.0.) PTSMX1=MIN(PTSMX1,CKIN(52)**2)
        PTSMN1=CKIN(51)**2
        PTS1=MAX(PTSMN1,(PMRS1+PTSMN1)*(PMRS1+PTSMX1)/
     &  (PMRS1+PTSMN1+PYR(0)*(PTSMX1-PTSMN1))-PMRS1)
        PTSMX2=((SHP-PM2**2-(PM1+PM3)**2)**2-(2.*PM2*(PM1+PM3))**2)/
     &  (4.*SHP)
        IF(CKIN(54).GT.0.) PTSMX2=MIN(PTSMX2,CKIN(54)**2)
        PTSMN2=CKIN(53)**2
        PTS2=MAX(PTSMN2,(PMRS2+PTSMN2)*(PMRS2+PTSMX2)/
     &  (PMRS2+PTSMN2+PYR(0)*(PTSMX2-PTSMN2))-PMRS2)
        PHI1=PARU(2)*PYR(0)
        PHI2=PARU(2)*PYR(0)
        PHIR=PHI2-PHI1
        PTS3=MAX(0.,PTS1+PTS2+2.*SQRT(PTS1*PTS2)*COS(PHIR))
        IF(PTS3.LT.CKIN(55)**2.OR.(CKIN(56).GT.0..AND.PTS3.GT.
     &  CKIN(56)**2)) THEN
          MINT(51)=1
          RETURN
        ENDIF

C...Calculate transverse masses and check phase space not closed.
        PMS1=PM1**2+PTS1
        PMS2=PM2**2+PTS2
        PMS3=PM3**2+PTS3
        PMT1=SQRT(PMS1)
        PMT2=SQRT(PMS2)
        PMT3=SQRT(PMS3)
        PM12=(PMT1+PMT2)**2
        IF(PMT1+PMT2+PMT3.GT.0.9999*SHPR) THEN
          MINT(51)=1
          RETURN
        ENDIF

C...Select rapidity for particle 3 and check phase space not closed.
        Y3MAX=LOG((SHP+PMS3-PM12+SQRT(MAX(0.,(SHP-PMS3-PM12)**2-
     &  4.*PMS3*PM12)))/(2.*SHPR*PMT3))
        IF(Y3MAX.LT.1E-6) THEN
          MINT(51)=1
          RETURN
        ENDIF
        Y3=(2.*PYR(0)-1.)*0.999999*Y3MAX
        PZ3=PMT3*SINH(Y3)
        PE3=PMT3*COSH(Y3)

C...Find momentum transfers in two mirror solutions (in 1-2 frame).
        PZ12=-PZ3
        PE12=SHPR-PE3
        PMS12=PE12**2-PZ12**2
        SQL12=SQRT(MAX(0.,(PMS12-PMS1-PMS2)**2-4.*PMS1*PMS2))
        IF(SQL12.LT.1E-6*SHP) THEN
          MINT(51)=1
          RETURN
        ENDIF
        PMM1=PMS12+PMS1-PMS2
        PMM2=PMS12+PMS2-PMS1
        TFAC=-SHPR/(2.*PMS12)
        T1P=TFAC*(PE12-PZ12)*(PMM1-SQL12)
        T1N=TFAC*(PE12-PZ12)*(PMM1+SQL12)
        T2P=TFAC*(PE12+PZ12)*(PMM2-SQL12)
        T2N=TFAC*(PE12+PZ12)*(PMM2+SQL12)

C...Construct relative mirror weights and make choice.
        WTPU=1./((T1P-PMRS1)*(T2P-PMRS2))**2
        WTNU=1./((T1N-PMRS1)*(T2N-PMRS2))**2
        WTP=WTPU/(WTPU+WTNU)
        WTN=WTNU/(WTPU+WTNU)
        EPS=1.
        IF(WTN.GT.PYR(0)) EPS=-1.

C...Store result of variable choice and associated weights.
        VINT(202)=PTS1
        VINT(207)=PTS2
        VINT(203)=PHI1
        VINT(208)=PHI2
        VINT(205)=(PMRS1+PTS1)**2*(PTSMX1-PTSMN1)/((PMRS1+PTSMN1)*
     &  (PMRS1+PTSMX1))
        VINT(210)=(PMRS2+PTS2)**2*(PTSMX2-PTSMN2)/((PMRS2+PTSMN2)*
     &  (PMRS2+PTSMX2))
        VINT(211)=Y3
        VINT(212)=Y3MAX
        VINT(213)=EPS
        IF(EPS.GT.0.) THEN
          VINT(214)=1./WTP
          VINT(215)=T1P
          VINT(216)=T2P
        ELSE
          VINT(214)=1./WTN
          VINT(215)=T1N
          VINT(216)=T2N
        ENDIF
        VINT(217)=-0.5*TFAC*(PE12-PZ12)*(PMM2+EPS*SQL12)
        VINT(218)=-0.5*TFAC*(PE12+PZ12)*(PMM1+EPS*SQL12)
        VINT(219)=0.5*(PMS12-PTS3)
        VINT(220)=SQL12
      ENDIF

      RETURN
      END

C***********************************************************************

      SUBROUTINE RYSIGH(NCHN,SIGS)

C...Differential matrix elements for all included subprocesses.
C...Note that what is coded is (disregarding the COMFAC factor)
C...1) for 2 -> 1 processes: s-hat/pi*d(sigma-hat), where,
C...when d(sigma-hat) is given in the zero-width limit, the delta
C...function in tau is replaced by a (modified) Breit-Wigner:
C...1/pi*s*H_res/((s*tau-m_res^2)^2+H_res^2),
C...where H_res = s-hat/m_res*Gamma_res(s-hat);
C...2) for 2 -> 2 processes: (s-hat)**2/pi*d(sigma-hat)/d(t-hat);
C...i.e., dimensionless quantities.
C...3) for 2 -> 3 processes: abs(M)^2, where the total cross-section is
C...Integral abs(M)^2/(2shat') * (prod_(i=1)^3 d^3p_i/((2pi)^3*2E_i)) *
C...(2pi)^4 delta^4(P - sum p_i).
C...COMFAC contains the factor pi/s (or equivalent) and
C...the conversion factor from GeV^-2 to mb.
      PARAMETER (KSZJ=4000)

      
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      integer N,NPAD,K
      double precision P,V
      SAVE /PYJETS/


      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/LUDAT3/MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5)
      COMMON/RYSUBS/MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200)
      COMMON/RYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/RYINT1/MINT(400),VINT(400)
      COMMON/RYINT2/ISET(200),KFPR(200,2),COEF(200,20),ICOL(40,4,2)
      COMMON/RYINT3/XSFX(2,-40:40),ISIG(1000,3),SIGH(1000)
      COMMON/RYINT4/WIDP(21:40,0:40),WIDE(21:40,0:40),WIDS(21:40,3)
      COMMON/RYINT5/NGEN(0:200,3),XSEC(0:200,3)
      SAVE /LUDAT1/,/LUDAT2/,/LUDAT3/
      SAVE /RYSUBS/,/RYPARS/,/RYINT1/,/RYINT2/,/RYINT3/,/RYINT4/,
     &/RYINT5/
      DIMENSION X(2),XPQ(-25:25),KFAC(2,-40:40),WDTP(0:40),
     &WDTE(0:40,0:5),HGZ(6,3),HL3(3),HR3(3),HL4(3),HR4(3)
      COMPLEX A004,A204,A114,A00U,A20U,A11U

C...The following gives an interface for process 131, gg -> Zqq,
C...to the matrix element package of Ronald Kleiss.
      COMMON/RKBBVC/RKMQ,RKMZ,RKGZ,RKVQ,RKAQ,RKVL,RKAL
      SAVE /RKBBVC/
      DIMENSION RKG1(0:3),RKG2(0:3),RKQ1(0:3),RKQ2(0:3),RKL1(0:3),
     &RKL2(0:3)

C...Reset number of channels and cross-section.
      NCHN=0
      SIGS=0.

C...Convert H' or A process into equivalent H one.
      ISUB=MINT(1)
      ISUBSV=ISUB
      IHIGG=1
      KFHIGG=25
      IF((ISUB.GE.151.AND.ISUB.LE.160).OR.(ISUB.GE.171.AND.
     &ISUB.LE.180)) THEN
        IHIGG=2
        IF(MOD(ISUB-1,10).GE.5) IHIGG=3
        KFHIGG=33+IHIGG
        IF(ISUB.EQ.151.OR.ISUB.EQ.156) ISUB=3
        IF(ISUB.EQ.152.OR.ISUB.EQ.157) ISUB=102
        IF(ISUB.EQ.153.OR.ISUB.EQ.158) ISUB=103
        IF(ISUB.EQ.171.OR.ISUB.EQ.176) ISUB=24
        IF(ISUB.EQ.172.OR.ISUB.EQ.177) ISUB=26
        IF(ISUB.EQ.173.OR.ISUB.EQ.178) ISUB=123
        IF(ISUB.EQ.174.OR.ISUB.EQ.179) ISUB=124
      ENDIF

C...Read kinematical variables and limits.
      ISTSB=ISET(ISUB)
      TAUMIN=VINT(11)
      YSTMIN=VINT(12)
      CTNMIN=VINT(13)
      CTPMIN=VINT(14)
      TAUPMN=VINT(16)
      TAU=VINT(21)
      YST=VINT(22)
      CTH=VINT(23)
      XT2=VINT(25)
      TAUP=VINT(26)
      TAUMAX=VINT(31)
      YSTMAX=VINT(32)
      CTNMAX=VINT(33)
      CTPMAX=VINT(34)
      TAUPMX=VINT(36)

C...Derive kinematical quantities.
      TAUE=TAU
      IF(ISTSB.GE.3.AND.ISTSB.LE.5) TAUE=TAUP
      X(1)=SQRT(TAUE)*EXP(YST)
      X(2)=SQRT(TAUE)*EXP(-YST)
      IF(MINT(45).EQ.2.AND.ISTSB.GE.1) THEN
        IF(X(1).GT.0.9999) RETURN
      ELSEIF(MINT(45).EQ.3) THEN
        X(1)=MIN(0.9999989,X(1))
      ENDIF
      IF(MINT(46).EQ.2.AND.ISTSB.GE.1) THEN
        IF(X(2).GT.0.9999) RETURN
      ELSEIF(MINT(46).EQ.3) THEN
        X(2)=MIN(0.9999989,X(2))
      ENDIF
      SH=TAU*VINT(2)
      SQM3=VINT(63)
      SQM4=VINT(64)
      RM3=SQM3/SH
      RM4=SQM4/SH
      BE34=SQRT(MAX(0.,(1.-RM3-RM4)**2-4.*RM3*RM4))
      RPTS=4.*VINT(71)**2/SH
      BE34L=SQRT(MAX(0.,(1.-RM3-RM4)**2-4.*RM3*RM4-RPTS))
      RM34=2.*RM3*RM4
      IF(ISUB.EQ.83) RM34=MAX(1E-20,RM34)
      RSQM=1.+RM34
      RTHM=(4.*RM3*RM4+RPTS)/(1.-RM3-RM4+BE34L)
      TH=-0.5*SH*MAX(RTHM,1.-RM3-RM4-BE34*CTH)
      UH=-0.5*SH*MAX(RTHM,1.-RM3-RM4+BE34*CTH)
      SQPTH=MAX(VINT(71)**2,0.25*SH*BE34**2*(1.-CTH**2))
      SH2=SH**2
      TH2=TH**2
      UH2=UH**2

C...Choice of Q2 scale.
      IF(ISTSB.EQ.1.OR.ISTSB.EQ.3.OR.ISTSB.EQ.5) THEN
        Q2=SH
      ELSEIF(MOD(ISTSB,2).EQ.0.OR.ISTSB.EQ.9) THEN
        IF(MSTP(32).EQ.1) THEN
          Q2=2.*SH*TH*UH/(SH**2+TH**2+UH**2)
        ELSEIF(MSTP(32).EQ.2) THEN
          Q2=SQPTH+0.5*(SQM3+SQM4)
        ELSEIF(MSTP(32).EQ.3) THEN
          Q2=MIN(-TH,-UH)
        ELSEIF(MSTP(32).EQ.4) THEN
          Q2=SH
        ENDIF
        IF(ISTSB.EQ.9.AND.MSTP(82).GE.2) Q2=Q2+PARP(82)**2
      ENDIF

C...Store derived kinematical quantities.
      VINT(41)=X(1)
      VINT(42)=X(2)
      VINT(44)=SH
      VINT(43)=SQRT(SH)
      VINT(45)=TH
      VINT(46)=UH
      VINT(48)=SQPTH
      VINT(47)=SQRT(SQPTH)
      VINT(50)=TAUP*VINT(2)
      VINT(49)=SQRT(MAX(0.,VINT(50)))
      VINT(52)=Q2
      VINT(51)=SQRT(Q2)

C...Calculate parton structure functions.
      IF(ISTSB.LE.0) GOTO 150
      IF(MINT(47).GE.2) THEN
        Q2SF=Q2
        IF(ISTSB.GE.3.AND.ISTSB.LE.5) THEN
          Q2SF=PMAS(23,1)**2
          IF(ISUB.EQ.8.OR.ISUB.EQ.76.OR.ISUB.EQ.77.OR.ISUB.EQ.124.OR.
     &    ISUB.EQ.174.OR.ISUB.EQ.179) Q2SF=PMAS(24,1)**2
        ENDIF
        DO 100 I=3-MIN(2,MINT(45)),MIN(2,MINT(46))
        XSF=X(I)
        IF(ISTSB.EQ.9) XSF=X(I)/VINT(142+I)
        CALL RYSTFU(MINT(10+I),XSF,Q2SF,XPQ)
        DO 100 KFL=-25,25
  100   XSFX(I,KFL)=XPQ(KFL)
      ENDIF

C...Calculate alpha_em, alpha_strong and K-factor.
      AEM=ULALEM(Q2)
      IF(MSTP(33).NE.3) AS=ULALPS(Q2)
      FACK=1.
      FACA=1.
      IF(MSTP(33).EQ.1) THEN
        FACK=PARP(31)
      ELSEIF(MSTP(33).EQ.2) THEN
        FACK=PARP(31)
        FACA=PARP(32)/PARP(31)
      ELSEIF(MSTP(33).EQ.3) THEN
        Q2AS=PARP(33)*Q2
        IF(ISTSB.EQ.9.AND.MSTP(82).GE.2) Q2AS=Q2AS+
     &  PARU(112)*PARP(82)
        AS=ULALPS(Q2AS)
      ENDIF
      VINT(138)=1.

C...Set flags for allowed reacting partons/leptons.
      DO 130 I=1,2
      DO 110 J=-25,25
  110 KFAC(I,J)=0
      IF(MINT(44+I).EQ.1) THEN
        KFAC(I,MINT(10+I))=1
      ELSEIF(MINT(40+I).EQ.1.AND.MSTP(12).EQ.0) THEN
        KFAC(I,MINT(10+I))=1
        KFAC(I,22)=1
        KFAC(I,24)=1
        KFAC(I,-24)=1
      ELSE
        DO 120 J=-25,25
        KFAC(I,J)=KFIN(I,J)
        IF(ABS(J).GT.MSTP(54).AND.J.LE.10) KFAC(I,J)=0
        IF(ABS(J).NE.21) THEN
          IF(XSFX(I,J).LT.1.E-10) KFAC(I,J)=0
        ELSE
          IF(XSFX(I,0).LT.1.E-10) KFAC(I,21)=0
        ENDIF
  120   CONTINUE
      ENDIF
  130 CONTINUE

C...Lower and upper limit for fermion flavour loops.
      MIN1=0
      MAX1=0
      MIN2=0
      MAX2=0
      DO 140 J=-20,20
      IF(KFAC(1,-J).EQ.1) MIN1=-J
      IF(KFAC(1,J).EQ.1) MAX1=J
      IF(KFAC(2,-J).EQ.1) MIN2=-J
      IF(KFAC(2,J).EQ.1) MAX2=J
  140 CONTINUE
      MINA=MIN(MIN1,MIN2)
      MAXA=MAX(MAX1,MAX2)

C...Common conversion factors (including Jacobian) for subprocesses.
      SQMZ=PMAS(23,1)**2
      SQMW=PMAS(24,1)**2
      SQMH=PMAS(KFHIGG,1)**2
      GMMH=PMAS(KFHIGG,1)*PMAS(KFHIGG,2)
      SQMZP=PMAS(32,1)**2
      SQMWP=PMAS(34,1)**2
      SQMHC=PMAS(37,1)**2
      SQMLQ=PMAS(39,1)**2
      SQMR=PMAS(40,1)**2
      XW=PARU(102)
      XWC=1./(16.*XW*(1.-XW))

C...Phase space integral in tau.
      COMFAC=PARU(1)*PARU(5)/VINT(2)
      IF(MINT(43).EQ.4) COMFAC=COMFAC*FACK
      IF((MINT(47).GE.2.OR.(ISTSB.GE.3.AND.ISTSB.LE.5)).AND.
     &ISTSB.NE.9) THEN
        ATAU1=LOG(TAUMAX/TAUMIN)
        ATAU2=(TAUMAX-TAUMIN)/(TAUMAX*TAUMIN)
        H1=COEF(ISUB,1)+(ATAU1/ATAU2)*COEF(ISUB,2)/TAU
        IF(MINT(72).GE.1) THEN
          TAUR1=VINT(73)
          GAMR1=VINT(74)
          ATAUD=LOG(TAUMAX/TAUMIN*(TAUMIN+TAUR1)/(TAUMAX+TAUR1))
          ATAU3=ATAUD/TAUR1
          IF(ATAUD.GT.1E-6) H1=H1+(ATAU1/ATAU3)*COEF(ISUB,3)/(TAU+TAUR1)
          ATAUD=ATAN((TAUMAX-TAUR1)/GAMR1)-ATAN((TAUMIN-TAUR1)/GAMR1)
          ATAU4=ATAUD/GAMR1
          IF(ATAUD.GT.1E-6) H1=H1+
     &    (ATAU1/ATAU4)*COEF(ISUB,4)*TAU/((TAU-TAUR1)**2+GAMR1**2)
        ENDIF
        IF(MINT(72).EQ.2) THEN
          TAUR2=VINT(75)
          GAMR2=VINT(76)
          ATAUD=LOG(TAUMAX/TAUMIN*(TAUMIN+TAUR2)/(TAUMAX+TAUR2))
          ATAU5=ATAUD/TAUR2
          IF(ATAUD.GT.1E-6) H1=H1+(ATAU1/ATAU5)*COEF(ISUB,5)/(TAU+TAUR2)
          ATAUD=ATAN((TAUMAX-TAUR2)/GAMR2)-ATAN((TAUMIN-TAUR2)/GAMR2)
          ATAU6=ATAUD/GAMR2
          IF(ATAUD.GT.1E-6) H1=H1+
     &    (ATAU1/ATAU6)*COEF(ISUB,6)*TAU/((TAU-TAUR2)**2+GAMR2**2)
        ENDIF
        IF(MINT(47).EQ.5.AND.(ISTSB.LE.2.OR.ISTSB.GE.6)) THEN
          ATAU7=LOG(MAX(2E-6,1.-TAUMIN)/MAX(2E-6,1.-TAUMAX))
          H1=H1+(ATAU1/ATAU7)*COEF(ISUB,7)*TAU/MAX(2E-6,1.-TAU)
        ENDIF
        COMFAC=COMFAC*ATAU1/(TAU*H1)
      ENDIF

C...Phase space integral in y*.
      IF(MINT(47).GE.4.AND.ISTSB.NE.9) THEN
        AYST0=YSTMAX-YSTMIN
        AYST1=0.5*(YSTMAX-YSTMIN)**2
        AYST2=AYST1
        AYST3=2.*(ATAN(EXP(YSTMAX))-ATAN(EXP(YSTMIN)))
        H2=(AYST0/AYST1)*COEF(ISUB,8)*(YST-YSTMIN)+
     &  (AYST0/AYST2)*COEF(ISUB,9)*(YSTMAX-YST)+
     &  (AYST0/AYST3)*COEF(ISUB,10)/COSH(YST)
        IF(MINT(45).EQ.3) THEN
          YST0=-0.5*LOG(TAUE)
          AYST4=LOG(MAX(1E-6,EXP(YST0-YSTMIN)-1.)/
     &    MAX(1E-6,EXP(YST0-YSTMAX)-1.))
          H2=H2+(AYST0/AYST4)*COEF(ISUB,11)/MAX(1E-6,1.-EXP(YST-YST0))
        ENDIF
        IF(MINT(46).EQ.3) THEN
          YST0=-0.5*LOG(TAUE)
          AYST5=LOG(MAX(1E-6,EXP(YST0+YSTMAX)-1.)/
     &    MAX(1E-6,EXP(YST0+YSTMIN)-1.))
          H2=H2+(AYST0/AYST5)*COEF(ISUB,12)/MAX(1E-6,1.-EXP(-YST-YST0))
        ENDIF
        COMFAC=COMFAC*AYST0/H2
      ENDIF

C...2 -> 1 processes: reduction in angular part of phase space integral
C...for case of decaying resonance.
      ACTH0=CTNMAX-CTNMIN+CTPMAX-CTPMIN
      IF((ISTSB.EQ.1.OR.ISTSB.EQ.3.OR.ISTSB.EQ.5).AND.
     &MDCY(KFPR(ISUBSV,1),1).EQ.1) THEN
        IF(KFPR(ISUB,1).EQ.25.OR.KFPR(ISUB,1).EQ.37.OR.KFPR(ISUB,1).EQ.
     &  39) THEN
          COMFAC=COMFAC*0.5*ACTH0
        ELSE
          COMFAC=COMFAC*0.125*(3.*ACTH0+CTNMAX**3-CTNMIN**3+
     &    CTPMAX**3-CTPMIN**3)
        ENDIF

C...2 -> 2 processes: angular part of phase space integral.
      ELSEIF(ISTSB.EQ.2.OR.ISTSB.EQ.4.OR.ISTSB.EQ.6) THEN
        ACTH1=LOG((MAX(RM34,RSQM-CTNMIN)*MAX(RM34,RSQM-CTPMIN))/
     &  (MAX(RM34,RSQM-CTNMAX)*MAX(RM34,RSQM-CTPMAX)))
        ACTH2=LOG((MAX(RM34,RSQM+CTNMAX)*MAX(RM34,RSQM+CTPMAX))/
     &  (MAX(RM34,RSQM+CTNMIN)*MAX(RM34,RSQM+CTPMIN)))
        ACTH3=1./MAX(RM34,RSQM-CTNMAX)-1./MAX(RM34,RSQM-CTNMIN)+
     &  1./MAX(RM34,RSQM-CTPMAX)-1./MAX(RM34,RSQM-CTPMIN)
        ACTH4=1./MAX(RM34,RSQM+CTNMIN)-1./MAX(RM34,RSQM+CTNMAX)+
     &  1./MAX(RM34,RSQM+CTPMIN)-1./MAX(RM34,RSQM+CTPMAX)
        H3=COEF(ISUB,13)+
     &  (ACTH0/ACTH1)*COEF(ISUB,14)/MAX(RM34,RSQM-CTH)+
     &  (ACTH0/ACTH2)*COEF(ISUB,15)/MAX(RM34,RSQM+CTH)+
     &  (ACTH0/ACTH3)*COEF(ISUB,16)/MAX(RM34,RSQM-CTH)**2+
     &  (ACTH0/ACTH4)*COEF(ISUB,17)/MAX(RM34,RSQM+CTH)**2
        COMFAC=COMFAC*ACTH0*0.5*BE34/H3

C...2 -> 2 processes: take into account final state Breit-Wigners.
        COMFAC=COMFAC*VINT(80)
      ENDIF

C...2 -> 3, 4 processes: phace space integral in tau'.
      IF(MINT(47).GE.2.AND.ISTSB.GE.3.AND.ISTSB.LE.5) THEN
        ATAUP1=LOG(TAUPMX/TAUPMN)
        ATAUP2=((1.-TAU/TAUPMX)**4-(1.-TAU/TAUPMN)**4)/(4.*TAU)
        H4=COEF(ISUB,18)+
     &  (ATAUP1/ATAUP2)*COEF(ISUB,19)*(1.-TAU/TAUP)**3/TAUP
        IF(MINT(47).EQ.5) THEN
          ATAUP3=LOG(MAX(2E-6,1.-TAUPMN)/MAX(2E-6,1.-TAUPMX))
          H4=H4+(ATAUP1/ATAUP3)*COEF(ISUB,20)*TAUP/MAX(2E-6,1.-TAUP)
        ENDIF
        COMFAC=COMFAC*ATAUP1/H4
      ENDIF

C...2 -> 3, 4 processes: effective W/Z structure functions.
      IF(ISTSB.EQ.3.OR.ISTSB.EQ.4) THEN
        IF(1.-TAU/TAUP.GT.1.E-4) THEN
          FZW=(1.+TAU/TAUP)*LOG(TAUP/TAU)-2.*(1.-TAU/TAUP)
        ELSE
          FZW=1./6.*(1.-TAU/TAUP)**3*TAU/TAUP
        ENDIF
        COMFAC=COMFAC*FZW
      ENDIF

C...2 -> 3 processes: phase space integrals for pT1, pT2, y3, mirror.
      IF(ISTSB.EQ.5) THEN
        COMFAC=COMFAC*VINT(205)*VINT(210)*VINT(212)*VINT(214)/
     &  (128.*PARU(1)**4*VINT(220))*(TAU**2/TAUP)
      ENDIF

C...Phase space integral for low-pT and multiple interactions.
      IF(ISTSB.EQ.9) THEN
        COMFAC=PARU(1)*PARU(5)*FACK*0.5*VINT(2)/SH2
        ATAU1=LOG(2.*(1.+SQRT(1.-XT2))/XT2-1.)
        ATAU2=2.*ATAN(1./XT2-1.)/SQRT(XT2)
        H1=COEF(ISUB,1)+(ATAU1/ATAU2)*COEF(ISUB,2)/SQRT(TAU)
        COMFAC=COMFAC*ATAU1/H1
        AYST0=YSTMAX-YSTMIN
        AYST1=0.5*(YSTMAX-YSTMIN)**2
        AYST3=2.*(ATAN(EXP(YSTMAX))-ATAN(EXP(YSTMIN)))
        H2=(AYST0/AYST1)*COEF(ISUB,8)*(YST-YSTMIN)+
     &  (AYST0/AYST1)*COEF(ISUB,9)*(YSTMAX-YST)+
     &  (AYST0/AYST3)*COEF(ISUB,10)/COSH(YST)
        COMFAC=COMFAC*AYST0/H2
        IF(MSTP(82).LE.1) COMFAC=COMFAC*XT2**2*(1./VINT(149)-1.)
C...For MSTP(82)>=2 an additional factor (xT2/(xT2+VINT(149))**2 is
C...introduced to make cross-section finite for xT2 -> 0.
        IF(MSTP(82).GE.2) COMFAC=COMFAC*XT2**2/(VINT(149)*
     &  (1.+VINT(149)))
      ENDIF

C...Strongly interacting Z_L/W_L model of Dobado, Herrero, Terron.
      IF((MSTP(46).GE.3.AND.MSTP(46).LE.6).AND.(ISUB.EQ.71.OR.ISUB.EQ.
     &72.OR.ISUB.EQ.73.OR.ISUB.EQ.76.OR.ISUB.EQ.77)) THEN
C...Calculate M_R and N_R functions for Higgs-like and QCD-like models.
        IF(MSTP(46).LE.4) THEN
          HDTLH=LOG(PMAS(25,1)/PARP(44))
          HDTMR=(4.5*PARU(1)/SQRT(3.)-74./9.)/8.+HDTLH/12.
          HDTNR=-1./18.+HDTLH/6.
        ELSE
          HDTNM=0.125*(1./(288.*PARU(1)**2)+(246./PARP(45))**2)
          HDTLQ=LOG(PARP(45)/PARP(44))
          HDTMR=-(4.*PARU(1))**2*0.5*HDTNM+HDTLQ/12.
          HDTNR=(4.*PARU(1))**2*HDTNM+HDTLQ/6.
        ENDIF

C...Calculate lowest and next-to-lowest order partial wave amplitudes.
        HDTV=1./(16.*PARU(1)*246.**2)
        A00L=HDTV*SH
        A20L=-0.5*A00L
        A11L=A00L/6.
        HDTLS=LOG(SH/PARP(44)**2)
        A004=(HDTV*SH)**2/(4.*PARU(1))*CMPLX((176.*HDTMR+112.*HDTNR)/3.+
     &  11./27.-(50./9.)*HDTLS,4.*PARU(1))
        A204=(HDTV*SH)**2/(4.*PARU(1))*CMPLX(32.*(HDTMR+2.*HDTNR)/3.+
     &  25./54.-(20./9.)*HDTLS,PARU(1))
        A114=(HDTV*SH)**2/(6.*PARU(1))*CMPLX(4.*(-2.*HDTMR+HDTNR)-
     &  1./18.,PARU(1)/6.)

C...Unitarize partial wave amplitudes with Pade or K-matrix method.
        IF(MSTP(46).EQ.3.OR.MSTP(46).EQ.5) THEN
          A00U=A00L/(1.-A004/A00L)
          A20U=A20L/(1.-A204/A20L)
          A11U=A11L/(1.-A114/A11L)
        ELSE
          A00U=(A00L+REAL(A004))/(1.-CMPLX(0.,A00L+REAL(A004)))
          A20U=(A20L+REAL(A204))/(1.-CMPLX(0.,A20L+REAL(A204)))
          A11U=(A11L+REAL(A114))/(1.-CMPLX(0.,A11L+REAL(A114)))
        ENDIF
      ENDIF

C...A: 2 -> 1, tree diagrams.

  150 IF(ISUB.LE.10) THEN
      IF(ISUB.EQ.1) THEN
C...f + f~ -> gamma*/Z0.
        MINT(61)=2
        CALL RYWIDT(23,SH,WDTP,WDTE)
        HP0=AEM/3.*SH
        HP1=AEM/3.*XWC*SH
        HS=HP1*WDTP(0)
        FACZ=4.*COMFAC*3.
        DO 160 I=MINA,MAXA
        IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 160
        EI=KCHG(IABS(I),1)/3.
        AI=SIGN(1.,EI)
        VI=AI-4.*EI*XW
        HI0=HP0
        IF(IABS(I).LE.10) HI0=HI0*FACA/3.
        HI1=HP1
        IF(IABS(I).LE.10) HI1=HI1*FACA/3.
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=-I
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACZ*(EI**2/SH2*HI0*HP0*VINT(111)+EI*VI*(1.-SQMZ/SH)/
     &  ((SH-SQMZ)**2+HS**2)*(HI0*HP1+HI1*HP0)*VINT(112)+
     &  (VI**2+AI**2)/((SH-SQMZ)**2+HS**2)*HI1*HP1*VINT(114))
  160   CONTINUE

      ELSEIF(ISUB.EQ.2) THEN
C...f + f~' -> W+/-.
        CALL RYWIDT(24,SH,WDTP,WDTE)
        HP=AEM/(24.*XW)*SH
        HS=HP*WDTP(0)
        FACBW=4.*COMFAC/((SH-SQMW)**2+HS**2)*3.
        DO 180 I=MIN1,MAX1
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 180
        IA=IABS(I)
        DO 170 J=MIN2,MAX2
        IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 170
        JA=IABS(J)
        IF(I*J.GT.0.OR.MOD(IA+JA,2).EQ.0) GOTO 170
        IF((IA.LE.10.AND.JA.GT.10).OR.(IA.GT.10.AND.JA.LE.10)) GOTO 170
        KCHW=(KCHG(IA,1)*ISIGN(1,I)+KCHG(JA,1)*ISIGN(1,J))/3
        HI=HP*2.
        IF(IA.LE.10) HI=HI*VCKM((IA+1)/2,(JA+1)/2)*FACA/3.
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=J
        ISIG(NCHN,3)=1
        HF=HP*(WDTE(0,1)+WDTE(0,(5-KCHW)/2)+WDTE(0,4))
        SIGH(NCHN)=HI*FACBW*HF
  170   CONTINUE
  180   CONTINUE

      ELSEIF(ISUB.EQ.3) THEN
C...f + f~ -> H0 (or H'0, or A0).
        CALL RYWIDT(KFHIGG,SH,WDTP,WDTE)
        HP=AEM/(8.*XW)*SH/SQMW*SH
        HS=HP*WDTP(0)
        HF=HP*(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))
        FACBW=4.*COMFAC/((SH-SQMH)**2+HS**2)
        IF(ABS(SH-SQMH).GT.100.*HS) FACBW=0.
        DO 190 I=MINA,MAXA
        IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 190
        IA=IABS(I)
        RMQ=PMAS(IA,1)**2/SH
        HI=HP*RMQ
        IF(IA.LE.10) HI=HP*RMQ*FACA/3.
        IF(IA.LE.10.AND.MSTP(37).EQ.1) HI=HI*
     &  (LOG(MAX(4.,PARP(37)**2*RMQ*SH/PARU(117)**2))/
     &  LOG(MAX(4.,SH/PARU(117)**2)))**(24./(33.-2.*MSTU(118)))
        IF(MSTP(4).GE.1.OR.IHIGG.GE.2) THEN
          IKFI=1
          IF(IA.LE.10.AND.MOD(IA,2).EQ.0) IKFI=2
          IF(IA.GT.10) IKFI=3
          HI=HI*PARU(150+10*IHIGG+IKFI)**2
        ENDIF
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=-I
        ISIG(NCHN,3)=1
        SIGH(NCHN)=HI*FACBW*HF
  190   CONTINUE

      ELSEIF(ISUB.EQ.4) THEN
C...gamma + W+/- -> W+/-.

      ELSEIF(ISUB.EQ.5) THEN
C...Z0 + Z0 -> H0.
        CALL RYWIDT(25,SH,WDTP,WDTE)
        HP=AEM/(8.*XW)*SH/SQMW*SH
        HS=HP*WDTP(0)
        HF=HP*(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))
        FACBW=4.*COMFAC/((SH-SQMH)**2+HS**2)
        IF(ABS(SH-SQMH).GT.100.*HS) FACBW=0.
        HI=HP/4.
        FACI=8./(PARU(1)**2*(1.-XW))*(AEM*XWC)**2
        DO 210 I=MIN1,MAX1
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 210
        DO 200 J=MIN2,MAX2
        IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 200
        EI=KCHG(IABS(I),1)/3.
        AI=SIGN(1.,EI)
        VI=AI-4.*EI*XW
        EJ=KCHG(IABS(J),1)/3.
        AJ=SIGN(1.,EJ)
        VJ=AJ-4.*EJ*XW
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=J
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACI*(VI**2+AI**2)*(VJ**2+AJ**2)*HI*FACBW*HF
  200   CONTINUE
  210   CONTINUE

      ELSEIF(ISUB.EQ.6) THEN
C...Z0 + W+/- -> W+/-.

      ELSEIF(ISUB.EQ.7) THEN
C...W+ + W- -> Z0.

      ELSEIF(ISUB.EQ.8) THEN
C...W+ + W- -> H0.
        CALL RYWIDT(25,SH,WDTP,WDTE)
        HP=AEM/(8.*XW)*SH/SQMW*SH
        HS=HP*WDTP(0)
        HF=HP*(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))
        FACBW=4.*COMFAC/((SH-SQMH)**2+HS**2)
        IF(ABS(SH-SQMH).GT.100.*HS) FACBW=0.
        HI=HP/2.
        FACI=1./(4.*PARU(1)**2)*(AEM/XW)**2
        DO 230 I=MIN1,MAX1
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 230
        EI=SIGN(1.,FLOAT(I))*KCHG(IABS(I),1)
        DO 220 J=MIN2,MAX2
        IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 220
        EJ=SIGN(1.,FLOAT(J))*KCHG(IABS(J),1)
        IF(EI*EJ.GT.0.) GOTO 220
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=J
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACI*VINT(180+I)*VINT(180+J)*HI*FACBW*HF
  220   CONTINUE
  230   CONTINUE
      ENDIF

C...B: 2 -> 2, tree diagrams.

      ELSEIF(ISUB.LE.20) THEN
      IF(ISUB.EQ.11) THEN
C...f + f' -> f + f'.
C...Gluon exchange (in t or s channel):
        FACQQ1=COMFAC*AS**2*4./9.*(SH2+UH2)/TH2
        FACQQB=COMFAC*AS**2*4./9.*((SH2+UH2)/TH2*FACA-
     &  MSTP(34)*2./3.*UH2/(SH*TH))
        FACQQ2=COMFAC*AS**2*4./9.*((SH2+TH2)/UH2-
     &  MSTP(34)*2./3.*SH2/(TH*UH))
C...gamma, gamma/Z, Z and W exchange (in t channel):
        FACGGF=COMFAC*AEM**2*2.*(SH2+UH2)/TH2
        FACGZF=COMFAC*AEM**2*XWC*4.*SH2/(TH*(TH-SQMZ))
        FACZZF=COMFAC*(AEM*XWC)**2*2.*SH2/(TH-SQMZ)**2
        FACWWF=COMFAC*(0.5*AEM/XW)**2*SH2/(TH-SQMW)**2
        DO 250 I=MIN1,MAX1
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 250
        IA=IABS(I)
        DO 240 J=MIN2,MAX2
        IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 240
        JA=IABS(J)

        IF(IA.LT.10.AND.JA.LT.10) THEN
C...Gluon exchange.
          NCHN=NCHN+1
          ISIG(NCHN,1)=I
          ISIG(NCHN,2)=J
          ISIG(NCHN,3)=1
          SIGH(NCHN)=FACQQ1
          IF(I.EQ.-J) SIGH(NCHN)=FACQQB
          IF(I.EQ.J) THEN
            SIGH(NCHN)=0.5*SIGH(NCHN)
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=J
            ISIG(NCHN,3)=2
            SIGH(NCHN)=0.5*FACQQ2
          ENDIF

        ELSE
C...Electroweak couplings.
          EI=KCHG(IA,1)*ISIGN(1,I)/3.
          AI=SIGN(1.,KCHG(IA,1)+0.5)*ISIGN(1,I)
          VI=AI-4.*EI*XW
          EJ=KCHG(JA,1)*ISIGN(1,J)/3.
          AJ=SIGN(1.,KCHG(JA,1)+0.5)*ISIGN(1,J)
          VJ=AJ-4.*EJ*XW
          EPSIJ=ISIGN(1,I*J)
C...gamma/Z exchange, only gamma exchange, or only Z exchange.
          IF(MSTP(21).GE.1.AND.MSTP(21).LE.4) THEN
            IF(MSTP(21).EQ.1.OR.MSTP(21).EQ.4) THEN
              FACNCF=FACGGF*EI**2*EJ**2+FACGZF*EI*EJ*
     &        (VI*VJ*(1.+UH2/SH2)+AI*AJ*EPSIJ*(1.-UH2/SH2))+
     &        FACZZF*((VI**2+AI**2)*(VJ**2+AJ**2)*(1.+UH2/SH2)+
     &        4.*VI*VJ*AI*AJ*EPSIJ*(1.-UH2/SH2))
            ELSEIF(MSTP(21).EQ.2) THEN
              FACNCF=FACGGF*EI**2*EJ**2
            ELSE
              FACNCF=FACZZF*((VI**2+AI**2)*(VJ**2+AJ**2)*(1.+UH2/SH2)+
     &        4.*VI*VJ*AI*AJ*EPSIJ*(1.-UH2/SH2))
            ENDIF
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=J
            ISIG(NCHN,3)=5
            SIGH(NCHN)=FACNCF
          ENDIF
C...W exchange.
          IF((MSTP(21).EQ.1.OR.MSTP(21).EQ.5).AND.AI*AJ.LT.0.) THEN
            FACCCF=FACWWF*VINT(180+I)*VINT(180+J)
            IF(EPSIJ.LT.0.) FACCCF=FACCCF*UH2/SH2
            IF(IA.GT.10.AND.MOD(IA,2).EQ.0) FACCCF=2.*FACCCF
            IF(JA.GT.10.AND.MOD(JA,2).EQ.0) FACCCF=2.*FACCCF
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=J
            ISIG(NCHN,3)=6
            SIGH(NCHN)=FACCCF
          ENDIF
        ENDIF
  240   CONTINUE
  250   CONTINUE

      ELSEIF(ISUB.EQ.12) THEN
C...f + f~ -> f' + f~' (q + q~ -> q' + q~' only).
        CALL RYWIDT(21,SH,WDTP,WDTE)
        FACQQB=COMFAC*AS**2*4./9.*(TH2+UH2)/SH2*(WDTE(0,1)+WDTE(0,2)+
     &  WDTE(0,3)+WDTE(0,4))
        DO 260 I=MINA,MAXA
        IF(I.EQ.0.OR.IABS(I).GT.MSTP(54).OR.
     &  KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 260
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=-I
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACQQB
  260   CONTINUE

      ELSEIF(ISUB.EQ.13) THEN
C...f + f~ -> g + g (q + q~ -> g + g only).
        FACGG1=COMFAC*AS**2*32./27.*(UH/TH-(2.+MSTP(34)*1./4.)*UH2/SH2)
        FACGG2=COMFAC*AS**2*32./27.*(TH/UH-(2.+MSTP(34)*1./4.)*TH2/SH2)
        DO 270 I=MINA,MAXA
        IF(I.EQ.0.OR.IABS(I).GT.MSTP(54).OR.
     &  KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 270
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=-I
        ISIG(NCHN,3)=1
        SIGH(NCHN)=0.5*FACGG1
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=-I
        ISIG(NCHN,3)=2
        SIGH(NCHN)=0.5*FACGG2
  270   CONTINUE

      ELSEIF(ISUB.EQ.14) THEN
C...f + f~ -> g + gamma (q + q~ -> g + gamma only).
        FACGG=COMFAC*AS*AEM*8./9.*(TH2+UH2)/(TH*UH)
        DO 280 I=MINA,MAXA
        IF(I.EQ.0.OR.IABS(I).GT.MSTP(54).OR.
     &  KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 280
        EI=KCHG(IABS(I),1)/3.
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=-I
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACGG*EI**2
  280   CONTINUE

      ELSEIF(ISUB.EQ.15) THEN
C...f + f~ -> g + Z0 (q + q~ -> g + Z0 only).
        FACZG=COMFAC*AS*AEM/(XW*(1.-XW))*1./18.*
     &  (TH2+UH2+2.*SQM4*SH)/(TH*UH)
        FACZG=FACZG*WIDS(23,2)
        DO 290 I=MINA,MAXA
        IF(I.EQ.0.OR.IABS(I).GT.MSTP(54).OR.
     &  KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 290
        EI=KCHG(IABS(I),1)/3.
        AI=SIGN(1.,EI)
        VI=AI-4.*EI*XW
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=-I
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACZG*(VI**2+AI**2)
  290   CONTINUE

      ELSEIF(ISUB.EQ.16) THEN
C...f + f~' -> g + W+/- (q + q~' -> g + W+/- only).
        FACWG=COMFAC*AS*AEM/XW*2./9.*(TH2+UH2+2.*SQM4*SH)/(TH*UH)
        DO 310 I=MIN1,MAX1
        IA=IABS(I)
        IF(I.EQ.0.OR.IA.GT.10.OR.KFAC(1,I).EQ.0) GOTO 310
        DO 300 J=MIN2,MAX2
        JA=IABS(J)
        IF(J.EQ.0.OR.JA.GT.10.OR.KFAC(2,J).EQ.0) GOTO 300
        IF(I*J.GT.0.OR.MOD(IA+JA,2).EQ.0) GOTO 300
        KCHW=(KCHG(IA,1)*ISIGN(1,I)+KCHG(JA,1)*ISIGN(1,J))/3
        FCKM=VCKM((IA+1)/2,(JA+1)/2)
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=J
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACWG*FCKM*WIDS(24,(5-KCHW)/2)
  300   CONTINUE
  310   CONTINUE

      ELSEIF(ISUB.EQ.17) THEN
C...f + f~ -> g + H0 (q + q~ -> g + H0 only).

      ELSEIF(ISUB.EQ.18) THEN
C...f + f~ -> gamma + gamma.
        FACGG=COMFAC*AEM**2*2.*(TH2+UH2)/(TH*UH)
        DO 320 I=MINA,MAXA
        IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 320
        EI=KCHG(IABS(I),1)/3.
        FCOI=1.
        IF(IABS(I).LE.10) FCOI=FACA/3.
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=-I
        ISIG(NCHN,3)=1
        SIGH(NCHN)=0.5*FACGG*FCOI*EI**4
  320   CONTINUE

      ELSEIF(ISUB.EQ.19) THEN
C...f + f~ -> gamma + Z0.
        FACGZ=COMFAC*2.*AEM**2*XWC*(TH2+UH2+2.*SQM4*SH)/(TH*UH)
        FACGZ=FACGZ*WIDS(23,2)
        DO 330 I=MINA,MAXA
        IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 330
        EI=KCHG(IABS(I),1)/3.
        AI=SIGN(1.,EI)
        VI=AI-4.*EI*XW
        FCOI=1.
        IF(IABS(I).LE.10) FCOI=FACA/3.
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=-I
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACGZ*FCOI*EI**2*(VI**2+AI**2)
  330   CONTINUE

      ELSEIF(ISUB.EQ.20) THEN
C...f + f~' -> gamma + W+/-.
        FACGW=COMFAC*0.5*AEM**2/XW*
     &  ((2.*UH-TH)/(3.*(SH-SQM4)))**2*(TH2+UH2+2.*SQM4*SH)/(TH*UH)
        DO 350 I=MIN1,MAX1
        IA=IABS(I)
        IF(I.EQ.0.OR.IA.GT.20.OR.KFAC(1,I).EQ.0) GOTO 350
        DO 340 J=MIN2,MAX2
        JA=IABS(J)
        IF(J.EQ.0.OR.JA.GT.20.OR.KFAC(2,J).EQ.0) GOTO 340
        IF(I*J.GT.0.OR.MOD(IA+JA,2).EQ.0) GOTO 340
        IF((IA.LE.10.AND.JA.GT.10).OR.(IA.GT.10.AND.JA.LE.10)) GOTO 340
        KCHW=(KCHG(IA,1)*ISIGN(1,I)+KCHG(JA,1)*ISIGN(1,J))/3
        FCKM=1.
        IF(IA.LE.10) FCKM=VCKM((IA+1)/2,(JA+1)/2)
        FCOI=1.
        IF(IA.LE.10) FCOI=FACA/3.
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=J
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACGW*FCOI*FCKM*WIDS(24,(5-KCHW)/2)
  340   CONTINUE
  350   CONTINUE
      ENDIF

      ELSEIF(ISUB.LE.30) THEN
      IF(ISUB.EQ.21) THEN
C...f + f~ -> gamma + H0.

      ELSEIF(ISUB.EQ.22) THEN
C...f + f~ -> (gamma*/Z0) + (gamma*/Z0).
C...Kinematics dependence.
        FACZZ=COMFAC*AEM**2*((TH2+UH2+2.*(SQM3+SQM4)*SH)/(TH*UH)-
     &  SQM3*SQM4*(1./TH2+1./UH2))
C...gamma, gamma/Z interference and Z couplings to final fermion pairs.
        DO 360 I=1,6
        DO 360 J=1,3
  360   HGZ(I,J)=0.
        HBW3=0.
        HBW4=0.
        RADC3=1.+ULALPS(SQM3)/PARU(1)
        RADC4=1.+ULALPS(SQM4)/PARU(1)
        DO 370 I=1,MIN(16,MDCY(23,3))
        IDC=I+MDCY(23,2)-1
        IF(MDME(IDC,1).LT.0) GOTO 370
        IMDM=0
        IF(MDME(IDC,1).EQ.1.OR.MDME(IDC,1).EQ.2) IMDM=1
        IF(MDME(IDC,1).EQ.4.OR.MDME(IDC,1).EQ.5) IMDM=MDME(IDC,1)-2
        IF(I.LE.8) THEN
          EF=KCHG(I,1)/3.
          AF=SIGN(1.,EF+0.1)
          VF=AF-4.*EF*XW
        ELSEIF(I.LE.16) THEN
          EF=KCHG(I+2,1)/3.
          AF=SIGN(1.,EF+0.1)
          VF=AF-4.*EF*XW
        ENDIF
        RM1=PMAS(IABS(KFDP(IDC,1)),1)**2/SQM3
        IF(4.*RM1.LT.1.) THEN
          FCOF=1.
          IF(I.LE.8) FCOF=3.*RADC3
          BE34=SQRT(MAX(0.,1.-4.*RM1))
          IF(IMDM.GE.1) THEN
            HGZ(1,IMDM)=HGZ(1,IMDM)+FCOF*EF**2*(1.+2.*RM1)*BE34
            HGZ(2,IMDM)=HGZ(2,IMDM)+FCOF*EF*VF*(1.+2.*RM1)*BE34
            HGZ(3,IMDM)=HGZ(3,IMDM)+FCOF*(VF**2*(1.+2.*RM1)+
     &      AF**2*(1.-4.*RM1))*BE34
          ENDIF
          HBW3=HBW3+FCOF*(VF**2*(1.+2.*RM1)+AF**2*(1.-4.*RM1))*BE34
        ENDIF
        RM1=PMAS(IABS(KFDP(IDC,1)),1)**2/SQM4
        IF(4.*RM1.LT.1.) THEN
          FCOF=1.
          IF(I.LE.8) FCOF=3.*RADC4
          BE34=SQRT(MAX(0.,1.-4.*RM1))
          IF(IMDM.GE.1) THEN
            HGZ(4,IMDM)=HGZ(4,IMDM)+FCOF*EF**2*(1.+2.*RM1)*BE34
            HGZ(5,IMDM)=HGZ(5,IMDM)+FCOF*EF*VF*(1.+2.*RM1)*BE34
            HGZ(6,IMDM)=HGZ(6,IMDM)+FCOF*(VF**2*(1.+2.*RM1)+
     &      AF**2*(1.-4.*RM1))*BE34
          ENDIF
          HBW4=HBW4+FCOF*(VF**2*(1.+2.*RM1)+AF**2*(1.-4.*RM1))*BE34
        ENDIF
  370   CONTINUE
C...Propagators: as simulated in RYOFSH and as desired.
        GMMZ=PMAS(23,1)*PMAS(23,2)
        HBW3=HBW3*XWC*SQMZ/((SQM3-SQMZ)**2+GMMZ**2)
        HBW4=HBW4*XWC*SQMZ/((SQM4-SQMZ)**2+GMMZ**2)
        MINT(15)=1
        MINT(61)=1
        CALL RYWIDT(23,SQM3,WDTP,WDTE)
        DO 380 J=1,3
        HGZ(1,J)=HGZ(1,J)*VINT(111)/SQM3
        HGZ(2,J)=HGZ(2,J)*VINT(112)/SQM3
  380   HGZ(3,J)=HGZ(3,J)*VINT(114)/SQM3
        MINT(61)=1
        CALL RYWIDT(23,SQM4,WDTP,WDTE)
        DO 390 J=1,3
        HGZ(4,J)=HGZ(4,J)*VINT(111)/SQM4
        HGZ(5,J)=HGZ(5,J)*VINT(112)/SQM4
  390   HGZ(6,J)=HGZ(6,J)*VINT(114)/SQM4
C...Loop over flavours; separate left- and right-handed couplings.
        DO 410 I=MINA,MAXA
        IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 410
        EI=KCHG(IABS(I),1)/3.
        AI=SIGN(1.,EI)
        VI=AI-4.*EI*XW
        VALI=VI-AI
        VARI=VI+AI
        FCOI=1.
        IF(IABS(I).LE.10) FCOI=FACA/3.
        DO 400 J=1,3
        HL3(J)=EI**2*HGZ(1,J)+EI*VALI*HGZ(2,J)+VALI**2*HGZ(3,J)
        HR3(J)=EI**2*HGZ(1,J)+EI*VARI*HGZ(2,J)+VARI**2*HGZ(3,J)
        HL4(J)=EI**2*HGZ(4,J)+EI*VALI*HGZ(5,J)+VALI**2*HGZ(6,J)
  400   HR4(J)=EI**2*HGZ(4,J)+EI*VARI*HGZ(5,J)+VARI**2*HGZ(6,J)
        FACLR=HL3(1)*HL4(1)+HL3(1)*(HL4(2)+HL4(3))+
     &  HL4(1)*(HL3(2)+HL3(3))+HL3(2)*HL4(3)+HL4(2)*HL3(3)+
     &  HR3(1)*HR4(1)+HR3(1)*(HR4(2)+HR4(3))+
     &  HR4(1)*(HR3(2)+HR3(3))+HR3(2)*HR4(3)+HR4(2)*HR3(3)
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=-I
        ISIG(NCHN,3)=1
        SIGH(NCHN)=0.5*FACZZ*FCOI*FACLR/(HBW3*HBW4)
  410   CONTINUE

      ELSEIF(ISUB.EQ.23) THEN
C...f + f~' -> Z0 + W+/-.
        FACZW=COMFAC*0.5*(AEM/XW)**2
        FACZW=FACZW*WIDS(23,2)
        THUH=MAX(TH*UH-SQM3*SQM4,SH*CKIN(3)**2)
        DO 430 I=MIN1,MAX1
        IA=IABS(I)
        IF(I.EQ.0.OR.IA.GT.20.OR.KFAC(1,I).EQ.0) GOTO 430
        DO 420 J=MIN2,MAX2
        JA=IABS(J)
        IF(J.EQ.0.OR.JA.GT.20.OR.KFAC(2,J).EQ.0) GOTO 420
        IF(I*J.GT.0.OR.MOD(IA+JA,2).EQ.0) GOTO 420
        IF((IA.LE.10.AND.JA.GT.10).OR.(IA.GT.10.AND.JA.LE.10)) GOTO 420
        KCHW=(KCHG(IA,1)*ISIGN(1,I)+KCHG(JA,1)*ISIGN(1,J))/3
        EI=KCHG(IA,1)/3.
        AI=SIGN(1.,EI)
        VI=AI-4.*EI*XW
        EJ=KCHG(JA,1)/3.
        AJ=SIGN(1.,EJ)
        VJ=AJ-4.*EJ*XW
        IF(VI+AI.GT.0) THEN
          VISAV=VI
          AISAV=AI
          VI=VJ
          AI=AJ
          VJ=VISAV
          AJ=AISAV
        ENDIF
        FCKM=1.
        IF(IA.LE.10) FCKM=VCKM((IA+1)/2,(JA+1)/2)
        FCOI=1.
        IF(IA.LE.10) FCOI=FACA/3.
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=J
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACZW*FCOI*FCKM*(1./(SH-SQMW)**2*
     &  ((9.-8.*XW)/4.*THUH+(8.*XW-6.)/4.*SH*(SQM3+SQM4))+
     &  (THUH-SH*(SQM3+SQM4))/(2.*(SH-SQMW))*((VJ+AJ)/TH-(VI+AI)/UH)+
     &  THUH/(16.*(1.-XW))*((VJ+AJ)**2/TH2+(VI+AI)**2/UH2)+
     &  SH*(SQM3+SQM4)/(8.*(1.-XW))*(VI+AI)*(VJ+AJ)/(TH*UH))*
     &  WIDS(24,(5-KCHW)/2)
  420   CONTINUE
  430   CONTINUE

      ELSEIF(ISUB.EQ.24) THEN
C...f + f~ -> Z0 + H0 (or H'0, or A0).
        THUH=MAX(TH*UH-SQM3*SQM4,SH*CKIN(3)**2)
        FACHZ=COMFAC*8.*(AEM*XWC)**2*
     &  (THUH+2.*SH*SQMZ)/((SH-SQMZ)**2+SQMZ*PMAS(23,2)**2)
        FACHZ=FACHZ*WIDS(23,2)*WIDS(KFHIGG,2)
        IF(MSTP(4).GE.1.OR.IHIGG.GE.2) FACHZ=FACHZ*
     &  PARU(154+10*IHIGG)**2
        DO 440 I=MINA,MAXA
        IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 440
        EI=KCHG(IABS(I),1)/3.
        AI=SIGN(1.,EI)
        VI=AI-4.*EI*XW
        FCOI=1.
        IF(IABS(I).LE.10) FCOI=FACA/3.
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=-I
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACHZ*FCOI*(VI**2+AI**2)
  440   CONTINUE

      ELSEIF(ISUB.EQ.25) THEN
C...f + f~ -> W+ + W-.
        THUH=MAX(TH*UH-SQM3*SQM4,SH*CKIN(3)**2)
        FACWW=COMFAC*0.25*(AEM/XW)**2
        FACWW=FACWW*WIDS(24,1)
        DO 450 I=MINA,MAXA
        IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 450
        EI=KCHG(IABS(I),1)/3.
        AI=SIGN(1.,EI)
        VI=AI-4.*EI*XW
        FCOI=1.
        IF(IABS(I).LE.10) FCOI=FACA/3.
        DSIGWW=THUH/SH2*(3.-(SH-3.*(SQM3+SQM4))/(SH-SQMZ)*
     &  (VI+AI)/(2.*AI*(1.-XW))+(SH/(SH-SQMZ))**2*
     &  (1.-2.*(SQM3+SQM4)/SH+12.*SQM3*SQM4/SH2)*(VI**2+AI**2)/
     &  (8.*(1.-XW)**2))-2.*SQMZ/(SH-SQMZ)*(VI+AI)/AI+
     &  SQMZ*SH/(SH-SQMZ)**2*(1.-2.*(SQM3+SQM4)/SH)*(VI**2+AI**2)/
     &  (2.*(1.-XW))
        IF(KCHG(IABS(I),1).LT.0) THEN
          DSIGWW=DSIGWW+2.*(1.+SQMZ/(SH-SQMZ)*(VI+AI)/(2.*AI))*
     &    (THUH/(SH*TH)-(SQM3+SQM4)/TH)+THUH/TH2
        ELSE
          DSIGWW=DSIGWW+2.*(1.+SQMZ/(SH-SQMZ)*(VI+AI)/(2.*AI))*
     &    (THUH/(SH*UH)-(SQM3+SQM4)/UH)+THUH/UH2
        ENDIF
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=-I
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACWW*FCOI*DSIGWW
  450   CONTINUE

      ELSEIF(ISUB.EQ.26) THEN
C...f + f~' -> W+/- + H0 (or H'0, or A0).
        THUH=MAX(TH*UH-SQM3*SQM4,SH*CKIN(3)**2)
        FACHW=COMFAC*0.125*(AEM/XW)**2*(THUH+2.*SH*SQMW)/
     &  ((SH-SQMW)**2+SQMW*PMAS(24,2)**2)
        FACHW=FACHW*WIDS(KFHIGG,2)
        IF(MSTP(4).GE.1.OR.IHIGG.GE.2) FACHW=FACHW*
     &  PARU(155+10*IHIGG)**2
        DO 470 I=MIN1,MAX1
        IA=IABS(I)
        IF(I.EQ.0.OR.IA.GT.20.OR.KFAC(1,I).EQ.0) GOTO 470
        DO 460 J=MIN2,MAX2
        JA=IABS(J)
        IF(J.EQ.0.OR.JA.GT.20.OR.KFAC(1,J).EQ.0) GOTO 460
        IF(I*J.GT.0.OR.MOD(IA+JA,2).EQ.0) GOTO 460
        IF((IA.LE.10.AND.JA.GT.10).OR.(IA.GT.10.AND.JA.LE.10)) GOTO 460
        KCHW=(KCHG(IA,1)*ISIGN(1,I)+KCHG(JA,1)*ISIGN(1,J))/3
        FCKM=1.
        IF(IA.LE.10) FCKM=VCKM((IA+1)/2,(JA+1)/2)
        FCOI=1.
        IF(IA.LE.10) FCOI=FACA/3.
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=J
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACHW*FCOI*FCKM*WIDS(24,(5-KCHW)/2)
  460   CONTINUE
  470   CONTINUE

      ELSEIF(ISUB.EQ.27) THEN
C...f + f~ -> H0 + H0.

      ELSEIF(ISUB.EQ.28) THEN
C...f + g -> f + g (q + g -> q + g only).
        FACQG1=COMFAC*AS**2*4./9.*((2.+MSTP(34)*1./4.)*UH2/TH2-UH/SH)*
     &  FACA
        FACQG2=COMFAC*AS**2*4./9.*((2.+MSTP(34)*1./4.)*SH2/TH2-SH/UH)
        DO 490 I=MINA,MAXA
        IF(I.EQ.0.OR.IABS(I).GT.10) GOTO 490
        DO 480 ISDE=1,2
        IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,21).EQ.0) GOTO 480
        IF(ISDE.EQ.2.AND.KFAC(1,21)*KFAC(2,I).EQ.0) GOTO 480
        NCHN=NCHN+1
        ISIG(NCHN,ISDE)=I
        ISIG(NCHN,3-ISDE)=21
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACQG1
        NCHN=NCHN+1
        ISIG(NCHN,ISDE)=I
        ISIG(NCHN,3-ISDE)=21
        ISIG(NCHN,3)=2
        SIGH(NCHN)=FACQG2
  480   CONTINUE
  490   CONTINUE

      ELSEIF(ISUB.EQ.29) THEN
C...f + g -> f + gamma (q + g -> q + gamma only).
        FGQ=COMFAC*FACA*AS*AEM*1./3.*(SH2+UH2)/(-SH*UH)
        DO 510 I=MINA,MAXA
        IF(I.EQ.0.OR.IABS(I).GT.MSTP(54)) GOTO 510
        EI=KCHG(IABS(I),1)/3.
        FACGQ=FGQ*EI**2
        DO 500 ISDE=1,2
        IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,21).EQ.0) GOTO 500
        IF(ISDE.EQ.2.AND.KFAC(1,21)*KFAC(2,I).EQ.0) GOTO 500
        NCHN=NCHN+1
        ISIG(NCHN,ISDE)=I
        ISIG(NCHN,3-ISDE)=21
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACGQ
  500   CONTINUE
  510   CONTINUE

      ELSEIF(ISUB.EQ.30) THEN
C...f + g -> f + Z0 (q + g -> q + Z0 only).
        FZQ=COMFAC*FACA*AS*AEM*XWC*1./3.*
     &  (SH2+UH2+2.*SQM4*TH)/(-SH*UH)
        FZQ=FZQ*WIDS(23,2)
        DO 530 I=MINA,MAXA
        IF(I.EQ.0.OR.IABS(I).GT.MSTP(54)) GOTO 530
        EI=KCHG(IABS(I),1)/3.
        AI=SIGN(1.,EI)
        VI=AI-4.*EI*XW
        FACZQ=FZQ*(VI**2+AI**2)
        DO 520 ISDE=1,2
        IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,21).EQ.0) GOTO 520
        IF(ISDE.EQ.2.AND.KFAC(1,21)*KFAC(2,I).EQ.0) GOTO 520
        NCHN=NCHN+1
        ISIG(NCHN,ISDE)=I
        ISIG(NCHN,3-ISDE)=21
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACZQ
  520   CONTINUE
  530   CONTINUE
      ENDIF

      ELSEIF(ISUB.LE.40) THEN
      IF(ISUB.EQ.31) THEN
C...f + g -> f' + W+/- (q + g -> q' + W+/- only).
        FACWQ=COMFAC*FACA*AS*AEM/XW*1./12.*
     &  (SH2+UH2+2.*SQM4*TH)/(-SH*UH)
        DO 550 I=MINA,MAXA
        IF(I.EQ.0.OR.IABS(I).GT.MSTP(54)) GOTO 550
        IA=IABS(I)
        KCHW=ISIGN(1,KCHG(IA,1)*ISIGN(1,I))
        DO 540 ISDE=1,2
        IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,21).EQ.0) GOTO 540
        IF(ISDE.EQ.2.AND.KFAC(1,21)*KFAC(2,I).EQ.0) GOTO 540
        NCHN=NCHN+1
        ISIG(NCHN,ISDE)=I
        ISIG(NCHN,3-ISDE)=21
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACWQ*VINT(180+I)*WIDS(24,(5-KCHW)/2)
  540   CONTINUE
  550   CONTINUE

      ELSEIF(ISUB.EQ.32) THEN
C...f + g -> f + H0 (q + g -> q + H0 only).

      ELSEIF(ISUB.EQ.33) THEN
C...f + gamma -> f + g (q + gamma -> q + g only).
        FGQ=COMFAC*AS*AEM*8./3.*(SH2+UH2)/(-SH*UH)
        DO 570 I=MINA,MAXA
        IF(I.EQ.0.OR.IABS(I).GT.MSTP(54)) GOTO 570
        EI=KCHG(IABS(I),1)/3.
        FACGQ=FGQ*EI**2
        DO 560 ISDE=1,2
        IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,22).EQ.0) GOTO 560
        IF(ISDE.EQ.2.AND.KFAC(1,22)*KFAC(2,I).EQ.0) GOTO 560
        NCHN=NCHN+1
        ISIG(NCHN,ISDE)=I
        ISIG(NCHN,3-ISDE)=22
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACGQ
  560   CONTINUE
  570   CONTINUE

      ELSEIF(ISUB.EQ.34) THEN
C...f + gamma -> f + gamma.
        FGQ=COMFAC*AEM**2*2.*(SH2+UH2)/(-SH*UH)
        DO 590 I=MINA,MAXA
        IF(I.EQ.0) GOTO 590
        EI=KCHG(IABS(I),1)/3.
        FACGQ=FGQ*EI**4
        DO 580 ISDE=1,2
        IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,22).EQ.0) GOTO 580
        IF(ISDE.EQ.2.AND.KFAC(1,22)*KFAC(2,I).EQ.0) GOTO 580
        NCHN=NCHN+1
        ISIG(NCHN,ISDE)=I
        ISIG(NCHN,3-ISDE)=22
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACGQ
  580   CONTINUE
  590   CONTINUE

      ELSEIF(ISUB.EQ.35) THEN
C...f + gamma -> f + Z0.
        FZQ=COMFAC*AEM**2*XWC*2.*
     &  (SH2+UH2+2.*SQM4*TH)/(SQPTH*SQM4-SH*UH)
        FZQ=FZQ*WIDS(23,2)
        DO 610 I=MINA,MAXA
        IF(I.EQ.0) GOTO 610
        EI=KCHG(IABS(I),1)/3.
        AI=SIGN(1.,EI)
        VI=AI-4.*EI*XW
        FACZQ=FZQ*EI**2*(VI**2+AI**2)
        DO 600 ISDE=1,2
        IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,22).EQ.0) GOTO 600
        IF(ISDE.EQ.2.AND.KFAC(1,22)*KFAC(2,I).EQ.0) GOTO 600
        NCHN=NCHN+1
        ISIG(NCHN,ISDE)=I
        ISIG(NCHN,3-ISDE)=22
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACZQ
  600   CONTINUE
  610   CONTINUE

      ELSEIF(ISUB.EQ.36) THEN
C...f + gamma -> f' + W+/-.
        FWQ=COMFAC*AEM**2/(2.*XW)*
     &  (SH2+UH2+2.*SQM4*TH)/(SQPTH*SQM4-SH*UH)
        DO 630 I=MINA,MAXA
        IF(I.EQ.0) GOTO 630
        IA=IABS(I)
        EIA=ABS(KCHG(IABS(I),1)/3.)
        FACWQ=FWQ*(EIA-SH/(SH+UH))**2
        KCHW=ISIGN(1,KCHG(IA,1)*ISIGN(1,I))
        DO 620 ISDE=1,2
        IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,22).EQ.0) GOTO 620
        IF(ISDE.EQ.2.AND.KFAC(1,22)*KFAC(2,I).EQ.0) GOTO 620
        NCHN=NCHN+1
        ISIG(NCHN,ISDE)=I
        ISIG(NCHN,3-ISDE)=22
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACWQ*VINT(180+I)*WIDS(24,(5-KCHW)/2)
  620   CONTINUE
  630   CONTINUE

      ELSEIF(ISUB.EQ.37) THEN
C...f + gamma -> f + H0.

      ELSEIF(ISUB.EQ.38) THEN
C...f + Z0 -> f + g (q + Z0 -> q + g only).

      ELSEIF(ISUB.EQ.39) THEN
C...f + Z0 -> f + gamma.

      ELSEIF(ISUB.EQ.40) THEN
C...f + Z0 -> f + Z0.
      ENDIF

      ELSEIF(ISUB.LE.50) THEN
      IF(ISUB.EQ.41) THEN
C...f + Z0 -> f' + W+/-.

      ELSEIF(ISUB.EQ.42) THEN
C...f + Z0 -> f + H0.

      ELSEIF(ISUB.EQ.43) THEN
C...f + W+/- -> f' + g (q + W+/- -> q' + g only).

      ELSEIF(ISUB.EQ.44) THEN
C...f + W+/- -> f' + gamma.

      ELSEIF(ISUB.EQ.45) THEN
C...f + W+/- -> f' + Z0.

      ELSEIF(ISUB.EQ.46) THEN
C...f + W+/- -> f' + W+/-.

      ELSEIF(ISUB.EQ.47) THEN
C...f + W+/- -> f' + H0.

      ELSEIF(ISUB.EQ.48) THEN
C...f + H0 -> f + g (q + H0 -> q + g only).

      ELSEIF(ISUB.EQ.49) THEN
C...f + H0 -> f + gamma.

      ELSEIF(ISUB.EQ.50) THEN
C...f + H0 -> f + Z0.
      ENDIF

      ELSEIF(ISUB.LE.60) THEN
      IF(ISUB.EQ.51) THEN
C...f + H0 -> f' + W+/-.

      ELSEIF(ISUB.EQ.52) THEN
C...f + H0 -> f + H0.

      ELSEIF(ISUB.EQ.53) THEN
C...g + g -> f + f~ (g + g -> q + q~ only).
        CALL RYWIDT(21,SH,WDTP,WDTE)
        FACQQ1=COMFAC*AS**2*1./6.*(UH/TH-(2.+MSTP(34)*1./4.)*UH2/SH2)*
     &  (WDTE(0,1)+WDTE(0,2)+WDTE(0,3)+WDTE(0,4))*FACA
        FACQQ2=COMFAC*AS**2*1./6.*(TH/UH-(2.+MSTP(34)*1./4.)*TH2/SH2)*
     &  (WDTE(0,1)+WDTE(0,2)+WDTE(0,3)+WDTE(0,4))*FACA
        IF(KFAC(1,21)*KFAC(2,21).EQ.0) GOTO 640
        NCHN=NCHN+1
        ISIG(NCHN,1)=21
        ISIG(NCHN,2)=21
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACQQ1
        NCHN=NCHN+1
        ISIG(NCHN,1)=21
        ISIG(NCHN,2)=21
        ISIG(NCHN,3)=2
        SIGH(NCHN)=FACQQ2
  640   CONTINUE

      ELSEIF(ISUB.EQ.54) THEN
C...g + gamma -> f + f~ (g + gamma -> q + q~ only).
        CALL RYWIDT(21,SH,WDTP,WDTE)
        WDTESU=0.
        DO 650 I=1,MIN(8,MDCY(21,3))
        EF=KCHG(I,1)/3.
  650   WDTESU=WDTESU+EF**2*(WDTE(I,1)+WDTE(I,2)+WDTE(I,3)+WDTE(I,4))
        FACQQ=COMFAC*AEM*AS*WDTESU*(TH2+UH2)/(TH*UH)
        IF(KFAC(1,21)*KFAC(2,22).NE.0) THEN
          NCHN=NCHN+1
          ISIG(NCHN,1)=21
          ISIG(NCHN,2)=22
          ISIG(NCHN,3)=1
          SIGH(NCHN)=FACQQ
        ENDIF
        IF(KFAC(1,22)*KFAC(2,21).NE.0) THEN
          NCHN=NCHN+1
          ISIG(NCHN,1)=22
          ISIG(NCHN,2)=21
          ISIG(NCHN,3)=1
          SIGH(NCHN)=FACQQ
        ENDIF

      ELSEIF(ISUB.EQ.55) THEN
C...g + Z -> f + f~ (g + Z -> q + q~ only).

      ELSEIF(ISUB.EQ.56) THEN
C...g + W -> f + f'~ (g + W -> q + q'~ only).

      ELSEIF(ISUB.EQ.57) THEN
C...g + H0 -> f + f~ (g + H0 -> q + q~ only).

      ELSEIF(ISUB.EQ.58) THEN
C...gamma + gamma -> f + f~.
        CALL RYWIDT(22,SH,WDTP,WDTE)
        WDTESU=0.
        DO 660 I=1,MIN(12,MDCY(22,3))
        IF(I.LE.8) EF= KCHG(I,1)/3.
        IF(I.GE.9) EF= KCHG(9+2*(I-8),1)/3.
  660   WDTESU=WDTESU+EF**2*(WDTE(I,1)+WDTE(I,2)+WDTE(I,3)+WDTE(I,4))
        FACFF=COMFAC*AEM**2*WDTESU*2.*(TH2+UH2)/(TH*UH)
        IF(KFAC(1,22)*KFAC(2,22).NE.0) THEN
          NCHN=NCHN+1
          ISIG(NCHN,1)=22
          ISIG(NCHN,2)=22
          ISIG(NCHN,3)=1
          SIGH(NCHN)=FACFF
        ENDIF

      ELSEIF(ISUB.EQ.59) THEN
C...gamma + Z0 -> f + f~.

      ELSEIF(ISUB.EQ.60) THEN
C...gamma + W+/- -> f + f~'.
      ENDIF

      ELSEIF(ISUB.LE.70) THEN
      IF(ISUB.EQ.61) THEN
C...gamma + H0 -> f + f~.

      ELSEIF(ISUB.EQ.62) THEN
C...Z0 + Z0 -> f + f~.

      ELSEIF(ISUB.EQ.63) THEN
C...Z0 + W+/- -> f + f~'.

      ELSEIF(ISUB.EQ.64) THEN
C...Z0 + H0 -> f + f~.

      ELSEIF(ISUB.EQ.65) THEN
C...W+ + W- -> f + f~.

      ELSEIF(ISUB.EQ.66) THEN
C...W+/- + H0 -> f + f~'.

      ELSEIF(ISUB.EQ.67) THEN
C...H0 + H0 -> f + f~.

      ELSEIF(ISUB.EQ.68) THEN
C...g + g -> g + g.
        FACGG1=COMFAC*AS**2*9./4.*(SH2/TH2+2.*SH/TH+3.+2.*TH/SH+
     &  TH2/SH2)*FACA
        FACGG2=COMFAC*AS**2*9./4.*(UH2/SH2+2.*UH/SH+3.+2.*SH/UH+
     &  SH2/UH2)*FACA
        FACGG3=COMFAC*AS**2*9./4.*(TH2/UH2+2.*TH/UH+3+2.*UH/TH+UH2/TH2)
        IF(KFAC(1,21)*KFAC(2,21).EQ.0) GOTO 670
        NCHN=NCHN+1
        ISIG(NCHN,1)=21
        ISIG(NCHN,2)=21
        ISIG(NCHN,3)=1
        SIGH(NCHN)=0.5*FACGG1
        NCHN=NCHN+1
        ISIG(NCHN,1)=21
        ISIG(NCHN,2)=21
        ISIG(NCHN,3)=2
        SIGH(NCHN)=0.5*FACGG2
        NCHN=NCHN+1
        ISIG(NCHN,1)=21
        ISIG(NCHN,2)=21
        ISIG(NCHN,3)=3
        SIGH(NCHN)=0.5*FACGG3
  670   CONTINUE

      ELSEIF(ISUB.EQ.69) THEN
C...gamma + gamma -> W+ + W-.
        SQMWE=MAX(0.5*SQMW,SQRT(SQM3*SQM4))
        FPROP=SH2/((SQMWE-TH)*(SQMWE-UH))
        FACWW=COMFAC*6.*AEM**2*(1.-FPROP*(4./3.+2.*SQMWE/SH)+
     &  FPROP**2*(2./3.+2.*(SQMWE/SH)**2))*WIDS(24,1)
        IF(KFAC(1,22)*KFAC(2,22).EQ.0) GOTO 680
        NCHN=NCHN+1
        ISIG(NCHN,1)=22
        ISIG(NCHN,2)=22
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACWW
  680   CONTINUE

      ELSEIF(ISUB.EQ.70) THEN
C...gamma + W+/- -> Z0 + W+/-.
        SQMWE=MAX(0.5*SQMW,SQRT(SQM3*SQM4))
        FPROP=(TH-SQMWE)**2/(-SH*(SQMWE-UH))
        FACZW=COMFAC*6.*AEM**2*((1.-XW)/XW)*
     &  (1.-FPROP*(4./3.+2.*SQMWE/(TH-SQMWE))+
     &  FPROP**2*(2./3.+2.*(SQMWE/(TH-SQMWE))**2))*WIDS(23,2)
        DO 700 KCHW=1,-1,-2
        DO 690 ISDE=1,2
        IF(KFAC(ISDE,22)*KFAC(3-ISDE,24*KCHW).EQ.0) GOTO 690
        NCHN=NCHN+1
        ISIG(NCHN,ISDE)=22
        ISIG(NCHN,3-ISDE)=24*KCHW
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACZW*WIDS(24,(5-KCHW)/2)
  690   CONTINUE
  700   CONTINUE
      ENDIF

      ELSEIF(ISUB.LE.80) THEN
      IF(ISUB.EQ.71) THEN
C...Z0 + Z0 -> Z0 + Z0.
        IF(SH.LE.4.01*SQMZ) GOTO 730

        IF(MSTP(46).LE.2) THEN
C...Exact scattering ME:s for on-mass-shell gauge bosons.
          BE2=1.-4.*SQMZ/SH
          TH=-0.5*SH*BE2*(1.-CTH)
          UH=-0.5*SH*BE2*(1.+CTH)
          IF(MAX(TH,UH).GT.-1.) GOTO 730
          SHANG=1./(1.-XW)*SQMW/SQMZ*(1.+BE2)**2
          ASHRE=(SH-SQMH)/((SH-SQMH)**2+GMMH**2)*SHANG
          ASHIM=-GMMH/((SH-SQMH)**2+GMMH**2)*SHANG
          THANG=1./(1.-XW)*SQMW/SQMZ*(BE2-CTH)**2
          ATHRE=(TH-SQMH)/((TH-SQMH)**2+GMMH**2)*THANG
          ATHIM=-GMMH/((TH-SQMH)**2+GMMH**2)*THANG
          UHANG=1./(1.-XW)*SQMW/SQMZ*(BE2+CTH)**2
          AUHRE=(UH-SQMH)/((UH-SQMH)**2+GMMH**2)*UHANG
          AUHIM=-GMMH/((UH-SQMH)**2+GMMH**2)*UHANG
          FACZZ=COMFAC*1./(4096.*PARU(1)**2*16.*(1.-XW)**2)*
     &    (AEM/XW)**4*(SH/SQMW)**2*(SQMZ/SQMW)*SH2
          IF(MSTP(46).LE.0) FACZZ=FACZZ*(ASHRE**2+ASHIM**2)
          IF(MSTP(46).EQ.1) FACZZ=FACZZ*((ASHRE+ATHRE+AUHRE)**2+
     &    (ASHIM+ATHIM+AUHIM)**2)
          IF(MSTP(46).EQ.2) FACZZ=0.

        ELSE
C...Strongly interacting Z_L/W_L model of Dobado, Herrero, Terron.
          FACZZ=COMFAC*(AEM/(16.*PARU(1)*XW*(1.-XW)))**2*(64./9.)*
     &    ABS(A00U+2.*A20U)**2
        ENDIF
        FACZZ=FACZZ*WIDS(23,1)

        DO 720 I=MIN1,MAX1
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 720
        EI=KCHG(IABS(I),1)/3.
        AI=SIGN(1.,EI)
        VI=AI-4.*EI*XW
        AVI=AI**2+VI**2
        DO 710 J=MIN2,MAX2
        IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 710
        EJ=KCHG(IABS(J),1)/3.
        AJ=SIGN(1.,EJ)
        VJ=AJ-4.*EJ*XW
        AVJ=AJ**2+VJ**2
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=J
        ISIG(NCHN,3)=1
        SIGH(NCHN)=0.5*FACZZ*AVI*AVJ
  710   CONTINUE
  720   CONTINUE
  730   CONTINUE

      ELSEIF(ISUB.EQ.72) THEN
C...Z0 + Z0 -> W+ + W-.
        IF(SH.LE.4.01*SQMZ) GOTO 760

        IF(MSTP(46).LE.2) THEN
C...Exact scattering ME:s for on-mass-shell gauge bosons.
          BE2=SQRT((1.-4.*SQMW/SH)*(1.-4.*SQMZ/SH))
          CTH2=CTH**2
          TH=-0.5*SH*(1.-2.*(SQMW+SQMZ)/SH-BE2*CTH)
          UH=-0.5*SH*(1.-2.*(SQMW+SQMZ)/SH+BE2*CTH)
          IF(MAX(TH,UH).GT.-1.) GOTO 760
          SHANG=4.*SQRT(SQMW/(SQMZ*(1.-XW)))*(1.-2.*SQMW/SH)*
     &    (1.-2.*SQMZ/SH)
          ASHRE=(SH-SQMH)/((SH-SQMH)**2+GMMH**2)*SHANG
          ASHIM=-GMMH/((SH-SQMH)**2+GMMH**2)*SHANG
          ATWRE=(1.-XW)/SQMZ*SH/(TH-SQMW)*((CTH-BE2)**2*(3./2.+BE2/2.*
     &    CTH-(SQMW+SQMZ)/SH+(SQMW-SQMZ)**2/(SH*SQMW))+4.*((SQMW+SQMZ)/
     &    SH*(1.-3.*CTH2)+8.*SQMW*SQMZ/SH2*(2.*CTH2-1.)+
     &    4.*(SQMW**2+SQMZ**2)/SH2*CTH2+2.*(SQMW+SQMZ)/SH*BE2*CTH))
          ATWIM=0.
          AUWRE=(1.-XW)/SQMZ*SH/(UH-SQMW)*((CTH+BE2)**2*(3./2.-BE2/2.*
     &    CTH-(SQMW+SQMZ)/SH+(SQMW-SQMZ)**2/(SH*SQMW))+4.*((SQMW+SQMZ)/
     &    SH*(1.-3.*CTH2)+8.*SQMW*SQMZ/SH2*(2.*CTH2-1.)+
     &    4.*(SQMW**2+SQMZ**2)/SH2*CTH2-2.*(SQMW+SQMZ)/SH*BE2*CTH))
          AUWIM=0.
          A4RE=2.*(1.-XW)/SQMZ*(3.-CTH2-4.*(SQMW+SQMZ)/SH)
          A4IM=0.
          FACWW=COMFAC*1./(4096.*PARU(1)**2*16.*(1.-XW)**2)*
     &    (AEM/XW)**4*(SH/SQMW)**2*(SQMZ/SQMW)*SH2
          IF(MSTP(46).LE.0) FACWW=FACWW*(ASHRE**2+ASHIM**2)
          IF(MSTP(46).EQ.1) FACWW=FACWW*((ASHRE+ATWRE+AUWRE+A4RE)**2+
     &    (ASHIM+ATWIM+AUWIM+A4IM)**2)
          IF(MSTP(46).EQ.2) FACWW=FACWW*((ATWRE+AUWRE+A4RE)**2+
     &    (ATWIM+AUWIM+A4IM)**2)

        ELSE
C...Strongly interacting Z_L/W_L model of Dobado, Herrero, Terron.
          FACWW=COMFAC*(AEM/(16.*PARU(1)*XW*(1.-XW)))**2*(64./9.)*
     &    ABS(A00U-A20U)**2
        ENDIF
        FACWW=FACWW*WIDS(24,1)

        DO 750 I=MIN1,MAX1
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 750
        EI=KCHG(IABS(I),1)/3.
        AI=SIGN(1.,EI)
        VI=AI-4.*EI*XW
        AVI=AI**2+VI**2
        DO 740 J=MIN2,MAX2
        IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 740
        EJ=KCHG(IABS(J),1)/3.
        AJ=SIGN(1.,EJ)
        VJ=AJ-4.*EJ*XW
        AVJ=AJ**2+VJ**2
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=J
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACWW*AVI*AVJ
  740   CONTINUE
  750   CONTINUE
  760   CONTINUE

      ELSEIF(ISUB.EQ.73) THEN
C...Z0 + W+/- -> Z0 + W+/-.
        IF(SH.LE.2.*SQMZ+2.*SQMW) GOTO 790

        IF(MSTP(46).LE.2) THEN
C...Exact scattering ME:s for on-mass-shell gauge bosons.
          BE2=1.-2.*(SQMZ+SQMW)/SH+((SQMZ-SQMW)/SH)**2
          EP1=1.-(SQMZ-SQMW)/SH
          EP2=1.+(SQMZ-SQMW)/SH
          TH=-0.5*SH*BE2*(1.-CTH)
          UH=(SQMZ-SQMW)**2/SH-0.5*SH*BE2*(1.+CTH)
          IF(MAX(TH,UH).GT.-1.) GOTO 790
          THANG=(BE2-EP1*CTH)*(BE2-EP2*CTH)
          ATHRE=(TH-SQMH)/((TH-SQMH)**2+GMMH**2)*THANG
          ATHIM=-GMMH/((TH-SQMH)**2+GMMH**2)*THANG
          ASWRE=-(1.-XW)/SQMZ*SH/(SH-SQMW)*(-BE2*(EP1+EP2)**4*CTH+
     &    1./4.*(BE2+EP1*EP2)**2*((EP1-EP2)**2-4.*BE2*CTH)+
     &    2.*BE2*(BE2+EP1*EP2)*(EP1+EP2)**2*CTH-
     &    1./16.*SH/SQMW*(EP1**2-EP2**2)**2*(BE2+EP1*EP2)**2)
          ASWIM=0.
          AUWRE=(1.-XW)/SQMZ*SH/(UH-SQMW)*(-BE2*(EP2+EP1*CTH)*
     &    (EP1+EP2*CTH)*(BE2+EP1*EP2)+BE2*(EP2+EP1*CTH)*
     &    (BE2+EP1*EP2*CTH)*(2.*EP2-EP2*CTH+EP1)-BE2*(EP2+EP1*CTH)**2*
     &    (BE2-EP2**2*CTH)-1./8.*(BE2+EP1*EP2*CTH)**2*((EP1+EP2)**2+
     &    2.*BE2*(1.-CTH))+1./32.*SH/SQMW*(BE2+EP1*EP2*CTH)**2*
     &    (EP1**2-EP2**2)**2-BE2*(EP1+EP2*CTH)*(EP2+EP1*CTH)*
     &    (BE2+EP1*EP2)+BE2*(EP1+EP2*CTH)*(BE2+EP1*EP2*CTH)*
     &    (2.*EP1-EP1*CTH+EP2)-BE2*(EP1+EP2*CTH)**2*(BE2-EP1**2*CTH)-
     &    1./8.*(BE2+EP1*EP2*CTH)**2*((EP1+EP2)**2+2.*BE2*(1.-CTH))+
     &    1./32.*SH/SQMW*(BE2+EP1*EP2*CTH)**2*(EP1**2-EP2**2)**2)
          AUWIM=0.
          A4RE=(1.-XW)/SQMZ*(EP1**2*EP2**2*(CTH**2-1.)-
     &    2.*BE2*(EP1**2+EP2**2+EP1*EP2)*CTH-2.*BE2*EP1*EP2)
          A4IM=0.
          FACZW=COMFAC*1./(4096.*PARU(1)**2*4.*(1.-XW))*(AEM/XW)**4*
     &    (SH/SQMW)**2*SQRT(SQMZ/SQMW)*SH2
          IF(MSTP(46).LE.0) FACZW=0.
          IF(MSTP(46).EQ.1) FACZW=FACZW*((ATHRE+ASWRE+AUWRE+A4RE)**2+
     &    (ATHIM+ASWIM+AUWIM+A4IM)**2)
          IF(MSTP(46).EQ.2) FACZW=FACZW*((ASWRE+AUWRE+A4RE)**2+
     &    (ASWIM+AUWIM+A4IM)**2)

        ELSE
C...Strongly interacting Z_L/W_L model of Dobado, Herrero, Terron.
          FACZW=COMFAC*AEM**2/(64.*PARU(1)**2*XW**2*(1.-XW))*16.*
     &    ABS(A20U+3.*A11U*CTH)**2
        ENDIF
        FACZW=FACZW*WIDS(23,2)

        DO 780 I=MIN1,MAX1
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 780
        EI=KCHG(IABS(I),1)/3.
        AI=SIGN(1.,EI)
        VI=AI-4.*EI*XW
        AVI=AI**2+VI**2
        KCHWI=ISIGN(1,KCHG(IABS(I),1)*ISIGN(1,I))
        DO 770 J=MIN2,MAX2
        IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 770
        EJ=KCHG(IABS(J),1)/3.
        AJ=SIGN(1.,EJ)
        VJ=AI-4.*EJ*XW
        AVJ=AJ**2+VJ**2
        KCHWJ=ISIGN(1,KCHG(IABS(J),1)*ISIGN(1,J))
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=J
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACZW*(AVI*VINT(180+J)*WIDS(24,(5-KCHWJ)/2)+
     &  VINT(180+I)*WIDS(24,(5-KCHWI)/2)*AVJ)
  770   CONTINUE
  780   CONTINUE
  790   CONTINUE

      ELSEIF(ISUB.EQ.75) THEN
C...W+ + W- -> gamma + gamma.

      ELSEIF(ISUB.EQ.76) THEN
C...W+ + W- -> Z0 + Z0.
        IF(SH.LE.4.01*SQMZ) GOTO 820

        IF(MSTP(46).LE.2) THEN
C...Exact scattering ME:s for on-mass-shell gauge bosons.
          BE2=SQRT((1.-4.*SQMW/SH)*(1.-4.*SQMZ/SH))
          CTH2=CTH**2
          TH=-0.5*SH*(1.-2.*(SQMW+SQMZ)/SH-BE2*CTH)
          UH=-0.5*SH*(1.-2.*(SQMW+SQMZ)/SH+BE2*CTH)
          IF(MAX(TH,UH).GT.-1.) GOTO 820
          SHANG=4.*SQRT(SQMW/(SQMZ*(1.-XW)))*(1.-2.*SQMW/SH)*
     &    (1.-2.*SQMZ/SH)
          ASHRE=(SH-SQMH)/((SH-SQMH)**2+GMMH**2)*SHANG
          ASHIM=-GMMH/((SH-SQMH)**2+GMMH**2)*SHANG
          ATWRE=(1.-XW)/SQMZ*SH/(TH-SQMW)*((CTH-BE2)**2*(3./2.+BE2/2.*
     &    CTH-(SQMW+SQMZ)/SH+(SQMW-SQMZ)**2/(SH*SQMW))+4.*((SQMW+SQMZ)/
     &    SH*(1.-3.*CTH2)+8.*SQMW*SQMZ/SH2*(2.*CTH2-1.)+
     &    4.*(SQMW**2+SQMZ**2)/SH2*CTH2+2.*(SQMW+SQMZ)/SH*BE2*CTH))
          ATWIM=0.
          AUWRE=(1.-XW)/SQMZ*SH/(UH-SQMW)*((CTH+BE2)**2*(3./2.-BE2/2.*
     &    CTH-(SQMW+SQMZ)/SH+(SQMW-SQMZ)**2/(SH*SQMW))+4.*((SQMW+SQMZ)/
     &    SH*(1.-3.*CTH2)+8.*SQMW*SQMZ/SH2*(2.*CTH2-1.)+
     &    4.*(SQMW**2+SQMZ**2)/SH2*CTH2-2.*(SQMW+SQMZ)/SH*BE2*CTH))
          AUWIM=0.
          A4RE=2.*(1.-XW)/SQMZ*(3.-CTH2-4.*(SQMW+SQMZ)/SH)
          A4IM=0.
          FACZZ=COMFAC*1./(4096.*PARU(1)**2)*(AEM/XW)**4*
     &    (SH/SQMW)**2*SH2
          IF(MSTP(46).LE.0) FACZZ=FACZZ*(ASHRE**2+ASHIM**2)
          IF(MSTP(46).EQ.1) FACZZ=FACZZ*((ASHRE+ATWRE+AUWRE+A4RE)**2+
     &    (ASHIM+ATWIM+AUWIM+A4IM)**2)
          IF(MSTP(46).EQ.2) FACZZ=FACZZ*((ATWRE+AUWRE+A4RE)**2+
     &    (ATWIM+AUWIM+A4IM)**2)

        ELSE
C...Strongly interacting Z_L/W_L model of Dobado, Herrero, Terron.
          FACZZ=COMFAC*(AEM/(4.*PARU(1)*XW))**2*(64./9.)*
     &    ABS(A00U-A20U)**2
        ENDIF
        FACZZ=FACZZ*WIDS(23,1)

        DO 810 I=MIN1,MAX1
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 810
        EI=SIGN(1.,FLOAT(I))*KCHG(IABS(I),1)
        DO 800 J=MIN2,MAX2
        IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 800
        EJ=SIGN(1.,FLOAT(J))*KCHG(IABS(J),1)
        IF(EI*EJ.GT.0.) GOTO 800
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=J
        ISIG(NCHN,3)=1
        SIGH(NCHN)=0.5*FACZZ*VINT(180+I)*VINT(180+J)
  800   CONTINUE
  810   CONTINUE
  820   CONTINUE

      ELSEIF(ISUB.EQ.77) THEN
C...W+/- + W+/- -> W+/- + W+/-.
        IF(SH.LE.4.01*SQMW) GOTO 850

        IF(MSTP(46).LE.2) THEN
C...Exact scattering ME:s for on-mass-shell gauge bosons.
          BE2=1.-4.*SQMW/SH
          BE4=BE2**2
          CTH2=CTH**2
          CTH3=CTH**3
          TH=-0.5*SH*BE2*(1.-CTH)
          UH=-0.5*SH*BE2*(1.+CTH)
          IF(MAX(TH,UH).GT.-1.) GOTO 850
          SHANG=(1.+BE2)**2
          ASHRE=(SH-SQMH)/((SH-SQMH)**2+GMMH**2)*SHANG
          ASHIM=-GMMH/((SH-SQMH)**2+GMMH**2)*SHANG
          THANG=(BE2-CTH)**2
          ATHRE=(TH-SQMH)/((TH-SQMH)**2+GMMH**2)*THANG
          ATHIM=-GMMH/((TH-SQMH)**2+GMMH**2)*THANG
          UHANG=(BE2+CTH)**2
          AUHRE=(UH-SQMH)/((UH-SQMH)**2+GMMH**2)*UHANG
          AUHIM=-GMMH/((UH-SQMH)**2+GMMH**2)*UHANG
          SGZANG=1./SQMW*BE2*(3.-BE2)**2*CTH
          ASGRE=XW*SGZANG
          ASGIM=0.
          ASZRE=(1.-XW)*SH/(SH-SQMZ)*SGZANG
          ASZIM=0.
          TGZANG=1./SQMW*(BE2*(4.-2.*BE2+BE4)+BE2*(4.-10.*BE2+BE4)*CTH+
     &    (2.-11.*BE2+10.*BE4)*CTH2+BE2*CTH3)
          ATGRE=0.5*XW*SH/TH*TGZANG
          ATGIM=0.
          ATZRE=0.5*(1.-XW)*SH/(TH-SQMZ)*TGZANG
          ATZIM=0.
          UGZANG=1./SQMW*(BE2*(4.-2.*BE2+BE4)-BE2*(4.-10.*BE2+BE4)*CTH+
     &    (2.-11.*BE2+10.*BE4)*CTH2-BE2*CTH3)
          AUGRE=0.5*XW*SH/UH*UGZANG
          AUGIM=0.
          AUZRE=0.5*(1.-XW)*SH/(UH-SQMZ)*UGZANG
          AUZIM=0.
          A4RE=1./SQMW*(1.+2.*BE2-6.*BE2*CTH-CTH2)
          A4IM=0.
          FWW=COMFAC*1./(4096.*PARU(1)**2)*(AEM/XW)**4*(SH/SQMW)**2*SH2
          IF(MSTP(46).LE.0) THEN
            AWWARE=ASHRE
            AWWAIM=ASHIM
            AWWSRE=0.
            AWWSIM=0.
          ELSEIF(MSTP(46).EQ.1) THEN
            AWWARE=ASHRE+ATHRE+ASGRE+ASZRE+ATGRE+ATZRE+A4RE
            AWWAIM=ASHIM+ATHIM+ASGIM+ASZIM+ATGIM+ATZIM+A4IM
            AWWSRE=ATHRE+AUHRE+ATGRE+ATZRE+AUGRE+AUZRE
            AWWSIM=ATHIM+AUHIM+ATGIM+ATZIM+AUGIM+AUZIM
          ELSE
            AWWARE=ASGRE+ASZRE+ATGRE+ATZRE+A4RE
            AWWAIM=ASGIM+ASZIM+ATGIM+ATZIM+A4IM
            AWWSRE=ATGRE+ATZRE+AUGRE+AUZRE
            AWWSIM=ATGIM+ATZIM+AUGIM+AUZIM
          ENDIF
          AWWA2=AWWARE**2+AWWAIM**2
          AWWS2=AWWSRE**2+AWWSIM**2

        ELSE
C...Strongly interacting Z_L/W_L model of Dobado, Herrero, Terron.
          FWWA=COMFAC*(AEM/(4.*PARU(1)*XW))**2*(64./9.)*
     &    ABS(A00U+0.5*A20U+4.5*A11U*CTH)**2
          FWWS=COMFAC*(AEM/(4.*PARU(1)*XW))**2*64.*ABS(A20U)**2
        ENDIF

        DO 840 I=MIN1,MAX1
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 840
        EI=SIGN(1.,FLOAT(I))*KCHG(IABS(I),1)
        DO 830 J=MIN2,MAX2
        IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 830
        EJ=SIGN(1.,FLOAT(J))*KCHG(IABS(J),1)
        IF(EI*EJ.LT.0.) THEN
C...W+W-
          IF(MSTP(45).EQ.1) GOTO 830
          IF(MSTP(46).LE.2) FACWW=FWW*AWWA2*WIDS(24,1)
          IF(MSTP(46).GE.3) FACWW=FWWA*WIDS(24,1)
        ELSE
C...W+W+/W-W-
          IF(MSTP(45).EQ.2) GOTO 830
          IF(MSTP(46).LE.2) FACWW=FWW*AWWS2
          IF(MSTP(46).GE.3) FACWW=FWWS
          IF(EI.GT.0.) FACWW=FACWW*VINT(91)
          IF(EI.LT.0.) FACWW=FACWW*VINT(92)
        ENDIF
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=J
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACWW*VINT(180+I)*VINT(180+J)
        IF(EI*EJ.GT.0.) SIGH(NCHN)=0.5*SIGH(NCHN)
  830   CONTINUE
  840   CONTINUE
  850   CONTINUE

      ELSEIF(ISUB.EQ.78) THEN
C...W+/- + H0 -> W+/- + H0.

      ELSEIF(ISUB.EQ.79) THEN
C...H0 + H0 -> H0 + H0.

      ENDIF

C...C: 2 -> 2, tree diagrams with masses.

      ELSEIF(ISUB.LE.90) THEN
      IF(ISUB.EQ.81) THEN
C...q + q~ -> Q + Q~.
        FACQQB=COMFAC*AS**2*4./9.*(((TH-SQM3)**2+
     &  (UH-SQM3)**2)/SH2+2.*SQM3/SH)
        IF(MSTP(35).GE.1) THEN
          IF(MSTP(35).EQ.1) THEN
            ALSSG=PARP(35)
          ELSE
            MST115=MSTU(115)
            MSTU(115)=MSTP(36)
            Q2BN=SQRT(SQM3*((SQRT(SH)-2.*SQRT(SQM3))**2+PARP(36)**2))
            ALSSG=ULALPS(Q2BN)
            MSTU(115)=MST115
          ENDIF
          XREPU=PARU(1)*ALSSG/(6.*SQRT(MAX(1E-20,1.-4.*SQM3/SH)))
          FREPU=XREPU/(EXP(MIN(87.,XREPU))-1.)
          VINT(138)=FREPU
          FACQQB=FACQQB*FREPU
        ENDIF
        DO 860 I=MINA,MAXA
        IF(I.EQ.0.OR.IABS(I).GT.MSTP(54).OR.
     &  KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 860
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=-I
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACQQB
  860   CONTINUE

      ELSEIF(ISUB.EQ.82) THEN
C...g + g -> Q + Q~.
        FACQQ1=COMFAC*FACA*AS**2*1./6.*((UH-SQM3)/(TH-SQM3)-
     &  2.*(UH-SQM3)**2/SH2+4.*SQM3/SH*(TH*UH-SQM3**2)/(TH-SQM3)**2)
        FACQQ2=COMFAC*FACA*AS**2*1./6.*((TH-SQM3)/(UH-SQM3)-
     &  2.*(TH-SQM3)**2/SH2+4.*SQM3/SH*(TH*UH-SQM3**2)/(UH-SQM3)**2)
        IF(MSTP(35).GE.1) THEN
          IF(MSTP(35).EQ.1) THEN
            ALSSG=PARP(35)
          ELSE
            MST115=MSTU(115)
            MSTU(115)=MSTP(36)
            Q2BN=SQRT(SQM3*((SQRT(SH)-2.*SQRT(SQM3))**2+PARP(36)**2))
            ALSSG=ULALPS(Q2BN)
            MSTU(115)=MST115
          ENDIF
          XATTR=4.*PARU(1)*ALSSG/(3.*SQRT(MAX(1E-20,1.-4.*SQM3/SH)))
          FATTR=XATTR/(1.-EXP(-MIN(87.,XATTR)))
          XREPU=PARU(1)*ALSSG/(6.*SQRT(MAX(1E-20,1.-4.*SQM3/SH)))
          FREPU=XREPU/(EXP(MIN(87.,XREPU))-1.)
          FATRE=(2.*FATTR+5.*FREPU)/7.
          VINT(138)=FATRE
          FACQQ1=FACQQ1*FATRE
          FACQQ2=FACQQ2*FATRE
        ENDIF
        IF(KFAC(1,21)*KFAC(2,21).EQ.0) GOTO 870
        NCHN=NCHN+1
        ISIG(NCHN,1)=21
        ISIG(NCHN,2)=21
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACQQ1
        NCHN=NCHN+1
        ISIG(NCHN,1)=21
        ISIG(NCHN,2)=21
        ISIG(NCHN,3)=2
        SIGH(NCHN)=FACQQ2
  870   CONTINUE

      ELSEIF(ISUB.EQ.83) THEN
C...f + q -> f' + Q.
        FACQQS=COMFAC*(0.5*AEM/XW)**2*SH*(SH-SQM3)/(SQMW-TH)**2
        FACQQU=COMFAC*(0.5*AEM/XW)**2*UH*(UH-SQM3)/(SQMW-TH)**2
        DO 890 I=MIN1,MAX1
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 890
        DO 880 J=MIN2,MAX2
        IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 880
        IF(I*J.GT.0.AND.MOD(IABS(I+J),2).EQ.0) GOTO 880
        IF(I*J.LT.0.AND.MOD(IABS(I+J),2).EQ.1) GOTO 880
        IF(IABS(I).LT.MINT(55).AND.MOD(IABS(I+MINT(55)),2).EQ.1) THEN
          NCHN=NCHN+1
          ISIG(NCHN,1)=I
          ISIG(NCHN,2)=J
          ISIG(NCHN,3)=1
          IF(MOD(MINT(55),2).EQ.0) FACCKM=VCKM(MINT(55)/2,
     &    (IABS(I)+1)/2)*VINT(180+J)
          IF(MOD(MINT(55),2).EQ.1) FACCKM=VCKM(IABS(I)/2,
     &    (MINT(55)+1)/2)*VINT(180+J)
          IF(I*J.GT.0) SIGH(NCHN)=FACQQS*FACCKM
          IF(I*J.LT.0) SIGH(NCHN)=FACQQU*FACCKM
        ENDIF
        IF(IABS(J).LT.MINT(55).AND.MOD(IABS(J+MINT(55)),2).EQ.1) THEN
          NCHN=NCHN+1
          ISIG(NCHN,1)=I
          ISIG(NCHN,2)=J
          ISIG(NCHN,3)=2
          IF(MOD(MINT(55),2).EQ.0) FACCKM=VCKM(MINT(55)/2,
     &    (IABS(J)+1)/2)*VINT(180+I)
          IF(MOD(MINT(55),2).EQ.1) FACCKM=VCKM(IABS(J)/2,
     &    (MINT(55)+1)/2)*VINT(180+I)
          IF(I*J.GT.0) SIGH(NCHN)=FACQQS*FACCKM
          IF(I*J.LT.0) SIGH(NCHN)=FACQQU*FACCKM
        ENDIF
  880   CONTINUE
  890   CONTINUE

      ELSEIF(ISUB.EQ.84) THEN
C...g + gamma -> Q + Q~.
        FMTU=SQM3/(SQM3-TH)+SQM3/(SQM3-UH)
        FACQQ=COMFAC*AS*AEM*(KCHG(IABS(MINT(55)),1)/3.)**2*
     &  ((SQM3-TH)/(SQM3-UH)+(SQM3-UH)/(SQM3-TH)+4.*FMTU*(1.-FMTU))
        IF(MSTP(35).GE.1) THEN
          IF(MSTP(35).EQ.1) THEN
            ALSSG=PARP(35)
          ELSE
            MST115=MSTU(115)
            MSTU(115)=MSTP(36)
            Q2BN=SQRT(SQM3*((SQRT(SH)-2.*SQRT(SQM3))**2+PARP(36)**2))
            ALSSG=ULALPS(Q2BN)
            MSTU(115)=MST115
          ENDIF
          XREPU=PARU(1)*ALSSG/(6.*SQRT(MAX(1E-20,1.-4.*SQM3/SH)))
          FREPU=XREPU/(EXP(MIN(87.,XREPU))-1.)
          VINT(138)=FREPU
          FACQQ=FACQQ*FREPU
        ENDIF
        IF(KFAC(1,21)*KFAC(2,22).NE.0) THEN
          NCHN=NCHN+1
          ISIG(NCHN,1)=21
          ISIG(NCHN,2)=22
          ISIG(NCHN,3)=1
          SIGH(NCHN)=FACQQ
        ENDIF
        IF(KFAC(1,22)*KFAC(2,21).NE.0) THEN
          NCHN=NCHN+1
          ISIG(NCHN,1)=22
          ISIG(NCHN,2)=21
          ISIG(NCHN,3)=1
          SIGH(NCHN)=FACQQ2
        ENDIF

      ELSEIF(ISUB.EQ.85) THEN
C...gamma + gamma -> F + F~ (heavy fermion, quark or lepton).
        FMTU=SQM3/(SQM3-TH)+SQM3/(SQM3-UH)
        FACFF=COMFAC*AEM**2*(KCHG(IABS(MINT(56)),1)/3.)**4*2.*
     &  ((SQM3-TH)/(SQM3-UH)+(SQM3-UH)/(SQM3-TH)+4.*FMTU*(1.-FMTU))
        IF(IABS(MINT(56)).LT.10) FACFF=3.*FACFF
        IF(IABS(MINT(56)).LT.10.AND.MSTP(35).GE.1) THEN
          IF(MSTP(35).EQ.1) THEN
            ALSSG=PARP(35)
          ELSE
            MST115=MSTU(115)
            MSTU(115)=MSTP(36)
            Q2BN=SQRT(SQM3*((SQRT(SH)-2.*SQRT(SQM3))**2+PARP(36)**2))
            ALSSG=ULALPS(Q2BN)
            MSTU(115)=MST115
          ENDIF
          XATTR=4.*PARU(1)*ALSSG/(3.*SQRT(MAX(1E-20,1.-4.*SQM3/SH)))
          FATTR=XATTR/(1.-EXP(-MIN(87.,XATTR)))
          VINT(138)=FATTR
          FACFF=FACFF*FATTR
        ENDIF
        IF(KFAC(1,22)*KFAC(2,22).NE.0) THEN
          NCHN=NCHN+1
          ISIG(NCHN,1)=22
          ISIG(NCHN,2)=22
          ISIG(NCHN,3)=1
          SIGH(NCHN)=FACFF
        ENDIF
      ENDIF

C...D: Mimimum bias processes.

      ELSEIF(ISUB.LE.100) THEN
      IF(ISUB.EQ.91) THEN
C...Elastic scattering.
        SIGS=XSEC(ISUB,1)

      ELSEIF(ISUB.EQ.92) THEN
C...Single diffractive scattering.
        SIGS=XSEC(ISUB,1)

      ELSEIF(ISUB.EQ.93) THEN
C...Double diffractive scattering.
        SIGS=XSEC(ISUB,1)

      ELSEIF(ISUB.EQ.94) THEN
C...Central diffractive scattering.
        SIGS=XSEC(ISUB,1)

      ELSEIF(ISUB.EQ.95) THEN
C...Low-pT scattering.
        SIGS=XSEC(ISUB,1)

      ELSEIF(ISUB.EQ.96) THEN
C...Multiple interactions: sum of QCD processes.
        CALL RYWIDT(21,SH,WDTP,WDTE)

C...q + q' -> q + q'.
        FACQQ1=COMFAC*AS**2*4./9.*(SH2+UH2)/TH2
        FACQQB=COMFAC*AS**2*4./9.*((SH2+UH2)/TH2*FACA-
     &  MSTP(34)*2./3.*UH2/(SH*TH))
        FACQQ2=COMFAC*AS**2*4./9.*((SH2+TH2)/UH2-
     &  MSTP(34)*2./3.*SH2/(TH*UH))
        DO 910 I=-3,3
        IF(I.EQ.0) GOTO 910
        DO 900 J=-3,3
        IF(J.EQ.0) GOTO 900
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=J
        ISIG(NCHN,3)=111
        SIGH(NCHN)=FACQQ1
        IF(I.EQ.-J) SIGH(NCHN)=FACQQB
        IF(I.EQ.J) THEN
          SIGH(NCHN)=0.5*SIGH(NCHN)
          NCHN=NCHN+1
          ISIG(NCHN,1)=I
          ISIG(NCHN,2)=J
          ISIG(NCHN,3)=112
          SIGH(NCHN)=0.5*FACQQ2
        ENDIF
  900   CONTINUE
  910   CONTINUE

C...q + q~ -> q' + q~' or g + g.
        FACQQB=COMFAC*AS**2*4./9.*(TH2+UH2)/SH2*(WDTE(0,1)+WDTE(0,2)+
     &  WDTE(0,3)+WDTE(0,4))
        FACGG1=COMFAC*AS**2*32./27.*(UH/TH-(2.+MSTP(34)*1./4.)*UH2/SH2)
        FACGG2=COMFAC*AS**2*32./27.*(TH/UH-(2.+MSTP(34)*1./4.)*TH2/SH2)
        DO 920 I=-3,3
        IF(I.EQ.0) GOTO 920
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=-I
        ISIG(NCHN,3)=121
        SIGH(NCHN)=FACQQB
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=-I
        ISIG(NCHN,3)=131
        SIGH(NCHN)=0.5*FACGG1
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=-I
        ISIG(NCHN,3)=132
        SIGH(NCHN)=0.5*FACGG2
  920   CONTINUE

C...q + g -> q + g.
        FACQG1=COMFAC*AS**2*4./9.*((2.+MSTP(34)*1./4.)*UH2/TH2-UH/SH)*
     &  FACA
        FACQG2=COMFAC*AS**2*4./9.*((2.+MSTP(34)*1./4.)*SH2/TH2-SH/UH)
        DO 940 I=-3,3
        IF(I.EQ.0) GOTO 940
        DO 930 ISDE=1,2
        NCHN=NCHN+1
        ISIG(NCHN,ISDE)=I
        ISIG(NCHN,3-ISDE)=21
        ISIG(NCHN,3)=281
        SIGH(NCHN)=FACQG1
        NCHN=NCHN+1
        ISIG(NCHN,ISDE)=I
        ISIG(NCHN,3-ISDE)=21
        ISIG(NCHN,3)=282
        SIGH(NCHN)=FACQG2
  930   CONTINUE
  940   CONTINUE

C...g + g -> q + q~ or g + g.
        FACQQ1=COMFAC*AS**2*1./6.*(UH/TH-(2.+MSTP(34)*1./4.)*UH2/SH2)*
     &  (WDTE(0,1)+WDTE(0,2)+WDTE(0,3)+WDTE(0,4))*FACA
        FACQQ2=COMFAC*AS**2*1./6.*(TH/UH-(2.+MSTP(34)*1./4.)*TH2/SH2)*
     &  (WDTE(0,1)+WDTE(0,2)+WDTE(0,3)+WDTE(0,4))*FACA
        FACGG1=COMFAC*AS**2*9./4.*(SH2/TH2+2.*SH/TH+3.+2.*TH/SH+
     &  TH2/SH2)*FACA
        FACGG2=COMFAC*AS**2*9./4.*(UH2/SH2+2.*UH/SH+3.+2.*SH/UH+
     &  SH2/UH2)*FACA
        FACGG3=COMFAC*AS**2*9./4.*(TH2/UH2+2.*TH/UH+3+2.*UH/TH+UH2/TH2)
        NCHN=NCHN+1
        ISIG(NCHN,1)=21
        ISIG(NCHN,2)=21
        ISIG(NCHN,3)=531
        SIGH(NCHN)=FACQQ1
        NCHN=NCHN+1
        ISIG(NCHN,1)=21
        ISIG(NCHN,2)=21
        ISIG(NCHN,3)=532
        SIGH(NCHN)=FACQQ2
        NCHN=NCHN+1
        ISIG(NCHN,1)=21
        ISIG(NCHN,2)=21
        ISIG(NCHN,3)=681
        SIGH(NCHN)=0.5*FACGG1
        NCHN=NCHN+1
        ISIG(NCHN,1)=21
        ISIG(NCHN,2)=21
        ISIG(NCHN,3)=682
        SIGH(NCHN)=0.5*FACGG2
        NCHN=NCHN+1
        ISIG(NCHN,1)=21
        ISIG(NCHN,2)=21
        ISIG(NCHN,3)=683
        SIGH(NCHN)=0.5*FACGG3
      ENDIF

C...E: 2 -> 1, loop diagrams.

      ELSEIF(ISUB.LE.110) THEN
      IF(ISUB.EQ.101) THEN
C...g + g -> gamma*/Z0.

      ELSEIF(ISUB.EQ.102) THEN
C...g + g -> H0 (or H'0, or A0).
        CALL RYWIDT(KFHIGG,SH,WDTP,WDTE)
        HP=AEM/(8.*XW)*SH/SQMW*SH
        HS=HP*WDTP(0)
        HF=HP*(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))
        FACBW=4.*COMFAC/((SH-SQMH)**2+HS**2)
        IF(ABS(SH-SQMH).GT.100.*HS) FACBW=0.
        HI=HP*WDTP(13)/32.
        IF(KFAC(1,21)*KFAC(2,21).EQ.0) GOTO 950
        NCHN=NCHN+1
        ISIG(NCHN,1)=21
        ISIG(NCHN,2)=21
        ISIG(NCHN,3)=1
        SIGH(NCHN)=HI*FACBW*HF
  950   CONTINUE

      ELSEIF(ISUB.EQ.103) THEN
C...gamma + gamma -> H0 (or H'0, or A0).
        CALL RYWIDT(KFHIGG,SH,WDTP,WDTE)
        HP=AEM/(8.*XW)*SH/SQMW*SH
        HS=HP*WDTP(0)
        HF=HP*(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))
        FACBW=4.*COMFAC/((SH-SQMH)**2+HS**2)
        IF(ABS(SH-SQMH).GT.100.*HS) FACBW=0.
        HI=HP*WDTP(14)*2.
        IF(KFAC(1,22)*KFAC(2,22).EQ.0) GOTO 960
        NCHN=NCHN+1
        ISIG(NCHN,1)=22
        ISIG(NCHN,2)=22
        ISIG(NCHN,3)=1
        SIGH(NCHN)=HI*FACBW*HF
  960   CONTINUE

      ENDIF

C...F: 2 -> 2, box diagrams.

      ELSEIF(ISUB.LE.120) THEN
      IF(ISUB.EQ.111) THEN
C...f + f~ -> g + H0 (q + q~ -> g + H0 only).
        A5STUR=0.
        A5STUI=0.
        DO 970 I=1,2*MSTP(1)
        SQMQ=PMAS(I,1)**2
        EPSS=4.*SQMQ/SH
        EPSH=4.*SQMQ/SQMH
        CALL RYWAUX(1,EPSS,W1SR,W1SI)
        CALL RYWAUX(1,EPSH,W1HR,W1HI)
        CALL RYWAUX(2,EPSS,W2SR,W2SI)
        CALL RYWAUX(2,EPSH,W2HR,W2HI)
        A5STUR=A5STUR+EPSH*(1.+SH/(TH+UH)*(W1SR-W1HR)+
     &  (0.25-SQMQ/(TH+UH))*(W2SR-W2HR))
        A5STUI=A5STUI+EPSH*(SH/(TH+UH)*(W1SI-W1HI)+
     &  (0.25-SQMQ/(TH+UH))*(W2SI-W2HI))
  970   CONTINUE
        FACGH=COMFAC*FACA/(144.*PARU(1)**2)*AEM/XW*AS**3*SQMH/SQMW*
     &  SQMH/SH*(UH**2+TH**2)/(UH+TH)**2*(A5STUR**2+A5STUI**2)
        FACGH=FACGH*WIDS(25,2)
        DO 980 I=MINA,MAXA
        IF(I.EQ.0.OR.IABS(I).GT.MSTP(54).OR.
     &  KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 980
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=-I
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACGH
  980   CONTINUE

      ELSEIF(ISUB.EQ.112) THEN
C...f + g -> f + H0 (q + g -> q + H0 only).
        A5TSUR=0.
        A5TSUI=0.
        DO 990 I=1,2*MSTP(1)
        SQMQ=PMAS(I,1)**2
        EPST=4.*SQMQ/TH
        EPSH=4.*SQMQ/SQMH
        CALL RYWAUX(1,EPST,W1TR,W1TI)
        CALL RYWAUX(1,EPSH,W1HR,W1HI)
        CALL RYWAUX(2,EPST,W2TR,W2TI)
        CALL RYWAUX(2,EPSH,W2HR,W2HI)
        A5TSUR=A5TSUR+EPSH*(1.+TH/(SH+UH)*(W1TR-W1HR)+
     &  (0.25-SQMQ/(SH+UH))*(W2TR-W2HR))
        A5TSUI=A5TSUI+EPSH*(TH/(SH+UH)*(W1TI-W1HI)+
     &  (0.25-SQMQ/(SH+UH))*(W2TI-W2HI))
  990   CONTINUE
        FACQH=COMFAC*FACA/(384.*PARU(1)**2)*AEM/XW*AS**3*SQMH/SQMW*
     &  SQMH/(-TH)*(UH**2+SH**2)/(UH+SH)**2*(A5TSUR**2+A5TSUI**2)
        FACQH=FACQH*WIDS(25,2)
        DO 1010 I=MINA,MAXA
        IF(I.EQ.0.OR.IABS(I).GT.MSTP(54)) GOTO 1010
        DO 1000 ISDE=1,2
        IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,21).EQ.0) GOTO 1000
        IF(ISDE.EQ.2.AND.KFAC(1,21)*KFAC(2,I).EQ.0) GOTO 1000
        NCHN=NCHN+1
        ISIG(NCHN,ISDE)=I
        ISIG(NCHN,3-ISDE)=21
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACQH
 1000   CONTINUE
 1010   CONTINUE

      ELSEIF(ISUB.EQ.113) THEN
C...g + g -> g + H0.
        A2STUR=0.
        A2STUI=0.
        A2USTR=0.
        A2USTI=0.
        A2TUSR=0.
        A2TUSI=0.
        A4STUR=0.
        A4STUI=0.
        DO 1020 I=1,2*MSTP(1)
        SQMQ=PMAS(I,1)**2
        EPSS=4.*SQMQ/SH
        EPST=4.*SQMQ/TH
        EPSU=4.*SQMQ/UH
        EPSH=4.*SQMQ/SQMH
        IF(EPSH.LT.1.E-6) GOTO 1020
        CALL RYWAUX(1,EPSS,W1SR,W1SI)
        CALL RYWAUX(1,EPST,W1TR,W1TI)
        CALL RYWAUX(1,EPSU,W1UR,W1UI)
        CALL RYWAUX(1,EPSH,W1HR,W1HI)
        CALL RYWAUX(2,EPSS,W2SR,W2SI)
        CALL RYWAUX(2,EPST,W2TR,W2TI)
        CALL RYWAUX(2,EPSU,W2UR,W2UI)
        CALL RYWAUX(2,EPSH,W2HR,W2HI)
        CALL RYI3AU(EPSS,TH/UH,Y3STUR,Y3STUI)
        CALL RYI3AU(EPSS,UH/TH,Y3SUTR,Y3SUTI)
        CALL RYI3AU(EPST,SH/UH,Y3TSUR,Y3TSUI)
        CALL RYI3AU(EPST,UH/SH,Y3TUSR,Y3TUSI)
        CALL RYI3AU(EPSU,SH/TH,Y3USTR,Y3USTI)
        CALL RYI3AU(EPSU,TH/SH,Y3UTSR,Y3UTSI)
        CALL RYI3AU(EPSH,SQMH/SH*TH/UH,YHSTUR,YHSTUI)
        CALL RYI3AU(EPSH,SQMH/SH*UH/TH,YHSUTR,YHSUTI)
        CALL RYI3AU(EPSH,SQMH/TH*SH/UH,YHTSUR,YHTSUI)
        CALL RYI3AU(EPSH,SQMH/TH*UH/SH,YHTUSR,YHTUSI)
        CALL RYI3AU(EPSH,SQMH/UH*SH/TH,YHUSTR,YHUSTI)
        CALL RYI3AU(EPSH,SQMH/UH*TH/SH,YHUTSR,YHUTSI)
        W3STUR=YHSTUR-Y3STUR-Y3UTSR
        W3STUI=YHSTUI-Y3STUI-Y3UTSI
        W3SUTR=YHSUTR-Y3SUTR-Y3TUSR
        W3SUTI=YHSUTI-Y3SUTI-Y3TUSI
        W3TSUR=YHTSUR-Y3TSUR-Y3USTR
        W3TSUI=YHTSUI-Y3TSUI-Y3USTI
        W3TUSR=YHTUSR-Y3TUSR-Y3SUTR
        W3TUSI=YHTUSI-Y3TUSI-Y3SUTI
        W3USTR=YHUSTR-Y3USTR-Y3TSUR
        W3USTI=YHUSTI-Y3USTI-Y3TSUI
        W3UTSR=YHUTSR-Y3UTSR-Y3STUR
        W3UTSI=YHUTSI-Y3UTSI-Y3STUI
        B2STUR=SQMQ/SQMH**2*(SH*(UH-SH)/(SH+UH)+2.*TH*UH*(UH+2.*SH)/
     &  (SH+UH)**2*(W1TR-W1HR)+(SQMQ-SH/4.)*(0.5*W2SR+0.5*W2HR-W2TR+
     &  W3STUR)+SH2*(2.*SQMQ/(SH+UH)**2-0.5/(SH+UH))*(W2TR-W2HR)+
     &  0.5*TH*UH/SH*(W2HR-2.*W2TR)+0.125*(SH-12.*SQMQ-4.*TH*UH/SH)*
     &  W3TSUR)
        B2STUI=SQMQ/SQMH**2*(2.*TH*UH*(UH+2.*SH)/(SH+UH)**2*
     &  (W1TI-W1HI)+(SQMQ-SH/4.)*(0.5*W2SI+0.5*W2HI-W2TI+W3STUI)+
     &  SH2*(2.*SQMQ/(SH+UH)**2-0.5/(SH+UH))*(W2TI-W2HI)+0.5*TH*UH/SH*
     &  (W2HI-2.*W2TI)+0.125*(SH-12.*SQMQ-4.*TH*UH/SH)*W3TSUI)
        B2SUTR=SQMQ/SQMH**2*(SH*(TH-SH)/(SH+TH)+2.*UH*TH*(TH+2.*SH)/
     &  (SH+TH)**2*(W1UR-W1HR)+(SQMQ-SH/4.)*(0.5*W2SR+0.5*W2HR-W2UR+
     &  W3SUTR)+SH2*(2.*SQMQ/(SH+TH)**2-0.5/(SH+TH))*(W2UR-W2HR)+
     &  0.5*UH*TH/SH*(W2HR-2.*W2UR)+0.125*(SH-12.*SQMQ-4.*UH*TH/SH)*
     &  W3USTR)
        B2SUTI=SQMQ/SQMH**2*(2.*UH*TH*(TH+2.*SH)/(SH+TH)**2*
     &  (W1UI-W1HI)+(SQMQ-SH/4.)*(0.5*W2SI+0.5*W2HI-W2UI+W3SUTI)+
     &  SH2*(2.*SQMQ/(SH+TH)**2-0.5/(SH+TH))*(W2UI-W2HI)+0.5*UH*TH/SH*
     &  (W2HI-2.*W2UI)+0.125*(SH-12.*SQMQ-4.*UH*TH/SH)*W3USTI)
        B2TSUR=SQMQ/SQMH**2*(TH*(UH-TH)/(TH+UH)+2.*SH*UH*(UH+2.*TH)/
     &  (TH+UH)**2*(W1SR-W1HR)+(SQMQ-TH/4.)*(0.5*W2TR+0.5*W2HR-W2SR+
     &  W3TSUR)+TH2*(2.*SQMQ/(TH+UH)**2-0.5/(TH+UH))*(W2SR-W2HR)+
     &  0.5*SH*UH/TH*(W2HR-2.*W2SR)+0.125*(TH-12.*SQMQ-4.*SH*UH/TH)*
     &  W3STUR)
        B2TSUI=SQMQ/SQMH**2*(2.*SH*UH*(UH+2.*TH)/(TH+UH)**2*
     &  (W1SI-W1HI)+(SQMQ-TH/4.)*(0.5*W2TI+0.5*W2HI-W2SI+W3TSUI)+
     &  TH2*(2.*SQMQ/(TH+UH)**2-0.5/(TH+UH))*(W2SI-W2HI)+0.5*SH*UH/TH*
     &  (W2HI-2.*W2SI)+0.125*(TH-12.*SQMQ-4.*SH*UH/TH)*W3STUI)
        B2TUSR=SQMQ/SQMH**2*(TH*(SH-TH)/(TH+SH)+2.*UH*SH*(SH+2.*TH)/
     &  (TH+SH)**2*(W1UR-W1HR)+(SQMQ-TH/4.)*(0.5*W2TR+0.5*W2HR-W2UR+
     &  W3TUSR)+TH2*(2.*SQMQ/(TH+SH)**2-0.5/(TH+SH))*(W2UR-W2HR)+
     &  0.5*UH*SH/TH*(W2HR-2.*W2UR)+0.125*(TH-12.*SQMQ-4.*UH*SH/TH)*
     &  W3UTSR)
        B2TUSI=SQMQ/SQMH**2*(2.*UH*SH*(SH+2.*TH)/(TH+SH)**2*
     &  (W1UI-W1HI)+(SQMQ-TH/4.)*(0.5*W2TI+0.5*W2HI-W2UI+W3TUSI)+
     &  TH2*(2.*SQMQ/(TH+SH)**2-0.5/(TH+SH))*(W2UI-W2HI)+0.5*UH*SH/TH*
     &  (W2HI-2.*W2UI)+0.125*(TH-12.*SQMQ-4.*UH*SH/TH)*W3UTSI)
        B2USTR=SQMQ/SQMH**2*(UH*(TH-UH)/(UH+TH)+2.*SH*TH*(TH+2.*UH)/
     &  (UH+TH)**2*(W1SR-W1HR)+(SQMQ-UH/4.)*(0.5*W2UR+0.5*W2HR-W2SR+
     &  W3USTR)+UH2*(2.*SQMQ/(UH+TH)**2-0.5/(UH+TH))*(W2SR-W2HR)+
     &  0.5*SH*TH/UH*(W2HR-2.*W2SR)+0.125*(UH-12.*SQMQ-4.*SH*TH/UH)*
     &  W3SUTR)
        B2USTI=SQMQ/SQMH**2*(2.*SH*TH*(TH+2.*UH)/(UH+TH)**2*
     &  (W1SI-W1HI)+(SQMQ-UH/4.)*(0.5*W2UI+0.5*W2HI-W2SI+W3USTI)+
     &  UH2*(2.*SQMQ/(UH+TH)**2-0.5/(UH+TH))*(W2SI-W2HI)+0.5*SH*TH/UH*
     &  (W2HI-2.*W2SI)+0.125*(UH-12.*SQMQ-4.*SH*TH/UH)*W3SUTI)
        B2UTSR=SQMQ/SQMH**2*(UH*(SH-UH)/(UH+SH)+2.*TH*SH*(SH+2.*UH)/
     &  (UH+SH)**2*(W1TR-W1HR)+(SQMQ-UH/4.)*(0.5*W2UR+0.5*W2HR-W2TR+
     &  W3UTSR)+UH2*(2.*SQMQ/(UH+SH)**2-0.5/(UH+SH))*(W2TR-W2HR)+
     &  0.5*TH*SH/UH*(W2HR-2.*W2TR)+0.125*(UH-12.*SQMQ-4.*TH*SH/UH)*
     &  W3TUSR)
        B2UTSI=SQMQ/SQMH**2*(2.*TH*SH*(SH+2.*UH)/(UH+SH)**2*
     &  (W1TI-W1HI)+(SQMQ-UH/4.)*(0.5*W2UI+0.5*W2HI-W2TI+W3UTSI)+
     &  UH2*(2.*SQMQ/(UH+SH)**2-0.5/(UH+SH))*(W2TI-W2HI)+0.5*TH*SH/UH*
     &  (W2HI-2.*W2TI)+0.125*(UH-12.*SQMQ-4.*TH*SH/UH)*W3TUSI)
        B4STUR=0.25*EPSH*(-2./3.+0.25*(EPSH-1.)*(W2SR-W2HR+W3STUR))
        B4STUI=0.25*EPSH*0.25*(EPSH-1.)*(W2SI-W2HI+W3STUI)
        B4TUSR=0.25*EPSH*(-2./3.+0.25*(EPSH-1.)*(W2TR-W2HR+W3TUSR))
        B4TUSI=0.25*EPSH*0.25*(EPSH-1.)*(W2TI-W2HI+W3TUSI)
        B4USTR=0.25*EPSH*(-2./3.+0.25*(EPSH-1.)*(W2UR-W2HR+W3USTR))
        B4USTI=0.25*EPSH*0.25*(EPSH-1.)*(W2UI-W2HI+W3USTI)
        A2STUR=A2STUR+B2STUR+B2SUTR
        A2STUI=A2STUI+B2STUI+B2SUTI
        A2USTR=A2USTR+B2USTR+B2UTSR
        A2USTI=A2USTI+B2USTI+B2UTSI
        A2TUSR=A2TUSR+B2TUSR+B2TSUR
        A2TUSI=A2TUSI+B2TUSI+B2TSUI
        A4STUR=A4STUR+B4STUR+B4USTR+B4TUSR
        A4STUI=A4STUI+B4STUI+B4USTI+B4TUSI
 1020   CONTINUE
        FACGH=COMFAC*FACA*3./(128.*PARU(1)**2)*AEM/XW*AS**3*
     &  SQMH/SQMW*SQMH**3/(SH*TH*UH)*(A2STUR**2+A2STUI**2+A2USTR**2+
     &  A2USTI**2+A2TUSR**2+A2TUSI**2+A4STUR**2+A4STUI**2)
        FACGH=FACGH*WIDS(25,2)
        IF(KFAC(1,21)*KFAC(2,21).EQ.0) GOTO 1030
        NCHN=NCHN+1
        ISIG(NCHN,1)=21
        ISIG(NCHN,2)=21
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACGH
 1030   CONTINUE

      ELSEIF(ISUB.EQ.114.OR.ISUB.EQ.115) THEN
C...g + g -> gamma + gamma or g + g -> g + gamma.
        A0STUR=0.
        A0STUI=0.
        A0TSUR=0.
        A0TSUI=0.
        A0UTSR=0.
        A0UTSI=0.
        A1STUR=0.
        A1STUI=0.
        A2STUR=0.
        A2STUI=0.
        ALST=LOG(-SH/TH)
        ALSU=LOG(-SH/UH)
        ALTU=LOG(TH/UH)
        IMAX=2*MSTP(1)
        IF(MSTP(38).GE.1.AND.MSTP(38).LE.8) IMAX=MSTP(38)
        DO 1040 I=1,IMAX
        EI=KCHG(IABS(I),1)/3.
        EIWT=EI**2
        IF(ISUB.EQ.115) EIWT=EI
        SQMQ=PMAS(I,1)**2
        EPSS=4.*SQMQ/SH
        EPST=4.*SQMQ/TH
        EPSU=4.*SQMQ/UH
        IF((MSTP(38).GE.1.AND.MSTP(38).LE.8).OR.EPSS.LT.1.E-4) THEN
          B0STUR=1.+(TH-UH)/SH*ALTU+0.5*(TH2+UH2)/SH2*(ALTU**2+
     &    PARU(1)**2)
          B0STUI=0.
          B0TSUR=1.+(SH-UH)/TH*ALSU+0.5*(SH2+UH2)/TH2*ALSU**2
          B0TSUI=-PARU(1)*((SH-UH)/TH+(SH2+UH2)/TH2*ALSU)
          B0UTSR=1.+(SH-TH)/UH*ALST+0.5*(SH2+TH2)/UH2*ALST**2
          B0UTSI=-PARU(1)*((SH-TH)/UH+(SH2+TH2)/UH2*ALST)
          B1STUR=-1.
          B1STUI=0.
          B2STUR=-1.
          B2STUI=0.
        ELSE
          CALL RYWAUX(1,EPSS,W1SR,W1SI)
          CALL RYWAUX(1,EPST,W1TR,W1TI)
          CALL RYWAUX(1,EPSU,W1UR,W1UI)
          CALL RYWAUX(2,EPSS,W2SR,W2SI)
          CALL RYWAUX(2,EPST,W2TR,W2TI)
          CALL RYWAUX(2,EPSU,W2UR,W2UI)
          CALL RYI3AU(EPSS,TH/UH,Y3STUR,Y3STUI)
          CALL RYI3AU(EPSS,UH/TH,Y3SUTR,Y3SUTI)
          CALL RYI3AU(EPST,SH/UH,Y3TSUR,Y3TSUI)
          CALL RYI3AU(EPST,UH/SH,Y3TUSR,Y3TUSI)
          CALL RYI3AU(EPSU,SH/TH,Y3USTR,Y3USTI)
          CALL RYI3AU(EPSU,TH/SH,Y3UTSR,Y3UTSI)
          B0STUR=1.+(1.+2.*TH/SH)*W1TR+(1.+2.*UH/SH)*W1UR+
     &    0.5*((TH2+UH2)/SH2-EPSS)*(W2TR+W2UR)-
     &    0.25*EPST*(1.-0.5*EPSS)*(Y3SUTR+Y3TUSR)-
     &    0.25*EPSU*(1.-0.5*EPSS)*(Y3STUR+Y3UTSR)+
     &    0.25*(-2.*(TH2+UH2)/SH2+4.*EPSS+EPST+EPSU+0.5*EPST*EPSU)*
     &    (Y3TSUR+Y3USTR)
          B0STUI=(1.+2.*TH/SH)*W1TI+(1.+2.*UH/SH)*W1UI+
     &    0.5*((TH2+UH2)/SH2-EPSS)*(W2TI+W2UI)-
     &    0.25*EPST*(1.-0.5*EPSS)*(Y3SUTI+Y3TUSI)-
     &    0.25*EPSU*(1.-0.5*EPSS)*(Y3STUI+Y3UTSI)+
     &    0.25*(-2.*(TH2+UH2)/SH2+4.*EPSS+EPST+EPSU+0.5*EPST*EPSU)*
     &    (Y3TSUI+Y3USTI)
          B0TSUR=1.+(1.+2.*SH/TH)*W1SR+(1.+2.*UH/TH)*W1UR+
     &    0.5*((SH2+UH2)/TH2-EPST)*(W2SR+W2UR)-
     &    0.25*EPSS*(1.-0.5*EPST)*(Y3TUSR+Y3SUTR)-
     &    0.25*EPSU*(1.-0.5*EPST)*(Y3TSUR+Y3USTR)+
     &    0.25*(-2.*(SH2+UH2)/TH2+4.*EPST+EPSS+EPSU+0.5*EPSS*EPSU)*
     &    (Y3STUR+Y3UTSR)
          B0TSUI=(1.+2.*SH/TH)*W1SI+(1.+2.*UH/TH)*W1UI+
     &    0.5*((SH2+UH2)/TH2-EPST)*(W2SI+W2UI)-
     &    0.25*EPSS*(1.-0.5*EPST)*(Y3TUSI+Y3SUTI)-
     &    0.25*EPSU*(1.-0.5*EPST)*(Y3TSUI+Y3USTI)+
     &    0.25*(-2.*(SH2+UH2)/TH2+4.*EPST+EPSS+EPSU+0.5*EPSS*EPSU)*
     &    (Y3STUI+Y3UTSI)
          B0UTSR=1.+(1.+2.*TH/UH)*W1TR+(1.+2.*SH/UH)*W1SR+
     &    0.5*((TH2+SH2)/UH2-EPSU)*(W2TR+W2SR)-
     &    0.25*EPST*(1.-0.5*EPSU)*(Y3USTR+Y3TSUR)-
     &    0.25*EPSS*(1.-0.5*EPSU)*(Y3UTSR+Y3STUR)+
     &    0.25*(-2.*(TH2+SH2)/UH2+4.*EPSU+EPST+EPSS+0.5*EPST*EPSS)*
     &    (Y3TUSR+Y3SUTR)
          B0UTSI=(1.+2.*TH/UH)*W1TI+(1.+2.*SH/UH)*W1SI+
     &    0.5*((TH2+SH2)/UH2-EPSU)*(W2TI+W2SI)-
     &    0.25*EPST*(1.-0.5*EPSU)*(Y3USTI+Y3TSUI)-
     &    0.25*EPSS*(1.-0.5*EPSU)*(Y3UTSI+Y3STUI)+
     &    0.25*(-2.*(TH2+SH2)/UH2+4.*EPSU+EPST+EPSS+0.5*EPST*EPSS)*
     &    (Y3TUSI+Y3SUTI)
          B1STUR=-1.-0.25*(EPSS+EPST+EPSU)*(W2SR+W2TR+W2UR)+
     &    0.25*(EPSU+0.5*EPSS*EPST)*(Y3SUTR+Y3TUSR)+
     &    0.25*(EPST+0.5*EPSS*EPSU)*(Y3STUR+Y3UTSR)+
     &    0.25*(EPSS+0.5*EPST*EPSU)*(Y3TSUR+Y3USTR)
          B1STUI=-0.25*(EPSS+EPST+EPSU)*(W2SI+W2TI+W2UI)+
     &    0.25*(EPSU+0.5*EPSS*EPST)*(Y3SUTI+Y3TUSI)+
     &    0.25*(EPST+0.5*EPSS*EPSU)*(Y3STUI+Y3UTSI)+
     &    0.25*(EPSS+0.5*EPST*EPSU)*(Y3TSUI+Y3USTI)
          B2STUR=-1.+0.125*EPSS*EPST*(Y3SUTR+Y3TUSR)+
     &    0.125*EPSS*EPSU*(Y3STUR+Y3UTSR)+
     &    0.125*EPST*EPSU*(Y3TSUR+Y3USTR)
          B2STUI=0.125*EPSS*EPST*(Y3SUTI+Y3TUSI)+
     &    0.125*EPSS*EPSU*(Y3STUI+Y3UTSI)+
     &    0.125*EPST*EPSU*(Y3TSUI+Y3USTI)
        ENDIF
        A0STUR=A0STUR+EIWT*B0STUR
        A0STUI=A0STUI+EIWT*B0STUI
        A0TSUR=A0TSUR+EIWT*B0TSUR
        A0TSUI=A0TSUI+EIWT*B0TSUI
        A0UTSR=A0UTSR+EIWT*B0UTSR
        A0UTSI=A0UTSI+EIWT*B0UTSI
        A1STUR=A1STUR+EIWT*B1STUR
        A1STUI=A1STUI+EIWT*B1STUI
        A2STUR=A2STUR+EIWT*B2STUR
        A2STUI=A2STUI+EIWT*B2STUI
 1040   CONTINUE
        ASQSUM=A0STUR**2+A0STUI**2+A0TSUR**2+A0TSUI**2+A0UTSR**2+
     &  A0UTSI**2+4.*A1STUR**2+4.*A1STUI**2+A2STUR**2+A2STUI**2
        FACGG=COMFAC*FACA/(16.*PARU(1)**2)*AS**2*AEM**2*ASQSUM
        FACGP=COMFAC*FACA*5./(192.*PARU(1)**2)*AS**3*AEM*ASQSUM
        IF(KFAC(1,21)*KFAC(2,21).EQ.0) GOTO 1050
        NCHN=NCHN+1
        ISIG(NCHN,1)=21
        ISIG(NCHN,2)=21
        ISIG(NCHN,3)=1
        IF(ISUB.EQ.114) SIGH(NCHN)=0.5*FACGG
        IF(ISUB.EQ.115) SIGH(NCHN)=FACGP
 1050   CONTINUE

      ELSEIF(ISUB.EQ.116) THEN
C...g + g -> gamma + Z0.

      ELSEIF(ISUB.EQ.117) THEN
C...g + g -> Z0 + Z0.

      ELSEIF(ISUB.EQ.118) THEN
C...g + g -> W+ + W-.

      ENDIF

C...G: 2 -> 3, tree diagrams.

      ELSEIF(ISUB.LE.140) THEN
      IF(ISUB.EQ.121) THEN
C...g + g -> f + f~ + H0 (f + f~ -> H0 as inner process).

      ELSEIF(ISUB.EQ.122) THEN
C...gamma + gamma -> f + f' + H0 (f + f~ -> H0 as inner process).

      ELSEIF(ISUB.EQ.123) THEN
C...f + f' -> f + f' + H0 (or H'0, or A0) (Z0 + Z0 -> H0 as
C...inner process).
        FACNOR=COMFAC*(4.*PARU(1)*AEM/(XW*(1.-XW)))**3*SQMZ/32.
        IF(MSTP(4).GE.1.OR.IHIGG.GE.2) FACNOR=FACNOR*
     &  PARU(154+10*IHIGG)**2
        FACPRP=1./((VINT(215)-VINT(204)**2)*(VINT(216)-VINT(209)**2))**2
        FACZZ1=FACNOR*FACPRP*(0.5*TAUP*VINT(2))*VINT(219)
        FACZZ2=FACNOR*FACPRP*VINT(217)*VINT(218)
        CALL RYWIDT(KFHIGG,SH,WDTP,WDTE)
        HP=AEM/(8.*XW)*SH/SQMW*SH
        HS=HP*WDTP(0)
        HF=HP*(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))
        FACBW=(1./PARU(1))*VINT(2)*HF/((SH-SQMH)**2+HS**2)
        IF(ABS(SH-SQMH).GT.100.*HS) FACBW=0.
        DO 1070 I=MIN1,MAX1
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 1070
        IA=IABS(I)
        DO 1060 J=MIN2,MAX2
        IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 1060
        JA=IABS(J)
        EI=KCHG(IA,1)*ISIGN(1,I)/3.
        AI=SIGN(1.,KCHG(IA,1)+0.5)*ISIGN(1,I)
        VI=AI-4.*EI*XW
        EJ=KCHG(JA,1)*ISIGN(1,J)/3.
        AJ=SIGN(1.,KCHG(JA,1)+0.5)*ISIGN(1,J)
        VJ=AJ-4.*EJ*XW
        FACLR1=(VI**2+AI**2)*(VJ**2+AJ**2)+4.*VI*AI*VJ*AJ
        FACLR2=(VI**2+AI**2)*(VJ**2+AJ**2)-4.*VI*AI*VJ*AJ
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=J
        ISIG(NCHN,3)=1
        SIGH(NCHN)=(FACLR1*FACZZ1+FACLR2*FACZZ2)*FACBW
 1060   CONTINUE
 1070   CONTINUE

      ELSEIF(ISUB.EQ.124) THEN
C...f + f' -> f" + f"' + H0 (or H'0, or A0) (W+ + W- -> H0 as
C...inner process).
        FACNOR=COMFAC*(4.*PARU(1)*AEM/XW)**3*SQMW
        IF(MSTP(4).GE.1.OR.IHIGG.GE.2) FACNOR=FACNOR*
     &  PARU(155+10*IHIGG)**2
        FACPRP=1./((VINT(215)-VINT(204)**2)*(VINT(216)-VINT(209)**2))**2
        FACWW=FACNOR*FACPRP*(0.5*TAUP*VINT(2))*VINT(219)
        CALL RYWIDT(KFHIGG,SH,WDTP,WDTE)
        HP=AEM/(8.*XW)*SH/SQMW*SH
        HS=HP*WDTP(0)
        HF=HP*(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))
        FACBW=(1./PARU(1))*VINT(2)*HF/((SH-SQMH)**2+HS**2)
        IF(ABS(SH-SQMH).GT.100.*HS) FACBW=0.
        DO 1090 I=MIN1,MAX1
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 1090
        EI=SIGN(1.,FLOAT(I))*KCHG(IABS(I),1)
        DO 1080 J=MIN2,MAX2
        IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 1080
        EJ=SIGN(1.,FLOAT(J))*KCHG(IABS(J),1)
        IF(EI*EJ.GT.0.) GOTO 1080
        FACLR=VINT(180+I)*VINT(180+J)
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=J
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACLR*FACWW*FACBW
 1080   CONTINUE
 1090   CONTINUE

      ELSEIF(ISUB.EQ.131) THEN
C...g + g -> Z0 + q + qbar.
        IF(KFAC(1,21)*KFAC(2,21).EQ.0) GOTO 1120

C...Read out information on flavours, masses, couplings.
        KFQ=KFPR(131,2)
        KFL=IABS(KFDP(MINT(35),1))
        PMH=SQRT(SH)
        PMQQ=SQRT(VINT(64))
        PMLL=SQRT(VINT(63))
        PMQ=PMAS(KFQ,1)
        QFQ=KCHG(KFQ,1)/3.
        AFQ=SIGN(1.,QFQ+0.1)
        VFQ=AFQ-4.*XW*QFQ
        QFL=KCHG(KFL,1)/3.
        AFL=SIGN(1.,QFL+0.1)
        VFL=AFL-4.*XW*QFL

C...Set line numbers for particles.
        IG1=MINT(84)+1
        IG2=MINT(84)+2
        IQ1=MINT(84)+3
        IQ2=MINT(84)+4
        IL1=MINT(84)+5
        IL2=MINT(84)+6
        IZ=MINT(84)+7

C...Reconstruct decay kinematics.
        DO 1100 I=MINT(84)+1,MINT(84)+7
        K(I,1)=1
        DO 1100 J=1,5
 1100   P(I,J)=0.
        P(IG1,4)=0.5*PMH
        P(IG1,3)=P(IG1,4)
        P(IG2,4)=P(IG1,4)
        P(IG2,3)=-P(IG1,3)
        P(IQ1,5)=PMQ
        P(IQ1,4)=0.5*PMQQ
        P(IQ1,3)=SQRT(MAX(0.,P(IQ1,4)**2-PMQ**2))
        P(IQ2,5)=PMQ
        P(IQ2,4)=P(IQ1,4)
        P(IQ2,3)=-P(IQ1,3)
        P(IL1,4)=0.5*PMLL
        P(IL1,3)=P(IL1,4)
        P(IL2,4)=P(IL1,4)
        P(IL2,3)=-P(IL1,3)
        P(IZ,5)=PMLL
        P(IZ,4)=0.5*(PMH+(PMLL**2-PMQQ**2)/PMH)
        P(IZ,3)=SQRT(MAX(0.,P(IZ,4)**2-PMLL**2))
        CALL LUDBRB(IQ1,IQ2,ACOS(VINT(83)),VINT(84),0D0,0D0,
     &  -DBLE(P(IZ,3)/(PMH-P(IZ,4))))
        CALL LUDBRB(IL1,IL2,ACOS(VINT(81)),VINT(82),0D0,0D0,
     &  DBLE(P(IZ,3)/P(IZ,4)))
        CALL LUDBRB(IQ1,IZ,ACOS(VINT(23)),VINT(24),0D0,0D0,0D0)

C...Interface information to program of Ronald Kleiss.
        RKMQ=PMQ
        RKMZ=PMAS(23,1)
        RKGZ=PMAS(23,2)
        RKVQ=VFQ
        RKAQ=AFQ
        RKVL=VFL
        RKAL=AFL
        RKG1(0)=P(IG1,4)
        RKG2(0)=P(IG2,4)
        RKQ1(0)=P(IQ1,4)
        RKQ2(0)=P(IQ2,4)
        RKL1(0)=P(IL1,4)
        RKL2(0)=P(IL2,4)
        DO 1110 J=1,3
        RKG1(J)=P(IG1,J)
        RKG2(J)=P(IG2,J)
        RKQ1(J)=P(IQ1,J)
        RKQ2(J)=P(IQ2,J)
        RKL1(J)=P(IL1,J)
        RKL2(J)=P(IL2,J)
 1110   CONTINUE
        CALL RKBBV(RKG1,RKG2,RKQ1,RKQ2,RKL1,RKL2,1,RKRES)

C...Multiply with normalization factors.
        WTMEP=1./(2.*SH*PARU(2)**8)
        WTCOU=AS**2*(4.*PARU(1)*AEM*XWC)**2
        WTZQQ=WTMEP*WTCOU*RKRES
        WTPHS=(PARU(1)/2.)**2*PMQQ**2*
     &  (PARU(1)*((PMLL**2-PMAS(23,1)**2)**2+(PMAS(23,1)*
     &  PMAS(23,2))**2)/(PMAS(23,1)*PMAS(23,2)))*0.5*SH
        NCHN=NCHN+1
        ISIG(NCHN,1)=21
        ISIG(NCHN,2)=21
        ISIG(NCHN,3)=INT(1.5+PYR(0))
        SIGH(NCHN)=COMFAC*WTPHS*WTZQQ
 1120   CONTINUE
      ENDIF

C...H: 2 -> 1, tree diagrams, non-standard model processes.

      ELSEIF(ISUB.LE.160) THEN
      IF(ISUB.EQ.141) THEN
C...f + f~ -> gamma*/Z0/Z'0.
        MINT(61)=2
        CALL RYWIDT(32,SH,WDTP,WDTE)
        HP0=AEM/3.*SH
        HP1=AEM/3.*XWC*SH
        HP2=HP1
        HS=HP1*VINT(117)
        HSP=HP2*WDTP(0)
        FACZP=4.*COMFAC*3.
        DO 1130 I=MINA,MAXA
        IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 1130
        EI=KCHG(IABS(I),1)/3.
        AI=SIGN(1.,EI)
        VI=AI-4.*EI*XW
        IF(IABS(I).LT.10) THEN
          VPI=PARU(123-2*MOD(IABS(I),2))
          API=PARU(124-2*MOD(IABS(I),2))
        ELSE
          VPI=PARU(127-2*MOD(IABS(I),2))
          API=PARU(128-2*MOD(IABS(I),2))
        ENDIF
        HI0=HP0
        IF(IABS(I).LE.10) HI0=HI0*FACA/3.
        HI1=HP1
        IF(IABS(I).LE.10) HI1=HI1*FACA/3.
        HI2=HP2
        IF(IABS(I).LE.10) HI2=HI2*FACA/3.
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=-I
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACZP*(EI**2/SH2*HI0*HP0*VINT(111)+EI*VI*
     &  (1.-SQMZ/SH)/((SH-SQMZ)**2+HS**2)*(HI0*HP1+HI1*HP0)*VINT(112)+
     &  EI*VPI*(1.-SQMZP/SH)/((SH-SQMZP)**2+HSP**2)*(HI0*HP2+HI2*HP0)*
     &  VINT(113)+(VI**2+AI**2)/((SH-SQMZ)**2+HS**2)*HI1*HP1*VINT(114)+
     &  (VI*VPI+AI*API)*((SH-SQMZ)*(SH-SQMZP)+HS*HSP)/(((SH-SQMZ)**2+
     &  HS**2)*((SH-SQMZP)**2+HSP**2))*(HI1*HP2+HI2*HP1)*VINT(115)+
     &  (VPI**2+API**2)/((SH-SQMZP)**2+HSP**2)*HI2*HP2*VINT(116))
 1130   CONTINUE

      ELSEIF(ISUB.EQ.142) THEN
C...f + f~' -> W'+/-.
        CALL RYWIDT(34,SH,WDTP,WDTE)
        HP=AEM/(24.*XW)*SH
        HS=HP*WDTP(0)
        FACBW=4.*COMFAC/((SH-SQMWP)**2+HS**2)*3.
        DO 1150 I=MIN1,MAX1
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 1150
        IA=IABS(I)
        DO 1140 J=MIN2,MAX2
        IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 1140
        JA=IABS(J)
        IF(I*J.GT.0.OR.MOD(IA+JA,2).EQ.0) GOTO 1140
        IF((IA.LE.10.AND.JA.GT.10).OR.(IA.GT.10.AND.JA.LE.10)) GOTO 1140
        KCHW=(KCHG(IA,1)*ISIGN(1,I)+KCHG(JA,1)*ISIGN(1,J))/3
        HI=HP*(PARU(133)**2+PARU(134)**2)
        IF(IA.LE.10) HI=HP*(PARU(131)**2+PARU(132)**2)*
     &  VCKM((IA+1)/2,(JA+1)/2)*FACA/3.
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=J
        ISIG(NCHN,3)=1
        HF=HP*(WDTE(0,1)+WDTE(0,(5-KCHW)/2)+WDTE(0,4))
        SIGH(NCHN)=HI*FACBW*HF
 1140   CONTINUE
 1150   CONTINUE

      ELSEIF(ISUB.EQ.143) THEN
C...f + f~' -> H+/-.
        CALL RYWIDT(37,SH,WDTP,WDTE)
        HP=AEM/(8.*XW)*SH/SQMW*SH
        HS=HP*WDTP(0)
        FACBW=4.*COMFAC/((SH-SQMHC)**2+HS**2)
        DO 1170 I=MIN1,MAX1
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 1170
        IA=IABS(I)
        IM=(MOD(IA,10)+1)/2
        DO 1160 J=MIN2,MAX2
        IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 1160
        JA=IABS(J)
        JM=(MOD(JA,10)+1)/2
        IF(I*J.GT.0.OR.IA.EQ.JA.OR.IM.NE.JM) GOTO 1160
        IF((IA.LE.10.AND.JA.GT.10).OR.(IA.GT.10.AND.JA.LE.10)) GOTO 1160
        IF(MOD(IA,2).EQ.0) THEN
          IU=IA
          IL=JA
        ELSE
          IU=JA
          IL=IA
        ENDIF
        RML=PMAS(IL,1)**2/SH
        RMU=PMAS(IU,1)**2/SH
        HI=HP*(RML*PARU(141)**2+RMU/PARU(141)**2)
        IF(IA.LE.10) HI=HI*FACA/3.
        KCHHC=(KCHG(IA,1)*ISIGN(1,I)+KCHG(JA,1)*ISIGN(1,J))/3
        HF=HP*(WDTE(0,1)+WDTE(0,(5-KCHHC)/2)+WDTE(0,4))
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=J
        ISIG(NCHN,3)=1
        SIGH(NCHN)=HI*FACBW*HF
 1160   CONTINUE
 1170   CONTINUE

      ELSEIF(ISUB.EQ.144) THEN
C...f + f~' -> R.
        CALL RYWIDT(40,SH,WDTP,WDTE)
        HP=AEM/(12.*XW)*SH
        HS=HP*WDTP(0)
        FACBW=4.*COMFAC/((SH-SQMR)**2+HS**2)*3.
        DO 1190 I=MIN1,MAX1
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 1190
        IA=IABS(I)
        DO 1180 J=MIN2,MAX2
        IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 1180
        JA=IABS(J)
        IF(I*J.GT.0.OR.IABS(IA-JA).NE.2) GOTO 1180
        HI=HP
        IF(IA.LE.10) HI=HI*FACA/3.
        HF=HP*(WDTE(0,1)+WDTE(0,(10-(I+J))/4)+WDTE(0,4))
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=J
        ISIG(NCHN,3)=1
        SIGH(NCHN)=HI*FACBW*HF
 1180   CONTINUE
 1190   CONTINUE

      ELSEIF(ISUB.EQ.145) THEN
C...q + l -> LQ (leptoquark).
        CALL RYWIDT(39,SH,WDTP,WDTE)
        HP=AEM/4.*SH
        HS=HP*WDTP(0)
        FACBW=4.*COMFAC/((SH-SQMLQ)**2+HS**2)
        IF(ABS(SH-SQMLQ).GT.100.*HS) FACBW=0.
        KFLQQ=KFDP(MDCY(39,2),1)
        KFLQL=KFDP(MDCY(39,2),2)
        DO 1210 I=MIN1,MAX1
        IF(KFAC(1,I).EQ.0) GOTO 1210
        IA=IABS(I)
        IF(IA.NE.KFLQQ.AND.IA.NE.KFLQL) GOTO 1210
        DO 1200 J=MIN2,MAX2
        IF(KFAC(2,J).EQ.0) GOTO 1200
        JA=IABS(J)
        IF(JA.NE.KFLQQ.AND.JA.NE.KFLQL) GOTO 1200
        IF(I*J.NE.KFLQQ*KFLQL) GOTO 1200
        IF(IA.EQ.KFLQQ) KCHLQ=ISIGN(1,I)
        IF(JA.EQ.KFLQQ) KCHLQ=ISIGN(1,J)
        HI=HP*PARU(151)
        HF=HP*(WDTE(0,1)+WDTE(0,(5-KCHLQ)/2)+WDTE(0,4))
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=J
        ISIG(NCHN,3)=1
        SIGH(NCHN)=HI*FACBW*HF
 1200   CONTINUE
 1210   CONTINUE

      ENDIF

C...I: 2 -> 2, tree diagrams, non-standard model processes.

      ELSE
      IF(ISUB.EQ.161) THEN
C...f + g -> f' + H+/- (b + g -> t + H+/- only)
C...(choice of only b and t to avoid kinematics problems).
        FHCQ=COMFAC*FACA*AS*AEM/XW*1./24
        DO 1230 I=MINA,MAXA
        IA=IABS(I)
        IF(IA.NE.5) GOTO 1230
        IUA=IA+MOD(IA,2)
        SQMQ=PMAS(IUA,1)**2
        FACHCQ=FHCQ/PARU(141)**2*SQMQ/SQMW*(SH/(SQMQ-UH)+
     &  2.*SQMQ*(SQMHC-UH)/(SQMQ-UH)**2+(SQMQ-UH)/SH+
     &  2.*SQMQ/(SQMQ-UH)+2.*(SQMHC-UH)/(SQMQ-UH)*(SQMHC-SQMQ-SH)/SH)
        KCHHC=ISIGN(1,KCHG(IA,1)*ISIGN(1,I))
        DO 1220 ISDE=1,2
        IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,21).EQ.0) GOTO 1220
        IF(ISDE.EQ.2.AND.KFAC(1,21)*KFAC(2,1).EQ.0) GOTO 1220
        NCHN=NCHN+1
        ISIG(NCHN,ISDE)=I
        ISIG(NCHN,3-ISDE)=21
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACHCQ*WIDS(37,(5-KCHHC)/2)
 1220   CONTINUE
 1230   CONTINUE

      ELSEIF(ISUB.EQ.162) THEN
C...q + g -> LQ + l~; LQ=leptoquark.
        FACLQ=COMFAC*FACA*PARU(151)*(AS*AEM/6.)*(-TH/SH)*
     &  (UH2+SQMLQ**2)/(UH-SQMLQ)**2
        KFLQQ=KFDP(MDCY(39,2),1)
        DO 1250 I=MINA,MAXA
        IF(IABS(I).NE.KFLQQ) GOTO 1250
        KCHLQ=ISIGN(1,I)
        DO 1240 ISDE=1,2
        IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,21).EQ.0) GOTO 1240
        IF(ISDE.EQ.2.AND.KFAC(1,21)*KFAC(2,I).EQ.0) GOTO 1240
        NCHN=NCHN+1
        ISIG(NCHN,ISDE)=I
        ISIG(NCHN,3-ISDE)=21
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACLQ*WIDS(39,(5-KCHLQ)/2)
 1240   CONTINUE
 1250   CONTINUE

      ELSEIF(ISUB.EQ.163) THEN
C...g + g -> LQ + LQ~; LQ=leptoquark.
        FACLQ=COMFAC*FACA*WIDS(39,1)*(AS**2/2.)*
     &  (7./48.+3.*(UH-TH)**2/(16.*SH2))*(1.+2.*SQMLQ*TH/(TH-SQMLQ)**2+
     &  2.*SQMLQ*UH/(UH-SQMLQ)**2+4.*SQMLQ**2/((TH-SQMLQ)*(UH-SQMLQ)))
        IF(KFAC(1,21)*KFAC(2,21).EQ.0) GOTO 1260
        NCHN=NCHN+1
        ISIG(NCHN,1)=21
        ISIG(NCHN,2)=21
C...Since don't know proper colour flow, randomize between alternatives.
        ISIG(NCHN,3)=INT(1.5+PYR(0))
        SIGH(NCHN)=FACLQ
 1260   CONTINUE

      ELSEIF(ISUB.EQ.164) THEN
C...q + q~ -> LQ + LQ~; LQ=leptoquark.
        FACLQA=COMFAC*WIDS(39,1)*(AS**2/9.)*
     &  (SH*(SH-4.*SQMLQ)-(UH-TH)**2)/SH2
        FACLQS=COMFAC*WIDS(39,1)*((PARU(151)**2*AEM**2/8.)*
     &  (-SH*TH-(SQMLQ-TH)**2)/TH2+(PARU(151)*AEM*AS/18.)*
     &  ((SQMLQ-TH)*(UH-TH)+SH*(SQMLQ+TH))/(SH*TH))
        KFLQQ=KFDP(MDCY(39,2),1)
        DO 1270 I=MINA,MAXA
        IF(I.EQ.0.OR.IABS(I).GT.MSTP(54).OR.
     &  KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 1270
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=-I
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACLQA
        IF(IABS(I).EQ.KFLQQ) SIGH(NCHN)=FACLQA+FACLQS
 1270   CONTINUE

      ENDIF
      ENDIF

C...Multiply with structure functions.
      IF(ISUB.LE.90.OR.ISUB.GE.96) THEN
        DO 1280 ICHN=1,NCHN
        IF(MINT(45).GE.2) THEN
          KFL1=ISIG(ICHN,1)
          IF(KFL1.EQ.21) KFL1=0
          SIGH(ICHN)=SIGH(ICHN)*XSFX(1,KFL1)
        ENDIF
        IF(MINT(46).GE.2) THEN
          KFL2=ISIG(ICHN,2)
          IF(KFL2.EQ.21) KFL2=0
          SIGH(ICHN)=SIGH(ICHN)*XSFX(2,KFL2)
        ENDIF
 1280   SIGS=SIGS+SIGH(ICHN)
      ENDIF

      RETURN
      END

C*********************************************************************

      SUBROUTINE RYSTFU(KF,X,Q2,XPQ)

C...Gives electron, photon, pi+, neutron and proton parton structure
C...functions according to a few different parametrizations. Note
C...that what is coded is x times the probability distribution,
C...i.e. xq(x,Q2) etc.
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/RYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/RYINT1/MINT(400),VINT(400)
      SAVE /LUDAT1/,/LUDAT2/
      SAVE /RYPARS/,/RYINT1/
      DIMENSION XPQ(-25:25),XQ(9),TX(6),TT(6),TS(6),NEHLQ(8,2),
     &CEHLQ(6,6,2,8,2),CDO(3,6,5,2),COW(3,5,4,2),CMT(0:3,0:2,9,4),
     &EXMT(0:3)

C...The following data lines are coefficients needed in the
C...Owens pion structure function parametrizations, see below.
C...Expansion coefficients for up and down valence quark distributions.
      DATA ((COW(IP,IS,1,1),IS=1,5),IP=1,3)/
     1  4.0000E-01,  7.0000E-01,  0.0000E+00,  0.0000E+00,  0.0000E+00,
     2 -6.2120E-02,  6.4780E-01,  0.0000E+00,  0.0000E+00,  0.0000E+00,
     3 -7.1090E-03,  1.3350E-02,  0.0000E+00,  0.0000E+00,  0.0000E+00/
      DATA ((COW(IP,IS,1,2),IS=1,5),IP=1,3)/
     1  4.0000E-01,  6.2800E-01,  0.0000E+00,  0.0000E+00,  0.0000E+00,
     2 -5.9090E-02,  6.4360E-01,  0.0000E+00,  0.0000E+00,  0.0000E+00,
     3 -6.5240E-03,  1.4510E-02,  0.0000E+00,  0.0000E+00,  0.0000E+00/
C...Expansion coefficients for gluon distribution.
      DATA ((COW(IP,IS,2,1),IS=1,5),IP=1,3)/
     1  8.8800E-01,  0.0000E+00,  3.1100E+00,  6.0000E+00,  0.0000E+00,
     2 -1.8020E+00, -1.5760E+00, -1.3170E-01,  2.8010E+00, -1.7280E+01,
     3  1.8120E+00,  1.2000E+00,  5.0680E-01, -1.2160E+01,  2.0490E+01/
      DATA ((COW(IP,IS,2,2),IS=1,5),IP=1,3)/
     1  7.9400E-01,  0.0000E+00,  2.8900E+00,  6.0000E+00,  0.0000E+00,
     2 -9.1440E-01, -1.2370E+00,  5.9660E-01, -3.6710E+00, -8.1910E+00,
     3  5.9660E-01,  6.5820E-01, -2.5500E-01, -2.3040E+00,  7.7580E+00/
C...Expansion coefficients for (up+down+strange) quark sea distribution.
      DATA ((COW(IP,IS,3,1),IS=1,5),IP=1,3)/
     1  9.0000E-01,  0.0000E+00,  5.0000E+00,  0.0000E+00,  0.0000E+00,
     2 -2.4280E-01, -2.1200E-01,  8.6730E-01,  1.2660E+00,  2.3820E+00,
     3  1.3860E-01,  3.6710E-03,  4.7470E-02, -2.2150E+00,  3.4820E-01/
      DATA ((COW(IP,IS,3,2),IS=1,5),IP=1,3)/
     1  9.0000E-01,  0.0000E+00,  5.0000E+00,  0.0000E+00,  0.0000E+00,
     2 -1.4170E-01, -1.6970E-01, -2.4740E+00, -2.5340E+00,  5.6210E-01,
     3 -1.7400E-01, -9.6230E-02,  1.5750E+00,  1.3780E+00, -2.7010E-01/
C...Expansion coefficients for charm quark sea distribution.
      DATA ((COW(IP,IS,4,1),IS=1,5),IP=1,3)/
     1  0.0000E+00, -2.2120E-02,  2.8940E+00,  0.0000E+00,  0.0000E+00,
     2  7.9280E-02, -3.7850E-01,  9.4330E+00,  5.2480E+00,  8.3880E+00,
     3 -6.1340E-02, -1.0880E-01, -1.0852E+01, -7.1870E+00, -1.1610E+01/
      DATA ((COW(IP,IS,4,2),IS=1,5),IP=1,3)/
     1  0.0000E+00, -8.8200E-02,  1.9240E+00,  0.0000E+00,  0.0000E+00,
     2  6.2290E-02, -2.8920E-01,  2.4240E-01, -4.4630E+00, -8.3670E-01,
     3 -4.0990E-02, -1.0820E-01,  2.0360E+00,  5.2090E+00, -4.8400E-02/

C...The following data lines are coefficients needed in the
C...Eichten, Hinchliffe, Lane, Quigg proton structure function
C...parametrizations, see below.
C...Powers of 1-x in different cases.
      DATA NEHLQ/3,4,7,5,7,7,7,7,3,4,7,6,7,7,7,7/
C...Expansion coefficients for up valence quark distribution.
      DATA (((CEHLQ(IX,IT,NX,1,1),IX=1,6),IT=1,6),NX=1,2)/
     1 7.677E-01,-2.087E-01,-3.303E-01,-2.517E-02,-1.570E-02,-1.000E-04,
     2-5.326E-01,-2.661E-01, 3.201E-01, 1.192E-01, 2.434E-02, 7.620E-03,
     3 2.162E-01, 1.881E-01,-8.375E-02,-6.515E-02,-1.743E-02,-5.040E-03,
     4-9.211E-02,-9.952E-02, 1.373E-02, 2.506E-02, 8.770E-03, 2.550E-03,
     5 3.670E-02, 4.409E-02, 9.600E-04,-7.960E-03,-3.420E-03,-1.050E-03,
     6-1.549E-02,-2.026E-02,-3.060E-03, 2.220E-03, 1.240E-03, 4.100E-04,
     1 2.395E-01, 2.905E-01, 9.778E-02, 2.149E-02, 3.440E-03, 5.000E-04,
     2 1.751E-02,-6.090E-03,-2.687E-02,-1.916E-02,-7.970E-03,-2.750E-03,
     3-5.760E-03,-5.040E-03, 1.080E-03, 2.490E-03, 1.530E-03, 7.500E-04,
     4 1.740E-03, 1.960E-03, 3.000E-04,-3.400E-04,-2.900E-04,-1.800E-04,
     5-5.300E-04,-6.400E-04,-1.700E-04, 4.000E-05, 6.000E-05, 4.000E-05,
     6 1.700E-04, 2.200E-04, 8.000E-05, 1.000E-05,-1.000E-05,-1.000E-05/
      DATA (((CEHLQ(IX,IT,NX,1,2),IX=1,6),IT=1,6),NX=1,2)/
     1 7.237E-01,-2.189E-01,-2.995E-01,-1.909E-02,-1.477E-02, 2.500E-04,
     2-5.314E-01,-2.425E-01, 3.283E-01, 1.119E-01, 2.223E-02, 7.070E-03,
     3 2.289E-01, 1.890E-01,-9.859E-02,-6.900E-02,-1.747E-02,-5.080E-03,
     4-1.041E-01,-1.084E-01, 2.108E-02, 2.975E-02, 9.830E-03, 2.830E-03,
     5 4.394E-02, 5.116E-02,-1.410E-03,-1.055E-02,-4.230E-03,-1.270E-03,
     6-1.991E-02,-2.539E-02,-2.780E-03, 3.430E-03, 1.720E-03, 5.500E-04,
     1 2.410E-01, 2.884E-01, 9.369E-02, 1.900E-02, 2.530E-03, 2.400E-04,
     2 1.765E-02,-9.220E-03,-3.037E-02,-2.085E-02,-8.440E-03,-2.810E-03,
     3-6.450E-03,-5.260E-03, 1.720E-03, 3.110E-03, 1.830E-03, 8.700E-04,
     4 2.120E-03, 2.320E-03, 2.600E-04,-4.900E-04,-3.900E-04,-2.300E-04,
     5-6.900E-04,-8.200E-04,-2.000E-04, 7.000E-05, 9.000E-05, 6.000E-05,
     6 2.400E-04, 3.100E-04, 1.100E-04, 0.000E+00,-2.000E-05,-2.000E-05/
C...Expansion coefficients for down valence quark distribution.
      DATA (((CEHLQ(IX,IT,NX,2,1),IX=1,6),IT=1,6),NX=1,2)/
     1 3.813E-01,-8.090E-02,-1.634E-01,-2.185E-02,-8.430E-03,-6.200E-04,
     2-2.948E-01,-1.435E-01, 1.665E-01, 6.638E-02, 1.473E-02, 4.080E-03,
     3 1.252E-01, 1.042E-01,-4.722E-02,-3.683E-02,-1.038E-02,-2.860E-03,
     4-5.478E-02,-5.678E-02, 8.900E-03, 1.484E-02, 5.340E-03, 1.520E-03,
     5 2.220E-02, 2.567E-02,-3.000E-05,-4.970E-03,-2.160E-03,-6.500E-04,
     6-9.530E-03,-1.204E-02,-1.510E-03, 1.510E-03, 8.300E-04, 2.700E-04,
     1 1.261E-01, 1.354E-01, 3.958E-02, 8.240E-03, 1.660E-03, 4.500E-04,
     2 3.890E-03,-1.159E-02,-1.625E-02,-9.610E-03,-3.710E-03,-1.260E-03,
     3-1.910E-03,-5.600E-04, 1.590E-03, 1.590E-03, 8.400E-04, 3.900E-04,
     4 6.400E-04, 4.900E-04,-1.500E-04,-2.900E-04,-1.800E-04,-1.000E-04,
     5-2.000E-04,-1.900E-04, 0.000E+00, 6.000E-05, 4.000E-05, 3.000E-05,
     6 7.000E-05, 8.000E-05, 2.000E-05,-1.000E-05,-1.000E-05,-1.000E-05/
      DATA (((CEHLQ(IX,IT,NX,2,2),IX=1,6),IT=1,6),NX=1,2)/
     1 3.578E-01,-8.622E-02,-1.480E-01,-1.840E-02,-7.820E-03,-4.500E-04,
     2-2.925E-01,-1.304E-01, 1.696E-01, 6.243E-02, 1.353E-02, 3.750E-03,
     3 1.318E-01, 1.041E-01,-5.486E-02,-3.872E-02,-1.038E-02,-2.850E-03,
     4-6.162E-02,-6.143E-02, 1.303E-02, 1.740E-02, 5.940E-03, 1.670E-03,
     5 2.643E-02, 2.957E-02,-1.490E-03,-6.450E-03,-2.630E-03,-7.700E-04,
     6-1.218E-02,-1.497E-02,-1.260E-03, 2.240E-03, 1.120E-03, 3.500E-04,
     1 1.263E-01, 1.334E-01, 3.732E-02, 7.070E-03, 1.260E-03, 3.400E-04,
     2 3.660E-03,-1.357E-02,-1.795E-02,-1.031E-02,-3.880E-03,-1.280E-03,
     3-2.100E-03,-3.600E-04, 2.050E-03, 1.920E-03, 9.800E-04, 4.400E-04,
     4 7.700E-04, 5.400E-04,-2.400E-04,-3.900E-04,-2.400E-04,-1.300E-04,
     5-2.600E-04,-2.300E-04, 2.000E-05, 9.000E-05, 6.000E-05, 4.000E-05,
     6 9.000E-05, 1.000E-04, 2.000E-05,-2.000E-05,-2.000E-05,-1.000E-05/
C...Expansion coefficients for up and down sea quark distributions.
      DATA (((CEHLQ(IX,IT,NX,3,1),IX=1,6),IT=1,6),NX=1,2)/
     1 6.870E-02,-6.861E-02, 2.973E-02,-5.400E-03, 3.780E-03,-9.700E-04,
     2-1.802E-02, 1.400E-04, 6.490E-03,-8.540E-03, 1.220E-03,-1.750E-03,
     3-4.650E-03, 1.480E-03,-5.930E-03, 6.000E-04,-1.030E-03,-8.000E-05,
     4 6.440E-03, 2.570E-03, 2.830E-03, 1.150E-03, 7.100E-04, 3.300E-04,
     5-3.930E-03,-2.540E-03,-1.160E-03,-7.700E-04,-3.600E-04,-1.900E-04,
     6 2.340E-03, 1.930E-03, 5.300E-04, 3.700E-04, 1.600E-04, 9.000E-05,
     1 1.014E+00,-1.106E+00, 3.374E-01,-7.444E-02, 8.850E-03,-8.700E-04,
     2 9.233E-01,-1.285E+00, 4.475E-01,-9.786E-02, 1.419E-02,-1.120E-03,
     3 4.888E-02,-1.271E-01, 8.606E-02,-2.608E-02, 4.780E-03,-6.000E-04,
     4-2.691E-02, 4.887E-02,-1.771E-02, 1.620E-03, 2.500E-04,-6.000E-05,
     5 7.040E-03,-1.113E-02, 1.590E-03, 7.000E-04,-2.000E-04, 0.000E+00,
     6-1.710E-03, 2.290E-03, 3.800E-04,-3.500E-04, 4.000E-05, 1.000E-05/
      DATA (((CEHLQ(IX,IT,NX,3,2),IX=1,6),IT=1,6),NX=1,2)/
     1 1.008E-01,-7.100E-02, 1.973E-02,-5.710E-03, 2.930E-03,-9.900E-04,
     2-5.271E-02,-1.823E-02, 1.792E-02,-6.580E-03, 1.750E-03,-1.550E-03,
     3 1.220E-02, 1.763E-02,-8.690E-03,-8.800E-04,-1.160E-03,-2.100E-04,
     4-1.190E-03,-7.180E-03, 2.360E-03, 1.890E-03, 7.700E-04, 4.100E-04,
     5-9.100E-04, 2.040E-03,-3.100E-04,-1.050E-03,-4.000E-04,-2.400E-04,
     6 1.190E-03,-1.700E-04,-2.000E-04, 4.200E-04, 1.700E-04, 1.000E-04,
     1 1.081E+00,-1.189E+00, 3.868E-01,-8.617E-02, 1.115E-02,-1.180E-03,
     2 9.917E-01,-1.396E+00, 4.998E-01,-1.159E-01, 1.674E-02,-1.720E-03,
     3 5.099E-02,-1.338E-01, 9.173E-02,-2.885E-02, 5.890E-03,-6.500E-04,
     4-3.178E-02, 5.703E-02,-2.070E-02, 2.440E-03, 1.100E-04,-9.000E-05,
     5 8.970E-03,-1.392E-02, 2.050E-03, 6.500E-04,-2.300E-04, 2.000E-05,
     6-2.340E-03, 3.010E-03, 5.000E-04,-3.900E-04, 6.000E-05, 1.000E-05/
C...Expansion coefficients for gluon distribution.
      DATA (((CEHLQ(IX,IT,NX,4,1),IX=1,6),IT=1,6),NX=1,2)/
     1 9.482E-01,-9.578E-01, 1.009E-01,-1.051E-01, 3.456E-02,-3.054E-02,
     2-9.627E-01, 5.379E-01, 3.368E-01,-9.525E-02, 1.488E-02,-2.051E-02,
     3 4.300E-01,-8.306E-02,-3.372E-01, 4.902E-02,-9.160E-03, 1.041E-02,
     4-1.925E-01,-1.790E-02, 2.183E-01, 7.490E-03, 4.140E-03,-1.860E-03,
     5 8.183E-02, 1.926E-02,-1.072E-01,-1.944E-02,-2.770E-03,-5.200E-04,
     6-3.884E-02,-1.234E-02, 5.410E-02, 1.879E-02, 3.350E-03, 1.040E-03,
     1 2.948E+01,-3.902E+01, 1.464E+01,-3.335E+00, 5.054E-01,-5.915E-02,
     2 2.559E+01,-3.955E+01, 1.661E+01,-4.299E+00, 6.904E-01,-8.243E-02,
     3-1.663E+00, 1.176E+00, 1.118E+00,-7.099E-01, 1.948E-01,-2.404E-02,
     4-2.168E-01, 8.170E-01,-7.169E-01, 1.851E-01,-1.924E-02,-3.250E-03,
     5 2.088E-01,-4.355E-01, 2.239E-01,-2.446E-02,-3.620E-03, 1.910E-03,
     6-9.097E-02, 1.601E-01,-5.681E-02,-2.500E-03, 2.580E-03,-4.700E-04/
      DATA (((CEHLQ(IX,IT,NX,4,2),IX=1,6),IT=1,6),NX=1,2)/
     1 2.367E+00, 4.453E-01, 3.660E-01, 9.467E-02, 1.341E-01, 1.661E-02,
     2-3.170E+00,-1.795E+00, 3.313E-02,-2.874E-01,-9.827E-02,-7.119E-02,
     3 1.823E+00, 1.457E+00,-2.465E-01, 3.739E-02, 6.090E-03, 1.814E-02,
     4-1.033E+00,-9.827E-01, 2.136E-01, 1.169E-01, 5.001E-02, 1.684E-02,
     5 5.133E-01, 5.259E-01,-1.173E-01,-1.139E-01,-4.988E-02,-2.021E-02,
     6-2.881E-01,-3.145E-01, 5.667E-02, 9.161E-02, 4.568E-02, 1.951E-02,
     1 3.036E+01,-4.062E+01, 1.578E+01,-3.699E+00, 6.020E-01,-7.031E-02,
     2 2.700E+01,-4.167E+01, 1.770E+01,-4.804E+00, 7.862E-01,-1.060E-01,
     3-1.909E+00, 1.357E+00, 1.127E+00,-7.181E-01, 2.232E-01,-2.481E-02,
     4-2.488E-01, 9.781E-01,-8.127E-01, 2.094E-01,-2.997E-02,-4.710E-03,
     5 2.506E-01,-5.427E-01, 2.672E-01,-3.103E-02,-1.800E-03, 2.870E-03,
     6-1.128E-01, 2.087E-01,-6.972E-02,-2.480E-03, 2.630E-03,-8.400E-04/
C...Expansion coefficients for strange sea quark distribution.
      DATA (((CEHLQ(IX,IT,NX,5,1),IX=1,6),IT=1,6),NX=1,2)/
     1 4.968E-02,-4.173E-02, 2.102E-02,-3.270E-03, 3.240E-03,-6.700E-04,
     2-6.150E-03,-1.294E-02, 6.740E-03,-6.890E-03, 9.000E-04,-1.510E-03,
     3-8.580E-03, 5.050E-03,-4.900E-03,-1.600E-04,-9.400E-04,-1.500E-04,
     4 7.840E-03, 1.510E-03, 2.220E-03, 1.400E-03, 7.000E-04, 3.500E-04,
     5-4.410E-03,-2.220E-03,-8.900E-04,-8.500E-04,-3.600E-04,-2.000E-04,
     6 2.520E-03, 1.840E-03, 4.100E-04, 3.900E-04, 1.600E-04, 9.000E-05,
     1 9.235E-01,-1.085E+00, 3.464E-01,-7.210E-02, 9.140E-03,-9.100E-04,
     2 9.315E-01,-1.274E+00, 4.512E-01,-9.775E-02, 1.380E-02,-1.310E-03,
     3 4.739E-02,-1.296E-01, 8.482E-02,-2.642E-02, 4.760E-03,-5.700E-04,
     4-2.653E-02, 4.953E-02,-1.735E-02, 1.750E-03, 2.800E-04,-6.000E-05,
     5 6.940E-03,-1.132E-02, 1.480E-03, 6.500E-04,-2.100E-04, 0.000E+00,
     6-1.680E-03, 2.340E-03, 4.200E-04,-3.400E-04, 5.000E-05, 1.000E-05/
      DATA (((CEHLQ(IX,IT,NX,5,2),IX=1,6),IT=1,6),NX=1,2)/
     1 6.478E-02,-4.537E-02, 1.643E-02,-3.490E-03, 2.710E-03,-6.700E-04,
     2-2.223E-02,-2.126E-02, 1.247E-02,-6.290E-03, 1.120E-03,-1.440E-03,
     3-1.340E-03, 1.362E-02,-6.130E-03,-7.900E-04,-9.000E-04,-2.000E-04,
     4 5.080E-03,-3.610E-03, 1.700E-03, 1.830E-03, 6.800E-04, 4.000E-04,
     5-3.580E-03, 6.000E-05,-2.600E-04,-1.050E-03,-3.800E-04,-2.300E-04,
     6 2.420E-03, 9.300E-04,-1.000E-04, 4.500E-04, 1.700E-04, 1.100E-04,
     1 9.868E-01,-1.171E+00, 3.940E-01,-8.459E-02, 1.124E-02,-1.250E-03,
     2 1.001E+00,-1.383E+00, 5.044E-01,-1.152E-01, 1.658E-02,-1.830E-03,
     3 4.928E-02,-1.368E-01, 9.021E-02,-2.935E-02, 5.800E-03,-6.600E-04,
     4-3.133E-02, 5.785E-02,-2.023E-02, 2.630E-03, 1.600E-04,-8.000E-05,
     5 8.840E-03,-1.416E-02, 1.900E-03, 5.800E-04,-2.500E-04, 1.000E-05,
     6-2.300E-03, 3.080E-03, 5.500E-04,-3.700E-04, 7.000E-05, 1.000E-05/
C...Expansion coefficients for charm sea quark distribution.
      DATA (((CEHLQ(IX,IT,NX,6,1),IX=1,6),IT=1,6),NX=1,2)/
     1 9.270E-03,-1.817E-02, 9.590E-03,-6.390E-03, 1.690E-03,-1.540E-03,
     2 5.710E-03,-1.188E-02, 6.090E-03,-4.650E-03, 1.240E-03,-1.310E-03,
     3-3.960E-03, 7.100E-03,-3.590E-03, 1.840E-03,-3.900E-04, 3.400E-04,
     4 1.120E-03,-1.960E-03, 1.120E-03,-4.800E-04, 1.000E-04,-4.000E-05,
     5 4.000E-05,-3.000E-05,-1.800E-04, 9.000E-05,-5.000E-05,-2.000E-05,
     6-4.200E-04, 7.300E-04,-1.600E-04, 5.000E-05, 5.000E-05, 5.000E-05,
     1 8.098E-01,-1.042E+00, 3.398E-01,-6.824E-02, 8.760E-03,-9.000E-04,
     2 8.961E-01,-1.217E+00, 4.339E-01,-9.287E-02, 1.304E-02,-1.290E-03,
     3 3.058E-02,-1.040E-01, 7.604E-02,-2.415E-02, 4.600E-03,-5.000E-04,
     4-2.451E-02, 4.432E-02,-1.651E-02, 1.430E-03, 1.200E-04,-1.000E-04,
     5 1.122E-02,-1.457E-02, 2.680E-03, 5.800E-04,-1.200E-04, 3.000E-05,
     6-7.730E-03, 7.330E-03,-7.600E-04,-2.400E-04, 1.000E-05, 0.000E+00/
      DATA (((CEHLQ(IX,IT,NX,6,2),IX=1,6),IT=1,6),NX=1,2)/
     1 9.980E-03,-1.945E-02, 1.055E-02,-6.870E-03, 1.860E-03,-1.560E-03,
     2 5.700E-03,-1.203E-02, 6.250E-03,-4.860E-03, 1.310E-03,-1.370E-03,
     3-4.490E-03, 7.990E-03,-4.170E-03, 2.050E-03,-4.400E-04, 3.300E-04,
     4 1.470E-03,-2.480E-03, 1.460E-03,-5.700E-04, 1.200E-04,-1.000E-05,
     5-9.000E-05, 1.500E-04,-3.200E-04, 1.200E-04,-6.000E-05,-4.000E-05,
     6-4.200E-04, 7.600E-04,-1.400E-04, 4.000E-05, 7.000E-05, 5.000E-05,
     1 8.698E-01,-1.131E+00, 3.836E-01,-8.111E-02, 1.048E-02,-1.300E-03,
     2 9.626E-01,-1.321E+00, 4.854E-01,-1.091E-01, 1.583E-02,-1.700E-03,
     3 3.057E-02,-1.088E-01, 8.022E-02,-2.676E-02, 5.590E-03,-5.600E-04,
     4-2.845E-02, 5.164E-02,-1.918E-02, 2.210E-03,-4.000E-05,-1.500E-04,
     5 1.311E-02,-1.751E-02, 3.310E-03, 5.100E-04,-1.200E-04, 5.000E-05,
     6-8.590E-03, 8.380E-03,-9.200E-04,-2.600E-04, 1.000E-05,-1.000E-05/
C...Expansion coefficients for bottom sea quark distribution.
      DATA (((CEHLQ(IX,IT,NX,7,1),IX=1,6),IT=1,6),NX=1,2)/
     1 9.010E-03,-1.401E-02, 7.150E-03,-4.130E-03, 1.260E-03,-1.040E-03,
     2 6.280E-03,-9.320E-03, 4.780E-03,-2.890E-03, 9.100E-04,-8.200E-04,
     3-2.930E-03, 4.090E-03,-1.890E-03, 7.600E-04,-2.300E-04, 1.400E-04,
     4 3.900E-04,-1.200E-03, 4.400E-04,-2.500E-04, 2.000E-05,-2.000E-05,
     5 2.600E-04, 1.400E-04,-8.000E-05, 1.000E-04, 1.000E-05, 1.000E-05,
     6-2.600E-04, 3.200E-04, 1.000E-05,-1.000E-05, 1.000E-05,-1.000E-05,
     1 8.029E-01,-1.075E+00, 3.792E-01,-7.843E-02, 1.007E-02,-1.090E-03,
     2 7.903E-01,-1.099E+00, 4.153E-01,-9.301E-02, 1.317E-02,-1.410E-03,
     3-1.704E-02,-1.130E-02, 2.882E-02,-1.341E-02, 3.040E-03,-3.600E-04,
     4-7.200E-04, 7.230E-03,-5.160E-03, 1.080E-03,-5.000E-05,-4.000E-05,
     5 3.050E-03,-4.610E-03, 1.660E-03,-1.300E-04,-1.000E-05, 1.000E-05,
     6-4.360E-03, 5.230E-03,-1.610E-03, 2.000E-04,-2.000E-05, 0.000E+00/
      DATA (((CEHLQ(IX,IT,NX,7,2),IX=1,6),IT=1,6),NX=1,2)/
     1 8.980E-03,-1.459E-02, 7.510E-03,-4.410E-03, 1.310E-03,-1.070E-03,
     2 5.970E-03,-9.440E-03, 4.800E-03,-3.020E-03, 9.100E-04,-8.500E-04,
     3-3.050E-03, 4.440E-03,-2.100E-03, 8.500E-04,-2.400E-04, 1.400E-04,
     4 5.300E-04,-1.300E-03, 5.600E-04,-2.700E-04, 3.000E-05,-2.000E-05,
     5 2.000E-04, 1.400E-04,-1.100E-04, 1.000E-04, 0.000E+00, 0.000E+00,
     6-2.600E-04, 3.200E-04, 0.000E+00,-3.000E-05, 1.000E-05,-1.000E-05,
     1 8.672E-01,-1.174E+00, 4.265E-01,-9.252E-02, 1.244E-02,-1.460E-03,
     2 8.500E-01,-1.194E+00, 4.630E-01,-1.083E-01, 1.614E-02,-1.830E-03,
     3-2.241E-02,-5.630E-03, 2.815E-02,-1.425E-02, 3.520E-03,-4.300E-04,
     4-7.300E-04, 8.030E-03,-5.780E-03, 1.380E-03,-1.300E-04,-4.000E-05,
     5 3.460E-03,-5.380E-03, 1.960E-03,-2.100E-04, 1.000E-05, 1.000E-05,
     6-4.850E-03, 5.950E-03,-1.890E-03, 2.600E-04,-3.000E-05, 0.000E+00/
C...Expansion coefficients for top sea quark distribution.
      DATA (((CEHLQ(IX,IT,NX,8,1),IX=1,6),IT=1,6),NX=1,2)/
     1 4.410E-03,-7.480E-03, 3.770E-03,-2.580E-03, 7.300E-04,-7.100E-04,
     2 3.840E-03,-6.050E-03, 3.030E-03,-2.030E-03, 5.800E-04,-5.900E-04,
     3-8.800E-04, 1.660E-03,-7.500E-04, 4.700E-04,-1.000E-04, 1.000E-04,
     4-8.000E-05,-1.500E-04, 1.200E-04,-9.000E-05, 3.000E-05, 0.000E+00,
     5 1.300E-04,-2.200E-04,-2.000E-05,-2.000E-05,-2.000E-05,-2.000E-05,
     6-7.000E-05, 1.900E-04,-4.000E-05, 2.000E-05, 0.000E+00, 0.000E+00,
     1 6.623E-01,-9.248E-01, 3.519E-01,-7.930E-02, 1.110E-02,-1.180E-03,
     2 6.380E-01,-9.062E-01, 3.582E-01,-8.479E-02, 1.265E-02,-1.390E-03,
     3-2.581E-02, 2.125E-02, 4.190E-03,-4.980E-03, 1.490E-03,-2.100E-04,
     4 7.100E-04, 5.300E-04,-1.270E-03, 3.900E-04,-5.000E-05,-1.000E-05,
     5 3.850E-03,-5.060E-03, 1.860E-03,-3.500E-04, 4.000E-05, 0.000E+00,
     6-3.530E-03, 4.460E-03,-1.500E-03, 2.700E-04,-3.000E-05, 0.000E+00/
      DATA (((CEHLQ(IX,IT,NX,8,2),IX=1,6),IT=1,6),NX=1,2)/
     1 4.260E-03,-7.530E-03, 3.830E-03,-2.680E-03, 7.600E-04,-7.300E-04,
     2 3.640E-03,-6.050E-03, 3.030E-03,-2.090E-03, 5.900E-04,-6.000E-04,
     3-9.200E-04, 1.710E-03,-8.200E-04, 5.000E-04,-1.200E-04, 1.000E-04,
     4-5.000E-05,-1.600E-04, 1.300E-04,-9.000E-05, 3.000E-05, 0.000E+00,
     5 1.300E-04,-2.100E-04,-1.000E-05,-2.000E-05,-2.000E-05,-1.000E-05,
     6-8.000E-05, 1.800E-04,-5.000E-05, 2.000E-05, 0.000E+00, 0.000E+00,
     1 7.146E-01,-1.007E+00, 3.932E-01,-9.246E-02, 1.366E-02,-1.540E-03,
     2 6.856E-01,-9.828E-01, 3.977E-01,-9.795E-02, 1.540E-02,-1.790E-03,
     3-3.053E-02, 2.758E-02, 2.150E-03,-4.880E-03, 1.640E-03,-2.500E-04,
     4 9.200E-04, 4.200E-04,-1.340E-03, 4.600E-04,-8.000E-05,-1.000E-05,
     5 4.230E-03,-5.660E-03, 2.140E-03,-4.300E-04, 6.000E-05, 0.000E+00,
     6-3.890E-03, 5.000E-03,-1.740E-03, 3.300E-04,-4.000E-05, 0.000E+00/

C...The following data lines are coefficients needed in the
C...Duke, Owens proton structure function parametrizations, see below.
C...Expansion coefficients for (up+down) valence quark distribution.
      DATA ((CDO(IP,IS,1,1),IS=1,6),IP=1,3)/
     1 4.190E-01, 3.460E+00, 4.400E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     2 4.000E-03, 7.240E-01,-4.860E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     3-7.000E-03,-6.600E-02, 1.330E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA ((CDO(IP,IS,1,2),IS=1,6),IP=1,3)/
     1 3.740E-01, 3.330E+00, 6.030E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     2 1.400E-02, 7.530E-01,-6.220E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     3 0.000E+00,-7.600E-02, 1.560E+00, 0.000E+00, 0.000E+00, 0.000E+00/
C...Expansion coefficients for down valence quark distribution.
      DATA ((CDO(IP,IS,2,1),IS=1,6),IP=1,3)/
     1 7.630E-01, 4.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     2-2.370E-01, 6.270E-01,-4.210E-01, 0.000E+00, 0.000E+00, 0.000E+00,
     3 2.600E-02,-1.900E-02, 3.300E-02, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA ((CDO(IP,IS,2,2),IS=1,6),IP=1,3)/
     1 7.610E-01, 3.830E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     2-2.320E-01, 6.270E-01,-4.180E-01, 0.000E+00, 0.000E+00, 0.000E+00,
     3 2.300E-02,-1.900E-02, 3.600E-02, 0.000E+00, 0.000E+00, 0.000E+00/
C...Expansion coefficients for (up+down+strange) sea quark distribution.
      DATA ((CDO(IP,IS,3,1),IS=1,6),IP=1,3)/
     1 1.265E+00, 0.000E+00, 8.050E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     2-1.132E+00,-3.720E-01, 1.590E+00, 6.310E+00,-1.050E+01, 1.470E+01,
     3 2.930E-01,-2.900E-02,-1.530E-01,-2.730E-01,-3.170E+00, 9.800E+00/
      DATA ((CDO(IP,IS,3,2),IS=1,6),IP=1,3)/
     1 1.670E+00, 0.000E+00, 9.150E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     2-1.920E+00,-2.730E-01, 5.300E-01, 1.570E+01,-1.010E+02, 2.230E+02,
     3 5.820E-01,-1.640E-01,-7.630E-01,-2.830E+00, 4.470E+01,-1.170E+02/
C...Expansion coefficients for charm sea quark distribution.
      DATA ((CDO(IP,IS,4,1),IS=1,6),IP=1,3)/
     1 0.000E+00,-3.600E-02, 6.350E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     2 1.350E-01,-2.220E-01, 3.260E+00,-3.030E+00, 1.740E+01,-1.790E+01,
     3-7.500E-02,-5.800E-02,-9.090E-01, 1.500E+00,-1.130E+01, 1.560E+01/
       DATA ((CDO(IP,IS,4,2),IS=1,6),IP=1,3)/
     1 0.000E+00,-1.200E-01, 3.510E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     2 6.700E-02,-2.330E-01, 3.660E+00,-4.740E-01, 9.500E+00,-1.660E+01,
     3-3.100E-02,-2.300E-02,-4.530E-01, 3.580E-01,-5.430E+00, 1.550E+01/
C...Expansion coefficients for gluon distribution.
      DATA ((CDO(IP,IS,5,1),IS=1,6),IP=1,3)/
     1 1.560E+00, 0.000E+00, 6.000E+00, 9.000E+00, 0.000E+00, 0.000E+00,
     2-1.710E+00,-9.490E-01, 1.440E+00,-7.190E+00,-1.650E+01, 1.530E+01,
     3 6.380E-01, 3.250E-01,-1.050E+00, 2.550E-01, 1.090E+01,-1.010E+01/
      DATA ((CDO(IP,IS,5,2),IS=1,6),IP=1,3)/
     1 8.790E-01, 0.000E+00, 4.000E+00, 9.000E+00, 0.000E+00, 0.000E+00,
     2-9.710E-01,-1.160E+00, 1.230E+00,-5.640E+00,-7.540E+00,-5.960E-01,
     3 4.340E-01, 4.760E-01,-2.540E-01,-8.170E-01, 5.500E+00, 1.260E-01/

C...The following data lines are coefficients needed in the
C...Morfin and Tung structure function parametrizations.
C...12 coefficients each for d(valence), u(valence), g, u(sea),
C...d(sea), s, c, b and t, in that order.
C...Expansion coefficients for set 1 (fit S1).
      DATA (((CMT(IEX,IPN,IFL,1),IFL=1,9),IPN=0,2),IEX=0,3)/
     &   1.30,  1.64,  1.86, -0.60, -0.45, -1.10, -3.87, -6.14,-12.53,
     &  -0.57, -0.33, -2.76, -1.68, -1.64, -1.66,  0.79,  2.65,  8.13,
     &  -0.08, -0.10,  0.10,  0.08,  0.05,  0.13, -0.70, -1.24, -2.64,
     &   0.18,  0.08, -0.17, -0.19, -0.18, -0.19, -0.03, -0.10, -0.38,
     &   0.16,  0.14, -0.07, -0.16, -0.19, -0.09, -0.17, -0.03,  0.34,
     &  -0.02, -0.01,  0.02,  0.04,  0.06,  0.01,  0.03, -0.02, -0.14,
     &   5.27,  3.74,  7.33,  9.31,  9.36,  9.07,  7.96,  6.90, 16.30,
     &   0.43,  0.54, -0.88, -1.17, -1.01, -1.39,  0.95,  1.52,-13.23,
     &   0.06,  0.03, -0.08,  0.29,  0.20,  0.47, -0.38, -0.50,  4.77,
     &  -1.85, -2.04, -0.88, -1.45, -1.48, -1.26,  0.60,  0.80, -0.57,
     &   1.08,  0.88,  2.47,  1.65,  1.49,  1.96,  0.60,  1.05,  3.58,
     &  -0.03,  0.02, -0.32, -0.20, -0.12, -0.36,  0.08, -0.14, -0.99/
C...Expansion coefficients for set 2 (fit B1).
      DATA (((CMT(IEX,IPN,IFL,2),IFL=1,9),IPN=0,2),IEX=0,3)/
     &   1.34,  1.62,  1.88, -0.99, -0.99, -0.99, -3.98, -6.28,-13.08,
     &  -0.57, -0.33, -2.78, -1.54, -1.54, -1.54,  0.72,  2.62,  8.54,
     &  -0.08, -0.10,  0.13,  0.10,  0.10,  0.10, -0.63, -1.18, -2.70,
     &   0.15,  0.11, -0.33, -0.33, -0.33, -0.33, -0.15, -0.18, -0.40,
     &   0.16,  0.14,  0.10,  0.03,  0.03,  0.03, -0.06,  0.02,  0.31,
     &  -0.02, -0.01, -0.04, -0.03, -0.03, -0.03,  0.00, -0.03, -0.12,
     &   5.30,  3.68,  7.52,  8.53,  8.53,  8.53,  7.46,  6.56, 15.35,
     &   0.43,  0.53, -1.13, -1.08, -1.08, -1.08,  0.96,  1.40,-11.83,
     &   0.06,  0.03,  0.04,  0.39,  0.39,  0.39, -0.30, -0.38,  4.16,
     &  -1.96, -1.94, -1.34, -1.55, -1.55, -1.55,  0.35,  0.65, -0.43,
     &   1.08,  0.87,  2.92,  2.02,  2.02,  2.02,  0.89,  1.13,  3.18,
     &  -0.03,  0.02, -0.49, -0.39, -0.39, -0.39, -0.04, -0.16, -0.82/
C...Expansion coefficients for set 3 (fit B2).
      DATA (((CMT(IEX,IPN,IFL,3),IFL=1,9),IPN=0,2),IEX=0,3)/
     &   1.38,  1.64,  1.52, -0.85, -0.85, -0.85, -3.74, -6.07,-12.08,
     &  -0.59, -0.33, -2.71, -1.43, -1.43, -1.43,  0.21,  2.33,  7.31,
     &  -0.08, -0.10,  0.15, -0.03, -0.03, -0.03, -0.50, -1.15, -2.35,
     &   0.18,  0.09, -0.72, -0.82, -0.82, -0.82, -0.58, -0.52, -0.73,
     &   0.16,  0.14,  0.45,  0.35,  0.35,  0.35,  0.24,  0.22,  0.54,
     &  -0.02, -0.01, -0.15, -0.09, -0.10, -0.10, -0.07, -0.07, -0.18,
     &   5.40,  3.74,  7.75,  9.19,  9.19,  9.19,  9.63,  8.33, 21.14,
     &   0.42,  0.54, -1.56, -0.92, -0.92, -0.92, -1.13,  0.28,-19.17,
     &   0.06,  0.03,  0.16,  0.12,  0.12,  0.12,  0.25, -0.28,  6.64,
     &  -1.91, -2.02, -2.18, -2.76, -2.76, -2.76, -1.09, -0.52, -1.92,
     &   1.11,  0.88,  3.75,  2.56,  2.56,  2.56,  2.10,  1.91,  4.59,
     &  -0.03,  0.02, -0.76, -0.40, -0.40, -0.40, -0.33, -0.31, -1.25/
C...Expansion coefficients for set 4 (fit E1).
      DATA (((CMT(IEX,IPN,IFL,4),IFL=1,9),IPN=0,2),IEX=0,3)/
     &   1.43,  1.69,  2.11, -0.84, -0.84, -0.84, -3.87, -6.09,-12.56,
     &  -0.65, -0.33, -3.01, -1.65, -1.65, -1.65,  0.85,  2.81,  8.69,
     &  -0.08, -0.11,  0.18,  0.12,  0.12,  0.12, -0.73, -1.34, -2.93,
     &   0.16,  0.11, -0.33, -0.32, -0.32, -0.32, -0.15, -0.17, -0.38,
     &   0.16,  0.14,  0.10,  0.02,  0.02,  0.02, -0.07,  0.01,  0.30,
     &  -0.02, -0.01, -0.04, -0.03, -0.03, -0.03,  0.00, -0.03, -0.12,
     &   6.17,  3.69,  7.93,  8.96,  8.96,  8.96,  7.83,  6.75, 14.62,
     &   0.43,  0.54, -1.40, -1.24, -1.24, -1.24,  1.00,  1.74,-11.27,
     &   0.06,  0.03,  0.09,  0.45,  0.45,  0.45, -0.36, -0.56,  4.29,
     &  -1.94, -1.99, -1.51, -1.70, -1.70, -1.70,  0.21,  0.54, -0.41,
     &   1.12,  0.90,  3.14,  2.15,  2.15,  2.15,  0.93,  1.15,  3.19,
     &  -0.02,  0.02, -0.55, -0.43, -0.43, -0.43, -0.03, -0.16, -0.87/

C...Euler's beta function, requires ordinary Gamma function
      EULBET(X,Y)=RYGAMM(X)*RYGAMM(Y)/RYGAMM(X+Y)

C...Reset structure functions, check x and hadron flavour.
      ALAM=0.
      DO 100 KFL=-25,25
  100 XPQ(KFL)=0.
      IF(X.LE.0..OR.X.GE.1.) THEN
        WRITE(MSTU(11),5000) X
        RETURN
      ENDIF
      KFA=IABS(KF)
      IF(KFA.NE.11.AND.KFA.NE.22.AND.KFA.NE.211.AND.KFA.NE.2112.AND.
     &KFA.NE.2212) THEN
        WRITE(MSTU(11),5100) KF
        RETURN
      ENDIF

C...Call user-supplied structure function.
      IF(MSTP(51).EQ.0.OR.MSTP(52).GE.2) THEN
        KFE=KFA
        IF(KFA.EQ.2112) KFE=2212
        CALL RYSTFE(KFE,X,Q2,XPQ)

      ELSEIF(KFA.EQ.11) THEN
C...Electron structure function.
        AEM=PARU(101)
        PME=PMAS(11,1)
        XL=LOG(MAX(1E-10,X))
        X1L=LOG(MAX(1E-10,1.-X))
        HLE=LOG(MAX(3.,Q2/PME**2))
        HBE=(2.*AEM/PARU(1))*(HLE-1.)

C...Electron inside electron, see R. Kleiss et al., in Z physics at
C...LEP 1, CERN 89-08, p. 34
        IF(MSTP(11).LE.1) THEN
          HDE=1.+(AEM/PARU(1))*(1.5*HLE+1.289868)+(AEM/PARU(1))**2*
     &    (-2.164868*HLE**2+9.840808*HLE-10.130464)
          HEE=0.5*HBE*(1.-X)**(0.5*HBE-1.)*SQRT(MAX(0.,HDE))-
     &    0.25*HBE*(1.+X)+HBE**2/32.*((1.+X)*(-4.*X1L+3.*XL)-
     &    4.*XL/(1.-X)-5.-X)
          HCB=0.5*HBE
        ELSE
          HCA=PARP(11)
          HCB=PARP(12)
          IF(MSTP(11).EQ.3) HCB=HCB+0.5*HBE
          HEE=X**HCA*(1.-X)**HCB/EULBET(1.+HCA,1.+HCB)
        ENDIF
        IF(X.GT.0.9999.AND.X.LE.0.999999) THEN
          HEE=HEE*100.**HCB/(100.**HCB-1.)
        ELSEIF(X.GT.0.999999) THEN
          HEE=0.
        ENDIF
        XPQ(11)=X*HEE

C...Photon and (transverse) W- inside electron.
        AEMP=ULALEM(PME*SQRT(MAX(0.,Q2)))/PARU(2)
        IF(MSTP(13).LE.1) THEN
          HLG=HLE
        ELSE
          HLG=LOG((PARP(13)/PME**2)*(1.-X)/X**2)
        ENDIF
        XPQ(22)=AEMP*HLG*(1.+(1.-X)**2)
        HLW=LOG(1.+Q2/PMAS(24,1)**2)/(4.*PARU(102))
        XPQ(-24)=AEMP*HLW*(1.+(1.-X)**2)

C..Quarks and gluons inside photon inside electron.
        IF(MSTP(12).EQ.1) THEN
          T=LOG(MIN(1E4,MAX(1.,Q2))/0.16)
          NF=3
          IF(Q2.GT.25.) NF=4
          IF(Q2.GT.300.) NF=5
          NFE=NF-2
          XL=LOG(MAX(1E-10,X))

C...Numerical integration of struncture function convolution.
          SXPGL=0.
          SXPQU=0.
          SXPQD=0.
          SUMXPP=0.
          ITER=-1
  110     ITER=ITER+1
          SUMXP=SUMXPP
          NSTP=2**(ITER-1)
          IF(ITER.EQ.0) NSTP=2
          SXPGL=0.5*SXPGL
          SXPQU=0.5*SXPQU
          SXPQD=0.5*SXPQD
          WTSTP=0.5/NSTP
          IF(ITER.EQ.0) WTSTP=0.5
          DO 120 ISTP=1,NSTP
          IF(ITER.EQ.0) THEN
            XLE=XL*(ISTP-1)
          ELSE
            XLE=XL*(ISTP-0.5)/NSTP
          ENDIF
          XE=EXP(XLE)
          XG=MIN(0.999999,X/XE)
          XPGA=1.+(1.-XE)**2
          CALL RYSTGA(NFE,XG,T,XPGL,XPQU,XPQD)
          SXPGL=SXPGL+WTSTP*XPGA*XPGL
          SXPQU=SXPQU+WTSTP*XPGA*XPQU
  120     SXPQD=SXPQD+WTSTP*XPGA*XPQD
          SUMXPP=SXPGL+SXPQU+SXPQD
          IF(ITER.LE.2.OR.(ITER.LE.7.AND.ABS(SUMXPP-SUMXP).GT.
     &    PARP(14)*(SUMXPP+SUMXP))) GOTO 110
          FCONV=AEMP*HLE*AEM*(-XL)

C...Put into output arrays.
          XPQ(0)=FCONV*SXPGL
          XPQ(1)=FCONV*SXPQD
          XPQ(-1)=XPQ(1)
          XPQ(2)=FCONV*SXPQU
          XPQ(-2)=XPQ(2)
          XPQ(3)=FCONV*SXPQD
          XPQ(-3)=XPQ(3)
          IF(NFE.GE.2) THEN
            XPQ(4)=FCONV*SXPQU
            XPQ(-4)=XPQ(4)
          ENDIF
          IF(NFE.EQ.3) THEN
            XPQ(5)=FCONV*SXPQD
            XPQ(-5)=XPQ(5)
          ENDIF
        ENDIF

      ELSEIF(KFA.EQ.22) THEN
C...Photon structure function from Drees and Grassie.
C...Allowed variable range: 1 GeV^2 < Q^2 < 10000 GeV^2.
        T=LOG(MIN(1E4,MAX(1.,Q2))/0.16)
        NF=3
        IF(Q2.GT.25.) NF=4
        IF(Q2.GT.300.) NF=5
        NFE=NF-2
        CALL RYSTGA(NFE,X,T,XPGL,XPQU,XPQD)
        AEM=PARU(101)

C...Put into output arrays.
        XPQ(0)=AEM*XPGL
        XPQ(1)=AEM*XPQD
        XPQ(-1)=XPQ(1)
        XPQ(2)=AEM*XPQU
        XPQ(-2)=XPQ(2)
        XPQ(3)=AEM*XPQD
        XPQ(-3)=XPQ(3)
        IF(NFE.GE.2) THEN
          XPQ(4)=AEM*XPQU
          XPQ(-4)=XPQ(4)
        ENDIF
        IF(NFE.EQ.3) THEN
          XPQ(5)=AEM*XPQD
          XPQ(-5)=XPQ(5)
        ENDIF

      ELSEIF(KFA.EQ.211) THEN
C...Pion structure functions from Owens.
C...Allowed variable range: 4 GeV^2 < Q^2 < approx 2000 GeV^2.

C...Determine set, Lambda and s expansion variable.
        NSET=1
        IF(MSTP(51).EQ.2.OR.MSTP(51).EQ.4.OR.MSTP(51).EQ.13) NSET=2
        IF(NSET.EQ.1) ALAM=0.2
        IF(NSET.EQ.2) ALAM=0.4
        IF(MSTP(52).LE.0) THEN
          SD=0.
        ELSE
          Q2IN=MIN(2E3,MAX(4.,Q2))
          SD=LOG(LOG(Q2IN/ALAM**2)/LOG(4./ALAM**2))
        ENDIF

C...Calculate structure functions.
        DO 140 KFL=1,4
        DO 130 IS=1,5
  130   TS(IS)=COW(1,IS,KFL,NSET)+COW(2,IS,KFL,NSET)*SD+
     &  COW(3,IS,KFL,NSET)*SD**2
        IF(KFL.EQ.1) THEN
          XQ(KFL)=X**TS(1)*(1.-X)**TS(2)/EULBET(TS(1),TS(2)+1.)
        ELSE
          XQ(KFL)=TS(1)*X**TS(2)*(1.-X)**TS(3)*(1.+TS(4)*X+TS(5)*X**2)
        ENDIF
  140   CONTINUE

C...Put into output arrays.
        XPQ(0)=XQ(2)
        XPQ(1)=XQ(3)/6.
        XPQ(2)=XQ(1)+XQ(3)/6.
        XPQ(3)=XQ(3)/6.
        XPQ(4)=XQ(4)
        XPQ(-1)=XQ(1)+XQ(3)/6.
        XPQ(-2)=XQ(3)/6.
        XPQ(-3)=XQ(3)/6.
        XPQ(-4)=XQ(4)

      ELSEIF(MSTP(51).EQ.1.OR.MSTP(51).EQ.2) THEN
C...Proton structure functions from Eichten, Hinchliffe, Lane, Quigg.
C...Allowed variable range: 5 GeV^2 < Q^2 < 1E8 GeV^2; 1E-4 < x < 1

C...Determine set, Lamdba and x and t expansion variables.
        NSET=MSTP(51)
        IF(NSET.EQ.1) ALAM=0.2
        IF(NSET.EQ.2) ALAM=0.29
        TMIN=LOG(5./ALAM**2)
        TMAX=LOG(1E8/ALAM**2)
        IF(MSTP(52).EQ.0) THEN
          T=TMIN
        ELSE
          T=LOG(MAX(1.,Q2/ALAM**2))
        ENDIF
        VT=MAX(-1.,MIN(1.,(2.*T-TMAX-TMIN)/(TMAX-TMIN)))
        NX=1
        IF(X.LE.0.1) NX=2
        IF(NX.EQ.1) VX=(2.*X-1.1)/0.9
        IF(NX.EQ.2) VX=MAX(-1.,(2.*LOG(X)+11.51293)/6.90776)
        CXS=1.
        IF(X.LT.1E-4.AND.ABS(PARP(51)-1.).GT.0.01) CXS=
     &  (1E-4/X)**(PARP(51)-1.)

C...Chebyshev polynomials for x and t expansion.
        TX(1)=1.
        TX(2)=VX
        TX(3)=2.*VX**2-1.
        TX(4)=4.*VX**3-3.*VX
        TX(5)=8.*VX**4-8.*VX**2+1.
        TX(6)=16.*VX**5-20.*VX**3+5.*VX
        TT(1)=1.
        TT(2)=VT
        TT(3)=2.*VT**2-1.
        TT(4)=4.*VT**3-3.*VT
        TT(5)=8.*VT**4-8.*VT**2+1.
        TT(6)=16.*VT**5-20.*VT**3+5.*VT

C...Calculate structure functions.
        DO 160 KFL=1,6
        XQSUM=0.
        DO 150 IT=1,6
        DO 150 IX=1,6
  150   XQSUM=XQSUM+CEHLQ(IX,IT,NX,KFL,NSET)*TX(IX)*TT(IT)
  160   XQ(KFL)=XQSUM*(1.-X)**NEHLQ(KFL,NSET)*CXS

C...Put into output array.
        XPQ(0)=XQ(4)
        XPQ(1)=XQ(2)+XQ(3)
        XPQ(2)=XQ(1)+XQ(3)
        XPQ(3)=XQ(5)
        XPQ(4)=XQ(6)
        XPQ(-1)=XQ(3)
        XPQ(-2)=XQ(3)
        XPQ(-3)=XQ(5)
        XPQ(-4)=XQ(6)

C...Special expansion for bottom (threshold effects).
        IF(MSTP(54).GE.5) THEN
          IF(NSET.EQ.1) TMIN=8.1905
          IF(NSET.EQ.2) TMIN=7.4474
          IF(T.LE.TMIN) GOTO 180
          VT=MAX(-1.,MIN(1.,(2.*T-TMAX-TMIN)/(TMAX-TMIN)))
          TT(1)=1.
          TT(2)=VT
          TT(3)=2.*VT**2-1.
          TT(4)=4.*VT**3-3.*VT
          TT(5)=8.*VT**4-8.*VT**2+1.
          TT(6)=16.*VT**5-20.*VT**3+5.*VT
          XQSUM=0.
          DO 170 IT=1,6
          DO 170 IX=1,6
  170     XQSUM=XQSUM+CEHLQ(IX,IT,NX,7,NSET)*TX(IX)*TT(IT)
          XPQ(5)=XQSUM*(1.-X)**NEHLQ(7,NSET)*CXS
          XPQ(-5)=XPQ(5)
  180     CONTINUE
        ENDIF

C...Special expansion for top (threshold effects).
        IF(MSTP(54).GE.6) THEN
          IF(NSET.EQ.1) TMIN=11.5528
          IF(NSET.EQ.2) TMIN=10.8097
          TMIN=TMIN+2.*LOG(PMAS(6,1)/30.)
          TMAX=TMAX+2.*LOG(PMAS(6,1)/30.)
          IF(T.LE.TMIN) GOTO 200
          VT=MAX(-1.,MIN(1.,(2.*T-TMAX-TMIN)/(TMAX-TMIN)))
          TT(1)=1.
          TT(2)=VT
          TT(3)=2.*VT**2-1.
          TT(4)=4.*VT**3-3.*VT
          TT(5)=8.*VT**4-8.*VT**2+1.
          TT(6)=16.*VT**5-20.*VT**3+5.*VT
          XQSUM=0.
          DO 190 IT=1,6
          DO 190 IX=1,6
  190     XQSUM=XQSUM+CEHLQ(IX,IT,NX,8,NSET)*TX(IX)*TT(IT)
          XPQ(6)=XQSUM*(1.-X)**NEHLQ(8,NSET)*CXS
          XPQ(-6)=XPQ(6)
  200     CONTINUE
        ENDIF

      ELSEIF(MSTP(51).EQ.3.OR.MSTP(51).EQ.4) THEN
C...Proton structure functions from Duke, Owens.
C...Allowed variable range: 4 GeV^2 < Q^2 < approx 1E6 GeV^2.

C...Determine set, Lambda and s expansion parameter.
        NSET=MSTP(51)-2
        IF(NSET.EQ.1) ALAM=0.2
        IF(NSET.EQ.2) ALAM=0.4
        IF(MSTP(52).LE.0) THEN
          SD=0.
        ELSE
          Q2IN=MIN(1E6,MAX(4.,Q2))
          SD=LOG(LOG(Q2IN/ALAM**2)/LOG(4./ALAM**2))
        ENDIF

C...Calculate structure functions.
        DO 220 KFL=1,5
        DO 210 IS=1,6
  210   TS(IS)=CDO(1,IS,KFL,NSET)+CDO(2,IS,KFL,NSET)*SD+
     &  CDO(3,IS,KFL,NSET)*SD**2
        IF(KFL.LE.2) THEN
          XQ(KFL)=X**TS(1)*(1.-X)**TS(2)*(1.+TS(3)*X)/(EULBET(TS(1),
     &    TS(2)+1.)*(1.+TS(3)*TS(1)/(TS(1)+TS(2)+1.)))
        ELSE
          XQ(KFL)=TS(1)*X**TS(2)*(1.-X)**TS(3)*(1.+TS(4)*X+TS(5)*X**2+
     &    TS(6)*X**3)
        ENDIF
  220   CONTINUE

C...Put into output arrays.
        XPQ(0)=XQ(5)
        XPQ(1)=XQ(2)+XQ(3)/6.
        XPQ(2)=3.*XQ(1)-XQ(2)+XQ(3)/6.
        XPQ(3)=XQ(3)/6.
        XPQ(4)=XQ(4)
        XPQ(-1)=XQ(3)/6.
        XPQ(-2)=XQ(3)/6.
        XPQ(-3)=XQ(3)/6.
        XPQ(-4)=XQ(4)

      ELSEIF(MSTP(51).GE.5.AND.MSTP(51).LE.8) THEN
C...Proton structure functions from Morfin and Tung.
C...Allowed variable range: 4 GeV^2 < Q^2 < 1E8 GeV^2, 0 < x < 1.

C...Calculate expansion parameters.
        NSET=MSTP(51)-4
        IF(NSET.EQ.1) ALAM=0.187
        IF(NSET.EQ.2) ALAM=0.212
        IF(NSET.EQ.3) ALAM=0.191
        IF(NSET.EQ.4) ALAM=0.155
        IF(MSTP(52).EQ.0) THEN
          SD=0.
        ELSE
          SD=LOG(LOG(MAX(4.,Q2)/ALAM**2)/LOG(4./ALAM**2))
        ENDIF
        XL=LOG(MAX(1E-10,X))
        X1L=LOG(MAX(1E-10,1.-X))
        XLL=LOG(MAX(1E-10,LOG(1.+1./MAX(1E-10,X))))

C...Calculate structure functions up to b.
        DO 240 IP=1,8
        DO 230 IEX=0,3
  230   EXMT(IEX)=CMT(IEX,0,IP,NSET)+CMT(IEX,1,IP,NSET)*SD+
     &  CMT(IEX,2,IP,NSET)*SD**2
        EXMTS=EXMT(0)+EXMT(1)*XL+EXMT(2)*X1L+EXMT(3)*XLL
        IF(EXMTS.LT.-50.) THEN
          XQ(IP)=0.
        ELSE
          XQ(IP)=EXP(EXMTS)
        ENDIF
  240   CONTINUE
        IF(Q2.LE.2.25) XQ(7)=0
        IF(Q2.LE.25.0) XQ(8)=0

C...Calculate t structure function, shifting effective Q scale for
C...nondefault t mass, Q_actual = Q_nominal * m_t_nominal/m_t_actual.
        IF(MSTP(52).EQ.0.OR.Q2.LE.PMAS(6,1)**2) THEN
          XQ(9)=0.
        ELSE
          SD=LOG(LOG(MAX(4.,Q2)/ALAM**2*(100./PMAS(6,1))**2)/
     &    LOG(4./ALAM**2))
          DO 250 IEX=0,3
  250     EXMT(IEX)=CMT(IEX,0,9,NSET)+CMT(IEX,1,9,NSET)*SD+
     &    CMT(IEX,2,9,NSET)*SD**2
          EXMTS=EXMT(0)+EXMT(1)*XL+EXMT(2)*X1L+EXMT(3)*XLL
          IF(EXMTS.LT.-50.) THEN
            XQ(IP)=0.
          ELSE
            XQ(IP)=EXP(EXMTS)
          ENDIF
        ENDIF

C...Put into output array.
        XPQ(0)=XQ(3)
        XPQ(1)=XQ(1)+XQ(5)
        XPQ(-1)=XQ(5)
        XPQ(2)=XQ(2)+XQ(4)
        XPQ(-2)=XQ(4)
        XPQ(3)=XQ(6)
        XPQ(-3)=XQ(6)
        XPQ(4)=XQ(7)
        XPQ(-4)=XQ(7)
        XPQ(5)=XQ(8)
        XPQ(-5)=XQ(8)
        XPQ(6)=XQ(9)
        XPQ(-6)=XQ(9)

      ELSEIF(MSTP(51).EQ.9) THEN
C...Lowest order parametrization of Gluck, Reya, Vogt.
C...Allowed variable range: 0.2 GeV^2 < Q2 < 1E6 GeV^2; 1E-4 < x < 1;
C...extended to 0.2 GeV^2 < Q2 < 1E8 GeV^2; 1E-6 < x < 1
C...after consultation with the authors.

C...Determine s and x.
        ALAM=0.25
        IF(MSTP(52).EQ.0) THEN
          SD=0.
        ELSE
          Q2IN=MIN(1E8,MAX(0.2,Q2))
          SD=LOG(LOG(Q2IN/ALAM**2)/LOG(0.2/ALAM**2))
        ENDIF
        XC=MAX(1E-6,X)
        XL=-LOG(XC)

C...Calculate structure functions.
        XQ(1)=(0.794+0.312*SD)*XC**(0.427-0.011*SD)*
     &  (1.+(6.887-2.227*SD)*XC+(-11.083+2.136*SD)*XC**2+
     &  (3.900+1.079*SD)*XC**3)*(1.-XC)**(1.037+1.077*SD)
        XQ(2)=(0.486+0.139*SD)*XC**(0.434-0.018*SD)*
     &  (1.+(7.716-2.684*SD)*XC+(-12.768+3.122*SD)*XC**2+
     &  (4.564+0.470*SD)*XC**3)*(1.-XC)**(1.889+1.129*SD)
        XQ(3)=(XC**(0.415+0.186*SD)*((0.786+0.942*SD)+
     &  (5.256-5.810*SD)*XC+(-4.599+5.842*SD)*XC**2)+SD**0.592*
     &  EXP(-(0.398+2.135*SD)+SQRT(3.779*SD**1.250*XL)))*
     &  (1.-XC)**(1.622+1.980*SD)
        XQ(4)=SD**0.448*(1.-XC)**(5.540-0.445*SD)*
     &  EXP(-(4.668+1.230*SD)+SQRT((13.173-1.361*SD)*SD**0.442*XL))/
     &  XL**(3.181-0.862*SD)
        XQ(5)=0.
        IF(SD.GT.1.125) XQ(5)=(SD-1.125)*(1.-XC)**(2.038+1.022*SD)*
     &  EXP(-(4.290+1.569*SD)+SQRT((2.982+1.452*SD)*SD**0.5*XL))
        XQ(6)=0.
        IF(SD.GT.1.603) XQ(6)=(SD-1.603)*(1.-XC)**(2.230+1.052*SD)*
     &  EXP(-(4.566+1.559*SD)+SQRT((4.147+1.268*SD)*SD**0.5*XL))

C...Put into output array - special factor for small x.
        CXS=1.
        IF(X.LT.1E-6.AND.ABS(PARP(51)-1.).GT.0.01)
     &  CXS=(1E-6/X)**(PARP(51)-1.)
        XPQ(0)=CXS*XQ(3)
        XPQ(1)=CXS*(XQ(2)+XQ(4))
        XPQ(-1)=CXS*XQ(4)
        XPQ(2)=CXS*(XQ(1)+XQ(4))
        XPQ(-2)=CXS*XQ(4)
        XPQ(3)=CXS*XQ(4)
        XPQ(-3)=CXS*XQ(4)
        XPQ(4)=CXS*XQ(5)
        XPQ(-4)=CXS*XQ(5)
        XPQ(5)=CXS*XQ(6)
        XPQ(-5)=CXS*XQ(6)

      ELSEIF(MSTP(51).EQ.10) THEN
C...Higher order parametrization of Gluck, Reya, Vogt.
C...Allowed variable range: 0.2 GeV^2 < Q2 < 1E6 GeV^2; 1E-4 < x < 1;
C...extended to 0.2 GeV^2 < Q2 < 1E8 GeV^2; 1E-6 < x < 1
C...after consultation with the authors.

C...Determine s and x.
        ALAM=0.20
        IF(MSTP(52).EQ.0) THEN
          SD=0.
        ELSE
          Q2IN=MIN(1E8,MAX(0.2,Q2))
          SD=LOG(LOG(Q2IN/ALAM**2)/LOG(0.2/ALAM**2))
        ENDIF
        SD2=SD**2
        XC=MAX(1E-6,X)
        XL=-LOG(XC)

C...Calculate structure functions.
        XQ(1)=(1.364+0.989*SD-0.236*SD2)*XC**(0.593-0.048*SD)*
     &  (1.+(8.912-6.092*SD+0.852*SD2)*XC+(-16.737+7.039*SD)*XC**2+
     &  (10.275+0.806*SD-2.000*SD2)*XC**3)*
     &  (1.-XC)**(2.043+1.408*SD-0.283*SD2)
        XQ(2)=(0.835+0.527*SD-0.144*SD2)*XC**(0.600-0.054*SD)*
     &  (1.+(10.245-7.821*SD+1.325*SD2)*XC+(-19.511+10.940*SD-
     &  1.133*SD2)*XC**2+(12.836-2.570*SD-1.041*SD2)*XC**3)*
     &  (1.-XC)**(3.083+1.382*SD-0.276*SD2)
        XQ(3)=(XC**(0.321-0.135*SD)*((10.51-2.299*SD)+
     &  (-17.28+0.755*SD)*XC+(8.242+2.543*SD)*XC**2)*
     &  XL**(-2.023-0.103*SD)+SD**1.044*
     &  EXP(-(-1.178+2.792*SD)+SQRT(2.318*SD**1.673*XL)))*
     &  (1.-XC)**(3.720+2.337*SD-0.199*SD2)
        XQ(4)=SD**0.761*(1.+(6.078-2.065*SD)*XC)*(1.-XC)**(4.654+
     &  0.603*SD-0.326*SD2)*EXP(-(4.231+1.036*SD)+SQRT(3.419*SD**0.316*
     &  XL))/XL**(0.897-0.618*SD)
        XQ(5)=0.
        IF(SD.GT.0.918) XQ(5)=(SD-0.918)*(1.-XC)**(3.328+0.859*SD)*
     &  EXP(-(3.837+1.504*SD)+SQRT((2.150+1.291*SD)*SD**0.5*XL))
        XQ(6)=0.
        IF(SD.GT.1.353) XQ(6)=(SD-1.353)*(1.-XC)**(3.382+0.909*SD)*
     &  EXP(-(4.130+1.486*SD)+SQRT((2.895+1.240*SD)*SD**0.5*XL))

C...Put into output array - special factor for small x.
        CXS=1.
        IF(X.LT.1E-6.AND.ABS(PARP(51)-1.).GT.0.01)
     &  CXS=(1E-6/X)**(PARP(51)-1.)
        XPQ(0)=CXS*XQ(3)
        XPQ(1)=CXS*(XQ(2)+XQ(4))
        XPQ(-1)=CXS*XQ(4)
        XPQ(2)=CXS*(XQ(1)+XQ(4))
        XPQ(-2)=CXS*XQ(4)
        XPQ(3)=CXS*XQ(4)
        XPQ(-3)=CXS*XQ(4)
        XPQ(4)=CXS*XQ(5)
        XPQ(-4)=CXS*XQ(5)
        XPQ(5)=CXS*XQ(6)
        XPQ(-5)=CXS*XQ(6)

C...Proton structure functions from Diemoz, Ferroni, Longo, Martinelli.
C...These are accessed via RYSTFE since the files needed may not always
C...available.
      ELSEIF(MSTP(51).GE.11.AND.MSTP(51).LE.13) THEN
        CALL RYSTFE(2212,X,Q2,XPQ)

C...Unknown proton parametrization.
      ELSE
        WRITE(MSTU(11),5200) MSTP(51)
      ENDIF

C...Isospin conjugation for neutron.
      IF(KFA.EQ.2112) THEN
        XPS=XPQ(1)
        XPQ(1)=XPQ(2)
        XPQ(2)=XPS
        XPS=XPQ(-1)
        XPQ(-1)=XPQ(-2)
        XPQ(-2)=XPS
      ENDIF

C...Charge conjugation for antiparticle.
      IF(KF.LT.0) THEN
        DO 260 KFL=1,25
        IF(KFL.EQ.21.OR.KFL.EQ.22.OR.KFL.EQ.23.OR.KFL.EQ.25) GOTO 260
        XPS=XPQ(KFL)
        XPQ(KFL)=XPQ(-KFL)
        XPQ(-KFL)=XPS
  260   CONTINUE
      ENDIF

C...Check positivity and reset above maximum allowed flavour.
      DO 270 KFL=-25,25
      XPQ(KFL)=MAX(0.,XPQ(KFL))
  270 IF(IABS(KFL).LE.8.AND.IABS(KFL).GT.MSTP(54)) XPQ(KFL)=0.

C...Formats for error printouts.
 5000 FORMAT(' Error: x value outside physical range, x =',1P,E12.3)
 5100 FORMAT(' Error: illegal particle code for structure function,',
     &' KF =',I5)
 5200 FORMAT(' Error: bad value of parameter MSTP(51) in RYSTFU,',
     &' MSTP(51) =',I5)

      RETURN
      END

C*********************************************************************

      SUBROUTINE RYSTGA(NFE,X,T,XPGL,XPQU,XPQD)

C...Gives photon structure function; external to RYSTFU since it
C...may be called several times for convolution of photon structure
C...functions with photons in electron structure function.
      DIMENSION DGAG(4,3),DGBG(4,3),DGCG(4,3),DGAN(4,3),DGBN(4,3),
     &DGCN(4,3),DGDN(4,3),DGEN(4,3),DGAS(4,3),DGBS(4,3),DGCS(4,3),
     &DGDS(4,3),DGES(4,3)

C...The following data lines are coefficients needed in the
C...Drees and Grassie photon structure function parametrizations.
      DATA DGAG/-.207E0,.6158E0,1.074E0,0.E0,.8926E-2,.6594E0,
     &.4766E0,.1975E-1,.03197E0,1.018E0,.2461E0,.2707E-1/
      DATA DGBG/-.1987E0,.6257E0,8.352E0,5.024E0,.5085E-1,.2774E0,
     &-.3906E0,-.3212E0,-.618E-2,.9476E0,-.6094E0,-.1067E-1/
      DATA DGCG/5.119E0,-.2752E0,-6.993E0,2.298E0,-.2313E0,.1382E0,
     &6.542E0,.5162E0,-.1216E0,.9047E0,2.653E0,.2003E-2/
      DATA DGAN/2.285E0,-.1526E-1,1330.E0,4.219E0,-.3711E0,1.061E0,
     &4.758E0,-.1503E-1,15.8E0,-.9464E0,-.5E0,-.2118E0/
      DATA DGBN/6.073E0,-.8132E0,-41.31E0,3.165E0,-.1717E0,.7815E0,
     &1.535E0,.7067E-2,2.742E0,-.7332E0,.7148E0,3.287E0/
      DATA DGCN/-.4202E0,.1778E-1,.9216E0,.18E0,.8766E-1,.2197E-1,
     &.1096E0,.204E0,.2917E-1,.4657E-1,.1785E0,.4811E-1/
      DATA DGDN/-.8083E-1,.6346E0,1.208E0,.203E0,-.8915E0,.2857E0,
     &2.973E0,.1185E0,-.342E-1,.7196E0,.7338E0,.8139E-1/
      DATA DGEN/.5526E-1,1.136E0,.9512E0,.1163E-1,-.1816E0,.5866E0,
     &2.421E0,.4059E0,-.2302E-1,.9229E0,.5873E0,-.79E-4/
      DATA DGAS/16.69E0,-.7916E0,1099.E0,4.428E0,-.1207E0,1.071E0,
     &1.977E0,-.8625E-2,6.734E0,-1.008E0,-.8594E-1,.7625E-1/
      DATA DGBS/.176E0,.4794E-1,1.047E0,.25E-1,25.E0,-1.648E0,
     &-.1563E-1,6.438E0,59.88E0,-2.983E0,4.48E0,.9686E0/
      DATA DGCS/-.208E-1,.3386E-2,4.853E0,.8404E0,-.123E-1,1.162E0,
     &.4824E0,-.11E-1,-.3226E-2,.8432E0,.3616E0,.1383E-2/
      DATA DGDS/-.1685E-1,1.353E0,1.426E0,1.239E0,-.9194E-1,.7912E0,
     &.6397E0,2.327E0,-.3321E-1,.9475E0,-.3198E0,.2132E-1/
      DATA DGES/-.1986E0,1.1E0,1.136E0,-.2779E0,.2015E-1,.9869E0,
     &-.7036E-1,.1694E-1,.1059E0,.6954E0,-.6663E0,.3683E0/

C...Photon structure function from Drees and Grassie.
C...Allowed variable range: 1 GeV^2 < Q^2 < 10000 GeV^2.
      X1=1.-X

C...Evaluate gluon content.
      DGA=DGAG(1,NFE)*T**DGAG(2,NFE)+DGAG(3,NFE)*T**(-DGAG(4,NFE))
      DGB=DGBG(1,NFE)*T**DGBG(2,NFE)+DGBG(3,NFE)*T**(-DGBG(4,NFE))
      DGC=DGCG(1,NFE)*T**DGCG(2,NFE)+DGCG(3,NFE)*T**(-DGCG(4,NFE))
      XPGL=DGA*X**DGB*X1**DGC

C...Evaluate up- and down-type quark content.
      DGA=DGAN(1,NFE)*T**DGAN(2,NFE)+DGAN(3,NFE)*T**(-DGAN(4,NFE))
      DGB=DGBN(1,NFE)*T**DGBN(2,NFE)+DGBN(3,NFE)*T**(-DGBN(4,NFE))
      DGC=DGCN(1,NFE)*T**DGCN(2,NFE)+DGCN(3,NFE)*T**(-DGCN(4,NFE))
      DGD=DGDN(1,NFE)*T**DGDN(2,NFE)+DGDN(3,NFE)*T**(-DGDN(4,NFE))
      DGE=DGEN(1,NFE)*T**DGEN(2,NFE)+DGEN(3,NFE)*T**(-DGEN(4,NFE))
      XPQN=X*(X**2+X1**2)/(DGA-DGB*LOG(X1))+DGC*X**DGD*X1**DGE
      DGA=DGAS(1,NFE)*T**DGAS(2,NFE)+DGAS(3,NFE)*T**(-DGAS(4,NFE))
      DGB=DGBS(1,NFE)*T**DGBS(2,NFE)+DGBS(3,NFE)*T**(-DGBS(4,NFE))
      DGC=DGCS(1,NFE)*T**DGCS(2,NFE)+DGCS(3,NFE)*T**(-DGCS(4,NFE))
      DGD=DGDS(1,NFE)*T**DGDS(2,NFE)+DGDS(3,NFE)*T**(-DGDS(4,NFE))
      DGE=DGES(1,NFE)*T**DGES(2,NFE)+DGES(3,NFE)*T**(-DGES(4,NFE))
      DGF=9.
      IF(NFE.EQ.2) DGF=10.
      IF(NFE.EQ.3) DGF=55./6.
      XPQS=DGF*X*(X**2+X1**2)/(DGA-DGB*LOG(X1))+DGC*X**DGD*X1**DGE
      IF(NFE.LE.1) THEN
        XPQU=(XPQS+9.*XPQN)/6.
        XPQD=(XPQS-4.5*XPQN)/6.
      ELSEIF(NFE.EQ.2) THEN
        XPQU=(XPQS+6.*XPQN)/8.
        XPQD=(XPQS-6.*XPQN)/8.
      ELSE
        XPQU=(XPQS+7.5*XPQN)/10.
        XPQD=(XPQS-5.*XPQN)/10.
      ENDIF

      RETURN
      END


C*********************************************************************

      SUBROUTINE RYSPLI(KF,KFLIN,KFLCH,KFLSP)

C...In case of a hadron remnant which is more complicated than just a
C...quark or a diquark, split it into two (partons or hadron + parton).
      DIMENSION KFL(3)

C...Preliminaries. Parton composition.
      KFA=IABS(KF)
      KFS=ISIGN(1,KF)
      KFL(1)=MOD(KFA/1000,10)
      KFL(2)=MOD(KFA/100,10)
      KFL(3)=MOD(KFA/10,10)
      IF(KFLIN.NE.21.AND.KFLIN.NE.22.AND.KFLIN.NE.23) THEN
        KFLR=KFLIN*KFS
      ELSE
        KFLR=KFLIN
      ENDIF
      KFLCH=0

C...Subdivide lepton.
      IF(KFA.GE.11.AND.KFA.LE.18) THEN
        IF(KFLR.EQ.KFA) THEN
          KFLSP=KFS*22
        ELSEIF(KFLR.EQ.22) THEN
          KFLSP=KFA
        ELSEIF(KFLR.EQ.-24.AND.MOD(KFA,2).EQ.1) THEN
          KFLSP=KFA+1
        ELSEIF(KFLR.EQ.24.AND.MOD(KFA,2).EQ.0) THEN
          KFLSP=KFA-1
        ELSEIF(KFLR.EQ.21) THEN
          KFLSP=KFA
          KFLCH=KFS*21
        ELSE
          KFLSP=KFA
          KFLCH=-KFLR
        ENDIF

C...Subdivide photon.
      ELSEIF(KFA.EQ.22) THEN
        IF(KFLR.NE.21) THEN
          KFLSP=-KFLR
        ELSE
          RAGR=0.75*PYR(0)
          KFLSP=1
          IF(RAGR.GT.0.125) KFLSP=2
          IF(RAGR.GT.0.625) KFLSP=3
          IF(PYR(0).GT.0.5) KFLSP=-KFLSP
          KFLCH=-KFLSP
        ENDIF

C...Subdivide meson.
      ELSEIF(KFL(1).EQ.0) THEN
        KFL(2)=KFL(2)*(-1)**KFL(2)
        KFL(3)=-KFL(3)*(-1)**IABS(KFL(2))
        IF(KFLR.EQ.KFL(2)) THEN
          KFLSP=KFL(3)
        ELSEIF(KFLR.EQ.KFL(3)) THEN
          KFLSP=KFL(2)
        ELSEIF(KFLR.EQ.21.AND.PYR(0).GT.0.5) THEN
          KFLSP=KFL(2)
          KFLCH=KFL(3)
        ELSEIF(KFLR.EQ.21) THEN
          KFLSP=KFL(3)
          KFLCH=KFL(2)
        ELSEIF(KFLR*KFL(2).GT.0) THEN
          CALL LUKFDI(-KFLR,KFL(2),KFDUMP,KFLCH)
          KFLSP=KFL(3)
        ELSE
          CALL LUKFDI(-KFLR,KFL(3),KFDUMP,KFLCH)
          KFLSP=KFL(2)
        ENDIF

C...Subdivide baryon.
      ELSE
        NAGR=0
        DO 100 J=1,3
  100   IF(KFLR.EQ.KFL(J)) NAGR=NAGR+1
        IF(NAGR.GE.1) THEN
          RAGR=0.00001+(NAGR-0.00002)*PYR(0)
          IAGR=0
          DO 110 J=1,3
          IF(KFLR.EQ.KFL(J)) RAGR=RAGR-1.
  110     IF(IAGR.EQ.0.AND.RAGR.LE.0.) IAGR=J
        ELSE
          IAGR=1.00001+2.99998*PYR(0)
        ENDIF
        ID1=1
        IF(IAGR.EQ.1) ID1=2
        IF(IAGR.EQ.1.AND.KFL(3).GT.KFL(2)) ID1=3
        ID2=6-IAGR-ID1
        KSP=3
        IF(MOD(KFA,10).EQ.2.AND.KFL(1).EQ.KFL(2)) THEN
          IF(IAGR.NE.3.AND.PYR(0).GT.0.25) KSP=1
        ELSEIF(MOD(KFA,10).EQ.2.AND.KFL(2).GE.KFL(3)) THEN
          IF(IAGR.NE.1.AND.PYR(0).GT.0.25) KSP=1
        ELSEIF(MOD(KFA,10).EQ.2) THEN
          IF(IAGR.EQ.1) KSP=1
          IF(IAGR.NE.1.AND.PYR(0).GT.0.75) KSP=1
        ENDIF
        KFLSP=1000*KFL(ID1)+100*KFL(ID2)+KSP
        IF(KFLR.EQ.21) THEN
          KFLCH=KFL(IAGR)
        ELSEIF(NAGR.EQ.0.AND.KFLR.GT.0) THEN
          CALL LUKFDI(-KFLR,KFL(IAGR),KFDUMP,KFLCH)
        ELSEIF(NAGR.EQ.0) THEN
          CALL LUKFDI(10000+KFLSP,-KFLR,KFDUMP,KFLCH)
          KFLSP=KFL(IAGR)
        ENDIF
      ENDIF

C...Add on correct sign for result.
      KFLCH=KFLCH*KFS
      KFLSP=KFLSP*KFS

      RETURN
      END

C*********************************************************************

      FUNCTION RYGAMM(X)

C...Gives ordinary Gamma function Gamma(x) for positive, real arguments;
C...see M. Abramowitz, I. A. Stegun: Handbook of Mathematical Functions
C...(Dover, 1965) 6.1.36.
      DIMENSION B(8)
      DATA B/-0.577191652,0.988205891,-0.897056937,0.918206857,
     &-0.756704078,0.482199394,-0.193527818,0.035868343/

      NX=INT(X)
      DX=X-NX

      RYGAMM=1.
      DXP=1.
      DO 100 I=1,8
      DXP=DXP*DX
  100 RYGAMM=RYGAMM+B(I)*DXP
      IF(X.LT.1.) THEN
        RYGAMM=RYGAMM/X
      ELSE
        DO 110 IX=1,NX-1
  110   RYGAMM=(X-IX)*RYGAMM
      ENDIF

      RETURN
      END

C***********************************************************************

      SUBROUTINE RYWAUX(IAUX,EPS,WRE,WIM)

C...Calculates real and imaginary parts of the auxiliary functions W1
C...and W2; see R. K. Ellis, I. Hinchliffe, M. Soldate and J. J. van
C...der Bij, Nucl. Phys. B297 (1988) 221.
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /LUDAT1/

      ASINH(X)=LOG(X+SQRT(X**2+1.))
      ACOSH(X)=LOG(X+SQRT(X**2-1.))

      IF(EPS.LT.0.) THEN
        IF(IAUX.EQ.1) WRE=2.*SQRT(1.-EPS)*ASINH(SQRT(-1./EPS))
        IF(IAUX.EQ.2) WRE=4.*(ASINH(SQRT(-1./EPS)))**2
        WIM=0.
      ELSEIF(EPS.LT.1.) THEN
        IF(IAUX.EQ.1) WRE=2.*SQRT(1.-EPS)*ACOSH(SQRT(1./EPS))
        IF(IAUX.EQ.2) WRE=4.*(ACOSH(SQRT(1./EPS)))**2-PARU(1)**2
        IF(IAUX.EQ.1) WIM=-PARU(1)*SQRT(1.-EPS)
        IF(IAUX.EQ.2) WIM=-4.*PARU(1)*ACOSH(SQRT(1./EPS))
      ELSE
        IF(IAUX.EQ.1) WRE=2.*SQRT(EPS-1.)*ASIN(SQRT(1./EPS))
        IF(IAUX.EQ.2) WRE=-4.*(ASIN(SQRT(1./EPS)))**2
        WIM=0.
      ENDIF

      RETURN
      END

C***********************************************************************

      SUBROUTINE RYI3AU(EPS,RAT,Y3RE,Y3IM)

C...Calculates real and imaginary parts of the auxiliary function I3;
C...see R. K. Ellis, I. Hinchliffe, M. Soldate and J. J. van der Bij,
C...Nucl. Phys. B297 (1988) 221.
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /LUDAT1/

      BE=0.5*(1.+SQRT(1.+RAT*EPS))
      IF(EPS.LT.1.) GA=0.5*(1.+SQRT(1.-EPS))

      IF(EPS.LT.0.) THEN
        IF(ABS(EPS).LT.1.E-4.AND.ABS(RAT*EPS).LT.1.E-4) THEN
          F3RE=RYSPEN(-0.25*EPS/(1.+0.25*(RAT-1.)*EPS),0.,1)-
     &    RYSPEN((1.-0.25*EPS)/(1.+0.25*(RAT-1.)*EPS),0.,1)+
     &    RYSPEN(0.25*(RAT+1.)*EPS/(1.+0.25*RAT*EPS),0.,1)-
     &    RYSPEN((RAT+1.)/RAT,0.,1)+0.5*(LOG(1.+0.25*RAT*EPS)**2-
     &    LOG(0.25*RAT*EPS)**2)+LOG(1.-0.25*EPS)*
     &    LOG((1.+0.25*(RAT-1.)*EPS)/(1.+0.25*RAT*EPS))+
     &    LOG(-0.25*EPS)*LOG(0.25*RAT*EPS/(1.+0.25*(RAT-1.)*EPS))
        ELSEIF(ABS(EPS).LT.1.E-4.AND.ABS(RAT*EPS).GE.1.E-4) THEN
          F3RE=RYSPEN(-0.25*EPS/(BE-0.25*EPS),0.,1)-
     &    RYSPEN((1.-0.25*EPS)/(BE-0.25*EPS),0.,1)+
     &    RYSPEN((BE-1.+0.25*EPS)/BE,0.,1)-
     &    RYSPEN((BE-1.+0.25*EPS)/(BE-1.),0.,1)+
     &    0.5*(LOG(BE)**2-LOG(BE-1.)**2)+
     &    LOG(1.-0.25*EPS)*LOG((BE-0.25*EPS)/BE)+
     &    LOG(-0.25*EPS)*LOG((BE-1.)/(BE-0.25*EPS))
        ELSEIF(ABS(EPS).GE.1.E-4.AND.ABS(RAT*EPS).LT.1.E-4) THEN
          F3RE=RYSPEN((GA-1.)/(GA+0.25*RAT*EPS),0.,1)-
     &    RYSPEN(GA/(GA+0.25*RAT*EPS),0.,1)+
     &    RYSPEN((1.+0.25*RAT*EPS-GA)/(1.+0.25*RAT*EPS),0.,1)-
     &    RYSPEN((1.+0.25*RAT*EPS-GA)/(0.25*RAT*EPS),0.,1)+
     &    0.5*(LOG(1.+0.25*RAT*EPS)**2-LOG(0.25*RAT*EPS)**2)+
     &    LOG(GA)*LOG((GA+0.25*RAT*EPS)/(1.+0.25*RAT*EPS))+
     &    LOG(GA-1.)*LOG(0.25*RAT*EPS/(GA+0.25*RAT*EPS))
        ELSE
          F3RE=RYSPEN((GA-1.)/(GA+BE-1.),0.,1)-
     &    RYSPEN(GA/(GA+BE-1.),0.,1)+RYSPEN((BE-GA)/BE,0.,1)-
     &    RYSPEN((BE-GA)/(BE-1.),0.,1)+0.5*(LOG(BE)**2-LOG(BE-1.)**2)+
     &    LOG(GA)*LOG((GA+BE-1.)/BE)+LOG(GA-1.)*LOG((BE-1.)/(GA+BE-1.))
        ENDIF
        F3IM=0.
      ELSEIF(EPS.LT.1.) THEN
        IF(ABS(EPS).LT.1.E-4.AND.ABS(RAT*EPS).LT.1.E-4) THEN
          F3RE=RYSPEN(-0.25*EPS/(1.+0.25*(RAT-1.)*EPS),0.,1)-
     &    RYSPEN((1.-0.25*EPS)/(1.+0.25*(RAT-1.)*EPS),0.,1)+
     &    RYSPEN((1.-0.25*EPS)/(-0.25*(RAT+1.)*EPS),0.,1)-
     &    RYSPEN(1./(RAT+1.),0.,1)+LOG((1.-0.25*EPS)/(0.25*EPS))*
     &    LOG((1.+0.25*(RAT-1.)*EPS)/(0.25*(RAT+1.)*EPS))
          F3IM=-PARU(1)*LOG((1.+0.25*(RAT-1.)*EPS)/(0.25*(RAT+1.)*EPS))
        ELSEIF(ABS(EPS).LT.1.E-4.AND.ABS(RAT*EPS).GE.1.E-4) THEN
          F3RE=RYSPEN(-0.25*EPS/(BE-0.25*EPS),0.,1)-
     &    RYSPEN((1.-0.25*EPS)/(BE-0.25*EPS),0.,1)+
     &    RYSPEN((1.-0.25*EPS)/(1.-0.25*EPS-BE),0.,1)-
     &    RYSPEN(-0.25*EPS/(1.-0.25*EPS-BE),0.,1)+
     &    LOG((1.-0.25*EPS)/(0.25*EPS))*
     &    LOG((BE-0.25*EPS)/(BE-1.+0.25*EPS))
          F3IM=-PARU(1)*LOG((BE-0.25*EPS)/(BE-1.+0.25*EPS))
        ELSEIF(ABS(EPS).GE.1.E-4.AND.ABS(RAT*EPS).LT.1.E-4) THEN
          F3RE=RYSPEN((GA-1.)/(GA+0.25*RAT*EPS),0.,1)-
     &    RYSPEN(GA/(GA+0.25*RAT*EPS),0.,1)+
     &    RYSPEN(GA/(GA-1.-0.25*RAT*EPS),0.,1)-
     &    RYSPEN((GA-1.)/(GA-1.-0.25*RAT*EPS),0.,1)+
     &    LOG(GA/(1.-GA))*LOG((GA+0.25*RAT*EPS)/(1.+0.25*RAT*EPS-GA))
          F3IM=-PARU(1)*LOG((GA+0.25*RAT*EPS)/(1.+0.25*RAT*EPS-GA))
        ELSE
          F3RE=RYSPEN((GA-1.)/(GA+BE-1.),0.,1)-
     &    RYSPEN(GA/(GA+BE-1.),0.,1)+RYSPEN(GA/(GA-BE),0.,1)-
     &    RYSPEN((GA-1.)/(GA-BE),0.,1)+LOG(GA/(1.-GA))*
     &    LOG((GA+BE-1.)/(BE-GA))
          F3IM=-PARU(1)*LOG((GA+BE-1.)/(BE-GA))
         ENDIF
      ELSE
        RSQ=EPS/(EPS-1.+(2.*BE-1.)**2)
        RCTHE=RSQ*(1.-2.*BE/EPS)
        RSTHE=SQRT(MAX(0.,RSQ-RCTHE**2))
        RCPHI=RSQ*(1.+2.*(BE-1.)/EPS)
        RSPHI=SQRT(MAX(0.,RSQ-RCPHI**2))
        R=SQRT(RSQ)
        THE=ACOS(RCTHE/R)
        PHI=ACOS(RCPHI/R)
        F3RE=RYSPEN(RCTHE,RSTHE,1)+RYSPEN(RCTHE,-RSTHE,1)-
     &  RYSPEN(RCPHI,RSPHI,1)-RYSPEN(RCPHI,-RSPHI,1)+
     &  (PHI-THE)*(PHI+THE-PARU(1))
        F3IM=RYSPEN(RCTHE,RSTHE,2)+RYSPEN(RCTHE,-RSTHE,2)-
     &  RYSPEN(RCPHI,RSPHI,2)-RYSPEN(RCPHI,-RSPHI,2)
      ENDIF

      Y3RE=2./(2.*BE-1.)*F3RE
      Y3IM=2./(2.*BE-1.)*F3IM

      RETURN
      END

C***********************************************************************

      FUNCTION RYSPEN(XREIN,XIMIN,IREIM)

C...Calculates real and imaginary part of Spence function; see
C...G. 't Hooft and M. Veltman, Nucl. Phys. B153 (1979) 365.
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /LUDAT1/
      DIMENSION B(0:14)

      DATA B/
     & 1.000000E+00,        -5.000000E-01,         1.666667E-01,
     & 0.000000E+00,        -3.333333E-02,         0.000000E+00,
     & 2.380952E-02,         0.000000E+00,        -3.333333E-02,
     & 0.000000E+00,         7.575757E-02,         0.000000E+00,
     &-2.531135E-01,         0.000000E+00,         1.166667E+00/

      XRE=XREIN
      XIM=XIMIN
      IF(ABS(1.-XRE).LT.1.E-6.AND.ABS(XIM).LT.1.E-6) THEN
        IF(IREIM.EQ.1) RYSPEN=PARU(1)**2/6.
        IF(IREIM.EQ.2) RYSPEN=0.
        RETURN
      ENDIF

      XMOD=SQRT(XRE**2+XIM**2)
      IF(XMOD.LT.1.E-6) THEN
        IF(IREIM.EQ.1) RYSPEN=0.
        IF(IREIM.EQ.2) RYSPEN=0.
        RETURN
      ENDIF

      XARG=SIGN(ACOS(XRE/XMOD),XIM)
      SP0RE=0.
      SP0IM=0.
      SGN=1.
      IF(XMOD.GT.1.) THEN
        ALGXRE=LOG(XMOD)
        ALGXIM=XARG-SIGN(PARU(1),XARG)
        SP0RE=-PARU(1)**2/6.-(ALGXRE**2-ALGXIM**2)/2.
        SP0IM=-ALGXRE*ALGXIM
        SGN=-1.
        XMOD=1./XMOD
        XARG=-XARG
        XRE=XMOD*COS(XARG)
        XIM=XMOD*SIN(XARG)
      ENDIF
      IF(XRE.GT.0.5) THEN
        ALGXRE=LOG(XMOD)
        ALGXIM=XARG
        XRE=1.-XRE
        XIM=-XIM
        XMOD=SQRT(XRE**2+XIM**2)
        XARG=SIGN(ACOS(XRE/XMOD),XIM)
        ALGYRE=LOG(XMOD)
        ALGYIM=XARG
        SP0RE=SP0RE+SGN*(PARU(1)**2/6.-(ALGXRE*ALGYRE-ALGXIM*ALGYIM))
        SP0IM=SP0IM-SGN*(ALGXRE*ALGYIM+ALGXIM*ALGYRE)
        SGN=-SGN
      ENDIF

      XRE=1.-XRE
      XIM=-XIM
      XMOD=SQRT(XRE**2+XIM**2)
      XARG=SIGN(ACOS(XRE/XMOD),XIM)
      ZRE=-LOG(XMOD)
      ZIM=-XARG

      SPRE=0.
      SPIM=0.
      SAVERE=1.
      SAVEIM=0.
      DO 100 I=0,14
      IF(MAX(ABS(SAVERE),ABS(SAVEIM)).LT.1E-30) GOTO 110
      TERMRE=(SAVERE*ZRE-SAVEIM*ZIM)/FLOAT(I+1)
      TERMIM=(SAVERE*ZIM+SAVEIM*ZRE)/FLOAT(I+1)
      SAVERE=TERMRE
      SAVEIM=TERMIM
      SPRE=SPRE+B(I)*TERMRE
  100 SPIM=SPIM+B(I)*TERMIM

  110 IF(IREIM.EQ.1) RYSPEN=SP0RE+SGN*SPRE
      IF(IREIM.EQ.2) RYSPEN=SP0IM+SGN*SPIM

      RETURN
      END

***********************************************************************

      SUBROUTINE RYTEST(MTEST)
      PARAMETER (KSZJ=4000)

C...Purpose: to provide a simple program (disguised as a subroutine) to
C...run at installation as a check that the program works as intended.
      
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      integer N,NPAD,K
      double precision P,V
      SAVE /PYJETS/

      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/LUDAT3/MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5)
      COMMON/RYSUBS/MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200)
      COMMON/RYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      SAVE /LUDAT1/,/LUDAT2/,/LUDAT3/
      SAVE /RYSUBS/,/RYPARS/

C...Common initial values. Loop over initiating conditions.
      MSTP(122)=MAX(0,MIN(2,MTEST))
      MDCY(LUCOMP(111),1)=0
      NERR=0
      DO 130 IPROC=1,8

C...Reset process type, kinematics cuts, and the flags used.
      MSEL=0
      DO 100 ISUB=1,200
  100 MSUB(ISUB)=0
      CKIN(1)=2.
      CKIN(3)=0.
      MSTP(2)=1
      MSTP(33)=0
      MSTP(81)=1
      MSTP(82)=1
      MSTP(111)=1
      MSTP(131)=0
      MSTP(133)=0
      PARP(131)=0.01

C...Prompt photon production at fixed target.
      IF(IPROC.EQ.1) THEN
        PZSUM=300.
        PESUM=SQRT(PZSUM**2+ULMASS(211)**2)+ULMASS(2212)
        PQSUM=2.
        MSEL=10
        CKIN(3)=5.
        CALL RYINIT('FIXT','pi+','p',PZSUM)

C...QCD processes at ISR energies.
      ELSEIF(IPROC.EQ.2) THEN
        PESUM=63.
        PZSUM=0.
        PQSUM=2.
        MSEL=1
        CKIN(3)=5.
        CALL RYINIT('CMS','p','p',PESUM)

C...W production + multiple interactions at CERN Collider.
      ELSEIF(IPROC.EQ.3) THEN
        PESUM=630.
        PZSUM=0.
        PQSUM=0.
        MSEL=12
        CKIN(1)=20.
        MSTP(82)=4
        MSTP(2)=2
        MSTP(33)=3
        CALL RYINIT('CMS','p','pbar',PESUM)

C...W/Z gauge boson pairs + pileup events at the Tevatron.
      ELSEIF(IPROC.EQ.4) THEN
        PESUM=1800.
        PZSUM=0.
        PQSUM=0.
        MSUB(22)=1
        MSUB(23)=1
        MSUB(25)=1
        CKIN(1)=200.
        MSTP(111)=0
        MSTP(131)=1
        MSTP(133)=2
        PARP(131)=0.04
        CALL RYINIT('CMS','p','pbar',PESUM)

C...Higgs production at LHC.
      ELSEIF(IPROC.EQ.5) THEN
        PESUM=17000.
        PZSUM=0.
        PQSUM=2.
        MSUB(3)=1
        MSUB(102)=1
        MSUB(123)=1
        MSUB(124)=1
        PMAS(25,1)=300.
        CKIN(1)=200.
        MSTP(81)=0
        MSTP(111)=0
        CALL RYINIT('CMS','p','p',PESUM)

C...Z' production at SSC.
      ELSEIF(IPROC.EQ.6) THEN
        PESUM=40000.
        PZSUM=0.
        PQSUM=2.
        MSEL=21
        PMAS(32,1)=600.
        CKIN(1)=400.
        MSTP(81)=0
        MSTP(111)=0
        CALL RYINIT('CMS','p','p',PESUM)

C...W pair production at 1 TeV e+e- collider.
      ELSEIF(IPROC.EQ.7) THEN
        PESUM=1000.
        PZSUM=0.
        PQSUM=0.
        MSUB(25)=1
        MSUB(69)=1
        CALL RYINIT('CMS','e+','e-',PESUM)

C...Deep inelastic scattering at a LEP+LHC ep collider.
      ELSEIF(IPROC.EQ.8) THEN
        P(1,1)=0.
        P(1,2)=0.
        P(1,3)=8000.
        P(2,1)=0.
        P(2,2)=0.
        P(2,3)=-80.
        PESUM=8080.
        PZSUM=7920.
        PQSUM=0.
        MSUB(11)=1
        CKIN(3)=50.
        MSTP(111)=0
        CALL RYINIT('USER','p','e-',PESUM)
      ENDIF

C...Generate 20 events of each required type.
      DO 120 IEV=1,20
      CALL RYEVNT
      PESUMM=PESUM
      IF(IPROC.EQ.4) PESUMM=MSTI(41)*PESUM

C...Check conservation of energy/momentum/flavour.
      MERR=0
      DEVE=ABS(PLU(0,4)-PESUMM)+ABS(PLU(0,3)-PZSUM)
      DEVT=ABS(PLU(0,1))+ABS(PLU(0,2))
      DEVQ=ABS(PLU(0,6)-PQSUM)
      IF(DEVE.GT.2E-3*PESUM.OR.DEVT.GT.MAX(0.01,1E-4*PESUM).OR.
     &DEVQ.GT.0.1) MERR=1
      IF(MERR.NE.0) WRITE(MSTU(11),5000) IPROC,IEV

C...Check that all KF codes are known ones, and that partons/particles
C...satisfy energy-momentum-mass relation.
      DO 110 I=1,N
      IF(K(I,1).GT.20) GOTO 110
      IF(LUCOMP(K(I,2)).EQ.0) THEN
        WRITE(MSTU(11),5100) I
        MERR=MERR+1
      ENDIF
      PD=P(I,4)**2-P(I,1)**2-P(I,2)**2-P(I,3)**2-P(I,5)**2*
     &SIGN(1.,P(I,5))
      IF(ABS(PD).GT.MAX(0.1,0.002*P(I,4)**2,0.002*P(I,5)**2).OR.
     &(P(I,5).GE.0..AND.P(I,4).LT.0.)) THEN
        WRITE(MSTU(11),5200) I
        MERR=MERR+1
      ENDIF
  110 CONTINUE

C...Listing of erroneous events, and first event of each type.
      IF(MERR.GE.1) NERR=NERR+1
      IF(NERR.GE.10) THEN
        WRITE(MSTU(11),5300)
        CALL LULIST(1)
        STOP
      ENDIF
      IF(MTEST.GE.1.AND.(MERR.GE.1.OR.IEV.EQ.1)) THEN
        IF(MERR.GE.1) WRITE(MSTU(11),5400)
        CALL LULIST(1)
      ENDIF
  120 CONTINUE

C...List statistics for each process type.
      IF(MTEST.GE.1) CALL RYSTAT(1)
  130 CONTINUE

C...Summarize result of run.
      IF(NERR.EQ.0) WRITE(MSTU(11),5500)
      IF(NERR.GT.0) WRITE(MSTU(11),5600) NERR
      RETURN

C...Formats for information.
 5000 FORMAT(/5X,'Energy/momentum/flavour nonconservation for process',
     &I2,', event',I4)
 5100 FORMAT(/5X,'Entry no.',I4,' in following event not known code')
 5200 FORMAT(/5X,'Entry no.',I4,' in following event has faulty ',
     &'kinematics')
 5300 FORMAT(/5X,'This is the tenth error experienced! Something is ',
     &'wrong.'/5X,'Execution will be stopped after listing of event.')
 5400 FORMAT(5X,'Faulty event follows:')
 5500 FORMAT(//5X,'End result of run: no errors detected.')
 5600 FORMAT(//5X,'End result of run:',I2,' errors detected.'/
     &5X,'This should not have happened!')
      END

C*********************************************************************

      BLOCK DATA RYDATA

C...Give sensible default values to all status codes and parameters.
      COMMON/RYSUBS/MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200)
      COMMON/RYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/RYINT1/MINT(400),VINT(400)
      COMMON/RYINT2/ISET(200),KFPR(200,2),COEF(200,20),ICOL(40,4,2)
      COMMON/RYINT3/XSFX(2,-40:40),ISIG(1000,3),SIGH(1000)
      COMMON/RYINT4/WIDP(21:40,0:40),WIDE(21:40,0:40),WIDS(21:40,3)
      COMMON/RYINT5/NGEN(0:200,3),XSEC(0:200,3)
      COMMON/RYINT6/PROC(0:200)
      CHARACTER PROC*28
      SAVE /RYSUBS/,/RYPARS/,/RYINT1/,/RYINT2/,/RYINT3/,/RYINT4/,
     &/RYINT5/,/RYINT6/

C...Default values for allowed processes and kinematics constraints.
      DATA MSEL/1/
      DATA MSUB/200*0/
      DATA ((KFIN(I,J),J=-40,40),I=1,2)/40*1,0,80*1,0,40*1/
      DATA CKIN/
     &   2.0, -1.0,  0.0, -1.0,  1.0,  1.0, -10.,  10., -10.,  10.,
     1  -10.,  10., -10.,  10., -10.,  10., -1.0,  1.0, -1.0,  1.0,
     2   0.0,  1.0,  0.0,  1.0, -1.0,  1.0, -1.0,  1.0,   0.,   0.,
     3   2.0, -1.0,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
     4  12.0, -1.0, 12.0, -1.0, 12.0, -1.0, 12.0, -1.0,   0.,   0.,
     5   0.0, -1.0,  0.0, -1.0,  0.0, -1.0,   0.,   0.,   0.,   0.,
     6   140*0./

C...Default values for main switches and parameters. Reset information.
      DATA (MSTP(I),I=1,100)/
     &     3,    1,    2,    0,    0,    0,    0,    0,    0,    0,
     1     0,    0,    1,    0,    0,    0,    0,    0,    0,    0,
     2     1,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     3     1,    2,    0,    0,    0,    2,    1,    5,    0,    0,
     4     1,    1,    3,    7,    2,    1,    1,    0,    0,    0,
     5     1,    1,   20,    6,    0,    0,    0,    0,    0,    0,
     6     1,    2,    2,    2,    1,    0,    0,    0,    0,    0,
     7     1,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     8     1,    1,  100,    0,    0,    0,    0,    0,    0,    0,
     9     1,    4,    0,    0,    0,    0,    0,    0,    0,    0/
      DATA (MSTP(I),I=101,200)/
     &     1,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     1     1,    1,    1,    0,    0,    0,    0,    0,    0,    0,
     2     0,    1,    2,    1,    1,   20,    1,    0,   10,    0,
     3     0,    4,    0,    1,    0,    0,    0,    0,    0,    0,
     4     0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     5     0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     6     0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     7     0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     8     5,    5, 1991,   03,   08,    0,    0,    0,    0,    0,
     9     0,    0,    0,    0,    0,    0,    0,    0,    0,    0/
      DATA (PARP(I),I=1,100)/
     &  0.25,  10.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
     1   0.0,-0.67, 25.0, 0.01,   0.,   0.,   0.,   0.,   0.,   0.,
     2    0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
     3   1.5,  2.0, 0.075,  0.,  0.2,   0.,  2.0,   0.,   0.,   0.,
     4  0.01,  2.0, 0.10, 1000., 2054., 0.,   0.,   0.,   0.,   0.,
     5   1.0, 2.26, 1.E4, 1.E-4,  0.,   0.,   0.,   0.,   0.,   0.,
     6  0.25,  1.0, 0.25,  1.0,  2.0, 1.E-3, 4.0,   0.,   0.,   0.,
     7   4.0,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
     8  1.45, 1.70,  0.5,  0.2, 0.33, 0.66,  0.7,  0.5,   0.,   0.,
     9  0.44, 0.44,  2.0,  1.0,   0.,  3.0,  1.0, 0.75,   0.,   0./
      DATA (PARP(I),I=101,200)/
     & -0.02,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
     1   2.0,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
     2   1.0,  0.4,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
     3  0.01,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
     4    0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
     5    0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
     6    0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
     7    0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
     8    0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
     9    0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0./
      DATA MSTI/200*0/
      DATA PARI/200*0./
      DATA MINT/400*0/
      DATA VINT/400*0./

C...Constants for the generation of the various processes.
      DATA (ISET(I),I=1,100)/
     &    1,    1,    1,   -1,    3,   -1,   -1,    3,   -2,   -2,
     1    2,    2,    2,    2,    2,    2,   -1,    2,    2,    2,
     2   -1,    2,    2,    2,    2,    2,   -1,    2,    2,    2,
     3    2,   -1,    2,    2,    2,    2,   -1,   -1,   -1,   -1,
     4   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
     5   -1,   -1,    2,    2,   -1,   -1,   -1,    2,   -1,   -1,
     6   -1,   -1,   -1,   -1,   -1,   -1,   -1,    2,    2,    2,
     7    4,    4,    4,   -1,   -1,    4,    4,   -1,   -1,   -2,
     8    2,    2,    2,    2,    2,   -2,   -2,   -2,   -2,   -2,
     9    0,    0,    0,   -1,    0,    9,   -2,   -2,   -2,   -2/
      DATA (ISET(I),I=101,200)/
     &   -1,    1,    1,   -2,   -2,   -2,   -2,   -2,   -2,   -2,
     1    2,    2,    2,    2,    2,   -1,   -1,   -1,   -2,   -2,
     2   -1,   -1,    5,    5,   -2,   -2,   -2,   -2,   -2,   -2,
     3    6,   -2,   -2,   -2,   -2,   -2,   -2,   -2,   -2,   -2,
     4    1,    1,    1,    1,    1,   -2,   -2,   -2,   -2,   -2,
     5    1,    1,    1,   -2,   -2,    1,    1,    1,   -2,   -2,
     6    2,    2,    2,    2,   -2,   -2,   -2,   -2,   -2,   -2,
     7    2,    2,    5,    5,   -2,    2,    2,    5,    5,   -2,
     8   -2,   -2,   -2,   -2,   -2,   -2,   -2,   -2,   -2,   -2,
     9   -2,   -2,   -2,   -2,   -2,   -2,   -2,   -2,   -2,   -2/
      DATA ((KFPR(I,J),J=1,2),I=1,50)/
     &   23,    0,   24,    0,   25,    0,   24,    0,   25,    0,
     &   24,    0,   23,    0,   25,    0,    0,    0,    0,    0,
     1    0,    0,    0,    0,   21,   21,   21,   22,   21,   23,
     1   21,   24,   21,   25,   22,   22,   22,   23,   22,   24,
     2   22,   25,   23,   23,   23,   24,   23,   25,   24,   24,
     2   24,   25,   25,   25,    0,   21,    0,   22,    0,   23,
     3    0,   24,    0,   25,    0,   21,    0,   22,    0,   23,
     3    0,   24,    0,   25,    0,   21,    0,   22,    0,   23,
     4    0,   24,    0,   25,    0,   21,    0,   22,    0,   23,
     4    0,   24,    0,   25,    0,   21,    0,   22,    0,   23/
      DATA ((KFPR(I,J),J=1,2),I=51,100)/
     5    0,   24,    0,   25,    0,    0,    0,    0,    0,    0,
     5    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     6    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     6    0,    0,    0,    0,   21,   21,   24,   24,   23,   24,
     7   23,   23,   24,   24,   23,   24,   23,   25,   22,   22,
     7   23,   23,   24,   24,   24,   25,   25,   25,    0,    0,
     8    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     8    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     9    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     9    0,    0,    0,    0,    0,    0,    0,    0,    0,    0/
      DATA ((KFPR(I,J),J=1,2),I=101,150)/
     &   23,    0,   25,    0,   25,    0,    0,    0,    0,    0,
     &    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     1   21,   25,    0,   25,   21,   25,   22,   22,   21,   22,
     1   22,   23,   23,   23,   24,   24,    0,    0,    0,    0,
     2   25,    0,   25,    0,   25,    0,   25,    0,    0,    0,
     2    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     3   23,    5,    0,    0,    0,    0,    0,    0,    0,    0,
     3    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     4   32,    0,   34,    0,   37,    0,   40,    0,   39,    0,
     4    0,    0,    0,    0,    0,    0,    0,    0,    0,    0/
      DATA ((KFPR(I,J),J=1,2),I=151,200)/
     5   35,    0,   35,    0,   35,    0,    0,    0,    0,    0,
     5   36,    0,   36,    0,   36,    0,    0,    0,    0,    0,
     6    6,   37,   39,    0,   39,   39,   39,   39,    0,    0,
     6    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     7   23,   35,   24,   35,   35,    0,   35,    0,    0,    0,
     7   23,   36,   24,   36,   36,    0,   36,    0,    0,    0,
     8    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     8    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     9    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     9    0,    0,    0,    0,    0,    0,    0,    0,    0,    0/
      DATA COEF/4000*0./
      DATA (((ICOL(I,J,K),K=1,2),J=1,4),I=1,40)/
     1 4,0,3,0,2,0,1,0,3,0,4,0,1,0,2,0,2,0,0,1,4,0,0,3,3,0,0,4,1,0,0,2,
     2 3,0,0,4,1,4,3,2,4,0,0,3,4,2,1,3,2,0,4,1,4,0,2,3,4,0,3,4,2,0,1,2,
     3 3,2,1,0,1,4,3,0,4,3,3,0,2,1,1,0,3,2,1,4,1,0,0,2,2,4,3,1,2,0,0,1,
     4 3,2,1,4,1,4,3,2,4,2,1,3,4,2,1,3,3,4,4,3,1,2,2,1,2,0,3,1,2,0,0,0,
     5 4,2,1,0,0,0,1,0,3,0,0,3,1,2,0,0,4,0,0,4,0,0,1,2,2,0,0,1,4,4,3,3,
     6 2,2,1,1,4,4,3,3,3,3,4,4,1,1,2,2,3,2,1,3,1,2,0,0,4,2,1,4,0,0,1,2,
     7 4,0,0,0,4,0,1,3,0,0,3,0,2,4,3,0,3,4,0,0,1,0,0,1,0,0,3,4,2,0,0,2,
     8 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     9 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     & 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/

C...Character constants: name of processes.
      DATA PROC(0)/                    'All included subprocesses   '/
      DATA (PROC(I),I=1,20)/
     1'f + f~ -> gamma*/Z0         ',  'f + f~'' -> W+/-             ',
     2'f + f~ -> H0                ',  'gamma + W+/- -> W+/-        ',
     3'Z0 + Z0 -> H0               ',  'Z0 + W+/- -> W+/-           ',
     4'                            ',  'W+ + W- -> H0               ',
     5'                            ',  '                            ',
     6'f + f'' -> f + f''            ','f + f~ -> f'' + f~''          ',
     7'f + f~ -> g + g             ',  'f + f~ -> g + gamma         ',
     8'f + f~ -> g + Z0            ',  'f + f~'' -> g + W+/-         ',
     9'f + f~ -> g + H0            ',  'f + f~ -> gamma + gamma     ',
     &'f + f~ -> gamma + Z0        ',  'f + f~'' -> gamma + W+/-     '/
      DATA (PROC(I),I=21,40)/
     1'f + f~ -> gamma + H0        ',  'f + f~ -> Z0 + Z0           ',
     2'f + f~'' -> Z0 + W+/-        ', 'f + f~ -> Z0 + H0           ',
     3'f + f~ -> W+ + W-           ',  'f + f~'' -> W+/- + H0        ',
     4'f + f~ -> H0 + H0           ',  'f + g -> f + g              ',
     5'f + g -> f + gamma          ',  'f + g -> f + Z0             ',
     6'f + g -> f'' + W+/-          ', 'f + g -> f + H0             ',
     7'f + gamma -> f + g          ',  'f + gamma -> f + gamma      ',
     8'f + gamma -> f + Z0         ',  'f + gamma -> f'' + W+/-      ',
     9'f + gamma -> f + H0         ',  'f + Z0 -> f + g             ',
     &'f + Z0 -> f + gamma         ',  'f + Z0 -> f + Z0            '/
      DATA (PROC(I),I=41,60)/
     1'f + Z0 -> f'' + W+/-         ', 'f + Z0 -> f + H0            ',
     2'f + W+/- -> f'' + g          ', 'f + W+/- -> f'' + gamma      ',
     3'f + W+/- -> f'' + Z0         ', 'f + W+/- -> f'' + W+/-       ',
     4'f + W+/- -> f'' + H0         ', 'f + H0 -> f + g             ',
     5'f + H0 -> f + gamma         ',  'f + H0 -> f + Z0            ',
     6'f + H0 -> f'' + W+/-         ', 'f + H0 -> f + H0            ',
     7'g + g -> f + f~             ',  'g + gamma -> f + f~         ',
     8'g + Z0 -> f + f~            ',  'g + W+/- -> f + f~''         ',
     9'g + H0 -> f + f~            ',  'gamma + gamma -> f + f~     ',
     &'gamma + Z0 -> f + f~        ',  'gamma + W+/- -> f + f~''     '/
      DATA (PROC(I),I=61,80)/
     1'gamma + H0 -> f + f~        ',  'Z0 + Z0 -> f + f~           ',
     2'Z0 + W+/- -> f + f~''        ', 'Z0 + H0 -> f + f~           ',
     3'W+ + W- -> f + f~           ',  'W+/- + H0 -> f + f~''        ',
     4'H0 + H0 -> f + f~           ',  'g + g -> g + g              ',
     5'gamma + gamma -> W+ + W-    ',  'gamma + W+/- -> Z0 + W+/-   ',
     6'Z0 + Z0 -> Z0 + Z0          ',  'Z0 + Z0 -> W+ + W-          ',
     7'Z0 + W+/- -> Z0 + W+/-      ',  'Z0 + Z0 -> Z0 + H0          ',
     8'W+ + W- -> gamma + gamma    ',  'W+ + W- -> Z0 + Z0          ',
     9'W+/- + W+/- -> W+/- + W+/-  ',  'W+/- + H0 -> W+/- + H0      ',
     &'H0 + H0 -> H0 + H0          ',  '                            '/
      DATA (PROC(I),I=81,100)/
     1'q + q~ -> Q + Q~, massive   ',  'g + g -> Q + Q~, massive    ',
     2'f + q -> f'' + Q, massive    ', 'g + gamma -> Q + Q~, massive',
     3'gamma + gamma -> F + F~, mas',  '                            ',
     4'                            ',  '                            ',
     5'                            ',  '                            ',
     6'Elastic scattering          ',  'Single diffractive          ',
     7'Double diffractive          ',  'Central diffractive         ',
     8'Low-pT scattering           ',  'Semihard QCD 2 -> 2         ',
     9'                            ',  '                            ',
     &'                            ',  '                            '/
      DATA (PROC(I),I=101,120)/
     1'g + g -> gamma*/Z0          ',  'g + g -> H0                 ',
     2'gamma + gamma -> H0         ',  '                            ',
     3'                            ',  '                            ',
     4'                            ',  '                            ',
     5'                            ',  '                            ',
     6'f + f~ -> g + H0            ',  'q + g -> q + H0             ',
     7'g + g -> g + H0             ',  'g + g -> gamma + gamma      ',
     8'g + g -> g + gamma          ',  'g + g -> gamma + Z0         ',
     9'g + g -> Z0 + Z0            ',  'g + g -> W+ + W-            ',
     &'                            ',  '                            '/
      DATA (PROC(I),I=121,140)/
     1'g + g -> f + f~ + H0        ',  'gamma + gamma -> f + f~ + H0',
     2'f + f'' -> f + f'' + H0       ',
     2'f + f'' -> f" + f"'' + H0     ',
     3'                            ',  '                            ',
     4'                            ',  '                            ',
     5'                            ',  '                            ',
     6'g + g -> Z0 + q + q~        ',  '                            ',
     7'                            ',  '                            ',
     8'                            ',  '                            ',
     9'                            ',  '                            ',
     &'                            ',  '                            '/
      DATA (PROC(I),I=141,160)/
     1'f + f~ -> gamma*/Z0/Z''0     ', 'f + f~'' -> W''+/-            ',
     2'f + f~'' -> H+/-             ', 'f + f~'' -> R                ',
     3'q + l -> LQ                 ',  '                            ',
     4'                            ',  '                            ',
     5'                            ',  '                            ',
     6'f + f~ -> H''0               ', 'g + g -> H''0                ',
     7'gamma + gamma -> H''0        ', '                            ',
     8'                            ',  'f + f~ -> A0                ',
     9'g + g -> A0                 ',  'gamma + gamma -> A0         ',
     &'                            ',  '                            '/
      DATA (PROC(I),I=161,180)/
     1'f + g -> f'' + H+/-          ', 'q + g -> LQ + l~            ',
     2'g + g -> LQ + LQ~           ',  'q + q~ -> LQ + LQ~          ',
     3'                            ',  '                            ',
     4'                            ',  '                            ',
     5'                            ',  '                            ',
     6'f + f~ -> Z0 + H''0          ', 'f + f~'' -> W+/- + H''0       ',
     7'f + f'' -> f + f'' + H''0      ',
     7'f + f'' -> f" + f"'' + H''0    ',
     8'                            ',  'f + f~ -> Z0 + A0           ',
     9'f + f~'' -> W+/- + A0        ',
     9'f + f'' -> f + f'' + A0       ',
     &'f + f'' -> f" + f"'' + A0     ',
     &'                            '/
      DATA (PROC(I),I=181,200)/     20*'                            '/

      END


C*********************************************************************

C...The following routines have been written by Ronald Kleiss,
C...to evaluate the matrix element for g + g -> Z + q + qbar,
C...with massive quarks (e.g. q = b).
C...They have been modified, so that all routines and commonblocks
C...have names beginning with RK, and so that some unnecessary
C...initialization information is not printed. Further, COMPLEX*16
C...has been changed to COMPLEX and REAL*8 to DOUBLE PRECISION
C...(in a few cases to REAL), so as to make the program better
C...transportable.

      SUBROUTINE RKBBV(AK1,AK2,AP1,AP2,ALEP1,ALEP2,IMC,RESULT)
* THE CROSS SECTION FOR
* G(K1) + G(K2) ---> Z(QV) + B(P1) + B_BAR(P2)
*                     |
*                     +---> L(LEP1) + LEP_BAR(LEP2)
* THE B QUARKS HAVE TO BE ON-SHELL, THE LEPTONS MASSLESS
* THE OPTION IMC=0 PERFORMS THE STANDARD SPIN SUM
* THE OPTION IMC=1 PERFORMS THE CALCULATION FOR 'NMC' RANDOMLY
* CHOSEN HELICITY STATES WHICH IMPROVES THE
* SPEED BY A FACTOR 32/NMC
      IMPLICIT NONE
      SAVE

      REAL AK1(0:3),AK2(0:3),AP1(0:3),AP2(0:3),ALEP1(0:3),ALEP2(0:3)
      DOUBLE PRECISION K1(0:4),K2(0:4),P1(0:4),P2(0:4),LEP1(0:4),
     &LEP2(0:4)
      REAL RMQ,RMV,RGV,GSTR,VB,AB,VL,AL
      INTEGER INIT
      INTEGER J1,J2,J3,J4,J5
      INTEGER K,IMC,KLOW,KUPP,NMC,OLDIMC
      DOUBLE PRECISION RKRAND,RKDOT,MULT,RMB
      INTEGER CHKGL1,CHKGL2
      DOUBLE PRECISION QV(0:4),R1(0:4),R2(0:4),Q1(0:4),Q2(0:4)
      DOUBLE PRECISION PP2(0:4)
      DOUBLE PRECISION CROSS
      INTEGER LG1,LG2,LV,L1,L2,HELIX,HELI
      COMPLEX ZFACV,ZFAC1,ZFAC2
      DOUBLE PRECISION ZFACS,ZFACB,ZFACBB,ZFACL
      COMPLEX RKZPR,RKZSF
      COMPLEX ZFAC
      DOUBLE PRECISION VPA,VMA
      DOUBLE PRECISION RR1(0:4),RR2(0:4)
      DOUBLE PRECISION ZD12V,ZD21V,ZD1V2,ZD2V1,ZDV12,ZDV21
      COMPLEX RKZF,ZN12V,ZN21V,ZN1V2,ZN2V1,ZNV12,ZNV21
      COMPLEX ZDIA1,ZDIA2,ZDIA3,ZDIA4,ZDIA5,ZDIA6,ZDIA7,ZDIA8
      COMPLEX ZC12V,ZC21V,ZCV12,ZCV21
      DOUBLE PRECISION S,ZD11,ZD22
      COMPLEX ZABEL,ZNABEL,ZNABEM
      REAL RESULT
      DOUBLE PRECISION AZABEL,CABEL,ANABEL,CNABEL,THIS1
      COMPLEX ANSS(-1:1,1:4,-1:1,1:4)
      INTEGER DONS(-1:1,1:4,-1:1,1:4)
      COMPLEX ANSF(-1:1,1:4,1:8,-1:1,1:4)
      INTEGER DONF(-1:1,1:4,1:8,-1:1,1:4)

      PARAMETER(CHKGL1=0,CHKGL2=0,NMC=1)

      COMMON / RKZSCO / ANSS,DONS
      COMMON / RKZFCO / ANSF,DONF
      COMMON / RKBBVC / RMQ,RMV,RGV,VB,AB,VL,AL
      DATA INIT/0/

* CHECK ON EITHER FIRST CALL OR CHANGE IN IMC
      IF(INIT.EQ.0.OR.IMC.NE.OLDIMC) THEN
        OLDIMC=IMC
        INIT=1
* REPRODUCE INPUT DATA
C       WRITE(6,*) ' ----------------------------------------'
C       WRITE(6,*) ' BBV: G G ---> B B_BAR Z, Z ---> L L_BAR'
C       WRITE(6,*) ' B QUARK MASS      = ',RMB,' GEV'
C       WRITE(6,*) ' BOSON MASS        = ',RMV,' GEV'
C       WRITE(6,*) ' BOSON WIDTH       = ',RGV,' GEV'
C       WRITE(6,*) ' B VECTOR C.       = ',VB
C       WRITE(6,*) ' B AXIAL C.        = ',AB
C       WRITE(6,*) ' LEPTON VECTOR C.  = ',VL
C       WRITE(6,*) ' LEPTON AXIAL C.   = ',AL
        RMB=RMQ
* ADJUST STRONG COUPLING SO AS TO GIVE EFFECTIVELY ALPHA_S=1
        GSTR=4D0*DSQRT(DATAN(1D0))
C       WRITE(6,*) ' QCD COUPLING      = ',GSTR
* SEE WETHER GAUGE CHECKS ARE REQUIRED
        IF(CHKGL1.EQ.1) THEN
          WRITE(6,*) ' GAUGE CHECK ON GLUON 1'
        ENDIF
        IF(CHKGL2.EQ.1) THEN
          WRITE(6,*) ' GAUGE CHECK ON GLUON 2'
        ENDIF
* SEE WETHER HELICITY MONTE CARLO IS REQUIRED
        IF(IMC.EQ.0) THEN
          KLOW=1
          KUPP=32
          MULT=1D0
          WRITE(6,*) ' SUM OVER HELICITIES SELECTED'
        ELSEIF(IMC.EQ.1) THEN
          KLOW=1
          KUPP=NMC
          MULT=32D0/(1D0*NMC)
C         WRITE(6,*) ' MONTE CARLO OVER HELICITES SELECTED'
C         WRITE(6,*) ' WITH ',NMC,' HELICITY TRIALS'
C         WRITE(6,*) ' RESULT THEN MULTIPLIED BY ',MULT
        ELSE
          WRITE(6,*) ' ERROR: WRONG OPTION IMC=',IMC
        ENDIF
C       WRITE(6,*) ' THE RESULT IS BASED ON ALPHA_S=1,',
C    .  ' MUST BE MULTIPLIED BY ALPHA_S**2'
C       WRITE(6,*) ' ----------------------------------------'
C       WRITE(6,800)'NO.','LG1','LG2','LV','L1','L2','AMP**2'
C 800   FORMAT(' ',6A4,A10)
      ENDIF

* INITIALIZE THE ARRAYS ANSS,DONS
      DO 130 J1=-1,1,2
        DO 120 J2=1,4
          DO 110 J3=-1,1,2
            DO 100 J4=1,4
              ANSS(J1,J2,J3,J4)=(0.,0.)
              DONS(J1,J2,J3,J4)=0
  100       CONTINUE
  110     CONTINUE
  120   CONTINUE
  130 CONTINUE

* INITIALIZE THE ARRAYS ANSF,DONF
      DO 180 J1=-1,1,2
        DO 170 J2=1,4
          DO 160 J3=1,8
            DO 150 J4=-1,1,2
              DO 140 J5=1,4
                 ANSF(J1,J2,J3,J4,J5)=(0.,0.)
                 DONF(J1,J2,J3,J4,J5)=0
  140         CONTINUE
  150       CONTINUE
  160     CONTINUE
  170   CONTINUE
  180 CONTINUE

* EQUATE THE (0:4) INTERNAL MOMENTA TO THE (0:3) ARGUMENTS MOMENTA
      DO 190 K=0,3
        K1(K)=AK1(K)
        K2(K)=AK2(K)
        P1(K)=AP1(K)
        P2(K)=AP2(K)
        LEP1(K)=ALEP1(K)
        LEP2(K)=ALEP2(K)
  190 CONTINUE

* ASSIGN LABELS TO THE MOMENTA FOR RECOGNITION
* THE MOMENTA K1,K2,LEP1,LEP2 (AND R1,R2) CAN OCCUR AS THE MASSLESS
* MOMENTA IN ARGUMENTS NO.2 AND 6 IN ZF, AND NO.2 AND 4 IN RKZSF
* R1,R2 AND Q1,Q2 ARE SOME OF THESE, AND CAN ALSO OCCUR
* AS ARGUMENTS NO.2 AND 6 IN ZF AND NO.2 AND 4 IN RKZSF
        K1(4)=1D0
        K2(4)=2D0
        LEP1(4)=3D0
        LEP2(4)=4D0
* THE OTHER MOMENTA P1,P2 AND THE VARIOUS RR1,RR2 CAN OCCUR ONLY
* AS ARGUMENT NO.3 IN ZF
        P1(4)=1D0
        P2(4)=2D0

* THE TOTAL BOSON MOMENTUM
* NO NEED TO ASSIGN 4TH COMPONENT LABEL SINCE IT IS NOT USED
      DO 200 K=0,3
        QV(K)=LEP1(K)+LEP2(K)
  200 CONTINUE


* DEFINE THE AUXILIARY VECTORS: THE RESULT SHOULD BE THE SAME
* FOR EVERY NON-SINGULAR CHOICE OF THE AUXILIARY VECTORS
* SINGULAR CHOICES ARE R1=K1 OR R2=K2
* THESE ARE OBTAINED BY PUTTING CHKGL1=1 OR CHKGL2=1

* AUXILIARY VECTOR FOR GLUON 1
* NEED TO ASSIGN ALSO 4TH COMPONENT LABELS HERE!
      IF(CHKGL1.EQ.1) THEN
        DO 210 K=0,4
          R1(K)=K1(K)
  210   CONTINUE
      ELSE
        DO 220 K=0,4
          R1(K)=K2(K)
  220   CONTINUE
      ENDIF

* AUXILIARY VECTOR FOR GLUON 2
      IF(CHKGL2.EQ.1) THEN
        DO 230 K=0,4
          R2(K)=K2(K)
  230   CONTINUE
      ELSE
        DO 240 K=0,4
          R2(K)=K1(K)
  240   CONTINUE
      ENDIF

* AUXILIARY VECTOR FOR THE B QUARK
      DO 250 K=0,4
        Q1(K)=LEP1(K)
  250 CONTINUE

* AUXILIARY VECTOR FOR THE B_BAR QUARK
      DO 260 K=0,4
        Q2(K)=LEP2(K)
  260 CONTINUE

* INITIALIZE THE CROSS SECTION TO ZERO
      CROSS=0D0

* SINCE P2 CORRESPONDS TO AN ANTIFERMION WE HAVE TO
* CHANGE ITS SIGN MOMENTARILY: PUT THE OLD RESULT IN PP2(0:3)
* BU MAKE SURE TO KEEP THE LABEL POSITIVE!
      DO 270 K=0,3
        PP2(K)=P2(K)
        P2(K)=-P2(K)
  270 CONTINUE

* COMPUTE OVERALL FACTORS: FOR EVERY SLASHED POLARIZATION THERE
* APPEARS A FACTOR OF 2 IN ADDITION TO THE NORMALIZATION
* FOLLOWING FROM THE CHISHOLM IDENTITY
* IN PRINCIPLE THE OVERALL FACTORS ARE DIFFERENT FOR EACH DIFFERENT
* HELICITY COMPBINATION BUT IN THIS CASE WE ARE ONLY INTERESTED IN
* THEIR ABSOLUTE VALUE (NO TRANSVERSE GLUON POLARIZATION ETC.)
* SO WE CAN TAKE THIS OUT OF THE LOOP, EXCEPT FOR THE NONTRIVIAL
* HELICITY DEPENDENCE IN 'ZFACV'

* OVERALL FACTOR FOR THE BOSON CURRENT, WITH BREIT-WIGNER
      ZFACV=2./CMPLX(SNGL(RKDOT(QV,QV))-RMV**2,RMV*RGV)

* OVERALL FACTOR FOR GLUON 1
      IF(CHKGL1.EQ.1) THEN
        ZFAC1=(1.,0.)
      ELSE
* ORIGINAL FORM: ZFAC1=2D0*LG1/(DSQRT(2D0)*RKZPR(-LG1,K1,R1))
        ZFAC1=DSQRT(2D0)/RKZSF(1,K1,-1,R1)
      ENDIF

* OVERALL FACTOR FOR GLUON 2
      IF(CHKGL2.EQ.1) THEN
        ZFAC2=1D0
      ELSE
* ORIGINAL FORM: ZFAC2=2D0*LG2/(DSQRT(2D0)*RKZPR(-LG2,K2,R2))
        ZFAC2=DSQRT(2D0)/RKZSF(1,K2,-1,R2)
      ENDIF

* OVERALL FACTOR FOR QCD COUPLINGS
      ZFACS=GSTR**2

* OVERALL FACTOR FOR THE B QUARK
      ZFACB=1/DSQRT(2D0*RKDOT(P1,Q1))

* OVERALL FACTOR FOR THE B_BAR QUARK
      ZFACBB=1D0/DSQRT(2D0*RKDOT(PP2,Q2))

* FINAL OVERALL FACTOR
      ZFAC=ZFACV*ZFAC1*ZFAC2*ZFACS*ZFACB*ZFACBB

* DO A BIG LOOP OVER ALL HELICITIES OR A RANDOM CHOICE OF HELICITIES
* NB: FUNNY INDENTATION HERE!
* ALSO INITIALIZE COUNTERS FOR RKZSF AND ZF

      DO 360 HELIX=KLOW,KUPP
      IF(IMC.EQ.0) THEN
        CALL RKHLPK(HELIX,LG1,LG2,LV,L1,L2)
      ELSE
        HELI=IDINT(32D0*RKRAND(HELIX))+1
        CALL RKHLPK(HELI,LG1,LG2,LV,L1,L2)
      ENDIF

* DETERMINE THE 'LEFT-' AND 'RIGHT-'HANDED COUPLINGS OF THE B TO THE Z
      VPA=VB+LV*AB
      VMA=VB-LV*AB
* AND THE LEPTON HELICITY FACTOR
      ZFACL=(VL-LV*AL)

* FIRST PART OF THE RESULT: THE ABELIAN TERMS
* COMPUTE THE NUMERATORS (ZN...) USING THE ZF FUNCTION
* AND THE DENOMINATORS (ZD...) THE STANDARD WAY
* THE INTERNAL FERMION MOMENTA ARE DIFFERENT IN EACH DIAGRAM
* AND ARE DENOTED BY RR1 AND RR2
* THE 4TH COMPONENT LABELS ARE NONTRIVIAL HERE: HAVING ALREADY
* P1(4)=1 AND P2(4)=2 WE ALSO DEFINE
* (P1-K1)(4)=3,
* (P1-K1-K2)(4)=(P1-K2-K1)(4)=4
* (P1-K2)(4)=5
* (P1-K1+QV)(4)=6
* (P1-K2+QV)(5)=7
* (P1+QV)(4)=8
* SO THAT IN THE VARIOUS DIAGRAMS WE HAVE
* IN ZN12V: RR1(4)=3, RR2(4)=4
* IN ZN21V: RR1(4)=5, RR2(4)=4
* IN ZN1V2: RR1(4)=3, RR2(4)=6
* IN ZN2V1: RR1(4)=5, RR2(4)=7
* IN ZNV12: RR1(4)=8, RR2(4)=6
* IN ZNV21: RR1(4)=8, RR2(4)=7

      DO 280 K=0,3
        RR1(K)=P1(K)-K1(K)
        RR2(K)=RR1(K)-K2(K)
  280 CONTINUE
      RR1(4)=3D0
      RR2(4)=4D0
      ZD12V=(RKDOT(RR1,RR1)-RMB**2)*(RKDOT(RR2,RR2)-RMB**2)
      ZN12V =
     . + RKZF(L1,Q1,P1,RMB,LG1,R1)     *RKZF(LG1,K1,RR1,RMB,LG2,R2)
     .  *RKZF(LG2,K2,RR2,RMB,LV,LEP2)  *RKZF(LV,LEP1,P2,RMB,L2,Q2)*VPA
     . + RKZF(L1,Q1,P1,RMB,LG1,R1)     *RKZF(LG1,K1,RR1,RMB,LG2,R2)
     .  *RKZF(LG2,K2,RR2,RMB,-LV,LEP1) *RKZF(-LV,LEP2,P2,RMB,L2,Q2)*VMA
     . + RKZF(L1,Q1,P1,RMB,LG1,R1)     *RKZF(LG1,K1,RR1,RMB,-LG2,K2)
     .  *RKZF(-LG2,R2,RR2,RMB,LV,LEP2) *RKZF(LV,LEP1,P2,RMB,L2,Q2)*VPA
     . + RKZF(L1,Q1,P1,RMB,LG1,R1)     *RKZF(LG1,K1,RR1,RMB,-LG2,K2)
     .  *RKZF(-LG2,R2,RR2,RMB,-LV,LEP1)*RKZF(-LV,LEP2,P2,RMB,L2,Q2)*VMA
     . + RKZF(L1,Q1,P1,RMB,-LG1,K1)    *RKZF(-LG1,R1,RR1,RMB,LG2,R2)
     .  *RKZF(LG2,K2,RR2,RMB,LV,LEP2)  *RKZF(LV,LEP1,P2,RMB,L2,Q2)*VPA
     . + RKZF(L1,Q1,P1,RMB,-LG1,K1)    *RKZF(-LG1,R1,RR1,RMB,LG2,R2)
     .  *RKZF(LG2,K2,RR2,RMB,-LV,LEP1) *RKZF(-LV,LEP2,P2,RMB,L2,Q2)*VMA
     . + RKZF(L1,Q1,P1,RMB,-LG1,K1)    *RKZF(-LG1,R1,RR1,RMB,-LG2,K2)
     .  *RKZF(-LG2,R2,RR2,RMB,LV,LEP2) *RKZF(LV,LEP1,P2,RMB,L2,Q2)*VPA
     . + RKZF(L1,Q1,P1,RMB,-LG1,K1)    *RKZF(-LG1,R1,RR1,RMB,-LG2,K2)
     .  *RKZF(-LG2,R2,RR2,RMB,-LV,LEP1)*RKZF(-LV,LEP2,P2,RMB,L2,Q2)*VMA

      DO 290 K=0,3
        RR1(K)=P1(K)-K2(K)
        RR2(K)=RR1(K)-K1(K)
  290 CONTINUE
      RR1(4)=5D0
      RR2(4)=4D0
      ZD21V=(RKDOT(RR1,RR1)-RMB**2)*(RKDOT(RR2,RR2)-RMB**2)
      ZN21V =
     .   RKZF(L1,Q1,P1,RMB,LG2,R2)     *RKZF(LG2,K2,RR1,RMB,LG1,R1)
     .  *RKZF(LG1,K1,RR2,RMB,LV,LEP2)  *RKZF(LV,LEP1,P2,RMB,L2,Q2)*VPA
     . + RKZF(L1,Q1,P1,RMB,LG2,R2)     *RKZF(LG2,K2,RR1,RMB,LG1,R1)
     .  *RKZF(LG1,K1,RR2,RMB,-LV,LEP1) *RKZF(-LV,LEP2,P2,RMB,L2,Q2)*VMA
     . + RKZF(L1,Q1,P1,RMB,LG2,R2)     *RKZF(LG2,K2,RR1,RMB,-LG1,K1)
     .  *RKZF(-LG1,R1,RR2,RMB,LV,LEP2) *RKZF(LV,LEP1,P2,RMB,L2,Q2)*VPA
     . + RKZF(L1,Q1,P1,RMB,LG2,R2)     *RKZF(LG2,K2,RR1,RMB,-LG1,K1)
     .  *RKZF(-LG1,R1,RR2,RMB,-LV,LEP1)*RKZF(-LV,LEP2,P2,RMB,L2,Q2)*VMA
     . + RKZF(L1,Q1,P1,RMB,-LG2,K2)    *RKZF(-LG2,R2,RR1,RMB,LG1,R1)
     .  *RKZF(LG1,K1,RR2,RMB,LV,LEP2)  *RKZF(LV,LEP1,P2,RMB,L2,Q2)*VPA
     . + RKZF(L1,Q1,P1,RMB,-LG2,K2)    *RKZF(-LG2,R2,RR1,RMB,LG1,R1)
     .  *RKZF(LG1,K1,RR2,RMB,-LV,LEP1) *RKZF(-LV,LEP2,P2,RMB,L2,Q2)*VMA
     . + RKZF(L1,Q1,P1,RMB,-LG2,K2)    *RKZF(-LG2,R2,RR1,RMB,-LG1,K1)
     .  *RKZF(-LG1,R1,RR2,RMB,LV,LEP2) *RKZF(LV,LEP1,P2,RMB,L2,Q2)*VPA
     . + RKZF(L1,Q1,P1,RMB,-LG2,K2)    *RKZF(-LG2,R2,RR1,RMB,-LG1,K1)
     .  *RKZF(-LG1,R1,RR2,RMB,-LV,LEP1)*RKZF(-LV,LEP2,P2,RMB,L2,Q2)*VMA

      DO 300 K=0,3
        RR1(K)=P1(K)-K1(K)
        RR2(K)=RR1(K)+QV(K)
  300 CONTINUE
      RR1(4)=3D0
      RR2(4)=6D0
      ZD1V2=(RKDOT(RR1,RR1)-RMB**2)*(RKDOT(RR2,RR2)-RMB**2)
      ZN1V2 =
     .   RKZF(L1,Q1,P1,RMB,LG1,R1)     *RKZF(LG1,K1,RR1,RMB,LV,LEP2)
     .  *RKZF(LV,LEP1,RR2,RMB,LG2,R2)  *RKZF(LG2,K2,P2,RMB,L2,Q2)*VPA
     . + RKZF(L1,Q1,P1,RMB,LG1,R1)     *RKZF(LG1,K1,RR1,RMB,LV,LEP2)
     .  *RKZF(LV,LEP1,RR2,RMB,-LG2,K2) *RKZF(-LG2,R2,P2,RMB,L2,Q2)*VPA
     . + RKZF(L1,Q1,P1,RMB,LG1,R1)     *RKZF(LG1,K1,RR1,RMB,-LV,LEP1)
     .  *RKZF(-LV,LEP2,RR2,RMB,LG2,R2) *RKZF(LG2,K2,P2,RMB,L2,Q2)*VMA
     . + RKZF(L1,Q1,P1,RMB,LG1,R1)     *RKZF(LG1,K1,RR1,RMB,-LV,LEP1)
     .  *RKZF(-LV,LEP2,RR2,RMB,-LG2,K2)*RKZF(-LG2,R2,P2,RMB,L2,Q2)*VMA
     . + RKZF(L1,Q1,P1,RMB,-LG1,K1)    *RKZF(-LG1,R1,RR1,RMB,LV,LEP2)
     .  *RKZF(LV,LEP1,RR2,RMB,LG2,R2)  *RKZF(LG2,K2,P2,RMB,L2,Q2)*VPA
     . + RKZF(L1,Q1,P1,RMB,-LG1,K1)    *RKZF(-LG1,R1,RR1,RMB,LV,LEP2)
     .  *RKZF(LV,LEP1,RR2,RMB,-LG2,K2) *RKZF(-LG2,R2,P2,RMB,L2,Q2)*VPA
     . + RKZF(L1,Q1,P1,RMB,-LG1,K1)    *RKZF(-LG1,R1,RR1,RMB,-LV,LEP1)
     .  *RKZF(-LV,LEP2,RR2,RMB,LG2,R2) *RKZF(LG2,K2,P2,RMB,L2,Q2)*VMA
     . + RKZF(L1,Q1,P1,RMB,-LG1,K1)    *RKZF(-LG1,R1,RR1,RMB,-LV,LEP1)
     .  *RKZF(-LV,LEP2,RR2,RMB,-LG2,K2)*RKZF(-LG2,R2,P2,RMB,L2,Q2)*VMA

      DO 310 K=0,3
        RR1(K)=P1(K)-K2(K)
        RR2(K)=RR1(K)+QV(K)
  310 CONTINUE
      RR1(4)=5D0
      RR2(4)=7D0
      ZD2V1=(RKDOT(RR1,RR1)-RMB**2)*(RKDOT(RR2,RR2)-RMB**2)
      ZN2V1 =
     .   RKZF(L1,Q1,P1,RMB,LG2,R2)     *RKZF(LG2,K2,RR1,RMB,LV,LEP2)
     .  *RKZF(LV,LEP1,RR2,RMB,LG1,R1)  *RKZF(LG1,K1,P2,RMB,L2,Q2)*VPA
     . + RKZF(L1,Q1,P1,RMB,LG2,R2)     *RKZF(LG2,K2,RR1,RMB,LV,LEP2)
     .  *RKZF(LV,LEP1,RR2,RMB,-LG1,K1) *RKZF(-LG1,R1,P2,RMB,L2,Q2)*VPA
     . + RKZF(L1,Q1,P1,RMB,LG2,R2)     *RKZF(LG2,K2,RR1,RMB,-LV,LEP1)
     .  *RKZF(-LV,LEP2,RR2,RMB,LG1,R1) *RKZF(LG1,K1,P2,RMB,L2,Q2)*VMA
     . + RKZF(L1,Q1,P1,RMB,LG2,R2)     *RKZF(LG2,K2,RR1,RMB,-LV,LEP1)
     .  *RKZF(-LV,LEP2,RR2,RMB,-LG1,K1)*RKZF(-LG1,R1,P2,RMB,L2,Q2)*VMA
     . + RKZF(L1,Q1,P1,RMB,-LG2,K2)    *RKZF(-LG2,R2,RR1,RMB,LV,LEP2)
     .  *RKZF(LV,LEP1,RR2,RMB,LG1,R1)  *RKZF(LG1,K1,P2,RMB,L2,Q2)*VPA
     . + RKZF(L1,Q1,P1,RMB,-LG2,K2)    *RKZF(-LG2,R2,RR1,RMB,LV,LEP2)
     .  *RKZF(LV,LEP1,RR2,RMB,-LG1,K1) *RKZF(-LG1,R1,P2,RMB,L2,Q2)*VPA
     . + RKZF(L1,Q1,P1,RMB,-LG2,K2)    *RKZF(-LG2,R2,RR1,RMB,-LV,LEP1)
     .  *RKZF(-LV,LEP2,RR2,RMB,LG1,R1) *RKZF(LG1,K1,P2,RMB,L2,Q2)*VMA
     . + RKZF(L1,Q1,P1,RMB,-LG2,K2)    *RKZF(-LG2,R2,RR1,RMB,-LV,LEP1)
     .  *RKZF(-LV,LEP2,RR2,RMB,-LG1,K1)*RKZF(-LG1,R1,P2,RMB,L2,Q2)*VMA

      DO 320 K=0,3
        RR1(K)=P1(K)+QV(K)
        RR2(K)=RR1(K)-K1(K)
  320 CONTINUE
      RR1(4)=8D0
      RR2(4)=6D0
      ZDV12=(RKDOT(RR1,RR1)-RMB**2)*(RKDOT(RR2,RR2)-RMB**2)
      ZNV12 =
     .   RKZF(L1,Q1,P1,RMB,LV,LEP2)   *RKZF(LV,LEP1,RR1,RMB,LG1,R1)
     .  *RKZF(LG1,K1,RR2,RMB,LG2,R2)  *RKZF(LG2,K2,P2,RMB,L2,Q2)*VPA
     . + RKZF(L1,Q1,P1,RMB,LV,LEP2)   *RKZF(LV,LEP1,RR1,RMB,LG1,R1)
     .  *RKZF(LG1,K1,RR2,RMB,-LG2,K2) *RKZF(-LG2,R2,P2,RMB,L2,Q2)*VPA
     . + RKZF(L1,Q1,P1,RMB,LV,LEP2)   *RKZF(LV,LEP1,RR1,RMB,-LG1,K1)
     .  *RKZF(-LG1,R1,RR2,RMB,LG2,R2) *RKZF(LG2,K2,P2,RMB,L2,Q2)*VPA
     . + RKZF(L1,Q1,P1,RMB,LV,LEP2)   *RKZF(LV,LEP1,RR1,RMB,-LG1,K1)
     .  *RKZF(-LG1,R1,RR2,RMB,-LG2,K2)*RKZF(-LG2,R2,P2,RMB,L2,Q2)*VPA
     . + RKZF(L1,Q1,P1,RMB,-LV,LEP1)  *RKZF(-LV,LEP2,RR1,RMB,LG1,R1)
     .  *RKZF(LG1,K1,RR2,RMB,LG2,R2)  *RKZF(LG2,K2,P2,RMB,L2,Q2)*VMA
     . + RKZF(L1,Q1,P1,RMB,-LV,LEP1)  *RKZF(-LV,LEP2,RR1,RMB,LG1,R1)
     .  *RKZF(LG1,K1,RR2,RMB,-LG2,K2) *RKZF(-LG2,R2,P2,RMB,L2,Q2)*VMA
     . + RKZF(L1,Q1,P1,RMB,-LV,LEP1)  *RKZF(-LV,LEP2,RR1,RMB,-LG1,K1)
     .  *RKZF(-LG1,R1,RR2,RMB,LG2,R2) *RKZF(LG2,K2,P2,RMB,L2,Q2)*VMA
     . + RKZF(L1,Q1,P1,RMB,-LV,LEP1)  *RKZF(-LV,LEP2,RR1,RMB,-LG1,K1)
     .  *RKZF(-LG1,R1,RR2,RMB,-LG2,K2)*RKZF(-LG2,R2,P2,RMB,L2,Q2)*VMA

      DO 330 K=0,3
        RR1(K)=P1(K)+QV(K)
        RR2(K)=RR1(K)-K2(K)
  330 CONTINUE
      RR1(4)=8D0
      RR2(4)=7D0
      ZDV21=(RKDOT(RR1,RR1)-RMB**2)*(RKDOT(RR2,RR2)-RMB**2)
      ZNV21 =
     .   RKZF(L1,Q1,P1,RMB,LV,LEP2)   *RKZF(LV,LEP1,RR1,RMB,LG2,R2)
     .  *RKZF(LG2,K2,RR2,RMB,LG1,R1)  *RKZF(LG1,K1,P2,RMB,L2,Q2)*VPA
     . + RKZF(L1,Q1,P1,RMB,LV,LEP2)   *RKZF(LV,LEP1,RR1,RMB,LG2,R2)
     .  *RKZF(LG2,K2,RR2,RMB,-LG1,K1) *RKZF(-LG1,R1,P2,RMB,L2,Q2)*VPA
     . + RKZF(L1,Q1,P1,RMB,LV,LEP2)   *RKZF(LV,LEP1,RR1,RMB,-LG2,K2)
     .  *RKZF(-LG2,R2,RR2,RMB,LG1,R1) *RKZF(LG1,K1,P2,RMB,L2,Q2)*VPA
     . + RKZF(L1,Q1,P1,RMB,LV,LEP2)   *RKZF(LV,LEP1,RR1,RMB,-LG2,K2)
     .  *RKZF(-LG2,R2,RR2,RMB,-LG1,K1)*RKZF(-LG1,R1,P2,RMB,L2,Q2)*VPA
     . + RKZF(L1,Q1,P1,RMB,-LV,LEP1)  *RKZF(-LV,LEP2,RR1,RMB,LG2,R2)
     .  *RKZF(LG2,K2,RR2,RMB,LG1,R1)  *RKZF(LG1,K1,P2,RMB,L2,Q2)*VMA
     . + RKZF(L1,Q1,P1,RMB,-LV,LEP1)  *RKZF(-LV,LEP2,RR1,RMB,LG2,R2)
     .  *RKZF(LG2,K2,RR2,RMB,-LG1,K1) *RKZF(-LG1,R1,P2,RMB,L2,Q2)*VMA
     . + RKZF(L1,Q1,P1,RMB,-LV,LEP1)  *RKZF(-LV,LEP2,RR1,RMB,-LG2,K2)
     .  *RKZF(-LG2,R2,RR2,RMB,LG1,R1) *RKZF(LG1,K1,P2,RMB,L2,Q2)*VMA
     . + RKZF(L1,Q1,P1,RMB,-LV,LEP1)  *RKZF(-LV,LEP2,RR1,RMB,-LG2,K2)
     .  *RKZF(-LG2,R2,RR2,RMB,-LG1,K1)*RKZF(-LG1,R1,P2,RMB,L2,Q2)*VMA

* COMPUTE THE DIAGRAMS SO FAR
      ZDIA1=ZN12V/ZD12V
      ZDIA2=ZN21V/ZD21V
      ZDIA3=ZN1V2/ZD1V2
      ZDIA4=ZN2V1/ZD2V1
      ZDIA5=ZNV12/ZDV12
      ZDIA6=ZNV21/ZDV21

* SECOND PART OF THE RESULT: THE NONABELIAN PART.
* THIS IS MADE UP PARTLY FROM THE ABELIAN PART AND PARTLY FROM
* NEW PIECES
* THE ASSIGNMENT OF THE 4TH COMPONENT LABELS IS NOW UNNECESSARY
* FOR RR1 SINCE IT DOES NOT OCCUR IN ANY ZF HERE

      S=2D0*RKDOT(K1,K2)

      DO 340 K=0,3
        RR1(K)=PP2(K)+QV(K)
  340 CONTINUE
      ZD11=S*(RKDOT(RR1,RR1)-RMB**2)

      ZC12V =
     . + RKZF(L1,Q1,P1,RMB,LG1,R1) *RKZSF(LG1,K1,LG2,R2)
     .  *RKZSF(LG2,K2,LV,LEP2)  *RKZF(LV,LEP1,P2,RMB,L2,Q2)*VPA
     . + RKZF(L1,Q1,P1,RMB,LG1,R1) *RKZSF(LG1,K1,LG2,R2)
     .  *RKZSF(LG2,K2,-LV,LEP1) *RKZF(-LV,LEP2,P2,RMB,L2,Q2)*VMA
     . + RKZF(L1,Q1,P1,RMB,LG1,R1) *RKZSF(LG1,K1,-LG2,K2)
     .  *RKZSF(-LG2,R2,LV,LEP2) *RKZF(LV,LEP1,P2,RMB,L2,Q2)*VPA
     . + RKZF(L1,Q1,P1,RMB,LG1,R1) *RKZSF(LG1,K1,-LG2,K2)
     .  *RKZSF(-LG2,R2,-LV,LEP1)*RKZF(-LV,LEP2,P2,RMB,L2,Q2)*VMA
     . + RKZF(L1,Q1,P1,RMB,-LG1,K1)*RKZSF(-LG1,R1,LG2,R2)
     .  *RKZSF(LG2,K2,LV,LEP2)  *RKZF(LV,LEP1,P2,RMB,L2,Q2)*VPA
     . + RKZF(L1,Q1,P1,RMB,-LG1,K1)*RKZSF(-LG1,R1,LG2,R2)
     .  *RKZSF(LG2,K2,-LV,LEP1) *RKZF(-LV,LEP2,P2,RMB,L2,Q2)*VMA
     . + RKZF(L1,Q1,P1,RMB,-LG1,K1)*RKZSF(-LG1,R1,-LG2,K2)
     .  *RKZSF(-LG2,R2,LV,LEP2) *RKZF(LV,LEP1,P2,RMB,L2,Q2)*VPA
     . + RKZF(L1,Q1,P1,RMB,-LG1,K1)*RKZSF(-LG1,R1,-LG2,K2)
     .  *RKZSF(-LG2,R2,-LV,LEP1)*RKZF(-LV,LEP2,P2,RMB,L2,Q2)*VMA

      ZC21V =
     . + RKZF(L1,Q1,P1,RMB,LG2,R2) *RKZSF(LG2,K2,LG1,R1)
     .  *RKZSF(LG1,K1,LV,LEP2)  *RKZF(LV,LEP1,P2,RMB,L2,Q2)*VPA
     . + RKZF(L1,Q1,P1,RMB,LG2,R2) *RKZSF(LG2,K2,LG1,R1)
     .  *RKZSF(LG1,K1,-LV,LEP1) *RKZF(-LV,LEP2,P2,RMB,L2,Q2)*VMA
     . + RKZF(L1,Q1,P1,RMB,LG2,R2) *RKZSF(LG2,K2,-LG1,K1)
     .  *RKZSF(-LG1,R1,LV,LEP2) *RKZF(LV,LEP1,P2,RMB,L2,Q2)*VPA
     . + RKZF(L1,Q1,P1,RMB,LG2,R2) *RKZSF(LG2,K2,-LG1,K1)
     .  *RKZSF(-LG1,R1,-LV,LEP1)*RKZF(-LV,LEP2,P2,RMB,L2,Q2)*VMA
     . + RKZF(L1,Q1,P1,RMB,-LG2,K2)*RKZSF(-LG2,R2,LG1,R1)
     .  *RKZSF(LG1,K1,LV,LEP2)  *RKZF(LV,LEP1,P2,RMB,L2,Q2)*VPA
     . + RKZF(L1,Q1,P1,RMB,-LG2,K2)*RKZSF(-LG2,R2,LG1,R1)
     .  *RKZSF(LG1,K1,-LV,LEP1) *RKZF(-LV,LEP2,P2,RMB,L2,Q2)*VMA
     . + RKZF(L1,Q1,P1,RMB,-LG2,K2)*RKZSF(-LG2,R2,-LG1,K1)
     .  *RKZSF(-LG1,R1,LV,LEP2) *RKZF(LV,LEP1,P2,RMB,L2,Q2)*VPA
     . + RKZF(L1,Q1,P1,RMB,-LG2,K2)*RKZSF(-LG2,R2,-LG1,K1)
     .  *RKZSF(-LG1,R1,-LV,LEP1)*RKZF(-LV,LEP2,P2,RMB,L2,Q2)*VMA
      ZDIA7=(-ZN12V+ZN21V)/ZD11-(ZC12V-ZC21V)/(2D0*S)

      DO 350 K=0,3
        RR1(K)=P1(K)+QV(K)
  350 CONTINUE
      ZD22=S*(RKDOT(RR1,RR1)-RMB**2)

      ZCV12 =
     . + RKZF(L1,Q1,P1,RMB,LV,LEP2) *RKZSF(LV,LEP1,LG1,R1)
     .  *RKZSF(LG1,K1,LG2,R2)    *RKZF(LG2,K2,P2,RMB,L2,Q2)*VPA
     . + RKZF(L1,Q1,P1,RMB,LV,LEP2) *RKZSF(LV,LEP1,LG1,R1)
     .  *RKZSF(LG1,K1,-LG2,K2)   *RKZF(-LG2,R2,P2,RMB,L2,Q2)*VPA
     . + RKZF(L1,Q1,P1,RMB,LV,LEP2) *RKZSF(LV,LEP1,-LG1,K1)
     .  *RKZSF(-LG1,R1,LG2,R2)   *RKZF(LG2,K2,P2,RMB,L2,Q2)*VPA
     . + RKZF(L1,Q1,P1,RMB,LV,LEP2) *RKZSF(LV,LEP1,-LG1,K1)
     .  *RKZSF(-LG1,R1,-LG2,K2)  *RKZF(-LG2,R2,P2,RMB,L2,Q2)*VPA
     . + RKZF(L1,Q1,P1,RMB,-LV,LEP1)*RKZSF(-LV,LEP2,LG1,R1)
     .  *RKZSF(LG1,K1,LG2,R2)    *RKZF(LG2,K2,P2,RMB,L2,Q2)*VMA
     . + RKZF(L1,Q1,P1,RMB,-LV,LEP1)*RKZSF(-LV,LEP2,LG1,R1)
     .  *RKZSF(LG1,K1,-LG2,K2)   *RKZF(-LG2,R2,P2,RMB,L2,Q2)*VMA
     . + RKZF(L1,Q1,P1,RMB,-LV,LEP1)*RKZSF(-LV,LEP2,-LG1,K1)
     .  *RKZSF(-LG1,R1,LG2,R2)   *RKZF(LG2,K2,P2,RMB,L2,Q2)*VMA
     . + RKZF(L1,Q1,P1,RMB,-LV,LEP1)*RKZSF(-LV,LEP2,-LG1,K1)
     .  *RKZSF(-LG1,R1,-LG2,K2)  *RKZF(-LG2,R2,P2,RMB,L2,Q2)*VMA

* THE FOURTH COMBINATION CAN BE GOTTEN FROM
* THE FIRST THREE USING DIRAC ALGEBRA:
* EPS1*EPS2*EPVS+EPS2*EPS1*EPSV = 2(EPS1.EPS2)*EPSV ETC.
      ZCV21=ZC12V+ZC21V-ZCV12

      ZDIA8=(-ZNV12+ZNV21)/ZD22-(ZCV12-ZCV21)/(2D0*S)

* CONSTRUCT THE ABELIAN AND NONABELIAN PART

      ZABEL= ZDIA1+ZDIA2+ZDIA3+ZDIA4+ZDIA5+ZDIA6
      ZNABEL=ZDIA1-ZDIA2+ZDIA3-ZDIA4+ZDIA5-ZDIA6
      ZNABEM=2D0*ZDIA7+2D0*ZDIA8
      ZNABEL=ZNABEL-ZNABEM
      ZABEL=ZABEL*ZFAC*ZFACL
      ZNABEL=ZNABEL*ZFAC*ZFACL

* INCLUDE COLOUR FACTORS:
* (N**2-1)*(N**2-2)/(8*N) = 7/3 FOR THE ABELIAN PART
* N*(N**2-1)/8 = 3 FOR THE NONABELIAN PART
* AND ADD THE RESULT TO THE CROSS SECTION
      THIS1=7D0/3D0*ABS(ZABEL)**2+3D0*ABS(ZNABEL)**2
CC    WRITE(6,801)HELIX,LG1,LG2,LV,L1,L2,THIS1
CC801 FORMAT(' ',6I4,D30.20)
      CROSS=CROSS+THIS1

* END OF THE BIG LOOP OVER HELICITIES
  360 CONTINUE

* DO NOT FORGET TO PUT P2 BACK TO ITS ORIGINAL VALUE IN PP2!
      DO 370 K=0,3
        P2(K)=PP2(K)
  370 CONTINUE

* ADD AVERAGING FACTORS:
* 1/2 FOR EACH GLUON SPIN, 1/8 FOR EACH GLUON COLOUR
      CROSS=CROSS/256D0

* TAKE INTO ACCOUNT A POSSIBLE FACTOR FOR THE HELICITY SUM OPTION
* AND RETURN THE FINAL RESULT
      IF(IMC.EQ.1) CROSS=CROSS*MULT
      RESULT=CROSS
      END

*==================================================================

      FUNCTION RKZF(L1,P1,Q,RMB,L2,P2)
* COMPUTES THE SCALAR STRUCTURE
* U_BAR(L1,P1)(SLASH(Q)+RMB)U(L2,P2)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMPLEX RKZF,RKZPR,RKZSF,ZFX
      COMPLEX ANSF(-1:1,1:4,1:8,-1:1,1:4)
      INTEGER DONF(-1:1,1:4,1:8,-1:1,1:4)
      COMPLEX CHECK
      COMMON / RKZFCO / ANSF,DONF
      DIMENSION P1(0:4),P2(0:4),Q(0:4),R(0:4)
* CHECK ON CORRECT LABEL INPUT
      IP1=IDINT(P1(4))
      IQ=IDINT(Q(4))
      IP2=IDINT(P2(4))
      IF(IABS(L1).NE.1.OR.IABS(L2).NE.1.OR.
     . IP1.LT.1.OR.IP1.GT.4            .OR.
     . IQ.LT.1.OR.IQ.GT.8              .OR.
     . IP2.LT.1.OR.IP2.GT.4) THEN
        WRITE(6,*) ' RKZF LABEL ERROR'
        WRITE(6,*) 'L1=',L1,' IP1=',IP1,' IQ=',IQ,
     .             ' L2=',L2,' IP2=',IP2
        STOP
      ENDIF
* CHECK WHETHER THIS ONE HAS BEEN CALCULATED ALREADY
      IF(DONF(L1,IP1,IQ,L2,IP2).EQ.0) THEN
* THIS ONE NOT DONE YET: DO IT AND STORE THE RESULT IN ARRAY 'ANSF'
        IF(L1.EQ.L2) THEN
          A=2D0*RKDOT(Q,P2)
C         IF(DABS(A).LT.(1D-10*P2(0)*Q(0))) THEN
C...The check above is extended to following.
          IF(ABS(A).LT.MAX(1D-8,ABS(1D-10*P2(0)*Q(0)))) THEN
            ANSF(L1,IP1,IQ,L2,IP2)=(0.,0.)
          ELSE
            A=RKDOT(Q,Q)/A
            DO 100 K=0,3
              R(K)=Q(K)-A*P2(K)
  100       CONTINUE
            IF(R(0).GT.0D0) THEN
              C=1D0
            ELSE
              DO 110 K=0,3
                R(K)=-R(K)
  110         CONTINUE
              C=-1D0
            ENDIF
            ANSF(L1,IP1,IQ,L2,IP2)=C*RKZPR(L1,P1,R)*RKZPR(-L1,R,P2)
          ENDIF
        ELSEIF(L1.EQ.-L2) THEN
          ANSF(L1,IP1,IQ,L2,IP2)=RMB*RKZSF(L1,P1,L2,P2)
        ELSE
          WRITE(6,*) ' ERROR IN RKZF: L1=',L1,'  L2=',L2
          STOP
        ENDIF
        RKZF=ANSF(L1,IP1,IQ,L2,IP2)
        DONF(L1,IP1,IQ,L2,IP2)=1
      ELSE
        RKZF=ANSF(L1,IP1,IQ,L2,IP2)
      ENDIF
      END

*==================================================================

      FUNCTION RKRAND(IDUMMY)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      DATA INIT/0/
      IF(INIT.EQ.0) THEN
        INIT=1
        X=DMOD(DSQRT(2D0),1D0)
        Y=DMOD(DSQRT(3D0),1D0)
        Z=DMOD(DSQRT(5D0),1D0)
      ELSE
        X=DMOD(X+Y+Z,1D0)
        Y=DMOD(X+Y+Z,1D0)
        Z=DMOD(X+Y+Z,1D0)
      ENDIF
      RKRAND=X
      END

*==================================================================

      FUNCTION RKDOT(P,Q)
      DOUBLE PRECISION P(0:4),Q(0:4),RKDOT
      RKDOT=P(0)*Q(0)-P(1)*Q(1)-P(2)*Q(2)-P(3)*Q(3)
      END

*==================================================================

      FUNCTION RKZPR(L,Q1,Q2)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMPLEX RKZPR
      DIMENSION Q1(0:4),Q2(0:4)
      IF(IABS(L).NE.1) THEN
        WRITE(6,*) ' RKZPR: ERROR   L=',L
        STOP
      ENDIF
C...Introduce cutoff to check that R1 and R2 not zero.
      R1=DSQRT(MAX(1D-10,Q1(0)-Q1(1)))
      R2=DSQRT(MAX(1D-10,Q2(0)-Q2(1)))
      RKZPR=CMPLX(SNGL(Q1(2)),SNGL(Q1(3)))*R2/R1
     .     -CMPLX(SNGL(Q2(2)),SNGL(Q2(3)))*R1/R2
      IF(L.EQ.-1) RKZPR=-CONJG(RKZPR)
      END

*==================================================================

      FUNCTION RKZSF(L1,P1,L2,P2)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMPLEX RKZSF,RKZPR,RKZSFK,CHECK
      COMPLEX ANSS(-1:1,1:4,-1:1,1:4)
      INTEGER DONS(-1:1,1:4,-1:1,1:4)
      COMMON / RKZSCO / ANSS,DONS
      DIMENSION P1(0:4),P2(0:4)
* CHECK ON CORRECT LABEL INPUT
      IP1=IDINT(P1(4))
      IP2=IDINT(P2(4))
      IF(IABS(L1).NE.1.OR.IABS(L2).NE.1.OR.
     . IP1.LT.1.OR.IP2.GT.4.OR.IP2.LT.1.OR.IP2.GT.4) THEN
       WRITE(6,*)
     .  ' RKZSF: ERROR L1=',L1,' L2=',L2,' IP1=',IP1,' IP2=',IP2
       STOP
      ENDIF
* CHECK WHETER THIS ONE WAS ALREADY COMPUTED
* DONS(,,,)=0: NOT YET COMPUTED, DONS(,,,)=1: ALREADY COMPUTED
* IF NOT YET COMPUTED: COMPUTE IT, AND STORE IN ARRAY 'ANSS'
* IF ALREADY COMPUTED: GET THE RESULT FROM ARRAY 'ANSS'
      IF(DONS(L1,IP1,L2,IP2).EQ.0) THEN
        IF(L1.EQ.L2) THEN
          ANSS(L1,IP1,L2,IP2)=(0.,0.)
        ELSE
          ANSS(L1,IP1,L2,IP2)=RKZPR(L1,P1,P2)
        ENDIF
        DONS(L1,IP1,L2,IP2)=1
      ENDIF
      RKZSF=ANSS(L1,IP1,L2,IP2)
      END

*==================================================================

      SUBROUTINE RKHLPK(NUM,LGL1,LGL2,LLV,LL1,LL2)
      IMPLICIT INTEGER(A-Z)
      SAVE
      DIMENSION CONFIG(32,6)
      DATA INIT/0/
      IF(INIT.EQ.0) THEN
        INIT=1
        MUM=0
        DO 140 GL1=1,-1,-2
          DO 130 GL2=1,-1,-2
            DO 120 LV=1,-1,-2
              DO 110 L1=1,-1,-2
                DO 100 L2=1,-1,-2
                  MUM=MUM+1
                  CONFIG(MUM,1)=GL1
                  CONFIG(MUM,2)=GL2
                  CONFIG(MUM,3)=LV
                  CONFIG(MUM,4)=L1
                  CONFIG(MUM,5)=L2
  100           CONTINUE
  110         CONTINUE
  120       CONTINUE
  130     CONTINUE
  140   CONTINUE
      ENDIF
      LGL1=CONFIG(NUM,1)
      LGL2=CONFIG(NUM,2)
      LLV =CONFIG(NUM,3)
      LL1 =CONFIG(NUM,4)
      LL2 =CONFIG(NUM,5)
      END

C*********************************************************************

      SUBROUTINE RYKCUT(MCUT)

C...Dummy routine, which the user can replace in order to make cuts on
C...the kinematics on the parton level before the matrix elements are
C...evaluated and the event is generated. The cross-section estimates
C...will automatically take these cuts into account, so the given
C...values are for the allowed phase space region only. MCUT=0 means
C...that the event has passed the cuts, MCUT=1 that it has failed.
      COMMON/RYINT1/MINT(400),VINT(400)
      COMMON/RYINT2/ISET(200),KFPR(200,2),COEF(200,20),ICOL(40,4,2)
      SAVE /RYINT1/,/RYINT2/

C...Set default value (accepting event) for MCUT.
      MCUT=0

C...Read out subprocess number.
      ISUB=MINT(1)
      ISTSB=ISET(ISUB)

C...Read out tau, y*, cos(theta), tau' (where defined, else =0).
      TAU=VINT(21)
      YST=VINT(22)
      CTH=0.
      IF(ISTSB.EQ.2.OR.ISTSB.EQ.4.OR.ISTSB.EQ.6) CTH=VINT(23)
      TAUP=0.
      IF(ISTSB.GE.3.AND.ISTSB.LE.5) TAUP=VINT(26)

C...Calculate x_1, x_2, x_F.
      IF(ISTSB.LE.2.OR.ISTSB.GE.6) THEN
        X1=SQRT(TAU)*EXP(YST)
        X2=SQRT(TAU)*EXP(-YST)
      ELSE
        X1=SQRT(TAUP)*EXP(YST)
        X2=SQRT(TAUP)*EXP(-YST)
      ENDIF
      XF=X1-X2

C...Calculate shat, that, uhat, p_T^2.
      SHAT=TAU*VINT(2)
      SQM3=VINT(63)
      SQM4=VINT(64)
      RM3=SQM3/SHAT
      RM4=SQM4/SHAT
      BE34=SQRT(MAX(0.,(1.-RM3-RM4)**2-4.*RM3*RM4))
      RPTS=4.*VINT(71)**2/SHAT
      BE34L=SQRT(MAX(0.,(1.-RM3-RM4)**2-4.*RM3*RM4-RPTS))
      RM34=2.*RM3*RM4
      RSQM=1.+RM34
      RTHM=(4.*RM3*RM4+RPTS)/(1.-RM3-RM4+BE34L)
      THAT=-0.5*SHAT*MAX(RTHM,1.-RM3-RM4-BE34*CTH)
      UHAT=-0.5*SHAT*MAX(RTHM,1.-RM3-RM4+BE34*CTH)
      PT2=MAX(VINT(71)**2,0.25*SHAT*BE34**2*(1.-CTH**2))

C...Decisions by user to be put here.

      RETURN
      END

C*********************************************************************

      SUBROUTINE RYSTFE(KF,X,Q2,XPQ)

C...This is a dummy routine, where the user can introduce an interface
C...to his own external structure function parametrization.
C...Arguments in:
C...KF : 11 for e-, 22 for gamma, 211 for pi+, 2212 for p.
C...    Isospin conjugation for n and charge conjugation for
C...    e+, pi-, pbar and nbar is performed by RYSTFU.
C...X : x value.
C...Q2 : Q^2 value.
C...Arguments out:
C...XPQ(-25:25) : x * f(x,Q^2), with index according to KF code,
C...    except that gluon is placed in 0. Thus XPQ(0) = xg,
C...    XPQ(1) = xd, XPQ(-1) = xdbar, XPQ(2) = xu, XPQ(-2) = xubar,
C...    XPQ(3) = xs, XPQ(-3) = xsbar, XPQ(4) = xc, XPQ(-4) = xcbar,
C...    XPQ(5) = xb, XPQ(-5) = xbbar, XPQ(6) = xt, XPQ(-6) = xtbar,
C...    XPQ(11) = xe-, XPQ(-11) = xe+, XPQ(22) = xgamma.
C...
C...One such interface, to the Diemos, Ferroni, Longo, Martinelli
C...proton structure functions, already comes with the package. What
C...the user needs here is external files with the three routines
C...FXG160, FXG260 and FXG360 of the authors above, plus the
C...interpolation routine FINT, which is part of the CERN library
C...KERNLIB package. To avoid problems with unresolved external
C...references, the external calls are commented in the current
C...version. To enable this option, remove the C* at the beginning
C...of the relevant lines.
C...
C...Alternatively, the routine can be used as an interface to the
C...structure function evolution program of Tung. This can be achieved
C...by removing C* at the beginning of some of the lines below.
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/RYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      SAVE /LUDAT1/,/LUDAT2/
      SAVE /RYPARS/
      DIMENSION XPQ(-25:25),XFDFLM(9)
      CHARACTER CHDFLM(9)*5,HEADER*40
      DATA CHDFLM/'UPVAL','DOVAL','GLUON','QBAR ','UBAR ','SBAR ',
     &'CBAR ','BBAR ','TBAR '/
      DATA HEADER/'Tung evolution package has been invoked'/
      DATA INIT/0/

C...Proton structure functions from Diemoz, Ferroni, Longo, Martinelli.
C...Allowed variable range 10 GeV2 < Q2 < 1E8 GeV2, 5E-5 < x < .95.
      IF(MSTP(51).GE.11.AND.MSTP(51).LE.13.AND.MSTP(52).LE.1) THEN
        XDFLM=MAX(0.51E-4,X)
        Q2DFLM=MAX(10.,MIN(1E8,Q2))
        IF(MSTP(52).EQ.0) Q2DFLM=10.
        DO 100 J=1,9
        IF(MSTP(52).EQ.1.AND.J.EQ.9) THEN
          Q2DFLM=Q2DFLM*(40./PMAS(6,1))**2
          Q2DFLM=MAX(10.,MIN(1E8,Q2))
        ENDIF
        XFDFLM(J)=0.
C...Remove C* on following three lines to enable the DFLM options.
C*      IF(MSTP(51).EQ.11) CALL FXG160(XDFLM,Q2DFLM,CHDFLM(J),XFDFLM(J))
C*      IF(MSTP(51).EQ.12) CALL FXG260(XDFLM,Q2DFLM,CHDFLM(J),XFDFLM(J))
C*      IF(MSTP(51).EQ.13) CALL FXG360(XDFLM,Q2DFLM,CHDFLM(J),XFDFLM(J))
  100   CONTINUE
        IF(X.LT.0.51E-4.AND.ABS(PARP(51)-1.).GT.0.01) THEN
          CXS=(0.51E-4/X)**(PARP(51)-1.)
          DO 110 J=1,7
  110     XFDFLM(J)=XFDFLM(J)*CXS
        ENDIF
        XPQ(0)=XFDFLM(3)
        XPQ(1)=XFDFLM(2)+XFDFLM(5)
        XPQ(2)=XFDFLM(1)+XFDFLM(5)
        XPQ(3)=XFDFLM(6)
        XPQ(4)=XFDFLM(7)
        XPQ(5)=XFDFLM(8)
        XPQ(6)=XFDFLM(9)
        XPQ(-1)=XFDFLM(5)
        XPQ(-2)=XFDFLM(5)
        XPQ(-3)=XFDFLM(6)
        XPQ(-4)=XFDFLM(7)
        XPQ(-5)=XFDFLM(8)
        XPQ(-6)=XFDFLM(9)

C...Proton structure function evolution from Wu-Ki Tung: parton
C...distribution functions incorporating heavy quark mass effects.
C...Allowed variable range: PARP(52) < Q < PARP(53); PARP(54) < x < 1.
      ELSE
        IF(INIT.EQ.0) THEN
          I1=0
          IF(MSTP(52).EQ.4) I1=1
          IHDRN=1
          NU=MSTP(53)
          I2=MSTP(51)
          IF(MSTP(51).GE.11) I2=MSTP(51)-3
          I3=0
          IF(MSTP(52).EQ.3) I3=1

C...Convert to Lambda in CWZ scheme (approximately linear relation).
          ALAM=0.75*PARP(1)
          TPMS=PMAS(6,1)
          QINI=PARP(52)
          QMAX=PARP(53)
          XMIN=PARP(54)

C...Initialize evolution (perform calculation or read results from
C...file).
C...Remove C* on following two lines to enable Tung initialization.
C*        CALL PDFSET(I1,IHDRN,ALAM,TPMS,QINI,QMAX,XMIN,NU,HEADER,
C*   &    I2,I3,IRET,IRR)
          INIT=1
        ENDIF

C...Put into output array.
        Q=SQRT(Q2)
        DO 120 I=-6,6
        FIXQ=0.
C...Remove C* on following line to enable structure function call.
C*      FIXQ=MAX(0.,PDF(10,1,I,X,Q,IR))
  120   XPQ(I)=X*FIXQ

C...Change order of u and d quarks from Tung to RYTHIA convention.
        XPS=XPQ(1)
        XPQ(1)=XPQ(2)
        XPQ(2)=XPS
        XPS=XPQ(-1)
        XPQ(-1)=XPQ(-2)
        XPQ(-2)=XPS
      ENDIF

      RETURN
      END

C************ THIS IS THE END OF RYTHIA PACKAGE ************************
