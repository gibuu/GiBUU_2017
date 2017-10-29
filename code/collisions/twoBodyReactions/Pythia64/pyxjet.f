 
C*********************************************************************
 
C...PYXJET
C...Selects number of jets in matrix element approach.
 
      SUBROUTINE PYXJET(ECM,NJET,CUT)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /PYDAT1/
C...Local array and data.
      DIMENSION ZHUT(5)
      DATA ZHUT/3.0922D0, 6.2291D0, 7.4782D0, 7.8440D0, 8.2560D0/
 
C...Trivial result for two-jets only, including parton shower.
      IF(MSTJ(101).EQ.0.OR.MSTJ(101).EQ.5) THEN
        CUT=0D0
 
C...QCD and Abelian vector gluon theory: Q^2 for jet rate and R.
      ELSEIF(MSTJ(109).EQ.0.OR.MSTJ(109).EQ.2) THEN
        CF=4D0/3D0
        IF(MSTJ(109).EQ.2) CF=1D0
        IF(MSTJ(111).EQ.0) THEN
          Q2=ECM**2
          Q2R=ECM**2
        ELSEIF(MSTU(111).EQ.0) THEN
          PARJ(169)=MIN(1D0,PARJ(129))
          Q2=PARJ(169)*ECM**2
          PARJ(168)=MIN(1D0,MAX(PARJ(128),EXP(-12D0*PARU(1)/
     &    ((33D0-2D0*MSTU(112))*PARU(111)))))
          Q2R=PARJ(168)*ECM**2
        ELSE
          PARJ(169)=MIN(1D0,MAX(PARJ(129),(2D0*PARU(112)/ECM)**2))
          Q2=PARJ(169)*ECM**2
          PARJ(168)=MIN(1D0,MAX(PARJ(128),PARU(112)/ECM,
     &    (2D0*PARU(112)/ECM)**2))
          Q2R=PARJ(168)*ECM**2
        ENDIF
 
C...alpha_strong for R and R itself.
        ALSPI=(3D0/4D0)*CF*PYALPS(Q2R)/PARU(1)
        IF(IABS(MSTJ(101)).EQ.1) THEN
          RQCD=1D0+ALSPI
        ELSEIF(MSTJ(109).EQ.0) THEN
          RQCD=1D0+ALSPI+(1.986D0-0.115D0*MSTU(118))*ALSPI**2
          IF(MSTJ(111).EQ.1) RQCD=MAX(1D0,RQCD+
     &    (33D0-2D0*MSTU(112))/12D0*LOG(PARJ(168))*ALSPI**2)
        ELSE
          RQCD=1D0+ALSPI-(3D0/32D0+0.519D0*MSTU(118))*(4D0*ALSPI/3D0)**2
        ENDIF
 
C...alpha_strong for jet rate. Initial value for y cut.
        ALSPI=(3D0/4D0)*CF*PYALPS(Q2)/PARU(1)
        CUT=MAX(0.001D0,PARJ(125),(PARJ(126)/ECM)**2)
        IF(IABS(MSTJ(101)).LE.1.OR.(MSTJ(109).EQ.0.AND.MSTJ(111).EQ.0))
     &  CUT=MAX(CUT,EXP(-SQRT(0.75D0/ALSPI))/2D0)
        IF(MSTJ(110).EQ.2) CUT=MAX(0.01D0,MIN(0.05D0,CUT))
 
C...Parametrization of first order three-jet cross-section.
  100   IF(MSTJ(101).EQ.0.OR.CUT.GE.0.25D0) THEN
          PARJ(152)=0D0
        ELSE
          PARJ(152)=(2D0*ALSPI/3D0)*((3D0-6D0*CUT+2D0*LOG(CUT))*
     &    LOG(CUT/(1D0-2D0*CUT))+(2.5D0+1.5D0*CUT-6.571D0)*
     &    (1D0-3D0*CUT)+5.833D0*(1D0-3D0*CUT)**2-3.894D0*
     &    (1D0-3D0*CUT)**3+1.342D0*(1D0-3D0*CUT)**4)/RQCD
          IF(MSTJ(109).EQ.2.AND.(MSTJ(101).EQ.2.OR.MSTJ(101).LE.-2))
     &    PARJ(152)=0D0
        ENDIF
 
C...Parametrization of second order three-jet cross-section.
        IF(IABS(MSTJ(101)).LE.1.OR.MSTJ(101).EQ.3.OR.MSTJ(109).EQ.2.OR.
     &  CUT.GE.0.25D0) THEN
          PARJ(153)=0D0
        ELSEIF(MSTJ(110).LE.1) THEN
          CT=LOG(1D0/CUT-2D0)
          PARJ(153)=ALSPI**2*CT**2*(2.419D0+0.5989D0*CT+0.6782D0*CT**2-
     &    0.2661D0*CT**3+0.01159D0*CT**4)/RQCD
 
C...Interpolation in second/first order ratio for Zhu parametrization.
        ELSEIF(MSTJ(110).EQ.2) THEN
          IZA=0
          DO 110 IY=1,5
            IF(ABS(CUT-0.01D0*IY).LT.0.0001D0) IZA=IY
  110     CONTINUE
          IF(IZA.NE.0) THEN
            ZHURAT=ZHUT(IZA)
          ELSE
            IZ=100D0*CUT
            ZHURAT=ZHUT(IZ)+(100D0*CUT-IZ)*(ZHUT(IZ+1)-ZHUT(IZ))
          ENDIF
          PARJ(153)=ALSPI*PARJ(152)*ZHURAT
        ENDIF
 
C...Shift in second order three-jet cross-section with optimized Q^2.
        IF(MSTJ(111).EQ.1.AND.IABS(MSTJ(101)).GE.2.AND.MSTJ(101).NE.3
     &  .AND.CUT.LT.0.25D0) PARJ(153)=PARJ(153)+
     &  (33D0-2D0*MSTU(112))/12D0*LOG(PARJ(169))*ALSPI*PARJ(152)
 
C...Parametrization of second order four-jet cross-section.
        IF(IABS(MSTJ(101)).LE.1.OR.CUT.GE.0.125D0) THEN
          PARJ(154)=0D0
        ELSE
          CT=LOG(1D0/CUT-5D0)
          IF(CUT.LE.0.018D0) THEN
            XQQGG=6.349D0-4.330D0*CT+0.8304D0*CT**2
            IF(MSTJ(109).EQ.2) XQQGG=(4D0/3D0)**2*(3.035D0-2.091D0*CT+
     &      0.4059D0*CT**2)
            XQQQQ=1.25D0*(-0.1080D0+0.01486D0*CT+0.009364D0*CT**2)
            IF(MSTJ(109).EQ.2) XQQQQ=8D0*XQQQQ
          ELSE
            XQQGG=-0.09773D0+0.2959D0*CT-0.2764D0*CT**2+0.08832D0*CT**3
            IF(MSTJ(109).EQ.2) XQQGG=(4D0/3D0)**2*(-0.04079D0+
     &      0.1340D0*CT-0.1326D0*CT**2+0.04365D0*CT**3)
            XQQQQ=1.25D0*(0.003661D0-0.004888D0*CT-0.001081D0*CT**2+
     &      0.002093D0*CT**3)
            IF(MSTJ(109).EQ.2) XQQQQ=8D0*XQQQQ
          ENDIF
          PARJ(154)=ALSPI**2*CT**2*(XQQGG+XQQQQ)/RQCD
          PARJ(155)=XQQQQ/(XQQGG+XQQQQ)
        ENDIF
 
C...If negative three-jet rate, change y' optimization parameter.
        IF(MSTJ(111).EQ.1.AND.PARJ(152)+PARJ(153).LT.0D0.AND.
     &  PARJ(169).LT.0.99D0) THEN
          PARJ(169)=MIN(1D0,1.2D0*PARJ(169))
          Q2=PARJ(169)*ECM**2
          ALSPI=(3D0/4D0)*CF*PYALPS(Q2)/PARU(1)
          GOTO 100
        ENDIF
 
C...If too high cross-section, use harder cuts, or fail.
        IF(PARJ(152)+PARJ(153)+PARJ(154).GE.1) THEN
          IF(MSTJ(110).EQ.2.AND.CUT.GT.0.0499D0.AND.MSTJ(111).EQ.1.AND.
     &    PARJ(169).LT.0.99D0) THEN
            PARJ(169)=MIN(1D0,1.2D0*PARJ(169))
            Q2=PARJ(169)*ECM**2
            ALSPI=(3D0/4D0)*CF*PYALPS(Q2)/PARU(1)
            GOTO 100
          ELSEIF(MSTJ(110).EQ.2.AND.CUT.GT.0.0499D0) THEN
            CALL PYERRM(26,
     &      '(PYXJET:) no allowed y cut value for Zhu parametrization')
          ENDIF
          CUT=0.26D0*(4D0*CUT)**(PARJ(152)+PARJ(153)+
     &    PARJ(154))**(-1D0/3D0)
          IF(MSTJ(110).EQ.2) CUT=MAX(0.01D0,MIN(0.05D0,CUT))
          GOTO 100
        ENDIF
 
C...Scalar gluon (first order only).
      ELSE
        ALSPI=PYALPS(ECM**2)/PARU(1)
        CUT=MAX(0.001D0,PARJ(125),(PARJ(126)/ECM)**2,EXP(-3D0/ALSPI))
        PARJ(152)=0D0
        IF(CUT.LT.0.25D0) PARJ(152)=(ALSPI/3D0)*((1D0-2D0*CUT)*
     &  LOG((1D0-2D0*CUT)/CUT)+0.5D0*(9D0*CUT**2-1D0))
        PARJ(153)=0D0
        PARJ(154)=0D0
      ENDIF
 
C...Select number of jets.
      PARJ(150)=CUT
      IF(MSTJ(101).EQ.0.OR.MSTJ(101).EQ.5) THEN
        NJET=2
      ELSEIF(MSTJ(101).LE.0) THEN
        NJET=MIN(4,2-MSTJ(101))
      ELSE
        RNJ=PYR(0)
        NJET=2
        IF(PARJ(152)+PARJ(153)+PARJ(154).GT.RNJ) NJET=3
        IF(PARJ(154).GT.RNJ) NJET=4
      ENDIF
 
      RETURN
      END
