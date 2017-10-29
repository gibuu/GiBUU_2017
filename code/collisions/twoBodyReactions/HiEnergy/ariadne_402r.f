C***********************************************************************
C in this version, all occurences of PY... are replaced by RY...
C***********************************************************************

C***********************************************************************
C $Id: arinit.F,v 0.17 1992/03/11 14:14:59 lonnblad Exp $
C**********************************************************************C
C                                                                      C
C                            A R I A D N E                             C
C                                                                      C
C           A Monte Carlo program for colour dipole radiation          C
C                                                                      C
C                        Version 4 revision 02                         C
C                  Latest date of change: 30 Apr 1992                  C
C                                                                      C
C                              Author :                                C
C                                                                      C
C                           Leif Lonnblad                              C
C                                                                      C
C                Deutsches Elektronen Synchrotron - DESY               C
C               Notkestrasse 85, 2000 Hamburg 50, Germany              C
C                                                                      C
C                       tel  int+49-4089982048                         C
C                       fax  int+49-4089982777                         C
C                                                                      C
C                   E-mail lonnblad@apollo3.desy.de                    C
C                                                                      C
C                   Copyright (C) 1992 Leif Lonnblad                   C
C                                                                      C
C                Please report any errors to the author                C
C                                                                      C
C**********************************************************************C

C**********************************************************************C
C     This program must be loaded together with JETSET 73              C
C     The model is described in Nucl. Phys. B306 (1988) 746,           C
C     Z. Phys. C43 (1989) 625, and Nucl. Phys. B339 (1990) 393.        C
C**********************************************************************C

C***********************************************************************
CCPH:   Double precision (B) used.   
CCPH:.. This version of ARIADNE is to be used with FRITIOF 7.1
CCPH:.. All Fritiof-associated changes can be searched under "CCPH"
CCPH:   The role of LUJETS common block is replaced by 
CCPH:      COMMON /ARJETX/ N,K(300,5),P(300,5),V(300,5)
CCPH:.. Since Fritiof requires the strings to be treated one at a time,
CCPH:   the original approach of having entire LUJETS to be treated
CCPH:   by Ariadne is not possible to use.  Therefore this approach
CCPH:   of replacing LUJETS by ARJETX is adopted.  After emission is done,
CCPH:   the partons can be copied back onto LUJETS in Fritiof.
CCPH:     
C***********************************************************************

      SUBROUTINE ARINIT(MODE)

C...ARiadne subroutine INITialize

C...Initializes Ariadne to run with other (Lund) MC programs


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)

      COMMON /ARPART/ BP(MAXPAR,5),IFL(MAXPAR),IEX(MAXPAR),QQ(MAXPAR),
     $                IDI(MAXPAR),IDO(MAXPAR),INO(MAXPAR),IPART
      SAVE /ARPART/

      COMMON /ARDIPS/ BX1(MAXDIP),BX3(MAXDIP),PT2IN(MAXDIP),
     $                SDIP(MAXDIP),IP1(MAXDIP),IP3(MAXDIP),
     $                AEX1(MAXDIP),AEX3(MAXDIP),QDONE(MAXDIP),
     $                QEM(MAXDIP),IRAD(MAXDIP),ISTR(MAXDIP),IDIPS
      SAVE /ARDIPS/

      COMMON /ARSTRS/ IPF(MAXSTR),IPL(MAXSTR),IFLOW(MAXSTR),
     $                PT2LST,IMF,IML,IO,QDUMP,ISTRS
      SAVE /ARSTRS/

      COMMON /ARDAT1/ PARA(40),MSTA(40)
      SAVE /ARDAT1/

      COMMON /ARDAT2/ PQMAS(10)
      SAVE /ARDAT2/

      COMMON /ARDAT3/ IWRN(40)
      SAVE /ARDAT3/

      COMMON /ARJETX/ N,K(300,5),P(300,5),V(300,5)
      SAVE /ARJETX/

      COMMON /LUDAT1/ MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /LUDAT1/

      COMMON /LUDAT2/ KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      SAVE /LUDAT2/

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,XQ2,U
      SAVE /LEPTOU/

      COMMON /RYPARS/ MSTP(200),PARP(200),MSTI(200),PARI(200)
      SAVE /RYPARS/

      COMMON /RYSUBS/ MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200)
      SAVE /RYSUBS/
      CHARACTER MODE*(*)


C...Set output files if not already done
      IF(MSTA(7).LT.0) MSTA(7)=MSTU(11)
      IF(MSTA(8).LT.0) MSTA(8)=MSTU(11)

C...Write out header
      WRITE(MSTA(7),1000)
      MSTA(2)=1

C...If Ariadne mode, do nothing special
      IF(MODE(1:7).EQ.'ARIADNE') THEN
        MSTA(1)=0

C...If JETSET mode, switch off cascade and fragmentation in JETSET
      ELSEIF(MODE(1:6).EQ.'JETSET') THEN
        MSTA(1)=1
        MSTA(5)=MIN(MAX(MSTJ(105),0),1)
        MSTJ(101)=5
        MSTJ(41)=0
        MSTJ(105)=0
        WRITE(MSTA(7),1010)

C...If PYTHIA mode, switch off cascades and fragmentation. Check that 
C...Ariadne can handle selected processes
      ELSEIF(MODE(1:6).EQ.'PYTHIA') THEN

        MSTA(1)=2
        WRITE(MSTA(7),1020)
        MSTA(5)=MIN(MAX(MSTP(111),0),1)
        MSTP(61)=0
        MSTP(71)=0
        MSTP(111)=0

C...If LEPTO mode, switch off cascades and fragmentation.
      ELSEIF(MODE(1:5).EQ.'LEPTO') THEN
        MSTA(1)=3
        WRITE(MSTA(7),1030)
        LST(8)=9
        MSTA(5)=MIN(MAX(LST(7),0),1)
        LST(7)=0
      ENDIF

C...Set quark masses
      IF(MSTA(24).GT.0) THEN
        DO 100 I=1,8
          PQMAS(I)=PMAS(I,1)
 100    CONTINUE
      ENDIF

      IF(MSTA(24).GE.2) THEN
        DO 110 I=1,5
          PQMAS(I)=PARF(100+I)
 110    CONTINUE
      ENDIF

      IF(MSTA(3).EQ.1) CALL ARTUNE('DELPHI')

 1000 FORMAT(/,14X,
     $     'The Lund Monte Carlo - Ariadne version 4 revision 02',/,
     $     23X,'Latest date of change: 30 Apr 1992')
 1010 FORMAT(18X,'Initialization done for running with JETSET')
 1020 FORMAT(18X,'Initialization done for running with PYTHIA')
 1030 FORMAT(18X,'Initialization done for running with LEPTO')


      RETURN

C**** END OF ARINIT ****************************************************
      END
C***********************************************************************

      BLOCK DATA ARDATA

C...ARiadne block DATA statements

C...Initialization of the common blocks used in Ariadne


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)

      COMMON /ARDAT1/ PARA(40),MSTA(40)
      SAVE /ARDAT1/

      COMMON /ARDAT2/ PQMAS(10)
      SAVE /ARDAT2/

      COMMON /ARDAT3/ IWRN(40)
      SAVE /ARDAT3/

C...Breif explanation of parameters and switches:
C...
C...
C...Parameters:
C...
C...PARA(1) (D=0.200) lambda_QCD
C...PARA(2) (D=0.200) Constant alpha_QCD for MSTA(12)=0
C...PARA(3) (D=1.000) Cutoff in invariant p_t for QCD emission
C...PARA(4) (D=1/137) Constant alpha_EM
C...PARA(5) (D=1.000) Cutoff in invariant p_t for EM emission
C...PARA(6) (D=-1.00) Maximum allowed invariant p_t (if >0)
C...PARA(7) (D=0.000) Maximum invariant mass (if >0)
C...PARA(8-9) not used
C...PARA(10)(D=1.000) Power in soft suppression (dimnsionality of
C...                  the extended source)
C...PARA(11)(D=0.938) Soft suppression parameter for code 1
C...PARA(12)(D=0.938) Soft suppression parameter for code 2
C...PARA(13)(D=0.938) Soft suppression parameter for code 3
C...PARA(14-19) not used
C...PARA(20)(D=1.000) Minimum p_t^2/Q^2 of q-qbar pair in boson-gluon
C...                  fusion to be allowed to be treated as such. Else
C...                  treated as sea-quark interaction. Only used when
C...                  running with LEPTO version 6.0 or higher
C...PARA(21-30) not used
C...PARA(31)(D=1.000) Maximum invariant p_t^2 for clustering three jets
C...                  into two in ARCLUS
C...PARA(32-38) not used
C...PARA(39)(D=0.001) Tolerance factor for momentum conservation
C...PARA(40)(D=1E32)  Maximum allowed floating point number ("minimum"
C...                  is 1/PARA(40)
C...
C...Switches:
C...
C...MSTA(1) (R)       Ariadne mode (set by ARINIT) for treatment of
C...                  incomming events.
C...         0 =>      No special treatment
C...         1 =>      as if produced by JETSET
C...         2 =>      as if produced by PYTHIA
C...         3 =>      as if produced by LEPTO
C...MSTA(2) (I)       Initialization done and headers written
C...MSTA(3) (D=0)     Setting of parameters in  Ariadne, JETSET, 
C...                  PYTHIA and LEPTO to suitable values.
C...         0 =>      off
C...         1 =>      on
C...MSTA(4) (I)       Number of calls made to AREXEC
C...MSTA(5) (D=0)     Perform fragmentation at the end of AREXEC
C...         0 =>      off
C...         1 =>      on
C...                  When running with JETSET, PYTHIA or LEPTO this
C...                  switch is set to the value of the corresponding
C...                  switch in these programs.
C...MSTA(6) (D=-1)    Maximum number of emission (per string) in a
C...                  AREXEC call (if <0 - no limit)
C...MSTA(7) (D=6)     File number for output (stdout) from Ariadne
C...                  set to MSTU(11) by ARINIT
C...MSTA(8) (D=6)     File number for error messages (stdout) from
C...                  Ariadne
C...MSTA(9) (D=1)     Debug mode
C...         0 =>      debug off
C...         1 =>      check that energy and momentum is conserved after
C...                   each call to AREXEC produce. Warns if change
C...                   in momentum is larger a factor PARA(39)
C...         2 =>      as for 1 but check every emission
C...         3 =>      as for 2 but dump string to /LUJETS/ after each 
C...                   emission
C...MSTA(10)(D=5)     Maximum number of warnings (of each kind) issued
C...                  by Ariadne
C...MSTA(11)(D=0)     Phase space restrictions. The maximum p_t of an 
C...                  emission is set to the p_t of the last emission
C...                  (otherwise no restrictions) for:
C...                    gluon  q-qbar  photon  emissions
C...         0 =>        yes     yes     yes
C...         1 =>        yes     yes     no
C...         2 =>        yes     no      yes
C...         3 =>        yes     no      no
C...         4 =>        no      no      no
C...MSTA(12)(D=1)     Running alpha_QCD
C...         0 =>      off
C...         1 =>      on
C...MSTA(13) (R)      Error experienced by Ariadne in last call to 
C...                  AREXECc. Reset to 0 at each call to AREXEC
C...MSTA(14)(D=1)     The maximum allowed p_t is set to the minimum
C...                  invariant p_t of all gluons in an incomming
C...                  string
C...         0 =>      off
C...         1 =>      on
C...MSTA(15)(D=5)     Number of flavours allowed in q-qbar emission
C...MSTA(16)(D=2)     Recoil treatment
C...         0 =>      minimize p_t1^2 + p_t3^2
C...         1 =>      as for 0 but pointlike string ends takes
C...                   all recoil
C...         2=>       as for 0 but also extended string ends which
C...                   have a>0 takes full recoil
C...MSTA(17)(D=2)     Recoil treatment for extended dipoles
C...         0 =>      no special treatment (but cf. MSTA(16))
C...         1 =>      emit recoil gluon (except if pointlike quark
C...                   in other dipole end for MSTA(16)=1)
C...         2 =>      emit recoilgluon according to new strategy
C...MSTA(18)(D=1)     P_t ordering of recoil gluons
C...         0 =>      off
C...         1 =>      on
C...MSTA(19)(D=1)     Correct or quick treatment of emissions from
C...                  heavy quarks
C...         0 =>      quick
C...         1 =>      correct
C...MSTA(20)(D=0)     Final state photon radiation
C...         0 =>      off
C...         1 =>      on
C...         2 =>      on but turned off at the first occurence of
C...                   q-qbar emission in a string.
C...MSTA(21)(D=0)     Photon radiation when run with PYTHIA or LEPTO
C...         0 =>      off
C...         1 =>      on
C...MSTA(22)(D=0)     Transfer of recoils in Drell-Yan processes
C...         0 =>      off
C...         1 =>      on
C...MSTA(23)(I)       Line number of particle to transfer recoil to 
C...                  for MSTA(22) > 0
C...MSTA(24)(D=2)     Quark masses to be used in q-qbar emissions
C...         0 =>      as specified in PMAS(1-8) in /ARDAT2/
C...         1 =>      "bare" quark masses as specified in PMAS(1-8)
C...                   in /LUDAT2/
C...         2 =>      "constituent" quark masses as specified in 
C...                   PARF(101-108) in /LUDAT2/
C...MSTA(25-29) not used
C...MSTA(30)(D=1)    various options for running with Lepto
C...         0 =>      Stuck quark point like, remnant extended with PARA(11)
C...         1 =>      as 0 but remnant extended with PARA(11)/(1-x)
C...         2 =>      as 1 bur struck quark extended with Q
C...MSTA(31)(D=1)    mass of extended partons
C...         0 =>      set to zero for backward compatibility
C...         1 =>      keeps value given
C...MSTA(32-40) not used
C...
C...End of description

      DATA PARA/0.2,0.2,1.0,0.007297353,1.0,-1.0,0.0,0.0,0.0,1.0,
     $          0.938,0.938,0.938,6*0.0,1.0,
     $          10*0.0,
     $          1.0,7*0.0,0.001,1.0E32/
      DATA MSTA/0,0,0,0,0,-1,6,6,1,5,
     $          0,1,0,1,5,2,2,1,1,0,
     $          0,0,0,2,5*0,1,
     $          1,9*0/
      DATA PQMAS/10*0.0/
      DATA IWRN/40*0/


C**** END OF ARDATA ****************************************************
      END
C***********************************************************************
C $Id: araddg.F,v 0.5 1992/01/31 16:14:59 lonnblad Exp $

      SUBROUTINE ARADDG(ID)

C...ARiadne subroutine ADD Gluon

C...Adds a gluon entry between the partons in dipole ID thus creating 
C...a new dipole


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)

      COMMON /ARPART/ BP(MAXPAR,5),IFL(MAXPAR),IEX(MAXPAR),QQ(MAXPAR),
     $                IDI(MAXPAR),IDO(MAXPAR),INO(MAXPAR),IPART
      SAVE /ARPART/

      COMMON /ARDIPS/ BX1(MAXDIP),BX3(MAXDIP),PT2IN(MAXDIP),
     $                SDIP(MAXDIP),IP1(MAXDIP),IP3(MAXDIP),
     $                AEX1(MAXDIP),AEX3(MAXDIP),QDONE(MAXDIP),
     $                QEM(MAXDIP),IRAD(MAXDIP),ISTR(MAXDIP),IDIPS
      SAVE /ARDIPS/

      COMMON /ARSTRS/ IPF(MAXSTR),IPL(MAXSTR),IFLOW(MAXSTR),
     $                PT2LST,IMF,IML,IO,QDUMP,ISTRS
      SAVE /ARSTRS/

      INXT(I)=IDO(IP3(I))
      IPRV(I)=IDI(IP1(I))


C...Allocate new gluon and new dipole at postitons IPART+1 and IDIPS+1
C...if there is space left.
      IPART=IPART+1
      IDIPS=IDIPS+1
      IF(IPART.GE.MAXPAR-1) CALL ARERRM('ARADDG',6,0)
      IF(IDIPS.GE.MAXDIP-1) CALL ARERRM('ARADDG',7,0)

C...Set properties of new gluon
      DO 100 I=1,5
        BP(IPART,I)=0.0
 100  CONTINUE
      IFL(IPART)=21
      IEX(IPART)=0
      QQ(IPART)=.FALSE.
      IDI(IPART)=ID
      IDO(IPART)=IDIPS

C...Set properties of new dipole
      IP1(IDIPS)=IPART
      IP3(IDIPS)=IP3(ID)
      QDONE(IDIPS)=.FALSE.
      QEM(IDIPS)=.FALSE.
      ISTR(IDIPS)=ISTR(ID)

C...Fix pointers for old dipole
      IP3(ID)=IPART
      IDI(IP3(IDIPS))=IDIPS
      IF(IPRV(ID).NE.0) QDONE(IPRV(ID))=.FALSE.
      QDONE(ID)=.FALSE.
      IF(INXT(IDIPS).NE.0) QDONE(INXT(IDIPS))=.FALSE.

      RETURN

C**** END OF ARADDG ****************************************************
      END
C***********************************************************************
C $Id: arangl.F,v 0.2 1991/09/26 12:41:40 lonnblad Exp $

      REAL FUNCTION ARANGL(I1,I2)

C...ARiadne function ANGLe

C...Returns the angle between paron I1 and I2


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)

      COMMON /ARPART/ BP(MAXPAR,5),IFL(MAXPAR),IEX(MAXPAR),QQ(MAXPAR),
     $                IDI(MAXPAR),IDO(MAXPAR),INO(MAXPAR),IPART
      SAVE /ARPART/


      ARANGL=BP(I1,1)*BP(I2,1)+BP(I1,2)*BP(I2,2)+BP(I1,3)*BP(I2,3)
      BP1=SQRT(BP(I1,1)**2+BP(I1,2)**2+BP(I1,3)**2)
      BP2=SQRT(BP(I2,1)**2+BP(I2,2)**2+BP(I2,3)**2)
      ARANGL=ACOS(ARANGL/(BP1*BP2))

      RETURN

C**** END OF ARANGL ****************************************************
      END
C***********************************************************************
C $Id: arbocm.F,v 0.5 1991/11/06 13:46:04 lonnblad Exp $

      SUBROUTINE ARBOCM(ID)

C...ARiadne subroutine BOost to Center of Mass

C...Boost the partons in dipole ID to the CMS of the dipole


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)

      COMMON /ARPART/ BP(MAXPAR,5),IFL(MAXPAR),IEX(MAXPAR),QQ(MAXPAR),
     $                IDI(MAXPAR),IDO(MAXPAR),INO(MAXPAR),IPART
      SAVE /ARPART/

      COMMON /ARDIPS/ BX1(MAXDIP),BX3(MAXDIP),PT2IN(MAXDIP),
     $                SDIP(MAXDIP),IP1(MAXDIP),IP3(MAXDIP),
     $                AEX1(MAXDIP),AEX3(MAXDIP),QDONE(MAXDIP),
     $                QEM(MAXDIP),IRAD(MAXDIP),ISTR(MAXDIP),IDIPS
      SAVE /ARDIPS/

      COMMON /ARSTRS/ IPF(MAXSTR),IPL(MAXSTR),IFLOW(MAXSTR),
     $                PT2LST,IMF,IML,IO,QDUMP,ISTRS
      SAVE /ARSTRS/

      COMMON /ARINT2/ DBEX,DBEY,DBEZ,PHI,THE
      SAVE /ARINT2/


C...Calculate boostvector and boost
      I1=IP1(ID)
      I3=IP3(ID)
      DPE1=BP(I1,4)
      DPE3=BP(I3,4)
      DPE=DPE1+DPE3
      DPX1=BP(I1,1)
      DPX3=BP(I3,1)
      DBEX=(DPX1+DPX3)/DPE
      DPY1=BP(I1,2)
      DPY3=BP(I3,2)
      DBEY=(DPY1+DPY3)/DPE
      DPZ1=BP(I1,3)
      DPZ3=BP(I3,3)
      DBEZ=(DPZ1+DPZ3)/DPE
      CALL AROBO2(0.0,0.0,-DBEX,-DBEY,-DBEZ,I1,I3)

C...Calculate rotation angles but no need for rotation yet
      PX=BP(I1,1)
      PY=BP(I1,2)
      PZ=BP(I1,3)
      PHI=ULANGL(PX,PY)
      THE=ULANGL(PZ,SQRT(PX**2+PY**2))

      RETURN

C**** END OF ARBOCM ****************************************************
      END
C***********************************************************************
C $Id: arcasc.F,v 0.7 1992/04/28 13:16:33 lonnblad Exp $

      SUBROUTINE ARCASC

C...ARiadne subroutine perform dipole CASCade

C...Performs a colour dipole cascade on string put in the ariadne
C...event record.


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)

      COMMON /ARPART/ BP(MAXPAR,5),IFL(MAXPAR),IEX(MAXPAR),QQ(MAXPAR),
     $                IDI(MAXPAR),IDO(MAXPAR),INO(MAXPAR),IPART
      SAVE /ARPART/

      COMMON /ARDIPS/ BX1(MAXDIP),BX3(MAXDIP),PT2IN(MAXDIP),
     $                SDIP(MAXDIP),IP1(MAXDIP),IP3(MAXDIP),
     $                AEX1(MAXDIP),AEX3(MAXDIP),QDONE(MAXDIP),
     $                QEM(MAXDIP),IRAD(MAXDIP),ISTR(MAXDIP),IDIPS
      SAVE /ARDIPS/

      COMMON /ARSTRS/ IPF(MAXSTR),IPL(MAXSTR),IFLOW(MAXSTR),
     $                PT2LST,IMF,IML,IO,QDUMP,ISTRS
      SAVE /ARSTRS/

      COMMON /ARDAT1/ PARA(40),MSTA(40)
      SAVE /ARDAT1/


C...Calculate total momentum of strings for debugging
      IF(MSTA(9).GT.0) CALL ARCHEM(1)

C...Reset counter
      IO=0

C...Loop over all dipole to find largest possible p_t^2
 100  ISEL=0
      PT2MAX=0.0
      DO 110 I=1,IDIPS
        PT2I=ARGPT2(I)
        IF(PT2I.GT.PT2MAX) THEN
          PT2MAX=PT2I
          ISEL=I
        ENDIF
 110  CONTINUE

C...Check that largest p_t^2 is above cuts.
      IF(ISEL.GT.0) THEN
        IF((QEM(ISEL).AND.PT2MAX.LE.PARA(5)**2).OR.
     $     ((.NOT.QEM(ISEL)).AND.PT2MAX.LE.PARA(3)**2)) ISEL=0
      ENDIF

      IF(MSTA(6).GE.0.AND.IO.GE.MSTA(6)) ISEL=0

C...Exit if below cuts or limit of number of emissions is reached
      IF(ISEL.EQ.0) THEN
        CALL ARDUMP
        IF(MSTA(9).GT.0) CALL ARCHEM(0)
        RETURN
      ENDIF

C...Perform the emission
      IO=IO+1
      PT2LST=PT2MAX
      CALL AREMIT(ISEL)
      QDUMP=.FALSE.

C...Check total momentum and dump according to debug mode
      IF(MSTA(9).GT.2) CALL ARDUMP
      IF(MSTA(9).GT.1) CALL ARCHEM(0)
      GOTO 100

C**** END OF ARCASC ****************************************************
      END
C***********************************************************************
C $Id: archem.F,v 0.6 1992/01/31 16:14:59 lonnblad Exp $

      SUBROUTINE ARCHEM(IMOD)

C...ARiadne subroutine CHEck Momentum conservation

C...Checks that momentum is conserved in ariadne


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)

      COMMON /ARPART/ BP(MAXPAR,5),IFL(MAXPAR),IEX(MAXPAR),QQ(MAXPAR),
     $                IDI(MAXPAR),IDO(MAXPAR),INO(MAXPAR),IPART
      SAVE /ARPART/

      COMMON /ARINT3/ DPTOT(5)
      SAVE /ARINT3/

      COMMON /ARDAT1/ PARA(40),MSTA(40)
      SAVE /ARDAT1/

      COMMON /ARJETX/ N,K(300,5),P(300,5),V(300,5)
      SAVE /ARJETX/

      DIMENSION DTOT(5)


C...Reset momentum counter. Include Drell-Yan produced particle
C...if present and check its momentum consistency.
      IF(MSTA(23).GT.0) THEN
        I=MSTA(23)
        DO 100 J=1,4
          DTOT(J)=P(I,J)
 100    CONTINUE
        IF(ABS(P(I,4)**2-P(I,3)**2-P(I,2)**2-P(I,1)**2-P(I,5)**2)
     $         .GT.PARA(39)*P(I,4)**2) CALL ARERRM('ARCHEM',10,I)
      ELSE
        DO 110 J=1,4
          DTOT(J)=0.0D0
 110    CONTINUE
      ENDIF

C...Sum all partons momentum and check their momentum concistency.
      DO 120 I=1,IPART
        DO 130 J=1,4
          DTOT(J)=DTOT(J)+BP(I,J)
 130    CONTINUE
        IF(ABS(BP(I,4)**2-BP(I,3)**2-BP(I,2)**2-BP(I,1)**2-BP(I,5)**2)
     $       .GT.PARA(39)*BP(I,4)**2.AND.MSTA(9).GE.2)
     $       CALL ARERRM('ARCHEM',10,I+N)
 120  CONTINUE
      DTOT(5)=DSQRT(DTOT(4)**2-DTOT(3)**2-DTOT(2)**2-DTOT(1)**2)

C...If IMOD=1 save total momentum for later use
      IF(IMOD.EQ.1) THEN
        DO 200 J=1,5
          DPTOT(J)=DTOT(J)
 200    CONTINUE
        RETURN
      ENDIF

C...IF IMOD=1 compare total momentum with old one
      DO 300 J=1,5
        IF(ABS(DTOT(J)-DPTOT(J)).GT.DPTOT(5)*PARA(39))
     $       CALL ARERRM('ARCHEM',9,0)
 300  CONTINUE

      RETURN

C**** END OF ARCHEM ****************************************************
      END
C***********************************************************************
C $Id: arclus.F,v 0.7 1992/03/13 13:36:14 lonnblad Exp $

      SUBROUTINE ARCLUS(NJET)

C...ARiadne subroutine jet CLUStering

C...Clusters particle in the /LUJETS/ common block into jets according
C...the dipole clustering algorithm.


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)

      COMMON /ARDAT1/ PARA(40),MSTA(40)
      SAVE /ARDAT1/

      COMMON /ARJETX/ N,K(300,5),P(300,5),V(300,5)
      SAVE /ARJETX/

      COMMON /LUDAT1/ MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /LUDAT1/


C...Reset error flag.
      MSTA(13)=0

C...Copy all particle to be considered as jet-initiators to the end
C...of the event record.
      IF(MSTU(48).EQ.0) CALL ARCOPJ

C...The total number of jetinitiators = current number of jets.
      NJET=MSTU(3)
      I1=0
      I3=0

C...Loop over all possible three-jets to find the three jets with
C...smallest invariant p_t^2
100   IF(NJET.LE.MAX(MSTU(47),2)) THEN
        CALL ARORDJ
        RETURN
      ENDIF

      J1=0
      J2=0
      J3=0
      PT2MIN=PARA(31)

      DO 110 I2=N+1,N+MSTU(3)
        IF(K(I2,5).LT.0) GOTO 110
        CALL ARUPDJ(I2,I1,I3)
        IF(V(I2,5).LT.PT2MIN) THEN
          J1=K(I2,3)
          J2=I2
          J3=K(I2,4)
          PT2MIN=V(I2,5)
        ENDIF
110   CONTINUE

C...Exit if smallest p_t^2 is above cutoff
      IF(J1.EQ.0) THEN
        CALL ARORDJ
        RETURN
      ENDIF

C...Else join the three jets into two and redo the procedure.
      CALL ARJOIN(J1,J2,J3)
      K(J2,5)=-1
      I1=J1
      I3=J3
      NJET=NJET-1

      GOTO 100

C**** END OF ARCLUS ****************************************************
      END
C***********************************************************************
C $Id: arcopa.F,v 0.4 1992/03/13 09:02:26 lonnblad Exp $

      SUBROUTINE ARCOPA(IJ,IP,ITYP)

C...ARiadne subroutine COpy PArton

C...Copies a parton from position IJ in /LUJETS/ common block to
C...Position IP in /ARPART/ common block.


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)

      COMMON /ARPART/ BP(MAXPAR,5),IFL(MAXPAR),IEX(MAXPAR),QQ(MAXPAR),
     $                IDI(MAXPAR),IDO(MAXPAR),INO(MAXPAR),IPART
      SAVE /ARPART/

      COMMON /ARJETX/ N,K(300,5),P(300,5),V(300,5)
      SAVE /ARJETX/


      DO 100 I=1,5
        BP(IP,I)=P(IJ,I)
 100  CONTINUE

      IFL(IP)=K(IJ,2)
      IEX(IP)=MOD(K(IJ,4),10)
      QQ(IP)=(ITYP.NE.2)
      INO(IP)=0
      RETURN

C**** END OF ARCOPA ****************************************************
      END
C***********************************************************************
C $Id: arcopj.F,v 0.6 1992/02/07 16:03:59 lonnblad Exp $

      SUBROUTINE ARCOPJ

C...ARiadne subroutine COPy Jet

C...Copies particles into jet initiators in /LUCLUS/


      COMMON /ARJETX/ N,K(300,5),P(300,5),V(300,5)
      SAVE /ARJETX/

      COMMON /LUDAT1/ MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /LUDAT1/

      COMMON /LUDAT2/ KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      SAVE /LUDAT2/


C...Reset jet counter
      MSTU(3)=0

C...Loop over all particles in the event record
      DO 100 I=1,N

C...Disregard all decayed particles and unknown entries
        IF(K(I,1).LE.0.OR.K(I,1).GE.10) GOTO 100

C...Disregard neutrinos and neutral particles according to MSTU(41)
        IF(MSTU(41).GE.2) THEN
          KC=LUCOMP(K(I,2))
          IF(KC.EQ.0.OR.KC.EQ.12.OR.KC.EQ.14.OR.KC.EQ.16.OR.KC.EQ.18) 
     $      GOTO 100
          IF(MSTU(41).GE.3.AND.KCHG(KC,2).EQ.0.AND.LUCHGE(K(I,2)).EQ.0)
     $      GOTO 100
        ENDIF

        IF(N+MSTU(3)+1.GT.MSTU(4)) THEN
          CALL LUERRM(11,'(ARCLUS:) no more memory left in LUJETS')
          MSTU(3)=-1
          RETURN
        ENDIF

C...Tag found jet-initiator
        MSTU(3)=MSTU(3)+1
        IJ=N+MSTU(3)
        DO 200 J=1,5
          P(IJ,J)=P(I,J)
 200    CONTINUE
        K(IJ,1)=31
        K(IJ,2)=97
        K(IJ,3)=0
        K(IJ,4)=0
        K(IJ,5)=0
        V(IJ,1)=P(I,4)**2-P(I,3)**2-P(I,2)**2-P(I,1)*2
        V(IJ,5)=0

 100  CONTINUE

      RETURN

C**** END OF ARCOPJ ****************************************************
      END
C***********************************************************************
C $Id: arcrdi.F,v 0.2 1991/09/26 12:42:05 lonnblad Exp $

      SUBROUTINE ARCRDI(ID,IPA1,IPA3,IS,QED)

C...ARiadne subroutine CReate DIpole

C...Creates a dipole from partons IPA1 and IPA3


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)

      COMMON /ARPART/ BP(MAXPAR,5),IFL(MAXPAR),IEX(MAXPAR),QQ(MAXPAR),
     $                IDI(MAXPAR),IDO(MAXPAR),INO(MAXPAR),IPART
      SAVE /ARPART/

      COMMON /ARDIPS/ BX1(MAXDIP),BX3(MAXDIP),PT2IN(MAXDIP),
     $                SDIP(MAXDIP),IP1(MAXDIP),IP3(MAXDIP),
     $                AEX1(MAXDIP),AEX3(MAXDIP),QDONE(MAXDIP),
     $                QEM(MAXDIP),IRAD(MAXDIP),ISTR(MAXDIP),IDIPS
      SAVE /ARDIPS/

      COMMON /ARSTRS/ IPF(MAXSTR),IPL(MAXSTR),IFLOW(MAXSTR),
     $                PT2LST,IMF,IML,IO,QDUMP,ISTRS
      SAVE /ARSTRS/


      IDO(IPA1)=ID
      IDI(IPA3)=ID
      IP1(ID)=IPA1
      IP3(ID)=IPA3
      ISTR(ID)=IS
      QDONE(ID)=.FALSE.
      QEM(ID)=QED

      RETURN

C**** END OF ARCRDI ****************************************************
      END
C***********************************************************************
C $Id: ardump.F,v 0.5 1992/02/07 16:01:36 lonnblad Exp $

      SUBROUTINE ARDUMP

C...ARiadne subroutine DUMP 

C...Dumps the entries in /ARPART/ into /LUJETS/


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)

      COMMON /ARPART/ BP(MAXPAR,5),IFL(MAXPAR),IEX(MAXPAR),QQ(MAXPAR),
     $                IDI(MAXPAR),IDO(MAXPAR),INO(MAXPAR),IPART
      SAVE /ARPART/

      COMMON /ARDIPS/ BX1(MAXDIP),BX3(MAXDIP),PT2IN(MAXDIP),
     $                SDIP(MAXDIP),IP1(MAXDIP),IP3(MAXDIP),
     $                AEX1(MAXDIP),AEX3(MAXDIP),QDONE(MAXDIP),
     $                QEM(MAXDIP),IRAD(MAXDIP),ISTR(MAXDIP),IDIPS
      SAVE /ARDIPS/

      COMMON /ARSTRS/ IPF(MAXSTR),IPL(MAXSTR),IFLOW(MAXSTR),
     $                PT2LST,IMF,IML,IO,QDUMP,ISTRS
      SAVE /ARSTRS/

      COMMON /ARJETX/ N,K(300,5),P(300,5),V(300,5)
      SAVE /ARJETX/

CCPH:..............................................................
      COMMON /ARDAT1/ PARA(40),MSTA(40)
      SAVE /ARDAT1/
CCPH:^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      INXT(I)=IP3(IDO(I))

C...Tag particles in old string with pointers to cascaded string 
      DO 100 I=IMF,IML
        K(I,1)=K(I,1)+10
        K(I,4)=N+1
        K(I,5)=N+IPART
 100  CONTINUE

C...Loop over all strings in dipole record
      DO 200 IS=1,ISTRS

C...Loop over all particles in each string
        I=IPF(IS)
 210    N=N+1

CCPH:..............................................................
        IF(N.GT.300) THEN
        WRITE(MSTA(8),*) '**',N,'  Extend ARJETX in Ariadne **'
        STOP
        ENDIF 
CCPH:^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        
        DO 220 J=1,5
          P(N,J)=BP(I,J)
          V(N,J)=V(IMF,J)
 220    CONTINUE
        K(N,2)=IFL(I)
        K(N,3)=IMF
        K(N,4)=IEX(I)
        K(N,5)=INO(I)
        IF(I.EQ.IPL(IS)) THEN
          K(N,1)=1
        ELSE
          K(N,1)=2
          I=INXT(I)
          GOTO 210
        ENDIF
 200  CONTINUE

C...Set pointers to cascaded string
      IMF=N+1-IPART
      IML=N
      QDUMP=.TRUE.

      RETURN

C**** END OF ARDUMP ****************************************************
      END
C***********************************************************************
C $Id: arduph.F,v 0.4 1992/02/07 16:00:36 lonnblad Exp $

      SUBROUTINE ARDUPH

C...ARiadne subroutine DUmp PHoton

C...Moves photon emitted by Ariadne to /LUJETS/


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)

      COMMON /ARPART/ BP(MAXPAR,5),IFL(MAXPAR),IEX(MAXPAR),QQ(MAXPAR),
     $                IDI(MAXPAR),IDO(MAXPAR),INO(MAXPAR),IPART
      SAVE /ARPART/

      COMMON /ARDIPS/ BX1(MAXDIP),BX3(MAXDIP),PT2IN(MAXDIP),
     $                SDIP(MAXDIP),IP1(MAXDIP),IP3(MAXDIP),
     $                AEX1(MAXDIP),AEX3(MAXDIP),QDONE(MAXDIP),
     $                QEM(MAXDIP),IRAD(MAXDIP),ISTR(MAXDIP),IDIPS
      SAVE /ARDIPS/

      COMMON /ARSTRS/ IPF(MAXSTR),IPL(MAXSTR),IFLOW(MAXSTR),
     $                PT2LST,IMF,IML,IO,QDUMP,ISTRS
      SAVE /ARSTRS/

      COMMON /ARINT3/ DPTOT(5)
      SAVE /ARINT3/

      COMMON /ARJETX/ N,K(300,5),P(300,5),V(300,5)
      SAVE /ARJETX/


      N=N+1
      DO 100 I=1,5
        P(N,I)=BP(IPART+1,I)
        DPTOT(I)=DPTOT(I)-BP(IPART+1,I)
        V(N,I)=V(IMF,I)
 100  CONTINUE

      DPTOT(5)=DSQRT(DPTOT(4)**2-DPTOT(3)**2-DPTOT(2)**2-DPTOT(1)**2)

      K(N,1)=1
      K(N,2)=22
      K(N,3)=IMF
      K(N,4)=0
      K(N,5)=IO

      RETURN

C**** END OF ARDUPH ****************************************************
      END
C***********************************************************************
C $Id: ardyre.F,v 0.7 1992/01/31 16:14:59 lonnblad Exp $

      SUBROUTINE ARDYRE(ID,*)

C...ARiadne subroutine Drell-Yan REcoil treatment

C...Transfers the recoil of an emission to a Drell-Yan produced
C...particle if the emission and the particle are in the same
C...phase space region.


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)

      double precision PYR

      COMMON /ARPART/ BP(MAXPAR,5),IFL(MAXPAR),IEX(MAXPAR),QQ(MAXPAR),
     $                IDI(MAXPAR),IDO(MAXPAR),INO(MAXPAR),IPART
      SAVE /ARPART/

      COMMON /ARDIPS/ BX1(MAXDIP),BX3(MAXDIP),PT2IN(MAXDIP),
     $                SDIP(MAXDIP),IP1(MAXDIP),IP3(MAXDIP),
     $                AEX1(MAXDIP),AEX3(MAXDIP),QDONE(MAXDIP),
     $                QEM(MAXDIP),IRAD(MAXDIP),ISTR(MAXDIP),IDIPS
      SAVE /ARDIPS/

      COMMON /ARSTRS/ IPF(MAXSTR),IPL(MAXSTR),IFLOW(MAXSTR),
     $                PT2LST,IMF,IML,IO,QDUMP,ISTRS
      SAVE /ARSTRS/

      COMMON /ARDAT1/ PARA(40),MSTA(40)
      SAVE /ARDAT1/

      COMMON /ARINT2/ DBEX,DBEY,DBEZ,PHI,THE
      SAVE /ARINT2/

      COMMON /ARJETX/ N,K(300,5),P(300,5),V(300,5)
      SAVE /ARJETX/

      COMMON /LUDAT1/ MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /LUDAT1/

      


C...Locate Drell-Yan produced particle (IDY) and boost it to CMS 
C...of dipole
      IDY=MSTA(23)
      CALL LUDBRB(IDY,IDY,0.0,0.0,-DBEX,-DBEY,-DBEZ)
      CALL LUDBRB(IDY,IDY,0.0,-PHI,0.0D0,0.0D0,0.0D0)
      CALL LUDBRB(IDY,IDY,-THE,0.0,0.0D0,0.0D0,0.0D0)

C...Calculate p_t and y for emitted gluon and light cone momenta for
C...IDY
      PTG=SQRT(PT2IN(ID))
      ZG=SQRT((1.0-BX1(ID))/(1.0-BX3(ID)))
      BPDY=P(IDY,4)+P(IDY,3)
      BMDY=P(IDY,4)-P(IDY,3)

C...If gluon is 'outside' IDYs phase-space, exit and perform normal 
C...emission
      IF(PTG.GT.BPDY*BMDY/(BMDY*ZG+BPDY/ZG)) THEN
        CALL LUDBRB(IDY,IDY,THE,PHI,DBEX,DBEY,DBEZ)
        RETURN
      ENDIF

C...Calculate positions of particles in imitting dipole and emitted
C...gluon
      I1=IP1(ID)
      I3=IP3(ID)
      CALL ARADDG(ID)
      IG=IP3(ID)

C...Set momenta for gluon
      BPG=PTG*ZG
      BMG=PTG/ZG

      W=SQRT(SDIP(ID))
      BPTOT=W+BPDY
      BMTOT=W+BMDY

      BET=PARU(2)*PYR(IDUM)

      BP(IG,1)=PTG*SIN(BET)
      BP(IG,2)=PTG*COS(BET)
      BP(IG,3)=0.5*(BPG-BMG)
      BP(IG,4)=0.5*(BPG+BMG)
      BP(IG,5)=0.0

C...Transfer transverse recoil to IDY and set new momenta for IDY
      P(IDY,1)=P(IDY,1)-BP(IG,1)
      P(IDY,2)=P(IDY,2)-BP(IG,2)
      XMT2=(P(IDY,1)**2+P(IDY,2)**2+P(IDY,5)**2)/(BPDY*BMDY)
      BPDY=BPDY*SQRT(XMT2)
      BMDY=BMDY*SQRT(XMT2)
      P(IDY,3)=0.5*(BPDY-BMDY)
      P(IDY,4)=0.5*(BPDY+BMDY)

      BPTOT=BPTOT-BPDY-BPG
      BMTOT=BMTOT-BMDY-BMG

C...Set new momenta for particles in emitting dipole and exit if
C...the recoil transfer is not kinematically allowed
      Y1=BP(I1,5)**2
      Y3=BP(I3,5)**2

      IF(BMTOT.LT.1.0E-20) CALL ARERRM('ARDYRE',11,0)

      BB=0.5*(BPTOT*BMTOT+Y1-Y3)/BMTOT
      BA=Y1*BPTOT/BMTOT

      IF(BB**2-BA.LT.0.0) CALL ARERRM('ARDYRE',11,0)

      BP1=BB+SQRT(BB**2-BA)

      IF(BP1.LE.SQRT(Y1)) CALL ARERRM('ARDYRE',11,0)
      BM1=Y1/BP1

      BM3=BMTOT-BM1

      IF(BM3.LE.SQRT(Y3)) CALL ARERRM('ARDYRE',11,0)
      BP3=Y3/BM3

      BP(I1,1)=0.0
      BP(I1,2)=0.0
      BP(I1,3)=0.5*(BP1-BM1)
      BP(I1,4)=0.5*(BP1+BM1)

      BP(I3,1)=0.0
      BP(I3,2)=0.0
      BP(I3,3)=0.5*(BP3-BM3)
      BP(I3,4)=0.5*(BP3+BM3)

C...Boost back all particles to original system
      CALL LUDBRB(IDY,IDY,THE,PHI,DBEX,DBEY,DBEZ)
      CALL AROBO3(THE,PHI,DBEX,DBEY,DBEZ,I1,IG,I3)

      RETURN 1

C**** END OF ARDYRE ****************************************************
      END
C***********************************************************************
C $Id: aremit.F,v 0.12 1992/03/12 12:02:30 lonnblad Exp $

      SUBROUTINE AREMIT(ID)

C...ARiadne subroutine EMIT

C...Administers the an emission from dipole ID


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)

      COMMON /ARPART/ BP(MAXPAR,5),IFL(MAXPAR),IEX(MAXPAR),QQ(MAXPAR),
     $                IDI(MAXPAR),IDO(MAXPAR),INO(MAXPAR),IPART
      SAVE /ARPART/

      COMMON /ARDIPS/ BX1(MAXDIP),BX3(MAXDIP),PT2IN(MAXDIP),
     $                SDIP(MAXDIP),IP1(MAXDIP),IP3(MAXDIP),
     $                AEX1(MAXDIP),AEX3(MAXDIP),QDONE(MAXDIP),
     $                QEM(MAXDIP),IRAD(MAXDIP),ISTR(MAXDIP),IDIPS
      SAVE /ARDIPS/

      COMMON /ARSTRS/ IPF(MAXSTR),IPL(MAXSTR),IFLOW(MAXSTR),
     $                PT2LST,IMF,IML,IO,QDUMP,ISTRS
      SAVE /ARSTRS/

      COMMON /ARDAT1/ PARA(40),MSTA(40)
      SAVE /ARDAT1/

      INXT(I)=IDO(IP3(I))


C...If FSR photon emission go a head
      IF(QEM(ID)) THEN
        CALL ARRADP(ID)
        RETURN

C...If q-qbar splitting go a head
      ELSEIF(IRAD(ID).NE.0) THEN
        CALL ARRADQ(ID)
        RETURN

C...If gluon emission from point-like dipole or if no p_t-ordered
C...recoil gluon, go a head
      ELSEIF((IEX(IP1(ID)).EQ.0.AND.IEX(IP3(ID)).EQ.0)
     $             .OR.MSTA(18).EQ.0) THEN
        CALL ARRADG(ID,0,SNR,PT21,PT23)
        RETURN
      ENDIF

C...If p_t-ordered recoil gluon, first save initial configuration
C...Then perform trial emission
      CALL ARSTOR(ID,IDS,IS1,IS3)
      CALL ARRADG(ID,0,SNR,PT21,PT23)

C...If no recoil gluon was produces keep trial emission
      IF(SNR.LE.1.0) RETURN

C...If two recoil gluons, tag the smallest one for p_t-ordering
      IF(AEX1(ID).LT.1.0.AND.AEX3(ID).LT.1.0) THEN
        INEWD=3
        IF(PT23.LT.PT21) THEN
          IGR=3
          IDR=INXT(INXT(ID))
        ELSE
          IGR=1
          IDR=ID
        ENDIF

C...If only one recoil gluon, tag it for p_t-ordering
      ELSEIF(AEX1(ID).LT.1.0.AND.AEX3(ID).GE.1.0) THEN
        IGR=1
        IDR=ID
        INEWD=2
      ELSEIF(AEX1(ID).GE.1.0.AND.AEX3(ID).LT.1.0) THEN
        IGR=3
        IDR=INXT(ID)
        INEWD=2
      ENDIF

      IDT=MAXDIP-1

C...Calculate the p_t^2 of a possibly earlier emission in place
C...of the recoil gluon. If this p_t^2 is lower than that of the
C...recoil gluon it could not have been emitted earlier and hence
C...the recoil gluon from the trial emission is kept.
      IF(IGR.EQ.1) THEN
        SY1=BP(IS1,5)/SQRT(SNR)
        CALL ARGQTE(IDT,SNR,PT2IN(IDS)/SNR,QQ(IS1),.FALSE.,
     $              IEX(IS1),0,SY1,0.0)
        IF(PT2IN(IDT).LT.PT21.AND.PT21.GT.PARA(3)**2
     $       .AND.PT21.GT.PARA(10+IEX(IS1))**2) RETURN
      ELSE
        SY3=BP(IS3,5)/SQRT(SNR)
        CALL ARGQTE(IDT,SNR,PT2IN(IDS)/SNR,.FALSE.,QQ(IS3),
     $              0,IEX(IS3),0.0,SY3)
        IF(PT2IN(IDT).LT.PT23.AND.PT23.GT.PARA(3)**2
     $       .AND.PT23.GT.PARA(10+IEX(IS3))**2) RETURN
      ENDIF

C...A gluon can be emittes in place of the recoil gluon at an earlier 
C...time. Recall the initial configuration and redo the emission without
C...recoil gluon
      CALL ARRECA(ID,IDS,IS1,IS3)

      IDIPS=IDIPS-INEWD
      IPART=IPART-INEWD
      CALL ARRADG(ID,IGR,SNREF,PT21,PT23)

C...Set p_t^2 for the emission in place of the recoil gluon
      IDS=ID
      IF(IGR.EQ.3) THEN
        IDS=INXT(ID)
        IF(INEWD.EQ.3) IDS=INXT(IDS)
      ENDIF

      CALL ARSTOR(IDS,IDSS,ISS1,ISS3)
      IP1(IDSS)=ISS1
      IP3(IDSS)=ISS3
      CALL ARBOCM(IDSS)

      QDONE(IDS)=.TRUE.
      SDIP(IDS)=ARMAS2(ISS1,ISS3)
      BX1(IDS)=BX1(IDT)
      BX3(IDS)=BX3(IDT)
      AEX1(IDS)=AEX1(IDT)
      AEX3(IDS)=AEX3(IDT)
      IRAD(IDS)=IRAD(IDT)
      PT2IN(IDS)=PT2IN(IDT)

      CALL ARCHKI(IDS,IOK)
      IF(IOK.EQ.0.AND.PT2IN(IDS).GT.PARA(3)**2) THEN
        QDONE(IDS)=.FALSE.
      ENDIF

      RETURN

C**** END OF AREMIT ****************************************************
      END
C***********************************************************************
C $Id: arerrm.F,v 0.12 1992/03/13 14:21:02 lonnblad Exp $

      SUBROUTINE ARERRM(SUB,IERR,ILINE)

C...ARiadne subroutine ERRor Message

C...Writes out an error message and optionally terminates the program


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)

      COMMON /ARSTRS/ IPF(MAXSTR),IPL(MAXSTR),IFLOW(MAXSTR),
     $                PT2LST,IMF,IML,IO,QDUMP,ISTRS
      SAVE /ARSTRS/

      COMMON /ARDAT1/ PARA(40),MSTA(40)
      SAVE /ARDAT1/

      COMMON /ARDAT3/ IWRN(40)
      SAVE /ARDAT3/

      COMMON /LUDAT1/ MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /LUDAT1/
      CHARACTER SUB*(*)

C...Write out common message
      IF(IWRN(IERR).LT.MSTA(10)) WRITE(MSTA(8),1000) SUB,IERR,MSTA(4)
      MSTA(13)=IERR
      IWRN(IERR)=IWRN(IERR)+1
      IFATAL=0
      IDUMP=0

C...Check error code and write appropriate message
      IF(IERR.EQ.1) THEN
        WRITE(MSTA(8),1010)
        WRITE(MSTA(8),1001) ILINE
        IFATAL=1
        IDUMP=1
      ELSEIF(IERR.EQ.2) THEN
        WRITE(MSTA(8),1020)
        WRITE(MSTA(8),1001) ILINE
        IFATAL=1
        IDUMP=1
      ELSEIF(IERR.EQ.3) THEN
        IF(IWRN(3).GT.MSTA(10)) RETURN
        IWRN(3)=IWRN(3)+1
        WRITE(MSTA(8),1030)
        IF(IWRN(3).EQ.MSTA(10)) THEN
          WRITE(MSTA(8),1001) ILINE
          IDUMP=1
        ENDIF
      ELSEIF(IERR.EQ.4) THEN
        WRITE(MSTA(8),1040)
        WRITE(MSTA(8),1001) ILINE
        IFATAL=1
        IDUMP=1
      ELSEIF(IERR.EQ.5) THEN
        WRITE(MSTA(8),1050)
        WRITE(MSTA(8),1001) ILINE
        IFATAL=1
        IDUMP=1
      ELSEIF(IERR.EQ.6) THEN
        WRITE(MSTA(8),1060) MAXPAR
        IFATAL=1
      ELSEIF(IERR.EQ.7) THEN
        WRITE(MSTA(8),1070) MAXDIP
        IFATAL=1
      ELSEIF(IERR.EQ.8) THEN
        WRITE(MSTA(8),1080) MAXSTR
        IFATAL=1
      ELSEIF(IERR.EQ.9) THEN
        IF(IWRN(9).GT.MSTA(10)) RETURN
        WRITE(MSTA(8),1090)
        IF(IWRN(9).EQ.MSTA(10)) IDUMP=1
      ELSEIF(IERR.EQ.10) THEN
        IF(IWRN(10).GT.MSTA(10)) RETURN
        WRITE(MSTA(8),1100)
      ELSEIF(IERR.EQ.11) THEN
        WRITE(MSTA(8),1110)
        IFATAL=1
        IDUMP=1
      ELSEIF(IERR.EQ.12) THEN
        WRITE(MSTA(8),1120)
        IFATAL=1
      ELSEIF(IERR.EQ.13) THEN
        IF(IWRN(13).GT.MSTA(10)) RETURN
        WRITE(MSTA(8),1130)
      ELSEIF(IERR.EQ.14) THEN
        WRITE(MSTA(8),1140)
        IFATAL=1
      ELSEIF(IERR.EQ.20) THEN
        IF(IWRN(20).GT.MSTA(10)) RETURN
        WRITE(MSTA(8),1200)
      ELSEIF(IERR.EQ.21) THEN
        IF(IWRN(21).GT.MSTA(10)) RETURN
        WRITE(MSTA(8),1210)
      ENDIF

C...Dump ariadne dipole record and list the event if necessary
      IF(IDUMP.GT.0) THEN
        IF(.NOT.QDUMP) CALL ARDUMP
        WRITE(MSTA(8),1002)
        CALL LULIST(2)
      ENDIF

C...Stop execution if necessary
      IF(IFATAL.GT.0) THEN
        WRITE(MSTA(8),1003)
        STOP 0
      ENDIF

 1000 FORMAT('*** ERROR Found by Ariadne ***'/'In routine ',A6,
     $     '. Error type =',I3,'. Ariadne call number:',I7)
 1001 FORMAT('Line number:',I4)
 1002 FORMAT('Dump of event follows:')
 1003 FORMAT('Error is fatal. Execution stopped.')

 1010 FORMAT('Found colour-singlet particle in string.')
 1020 FORMAT('Found colour-triplet particle in string.')
 1030 FORMAT('Found colour-singlet particle in string.',
     $       ' Will try to cope...')
 1040 FORMAT('Found colour-triplet particle in purely gluonic string.')
 1050 FORMAT('Inconsistent colour flow in string.')
 1060 FORMAT('Maximum number of partons (',I5,') exceeded. See manual.')
 1070 FORMAT('Maximum number of dipoles (',I5,') exceeded. See manual.')
 1080 FORMAT('Maximum number of strings (',I5,') exceeded. See manual.')
 1090 FORMAT('Four-momentum was not conserved.')
 1100 FORMAT('Particle has inconsistent four-momentum. ',
     $     'Will try to cope...')
 1110 FORMAT('Recoil transfer for Drell-Yan process was not',
     $       ' kinematically allowed.')
 1120 FORMAT('AREXEC called before initialization with ARINIT.')
 1130 FORMAT('Dipole has inconsistent mass. Will try to cope...')
 1140 FORMAT('Unphysical boost vector.',/,
     $     'Try switching to double precision - see manual')
 1200 FORMAT('Selected sub-process in PYTHIA is not suported by',
     $  ' Ariadne.',/,
     $  '(only processes 11,12,13,28,53,68 are currently supported)',/,
     $  'Will try to continue but results may not be meaningful.')
 1210 FORMAT('Too many jets. ARCLUS not able to order jets in energy.')

      RETURN

C**** END OF ARERRM ****************************************************
      END
C***********************************************************************
C $Id: arexec.F,v 0.9 1992/03/13 15:34:17 lonnblad Exp $

      SUBROUTINE AREXEC

C...ARiadne subroutine EXECute ariadne

C...The Main driver routine in Ariadne.


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)

      COMMON /ARPART/ BP(MAXPAR,5),IFL(MAXPAR),IEX(MAXPAR),QQ(MAXPAR),
     $                IDI(MAXPAR),IDO(MAXPAR),INO(MAXPAR),IPART
      SAVE /ARPART/

      COMMON /ARDIPS/ BX1(MAXDIP),BX3(MAXDIP),PT2IN(MAXDIP),
     $                SDIP(MAXDIP),IP1(MAXDIP),IP3(MAXDIP),
     $                AEX1(MAXDIP),AEX3(MAXDIP),QDONE(MAXDIP),
     $                QEM(MAXDIP),IRAD(MAXDIP),ISTR(MAXDIP),IDIPS
      SAVE /ARDIPS/

      COMMON /ARSTRS/ IPF(MAXSTR),IPL(MAXSTR),IFLOW(MAXSTR),
     $                PT2LST,IMF,IML,IO,QDUMP,ISTRS
      SAVE /ARSTRS/

      COMMON /ARDAT1/ PARA(40),MSTA(40)
      SAVE /ARDAT1/

      COMMON /ARDAT2/ PQMAS(10)
      SAVE /ARDAT2/

      COMMON /ARDAT3/ IWRN(40)
      SAVE /ARDAT3/

      COMMON /ARINT2/ DBEX,DBEY,DBEZ,PHI,THE
      SAVE /ARINT2/

      COMMON /ARJETX/ N,K(300,5),P(300,5),V(300,5)
      SAVE /ARJETX/

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,XQ2,U
      SAVE /LEPTOU/

      COMMON /RYPARS/ MSTP(200),PARP(200),MSTI(200),PARI(200)
      SAVE /RYPARS/

      COMMON /RYSUBS/ MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200)
      SAVE /RYSUBS/

      COMMON /RYINT1/ MINT(400),VINT(400)
      SAVE /RYINT1/

C...Step counter
      MSTA(4)=MSTA(4)+1

C...Reset error log
      MSTA(13)=0

C...Error if ARINIT has not been called
      IF(MSTA(2).EQ.0) CALL ARERRM('AREXEC',12,0)

C...If ariadne mode just pass event through to ARPARS
      IF(MSTA(1).EQ.0) THEN

      CALL ARPARS(1,N)

C...If JETSET mode should work by just passing event on to ARPARS
      ELSEIF(MSTA(1).EQ.1) THEN
        CALL ARPARS(1,N)

C...If PYTHIA mode tag extended partons etc.
      ELSEIF(MSTA(1).EQ.2) THEN

        ISUB=MINT(1)
        IF(ISUB.NE.11.AND.ISUB.NE.12.AND.ISUB.NE.13.AND.
     $     ISUB.NE.28.AND.ISUB.NE.53.AND.ISUB.NE.68)
     $       CALL ARERRM('AREXEC',20,0)

        IFIRST=1
        ILAST=N

        DO 100 I=IFIRST,ILAST
          IF(K(I,1).GT.2) GOTO 100
          CALL ARGTYP(I,ITYP)
          IF(ITYP.EQ.0) GOTO 100
          IF(K(I,3).EQ.1.OR.K(I,3).EQ.2) THEN
            K(I,4)=1
          ELSE
            K(I,4)=0
          ENDIF
 100    CONTINUE

        CALL ARPARS(IFIRST,ILAST)

C...If LEPTO mode tag extended partons
      ELSEIF(MSTA(1).EQ.3) THEN
        IF(LST(24).EQ.1) THEN

C...Boost to hadronic cm to avoid precision problems
          DEL=DBLE(P(5,4))+DBLE(P(6,4))
          DBXL=(DBLE(P(5,1))+DBLE(P(6,1)))/DEL
          DBYL=(DBLE(P(5,2))+DBLE(P(6,2)))/DEL
          DBZL=(DBLE(P(5,3))+DBLE(P(6,3)))/DEL
          CALL LUDBRB(5,N,0.0,0.0,-DBXL,-DBYL,-DBZL)

          IF(MSTA(30).LT.2) THEN
            K(5,4)=0
          ELSE
            K(5,4)=3
            PARA(13)=SQRT(XQ2)
          ENDIF
          IF(MSTA(30).EQ.0) THEN
            K(6,4)=1
          ELSE
            K(6,4)=2
            PARA(12)=PARA(11)/(1.0-X)
          ENDIF
          CALL ARPARS(5,6)
          CALL LUDBRB(5,N,0.0,0.0,DBXL,DBYL,DBZL)
        ELSEIF(LST(24).EQ.3) THEN

C...Boost to hadronic cm to avoid precision problems
          DEL=DBLE(P(5,4))+DBLE(P(6,4))+DBLE(P(7,4))+DBLE(P(8,4))
          DBXL=(DBLE(P(5,1))+DBLE(P(6,1))+
     $         DBLE(P(7,1))+DBLE(P(8,1)))/DEL
          DBYL=(DBLE(P(5,2))+DBLE(P(6,2))+
     $         DBLE(P(7,2))+DBLE(P(8,2)))/DEL
          DBZL=(DBLE(P(5,3))+DBLE(P(6,3))+
     $         DBLE(P(7,3))+DBLE(P(8,4)))/DEL
          CALL LUDBRB(5,N,0.0,0.0,-DBXL,-DBYL,-DBZL)

          IF(MSTA(30).LT.2) THEN
            K(5,4)=0
          ELSE
            K(5,4)=3
            PARA(13)=SQRT(XQ2)
          ENDIF
          IF(MSTA(30).EQ.0) THEN
            K(6,4)=1
          ELSE
            K(6,4)=2
            PARA(12)=PARA(11)/(1.0-X)
          ENDIF
          CALL ARPARS(5,6)
          IF(MSTA(30).LT.2) THEN
            K(7,4)=0
          ELSE
            K(7,4)=3
            PARA(13)=SQRT(XQ2)
          ENDIF
          IF(MSTA(30).EQ.0) THEN
            K(8,4)=1
          ELSE
            K(8,4)=2
            PARA(12)=PARA(11)/(1.0-X)
          ENDIF
          CALL ARPARS(7,8)
          CALL LUDBRB(5,N,0.0,0.0,DBXL,DBYL,DBZL)
        ENDIF
      ENDIF

C...Perform fragmentation if requested
      IF(MSTA(5).EQ.1) CALL LUEXEC

      RETURN

C**** END OF AREXEC ****************************************************
      END
C***********************************************************************
C $Id: argpt2.F,v 0.4 1992/04/28 13:19:22 lonnblad Exp $

      REAL FUNCTION ARGPT2(ID)

C...ARiadne function Generate PT2

C...Returns the p_t^2 for a possible emission from dipole ID.


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)

      COMMON /ARDAT1/ PARA(40),MSTA(40)
      SAVE /ARDAT1/

      COMMON /ARPART/ BP(MAXPAR,5),IFL(MAXPAR),IEX(MAXPAR),QQ(MAXPAR),
     $                IDI(MAXPAR),IDO(MAXPAR),INO(MAXPAR),IPART
      SAVE /ARPART/

      COMMON /ARDIPS/ BX1(MAXDIP),BX3(MAXDIP),PT2IN(MAXDIP),
     $                SDIP(MAXDIP),IP1(MAXDIP),IP3(MAXDIP),
     $                AEX1(MAXDIP),AEX3(MAXDIP),QDONE(MAXDIP),
     $                QEM(MAXDIP),IRAD(MAXDIP),ISTR(MAXDIP),IDIPS
      SAVE /ARDIPS/

      COMMON /ARSTRS/ IPF(MAXSTR),IPL(MAXSTR),IFLOW(MAXSTR),
     $                PT2LST,IMF,IML,IO,QDUMP,ISTRS
      SAVE /ARSTRS/


C...Set invariant mass squared in the dipole and generate a p_t^2
C...with the appropriate Monte Carlo subroutine
      IF(QEM(ID).AND.MSTA(20).GE.2.AND.ISTRS.GE.2) THEN
        PT2IN(ID)=0.0
        QDONE(ID)=.TRUE.
      ENDIF
      IF(.NOT.QDONE(ID)) THEN
        SDIP(ID)=ARMAS2(IP1(ID),IP3(ID))
        IF(QEM(ID)) THEN
          CALL ARGQED(ID)
        ELSE
          CALL ARGQCD(ID)
        ENDIF
        QDONE(ID)=.TRUE.
      ENDIF

      ARGPT2=PT2IN(ID)

      RETURN

C**** END OF ARGPT2 ****************************************************
      END
C***********************************************************************
C $Id: argqcd.F,v 0.12 1992/03/11 08:56:44 lonnblad Exp $

      SUBROUTINE ARGQCD(ID)

C...ARiadne subroutine Generate pt2 for QCD emission.

C...Generates a p_t^2 for a possible QCD emission from dipole ID


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)
      
      double precision PYR

      COMMON /ARPART/ BP(MAXPAR,5),IFL(MAXPAR),IEX(MAXPAR),QQ(MAXPAR),
     $                IDI(MAXPAR),IDO(MAXPAR),INO(MAXPAR),IPART
      SAVE /ARPART/

      COMMON /ARDIPS/ BX1(MAXDIP),BX3(MAXDIP),PT2IN(MAXDIP),
     $                SDIP(MAXDIP),IP1(MAXDIP),IP3(MAXDIP),
     $                AEX1(MAXDIP),AEX3(MAXDIP),QDONE(MAXDIP),
     $                QEM(MAXDIP),IRAD(MAXDIP),ISTR(MAXDIP),IDIPS
      SAVE /ARDIPS/

      COMMON /ARSTRS/ IPF(MAXSTR),IPL(MAXSTR),IFLOW(MAXSTR),
     $                PT2LST,IMF,IML,IO,QDUMP,ISTRS
      SAVE /ARSTRS/

      COMMON /ARDAT1/ PARA(40),MSTA(40)
      SAVE /ARDAT1/

      COMMON /ARDAT2/ PQMAS(10)
      SAVE /ARDAT2/

      COMMON /ARINT1/ BC1,BC3,BZM,BZP,BP1,BM1,BP3,BM3,
     $                B1,B2,B3,XT2,XT,Y,QQ1,QQ3,NE1,NE3,
     $                S,W,C,CN,ALPHA0,XLAM2,IFLG,
     $                XT2MP,XT2ME,XT2M,XT2C,XTS,XT3,XT1,
     $                YINT,YMAX,YMIN,
     $                Y1,Y2,Y3,SY1,SY2,SY3,SSY,
     $                AE1,AE3,NXP1,NXP3,FQ1,FQ3
      SAVE /ARINT1/

      COMMON /LUDAT1/ MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /LUDAT1/
      EXTERNAL ARNDX1,ARNDX2,ARNDX3,ARNDY1,ARNDY2,ARNDY3,ARNDY4,
     $         ARVET3,ARVET4,ARVET5
      REAL ARNDX1,ARNDX2,ARNDX3,ARNDY1,ARNDY2,ARNDY3,ARNDY4,
     $         ARVET3,ARVET4,ARVET5


C...Copy some information from dipole record
C...S     = the invariant mass squared
C...W     = total energy in dipole
C...XT2MP = maximum allowed fractional p_t^2 (x_t^2) for restricted  
C...        phase space option
C...QQ1(3)= Boolean variable 'is quark' for parton 1(3)
C...NE1(3)= integer determining extention of parton 1(3) (0=pointlike)
C...SY1(3)= fractional mass of parton 1(3)
      PT2IN(ID)=0.0
      S=SDIP(ID)
      IF(S.LE.4.0*PARA(3)**2) RETURN
      W=SQRT(S)
      XT2MP=PT2LST/S
      QQ1=QQ(IP1(ID))
      QQ3=QQ(IP3(ID))
      NE1=IEX(IP1(ID))
      NE3=IEX(IP3(ID))
      SY1=BP(IP1(ID),5)/W
      SY3=BP(IP3(ID),5)/W

      GOTO 100

C...Special entry for checking p_t-ordering of recoil gluons
      ENTRY ARGQTE(ID,SI,XT2MPI,QQ1I,QQ3I,NE1I,NE3I,SY1I,SY3I)
            PT2IN(ID)=0.0
            S=SI
            IF(S.LE.4.0*PARA(3)**2) RETURN
            W=SQRT(S)
            XT2MP=XT2MPI
            QQ1=QQ1I
            QQ3=QQ3I
            NE1=NE1I
            NE3=NE3I
            SY1=SY1I
            SY3=SY3I

 100  SY2=0.0

C...Calculate maximum x_t^2 for extended dipole
      IF(NE1.GT.0.AND.NE3.EQ.0) XT2ME=((0.25*S*(PARA(10+NE1)**
     $                          PARA(10)))**(2.0/(2.0+PARA(10))))/S
      IF(NE1.EQ.0.AND.NE3.GT.0) XT2ME=((0.25*S*(PARA(10+NE3)**
     $                          PARA(10)))**(2.0/(2.0+PARA(10))))/S
      IF(NE1.GT.0.AND.NE3.GT.0) XT2ME=((0.25*S*((PARA(10+NE1)*
     $           PARA(10+NE3))**PARA(10)))**(1.0/(1.0+PARA(10))))/S

C...XLAM = scaled lambda_QCD squared
      XLAM2=PARA(1)**2/S

C...C = colour factors etc. in cross section
      C=6.0/(4.0*PARU(1))
      IF(QQ1.AND.QQ3) C=4.0/(3.0*PARU(1))

C...alpha_0 for alpha_QCD = alpha_0/ln(p_t^2/lambda_QCD^2)
      ALPHA0=12.0*PARU(1)/(33.0-2.0*MAX(ARNOFL(W,MAX(5,MSTA(15))),3.0))

C...Set exponents in cross section
      NXP1=3
      NXP3=3
      IF(QQ1) NXP1=2
      IF(QQ3) NXP3=2

C...Flavour of this emission 0 = gluon emission
      IFLG=0

C...Minimum x_t^2
      XT2C=PARA(3)**2/S

C...Calculate mass dependent parameters
      CALL ARMADE

C...Set maximum x_t^2
      IF(MSTA(11).LT.4) XT2M=MIN(XT2M,XT2MP)
      IF(NE1.GT.0.OR.NE3.GT.0) XT2M=MIN(XT2M,XT2ME)

      IF(XT2M.LE.XT2C) THEN
        PT2IN(ID)=0.0
        RETURN
      ENDIF

C...Set additional parameters and call the veto algorith with
C...Suitable random functions
      IF(MSTA(12).GT.0) THEN
C.......Running alpha_QDC
        YINT=2.0*LOG(0.5/SQRT(XLAM2)+SQRT(0.25/XLAM2-1.0))
        CN=1.0/(YINT*C*ALPHA0)
        IF(NE1.GT.0.OR.NE3.GT.0) THEN
C.........Extended dipole
          CALL ARMCDI(ARNDX1,ARNDY2,ARVET4)
        ELSE
C.........Pointlike dipole
          CALL ARMCDI(ARNDX1,ARNDY1,ARVET4)
        ENDIF
      ELSE
C.......Constant alpha_QCD
        YINT=1.0
        CN=2.0/(C*PARA(2))
        IF(NE1.GT.0.OR.NE3.GT.0) THEN
C.........Extended dipole
          CALL ARMCDI(ARNDX2,ARNDY2,ARVET3)
        ELSE
C.........Pointlike dipole
          CALL ARMCDI(ARNDX2,ARNDY1,ARVET3)
        ENDIF
      ENDIF

C...Save the generated values of p_t^2, x1, x3, a1 and a3
      PT2IN(ID)=XT2*S
      BX1(ID)=B1
      BX3(ID)=B3
      AEX1(ID)=AE1
      AEX3(ID)=AE3
      IRAD(ID)=0

C...Exit if no q-qbar emission
      IF(MSTA(15).LE.0) RETURN
      QG1=((.NOT.QQ1).AND.NE1.EQ.0)
      QG3=((.NOT.QQ3).AND.NE3.EQ.0)
      IF((.NOT.QG1).AND.(.NOT.QG3)) RETURN

C...Colour factors and things in cross section. If g-g dipole
C...q-qbar splitting only calculated forone gluon but double
C...cross section
      C=1.0/(8.0*PARU(1))
      IF(QG1.AND.QG3) C=C*2.0

C...Parton 3 is always assumed to be split
      IF(QG1) THEN
        SY1=SY3
        NE1=NE3
        NE3=0
      ENDIF
C...set 'minimum' XT2 to the XT2 of the gluon emission. XT2s below that
C...are not relevant
      XT2C=MAX(XT2,XT2C)

C...Loop over allowed flavours
      DO 200 IFLG=1,MSTA(15)

C...Set mass dependent parameters
        SY2=PQMAS(IFLG)/W
        SY3=SY2
        CALL ARMADE

C...Set phase space restrictions
        IF(MSTA(11).LT.2) XT2M=MIN(XT2M,XT2MP)
        IF(NE1.GT.0.OR.NE3.GT.0) XT2M=MIN(XT2M,XT2ME)

C...Exit if not enough energy
        IF(XT2M.LE.XT2C.OR.SSY.GE.1.0) GOTO 300

C...Set additional parameters and call the veto algorith with
C...Suitable random functions
        YINT=2.0*SQRT(S)
C.......Running alpha_QCD
        IF(MSTA(12).GT.0) THEN
          CN=1.0/(YINT*C*ALPHA0)
          IF(NE1.GT.0.OR.NE3.GT.0) THEN
C...........pointlike dipole
            CALL ARMCDI(ARNDX1,ARNDY4,ARVET5)
          ELSE
C...........extended dipole
            CALL ARMCDI(ARNDX1,ARNDY3,ARVET5)
          ENDIF
        ELSE
C.........Constant alpha_QCD
          CN=2.0/(YINT*C*PARA(2))
          CN=1.0/(YINT*C*ALPHA0)
          IF(NE1.GT.0.OR.NE3.GT.0) THEN
C...........pointlike dipole
            CALL ARMCDI(ARNDX3,ARNDY4,ARVET5)
          ELSE
C...........extended dipole
            CALL ARMCDI(ARNDX3,ARNDY3,ARVET5)
          ENDIF
        ENDIF

C...If Generated XT2 is larger than previous XT2 accept this and save
C...the generated values of p_t^2, x1, x3, a1 and a3
        IF(XT2.GT.XT2C) THEN
          PT2IN(ID)=XT2*S
          BX1(ID)=B1
          BX3(ID)=B3
          AEX1(ID)=AE1
          AEX3(ID)=AE3
          IRAD(ID)=IFLG
          XT2C=XT2
        ENDIF

 200  CONTINUE

C...Exit if gluon emission was chosen 
 300  IF(IRAD(ID).EQ.0) RETURN

C...Select wich gluon to split
      IF((.NOT.QG3).OR.(QG1.AND.PYR(IDUM).GT.0.5)) THEN
        IRAD(ID)=-IRAD(ID)
        B1=BX1(ID)
        BX1(ID)=BX3(ID)
        BX3(ID)=B1
        AE1=AEX1(ID)
        AEX1(ID)=AEX3(ID)
        AEX3(ID)=AE1
      ENDIF

      RETURN

C**** END OF ARGQCD ****************************************************
      END
C***********************************************************************
C $Id: argqed.F,v 0.10 1992/04/30 10:58:34 lonnblad Exp $

      SUBROUTINE ARGQED(ID)

C...ARiadne subroutine Generate pt2 for QED emission

C...Generates a p-t^2 for a possible QED emission from dipole ID.


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)

      COMMON /ARPART/ BP(MAXPAR,5),IFL(MAXPAR),IEX(MAXPAR),QQ(MAXPAR),
     $                IDI(MAXPAR),IDO(MAXPAR),INO(MAXPAR),IPART
      SAVE /ARPART/

      COMMON /ARDIPS/ BX1(MAXDIP),BX3(MAXDIP),PT2IN(MAXDIP),
     $                SDIP(MAXDIP),IP1(MAXDIP),IP3(MAXDIP),
     $                AEX1(MAXDIP),AEX3(MAXDIP),QDONE(MAXDIP),
     $                QEM(MAXDIP),IRAD(MAXDIP),ISTR(MAXDIP),IDIPS
      SAVE /ARDIPS/

      COMMON /ARSTRS/ IPF(MAXSTR),IPL(MAXSTR),IFLOW(MAXSTR),
     $                PT2LST,IMF,IML,IO,QDUMP,ISTRS
      SAVE /ARSTRS/

      COMMON /ARDAT1/ PARA(40),MSTA(40)
      SAVE /ARDAT1/

      COMMON /ARINT1/ BC1,BC3,BZM,BZP,BP1,BM1,BP3,BM3,
     $                B1,B2,B3,XT2,XT,Y,QQ1,QQ3,NE1,NE3,
     $                S,W,C,CN,ALPHA0,XLAM2,IFLG,
     $                XT2MP,XT2ME,XT2M,XT2C,XTS,XT3,XT1,
     $                YINT,YMAX,YMIN,
     $                Y1,Y2,Y3,SY1,SY2,SY3,SSY,
     $                AE1,AE3,NXP1,NXP3,FQ1,FQ3
      SAVE /ARINT1/

      COMMON /LUDAT1/ MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /LUDAT1/
      EXTERNAL ARNDX2,ARNDY1,ARVET1,ARVET2
      REAL ARNDX2,ARNDY1,ARVET1,ARVET2


C...Copy information about partons in dipole (for explanation see
C...subroutine ARGQCD
      PT2IN(ID)=0.0
      S=SDIP(ID)
      IF(S.LE.4.0*PARA(5)**2) RETURN
      IF(MSTA(20).GE.2.AND.ISTRS.GE.2) RETURN
      W=SQRT(S)
      XT2MP=PT2LST/S
      QQ1=QQ(IP1(ID))
      QQ3=QQ(IP3(ID))
      NE1=IEX(IP1(ID))
      NE3=IEX(IP3(ID))

      SY1=BP(IP1(ID),5)/W
      SY2=0.0
      SY3=BP(IP3(ID),5)/W

      XT2C=PARA(5)**2/S
      NXP1=2
      NXP3=2

C...Set charges of emitting quarks and set constant in cross section
      IQ1=LUCHGE(IFL(IP1(ID)))
      IQ3=LUCHGE(IFL(IP3(ID)))
      FQMAX=FLOAT(MAX(ABS(IQ1),ABS(IQ3)))
      FQ1=FLOAT(IQ1)/FQMAX
      FQ3=FLOAT(IQ3)/FQMAX
      C=(FQMAX**2)/(9.0*PARU(1))
      IFLG=-1

C...Set mass dependent parameters
      CALL ARMADE

C...Restrict phase space if demanded
      IF(MSTA(11).EQ.0.OR.MSTA(11).EQ.2) XT2M=MIN(XT2M,XT2MP)

C...Set some further parameters and call the veto algorithm with
C...suitable random functions for constant alpha_EM.
      YINT=1.0
      CN=2.0/(C*PARA(4))
      CALL ARMCDI(ARNDX2,ARNDY1,ARVET1)

C...Save information about emission
      PT2IN(ID)=XT2*S
      BX1(ID)=B1
      BX3(ID)=B3
      AEX1(ID)=AE1
      AEX3(ID)=AE3

      RETURN

C**** END OF ARGQED ****************************************************
      END
C***********************************************************************
C $Id: argtyp.F,v 0.3 1991/09/26 12:42:47 lonnblad Exp $

      SUBROUTINE ARGTYP(I,ITYP)

C...ARiadne subroutine Get TYpe of Particle

C...Determines the type of particle I according to ITYP=2: gluon,
C...ITYP=1: quark or anti-di-quark, ITYP=-1: anti-quark or di-quark,
C...ITYP=0: other.


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)

      COMMON /ARJETX/ N,K(300,5),P(300,5),V(300,5)
      SAVE /ARJETX/

      COMMON /LUDAT2/ KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      SAVE /LUDAT2/


      if (LUCOMP(K(I,2)).le.0) then
         ityp = 0
      else
         ITYP=KCHG(LUCOMP(K(I,2)),2)*ISIGN(1,K(I,2))
      endif

      RETURN

C**** END OF ARGTYP ****************************************************
      END
C***********************************************************************
C $Id: aript2.F,v 0.4 1992/01/27 16:03:19 lonnblad Exp $

      REAL FUNCTION ARIPT2(I1,I2,I3)

C...ARiadne function Invariant PT2

C...Returns the invariant p_t^2 of three partons


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)

      COMMON /ARPART/ BP(MAXPAR,5),IFL(MAXPAR),IEX(MAXPAR),QQ(MAXPAR),
     $                IDI(MAXPAR),IDO(MAXPAR),INO(MAXPAR),IPART
      SAVE /ARPART/


      ARIPT2=(ARMAS2(I1,I2)-(BP(I1,5)+BP(I2,5))**2)*
     $       (ARMAS2(I2,I3)-(BP(I2,5)+BP(I3,5))**2)/
     $        ARMAS3(I1,I2,I3)

      RETURN

C**** END OF ARIPT2 ****************************************************
      END
C***********************************************************************
C $Id: arjoin.F,v 0.6 1992/01/27 16:03:19 lonnblad Exp $

      SUBROUTINE ARJOIN(J1,J2,J3)

C...ARiadne subroutine JOIN jets

C...Join three jets into two


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)

      COMMON /ARPART/ BP(MAXPAR,5),IFL(MAXPAR),IEX(MAXPAR),QQ(MAXPAR),
     $                IDI(MAXPAR),IDO(MAXPAR),INO(MAXPAR),IPART
      SAVE /ARPART/

      COMMON /ARJETX/ N,K(300,5),P(300,5),V(300,5)
      SAVE /ARJETX/

      COMMON /LUDAT1/ MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /LUDAT1/


C...Copy jets into ariadne dipole record
      CALL ARCOPA(J1,1,1)
      CALL ARCOPA(J2,2,1)
      CALL ARCOPA(J3,3,1)

C...Boost to CMS of jets
      DE=BP(1,4)+BP(2,4)+BP(3,4)
      DBEX=(BP(1,1)+BP(2,1)+BP(3,1))/DE
      DBEY=(BP(1,2)+BP(2,2)+BP(3,2))/DE
      DBEZ=(BP(1,3)+BP(2,3)+BP(3,3))/DE

      CALL AROBO3(0.0,0.0,-DBEX,-DBEY,-DBEZ,1,2,3)

C...Rotate Jet 1 to z-axis and jet 2 to xz plane
      PX=BP(1,1)
      PY=BP(1,2)
      PZ=BP(1,3)
      PHI=ULANGL(PX,PY)
      CALL AROBO3(0.0,-PHI,0.0D0,0.0D0,0.0D0,1,2,3)
      THE=ULANGL(PZ,SQRT(PX**2+PY**2))
      CALL AROBO3(-THE,0.0,0.0D0,0.0D0,0.0D0,1,2,3)
      PX=BP(2,1)
      PY=BP(2,2)
      PHI2=ULANGL(PX,PY)
      CALL AROBO3(0.0,-PHI2,0.0D0,0.0D0,0.0D0,1,2,3)

C...Calculate energy fractions
      BE=BP(1,4)+BP(2,4)+BP(3,4)
      B1=2.0*BP(1,4)/BE
      B3=2.0*BP(3,4)/BE

C...Determine recoil angle
      BET=ARANGL(1,3)
      PSI=(PARU(1)-BET)*(B3**2)/(B1**2+B3**2)
      BP(1,1)=0.0
      BP(1,2)=0.0
      BP(1,3)=BE*0.5
      BP(1,4)=BE*0.5
      BP(1,5)=0.0
      BP(2,1)=0.0
      BP(2,2)=0.0
      BP(2,3)=-BE*0.5
      BP(2,4)=BE*0.5
      BP(2,5)=0.0

C...Rotate and boost back
      CALL AROBO2(PSI,0.0,0.0D0,0.0D0,0.0D0,1,2)
      CALL AROBO2(0.0,PHI2,0.0D0,0.0D0,0.0D0,1,2)
      CALL AROBO2(THE,PHI,DBEX,DBEY,DBEZ,1,2)

C...Copy jets to /LUJETS/
      DO 100 J=1,5
        P(J1,J)=BP(1,J)
        P(J3,J)=BP(2,J)
 100  CONTINUE
      V(J1,1)=0.0
      V(J3,1)=0.0

      RETURN

C**** END OF ARJOIN ****************************************************
      END
C***********************************************************************
C $Id: armade.F,v 0.4 1992/01/27 16:03:19 lonnblad Exp $

      SUBROUTINE ARMADE

C...ARiadne subroutine set MAss DEpendencies

C...Sets some mass dependencies needed for ARMCDI


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)

      COMMON /ARINT1/ BC1,BC3,BZM,BZP,BP1,BM1,BP3,BM3,
     $                B1,B2,B3,XT2,XT,Y,QQ1,QQ3,NE1,NE3,
     $                S,W,C,CN,ALPHA0,XLAM2,IFLG,
     $                XT2MP,XT2ME,XT2M,XT2C,XTS,XT3,XT1,
     $                YINT,YMAX,YMIN,
     $                Y1,Y2,Y3,SY1,SY2,SY3,SSY,
     $                AE1,AE3,NXP1,NXP3,FQ1,FQ3
      SAVE /ARINT1/


      SSY=SY1+SY2+SY3
      Y1=SY1**2
      Y2=SY2**2
      Y3=SY3**2

      BC1=-DBLE(Y1)-1.0D0+DBLE(SY2+SY3)**2
      BC3=-DBLE(Y3)-1.0D0+DBLE(SY2+SY1)**2
      XT2M=0.0
      IF(SQRT(0.25+Y2)-1.0-(BC1+BC3)/2.0.LT.0.0) RETURN
      XTS=(SQRT(0.25+Y2)-1.0-(BC1+BC3)/2.0)**2
      XT1=-2.0*SY1-BC1
      XT3=-2.0*SY3-BC3
      IF(XT1.LT.0.0) RETURN
      IF(XT3.LT.0.0) RETURN
      XT2M=MIN(XTS,XT1*XT3)

      BZP=0.5*(1.0+Y1-Y3+SQRT(1.0+(Y1-Y3)**2-2.0*(Y1+Y3)))
      BZM=0.5*(1.0+Y3-Y1+SQRT(1.0+(Y1-Y3)**2-2.0*(Y1+Y3)))

      RETURN

C**** END OF ARMADE ****************************************************
      END
C***********************************************************************
C $Id: armass.F,v 0.6 1992/01/27 16:03:19 lonnblad Exp $

      REAL FUNCTION ARMAS2(I1,I2)

C...ARiadne function invariant MASs of 2 partons

C...Returns the invariant mass^2 of partons I1 and I2

      DIMENSION I(2)


      I(1)=I1
      I(2)=I2

      ARMAS2=ARMASS(2,I)

      RETURN

C**** END OF ARMAS2 ****************************************************
      END
C***********************************************************************
C $Id: armass.F,v 0.6 1992/01/27 16:03:19 lonnblad Exp $

      REAL FUNCTION ARMAS3(I1,I2,I3)

C...ARiadne function invariant MASs of 3 partons

C...Returns the invariant mass^2 of partons I1, I2 and I3

      DIMENSION I(3)


      I(1)=I1
      I(2)=I2
      I(3)=I3
      
      ARMAS3=ARMASS(3,I)

      RETURN

C**** END OF ARMAS3 ****************************************************
      END
C***********************************************************************
C $Id: armass.F,v 0.6 1992/01/27 16:03:19 lonnblad Exp $

      REAL FUNCTION ARMASS(N,I)

C...ARiadne function invariant MASS of partons

C...Returns the total invariant mass^2 of N partons


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)

      COMMON /ARPART/ BP(MAXPAR,5),IFL(MAXPAR),IEX(MAXPAR),QQ(MAXPAR),
     $                IDI(MAXPAR),IDO(MAXPAR),INO(MAXPAR),IPART
      SAVE /ARPART/
      DIMENSION I(N),DPS(4)


      DO 100 IK=1,4
        DPS(IK)=0.0D0
        DO 200 IJ=1,N
          DPS(IK)=DPS(IK)+BP(I(IJ),IK)
 200    CONTINUE
 100  CONTINUE

      DMASS=DPS(4)**2-DPS(3)**2-DPS(2)**2-DPS(1)**2
      ARMASS=MAX(DMASS,0.0D0)

      RETURN

C**** END OF ARMASS ****************************************************
      END
C***********************************************************************
C $Id: armcdi.F,v 0.9 1992/01/31 16:14:59 lonnblad Exp $

      SUBROUTINE ARMCDI(ARRNDX,ARRNDY,ARVETO)

C...ARiadne subroutine Monte Carlo DIstribution

C...Generates x_1 and x_3 for a radiating dipole


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)

      double precision PYR

      COMMON /ARINT1/ BC1,BC3,BZM,BZP,BP1,BM1,BP3,BM3,
     $                B1,B2,B3,XT2,XT,Y,QQ1,QQ3,NE1,NE3,
     $                S,W,C,CN,ALPHA0,XLAM2,IFLG,
     $                XT2MP,XT2ME,XT2M,XT2C,XTS,XT3,XT1,
     $                YINT,YMAX,YMIN,
     $                Y1,Y2,Y3,SY1,SY2,SY3,SSY,
     $                AE1,AE3,NXP1,NXP3,FQ1,FQ3
      SAVE /ARINT1/


C...Exit if below cut
 100  IF(XT2M.LT.XT2C) GOTO 900

C...Generate random XT2
      XT2=ARRNDX()
      IF(XT2.LT.XT2C) GOTO 900
      XT=SQRT(XT2)

C...Generate rapidity Y
      Y=ARRNDY()

C...Calculate energy fractions
      B1=-XT*EXP(Y)-BC1
      B3=-XT*EXP(-Y)-BC3
      B2=2.0-B1-B3

C...Set maximum XT2 for possible next random call (VETO algorithm)
      XT2M=XT2

C...Redo random calls according to veto-algorithm
      IF(ARVETO().LT.PYR(IDUM)) GOTO 100

C...Check that Current values are kinematically allowed
      CALL ARCHKI(0,IOK)
      IF(IOK.EQ.0) GOTO 100

      RETURN

C...If below cuts set XT2 to 0
 900  B1=-BC1
      B3=-BC3
      XT2=0.0

      RETURN

C**** END OF ARMCDI ****************************************************
      END
C***********************************************************************
C $Id: armipt.F,v 0.3 1991/09/26 12:43:33 lonnblad Exp $

      REAL FUNCTION ARMIPT(IF,IL)

C...ARiadne function determine MInimum PT2

C...Determines the minimum p_t^2 of any gluon between positions 
C...IF and IL.


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)

      COMMON /ARPART/ BP(MAXPAR,5),IFL(MAXPAR),IEX(MAXPAR),QQ(MAXPAR),
     $                IDI(MAXPAR),IDO(MAXPAR),INO(MAXPAR),IPART
      SAVE /ARPART/

      COMMON /ARDIPS/ BX1(MAXDIP),BX3(MAXDIP),PT2IN(MAXDIP),
     $                SDIP(MAXDIP),IP1(MAXDIP),IP3(MAXDIP),
     $                AEX1(MAXDIP),AEX3(MAXDIP),QDONE(MAXDIP),
     $                QEM(MAXDIP),IRAD(MAXDIP),ISTR(MAXDIP),IDIPS
      SAVE /ARDIPS/

      COMMON /ARSTRS/ IPF(MAXSTR),IPL(MAXSTR),IFLOW(MAXSTR),
     $                PT2LST,IMF,IML,IO,QDUMP,ISTRS
      SAVE /ARSTRS/

      COMMON /ARDAT1/ PARA(40),MSTA(40)
      SAVE /ARDAT1/

      INXT(IP)=IP3(IDO(IP))
      IPRV(IP)=IP1(IDI(IP))


      ARMIPT=PARA(40)
      DO 100 I=IF,IL
        IF(.NOT.QQ(I)) ARMIPT=MIN(ARMIPT,ARIPT2(IPRV(I),I,INXT(I)))
 100  CONTINUE

      RETURN

C**** END OF ARMIPT ****************************************************
      END
C***********************************************************************
C $Id: arnofl.F,v 0.3 1991/09/26 12:43:36 lonnblad Exp $

      REAL FUNCTION ARNOFL(W,MNOFL)

C...ARiadne function Number Of FLavours 

C...Returns the number of flavourspossible at energy W


      COMMON /ARDAT2/ PQMAS(10)
      SAVE /ARDAT2/


      ARNOFL=0.0
      DO 100 I=1,MNOFL
        IF(W.LT.2.0*PQMAS(I)) RETURN
        ARNOFL=FLOAT(I)
 100  CONTINUE

      RETURN

C**** END OF ARNOFL ****************************************************
      END
C***********************************************************************
C $Id: arordj.F,v 0.6 1992/01/31 16:14:59 lonnblad Exp $

      SUBROUTINE ARORDJ

C...ARiadne subroutine ORDer Jets

C...Orders jets in falling energy


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)

      COMMON /ARPART/ BP(MAXPAR,5),IFL(MAXPAR),IEX(MAXPAR),QQ(MAXPAR),
     $                IDI(MAXPAR),IDO(MAXPAR),INO(MAXPAR),IPART
      SAVE /ARPART/

      COMMON /ARJETX/ N,K(300,5),P(300,5),V(300,5)
      SAVE /ARJETX/

      COMMON /LUDAT1/ MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /LUDAT1/


C...Error if no space left in /ARPART/
      IF(MSTU(3).GT.MAXPAR) CALL ARERRM('ARORDJ',10,0)

C...Copy jets into /ARPART/ and link them with pointers
      IPF=1
      NJET=0
      DO 100 I=1,MSTU(3)
        IF(K(N+I,5).LT.0) GOTO 100
        NJET=NJET+1
        IDO(NJET)=NJET+1
        IDI(NJET)=NJET-1
        DO 110 J=1,5
          BP(NJET,J)=P(N+I,J)
 110    CONTINUE
 100  CONTINUE
      IDI(1)=0
      IDO(NJET)=0

C...Copy back jets to /LUJETS/ in falling order in energy
      MSTU(3)=NJET
      DO 200 I=1,MSTU(3)
        EMAX=0.0
        IM=0
        IP=IPF
 210    IF(BP(IP,4).GT.EMAX) THEN
          EMAX=BP(IP,4)
          IM=IP
        ENDIF
        IF(IDO(IP).NE.0) THEN
          IP=IDO(IP)
          GOTO 210
        ENDIF

        DO 220 J=1,5
          P(N+I,J)=BP(IM,J)
 220    CONTINUE

        IF(IM.EQ.IPF) THEN
          IPF=IDO(IM)
        ELSE
          IDO(IDI(IM))=IDO(IM)
          IF(IDO(IM).NE.0) IDI(IDO(IM))=IDI(IM)
        ENDIF

 200  CONTINUE

      RETURN

C**** END OF ARORDJ ****************************************************
      END
C***********************************************************************
C $Id: arorie.F,v 0.18 1992/03/12 11:23:19 lonnblad Exp $

      SUBROUTINE ARORIE(I1,I2,I3,BS,B1,B3,QR1,QR3,PT21,PT23)

C...ARiadne subroutine ORIEnt

C...Orients three partons according to recoil strategy determined
C...by QR1 and QR3


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)

      double precision PYR

      COMMON /ARPART/ BP(MAXPAR,5),IFL(MAXPAR),IEX(MAXPAR),QQ(MAXPAR),
     $                IDI(MAXPAR),IDO(MAXPAR),INO(MAXPAR),IPART
      SAVE /ARPART/

      COMMON /ARDIPS/ BX1(MAXDIP),BX3(MAXDIP),PT2IN(MAXDIP),
     $                SDIP(MAXDIP),IP1(MAXDIP),IP3(MAXDIP),
     $                AEX1(MAXDIP),AEX3(MAXDIP),QDONE(MAXDIP),
     $                QEM(MAXDIP),IRAD(MAXDIP),ISTR(MAXDIP),IDIPS
      SAVE /ARDIPS/

      COMMON /ARSTRS/ IPF(MAXSTR),IPL(MAXSTR),IFLOW(MAXSTR),
     $                PT2LST,IMF,IML,IO,QDUMP,ISTRS
      SAVE /ARSTRS/

      COMMON /LUDAT1/ MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /LUDAT1/

      INXT(I)=IP3(IDO(I))
      IPRV(I)=IP1(IDI(I))


C...Set parton energies and momentum and total energy
      BW=SQRT(BS)
      IF(B1.LE.0.0) CALL ARERRM('ARORIE',9,0)
      DE1=0.5*B1*BW
      IF(B3.LE.0.0) CALL ARERRM('ARORIE',9,0)
      DE3=0.5*B3*BW
      DE2=BW-DE1-DE3
      IF(DE2.LT.BP(I2,5)) CALL ARERRM('ARORIE',9,0)

CCPH:  Negativity may arise here occasionally.  Remove the problem brutally.
C      DP1=SQRT(DE1**2-BP(I1,5)**2)
C      DP2=SQRT(DE2**2-BP(I2,5)**2)
C      DP3=SQRT(DE3**2-BP(I3,5)**2)
      DP1=DE1**2-BP(I1,5)**2
      DP2=DE2**2-BP(I2,5)**2
      DP3=DE3**2-BP(I3,5)**2
      if(DP1.LT.0.) THEN
       DP1 = 0.D0
       BP(I1,5) = DE1
      ENDIF
      if(DP2.LT.0.) THEN
       DP2 = 0.D0
       BP(I2,5) = DE2
      ENDIF
      if(DP3.LT.0.) THEN
       DP3 = 0.D0
       BP(I3,5) = DE3
      ENDIF
      DP1=DSQRT(DP1)
      DP2=DSQRT(DP2)
      DP3=DSQRT(DP3)

C...If both partons 1 and 3 can take full recoil choose one according to
C...Kleiss
      IF(QR1.AND.QR3) THEN
        IF(B1**2.LT.(B1**2+B3**2)*PYR(IDUM)) THEN
          QR1=.FALSE.
        ELSE
          QR3=.FALSE.
        ENDIF
      ENDIF

C...Calculate angle between partons 1 and 3
      BCALP=1.0
      IF(DP1.GT.0.0.AND.DP3.GT.0.0) THEN
        BCALP=(DP2**2-DP1**2-DP3**2)/(2.0*DP1*DP3)
      ELSE
        CALL ARERRM('ARORIE',9,0)
      ENDIF
      IF(ABS(BCALP).GT.1.0) CALL ARERRM('ARORIE',9,0)
      BCALP=MAX(-1.0D0,MIN(1.0D0,DBLE(BCALP)))
      BALP=ACOS(BCALP)

C...Determine angle between parton 1 and z-axis
      IF(QR1.AND.PT21.LE.0.0.AND.PT23.LE.0.0) THEN
        BPSI=PARU(1)-BALP
      ELSEIF(QR3.AND.PT21.LE.0.0.AND.PT23.LE.0.0) THEN
        BPSI=0.0
      ELSE
        BPSI=(PARU(1)-BALP)*(B3**2)/(B1**2+B3**2)

C...New recoil strategy
        IF(PT21.GT.0.0.AND.PT21.GE.PT23) THEN
          I0=IPRV(I1)
          BPSI=ARECOI(BP(I0,4),DE1,DE2,DE3,ABS(BP(I0,3)),DP1,DP2,DP3,
     $         BALP,PT21)
        ELSEIF(PT23.GT.0.0.AND.PT23.GT.PT21) THEN
          I4=INXT(I3)
          BPSI=PARU(1)-BALP-
     $         ARECOI(BP(I4,4),DE3,DE2,DE1,ABS(BP(I4,3)),
     $         DP3,DP2,DP1,BALP,PT23)
        ENDIF
      ENDIF

C...Set random azimuth angle
      BGAM=PARU(2)*PYR(IDUM)
      BSGAM=SIN(BGAM)
      BCGAM=COS(BGAM)
      BSPSI=SIN(BPSI)
      BCPSI=COS(BPSI)
      BSPSA=SIN(BPSI+BALP)
      BCPSA=COS(BPSI+BALP)

C...Set fourmomentum of partons
      BP(I1,1)=DP1*BSPSI*BSGAM
      BP(I1,2)=-DP1*BSPSI*BCGAM
      BP(I1,3)=DP1*BCPSI
      BP(I1,4)=DE1

      BP(I3,1)=DP3*BSPSA*BSGAM
      BP(I3,2)=-DP3*BSPSA*BCGAM
      BP(I3,3)=DP3*BCPSA
      BP(I3,4)=DE3

      DZ2=-DP1*BCPSI-DP3*BCPSA
      DT2=DSQRT(MAX(DP2**2-DZ2**2,0.0D0))
      BP(I2,1)=-DT2*BSGAM
      BP(I2,2)=DT2*BCGAM
      BP(I2,3)=DZ2
      BP(I2,4)=DE2

      RETURN

C**** END OF ARORIE ****************************************************
      END
C***********************************************************************
C $Id: arorie.F,v 0.18 1992/03/12 11:23:19 lonnblad Exp $

      REAL FUNCTION ARECOI(BE0,DE1,DE2,DE3,BP0,DP1,DP2,DP3,BALP,PT12)

C...Ariadne function RECOIl

C...Calculates the angle of a recoil gluon according to the new
C...Recoil strategy: p_t1^2*exp(-y_1)=p_t2^2*exp(-y_2)


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)

      COMMON /LUDAT1/ MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /LUDAT1/

      COMMON /ARDAT1/ PARA(40),MSTA(40)
      SAVE /ARDAT1/

C...Calculate the maximum and minimum angle
      PHIL=0.0
      PHIU=PARU(1)-BALP

C...Calculate angle of recoil gluon
      BW=DE1+DE2+DE3
      BS=BW**2
      BM3=DE3**2-DP3**2
      S0123=(BW+BE0)**2-BP0**2
      S12=BS-2.0*BW*DE3+BM3
      S23=BS-2.0*BW*DE1
      S13=BS-2.0*BW*DE2
      D01=2.0*S12*DE1*BE0
      D02=2.0*S12*DP1*BP0
      D03=PT12*(S0123-S13-S23+BM3-2.0*DE3*BE0)
      D04=2*DP3*BP0*PT12*COS(BALP)
      D05=2*DP3*BP0*PT12*SIN(BALP)
      D11=D01-D03
      D12=D05
      D13=D04+D02
      D21=(D11**2-D13**2)/(D12**2+D13**2)
      D22=D12*D11/(D12**2+D13**2)
      D31=D22**2-D21
      DSPHI=SQRT(MAX(D31,0.0D0))-D22
      IF(DSPHI.LT.0.0D0) THEN
        PHI=PHIL
      ELSEIF(DSPHI.GE.1.0D0) THEN
        PHI=PHIU
      ELSE
        PHI=MIN(ASIN(DSPHI),DBLE(PARU(1))-BALP)
      ENDIF

      ARECOI=PHI

      RETURN


C**** END OF ARECOI ****************************************************
      END
C***********************************************************************
C $Id: arpars.F,v 0.10 1992/03/04 19:36:20 lonnblad Exp $

      SUBROUTINE ARPARS(NSTART,NEND)

C...ARiadne subroutine PARSe the event record

C...Parse through the /LUJETS/ event record to find un-cascaded
C...strings. Performs dipole cascade on each found.


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)

      COMMON /ARPART/ BP(MAXPAR,5),IFL(MAXPAR),IEX(MAXPAR),QQ(MAXPAR),
     $                IDI(MAXPAR),IDO(MAXPAR),INO(MAXPAR),IPART
      SAVE /ARPART/

      COMMON /ARDIPS/ BX1(MAXDIP),BX3(MAXDIP),PT2IN(MAXDIP),
     $                SDIP(MAXDIP),IP1(MAXDIP),IP3(MAXDIP),
     $                AEX1(MAXDIP),AEX3(MAXDIP),QDONE(MAXDIP),
     $                QEM(MAXDIP),IRAD(MAXDIP),ISTR(MAXDIP),IDIPS
      SAVE /ARDIPS/

      COMMON /ARSTRS/ IPF(MAXSTR),IPL(MAXSTR),IFLOW(MAXSTR),
     $                PT2LST,IMF,IML,IO,QDUMP,ISTRS
      SAVE /ARSTRS/

      COMMON /ARDAT1/ PARA(40),MSTA(40)
      SAVE /ARDAT1/

      COMMON /ARJETX/ N,K(300,5),P(300,5),V(300,5)
      SAVE /ARJETX/


      IDIR=0

C...Loop over entries in /LUJETS/ to be considered
      DO 100 I=NSTART,NEND

C...If IDIR=0 there is no current string so skip all entries which
C...are not the begining of a string (K(I,1)=2) otherwise copy
C...parton to dipole record
        IF(IDIR.EQ.0) THEN
          IF(K(I,1).NE.2) GOTO 100
          CALL ARGTYP(I,ITYP)
          IF(ITYP.EQ.0) CALL ARERRM('ARPARa',1,I)
          IDIR=ITYP
          IMF=I
          IPART=1
          IDIPS=0
          CALL ARCOPA(I,IPART,ITYP)
        ELSE

C...If in a string, copy parton and create a dipole. Error if
C...colour singlets of triplets are found
          IF(K(I,1).EQ.2) THEN
            CALL ARGTYP(I,ITYP)
            IF(ABS(ITYP).EQ.1) CALL ARERRM('ARPARS',2,I)
            IF(ABS(ITYP).EQ.0) CALL ARERRM('ARPARb',1,I)
            IPART=IPART+1
            IDIPS=IDIPS+1
            CALL ARCOPA(I,IPART,ITYP)
            CALL ARCRDI(IDIPS,IPART-1,IPART,1,.FALSE.)

C...If the end of a string check colour flow and consistency
          ELSEIF(K(I,1).EQ.1) THEN
            CALL ARGTYP(I,ITYP)
            IF(ITYP.EQ.0) CALL ARERRM('ARPARc',1,I)
            IML=I
            IPART=IPART+1
            IDIPS=IDIPS+1
            CALL ARCOPA(I,IPART,ITYP)
            CALL ARCRDI(IDIPS,IPART-1,IPART,1,.FALSE.)
C...........If purely gluonic string create extra dipole
            IF(ITYP.EQ.2) THEN
              IF(IDIR.NE.2) CALL ARERRM('ARPARS',4,I)
              IDIPS=IDIPS+1
              CALL ARCRDI(IDIPS,IPART,1,1,.FALSE.)
C...........If ordinary string create EM-dipole
            ELSE
              IF(ITYP.NE.-IDIR) CALL ARERRM('ARPARS',5,I)
              IF(MSTA(20).GT.0.AND.IDIPS.EQ.1.AND.
     $               IEX(1).EQ.0.AND.IEX(IPART).EQ.0) THEN
                IDIPS=IDIPS+1
                CALL ARCRDI(IDIPS,IPART,1,1,.TRUE.)
              ENDIF
            ENDIF

C...Initialize string variables in dipole record and perform cascade
            PT2LST=PARA(40)
            IF(MSTA(14).GT.1.AND.IPART.GT.2) PT2LST=ARMIPT(1,IPART)
            IF(PARA(6).GT.0.0) PT2LST=MIN(PT2LST,PARA(6)**2)
C...Don't cascade purelu gluonic string
            IF(IDIR.NE.2) THEN
              IPF(1)=1
              IPL(1)=IPART
              ISTRS=1
              IFLOW(1)=IDIR
              CALL AREXMA(1,IPART)
              CALL ARCASC
              IDIR=0
            ENDIF
          ENDIF
        ENDIF
 100  CONTINUE


      RETURN

C**** END OF ARPARS ****************************************************
      END
C***********************************************************************
C $Id: arradg.F,v 0.15 1992/03/12 12:08:18 lonnblad Exp $

      SUBROUTINE ARRADG(ID,NREM,SNR,PT21,PT23)

C...ARiadne subroutine RADiate Gluon

C...Performs the radiation of a gluon from dipole ID


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)

      COMMON /ARPART/ BP(MAXPAR,5),IFL(MAXPAR),IEX(MAXPAR),QQ(MAXPAR),
     $                IDI(MAXPAR),IDO(MAXPAR),INO(MAXPAR),IPART
      SAVE /ARPART/

      COMMON /ARDIPS/ BX1(MAXDIP),BX3(MAXDIP),PT2IN(MAXDIP),
     $                SDIP(MAXDIP),IP1(MAXDIP),IP3(MAXDIP),
     $                AEX1(MAXDIP),AEX3(MAXDIP),QDONE(MAXDIP),
     $                QEM(MAXDIP),IRAD(MAXDIP),ISTR(MAXDIP),IDIPS
      SAVE /ARDIPS/

      COMMON /ARSTRS/ IPF(MAXSTR),IPL(MAXSTR),IFLOW(MAXSTR),
     $                PT2LST,IMF,IML,IO,QDUMP,ISTRS
      SAVE /ARSTRS/

      COMMON /ARDAT1/ PARA(40),MSTA(40)
      SAVE /ARDAT1/

      COMMON /ARINT2/ DBEX,DBEY,DBEZ,PHI,THE
      SAVE /ARINT2/

      INXT(I)=IDO(IP3(I))


C...Boost dipole to its CMS
      CALL ARBOCM(ID)

C...Copy some information about dipole
      BS=ARMAS2(IP1(ID),IP3(ID))
      IF(ABS(BS-SDIP(ID)).GT.(BS+SDIP(ID))*PARA(39).AND.
     $     MSTA(9).GE.2) CALL ARERRM('ARRADG',13,0)

      BW=SQRT(BS)
      B1=BX1(ID)
      B3=BX3(ID)
      NE1=IEX(IP1(ID))
      NE3=IEX(IP3(ID))

C...If parton not extended - no recoil gluon (trivial)
      IF(NE1.EQ.0) AEX1(ID)=2.0
      IF(NE3.EQ.0) AEX3(ID)=2.0

C...No recoil gluon if reemission
      IF(NREM.EQ.1) AEX1(ID)=2.0
      IF(NREM.EQ.3) AEX3(ID)=2.0

C...If AEX1(3) >= 1 then no recoil gluon
      IF(MSTA(17).EQ.0) THEN
        AEX1(ID)=2.0
        AEX3(ID)=2.0
      ENDIF

C...No recoil gluons if not enough energy left for original parton
      IF(AEX1(ID).LT.1.0.OR.AEX3(ID).LT.1.0) THEN
        BY1=BP(IP1(ID),5)**2/BS
        BY3=BP(IP3(ID),5)**2/BS
        BPT=0.5*SQRT(BS)*
     $       (1.0+BY1-BY3+SQRT(1.0+(BY1-BY3)**2-2.0*(BY1+BY3)))
        B1P=(1.0-AEX1(ID))*BPT
        IF(B1P.LT.BP(IP1(ID),5)) THEN
          AEX1(ID)=2.0
          B1P=0.0
          B1M=0.0
        ELSE
          B1M=BS*BY1/B1P
        ENDIF
        BMT=0.5*SQRT(BS)*
     $       (1.0+BY3-BY1+SQRT(1.0+(BY1-BY3)**2-2.0*(BY1+BY3)))
        B3M=(1.0-AEX3(ID))*BMT
        IF(B3M.LT.BP(IP3(ID),5)) THEN
          AEX3(ID)=2.0
          B3P=0.0
          B3M=0.0
        ELSE
          B3P=BS*BY3/B3M
        ENDIF
      ENDIF

C...Check if any parton can take full recoil.
      QR1=(QQ(IP1(ID)).AND.MSTA(16).GE.1.AND.(IEX(IP1(ID)).EQ.0.OR.
     $     (IEX(IP1(ID)).GE.1.AND.MSTA(16).EQ.2.AND.AEX1(ID).GE.1.0)))
      QR3=(QQ(IP3(ID)).AND.MSTA(16).GE.1.AND.(IEX(IP3(ID)).EQ.0.OR.
     $     (IEX(IP3(ID)).GE.1.AND.MSTA(16).EQ.2.AND.AEX3(ID).GE.1.0)))

C...Special treatment for Drell-Yan produced particles
      IF(MSTA(23).GT.0) CALL ARDYRE(ID,*100)

C...No recoil gluons if one parton can take full recoil
      IF((AEX1(ID).LT.1.0.OR.AEX3(ID).LT.1.0).AND.MSTA(17).EQ.1) THEN
        IF(QR3) AEX1(ID)=2.0
        IF(QR1) AEX3(ID)=2.0
      ENDIF

      QRG1=(AEX1(ID).LT.1.0)
      QRG3=(AEX3(ID).LT.1.0)

      IDE=ID

C...Add recoil gluon for parton 1
      IF(QRG1) THEN
        CALL ARADDG(ID)
        IDE=INXT(ID)
        BP(IP1(ID),1)=0.0
        BP(IP1(ID),2)=0.0
        BP(IP1(ID),3)=0.5*(B1P-B1M)
        BP(IP1(ID),4)=0.5*(B1P+B1M)
        INO(IP3(ID))=-IO
      ENDIF

C...Add emitted gluon
      CALL ARADDG(IDE)
      INO(IP3(IDE))=IO

C...Add recoil gluon for parton 3
      IF(QRG3) THEN
        IDL=INXT(IDE)
        CALL ARADDG(IDL)
        IDL=INXT(IDL)
        BP(IP3(IDL),1)=0.0
        BP(IP3(IDL),2)=0.0
        BP(IP3(IDL),3)=0.5*(B3P-B3M)
        BP(IP3(IDL),4)=0.5*(B3P+B3M)
        INO(IP1(IDL))=-IO
      ENDIF

      IF(NREM.EQ.0) THEN
        IF(QRG1.AND.QRG3) THEN
          SNR3=BS*((BW-B1M)*(1.0-B1+BY1-BY3)/BW+BY3)
          SNR1=BS*((BW-B3P)*(1.0-B3+BY3-BY1)/BW+BY1)
        ELSEIF(QRG1) THEN
          SNR=BS*(1.0-B3+BY3)
        ELSEIF(QRG3) THEN
          SNR=BS*(1.0-B1+BY1)
        ELSE
          SNR=0.0
        ENDIF
      ENDIF

      PT21=0.0
      PT23=0.0
      IF(QRG1.OR.QRG3) THEN
        B2M=(1.0-B3+BY3-BY1)*BW
        B2P=(1.0-B1+BY1-BY3)*BW
        IF(QRG1.AND.MSTA(17).GE.2) PT21=(B2M*B2P**3)/(BW-B1P-B2P)**2
        IF(QRG3.AND.MSTA(17).GE.2) PT23=(B2P*B2M**3)/(BW-B3M-B2M)**2
        DA=(BW-B1P-B3P)/(BW-B1M-B3M)
        SA=(BW-B1P-B3P)*(BW-B1M-B3M)/BS
        DB=(DA-1.0D0)/(DA+1.0D0)
        BY1A=BY1/SA
        IF(QRG1) BY1A=0.0
        BY3A=BY3/SA
        IF(QRG3) BY3A=0.0
        BS=BS*SA
        B1=1.0-(1.0-B1+BY1-BY3)/SQRT(SA*DA)+BY1A-BY3A
        B3=1.0-(1.0-B3+BY3-BY1)/SQRT(SA/DA)+BY3A-BY1A

        IF(QRG1) CALL AROBO1(0.0,0.0,0.0D0,0.0D0,-DB,IP1(ID))
        IF(QRG3) CALL AROBO1(0.0,0.0,0.0D0,0.0D0,-DB,IP3(IDL))
      ENDIF

C...Disable Kleiss orientation if extended partons
      IF(QR1.AND.QR3.AND.NE1+NE3.NE.0) THEN
        QR1=.FALSE.
        QR3=.FALSE.
      ENDIF

C...Orientate the emitted partons
      IF(NREM.EQ.0) THEN
        CALL ARORIE(IP1(IDE),IP3(IDE),IP3(INXT(IDE)),BS,B1,B3,QR1,QR3,
     $       PT21,PT23)
      ELSEIF(NREM.EQ.1) THEN
        QR1=.FALSE.
        QR3=.TRUE.
        CALL ARORIE(IP1(IDE),IP3(IDE),IP3(INXT(IDE)),BS,B1,B3,
     $       QR1,QR3,0.0,0.0)
      ELSEIF(NREM.EQ.3) THEN
        QR1=.TRUE.
        QR3=.FALSE.
        CALL ARORIE(IP1(IDE),IP3(IDE),IP3(INXT(IDE)),BS,B1,B3,
     $       QR1,QR3,0.0,0.0)
      ENDIF
        
C...Boost created dipoles back to original CMS
      IF((.NOT.QRG1).AND.(.NOT.QRG3)) THEN
        CALL AROBO3(THE,PHI,DBEX,DBEY,DBEZ,
     $              IP1(IDE),IP3(IDE),IP3(INXT(IDE)))
      ELSEIF(QRG1.AND.(.NOT.QRG3)) THEN
        IF(MSTA(17).LT.2) PT21=ARIPT2(IP1(ID),IP1(IDE),IP3(IDE))
        CALL AROBO4(0.0,0.0,0.0D0,0.0D0,DB,
     $              IP1(ID),IP1(IDE),IP3(IDE),IP3(INXT(IDE)))
        CALL AROBO4(THE,PHI,DBEX,DBEY,DBEZ,
     $              IP1(ID),IP1(IDE),IP3(IDE),IP3(INXT(IDE)))
      ELSEIF((.NOT.QRG1).AND.QRG3) THEN
        IF(MSTA(17).LT.2) PT23=ARIPT2(IP1(IDE),IP3(IDE),IP3(INXT(IDE)))
        CALL AROBO4(0.0,0.0,0.0D0,0.0D0,DB,
     $              IP1(IDE),IP3(IDE),IP3(INXT(IDE)),IP3(IDL))
        CALL AROBO4(THE,PHI,DBEX,DBEY,DBEZ,
     $              IP1(IDE),IP3(IDE),IP3(INXT(IDE)),IP3(IDL))
      ELSEIF(QRG1.AND.QRG3) THEN
        IF(MSTA(17).LT.2) THEN
          PT21=ARIPT2(IP1(ID),IP1(IDE),IP3(IDE))
          PT23=ARIPT2(IP3(IDE),IP3(INXT(IDE)),IP3(IDL))
        ENDIF
        IF(PT21.GE.PT23) THEN
          SNR=SNR3
        ELSE
          SNR=SNR1
        ENDIF
        CALL AROBO5(0.0,0.0,0.0D0,0.0D0,DB,
     $              IP1(ID),IP1(IDE),IP3(IDE),IP3(INXT(IDE)),IP3(IDL))
        CALL AROBO5(THE,PHI,DBEX,DBEY,DBEZ,
     $              IP1(ID),IP1(IDE),IP3(IDE),IP3(INXT(IDE)),IP3(IDL))
      ENDIF

 100  CONTINUE

      RETURN

C**** END OF ARRADG ****************************************************
      END
C***********************************************************************
C $Id: arradp.F,v 0.9 1992/02/05 14:48:41 lonnblad Exp $

      SUBROUTINE ARRADP(ID)

C...ARiadne subroutine RADiate Photon

C...Performs the radiation of a photon from EM-dipole ID


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)

      COMMON /ARPART/ BP(MAXPAR,5),IFL(MAXPAR),IEX(MAXPAR),QQ(MAXPAR),
     $                IDI(MAXPAR),IDO(MAXPAR),INO(MAXPAR),IPART
      SAVE /ARPART/

      COMMON /ARDIPS/ BX1(MAXDIP),BX3(MAXDIP),PT2IN(MAXDIP),
     $                SDIP(MAXDIP),IP1(MAXDIP),IP3(MAXDIP),
     $                AEX1(MAXDIP),AEX3(MAXDIP),QDONE(MAXDIP),
     $                QEM(MAXDIP),IRAD(MAXDIP),ISTR(MAXDIP),IDIPS
      SAVE /ARDIPS/

      COMMON /ARSTRS/ IPF(MAXSTR),IPL(MAXSTR),IFLOW(MAXSTR),
     $                PT2LST,IMF,IML,IO,QDUMP,ISTRS
      SAVE /ARSTRS/

      COMMON /ARDAT1/ PARA(40),MSTA(40)
      SAVE /ARDAT1/

      COMMON /ARINT2/ DBEX,DBEY,DBEZ,PHI,THE
      SAVE /ARINT2/

      INXT(I)=IDO(IP3(I))
      IPRV(I)=IDI(IP1(I))

C...Boost dipole to its CMS, and get its invaiant mass^2
      CALL ARBOCM(ID)
      BS=ARMAS2(IP1(ID),IP3(ID))
      IF(ABS(BS-SDIP(ID)).GT.(BS+SDIP(ID))*PARA(39).AND.
     $     MSTA(9).GE.2) CALL ARERRM('ARRADG',13,0)

      QR1=.TRUE.
      QR3=.TRUE.
C...Use position IPART+1 temporarily for the photon and orientate
C...the particles/partons
      BP(IPART+1,5)=0.0
      CALL ARORIE(IP1(ID),IPART+1,IP3(ID),BS,BX1(ID),BX3(ID),
     $            QR1,QR3,0.0,0.0)

C...Boost back to original CMS
      CALL AROBO3(THE,PHI,DBEX,DBEY,DBEZ,
     $            IP1(ID),IPART+1,IP3(ID))
C...Copy photon information to /LUJETS/
      CALL ARDUPH

C...Flagg dipoles that were affected by the emission
      QDONE(INXT(ID))=.FALSE.
      QDONE(IPRV(ID))=.FALSE.
      QDONE(ID)=.FALSE.

      RETURN

C**** END OF ARRADP ****************************************************
      END
C***********************************************************************
C $Id: arradq.F,v 0.8 1992/01/31 16:14:59 lonnblad Exp $

      SUBROUTINE ARRADQ(ID)

C...ARiadne subroutine RADiate Q-qbar pair

C...Performs the emission of a q-qbar pair from gluon in dipole ID


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)

      COMMON /ARPART/ BP(MAXPAR,5),IFL(MAXPAR),IEX(MAXPAR),QQ(MAXPAR),
     $                IDI(MAXPAR),IDO(MAXPAR),INO(MAXPAR),IPART
      SAVE /ARPART/

      COMMON /ARDIPS/ BX1(MAXDIP),BX3(MAXDIP),PT2IN(MAXDIP),
     $                SDIP(MAXDIP),IP1(MAXDIP),IP3(MAXDIP),
     $                AEX1(MAXDIP),AEX3(MAXDIP),QDONE(MAXDIP),
     $                QEM(MAXDIP),IRAD(MAXDIP),ISTR(MAXDIP),IDIPS
      SAVE /ARDIPS/

      COMMON /ARSTRS/ IPF(MAXSTR),IPL(MAXSTR),IFLOW(MAXSTR),
     $                PT2LST,IMF,IML,IO,QDUMP,ISTRS
      SAVE /ARSTRS/

      COMMON /ARDAT1/ PARA(40),MSTA(40)
      SAVE /ARDAT1/

      COMMON /ARINT2/ DBEX,DBEY,DBEZ,PHI,THE
      SAVE /ARINT2/

      INXT(I)=IDO(IP3(I))
      IPRV(I)=IDI(IP1(I))

C...Boost dipole to its CMS and copy its invariant mass^2
      CALL ARBOCM(ID)
      BS=ARMAS2(IP1(ID),IP3(ID))
      IF(ABS(BS-SDIP(ID)).GT.(BS+SDIP(ID))*PARA(39).AND.
     $     MSTA(9).GE.2) CALL ARERRM('ARRADQ',13,0)

C...Check which gluon to split
      IF(IRAD(ID).LT.0) THEN
C.......Determine patons ability to recoil, save pointers and flag
C.......affected dipoles
        QR1=.TRUE.
        QR3=(QQ(IP3(ID)).AND.MSTA(16).GT.0)
        IPG=IP1(ID)
        IDN=ID
        IDP=IPRV(ID)
        IF(INXT(ID).NE.0) QDONE(INXT(ID))=.FALSE.
C.......Split the gluon entry, orientate the partons, and boost back
        CALL ARSPLG(IPG,ABS(IRAD(ID)))
        CALL ARORIE(IP3(IDP),IP1(IDN),IP3(IDN),BS,BX1(ID),BX3(ID),
     $              QR1,QR3,0.0,0.0)
        CALL AROBO3(THE,PHI,DBEX,DBEY,DBEZ,
     $              IP3(IDP),IP1(IDN),IP3(IDN))
      ELSE
C.......Determine patons ability to recoil, save pointers and flag
C.......affected dipoles
        QR3=.TRUE.
        QR1=(QQ(IP1(ID)).AND.MSTA(16).GT.0)
        IPG=IP3(ID)
        IDP=ID
        IDN=INXT(ID)
        IF(IPRV(ID).NE.0) QDONE(IPRV(ID))=.FALSE.
C.......Split the gluon entry, orientate the partons, and boost back
        CALL ARSPLG(IPG,ABS(IRAD(ID)))
        CALL ARORIE(IP1(IDP),IP3(IDP),IP1(IDN),BS,BX1(ID),BX3(ID),
     $              QR1,QR3,0.0,0.0)
        CALL AROBO3(THE,PHI,DBEX,DBEY,DBEZ,
     $              IP1(IDP),IP3(IDP),IP1(IDN))
      ENDIF



      RETURN

C**** END OF ARRADQ ****************************************************
      END
C***********************************************************************
C $Id: arreca.F,v 0.4 1991/11/06 14:14:58 lonnblad Exp $

      SUBROUTINE ARRECA(ID,IDS,IS1,IS3)

C...ARiadne function RECAll

C...Recalls a dipole entry stored by ARSTOR


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)

      COMMON /ARPART/ BP(MAXPAR,5),IFL(MAXPAR),IEX(MAXPAR),QQ(MAXPAR),
     $                IDI(MAXPAR),IDO(MAXPAR),INO(MAXPAR),IPART
      SAVE /ARPART/

      COMMON /ARDIPS/ BX1(MAXDIP),BX3(MAXDIP),PT2IN(MAXDIP),
     $                SDIP(MAXDIP),IP1(MAXDIP),IP3(MAXDIP),
     $                AEX1(MAXDIP),AEX3(MAXDIP),QDONE(MAXDIP),
     $                QEM(MAXDIP),IRAD(MAXDIP),ISTR(MAXDIP),IDIPS
      SAVE /ARDIPS/

      COMMON /ARSTRS/ IPF(MAXSTR),IPL(MAXSTR),IFLOW(MAXSTR),
     $                PT2LST,IMF,IML,IO,QDUMP,ISTRS
      SAVE /ARSTRS/



      SDIP(ID)=SDIP(IDS)
      IP1(ID)=IP1(IDS)
      IP3(ID)=IP3(IDS)
      BX1(ID)=BX1(IDS)
      BX3(ID)=BX3(IDS)
      PT2IN(ID)=PT2IN(IDS)
      AEX1(ID)=AEX1(IDS)
      AEX3(ID)=AEX3(IDS)
      QDONE(ID)=QDONE(IDS)
      QEM(ID)=QEM(IDS)
      IRAD(ID)=IRAD(IDS)
      ISTR(ID)=ISTR(IDS)

      I1=IP1(ID)
      I3=IP3(ID)

      DO 100 I=1,5
        BP(I1,I)=BP(IS1,I)
        BP(I3,I)=BP(IS3,I)
 100  CONTINUE
      IFL(I1)=IFL(IS1)
      IFL(I3)=IFL(IS3)
      IEX(I1)=IEX(IS1)
      IEX(I3)=IEX(IS3)
      QQ(I1)=QQ(IS1)
      QQ(I3)=QQ(IS3)
      IDI(I1)=IDI(IS1)
      IDI(I3)=IDI(IS3)
      IDO(I1)=IDO(IS1)
      IDO(I3)=IDO(IS3)
      INO(I1)=INO(IS1)
      INO(I3)=INO(IS3)

      RETURN

C**** END OF ARRECA ****************************************************
      END
C***********************************************************************
C $Id: arupdj.F,v 0.1 1992/01/27 16:03:19 lonnblad Exp $

      SUBROUTINE ARUPDJ(I2,I1,I3)

C...ARiadne subroutine UPDate Jet entry

C...Takes a jet entry I2 and determines its minimum invariant pt wrt.
C...all other jets. I1 and I3 indicates which jets have been changed
C...since last call


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)

      COMMON /ARDAT1/ PARA(40),MSTA(40)
      SAVE /ARDAT1/

      COMMON /ARJETX/ N,K(300,5),P(300,5),V(300,5)
      SAVE /ARJETX/

      COMMON /LUDAT1/ MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /LUDAT1/


      IF(K(I2,5).LT.0) RETURN
      IF(I1.EQ.0.OR.K(K(I2,3),5).LT.0.OR.K(K(I2,4),5).LT.0.OR.
     $     I2.EQ.I1.OR.I2.EQ.I3.OR.K(I2,3).EQ.I1.OR.K(I2,4).EQ.I3) THEN
        V(I2,5)=PARA(40)
        DO 100 J1=N+1,N+MSTU(3)-1
          IF(K(J1,5).LT.0) GOTO 100
          IF(J1.EQ.I2) GOTO 100
          DO 110 J3=J1+1,N+MSTU(3)
            IF(K(J3,5).LT.0) GOTO 110
            IF(J3.EQ.I2) GOTO 110
            CALL ARSMPT(J1,I2,J3)
 110      CONTINUE
 100    CONTINUE
      ELSE
        DO 200 J=N+1,N+MSTU(3)
          IF(J.EQ.I2.OR.K(J,5).LT.0) GOTO 200
          IF(J.GT.I1) CALL ARSMPT(I1,I2,J)
          IF(J.LT.I3) CALL ARSMPT(J,I2,I3)
          IF(J.LT.I1) CALL ARSMPT(J,I2,I1)
          IF(J.GT.I3) CALL ARSMPT(I3,I2,J)
 200    CONTINUE
      ENDIF

      RETURN

C**** END OF ARUPDJ ****************************************************
      END

      SUBROUTINE ARSMPT(I1,I2,I3)


      COMMON /ARJETX/ N,K(300,5),P(300,5),V(300,5)
      SAVE /ARJETX/

      P12=P(I1,4)*P(I2,4)-P(I1,3)*P(I2,3)-
     $    P(I1,2)*P(I2,2)-P(I1,1)*P(I2,1)
      P23=P(I3,4)*P(I2,4)-P(I3,3)*P(I2,3)-
     $    P(I3,2)*P(I2,2)-P(I3,1)*P(I2,1)
      P31=P(I1,4)*P(I3,4)-P(I1,3)*P(I3,3)-
     $    P(I1,2)*P(I3,2)-P(I1,1)*P(I3,1)
      PT2I=(V(I1,1)+V(I2,1)+P12)*(V(I2,1)+V(I3,1)+P23)/
     $     (V(I1,1)+V(I2,1)+V(I3,1)+P12+P23+P31)

      IF(PT2I.GE.V(I2,5)) RETURN
      V(I2,5)=PT2I
      K(I2,3)=I1
      K(I2,4)=I3  

      RETURN
      END
C***********************************************************************
C $Id: arrndx.F,v 0.4 1992/01/31 16:14:59 lonnblad Exp $

      REAL FUNCTION ARNDX1()

C...Ariadne function RNDom Xt2 version 1

C...Generate an x_t^2 according to a Sudakov suppressed distribution.
C...Suitable for running alpha_QCD


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)

      double precision PYR

      COMMON /ARDAT1/ PARA(40),MSTA(40)
      SAVE /ARDAT1/

      COMMON /ARINT1/ BC1,BC3,BZM,BZP,BP1,BM1,BP3,BM3,
     $                B1,B2,B3,XT2,XT,Y,QQ1,QQ3,NE1,NE3,
     $                S,W,C,CN,ALPHA0,XLAM2,IFLG,
     $                XT2MP,XT2ME,XT2M,XT2C,XTS,XT3,XT1,
     $                YINT,YMAX,YMIN,
     $                Y1,Y2,Y3,SY1,SY2,SY3,SSY,
     $                AE1,AE3,NXP1,NXP3,FQ1,FQ3
      SAVE /ARINT1/


      ARNDX1=0.0
      ARG=PYR(IDUM)
      IF(LOG(ARG)*CN.LT.LOG(LOG(XT2C/XLAM2)/LOG(XT2M/XLAM2))) RETURN
      ARNDX1=XLAM2*(XT2M/XLAM2)**(ARG**CN)

      RETURN

C**** END OF ARNDX1 ****************************************************
      END
C***********************************************************************

      REAL FUNCTION ARNDX2()

C...Ariadne function RNDom Xt2 version 2

C...Generate an x_t^2 according to a Sudakov suppressed distribution.
C...Suitable for constant alpha_QCD and QED emission


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)

      double precision PYR

      COMMON /ARINT1/ BC1,BC3,BZM,BZP,BP1,BM1,BP3,BM3,
     $                B1,B2,B3,XT2,XT,Y,QQ1,QQ3,NE1,NE3,
     $                S,W,C,CN,ALPHA0,XLAM2,IFLG,
     $                XT2MP,XT2ME,XT2M,XT2C,XTS,XT3,XT1,
     $                YINT,YMAX,YMIN,
     $                Y1,Y2,Y3,SY1,SY2,SY3,SSY,
     $                AE1,AE3,NXP1,NXP3,FQ1,FQ3
      SAVE /ARINT1/


      ARNDX2=0.0
      ARG=PYR(IDUM)
      IF(CN*LOG(ARG).LT.(LOG(XT2M))**2-(LOG(XT2C))**2) RETURN
      ARNDX2=EXP(-SQRT((LOG(XT2M))**2-LOG(ARG)*CN))

      RETURN

C**** END OF ARNDX2 ****************************************************
      END
C***********************************************************************

      REAL FUNCTION ARNDX3()

C...Ariadne function RNDom Xt2 version 3

C...Generate an x_t^2 according to a Sudakov suppressed distribution.
C...Suitable for constant alpha_QCD q-qbar emission


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)

      double precision PYR

      COMMON /ARDAT1/ PARA(40),MSTA(40)
      SAVE /ARDAT1/

      COMMON /ARINT1/ BC1,BC3,BZM,BZP,BP1,BM1,BP3,BM3,
     $                B1,B2,B3,XT2,XT,Y,QQ1,QQ3,NE1,NE3,
     $                S,W,C,CN,ALPHA0,XLAM2,IFLG,
     $                XT2MP,XT2ME,XT2M,XT2C,XTS,XT3,XT1,
     $                YINT,YMAX,YMIN,
     $                Y1,Y2,Y3,SY1,SY2,SY3,SSY,
     $                AE1,AE3,NXP1,NXP3,FQ1,FQ3
      SAVE /ARINT1/


      ARNDX3=0.0
      ARG=PYR(IDUM)
      IF(LOG(ARG)*CN.LT.LOG(XT2C/XT2M)) RETURN
      ARNDX3=XT2M*(ARG**CN)

      RETURN

C**** END OF ARNDX3 ****************************************************
      END
C***********************************************************************
C $Id: arrndy.F,v 0.6 1992/02/03 09:48:57 lonnblad Exp $

      REAL FUNCTION ARNDY1()

C...Ariadne function RaNDom Y version 1

C...Generates a properly distributed Y
C...Suitable for gluon and photon emission


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)

      double precision PYR

      COMMON /ARINT1/ BC1,BC3,BZM,BZP,BP1,BM1,BP3,BM3,
     $                B1,B2,B3,XT2,XT,Y,QQ1,QQ3,NE1,NE3,
     $                S,W,C,CN,ALPHA0,XLAM2,IFLG,
     $                XT2MP,XT2ME,XT2M,XT2C,XTS,XT3,XT1,
     $                YINT,YMAX,YMIN,
     $                Y1,Y2,Y3,SY1,SY2,SY3,SSY,
     $                AE1,AE3,NXP1,NXP3,FQ1,FQ3
      SAVE /ARINT1/


      ZMAX=SQRT(XTS/XT2)+SQRT(MAX(XTS/XT2-1.0,0.0))
      YMAX=LOG(MIN(ZMAX,XT3/XT))
      YMIN=-LOG(MIN(ZMAX,XT1/XT))

      ARNDY1=YMIN+PYR(IDUM)*(YMAX-YMIN)

      RETURN

C**** END OF ARNDY1 ****************************************************
      END
C***********************************************************************

      REAL FUNCTION ARNDY2()

C...Ariadne function RaNDom Y version 2

C...Generates a properly distributed Y
C...Suitable for gluon emission from extended dipole


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)

      double precision PYR

      COMMON /ARDAT1/ PARA(40),MSTA(40)
      SAVE /ARDAT1/

      COMMON /ARINT1/ BC1,BC3,BZM,BZP,BP1,BM1,BP3,BM3,
     $                B1,B2,B3,XT2,XT,Y,QQ1,QQ3,NE1,NE3,
     $                S,W,C,CN,ALPHA0,XLAM2,IFLG,
     $                XT2MP,XT2ME,XT2M,XT2C,XTS,XT3,XT1,
     $                YINT,YMAX,YMIN,
     $                Y1,Y2,Y3,SY1,SY2,SY3,SSY,
     $                AE1,AE3,NXP1,NXP3,FQ1,FQ3
      SAVE /ARINT1/


      ZMAX=SQRT(XTS/XT2)+SQRT(MAX(XTS/XT2-1.0,0.0))

      AE1=1.0
      AE3=1.0
      IF(NE1.GT.0) AE1=(PARA(10+NE1)/(XT*W))**PARA(10)
      IF(NE3.GT.0) AE3=(PARA(10+NE3)/(XT*W))**PARA(10)
      BP1=(1.0-AE1)*BZP
      IF(BP1.LE.SY1) THEN
        BP1=0.0
        BM1=0.0
      ELSE
        BM1=Y1/BP1
      ENDIF
      BM3=(1.0-AE3)*BZM
      IF(BM3.LE.SY3) THEN
        BM3=0.0
        BP3=0.0
      ELSE
        BP3=Y3/BM3
      ENDIF
      AZ1=1.0-BP1-BP3
      AZ3=1.0-BM1-BM3
      A=(0.5+SQRT(MAX(0.25-XT2/(AZ1*AZ3),0.0)))/XT

      YMAX=LOG(MIN(ZMAX,MIN(XT3/XT,ABS(AZ1)*A)))
      YMIN=-LOG(MIN(ZMAX,MIN(XT1/XT,ABS(AZ3)*A)))

      ARNDY2=YMIN+PYR(IDUM)*(YMAX-YMIN)

      RETURN

C**** END OF ARNDY2 ****************************************************
      END
C***********************************************************************

      REAL FUNCTION ARNDY3()

C...Ariadne function RaNDom Y version 3

C...Generates a properly distributed Y
C...Suitable for q-qbar emission


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)

      double precision PYR

      COMMON /ARINT1/ BC1,BC3,BZM,BZP,BP1,BM1,BP3,BM3,
     $                B1,B2,B3,XT2,XT,Y,QQ1,QQ3,NE1,NE3,
     $                S,W,C,CN,ALPHA0,XLAM2,IFLG,
     $                XT2MP,XT2ME,XT2M,XT2C,XTS,XT3,XT1,
     $                YINT,YMAX,YMIN,
     $                Y1,Y2,Y3,SY1,SY2,SY3,SSY,
     $                AE1,AE3,NXP1,NXP3,FQ1,FQ3
      SAVE /ARINT1/


      ZMAX=SQRT(XTS/XT2)+SQRT(MAX(XTS/XT2-1.0,0.0))
      ZMIN=MIN(ZMAX,XT1/XT)
      ZMAX=MIN(ZMAX,XT3/XT)

      YMAX=LOG(ZMAX)
      YMIN=-LOG(ZMIN)

      ARNDY3=-LOG(1.0/ZMAX+PYR(IDUM)*(ZMIN-1.0/ZMAX))

      RETURN

C**** END OF ARNDY3 ****************************************************
      END
C***********************************************************************

      REAL FUNCTION ARNDY4()

C...Ariadne function RaNDom Y version 4

C...Generates a properly distributed Y
C...Suitable for q-qbar emission from extended dipole


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)

      double precision PYR

      COMMON /ARINT1/ BC1,BC3,BZM,BZP,BP1,BM1,BP3,BM3,
     $                B1,B2,B3,XT2,XT,Y,QQ1,QQ3,NE1,NE3,
     $                S,W,C,CN,ALPHA0,XLAM2,IFLG,
     $                XT2MP,XT2ME,XT2M,XT2C,XTS,XT3,XT1,
     $                YINT,YMAX,YMIN,
     $                Y1,Y2,Y3,SY1,SY2,SY3,SSY,
     $                AE1,AE3,NXP1,NXP3,FQ1,FQ3
      SAVE /ARINT1/

      COMMON /ARDAT1/ PARA(40),MSTA(40)
      SAVE /ARDAT1/


      ZMAX=SQRT(XTS/XT2)+SQRT(MAX(XTS/XT2-1.0,0.0))
      ZMIN=MIN(ZMAX,XT1/XT)
      ZMAX=MIN(ZMAX,XT3/XT)

      AE1=1.0
      AE3=1.0
      IF(NE1.GT.0) AE1=(PARA(10+NE1)/(XT*W))**PARA(10)
      IF(NE3.GT.0) AE3=(PARA(10+NE3)/(XT*W))**PARA(10)
      BP1=(1.0-AE1)*BZP
      IF(BP1.LE.SY1) THEN
        BP1=0.0
        BM1=0.0
      ELSE
        BM1=Y1/BP1
      ENDIF
      BM3=(1.0-AE3)*BZM
      IF(BM3.LE.SY3) THEN
        BM3=0.0
        BP3=0.0
      ELSE
        BP3=Y3/BM3
      ENDIF
      AZ1=1.0-BP1-BP3
      AZ3=1.0-BM1-BM3
      A=(0.5+SQRT(MAX(0.25-XT2/(AZ1*AZ3),0.0)))/XT

      ZMAX=MIN(ZMAX,ABS(AZ1)*A)
      ZMIN=MIN(ZMIN,ABS(AZ3)*A)

      YMAX=LOG(ZMAX)
      YMIN=-LOG(ZMIN)

      ARNDY4=-LOG(1.0/ZMAX+PYR(IDUM)*(ZMIN-1.0/ZMAX))

      RETURN

C**** END OF ARNDY4 ****************************************************
      END
C***********************************************************************
C $Id: arrobo.F,v 0.5 1992/03/12 12:10:13 lonnblad Exp $

      SUBROUTINE AROBO1(THE,PHI,DBEX,DBEY,DBEZ,I1)

C...Ariadne subroutine ROtate BOost 1 parton

C...Rotates and boosts 1 parton in /ARPART/


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)
      DIMENSION I(1)


      I(1)=I1
      CALL ARROBO(THE,PHI,DBEX,DBEY,DBEZ,1,I)

      RETURN

C**** END OF AROBO1 ****************************************************
      END
C***********************************************************************
C $Id: arrobo.F,v 0.5 1992/03/12 12:10:13 lonnblad Exp $

      SUBROUTINE AROBO2(THE,PHI,DBEX,DBEY,DBEZ,I1,I2)

C...Ariadne subroutine ROtate BOost 2 partons

C...Rotates and boosts 2 partons in /ARPART/


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)
      DIMENSION I(2)


      I(1)=I1
      I(2)=I2
      CALL ARROBO(THE,PHI,DBEX,DBEY,DBEZ,2,I)

      RETURN

C**** END OF AROBO2 ****************************************************
      END
C***********************************************************************
C $Id: arrobo.F,v 0.5 1992/03/12 12:10:13 lonnblad Exp $

      SUBROUTINE AROBO3(THE,PHI,DBEX,DBEY,DBEZ,I1,I2,I3)

C...Ariadne subroutine ROtate BOost 3 partons

C...Rotates and boosts 3 partons in /ARPART/


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)
      DIMENSION I(3)


      I(1)=I1
      I(2)=I2
      I(3)=I3
      CALL ARROBO(THE,PHI,DBEX,DBEY,DBEZ,3,I)

      RETURN

C**** END OF AROBO3 ****************************************************
      END
C***********************************************************************
C $Id: arrobo.F,v 0.5 1992/03/12 12:10:13 lonnblad Exp $

      SUBROUTINE AROBO4(THE,PHI,DBEX,DBEY,DBEZ,I1,I2,I3,I4)

C...Ariadne subroutine ROtate BOost 4 partons

C...Rotates and boosts 4 partons in /ARPART/


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)
      DIMENSION I(4)


      I(1)=I1
      I(2)=I2
      I(3)=I3
      I(4)=I4
      CALL ARROBO(THE,PHI,DBEX,DBEY,DBEZ,4,I)

      RETURN

C**** END OF AROBO4 ****************************************************
      END
C***********************************************************************
C $Id: arrobo.F,v 0.5 1992/03/12 12:10:13 lonnblad Exp $

      SUBROUTINE AROBO5(THE,PHI,DBEX,DBEY,DBEZ,I1,I2,I3,I4,I5)

C...Ariadne subroutine ROtate BOost 5 partons

C...Rotates and boosts 5 partons in /ARPART/


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)
      DIMENSION I(5)


      I(1)=I1
      I(2)=I2
      I(3)=I3
      I(4)=I4
      I(5)=I5
      CALL ARROBO(THE,PHI,DBEX,DBEY,DBEZ,5,I)

      RETURN

C**** END OF AROBO5 ****************************************************
      END
C***********************************************************************
C $Id: arrobo.F,v 0.5 1992/03/12 12:10:13 lonnblad Exp $

      SUBROUTINE ARROBO(THE,PHI,DBEX,DBEY,DBEZ,NI,I)

C...ARiadne subroutine ROtate and BOost

C...Rotates and boost NI particles in /ARPART/


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)

      COMMON /ARPART/ BP(MAXPAR,5),IFL(MAXPAR),IEX(MAXPAR),QQ(MAXPAR),
     $                IDI(MAXPAR),IDO(MAXPAR),INO(MAXPAR),IPART
      SAVE /ARPART/

      DIMENSION I(NI),BR(3,3),BV(3),DP(4)


      IF(THE**2+PHI**2.GT.1.0E-20) THEN

C...Rotate (typically from z axis to direction theta,phi)

        SP=SIN(PHI)
        CP=COS(PHI)
        ST=SIN(THE)
        CT=COS(THE)

        BR(1,1)=CT*CP
        BR(1,2)=-SP
        BR(1,3)=ST*CP
        BR(2,1)=CT*SP
        BR(2,2)=CP
        BR(2,3)=ST*SP
        BR(3,1)=-ST
        BR(3,2)=0.0
        BR(3,3)=CT

        DO 100 IJ=1,NI
          DO 110 J=1,3
            BV(J)=BP(I(IJ),J)
 110      CONTINUE
          DO 120 J=1,3
            BP(I(IJ),J)=BR(J,1)*BV(1)+BR(J,2)*BV(2)+BR(J,3)*BV(3)
 120      CONTINUE
 100    CONTINUE

      ENDIF

      DBTOT2=DBEX**2+DBEY**2+DBEZ**2
      IF(DBTOT2.GT.1.0D-20) THEN
        IF(DBTOT2.GE.1.0D0) CALL ARERRM('ARROBO',14,0)
        DGA=1.0D0/DSQRT(1.0D0-DBTOT2)

        DO 200 IJ=1,NI
          DO 210 J=1,4
            DP(J)=BP(I(IJ),J)
 210      CONTINUE
          DBEP=DBEX*DP(1)+DBEY*DP(2)+DBEZ*DP(3)
          DGABEP=DGA*(DGA*DBEP/(1.0D0+DGA)+DP(4))

          BP(I(IJ),1)=DP(1)+DGABEP*DBEX
          BP(I(IJ),2)=DP(2)+DGABEP*DBEY
          BP(I(IJ),3)=DP(3)+DGABEP*DBEZ
          BP(I(IJ),4)=DGA*(DP(4)+DBEP)

 200    CONTINUE

      ENDIF

      RETURN

C**** END OF ARROBO ****************************************************
      END
C***********************************************************************
C $Id: arsplg.F,v 0.10 1992/03/04 19:36:57 lonnblad Exp $

      SUBROUTINE ARSPLG(IG,IFLAV)

C...ARiadne subroutine SPLit Gluon

C...Splits a gluon entry into a q and a q-bar entry


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)

      double precision PYR

      COMMON /ARPART/ BP(MAXPAR,5),IFL(MAXPAR),IEX(MAXPAR),QQ(MAXPAR),
     $                IDI(MAXPAR),IDO(MAXPAR),INO(MAXPAR),IPART
      SAVE /ARPART/

      COMMON /ARDIPS/ BX1(MAXDIP),BX3(MAXDIP),PT2IN(MAXDIP),
     $                SDIP(MAXDIP),IP1(MAXDIP),IP3(MAXDIP),
     $                AEX1(MAXDIP),AEX3(MAXDIP),QDONE(MAXDIP),
     $                QEM(MAXDIP),IRAD(MAXDIP),ISTR(MAXDIP),IDIPS
      SAVE /ARDIPS/

      COMMON /ARSTRS/ IPF(MAXSTR),IPL(MAXSTR),IFLOW(MAXSTR),
     $                PT2LST,IMF,IML,IO,QDUMP,ISTRS
      SAVE /ARSTRS/

      COMMON /ARDAT2/ PQMAS(10)
      SAVE /ARDAT2/

      INXT(I)=IDO(IP3(I))


C...Allocate space for new parton and new string if there is room
      IPART=IPART+1
      ISTRS=ISTRS+1

      IF(IPART.GE.MAXPAR-1) CALL ARERRM('ARSPLG',6,0)
      IF(ISTRS.GT.MAXSTR) CALL ARERRM('ARSPLG',8,0)

C...Set new pointers
      IDP=IDI(IG)
      IDN=IDO(IG)
      IDO(IG)=0
      IDI(IPART)=0
      IDO(IPART)=IDN

      IP1(IDN)=IPART

      IS=ISTR(IDP)

C...If closed gluonic string, no new string is created. The colour flow
C...which was previously undefined is set randomly
      IF(IFLOW(IS).EQ.2) THEN
        ISTRS=ISTRS-1
        IFLOW(IS)=1
        IPF(IS)=IPART
        IPL(IS)=IG
        IF(PYR(IDUM).GT.0.5) IFLOW(IS)=-1
        IFL(IG)=IFLAV*IFLOW(IS)
        IFL(IPART)=-IFL(IG)

C...If new string is created set pointers for its partons
      ELSE
        IFLOW(ISTRS)=IFLOW(IS)
        IPF(ISTRS)=IPART
        IPL(ISTRS)=IPL(IS)
        IPL(IS)=IG
        IFL(IPART)=IFLAV*IFLOW(IS)
        IFL(IG)=-IFL(IPART)
        IDNI=IDN
 100    ISTR(IDNI)=ISTRS
        IF(.NOT.QQ(IP3(IDNI))) THEN
          IDNI=INXT(IDNI)
          GOTO 100
        ENDIF
      ENDIF

C...Reset momenta for created quarks and flag affected dipoles
      DO 200 I=1,4
        BP(IG,I)=0.0
        BP(IPART,I)=0.0
 200  CONTINUE
      BP(IG,5)=PQMAS(IFLAV)
      BP(IPART,5)=PQMAS(IFLAV)
      IEX(IG)=0
      IEX(IPART)=0
      QQ(IG)=.TRUE.
      QQ(IPART)=.TRUE.
      QDONE(IDP)=.FALSE.
      QDONE(IDN)=.FALSE.
      INO(IG)=SIGN(1000*ABS(INO(IG))+IO,INO(IG))
      INO(IPART)=INO(IG)

      RETURN

C**** END OF ARSPLG ****************************************************
      END
C***********************************************************************
C $Id: arstor.F,v 0.3 1991/11/06 14:15:07 lonnblad Exp $

      SUBROUTINE ARSTOR(ID,IDS,IS1,IS3)

C...ARiadne subroutine STORe 

C...Stores a dipole entry for later use


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)

      COMMON /ARPART/ BP(MAXPAR,5),IFL(MAXPAR),IEX(MAXPAR),QQ(MAXPAR),
     $                IDI(MAXPAR),IDO(MAXPAR),INO(MAXPAR),IPART
      SAVE /ARPART/

      COMMON /ARDIPS/ BX1(MAXDIP),BX3(MAXDIP),PT2IN(MAXDIP),
     $                SDIP(MAXDIP),IP1(MAXDIP),IP3(MAXDIP),
     $                AEX1(MAXDIP),AEX3(MAXDIP),QDONE(MAXDIP),
     $                QEM(MAXDIP),IRAD(MAXDIP),ISTR(MAXDIP),IDIPS
      SAVE /ARDIPS/

      COMMON /ARSTRS/ IPF(MAXSTR),IPL(MAXSTR),IFLOW(MAXSTR),
     $                PT2LST,IMF,IML,IO,QDUMP,ISTRS
      SAVE /ARSTRS/

      COMMON /ARINT2/ DBEX,DBEY,DBEZ,PHI,THE
      SAVE /ARINT2/


      IDS=MAXDIP
      SDIP(IDS)=SDIP(ID)
      IP1(IDS)=IP1(ID)
      IP3(IDS)=IP3(ID)
      BX1(IDS)=BX1(ID)
      BX3(IDS)=BX3(ID)
      PT2IN(IDS)=PT2IN(ID)
      AEX1(IDS)=AEX1(ID)
      AEX3(IDS)=AEX3(ID)
      QDONE(IDS)=QDONE(ID)
      QEM(IDS)=QEM(ID)
      IRAD(IDS)=IRAD(ID)
      ISTR(IDS)=ISTR(ID)

      I1=IP1(ID)
      I3=IP3(ID)
      IS1=MAXPAR-1
      IS3=MAXPAR
      DO 100 I=1,5
        BP(IS1,I)=BP(I1,I)
        BP(IS3,I)=BP(I3,I)
 100  CONTINUE
      IFL(IS1)=IFL(I1)
      IFL(IS3)=IFL(I3)
      IEX(IS1)=IEX(I1)
      IEX(IS3)=IEX(I3)
      QQ(IS1)=QQ(I1)
      QQ(IS3)=QQ(I3)
      IDI(IS1)=IDI(I1)
      IDI(IS3)=IDI(I3)
      IDO(IS1)=IDO(I1)
      IDO(IS3)=IDO(I3)
      INO(IS1)=INO(I1)
      INO(IS3)=INO(I3)

      RETURN

C**** END OF ARSTOR ****************************************************
      END
C***********************************************************************
C $Id: artest.F,v 0.12 1992/03/11 09:36:15 lonnblad Exp $

      SUBROUTINE ARTEST(IPRINT)

C...ARiadne subroutine TEST

C...Performs various tests on Ariadne


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)

      double precision PYR

      COMMON /ARPART/ BP(MAXPAR,5),IFL(MAXPAR),IEX(MAXPAR),QQ(MAXPAR),
     $                IDI(MAXPAR),IDO(MAXPAR),INO(MAXPAR),IPART
      SAVE /ARPART/

      COMMON /ARDIPS/ BX1(MAXDIP),BX3(MAXDIP),PT2IN(MAXDIP),
     $                SDIP(MAXDIP),IP1(MAXDIP),IP3(MAXDIP),
     $                AEX1(MAXDIP),AEX3(MAXDIP),QDONE(MAXDIP),
     $                QEM(MAXDIP),IRAD(MAXDIP),ISTR(MAXDIP),IDIPS
      SAVE /ARDIPS/

      COMMON /ARSTRS/ IPF(MAXSTR),IPL(MAXSTR),IFLOW(MAXSTR),
     $                PT2LST,IMF,IML,IO,QDUMP,ISTRS
      SAVE /ARSTRS/

      COMMON /ARDAT1/ PARA(40),MSTA(40)
      SAVE /ARDAT1/

      COMMON /ARDAT2/ PQMAS(10)
      SAVE /ARDAT2/

      COMMON /ARDAT3/ IWRN(40)
      SAVE /ARDAT3/

      COMMON /ARINT1/ BC1,BC3,BZM,BZP,BP1,BM1,BP3,BM3,
     $                B1,B2,B3,XT2,XT,Y,QQ1,QQ3,NE1,NE3,
     $                S,W,C,CN,ALPHA0,XLAM2,IFLG,
     $                XT2MP,XT2ME,XT2M,XT2C,XTS,XT3,XT1,
     $                YINT,YMAX,YMIN,
     $                Y1,Y2,Y3,SY1,SY2,SY3,SSY,
     $                AE1,AE3,NXP1,NXP3,FQ1,FQ3
      SAVE /ARINT1/

      COMMON /ARINT2/ DBEX,DBEY,DBEZ,PHI,THE
      SAVE /ARINT2/

      COMMON /ARJETX/ N,K(300,5),P(300,5),V(300,5)
      SAVE /ARJETX/

      COMMON /LUDAT1/ MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /LUDAT1/

      MSTA(9)=1
      MSTA(6)=-1
      MSTA(20)=1

      MSTJ(21)=0

      PARA(12)=2.0
      PARA(13)=10.0

      CALL ARINIT('ARIADNE')

      DO 110 I=1,10000

        W=10.0*EXP(PYR(IDUM)*LOG(1000.0))
 100    SM1=PYR(IDUM)*20.0
        SM2=PYR(IDUM)*20.0
        E1=0.5*(W**2+SM1**2-SM2**2)/W
        E2=W-E1
        IF(E1.LT.SM1) GOTO 100
        IF(E2.LT.SM2) GOTO 100
        NE1=INT(PYR(IDUM)*4.0)
        NE2=INT(PYR(IDUM)*4.0)
        N=2
        P(1,1)=0.0
        P(1,2)=0.0
        P(1,3)=-SQRT(E1**2-SM1**2)
        P(1,4)=E1
        P(1,5)=SM1
        K(1,1)=2
        K(1,2)=1
        K(1,3)=0
        K(1,4)=NE1
        K(1,5)=0
        P(2,1)=0.0
        P(2,2)=0.0
        P(2,3)=SQRT(E2**2-SM2**2)
        P(2,4)=E2
        P(2,5)=SM2
        K(2,1)=1
        K(2,2)=-1
        K(2,3)=0
        K(2,4)=NE2
        K(2,5)=0

        CALL AREXEC
        CALL LUEXEC

        IF(IPRINT.GT.0.AND.MOD(I,100).EQ.0) CALL LULIST(2)

 110  CONTINUE

      NERRA=0
      DO 200 I=1,40
        NERRA=NERRA+IWRN(I)
 200  CONTINUE

      NWRNA=IWRN(13)+IWRN(10)
      NERRA=NERRA-NWRNA
      IF(NERRA.EQ.0) THEN
        WRITE(MSTA(7),1000)
      ELSE
        WRITE(MSTA(7),1010) NERRA
      ENDIF

      IF(NWRNA.GT.0) WRITE(MSTA(7),1020) NWRNA

      NWRNJ=MSTU(27)
      NERRJ=MSTU(23)

      IF(NWRNJ+NERRJ.NE.0) WRITE(MSTA(7),1030) NWRNJ,NERRJ

 1000 FORMAT('No errors experienced by Ariadne.')
 1010 FORMAT(I5,' errors occurred in Ariadne.')
 1020 FORMAT(I5,' Non-serious warnings issued by Ariadne')
 1030 FORMAT(I5,' warnings and',I5,' errors occured in JETSET when ',
     $     'attempting to fragment',/
     $     ,' parton state produced by Ariadne.')

      RETURN

C**** END OF ARTEST ****************************************************
      END
C***********************************************************************
C $Id: arveto.F,v 0.8 1992/01/27 16:03:19 lonnblad Exp $

      REAL FUNCTION ARVET1()

C...ARiadne function VETo factor version 1

C...Determine the acceptance factor for chosen x_t^2 and y
C...Suitable for photon emission with constant alpha_EM


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)

      COMMON /ARDAT1/ PARA(40),MSTA(40)
      SAVE /ARDAT1/

      COMMON /ARINT1/ BC1,BC3,BZM,BZP,BP1,BM1,BP3,BM3,
     $                B1,B2,B3,XT2,XT,Y,QQ1,QQ3,NE1,NE3,
     $                S,W,C,CN,ALPHA0,XLAM2,IFLG,
     $                XT2MP,XT2ME,XT2M,XT2C,XTS,XT3,XT1,
     $                YINT,YMAX,YMIN,
     $                Y1,Y2,Y3,SY1,SY2,SY3,SSY,
     $                AE1,AE3,NXP1,NXP3,FQ1,FQ3
      SAVE /ARINT1/


      ARVET1=0.0
      IF(B2.LE.0) RETURN
      ARVET1=-((FQ1*(1.0-B1)/B2-FQ3*(1.0-B3)/B2)**2)*
     $         (B1**NXP1+B3**NXP3)*(YMAX-YMIN)*0.5/LOG(XT2)

      IF(MSTA(19).EQ.0) RETURN

      ARVET1=ARVET1*ARVETH()

      RETURN

C**** END OF ARVET1 ****************************************************
      END
C***********************************************************************

      REAL FUNCTION ARVET2()

C...ARiadne function VETo factor version 2

C...Determine the acceptance factor for chosen x_t^2 and y
C...Suitable for photon emission with running alpha_EM


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)

      COMMON /ARINT1/ BC1,BC3,BZM,BZP,BP1,BM1,BP3,BM3,
     $                B1,B2,B3,XT2,XT,Y,QQ1,QQ3,NE1,NE3,
     $                S,W,C,CN,ALPHA0,XLAM2,IFLG,
     $                XT2MP,XT2ME,XT2M,XT2C,XTS,XT3,XT1,
     $                YINT,YMAX,YMIN,
     $                Y1,Y2,Y3,SY1,SY2,SY3,SSY,
     $                AE1,AE3,NXP1,NXP3,FQ1,FQ3
      SAVE /ARINT1/


      ARVET2=ARVET1()*ULALEM(XT2*S)/ULALEM(0.25*S)

      RETURN

C**** END OF ARVET2 ****************************************************
      END
C***********************************************************************

      REAL FUNCTION ARVET3()

C...ARiadne function VETo factor version 3

C...Determine the acceptance factor for chosen x_t^2 and y
C...Suitable for gluon emission with constant alpha_QCD


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)

      COMMON /ARDAT1/ PARA(40),MSTA(40)
      SAVE /ARDAT1/

      COMMON /ARINT1/ BC1,BC3,BZM,BZP,BP1,BM1,BP3,BM3,
     $                B1,B2,B3,XT2,XT,Y,QQ1,QQ3,NE1,NE3,
     $                S,W,C,CN,ALPHA0,XLAM2,IFLG,
     $                XT2MP,XT2ME,XT2M,XT2C,XTS,XT3,XT1,
     $                YINT,YMAX,YMIN,
     $                Y1,Y2,Y3,SY1,SY2,SY3,SSY,
     $                AE1,AE3,NXP1,NXP3,FQ1,FQ3
      SAVE /ARINT1/


      ARVET3=-(B1**NXP1+B3**NXP3)*(YMAX-YMIN)*0.5/LOG(XT2)

      IF(MSTA(19).EQ.0) RETURN

      ARVET3=ARVET3*ARVETH()


      RETURN

C**** END OF ARVET3 ****************************************************
      END
C***********************************************************************

      REAL FUNCTION ARVET4()

C...ARiadne function VETo factor version 4

C...Determine the acceptance factor for chosen x_t^2 and y
C...Suitable for gluon emission with running alpha_QCD


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)

      COMMON /ARDAT1/ PARA(40),MSTA(40)
      SAVE /ARDAT1/

      COMMON /ARINT1/ BC1,BC3,BZM,BZP,BP1,BM1,BP3,BM3,
     $                B1,B2,B3,XT2,XT,Y,QQ1,QQ3,NE1,NE3,
     $                S,W,C,CN,ALPHA0,XLAM2,IFLG,
     $                XT2MP,XT2ME,XT2M,XT2C,XTS,XT3,XT1,
     $                YINT,YMAX,YMIN,
     $                Y1,Y2,Y3,SY1,SY2,SY3,SSY,
     $                AE1,AE3,NXP1,NXP3,FQ1,FQ3
      SAVE /ARINT1/


      ARVET4=(B1**NXP1+B3**NXP3)*(YMAX-YMIN)*0.5/YINT

      IF(MSTA(19).EQ.0) RETURN

      ARVET4=ARVET4*ARVETH()


      RETURN

C**** END OF ARVET4 ****************************************************
      END
C***********************************************************************

      REAL FUNCTION ARVET5()

C...ARiadne function VETo factor version 5

C...Determine the acceptance factor for chosen x_t^2 and y
C...Suitable for q-qbar emission


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)

      COMMON /ARINT1/ BC1,BC3,BZM,BZP,BP1,BM1,BP3,BM3,
     $                B1,B2,B3,XT2,XT,Y,QQ1,QQ3,NE1,NE3,
     $                S,W,C,CN,ALPHA0,XLAM2,IFLG,
     $                XT2MP,XT2ME,XT2M,XT2C,XTS,XT3,XT1,
     $                YINT,YMAX,YMIN,
     $                Y1,Y2,Y3,SY1,SY2,SY3,SSY,
     $                AE1,AE3,NXP1,NXP3,FQ1,FQ3
      SAVE /ARINT1/


      ARVET5=((1.0D0-B3+Y3)**2+(1.0D0-B2+Y2)**2)*XT*
     $         (EXP(-YMIN)-EXP(-YMAX))/YINT

      RETURN

C**** END OF ARVET5 ****************************************************
      END
C***********************************************************************

      REAL FUNCTION ARVETH()

C...ARiadne function Heavy VETo factor

C...Extra acceptance factor for heavy dipoles


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)

      COMMON /ARINT1/ BC1,BC3,BZM,BZP,BP1,BM1,BP3,BM3,
     $                B1,B2,B3,XT2,XT,Y,QQ1,QQ3,NE1,NE3,
     $                S,W,C,CN,ALPHA0,XLAM2,IFLG,
     $                XT2MP,XT2ME,XT2M,XT2C,XTS,XT3,XT1,
     $                YINT,YMAX,YMIN,
     $                Y1,Y2,Y3,SY1,SY2,SY3,SSY,
     $                AE1,AE3,NXP1,NXP3,FQ1,FQ3
      SAVE /ARINT1/


      ARVETH=0.0
      BX1=1.0-B1+Y1-Y3
      BX3=1.0-B3+Y3-Y1
      IF(B2.GE.1.0.OR.BX1.LE.0.OR.BX3.LE.0) RETURN
      BXM=BX1/BX3
      ARVETH=1.0-(Y1*BXM+Y3/BXM)/(1.0-B2)

      RETURN

C**** END OF ARVETH ****************************************************
      END


C***********************************************************************
C $Id: artune.F,v 0.4 1992/02/14 11:17:37 lonnblad Exp $

      SUBROUTINE ARTUNE(SET)

C...ARiadne subroutine TUNE

C...Sets parameters and switches in Ariadne and other programs which
C...Ariadne runs with.


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)

      COMMON /ARDAT1/ PARA(40),MSTA(40)
      SAVE /ARDAT1/

      COMMON /LUDAT1/ MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /LUDAT1/

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,XQ2,U
      SAVE /LEPTOU/
      CHARACTER SET*(*)


      IF(SET.EQ.'DELPHI') THEN
        PARA(1)=0.22
        PARA(3)=0.6
        PARA(5)=0.6
        MSTJ(11)=1
        PARJ(41)=0.23
        PARJ(42)=0.34
        PARJ(21)=0.405
        WRITE(MSTA(7),1000) SET
      ELSEIF(SET.EQ.'OPAL') THEN
        PARA(1)=0.20
        PARA(3)=1.0
        PARA(5)=1.0
        PARJ(41)=0.18
        PARJ(42)=0.34
        PARJ(21)=0.37
        WRITE(MSTA(7),1000) SET
      ELSE
        WRITE(MSTA(7),1010) SET
      ENDIF

 1000 FORMAT('Parameters and switches initialized using the "',A,
     $     '" tuning set')
 1010 FORMAT('Tuning set "',A,'" does not exist. Parameters and',
     $     ' switches retains their default value')

      RETURN

C**** END OF ARTUNE ****************************************************
      END
C***********************************************************************
C $Id: arprda.F,v 0.2 1992/01/27 16:03:19 lonnblad Exp $

      SUBROUTINE ARPRDA

C...ARiadne subroutine PRint DAta

C...Prints out parameters and switches used in Ariadne.


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)

      COMMON /ARDAT1/ PARA(40),MSTA(40)
      SAVE /ARDAT1/


      WRITE(MSTA(7),*)
      WRITE(MSTA(7),1000)
      DO 100 I=1,20
        WRITE(MSTA(7),1010) I,MSTA(I),MSTA(I+20),PARA(I),PARA(I+20)
 100  CONTINUE
      WRITE(MSTA(7),*)

 1000 FORMAT(10X,'Parameters and switches used by Ariadne:',/,/,
     $     '         I   MSTA(I) MSTA(I+20)   PARA(I) PARA(I+20)',/)
 1010 FORMAT(2I10,I11,3X,2G11.3)

      RETURN

C**** END OF ARPRDA ****************************************************
      END
C***********************************************************************
C $Id: archki.F,v 0.6 1992/03/12 12:09:11 lonnblad Exp $

      SUBROUTINE ARCHKI(ID,IOK)

C...ARiadne subroutine CHeck KInematics

C...Checks if the generated emission for dipole ID (or current dipole
C...if ID=0) is kinematically allowed.


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)

      COMMON /ARDAT1/ PARA(40),MSTA(40)
      SAVE /ARDAT1/

      COMMON /ARINT1/ BC1,BC3,BZM,BZP,BP1,BM1,BP3,BM3,
     $                B1,B2,B3,XT2,XT,Y,QQ1,QQ3,NE1,NE3,
     $                S,W,C,CN,ALPHA0,XLAM2,IFLG,
     $                XT2MP,XT2ME,XT2M,XT2C,XTS,XT3,XT1,
     $                YINT,YMAX,YMIN,
     $                Y1,Y2,Y3,SY1,SY2,SY3,SSY,
     $                AE1,AE3,NXP1,NXP3,FQ1,FQ3
      SAVE /ARINT1/

      COMMON /ARPART/ BP(MAXPAR,5),IFL(MAXPAR),IEX(MAXPAR),QQ(MAXPAR),
     $                IDI(MAXPAR),IDO(MAXPAR),INO(MAXPAR),IPART
      SAVE /ARPART/

      COMMON /ARDIPS/ BX1(MAXDIP),BX3(MAXDIP),PT2IN(MAXDIP),
     $                SDIP(MAXDIP),IP1(MAXDIP),IP3(MAXDIP),
     $                AEX1(MAXDIP),AEX3(MAXDIP),QDONE(MAXDIP),
     $                QEM(MAXDIP),IRAD(MAXDIP),ISTR(MAXDIP),IDIPS
      SAVE /ARDIPS/

      COMMON /ARDAT2/ PQMAS(10)
      SAVE /ARDAT2/

      IOK=0
      IF(ID.NE.0) THEN
        IFLG=IRAD(ID)
        NE1=IEX(IP1(ID))
        NE3=IEX(IP3(ID))
        QQ1=QQ(IP1(ID))
        QQ3=QQ(IP3(ID))
        S=SDIP(ID)
        W=SQRT(S)
        SY1=BP(IP1(ID),5)/W
        SY3=BP(IP3(ID),5)/W
        SY2=0.0
        IF(IFLG.NE.0) SY2=PQMAS(ABS(IFLG))/W
        SSY=SY1+SY2+SY3
        Y1=SY1**2
        Y2=SY2**2
        Y3=SY3**2
        BZP=0.5*(1.0+Y1-Y3+SQRT(1.0+(Y1-Y3)**2-2.0*(Y1+Y3)))
        BZM=0.5*(1.0+Y3-Y1+SQRT(1.0+(Y1-Y3)**2-2.0*(Y1+Y3)))
        B1=BX1(ID)
        B3=BX3(ID)
        AE1=AEX1(ID)
        AE3=AEX3(ID)
        BP1=(1.0-AE1)*BZP
        IF(BP1.LE.SY1) THEN
          BP1=0.0
          BM1=0.0
        ELSE
          BM1=Y1/BP1
        ENDIF
        BM3=(1.0-AE3)*BZM
        IF(BM3.LE.SY3) THEN
          BM3=0.0
          BP3=0.0
        ELSE
          BP3=Y3/BM3
        ENDIF
      ENDIF

      QRG=(MSTA(17).NE.0)
      QNG=(MSTA(17).EQ.1.AND.MSTA(16).GT.0.AND.
     $     ((QQ1.AND.NE1.EQ.0).OR.(QQ3.AND.NE3.EQ.0)))
      QRG=(QRG.AND.(.NOT.QNG))

      IF(IFLG.EQ.0.AND.(NE1.GT.0.OR.NE3.GT.0).AND.
     $     (BP1.GT.0.0.OR.BM3.GT.0.0).AND.QRG) THEN
        AZ1=1.0-BP1-BP3
        AZ3=1.0-BM1-BM3
        IF(AZ1.LE.0.0) RETURN
        IF(AZ3.LE.0.0) RETURN
        Y1A=Y1/(AZ1*AZ3)
        IF(BP1.GT.0.0) Y1A=0.0
        Y3A=Y3/(AZ1*AZ3)
        IF(BM3.GT.0.0) Y3A=0.0
        BE1=0.5*(1.0-(1.0-B1+Y1-Y3)/AZ1+Y1A-Y3A)
        BE3=0.5*(1.0-(1.0-B3+Y3-Y1)/AZ3+Y3A-Y1A)
        BE2=1.0-BE1-BE3
        BP1A=BE1**2-Y1A
        BP2A=BE2**2-Y2
        BP3A=BE3**2-Y3A
        IF(2.0*(BP1A*BP2A+BP2A*BP3A+BP3A*BP1A).LE.
     $       (BP1A**2+BP2A**2+BP3A**2)) RETURN
        IF(BE1.LT.SQRT(Y1A)) RETURN
        IF(BE2.LT.SY2) RETURN
        IF(BE3.LT.SQRT(Y3A)) RETURN
        PT21=B3
        PT23=B1
        IF(BP1.GT.0.0.AND.BM3.GT.0.0.AND.MSTA(17).GE.2) THEN
          BP2=1.0-B1+Y1-Y3
          BM2=1.0-B3+Y3-Y1
          PT21=(BM2*BP2**3)/(1.0-BP1-BP2)**2
          PT23=(BP2*BM2**3)/(1.0-BM3-BM2)**2
        ENDIF
        IF((BP1.GT.0.0.AND.BM3.GT.0.0.AND.PT21.GE.PT23).OR.
     $       (BM3.GT.0.0.AND.BP1.LE.0.0)) THEN
          BM3=0.0
          BP3=0.0
        ELSE
          BP1=0.0
          BM1=0.0
        ENDIF
        AZ1=1.0-BP1-BP3
        AZ3=1.0-BM1-BM3
        IF(AZ1.LE.0.0) RETURN
        IF(AZ3.LE.0.0) RETURN
        Y1A=Y1/(AZ1*AZ3)
        IF(BP1.GT.0.0) Y1A=0.0
        Y3A=Y3/(AZ1*AZ3)
        IF(BM3.GT.0.0) Y3A=0.0
        BE1=0.5*(1.0-(1.0-B1+Y1-Y3)/AZ1+Y1A-Y3A)
        BE3=0.5*(1.0-(1.0-B3+Y3-Y1)/AZ3+Y3A-Y1A)
        BE2=1.0-BE1-BE3
        BP1A=BE1**2-Y1A
        BP2A=BE2**2-Y2
        BP3A=BE3**2-Y3A
        IF(2.0*(BP1A*BP2A+BP2A*BP3A+BP3A*BP1A).LE.
     $       (BP1A**2+BP2A**2+BP3A**2)) RETURN
        IF(BE1.LT.SQRT(Y1A)) RETURN
        IF(BE2.LT.SY2) RETURN
        IF(BE3.LT.SQRT(Y3A)) RETURN
      ELSE
        BE1=0.5*B1
        BE3=0.5*B3
        BE2=1.0-BE1-BE3
        BP1A=BE1**2-Y1
        BP2A=BE2**2-Y2
        BP3A=BE3**2-Y3
        IF(2.0*(BP1A*BP2A+BP2A*BP3A+BP3A*BP1A).LE.
     $       (BP1A**2+BP2A**2+BP3A**2)) RETURN
        IF(BE1.LT.SY1) RETURN
        IF(BE2.LT.SY2) RETURN
        IF(BE3.LT.SY3) RETURN
      ENDIF

      IOK=1

      RETURN

C**** END OF ARCHKI ****************************************************
      END
C***********************************************************************
C $Id: arexma.F,v 0.1 1992/02/11 10:36:11 lonnblad Exp $

      SUBROUTINE AREXMA(I1,I3)

C...ARiadne subroutine make EXtended partons MAssless

C...Makes extended partons massless.


      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)

      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)

      COMMON /ARPART/ BP(MAXPAR,5),IFL(MAXPAR),IEX(MAXPAR),QQ(MAXPAR),
     $                IDI(MAXPAR),IDO(MAXPAR),INO(MAXPAR),IPART
      SAVE /ARPART/

      COMMON /ARDAT1/ PARA(40),MSTA(40)
      SAVE /ARDAT1/


      IF(MSTA(31).GT.0) RETURN
      IF(IEX(I1).EQ.0.AND.IEX(I3).EQ.0) RETURN
      DPE1=BP(I1,4)
      DPE3=BP(I3,4)
      DPE=DPE1+DPE3
      DPX1=BP(I1,1)
      DPX3=BP(I3,1)
      DBEX=(DPX1+DPX3)/DPE
      DPY1=BP(I1,2)
      DPY3=BP(I3,2)
      DBEY=(DPY1+DPY3)/DPE
      DPZ1=BP(I1,3)
      DPZ3=BP(I3,3)
      DBEZ=(DPZ1+DPZ3)/DPE
      CALL AROBO2(0.0,0.0,-DBEX,-DBEY,-DBEZ,I1,I3)

C...Calculate rotation angles but no need for rotation yet
      PX=BP(I1,1)
      PY=BP(I1,2)
      PZ=BP(I1,3)
      PHI=ULANGL(PX,PY)
      THE=ULANGL(PZ,SQRT(PX**2+PY**2))
      CALL AROBO2(0.0,-PHI,0.0D0,0.0D0,0.0D0,I1,I3)
      CALL AROBO2(-THE,0.0,0.0D0,0.0D0,0.0D0,I1,I3)
      IF(IEX(I1).GT.0) BP(I1,5)=0.0
      IF(IEX(I3).GT.0) BP(I3,5)=0.0
      BE=BP(I1,4)+BP(I3,4)
      BP(I1,4)=0.5*(BE**2+BP(I1,5)**2-BP(I3,5)**2)/BE
      BP(I3,4)=BE-BP(I1,4)
      BP(I1,3)=SQRT(BP(I1,4)**2-BP(I1,5)**2)
      BP(I3,3)=-BP(I1,3)
      BP(I1,2)=0.0
      BP(I3,2)=0.0
      BP(I1,1)=0.0
      BP(I3,1)=0.0

      CALL AROBO2(THE,PHI,DBEX,DBEY,DBEZ,I1,I3)

      RETURN

C**** END OF AREXMA ****************************************************
      END
