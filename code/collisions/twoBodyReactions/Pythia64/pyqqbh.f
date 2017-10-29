 
C***********************************************************************
 
C...PYQQBH
C...Calculates the matrix element for the processes
C...g + g or q + qbar -> Q + Qbar + H (normally with Q = t).
C...REDUCE output and part of the rest courtesy Z. Kunszt, see
C...Z. Kunszt, Nucl. Phys. B247 (1984) 339.
 
      SUBROUTINE PYQQBH(WTQQBH)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      COMMON/PYINT2/ISET(500),KFPR(500,2),COEF(500,20),ICOL(40,4,2)
      SAVE /PYDAT1/,/PYDAT2/,/PYPARS/,/PYINT1/,/PYINT2/
C...Local arrays and function.
      DIMENSION PP(15,4),CLR(8,8),FM(10,10),RM(8,8),DX(8)
      DOT(I,J)=PP(I,4)*PP(J,4)-PP(I,1)*PP(J,1)-PP(I,2)*PP(J,2)-
     &PP(I,3)*PP(J,3)
 
C...Mass parameters.
      WTQQBH=0D0
      ISUB=MINT(1)
      SHPR=SQRT(VINT(26))*VINT(1)
      PQ=PMAS(PYCOMP(KFPR(ISUB,2)),1)
      PH=SQRT(VINT(21))*VINT(1)
      SPQ=PQ**2
      SPH=PH**2
 
C...Set up outgoing kinematics: 1=t, 2=tbar, 3=H.
      DO 100 I=1,2
        PT=SQRT(MAX(0D0,VINT(197+5*I)))
        PP(I,1)=PT*COS(VINT(198+5*I))
        PP(I,2)=PT*SIN(VINT(198+5*I))
  100 CONTINUE
      PP(3,1)=-PP(1,1)-PP(2,1)
      PP(3,2)=-PP(1,2)-PP(2,2)
      PMS1=SPQ+PP(1,1)**2+PP(1,2)**2
      PMS2=SPQ+PP(2,1)**2+PP(2,2)**2
      PMS3=SPH+PP(3,1)**2+PP(3,2)**2
      PMT3=SQRT(PMS3)
      PP(3,3)=PMT3*SINH(VINT(211))
      PP(3,4)=PMT3*COSH(VINT(211))
      PMS12=(SHPR-PP(3,4))**2-PP(3,3)**2
      PP(1,3)=(-PP(3,3)*(PMS12+PMS1-PMS2)+
     &VINT(213)*(SHPR-PP(3,4))*VINT(220))/(2D0*PMS12)
      PP(2,3)=-PP(1,3)-PP(3,3)
      PP(1,4)=SQRT(PMS1+PP(1,3)**2)
      PP(2,4)=SQRT(PMS2+PP(2,3)**2)
 
C...Set up incoming kinematics and derived momentum combinations.
      DO 110 I=4,5
        PP(I,1)=0D0
        PP(I,2)=0D0
        PP(I,3)=-0.5D0*SHPR*(-1)**I
        PP(I,4)=-0.5D0*SHPR
  110 CONTINUE
      DO 120 J=1,4
        PP(6,J)=PP(1,J)+PP(2,J)
        PP(7,J)=PP(1,J)+PP(3,J)
        PP(8,J)=PP(1,J)+PP(4,J)
        PP(9,J)=PP(1,J)+PP(5,J)
        PP(10,J)=-PP(2,J)-PP(3,J)
        PP(11,J)=-PP(2,J)-PP(4,J)
        PP(12,J)=-PP(2,J)-PP(5,J)
        PP(13,J)=-PP(4,J)-PP(5,J)
  120 CONTINUE
 
C...Derived kinematics invariants.
      X1=DOT(1,2)
      X2=DOT(1,3)
      X3=DOT(1,4)
      X4=DOT(1,5)
      X5=DOT(2,3)
      X6=DOT(2,4)
      X7=DOT(2,5)
      X8=DOT(3,4)
      X9=DOT(3,5)
      X10=DOT(4,5)
 
C...Propagators.
      SS1=DOT(7,7)-SPQ
      SS2=DOT(8,8)-SPQ
      SS3=DOT(9,9)-SPQ
      SS4=DOT(10,10)-SPQ
      SS5=DOT(11,11)-SPQ
      SS6=DOT(12,12)-SPQ
      SS7=DOT(13,13)
      DX(1)=SS1*SS6
      DX(2)=SS2*SS6
      DX(3)=SS2*SS4
      DX(4)=SS1*SS5
      DX(5)=SS3*SS5
      DX(6)=SS3*SS4
      DX(7)=SS7*SS1
      DX(8)=SS7*SS4
 
C...Define colour coefficients for g + g -> Q + Qbar + H.
      IF(ISUB.EQ.121.OR.ISUB.EQ.181.OR.ISUB.EQ.186) THEN
        DO 140 I=1,3
          DO 130 J=1,3
            CLR(I,J)=16D0/3D0
            CLR(I+3,J+3)=16D0/3D0
            CLR(I,J+3)=-2D0/3D0
            CLR(I+3,J)=-2D0/3D0
  130     CONTINUE
  140   CONTINUE
        DO 160 L=1,2
          DO 150 I=1,3
            CLR(I,6+L)=-6D0
            CLR(I+3,6+L)=6D0
            CLR(6+L,I)=-6D0
            CLR(6+L,I+3)=6D0
  150     CONTINUE
  160   CONTINUE
        DO 180 K1=1,2
          DO 170 K2=1,2
            CLR(6+K1,6+K2)=12D0
  170     CONTINUE
  180   CONTINUE
 
C...Evaluate matrix elements for g + g -> Q + Qbar + H.
        FM(1,1)=64*PQ**6+16*PQ**4*PH**2+32*PQ**4*(X1+2*X2+X4+X9+2*
     &  X7+X5)+8*PQ**2*PH**2*(-X1-X4+2*X7)+16*PQ**2*(X2*X9+4*X2*
     &  X7+X2*X5-2*X4*X7-2*X9*X7)+8*PH**2*X4*X7-16*X2*X9*X7
        FM(1,2)=16*PQ**6+8*PQ**4*(-2*X1+X2-2*X3-2*X4-4*X10+X9-X8+2
     &  *X7-4*X6+X5)+8*PQ**2*(-2*X1*X2-2*X2*X4-2*X2*X10+X2*X7-2*
     &  X2*X6-2*X3*X7+2*X4*X7+4*X10*X7-X9*X7-X8*X7)+16*X2*X7*(X4+
     &  X10)
        FM(1,3)=16*PQ**6-4*PQ**4*PH**2+8*PQ**4*(-2*X1+2*X2-2*X3-4*
     &  X4-8*X10+X9+X8-2*X7-4*X6+2*X5)-(4*PQ**2*PH**2)*(X1+X4+X10
     &  +X6)+8*PQ**2*(-2*X1*X2-2*X1*X10+X1*X9+X1*X8-2*X1*X5+X2**2
     &  -4*X2*X4-5*X2*X10+X2*X8-X2*X7-3*X2*X6+X2*X5+X3*X9+2*X3*X7
     &  -X3*X5+X4*X8+2*X4*X6-3*X4*X5-5*X10*X5+X9*X8+X9*X6+X9*X5+
     &  X8*X7-4*X6*X5+X5**2)-(16*X2*X5)*(X1+X4+X10+X6)
        FM(1,4)=16*PQ**6+4*PQ**4*PH**2+16*PQ**4*(-X1+X2-X3-X4+X10-
     &  X9-X8+2*X7+2*X6-X5)+4*PQ**2*PH**2*(X1+X3+X4+X10+2*X7+2*X6
     &  )+8*PQ**2*(4*X1*X10+4*X1*X7+4*X1*X6+2*X2*X10-X2*X9-X2*X8+
     &  4*X2*X7+4*X2*X6-X2*X5+4*X10*X5+4*X7*X5+4*X6*X5)-(8*PH**2*
     &  X1)*(X10+X7+X6)+16*X2*X5*(X10+X7+X6)
        FM(1,5)=8*PQ**4*(-2*X1-2*X4+X10-X9)+4*PQ**2*(4*X1**2-2*X1*
     &  X2+8*X1*X3+6*X1*X10-2*X1*X9+4*X1*X8+4*X1*X7+4*X1*X6+2*X1*
     &  X5+X2*X10+4*X3*X4-X3*X9+2*X3*X7+3*X4*X8-2*X4*X6+2*X4*X5-4
     &  *X10*X7+3*X10*X5-3*X9*X6+3*X8*X7-4*X7**2+4*X7*X5)+8*(X1**
     &  2*X9-X1**2*X8-X1*X2*X7+X1*X2*X6+X1*X3*X9+X1*X3*X5-X1*X4*
     &  X8-X1*X4*X5+X1*X10*X9+X1*X9*X7+X1*X9*X6-X1*X8*X7-X2*X3*X7
     &  +X2*X4*X6-X2*X10*X7-X2*X7**2+X3*X7*X5-X4*X10*X5-X4*X7*X5-
     &  X4*X6*X5)
        FM(1,6)=16*PQ**4*(-4*X1-X4+X9-X7)+4*PQ**2*PH**2*(-2*X1-X4-
     &  X7)+16*PQ**2*(-2*X1**2-3*X1*X2-2*X1*X4-3*X1*X9-2*X1*X7-3*
     &  X1*X5-2*X2*X4-2*X7*X5)-8*PH**2*X4*X7+8*(-X1*X2*X9-2*X1*X2
     &  *X5-X1*X9**2-X1*X9*X5+X2**2*X7-X2*X4*X5+X2*X9*X7-X2*X7*X5
     &  +X4*X9*X5+X4*X5**2)
        FM(1,7)=8*PQ**4*(2*X3+X4+3*X10+X9+2*X8+3*X7+6*X6)+2*PQ**2*
     &  PH**2*(-2*X3-X4+3*X10+3*X7+6*X6)+4*PQ**2*(4*X1*X10+4*X1*
     &  X7+8*X1*X6+6*X2*X10+X2*X9+2*X2*X8+6*X2*X7+12*X2*X6-8*X3*
     &  X7+4*X4*X7+4*X4*X6+4*X10*X5+4*X9*X7+4*X9*X6-8*X8*X7+4*X7*
     &  X5+8*X6*X5)+4*PH**2*(-X1*X10-X1*X7-2*X1*X6+2*X3*X7-X4*X7-
     &  X4*X6)+8*X2*(X10*X5+X9*X7+X9*X6-2*X8*X7+X7*X5+2*X6*X5)
        FM(1,8)=8*PQ**4*(2*X3+X4+3*X10+2*X9+X8+3*X7+6*X6)+2*PQ**2*
     &  PH**2*(-2*X3-X4+2*X10+X7+2*X6)+4*PQ**2*(4*X1*X10-2*X1*X9+
     &  2*X1*X8+4*X1*X7+8*X1*X6+5*X2*X10+2*X2*X9+X2*X8+4*X2*X7+8*
     &  X2*X6-X3*X9-8*X3*X7+2*X3*X5+2*X4*X9-X4*X8+4*X4*X7+4*X4*X6
     &  +4*X4*X5+5*X10*X5+X9**2-X9*X8+2*X9*X7+5*X9*X6+X9*X5-7*X8*
     &  X7+2*X8*X5+2*X7*X5+10*X6*X5)+2*PH**2*(-X1*X10+X3*X7-2*X4*
     &  X7+X4*X6)+4*(-X1*X9**2+X1*X9*X8-2*X1*X9*X5-X1*X8*X5+2*X2*
     &  X10*X5+X2*X9*X7+X2*X9*X6-2*X2*X8*X7+3*X2*X6*X5+X3*X9*X5+
     &  X3*X5**2+X4*X9*X5-2*X4*X8*X5+2*X4*X5**2)
        FM(2,2)=16*PQ**6+16*PQ**4*(-X1+X3-X4-X10+X7-X6)+16*PQ**2*(
     &  X3*X10+X3*X7+X3*X6+X4*X7+X10*X7)-16*X3*X10*X7
        FM(2,3)=16*PQ**6+8*PQ**4*(-2*X1+X2+2*X3-4*X4-4*X10-X9+X8-2
     &  *X7-2*X6+X5)+8*PQ**2*(-2*X1*X5+4*X3*X10-X3*X9-X3*X8-2*X3*
     &  X7+2*X3*X6+X3*X5-2*X4*X5-2*X10*X5-2*X6*X5)+16*X3*X5*(X10+
     &  X6)
        FM(2,4)=8*PQ**4*(-2*X1-2*X3+X10-X8)+4*PQ**2*(4*X1**2-2*X1*
     &  X2+8*X1*X4+6*X1*X10+4*X1*X9-2*X1*X8+4*X1*X7+4*X1*X6+2*X1*
     &  X5+X2*X10+4*X3*X4+3*X3*X9-2*X3*X7+2*X3*X5-X4*X8+2*X4*X6-4
     &  *X10*X6+3*X10*X5+3*X9*X6-3*X8*X7-4*X6**2+4*X6*X5)+8*(-X1
     &  **2*X9+X1**2*X8+X1*X2*X7-X1*X2*X6-X1*X3*X9-X1*X3*X5+X1*X4
     &  *X8+X1*X4*X5+X1*X10*X8-X1*X9*X6+X1*X8*X7+X1*X8*X6+X2*X3*
     &  X7-X2*X4*X6-X2*X10*X6-X2*X6**2-X3*X10*X5-X3*X7*X5-X3*X6*
     &  X5+X4*X6*X5)
        FM(2,5)=16*PQ**4*X10+8*PQ**2*(2*X1**2+2*X1*X3+2*X1*X4+2*X1
     &  *X10+2*X1*X7+2*X1*X6+X3*X7+X4*X6)+8*(-2*X1**3-2*X1**2*X3-
     &  2*X1**2*X4-2*X1**2*X10-2*X1**2*X7-2*X1**2*X6-2*X1*X3*X4-
     &  X1*X3*X10-2*X1*X3*X6-X1*X4*X10-2*X1*X4*X7-X1*X10**2-X1*
     &  X10*X7-X1*X10*X6-2*X1*X7*X6+X3**2*X7-X3*X4*X7-X3*X4*X6+X3
     &  *X10*X7+X3*X7**2-X3*X7*X6+X4**2*X6+X4*X10*X6-X4*X7*X6+X4*
     &  X6**2)
        FM(2,6)=8*PQ**4*(-2*X1+X10-X9-2*X7)+4*PQ**2*(4*X1**2+2*X1*
     &  X2+4*X1*X3+4*X1*X4+6*X1*X10-2*X1*X9+4*X1*X8+8*X1*X6-2*X1*
     &  X5+4*X2*X4+3*X2*X10+2*X2*X7-3*X3*X9-2*X3*X7-4*X4**2-4*X4*
     &  X10+3*X4*X8+2*X4*X6+X10*X5-X9*X6+3*X8*X7+4*X7*X6)+8*(X1**
     &  2*X9-X1**2*X8-X1*X2*X7+X1*X2*X6+X1*X3*X9+X1*X3*X5+X1*X4*
     &  X9-X1*X4*X8-X1*X4*X5+X1*X10*X9+X1*X9*X6-X1*X8*X7-X2*X3*X7
     &  -X2*X4*X7+X2*X4*X6-X2*X10*X7+X3*X7*X5-X4**2*X5-X4*X10*X5-
     &  X4*X6*X5)
        FM(2,7)=8*PQ**4*(X3+2*X4+3*X10+X7+2*X6)+4*PQ**2*(-4*X1*X3-
     &  2*X1*X4-2*X1*X10+X1*X9-X1*X8-4*X1*X7-2*X1*X6+X2*X3+2*X2*
     &  X4+3*X2*X10+X2*X7+2*X2*X6-6*X3*X4-6*X3*X10-2*X3*X9-2*X3*
     &  X7-4*X3*X6-X3*X5-6*X4**2-6*X4*X10-3*X4*X9-X4*X8-4*X4*X7-2
     &  *X4*X6-2*X4*X5-3*X10*X9-3*X10*X8-6*X10*X7-6*X10*X6+X10*X5
     &  +X9*X7-2*X8*X7-2*X8*X6-6*X7*X6+X7*X5-6*X6**2+2*X6*X5)+4*(
     &  -X1**2*X9+X1**2*X8-2*X1*X2*X10-3*X1*X2*X7-3*X1*X2*X6+X1*
     &  X3*X9-X1*X3*X5+X1*X4*X9+X1*X4*X8+X1*X4*X5+X1*X10*X9+X1*
     &  X10*X8-X1*X9*X6+X1*X8*X6+X2*X3*X7-3*X2*X4*X7-X2*X4*X6-3*
     &  X2*X10*X7-3*X2*X10*X6-3*X2*X7*X6-3*X2*X6**2-2*X3*X4*X5-X3
     &  *X10*X5-X3*X6*X5-X4**2*X5-X4*X10*X5+X4*X6*X5)
        FM(2,8)=8*PQ**4*(X3+2*X4+3*X10+X7+2*X6)+4*PQ**2*(-4*X1*X3-
     &  2*X1*X4-2*X1*X10-X1*X9+X1*X8-4*X1*X7-2*X1*X6+X2*X3+2*X2*
     &  X4+X2*X10-X2*X7-2*X2*X6-6*X3*X4-6*X3*X10-2*X3*X9+X3*X8-2*
     &  X3*X7-4*X3*X6+X3*X5-6*X4**2-6*X4*X10-2*X4*X9-4*X4*X7-2*X4
     &  *X6+2*X4*X5-3*X10*X9-3*X10*X8-6*X10*X7-6*X10*X6+3*X10*X5-
     &  X9*X6-2*X8*X7-3*X8*X6-6*X7*X6+X7*X5-6*X6**2+2*X6*X5)+4*(
     &  X1**2*X9-X1**2*X8-X1*X2*X7+X1*X2*X6-3*X1*X3*X5+X1*X4*X9-
     &  X1*X4*X8-3*X1*X4*X5+X1*X10*X9+X1*X10*X8-2*X1*X10*X5+X1*X9
     &  *X6+X1*X8*X7+X1*X8*X6-X2*X4*X7+X2*X4*X6-X2*X10*X7-X2*X10*
     &  X6-2*X2*X7*X6-X2*X6**2-3*X3*X4*X5-3*X3*X10*X5+X3*X7*X5-3*
     &  X3*X6*X5-3*X4**2*X5-3*X4*X10*X5-X4*X6*X5)
        FM(3,3)=64*PQ**6+16*PQ**4*PH**2+32*PQ**4*(X1+X2+2*X3+X8+X6
     &  +2*X5)+8*PQ**2*PH**2*(-X1+2*X3-X6)+16*PQ**2*(X2*X5-2*X3*
     &  X8-2*X3*X6+4*X3*X5+X8*X5)+8*PH**2*X3*X6-16*X3*X8*X5
        FM(3,4)=16*PQ**4*(-4*X1-X3+X8-X6)+4*PQ**2*PH**2*(-2*X1-X3-
     &  X6)+16*PQ**2*(-2*X1**2-3*X1*X2-2*X1*X3-3*X1*X8-2*X1*X6-3*
     &  X1*X5-2*X2*X3-2*X6*X5)-8*PH**2*X3*X6+8*(-X1*X2*X8-2*X1*X2
     &  *X5-X1*X8**2-X1*X8*X5+X2**2*X6-X2*X3*X5+X2*X8*X6-X2*X6*X5
     &  +X3*X8*X5+X3*X5**2)
        FM(3,5)=8*PQ**4*(-2*X1+X10-X8-2*X6)+4*PQ**2*(4*X1**2+2*X1*
     &  X2+4*X1*X3+4*X1*X4+6*X1*X10+4*X1*X9-2*X1*X8+8*X1*X7-2*X1*
     &  X5+4*X2*X3+3*X2*X10+2*X2*X6-4*X3**2-4*X3*X10+3*X3*X9+2*X3
     &  *X7-3*X4*X8-2*X4*X6+X10*X5+3*X9*X6-X8*X7+4*X7*X6)+8*(-X1
     &  **2*X9+X1**2*X8+X1*X2*X7-X1*X2*X6-X1*X3*X9+X1*X3*X8-X1*X3
     &  *X5+X1*X4*X8+X1*X4*X5+X1*X10*X8-X1*X9*X6+X1*X8*X7+X2*X3*
     &  X7-X2*X3*X6-X2*X4*X6-X2*X10*X6-X3**2*X5-X3*X10*X5-X3*X7*
     &  X5+X4*X6*X5)
        FM(3,6)=16*PQ**6+4*PQ**4*PH**2+16*PQ**4*(-X1-X2+2*X3+2*X4+
     &  X10-X9-X8-X7-X6+X5)+4*PQ**2*PH**2*(X1+2*X3+2*X4+X10+X7+X6
     &  )+8*PQ**2*(4*X1*X3+4*X1*X4+4*X1*X10+4*X2*X3+4*X2*X4+4*X2*
     &  X10-X2*X5+4*X3*X5+4*X4*X5+2*X10*X5-X9*X5-X8*X5)-(8*PH**2*
     &  X1)*(X3+X4+X10)+16*X2*X5*(X3+X4+X10)
        FM(3,7)=8*PQ**4*(3*X3+6*X4+3*X10+X9+2*X8+2*X7+X6)+2*PQ**2*
     &  PH**2*(X3+2*X4+2*X10-2*X7-X6)+4*PQ**2*(4*X1*X3+8*X1*X4+4*
     &  X1*X10+2*X1*X9-2*X1*X8+2*X2*X3+10*X2*X4+5*X2*X10+2*X2*X9+
     &  X2*X8+2*X2*X7+4*X2*X6-7*X3*X9+2*X3*X8-8*X3*X7+4*X3*X6+4*
     &  X3*X5+5*X4*X8+4*X4*X6+8*X4*X5+5*X10*X5-X9*X8-X9*X6+X9*X5+
     &  X8**2-X8*X7+2*X8*X6+2*X8*X5)+2*PH**2*(-X1*X10+X3*X7-2*X3*
     &  X6+X4*X6)+4*(-X1*X2*X9-2*X1*X2*X8+X1*X9*X8-X1*X8**2+X2**2
     &  *X7+2*X2**2*X6+3*X2*X4*X5+2*X2*X10*X5-2*X2*X9*X6+X2*X8*X7
     &  +X2*X8*X6-2*X3*X9*X5+X3*X8*X5+X4*X8*X5)
        FM(3,8)=8*PQ**4*(3*X3+6*X4+3*X10+2*X9+X8+2*X7+X6)+2*PQ**2*
     &  PH**2*(3*X3+6*X4+3*X10-2*X7-X6)+4*PQ**2*(4*X1*X3+8*X1*X4+
     &  4*X1*X10+4*X2*X3+8*X2*X4+4*X2*X10-8*X3*X9+4*X3*X8-8*X3*X7
     &  +4*X3*X6+6*X3*X5+4*X4*X8+4*X4*X6+12*X4*X5+6*X10*X5+2*X9*
     &  X5+X8*X5)+4*PH**2*(-X1*X3-2*X1*X4-X1*X10+2*X3*X7-X3*X6-X4
     &  *X6)+8*X5*(X2*X3+2*X2*X4+X2*X10-2*X3*X9+X3*X8+X4*X8)
        FM(4,4)=64*PQ**6+16*PQ**4*PH**2+32*PQ**4*(X1+2*X2+X3+X8+2*
     &  X6+X5)+8*PQ**2*PH**2*(-X1-X3+2*X6)+16*PQ**2*(X2*X8+4*X2*
     &  X6+X2*X5-2*X3*X6-2*X8*X6)+8*PH**2*X3*X6-16*X2*X8*X6
        FM(4,5)=16*PQ**6+8*PQ**4*(-2*X1+X2-2*X3-2*X4-4*X10-X9+X8-4
     &  *X7+2*X6+X5)+8*PQ**2*(-2*X1*X2-2*X2*X3-2*X2*X10-2*X2*X7+
     &  X2*X6+2*X3*X6-2*X4*X6+4*X10*X6-X9*X6-X8*X6)+16*X2*X6*(X3+
     &  X10)
        FM(4,6)=16*PQ**6-4*PQ**4*PH**2+8*PQ**4*(-2*X1+2*X2-4*X3-2*
     &  X4-8*X10+X9+X8-4*X7-2*X6+2*X5)-(4*PQ**2*PH**2)*(X1+X3+X10
     &  +X7)+8*PQ**2*(-2*X1*X2-2*X1*X10+X1*X9+X1*X8-2*X1*X5+X2**2
     &  -4*X2*X3-5*X2*X10+X2*X9-3*X2*X7-X2*X6+X2*X5+X3*X9+2*X3*X7
     &  -3*X3*X5+X4*X8+2*X4*X6-X4*X5-5*X10*X5+X9*X8+X9*X6+X8*X7+
     &  X8*X5-4*X7*X5+X5**2)-(16*X2*X5)*(X1+X3+X10+X7)
        FM(4,7)=8*PQ**4*(-X3-2*X4-3*X10-2*X9-X8-6*X7-3*X6)+2*PQ**2
     &  *PH**2*(X3+2*X4-3*X10-6*X7-3*X6)+4*PQ**2*(-4*X1*X10-8*X1*
     &  X7-4*X1*X6-6*X2*X10-2*X2*X9-X2*X8-12*X2*X7-6*X2*X6-4*X3*
     &  X7-4*X3*X6+8*X4*X6-4*X10*X5+8*X9*X6-4*X8*X7-4*X8*X6-8*X7*
     &  X5-4*X6*X5)+4*PH**2*(X1*X10+2*X1*X7+X1*X6+X3*X7+X3*X6-2*
     &  X4*X6)+8*X2*(-X10*X5+2*X9*X6-X8*X7-X8*X6-2*X7*X5-X6*X5)
        FM(4,8)=8*PQ**4*(-X3-2*X4-3*X10-X9-2*X8-6*X7-3*X6)+2*PQ**2
     &  *PH**2*(X3+2*X4-2*X10-2*X7-X6)+4*PQ**2*(-4*X1*X10-2*X1*X9
     &  +2*X1*X8-8*X1*X7-4*X1*X6-5*X2*X10-X2*X9-2*X2*X8-8*X2*X7-4
     &  *X2*X6+X3*X9-2*X3*X8-4*X3*X7-4*X3*X6-4*X3*X5+X4*X8+8*X4*
     &  X6-2*X4*X5-5*X10*X5+X9*X8+7*X9*X6-2*X9*X5-X8**2-5*X8*X7-2
     &  *X8*X6-X8*X5-10*X7*X5-2*X6*X5)+2*PH**2*(X1*X10-X3*X7+2*X3
     &  *X6-X4*X6)+4*(-X1*X9*X8+X1*X9*X5+X1*X8**2+2*X1*X8*X5-2*X2
     &  *X10*X5+2*X2*X9*X6-X2*X8*X7-X2*X8*X6-3*X2*X7*X5+2*X3*X9*
     &  X5-X3*X8*X5-2*X3*X5**2-X4*X8*X5-X4*X5**2)
        FM(5,5)=16*PQ**6+16*PQ**4*(-X1-X3+X4-X10-X7+X6)+16*PQ**2*(
     &  X3*X6+X4*X10+X4*X7+X4*X6+X10*X6)-16*X4*X10*X6
        FM(5,6)=16*PQ**6+8*PQ**4*(-2*X1+X2-4*X3+2*X4-4*X10+X9-X8-2
     &  *X7-2*X6+X5)+8*PQ**2*(-2*X1*X5-2*X3*X5+4*X4*X10-X4*X9-X4*
     &  X8+2*X4*X7-2*X4*X6+X4*X5-2*X10*X5-2*X7*X5)+16*X4*X5*(X10+
     &  X7)
        FM(5,7)=8*PQ**4*(-2*X3-X4-3*X10-2*X7-X6)+4*PQ**2*(2*X1*X3+
     &  4*X1*X4+2*X1*X10+X1*X9-X1*X8+2*X1*X7+4*X1*X6-2*X2*X3-X2*
     &  X4-3*X2*X10-2*X2*X7-X2*X6+6*X3**2+6*X3*X4+6*X3*X10+X3*X9+
     &  3*X3*X8+2*X3*X7+4*X3*X6+2*X3*X5+6*X4*X10+2*X4*X8+4*X4*X7+
     &  2*X4*X6+X4*X5+3*X10*X9+3*X10*X8+6*X10*X7+6*X10*X6-X10*X5+
     &  2*X9*X7+2*X9*X6-X8*X6+6*X7**2+6*X7*X6-2*X7*X5-X6*X5)+4*(-
     &  X1**2*X9+X1**2*X8+2*X1*X2*X10+3*X1*X2*X7+3*X1*X2*X6-X1*X3
     &  *X9-X1*X3*X8-X1*X3*X5-X1*X4*X8+X1*X4*X5-X1*X10*X9-X1*X10*
     &  X8-X1*X9*X7+X1*X8*X7+X2*X3*X7+3*X2*X3*X6-X2*X4*X6+3*X2*
     &  X10*X7+3*X2*X10*X6+3*X2*X7**2+3*X2*X7*X6+X3**2*X5+2*X3*X4
     &  *X5+X3*X10*X5-X3*X7*X5+X4*X10*X5+X4*X7*X5)
        FM(5,8)=8*PQ**4*(-2*X3-X4-3*X10-2*X7-X6)+4*PQ**2*(2*X1*X3+
     &  4*X1*X4+2*X1*X10-X1*X9+X1*X8+2*X1*X7+4*X1*X6-2*X2*X3-X2*
     &  X4-X2*X10+2*X2*X7+X2*X6+6*X3**2+6*X3*X4+6*X3*X10+2*X3*X8+
     &  2*X3*X7+4*X3*X6-2*X3*X5+6*X4*X10-X4*X9+2*X4*X8+4*X4*X7+2*
     &  X4*X6-X4*X5+3*X10*X9+3*X10*X8+6*X10*X7+6*X10*X6-3*X10*X5+
     &  3*X9*X7+2*X9*X6+X8*X7+6*X7**2+6*X7*X6-2*X7*X5-X6*X5)+4*(
     &  X1**2*X9-X1**2*X8-X1*X2*X7+X1*X2*X6+X1*X3*X9-X1*X3*X8+3*
     &  X1*X3*X5+3*X1*X4*X5-X1*X10*X9-X1*X10*X8+2*X1*X10*X5-X1*X9
     &  *X7-X1*X9*X6-X1*X8*X7-X2*X3*X7+X2*X3*X6+X2*X10*X7+X2*X10*
     &  X6+X2*X7**2+2*X2*X7*X6+3*X3**2*X5+3*X3*X4*X5+3*X3*X10*X5+
     &  X3*X7*X5+3*X4*X10*X5+3*X4*X7*X5-X4*X6*X5)
        FM(6,6)=64*PQ**6+16*PQ**4*PH**2+32*PQ**4*(X1+X2+2*X4+X9+X7
     &  +2*X5)+8*PQ**2*PH**2*(-X1+2*X4-X7)+16*PQ**2*(X2*X5-2*X4*
     &  X9-2*X4*X7+4*X4*X5+X9*X5)+8*PH**2*X4*X7-16*X4*X9*X5
        FM(6,7)=8*PQ**4*(-6*X3-3*X4-3*X10-2*X9-X8-X7-2*X6)+2*PQ**2
     &  *PH**2*(-2*X3-X4-2*X10+X7+2*X6)+4*PQ**2*(-8*X1*X3-4*X1*X4
     &  -4*X1*X10+2*X1*X9-2*X1*X8-10*X2*X3-2*X2*X4-5*X2*X10-X2*X9
     &  -2*X2*X8-4*X2*X7-2*X2*X6-5*X3*X9-4*X3*X7-8*X3*X5-2*X4*X9+
     &  7*X4*X8-4*X4*X7+8*X4*X6-4*X4*X5-5*X10*X5-X9**2+X9*X8-2*X9
     &  *X7+X9*X6-2*X9*X5+X8*X7-X8*X5)+2*PH**2*(X1*X10-X3*X7+2*X4
     &  *X7-X4*X6)+4*(2*X1*X2*X9+X1*X2*X8+X1*X9**2-X1*X9*X8-2*X2
     &  **2*X7-X2**2*X6-3*X2*X3*X5-2*X2*X10*X5-X2*X9*X7-X2*X9*X6+
     &  2*X2*X8*X7-X3*X9*X5-X4*X9*X5+2*X4*X8*X5)
        FM(6,8)=8*PQ**4*(-6*X3-3*X4-3*X10-X9-2*X8-X7-2*X6)+2*PQ**2
     &  *PH**2*(-6*X3-3*X4-3*X10+X7+2*X6)+4*PQ**2*(-8*X1*X3-4*X1*
     &  X4-4*X1*X10-8*X2*X3-4*X2*X4-4*X2*X10-4*X3*X9-4*X3*X7-12*
     &  X3*X5-4*X4*X9+8*X4*X8-4*X4*X7+8*X4*X6-6*X4*X5-6*X10*X5-X9
     &  *X5-2*X8*X5)+4*PH**2*(2*X1*X3+X1*X4+X1*X10+X3*X7+X4*X7-2*
     &  X4*X6)+8*X5*(-2*X2*X3-X2*X4-X2*X10-X3*X9-X4*X9+2*X4*X8)
        FM(7,7)=72*PQ**4*X10+18*PQ**2*PH**2*X10+8*PQ**2*(X1*X10+9*
     &  X2*X10+7*X3*X7+2*X3*X6+2*X4*X7+7*X4*X6+X10*X5+2*X9*X7+7*
     &  X9*X6+7*X8*X7+2*X8*X6)+2*PH**2*(-X1*X10-7*X3*X7-2*X3*X6-2
     &  *X4*X7-7*X4*X6)+4*X2*(X10*X5+2*X9*X7+7*X9*X6+7*X8*X7+2*X8
     &  *X6)
        FM(7,8)=72*PQ**4*X10+2*PQ**2*PH**2*X10+4*PQ**2*(2*X1*X10+
     &  10*X2*X10+7*X3*X9+2*X3*X8+14*X3*X7+4*X3*X6+2*X4*X9+7*X4*
     &  X8+4*X4*X7+14*X4*X6+10*X10*X5+X9**2+7*X9*X8+2*X9*X7+7*X9*
     &  X6+X8**2+7*X8*X7+2*X8*X6)+2*PH**2*(7*X1*X10-7*X3*X7-2*X3*
     &  X6-2*X4*X7-7*X4*X6)+2*(-2*X1*X9**2-14*X1*X9*X8-2*X1*X8**2
     &  +2*X2*X10*X5+2*X2*X9*X7+7*X2*X9*X6+7*X2*X8*X7+2*X2*X8*X6+
     &  7*X3*X9*X5+2*X3*X8*X5+2*X4*X9*X5+7*X4*X8*X5)
        FM(8,8)=72*PQ**4*X10+18*PQ**2*PH**2*X10+8*PQ**2*(X1*X10+X2
     &  *X10+7*X3*X9+2*X3*X8+7*X3*X7+2*X3*X6+2*X4*X9+7*X4*X8+2*X4
     &  *X7+7*X4*X6+9*X10*X5)+2*PH**2*(-X1*X10-7*X3*X7-2*X3*X6-2*
     &  X4*X7-7*X4*X6)+4*X5*(X2*X10+7*X3*X9+2*X3*X8+2*X4*X9+7*X4*
     &  X8)
        FM(9,9)=-4*PQ**4*X10-PQ**2*PH**2*X10+4*PQ**2*(-X1*X10-X2*X10+
     &  X3*X7+X4*X6-X10*X5+X9*X6+X8*X7)+PH**2*(X1*X10-X3*X7-X4*X6
     &  )+2*X2*(-X10*X5+X9*X6+X8*X7)
        FM(9,10)=-4*PQ**4*X10-PQ**2*PH**2*X10+2*PQ**2*(-2*X1*X10-2*X2*
     &  X10+2*X3*X9+2*X3*X7+2*X4*X6-2*X10*X5+X9*X8+2*X8*X7)+PH**2
     &  *(X1*X10-X3*X7-X4*X6)+2*(-X1*X9*X8-X2*X10*X5+X2*X8*X7+X3*
     &  X9*X5)
        FMXX=-4*PQ**4*X10-PQ**2*PH**2*X10+2*PQ**2*(-2*X1*X10-2*X2*
     &  X10+2*X4*X8+2*X4*X6+2*X3*X7-2*X10*X5+X9*X8+2*X9*X6)+PH**2
     &  *(X1*X10-X3*X7-X4*X6)+2*(-X1*X9*X8-X2*X10*X5+X2*X9*X6+X4*
     &  X8*X5)
        FM(9,10)=0.5D0*(FMXX+FM(9,10))
        FM(10,10)=-4*PQ**4*X10-PQ**2*PH**2*X10+4*PQ**2*(-X1*X10-X2*X10+
     &  X3*X7+X4*X6-X10*X5+X9*X3+X8*X4)+PH**2*(X1*X10-X3*X7-X4*X6
     &  )+2*X5*(-X10*X2+X9*X3+X8*X4)
 
C...Repackage matrix elements.
        DO 200 I=1,8
          DO 190 J=I,8
            RM(I,J)=FM(I,J)
  190     CONTINUE
  200   CONTINUE
        RM(7,7)=FM(7,7)-2D0*FM(9,9)
        RM(7,8)=FM(7,8)-2D0*FM(9,10)
        RM(8,8)=FM(8,8)-2D0*FM(10,10)
 
C...Produce final result: matrix elements * colours * propagators.
        DO 220 I=1,8
          DO 210 J=I,8
            FAC=8D0
            IF(I.EQ.J)FAC=4D0
            WTQQBH=WTQQBH+RM(I,J)*FAC*CLR(I,J)/(DX(I)*DX(J))
  210     CONTINUE
  220   CONTINUE
        WTQQBH=-WTQQBH/256D0
 
      ELSE
C...Evaluate matrix elements for q + qbar -> Q + Qbar + H.
        A11=-8D0*PQ**4*X10-2D0*PQ**2*PH**2*X10-(8D0*PQ**2)*(X2*X10+X3
     &  *X7+X4*X6+X9*X6+X8*X7)+2D0*PH**2*(X3*X7+X4*X6)-(4D0*X2)*(X9
     &  *X6+X8*X7)
        A12=-8D0*PQ**4*X10+4D0*PQ**2*(-X2*X10-X3*X9-2D0*X3*X7-X4*X8-
     &  2D0*X4*X6-X10*X5-X9*X8-X9*X6-X8*X7)+2D0*PH**2*(-X1*X10+X3*X7
     &  +X4*X6)+2D0*(2D0*X1*X9*X8-X2*X9*X6-X2*X8*X7-X3*X9*X5-X4*X8*
     &  X5)
        A22=-8D0*PQ**4*X10-2D0*PQ**2*PH**2*X10-(8D0*PQ**2)*(X3*X9+X3*
     &  X7+X4*X8+X4*X6+X10*X5)+2D0*PH**2*(X3*X7+X4*X6)-(4D0*X5)*(X3
     &  *X9+X4*X8)
 
C...Produce final result: matrix elements * propagators.
        A11=A11/DX(7)**2
        A12=A12/(DX(7)*DX(8))
        A22=A22/DX(8)**2
        WTQQBH=-(A11+A22+2D0*A12)*8D0/9D0
      ENDIF
 
      RETURN
      END
