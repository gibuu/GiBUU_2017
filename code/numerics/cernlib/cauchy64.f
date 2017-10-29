***************************************************************************************************************************
* This subroutine is part of CERNLIB
* (C) Copyright CERN except where explicitly stated otherwise. 
* Permission to use and/or redistribute this work is granted under the terms of the GNU General Public License. 
***************************************************************************************************************************
*
* $Id: cauchy64.F,v 1.1.1.1 1996/04/01 15:02:14 mclareni Exp $
*
* $Log: cauchy64.F,v $
* Revision 1.1.1.1  1996/04/01 15:02:14  mclareni
* Mathlib gen
*
*
      FUNCTION DCAUCH(F,A,B,S,EPS)
*
* $Id: imp64.inc,v 1.1.1.1 1996/04/01 15:02:59 mclareni Exp $
*
* $Log: imp64.inc,v $
* Revision 1.1.1.1  1996/04/01 15:02:59  mclareni
* Mathlib gen
*
*
* imp64.inc
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER NAME*(*)
      CHARACTER*80 ERRTXT
      PARAMETER (NAME = 'DCAUCH')
      EXTERNAL F
      DIMENSION X(12),W(12)
      PARAMETER (Z1 = 1, HF = Z1/2, CST = 5*Z1/1000)
      DATA X( 1) /9.6028985649753623D-1/, W( 1) /1.0122853629037626D-1/
      DATA X( 2) /7.9666647741362674D-1/, W( 2) /2.2238103445337447D-1/
      DATA X( 3) /5.2553240991632899D-1/, W( 3) /3.1370664587788729D-1/
      DATA X( 4) /1.8343464249564980D-1/, W( 4) /3.6268378337836198D-1/
      DATA X( 5) /9.8940093499164993D-1/, W( 5) /2.7152459411754095D-2/
      DATA X( 6) /9.4457502307323258D-1/, W( 6) /6.2253523938647893D-2/
      DATA X( 7) /8.6563120238783174D-1/, W( 7) /9.5158511682492785D-2/
      DATA X( 8) /7.5540440835500303D-1/, W( 8) /1.2462897125553387D-1/
      DATA X( 9) /6.1787624440264375D-1/, W( 9) /1.4959598881657673D-1/
      DATA X(10) /4.5801677765722739D-1/, W(10) /1.6915651939500254D-1/
      DATA X(11) /2.8160355077925891D-1/, W(11) /1.8260341504492359D-1/
      DATA X(12) /9.5012509837637440D-2/, W(12) /1.8945061045506850D-1/
      IF(S .EQ. A .OR. S .EQ. B) THEN
       H=0
       WRITE(ERRTXT,101) S
       CALL MTLPRT(NAME,'D104.1',ERRTXT)
      ELSEIF(S .LT. MIN(A,B) .OR. S .GT. MAX(A,B)) THEN
       H=DGAUSS(F,A,B,EPS)
      ELSE
       IF(2*S .LE. A+B) THEN
        H=DGAUSS(F,2*S-A,B,EPS)
        B0=S-A
       ELSE
        H=DGAUSS(F,A,2*S-B,EPS)
        B0=B-S
       ENDIF
       C=CST/B0
       BB=0
    1  AA=BB
       BB=B0
    2  C1=HF*(BB+AA)
       C2=HF*(BB-AA)
       C3=S+C1
       C4=S-C1
       S8=0
       DO 3 I = 1,4
       U=C2*X(I)
    3  S8=S8+W(I)*((F(C3+U)+F(C4-U))+(F(C3-U)+F(C4+U)))
       S8=C2*S8
       S16=0
       DO 4 I = 5,12
       U=C2*X(I)
    4  S16=S16+W(I)*((F(C3+U)+F(C4-U))+(F(C3-U)+F(C4+U)))
       S16=C2*S16
       IF(ABS(S16-S8) .LE. EPS*(1+ABS(S16))) GO TO 5
       BB=C1
       IF(1+ABS(C*C2) .NE. 1) GO TO 2
       H=0
       CALL MTLPRT(NAME,'D104.2','TOO HIGH ACCURACY REQUIRED')
       GO TO 9
    5  H=H+S16
       IF(BB .NE. B0) GO TO 1
      END IF
    9 DCAUCH=H
      RETURN
  101 FORMAT('SINGULARITY  S = ',D15.8,' AT END-POINT OF INTERVAL')
      END
