C     -----**************************************************-----
C     -----*                                                *-----
C     -----*     STATISTICAL MODEL WITH UNITARY SYMMETRY    *-----
C     -----*       MAIN PROGRAM FOR  ANTI-N + N  AND        *-----
C     -----*          ANTI-N + (NN)  ANNIHILATION           *-----
C     -----*     BELOW 10 GEV/C FOR MULTIPLICITIES 2<N<6    *-----
C     -----*          AUTHOR: I.A.PSHENICHNOV (INR,MOSCOW)  *-----
C     -----*          VERSION: JULY, 1992                   *-----
C     -----*     Reference: Nucl. Phys. A537, 393 (1992)    *-----
C     -----*     Current address: Frankfurt Institute for   *-----
C     -----*                      Advanced Studies (FIAS)   *-----
C     -----*     Contact tel.: ++069/798-47624              *-----
C     -----*     e-mail:  pshenich@fias.uni-frankfurt.de    *-----
C     -----**************************************************-----

C     input:  PLAB IN GEV/C  ------------------------


      SUBROUTINE UM_TN(PLAB)

!     CHARACTER*7 MARK_LINE, CURRENT_LINE
!     CHARACTER*43 MEANING
      CHARACTER*5 GROUP_NAME

      LOGICAL TWO_NUCLEON, FITTING
      logical verbose
      CHARACTER*4 SYMB,SYC

      INTEGER NCALLS
      DATA NCALLS/0/

      COMMON /MES/SYMB(44),SMM(44),ISM(44),ISM3(44),IY(44),
     &     JSP(44),IG(44),IP(44),IC(44),ICH(44),SYC(44)
      COMMON /CHGR/ICG(6,6000,3),NSUPL(3)
      COMMON /WFFW/FWP(6000),FWP_N(6000),FWN(6000)
      COMMON /PARSWT/IPSWT
      COMMON /SWITCH/TWO_NUCLEON, FITTING
      COMMON /ASB/ASB(9)
      COMMON /VERB/ verbose

      SAVE  /MES/, /CHGR/, /WFFW/, /PARSWT/, /SWITCH/, /ASB/, /VERB/
      SAVE GROUP_NAME,NCALLS,CMAPI,CLAMSP,CLAMSN

      NCALLS=NCALLS+1

      IF(NCALLS.EQ.1) THEN

C     ---------   CMAPI - INTERACTION VOLUME CONSTANT IN GEV  -------

         CMAPI=0.540

C     ----------------------  TWO-NUCLEON TOGGLE --------------------

         TWO_NUCLEON=.FALSE.

C****
         FITTING=.FALSE.

C     ------------------ SYMMETRY BREAKING CONSTANTS ----------------

c**** Pion:
         ASB(1)=1.
c**** Kaon:
         ASB(2)=0.07
c**** Eta:
         ASB(3)=0.13
c**** Rho:
         ASB(4)=0.24
c**** K*:
         ASB(5)=0.05
c**** Omega:
         ASB(6)=0.18
c**** Nucleon:
         ASB(7)=1.00
c**** Sigma,Lambda:
         ASB(8)=1.00
c**** Ksi:
         ASB(9)=1.00

C     ------------- NAME FOR A GROUP OF THE OUTPUT FILES  -----------
         GROUP_NAME='um0_0'

         if(verbose) OPEN(13,FILE=GROUP_NAME//'.out1',STATUS='UNKNOWN')
         IF (TWO_NUCLEON) THEN
           if(verbose) OPEN(16,FILE=GROUP_NAME//'.ppp',STATUS='UNKNOWN')
           if(verbose) OPEN(17,FILE=GROUP_NAME//'.ppn',STATUS='UNKNOWN')
           if(verbose) OPEN(18,FILE=GROUP_NAME//'.pnn',STATUS='UNKNOWN')
         ELSE
           if(verbose) OPEN(14,FILE=GROUP_NAME//'.fnp',STATUS='UNKNOWN')
           if(verbose) OPEN(15,FILE=GROUP_NAME//'.fnn',STATUS='UNKNOWN')
         END IF

C     ---------   IPSWT - NUMBER OF VARIANT FOR INITIAL SET UP ------
C     OF PARAMETERS, WHEN FITTING=.T.   (IPSWT.LE.10)
         IPSWT=2
C     ---------------------------------------------------------------

         CALL MESUS

         CALL USUN

         close(13)

      END IF

C     ----- CLAMSP, CLAMSN - STRANGNESS SUPRESSION FACTORS        -------
c     AL         FOR PBAR+P AND PBAR+N

c     AL***       Standard values: *************************************
c     CLAMSP=1.
c     CLAMSN=1.
c     AL***      Momentum dependent strangeness suppression factors:
      CLAMSP=0.50*exp(-0.16*PLAB)+1.
      CLAMSN=CLAMSP
      CLAMSP=4.2*exp(-1.4*PLAB) + CLAMSP
      CLAMSN=2.0*exp(-0.7*PLAB) + CLAMSN
c******************************************************************


      if(verbose) OPEN(13,FILE=GROUP_NAME//'.out2',STATUS='UNKNOWN')
      if(verbose) WRITE(13,1)PLAB,CMAPI,CLAMSP,CLAMSN
 1    FORMAT(72('*')//8X,'UNITARY-SYMMETRICAL THEORY OF ANTI-N+N',
     &        ' AND ANTI-N+(NN) '/
     &        8X,'ANNIHILATION, VERSION: JULY, 1992'//72('*')///20X,
     &        'MAIN CONSTANTS ARE:'/15X,'PLAB=',F8.3,' GEV/C'/
     &        15X,'INTERACTION VOLUME CONSTANT, CMAPI=',F8.3,' GEV'/
     &        15X,'STRANGENESS SUPRESSION FACTORS, CLAMSP=',F8.3,
     &        ', CLAMSN=',F8.3//)

      CALL ANNIHIL(PLAB,CMAPI,CLAMSP,CLAMSN)

      close(13)

      CALL RESET

      IF(NCALLS.EQ.1) THEN

         IF (TWO_NUCLEON) THEN

            if(verbose) WRITE(16,1)PLAB,CMAPI,CLAMSP,CLAMSN
            if(verbose) WRITE(17,1)PLAB,CMAPI,CLAMSP,CLAMSN
            if(verbose) WRITE(18,1)PLAB,CMAPI,CLAMSP,CLAMSN
            if(verbose) WRITE(16,2)ASB
            if(verbose) WRITE(17,2)ASB
            if(verbose) WRITE(18,2)ASB
 2       FORMAT(/15X,'SYMMETRY BREAKING CONSTANTS ARE:'/
     &       1X,'   PI   ','    K    ','   ETA  ','   RO   ','   K*   ',
     &       1X,'  OMEGA ',' NUCLEON ','SIG,LAMB','  KSI   '/1X,9F8.3)
C     -------------------------   FINAL RESULTS ----------------------
            if(verbose) WRITE(16,16)
 16   FORMAT(///73('*')/5X,'COMPLETE LIST OF ANNIHILATION CHANNELS',
     &     '  WITH  PP-PAIR (%)'/
     &     73('*')/1X,'NUMBER',9X,'IDENTIFIERS',12X,'S Y M B O L S',11X,
     &     'WITH PP')
            NUMERO=0
            DO I=1,NSUPL(1)
               IF (FWP(I).GT.9.999999E-05) THEN
                  NUMERO=NUMERO+1
                  if(verbose) WRITE(16,405)NUMERO,(ICG(L1,I,1),L1=1,6),
     &                 (SYMB(ICG(L2,I,1)),L2=1,6),FWP(I)
               END IF
            END DO
            CLOSE(16)

            if(verbose) WRITE(17,17)
 17   FORMAT(///73('*')/5X,'COMPLETE LIST OF ANNIHILATION CHANNELS',
     &  '  WITH  PN-PAIR (%)'/
     &  73('*')/1X,'NUMBER',9X,'IDENTIFIERS',12X,'S Y M B O L S',11X,
     &  'WITH PN')
            NUMERO=0
            DO I=1,NSUPL(2)
               IF (FWP_N(I).GT.9.999999E-05) THEN
                  NUMERO=NUMERO+1
                  if(verbose) WRITE(17,405)NUMERO,(ICG(L1,I,2),L1=1,6),
     &                 (SYMB(ICG(L2,I,2)),L2=1,6),FWP_N(I)
               END IF
            END DO
            CLOSE(17)

            if(verbose) WRITE(18,18)
 18   FORMAT(///73('*')/5X,'COMPLETE LIST OF ANNIHILATION CHANNELS',
     &           '  WITH  NN-PAIR (%)'/
     &           73('*')/1X,'NUMBER',9X,
     &           'IDENTIFIERS',12X,'S Y M B O L S',11X,'WITH NN')
            NUMERO=0
            DO I=1,NSUPL(3)
               IF (FWN(I).GT.9.999999E-05) THEN
                  NUMERO=NUMERO+1
                  if(verbose) WRITE(18,405)NUMERO,(ICG(L3,I,3),L3=1,6),
     &                 (SYMB(ICG(L4,I,3)),L4=1,6),FWN(I)
               END IF
            END DO
            CLOSE(18)

         ELSE

            if(verbose) WRITE(14,1)PLAB,CMAPI,CLAMSP,CLAMSN
            if(verbose) WRITE(15,1)PLAB,CMAPI,CLAMSP,CLAMSN
            if(verbose) WRITE(14,2)ASB
            if(verbose) WRITE(15,2)ASB

C     -------------------------   FINAL RESULTS ----------------------
            if(verbose) WRITE(14,703)
 703  FORMAT(///73('*')/5X,'COMPLETE LIST OF ANNIHILATION CHANNELS',
     & '  WITH  P R O T O N (%)'/
     & 73('*')/1X,'NUMBER',9X,'IDENTIFIERS',12X,'S Y M B O L S',11X,
     & 'WITH P')
            NUMERO=0
            DO I=1,NSUPL(1)
               IF (FWP(I).GT.9.999999E-05) THEN
                  NUMERO=NUMERO+1
                  if(verbose) WRITE(14,405)NUMERO,(ICG(L1,I,1),L1=1,6),
     &                 (SYMB(ICG(L2,I,1)),L2=1,6),FWP(I)
 405              FORMAT(1x,I5,1x,6(1X,I3),1x,6(1X,A4),1x,F8.4)
               END IF
            END DO
            CLOSE(14)

            if(verbose) WRITE(15,406)
 406  FORMAT(///73('*')/5X,'COMPLETE LIST OF ANNIHILATION CHANNELS',
     &     '  WITH  N E U T R O N (%)'/
     &     73('*')/1X,'NUMBER',9X,
     &     'IDENTIFIERS',12X,'S Y M B O L S',11X,'WITH N')
            NUMERO=0
            DO I=1,NSUPL(3)
               IF (FWN(I).GT.9.999999E-05) THEN
                  NUMERO=NUMERO+1
                  if(verbose) WRITE(15,405)NUMERO,(ICG(L3,I,3),L3=1,6),
     &                 (SYMB(ICG(L4,I,3)),L4=1,6),FWN(I)
               END IF
            END DO
            CLOSE(15)

         END IF

      END IF

      END


      BLOCK DATA ONE
C     ---------- INITIAL SETS OF PARAMETERS: FROM 1 TO 10 -----------
      COMMON /GRPAR/SA(10,6),SPL(10,6),SAMX(10,6),SAMN(10,6)
      SAVE /GRPAR/
C     __1__ __2__ __3__ __4__ __5__ __6__ __7__ __8__ __9__ __10_
      DATA   SA/
     1     1.00, 1.00, 0.95, 0.68, 0.75, 0.80, 1.00, 1.00, 1.00, 1.00,
     2     1.00, 0.00, 0.18, 0.10, 0.10, 0.12, 1.00, 0.40, 1.00, 0.40,
     3     1.00, 0.15, 0.13, 0.10, 0.05, 0.07, 1.00, 1.00, 1.00, 1.00,
     4     1.00, 0.40, 0.35, 0.27, 0.30, 0.40, 0.00, 2.00, 2.00, 2.00,
     5     1.00, 0.00, 0.10, 0.10, 0.10, 0.12, 1.10, 1.00, 1.00, 1.00,
     6     1.00, 0.24, 0.21, 0.16, 0.20, 0.22, 0.00, 0.80, 0.80, 0.80/
C     __1__ __2__ __3__ __4__ __5__ __6__ __7__ __8__ __9__ __10_
      DATA  SPL/
     1     0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     2     0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     3     0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     4     0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     5     0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     6     0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00/
C     __1__ __2__ __3__ __4__ __5__ __6__ __7__ __8__ __9__ __10_
      DATA SAMX/
     1     99.00,99.00,99.00,99.00,99.00,99.00,99.00,99.00,99.00,99.00,
     2     99.00,99.00,99.00,99.00,99.00,99.00,99.00,99.00,99.00,99.00,
     3     99.00,99.00,99.00,99.00,99.00,99.00,99.00,99.00,99.00,99.00,
     4     99.00,99.00,99.00,99.00,99.00,99.00,99.00,99.00,99.00,99.00,
     5     99.00,99.00,99.00,99.00,99.00,99.00,99.00,99.00,99.00,99.00,
     6     99.00,99.00,99.00,99.00,99.00,99.00,99.00,99.00,99.00,99.00/
C     __1__ __2__ __3__ __4__ __5__ __6__ __7__ __8__ __9__ __10_
      DATA SAMN/
     1     0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     2     0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     3     0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     4     0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     5     0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     6     0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00/
C     __1__ __2__ __3__ __4__ __5__ __6__ __7__ __8__ __9__ __10_
      END


      SUBROUTINE RESET
C     ------ THE WORK WAS COMPLETED. TO CALCULATE ALL PRODUCTS OF THE  ---
C     - WEIGHTS AND PARAMETERS AND TOTAL SUM OF THEM. TO NORMALISE ... ---
      LOGICAL TWO_NUCLEON, FITTING
      COMMON /CHGR/ICG(6,6000,3),NSUPL(3)
      COMMON /WFFW/FWP(6000),FWP_N(6000),FWN(6000)
      COMMON /SWITCH/TWO_NUCLEON, FITTING
      SAVE /CHGR/, /WFFW/, /SWITCH/

      CALL ASS
      IF (TWO_NUCLEON) THEN
         K=1
      ELSE
         K=2
      END IF
      DO ISWT=1,3,K
         AIN=0.
         DO I=1,NSUPL(ISWT)
            IF (ISWT.EQ.1) THEN
               FWP(I)=WMUL(I,ISWT)
               AIN=AIN+FWP(I)
            ELSE IF (ISWT.EQ.2) THEN
               FWP_N(I)=WMUL(I,ISWT)
               AIN=AIN+FWP_N(I)
            ELSE
               FWN(I)=WMUL(I,ISWT)
               AIN=AIN+FWN(I)
            END IF
         END DO

         DO I=1,NSUPL(ISWT)
            IF (ISWT.EQ.1) THEN
               FWP(I)=100.*FWP(I)/AIN
            ELSE IF (ISWT.EQ.2) THEN
               FWP_N(I)=100.*FWP_N(I)/AIN
            ELSE
               FWN(I)=100.*FWN(I)/AIN
            END IF
         END DO
      END DO

      RETURN
      END


      SUBROUTINE ASS
C     ------------------- TO FILL THE  FTLAM ARRAY ---------------
      DOUBLE PRECISION A(100)
      LOGICAL TWO_NUCLEON, FITTING
      COMMON /A/A
      COMMON /FTL/FTLAM(44)
      COMMON /SBST/ISBST(44)
      COMMON /SWITCH/TWO_NUCLEON, FITTING
      COMMON /ASB/ASB(9)
      SAVE /A/, /FTL/, /SBST/, /SWITCH/, /ASB/
      IF (FITTING) THEN
         DO I=1,44
            FTLAM(I)=A(ISBST(I))
         END DO
         DO K=1,9
            ASB(K)=A(K)
         END DO
      ELSE
         DO I=1,44
            FTLAM(I)=ASB(ISBST(I))
         END DO
      END IF
C     ---------------------- DUMMY PARAMETER ----------------------
      FTLAM(19)=1.
C     ----------------- EXCEPTION FOR ETA-P AND PHI ---------------
      FTLAM(9)=1.
      FTLAM(18)=1.
      RETURN
      END


      BLOCK DATA TWO
C     ----------- TO ASSIGN CORRESPONDENSE: MESON AND PARAMETER ------
      COMMON /SBST/ISBST(44)
      SAVE /SBST/
      DATA ISBST/1, 1, 2, 2, 2, 2, 1, 3, 1,
     &     4, 4, 5, 5, 5, 5, 4, 6, 1,
     &     1, 17*1,
     &     7, 7, 8, 8, 8, 9, 9, 8/
      END

      REAL FUNCTION WMUL(NBRN,ISWT)
C     ----- CALCULATION THE PRODUCT OF WEIGHTS AND PARAMETERS -----
      LOGICAL TWO_NUCLEON, FITTING
      COMMON /FTL/FTLAM(44)
      COMMON /WEITHS/WCP(4,1000),WCN(1000),WPG(4,6000),WNG(6000)
      COMMON /CHGR/ICG(6,6000,3),NSUPL(3)
      COMMON /SWITCH/TWO_NUCLEON, FITTING
      SAVE /FTL/, /WEITHS/, /CHGR/, /SWITCH/

      IF (ISWT.EQ.1) THEN
         W=WPG(1,NBRN)+WPG(2,NBRN)
      ELSE IF (ISWT.EQ.2) THEN
         W=WPG(3,NBRN)+WPG(4,NBRN)
      ELSE
         W=WNG(NBRN)
      END IF

      DO I=1,6
         W=W*FTLAM(ICG(I,NBRN,ISWT))
      END DO
      WMUL=W
      RETURN
      END


      SUBROUTINE ANNIHIL(PLAB,CMAPI,CLAMSP,CLAMSN)

      use um1, only: SN

      CHARACTER*4 SYMB,SYC,SYMPAT
      DIMENSION VN(8),AMI(21),INDPAT(8),SYMPAT(8),PVOL(1000),SNS(1000),
     &     ISTART(8),IFINISH(8),SUMTOP(3,8),IPT6(6)
      LOGICAL TWO_NUCLEON, FITTING
      logical verbose
      INTEGER BM2, BM3, BM4, BM5, BM6
      INTEGER BARIONM2, BARIONM3, BARIONM4, BARIONM5, BARIONM6
      COMMON /WEITHS/WCP(4,1000),WCN(1000),WPG(4,6000),WNG(6000)
      COMMON /CHGR/ICG(6,6000,3),NSUPL(3)
      COMMON /MES/SYMB(44),SMM(44),ISM(44),ISM3(44),IY(44),JSP(44),
     &     IG(44),IP(44),IC(44),ICH(44),SYC(44)
      COMMON /UN/ I2(2,4,2), IY2(2,4,2), U2(4,4,2),
     &     I3(3,6,2), IY3(3,6,2), U3(4,6,2),
     &     I4(4,9,2), IY4(4,9,2), U4(4,9,2),
     &     I5(5,12,2),IY5(5,12,2),U5(4,12,2),
     &     I6(6,16,2),IY6(6,16,2),U6(4,16,2),
     &     USUM(4,2,6), ISOS(2,6), IYSTR(2,6), NUMLIN(6)
      COMMON /MULB/ BM2(2,4),BM3(3,8),BM4(4,16),BM5(5,32),BM6(6,64)
      COMMON /AK/AKPQ(6,8),AMAN,ICHN(6,1000),NCHAN,
     &     IEYE(44,600,3),NEYE(3)
      COMMON /FILIN/FILIN(5),NII(3)
      COMMON /WFFW/FWP(6000),FWP_N(6000),FWN(6000)
      COMMON /SWITCH/TWO_NUCLEON,FITTING
      COMMON /MULBAR/ BARIONM2(2,4),BARIONM3(3,12),BARIONM4(4,32),
     &     BARIONM5(5,80),BARIONM6(6,192)
      COMMON /VERB/ verbose

      SAVE /WEITHS/, /CHGR/, /MES/, /UN/, /MULB/, /AK/
      SAVE /FILIN/, /WFFW/, /SWITCH/, /MULBAR/, /VERB/

      V0=4.*3.1415926*0.333333*(1./CMAPI)**3

      IF (TWO_NUCLEON) THEN
         E0=SQRT(5*AMAN**2+4*SQRT(PLAB**2+AMAN**2)*AMAN)
         S_AT_REST=9*AMAN**2
         BLAM=4.*AMAN*E0/(E0**2+3.*AMAN**2)
      ELSE
         E0=SQRT(2*AMAN**2+2*SQRT(PLAB**2+AMAN**2)*AMAN)
         S_AT_REST=4*AMAN**2
         BLAM=2.*AMAN/E0
      END IF
      S=E0**2
      VT0=V0*(S_AT_REST*ALOG(S))/(S*ALOG(S_AT_REST))

      if(verbose) WRITE(13,1)E0,V0,VT0,BLAM
 1    FORMAT(////3X,'TOTAL ENERGY=',F8.3,' GEV'/5X,
     &     'V0=',F8.3,' 1/GEV**3',5X,'VT0=',F8.3,' 1/GEV**3',
     &     3X,'BLAM=',F8.3)

      V0=VT0
      SUMP=0.
      SUMP_N=0.
      SUMN=0.
      NCHAN=0

      DO 90 I=1,3
         NSUPL(I)=0
         DO 91 J=1,8
            SUMTOP(I,J)=0.
  91        CONTINUE
  90      CONTINUE


          DO 702 K=1,3
           DO 701 J=1,6000
            DO 700 I=1,6
              ICG(I,J,K)=19
 700        CONTINUE
 701       CONTINUE
 702      CONTINUE

           DO 801 J=1,1000
            DO 800 I=1,6
              ICHN(I,J)=19
 800        CONTINUE
 801       CONTINUE

C
C -------------- TWO PARTICLES IN THE FINAL STATE ----------------
C
         N=2
         VN(N)=(BLAM*V0)**(N-1)
      if(verbose) WRITE(13,20)VN(2)
  20  FORMAT(72('*')/6X,'ANNIHILATION IN TWO MESONS',6X,'VN(2)=',E10.3,
     &' 1/GEV**3'/72('*'))
      IF (TWO_NUCLEON) THEN
         NCOMB=N*(2**(N-1))
      ELSE
         NCOMB=2**N
      END IF
         ISTART(N)=1
         CALL IENUL
         NLI=NUMLIN(N)
      DO 21 I=1,NLI
         DO 22 ICM=1,NCOMB
         SPN=1.
         PRODMAS=1.
         NKAON=0
            DO 23 IPAT=1,N
              ISMT=I2(IPAT,I,1)
              IYT=IY2(IPAT,I,1)
         IF (TWO_NUCLEON) THEN
              JSPT=BARIONM2(IPAT,ICM)
         ELSE
              JSPT=BM2(IPAT,ICM)
         END IF
                CALL ASG(ISMT,IYT,JSPT,NUMBP)
                INDPAT(IPAT)=NUMBP
                SYMPAT(IPAT)=SYC(NUMBP)
                AMI(IPAT)=SMM(NUMBP)
                PRODMAS=PRODMAS*2.*AMI(IPAT)
            IF (JSPT.GT.4) THEN
                SPN=SPN*(0.2*JSPT+1.)
                NKAON=NKAON+1-IY(NUMBP)
            ELSE
                SPN=SPN*(2.*JSPT+1.)
                NKAON=NKAON+IABS(IY(NUMBP))
            END IF
 23         CONTINUE
             IF(NKAON.EQ.0) THEN
                FFKP=1.
                FFKN=1.
             ELSE
                FFKP=CLAMSP**NKAON
                FFKN=CLAMSN**NKAON
             END IF
          IF(FFKP.LE.1.E-14 .AND. FFKN.LE.1.E-14) GOTO 22
          PVN=SN(E0,N,AMI)*PRODMAS
          IF(PVN.LE.1.E-06) GOTO 22
          CALL GIND(INDPAT,N,GN)
          IF(GN.LE.0.) GOTO 22
                    NCHAN=NCHAN+1
             PVOL(NCHAN)=PVN
             SNS(NCHAN)=SPN
                    DO 24 IPAT=1,N
                      ICHN(IPAT,NCHAN)=INDPAT(IPAT)
 24                 CONTINUE

                       IF (TWO_NUCLEON) THEN
C --------  SUMMATION OVER SU(3), ANNIHILATION WITH PP AND PN ---------
                          DO ISST=1,4
C Supply the remainder of two integer parameters, arg1/arg2
                      MODULUS=MOD(ISST,2)
      WCP(ISST,NCHAN)=
     &      ((AKPQ(2,ISST+3)+AKPQ(3,ISST+3))*U2(1,I,2-MODULUS)
     &      + AKPQ(4,ISST+3)*U2(2,I,2-MODULUS)
     &      + AKPQ(5,ISST+3)*U2(3,I,2-MODULUS)
     &      + AKPQ(6,ISST+3)*U2(4,I,2-MODULUS))*(MODULUS+1)/3.
      WCP(ISST,NCHAN)=VN(N)*SNS(NCHAN)*FFKP*PVOL(NCHAN)*WCP(ISST,NCHAN)
                          END DO
              SUMP = SUMP + WCP(1,NCHAN) + WCP(2,NCHAN)
              SUMP_N = SUMP_N + WCP(3,NCHAN) + WCP(4,NCHAN)
              SUMTOP(1,N)=SUMTOP(1,N) + WCP(1,NCHAN) + WCP(2,NCHAN)
              SUMTOP(2,N)=SUMTOP(2,N) + WCP(3,NCHAN) + WCP(4,NCHAN)

C ------------ SUMMATION OVER SU(3), ANNIHILATION WITH NN -----------
        WCN(NCHAN)=(AKPQ(2,8)+AKPQ(3,8))*U2(1,I,2)
     &  +AKPQ(4,8)*U2(2,I,2)+AKPQ(5,8)*U2(3,I,2)
     &             +AKPQ(6,8)*U2(4,I,2)
      WCN(NCHAN)=VN(N)*SNS(NCHAN)*FFKN*PVOL(NCHAN)*WCN(NCHAN)
                   SUMN=SUMN+WCN(NCHAN)
                   SUMTOP(3,N)=SUMTOP(3,N)+WCN(NCHAN)

                               ELSE
C ---- SUMMATION OVER SU(3) REPRESENTATIONS, ANNIHILATION WITH PROTON ---
       DO 25 ISST=1,2
        WCP(ISST,NCHAN)=AKPQ(1,ISST)*U2(1,I,ISST)
     &             +(AKPQ(2,ISST)+AKPQ(3,ISST))*U2(2,I,ISST)
     &             +(AKPQ(4,ISST)+AKPQ(5,ISST))*U2(3,I,ISST)
     &             +AKPQ(6,ISST)*U2(4,I,ISST)
        WCP(ISST,NCHAN)=VN(N)*SNS(NCHAN)*FFKP*PVOL(NCHAN)*
     &           WCP(ISST,NCHAN)
                   SUMP=SUMP+WCP(ISST,NCHAN)
                   SUMTOP(1,N)=SUMTOP(1,N)+WCP(ISST,NCHAN)
 25    CONTINUE
C -- SUMMATION OVER SU(3) REPRESENTATIONS, ANNIHILATION WITH NEUTERON ---
        WCN(NCHAN)=AKPQ(1,3)*U2(1,I,2)
     &             +(AKPQ(2,3)+AKPQ(3,3))*U2(2,I,2)
     &             +(AKPQ(4,3)+AKPQ(5,3))*U2(3,I,2)
     &             +AKPQ(6,3)*U2(4,I,2)
      WCN(NCHAN)=VN(N)*SNS(NCHAN)*FFKN*PVOL(NCHAN)*WCN(NCHAN)
                   SUMN=SUMN+WCN(NCHAN)
                   SUMTOP(2,N)=SUMTOP(2,N)+WCN(NCHAN)
                                END IF

 22    CONTINUE
 21    CONTINUE
                IFINISH(N)=NCHAN
C
C -------------- THREE PARTICLES IN THE FINAL STATE ----------------
C
         N=3
         VN(N)=(BLAM*V0)**(N-1)
      if(verbose) WRITE(13,30)VN(3)
  30  FORMAT(72('*')/6X,'ANNIHILATION IN THREE MESONS',6X,'VN(3)=',
     &E10.3,' 1/GEV**3'/72('*'))
                         IF (TWO_NUCLEON) THEN
                            NCOMB=N*(2**(N-1))
                         ELSE
                            NCOMB=2**N
                         END IF
         ISTART(N)=NCHAN+1
         CALL IENUL
         NLI=NUMLIN(N)
      DO 31 I=1,NLI
         DO 32 ICM=1,NCOMB
         SPN=1.
         PRODMAS=1.
         NKAON=0
            DO 33 IPAT=1,N
              ISMT=I3(IPAT,I,1)
              IYT=IY3(IPAT,I,1)
         IF (TWO_NUCLEON) THEN
              JSPT=BARIONM3(IPAT,ICM)
         ELSE
              JSPT=BM3(IPAT,ICM)
         END IF
                CALL ASG(ISMT,IYT,JSPT,NUMBP)
                INDPAT(IPAT)=NUMBP
                SYMPAT(IPAT)=SYC(NUMBP)
                AMI(IPAT)=SMM(NUMBP)
                PRODMAS=PRODMAS*2.*AMI(IPAT)
            IF (JSPT.GT.4) THEN
                SPN=SPN*(0.2*JSPT+1.)
                NKAON=NKAON+1-IY(NUMBP)
            ELSE
                SPN=SPN*(2.*JSPT+1.)
                NKAON=NKAON+IABS(IY(NUMBP))
            END IF
 33         CONTINUE
             IF(NKAON.EQ.0) THEN
                FFKP=1.
                FFKN=1.
             ELSE
                FFKP=CLAMSP**NKAON
                FFKN=CLAMSN**NKAON
             END IF
          IF(FFKP.LE.1.E-14 .AND. FFKN.LE.1.E-14) GOTO 32
          PVN=SN(E0,N,AMI)*PRODMAS
          IF(PVN.LE.1.E-06) GOTO 32
          CALL GIND(INDPAT,N,GN)
          IF(GN.LE.0.) GOTO 32
                    NCHAN=NCHAN+1
             PVOL(NCHAN)=PVN
             SNS(NCHAN)=SPN
                    DO 34 IPAT=1,N
                      ICHN(IPAT,NCHAN)=INDPAT(IPAT)
 34                 CONTINUE

                       IF (TWO_NUCLEON) THEN
C --------  SUMMATION OVER SU(3), ANNIHILATION WITH PP AND PN ---------
                          DO ISST=1,4
C Supply the remainder of two integer parameters, arg1/arg2
                      MODULUS=MOD(ISST,2)
      WCP(ISST,NCHAN)=
     &      ((AKPQ(2,ISST+3)+AKPQ(3,ISST+3))*U3(1,I,2-MODULUS)
     &      + AKPQ(4,ISST+3)*U3(2,I,2-MODULUS)
     &      + AKPQ(5,ISST+3)*U3(3,I,2-MODULUS)
     &      + AKPQ(6,ISST+3)*U3(4,I,2-MODULUS))*(MODULUS+1)/3.
      WCP(ISST,NCHAN)=VN(N)*SNS(NCHAN)*FFKP*PVOL(NCHAN)*WCP(ISST,NCHAN)
                          END DO
              SUMP = SUMP + WCP(1,NCHAN) + WCP(2,NCHAN)
              SUMP_N = SUMP_N + WCP(3,NCHAN) + WCP(4,NCHAN)
              SUMTOP(1,N)=SUMTOP(1,N) + WCP(1,NCHAN) + WCP(2,NCHAN)
              SUMTOP(2,N)=SUMTOP(2,N) + WCP(3,NCHAN) + WCP(4,NCHAN)

C ------------ SUMMATION OVER SU(3), ANNIHILATION WITH NN -----------
        WCN(NCHAN)=(AKPQ(2,8)+AKPQ(3,8))*U3(1,I,2)
     &  +AKPQ(4,8)*U3(2,I,2)+AKPQ(5,8)*U3(3,I,2)
     &             +AKPQ(6,8)*U3(4,I,2)
      WCN(NCHAN)=VN(N)*SNS(NCHAN)*FFKN*PVOL(NCHAN)*WCN(NCHAN)
                   SUMN=SUMN+WCN(NCHAN)
                   SUMTOP(3,N)=SUMTOP(3,N)+WCN(NCHAN)

                               ELSE
C ---- SUMMATION OVER SU(3) REPRESENTATIONS, ANNIHILATION WITH PROTON ---
       DO 35 ISST=1,2
        WCP(ISST,NCHAN)=AKPQ(1,ISST)*U3(1,I,ISST)
     &             +(AKPQ(2,ISST)+AKPQ(3,ISST))*U3(2,I,ISST)
     &             +(AKPQ(4,ISST)+AKPQ(5,ISST))*U3(3,I,ISST)
     &             +AKPQ(6,ISST)*U3(4,I,ISST)
        WCP(ISST,NCHAN)=VN(N)*SNS(NCHAN)*FFKP*PVOL(NCHAN)*
     &           WCP(ISST,NCHAN)
                   SUMP=SUMP+WCP(ISST,NCHAN)
                   SUMTOP(1,N)=SUMTOP(1,N)+WCP(ISST,NCHAN)
 35    CONTINUE
C -- SUMMATION OVER SU(3) REPRESENTATIONS, ANNIHILATION WITH NEUTERON ---
        WCN(NCHAN)=AKPQ(1,3)*U3(1,I,2)
     &             +(AKPQ(2,3)+AKPQ(3,3))*U3(2,I,2)
     &             +(AKPQ(4,3)+AKPQ(5,3))*U3(3,I,2)
     &             +AKPQ(6,3)*U3(4,I,2)
      WCN(NCHAN)=VN(N)*SNS(NCHAN)*FFKN*PVOL(NCHAN)*WCN(NCHAN)
                   SUMN=SUMN+WCN(NCHAN)
                   SUMTOP(2,N)=SUMTOP(2,N)+WCN(NCHAN)
                                END IF

 32    CONTINUE
 31    CONTINUE
                IFINISH(N)=NCHAN
C
C -------------- FOUR PARTICLES IN THE FINAL STATE ----------------
C
         N=4
         VN(N)=(BLAM*V0)**(N-1)
      if(verbose) WRITE(13,40)VN(4)
  40  FORMAT(72('*')/6X,'ANNIHILATION IN FOUR MESONS',6X,'VN(4)=',E10.3,
     &' 1/GEV**3'/72('*'))
                         IF (TWO_NUCLEON) THEN
                            NCOMB=N*(2**(N-1))
                         ELSE
                            NCOMB=2**N
                         END IF
         ISTART(N)=NCHAN+1
         CALL IENUL
         NLI=NUMLIN(N)
      DO 41 I=1,NLI
         DO 42 ICM=1,NCOMB
         SPN=1.
         PRODMAS=1.
         NKAON=0
            DO 43 IPAT=1,N
              ISMT=I4(IPAT,I,1)
              IYT=IY4(IPAT,I,1)
         IF (TWO_NUCLEON) THEN
              JSPT=BARIONM4(IPAT,ICM)
         ELSE
              JSPT=BM4(IPAT,ICM)
         END IF
                CALL ASG(ISMT,IYT,JSPT,NUMBP)
                INDPAT(IPAT)=NUMBP
                SYMPAT(IPAT)=SYC(NUMBP)
                AMI(IPAT)=SMM(NUMBP)
                PRODMAS=PRODMAS*2.*AMI(IPAT)
            IF (JSPT.GT.4) THEN
                SPN=SPN*(0.2*JSPT+1.)
                NKAON=NKAON+1-IY(NUMBP)
            ELSE
                SPN=SPN*(2.*JSPT+1.)
                NKAON=NKAON+IABS(IY(NUMBP))
            END IF
 43         CONTINUE
             IF(NKAON.EQ.0) THEN
                FFKP=1.
                FFKN=1.
             ELSE
                FFKP=CLAMSP**NKAON
                FFKN=CLAMSN**NKAON
             END IF
          IF(FFKP.LE.1.E-14 .AND. FFKN.LE.1.E-14) GOTO 42
          PVN=SN(E0,N,AMI)*PRODMAS
          IF(PVN.LE.1.E-06) GOTO 42
          CALL GIND(INDPAT,N,GN)
          IF(GN.LE.0.) GOTO 42
                    NCHAN=NCHAN+1
             PVOL(NCHAN)=PVN
             SNS(NCHAN)=SPN
                    DO 44 IPAT=1,N
                      ICHN(IPAT,NCHAN)=INDPAT(IPAT)
 44                 CONTINUE

                       IF (TWO_NUCLEON) THEN
C --------  SUMMATION OVER SU(3), ANNIHILATION WITH PP AND PN ---------
                          DO ISST=1,4
C Supply the remainder of two integer parameters, arg1/arg2
                      MODULUS=MOD(ISST,2)
      WCP(ISST,NCHAN)=
     &      ((AKPQ(2,ISST+3)+AKPQ(3,ISST+3))*U4(1,I,2-MODULUS)
     &      + AKPQ(4,ISST+3)*U4(2,I,2-MODULUS)
     &      + AKPQ(5,ISST+3)*U4(3,I,2-MODULUS)
     &      + AKPQ(6,ISST+3)*U4(4,I,2-MODULUS))*(MODULUS+1)/3.
      WCP(ISST,NCHAN)=VN(N)*SNS(NCHAN)*FFKP*PVOL(NCHAN)*WCP(ISST,NCHAN)
                          END DO
              SUMP = SUMP + WCP(1,NCHAN) + WCP(2,NCHAN)
              SUMP_N = SUMP_N + WCP(3,NCHAN) + WCP(4,NCHAN)
              SUMTOP(1,N)=SUMTOP(1,N) + WCP(1,NCHAN) + WCP(2,NCHAN)
              SUMTOP(2,N)=SUMTOP(2,N) + WCP(3,NCHAN) + WCP(4,NCHAN)

C ------------ SUMMATION OVER SU(3), ANNIHILATION WITH NN -----------
        WCN(NCHAN)=(AKPQ(2,8)+AKPQ(3,8))*U4(1,I,2)
     &  +AKPQ(4,8)*U4(2,I,2)+AKPQ(5,8)*U4(3,I,2)
     &             +AKPQ(6,8)*U4(4,I,2)
      WCN(NCHAN)=VN(N)*SNS(NCHAN)*FFKN*PVOL(NCHAN)*WCN(NCHAN)
                   SUMN=SUMN+WCN(NCHAN)
                   SUMTOP(3,N)=SUMTOP(3,N)+WCN(NCHAN)

                               ELSE
C ---- SUMMATION OVER SU(3) REPRESENTATIONS, ANNIHILATION WITH PROTON ---
       DO 45 ISST=1,2
        WCP(ISST,NCHAN)=AKPQ(1,ISST)*U4(1,I,ISST)
     &             +(AKPQ(2,ISST)+AKPQ(3,ISST))*U4(2,I,ISST)
     &             +(AKPQ(4,ISST)+AKPQ(5,ISST))*U4(3,I,ISST)
     &             +AKPQ(6,ISST)*U4(4,I,ISST)
        WCP(ISST,NCHAN)=VN(N)*SNS(NCHAN)*FFKP*PVOL(NCHAN)*
     &           WCP(ISST,NCHAN)
                   SUMP=SUMP+WCP(ISST,NCHAN)
                   SUMTOP(1,N)=SUMTOP(1,N)+WCP(ISST,NCHAN)
 45    CONTINUE
C -- SUMMATION OVER SU(3) REPRESENTATIONS, ANNIHILATION WITH NEUTERON ---
        WCN(NCHAN)=AKPQ(1,3)*U4(1,I,2)
     &             +(AKPQ(2,3)+AKPQ(3,3))*U4(2,I,2)
     &             +(AKPQ(4,3)+AKPQ(5,3))*U4(3,I,2)
     &             +AKPQ(6,3)*U4(4,I,2)
      WCN(NCHAN)=VN(N)*SNS(NCHAN)*FFKN*PVOL(NCHAN)*WCN(NCHAN)
                   SUMN=SUMN+WCN(NCHAN)
                   SUMTOP(2,N)=SUMTOP(2,N)+WCN(NCHAN)
                                END IF

 42    CONTINUE
 41    CONTINUE
                IFINISH(N)=NCHAN
C
C -------------- FIVE PARTICLES IN THE FINAL STATE ----------------
C
         N=5
         VN(N)=(BLAM*V0)**(N-1)
      if(verbose) WRITE(13,50)VN(5)
  50  FORMAT(72('*')/6X,'ANNIHILATION IN FIVE MESONS',6X,'VN(5)=',E10.3,
     &' 1/GEV**3'/72('*'))
                         IF (TWO_NUCLEON) THEN
                            NCOMB=N*(2**(N-1))
                         ELSE
                            NCOMB=2**N
                         END IF
         ISTART(N)=NCHAN+1
         CALL IENUL
         NLI=NUMLIN(N)
      DO 51 I=1,NLI
         DO 52 ICM=1,NCOMB
         SPN=1.
         PRODMAS=1.
         NKAON=0
            DO 53 IPAT=1,N
              ISMT=I5(IPAT,I,1)
              IYT=IY5(IPAT,I,1)
         IF (TWO_NUCLEON) THEN
              JSPT=BARIONM5(IPAT,ICM)
         ELSE
              JSPT=BM5(IPAT,ICM)
         END IF
                CALL ASG(ISMT,IYT,JSPT,NUMBP)
                INDPAT(IPAT)=NUMBP
                SYMPAT(IPAT)=SYC(NUMBP)
                AMI(IPAT)=SMM(NUMBP)
                PRODMAS=PRODMAS*2.*AMI(IPAT)
            IF (JSPT.GT.4) THEN
                SPN=SPN*(0.2*JSPT+1.)
                NKAON=NKAON+1-IY(NUMBP)
            ELSE
                SPN=SPN*(2.*JSPT+1.)
                NKAON=NKAON+IABS(IY(NUMBP))
            END IF
 53         CONTINUE
             IF(NKAON.EQ.0) THEN
                FFKP=1.
                FFKN=1.
             ELSE
                FFKP=CLAMSP**NKAON
                FFKN=CLAMSN**NKAON
             END IF
          IF(FFKP.LE.1.E-14 .AND. FFKN.LE.1.E-14) GOTO 52
          PVN=SN(E0,N,AMI)*PRODMAS
          IF(PVN.LE.1.E-06) GOTO 52
          CALL GIND(INDPAT,N,GN)
          IF(GN.LE.0.) GOTO 52
                    NCHAN=NCHAN+1
             PVOL(NCHAN)=PVN
             SNS(NCHAN)=SPN
                    DO 54 IPAT=1,N
                      ICHN(IPAT,NCHAN)=INDPAT(IPAT)
 54                 CONTINUE

                       IF (TWO_NUCLEON) THEN
C --------  SUMMATION OVER SU(3), ANNIHILATION WITH PP AND PN ---------
                          DO ISST=1,4
C Supply the remainder of two integer parameters, arg1/arg2
                      MODULUS=MOD(ISST,2)
      WCP(ISST,NCHAN)=
     &      ((AKPQ(2,ISST+3)+AKPQ(3,ISST+3))*U5(1,I,2-MODULUS)
     &      + AKPQ(4,ISST+3)*U5(2,I,2-MODULUS)
     &      + AKPQ(5,ISST+3)*U5(3,I,2-MODULUS)
     &      + AKPQ(6,ISST+3)*U5(4,I,2-MODULUS))*(MODULUS+1)/3.
      WCP(ISST,NCHAN)=VN(N)*SNS(NCHAN)*FFKP*PVOL(NCHAN)*WCP(ISST,NCHAN)
                          END DO
              SUMP = SUMP + WCP(1,NCHAN) + WCP(2,NCHAN)
              SUMP_N = SUMP_N + WCP(3,NCHAN) + WCP(4,NCHAN)
              SUMTOP(1,N)=SUMTOP(1,N) + WCP(1,NCHAN) + WCP(2,NCHAN)
              SUMTOP(2,N)=SUMTOP(2,N) + WCP(3,NCHAN) + WCP(4,NCHAN)

C ------------ SUMMATION OVER SU(3), ANNIHILATION WITH NN -----------
        WCN(NCHAN)=(AKPQ(2,8)+AKPQ(3,8))*U5(1,I,2)
     &  +AKPQ(4,8)*U5(2,I,2)+AKPQ(5,8)*U5(3,I,2)
     &             +AKPQ(6,8)*U5(4,I,2)
      WCN(NCHAN)=VN(N)*SNS(NCHAN)*FFKN*PVOL(NCHAN)*WCN(NCHAN)
                   SUMN=SUMN+WCN(NCHAN)
                   SUMTOP(3,N)=SUMTOP(3,N)+WCN(NCHAN)

                               ELSE
C ---- SUMMATION OVER SU(3) REPRESENTATIONS, ANNIHILATION WITH PROTON ---
       DO 55 ISST=1,2
        WCP(ISST,NCHAN)=AKPQ(1,ISST)*U5(1,I,ISST)
     &             +(AKPQ(2,ISST)+AKPQ(3,ISST))*U5(2,I,ISST)
     &             +(AKPQ(4,ISST)+AKPQ(5,ISST))*U5(3,I,ISST)
     &             +AKPQ(6,ISST)*U5(4,I,ISST)
        WCP(ISST,NCHAN)=VN(N)*SNS(NCHAN)*FFKP*PVOL(NCHAN)*
     &           WCP(ISST,NCHAN)
                   SUMP=SUMP+WCP(ISST,NCHAN)
                   SUMTOP(1,N)=SUMTOP(1,N)+WCP(ISST,NCHAN)
 55    CONTINUE
C -- SUMMATION OVER SU(3) REPRESENTATIONS, ANNIHILATION WITH NEUTERON ---
        WCN(NCHAN)=AKPQ(1,3)*U5(1,I,2)
     &             +(AKPQ(2,3)+AKPQ(3,3))*U5(2,I,2)
     &             +(AKPQ(4,3)+AKPQ(5,3))*U5(3,I,2)
     &             +AKPQ(6,3)*U5(4,I,2)
      WCN(NCHAN)=VN(N)*SNS(NCHAN)*FFKN*PVOL(NCHAN)*WCN(NCHAN)
                   SUMN=SUMN+WCN(NCHAN)
                   SUMTOP(2,N)=SUMTOP(2,N)+WCN(NCHAN)
                                END IF

 52    CONTINUE
 51    CONTINUE
                IFINISH(N)=NCHAN
C
C -------------- SIX PARTICLES IN THE FINAL STATE ----------------
C
         N=6
         VN(N)=(BLAM*V0)**(N-1)
      if(verbose) WRITE(13,60)VN(6)
  60  FORMAT(72('*')/6X,'ANNIHILATION IN SIX MESONS',6X,'VN(6)=',E10.3,
     &' 1/GEV**3'/72('*'))
                         IF (TWO_NUCLEON) THEN
                            NCOMB=N*(2**(N-1))
                         ELSE
                            NCOMB=2**N
                         END IF
         ISTART(N)=NCHAN+1
         CALL IENUL
         NLI=NUMLIN(N)
      DO 61 I=1,NLI
         DO 62 ICM=1,NCOMB
         SPN=1.
         PRODMAS=1.
         NKAON=0
            DO 63 IPAT=1,N
              ISMT=I6(IPAT,I,1)
              IYT=IY6(IPAT,I,1)
         IF (TWO_NUCLEON) THEN
              JSPT=BARIONM6(IPAT,ICM)
         ELSE
              JSPT=BM6(IPAT,ICM)
         END IF
                CALL ASG(ISMT,IYT,JSPT,NUMBP)
                INDPAT(IPAT)=NUMBP
                SYMPAT(IPAT)=SYC(NUMBP)
                AMI(IPAT)=SMM(NUMBP)
                PRODMAS=PRODMAS*2.*AMI(IPAT)
            IF (JSPT.GT.4) THEN
                SPN=SPN*(0.2*JSPT+1.)
                NKAON=NKAON+1-IY(NUMBP)
            ELSE
                SPN=SPN*(2.*JSPT+1.)
                NKAON=NKAON+IABS(IY(NUMBP))
            END IF
 63         CONTINUE
             IF(NKAON.EQ.0) THEN
                FFKP=1.
                FFKN=1.
             ELSE
                FFKP=CLAMSP**NKAON
                FFKN=CLAMSN**NKAON
             END IF
          IF(FFKP.LE.1.E-14 .AND. FFKN.LE.1.E-14) GOTO 62
          PVN=SN(E0,N,AMI)*PRODMAS
          IF(PVN.LE.1.E-06) GOTO 62
          CALL GIND(INDPAT,N,GN)
          IF(GN.LE.0.) GOTO 62
                    NCHAN=NCHAN+1
             PVOL(NCHAN)=PVN
             SNS(NCHAN)=SPN
                    DO 64 IPAT=1,N
                      ICHN(IPAT,NCHAN)=INDPAT(IPAT)
 64                 CONTINUE

                       IF (TWO_NUCLEON) THEN
C --------  SUMMATION OVER SU(3), ANNIHILATION WITH PP AND PN ---------
                          DO ISST=1,4
C Supply the remainder of two integer parameters, arg1/arg2
                      MODULUS=MOD(ISST,2)
      WCP(ISST,NCHAN)=
     &      ((AKPQ(2,ISST+3)+AKPQ(3,ISST+3))*U6(1,I,2-MODULUS)
     &      + AKPQ(4,ISST+3)*U6(2,I,2-MODULUS)
     &      + AKPQ(5,ISST+3)*U6(3,I,2-MODULUS)
     &      + AKPQ(6,ISST+3)*U6(4,I,2-MODULUS))*(MODULUS+1)/3.
      WCP(ISST,NCHAN)=VN(N)*SNS(NCHAN)*FFKP*PVOL(NCHAN)*WCP(ISST,NCHAN)
                          END DO
              SUMP = SUMP + WCP(1,NCHAN) + WCP(2,NCHAN)
              SUMP_N = SUMP_N + WCP(3,NCHAN) + WCP(4,NCHAN)
              SUMTOP(1,N)=SUMTOP(1,N) + WCP(1,NCHAN) + WCP(2,NCHAN)
              SUMTOP(2,N)=SUMTOP(2,N) + WCP(3,NCHAN) + WCP(4,NCHAN)

C ------------ SUMMATION OVER SU(3), ANNIHILATION WITH NN -----------
        WCN(NCHAN)=(AKPQ(2,8)+AKPQ(3,8))*U6(1,I,2)
     &  +AKPQ(4,8)*U6(2,I,2)+AKPQ(5,8)*U6(3,I,2)
     &             +AKPQ(6,8)*U6(4,I,2)
      WCN(NCHAN)=VN(N)*SNS(NCHAN)*FFKN*PVOL(NCHAN)*WCN(NCHAN)
                   SUMN=SUMN+WCN(NCHAN)
                   SUMTOP(3,N)=SUMTOP(3,N)+WCN(NCHAN)

                               ELSE
C ---- SUMMATION OVER SU(3) REPRESENTATIONS, ANNIHILATION WITH PROTON ---
       DO 65 ISST=1,2
        WCP(ISST,NCHAN)=AKPQ(1,ISST)*U6(1,I,ISST)
     &             +(AKPQ(2,ISST)+AKPQ(3,ISST))*U6(2,I,ISST)
     &             +(AKPQ(4,ISST)+AKPQ(5,ISST))*U6(3,I,ISST)
     &             +AKPQ(6,ISST)*U6(4,I,ISST)
        WCP(ISST,NCHAN)=VN(N)*SNS(NCHAN)*FFKP*PVOL(NCHAN)*
     &           WCP(ISST,NCHAN)
                   SUMP=SUMP+WCP(ISST,NCHAN)
                   SUMTOP(1,N)=SUMTOP(1,N)+WCP(ISST,NCHAN)
 65    CONTINUE
C -- SUMMATION OVER SU(3) REPRESENTATIONS, ANNIHILATION WITH NEUTERON ---
        WCN(NCHAN)=AKPQ(1,3)*U6(1,I,2)
     &             +(AKPQ(2,3)+AKPQ(3,3))*U6(2,I,2)
     &             +(AKPQ(4,3)+AKPQ(5,3))*U6(3,I,2)
     &             +AKPQ(6,3)*U6(4,I,2)
      WCN(NCHAN)=VN(N)*SNS(NCHAN)*FFKN*PVOL(NCHAN)*WCN(NCHAN)
                   SUMN=SUMN+WCN(NCHAN)
                   SUMTOP(2,N)=SUMTOP(2,N)+WCN(NCHAN)
                                END IF

 62    CONTINUE
 61    CONTINUE
                IFINISH(N)=NCHAN
C
C    ----------  NORMALISATION PROCEDURE AND PRINTING    ------------
C
                         IF (TWO_NUCLEON) THEN
      if(verbose) WRITE(13,510)
 510  FORMAT(///100('*')/5X,'NCHAN',8X,'S Y M B O L S',
     &14X,'P.V.*MASS ',5X,'SPIN FAC. ',5X,' WITH PP ',5X,' WITH PN ',
     &5X,' WITH NN ')
       DO I=1,NCHAN
            WCP(1,I)=WCP(1,I)/SUMP
            WCP(2,I)=WCP(2,I)/SUMP
                WCP00=WCP(1,I)+WCP(2,I)
            WCP(3,I)=WCP(3,I)/SUMP_N
            WCP(4,I)=WCP(4,I)/SUMP_N
                WCP_N00=WCP(3,I)+WCP(4,I)
            WCN(I)=WCN(I)/SUMN
            if(verbose) WRITE(13,511)I,(SYC(ICHN(L,I)),L=1,6),PVOL(I),
     &        SNS(I),WCP00,WCP_N00,WCN(I)
 511   FORMAT(4X,I4,2X,6(A4,1X),5X,E10.2,5X,F5.1,5X,3(5X,F6.5,5X))
       END DO
C            -----  CROSS-SECTIONS FOR DIFFERENT N ----
                 DO I=2,6
                    SUMTOP(1,I)=SUMTOP(1,I)/SUMP
                    SUMTOP(2,I)=SUMTOP(2,I)/SUMP_N
                    SUMTOP(3,I)=SUMTOP(3,I)/SUMN
                 END DO
           if(verbose) WRITE(13,512)
 512  FORMAT(//'    TOPOLOGICAL PROBABILITIES WITH PP, PN AND NN'/)
                   DO K=1,3
           if(verbose) WRITE(13,505)(M,SUMTOP(K,M),M=2,6)
 505       FORMAT(5X,5('    N=',I2,' SUM=',F8.4)/)
                   END DO

                               ELSE
      if(verbose) WRITE(13,500)
 500  FORMAT(///100('*')/5X,'NCHAN',8X,'S Y M B O L S',
     &14X,'PHASE VOL.',5X,'SPIN FAC. ',10X,' WITH P  ',
     &10X,' WITH N  ')
       DO 501 I=1,NCHAN
            WCP(1,I)=WCP(1,I)/SUMP
            WCP(2,I)=WCP(2,I)/SUMP
                WCP00=WCP(1,I)+WCP(2,I)
            WCN(I)=WCN(I)/SUMN
            if(verbose) WRITE(13,502)I,(SYC(ICHN(L,I)),L=1,6),PVOL(I),
     &       SNS(I),WCP00,WCN(I)
 502   FORMAT(4X,I4,2X,6(A4,1X),5X,E10.2,5X,F5.1,5X,2(10X,F6.5,5X))
 501   CONTINUE
C            -----  CROSS-SECTIONS FOR DIFFERENT N ----
                 DO 503 I=2,6
                    SUMTOP(1,I)=SUMTOP(1,I)/SUMP
                    SUMTOP(2,I)=SUMTOP(2,I)/SUMN
 503             CONTINUE
           if(verbose) WRITE(13,506)
 506  FORMAT(///'    TOPOLOGICAL PROBABILITIES WITH P AND WITH N'///)
                   DO 504 K=1,2
           if(verbose) WRITE(13,505)(M,SUMTOP(K,M),M=2,6)
 504               CONTINUE
                               END IF
C
C     -- CHANNELS WITH DEFINITE CHARGE: COME BACK TO SU(2)-ISOSPIN --
C     --              TO FILL THE ICG(6,6000,3)                    --
C
            DO 600 N=2,6
        NVR1=ISTART(N)
        NVR2=IFINISH(N)
              DO 601 NCVAR=NVR1,NVR2
                    CALL IENUL
                     DO 602 I60=1,6
                  IPT6(I60)=ICHN(I60,NCVAR)
 602                 CONTINUE
                 NCHAN=NCVAR
                     DO 603 IL=1,5
                     FILIN(IL)=0.
 603                 CONTINUE
                 CALL CHACOR(IPT6,N)
 601          CONTINUE
 600       CONTINUE

       RETURN
       END



                           SUBROUTINE MESUS
C  ---------  QUANTUM NUMBERS AND NAMES OF THE MESON'S AND -----------
C  ------------------  BARION'S ZOO IN PRINTING ----------------------
      CHARACTER*4 SYMB,SYC
      logical verbose
      COMMON /MES/SYMB(44),SMM(44),ISM(44),ISM3(44),IY(44),JSP(44),
     &IG(44),IP(44),IC(44),ICH(44),SYC(44)
      COMMON /VERB/ verbose
      SAVE /MES/, /VERB/

           if(verbose) WRITE(13,1)
 1    FORMAT(///8X,'N.B. I=5 MEANS I=1/2'/3X,72('*')/3X,
     &'!NUMBER!','SYMBOL!',' M,GEV!','   I   !','  I3  !','  Y   !',
     &'  J   !','  G   !','  P   !','  CH  !'/3X,72('*'))
      DO I=1,44
                   IF (SMM(I).GE.1.E-06) THEN
      if(verbose) WRITE(13,2)I,SYMB(I),SMM(I),ISM(I),ISM3(I),IY(I),
     & JSP(I),IG(I),IP(I),ICH(I)
 2    FORMAT(3X,'!',I3,3X,'! ',A4,' !',F6.3,'! ',8(I3,3X,'!'))
                   END IF
      END DO
           if(verbose) WRITE(13,4)
 4    FORMAT(///8X,'TRUNCATED TABLE (WITHOUT CHARGE)')
           if(verbose) WRITE(13,10)
 10   FORMAT(///8X,'N.B. I=5 MEANS I=1/2'/3X,78('*')/3X,
     &'!NUMBER!','SYMBOL!',' M,GEV!','   I   !','      !','  Y   !',
     &'  J   !','      !','      !'/3X,78('*'))
      DO I=1,44
                   IF (SMM(I).GE.1.E-06) THEN
      if(verbose) WRITE(13,20)I,SYC(I),SMM(I),ISM(I),IY(I),JSP(I)
 20   FORMAT(3X,'!',I3,3X,'! ',A4,' !',F6.3,'! ',I3,3X,'!',6X,'!',
     &2(I3,3X,'!'),3(6X,'!'))
                   END IF
      END DO
      RETURN
      END

                             BLOCK DATA THREE
C          ---- QUANTUM NUMBERS AND NAMES OF PARTICLES ----
      CHARACTER*4 SYMB,SYC
      COMMON /MES/SYMB(44),SMM(44),ISM(44),ISM3(44),IY(44),JSP(44),
     &IG(44),IP(44),IC(44),ICH(44),SYC(44)
      COMMON /WHO/NUMWHO(44),NUMIS(3,44)
      COMMON /CHGR/ICG(6,6000,3),NSUPL(3)
      COMMON /VLF/VALFAC(6)
      SAVE /MES/, /WHO/, /CHGR/, /VLF/
      DATA SYMB/' PI+',' PI-',' K+ ',' K- ',' K0 ',' AK0',' PI0',' ETA',
     &          'ETAP',' RO+',' RO-',' K*+',' K*-',' K*0','AK*0',' RO0',
     &          ' OMG',' PHI','    ',    17*'    ',
     &          '  P ','  N ','SIG-','SIG+','SIG0','KSI-','KSI0','LAMB'/
      DATA SMM/  0.138, 0.138, 0.496, 0.496, 0.496, 0.496, 0.138, 0.549,
     &           0.958, 0.770, 0.770, 0.892, 0.892, 0.892, 0.892, 0.770,
     &           0.783, 1.019, 0.000,    17*0.000,
     &           0.939, 0.939, 1.193, 1.193, 1.193, 1.318, 1.318, 1.116/
      DATA ISM/      1,     1,     5,     5,     5,     5,     1,     0,
     &               0,     1,     1,     5,     5,     5,     5,     1,
     &               0,     0,     0,     17*0,
     &               5,     5,     1,     1,     1,     5,     5,     0/
      DATA ISM3/     1,    -1,     5,    -5,    -5,     5,     0,     0,
     &               0,     1,    -1,     5,    -5,    -5,     5,     0,
     &               0,     0,     0,     17*0,
     &               5,    -5,    -1,     1,     0,    -5,     5,     0/
      DATA IY/       0,     0,     1,    -1,     1,    -1,     0,     0,
     &               0,     0,     0,     1,    -1,     1,    -1,     0,
     &               0,     0,     0,     17*0,
     &               1,     1,     0,     0,     0,    -1,    -1,     0/
      DATA JSP/     9*0,   9*1,    0,     17*0,
     &              8*5/
      DATA IG/      -1,    -1,     0,     0,     0,     0,    -1,     1,
     &               1,     1,     1,     0,     0,     0,     0,     1,
     &              -1,    -1,     0,     17*0,
     &              8*0/
      DATA IP/     18*-1,          0,     17*0,
     &              8*1/
      DATA IC/       1,     1,     0,     0,     0,     0,     1,     1,
     &               1,    -1,    -1,     0,     0,     0,     0,    -1,
     &              -1,    -1,     0,     17*0,
     &              8*0/
      DATA ICH/      1,    -1,     1,    -1,     0,     0,     0,     0,
     &               0,     1,    -1,     1,    -1,     0,     0,     0,
     &               0,     0,     0,     17*0,
     &               1,     0,    -1,     1,     0,    -1,     0,     0/
      DATA SYC/ ' PI ','    ',' K3 ',' K4 ','    ','    ','    ',' ETA',
     &          'ETAP',' RO ','    ',' K12',' K13','    ','    ','    ',
     &          ' OMG',' PHI','    ',     17*'    ',
     &          'NUCL','    ','SIGM','    ','    ',' KSI','    ','LAMB'/
      DATA NUMWHO/   3,     0,     2,     2,     0,     0,     0,     1,
     &               1,     3,     0,     2,     2,     0,     0,     0,
     &               1,     1,     1,    17*1,
     &               2,     0,     3,     0,     0,     2,     0,     1/
      DATA NUMIS/1,2,7,        0,0,0,        3,5,0,         4,6,0,
     &           0,0,0,        0,0,0,        0,0,0,         8,0,0,
     &           9,0,0,     10,11,16,        0,0,0,       12,14,0,
     &         13,15,0,        0,0,0,        0,0,0,         0,0,0,
     &          17,0,0,       18,0,0,       19,0,0,     17*0,17*0,17*0,
     &         37,38,0,        0,0,0,     39,40,41,         0,0,0,
     &           0,0,0,      42,43,0,        0,0,0,        44,0,0/
      DATA VALFAC/1.,2.,6.,24.,120.,720./
           END
                          BLOCK DATA FOUR
C          ------- TABLE FOR K**2 IN INITIAL STATE ----------
      COMMON /AK/AKPQ(6,8),AMAN,ICHN(6,1000),NCHAN,
     &IEYE(44,600,3),NEYE(3)
      SAVE /AK/
      DATA AKPQ/0.25000, 0.10000, 0.50000, 0.00000, 0.00000, 0.15000,
     &          0.00000, 0.30000, 0.16666, 0.16666, 0.16666, 0.20000,
     &          0.00000, 0.30000, 0.16666, 0.16666, 0.16666, 0.20000,
     &          0.00000, 0.45000, 0.25000, 0.00000, 0.25000, 0.05000,
     &          0.00000, 0.00000, 0.00000, 0.50000, 0.00000, 0.50000,
     &          0.00000, 0.45000, 0.25000, 0.00000, 0.25000, 0.05000,
     &          0.00000, 0.00000, 0.00000, 0.50000, 0.00000, 0.50000,
     &          0.00000, 0.00000, 0.00000, 0.50000, 0.00000, 0.50000/
      DATA AMAN/0.940/,NCHAN/0/
           END
                          BLOCK DATA FIVE
C             --------  PERMUTATION OF SPINS  ----------
      INTEGER BM2, BM3, BM4, BM5, BM6
      INTEGER BARIONM2, BARIONM3, BARIONM4, BARIONM5, BARIONM6
      COMMON /MULB/ BM2(2,4),BM3(3,8),BM4(4,16),BM5(5,32),BM6(6,64)
      COMMON /MULBAR/ BARIONM2(2,4),BARIONM3(3,12),BARIONM4(4,32),
     &                BARIONM5(5,80),BARIONM6(6,192)
      SAVE /MULB/, /MULBAR/
      DATA BM2/0,0, 0,1, 1,0, 1,1/
      DATA BM3/0,0,0, 0,0,1, 0,1,0, 0,1,1,
     &         1,0,0, 1,0,1, 1,1,0, 1,1,1/
      DATA BM4/0,0,0,0, 0,0,0,1, 0,0,1,0, 0,0,1,1,
     &         0,1,0,0, 0,1,0,1, 0,1,1,0, 0,1,1,1,
     &         1,0,0,0, 1,0,0,1, 1,0,1,0, 1,0,1,1,
     &         1,1,0,0, 1,1,0,1, 1,1,1,0, 1,1,1,1/
      DATA BM5/0,0,0,0,0, 0,0,0,0,1, 0,0,0,1,0, 0,0,0,1,1,
     &         0,0,1,0,0, 0,0,1,0,1, 0,0,1,1,0, 0,0,1,1,1,
     &         0,1,0,0,0, 0,1,0,0,1, 0,1,0,1,0, 0,1,0,1,1,
     &         0,1,1,0,0, 0,1,1,0,1, 0,1,1,1,0, 0,1,1,1,1,
     &         1,0,0,0,0, 1,0,0,0,1, 1,0,0,1,0, 1,0,0,1,1,
     &         1,0,1,0,0, 1,0,1,0,1, 1,0,1,1,0, 1,0,1,1,1,
     &         1,1,0,0,0, 1,1,0,0,1, 1,1,0,1,0, 1,1,0,1,1,
     &         1,1,1,0,0, 1,1,1,0,1, 1,1,1,1,0, 1,1,1,1,1/
      DATA BM6/0,0,0,0,0,0, 0,0,0,0,0,1, 0,0,0,0,1,0, 0,0,0,0,1,1,
     &         0,0,0,1,0,0, 0,0,0,1,0,1, 0,0,0,1,1,0, 0,0,0,1,1,1,
     &         0,0,1,0,0,0, 0,0,1,0,0,1, 0,0,1,0,1,0, 0,0,1,0,1,1,
     &         0,0,1,1,0,0, 0,0,1,1,0,1, 0,0,1,1,1,0, 0,0,1,1,1,1,
     &         0,1,0,0,0,0, 0,1,0,0,0,1, 0,1,0,0,1,0, 0,1,0,0,1,1,
     &         0,1,0,1,0,0, 0,1,0,1,0,1, 0,1,0,1,1,0, 0,1,0,1,1,1,
     &         0,1,1,0,0,0, 0,1,1,0,0,1, 0,1,1,0,1,0, 0,1,1,0,1,1,
     &         0,1,1,1,0,0, 0,1,1,1,0,1, 0,1,1,1,1,0, 0,1,1,1,1,1,
     &         1,0,0,0,0,0, 1,0,0,0,0,1, 1,0,0,0,1,0, 1,0,0,0,1,1,
     &         1,0,0,1,0,0, 1,0,0,1,0,1, 1,0,0,1,1,0, 1,0,0,1,1,1,
     &         1,0,1,0,0,0, 1,0,1,0,0,1, 1,0,1,0,1,0, 1,0,1,0,1,1,
     &         1,0,1,1,0,0, 1,0,1,1,0,1, 1,0,1,1,1,0, 1,0,1,1,1,1,
     &         1,1,0,0,0,0, 1,1,0,0,0,1, 1,1,0,0,1,0, 1,1,0,0,1,1,
     &         1,1,0,1,0,0, 1,1,0,1,0,1, 1,1,0,1,1,0, 1,1,0,1,1,1,
     &         1,1,1,0,0,0, 1,1,1,0,0,1, 1,1,1,0,1,0, 1,1,1,0,1,1,
     &         1,1,1,1,0,0, 1,1,1,1,0,1, 1,1,1,1,1,0, 1,1,1,1,1,1/

      DATA BARIONM2/5,0, 5,1, 0,5, 1,5/
      DATA BARIONM3/5,0,0, 5,0,1, 5,1,0, 5,1,1,
     &              0,5,0, 0,5,1, 1,5,0, 1,5,1,
     &              0,0,5, 0,1,5, 1,0,5, 1,1,5/
      DATA BARIONM4/5,0,0,0, 5,0,0,1, 5,0,1,0, 5,0,1,1,
     &              5,1,0,0, 5,1,0,1, 5,1,1,0, 5,1,1,1,
     &              0,5,0,0, 0,5,0,1, 0,5,1,0, 0,5,1,1,
     &              1,5,0,0, 1,5,0,1, 1,5,1,0, 1,5,1,1,
     &              0,0,5,0, 0,0,5,1, 0,1,5,0, 0,1,5,1,
     &              1,0,5,0, 1,0,5,1, 1,1,5,0, 1,1,5,1,
     &              0,0,0,5, 0,0,1,5, 0,1,0,5, 0,1,1,5,
     &              1,0,0,5, 1,0,1,5, 1,1,0,5, 1,1,1,5/
      DATA BARIONM5/5,0,0,0,0, 5,0,0,0,1, 5,0,0,1,0, 5,0,0,1,1,
     &              5,0,1,0,0, 5,0,1,0,1, 5,0,1,1,0, 5,0,1,1,1,
     &              5,1,0,0,0, 5,1,0,0,1, 5,1,0,1,0, 5,1,0,1,1,
     &              5,1,1,0,0, 5,1,1,0,1, 5,1,1,1,0, 5,1,1,1,1,
     &              0,5,0,0,0, 0,5,0,0,1, 0,5,0,1,0, 0,5,0,1,1,
     &              0,5,1,0,0, 0,5,1,0,1, 0,5,1,1,0, 0,5,1,1,1,
     &              1,5,0,0,0, 1,5,0,0,1, 1,5,0,1,0, 1,5,0,1,1,
     &              1,5,1,0,0, 1,5,1,0,1, 1,5,1,1,0, 1,5,1,1,1,
     &              0,0,5,0,0, 0,0,5,0,1, 0,0,5,1,0, 0,0,5,1,1,
     &              0,1,5,0,0, 0,1,5,0,1, 0,1,5,1,0, 0,1,5,1,1,
     &              1,0,5,0,0, 1,0,5,0,1, 1,0,5,1,0, 1,0,5,1,1,
     &              1,1,5,0,0, 1,1,5,0,1, 1,1,5,1,0, 1,1,5,1,1,
     &              0,0,0,5,0, 0,0,5,0,1, 0,0,5,1,0, 0,0,5,1,1,
     &              0,1,0,5,0, 0,1,5,0,1, 0,1,5,1,0, 0,1,5,1,1,
     &              1,0,0,5,0, 1,0,5,0,1, 1,0,5,1,0, 1,0,5,1,1,
     &              1,1,0,5,0, 1,1,5,0,1, 1,1,5,1,0, 1,1,5,1,1,
     &              0,0,0,0,5, 0,0,0,1,5, 0,0,1,0,5, 0,0,1,1,5,
     &              0,1,0,0,5, 0,1,0,1,5, 0,1,1,0,5, 0,1,1,1,5,
     &              1,0,0,0,5, 1,0,0,1,5, 1,0,1,0,5, 1,0,1,1,5,
     &              1,1,0,0,5, 1,1,0,1,5, 1,1,1,0,5, 1,1,1,1,5/
      DATA BARIONM6/5,0,0,0,0,0, 5,0,0,0,0,1, 5,0,0,0,1,0, 5,0,0,0,1,1,
     &              5,0,0,1,0,0, 5,0,0,1,0,1, 5,0,0,1,1,0, 5,0,0,1,1,1,
     &              5,0,1,0,0,0, 5,0,1,0,0,1, 5,0,1,0,1,0, 5,0,1,0,1,1,
     &              5,0,1,1,0,0, 5,0,1,1,0,1, 5,0,1,1,1,0, 5,0,1,1,1,1,
     &              5,1,0,0,0,0, 5,1,0,0,0,1, 5,1,0,0,1,0, 5,1,0,0,1,1,
     &              5,1,0,1,0,0, 5,1,0,1,0,1, 5,1,0,1,1,0, 5,1,0,1,1,1,
     &              5,1,1,0,0,0, 5,1,1,0,0,1, 5,1,1,0,1,0, 5,1,1,0,1,1,
     &              5,1,1,1,0,0, 5,1,1,1,0,1, 5,1,1,1,1,0, 5,1,1,1,1,1,
     &              0,5,0,0,0,0, 0,5,0,0,0,1, 0,5,0,0,1,0, 0,5,0,0,1,1,
     &              0,5,0,1,0,0, 0,5,0,1,0,1, 0,5,0,1,1,0, 0,5,0,1,1,1,
     &              0,5,1,0,0,0, 0,5,1,0,0,1, 0,5,1,0,1,0, 0,5,1,0,1,1,
     &              0,5,1,1,0,0, 0,5,1,1,0,1, 0,5,1,1,1,0, 0,5,1,1,1,1,
     &              1,5,0,0,0,0, 1,5,0,0,0,1, 1,5,0,0,1,0, 1,5,0,0,1,1,
     &              1,5,0,1,0,0, 1,5,0,1,0,1, 1,5,0,1,1,0, 1,5,0,1,1,1,
     &              1,5,1,0,0,0, 1,5,1,0,0,1, 1,5,1,0,1,0, 1,5,1,0,1,1,
     &              1,5,1,1,0,0, 1,5,1,1,0,1, 1,5,1,1,1,0, 1,5,1,1,1,1,
     &              0,0,5,0,0,0, 0,0,5,0,0,1, 0,0,5,0,1,0, 0,0,5,0,1,1,
     &              0,0,5,1,0,0, 0,0,5,1,0,1, 0,0,5,1,1,0, 0,0,5,1,1,1,
     &              0,1,5,0,0,0, 0,1,5,0,0,1, 0,1,5,0,1,0, 0,1,5,0,1,1,
     &              0,1,5,1,0,0, 0,1,5,1,0,1, 0,1,5,1,1,0, 0,1,5,1,1,1,
     &              1,0,5,0,0,0, 1,0,5,0,0,1, 1,0,5,0,1,0, 1,0,5,0,1,1,
     &              1,0,5,1,0,0, 1,0,5,1,0,1, 1,0,5,1,1,0, 1,0,5,1,1,1,
     &              1,1,5,0,0,0, 1,1,5,0,0,1, 1,1,5,0,1,0, 1,1,5,0,1,1,
     &              1,1,5,1,0,0, 1,1,5,1,0,1, 1,1,5,1,1,0, 1,1,5,1,1,1,
     &              0,0,0,5,0,0, 0,0,0,5,0,1, 0,0,0,5,1,0, 0,0,0,5,1,1,
     &              0,0,1,5,0,0, 0,0,1,5,0,1, 0,0,1,5,1,0, 0,0,1,5,1,1,
     &              0,1,0,5,0,0, 0,1,0,5,0,1, 0,1,0,5,1,0, 0,1,0,5,1,1,
     &              0,1,1,5,0,0, 0,1,1,5,0,1, 0,1,1,5,1,0, 0,1,1,5,1,1,
     &              1,0,0,5,0,0, 1,0,0,5,0,1, 1,0,0,5,1,0, 1,0,0,5,1,1,
     &              1,0,1,5,0,0, 1,0,1,5,0,1, 1,0,1,5,1,0, 1,0,1,5,1,1,
     &              1,1,0,5,0,0, 1,1,0,5,0,1, 1,1,0,5,1,0, 1,1,0,5,1,1,
     &              1,1,1,5,0,0, 1,1,1,5,0,1, 1,1,1,5,1,0, 1,1,1,5,1,1,
     &              0,0,0,0,5,0, 0,0,0,0,5,1, 0,0,0,1,5,0, 0,0,0,1,5,1,
     &              0,0,1,0,5,0, 0,0,1,0,5,1, 0,0,1,1,5,0, 0,0,1,1,5,1,
     &              0,1,0,0,5,0, 0,1,0,0,5,1, 0,1,0,1,5,0, 0,1,0,1,5,1,
     &              0,1,1,0,5,0, 0,1,1,0,5,1, 0,1,1,1,5,0, 0,1,1,1,5,1,
     &              1,0,0,0,5,0, 1,0,0,0,5,1, 1,0,0,1,5,0, 1,0,0,1,5,1,
     &              1,0,1,0,5,0, 1,0,1,0,5,1, 1,0,1,1,5,0, 1,0,1,1,5,1,
     &              1,1,0,0,5,0, 1,1,0,0,5,1, 1,1,0,1,5,0, 1,1,0,1,5,1,
     &              1,1,1,0,5,0, 1,1,1,0,5,1, 1,1,1,1,5,0, 1,1,1,1,5,1,
     &              0,0,0,0,0,5, 0,0,0,0,1,5, 0,0,0,1,0,5, 0,0,0,1,1,5,
     &              0,0,1,0,0,5, 0,0,1,0,1,5, 0,0,1,1,0,5, 0,0,1,1,1,5,
     &              0,1,0,0,0,5, 0,1,0,0,1,5, 0,1,0,1,0,5, 0,1,0,1,1,5,
     &              0,1,1,0,0,5, 0,1,1,0,1,5, 0,1,1,1,0,5, 0,1,1,1,1,5,
     &              1,0,0,0,0,5, 1,0,0,0,1,5, 1,0,0,1,0,5, 1,0,0,1,1,5,
     &              1,0,1,0,0,5, 1,0,1,0,1,5, 1,0,1,1,0,5, 1,0,1,1,1,5,
     &              1,1,0,0,0,5, 1,1,0,0,1,5, 1,1,0,1,0,5, 1,1,0,1,1,5,
     &              1,1,1,0,0,5, 1,1,1,0,1,5, 1,1,1,1,0,5, 1,1,1,1,1,5/

                     END


                   SUBROUTINE ASG(ISMT,IYT,JSPT,NUMBP)
C        --------   TO FIND THE PARTICLE NUMBER  ----------
      CHARACTER*4 SYMB,SYC
      logical verbose
      COMMON /MES/SYMB(44),SMM(44),ISM(44),ISM3(44),IY(44),JSP(44),
     &IG(44),IP(44),IC(44),ICH(44),SYC(44)
      COMMON /VERB/ verbose
      SAVE /MES/, /VERB/
      DO 1 I=1,44
        IF(ISMT.NE.ISM(I)) GOTO 1
        IF (IYT.NE.IY(I))  GOTO 1
        IF(JSPT.NE.JSP(I)) GOTO 1
                  NUMBP=I
          GOTO 2
 1    CONTINUE
             if(verbose) WRITE(13,3)ISMT,IYT,JSPT
 3    FORMAT(5X,'MESSAGE FROM SUBROUTINE ASG: IMPOSSIBLE COMBINATION ',
     &'ISMT=',I4,'  IYT=',I4,'  JSPT=',I4)
                  NUMBP=19
 2    RETURN
      END

                         SUBROUTINE IENUL
C --------------------- FILL IEYE AND NEYE BY ZEROS -------------------
      COMMON /AK/AKPQ(6,8),AMAN,ICHN(6,1000),NCHAN,
     &IEYE(44,600,3),NEYE(3)
      SAVE /AK/
          DO K=1,3
              NEYE(K)=0
              DO M=1,600
                 DO L=1,44
                     IEYE(L,M,K)=0
                 END DO
              END DO
          END DO
      RETURN
      END




                      SUBROUTINE GIND(INDPAT,N,GN)
C  ---------  IDENTITY FACTOR CALCULATION AND CANCELATION OF THE  -----
C  ---------         CHANNELS, WHICH WERE TREATED BEFOR           -----
      DIMENSION INDPAT(8),ISEARCH(44),IFACT(44)
      COMMON /AK/AKPQ(6,8),AMAN,ICHN(6,1000),NCHAN,
     &IEYE(44,600,3),NEYE(3)
      SAVE /AK/
             DO 1 I=1,44
             ISEARCH(I)=0
             IFACT(I)=1
 1           CONTINUE
      DO 2 K=1,N
      ISEARCH(INDPAT(K))=ISEARCH(INDPAT(K))+1
      IFACT(INDPAT(K))=IFACT(INDPAT(K))*ISEARCH(INDPAT(K))
 2    CONTINUE
                IF(NCHAN.EQ.0) GOTO 5
      DO 3 IN=1,NEYE(1)
      DO 4 I=1,44
      IF(IEYE(I,IN,1).NE.ISEARCH(I)) GOTO 3
 4    CONTINUE
                    GN=-1.
      GOTO 10
 3    CONTINUE
 5                  GN=1.
           NEYE(1)=NEYE(1)+1
      DO 6 K=1,44
      IEYE(K,NEYE(1),1)=ISEARCH(K)
      GN=GN*IFACT(K)
 6    CONTINUE
 10   RETURN
      END

                   SUBROUTINE CHACOR(IPT6,N)
C         ---- TO DEFINE CHARGE COMBINATION AND ITS WEIGHT ----
      DIMENSION IPT6(6),IWORK(6),IWORKI(6),IWORKI3(6),IFIND(44),
     &IFACT(44)
      LOGICAL TWO_NUCLEON, FITTING
      CHARACTER*4 SYMB,SYC
      COMMON /WEITHS/WCP(4,1000),WCN(1000),WPG(4,6000),WNG(6000)
      COMMON /CHGR/ICG(6,6000,3),NSUPL(3)
      COMMON /MES/SYMB(44),SMM(44),ISM(44),ISM3(44),IY(44),JSP(44),
     &IG(44),IP(44),IC(44),ICH(44),SYC(44)
      COMMON /AK/AKPQ(6,8),AMAN,ICHN(6,1000),NCHAN,
     &IEYE(44,600,3),NEYE(3)
      COMMON /WHO/NUMWHO(44),NUMIS(3,44)
      COMMON /VLF/VALFAC(6)
      COMMON /FILIN/FILIN(5),NII(3)
      COMMON /SWITCH/TWO_NUCLEON, FITTING
      SAVE /WEITHS/, /CHGR/, /MES/, /AK/, /WHO/
      SAVE /VLF/, /FILIN/, /SWITCH/
           DO 100 M=1,3
           NII(M)=NSUPL(M)+1
 100       CONTINUE
            MX1=NUMWHO(IPT6(1))
      DO 1 I1=1,MX1
       IWORK(1)=NUMIS(I1,IPT6(1))
       MX2=NUMWHO(IPT6(2))
        DO 2 I2=1,MX2
         IWORK(2)=NUMIS(I2,IPT6(2))
         MX3=NUMWHO(IPT6(3))
          DO 3 I3=1,MX3
           IWORK(3)=NUMIS(I3,IPT6(3))
           MX4=NUMWHO(IPT6(4))
            DO 4 I4=1,MX4
             IWORK(4)=NUMIS(I4,IPT6(4))
             MX5=NUMWHO(IPT6(5))
              DO 5 I5=1,MX5
               IWORK(5)=NUMIS(I5,IPT6(5))
               MX6=NUMWHO(IPT6(6))
                DO 6 I6=1,MX6
                 IWORK(6)=NUMIS(I6,IPT6(6))
                 ICHGSUM=0
         DO I60=1,N
          ICHGSUM=ICHGSUM+ICH(IWORK(I60))
          IWORKI(I60)=ISM(IWORK(I60))
          IWORKI3(I60)=ISM3(IWORK(I60))
         END DO

       IF (  ( (IABS(ICHGSUM).LE.1).AND.TWO_NUCLEON )
     &                      .OR.
     &(((ICHGSUM.EQ.-1).OR.(ICHGSUM.EQ.0.)).AND.(.NOT.TWO_NUCLEON)) )
     &                      THEN
C
C ------------- TO AVOID THE FAMILIAR COMBINATIONS ------------
C
         DO IFD=1,44
           IFIND(IFD)=0
           IFACT(IFD)=1
         END DO

            DO K=1,N
               IFIND(IWORK(K))=IFIND(IWORK(K))+1
               IFACT(IWORK(K))=IFACT(IWORK(K))*IFIND(IWORK(K))
            END DO

                   IF (TWO_NUCLEON) THEN

          IF ( (ICHGSUM.EQ.1).OR.(ICHGSUM.EQ.0) ) THEN
C
C ------------     CASE OF ANNIHILATION WITH PP OR PN --------------
C
                  IPLACE=2-ICHGSUM

                  IF(NEYE(IPLACE).EQ.0) GOTO 68
                  DO 163 IN=1,NEYE(IPLACE)
                  DO I=1,44
                    IF(IEYE(I,IN,IPLACE).NE.IFIND(I)) GOTO 163
                  END DO
                          GOTO 6
 163              CONTINUE
 68       GNP=1.
          NEYE(IPLACE)=NEYE(IPLACE)+1
          NSUPL(IPLACE)=NSUPL(IPLACE)+1
               DO K=1,N
                ICG(K,NSUPL(IPLACE),IPLACE)=IWORK(K)
               END DO
            DO K=1,44
             IEYE(K,NEYE(IPLACE),IPLACE)=IFIND(K)
             GNP=GNP*IFACT(K)
            END DO
                  DO J=1,2
                   ISST=J+2*(IPLACE-1)
C Supply the remainder of two integer parameters, arg1/arg2
                   C=1.5-MOD(ISST,2)
                   Z=FLOAT(ICHGSUM)-0.5
               CALL SUBCLBH(N,IWORKI,IWORKI3,C,Z,UISO)
                   WPG(ISST,NSUPL(IPLACE))=UISO*VALFAC(N)/GNP
                   FILIN(ISST)=FILIN(ISST)+WPG(ISST,NSUPL(IPLACE))
                  END DO

          ELSE
C
C ------------      CASE OF ANNIHILATION WITH NN   --------------
C
          IF(NEYE(3).EQ.0) GOTO 78
              DO 73 IN=1,NEYE(3)
              DO I=1,44
           IF(IEYE(I,IN,3).NE.IFIND(I)) GOTO 73
              END DO
         GOTO 6
 73           CONTINUE
 78      GNN=1.
         NEYE(3)=NEYE(3)+1
         NSUPL(3)=NSUPL(3)+1
                DO K=1,N
                 ICG(K,NSUPL(3),3)=IWORK(K)
                END DO
             DO K=1,44
              IEYE(K,NEYE(3),3)=IFIND(K)
              GNN=GNN*IFACT(K)
             END DO
           C=1.5
           Z=-1.5
               CALL SUBCLBH(N,IWORKI,IWORKI3,C,Z,UISO)
                   WNG(NSUPL(3))=UISO*VALFAC(N)/GNN
                   FILIN(5)=FILIN(5)+WNG(NSUPL(3))
          END IF
                   ELSE

          IF (ICHGSUM.EQ.0) THEN
C
C ------------     CASE OF ANNIHILATION WITH P    --------------
C
                  IF(NEYE(1).EQ.0) GOTO 168
                  DO 63 IN=1,NEYE(1)
                  DO I=1,44
                    IF(IEYE(I,IN,1).NE.IFIND(I)) GOTO 63
                  END DO
                          GOTO 6
 63               CONTINUE
 168      GNP=1.
          NEYE(1)=NEYE(1)+1
          NSUPL(1)=NSUPL(1)+1
               DO K=1,N
                ICG(K,NSUPL(1),1)=IWORK(K)
               END DO
            DO K=1,44
             IEYE(K,NEYE(1),1)=IFIND(K)
             GNP=GNP*IFACT(K)
            END DO
                  DO ISST=1,2
                   C=FLOAT(ISST-1)
                   Z=FLOAT(ICHGSUM)
               CALL SUBCLBH(N,IWORKI,IWORKI3,C,Z,UISO)
                   WPG(ISST,NSUPL(1))=UISO*VALFAC(N)/GNP
                   FILIN(ISST)=FILIN(ISST)+WPG(ISST,NSUPL(1))
                  END DO

          ELSE
C
C ------------      CASE OF ANNIHILATION WITH N    --------------
C
          IF(NEYE(3).EQ.0) GOTO 178
              DO 173 IN=1,NEYE(3)
              DO I=1,44
           IF(IEYE(I,IN,3).NE.IFIND(I)) GOTO 173
              END DO
         GOTO 6
 173           CONTINUE
 178     GNN=1.
         NEYE(3)=NEYE(3)+1
         NSUPL(3)=NSUPL(3)+1
                DO K=1,N
                 ICG(K,NSUPL(3),3)=IWORK(K)
                END DO
             DO K=1,44
              IEYE(K,NEYE(3),3)=IFIND(K)
              GNN=GNN*IFACT(K)
             END DO
           C=1.
           Z=-1.
               CALL SUBCLBH(N,IWORKI,IWORKI3,C,Z,UISO)
                   WNG(NSUPL(3))=UISO*VALFAC(N)/GNN
                   FILIN(5)=FILIN(5)+WNG(NSUPL(3))
          END IF
                   END IF
       END IF
 6         CONTINUE
 5          CONTINUE
 4           CONTINUE
 3            CONTINUE
 2             CONTINUE
 1              CONTINUE
C
C     ----  THE SUM OF CHARGE-CORRELATION MUST BE EQUAL TO ONE ----
C

            DO ISST=1,2
             IF(FILIN(ISST).GT.1.E-20) THEN
               DO K=NII(1),NSUPL(1)
         WPG(ISST,K)=100.*WCP(ISST,NCHAN)*WPG(ISST,K)/FILIN(ISST)
               END DO
             END IF
            END DO

           DO ISST=3,4
             IF(FILIN(ISST).GT.1.E-20) THEN
               DO K=NII(2),NSUPL(2)
         WPG(ISST,K)=100.*WCP(ISST,NCHAN)*WPG(ISST,K)/FILIN(ISST)
               END DO
             END IF
           END DO


                     IF(FILIN(5).GT.1.E-20) THEN
         DO K=NII(3),NSUPL(3)
         WNG(K)=100.*WCN(NCHAN)*WNG(K)/FILIN(5)
         END DO
                     END IF

             RETURN
             END



           SUBROUTINE SUBCLBH(NP,IW,IW3,C,Z,U)
              use ClebschGordan, only : clebschSquared
              DIMENSION IW(6),IW3(6),CW(6),CW3(6)
       DO 100 I=1,NP
          CW(I)=FLOAT(IW(I))
          CW3(I)=FLOAT(IW3(I))
         IF(IW(I).GT.4) CW(I)=0.1*CW(I)
         IF(IABS(IW3(I)).GT.4) CW3(I)=0.1*CW3(I)
 100   CONTINUE
       N=NP
       GOTO (1,2,3,4,5,6)N
 1     WRITE(13,200)
 200   FORMAT(5X,'ERROR REPORTED BY SUBCLBH: NP MUST BE GREATER THAN 1')
       RETURN
 2      U=ClebschSquared(CW(1),CW(2),C,CW3(1),CW3(2),Z)
       RETURN
 3      U=US2(CW(1),CW3(1),CW(2),CW3(2),CW(3),CW3(3),C,Z)
       RETURN
 4      CPM=ABS(CW(1)-CW(2))
        CPX=CW(1)+CW(2)
        CALL HALF(CPX,SPT)
        XP=CW3(1)+CW3(2)
        CPM=AMAX1(ABS(XP),CPM)
            UV=0.
            CP=CPM
 41       UV=UV+ClebschSquared(CW(1),CW(2),CP,CW3(1),CW3(2),XP)*
     &                    US2(CP,XP,CW(3),CW3(3),CW(4),CW3(4),C,Z)
          CP=CP+SPT
          DF=CP-CPX
          IF(DF.GT.1.E-04) GOTO 40
          GOTO 41
 40     U=UV
       RETURN
 5     XPP=CW3(1)+CW3(2)+CW3(3)
       CPPM=ABS(XPP)
        CPPX=AMIN1(CW(1)+CW(2)+CW(3),CW(4)+CW(5)+C)
        CALL HALF(CPPX,SPPT)
            UV=0.
            IF(CPPM.GT.(CPPX+0.01)) GOTO 50
            CPP=CPPM
 51       UV=UV+US2(CW(1),CW3(1),CW(2),CW3(2),CW(3),CW3(3),CPP,XPP)*
     &                    US2(CPP,XPP,CW(4),CW3(4),CW(5),CW3(5),C,Z)
          CPP=CPP+SPPT
          DF=CPP-CPPX
          IF(DF.GT.1.E-04) GOTO 50
          GOTO 51
 50     U=UV
       RETURN
 6      CPM=ABS(CW(1)-CW(2))
        CPX=CW(1)+CW(2)
        CALL HALF(CPX,SPT)
        XP=CW3(1)+CW3(2)
        CPM=AMAX1(ABS(XP),CPM)
        CP=CPM
               XPPP=CW3(1)+CW3(2)+CW3(3)+CW3(4)
               CPPPM=ABS(XPPP)
               CPPPX=AMIN1(CW(1)+CW(2)+CW(3)+CW(4),CW(5)+CW(6)+C)
               CALL HALF(CPPPX,SPPPT)
               UV=0.
               IF(CPPPM.GT.(CPPPX+0.01)) GOTO 60
 63            CPPP=CPPPM
 61       UV=UV+ClebschSquared(CW(1),CW(2),CP,CW3(1),CW3(2),XP)*
     &                  US2(CP,XP,CW(3),CW3(3),CW(4),CW3(4),CPPP,XPPP)*
     &                  US2(CPPP,XPPP,CW(5),CW3(5),CW(6),CW3(6),C,Z)
               CPPP=CPPP+SPPPT
               DF3=CPPP-CPPPX
               IF(DF3.GT.1.E-04)  GOTO 62
               GOTO 61
 62     CP=CP+SPT
               DF1=CP-CPX
               IF(DF1.GT.1.E-04)  GOTO 60
               GOTO 63

 60    U=UV
      RETURN
      END


                SUBROUTINE HALF(CPX,SPT)
C
C -------- TO DEFINE: SPT IS HALF-INTEGER IF CPX IS HALF-INTEGER  -------
C
                SPT=1.
                I2=IFIX(2.*CPX)
C Supply the remainder of two integer parameters, arg1/arg2
                I=MOD(I2,2)
                IF(I.EQ.1) SPT=0.5
                RETURN
                END



            FUNCTION US2(A1,X1,A2,X2,A3,X3,C,Z)
              use ClebschGordan, only : CG
               CPM=ABS(A1-A2)
               CPX=A1+A2
               CALL HALF(CPX,SPT)
               XP=X1+X2
               CPM=AMAX1(ABS(XP),CPM)
               UV=0.
               IF(CPM.GT.(CPX+0.01)) GOTO 30
               NLOOP=0
            CP=CPM
 31         UV=UV+(CG(nint(A1*2),nint(A2*2),nint(CP*2),
     &                nint(X1*2),nint(X2*2),nint(XP*2))*
     &             CG(nint(CP*2),nint(A3*2),nint(C*2),
     &                nint(XP*2),nint(X3*2),nint(Z*2)))**2
            CP=CP+SPT
            DF=CP-CPX
            NLOOP=NLOOP+1
            IF(NLOOP.GT.40) STOP
            IF(DF.GT.1.E-04) GOTO 30
            GOTO 31
 30   US2=UV
      RETURN
      END
