C***********************************************************************
C in this version, all occurences of PY... are replaced by RY...
C***********************************************************************

C*********************************************************************
C This is am modified version of the "original" FRITIOF v7.02
C Kai Gallmeister:
C * LUSTRF is taken out
C * Since this code has to be compiled with the automatic conversion
C   real -> double, all occurences of /LUJETS/ are replaced by /PYJETS/,
C   the corresponding common block of PYTHIAv6.2
C   Due to this, many routines do not need to know any more, whether
C   FRITIOF or PYTHIA generated the output
C
c===================================================================
c FRITIOF 7.02
c with (minor) changes by Kai Gallmeister [KG], 27.1.2003..27.1.2003
c===================================================================
c
c - There is no direct call to ULMASS in any FRITIOF part any more.
c   Instead everything is redirected to ULGAMASS. This ensures, that
c   for a parton-parton-event the initial partons have the same mass
c   during the whole event!
c   If events including nuclei are still working is not checked!
c   If this works correctly for hadronic events is not checked!
c
c
c===================================================================

c===================================================================
c...Start a Hadron-Hadron-Event
c...Use this instead of FREVENT(...)!
c
c KCD(1|2),NPROT(1|2), EXMA(1|2,1..2) have to be set up correctly.
c

      subroutine FRHADHAD(FRAME,BEAM,TARGET,WIN, BeamMass, Targetmass)

      CHARACTER*(*) FRAME,BEAM,TARGET

c internal common blocks:

      common /DATAULGA/ GAMASS(2)
      save /DATAULGA/

      GAMASS(1) = BeamMass
      GAMASS(2) = TargetMass

      call FREVENT(FRAME,BEAM,TARGET,WIN)

      end

c===================================================================
c...wrapper for function calls according ULMASS(IDN(L,N))
c
c works only for a hadron-hadron event, not for nuclei!
c
c The Masses in GAMASS have to be set externally before any
c calling of FREVENT.
c
c No checks, no warranties, all responsibilities to the user...
c

      function ULGAMASS(L,N)

      common /DATAULGA/ GAMASS(2)
      save /DATAULGA/

      PARAMETER (KSZ2=300)
      COMMON/FRINTN3/IDN(2,KSZ2),FMN(2,KSZ2),NUC(2,3000)
      save /FRINTN3/

c$$$      write(*,*) 'ULGAMASS:',L,N,ULMASS(IDN(L,N)),GAMASS(L)

      if (N.le.1) then
         ULGAMASS = GAMASS(L)
      else
         write(*,*) 'ULGAMASS: N>1:',N,' STOP!!!'
         stop
         ULGAMASS = ULMASS(IDN(L,N))
      endif
      return

      end

c===================================================================
c HERE STARTS THE 'ORIGINAL' FRITIOF PACKAGE...
c changes by Kai Gallmester (all are marked by [KG] !!!):
c - function calls to ULMASS() are redirected to ULGAMASS
c
c===================================================================

C*************************************************************************
C*                                                                       *
C*   FRITIOF 7.02 subroutine packages                                     *
C*                                                                       *
C------------------------------------------------------------------------*
C*************************************************************************

C**************************** FREDITD ***********************************

      SUBROUTINE FREDITD()

C...This is a dummy subroutine in connection to option KFR(13)>=4.
C...User may elect to write his own special purpose codes here that
C...edits and compresses the event record LUJETS.  There may also be
C...times the user wish to keep a trace on certain decay products by
C...assigning them a special codes here.

      PARAMETER (KSZJ=4000,KSZ1=20)
      COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)
      SAVE /FRPARA1/

      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      integer N,NPAD,K
      double precision P,V
      SAVE /PYJETS/

C.....DUMMY
      RETURN
      END

C**************************** END FREDITD *******************************


C**************************** FRINITA ************************************

      SUBROUTINE FRINITA(CFRAME,CBEAM,CTARG,WIN)

C Purpose: identifies the particles involved and fills common block
C Calculates initial momenta and fills common block FRINTN0.
C Write
C The program is stopped if the BEAM or TARGET particles or the frame
C are not recognized.
C This routine calls FRHILDN for setting the particle codes and masses.

      PARAMETER (KSZ1=20, KSZ2=300)
      CHARACTER CFRAME*4,CBEAM*4,CTARG*4, PARTIC*4,PACD*4
      CHARACTER INIT*42, CGDATE*11
      COMMON/FRINTN0/PLI0(2,4),AOP(KSZ1),IOP(KSZ1),NFR(KSZ1)
      COMMON/FRINTN3/IDN(2,KSZ2),FMN(2,KSZ2),NUC(2,3000)
      COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)
      COMMON/FRCODES/IPT(2),PACD(27),NNUC(27),NPROT(27),KCD(27)
     >           ,RO1(27,2),EXMA(9,2)
      COMMON/FRGEOMC/NFLG,NUMROP,NUMROT,NUMREP
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)

      SAVE VFR10
      SAVE /FRINTN0/,/FRINTN3/,/FRPARA1/,/FRCODES/,/FRGEOMC/,/LUDAT1/
      DATA VFR10 /0./

      DATA IQFST /0/
      DATA CGDATE /'14 NOV 1993'/

      IWR = 0
      IF(KFR(11).LT.0.OR.IQFST.LE.KFR(11)) IWR=1
C.....................Identify the particles involved...................
c     .

c.....IPT(i) gives for example 15 for Oxygen, 1 for NEW1 ....
c.....IOP parameters see H. Pi paper page 186
      DO 110 L=1, 2
         IPT(L) = 0
         IF(L.EQ.1) PARTIC=CBEAM
         IF(L.EQ.2) PARTIC=CTARG
         DO 100 J=1,27
            IF(PARTIC.EQ.PACD(J)) THEN
               IPT(L)=J
               IOP(3+2*(L-1))=NNUC(IPT(L))
               IOP(4+2*(L-1))=NPROT(IPT(L))
               IOP(6+L)=KCD(IPT(L))
               GOTO 110
            ENDIF
 100     CONTINUE
 110  CONTINUE

C................................................write headers........
      IF(MSTU(12).EQ.1) THEN
         WRITE(MSTU(11),1000) CGDATE
         CALL LULIST(0)
         WRITE(MSTU(11),*)
         WRITE(MSTU(11),*)
      ENDIF

c......ipt(1) or ipt(2) .eq.0 then no target or projectile is specified
      IF(IPT(1).EQ.0) WRITE(MSTU(11),1100) CBEAM
      IF(IPT(2).EQ.0) WRITE(MSTU(11),1200) CTARG
      IF(IPT(1).EQ.0.OR.IPT(2).EQ.0) STOP 99901

c......iop(3,5).eq.0 => no proj/target nucleons, then stop
      IF(IOP(3).EQ.0.OR.IOP(5).EQ.0) THEN
         WRITE(MSTU(11),1300)
         STOP 99902
      ENDIF

c......test for the number of nucleons in the target/projectile:
      IF(IOP(3).GT.KSZ2) WRITE(MSTU(11),1310) IOP(3)
      IF(IOP(5).GT.KSZ2) WRITE(MSTU(11),1320) IOP(5)
      IF(IOP(3).GT.KSZ2.OR.IOP(5).GT.KSZ2) STOP 99903

C.....set minimum mass and minimum diffractive mass

      DO 112 L=1, 2
         IF(NNUC(IPT(L)).LE.1.and.IPT(L).LE.9) THEN
            AOP(8+L) = EXMA(IPT(L),1)
c..jochen: habe diffrctive masse auf min mass gesetzt!!!!!!!!!
            AOP(10+L) = EXMA(IPT(L),1)
c            AOP(10+L) = EXMA(IPT(L),2)
         ELSEIF(NNUC(IPT(L)).GE.2) THEN
            AOP(8+L) = EXMA(8,1)
            AOP(10+L) = EXMA(8,2)
         ENDIF
         IF(KCD(IPT(L)).NE.0) THEN
c[KG]       PMAS = ULMASS(KCD(IPT(L)))
            PMAS = ULGAMASS(L,1) ! [KG]
            IF(AOP(8+L).LT.PMAS) AOP(8+L)=PMAS
         ENDIF
         IF(AOP(8+L).LE.0.0001) AOP(8+L) = EXMA(8,1)
         IF(AOP(10+L).LT.AOP(8+L)) AOP(10+L) = AOP(8+L)
 112  CONTINUE


c      write(*,*) '===',ULGAMASS(1,1),AOP(9),AOP(11),IPT(1)
c      write(*,*) '===',ULGAMASS(2,1),AOP(10),AOP(12),IPT(2)



      CALL FRHILDN

      DO 115 L=1,2
         DO 115 LO=1,2
 115  PLI0(L,LO) = 0.

      IF(CFRAME.EQ.'CMS') THEN
         INIT=CBEAM//'-'//CTARG//' COLLIDER'//' '
         IF(IWR.EQ.1) WRITE(MSTU(11),1400) INIT
         IF(IWR.EQ.1) WRITE(MSTU(11),1500) WIN

C...IN CASE OF NUCLEUS, HERE WE HAVE NEGLECTED THE MASS DIFFERECE
C...BETWEEN THE NUCLEONS.........................................
         S0=WIN**2
         FP = FMN(1,1)**2
         FT = FMN(2,1)**2
         PL2=S0/4-(FP+FT)/2+(FP-FT)**2/(4*S0)
         ESEN=FRSQR(PL2 + FP, 'ESENPL' )
         PLI0(1,4) =ESEN+FRSQR(PL2, 'euiron' )
         PLI0(1,3) =FP/PLI0(1,4)

         ESEN=FRSQR(PL2 + FT, 'ESENFT')
         PLI0(2,3) =ESEN+FRSQR(PL2,'iopji1' )
         PLI0(2,4) =FT/PLI0(2,3)
         AOP(1) = WIN

      ELSEIF(CFRAME.EQ.'FIXT') THEN

         PLAB = WIN
         ELAB = SQRT(PLAB**2+FMN(1,1)**2)
         PLI0(1,4) =ELAB+PLAB
         PLI0(1,3) =FMN(1,1)**2/PLI0(1,4)
         PLI0(2,4) =FMN(2,1)
         PLI0(2,3) =FMN(2,1)
         S0 = (PLI0(1,4)+PLI0(2,4))*(PLI0(1,3)+PLI0(2,3))
         AOP(1) = SQRT(S0)
         INIT=CBEAM//' ON '//CTARG//' FIXED TARGET'//' '
         IF(IWR.EQ.1) THEN
            WRITE(MSTU(11),1400) INIT
            IF(IOP(3).GT.1) WRITE(MSTU(11),1600) WIN, AOP(1)
            IF(IOP(3).EQ.1) WRITE(MSTU(11),1610) WIN, AOP(1)
         ENDIF
      ELSE
         WRITE(MSTU(11),2000) CFRAME
         STOP 99904
      ENDIF

      IF(IWR.EQ.1) THEN
         WRITE(MSTU(11),2005) IOP(3),IOP(4),IOP(5),IOP(6)
         IF(IOP(3)+IOP(5).GT.2) WRITE(MSTU(11),2007)
      ENDIF


C.....Evaluate cross sections :

      IF(NFR(1).GT.0.AND.ABS(VFR(10)-VFR10).LT.0.001) THEN
         VFR(10) = 0.
         VFR(11) = 0.
      ENDIF
      CALL FRQPROB(IDN(1,1),IDN(2,1),IWR)
      VFR10 = VFR(10)

C.....Set up a few control parameters for the geometry package:
C.....NFLG is the entry control of subroutines FRPACOL & FRAACOL

      NFLG=0
      NUMROP=1
      if (IOP(3).gt.25) NUMROP=3
      NUMROT=1
      if (IOP(5).gt.25) NUMROT=3
      if (KFR(6).eq.1.or.(KFR(6).eq.2.and.IOP(5).gt.79) ) NUMROT=200
      NUMREP=1
      if (IOP(3).gt. 50.or.IOP(5).gt. 50) NUMREP=3
      if (IOP(3).gt.100.or.IOP(5).gt.100) NUMREP=6
      if (IOP(3).gt.200.or.IOP(5).gt.200) NUMREP=10


      IF(NFR(1).EQ.0) THEN
         IF(KFR(2).EQ.0) THEN
            WRITE(MSTU(11),2100)
         ELSEIF(KFR(2).EQ.1) THEN
            WRITE(MSTU(11),2110)
            WRITE(MSTU(11),2120) VFR(8),VFR(9)
         ENDIF
         IF(KFR(1).EQ.0) THEN
            WRITE(MSTU(11),2130)
         ELSEIF(KFR(1).EQ.1) THEN
            WRITE(MSTU(11),2140)
         ENDIF
         IF(IOP(5).GT.1) THEN
            IF(KFR(3).EQ.0) THEN
               WRITE(MSTU(11),2200)
            ELSEIF(KFR(3).EQ.1.or.KFR(3).EQ.3) THEN
               WRITE(MSTU(11),2210)
            ELSEIF(KFR(3).EQ.2.or.KFR(3).EQ.3) THEN
               WRITE(MSTU(11),2220) VFR(1),VFR(2)
            ENDIF
            IF(KFR(4).EQ.0) THEN
               WRITE(MSTU(11),2230)
            ELSEIF(KFR(4).EQ.1) THEN
               WRITE(MSTU(11),2240)
            ENDIF
            IF(KFR(5).EQ.0) THEN
               WRITE(MSTU(11),2250)
            ELSEIF(KFR(5).EQ.1) THEN
               WRITE(MSTU(11),2254)
            ELSEIF(KFR(5).EQ.2) THEN
               WRITE(MSTU(11),2258)
            ENDIF
            IF(KFR(6).EQ.0.or.(KFR(6).eq.2.and.IOP(5).le.79)) THEN
               WRITE(MSTU(11),2260)
            ELSE
               WRITE(MSTU(11),2270) VFR(4), VFR(5)
            ENDIF
         ENDIF
         IF(KFR(7).eq.0) THEN
            WRITE(MSTU(11),2280)
         ELSEIF(KFR(7).eq.1) THEN
            WRITE(MSTU(11),2290)
         ELSEIF(KFR(7).eq.2) THEN
            WRITE(MSTU(11),3000)
         ENDIF
         WRITE(MSTU(11),*) '  '
      ENDIF

      IQFST= IQFST+1

C........................FORMATS FOR INITIALIZATION AND ERROR INFORMATION
 1000 FORMAT(//20X,'THE LUND MONTE CARLO - FRITIOF VERSION 7.02'/
     *         20X,'LAST DATE OF CHANGE/BUG FIXING: ', A11)
 1100 FORMAT(1X,'ERROR: UNRECOGNIZED BEAM PARTICLE ''',A,
     * '''. EXECUTION STOPPED.')
 1200 FORMAT(1X,'ERROR: UNRECOGNIZED TARGET PARTICLE ''',A,
     * '''. EXECUTION STOPPED.')
 1300 FORMAT(1X,'ERROR: PARTICLES NOT WELL DEFINED. EXECUTION STOPPED.')
 1310 FORMAT(1X,'ERROR: TOO LARGE PROJECTILE, IOP(3)= ',I5,
     * '. EXECUTION STOPPED.')
 1320 FORMAT(1X,'ERROR: TOO LARGE TARGET,IOP(5)= ',I5,
     * '. EXECUTION STOPPED.')
 1400 FORMAT(/1X,77('=')/1X,'|',75X,'|'/1X,'|',8X,'FRITIOF WILL BE ',
     *'INITIALIZED FOR',1X,A34,1X,'|')
 1500 FORMAT(1X,'|',15X,'AT',1X,F10.3,1X,'GEV CENTER-OF-MASS ENERGY',
     *21X,'|'/1X,'|',75X,'|'/1X,77('='))
 1600 FORMAT(1X,'|',12X,'AT',1X,F10.3,1X,'GEV/C LAB-MOMENTUM PER NUCLEON
     *',19X,'|'/1x,'|',18X,'Equivalent CMS energy W= ',F9.4,1X,'GeV',
     > 19X,'|'/1X,'|',75X,'|'/1X,77('='))
 1610 FORMAT(1X,'|',23X,'AT',1X,F10.3,1X,'GEV/C LAB-MOMENTUM
     *',8X,'|'/1x,'|',18X,'Equivalent CMS energy W= ',F9.4,1X,'GeV',
     > 19X,'|'/1X,'|',75X,'|'/1X,77('='))
 2000 FORMAT(1X,'ERROR: UNRECOGNIZED COORDINATE FRAME ''',A,
     *'''. EXECUTION STOPPED.')
 2005 FORMAT(1X,5X,'PROJECTILE (A,Z)= ','(',I3,',',I3,')',
     > 4X,'TARGET (A,Z)= ','(',I3,',',I3,')',/)
 2007 FORMAT(/,/,4X,'A REMINDER: if the event is listed by LULIST,',
     > 1x,'the lines without' /,4x,'character names represent the',
     > 1x,'spectator nuclei.')
2100   FORMAT(/,/,4X,'No gluon radiation')
2110   FORMAT(    4X,'Gluon radiation included')
2120   FORMAT(4X,'Mu (projectile) =',F5.2,1x,'GeV;',3x,
     >           'Mu (target) =',F5.2,1x,'GeV')
2130   FORMAT(    4X,'Fragmentation not performed')
2140   FORMAT(    4X,'Fragmentation performed')
2200   FORMAT(    4X,'All interactions recorded')
2210   FORMAT(    4X,'Spectator veto')
2220   FORMAT(4X,'Impact parameter restricted in',
     >                                 1x,f7.3,'-',f7.3,' fm')
2230   FORMAT(    4X,'NO Fermi motion')
2240   FORMAT(    4X,'Fermi motion included')
2250   FORMAT(    4X,'Overlap function: eikonal')
2254   FORMAT(    4X,'Overlap function: gaussian')
2258   FORMAT(    4X,'Overlap function: gray disc')
2260   FORMAT(    4X,'No target nucleus deformation')
2270   FORMAT(    4X,'Target nucleus deformation applied.',1x,
     >            'Dipole and quadrupole coeff: ',f7.3,',',2x,f7.3)
2280   FORMAT(    4X,'No QCD parton scattering')
2290   FORMAT(    4X,'Hardest Rutherford parton scattering included')
3000   FORMAT(    4X,'Multiple parton scattering included')

      RETURN
      END

C******************************** END FRINITA ***************************


C******************************** FRHILDN *******************************

      SUBROUTINE FRHILDN

C...This routine sets particle codes and masses. Fills common block
C...FRINTN3-IDN,FMN; randomly order the neutrons and protons.

      PARAMETER (KSZ1=20,KSZ2=300)
      COMMON/FRINTN0/PLI0(2,4),AOP(KSZ1),IOP(KSZ1),NFR(KSZ1)
      COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)
      COMMON/FRINTN3/IDN(2,KSZ2),FMN(2,KSZ2),NUC(2,3000)
      SAVE /FRINTN0/,/FRPARA1/,/FRINTN3/

      DO 100 L=1,2

         IF(IOP(3+2*(L-1)).LE.1) THEN
c..........................iop(3,5)=number of projectile,target nucleons
            IDN(L,1)=IOP(6+L)
c..........................iop(7,8)=KF-Code of projectile,target
c[KG]:      FMN(L,1)=ULMASS(IDN(L,1))
            FMN(L,1)=ULGAMASS(L,1) ! [KG]
c..........................FMN gives the corresponding mass
         ELSE

            IPR=0
            INU=0
            IZ0=IOP(4+2*(L-1))
            IA0=IOP(3+2*(L-1))
            DO 30 I=1, IA0
               S=PYR(0)
               Q=FLOAT(IZ0-IPR)/(IA0-IPR-INU)
               IF (IPR.LT.IZ0.AND.(S.LT.Q.OR.INU.EQ.IA0-IZ0)) THEN
                  IDN(L,I)=2212
                  IPR= IPR+1
               ELSE
                  IDN(L,I)=2112
                  INU= INU+1
               ENDIF
c[KG]:         FMN(L,I)=ULMASS(IDN(L,I))
               FMN(L,I)=ULGAMASS(L,I) ! [KG]
 30         CONTINUE
            IF(IPR.NE.IZ0.OR.INU.NE.IA0-IZ0) CALL FRMGOUT(0,0,
     >           ' Proton or Neutron numbers incorrect!',float(IA0)
     $           ,float(IZ0),float(IPR),float(INU),0.)
         ENDIF
 100  CONTINUE

      RETURN
      END

C******************************** END FRHILDN ****************************

C******************************** FRINGEB  *******************************

      SUBROUTINE FRINGEB

C........................This routine administrates one complete event

      PARAMETER (KSZJ=4000,KSZ1=20,KSZ2=300)
      COMMON/FRINTN0/PLI0(2,4),AOP(KSZ1),IOP(KSZ1),NFR(KSZ1)
      COMMON/FRINTN1/PPS(2,KSZ2,5),PPH(2,KSZ2,5),PPSY(2,KSZ2,5),PPA(2,5)
      COMMON/FRINTN3/IDN(2,KSZ2),FMN(2,KSZ2),NUC(2,3000)
      COMMON/FRINTN4/KFEND(2,KSZ2,2)
      COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)
      COMMON/FRGEOMC/NFLG,NUMROP,NUMROT,NUMREP
      COMMON/FRCONT2/ICT(10),ICTT(10)

      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      integer N,NPAD,K
      double precision P,V
      SAVE /PYJETS/

      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)

      SAVE /FRINTN0/,/FRINTN1/,/FRINTN3/,/FRINTN4/,/FRPARA1/,/FRGEOMC/,
     >     /FRCONT2/,/LUDAT1/

      NFR(1) = NFR(1) + 1
  2   IOP(16)=0

C.....Randomly order protons and neutrons:
      CALL FRHILDN

C.....Create the nuclei and calculate number of collisions
      CALL FRANGAN
C.....Fix the end flavors for wounded strings
      DO L=1,2
        DO II=1, IOP(9+(L-1))
c............iop(9/19)=number of wounded projectile/target nucleons
          CALL FRBELEO(KFEND(L,II,1),KFEND(L,II,2),IDN(L,II))
        enddo
      enddo
C.....Generate the masses and momenta after the collisions
 10   CALL FRRINGO

      MSTU24=0
      NMEM=0
C.. Fill the strings and emit gluons:
      N=0
      NSTR = 0
      IOP(15) = 0
      DO L=1, 2
        DO J=1,IOP(8+L)
          CALL FRTORST(L,J)
          NSTR = NSTR+1
          NQG = N-NMEM

          IF( (KFR(1).EQ.1.AND.KFR(13).GE.1).AND.
     >         (NSTR.GT.100.OR.NQG.GT.(KSZJ-N)/10) ) THEN
            CALL  LUEXEC
c................administrates the fragmentation and decay chain
            IF(MSTU(24).EQ.4) MSTU24=1
c................default: mstu(24)=0
            IF(KFR(13).LE.3) THEN
              CALL LUEDIT(KFR(13))
            ELSEIF(KFR(13).GE.4) THEN
              CALL FREDITD()
            ENDIF
            NMEM=N
            NSTR=0
          ENDIF
        enddo
      enddo

C....TO ADD ONTO LUJETS THE COLOUR NEUTRAL PARTICLES THAT MAY
C....HAVE BEEN PRODUCED FROM PARTON-PARTON PROCESSES:

      CALL FRFILHW

      IF(N.GE.KSZJ-2) CALL FRMGOUT(0,1,
     > 'LUJETS array size KSZJ must be expanded',float(N),float(KSZJ),
     >  0.,0.,0.)

      IF(KFR(1).EQ.1.AND.N.GT.NMEM) THEN
        CALL  LUEXEC
c..........administrates the fragmentation and decay chain
        IF(MSTU(24).EQ.4) MSTU24=1
c..........default: mstu(24)=0
        IF(KFR(13).GE.1.AND.KFR(13).LE.3) THEN
c..........compresses data in lujets
          CALL LUEDIT(KFR(13))
        ELSEIF(KFR(13).GE.4) THEN
c.........user defined compresseion
          CALL FREDITD()
        ENDIF
      ENDIF

C... Regenerate event in case of 'infinite-loop error' in Jetset:
*      IF(MSTU24.GT.0) GOTO 2
      IF(MSTU24.GT.0) return

C....Record the number of N-N collisions:
      IOP(2) = IOP(2)-ICT(3)-ICT(8)-ICT(10)
      NFR(3) = NFR(3) + IOP(2)

       NFR(4) = NFR(4) + IOP(13)
       IF(IOP(13).GE.1) NFR(5) = NFR(5)+1

C-----RECORDING OF IMPACT PARAMETER AND COUNTING OF SPECTATOR PROTONS

      IOP(11)=IOP(4)
      IOP(12)=IOP(6)
      DO L=1,2
        DO J=1,IOP(9+L-1)
          IF(IDN(L,J).EQ.2212) IOP(10+L)=IOP(10+L)-1
        enddo
      enddo

      IF(IDN(1,1).NE.2212.AND.IDN(1,1).NE.2112) IOP(11)=0

C.....Add the nuclei spectators onto the event record:
      DO L=1,2
        IF(IOP(3+2*(L-1))-IOP(8+L).GE.1) THEN
          N = N+1
          P(N,1)=PPA(L,1)
          P(N,2)=PPA(L,2)
          P(N,3)=0.5*(PPA(L,4)-PPA(L,3))
          P(N,4)=0.5*(PPA(L,4)+PPA(L,3))
          P(N,5)=PPA(L,5)
          K(N,1) = 4
          K(N,2) = (10000+IOP(10+L))*(-1)**(L-1)
          K(N,3) = 0
          K(N,4) = 0
          K(N,5) = 0
        ENDIF
c300  CONTINUE
      enddo
C....Check Energy, momentum and charge conservation:
      IF(KFR(14).EQ.1) CALL FRCHKEP(1)

C....Dump out the event for inspection when error occurs in FRMGOUT:
      IF(IOP(16).GE.1) CALL LULIST(2)

      RETURN
      END

C******************************** END FRINGEB  ***************************


C*************************************************************************
C*************************************************************************
C                                                                        *
C      This is the routine package for nuclear geometry                  *
C      Input FRINTN0, output FRINTN3_NUC                                  *
C                                                                        *
C*************************************************************************


C********************************* FRANGAN ****************************

      SUBROUTINE FRANGAN

      PARAMETER (KSZ1=20,KSZ2=300)
      COMMON/FRGEOMC/NFLG,NUMROP,NUMROT,NUMREP
      COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)
      COMMON/FRINTN0/PLI0(2,4),AOP(KSZ1),IOP(KSZ1),NFR(KSZ1)
      COMMON/FRINTN3/IDN(2,KSZ2),FMN(2,KSZ2),NUC(2,3000)
      SAVE /FRGEOMC/,/FRPARA1/,/FRINTN0/,/FRINTN3/

c....iop(3,5)=number of projectile/target nucleons
      IF (IOP(3).LE.1.and.IOP(5).LE.1) THEN
         CALL FRPPCOL
      ELSE
         CALL FROVLAP

         IF(IOP(3).LE.1.or.IOP(5).LE.1) then
            CALL FRPACOL
         ELSE
            CALL FRAACOL
         ENDIF

      ENDIF

      RETURN
      END

C********************************* END FRANGAN ****************************

C********************************* FRPPCOL *****************************

          SUBROUTINE FRPPCOL

C --- this routine takes care of hadron-hadron collisions

      PARAMETER (KSZ1=20,KSZ2=300)
      COMMON/FRINTN0/PLI0(2,4),AOP(KSZ1),IOP(KSZ1),NFR(KSZ1)
      COMMON/FRINTN3/IDN(2,KSZ2),FMN(2,KSZ2),NUC(2,3000)
      SAVE /FRINTN0/,/FRINTN3/

      IOP(9)=1
      IOP(10)=1
      IOP(2)=1
      NUC(1,1)=1
      NUC(2,1)=1

      RETURN
      END

C********************************* END FRPPCOL *************************

C********************************* FRPACOL *****************************

          SUBROUTINE FRPACOL

C --- this routine deals with p-A collisions working in the rest frame
C      of the target center and taking z axis parallel to the projectile
C      incident direction

      PARAMETER (KSZ1=20,KSZ2=300)
      COMMON/FRGEOMC/NFLG,NUMROP,NUMROT,NUMREP
      COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)
      COMMON/FRINTN0/PLI0(2,4),AOP(KSZ1),IOP(KSZ1),NFR(KSZ1)
      COMMON/FRINTN3/IDN(2,KSZ2),FMN(2,KSZ2),NUC(2,3000)
      SAVE /FRGEOMC/,/FRPARA1/,/FRINTN0/,/FRINTN3/
      DIMENSION CORT(KSZ2,3)

50    IF (NFLG.NE.0) GOTO 100
C ==> first entry
      NFLG=1
      NFLG2=0

C --- initialization
C.....RMIN is the minimum distance required between two nucleons.

      RMIN=VFR(3)
      RMIN2=RMIN*RMIN
      CUTOFF=AOP(7)*3.
C --- parameters of the nucleon density distribution
      NA=IOP(5)
        IF (KFR(6).EQ.0.OR.(KFR(6).EQ.2.AND.IOP(5).LE.79) ) THEN
         CALL FRSEARC(2,FMMT,RMMT)
        ELSE
         CALL FRNUCDF(NA,A0,A2,A4,RMAX3)
        ENDIF

100   IF (NFLG2.NE.0) GOTO 120
C ==> second entry
      NFLG2=1
      NFLG3=0
      NROT=0

C --- determine the coordinates of target nucleons
      NA=IOP(5)
        IF (KFR(6).EQ.0.OR.(KFR(6).EQ.2.AND.IOP(5).LE.79) ) THEN
          CALL FRNUCOR(2,NA,RMIN2,FMMT,RMMT,CORT)
        ELSE
          CALL FRNUCOD(NA,RMIN2,A0,A2,A4,RMAX3,CORT)
        ENDIF

120   IF (NFLG3.NE.0) GOTO 150
C ==> third entry
      NFLG3=1
      NREP=0
      NROT=NROT+1

        IF (KFR(6).EQ.0.OR.(KFR(6).EQ.2.AND.IOP(5).LE.79) ) THEN
C --- rotate the target nucleus 90 degrees
      DO 130 I=1,IOP(5)
      W=CORT(I,1)
      CORT(I,1)=CORT(I,2)
      CORT(I,2)=CORT(I,3)
      CORT(I,3)=W
130   CONTINUE
          ELSE
C --- rotate target nucleus randomly along each axis with probabity 1/3
C     it is checked that this rotation gives an even solid angle dist.
      CDELTA=-1.+2.*PYR(0)
      SDELTA=FRSQR(MAX(0.,1.-CDELTA**2), 'SDELCA' )
      RA=PYR(0)
        IF (RA.GT.0.6666667) THEN
C --- rotate aronud z axis
      DO 132 I=1,IOP(5)
      X=CORT(I,1)*CDELTA-CORT(I,2)*SDELTA
      CORT(I,2)=CORT(I,1)*SDELTA+CORT(I,2)*CDELTA
      CORT(I,1)=X
132   CONTINUE
        ELSE IF (RA.GT.0.3333333) THEN
C --- rotate aronud x axis
      DO 134 I=1,IOP(5)
      Y=CORT(I,2)*CDELTA-CORT(I,3)*SDELTA
      CORT(I,3)=CORT(I,2)*SDELTA+CORT(I,3)*CDELTA
      CORT(I,2)=Y
134   CONTINUE
        ELSE
C --- rotate aronud y axis
      DO 136 I=1,IOP(5)
      Z=CORT(I,3)*CDELTA-CORT(I,1)*SDELTA
      CORT(I,1)=CORT(I,3)*SDELTA+CORT(I,1)*CDELTA
      CORT(I,3)=Z
136   CONTINUE
        ENDIF
          ENDIF

C --- find out the scope in X-Y plane of target nucleus
      XMAXT=CORT(1,1)
      XMINT=XMAXT
      YMAXT=CORT(1,2)
      YMINT=YMAXT
        DO 140 I=2,IOP(5)
      IF (CORT(I,1).GE.XMAXT) XMAXT=CORT(I,1)
      IF (CORT(I,1).LE.XMINT) XMINT=CORT(I,1)
      IF (CORT(I,2).GE.YMAXT) YMAXT=CORT(I,2)
      IF (CORT(I,2).LE.YMINT) YMINT=CORT(I,2)
140     CONTINUE

C --- target area in X-Y plane to be shooted
      XMAX=XMAXT+CUTOFF
      XMIN=XMINT-CUTOFF
      YMAX=YMAXT+CUTOFF
      YMIN=YMINT-CUTOFF

      IF (NROT.EQ.NUMROT) NFLG2=0

C ==> fourth entry
150   NREP=NREP+1

C --- sample impact, (XP,YP), of projectile
        IF (KFR(3).EQ.2.OR.KFR(3).EQ.3) THEN
      BPRO=FRSQR(PYR(0)*(VFR(2)*VFR(2)-VFR(1)*VFR(1))
     >                    +VFR(1)*VFR(1),'BPRO09')
      BPHI=6.2832*PYR(0)
      XP=BPRO*COS(BPHI)
      YP=BPRO*SIN(BPHI)
        ELSE
      XP=(XMAX-XMIN)*PYR(0)+XMIN
      YP=(YMAX-YMIN)*PYR(0)+YMIN
        ENDIF
      AOP(2)=FRSQR(XP**2+YP**2,'bipa22')
      IOP(2)=0
      IOP(9)=0
      IOP(10)=0

          DO 200 I=1,IOP(5)
C --- distance between the projectile proton and a target nucleon
      R2=(XP-CORT(I,1))**2+(YP-CORT(I,2))**2
C --- judge if a binary collision takes place
      PP = FRVOV(R2)
        IF (PYR(0).LT.PP) THEN
      IOP(2)=IOP(2)+1
c.............iop(2)=number of n-n-subcollisions
      IOP(10)=IOP(10)+1
c.............iop(10)=number of wounded targetnucleons
      IF (IOP(2).GT.3000) CALL FRMGOUT(0,0,'Array NUC() needs to be
     > expanded (3000 not enough)', 0.,0.,0.,0.,0.)
      NUC(1,IOP(2))=1
      NUC(2,IOP(2))=I
        ENDIF
200   CONTINUE

      IF (IOP(10).GT.0) IOP(9)=1

      IF (NREP.EQ.NUMREP) NFLG3=0

      IF ((KFR(3).EQ.1.OR.KFR(3).EQ.3).AND.IOP(9).LT.IOP(3)) GOTO 50
      IF (IOP(9).EQ.0) GOTO 50

C --- make order numbers of wounded target nucleons tightly
      DO 770 I=1,IOP(5)
c...............iop(5)=number of target nucleons
        II=999
        DO 750 J=1,IOP(2)
c.................iop(2)=number of nn-subcollisions
        IF (NUC(2,J).GE.I.AND.NUC(2,J).LT.II) II=NUC(2,J)
750     CONTINUE
        IF (II.EQ.I) GOTO 770
        IF (II.EQ.999) GOTO 780
        DO 760 J=1,IOP(2)
        IF (NUC(2,J).EQ.II) THEN
          IDN(2,I)= IDN(2,II)
          FMN(2,I)= FMN(2,II)
          NUC(2,J)=I
        ENDIF
760     CONTINUE
770   CONTINUE
780   CONTINUE

      RETURN
      END

C********************************* END FRPACOL ************************

C********************************* FRAACOL ******************************

       SUBROUTINE FRAACOL

C --- this subroutine deals with A-A collisions working in the rest
C       frame of the target center and taking Z axis parallel to the
C       projectile incident direction

      PARAMETER (KSZ1=20,KSZ2=300)
      COMMON/FRGEOMC/NFLG,NUMROP,NUMROT,NUMREP
      COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)
      COMMON/FRINTN0/PLI0(2,4),AOP(KSZ1),IOP(KSZ1),NFR(KSZ1)
      COMMON/FRINTN3/IDN(2,KSZ2),FMN(2,KSZ2),NUC(2,3000)
      SAVE /FRGEOMC/,/FRPARA1/,/FRINTN0/,/FRINTN3/
      DIMENSION CORP(KSZ2,3),CORT(KSZ2,3),MARKP(KSZ2),MARKT(KSZ2)

50    IF (NFLG.NE.0) GOTO 100
C ==> first entry
      NFLG=1
      NFLG2=0

C --- initialization
      RMIN=VFR(3)
      RMIN2=RMIN*RMIN
      CUTOFF=AOP(7)*3.
C --- parameters of nucleon density distribution
      NA=IOP(3)
      CALL FRSEARC(1,FMMP,RMMP)
      NA=IOP(5)
        IF (KFR(6).EQ.0.OR.(KFR(6).EQ.2.AND.IOP(5).LE.79) ) THEN
            CALL FRSEARC(2,FMMT,RMMT)
        ELSE
        CALL FRNUCDF(NA,A0,A2,A4,RMAX3)
        ENDIF

100   IF (NFLG2.NE.0) GOTO 200
C ==> second entry
      NFLG2=1
      NFLG3=0
      NROTT=0

C --- determine the coordinates of target nucleons
      NA=IOP(5)
        IF (KFR(6).EQ.0.OR.(KFR(6).EQ.2.AND.IOP(5).LE.79) ) THEN
      CALL FRNUCOR(2,NA,RMIN2,FMMT,RMMT,CORT)
        ELSE
      CALL FRNUCOD(NA,RMIN2,A0,A2,A4,RMAX3,CORT)
        ENDIF

200   IF (NFLG3.NE.0) GOTO 300
C ==> third entry
      NFLG3=1
      NFLG4=0
      NROTT=NROTT+1

        IF (KFR(6).EQ.0.OR.(KFR(6).EQ.2.AND.IOP(5).LE.79) ) THEN
C --- rotate the target nucleus
        DO 210 I=1,IOP(5)
      WW=CORT(I,1)
      CORT(I,1)=CORT(I,2)
      CORT(I,2)=CORT(I,3)
      CORT(I,3)=WW
210     CONTINUE
          ELSE
C --- rotate target nucleus randomly along an axis with probabity 1/3
C     it is checked that this rotation gives an even solid angle dist.
      CDELTA=-1.+2.*PYR(0)
      SDELTA=FRSQR(MAX(0.,1.-CDELTA**2), 'SDE9A6')
      RA=PYR(0)
        IF (RA.GT.0.6666667) THEN
C --- rotate aronud z axis
      DO 212 I=1,IOP(5)
      X=CORT(I,1)*CDELTA-CORT(I,2)*SDELTA
      CORT(I,2)=CORT(I,1)*SDELTA+CORT(I,2)*CDELTA
      CORT(I,1)=X
212   CONTINUE
        ELSE IF (RA.GT.0.3333333) THEN
C --- rotate aronud x axis
      DO 214 I=1,IOP(5)
      Y=CORT(I,2)*CDELTA-CORT(I,3)*SDELTA
      CORT(I,3)=CORT(I,2)*SDELTA+CORT(I,3)*CDELTA
      CORT(I,2)=Y
214   CONTINUE
        ELSE
C --- rotate aronud y axis
      DO 216 I=1,IOP(5)
      Z=CORT(I,3)*CDELTA-CORT(I,1)*SDELTA
      CORT(I,1)=CORT(I,3)*SDELTA+CORT(I,1)*CDELTA
      CORT(I,3)=Z
216   CONTINUE
        ENDIF
          ENDIF

C --- find out the scope in X-Y plane of sampled target nucleons
      XMAXT=CORT(1,1)
      XMINT=XMAXT
      YMAXT=CORT(1,2)
      YMINT=YMAXT
        DO 220 I=2,IOP(5)
      IF (CORT(I,1).GE.XMAXT) XMAXT=CORT(I,1)
      IF (CORT(I,1).LE.XMINT) XMINT=CORT(I,1)
      IF (CORT(I,2).GE.YMAXT) YMAXT=CORT(I,2)
      IF (CORT(I,2).LE.YMINT) YMINT=CORT(I,2)
220     CONTINUE

      NROTP=0
C --- determine the coordinates of projectile nucleons
C       with respect to the rest frame of the projectile center
C       (Z axes of two frames are assumed to be parallel each other)
      NA=IOP(3)
      CALL FRNUCOR(1,NA,RMIN2,FMMP,RMMP,CORP)

      IF (NROTT.EQ.NUMROT) NFLG2=0

300   IF (NFLG4.NE.0) GOTO 400
C ==> fourth entry
      NFLG4=1
      NREP=0
      NROTP=NROTP+1

C --- rotate the projectile nucleus 90 degrees
        DO 310 I=1,IOP(3)
      WW=CORP(I,1)
      CORP(I,1)=CORP(I,2)
      CORP(I,2)=CORP(I,3)
      CORP(I,3)=WW
310     CONTINUE

C --- find out the scope of sampled projectile nucleons in thw X-Y plane
C       of the projectile rest frame
      XMAXP=CORP(1,1)
      XMINP=XMAXP
      YMAXP=CORP(1,2)
      YMINP=XMAXT
        DO 320 I=2,IOP(3)
      IF (CORP(I,1).GE.XMAXP) XMAXP=CORP(I,1)
      IF (CORP(I,1).LE.XMINP) XMINP=CORP(I,1)
      IF (CORP(I,2).GE.YMAXP) YMAXP=CORP(I,2)
      IF (CORP(I,2).LE.YMINP) YMINP=CORP(I,2)
320     CONTINUE

C --- start to treat the nucleus-nucleus collision
C --- first determine the area of the projectile nucleus shooting
C      with respect to the rest frame of the target center
      XMAX=XMAXT-XMINP+CUTOFF
      XMIN=XMINT-XMAXP-CUTOFF
      YMAX=YMAXT-YMINP+CUTOFF
      YMIN=YMINT-YMAXP-CUTOFF

      IF (NROTP.EQ.NUMROP) NFLG3=0

C ==> fifth entry
400   NREP=NREP+1

C --- sample impact of projectile with respect to the target rest frame
        IF (KFR(3).EQ.2.OR.KFR(3).EQ.3) THEN
      BPRO=FRSQR(PYR(0)*(VFR(2)*VFR(2)-VFR(1)*VFR(1))+
     >                      VFR(1)*VFR(1),'BPRO12')
      BPHI=6.2832*PYR(0)
      XPRO=BPRO*COS(BPHI)
      YPRO=BPRO*SIN(BPHI)
        ELSE
      XPRO=(XMAX-XMIN)*PYR(0)+XMIN
      YPRO=(YMAX-YMIN)*PYR(0)+YMIN
        ENDIF
      AOP(2)=FRSQR(XPRO**2+YPRO**2,'bipa222')

      IOP(2)=0
      IOP(9)=0
      IOP(10)=0

      DO 410 I=1,IOP(3)
      MARKP(I)=0
410   CONTINUE
      DO 420 I=1,IOP(5)
      MARKT(I)=0
420   CONTINUE

C --- treat the collisions between two nucleons
          DO 600 NUP=1,IOP(3)
C --- coordinates of projectile nucleon, (XP,YP)
C       with respect to the rest frame of the target center
      XP=CORP(NUP,1)+XPRO
      YP=CORP(NUP,2)+YPRO
      IF (XP.GT.XMAXT+CUTOFF) GOTO 600
      IF (XP.LT.XMINT-CUTOFF) GOTO 600
      IF (YP.GT.YMAXT+CUTOFF) GOTO 600
      IF (YP.LT.YMINT-CUTOFF) GOTO 600

          DO 500 NUT=1,IOP(5)
C --- distance between a projectile nucleon and a target nucleon
      R2=(XP-CORT(NUT,1))**2+(YP-CORT(NUT,2))**2
C --- judge if a binary collision takes place

      PP = FRVOV(R2)

      IF (PYR(0).LT.PP) THEN
         IOP(2)=IOP(2)+1
        IF (IOP(2).GT.3000) CALL FRMGOUT(0,0,'Array NUCT needs to be
     >  expanded (3000 insufficient)',0.,0.,0.,0.,0.)
         NUC(1,IOP(2))=NUP
         NUC(2,IOP(2))=NUT
       IF (MARKP(NUP).EQ.0) THEN
         MARKP(NUP)=1
         IOP(9)=IOP(9)+1
       ENDIF
       IF (MARKT(NUT).EQ.0) THEN
         MARKT(NUT)=1
         IOP(10)=IOP(10)+1
       ENDIF
      ENDIF
500   CONTINUE
600   CONTINUE

      IF (NREP.EQ.NUMREP) NFLG4=0

      IF ((KFR(3).EQ.1.OR.KFR(3).EQ.3).AND.IOP(9).LT.IOP(3)) GOTO 50
      IF (IOP(9).EQ.0) GOTO 50

C --- Order and pack the wounded nucleons in the front:
      DO 700 L=1, 2
      DO 670 I=1,IOP(3+2*(L-1))
        II=999
        DO 650 J=1,IOP(2)
650     IF (NUC(L,J).GE.I.AND.NUC(L,J).LT.II) II=NUC(L,J)

        IF (II.EQ.I) GOTO 670
        IF (II.EQ.999) GOTO 700
        DO 660 J=1,IOP(2)
        IF (NUC(L,J).EQ.II) THEN
          IDN(L,I)= IDN(L,II)
          FMN(L,I)= FMN(L,II)
          NUC(L,J)=I
        ENDIF
660     CONTINUE
670   CONTINUE
700   CONTINUE


      RETURN
      END

C********************************* END FRAACOL ***************************


C********************************* FRNUCOR *******************************

      SUBROUTINE FRNUCOR(L,NA,RMIN2,FMM,RMM,COR)

C --- this subroutine determines nucleon coordinates inside a nucleus
C     and recenter the sampled nucleons with respect to the rest frame of
C     the nucleus center. L=1 for proj, L=2 for target.

      PARAMETER (KSZ2=300)
      DIMENSION COR(KSZ2,3),SUM(3)

      DO 150 J=1,NA
C --- sample a nucleon from the nucleus
C --- first, sample r
100   RR1=RMM*PYR(0)
      RR2=FMM*PYR(0)

      FR = FRROR(L,RR1)

        IF (RR2.LT.FR) THEN
      R=RR1
        ELSE
      GOTO 100
        ENDIF

C --- then sample COS(sita) & fai
        DO 140 NUM=1,10
      CTHITA=1.-2.*PYR(0)
      STHITA=FRSQR(MAX(0.,1.-CTHITA**2), 'SITFAI')
      FAI=6.2832*PYR(0)
      COR(J,1)=R*STHITA*COS(FAI)
      COR(J,2)=R*STHITA*SIN(FAI)
      COR(J,3)=R*CTHITA
        IF (J.EQ.1) GOTO 150
C --- check if there are two nucleons too close each other
        DO 130 J1=1,J-1
      DICX=COR(J,1)-COR(J1,1)
      DICY=COR(J,2)-COR(J1,2)
      DICZ=COR(J,3)-COR(J1,3)
      DIC2=DICX*DICX+DICY*DICY+DICZ*DICZ
C --- if two nucleons too close each other, sample THITA & FAI once more
      IF (DIC2.LT.RMIN2) GOTO 140
130   CONTINUE
      GOTO 150
140   CONTINUE
C --- if 10 times of repeated saplings don't help, then sample R again
      GOTO 100
150   CONTINUE

C --- recenter the sampled nucleons within a nucleus
      SUM(1)=0.
      SUM(2)=0.
      SUM(3)=0.
        DO 170 J=1,NA
      SUM(1)=SUM(1)+COR(J,1)
      SUM(2)=SUM(2)+COR(J,2)
      SUM(3)=SUM(3)+COR(J,3)
170     CONTINUE
        DO 180 J=1,NA
      COR(J,1)=COR(J,1)-SUM(1)/NA
      COR(J,2)=COR(J,2)-SUM(2)/NA
      COR(J,3)=COR(J,3)-SUM(3)/NA
180     CONTINUE

C --- order the nucleons on increasing z-coordinates
      DO 220 J1=2,NA
        ZCO=COR(J1,3)
        DO 210 J2=1,J1-1
        IF (ZCO.LT.COR(J2,3)) THEN
          XCO=COR(J1,1)
          YCO=COR(J1,2)
          DO 200 J3=1,J1-J2
            COR(J1+1-J3,1)=COR(J1-J3,1)
            COR(J1+1-J3,2)=COR(J1-J3,2)
            COR(J1+1-J3,3)=COR(J1-J3,3)
200       CONTINUE
          COR(J2,1)=XCO
          COR(J2,2)=YCO
          COR(J2,3)=ZCO
          GOTO 220
          ENDIF
210     CONTINUE
220   CONTINUE
      RETURN
      END

C********************************* END FRNUCOR ***************************

C********************************* FRNUCOD *******************************

      SUBROUTINE FRNUCOD(NA,RMIN2,A0,A2,A4,RMAX3,COR)

C --- this subroutine determines nucleon coordinates inside a deformed
C     nucleus and recenter the sampled nucleons
C     with respect to the rest frame of the nucleus center

      PARAMETER (KSZ2=300)
      DIMENSION COR(KSZ2,3),R0(38),FM(38),RM(38),SUM(3)
      DATA R0/5.8,5.9,6.0,6.1,6.2,6.3,6.4,6.5,6.6,6.7,6.8,6.9,7.0,7.1,
     $    7.2,7.3,7.4,7.5,7.6,7.7,7.8,7.9,8.0,8.1,8.2,8.3,8.4,8.5,8.6,
     $    8.7,8.8,8.9,9.0,9.1,9.2,9.3,9.4,9.5/
      DATA FM/  20.32,21.13,21.96,22.80,23.66,24.52,25.43,26.35,27.28,
     $    28.23,29.19,30.18,31.18,32.20,33.24,34.30,35.37,36.47,37.58,
     $    38.71,39.86,41.02,42.21,43.41,44.63,45.87,47.13,48.41,49.70,
     $    51.02,52.35,53.70,55.07,56.46,57.86,59.29,60.73,62.19/
      DATA RM/  16.49,16.59,16.70,16.81,16.91,17.02,17.13,17.23,17.34,
     $    17.45,17.56,17.66,17.77,17.87,17.98,18.09,18.19,18.30,18.41,
     $    18.51,18.62,18.73,18.83,18.94,19.04,19.15,19.26,19.36,19.47,
     $    19.57,19.68,19.79,19.89,20.00,20.10,20.21,20.32,20.42/

      DO 150 J=1,NA
C --- sample a nucleon from the target
C --- first sample sita from R(sita)**3. W1:cos(sita)
50    W1=-1.+2.*PYR(0)
      W2=RMAX3*PYR(0)
      W12=W1*W1
      W14=W12*W12
      RSITA=A0+A2*W12+A4*W14
      RSITA3=RSITA*RSITA*RSITA
      IF (RSITA3.LT.W2) GOTO 50
      CTHITA=W1
      STHITA=FRSQR(MAX(0.,1.-W12), 'STHW12')
      FAI=6.2832*PYR(0)

C --- then sample r

      RT0=RSITA
      FMM=FRINT(FM,R0,38,RT0)
      RMM=FRINT(RM,R0,38,RT0)

        DO 140 NUM=1,10
100   RR1=RMM*PYR(0)
      RR2=FMM*PYR(0)
      RR1S=RR1*RR1
      FR=RR1S/(1.+FRREX((RR1-RT0)/.55))
        IF (RR2.LT.FR) THEN
      R=RR1
        ELSE
      GOTO 100
        ENDIF

      COR(J,1)=R*STHITA*COS(FAI)
      COR(J,2)=R*STHITA*SIN(FAI)
      COR(J,3)=R*CTHITA
        IF (J.EQ.1) GOTO 150
C --- check if there are two nucleons too close each other
      DO 130 J1=1,J-1
      DICX=COR(J,1)-COR(J1,1)
      DICY=COR(J,2)-COR(J1,2)
      DICZ=COR(J,3)-COR(J1,3)
      DIC2=DICX*DICX+DICY*DICY+DICZ*DICZ
C --- if two nucleons too close each other, sample R once more
      IF (DIC2.LT.RMIN2) GOTO 140
130   CONTINUE
      GOTO 150
140   CONTINUE
C --- if 10 times repeated don't help, then sample THITA and FAI again
      GOTO 50
150   CONTINUE

C --- recenter the sampled nucleons within a nucleus
      SUM(1)=0.
      SUM(2)=0.
      SUM(3)=0.
        DO 160 J=1,NA
      SUM(1)=SUM(1)+COR(J,1)
      SUM(2)=SUM(2)+COR(J,2)
      SUM(3)=SUM(3)+COR(J,3)
160     CONTINUE
        DO 170 J=1,NA
      COR(J,1)=COR(J,1)-SUM(1)/NA
      COR(J,2)=COR(J,2)-SUM(2)/NA
      COR(J,3)=COR(J,3)-SUM(3)/NA
170     CONTINUE

C --- order the nucleons on increasing z-coordinates
      DO 220 J1=2,NA
        ZCO=COR(J1,3)
        DO 210 J2=1,J1-1
          IF (ZCO.LT.COR(J2,3)) THEN
          XCO=COR(J1,1)
          YCO=COR(J1,2)
          DO 200 J3=1,J1-J2
            COR(J1+1-J3,1)=COR(J1-J3,1)
            COR(J1+1-J3,2)=COR(J1-J3,2)
            COR(J1+1-J3,3)=COR(J1-J3,3)
200       CONTINUE
          COR(J2,1)=XCO
          COR(J2,2)=YCO
          COR(J2,3)=ZCO
          GOTO 220
          ENDIF
210     CONTINUE
220   CONTINUE


      RETURN
      END

C********************************* END FRNUCOD ***************************

C********************************* FRNUCDF ******************************

      SUBROUTINE FRNUCDF(NA,A0,A2,A4,RMAX3)

C --- this subroutine determines parameters of nucleon density
C       distribution of deformed nuclei

      PARAMETER (KSZ1=20)
      COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)
      SAVE /FRPARA1/

      BETA2S=VFR(4)*VFR(4)
      BETA4S=VFR(5)*VFR(5)
      R0=1.16*NA**.3333333-1.35*NA**(-.3333333)
      A0=R0*(1.-.0795774*(BETA2S+BETA4S)-.446031*VFR(4)+.44881*VFR(5))
      A2=R0*(1.338093*VFR(4)-4.4881*VFR(5))
      A4=R0*(5.236117*VFR(5))
      RMAX3=(A0+A2+A4)**3
      RETURN
      END

C********************************* END FRNUCDF **************************


C********************************* FRSEARC *******************************

       SUBROUTINE FRSEARC(L,fmax,xmin)

C...this is a routine for finding maximum and 'sufficient minimum'
C...of Woods-Saxon or harmonic oscillator.  xm, fm are outputs.
C...Input L=1 for proj, L=2 for target.

       dimension xx(4)
       xx(1)=101.
       do 10 j=1,3
        fmax=0.
        do 20 i=1,200
         x=(.01**(j-1))*(i-1)+xx(j)-100*(.01**(j-1))
       F = FRROR(L,X)
        if(f.gt.fmax) then
         fmax=f
         xx(j+1)=x
        endif
  20   continue
  10   continue

        x1=xx(4)
        fmin=fmax
       do 30 i=int(x1),50
        x=float(i)
       F = FRROR(L,X)
           if(f.lt.fmin) then
      fmin=f
      xy=x
         endif
        if(fmin.lt.1.E-3*fmax) goto 35
  30    continue

  35    do 40 i=1,200
         x=.01*float(i-1)+xy-1.
       F = FRROR(L,X)
            if(f.lt.fmin) then
       fmin=f
         xmin=x
            endif
  40    continue

        return
       end

C********************************* END FRSEARC ***************************

       REAL FUNCTION FRROR(L,R)

C......Gives value of nuclear density.  L=1,2 for projectile and target
C......the dinsity FRROR is given by
C...... for A<,=16, harmonic oscilator potential shell model density;
C...... for A> 16, (Woods-Saxon) Fermi distribution.

      PARAMETER (KSZ1=20)
      COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)
      CHARACTER*4 PACD
      COMMON/FRCODES/IPT(2),PACD(27),NNUC(27),NPROT(27),KCD(27)
     >           ,RO1(27,2),EXMA(9,2)
      SAVE /FRPARA1/,/FRCODES/

      A = NNUC(IPT(L))

      IF(A.LE.16) THEN
      RCH = RO1(IPT(L),1)
      D2 = (2.5-4./A)**(-1) * (RCH**2 - 0.81**2)
C.......(Proton radius 0.81 was subtracted off from the charge radius.)
      FRROR = (1.+((A-4.)/6.)* R**2/D2) * FRREX(-R**2/D2)

      ELSE
C                               (Woods-Saxon) distribution

       R0 = RO1(IPT(L),1)
       C =  RO1(IPT(L),2)
        AP=A**(1./3.)
        IF(R0.LE.0.) r0=1.16*(1.-1.16/ap**2)
        IF(C.LE.0.) C =0.5
       ARG = (R - R0*AP)/C
       FRROR=( 1.+FRREX(ARG) )**(-1)

       ENDIF

       FRROR = r**2 * FRROR

       return
       end


C********************************* FROVLAP ******************************


      SUBROUTINE FROVLAP

C...This subroutine determines the parameters of overlaping functions
C...for given cross sections SIGtot and SIGel.  Output is
C   AOP(3-4): OMEGA0,BETA ---- eikonal ovlap=1 - exp(-2Omega0*exp(-beta*b2))
C   AOP(5-6): GAUA,GAUB ---- Gaussian ovlap = 1-(1-gaexp(-gb*b2))**2
C   AOP(7-8) REC0,ALFA ---- gray disk = alfa  if b<REC0.
C   Fitting such that SIGtot = 2 integral d^b (1-sqrt(1-OVLAP));
C                       SIGinel =  integral d^b (OVLAP);
C...Units: AOP(4),AOP(6) are in Fm^-2; AOP(7) is in Fermi.


      PARAMETER (KSZ1=20,PI=3.1415926)
      COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)
      COMMON/FRINTN0/PLI0(2,4),AOP(KSZ1),IOP(KSZ1),NFR(KSZ1)
      dimension OMEG(2)
      SAVE /FRPARA1/,/FRINTN0/

      SGTOT = VFR(10)
      SGEL = VFR(11)
      SGTOT = SGTOT/0.389
      SGEL = SGEL/0.389
      SGINEL = SGTOT - SGEL
      RATIO = (SGINEL)/SGTOT

C............Eikonal parameters....................................

      OMEG(1) = 0.
      OMEG(2) = 2.
5      FOMEG = RATIO - 0.5*FRSUM(-2.*OMEG(2))/FRSUM(-OMEG(2))
      IF(ABS(FOMEG).LT.0.001) THEN
      OMEGTY = OMEG(2)
      GOTO 100
      ELSEIF(FOMEG.LT.0.) THEN
      OMEG(2) = OMEG(2) + 0.5
      GOTO 5
      ENDIF

      I = 0
10    I = I+1

      OMEGTY = (OMEG(1)+OMEG(2))/2.
      FOMEG = RATIO - 0.5*FRSUM(-2.*OMEGTY)/FRSUM(-OMEGTY)
      IF(ABS(FOMEG).LT.0.001) GOTO 100
      IF(FOMEG.LT.0.)  OMEG(1) = OMEGTY
      IF(FOMEG.GT.0.)  OMEG(2) = OMEGTY

      GOTO 10

100   BETA = (-2.*PI/SGTOT)* FRSUM(-OMEGTY)
      AOP(4) = BETA/(0.197)**2

      AOP(3) = OMEGTY

C............Gaussian parameters....................................

      AOP(5) = 4.*SGEL/SGTOT
      AOP(6) = 2.*PI*AOP(5)/SGTOT/(0.197)**2

C.............GRAY DISK PARAMETERS................................

      AOP(8) = 4.* (SGINEL)*SGEL/SGTOT**2
      REC02 = SGINEL/(PI*AOP(8))
      AOP(7) = SQRT(REC02) *0.197

      return

      END

C................................................

      FUNCTION FRSUM(X)

C...Summation FRSUM=SUM (x^n/n*n!), used for integrating eikonal overlap func.

      FRSUM = 0.
      IF(X.EQ.0.) RETURN

      I = 0
      TERM = 1.
10    I = I+1
      TERM = TERM * (X) * MAX(1,I-1)/FLOAT(I**2)
      FRSUM = FRSUM + TERM
      IF(ABS(TERM).LT.MIN(1.E-6,1.E-6*ABS(X))) GOTO 100

      CALL FRLOOPU(*10,I,2000,'LOPFRSUM')

100   CONTINUE

      RETURN
      END

C********************************* END FROVLAP **************************

C************************  FUNCTION FRVOV ********************************

      FUNCTION FRVOV(R2)

C...... Gives the value of overlap function at b^2=R2.
C...... KFR(5)=IQ=0, eikonal, IQ=1,Gaussian, IQ=2, grey disk.

      PARAMETER (KSZ1=20)
      COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)
      COMMON/FRINTN0/PLI0(2,4),AOP(KSZ1),IOP(KSZ1),NFR(KSZ1)
      SAVE /FRPARA1/,/FRINTN0/

      IF(KFR(5).EQ.0) THEN
      OMEGA0 = AOP(3)
      BETA = AOP(4)
      FRVOV = 1.-FRREX( -2.*OMEGA0*FRREX(-BETA*R2))
            ELSE IF (KFR(5).EQ.1) THEN
      GAUA = AOP(5)
      GAUB = AOP(6)
            FRVOV =1.-(1.-GAUA*FRREX(-GAUB*R2))**2
            ELSE IF (KFR(5).EQ.2) THEN
      REC02 = AOP(7)**2
      ALFA = AOP(8)
         IF (R2.LT.REC02) THEN
         FRVOV=ALFA
         ELSE
         FRVOV=0.
         ENDIF
      ENDIF

      RETURN
      END

C************************ END FUNCTION FRVOV ********************************

C*************************************************************************

      FUNCTION FRINT(FF,FX,NN,X)

CC........this function deals with the linear interpolation

      DIMENSION FF(NN),FX(NN)

      IF (X.LT.FX(1).OR.X.GT.FX(NN)) THEN
         CALL FRMGOUT(0,0,
     >   'X out of range! Be advised to set KFR(6)=0.',
     >    X,FX(1),FX(2),0.,0.)
      ENDIF
        DO 100 I=1,NN-1
      IF (X.GT.FX(I).AND.X.LE.FX(I+1)) THEN
      FRINT=FF(I)+(FF(I+1)-FF(I))*(X-FX(I))/(FX(I+1)-FX(I))
      goto 200
      ENDIF
100     CONTINUE

200   RETURN
      END

C*************************************************************************

C******************* END OF PACKAGE FOR NUCLEAR GEOMETRY *****************




C*************************************************************************
C                                                                        *
C This is the routine package for generating nucleon-nucleon collisions  *
C                                                                        *
C*************************************************************************

C********************************* FRRINGO *******************************

      SUBROUTINE FRRINGO

C ....................    THE ROUTINE GIVES MASSES TO THE EXCITED NUCLEONS

      IMPLICIT DOUBLE PRECISION (D)
      PARAMETER (KSZ1=20,KSZ2=300,KSZJ=4000)
      COMMON/FRINTN0/PLI0(2,4),AOP(KSZ1),IOP(KSZ1),NFR(KSZ1)
      COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)
      COMMON/FRINTN3/IDN(2,KSZ2),FMN(2,KSZ2),NUC(2,3000)
      COMMON/FRINTN1/PPS(2,KSZ2,5),PPH(2,KSZ2,5),PPSY(2,KSZ2,5),PPA(2,5)
      COMMON/FRINTN2/NHP(2),IHQP(2,KSZ2),KHP(2,KSZ2,100,5),
     >   PHP(2,KSZ2,100,5)
      COMMON/FRCNUT/NR,KR(10,5),PR(10,5),NR0
      COMMON/FRJETS/NJ,KJ(100,5),PJ(100,5)
      COMMON/RYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)

      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      integer N,NPAD,K
      double precision P,V
      SAVE /PYJETS/

      COMMON/RYSUBS/MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200)

      DIMENSION DPTSQ(2),DP(2,5),DN(2,5),DSM(2,5),DSYM(2,5),DBETA(3),
     > PJK(2,4),PPSYU(2,5),FMN0(2),RM(2,2)
      SAVE PXGAU,PYGAU,RM
      SAVE /FRINTN0/,/FRPARA1/,/FRINTN3/,/FRINTN1/,/FRINTN2/,/FRCNUT/,
     >     /FRJETS/,/RYPARS/,/RYSUBS/
      DATA RM /0.,0.,0.,0./
      CALL FRDOICT(-1)

C..Monitors ICT(10) (for current event) and ICTT (for all events) are set up:
C...   ICT(kfel) counts the number of times error KFEL occurs.
C...   ICT(7) is the number of ``single diff'' collisions.
C...   ICT(8) the number of collisions bypassed due to repeated errors.
C...  KFEL=2,3,8,10 result in the collision being skipped.
C...  KFEL used: 1-10

C...........................................set initial momenta

 10   KFEL=0
      DO  L=1,2
        DO  I=1,IOP(3+2*(L-1))
c..................iop(3/5)=number of projectile/target nucleons
          DO  LO=1,4
            PPS(L,I,LO)= PLI0(L,LO)
c.....pli0 stores the initial momenta of the colliding hadron/nucleon
c.....pps stores th initial momenta of each nucleon.
          enddo
          PPS(L,I,5)= FRSQR(PPS(L,I,4)*PPS(L,I,3)-PPS(L,I,1)**2
     >         -PPS(L,I,2)**2, 'reyu90')
          IHQP(L,I) = 0
          kHP(L,I,1,5)=0
          PHP(L,I,1,5)=-1.
          PHP(L,I,2,5)=-1.
          DO  LO=1,5
            PPH(L,I,LO) = 0.
            PPSY(L,I,LO) = PPS(L,I,LO)
          enddo
        enddo
      enddo
C....Spectartor momenta and FERMI-MOTION  ............

      CALL FRHELGE

C.....loop over all binary collisions

      IOP(13) = 0
      IOP(14) = 0
      IOP(15) =0
      NHP(1) = 0
      NHP(2) = 0
      NR = 0
      IF(IOP(2).GT.3000) CALL FRMGOUT(0,0,
     >  '** NUC(2,3000) must be increased **',0.,0.,0.,0.,0.)

      DO 1000 I=1,IOP(2)

      IOP(1) = I
      NU1 = NUC(1,I)
      NU2 = NUC(2,I)
      KFELS=0

C.....FRVECTC SETS FOUR-VECTORS DP(2,4) TO PPS .....
      CALL FRVECTC(I, 1, DP)
      CALL FRVECTC(I, 1, DSM)
      CALL FRVECTC(I, 2, DSYM)
C...............................TOTAL E-P BEFORE THE COLLISION

      DO  LO = 1, 4
        PPSYU(1,LO) = PPSY(1,NU1,LO)+PPSY(2,NU2,LO)
      enddo

C.....BOOST to remnent-remnent CMS ......................

      CALL FRTOCMS(1, 1, DP, DBETA)

C....Skip the collision when strings moving backwards in CMS:
      IF(DP(1,4)-DP(1,3).LE.0.D0) THEN
C     KFEL=10
      CALL FRDOICT(10)
      GOTO 1000
      ENDIF

C.... ANGLES OF THE MOMENTUM VECTORS ................................

      CALL FRPOLAR(DTHE,DPHI, DP)

C...  ROTATE SO P GOES TO THE Z-AXE ...................

      CALL FRROTAR(DTHE,DPHI,1, DP)
C....
      NR0 = NR+1

      SMP= REAL( DP(1,3)*DP(1,4)-DP(1,1)**2-DP(1,2)**2)
      SMT= REAL( DP(2,3)*DP(2,4)-DP(2,1)**2-DP(2,2)**2)
c......smp/smt=masses of projectile/target

      W=REAL( DP(1,4)+DP(2,4) )
      PK2= FRKVM(W,SMP,SMT)

50    IF(IOP(15).EQ.0) THEN
        DO  LL=1,2
          DO  LO=1,4
            PJK(LL,LO) = 0.
          enddo
        enddo
      ENDIF

C............................................GENERATE HARD PARTONS.......

      IHAV = 0
      NJ = 0
      IF(IOP(18).GE.1.AND.IOP(15).EQ.0) then

c [KG] :        FMN0(1) = ULMASS(IDN(1,NU1))
c [KG] :        FMN0(2) = ULMASS(IDN(2,NU2))
         FMN0(1) = ULGAMASS(1,NU1) ! [KG]
         FMN0(2) = ULGAMASS(2,NU2) ! [KG]

C...................................effective energy for hard scattering
        P13 = FMN0(1)**2/REAL(DP(1,4))
        P24 = FMN0(2)**2/REAL(DP(2,3))
c..p13=m^1**2/p_1+
c..p24=m^2**2/p_2-
        E1 = 0.5*( REAL(DP(1,4))+P13 )
        P1 = 0.5*( REAL(DP(1,4))-P13 )
        E2 = 0.5*( REAL(DP(2,3))+P24 )
        P2 = 0.5*(-REAL(DP(2,3))+P24 )
c..E1=0.5*( (p_1+) + m^1**2/p_1+)
c..p1=0.5*( (p_1+) - m^1**2/p_1+)
c..E2=0.5*( (p_2-) + m^2**2/p_2-)
c..p2=0.5*(-(p_2-) - m^2**2/p_2-)

        CKIN(22)=MIN(1.,REAL(DP(1,3))/P13)
        CKIN(24)=MIN(1.,REAL(DP(2,4))/P24)

        WEF =FRSQR(FMN0(1)**2+FMN0(2)**2+2.*(E1*E2-P1*P2),' WEFIS1')
        IF(WEF.LT.PARP(2)) GOTO 495

        N=0

        P(1,1) = 0.
        P(1,2) = 0.
        P(2,1) = 0.
        P(2,2) = 0.
        P(1,3) = P1
        P(2,3) = P2
        P(1,4) = E1
        P(2,4) = E2

C.....PYTHIA is reinitialized whenever the collision enviroment is changed.
C.....if momenta unchanged, no need to reinitialize PYTHIA, save time.....
        INI = 2
        DO  LO=1,2
          DO  J=3,4
            IF( ABS(P(LO,J)-RM(LO,J-2)).GT.0.01*RM(LO,J-2) ) INI=1
          enddo
        enddo
        IF(INI.EQ.1) THEN
          DO LO=1,2
            DO J=1,2
              RM(LO,J) = P(LO,J+2)
            enddo
          enddo
C.....resetparametersbefore every new collision:
          CALL FRSETRY(1)
        ENDIF

C.....PYTHIA is reinitialized here only when the collision energy changes, but
C.....not when the particle changes (neutron instead of a proton, for exemple).

        NFR(2) = NFR(2)+1
        CALL FRHARDP(IDN(1,NU1),IDN(2,NU2),Wef,IHAV,INI)
        IF(IHAV.EQ.0)     GOTO 495
        DO J=1, NJ
          PJK(ABS(KJ(J,3)),1) = PJK(ABS(KJ(J,3)),1)+ PJ(J,1)
          PJK(ABS(KJ(J,3)),2) = PJK(ABS(KJ(J,3)),2)+ PJ(J,2)
          PJK(ABS(KJ(J,3)),3) = PJK(ABS(KJ(J,3)),3)+ PJ(J,4) - PJ(J,3)
          PJK(ABS(KJ(J,3)),4) = PJK(ABS(KJ(J,3)),4)+ PJ(J,4) + PJ(J,3)
 220    enddo
      ENDIF
C........................................................................


C..........generate soft PT in REMNENT-REMNENT CMS frame ......

495   ICPK=0
      Irep = 0

      AHT=0.
      IF(IOP(15).GE.2.AND.KFR(9).EQ.1) AHT=1.0

600   PXGAU = 0.
      PYGAU = 0.

      IF(ICPK.LE.100) THEN
      PX0 = AHT*PJK(1,1)
      PY0 = AHT*PJK(1,2)
      CALL FRCOLPT(PK2,PXGAU,PYGAU,PX0,PY0)
      ICPK = ICPK+1
        IF(IOP(15).EQ.0.AND.NJ.GT.0) THEN
        PTTRY = (PXGAU+PJK(1,1))**2+ (PYGAU+PJK(1,2))**2
        IF(PTTRY.GT.PK2) GOTO 600
        ENDIF
      ENDIF

500   DN(1,1) = DBLE(PXGAU)
      DN(1,2) = DBLE(PYGAU)
      DN(2,1) = -DBLE(PXGAU)
      DN(2,2) = -DBLE(PYGAU)

      DPTSQ(1) = DN(1,1)**2 + DN(1,2)**2
      DPTSQ(2) = DN(2,1)**2 + DN(2,2)**2

C.............................................................

      IF(KFEL.GT.0) THEN
      KFELS=KFELS+1
      CALL FRDOICT(KFEL)
      ENDIF

C  Bypass the collision if too many errors:
      IF (KFELS.GT.20) THEN
      IF(I.EQ.1) GOTO 10
      KFEL=8
      GOTO 999
      ENDIF

C............... SOFT MOMENTUM TRANSFERS ...........................

610   CALL FRPSOFT(I,PJK,DP,DN,KFEL)
      IF(KFEL.GT.0) GOTO 999

C.....excited masses for soft remnents
      DWI =DN(1,4)*DN(1,3)
      DWT =DN(2,4)*DN(2,3)
      IF(DWI.LT.DPTSQ(1).OR.DWT.LT.DPTSQ(2)) THEN
      kfel=1
      GOTO 500
      ENDIF

      DN(1,5)=DFRSQR(DWI-DPTSQ(1), 'DWI123')
      DN(2,5)=DFRSQR(DWT-DPTSQ(2), 'DNW935')


C.....transform DN & PJ() & PR() back to original frame

      IF(IHAV.EQ.1.AND.IREP.EQ.0) THEN
      IQL = -2
      ELSE
      IQL = -1
      ENDIF

       CALL FRROTAR(DTHE,DPHI,IQL, DN)
       CALL FRTOCMS(1, IQL, DN, DBETA)
       Irep = 1

C......... UPDATE THE PPS(1-2,I,4) ARRAY:
      CALL FRVECTC(I, -1, DN)

      PPS(1,NU1,5) =REAL(DN(1,5))
      PPS(2,NU2,5) =REAL(DN(2,5))

C..................................update PPSY ..........
      CALL FRFILHD(I,0,KFEL)
      IF(KFEL.GT.0) GOTO 500

C.....treatment of diffractive collision

      CALL FRSETDM(I,KFEL)
C     ...single diffractive
      IF(KFEL.EQ.-1) CALL FRDOICT(7)
      IF(KFEL.GE.5)  GOTO 500

C....................STORE THE HARD PARTONS TO FRINTN2 ..........

      IF(IOP(15).EQ.0) CALL FRFILHD(I,1,KFEL)

C...TO test the Rutherford scattering against the background:

       N=0
       IF(NJ.GT.0.AND.IOP(15).EQ.0)
     >       CALL FRCHEXG(*50,I)

       IF(IOP(15).LE.1.AND.NJ.GT.0) THEN
         IOP(13) = IOP(13) + 1
         IOP(14) = IOP(14) + (NJ+1)/2
       ENDIF

 990   KFEL=0
       IOP(15)=0
       NJ=0

C.........CHECK E-P CONSERVATION .....................
      DO 910 LO = 1,4
910   PPSYU(2,LO) = PPSY(1,NU1,LO)+PPSY(2,NU2,LO)

      DO 920 LO = 4, 1, -1
      PERR = ABS(PPSYU(2,LO)-PPSYU(1,LO))
      IF(PERR.GE.MAX(0.1+0.5*(LO/3),0.05*PPSYU(1,LO)) ) THEN
      CALL FRMGOUT(1,1,'FRRINGO E-P not conserved:',
     >  PPSYU(2,1)-PPSYU(1,1),PPSYU(2,2)-PPSYU(1,2),
     >  PPSYU(2,3)-PPSYU(1,3),PPSYU(2,4)-PPSYU(1,4),PPSYU(1,4))
      GOTO 1000
      ENDIF
920   CONTINUE

      GOTO 1000
999   CALL FRVECTC(I, -1, DSM)
      CALL FRVECTC(I, -2, DSYM)
      CALL FRDOICT(KFEL)
1000  CONTINUE

      RETURN
      END


C-------------------------------------------------------------------

      FUNCTION FRKVM(W,AM1,AM2)

C..FRKVM=3-MOMENTUM K^2 IN A CMS FRAME OF TWO PARTICLES WITH MASS AM1 & AM2:

      FRKVM = ((W**2-AM1-AM2)**2-4.*AM1*AM2)/(4.*W**2)
      RETURN
      END


C********************************* END FRRINGO ***************************


C********************************* FRPSOFT *******************************

      SUBROUTINE FRPSOFT(I,PJK,DP,DN,KFEL)

C.....To generate momenta for the soft remnents
C.....  IQ=0: generated according to dQ/Q
C.....    =1: dQ/Q+k
C.....  kfel>0: no phase space for the collision

      PARAMETER (KSZ1=20,KSZ2=300)
      IMPLICIT DOUBLE PRECISION (D)
      COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)
      COMMON/FRINTN3/IDN(2,KSZ2),FMN(2,KSZ2),NUC(2,3000)
      COMMON/FRINTN0/PLI0(2,4),AOP(KSZ1),IOP(KSZ1),NFR(KSZ1)
      COMMON/FRINTN1/PPS(2,KSZ2,5),PPH(2,KSZ2,5),PPSY(2,KSZ2,5),PPA(2,5)
      COMMON/FRJETS/NJ,KJ(100,5),PJ(100,5)
      DIMENSION DP(2,5),DN(2,5),PJK(2,4)
      SAVE /FRPARA1/,/FRINTN3/,/FRINTN0/,/FRINTN1/,/FRJETS/

      KFEL=0
      IK9=0
      IQ = 1
      IF(IOP(15).GE.2) IQ =0
      NU1 = NUC(1,I)
      NU2 = NUC(2,I)

      DWP = DP(1,4)+DP(2,4)+PPH(1,NU1,4)+PPH(2,NU2,4)
      DWM = DP(1,3)+DP(2,3)+PPH(1,NU1,3)+PPH(2,NU2,3)
      DA =(DN(1,1)+PPH(1,NU1,1)+DBLE(IQ*PJK(1,1)))**2+
     >    (DN(1,2)+PPH(1,NU1,2)+DBLE(IQ*PJK(1,2)))**2
      DB =(DN(2,1)+PPH(2,NU2,1)+DBLE(IQ*PJK(2,1)))**2+
     >    (DN(2,2)+PPH(2,NU2,2)+DBLE(IQ*PJK(2,2)))**2
      DA = DA + AOP(9)**2
      DB = DB + AOP(10)**2

10    CALL FRPLIMT(DWP,DWM,DA,DB,DPLO3,DPHI3,DPLO4,DPHI4,KFEL)
      IF(KFEL.GT.0) RETURN

      DPLO3 = DMAX1(DPLO3, DP(1,3)+PPH(1,NU1,3))
      DPLO4 = DMAX1(DPLO4, DP(2,4)+PPH(2,NU2,4))
      DPHI3 = DMIN1(DPHI3, DP(1,4)+PPH(1,NU1,4))
      DPHI4 = DMIN1(DPHI4, DP(2,3)+PPH(2,NU2,3))
      IF(IQ.EQ.1) THEN
      DPLO3 = DMAX1(DPLO3, DBLE(PJK(1,3)+PPH(1,NU1,3)) )
      DPLO4 = DMAX1(DPLO4, DBLE(PJK(2,4)+PPH(2,NU2,4)) )
      ENDIF

      P0 = 0.

      IF(PYR(0).LT.0.500) THEN

            IF(DPHI3.LT.DPLO3) GOTO 99
      P0 = -PPH(1,NU1,3)
      DN(1,3) = DFRDPQ(DPLO3,DPHI3,P0)

      DPLO4 = DMAX1(DPLO4,DB/(DWM-DN(1,3)))
      DPHI4 = DMIN1(DPHI4, DWP - DA/DN(1,3) )
            IF(DPHI4.LT.DPLO4) GOTO 99
      P0 = -PPH(2,NU2,4)
      DN(2,4) = DFRDPQ(DPLO4,DPHI4,P0)

      ELSE

            IF(DPHI4.LT.DPLO4) GOTO 99
      P0 = -PPH(2,NU2,4)
      DN(2,4) = DFRDPQ(DPLO4,DPHI4,P0)

      DPLO3 = DMAX1(DPLO3,DA/(DWP-DN(2,4)))
      DPHI3 = DMIN1(DPHI3, DWM - DB/DN(2,4) )
           IF(DPHI3.LT.DPLO3) GOTO 99
      P0 = -PPH(1,NU1,3)
      DN(1,3) = DFRDPQ(DPLO3,DPHI3,P0)

      ENDIF

      DN(1,4) = DWP - DN(2,4)
      DN(2,3) = DWM - DN(1,3)

      DO 70 L=1,2
      DO 70 LO=3,4
      DN(L,LO) = DN(L,LO) - DBLE(PPH(L,NUC(L,I),LO))
      IF(IQ.EQ.1) DN(L,LO) = DN(L,LO) - DBLE(PJK(L,LO))
        IF(DN(L,LO).LE.0.D0) THEN
      IK9=IK9+1
      IF(IK9.LT.10) GOTO 10
      KFEL=9
        ENDIF
70    CONTINUE

      GOTO 100
99    KFEL=2
100   RETURN
      END

C********************************* END FRPSOFT ***************************

C********************************* FUNCTION FRDPQ *************************

      DOUBLE PRECISION FUNCTION DFRDPQ(DPMIN,DPMAX,P0)

C......generate PQ (P- or P+ ) according to dPQ/PQ+P0,
C......with Pmin< PQ <Pmax.

      IMPLICIT DOUBLE PRECISION (D)
      PARAMETER (KSZ1=20)
      COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)
      SAVE /FRPARA1/

      DP0 = DBLE(P0)
      IF(DPMIN+DP0.LE.0.) THEN
      DFRDPQ = 0.
      RETURN
      ENDIF

        IF(DPMAX.LE.DPMIN) THEN
      DFRDPQ = DPMIN
      ELSE
      DR = DBLE(PYR(0))
c....... longitudinal momentum transfer
c....changed by jochen geiss
      DFRDPQ = (DPMIN+DP0) *((DPMAX+DP0)/(DPMIN+DP0))**DR -DP0
c      DFRDPQ = (DPMIN+DP0) *((DPMAX+DP0)/(DPMIN+DP0))**(DR**1.5)-DP0
      ENDIF

      RETURN
      END


C****************************** END FUNCTION FRDPQ ***********************


C********************************* FRPLIMT ******************************

C....TO GIVE THE UPPER AND LOWER LIMITS FOR FINAL MOMENTA
C....  DP(2,4)  -  initial momenta
C...   DA=MINIMUM OF M_3T**2; DB=MINIMUM OF M_4T**2.

      SUBROUTINE FRPLIMT(DWP,DWM,DA,DB,DPLO3,DPHI3,DPLO4,DPHI4,KFEL)

      IMPLICIT DOUBLE PRECISION (D)

      KFEL=0
      DS = DWP*DWM
      DTM1 = (DS + DA - DB)
      DTM2 = (DTM1**2-4.D0*DA*DS)
            IF(DTM2.LT.0) THEN
      KFEL=3
      RETURN
            ENDIF
      DTM2 = DSQRT( DTM2 )
      DPLO3 = (DTM1 - DTM2)/(2.D0*DWP)
      DPHI3 = (DTM1 + DTM2)/(2.D0*DWP)
      DPLO4 = (DB-DA+ DWP*DPLO3)/DWM
      DPHI4 = (DB-DA+ DWP*DPHI3)/DWM

      RETURN
      END

C********************************* END FRPLIMT **************************

C********************************* FRCHEXG ***************************

      SUBROUTINE FRCHEXG(*,I)

C.....To administrate the preliminary checking and testing of RPS
C.....This routine will manages the final gluon emission if there is only
C.....one collision, else the gluon emission will be managed in FRINGEB.
C...  IOP(15)>=2 signals hard partons are drawned.

      IMPLICIT DOUBLE PRECISION (D)
      PARAMETER (KSZJ=4000,KSZ1=20,KSZ2=300)
      COMMON/FRINTN0/PLI0(2,4),AOP(KSZ1),IOP(KSZ1),NFR(KSZ1)
      COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)
      COMMON/FRINTN3/IDN(2,KSZ2),FMN(2,KSZ2),NUC(2,3000)
      COMMON/FRJETS/NJ,KJ(100,5),PJ(100,5)

      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      integer N,NPAD,K
      double precision P,V
      SAVE /PYJETS/

      SAVE /FRINTN0/,/FRPARA1/,/FRINTN3/,/FRJETS/

       IF(IOP(15).EQ.0.AND.NJ.GT.0) IOP(15)=1

       CALL FRTORST(1,NUC(1,I))
       IF(IOP(15).NE.2) CALL FRTORST(2,NUC(2,I))

       IF(IOP(15).GE.2) THEN

       CALL FRFILHD(I,-1,KFEL)

       RETURN 1

       ENDIF

       RETURN
       END

C************************ END FRCHEXG ***********************************


C********************************* FRFILHD ******************************

      SUBROUTINE FRFILHD(I,IQ,KFEL)

C....TO INCORPORATE THE HARD PARTON MOMENTA INTO THE NUCLEON SYSTEM AND
C....TO EVELUATE THE SYSTEM MASS PPSY(,,5).
C......IQ = 0, PPSY is updated but hard partons not stored (PPH not updated);
C......IQ = 1, no effect to PPSY but PPH, PHP are updated
C......IQ <0: hard partons stripped off from the record.
C......Output flag for IQ=0: (kfel is dummy for ABS(IQ)=1)
C......Kfel=0, no problem;
C......    =4, system mass smaller than minimum, FRFILHD aborts;

      PARAMETER (KSZ1=20,KSZ2=300)
      COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)
      COMMON/FRINTN0/PLI0(2,4),AOP(KSZ1),IOP(KSZ1),NFR(KSZ1)
      COMMON/FRINTN3/IDN(2,KSZ2),FMN(2,KSZ2),NUC(2,3000)
      COMMON/FRJETS/NJ,KJ(100,5),PJ(100,5)
      COMMON/FRINTN2/NHP(2),IHQP(2,KSZ2),KHP(2,KSZ2,100,5),
     >   PHP(2,KSZ2,100,5)
      COMMON/FRINTN1/PPS(2,KSZ2,5),PPH(2,KSZ2,5),PPSY(2,KSZ2,5),PPA(2,5)

      DIMENSION PPHN(2,4),IHQPM(2)
      SAVE PPHN
      SAVE /FRPARA1/,/FRINTN0/,/FRINTN3/,/FRJETS/,/FRINTN2/,/FRINTN1/

      kfel=0

      IF(IQ.EQ.0) THEN

      DO 11 L=1,2
      DO 11 J=1,4
11    PPHN(L,J) = 0.

      IF(NJ.GE.1) THEN
      DO 21 LO=1, NJ
      L = ABS(KJ(LO,3))
      PPHN(L,1) = PPHN(L,1)+ PJ(LO,1)
      PPHN(L,2) = PPHN(L,2)+ PJ(LO,2)
      PPHN(L,3) = PPHN(L,3)+ PJ(LO,4)-PJ(LO,3)
21    PPHN(L,4) = PPHN(L,4)+ PJ(LO,4)+PJ(LO,3)
      ENDIF

      DO 30 L=1, 2
      DO 30 LO=1, 4
30    PPSY(L,NUC(L,I),LO) = PPS(L,NUC(L,I),LO)+
     >               PPH(L,NUC(L,I),LO)+ PPHN(L,LO)

      DO 570 L=1,2
      SMSY2 = PPSY(L,NUC(L,I),4)*PPSY(L,NUC(L,I),3)-
     >    PPSY(L,NUC(L,I),1)**2-PPSY(L,NUC(L,I),2)**2
      IF(SMSY2.LT.AOP(8+L)**2) THEN
      KFEL=4
      RETURN
      ENDIF
570   PPSY(L,NUC(L,I),5) = SQRT(SMSY2 )

      ELSEIF(IQ.EQ.1.and.NJ.GT.0) THEN
C.....................................STORE THE HARD PARTONS TO FRINTN2

       IHQPM(1) = IHQP(1,NUC(1,I))
       IHQPM(2) = IHQP(2,NUC(2,I))
       DO 512 LO = 1, NJ
       ISIDE = ABS(KJ(LO,3))
       INUC = NUC(ISIDE,I)
       IHQP(ISIDE,INUC) = IHQP(ISIDE,INUC) +1
       NHP(ISIDE) = NHP(ISIDE) + 1
       DO 510 L=1,4
       PHP(ISIDE,INUC,IHQP(ISIDE,INUC),L)= PJ(LO,L)
510    KHP(ISIDE,INUC,IHQP(ISIDE,INUC),L)= KJ(LO,L)
       IF(kfr(9).NE.0) KHP(ISIDE,INUC,IHQP(ISIDE,INUC),4)=IHQPM(ISIDE)
512    CONTINUE

       DO 90 L = 1, 2
       DO 90 L2 = 1,4
 90    PPH(L,NUC(L,I),L2) = PPH(L,NUC(L,I),L2)+ PPHN(L,L2)

       ELSEIF(IQ.LE.-1.and.NJ.GT.0) THEN
C....................................STRIP OFF THE HARD PARTONS

        DO 520 LO = 1, NJ
        ISIDE = ABS(KJ(LO,3))
        INUC = NUC(ISIDE,I)
         DO 517 L=1,4
         PHP(ISIDE,INUC,IHQP(ISIDE,INUC),L)= 0.0
517      KHP(ISIDE,INUC,IHQP(ISIDE,INUC),L)= 0.0
        NHP(ISIDE) = NHP(ISIDE) - 1
520     IHQP(ISIDE,INUC) = IHQP(ISIDE,INUC) -1

        DO 95 L = 1, 2
        DO 95 L2 = 1,4
 95     PPH(L,NUC(L,I),L2) = PPH(L,NUC(L,I),L2)- PPHN(L,L2)

      NJ = 0

      ENDIF

      RETURN
      END

C********************************* END FRFILHD **************************

C********************************* FRSETDM ******************************

      SUBROUTINE FRSETDM(I,KFEL)

C...TO RESET THE MASS IF FOUND DIFFRACTIVE, M<AOP(10+L).
C...I is the index for the collision
C...kfel = 0: not a "diffractive" event, no reset necessary
C....... = 5: fail to properly reset, regenerate the event may be needed.
C....... = 6: double diffractive
C....... = -1: reset successfully

      IMPLICIT DOUBLE PRECISION (D)
      PARAMETER (KSZ1=20,KSZ2=300)
      COMMON/FRINTN0/PLI0(2,4),AOP(KSZ1),IOP(KSZ1),NFR(KSZ1)
      COMMON/FRINTN3/IDN(2,KSZ2),FMN(2,KSZ2),NUC(2,3000)
      COMMON/FRINTN1/PPS(2,KSZ2,5),PPH(2,KSZ2,5),PPSY(2,KSZ2,5),PPA(2,5)
      DIMENSION PTSQ(2),PTSQS(2),PNEW(2)
      SAVE /FRINTN0/,/FRINTN3/,/FRINTN1/

      IF(PPSY(1,NUC(1,I),5).LT.AOP(11).and.
     >                 PPSY(2,NUC(2,I),5).LT.AOP(12)) THEN
C      "double diffractive"
      KFEL=6
      RETURN
      ELSEIF(PPSY(1,NUC(1,I),5).GE.AOP(11).and.
     >                 PPSY(2,NUC(2,I),5).GE.AOP(12)) THEN
      KFEL=0
      RETURN
      ELSEIF(PPSY(1,NUC(1,I),5).LT.AOP(11)) THEN
      KFEL=-1
             L = 1
      ELSEIF(PPSY(2,NUC(2,I),5).LT.AOP(12)) THEN
      KFEL=-1
             L = 2
      ENDIF

C-------RESET THE MASS -------------------------------------------------
C.......Keep Pt fixed and choose P_large between pmin and pmax acc dP/P,
C.......which is equivalent to a uniform distributn of y.
C.......(A bad choice here may result in an ugly peak in dn/dy_proton)

      NV = NUC(L,I)
      NVV = NUC(3-L,I)
c [KG] :      FM = ULMASS(IDN(L,NUC(L,I)))
      FM = ULGAMASS(L,NUC(L,I)) ! [KG]
      LG=4
      IF( PPSY(L,NV,3).GT.PPSY(L,NV,4) ) LG=3

      PTSQ(L) = PPSY(L,NV,1)**2+PPSY(L,NV,2)**2
      PTSQS(L) = PPS(L,NV,1)**2+PPS(L,NV,2)**2
      TMP20 = FM**2 + PTSQ(L)
      DPMAX = PPSY(L,NV,LG)
      DPMIN = TMP20/PPSY(L,NV,7-LG)
         PNEW(LG-2) = DFRDPQ(DPMIN,DPMAX,0.)
         PNEW(5-LG) = TMP20/PNEW(LG-2)
      ADELP = PNEW(2) - PPSY(L,NV,4)
      ADELM = PNEW(1) - PPSY(L,NV,3)
      PPSY(L,NV,4) = PNEW(2)
      PPSY(L,NV,3) = PNEW(1)
      PPS(L,NV,4) = PPS(L,NV,4)  + ADELP
      PPS(L,NV,3) = PPS(L,NV,3)  + ADELM
      RMS34S = PPS(L,NV,3)*PPS(L,NV,4)
          IF(RMS34S.LE.PTSQS(L)) GOTO 500
      PPSY(L,NV,5) =FM
      PPS(L,NV,5) =FRSQR(RMS34S-PTSQS(L),'ppslnv5')

C   Rebalance the energy-momentum

      PTSQ(3-L) = PPSY(3-L,NVV,1)**2+PPSY(3-L,NVV,2)**2
      PTSQS(3-L) = PPS(3-L,NVV,1)**2+PPS(3-L,NVV,2)**2
      PPSY(3-L,NVV,4) = PPSY(3-L,NVV,4) - ADELP
      PPSY(3-L,NVV,3) = PPSY(3-L,NVV,3) - ADELM
      PPS(3-L,NVV,4) = PPS(3-L,NVV,4) - ADELP
      PPS(3-L,NVV,3) = PPS(3-L,NVV,3) - ADELM
      RMSPT = AOP(12-L+1)**2 + PTSQ(3-L)
      RMS34 = PPSY(3-L,NVV,4)*PPSY(3-L,NVV,3)
      RMS34S = PPS(3-L,NVV,4)*PPS(3-L,NVV,3)
          IF(RMS34.LT.RMSPT.OR.RMS34S.LE.PTSQS(3-L)) GOTO 500
      PPSY(3-L,NVV,5)=FRSQR(RMS34-PTSQ(3-L), 'DN25IS')
      PPS(3-L,NVV,5)=FRSQR(RMS34S-PTSQS(3-L), 'DN25SS')

      GOTO 600
500   KFEL=5
600   RETURN
      END


C********************************* END FRSETDM **************************


C********************************* FRVECTC ******************************


      SUBROUTINE FRVECTC(I,IQ,DP)

C......CONVERSIONS BETWEEN PPS OR PPSY ARRAYS AND DP(2,4):
C.........TO SET PROJECTILE AND TARGET FOUR VECTORS TO A FORM:
C.........DP(L,1-4)=(PX,PY,P_,P+) ............................
C.........IQ=1, DP=PPS;
C.........IQ=2, DP=PPSY;
C.........FOR IQ < 0, THE REVERSE IS DONE.

      PARAMETER (KSZ2=300)
      COMMON/FRINTN1/PPS(2,KSZ2,5),PPH(2,KSZ2,5),PPSY(2,KSZ2,5),PPA(2,5)
      COMMON/FRINTN3/IDN(2,KSZ2),FMN(2,KSZ2),NUC(2,3000)
      DOUBLE PRECISION DP(2,5)
      SAVE /FRINTN1/,/FRINTN3/

      IF(IQ.EQ.1) THEN
        DO  L = 1, 2
          DO  J = 1, 4
            DP(L,J) =DBLE(PPS(L,NUC(L,I),J))
          enddo
        enddo
      ELSEIF(IQ.EQ.2) THEN
        DO  L = 1, 2
          DO  J = 1, 4
            DP(L,J) =DBLE(PPSY(L,NUC(L,I),J))
          enddo
        enddo
      ELSEIF(IQ.EQ.-1) THEN
        DO  L = 1, 2
          DO  J = 1, 4
            PPS(L,NUC(L,I),J)=REAL(DP(L,J))
          enddo
        enddo
      ELSEIF(IQ.EQ.-2) THEN
        DO  L = 1, 2
          DO  J = 1, 4
            PPSY(L,NUC(L,I),J)=REAL(DP(L,J))
          enddo
        enddo

      ENDIF

      RETURN
      END

C********************************* END FRVECTC **************************


C********************************* FRCOLPT ******************************

      SUBROUTINE FRCOLPT(PK2M,PX,PY,PX0,PY0)

C-----------------------------------------------------------------------
C     GIVING THE EXCITED NUCLEONS gaussian PT:
C      dP ~ exp(-(Px-Px0)**2/sig)*exp(-(Py-Py0)**2/sig)
C     PK2M: the upper cut off for PT^2.
C-----------------------------------------------------------------------

      PARAMETER (KSZ1=20)
      implicit DOUBLE PRECISION (D)
      COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)
      SAVE /FRPARA1/

      PX = 0.
      PY = 0.
      PK0 = PX0**2+PY0**2
      P2MX = PK2M-PK0
      IF(P2MX.LE.1.E-5) RETURN

      IF(VFR(6).LE.0.000001) THEN
      PX = PX0
      PY = PY0

      ELSE

      ITRY = 0
10    CALL FRGAUSS(P, VFR(6), P2MX)

      Adelpt=FRSQR(p, 'PIO086')

      AFI=2*3.1415926*PYR(0)
      PX=PX0+ADELPT*COS(AFI)
      PY=PY0+ADELPT*SIN(AFI)
        IF(PX**2+PY**2.GT.PK2M) THEN
        ITRY = ITRY+1
        CALL FRLOOPU(*10,ITRY,100,'LPfrcolpt')
        ENDIF

      ENDIF

      RETURN
      END


C********************************* END FRCOLPT **************************


C********************************* FRHELGE ******************************

      SUBROUTINE FRHELGE

C.... HELGE calculates the energy momenta for the nuclear spectator remnent
C.... and GIVES FERMI-MOTION TO THE NUCLEONS IN THE NUCLEI, here the input
C.... nucleons are assumed to be moving along the Z-axis.
C.... By giving the rest nucleons a fermi momentum, then binding energy
C.... has to be included to ensure energy-momentum conservation. This
C.... is achieved here by puting the nucleons off shell.

      PARAMETER (KSZ1=20, KSZ2=300)
      COMMON/FRINTN0/PLI0(2,4),AOP(KSZ1),IOP(KSZ1),NFR(KSZ1)
      COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)
      COMMON/FRINTN1/PPS(2,KSZ2,5),PPH(2,KSZ2,5),PPSY(2,KSZ2,5),PPA(2,5)
      COMMON/FRINTN3/IDN(2,KSZ2),FMN(2,KSZ2),NUC(2,3000)
      DIMENSION FVECT(3,KSZ2)
      SAVE /FRINTN0/,/FRPARA1/,/FRINTN1/,/FRINTN3/
      DATA PI /3.1415926/

      DO 7 L=1,2
      NWD = IOP(3+2*(L-1))-IOP(8+L)
      DO 9 LO=1,4
  9      PPA(L,LO) = FLOAT(NWD)* PLI0(L,LO)
  7      PPA(L,5)=FRSQR(PPA(L,4)*PPA(L,3),'PPA')

      IF(KFR(4).EQ.0) RETURN

C......................................................................
      DO 700 L=1, 2

      IF (IOP(3+2*(L-1)).LE.1) GOTO 700

   15    DO 18 LO=1,5
   18    PPA(L,LO) = 0.0

       DO 10 I=1,3

        SUM=0.
        DO 20 J=1,IOP(3+2*(L-1))
          SL1=PYR(0)
          SL2=PYR(0)
          SL1=FRSQR(-2.*LOG(MAX(1.E-15,SL1))*IOP(3+2*(L-1))
     >                      /(IOP(3+2*(L-1))-1), 'SL1KL7')
          SL2=COS(2.*PI*SL2)
          FVECT(I,J)=SL1*SL2*.1
   20     SUM=SUM+FVECT(I,J)

        DO 12 J=1,IOP(3+2*(L-1))
   12     FVECT(I,J)=FVECT(I,J)-SUM/FLOAT(IOP(3+2*(L-1)))

   10 CONTINUE

      DO 30 J=1,IOP(3+2*(L-1))
        SUM2=0.
        DO 32 I=1,3
   32   SUM2=SUM2+FVECT(I,J)**2
        IF (SUM2.GE..3) GOTO 15
   30 CONTINUE

      DO 50 J=1,IOP(3+2*(L-1))
        PPS(L,J,1)=FVECT(1,J)
        PPS(L,J,2)=FVECT(2,J)
        BOSFAC=FRSQR(PPS(L,J,4)/PPS(L,J,3), 'BOS590')
C.......In the rest frame, keeping E=M unchanged,
C.......nucleon off shell, FMN(L,J) changed from nucleon rest masses.
      EE = PPS(L,J,5)
        PPLUS=EE+FVECT(3,J)
        PMINUS=EE-FVECT(3,J)
        EM2=PPLUS*PMINUS-FVECT(1,J)**2-FVECT(2,J)**2

      IF(EM2.LE.0.) GOTO 15

      PPS(L,J,5) = SQRT(EM2)
      FMN(L,J) = PPS(L,J,5)
        PPS(L,J,4)=PPLUS*BOSFAC
        PPS(L,J,3)=PMINUS/BOSFAC
      IF(J.GT.IOP(8+L)) THEN
      DO 52 LO=1,4
   52 PPA(L,LO)= PPA(L,LO)+ PPS(L,J,LO)
      ENDIF
       DO 55 LO=1,4
   55  PPSY(L,J,LO) = PPS(L,J,LO)
       PPSY(L,J,5) = FRSQR(PPSY(L,J,4)*PPSY(L,J,3)-PPSY(L,J,1)**2
     >               -PPSY(L,J,2)**2, 'PPSYHG')
   50  CONTINUE

       PPA(L,5)=FRSQR(PPA(L,4)*PPA(L,3)-PPA(L,1)**2-PPA(L,2)**2,'PA1')

 700   CONTINUE

      RETURN
      END

C********************************* END FRHELGE **************************



C*************************************************************************
C      This is the routine package that sets up strings for Ariadne      *
C*************************************************************************

C********************************* FRTORST ******************************

      SUBROUTINE FRTORST(L,J)

C------------------------------------------------------------------
C Purpose: to set parton codes and momenta before entering ARIADNE and
C JETSET, and to take care of the diffractive hadrons. Fills common block
C LUJETS. J is the nucleon label in nuclus. L=1 for projectile and =2
C for the target.
C From here FRANGUR are called that handles the diffractive particles,
C FRATLEO that sets parton momenta and calls ARIADNE.
C------------------------------------------------------------------

      IMPLICIT DOUBLE PRECISION (D)
      PARAMETER (KSZ1=20,KSZ2=300)
      COMMON/FRINTN0/PLI0(2,4),AOP(KSZ1),IOP(KSZ1),NFR(KSZ1)
      COMMON/FRINTN1/PPS(2,KSZ2,5),PPH(2,KSZ2,5),PPSY(2,KSZ2,5),PPA(2,5)
      COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)
      SAVE /FRINTN0/,/FRINTN1/,/FRPARA1/

      IF(IOP(15).EQ.0.AND.PPSY(L,J,5).LT.AOP(10+L)) THEN
         CALL FRANGUR(L,J)
      ELSE
         CALL FRATLEO(L,J)
      ENDIF

      RETURN
      END

C********************************* END FRTORST **************************

C********************************* FRBELEO ******************************


      SUBROUTINE FRBELEO(IFLA,IFLB,KF)

C-----------------------------------------------------------------------
C     GIVING SPIN AND QUARKFLAVOUR TO THE ENDS OF THE EXCITED STRINGS.
C     FOR MESONS, THE ORDER OF THE END FLAVORS IS RANDOMLY GIVEN;
C     FOR BARYONS, where a quark-diquark combination, the diquark is
C     always assigned to IFLB.
C     IFLA AND IFLB ARE ADAPTED TO THE STANDARD KF CODES
C-----------------------------------------------------------------------

      PARAMETER (KSZ1=20)
      COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)
      SAVE /FRPARA1/
      INTEGER IFRKFC

      SPIN=PYR(0)
      J = ABS(KF)

      IF(J.LT.1000) THEN
C...   identify the quark and antiquark in mesons:
        J100=  J/100
        J10 = (J-J100*100)/10
        ISGN = (-1)**MAX(J100, J10)
        IF(KF.LT.0) ISGN = -ISGN
        IF(ISGN.GT.0) J10 = -J10
        IF(ISGN.LT.0) J100 = -J100

        IF(SPIN.LT..5) THEN
          IFLA=J100
          IFLB=J10
        ELSE
          IFLA=J10
          IFLB=J100
        ENDIF

      ELSEIF(J.LT.10000) THEN
        J1000=  J/1000
        J100 = (J-J1000*1000)/100
        J10 = (J-J1000*1000-J100*100)/10
        IF(KF.LT.0) THEN
        J1000=  -J1000
        J100 = -J100
        J10 = -J10
        ENDIF
        IF(SPIN.LT.VFR(13)) THEN
          IFLA=J1000
          IFLB=IFRKFC(J100,J10,0,1.)
        ELSEIF(SPIN.LT.VFR(13)+VFR(14)) THEN
          IFLA=J10
          IFLB=IFRKFC(J1000,J100,0,1.)
        ELSEIF(SPIN.LT.VFR(13)+VFR(14)+VFR(15)) THEN
          IFLA=J100
          IFLB=IFRKFC(J1000,J10,0,0.)
        ENDIF
C...Certain Lambda-like hadrons have two lightest quarks in spin-0:
          IF(ABS(J100).LT.ABS(J10)) THEN
        IF(SPIN.LT.VFR(13)) THEN
          IFLA=J1000
          IFLB=IFRKFC(J100,J10,0,0.)
        ELSEIF(SPIN.LT.VFR(13)+VFR(14)) THEN
          IFLA=J10
          IFLB=IFRKFC(J1000,J100,0,1.)
        ELSEIF(SPIN.LT.VFR(13)+VFR(14)+VFR(15)) THEN
          IFLA=J100
          IFLB=IFRKFC(J1000,J10,0,1.)
        ENDIF
          ENDIF

      ELSE

        CALL FRMGOUT(0,0,'Unrecognized particle KF code',
     >      real(KF),0.,0.,0.,0.)

      ENDIF

      RETURN
      END

C********************************* END FRBELEO **************************

C********************************* IFRKFC ******************************

C... THE KF CODE FOR A 2- OR 3-QUARK SYSTEM OF SPIN S COMPOSED BY
C... FLAVOR IA, IB, IC: (THE SYSTEM MUST BE qq or qqq, not qqbar, etc).
C... IT CORRESPONDS TO A DIQUARK SYSTEM IF IC=0.........................

      INTEGER FUNCTION IFRKFC(IA,IB,IC,S)

      IA0 = MAX( IABS(IA), MAX(IABS(IB),IABS(IC)))
      IC0 = MIN( IABS(IA), MIN(IABS(IB),IABS(IC)))
      IB0 = IABS(IA+IB+IC)-IA0-IC0
      IFRKFC = 1000*IA0 + 100*IB0 + 10*IC0 + INT(2.*(S+0.2))+ 1
      IF(IA.NE.IABS(IA).OR.IB.NE.IABS(IB)) IFRKFC = -IFRKFC
      RETURN
      END

C********************************* END IFRKFC **************************


C********************************* FRANGUR ******************************

      SUBROUTINE FRANGUR(L,J)

C-----------------------------------------------------------------------
C     ADD THE DIFFRACTIVE PARTICLES TO THE EVENT RECORD
C-----------------------------------------------------------------------

      PARAMETER (KSZJ=4000,KSZ1=20,KSZ2=300)
      COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)
      COMMON/FRINTN1/PPS(2,KSZ2,5),PPH(2,KSZ2,5),PPSY(2,KSZ2,5),PPA(2,5)
      COMMON/FRINTN3/IDN(2,KSZ2),FMN(2,KSZ2),NUC(2,3000)

      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      integer N,NPAD,K
      double precision P,V
      SAVE /PYJETS/

      SAVE /FRPARA1/,/FRINTN1/,/FRINTN3/

      N=N+1
      K(N,1)=1
      K(N,2)=IDN(L,J)
      K(N,3)=0
      K(N,4)=0
      K(N,5)=0
      P(N,1)=PPSY(L,J,1)
      P(N,2)=PPSY(L,J,2)
      P(N,3)=(PPSY(L,J,4)-PPSY(L,J,3))/2.
      P(N,4)=(PPSY(L,J,4)+PPSY(L,J,3))/2.
c [KG]:      P(N,5)=ULMASS(IDN(L,J))
      P(N,5)=ULGAMASS(L,J)      ! [KG]

      RETURN
      END

C********************************* END FRANGUR **************************

C********************************* FRATLEO ******************************

      SUBROUTINE FRATLEO(L,J)

      IMPLICIT DOUBLE PRECISION (D)
      PARAMETER (KSZJ=4000,KSZ1=20,KSZ2=300,MAXSTR=100)

      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      integer N,NPAD,K
      double precision P,V
      SAVE /PYJETS/

      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/RYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/ARDAT1/PARA(40),MSTA(40)
      COMMON/ARSTRS/IPF(MAXSTR),IPL(MAXSTR),IFLOW(MAXSTR),
     $                PT2LST,IMF,IML,IO,QDUMP,ISTRS
      logical QDUMP
      COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)
      COMMON/FRINTN0/PLI0(2,4),AOP(KSZ1),IOP(KSZ1),NFR(KSZ1)
      COMMON/FRINTN1/PPS(2,KSZ2,5),PPH(2,KSZ2,5),PPSY(2,KSZ2,5),PPA(2,5)
      COMMON/FRINTN2/NHP(2),IHQP(2,KSZ2),KHP(2,KSZ2,100,5),
     >               PHP(2,KSZ2,100,5)
      COMMON/FRINTN3/IDN(2,KSZ2),FMN(2,KSZ2),NUC(2,3000)
      COMMON/FRINTN4/KFEND(2,KSZ2,2)
      COMMON/FRJETS/NJ,KJ(100,5),PJ(100,5)

      COMMON/FRATLE1/NA1,KA1(KSZ2,5),PA1(KSZ2,5)

      DIMENSION DPV1(4),DPV2(4),PPSR(4),RFA(2,2),XPQ(-25:25)
     >        ,PGL(2)

      SAVE RFA, J1M,IFL1M,IFL2M
      SAVE /LUDAT1/,/RYPARS/,/ARDAT1/,/ARSTRS/,/FRPARA1/,
     >   /FRINTN0/,/FRINTN1/,/FRINTN2/,/FRINTN3/,/FRINTN4/,
     >   /FRJETS/,/FRATLE1/

C-----------------------------------------------------------------------
C to set parton codes and momenta for the J-th hadron string, with ends
C flavours IFL1 and IFL2.  Calls ARIADNE for dipole shower.
C For baryons, IFL2 always corresponds to the diquark end.
C This routine is entered twice (for each collision), once in test mode
C    and once for the final gluon radiation.
C
C Codes used for partons:
C   K(J,5)=111 for a hard valence quark and its accompaning gluon kink;
C         =221 for a sea quark that has converted into a gluon.
C         =222 for the soft gluon kink accompaning a hard gluon.
C         =100 for the hard gluon.
C
C L =1 FOR PROJECTILE; AND L=2 FOR TARGET.
C    IOP(15)=0 normal mode, arrange the partons and do bremsstrahlung
C    IOP(15)=1 test mode, look for valence quarks and testing if RPS is
C              drowned. Output flag:
C           =1+L: gluon on L drowned;
C PHP(L,J,K,5) will be used to memorize the fractions.
C KHP(L,J,1,5)=I   labels the valence quark
C-----------------------------------------------------------------------

      IF(L.EQ.1) J1M= 0
                 J2M= 0
      IFL1= KFEND(L,J,1)
      IFL2= KFEND(L,J,2)

        IF(IOP(15).GT.0) THEN
      IFL2A = IFL2
      IFL2B = 0
      SPIN0 = 0.
      IF(ABS(IFL2).GT.1000) THEN
      IFL2A = (IFL2/ABS(IFL2))* (ABS(IFL2)/1000)
      IFL2B = (IFL2/ABS(IFL2))* ((ABS(IFL2)-ABS(IFL2A)*1000)/100)
      SPIN0 = FLOAT((ABS(IFL2)-ABS(IFL2A)*1000-ABS(IFL2B)*100)/2)
      ENDIF
        ENDIF

      IQQK = 0
      IQGL = 0
      NA1 = 0

      IF(IHQP(L,J).LE.0) GOTO 150

      I0 = KHP(L,J,IHQP(L,J),4)+1

      DO 1000 I = I0, IHQP(L,J)

      KHP(L,J,I,4) = 0

C   Check if a quark is possibly a valence:
C   The processes in which quark flavorS changed are treated as sea quarks:
C   Note since PYTHIA only has pions and protons, only u,d can be valence.

      IF(IOP(15).GT.0.AND.NJ.GT.0
     >   .AND.KHP(L,J,1,5).eq.0.AND.KHP(L,J,I,2).NE.21
     >   .AND.(MSTI(1).NE.12.AND.MSTI(1).NE.53) ) THEN

      IU = KHP(L,J,I,2)
      IUSN = (IU/IABS(IU))

        IF(IU.EQ.IFL1.OR.IUSN*(IABS(IU)+2).EQ.IFL1.OR.
     >     IU.EQ.IFL2A.OR.IUSN*(IABS(IU)+2).EQ.IFL2A.OR.
     >     IU.EQ.IFL2B.OR.IUSN*(IABS(IU)+2).EQ.IFL2B) THEN
C                          weighted by the structure functions
        CALL RYSTFU(MSTI(11+L-1),PARI(33+L-1),PARI(21)**2,XPQ)

        RVAL = ABS(XPQ(IU)-XPQ(-IU))/MAX(XPQ(IU),XPQ(-IU))

          IF(PYR(0).LE.RVAL) THEN
          KHP(L,J,1,5) = I
           IF(IU.EQ.IFL1.OR.IUSN*(IABS(IU)+2).EQ.IFL1) THEN
             KHP(L,J,I,2) = IFL1
           ELSEIF(IU.EQ.IFL2A.OR.IUSN*(IABS(IU)+2).EQ.IFL2A) THEN
             KHP(L,J,I,2) = IFL2A
             IF(IFL2B.EQ.0) THEN
               IFL2 = IFL1
               IFL1 = IFL2A
             ELSE
               SPIN = SPIN0
               IF(IFL1.EQ.IFL2B) SPIN = 1.
               IFL2 = IFRKFC(IFL1,IFL2B,0,SPIN)
               IFL1 = IFL2A
             ENDIF
           ELSEIF(IU.EQ.IFL2B.OR.IUSN*(IABS(IU)+2).EQ.IFL2B) THEN
             KHP(L,J,I,2) = IFL2B
             SPIN = SPIN0
             IF(IFL1.EQ.IFL2A) SPIN = 1.
             IFL2 = IFRKFC(IFL1,IFL2A,0,SPIN)
             IFL1 = IFL2B
           ELSE
             KHP(L,J,1,5) = 0
           ENDIF
          ENDIF
        ENDIF
      ENDIF


      IF(I.EQ.KHP(L,J,1,5).AND.IQQK.EQ.0) THEN
      DO 131 LO=1,4
      KA1(1,LO) = KHP(L,J,I,LO)
131   PA1(1,LO) = PHP(L,J,I,LO)
      PA1(1,5) = 0.0
        KA1(1,1)= 2
        KA1(1,3)= 0
        KA1(1,4)= 0
        KA1(1,5)= 111
        IQQK=1

C..       Save the end configurations:
        IF(IOP(15).GT.0) THEN
          IF(L.EQ.1) THEN
          J1M= J
          IFL1M= IFL1
          IFL2M= IFL2
          ELSE
          J2M= J
          ENDIF
        ENDIF

      ELSE
C         GLUON entry started at NA1=2
      NA1 = MAX(NA1+1,2)
      DO 138 LO=1,4
      PA1(NA1,LO) = PHP(L,J,I,LO)
138   KA1(NA1,LO) = KHP(L,J,I,LO)
      PA1(NA1,5) = 0.0
      KA1(NA1,1)=2
      KA1(NA1,5)=0
      IF(KA1(NA1,2).NE.21) KA1(NA1,5) = 221
      KA1(NA1,2)=21
      KA1(NA1,3)=0
      KA1(NA1,4)=0
      KA1(NA1,5)=0
      IQGL= IQGL + 1

      ENDIF

1000  CONTINUE

C...ORDER THE GLUONS ACCORDING TO PT..............
      IF(NA1.GT.2) CALL FRORDER(L,2,NA1)

C...DIPOLE RADIATION.......................

150   IOP(17) = N+1

      IF(IQQK.EQ.1) THEN
      N = N+1
      DO 190 LO=1,5
      K(N,LO) = KA1(1,LO)
190   P(N,LO) = PA1(1,LO)
      P(N,5) = 0.0
      P(N,4) = SQRT(P(N,1)**2+P(N,2)**2+P(N,3)**2)
      ENDIF

C...............Set up the two fractions:

      RFA(L,1) = 1.0
      IF(KFR(10).NE.0.AND.IQGL.GT.0) THEN
        IF(PHP(L,J,1,5).GT.0) THEN
        RFA(L,1) = PHP(L,J,1,5)
        ELSE
          IF(KFR(10).EQ.1) THEN
          RFA(L,1) = VFR(16)
          ELSEIF(KFR(10).EQ.2) THEN
          RFA(L,1) = PYR(0)
          ENDIF
        PHP(L,J,1,5)= RFA(L,1)
        ENDIF
      ENDIF
      IF(RFA(L,1).GT.1.0) CALL FRMGOUT(0,0,'VFR(16)>1 NOT ALLOWED!',
     >   VFR(16),RFA(L,1),0.,0.,0.)

      RFA(L,2) = 1.0
      IF(KFR(8).GE.1.AND.IQGL.GT.0) THEN
        IF(PHP(L,J,2,5).GT.0) THEN
        RFA(L,2) = PHP(L,J,2,5)
        ELSE
      PGL(1) = PA1(IQGL+1,4)-PA1(IQGL+1,3)
      PGL(2) = PA1(IQGL+1,4)+PA1(IQGL+1,3)
      AP = SQRT( PGL(L)/PPSY(L,J,2+L))
C       Rarely PGL can become zero due to inaccuracy when Pz is very large.
      RKK = 0.
      RKKMX= 1.0-(PPS(L,J,1)**2+PPS(L,J,2)**2)/(PPS(L,J,3)*PPS(L,J,4))
      RKKMX = 0.99*RKKMX
593   IF(AP.GT.0.) RKK=AP*((RKKMX+AP)/AP)**PYR(0)-AP
      IF(RKK.LE.0..OR.RKK.GE.RKKMX) GOTO 593
      RFA(L,2) = 1.- RKK
      PHP(L,J,2,5)= RFA(L,2)
        ENDIF
      ENDIF

      DO 588 LO=1,4
588   PPSR(LO) = PPS(L,J,LO)

C...for the kink:
      IF(RFA(L,2).GT.0.and.RFA(L,2).LT.1.0) THEN
      PGL(L) = (1.-RFA(L,2))* PPS(L,J,2+L)
      PA1(1,1)= 0.0
      PA1(1,2)= 0.0
      PGL(3-L) = 0.0
      PA1(1,3)=0.5*(PGL(2)-PGL(1))
      PA1(1,4)=0.5*(PGL(2)+PGL(1))
      KA1(1,1)=2
      KA1(1,2)=21
      KA1(1,4)=3
      KA1(1,5)=222
      PPSR(3) = PPSR(3) - PGL(1)
      PPSR(4) = PPSR(4) - PGL(2)
      ENDIF

590   CALL FRPPART(L,PPSR,DPV1,DPV2)

      N = N+1
      DO 599 LO=1,4
599   P(N,LO) = DPV1(LO)
      P(N,5) = 0.0
      P(N,4) = SQRT(P(N,1)**2+P(N,2)**2+P(N,3)**2)
      K(N,1) = 2
      K(N,3) = 0
      K(N,4) = 1
      IF(IQQK.EQ.0) THEN
      K(N,2) = IFL1
      K(N,5) = 0
      ELSEIF(IQQK.EQ.1) THEN
      K(N,2) = 21
      K(N,5) = 111
      ENDIF

      N = N+1
      DO 600 LO=1,4
600   P(N,LO) = DPV2(LO)
      P(N,5) = 0.0
      P(N,4) = SQRT(P(N,1)**2+P(N,2)**2+P(N,3)**2)
      K(N,1) = 1
      K(N,2) = IFL2
      K(N,3) = 0
      K(N,4) = 2
      K(N,5) = 0

C...  Check mass and remove negative mass arising from numerical
C...  imprecisions.
      RMS0= ppsr(3)*ppsr(4)- ppsr(1)**2-ppsr(2)**2
      JN=N
      IF(P(N-1,4).LT.P(N,4)) JN=N-1
      PP2=(P(N,3)+P(N-1,3))**2+(P(N,2)+P(N-1,2))**2
     >     +(P(N,1)+P(N-1,1))**2
      XRMS20=(P(N,4)+P(N-1,4))**2-PP2 -RMS0

      IF(XRMS20.LT.0.) THEN
      XADD0=0.0
      XADD1=0.0
700   XADD1=XADD1+0.1
      PJNV= P(JN,4)+XADD1
      XRMS2=(P(2*N-1-JN,4)+PJNV)**2-PP2-RMS0
      IF(XRMS2*XRMS20.GT.0.) GOTO 700

      NTRY=0
710   XADD= (XADD0+XADD1)/2.0
       XRMS2M=XRMS2
       PJNV= P(JN,4)+XADD
       XRMS2=(P(2*N-1-JN,4)+PJNV)**2-PP2-RMS0
       IF(XRMS2*XRMS20.GT.0.) THEN
        XADD0=XADD
       ELSE
        XADD1=XADD
       ENDIF
       IF(XRMS2.EQ.XRMS2M) THEN
         NTRY=NTRY+1
       ELSE
         NTRY=0
       ENDIF

       IF(NTRY.GE.5) THEN
         IF(XRMS2.LE.-0.5*RMS0) XADD=XADD1
         GOTO 720
       ENDIF

       GOTO 710

720    P(JN,4) = P(JN,4)+ XADD
      ENDIF

      MSTA(11) = 0

C..Test hard partons against bremsstrahlung
      IF(KFR(9).NE.0.AND.IOP(15).EQ.1.AND.IQGL+IQQK.GT.0) THEN
        IF(IQGL.GT.0.AND.RFA(L,1).LT.1.0) THEN
        PARA(11) = VFR(7+L)/(RFA(L,1))
        PARA(12) = VFR(7+L)/(1.-RFA(L,1))
        ENDIF
        CALL FRTESTG(L,IQQK,IQGL,IOK,RFA)
        IF(IOK.EQ.0) THEN
        IOP(15) = IOP(15) +L
          IF(J1M.GT.0) KHP(1,J1M,1,5)= 0
          IF(J2M.GT.0) KHP(2,J2M,1,5)= 0
        GOTO 999
        ELSEIF(L.EQ.2.AND.IQQK+J1M.GT.0) THEN
C..            If valence quarks survived, keep the new ends:
        KFEND(L,J,1) = IFL1
        KFEND(L,J,2) = IFL2
          IF(J1M.GT.0) THEN
        KFEND(1,J1M,1) = IFL1M
        KFEND(1,J1M,2) = IFL2M
          ENDIF
        ENDIF
      ENDIF

C.. Insert hard gluons one by one, and do emission.  The Pt of emission is
C.. restricted such that   Pt_next gluon < Pt < Previous emmision.

                         IF(IOP(15).EQ.0) THEN
      PARA(11) = VFR(7+L)
      PARA(12) = PARA(11)
      MSTA(11) = 0
      IF(RFA(L,1).LT.1.0) THEN
       PARA(11) = VFR(7+L)/(RFA(L,1))
       PARA(12) = VFR(7+L)/(1.-RFA(L,1))
      ENDIF

      IARI = 0
      PARA3=PARA(3)
      PARA6=PARA(6)

      DO 900 I = 1, MAX(IQGL,1)

        INUP= 0
        IF(I.EQ.1) INUP= 0
        CALL FRMXGPT(IOP(17),INUP,IMX,VRPTNX,1)
        IF(IMX.GT.0) THEN
          CALL FRINSET(IMX,IOP(17),N,NOG,1)
          IF(I.EQ.1.AND.KFR(8).GT.0) THEN
            CALL FRINSKK(RFA(L,2),NKK)
            IF(NKK.GT.0) THEN
            PARA(13) = PARA(11)/(1.-RFA(L,2))
            PARA(11) = PARA(11)/RFA(L,2)
            ENDIF
          ENDIF
        ENDIF

        IF(KFR(2).EQ.1) THEN
        PARA(3) = MAX(PARA3,VRPTNX)
        IF(IARI.GT.0) PARA(6) = FRSQR(PT2LST,'PT2LSTAR')
        NMEM = N

        CALL FRARIAD
        IARI= N-NMEM
        ENDIF

900   CONTINUE
      PARA(3)=PARA3
      PARA(6)=PARA6
                                 ENDIF

999   RETURN
      END

C************************************************************************

      REAL FUNCTION FRIPT(I,N1,N2,IQ)

C....Evaluate the invariant P_T**2 of parton I as if I is put between
C....partons N1 and N2.
C....IQ=1: Definition 1: s12*s13/s123
C....  =2: Definition 2: s12*s13/(s123-s12-s13) (true P_T^2 in CMS of N1 N2)
C....  <0: use the N1 and N2 of previously memorized (N1,N2 dummy here).

      IMPLICIT DOUBLE PRECISION (D)
      PARAMETER (KSZJ=4000)

      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      integer N,NPAD,K
      double precision P,V
      SAVE /PYJETS/


      DIMENSION DM(2,5), DI(4), DS(3)
      SAVE DM

      IF(IQ.GT.0) THEN
      DO 20 LO=1,5
      DM(1,LO) = DBLE( P(N1,LO))
20    DM(2,LO) = DBLE( P(N2,LO))
            IF(I.EQ.N1.OR.I.EQ.N2) THEN
      FRIPT = 0.
      RETURN
            ENDIF
      ENDIF

      DO 30 L=1,3
      DO 35 LO=1,4
      IF(L.LE.2) DI(LO) = DM(L,LO)+ DBLE( P(I,LO))
      IF(L.EQ.3) DI(LO) = DM(1,LO)+ DI(LO)
35    CONTINUE
30    DS(L) = DI(4)**2- DI(3)**2-DI(2)**2-DI(1)**2

      DS(1) = DS(1) - (DM(1,5)+DBLE(P(I,5)))**2
      DS(2) = DS(2) - (DM(2,5)+DBLE(P(I,5)))**2
      IF(IABS(IQ).EQ.2) DS(3) = DS(3) - DS(1) - DS(2)
      IF(DABS(DS(3)).LE.1.D-5) DS(3) = 1.D-5

      FRIPT = SNGL (DS(1)*DS(2)/DS(3) )

      RETURN
      END

C************************************************************************
C********************************* END FRATLEO **************************

C*********************** SUBROUTINE FRTESTG *****************************

      SUBROUTINE FRTESTG(L,IQQK,IQGL,IOK,RFA)

C.......Test a gluon against the bremsstrahlung background
C....... IQGL = the number of gluons
C....... IOK =1 - gluon sticks out ;  IOK=0 - drawned.

      PARAMETER (KSZJ=4000,KSZ1=20,KSZ2=300)
      COMMON/FRINTN1/PPS(2,KSZ2,5),PPH(2,KSZ2,5),PPSY(2,KSZ2,5),PPA(2,5)
      COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)
      COMMON/FRINTN0/PLI0(2,4),AOP(KSZ1),IOP(KSZ1),NFR(KSZ1)

      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      integer N,NPAD,K
      double precision P,V
      SAVE /PYJETS/

      COMMON/ARDAT1/PARA(40),MSTA(40)
      SAVE /FRINTN1/,/FRPARA1/,/FRINTN0/,/ARDAT1/
      DIMENSION RFA(2,2)


      IOK = 0

      MSTA6 = MSTA(6)
      MSTA(6) = 1
      PARA11 = PARA(11)
      PARA3 = PARA(3)

      CALL FRSAVEN(IOP(17),0)

      NOG = 0
      PTGL2 = 0.
      IF(IQQK.GT.0) PTGL2=MAX(PTGL2,FRIPT(IOP(17),IOP(17)+1,
     >                        IOP(17)+2,1))
       IF(IQGL.GT.0) THEN
       CALL FRMXGPT(IOP(17),N,IMX,VRPTNX,0)
       CALL FRINSET(IMX,IOP(17),N,NOG,1)
        CALL FRINSKK(RFA(L,2),NKK)
        IF(NKK.GT.0) THEN
        PARA(13) = PARA(11)/(1.-RFA(L,2))
        PARA(11) = PARA(11)/RFA(L,2)
        ENDIF
       IF(NOG.GT.0) PTGL2 = MAX(PTGL2,FRIPT(NOG,NOG-1,NOG+1,1) )
       ENDIF

      IF(PTGL2.GE.PARA3**2) THEN
C        PARA(3) = SQRT(PTGL2)
C        NM = N
C        CALL FRARIAD
C        IF(N-NM.EQ.0) IOK= 1
C
      PTAR = 0.
      NM = N
      CALL FRARIAD
       IF(N -NM.GE.1) THEN
        DO 11 II=IOP(17)+1, N-1
          IF(K(II,2).EQ.21.AND.K(II,5).EQ.1) THEN
          PTAR = FRIPT(II,II-1,II+1,1)
          GOTO 100
          ENDIF
 11     CONTINUE
       ENDIF
100   IF(PTAR.LE.PTGL2) IOK= 1

      ENDIF

C.......Restore the configuration:

      CALL FRSAVEN(IOP(17),1)
      MSTA(6) = MSTA6
      PARA(3) = PARA3
      PARA(11) = PARA11

      RETURN
      END


C************************END FRTESTG ************************************

C*********************** SUBROUTINE FRARIAD ***********************

      SUBROUTINE FRARIAD

C..Fritiof interface to Ariadne_4.02r.  LUJETS entries from IOP(17) to N
C..are copied to Ariadne event record ARJETX, and after emission is done
C..partons are copied back onto LUJETS.


      PARAMETER (KSZJ=4000,KSZ1=20)
      COMMON/FRINTN0/PLI0(2,4),AOP(KSZ1),IOP(KSZ1),NFR(KSZ1)

      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      integer N,NPAD,K
      double precision P,V
      SAVE /PYJETS/

      COMMON/ARJETX/NO,KO(300,5),PO(300,5),VO(300,5)
      SAVE /FRINTN0/,/ARJETX/

      NO=0
      DO 100 I=IOP(17),N
      NO=NO+1
      DO 100 LO=1,5
      KO(NO,LO) = K(I,LO)
100   PO(NO,LO) = P(I,LO)

      CALL AREXEC

      N=IOP(17)-1
      DO 200 IO=1,NO
      IF(KO(IO,1).GE.11) GOTO 200
      N=N+1
      DO 250 LO=1,5
      K(N,LO) = KO(IO,LO)
250   P(N,LO) = PO(IO,LO)
200   CONTINUE

      RETURN
      END

C*********************** END FRARIAD *****************************

C******************************** FRMXGPT ********************************

      SUBROUTINE FRMXGPT(N1,N2,IMX,VRPTNX,IDROP)

C..   To find the gluon with maximu inv-pt among those on FRATLE1:2-NA1.
C..   IDROP=1: The gluon will have KA1(I,3)=-2001 once it has been used.
C..        =0: Gluons will not be marked (used when called from FRTESTG.
C..   The inv pt is calculated assuming the gluon is between N1 and N2.
C..   IMX=0 - no more gluons are available;
C..   IMX>=2 - index of the gluon with largest inv-pt.
C..   VRPTNX gives value of the NEXT largest pt.
C..   When N1 or N2 <= 0, previous memorized p(N1,)and p(N2,) will be used.

      PARAMETER (KSZJ=4000,KSZ2=300)

      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      integer N,NPAD,K
      double precision P,V
      SAVE /PYJETS/

      COMMON/FRATLE1/NA1,KA1(KSZ2,5),PA1(KSZ2,5)
      SAVE /FRATLE1/

      VRPTNX=0.
      IMX=0

      IF(NA1.EQ.2) THEN
      IMX=2

      ELSEIF(NA1.GT.2) THEN

      IQ=1
      IF(N1.LE.0.OR.N2.LE.0) IQ=-1
      VRPT2M=-1.E4
      DO 100 I=2, NA1
      IF(KA1(I,3).EQ.-2001) GOTO 100
        DO 120 LO=1,5
120     P(N+1,LO)=PA1(I,LO)
      VRPT2 = FRIPT(N+1,N1,N2,IQ)
      IQ=-1

      IF(VRPT2.GT.VRPT2M.OR.IMX.EQ.0) THEN
      IMX = I
      VRPTNX=SQRT(MAX(VRPT2M,0.))
      VRPT2M=VRPT2
      ENDIF

100   CONTINUE

      ENDIF

      IF(IDROP.GT.0.AND.IMX.GT.0) KA1(IMX,3)=-2001

      RETURN
      END

C******************************** END FRMXGPT ********************************

C******************************** FRINSET ********************************

      SUBROUTINE FRINSET(NA,N1,N2,NOG,IQ)

C.......To place a gluon specified by NA on FRATLE1 block onto LUJETS
C.......between N1,N2, where N1 is assumed to be quark N2G could be
C.......the gluon kink or the end diquark.
C.......The placing is based on rapidity
C.......ordering, and if more than one rapidity-ordered spots are
C.......found then the one giving maximum invariant mass is used.
C.......If no such place is found then NA is placed near N1.
C....... IQ:=0 the actual insertion will not take place, N1,N2,N unchanged;
C........   =1 the insertion takes place, and the entire LUJETS N>NOG
C........      is shifted by 1.
C....... output NOG: the actual place NA was placed;
C....... The inserted gluon has a code: K(NOG,5)=100.
C....... The gluon will not be placed at the string ends.  So in case of
C....... q-qbar pairs from gluon splitting in Ariadne, the gluon
C....... will not be placed between the pair.


      PARAMETER (KSZJ=4000,KSZ2=300)

      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      integer N,NPAD,K
      double precision P,V
      SAVE /PYJETS/

      COMMON/FRATLE1/NA1,KA1(KSZ2,5),PA1(KSZ2,5)
      DIMENSION Y(2)
      SAVE /FRATLE1/

      NOG = 0
      IF(N2.LE.N1+1) GOTO 250

      YGL = PA1(NA,4) - PA1(NA,3)
      YGU = PA1(NA,4) + PA1(NA,3)
      IF(YGU.LE.0.) THEN
      YG = -.9E10
      ELSEIF(YGL.LE.0.) THEN
      YG = +.9E10
      ELSE
      YG = .5* LOG(YGU/YGL)
      ENDIF

      VM2X = 0.
      DO 200 I = N1, N2-1
      IF(K(I,5).EQ.222.AND.K(I,2).EQ.21) GOTO 250
      NUM=2
      IF(I.NE.N1) THEN
      NUM=1
      Y(1) = Y(2)
      ENDIF
       DO 100 II=1,NUM
       IR = I+2-II
       YL = P(IR,4) - P(IR,3)
       YU = P(IR,4) + P(IR,3)
CC       IF(YU.LE.0.AND.YL.LE.0.) CALL FRMGOUT(0,1,' 0 MOMENTA!',
CC     >         FLOAT(IR),P(IR,1),P(IR,2),P(IR,3),P(IR,4))
       IF(YU.LE.0.) THEN
       YR = -1.E10
       ELSEIF(YL.LE.0.) THEN
       YR = +1.E10
       ELSE
       YR = .5* LOG(YU/YL)
       ENDIF
       Y(3-II) = YR
100    CONTINUE

      IF( K(I,1).EQ.2 .AND.
     >  (YG.GE.MIN(Y(1),Y(2)).AND.YG.LE.MAX(Y(1),Y(2))) ) THEN
      VM2=(PA1(NA,4)+P(I,4)+P(I+1,4))**2-(PA1(NA,3)
     >     +P(I,3)+P(I+1,3))**2-(PA1(NA,2)+P(I,2)+P(I+1,2))**2-
     >     (PA1(NA,1)+P(I,1)+P(I+1,1))**2
       IF(VM2.GT.VM2X) THEN
       VM2X = VM2
       NOG = I+1
       ENDIF
      ENDIF
200   CONTINUE

250    NOG = MAX(N1+1, NOG)

      IF (IQ.GE.1) THEN
      DO 150 I = N, NOG, -1
      DO 150 LO = 1,5
      K(I+1,LO) = K(I,LO)
150   P(I+1,LO) = P(I,LO)

      DO 160 LO = 1,5
      K(NOG,LO) = KA1(NA,LO)
160   P(NOG,LO) = PA1(NA,LO)
      P(NOG,4) = SQRT(P(NOG,1)**2+P(NOG,2)**2+P(NOG,3)**2)
      P(NOG,5) = 0.0
      K(NOG,3) = 0
      K(NOG,4) = 0
      K(NOG,5) = 100
      N = N+1
      ENDIF

      RETURN
      END

C******************************** END FRINSET ********************************

C****************************** FRPPART *********************************

      SUBROUTINE FRPPART(L,PPSR,DPV1,DPV2)

C........If a system has a lightcone-momenta PPSR(4), it is partitioned
C........into two momenta corresponding to the quark (DPV1) and diquark(DPV2)
C........end.  Note the input PPSR is light-cone: Px,Py,P-,P+, and the output
C........DPVs are normal 4-vector: Px,Py,Pz,E.

      PARAMETER (KSZ1=20,PI = 3.1415926)
      IMPLICIT DOUBLE PRECISION (D)
      DIMENSION PPSR(4),DPV1(4),DPV2(4)
      COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)
      SAVE /FRPARA1/

      PPST2 = PPSR(1)**2+PPSR(2)**2
      WT2 = PPSR(3)*PPSR(4)
      SM = WT2 - PPST2
      GP2MX = 0.25* (SQRT(WT2)-SQRT(PPST2))**2

      NTRY= 0
10    NTRY=NTRY+1
      IF(NTRY.GT.200) CALL FRMGOUT(0,1,'NTRY runaway loop:',
     >      FLOAT(NTRY),0.,0.,0.,0.)

C.... ....TO GENERATE A GAUSSIAN PT FOR THE DIQUARK ..............

      CALL FRGAUSS(GP2, VFR(7), GP2MX)
        Gpt = FRSQR(GP2, 'PT2HGF')
        phi= 2.*PI*PYR(0)
        DPV2(1) = DBLE(Gpt*cos(phi))
        DPV2(2) = DBLE(Gpt*sin(phi))

      DPV1(1) = -DPV2(1) + DBLE(PPSR(1))
      DPV1(2) = -DPV2(2) + DBLE(PPSR(2))

C........DPV1 CORRES. TO THE QUARK END, DPV2 CORRES. TO THE DIQUARK END....

      DGPT1 = DPV1(1)**2 + DPV1(2)**2
      DGPT2 = DPV2(1)**2 + DPV2(2)**2
      DTM1 =DBLE(PPSR(3))+ (DGPT2-DGPT1)/DBLE(PPSR(4))
      DTM2=DTM1**2-4.D0*DGPT2*DBLE(PPSR(3))/DBLE(PPSR(4))
      IF(DTM2.LT.-0.1) CALL FRMGOUT(0,1,'CHECK DTM2',WT2,PPST2,
     >    REAL(DGPT1),REAL(DGPT2),GP2MX)
      DTM2=DFRSQR(DMAX1(0.D0,DTM2),'DTM2$')

      DGPV2M = (DTM1-(-1)**(L-1)*DTM2)/2.D0
      DGPV1M = DBLE(PPSR(3))-DGPV2M
      DGPV1P = DGPT1/DGPV1M
      DGPV2P = DBLE(PPSR(4))-DGPV1P
CC      DGPV1P = (DBLE(PPSR(4))*DGPV2M-DGPT2+DGPT1)/DBLE(PPSR(3))

      DPV1(3) = 0.5D0*(DGPV1P-DGPV1M)
      DPV1(4) = 0.5D0*(DGPV1P+DGPV1M)
      DPV2(3) = 0.5D0*(DGPV2P-DGPV2M)
      DPV2(4) = 0.5D0*(DGPV2P+DGPV2M)

      RETURN
      END

C****************************** END FRPPART *****************************

C************************************ FRINSKK ***************************

      SUBROUTINE FRINSKK(XF,NKK)

C...To insert the soft gluon kink stored at NA1=1 to LUJETS_N.
C...OUTPUT: NkK=location of the kink.
C.......   NKK=0 - gluon kink is not inserted;

      PARAMETER (KSZJ=4000,KSZ1=20,KSZ2=300)
      COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)
      COMMON/FRINTN1/PPS(2,KSZ2,5),PPH(2,KSZ2,5),PPSY(2,KSZ2,5),PPA(2,5)
      COMMON/FRATLE1/NA1,KA1(KSZ2,5),PA1(KSZ2,5)

      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      integer N,NPAD,K
      double precision P,V
      SAVE /PYJETS/

      SAVE /FRPARA1/,/FRINTN1/,/FRATLE1/

      IF(XF.LE.1E-5.OR.XF.GE.0.99999.or.KA1(1,2).NE.21) THEN
      NKK = 0
      RETURN
      ENDIF

      NKK = N
        DO 200 LO = 1, 5
        P(N+1,LO) = P(N,LO)
200     K(N+1,LO) = K(N,LO)

        DO 205 LO = 1, 5
        K(NKK,LO) = KA1(1,LO)
205     P(NKK,LO) = PA1(1,LO)
        P(NKK,4) = SQRT(P(NKK,1)**2+P(NKK,2)**2+P(NKK,3)**2)
        P(NKK,5) = 0.0

        N = N+1

      RETURN
      END

C********************************* END FRINSKK ***************************

C*************************** FRSAVEN ********************************


      SUBROUTINE FRSAVEN(N1,IQ)

C..IQ=0: To save a lujets configuration (with 3 partons) N1 and N temporarily
C..IQ=\=0: restore the configuration at N1 and N1+1=N

      PARAMETER (KSZJ=4000)

      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      integer N,NPAD,K
      double precision P,V
      SAVE /PYJETS/

      DIMENSION KM(3,5),PM(3,5)
      SAVE KM, PM, NUM

      IF(IQ.EQ.0) THEN
      NUM = N-N1+1
      IF(NUM.GT.3) CALL FRMGOUT(0,0,'More than 3 partons in
     > FRSAVEN',float(n1),float(N),0.,0.,0.)

      DO 100 I=N1,N
      IM = I-N1+1
      DO 100 LO=1,5
      PM(IM,LO) = P(I,LO)
100   KM(IM,LO) = K(I,LO)

      ELSE

      DO 200 I=1,NUM
      NI = N1+I-1
      DO 200 LO=1,5
      P(NI,LO) = PM(I,LO)
200   K(NI,LO) = KM(I,LO)

      N = N1+NUM-1

      ENDIF

      RETURN
      END

C*************************** END FRSAVEN ********************************

C********************************* FRFILHW ******************************

      SUBROUTINE FRFILHW

C     TO ADD THOSE VECTOR BOSONS,HIGGS ETC (if they are produced from
C     parton subprocesses) TO LUJETS

      PARAMETER (KSZJ=4000)

      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      integer N,NPAD,K
      double precision P,V
      SAVE /PYJETS/

      COMMON/FRCNUT/NR,KR(10,5),PR(10,5),NR0
      SAVE /FRCNUT/

      IF(NR.LE.0) RETURN
      DO 10 I = 1, NR
       N= N+1
       DO 20 J=1, 5
       K(N, J) = KR(I,J)
20     P(N, J) = PR(I,J)
        K(N,1)=1
        K(N,3)=0
10    CONTINUE

      RETURN
      END

C********************************* END FRFILHW **************************

C********************************* FRORDER ******************************

      SUBROUTINE FRORDER(L,NS,NE)

C......TO ORDER PARTICLES (gluons) ACCORDING TO
C........FOR KFR12=1, ASCENDING RAPIDITY FOR PROJECTILE
C........              DESCENDING RAPIDITY FOR TARGET
C........FOR KFR12>,=2, ASCENDING PT

      PARAMETER (KSZ2=300,KSZ1=20)
      COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)
      COMMON/FRATLE1/NA1,KA1(KSZ2,5),PA1(KSZ2,5)
      COMMON/FRDUMYOR/KO(KSZ2,5),PO(KSZ2,5),Y(KSZ2),II(KSZ2)
      SAVE /FRPARA1/,/FRATLE1/,/FRDUMYOR/

      KFR12 = 2
      IF(NE.LE.NS) RETURN
      IF(NE-NS+1.GT.KSZ2) CALL FRMGOUT(0,1,'FRDUMYOR array size
     > insufficient', FLOAT(NS),FLOAT(NE),0.,0.,0.)

      SML = 1.E-20
      if(L.EQ.1.OR.KFR12.ge.2) THEN
C             !ascending
      IQ = -1
      ELSE
C              !descending
      IQ = 1
      ENDIF

      IOR=0
      DO 31 IO = NS, NE
      IR = IO - NS + 1
      DO 33 J=1, 5
      KO(IR,J) = KA1(IO,J)
33    PO(IR,J) = PA1(IO,J)
      PPLS = PA1(IO,4) + PA1(IO,3)
      PMIS = PA1(IO,4) - PA1(IO,3)
      IF(KFR12.EQ.1) THEN
      Y(IR) = .5* LOG( MAX(PPLS,SML)/MAX(PMIS,SML))
      ELSEIF(KFR12.GE.2) THEN
      Y(IR) = PA1(IO,1)**2+PA1(IO,2)**2
      ENDIF
      IF(IOR.EQ.0.AND.IR.GE.2) THEN
       IF(IQ.EQ.1.AND.Y(IR).GT.YLIM) IOR=1
       IF(IQ.EQ.-1.AND.Y(IR).LT.YLIM) IOR=1
      ENDIF
      IF(IR.EQ.1) THEN
      YLIM = Y(1)
      ELSEIF(IR.GT.1.AND.IQ.EQ.1) THEN
      YLIM = MIN(YLIM,Y(IR))
      ELSEIF(IQ.EQ.-1) THEN
      YLIM = MAX(YLIM,Y(IR))
      ENDIF

31    CONTINUE

      IF(IOR.EQ.0) RETURN

      CALL FRORD01(Y, II, NE-NS+1, IQ)

      DO 35 IO = NS,NE
      IR = IO - NS + 1
       DO 35 J = 1, 5
      KA1(IO,J) = KO(II(IR),J)
35    PA1(IO,J) = PO(II(IR),J)

      RETURN
      END

C************************************ FRORD01 **************************

      SUBROUTINE FRORD01(P,II,N,IQ)

C           Routine to arrange II so that
C           P(ii(1)) >= P(ii(2)) >= ... >= P(ii(n)), IF IQ>,=0
C           P(ii(1)) <= P(ii(2)) <= ... <= P(ii(n)), IF IQ<0
C
      real p(1)
      integer ii(1)
      logical done

      do 101 k=1,n
101   ii(k)=k

      do 110 nlim = n-1,1,-1
            done = .true.
              IF(IQ.GE.0) THEN
            do 120 k = 1,nlim
                  if ( p(ii(k)) .lt. p(ii(k+1)) ) then
                        done=.false.
                        itemp=ii(k)
                        ii(k)=ii(k+1)
                        ii(k+1)=itemp
                        end if
 120        continue
            if (done) return
              ELSE
            do 130 k = 1,nlim
                  if ( p(ii(k)) .gt. p(ii(k+1)) ) then
                        done=.false.
                        itemp=ii(k)
                        ii(k)=ii(k+1)
                        ii(k+1)=itemp
                        end if
 130        continue
            if (done) return
             ENDIF
110   continue

      RETURN
      END


C********************************* END FRORD01 **************************



C*************************************************************************
C**                                                                     **
C**   This package interfaces with RYTHIA and handles the generated     **
C**   hard partons                                                      **
C**                                                                     **
C*************************************************************************

C********************************* FRQPROB ******************************

      SUBROUTINE FRQPROB(KFI,KFT,IQ)

C.... TO ESTIMATE cross sections.
C......KFI, KFT - The KF codes of the incident and target particle (nucleon).
C......IQ=0 will suppress all the write out.


      PARAMETER (KSZJ=4000,KSZ1=20)
      COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)
      COMMON/FRCODES/IPT(2),PACD(27),NNUC(27),NPROT(27),KCD(27)
     >           ,RO1(27,2),EXMA(9,2)
      COMMON/FRINTN0/PLI0(2,4),AOP(KSZ1),IOP(KSZ1),NFR(KSZ1)

      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      integer N,NPAD,K
      double precision P,V
      SAVE /PYJETS/

      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/RYSUBS/MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200)
      COMMON/RYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/RYINT1/MINT(400),VINT(400)
      COMMON/RYINT5/NGEN(0:200,3),XSEC(0:200,3)
      CHARACTER PACD*4,NAM1*16,NAM2*16,NAM1S*8,NAM2S*8,WORDS*21
      SAVE /FRPARA1/,/FRCODES/,/FRINTN0/,/LUDAT1/,/RYSUBS/,
     >     /RYPARS/,/RYINT1/,/RYINT5/

      N=0
      DO 100 L = 1,2
       DO 120 J = 1, 2
120    P(L,J) = PLI0(L,J)
      P(L,3) = (PLI0(L,4)-PLI0(L,3))/2.
100   P(L,4) = (PLI0(L,4)+PLI0(L,3))/2.

!      write(*,*) 'FRITIOF: KFR(7),PARP(2) = ',KFR(7),PARP(2) ! [KG]

      W = AOP(1)
c.......in rythia: parp(2)=10
      IF(W.LE.PARP(2)) THEN
!        IF(KFR(7).GT.0) WRITE(MSTU(11),3000) W
        IF(KFR(7).GT.0.and.IQ.gt.0) WRITE(MSTU(11),3000) W ! [KG]
        KFR(7)= 0
        PARP(2)= 0.9*W
      ENDIF
      IOP(18)= KFR(7)
      CALL LUNAME(KFI,NAM1)
      CALL LUNAME(KFT,NAM2)
      NAM1S = NAM1
      NAM2S = NAM2

C.....................................
      CALL FRSETRY(1)
c.........jochen:
      if(kfr(7).ne.0) CALL FRHARDP(KFI, KFT, W, IHAV,-1)
c      CALL FRHARDP(KFI, KFT, W, IHAV,-1)

Cc      XQCD = XSEC(11,3)+XSEC(12,3)+XSEC(13,3)+
Cc     >         XSEC(28,3)+XSEC(53,3)+XSEC(68,3)+
Cc     >         XSEC(81,3)+XSEC(82,3)+XSEC(83,3)

C......VINT(103) = SIGL DIFF CROSS SECTION; VINT(106) = NON-DIFF INELASTIC
C......If target is a nuclei, xsections are taken as the average of N,P:

990   IF(IQ.GT.0) WRITE(MSTU(11),999)
999   FORMAT(/1x,79('-')/
     >   4x,'FRITIOF-FRQPROB reporting: ',/)

      XINEL = VINT(106) + VINT(103)
      XTOT = VINT(101)
      XEL = VINT(102)

      IF(NNUC(IPT(2)).GT.1) THEN

c........jochen
      if(kfr(7).ne.0)CALL FRHARDP(KFI, 4324-IABS(KFT), W, IHAV,-1)
c      CALL FRHARDP(KFI, 4324-IABS(KFT), W, IHAV,-1)
      XTOT = (XTOT+VINT(101))/2.
      XEL = (XEL+VINT(102))/2.
      XINEL = (XINEL+VINT(106)+VINT(103))/2.
      NAM2 = 'Nucleon'//' '
      ENDIF

      WORDS = '(from the input)'
      IF(VFR(10).LE.0.OR.VFR(11).LE.0.OR.IQ.GT.0)THEN
       IF(VFR(10).LE.0.) VFR(10) = XTOT
       IF(VFR(11).LE.0.) VFR(11) = XEL
       words = '(from Block-Cahn fit)'
      ENDIF

       IF(IQ.GT.0) THEN
       WRITE(MSTU(11),2100)
     >       NAM1S,NAM2S,WORDS, VFR(10),VFR(11),XINEL
        WRITE(MSTU(11),3001)
       ENDIF

2100  FORMAT(6X,'Cross sections for ',A8,'-- ',A8,' are',1X,A21,':',/
     >        8x,'Total cross section=', F10.3, ' mb',/
     >        8x,'Elastic cross section=', F10.3, ' mb',/
     >   6X,'Non-double diffractive inelastic xsection= ',F8.3,' mb')
3000  FORMAT(/4x,'Warning! W_CMS=',F6.2,'-- W too small for ',
     > 'hard scattering!',
     > /4x,'Excecution continues with RPS switched off!',
     > /4x,'(Please refer to the RYTHIA parameter PARP(2))',/ )
3001  FORMAT(1x,79('-')/)

      RETURN
      END

C********************************* END FRQPROB **************************


C********************************* FRHARDP *******************************

      SUBROUTINE FRHARDP(KFI,KFT,W,IHAV,IQ)

C....GIVEN PARTICLE KF CODES FOR THE PROJECTILE KFI AND THE TARGET KFT,
C....AND TOTAL CMS ENERGY W, THIS ROUTINE WILL GENERATE PARTON-PARTON
C....PROCESSES (INCLUDING QCD 2->2 PROCESSES, VECTOR BOSONS OR HIGGS
C....PRODUCTIONS, HEAVY QUARKS ECT).  THE GENERATED COLORED OBJECTS
C....ARE TRANSFERED FROM LUJETS TO BLOCK FRJETS, AND THE COLOR-NEUTRAL
C....PARTICLES ARE STORED IN BLOCK FRCNUT. AFTERWARDS THE N IN LUJETS
C....IS RESET TO ZERO.
C......Ihav=0 - non-hard events; Ihav=1, event containing hard process.
C......IQ=-1 OR 1: PYINIT - PYTHIA INITIALIZATION IS MADE;
C......IQ=OTHERS: NO PYTHIA INITIALIZATION
C......FOR IQ<0, Only PYEVENT IS CALLED, FREDIPY will not be callled.

      PARAMETER (KSZJ=4000,KSZ1=20)
      CHARACTER*6 BEAM(2), PARCDE7(15)

      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      integer N,NPAD,K
      double precision P,V
      SAVE /PYJETS/

      COMMON/RYSUBS/MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200)
      COMMON/RYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/FRINTN0/PLI0(2,4),AOP(KSZ1),IOP(KSZ1),NFR(KSZ1)
      COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)
      COMMON/FRCNUT/NR,KR(10,5),PR(10,5),NR0
      COMMON/FRJETS/NJ,KJ(100,5),PJ(100,5)
      COMMON/FRPICKJ/NH,KP(100,5),PP(100,5)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)

      COMMON/FRTORYXTO/KF0(2)

      DIMENSION KCD7(15),KF(2),IKF(2)
      SAVE BEAM
      SAVE /RYSUBS/,/RYPARS/,/FRINTN0/,/FRPARA1/,/FRCNUT/,
     >     /FRJETS/,/FRPICKJ/,/FRTORYXTO/,/LUDAT1/

C....The following is a list of particles available in RYTHIA 5.5.
C....Note for generating hard scattering, all mesons not on the list
C....will be reduced to corresponding pions, and baryons to protons or
C....neutrons. Obviously this can only be a good approximation
C....at high energies when the gluon distribution dominents.

        DATA PARCDE7/
     >  'p+    ','p~-   ','n0    ','n~0   ','pi+   ','pi-   ',
     >  'e-    ','e+    ','nu_e  ','nu_e~ ','mu-   ','mu+   ',
     >  'nu_mu ','nu_mu~','gamma '/
        DATA KCD7/2212,-2212,2112,-2112,211,-211,
     >            11,  -11,  12,  -12,  13, -13,
     >            14,  -14,  22/

      N=0
C...            Identify the beam particles:
c...and store them in the strings beam(1),beam(2)
      KF(1) = KFI
      KF(2) = KFT

      DO L=1, 2
        IF( KF(L).NE.KF0(L) ) THEN
          KF0(L) = KF(L)
          IKF(L) = 0

 5        DO  I=1, 15
            IF(KF(L).EQ.KCD7(I)) THEN
              BEAM(L) = PARCDE7(I)
              IF(IKF(L).EQ.0) IKF(L) = 1
              GOTO 14
            ENDIF
          enddo

 14       IF(IKF(L).EQ.0) THEN
            KFV = IABS(KF(L))
            IF(KFV.GT.1000.AND.KFV.LT.9999) THEN
              KCG=0
              K1 = KFV/1000
              K1R = 2 - MOD(K1,2)
              K2 = (KFV-K1*1000)/100
              K2R = 2 - MOD(K2,2)
              K3 = (KFV-K1*1000-K2*100)/10
              K3R = 2 - MOD(K3,2)
              KLG= MAX(K1R,MAX(K2R,K3R))
              KSM= MIN(K1R,MIN(K2R,K3R))
              KME= K1R+K2R+K3R -KLG-KSM
              KFVR= 1000*KLG+ 100*KME+ 10*KSM+ 2
              IF(KFVR.NE.2212.AND.KFVR.NE.2112) KFVR= 2212
              KF(L) = KFVR *(KF(L)/KFV)
            ELSEIF(KFV.GT.100.AND.KFV.LT.999) THEN
              K1= KFV/100
              K2= (KFV-K1*100)/10
              K3= (KFV-K1*100-K2*10)
              K1R= 2 - MOD(K1,2)
              K2R= 2 - MOD(K2,2)
              KF(L)=(MAX(K1R,K2R)*100+MIN(K2R,K1R)*10+K3)
     >             *(-1)**MAX(K1R,K2R)
              IF(K1R.LT.K2R) KF(L)= -KF(L)
              IF(ABS(KF(L)).NE.211) KF(L) = (KF(L)/IABS(KF(L)))*211
            ENDIF
            IKF(L) = -1
            GOTO 5
          ENDIF

          IF(IKF(L).EQ.0) THEN
            CALL FRMGOUT(0,1,'RYTHIA unrecognized particle',
     >           FLOAT(L),FLOAT(KF(1)),FLOAT(KF(2)),0.,0.)
          ELSEIF(IKF(L).EQ.-1) THEN
            WRITE(MSTU(11),1010)
            IF(L.EQ.1) WRITE(MSTU(11),1012) 'The projectile',
     >           KF0(L), BEAM(L),KF(L)
            IF(L.EQ.2) WRITE(MSTU(11),1012) 'The target    ',
     >           KF0(L), BEAM(L),KF(L)
          ENDIF
        ENDIF
      enddo
 20   CONTINUE

C................................
      NJ = 0
      NH = 0
      MSTI(2) = 0
      MSTI(3) = 0
      MSTI(31) = 0
      IF(ABS(IQ).EQ.1) CALL FRRYINI('USER',BEAM(1),BEAM(2),W)

      CALL RYEVNT

      IF(MSTI(1).GE.91.AND.MSTI(1).LE.95) THEN
      IHAV = 0
      ELSE
      IHAV = 1
      ENDIF

      if(IQ.LT.0.OR.IHAV.EQ.0) RETURN

      CALL FREDIRY

      DO LO = 1, NH
        NJ = NJ+1
        CALL FRVECRC(NJ,LO,1)
      enddo

1010  FORMAT(/,6x,
     >'For the purpose of interfacing with RYTHIA:')
1012  FORMAT(8x,A14,' (code=',I5,')',' is treated as ',A6,
     >       '(code=',I5,')',/ )

      RETURN
      END

C********************************* END FRHARDP ***************************


C********************************* FRSETRY *******************************

      SUBROUTINE FRSETRY(IQ)

C.......TO SET PT_MIN AND CERTAIN SWITCHES FOR PYTHIA................
C....... IQ = -1, SET ON ONLY QCD PROCESSES;
C........   = +1, SET ON QCD + LOW_PT PROCESSES
C:

      PARAMETER (KSZ1=20)
      COMMON/FRINTN0/PLI0(2,4),AOP(KSZ1),IOP(KSZ1),NFR(KSZ1)
      COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)
      COMMON/RYSUBS/MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200)
      COMMON/RYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      SAVE /FRINTN0/,/FRPARA1/,/RYSUBS/,/RYPARS/

      MSTP(61) = 0
      MSTP(63) = 0
      MSTP(65) = 0
      MSTP(71) = 0
      MSTP(91) = 0
      MSTP(111) = 0
      MSTP(122) = 0
      MSTP(31) =5

C....since sea quarks are treated as gluons in current model, heavy
C....quarks are not included...........
      MSUB(81) = 0
      MSUB(82) = 0
      MSUB(83) = 0

      IF(IQ.EQ.-1) THEN
      MSEL = 0
      MSUB(11) = 1
      MSUB(12) = 1
      MSUB(13) = 1
      MSUB(28) = 1
      MSUB(53) = 1
      MSUB(68) = 1
      MSUB(92) = 0
      MSUB(95) = 0
      CKIN(3) = VFR(12)
      CKIN(5) = VFR(12)
      CKIN(6) = VFR(12)
      ELSEIF(IQ.EQ.1) THEN
Cc      MSUB(92) = 1
Cc      MSUB(95) = 1
      MSEL = 1

      ENDIF


C.........multiple interactions....................
      PARP(81) = VFR(12)

      RETURN
      END


C********************************* END FRSETRY ***************************

C********************************* FRVECRC *******************************

C...TO SET VECTORS BETWEEN PJ(JF,), PR(JF,), AND PP(L,)

      SUBROUTINE FRVECRC(JF,L,IQ)

      COMMON/FRJETS/NJ,KJ(100,5),PJ(100,5)
      COMMON/FRCNUT/NR,KR(10,5),PR(10,5),NR0
      COMMON/FRPICKJ/NH,KP(100,5),PP(100,5)
      SAVE /FRJETS/,/FRCNUT/,/FRPICKJ/

      IF(IQ.EQ.1) THEN
       DO 11 J=1,5
       KJ(JF,J) = KP(L,J)
11     PJ(JF,J) = PP(L,J)
      ELSEIF(IQ.EQ.-1) THEN
       DO 21 J=1,5
      KP(L,J) = KJ(JF,J)
21    PP(L,J) = PJ(JF,J)
      ELSEIF(IQ.EQ.2) THEN
       DO 31 J=1,5
       KR(JF,J) = KP(L,J)
31     PR(JF,J) = PP(L,J)
      ELSEIF(IQ.EQ.-2) THEN
       DO 41 J=1,5
      KP(L,J) = KR(JF,J)
41    PP(L,J) = PR(JF,J)
      ENDIF

      RETURN
      END

C********************************* END FRVECRC ***************************


C********************************* FREDIRY *******************************

      SUBROUTINE FREDIRY

C....TO PICK OUT THE SCATTERED PARTONS OUT OF PYTHIA'S EVENT RECORD.
C....NOTE ONE MUST HAVE SET MSTP(61,63,65,71,81,111) ALL TO 0 .
C....THE PARTONS PICKED ARE STORED IN BLOCK FRPICKJ.

C....KFR19 controls the assignment of the partons to the original nucleons:
C..  (currently in effect: KFR19=1)
C..   =1 The hardest pairs are assigned according to PYTHIA, the rest randomly.
C..   =2 Both pairs of partons are assigned randomly to one of the nucleon;
C..   =0 assignment to be made by inspecting the Feynman diagrams in FRHQSGN.

C...  Side-1: KP(J,3) = 1;  Side-2: KP(J,3) = 2  - for the hardest pair;
C...  Side-1: KP(J,3) = -1;  Side-2: KP(J,3) = -2 - for the softer pairs.

      PARAMETER (KSZJ=4000, KSZ1=20)
      COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)

      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      integer N,NPAD,K
      double precision P,V
      SAVE /PYJETS/

      COMMON/RYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/RYSUBS/MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200)
      COMMON/FRPICKJ/NH,KP(100,5),PP(100,5)
      SAVE /FRPARA1/,/RYPARS/,/RYSUBS/,/FRPICKJ/

      kfr19 = 1

      IFEL = 0
      DO 10 L = 1, MSTI(3)
      IF(L.EQ.1) THEN
      LINE = MSTI(7)
      IF(LINE.EQ.0) LINE = MSTI(8)
      ELSEIF(L.EQ.2) THEN
      LINE = MSTI(8)
      ELSE
      IFEL = 1
      ENDIF

      IF(LINE.EQ.0) IFEL = 1
      IF(IFEL.EQ.1) CALL FRMGOUT(0,1,'Error in FREDITRY:',
     >      FLOAT(MSTI(3)),FLOAT(MSTI(7)),FLOAT(MSTI(8)),0.,0.)
      DO 30 J=1,5
      PP(L,J) = P(LINE,J)
30    KP(L,J) = K(LINE,J)

      IF(L.EQ.1.AND.KFR19.EQ.2) THEN
      KP(1,3) = INT(1.5+PYR(0))
      ELSEIF(KFR19.EQ.2) THEN
      KP(L,3) = KP(1,3)
      ENDIF

10    CONTINUE


            IF(KFR19.EQ.1) THEN
      THAT = (PARI(33)*PARI(11)/2.-PP(1,4))**2-
     >   (PARI(33)*PARI(11)/2.-PP(1,3))**2-PP(1,2)**2-PP(1,1)**2
      T12 = ABS( THAT - PARI(15))
      U12 = ABS( THAT - PARI(16))

        if(T12.LE.U12) THEN
      KP(1,3) = 1
      KP(2,3) = 2
      ELSE
      KP(1,3) = 2
      KP(2,3) = 1
        ENDIF
            ENDIF


      NH = MSTI(3)

      IF(MSTP(81).EQ.0.or.KFR(7).NE.2) RETURN

      DO 39 L = MSTI(8)+1, N
      IF(K(L,3).EQ.0) THEN

      NH = NH + 1
      IF(NH.GT.100)CALL FRMGOUT(0,0,'Extend blocks FRPICKJ and FRJETS',
     >                         float(NH),0.,0.,0.,0.)
      ISIDE1 = - INT(1.5+ PYR(0))
      ISIDE2 = -3 - ISIDE1
       DO 35 LI = L+1, N
       XB = ABS(P(L,1)+P(LI,1))
       YB = ABS(P(L,2)+P(LI,2))
       IF(XB.GT.0.0001.OR.YB.GT.0.0001) GOTO 35
       NH = NH + 1
       DO 40 J=1,5
       PP(NH-1,J) = P(L,J)
       KP(NH-1,J) = K(L,J)
       PP(NH,J) = P(LI,J)
40     KP(NH,J) = K(LI,J)
       KP(NH-1,3) = ISIDE1
       IF(KFR19.EQ.2) THEN
       KP(NH,3) = ISIDE1
       K(L,3) = ISIDE1
       K(LI,3) = ISIDE1
       ELSE
       KP(NH,3) = ISIDE2
       K(L,3) = ISIDE1
       K(LI,3) = ISIDE2
       ENDIF
       GOTO 39
35     CONTINUE
      ENDIF

39    CONTINUE

      RETURN
      END

C********************************* END FREDIRY ***************************

C********************************* FRHPLIS *******************************

      SUBROUTINE FRHPLIS

C.........TO LIST THE HARD PARTONS EXTRACTED FROM PYTHIA EVENT RECORD
C.........in case one wants to examing it                       ......

      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/RYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/FRCNUT/NR,KR(10,5),PR(10,5),NR0
      COMMON/FRJETS/NJ,KJ(100,5),PJ(100,5)
      COMMON/FRPICKJ/NH,KP(100,5),PP(100,5)
      SAVE /LUDAT1/,/RYPARS/,/FRCNUT/,/FRJETS/,/FRPICKJ/

      WRITE(MSTU(11),*) '===================================='
      WRITE(MSTU(11),10) MSTI(5)
      WRITE(MSTU(11),15) MSTI(11),MSTI(12)
      WRITE(MSTU(11),20) MSTI(1), MSTI(2)
      WRITE(MSTU(11),30) MSTI(3),MSTI(31)
      WRITE(MSTU(11),50) MSTI(7),MSTI(8)
      WRITE(MSTU(11),60) MSTI(13),MSTI(14)
      WRITE(MSTU(11),70) MSTI(15),MSTI(16),(MSTI(L),l=21,24)
      WRITE(MSTU(11),80) PARI(33), PARI(34)

      WRITE(MSTU(11),*) ' NH=',NH, '    -- FRPICKJ '
      IF(NH.GT.0) THEN
      DO 301 J=1, NH
301   WRITE(MSTU(11),3401) J, (KP(J,L),L=1,5),(PP(J,L),L=1,5)
      ENDIF
      WRITE(MSTU(11),*) ' NJ=',NJ, '    -- FRJETS '
      IF(NJ.GT.0) THEN
      DO 303 J=1, NJ
303   WRITE(MSTU(11),3401) J, (KJ(J,L),L=1,5),(PJ(J,L),L=1,5)
      ENDIF

      WRITE(MSTU(11),*) ' NR=',NR, '    -- FRCNUT '
      IF(NR.GT.0) THEN
      DO 305 J=1, NR
305   WRITE(MSTU(11),3401) J, (KR(J,L),L=1,5),(PR(J,L),L=1,5)
      ENDIF

10    FORMAT( /' PARTON LIST AT ',I6,2X,'-th CALL TO RYTHIA',/)
15    FORMAT( ' Collision between ',I5,' & ',I5,/)
20    FORMAT( ' Subprocess type - MSTI(1,2): ',2I6 )
30    FORMAT( ' No. of partons produced: ',I6, 2X,
     >          ' No. of interactions: ', I3 )
50    FORMAT( ' Parton line number - MSTI(7,8): ',2I6 )
60    FORMAT( ' Initial shower initiaters - MSTI(13,14): ',2I5 )
70    FORMAT( ' The process is ',I4,' +',I4,' ->',I4,' +',I4,
     >     I4,' +',I4 )
80    FORMAT( ' X_1, X_2 = ',2G13.6, / )

3401  FORMAT(1X,I3,';',2x,5I6,2X,5F10.4, ' -- K, P ' )
      RETURN
      END

C********************************* END FRHPLIS ***************************



C*************************************************************************
C**                                                                     **
C**   This is the package for auxililary subroutines                    **
C**                                                                     **
C*************************************************************************

C********************************* FRGAUSS ****************************

      SUBROUTINE FRGAUSS(P2,V,PMAX)

C.... TO RETURN A VALUE P2 WHICH HAS A MAXIMUM SET BY PMAX, AND A
C.... 2-D GAUSSIAN DISTRIBUTION WITH WIDTH V, i.e., e^(-P2/V)dP2, 0<P2<PMAX.
C.... set PMAX < 0 if PMAX should be infinity.

      P2 = 0
      IF(V.LE.1.E-8) RETURN

      IF(PMAX.LT.0) THEN
      A = 1.
      ELSEIF(PMAX.LT.1.E-9) THEN
      RETURN
      ELSE
      A = 1. - FRREX(-PMAX/V)
      ENDIF

10    P2 = -V* LOG(MAX(1.E-20,1. - A*PYR(0)))
      IF(P2.LT.0.) GOTO 10

      RETURN
      END

C********************************* END FRGAUSS ************************

C********************************* FRBETAV ****************************

      SUBROUTINE FRBETAV(ID,DBETA,DP)

C...FOR GIVEN PAIR OF MEMENTA DP(2,4), THIS IS TO FILL THE ARRAY
C...DBETA(3) WHICH ARE THE BETA FACTORS FOR THE CMS FRAME.
C.. ID = 0 FOR NORMAL 4-VECTORS
C.. ID = 1 IF THE VECTORS ARE LIGHT-CONE FORM: PX,PY,P-,P+

      IMPLICIT DOUBLE PRECISION (D)
      DIMENSION DBETA(3),DP(2,5),DR(2,4)

      DO 10 I = 1, 2
      DR(I,1) = DP(I,1)
      DR(I,2) = DP(I,2)
      IF(ID.EQ.0) THEN
      DR(I,3) = DP(I,3)
      DR(I,4) = DP(I,4)
      ELSE
      DR(I,3) = (DP(I,4)-DP(I,3))/2.D0
      DR(I,4) = (DP(I,4)+DP(I,3))/2.D0
      ENDIF
10    CONTINUE

C..... BETA_X, BETA_Y, BETA_Z
      DESUM = (DR(1,4)+DR(2,4))
      IF(DESUM.LE.0.) STOP 'FRBETA: E 0'
      DO 15 I = 1, 3
15    DBETA(I) = (DR(1,I)+DR(2,I))/DESUM

      RETURN
      END

C********************************* END FRBETAV ************************

C********************************* FRTOCMS ****************************

        SUBROUTINE FRTOCMS(ID,IQ,DP,DBETAO)

C... IN DOUBLE PRECISION.
C... GIVING A PAIR OF MEMONTA DP(2,4), A CALL TO THIS ROUTINE WILL
C... TRANSFORM THE MOMENTA INTO THEIR CMS FRAME.
C... IF ID = 0, DP IS ASSUMED TO BE THE ORDINARY 4-VECTOR;
C...    ID = 1, DP IS ASSUMED TO BE THE LIGHT-CONE FORM: PX,PY,P-,P+.
C... IF IQ > 0, A BOOST IS PERFORMED TO THE CMS FRAM; DBETA is OUTPUT of
C...             the beta factor used;
C...    IQ = 0, A BOOST is done on DP with a known Dbeta
C...    IQ < 0, A BOOST is done on DP with Dbeta = -DBETA, i.e., inverse boost,
C...             here DBETA must be given as INPUT.
C... SO, FOR EXEMPLE,
C: ...  CALL FRTOCMS(ID, 1, DP, DBETA)    -- TO THE CMS FRAME;
C: ...  CALL FRTOCMS(ID, -1, DP, DBETA)   -- BACK TO THE ORGINAL FRAME.
C:
C: ...SPECIALLY FOR |IQ|=2, PJ IN FRJETS BLOCK AND PR(NR0-NR)
C: ...IN FRCNUT BLOCK ARE ALSO BOSTED.
C: ...NR-NR0<1 IS ASSUMED.

      IMPLICIT DOUBLE PRECISION (D)
      COMMON/FRCNUT/NR,KR(10,5),PR(10,5),NR0
      COMMON/FRJETS/NJ,KJ(100,5),PJ(100,5)
      DIMENSION DP(2,5),DBETAO(3),DBETA(3),DRV(4),DBP(100)
      SAVE /FRCNUT/,/FRJETS/

C..... BETA_X, BETA_Y, BETA_Z
      DBET2 = 0.

      IF(ID.EQ.0) THEN
      DESUM = (DP(1,4)+DP(2,4))
      DZSUM = DP(1,3)+DP(2,3)
      else
      DESUM = (DP(1,4)+DP(1,3))/2.D0 +(DP(2,4)+DP(2,3))/2.D0
      DZSUM = (DP(1,4)-DP(1,3))/2.D0 +(DP(2,4)-DP(2,3))/2.D0
      ENDIF

      IF(DESUM.LE.0.) STOP 'FRCMS: E LESS THAN 0'

      DO 5 I = 1, 3
      IF(IQ.GT.0) THEN
       IF(I.LE.2) DBETA(I) = (DP(1,I)+DP(2,I))/DESUM
       IF(I.EQ.3) DBETA(I) = DZSUM/DESUM
       DBETAO(I) = DBETA(I)
      ELSEIF(IQ.LT.0) THEN
       DBETA(I) = - DBETAO(I)
      ELSEIF(IQ.EQ.0) THEN
       DBETA(I) = DBETAO(I)
      ENDIF
      DBET2 = DBET2 + DBETA(I)**2
5     CONTINUE

      IF(DBET2.GE.1.D0) STOP 'FRCMS: BETA = 1'
      IF(DBET2.LT.1.D-8) RETURN
      DGAMA = 1.D0/DFRSQR(1.D0-DBET2, 'DGA678')
      DEFF = DGAMA/(1.D0+DGAMA)


      K=0
9     K=K+1
      IF(K.EQ.1) THEN
      Iup=2
      ELSEIF(K.EQ.2.OR.K.EQ.4) THEN
      Iup=NJ
      ELSEIF(K.EQ.3) THEN
         IF(NR-NR0.GE.2) CALL FRMGOUT(0,1,'Check NR,NR0:',
     >                 FLOAT(NR),FLOAT(NR0),0.,0.,0.)
      Iup=NR-NR0+1
      ENDIF

      DO 888 I = 1, IUP
       IF(K.EQ.1) THEN
      DRV(1) = DP(I,1)
      DRV(2) = DP(I,2)
         IF(ID.EQ.0) THEN
      DRV(3) = DP(I,3)
      DRV(4) = DP(I,4)
         ELSE
      DRV(3) = (DP(I,4)-DP(I,3))/2.D0
      DRV(4) = (DP(I,4)+DP(I,3))/2.D0
         ENDIF
      ELSEIF(K.EQ.2) THEN
       DO 11 J=1,4
11     DRV(J) = PJ(I,J)
      ELSEIF(K.EQ.3) THEN
      IR = NR0+I-1
       DO 12 J=1,4
12     DRV(J) = PR(IR,J)
      ENDIF

      DBP(I) = 0.
       DO 30 J = 1, 3
30     DBP(I) = DBP(I) + DBETA(J)* DRV(J)
       DO 35 J = 1, 3
35     DRV(J) = DRV(J)+ (DEFF*DBP(I)-DRV(4))*DGAMA*DBETA(J)
40    DRV(4) = DGAMA* (DRV(4)-DBP(I))

       IF(K.EQ.1) THEN
      DP(I,1) = DRV(1)
      DP(I,2) = DRV(2)
        IF(ID.EQ.0) THEN
      DP(I,3) = DRV(3)
      DP(I,4) = DRV(4)
        ELSE
      DP(I,4) = DRV(4) + DRV(3)
      DP(I,3) = DRV(4) - DRV(3)
        ENDIF
      ELSEIF(K.EQ.2) THEN
      DO 41 J=1,4
41    PJ(I,J) = DRV(J)
      ELSEIF(K.EQ.3) THEN
      DO 42 J=1,4
42    PR(NR0+I-1,J) = DRV(J)
       ENDIF
888   CONTINUE

      IF(ABS(IQ).EQ.2.AND.NJ.GT.0.AND.K.LE.1) GO TO 9
      IF(ABS(IQ).EQ.2.and.NR-NR0.GE.0.AND.K.LE.2) GO TO 9

      RETURN
      END

C********************************* END FRTOCMS ************************

C********************************* FRPOLAR ****************************

      SUBROUTINE FRPOLAR(DTHE,DPHI,DP)

C....FOR A GIVEN PAIR OF MEMENTA DP(2,4), WHICH ARE IN LIGHT-CONE FORM
C....AND IN THEIR CMS FRAME, THIS ROUTINE IS TO FIND THE POLAR ANGEL
C....(THETA,PHI) WITH WHICH A FOLLOW UP CALL TO FRROTAR (ROTATION)
C....WILL MAKE THE PAIRWISE MOMENTA LIE ON THE NEW Z-AXES.
C.......................................................................

      IMPLICIT DOUBLE PRECISION (D)
      DIMENSION DP(2,5)
C............................ANGLES OF THE MOMENTUM VECTOR ...............
      DTHE = 0D0
      DPHI = 0D0
      DBP3=.5D0*(DP(1,4)-DP(1,3))
c..........dbp3=p_z
      DXY=(DP(1,1)**2+DP(1,2)**2)
      DXYZ = DXY + DBP3**2
      IF(DXYZ.EQ.0.D0) RETURN

      DCTH=DBP3/DFRSQR(DXYZ, 'dcthiowe')
      DCTH=DMAX1(-1.D0,DMIN1(DCTH,1.D0))
      DTHE=DACOS(DCTH)
      IF(DTHE.EQ.0.D0.OR.DXY.EQ.0D0) RETURN
      DCPH = DP(1,1)/DFRSQR(DXY, 'dxyeuw')
      DPHI = DACOS(DCPH)
      IF(DP(1,2).LT.0.D0) DPHI = - DPHI

      RETURN
      END

C********************************* END FRPOLAR ************************

C********************************* FRROTAR ****************************

       SUBROUTINE FRROTAR(DTHE,DPHI,IQ,DP)

C...ROTATIONS
C:  THE ROTATION MEANS:
C:  IF IQ >0, THE COORDINATES ARE FIRST ROTATED DTHE ANGLE ABOUT
C:  Y-AXES AND THEN AN ANGLE DPHI AROUND Z-ZXES. THE EFFECT IS THAT
C:  A VECTOR ORIGINALLY AT (DTHE,DPHI) POLAR ANGLE WILL BE MOVED TO
C:  THE NEW Z-AXE.
C:  IF IQ <0, THE COORDINATES ARE FIRST ROTATED -DPHI ANGLE ABOUT
C:  Z-AXES AND THEN AN ANGLE -DTHE AROUND Y-ZXES.
C:  THIS IS TO COUNTER THE EFFECT OF A PREVIOUS "IQ>0" ROTATION.
C:  THE MORAL IS, IF CALLED TWICE LIKE THIS,
C: ...  CALL FRROTAR(DTHE,DPHI,1,DP)
C: ...  CALL FRROTAR(DTHE,DPHI,-1,DP)
C:  WE ARE BACK TO THE ORIGINAL FRAME AND NOTHING IS CHANGED.
C:
C: ...SPECIALLY FOR |IQ|=2, PJ,PJI IN FRJETS BLOCK AND PR(NR0-NR)
C: ...IN FRCNUT BLOCK ARE ALSO ROTATED.

      IMPLICIT DOUBLE PRECISION (D)
      COMMON/FRCNUT/NR,KR(10,5),PR(10,5),NR0
      COMMON/FRJETS/NJ,KJ(100,5),PJ(100,5)

      COMMON/DUMYFROTA/ DPV(4)
      DIMENSION DP(2,5)
      SAVE /FRCNUT/,/FRJETS/,/DUMYFROTA/

      IF(DTHE**2+DPHI**2.LT.1D-20) RETURN

      K=0
9     K=K+1
      IF(K.EQ.1) THEN
      IUP=2
      ELSEIF(K.EQ.2.OR.K.EQ.4) THEN
      IUP=NJ
      ELSEIF(K.EQ.3) THEN
      IUP=NR-NR0+1
      ENDIF

      DO 120 I=1,IUP
        IF(K.EQ.1) THEN
          DO 100 J=1,2
100     DPV(J) = DP(I,J)
          DPV(3)=.5D0*(DP(I,4)-DP(I,3))
          DPV(4) =.5D0*(DP(I,4)+ DP(I,3))
        ELSEIF(K.EQ.2) THEN
          DO 101 J=1,4
101     DPV(J) = PJ(I,J)
        ELSEIF(K.EQ.3) THEN
         IF(NR-NR0.GE.2) CALL FRMGOUT(0,1,'CHECK here NR,NR0:',
     >           FLOAT(NR),FLOAT(NR0),0.,0.,0.)
          DO 103 J=1,4
103     DPV(J) = PR(NR0+I-1,J)
       ENDIF

       IF(IQ.GT.0) THEN
       CALL FRROTAZ(DPHI,DPV)
       CALL FRROTAY(DTHE,DPV)
       ELSEIF(IQ.LT.0) THEN
       CALL FRROTAY(-DTHE,DPV)
       CALL FRROTAZ(-DPHI,DPV)
       ENDIF

        IF(K.EQ.1) THEN
        DP(I,1) = DPV(1)
        DP(I,2) = DPV(2)
          DP(I,3)=DPV(4) - DPV(3)
          DP(I,4)=DPV(4) + DPV(3)
        ELSEIF(K.EQ.2) THEN
          DO 111 J=1,4
111     PJ(I,J) = DPV(J)
        ELSEIF(K.EQ.3) THEN
          DO 113 J=1,4
113     PR(NR0+I-1,J) = DPV(J)
      ENDIF

120     CONTINUE

      IF(ABS(IQ).EQ.2.AND.NJ.GT.0.AND.K.LE.1) GO TO 9
      IF(ABS(IQ).EQ.2.and.NR-NR0.GE.0.AND.K.LE.2) GO TO 9

      RETURN
      END

C********************************* END FRROTAR ************************

C********************************* FRROTAY ****************************

      SUBROUTINE FRROTAY(DTHE, DPV)

C:  ROTATE COORDINATES AROUND Y-AXIS BY AN ANGLE DTHE
C:  DPV(3) GIVES THE SPACE COMPONENTS OF A VECTOR.

      IMPLICIT DOUBLE PRECISION (D)
      DIMENSION DPV(4)
      IF(DTHE**2.LT.1D-20) RETURN

      DPVX=DPV(1)*DCOS(DTHE)-DPV(3)*DSIN(DTHE)
      DPVY = DPV(2)
      DPVZ=DPV(1)*DSIN(DTHE)+DPV(3)*DCOS(DTHE)

      DPV(1) = DPVX
      DPV(2) = DPVY
      DPV(3) = DPVZ

      RETURN
      END

C********************************* END FRROTAY ************************

C********************************* FRROTAZ ****************************

      SUBROUTINE FRROTAZ(DPHI, DPV)

C:  ROTATE COORDINATES AROUND Z-AXES BY AN ANGLE DPHI
C:  DPV(3) GIVES THE SPACE COMPONENTS OF A VECTOR.

      IMPLICIT DOUBLE PRECISION (D)
      DIMENSION DPV(4)
      IF(DPHI**2.LT.1D-20) RETURN

      DPVX=DPV(1)*DCOS(DPHI)+DPV(2)*DSIN(DPHI)
      DPVZ = DPV(3)
      DPVY=-DPV(1)*DSIN(DPHI)+DPV(2)*DCOS(DPHI)

      DPV(1) = DPVX
      DPV(2) = DPVY
      DPV(3) = DPVZ

      RETURN
      END

C********************************* END FRROTAZ ************************

C********************************* FRBOOT1 ****************************

      SUBROUTINE FRBOOT1(ID,DPV,DBETA)

C... TO BOOST AN SINGLE MOMENTA BY A DBETA(3) FACTOR.
C... ID =0, DPV(1-4)=P_X, P_Y, P_Z, E;
C... ID =1, DPV(1-4)=P_X, P_Y, P_, P+

      IMPLICIT DOUBLE PRECISION (D)
      DIMENSION DBETA(3), DPV(4), DRV(4)

      DBET2 = 0.D0
      DO 10 J = 1, 3
10    DBET2 = DBET2 + DBETA(J)**2
      IF(DBET2.LT.1.D-10) RETURN
      IF(DBET2.GT.1.D0) THEN
      CALL FRMGOUT(0,1,' FRBOOT1: CHECK BETA > 1',
     >    REAL(DBETA(1)),REAL(DBETA(2)),REAL(DBETA(3)),REAL(DBET2),0.)
      ENDIF
      DBET = DSQRT(DBET2)

        IF(DBET.GT.0.99999999D0) THEN
      DO 13 J=1, 3
13    DBETA(J) = DBETA(J)* 0.99999999D0/DBET
        DBET=0.99999999D0
      DBET2 = DBET**2
        ENDIF

      DGAMA = 1.D0/DFRSQR(1.D0-DBET2, 'UIOP09')
        DEFF = DGAMA/(1.D0+DGAMA)

      DRV(1) = DPV(1)
      DRV(2) = DPV(2)
      IF(ID.EQ.0) THEN
      DRV(3) = DPV(3)
      DRV(4) = DPV(4)
      ELSE
      DRV(3) = (DPV(4) - DPV(3))/2D0
      DRV(4) = (DPV(4) + DPV(3))/2D0
      ENDIF

       DBP = 0.
       DO 25 I = 1, 3
25     DBP = DBP + DBETA(I)* DRV(I)

C..........................................BOOST..............
      DO 30 J = 1, 3
  30    DRV(J) = DRV(J) +(DEFF*DBP -DRV(4))*DGAMA*DBETA(J)
      DRV(4) = DGAMA* (DRV(4)-DBP)

      DPV(1) = DRV(1)
      DPV(2) = DRV(2)
      IF(ID.EQ.0) THEN
      DPV(3) = DRV(3)
      DPV(4) = DRV(4)
      ELSE
      DPV(3) = (DRV(4) - DRV(3))
      DPV(4) = (DRV(4) + DRV(3))
      ENDIF


      RETURN
      END

C********************************* END FRBOOT1 ************************


C********************************* FRMGOUT *****************************

      SUBROUTINE FRMGOUT(ID,ILIST,MESG,A,B,C,D,E)

C...FOR GENERAL MESSAGE PRINT OUT ....................................
C...ID:  -50 - 50  ID number for the error:
C.....   If ID=0, the execution will be stopped upon MESG printout;
C.....   If ID>0, the execution will continue but the MESG printout
C.....     is limited to NTERM times;
C.....   If ID<0, the execution will continue with the MESG printout
C......     but execution stops if it repeats NTERM times.
C......A,B,C,D,E some variables to be printed out for inspection,
C.........          however, the will not be printed if they are all 0.
C...ILIST=1: a full list including FR-status, event list, hard parton list
C........          is given.
C........      =0: full list suppressed.
C..
C...MGO(1) - energy nonconservation in FRRINGO;
C...MGO(2) - energy nonconservation in FRPPART;


      PARAMETER (KSZ1=20)
      COMMON/FRINTN0/PLI0(2,4),AOP(KSZ1),IOP(KSZ1),NFR(KSZ1)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/FRJETS/NJ,KJ(100,5),PJ(100,5)

      COMMON/FRMGOV/MGO(-5:10)

      CHARACTER *(*) MESG

      SAVE /FRINTN0/,/LUDAT1/,/FRJETS/,/FRMGOV/
      DATA NTERM/10/

   2  FORMAT(2X, A,3X,'Error ID=',I3,1x,'Count=',I3, /)

      MGO(ID) = MGO(ID)+1
      IF(MGO(ID).LE.NTERM) IOP(16)=1
      IF(ID.GT.0.AND.MGO(ID).GT.NTERM) RETURN

      IQAB = 1
      IF(ABS(A)+ABS(B)+ABS(C)+ABS(D)+ABS(E).EQ.0.) THEN
      IQAB = -1
      ENDIF
      WRITE(MSTU(11),10) NFR(1), IOP(1)

      WRITE(MSTU(11),2) MESG, ID, MGO(ID)

      IF(IQAB.GT.0) WRITE(MSTU(11),*) A,B,C,D,E

      CALL FRVALUE(0)

      IF(ILIST.EQ.1) THEN
      CALL LULIST(2)
      IF(NJ.GT.0) CALL FRHPLIS
      ENDIF

      WRITE(MSTU(11),20)

      IF(ID.EQ.0.OR.(ID.LT.0.AND.MGO(ID).GE.NTERM)) THEN
      WRITE(MSTU(11),*) ' Severe! EXECUTION STOPPED BY FRMGOUT'
      STOP 'FRMGOUT:'
      ENDIF

10    FORMAT(/72('*') /72('?')/,' POSSIBLY AN ERROR! AT EVENT NO. ',I7,
     > 3X,'SUBCOLLISION ',I4)
20    format(72('|')/72('*'),/)

      RETURN
      END

C********************************* END FRMGOUT *************************

C********************************* FRDOICT *******************************

C.... TO MANAGE COUNTING ......................................
C.... 1 IS ADDED TO ICT(I) IF the argument is > 0.
C.... The record in ICT is cleared if I <= -1; and
C.... (This routine is not listed in the manual)

      SUBROUTINE FRDOICT(I)

      COMMON/FRCONT2/ICT(10),ICTT(10)
      SAVE /FRCONT2/

      IF(I.LE.-1) THEN
      DO 10 L=1, 10
10    ICT(L) = 0
      RETURN
      ENDIF

      IF(I.GT.0) THEN
      ICT(I) = ICT(I) + 1
      ICTT(I) = ICTT(I) + 1
      ENDIF

      RETURN
      END

C********************************* END FRDOICT ***************************

C******************************** AUXILLIARY ROUTINES *******************

      REAL FUNCTION FRSQR( X, MESSAGE )

C....Optional character 'MESSAGE' helps to identify the source of error.
C....Allow a little numerical error margin...

      CHARACTER*(*) MESSAGE
      IFLAG = 0
      IF(X.LT.-0.001) THEN
      IFLAG = 1
      WRITE(6,*) X, ' --SQRT-NEGATIVE VALUE '
      write(6,100) MESSAGE
100   FORMAT( A )
*martin
*      STOP 'FRSQR: `NEG-ROOT'''
      ENDIF
      FRSQR = SQRT(MAX(X,0.))
      RETURN
      END

      DOUBLE PRECISION FUNCTION DFRSQR( DX, MESSAGE )
      IMPLICIT DOUBLE PRECISION (D)
      CHARACTER*(*) MESSAGE
      IFLAG = 0
      IF(DX.LT.-0.001D0) THEN
      IFLAG = 1
      WRITE(6,*) DX, ' --SQRT-NEGATIVE VALUE '
      write(6,100) MESSAGE
100   FORMAT( A )
      STOP 'DFRSQR: `NEG-D_ROOT'''
      ENDIF
      DFRSQR = DSQRT(DMAX1(DX,0.D0))
      RETURN
      END

C..........................................................

      FUNCTION FRREX(X)

C.......To take care the over_under_flow problem in large exponentials

      ARG = MIN(ABS(X), 80.0)
      IF(X.LT.0.) ARG = -ARG

      FRREX = EXP(ARG)

      RETURN
      END



C********************************* FRUPCAS ***************************

      SUBROUTINE FRUPCAS(STR)

C.... Convert a string character (length<20) into UPPER CASE one
C.... STR must not contain more than 3 spaces between the characters.

      CHARACTER STR*(*),CHALP(2)*26,STR0*23
      DATA CHALP/'abcdefghijklmnopqrstuvwxyz',
     &'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/

      I=0
      IPH=0
      STR0 = STR
10    I=I+1
      IF(STR0(I:I).EQ.' ') IPH=IPH+1
      IF(IPH.GT.3) GOTO 30
      IF (STR0(I:I) .GE. 'a'  .AND.  STR0(I:I) .LE. 'z')  THEN
CC          STR(I:I) = CHAR(ICHAR(STR(I:I)) - '20'X)
      DO 20 J=1,26
20      IF(STR(I:I).EQ.CHALP(1)(J:J)) STR(I:I)= CHALP(2)(J:J)
      ENDIF
      GOTO 10

30    RETURN
      END

C********************************* END FRUPCAS ***********************

C********************************* FRLOOPU ***************************

      SUBROUTINE FRLOOPU(*,I,IMAX,MESSAGE)

C....Handle loops to avoid infinite loop.
C....I=loop index;  IMAX=maximum number of looping;
C....MESSAGE = character code to mark the loop.

      PARAMETER (KSZ1=20)
      CHARACTER*(*) MESSAGE
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/FRINTN0/PLI0(2,4),AOP(KSZ1),IOP(KSZ1),NFR(KSZ1)
      SAVE /LUDAT1/,/FRINTN0/

      IF(I.LE.IMAX) THEN
      RETURN 1
      ELSE
      WRITE(MSTU(11),110) NFR(1),IOP(1),MESSAGE
      RETURN
      ENDIF

110   FORMAT(/,15('?'),2x,'Event-subcollision number:',I6,'_',I4,
     >      10('?')/,4X,'Loop aborted at ', A,/)

      END

C***************************** END FRLOOPU ***************************

C*************************** FRVALUE *********************************

C...To output the values of FRITIOF parameters ...................
C... IF IQ=0, THE OUTPUT IS WRITTEN ON MSTU(11);
C... IF IQ>0, IT WILL WRITES ON AN AUX FILE NAMED 'Oxchk', which
C... is refreshed at an interval of every INUM events.
C... IF IQ<0, as in IQ>0 but "Execution completed" will ALSO be printed.

      SUBROUTINE FRVALUE(iq)

      parameter (KSZ1=20)
      COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)
      COMMON/FRINTN0/PLI0(2,4),AOP(KSZ1),IOP(KSZ1),NFR(KSZ1)
      COMMON/FRCONT2/ICT(10),ICTT(10)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/FRMGOV/MGO(-5:10)

      SAVE ICO
      SAVE /FRPARA1/,/FRINTN0/,/FRCONT2/,/LUDAT1/,/FRMGOV/
      DATA ICO,INUM /0,20/

      IF(IQ.EQ.0) THEN
      IFILE = MSTU(11)
      ELSE
      ICO = ICO + 1
      IFILE = 111
       IF(ICO.EQ.1) THEN
       OPEN(IFILE,FILE='Oxchk',STATUS='UNKNOWN')
       ELSEIF(MOD(ICO,INUM).EQ.1) THEN
       OPEN(IFILE,FILE='Oxchk',STATUS='OLD')
       ENDIF
      ENDIF

      WRITE(IFILE,100)
100      FORMAT(/,79('='),/2x,'Fritiof Status Report',/)

      WRITE(IFILE,106)
      WRITE(IFILE,166)(KFR(L),L=1,15),(VFR(L),L=1,15),(IOP(L),L=1,15)

      WRITE(IFILE,103)(NFR(L),L=1,10),(ICT(L),L=1,10),(ICTT(L),L=1,10)
103   FORMAT(79('-'),/,
     >   4X,7X,'1|',6X,'2|',6X,'3|',6X,'4|',6X,'5|',5X,'6|',5X,'7|',
     >     5x,'8|',5X,'9|',4X,'10|',/,79('-'),/,
     >   1X,'NFR',5I8,5I7,/,
     >   1X,'ICT',5I8,5I7,/,'ICTT',5I8,5I7,/,79('-'),/ )

          ratio1 = 0.
          ratio2 = 0.
          if(NFR(3).GT.0) RATIO1 = FLOAT(NFR(4))/FLOAT(NFR(3))
          if(NFR(1).GT.0) RATIO2 = FLOAT(NFR(5))/FLOAT(NFR(1))

      WRITE(IFILE,*)' Percent of collisions having hard scattering: ',
     >     ratio1
      WRITE(IFILE,*)' Percent of events having hard scattering: ',
     >     ratio2

      IF(MGO(1).GT.0) WRITE(IFILE, 171) MGO(1)
      IF(MGO(2).GT.0) WRITE(IFILE, 172) MGO(2)

CC      call datecpu(IFILE)

      IF(IQ.LT.0) WRITE(IFILE, *) ' &&& EXECUTION COMPLETES &&&&'

        IF(IQ.NE.0.AND.MOD(ICO,INUM).EQ.0) CLOSE(IFILE)

106   FORMAT(20X,'KFR(L), VFR(L), IOP(L)  ',30X,/,78('-') )
166   FORMAT(
     >   7X,'1|',3X,'2|',3X,'3|',3X,'4|',3X,'5|',3X,'6|',3X,'7|',
     >     3x,'8|',3X,'9|',2X,'10|',2X,'11|',2X,'12|',2X,'13|',
     >     2X,'14|',2X,'15|',/,78('-'),/,
     >   'KFR', 15I5,/,'VFR', 15F5.2,/,'IOP', 15I5,/,78('-'),/ )

171   FORMAT(/,2X,'No. of errors - energy nonconserv. in FRRINGO:',I6)
172   FORMAT(/,2X,'No. of errors - energy nonconserv. in FRPPART:',I6,/)

      return
      end

C***************************** END FRVALUE *****************************

C***********************************************************************
C..............Subroutines originally from PYTHIA.......................
C... modified slightly for FRITIOF to accomodate some changes in the way
C... meson-nucleon cross sections are handled.
C.......................................................................

      SUBROUTINE FRRYINI(FRAME,BEAM,TARGET,WIN)

C...This routine is identical to PYTHIA's PYINIT except the routine
C...for xsections PYXTOT is replaced.

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

CC    CALL RYDATA

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
      IF(MINT(44).EQ.4) CALL FRRYXTO

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

      SUBROUTINE FRRYXTO

C...This routine is borrowed from PYTHIA's PYXTOT.  Modifications:
C... The Block and Cahn fit No.2 for slope parameters are changed according
C    to their paper in Physics Simulations at High Energy
C    (edited by V.Barger,etc), which gives better fit.
C...  MSTU(31)=6: Block-Cahn fit 8 for total xsection, and fit 1 for slope;
C...  The pion scaling factor is slightly modified to better accomodate
C...  the low energy data.
C...  Kaons added.
C
C...Parametrizes total, double diffractive, single diffractive and
C...elastic cross-sections for different energies and beams.
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/RYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/RYINT1/MINT(400),VINT(400)
      COMMON/RYINT5/NGEN(0:200,3),XSEC(0:200,3)

      COMMON/FRTORYXTO/KF0(2)

      SAVE /LUDAT1/,/RYPARS/,/RYINT1/,/RYINT5/,/FRTORYXTO/
      DIMENSION BCS(5,8),BCB(2,5),BCC(3),SCALE(2)

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
     2  9.92, 0.27, 0.013, 18.9, 1.93/
CC   2  9.92, -0.027, 0.013, 18.9, 1.07/
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
      IF(MSTP(31).GE.4.and.MSTP(31).LE.5) NFIT=2
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

C...Rescale for MESONS.
      SCALPI = (2./3.-1.13/VINT(1))
      SCALK = (2./3.-3.27/VINT(1))

      DO 110 LO = 1, 2
      IF(IABS(KF0(LO)).EQ.321) THEN
      SCALE(LO) = SCALK
      ELSEIF(IABS(KF0(LO)).LE.999) THEN
      SCALE(LO) = SCALPI
      ELSE
      SCALE(LO) = 1.0
      ENDIF
110   CONTINUE

      SIGMA=SIGMA *SCALE(1)*SCALE(2)
      SIGDD=SIGDD *SCALE(1)*SCALE(2)
      SIGSD=SIGSD *SCALE(1)*SCALE(2)
      SIGEL=SIGEL *SCALE(1)*SCALE(2)
      SIGND=SIGND *SCALE(1)*SCALE(2)

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
C********************************************************************

        SUBROUTINE FRCHKEP(IQ)

C... CHECK THE TOTAL CHARGE, ENERGY, MOMENTUM CONSERVATION
C... IQ=0: Just check the sums and take no further steps
C...   =1: Monitor the number of errors and signal (IOP(16)) for printout
C...       via FRMGOUT.

      PARAMETER (KSZJ=4000,KSZ1=20, KSZ2=300)
      COMMON/FRINTN0/PLI0(2,4),AOP(KSZ1),IOP(KSZ1),NFR(KSZ1)

      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      integer N,NPAD,K
      double precision P,V
      SAVE /PYJETS/

      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)

      Data ifst /0/
      Save Etot, Ptot, Cgtot
      SAVE /FRINTN0/,/LUDAT1/

      if(Ifst.eq.0) then
C  TOTAL Beam energy, momentum:
        EBM = 0.5*(PLI0(1,4)+PLI0(1,3))* IOP(3)
        PBM = 0.5*(PLI0(1,4)-PLI0(1,3))* IOP(3)
C  TOTAL target energy, momentum:
        ETG = 0.5*(PLI0(2,4)+PLI0(2,3))* IOP(5)
        PTG = 0.5*(PLI0(2,4)-PLI0(2,3))* IOP(5)

        Etot = EBM + ETG
        Ptot = PBM + PTG
        Cgtot = IOP(4)+ IOP(6)
        Ifst=1
      endif

      charg =0.0
      EE=0.0
      PPz=0.0
      PPx=0.0
      PPy=0.0
      DO 100 J=1, N
        IF(ABS(K(J,2)).GE.10000) then
          charg = charg+ ABS(K(J,2))-10000
          EE = EE+ P(j,4)
          PPz = PPz+ P(j,3)
          PPy = PPy+ P(j,2)
          PPx = PPx+ P(j,1)
        elseif(K(J,1).ge.1.and.K(J,1).le.5) then
          charg = charg+ PLU(j,6)
          EE = EE+ P(j,4)
          PPz = PPz+ P(j,3)
          PPy = PPy+ P(j,2)
          PPx = PPx+ P(j,1)
        endif
100   CONTINUE

      Ifel=0
      if(abs(PPx).gt.0.5.or.abs(PPy).gt.0.5) Ifel =1
      if(abs(PPZ-Ptot).gt.MAX(0.01*Ptot,0.5)) Ifel =1
      if(abs(EE-Etot).gt.MAX(0.01*Etot,0.5)) Ifel =1
      if(abs(charg-Cgtot).gt.0.01) Ifel =1

      if(ifel.eq.1) then
        write(MSTU(11), 1000)  NFR(1)
        IF(IQ.EQ.1)
     >     CALL FRMGOUT(-1,0,'Charge or energy non-conservation',
     >                   0.,0.,0.,0.,0.)
        write(MSTU(11), 1010) Ptot, Etot, Cgtot
        write(MSTU(11), 1020) PPx,PPy, PPz, EE, charg
      endif

1000  format( /,'???????????????????????????????????????????'
     >      /,'  Charge or energy non-conservation at event:', I6 )
1010  format('  Original Pz, E, Charge: ', 26x, 2E13.4, F6.1)
1020  format(' Total Px, Py, Pz, E, Cg: ', 4E13.4, F6.1 )
      RETURN
      END

C*********************************************************************

C*********************************************************************
C********************************* DATA FRDATA ***********************

      BLOCK DATA FRDATA
      PARAMETER (KSZ1=20)
      CHARACTER*4 PACD
      COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)
      COMMON/FRCODES/IPT(2),PACD(27),NNUC(27),NPROT(27),KCD(27)
     >           ,RO1(27,2),EXMA(9,2)
      SAVE /FRPARA1/,/FRCODES/

      DATA KFR/1,1,0,0,0,0,1,1,1,1, 4,2,0,0,0, 5*0/
c      DATA KFR/1,1,0,0,0,0,0,1,1,1, 4,2,0,0,0, 5*0/

c.jochen
      DATA VFR/0.,0.2,0.8,0.2,0.1, 0.01, 1.3,.75,.75,0.,
     >         0.,1.0,0.1,.6,.3, .5,4*0./
c       DATA VFR/0.,0.2,0.8,0.2,0.1, 0.01,1.0,.75,.75,0.,
c     >         0.,1.0,0.167,.333,.5, .5,4*0./

C.....The following are particles in store.  In particulare RO1 are the
C.....parameters for nuclei density.  For A<,=16, shell model harmonic
C.....oscilator density, RO1(j,1) gives the nuclear root-mean-square-
C.....(charge) radius.  For A>16 RO1(j,2) gives the two parameters r0 and C
C.....to the Woods-Saxon density.
C.....EXMA(J,1) AND EXMA(J,2) correspond to the "minimum excitation mass"
C.....and "diffractive mass" respectively.

      DATA PACD/'NEW1','NEW2','PI+ ','PI- ','K+  ','K-  ','N   ','P   '
     >         ,'PBAR','D   ','HE  ','BE  ','B   ','C   ','O   ','AL  '
     >         ,'SI  ','S   ','AR  ','CA  ','CU  ','AG  ','XE  ','W   '
     >         ,'AU  ','PB  ','U   '/
      DATA NNUC/  1,     1,     1,     1,     1,     1,     1,     1,
     >            1,     2,     4,     9,    11,    12,    16,    27,
     >           28,    32,    40,    40,    64,   108,   131,   184,
     >          197,   207,   238/
      DATA NPROT/ 0,    0,      1,    -1,     1,    -1,     0,     1,
     >           -1,    1,      2,     4,     5,     6,     8,    13,
     >           14,   16,     18,    20,    29,    47,    54,    74,
     >           79,    82,    92 /
c....see also MONTE CARLO PARTICLE NUMBERING SCHEME in Phys.Lett. B239 (1990)
c....(PDG) on page 67.III
      DATA KCD /  0,    0,    211,  -211,   321,   -321, 2112,  2212,
     >        -2212,  18*0 /

C.....Source for RO1: RO(j,1) for A<17 is from ref. BJ; and RO(j,1-2) for
C.....A>17 are taken from FRITIOF 6.0, where RO(j,1)=r0=1.16*(1.-1.16/A**(2/3)):
c.....nuclear density distribution:
c.....ro1(i,1) gives the radius parameter r_0, ro1(i,2) the edge
c.....thickness parameter in the Woods-Saxon density. All numbers are
c.....given in fermi.

      DATA RO1 / 0.,    0.,    0.,    0.,     0.,    0.,    0.,    0.,
     1           0., 2.095,  1.74, 2.519,  2.37, 2.446, 2.724,   1.01,
     1        1.014, 1.027, 1.045, 1.045, 1.076, 1.101, 1.108,  1.118,
     1        1.120, 1.122, 1.125,
     2           0.,    0.,    0.,    0.,     0.,    0.,    0.,    0.,
     2           0.,    0.,    0.,    0.,     0.,    0.,    0., 0.478,
     2        0.480, 0.490, 0.490, 0.490,  0.490, 0.495,  0.52, 0.530,
     2        0.540, 0.545,  0.55 /

C.....RO1 follows are taken from experimental measurements from ref. BJ,
C.....and some unfound in the book are hand-picked. Not going to be used.
C
C      DATA RO1 / 0.,    0.,    0.,    0.,     0.,    0.,    0.,    0.,
C     1           0., 2.095,  1.74, 2.519,  2.37, 2.446, 2.724,  0.947,
C     1        1.035, 1.016,  1.00, 1.023, 1.068, 1.119, 1.114,  1.116,
C     1        1.109, 1.121, 1.123,
C     2           0.,    0.,    0.,    0.,     0.,    0.,    0.,    0.,
C     2           0.,    0.,    0.,    0.,     0.,    0.,    0., 0.569,
C     2        0.537, 0.540, 0.543,  0.55,  0.579, 0.523,  0.52, 0.525,
C     2        0.535, 0.535,  0.55 /

      DATA EXMA/0.00,  0.00,   0.14,   0.14,  0.50,  0.50,  0.94,  0.94,
     1          0.94,
     2          0.00,  0.00,   0.40,   0.40,  0.75,  0.75,  1.20,  1.20,
     2          1.20 /

C.....Reference. BJ:  R.C. Barrett and D.F.Jackson,
C.....Nuclear Sizes and Structure. However Ca, S, Ar, Xe, W can not be found
C.....in the book.

      END


CC****************************** END FRDATA ***********************

C******************************************************************
C...........Main switches and parameters...........................

C\item[KFR(1)] (D=1) Fragmentation
C       \item[=0] Off.
C       \item[=1] On.
C\item[KFR(2)] (D=1) Multiple gluon emission (dipole radiation)
C       \item[=0] Off.
C       \item[=1] On.
C\item[KFR(3)] (D=0) Event selection for collisions with a nucleus
C       \item[=0] Generate minimum bias events (all interactions recorded).
C       \item[=1] Generate only events with all projectile nucleons
C                  participated.
C       \item[=2] Generate only events with impact parameter between
C                  $b_{min}$=VFR(1) and $b_{max}$=VFR(2).
C       \item[=3] Apply both requirements in 1 and 2.
C\item[KFR(4)] (D=1) Fermi motion in nuclei
C       \item[=0] Neglected.
C       \item[=1] Included.
C\item[KFR(5)] (D=0) Nucleon-nucleon overlap function
C       \item[=0] Eikonal.
C       \item[=1] Gaussian.
C       \item[=2] Gray disc.
C\item[KFR(6)] (D=2) Target Nucleus deformation
C       \item[=0] No deformation.
C       \item[=1] Deformed target nucleus.
C       \item[=2] Apply deformation only if the target atomic number $A\geq 80$.
C\item[KFR(7)] (D=1) Rutherford parton scattering processes
C       \item[=0] Off.
C       \item[=1] On. Here only the hardest RPS is used in FRITIOF.
C       \item[=2] On. The full multiple hard scattering scenario of PYTHIA
C               is used.
C\item[KFR(8)](D=1) Hard gluons cause a corner (soft gluon kink) on the string
C       \item[=0] No kink is formed.
C       \item[=1] Gluon kink is formed.
C\item[KFR(9)] (D=1) `Drowning' of Rutherford parton scattering
C       \item[=0] Off. Accept all RPS events.
C       \item[=-1] On. Throw away the drowned RPS event completely
C                  and replace it by a purely soft event.
C       \item[=1] As in -1, but the
C            transverse momentum transfer of the soft collision is superimposed
C            by the $q_T$ of the drowned RPS.
C\item[KFR(10)] (D=1)
C        SRM parameters in RPS events: $\mu_1=\mu_0/r$, $\mu_2=\mu_0/(1-r)$
C   \item[=0] $\mu$ remains the same as in a soft event: $\mu_1=\mu_2=\mu_0$.
C      \item[=1] $r$\,=\,VFR(16).
C      \item[=2] $r$ takes a uniform distribution in (0,1).
C\item[KFR(11)] (D=4)
C      Write out of a message when the arguments in FREVENT is changed.
C      \item[=-1] Write it out every time the change occurs.
C      \item[=$n$ ($n\geq 0$)] The write out is limited to $n$ times.
C\item[KFR(12)] (D=2)
C        Set up of the dipole cascade and string fragmentation parameters.
C      \item[=0] No set up.  The default values are used.
C      \item[=1] Set to the values optimised by
C       OPAL collaboration \cite{opal}:
C      PARA(1)=0.20, PARA(3)=1.0, PARJ(21)=0.37, PARJ(41)=0.18, PARJ(42)=0.34.
C      \item[=2] Set to the values optimised by DELPHI collaboration:
C      PARA(1)=0.22, PARA(3)=0.6, PARJ(21)=0.405, PARJ(41)=0.23, PARJ(42)=0.34.
C\item[KFR(13)] (D=0)
C       Compresses the event record to save space in LUJETS.  This switch
C       is particularly needed for heavy ion collisions at high energy
C       where LUJETS must be compressed before it gets overfilled.
C      \item[=0] Do not compress LUJETS.
C      \item[=1-3] LUEDIT(KFR(13)) is called and LUJETS is compressed.
C        Specifically, for KFR(13)=1 fragmented jets and
C        decayed particles are removed, for KFR(13)=2 neutrinos and
C        unknown particles are also removed, and for KFR(13)=3
C        neutral particles are further excluded.
C      \item[=4] A dummy subroutine FREDITD() is provided as an interface
C        in which a user may write his own special purpose codes to edit
C        and compress LUJETS.
C\item[KFR(14)] (D=0)
C       If set to 1, the outcome of each event will be checked for
C       charge and energy-momentum conservation.
C%%
C\item[VFR(1)] (D=0.0 fm)
C       Minimum impact parameter for options KFR(3)=2 or 3.
C\item[VFR(2)] (D=0.2 fm)
C       Maximum impact parameter for options KFR(3)=2 or 3.
C\item[VFR(3)] (D=0.8 fm)
C  The minimum allowable distance $R_{min}$ between nucleons in a nucleus.
C\item[VFR(4-5)] (D=0.2, 0.1)
C  Dipole and quadrupole deformation coefficients for deformed target nucleus.
C\item[VFR(6)] (D=0.01 GeV$^2/c^2$)
C       The $<Q_T^2>$ for the Gaussian distribution of soft transverse
C       momentum transfer.
C\item[VFR(7)] (D=0.30 GeV$^2/c^2$)
C       The $<Q^2_{2T}>$ for the Gaussian distribution of primordial transverse
C       momenta on the string ends.
C\item[VFR(8)] (D=0.75 GeV)
C       Soft radiation coherence parameter $\mu_0$ for projectile hadron or
C       nucleon.
C\item[VFR(9)] (D=0.75 GeV)
C   Soft radiation coherence parameter $\mu_0$ for target hadron or nucleon.
C\item[VFR(10-11)] (D=0.0, 0.0 mb)
C       Projectile-target nucleon total and elastic cross
C       sections, respectively.
C       By default, they are taken from the parametrization
C       of Block and Cahn \cite{block} (MSTP(31)=5 in PYTHIA).
C       The meson-nucleon cross sections are obtained
C       simply by scaling down the Block-Cahn fit. The scale factor is
C       $(2/3-a/\sqrt s)$, where $a=1.13$ GeV for pions and $a=3.27$ GeV for
C       kaons are chosen to reproduce the low energy experimental data. For
C       all the other baryons, it is treated as a pion if it is a meson
C       and it is treated as a proton if it is a baryon.
C       User may override the default by setting VFR(17-18) to positive values.
C       However, the user assigned cross sections will only affect the
C       N-N interaction probability in nucleus collisions. The
C       probability for Rutherford parton scattering is not affected.
C\item[VFR(12)] (D=1.0 GeV/$c$)
C       The $q_{Tmin}$ for Rutherford parton scattering.
C\item[VFR(13-15)] (D=1/6, 1/3, 1/2)
C       The probabilities for assigning various spins and flavours to the
C      diquark end of the string.  For example in a proton, VFR(13-15) are the
C       probabilities of finding a $ud$ diquark of spin 1,
C       a $uu$ diquark of spin 1, and a $ud$ diquark of spin 0, respectively.
C\item[VFR(16)] (D=0.5)
C       The fraction $r$ in option KFR(10)=1.
C
CC**************************************************************************
C********************************* END OF FRITIOF PACKAGE ******************



C**************************** FREVENT ***********************************

      SUBROUTINE FREVENT(FRAME,BEAM,TARGET,WIN)

C...This is the main routine that initializes and call FRINGEB to
C...administrate the event generation

      PARAMETER (KSZJ=4000,KSZ1=20,KSZ2=300)
      CHARACTER*(*) FRAME,BEAM,TARGET
      CHARACTER CFRAME*4,CBEAM*4,CTARG*4,CMEM*4
      CHARACTER PARTIC*4,PACD*4
      COMMON/FRCODES/IPT(2),PACD(27),NNUC(27),NPROT(27),KCD(27)
     >           ,RO1(27,2),EXMA(9,2)
      COMMON/FRINTN0/PLI0(2,4),AOP(KSZ1),IOP(KSZ1),NFR(KSZ1)
      COMMON/FRINTN3/IDN(2,KSZ2),FMN(2,KSZ2),NUC(2,3000)
      COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/ARDAT1/PARA(40),MSTA(40)

      DIMENSION CMEM(3)
      SAVE CMEM,WINM,IFST
      SAVE /FRINTN0/,/FRINTN3/,/FRPARA1/,/LUDAT1/,/ARDAT1/
      DATA WINM,IFST /0.,0/

C.....check particles and frame, set initial momenta and write header

      MSTU(21)=1

       CBEAM=BEAM
       CTARG=TARGET
       CFRAME=FRAME

       CALL FRUPCAS(CFRAME)
       CALL FRUPCAS(CBEAM)
       CALL FRUPCAS(CTARG)
       INEW = 0
c......cmem(i)='   ', i=1..3

*5.7.97 M.E.
       inew=1

*       IF((CFRAME.NE.CMEM(1)).or.(CBEAM.NE.CMEM(2)).or.
*     >    (CTARG.NE.CMEM(3)).or.(ABS(WIN-WINM).GT.0.001)) INEW=1

       IF(INEW.GT.0) THEN
*          write(*,*)'in initialization'
         CMEM(1) = CFRAME
         CMEM(2) = CBEAM
         CMEM(3) = CTARG
         WINM = WIN
         CALL FRINITA(CFRAME,CBEAM,CTARG,WIN)

         IF(KSZJ.NE.MSTU(4)) THEN
         CALL FRMGOUT(0,0,'LUJETS not compatible in FRITIOF and JETSET'
     >                ,float(KSZJ),float(MSTU(4)),0.,0.,0.)
         ENDIF
      ENDIF

      IF(IFST.EQ.0) THEN
      DO 12 LO=1, KSZ1
 12     NFR(LO)= 0
C.....     set some control parameters for Ariadne and jetset:
      MSTA(7) = MSTU(11)
      MSTA(8) = MSTU(11)
      MSTA(9) = 0
      MSTA(14) = 0
C.....extended partons made massless, for compatibility with AR3.03
      MSTA(31) = 0
        IF(KFR(12).EQ.1) THEN
C             Fragmentation parameters set to the Opal optimized values:
        CALL ARTUNE('OPAL')
        ELSEIF(KFR(12).GE.2) THEN
C                               DELPHI optimized values
        CALL ARTUNE('DELPHI')
        ENDIF
      CALL ARINIT('ARIADNE')
      IFST= 1
      ENDIF

C.....administrate one event.........................................

      CALL FRINGEB

      RETURN
      END

C**************************** END FRITIOF *******************************
