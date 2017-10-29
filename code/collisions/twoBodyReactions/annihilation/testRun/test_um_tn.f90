      program test_um_tn

      use inputGeneral

      implicit none


      COMMON /MES/SYMB(44),SMM(44),ISM(44),ISM3(44),IY(44),&
     &            JSP(44),IG(44),IP(44),IC(44),ICH(44),SYC(44)
      CHARACTER*4 SYMB,SYC
      REAL SMM
      INTEGER ISM,ISM3,IY,JSP,IG,IP,IC,ICH

      COMMON /CHGR/ICG(6,6000,3),NSUPL(3)
      integer ICG,NSUPL
      COMMON /WFFW/FWP(6000),FWP_N(6000),FWN(6000)
      real FWP,FWP_N,FWN
      SAVE /CHGR/, /WFFW/

      integer j,i,NUMERO,L1,L2,L3,L4
      real PLAB

      call readInputGeneral


      open(1,file='chk_p.out',status='unknown')
      open(2,file='chk_n.out',status='unknown')


      do j=1,1

        if (j.eq.1) PLAB=0.
        if (j.eq.2) PLAB=1.
        if (j.eq.3) PLAB=10.
        if (j.eq.4) PLAB=5.
        if (j.eq.5) PLAB=0.
        
        call UM_TN(PLAB)

        write(*,*)' Number of pbarp and pbarn annihilation channels :',&
     &            NSUPL(1),NSUPL(3)

                write(1,*)' j= ', j
 
               WRITE(1,703)
 703  FORMAT(///73('*')/5X,'COMPLETE LIST OF ANNIHILATION CHANNELS',&
     &'  WITH  P R O T O N (%)'/&
     &73('*')/1X,'NUMBER',9X,'IDENTIFIERS',12X,'S Y M B O L S',11X,&
     &'WITH P')
                NUMERO=0
                DO I=1,NSUPL(1)
                  IF (FWP(I).GT.9.999999E-05) THEN
                    NUMERO=NUMERO+1 
                    WRITE(1,405)NUMERO,(ICG(L1,I,1),L1=1,6),&
     &                           (SYMB(ICG(L2,I,1)),L2=1,6),FWP(I)
 405  FORMAT(' ',I5,' ',6(1X,I3),' ',6(1X,A4),' ',F8.4)
                  END IF 
                END DO 

                write(2,*)' j= ', j
                WRITE(2,406)
 406  FORMAT(///73('*')/5X,'COMPLETE LIST OF ANNIHILATION CHANNELS',&
     &'  WITH  N E U T R O N (%)'/&
     &73('*')/1X,'NUMBER',9X,&
     &'IDENTIFIERS',12X,'S Y M B O L S',11X,'WITH N')
                NUMERO=0
                DO I=1,NSUPL(3)
                  IF (FWN(I).GT.9.999999E-05) THEN
                    NUMERO=NUMERO+1 
                    WRITE(2,405)NUMERO,(ICG(L3,I,3),L3=1,6),&
     &                           (SYMB(ICG(L4,I,3)),L4=1,6),FWN(I)
                  END IF 
                END DO 

      end do

      close(1)
      close(2)

      end


