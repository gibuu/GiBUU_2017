c===================================================================
c HandleSFREP.F: Kai Gallmeister, 6.5.2004...30.1.2006
c
c Routines to handle the common block /DataSFREP/.
c Needed for calculation of production ans formation times from the
c (modified) JETSET-Fragmentation routine "PYSTRF.FT.F".
c
c 30.9.2004:
c 30.1.2006: SFREP_Reset optimized
c
c===================================================================

c===================================================================
c Reset the entries in /DataSFREP/ for lines iline1..iline1.
c
c if iline2<iline2: all lines 1..4000 are reset
c
c
c ***** used by PYSTRFT.FT.F *****
c
      subroutine SFREP_Reset(iLine1,iLine2)
      IMPLICIT NONE
      integer iLine1,iLine2

      integer SFR_NvecMax
      parameter (SFR_NvecMax=100) ! former value: 40

      common /DataSFREP/
     $     SFR_Gamma(4000),     ! reported Gamma
     $     SFR_Z(4000),         ! reported Z 
     $     SFR_mT2(4000),       ! reported mT2
     $     SFR_Flag(4000),      ! Flag: JS*(1..4)
     $     SFR_Flag1(4000),     ! Flag-Counter: added 1,if Flag is set
     $     SFR_REV(4000),       ! reverse ordering (unused)
     $     SFR_REV2(4000),      ! reverse ordering (unused)
     $     SFR_Nvec(4000),      ! number of vectors per particle
     $     SFR_Coeff(4000,SFR_NvecMax), ! coeff of vector
     $     SFR_Vec(4,4000,SFR_NvecMax), ! vector
     $     SFR_VecID(4000,SFR_NvecMax), ! ID of vector (cf. SFREPS)
     $     SFR_JJ(2,4000),      ! used for string-break ordering
     $     SFR_KK(2,4000),      !     -"-
     $     SFR_ID0(4000),       !     -"-
     $     SFR_DHM(4,4000),     ! reported coeffs for M-calc
     $     SFR_DHG(4,4000),     ! reported coeffs for Gamma-calc
     $     SFR_FB(4000),        ! reported parameter FB (???)
     $     SFR_GamIN(2,4000),   ! reported IN(1),IN(2)
     $     SFR_GamC(20,4000),   ! reported Coeffs for Gamma
     $     SFR_Ntot,            ! #(used entries) = #particles
     $     SFR_Nused(2)         ! indizes of lines which has been used
                                !                (min,max)

      integer SFR_Ntot
      integer SFR_Nused
      double precision SFR_Gamma,SFR_Z, SFR_mT2
      integer SFR_Flag, SFR_Flag1, SFR_REV, SFR_REV2
      integer SFR_Nvec, SFR_VecID
      double precision SFR_Coeff, SFR_Vec
      integer SFR_JJ,SFR_KK,SFR_ID0
      double precision SFR_DHM, SFR_DHG
      double precision SFR_FB
      integer SFR_GamIN
      double precision SFR_GamC

      save /DataSFREP/
      

      integer iL,iL1,iL2,j

c      write(*,*) '==> SFREP_Reset(',iLine1,iLine2,')'

      if (iLine2.lt.iLine1) then
         iL1 = 1
         iL2 = 4000
      else
         iL1 = iLine1
         iL2 = iLine2
      endif

      call SFREP_CalcToReset(iL1,iL2)

c      write(*,*) '==> SFREP_Reset(',iLine1,iLine2,')',iL1,iL2

      do iL=iL1,iL2
         SFR_Gamma(iL)= 0d0
         SFR_Z(iL)    = 0d0
         SFR_mT2(iL)  = 0d0
         SFR_Flag(iL) = 0
         SFR_Flag1(iL)= 0
         SFR_REV(iL)  = 0
         SFR_REV2(iL) = 0
         SFR_NVec(iL) = 0
         do j=1,2
            SFR_KK(j,iL) = 0
            SFR_JJ(j,iL) = 0
         enddo
         SFR_ID0(iL)  = 0
         do j=1,4
            SFR_DHM(j,iL) = 0d0
            SFR_DHG(j,iL) = 0d0
         enddo
         SFR_FB(iL)   = 0d0

         SFR_GamIN(1,iL) = 0
         SFR_GamIN(2,iL) = 0
         do j=1,10
            SFR_GamC(j,iL) = 0d0
         enddo
      enddo

      end
c====================================================
c Set the number of produced particles

      subroutine SFREP_SetN()
      IMPLICIT NONE

      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      integer N,NPAD,K
      double precision P,V
      SAVE /PYJETS/

      integer SFR_NvecMax
      parameter (SFR_NvecMax=100) ! former value: 40

      common /DataSFREP/
     $     SFR_Gamma(4000),     ! reported Gamma
     $     SFR_Z(4000),         ! reported Z 
     $     SFR_mT2(4000),       ! reported mT2
     $     SFR_Flag(4000),      ! Flag: JS*(1..4)
     $     SFR_Flag1(4000),     ! Flag-Counter: added 1,if Flag is set
     $     SFR_REV(4000),       ! reverse ordering (unused)
     $     SFR_REV2(4000),      ! reverse ordering (unused)
     $     SFR_Nvec(4000),      ! number of vectors per particle
     $     SFR_Coeff(4000,SFR_NvecMax), ! coeff of vector
     $     SFR_Vec(4,4000,SFR_NvecMax), ! vector
     $     SFR_VecID(4000,SFR_NvecMax), ! ID of vector (cf. SFREPS)
     $     SFR_JJ(2,4000),      ! used for string-break ordering
     $     SFR_KK(2,4000),      !     -"-
     $     SFR_ID0(4000),       !     -"-
     $     SFR_DHM(4,4000),     ! reported coeffs for M-calc
     $     SFR_DHG(4,4000),     ! reported coeffs for Gamma-calc
     $     SFR_FB(4000),        ! reported parameter FB (???)
     $     SFR_GamIN(2,4000),   ! reported IN(1),IN(2)
     $     SFR_GamC(20,4000),   ! reported Coeffs for Gamma
     $     SFR_Ntot,            ! #(used entries) = #particles
     $     SFR_Nused(2)         ! indizes of lines which has been used
                                !                (min,max)

      integer SFR_Ntot
      integer SFR_Nused
      double precision SFR_Gamma,SFR_Z, SFR_mT2
      integer SFR_Flag, SFR_Flag1, SFR_REV, SFR_REV2
      integer SFR_Nvec, SFR_VecID
      double precision SFR_Coeff, SFR_Vec
      integer SFR_JJ,SFR_KK,SFR_ID0
      double precision SFR_DHM, SFR_DHG
      double precision SFR_FB
      integer SFR_GamIN
      double precision SFR_GamC

      save /DataSFREP/
      

      SFR_Ntot = N
      end

c====================================================
c Set the value SFR_Gamma(iLine) to G.

      subroutine SFREP_SetGamma(iLine,G)
      IMPLICIT NONE
      integer iLine
      double precision G

      integer SFR_NvecMax
      parameter (SFR_NvecMax=100) ! former value: 40

      common /DataSFREP/
     $     SFR_Gamma(4000),     ! reported Gamma
     $     SFR_Z(4000),         ! reported Z 
     $     SFR_mT2(4000),       ! reported mT2
     $     SFR_Flag(4000),      ! Flag: JS*(1..4)
     $     SFR_Flag1(4000),     ! Flag-Counter: added 1,if Flag is set
     $     SFR_REV(4000),       ! reverse ordering (unused)
     $     SFR_REV2(4000),      ! reverse ordering (unused)
     $     SFR_Nvec(4000),      ! number of vectors per particle
     $     SFR_Coeff(4000,SFR_NvecMax), ! coeff of vector
     $     SFR_Vec(4,4000,SFR_NvecMax), ! vector
     $     SFR_VecID(4000,SFR_NvecMax), ! ID of vector (cf. SFREPS)
     $     SFR_JJ(2,4000),      ! used for string-break ordering
     $     SFR_KK(2,4000),      !     -"-
     $     SFR_ID0(4000),       !     -"-
     $     SFR_DHM(4,4000),     ! reported coeffs for M-calc
     $     SFR_DHG(4,4000),     ! reported coeffs for Gamma-calc
     $     SFR_FB(4000),        ! reported parameter FB (???)
     $     SFR_GamIN(2,4000),   ! reported IN(1),IN(2)
     $     SFR_GamC(20,4000),   ! reported Coeffs for Gamma
     $     SFR_Ntot,            ! #(used entries) = #particles
     $     SFR_Nused(2)         ! indizes of lines which has been used
                                !                (min,max)

      integer SFR_Ntot
      integer SFR_Nused
      double precision SFR_Gamma,SFR_Z, SFR_mT2
      integer SFR_Flag, SFR_Flag1, SFR_REV, SFR_REV2
      integer SFR_Nvec, SFR_VecID
      double precision SFR_Coeff, SFR_Vec
      integer SFR_JJ,SFR_KK,SFR_ID0
      double precision SFR_DHM, SFR_DHG
      double precision SFR_FB
      integer SFR_GamIN
      double precision SFR_GamC

      save /DataSFREP/
      

      call SFREP_SetUsed(iLine)

      SFR_Gamma(iLine) = G
      end

c====================================================
c Set the values SFR_Gamma(iLine),SFR_Z(iLine),SFR_mT2(iLine) to G,z,m.
c
c ***** used by PYSTRFT.FT.F *****
c
      subroutine SFREP_SetGZM(iLine,G,z,m)
      IMPLICIT NONE
      integer iLine
      double precision G,z,m

      integer SFR_NvecMax
      parameter (SFR_NvecMax=100) ! former value: 40

      common /DataSFREP/
     $     SFR_Gamma(4000),     ! reported Gamma
     $     SFR_Z(4000),         ! reported Z 
     $     SFR_mT2(4000),       ! reported mT2
     $     SFR_Flag(4000),      ! Flag: JS*(1..4)
     $     SFR_Flag1(4000),     ! Flag-Counter: added 1,if Flag is set
     $     SFR_REV(4000),       ! reverse ordering (unused)
     $     SFR_REV2(4000),      ! reverse ordering (unused)
     $     SFR_Nvec(4000),      ! number of vectors per particle
     $     SFR_Coeff(4000,SFR_NvecMax), ! coeff of vector
     $     SFR_Vec(4,4000,SFR_NvecMax), ! vector
     $     SFR_VecID(4000,SFR_NvecMax), ! ID of vector (cf. SFREPS)
     $     SFR_JJ(2,4000),      ! used for string-break ordering
     $     SFR_KK(2,4000),      !     -"-
     $     SFR_ID0(4000),       !     -"-
     $     SFR_DHM(4,4000),     ! reported coeffs for M-calc
     $     SFR_DHG(4,4000),     ! reported coeffs for Gamma-calc
     $     SFR_FB(4000),        ! reported parameter FB (???)
     $     SFR_GamIN(2,4000),   ! reported IN(1),IN(2)
     $     SFR_GamC(20,4000),   ! reported Coeffs for Gamma
     $     SFR_Ntot,            ! #(used entries) = #particles
     $     SFR_Nused(2)         ! indizes of lines which has been used
                                !                (min,max)

      integer SFR_Ntot
      integer SFR_Nused
      double precision SFR_Gamma,SFR_Z, SFR_mT2
      integer SFR_Flag, SFR_Flag1, SFR_REV, SFR_REV2
      integer SFR_Nvec, SFR_VecID
      double precision SFR_Coeff, SFR_Vec
      integer SFR_JJ,SFR_KK,SFR_ID0
      double precision SFR_DHM, SFR_DHG
      double precision SFR_FB
      integer SFR_GamIN
      double precision SFR_GamC

      save /DataSFREP/
      

c      write(*,*) '==> SFREP_SetGZM(',iLine,G,z,m,')'

      call SFREP_SetUsed(iLine)

      SFR_Gamma(iLine) = G
      SFR_Z(iLine) = z
      SFR_mT2(iLine) = m
      end

c====================================================
c Set the values SFR_DHM()
c
c ***** used by PYSTRFT.FT.F *****
c
      subroutine SFREP_SetDHM(iLine,DHM)
      IMPLICIT NONE
      integer iLine
      double precision DHM(4)

      integer SFR_NvecMax
      parameter (SFR_NvecMax=100) ! former value: 40

      common /DataSFREP/
     $     SFR_Gamma(4000),     ! reported Gamma
     $     SFR_Z(4000),         ! reported Z 
     $     SFR_mT2(4000),       ! reported mT2
     $     SFR_Flag(4000),      ! Flag: JS*(1..4)
     $     SFR_Flag1(4000),     ! Flag-Counter: added 1,if Flag is set
     $     SFR_REV(4000),       ! reverse ordering (unused)
     $     SFR_REV2(4000),      ! reverse ordering (unused)
     $     SFR_Nvec(4000),      ! number of vectors per particle
     $     SFR_Coeff(4000,SFR_NvecMax), ! coeff of vector
     $     SFR_Vec(4,4000,SFR_NvecMax), ! vector
     $     SFR_VecID(4000,SFR_NvecMax), ! ID of vector (cf. SFREPS)
     $     SFR_JJ(2,4000),      ! used for string-break ordering
     $     SFR_KK(2,4000),      !     -"-
     $     SFR_ID0(4000),       !     -"-
     $     SFR_DHM(4,4000),     ! reported coeffs for M-calc
     $     SFR_DHG(4,4000),     ! reported coeffs for Gamma-calc
     $     SFR_FB(4000),        ! reported parameter FB (???)
     $     SFR_GamIN(2,4000),   ! reported IN(1),IN(2)
     $     SFR_GamC(20,4000),   ! reported Coeffs for Gamma
     $     SFR_Ntot,            ! #(used entries) = #particles
     $     SFR_Nused(2)         ! indizes of lines which has been used
                                !                (min,max)

      integer SFR_Ntot
      integer SFR_Nused
      double precision SFR_Gamma,SFR_Z, SFR_mT2
      integer SFR_Flag, SFR_Flag1, SFR_REV, SFR_REV2
      integer SFR_Nvec, SFR_VecID
      double precision SFR_Coeff, SFR_Vec
      integer SFR_JJ,SFR_KK,SFR_ID0
      double precision SFR_DHM, SFR_DHG
      double precision SFR_FB
      integer SFR_GamIN
      double precision SFR_GamC

      save /DataSFREP/
      

      integer j

c      write(*,*) '==> SFREP_SetDHM(',iLine,(DHM(j),j=1,4),')'

      call SFREP_SetUsed(iLine)

      do j=1,4
         SFR_DHM(j,iLine) = DHM(j)
      enddo
      end

c====================================================
c Set the values SFR_DHG()
c
c ***** used by PYSTRFT.FT.F *****
c
      subroutine SFREP_SetDHG(iLine,DHG)
      IMPLICIT NONE
      integer iLine
      double precision DHG(4)

      integer SFR_NvecMax
      parameter (SFR_NvecMax=100) ! former value: 40

      common /DataSFREP/
     $     SFR_Gamma(4000),     ! reported Gamma
     $     SFR_Z(4000),         ! reported Z 
     $     SFR_mT2(4000),       ! reported mT2
     $     SFR_Flag(4000),      ! Flag: JS*(1..4)
     $     SFR_Flag1(4000),     ! Flag-Counter: added 1,if Flag is set
     $     SFR_REV(4000),       ! reverse ordering (unused)
     $     SFR_REV2(4000),      ! reverse ordering (unused)
     $     SFR_Nvec(4000),      ! number of vectors per particle
     $     SFR_Coeff(4000,SFR_NvecMax), ! coeff of vector
     $     SFR_Vec(4,4000,SFR_NvecMax), ! vector
     $     SFR_VecID(4000,SFR_NvecMax), ! ID of vector (cf. SFREPS)
     $     SFR_JJ(2,4000),      ! used for string-break ordering
     $     SFR_KK(2,4000),      !     -"-
     $     SFR_ID0(4000),       !     -"-
     $     SFR_DHM(4,4000),     ! reported coeffs for M-calc
     $     SFR_DHG(4,4000),     ! reported coeffs for Gamma-calc
     $     SFR_FB(4000),        ! reported parameter FB (???)
     $     SFR_GamIN(2,4000),   ! reported IN(1),IN(2)
     $     SFR_GamC(20,4000),   ! reported Coeffs for Gamma
     $     SFR_Ntot,            ! #(used entries) = #particles
     $     SFR_Nused(2)         ! indizes of lines which has been used
                                !                (min,max)

      integer SFR_Ntot
      integer SFR_Nused
      double precision SFR_Gamma,SFR_Z, SFR_mT2
      integer SFR_Flag, SFR_Flag1, SFR_REV, SFR_REV2
      integer SFR_Nvec, SFR_VecID
      double precision SFR_Coeff, SFR_Vec
      integer SFR_JJ,SFR_KK,SFR_ID0
      double precision SFR_DHM, SFR_DHG
      double precision SFR_FB
      integer SFR_GamIN
      double precision SFR_GamC

      save /DataSFREP/
      

      integer j

c      write(*,*) '==> SFREP_SetDHG(',iLine,(DHG(j),j=1,4),')'

      call SFREP_SetUsed(iLine)

      do j=1,4
         SFR_DHG(j,iLine) = DHG(j)
      enddo
      end
c====================================================
c Set the value SFR_FB(iLine) to FB.
c
c ***** used by PYSTRFT.FT.F *****
c
      subroutine SFREP_SetFB(iLine,FB)
      IMPLICIT NONE
      integer iLine
      double precision FB

      integer SFR_NvecMax
      parameter (SFR_NvecMax=100) ! former value: 40

      common /DataSFREP/
     $     SFR_Gamma(4000),     ! reported Gamma
     $     SFR_Z(4000),         ! reported Z 
     $     SFR_mT2(4000),       ! reported mT2
     $     SFR_Flag(4000),      ! Flag: JS*(1..4)
     $     SFR_Flag1(4000),     ! Flag-Counter: added 1,if Flag is set
     $     SFR_REV(4000),       ! reverse ordering (unused)
     $     SFR_REV2(4000),      ! reverse ordering (unused)
     $     SFR_Nvec(4000),      ! number of vectors per particle
     $     SFR_Coeff(4000,SFR_NvecMax), ! coeff of vector
     $     SFR_Vec(4,4000,SFR_NvecMax), ! vector
     $     SFR_VecID(4000,SFR_NvecMax), ! ID of vector (cf. SFREPS)
     $     SFR_JJ(2,4000),      ! used for string-break ordering
     $     SFR_KK(2,4000),      !     -"-
     $     SFR_ID0(4000),       !     -"-
     $     SFR_DHM(4,4000),     ! reported coeffs for M-calc
     $     SFR_DHG(4,4000),     ! reported coeffs for Gamma-calc
     $     SFR_FB(4000),        ! reported parameter FB (???)
     $     SFR_GamIN(2,4000),   ! reported IN(1),IN(2)
     $     SFR_GamC(20,4000),   ! reported Coeffs for Gamma
     $     SFR_Ntot,            ! #(used entries) = #particles
     $     SFR_Nused(2)         ! indizes of lines which has been used
                                !                (min,max)

      integer SFR_Ntot
      integer SFR_Nused
      double precision SFR_Gamma,SFR_Z, SFR_mT2
      integer SFR_Flag, SFR_Flag1, SFR_REV, SFR_REV2
      integer SFR_Nvec, SFR_VecID
      double precision SFR_Coeff, SFR_Vec
      integer SFR_JJ,SFR_KK,SFR_ID0
      double precision SFR_DHM, SFR_DHG
      double precision SFR_FB
      integer SFR_GamIN
      double precision SFR_GamC

      save /DataSFREP/
      

c      write(*,*) '==> SFREP_SetFB(',iLine,FB,')'

      call SFREP_SetUsed(iLine)

      SFR_FB(iLine) = FB
      end
c====================================================
c Set the flag SFR_Flag(iLine) to F.
c
c increases an internal counter, which is stored in SFR_Flag1(iLine)
c
c ***** used by PYSTRFT.FT.F *****
c
      subroutine SFREP_SetFlag(iLine,F)
      IMPLICIT NONE
      integer iLine, F

      integer SFR_NvecMax
      parameter (SFR_NvecMax=100) ! former value: 40

      common /DataSFREP/
     $     SFR_Gamma(4000),     ! reported Gamma
     $     SFR_Z(4000),         ! reported Z 
     $     SFR_mT2(4000),       ! reported mT2
     $     SFR_Flag(4000),      ! Flag: JS*(1..4)
     $     SFR_Flag1(4000),     ! Flag-Counter: added 1,if Flag is set
     $     SFR_REV(4000),       ! reverse ordering (unused)
     $     SFR_REV2(4000),      ! reverse ordering (unused)
     $     SFR_Nvec(4000),      ! number of vectors per particle
     $     SFR_Coeff(4000,SFR_NvecMax), ! coeff of vector
     $     SFR_Vec(4,4000,SFR_NvecMax), ! vector
     $     SFR_VecID(4000,SFR_NvecMax), ! ID of vector (cf. SFREPS)
     $     SFR_JJ(2,4000),      ! used for string-break ordering
     $     SFR_KK(2,4000),      !     -"-
     $     SFR_ID0(4000),       !     -"-
     $     SFR_DHM(4,4000),     ! reported coeffs for M-calc
     $     SFR_DHG(4,4000),     ! reported coeffs for Gamma-calc
     $     SFR_FB(4000),        ! reported parameter FB (???)
     $     SFR_GamIN(2,4000),   ! reported IN(1),IN(2)
     $     SFR_GamC(20,4000),   ! reported Coeffs for Gamma
     $     SFR_Ntot,            ! #(used entries) = #particles
     $     SFR_Nused(2)         ! indizes of lines which has been used
                                !                (min,max)

      integer SFR_Ntot
      integer SFR_Nused
      double precision SFR_Gamma,SFR_Z, SFR_mT2
      integer SFR_Flag, SFR_Flag1, SFR_REV, SFR_REV2
      integer SFR_Nvec, SFR_VecID
      double precision SFR_Coeff, SFR_Vec
      integer SFR_JJ,SFR_KK,SFR_ID0
      double precision SFR_DHM, SFR_DHG
      double precision SFR_FB
      integer SFR_GamIN
      double precision SFR_GamC

      save /DataSFREP/
      

      integer FF1
      save FF1

      data FF1 /0/

c      write(*,*) '==> SFREP_SetFlag(',iLine,F,')'

      call SFREP_SetUsed(iLine)

      SFR_Flag(iLine) = F

      FF1 = FF1+1
      SFR_Flag1(iLine) = FF1
      end

c====================================================
c Set the value SFR_REV(iLine) to Q.

      subroutine SFREP_SetREV(iLine,Q)
      IMPLICIT NONE
      integer iLine
      integer Q

      integer SFR_NvecMax
      parameter (SFR_NvecMax=100) ! former value: 40

      common /DataSFREP/
     $     SFR_Gamma(4000),     ! reported Gamma
     $     SFR_Z(4000),         ! reported Z 
     $     SFR_mT2(4000),       ! reported mT2
     $     SFR_Flag(4000),      ! Flag: JS*(1..4)
     $     SFR_Flag1(4000),     ! Flag-Counter: added 1,if Flag is set
     $     SFR_REV(4000),       ! reverse ordering (unused)
     $     SFR_REV2(4000),      ! reverse ordering (unused)
     $     SFR_Nvec(4000),      ! number of vectors per particle
     $     SFR_Coeff(4000,SFR_NvecMax), ! coeff of vector
     $     SFR_Vec(4,4000,SFR_NvecMax), ! vector
     $     SFR_VecID(4000,SFR_NvecMax), ! ID of vector (cf. SFREPS)
     $     SFR_JJ(2,4000),      ! used for string-break ordering
     $     SFR_KK(2,4000),      !     -"-
     $     SFR_ID0(4000),       !     -"-
     $     SFR_DHM(4,4000),     ! reported coeffs for M-calc
     $     SFR_DHG(4,4000),     ! reported coeffs for Gamma-calc
     $     SFR_FB(4000),        ! reported parameter FB (???)
     $     SFR_GamIN(2,4000),   ! reported IN(1),IN(2)
     $     SFR_GamC(20,4000),   ! reported Coeffs for Gamma
     $     SFR_Ntot,            ! #(used entries) = #particles
     $     SFR_Nused(2)         ! indizes of lines which has been used
                                !                (min,max)

      integer SFR_Ntot
      integer SFR_Nused
      double precision SFR_Gamma,SFR_Z, SFR_mT2
      integer SFR_Flag, SFR_Flag1, SFR_REV, SFR_REV2
      integer SFR_Nvec, SFR_VecID
      double precision SFR_Coeff, SFR_Vec
      integer SFR_JJ,SFR_KK,SFR_ID0
      double precision SFR_DHM, SFR_DHG
      double precision SFR_FB
      integer SFR_GamIN
      double precision SFR_GamC

      save /DataSFREP/
      

      call SFREP_SetUsed(iLine)

      SFR_REV(iLine) = Q
      end

c====================================================
c Set the value SFR_REV2(iLine) to Q.
c
c ***** used by PYSTRFT.FT.F *****
c
      subroutine SFREP_SetREV2(iLine,Q)
      IMPLICIT NONE
      integer iLine
      integer Q

      integer SFR_NvecMax
      parameter (SFR_NvecMax=100) ! former value: 40

      common /DataSFREP/
     $     SFR_Gamma(4000),     ! reported Gamma
     $     SFR_Z(4000),         ! reported Z 
     $     SFR_mT2(4000),       ! reported mT2
     $     SFR_Flag(4000),      ! Flag: JS*(1..4)
     $     SFR_Flag1(4000),     ! Flag-Counter: added 1,if Flag is set
     $     SFR_REV(4000),       ! reverse ordering (unused)
     $     SFR_REV2(4000),      ! reverse ordering (unused)
     $     SFR_Nvec(4000),      ! number of vectors per particle
     $     SFR_Coeff(4000,SFR_NvecMax), ! coeff of vector
     $     SFR_Vec(4,4000,SFR_NvecMax), ! vector
     $     SFR_VecID(4000,SFR_NvecMax), ! ID of vector (cf. SFREPS)
     $     SFR_JJ(2,4000),      ! used for string-break ordering
     $     SFR_KK(2,4000),      !     -"-
     $     SFR_ID0(4000),       !     -"-
     $     SFR_DHM(4,4000),     ! reported coeffs for M-calc
     $     SFR_DHG(4,4000),     ! reported coeffs for Gamma-calc
     $     SFR_FB(4000),        ! reported parameter FB (???)
     $     SFR_GamIN(2,4000),   ! reported IN(1),IN(2)
     $     SFR_GamC(20,4000),   ! reported Coeffs for Gamma
     $     SFR_Ntot,            ! #(used entries) = #particles
     $     SFR_Nused(2)         ! indizes of lines which has been used
                                !                (min,max)

      integer SFR_Ntot
      integer SFR_Nused
      double precision SFR_Gamma,SFR_Z, SFR_mT2
      integer SFR_Flag, SFR_Flag1, SFR_REV, SFR_REV2
      integer SFR_Nvec, SFR_VecID
      double precision SFR_Coeff, SFR_Vec
      integer SFR_JJ,SFR_KK,SFR_ID0
      double precision SFR_DHM, SFR_DHG
      double precision SFR_FB
      integer SFR_GamIN
      double precision SFR_GamC

      save /DataSFREP/
      

c      write(*,*) '==> SFREP_SetREV2(',iLine,Q,')'

      call SFREP_SetUsed(iLine)

      SFR_REV2(iLine) = Q
      end

c====================================================
c Reset the Vector iVec of iLine
c
      subroutine SFREP_ResetVec(iLine,iVec)
      IMPLICIT NONE
      integer iLine,iVec

      integer SFR_NvecMax
      parameter (SFR_NvecMax=100) ! former value: 40

      common /DataSFREP/
     $     SFR_Gamma(4000),     ! reported Gamma
     $     SFR_Z(4000),         ! reported Z 
     $     SFR_mT2(4000),       ! reported mT2
     $     SFR_Flag(4000),      ! Flag: JS*(1..4)
     $     SFR_Flag1(4000),     ! Flag-Counter: added 1,if Flag is set
     $     SFR_REV(4000),       ! reverse ordering (unused)
     $     SFR_REV2(4000),      ! reverse ordering (unused)
     $     SFR_Nvec(4000),      ! number of vectors per particle
     $     SFR_Coeff(4000,SFR_NvecMax), ! coeff of vector
     $     SFR_Vec(4,4000,SFR_NvecMax), ! vector
     $     SFR_VecID(4000,SFR_NvecMax), ! ID of vector (cf. SFREPS)
     $     SFR_JJ(2,4000),      ! used for string-break ordering
     $     SFR_KK(2,4000),      !     -"-
     $     SFR_ID0(4000),       !     -"-
     $     SFR_DHM(4,4000),     ! reported coeffs for M-calc
     $     SFR_DHG(4,4000),     ! reported coeffs for Gamma-calc
     $     SFR_FB(4000),        ! reported parameter FB (???)
     $     SFR_GamIN(2,4000),   ! reported IN(1),IN(2)
     $     SFR_GamC(20,4000),   ! reported Coeffs for Gamma
     $     SFR_Ntot,            ! #(used entries) = #particles
     $     SFR_Nused(2)         ! indizes of lines which has been used
                                !                (min,max)

      integer SFR_Ntot
      integer SFR_Nused
      double precision SFR_Gamma,SFR_Z, SFR_mT2
      integer SFR_Flag, SFR_Flag1, SFR_REV, SFR_REV2
      integer SFR_Nvec, SFR_VecID
      double precision SFR_Coeff, SFR_Vec
      integer SFR_JJ,SFR_KK,SFR_ID0
      double precision SFR_DHM, SFR_DHG
      double precision SFR_FB
      integer SFR_GamIN
      double precision SFR_GamC

      save /DataSFREP/
      

      integer j

      if (iVec.gt.SFR_NVecMax) then
         call SFREP_ERROR(1,'SFREP_ResetVec: iVec > NVecMax!')
         return
      endif

      SFR_Coeff(iLine,iVec) = 0d0
      do j=1,4
         SFR_Vec(j,iLine,iVec) = 0d0
      enddo
      SFR_VecID(iLine,iVec) = 0
      
      end

c====================================================
c Add one momentum-vector for the hadron momentum to the list
c
c Adds the momentum vector given by P(iP,..) to the actual list af
c momenta building up the vector of hadron in line number iLine.
c Parameter c gives the coefficient how to multiply this vector in the
c sum.
c
c Every call adds one vector, until the maximum in the list is reached.
c
c ***** used by PYSTRFT.FT.F *****
c
      subroutine SFREP_AddVec(iLine, c, iP)
      IMPLICIT NONE
      integer iLine, iP
      double precision c

      integer SFR_NvecMax
      parameter (SFR_NvecMax=100) ! former value: 40

      common /DataSFREP/
     $     SFR_Gamma(4000),     ! reported Gamma
     $     SFR_Z(4000),         ! reported Z 
     $     SFR_mT2(4000),       ! reported mT2
     $     SFR_Flag(4000),      ! Flag: JS*(1..4)
     $     SFR_Flag1(4000),     ! Flag-Counter: added 1,if Flag is set
     $     SFR_REV(4000),       ! reverse ordering (unused)
     $     SFR_REV2(4000),      ! reverse ordering (unused)
     $     SFR_Nvec(4000),      ! number of vectors per particle
     $     SFR_Coeff(4000,SFR_NvecMax), ! coeff of vector
     $     SFR_Vec(4,4000,SFR_NvecMax), ! vector
     $     SFR_VecID(4000,SFR_NvecMax), ! ID of vector (cf. SFREPS)
     $     SFR_JJ(2,4000),      ! used for string-break ordering
     $     SFR_KK(2,4000),      !     -"-
     $     SFR_ID0(4000),       !     -"-
     $     SFR_DHM(4,4000),     ! reported coeffs for M-calc
     $     SFR_DHG(4,4000),     ! reported coeffs for Gamma-calc
     $     SFR_FB(4000),        ! reported parameter FB (???)
     $     SFR_GamIN(2,4000),   ! reported IN(1),IN(2)
     $     SFR_GamC(20,4000),   ! reported Coeffs for Gamma
     $     SFR_Ntot,            ! #(used entries) = #particles
     $     SFR_Nused(2)         ! indizes of lines which has been used
                                !                (min,max)

      integer SFR_Ntot
      integer SFR_Nused
      double precision SFR_Gamma,SFR_Z, SFR_mT2
      integer SFR_Flag, SFR_Flag1, SFR_REV, SFR_REV2
      integer SFR_Nvec, SFR_VecID
      double precision SFR_Coeff, SFR_Vec
      integer SFR_JJ,SFR_KK,SFR_ID0
      double precision SFR_DHM, SFR_DHG
      double precision SFR_FB
      integer SFR_GamIN
      double precision SFR_GamC

      save /DataSFREP/
      

      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      integer N,NPAD,K
      double precision P,V
      SAVE /PYJETS/

      integer iV,j

c      write(*,*) '==> SFREP_AddVec(',iLine,c,iP,')'

      call SFREP_SetUsed(iLine)

      if (SFR_NVec(iLine).eq.SFR_NvecMax) then
         call SFREP_ERROR(2,'SFREP_AddVec: NVec exhausted!')
         return
      endif

      SFR_NVec(iLine)=SFR_NVec(iLine)+1
      iV = SFR_NVec(iLine)

      SFR_Coeff(iLine,iV) = c
      do j=1,4
         SFR_Vec(j,iLine,iV) = P(iP,j)
      enddo
      SFR_VecID(iLine,iV) = iP
      end

c====================================================
c Add one momentum-vector for the hadron momentum to the list
c
c Adds the momentum vector given by c,v(4),ID to the actual list af
c momenta building up the vector of hadron in line number iLine.
c Parameter c gives the coefficient how to multiply this vector in the
c sum.
c
c Every call adds one vector, until the maximum in the list is reached.
c
      subroutine SFREP_AddVec0(iLine, c, V, ID)
      IMPLICIT NONE
      integer iLine, ID
      double precision c,V(4)

      integer SFR_NvecMax
      parameter (SFR_NvecMax=100) ! former value: 40

      common /DataSFREP/
     $     SFR_Gamma(4000),     ! reported Gamma
     $     SFR_Z(4000),         ! reported Z 
     $     SFR_mT2(4000),       ! reported mT2
     $     SFR_Flag(4000),      ! Flag: JS*(1..4)
     $     SFR_Flag1(4000),     ! Flag-Counter: added 1,if Flag is set
     $     SFR_REV(4000),       ! reverse ordering (unused)
     $     SFR_REV2(4000),      ! reverse ordering (unused)
     $     SFR_Nvec(4000),      ! number of vectors per particle
     $     SFR_Coeff(4000,SFR_NvecMax), ! coeff of vector
     $     SFR_Vec(4,4000,SFR_NvecMax), ! vector
     $     SFR_VecID(4000,SFR_NvecMax), ! ID of vector (cf. SFREPS)
     $     SFR_JJ(2,4000),      ! used for string-break ordering
     $     SFR_KK(2,4000),      !     -"-
     $     SFR_ID0(4000),       !     -"-
     $     SFR_DHM(4,4000),     ! reported coeffs for M-calc
     $     SFR_DHG(4,4000),     ! reported coeffs for Gamma-calc
     $     SFR_FB(4000),        ! reported parameter FB (???)
     $     SFR_GamIN(2,4000),   ! reported IN(1),IN(2)
     $     SFR_GamC(20,4000),   ! reported Coeffs for Gamma
     $     SFR_Ntot,            ! #(used entries) = #particles
     $     SFR_Nused(2)         ! indizes of lines which has been used
                                !                (min,max)

      integer SFR_Ntot
      integer SFR_Nused
      double precision SFR_Gamma,SFR_Z, SFR_mT2
      integer SFR_Flag, SFR_Flag1, SFR_REV, SFR_REV2
      integer SFR_Nvec, SFR_VecID
      double precision SFR_Coeff, SFR_Vec
      integer SFR_JJ,SFR_KK,SFR_ID0
      double precision SFR_DHM, SFR_DHG
      double precision SFR_FB
      integer SFR_GamIN
      double precision SFR_GamC

      save /DataSFREP/
      

      integer iV,j

      call SFREP_SetUsed(iLine)

      if (SFR_NVec(iLine).eq.SFR_NvecMax) then
         call SFREP_ERROR(3,'SFREP_AddVec0: NVec exhausted!')
         return
      endif

      SFR_NVec(iLine)=SFR_NVec(iLine)+1
      iV = SFR_NVec(iLine)

      SFR_Coeff(iLine,iV) = c
      do j=1,4
         SFR_Vec(j,iLine,iV) = V(j)
      enddo
      SFR_VecID(iLine,iV) = ID
      end
      
c====================================================
c Sets all entries in line iLine1 as in iLine2
c
c (needed for the final rearrangement of hadrons)
c
c ***** used by PYSTRFT.FT.F *****
c
      subroutine SFREP_Copy(iLine1, iLine2)
      IMPLICIT NONE
      integer iLine1, iLine2

      integer SFR_NvecMax
      parameter (SFR_NvecMax=100) ! former value: 40

      common /DataSFREP/
     $     SFR_Gamma(4000),     ! reported Gamma
     $     SFR_Z(4000),         ! reported Z 
     $     SFR_mT2(4000),       ! reported mT2
     $     SFR_Flag(4000),      ! Flag: JS*(1..4)
     $     SFR_Flag1(4000),     ! Flag-Counter: added 1,if Flag is set
     $     SFR_REV(4000),       ! reverse ordering (unused)
     $     SFR_REV2(4000),      ! reverse ordering (unused)
     $     SFR_Nvec(4000),      ! number of vectors per particle
     $     SFR_Coeff(4000,SFR_NvecMax), ! coeff of vector
     $     SFR_Vec(4,4000,SFR_NvecMax), ! vector
     $     SFR_VecID(4000,SFR_NvecMax), ! ID of vector (cf. SFREPS)
     $     SFR_JJ(2,4000),      ! used for string-break ordering
     $     SFR_KK(2,4000),      !     -"-
     $     SFR_ID0(4000),       !     -"-
     $     SFR_DHM(4,4000),     ! reported coeffs for M-calc
     $     SFR_DHG(4,4000),     ! reported coeffs for Gamma-calc
     $     SFR_FB(4000),        ! reported parameter FB (???)
     $     SFR_GamIN(2,4000),   ! reported IN(1),IN(2)
     $     SFR_GamC(20,4000),   ! reported Coeffs for Gamma
     $     SFR_Ntot,            ! #(used entries) = #particles
     $     SFR_Nused(2)         ! indizes of lines which has been used
                                !                (min,max)

      integer SFR_Ntot
      integer SFR_Nused
      double precision SFR_Gamma,SFR_Z, SFR_mT2
      integer SFR_Flag, SFR_Flag1, SFR_REV, SFR_REV2
      integer SFR_Nvec, SFR_VecID
      double precision SFR_Coeff, SFR_Vec
      integer SFR_JJ,SFR_KK,SFR_ID0
      double precision SFR_DHM, SFR_DHG
      double precision SFR_FB
      integer SFR_GamIN
      double precision SFR_GamC

      save /DataSFREP/
      

      integer iV,nV,j

c      write(*,*) '==> SFREP_Copy(',iLine1,iLine2,')'

      call SFREP_SetUsed(iLine1)

      SFR_Gamma(iLine1) = SFR_Gamma(iLine2)
      SFR_Z(iLine1)     = SFR_Z(iLine2)
      SFR_mT2(iLine1)   = SFR_mT2(iLine2)
      
      SFR_Flag(iLine1)  = SFR_Flag(iLine2)
      SFR_Flag1(iLine1) = SFR_Flag1(iLine2)
      SFR_REV(iLine1)   = SFR_REV(iLine2) 
      SFR_REV2(iLine1)  = SFR_REV2(iLine2) 

      SFR_Nvec(iLine1)  = SFR_Nvec(iLine2)
      nV = SFR_NVec(iLine1)
      do iV=1,nV
         SFR_Coeff(iLine1,iV) = SFR_Coeff(iLine2,iV)
         do j=1,4
            SFR_Vec(j,iLine1,iV) = SFR_Vec(j,iLine2,iV)
         enddo
         SFR_VecID(iLine1,iV) = SFR_VecID(iLine2,iV)
      enddo

      do j=1,2
         SFR_JJ(j,iLine1) = SFR_JJ(j,iLine2)
         SFR_KK(j,iLine1) = SFR_KK(j,iLine2)
      enddo
      SFR_ID0(iLine1) = SFR_ID0(iLine2)
      do j=1,4
         SFR_DHM(j,iLine1) = SFR_DHM(j,iLine2)
         SFR_DHG(j,iLine1) = SFR_DHG(j,iLine2)
      enddo

      SFR_FB(iLine1) = SFR_FB(iLine2)

      do j=1,2
         SFR_GamIN(j,iLine1) = SFR_GamIN(j,iLine2)
      enddo
      do j=1,10
         SFR_GamC(j,iLine1) = SFR_GamC(j,iLine2)
      enddo

      
      call SFREP_Reset(iLine2,iLine2)

      end

c====================================================
c Print out the values from iLine1 to iLine2
c
c (just for internal purposes to have a print routine)
c Uses routines from "HandleSFREP_S.F" to find vectors in the global
c list. 
c
      subroutine SFREP_Write(iLine1,iLine2)
      IMPLICIT NONE
      integer iLine1,iLine2

      integer SFR_NvecMax
      parameter (SFR_NvecMax=100) ! former value: 40

      common /DataSFREP/
     $     SFR_Gamma(4000),     ! reported Gamma
     $     SFR_Z(4000),         ! reported Z 
     $     SFR_mT2(4000),       ! reported mT2
     $     SFR_Flag(4000),      ! Flag: JS*(1..4)
     $     SFR_Flag1(4000),     ! Flag-Counter: added 1,if Flag is set
     $     SFR_REV(4000),       ! reverse ordering (unused)
     $     SFR_REV2(4000),      ! reverse ordering (unused)
     $     SFR_Nvec(4000),      ! number of vectors per particle
     $     SFR_Coeff(4000,SFR_NvecMax), ! coeff of vector
     $     SFR_Vec(4,4000,SFR_NvecMax), ! vector
     $     SFR_VecID(4000,SFR_NvecMax), ! ID of vector (cf. SFREPS)
     $     SFR_JJ(2,4000),      ! used for string-break ordering
     $     SFR_KK(2,4000),      !     -"-
     $     SFR_ID0(4000),       !     -"-
     $     SFR_DHM(4,4000),     ! reported coeffs for M-calc
     $     SFR_DHG(4,4000),     ! reported coeffs for Gamma-calc
     $     SFR_FB(4000),        ! reported parameter FB (???)
     $     SFR_GamIN(2,4000),   ! reported IN(1),IN(2)
     $     SFR_GamC(20,4000),   ! reported Coeffs for Gamma
     $     SFR_Ntot,            ! #(used entries) = #particles
     $     SFR_Nused(2)         ! indizes of lines which has been used
                                !                (min,max)

      integer SFR_Ntot
      integer SFR_Nused
      double precision SFR_Gamma,SFR_Z, SFR_mT2
      integer SFR_Flag, SFR_Flag1, SFR_REV, SFR_REV2
      integer SFR_Nvec, SFR_VecID
      double precision SFR_Coeff, SFR_Vec
      integer SFR_JJ,SFR_KK,SFR_ID0
      double precision SFR_DHM, SFR_DHG
      double precision SFR_FB
      integer SFR_GamIN
      double precision SFR_GamC

      save /DataSFREP/
      

      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      integer N,NPAD,K
      double precision P,V
      SAVE /PYJETS/

      double precision xiP(1000),xiM(1000)

      integer iL,iL1,iL2,iV,i,j,ID,IDfirst!,h
      integer iString, iVec
      data iString /0/

      double precision Vsum(4)
      double precision VV(4), Vert1(4),Vert2(4), xh1,xh2,xh,hFak
      integer hID, hi

c$$$      integer ii1,ii2
c$$$      double precision FOUR,FOURVS
c$$$      FOUR(ii1,ii2)=-SFR_Vec(1,iL,ii1)*SFR_Vec(1,iL,ii2)
c$$$     $     -SFR_Vec(2,iL,ii1)*SFR_Vec(2,iL,ii2)
c$$$     $     -SFR_Vec(3,iL,ii1)*SFR_Vec(3,iL,ii2)
c$$$     $     +SFR_Vec(4,iL,ii1)*SFR_Vec(4,iL,ii2)
c$$$      FOURVS(ii1)=-SFR_Vec(1,iL,ii1)*Vsum(1)
c$$$     $     -SFR_Vec(2,iL,ii1)*Vsum(2)
c$$$     $     -SFR_Vec(3,iL,ii1)*Vsum(3)
c$$$     $     +SFR_Vec(4,iL,ii1)*Vsum(4)

      if (iLine2.lt.iLine1) then
         iL1 = 1
         iL2 = SFR_Ntot
      else
         iL1 = iLine1
         iL2 = iLine2
      endif

      do i=1,1000
         xiP(i) = 1d0
         xiM(i) = 0d0
      enddo

      do iL=iL1,iL2
         if (SFR_Flag(iL).eq.0) then 
            write(*,1000) iL
            goto 900
         endif

         do j=1,4
            Vsum(j) = 0d0
            Vert1(j) = 0d0
            Vert2(j) = 0d0
         enddo

         write(*,1001) iL,'Gamma^2,z,mT^2,pT^2, m',
     $        SFR_Gamma(iL),SFR_Z(iL),SFR_mT2(iL),
     $        SFR_mT2(iL)-P(iL,5)**2,P(iL,5)
         write(*,1002) iL,'Flag ',
     $        SFR_Flag(iL),SFR_Flag1(iL),SFR_REV(iL),SFR_REV2(iL)
         write(*,1003) iL,SFR_JJ(1,iL),SFR_JJ(2,iL),
     $        SFR_KK(1,iL),SFR_KK(2,iL),
     $        (SFR_JJ(1,iL)-SFR_ID0(iL))/4+1,
     $        (SFR_KK(1,iL)-1-SFR_ID0(iL))/4+1,
     $        (SFR_JJ(2,iL)-SFR_ID0(iL))/4+1,
     $        (SFR_KK(2,iL)-1-SFR_ID0(iL))/4+1

c         if (abs(SFR_Flag(iL)).ge.3) then
            write(*,1001) iL,'** DHM',(SFR_DHM(j,iL),j=1,4)
            write(*,1001) iL,'** DHG',(SFR_DHG(j,iL),j=1,4)
            write(*,1001) iL,'** FB ',SFR_FB(iL)
c         endif

         write(*,1050) iL,'** GammaC', SFR_GamIN(1,iL),SFR_GamIN(2,iL),
     $        (SFR_GamIN(1,iL)+j*4,SFR_GamC(j*2+1,iL),
     $        SFR_GamIN(1,iL)+j*4+1,SFR_GamC(j*2+2,iL),
     $        j=0,(SFR_GamIN(2,iL)-1-SFR_GamIN(1,iL))/4)
            

c$$$         do iV=1,SFR_NVec(iL)
c$$$            do h=1,iV
c$$$               write(*,*) 'Four:',iV,h,':',
c$$$     $              SFR_VecID(iL,IV),SFR_VecID(iL,h),
c$$$     $              FOUR(iV,h)
c$$$            enddo
c$$$         enddo

         do iV=1,SFR_NVec(iL)

            do j=1,4
               Vsum(j) = Vsum(j)+SFR_Coeff(iL,iV)*SFR_Vec(j,iL,iV)
            enddo

            call SFREPS_FindVecI(2, iL,iV, iString,iVec)

            if (iString.ne.0) then
               write(*,1011) iL,SFR_Coeff(iL,iV),
     $              (SFR_Vec(j,iL,iV),j=1,4),SFR_VecID(iL,iV),
     $              2,iString,iVec
            else
               write(*,1010) iL,SFR_Coeff(iL,iV),
     $              (SFR_Vec(j,iL,iV),j=1,4),SFR_VecID(iL,iV)
            endif

            if (iString.ne.0) then

               ID = SFR_VecID(iL,iV)
               call SFREPS_GetVec(2,iString,1, VV,IDfirst)

c               write(*,*) 'iV,ID,+/-,xiP,xiM:',
c     $              iV,ID,mod((ID-IDfirst),2),xiP(ID),xiM(ID)

               if (mod((ID-IDfirst),2).eq.0) then 
c... (+)-Vector:
                  do j=1,4
                     Vert1(j) = Vert1(j) + xiP(ID)*SFR_Vec(j,iL,iV)
                  enddo
                  xiP(ID) = xiP(ID)-SFR_Coeff(iL,iV)
                  do j=1,4
                     Vert2(j) = Vert2(j) + xiP(ID)*SFR_Vec(j,iL,iV)
                  enddo

               else             
c... (-)-Vector: 
                  do j=1,4
                     Vert1(j) = Vert1(j) + xiM(ID)*SFR_Vec(j,iL,iV)
                  enddo
                  xiM(ID) = xiM(ID)+SFR_Coeff(iL,iV)
                  do j=1,4
                     Vert2(j) = Vert2(j) + xiM(ID)*SFR_Vec(j,iL,iV)
                  enddo
               endif
            endif

         enddo

         if (SFR_NVec(iL).gt.0) then

            call SFREPS_GetVecID(3,iString,iVec, VV,iL)

c            call SFREPS_FindVec(3,Vsum, iString,iVec)

            if (iVec.ne.0) then
               write(*,1021) iL,(VV(j),j=1,4),3,iString,iVec
            endif
            write(*,1022) iL,(Vsum(j),j=1,4),'Sum'
         endif

         xh1 = Vert1(4)**2-Vert1(1)**2-Vert1(2)**2-Vert1(3)**2
         xh2 = sqrt(max(xh1,0d0))
         write(*,1040) iL,'Vertex 1',(Vert1(j),j=1,4),
     $        xh1,xh2,xh1-SFR_Gamma(iL)
         xh1 = Vert2(4)**2-Vert2(1)**2-Vert2(2)**2-Vert2(3)**2
         xh2 = sqrt(max(xh1,0d0))
         write(*,1040) iL,'Vertex 2',(Vert2(j),j=1,4),
     $        xh1,xh2,xh1-SFR_Gamma(iL)
         write(*,*)'    '

c... write checks of Gamma, z, mT^2:

         if (SFR_NVec(iL).ge.6) then

c... check Gamma-calculation:
            
            xh1 = 0d0
            xh2 = 0d0
            do iV=5,SFR_NVec(iL)
               if (SFR_VecID(iL,iV).eq.SFR_GamIN(1,iL))
     $              xh1=SFR_Coeff(iL,iV)
               if (SFR_VecID(iL,iV).eq.SFR_GamIN(2,iL))
     $              xh2=SFR_Coeff(iL,iV)
            enddo

            xh = SFR_DHG(1,iL)
     $           + SFR_DHG(2,iL)*xh1
     $           + SFR_DHG(3,iL)*xh2
     $           + SFR_DHG(4,iL)*xh1*xh2
            
            write(*,1001) iL,' recalc Gamma',xh,SFR_Gamma(iL),
     $           xh-SFR_Gamma(iL)
            
c... calculate Vertex-Vector:
c... assume to be set: xh1,xh2, iString

            if (SFR_GamIN(1,iL).eq.0) goto 100 ! skip reconstruct

            write(*,1001) iL,' reconstruct V'

            call SetArrayD(4,Vsum,0d0)

            do hi=0,(SFR_GamIN(2,iL)-1-SFR_GamIN(1,iL))/4
               hID = SFR_GamIN(1,iL)+hi*4
               call SFREPS_GetVecID(2,iString,iVec, VV, hID)
               if (hID.eq.SFR_GamIN(1,iL)) then
c                  hFak=SFR_GamC(2*hi+1,iL)-xh1
c                  hFak=SFR_GamC(2*hi+1,iL)-sign(xh1,dble(SFR_Flag(iL)))
                  hFak=SFR_GamC(2*hi+1,iL)-sign(1,SFR_Flag(iL))*xh1
               else
                  hFak=SFR_GamC(2*hi+1,iL)/2
c                 hFak=SFR_GamC(2*hi+1,iL)*
c    $                 3.71515E+00/(3.71515E+00+3.15079E+00)
               endif
               call AddVert(Vsum, hFak, VV)
               write(*,1010) iL,hFak,(VV(j),j=1,4),hID

               hID=hID+1
               call SFREPS_GetVecID(2,iString,iVec, VV, hID)
               if (hID.eq.SFR_GamIN(2,iL)) then
c                  hFak=SFR_GamC(2*hi+2,iL)+xh2
c                  hFak=SFR_GamC(2*hi+2,iL)+sign(xh2,dble(SFR_Flag(iL)))
                  hFak=SFR_GamC(2*hi+2,iL)+sign(1,SFR_Flag(iL))*xh2
               else
                  hFak=SFR_GamC(2*hi+2,iL)/2
c                 hFak=SFR_GamC(2*hi+2,iL)*
c    $              3.15079E+00/(3.71515E+00+3.15079E+00)   
               endif
               call AddVert(Vsum, hFak, VV)
               write(*,1010) iL,hFak,(VV(j),j=1,4),hID
               
            enddo
            write(*,1022) iL,(Vsum(j),j=1,4),'Sum'
            call SetArrayArr(4,VV,Vsum)

            if (SFR_FLag(iL).gt.0) then 
               call AddVert(VV,-1d0,Vert2)
               write(*,1022) iL,(VV(j),j=1,4),'Sum-Vert2'
            else
               call AddVert(VV,-1d0,Vert1)
               write(*,1022) iL,(VV(j),j=1,4),'Sum-Vert1'
            endif

            xh1= Vsum(4)**2-Vsum(1)**2-Vsum(2)**2-Vsum(3)**2
            xh2= sqrt(max(xh1,0d0))

            write(*,1040) iL,'recon. V',(Vsum(j),j=1,4),
     $           xh1,xh2,xh1-SFR_Gamma(iL)

 100       continue

c$$$            value=Vsum(1)**2+Vsum(2)**2+VSum(3)**2-Vsum(4)**2
c$$$            write(*,1005) iL,SFR_NVec(iL)-1,
c$$$     $           0.5d0*value/FOURVS(SFR_NVec(iL)-1)
c$$$            write(*,1005) iL,SFR_NVec(iL),
c$$$     $           0.5D0*value/FOURVS(SFR_NVec(iL))

c$$$            value = (SFR_Coeff(iL,1)+SFR_Coeff(iL,2))**2
c$$$     $           +(SFR_Coeff(iL,3)+SFR_Coeff(iL,4))**2
c$$$     $           +P(iL,5)**2
c$$$            write(*,1001) iL,': check: pT^2 (A): ',
c$$$     $           value,value-SFR_mT2(iL)

c$$$            value = SFR_Coeff(iL,5)*SFR_Vec(4,iL,5)
c$$$     $           *SFR_Coeff(iL,6)*SFR_Vec(4,iL,6) * 4
c$$$c     $           +P(iL,5)**2
c$$$            write(*,1001) iL,': check: pT^2 (B): ',
c$$$     $           value, value-SFR_mT2(iL)       

c$$$            do j=1,4
c$$$               Vsum(j) = 0d0
c$$$            enddo
c$$$            do iV=1,SFR_NVec(iL)-2
c$$$               do j=1,4
c$$$                  Vsum(j) = Vsum(j) * SFR_Coeff(iL,iV)*SFR_Vec(j,iL,iV)
c$$$               enddo
c$$$            enddo
c$$$            value = Vsum(1)**2+Vsum(2)**2+Vsum(3)**2-Vsum(4)**2
c$$$
c$$$            write(*,*) ': check: 0    (C): ',
c$$$     $           value, 0d0

c            do i=1,4
c               write(*,1030), i,(FOUR(i,j),j=1,4)
c            enddo

c$$$            value = SFR_Coeff(iL,1)**2+SFR_Coeff(iL,2)**2
c$$$     $           +SFR_Coeff(iL,3)**2+SFR_Coeff(iL,4)**2
c$$$     $           +2*(SFR_Coeff(iL,1)*SFR_Coeff(iL,2)*FOUR(1,2)
c$$$     $           +SFR_Coeff(iL,1)*SFR_Coeff(iL,4)*FOUR(1,4)
c$$$     $           +SFR_Coeff(iL,2)*SFR_Coeff(iL,3)*FOUR(2,3)
c$$$     $           +SFR_Coeff(iL,3)*SFR_Coeff(iL,4)*FOUR(4,5) )
c$$$     $           +P(iL,5)**2
c$$$            write(*,1001) iL,': check: pT^2 (D): ',
c$$$     $           value,value-SFR_mT2(iL)
c$$$
c$$$            value =SFR_DHM(1,iL)
c$$$     $           +SFR_DHM(2,iL)*SFR_Coeff(iL,SFR_NVec(iL)-1)
c$$$     $           +SFR_DHM(3,iL)*SFR_Coeff(iL,SFR_NVec(iL))
c$$$     $           +SFR_DHM(4,iL)*SFR_Coeff(iL,SFR_NVec(iL)-1)
c$$$     $           *SFR_Coeff(iL,SFR_NVec(iL))
c$$$            write(*,1001) iL,': check: m^2  (E): ',
c$$$     $           value,value-P(iL,5)**2
c$$$
c$$$            value =SFR_DHG(1,iL)
c$$$     $           +SFR_DHG(2,iL)*SFR_Coeff(iL,SFR_NVec(iL)-1)
c$$$     $           +SFR_DHG(3,iL)*SFR_Coeff(iL,SFR_NVec(iL))
c$$$     $           +SFR_DHG(4,iL)*SFR_Coeff(iL,SFR_NVec(iL)-1)
c$$$     $           *SFR_Coeff(iL,SFR_NVec(iL))
c$$$            write(*,1001) iL,': check: G    (F): ',
c$$$     $           value,value-SFR_Gamma(iL)

         endif

 900     write(*,*)

      enddo

 1000 FORMAT(0P,1i4,': ==========')
 1001 FORMAT(0P,1i4,': ',A,'=',1P,6e13.5)
 1002 FORMAT(0P,1i4,': ',A,'=',4i7)
 1003 FORMAT(0P,1i4,': StringBreak: ',4i4,' -> (',2i3,')(',2i3,')')
 1005 FORMAT(0P,1i4,': FOUR(V,',i1,') = ',1P,e12.5)
 1006 FORMAT(0P,1i4,': ',A)

 1010 FORMAT(0P,1i4,':    ',1P,e12.5,' x(',3(e12.5,','),e12.5,') :['
     $     ,i3,']')
 1011 FORMAT(0P,1i4,':    ',1P,e12.5,' x(',3(e12.5,','),e12.5,') :['
     $     ,i3,'] =<<SFREP_S,',i1,'>>[',i2,',',i2,']')
 1020 FORMAT(0P,1i4,':     ==           (',1P,3(e12.5,','),e12.5,')')
 1021 FORMAT(0P,1i4,':     ==           (',1P,3(e12.5,','),e12.5,') =',
     $     '<<SFREP_S,',i1,'>>[',i2,',',i2,']')
 1022 FORMAT(0P,1i4,':     ==           (',1P,3(e12.5,','),e12.5,') = ',
     $     A)   
 1030 FORMAT(0P,1i4,': ',1P,4e13.5)
 1040 FORMAT(0P,1i4,': ',A,':',1P,4e13.5,' ; ',3e13.5)

 1050 FORMAT(0P,1i4,': ',A,':  IN(1)=',i3,' IN(2)=',i3,' ',1P,
     $     10(' C(',i3,')=',e12.5))
      end

c====================================================
c Convert the String Breaks into Gluons
c
      subroutine SFREP_ConvertToGluons(iLine1,iLine2,iFile)
      IMPLICIT NONE
      integer iLine1,iLine2,iFile

      integer SFR_NvecMax
      parameter (SFR_NvecMax=100) ! former value: 40

      common /DataSFREP/
     $     SFR_Gamma(4000),     ! reported Gamma
     $     SFR_Z(4000),         ! reported Z 
     $     SFR_mT2(4000),       ! reported mT2
     $     SFR_Flag(4000),      ! Flag: JS*(1..4)
     $     SFR_Flag1(4000),     ! Flag-Counter: added 1,if Flag is set
     $     SFR_REV(4000),       ! reverse ordering (unused)
     $     SFR_REV2(4000),      ! reverse ordering (unused)
     $     SFR_Nvec(4000),      ! number of vectors per particle
     $     SFR_Coeff(4000,SFR_NvecMax), ! coeff of vector
     $     SFR_Vec(4,4000,SFR_NvecMax), ! vector
     $     SFR_VecID(4000,SFR_NvecMax), ! ID of vector (cf. SFREPS)
     $     SFR_JJ(2,4000),      ! used for string-break ordering
     $     SFR_KK(2,4000),      !     -"-
     $     SFR_ID0(4000),       !     -"-
     $     SFR_DHM(4,4000),     ! reported coeffs for M-calc
     $     SFR_DHG(4,4000),     ! reported coeffs for Gamma-calc
     $     SFR_FB(4000),        ! reported parameter FB (???)
     $     SFR_GamIN(2,4000),   ! reported IN(1),IN(2)
     $     SFR_GamC(20,4000),   ! reported Coeffs for Gamma
     $     SFR_Ntot,            ! #(used entries) = #particles
     $     SFR_Nused(2)         ! indizes of lines which has been used
                                !                (min,max)

      integer SFR_Ntot
      integer SFR_Nused
      double precision SFR_Gamma,SFR_Z, SFR_mT2
      integer SFR_Flag, SFR_Flag1, SFR_REV, SFR_REV2
      integer SFR_Nvec, SFR_VecID
      double precision SFR_Coeff, SFR_Vec
      integer SFR_JJ,SFR_KK,SFR_ID0
      double precision SFR_DHM, SFR_DHG
      double precision SFR_FB
      integer SFR_GamIN
      double precision SFR_GamC

      save /DataSFREP/

      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      integer N,NPAD,K
      double precision P,V
      SAVE /PYJETS/

      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      integer MSTU,MSTJ
      double precision PARU,PARJ
      SAVE /PYDAT1/

      integer iL,iL1,iL2,iV,i,j,ID,IDfirst,iiL!,h
      integer iString, iVec
      data iString /0/

      double precision VV(4),Vsum(4),meff

      double precision RRR(1000)
      logical LLL(1000)
      integer iLnew,iLnew1
      logical FirstParticle
      integer RevAdd
      double precision DPS(5)
      double precision HHBZ,HHPMT,HHPEZ


      
      if (iLine2.lt.iLine1) then
         iL1 = 1
         iL2 = SFR_Ntot
      else
         iL1 = iLine1
         iL2 = iLine2
      endif

      call SFREP_Reset(iL2+1,4000)
      iLnew = SFR_Ntot

      FirstParticle = .true.
      RevAdd = 1

      do iL=iL1,iL2

         if (SFR_Flag(iL).eq.0) then 
            if(iFile.gt.0) write(iFile,*) 'Sums of this string~~~~~~~~'
            
            do iiL=1,1000
               if (LLL(iiL)) then
                   if(iFile.gt.0) then
             write(iFile,'(A,i3,A,1P,6e13.5)') '.....',iiL,'=',RRR(iiL)
                   endif
               endif

               LLL(iiL)=.false.
               RRR(iiL)=0d0
            end do
            if(iFile.gt.0) write(iFile,1000) iL

            FirstParticle = .true.

            goto 900
         endif

         iLnew = iLnew+1
         SFR_ID0(iLnew)=SFR_ID0(iL)
         if (FirstParticle) then
            iLnew = iLnew+1
            SFR_ID0(iLnew)=SFR_ID0(iL)
         endif
         

         if(iFile.gt.0) write(iFile,1003)
     $        iL,SFR_JJ(1,iL),SFR_JJ(2,iL),
     $        SFR_KK(1,iL),SFR_KK(2,iL),
     $        (SFR_JJ(1,iL)-SFR_ID0(iL))/4+1,
     $        (SFR_KK(1,iL)-1-SFR_ID0(iL))/4+1,
     $        (SFR_JJ(2,iL)-SFR_ID0(iL))/4+1,
     $        (SFR_KK(2,iL)-1-SFR_ID0(iL))/4+1

         do j=1,4
            Vsum(j) = 0d0
         enddo

         iLnew1 = iLnew
         if (SFR_REV(iL).gt.0) then
            iLnew1=iLnew1+RevAdd
            RevAdd = -RevAdd
         else
            RevAdd = 1
         endif

         do iV=1,SFR_NVec(iL)

            do j=1,4
               Vsum(j) = Vsum(j)+SFR_Coeff(iL,iV)*SFR_Vec(j,iL,iV)
            enddo

            call SFREPS_FindVecI(2, iL,iV, iString,iVec)
            if (iString.ne.0) then
               if(iFile.gt.0) write(iFile,1011)
     $              iL,SFR_Coeff(iL,iV),
     $              (SFR_Vec(j,iL,iV),j=1,4),SFR_VecID(iL,iV),
     $              2,iString,iVec
             LLL(SFR_VecID(iL,iV)) = .TRUE.
             RRR(SFR_VecID(iL,iV)) =
     $           RRR(SFR_VecID(iL,iV))+SFR_Coeff(iL,iV)  

             ID = SFR_VecID(iL,iV)
             call SFREPS_GetVec(2,iString,1, VV,IDfirst)
             do j=1,4
                VV(j)=SFR_Vec(j,iL,iV)
             enddo

             if (mod((ID-IDfirst),2).eq.0) then 
c...  (+)-Vector:
c                write(*,*) 'ID,+',ID
                call SFREP_AddVec0(iLnew1-1,SFR_Coeff(iL,iV),VV,
     $               SFR_VecID(iL,iV))
                SFR_KK(2,iLnew1-1)=SFR_KK(1,iL)
                SFR_JJ(2,iLnew1-1)=SFR_JJ(1,iL)
                SFR_Flag1(iLnew1-1) = iL ! abuse of SFR_Flag1

               else             
c... (-)-Vector: 
c                  write(*,*) 'ID,-',ID
                  call SFREP_AddVec0(iLnew1,SFR_Coeff(iL,iV),VV,
     $               SFR_VecID(iL,iV))
                  SFR_KK(1,iLnew1)=SFR_KK(2,iL)
                  SFR_JJ(1,iLnew1)=SFR_JJ(2,iL)
                  SFR_Flag1(iLnew1) = iL ! abuse of SFR_Flag1
               endif

            else
c               write(*,1010) iL,SFR_Coeff(iL,iV),
c     $              (SFR_Vec(j,iL,iV),j=1,4),SFR_VecID(iL,iV)
            endif

c            LLL(SFR_VecID(iL,iV)) = .TRUE.
c            RRR(SFR_VecID(iL,iV)) =
c     $           RRR(SFR_VecID(iL,iV))+SFR_Coeff(iL,iV)

         enddo


         FirstParticle = .false.

 900     if(iFile.gt.0) write(iFile,*)
      enddo


      if(iFile.gt.0) write(iFile,*) 'Sums of this string~~~~~~~~'
      do iiL=1,1000
         if (LLL(iiL)) then
            if(iFile.gt.0) then
          write(iFile,'(A,i3,A,1P,6e13.5)') '.....',iiL,'=',RRR(iiL)
            endif
         endif

         LLL(iiL)=.false.
         RRR(iiL)=0d0
      end do
      if(iFile.gt.0) write(iFile,1000) iL

cccccccccccccccccccccc

      do iL=1,iL2
         if (K(iL,1).lt.10) K(iL,1)=K(iL,1)+10
      enddo

cccccccccccccccccccccc

      do iL=iL2+1,iLnew
         if(iFile.gt.0) write(iFile,1003)
     $        iL,SFR_JJ(1,iL),SFR_JJ(2,iL),
     $        SFR_KK(1,iL),SFR_KK(2,iL),
     $        (SFR_JJ(1,iL)-SFR_ID0(iL))/4+1,
     $        (SFR_KK(1,iL)-1-SFR_ID0(iL))/4+1,
     $        (SFR_JJ(2,iL)-SFR_ID0(iL))/4+1,
     $        (SFR_KK(2,iL)-1-SFR_ID0(iL))/4+1
         do j=1,4
            Vsum(j) = 0d0
         enddo
         do iV=1,SFR_NVec(iL)

            do j=1,4
               Vsum(j) = Vsum(j)+SFR_Coeff(iL,iV)*SFR_Vec(j,iL,iV)
            enddo

            call SFREPS_FindVecI(2, iL,iV, iString,iVec)
            if (iString.ne.0) then
               if(iFile.gt.0) write(iFile,1011)
     $              iL,SFR_Coeff(iL,iV),
     $              (SFR_Vec(j,iL,iV),j=1,4),SFR_VecID(iL,iV),
     $              2,iString,iVec
            endif
         enddo
         meff=Vsum(4)**2-Vsum(1)**2-Vsum(2)**2-Vsum(3)**2
         if (meff.gt.0d0) then
            meff=sqrt(meff)
         else
            meff=-sqrt(-meff)
         endif

         if(iFile.gt.0) write(iFile,1023) iL,(Vsum(j),j=1,4),meff
         if(iFile.gt.0) write(iFile,*)

         K(iL,1) = 2
         if (SFR_JJ(1,iL).eq.0) then ! first particle
            K(iL,2) = K(K(SFR_Flag1(iL),3),2)
         else if (SFR_JJ(2,iL).eq.0) then ! last particle
            K(iL,2) = K(K(SFR_Flag1(iL),3),2)
            K(iL,1) = 1
         else
            K(iL,2) = 21
         endif
         K(iL,3) = SFR_Flag1(iL)
         K(iL,4) = 0
         K(iL,5) = 0

         do j=1,4
            P(iL,j) = Vsum(j)
         enddo
         P(iL,5) = meff

      enddo

      N = iLnew
      MSTU(70)=MSTU(70)+1
      MSTU(70+MSTU(70))=iLnew

c      call PYLIST(2)

      iString = 1
      iL1 = iL2+1
      do iL=iL1,iLnew
         if (K(iL,1).eq.1) then
            call SFREPS_GetBoost(iString,DPS)
            if (ABS(DPS(3)).LT.0.99D0*DPS(4)) then
               call PYROBO(iL1,iL,0D0,0D0,
     $              DPS(1)/DPS(4),DPS(2)/DPS(4),DPS(3)/DPS(4))
            else
         HHBZ=SQRT(MAX(1D-6,DPS(4)+DPS(3))/MAX(1D-6,DPS(4)-DPS(3)))
               do iiL=iL1,iL
                  HHPMT=P(IIL,1)**2+P(IIL,2)**2+P(IIL,5)**2
                  IF(P(IIL,3).GT.0D0) THEN
                     HHPEZ=(P(IIL,4)+P(IIL,3))*HHBZ
                     P(IIL,3)=0.5D0*(HHPEZ-HHPMT/HHPEZ)
                     P(IIL,4)=0.5D0*(HHPEZ+HHPMT/HHPEZ)
                  ELSE
                     HHPEZ=(P(IIL,4)-P(IIL,3))/HHBZ
                     P(IIL,3)=-0.5D0*(HHPEZ-HHPMT/HHPEZ)
                     P(IIL,4)=0.5D0*(HHPEZ+HHPMT/HHPEZ)
                  ENDIF
               enddo

            endif
            do j=1,4
               Vsum(j) = 0d0
            enddo
            do iiL=iL1,iL
               do j=1,4
                  Vsum(j) = Vsum(j) + P(iiL,j)
               enddo
            enddo
            if (iFile.gt.0) write(iFile,1030) iString,(Vsum(j),j=1,4)

            iString=iString+1
            iL1=iL+1
         endif
      enddo

c      call PYLIST(2)


 1000 FORMAT(0P,1i4,': ==========')
 1001 FORMAT(0P,1i4,': ',A,'=',1P,6e13.5)
 1002 FORMAT(0P,1i4,': ',A,'=',4i7)
 1003 FORMAT(0P,1i4,': StringBreak: ',4i4,' -> (',2i3,')(',2i3,')')
 1005 FORMAT(0P,1i4,': FOUR(V,',i1,') = ',1P,e12.5)
 1006 FORMAT(0P,1i4,': ',A)

 1010 FORMAT(0P,1i4,':    ',1P,e12.5,' x(',3(e12.5,','),e12.5,') :['
     $     ,i3,']')
 1011 FORMAT(0P,1i4,':    ',1P,e12.5,' x(',3(e12.5,','),e12.5,') :['
     $     ,i3,'] =<<SFREP_S,',i1,'>>[',i2,',',i2,']')
 1020 FORMAT(0P,1i4,':     ==           (',1P,3(e12.5,','),e12.5,')')
 1021 FORMAT(0P,1i4,':     ==           (',1P,3(e12.5,','),e12.5,') =',
     $     '<<SFREP_S,',i1,'>>[',i2,',',i2,']')
 1022 FORMAT(0P,1i4,':     ==           (',1P,3(e12.5,','),e12.5,') = ',
     $     A) 
 1023 FORMAT(0P,1i4,':     ==           (',1P,3(e12.5,','),e12.5,') = ',
     $     e12.5) 
 1030 FORMAT(0P,1i4,': ',1P,4e13.5)
 1040 FORMAT(0P,1i4,': ',A,':',1P,4e13.5,' ; ',3e13.5)

 1050 FORMAT(0P,1i4,': ',A,':  IN(1)=',i3,' IN(2)=',i3,' ',1P,
     $     10(' C(',i3,')=',e12.5))
      end

c====================================================
c... Reconstruct "last hadron" for a (1+1) String
c
c Uses the reported values z_(h,i)^+ and z_(h,i)^- to set up the vectors
c for the final hadron (flag=4).
c
c Please note: "transverse momenta"
c A "normal" particle with SFR_Flag(iL)>0 would have transverse momenta
c according the first column, but a particle with SFR_Flag(iLX)==4 lies
c at the border, therefore by trivial string breaking order, colums 2 is
c valid.
c   (iL,1)  =(iL-1,2)  =(iL-1,2)
c   (iL,2)  =(iL+1,1)  =(iL+1,2)
c   (iL,3)  =(iL-1,4)  =(iL-1,4)
c   (iL,4)  =(iL+1,3)  =(iL+1,4)
c

      subroutine SFREP_ReconLastQQ(iLine1,iLine2)
      IMPLICIT NONE
      integer iLine1,iLine2

      integer SFR_NvecMax
      parameter (SFR_NvecMax=100) ! former value: 40

      common /DataSFREP/
     $     SFR_Gamma(4000),     ! reported Gamma
     $     SFR_Z(4000),         ! reported Z 
     $     SFR_mT2(4000),       ! reported mT2
     $     SFR_Flag(4000),      ! Flag: JS*(1..4)
     $     SFR_Flag1(4000),     ! Flag-Counter: added 1,if Flag is set
     $     SFR_REV(4000),       ! reverse ordering (unused)
     $     SFR_REV2(4000),      ! reverse ordering (unused)
     $     SFR_Nvec(4000),      ! number of vectors per particle
     $     SFR_Coeff(4000,SFR_NvecMax), ! coeff of vector
     $     SFR_Vec(4,4000,SFR_NvecMax), ! vector
     $     SFR_VecID(4000,SFR_NvecMax), ! ID of vector (cf. SFREPS)
     $     SFR_JJ(2,4000),      ! used for string-break ordering
     $     SFR_KK(2,4000),      !     -"-
     $     SFR_ID0(4000),       !     -"-
     $     SFR_DHM(4,4000),     ! reported coeffs for M-calc
     $     SFR_DHG(4,4000),     ! reported coeffs for Gamma-calc
     $     SFR_FB(4000),        ! reported parameter FB (???)
     $     SFR_GamIN(2,4000),   ! reported IN(1),IN(2)
     $     SFR_GamC(20,4000),   ! reported Coeffs for Gamma
     $     SFR_Ntot,            ! #(used entries) = #particles
     $     SFR_Nused(2)         ! indizes of lines which has been used
                                !                (min,max)

      integer SFR_Ntot
      integer SFR_Nused
      double precision SFR_Gamma,SFR_Z, SFR_mT2
      integer SFR_Flag, SFR_Flag1, SFR_REV, SFR_REV2
      integer SFR_Nvec, SFR_VecID
      double precision SFR_Coeff, SFR_Vec
      integer SFR_JJ,SFR_KK,SFR_ID0
      double precision SFR_DHM, SFR_DHG
      double precision SFR_FB
      integer SFR_GamIN
      double precision SFR_GamC

      save /DataSFREP/
      

      integer iL,iL4,iV,j,iL0
      double precision c(6)
      integer iString, iVec
      integer ID,ID0
      double precision VV(4)
      integer JJ(2),KK(2)
      logical ccOK
c      double precision Vsum(4),h
      

      do j=1,6
         c(j) = 0d0
      enddo

      do iL=iLine1,iLine2
         if (abs(SFR_Flag(iL)).eq.4) then
            iL4 = iL
         else
            do iV=5,6
               c(iV) = c(iV) + SFR_Coeff(iL,iV)
            enddo
         endif
      enddo

c      if (iL4.le.iLine1.or.iL4.ge.iLine2) then
c         write(*,*) 'SFREP_ReconLast1: Problem iL4,iL1,iL2=',
c     $        iL4,iLine1,iLine2
c      endif
      
      iL0 = iLine1
      if (iL4.eq.iL0) iL0 = iLine2

      SFR_Nvec(iL4)  = SFR_Nvec(iL0)

c... Transverse momenta:

      if (iL4.eq.iLine1) then
         if (SFR_FLag(iL4).eq.4) then
            call SFREP_ERROR(4,'SFREP_ReconLastQQ: iL4=iL1,Flag==4.')
            return
         else
            call SFREP_CopyVec(iL4,1, iL4+1,2, -1d0)
            call SFREP_CopyVec(iL4,2, iL4+1,2,  0d0)
            call SFREP_CopyVec(iL4,3, iL4+1,4, -1d0)
            call SFREP_CopyVec(iL4,4, iL4+1,4,  0d0)
         endif

      elseif (iL4.eq.iLine2) then
         if (SFR_FLag(iL4).eq.-4) then
            call SFREP_ERROR(5,'SFREP_ReconLastQQ: iL4=iL2,Flag==-4.')
            return
         else
            call SFREP_CopyVec(iL4,1, iL4-1,2, -1d0)
            call SFREP_CopyVec(iL4,2, iL4-1,2,  0d0)
            call SFREP_CopyVec(iL4,3, iL4-1,4, -1d0)
            call SFREP_CopyVec(iL4,4, iL4-1,4,  0d0)
         endif

      else
         if (SFR_FLag(iL4).eq.4) then
            call SFREP_CopyVec(iL4,1, iL4-1,2, -1d0)
            call SFREP_CopyVec(iL4,2, iL4+1,2, -1d0)
            call SFREP_CopyVec(iL4,3, iL4-1,4, -1d0)
            call SFREP_CopyVec(iL4,4, iL4+1,4, -1d0)
         else
            call SFREP_CopyVec(iL4,1, iL4+1,2, -1d0)
            call SFREP_CopyVec(iL4,2, iL4-1,2, -1d0)
            call SFREP_CopyVec(iL4,3, iL4+1,4, -1d0)
            call SFREP_CopyVec(iL4,4, iL4-1,4, -1d0)
         endif
      endif

c... Longitudinal momenta:

      do iV=5,6
         SFR_Coeff(iL4,iV) = (1-c(iV))
         do j=1,4
            SFR_Vec(j,iL4,iV) = SFR_Vec(j,iL0,iV)
         enddo 
         SFR_VecID(iL4,iV) = SFR_VecID(iL0,iV)
      enddo

c... Set StringBreak:

      call SFREPS_FindVecI(2, iLine1,5, iString,iVec) ! to set iString
      call SFREPS_GetVec(2,iString,1, VV,ID0) ! to set ID0

      do iL=iLine1,iLine2
         JJ(1)=20000
         JJ(2)=-1
         KK(1)=20000
         KK(2)=-1

         do iV=5,SFR_NVec(iL)
            ID = SFR_VecID(iL,iV)
            if (mod((ID-ID0),2).eq.1) then
               if (ID.lt.KK(1)) KK(1)=SFR_VecID(iL,iV)
               if (ID.gt.KK(2)) KK(2)=SFR_VecID(iL,iV)
            else
               if (ID.lt.JJ(1)) JJ(1)=SFR_VecID(iL,iV)
               if (ID.gt.JJ(2)) JJ(2)=SFR_VecID(iL,iV)
            endif
         enddo
         SFR_ID0(iL) = ID0
         do j=1,2
            SFR_KK(j,iL) = KK(j)
            SFR_JJ(j,iL) = JJ(j)
         enddo
      enddo

c      if (iL4.le.iLine1.or.iL4.ge.iLine2) then
c         write(*,*) 'SFREP_ReconLast1: Solved, continue...'
cc         call PYLIST(2)
cc         call SFREP_Write(iLine1,iLine2)
c      endif

c... Check Correction:

      call SFREP_CheckCorr('SFREP_ReconLastQQ',
     $     iLine1,iLine2,iL4,iString,ccOK)

      end

c====================================================
c... Reconstruct "last hadron" for a String with Gluons
c
c Uses the reported values z_(h,i)^+ and z_(h,i)^- to set up the vectors
c for the final hadron (flag=4).
c
c ... Implemntation not fully checked !!! ...
c
      subroutine SFREP_ReconLastQgQ(iLine1,iLine2)
      IMPLICIT NONE
      integer iLine1,iLine2

      integer SFR_NvecMax
      parameter (SFR_NvecMax=100) ! former value: 40

      common /DataSFREP/
     $     SFR_Gamma(4000),     ! reported Gamma
     $     SFR_Z(4000),         ! reported Z 
     $     SFR_mT2(4000),       ! reported mT2
     $     SFR_Flag(4000),      ! Flag: JS*(1..4)
     $     SFR_Flag1(4000),     ! Flag-Counter: added 1,if Flag is set
     $     SFR_REV(4000),       ! reverse ordering (unused)
     $     SFR_REV2(4000),      ! reverse ordering (unused)
     $     SFR_Nvec(4000),      ! number of vectors per particle
     $     SFR_Coeff(4000,SFR_NvecMax), ! coeff of vector
     $     SFR_Vec(4,4000,SFR_NvecMax), ! vector
     $     SFR_VecID(4000,SFR_NvecMax), ! ID of vector (cf. SFREPS)
     $     SFR_JJ(2,4000),      ! used for string-break ordering
     $     SFR_KK(2,4000),      !     -"-
     $     SFR_ID0(4000),       !     -"-
     $     SFR_DHM(4,4000),     ! reported coeffs for M-calc
     $     SFR_DHG(4,4000),     ! reported coeffs for Gamma-calc
     $     SFR_FB(4000),        ! reported parameter FB (???)
     $     SFR_GamIN(2,4000),   ! reported IN(1),IN(2)
     $     SFR_GamC(20,4000),   ! reported Coeffs for Gamma
     $     SFR_Ntot,            ! #(used entries) = #particles
     $     SFR_Nused(2)         ! indizes of lines which has been used
                                !                (min,max)

      integer SFR_Ntot
      integer SFR_Nused
      double precision SFR_Gamma,SFR_Z, SFR_mT2
      integer SFR_Flag, SFR_Flag1, SFR_REV, SFR_REV2
      integer SFR_Nvec, SFR_VecID
      double precision SFR_Coeff, SFR_Vec
      integer SFR_JJ,SFR_KK,SFR_ID0
      double precision SFR_DHM, SFR_DHG
      double precision SFR_FB
      integer SFR_GamIN
      double precision SFR_GamC

      save /DataSFREP/
      

      integer iL, iLL, iL4, iV, j
      double precision CC(SFR_NvecMax)
      integer iString, iVec
      data iString /0/
      integer ID,ID0
      double precision VV(4)
      integer JJ(2),KK(2), ILSB1, ILSB2
      logical ccOK

      integer SFREPS_NVec       ! prototype

c      write(*,*) '>> SFREP_ReconLastQgQ(',iLine1,iLine2,')'

      do iL=iLine1,iLine2
         if (abs(SFR_Flag(iL)).eq.4) iL4 = iL
      enddo

c... Transverse momenta: (save room for momenta)

      SFR_Nvec(iL4) = 4
      call SFREP_ResetVec(iL4,1)
      call SFREP_ResetVec(iL4,2)
      call SFREP_ResetVec(iL4,3)
      call SFREP_ResetVec(iL4,4)

c... Longitudinal momenta:
c... Get Coefficients:

      do iV=1,SFR_NvecMax
         CC(iV) = 0d0
      enddo

      do iL=iLine1,iLine2
         do iV=5,SFR_NVec(iL)
            call SFREPS_FindVecI(2, iL,iV, iString,iVec)
            if (iString.ne.0) then
               CC(iVec) = CC(iVec)+SFR_Coeff(iL,iV)
            endif
         enddo
      enddo

c... Add Longitudinal Momenta:

      do iV=1,SFREPS_NVec(2,iString)
         if (abs(CC(iV)-1d0).gt.1d-10) then
            call SFREPS_GetVec(2,iString,iV, VV,ID)
            call SFREP_AddVec0(iL4, 1d0-CC(iV), VV,ID)
         endif
      enddo

c... Set StringBreak:

      iLL = iLine1
      if (iL4.eq.iLL) iLL = iLine2
      call SFREPS_FindVecI(2, iLL,5, iString,iVec) ! to set iString
      call SFREPS_GetVec(2,iString,1, VV,ID0) ! to set ID0

      do iL=iLine1,iLine2
         JJ(1)=20000
         JJ(2)=-1
         KK(1)=20000
         KK(2)=-1

         do iV=5,SFR_NVec(iL)
            ID = SFR_VecID(iL,iV)
            if (mod((ID-ID0),2).eq.1) then
               if (ID.lt.KK(1)) KK(1)=SFR_VecID(iL,iV)
               if (ID.gt.KK(2)) KK(2)=SFR_VecID(iL,iV)
            else
               if (ID.lt.JJ(1)) JJ(1)=SFR_VecID(iL,iV)
               if (ID.gt.JJ(2)) JJ(2)=SFR_VecID(iL,iV)
            endif
         enddo
         SFR_ID0(iL) = ID0
         do j=1,2
            SFR_KK(j,iL) = KK(j)
            SFR_JJ(j,iL) = JJ(j)
         enddo
      enddo

c... Check String Breaking ordering:

      ILSB1 = 0                 ! from Rand, Default
      if (iL4.eq.iLine1) goto 701

      ILSB1 = iL4-1             ! SB1 == SB2(letzte Zeile)
      if (SFR_KK(1,iL4).eq.SFR_KK(2,iLSB1) .and.
     $     SFR_JJ(1,iL4).eq.SFR_JJ(2,iLSB1)) goto 701

      if (iL4.gt.iLine1+1) then
         ILSB1 = iL4-2          ! SB1 == SB2(vorletzte Zeile)
         if (SFR_KK(1,iL4).eq.SFR_KK(2,iLSB1) .and.
     $        SFR_JJ(1,iL4).eq.SFR_JJ(2,iLSB1)) goto 701
      else
         ILSB1 = 0              ! from Rand
         goto 701
      endif
      call SFREP_ERROR(6,'SFREP_ReconLastQgQ:')
      write(*,1000) iL4,'OOOps, SB1.'
      call SFREP_Write(iLine1,iLine2)
      return

 701  continue

      ILSB2 = 0                 ! from Rand, Default
      if (iL4.eq.iLine2) goto 702
      
      iLSB2 = iL4+1             ! SB2 == SB1(nächste Zeile)
      if (SFR_KK(2,iL4).eq.SFR_KK(1,iLSB2) .and.
     $     SFR_JJ(2,iL4).eq.SFR_JJ(1,iLSB2)) goto 702

      if (iL4.lt.iLine2-1) then
         iLSB2 = iL4+2          ! SB2 == SB1(übernächste Zeile)
         if (SFR_KK(2,iL4).eq.SFR_KK(1,iLSB2) .and.
     $        SFR_JJ(2,iL4).eq.SFR_JJ(1,iLSB2)) goto 702
      else
         ILSB2 = 0              ! from Rand
         goto 702
      endif
      call SFREP_ERROR(7,'SFREP_ReconLastQgQ:')
      write(*,1000) iL4,'OOOps, SB2.'
      call SFREP_Write(iLine1,iLine2)
      return

 702  continue

c... Set Reverse Ordering Flag:

c... may happen for: 2 Childs, reverse ordering:
      if (iLSB1.eq.0 .and. iLSB2.eq.0) then 
         SFR_REV(iL4) = 1
      endif

      if (iLSB1.eq.iL4-2 .or. iLSB2.eq.iL4+2) then
         SFR_REV(iL4) = 1
         if (abs(SFR_Flag(iL4-1)).eq.3) SFR_REV(iL4-1) = 1
         if (abs(SFR_Flag(iL4+1)).eq.3) SFR_REV(iL4+1) = 1
      endif

c... NOW: Transverse momenta:

      if (SFR_REV(iL4) .eq. 1) then
         if (iL4.eq.iLine2) ILSB2 = iL4-1
         if (iL4.eq.iLine1) ILSB1 = iL4+1
      endif

      
      if (SFR_FLag(iL4).eq.4) then
         if (iLSB1.ne.0) then
            call SFREP_CopyVec(iL4,1, iLSB1,2, -1d0)
            call SFREP_CopyVec(iL4,3, iLSB1,4, -1d0) 
         endif
         
         if (iLSB2.ne.0) then
            call SFREP_CopyVec(iL4,2, iLSB2,2, -1d0)
            call SFREP_CopyVec(iL4,4, iLSB2,4, -1d0)
         endif

         if (iLSB1.ne.iL4-1) then
            if (iLSB2.ne.iL4-1) then
               SFR_Coeff(iL4,2) = SFR_Coeff(iL4,2)
     $              -SFR_Coeff(iL4-1,2)
               SFR_Coeff(iL4,4) = SFR_Coeff(iL4,4)
     $              -SFR_Coeff(iL4-1,4)
            endif

            SFR_Coeff(iL4,2) = SFR_Coeff(iL4,2)
     $           -SFR_Coeff(iL4-1,1)
            SFR_Coeff(iL4,4) = SFR_Coeff(iL4,4)
     $           -SFR_Coeff(iL4-1,3)
         endif
      else
         if (iLSB2.ne.0) then
            call SFREP_CopyVec(iL4,1, iLSB2,2, -1d0)
            call SFREP_CopyVec(iL4,3, iLSB2,4, -1d0)
         endif

         if (iLSB1.ne.0) then
            call SFREP_CopyVec(iL4,2, iLSB1,2, -1d0)
            call SFREP_CopyVec(iL4,4, iLSB1,4, -1d0)
         endif

         if (iLSB2.ne.iL4+1) then
            if (iLSB1.ne.iL4+1) then
               SFR_Coeff(iL4,2) = SFR_Coeff(iL4,2)
     $              -SFR_Coeff(iL4+1,2)
               SFR_Coeff(iL4,4) = SFR_Coeff(iL4,4)
     $              -SFR_Coeff(iL4+1,4)
            endif
            SFR_Coeff(iL4,2) = SFR_Coeff(iL4,2)
     $           -SFR_Coeff(iL4+1,1)
            SFR_Coeff(iL4,4) = SFR_Coeff(iL4,4)
     $           -SFR_Coeff(iL4+1,3)
         endif
      endif

 900  continue

c... Check Correction:

      call SFREP_CheckCorr('SFREP_ReconLastQgQ',
     $     iLine1,iLine2,iL4,iString,ccOK)

      if (ccOK) then
      else
         call SFREP_ERROR(8,'SFREP_ReconLastQgQ:')
         write(*,1001) iL4,'Problem ', iL4,iLine1,iLine2,
     $        iLSB1,iLSB2,SFR_FLag(iL4),SFR_REV(iL4)
         call SFREP_Write(iLine1,iLine2)
      endif

 1000 FORMAT('SFREP_ReconLastQgQ:',i3,' ',A,1P,4e13.5)
 1001 FORMAT('SFREP_ReconLastQgQ:',i3,' ',A,8i5)

      end      

c====================================================
c... Reconstruct "last hadron" for a GG-String
c
c Uses the reported values z_(h,i)^+ and z_(h,i)^- to set up the vectors
c for the final hadron (flag=4).
c
c SORRY: Not possible.
c only Flags are set correctly, last hadron ist still not
c reconstructed!!!
c
c Erklärung des Problems:
c Ein String aus 2 Gluonen wird in 6 (!) Vektoren aufgesplittet, die die
c String-Regionen aufspannen. Summiere ich die Vorfaktoren dieser
c Vektoren, die alle bis auf das letzte Teilchen erzeugen, sollten
c bis auf maximal 4 Faktoren alle gleich 1 sein. Diese Vektoren
c erzeugen dann das letzte Teilchen. (Analoges sollte für die
c Transversalimpulse gelten.)
c Dies ist hier aber nicht der Fall: Wenn 6 der 6 Vorfaktoren != 1
c sind, kann ich nicht ersehen, wo die beiden Stringbreaks
c stattgefunden haben. (Ich weiß es zwar, da ich die Stringbreaks
c des linken und rechten Teilchens kenne, eine Rekonstruktion der
c Impulse ist so aber nicht möglich!)

      subroutine SFREP_ReconLastGG(iLine1,iLine2)
      IMPLICIT NONE
      integer iLine1,iLine2

      integer SFR_NvecMax
      parameter (SFR_NvecMax=100) ! former value: 40

      common /DataSFREP/
     $     SFR_Gamma(4000),     ! reported Gamma
     $     SFR_Z(4000),         ! reported Z 
     $     SFR_mT2(4000),       ! reported mT2
     $     SFR_Flag(4000),      ! Flag: JS*(1..4)
     $     SFR_Flag1(4000),     ! Flag-Counter: added 1,if Flag is set
     $     SFR_REV(4000),       ! reverse ordering (unused)
     $     SFR_REV2(4000),      ! reverse ordering (unused)
     $     SFR_Nvec(4000),      ! number of vectors per particle
     $     SFR_Coeff(4000,SFR_NvecMax), ! coeff of vector
     $     SFR_Vec(4,4000,SFR_NvecMax), ! vector
     $     SFR_VecID(4000,SFR_NvecMax), ! ID of vector (cf. SFREPS)
     $     SFR_JJ(2,4000),      ! used for string-break ordering
     $     SFR_KK(2,4000),      !     -"-
     $     SFR_ID0(4000),       !     -"-
     $     SFR_DHM(4,4000),     ! reported coeffs for M-calc
     $     SFR_DHG(4,4000),     ! reported coeffs for Gamma-calc
     $     SFR_FB(4000),        ! reported parameter FB (???)
     $     SFR_GamIN(2,4000),   ! reported IN(1),IN(2)
     $     SFR_GamC(20,4000),   ! reported Coeffs for Gamma
     $     SFR_Ntot,            ! #(used entries) = #particles
     $     SFR_Nused(2)         ! indizes of lines which has been used
                                !                (min,max)

      integer SFR_Ntot
      integer SFR_Nused
      double precision SFR_Gamma,SFR_Z, SFR_mT2
      integer SFR_Flag, SFR_Flag1, SFR_REV, SFR_REV2
      integer SFR_Nvec, SFR_VecID
      double precision SFR_Coeff, SFR_Vec
      integer SFR_JJ,SFR_KK,SFR_ID0
      double precision SFR_DHM, SFR_DHG
      double precision SFR_FB
      integer SFR_GamIN
      double precision SFR_GamC

      save /DataSFREP/
      

      integer iL, iLL, iL4, iV, j
      double precision CC(SFR_NvecMax)
      integer iString, iVec
      data iString /0/
      integer ID,ID0
      double precision VV(4)
      integer JJ(2),KK(2), ILSB1, ILSB2
      logical ccOK

      integer SFREPS_NVec       ! prototype

c$$$      write(*,*) '>> SFREP_ReconLastGG(',iLine1,iLine2,')'
c$$$      call PYLIST(2)
c$$$      call SFREPS_Write(1)
c$$$      call SFREPS_Write(2)
c$$$      call SFREPS_Write(3)
c$$$c      call SFREP_Write(iLine1,iLine2)

      do iL=iLine1,iLine2
         if (abs(SFR_Flag(iL)).eq.4) iL4 = iL
      enddo

c... Set StringBreak:

      iLL = iLine1
      if (iL4.eq.iLL) iLL = iLine2
      call SFREPS_FindVecI(2, iLL,5, iString,iVec) ! to set iString
      call SFREPS_GetVec(2,iString,1, VV,ID0)      ! to set ID0

      do iL=iLine1,iLine2
         if (iL.eq.iL4) goto 100
         JJ(1)=20000
         JJ(2)=-1
         KK(1)=20000
         KK(2)=-1

         do iV=5,SFR_NVec(iL)
            ID = SFR_VecID(iL,iV)
            if (mod((ID-ID0),2).eq.1) then
               if (ID.lt.KK(1)) KK(1)=SFR_VecID(iL,iV)
               if (ID.gt.KK(2)) KK(2)=SFR_VecID(iL,iV)
            else
               if (ID.lt.JJ(1)) JJ(1)=SFR_VecID(iL,iV)
               if (ID.gt.JJ(2)) JJ(2)=SFR_VecID(iL,iV)
            endif
         enddo
         SFR_ID0(iL) = ID0
         do j=1,2
            SFR_KK(j,iL) = KK(j)
            SFR_JJ(j,iL) = JJ(j)
         enddo
 100     continue
      enddo

c$$$      call SFREP_Write(iLine1,iLine2)
c$$$      stop

      return

c...
c... WHAT FOLLOWS IS NOT USED
c...

c... Transverse momenta: (save room for momenta)

      SFR_Nvec(iL4) = 4
      call SFREP_ResetVec(iL4,1)
      call SFREP_ResetVec(iL4,2)
      call SFREP_ResetVec(iL4,3)
      call SFREP_ResetVec(iL4,4)

c... Longitudinal momenta:
c... Get Coefficients:

      do iV=1,SFR_NvecMax
         CC(iV) = 0d0
      enddo

      do iL=iLine1,iLine2
         do iV=5,SFR_NVec(iL)
            call SFREPS_FindVecI(2, iL,iV, iString,iVec)
            if (iString.ne.0) then
               CC(iVec) = CC(iVec)+SFR_Coeff(iL,iV)
            endif
         enddo
      enddo

      do iV=1,SFR_NvecMax
         write(*,*) iV,CC(iV)
      enddo

      stop

c... Add Longitudinal Momenta:

      do iV=1,SFREPS_NVec(2,iString)
         if (abs(CC(iV)-1d0).gt.1d-10) then
            call SFREPS_GetVec(2,iString,iV, VV,ID)
            call SFREP_AddVec0(iL4, 1d0-CC(iV), VV,ID)
         endif
      enddo

c... Set StringBreak:

      iLL = iLine1
      if (iL4.eq.iLL) iLL = iLine2
      call SFREPS_FindVecI(2, iLL,5, iString,iVec) ! to set iString
      call SFREPS_GetVec(2,iString,1, VV,ID0)      ! to set ID0

      do iL=iLine1,iLine2
         JJ(1)=20000
         JJ(2)=-1
         KK(1)=20000
         KK(2)=-1

         do iV=5,SFR_NVec(iL)
            ID = SFR_VecID(iL,iV)
            if (mod((ID-ID0),2).eq.1) then
               if (ID.lt.KK(1)) KK(1)=SFR_VecID(iL,iV)
               if (ID.gt.KK(2)) KK(2)=SFR_VecID(iL,iV)
            else
               if (ID.lt.JJ(1)) JJ(1)=SFR_VecID(iL,iV)
               if (ID.gt.JJ(2)) JJ(2)=SFR_VecID(iL,iV)
            endif
         enddo
         SFR_ID0(iL) = ID0
         do j=1,2
            SFR_KK(j,iL) = KK(j)
            SFR_JJ(j,iL) = JJ(j)
         enddo
      enddo

c... Check String Breaking ordering:

      ILSB1 = 0                 ! from Rand, Default
      if (iL4.eq.iLine1) goto 701

      ILSB1 = iL4-1             ! SB1 == SB2(letzte Zeile)
      if (SFR_KK(1,iL4).eq.SFR_KK(2,iLSB1) .and.
     $     SFR_JJ(1,iL4).eq.SFR_JJ(2,iLSB1)) goto 701

      if (iL4.gt.iLine1+1) then
         ILSB1 = iL4-2          ! SB1 == SB2(vorletzte Zeile)
         if (SFR_KK(1,iL4).eq.SFR_KK(2,iLSB1) .and.
     $        SFR_JJ(1,iL4).eq.SFR_JJ(2,iLSB1)) goto 701
      else
         ILSB1 = 0              ! from Rand
         goto 701
      endif
c      call SFREP_ERROR(6,'SFREP_ReconLastGG:')
c      write(*,1000) iL4,'OOOps, SB1.'
c      call SFREP_Write(iLine1,iLine2)
      return

 701  continue

      ILSB2 = 0                 ! from Rand, Default
      if (iL4.eq.iLine2) goto 702
      
      iLSB2 = iL4+1             ! SB2 == SB1(nächste Zeile)
      if (SFR_KK(2,iL4).eq.SFR_KK(1,iLSB2) .and.
     $     SFR_JJ(2,iL4).eq.SFR_JJ(1,iLSB2)) goto 702

      if (iL4.lt.iLine2-1) then
         iLSB2 = iL4+2          ! SB2 == SB1(übernächste Zeile)
         if (SFR_KK(2,iL4).eq.SFR_KK(1,iLSB2) .and.
     $        SFR_JJ(2,iL4).eq.SFR_JJ(1,iLSB2)) goto 702
      else
         ILSB2 = 0              ! from Rand
         goto 702
      endif
c      call SFREP_ERROR(7,'SFREP_ReconLastGG:')
c      write(*,1000) iL4,'OOOps, SB2.'
c      call SFREP_Write(iLine1,iLine2)
      return

 702  continue

c... Set Reverse Ordering Flag:

c... may happen for: 2 Childs, reverse ordering:
      if (iLSB1.eq.0 .and. iLSB2.eq.0) then 
         SFR_REV(iL4) = 1
      endif

      if (iLSB1.eq.iL4-2 .or. iLSB2.eq.iL4+2) then
         SFR_REV(iL4) = 1
         if (abs(SFR_Flag(iL4-1)).eq.3) SFR_REV(iL4-1) = 1
         if (abs(SFR_Flag(iL4+1)).eq.3) SFR_REV(iL4+1) = 1
      endif

c... NOW: Transverse momenta:

      if (SFR_REV(iL4) .eq. 1) then
         if (iL4.eq.iLine2) ILSB2 = iL4-1
         if (iL4.eq.iLine1) ILSB1 = iL4+1
      endif

      
      if (SFR_FLag(iL4).eq.4) then
         if (iLSB1.ne.0) then
            call SFREP_CopyVec(iL4,1, iLSB1,2, -1d0)
            call SFREP_CopyVec(iL4,3, iLSB1,4, -1d0) 
         endif
         
         if (iLSB2.ne.0) then
            call SFREP_CopyVec(iL4,2, iLSB2,2, -1d0)
            call SFREP_CopyVec(iL4,4, iLSB2,4, -1d0)
         endif

         if (iLSB1.ne.iL4-1) then
            if (iLSB2.ne.iL4-1) then
               SFR_Coeff(iL4,2) = SFR_Coeff(iL4,2)
     $              -SFR_Coeff(iL4-1,2)
               SFR_Coeff(iL4,4) = SFR_Coeff(iL4,4)
     $              -SFR_Coeff(iL4-1,4)
            endif

            SFR_Coeff(iL4,2) = SFR_Coeff(iL4,2)
     $           -SFR_Coeff(iL4-1,1)
            SFR_Coeff(iL4,4) = SFR_Coeff(iL4,4)
     $           -SFR_Coeff(iL4-1,3)
         endif
      else
         if (iLSB2.ne.0) then
            call SFREP_CopyVec(iL4,1, iLSB2,2, -1d0)
            call SFREP_CopyVec(iL4,3, iLSB2,4, -1d0)
         endif

         if (iLSB1.ne.0) then
            call SFREP_CopyVec(iL4,2, iLSB1,2, -1d0)
            call SFREP_CopyVec(iL4,4, iLSB1,4, -1d0)
         endif

         if (iLSB2.ne.iL4+1) then
            if (iLSB1.ne.iL4+1) then
               SFR_Coeff(iL4,2) = SFR_Coeff(iL4,2)
     $              -SFR_Coeff(iL4+1,2)
               SFR_Coeff(iL4,4) = SFR_Coeff(iL4,4)
     $              -SFR_Coeff(iL4+1,4)
            endif
            SFR_Coeff(iL4,2) = SFR_Coeff(iL4,2)
     $           -SFR_Coeff(iL4+1,1)
            SFR_Coeff(iL4,4) = SFR_Coeff(iL4,4)
     $           -SFR_Coeff(iL4+1,3)
         endif
      endif

 900  continue

c... Check Correction:

      call SFREP_CheckCorr('SFREP_ReconLastGG',
     $     iLine1,iLine2,iL4,iString,ccOK)

      if (ccOK) then
      else
c         call SFREP_ERROR(8,'SFREP_ReconLastGG:')
c         write(*,1001) iL4,'Problem ', iL4,iLine1,iLine2,
c     $        iLSB1,iLSB2,SFR_FLag(iL4),SFR_REV(iL4)
c         call SFREP_Write(iLine1,iLine2)
      endif

 1000 FORMAT('SFREP_ReconLastGG:',i3,' ',A,1P,4e13.5)
 1001 FORMAT('SFREP_ReconLastGG:',i3,' ',A,8i5)

      end      

c====================================================
c check the correction
c
      subroutine SFREP_CheckCorr(S,iLine1,iLine2,iLL,iString, ccOK)
      IMPLICIT NONE
      character S*(*)
      integer iLine1,iLine2,iLL,iString
      logical ccOK

      integer SFR_NvecMax
      parameter (SFR_NvecMax=100) ! former value: 40

      common /DataSFREP/
     $     SFR_Gamma(4000),     ! reported Gamma
     $     SFR_Z(4000),         ! reported Z 
     $     SFR_mT2(4000),       ! reported mT2
     $     SFR_Flag(4000),      ! Flag: JS*(1..4)
     $     SFR_Flag1(4000),     ! Flag-Counter: added 1,if Flag is set
     $     SFR_REV(4000),       ! reverse ordering (unused)
     $     SFR_REV2(4000),      ! reverse ordering (unused)
     $     SFR_Nvec(4000),      ! number of vectors per particle
     $     SFR_Coeff(4000,SFR_NvecMax), ! coeff of vector
     $     SFR_Vec(4,4000,SFR_NvecMax), ! vector
     $     SFR_VecID(4000,SFR_NvecMax), ! ID of vector (cf. SFREPS)
     $     SFR_JJ(2,4000),      ! used for string-break ordering
     $     SFR_KK(2,4000),      !     -"-
     $     SFR_ID0(4000),       !     -"-
     $     SFR_DHM(4,4000),     ! reported coeffs for M-calc
     $     SFR_DHG(4,4000),     ! reported coeffs for Gamma-calc
     $     SFR_FB(4000),        ! reported parameter FB (???)
     $     SFR_GamIN(2,4000),   ! reported IN(1),IN(2)
     $     SFR_GamC(20,4000),   ! reported Coeffs for Gamma
     $     SFR_Ntot,            ! #(used entries) = #particles
     $     SFR_Nused(2)         ! indizes of lines which has been used
                                !                (min,max)

      integer SFR_Ntot
      integer SFR_Nused
      double precision SFR_Gamma,SFR_Z, SFR_mT2
      integer SFR_Flag, SFR_Flag1, SFR_REV, SFR_REV2
      integer SFR_Nvec, SFR_VecID
      double precision SFR_Coeff, SFR_Vec
      integer SFR_JJ,SFR_KK,SFR_ID0
      double precision SFR_DHM, SFR_DHG
      double precision SFR_FB
      integer SFR_GamIN
      double precision SFR_GamC

      save /DataSFREP/
      

      integer iL,iV, j
      double precision V(4),VV(4),h

c... Check Correction: V_T

      ccOK = .TRUE.

      do j=1,4
         V(j) = 0d0
      enddo

      do iL=iLine1,iLine2
         do iV=1,4              ! only transverse momenta
            do j=1,4
               V(j) = V(j)+SFR_Coeff(iL,iV)*SFR_Vec(j,iL,iV)
            enddo
         enddo
      enddo

      if (V(1)**2+V(2)**2+V(3)**2+V(4)**2.gt.1d-15) then
                                ! not FOUR-Product !
cc         write(*,1000) S,iLL,'V_T=',(V(j),j=1,4)
         ccOK = .FALSE.
      endif

c... Check Correction: V

      do iL=iLine1,iLine2
         do j=1,4
            V(j) = 0d0
         enddo
         do iV=1,SFR_NVec(iL)              ! all momenta
            do j=1,4
               V(j) = V(j)+SFR_Coeff(iL,iV)*SFR_Vec(j,iL,iV)
            enddo
         enddo

         call SFREPS_GetVecID(3,iString,iV, VV,iL)
         

         h = 0d0
         do j=1,4
            h = h + (V(j)-VV(j))**2 ! not Four-Product
         enddo

         if (h.gt.1d-15) then
cc            call SFREP_ERROR(9,'SFREP_CheckCorr:')
c            write(*,1001) S,iL,'',iLine1,iLine2,iLL,iString
cc            write(*,1000) S,iL,'Problem: V =',(V(j),j=1,4)
cc            write(*,1000) S,iL,'Problem: VV=',(VV(j),j=1,4)
c            call SFREP_Write(iLine1,iLine2)
c            call PYLIST(2)
c            write(*,*) '====='
            ccOK = .FALSE.
         endif
      enddo

 1000 FORMAT(A,':',i3,' [CC]: ',A,1P,4e13.5)
 1001 FORMAT(A,':',i3,' [CC]: ',A,4i5)

      end

c====================================================
c Set vector iV1 from iLine as vector iV2 for iLine2
c
c sets the coefficient of (iV1,iLine1) Coeff(1,1) to c*Coeff(2,2)
c
      subroutine SFREP_CopyVec(iL1,iV1, iL2,iV2, c)
      IMPLICIT NONE
      integer iL1,iV1, iL2,iV2
      double precision c

      integer SFR_NvecMax
      parameter (SFR_NvecMax=100) ! former value: 40

      common /DataSFREP/
     $     SFR_Gamma(4000),     ! reported Gamma
     $     SFR_Z(4000),         ! reported Z 
     $     SFR_mT2(4000),       ! reported mT2
     $     SFR_Flag(4000),      ! Flag: JS*(1..4)
     $     SFR_Flag1(4000),     ! Flag-Counter: added 1,if Flag is set
     $     SFR_REV(4000),       ! reverse ordering (unused)
     $     SFR_REV2(4000),      ! reverse ordering (unused)
     $     SFR_Nvec(4000),      ! number of vectors per particle
     $     SFR_Coeff(4000,SFR_NvecMax), ! coeff of vector
     $     SFR_Vec(4,4000,SFR_NvecMax), ! vector
     $     SFR_VecID(4000,SFR_NvecMax), ! ID of vector (cf. SFREPS)
     $     SFR_JJ(2,4000),      ! used for string-break ordering
     $     SFR_KK(2,4000),      !     -"-
     $     SFR_ID0(4000),       !     -"-
     $     SFR_DHM(4,4000),     ! reported coeffs for M-calc
     $     SFR_DHG(4,4000),     ! reported coeffs for Gamma-calc
     $     SFR_FB(4000),        ! reported parameter FB (???)
     $     SFR_GamIN(2,4000),   ! reported IN(1),IN(2)
     $     SFR_GamC(20,4000),   ! reported Coeffs for Gamma
     $     SFR_Ntot,            ! #(used entries) = #particles
     $     SFR_Nused(2)         ! indizes of lines which has been used
                                !                (min,max)

      integer SFR_Ntot
      integer SFR_Nused
      double precision SFR_Gamma,SFR_Z, SFR_mT2
      integer SFR_Flag, SFR_Flag1, SFR_REV, SFR_REV2
      integer SFR_Nvec, SFR_VecID
      double precision SFR_Coeff, SFR_Vec
      integer SFR_JJ,SFR_KK,SFR_ID0
      double precision SFR_DHM, SFR_DHG
      double precision SFR_FB
      integer SFR_GamIN
      double precision SFR_GamC

      save /DataSFREP/
      

      integer j

      SFR_Coeff(iL1,iV1) = c*SFR_Coeff(iL2,iV2)
      do j=1,4
         SFR_Vec(j,iL1,iV1) = SFR_Vec(j,iL2,iV2)
      enddo
      SFR_VecID(iL1,iV1) = SFR_VecID(iL2,iV2)
      end

c====================================================
c Set the parameters necessary to recalculate DHG
c
c ***** used by PYSTRFT.FT.F *****
c
      subroutine SFREP_SetGammaC(iLine,IN1,IN2,JT)
      IMPLICIT NONE
      integer iLine, IN1, IN2, JT

      integer SFR_NvecMax
      parameter (SFR_NvecMax=100) ! former value: 40

      common /DataSFREP/
     $     SFR_Gamma(4000),     ! reported Gamma
     $     SFR_Z(4000),         ! reported Z 
     $     SFR_mT2(4000),       ! reported mT2
     $     SFR_Flag(4000),      ! Flag: JS*(1..4)
     $     SFR_Flag1(4000),     ! Flag-Counter: added 1,if Flag is set
     $     SFR_REV(4000),       ! reverse ordering (unused)
     $     SFR_REV2(4000),      ! reverse ordering (unused)
     $     SFR_Nvec(4000),      ! number of vectors per particle
     $     SFR_Coeff(4000,SFR_NvecMax), ! coeff of vector
     $     SFR_Vec(4,4000,SFR_NvecMax), ! vector
     $     SFR_VecID(4000,SFR_NvecMax), ! ID of vector (cf. SFREPS)
     $     SFR_JJ(2,4000),      ! used for string-break ordering
     $     SFR_KK(2,4000),      !     -"-
     $     SFR_ID0(4000),       !     -"-
     $     SFR_DHM(4,4000),     ! reported coeffs for M-calc
     $     SFR_DHG(4,4000),     ! reported coeffs for Gamma-calc
     $     SFR_FB(4000),        ! reported parameter FB (???)
     $     SFR_GamIN(2,4000),   ! reported IN(1),IN(2)
     $     SFR_GamC(20,4000),   ! reported Coeffs for Gamma
     $     SFR_Ntot,            ! #(used entries) = #particles
     $     SFR_Nused(2)         ! indizes of lines which has been used
                                !                (min,max)

      integer SFR_Ntot
      integer SFR_Nused
      double precision SFR_Gamma,SFR_Z, SFR_mT2
      integer SFR_Flag, SFR_Flag1, SFR_REV, SFR_REV2
      integer SFR_Nvec, SFR_VecID
      double precision SFR_Coeff, SFR_Vec
      integer SFR_JJ,SFR_KK,SFR_ID0
      double precision SFR_DHM, SFR_DHG
      double precision SFR_FB
      integer SFR_GamIN
      double precision SFR_GamC

      save /DataSFREP/
      

      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      integer N,NPAD,K
      double precision P,V
      SAVE /PYJETS/

      integer IN,i

c      write(*,*) '==> SFREP_SetGammaC(',iLine,IN1,IN2,JT,')'

      SFR_GamIN(1,iLine) = IN1
      SFR_GamIN(2,iLine) = IN2

      do i=1,10
         SFR_GamC(i,iLine) = 0d0
      enddo

      if ((IN2-1-IN1)/4+1.gt.10) then
         call SFREP_ERROR(2,'SFREP_SetGammaC: more than 20 Coeff!')
         return
      endif
      
      i=0
      do IN=IN1,IN2-1,4
c         write(*,*) 'Setting...',i+1,IN,P(IN+2,JT)
c         write(*,*) 'Setting...',i+2,IN+1,P(IN+3,JT)

         SFR_GamC(i+1,iLine) = P(IN+2,JT)
         SFR_GamC(i+2,iLine) = P(IN+3,JT)
         i=i+2
      enddo

      end

c====================================================
c=== LOW LEVEL ROUTINES =============================
c====================================================

c===================================================================
c Report in Nused, which Lines have been used
c
      subroutine SFREP_SetUsed(iLine)
      IMPLICIT NONE
      integer iLine
      

      integer SFR_NvecMax
      parameter (SFR_NvecMax=100) ! former value: 40

      common /DataSFREP/
     $     SFR_Gamma(4000),     ! reported Gamma
     $     SFR_Z(4000),         ! reported Z 
     $     SFR_mT2(4000),       ! reported mT2
     $     SFR_Flag(4000),      ! Flag: JS*(1..4)
     $     SFR_Flag1(4000),     ! Flag-Counter: added 1,if Flag is set
     $     SFR_REV(4000),       ! reverse ordering (unused)
     $     SFR_REV2(4000),      ! reverse ordering (unused)
     $     SFR_Nvec(4000),      ! number of vectors per particle
     $     SFR_Coeff(4000,SFR_NvecMax), ! coeff of vector
     $     SFR_Vec(4,4000,SFR_NvecMax), ! vector
     $     SFR_VecID(4000,SFR_NvecMax), ! ID of vector (cf. SFREPS)
     $     SFR_JJ(2,4000),      ! used for string-break ordering
     $     SFR_KK(2,4000),      !     -"-
     $     SFR_ID0(4000),       !     -"-
     $     SFR_DHM(4,4000),     ! reported coeffs for M-calc
     $     SFR_DHG(4,4000),     ! reported coeffs for Gamma-calc
     $     SFR_FB(4000),        ! reported parameter FB (???)
     $     SFR_GamIN(2,4000),   ! reported IN(1),IN(2)
     $     SFR_GamC(20,4000),   ! reported Coeffs for Gamma
     $     SFR_Ntot,            ! #(used entries) = #particles
     $     SFR_Nused(2)         ! indizes of lines which has been used
                                !                (min,max)

      integer SFR_Ntot
      integer SFR_Nused
      double precision SFR_Gamma,SFR_Z, SFR_mT2
      integer SFR_Flag, SFR_Flag1, SFR_REV, SFR_REV2
      integer SFR_Nvec, SFR_VecID
      double precision SFR_Coeff, SFR_Vec
      integer SFR_JJ,SFR_KK,SFR_ID0
      double precision SFR_DHM, SFR_DHG
      double precision SFR_FB
      integer SFR_GamIN
      double precision SFR_GamC

      save /DataSFREP/
      

      if (iLine .lt. SFR_Nused(1)) SFR_Nused(1)=iLine
      if (iLine .gt. SFR_Nused(2)) SFR_Nused(2)=iLine

      end
c-------------------------------------------------------------------
c Reset the "Report in Nused,..."
c
c If you want mark all lines as empty, i.e. as unused 
c ... iL1 = 4001, iL2 = 0
c If you want to mark all lines to be deleted by the next call to
c SFREP_Reset, then set
c ... il1 = 1, iL2 = 4000

      subroutine SFREP_SetUsed0(allUsed)
      IMPLICIT NONE
      logical allUsed

      integer SFR_NvecMax
      parameter (SFR_NvecMax=100) ! former value: 40

      common /DataSFREP/
     $     SFR_Gamma(4000),     ! reported Gamma
     $     SFR_Z(4000),         ! reported Z 
     $     SFR_mT2(4000),       ! reported mT2
     $     SFR_Flag(4000),      ! Flag: JS*(1..4)
     $     SFR_Flag1(4000),     ! Flag-Counter: added 1,if Flag is set
     $     SFR_REV(4000),       ! reverse ordering (unused)
     $     SFR_REV2(4000),      ! reverse ordering (unused)
     $     SFR_Nvec(4000),      ! number of vectors per particle
     $     SFR_Coeff(4000,SFR_NvecMax), ! coeff of vector
     $     SFR_Vec(4,4000,SFR_NvecMax), ! vector
     $     SFR_VecID(4000,SFR_NvecMax), ! ID of vector (cf. SFREPS)
     $     SFR_JJ(2,4000),      ! used for string-break ordering
     $     SFR_KK(2,4000),      !     -"-
     $     SFR_ID0(4000),       !     -"-
     $     SFR_DHM(4,4000),     ! reported coeffs for M-calc
     $     SFR_DHG(4,4000),     ! reported coeffs for Gamma-calc
     $     SFR_FB(4000),        ! reported parameter FB (???)
     $     SFR_GamIN(2,4000),   ! reported IN(1),IN(2)
     $     SFR_GamC(20,4000),   ! reported Coeffs for Gamma
     $     SFR_Ntot,            ! #(used entries) = #particles
     $     SFR_Nused(2)         ! indizes of lines which has been used
                                !                (min,max)

      integer SFR_Ntot
      integer SFR_Nused
      double precision SFR_Gamma,SFR_Z, SFR_mT2
      integer SFR_Flag, SFR_Flag1, SFR_REV, SFR_REV2
      integer SFR_Nvec, SFR_VecID
      double precision SFR_Coeff, SFR_Vec
      integer SFR_JJ,SFR_KK,SFR_ID0
      double precision SFR_DHM, SFR_DHG
      double precision SFR_FB
      integer SFR_GamIN
      double precision SFR_GamC

      save /DataSFREP/
      

c      write(*,*) '==> SFREP_SetUsed0(',allUsed,')'

      if (allUsed) then
         SFR_Nused(1) = 1       ! = Minimum !!!
         SFR_Nused(2) = 4000    ! = Maximum !!!
      else
         SFR_Nused(1) = 4001    ! = Minimum !!!
         SFR_Nused(2) = 0       ! = Maximum !!!
      end if
      end
c-------------------------------------------------------------------
c Calculate, which lines have to be reset REALLY
c
c INPUTS
c * iL1,iL2: line range to be reset
c OUTPUT
c * iL1,iL2: line range to be reset
c * SFR_Nused(i) is reset
c 
      subroutine SFREP_CalcToReset(iL1,iL2)
      IMPLICIT NONE

      integer SFR_NvecMax
      parameter (SFR_NvecMax=100) ! former value: 40

      common /DataSFREP/
     $     SFR_Gamma(4000),     ! reported Gamma
     $     SFR_Z(4000),         ! reported Z 
     $     SFR_mT2(4000),       ! reported mT2
     $     SFR_Flag(4000),      ! Flag: JS*(1..4)
     $     SFR_Flag1(4000),     ! Flag-Counter: added 1,if Flag is set
     $     SFR_REV(4000),       ! reverse ordering (unused)
     $     SFR_REV2(4000),      ! reverse ordering (unused)
     $     SFR_Nvec(4000),      ! number of vectors per particle
     $     SFR_Coeff(4000,SFR_NvecMax), ! coeff of vector
     $     SFR_Vec(4,4000,SFR_NvecMax), ! vector
     $     SFR_VecID(4000,SFR_NvecMax), ! ID of vector (cf. SFREPS)
     $     SFR_JJ(2,4000),      ! used for string-break ordering
     $     SFR_KK(2,4000),      !     -"-
     $     SFR_ID0(4000),       !     -"-
     $     SFR_DHM(4,4000),     ! reported coeffs for M-calc
     $     SFR_DHG(4,4000),     ! reported coeffs for Gamma-calc
     $     SFR_FB(4000),        ! reported parameter FB (???)
     $     SFR_GamIN(2,4000),   ! reported IN(1),IN(2)
     $     SFR_GamC(20,4000),   ! reported Coeffs for Gamma
     $     SFR_Ntot,            ! #(used entries) = #particles
     $     SFR_Nused(2)         ! indizes of lines which has been used
                                !                (min,max)

      integer SFR_Ntot
      integer SFR_Nused
      double precision SFR_Gamma,SFR_Z, SFR_mT2
      integer SFR_Flag, SFR_Flag1, SFR_REV, SFR_REV2
      integer SFR_Nvec, SFR_VecID
      double precision SFR_Coeff, SFR_Vec
      integer SFR_JJ,SFR_KK,SFR_ID0
      double precision SFR_DHM, SFR_DHG
      double precision SFR_FB
      integer SFR_GamIN
      double precision SFR_GamC

      save /DataSFREP/
      

      integer iL1,iL2

      integer iiL1,iiL2, iiN1,iiN2
      
      iiL1 = iL1
      iiL2 = iL2
      iiN1 = SFR_Nused(1)
      iiN2 = SFR_Nused(2)

!      write(*,*) '==> SFREP_CalcToReset: ',iiL1,iiL2,iiN1,iiN2

!      return ! FOR TEMPORARY USE !!!!
      
      iL1 = 4001
      iL2 = 0

      
      if (iiL1.le.iiN1) then
         if (iiL2.lt.iiN1) return
         if (iiL2.le.iiN2) then
            iL1 = iiN1
            iL2 = iiL2
            SFR_Nused(1) = iiL2+1
         else
            iL1 = iiN1
            iL2 = iiN2
            SFR_Nused(1) = iiL2+1
            SFR_Nused(2) = iiL1-1
         endif
      else if (iiL1.le.iiN2) then
         if (iiL2.lt.iiN2) then
            iL1 = iiL1
            iL2 = iiL2
         else
            iL1 = iiL1
            iL2 = iiN2
            SFR_Nused(2) = iiL1-1
         endif
      endif
      
      if (SFR_Nused(1).gt.SFR_Nused(2)) then
         SFR_Nused(1)=4001
         SFR_Nused(2)=0
      endif
      
      
      end
c===================================================================

c===================================================================
c-------------------------------------------------------------------
c LOW LEVEL: Set an Array of doubles (size nA) to Val
c
      subroutine SetArrayD(nA, A, val)
      IMPLICIT NONE
      integer nA
      double precision A(*), val

      integer j

      do j=1,nA
         A(j) = val
      enddo
      end
c-------------------------------------------------------------------
c LOW LEVEL: Set an Array of integers (size nA) to Val
c
      subroutine SetArrayI(nA, A, val)
      IMPLICIT NONE
      integer nA
      integer A(*), val

      integer j

      do j=1,nA
         A(j) = val
      enddo
      end
c-------------------------------------------------------------------
c LOW LEVEL: Set an Array of size nA to Values of Array A2
c
      subroutine SetArrayArr(nA, A, A2)
      IMPLICIT NONE
      integer nA
      double precision A(*), A2(*)

      integer j

      do j=1,nA
         A(j) = A2(j)
      enddo
      end
c-------------------------------------------------------------------
c LOW LEVEL : Add to a Vector fak*VV
c
      subroutine AddVert(V,fak,VV)
      IMPLICIT NONE
      double precision V(4), fak, VV(4)

      integer j

      do j=1,4
         V(j) = V(j) + fak*VV(j)
      enddo
      end
c-------------------------------------------------------------------
c LOW LEVEL : Add to a Vector fak*SFR_Vec(...,iL,iV)
c
      subroutine AddVertSFR(V,fak,iL,iV)
      IMPLICIT NONE
      double precision V(4), fak
      integer iL,iV

      integer SFR_NvecMax
      parameter (SFR_NvecMax=100) ! former value: 40

      common /DataSFREP/
     $     SFR_Gamma(4000),     ! reported Gamma
     $     SFR_Z(4000),         ! reported Z 
     $     SFR_mT2(4000),       ! reported mT2
     $     SFR_Flag(4000),      ! Flag: JS*(1..4)
     $     SFR_Flag1(4000),     ! Flag-Counter: added 1,if Flag is set
     $     SFR_REV(4000),       ! reverse ordering (unused)
     $     SFR_REV2(4000),      ! reverse ordering (unused)
     $     SFR_Nvec(4000),      ! number of vectors per particle
     $     SFR_Coeff(4000,SFR_NvecMax), ! coeff of vector
     $     SFR_Vec(4,4000,SFR_NvecMax), ! vector
     $     SFR_VecID(4000,SFR_NvecMax), ! ID of vector (cf. SFREPS)
     $     SFR_JJ(2,4000),      ! used for string-break ordering
     $     SFR_KK(2,4000),      !     -"-
     $     SFR_ID0(4000),       !     -"-
     $     SFR_DHM(4,4000),     ! reported coeffs for M-calc
     $     SFR_DHG(4,4000),     ! reported coeffs for Gamma-calc
     $     SFR_FB(4000),        ! reported parameter FB (???)
     $     SFR_GamIN(2,4000),   ! reported IN(1),IN(2)
     $     SFR_GamC(20,4000),   ! reported Coeffs for Gamma
     $     SFR_Ntot,            ! #(used entries) = #particles
     $     SFR_Nused(2)         ! indizes of lines which has been used
                                !                (min,max)

      integer SFR_Ntot
      integer SFR_Nused
      double precision SFR_Gamma,SFR_Z, SFR_mT2
      integer SFR_Flag, SFR_Flag1, SFR_REV, SFR_REV2
      integer SFR_Nvec, SFR_VecID
      double precision SFR_Coeff, SFR_Vec
      integer SFR_JJ,SFR_KK,SFR_ID0
      double precision SFR_DHM, SFR_DHG
      double precision SFR_FB
      integer SFR_GamIN
      double precision SFR_GamC

      save /DataSFREP/
      

      integer j

      do j=1,4
         V(j) = V(j) + fak*SFR_Vec(j,iL,iV)
      enddo
      end
c-------------------------------------------------------------------
c LOW LEVEL : Length of Vektor
c
      subroutine VertLength(V,L)
      IMPLICIT NONE
      double precision V(4),L

      L = V(4)**2-V(1)**2-V(2)**2-V(3)**2
      return
      end

c===================================================================
