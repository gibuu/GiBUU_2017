c===================================================================
c HandleSFREP_S.F: Kai Gallmeister, 6.5.2004...30.01.2006
c
c Routines to handle the common block /DataSFREP_S/.
c Needed for calculation of production ans formation times from the
c (modified) JETSET-Fragmentation routine "PYSTRF.FT.F".
c
c Here overall (momentum) vectors are handled.
c
c===================================================================

c===================================================================
c Reset the full list for a given string.
c
c if the number of the string < 0, all strings are deleted.
c
      subroutine SFREPS_Reset(iString)
      IMPLICIT NONE
      integer iString

      integer SFRS_Nmax, SFRS_NVecmax
      parameter (SFRS_Nmax=100,SFRS_NVecmax=100) ! former value: 10,50

      common /DataSFREP_S/
     $     SFRS_N,               ! #Strings
     $     SFRS_NVec(3,SFRS_Nmax), ! #(InVectors) of String i
     $     SFRS_Vec(4,3,SFRS_Nmax, SFRS_NVecmax),
     $     SFRS_VecID(3,SFRS_Nmax, SFRS_NVecmax),
     $     SFRS_DPS(5,SFRS_Nmax) ! Boost-Vector of String

      integer SFRS_N, SFRS_NVec, SFRS_VecID
      double precision SFRS_Vec
      double precision SFRS_DPS

      save /DataSFREP_S/

      integer iS,iA

      if (iString.le.0) then
         SFRS_N = 0
         do iA=1,3
            do iS=1,SFRS_Nmax
               SFRS_NVec(iA,iS) = 0
            enddo
         enddo
      else
         do iA=1,3
            SFRS_NVec(iA,iString) = 0
         enddo
      endif
      end

c===================================================================
c Make up a new String. 
c
c Increase the internal counter of strings and reset all the values of
c the "new" string
c
c ***** used by PYSTRFT.FT.F *****
c
      subroutine SFREPS_NewString()
      IMPLICIT NONE

      integer SFRS_Nmax, SFRS_NVecmax
      parameter (SFRS_Nmax=100,SFRS_NVecmax=100) ! former value: 10,50

      common /DataSFREP_S/
     $     SFRS_N,               ! #Strings
     $     SFRS_NVec(3,SFRS_Nmax), ! #(InVectors) of String i
     $     SFRS_Vec(4,3,SFRS_Nmax, SFRS_NVecmax),
     $     SFRS_VecID(3,SFRS_Nmax, SFRS_NVecmax),
     $     SFRS_DPS(5,SFRS_Nmax) ! Boost-Vector of String

      integer SFRS_N, SFRS_NVec, SFRS_VecID
      double precision SFRS_Vec
      double precision SFRS_DPS

      save /DataSFREP_S/

      integer iA 

c      write(*,*) '==> SFREPS_NewString(',')'

      if (SFRS_N.eq.SFRS_Nmax) then
         call SFREP_ERROR(-1,'SFREPS_NewString: N exhausted!')
         stop
         return
      endif
      SFRS_N = SFRS_N + 1

      do iA=1,3
         SFRS_NVec(iA,SFRS_N) = 0
      enddo
      end

c===================================================================
c Add one momentum vector to the actual string.
c
c iA=1,2,3 indicates the internal sublist
c
      subroutine SFREPS_AddVec(iA,iP)
      IMPLICIT NONE
      integer iA,iP
      

      integer SFRS_Nmax, SFRS_NVecmax
      parameter (SFRS_Nmax=100,SFRS_NVecmax=100) ! former value: 10,50

      common /DataSFREP_S/
     $     SFRS_N,               ! #Strings
     $     SFRS_NVec(3,SFRS_Nmax), ! #(InVectors) of String i
     $     SFRS_Vec(4,3,SFRS_Nmax, SFRS_NVecmax),
     $     SFRS_VecID(3,SFRS_Nmax, SFRS_NVecmax),
     $     SFRS_DPS(5,SFRS_Nmax) ! Boost-Vector of String

      integer SFRS_N, SFRS_NVec, SFRS_VecID
      double precision SFRS_Vec
      double precision SFRS_DPS

      save /DataSFREP_S/

      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      integer N,NPAD,K
      double precision P,V
      SAVE /PYJETS/

      integer iV,j

      if (SFRS_N.lt.1) then
         call SFREP_ERROR(-2,'SFREPS_AddVec: SFRS_N < 1.')
         return
      endif

      if (SFRS_NVec(iA,SFRS_N).eq.SFRS_NVecmax) then
         call SFREP_ERROR(-3,'SFREPS_AddVec: Nvec() exhausted.')
         write(*,*) 'OOOps: iA,SFRS_N=',iA,SFRS_N
         return
      endif

      SFRS_NVec(iA,SFRS_N) = SFRS_NVec(iA,SFRS_N)+1
      iV = SFRS_NVec(iA,SFRS_N)

      do j=1,4
         SFRS_Vec(j,iA,SFRS_N,iV) = P(iP,j)
      enddo
      SFRS_VecID(iA,SFRS_N,iV) = iP

      end
c===================================================================
c Get the Number of stored vectors
c
c iA=1,2,3 indicates the internal sublist
c iS indicates the string number
c
      integer function SFREPS_NVec(iA,iS)
      IMPLICIT NONE
      integer iA,iS
      

      integer SFRS_Nmax, SFRS_NVecmax
      parameter (SFRS_Nmax=100,SFRS_NVecmax=100) ! former value: 10,50

      common /DataSFREP_S/
     $     SFRS_N,               ! #Strings
     $     SFRS_NVec(3,SFRS_Nmax), ! #(InVectors) of String i
     $     SFRS_Vec(4,3,SFRS_Nmax, SFRS_NVecmax),
     $     SFRS_VecID(3,SFRS_Nmax, SFRS_NVecmax),
     $     SFRS_DPS(5,SFRS_Nmax) ! Boost-Vector of String

      integer SFRS_N, SFRS_NVec, SFRS_VecID
      double precision SFRS_Vec
      double precision SFRS_DPS

      save /DataSFREP_S/

      SFREPS_NVec = SFRS_NVec(iA,iS)
      return
      end

c===================================================================
c Add a manifold of vectors (iP1..iP2) to the actual string.
c
c iA=1,2,3 indicates the internal sublist
c
c ***** used by PYSTRFT.FT.F *****
c
      subroutine SFREPS_AddVecs(iA,iP1,iP2)
      IMPLICIT NONE
      integer iA, iP1,IP2,iP

c      write(*,*) '==> SFREPS_AddVecs(',iA,iP1,iP2,')'
      
      do iP=iP1,iP2
         call SFREPS_AddVec(iA,iP)
      enddo
      end

c===================================================================
c Set the Boost-Vector of the actual string
c
c ***** used by PYSTRFT.FT.F *****
c
      subroutine SFREPS_SetBoost(DPS)
      IMPLICIT NONE
      double precision DPS(5)

      integer SFRS_Nmax, SFRS_NVecmax
      parameter (SFRS_Nmax=100,SFRS_NVecmax=100) ! former value: 10,50

      common /DataSFREP_S/
     $     SFRS_N,               ! #Strings
     $     SFRS_NVec(3,SFRS_Nmax), ! #(InVectors) of String i
     $     SFRS_Vec(4,3,SFRS_Nmax, SFRS_NVecmax),
     $     SFRS_VecID(3,SFRS_Nmax, SFRS_NVecmax),
     $     SFRS_DPS(5,SFRS_Nmax) ! Boost-Vector of String

      integer SFRS_N, SFRS_NVec, SFRS_VecID
      double precision SFRS_Vec
      double precision SFRS_DPS

      save /DataSFREP_S/

      integer j

c      write(*,*) '==> SFREPS_SetBoost(',(DPS(j),j=1,5),')'

      do j=1,5
         SFRS_DPS(j,SFRS_N) = DPS(j)
      enddo

      end

c===================================================================
c Get the Boost-Vector of string iS
c

      subroutine SFREPS_GetBoost(iS,DPS)
      IMPLICIT NONE
      integer iS
      double precision DPS(5)

      integer SFRS_Nmax, SFRS_NVecmax
      parameter (SFRS_Nmax=100,SFRS_NVecmax=100) ! former value: 10,50

      common /DataSFREP_S/
     $     SFRS_N,               ! #Strings
     $     SFRS_NVec(3,SFRS_Nmax), ! #(InVectors) of String i
     $     SFRS_Vec(4,3,SFRS_Nmax, SFRS_NVecmax),
     $     SFRS_VecID(3,SFRS_Nmax, SFRS_NVecmax),
     $     SFRS_DPS(5,SFRS_Nmax) ! Boost-Vector of String

      integer SFRS_N, SFRS_NVec, SFRS_VecID
      double precision SFRS_Vec
      double precision SFRS_DPS

      save /DataSFREP_S/

      integer j

      do j=1,5
         DPS(j) = SFRS_DPS(j,iS)
      enddo

      end

c===================================================================     
c Write out the given values.
c
      subroutine SFREPS_Write(iA)
      IMPLICIT NONE
      integer iA

      integer SFRS_Nmax, SFRS_NVecmax
      parameter (SFRS_Nmax=100,SFRS_NVecmax=100) ! former value: 10,50

      common /DataSFREP_S/
     $     SFRS_N,               ! #Strings
     $     SFRS_NVec(3,SFRS_Nmax), ! #(InVectors) of String i
     $     SFRS_Vec(4,3,SFRS_Nmax, SFRS_NVecmax),
     $     SFRS_VecID(3,SFRS_Nmax, SFRS_NVecmax),
     $     SFRS_DPS(5,SFRS_Nmax) ! Boost-Vector of String

      integer SFRS_N, SFRS_NVec, SFRS_VecID
      double precision SFRS_Vec
      double precision SFRS_DPS

      save /DataSFREP_S/

      double precision SFREPS_FOUR ! prototype 

      integer iS,iV,j, iV1
      double precision S(4)

      write(*,2001)

      if (iA.eq.0) then
c... write Boost-Vector:

         do iS=1,SFRS_N
            write(*,1011) iS,(SFRS_DPS(j,iS),j=1,5),
     $           ABS(SFRS_DPS(3,iS)).LT.0.99D0*SFRS_DPS(4,iS)
         enddo
         goto 999
      endif

c... write all other lists:

      do iS=1,SFRS_N
         do j=1,4
            S(j) = 0d0
         enddo
         do iV=1,SFRS_NVec(iA,iS)
            write(*,1000) iA,iS,iV,(SFRS_Vec(j,iA,iS,iV),j=1,4),
     $        SFRS_VecID(iA,iS,iV)   
            do j=1,4
               S(j) = S(j)+SFRS_Vec(j,iA,iS,iV)
            enddo
         enddo

         do j=1,3
            if (S(j).lt.1d-10) S(j) = 0d0
         enddo

         write(*,1001) (S(j),j=1,4)

c... write 4-Vector Products:

 990     continue
         
         if (iA.lt.3) then 
            write(*,*)
            write(*,1201) iA,(iS,iV1,SFRS_VecID(iA,iS,iV1),
     $           iV1=1,SFRS_NVec(iA,iS))
            do iV=1,SFRS_NVec(iA,iS)
               write(*,1210) iA,iS,iV,SFRS_VecID(iA,iS,iV),
     $              (SFREPS_FOUR(iA,iS,iV,iV1,1d-10,0d0),iV1=1,iV)
            enddo
c            write(*,2002)
         endif

c... end writing all other lines

         if (iS.ne.SFRS_N) write(*,2002)

      enddo

c... end of writing:

 999  write(*,2001)
      write(*,*)

 1000 FORMAT('<<SFREP_S,',i1,'>>{',i2,',',i2,'} :(',
     $     1P,3(e12.5,','),e12.5,') :[',i3,']')
 1001 FORMAT('                ==   :(',
     $     1P,3(e12.5,','),e12.5,')')

 1011 FORMAT('<<SFREP_S,B>>{',i2,'} :(',
     $     1P,4(e12.5,','),e12.5,') :',L1)

 1201 FORMAT('<<SFREP_S,',i1,
     $     '>> FourProd     ',
     $     20(' {',i1,',',i2,'}[',i3,']'))
 1210 FORMAT('<<SFREP_S,',i1,
     $     '>> {',i1,',',i2,'}[',i3,']',
     $     ' :',1P,20e12.4)

 2001 FORMAT(115('='))
 2002 FORMAT(115('-'))

      end

c===================================================================
c Find a vector V(1..4) in the sublist iA.
c
c INPUT:  iA      : Sublist to search
c         V(1..4) : 4-Vector to search
c OUTPUT: iString, iVec : matching String and matching vector of String
c                         (both =0 -> failure)  
c
c Failure: - vector not found
c
      subroutine SFREPS_FindVec(iA, V, iString,iVec)
      IMPLICIT NONE
      double precision V(4)
      integer iA, iString,iVec

      integer SFRS_Nmax, SFRS_NVecmax
      parameter (SFRS_Nmax=100,SFRS_NVecmax=100) ! former value: 10,50

      common /DataSFREP_S/
     $     SFRS_N,               ! #Strings
     $     SFRS_NVec(3,SFRS_Nmax), ! #(InVectors) of String i
     $     SFRS_Vec(4,3,SFRS_Nmax, SFRS_NVecmax),
     $     SFRS_VecID(3,SFRS_Nmax, SFRS_NVecmax),
     $     SFRS_DPS(5,SFRS_Nmax) ! Boost-Vector of String

      integer SFRS_N, SFRS_NVec, SFRS_VecID
      double precision SFRS_Vec
      double precision SFRS_DPS

      save /DataSFREP_S/

      integer iS, iV,nV,j

      do iS=1,SFRS_N
         nV = SFRS_NVec(iA,iS)
         do iV=1,nV
            do j=1,4
               if (V(j).ne.SFRS_Vec(j,iA,iS,iV)) goto 100
            enddo
            iString = iS
            iVec = iV
            return
 100     continue
         enddo
      enddo

      iString = 0
      iVec = 0
      end

c===================================================================
c Find a vector SFR_Vec(:,iL,iV) in the sublist iA.
c
c As SFREPS_FindVec(), but now the 4-Vector is taken from the 
c array SFR_Vec(..,..).
c Attention: now also the (internal) ID is compared!
c
c INPUT:  iA      : Sublist to search
c         iL, iV  : indizes of the vector to search
c OUTPUT: iString, iVec : matching String and matching vector of String
c                         (both =0 -> failure) 
c Failure: - no matching ID found
c          - Vector components do not match (no vector found)
c
      subroutine SFREPS_FindVecI(iA, iL,iV, iString,iVec)
      IMPLICIT NONE
      integer iA, iL,iV, iString,iVec

      integer SFRS_Nmax, SFRS_NVecmax
      parameter (SFRS_Nmax=100,SFRS_NVecmax=100) ! former value: 10,50

      common /DataSFREP_S/
     $     SFRS_N,               ! #Strings
     $     SFRS_NVec(3,SFRS_Nmax), ! #(InVectors) of String i
     $     SFRS_Vec(4,3,SFRS_Nmax, SFRS_NVecmax),
     $     SFRS_VecID(3,SFRS_Nmax, SFRS_NVecmax),
     $     SFRS_DPS(5,SFRS_Nmax) ! Boost-Vector of String

      integer SFRS_N, SFRS_NVec, SFRS_VecID
      double precision SFRS_Vec
      double precision SFRS_DPS

      save /DataSFREP_S/

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
      

      integer iS, iVs,nV,j

      do iS=1,SFRS_N
         nV = SFRS_NVec(iA,iS)
         do iVs=1,nV
            if (SFR_VecID(iL,iV).ne.SFRS_VecID(iA,iS,iVs)) goto 100
            do j=1,4
               if (SFR_Vec(j,iL,iV).ne.SFRS_Vec(j,iA,iS,iVs)) goto 100
            enddo
            iString = iS
            iVec = iVs
            return
 100     continue
         enddo
      enddo

      iString = 0
      iVec = 0
      end
c===================================================================
c Return a vector SFRS_Vec(:) from the sublist iA.
c
c INPUT:  iA      : Sublist to search
c         iString, iVec : indizes of vector
c OUTPUT: V(1..4) : 4-Vector
c         ID      : internal ID (<0 -> failure)
c
c Failure: - iVec larger than stored list 
c
      subroutine SFREPS_GetVec(iA,iString,iVec, V,ID)
      IMPLICIT NONE
      integer iA,iString,iVec, ID
      double precision V(4)

      integer SFRS_Nmax, SFRS_NVecmax
      parameter (SFRS_Nmax=100,SFRS_NVecmax=100) ! former value: 10,50

      common /DataSFREP_S/
     $     SFRS_N,               ! #Strings
     $     SFRS_NVec(3,SFRS_Nmax), ! #(InVectors) of String i
     $     SFRS_Vec(4,3,SFRS_Nmax, SFRS_NVecmax),
     $     SFRS_VecID(3,SFRS_Nmax, SFRS_NVecmax),
     $     SFRS_DPS(5,SFRS_Nmax) ! Boost-Vector of String

      integer SFRS_N, SFRS_NVec, SFRS_VecID
      double precision SFRS_Vec
      double precision SFRS_DPS

      save /DataSFREP_S/

      integer j

      if (iVec.gt.SFRS_NVec(iA,iString)) then
         ID = -1
         return
      endif

      do j=1,4
         V(j) = SFRS_Vec(j,iA,iString,iVec)
      enddo
      ID = SFRS_VecID(iA,iString,iVec)

      end
c===================================================================
c Return the vector SFRS_Vec(:) in the sublist iA with ID.
c
c INPUT:  iA      : Sublist to search
c         iString : string to search
c         ID      : (internal) ID of vector
c OUTPUT: iVec    : matching vector index (failure -> 0)
c         V(1..4) : matching vector components (failure -> 0)
c
c Failure: - ID not found
c 
      subroutine SFREPS_GetVecID(iA,iString,iVec, V, ID)
      IMPLICIT NONE
      integer iA,iString, iVec, ID
      double precision V(4)

      integer SFRS_Nmax, SFRS_NVecmax
      parameter (SFRS_Nmax=100,SFRS_NVecmax=100) ! former value: 10,50

      common /DataSFREP_S/
     $     SFRS_N,               ! #Strings
     $     SFRS_NVec(3,SFRS_Nmax), ! #(InVectors) of String i
     $     SFRS_Vec(4,3,SFRS_Nmax, SFRS_NVecmax),
     $     SFRS_VecID(3,SFRS_Nmax, SFRS_NVecmax),
     $     SFRS_DPS(5,SFRS_Nmax) ! Boost-Vector of String

      integer SFRS_N, SFRS_NVec, SFRS_VecID
      double precision SFRS_Vec
      double precision SFRS_DPS

      save /DataSFREP_S/

      integer j

      do j=1,SFRS_NVec(iA,iString)
         if (SFRS_VecID(iA,iString,j).eq.ID) then
            iVec = j
            goto 100
         endif
      enddo

      call SFREP_ERROR(-4,'SFREPS_GetVecID: ID not found.')
      write(*,*) 'iA, iString, ID=',iA,iString,ID

      do j=1,4
         V(j) = 0d0
      enddo
      iVec = 0
      return

 100  continue

      do j=1,4
         V(j) = SFRS_Vec(j,iA,iString,iVec)
      enddo

      end

c===================================================================
c Write an Error Message and set a error flag

      subroutine SFREP_ERROR(iErr,S)
      IMPLICIT NONE
      integer iErr
      character S*(*)

      common /DataSFREP_Err/ SFR_Err
      integer SFR_Err
      save /DataSFREP_Err/
      
      SFR_Err = iErr

      write(*,*) 'SFREP_ERROR:',S,' [ErrCode:',iErr,']'
      
      end
c===================================================================
c Check the vectors saved

      subroutine SFREPS_Check
      IMPLICIT NONE

      integer SFRS_Nmax, SFRS_NVecmax
      parameter (SFRS_Nmax=100,SFRS_NVecmax=100) ! former value: 10,50

      common /DataSFREP_S/
     $     SFRS_N,               ! #Strings
     $     SFRS_NVec(3,SFRS_Nmax), ! #(InVectors) of String i
     $     SFRS_Vec(4,3,SFRS_Nmax, SFRS_NVecmax),
     $     SFRS_VecID(3,SFRS_Nmax, SFRS_NVecmax),
     $     SFRS_DPS(5,SFRS_Nmax) ! Boost-Vector of String

      integer SFRS_N, SFRS_NVec, SFRS_VecID
      double precision SFRS_Vec
      double precision SFRS_DPS

      save /DataSFREP_S/

      integer iS

      do iS=1,SFRS_N
         if (SFRS_NVec(3,iS).eq.0) then
c$$$            write(*,*) 'SFREPS_Check: Empty String ',iS
c$$$
c$$$            call PYLIST(2)
c$$$            call SFREPS_Write(1)
c$$$            call SFREPS_Write(2)
c$$$            call SFREPS_Write(3)
c$$$            call SFREP_Write(1,0)
c$$$c            stop

            SFRS_NVec(1,iS) = 0
            SFRS_NVec(2,iS) = 0
         endif
      enddo
      end

c===================================================================
c Four Vector Product of two stored vectors in list (iA,iS)
c
c There is an additional check for the return value: (C-syntax)
c    Value = (Value<DefB ? DefV : Value)

      double precision function SFREPS_FOUR(iA,iS,iV1,iV2,DefB,DefV)
      IMPLICIT NONE
      integer iA,iS,iV1,iV2
      double precision DefB, DefV
     

      integer SFRS_Nmax, SFRS_NVecmax
      parameter (SFRS_Nmax=100,SFRS_NVecmax=100) ! former value: 10,50

      common /DataSFREP_S/
     $     SFRS_N,               ! #Strings
     $     SFRS_NVec(3,SFRS_Nmax), ! #(InVectors) of String i
     $     SFRS_Vec(4,3,SFRS_Nmax, SFRS_NVecmax),
     $     SFRS_VecID(3,SFRS_Nmax, SFRS_NVecmax),
     $     SFRS_DPS(5,SFRS_Nmax) ! Boost-Vector of String

      integer SFRS_N, SFRS_NVec, SFRS_VecID
      double precision SFRS_Vec
      double precision SFRS_DPS

      save /DataSFREP_S/

      SFREPS_FOUR = SFRS_Vec(4,iA,iS,iV1)*SFRS_Vec(4,iA,iS,iV2)
     $     -SFRS_Vec(1,iA,iS,iV1)*SFRS_Vec(1,iA,iS,iV2)
     $     -SFRS_Vec(2,iA,iS,iV1)*SFRS_Vec(2,iA,iS,iV2)
     $     -SFRS_Vec(3,iA,iS,iV1)*SFRS_Vec(3,iA,iS,iV2)

      if (SFREPS_FOUR.lt.DefB) SFREPS_FOUR=DefV

      return
      end
c===================================================================

