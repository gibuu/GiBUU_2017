c*****************************************************************==
c****m* /GetJetsetVec
c NAME
c GetJetsetVec.F
c PURPOSE
c return a list of 'production--' and 'formation--vertices'
c of the hadrons for a given PYTHIA v6.2 event.
c
c HISTORY
c Kai Gallmeister, 09.11.2004...30.01.2006
c
c This relies on (and replaces):
c * GetFormVec.F:  Kai Gallmeister, 24.8.2004...30.9.2004
c * GetFormTime.F: Kai Gallmeister, 23.3.2004...17.8.2004
c 
c some dates:
c * 31.01.05: Output via common block /DataGJV/, 
c   cluster->2: eArr(1,..)=3,4: last/next-to-last hadron
c * 02.03.05: cluster->2: Times correct, eArr(1,..)=0 again
c * 30.03.05: Pythia v6.225
c
c NOTES
c
c The Pythia-subroutine PYSTR has to be expanded to report the history
c of the String-Fragmentation through calls to routines in 
c "HandleSFREP.F", "HandleSFREP_S.F", which use internal common blocks.
c
c Link your code with "PYSTRF.Modif.6225.o" (or "PYSTRF.Modif.6225.o").
c
c Please note:
c A lot of the work is done in [HandleSFREP.F] via: 
c SFREP_ReconLastQQ, SFREP_ReconLastQgQ, (SFREP_ReconLastGG)
c
c The Input/Output is done via the common block /DataGJV/
c*****************************************************************==


c*****************************************************************==
c
c USER CALLABLE ROUTINES:
c =======================
c
c subroutine GetJetsetVecINIT()
c   ... initialize calculation of the vertices
c
c subroutine GetJetsetVec()
c   ... calculate the vertices from the Jetset-Fragmentation
c
c subroutine GetJetsetVecPYEDIT()
c   ... reorder Arr,EArr like a PYEDIT(1) call
c
c subroutine GetJetsetVecPYROBO(THE,PHI,BEX,BEY,BEZ)
c   ... boost like PYROBO would do
c
c subroutine GetJetsetVecCheckT(GammaLevel)
c   ... check times and reset EArr-Flags
c
c Internal Routines:
c ------------------
c
c - GJV_Cluster, GJV_StringQQ, GJV_StringQgQ, GJV_StringGG,
c   GJV_String_0
c
c - Low Level: 
c   AddVert_Rep, SetEArr, BoostVert, TransformLL
c   SetArrayD, SetArrayI, SetArrayArr
c   AddVert, AddVertSFR, VertLength
c
c   GJV_ERR, GJV_MSG, GJV_MSG_L, GJV_ListFinalEvent
c
c - SFREP_*, SFREPS_* : Handle /DataSFREP/, /DataSFREP_S/
c
c - Dummies for LULIST and PYLIST are included.
c
c-------------------------------------------------------------------
c
c GetJetsetVec(nArrMax, Arr,EArr,verb,AtOrigin):
c ============================================
c
c The array Arr gives as output the two production times and
c the formation time.
c The array EArr reports possible errors (see below) and the actual
c string number.
c nArrmax (as Array size) should be choosen 4000.
c If verb is set to a nonzero value, some reports are written. 
c The Flag AtOrigin changes the treatment of the outmost production
c points: if AtOrigin is TRUE, then the production point is the origin
c (in the cms of the string/cluster), otherwise it is the turning point.
c
c All Output is suitable to the Jetset Arrays BEFORE any PYEDIT(1) calls.
c The routine GetJetsetVecPYEDIT shrinks the arrays as the PYEDIT routine.
c
c You have to ensure, that before every call of PYEVNT
c     call SFREP_Reset(1,0)
c     call SFREPS_Reset(0)
c is given. Then, after
c     call PYEVNT
c and
c     call GetJetsetVec(...,Arr,EArr,...)
c the array Arr holds the production and formation 4-vectors 
c and the array EArr indicates possible errors and the string number.
c Now 
c     call GetJetsetVecPYEDIT(...,Arr,EArr)
c shrinks the output arrays, to be compatible with the momenta after
c     call PYEDIT(...)
c in the PYTHIA-Arrays.
c
c
c The coding of the flag EArr(...) is as follows:
c The number of the string is given in EArr(4,...)
c The Rank of the particle is given in EArr(5,...), like it was
c calculated, EArr(6,...) gives the rank as the minmal number of
c particles to the left or right (plus 1). (reverse ordering is considered)
c
c The meaning of the 3 lowest entries EArr(1,...) to EArr(3,...) is:
c
c   3 2 1
c  -------
c   x x 0 : hadron has not been processed
c   x x 1 : everything should be ok (in principle)
c   x x 3 : next to last hadron in fragmentation
c   x x 4 : last hadron in fragmentation
c
c   x 1 x : problems with ProdTime 1 ; additive
c   x 2 x : problems with ProdTime 2 ; additive
c   x 4 x : problems with FormTime   ; additive
c
c   1 x x : hadron from QQ-String
c   2 x x : hadron from QQ-String with internal Gluons
c   3 x x : hadron from pure gluonic GG-string
c   4 x x : hadron from Cluster-Decay -> 1 hadron
c   5 x x : hadron from Cluster-Decay -> 2 hadrons
c   6 x x : hadron from Doku-Line
c
c Values of EArr(2,...) can be tested by the intrinsic f77-function
c IAnd(...,1|2|4).
c 
c===================================================================
c (Known) Problems of the current version:
c
c - calculated Vectors in GG-Strings are more or less not understood
c
c - Reverse ordering is implemented in an unsatisfactory way: only if
c   StringBreak ordering is unordered (i.e. only for special cases in
c   QgQ-Strings), reverse ordering is taken into account.
c
c - The formation time in QgQ-Strings needs some reaximination.
c   (Occurence of tau_F<tau_P! Probably due to some xi<0.)
c
c - Addition of vectors in QgQ-Strings should be more sophisticated,
c   i.e. the use of /DataRepString0/ and AddVert_Rep should be
c   eliminated.
c
c - A lot of source code is just for reporting or some old stuff.
c   The code could be much more stringent!
c
c===================================================================
c! BE AWARE:
c! =========
c! The Code returns values for the prod- and form-vectors. There is no
c! guarantee, that these values make sense! In detail:
c!
c! - the errorflag reports errors that are catched.
c!   There may be errors (eg. NAN) that are missed
c!   Also, form-times may be smaller than prod-times.
c!
c! - Ordering along the string:
c!   Ordering of the hadrons in the same direction as the string building
c!   quarks is assumed, but this may be reversed. 
c!   Inverse ordering could cure the problem; 
c!   but "Why only there?", should one implement the possibility of
c!   inverse ordering everywhere?
c!
c! - Pure gluonic strings (GG) are handled at the moment in a somehow
c!   crude way. Handling like QgQ-strings is anticipated, but first tries
c!   showed lot of error messages.
c! 
c! - Jetset relies on m_q==0 in QgQ-Strings: formulae assume, that
c!   V(i)==V(i+3), i.e. that all p_+/- from a gluon are identical. 
c!   Setting the Q- and Qbar-momenta to lightlike (cf. TransformLL) 
c!   changes 'gluon'-p_+/-, which is not taken care of in the code.   
c!
c! - Reported Gamma values in QQ-Strings (and for QQ-like string parts in
c!   QgQ strings) differ dramatically from calculated values. This is due
c!   to the fact that here Jetset uses mT to calculate Gamma, while
c!   the presented "reconstructing" procedures neglect pT and therefore
c!   use only the invariant mass m. 
c!
c===================================================================

c*****************************************************************==
c Initialize calculation of the vertices
c

      subroutine GetJetsetVecINIT()
      IMPLICIT NONE

      integer nArrMax
      parameter (nArrMax=200)   ! maximum size of arrays

      common /DataGJV/
     $     Arr(3,4,nArrMax),    ! 3* 4D-Vertizes
     $     EArr(6,nArrMax),     ! errFlag, rank
     $     verb,                ! verbosity
     $     AtOrigin             ! treatment of outmost prod points

      double precision Arr
      integer EArr
      integer verb
      logical AtOrigin

      save /DataGJV/

c... Prepare for reporting:

      call SFREP_Reset(1,0)
      call SFREPS_Reset(0)

c... Reset the output arrays:

      Arr = 0d0
      EArr = 0

      end

c*****************************************************************==
c****s* GetJetsetVec/GetJetsetVec
c NAME
c subroutine GetJetsetVec(withPythia)
c PURPOSE
c Calculate the vertices from the Jetset-Fragmentation;
c Return Arrays of Production and Formation Vertices of an Event
c
c This is the (only) routine to be called by the user.
c
c Arr is an array of dimension (3,4,nArrMax), which holds the output:
c   1,2: production time 1 and 2, 3: formation time
c
c EArr is an array of dimension (4,nArrMax), which holds error flags
c and the string number the particle belongs to
c
c If verb =1,2,... some output is produced.
c
c*****************************************************************==
      subroutine GetJetsetVec(withPythia)
      IMPLICIT NONE

      logical withPythia
      
c...common blocks:

      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      integer N,NPAD,K
      double precision P,V
      SAVE /PYJETS/

      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      integer MSTP,MSTI
      double precision PARP,PARI
      SAVE /PYPARS/


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
      

      common /DataSFREP_Err/ SFR_Err
      integer SFR_Err
      save /DataSFREP_Err/

      integer nArrMax
      parameter (nArrMax=200)   ! maximum size of arrays

      common /DataGJV/
     $     Arr(3,4,nArrMax),    ! 3* 4D-Vertizes
     $     EArr(6,nArrMax),     ! errFlag, rank
     $     verb,                ! verbosity
     $     AtOrigin             ! treatment of outmost prod points

      double precision Arr
      integer EArr
      integer verb
      logical AtOrigin

      save /DataGJV/

      common /DataGJVerror/ UseWithPYTHIA, ListFinalEvent
      integer UseWithPYTHIA, ListFinalEvent
      save /DataGJVerror/

c      data UseWithPYTHIA /1/

      integer iL                ! actual line
      integer sL                ! line where string is listed
      integer pL1, pL2          ! parent lines (see below)
      integer cL1, cL2          ! child lines

      integer pL1C              ! saves K(pL1,1) = K(pL1,1)

      integer iString           ! nr. of current string

      integer MSTI4

c...helpful variables:

      logical HasDokuLine

      integer j

      HasDokuLine = .FALSE.

      if (withPythia) then
         UseWithPYTHIA = 1
         MSTI4 = MSTI(4)
      else
         UseWithPythia = 0
         MSTI4 = 0
      endif


      SFR_Err = 0               ! all ok
      call SFREP_SetN           ! SFR_Ntot = N

      call SFREPS_Check         ! Check the vectors saved

c      iL = MSTI(4)
      iL=0                      ! WARUM DAS HIER ?????
      iString = 0

 100  iL = iL+1
      if (K(iL,1).ge.10) goto 999 ! not stable, skip ahead

      pL1 = K(iL,3)
      if (pL1.le.0) then        ! check this for FRITIOF
         pL1C = 21
      else
         pL1C = K(pL1,1)
      endif

      if ((pL1C.eq.12).or.(pL1C.eq.51)) then 

c...beginning of q-string:
c...======================
         
c     ...pL1,pL2 = begin & end of quark string (parent)
c     ...cL1,cL2 = begin & end of resulting particles (child)

         cL1 = K(pL1,4)
         cL2 = K(pL1,5)

         if (cL1.ne.iL) call GJV_ERR('cL1.ne.iL',0,iL)

         do j=pL1+1,N
            if (K(j,1).eq.11) then
               pL2 = j
               goto 111
            endif
         enddo
         call GJV_ERR('end of q-string not found',0,iL)
 111     continue

c... search line with 'string':
         do j=pL2+1,N
            if ((K(j,4).eq.cL1).and.(K(j,5).eq.cL2)) then
               sL = j
               goto 112
            endif
         enddo
         call GJV_ERR('string line not found',0,iL)
 112     continue

         iString = iString+1

c... set first line to be used in output array
c... and check output array size: [ED: do not delete !!!]
         
         if (nArrMax<cL2) then
            write(*,*) 'GetJetsetVec: FATAL ERROR: nArrMax < ', cL2
            call PYLIST(2)
            stop
         endif

c... calculate values:

         if (K(sL,2).eq.91) then ! This was a cluster decay !
            call GJV_Cluster(pL1,pL2,sL,cL1,cL2,iString)
         else                   ! This was a string decay !
            if (K(pL1,2).eq.21) then
               call GJV_StringGG(pL1,pL2,sL,cL1,cL2,iString)
            else
               if (pL2-pL1.eq.1) then
                  call GJV_StringQQ(pL1,pL2,sL,cL1,cL2,iString)
               else
                  call GJV_StringQgQ(pL1,pL2,sL,cL1,cL2,iString)
               endif
            endif
         endif

      else if (pL1C.eq.21) then 

c...from decayed particle:
c...======================

         HasDokuLine = .TRUE.
         
         if (pL1.le.MSTI4) then ! from docu
            call SetEArr(iL, 1,1, 0, 6, 0, 1)
         else
            call GJV_ERR('where does this come from?',9,iL)
         endif
         cL2 = iL         ! last line worked on
         
      else
         
c...  something unknown:
c...==================

         call SFREP_ERROR(-1,'SEVERE PROBLEM: unknown Code pL1c')
         write(*,*) 'pL1c = ',pL1C,' iL = ',iL
         return
         
c         call GJV_ERR('unknown code',0,iL)

      endif
      
      iL = cL2

 999  if (iL.lt.N) goto 100  ! loop: find next string 

      if (HasDokuLine) call LookForDokuLine

      end

c-----------------------------------------------------------------

      subroutine LookForDokuLine
      IMPLICIT NONE

c...common blocks:

      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      integer N,NPAD,K
      double precision P,V
      SAVE /PYJETS/

      integer nArrMax
      parameter (nArrMax=200)   ! maximum size of arrays

      common /DataGJV/
     $     Arr(3,4,nArrMax),    ! 3* 4D-Vertizes
     $     EArr(6,nArrMax),     ! errFlag, rank
     $     verb,                ! verbosity
     $     AtOrigin             ! treatment of outmost prod points

      double precision Arr
      integer EArr
      integer verb
      logical AtOrigin

      save /DataGJV/


      integer iL1,iL2, iPL

      logical IdenticalLines ! PROTOTYPE

c$$$      write(*,*)    
c$$$      write(*,*)
c$$$      write(*,*) '#############################################'
c$$$      call PYGIVE('MSTP(14)=')
c$$$      call PYGIVE('MSTI(1)=')
c$$$      call PYGIVE('MINT(1)=')
c$$$      call PYGIVE('MSTI(9)=')
c$$$      call PYGIVE('MINT(103)=')
c$$$      call PYGIVE('MINT(107)=')
c$$$      call PYGIVE('MINT(122)=')
c$$$      call PYGIVE('MINT(123)=')
c$$$      call PYGIVE('MINT(17)=')
c$$$      call PYGIVE('MINT(18)=')
c$$$      call PYLIST(2)
c$$$!      call GetJetSetVec_List(6,1,N)

      iL1 = 0
 100  iL1 = iL1+1
      if (iL1>=N) goto 900
      if (iL1>=nArrMax) goto 900
      if (EArr(3,iL1).ne.6) goto 100
      if (K(iL1,2).eq.22) goto 100 ! skip photons

      iPL = K(iL1,3)
      iL2 = iL1

 200  iL2=iL2+1
      if ((EArr(3,iL2)==6).and.(K(iL2,3)==iPL)) goto 200

      iL2 = iL2-1

! we found something:


      if (iL2-iL1 == 0) then
         if (IdenticalLines(iL1,iPL)) then
c... This is scattered "elastically". Therefore: tP1 = tP2 = tF = 0
c... (nothing to do)

         else
c... This is like Cluster->1
            Arr(1,1:4,iL1) = 0d0
            Arr(2,1:4,iL1) = 0d0
            Arr(3,1:4,iL1) = P(iL1,1:4)
            EArr(2,iL1) = 0     ! all ok
         endif

      else if (iL2-iL1 == 1) then
!         write(*,*) 'found:', iL1,iL2,' -> ', iPL,' ### 2'

         call Guess2BodyTimes(iPL,iL1,iL2)
!         call GetJetSetVec_List(6,1,N)
!         stop
      else
         write(*,*) 'found:', iL1,iL2,' -> ', iPL,' ### more'
         call PYLIST(2)
         stop
      endif




      iL1 = iL2
      goto 100


 900  continue

!      stop

      end

c*****************************************************************==

      logical function IdenticalLines(iL1,iL2)
      IMPLICIT NONE
      integer iL1, iL2

      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      integer N,NPAD,K
      double precision P,V
      SAVE /PYJETS/

      integer i

      IdenticalLines = .FALSE.

      if ((iL1==0).or.(iL2==0)) return

      if (K(iL1,2) .ne. K(iL2,2)) return
!      do i=1,5
      do i=5,5
         if (P(iL1,i) .ne. P(iL2,i)) return
      enddo

      IdenticalLines = .TRUE.

      end

c*****************************************************************==
c****s* GetJetsetVec/Guess2BodyTimes
c NAME
c subroutine Guess2BodyTimes(pL,cL1,cL2)
c PURPOSE
c 
c*****************************************************************==
      subroutine Guess2BodyTimes(pL,cL1,cL2)
      IMPLICIT NONE
      integer pL,cL1,cL2

c... Common Blocks

      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      integer N,NPAD,K
      double precision P,V
      SAVE /PYJETS/

      integer nArrMax
      parameter (nArrMax=200)   ! maximum size of arrays

      common /DataGJV/
     $     Arr(3,4,nArrMax),    ! 3* 4D-Vertizes
     $     EArr(6,nArrMax),     ! errFlag, rank
     $     verb,                ! verbosity
     $     AtOrigin             ! treatment of outmost prod points

      double precision Arr
      integer EArr
      integer verb
      logical AtOrigin

      save /DataGJV/

c... internal variables/functions

      double precision PLL(2,4)

      double precision FLL12,xiP(2),xiM(2)
      double precision sqrts,m1,m2

      double precision sqrtsR, thetaR, phiR, betaR(3)

      integer j

      double precision MP_ScalProd4 ! prototype
      double precision MP_P ! prototype

      call MP_Set4(1, P(cL1,5), P(cL1,1),P(cL1,2),P(cL1,3),P(cL1,4))
      call MP_Set4(2, P(cL2,5), P(cL2,1),P(cL2,2),P(cL2,3),P(cL2,4))
!      call MP_Set4(3, 0d0, 0d0,0d0, sqrts/2,sqrts/2)
!      call MP_Set4(4, 0d0, 0d0,0d0,-sqrts/2,sqrts/2)
      
      call MP_CalcROBO(1,2, sqrtsR, thetaR, phiR, betaR)

      sqrts = sqrtsR
      FLL12 = sqrts**2/2
      m1 = P(cL1,5)
      m2 = P(cL2,5)

      if (sqrts < m1+m2) then 
         write(*,*) 'ooops: sqrts < m1+m2 ::', sqrts,m1,m2
         write(*,*) '       cL1,cL2,pL = ',cL1,cL2,pL
         call PYLIST(2)
      endif

!      p0 = sqrt( (sqrts**2-m1**2-m2**2)**2-4*m1**2*m2**2 )/(2*sqrts)



      call MP_Set4(3, 0d0, 0d0,0d0, sqrts/2,sqrts/2)
      call MP_Set4(4, 0d0, 0d0,0d0,-sqrts/2,sqrts/2)

!      call MP_Write(6,1,4)

!      call MP_ROBO_INV(1,2,thetaR,phiR,betaR(1),betaR(2),betaR(3))
      call MP_ROBO(3,4,thetaR,phiR,betaR(1),betaR(2),betaR(3))


!      write(*,*) '### CL=',cL1,cL2, pL
!      call MP_Write(6,1,4)
!      call PYLIST(2)
!      call GetJetSetVec_List(6,0,0)

      do j=1,4
         PLL(1,j) = MP_P(3,j) ! <- p_+
         PLL(2,j) = MP_P(4,j) ! <- p_-
      enddo


!      write(*,*) sqrt(2*MP_ScalProd4(3,4))
!      write(*,*) MP_SqrtS(1,2),MP_SqrtS(3,4)


      xiP(1) = MP_ScalProd4(1,4)/FLL12
      xiP(2) = MP_ScalProd4(2,4)/FLL12
      xiM(1) = MP_ScalProd4(1,3)/FLL12
      xiM(2) = MP_ScalProd4(2,3)/FLL12
      

!      write(*,*) 'P:',xiP(1),xiP(2)
!      write(*,*) 'M:',xiM(1),xiM(2)

      if (xiP(2)*xiM(1) < xiP(1)*xiM(2)) then
!         write(*,*) 'links'

         Arr(1,1:4,cL1) = PLL(1,1:4)
         Arr(2,1:4,cL1) = PLL(1,1:4)*xiP(2) + PLL(2,1:4)*xiM(1)
         Arr(3,1:4,cL1) = PLL(1,1:4)        + PLL(2,1:4)*xiM(1)

         Arr(1,1:4,cL2) = Arr(2,1:4,cL1)
         Arr(2,1:4,cL2) =                     PLL(2,1:4)
         Arr(3,1:4,cL2) = PLL(1,1:4)*xiP(2) + PLL(2,1:4)

      else
         write(*,*) 'Guess2BodyTimes: rechts! STOP'
         call PYLIST(2)
         stop
      endif

      if (AtOrigin) then
         Arr(1,1:4,cL1) = 0d0
         Arr(2,1:4,cL2) = 0d0
      endif 

      EArr(2,cL1) = 0           ! all ok
      EArr(2,cL2) = 0           ! all ok

!      call GetJetSetVec_List(6,0,0)

!      if (pl <= 0) stop

!      stop


      end
c*****************************************************************==
c Handle a cluster decay
c
c INPUT:
c pL1,pL2,sL,cL1,cL2: lines in the Event record
c iStrNr:             actual String number
c
c OUTPUT:
c Arr, EArr:          some entries changed
c
c-----------------------------------------------------------------
c cluster -> 1 hadron:
c --------------------
c Produktion-Punkte werden auf {0} gesetzt. (AtOrigin hat keine Wirkung)
c Fomation-Punkt: entspricht Hadron-Masse (nicht Cluster-Masse)
c ... Die Cluster-Masse wird in Jetset auf eine Hadron-Masse ge-map-t,
c     wobei ein Gluon in den naechsten String/Cluster eingebaut wird, um
c     dies zu kompensieren.
c
c cluster -> 2 hadrons:
c ---------------------
c Hier gibt es das beruehmt/beruechtigte "Inverse Ordering" Problem. 
c (vgl. Year04/W25), Year05/W10)
c
c Als erstes projeziere ich die Cluster-building Impulse auf masselose
c partonen: call TransformLL
c Trotzdem gibt es zwei Zuordnungsmoeglichkeiten...
c Ich habe im Folgenden eine bestimmte Loesung gewaehlt. Die alternative
c Loesung ('reverse ordering') vertauscht die Reihenfolge der Quark-Impulse.
c Gemaess 'area law' sind beide Loesungen moeglich, waehrend die eine
c statistisch unterdrueckt ist. 
c Dies ist allerdings schon im JetSet-Output beruecksichtigt:
c Berechne ich \Delta = \Gamma^{rechts)-\Gamma^{links}, so ist dieser
c Wert negativ, wenn Jetset "inverse ordering" angewandt hat.
c Somit kann ich durch Benutzung einer Berechnungs-Variante Zeiten
c ausrechnen, die manchmal der groesseren, manchal der kleineren
c Vertex-Zeit entsprechen. Die Wahl hat jedoch schon Jetset getroffen.
c In ca. 60% der Faelle ist Delta<0, d.h. die 'rechte' Kombination
c liefert kleinere Produktionszeiten als die 'linke'. Da die kleinere
c Produktionszeit  gemaess AreaLaw wahrscheinlicher sein soll, ist also
c die RECHTE Kombination (fuer alle Berechnungen) die richtige, auch wenn
c sie manchmal groessere Produktionszeiten liefert!
c
c (Leider gab es in alten Code-Varianten einen Vertauschung/Fehler bei
c der Berechnung von xiP(j),xiM(j) !!!)
c
c-----------------------------------------------------------------

      subroutine GJV_Cluster(pL1,pL2,sL,cL1,cL2, iStrNr)
      IMPLICIT NONE

      integer pL1,pL2,sL,cL1,cL2
      integer iStrNr

      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      integer N,NPAD,K
      double precision P,V
      SAVE /PYJETS/

      integer nArrMax
      parameter (nArrMax=200)   ! maximum size of arrays

      common /DataGJV/
     $     Arr(3,4,nArrMax),    ! 3* 4D-Vertizes
     $     EArr(6,nArrMax),     ! errFlag, rank
     $     verb,                ! verbosity
     $     AtOrigin             ! treatment of outmost prod points

      double precision Arr
      integer EArr
      integer verb
      logical AtOrigin

      save /DataGJV/

      integer cL,iCode
      integer cL12, j

      double precision PLL(2,4)

      integer iF1,iF2
      double precision FOURL, FOURLL
      FOURL(iF1,iF2) = P(iF1,4)*PLL(iF2,4)
     $     -P(iF1,1)*PLL(iF2,1)
     $     -P(iF1,2)*PLL(iF2,2)
     $     -P(iF1,3)*PLL(iF2,3)
      FOURLL(iF1,IF2) = PLL(iF1,4)*PLL(iF2,4)
     $     -PLL(iF1,1)*PLL(iF2,1)
     $     -PLL(iF1,2)*PLL(iF2,2)
     $     -PLL(iF1,3)*PLL(iF2,3)

      double precision FLL12,xiP(2),xiM(2)


      integer cDelta(2)
      data cDelta /0,0/
      save cDelta

      cL12 = cL2-cL1+1
      if (cL12.eq.1) then
         iCode = 4
      else if (cL12.eq.2) then
         iCode = 5
      else
         write(*,*) 'GJV_Cluster: FATAL ERROR: cL12 == ',
     $        cL12,' STOP!!!'
         stop
      endif

c... set default values:

      do cL=cL1,cL2
         call SetEArr(cL, 1,1, iStrNr,iCode,7,1)
      enddo

c... Cluster -> 1 hadron

      if (iCode.eq.4) then
         Arr(1,1:4,cL1) = 0d0
         Arr(2,1:4,cL1) = 0d0
         Arr(3,1:4,cL1) = P(cL1,1:4)

         EArr(2,cL1) = 0        ! all ok
         return                 ! -> ENDE
      endif

c... Cluster build up with Gluons:

      if (pL2-pL1.gt.1) then    ! this is a Cluster with GLUON
         return                 ! error==7 is set !
      endif

c... Cluster -> 2 hadrons
c... SEE EXPLANATION TO THIS ROUTINE FOR THE FOLLOWING!!!!!
c...
c... p_- = PLL(1), p_+ = PLL(2), p^G = P(cL1), p^R = P(cL2)
c... G == 1 (cL1) , R == 2 (cL2) 

      call TransformLL(pL1,pL2,PLL)

      FLL12 = FOURLL(1,2)
      xiP(1) = FOURL(cL1,1)/FLL12 ! (p^G p_-) / (p_+ p_-), G
      xiP(2) = FOURL(cL2,1)/FLL12 ! (p^R p_-) / (p_+ p_-), R
      xiM(1) = FOURL(cL1,2)/FLL12 ! (p^G p_+) / (p_+ p_-), G
      xiM(2) = FOURL(cL2,2)/FLL12 ! (p^R p_+) / (p_+ p_-), R

      do j=1,4
         Arr(1,j,cL1) =                          PLL(1,j)
         Arr(2,j,cL1) = xiP(1)*PLL(2,j) + xiM(2)*PLL(1,j)
         Arr(3,j,cL1) = xiP(1)*PLL(2,j) +        PLL(1,j)

         Arr(1,j,cL2) = Arr(2,j,cL1) ! gleicher Vertex wie oben
         Arr(2,j,cL2) =        PLL(2,j)
         Arr(3,j,cL2) =        PLL(2,j) + xiM(2)*PLL(1,j)
      enddo 
    
      EArr(2,cL1) = 0
      EArr(2,cL2) = 0
      
      if (AtOrigin) then
         do j=1,4
            Arr(1,j,cL1) = 0d0
            Arr(2,j,cL2) = 0d0
         enddo
      endif 

      return

      end

c*****************************************************************==
c Handle a qq-string decay (i.e. a (1+1)-String)
c
c INPUT, OUTPUT: as GJV_Cluster
c
      subroutine GJV_StringQQ(pL1,pL2,sL,cL1,cL2, iStrNr)
      IMPLICIT NONE

      integer pL1,pL2,sL,cL1,cL2
      integer iStrNr

      common /DataSFREP_Err/ SFR_Err
      integer SFR_Err
      save /DataSFREP_Err/

      integer nArrMax
      parameter (nArrMax=200)   ! maximum size of arrays

      common /DataGJV/
     $     Arr(3,4,nArrMax),    ! 3* 4D-Vertizes
     $     EArr(6,nArrMax),     ! errFlag, rank
     $     verb,                ! verbosity
     $     AtOrigin             ! treatment of outmost prod points

      double precision Arr
      integer EArr
      integer verb
      logical AtOrigin

      save /DataGJV/

      integer cL,iCode
      parameter (iCode=1)

c... set default values:

      do cL=cL1,cL2
         call SetEArr(cL, 0,0, iStrNr,iCode,7,1)
      enddo

c... Reconstruct last and next to last vector

      call SFREP_ReconLastQQ(cL1,cL2)
      if (SFR_Err.ne.0) return

      call GJV_String_0(cL1,cL2, iCode)

      end

c*****************************************************************==
c Handle a q...g...q-string decay
c
c INPUT, OUTPUT: as GJV_Cluster
c
      subroutine GJV_StringQgQ(pL1,pL2,sL,cL1,cL2, iStrNr)
      IMPLICIT NONE

      integer pL1,pL2,sL,cL1,cL2
      integer iStrNr

      common /DataSFREP_Err/ SFR_Err
      integer SFR_Err
      save /DataSFREP_Err/

      integer nArrMax
      parameter (nArrMax=200)   ! maximum size of arrays

      common /DataGJV/
     $     Arr(3,4,nArrMax),    ! 3* 4D-Vertizes
     $     EArr(6,nArrMax),     ! errFlag, rank
     $     verb,                ! verbosity
     $     AtOrigin             ! treatment of outmost prod points

      double precision Arr
      integer EArr
      integer verb
      logical AtOrigin

      save /DataGJV/

      integer cL,iCode
      parameter (iCode=2)

c... set default values:

      do cL=cL1,cL2
         call SetEArr(cL, 0,0,iStrNr,iCode,7,1)
      enddo

c... Reconstruct last and next to last vector

      call SFREP_ReconLastQgQ(cL1,cL2)
      if (SFR_Err.ne.0) return

      call GJV_String_0(cL1,cL2, iCode)

      end
c*****************************************************************==
c Handle a g...g...g-string decay (i.e. a pure gluonic String)
c
c INPUT, OUTPUT: as GJV_Cluster
c
      subroutine GJV_StringGG(pL1,pL2,sL,cL1,cL2, iStrNr)
      IMPLICIT NONE

      integer pL1,pL2,sL,cL1,cL2
      integer iStrNr

      common /DataSFREP_Err/ SFR_Err
      integer SFR_Err
      save /DataSFREP_Err/

      integer nArrMax
      parameter (nArrMax=200)   ! maximum size of arrays

      common /DataGJV/
     $     Arr(3,4,nArrMax),    ! 3* 4D-Vertizes
     $     EArr(6,nArrMax),     ! errFlag, rank
     $     verb,                ! verbosity
     $     AtOrigin             ! treatment of outmost prod points

      double precision Arr
      integer EArr
      integer verb
      logical AtOrigin

      save /DataGJV/

      integer cL,iCode
      parameter (iCode=3)

c... set default values:

      do cL=cL1,cL2
         call SetEArr(cL, 0,0, iStrNr,iCode,7,1)
      enddo

c... Reconstruct last and next to last vector

      call SFREP_ReconLastGG(cL1,cL2)
      if (SFR_Err.ne.0) return

      call GJV_String_0(cL1,cL2, iCode)

      end

c*****************************************************************==
c Handle a String decay after reconstructiong all vectors
c
c (cf. SFREP_Write)
c
c-------------------------------------------------------------------
c Dies ist eine "brute force"/"dirty hack" routine !!!
c Bessere Loesungen sind willkommen!!!
c-------------------------------------------------------------------
c
c

      subroutine GJV_String_0(cL1,cL2, iCode)
      IMPLICIT NONE
      integer cL1,cL2
      integer iCode

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
      

      integer nArrMax
      parameter (nArrMax=200)   ! maximum size of arrays

      common /DataGJV/
     $     Arr(3,4,nArrMax),    ! 3* 4D-Vertizes
     $     EArr(6,nArrMax),     ! errFlag, rank
     $     verb,                ! verbosity
     $     AtOrigin             ! treatment of outmost prod points

      double precision Arr
      integer EArr
      integer verb
      logical AtOrigin

      save /DataGJV/

      integer nxi
      parameter (nxi=1000)

      double precision xiP(nxi),xiM(nxi)
      integer i,j,cL,iV,iString,iVec, ID,ID0
      double precision VV(4),
     $     Vert1(4),Vert2(4), VertP1(4),VertP2(4)
      double precision DPS(5)

      integer cL3,cL4
      data cL3,cL4 /0,0/
c      integer iL1,iL2,iL3,iL4
      integer iiCL,iH,ii

      common /DataRepString0/ RepFak,RepN,RepID,RepUse
      double precision RepFak(4,20)
      integer RepN(4), RepID(4,20), RepUse(4000)
      save /DataRepString0/

      integer OutmostCL(2)

      common /DataICL/ NiCL,iCL,NiCL4
      integer NiCL(2),iCL(2,100),NiCL4
      save /DataICL/

      logical lReverse, lInclude4

      integer dID
      integer j1,j2,k1,k2
      double precision hG

      integer JJold,KKold

      common /DataRevOrd/ ArrRO
      integer ArrRO(6,2,-1:1)
      save /DataRevOrd/

c      write(*,*) 'GJV_String_0 ', cL1,cL2

c... find last and next-to-last line:

      do cL=cL1,cL2
         if (abs(SFR_Flag(cL)).eq.3) cL3 = cL
         if (abs(SFR_Flag(cL)).eq.4) cL4 = cL
         if (SFR_Flag(cL).eq.0) then 
            write(*,*) 'GJV_String_0: FATAL ERROR, Line ',cL,
     $           ' should not be here! STOP!!!'
            call PYLIST(2)
            stop
         endif
      enddo

      EArr(1,cL3) = 3
      EArr(1,cL4) = 4

c... Check String Break Ordering, Set Reverse-Flag:
c
c Ich gehe die Liste durch und untersuche, ob das String-Breaking durch
c SFR_JJ, SFR_KK reported in der richtigen Reihenfolge ist. Das geht 
c natuerlich nur fuer QgQ-Strings (und GG-Strings, aber dort...), und dann
c auch nur in speziellen Faellen. 
c Ist evtl. SFR_REV oder SFR_REV2 eine Hilfe, das allgemeiner zu machen?

      lReverse = .FALSE.

      ID0 = SFR_ID0(cL1)

      JJold = ID0
      KKold = ID0+1

      do cL=cL1,cL2
         if (SFR_JJ(1,cL).eq.JJold
     $        .and. SFR_KK(1,cL).eq.KKold) then
            JJold = SFR_JJ(2,cL)
            KKold = SFR_KK(2,cL)
         else
            lReverse =.TRUE.
         endif
      enddo

      if (iCode.eq.3) lReverse = .FALSE. ! Gluonic String

c...FOLLOWING IS JUST FOR INTERNAL PURPOSES:

      if (lReverse) then

c$$$         call SFREP_ERROR(cl4,'Reverse!!!:')
c$$$         
c$$$
c$$$         if (SFR_REV2(cL4).eq.-1) then
c$$$            call SFREP_ERROR(cl4,'Reverse???:')
c$$$         endif

c$$$         if (SFR_Flag(cL4).lt.0) then
c$$$            if (cL4.lt.cL3) then 
c$$$               ArrRO(iCode,2,-1) = ArrRO(iCode,2,-1)+1
c$$$            else
c$$$               ArrRO(iCode,2,1) = ArrRO(iCode,2,1)+1
c$$$            endif
c$$$         else
c$$$            if (cL4.lt.cL3) then 
c$$$               ArrRO(iCode,2,1) = ArrRO(iCode,2,1)+1
c$$$            else
c$$$               ArrRO(iCode,2,-1) = ArrRO(iCode,2,-1)+1
c$$$            endif
c$$$         endif

c$$$         ArrRO(iCode,2,SFR_REV(cL4))
c$$$     $        = ArrRO(iCode,2,SFR_REV(cL4))+1

c$$$
c$$$         if (SFR_Flag(cL4).lt.SFR_Flag(cL3)) then
c$$$            ArrRO(iCode,2,-1)
c$$$     $           = ArrRO(iCode,2,-1)+1
c$$$         else
c$$$            ArrRO(iCode,2,1)
c$$$     $           = ArrRO(iCode,2,1)+1
c$$$         endif

         ArrRO(iCode,2,SFR_REV2(cL4))
     $        = ArrRO(iCode,2,SFR_REV2(cL4))+1
         
c$$$         if (SFR_Flag(cL4).lt.0) then
c$$$            ArrRO(iCode,2,-SFR_REV2(cL4))
c$$$     $           = ArrRO(iCode,2,-SFR_REV2(cL4))+1
c$$$         else
c$$$            ArrRO(iCode,2,SFR_REV2(cL4))
c$$$     $           = ArrRO(iCode,2,SFR_REV2(cL4))+1
c$$$         endif
      else

c$$$         if (SFR_Flag(cL4).lt.0) then
c$$$            if (cL4.lt.cL3) then 
c$$$               ArrRO(iCode,1,-1) = ArrRO(iCode,1,-1)+1
c$$$            else
c$$$               ArrRO(iCode,1,1) = ArrRO(iCode,1,1)+1
c$$$            endif
c$$$         else
c$$$            if (cL4.lt.cL3) then 
c$$$               ArrRO(iCode,1,1) = ArrRO(iCode,1,1)+1
c$$$            else
c$$$               ArrRO(iCode,1,-1) = ArrRO(iCode,1,-1)+1
c$$$            endif
c$$$         endif

c$$$         ArrRO(iCode,1,SFR_REV(cL4))
c$$$     $        = ArrRO(iCode,1,SFR_REV(cL4))+1

c$$$         if (SFR_Flag(cL4).lt.SFR_Flag(cL3)) then
c$$$            ArrRO(iCode,1,-1)
c$$$     $           = ArrRO(iCode,1,-1)+1
c$$$         else
c$$$            ArrRO(iCode,1,1)
c$$$     $           = ArrRO(iCode,1,1)+1
c$$$         endif

         ArrRO(iCode,1,SFR_REV2(cL4))
     $        = ArrRO(iCode,1,SFR_REV2(cL4))+1

c$$$         if (SFR_Flag(cL4).lt.0) then
c$$$            ArrRO(iCode,1,-SFR_REV2(cL4))
c$$$     $           = ArrRO(iCode,1,-SFR_REV2(cL4))+1
c$$$         else
c$$$            ArrRO(iCode,1,SFR_REV2(cL4))
c$$$     $           = ArrRO(iCode,1,SFR_REV2(cL4))+1
c$$$         endif
      endif

c... Exclude Particle 4 for GG-Strings

      lInclude4 = .TRUE.
      if (iCode.eq.3) lInclude4 = .FALSE.

c... find lines to do from left and to do from right:

      NiCL(1) = 0
      NiCL(2) = 0

      do cL=cL1,min(cL3,cL4)-1
         NiCL(1) = NiCL(1)+1
         iCL(1,NiCL(1)) = cL
      enddo
      do cL=cL2,max(cL3,cL4)+1,-1
         NiCL(2) = NiCL(2)+1
         iCL(2,NiCL(2)) = cL
      enddo

c... Set an index to get iString:
      if (NiCL(1).gt.0) then
         iiCL = iCL(1,1)
      else if (NiCL(2).gt.0) then
         iiCL = iCL(2,1)
      else
c... (ignore the problems, set iiCL=cL3)                                
c         call SFREP_ERROR(90,'Set an index to get iString')
c         write (*,*) cL1,CL2,cL3,cL4,NiCL(1),NiCL(2)
c         return
         iiCL = cL3
      endif
      call SFREPS_FindVecI(2, iiCL, 5, iString,iVec)

c... now add Line with 3 and 4

      iH = 1
      if (cL3.gt.cL4) iH = 2

      if (lReverse) then
         NiCL(iH) = NiCL(iH)+1
         NiCL4 = NiCL(iH)
         iCL(iH,NiCL(iH)) = cL4
      endif 
      NiCL(iH) = NiCL(iH)+1
      iCL(iH,NiCL(iH)) = cL3
      if (.not.lReverse) then
         NiCL(iH) = NiCL(iH)+1
         NiCL4 = NiCL(iH)
         iCL(iH,NiCL(iH)) = cL4
      endif 

c... Set the Outmost particles:

      if (NiCL(1).gt.0 .and. NiCL(2).gt.0) then
         OutmostCL(1) = iCL(1,1)       
         OutmostCL(2) = iCL(2,1)       
      elseif (NiCL(1).gt.0) then
         OutmostCL(1) = iCL(1,1)
         OutmostCL(2) = iCL(1,NiCL(1))
      elseif (NiCL(2).gt.0) then
         OutmostCL(1) = iCL(2,NiCL(2))
         OutmostCL(2) = iCL(2,1)  
      else
         ! this should never be possible
      endif

c... Delete Particle 4 from the list:

      if (.not.lInclude4) then
         do i=NiCL4,NiCL(iH)
            iCL(iH,i) = iCL(iH,i+1)
         enddo
         NiCL(iH) = NiCL(iH)-1
      endif

c... Set the rank of the particles:

      do ii=1,2
         iH = NiCL(3-ii)
         if (.not. lInclude4) iH = iH+1
         do iiCL = 1,NiCL(ii)
            cL = iCL(ii,iiCL)
            EArr(5,cL) = iiCL
            EArr(6,cL) = min(iiCL,NiCL(ii)-iiCL+1+iH)
         enddo
      enddo
      

c... LOOP OVER THE 'LEFT' LINES:

      call SetArrayD(nxi,xiP,1d0) ! reset all xi-Values
      call SetArrayD(nxi,xiM,0d0)

      do iiCL = 1,NiCL(1)
         cL = iCL(1,iiCL)

         call SetArrayD(4,Vert1, 0d0)
         call SetArrayD(4,Vert2, 0d0)
         call SetArrayD(4,VertP1, 0d0)
         call SetArrayD(4,VertP2, 0d0)

         call SetArrayI(4,RepN,0)
         call SetArrayI(4000,RepUse,0)
         do j=SFR_JJ(1,cL),SFR_KK(1,cL)-1,4
            RepUse(j) = 1
            RepUse(j+1) = 1
         enddo
         do j=SFR_JJ(2,cL),SFR_KK(2,cL)-1,4
            RepUse(j) = 1
            RepUse(j+1) = 1
         enddo

         ID0 = SFR_ID0(cL)

         if (verb.gt.0) then
            write(*,3001) (j,xiP(j),j=SFR_JJ(1,cL),SFR_KK(2,cL)-1,4)
            write(*,3002) (j,xiM(j),j=SFR_JJ(1,cL)+1,SFR_KK(2,cL),4)
         endif

         do iV=5,SFR_NVec(cL)
            ID = SFR_VecID(cL,iV)
            if (mod((ID-ID0),2).eq.0) then 
c... (+)-Vector:

               if (ID.le.SFR_KK(1,cL)) then
                  call AddVertSFR(Vert1,xiP(ID),cL,iV)
                  call AddVert_Rep(1,ID,xiP(ID))
c$$$               else
c$$$                  write(*,2000) 'Wegfall 1+',xiP(ID),ID
               endif

               xiP(ID) = xiP(ID)-SFR_Coeff(cL,iV)

               if (ID.ge.SFR_JJ(2,cL)) then
                  call AddVertSFR(Vert2,xiP(ID),cL,iV)
                  call AddVert_Rep(2,ID,xiP(ID))
c$$$               else
c$$$                  write(*,2000) 'Wegfall 2+',xiP(ID),ID
               endif
               call AddVertSFR(VertP1,SFR_Coeff(cL,iV),cL,iV)
               call AddVert_Rep(3,ID,SFR_Coeff(cL,iV))

            else             
c... (-)-Vector: 

               if (ID.le.SFR_KK(1,cL)) then
                  call AddVertSFR(Vert1,xiM(ID),cL,iV)
                  call AddVert_Rep(1,ID,xiM(ID))
c$$$               else
c$$$                  write(*,2000) 'Wegfall 1-',xiM(ID),ID
               endif

               xiM(ID) = xiM(ID)+SFR_Coeff(cL,iV)

               if (ID.ge.SFR_JJ(2,cL)) then
                  call AddVertSFR(Vert2,xiM(ID),cL,iV)
                  call AddVert_Rep(2,ID,xiM(ID))
c$$$               else
c$$$                  write(*,2000) 'Wegfall 2-',xiM(ID),ID
               endif
               call AddVertSFR(VertP2,SFR_Coeff(cL,iV),cL,iV)
               call AddVert_Rep(4,ID,SFR_Coeff(cL,iV))

            endif
         enddo ! iV

c... Some (brute force) corrections:

         if (iCode.eq.3) goto 199 ! no correct for GG-String (not understood)

c... (zu wenige Vektoren:) [z.B. in (1,2)(1,2)]

         j1 = (SFR_JJ(1,cL)-ID0)/4+1
         k1 = (SFR_KK(1,cL)-1-ID0)/4+1
         j2 = (SFR_JJ(2,cL)-ID0)/4+1
         k2 = (SFR_KK(2,cL)-1-ID0)/4+1
         if (j1.eq.j2 .and. k1.eq.k2 .and. k1.gt.j1) then
            ID = SFR_JJ(1,cL)+1
            dID = 3

 110        continue
            call SFREPS_GetVecID(2,iString,iVec,VV,ID)

            call AddVert_Rep(1,ID,0.5d0)
            call AddVert_Rep(2,ID,0.5d0)

            call AddVert(Vert1,0.5d0,VV)
            call AddVert(Vert2,0.5d0,VV)

            ID = ID+dID
            dID = 4-dID ! 3,1,3,1,...
            if (ID.lt.SFR_KK(1,cL)) goto 110
         endif

c... (was noch nicht erschlagen ist [zu wenige Vektoren], kommt jetzt:)

         j1 = min(SFR_JJ(1,cL),SFR_JJ(2,cL))
         j2 = max(SFR_KK(1,cL),SFR_KK(2,cL))
         do ID=j1,j2
            if (RepUse(ID).ne.0) then
               dID = 0
               if (RepUse(ID+3).ne.0) dID=ID+3
               call SFREPS_GetVecID(2,iString,iVec,VV,ID)

               call AddVert_Rep(1,ID,0.5d0)
               call AddVert_Rep(2,ID,0.5d0)

               call AddVert(Vert1,0.5d0,VV)
               call AddVert(Vert2,0.5d0,VV)
               if (dID.ne.0) RepUse(dID) = 1
            endif
         enddo

 199     continue

c... Set the output arrays:

         do j=1,4
            Arr(1,j,cL) = Vert1(j)
            Arr(2,j,cL) = Vert2(j)
            Arr(3,j,cL) = VertP1(j)+Vert2(j)
         enddo
         EArr(2,cL) = 0         ! Error Flag

c---------------------------
         if (verb.gt.0) then

         write(*,*) '****',cL,' Left <<<<<<<<<<<<<<'

         write(*,*) '**** ... ',SFR_JJ(1,cL),SFR_KK(1,cL),
     $        SFR_JJ(2,cL),SFR_KK(2,cL),'  ID0=',ID0

         write(*,3001) (j,xiP(j),j=SFR_JJ(1,cL),SFR_KK(2,cL)-1,4)
         write(*,3002) (j,xiM(j),j=SFR_JJ(1,cL)+1,SFR_KK(2,cL),4)

         call VertLength(Vert1,hG)
         write(*,1000) 'Vert1   ',(Vert1(j),j=1,4),hG
         call VertLength(Vert2,hG)
         if (SFR_Gamma(cL).gt.0d0) then
         write(*,1000) 'Vert2   ',(Vert2(j),j=1,4),hG,hG-SFR_Gamma(cL)
         else
         write(*,1000) 'Vert2   ',(Vert2(j),j=1,4),hG
         endif
         write(*,1000) 'Vert1+P2',(Vert1(j)+VertP2(j),j=1,4)
         write(*,1000) 'Vert2+P1',(Vert2(j)+VertP1(j),j=1,4)
         write(*,*)

         write(*,2000) 'Vert1 ',(RepFak(1,j),RepID(1,j),j=1,RepN(1))
         write(*,2000) 'Vert2 ',(RepFak(2,j),RepID(2,j),j=1,RepN(2))
         write(*,2000) 'VertP1',(RepFak(3,j),RepID(3,j),j=1,RepN(3))
         write(*,2000) 'VertP2',(RepFak(4,j),RepID(4,j),j=1,RepN(4))
         write(*,*)

         endif ! verb
c---------------------------
      enddo ! cL

c... LOOP OVER THE 'RIGHT' LINES:

      call SetArrayD(nxi,xiP,0d0) ! reset all xi-Values
      call SetArrayD(nxi,xiM,1d0)

      do iiCL = 1,NiCL(2)
         cL = iCL(2,iiCL)

         call SetArrayD(4,Vert1, 0d0)
         call SetArrayD(4,Vert2, 0d0)
         call SetArrayD(4,VertP1, 0d0)
         call SetArrayD(4,VertP2, 0d0)

         call SetArrayI(4,RepN,0)
         call SetArrayI(4000,RepUse,0)
         do j=SFR_JJ(1,cL),SFR_KK(1,cL)-1,4
            RepUse(j) = 1
            RepUse(j+1) = 1
         enddo
         do j=SFR_JJ(2,cL),SFR_KK(2,cL)-1,4
            RepUse(j) = 1
            RepUse(j+1) = 1
         enddo

         ID0 = SFR_ID0(cL)

         if (verb.gt.0) then
            write(*,3001) (j,xiP(j),j=SFR_JJ(1,cL),SFR_KK(2,cL)-1,4)
            write(*,3002) (j,xiM(j),j=SFR_JJ(1,cL)+1,SFR_KK(2,cL),4)
         endif

         do iV=5,SFR_NVec(cL)
            ID = SFR_VecID(cL,iV)
            if (mod((ID-ID0),2).eq.0) then 
c... (+)-Vector:

               if (ID.ge.SFR_JJ(2,cL)) then
                  call AddVertSFR(Vert1,xiP(ID),cL,iV)
                  call AddVert_Rep(1,ID,xiP(ID))
c$$$               else
c$$$                  write(*,2000) 'WEGFALL 1+',xiP(ID),ID
               endif

               xiP(ID) = xiP(ID)+SFR_Coeff(cL,iV)

               if (ID.le.SFR_KK(1,cL)) then
                  call AddVertSFR(Vert2,xiP(ID),cL,iV)
                  call AddVert_Rep(2,ID,xiP(ID))
c$$$               else
c$$$                  write(*,2000) 'WEGFALL 2+',xiP(ID),ID
               endif
               call AddVertSFR(VertP1,SFR_Coeff(cL,iV),cL,iV)
               call AddVert_Rep(3,ID,SFR_Coeff(cL,iV))
               
            else             
c... (-)-Vector: 

               if (ID.ge.SFR_JJ(2,cL)) then
                  call AddVertSFR(Vert1,xiM(ID),cL,iV)
                  call AddVert_Rep(1,ID,xiM(ID))
c$$$               else
c$$$                  write(*,2000) 'WEGFALL 1-',xiM(ID),ID
               endif

               xiM(ID) = xiM(ID)-SFR_Coeff(cL,iV)

               if (ID.le.SFR_KK(1,cL)) then
                  call AddVertSFR(Vert2,xiM(ID),cL,iV)
                  call AddVert_Rep(2,ID,xiM(ID))
c$$$               else
c$$$                  write(*,2000) 'WEGFALL 2-',xiM(ID),ID
               endif
               call AddVertSFR(VertP2,SFR_Coeff(cL,iV),cL,iV)
               call AddVert_Rep(4,ID,SFR_Coeff(cL,iV))
               
            endif
         enddo ! iV

c... Some (brute force) corrections:

         if (iCode.eq.3) goto 299 ! no correct for GG-String (not understood)

c... (zu wenige Vektoren:) [z.B. in (1,2)(1,2)]

         j1 = (SFR_JJ(1,cL)-ID0)/4+1
         k1 = (SFR_KK(1,cL)-1-ID0)/4+1
         j2 = (SFR_JJ(2,cL)-ID0)/4+1
         k2 = (SFR_KK(2,cL)-1-ID0)/4+1
         if (j1.eq.j2 .and. k1.eq.k2 .and. k1.gt.j1) then
            ID = SFR_JJ(1,cL)+1
            dID = 3

 210        continue
            call SFREPS_GetVecID(2,iString,iVec,VV,ID)

            call AddVert_Rep(1,ID,0.5d0)
            call AddVert_Rep(2,ID,0.5d0)

            call AddVert(Vert1,0.5d0,VV)
            call AddVert(Vert2,0.5d0,VV)

            ID = ID+dID
            dID = 4-dID ! 3,1,3,1,...
            if (ID.lt.SFR_KK(1,cL)) goto 210
         endif

c... (was noch nicht erschlagen ist [zu wenige Vektoren], kommt jetzt:)

         j1 = min(SFR_JJ(1,cL),SFR_JJ(2,cL))
         j2 = max(SFR_KK(1,cL),SFR_KK(2,cL))
         do ID=j1,j2
            if (RepUse(ID).ne.0) then
               dID = 0
               if (RepUse(ID+3).ne.0) dID=ID+3
               call SFREPS_GetVecID(2,iString,iVec,VV,ID)

               call AddVert_Rep(1,ID,0.5d0)
               call AddVert_Rep(2,ID,0.5d0)

               call AddVert(Vert1,0.5d0,VV)
               call AddVert(Vert2,0.5d0,VV)
               if (dID.ne.0) RepUse(dID) = 1
            endif
         enddo

 299     continue

c... Set the output arrays:

         do j=1,4
            Arr(1,j,cL) = Vert2(j)
            Arr(2,j,cL) = Vert1(j)
            Arr(3,j,cL) = VertP1(j)+Vert1(j)
         enddo
         EArr(2,cL) = 0         ! Error Flag

c---------------------------
         if (verb.gt.0) then

         write(*,*) '****',cL,' Right >>>>>>>>>>>>>'

         write(*,*) '**** ... ',SFR_JJ(1,cL),SFR_KK(1,cL),
     $        SFR_JJ(2,cL),SFR_KK(2,cL),'  ID0=',ID0

         if (mod(SFR_JJ(1,cL)-ID0,2).eq.1) write(*,*) 'oops 1'
         if (mod(SFR_JJ(2,cL)-ID0,2).eq.1) write(*,*) 'oops 2'

         write(*,3001) (j,xiP(j),
     $        j=SFR_JJ(1,cL),SFR_KK(2,cL)-1,4)
         write(*,3002) (j,xiM(j),
     $        j=SFR_JJ(1,cL)+1,SFR_KK(2,cL),4)

         call VertLength(Vert1,hG)
         write(*,1000) 'Vert1   ',(Vert1(j),j=1,4),hG
         call VertLength(Vert2,hG)
         if (SFR_Gamma(cL).gt.0d0) then
         write(*,1000) 'Vert2   ',(Vert2(j),j=1,4),hG,hG-SFR_Gamma(cL)
         else
         write(*,1000) 'Vert2   ',(Vert2(j),j=1,4),hG
         endif
         write(*,1000) 'Vert1+P1',(Vert1(j)+VertP1(j),j=1,4)
         write(*,1000) 'Vert2+P2',(Vert2(j)+VertP2(j),j=1,4)
         write(*,*)

         write(*,2000) 'Vert1 ',(RepFak(1,j),RepID(1,j),j=1,RepN(1))
         write(*,2000) 'Vert2 ',(RepFak(2,j),RepID(2,j),j=1,RepN(2))
         write(*,2000) 'VertP1',(RepFak(3,j),RepID(3,j),j=1,RepN(3))
         write(*,2000) 'VertP2',(RepFak(4,j),RepID(4,j),j=1,RepN(4))
         write(*,*)

         endif ! verb
c---------------------------

      enddo ! cL

c... Set the outmost production points:

      if (AtOrigin) then
         do j=1,4
            Arr(1,j,OutmostCL(1)) = 0d0
            Arr(2,j,OutmostCL(2)) = 0d0
         enddo
      endif

c... now boost all the vectors back:

c... COMMENT THIS FOR TEMPORARY USE ONLY !!!!!!!!

      call SFREPS_GetBoost(iString,DPS)
      call BoostVert(DPS,cL1,cL2)

c... ENDE

      return                    ! remove this for the tests below

c-------------------------------------------------------------------
c... Fuer Tests:
c... EVENTS, BEI DENEN ABGEBROCHEN WERDEN SOLLTE:

      if (iCode.eq.1) return    ! dont not check (1+1)-strings
      if (iCode.eq.3) return    ! dont not check GG-strings

c...1) Hadron 1 auf (1,1)(2,2)

      cL = iCL(1,1)
      j1 = (SFR_JJ(1,cL)-SFR_ID0(cL))/4+1
      k1 = (SFR_KK(1,cL)-1-SFR_ID0(cL))/4+1
      j2 = (SFR_JJ(2,cL)-SFR_ID0(cL))/4+1
      k2 = (SFR_KK(2,cL)-1-SFR_ID0(cL))/4+1

      if(j1.eq.1 .and. k1.eq.1 .and. j2.eq.2 .and. k2.eq.2) then
         call SFREP_ERROR(99,'Hadron 1 auf (1,1)(2,2)')
      endif

c...2) Hadron auf String (i,i+1)(i,i+1)

c$$$      do cL=iL1,iL2
c$$$         j1 = (SFR_JJ(1,cL)-SFR_ID0(cL))/4+1
c$$$         k1 = (SFR_KK(1,cL)-1-SFR_ID0(cL))/4+1
c$$$         j2 = (SFR_JJ(2,cL)-SFR_ID0(cL))/4+1
c$$$         k2 = (SFR_KK(2,cL)-1-SFR_ID0(cL))/4+1
c$$$         if (j1.eq.j2 .and. k1.eq.k2 .and. k1.gt.j1) then
c$$$            call SFREP_ERROR(99,
c$$$     $           'Hadron (von links) auf String (i,i+..)(i,i+..)')
c$$$            write(*,*) 'Hadron = ',cL
c$$$         endif
c$$$      enddo
c$$$      do cL=iL3,iL4
c$$$         j1 = (SFR_JJ(1,cL)-SFR_ID0(cL))/4+1
c$$$         k1 = (SFR_KK(1,cL)-1-SFR_ID0(cL))/4+1
c$$$         j2 = (SFR_JJ(2,cL)-SFR_ID0(cL))/4+1
c$$$         k2 = (SFR_KK(2,cL)-1-SFR_ID0(cL))/4+1
c$$$         if (j1.eq.j2 .and. k1.eq.k2 .and. k1.gt.j1) then
c$$$            call SFREP_ERROR(99,
c$$$     $           'Hadron (von rechts) auf String (i,i+..)(i,i+..)')
c$$$            write(*,*) 'Hadron = ',cL
c$$$         endif
c$$$      enddo

c... 3) |Gamma - Gamma_reported| > 1e-1
      do iiCL = 1,NiCL(1)
         cL = iCL(1,iiCL)

         j1 = (SFR_JJ(1,cL)-SFR_ID0(cL))/4+1
         k1 = (SFR_KK(1,cL)-1-SFR_ID0(cL))/4+1
         j2 = (SFR_JJ(2,cL)-SFR_ID0(cL))/4+1
         k2 = (SFR_KK(2,cL)-1-SFR_ID0(cL))/4+1
         if (j1.eq.j2.and.k1.eq.k2.and.j1.eq.k1) then
            ! this is like (1+1) string, dont check                   !
         else
            hG = Arr(2,4,cL)**2
     $           -Arr(2,1,cL)**2-Arr(2,2,cL)**2-Arr(2,3,cL)**2
            if (SFR_Gamma(cL).ne.0
     $           .and. abs(hG-SFR_Gamma(cL)).gt.1d-1) then
               call SFREP_ERROR(99,
     $              '|Gamma - Gamma_reported| > 1e-1 (a)')
               write(*,*) 'Hadron = ',cL
            endif
         endif
      enddo
      do iiCL = 1,NiCL(2)
         cL = iCL(2,iiCL)

         j1 = (SFR_JJ(1,cL)-SFR_ID0(cL))/4+1
         k1 = (SFR_KK(1,cL)-1-SFR_ID0(cL))/4+1
         j2 = (SFR_JJ(2,cL)-SFR_ID0(cL))/4+1
         k2 = (SFR_KK(2,cL)-1-SFR_ID0(cL))/4+1
         if (j1.eq.j2.and.k1.eq.k2.and.j1.eq.k1) then
            ! this is like (1+1) string, dont check                   !
         else
            hG = Arr(1,4,cL)**2
     $           -Arr(1,1,cL)**2-Arr(1,2,cL)**2-Arr(1,3,cL)**2
            if (SFR_Gamma(cL).ne.0
     $           .and. abs(hG-SFR_Gamma(cL)).gt.1d-1) then
               call SFREP_ERROR(99,
     $              '|Gamma - Gamma_reported| > 1e-1 (b)')
               write(*,*) 'Hadron = ',cL
            endif
         endif
      enddo

 1000 FORMAT('GJV_String_0',':',A,'= ',1P,4e13.5,' : ',2e13.5)
 2000 FORMAT('GJV_String_0',':',A,'= ',1P,10(e13.5,'{',i2,'}'))

 3001 FORMAT('***',1P,10(' xiP(',i2,')=',e13.5))
 3002 FORMAT('***',1P,10(' xiM(',i2,')=',e13.5))    

      end
c*****************************************************************==
c Do the "Shrinking" a PYEDIT call would do
c
c This has to be done, BEFORE a call to PYEDIT is done!
c
      subroutine GetJetsetVecPYEDIT()
      IMPLICIT NONE

      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      integer N,NPAD,K
      double precision P,V
      SAVE /PYJETS/

      integer nArrMax
      parameter (nArrMax=200)   ! maximum size of arrays

      common /DataGJV/
     $     Arr(3,4,nArrMax),    ! 3* 4D-Vertizes
     $     EArr(6,nArrMax),     ! errFlag, rank
     $     verb,                ! verbosity
     $     AtOrigin             ! treatment of outmost prod points

      double precision Arr
      integer EArr
      integer verb
      logical AtOrigin

      save /DataGJV/

      integer iOld,iNew, i,j

      iNew = 0

      do iOld=1,N
         if (K(iOld,1).lt.10) then
            iNew=iNew+1

            do j=1,6
               EArr(j,iNew) = EArr(j,iOld)
            enddo
            do j=1,4
               do i=1,3
                  Arr(i,j,iNew) = Arr(i,j,iOld)
               enddo
            enddo
            
         endif
      enddo

      end
c*****************************************************************==
c Boost the output array like PYROBO
c
c just uses entries N+1 to N+... and performs a PYROBO
c
      subroutine GetJetsetVecPYROBO(THE,PHI,BEX,BEY,BEZ)
      IMPLICIT NONE
      double precision THE,PHI,BEX,BEY,BEZ

      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      integer N,NPAD,K
      double precision P,V
      SAVE /PYJETS/

      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      integer MSTU,MSTJ
      double precision PARU,PARJ
      SAVE /PYDAT1/

      COMMON/LUDAT1/MSTUf(200),PARUf(200),MSTJf(200),PARJf(200)
      integer MSTUf,MSTJf
      real PARUf,PARJf
      save/LUDAT1/

      integer nArrMax
      parameter (nArrMax=200)   ! maximum size of arrays

      common /DataGJV/
     $     Arr(3,4,nArrMax),    ! 3* 4D-Vertizes
     $     EArr(6,nArrMax),     ! errFlag, rank
     $     verb,                ! verbosity
     $     AtOrigin             ! treatment of outmost prod points

      double precision Arr
      integer EArr
      integer verb
      logical AtOrigin

      save /DataGJV/

      common /DataGJVerror/ UseWithPYTHIA, ListFinalEvent
      integer UseWithPYTHIA, ListFinalEvent
      save /DataGJVerror/

      integer iLA, iLJ, j

      if (UseWithPythia > 0) then
         MSTU(33) = 0           ! do not reset V
      else
         MSTUf(33) = 0           ! do not reset V
      endif

c... copy vectors to Jetset-arrays

      iLJ = N
      do iLA=1,N
         K(iLJ+1,1) = 1         ! assume to be boosted!
         K(iLJ+2,1) = 1
         do j=1,4
            P(iLJ+1,j) = Arr(1,j,iLA)
            V(iLJ+1,j) = Arr(2,j,iLA)
            P(iLJ+2,j) = Arr(3,j,iLA)
         enddo
         iLJ = iLJ+2
      enddo

c... perform boost
      if (UseWithPythia > 0) then
         call PYROBO(N+1,iLJ, THE,PHI,BEX,BEY,BEZ)
      else
         call LUDBRB(N+1,iLJ, THE,PHI,BEX,BEY,BEZ)
      endif

c... copy vectors back

      iLJ = N
      do iLA=1,N
         do j=1,4
            Arr(1,j,iLA) = P(iLJ+1,j)
            Arr(2,j,iLA) = V(iLJ+1,j)
            Arr(3,j,iLA) = P(iLJ+2,j)
         enddo
         iLJ = iLJ+2
      enddo

      end

c*****************************************************************==
c Check whether returned time makes sense, adjust the error flag array:
c 1) Error, if Gamma(1...3) < GammaLevel (approx. -1d-10)
c 2) Error, if Gamma(3) < max(Gamma(1).Gamma(2))
c
c Tests are only done for particles without any error flag set at input.
c
      subroutine GetJetsetVecCheckT(GammaLevel)
      IMPLICIT NONE
      double precision GammaLevel

      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      integer N,NPAD,K
      double precision P,V
      SAVE /PYJETS/

      integer nArrMax
      parameter (nArrMax=200)   ! maximum size of arrays

      common /DataGJV/
     $     Arr(3,4,nArrMax),    ! 3* 4D-Vertizes
     $     EArr(6,nArrMax),     ! errFlag, rank
     $     verb,                ! verbosity
     $     AtOrigin             ! treatment of outmost prod points

      double precision Arr
      integer EArr
      integer verb
      logical AtOrigin

      save /DataGJV/

      integer i,j
      double precision Gamma(3),GammaL
      integer iFl(3)

      GammaL = -1d-5
      if(GammaLevel.gt.GammaL) then
         write(*,*) 'GetJetsetVecCheckT: GammaLevel:',
     $        GammaLevel,'>=',GammaL
         write(*,*) '...please check code! STOP!'
         stop
      endif
      GammaL=GammaLevel

      do i=1,N
         if (EArr(1,i).eq.0) goto 199 ! only processed 
         if (EArr(2,i).ne.0) goto 199 ! already problems with times

         do j=1,3
            iFl(j) = 0
         enddo
            
         do j=1,3
            Gamma(j) = 
     $           Arr(j,4,i)**2-Arr(j,1,i)**2-Arr(j,2,i)**2-Arr(j,3,i)**2
            if (Gamma(j).lt.GammaL) iFl(j) = 1
         enddo

         if (Gamma(3).lt.Gamma(1).or.Gamma(3).lt.Gamma(2))
     $        iFl(3) = 1
         
         EArr(2,i) = 4*iFl(3)+2*iFl(2)+iFl(1)

 199     continue
      enddo

      end

c*****************************************************************==
c List the information 
c

      subroutine GetJetSetVec_List(iFile, i1,i2)
      IMPLICIT NONE

      integer iFile,i1,i2

      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      integer N,NPAD,K
      double precision P,V
      SAVE /PYJETS/

      integer nArrMax
      parameter (nArrMax=200)   ! maximum size of arrays

      common /DataGJV/
     $     Arr(3,4,nArrMax),    ! 3* 4D-Vertizes
     $     EArr(6,nArrMax),     ! errFlag, rank
     $     verb,                ! verbosity
     $     AtOrigin             ! treatment of outmost prod points

      double precision Arr
      integer EArr
      integer verb
      logical AtOrigin

      save /DataGJV/

      integer iArr,ii,jj,ii1,ii2,j
      double precision Gamma(3), sGamma(3)
      

      ii1 = i1
      ii2 = i2

      if (i1.le.0) ii1 = 1
      if (i2.le.0) ii2 = N

      do iArr=ii1,ii2
         write(iFile,1000) iArr,((Arr(ii,jj,iArr),jj=1,4),ii=1,3),
     $        (eArr(ii,iArr),ii=1,6)

         do j=1,3
            Gamma(j) = 
     $           Arr(j,4,iArr)**2-Arr(j,1,iArr)**2
     $           -Arr(j,2,iArr)**2-Arr(j,3,iArr)**2
            

            if (Gamma(j).gt.0.0) then
               sGamma(j) = sqrt(Gamma(j))
            elseif (Gamma(j).gt.-1d-5) then
               sGamma(j) = 0.0
            else
               sGamma(j) = -99d99
            endif
         enddo

         write(iFile,1020) ((/Gamma(ii),sGamma(ii)/),ii=1,3)

      enddo
      
 1000 FORMAT (i3,3('(',4f12.3,')'),6i3)
 1020 FORMAT ('   ',3('( --> (',f14.5,')^1/2 = ',f14.5,6(' '),')'))
      end      

c===================================================================
c-------------------------------------------------------------------
c LOW LEVEL: Add a Vector to DataRepString0

      subroutine AddVert_Rep(iV, ID, Fak)
      IMPLICIT NONE
      integer iV, ID
      double precision Fak

      double precision RepFak(4,20)
      integer RepN(4), RepID(4,20), RepUse(4000)
      common /DataRepString0/ RepFak,RepN,RepID,RepUse
      save /DataRepString0/

c      if (Fak.eq.0d0) return

      RepN(iV) = RepN(iV)+1
      if (RepN(iV).gt.20) then 
         call SFREP_Error(10,'AddVert_Rep: too many vectors')
         return
      endif
      RepID(iV,RepN(iV)) = ID
      RepFak(iV,RepN(iV)) = Fak
      RepUse(ID) = 0
      RepUse(ID+3) = 0
      RepUse(ID-3) = 0
      end 
     
c-------------------------------------------------------------------
c LOW LEVEL: Set the Error Array
c
      subroutine SetEArr(cL, i6,i5,i4,i3,i2,i1)
      IMPLICIT NONE
      integer cL,i6,i5,i4,i3,i2,i1

      integer nArrMax
      parameter (nArrMax=200)   ! maximum size of arrays

      common /DataGJV/
     $     Arr(3,4,nArrMax),    ! 3* 4D-Vertizes
     $     EArr(6,nArrMax),     ! errFlag, rank
     $     verb,                ! verbosity
     $     AtOrigin             ! treatment of outmost prod points

      double precision Arr
      integer EArr
      integer verb
      logical AtOrigin

      save /DataGJV/

      EArr(6,cL) = i6
      EArr(5,cL) = i5
      EArr(4,cL) = i4
      EArr(3,cL) = i3
      EArr(2,cL) = i2
      EArr(1,cL) = i1
      
      end

c-------------------------------------------------------------------
c LOW LEVEL: Boost the Vertices as Jetset would do
c
      subroutine BoostVert(DPS,cL1,cL2)
      IMPLICIT NONE
      double precision DPS(5)
      integer cL1,cL2

      integer nArrMax
      parameter (nArrMax=200)   ! maximum size of arrays

      common /DataGJV/
     $     Arr(3,4,nArrMax),    ! 3* 4D-Vertizes
     $     EArr(6,nArrMax),     ! errFlag, rank
     $     verb,                ! verbosity
     $     AtOrigin             ! treatment of outmost prod points

      double precision Arr
      integer EArr
      integer verb
      logical AtOrigin

      save /DataGJV/

      double precision DBX,DBY,DBZ,DB,EPS1,DGA,DV(4),DBV,DGABV
      double precision HHBZ,HHPMT,HHPEZ

      integer cL,j, iT

      if (ABS(DPS(3)).LT.0.99D0*DPS(4)) then

         DBX = DPS(1)/DPS(4)
         DBY = DPS(2)/DPS(4)
         DBZ = DPS(3)/DPS(4)

         IF(DBX**2+DBY**2+DBZ**2.GT.1D-20) THEN

c like PYROBO(...,...,0D0,0D0,DPS(1)/DPS(4),DPS(2)/DPS(4),DPS(3)/DPS(4))

            DB=SQRT(DBX**2+DBY**2+DBZ**2)
            EPS1=1D0-1D-12
            IF(DB.GT.EPS1) THEN
               CALL PYERRM(3,'(BoostVert:) boost vector too large')
               DBX=DBX*(EPS1/DB)
               DBY=DBY*(EPS1/DB)
               DBZ=DBZ*(EPS1/DB)
               DB=EPS1
            ENDIF
            DGA=1D0/SQRT(1D0-DB**2)

            do cL=cL1,cL2
               do iT=1,3
                  do j=1,4
                     DV(j)=Arr(iT,j,cL)
                  enddo
                  DBV=DBX*DV(1)+DBY*DV(2)+DBZ*DV(3)
                  DGABV=DGA*(DGA*DBV/(1D0+DGA)+DV(4))
                  Arr(iT,1,cL) = DV(1)+DGABV*DBX
                  Arr(iT,2,cL) = DV(2)+DGABV*DBY
                  Arr(iT,3,cL) = DV(3)+DGABV*DBZ
                  Arr(iT,4,cL) = DGA*(DV(4)+DBV)
               enddo
            enddo
         ENDIF

      else

c make shift in z-direction (i.e. recalculate rapidity y)

         HHBZ=SQRT(MAX(1D-6,DPS(4)+DPS(3))/MAX(1D-6,DPS(4)-DPS(3)))
         do cL=cL1,cL2
            do iT=1,3
               HHPMT = max(1d-10,Arr(iT,4,cL)**2-Arr(iT,3,cL)**2)
               IF(Arr(iT,3,cL).GT.0D0) THEN
                  HHPEZ=(Arr(iT,4,cL)+Arr(iT,3,cL))*HHBZ
                  if (HHPEZ.lt.1d-10) goto 100
                  Arr(iT,3,cL)=0.5D0*(HHPEZ-HHPMT/HHPEZ)
                  Arr(iT,4,cL)=0.5D0*(HHPEZ+HHPMT/HHPEZ)
               ELSE
                  HHPEZ=(Arr(iT,4,cL)-Arr(iT,3,cL))/HHBZ
                  if (HHPEZ.lt.1d-10) goto 100
                  Arr(iT,3,cL)=-0.5D0*(HHPEZ-HHPMT/HHPEZ)
                  Arr(iT,4,cL)=0.5D0*(HHPEZ+HHPMT/HHPEZ)
               ENDIF
 100           continue
            enddo
         enddo

         return
      endif

      end

c-------------------------------------------------------------------
c LOW LEVEL: Get to LightLike Vectors
c
c Transforms 2 massive vectors (from PYJETS) and returns two 
c lightlike vectors
c (Cut out from Pythia 6.208, PYSTRF, loop with label 670)
c
      subroutine TransformLL(iL1,iL2,PLL)
      IMPLICIT NONE
      integer iL1, iL2
      double precision PLL(2,4)

      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      integer N,NPAD,K
      double precision P,V
      SAVE /PYJETS/

      double precision DHKC, DHKS, DHK1,DHK2, DHM1,DHM2
      integer j

      DHM1 = P(iL1,5)**2
      DHM2 = P(iL2,5)**2
      DHKC = P(iL1,4)*P(iL2,4)
     $     -P(iL1,1)*P(iL2,1)-P(iL1,2)*P(iL2,2)-P(iL1,3)*P(iL2,3)
      DHKS = sqrt(DHKC**2-DHM1*DHM2)
      DHK1 = 0.5D0*((DHM2+DHKC)/DHKS-1D0)
      DHK2 = 0.5D0*((DHM1+DHKC)/DHKS-1D0)
      do j=1,4
         PLL(1,j) = (1d0+DHK1)* P(iL1,j)-DHK2* P(iL2,j)
         PLL(2,j) = (1d0+DHK2)* P(iL2,j)-DHK1* P(iL1,j)
      enddo

c$$$      write(*,*) '=== TransformLL ==='
c$$$      write(*,1000) 'P1  ',(P(iL1,j),j=1,5) 
c$$$      write(*,1000) 'P2  ',(P(iL2,j),j=1,5) 
c$$$      write(*,1000) 'Sum ',(P(iL1,j)+P(iL2,j),j=1,4)
c$$$
c$$$      write(*,1000) 'PLL1',(PLL(1,j),j=1,4),
c$$$     $     PLL(1,4)**2-PLL(1,1)**2-PLL(1,2)**2-PLL(1,3)**2
c$$$      write(*,1000) 'PLL2',(PLL(2,j),j=1,4),
c$$$     $     PLL(2,4)**2-PLL(2,1)**2-PLL(2,2)**2-PLL(2,3)**2
c$$$      write(*,1000) 'Sum ',(PLL(1,j)+PLL(2,j),j=1,4)
c$$$      stop

 1000 FORMAT(A,' : ',1P,5e13.5)
      end

c===================================================================

c=================================================================
c ERROR AND MESSAGING ROUTINES
c=================================================================
c
c emit a 'GetJetsetVec:Error'--Message, list the event and stop.
c
      subroutine GJV_ERR(Text,iLevel,iL)
      IMPLICIT NONE
      character*(*) Text
      integer iLevel,iL

      common /DataGJVerror/ UseWithPYTHIA, ListFinalEvent
      integer UseWithPYTHIA, ListFinalEvent
      save /DataGJVerror/

      write(*,1001) iLevel,iL,Text
      if (UseWithPYTHIA.eq.1) then
         call PYLIST(1)
      else
         call LULIST(1)
      endif
      stop

 1001 FORMAT('GetJetsetVec:ERROR[',i1,'] (L',i2,'): ',A,' STOP!')

      end

c-----------------------------------------------------------------
c
c emit a 'GetJetsetVec:Mesage'--Message.
c
c (remember event with ListLevel=1 for a final GJV_ListFinalEvent).
c
      subroutine GJV_MSG(Text,iLevel,iL)
      IMPLICIT NONE
      character*(*) Text
      integer iLevel,iL

      common /DataGJVerror/ UseWithPYTHIA, ListFinalEvent
      integer UseWithPYTHIA, ListFinalEvent
      save /DataGJVerror/

      write(*,1001) iLevel,iL,Text
      ListFinalEvent = 1

 1001 FORMAT('GetJetsetVec:Message[',i1,'] (L',i2,'): ',A)

      end

c-----------------------------------------------------------------
c
c emit a 'GetJetsetVec:Mesage'--Message and list the event
c
c (remember event with ListLevel=2 for a final GJV_ListFinalEvent).
c
      subroutine GJV_MSG_L(Text,iLevel,iL)
      IMPLICIT NONE
      character*(*) Text
      integer iLevel,iL

      common /DataGJVerror/ UseWithPYTHIA, ListFinalEvent
      integer UseWithPYTHIA, ListFinalEvent
      save /DataGJVerror/

      write(*,1001) iLevel,iL,Text

      if (UseWithPYTHIA.eq.1) then
         call PYLIST(1)
      else
         call LULIST(1)
      endif
      ListFinalEvent = 2

 1001 FORMAT('GetJetsetVec:Message[',i1,'] (L',i2,'): ',A)

      end

c-----------------------------------------------------------------
c
c List the final event according previous GJV_MSG* calls
c
c A call with 'ListLevel'==0 initializes this code.
c A call with 'ListLevel'==1 lists the current event (after a call to
c     LU/PYEDIT), if there has been a call to GJV_MSG(...) inbetween.
c A call with 'ListLevel'==2 lists the current event (after a call to
c     LU/PYEDIT), if there has been a call to GJV_MSG_L(...) inbetween.
c
      subroutine GJV_ListFinalEvent(ListLevel)
      IMPLICIT NONE
      integer ListLevel

      common /DataGJVerror/ UseWithPYTHIA, ListFinalEvent
      integer UseWithPYTHIA, ListFinalEvent
      save /DataGJVerror/

      if (ListLevel.eq.0) then ! Initialisation
         ListFinalEvent = 0
      else
         if (ListLevel.le.ListFinalEvent) then
            if (UseWithPYTHIA.eq.1) then
               call PYLIST(1)
            else
               call LULIST(1)
            endif  
         endif
         ListFinalEvent = 0
      endif
      end

c=================================================================
c
c Here we add some dummies...
c
c=================================================================

c=================================================================

c===================================================================
