!******************************************************************************
!****m* /CollGetLeading
! PURPOSE
! Recontruct the number of leading quraks for a JetSet fragmentation
!******************************************************************************
module CollGetLeading

  IMPLICIT NONE
  private

  public :: GetLeading_PY, GetLeading_FR


  logical :: UseWithPYTHIA = .true.

contains

  !****************************************************************************
  !****s* CollGetLeading/GetLeading_PY
  ! NAME
  ! subroutine GetLeading_PY
  !
  ! PURPOSE
  ! Try to reconstruct for every outgoing partice of a Jetset fragmentation
  ! the number of leading quarks of this particle. i.e. the number of
  ! quarks contained in the final partickle, which already spanned the
  ! fragmented string.
  !
  ! INPUTS
  ! The Pythia common blocks have to be set (no PYEDIT!)
  !
  ! OUTPUT
  ! in common block /PYJETS/:
  ! * K(i,4): 0,1,2,3 according number of leading quarks
  ! * K(i,5): number of fragmented string in this event
  !****************************************************************************

  subroutine GetLeading_PY

    !...common blocks:

    COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
    integer N,NPAD,K
    double precision P,V
    SAVE /PYJETS/

    COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
    integer MSTP,MSTI
    double precision PARP,PARI
    SAVE /PYPARS/

    COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
    integer MSTU,MSTJ
    double precision PARU,PARJ
    SAVE /PYDAT1/

    !...overall variables:

    integer iL                ! actual line
    integer pL1, pL2          ! parent lines (see below)
    integer cL1, cL2          ! child lines  (see below)
    integer cL12              ! #(child lines)= cL2-cL1+1

    integer pQ1, pQ2          ! #quarks in parent 1,2

    !...helpful variables:

    integer i,j
    integer iQ1,iQ2,iQ3,iQA1,iQA2
    integer nQ1, nQ2
    integer pL1C              ! saves CODE(pL1) = K(pL1,1)
    integer iS

    UseWithPYTHIA = .TRUE.

    if (MSTU(16).ne.2) call GL_ERR('GetLeading_PY requires MSTU(16)=2 !!!',0,0)


    iL = MSTI(4)
    iS = 0

100 iL = iL+1
    if (K(iL,1).ge.10) goto 999 ! not stable, skip ahead

    pL1 = K(iL,3)
    if (pL1.le.0) then        ! check this for FRITIOF
       pL1C = 21
    else
       pL1C = K(pL1,1)
    end if
    iS = iS+1

    select case (pL1C)

    case (12,51) !...pL1C->beginning of q-string:
       !...============================
       !     ...pL1,pL2 = begin & end of quark string (parent)
       !     ...cL1,cL2 = begin & end of resulting particles (child)

       cL1 = K(pL1,4)
       cL2 = K(pL1,5)
       cL12 = cL2-cL1+1

       if (cL1.ne.iL) call GL_ERR('cL1.ne.iL',0,iL)

       do j=pL1+1,N
          if (K(j,1).eq.11) then
             pL2 = j
             goto 111
          end if
       end do
       call GL_ERR('end of q-string not found',0,iL)
111    continue

!!!         pL2 = K(cL2,3) ???

       !... reserve Entries

       do i=cL1,cL2
          K(i,4) = 0
          K(i,5) = iS
       end do

       select case (cL12) !... check number of child particles

       case (1)  !----- 1 particle, collapse

          pQ1 = CalcParQuarks(K(pL1,2),iQ1,iQ2)
          pQ2 = CalcParQuarks(K(pL2,2),iQA1,iQA2)

          if (pQ1.eq.0) goto 990
          if (pQ2.eq.0) goto 990

          if (pQ1.eq.1) then
             iQ2 = iQA1
             iQ3 = iQA2
          else
             iQ3 = iQA1
          end if
          if (QQQeP(K(cL1,2),iQ1,iQ2,iQ3)) then
             K(cL1,4) = pQ1+pQ2
          else
             call GL_ERR('must fit',1,iL)
          end if

       case (2)  !----- 2 particles

          pQ1 = CalcParQuarks(K(pL1,2),iQ1,iQ2)
          pQ2 = CalcParQuarks(K(pL2,2),iQA1,iQA2)

          if (pQ1.eq.0) goto 990
          if (pQ2.eq.0) goto 990

          if (QQQePP(K(cL1,2),K(cL2,2), pQ1,pQ2, iQ1,iQ2, iQA1,iQA2,nQ1,nQ2)) then
             K(cL1,4) = nQ1
             K(cL2,4) = nQ2
          else
             call GL_ERR('2 particles, urgh!!!',2,iL)
          end if

       case (3:)  !----- 3 or more particles

          call Check3(pL1,cL1, 1)
          call Check3(pL2,cL2,-1)

       end select


    case (21) !...pL1C->from decayed particle:
             !...============================

       if (pL1.le.MSTI(4)) then ! from docu
          K(iL,4) = CountQ(K(iL,2))
       else
          call GL_ERR('where does this come from?',9,iL)
       end if
       cL2 = iL         ! last line worked on


    case default  !...pL1C->something unknown:
                  !...========================

       call GL_ERR('unknown code',0,iL)

    end select

990 continue
    iL = cL2

999 if (iL.lt.N) goto 100  ! loop: find next string


  contains
    !==========================================================================
    ! Check, whether particle in line 'iL' and/or the following/preceeding
    ! lines (depending on 'dir') may contain the quark content given by the
    ! parent line 'parLine'.
    ! This routine treats the case that we have 3 or more output particles.

    subroutine Check3(parLine, iL, dir)

      integer parLine,iL,dir

      !...internal variables:
      integer pQ, iQ1, iQ2
      integer, parameter :: iDiWeight=1  ! splitted diquark counts/counts not



      pQ = CalcParQuarks(K(parLine,2),iQ1,iQ2)
      if (pQ.eq.0) return

      if (pQ.eq.1) then ! quark
         if (QQQeP(K(iL,2),iQ1,0,0)) then
            K(iL,4) = 1
         else
            call GL_ERR('must fit (1)',0,iL)
         end if
      else                      ! diquark
         if (QQQeP(K(iL,2),iQ1,iQ2,0)) then
            K(iL,4) = 2
         else
            !if (K(iL,3).ne.K(iL+dir,3)) then
            !        something strange, but outcome should be correct !
            !                  call GL_MSG_L('parents not equal (1)',0,iL)
            !end if
            if (QQQeP(K(iL,2),iQ1,0,0)) then
               if (QQQeP(K(iL+dir,2),iQ2,0,0)) then
                  K(iL,4)     = iDiWeight
                  K(iL+dir,4) = iDiWeight
               else
                  call GL_ERR('must fit (1a)',0,iL)
               end if
            else if (QQQeP(K(iL,2),iQ2,0,0)) then
               if (QQQeP(K(iL+dir,2),iQ1,0,0)) then
                  K(iL,4)     = iDiWeight
                  K(iL+dir,4) = iDiWeight
               else
                  call GL_ERR('must fit (1b)',0,iL)
               end if
            else
               call GL_ERR('must fit (1c)',0,iL)
            end if
         end if

      end if
    end subroutine Check3

  end subroutine GetLeading_PY

  !****************************************************************************
  ! This is more or less a copy of GetLeading_PY !


  subroutine GetLeading_FR

    !...common blocks:

    COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
    integer N,NPAD,K
    double precision P,V
    SAVE /PYJETS/

    COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
    integer MSTU,MSTJ
    real PARU,PARJ
    SAVE /LUDAT1/

    !...overall variables:

    integer iL                ! actual line
    integer pL1, pL2          ! parent lines (see below)
    integer cL1, cL2          ! child lines  (see below)
    integer cL12              ! #(child lines)= cL2-cL1+1

    integer pQ1, pQ2          ! #quarks in parent 1,2

    !...helpful variables:

    integer i,j
    integer iQ1,iQ2,iQ3,iQA1,iQA2
    integer nQ1, nQ2
    integer pL1C              ! saves CODE(pL1) = K(pL1,1)
    integer iS
    integer, parameter :: MSTI4 = 0

    UseWithPYTHIA = .FALSE.

    if (MSTU(16).ne.2) call GL_ERR('GetLeading_FR requires MSTU(16)=2 !!!',0,0)


    iL = MSTI4
    iS = 0

!    call LULIST(2)


100 iL = iL+1
    if (K(iL,1).ge.10) goto 999 ! not stable, skip ahead

    pL1 = K(iL,3)
    if (pL1.le.0) then        ! check this for FRITIOF
       pL1C = 21
    else
       pL1C = K(pL1,1)
    end if
    iS = iS+1

    select case (pL1C)

    case (12) !...pL1C->beginning of q-string:
       !...============================
       !     ...pL1,pL2 = begin & end of quark string (parent)
       !     ...cL1,cL2 = begin & end of resulting particles (child)

       cL1 = K(pL1,4)
       cL2 = K(pL1,5)
       cL12 = cL2-cL1+1

!       write(*,*) iL, pL1,cL1,cL2

       if (cL1.ne.iL) call GL_ERR('cL1.ne.iL',0,iL)

       do j=pL1+1,N
          if (K(j,1).eq.11) then
             pL2 = j
             goto 111
          end if
       end do
       call GL_ERR('end of q-string not found',0,iL)
111    continue

!!!         pL2 = K(cL2,3) ???

       !... reserve Entries

       do i=cL1,cL2
          K(i,4) = 0
          K(i,5) = iS
       end do

       select case (cL12) !... check number of child particles

       case (1)  !----- 1 particle, collapse

          pQ1 = CalcParQuarks(K(pL1,2),iQ1,iQ2)
          pQ2 = CalcParQuarks(K(pL2,2),iQA1,iQA2)

          if (pQ1.eq.0) goto 990
          if (pQ2.eq.0) goto 990

          if (pQ1.eq.1) then
             iQ2 = iQA1
             iQ3 = iQA2
          else
             iQ3 = iQA1
          end if
          if (QQQeP(K(cL1,2),iQ1,iQ2,iQ3)) then
             K(cL1,4) = pQ1+pQ2
          else
             call GL_ERR('must fit',1,iL)
          end if

       case (2)  !----- 2 particles

          pQ1 = CalcParQuarks(K(pL1,2),iQ1,iQ2)
          pQ2 = CalcParQuarks(K(pL2,2),iQA1,iQA2)

          if (pQ1.eq.0) goto 990
          if (pQ2.eq.0) goto 990

          if (QQQePP(K(cL1,2),K(cL2,2), pQ1,pQ2, iQ1,iQ2, iQA1,iQA2,nQ1,nQ2)) then
             K(cL1,4) = nQ1
             K(cL2,4) = nQ2
          else
             call GL_ERR('2 particles, urgh!!!',2,iL)
          end if

       case (3:)  !----- 3 or more particles

          call Check3(pL1,cL1, 1)
          call Check3(pL2,cL2,-1)

       end select



    case (21) !...pL1C->from decayed particle:
       !...============================

       if (pL1.le.MSTI4) then ! from docu
          K(iL,4) = CountQ(K(iL,2))
       else
          call GL_ERR('where does this come from?',9,iL)
       end if
       cL2 = iL         ! last line worked on


    case default  !...pL1C->something unknown:
       !...========================

       call GL_ERR('unknown code',0,iL)

    end select

990 continue
    iL = cL2

999 if (iL.lt.N) goto 100  ! loop: find next string


  contains
    !==========================================================================
    ! Check, whether particle in line 'iL' and/or the following/preceeding
    ! lines (depending on 'dir') may contain the quark content given by the
    ! parent line 'parLine'.
    ! This routine treats the case that we have 3 or more output particles.

    subroutine Check3(parLine, iL, dir)

      integer parLine,iL,dir

      !...internal variables:
      integer pQ, iQ1, iQ2
      integer, parameter :: iDiWeight=1  ! splitted diquark counts/counts not


      pQ = CalcParQuarks(K(parLine,2),iQ1,iQ2)
      if (pQ.eq.0) return

      if (pQ.eq.1) then ! quark
         if (QQQeP(K(iL,2),iQ1,0,0)) then
            K(iL,4) = 1
         else
            call GL_ERR('must fit (1)',0,iL)
         end if
      else                      ! diquark
         if (QQQeP(K(iL,2),iQ1,iQ2,0)) then
            K(iL,4) = 2
         else
            !if (K(iL,3).ne.K(iL+dir,3)) then
            !        something strange, but outcome should be correct !
            !                  call GL_MSG_L('parents not equal (1)',0,iL)
            !endif
            if (QQQeP(K(iL,2),iQ1,0,0)) then
               if (QQQeP(K(iL+dir,2),iQ2,0,0)) then
                  K(iL,4)     = iDiWeight
                  K(iL+dir,4) = iDiWeight
               else
                  call GL_ERR('must fit (1a)',0,iL)
               end if
            else if (QQQeP(K(iL,2),iQ2,0,0)) then
               if (QQQeP(K(iL+dir,2),iQ1,0,0)) then
                  K(iL,4)     = iDiWeight
                  K(iL+dir,4) = iDiWeight
               else
                  call GL_ERR('must fit (1b)',0,iL)
               end if
            else
               call GL_ERR('must fit (1c)',0,iL)
            end if
         end if

      end if
    end subroutine Check3

  end subroutine GetLeading_FR

  !****************************************************************************
  !****************************************************************************
  !****************************************************************************

  integer function CalcParQuarks(parKF, parQ1, parQ2)

    use ID_translation, only: SplitDiquark

    integer parKF
    integer,intent(out):: parQ1, parQ2

    CalcParQuarks = 0
    if (abs(parKF).lt.10) then
       CalcParQuarks = 1      ! single quark, 1 leading quark
       parQ1 = parKF
       parQ2 = 0
    else if (mod(abs(parKF)/10,10).eq.0) then
       CalcParQuarks = 2      ! diquark, max 2 leading quarks
       call SplitDiquark(parKF,parQ1,parQ2)
    else if (parKF.eq.21) then
       CalcParQuarks = 0      ! GLUON !
       !         call GL_MSG_L('parent is gluon',0,0)
    else if (parKF.eq.88) then
       CalcParQuarks = 0      ! Junction !
    else
       write(*,*) 'CalcParQuarks: unknown quark content of ',parKF,'. STOP.'
       call GL_ERR('',0,0)
    end if
    return
  end function CalcParQuarks

  !============================================================================
  !
  ! neglecting: iKF==0; Reggeon, Pomeron etc.; iKF>10000

  integer function CountQ(iKF)

    integer iKF

    integer aKF
    aKF = abs(iKF)

    if (aKF.lt.10) then
       CountQ = 1             ! quarks
    else if (aKF.lt.100) then
       CountQ = 0             ! leptons etc, special codes
    else if (aKF.lt.1000) then
       CountQ = 2             ! mesons
    else if (aKF.lt.10000) then
       if (mod(aKF/10,10).eq.0) then
          CountQ = 2          ! diquark
       else
          CountQ = 3          ! baryons
       end if
    else
       write(*,*) 'CountQ: not prepared fo KF=|',iKF,'| >10000. STOP.'
       stop
    end if
  end function CountQ

  !============================================================================

  logical function QQQeP(KFP, KFQ1,KFQ2,KFQ3)

    use ID_translation, only: splitmeson, splitbaryon

    integer KFP, KFQ1, KFQ2, KFQ3

    integer aKFP, aKFQ1, aKFQ2, aKFQ3, iQ1,iQ2,iQ3

    if (KFP.lt.0) then
       aKFP = -KFP
       aKFQ1 = -KFQ1
       aKFQ2 = -KFQ2
       aKFQ3 = -KFQ3
    else
       aKFP = KFP
       aKFQ1 = KFQ1
       aKFQ2 = KFQ2
       aKFQ3 = KFQ3
    end if

    QQQeP = .FALSE.

    aKFP = mod(aKFP,10000)

    if (aKFP.ge.1000) goto 100 ! -> Baryons

    !----- Mesons: ---------

    if (aKFQ3.ne.0) return

    select case (aKFP)

    case (111,113,115,223) !... pi0,rho0,omega,115=a_20:

       if (abs(aKFQ1).gt.2) return
       if ((aKFQ2.ne.0).and.(aKFQ1.ne.-aKFQ2)) return


    case (221,331)     !... eta,etaprime:

       if (abs(aKFQ1).gt.3) return
       if ((aKFQ2.ne.0).and.(aKFQ1.ne.-aKFQ2)) return

    case default

       call SplitMeson(aKFP, iQ1, iQ2)
       if (aKFQ1.eq.iQ1) then
          if ((aKFQ2.ne.0).and.(aKFQ2.ne.iQ2)) return
       else if (aKFQ1.eq.iQ2) then
          if ((aKFQ2.ne.0).and.(aKFQ2.ne.iQ1)) return
       else
          return
       end if

    end select


    QQQeP = .TRUE.
    return

    !----- Baryons: ---------
100 continue

    if (aKFQ1.lt.0) return   ! no antiquark in baryon
    if (aKFQ2.lt.0) return   ! no antiquark in baryon
    if (aKFQ3.lt.0) return   ! no antiquark in baryon

    call SplitBaryon(aKFP,iQ1,iQ2,iQ3)

    if (aKFQ1.eq.iQ1) then
       iQ1 = 0
    else if (aKFQ1.eq.iQ2) then
       iQ2 = 0
    else if (aKFQ1.eq.iQ3) then
       iQ3 = 0
    else
       return
    end if

    ! beachte: keine Abfrage (aKFQ2.ne.0) nötig!
    if (aKFQ2.eq.iQ1) then
       iQ1 = 0
    else if (aKFQ2.eq.iQ2) then
       iQ2 = 0
    else if (aKFQ2.eq.iQ3) then
       iQ3 = 0
    else
       return
    end if

    if (aKFQ3.eq.iQ1) then
       iQ1 = 0
    else if (aKFQ3.eq.iQ2) then
       iQ2 = 0
    else if (aKFQ3.eq.iQ3) then
       iQ3 = 0
    else
       return
    end if

    QQQeP = .TRUE.

    return
  end function QQQeP

  !============================================================================
  !
  !
  ! INPUT:
  !     P1,  P2 : KF codes of particle 1, 2
  !     nA,  nB : number of quarks from parent A, B
  !     QA1, QA2: quarks from parent A
  !     QB1, QB2: quarks from parent B
  ! OUTPUT:
  !     n1,  n2 : number of leading quarks in particle 1, 2
  !
  ! assumed: nA,nB > 0; ...

  logical function QQQePP(P1,P2, nA,nB, QA1,QA2, QB1,QB2, n1,n2)

    integer P1,P2, nA,nB, QA1,QA2, QB1,QB2
    integer,intent(out):: n1,n2

    integer n, Q1, Q2, Q3, Q4

    n = nA+nB

    !$$$      write (*,*) 'calling QQQePP: ',
    !$$$     $     'P1,P2:',P1,P2,'  nA,nB:',nA,nB,
    !$$$     $     '  QA1,QA2, QB1,QB2:',QA1,QA2, QB1,QB2

    select case (n)

    case (2)

       n1 = 1
       n2 = 1
       Q1 = QA1
       Q2 = QB1

       if ((QQQeP(P1,Q1,0,0).and.QQQeP(P2,Q2,0,0)).or. &
            (QQQeP(P1,Q2,0,0).and.QQQeP(P2,Q1,0,0))) then
          ! nothing to do
       else
          write(*,*) 'QQQePP: 2: more complicated. STOP.'
          stop
       end if

    case (3)

       Q1 = QA1
       if (nA.eq.1) then      ! nA = 1, nB = 2
          Q2 = QB1
          Q3 = QB2
       else                   ! nA = 2, nB = 1
          Q2 = QA2
          Q3 = QB1
       end if

       if ( (QQQeP(P1,Q1,0,0).and.QQQeP(P2,Q2,Q3,0)).or. &
            (QQQeP(P1,Q2,0,0).and.QQQeP(P2,Q3,Q1,0)).or. &
            (QQQeP(P1,Q3,0,0).and.QQQeP(P2,Q1,Q2,0))) then
          n1=1
          n2=2
       else if ((QQQeP(P2,Q1,0,0).and.QQQeP(P1,Q2,Q3,0)).or. &
            (QQQeP(P2,Q2,0,0).and.QQQeP(P1,Q3,Q1,0)).or. &
            (QQQeP(P2,Q3,0,0).and.QQQeP(P1,Q1,Q2,0))) then
          n1=2
          n2=1
       else
          write(*,*) 'QQQePP: 3(a): more complicated. STOP.'
          write(*,*) ': P1,P2: ',P1,P2
          write(*,*) ': QA1,QA2, QB1,QB2: ',QA1,QA2, QB1,QB2
          write(*,*) ': Q1, Q2, Q3: ', Q1, Q2, Q3
          stop
       end if

    case (4)

!!$       write(*,*) 'QQQePP: nA+nB = 4'
!!$       write (*,*)     'P1,P2:',P1,P2,'  nA,nB:',nA,nB
!!$       write (*,*)     '  QA1,QA2, QB1,QB2:',QA1,QA2, QB1,QB2
!!$       call PYLIST(2)

       if (nA.eq.2) then
          Q1 = QA1
          Q2 = QA2
          Q3 = QB1
          Q4 = QB2
       else
          write(*,*) 'QQQePP: nA != 2. STOP'
          stop
       end if

       if ( (QQQeP(P1,Q1,0,0).and.QQQeP(P2,Q2,Q3,Q4)).or. &
            (QQQeP(P1,Q2,0,0).and.QQQeP(P2,Q1,Q3,Q4)).or. &
            (QQQeP(P1,Q3,0,0).and.QQQeP(P2,Q1,Q2,Q4)).or. &
            (QQQeP(P1,Q4,0,0).and.QQQeP(P2,Q1,Q2,Q3)) ) then
          n1 = 1
          n2 = 3
       else if ( (QQQeP(P1,Q1,Q2,0).and.QQQeP(P2,Q3,Q4,0)).or. &
            (QQQeP(P1,Q1,Q3,0).and.QQQeP(P2,Q2,Q4,0)).or. &
            (QQQeP(P1,Q1,Q4,0).and.QQQeP(P2,Q2,Q3,0)).or. &
            (QQQeP(P1,Q2,Q3,0).and.QQQeP(P2,Q1,Q4,0)).or. &
            (QQQeP(P1,Q2,Q4,0).and.QQQeP(P2,Q1,Q3,0)).or. &
            (QQQeP(P1,Q3,Q4,0).and.QQQeP(P2,Q1,Q2,0)) ) then
          n1 = 2
          n2 = 2
       else if ( (QQQeP(P2,Q1,0,0).and.QQQeP(P1,Q2,Q3,Q4)).or. &
            (QQQeP(P2,Q2,0,0).and.QQQeP(P1,Q1,Q3,Q4)).or. &
            (QQQeP(P2,Q3,0,0).and.QQQeP(P1,Q1,Q2,Q4)).or. &
            (QQQeP(P2,Q4,0,0).and.QQQeP(P1,Q1,Q2,Q3)) ) then
          n2 = 1
          n1 = 3
       else
          write(*,*) 'QQQePP: 4(a): more complicated. STOP.'
          write(*,*) ': P1,P2: ',P1,P2
          write(*,*) ': QA1,QA2, QB1,QB2: ',QA1,QA2, QB1,QB2
          write(*,*) ': Q1, Q2, Q3, Q4: ', Q1, Q2, Q3, Q4
          stop
       end if

!!$       write(*,*) 'n1,n2:',n1,n2
!!$       stop

    case default

       write(*,*) 'QQQePP: nA+nB != 2,3,4. STOP'

       write(*,*) 'P1,P2:',P1,P2,'  nA,nB:',nA,nB
       write(*,*) '  QA1,QA2, QB1,QB2:',QA1,QA2, QB1,QB2

       call PYLIST(2)

       stop
    end select

    QQQePP = .TRUE.
    return
  end function QQQePP


  !============================================================================
  !
  ! emit a 'GetLeading:Error'--Message, list the event and stop.

  subroutine GL_ERR(Text,iLevel,iL)

    character*(*) Text
    integer iLevel,iL

    write(*,1001) iLevel,iL,Text
    if (UseWithPYTHIA) then
       call PYLIST(2)
    else
       call LULIST(2)
    end if
    stop

1001 FORMAT('GetLeading:ERROR[',i1,'] (L',i2,'): ',A,' STOP!')

  end subroutine GL_ERR

  !-----------------------------------------------------------------
  !
  ! emit a 'GetLeading:Mesage'--Message and list the event

!   subroutine GL_MSG_L(Text,iLevel,iL)
!
!     character*(*) Text
!     integer iLevel,iL
!
!     write(*,1001) iLevel,iL,Text
!
!     if (UseWithPYTHIA) then
!        call PYLIST(1)
!     else
!        call LULIST(1)
!     endif
!
! 1001 FORMAT('GetLeading:Message[',i1,'] (L',i2,'): ',A)
!
!   end subroutine GL_MSG_L

  !============================================================================




end module CollGetLeading
