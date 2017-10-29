!******************************************************************************
!****m* /ID_translation
! NAME
! module ID_translation
! PURPOSE
! Translation of ID numbers between BUU-codes and JETSET-codes.
! NOTES
! These routines are collected from the original files
! "fritzi.f" and "fritiof.f".
!******************************************************************************
module ID_translation

  implicit none
  private

  integer, parameter :: zero156 (156) = 0

  integer, parameter :: kfbuu (122,-1:2) = reshape ( &   ! array of KF-codes for BUU |ID| and IZ
    !--- BARYONS:
    (/     0,  2112,  2212,     0, &  !  1: N(938)
        1114,  2114,  2214,  2224, &  !  2: P33(1232) ****
           0,  2112,  2212,     0, &  !  3: P11(1440) ****
           0,  2112,  2212,     0, &  !  4: S11(1535) ***
           0,  2112,  2212,     0, &  !  5: S11(1650) ****
           0,  2112,  2212,     0, &  !  6: S11(2090) *
           0,  2112,  2212,     0, &  !  7: D13(1520) ****
           0,  2112,  2212,     0, &  !  8: D13(1700) *
           0,  2112,  2212,     0, &  !  9: D13(2080) *
           0,  2112,  2212,     0, &  ! 10: D15(1675) ****
           0,  2112,  2212,     0, &  ! 11: G17(2190) ****
           0,  2112,  2212,     0, &  ! 12: P11(1710) *
           0,  2112,  2212,     0, &  ! 13: P11(2100) *
           0,  2112,  2212,     0, &  ! 14: P13(1720) *
           0,  2112,  2212,     0, &  ! 15: P13(1900) ***
           0,  2112,  2212,     0, &  ! 16: F15(1680) ****
           0,  2112,  2212,     0, &  ! 17: F15(2000) *
           0,  2112,  2212,     0, &  ! 18: F17(1990) **
        1114,  2114,  2214,  2224, &  ! 19: S31(1620) **
        1114,  2114,  2214,  2224, &  ! 20: S31(1900) ***
        1114,  2114,  2214,  2224, &  ! 21: D33(1700) *
        1114,  2114,  2214,  2224, &  ! 22: D33(1940) *
        1114,  2114,  2214,  2224, &  ! 23: D35(1930) **
        1114,  2114,  2214,  2224, &  ! 24: D35(2350) **
        1114,  2114,  2214,  2224, &  ! 25: P31(1750) *
        1114,  2114,  2214,  2224, &  ! 26: P31(1910) ****
        1114,  2114,  2214,  2224, &  ! 27: P33(1600) ***
        1114,  2114,  2214,  2224, &  ! 28: P33(1920) *
        1114,  2114,  2214,  2224, &  ! 29: F35(1750) *
        1114,  2114,  2214,  2224, &  ! 30: F35(1905) ***
        1114,  2114,  2214,  2224, &  ! 31: F37(1950) ****
    !---
           0,  3122,     0,     0, &  ! 32: Lambda(1116)
        3112,  3212,  3222,     0, &  ! 33: Sigma(1189)
        3114,  3214,  3224,     0, &  ! 34: Sigma P13(1385)
           0,  3122,     0,     0, &  ! 35: Lambda S01(1405)
           0,  3122,     0,     0, &  ! 36: Lambda D03(1520)
           0,  3122,     0,     0, &  ! 37: Lambda P01(1600)
           0,  3122,     0,     0, &  ! 38: Lambda P11(1660)
           0,  3122,     0,     0, &  ! 39: Lambda S01(1670)
           0,  3122,     0,     0, &  ! 40: Lambda D13(1670)
           0,  3122,     0,     0, &  ! 41: Lambda D03(1690)
           0,  3122,     0,     0, &  ! 42: Lambda S11(1750)
        3114,  3214,  3224,     0, &  ! 43: Sigma D15(1775)
        3114,  3214,  3224,     0, &  ! 44: Sigma S01(1800)
        3114,  3214,  3224,     0, &  ! 45: Sigma P01(1810)
           0,  3122,     0,     0, &  ! 46: Lambda F05(1820)
           0,  3122,     0,     0, &  ! 47: Lambda D05(1830)
           0,  3122,     0,     0, &  ! 48: Lambda P03(1890)
           0,  3122,     0,     0, &  ! 49: Lambda F15(1915)
        3114,  3214,  3224,     0, &  ! 50: Sigma F17(2030)
        3114,  3214,  3224,     0, &  ! 51: Sigma G07(2100)
        3114,  3214,  3224,     0, &  ! 52: Sigma F05(2110)
    !---
        3312,  3322,     0,     0, &  ! 53: Xi
        3314,  3324,     0,     0, &  ! 54: Xi*
    !---
        3334,     0,     0,     0, &  ! 55: Omega
    !---
           0,     0,  4122,     0, &  ! 56: Lambda_c
           0,  4112,  4212,  4222, &  ! 57: Sigma_c
           0,  4114,  4214,  4224, &  ! 58: Sigma*_c
           0,  4132,  4232,     0, &  ! 59: Xi_c
           0,  4314,  4324,     0, &  ! 60: Xi*_c
           0,  4332,     0,     0, &  ! 61: Omega_c
    !--- EMPTY ENTRIES:
           zero156, &
    !--- MESONS:
        -211,   111,   211,     0, &  ! 101: pi
           0,   221,     0,     0, &  ! 102: eta
        -213,   113,   213,     0, &  ! 103: rho
           0, 60221,     0,     0, &  ! 104: sigma
           0,   223,     0,     0, &  ! 105: omega
           0,   331,     0,     0, &  ! 106: eta prime
           0,   333,     0,     0, &  ! 107: phi
           0,   441,     0,     0, &  ! 108: eta_c
           0,   443,     0,     0, &  ! 109: J/psi
           0,   311,   321,     0, &  ! 110: K+,K0
        -321,  -311,     0,     0, &  ! 111: K-,K~0
           0,   313,   323,     0, &  ! 112: K*
        -323,  -313,     0,     0, &  ! 113: K~*
           0,   421,   411,     0, &  ! 114: D
        -411,  -421,     0,     0, &  ! 115: D~
           0,   423,   413,     0, &  ! 116: D~
        -413,  -423,     0,     0, &  ! 117: D~*
           0,     0,   431,     0, &  ! 118: D_s+
        -431,     0,     0,     0, &  ! 119: D_s-
           0,     0,   433,     0, &  ! 120: D_s*+
        -433,     0,     0,     0, &  ! 121: D_s*-
           0,   225,     0,     0/) & ! 122: f2(1270)
    ,(/122,4/), order=(/2,1/) )


  integer, parameter :: kfbuuC (36,-1:3) = reshape ( &         ! compressed for reverse search
    !--- BARYONS: (14 Entries)
    (/     0,  2112,  2212,     0,   1, &  !  N(938)
        1114,  2114,  2214,  2224,   2, &  !  P33(1232)
    !---
           0,  3122,     0,     0,  32, &  !  Lambda(1116)
        3112,  3212,  3222,     0,  33, &  !  Sigma(1189)
        3114,  3214,  3224,     0,  34, &  !  Sigma P13(1385)
    !---
        3312,  3322,     0,     0,  53, &  !  Xi
        3314,  3324,     0,     0,  54, &  !  Xi*
    !---
        3334,     0,     0,     0,  55, &  !  Omega
    !---
           0,     0,  4122,     0,  56, &  !  Lambda_c
           0,  4112,  4212,  4222,  57, &  !  Sigma_c
           0,  4114,  4214,  4224,  58, &  !  Sigma*_c
           0,  4132,  4232,     0,  59, &  !  Xi_c
           0,  4314,  4324,     0,  60, &  !  Xi*_c
           0,  4332,     0,     0,  61, &  !  Omega_c
    !--- MESONS: (19 Entries)
        -211,   111,   211,     0,  101, &  !  pi
           0,   221,     0,     0,  102, &  !  eta
        -213,   113,   213,     0,  103, &  !  rho
           0, 60221,     0,     0,  104, &  !  sigma
           0,   223,     0,     0,  105, &  !  omega
           0,   331,     0,     0,  106, &  !  eta prime
           0,   333,     0,     0,  107, &  !  phi
           0,   441,     0,     0,  108, &  !  eta_c
           0,   443,     0,     0,  109, &  !  J/psi
           0,   311,   321,     0,  110, &  !  K+,K0
        -321,  -311,     0,     0,  111, &  !  K-,K~0
           0,   313,   323,     0,  112, &  !  K*
        -323,  -313,     0,     0,  113, &  !  K~*
           0,   421,   411,     0,  114, &  !  D
        -411,  -421,     0,     0,  115, &  !  D~
           0,   423,   413,     0,  116, &  !  D~
        -413,  -423,     0,     0,  117, &  !  D~*
           0,     0,   431,     0,  118, &  !  D_s+
        -431,     0,     0,     0,  119, &  !  D_s-
           0,     0,   433,     0,  120, &  !  D_s*+
        -433,     0,     0,     0,  121, &  !  D_s*-
           0,   225,     0,     0,  122/) & !  f2(1270)
    ,(/36,5/), order=(/2,1/) )


  type tMassEntry
     integer :: KF, KC_FR,KC_PY
     real    :: massGiBUU
  end type tMassEntry

  integer, save :: nMassEntry = 0
  integer, parameter :: nMassEntryMax = 100
  type(tMassEntry), dimension(nMassEntryMax) :: MassEntry

  interface KFfromBUU
    module procedure KFfromBUU1, KFfromBUU2
  end interface

  public :: KFfromBUU, KFtoBUU, BUUKFDeltaQ
  public :: SetBruteForceMasses_PY, SetBruteForceMasses_FR

  public :: SplitMeson, SplitDiquark, SplitBaryon, SplitKFtoQ

contains

  !****************************************************************************
  !****f* ID_translation/KFfromBUU
  ! NAME
  ! function KFfromBUU1 (ID, IZ) result (KF)
  ! PURPOSE
  ! translate from BUU code (ID,IZ) to Jetset code (KF)
  !****************************************************************************
  pure function KFfromBUU1 (ID, IZ) result (KF)
    integer, intent(in) :: ID, IZ
    integer :: KF

    if (ID>0) then
      KF = kfbuu(ID,IZ)
    else
      KF = -kfbuu(-ID,-IZ)
    end if
  end function KFfromBUU1


  !****************************************************************************
  !****f* ID_translation/KFfromBUU2
  ! NAME
  ! function KFfromBUU2 (part) result (KF)
  ! PURPOSE
  ! translate from BUU code to Jetset code (KF)
  !****************************************************************************
  pure function KFfromBUU2 (part) result (KF)
    use particleDefinition
    use IDTable, only: isHadron, photon
    type(particle), intent(in) :: part
    integer :: KF

    if (isHadron(part%ID)) then
      if (part%Antiparticle) then
        KF = -kfbuu(part%ID,-part%charge)
      else
        KF = kfbuu(part%ID,part%charge)
      end if
    else if (part%ID==photon) then
      KF = 22
    else
      KF = -99999
    end if
  end function KFfromBUU2


  !****************************************************************************
  !****s* ID_translation/KFtoBUU
  ! NAME
  ! subroutine KFtoBUU (KF, ID, IZ)
  ! PURPOSE
  ! Translate from Jetset code (KF) to BUU code (ID,IZ).
  !****************************************************************************
  pure subroutine KFtoBUU (KF, ID, IZ)
    use idTable, only: photon, electron
    integer, intent(in) :: KF
    integer, intent(out) :: ID, IZ

    integer :: i,j,hKF,aKF

    hKF = KF ! we need this due to 'intent(in)'
    aKF = abs(hKF)

    if (akF==22) then

      ID=photon
      IZ=0
      return

    else if (aKF==11) then

      ID = electron
      IZ = -1
      if (hKF<0) then
         ID = -ID
         IZ = -IZ
      end if
      return

    else if (aKF<100) then      ! LEPTONS or something else

!      write(*,*) '------ Warning: KFtoBUU : notFound(lep):', hKF

    else if (aKF<1000) then ! MESONS
      if ((hKF==130).or.(hKF==310)) then
         ! convert K0_L / K0_s to K0
         hKF = 311
         aKF = 311
      end if

      do i=15,35
         do j=-1,1
            if (hKF==kfbuuC(i,j)) then
               ID = kfbuuC(i,3)
               IZ = j
               return
            end if
         end do
      end do

!      write(*,*) '------ Warning: KFtoBUU : notFound(MES):', hKF

    else if (aKF<10000) then   ! BARYONS

      do i=1,14
         do j=-1,2
            if (aKF==kfbuuC(i,j)) then
               ID = kfbuuC(i,3)
               IZ = j
               if (hKF<0) then
                  ID = -ID
                  IZ = -IZ
               end if
               return
            end if
         end do
      end do

!      write(*,*) '------ Warning: KFtoBUU : notFound(BAR):', hKF

    else

!      write(*,*) '------ Warning: KFtoBUU : notFound(???):', hKF

    end if

    ID = 0
    IZ = 0

  end subroutine KFtoBUU


  !****************************************************************************
  !****f* ID_translation/BUUKFDeltaQ
  ! NAME
  ! integer function BUUKFDeltaQ (deltaQ, ID, IZ) result(iQ)
  ! PURPOSE
  ! return the maximal possible charge correction (iQ) into the direction
  ! (deltaQ) for particle (ID,IZ)
  !****************************************************************************
  integer function BUUKFDeltaQ (deltaQ, ID, IZ) result(iQ)
    use idTable, only: photon, isBaryon, isMeson
    use CallStack, only: traceback
    integer, intent(in) :: deltaQ, ID, IZ
    integer, save :: kfbuuDeltaQ(122,-1:2,2)  ! tabulated values
    logical :: first = .true.

    if (first) call init

    if (isBaryon(abs(ID))) then  ! BARYONS
      if (deltaQ<0) then
        if (ID>0) then
          iQ = kfbuuDeltaQ(ID,IZ,1)
        else
          iQ = kfbuuDeltaQ(-ID,-IZ,2)
        end if
      else
        if (ID>0) then
          iQ = kfbuuDeltaQ(ID,IZ,2)
        else
          iQ = kfbuuDeltaQ(-ID,-IZ,1)
        end if
      end if
    else if (isMeson(abs(ID))) then        ! MESONS
      if (deltaQ<0) then
        iQ = kfbuuDeltaQ(ID,IZ,1)
      else
        iQ = kfbuuDeltaQ(ID,IZ,2)
      end if
    else if (ID==photon) then
      iQ = 0
    else
      write(*,*) "funny ID in BUUKFDeltaQ: ", ID, IZ, deltaQ
      call traceback("error in BUUKFDeltaQ")
    end if

  contains

    subroutine init
      integer :: i, j, j1
      ! Calculate maximal possible charge corretion array
      do i=1,122
        do j=-1,2
          kfbuuDeltaQ(i,j,1:2) = 0
          if (kfbuu(i,j)/=0) then
            do j1=-1,j-1
              if (kfbuu(i,j1)/=0) kfbuuDeltaQ(i,j,1) = kfbuuDeltaQ(i,j,1) - 1
            end do
            do j1=j+1,2
              if (kfbuu(i,j1)/=0) kfbuuDeltaQ(i,j,2) = kfbuuDeltaQ(i,j,2) + 1
            end do
          end if
        end do ! j
      end do ! i
      first = .false.
    end subroutine init

  end function BUUKFDeltaQ


  !****************************************************************************
  !****s* ID_translation/SetBruteForceMasses_INIT
  ! NAME
  ! subroutine SetBruteForceMasses_INIT
  !
  ! PURPOSE
  ! initialize the "SetBruteForceMasses" routines
  !****************************************************************************
  subroutine SetBruteForceMasses_INIT
    use ParticleProperties, only: hadron

    integer LUCOMP ! prototype
    integer PYCOMP ! prototype

    integer :: i,j,aKF
    integer :: iMassEntry

    if (nMassEntry>0) return

    do i=1,14  ! --- BARYON ---
       jLoop1: do j=-1,2
          aKF = kfbuuC(i,j)
          if (aKF.eq.0) cycle jLoop1
          do iMassEntry=1,nMassEntry
             if (MassEntry(iMassEntry)%KF .eq. aKF) cycle jLoop1
          end do

          nMassEntry = nMassEntry+1
          MassEntry(nMassEntry)%KF = aKF
          MassEntry(nMassEntry)%KC_FR = LUCOMP(aKF)
          MassEntry(nMassEntry)%KC_PY = PYCOMP(aKF)
          MassEntry(nMassEntry)%massGiBUU = hadron(kfbuuC(i,3))%mass

       end do jLoop1
    end do
    do i=15,36 ! --- MESON ---
       jLoop2: do j=-1,1
          aKF = abs(kfbuuC(i,j))
          if (aKF.eq.0) cycle jLoop2
          if (aKF.gt.10000) cycle jLoop2
          do iMassEntry=1,nMassEntry
             if (MassEntry(iMassEntry)%KF .eq. aKF) cycle jLoop2
          end do

          nMassEntry = nMassEntry+1
          MassEntry(nMassEntry)%KF = aKF
          MassEntry(nMassEntry)%KC_FR = LUCOMP(aKF)
          MassEntry(nMassEntry)%KC_PY = PYCOMP(aKF)
          MassEntry(nMassEntry)%massGiBUU = hadron(kfbuuC(i,3))%mass

       end do jLoop2
    end do


!!$    write(*,*) '------------------------------------'
!!$    write(*,*) 'SetBruteForceMasses_INIT:'
!!$    write(*,*) '  i    aKF  KC_FR  KC_PY        mass'
!!$    do i=1,nMassEntry
!!$       write(*,'(i4,3i7,f12.3)') i,&
!!$            &MassEntry(i)%KF,MassEntry(i)%KC_FR,MassEntry(i)%KC_PY,&
!!$            &MassEntry(i)%massGiBUU
!!$    enddo
!!$    write(*,*) '------------------------------------'

  end subroutine SetBruteForceMasses_INIT


  !****************************************************************************
  !****s* ID_translation/SetBruteForceMasses_PY
  ! NAME
  ! subroutine SetBruteForceMasses_PY
  !
  ! PURPOSE
  ! set all PYTHIA particle masses to the GiBUU value
  !
  !****************************************************************************
  subroutine SetBruteForceMasses_PY

    COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
    integer KCHG
    double precision PMAS,PARF,VCKM
    SAVE /PYDAT2/

    integer :: i

    if (nMassEntry.eq.0) call SetBruteForceMasses_INIT

    do i=1,nMassEntry
       PMAS(MassEntry(i)%KC_PY,1)=MassEntry(i)%massGiBUU
    end do

  end subroutine SetBruteForceMasses_PY


  !****************************************************************************
  !****s* ID_translation/SetBruteForceMasses_FR
  ! NAME
  ! subroutine SetBruteForceMasses_FR
  !
  ! PURPOSE
  ! set all FRITIOF particle masses to the GiBUU value
  !
  !****************************************************************************
  subroutine SetBruteForceMasses_FR

    COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
    integer KCHG
    real PMAS,PARF,VCKM
    SAVE /LUDAT2/

    integer :: i

    if (nMassEntry.eq.0) call SetBruteForceMasses_INIT

    do i=1,nMassEntry
       PMAS(MassEntry(i)%KC_FR,1)=MassEntry(i)%massGiBUU
    end do

  end subroutine SetBruteForceMasses_FR


  !****************************************************************************
  !****s* ID_translation/SplitMeson(KF, Q1, Q2)
  ! NAME
  ! subroutine SplitMeson(KF, Q1, Q2)
  ! PURPOSE
  ! Split the KF code of a meson into its quark content
  !****************************************************************************
  subroutine SplitMeson(KF, Q1, Q2)

    integer, intent(in) :: KF
    integer, intent(out):: Q1, Q2

    Q1 = abs(KF)/100
    Q2 = -(abs(KF)-100*Q1)/10

    if ((Q1.eq.1).or.(Q1.eq.3).or.(Q1.eq.5)) then
       Q1 = -Q1
       Q2 = -Q2
    end if

    if (KF.lt.0) then
       Q1 = -Q1
       Q2 = -Q2
    end if

  end subroutine SplitMeson

  !****************************************************************************
  !****s* ID_translation/SplitDiquark(KF, Q1, Q2)
  ! NAME
  ! subroutine SplitDiquark(KF, Q1, Q2)
  ! PURPOSE
  ! Split the KF code of a diquark into its quark content
  !****************************************************************************
  subroutine SplitDiquark(KF, Q1, Q2)

    integer, intent(in) :: KF
    integer, intent(out):: Q1, Q2

    Q1 = abs(KF)/1000
    Q2 = (abs(KF) - 1000*Q1)/100
    if (KF.lt.0) then
       Q1 = -Q1
       Q2 = -Q2
    end if

  end subroutine SplitDiquark

  !****************************************************************************
  !****s* ID_translation/SplitBaryon(KF, Q1, Q2, Q3)
  ! NAME
  ! subroutine SplitBaryon(KF, Q1, Q2, Q3)
  ! PURPOSE
  ! Split the KF code of a baryon into its quark content
  !****************************************************************************
  subroutine SplitBaryon(KF, Q1, Q2, Q3)

    integer, intent(in) :: KF
    integer, intent(out):: Q1, Q2, Q3

    Q1 = abs(KF)/1000
    Q2 = (abs(KF)-1000*Q1)/100
    Q3 = (abs(KF)-1000*Q1-100*Q2)/10

    if (KF.lt.0) then
       Q1 = -Q1
       Q2 = -Q2
       Q3 = -Q3
    end if
  end subroutine SplitBaryon

  !****************************************************************************
  !****s* ID_translation/SplitKFtoQ(KF, Q1, Q2, Q3)
  ! NAME
  ! subroutine SplitKFtoQ(KF, Q1, Q2, Q3)
  ! PURPOSE
  ! Split the KF code of any particle into its quark content
  !****************************************************************************
  subroutine SplitKFtoQ(KF, Q1, Q2, Q3)

    integer, intent(in) :: KF
    integer, intent(out):: Q1, Q2, Q3

    integer :: aKF

    aKF = abs(KF)

    Q1 = 0
    Q2 = 0
    Q3 = 0

    if (aKF < 10) then
       Q3 = KF
    else if (aKF < 100) then
       !
    else if (aKF < 1000) then
       call SplitMeson(KF, Q2, Q3)
    else if (aKF < 10000) then
       call SplitBaryon(KF, Q1, Q2, Q3)

       if (Q3 == 0) then ! it was a diquark
          Q3 = Q2
          Q2 = Q1
          Q1 = 0
       end if

    else
       write(*,*) 'SplitKFtoQ: KF = ',KF
    end if

  end subroutine SplitKFtoQ


end module ID_translation
