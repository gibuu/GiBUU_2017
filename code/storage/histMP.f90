!******************************************************************************
!****m* /histMPf90
! NAME
! module histMPf90
! PURPOSE
! Encapsulate all routines and datas for multi-particle Histogramms.
!
! Provides all features as histf90, but here the information is
! stored for many kinds of particles, as e.g. all pions, kaons etc
! at the same time.
! Output is in a gnuplot like multicolumn format
!
! Features of Histograms provided by this module:
! - store paramaters of the x-binning
! - enable two y-values (y and y2)
! - track keeping of under-/over-score the given extreme values of x.
! - provide simple-to-understand output routines
! - provide simple histogram arithmetic
! - many others...
!
! NOTES
! Programming is fairly similar to "module histf90"
!
! There are different sets of particles prepared:
! * 1: pi, K, N
! * 2: pi, rho, (div. mes.), K, N, Delta
! * 3: pi, eta, omega, phi, J/psi, K, D
! * 4: "all mesons"
! * 5: "all stable hadrons": pi,K,N,Lambda,Sigma,Xi,Omega
!
! INPUTS
! ...(This module needs no input)
!******************************************************************************
module histMPf90

  use histf90

  implicit none

  public :: histogramMP

  !****************************************************************************
  !****t* histMPf90/histogramMP
  ! NAME
  ! type histogramMP
  ! PURPOSE
  ! Type definition to store all information for a multiparticle 1D Histogram.
  ! SOURCE
  !
  type histogramMP
     real             :: xMin        ! smallest x-value
     real             :: xMax        ! largest x-value
     real             :: xBin        ! width of x-Bin
     real,allocatable :: xExtreme(:,:) ! extremes of used x-values
                                     ! (id), (min,max)
     character*(NameLength) :: Name  ! name to be written
     real,allocatable :: yVal(:,:,:) ! histogramm values:
                                     ! (id), (x), (yy1,yy2,yy3...)
                                     ! x-bin 0,-1 : extreme values !!!
     integer          :: iSet        ! particle set
  end type histogramMP
  !****************************************************************************

  private

  integer, parameter :: nSet = 5

  logical, save :: initFlag=.true.

  integer, dimension(nSet,22,-1:1), save :: Map2Hist_Meson
  integer, dimension(nSet,61,-1:2), save :: Map2Hist_Baryon
  integer, dimension(nSet,61,-2:1), save :: Map2Hist_ABaryon

  integer, dimension(nSet), save :: Map2Hist_N

  character*12,dimension(nSet,50), save :: Map2Hist_TIT

  public :: histMPf90_Init
  public :: CreateHistMP
  public :: AddHistMP
  public :: WriteHistMP, WriteHistMP_Names
  public :: ClearHistMP

  public :: FetchHistMP

  public :: Map2HistMP
  public :: Map2HistMP_getN

  !****************************************************************************
  !****s* histMPf90/AddHistMP
  ! NAME
  ! subroutine AddHistMP(H, ID,IC,anti, x,y,y2)
  !
  ! subroutine AddHistMP(H, Part, x,y,y2)
  !
  ! subroutine AddHistMP(H, KF, x,y,y2)
  ! PURPOSE
  ! Add to the given histogram at the given x-value the weight y (y2).
  ! If y2 is given, also the second column is filled.
  !
  ! INPUTS
  ! * type(histogramMP) :: H  -- Histogramm to be used
  ! * integer           :: ID,IC -- particle (ID and charge)
  ! * logical           :: anti  -- particle (antiparticle)
  ! * real              :: x  -- x-value
  ! * real              :: y  -- weight to be added
  ! * real              :: y2 -- second weight to be added [OPTIONAL]
  ! or:
  ! * type(histogramMP) :: H  -- Histogramm to be used
  ! * type(particle)    :: Part -- particle
  ! * real              :: x  -- x-value
  ! * real              :: y  -- weight to be added
  ! * real              :: y2 -- second weight to be added [OPTIONAL]
  ! or:
  ! * type(histogramMP) :: H  -- Histogramm to be used
  ! * integer           :: KF -- particle (KF code)
  ! * real              :: x  -- x-value
  ! * real              :: y  -- weight to be added
  ! * real              :: y2 -- second weight to be added [OPTIONAL]
  ! OUTPUT
  ! H is changed
  !****************************************************************************
  interface AddHistMP
     module procedure AddHistMP_I, AddHistMP_P,AddHistMP_KF
  end interface

  !****************************************************************************
  !****f* histMPf90/Map2HistMP
  ! NAME
  ! integer function Map2HistMP(ID,IC,anti, iSet)
  !
  ! integer function Map2HistMP(Part, iSet)
  !
  ! integer function Map2HistMP(KF, iSet)
  !
  ! PURPOSE
  ! calculate the number of the histogram
  ! INPUTS:
  ! * integer           :: ID,IC -- particle (ID and charge)
  ! * logical           :: anti  -- particle (antiparticle)
  ! or:
  ! * type(particle)    :: Part -- particle
  ! or:
  ! * integer           :: KF -- particle (KF code)
  !****************************************************************************
  interface Map2HistMP
     module procedure Map2HistMP_I, Map2HistMP_P, Map2HistMP_KF
  end interface Map2HistMP

  !****************************************************************************
  !****s* histMPf90/WriteHistMP_Names
  ! NAME
  ! subroutine WriteHistMP_Names(H,iFile)
  !
  ! subroutine WriteHistMP_Names(iSet,iFile)
  ! PURPOSE
  ! Write a Header line with the names of the particles
  ! INPUTS
  ! * type(histogramMP) :: H -- histogram to use
  ! * integer :: iFile -- number of file to write to
  ! or:
  ! * integer :: iSet -- set of histogram
  ! * integer :: iFile -- number of file to write to
  !****************************************************************************
  interface WriteHistMP_Names
     module procedure WriteHistMP_Names_H, WriteHistMP_Names_I
  end interface WriteHistMP_Names

contains
  !****************************************************************************
  !****s* histMPf90/histMPf90_Init
  ! NAME
  ! subroutine histMPf90_Init
  ! PURPOSE
  ! Init the module
  !****************************************************************************

  subroutine histMPf90_Init

    integer :: iS,n!,iI,iC

    if (.not.initFlag) return

    ! (1) init maps:

    Map2Hist_Meson = 0
    Map2Hist_Baryon = 0
    Map2Hist_ABaryon = 0

    write(*,*) 'histMPf90_Init'

100 continue ! we are looping twice over these definitions !!!!!

    ! (1a) map of iPartSet==1:

    Map2Hist_TIT(1, 1) = 'pi+'
    Map2Hist_TIT(1, 2) = 'pi0'
    Map2Hist_TIT(1, 3) = 'pi-'
    Map2Hist_TIT(1, 4) = 'K+'
    Map2Hist_TIT(1, 5) = 'K0'
    Map2Hist_TIT(1, 6) = 'K~0'
    Map2Hist_TIT(1, 7) = 'K-'
    Map2Hist_TIT(1, 8) = 'p'
    Map2Hist_TIT(1, 9) = 'n'
    Map2Hist_TIT(1,10) = 'n~'
    Map2Hist_TIT(1,11) = 'p~'

    Map2Hist_Meson(1, 1, 1) =  1 ! pi+
    Map2Hist_Meson(1, 1, 0) =  2 ! pi0
    Map2Hist_Meson(1, 1,-1) =  3 ! pi-

    Map2Hist_Meson(1,10, 1) =  4 ! K+
    Map2Hist_Meson(1,10, 0) =  5 ! K0
    Map2Hist_Meson(1,11, 0) =  6 ! K~0
    Map2Hist_Meson(1,11,-1) =  7 ! K-

    Map2Hist_Baryon(1, 1, 1) =  8 ! p
    Map2Hist_Baryon(1, 1, 0) =  9 ! n

    Map2Hist_ABaryon(1, 1, 0) = 10 ! n~
    Map2Hist_ABaryon(1, 1,-1) = 11 ! p~


    ! (1b) map of iPartSet==2:

    Map2Hist_TIT(2, 1) = 'pi+'
    Map2Hist_TIT(2, 2) = 'pi0'
    Map2Hist_TIT(2, 3) = 'pi-'
    Map2Hist_TIT(2, 4) = 'rho+'
    Map2Hist_TIT(2, 5) = 'rho0'
    Map2Hist_TIT(2, 6) = 'rho-'
    Map2Hist_TIT(2, 7) = 'divMes'
    Map2Hist_TIT(2, 8) = 'K+'
    Map2Hist_TIT(2, 9) = 'K0'
    Map2Hist_TIT(2,10) = 'K~0'
    Map2Hist_TIT(2,11) = 'K-'
    Map2Hist_TIT(2,12) = 'p'
    Map2Hist_TIT(2,13) = 'n'
    Map2Hist_TIT(2,14) = 'n~'
    Map2Hist_TIT(2,15) = 'p~'
    Map2Hist_TIT(2,16) = 'Delta++'
    Map2Hist_TIT(2,17) = 'Delta+'
    Map2Hist_TIT(2,18) = 'Delta0'
    Map2Hist_TIT(2,19) = 'Delta-'


    Map2Hist_Meson(2, 1, 1) =  1 ! pi+
    Map2Hist_Meson(2, 1, 0) =  2 ! pi0
    Map2Hist_Meson(2, 1,-1) =  3 ! pi-

    Map2Hist_Meson(2, 3, 1) =  4 ! rho+
    Map2Hist_Meson(2, 3, 0) =  5 ! rho0
    Map2Hist_Meson(2, 3,-1) =  6 ! rho-

    Map2Hist_Meson(2, 2, 0) =  7 ! eta
    Map2Hist_Meson(2, 4, 0) =  7 ! sigma
    Map2Hist_Meson(2, 5, 0) =  7 ! omega
    Map2Hist_Meson(2, 6, 0) =  7 ! etaP
    Map2Hist_Meson(2, 7, 0) =  7 ! phi

    Map2Hist_Meson(2,10, 1) =  8 ! K+
    Map2Hist_Meson(2,10, 0) =  9 ! K0
    Map2Hist_Meson(2,11, 0) = 10 ! K~0
    Map2Hist_Meson(2,11,-1) = 11 ! K-

    Map2Hist_Baryon(2, 1, 1) = 12 ! p
    Map2Hist_Baryon(2, 1, 0) = 13 ! n

    Map2Hist_ABaryon(2, 1, 0) = 14 ! n~
    Map2Hist_ABaryon(2, 1,-1) = 15 ! p~

    Map2Hist_Baryon(2, 2, 2) = 16 ! Delta++
    Map2Hist_Baryon(2, 2, 1) = 17 ! Delta+
    Map2Hist_Baryon(2, 2, 0) = 18 ! Delta0
    Map2Hist_Baryon(2, 2,-1) = 19 ! Delta-

    ! (1c) map of iPartSet==3:

    Map2Hist_TIT(3, 1) = 'pi+'
    Map2Hist_TIT(3, 2) = 'pi0'
    Map2Hist_TIT(3, 3) = 'pi-'
    Map2Hist_TIT(3, 4) = 'eta'
    Map2Hist_TIT(3, 5) = 'omega'
    Map2Hist_TIT(3, 6) = 'phi'
    Map2Hist_TIT(3, 7) = 'J/psi'
    Map2Hist_TIT(3, 8) = 'K+'
    Map2Hist_TIT(3, 9) = 'K0'
    Map2Hist_TIT(3,10) = 'K~0'
    Map2Hist_TIT(3,11) = 'K-'
    Map2Hist_TIT(3,12) = 'D+'
    Map2Hist_TIT(3,13) = 'D0'
    Map2Hist_TIT(3,14) = 'D~0'
    Map2Hist_TIT(3,15) = 'D-'


    Map2Hist_Meson(3, 1, 1) =  1 ! pi+
    Map2Hist_Meson(3, 1, 0) =  2 ! pi0
    Map2Hist_Meson(3, 1,-1) =  3 ! pi-

    Map2Hist_Meson(3, 2, 0) =  4 ! eta
    Map2Hist_Meson(3, 5, 0) =  5 ! omega
    Map2Hist_Meson(3, 7, 0) =  6 ! phi
    Map2Hist_Meson(3, 9, 0) =  7 ! J/psi

    Map2Hist_Meson(3,10, 1) =  8 ! K+
    Map2Hist_Meson(3,10, 0) =  9 ! K0
    Map2Hist_Meson(3,11, 0) = 10 ! K~0
    Map2Hist_Meson(3,11,-1) = 11 ! K-

    Map2Hist_Meson(3,14, 1) = 12 ! D+
    Map2Hist_Meson(3,14, 0) = 13 ! D0
    Map2Hist_Meson(3,15, 0) = 14 ! D~0
    Map2Hist_Meson(3,15,-1) = 15 ! D-


    ! (1d) map of iPartSet==4:

    Map2Hist_TIT(4, 1) = 'pi+'
    Map2Hist_TIT(4, 2) = 'pi0'
    Map2Hist_TIT(4, 3) = 'pi-'
    Map2Hist_TIT(4, 4) = 'rho+'
    Map2Hist_TIT(4, 5) = 'rho0'
    Map2Hist_TIT(4, 6) = 'rho-'
    Map2Hist_TIT(4, 7) = 'eta'
    Map2Hist_TIT(4, 8) = 'sigma'
    Map2Hist_TIT(4, 9) = 'omega'
    Map2Hist_TIT(4,10) = 'etaP'
    Map2Hist_TIT(4,11) = 'phi'
    Map2Hist_TIT(4,12) = 'K+'
    Map2Hist_TIT(4,13) = 'K0'
    Map2Hist_TIT(4,14) = 'K~0'
    Map2Hist_TIT(4,15) = 'K-'
    Map2Hist_TIT(4,16) = 'K+*'
    Map2Hist_TIT(4,17) = 'K0*'
    Map2Hist_TIT(4,18) = 'K~0*'
    Map2Hist_TIT(4,19) = 'K-*'
    Map2Hist_TIT(4,20) = 'f_2'

    Map2Hist_Meson(4, 1, 1) =  1 ! pi+
    Map2Hist_Meson(4, 1, 0) =  2 ! pi0
    Map2Hist_Meson(4, 1,-1) =  3 ! pi-

    Map2Hist_Meson(4, 3, 1) =  4 ! rho+
    Map2Hist_Meson(4, 3, 0) =  5 ! rho0
    Map2Hist_Meson(4, 3,-1) =  6 ! rho-

    Map2Hist_Meson(4, 2, 0) =  7 ! eta

    Map2Hist_Meson(4, 4, 0) =  8 ! sigma

    Map2Hist_Meson(4, 5, 0) =  9 ! omega
    Map2Hist_Meson(4, 6, 0) = 10 ! etaPrime
    Map2Hist_Meson(4, 7, 0) = 11 ! phi

    Map2Hist_Meson(4,10, 1) = 12 ! K+
    Map2Hist_Meson(4,10, 0) = 13 ! K0
    Map2Hist_Meson(4,11, 0) = 14 ! K~0
    Map2Hist_Meson(4,11,-1) = 15 ! K-

    Map2Hist_Meson(4,12, 1) = 16 ! K+*
    Map2Hist_Meson(4,12, 0) = 17 ! K0*
    Map2Hist_Meson(4,13, 0) = 18 ! K~0*
    Map2Hist_Meson(4,13,-1) = 19 ! K-*

    Map2Hist_Meson(4,22, 0) = 20 ! f_2(1270)


    ! (1e) map of iPartSet==5:

    Map2Hist_TIT(5, 1) = 'pi+'
    Map2Hist_TIT(5, 2) = 'pi0'
    Map2Hist_TIT(5, 3) = 'pi-'
    Map2Hist_TIT(5, 4) = 'K+'
    Map2Hist_TIT(5, 5) = 'K0'
    Map2Hist_TIT(5, 6) = 'K~0'
    Map2Hist_TIT(5, 7) = 'K-'
    Map2Hist_TIT(5, 8) = 'p'
    Map2Hist_TIT(5, 9) = 'n'
    Map2Hist_TIT(5,10) = 'n~'
    Map2Hist_TIT(5,11) = 'p~'
    Map2Hist_TIT(5,12) = 'Lambda0'
    Map2Hist_TIT(5,13) = 'Sigma+'
    Map2Hist_TIT(5,14) = 'Sigma0'
    Map2Hist_TIT(5,15) = 'Sigma-'
    Map2Hist_TIT(5,16) = 'Xi0'
    Map2Hist_TIT(5,17) = 'Xi-'
    Map2Hist_TIT(5,18) = 'Omega-'
    ! strange antibaryons are omitted...


    Map2Hist_Meson(5, 1, 1) =  1 ! pi+
    Map2Hist_Meson(5, 1, 0) =  2 ! pi0
    Map2Hist_Meson(5, 1,-1) =  3 ! pi-

    Map2Hist_Meson(5,10, 1) =  4 ! K+
    Map2Hist_Meson(5,10, 0) =  5 ! K0
    Map2Hist_Meson(5,11, 0) =  6 ! K~0
    Map2Hist_Meson(5,11,-1) =  7 ! K-

    Map2Hist_Baryon(5, 1, 1) =  8 ! p
    Map2Hist_Baryon(5, 1, 0) =  9 ! n

    Map2Hist_ABaryon(5, 1, 0) = 10 ! n~
    Map2Hist_ABaryon(5, 1,-1) = 11 ! p~

    Map2Hist_Baryon(5,32, 0) = 12 ! Lambda0
    Map2Hist_Baryon(5,33, 1) = 13 ! Sigma+
    Map2Hist_Baryon(5,33, 0) = 14 ! Sigma0
    Map2Hist_Baryon(5,33,-1) = 15 ! Sigma-
    Map2Hist_Baryon(5,53, 0) = 16 ! Xi0
    Map2Hist_Baryon(5,53,-1) = 17 ! Xi-
    Map2Hist_Baryon(5,55,-1) = 18 ! Omega-


    if (.not.initFlag) return ! return in the second go


    ! (2) count entries of map

    do iS=1,nSet
       n = max(MAXVAL(Map2Hist_Meson(iS,:,:)),&
            & MAXVAL(Map2Hist_Baryon(iS,:,:)),&
            & MAXVAL(Map2Hist_ABaryon(iS,:,:)))

       Map2Hist_Meson(iS,:,:)   = n+1
       Map2Hist_Baryon(iS,:,:)  = n+2
       Map2Hist_ABaryon(iS,:,:) = n+3

       Map2Hist_TIT(iS,n+1) = 'other M'
       Map2Hist_TIT(iS,n+2) = 'other B'
       Map2Hist_TIT(iS,n+3) = 'other B~'

       Map2Hist_N(iS) = n+3
    end do

    initFlag=.false.

    goto 100 ! do it again

  end subroutine histMPf90_Init

  !****************************************************************************
  !****s* histMPf90/CreateHistMP
  ! NAME
  ! subroutine CreateHistMP(H, HName,x1,x2,bin, iPartSet)
  ! PURPOSE
  ! This is the Constructor of a multi-particle 1D-Histogram!
  ! Allocate Memory for the entries and put additional variables to their
  ! default.
  ! The parameter iPartSet specifies, which particles should be included.
  !
  ! INPUTS
  ! * type(histogramMP) :: H         -- Histogramm to be created
  ! * character*(*)     :: HName     -- Name of Histogram
  ! * real              :: x1,x2,bin -- Minimal/maximal value for x-coordinate
  !                                     to be considered, bin-width
  ! * integer           :: iPartSet  -- particle set to consider
  !
  ! possible values for "iPartSet" are:
  ! * 1: pions, kaons, nucleons (i.e. "stable" particles)
  ! * 2: most mesons (except charm), nucleons, Deltas
  ! * 3: stable/reconstructable mesons, no baryons
  ! * 4: all non-strange mesons
  ! * ... to be continued
  !
  !
  ! OUTPUT
  ! ---
  !****************************************************************************
  subroutine CreateHistMP(H, HName,x1,x2,bin, iPartSet)

    type(histogramMP),intent(inout) :: H
    character*(*),intent(in) :: HName
    real,intent(in) :: x1
    real,intent(in) :: x2
    real,intent(in) :: bin
    integer, intent(in) :: iPartSet

    integer :: L,nHist

    if (initFlag) call histMPf90_Init

    nHist = Map2HistMP_getN(iPartSet)
    if (nHist==0) then
       write(*,*) 'CreateHistMP: iPartSet invalid:',iPartSet
       stop
    end if
    H%iSet = iPartSet

    H%xMin = x1
    H%xMax = x2
    H%xBin = bin

    if (allocated(H%xExtreme)) deallocate(H%xExtreme)
    if (allocated(H%yVal)) deallocate(H%yVal)


    allocate(H%xExtreme(nHist,2))
    H%xExtreme(:,1) =  99e9
    H%xExtreme(:,2) = -99e9

    if (len(HName) > NameLength) then
       H%Name = HName(1:NameLength)
    else
       H%Name = HName
    end if

    L = nint( (x2-x1)/bin )+1
    allocate(H%yVal(nHist,-1:L,3))

    write(*,*) 'CreateHistMP:',nHist,L

    H%yVal = 0.

  end subroutine CreateHistMP

  !****************************************************************************
  !****s* histMPf90/ClearHistMP
  ! NAME
  ! subroutine ClearHistMP(H)
  ! PURPOSE
  ! Sets the histogram to zero again
  ! INPUTS
  ! * type(histogramMP) :: H         -- Histogram to be cleared
  ! OUTPUT
  ! H is changed
  !****************************************************************************
  subroutine ClearHistMP(H)
    type(histogramMP),intent(inout) :: H

    H%xExtreme(:,1) =  99e9
    H%xExtreme(:,2) = -99e9
    H%yVal = 0.
  end subroutine ClearHistMP

  !****************************************************************************
  ! cf. interface AddHistMP
  !****************************************************************************
  pure subroutine AddHistMP_I(H, ID,IC,anti, x,y,y2)
    type(histogramMP),intent(inout) :: H
    integer,intent(in)       :: ID,IC
    logical,intent(in)       :: anti
    real,intent(in)          :: x
    real,intent(in)          :: y
    real,intent(in),optional :: y2

    integer :: iH
    real :: yy

    yy = 0.
    if (present(y2)) yy = y2

    iH = Map2HistMP_I(ID,IC,anti, H%iSet)
    if (iH < 1) return

    call AddHistMP_iH(H,iH,x,y,yy)

  end subroutine AddHistMP_I
  !----------------------------------------------------------------------------
  pure subroutine AddHistMP_P(H, Part, x,y,y2)
    use particleDefinition

    type(histogramMP),intent(inout) :: H
    type(particle),intent(in) :: Part
    real,intent(in)          :: x
    real,intent(in)          :: y
    real,intent(in),optional :: y2

    call AddHistMP_I(H, Part%ID,Part%charge,Part%antiparticle, x,y,y2)

  end subroutine AddHistMP_P
  !----------------------------------------------------------------------------
  pure subroutine AddHistMP_KF(H, KF, x,y,y2)

    type(histogramMP),intent(inout) :: H
    integer,intent(in)       :: KF
    real,intent(in)          :: x
    real,intent(in)          :: y
    real,intent(in),optional :: y2

    integer :: iH
    real :: yy

    yy = 0.
    if (present(y2)) yy = y2

    iH = Map2HistMP_KF(KF, H%iSet)
    if (iH < 1) return

    call AddHistMP_iH(H,iH,x,y,yy)

  end subroutine AddHistMP_KF

  !****************************************************************************
  !****s* histMPf90/AddHistMP_IH
  ! NAME
  ! subroutine AddHistMP_IH(H, iH, x,y,y2)
  ! PURPOSE
  ! Add to the given histogram at the given x-value the weight y (y2).
  ! If y2 is given, also the second column is filled.
  ! ID and IC specify the particle.
  !
  ! INPUTS
  ! * type(histogramMP) :: H  -- Histogramm to be used
  ! * integer           :: iH -- nr of Histogram
  ! * real              :: x  -- x-value
  ! * real              :: y  -- weight to be added
  ! * real              :: y2 -- second weight to be added [OPTIONAL]
  ! OUTPUT
  ! H is changed
  ! NOTES
  ! this is an internal routine
  !****************************************************************************
  pure subroutine AddHistMP_IH(H, iH, x,y,y2)

    type(histogramMP),intent(inout) :: H
    integer, intent(in)      :: iH
    real,intent(in)          :: x
    real,intent(in)          :: y
    real,intent(in),optional :: y2

    integer :: iBin
    real :: yy

    if (iH < 1) return

    yy = 0.
    if (present(y2)) yy = y2


    if (x.lt.H%xExtreme(iH,1)) H%xExtreme(iH,1)=x
    if (x.gt.H%xExtreme(iH,2)) H%xExtreme(iH,2)=x

    if (x < H%xMin) then
       iBin = -1
    else if (x >= H%xMax) then
       iBin = 0
    else
       iBin = int( (x-H%xMin)/H%xBin )+1
    end if

    H%yVal(iH,iBin,1) = H%yVal(iH,iBin,1)+y
    H%yVal(iH,iBin,2) = H%yVal(iH,iBin,2)+1.
    H%yVal(iH,iBin,3) = H%yVal(iH,iBin,3)+yy

  end subroutine AddHistMP_IH

  !****************************************************************************
  !****f* histMPf90/Map2HistMP_getN
  ! NAME
  ! integer function Map2HistMP_getN(iSet)
  ! PURPOSE
  ! return the number of histograms for the given particle set.
  ! returns 0, if the particle set is invalid.
  !****************************************************************************
  integer function Map2HistMP_getN(iSet)

    integer, intent(in):: iSet

    Map2HistMP_getN = 0

    if (initFlag) call histMPf90_Init

    if (iSet<0) return
    if (iSet>nSet) return

    Map2HistMP_getN = Map2Hist_N(iSet)
  end function Map2HistMP_getN

  !****************************************************************************
  ! cf. interface Map2HistMP
  !****************************************************************************
  pure integer function Map2HistMP_I(ID,IC,anti, iSet)
    use particleDefinition

    integer, intent(in):: ID,IC
    logical, intent(in):: anti
    integer, intent(in):: iSet

    !
    ! ATTENTION: no check of array boundaries !!!
    !

    Map2HistMP_I = 0
    if (ID<1) return
    if (ID>122) return

    if (ID>100) then ! === MESON
       Map2HistMP_I = Map2Hist_Meson(iSet,ID-100,IC)
    else             ! === BARYON
       if (anti) then
          Map2HistMP_I = Map2Hist_ABaryon(iSet,ID,IC)
       else
          Map2HistMP_I = Map2Hist_Baryon(iSet,ID,IC)
       end if
    end if

  end function Map2HistMP_I
  !----------------------------------------------------------------------------
  pure integer function Map2HistMP_P(Part, iSet)
    use particleDefinition
    type(particle),intent(in) :: Part
    integer, intent(in):: iSet

    Map2HistMP_P = Map2HistMP_I(Part%ID,Part%charge,Part%antiparticle, iSet)
  end function Map2HistMP_P
  !----------------------------------------------------------------------------
  pure integer function Map2HistMP_KF(KF, iSet)
    use particleDefinition
    use ID_translation, only: KFtoBUU

    integer, intent(in):: KF
    integer, intent(in):: iSet

    !
    ! ATTENTION: no check of array boundaries !!!
    !

    integer:: ID,IC
    logical:: anti

    Map2HistMP_KF = 0

    call KFtoBUU(KF, ID,IC)
    if (ID<0) then
       anti = .true.
       ID = -ID
    else
       anti = .false.
    end if
    if (ID<1) return

    if (ID>100) then ! === MESON
       Map2HistMP_KF = Map2Hist_Meson(iSet,ID-100,IC)
    else             ! === BARYON
       if (anti) then
          Map2HistMP_KF = Map2Hist_ABaryon(iSet,ID,IC)
       else
          Map2HistMP_KF = Map2Hist_Baryon(iSet,ID,IC)
       end if
    end if

  end function Map2HistMP_KF

  !****************************************************************************
  !****s* histMPf90/WriteHistMP
  ! NAME
  ! subroutine WriteHistMP(H,iFile,add,mul,iColumn,DoAve,maxVal,H2,file,dump)
  ! PURPOSE
  ! Write out the histogram.
  !
  ! The entries are multiplied by 'mul' and 'add' is added.
  ! INPUTS
  ! * type(histogramMP) :: H     -- Histogramm to be used
  ! * integer           :: iFile -- File number output to redirect [OPTIONAL]
  ! * real              :: add   -- factor to add      [OPTIONAL]
  ! * real              :: mul   -- factor to multiply [OPTIONAL]
  ! * integer           :: iColumn -- column to write [OPTIONAL]
  ! * logical           :: DoAve -- write "average" (cf Notes) [OPTIONAL]
  ! * real              :: maxVal -- value used instead "divide by zero" [OPTIONAL]
  ! * type(histogram)   :: H2    -- Histogramm to divide by [OPTIONAL]
  ! * character*(*)     :: file  -- name of file to open and close [OPTIONAL]
  ! * logical           :: dump  -- if true, also dump it binary to file [OPTIONAL]
  ! OUTPUT
  ! write to file number
  ! NOTES
  ! The Histogram Data is not affected!!!
  !
  ! If you select "DoAve", the average value (built by dividing column 3 by
  ! column 1) is printed instead of the value of column 3 or column 1.
  ! This behaviour differs from the behaviour in module histf90 and hist2Df90,
  ! where on additional column is inserted in teh output.
  ! Setting this flag overrides settings according iColumn, add and mul.
  !****************************************************************************
  subroutine WriteHistMP(H,iFile_in,add,mul,iColumn,DoAve,maxVal,H2,file,dump)
    use CALLSTACK, only: TRACEBACK

    type(histogramMP),intent(in)          :: H
    integer,        intent(in),optional :: iFile_in
    real,           intent(in),optional :: add
    real,           intent(in),optional :: mul
    integer,        intent(in),optional :: iColumn
    logical,        intent(in),optional :: DoAve
    real,           intent(in),optional :: maxVal
    type(histogram),intent(in),optional :: H2
    character*(*),  intent(in),optional :: file
    logical,        intent(in),optional :: dump

    real,dimension(:,:),allocatable :: S
    real,dimension(:),  allocatable :: z
    integer             :: iBin, iBinMax,i,iC,nP, iFile

    real :: addFak
    real :: mulFak
    real :: maxZ!,Z
    logical :: writeZ

    real,allocatable :: yVal(:,:,:)

    addFak = 0.
    mulFak = 1.
    iC = 1
    maxZ = 99.0
    writeZ = .false.
    iFile = 65

    if (present(add)) addFak = add
    if (present(mul)) mulFak = mul
    if (present(iColumn)) iC = iColumn
    if (present(maxVal)) maxZ = maxVal
    if (present(DoAve)) writeZ = DoAve
    if (present(iFile_in)) iFile = iFile_in
    if (writeZ) iC = 1

    if (present(file)) then
       open(iFile, file=file, status='unknown')
       rewind(iFile)
    end if

    if (.not.allocated(H%yVal)) return

    nP = Map2Hist_N(H%iSet)
    iBinMax = ubound(H%yVal,dim=2)

    allocate(S(nP,3))
    allocate(z(nP))

    allocate(yVal(nP,-1:iBinMax,3))
    yVal = H%yVal

    if (present(H2)) then
       if (ubound(H2%yVal,dim=1).ne.iBinMax) then
          write(*,*) '  histogram=',H%Name
          write(*,*) ubound(H2%yVal,dim=1),iBinMax
          call TRACEBACK("ERROR, dim not equal! STOP")
          stop
       end if

       do iBin=-1,iBinMax
          if (H2%yVal(iBin,1) .eq. 0.0) then
             do i=1,nP
                if (yVal(i,iBin,1) .ne. 0.0) yVal(i,iBin,1) = maxZ
                if (yVal(i,iBin,3) .ne. 0.0) yVal(i,iBin,3) = maxZ
             end do
          else
             do i=1,nP
                yVal(i,iBin,1) = yVal(i,iBin,1)/H2%yVal(iBin,1)
                yVal(i,iBin,3) = yVal(i,iBin,3)/H2%yVal(iBin,1) ! we divide by column 1!
             end do
          end if
       end do
    end if

    S = SUM(yVal,dim=2) - yVal(:,-1,1:3) - yVal(:,0,1:3)


1000 format (1X,'#####',/,1X,'##### Histogram ',A40,/,1X,'#####')
1001 format (1X,'#',A11,1P,100E12.4,0P)

2000 format (1X,1P,E12.4,100E12.4,0P)

    write(iFile,1000) H%Name

    write(iFile,1001) 'Underflow',(yVal(i,-1,iC)*mulFak+addFak,i=1,nP)
    write(iFile,1001) 'Entries  ',(S(i,iC)*mulFak+addFak      ,i=1,nP)
    write(iFile,1001) 'Overflow ',(yVal(i, 0,iC)*mulFak+addFak,i=1,nP)
    write(iFile,1001)
    write(iFile,1001) 'minimum',  (H%xExtreme(i,1),i=1,nP)
    write(iFile,1001) 'maximum',  (H%xExtreme(i,2),i=1,nP)

    call WriteHistMP_Names(H, iFile)

    if (writeZ) then
       do iBin=1,iBinMax
          do i=1,nP
             if (yVal(i,iBin,2) > 0) then
                if (yVal(i,iBin,1) .eq. 0.0) then
                   z(i) = maxz
                else
                   z(i) = yVal(i,iBin,3)/yVal(i,iBin,1)
                end if
             else
                z(i) = maxz
             end if
          end do
          write(iFile,2000) &
               & H%xMin+(real(iBin)-0.5)*H%xBin,&
               & (z(i), i=1,nP)

       end do
    else
       do iBin=1,iBinMax
          write(iFile,2000) &
               & H%xMin+(real(iBin)-0.5)*H%xBin,&
               & (yVal(i,iBin,iC)/H%xBin*mulFak+addFak, i=1,nP)
       end do
    end if


    if (present(file)) then
       close(iFile)
    end if

    if ((present(dump)).and.(dump)) then
       if (present(file)) then
          call DumpHistMP(H,file//".bin",iFile)
       else
          call DumpHistMP(H,"DumpHist.bin",iFile)
       end if
    end if

  end subroutine WriteHistMP

  !****************************************************************************
  ! cf. interface WriteHistMP_Names
  !****************************************************************************
  subroutine WriteHistMP_Names_I(iSet, iFile)

    integer, intent(in) :: iSet
    integer, intent(in) :: iFile

    integer :: i,nP

1002 format (1X,'#####',/,1X,'#####         ',100(i2,':',A9))

    nP = Map2Hist_N(iSet)
    write(iFile,1002) ( i+1,Map2Hist_TIT(iSet,i),i=1,nP)

  end subroutine WriteHistMP_Names_I

  subroutine WriteHistMP_Names_H(H, iFile)

    type(histogramMP),intent(in) :: H
    integer,          intent(in) :: iFile

    call WriteHistMP_Names_I(H%iSet, iFile)

  end subroutine WriteHistMP_Names_H


  !****************************************************************************
  !****s* histMPf90/DumpHistMP
  ! NAME
  ! subroutine DumpHistMP(H,file,iFile)
  ! PURPOSE
  ! Write all the histogram information unformatted (i.e. binary) to a file
  !
  ! INPUTS
  ! * type(histogramMP) :: H     -- Histogramm to be used
  ! * character*(*)     :: file  -- name of file to open and close
  ! * integer,OPTIONAL  :: iFile -- File number output to redirect [OPTIONAL]
  ! OUTPUT
  ! H is written UNFORMATTED to the given file
  !
  !****************************************************************************
  subroutine DumpHistMP(H,file,iFile)

    type(histogramMP),intent(in)          :: H
    character*(*),    intent(in)          :: file
    integer,          intent(in),optional :: iFile

    integer :: iF

    iF=121
    if (present(iFile)) iF = iFile

    open(iF,file=file,status='UNKNOWN',form='UNFORMATTED')
    rewind(iF)

    write(iF) H%iSet
    write(iF) H%xMin,H%xMax,H%xBin,H%xExtreme
    write(iF) H%Name
    write(iF) H%yVal

    close(iF)

  end subroutine DumpHistMP

  !****************************************************************************
  !****s* histMPf90/FetchHistMP
  ! NAME
  ! subroutine FetchHistMP(H,file,iFile,flagOK)
  ! PURPOSE
  ! Read in all the histogram information previously dumped unformatted
  ! (i.e. binary) to a file
  !
  ! INPUTS
  ! * character*(*)   :: file  -- name of file to open and close
  ! * integer,OPTIONAL:: iFile -- File number input to redirect [OPTIONAL]
  ! OUTPUT
  ! * type(histogramMP) :: H     -- Histogramm to be used
  ! * logical           :: flagOK -- flag, if reading was okay [OPTIONAL]
  !
  ! H is read UNFORMATTED from the given file. Sizes are calculated as in
  ! CreateHist, also memory is allocated.
  !
  ! NOTES
  ! No checks about input are performed!
  !****************************************************************************
  subroutine FetchHistMP(H,file,iFile,flagOK)

    type(histogramMP),intent(inout)       :: H
    character*(*),  intent(in)          :: file
    integer,        intent(in),optional :: iFile
    logical,        intent(out),optional:: flagOK

    integer :: iF
    integer :: L
    integer :: ios

    if (initFlag) call histMPf90_Init

    iF=121
    if (present(iFile)) iF = iFile

    open(iF,file=file,status='UNKNOWN',form='UNFORMATTED',iostat=ios)
    if (ios.ne.0) then
       close(iF)
       if (present(flagOK)) flagOK=.false.
       return
    end if
    rewind(iF)

    read(iF,iostat=ios) H%iSet
    if (ios.ne.0) then
       close(iF)
       if (present(flagOK)) flagOK=.false.
       return
    end if

    if (allocated(H%xExtreme)) deallocate(H%xExtreme)
    if (allocated(H%yVal)) deallocate(H%yVal)

    allocate(H%xExtreme(Map2Hist_N(H%iSet),2))

    read(iF) H%xMin,H%xMax,H%xBin,H%xExtreme
    read(iF) H%Name

    L = nint( (H%xMax-H%xMin)/H%xBin )+1

    allocate(H%yVal(Map2Hist_N(H%iSet),-1:L,3))

    read(iF) H%yVal

    close(iF)
    if (present(flagOK)) flagOK=.true.

  end subroutine FetchHistMP

  !****************************************************************************

end module histMPf90
