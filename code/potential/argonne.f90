!******************************************************************************
!****m* /argonnev18
! NAME
! module argonnev18
!
! PURPOSE
! This module implements the Argonne V18 NN-potential for Deuterium and Deuterium wave functions
!******************************************************************************

module argonnev18

  implicit none

  private

  real, save, public :: argonne_max_kSpace = 0.  ! maximum of |psi(k)|^2 * k**2 in units of GeV^-1

  public :: argonne_WF_kSpace, argonne_WF_rSpace, argonne_deuteriumPot

  logical, parameter :: debug = .false.

contains

  !****************************************************************************
  !****f* argonnev18/argonne_WF_kSpace
  ! NAME
  ! function argonne_WF_kSpace (k_GeV, noPWave_in) result(wf)
  !
  ! PURPOSE
  ! This function calculates the  square of the momentum space wave function of a deuteron
  ! according to the Argonne V18 potential.
  !
  ! INPUTS
  ! * real    :: k_GeV -- momentum in units of GeV
  !
  ! OUTPUT
  ! * real    :: wf --momentum space wave function**2 in units of [GeV^-3]
  !
  ! NOTES
  ! * Normalization : integral(wf(k)*k**2 dk)=1
  !
  !****************************************************************************
  function argonne_WF_kSpace (k_GeV, noPWave_in) result(wf)
    use inputGeneral, only: path_to_input
    use constants, only: hbarc

    real, intent(in) :: k_Gev
    logical, intent(in), optional :: noPWave_in
    real :: wf

    logical, save :: initFlag=.true.
    real, dimension(1:160,1:3),save :: table
    integer :: index, i
    real :: k, max_check
    real, dimension(1:3) :: readin
    real, save :: maxK
    character(100) :: dummy
    logical :: noPWave

    if (present(noPWave_in)) then
       noPWave=noPWave_in
    else
       noPWave=.false.
    end if

    k=k_GEV/hbarc

    if (initFlag) then
       open(10,file=trim(path_to_input)//"/argonnePotential/deut.wfk")
       do i=1,1
          read(10,*) dummy
          if (debug) write(*,*) i, dummy
       end do
       do i=lbound(table,dim=1),ubound(table,dim=1)
          read(10,*) readin
          if (debug) write(*,*) i,' ',readin
          table(i,:)=readin
          maxK=table(i,1)
          ! Remember maximum of |psi(k)|^2 * k**2
          max_check=(table(i,2)**2+table(i,3)**2)*table(i,1)**2
          if (max_check.gt.argonne_max_kSpace) argonne_max_kSpace=max_check
       end do
       argonne_max_kSpace=argonne_max_kSpace/hbarc
       write(*,*) 'argonne_max_kSpace=',argonne_max_kSpace
       initflag=.false.
       close(10)
    end if

    if (k.gt.maxk) then
       wf=0.
       return
    else
       index=min(ubound(table,dim=1),max(1,1+nint(k/0.05)))
       if (noPWave) then
          wf=table(index,2)**2
       else
          wf=table(index,2)**2+table(index,3)**2
       end if
    end if
    wf=wf/hbarc**3

  end function argonne_WF_kSpace



  !****************************************************************************
  !****f* argonnev18/argonne_WF_rSpace
  ! NAME
  ! function argonne_WF_rSpace (r, noPWave_in) result(wf)
  !
  ! PURPOSE
  ! This function calculates the  square of the position space wave function of a deuteron
  ! according to the Argonne V18 potential.
  !
  ! INPUTS
  ! * real    :: r -- relative distance in units of [fm]
  !
  ! OUTPUT
  ! * real    :: wf -- position space wave function**2 in units of [fm^-1]
  !
  ! NOTES
  ! * Normalization : integral(wf(r) dr)=1
  !
  !****************************************************************************
  function argonne_WF_rSpace (r, noPWave_in) result(wf)
    use inputGeneral, only: path_to_input
    implicit none
    real, intent(in) :: r
    logical, intent(in), optional :: noPWave_in
    real :: wf

    logical, save :: initFlag=.true.
    real, dimension(1:1490,1:5),save :: table
    integer :: index, i
    real, save :: maxR
    real, dimension(1:5) :: readin
    character(100) :: dummy
    logical  :: noPWave

    if (present(noPWave_in)) then
       ! No P-Wave contribution
       noPWave=noPWave_in
    else
       noPWave=.false.
    end if

    if (initFlag) then
       open(10,file=trim(path_to_input)//"/argonnePotential/deut.wf")
       do i=1,7
          read(10,*) dummy
          if (debug)  write(*,'(I8,2A)') i, ' ' ,dummy
       end do
       do i=lbound(table,dim=1),ubound(table,dim=1)
          read(10,*) readin
          table(i,:)=readin
          maxR=table(i,1)
          if (debug)  write(*,*) i, table(i,:)
       end do
       initflag=.false.
       close(10)
    end if

    if (r.gt.maxR) then
       wf=0.
       return
    else
       index=min(ubound(table,dim=1),max(1,nint(r/0.01)))
       if (debug) write(*,*) r, max(1,nint(r/0.01)), index
       if (noPWave) then
          ! No P-Wave contribution
          WF=table(index,2)**2
       else
          WF=table(index,2)**2+table(index,4)**2
       end if
    end if

  end function argonne_WF_rSpace



  !****************************************************************************
  !****s* argonnev18/argonne_deuteriumPot
  ! NAME
  ! real function argonne_deuteriumPot(r)
  !
  ! PURPOSE
  ! * This function returns the AV18 potential for deuterium.
  !
  ! INPUTS
  ! * real,intent(in) :: r        ! distance berween of proton and neutron  [fm]
  !
  ! OUTPUT
  ! * Potential in [GeV]
  !****************************************************************************
  real function argonne_deuteriumPot(r)
    use av18pot, only: av18pw

    implicit none
    real,intent(in) :: r        ! distance berween of proton and neutron

    integer,parameter :: lpot=1 ! =av18
    integer,parameter :: l=0    ! S-Wave l=0
    integer,parameter :: s=1    ! Spin=1
    integer,parameter :: j=1    ! total angular momentum J=1
    integer,parameter :: t=0    ! total isospin
    integer,parameter :: t1z=1  ! z-Isospin of proton
    integer,parameter :: t2z=-1 ! z-Isospin of neutron
    real(8),dimension(2,2) :: vpw

    call av18pw(lpot,l,s,j,t,t1z,t2z,real(r,8),vpw)
    argonne_deuteriumPot=vpw(1,1)/1000. !Converting in GeV

  end function argonne_deuteriumPot




end module argonnev18
