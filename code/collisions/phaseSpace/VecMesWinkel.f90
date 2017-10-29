!******************************************************************************
!****m* /VecMesWinkel
! NAME
! module VecMesWinkel
! PURPOSE
! Incorporates the functions used for angular distributions of final states in
! gamma N -> V N and gamma N -> V Delta.
!******************************************************************************
module VecMesWinkel
  implicit none
  private

  public :: vecmesa, vecdelta

contains

  !****************************************************************************
  !****s*  VecMesWinkel/vecmesa
  ! NAME
  ! subroutine vecmesa(sqrts,idM,massM,massIn2,cost,iParam,successflag)
  ! PURPOSE
  ! determines scattering angle for gamma N --> V N
  ! INPUTS
  ! * real    :: sqrts -- sqrt(s)
  ! * integer :: idM -- ID of vector meson
  ! * real    :: massM -- mass of outgoing vector meson
  ! * real    :: massIn2 -- mass squared of incoming gamma
  !   (negative: Q^2!)
  ! * integer :: iParam -- select the distribution (see below)
  !
  ! OUTPUT
  ! * real              :: cost -- cos(theta) of scattering angle
  ! * logical, OPTIONAL :: successflag --
  !   false, if problems with in-medium masses
  ! NOTES
  ! * massIn2 may be negative: Q^2 !
  ! * The slope parameter B and the total cross section should be
  !   consistent!
  !
  ! The parameter 'iParam' selects the used distribution:
  ! * =0: 'old' parametrisation for gammaN->VN (cf. Effenberger PhD):
  !   Slope paremeter B according ABBHHM collab, PR 175, 1669 (1968)
  ! * =1: Pythia parametrisation:
  !   Slope parameter B=2*b_p+2*b_V+4*s**eps-4.2
  ! * =2: 'Donnachie, Landshoff'
  !   Select t according dsig/dt as given by VecMesWinkel/dsigdt, not
  !   by a given slope parameter
  ! * =3: as 1, but for rho and W<~6GeV slope parameter adjusted
  !   according CLAS experimental data [Morrow et al, EPJ A39, 5 (2009)]
  ! * =4: old Muehlich routines
  ! * =5: as 1, but modified for low W (used in rho0 Toy Init)
  ! * =6: W and Q2 dependent slope as fitted to VMD distribution of
  !   PYTHIA. Only valid for W<3.3GeV. (used in  rho0 Toy Init)
  !
  ! The default value for iParam is 3 (until Jan 2010 the default was 2)
  !****************************************************************************
  subroutine vecmesa(sqrts,idM,massM,massIn2,cost,iParam_,successflag)
    use random, only: rn, rnFlat, rnExp
    use baryonWidthMedium, only: get_MediumSwitch_coll
    use mesonWidthMedium, only: get_MediumSwitchMesons
    use VecMesWinkel_Muehlich, only: SelectCost
    use constants, only: mN

    real,    intent(in)  :: sqrts
    integer, intent(in)  :: idM
    real,    intent(in)  :: massM,MassIn2
    real,    intent(out) :: cost
    integer, intent(in)  :: iParam_
    logical, intent(out), optional :: successflag


    integer, parameter :: nvec = 4
    integer, parameter :: mbi = 6
    integer            :: iParam

    real, parameter :: extmin = -2.5
    ! min. t in V photoproduction before dsig/dt becomes flat
    ! (has also to be changed in vecmes)

    real, dimension(mbi,2), parameter :: bi = RESHAPE( (/&
         &1.8,  2.5,  3.5,  4.5,  5.8, 1e20,&     ! e_gamma
         &5.75, 5.43, 6.92, 8.1 , 7.9, 7.9 /) ,&  ! slope
         (/mbi,2/) )

    real, parameter :: bp = 2.3
    real, parameter :: eps=0.0808
    real, dimension(nvec), parameter :: bv = (/&
         & 1.4, 1.4, 1.4, 0.23 /)

    integer :: i,iv
    real :: pini2, pfin2
    real :: t,tmin, tmax
    real :: egamma,s,bel,bel1
    real :: sigmax,sigtest,sigtest2
    real :: mN2
    real :: h1,h2

    if (present(successflag)) successflag=.true.

    iParam = iParam_
    if (iParam==4 .and. idM/=105 .and. idM/=107) iParam=0   ! use iParam==4 only for omega and phi

    select case (idM)
    case (103)
       iv = 1
    case (105)
       iv = 2
    case (107)
       iv = 3
    case (109)
       iv =4
    case default
       write(*,*) 'problems in t-dependence, idM=',idM
       stop
    end select

    s = sqrts**2
    mN2 = mN**2

    egamma=(s-mN2-massIn2)/(2*mN)

    pini2 = ((s-mN2-massIn2)**2-4*massIn2*mN2)/(4*s)
    pfin2 = (s-(mN+massM)**2)*(s-(mN-massM)**2)/(4*s)
    if (pfin2.le.0..or.pini2.le.0.) then
       write(*,*) 'problems in vecmesa:',sqrts,massM,massIn2,pini2,pfin2
       if (present(successflag).and.(get_MediumSwitch_coll().or.get_MediumSwitchMesons())) then
          successflag=.false.
          return
       else
          stop
       end if
    end if

    h1 = sqrt((mN2+pini2)*(mN2+pfin2))
    h2 = sqrt(pini2*pfin2)

    tmax=2*(mN2-(h1-h2))
    tmin=2*(mN2-(h1+h2))

    if (tmax.le.tmin) then
       t=tmin
    else
       select case (iParam)
       case (0) !##### 'old' parametrisation (Effenberger) #####

          bel = 0.0
          do i=1,mbi
             if (egamma.lt.bi(i,1)) then
                bel=bi(i,2)
                exit
             end if
          end do
          t = rnExp(bel,tmax,tmin)

       case (1) !##### Pythia parametrisation #####

          bel=2*bp+2*bv(iv)+4*s**eps-4.2
          t = rnExp(bel,tmax,tmin)

       case (2) !##### 'Donnachie, Landshoff' ######

          if (iv.eq.4.and.sqrts.lt.7.6525) then
             t = rnExp(1.13,tmax,tmin)
          else
             sigmax = dsigdt(sqrts,tmax,iv)
             do
                t = rnFlat (tmin, tmax)
                sigtest  = rn()*sigmax
                sigtest2 = dsigdt(sqrts,max(extmin,t),iv)
                if (sigtest2.ge.sigtest) exit
             end do
          end if

       case (3) !##### Pythia, enhanced ######

          select case (iv)
          case (4)  ! J/Psi N -> J/Psi N
             bel=-1.64+0.83*log(s)  ! the parameterization from
                                    ! D. Kharzeev et al., EPJC 9, 459 (1999)
          case (1)
             bel  = 2*bp+2*bv(iv)+4*s**eps-4.2
             bel1 = (log(sqrts)-log(2.0))/(log(6.0)-log(2.0))*6.4+0.6
             bel = min(bel,bel1)

          case default
             bel  = 2*bp+2*bv(iv)+4*s**eps-4.2

          end select
          t = rnExp(bel,tmax,tmin)

       case (4) !##### Muehlich #####

          cost = SelectCost(sqrts,idM,1,1,massM,mN)
          return

       case (5) !##### Rho0 Toy Init #####

          bel=2*bp+2*bv(iv)+4*s**eps-4.2

          if (iv.eq.1) then
             bel1 = 2.2+3.0*(sqrts-2.0)
             bel = min(bel,bel1)
          end if

          t = rnExp(bel,tmax,tmin)

       case (6) !##### Rho0 Toy Init: Fit to PYTHIA-VMD  #####

          if (sqrts>3.3) then
             write(*,*) 'vecmesa: sqrts=',sqrts,' too large. STOP!'
             stop
          end if

          ! note: massIn2 = -Q2 !!
          bel = (4.576+1.101*massIn2+0.1656*massIn2**2) &
               & * (-3.968+3.364*sqrts-0.5095*sqrts**2)
          t = rnExp(bel,tmax,tmin)

       case (7)  !##### Simple Flat Distribution #####

          t = rnFlat (tmin, tmax)

       case default

          write(*,*) 'vecmesa: iParam=',iParam,' not possible. STOP!'
          stop

       end select

    end if

    cost=(t-2*mN2+2*h1)/(2*h2)

    if (abs(cost).gt.1) then
       write(*,*) 'problems in vecmesa cost',cost
       write(*,*) 'vecmesa:',iv,massM,sqrts,tmax,tmin,t
       cost=sign(1.,cost)
    end if

  end subroutine vecmesa

  !****************************************************************************
  !****if*  VecMesWinkel/dsigdt
  ! NAME
  ! real function dsigdt(sqrts,t,iv)
  ! PURPOSE
  ! Calculate the cross section dsigma/dt for gamma N -> V N according
  ! a parametrization of scattering amplitude from
  ! Donnachie and Landshoff
  !
  ! INPUTS
  ! * real    :: sqrts -- sqrt(s)
  ! * real    :: t     -- Mandelstam-t ( <0!)
  ! * integer :: iv    -- ID of vector meson
  ! OUTPUT
  ! function value: dsig/dt in mub/GeV^2
  !
  ! NOTES
  ! Reference:
  ! * A. Donnachie and P. V. Landshoff,
  !   ``Exclusive vector photoproduction: Confirmation of Regge theory,''
  !   Phys. Lett.  B 478, 146 (2000).
  !
  ! The parametrization is for photoproduction.
  !
  ! Range of validity: sqrt(s)=6.9...200 GeV.
  !
  ! Maybe one should upgrade to the newest work of the authors:
  ! * A.~Donnachie and P.~V.~Landshoff,
  !   ``Successful description of exclusive vector meson electroproduction,''
  !   arXiv:0803.0686 [hep-ph].
  !****************************************************************************

  real function dsigdt(sqrts,t,iv)
    use constants, only: pi, mN

    real, intent(in) :: sqrts
    real, intent(in) :: t
    integer, intent(in) :: iv

    complex :: ampli, ic
    real :: gv,fn,ap1,ar,ap0,alp1p,alp0p,alrp
    real :: expo1,expo0,expor,s,const

    const = 0.

    ic=(0.,1.)
    s=sqrts**2

    alp1p=0.25
    alp0p=0.1
    alrp=0.93

    expo1=1.08+alp1p*t-1.
    expor=0.55+alrp*t-1.
    expo0=1.44+alp0p*t-1.

    select case (iv)
    case (1)
       ap1=6.
       ar=15.9
       ap0=0.036
       gv=1./(1.-t/0.71)
    case (2)
       ! omega: modification by Muehlich, mentioned in EPJ A20 (2004) 499-508
       const=2532.*exp(-(s-mN**2)/(2.*0.21928*mN))+0.1532
       ap1=3.5
       ar=7.
       ap0=0.036
       gv=1./(1.-t/0.71)
    case (3)
       ap1=1.49
       ar=0.
       ap0=0.014
       gv=1./(1.-t/1.5)
    case (4)
       ap1=0.17
       ar=0.
       ap0=0.016
       gv=1.
    end select

    fn=(4.*mN**2.-2.79*t)/(4*mN**2-t)*(1./(1.-t/0.71))**2

    ampli=ic*fn*gv &
         &        *(ap1*(alp1p*s)**(expo1)*cexp(-0.5*ic*pi*expo1) &
         &        +ar*(alrp*s)**(expor)*cexp(-0.5*ic*pi*expor) &
         &        +ap0*(alp0p*s)**(expo0)*cexp(-0.5*ic*pi*expo0))

    dsigdt=(real(ampli)**2+aimag(ampli)**2)+const !mub/GeV^2

  end function dsigdt



  !****************************************************************************
  !****s*  VecMesWinkel/vecdelta
  ! NAME
  ! subroutine vecdelta(sqrts,idM,massD,massM,massIn2,cost,successflag)
  ! PURPOSE
  ! determines scattering angle for gamma N --> V Delta
  ! INPUTS
  ! * real    :: sqrts -- sqrt(s)
  ! * integer :: idM -- ID of vector meson
  ! * real    :: massD -- mass of outgoing baryon (Delta)
  ! * real    :: massM -- mass of outgoing vector meson
  ! * real    :: massIn2 -- mass squared of incoming gamma (negative: Q^2!)
  ! OUTPUT
  ! * real              :: cost -- cos(theta) of scattering angle
  ! * logical, OPTIONAL :: successflag -- false, if problems with in-medium masses
  ! NOTES
  ! * massIn2 may be negative: Q^2 !
  !****************************************************************************

  subroutine vecdelta(sqrts,idM,massD,massM,massIn2,cost,successflag)
    use random, only: rn, rnExp
    use constants, only: mN, mPi
    use baryonWidthMedium, only: get_MediumSwitch_coll
    use mesonWidthMedium, only: get_MediumSwitchMesons

    real, intent(in) :: sqrts
    integer, intent(in) :: idM
    real, intent(in) :: massD,massM,MassIn2
    real, intent(out) :: cost
    logical, intent(out), optional :: successflag

    real :: pini2, pfin2
    real :: t,tmin, tmax, bb
    real :: sigma,sigmax,sigtest
    real :: s, massN2, massD2
    real :: h1,h2

    if (present(successflag)) successflag=.true.

    s = sqrts**2
    massN2 = mN**2
    massD2 = massD**2

    pini2 = ((s-massN2-massIn2)**2-4*massIn2*massN2)/(4*s)
    pfin2 = (s-(massD+massM)**2)*(s-(massD-massM)**2)/(4*s)

    if (pfin2.le.0.or.pini2.le.0) then
       write(*,*) 'problems in vecdelta',sqrts,massD,massM,massIn2,pini2,pfin2
       if (present(successflag).and.(get_MediumSwitch_coll().or. get_MediumSwitchMesons())) then
          successflag=.false.
          return
       else
          stop
       end if
    end if

    h1 = sqrt((massN2+pini2)*(massD2+pfin2))
    h2 = sqrt(pini2*pfin2)

    tmax=massN2+massD2 - 2*(h1-h2)
    tmin=massN2+massD2 - 2*(h1+h2)

    if (tmax.le.tmin) then
       t=tmin
    else
       if (idM.eq.107) then
          bb=6.5
          sigmax=exp(bb*tmax)/(mPi**2-tmax)**2
          do
             t=rn()*(tmax-tmin)+tmin
             sigma=exp(bb*t)/(mPi**2-t)**2
             sigtest=rn()*sigmax
             if (sigma.gt.sigtest) exit
          end do
       else
          bb=6.
          t = rnExp(-bb,tmax,tmin)
       end if
    end if

    cost=(t-massN2-massD2+2*h1)/(2*h2)

    if (abs(cost).gt.1) then
       write(*,*) 'problems in vecdelta: cost=',cost
       write(*,*) 'vecdelta:',idM,massM,sqrts,tmax,tmin,t
       cost=sign(1.,cost)
    end if

  end subroutine vecdelta


end module VecMesWinkel
