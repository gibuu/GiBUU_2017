!******************************************************************************
!****m* /photonXSections
! NAME
! module photonXSections
!
! PURPOSE
! This module collects some routines for cross-section parametrizations
! for the photon-induced reactions (gamma+N).
!******************************************************************************
module photonXSections

  use mediumDefinition

  implicit none

  private

  public :: calcXS_gammaN2VN
  public :: calcXS_gammaN2VDelta
  public :: calcXS_gammaN2strange
  public :: setIParam
  !, stotphi, vmdphi


  !****************************************************************************
  !****g* photonXSections/iParam
  ! SOURCE
  !
  integer, save :: iParam = 2
  ! PURPOSE
  ! Switch to select the kind of parametrization for gamma N -> V N:
  ! * 1: "old parametrization", fit to experimental data, cf. Effenberger PhD, p.53
  ! * 2: Pythia, cf. Friberg/Sjöstrand hep-ph/0007314
  ! * 3: Donnachie, Landshoff [citation needed]
  !****************************************************************************


  !****************************************************************************
  !****g* photonXSections/omega_saphir
  ! SOURCE
  !
  logical, save :: omega_saphir = .true.
  ! PURPOSE
  ! If .true. an improved fit (to SAPHIR data) will be used for gamma N -> omega N.
  ! cf. "calcXS_omega_saphir"
  !****************************************************************************


  ! global variables used by calcXS_omega_saphir
  real, dimension(0:100) :: x,y
  real :: m_omega,sqrts,spot_omega
  type(medium), save :: med


contains

  !****************************************************************************
  !****s* photonXSections/setIParam
  ! NAME
  ! subroutine setIParam(i,os)
  ! PURPOSE
  ! switching by hand between some parametrisations
  ! INPUTS
  ! * integer           :: i  -- new value of "iParam"
  ! * integer, OPTIONAL :: os -- new value of "omega_saphir"
  !****************************************************************************
  subroutine setIParam(i,os)
    integer,intent(in) :: i
    logical,intent(in),optional :: os
    iParam = i
    if (present(os)) omega_saphir = os
  end subroutine


  !****************************************************************************
  !****s* photonXSections/readInput
  ! NAME
  ! subroutine readInput
  ! PURPOSE
  ! Reads input out of jobcard. Namelist 'photonXS'.
  !****************************************************************************
  subroutine readInput
    use output

    integer :: ios

    !**************************************************************************
    !****n* photonXSections/photonXS
    ! NAME
    ! NAMELIST photonXS
    ! PURPOSE
    ! Includes the parameters:
    ! * iParam
    ! * omega_saphir
    !**************************************************************************

    NAMELIST /photonXS/ iParam,omega_saphir

    call Write_ReadingInput('photonXS',0)
    rewind(5)
    read(5,nml=photonXS,IOSTAT=ios)
    call Write_ReadingInput('low_photo_induced',0,ios)

    write(*,*) "iParam:       ",iParam
    write(*,*) "omega_saphir: ",omega_saphir

    call Write_ReadingInput('photonXS',1)

  end subroutine readInput



  !****************************************************************************
  !****s* photonXSections/gammaN2VN_matrix
  ! NAME
  ! subroutine gammaN2VN_matrix(srts,i,spot,matrix)
  ! PURPOSE
  ! Calculates squared matrix elements for gamma N -> V N
  ! with different parametrizations, cf. 'iParam'.
  ! INPUTS
  ! * real, intent(in)    :: srts   -- SQRT(s)
  ! * integer, intent(in) :: i      -- select vector meson V: 1=rho, 2=omega, 3=phi, 4=J/Psi
  ! * real, intent(in)    :: spot   -- scalar potential for vector meson in final state
  ! OUTPUT
  ! * real, intent(out)          :: matrix  -- squared matrix element for gamma N -> V N
  ! NOTES
  ! For rho, omega and phi this routine returns matrix elements,
  ! for J/Psi it directly returns the cross section.
  !
  ! The used coupling constants gV are not modified due to shadowing!
  !****************************************************************************
  subroutine gammaN2VN_matrix(srts,i,spot,matrix)

    use constants, only: pi, alphaQED, mN
    use particleProperties, only: hadron
    use IdTable, only: omegaMeson

    real, intent(in) :: srts
    integer, intent(in) :: i    ! 1=rho, 2=omega, 3=phi, 4=J/Psi
    real, intent(in) :: spot
    real, intent(out) :: matrix

    logical, save :: first = .true.
    integer, parameter :: nvec = 4
    real :: b1(nvec),b2(nvec),b3(nvec),bv(2),wcm2,tmax,tmin
    real :: a1(nvec),a2(nvec),a3(nvec),sigv,bel,bp
    real :: egam,hhh,pfinal2
    integer :: para

    real, parameter :: bpsi=1.13, npsi=4.300, mpsi=3.097
    ! pomeron parameter :
    real, parameter :: eps = 0.0808   ! Phytia
    real, parameter :: eta = -0.4525
    real, parameter :: eps2 = 0.44    ! Additional Donnachie parameter
    !  VMD coupling constants
    real, parameter :: gv(nvec) = (/ 2.2 , 23.6 , 18.4 , 11.5 /)
    ! thresholds
    real, parameter :: thres2(nvec) = (/ 9.856, 14.78, 11.626, 16.2812 /)
    real, parameter :: thres3(nvec) = (/ 11.7210, 16.559, 8.18, 16.2812 /)

    if (first) then
      call readInput()
      first = .false.
    end if

    wcm2=srts**2

    ! low energies: fall back to old parametrization
    para = iParam
    select case (IParam)
    case (2)  ! Pythia
      if (wcm2<thres2(i)) para = 1
    case (3)  ! Donnachie
      if (wcm2<thres3(i)) then
        para = 1
      else if (i==4) then
        para = 2
      end if
    end select


    select case (para)
    case (1) !************ 1. ) old parametrization (Effenberger)

       select case (i)
       case (1)
          matrix = 0.16
       case (2)
          ! Final momentum of omega Meson and nucleon in CM
          pfinal2 = (wcm2-(mN+hadron(omegaMeson)%mass+spot)**2) &
                  * (wcm2-(mN-hadron(omegaMeson)%mass-spot)**2) / (4*wcm2)
          if (pfinal2.lt.0) pfinal2 = 0.
          matrix = 0.08*pfinal2/(2.*(srts-1.73)**2+pfinal2)
       case (3)
          matrix = 0.004
       case (4)
          matrix = 0.
       end select

    case (2) !************ 2.) PYTHIA parametrization

       ! notes:
       ! 1) sigma_el ~ sigma_tot^2/B_el
       ! 2) The used param for B_el is different from Falter PhD.
       ! 3) (hbar c)^2 = 0.3894 GeV^2 mb

       bp = 2.3
       bv(1:2) = (/ 1.4, 0.23 /)
       a1(1:3) = (/ 13.63, 13.63, 10.01 /)
       a2(1:3) = (/ 31.79, 31.79, -1.52 /)

      select case (i)
      case (1,2,3)
        !*     rho
        sigv=a1(i)*wcm2**eps+a2(i)*wcm2**eta
!        bel=2*bp+2*bv(1)+2.*0.25*log(0.25*wcm2)
        bel=2*bp+2*bv(1)+4*wcm2**eps-4.2
        matrix = wcm2*alphaQED*sigv**2/(gv(i)*16*pi*bel*0.389)
      case (4)
        !*     J/Psi
        a1(4)=a1(3)/10.
        a2(4)=a2(3)/10.
        egam=(wcm2-mN**2)/(2.*mN)       ! photon energy
        hhh=sqrt(max(0.,(wcm2-mpsi**2-mN**2)**2-4*mN**2*mpsi**2))
        tmin = mpsi**2-(wcm2-mN**2)/(2.*wcm2)*((wcm2+mpsi**2-mN**2)+hhh)
        tmax = mpsi**2-(wcm2-mN**2)/(2.*wcm2)*((wcm2+mpsi**2-mN**2)-hhh)
        if (wcm2.lt.23.7995) then
          matrix = npsi/mpsi**4/(16.*bpsi*pi) * (exp(bpsi*tmax)-exp(bpsi*tmin))
        else if (wcm2.lt.58.5612) then
          matrix = (4.5*(1-cos(0.15*(egam-11.)))+0.52)/1000.
        else
          matrix = 0.0020*srts**0.77
        end if
      end select

    case (3) !************ 3.) Donnachie parametrization

       !*     rho
      a1(1)=51.3083
      a2(1)=7.75827
      a3(1)=0.
      b1(1)=5.20554
      b2(1)=13.8542
      b3(1)=0.0200568
      !*     omega
      a1(2)=a1(1)*gv(1)/gv(2)
      a2(2)=a2(1)*gv(1)/gv(2)
      a3(2)=a3(1)*gv(1)/gv(2)
      b1(2)=b1(1)*gv(1)/gv(2)
      b2(2)=b2(1)*gv(1)/gv(2)
      b3(2)=b3(1)*gv(1)/gv(2)
      !*     phi
      a1(3)=0.370166
      a2(3)=0.122403
      a3(3)=0.00121056
      b1(3)=0.377458
      b2(3)=0.0471389
      b3(3)=0.00100955

      select case (i)
      case (1,2)
        !*     rho, omega
        if (wcm2.lt.151.721) then
          sigv=a1(i)*wcm2**(-0.917546)+a2(i)*wcm2**0.026853
        else
          sigv=b1(i)*wcm2**eps+b2(i)*wcm2**eta+b3(i)*wcm2**eps2
        end if
      case (3)
        !*     phi
        if (wcm2.lt.229.607) then
          sigv=a1(3)*wcm2**eps+a2(3)*wcm2**eta+a3(3)*wcm2**eps2
        else
          sigv=b1(3)*wcm2**eps+b2(3)*wcm2**eta+b3(3)*wcm2**eps2
        end if
      end select

      matrix = wcm2*sigv/1000.

    case default

       write(*,*) 'problems in calcXS_gammaN2VN, iParam > 2'
       stop

    end select

  end subroutine



  !****************************************************************************
  !****s* photonXSections/calcXS_gammaN2VN
  ! NAME
  ! subroutine calcXS_gammaN2VN(srts,media,sig,sigi)
  ! PURPOSE
  ! Produces cross section for
  !     gamma N -> V  N
  ! and
  !    gamma N -> V pi N
  ! INPUTS
  ! * real, intent(in)         :: srts   -- SQRT(s)
  ! * type(medium),intent(in)  :: media  -- medium at position
  ! OUTPUT
  ! * real, dimension(1:4)           :: sig  -- cross section for gamma N -> (rho,omega,phi,JPsi) N
  ! * real, dimension(1:4), OPTIONAL :: sigi -- cross section for gamma N -> (rho,omega,phi,JPsi) N Pion
  ! NOTES
  ! Returned cross sections are in microbarn.
  !****************************************************************************
  subroutine calcXS_gammaN2VN(srts,media,sig,sigi)

    use IdTable, only: jPsi, omegaMeson, phi, rho
    use twoBodyPhaseSpace, only: Integrate_2bodyPS_resonance
    use threeBodyPhaseSpace, only: Integrate_3bodyPS_resonance
    use mediumDefinition
    use mesonPotentialModule, only: vecMes_massShift
    use constants, only: mN, mPi

    type(medium),intent(in)  :: media
    real, intent(in) :: srts
    real, dimension(1:4),intent(out) :: sig
    real, dimension(1:4),intent(out),optional :: sigi

    integer, parameter ::  nvec=4
    real :: matrix,matrix2(nvec),pinitial,ps(5),ps2(2),spot(3)
    integer :: i

    integer, parameter, dimension(1:nvec) :: ires = (/rho,omegaMeson,phi,jPsi/)

    ! Medium Corrections, simple medium mass shifts of the vector mesons
    if (media%useMedium) then
      do i=1,3
        spot(i) = vecMes_massShift(ires(i),media%density)
      end do
    else
      spot = 0.
    end if

    pinitial=(srts**2-mN**2)/(2*srts)

    ! **************************************************************************
    !  gamma N -> V N
    ! **************************************************************************

    do i=1,3 ! rho, omega, phi
      call gammaN2VN_matrix(srts,i,spot(i),matrix)                    ! get matrix element
      ps = Integrate_2bodyPS_resonance (ires(i), srts, mN, spot(i))   ! get phase space
      sig(i) = 1000.*matrix*ps(1)/(pinitial*srts**2)
    end do

    ! for the J/Psi, gammaN2VN_matrix returns the cross section, not the matrix element
    call gammaN2VN_matrix(srts,4,0.,matrix)
    sig(4) = matrix

    ! omega meson: fit to SAPHIR data
    spot_omega = spot(2)
    med = media
    if (omega_saphir .and. srts<2.45) call calcXS_omega_saphir(srts,sig(2))

    do i=1,4
      if (sig(i) < 1e-3) sig(i) =0.0
    end do

    ! **************************************************************************
    !  gamma N -> V pi N
    ! **************************************************************************

    if (present(sigi)) then
      matrix2(1:3) = (/ 0.9548, 0.2087, 0. /)
      sigi(4)=0.

      do i=1,3
        ps2 = Integrate_3bodyPS_Resonance (srts, mN, mPi, ires(i), spot(i))
        sigi(i)=1000.*matrix2(i)*ps2(1)/(pinitial*srts)
      end do

      do i=1,4
        if (sigi(i) < 1e-3) sigi(i)=0.0
      end do
    end if

  end subroutine calcXS_gammaN2VN


  !****************************************************************************
  !****s* photonXSections/calcXS_gammaN2VDelta
  ! NAME
  ! subroutine calcXS_gammaN2VDelta(srts,sigi,media)
  !
  ! PURPOSE
  ! Calculates the cross section for gamma N -> V Delta,
  ! where V is a vector meson.
  !
  ! INPUTS
  ! * real, intent(in)         :: srts   ! SQRT(s)
  ! * type(medium),intent(in)  :: media  ! medium at position
  !
  ! OUTPUT
  ! * real, dimension(1:4) :: sigi  ! cross sections
  ! * sigi(1) -> rho Delta
  ! * sigi(2) -> omega Delta
  ! * sigi(3) -> phi Delta
  ! * sigi(4) -> JPsi Delta
  !
  ! NOTES
  ! Units of cross sections: microbarn
  !****************************************************************************
  subroutine calcXS_gammaN2VDelta(srts,sigi,media)
    use particleProperties, only: hadron
    use baryonWidth, only: FullWidthBaryon
    use twoBodyPhaseSpace, only: Integrate_2bodyPS_resonance
    use idTable, only: delta, rho, omegaMeson,JPSI,phi
    use constants, only: pi, GeVSquared_times_mb, mN, mPi
    use mediumDefinition
    use mesonPotentialModule, only: vecMes_massShift

    real, intent(in) :: srts
    real, intent(out),dimension(1:4) :: sigi
    type(medium),intent(in)  :: media

    integer, parameter :: nvec = 4
    integer, dimension(1:nvec), parameter :: ires = (/ rho, omegaMeson, phi, JPSI /)
    !real qsq2
    real matrix(nvec),pinitial,ps(5),integral!,pfinal2
    real spot(nvec)!,ampli
    integer i,j,nm!,idelta
    real minmass,maxmass,mres0,gamres0,spectral,intfac,y,ymax,ymin
    real dya,mass,gamtot!,pfinal,schwelle

    real, parameter :: a1=47.3                   ! mubGeV2
    real, parameter :: m1=2.3                    ! GeV
    real, parameter :: g1=1.8                    ! GeV
    real, parameter :: dy = 2*pi/10.
    real, dimension(1:nvec), parameter :: gv = (/ 20.8/2.2, 20.8/23.6*11., 12./18.4, 2.2/11.5 /)   ! only #3 and #4 are used

    real, dimension(1:nvec) :: mesmass

    mesmass = (/ 2*mPi, 2*mPi, 3*mPi,  hadron(JPsi)%mass /)  ! minimal meson masses

    if (media%useMedium) then
      do i=1,nvec
        spot(i) = vecMes_massShift(ires(i),media%density)
      end do
    else
       spot = 0.
    end if

    !*      ampli=0.0866**2/gv(3)
    !*      do i=1,nvec
    !*         matrix(i)=gv(i)*ampli
    !*      end do

    do i=1,2
       matrix(i) = a1/((srts-m1)**2+g1**2/4.) * GeVSquared_times_mb/1000. ! Muehlich omega paper: EPJ A20 (2004)
    end do
    matrix(3)=0.0866**2
    matrix(4)=0.0866**2 * gv(4)/gv(3)

    pinitial=(srts**2-mN**2)/2./srts

    minmass=mN+mPi
    mres0=hadron(delta)%mass
    gamres0=hadron(delta)%width

    do i=1,nvec

       maxmass=srts-mesmass(i)-spot(i)

       if (maxmass.le.minmass) then
          integral=0.
       else ! Integration over the Delta's mass??
          integral=0.
          ymax=2.*atan((maxmass-mres0)/gamres0*2.)
          ymin=2.*atan((minmass-mres0)/gamres0*2.)
          nm=max(int((ymax-ymin)/dy),1)
          dya=(ymax-ymin)/float(nm)
          do j=1,nm
             y=ymin+(float(j)-0.5)*dya
             mass=.5*tan(y/2.)*gamres0+mres0
             mass=min(max(mass,minmass),maxmass)

             gamtot=FullWidthBaryon(delta,mass) ! Vacuum width

             spectral=2./pi*mass**2*gamtot/((mass**2-   mres0**2)**2+gamtot**2*mass**2)
             intfac=gamres0/((mass-mres0)**2+gamres0**2/4.)
             ps = Integrate_2bodyPS_resonance (ires(i), srts, mass, spot(i))

             integral=integral+ps(1)*spectral*dya/intfac
          end do
       end if

       sigi(i) = matrix(i)*integral/pinitial/srts**2 * 1000./GeVSquared_times_mb  ! mub
       if (sigi(i)<1e-3) sigi(i)=0.

    end do

  end subroutine calcXS_gammaN2VDelta




!!$
!!$  !***************************************************************
!!$  !***************************************************************
!!$  !***************************************************************
!!$
!!$
!!$  subroutine calcXS_gammaN2VNa(srts,idv,massv,massv2,cost)
!!$    implicit none
!!$    include "common"
!!$    include "cominput"
!!$    integer idv,nvec,i,j,mbi,iv,iexv
!!$    parameter(nvec=4,mbi=5)
!!$    real bi(mbi,2),c1,x0,x,egamma,extmin
!!$    real srts,cost,pfinal,tmax,tmin,pfinal2,pini,massv,t,rn
!!$    real sigmax,sigtest,sigtest2,pini2,massv2,wcm2,bel,bp,bv(nvec)
!!$    logical tflag,flag
!!$    data ((bi(i,j),j=1,2),i=1,mbi)
!!$    *       e_gamma  slope
!!$    &     /1.8, 5.75,
!!$    &     2.5, 5.43,
!!$    &     3.5, 6.92,
!!$    &     4.5, 8.1,
!!$    &     5.8, 7.9/
!!$    !*********************************
!!$    iexv=2                    ! 0-> old parametrization for gammaN->VN
!!$    ! 1-> PYTHIA parametrization
!!$    ! 2-> Donnachie, Landshoff
!!$    ! WARNING: has also to be changed in calcXS_gammaN2VN
!!$    extmin=-2.5               ! min. t in V photoproduction before dsig/dt
!!$    ! becomes flat
!!$    !*********************************
!!$    bp=2.3
!!$    do i=1,nvec-1
!!$       bv(i)=1.4
!!$    end do
!!$    bv(4)=0.23
!!$
!!$    if(idv.eq.103) then
!!$       iv=1
!!$    else if(idv.eq.105) then
!!$       iv=2
!!$    else if(idv.eq.107) then
!!$       iv=3
!!$    else if(idv.eq.109) then
!!$       iv=4
!!$    else
!!$       write(*,*) 'problems in t-dependence, idv=',idv
!!$    end if
!!$
!!$
!!$    egamma=(srts**2-rmass**2-massv2)/2./rmass
!!$    flag=.true.
!!$    i=0
!!$    do while(flag)
!!$       i=i+1
!!$       if(i.eq.mbi.or.egamma.lt.bi(i,1)) flag=.false.
!!$    end do
!!$    bel=bi(i,2)
!!$
!!$    if(iexv.eq.1) then
!!$       wcm2=srts**2
!!$       bel=2*bp+2*bv(iv)+2.*0.25*log(0.25*wcm2)
!!$    end if
!!$
!!$    pini2=((srts**2-rmass**2-massv2)**2-4*massv2*rmass**2)/
!!$    &     (4.*srts**2)
!!$    pfinal2=(srts**2-(rmass+massv)**2)*(srts**2-(rmass-massv)**2)/
!!$    &     4./srts**2
!!$    if(pfinal2.le.0.or.pini2.le.0) then
!!$       write(*,*)'problems in calcXS_gammaN2VNa',srts,rmass,massv,pini2,pfinal2
!!$       stop
!!$    end if
!!$    pfinal=sqrt(pfinal2)
!!$    pini=sqrt(pini2)
!!$
!!$    tmax=2.*rmass**2-2.*(sqrt(rmass**2+pini**2)*sqrt(rmass**2+
!!$    &     pfinal2)-pini*pfinal)
!!$    tmin=2.*rmass**2-2.*(sqrt(rmass**2+pini**2)*sqrt(rmass**2+
!!$    &     pfinal2)+pini*pfinal)
!!$
!!$    if(tmax.le.tmin) then
!!$       t=tmin
!!$    else if(iexv.le.1) then
!!$       c1=bel/(exp(bel*tmax)-exp(bel*tmin))
!!$       x0=-c1/bel*exp(bel*tmin)
!!$       x=rn(iseed)
!!$       t=log(bel/c1*(x-x0))/bel
!!$    else
!!$       if(iv.eq.4.and.srts.lt.7.6525) then
!!$          c1=1.13/(exp(1.13*tmax)-exp(1.13*tmin))
!!$          x0=-c1/1.13*exp(1.13*tmin)
!!$          x=rn(iseed)
!!$          t=log(1.13/c1*(x-x0))/1.13
!!$       else
!!$          *            if(tmin.lt.extmin) tmin=extmin
!!$          call dsig(srts,tmax,sigmax,iv)
!!$          tflag=.true.
!!$          do while(tflag)
!!$             t=rn(iseed)*(tmax-tmin)+tmin
!!$             sigtest=rn(iseed)*sigmax
!!$             if(t.lt.extmin) then
!!$                call dsig(srts,extmin,sigtest2,iv)
!!$             else
!!$                call dsig(srts,t,sigtest2,iv)
!!$             end if
!!$             if(sigtest2.ge.sigtest) tflag=.false.
!!$          end do
!!$       end if
!!$    end if
!!$
!!$    cost=(t-2*rmass**2+2*sqrt(rmass**2+pini**2)*sqrt(rmass**2+
!!$    &     pfinal2))/2./pini/pfinal
!!$
!!$    if(abs(cost).gt.1) then
!!$       write(*,*)'problems in calcXS_gammaN2VNa cost',cost
!!$       write(*,*)'calcXS_gammaN2VNa:',iv,massv,srts,tmax,tmin,t
!!$       cost=sign(1.,cost)
!!$    end if
!!$
!!$    return
!!$  end subroutine calcXS_gammaN2VNa
!!$
!!$
!!$
!!$  !*******************************************************
!!$  subroutine dsig(srts,t,dsigdt,iv)
!!$    implicit none
!!$    integer nvec,i,j,iv
!!$    parameter(nvec=4)
!!$    complex ampli,ic
!!$    real t,gv,fn,mn,ap1,ar,ap0,alp1p,alp0p,alrp
!!$    real srts,pi,expo1,expo0,expor,dsigdt,s
!!$    parameter(mn=0.938,pi=3.141592654)
!!$
!!$    ic=(0.,1.)
!!$    s=srts**2
!!$
!!$    alp1p=0.25
!!$    alp0p=0.1
!!$    alrp=0.93
!!$
!!$    expo1=1.08+alp1p*t-1.
!!$    expor=0.55+alrp*t-1.
!!$    expo0=1.44+alp0p*t-1.
!!$
!!$    if(iv.eq.1.or.iv.eq.2) then
!!$       ap1=6.
!!$       ar=15.9
!!$       ap0=0.036
!!$       gv=1./(1.-t/0.71)
!!$    else if(iv.eq.3) then
!!$       ap1=1.49
!!$       ar=0.
!!$       ap0=0.014
!!$       gv=1./(1.-t/1.5)
!!$    else if(iv.eq.4) then
!!$       ap1=0.17
!!$       ar=0.
!!$       ap0=0.016
!!$       gv=1.
!!$    end if
!!$
!!$    fn=(4.*mn**2.-2.79*t)/(4*mn**2-t)*(1./(1.-t/0.71))**2
!!$
!!$    ampli=ic*fn*gv                                                                          &
!!$         &        *(ap1*(alp1p*s)**(expo1)*cexp(-0.5*ic*pi*expo1)        &
!!$         &        +ar*(alrp*s)**(expor)*cexp(-0.5*ic*pi*expor)                &
!!$         &        +ap0*(alp0p*s)**(expo0)*cexp(-0.5*ic*pi*expo0))
!!$    dsigdt=(real(ampli)**2+aimag(ampli)**2) !mub/GeV^2
!!$
!!$    return
!!$  end subroutine dsig
!!$  !**********************************************************************
!!$  subroutine calcXS_gammaN2VNm(srts,minmass,mass1,mass2,id1,mass,rhor,pr, spot)
!!$    implicit none
!!$!    include "common"
!!$!    include "cominput"
!!$!    include "comwidth"
!!$    real srts,minmass,mass1,mass2,mass,massmax2,psam,x,xx,gamtot,rn,  &
!!$         &     spectral,psa,ratio(nmesch2+nmesch3),bwmes,rhor,pr,spot
!!$    integer id1,idm
!!$    logical flag
!!$
!!$    idm=id1-idbmax
!!$    if(mesprop1(idm,2).lt.1e-03) then
!!$       mass=mesprop1(idm,1)
!!$    else
!!$       massmax2=srts-mass1-mass2-spot
!!$       if(massmax2.le.minmass) then
!!$          write(*,*)'problems in calcXS_gammaN2VNm masses',  minmass,massmax2
!!$          stop
!!$       end if
!!$
!!$       call bops3(psam,srts,mass1,mass2,minmass,629)
!!$
!!$       flag=.true.
!!$       do while(flag)
!!$          x=rn(iseed)
!!$          xx=rn(iseed)
!!$          mass=minmass+(massmax2-minmass)*x
!!$
!!$          if((iwsc2.ge.1.and.idm.eq.3).or.   (iwsc2.ge.2.and.idm.eq.5)) then
!!$             rhores=rhor
!!$             pres=pr
!!$             gamtot=bwmes((mass+spot)**2,idm,0,0,0,ratio,1,629)
!!$          else
!!$             gamtot=bwmes((mass+spot)**2,idm,0,0,0,ratio,0,629)
!!$          end if
!!$          spectral=mass**2*gamtot/((mass**2-   mesprop1(idm,1)**2)**2+gamtot**2*mass**2)
!!$
!!$          call bops3(psa,srts,mass1,mass2,mass+spot,629)
!!$
!!$          if(spectral*psa.gt.xx*psam/mesprop1(idm,2)) then
!!$             flag=.false.
!!$          end if
!!$       end do
!!$    end if
!!$    return
!!$  end subroutine calcXS_gammaN2VNm
!!$  !*****************************************************************************




  !****************************************************************************
!   real function stotphi(s)
!     use constants, only: mN, pi
!     real :: s,tmin,tmax,t,qcm,sigtot,dsdop,dsdon,sth
!     real, parameter :: mp = 1.02
!     integer :: i
!     integer, parameter :: imax = 1000
!
!     sth=(mp+mn)**2
!     sigtot=0.
!     if(s.gt.sth)then
!        qcm=sqrt((s-(mn+mp)**2)*(s-(mn-mp)**2))/2./sqrt(s)
!        tmin=0.
!        tmax=-4.*qcm**2
!        do i=0,imax
!           t=float(i)*(tmax-tmin)/float(imax)+tmin
!           dsdop=vmdphi(s,tmin,tmax,t,1)
!           dsdon=vmdphi(s,tmin,tmax,t,0)
!           sigtot=sigtot+(dsdop+dsdon)/2.*2.*pi*2./float(imax)
!        end do
!     else
!        qcm=sqrt((sth-(mn+mp)**2)*(sth-(mn-mp)**2))/2./sqrt(sth)
!        tmin=0.
!        tmax=-4.*qcm**2
!        do i=0,imax
!           t=float(i)*(tmax-tmin)/float(imax)+tmin
!           dsdop=vmdphi(sth,tmin,tmax,t,1)
!           dsdon=vmdphi(sth,tmin,tmax,t,0)
!           sigtot=sigtot+(dsdop+dsdon)/2.*2.*pi*2./float(imax)
!        end do
!     end if
!     stotphi=sigtot
!
!     return
!   contains
!
!     !************************************************************************
!     !*     A.I.TITOV, T.S.H.LEE, H.TOKI, O.STRELTSOVA
!     !*     PHYS.REV.C 60, 035205 (1999)
!     !*     "STRCTURE OF THE PHI PHOTOPRODUCTION AMPLITUDE AT A FEW GEV"
!
!     real function vmdphi(s,tmin,tmax,t,fn)
!       use constants, only: mN, pi
!       integer imax,ichannel,nmax,fn!,i,j,l,ifail,n
!       parameter(imax=5000,ichannel=8,nmax=20)
!       real tmin,tmax
!       !real p,pb,egamma
!       real s,t,u,mp,eqed,s1,s2!,kcm,qcm
!       parameter(mp=1.02,eqed=0.303)
!       real alph10,alph20,alphb,alph1,msqr!,alph2,dsig,sth
!       parameter(alph10=1.08,alph20=-0.75,alphb=0.25)
!       complex ic,x1,x1d,x3,x4,x5,x3d,x5d,sum(ichannel)!,x2,x2d,x6,x6d,x4d
!       real mpi,gpinn,gphigapi,meta,getann,gphigaeta
!       parameter(mpi=0.136,gphigapi=-0.141,meta=0.547)
!       parameter(getann=-3.527,gphigaeta=-0.707)
!       real c1,c2,f,f1,bslope!,f2
!       parameter(c1=2.09,s1=1.6,bslope=1.,c2=0,s2=1)
!       !real kq,pq,theta,kp
!       complex y1,y2,y3,y4,gphinn,kapan,kapaphi!,sum8
!       parameter(gphinn=-0.24,kapaphi=0.2)
!       real fs,fu,lambda
!       parameter(lambda=1.87)
!       integer iform
!       real a1,a2!,mp1,f7
!       parameter(a1=0.9,a2=0.1)
!       real dsdo,gphi!,cost,pmin,pmax,plab,sigtot
!       parameter(gphi=12.88)
!
!       ic=cmplx(0.,1.)
!
!       if(.not.(a1+a2.eq.1.)) then
!          write(*,*)'problems norm a1+a2 not eq 1'
!          stop
!       end if
!
!       iform=1
!       if(abs(fn).eq.1)then
!          kapan=1.79
!          gpinn=-13.26
!       else if(fn.eq.0)then
!          kapan=-1.91
!          gpinn=13.26
!       end if
!
!       !sth=(mp+mn)**2
!       !qcm=sqrt((s-(mn+mp)**2)*(s-(mn-mp)**2))/2./sqrt(s)
!       !kcm=qcm
!       u=2.*mn**2+2.*mp**2-s-t
!       !kp=(s-mn**2-mp**2)/2.
!       !kq=mp**2-t/2.
!       !pq=(mn**2+mp**2-u)/2.
!       alph1=alph10+alphb*t
!       f=(4.*mn**2-2.8*t)/((4.*mn**2-t)*(1.-(t/0.7))**2)     &
!            &     *exp(bslope*(t-tmin)/2.)
!       f1=1.            &
!            &     /sqrt((2*mn**4*(2*mp**2 + mp**2 - t) +                                        &
!            &     mp**4*(4*mp**2 - 4*s - t) +                                                            &
!            &     2*mp**2*(mp**4 - 3*mp**2*s + 2*s**2 + 3*s*t+t**2)+                 &
!            &     2*mn**2*(mp**4 + mp**2*(5*mp**2 - 4*s - t) +                            &
!            &     2*s*(-mp**2 + t)) -                                                                &
!            &     (mp**2 - t)*(-2*s**2 - 2*s*t-t**2+mp**2*(2*s+t)))/               &
!            &     (4.*mp**2))
!       x1=c1*f*f1*cexp(-ic*pi/2.*alph1)*((s-s1)*alphb)**(alph1)
!       x1d=conjg(x1)
!       x3=-ic/(t-mpi**2)*eqed/mp*gpinn*gphigapi   *(0.7**2-mpi**2)/(0.7**2-t)*(0.77**2-mpi**2)/(0.77**2-t)
!       x3d=conjg(x3)
!       x4=-ic/(t-meta**2)*eqed/mp*getann*gphigaeta   *(1.**2-meta**2)/(1.**2-t)*(0.9**2-meta**2)/(0.9**2-t)
!       !x4d=conjg(x4)
!       x5=x3+x4
!       x5d=conjg(x5)
!       fs=lambda**4/(lambda**4+(s-mn**2)**2)
!       fu=lambda**4/(lambda**4+(u-mn**2)**2)
!       if(iform.eq.0)then
!          fs=1.
!          fu=1.
!       end if
!       y1=ic*kapan/2./mn
!       y2=ic*kapaphi/2./mn
!       y3=eqed*gphinn/(s-mn**2)*(a1*fs+a2*fu)
!       y4=eqed*gphinn/(u-mn**2)*(a1*fs+a2*fu)
!
!       !*     SCALAR MESON + POMERON I EXCHANGE
!
!       x3=x5
!       x3d=x5d
!       sum(6)=(6*mp**6*x1*x1d + 2*mn**4*(3*mp**2 - t)*x1*x1d -           &
!            &     t*(2*s**2 + 2*s*t + t**2)*x1*x1d +                                              &
!            &     2*mn**2*(6*mp**4 + 2*s*t - mp**2*(6*s + t))*x1*x1d -             &
!            &     2*mp**4*(6*s*x1*x1d + t*x1*x1d - 2*t**2*x3*x3d) +                 &
!            &     mp**2*(6*s**2*x1*x1d + 10*s*t*x1*x1d + 4*t**2*x1*x1d -       &
!            &     t**3*x3*x3d))/mp**2
!
!       !*     S- AND U-CHANNEL NUCLEON EXCHANGE (PHI RADIATION)
!
!       sum(7)=vsum1(s,t,u,y1,y2,y3,y4,fn)
!
!       !*     ALL TOGETHER FOR MODEL A (POMERON I, PI, ETA, S-CAHNNEL, U-CHANNEL)
!
!       sum(8)=sum(6)+sum(7)+vsum2(s,t,u,y1,y2,y3,y4,x1,x5,fn)
!       msqr=real(sum(8))
!       dsdo=1./(64.*pi**2*4.*s)*msqr*0.389*(gphi/eqed)**2*2./3.
!       vmdphi=dsdo
!
!       return
!     end function vmdphi
!     !*******************************************************************************
!
!     complex function vsum1(s,t,u,y1,y2,y3,y4,qn)
!       use constants, only: mN
!       integer qn
!       real s,t,u,fn
!       real, parameter :: mp = 1.02
!       complex y1,y3,y2,y4!,y1d,y2d,y3d,y4d
!       complex part1,part2!,part3
!
!       fn=real(qn)
!       !y1d=conjg(y1)
!       !y2d=conjg(y2)
!       !y3d=conjg(y3)
!       !y4d=conjg(y4)
!
!       part1=2*y3*(2*fn**2*(2*s*(s + t + s*t*y2**2)*y3 +                &
!            &     mp**6*y2**2*(-y3 + y4) + (0,6)*mn**5*y2*(y3 + y4) +    &
!            &     4*mn**6*y2**2*(y3 + y4) +        &
!            &     mp**4*(2*(1 + s*y2**2)*y3 + (-2*s + t)*y2**2*y4) -         &
!            &     mn**4*((-6 + 9*mp**2*y2**2 + 8*s*y2**2 - 2*t*y2**2)*y3 +         &
!            &     (7*mp**2*y2**2 + 4*(-2 + 2*s*y2**2 + t*y2**2))*y4) -         &
!            &     (0,6)*mn**3*y2*(t*(-y3 + y4) + 2*s*(y3 + y4) +         &
!            &     mp**2*(5*y3 + 2*y4)) -         &
!            &     (0,1)*mn*y2*(6*s*t*(y3 - y4) + 6*mp**4*(2*y3 - y4) +         &
!            &     t**2*y4 - 6*s**2*(y3 + y4) +         &
!            &     mp**2*(-6*s*(y3 - 2*y4) - 13*t*y4)) -         &
!            &     mp**2*(s**2*y2**2*(y3 - y4) + 4*t*y4 +         &
!            &     s*((4 + t*y2**2)*y3 - t*y2**2*y4)) +         &
!            &     mn**2*(-4*mp**4*y2**2*(3*y3 - 2*y4) +         &
!            &     4*s**2*y2**2*(y3 + y4) -         &
!            &     2*t*(y3 + 2*s*y2**2*y3 + y4 - 2*s*y2**2*y4) +         &
!            &     mp**2*((12 - 6*s*y2**2 + t*y2**2)*y3 +         &
!            &     (4 - 10*s*y2**2 + 3*t*y2**2)*y4))) +         &
!            &     fn*y1*(18*mn**6*y2*y3 -         &
!            &     (0,4)*mn**5*(3*(-1 + 2*mp**2*y2**2 + t*y2**2)*y3 +         &
!            &     (-3 + 2*mp**2*y2**2 - t*y2**2)*y4) +         &
!            &     mn**4*y2*(-54*s*y3 + t*(36*y3 + y4) +         &
!            &     mp**2*(-36*y3 + 10*y4)) -         &
!            &     mn**2*y2*(-54*s**2*y3 + 4*mp**2*(5*s - 12*t)*y4 +         &
!            &     15*t**2*y4 + 2*s*t*(18*y3 + y4) + 2*mp**4*(63*y3 + 10*y4))        &
!            &     + y2*(10*mp**6*y4 + 2*mp**2*s*(18*s*y3 + 5*s*y4 + 4*t*y4) +         &
!            &     mp**4*(-18*s*y3 - 20*s*y4 + 11*t*y4) +         &
!            &     s*(-18*s**2*y3 + s*t*y4 + t**2*y4)) +         &
!            &     (0,2)*mn**3*(2*mp**4*y2**2*(21*y3 + 4*y4) +         &
!            &     mp**2*((-30 + 24*s*y2**2 - 3*t*y2**2)*y3 +         &
!            &     4*(-3 + 2*s*y2**2 - 2*t*y2**2)*y4) +         &
!            &     2*(t*(3*y3 + (-3 + t*y2**2)*y4) +         &
!            &     s*(6*(-1 + t*y2**2)*y3 - 2*(3 + t*y2**2)*y4))) +         &
!            &     (0,1)*mn*(4*mp**6*y2**2*(3*y3 - 2*y4) +         &
!            &     4*mp**4*(3*(-2 + s*y2**2)*y3 +         &
!            &     (3 + 4*s*y2**2 - 4*t*y2**2)*y4) -         &
!            &     2*(t**2*y4 + s**2*        &
!            &     (6*(-1 + t*y2**2)*y3 - 2*(3 + t*y2**2)*y4) -         &
!            &     2*s*t*(-3*y3 + (3 + t*y2**2)*y4)) +         &
!            &     mp**2*(t*(26 + 3*t*y2**2)*y4 - 8*s**2*y2**2*(3*y3 + y4) +         &
!            &     2*s*(3*(2 + t*y2**2)*y3 - 4*(3 + 2*t*y2**2)*y4)))) +         &
!            &     y1**2*(4*mn**8*y2**2*y3 +         &
!            &     4*s**2*(t + s**2*y2**2 + s*t*y2**2)*y3 -         &
!            &     2*mp**2*s*(4*s**2*y2**2*y3 + t*(y3 - y4) +         &
!            &     s*(y3 + 2*t*y2**2*y3 - y4)) + mp**8*y2**2*(y3 - y4) +         &
!            &     mp**4*(s*(4 + t*y2**2)*(y3 - y4) + s**2*y2**2*(5*y3 - y4) +         &
!            &     2*t*y4) - 2*mp**6*        &
!            &     (y3 + s*y2**2*y3 - (1 + s*y2**2 - 2*t*y2**2)*y4) -         &
!            &     (0,4)*mn**5*y2*(t*(3*y3 - y4) + 2*mp**2*(3*y3 + y4)) +         &
!            &     (0,2)*mn**3*y2*(2*mp**4*(21*y3 + 4*y4) +         &
!            &     mp**2*(24*s*y3 - 3*t*y3 + 8*s*y4 - 8*t*y4) +         &
!            &     2*t*(6*s*y3 - 2*s*y4 + t*y4)) -         &
!            &     4*mn**6*((-2 + 6*mp**2*y2**2 + 4*s*y2**2 + t*y2**2)*y3 -         &
!            &     2*(y4 + t*y2**2*y4)) +         &
!            &     (0,1)*mn*y2*(4*mp**6*(3*y3 - 2*y4) +         &
!            &     4*mp**4*(3*s*y3 + 4*s*y4 - 4*t*y4) +         &
!            &     4*s*t*(t*y4 + s*(-3*y3 + y4)) +         &
!            &     mp**2*(s*t*(6*y3 - 16*y4) + 3*t**2*y4 - 8*s**2*(3*y3 + y4))        &
!            &     ) + mn**4*(mp**4*y2**2*(37*y3 - y4) +        &
!            &     2*mp**2*((-9 + 20*s*y2**2 - 2*t*y2**2)*y3 -         &
!            &     (7 + 8*t*y2**2)*y4) +         &
!            &     4*(6*s**2*y2**2*y3 + t*(y3 - 2*y4) +                                 &
!            &     s*((-4 + 3*t*y2**2)*y3 - 4*(y4 + t*y2**2*y4)))) +                                 &
!            &     mn**2*(2*mp**6*y2**2*(9*y3 + 7*y4) +                                 &
!            &     mp**4*((-24 + 22*s*y2**2 - t*y2**2)*y3 +                                 &
!            &     (16 + 2*s*y2**2 - 7*t*y2**2)*y4) -                                 &
!            &     2*mp**2*(4*s**2*y2**2*y3 - t*(y3 + (3 + 2*t*y2**2)*y4) +                                 &
!            &     s*((6 - 4*t*y2**2)*y3 + 2*(5 + 4*t*y2**2)*y4)) -                                 &
!            &     4*s*(4*s**2*y2**2*y3 - 2*t*(-y3 + y4 + t*y2**2*y4) +                                 &
!            &     s*((-2 + 3*t*y2**2)*y3 - 2*(y4 + t*y2**2*y4))))))
!
!       part2=2*y4*(2*fn**2*(mp**6*y2**2*(y3 - y4) +                                 &
!            &     2*(s + t)*(s + s*t*y2**2 + t**2*y2**2)*y4 +                                 &
!            &     (0,6)*mn**5*y2*(y3 + y4) + 4*mn**6*y2**2*(y3 + y4) -                                 &
!            &     mn**4*((7*mp**2*y2**2 + 4*(-2 + 2*s*y2**2 + t*y2**2))*y3 +                                 &
!            &     (-14 + 9*mp**2*y2**2 + 8*s*y2**2 + 6*t*y2**2)*y4) -                                 &
!            &     (0,6)*mn**3*y2*(mp**2*(2*y3 - y4) + 2*s*(y3 + y4) +                                 &
!            &     t*(y3 + 3*y4)) +                                 &
!            &     mp**4*(2*y4 + 2*s*y2**2*(-y3 + y4) + t*y2**2*(y3 + 8*y4)) -                                &
!            &     mp**2*(8*t**2*y2**2*y4 +                                 &
!            &     t*((4 - s*y2**2)*y3 + 9*s*y2**2*y4) +                                 &
!            &     s*(4*y4 + s*y2**2*(-y3 + y4))) +                                 &
!            &     mn**2*(4*mp**4*y2**2*(2*y3 - 3*y4) +                                 &
!            &     t*(-2*y3 + 4*s*y2**2*y3 - 6*y4 + 4*s*y2**2*y4) +                                 &
!            &     mp**2*((4 - 10*s*y2**2 + 3*t*y2**2)*y3 +                                 &
!            &     (20 - 6*s*y2**2 + t*y2**2)*y4) +                                 &
!            &     4*s*(-2*y4 + s*y2**2*(y3 + y4))) +                                 &
!            &     (0,1)*mn*y2*(-(t**2*(y3 - 12*y4)) + 6*s**2*(y3 + y4) +                                 &
!            &     6*s*t*(y3 + 3*y4) + 6*mp**4*(y3 + 4*y4) +                                 &
!            &     mp**2*(t*(13*y3 - 42*y4) - 6*s*(2*y3 + 5*y4)))) +                                 &
!            &     y1**2*(-(mp**8*y2**2*(y3 - 17*y4)) + 4*mn**8*y2**2*y4 +                                 &
!            &     4*(s + t)**2*(t + s**2*y2**2 + s*t*y2**2)*y4 +                                 &
!            &     4*mn**6*(2*(1 + t*y2**2)*y3 -                                 &
!            &     (-2 + 2*mp**2*y2**2 + 4*s*y2**2 + 3*t*y2**2)*y4) +                                 &
!            &     2*mp**6*((1 + s*y2**2 - 2*t*y2**2)*y3 -                                 &
!            &     (1 + 25*s*y2**2 + 16*t*y2**2)*y4) -                                 &
!            &     (0,4)*mn**5*y2*(-(t*(y3 - 3*y4)) + 2*mp**2*(y3 + 3*y4)) -                                 &
!            &     2*mp**2*(s + t)*(12*s**2*y2**2*y4 + 2*t*(4 + t*y2**2)*y4 +                                 &
!            &     s*(-y3 + y4 + 14*t*y2**2*y4)) +                                 &
!            &     mp**4*(-(s**2*y2**2*(y3 - 53*y4)) +                                 &
!            &     2*t*(y3 + 2*(4 + 5*t*y2**2)*y4) +                                 &
!            &     s*(-((4 + t*y2**2)*y3) + (4 + 73*t*y2**2)*y4)) +                                 &
!            &     (0,2)*mn**3*y2*(2*mp**4*(4*y3 + 3*y4) +                                 &
!            &     2*t*(-2*s*y3 + t*y3 + 6*s*y4 + 6*t*y4) +                                 &
!            &     mp**2*(-8*t*y3 + 3*t*y4 + 8*s*(y3 + 3*y4))) -                                 &
!            &     (0,1)*mn*y2*(4*mp**6*(2*y3 + 15*y4) +                                 &
!            &     4*t*(s + t)*(-(s*y3) + 3*s*y4 + 3*t*y4) -                                 &
!            &     4*mp**4*(4*s*y3 - 4*t*y3 + 21*s*y4 + 12*t*y4) +                                 &
!            &     mp**2*(8*s**2*(y3 + 3*y4) - 3*t**2*(y3 + 6*y4) +                                 &
!            &     s*t*(16*y3 + 6*y4))) +                                 &
!            &     mn**2*(2*mp**6*y2**2*(7*y3 + y4) +                                 &
!            &     mp**4*((16 + 2*s*y2**2 - 7*t*y2**2)*y3 +                                 &
!            &     (-24 - 10*s*y2**2 + 23*t*y2**2)*y4) +                                 &
!            &     2*mp**2*(20*s**2*y2**2*y4 + t*(3*y3 + 2*t*y2**2*y3 + y4) -                                 &
!            &     2*s*(5*y3 + 4*t*y2**2*y3 + 3*y4 - 10*t*y2**2*y4)) -                                 &
!            &     4*(s + t)*(4*s**2*y2**2*y4 + t**2*y2**2*y4 -                                 &
!            &     s*(2*y3 + 2*t*y2**2*y3 + 2*y4 - 5*t*y2**2*y4))) -                                 &
!            &     mn**4*(mp**4*y2**2*(y3 - 21*y4) +                                 &
!            &     2*mp**2*((7 + 8*t*y2**2)*y3 +                                 &
!            &     (9 + 4*s*y2**2 - 6*t*y2**2)*y4) -                                 &
!            &     4*(6*s**2*y2**2*y4 -                                 &
!            &     s*(4*y3 + 4*t*y2**2*y3 + 4*y4 - 9*t*y2**2*y4) +                                 &
!            &     t*(-2*y3 - 3*y4 + 3*t*y2**2*y4)))) +                                 &
!            &     fn*y1*(-18*mn**6*y2*y4 -                                 &
!            &     (0,4)*mn**5*((-3 + 2*mp**2*y2**2 - t*y2**2)*y3 +                                 &
!            &     3*(-1 + 2*mp**2*y2**2 + t*y2**2)*y4) +                                 &
!            &     mn**4*y2*(10*mp**2*y3 + 54*s*y4 + t*(y3 + 18*y4)) -                                 &
!            &     mn**2*y2*(54*s**2*y4 + 3*t**2*(5*y3 + 6*y4) +                                 &
!            &     10*mp**4*(2*y3 + 9*y4) + 2*s*t*(y3 + 36*y4) +                                 &
!            &     4*mp**2*(5*s*y3 - 12*t*y3 - 18*s*y4)) +                                 &
!            &     y2*(2*mp**6*(5*y3 - 18*y4) +                                 &
!            &     mp**4*(-20*s*y3 + 11*t*y3 + 90*s*y4 + 90*t*y4) +                                 &
!            &     2*mp**2*(s**2*(5*y3 - 36*y4) + 4*s*t*(y3 - 18*y4) -                                 &
!            &     36*t**2*y4) +                                 &
!            &     (s + t)*(18*s**2*y4 + 18*t**2*y4 + s*t*(y3 + 36*y4))) +                                  &
!            &     (0,2)*mn**3*(2*mp**4*y2**2*(4*y3 + 3*y4) +                                 &
!            &     mp**2*(4*(-3 + 2*s*y2**2 - 2*t*y2**2)*y3 +                                 &
!            &     3*(2 + 8*s*y2**2 + t*y2**2)*y4) -                                 &
!            &     2*(t*(3*y3 - t*y2**2*y3 + 9*y4 - 6*t*y2**2*y4) +                                 &
!            &     2*s*(3*y3 + t*y2**2*y3 + 3*y4 - 3*t*y2**2*y4))) -                                 &
!            &     (0,1)*mn*(4*mp**6*y2**2*(2*y3 + 15*y4) -                                 &
!            &     4*mp**4*((3 + 4*s*y2**2 - 4*t*y2**2)*y3 +                                 &
!            &     3*(4 + 7*s*y2**2 + 4*t*y2**2)*y4) -                                 &
!            &     2*(2*s**2*(3*y3 + t*y2**2*y3 + 3*y4 - 3*t*y2**2*y4) +                                 &
!            &     2*s*t*((3 + t*y2**2)*y3 + 3*(3 - 2*t*y2**2)*y4) -                                 &
!            &     t**2*(y3 + 6*(-2 + t*y2**2)*y4)) +                                 &
!            &     mp**2*(8*s**2*y2**2*(y3 + 3*y4) +                                 &
!            &     s*(8*(3 + 2*t*y2**2)*y3 + 6*(10 + t*y2**2)*y4) -                                 &
!            &     t*((26 + 3*t*y2**2)*y3 + 6*(-14 + 3*t*y2**2)*y4)))))
!
!       vsum1=part1+part2
!
!       return
!     end function vsum1
!
!     !*******************************************************************************
!
!     complex function vsum2(s,t,u,y1,y2,y3,y4,x1,x3,qn)
!       use constants, only: mN
!       integer qn
!       real s,t,u,fn
!       real, parameter :: mp = 1.02
!       complex y1,y2,y3,y4
!       complex x1,x3,x1d,x3d!,y1d,y2d,y3d,y4d
!       complex part1,part2,part3,part4!,part5
!
!       fn=real(qn)
!       x1d=conjg(x1)
!       x3d=conjg(x3)
!       !y1d=conjg(y1)
!       !y2d=conjg(y2)
!       !y3d=conjg(y3)
!       !y4d=conjg(y4)
!
!       part1=2*t*x3*(fn*((0,-1)*mn*(4*mp**2 - t)*(y3 + y4) +                                 &
!            &     mn**4*y2*(y3 + y4) -                                 &
!            &     mn**2*y2*(4*mp**2*y3 + t*(-y3 + y4) + 2*s*(y3 + y4)) +                                &
!            &     y2*(-(mp**4*(y3 - 3*y4)) + 2*s*t*y4 + t**2*y4 -                                 &
!            &     4*mp**2*(s + t)*y4 + s**2*(y3 + y4))) +                                 &
!            &     y1*(s**2*y3 - mp**4*(y3 - 3*y4) - (0,1)*mn**3*t*y2*(y3 - y4)+                                 &
!            &     s**2*y4 + 2*s*t*y4 + t**2*y4 - 4*mp**2*(s + t)*y4 +                                 &
!            &     mn**4*(y3 + y4) -                                 &
!            &     mn**2*(4*mp**2*y3 + t*(-y3 + y4) + 2*s*(y3 + y4)) +                                &
!            &     (0,1)*mn*y2*(2*mp**2*t*y4 + 4*mp**4*(y3 + y4) -                                 &
!            &     t*(t*y4 + s*(-y3 + y4)))))
!
!       part2=-2*t*x3d*(fn*((0,-1)*mn*(4*mp**2 - t)*(y3 + y4) +                                 &
!            &     mn**4*y2*(y3 + y4) -                                 &
!            &     mn**2*y2*(4*mp**2*y3 + t*(-y3 + y4) + 2*s*(y3 + y4)) +                                &
!            &     y2*(-(mp**4*(y3 - 3*y4)) + 2*s*t*y4 + t**2*y4 -                                 &
!            &     4*mp**2*(s + t)*y4 + s**2*(y3 + y4))) +                                 &
!            &     y1*(s**2*y3 - mp**4*(y3 - 3*y4) - (0,1)*mn**3*t*y2*(y3 - y4)+                                 &
!            &     s**2*y4 + 2*s*t*y4 + t**2*y4 - 4*mp**2*(s + t)*y4 +                                 &
!            &     mn**4*(y3 + y4) -                                 &
!            &     mn**2*(4*mp**2*y3 + t*(-y3 + y4) + 2*s*(y3 + y4)) +                                &
!            &     (0,1)*mn*y2*(2*mp**2*t*y4 + 4*mp**4*(y3 + y4) -                                 &
!            &     t*(t*y4 + s*(-y3 + y4)))))
!
!       part3=(-2*x1d*(fn*(-2*mp**6*(y3 - 2*y4) +                                 &
!            &     2*mp**2*(3*s**2 + 4*s*t + 2*t**2)*y4 + mn**6*(y3 + y4) +                                 &
!            &     (0,3)*mn**3*mp**2*(2*mp**2 - t)*y2*(y3 + y4) -                                 &
!            &     (s**2 + s*t + t**2)*(t*y4 + s*(y3 + y4)) +                                 &
!            &     mn**4*(t*y3 - 3*s*(y3 + y4) + mp**2*(-4*y3 + 2*y4)) +                                 &
!            &     mp**4*(3*s*(y3 - 3*y4) - t*(2*y3 + 3*y4)) +                                 &
!            &     mn**2*(t**2*y3 + 4*mp**2*s*(y3 - 2*y4) + 2*s*t*y4 +                                 &
!            &     3*s**2*(y3 + y4) + mp**4*(-7*y3 + 5*y4)) +                                 &
!            &     (0,1)*mn*mp**2*y2*                                &
!            &     (6*mp**4*(y3 + y4) -                                 &
!            &     2*mp**2*(-2*t*(y3 - 4*y4) + 3*s*(y3 + y4)) +                                 &
!            &     t*(-(t*(y3 - 4*y4)) + 3*s*(y3 + y4)))) +                                 &
!            &     y1*((0,1)*mn**3*(6*mp**4 - 3*mp**2*t - t**2)*(y3 + y4) -                                 &
!            &     2*mn**6*mp**2*y2*(y3 + y4) +                                 &
!            &     mn**4*mp**2*y2*                                &
!            &     (mp**2*(5*y3 - y4) + t*(-3*y3 + y4) + 6*s*(y3 + y4)) +                                  &
!            &     mn**2*mp**2*y2*                                &
!            &     (2*s*t*(y3 - 3*y4) + 4*mp**4*(2*y3 - y4) +                                 &
!            &     t**2*(-y3 + y4) - 6*s**2*(y3 + y4) +                                 &
!            &     mp**2*(-2*s*y3 + 3*t*y3 + 10*s*y4 - 7*t*y4)) +                                 &
!            &     (0,1)*mn*(6*mp**6*(y3 + y4) + t**2*(t*y4 + s*(y3 + y4)) -                                 &
!            &     2*mp**4*(-2*t*(y3 - 4*y4) + 3*s*(y3 + y4)) +                                 &
!            &     mp**2*t*(t*y4 + 3*s*(y3 + y4))) +                                 &
!            &     mp**2*y2*(mp**6*(y3 - 5*y4) +                                 &
!            &     (2*s + t)*(2*s*t*y4 + t**2*y4 + s**2*(y3 + y4)) -                                 &
!            &     mp**2*(-(s*t*(y3 - 13*y4)) + 4*t**2*y4 +                                 &
!            &     3*s**2*(y3 + 3*y4)) + mp**4*(12*s*y4 + t*(y3 + 5*y4)))                                &
!            &     )))/mp**2
!
!       part4=(-2*x1*(fn*(-2*mp**6*(y3 - 2*y4) +                                 &
!            &     2*mp**2*(3*s**2 + 4*s*t + 2*t**2)*y4 + mn**6*(y3 + y4) +                                 &
!            &     (0,3)*mn**3*mp**2*(2*mp**2 - t)*y2*(y3 + y4) -                                 &
!            &     (s**2 + s*t + t**2)*(t*y4 + s*(y3 + y4)) +                                 &
!            &     mn**4*(t*y3 - 3*s*(y3 + y4) + mp**2*(-4*y3 + 2*y4)) +                                 &
!            &     mp**4*(3*s*(y3 - 3*y4) - t*(2*y3 + 3*y4)) +                                 &
!            &     mn**2*(t**2*y3 + 4*mp**2*s*(y3 - 2*y4) + 2*s*t*y4 +                                 &
!            &     3*s**2*(y3 + y4) + mp**4*(-7*y3 + 5*y4)) +                                 &
!            &     (0,1)*mn*mp**2*y2*                                &
!            &     (6*mp**4*(y3 + y4) -                                 &
!            &     2*mp**2*(-2*t*(y3 - 4*y4) + 3*s*(y3 + y4)) +                                 &
!            &     t*(-(t*(y3 - 4*y4)) + 3*s*(y3 + y4)))) +                                 &
!            &     y1*((0,1)*mn**3*(6*mp**4 - 3*mp**2*t - t**2)*(y3 + y4) -                                 &
!            &     2*mn**6*mp**2*y2*(y3 + y4) +                                 &
!            &     mn**4*mp**2*y2*                                &
!            &     (mp**2*(5*y3 - y4) + t*(-3*y3 + y4) + 6*s*(y3 + y4)) +                                  &
!            &     mn**2*mp**2*y2*                                &
!            &     (2*s*t*(y3 - 3*y4) + 4*mp**4*(2*y3 - y4) +                                 &
!            &     t**2*(-y3 + y4) - 6*s**2*(y3 + y4) +                                 &
!            &     mp**2*(-2*s*y3 + 3*t*y3 + 10*s*y4 - 7*t*y4)) +                                 &
!            &     (0,1)*mn*(6*mp**6*(y3 + y4) + t**2*(t*y4 + s*(y3 + y4)) -                                 &
!            &     2*mp**4*(-2*t*(y3 - 4*y4) + 3*s*(y3 + y4)) +                                 &
!            &     mp**2*t*(t*y4 + 3*s*(y3 + y4))) +                                 &
!            &     mp**2*y2*(mp**6*(y3 - 5*y4) +                                 &
!            &     (2*s + t)*(2*s*t*y4 + t**2*y4 + s**2*(y3 + y4)) -                                 &
!            &     mp**2*(-(s*t*(y3 - 13*y4)) + 4*t**2*y4 +                                 &
!            &     3*s**2*(y3 + 3*y4)) + mp**4*(12*s*y4 + t*(y3 + 5*y4)))                                &
!            &     )))/mp**2
!
!       vsum2=part1+part2+part3+part4
!
!       return
!     end function vsum2
!
!   end function stotphi





  !****************************************************************************
  !****s* photonXSections/calcXS_gammaN2strange
  ! NAME
  ! subroutine calcXS_gammaN2strange(srts,sigma)
  ! PURPOSE
  ! Produces cross section for
  !     gamma N -> Lambda K, Sigma K, N K Kbar
  ! OUTPUT
  ! * real, dimension(1:3 :: sigma ! cross section
  ! * sigma(1) :: Lambda K
  ! * sigma(2) :: Sigma K
  ! * sigma(3) :: N K Kbar
  ! NOTES
  ! UNITS of cross sections????
  !****************************************************************************
  subroutine calcXS_gammaN2strange(srts,sigma)

      use IdTable, only: lambda,SigmaResonance
      use particleProperties, only: hadron
      use threeBodyPhaseSpace, only: Integrate_3bodyPS
      use constants, only: mN, mK

      real, intent(in) ::  srts
      real, intent(out),dimension(1:3) :: sigma

      real m2(3),lambdaField(3)
      real mass1,mass2,mass3
      real pfinal,pini,srts0,ps
      integer i

      data (m2(i),i=1,3) /13., 15.,12./
      data (lambdaField(i),i=1,3) /0.5, 0.4, 0.7/

      do i=1,3
         sigma(i)=0.
      end do
      pini=(srts**2-mN**2)/2./srts

      !***** Lambda K
      mass1=mK
      mass2=hadron(Lambda)%mass
      srts0=mass1+mass2

      if (srts.gt.mass1+mass2) then
         pfinal=sqrt((srts**2-(mass1+mass2)**2)*  (srts**2-(mass1-mass2)**2)/4./srts**2)
         sigma(1)=m2(1)*pfinal/pini/srts**2*lambdaField(1)**2/ (lambdaField(1)**2+(srts-srts0)**2)
      else
         sigma(1)=0.
      end if

      !***** Sigma K
      mass1=mK
      mass2=hadron(SigmaResonance)%mass
      !      mass1=mesprop1(10,1)
      !      mass2=sbarprop1(2,1)
      srts0=mass1+mass2
      if (srts.gt.mass1+mass2) then
         pfinal=sqrt((srts**2-(mass1+mass2)**2)*  (srts**2-(mass1-mass2)**2)/4./srts**2)
         sigma(2)=m2(2)*pfinal/pini/srts**2*lambdaField(2)**2/ (lambdaField(2)**2+(srts-srts0)**2)
      else
         sigma(2)=0.
      end if
      ! use the same cross section for both isospin channels
      sigma(2)=2.*sigma(2)

      !*****N K+ K-

      mass1=mK
      mass2=mK
      mass3=mN
      srts0=mass1+mass2+mass3
      if (srts.gt.mass1+mass2+mass3) then
         !   call bops3(ps,srts,mass1,mass2,mass3)
         ps = Integrate_3bodyPS (srts, mass1, mass2, mass3)
         sigma(3)=m2(3)*ps/srts/pini*lambdaField(3)**2/ (lambdaField(3)**2+(srts-srts0)**2)
      else
         sigma(3)=0.
      end if
      ! use the same cross section for all 3 isospin channels
      sigma(3)=3.*sigma(3)
      return

  end subroutine calcXS_gammaN2strange



  !****************************************************************************
  !****s* photonXSections/calcXS_omega_saphir
  ! NAME
  ! subroutine calcXS_omega_saphir(srts,cs)
  ! PURPOSE
  ! Calculates the cross section for gamma + N -> omega + N (fit to SAPHIR data).
  ! See P. Muehlich, diss., chapter 9.3.3.
  ! The data can be found in J.Barth, Eur. Phys. J. A18 (2003) 117-127
  ! INPUTS
  ! * real,intent(in):: srts --- sqrt(s) for gamma + nucleon (in GeV)
  ! OUTPUT
  ! * real,intent(out):: cs --- cross section in microbarn
  ! NOTES
  ! The matrix element for the process is tabulated in the file
  ! "gammaN_omegaN_ME_saphir.dat". It is computed from the (splined)
  ! SAPHIR cross sections.
  !****************************************************************************
  subroutine calcXS_omega_saphir(srts,cs)
    use IdTable, only: omegaMeson
    use particleProperties, only: hadron
    use quadpack, only: qags
    use constants, only: mN

    real,intent(in):: srts
    real,intent(out):: cs
    real:: a,b,p_i
    real,   parameter :: qag_absError=0.      ! absolute error=0  -> relative error gives the accuracy
    real,   parameter :: qag_relError=0.001
    real              :: qag_error_estimation=0.0
    integer           :: qag_neval=0,qag_error_code=0
    logical,save :: first = .true.

    if (first) then
      m_omega=hadron(omegaMeson)%mass
      call read_data()
      first=.false.
    end if

    sqrts=srts

    p_i=(srts**2-mN**2)/(2*srts)
    cs=0.
    a=0.
    b=srts-mN
    qag_error_estimation=0.0
    qag_neval=0
    qag_error_code=0

    call qags(integrand_saphir,a,b,qag_absError,qag_relError,cs,qag_error_estimation,qag_neval,qag_error_code)
    if (qag_error_code.ne.0 .and. qag_error_estimation/cs>qag_relError*3.) then
      write(*,*) 'Error with QAGS in calcXS_omega_saphir:',qag_error_code
      write(*,*) 'result: ',cs,' +- ',qag_error_estimation,' (',qag_error_estimation/cs*100.,'%)'
      write(*,*) 'neval=',qag_neval
      write(*,*) 'srts=',srts
    end if

    cs=cs/(p_i*srts**2)

  contains

    subroutine read_data()
      use inputGeneral, only: path_To_Input
      integer :: i
      open(202,file=trim(path_To_Input)//"/gammaN_omegaN_ME_saphir.dat")
      do i=lbound(y,dim=1),ubound(y,dim=1)
        read (202,*) x(i),y(i)
      end do
      close(202)
    end subroutine

  end subroutine

  real function integrand_saphir(m)
    use constants, only: pi, mN
    use IdTable, only: omegaMeson
    use mesonWidthMedium, only: WidthMesonMedium
    use spline, only: bsplint2
    real,intent(in):: m ! invariant mass of vector meson
    real :: p_f,Q,ME,w
    p_f = sqrt((sqrts**2-(mN+m+spot_omega)**2)*(sqrts**2-(mN-m-spot_omega)**2))/(2.*sqrts)
    Q = sqrts + m_omega - m
    if (Q>1.72) then
      ! get matrix element
      ME = bsplint2(x,y,Q)
    else
      ME = 0.
    end if
    integrand_saphir = p_f * ME
    ! spectral function with in-medium width (ATTENTION: neglecting momentum dependence)
    w = WidthMesonMedium(omegaMeson,m,(/0.,0.,0.,0./),med)
    integrand_saphir = integrand_saphir * 2./pi*m**2*w / ((m**2-(m_omega+spot_omega)**2)**2+m**2*w**2)
  end function


end module photonXSections
