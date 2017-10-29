!******************************************************************************
!****m* /mesonMeson_Xsections
! NAME
! module mesonMeson_Xsections
! PURPOSE
! Provide parametrizations for (some) meson-meson cross sections
!******************************************************************************
module mesonMeson_Xsections

  use IdTable

  implicit none
  private

  public :: sig_pipi2kkbar
  public :: kkbar_cross
  public :: kkbar_out
  public :: kstarkbar_cross
  public :: mesMes_Tabulate

  logical, save :: debug = .false.

  integer, parameter, dimension(2,21) :: channelID = reshape((/ &
       pion, pion, & ! 1
       rho, rho, & ! 2
       pion, rho, & ! 3
       pion, eta, & ! 4
       pion, sigmaMeson, & ! 5
       pion, omegaMeson, & ! 6
       pion, etaPrime, & ! 7
       eta, eta, & ! 8
       eta, rho, & ! 9
       eta, sigmaMeson, & ! 10
       eta, omegaMeson, & ! 11
       eta, etaPrime, & ! 12
       rho, sigmaMeson, & ! 13
       rho, omegaMeson, & ! 14
       rho, etaPrime, & ! 15
       sigmaMeson, sigmaMeson, & ! 16
       sigmaMeson, omegaMeson, & ! 17
       sigmaMeson, etaPrime, & ! 18
       omegaMeson, omegaMeson, & ! 19
       omegaMeson, etaPrime, & ! 20
       etaPrime, etaPrime /), & ! 21
       (/2,21/))

  integer, parameter :: nSrts = 500
  real, parameter :: dSrts = 0.05
  real, dimension(21,0:nSrts) :: arrPCM
  real, dimension(21) :: arrSrts0 ! the thresholds
  real, dimension(21) :: arrFakInt ! the spectral function norms

contains

  !****************************************************************************
  !****f* mesonMeson_Xsections/sig_pipi2kkbar
  ! NAME
  ! real function sig_pipi2kkbar(srts,srts0)
  ! PURPOSE
  ! Cross section of pi pi --> K Kbar averaged over isospins of
  ! incoming pions and summed over isospins of outgoing pions
  !
  ! INPUTS
  ! * real:: srts  -- c.m. energy of incoming particles (GeV),
  ! * real:: srts0 -- threshold c.m. energy (GeV),
  ! OUTPUT
  ! * sig_pipi2kkbar -- cross section (mbarn)
  !****************************************************************************

  real function sig_pipi2kkbar(srts,srts0)
    real, intent(in) :: srts,srts0
    if (srts.le.srts0) then
       sig_pipi2kkbar = 0.
       return
    end if
    !   Parameterization from W. Cassing et al., NPA 614 (1997) 415:
    sig_pipi2kkbar = 2.7*(1.-(srts0/srts)**2)**0.76
  end function sig_pipi2kkbar


  !****************************************************************************
  !****s* mesonMeson_Xsections/kkbar_cross
  ! NAME
  ! subroutine kkbar_cross(iQ,srts,pinitial2,const,msigbg,useWidth)
  ! PURPOSE
  ! calculate the cross-sections for the KKbar incoming channel
  !
  ! INPUTS
  ! * integer:: iQ -- total charge of K and Kbar
  ! * real   :: srts -- c.m. energy (GeV)
  ! * real   :: pinitial2 -- c.m. momentum squared of incoming particles (GeV/c)**2
  ! * real   :: const -- cross section of nonstrange+nonstrange --> K Kbar
  !   or --> Kstar Kbar, K Kbarstar (mbarn)
  ! * logical:: useWidth -- flag whether to use the width or just pole masses
  ! OUTPUT
  ! * real, dimension(21) :: msigbg(i), i=1,2,..,21 --- partial cross sections for different outgoing channels (mbarn)
  !****************************************************************************
  subroutine kkbar_cross(iQ,srts,pinitial2,const,msigbg,useWidth)

    use constants, only: mK

    integer, intent(in) :: iQ
    real, intent(in) :: srts,pinitial2,const
    real, dimension(21), intent(out) :: msigbg
    logical, intent(in) :: useWidth


    ! This collects the factors iso, spin, and 1/2 for identical particles
    ! this is done for iQ = 0, if-statements may change this
    real, parameter, dimension(21) :: fak = (/ &
       0.5 * 1.0 * 1.0, &  ! 01: pion, pion
       0.5 * 9.0 * 1.0, &  ! 02: rho, rho
       0.0, &              ! 03: pion, rho
       0.5 * 1.0 * 1.0, &  ! 04: pion, eta  --> iQ==0: 0.5, else: 1.0
       0.0, &              ! 05: pion, sigmaMeson
       0.0, &              ! 06: pion, omegaMeson
       0.5 * 1.0 * 1.0, &  ! 07: pion, etaPrime --> iQ==0: 0.5, else: 1.0
       0.5 * 1.0 * 0.5, &  ! 08: eta, eta --> iQ==0: 0.5, else: 0.0 !!!
       0.0, &              ! 09: eta, rho
       0.0, &              ! 10: eta, sigmaMeson
       0.0, &              ! 11: eta, omegaMeson
       0.5 * 1.0 * 1.0, &  ! 12: eta, etaPrime --> iQ==0: 0.5, else: 0.0 !!!
       0.0, &              ! 13: rho, sigmaMeson
       0.5 * 9.0 * 1.0, &  ! 14: rho, omegaMeson --> iQ==0: 0.5, else: 1.0
       0.0, &              ! 15: rho, etaPrime
       0.5 * 1.0 * 0.5, &  ! 16: sigmaMeson, sigmaMeson --> iQ==0: 0.5, else: 0.0 !!!
       0.0, &              ! 17: sigmaMeson, omegaMeson
       0.0, &              ! 18: sigmaMeson, etaPrime
       0.5 * 9.0 * 0.5, &  ! 19: omegaMeson, omegaMeson --> iQ==0: 9.0, else: 0.0 !!!
       0.0, &              ! 20: omegaMeson, etaPrime
       0.5 * 1.0 * 0.5 /)  ! 21: etaPrime, etaPrime --> iQ==0: 0.5, else: 0.0 !!!



    real :: srts0,pfinal2

    if (debug) write(*,*) 'In kkbar_cross:', srts

    !     fritiof might produce lower mass kaons
    if (srts.lt.2.*mK) then
       srts0 = 2.*0.493
    else
       srts0 = 2.*mK
    end if

    if (srts0.gt.srts) then
       write(*,*) 'srts0 gt srts in mesmes',srts0,srts
       stop
    end if

    ! set all processes to zero:
    msigbg = 0.0

    ! === KKbar -> pipi:
    pfinal2 = pcm2( srts, 1, useWidth )
    msigbg(1) = fak(1)*9./4.*sig_pipi2kkbar(srts,srts0)*pfinal2/pinitial2

    ! === KKbar -> rhorho:
    ! (same isofac as before but additional spinfactor)
    pfinal2 = pcm2( srts, 2, useWidth )
    msigbg(2) = fak(2)*9./4.*sig_pipi2kkbar(srts,srts0)*pfinal2/pinitial2

    ! === KKbar -> pi rho:
    ! (this is possible only in p-wave, hence for simplicity put to zero)
    !    msigbg(3) = 0.

    ! === KKbar -> pi eta:
    pfinal2 = pcm2( srts, 4, useWidth )
    msigbg(4) = fak(4)*const*pfinal2/pinitial2

    ! === KKbar -> pi sigma, because of parity:
    !    msigbg(5) = 0.

    ! === KKbar -> pi omega, because of p-wave:
    !    msigbg(6) = 0.

    ! === KKbar -> pi etap:
    pfinal2 = pcm2( srts, 7, useWidth )
    msigbg(7) = fak(7)*const*pfinal2/pinitial2

    ! === KKbar -> eta eta:
    pfinal2 = pcm2( srts, 8, useWidth )
    msigbg(8) = fak(8)*const*pfinal2/pinitial2

    ! === KKbar -> eta rho (p-wave):
    !    msigbg(9) = 0.

    ! === KKbar -> eta sigma, because of parity:
    !    msigbg(10) = 0.

    ! === KKbar -> eta omega (p-wave):
    !    msigbg(11) = 0.

    ! === KKbar -> eta etap:
    pfinal2 = pcm2( srts, 12, useWidth )
    msigbg(12) = fak(12)*const*pfinal2/pinitial2

    ! === KKbar -> rho sigma (p-wave):
    !    msigbg(13) = 0.

    ! === KKbar -> rho omega:
    pfinal2 = pcm2( srts, 14, useWidth )
    msigbg(14) = fak(14)*const*pfinal2/pinitial2

    ! === KKbar -> rho etap (p-wave):
    !    msigbg(15) = 0.

    ! === KKbar -> sigma sigma:
    pfinal2 = pcm2( srts, 16, useWidth )
    msigbg(16) = fak(16)*const*pfinal2/pinitial2

    ! === KKbar -> sigma omega (p-wave):
    !    msigbg(17) = 0.

    ! === KKbar -> sigma etap (parity):
    !    msigbg(18) = 0.

    ! === KKbar -> omega omega:
    pfinal2 = pcm2( srts, 19, useWidth )
    msigbg(19) = fak(19)*const*pfinal2/pinitial2

    ! === KKbar -> omega etap (p-wave):
    !    msigbg(20) = 0.

    ! === KKbar -> etap etap:
    pfinal2 = pcm2( srts, 21, useWidth )
    msigbg(21) = fak(21)*const*pfinal2/pinitial2


    ! now we are doing a charge-correction to many channels.

    if (iQ.ne.0) then
       msigbg( 4) = msigbg( 4) * 2
       msigbg( 7) = msigbg( 7) * 2
       msigbg( 8) = 0.0
       msigbg(12) = 0.0
       msigbg(14) = msigbg(14) * 2
       msigbg(16) = 0.0
       msigbg(19) = 0.0
       msigbg(21) = 0.0
    end if

  end subroutine kkbar_cross


  !****************************************************************************
  !****s* mesonMeson_Xsections/kstarkbar_cross
  ! NAME
  ! subroutine kstarkbar_cross(iQ,srts,pinitial2,const,msigbg,useWidth)
  ! PURPOSE
  ! calculate cross-sections for Kstar Kbar or K Kstarbar incoming channels
  !
  ! INPUTS
  ! * integer :: iQ --- total charge of Kstar and Kbar,
  ! * real :: srts --- c.m. energy (GeV),
  ! * real :: pinitial2 --- c.m. momentum squared of incoming particles (GeV/c)**2,
  ! * real :: const --- cross section of nonstrange+nonstrange --> K Kbar
  !   or --> Kstar Kbar, K Kbarstar (mbarn),
  ! * logical:: useWidth -- flag whether to use the width or just pole masses
  ! OUTPUT
  ! * real, dimension(:) :: msigbg(i), i=1,2,..,21 --- partial cross sections
  !   for different outgoing channels (mbarn)
  !****************************************************************************
  subroutine kstarkbar_cross(iQ,srts,pinitial2,const,msigbg,useWidth)

    use constants, only: mK
    use IdTable, only: kaonStar
    use particleProperties, only: hadron

    integer, intent(in) :: iQ
    real, intent(in) :: srts,pinitial2,const
    real, dimension(21), intent(out) :: msigbg
    logical, intent(in) :: useWidth

    real :: srts0,pfinal2

    ! This collects the factors iso, spin, and 1/2 for identical particles
    ! this is done for iQ = 0, if-statements may change this
    real, parameter, dimension(21) :: fak = (/ &
       0.0, &               ! 01: pion, pion
       0.0, &               ! 02: rho, rho
       1.0 * 1.0 * 1.0, &   ! 03: pion, rho
       0.0, &               ! 04: pion, eta
       0.0, &               ! 05: pion, sigmaMeson
       0.25 * 1.0 * 1.0, &  ! 06: pion, omegaMeson --> iQ==0: 0.25, else: 0.5
       0.0, &               ! 07: pion, etaPrime
       0.0, &               ! 08: eta, eta
       0.25 * 1.0 * 1.0, &  ! 09: eta, rho --> iQ==0: 0.25, else: 0.5
       0.0, &               ! 10: eta, sigmaMeson
       0.25 * 1.0 * 1.0, &  ! 11: eta, omegaMeson --> iQ==0: 0.25, else: 0.0 !!!!
       0.0, &               ! 12: eta, etaPrime
       0.0, &               ! 13: rho, sigmaMeson
       0.0, &               ! 14: rho, omegaMeson
       0.25 * 1.0 * 1.0, &  ! 15: rho, etaPrime --> iQ==0: 0.25, else: 0.5
       0.0, &               ! 16: sigmaMeson, sigmaMeson
       0.0, &               ! 17: sigmaMeson, omegaMeson
       0.0, &               ! 18: sigmaMeson, etaPrime
       0.0, &               ! 19: omegaMeson, omegaMeson
       0.25 * 1.0 * 1.0, &  ! 20: omegaMeson, etaPrime --> iQ==0: 0.25, else: 0.0 !!!!
       0.0 /)               ! 21: etaPrime, etaPrime



    if (debug) write(*,*) 'In kstarkbar_cross:', srts

    srts0 = mK+hadron(kaonStar)%minmass

    if (srts0.gt.srts) then
       write(*,*) 'srts0 gt srts in mesmes',srts0,srts
       stop
    end if

    ! set all processes to zero:
    msigbg = 0.0

    ! === Kbar Kstar-> pi rho:
    pfinal2 = pcm2( srts,  3, useWidth )
    msigbg(3) = fak(3)*9./4.*sig_pipi2kkbar(srts,srts0)*pfinal2/pinitial2

    ! === Kbar Kstar-> pi omega:
    pfinal2 = pcm2( srts,  6, useWidth )
    msigbg(6) = fak(6)*const*pfinal2/pinitial2

    ! === Kbar Kstar-> eta rho:
    pfinal2 = pcm2( srts,  9, useWidth )
    msigbg(9) = fak(9)*const*pfinal2/pinitial2

    ! === Kbar Kstar-> eta omega:
    pfinal2 = pcm2( srts, 11, useWidth )
    msigbg(11) = fak(11)*const*pfinal2/pinitial2

    ! === Kbar Kstar-> rho sigma (parity):
    !    msigbg(13) = 0.

    ! === Kbar Kstar-> rho etap:
    pfinal2 = pcm2( srts, 15, useWidth )
    msigbg(15) = fak(15)*const*pfinal2/pinitial2

    ! === Kbar Kstar-> sigma omega (parity):
    !    msigbg(17) = 0.

    ! === Kbar Kstar-> etap omega:
    pfinal2 = pcm2( srts, 20, useWidth )
    msigbg(20) = fak(20)*const*pfinal2/pinitial2



    ! now we are doing a iQ-correction to many channels.
    ! please check, whether this has to be done and everything is justified!

    if (iQ.ne.0) then
       msigbg( 6) = msigbg( 6) * 2
       msigbg( 9) = msigbg( 9) * 2
       msigbg(11) = 0.0
       msigbg(15) = msigbg(15) * 2
       msigbg(20) = 0.0
    end if

    msigbg = msigbg * 2.0

  end subroutine kstarkbar_cross

  !****************************************************************************
  !****s* mesonMeson_Xsections/kkbar_out
  ! NAME
  ! subroutine kkbar_out(iQ,msigbg,sigbgt,teilchenOut)
  ! PURPOSE
  ! choose outgoing state for the K Kbar annihilation
  ! INPUTS
  ! * integer :: iQ --- total charge of K and Kbar,
  ! * real, dimension(:) ::msigbg(i), i=1,2,...,21 --- partial cross sections
  !   for different outgoing channels (mbarn),
  ! * real ::sigbgt --- total background cross section (mbarn)
  ! OUTPUT
  ! * type(preEvent),dimension(21) :: teilchenOut   ---   outgoing particles
  !
  ! NOTE
  ! * This routine handles also the kstarkbar case
  !****************************************************************************
  subroutine kkbar_out(iQ,msigbg,sigbgt,partOut)

    use random, only: rn
    use preEventDefinition

    integer, intent(in) :: iQ
    real, dimension(21), intent(in) :: msigbg
    real, intent(in) :: sigbgt
    type(preEvent),dimension(:), intent(out) :: partOut

    logical :: flag
    real :: msig,x
    integer :: mch,n1,n2

    !     determine outgoing channel
    flag = .true.
    mch = 0
    msig = 0.
    x = rn()*sigbgt
    do while (flag)
       mch = mch + 1
       if (mch.gt.21) then
          write(*,*) 'problems in mesmes (kkbar_out)'
          stop
       end if
       msig = msig + msigbg(mch)
       if (x.le.msig) flag = .false.
    end do

    partOut(1:2)%ID = channelID(1:2,mch)

    !     Determine charges:

    select case (mch)
    case (1,2,3)
       !       outgoing pipi, rhorho, pirho:
       !       in order to distribute charge equally : random decision
       if (rn().lt.0.5) then
          n1 = 1
          n2 = 2
       else
          n2 = 1
          n1 = 2
       end if

       if (iQ.eq.0) then
          if (rn().lt.5./6.) then
             partOut(n1)%charge = 1
             partOut(n2)%charge = -1
          else
             partOut(n1)%charge = 0
             partOut(n2)%charge = 0
          end if
       else
          partOut(n1)%charge = iQ
          partOut(n2)%charge = 0
       end if

    case (9)
       partOut(1)%charge = 0
       partOut(2)%charge = iQ

    case default
       partOut(1)%charge = iQ
       partOut(2)%charge = 0

    end select

  end subroutine kkbar_out

  !****************************************************************************
  !****if* mesonMeson_Xsections/pcm2Integrated
  ! NAME
  ! real function pcm2Integrated( srts, iCh )
  ! PURPOSE
  ! calculate the phase space (squared) by integrating over the spectral
  ! functions of the particles of the given channel.
  !
  ! Threshold is treated correctly by using the minimal masses.
  !
  ! INPUTS
  ! * real :: srts -- c.m. energy (GeV)
  ! * integer :: iCh -- number of channel to consider
  ! OUTPUT
  ! * phase space squared in GeV^2
  !****************************************************************************
  real function pcm2Integrated( srts, iCh )

    use CallStack, only: Traceback
    use constants, only: mK

    real, intent(in) :: srts
    integer, intent(in) :: iCh

    integer :: iSrts
    real :: h

    pcm2Integrated = 0.0 ! default return value

    if (srts < arrSrts0(iCh)) return
    if (srts < 2*mK) return ! no tabulation below

    ! we use 'int', since we really want the lower boundary
    iSrts = int((srts-2*mK)/dSrts)
!    if (iSrts < 1) then
!       ! here we have to consider the threshold
!       h = (srts-2*mK)/dSrts-iSrts !(srts-2*mK)-iSrts*dSrts
!       write(*,*) srts,h
!
!
!       pcm2Integrated = h * arrPCM(iCh,iSrts+1)
!    else if (iSrts >= nSrts) then
    if (iSrts >= nSrts) then
       ! this is beyond the tabulation, do the calculation
       call traceback("(iSrts >= nSrts)")
    else
       ! here we are now in the tabulated region
       ! do a linear look-up
       h = (srts-2*mK)/dSrts-iSrts !(srts-2*mK)-iSrts*dSrts
       pcm2Integrated = (1-h) * arrPCM(iCh,iSrts) + h * arrPCM(iCh,iSrts+1)
    end if


  end function pcm2Integrated

  !****************************************************************************
  !****if* mesonMeson_Xsections/pcm2Pole
  ! NAME
  ! real function pcm2Pole( srts, iCh )
  ! PURPOSE
  ! calculate the phase space (squared) by calculating the pCM assuming
  ! pole masses of the particles in the given channel.
  !
  ! Threshold is treated correctly by using the pole masses.
  !
  ! INPUTS
  ! * real :: srts -- c.m. energy (GeV)
  ! * integer :: iCh -- number of channel to consider
  ! OUTPUT
  ! * phase space (squared) in GeV^2
  !****************************************************************************
  real function pcm2Pole( srts, iCh )

    use particleProperties, only: hadron
    use twoBodyTools, only: pCM_sqr

    real, intent(in) :: srts
    integer, intent(in) :: iCh

    pcm2Pole = max(0.,pCM_sqr(srts**2, &
            hadron(channelID(1,iCh))%mass**2, &
            hadron(channelID(2,iCh))%mass**2) )
  end function pcm2Pole

  !****************************************************************************
  !****if* mesonMeson_Xsections/pcm2
  ! NAME
  ! real function pcm2( srts, iCh, useWidth )
  ! PURPOSE
  ! Shorthand for calling pcm2Integrated or pcm2Pole
  !****************************************************************************
  real function pcm2( srts, iCh, useWidth )

    real, intent(in) :: srts
    integer, intent(in) :: iCh
    logical, intent(in) :: useWidth

    if (useWidth) then
       pcm2 = pcm2Integrated(srts, iCh )
    else
       pcm2 = pcm2Pole(srts, iCh )
    end if
  end function pcm2

  !****************************************************************************
  !****s* mesonMeson_Xsections/mesMes_Tabulate
  ! NAME
  ! subroutine mesMes_Tabulate
  ! PURPOSE
  ! Initialize the cross sections by tabulating them as function of sqrt(s)
  !****************************************************************************
  subroutine mesMes_Tabulate

    use constants, only: mK
    use mesonWidth, only: GetSpectralIntegral
    use output, only: Write_InitStatus
    use particleProperties, only: hadron
    use twoBodyTools, only: pCM

    integer :: iSrts, iCh
    real :: srts, srts0
    real :: valI, fakInt
    logical,dimension(1:2) :: stable

    call Write_InitStatus("mesonMeson_Xsections/mesMes_Tabulate",0)

    do iCh=1,21
       stable(1:2) = (hadron(channelID(1:2,iCh))%width < 1e-3)
       srts0 = hadron(channelID(1,iCh))%minmass + hadron(channelID(2,iCh))%minmass
       fakInt = 1.0/(GetSpectralIntegral(channelID(1,iCh))*GetSpectralIntegral(channelID(2,iCh)))

       arrSrts0(iCh) = srts0
       arrFakInt(iCh) = fakInt

       write(*,*) 'i=',iCh,channelID(1:2,iCh),stable(1:2),srts0,fakInt

       do iSrts=0,nSrts
          srts= iSrts*dSrts + 2*mK
          if (srts < srts0) cycle

          if (stable(1).and.stable(2)) then
             valI = pCM( srts, &
                  hadron(channelID(1,iCh))%mass, &
                  hadron(channelID(2,iCh))%mass )**2
          else if (stable(1)) then
             valI = calculate1(srts,channelID(1,iCh),channelID(2,iCh))
          else if (stable(2)) then
             valI = calculate1(srts,channelID(2,iCh),channelID(1,iCh))
          else
             valI = calculate2(srts,channelID(1,iCh),channelID(2,iCh))
          end if

          arrPCM(iCh,iSrts) = valI * fakInt

       end do
    end do

    call Write_InitStatus("mesonMeson_Xsections/mesMes_Tabulate",1)

  end subroutine mesMes_Tabulate


  !****************************************************************************
  !****if* mesonMeson_Xsections/calculate1
  ! NAME
  ! real function calculate1(srts, id1, id2)
  ! PURPOSE
  ! Calculates the integral of pCM**2 for id1 stable, id2 instable
  !****************************************************************************
  function calculate1(srts, id1, id2) result(integral)

    use constants, only: pi
    use mesonWidth, only: fullWidthMeson
    use particleProperties, only: hadron
    use twoBodyTools, only: pCM

    real,    intent(in)  :: srts
    integer, intent(in)  :: id1,id2
    real                 :: integral

    real :: mass1,mass2,mu1,mu2,pfinal
    integer :: nmu2,j
    real :: gamma1,gamma2,gamma2Tot
    real, parameter :: dy=pi/100.
    real :: dy2,minmu2,maxmu2,ymax2,ymin2,y2
    real :: spectral2,intfac2

    integral=0.

    mass1=hadron(id1)%mass
    mass2=hadron(id2)%mass
    gamma1=hadron(id1)%width
    gamma2=hadron(id2)%width

    !    write(*,*) "xx",mass1,mass2,gamma1,gamma2

    mu1 = mass1

    minmu2=hadron(id2)%minmass
    maxmu2=srts-mu1

    if (maxmu2 < minmu2) return

    ymax2=2.*atan((maxmu2-mass2) / gamma2*2.)
    ymin2=2.*atan((minmu2-mass2) / gamma2*2.)

    nmu2=max(int((ymax2-ymin2)/dy),1)
    dy2=(ymax2-ymin2)/float(nmu2)

    do j=1,nmu2 ! loop over second particle's mass
       y2=ymin2+(float(j)-0.5)*dy2
       mu2=.5*tan(y2/2.)*gamma2+mass2
       mu2=min(max(mu2,minmu2),maxmu2)

       pfinal = pCM(srts,mu1,mu2)**2
       gamma2Tot = fullWidthMeson(ID2, mu2)
       spectral2 = 2./pi * mu2**2 * gamma2tot / ((mu2**2-mass2**2)**2+mu2**2*gamma2tot**2)
       intfac2 = gamma2 / ((mu2-mass2)**2+gamma2**2/4.)

       !       write(*,*) spectral2,intfac2

       integral=integral+pfinal*spectral2/intfac2*dy2
    end do

  end function calculate1

  !****************************************************************************
  !****if* mesonMeson_Xsections/calculate2
  ! NAME
  ! real function calculate2(srts, id1, id2)
  ! PURPOSE
  ! Calculates the integral of pCM for id1 instable, id2 instable
  !****************************************************************************
  function calculate2(srts, id1, id2) result(integral)

    use constants, only: pi
    use mesonWidth, only: fullWidthMeson
    use particleProperties, only: hadron
    use twoBodyTools, only: pCM

    real,    intent(in)  :: srts
    integer, intent(in)  :: id1,id2
    real                 :: integral

    real :: mass1,mass2,mu1,mu2,pfinal
    integer :: nmu1,nmu2,i,j
    real :: gamma1,gamma2,gamma1Tot,gamma2Tot
    real, parameter :: dy=pi/100.
    real :: dy1,dy2,minmu1,maxmu1,minmu2,maxmu2,ymax1,ymin1,y1,y2
    real :: spectral1,spectral2,ymax2,ymin2,intfac1,intfac2

    integral=0.

    mass1=hadron(id1)%mass
    mass2=hadron(id2)%mass
    gamma1=hadron(id1)%width
    gamma2=hadron(id2)%width

    minmu1=hadron(id1)%minmass
    maxmu1=srts-hadron(id2)%minmass

    if (maxmu1 < minmu1) return

    ymax1=2.*atan((maxmu1-mass1) / gamma1*2.)
    ymin1=2.*atan((minmu1-mass1) / gamma1*2.)

    nmu1=max(int((ymax1-ymin1)/dy),1)
    dy1=(ymax1-ymin1)/float(nmu1)

    do i=1,nmu1 ! loop over first particle's mass
       y1=ymin1+(float(i)-0.5)*dy1
       mu1=.5*tan(y1/2.)*gamma1+mass1
       mu1=min(max(mu1,minmu1),maxmu1)

       minmu2=hadron(id2)%minmass
       maxmu2=srts-mu1

       ymax2=2.*atan((maxmu2-mass2) / gamma2*2.)
       ymin2=2.*atan((minmu2-mass2) / gamma2*2.)

       nmu2=max(int((ymax2-ymin2)/dy),1)
       dy2=(ymax2-ymin2)/float(nmu2)
       gamma1Tot = fullWidthMeson(ID1, mu1)

       spectral1 = 2./pi * mu1**2 * gamma1Tot / ((mu1**2-  mass1**2)**2+mu1**2*gamma1Tot**2)
       intfac1 = gamma1 / ((mu1-mass1)**2+gamma1**2/4.)

       do j=1,nmu2 ! loop over second particle's mass
          y2=ymin2+(float(j)-0.5)*dy2
          mu2=.5*tan(y2/2.)*gamma2+mass2
          mu2=min(max(mu2,minmu2),maxmu2)

          pfinal = pCM(srts,mu1,mu2)**2
          gamma2Tot = fullWidthMeson(ID2, mu2)

          spectral2 = 2./pi * mu2**2 * gamma2tot / ((mu2**2-mass2**2)**2+mu2**2*gamma2tot**2)
          intfac2 = gamma2 / ((mu2-mass2)**2+gamma2**2/4.)

          integral=integral+pfinal*spectral1*spectral2/(intfac1*intfac2)*dy1*dy2
       end do
    end do

  end function calculate2



end module mesonMeson_Xsections
