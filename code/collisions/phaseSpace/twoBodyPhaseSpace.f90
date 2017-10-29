!******************************************************************************
!****m*  /twoBodyPhaseSpace
! NAME
! module twoBodyPhaseSpace
! NOTES
! Includes all the routines which are necessary for the two-body phase-space.
!******************************************************************************
module twoBodyPhaseSpace

  implicit none
  private

  public :: Integrate_2bodyPS_resonance
  public :: nnRR
  public :: setMaxSqrts

  logical, parameter :: DoWrite = .false. ! switch on/off printig of data files

  ! the following is not a parameter, but will be overwritten by 'setMaxSqrts',
  ! according to the high energy thresholds in master_2body
  integer, save :: srts_maxIndex = 220  ! this

  real, parameter :: delta_Srts = 0.01
  real, parameter :: srtsMin = 1.8

contains

  !****************************************************************************
  !****s*  twoBodyPhaseSpace/setMaxSqrts
  ! NAME
  ! subroutine setMaxSqrts(srtsMax)
  ! PURPOSE
  ! set the upper bound of the tabulation
  !****************************************************************************
  subroutine setMaxSqrts(srtsMax)
    real, intent(in) :: srtsMax
    srts_maxIndex = ceiling((srtsMax - srtsMin) / delta_Srts)
  end subroutine setMaxSqrts


  !****************************************************************************
  !****f*  twoBodyPhaseSpace/Integrate_2bodyPS_resonance
  ! NAME
  ! function Integrate_2bodyPS_resonance(resID, srts, massStable, scalarPotential) result (ps)
  ! PURPOSE
  ! Returns Integral over the CM-momentum with a resonance among the two
  ! particles. Therefore one has to Integrate over the
  ! mass of the resonance as well:
  !    ps=Integral { p_final(m_R) * Spectral function(m_R} d(m_R)
  ! with
  !    m_R = massResonance
  !
  ! It uses an internal routine "Calculate" which is either
  ! called directly or we store its values to a field.
  !
  ! This routine is not giving the real Phase space !!!
  ! INPUTS
  ! * integer, intent(in)  :: resID -- Id of the resonance
  ! * real, intent(in)     :: srts -- sqrt(s) in the problem
  ! * real, intent(in)     :: massStable -- mass of the stable particle
  ! * real, intent(in)     :: scalarPotential -- scalarPotential of the resonance
  ! OUTPUT
  ! * real,dimension(1:5)  :: ps --  Integral as given above
  !
  ! Here the different components are:
  ! * ps(1): With full width evaluated in the nominator of spectral function
  ! * ps(2): With partial width (pion N)  evaluated in the nominator of spectral function
  ! * ps(3): With partial width (eta N)   evaluated in the nominator of spectral function
  ! * ps(4): With partial width (rho N)   evaluated in the nominator of spectral function
  ! * ps(5): With partial width (omega N) evaluated in the nominator of spectral function
  ! NOTES
  ! Formerly known as "massInt2"
  !****************************************************************************
  function Integrate_2bodyPS_resonance(resID, srts, massStable, scalarPotential) result (ps)

    use particleProperties, only: hadron
    use idTable
    use constants, only: mN

    integer,intent(in) :: resID
    real, intent(in) :: srts, massStable, scalarPotential
    real, dimension(1:5) ::  ps

    logical, save :: initFlag=.true.
    ! To tabulate the output
    integer, parameter :: resID_min = nucleon
    integer, parameter :: resID_max = F37_1950
    real, save, allocatable, dimension(:,:,:,:) :: ps_Field
    integer :: down
    real :: srts_down, weight

    if ((abs(massStable-mN) < 0.001 .or. &
         abs(massStable-hadron(eta)%mass) < 0.001) .and.(abs(scalarPotential) < 0.001)) then

       if (initFlag) call initialize

       down = floor((srts-srtsMin)/delta_srts)

       if (down<0 .or. down>srts_maxIndex-1 .or. resID<resID_min .or. resID>resID_max) then
          ps = Calculate(resID, srts, massStable, scalarPotential)
       else
          ! use tabulated fields for eta and nucleon as stable particles (with linear interpolation)
          srts_down = srtsMin+float(down)*delta_srts
          weight = (srts-srts_down)/delta_srts ! weight for the interpolation
          if (abs(massStable-mN) < 0.001) then
             ps(1:5) = ps_Field(1, down, resID, 1:5) * (1.-weight) + ps_Field(1, down+1, resID, 1:5) * weight
          else
             ps(1:5) = ps_Field(2, down, resID, 1:5) * (1.-weight) + ps_Field(2, down+1, resID, 1:5) * weight
          end if
       end if

    else

       ps = Calculate(resID, srts, massStable, scalarPotential)

    end if


  contains

    !**************************************************************************
    ! Save the output for massStable=nucleon mass and scalarPotential=0.
    subroutine initialize

      use output, only: DoPR,Write_InitStatus,intToChar

      integer :: iSrts,iRes,j
      real :: srtsDummy
      real, dimension(1:5) :: psDummy
      real, dimension(1:2) :: stableDummy
      logical :: DoW


      DoW = DoWrite .and. DoPR(1)

      call Write_InitStatus("Integrate_2bodyPS_resonance",0)

      allocate(ps_Field(1:2,0:srts_maxIndex,resID_min:resID_max,1:5))

      stableDummy = (/ mN, hadron(eta)%mass /)


!$OMP PARALLEL DO COLLAPSE(2) &
!$OMP DEFAULT(NONE) &
!$OMP SHARED(ps_Field,srts_maxIndex,stableDummy) &
!$OMP PRIVATE(j,iRes,iSrts,srtsDummy)
      do j=1,2
         do iRes = resID_min, resID_max

            do iSrts=0, srts_maxIndex
               srtsDummy=float(iSrts)*delta_Srts+srtsMin
               ps_Field(j,iSrts,iRes,1:5) = &
                    Calculate(iRes, srtsDummy, stableDummy(j), 0.)

            end do
         end do
      end do
!$OMP END PARALLEL DO

      if (DoW) then
         write(*,*) 'Writing Output to TwobodyPS_R_*.dat. where * is the ID of the resonance'

         do iRes = resID_min, resID_max
            open(100,file='TwobodyPS_R_'//intToChar(iRes)//'.dat')
            write(100,*) '#sqrt(s),integral,...'

            do iSrts=0, srts_maxIndex
               srtsDummy=float(iSrts)*delta_Srts+srtsMin
               write(100,'(8g12.5)') srtsDummy, &
                    ps_Field(:,iSrts,iRes,1), &
                    stableDummy(:)
            end do
            close(100)

         end do
      end if

      call Write_InitStatus("Integrate_2bodyPS_resonance",1)

      initFlag = .false.

    end subroutine initialize

    !**************************************************************************
    function Calculate(resID, srts, massStable, scalarPotential) result(ps)

      use IdTable, only: isMeson
      use particleProperties, only: hadron
      use mesonWidth, only: fullWidthMeson
      use baryonWidth, only: fullWidthBaryon, partialWidthBaryon
      use constants, only: pi

      integer, intent(in) :: resID
      real, intent(in) :: srts, massStable, scalarPotential
      real, dimension(1:5) ::  ps

      real :: minmass,maxmass,mass, gamtot,spectral(1:5),pfinal
      real :: intfac,y,ymax,ymin,mres0,gamres0,dya
      integer :: i,nm
      real, parameter :: dy=2*pi/1000.

      logical :: mesonFlag

      ! Standard output:
      ps=0.

      ! Check Input, check wether resonance is meson or baryon
      mres0   = hadron(resID)%mass     ! pole mass
      gamres0 = hadron(resID)%width    ! pole width
      minmass = hadron(resID)%minmass
      mesonFlag = isMeson(resID)

      maxmass=srts-massStable-scalarPotential

      ! integrate over mass distribution

      if (maxmass <= minmass) return

      if (gamres0 < 1e-05) then

         if (srts > mres0+massStable) then
            ! Then ps=pCM of stable and unstable particle
            ps(1)=sqrt((srts**2-(mres0+scalarPotential+massStable)**2)*(srts**2-(mres0+scalarPotential-massStable)**2)/4./srts**2)
         end if

      else

         ymax = 2. * atan((maxmass-mres0)/gamres0*2.)
         ymin = 2. * atan((minmass-mres0)/gamres0*2.)
         nm   = max(int((ymax-ymin)/dy), 1)
         dya  = (ymax-ymin) / float(nm)

         do i=1, nm
            y=ymin+(float(i)-0.5)*dya

            mass=.5*tan(y/2.)*gamres0+mres0
            mass=min(max(mass,minmass),maxmass)

            pfinal = sqrt((srts**2-(mass+scalarPotential+massStable)**2)*(srts**2-(mass+scalarPotential-massStable)**2)/4./srts**2)

            if (mesonFlag) then
               gamtot = fullWidthMeson(resID, mass)
            else
               gamtot = fullWidthBaryon(resID, mass)
            end if

            spectral(1) = 2./pi * mass**2 * gamtot / ((mass**2-mres0**2)**2+gamtot**2*mass**2)
            if (mesonFlag .or. gamtot<1E-6) then
              spectral(2:5) = 0.
            else
              spectral(2) = spectral(1) * partialWidthBaryon(resID, mass, .false., pion,       nucleon) / gamtot
              spectral(3) = spectral(1) * partialWidthBaryon(resID, mass, .false., eta,        nucleon) / gamtot
              spectral(4) = spectral(1) * partialWidthBaryon(resID, mass, .false., rho,        nucleon) / gamtot
              spectral(5) = spectral(1) * partialWidthBaryon(resID, mass, .false., omegaMeson, nucleon) / gamtot
            end if

            intfac=gamres0/((mass-mres0)**2+gamres0**2/4.)
            ps(:) = ps(:) + pfinal*spectral(:)*dya/intfac
         end do

      end if

    end function Calculate


  end function Integrate_2bodyPS_resonance



  !****************************************************************************
  !****s*  twoBodyPhaseSpace/nnRR
  ! NAME
  ! real function nnRR(srts, ID)
  ! PURPOSE
  ! Procedure for calculation of integrals of Resonance Resonance'
  ! CM-momentum
  !
  ! Evaluates Integral over two body phase space in vacuum for two baryon
  ! resonances. Therefore one has to Integrate over the mass of the two
  ! resonances as well.
  !    ps=Integral d(mass_A )  d(mass_B ) p_AB * Spectralfunction_A Spectralfunction_B
  ! INPUTS
  ! * integer :: ID(1:2)  -- Ids of the resonances
  ! * real    :: srts     -- sqrt(s) of the reaction
  ! OUTPUT
  ! * Integral as given above
  ! NOTES
  ! To increase speed, we tabulate the possible output for Delta+R.
  !****************************************************************************
  real function nnRR(srts, ID)

    use IdTable, only: Delta, F37_1950
    use constants, only: mN, mPi

    real,    intent(in)  :: srts
    integer, intent(in)  :: ID(1:2)

    ! Tabulation parameters
    real, parameter :: srtsMin = 2 * (mN + mPi)
    real, dimension(:,:), allocatable, save :: integral_Field
    integer :: down
    real :: srts_down, weight
    logical, save :: initFlag=.true.
    logical,parameter :: speedFlag=.true. ! Switch tabulating on/off

    if (id(1)==Delta .and. speedFlag) then

       if (initFlag) call initialize

       down = floor((srts-srtsMin)/delta_srts)

       if (down<0 .or. down>srts_maxIndex-1) then
          !          write(*,*) 'WARNING : Slow speed in nnRR'
          nnRR = calculate(srts,id(1),id(2))
       else
          ! use tabulated fields (with linear interpolation)
          srts_down = srtsMin + float(down) * delta_srts
          weight = (srts-srts_down)/delta_srts
          nnRR = integral_Field(id(2),down  ) * (1.-weight) &
               + integral_Field(id(2),down+1) * weight
       end if
    else
       nnRR = calculate(srts,id(1),id(2))
    end if

  contains

    !**************************************************************************
    subroutine initialize

      use output, only: DoPR, Write_InitStatus

      integer :: i,id
      real :: srtsDummy
      logical :: DoW

      allocate(integral_Field(Delta:F37_1950,0:srts_maxIndex))

      DoW = DoWrite .and. DoPR(1)

      call Write_InitStatus("NNRR",0)

!$OMP PARALLEL DO &
!$OMP DEFAULT(NONE) &
!$OMP SHARED(integral_Field,srts_maxIndex) &
!$OMP PRIVATE(id,i,srtsdummy)
      do id = lbound(integral_Field,dim=1), ubound(integral_Field,dim=1)
         do i = 0, srts_maxIndex
            srtsDummy = float(i)*delta_Srts + srtsMin
            integral_Field(id,i) = calculate(srtsDummy, Delta, id)
         end do
      end do
!$OMP END PARALLEL DO

      if (DoW) then
         Open(100,file='NNRR.dat')
         write(100,*) '#sqrt(s),integral, maximal contribution,i'
         do id = lBound(integral_Field,dim=1), uBound(integral_Field,dim=1)
            do i = 0, srts_maxIndex
               srtsDummy = float(i)*delta_Srts + srtsMin
               write(100,'(F7.3,2I4,F7.3)') srtsDummy,id,i,integral_Field(id,i)
            end do
         end do
         close(100)
      end if
      call Write_InitStatus("NNRR",1)

      initFlag = .false.

    end subroutine initialize

    !**************************************************************************
    function calculate(srts, id1, id2) result(integral)

      use IdTable, only: isBaryon
      use particleProperties, only: hadron
      use baryonWidth, only: fullWidthBaryon
      use constants, only: pi

      real,    intent(in)  :: srts
      integer, intent(in)  :: id1,id2
      real                 :: integral

      real mass1,mass2,mu1,mu2,pfinal2,pfinal
      integer nmu1,nmu2,i,j
      real gamma1,gamma2,gamma1Tot,gamma2Tot
      real , parameter :: dy=pi/100.
      real dy1,dy2,minmu1,maxmu1,minmu2,maxmu2,ymax1,ymin1,y1,y2
      real spectral1,spectral2,ymax2,ymin2,intfac1,intfac2

      integral=0.

      if (isBaryon(id1) .and. isBaryon(id2)) then
         mass1=hadron(id1)%mass
         mass2=hadron(id2)%mass
         gamma1=hadron(id1)%width
         gamma2=hadron(id2)%width
      else
         write(*,*) 'Problem in nnRR : Particles are no baryons :' , id1, id2
         write(*,*) 'Routine is only meant to work for baryons'
         write(*,*) 'Critical error! stop!'
         stop
      end if

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

         do j=1,nmu2 ! loop over second particle's mass
            y2=ymin2+(float(j)-0.5)*dy2
            mu2=.5*tan(y2/2.)*gamma2+mass2
            mu2=min(max(mu2,minmu2),maxmu2)

            pfinal2=(srts**2-(mu1+mu2)**2)*(srts**2-(mu1-mu2)**2)/ (4.*srts**2)
            if (pfinal2.lt.0) then
               write(*,*) 'pfinal2 lt 0',pfinal2
               pfinal=0.
            else
               pfinal=sqrt(pfinal2)
            end if

            ! Evaluate widht of both particles
            gamma1Tot = fullWidthBaryon(ID1, mu1)
            gamma2Tot = fullWidthBaryon(ID2, mu2)

            spectral1 = 2./pi * mu1**2 * gamma1Tot / ((mu1**2-  mass1**2)**2+mu1**2*gamma1Tot**2)
            intfac1 = gamma1 / ((mu1-mass1)**2+gamma1**2/4.)

            spectral2 = 2./pi * mu2**2 * gamma2tot / ((mu2**2-mass2**2)**2+mu2**2*gamma2tot**2)
            intfac2 = gamma2 / ((mu2-mass2)**2+gamma2**2/4.)

            integral=integral+pfinal*spectral1*spectral2/intfac1/intfac2*dy1*dy2
         end do
      end do

    end function calculate


  end function nnRR



end module twoBodyPhaseSpace
