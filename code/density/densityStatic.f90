!******************************************************************************
!****m* /densityStatic
! NAME
! module densityStatic
!
! PURPOSE
! Collect routines for STATIC density calculations.
!******************************************************************************
module densityStatic

  use nucleusDefinition
  use dichteDefinition
  use Callstack, only: Traceback

  implicit none


  private

  public :: staticDensity,staticDensityInit,densityLuis
  public :: ReAdjust


contains

  !****************************************************************************
  !****s* densityStatic/staticDensityInit
  ! NAME
  ! subroutine staticDensityInit (nuc)
  ! PURPOSE
  ! decide, which density parametrisation is used. Then tabulate this and
  ! also set the extreme values for the MC decision.
  ! INPUTS
  ! * type(tNucleus) :: nuc    -- nucleus which is regarded
  ! OUTPUT
  ! * type(tNucleus) :: nuc
  !****************************************************************************
  subroutine staticDensityInit (nuc)

    type(tNucleus),pointer :: nuc

    if (nuc%mass <= 2) then
       if (nuc%densityswitch_static > 0) then
          write(*,*)
          write(*,'(A)') '  WARNING: static density not possible with Mass<=2!'
          write(*,'(A)') '  setting densityswitch_static to 0 for this run!!!'
          write(*,*)
          nuc%densityswitch_static = 0
       end if
       return
    end if


    select case (nuc%densitySwitch_static)
    case (0)
       ! Density=0.
       call TabulateZero(nuc)

    case (1) ! Static Woods-Saxon distribution
       call TabulateDensityWoodsSaxon(nuc)

    case (2) ! Static density prescription implemented by L.Alvarez-Russo
            ! corresponds to Oset's papers, e.g. NPA 554

       call TabulateDensityLuis(nuc)

    case (3) ! Static density prescription according to Horst Lenske
            ! Implements different radii for neutrons and protons

       call TabulateDensityLenske(nuc)

    case (4) ! Static Harmonic Oscilator Shell model

       call TabulateDensityHarmOsc(nuc)

    case (5) ! Fermi gas model with no surface term

       call TabulateSphere(nuc)

    case (6) !Static Density based on LDA, implemented by Birger Steinmueller

       call TabulateDensityBirger(nuc)

    case (7) !Static Density based on LDA + Welke potential

       call TabulateDensityBirgerWelke(nuc)

    case (8) !Static Density prescription according Relativistic Thomas-Fermi
            !Valid only in RMF-mode
       call TabulateDensityExRTF(nuc)

    case default

       write(*,*) 'Error in static density: DensitySwitch_static is not well defined', nuc%densitySwitch_static
       call Traceback('Severe Error : STOP')
    end select

    call SearchMaxVals(nuc)
    call NucleusAverageDensity(nuc)

  end subroutine staticDensityInit


  !****************************************************************************
  !****s* densityStatic/TabulateZero
  ! NAME
  ! subroutine TabulateZero(nuc)
  ! PURPOSE
  ! Set the density table to zero
  ! INPUTS
  ! * type(tNucleus) :: nuc
  !
  ! OUTPUT
  ! * type(tNucleus) :: nuc
  !
  !****************************************************************************
  subroutine TabulateZero(nuc)
    use output

    type(tNucleus),pointer :: nuc

    integer :: i

    call Write_InitStatus('density tabulation (density=0.)',0)
    write(*,'(A,F8.4)') '    dr   =',nuc%dx
    write(*,'(A,F8.4)') '    r_max=',nuc%dx*float(nuc%maxIndex)

    do i=0,nuc%MaxIndex
       nuc%densTab(i,1:2) = 0.
    end do
    call Write_InitStatus('density tabulation (density=0.)',1)

  end subroutine TabulateZero



  !****************************************************************************
  !****s* densityStatic/TabulateSphere
  ! NAME
  ! subroutine TabulateFermiGas(nuc)
  ! PURPOSE
  ! Tabulate a sphere with constant density.
  !
  ! Tabulates the static density to make it available faster for later use
  ! INPUTS
  ! * type(tNucleus) :: nuc
  !
  ! OUTPUT
  ! * type(tNucleus) :: nuc
  !
  !****************************************************************************
  subroutine TabulateSphere(nuc)
    use output

    type(tNucleus),pointer :: nuc

    real :: x,ratio!,h
    integer :: i
    real, parameter :: epsilon=0.001

    call Write_InitStatus('density tabulation (Sphere)',0)
    write(*,'(A,F8.4)') '    dr   =',nuc%dx
    write(*,'(A,F8.4)') '    r_max=',nuc%dx*float(nuc%maxIndex)

    if (nuc%radius.lt.epsilon) then
       write(*,*) 'SEVERE ERROR: This nucleus is not well defined.'
       write(*,*) 'Radius=',nuc%radius
       write(*,*) 'A=',nuc%mass,'  Z= ',nuc%charge
       write(*,*) 'Not possible to use that nucleus for sphere density.'
       call Traceback('Stop')
    end if

    ratio = float(nuc%charge)/float(nuc%mass)

    do i=0,nuc%MaxIndex
       x = i*nuc%dx
       if (x.le.nuc%radius) then
          nuc%densTab(i,1) = nuc%density * ratio
          nuc%densTab(i,2) = nuc%density * (1-ratio)
       else
          nuc%densTab(i,1:2) = 0.
       end if
    end do
    call Write_InitStatus('density tabulation (Sphere)',1)

  end subroutine TabulateSphere


  !****************************************************************************
  !****s* densityStatic/TabulateDensityWoodsSaxon
  ! NAME
  ! subroutine TabulateDensityWoodsSaxon (nuc)
  ! PURPOSE
  ! Tabulate the Woods-Saxon distribution. Tabulates the static density
  ! to make it available faster for later use.
  ! INPUTS
  ! * type(tNucleus) :: nuc
  ! OUTPUT
  ! * type(tNucleus) :: nuc
  !****************************************************************************
  subroutine TabulateDensityWoodsSaxon (nuc)
    use output, only: write_initstatus

    type(tNucleus),pointer :: nuc

    real :: x,h,ratio
    integer :: i
    real, parameter :: epsilon=0.001

    call Write_InitStatus('density tabulation (Woods-Saxon)',0)
    write(*,'(A,F8.4)') '    dr   =',nuc%dx
    write(*,'(A,F8.4)') '    r_max=',nuc%dx*float(nuc%maxIndex)

    if (nuc%surface<epsilon .or. nuc%radius<epsilon) then
       write(*,*) 'SEVERE ERROR: This nucleus is not well defined.'
       write(*,*) 'Radius=',nuc%radius,'   Surface=',nuc%surface
       write(*,*) 'A=',nuc%mass,'  Z= ',nuc%charge
       write(*,*) 'Not possible to use that nucleus for Woods-Saxon density.'
       call Traceback('Stop')
    end if

    ratio = float(nuc%charge)/float(nuc%mass)

    do i=0,nuc%MaxIndex
       x = i*nuc%dx
       h = nuc%density /(1.+exp((x-nuc%radius)/nuc%surface))
       nuc%densTab(i,1) = h * ratio
       nuc%densTab(i,2) = h * (1-ratio)
    end do
    call Write_InitStatus('density tabulation (Woods-Saxon)',1)

  end subroutine TabulateDensityWoodsSaxon


  !****************************************************************************
  !****s* densityStatic/TabulateDensityHarmOsc
  ! NAME
  ! subroutine TabulateDensityHarmOsc(nuc)
  ! PURPOSE
  ! Tabulate the density distribution according harmonic oscillator shell
  ! modell.
  !
  ! Tabulates the static density to make it available faster for later use
  ! INPUTS
  ! * type(tNucleus) :: nuc
  !
  ! OUTPUT
  ! * type(tNucleus) :: nuc
  !
  ! NOTES
  ! The parameter are taken from the FRITIOF package.
  !****************************************************************************
  subroutine TabulateDensityHarmOsc(nuc)
    use output
    use constants

    type(tNucleus),pointer :: nuc

    real :: x,h,ratio, h1,h2,h3
    integer :: i
    real, parameter :: epsilon=0.001
    logical :: okay

    integer, parameter:: PossibleA(2:8) = (/   4, -1,     9,   11,    12, -1,    16/)
    real,    parameter:: rCh(2:8) =       (/1.74, -1., 2.519, 2.37, 2.446, -1., 2.724/)
    real,    parameter:: d2(2:8)   = (rCh**2-0.81**2) / (2.5-4./PossibleA)


    call Write_InitStatus('density tabulation (Harm. Osc.)',0)
    write(*,'(A,F8.4)') '    dr   =',nuc%dx
    write(*,'(A,F8.4)') '    r_max=',nuc%dx*float(nuc%maxIndex)

    okay = .true.
    if ( (nuc%charge<2) .or. (nuc%charge>8) ) then
       okay = .false.
    else
       if (nuc%mass.ne.PossibleA(nuc%charge)) okay = .false.
    end if

    if (.not.okay) then
       write(*,*) 'SEVERE ERROR:'
       write(*,*) 'Radius=',nuc%radius,'   Surface=',nuc%surface
       write(*,*) 'A=',nuc%mass,'  Z= ',nuc%charge
       write(*,*) 'Not possible to use that nucleus for harm. osc. density.'
       call Traceback('Stop')
    end if

    write(*,*) 'Parameters: r_ch = ',rCh(nuc%charge),' d2 = ',d2(nuc%charge)

    ratio = float(nuc%charge)/float(nuc%mass)

    h1 = 4/(pi * d2(nuc%charge))**(3./2.)
    h2 = (nuc%mass-4)/6.

    do i=0,nuc%MaxIndex
       x = i*nuc%dx
       h3 = x**2/d2(nuc%charge)
       h = h1 * (1+h2*h3) * exp(-h3)

       nuc%densTab(i,1) = h * ratio
       nuc%densTab(i,2) = h * (1-ratio)
    end do
    call Write_InitStatus('density tabulation (Harm. Osc.)',1)

  end subroutine TabulateDensityHarmOsc


  !****************************************************************************
  !****s* densityStatic/TabulateDensityLenske
  ! NAME
  ! subroutine TabulateDensityLenske (nuc)
  ! PURPOSE
  ! Tabulate the density distribution according to Woods-Saxon distribution
  ! but with refined charge radii for proton and neutron according to H. Lenske.
  ! Tabulates the static density to make it available faster for later use.
  ! INPUTS
  ! * type(tNucleus) :: nuc
  ! OUTPUT
  ! * type(tNucleus) :: nuc
  ! NOTES
  ! Everything in fm.
  !****************************************************************************
  subroutine TabulateDensityLenske (nuc)
    use nucD, only: dfs
    use output, only: write_initstatus

    type(tNucleus),pointer :: nuc

    real, dimension(1:3) :: radius, surface, rhoMax
    real :: x
    integer :: i

    call Write_InitStatus('density tabulation (Lenske & Woods-Saxon)',0)
    write(*,'(A,F8.4)') '    dr   =',nuc%dx
    write(*,'(A,F8.4)') '    r_max=',nuc%dx*float(nuc%maxIndex)

    call DFS(nuc%mass,nuc%charge,radius,surface,rhoMax)

    do i=0,nuc%MaxIndex
       x = i*nuc%dx
       nuc%densTab(i,1) = rhoMax(1)/(1+exp((x-radius(1))/surface(1)))
       nuc%densTab(i,2) = rhoMax(2)/(1+exp((x-radius(2))/surface(2)))
    end do
    call Write_InitStatus('density tabulation (Lenske & Woods-Saxon)',1)

  end subroutine TabulateDensityLenske


  !****************************************************************************
  !****s* densityStatic/TabulateDensityLuis
  ! NAME
  ! subroutine TabulateDensityLuis(nuc)
  ! PURPOSE
  ! Tabulate the density distribution of matter (p and n)
  ! and the density of centers ( p and n number densities )
  ! following J. Nieves et al., Pionic atoms... NPA554
  !
  ! Tabulates the static density to make it available faster for later use
  ! INPUTS
  ! * type(tNucleus) :: nuc
  !
  ! OUTPUT
  ! * type(tNucleus) :: nuc
  !
  ! NOTES
  ! everything in fm
  !****************************************************************************
  subroutine TabulateDensityLuis(nuc)
    use output

    type(tNucleus),pointer :: nuc

    real :: x
    real :: rp,ap,rho0p,rn,an,rho0n
    integer :: i

    call Write_InitStatus('density tabulation (NPA554)',0)
    write(*,'(A,F8.4)') '    dr   =',nuc%dx
    write(*,'(A,F8.4)') '    r_max=',nuc%dx*float(nuc%maxIndex)

    do i=0,nuc%MaxIndex
       x = i*nuc%dx
       call densityLuis(x,nuc%charge,nuc%mass,nuc%densTab(i,1),nuc%densTab(i,2),rp,ap,rho0p,rn,an,rho0n)
    end do

    write(*,paragraph) ' Parameters: '
    write(*,'(A,3(1x,F8.4))') ' rp, ap, rho0p : ', rp,ap,rho0p
    write(*,'(A,3(1x,F8.4))') ' rn, an, rho0n : ', rn,an,rho0n
    call Write_InitStatus('density tabulation (NPA554)',1)

  end subroutine TabulateDensityLuis

  !****************************************************************************
  !****s* densityStatic/TabulateDensityExRTF
  ! NAME
  ! subroutine TabulateDensityExRTF(nuc)
  ! PURPOSE
  ! Tabulate the density distribution according to Relativistic
  ! Thomas-Fermi model code from Horst Lenske.
  !
  ! Tabulates the static density to make it available faster for later use
  ! INPUTS
  ! * type(tNucleus) :: nuc
  !
  ! OUTPUT
  ! * type(tNucleus) :: nuc
  !
  ! NOTES
  ! everything in fm
  !****************************************************************************
  subroutine TabulateDensityExRTF(nuc)
    use NucExRTF, only: NucExRTF_Main
    use output

    type(tNucleus),pointer :: nuc

    integer :: i
    real, dimension(0:2000, 1:2) :: RTF_dens

    call Write_InitStatus('density tabulation (ExRTF Lenske)',0)
    write(*,'(A,F8.4)') '    dr   =',nuc%dx
    write(*,'(A,F8.4)') '    r_max=',nuc%dx*float(nuc%maxIndex)

    call NucExRTF_Main(nuc%mass,nuc%charge,RTF_dens)

    do i=0,nuc%MaxIndex
       nuc%densTab(i,1) = RTF_dens(i,1)
       nuc%densTab(i,2) = RTF_dens(i,2)
    end do
    call Write_InitStatus('density tabulation (ExRTF Lenske)',1)

  end subroutine TabulateDensityExRTF


  !****************************************************************************
  !****s* densityStatic/SearchMaxVals
  ! NAME
  ! subroutine SearchMaxVals(nuc)
  !
  ! PURPOSE
  ! go through the tabulated distributions to search for the extrema
  !
  !****************************************************************************
  subroutine SearchMaxVals(nuc)
    use constants, only:pi

    type(tNucleus),pointer :: nuc

    integer :: i,j
    real :: x, maxV(2), maxX
    real :: Sum(0:3,3)

    maxV = 0.
    maxX = 0.
    Sum = 0.


    do i=0,nuc%MaxIndex
       x = i*nuc%dx
       if (nuc%densTab(i,1) > maxV(1)) maxV(1) = nuc%densTab(i,1)
       if (nuc%densTab(i,2) > maxV(2)) maxV(2) = nuc%densTab(i,2)
       if (nuc%densTab(i,1)+nuc%densTab(i,2) > 1e-6) maxX = x

       do j=0,2
          Sum(j,1) = Sum(j,1)+x**(j+2)*nuc%densTab(i,1)
          Sum(j,2) = Sum(j,2)+x**(j+2)*nuc%densTab(i,2)
          Sum(j,3) = Sum(j,3)+x**(j+2)*(nuc%densTab(i,1)+nuc%densTab(i,2))
       end do

       Sum(3,1) = Sum(3,1)+x**2*nuc%densTab(i,1)**2
       Sum(3,2) = Sum(3,2)+x**2*nuc%densTab(i,2)**2
       Sum(3,3) = Sum(3,3)+x**2*(nuc%densTab(i,1)+nuc%densTab(i,2))**2


    end do

    Sum = Sum*nuc%dx

    nuc%MaxDist = maxX
    nuc%MaxDens = maxV

    write(*,*) 'Extrema for MC:'
    write(*,*) '    MaxDist = ',nuc%MaxDist
    write(*,*) '    MaxDens = ',nuc%MaxDens

    write(*,*) 'Integrations:'
    write(*,'("   ",A8,4A13)')  ' ', 'rho d^3r', 'rho r d^3r', 'rho r^2 d^3r', 'rho^2 d^3r'
    write(*,'("   ",A8,4f13.3)') 'Proton:  ',Sum(0:3,1)*4.*pi
    write(*,'("   ",A8,4f13.3)') 'Neutron: ',Sum(0:3,2)*4.*pi
    write(*,'("   ",A8,4f13.3)') 'Baryon:  ',Sum(0:3,3)*4.*pi
    write(*,*)

    if (nuc%radius < 0.001) then
       nuc%radius = sqrt(Sum(2,3)/Sum(0,3))
       write(*,*) 'Attention! Setting radius to ',nuc%radius
    end if

  end subroutine SearchMaxVals


  !****************************************************************************
  !****f* densityStatic/staticDensity
  ! NAME
  ! type(dichte) function staticDensity(r,nucl)
  ! PURPOSE
  ! gives density in the restframe of the nucleus "nucl" at position "r"
  ! INPUTS
  ! * real, dimension(1:3),intent(in) :: r -- position where density should be calculated
  ! * type(tNucleus),pointer    :: nucl    -- nucleus which is regarded
  ! USAGE
  ! (dichte)=staticDensity(...)
  !****************************************************************************
  type(dichte) function staticDensity(r,nucl)

    real, dimension(1:3),intent(in) :: r
    type(tNucleus),pointer       :: nucl

    real :: sqrtR
    integer :: i

    if (nucl%DoInit) then
       call Traceback('nucleus not initialized! stop')
    end if


    select case (nucl%densitySwitch_static)
    case (0)
       staticdensity%proton  = 0.
       staticdensity%neutron = 0.
       staticdensity%baryon  = 0.

    case (1:8) ! Tabulated density distributions

       sqrtR=sqrt(r(1)**2+r(2)**2+r(3)**2)

       if (sqrtR>nucl%MaxIndex*nucl%dx) then
          staticdensity%proton(0)  = 0.
          staticdensity%neutron(0) = 0.
          staticdensity%baryon(0)  = 0.
       else
          i = min(nint(sqrtR/nucl%dx),nucl%maxIndex)
          staticdensity%proton(0)  = nucl%densTab(i,1)
          staticdensity%neutron(0) = nucl%densTab(i,2)
          staticdensity%baryon(0)  = nucl%densTab(i,1)+nucl%densTab(i,2)
       end if

       ! set flux to zero (lrf):
       staticdensity%proton(1:3)  = 0.
       staticdensity%neutron(1:3) = 0.
       staticdensity%baryon(1:3)  = 0.

    case default

       write(*,*) 'Error in static density: DensitySwitch_static is not well defined',&
            & nucl%densitySwitch_static
       call Traceback('Severe Error : STOP')
    end select

  end function staticDensity


  !****************************************************************************
  !****s* densityStatic/densityLuis
  ! NAME
  ! subroutine densityLuis(r,z,a,rhop_mat,rhon_mat,center_in)
  ! PURPOSE
  ! * This routine calculates the density of matter (p and n)
  !   and the density of centers ( p and n number densities )
  ! * following J. Nieves et al., Pionic atoms... NPA554
  ! * everything in fm
  ! * returns per default the density of matter, use "center_in" to switch to density of centers
  ! INPUTS
  ! * real, intent(in)::r  -- radius
  ! * integer, intent(in)::z  -- charge of the nucleus
  ! * integer, intent(in)::a  -- atomic number
  ! * logical, optional, intent(in) :: center_in -- if true then density of centers is given in the output
  ! RESULT
  ! * real, intent(out):: rhop,rhon -- Proton and neutron densities at r
  ! * real, intent(out):: rp,ap,rho0p,rn,an,rho0n -- parameters of the density distributions
  !****************************************************************************
  subroutine densityLuis(r,z,a,rhop,rhon,rp,ap,rho0p,rn,an,rho0n,center_in)
    use constants, only: pi
    use gauss_integration, only: sg20r, rg20r

    real, intent(in)::r
    integer, intent(in)::z
    integer, intent(in)::a
    real, intent(out):: rhop,rhon
    real, intent(out) :: rp,ap,rho0p,rn,an,rho0n
    logical, optional, intent(in) :: center_in
    logical :: center

    real::x,rpc,apc,rnc,anc,rmax,rin
    real::resu1,resu2,resu3,resu4
    real,dimension(:),allocatable::absi,orde1,orde2,orde3,orde4
    integer::nin,nins,i
    real, parameter::r2=0.69

    ! Choose between density of matter or density of centers
    if (present(center_in)) then
       center=center_in
    else
       center=.false.
    end if

    call denspar(z,a,rp,ap,rn,an)

    if (a.le.18) then
       !       For light nuclei, harmonic oscilator densities are used
       !       protons
       rpc=sqrt(rp**2-2./3.*r2)
       x=ap*rp**2/(1.+3./2.*ap)/rpc**2
       apc=2.*x/(2.-3.*x)
       !       neutrons
       rnc=sqrt(rn**2-2./3.*r2)
       x=an*rn**2/(1.+3./2.*an)/rnc**2
       anc=2.*x/(2.-3.*x)

    else
       !       Fermi liquid type
       !       protons
       rpc=rp+5.*r2*rp/(15.*rp**2+7.*pi**2*ap**2)
       apc=sqrt((rp**3+pi**2*ap**2*rpc-rpc**3)/pi**2/rpc)
       !       neutrons
       rnc=rn+5.*r2*rn/(15.*rn**2+7.*pi**2*an**2)
       anc=sqrt((rn**3+pi**2*an**2*rnc-rnc**3)/pi**2/rnc)

    end if

    !       normalization
    nin=20
    rmax=20.
    allocate (absi(20*nin))
    allocate (orde1(20*nin))
    allocate (orde2(20*nin))
    allocate (orde3(20*nin))
    allocate (orde4(20*nin))
    call sg20r(0.,rmax,nin,absi,nins)
    do i=1,nins
       rin=absi(i)
       orde1(i)=4.*pi*rin**2*antz(rin,rp,ap)
       orde2(i)=4.*pi*rin**2*antz(rin,rn,an)
       orde3(i)=4.*pi*rin**2*antz(rin,rpc,apc)
       orde4(i)=4.*pi*rin**2*antz(rin,rnc,anc)
    end do
    call rg20r(0.,rmax,nin,orde1,resu1)
    call rg20r(0.,rmax,nin,orde2,resu2)
    call rg20r(0.,rmax,nin,orde3,resu3)
    call rg20r(0.,rmax,nin,orde4,resu4)


    if (center) then
       ! use the distributions of nucleon centers:
       rho0p=z/resu3
       rho0n=(a-z)/resu4
       rp=rpc
       ap=apc
       rn=rnc
       an=anc
    else
       ! use matter distributions:
       rho0p=z/resu1
       rho0n=(a-z)/resu2
    end if

    rhop=rho0p*antz(r,rp,ap)
    rhon=rho0n*antz(r,rn,an)

    deallocate(absi,orde1,orde2,orde3,orde4)

  contains

    function antz(rin,rg,ag)
      implicit none
      real, intent(in)::rin,rg,ag
      real:: antz

      if (a.le.18) then
         !       harmonic oscilator
         antz=(1.+ag*(rin/rg)**2)*exp(-(rin/rg)**2)
      else
         !       Fermi liquid
         antz=1./(1.+exp((rin-rg)/ag))
      end if
    end function antz

    subroutine denspar(z,a,rp,ap,rn,an)
      implicit none
      integer, intent(in)::z,a
      real, intent(out)::rp,ap,rn,an

      select case (z)
      case (4)  ! Be (9)
         ! proton values from DeJager et al.,
         ! At. Data and Nucl. Data Tables 14, 479 (1974)
         ! neutron values from Koptev et al., Yad. Fiz. 31, 1501 (1980)
         rp=1.78
         ap=0.631
         rn=2.11
         an=1.000
      case (6)   ! C(12)
         rp=1.692
         ap=1.082
         rn=rp
         an=ap
         !Print *,"Im Kohlenstoff"
         !Stop
      case (8)
         if (a.eq.16) then
            !             O(16)
            rp=1.833
            ap=1.544
            rn=1.815
            an=1.529
         else
            !             O(18)
            rp=1.881
            ap=1.544
            rn=1.975
            an=2.048
         end if
      case (13) ! Al(27)
         rp=3.05
         ap=0.535
         rn=rp
         an=ap
      case (20)
         if (a.eq.40) then
            !             Ca(40)
            rp=3.51
            ap=0.563
            rn=3.43
            an=ap
         else
            !             Ca(44)
            rp=3.573
            ap=0.563
            rn=3.714
            an=ap
         end if
      case (26) ! Fe(56)
         rp=3.971
         ap=0.5935
         rn=4.05
         an=ap
      case (29) ! Cu(63)
         rp=4.214
         ap=0.586
         rn=4.31
         an=ap
      case (33) ! As(75)
         rp=4.492
         ap=0.58
         rn=4.64
         an=ap
      case (58) ! Ce(142)
         rp=5.76
         ap=0.535
         rn=5.98
         an=ap
      case (50) ! Sn-isotopes:
         ! Values taken from R. Schmidt et al,
         ! PRC 67, 044308 (2003)
         ! --- p and n center distribution parameters,
         ! should not be further corrected
         select case (A)
         case (112)
            rp=5.416
            ap=0.497
            rn=rp
            an=0.543
         case (116)
            rp=5.399
            ap=0.486
            rn=rp
            an=0.552
         case (120)
            rp=5.356
            ap=0.515
            rn=rp
            an=0.565
         case (124)
            rp=5.530
            ap=0.467
            rn=rp
            an=0.558
         case default
            call Traceback('There is no init for this Sn isotope')
         end select
      case (73) ! Ta(181)
         ! proton values from DeJager et al
         ! neutron values from Koptev et al
         rp=6.38
         ap=0.64
         rn=6.42
         an=0.64
      case (79) ! Au(197)
         rp=6.55
         ap=0.522
         rn=6.79
         an=ap
      case (82)  ! Pb(208)
         rp=6.624
         ap=0.549
         rn=6.89
         an=0.549
      case default
         write(*,*) "For this core the density distribution according to NPA554 is not yet implemented!!"
         write(*,*) z

         call Traceback('Stop')
      end select
    end subroutine denspar

  end subroutine densityLuis


  !****************************************************************************
  !****s* densityStatic/TabulateDensityBirger
  ! NAME
  ! subroutine TabulateDensityBirger(nuc)
  ! PURPOSE
  ! Tabulate the density distribution based on a local density approximation
  ! first described by Brueckner et al.
  !
  ! INPUTS
  ! * type(tNucleus) :: nuc
  !
  ! OUTPUT
  ! * type(tNucleus) :: nuc
  !
  !****************************************************************************
  subroutine TabulateDensityBirger(nuc)
    use nucDLDA, only: DFLDA
    use output

    type(tNucleus),pointer :: nuc

    write(*,paragraph) 'Initializing density tabulation, LDA by Birger'
    write(*,'(A,F8.4)') '    dr   =',nuc%dx
    write(*,'(A,F8.4)') '    r_max=',nuc%dx*float(nuc%maxIndex)

    call DFLDA(nuc)

    write(*,*) 'Chemical Potential', nuc%chemPot

    write(*,paragraph) 'Finished initializing density tabulation'
    write(*,*)


  end subroutine TabulateDensityBirger

  !****************************************************************************
  !****s* densityStatic/TabulateDensityBirgerWelke
  ! NAME
  ! subroutine TabulateDensityBirgerWelke(nuc)
  ! PURPOSE
  ! Tabulate the density distribution based on a local density approximation
  ! first described by Brueckner et al. and a momentum-dependent potential
  !
  ! INPUTS
  ! * type(tNucleus) :: nuc
  !
  ! OUTPUT
  ! * type(tNucleus) :: nuc
  !
  !****************************************************************************
  subroutine TabulateDensityBirgerWelke(nuc)
    use nucDLDA, only: DFLDAWelke
    use output

    type(tNucleus),pointer :: nuc

    write(*,paragraph) 'Initializing density tabulation, LDA and Welke potential by Birger'
    write(*,'(A,F8.4)') '    dr   =',nuc%dx
    write(*,'(A,F8.4)') '    r_max=',nuc%dx*float(nuc%maxIndex)

    call DFLDAWelke(nuc)

    write(*,*) 'Chemical Potential', nuc%chemPot

    write(*,paragraph) 'Finished initializing density tabulation'
    write(*,*)


  end subroutine TabulateDensityBirgerWelke

  !****************************************************************************
  !****s* densityStatic/ReAdjust
  ! NAME
  ! subroutine ReAdjust(nuc, potP, potN, potC)
  !
  ! PURPOSE
  ! This routine recalculates the density distributions for protons and
  ! neutrons by considering the given potentials as static and fulfill the
  ! condition
  !   sqrt(p_F^2+m_N^2) + U - m_N == E_sep ~ -8MeV
  ! With the Local-Thomas-Fermi, we connect the resulting fermi momentum to
  ! a density,
  !   rho = p_F^3/(3pi^2)
  ! Since the potentials are given as function of r, we calculate rho(r).
  !
  ! Thus, given proton and nucleon baryon potential (for fixed momentum) and
  ! the coulomb potential, the parametrization of the nuclear density is
  ! readjusted.
  !
  ! This routine is called by
  ! baryonPotentialModule/HandPotentialToDensityStatic
  !
  !
  ! INPUTS
  ! * type(tNucleus),pointer :: nuc -- the nucleus to consider
  ! * real, dimension(0:) :: potP, potN -- the proton,neutron potentials
  !   with p=pF. The dimension has to be identical to nuc%densTab(0: ,1:2).
  ! * real, dimension(0:) :: potC -- The Coulomb potential (>0, in GeV)
  !
  ! OUTPUT
  ! * nuc%densTab(0: ,1:2) is changed
  !****************************************************************************
  subroutine ReAdjust(nuc, potP, potN, potC)
    use constants
    use particleDefinition
    use output, only: IntToChar

    logical, parameter :: verbose = .false.

    type(tNucleus),pointer :: nuc
    real, dimension(0:), intent(in) :: potP, potN, potC

    integer :: i
    real :: x, rho, pF

    real :: SumP,SumN, facP,facN, pFN, pFP
    integer, save :: nCall = 0

    nCall = nCall+1

    if (verbose) then
       write(*,*) 'in Readjust...'

       open(113,file='ReAdjust.orig.'//IntToChar(nCall)//'.dat', status='unknown')
       rewind(113)
       write(113,'("#",A12,10A13)') 'x','rhoN','rhoP','potN','potP','potC','pF_N','pF_P','pF_B'
       do i=0,nuc%MaxIndex
          x = i*nuc%dx
          pFN = (3*pi**2*(nuc%densTab(i,2)))**(1./3.)*hbarc
          pFP = (3*pi**2*(nuc%densTab(i,1)))**(1./3.)*hbarc
          pF  = (3*pi**2*(nuc%densTab(i,1)+nuc%densTab(i,2))/2)**(1./3.)*hbarc
          write(113,'(10f13.5)') x,nuc%densTab(i,2),nuc%densTab(i,1),&
               PotN(i),PotP(i),PotC(i), &
               pFN, pFP, pF
       end do
       close(113)
    end if

    SumP = 0
    SumN = 0
    do i=0,nuc%MaxIndex
       x = i*nuc%dx

!       pF = 2*mN*(-potP(i)+nuc%ConstBinding)-potP(i)**2+nuc%ConstBinding**2
       pF = (-(potP(i)+potC(i))+nuc%ConstBinding+mN)**2 - mN**2
       if (pF.lt.0.0) then
          nuc%densTab(i,1) = 0.0
       else
          pF = sqrt(pF)
          rho = (pF/hbarc)**3/(3*pi**2)
          nuc%densTab(i,1) = rho
          SumP = SumP + x**2*rho

!          write(*,*) x,pF,rho

       end if


!       pF = 2*mN*(-potN(i)+nuc%ConstBinding)-potN(i)**2+nuc%ConstBinding**2
       pF = (-potN(i)+nuc%ConstBinding+mN)**2 - mN**2
       if (pF.lt.0.0) then
          nuc%densTab(i,2) = 0.0
       else
          pF = sqrt(pF)
          rho = (pF/hbarc)**3/(3*pi**2)
          nuc%densTab(i,2) = rho
          SumN = SumN + x**2*rho
       end if
    end do

    if (verbose) then
       write(*,*) 'P: ',SumP*4*pi*nuc%dx,nuc%charge
       write(*,*) 'N: ',SumN*4*pi*nuc%dx,nuc%mass-nuc%charge
    end if

    if ((SumP.eq.0).or.(SumN.eq.0)) then
       call Traceback('Failure!')
    end if

    ! Ensure normalization to mass number:

    facN=(nuc%mass-nuc%charge)/(SumN*4*pi*nuc%dx)
    facP=(nuc%charge)/(SumP*4*pi*nuc%dx)

    ! method 1: just set the normalization

    if (verbose) then
       write(*,*) 'Scaling rhoN by ',facN
       write(*,*) 'Scaling rhoP by ',facP
    end if

    nuc%densTab(:,2)=nuc%densTab(:,2)*facN
    nuc%densTab(:,1)=nuc%densTab(:,1)*facP

    nuc%facN = facN
    nuc%facP = facP

    ! method 2: rescale the radius

!    write(*,*) 'Scaling x by ',facN**(1./3.)
!    nuc%dx = nuc%dx*facN**(1./3.)


    call SearchMaxVals(nuc)
    call NucleusAverageDensity(nuc)

    if (verbose) then
       open(113,file='ReAdjust.new.'//IntToChar(nCall)//'.dat', status='unknown')
       rewind(113)
       write(113,'("#",A12,10A13)') 'x','rhoN','rhoP','potN','potP','potC','pF_N','pF_P','pF_B'
       do i=0,nuc%MaxIndex
          x = i*nuc%dx
          pFN = (3*pi**2*(nuc%densTab(i,2)))**(1./3.)*hbarc
          pFP = (3*pi**2*(nuc%densTab(i,1)))**(1./3.)*hbarc
          pF  = (3*pi**2*(nuc%densTab(i,1)+nuc%densTab(i,2))/2)**(1./3.)*hbarc
          write(113,'(10f13.5)') x,nuc%densTab(i,2),nuc%densTab(i,1),&
               PotN(i),PotP(i),PotC(i), &
               pFN, pFP, pF
       end do
       close(113)
    end if


  end subroutine ReAdjust



end module densityStatic
