program OpticalGlauber

  use inputGeneral
  use version, only: PrintVersion
  use particleProperties, only: initParticleProperties
  use nucleusDefinition
  use nucleus, only : getProjectile,getTarget
  use dichteDefinition
  use densityStatic
  use output, only: Write_ReadingInput
  use CallStack, only: traceBack

  implicit none


  type(tNucleus),pointer :: NucA,NucB
  integer :: ios

  !***************************************************************************
  !****g* OpticalGlauber/impact
  ! SOURCE
  !
  real, save :: impact = 0.
  ! PURPOSE
  ! Impact parameter b [fm]. There are three options:
  ! * b>=0: The impact parameter is fixed to the given value.
  ! * -100<b<0: The impact parameter will be chosen randomly in each run between 0 and abs(b).
  ! * b<=-100: "Minimum bias". The impact parameter will be chosen randomly in each run (maximum = sum of radii plus twice the sum of surfaces).
  !***************************************************************************


 !***************************************************************************
  !****g* OpticalGlauber/impact_profile
  ! SOURCE
  !
  integer, save :: impact_profile = 0
  ! PURPOSE
  ! This switch provides impact-parameter distributions for trigger-biased
  ! setups. Only used for impact_parameter < 0.
  ! Possible values:
  ! * 0: minimum bias (default)
  ! * 1: HADES C+C    at 1.00 AGeV
  ! * 2: HADES C+C    at 2.00 AGeV
  ! * 3: HADES Ar+KCl at 1.76 AGeV
  ! * 4: HADES Au+Au  at 1.23 AGeV (all)
  ! * 5: HADES Au+Au  at 1.23 AGeV ( 0-10% central)
  ! * 6: HADES Au+Au  at 1.23 AGeV (10-20% central)
  ! * 7: HADES Au+Au  at 1.23 AGeV (20-30% central)
  ! * 8: HADES Au+Au  at 1.23 AGeV (30-40% central)
  !***************************************************************************

  real, save :: sigma = 30.0

  NAMELIST /Glauber/ impact, impact_profile, sigma

  real, parameter :: b0(6:8) = (/ 5.67, 7.27, 8.59 /)
  real, parameter :: db(6:8) = (/ 0.99, 0.85, 0.79 /)

  real, dimension(:,:), allocatable :: densTransA, densTransB
  real :: nColl, nPart
  real :: b = 0, bmax = 0, weight, deltab
  real, dimension(0:2) :: s
  integer :: ib
  integer, parameter :: nb = 100

  call PrintVersion

  call readInputGeneral
  call initParticleProperties

  NucA => getProjectile()
  NucB => getTarget()

  call Write_ReadingInput('Glauber',0)
  rewind(5)
  read(5,nml=Glauber,iostat=ios)
  call Write_ReadingInput('Glauber',0,ios)

  write(*,*) 'impact = ',impact,' fm'
  write(*,*) 'profile: ',impact_profile
  write(*,*) 'sigma  = ',sigma,' mb'

  call Write_ReadingInput('Glauber',1)

  allocate(densTransA(0:NucA%maxIndex,2))
  allocate(densTransB(0:NucB%maxIndex,2))

  call calcDensTrans(NucA, densTransA)
  call calcDensTrans(NucB, densTransB)

  call integrateDensTrans(NucA%dx,densTransA)
  call integrateDensTrans(NucB%dx,densTransB)

  if (impact >= 0.) then ! ===== impact parameter is fixed
     b = impact
     call doCalc(b, nColl, nPart)
  else ! ===== impact parameter integration
     if (impact > -100.) then
        bmax = abs(impact)
      else
        bmax = NucA%radius + NucB%radius + 2*(NucA%surface+NucB%surface)
      end if

      deltab=bmax/nb
      s = 0
      do ib=0,nb
         b=ib*deltab
         call doCalc(b, nColl, nPart)

         select case (impact_profile)
         case (0)
            ! default: minimum bias
            weight = 1.0
         case (1,2,3,4,5)
            ! use HADES trigger-biased centrality distributions with Woods-Saxon shape
            weight = impact_HADES(b,impact_profile)
         case (6,7,8)
!        ! use HADES trigger-biased centrality distributions with Gaussian shape
            !        b = rnGauss(db(impact_profile), b0(impact_profile))
            call TRACEBACK('not yet implemented')
         case default
            write(*,*) 'bad value for impact_profile: ', impact_profile
            call TRACEBACK('stop')
         end select

         s = s + weight*b*(/1.,nColl,nPart/)
!         write(*,*) b,weight
      end do

      nColl = s(1)/s(0)
      nPart = s(2)/s(0)

  end if

  write(*,*) 'nColl, nPart = ', nColl, nPart

contains

  subroutine calcDensTrans(Nuc, D)

    type(tNucleus),intent(in) :: Nuc
    real, dimension(0:,:), intent(inout) :: D

    integer :: iX,iZ,i
    real :: rX, rZ, r
    real, dimension(2) :: S, HH

    HH = (/real(Nuc%charge), real(Nuc%mass-Nuc%charge)/)

!    write(*,*) 'tabulating transversal density...'

    if (nuc%densityswitch_static == 0) then

       write(*,*) '(nuc%densityswitch_static == 0)'
       write(*,*) Nuc%dx,HH/(Nuc%dx**2)

       ! set everything to zero and insert some spike
       D = 0.0

       D(0,:) = HH/(Nuc%dx**2)

       return

    end if

!!$    S = 0.0
!!$    do iX=0,Nuc%maxIndex
!!$       rX = iX * Nuc%dx
!!$       S = S + rX**2 * nuc%densTab(iX,:)
!!$    end do
!!$    write(*,*) 'S = ', S*Nuc%dx*4*3.14159     ! --> Z, (A-Z)
!!$    write(*,*) '  = ', S*Nuc%dx*4*3.14159/HH  ! --> 1.0, 1.0


    do iX=0,Nuc%maxIndex
       rX = iX * Nuc%dx
       S = 0.0
       do iZ=0,Nuc%maxIndex
          rZ = iZ * Nuc%dx
          r = sqrt(rX**2 + rZ**2)
          i = nint(r/nuc%dx)
          if (i < nuc%maxIndex) then
             S = S + nuc%densTab(i,:)
          end if
       end do
       D(iX,:) = S * 2 * Nuc%dx
    end do

!!$    write(*,*) 'D(0): ',D(0,:)

!    write(*,*) 'tabulating transversal density... done.'
!    write(*,*)

  end subroutine calcDensTrans

  real function getDensTrans(x,y, dr, D)
    real, intent(in) :: x,y,dr
    real, dimension(0:,:), intent(in) :: D

    integer :: ir

    getDensTrans = 0.0
    ir = nint(sqrt(x**2+y**2)/dr)
    if (ir <= ubound(D,dim=1)) getDensTrans = sum(D(ir,:))

  end function getDensTrans

  subroutine integrateDensTrans(dr, D)
    real, intent(in) :: dr
    real, dimension(0:,:), intent(in) :: D
    integer :: ix,iy, nBin
    real :: x,y, S, S0, val

    nBin = ubound(D,dim=1)
    S = 0
    S0 = 0

!!$    do ix=0,nBin
!!$       write(101,*) ix*dr,D(ix,:)
!!$    end do

    do ix=-nBin,nBin
       x = ix * dr
       do iy=-nBin,nBin
          y = iy * dr

          val = getDensTrans(x,y, dr, D)
          S = S + val
          if (val>1e-3) S0 = S0 + 1

       end do
    end do

    write(*,*) 'Integral = ',S*dr**2


  end subroutine integrateDensTrans


  subroutine doCalc(b, nColl, nPart)

    real, intent(in) :: b
    real, intent(out) :: nColl, nPart

    real :: RA, RB
    integer :: NA,NB

    integer, parameter :: nBin = 100
    real :: dBin
    real, parameter :: valEps = 1e-10
    real :: xA, yA, valA, fakA
    real :: xB, yB, valB, fakB
    real :: SumColl, sumPart
    integer :: i,j

    RA = NucA%dx*NucA%maxIndex
    RB = NucB%dx*NucB%maxIndex

    NA = nucA%mass
    NB = nucB%mass

!!$    write(*,*) 'R: ',RA,RB
!!$    write(*,*) 'N: ',NA,NB

    dBin = RB/nBin

    fakA = sigma/NA*0.1 ! in fm^2
    fakB = sigma/NB*0.1 ! in fm^2

    SumColl = 0.0
    SumPart = 0.0

    do i=-nBin,nBin
       xB = i * dBin
       do j=-nBin,nBin
          yB = j * dBin

          valB = getDensTrans(xB,yB, NucB%dx, densTransB)
          if (valB < valEps) cycle

          xA = xB + b
          yA = yB

          valA = getDensTrans(xA,yA, NucA%dx, densTransA)
          if (valA < valEps) cycle

          SumColl = SumColl + valA*valB

          SumPart = SumPart &
               + valA * ( 1.0 - (1.0 - valB*fakB)**nB ) &
               + valB * ( 1.0 - (1.0 - valA*fakA)**nA )


       end do
    end do

    nColl = SumColl * dBin*dBin * sigma*0.1
    nPart = SumPart * dBin*dBin

  end subroutine doCalc

  real function impact_HADES (b, i)
    ! profile for the impact parameter distributions fitted to the HADES trigger
    ! for C+C see:    Agakishiev et al., EPJ A 40 (2009) 45, fig. 6
    ! for Ar+KCl see: Agakishiev et al., EPJ A 47 (2011) 21, fig. 1
    ! AuAu: Tetyana Galatyuk, private communications
    use distributions, only: woods_saxon
    use CallStack, only: traceBack

    real, intent(in) :: b      ! impact parameter in fm
    integer, intent(in) :: i   ! select reaction: 1 = CC1, 2 = CC2, 3 = ArKCl, 4 = AuAu (min.bias), 5 = AuAu (central)

    real, parameter :: b0(1:5) = (/ 3.705, 4.034, 4.922, 10.00, 4.506 /)  ! Woods-Saxon parameters for HADES centrality selection
    real, parameter :: db(1:5) = (/ 0.767, 0.749, 0.601, 0.800, 0.578 /)

    if (i<1 .or. i>5) then
       print *,"error in impact_HADES: ", i
       call traceBack()
    end if

    impact_HADES = woods_saxon (b, 1., b0(i), db(i))

  end function impact_HADES


end program OpticalGlauber
