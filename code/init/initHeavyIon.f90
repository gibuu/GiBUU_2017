!******************************************************************************
!****m* /initHeavyIon
! NAME
! module initHeavyIon
!
! PURPOSE
! Initialization of nucleus-nucleus collisions.
!
! INPUTS
! The namelist "heavyIon".
!******************************************************************************
module initHeavyIon

  implicit none
  private

  !****************************************************************************
  !****g* initHeavyIon/impact_parameter
  ! SOURCE
  !
  real, save :: impact_parameter = 0.
  ! PURPOSE
  ! Impact parameter b [fm]. There are three options:
  ! * b>=0: The impact parameter is fixed to the given value.
  ! * -100<b<0: The impact parameter will be chosen randomly in each run between 0 and abs(b).
  ! * b<=-100: "Minimum bias". The impact parameter will be chosen randomly in each run (maximum = sum of radii plus twice the sum of surfaces).
  !****************************************************************************


  !****************************************************************************
  !****g* initHeavyIon/impact_profile
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
  !****************************************************************************


  !****************************************************************************
  !****g* initHeavyIon/distance
  ! SOURCE
  !
  real, save :: distance = 0.
  ! PURPOSE
  ! Distance between centers of nuclei along z (i.e. beam)-direction [fm].
  ! This will be readjusted automatically in case it is too small.
  !****************************************************************************

  !****************************************************************************
  !****g* initHeavyIon/coulomb
  ! SOURCE
  !
  logical, save :: coulomb=.false.
  ! PURPOSE
  ! If .true., then a Coulomb propagation from coulombDistance = 10000 fm
  ! to distance is performed.
  !****************************************************************************

  !****************************************************************************
  !****g* initHeavyIon/ekin_lab_Projectile
  ! SOURCE
  !
  real, public, save :: ekin_lab_Projectile = 0.
  ! PURPOSE
  ! Kinetic energy per nucleon of projectile nucleus in lab frame [GeV].
  !****************************************************************************

  !****************************************************************************
  !****g* initHeavyIon/ekin_lab_Target
  ! SOURCE
  !
  real, save :: ekin_lab_Target = 0.
  ! PURPOSE
  ! Kinetic energy per nucleon of target nucleus in lab frame [GeV].
  !****************************************************************************

  !****************************************************************************
  !****g* initHeavyIon/adjustGridFlag
  ! SOURCE
  !
  logical, save ::  adjustGridFlag=.false.
  ! PURPOSE
  ! If .true., the grid spacing in z-direction will be readjusted.
  !****************************************************************************

  !****************************************************************************
  !****g* initHeavyIon/cmsFlag
  ! SOURCE
  !
  logical, save ::  cmsFlag=.true.
  ! PURPOSE
  ! If .true.,  the collision takes place in the CM frame of the two nuclei (default option).
  ! If .false., the collision takes place in the LAB frame (target at rest).
  !****************************************************************************

  real, public, save :: b   ! actual impact parameter

  public :: initHeavyIonCollision


contains


  !****************************************************************************
  !****s* initHeavyIon/initHeavyIonCollision
  ! NAME
  ! subroutine initHeavyIonCollision (targetNuc, projectileNuc)
  ! PURPOSE
  ! Provides initial configuration of a heavy ion event in the CM-frame
  ! INPUTS
  ! * type(tNucleus), pointer :: targetNuc     --- Target nucleus
  ! * type(tNucleus), pointer :: projectileNuc --- Projectile nucleus
  ! NOTES
  ! * Sets the kinematics of targetNuc and projectileNuc in their CM-frame
  ! * Expects that the masses (A,Z) of the nuclei are already well defined.
  !****************************************************************************
  subroutine initHeavyIonCollision (targetNuc, projectileNuc)
    use nucleusDefinition
    use Dilepton_Analysis, only: Dilep_Init
    use random, only: rn, rnGauss

    type(tNucleus), pointer :: targetNuc, projectileNuc

    real, dimension(1:3) :: betaCM2Lab           ! boost vector from CM to Lab frame
    real :: energyTargetCM,energyProjectileCM    ! energies in CM system
    real :: gammaProj, gammaTarg                 ! gamma factors in CM frame of the projectile and target
    real :: momentumCM, bmax
    logical, parameter :: debug=.true.
    logical, save :: first = .true.
    real, parameter :: b0(6:8) = (/ 5.67, 7.27, 8.59 /)
    real, parameter :: db(6:8) = (/ 0.99, 0.85, 0.79 /)

    if (first) call initInput

    if (impact_parameter >= 0.) then
      ! impact parameter is fixed
      b = impact_Parameter
    else
      ! choose random impact parameter in a given interval
      if (impact_parameter > -100.) then
        bmax = abs(impact_parameter)
      else
        bmax = targetNuc%radius + projectileNuc%radius + 2*(targetNuc%surface+projectileNuc%surface)
      end if

      select case (impact_profile)
      case (0)
        ! default: minimum bias
        b = bmax * sqrt(rn())
      case (1,2,3,4,5)
        ! use HADES trigger-biased centrality distributions with Woods-Saxon shape
        do
          b = bmax * sqrt(rn())
          if (rn() < impact_HADES(b,impact_profile)) exit  ! accept
        end do
      case (6,7,8)
        ! use HADES trigger-biased centrality distributions with Gaussian shape
        b = rnGauss(db(impact_profile), b0(impact_profile))
      case default
        write(*,*) 'initHeavyIonCollision: bad value for impact_profile: ', impact_profile
        stop
      end select

      write(*,*) 'initHeavyIonCollision: Impact Parameter = ',b
    end if

    call setVelocity
    call setPosition

    if (adjustGridFlag) call adjustGrid
    if (coulomb) call coulombCorrect

    call Dilep_Init (ekin_lab_Projectile + ekin_lab_Target, cms=cmsFlag, beta=betaCM2Lab)

  contains


    !**************************************************************************
    !****s* initHeavyIon/initInput
    ! NAME
    ! subroutine initInput
    ! PURPOSE
    ! Reads input out of jobcard. Namelist 'heavyIon'.
    !**************************************************************************
    subroutine initInput

      use output, only: Write_ReadingInput

      integer :: ios

      !************************************************************************
      !****n* initHeavyIon/heavyIon
      ! NAME
      ! NAMELIST heavyIon
      ! PURPOSE
      ! Includes parameters for initialization of a nucleus-nucleus collision:
      ! * impact_parameter    -- Impact Parameter [fm]
      ! * impact_profile      -- choose impact parameter distribution
      ! * distance            -- Distance between centers of nuclei along z (i.e. beam)-direction [fm]
      ! * coulomb             -- If .true., then a Coulomb propagation from coulombDistance = 10000 fm to distance is performed.
      ! * ekin_lab_Target     -- Kinetic energy per nucleon of target in lab frame [GeV]
      ! * ekin_lab_Projectile -- Kinetic energy per nucleon of projectile in lab frame [GeV]
      ! * adjustGridFlag      -- If .true., the grid spacing in z-direction will be readjusted.
      ! * cmsFlag             -- Perform the collision in the center-of-mass frame?
      !************************************************************************
      NAMELIST /heavyIon/ impact_parameter, impact_profile, distance, coulomb, &
                          ekin_lab_Target, ekin_lab_Projectile, adjustGridFlag, cmsFlag

      call Write_ReadingInput('heavyIon',0)
      rewind(5)
      read(5,nml=heavyIon,IOSTAT=ios)
      call Write_ReadingInput('heavyIon',0,ios)

      write(*,*) '  Impact Parameter=',impact_parameter
      write(*,*) '  Impact Profile  =',impact_profile
      write(*,*) '  Distance=',distance
      write(*,*) '  Coulomb correction for trajectories=',coulomb
      write(*,*) '  Kinetic Energy/nucleon of Target in lab frame    =',ekin_lab_Target
      write(*,*) '  Kinetic Energy/nucleon of Projectile in lab frame=',ekin_lab_Projectile
      write(*,*) '  Grid in z-direction will be adjusted ?:', adjustGridFlag
      if (cmsFlag) then
         write(*,*) '  calculational frame: CM'
      else
         write(*,*) '  calculational frame: LAB'
      end if

      call Write_ReadingInput('heavyIon',1)

      first = .false.

    end subroutine initInput


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
        write(*,*) "error in impact_HADES: ", i
        call traceBack()
      end if

      impact_HADES = woods_saxon (b, 1., b0(i), db(i))

    end function impact_HADES


    !**************************************************************************
    !****s* initHeavyIon/setVelocity
    ! NAME
    ! subroutine setVelocity
    ! PURPOSE
    ! Sets initial kinematics of the target and projectile nucleus (in the CM frame).
    !**************************************************************************
    subroutine setVelocity
      use constants, only: mN

      real :: energyTargetLab,energyProjectileLab  !energies in lab system
      real :: s  ! Mandelstam s
      real :: massTarget,massProjectile,momentumTargetLab,momentumProjectileLab


      s = (mN + ekin_lab_target + mN + ekin_lab_projectile)**2 &
           - ( sqrt((mN + ekin_lab_target)**2 - mN**2) &
           &   - sqrt((mN + ekin_lab_projectile)**2 -mN**2) )**2
      write(*,*) 'sqrt(s)_pp = ', sqrt(s)


      massTarget     = float(targetNuc%mass)     * mN
      massProjectile = float(projectileNuc%mass) * mN

      EnergyTargetLab     = float(targetNuc%mass)     * ekin_lab_target     + massTarget
      EnergyProjectileLab = float(projectileNuc%mass) * ekin_lab_projectile + massProjectile
      momentumTargetLab     = sqrt(EnergyTargetLab**2-massTarget**2)
      momentumProjectileLab = sqrt(EnergyProjectileLab**2-massProjectile**2)

      if (cmsFlag) then

         ! Projectile and Target are moving collinear, but in opposite direction, therefore:
         s = (EnergyTargetLab+EnergyProjectileLab)**2 &
              - (sqrt(EnergyTargetLab**2-massTarget**2) - sqrt(EnergyProjectileLab**2-massProjectile**2))**2

         momentumCM=Sqrt((s-(massTarget-massProjectile)**2)*(s-(massTarget+massProjectile)**2)/s)/2.
         if (debug) write(*,*) 'Momentum in CM = ', momentumCM
         EnergyTargetCM=SQRT(momentumCM**2+massTarget**2)
         EnergyProjectileCM=SQRT(momentumCM**2+massProjectile**2)
         targetNuc%velocity     = (/ 0., 0., -momentumCM/EnergyTargetCM     /)
         projectileNuc%velocity = (/ 0., 0.,  momentumCM/EnergyProjectileCM /)

      else

         if (debug) then
            write(*,*) 'Projectile momentum in LAB = ',momentumProjectileLab
            write(*,*) 'Target momentum in LAB     = ',momentumTargetLab
         end if
         targetNuc%velocity=(/ 0.,0.,-momentumTargetLab/EnergyTargetLab/)
         projectileNuc%velocity=(/0.,0.,momentumProjectileLab/EnergyProjectileLab/)

      end if

      betaCM2Lab = - (momentumProjectileLab-momentumTargetLab) / (EnergyProjectileLab+EnergyTargetLab) * (/0.,0.,1./)

      if (debug) then
        write(*,'(A,3G15.8)') ' Velocity Target     =', targetNuc%velocity
        write(*,'(A,3G15.8)') ' Velocity Projectile =', projectileNuc%velocity
        write(*,'(A,3G15.8)') ' betaCM2Lab          =', betaCM2Lab
      end if

      gammaTarg = 1./sqrt( 1. - targetNuc%velocity(3)**2 )
      gammaProj = 1./sqrt( 1. - projectileNuc%velocity(3)**2 )

    end subroutine setVelocity


    !**************************************************************************
    !****s* initHeavyIon/setPosition
    ! NAME
    ! subroutine setPosition
    ! PURPOSE
    ! Sets intitial positions of the target and projectile nucleus.
    ! NOTES
    ! In case of partial overlap of target and projectile,
    ! the input distance between the two nuclei is readjusted.
    !**************************************************************************
    subroutine setPosition

      real :: distance_touch, vP, vT

      ! check whether input distance is reasonable and readjust it if needed:
      distance_touch = (5.*targetNuc%surface+targetNuc%radius)/gammaTarg &
                     + (5.*projectileNuc%surface+projectileNuc%radius)/gammaProj + 1.

      write(*,*) 'Distance between centers along z-axis =', distance
      if (distance < distance_touch) then
         write(*,*)
         write(*,*) 'WARNING: Heavy Ions are sitting too close at initilization'
         write(*,*) 'Surface of target     =', targetNuc%surface
         write(*,*) 'Surface of projectile =', projectileNuc%surface
         write(*,*) '(Radius of target)/gammaTarg     =', targetNuc%radius / gammaTarg
         write(*,*) '(Radius of projectile)/gammaProj =', projectileNuc%radius / gammaProj
         distance = distance_touch
         write(*,*) 'Readjusted distance between centers along z-axis =', distance
      end if
      write(*,*) 'Overlap at t = ',distance/abs(targetNuc%velocity-projectileNuc%velocity)

      if (cmsFlag) then

         !Set position such that nuclei collide at center of grid (CM frame):
         vP = abs(projectileNuc%velocity(3))
         vT = abs(targetNuc%velocity(3))
         targetNuc%position     = (/ -b/2., 0.,  distance*vT/(vP+vT) /)
         projectileNuc%position = (/  b/2., 0., -distance*vP/(vP+vT) /)

      else

         !Set position such that target nucleus is at rest (LAB frame):
         targetNuc%position=(/ 0.,0.,0./)
         projectileNuc%position=(/ b,0.,-distance/)

      end if

    end subroutine setPosition


    !**************************************************************************
    !****s* initHeavyIon/adjustGrid
    ! NAME
    ! subroutine adjustGrid
    ! PURPOSE
    ! Adjustes the grid spacings to resolve the Lorentz-contracted radii along z-axis.
    ! NOTES
    ! The size of the grid is fixed.
    !**************************************************************************
    subroutine adjustGrid
      use densitymodule, only: acceptGrid

      real :: Rlong_min, Rtr_min!, Rlong_max, Rtr_max

      Rlong_min = min(targetNuc%radius/gammaTarg,projectileNuc%radius/gammaProj)
      !Rlong_max = max(targetNuc%radius/gammaTarg,projectileNuc%radius/gammaProj)

      Rtr_min = min(targetNuc%radius,projectileNuc%radius)
      !Rtr_max = max(targetNuc%radius,projectileNuc%radius)

      call acceptGrid(Rtr_min/6., Rtr_min/6., Rlong_min/6.)

    end subroutine adjustGrid


    !**************************************************************************
    !****s* initHeavyIon/coulombCorrect
    ! NAME
    ! subroutine coulombCorrect
    ! PURPOSE
    ! Do correction of the trajectories of the two nuclei due to Coulomb forces
    ! Treat the two nuclei like point particles
    !**************************************************************************
    subroutine coulombCorrect
      use coulombKorrektur, only: CoulpropaTwo
      use constants, only: mN

      real,parameter :: coulombDistance=10000. !Distance which is assumed to be infinity
      real,dimension(3) :: momentumTarget,momentumProjectile
      real :: massTarget, massProjectile, vP, vT

      write(*,*)
      write(*,*) '**Doing Coulomb correction to velocity and position of target and projectile.'
      write(*,*) ' We propagate them from',coulombdistance,'to', distance
      write(*,*) ' Positions before correction:'
      write(*,*) ' Target:',targetNuc%position
      write(*,*) ' Projectile:',projectileNuc%position
      write(*,*) ' Velocities before correction:'
      write(*,*) ' Target:',targetNuc%velocity
      write(*,*) ' Projectile:',projectileNuc%velocity

      massTarget     = float(targetNuc%mass)     * mN
      massProjectile = float(projectileNuc%mass) * mN

      momentumTarget     = (/ 0., 0., -momentumCM /)
      momentumProjectile = (/ 0., 0., +momentumCM /)
      distance           = Sqrt(distance**2+b**2)

      !Set total distance of both nuclei to "coulombDistance"
      vP = abs(projectileNuc%velocity(3))
      vT = abs(targetNuc%velocity(3))
      targetNuc%position     = (/  b/2., 0.,  coulombDistance*vT/(vP+vT) /)
      projectileNuc%position = (/ -b/2., 0., -coulombDistance*vP/(vP+vT) /)
      !Propagate both nuclei until they reach "distance
      call CoulpropaTwo (targetNuc%position, momentumTarget, targetNuc%charge, massTarget, &
                         projectileNuc%position, momentumProjectile, projectileNuc%charge, massProjectile, &
                         distance)

      EnergyTargetCM=SQRT(Dot_Product(momentumTarget,momentumTarget)+massTarget**2)
      EnergyProjectileCM=SQRT(Dot_Product(momentumProjectile,momentumProjectile)+massProjectile**2)
      targetNuc%velocity=momentumTarget/EnergyTargetCM
      projectileNuc%velocity=momentumProjectile/EnergyProjectileCM

      write(*,*)
      write(*,*) ' Positions after correction:'
      write(*,*) ' Target:',targetNuc%position
      write(*,*) ' Projectile:',projectileNuc%position
      write(*,*) ' Velocities after correction:'
      write(*,*) ' Target:',targetNuc%velocity
      write(*,*) ' Projectile:',projectileNuc%velocity

      write(*,*) '** Finished Coulomb correction to velocity and position of target and projectile.'


    end subroutine coulombCorrect

end subroutine initHeavyIonCollision

end module initHeavyIon
