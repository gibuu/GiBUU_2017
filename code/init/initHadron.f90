!******************************************************************************
!****m* /initHadron
! NAME
! module initHadron
!
! PURPOSE
! Initialises the hadron projectile for a hadron-nucleus collision.
! The eventtype 'hadron' always uses real particles.
!
! NOTES
! The target nucleus has to be initialised separately
! (by using initNucPhaseSpace routine) before calling the
! initHadronInduced routine.
!******************************************************************************
module initHadron

  implicit none
  private


  !****************************************************************************
  !****g* initHadron/impactParameter
  ! SOURCE
  !
  real, save :: impactParameter=0.
  ! PURPOSE
  ! smaller 0: Impact parameter will be chosen
  ! randomly in the interval [0;abs(impactParameter)]
  ! (see subroutine setGeometry).
  ! It is recommended to take very large negative value of impactParameter
  ! in order to have good automatic random choice, e.g. impactParameter=-100.
  !****************************************************************************


  !****************************************************************************
  !****g* initHadron/bRaiseFlag
  ! SOURCE
  !
  logical, save :: bRaiseFlag=.false.
  ! PURPOSE
  ! if .true.: actual impact parameter will be raised by deltaB after nRunsPerB
  ! subsequent runs. Starting value is given by the impactParameter variable.
  !****************************************************************************


  !****************************************************************************
  !****g* initHadron/deltaB
  ! SOURCE
  !
  real, save :: deltaB=0.
  ! PURPOSE
  ! impact parameter step (relevant if bRaiseFlag=.true.)
  !****************************************************************************


  !****************************************************************************
  !****g* initHadron/nRunsPerB
  ! SOURCE
  !
  integer, save :: nRunsPerB=1
  ! PURPOSE
  ! number of subsequent runs per impact parameter
  ! (relevant if bRaiseFlag=.true.)
  !****************************************************************************


  !****************************************************************************
  !****g* initHadron/perturbative
  ! SOURCE
  !
  logical, save :: perturbative=.false.
  ! PURPOSE
  ! if .true. then the hadron is a perturbative particle
  !****************************************************************************

  !****************************************************************************
  !****g* initHadron/numberParticles
  ! SOURCE
  !
  integer,save   :: numberParticles = 200
  !
  ! PURPOSE
  ! Number of projectile testparticles per ensemble
  ! in the case of a perturbative treatment
  !****************************************************************************


  !****************************************************************************
  !****g* initHadron/particleId
  ! SOURCE
  !
  integer, save :: particleId=1
  ! PURPOSE
  ! Identity of the projectile hadron.
  !****************************************************************************

  !****************************************************************************
  !****g* initHadron/antiParticle
  ! SOURCE
  !
  logical, save :: antiParticle=.false.
  ! PURPOSE
  ! if .true. then the hadron is an antiparticle
  !****************************************************************************


  !****************************************************************************
  !****g* initHadron/particleCharge
  ! SOURCE
  !
  integer, save :: particleCharge=0
  ! PURPOSE
  ! Charge of the hadron
  !****************************************************************************

  !****************************************************************************
  !****g* initHadron/ekin_lab
  ! SOURCE
  !
  real,    save :: ekin_lab=1.
  ! PURPOSE
  ! kinetic energy of the hadron in the rest frame of the target nucleus (GeV)
  ! NOTES
  ! If ekin_lab < 0.  --- initialization according to the binding energy
  !****************************************************************************

  !****************************************************************************
  !****g* initHadron/E_bind
  ! SOURCE
  !
  real,    save :: E_bind=0.
  ! PURPOSE
  ! binding energy of initialized hadron (GeV)
  ! NOTES
  ! Active for iniType= 0,2 if  ekin_lab < 0. is set.
  !****************************************************************************


  !****************************************************************************
  !****g* initHadron/iniType
  ! SOURCE
  !
  integer, save :: iniType=0
  ! PURPOSE
  ! * 0:   usual initialization for the hadron-nucleus collision
  ! * 1:   position and momentum of the hadron is chosen according to
  !        the Gaussians centered, resp., at the centre of the nucleus and at zero momentum
  !        (impactParameter, distance and ekin_lab have no effect in this case)
  ! * 2:   gaussian in coordinate space, but usual sharp momentum choice
  !        (impactParameter, distance and ekin_lab work as usual)
  !****************************************************************************

  !****************************************************************************
  !****g* initHadron/zChoice
  ! SOURCE
  !
  integer, save :: zChoice=1
  ! PURPOSE
  ! * 1:    hadron is initialised at fixed distance delta from nuclear surface
  ! * 2:    hadron is initialised at fixed z
  ! Relevant for iniType=0 or iniType=2.
  !****************************************************************************



  !****************************************************************************
  !****g* initHadron/delta
  ! SOURCE
  !
  real, save :: delta=0.5
  ! PURPOSE
  ! * for zChoice=1: distance from nuclear surface, at which the hadron is
  !   initialised [fm]
  ! * for zChoice=2: maximum distance from the edge of nucleus in transverse
  !   direction which restricts the choice of actual impact parameter for
  !   impactParameter < 0 (for impactParameter > 0 no restriction)
  ! Relevant for iniType=0 or  iniType=2.
  !****************************************************************************


  !****************************************************************************
  !****g* initHadron/deltaZ
  ! SOURCE
  !
  real, save :: deltaZ=5.
  ! PURPOSE
  ! z = -deltaZ - R_nucleus, where z is z-coordinate of the hadron
  ! Relevant for iniType=0,2 and zChoice=2.
  !****************************************************************************


  !****************************************************************************
  !****g* initHadron/width
  ! SOURCE
  !
  real, save ::    width=1.
  ! PURPOSE
  ! Width of a gaussian density profile [fm].
  ! Only relevant for iniType= 1 and 2.
  !****************************************************************************

  !****************************************************************************
  !****g* initHadron/debug
  ! SOURCE
  logical, parameter :: debug=.false.
  ! PURPOSE
  ! if .true. then additional printouts are done for debugging
  !****************************************************************************


  ! Working variables:
  real, save, public :: b, phi, z, p_lab

  public :: initHadronInduced
  public :: particleId,antiparticle,particleCharge,perturbative,E_bind

contains

  !****************************************************************************
  !****s* initHadron/initHadronInduced
  ! NAME
  ! subroutine initHadronInduced (teilchenReal, teilchenPert)
  ! PURPOSE
  ! Provides initial conditions for a hadron.
  ! INPUTS
  ! * type(particle),dimension(:,:),intent(inout) ::  teilchenReal ---
  !   array of real particles to store the hadron if it is real
  ! * type(particle),dimension(:,:),intent(inout) ::  teilchenPert ---
  !   array of perturbative particles to store the hadron if it is perturbative
  ! NOTES
  ! * The user has to provide the namelist 'hadron', which contains
  !   the impact parameter, lab. energy, id, charge and antiparticle-
  !   parameters of the hadron.
  !   If the input impact parameter is negative, then
  !   the actual impact parameter is chosen by Monte-Carlo between
  !   0 and abs(input impact parameter).
  ! * The target nucleus has to be initialised before calling this routine.
  ! * The program initialises the real (not perturbative) particle.
  !****************************************************************************
  subroutine initHadronInduced (teilchenReal, teilchenPert)

  use particleDefinition
  use nucleusDefinition
  use nucleus, only: getTarget, getProjectile
  use collisionNumbering, only: pert_Numbering

  type(particle), dimension(:,:), intent(inout) :: teilchenReal, teilchenPert

  type(tNucleus),pointer :: targetNuc                      ! Target nucleus
  type(tNucleus),pointer :: proj                           ! projectile (for storing velocity only)

  ! Working variables:
  type(particle) :: teilchen
  logical, save :: initFlag=.true.
  integer, save :: nRun=0
  real, save, dimension(1:3) :: beta
  integer :: i, j, index

  write(*,*)
  write(*,*) '**Initializing  hadron'

  if (initFlag) call initInput

  targetNuc => getTarget()
  proj => getProjectile()

  nRun=nRun+1

  if (iniType.eq.0 .or.iniType.eq.2) call setGeometry(.true.)

  if (.not.perturbative) then

      do j=1,size(teilchenReal,dim=1) ! Loop over all Ensembles

           index=1
           do !Search for empty space in particle vector
              if (index.gt.size(teilchenReal,dim=2)) then
                 write(*,*) 'Real particle vector too small. Stop in initHadronInduced.'
                 stop
              else if (teilchenReal(j,index)%ID > 0) then
                 index=index+1
              else
                 exit
              end if
           end do

           if (debug) write(*,*) ' In initHadronInduced real, index=', index

           call setToDefault(teilchen) !set teilchen to its default values
           call setPosition
           call setKinematics
           if (index.gt.1) then
              teilchen%event=teilchenReal(j,index-1)%event + 1
           else
              teilchen%event=1
           end if

           ! Give the particle its unique number
           call setNumber(teilchen)

           teilchenReal(j,index)=teilchen

      end do

  else

      do j=1,size(teilchenPert,dim=1) ! Loop over all Ensembles

           index=1
           do i=1,numberParticles
               do !Search for empty space in particle vector
                  if (index.gt.size(teilchenPert,dim=2)) then
                     write(*,*) 'Pert particle vector too small. Stop in initHadronInduced.'
                     stop
                  else if (teilchenPert(j,index)%ID > 0) then
                     index=index+1
                  else
                     exit
                  end if
               end do

               if (debug) write(*,*) ' In initHadronInduced pert, index=', index

               call setToDefault(teilchen) !set teilchen to its default values
               if (impactParameter.lt.0.) call setGeometry(.false.)
               call setPosition
               call setKinematics

               call setNumber(teilchen)
               teilchen%event(1:2)=pert_numbering()
               teilchen%perturbative=.true.
               teilchen%perweight=1./float(numberParticles)
               teilchenPert(j,index)=teilchen

           end do

      end do

  end if

  proj%velocity=beta

  contains

    !**************************************************************************
    !****s* initHadronInduced/initInput
    ! NAME
    ! subroutine initInput
    ! PURPOSE
    ! Reads input out of jobcard. Namelist 'hadron'.
    !**************************************************************************
    subroutine initInput

      use output, only: Write_ReadingInput
      use Dilepton_Analysis, only: Dilep_Init
      use checks, only: ChecksSwitchRealRun

      !************************************************************************
      !****n* initHadron/hadron
      ! NAME
      ! NAMELIST hadron
      ! PURPOSE
      ! Includes parameters for initialization of a hadron in the case of
      ! the hadron-nucleus collision:
      ! * impactParameter
      ! * bRaiseFlag
      ! * deltaB
      ! * nRunsPerB
      ! * perturbative
      ! * numberParticles
      ! * particleId
      ! * antiParticle
      ! * particleCharge
      ! * ekin_lab
      ! * E_bind
      ! * iniType
      ! * zChoice
      ! * delta
      ! * deltaZ
      ! * width
      !************************************************************************
      NAMELIST /hadron/ impactParameter, bRaiseFlag, deltaB, nRunsPerB, &
                        perturbative, numberParticles,                  &
                        particleId, antiparticle, particleCharge,       &
                        ekin_lab, E_bind, iniType, zChoice, delta,      &
                        deltaZ, width

      integer :: ios

      call Write_ReadingInput('hadron',0)
      rewind(5)
      read(5,nml=hadron,IOSTAT=ios)
      call Write_ReadingInput('hadron',0,ios)

      write(*,*) '  Impact Parameter=',impactParameter
      write(*,*) '  bRaiseFlag=',bRaiseFlag
      write(*,*) '  deltaB=',deltaB
      write(*,*) '  nRunsPerB=',nRunsPerB
      write(*,*) '  perturbative=',perturbative
      write(*,*) '  numberParticles=',numberParticles
      write(*,*) '  Id of a hadron=',particleId
      write(*,*) '  antiparticle=',antiparticle
      write(*,*) '  charge of a hadron=',particleCharge
      write(*,*) '  Kinetic Energy of a hadron in lab frame=',ekin_lab
      if (ekin_lab<0.) then
         write(*,*) '  *** Attention: initialization according to the binding energy of a hadron:'
         write(*,*) '  E_bind=', E_bind
      end if
      write(*,*) '  Initialization type, 0 -- usual, 1,2 -- gaussians =',iniType
      write(*,*) '  zChoice type, 1 -- fixed dist. from surf., 2 --- fixed z =', zChoice
      write(*,*) '  distance from surface =', delta
      write(*,*) '  shift from nuclear surface in -z direction =', deltaZ
      write(*,*) '  Gaussian density profile width of a hadron=',width

      call ChecksSwitchRealRun(.not.perturbative)

      call Write_ReadingInput('hadron',1)

      call Dilep_Init (ekin_lab)

      initFlag=.false.

    end subroutine initInput


    !**************************************************************************
    !****s* initHadronInduced/setKinematics
    ! NAME
    ! subroutine setKinematics
    ! PURPOSE
    ! Sets basic kinematics of a hadron colliding with a nucleus.
    !**************************************************************************
    subroutine setKinematics

      use particleProperties, only: hadron
      use constants, only: hbarc
      use random, only: rn
      use RMF, only: getRMF_flag, ModificationFactor, g_rho, g_omega, g_sigma
      use densitymodule, only: gridSpacing, gridPoints, rhoField, omegaField, sigmaField
      use IdTable, only: nucleon
      use coulomb, only: emfoca

      real :: E, width_p, pAbs, P_cut, fact, isofact, mstar, V0, Estar
      integer :: I1,I2,I3

      teilchen%ID=particleId
      teilchen%antiparticle=antiparticle
      teilchen%charge=particleCharge
      teilchen%mass=hadron(particleId)%mass

      select case (iniType)

      case (0,2)    ! Usual initialisation in momentum space:

         mstar = teilchen%mass
         V0 = 0.

         if ( getRMF_flag() ) then ! Determine vector and scalar potentials

            fact = ModificationFactor(teilchen%Id,teilchen%antiparticle)

            ! position in large grid:
            I1=NINT(teilchen%position(1)/gridSpacing(1))
            I2=NINT(teilchen%position(2)/gridSpacing(2))
            I3=NINT(teilchen%position(3)/gridSpacing(3))

            if (       abs(I1).le.gridPoints(1) &
              & .and. abs(I2).le.gridPoints(2) &
              & .and. abs(I3).le.gridPoints(3)  ) then

               mstar = teilchen%mass &
                    &+ fact*g_sigma*sigmaField(I1,I2,I3)
               if (mstar.le.0.) then
                 write(*,*) ' problems in initHadron, mstar: ', mstar
                 stop
               end if

               if (teilchen%ID==nucleon) then
                  if (teilchen%charge==0) then
                     isofact=-1.
                  else
                     isofact=1.
                  end if
               else
                  isofact=0. !isospin sector presently only only for protons and neutrons
               end if

               if ( teilchen%antiparticle ) fact=-fact

               V0 = g_omega*omegaField(I1,I2,I3,0)*fact &
                 &+ isofact*fact*g_rho*rhoField(I1,I2,I3,0)

            end if

            V0 = V0 + emfoca(teilchen%position,(/0.,0.,0./),teilchen%charge,teilchen%ID)
            !write(*,*) 'initHadronInduced/setKinematics V0, S: ', &
            !         & V0, mstar - teilchen%mass

         end if

         if (ekin_lab.ge.0.) then

            p_lab=Sqrt( (ekin_lab+teilchen%mass)**2 &
                      &-teilchen%mass**2)
            Estar=sqrt(p_lab**2+mstar**2)
            E=Estar+V0
            E_bind=teilchen%mass-E

         else ! Determine particle momentum according to the binding energy

            if ( .not.getRMF_flag() ) then
              write(*,*) ' ekin_lab < 0 is possible only in RMF mode currently'
              stop
            end if
            Estar=teilchen%mass-E_bind-V0
            if (Estar.le.mstar) then
               write(*,*) ' problems in initHadron: too strong binding'
               write(*,*) ' Estar, mstar: ', Estar, mstar
               stop
            end if
            p_lab=sqrt( Estar**2 - mstar**2 )
            E=sqrt(p_lab**2+mstar**2)+V0

         end if

         !write(*,*) 'initHadronInduced/setKinematics p_lab:', p_lab
         !write(*,*) 'initHadronInduced/setKinematics E_bind:', E_bind

         if (Estar.le.0.) then
             write(*,*) ' Problem in initHadron/setKinematics, Estar=',Estar
             stop
         end if

         teilchen%momentum(0)=Estar
         ! the hadron is initialized moving in positive z-direction:
         teilchen%momentum(1:3)=(/0.,0.,p_lab/)
         beta(1:3)=(/0.,0.,p_lab/Estar/)
         if (abs(beta(3)).gt.1.) then
            write(*,*) ' Problem in initHadron/setKinematics, beta(3)=',beta(3)
            stop
         end if

         if (iniType==2) &  ! Lorentz-contracted Gaussian wave packet:
           teilchen%position(3) = (teilchen%position(3)-z)*sqrt(1.-beta(3)**2) + z

      case (1)  ! Initialisation using Gaussian wave packet in momentum space:

         ! Monte Carlo sampling according to the coherent state momentum distribution:
         width_p=0.5*hbarc/width
         P_cut=5.*width_p
         do
           do
             teilchen%momentum(1)=(1.-2.*rn())*P_cut
             teilchen%momentum(2)=(1.-2.*rn())*P_cut
             teilchen%momentum(3)=(1.-2.*rn())*P_cut
             pAbs=Sqrt( dot_product(teilchen%momentum(1:3),teilchen%momentum(1:3)) )
             if (pAbs.le.P_cut) exit
           end do
           if ( rn() <= exp(-(pAbs/width_p)**2/2.) ) exit
         end do

         teilchen%momentum(0)=sqrt(teilchen%mass**2+pAbs**2)
         beta=0.

      case default

         write(*,*) 'Not valid initialisation type:',iniType,'STOP'
         Stop

      end select

      teilchen%velocity(1:3)=teilchen%momentum(1:3)/teilchen%momentum(0)

      if (debug) then
        write(*,*) 'Id of a hadron:',teilchen%Id
        write(*,*) 'Charge:',teilchen%charge
        write(*,*) 'Antiparticle ?:',teilchen%antiparticle
        write(*,*) 'Mass:',teilchen%mass
        write(*,*) '4-momentum:',teilchen%momentum
        write(*,*) 'Velocity:',teilchen%velocity
        write(*,*) 'Boost velocity:', beta
        write(*,*) 'Gamma factor:', 1./sqrt(1.-beta(3)**2)
      end if

    end subroutine setKinematics


    !**************************************************************************
    !****s* initHadronInduced/setPosition
    ! NAME
    ! subroutine setPosition
    ! PURPOSE
    ! Sets position of the hadron.
    ! NOTES
    ! In the case of usual initialisation (initype=0),
    ! the actual impact parameter b and the coordinate z has to be defined
    ! before by  calling the subroutine setGeometry.
    !**************************************************************************
    subroutine setPosition

      use random, only: rn

      real :: rAbs, R_cut            ! for initialisation using Gaussian (iniType=1,2)
      real, dimension(1:3) :: place


      select case (iniType)

      case (0)  ! Usual initialisation in coordinate space:

          teilchen%position=(/b*cos(phi),b*sin(phi),z/)

          if (debug) write(*,*) 'Position of hadron:',teilchen%position

      case (1,2) ! Initialisation using Gaussian wave packet in coordinate space:

          ! Monte Carlo distribution according to the gaussian density profile:
          R_cut=5.*width
          do
            do
              teilchen%position(1)=(1.-2.*rn())*R_cut
              teilchen%position(2)=(1.-2.*rn())*R_cut
              teilchen%position(3)=(1.-2.*rn())*R_cut
              rAbs=Sqrt( dot_product(teilchen%position(1:3),teilchen%position(1:3)) )
              if (rAbs.le.R_cut) exit
            end do
            if ( rn() <= exp(-(rAbs/width)**2/2.) ) exit
          end do

          if (iniType.eq.2) then
             teilchen%position(3)=teilchen%position(3)+z
             teilchen%position(1)=teilchen%position(1)+b
             place=teilchen%position
             teilchen%position(1)=place(1)*cos(phi)-place(2)*sin(phi)
             teilchen%position(2)=place(1)*sin(phi)+place(2)*cos(phi)
          end if

      case default

         write(*,*) 'Not valid initialisation type:',iniType,'STOP'
         Stop

      end select

    end subroutine setPosition


    !**************************************************************************
    !****s* initHadronInduced/setGeometry
    ! NAME
    ! subroutine setGeometry (flagPrint)
    ! PURPOSE
    ! Sets the impact parameter b and coordinate z of the hadron.
    ! INPUTS
    ! * logical, intent(in) :: flagPrint --  if .true., then print output parameters
    !   to standard output file
    ! NOTES
    ! If the input impactParameter is less than zero, than the
    ! actual impact parameter is choosen by a Monte-Carlo-decision.
    !**************************************************************************
    subroutine setGeometry (flagPrint)

      use random, only: rn
      use constants, only: pi

      logical, intent(in) :: flagPrint     !  if .true., then print output parameters
                                           !  into standart output file

      real :: Radial_distance

      Radial_distance = targetNuc%radius + delta

      if (impactParameter.ge.0.) then
        b = impactParameter
        if (bRaiseFlag .and. mod(nRun,nRunsPerB).eq.0) impactParameter= impactParameter + deltaB
        phi=0.
      else
        b = min(abs(impactParameter),Radial_distance)*sqrt(rn())
        phi=2.*pi*rn()
      end if

      select case (zChoice)

      case (1) ! The hadron is initialised at fixed distance delta from the nuclear surface:

         if (b.gt.Radial_distance) then
            write(*,*) ' Warning in setGeometry: too big impactParameter. b, Radial_distance: ',&
                      & b, Radial_distance
            b=Radial_distance-0.01
         end if
         z = -sqrt( max(0.,Radial_distance**2-b**2) )

      case (2) ! The hadron is initialised at fixed z:

         z = -targetNuc%radius - deltaZ

      case default

         write(*,*) ' in setGeometry: wrong zChoice=', zChoice
         stop

      end select

      if (flagPrint) then
          write(*,*) 'initHadronInduced/setGeometry Radial_distance [fm]:', Radial_distance
          write(*,*) 'initHadronInduced/setGeometry b, z [fm]:', b, z
      end if

    end subroutine setGeometry


  end subroutine initHadronInduced


end module initHadron
