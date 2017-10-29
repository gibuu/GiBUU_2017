!******************************************************************************
!****m* /initPion
! NAME
! module initPion
! PURPOSE
! Includes the initialization of pions for pion induced events.
!******************************************************************************
module initPion

  implicit none

  private
  public :: initPionInduced, getEkin, getTotalPerweight


  !Input variables out of jobcard:

  !****************************************************************************
  !****g* pionNucleus/UseCoulomb
  ! SOURCE
  logical,save   :: UseCoulomb=.false.
  ! PURPOSE
  ! if .true. then a Coulomb propagation from CoulombDistance to distance is
  ! performed
  !****************************************************************************

  !****************************************************************************
  !****g* pionNucleus/CoulombDistance
  ! SOURCE
  real,save      :: CoulombDistance=200. ! [fm]
  ! PURPOSE
  ! distance from where the Coulomb propagation starts
  !****************************************************************************

  !****************************************************************************
  !****g* pionNucleus/distance
  ! SOURCE
  real,save      :: distance=15. ! [fm]
  ! PURPOSE
  ! initialization distance
  !****************************************************************************

  !****************************************************************************
  !****g* pionNucleus/impact_parameter
  ! SOURCE
  real,save      :: impact_parameter=0. ! [fm]
  ! PURPOSE
  ! impact parameter.
  ! If less than 0, than an impact parameter integration is performed
  !****************************************************************************

  !****************************************************************************
  !****g* pionNucleus/charge
  ! SOURCE
  integer,save   :: charge=0
  ! PURPOSE
  ! charge of pion
  !****************************************************************************

  !****************************************************************************
  !****g* pionNucleus/numberPions
  ! SOURCE
  integer,save   :: numberPions=200
  ! PURPOSE
  ! number of initialized pions per ensemble
  !****************************************************************************

  !****************************************************************************
  !****g* pionNucleus/ekin_lab
  ! SOURCE
  real,save      :: ekin_lab=0.
  !
  ! PURPOSE
  ! kinetic energies of pions in lab frame.
  !****************************************************************************

  !****************************************************************************
  !****g* pionNucleus/delta_ekin_lab
  ! SOURCE
  real,save      :: delta_ekin_lab=0.01  ! [GeV]
  ! PURPOSE
  ! step size for kinetic energies in energy scans
  !****************************************************************************

  !****************************************************************************
  !****g* pionNucleus/debug
  ! SOURCE
  logical, parameter :: debug=.false.
  ! PURPOSE
  ! Switch on debug mode
  !****************************************************************************

  real,save      :: totalPerweight=0.


contains

  !****************************************************************************
  !****f* initPion/getEkin
  ! NAME
  ! real function getEkin
  ! PURPOSE
  ! This function returns the kinetic energy of the pions as they were
  ! initialized at the last call of initPionInduced.
  !****************************************************************************
  real function getEkin()
    getEkin=ekin_lab
  end function getEkin


  !****************************************************************************
  !****f* initPion/getTotalPerweight
  ! NAME
  ! real function getTotalPerweight
  ! PURPOSE
  ! This function returns the total perweight of the pions as they were
  ! initialized at the last call of initPionInduced.
  !****************************************************************************
  real function getTotalPerweight()
    getTotalPerweight=totalPerweight
  end function getTotalPerweight

  !****************************************************************************
  !****s* initPion/initPionInduced
  ! NAME
  ! subroutine initPionInduced(teilchen,raiseEnergyFlag,targetNuc)
  ! PURPOSE
  ! This routine initializes pions in pion-nucleus scattering.
  ! INPUTS
  ! * type(particle), intent(inout), dimension(:,:) :: teilchen
  !   -- vector to store pions in
  ! * logical, intent(in) :: raiseEnergyFlag
  !   -- if .true. energy of initialized pions is raised by delta_ekin_lab
  ! * type(tNucleus),pointer,optional :: targetNuc -- Target nucleus
  !****************************************************************************
  subroutine initPionInduced(teilchen,raiseEnergyFlag,targetNuc)
    use idTable
    use particleDefinition
    use nucleusDefinition
    use collisionNumbering, only:pert_Numbering ! Numbering of %event of the perturbative particles
    use output

    type(particle), intent(inout), dimension(:,:) :: teilchen
    logical, intent(in) :: raiseEnergyFlag
    type(tNucleus),pointer,optional :: targetNuc

    integer  :: index !of pion in vector teilchen
    integer  :: i,j   !loop indizes
    logical,save  :: outside
    logical  :: successFlag


    !integer :: numEnsembles
    logical, save :: initFlag=.true.

    !numEnsembles=size(teilchen,dim=1)

    write(*,*)
    write(*,subchapter) 'Initializing  pion induced events'

    if (initFlag) then
       !Read input and check whether pion is initialized in- or outside nucleus
       call initInput(outside)
       initFlag=.false.
    end if

    if (raiseEnergyFlag) ekin_lab=ekin_lab+delta_ekin_lab
    write(*,*) ' Kinetic energy of pions in lab frame=', ekin_lab
    write(*,*) ' Outside-Flag=', outside
    totalPerweight=0.

    do j=1,size(teilchen,dim=1) ! Loop over all Ensembles
       index=1
       i=0
       do
          do while(teilchen(j,index)%Id > 0) !Find free place in particle vector
             index=index+1
             if (index.gt.size(teilchen,dim=2)) then
                write(*,*) 'Particle vector too small in initPion'
                write(*,*) 'Size=',size(teilchen,dim=2)
                write(*,*) 'Ensemble: ',j
                write(*,*) 'Number pions per ensemble=',numberPions
                stop
             end if
          end do

          call setToDefault(teilchen(j,index))
          call setKinematics
          call setPosition
          ! Give the particle its unique number
          call setNumber(teilchen(j,index))

          !Do correction of kinematics if pion is initialized outside the nucleus:
          if (outside.and.UseCoulomb) call CoulombCorrect

          !Do correction of momentum if pion is initialized inside the nucleus:
          if (.not.outside) then
             call momentumCorrect(successFlag)
             if (.not.successFlag) then
                teilchen(j,index)%ID=0
                write(*,*) 'Generating event not succesful:',j,i
                cycle ! New Event!
             end if
          end if
          i=i+1
          teilchen(j,index)%firstevent=i+size(teilchen,dim=1)*j
          if (i.eq.numberpions) exit
       end do
    end do

    write(*,*) '**Finished Initializing pions for pion induced events'
    write(*,*) '**Total perweight=', totalPerweight

    write(*,*)

  contains

    !**************************************************************************
    !****s* initPionInduced/initInput
    ! NAME
    ! subroutine initInput(outside)
    ! PURPOSE
    ! Reads input out of jobcard. Namelist 'pionNucleus'.
    ! RESULT
    ! logical outside !whether pions are initialized in or outside the nucleus
    ! NOTES
    ! Checks wether pion is initialized in- or outside the nucleus
    !**************************************************************************
    subroutine initInput(outside)
      use output

      logical, intent(out) :: outside

      !************************************************************************
      !****n* initPion/pionNucleus
      ! NAME
      ! NAMELIST pionNucleus
      ! PURPOSE
      ! Includes parameters of pion initialization:
      ! * UseCoulomb
      ! * CoulombDistance
      ! * distance
      ! * impact_parameter
      ! * charge
      ! * numberPions
      ! * ekin_lab
      ! * delta_ekin_lab
      !************************************************************************

      NAMELIST /pionNucleus/ UseCoulomb,CoulombDistance,distance,&
           & impact_parameter,&
           & charge,numberPions,ekin_lab,delta_ekin_lab

      call Write_ReadingInput('pionNucleus',0)
      rewind(5)
      read(5,nml=pionNucleus)
      write(*,*) '  Impact Parameter=',impact_parameter
      write(*,*) '  Distance=',distance
      write(*,*) '  Coulomb correction for trajectories=',UseCoulomb
      write(*,*) '  Distance for the Coulomb correction=',CoulombDistance
      write(*,*) '  Kinetic Energy of pions in lab frame=',ekin_lab
      write(*,*) '  Delta(Energy) for energy scans =',delta_ekin_lab
      write(*,*) '  Number of pions per ensemble=',numberPions
      write(*,*) '  Charge of pions=',charge
      write(*,*)
      if (present(targetNuc)) then
         if (Sqrt(distance**2+impact_parameter**2).lt.(targetNuc%radius+targetNuc%surface)) then
            outside=.false.
            write(*,*) 'Pions are initialized inside the nucleus'
         else
            outside=.true.
            write(*,*) 'Pions are initialized outside the nucleus'
         end if
         write(*,*)
      else
         outside=.false.
         write(*,*) 'Pions are initialized inside the nucleus'
      end if
      call Write_ReadingInput('pionNucleus',1)

    end subroutine initInput

    !**************************************************************************
    !****s* initPionInduced/setKinematics
    ! NAME
    ! subroutine setKinematics
    ! PURPOSE
    ! Sets basic kinematics of the pions.
    !**************************************************************************
    subroutine setKinematics

      use constants, only: mPi

      Teilchen(j,index)%ID=pion
      Teilchen(j,index)%charge=charge
      Teilchen(j,index)%antiparticle=.false.
      Teilchen(j,index)%perturbative=.true.
      Teilchen(j,index)%productionTime=0.
      Teilchen(j,index)%mass=mPi
      Teilchen(j,index)%momentum(0)=ekin_lab+Teilchen(j,index)%mass
      !pion is initialized moving in positive z-direction:
      Teilchen(j,index)%momentum(1:3)=(/0.,0.,Sqrt(Teilchen(j,index)%momentum(0)**2-Teilchen(j,index)%mass**2)/)
      !assume vacuum dispersion relation:
      Teilchen(j,index)%velocity(1:3)=teilchen(j,index)%momentum(1:3)/teilchen(j,index)%momentum(0)
      teilchen(j,index)%event(1:2)=pert_numbering()
      if (debug) write(*,*) 'Masse=',teilchen(j,index)%mass
    end subroutine setKinematics

    !**************************************************************************
    !****s* initPion/setPosition
    ! NAME
    ! subroutine setPosition
    ! PURPOSE
    ! Sets positions of the pions.
    ! NOTES
    ! If Impact_Parameter is choosen to be less than zero, than the impact
    ! parameter is choosen by a Monte-Carlo-decision. This is made such that
    ! the pion is initialized on a disk of radius "bmax_Innerdisk" or on a
    ! ring which surrounds the inner disk and has an outer radius of
    ! "bmaxOuterRing".
    ! The probability to be on the inner disk is given by "pInnerDisk".
    ! The inner disk and the outer ring are separetely populated by a constant
    ! number density of pions.
    ! One distinguishes between inner disk and outer ring to have the
    ! possibility to have different
    ! population densities. Assumed one would only have one disk, then most of
    ! the particles would be
    ! initialized with high impact-parameter where only few reactions take
    ! place.
    !
    ! The perweight is given in units of mb for impact parameter integration.
    !**************************************************************************
    subroutine setPosition

      use random
      use inputGeneral
      use constants

      real :: bmax_OuterRing              !maximal Radius of outer ring
      real :: bmax_InnerDisk              !Radius of inner ring
      real, parameter :: pInnerDisk=0.7   !probability for initialization on inner ring
      real :: minimalDistance=2.52        !=pirp in old BUU
      !SQRT(maximal crossection of pion and nucleus/pi)
      real :: phi

      real, parameter :: ratioRadius=1.8 ! bmax_Outerring=ratioRadius*nuclearRadius+...

      integer :: totalNumPions
      real :: randomNumber
      real :: impact
      logical, save :: flag = .true.


      totalNumPions=numberPions*numEnsembles


      if (impact_parameter.ge.0.) then
         teilchen(j,index)%position=(/impact_Parameter,0.,-distance/)

         !         if(fullensemble) then
         !            teilchen(j,index)%perweight=1./float(numberPions)
         !         else
         teilchen(j,index)%perweight=1./float(totalNumPions)
         !         end if
      else  !Monte Carlo decision to have impact parameter integration in the end
         !maximum impact parameter of outer ring:
         if (fullEnsemble) then       !Full ensemble
            minimalDistance=minimalDistance/sqrt(float(numEnsembles))
         end if
         if (present(targetNuc)) then
            if (targetNuc%radius.gt.0.001) then  !No elementary event
               !maximum impact parameter of outer ring:
               bmax_OuterRing=ratioRadius*targetNuc%radius+minimalDistance

               !maximum impact parameter of inner disk:
               bmax_InnerDisk=targetNuc%radius+minimaldistance
            else
               !maximum impact parameter of outer ring:
               bmax_OuterRing=3.

               !maximum impact parameter of inner disk:
               bmax_InnerDisk=2.
               if (fullEnsemble) then       !Full ensemble
                  bmax_InnerDisk=bmax_InnerDisk  /sqrt(float(numEnsembles))
                  bmax_OuterRing=bmax_OuterRing  /sqrt(float(numEnsembles))
               end if
            end if
         else
            !maximum impact parameter of outer ring:
            bmax_OuterRing=3.

            !maximum impact parameter of inner disk:
            bmax_InnerDisk=2.
         end if
         if (flag) then
            write(*,*) '  Radius of outer ring:'  ,bmax_OuterRing
            write(*,*) '  Radius of inner circle:',bmax_InnerDisk
            write(*,*) '  perweight for pion in inner circle in fm^2:', (pi*bmax_InnerDisk**2)/pInnerDisk/float(totalNumPions)

            write(*,*) '  perweight for pion in outer ring in fm^2:', &
                 &        pi*(bmax_OuterRing**2-bmax_InnerDisk**2)/(1.-pInnerDisk)/float(totalNumPions)
            flag=.false.
         end if

         randomNumber=rn()
         phi=rn()*2*pi
         if (randomNumber.le.pInnerDisk) then ! impact parameter within nuclear radius
            impact=rn()
            teilchen(j,index)%position(1)=sqrt(impact)*bmax_InnerDisk*cos(phi)
            teilchen(j,index)%position(2)=sqrt(impact)*bmax_InnerDisk*sin(phi)
            teilchen(j,index)%position(3)=-distance
            teilchen(j,index)%perweight=(pi*bmax_InnerDisk**2)/pInnerDisk /float(totalNumPions)*10 ! in mB (factor 10 due to fm**2 to mb conversion)
         else                          !impact parameter not within nuclear radius
            impact=rn()
            teilchen(j,index)%position(1)=sqrt(impact*(bmax_OuterRing**2-bmax_InnerDisk**2)&
                 &                                      +bmax_InnerDisk**2)*cos(phi)
            teilchen(j,index)%position(2)=sqrt(impact*(bmax_OuterRing**2-bmax_InnerDisk**2)&
                 &                                      +bmax_InnerDisk**2)*sin(phi)
            teilchen(j,index)%position(3)=-distance

            teilchen(j,index)%perweight=pi*(bmax_OuterRing**2-bmax_InnerDisk**2)/(1.-pInnerDisk) &
                 &       /float(totalNumPions)*10  ! in mB (factor 10 due to fm**2 to mb conversion)
         end if
         !        if(fullensemble) teilchen(j,index)%perweight=teilchen(j,index)%perweight*float(numEnsembles)
      end if
      totalPerweight=totalPerweight+teilchen(j,index)%perweight
    end subroutine setPosition

    !**************************************************************************
    !****s* initPion/CoulombCorrect
    ! NAME
    ! subroutine CoulombCorrect
    ! PURPOSE
    ! Corrects the trajectory according to Coulomb forces.
    !**************************************************************************
    subroutine CoulombCorrect
      use CoulombKorrektur, only: Coulpropa

      if (DEBUG) then
         write(*,*) ' Before Coulomb correction of trajectory'
         write(*,*) 'position=', teilchen(j,index)%position
         write(*,*) '4-momentum=', teilchen(j,index)%momentum
      end if

      teilchen(j,index)%position(3)=-coulombdistance

      call Coulpropa (teilchen(j,index)%position(1:3), teilchen(j,index)%momentum(1:3), &
                      teilchen(j,index)%charge, teilchen(j,index)%mass, &
                      targetNuc%charge, distance)

      !Assume vacuum dispersion relation:
      teilchen(j,index)%velocity(1:3)=teilchen(j,index)%momentum(1:3)/FreeEnergy(teilchen(j,index))

      if (DEBUG) then
         write(*,*) ' After Coulomb correction of trajectory'
         write(*,*) 'position=', teilchen(j,index)%position
         write(*,*) '4-momentum=', teilchen(j,index)%momentum
         write(*,*)
      end if
    end subroutine coulombCorrect

    !**************************************************************************

    subroutine momentumCorrect(successFlag)
      use coulomb, only: emfoca

      real :: cPot
      logical, intent(out) :: successFlag


      if (debug) then
         write(*,*) 'Correct for in medium potentials'
         write(*,*) 'Vacuum kinetic energy:', ekin_lab
      end if

      ! Evaluate coulomb potential and correct for coulomb potential:
      cpot=0.
      if (UseCoulomb) then
         cpot = emfoca(teilchen(j,index)%position,teilchen(j,index)%momentum(1:3),teilchen(j,index)%charge,teilchen(j,index)%ID)

         if (kineticEnergy(teilchen(j,index))-cPot.lt.0) then
            write(*,*) "Error in Initpion"
            write(*,*) "Energy too small:",teilchen(j,index)%momentum
            write(*,*) "Cannot initiliaze pions at",teilchen(j,index)%position
            stop
         end if

         !Correct momentum  p**2+m**2=(E-V_coulomb)**2 :
         teilchen(j,index)%momentum(1:3)=(/0.,0.,SQRT((ekin_lab+teilchen(j,index)%mass&
              &                         -cPot)**2-teilchen(j,index)%mass**2)/)
         teilchen(j,index)%momentum(0)=ekin_lab+teilchen(j,index)%mass-cpot

         if (debug) write(*,*) 'Coulomb potential:',cpot
      end if

      !Correct for hadronic potential
      call RechneImpuls(teilchen(j,index),ekin_Lab+teilchen(j,index)%mass-cpot, successFlag)
      if (debug) then
         write(*,*) "In medium kinetic energy:",kineticEnergy(teilchen(j,index))
         write(*,'(A,4F9.6)') "In medium momentum:",teilchen(j,index)%momentum
         write(*,*)
      end if

    end subroutine momentumCorrect

    !!*******************************************************

    subroutine RechneImpuls(partIn,energy,success)
      ! Newton-Routine um Impuls des Pions zu bestimmen
      ! Lse Gleichung Wurzel(m**2+p**2)+V(p_in_LRF)=energy
      use particleDefinition
      use potentialModule
      use energyCalc

      real, intent(in) :: energy
      type(particle),intent(inOut) :: partIn
      logical, intent(out) :: success
      !local
      integer, parameter :: maxSteps=100
      type(particle) :: teilchen
      real, dimension(-1:1) ::f
      real, parameter :: dp=0.01
      real :: grad
      integer :: i,j


      teilchen=partIn

      do j=0,maxSteps
         ! Evaluate derivative d(Energy_of_TeilchenIn(p)-Energy)/dp
         do i=-1,1
            teilchen%momentum(1:2)=0.
            teilchen%momentum(3)=partIn%momentum(3)+i*dp
            call energyDetermination(teilchen)
            f(i)=teilchen%momentum(0)-energy
         end do
         if (abs(f(0)).lt.0.003) then
            if (debug) then
               write(*,*) 'Kinetic Energy in medium:',kineticEnergy(partIn)
               write(*,*) 'Kinetic Energy in vacuum:',ekin_lab
               write(*,*) 'Position=',partIn%position
            end if
            success=.true.
            return
         end if
         if (debug) write(*,*) 'f=',f
         grad=(f(1)-f(-1))/2./dp
         if (abs(grad).lt.0.0001) then
            write(*,*) "Gradient zero in RechneImpuls von initPion", grad
            write(*,*) "Energy", energy
            write(*,*) "Momentum", teilchen%momentum
            write(*,*) "Step ", j
            write(*,*) f(-1),f(0),f(1)
            success=.false.
            return
         end if
         partIn%momentum(1:2)=0.
         partIn%momentum(3)=partIn%momentum(3)-f(0)/grad
         call energyDetermination(teilchen)
      end do
      write(*,*) "Fehler in initpion.f90,RechneImpuls",f(0)
      success=.false.
    end subroutine RechneImpuls





  end subroutine initpionInduced


end module initPion
