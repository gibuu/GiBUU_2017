!******************************************************************************
!****m* /collisionCriteria
! NAME
! module collisionCriteria
!
! PURPOSE
! Includes routines to evaluate whether a 2-body collision
! among two particles is taking place.
!
! NOTES
! Available are two different collision criteria:
! * Local collision criterion: subroutine "localCollisionCriteria"
! * Kodama or "geometric collision criterion":
!   subroutines "kodama_position" and "kodama_time"
!******************************************************************************
module collisionCriteria
  implicit none
  private

  !****************************************************************************
  !****g*  collisionCriteria/debug
  ! SOURCE
  !
  logical, parameter :: debug = .false.
  ! PURPOSE
  ! Set to .true., this logical will cause debug information.
  !****************************************************************************


  !****************************************************************************
  !****g*  collisionCriteria/kodama_evalFrame
  ! SOURCE
  !
  logical, save :: kodama_evalFrame = .false.
  ! PURPOSE
  ! Set to .true., this logical will cause the kodama criterion to be
  ! evaluated in the laboratory/evaluation frame, not CM frame.
  !****************************************************************************


  logical,save :: initFlag = .true.

  public :: kodama_time,kodama_position,localCollisionCriteria

contains

  !****************************************************************************
  !****s* collisionCriteria/readInput_coll
  ! NAME
  ! subroutine readInput_coll
  !
  ! PURPOSE
  ! Reads input out of jobcard "collCriteria".
  !
  ! INPUTS
  ! none
  !
  ! OUTPUT
  ! Modifies global module variables.
  !****************************************************************************
  subroutine readInput_coll
    use output, only: Write_ReadingInput

    integer :: ios

    !**************************************************************************
    !****n*  collisionCriteria/collCriteria
    ! NAME
    ! NAMELIST collCriteria
    ! PURPOSE
    ! Includes the input switches:
    ! * kodama_evalFrame
    !**************************************************************************
    NAMELIST /collCriteria/ kodama_evalFrame

    initFlag=.false.
    call Write_ReadingInput("collCriteria",0)
    rewind(5)
    read(5,nml=collCriteria,IOSTAT=ios)
    call Write_ReadingInput("collCriteria",0,ios)

    if (kodama_evalFrame) write(*,*) 'WARNING : Kodama criterion is evaluated in the evaluation frame'

    call Write_ReadingInput("collCriteria",1)

  end subroutine readInput_coll


  !****************************************************************************
  !****f* collisionCriteria/localCollisionCriteria
  ! NAME
  ! logical function localCollisionCriteria(pair,sigmaTot,weightLocal,numEnsembles,deltaT)
  !
  ! INPUTS
  ! * type(particle),intent(in),dimension(1:2) :: pair -- Particles which should be checked for a possibility to collide
  ! * real, intent(in) :: delta_T --  Time step size
  ! * real, intent(in) :: sigmaTot     -- total cross section in mb
  ! * integer, intent(in)  :: weightLocal  -- rescaling factor for probability: p=p*weightLocal
  ! * integer, intent(in)  :: numEnsembles -- number of Ensembles
  !
  ! OUTPUT
  ! * The function is true if the criterion is fullfilled, and false if not!
  !
  ! NOTES
  ! This is based upon the prescription of Lang et al.,
  ! i.e. the "local collision criterion".
  !****************************************************************************
  logical function localCollisionCriteria(pair,sigmaTot,weightLocal,numEnsembles,deltaT)
    use particleDefinition
    use VolumeElements, only: VolumeElements_boxSize
    use random, only: rn
    use output, only: DoPr
    use lorentzTrafo, only: eval_sigmaBoost

    type(particle),intent(in),dimension(1:2) :: pair         ! in calculation frame
    real         , intent(in)                :: sigmaTot     ! in mb
    integer      , intent(in)                :: weightLocal  ! rescaling factor for probability: p=p*weightLocal
    integer      , intent(in)                :: numEnsembles ! number of Ensembles
    real         , intent(in)                :: deltaT       ! time step size [fm]

    real :: probability, dThreex, vrel, sigmaBoost
    real, dimension(1:3) :: vrel_vector ! relative velocity

    if (initFlag) call readInput_coll

    ! 1.) Get size of Volume boxes
    dThreeX=VolumeElements_boxSize()

    ! 2.) Evaluate relative velocity
    vrel_vector=pair(1)%velocity-pair(2)%velocity
    vrel=sqrt(Dot_product(vrel_vector,vrel_vector))

    ! 2a) Boost sigma to frame where both A and B are moving
    sigmaBOOST=eval_sigmaBoost(pair(1)%momentum,pair(2)%momentum)

    ! 3.) Evaluate collision probability (see e.g. Phd-thesis of Lang)
    if (dThreeX>1E-9) then
       probability=min(350.,sigmaTot*sigmaBOOST)/10./float(numensembles) * vRel* deltaT /dThreeX *float(weightLocal)
    else
       write(*,*) 'get_boxSize()=0. Critical Error! Stop!'
       stop
    end if

    ! 4.) Check collision probability and make Monte-Carlo-Decision
    if (probability>1) then
       if (DoPr(2)) then
          write(*,*) '****'
          write(*,*) 'Error in localCollisionCriteria: probability>1'
          write(*,*) 'pair = ',pair(1)%ID,pair(2)%ID
          write(*,*) 'p=',probability," sigma=",sigmaTot," srts=",sqrts(pair)
          write(*,*) "vrel= ", vrel, " deltaT=",deltaT," d^3x=",dThreeX," weight=", weightLocal
          write(*,*) '****'
       end if
       localCollisionCriteria=.true.
    else if (probability>rn()) then
       localCollisionCriteria=.true.
    else
       localCollisionCriteria=.false.
    end if

    ! 5.) Protocol the probability for debuging purposes
    if (debug) call protocol_prob

  contains

      !************************************************************************
      !****f* localCollisionCriteria/protocol_prob
      ! NAME
      ! subroutine protocol_prob
      !
      ! PURPOSE
      ! This subroutine is called if "debug=.true.". It stores all calculated
      ! probabilities in a histogram. After 1000 savings, it is writing the
      ! current histogram to file "CollisionProb.dat".
      !
      ! OUTPUT
      ! * file "CollisionProb.dat"
      !************************************************************************
      subroutine protocol_prob
        use histf90

        integer,save :: counter=0
        type(histogram) ,save :: p

        if (counter==0) call createHist(p,"Collision probability",0.,2.,0.05)

        call AddHist(p,probability,1.)
        counter=counter+1

        if (mod(counter,1000)==0) then
           open(133,File='CollisionProb.dat')
           call writeHist(p,133)
           close(133)
        end if

      end subroutine protocol_prob

  end function localCollisionCriteria


  !****************************************************************************
  !****f* collisionCriteria/kodama_position
  ! NAME
  ! logical function kodama_position(teilchen,bmax)
  !
  ! INPUTS
  ! Particles which should be checked for a possibility to collide:
  ! * type(particle),intent(in),dimension(1:2) :: teilchen
  !
  ! Maximal impact parameter:
  ! * real, intent(in) :: bMax
  !
  ! NOTES
  ! The particles' velocities and positions need to be properly defined!
  ! The notation is according to Effenberger's Dr. thesis  page 251. This
  ! implements the criterion that the relative distance in space of two
  ! particles should be smaller than bMax in the CM-Frame of both particles.
  !
  ! OUTPUT
  ! The function is true if the criterion is fullfilled, and false if not!
  !****************************************************************************
  logical function kodama_position(teilchen,bMax)
    use lorentzTrafo, only: lorentz
    use particleDefinition

    type(particle),intent(in),dimension(1:2) :: teilchen
    real,          intent(in)                :: bMax

    real, dimension(0:3,1:2) :: position,velocity ! first Index: fourVector, second: particle
    real, dimension(1:3) :: beta_12, x_12         ! according to Effe's notation, relative position and velocity
    real, dimension(1:3) :: beta                  ! velocity of CM-frame
    real :: bMin                                  ! Minimal distance of both particles

    if (initFlag) call readInput_coll

    ! Set 4-positions of both particles:
    ! Both particles at equal times in lab frame. Assume t=0:
    position(0,1:2)=0.

    position(1:3,1)=teilchen(1)%position(1:3)
    position(1:3,2)=teilchen(2)%position(1:3)

    ! Define 4-velocities of both particles
    ! First particle:
    velocity(0,1)=1./Sqrt(1.-Dot_Product(teilchen(1)%velocity(1:3),teilchen(1)%velocity(1:3)))
    velocity(1:3,1)=teilchen(1)%velocity(1:3)*velocity(0,1)
    ! Second particle:
    velocity(0,2)=1./Sqrt(1.-Dot_Product(teilchen(2)%velocity(1:3),teilchen(2)%velocity(1:3)))
    velocity(1:3,2)=teilchen(2)%velocity(1:3)*velocity(0,2)

    ! Evaluate the velocity of the CM-frame of both particles:
    beta(1:3)=(teilchen(1)%momentum(1:3)+teilchen(2)%momentum(1:3))/ &
              (teilchen(1)%momentum(0)+teilchen(2)%momentum(0))

    if (.not.kodama_evalFrame) then
       ! Boost everything to the CM-frame
       call lorentz (beta, position(0:3,1), 'kodama_position(1)')
       call lorentz (beta, position(0:3,2), 'kodama_position(2)')
       call lorentz (beta, velocity(0:3,1), 'kodama_position(3)')
       call lorentz (beta, velocity(0:3,2), 'kodama_position(4)')
    end if

    ! Evaluate the minimal distance bMin, assuming straight trajectories
    x_12(1:3)=position(1:3,1)-position(1:3,2) &
         -velocity(1:3,1)/velocity(0,1)*position(0,1)&
         +velocity(1:3,2)/velocity(0,2)*position(0,2)

    beta_12(1:3)=velocity(1:3,1)/velocity(0,1)-velocity(1:3,2)/velocity(0,2)

    if (Dot_Product(beta_12,beta_12) > 1E-9) then
       bMin=Dot_Product(x_12,x_12)-(Dot_Product(x_12,beta_12))**2/Dot_Product(beta_12,beta_12)
       if (bmin>0.) then
          bmin=sqrt(bmin)
       else if (bmin>-1E-4) then
          bmin=0.
       else
          write(*,*) 'Problem with bmin in kodama_position!'
          write(*,*) 'Warning: There are particles with bmin^2<0:'
          write(*,*) 'x_12=',x_12
          write(*,*) 'beta_12=',beta_12
          write(*,*) 'bmin=',bmin
          write(*,*) 'particle 1:',teilchen(1)%ID,'pos=',teilchen(1)%position, 'mom=',teilchen(1)%momentum
          write(*,*) 'particle 2:',teilchen(2)%ID,'pos=',teilchen(2)%position, 'mom=',teilchen(2)%momentum
          bmin=0 ! Set bmin to 0, otherwise the result is invalid since complex
       end if
    else
       bMin=sqrt(Dot_Product(x_12,x_12))
       write( *,*) 'Warning: There are particles which are in relative rest'
       write(*,*) teilchen(1)%ID, teilchen(1)%position, &
            & '###', teilchen(1)%momentum,'###',velocity(1:3,1), &
            & '###',teilchen(1)%lastcollisionTime, teilchen(1)%number,teilchen(1)%event
       write(*,*) teilchen(2)%ID, teilchen(2)%position, &
            & '###', teilchen(2)%momentum,'###',velocity(1:3,2), &
            & '###',teilchen(2)%lastcollisionTime, teilchen(2)%number,teilchen(2)%event
    end if

    kodama_position = (bMin<bmax)
    !    write(*,*) 'In kodama_position',bmin,bmax

  end function kodama_position



  !****************************************************************************
  !****f* collisionCriteria/kodama_time
  ! NAME
  ! logical function kodama_time(teilchen,delta_T,collision_time)
  !
  ! INPUTS
  ! Particles which should be checked for a possibility to collide:
  ! * type(particle),intent(in),dimension(1:2) :: teilchen
  !
  ! Time step size:
  ! * real, intent(in) :: delta_T
  !
  ! OUTPUT
  ! * Time instant in the computational frame when the two particles collide
  !   (with respect to the current "BUU time" used for time stepping):
  !   real, optional, intent(out) :: collision_time
  ! * The function is true if the criterion for tau_1 and tau_2 is fullfilled,
  !   and false if not!
  !
  ! NOTES
  ! The notation is according to Effenberger's Dr. thesis  pages 251-252.
  ! This function implements the criterion that the relative distance in time
  ! of two particles should fulfill Abs(tau_1+tau_2)<delta_T. This was first
  ! implemented by G.Wolf. The meanings of tau_1 and tau_2 are explained in
  ! Effenbergers thesis pages 251ff.
  !****************************************************************************
  logical function kodama_time(teilchen,delta_T,collision_time)

    use lorentzTrafo, only: lorentz
    use particleDefinition
    use output, only: WriteParticle, WriteParticle_debug
    use callstack, only: traceBack

    type(particle),intent(in),dimension(1:2) :: teilchen
    real, intent(in) :: delta_T
    real, optional, intent(out) :: collision_time

    real,dimension(1:3) :: beta, beta_ij, x_ij, beta_12, x_12
    real,dimension(0:3) :: velocity_j, position_i, position_j
    integer :: i,j
    real :: gamma, t_min, gamma_j, beta_12_squared
    real, dimension(1:2) :: tau

    if (initFlag) call readInput_coll

    ! Evaluate the minimal distance bMin, assuming straight trajectories:

    if (kodama_evalFrame) then
       ! ** Decision is performed in the evaluation frame
       x_12(1:3)=teilchen(1)%position(1:3)-teilchen(2)%position(1:3)
       beta_12(1:3)=teilchen(1)%velocity(1:3)-teilchen(2)%velocity(1:3)

       beta_12_squared=Dot_product(beta_12,beta_12)
       if (beta_12_squared>0.) then
          ! See Effenberger Dissertation, page 251
          ! Using the d^2(t) of formula B2 and demanding for  d(d^2(t))/dt=0  we get:
          t_min=-Dot_product(x_12,beta_12)/Dot_product(beta_12,beta_12)
       else
          write(*,*) 'Warning in Kodama Time criterion: beta_12_squared.le.0', beta_12_squared
          write(*,*) beta_12
          t_min=0.
       end if

       kodama_time = (abs(t_min)<=Delta_T/2.)

       if (present(collision_time)) collision_time=t_min

    else
       ! ** Decision is performed in the CM Frame (Wolff et al.) :
       ! We evaluate for each particle in its own rest frame the time tau, which it takes
       ! until it has reached closest distance to the other particle.
       ! This time tau might as well be negative since this point in time might lie in the past
       ! in this special rest frame.

       ! Loop over both particles' rest-frames:

       if (debug) write(*,*) 'velos 1=', teilchen(1)%velocity
       if (debug) write(*,*) 'velos 2=', teilchen(2)%velocity


       do i=1,2 !i is the index of the particle whose restframe we are using to evaluate tau
          ! j is index of second particle
          j=mod(i,2)+1
          if (debug) write(*,*) "i,j:" ,i,j
          !In lab frame at equal times, assume t=0 :
          position_j(0)=0.
          position_i(0)=0.

          !Set positions in lab frame:
          position_j(1:3)=teilchen(j)%position(1:3)
          position_i(1:3)=teilchen(i)%position(1:3)

          !Boost Parameters for boost to restframe of i-th particle:
          beta(1:3)=teilchen(i)%velocity(1:3)
          gamma = 1-Dot_Product(beta,beta)
          if (gamma<=0.) then
             write(*,*) 'kodama_time: wrong boost vector!'
             call WriteParticle(6,99,i, teilchen(i))
             stop
          end if
          gamma = 1./sqrt(gamma)

          if (debug) write(*,*) beta
          !Velocity of particle j:
          gamma_j = 1. - Dot_Product(teilchen(j)%velocity(1:3),teilchen(j)%velocity(1:3))
          if (gamma_j<=0.) then
             write(*,*) 'Severe ERROR!!!! You must debug it!!'
             write(*,*) 'kodama_time: wrong velocity_j vector!'
             write(*,*) 'velocity=',teilchen(j)%velocity
             write(*,*) 'v^2=',Dot_Product(teilchen(j)%velocity(1:3),teilchen(j)%velocity(1:3))
             write(*,*) 'gamma_j=', gamma_j
             call WriteParticle(6,99,j, teilchen(j))
             call traceBack('problem in kodama_time: bad velocity_j!',-1)
             stop
          end if

          velocity_j(0)=1./sqrt(gamma_j)
          velocity_j(1:3)=teilchen(j)%velocity(1:3)*velocity_j(0)

          if (debug) write(*,*) 'velo=', velocity_j

          ! Boost everything to restframe of particle i
          call lorentz (beta, position_j, 'kodama_time(1)')
          call lorentz (beta, position_i, 'kodama_time(2)')
          call lorentz (beta, velocity_j, 'kodama_time(3)')

          ! Since velocity of particle i is zero in its own restframe we get:
          x_ij(1:3)=position_i(1:3)-position_j(1:3)+velocity_j(1:3)/velocity_j(0)*position_j(0)
          beta_ij(1:3)=-velocity_j(1:3)/velocity_j(0)

          if (debug) then
             write(*,*) 'beta=',beta_ij,'x_ij=', x_ij
             call WriteParticle_debug(teilchen(1))
             call WriteParticle_debug(teilchen(2))
             write(*,*) -Dot_Product(x_ij,beta_ij),Dot_Product(beta_ij,beta_ij)
          end if
          t_min=-Dot_Product(x_ij,beta_ij)/Dot_Product(beta_ij,beta_ij)
          tau(i)=(t_min-position_i(0))*gamma
       end do

       ! The collision criterion: Here we decide wether the collision can
       ! happen in the timestep of size delta_T:
       kodama_time = (abs(tau(1)+tau(2))<=Delta_T)

       if (present(collision_time)) collision_time=(tau(1)+tau(2))/2.
    end if

  end function kodama_time


end module collisionCriteria
