!******************************************************************************
!****m* /coulombKorrektur
! NAME
! module CoulombKorrektur
! PURPOSE
! Includes routine for Coulomb correction of particle Trajectories.
!******************************************************************************
module CoulombKorrektur

  implicit none
  private

  public :: coulombPropagation, Coulpropa, CoulpropaTwo

contains


  !****************************************************************************
  !****s* coulombKorrektur/coulombPropagatioon
  ! NAME
  ! subroutine coulombPropagation(teilchen)
  ! PURPOSE
  ! This routine does a simple coulomb propagation for a particle vector in the field of the target nucleus.
  ! USES
  ! CoulPropa
  !****************************************************************************
  subroutine coulombPropagation (teilchen)
    use particleDefinition
    use nucleusDefinition
    use nucleus, only: getTarget
    use inputGeneral, only: PrintParticleVectors
    use output, only: writeParticleVector

    type(particle), intent(inOUT), dimension(:,:)  ::  teilchen
    integer i, j
    real, dimension(1:3) :: mom
    type(tNucleus),pointer :: targetNucX
    real, parameter :: propagationTime = 200.    ! duration of the Coulomb propagation in units of fm/c.

    if (PrintParticleVectors) then
       ! Return the Particles before the Coulombcorrection
       call writeParticleVector("PertParticles_BeforeCoulomb", teilchen)
    end if

    ! Do the propagation
    targetNucX => getTarget()

    do i=lbound(teilchen, dim=1),ubound(teilchen, dim=1)
       do j=lbound(teilchen, dim=2),ubound(teilchen, dim=2)
          mom=teilchen(i,j)%momentum(1:3)
          call Coulpropa(teilchen(i,j)%position,mom,teilchen(i,j)%charge,teilchen(i,j)%mass,targetNucX%charge,3.,propagationTime)
          teilchen(i,j)%momentum(1:3)=mom
          teilchen(i,j)%momentum(0)=SQRT(teilchen(i,j)%mass**2+Dot_Product(mom,mom))
       end do
    end do
  end subroutine coulombPropagation



  !****************************************************************************
  !****s* coulombKorrektur/coulPropa
  ! NAME
  ! subroutine Coulpropa(r,p,charge,mass,externalCharge,distance,propagationTime)
  ! PURPOSE
  ! This routine does a simple coulomb propagation for a particle with charge "charge".
  ! A point charge of strength "externalCharge" is assumed to sit in the origin.
  ! The propagation is terminated when
  ! ... "propagtionTime" is not present as input given and the particle is closer
  ! than "distance" to the origin or when it's starting to move away from the origin. This mode is
  ! used to correct the trajectories of incoming particles.
  ! ... "propagationTime" is present as input and the propagationTime is reached. This mode is used
  ! to correct the trajectories of outgoing particles.
  ! Relativistic kinematics!
  ! INPUTS
  ! real mass   : mass of particle
  ! integer charge : charge of particle
  ! integer ExternalCharge  : fixed charge at origin
  ! real distance  : minimal distance to which the particle is propagated
  !             towards the origin where the point charge is assumed to sit
  ! real propagationTime : optional argument
  ! real (1:3) r,p  : position and momentum of particle at start
  ! OUTPUT
  ! real (1:3) r,p : position and momentum of particle after propagation
  !****************************************************************************
  subroutine Coulpropa (r, p, charge, mass, externalCharge, distance, propagationTime, printout_in)
    use constants, only: alphaQED, hbarc
    real,dimension(3), intent(inout) :: r,p !initial position and momentum of propagated particle
    real, intent(in)    :: mass             !mass of particle
    integer, intent(in) :: charge           !charge of particle
    integer, intent(in) :: ExternalCharge   !fixed charge at origin
    real, intent(in)    ::distance          !minimal distance to which the particle is propagated towards the origin where the point charge is assumed to sit
    real, intent(in),optional :: propagationTime     !total time to propagate the particle
    logical, optional :: printout_in
    ! Propagation :
    integer :: nt=100000000             ! maximal number of time steps
    real, parameter    :: dt=0.2        ! time step size in fm/c
    real,dimension(3) :: grad           ! Gradient of Coulomb potential
    real :: E                           ! Energy
    !Predictor-Corrector method:
    logical, parameter :: predCor=.true.   ! Switch on predictor/corrector method
    real,dimension(3)  :: rp,pp            ! predicted position and momentum
    real,dimension(3)  :: gradPredictor    ! predicted gradient of Coulomb Potential
    real               :: EPredictor       ! predicted Energy
    logical            :: printout
    !Printout :
    integer, parameter :: printSteps=10  !printSteps*dt is difference in time for the printout
    integer :: i
    real :: rsave                         !variable used to save the distance to the origin for the next time step

    if (present(printout_in)) then
       printout=printout_in
    else
       printout=.false.
    end if

    if (present(propagationTime)) then
       nt=NINT(propagationTime/dt)
    end if

    if (printout) then
       open(99,File='CoulombPropagation.dat',position='append')
       open(100,File='CoulombPropagation_EnergyCheck.dat',position='append')
    end if
    rSave=2.*Dot_Product(r,r)
    do i=1, nt
       E=SQRT(Dot_Product(p,p)+mass**2)
       if (Sqrt(Dot_Product(r,r)).lt.0.0001) return
       if (PredCor) then
          grad=-r(:)/(Sqrt(Dot_Product(r,r)))**3*externalCharge*alphaQED*hbarc
          !predictor step
          pp=p-dt*grad*charge
          rp=r+p/E*dt

          if (Sqrt(Dot_Product(rp,rp)).lt.0.0001) return

          gradPredictor=-rp(:)/(Sqrt(Dot_Product(rp,rp)))**3*externalCharge*alphaQED*hbarc
          !corrector step
          EPredictor=SQRT(Dot_Product(pp,pp)+mass**2)
          r=r+dt/2.*(p/E+pp/EPredictor)
          p=p-dt/2.*(grad+gradPredictor)*charge
       else
          if (Sqrt(Dot_Product(r,r)).lt.0.0001) return
          grad=-r(:)/(Sqrt(Dot_Product(r,r)))**3*externalCharge*alphaQED*hbarc
          r=r+p/E*dt
          p=p-dt*grad*charge
       end if

       ! The propagation is terminated if the particle is closer
       ! than "distance" to
       ! the origin or when it's starting to move away from the origin.

       if (.not.present(PropagationTime)) then
          if (Sqrt(Dot_product(r,r)).lt.distance &
               & .or.(Dot_Product(r,r).gt.rSave)) exit
          rSave=Dot_Product(r,r)
       end if

       if (printOut) then !Print trajectories
          if (Mod(i,printSteps).eq.0.and.(sum(abs(r)).lt.40)) then
             !Print position
             write(99,*) r
             !Print full energy
             write(100,*) i*dt, sqrt(sum(p*p)+mass**2) &
                  & +charge*externalCharge/137./sqrt(sum(r*r))*0.2
          end if

       end if

    end do
    if (printout) then
       close(99)
       close(100)
    end if
  end subroutine Coulpropa



  !****************************************************************************
  !****s* coulombKorrektur/coulpropaTwo
  ! NAME
  ! subroutine CoulpropaTwo(r1,p1,charge1,mass1,r2,p2,charge2,mass2,distance,propagationTime)
  ! PURPOSE
  ! This routine does a simple propagation for two point charges which interact
  ! with each other via the coulomb force. Relativistic kinematics!
  ! The propagation is terminated when
  ! ... the particles are closer than "distance" to each other or when they are starting
  ! to move away from each other.
  ! ... or when the propagationTime is over.
  ! INPUTS
  ! real mass1,mass2   : masses of particle
  ! integer charge1,charge2 : charges of particle
  ! real distance  : minimal distance to which the particle is propagated
  !             towards the origin where the point charge is assumed to sit
  ! real propagationTime : optional argument
  ! real (1:3) r1,p1,r2,p2  : position and momentum of particles before propagation
  ! OUTPUT
  ! real (1:3) r1,p1,r2,p2  : position and momentum of particles after propagation
  !****************************************************************************
  subroutine CoulpropaTwo (r1, p1, charge1, mass1, r2, p2, charge2, mass2, distance)!,propagationTime)
    use constants, only: alphaQED, hbarc
    real, intent(in)                  :: distance          ! minimal distance of approach of both point charges
    real, dimension(3), intent(inout) :: r1,p1             ! initial position and momentum of particle 1
    real, dimension(3), intent(inout) :: r2,p2             ! initial position and momentum of particle 2
    real, intent(in)                  :: mass1,mass2       ! masses of charges  1&2
    integer, intent(in)               :: charge1,charge2   ! charges of charges 1&2
    !real, intent(in), optional        :: propagationTime   ! total time to propagate the particle
    !propagation
    integer, parameter :: nt=100000000 ! maximal number of time steps
    real, parameter    :: dt=0.05      ! in fm/c
    real,dimension(3) :: grad1,grad2  ! Gradient of Coulomb potential for charge 1 and 2
    real :: E1,E2                     ! Energies of particle 1&2
    !predictor-corrector-method
    logical, parameter :: predCor=.false.              ! Switch on predictor/corrector method
    real,dimension(3) :: pp1,rp1,pp2,rp2               ! predicted momenta and positions
    real,dimension(3) :: gradPredictor1,gradPredictor2 ! predicted gradients of Coulomb potential
    real :: EPredictor1,EPredictor2                    ! predicted energies
    !printout
    logical, parameter :: printout=.false. ! switch to turn off printout
    integer, parameter :: steps=10         ! steps*dt defines when to print out the particles' positions
    integer :: i
    real :: rsave

    if (printout) then
       open(11,File='particle_1_trajectory.dat',position='append')
       open(22,File='particle_2_trajectory.dat',position='append')
       open(222,File='particles_energyCheck.dat',position='append')
       open(111,File='particles_momentumCheck.dat',position='append')
    end if
    rSave=Dot_Product(r1-r2,r1-r2)*2.

!     If (Present(propagationTime)) then
!        nt=NINT(propagationTime/dt)
!     end if

    do i=1, nt
       if (PredCor) then
          !predictor step

          grad1=-(r1-r2)/Sqrt(dot_Product(r1-r2,r1-r2))**3*charge1*charge2*alphaQED*hbarc
          grad2=-(r2-r1)/Sqrt(dot_Product(r1-r2,r1-r2))**3*charge1*charge2*alphaQED*hbarc

          !Predict values for charge 1:
          E1=SQRT(Dot_Product(p1,p1)+mass1**2)
          rp1(:)=r1(:)+p1(:)/E1*dt
          pp1(:)=p1(:)-dt*grad1(:)
          EPredictor1=SQRT(Dot_Product(pp1,pp1)+mass1**2)

          !Predict values for charge 2:
          E2=SQRT(Dot_Product(p2,p2)+mass2**2)
          rp2(:)=r2(:)+p2(:)/E2*dt
          pp2(:)=p2(:)-dt*grad2(:)
          EPredictor2=SQRT(Dot_Product(pp2,pp2)+mass2**2)

          !Evaluate gradients of potential at predicted positions:
          gradPredictor1=-(rp1-rp2)/Sqrt(dot_Product(rp1-rp2,rp1-rp2))**3 &
               &          *charge1*charge2*alphaQED*hbarc
          gradPredictor2=-(rp2-rp1)/Sqrt(dot_Product(rp1-rp2,rp1-rp2))**3 &
               &          *charge1*charge2*alphaQED*hbarc

          !Corrector step:
          r1(:)=r1(:)+dt/2.*(p1(:)/E1+pp1(:)/EPredictor1)
          p1(:)=p1(:)-dt/2.*(grad1(:)+gradPredictor1(:))

          r2(:)=r2(:)+dt/2.*(p2(:)/E2+pp2(:)/EPredictor2)
          p2(:)=p2(:)-dt/2.*(grad2(:)+gradPredictor2(:))
       else
          !dV_Coulomb/dr1
          grad1=-(r1-r2)/Sqrt(dot_Product(r1-r2,r1-r2))**3*charge1*charge2*alphaQED*hbarc
          !dV_Coulomb/dr2
          grad2=-(r2-r1)/Sqrt(dot_Product(r1-r2,r1-r2))**3*charge1*charge2*alphaQED*hbarc

          !charge 1:
          E1=SQRT(Dot_Product(p1,p1)+mass1**2)
          r1(:)=r1(:)+p1(:)/E1*dt
          p1(:)=p1(:)-dt*grad1(:)

          !charge 2:
          E2=SQRT(Dot_Product(p2,p2)+mass2**2)
          r2(:)=r2(:)+p2(:)/E2*dt
          p2(:)=p2(:)-dt*grad2(:)
       end if

       if (Sum(abs(r1-r2)).lt.1E-08) then
          write(*,*) 'The point charges are getting to close!'
          write(*,*) 'This causes numerical difficulties in Coulomb potential (due to 1/r)!!!'
          write(*,*) 'Positions:',r1,'and',r2
          stop
       end if

       if (printout) then
          if (Mod(i,steps).eq.0) then
             !               Print *, nt,steps
             !Print position of particle 1
             write(11,*) i*dt,r1

             !Print position and of particle 2
             write(22,*) i*dt,r2
             !Print full energy of particles
             write(222,*) i*dt, E1+E2  &
                  & +charge1*charge2*alphaQED/sqrt(Dot_Product((r1-r2),(r1-r2)))*hbarc
             write(111,*) i*dt,p1+p2
          end if
       end if
!        If (.not.Present(PropagationTime)) then
       if ((Sqrt(Dot_product(r1-r2,r1-r2))<distance) .or. (rSave<Dot_Product(r1-r2,r1-r2))) exit
       rSave = Dot_Product(r1-r2,r1-r2)
!        end if
    end do
    if (printOut) then
       close(11)
       close(22)
       close(222)
       close(111)
    end if

    return
  end subroutine CoulpropaTwo



!   real function CoulombPotential(r,charge)
!     use constants, only : alphaQED, hbarc
!     real,intent(in)    :: r        ! abs(r) in fm of probe
!     integer,intent(in) :: charge   ! charge of probe
!
!     logical, save :: initFlag=.true.
!     integer, parameter :: numPoints=100
!     real, dimension(0:numPoints) :: pot
!     real, save :: dr
!     real :: r_cutOff
!     real :: A, diff
!     integer :: index_r_low,index_r_up
!
!     if(initFlag) then
!        ! Tabulate potential
!        call tabulate_CoulPot(pot,r_cutOff,A,dr)
!        initFlag=.false.
!     end if
!     if(r.ge.r_cutOff) then
!        coulombPotential=1/r*A*charge*alphaQED*hbarc
!     else
!        ! linear interpolation between grid points
!        index_r_low=int(r/dr)
!        index_r_up=int(r/dr)+1
!        diff=r-(index_r_low)*dr
!        coulombPotential=((1.-diff)*pot(index_r_low)+diff*pot(index_r_up))*charge
!     end if
!   end function CoulombPotential



!   subroutine tabulate_CoulPot(pot,r_cutOff,A,dr)
!     use densityModule, only :  densityAt
!     use constants, only : pi
!     use dichteDefinition
!     use minkowski, only : abs4
!     real, intent(out),  dimension(0:) :: pot
!     real, intent(out) :: dr
!     real, intent(out) :: r_cutOff
!     real, intent(out) :: A
!     real :: r,protonDens
!     type(dichte) :: dens
!     real, parameter :: dens_cutOff=1E-5
!     integer :: index
!     ! Find useful cutoff first
!     r=0.
!     A=0.
!     find_cutOff_loop: do
!        r=r+0.01
!        dens=densityAt((/r,0.,0./))
!        protonDens=abs4(dens%proton)
!        A=A+r**2*protonDens
!        if(protonDens.lt.dens_cutOff) then
!           r_cutOff=r
!           exit find_cutOff_loop
!        end if
!        if(r.gt.1000) then
!           write(*,*) 'Error in tabulate_CoulPot!'
!           write(*,*) 'Density does not approach 0!! STOP !!'
!           stop
!        end if
!     end do find_cutOff_loop
!     A=A*4.*pi
!     dr=r_cutoff/float(size(pot,dim=1)-1)
!
!
!     do index=lbound(pot,dim=1),ubound(pot,dim=1)
!        r=float(index)*dr
!        pot(index)=eval_pot(r)
!     end do
!
!     contains
!
!       real function eval_pot(r)
!         real, intent(in) :: r
!         eval_pot=0. ! Not yet finished
!       end function eval_pot
!
!   end subroutine tabulate_CoulPot


end module CoulombKorrektur
