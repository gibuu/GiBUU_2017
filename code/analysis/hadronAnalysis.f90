!******************************************************************************
!****m* /hadronAnalysis
! NAME
! module hadronAnalysis
!
! PURPOSE
! Contains output routines for the hadron induced reactions
! NOTES
! The routines of this module do not influence the physical processes.
! For some purposes also routines from the module HeavyIonAnalysis
! can be useful.
!******************************************************************************
module hadronAnalysis
  implicit none
  private

  !****************************************************************************
  !****n* hadronAnalysis/HadronAnalysis
  ! PURPOSE
  ! Namelist for hadronAnalysis includes:
  ! * flagAnalysis
  !****************************************************************************

  !****************************************************************************
  !****g* hadronAnalysis/flagAnalysis
  ! SOURCE
  !
  logical, save :: flagAnalysis = .false.
  ! PURPOSE
  ! If true, perform the output of a hadron at the latest
  ! time step before the hadron disappeared (file DoHadronAnalysisTime.dat).
  ! The output hadron has the same baryon/meson type and antiparticle-flag
  ! as the beam particle. In case if the hadron did not disappear, the output
  ! is done at the end of the time evolution.
  ! The output for the hadron is also done in three another files if its momentum
  ! becomes for the first time less than the cut values pCut1 and pCut2
  ! (files DoHadronAnalysisTime1.dat and  DoHadronAnalysisTime2.dat)
  ! and if it becomes bound (DoHadronAnalysisTime3.dat)
  !
  ! NOTES
  ! Presently feasible only for real particle simulations.
  !****************************************************************************

  !****************************************************************************
  !****t* hadronAnalysis/hadron
  ! NAME
  ! type particle
  ! PURPOSE
  ! The type definition only for the output purposes.
  ! SOURCE
  !
  Type hadron
     integer               :: ID=0
     integer               :: charge=0
     integer               :: event=0
     real                  :: mass=0.
     real, dimension (1:3) :: position=0.
     real, dimension (1:3) :: momentum=0.
  End Type hadron
  !****************************************************************************

  ! Working fields:
  real, save, dimension(:,:), allocatable :: times
  type(hadron), save, dimension(:,:), allocatable :: hadrons


  public :: DoHadronAnalysisTime,cleanUp


contains


  !****************************************************************************
  !****s* hadronAnalysis/readInput
  ! NAME
  ! subroutine readInput
  ! PURPOSE
  ! Reads input in jobcard out of namelist "HadronAnalysis"
  !****************************************************************************
  subroutine readInput

    use output
    integer :: ios

    NAMELIST /HadronAnalysis/ flagAnalysis

    rewind(5)
    call Write_ReadingInput('HadronAnalysis',0)
    read(5,nml=HadronAnalysis,iostat=ios)
    call Write_ReadingInput('HadronAnalysis',0,ios)
    if (flagAnalysis) then
      write(*,*) ' Analysing output for hadron-induced reaction will be performed:'
    end if
    call Write_ReadingInput('HadronAnalysis',1)

  end subroutine readInput


  !****************************************************************************
  !****s* hadronAnalysis/DoHadronAnalysisTime
  ! NAME
  ! subroutine DoHadronAnalysisTime (realparticles, time, finalFlag)
  !
  ! PURPOSE
  ! Makes the output of the particle of the same type as incoming hadron
  ! (same isBaryon flag and same antiparticle flag) on the disk.
  !
  ! INPUTS
  ! * type(particle), dimension(:,:), intent(in)  :: realparticles   ! real particle vector
  ! * real, intent(in) :: time                                       ! current time step
  ! * logical, intent (in)     :: finalFlag   ! .true. if it is the last call for one specific
  !   energy, therefore final output must be made.
  ! RESULT
  ! See the ouput file DoHadronAnalysisTime.dat for details.
  !****************************************************************************
  subroutine DoHadronAnalysisTime (realparticles, time, finalFlag)

    use particleDefinition
    use initHadron, only: b,particleId,antiparticle !,particleCharge
    use IdTable, only: isBaryon
    use densitymodule, only: Particle4Momentum_RMF
    use coulomb, only: emfoca
    use RMF, only:  getRMF_flag

    type(particle), dimension(:,:), intent(in) :: realParticles
    real, intent(in) :: time
    logical, intent (in) :: finalFlag

    real, parameter :: pCut1=0.5, pCut2=0.3 ! Cut values of the hadron momentum (GeV/c)
    real, parameter :: rCut=3.              ! Cut value of the hadron radial position (fm)

    integer :: i,j
    real :: p,r,E
    real, dimension(0:3) :: momentum
    real, dimension(1:3) :: place
    integer, save :: numEnsembles
    integer, save :: isu=0     ! counter of subsequent runs
    logical, save :: flagIni=.true.

    if (flagIni) then
       call readInput  ! read-in flagAnalysis
       if (flagAnalysis) then
          numEnsembles=size(realParticles,dim=1)
          allocate(hadrons(1:numEnsembles,0:3),times(1:numEnsembles,0:3))
       end if
       flagIni=.false.
    end if

    if (.not.flagAnalysis) return

    if (finalFlag) then
       isu=isu+1
       open(30,file='DoHadronAnalysisTime.dat',position='Append')
       do i = 1,numEnsembles
         if (hadrons(i,0)%Id.ne.0) then
            write(30,5) hadrons(i,0)%Id,hadrons(i,0)%charge,&
                       &hadrons(i,0)%event,&
                       &hadrons(i,0)%mass,hadrons(i,0)%position,&
                       &hadrons(i,0)%momentum,times(i,0),i,isu,b
5           format(i4,1x,i2,1x,i8,1x,f5.3,3(1x,f8.3),3(1x,f8.3),1x,f5.1,1x,i5,1x,i4,1x,f5.2)
         end if
       end do
       open(30,file='DoHadronAnalysisTime1.dat',position='Append')
       do i = 1,numEnsembles
         if (hadrons(i,1)%Id.ne.0) then
            write(30,5) hadrons(i,1)%Id,hadrons(i,1)%charge,&
                       &hadrons(i,1)%event,&
                       &hadrons(i,1)%mass,hadrons(i,1)%position,&
                       &hadrons(i,1)%momentum,times(i,1),i,isu,b
         end if
       end do
       open(30,file='DoHadronAnalysisTime2.dat',position='Append')
       do i = 1,numEnsembles
         if (hadrons(i,2)%Id.ne.0) then
            write(30,5) hadrons(i,2)%Id,hadrons(i,2)%charge,&
                       &hadrons(i,2)%event,&
                       &hadrons(i,2)%mass,hadrons(i,2)%position,&
                       &hadrons(i,2)%momentum,times(i,2),i,isu,b
         end if
       end do
       if (getRMF_flag()) then
          open(30,file='DoHadronAnalysisTime3.dat',position='Append')
          do i = 1,numEnsembles
            if (hadrons(i,3)%Id.ne.0) then
               write(30,5) hadrons(i,3)%Id,hadrons(i,3)%charge,&
                          &hadrons(i,3)%event,&
                          &hadrons(i,3)%mass,hadrons(i,3)%position,&
                          &hadrons(i,3)%momentum,times(i,3),i,isu,b
            end if
          end do
       end if
       hadrons(:,:)%Id=0
       return
    end if


    do i = 1,numEnsembles

       ParticleLoop1 : do j = 1,size(realParticles,dim=2)
         if (realParticles(i,j)%Id == 0) then
             cycle ParticleLoop1
         else if (realParticles(i,j)%Id < 0) then
             exit ParticleLoop1
         end if
!          if(        realParticles(i,j)%Id.eq.particleId &
!              &.and. (realParticles(i,j)%antiparticle.eqv.antiparticle) &
!              &.and. realParticles(i,j)%charge.eq.particleCharge ) then
          if (        (isBaryon(realParticles(i,j)%Id).eqv.isBaryon(particleId)) &
              &.and. (realParticles(i,j)%antiparticle.eqv.antiparticle) ) then
            !**** Select the particle if it exists:
            hadrons(i,0)%Id=realParticles(i,j)%Id
            if (realParticles(i,j)%antiparticle) hadrons(i,0)%Id=-hadrons(i,0)%Id
            hadrons(i,0)%charge=realParticles(i,j)%charge
            hadrons(i,0)%event=realParticles(i,j)%event(1)
            hadrons(i,0)%mass=realParticles(i,j)%mass
            hadrons(i,0)%position=realParticles(i,j)%position
            hadrons(i,0)%momentum(1:3)=realParticles(i,j)%momentum(1:3)
            times(i,0)=time
            p=sqrt(dot_product(realParticles(i,j)%momentum(1:3),&
                              &realParticles(i,j)%momentum(1:3)))
            r=sqrt(dot_product(realParticles(i,j)%position(1:3),&
                              &realParticles(i,j)%position(1:3)))
            !**** Select the particle if its momentum gets less than pCut1:
            if (p.lt.pCut1 .and. r.lt.rCut .and. hadrons(i,1)%Id.eq.0) then
                hadrons(i,1)%Id=realParticles(i,j)%Id
                if (realParticles(i,j)%antiparticle) hadrons(i,1)%Id=-hadrons(i,1)%Id
                hadrons(i,1)%charge=realParticles(i,j)%charge
                hadrons(i,1)%event=realParticles(i,j)%event(1)
                hadrons(i,1)%mass=realParticles(i,j)%mass
                hadrons(i,1)%position=realParticles(i,j)%position
                hadrons(i,1)%momentum(1:3)=realParticles(i,j)%momentum(1:3)
                times(i,1)=time
            end if
            !**** Select the particle if its momentum gets less than pCut2:
            if (p.lt.pCut2 .and. r.lt.rCut .and. hadrons(i,2)%Id.eq.0) then
                hadrons(i,2)%Id=realParticles(i,j)%Id
                if (realParticles(i,j)%antiparticle) hadrons(i,2)%Id=-hadrons(i,2)%Id
                hadrons(i,2)%charge=realParticles(i,j)%charge
                hadrons(i,2)%event=realParticles(i,j)%event(1)
                hadrons(i,2)%mass=realParticles(i,j)%mass
                hadrons(i,2)%position=realParticles(i,j)%position
                hadrons(i,2)%momentum(1:3)=realParticles(i,j)%momentum(1:3)
                times(i,2)=time
            end if
            !**** Select the particle if its binding energy becomes positive:
            if (getRMF_flag() .and. hadrons(i,3)%Id.eq.0) then
               call Particle4Momentum_RMF(realParticles(i,j),momentum)
               place=realParticles(i,j)%position
               E = momentum(0) + emfoca(place,(/0.,0.,0./),realParticles(i,j)%charge,realParticles(i,j)%ID)
               if (E.lt.realParticles(i,j)%mass) then
                   hadrons(i,3)%Id=realParticles(i,j)%Id
                   if (realParticles(i,j)%antiparticle) hadrons(i,3)%Id=-hadrons(i,3)%Id
                   hadrons(i,3)%charge=realParticles(i,j)%charge
                   hadrons(i,3)%event=realParticles(i,j)%event(1)
                   hadrons(i,3)%mass=realParticles(i,j)%mass
                   hadrons(i,3)%position=realParticles(i,j)%position
                   hadrons(i,3)%momentum(1:3)=realParticles(i,j)%momentum(1:3)
                   times(i,3)=time
               end if
            end if
          exit
          end if
       end do ParticleLoop1

    end do

  end subroutine DoHadronAnalysisTime

  subroutine cleanUp
    if (allocated(times)) deallocate(times)
    if (allocated(hadrons)) deallocate(hadrons)
  end subroutine cleanUp

end module hadronAnalysis
