!******************************************************************************
!****m* /statistics
! NAME
! module statistics
! PURPOSE
! This module is used to store statistical information about the particles.
! Calling "saveInfo" saves the information of
! the particle with %number and ensemble being the ordering scheme.
! Calling "getInfo(ensemble,number,...)" returns the information.
! This info is not included in the particleDefinition, since it would blow up
! the memory usage. If this module is used, and only then,
! memory is allocated for this information.
!******************************************************************************
module statistics
  private

  type info
     integer :: number=0 ! number of particle
     real, dimension(1:3)      :: productionPlace=0.
     integer, dimension(1:3)   :: producingParticles=0
  end type info

  integer,save :: dim_InfoVec = 100
  integer,save :: maxSize = 900000 ! = 54 MB

  type(info),allocatable, save, dimension(:,:) :: infoVec

  logical, save :: initFlag=.true.
  logical, save :: debug=.false.
  logical, save :: initSave=.true.

  public:: saveInfo, getInfo,splashInfo

contains

  !****************************************************************************
  !****s* statistics/splashInfo
  ! NAME
  ! subroutine splashInfo
  ! INPUTS
  ! NONE
  ! PURPOSE
  ! This subroutine sets the module back to normal and deallocates the information vector, if allocated.
  !****************************************************************************
  subroutine splashInfo
    implicit none
    initSave=.true.
    initFlag=.true.
    if (allocated(infoVec)) deallocate(infoVec)
  end subroutine splashInfo


  !****************************************************************************
  !****s* statistics/saveInfo
  ! NAME
  ! subroutine saveInfo(number,ensemble,place, producingParticle1,producingParticle2,producingParticle3)
  ! PURPOSE
  ! This subroutine saves given information into a vector. The information is stored by number and ensemble of the particle.
  ! INPUTS
  ! integer, intent(in) :: number ! number of particle
  ! integer, intent(in) :: ensemble ! ensemble of particle
  ! integer, intent(in) :: producingParticle1,producingParticle2,producingParticle3 ! particle ID's which produced the particle
  ! real, dimension(1:3) , intent(in) :: place ! place where particle was produced
  ! AUTHOR
  ! Oliver Buss
  !****************************************************************************
  subroutine saveInfo(number,ensemble,place, producingParticle1,producingParticle2,producingParticle3)
    use inputgeneral, only: numensembles
    implicit none
    integer, intent(in) :: number ! number of particle
    integer, intent(in) :: ensemble ! ensemble of particle
    integer, intent(in) :: producingParticle1,producingParticle2,producingParticle3 ! particle ID's which produced the particle
    real, dimension(1:3) , intent(in) :: place ! place where particle was produced
    integer, save, dimension(:),Allocatable :: index

    if (initSave) then
       if (.not.allocated(index)) allocate(index(1:numEnsembles))
       index=1
       initSave=.false.
    end if

    if (debug) write(*,*) 'In Saving'
    ! Check that infoVec is allocated
    if (.not.allocated(infoVec)) call enlargeInfoVec

    ! Write input information to infoVec
    do
       if (ubound(infoVec,dim=2).ge.index(ensemble)) then
          infoVec(ensemble,index(ensemble))%number=number
          infoVec(ensemble,index(ensemble))%productionPlace=place
          infoVec(ensemble,index(ensemble))%producingParticles=(/producingParticle1,producingParticle2,producingParticle3/)
          index(ensemble)=index(ensemble)+1
          exit
       else
          call enlargeInfoVec
       end if
    end do
    if (debug) write(*,*) 'Done Saving'
  end subroutine saveInfo


  !****************************************************************************
  !****s* statistics/getInfo
  ! NAME
  ! subroutine getInfo(number,ensemble,prodPlace,prodParticles,successFlag)
  ! PURPOSE
  ! This subroutine returns given information out of a vector. The information is stored by number and ensemble of the particle.
  ! INPUTS
  ! integer, intent(in) :: number ! number of particle
  ! integer, intent(in) :: ensemble ! ensemble of particle
  ! integer, dimension(1:3), intent(out) :: producingParticles ! particles which produced the particle
  ! real, dimension(1:3) , intent(in) :: prodPlace ! place where particle was produced
  ! AUTHOR
  ! Oliver Buss
  !****************************************************************************
  subroutine getInfo(number,ensemble,prodPlace,prodParticles,successFlag)
    ! returns info for a given particle with %number "number" and ensemble "ensemble"
    use inputGeneral, only: fullensemble
    implicit none
    integer, intent(in)     :: number
    integer, intent(in)     :: ensemble
    real, dimension(1:3) ,intent(out)    :: prodPlace
    integer, dimension(1:3),intent(out)  :: prodParticles
    logical,  intent(out)                :: successFlag
    logical, save :: onceFlag=.true.
    integer :: i, ens

    successFlag=.false.
    prodPlace=0.
    prodParticles=0

    if (debug) write(*,*) 'in getInfo'

    if (.not.allocated(infoVec)) then
       if (onceFlag) write(*,*) 'Error in getInfo: No info stored'
       onceFlag=.false.
    else
       if (.not.fullEnsemble) then
          do i=lbound(infoVec,dim=2),ubound(infoVec,dim=2)
             if (infoVec(ensemble,i)%number.eq.number) then
                prodPlace=infoVec(ensemble,i)%productionPlace
                prodParticles=infoVec(ensemble,i)%producingParticles
                successFlag=.true.
             end if
          end do
       else
          ! search everywhere, since in fullensemble mode, particle can be stored everywhere
          do ens=lbound(infoVec,dim=1),ubound(infoVec,dim=1)
             do i=lbound(infoVec,dim=2),ubound(infoVec,dim=2)
                if (infoVec(ens,i)%number.eq.number) then
                   prodPlace=infoVec(ens,i)%productionPlace
                   prodParticles=infoVec(ens,i)%producingParticles
                   successFlag=.true.
                end if
             end do
          end do
       end if
    end if

    if (debug) write(*,*) 'done getInfo'

  end subroutine getInfo

  !****************************************************************************
  !****s* statistics/enlargeInfoVec
  ! NAME
  ! subroutine enlargeInfoVec
  ! PURPOSE
  ! This subroutine enlarges the vector infoVec if it got too small.
  ! If the size is already greater than maxsize, then the code stops to prevent memory overflow
  ! AUTHOR
  ! Oliver Buss
  !****************************************************************************
  subroutine enlargeInfoVec

    use inputGeneral, only: numEnsembles
    implicit none
    type(info),allocatable, dimension(:,:) :: save_infoVec
    integer :: i


    if (initFlag) then
       allocate(infoVec(1:numEnsembles,1:dim_InfoVec))
       initFlag=.false.
       write(*,*) '############################################################################################################'
       write(*,'(A,F9.3,A)') 'WARNING in  enlargeInfoVec: Enlarging infoVec. New size=', &
            & numensembles*dim_InfoVec*(16*3+4*4)/(1024.**2),'MB'
       write(*,*) '############################################################################################################'
       if (debug) write(*,*) 'enlargeInfoVec: Initialising'
       return
    end if


    if (debug) write(*,*) 'enlargeInfoVec: Enlarging'


    ! Check that size is not already too big
    if (size(infoVec,dim=2)*size(infoVec,dim=1).gt.maxSize) then
       write(*,*) 'Vector infoVec gets too big in statistics', &
            & lbound(infoVec,dim=2),ubound(infoVec,dim=2), lBound(infoVec,dim=1),uBound(infoVec,dim=1)
       write(*,*) 'dim_InfoVec=',dim_infoVec
       write(*,*) 'maxsize=', maxsize
       stop
    end if




    ! Save the old vector to another vector
    allocate(save_infoVec(1:numensembles,lbound(infoVec,dim=2):uBound(infoVec,dim=2)))
    do i=lbound(infoVec,dim=1),ubound(infoVec,dim=1)
       save_infoVec(i,:)=infoVec(i,:)
    end do

    ! Enlarge oldVector
    deallocate(infoVec)
    dim_InfoVec=dim_InfoVec*3
    allocate(infoVec(1:numEnsembles,1:dim_InfoVec))
    write(*,*) '############################################################################################################'
    write(*,'(A,F9.3,A)') 'WARNING in  enlargeInfoVec: Enlarging infoVec. New size=', &
         & numensembles*dim_InfoVec*(16*3+4*4)/(1024.**2),'MB'
    write(*,*) 'dim_infovec=', dim_infoVec
    write(*,*) '############################################################################################################'

    ! Write Information back
    do i=lbound(save_infoVec,dim=1),ubound(save_infoVec,dim=1)
       infoVec(i,lbound(save_infoVec,dim=2):uBound(save_infoVec,dim=2))=save_infoVec(i,:)
    end do
  end subroutine enlargeInfoVec




end module statistics
