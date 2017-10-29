!******************************************************************************
!****m* /InABoxAnalysisPion
! NAME
! module InABoxAnalysisPion
! PURPOSE
! Modules includes routines which are used to evaluate the decay width (gamma)
! of the pion.
!******************************************************************************
module InABoxAnalysisPion
  private

  character(40),save  :: name='dummy'

  public :: InABoxAnalysisPion_count
  public :: InABoxAnalysisPion_eval

contains

  !****************************************************************************
  !****s* InABoxAnalysisPion/InABoxAnalysisPion_count
  ! NAME
  ! subroutine InABoxAnalysisPion_count(particleVector,time)
  ! INPUTS
  ! * type(particle), dimension(:,:),intent(in) :: particleVector
  ! * real,intent(in) :: time
  ! PURPOSE
  ! Counts number of pions in the particleVector and prints it to file "pionNumbers_*.dat" where * is the kinetic energy of the
  ! incoming pions
  !****************************************************************************

  subroutine InABoxAnalysisPion_count(particleVector,time)

    use particleDefinition
    use IDTAble, only: pion
    use initPion, only: getEkin
    use output, only: realToChar
    use collisionNumbering, only: pert_numbering
    implicit none
    type(particle), dimension(:,:),intent(in) :: particleVector
    real,intent(in) :: time
    real, parameter :: epsilon=0.001
    real :: velocity
    integer :: numPions, numPions_notScattered
    integer :: i,j
    logical, save:: initFlag=.true.
    real, save :: ekin=0.

    if (abs(ekin-getEkin()).gt.0.0001) then ! new energy
       ekin=getEkin()
       initFlag=.true.
       name='pion_Numbers_'//realToChar(ekin*1000.)//'.dat'
    end if

    numPions=0
    numPions_notScattered=0
    velocity=0.
    do i=lbound(particleVector,dim=1), ubound(particleVector,dim=1)
       do j=lbound(particleVector,dim=2), ubound(particleVector,dim=2)
          if (particleVector(i,j)%ID.eq.pion) then
             numPions=numPions+1
             if (particleVector(i,j)%event(1).eq.pert_numbering()) then  ! Initial event in initPion is pert_numbering()
                numPions_notScattered=numPions_notScattered+1
                velocity=sqrt(Dot_Product(particleVector(i,j)%velocity(1:3),particleVector(i,j)%velocity(1:3)))
             end if
          end if
       end do
    end do

    if (initFlag) then
       open(321,File=name)
       initFlag=.false.
    else
       open(321,File=name,position='append')
    end if

    write(321,'(F12.6,2I6,F12.6)') time/197.,numPions, numPions_notScattered,velocity
    close(321)
  end subroutine InABoxAnalysisPion_count


  !****************************************************************************
  !****s* InABoxAnalysisPion/InABoxAnalysisPion_eval
  ! NAME
  ! subroutine InABoxAnalysisPion_eval()
  ! INPUTS
  ! NONE
  ! PURPOSE
  ! Reads in the number of pions out of the files created by
  ! InABoxAnalysisPion_count
  ! and evaluates the decay width.
  ! Prints results to "pion_gamma.dat" and to "pion_gamma_mean.dat".
  !****************************************************************************
  subroutine InABoxAnalysisPion_eval()

    use initPion, only: getEkin

    implicit none
    integer :: ios,number
    real :: gammaTotal, gammaAbs, gammaMean
    integer, dimension(1:2) :: num, num_notScattered
    real, dimension(1:2) :: time
    logical,save :: initFlag=.true.
    real :: velocity,meanFreePath,veloMean,GammaError
    real, parameter :: fmMeV=197.

    if (initFlag) then
       open(11,File='pion_gamma.dat')
       open(12,File='pion_gamma_mean.dat')
       initFlag=.false.
    else
       open(11,File='pion_gamma.dat',position='append')
       open(12,File='pion_gamma_mean.dat',position='append')
    end if

    open(13,File=name, status='old',iostat=ios)
    if (ios.ne.0) then
       ! ios.ne.0 = File not available or some other error in reading it.
       write(*,*) 'File pion_Numbers.dat is not available!  Critical error in InABoxAnalysisPion_eval. Return!'
       return
    else
       number=0
       gammaMean=0.

       num=0
       num_notScattered=0
       veloMean=0.
       gammaError=0.
       do
          num(1)=num(2)
          num_notScattered(1)=num_notScattered(2)
          time(1)=time(2)
          read(13,'(F12.6,2I6,F12.6)',ioStat=ios) time(2),num(2), num_notScattered(2),velocity
          if (ios.ne.0) exit
          ! Factor 197. due to delta_T in units of fm and gamma in units of MEV
          if (  (num(2).gt.0).and.(num(1).gt.0).and.(num_notScattered(2).gt.0).and.(num_notScattered(1).gt.0)&
               &.and.(time(2)-time(1) .gt.0.00001)) then
             gammaAbs=-(log(float(num(2)))-log(float(num(1))))/(time(2)-time(1))
             gammaTotal=-(log(float(num_notScattered(2)))-log(float(num_notScattered(1))))/(time(2)-time(1))
             write(11,*) gammaAbs, gammaTotal
             gammaMean=gammaMean+gammaTotal
             gammaError=gammaError+gammaTotal**2
             number=number+1
             veloMean=veloMean+velocity
          end if
       end do
       if (gammaMean.gt.0) then
          meanFreePath=veloMean/gammaMean*fmMeV
       else
          meanFreePath=10000000000000.
       end if
       ! Error of the mean value:
       if (number.gt.1) then
          gammaError=sqrt(max(((gammaError-gammaMean**2/float(Number))/float(number-1)/float(number)),0.))
       else
          gammaError=999.
       end if
       write(11,'(A,F10.5,A,3F10.5)') '# Elab=,',getekin(),&
            & ' Mean value of gamma, velocity, meanfreepath; error of gamma, error of mfp'&
            &, gammaMean/float(number) , veloMean/float(number), meanFreePath
       write(12,'(6F10.5)') getekin(),gammaMean/float(number),veloMean/float(number),meanfreepath, &
            & gammaError, GammaError/gammaMean*meanfreepath
       end if
    close(11)
    close(12)
    close(13)


  end subroutine InABoxAnalysisPion_eval



end module InABoxAnalysisPion
