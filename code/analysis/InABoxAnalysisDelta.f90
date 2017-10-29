!******************************************************************************
!****m* /InABoxAnalysisDelta
! NAME
! module InABoxAnalysisDelta
! PURPOSE
! Modules includes routines which are used to evaluate the decay width (gamma)
! of the pion.
!******************************************************************************
module InABoxAnalysisDelta
  private

  public :: InABoxAnalysisDelta_count
  public :: InABoxAnalysisDelta_eval

contains

  !****************************************************************************
  !****s* InABoxAnalysisDelta/InABoxAnalysisDelta_count
  ! NAME
  ! subroutine InABoxAnalysisDelta_count(particleVector,time)
  ! INPUTS
  ! * type(particle), dimension(:,:),intent(in) :: particleVector
  ! * real,intent(in) :: time
  ! PURPOSE
  ! Counts number of pions in the particleVector and prints it to file "pionNumbers_*.dat" where * is the kinetic energy of the
  ! incoming pions
  !****************************************************************************

  subroutine InABoxAnalysisDelta_count(part,time)

        ! Zählt alle Deltas, die noch nicht gestossen haben und schreibt die Zahl nach "numDeltas.dat".
    use particleDefinition
    use random
    use idTable, only: delta
    implicit none
    integer :: i,j
    type(particle), dimension(:,:) :: part
    real :: time
    integer :: numDeltas

    numDeltas=0
    do i=lbound(part,dim=1),ubound(part,dim=1)
       do j=lbound(part,dim=2),ubound(part,dim=2)
          if (part(i,j)%ID.eq.delta) then
             if (part(i,j)%event(1).eq.0) then
                numDeltas=numDeltas+1
             end if
          end if
       end do
    end do
    open(120,file="numDeltas.dat",position='append')
    write(120,*) time, numdeltas
    close(120)



  end subroutine InABoxAnalysisDelta_count


  !****************************************************************************
  !****s* InABoxAnalysisDelta/InABoxAnalysisDelta_eval
  ! NAME
  ! subroutine InABoxAnalysisDelta_eval()
  ! INPUTS
  !
  ! PURPOSE
  !
  !****************************************************************************
  subroutine InABoxAnalysisDelta_eval()



  end subroutine InABoxAnalysisDelta_eval



end module InABoxAnalysisDelta
