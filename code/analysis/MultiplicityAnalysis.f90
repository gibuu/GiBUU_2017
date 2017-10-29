!******************************************************************************
!****m* /MultiplicityAnalysis
! NAME
! module MultiplicityAnalysis
!
! PURPOSE
! This module contains collision independent analysis routines concerning
! multiplicity of particles
!
! INPUTS
! (none)
!
!******************************************************************************
module MultiplicityAnalysis

  use AnaEventDefinition
  implicit none

  private

  integer,parameter :: nNBins = 37
  integer,parameter :: nPBins = 18

  real,dimension(nPBins,0:nNBins),save :: MultArr
  real,dimension(nPBins,0:2),save      :: MultAveArr

  character*(13),dimension(nPBins) :: BinName = (/ &
       & "ALL         ", &
       & "ALL charged ", &
       & "charge -    ", &
       & "charge 0    ", &
       & "charge +    ", &
       & "ALL pions   ", &
       & "pi-         ", &
       & "pi0         ", &
       & "pi+         ", &
       & "ALL baryons ", &
       & "neutrons    ", &
       & "protons     ", &
       & "ALL kaons   ", & ! and antikaons
       & "K-          ", &
       & "K0bar       ", &
       & "K0          ", &
       & "K+          ", &
       & "Lambda      " /)

  public :: Multiplicity_Reset
  public :: Multiplicity_AddEvent
  public :: Multiplicity_Write

contains
  !****************************************************************************
  !****s* MultiplicityAnalysis/Multiplicity_Reset
  ! NAME
  ! subroutine Multiplicity_Reset
  !
  ! PURPOSE
  !
  !****************************************************************************
  subroutine Multiplicity_Reset
    MultArr = 0
    MultAveArr = 0
  end subroutine Multiplicity_Reset


  !****************************************************************************
  !****s* MultiplicityAnalysis/Multiplicity_AddEvent
  ! NAME
  ! subroutine Multiplicity_AddEvent(E)
  !
  ! PURPOSE
  !
  !****************************************************************************
  subroutine Multiplicity_AddEvent(E)
    use IdTable, only: isBaryon
    use particlePointerListDefinition
    use particlePointerList

    type(tAnaEvent),intent(in) :: E

    type(tParticleListNode),Pointer  :: pNode
    integer, dimension(nPBins) :: nPart
    integer :: i
    real :: weight

    nPart = 0
    weight = 0

    pNode => E%particleList%first
    do
       if (.not. associated(pNode)) exit
       if (pNode%V%ID.gt.0) then

          weight = pNode%V%perweight

          nPart(1) = nPart(1)+1   ! == All
          if (pNode%V%charge.lt.0) then
             nPart(2)=nPart(2)+1  ! == all charged
             nPart(3)=nPart(3)+1  ! == negative charged
          else if (pNode%V%charge.gt.0) then
             nPart(2)=nPart(2)+1  ! == all charged
             nPart(5)=nPart(5)+1  ! == positive charged
          else
             nPart(4)=nPart(4)+1  ! == neutral
          end if

          if (pNode%V%ID.eq.101) then
             nPart(6) = nPart(6)+1   ! == All pions
             if (pNode%V%charge.lt.0) then
                nPart(7)=nPart(7)+1  ! == negative charged pions
             else if (pNode%V%charge.gt.0) then
                nPart(9)=nPart(9)+1  ! == positive charged pions
             else
                nPart(8)=nPart(8)+1  ! == neutral pions
             end if
          end if


          if (pNode%V%ID.eq.111) then
             nPart(13) = nPart(13)+1   ! == All kaons and __antikaons__
             if (pNode%V%charge.lt.0) then
                nPart(14)=nPart(14)+1  ! == negative charged antikaons
             else if (pNode%V%charge.gt.0) then
               write(*,*) 'Problem in MultiplicityAnalysis.f90;  positively charged antikaon found'
               stop
             else
                nPart(15)=nPart(15)+1  ! == neutral antikaons
             end if
          end if
          if (pNode%V%ID.eq.110) then
             nPart(13) = nPart(13)+1   ! == All __kaons__ and antikaons
             if (pNode%V%charge.lt.0) then
               write(*,*) 'Problem in MultiplicityAnalysis.f90;  negatively charged kaon found'
               stop
             else if (pNode%V%charge.gt.0) then
                nPart(17)=nPart(17)+1  ! == positive charged kaons
             else
                nPart(16)=nPart(16)+1  ! == neutral kaons
             end if
          end if


          if (isBaryon(pNode%V%ID) .and..not.pNode%V%antiParticle) then
             nPart(10) = nPart(10)+1   ! == All baryons
             if (pNode%V%ID.eq.32) nPart(18)=nPart(18)+1  ! Lambda baryons
             if (pNode%V%ID.eq.1) then
                if (pNode%V%charge.eq.0) then
                   nPart(11)=nPart(11)+1  ! == neutrons
                else if (pNode%V%charge.eq.1) then
                   nPart(12)=nPart(12)+1  ! == protons
                end if
             end if
          end if

       end if
       pNode => pNode%next
    end do

!    write(*,*) nPart
!    write(*,*)

    do i=1,nPBins
       MultAveArr(i,0)=MultAveArr(i,0)+weight
       MultAveArr(i,1)=MultAveArr(i,1)+weight*nPart(i)
       MultAveArr(i,2)=MultAveArr(i,2)+weight*(nPart(i))**2
       if (nPart(i).le.nNBins) MultArr(i,nPart(i))=MultArr(i,nPart(i))+weight
    end do

  end subroutine Multiplicity_AddEvent

  !****************************************************************************
  !****s* MultiplicityAnalysis/Multiplicity_Write
  ! NAME
  ! subroutine Multiplicity_Write(Prefix,iFile)
  !
  ! PURPOSE
  !
  !****************************************************************************
  subroutine Multiplicity_Write(Prefix,iFile)
    use output

    character*(*),intent(in)   :: Prefix
    integer,OPTIONAL,intent(in):: iFile

    integer :: iF,iL,i,j
    logical, save :: DoInit = .true.

    iF = 732 ! DUMMY value
    if (present(iFile)) iF = iFile

    iL = len(Prefix)

    !--- Write average and dispersion ---

    if (DoInit) then
       open(iF,file='Multiplicity_Ave.dat',status='unknown')
       write(iF,'("#",A,100("   ",A13,10("-")))') &
            & line(1:iL),BinName
    else
       open(iF,file='Multiplicity_Ave.dat',position='append')
    end if

    write(iF,'(" ",A,1P,100(e13.4))') Prefix,(f1(i),f2(i),i=1,nPBins)
    close(iF)

    !--- Write distribution ---

    if (DoInit) then
       open(iF,file='Multiplicity_Distr.dat',status='unknown')
       write(iF,'("#     ",100(A))') BinName
       write(iF,'("#")')
    else
       open(iF,file='Multiplicity_Distr.dat',position='append')
    end if

    write(iF,'("# ",A)') Prefix
    do i=0,nNBins
       write(iF,'(i3,1P,100(e13.4))') i,(g1(j,i),j=1,nPBins),&
            & (f1(j),j=1,nPBins)
    end do
    write(iF,*)
    write(iF,*)
    close(iF)

    DoInit = .false.

  contains

    real function f1(i)
      integer :: i
      if (MultAveArr(i,0).ne.0) then
         f1 = MultAveArr(i,1)/MultAveArr(i,0)
      else
         f1 = 99.9
      end if
      return
    end function f1

    real function f2(i)
      integer :: i
      if (MultAveArr(i,0).ne.0) then
         f2 = MultAveArr(i,2)/MultAveArr(i,0) - &
              & (MultAveArr(i,1)/MultAveArr(i,0))**2
         f2 = sqrt(max(0.0,f2))
      else
         f2 = 99.9
      end if
      return
    end function f2

    real function g1(i,j)
      integer :: i,j
      if (MultAveArr(i,0).ne.0) then
         g1 = MultArr(i,j)/MultAveArr(i,0)
      else
         g1 = 99.9
      end if
      return
    end function g1

  end subroutine Multiplicity_Write


       !***********************************************************************
       !****o* MultiplicityAnalysis/Multiplicity_Ave.dat
       ! NAME
       ! file Multiplicity_Ave.dat
       !
       ! PURPOSE
       ! The file shows the average multiplicity and its standard deviation for the various outging  hadrons
       !
       ! Columns:
       ! * #1: running variable (see Multiplicity_Distr.dat for explanation)
       ! * #2:  average multiplicity  of  "ALL"           particles
       ! * #3:  stadard deviation (of average multiplicity)  of  "ALL"           particles
       ! * #4:  average multiplicity  of  "ALL CHARGED"   particles
       ! * #5:  stadard deviation (of average multiplicity)  of  "ALL CHARGED"   particles
       ! * #6:  average multiplicity  of  "charge -1"     particles
       ! * #7:  stadard deviation (of average multiplicity)  of  "charge -1"     particles
       ! * #8:  average multiplicity  of  "charge 0"      particles
       ! * #9:  stadard deviation (of average multiplicity)  of  "charge 0"      particles
       ! * #10:  average multiplicity  of  "charge +1"     particles
       ! * #11:  stadard deviation (of average multiplicity) of  "charge +1"     particles
       ! * #12:  average multiplicity  of  "ALL pions"
       ! * #13:  stadard deviation (of average multiplicity) of  "ALL pions"
       ! * #14:  average multiplicity  of  "pion -1"
       ! * #15:  stadard deviation (of average multiplicity) of  "pion -1"
       ! * #16:  average multiplicity  of  "pion 0"
       ! * #17:  stadard deviation (of average multiplicity) of  "pion 0"
       ! * #18:  average multiplicity  of  "pion +1"
       ! * #19:  stadard deviation (of average multiplicity) of  "pion +1"
       ! * #20:  average multiplicity  of  "ALL nucleons"
       ! * #21:  stadard deviation (of average multiplicity) of  "ALL nucleons"
       ! * #22:  average multiplicity  of  "neutrons"
       ! * #23:  stadard deviation (of average multiplicity) of  "neutrons"
       ! * #24:  average multiplicity  of  "protons"
       ! * #25:  stadard deviation (of average multiplicity) of  "protons"
       ! * #26:  average multiplicity  of  "ALL kaons"
       ! * #27:  stadard deviation (of average multiplicity) of  "ALL kaons"
       ! * #28:  average multiplicity  of  "Kbar -1"
       ! * #29:  stadard deviation (of average multiplicity) of  "Kbar -1"
       ! * #30:  average multiplicity  of  "Kbar 0"
       ! * #31:  stadard deviation (of average multiplicity) of  "Kbar 0"
       ! * #32:  average multiplicity  of  "K 0"
       ! * #33:  stadard deviation (of average multiplicity) of  "K 0"
       ! * #34:  average multiplicity  of  "K +1"
       ! * #35:  stadard deviation (of average multiplicity) of  "K +1"
       ! * #36:  average multiplicity  of  "Lambda"
       ! * #37:  stadard deviation (of average multiplicity) of  "Lambda"
       !
       !***********************************************************************








       !***********************************************************************
       !****o* MultiplicityAnalysis/Multiplicity_Distr.dat
       ! NAME
       ! file Multiplicity_Distr.dat
       ! PURPOSE
       ! The file shows the multiplicity distributions of various hadrons in the final state
       !
       ! The commented line with  one number only (looks like #   1.000) is  the "running variable" from your jobcard,
       ! For example, for neutrino mode runs versus neutrino energy (mode 0 or mode 6) this is neutrino energy
       !  for neutrino mode runs versus Q2 (mode 3) this is Q2
       ! For each value of "running variable" you have several (="num_runs_sameEnergy") blocks of output separated by two empty lines
       !
       ! Columns:
       ! * #1: number of outgoing particles (integer)
       ! * #2:  percentage of events with number of "ALL"         particles equal to #1
       ! * #3:  percentage of events with number of "ALL CHARGED" particles equal to #1
       ! * #4:  percentage of events with number of "charge -1"   particles equal to #1
       ! * #5:  percentage of events with number of "charge 0"    particles equal to #1
       ! * #6:  percentage of events with number of "charge +1"   particles equal to #1
       ! * #7:  percentage of events with number of "ALL pions"    equal to #1
       ! * #8:  percentage of events with number of "pion -1"      equal to #1
       ! * #9:  percentage of events with number of "pion  0"      equal to #1
       ! * #10: percentage of events with number of "pion +1"      equal to #1
       ! * #11: percentage of events with number of "ALL baryons"  equal to #1
       ! * #12: percentage of events with number of "neutrons"     equal to #1
       ! * #13: percentage of events with number of "protons"      equal to #1
       ! * #14: percentage of events with number of "ALL kaons"   (and antikaons)  equal to #1
       ! * #15: percentage of events with number of "Kbar -1"      equal to #1
       ! * #16: percentage of events with number of "Kbar 0"       equal to #1
       ! * #17: percentage of events with number of "K0"           equal to #1
       ! * #18: percentage of events with number of "K+"           equal to #1
       ! * #19: percentage of events with number of "Lambda"       equal to #1
       ! * #20: average multiplicity of  "ALL"           (here column 1 makes not sense, so the values are the same)
       ! * #21: average multiplicity of  "ALL CHARGED"   (here column 1 makes not sense, so the values are the same)
       ! * #22: average multiplicity of  "charge -1"     (here column 1 makes not sense, so the values are the same)
       ! * #23: average multiplicity of  "charge 0"      (here column 1 makes not sense, so the values are the same)
       ! * #24: average multiplicity of  "charge +1"     (here column 1 makes not sense, so the values are the same)
       ! * #25: average multiplicity of  "ALL pions"     (here column 1 makes not sense, so the values are the same)
       ! * #26: average multiplicity of  "pion -1"       (here column 1 makes not sense, so the values are the same)
       ! * #27: average multiplicity of  "pion  0"       (here column 1 makes not sense, so the values are the same)
       ! * #28: average multiplicity of  "pion +1"       (here column 1 makes not sense, so the values are the same)
       ! * #29: average multiplicity of  "ALL baryons"   (here column 1 makes not sense, so the values are the same)
       ! * #30: average multiplicity of  "neutrons"      (here column 1 makes not sense, so the values are the same)
       ! * #31: average multiplicity of  "protons"       (here column 1 makes not sense, so the values are the same)
       ! * #32: average multiplicity of  "ALL kaons"     (here column 1 makes not sense, so the values are the same)
       ! * #33: average multiplicity of  "Kbar -1"       (here column 1 makes not sense, so the values are the same)
       ! * #34: average multiplicity of  "Kbar 0"        (here column 1 makes not sense, so the values are the same)
       ! * #35: average multiplicity of  "K0"            (here column 1 makes not sense, so the values are the same)
       ! * #36: average multiplicity of  "K+"            (here column 1 makes not sense, so the values are the same)
       ! * #37: average multiplicity of  "Lambda"        (here column 1 makes not sense, so the values are the same)
       !
       ! NOTES
       ! The numbers in one raw should NOT add-up in any way
       !
       ! Example: consider an event with 1 neutron, 1 pi+ and 1pi0 in the final state;
       ! it will contribute to  "3 ALL", "1 ALL CHARGED", "0 charge -1", "2 charge 0", "1 charge +1",
       !                        "2 ALL pions", "0 pion -1", "1 pion 0", "1 pion+", "1 ALL baryons", "1 neutrons", "0 protons"
       !
       ! In each column (in each block) the probabilities should sum up to 1
       !***********************************************************************

end module MultiplicityAnalysis
