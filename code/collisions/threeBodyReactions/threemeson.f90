!******************************************************************************
!****m* /ThreeMeson
! NAME
! module ThreeMeson
!
! PURPOSE
! This module implements the back reaction 3 -> 1 of mesonic decay channels
!
! INPUTS
! (none)
!******************************************************************************
module ThreeMeson

  use callstack, only: traceBack

  implicit none

  private

  public :: DoThreeMeson
  public :: calcPhi3

  logical, save :: initFlag=.true.

  integer, dimension(1:4,1:20,1:2) :: ArrChannel
  integer, dimension(1:4) :: nChannel

  ! indicate same particles in every decay channel.
  ! bool coding: 1,2,4
  ! * part 1 and part 2: '1 or 2' = 3
  ! * part 1 and part 2 and part 3: '1 or 2 or 4' = 7
  ! Please keep up to date if changed in decayChannels
  integer, parameter, dimension(1:4) :: isSameBoolArr = (/ 3, 0, 7, 0 /)

  integer, dimension(1:4,1:20,1:2) :: StatChannel ! for statistics

  character(11),dimension(1:4), parameter :: sChannel = (/ &
       'pi0 pi0 eta', &
       'pi- pi+ pi0', &
       'pi0 pi0 pi0', &
       'pi- pi+ eta' /)

contains
  !****************************************************************************
  !****s* ThreeMeson/DoThreeMeson
  ! NAME
  ! subroutine DoThreeMeson(time)
  ! PURPOSE
  ! process all the finding of possibilities for 3->1 resonance decay back
  ! reactions. Actually also do the collision.
  !****************************************************************************
  subroutine DoThreeMeson(time)
    use constants, only: hbarc,mPi !,pi
    use decayChannels
    use IdTable, only: rho, omegaMeson, phi
    use inputGeneral, only: delta_T, numEnsembles
    use mediumDefinition
    use mesonPotentialModule, only: vecMes_massShift
    use mesonWidthMedium, only: decayWidthMesonMedium, WidthMesonMedium
    use particleDefinition
    use ParticleProperties
!    use twoBodyTools, only: pCM
    use VolumeElements, only: VolumeElements_boxSize, VolumeElements_3Body
    use random, only: rn

    real, intent(in) :: time

    integer :: ID,dID,iChannel

    logical :: doInit
    integer, dimension(1:3) :: iEns, iInd
    type(particle), POINTER :: Part1, Part2, Part3
    real :: srts, s, gammaIn, gammaTot, Phi3, m0
    real, dimension(nDecays) :: decayWidth
    logical :: pauliFlag
    type(medium) :: mediumATcoll
    real, dimension(0:3) :: momLRF
    real :: spectral, preFak, scaleFak, prob31

    if (initFlag) call ThreeMeson_Init

!!$    do i=1,200
!!$       srts = 2*mPi + 0.02d0*i
!!$       Phi3 = calcPhi3(srts,mPi)
!!$       write(*,*) srts,phi3,pCM(srts,mPi,mPi)/(4*pi*srts)
!!$    enddo
!!$    stop


    StatChannel = 0

    ! The probability for the 3->1 is given by:
    ! P_{31} = \frac{\Delta t}{(\Delta^3x)^2} \frac{1}{8E_1'E_2'E_3'}
    !          \frac{2\sqrt{s}\Gamma}{\Phi_3} 2\pi{\cal A}
    ! While the decay probability is independent on the number of
    ! test particles/ensembles, the 2 body interaction probability
    ! scales as the cross section linearly with this number,
    ! and the 3 to 1 probability scales with N_test squared.
    preFak = delta_T/(numEnsembles*VolumeElements_boxSize())**2 * hbarc**5

    do dID=1,nDecay3bodyMeson
!!$       write(*,*) '================= dID = ',dID


       doInit = .true.
       do
          if (.not.VolumeElements_3Body(Decay3BodyMeson(dID), isSameBoolArr(dID), &
               doInit, iEns,iInd, Part1,Part2,Part3, scaleFak, 1)) exit

          srts = sqrtS(Part1,Part2,Part3)
          s = srts**2
          select case (dID)
          case (2,3)
             Phi3 = CalcPhi3(srts,mPi)
          case (1,4)
             Phi3 = CalcPhi3(srts,hadron(102)%mass)
          end select

          do iChannel=1,nChannel(dID)
             ID = ArrChannel(dID,iChannel,1) ! ID of mother particle
!             write(*,*) 'ID = ',ID

             decayWidth = decayWidthMesonMedium(ID,srts,0,pauliFlag)
             gammaIn  = decayWidth(ArrChannel(dID,iChannel,2))

             gammaTot = WidthMesonMedium(ID,srts, momLRF, mediumATcoll)
             m0 = hadron(ID)%mass
             if (ID==rho .or. ID==omegaMeson .or. ID==phi) then
                ! take into account a possible in-medium mass shift
                m0 = m0 + vecMes_massShift(ID, mediumATcoll%density)
             end if
             ! 2\pi{/cal A} =
             spectral = 2*srts * gammaTot / ((s-m0**2)**2+s*gammaTot**2)

             ! P_{31} =
             prob31 = preFak &
                  * 1.0/( 8.0*Part1%momentum(0)*Part2%momentum(0)*Part3%momentum(0) ) &
                  * 2.0 * srts * gammaIn / Phi3 &
                  * spectral &
                  * scaleFak * (hadron(ID)%Spin*2+1)

             if (prob31 > rn()) then
!!$                write(*,*) 'Doing it...'
!!$                write(*,*) ID,gammaTot,gammaIn,prob31
                call ThreeMeson_DoColl(Part1,Part2,Part3, ID, time)

                StatChannel(dID,iChannel,1) = StatChannel(dID,iChannel,1) + 1
                exit ! this triple can not scatter any more !
             end if

          end do

       end do

    end do ! iDecay

    write(*,*) 'decay channels statistics:'
    do dID=1,nDecay3bodyMeson
       write(*,'(i4,A13,i4,50i8)') dID,sChannel(dID),nChannel(dID),ArrChannel(dID,1:nChannel(dID),1)
       write(*,'(A4,A13,A4,50i8)') ' ',' ',' ',          StatChannel(dID,1:nChannel(dID),1)
    end do

  end subroutine DoThreeMeson

  !****************************************************************************
  !****if* ThreeMeson/calcPhi3
  ! NAME
  ! real function calcPhi3(srts,m3)
  ! PURPOSE
  ! Calculate the 3-particle phase space. Particle 1 and 2 are pions.
  !
  ! INPUTS
  ! * real :: srts     --- sqrt(s) = mass of decaying particle
  ! * real :: m3       --- mass of third particle
  ! OUTPUT
  ! return value in GeV^2
  !
  ! NOTES
  ! there is threeBodyPhaseSpace/Integrate_3bodyPS, which does the same,
  ! but just by simply adding 100 values.
  ! This here is closer to mesonWidthVacuum/threepi
  !
  ! We calculate the integral of
  !  \dd\Phi_n(P;p_1,\dots,p_n) = (2\pi)^4\delta^{(4)}(P-p_1-\dots-p_n)
  !  \prod_i{\frac{\dd^3p_i}{(2\pi)^3\,2E_i}}
  !  =
  !  (2\pi)^4\delta^{(4)}(P-p_1-\dots-p_n)
  !   \prod_i{\frac{\dd^4p_i}{(2\pi)^4}}
  !   (2\pi)\delta(p_i^2-m_i^2)\Theta(p_i^0)
  !
  ! Here the definition of the n-body phase space is slightly (by a factor
  ! $(2\pi)^4$) modified compared to the PDG\cite{Beringer:1900zz} version.
  !****************************************************************************
  real function calcPhi3(srts,m3)

    use gauss_integration, only: sg20r, rg20r
    use constants, only: mPi,pi

    real, intent(in) :: srts
    real, intent(in) :: m3
    integer :: ns
    integer, parameter :: n=3
    real, dimension(1:20*n) :: absi,orde
    real, parameter :: mm = mPI**2
    real :: mm3,s,resu

    if (srts>2.*mPi+m3) then

       call sg20r(2.*mPi,srts-m3,n,absi,ns)
       absi = absi**2
       mm3 = m3**2
       s = srts**2

       orde = sqrt( ((absi-s-mm3)**2-4*s*mm3) * (absi-4*mm) )

       call rg20r(2.*mPi,srts-m3,n,orde,resu)
       calcPhi3 = resu/(64*pi**3*s)

    else
       calcPhi3 = 0
    end if

  end function calcPhi3

  !****************************************************************************
  !****is* ThreeMeson/ThreeMeson_INIT
  ! NAME
  ! subroutine ThreeMeson_INIT()
  ! PURPOSE
  ! This routine initializes the 3-meson routines.
  ! It sets up internal arrays.
  !****************************************************************************
  subroutine ThreeMeson_INIT()

    use decayChannels, only: nDecay3bodyMeson
    use IdTable, only: nMes, pion
    use ParticleProperties, only: nDecays, hadron
    use output, only: Write_InitStatus

    integer :: ID, iDecay, dID

    initFlag=.false. ! to remember that routine has been called

    call Write_InitStatus("ThreeMeson", 0)

    if (nDecay3bodyMeson /= 4) then
       write(*,*) 'nDecay3bodyMeson=',nDecay3bodyMeson
       call traceback('nDecay3bodyMeson /= 4')
    end if

    ArrChannel = 0
    nChannel = 0

    do ID=pion,pion+nMes-1
       do iDecay=1,nDecays
          dID=-hadron(ID)%decaysID(iDecay) ! 3Body are negative
          if (dID>0) then
             if (hadron(ID)%decays(iDecay) > 0.0) then
                nChannel(dID) = nChannel(dID)+1
                ArrChannel(dID,nChannel(dID),1:2) = (/ ID, iDecay /)
             end if
          end if
       end do
    end do

    write(*,*) 'decay channels to consider:'
    do dID=1,nDecay3bodyMeson
       write(*,'(i4,A13,i4,50i4)') dID,sChannel(dID),nChannel(dID),ArrChannel(dID,1:nChannel(dID),1)
       write(*,'(A4,A13,A4,50i4)') ' ',' ',' ',      ArrChannel(dID,1:nChannel(dID),2)
    end do

    call Write_InitStatus("ThreeMeson", 1)
  end subroutine ThreeMeson_INIT

  !****************************************************************************
  !****is* ThreeMeson/ThreeMeson_DoColl
  ! NAME
  ! subroutine ThreeMeson_DoColl(Part1,Part2,Part3, ID, time)
  ! PURPOSE
  ! This routine does the actual collision.
  ! It replaces Part1 with the mother particle and erases Part2 and Part3
  !****************************************************************************
  subroutine ThreeMeson_DoColl(Part1,Part2,Part3, ID, time)
    use collisionNumbering, only: real_numbering, real_firstnumbering, ReportEventNumber
    use history, only: setHistory
    use particleDefinition
    use twoBodyStatistics, only: rate

    type(particle), POINTER, intent(inout) :: Part1, Part2, Part3
    integer, intent(in) :: ID
    real, intent(in) :: time

    type(particle) :: PartNew

    !===== Create new particle =====

    PartNew%ID = ID
!    PartNew%antiparticle = .false.

    call setNumber(PartNew) ! we give it a new number

    PartNew%charge = Part1%charge + Part2%charge + Part3%charge

    PartNew%position = (Part1%position+Part2%position+Part3%position)/3
    PartNew%momentum = (Part1%momentum+Part2%momentum+Part3%momentum)
    PartNew%velocity = Part1%momentum(1:3)/Part1%momentum(0)

    PartNew%mass = sqrtS(Part1, Part2, Part3)

    call setHistory(Part1,Part2,Part3, (/PartNew/) )

    ! no formation time effects:
    PartNew%lastCollisionTime = time
    PartNew%productionTime    = time
    PartNew%formationTime     = time
!    PartNew%scaleCS=1.
!    PartNew%in_Formation=.false.

    ! particles are not perturbative:
!    PartNew%perturbative = .false.
!    PartNew%perWeight = 1.0

    ! particles are 'on-shell':
!    PartNew%offshellParameter = 0.

    ! now the tricky part starts:
    PartNew%event(1:2) = real_numbering()
    PartNew%firstEvent = real_firstnumbering()

    !===== Set some global counters =====

    call reportEventNumber( (/Part1,Part2,Part3/), (/PartNew/), PartNew%event, time, 3111)

    call rate( (/Part1,Part2,Part3/), (/PartNew/), time)

    !===== Set particle 1, Delete particle 2 & 3 =====

    Part1 = PartNew
    call setToDefault(Part2)
    call setToDefault(Part3)

  end subroutine ThreeMeson_DoColl

end module ThreeMeson
