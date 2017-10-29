!******************************************************************************
!****m* /annihilation_channels
! NAME
! module annihilation_channels
! PURPOSE
! Includes the routines for selecting annihilation channel.
!******************************************************************************

module annihilation_channels

  implicit none

  private

  !****************************************************************************
  !****t* annihilation_channels/anni_channel
  ! SOURCE
  !
  Type anni_channel
     integer,dimension(1:6) :: id=0      ! Id's of the final state particles
     integer,dimension(1:6) :: charge=0  ! charges of the final state particles
     real :: probability=0.              ! channel probability
  end type anni_channel
  !****************************************************************************


  !****************************************************************************
  !****g* annihilation_channels/pbarp_channels_atRest
  ! SOURCE
  !
  Type(anni_channel), Allocatable, dimension(:), save :: pbarp_channels_atRest
  ! PURPOSE
  ! Stores the information for all possible channels of antiproton-proton
  ! annihilation at rest.
  !****************************************************************************

  !****************************************************************************
  !****g* annihilation_channels/pbarp_channels_eDep
  ! SOURCE
  !
  Type(anni_channel), Allocatable, dimension(:), save :: pbarp_channels_eDep
  ! PURPOSE
  ! Stores the information for all possible channels of antiproton-proton
  ! annihilation at the given beam energy.
  !****************************************************************************


  !****************************************************************************
  !****g* annihilation_channels/pbarn_channels_atRest
  ! SOURCE
  !
  Type(anni_channel), Allocatable, dimension(:), save :: pbarn_channels_atRest
  ! PURPOSE
  ! Stores the information for all possible channels of antiproton-neutron
  ! annihilation at rest.
  !****************************************************************************


  !****************************************************************************
  !****g* annihilation_channels/pbarn_channels_eDep
  ! SOURCE
  !
  Type(anni_channel), Allocatable, dimension(:), save :: pbarn_channels_eDep
  ! PURPOSE
  ! Stores the information for all possible channels of antiproton-neutron
  ! annihilation at the given beam energy.
  !****************************************************************************

  integer, save :: n_pbarp_atRest, n_pbarn_atRest   ! Sizes of arrays pbarp_channels_atRest
                                                    ! and pbarn_channels_atRest respectively.
  integer, save :: n_pbarp_eDep, n_pbarn_eDep       ! Sizes of arrays pbarp_channels_eDep
                                                    ! and pbarn_channels_eDep respectively.

  public :: choose_channel

  logical, parameter :: verbose = .false.  ! overall verbosity flag

contains

  !****************************************************************************
  !****s* annihilation_channels/choose_channel
  ! NAME
  ! subroutine choose_channel(srts,antibar,bar,time,finalState,success)
  ! PURPOSE
  ! Selects randomly the outgoing channel (id's and charges only) for
  ! a given antibaryon-baryon annihilating pair.
  ! INPUTS
  ! * real, intent(in)           ::  srts                         ! inv. energy
  ! * type(particle), intent(in) ::  antibar                      ! antibaryon particle
  ! * type(particle), intent(in) ::  bar                          ! baryon particle
  ! * real, intent(in)           ::  time                         ! current time step (fm/c)
  ! OUTPUT
  ! * type(particle), intent(inout), dimension(:) :: finalState   ! Array of the final state particles
  !                                                               ! (only id's and charges are determined)
  ! * logical, intent(out) ::        success                      ! .true. if the channel choice was
  !                                                               ! successfull
  ! NOTES
  !****************************************************************************
    subroutine choose_channel(srts,antibar,bar,time,finalState,success)

    use particleDefinition, only: particle
    use IdTable, only: getAntiMeson
    use ParticleProperties, only: hadron
    use random
    use inputGeneral, only: eventtype
    use eventtypes, only: elementary

    real, intent(in)           ::  srts                         ! inv. energy
    type(particle), intent(in) ::  antibar                      ! antibaryon particle
    type(particle), intent(in) ::  bar                          ! baryon particle
    real, intent(in)           ::  time                         ! current time step (fm/c)
    type(particle), intent(inout), dimension(:) :: finalState   ! Array of the final state particles
                                                                ! (only id's and charges are determined)
    logical, intent(out) ::        success                      ! .true. if the channel choice was
                                                                ! successfull

    integer :: totalCharm, totalStrangeness, totalCharge, i, k, antiID, antiCharge
    real :: Prob, x, sum

    Type(anni_channel), Allocatable, dimension(:) :: pbarp_channels, pbarn_channels ! Working arrays

    integer :: n_pbarp, n_pbarn                                 ! Sizes of arrays pbarp_channels
                                                                ! and pbarn_channels respectively.

    logical, parameter :: flagEnergyDependent=.true.            ! If .true. -- use srts-dependent tables
                                                                ! of channels, .false. -- use tables for
                                                                ! annihilation at rest.
    real, parameter :: srts_min=1.876                           ! Minimum possible value of srts (GeV)
                                                                ! in annihilation.
    real, parameter :: srts_max=2.6                             ! Maximum value of srts up to which
                                                                ! channel tables at rest are still contributing
                                                                ! (above srts_max pure Unitary Model tables
                                                                !  are used).

    real, save :: srtsAverage=0., srtsAverage_previous=-100., time_previous=-100.
    integer, save :: nCalls=0
    logical, save :: flagIni=.true.

    totalCharm= hadron(bar%Id)%charm - hadron(antibar%Id)%charm

    if (totalCharm.ne.0) then  ! Annihilation with nonvanishing total charm
                              ! is not implemented
       success=.false.
       return
    end if

    totalStrangeness=  hadron(bar%Id)%strangeness - hadron(antibar%Id)%strangeness

    if (totalStrangeness.ne.0) then  ! Annihilation with nonvanishing total strangeness
                                    ! is not implemented yet here (see, however, DoColl_BaB)
       success=.false.
       return
    end if

    totalCharge= antibar%charge + bar%charge

    if (abs(totalCharge).gt.1) then
       success=.false.
       return
    end if

    if (flagIni) then
       call readInput
       flagIni=.false.
    end if

    if (allocated(pbarp_channels)) deallocate(pbarp_channels)
    if (allocated(pbarn_channels)) deallocate(pbarn_channels)

    if (.not.flagEnergyDependent) then  ! Use channel tables at rest

       n_pbarp=n_pbarp_atRest
       n_pbarn=n_pbarn_atRest
       allocate(pbarp_channels(1:n_pbarp),pbarn_channels(1:n_pbarn))
       pbarp_channels=pbarp_channels_atRest
       pbarn_channels=pbarn_channels_atRest

    else  !  Choose randomly between tables at rest and tables generated by Unitary Model

       if (time.ne.time_previous .or. eventtype.eq.elementary) then
          if (nCalls.gt.0) then
             srtsAverage=srtsAverage/float(nCalls)
          else
             srtsAverage=srts
          end if
          !write(*,'(1x,A26,1x,3(e13.6,1x))')' time, srts, srtsAverage: ', time, srts, srtsAverage
          if (srtsAverage.ne.srtsAverage_previous) call generate_table(srtsAverage)
          srtsAverage_previous=srtsAverage
          srtsAverage=0.
          nCalls=0
       end if
       time_previous=time

       nCalls=nCalls+1
       srtsAverage=srtsAverage+srts

       !write(*,'(1x,A21,1x,e13.6,1x,i5,1x,e13.6)') ' time, nCalls, srts: ', time, nCalls, srts

       ! "Probability" for the tables at rest:
       Prob=1.-(srts-srts_min)/(srts_max-srts_min)

       if (rn().lt.Prob) then
          n_pbarp=n_pbarp_atRest
          n_pbarn=n_pbarn_atRest
          allocate(pbarp_channels(1:n_pbarp),pbarn_channels(1:n_pbarn))
          pbarp_channels=pbarp_channels_atRest
          pbarn_channels=pbarn_channels_atRest
       else
          n_pbarp=n_pbarp_eDep
          n_pbarn=n_pbarn_eDep
          allocate(pbarp_channels(1:n_pbarp),pbarn_channels(1:n_pbarn))
          pbarp_channels=pbarp_channels_eDep
          pbarn_channels=pbarn_channels_eDep
       end if

    end if

    ! Monte-Carlo decision:

    x=rn()

    if (totalCharge.eq.0) then

       sum=0.
       i=0
       do
          i=i+1
          if (i.gt.n_pbarp) then
             write(*,*) 'In choose_channel pbarp: no channel found'
             stop
          end if
          sum=sum+pbarp_channels(i)%probability
          if (sum.gt.x) exit
       end do

       do k=1,6
          if (pbarp_channels(i)%Id(k).gt.0) then
             finalState(k)%Id=pbarp_channels(i)%Id(k)
             finalState(k)%charge=pbarp_channels(i)%charge(k)
          else
             exit
          end if
       end do

       success=.true.
       return

    else

       sum=0.
       i=0
       do
          i=i+1
          if (i.gt.n_pbarn) then
             write(*,*) 'In choose_channel pbarn: no channel found'
             stop
          end if
          sum=sum+pbarn_channels(i)%probability
          if (sum.gt.x) exit
       end do

       do k=1,6
          if (pbarn_channels(i)%Id(k).gt.0) then
             finalState(k)%Id=pbarn_channels(i)%Id(k)
             finalState(k)%charge=pbarn_channels(i)%charge(k)
          else
             exit
          end if
       end do

       if (totalCharge.eq.1) then
          do k=1,6
             if (finalState(k)%Id.gt.0) then
                call getAntiMeson(finalState(k)%Id,finalState(k)%charge,antiID,antiCharge)
                finalState(k)%Id=antiID
                finalState(k)%charge=antiCharge
             else
                exit
             end if
          end do
       end if

       success=.true.
       return

    end if

    end subroutine choose_channel


  !****************************************************************************
  !****s* annihilation_channels/generate_table
  ! NAME
  ! subroutine generate_table(srts)
  ! PURPOSE
  ! Interface subroutine to Igor Pshenichnov code which  generates the tables
  ! of outgoing channels for PbarP and PbarN annihilation.
  ! INPUTS
  ! * real, intent(in)           ::  srts                         ! inv. energy
  !****************************************************************************
    subroutine generate_table(srts)

    use twoBodyTools, only: p_lab
    use constants, only: mN
    use output, only: Write_InitStatus

    real, intent(in)           ::  srts                         ! inv. energy

    real :: PLAB,fnorm
    integer :: i,j,k,id,iz

    COMMON /CHGR/ICG(6,6000,3),NSUPL(3)
    integer :: ICG,NSUPL
    COMMON /WFFW/FWP(6000),FWP_N(6000),FWN(6000)
    real :: FWP,FWP_N,FWN
    SAVE /CHGR/, /WFFW/

    if (verbose) call Write_InitStatus('annihilation_channels/generate_table',0)

    PLAB=p_lab(srts,mN,mN)

    call UM_TN(PLAB)

    n_pbarp_eDep=NSUPL(1)
    n_pbarn_eDep=NSUPL(3)

    if (verbose) write(*,*) 'srts, plab:', srts, PLAB
    if (verbose) write(*,*) 'number of pbarp and pbarn annihilation channels :', &
              & n_pbarp_eDep,n_pbarn_eDep

    if (allocated(pbarp_channels_eDep)) deallocate(pbarp_channels_eDep)
    if (allocated(pbarn_channels_eDep)) deallocate(pbarn_channels_eDep)

    allocate(pbarp_channels_eDep(1:n_pbarp_eDep),pbarn_channels_eDep(1:n_pbarn_eDep))

    j=0
    do i=1,n_pbarp_eDep
       if (FWP(I).GT.9.999999E-05) then
         j=j+1
         do k=1,6
           call decode_to_GiBUU(ICG(k,i,1),id,iz)
           pbarp_channels_eDep(j)%Id(k)=id
           pbarp_channels_eDep(j)%charge(k)=iz
         end do
         pbarp_channels_eDep(j)%probability=FWP(i)
       end if
    end do

    fnorm= sum(pbarp_channels_eDep(1:n_pbarp_eDep)%probability)
    if (verbose) write(*,*) 'Norma pbarp: ', fnorm
    pbarp_channels_eDep(:)%probability= pbarp_channels_eDep(:)%probability/fnorm

    j=0
    do i=1,n_pbarn_eDep
       if (FWN(I).GT.9.999999E-05) then
         j=j+1
         do k=1,6
           call decode_to_GiBUU(ICG(k,i,3),id,iz)
           pbarn_channels_eDep(j)%Id(k)=id
           pbarn_channels_eDep(j)%charge(k)=iz
         end do
         pbarn_channels_eDep(j)%probability=FWN(i)
       end if
    end do

    fnorm= sum(pbarn_channels_eDep(1:n_pbarn_eDep)%probability)
    if (verbose) write(*,*) 'Norma pbarn: ', fnorm
    pbarn_channels_eDep(:)%probability= pbarn_channels_eDep(:)%probability/fnorm

    if (verbose) call Write_InitStatus('annihilation_channels/generate_table',1)

    end subroutine generate_table



  !****************************************************************************
  !****s* annihilation_channels/init
  ! NAME
  ! subroutine init
  ! PURPOSE
  ! Initializes the tables of outgoing channels for PbarP and PbarN annihilation
  ! at rest according to I.A. Pshenichnov (private communication, see also
  ! E.S. Golubeva et al., NPA 537, 393 (1992))
  !****************************************************************************
    subroutine readInput

    use inputGeneral, only: path_To_Input
    use output, only: Write_ReadingInput, Write_InitStatus

    integer :: ios, i, k, num, id_anni(1:6), iz, id
    CHARACTER*4 :: SYMB(1:6)
    real :: w, fnorm

    call Write_InitStatus('annihilation_channels',0)

    if (allocated(pbarp_channels_atRest)) deallocate(pbarp_channels_atRest)
    if (allocated(pbarn_channels_atRest)) deallocate(pbarn_channels_atRest)

    n_pbarp_atRest=200
    n_pbarn_atRest=200

    allocate(pbarp_channels_atRest(1:n_pbarp_atRest),pbarn_channels_atRest(1:n_pbarn_atRest))

    call Write_ReadingInput('PbarP annihilation table',0)
    ios=0
!    open(1,file=trim(path_to_input)//'/annihilation/um0_0.fnp',status='old',iostat=ios)
!    if (verbose) open(2,file='um0_0_chk.fnp',status='new')
    open(1,file=trim(path_to_input)//'/annihilation/empir_table.fnp',status='old',iostat=ios)
    if (verbose) open(2,file='empir_table_chk.fnp',status='unknown')

    i=0
    do
       read(1,405,end=5,err=10) num, id_anni(1:6), symb(1:6), w
       i=i+1
!       if (verbose) write(2,405) num, id_anni(1:6), symb(1:6), w
 405   FORMAT(1x,I5,1x,6(1X,I3),1x,6(1X,A4),1x,F8.4)
       do k=1,6
         call decode_to_GiBUU(id_anni(k),id,iz)
         pbarp_channels_atRest(i)%Id(k)=id
         pbarp_channels_atRest(i)%charge(k)=iz
       end do
       pbarp_channels_atRest(i)%probability=w
       cycle
5      exit
10     write(*,*) 'In annihilation_channels/readInput: error in reading empir_table.fnp'
       stop
    end do
    close(1)
    if (verbose) write(2,*)' Number of channels: ', i
    fnorm= sum(pbarp_channels_atRest(1:i)%probability)
    if (verbose) write(2,*)' Norma: ', fnorm
    pbarp_channels_atRest(:)%probability= pbarp_channels_atRest(:)%probability/fnorm
    if (verbose) close(2)
    call Write_ReadingInput('PbarP annihilation table',1)

    call Write_ReadingInput('PbarN annihilation table',0)
    ios=0
!    open(1,file=trim(path_to_input)//'/annihilation/um0_0.fnn',status='old',iostat=ios)
!    if (verbose) open(2,file='um0_0_chk.fnn',status='new')
    open(1,file=trim(path_to_input)//'/annihilation/empir_table.fnn',status='old',iostat=ios)
    if (verbose) open(2,file='empir_table_chk.fnn',status='unknown')
    i=0
    do
       read(1,405,end=15,err=20) num, id_anni(1:6), symb(1:6), w
       i=i+1
       if (verbose) write(2,405) num, id_anni(1:6), symb(1:6), w
       do k=1,6
         call decode_to_GiBUU(id_anni(k),id,iz)
         pbarn_channels_atRest(i)%Id(k)=id
         pbarn_channels_atRest(i)%charge(k)=iz
       end do
       pbarn_channels_atRest(i)%probability=w
       cycle
15     exit
20     write(*,*) 'In annihilation_channels/readInput: error in reading empir_table.fnn'
       stop
    end do
    close(1)
    if (verbose) write(2,*)' Number of channels: ', i
    fnorm= sum(pbarn_channels_atRest(1:i)%probability)
    if (verbose) write(2,*)' Norma: ', fnorm
    pbarn_channels_atRest(:)%probability= pbarn_channels_atRest(:)%probability/fnorm
    if (verbose) close(2)
    call Write_ReadingInput('PbarN annihilation table',1)

    call Write_InitStatus('annihilation_channels',1)

    end subroutine readInput



  !****************************************************************************
  !****s* annihilation_channels/decode_to_GiBUU
  ! NAME
  ! subroutine decode_to_GiBUU(id_input,id,iz)
  ! PURPOSE
  ! Transform the internal particle coding of Igor Pshenichnov's annihilation code
  ! to GiBUU coding.
  ! INPUTS
  ! * integer, intent(in)         ::  id_input             ! Id in Pshenichnov's coding
  ! OUTPUT
  ! * integer, intent(out)        ::  id                   ! Id in GiBUU coding
  ! * integer, intent(out)        ::  iz                   ! charge of particle
  ! NOTES
  ! Only mesons are considered.
  !****************************************************************************
    subroutine decode_to_GiBUU(id_input,id,iz)

    use IdTable

    integer, intent(in)         ::  id_input             ! Id in Pshenichnov's coding
    integer, intent(out)        ::  id                   ! Id in GiBUU coding
    integer, intent(out)        ::  iz                   ! charge of particle

    select case (id_input)

    case (1)
      id=pion
      iz=1
    case (2)
      id=pion
      iz=-1
    case (3)
      id=kaon
      iz=1
    case (4)
      id=kaonBar
      iz=-1
    case (5)
      id=kaon
      iz=0
    case (6)
      id=kaonBar
      iz=0
    case (7)
      id=pion
      iz=0
    case (8)
      id=eta
      iz=0
    case (9)
      id=etaPrime
      iz=0
    case (10)
      id=rho
      iz=1
    case (11)
      id=rho
      iz=-1
    case (12)
      id=kaonStar
      iz=1
    case (13)
      id=kaonStarBar
      iz=-1
    case (14)
      id=kaonStar
      iz=0
    case (15)
      id=kaonStarBar
      iz=0
    case (16)
      id=rho
      iz=0
    case (17)
      id=omegaMeson
      iz=0
    case (18)
      id=phi
      iz=0
    case (19)
      id=0
      iz=0
    case default
      write(*,*) 'In decode_to_GiBUU: wrong input channel'
    end select

    end subroutine decode_to_GiBUU

end module annihilation_channels
