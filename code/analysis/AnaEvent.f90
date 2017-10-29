!******************************************************************************
!****m* /AnaEvent
! NAME
! module AnaEvent
!
! PURPOSE
! This module handles the type "tAnaEvent".
! It can be used for analyzing your particle vector. See e.g.
! lowElectronAnalysis.f90 for usage.
!
! This module includes several subroutines to print out cross sections
! when given an array of events.
!
! The routines given in this module belong to two categories:
! The first category are the more low-level routines, which act on building
! up an single tAnaEvent or dumping a list of tAnaEvent types, as e.g.:
! * event_init
! * event_clear
! * event_add
! * event_dump
! On the other hand we have routines as:
! * event_sigma
! * event_dSigma_dE
! * event_dSigma_dOmega
! * event_dSigma_dEcostheta
! which fill histograms by looping over a given list of tAnaEvent types. These
! are very complex.
!
! INPUTS
! No Namelist available.
!
! NOTES
!
! %perweight should not include any dOmega or the such.
! The sum of all perweight should simply return the total cross section.
!
!******************************************************************************
module AnaEvent

  use AnaEventDefinition
  use CallStack, only: TraceBack

  implicit none

  private

  logical :: debug=.false.
  logical :: exclusive_hadron
  real :: nuc_mom_pcut
  logical :: Q2p

  integer, parameter, public:: dimSigma=120      ! number of possible channels

  !****************************************************************************
  !****ig* AnaEvent/particleIDs_flag
  ! SOURCE
  !
  logical, dimension(1:numStableParts) :: particleIDs_flag = .true.
  ! PURPOSE
  ! Flags to switch on/off output of special ID's,
  ! e.g. particleIDs_flag(1)=.false. means that there is no output for pions
  ! The array particleIDs is set in module AnaEventDefinition.f90
  !****************************************************************************

  public :: set_particleIDs_flag
  public :: event_init,event_add,event_clear, event_dump
  public :: event_sigma
  public :: event_dSigma_dE,event_dSigma_dOmega,event_dSigma_dEcostheta,&
            &event_dSigma_dInvMass
  public :: event_pairPhotons
  public :: event_hadronicEnergy, event_dSigma_dLeptonVariables
  public :: event_GetParticle
  public :: makeerror_hist, seterror_hist ! used in LArAnalysis.f90
  public :: SpecificEvent_Name, IfPass_SpecificEvent
  public :: set_Exclusive 
  public :: set_QElike


contains

  !****************************************************************************
  !****s* AnaEvent/set_particleIDs_flag
  ! NAME
  ! subroutine set_particleIDs_flag(flags)
  !
  ! PURPOSE
  ! Set the variable set_particleIDs_flag to a new value
  !
  ! INPUTS
  ! * logical, dimension(1:numStableParts) :: flags
  !
  ! OUTPUT
  ! * sets     particleIds_flag=flags
  !****************************************************************************
  subroutine set_particleIDs_flag(flags)
    use IdTable
    logical, dimension(1:numStableParts), intent(in) :: flags
    particleIds_flag=flags
  end subroutine set_particleIDs_flag


!   !*************************************************************************
!   !****s* AnaEvent/event_getMultiplicities
!   ! NAME
!   ! logical function event_getMultiplicities(L,ID,Charge,N)
!   !
!   ! PURPOSE
!   ! Output the multiplicities of an event.
!   !
!   ! INPUTS
!   ! * integer      :: ID      ! ID of particle
!   ! * integer      :: Charge  ! Charge of particle
!   ! * type(tAnaEvent) :: L       ! The event
!   !
!   ! OUTPUT
!   ! * integer :: N  -- Multiplicity of particle species "ID" with charge
!   !   "charge" in the event L
!   ! * function value --
!   !   true if "ID" was counted seperately,
!   !   false if "ID"-Multiplicities were included in the general
!   !   multiplicities, then "N" ist the number of all particles not being
!   !   defined in the array particleIDs (see above).
!   !
!   ! NOTES
!   ! * This routine is not suited to return the multiplicities of
!   !   anti-particles.
!   ! * We count explicitly all "particles" (NOT ANTIPARTICLES) with IDs given
!   !   in the array "particleIDs".
!   !*************************************************************************
!   logical function event_getMultiplicities(L,ID,Charge,N)
!     use IdTable
!
!     integer, intent(in)     :: ID
!     integer, intent(in)     :: Charge
!     type(tAnaEvent),intent(in) :: L
!
!     integer, intent(out) :: N
!
!     integer :: i
!
!     Select Case(ID)
!     case(pion,eta,kaon,kaonBar,dMeson,dBar,nucleon,lambda,sigmaResonance,Xi, &
!        &  OmegaResonance, &
!        &  ds_plus,ds_minus)
!        do i=lBound(particleIDs,dim=1),uBound(particleIDs,dim=1)
!           if(particleIDs(i).eq.ID) exit
!        end do
!        N=L%numberParticles(i,charge)
!        event_getMultiplicities=.true.
!     case Default
!        N=L%numberParticles(numStableParts+1,charge)
!        event_getMultiplicities=.false.
!     end select
!   end function event_getMultiplicities




  !****************************************************************************
  !****s* AnaEvent/event_sigma
  ! NAME
  ! subroutine event_sigma(E,sigma,newInit,runNumber,identifier)
  !
  ! PURPOSE
  ! * Output for single particle production  sigma(total), e.g. pi+, pi-
  ! * Output for double particle production  sigma(total), e.g. 2pi, 2N
  ! * Output for three particle production  sigma(total), e.g. 3pi, 3N
  ! * Output for four particle production  sigma(total), e.g. 4pi, 4N
  !
  !
  ! INPUTS
  ! * type(tAnaEvent), dimension(:) :: E -- List of events
  ! * real, dimension(1:dimSigma,1:2),intent(inout),optional :: sigma --
  !   List of Xsection which is being updated in this routine
  !   First index: : channel
  !   Second index : 1=sum over runs( Xsection of each run),
  !   2=sum over runs(( Xsection of each run)**2)
  ! * logical,intent(inout),optional :: newInit
  ! * integer,intent(inout),optional :: runNumber
  ! * real, optional :: identifier --
  !   Plotted in column 1 of the output file,
  !   if not present than we plot runNumber there
  !
  ! The optional input variable sigma returns the evaluated cross sections.
  !
  ! If newInit=.true., then the input sigma is initialized again (=0).
  ! If newInit=.false., then the input is used as a result of earlier runs and
  ! the results of the present run are just added.
  !
  ! If .true. then we assume that runNumber is the number to use as a
  ! normalization for the histogram results;
  ! We assume that the Xsection contain the results of "runNumber" runs.
  !
  !
  ! NOTES
  ! * If there are more than one particle of the same species produced,
  !   then the event is not considered for the "single particle cross section".
  !   e.g. a pi^+ pi^0 event does not contribute to any of those cross
  !   sections.
  ! * If there are more than two particles of the same species produced,
  !   then the event is not considered for the "double particle cross section".
  ! * If there are more than three particles of the same species produced,
  !   then the event is not considered for the "three particle cross section".
  ! * If there are more than four particles of the same species produced,
  !   then the event is not considered for the "four particle cross section".
  !
  ! Note that e.g. 2pi/3pi means that two/three pions, no matter which
  ! charge, where found.
  !
  ! Note also that, e.g., "pi+ X" means that one piPlus was found. And if
  ! there is another particle of the same species, but with another charge,
  ! then we still count the event.
  ! E.g. for pi^+ pi^0 the event contributes both to "pi+ + X" and "pi0 + X"
  !
  ! OUTPUT
  ! * filename 'sigma.dat'
  !****************************************************************************
  subroutine event_sigma(E,sigmaOut,newInit,runNumber,identifier)
    use IdTable, only: isHadron
    use particleProperties, only: hadron, validCharge_ID, PartName
    use output, only: intToChar
    use particleDefinition
    use initNeutrino, only: get_init_namelist


    ! Input
    real, optional :: identifier
    type(tAnaEvent), intent(in), dimension(:) :: E
    real, dimension(1:dimSigma,1:2),intent(inout),optional :: sigmaOut
    logical,intent(in),optional :: newInit
    integer,intent(in),optional :: runNumber
    real, dimension(1:4) :: perweight

    ! Field to store the Xsections in
    real, dimension(1:dimSigma,1:2) :: sigma

    character(30), dimension(1:dimSigma,1:2) :: sigma_name
    ! Denotes the different output channels
    integer :: channel, charge, i,j
    type(particle) :: part
    character(20) :: name
    character(20) :: formatOut
    character(100) :: fname
    integer :: foundIds
    logical :: initSigma=.true.
    integer :: numberParticles, num

    integer :: outLeptonID, outLeptoncharge

    call get_init_namelist(outLepton_ID=outLeptonID, &
         & outLepton_charge=outLeptoncharge)

    channel=0
    sigma_name(:,1)=' '
    sigma_name(:,2)=' '
    if (present(sigmaOut).and.present(runNumber).and.present(newInit)) then
       if (newInit) then
          sigma(:,1)=0.
          sigma(:,2)=0.
          sigmaOut(:,1)=0.
          sigmaOut(:,2)=0.
       else
          sigma(:,1)=sigmaOut(:,1)
          sigma(:,2)=sigmaOut(:,2)
       end if
       channel=channel+1
       sigma_name(channel,1)=intTochar(channel)//': run'
       sigma_name(channel,2)=intTochar(channel)//'     '
       if (present(identifier)) then
          sigma(channel,1)=identifier
       else
          sigma(channel,1)=runNumber
       end if
       sigma(channel,2)=0
    else
       sigma(:,1)=0.
       sigma(:,2)=0.
    end if
    if (channel>dimSigma) then
       call TraceBack('index value channel > bounds')
    end if

    !**************************************************************************
    !
    ! One particle sigmas
    ! - no sum over the different charge states
    !
    !**************************************************************************
    idLoop: do i=lbound(particleIDs,dim=1),ubound(particleIDs,dim=1)
       chargeLoop: do charge=-2,2
          if (.not.validCharge_ID(particleIDs(i),charge)) cycle
          ! Check that charge is valid
          channel=channel+1
          if (channel>dimSigma) then
             call TraceBack('index value channel > bounds')
          end if

          ! (0) Create the name of the cross section channel to be stored
          ! in sigma_name
          name = PartName(particleIDs(i),charge,.false.)
          if (i <= numStableMesons .and. exclusive_hadron) &
              &  name = trim('excl'//name)
          sigma_name(channel,1)=intTochar(channel)//': '// trim(name)
          sigma_name(channel,2)=intTochar(channel+dimSigma)//': error '&
              & // trim(name)

             ! (1) Sum over all events
          eventLoop: do j=lbound(E,dim=1),ubound(E,dim=1)
             ! (2) Count contributions of any particles in the event with
             !     charge="charge"
             !     and ID="particleIDs(i):

           if (( E(j)%numberParticles(i,charge).eq.1)  .and.  &
                & ( sum(E(j)%numberParticles(1:numStableMesons,:)).eq.1)) then

                if ( event_GetParticle(E(j),particleIDs(i),charge,1,part)) then
                   sigma(channel,1)=sigma(channel,1)+part%perweight
                else
                   write(*,*) j,i,charge, E(j)%numberParticles(i,charge)
                   write(*,*) particleIDs(i)
                   write(*,*) ParticleIDs
                   fname = 'SpecialEvent'
                   call event_dump(1,E,fname)
                   call Traceback()
                end if
           end if
          end do eventLoop
       end do chargeLoop
    end do idLoop

    !**************************************************************************
    ! One,Two,Three and Four particle production sigmas, e.g. pi pi, n n, ...
    ! - Here we sum over all charge channels
    !**************************************************************************

    ! Loop over the number of particles:
    do numberParticles=1,4
       idLoop1: do i=lbound(particleIDs,dim=1),ubound(particleIDs,dim=1)
          channel=channel+1
          if (channel>dimSigma) then
             call TraceBack('index value channel > bounds')
          end if

          ! (0) Create the name of the cross section channel to be stored in
          !     sigma_name
          name = PartName(particleIDs(i))
          sigma_name(channel,1)=intTochar(channel)//  &
                   &': '//intTochar(numberParticles)//trim(name)
          sigma_name(channel,2)=intTochar(channel+dimSigma)//': error '  &
             &  //intTochar(numberParticles)//trim(name)

          perweight=0.
          ! (1) Sum over all events
          eventLoop1: do j=lbound(E,dim=1),ubound(E,dim=1)

             ! Check that there are the wished number of particles in the event:
             if ( Sum(E(j)%numberParticles(i,:)).ne.numberParticles) cycle

             foundIDs=0 ! =number of found particles in a special event
             ! (2) Sum over all charges
             ! Count contributions of any particles in the event with any
             ! charge="charge" and ID="particleIDs(i)":
             ch_loop1: do charge=-2,2
                if (.not.validCharge_ID(particleIDs(i),charge)) cycle
                ! Check if charge is valid
                numInEventLoop: do num=1,numberParticles
                   ! if there are less than num particles in this event
                   ! with that charge:
                   if (.not.event_GetParticle(E(j),particleIDs(i),charge,num, &
                       &part)) exit numInEventLoop

                   foundIDs=foundIds+1
                   perweight(foundIds)=part%perweight
                   if (foundIds.eq.numberParticles) exit ch_loop1
                end do numInEventLoop
             end do ch_loop1
             if (foundIds.ne.numberParticles) then
                write(*,*) 'serious error in event_sigma', foundIds,&
                           & numberParticles
                call TraceBack()
             else
                ! (3) Add up cross section.
                ! We add the minimum of all perweights
                ! since the minimum determines the right probability
                sigma(channel,1)=sigma(channel,1)  &
                                 & +minval(perweight(1:numberParticles))
             end if
          end do eventLoop1
       end do idLoop1
    end do

    !**************************************************************************
    ! Two particle production sigmas
    ! - Special channels
    !**************************************************************************
    channel=channel+1
    if (channel>dimSigma) then
       call TraceBack('index value channel > bounds')
    end if
    sigma_name(channel,1)=intTochar(channel)//': pi N'
    sigma_name(channel,2)=intTochar(channel+dimSigma)//': error  pi N'

    ! Count contributions of any particles in the event with any charge="charge"
    ! and ID="particleIDs(i)":
    eventLoop2: do j=lbound(E,dim=1),ubound(E,dim=1)
    ! 1 particle must be a pion (i=1), another one a nucleon (i=numStableMesons+1)
       if (( Sum(E(j)%numberParticles(1,:)).eq.1) .and.  &
          &  ( Sum(E(j)%numberParticles(numStableMesons + 1,:)).eq.1)) then
          foundIDs=0
          ! Search pion
          ch_loop2: do charge=-1,1
             if ( event_GetParticle(E(j),pion,charge,1,part)) then
                foundIDs=foundIds+1
                perweight(foundIds)=part%perweight
                exit ch_loop2
             end if
          end do ch_loop2
          if (foundIds.ne.1) then
             write(*,*) 'foundIds =', foundIds
             call TraceBack()
          end if
          ! Search nucleon
          ch_loop2n: do charge=0,1
             if ( event_GetParticle(E(j),nucleon,charge,1,part)) then
                foundIDs=foundIds+1
                perweight(foundIds)=part%perweight
                exit ch_loop2n
             end if
          end do ch_loop2n
          if (foundIds.ne.2) then
             write(*,*) 'foundIds =', foundIds
             call TraceBack()
          else
             sigma(channel,1)=sigma(channel,1)+min(perweight(1),perweight(2))
          end if
       end if
    end do eventLoop2

    !**************************************************************************
    ! Special channels: only one proton, but no pions
    !**************************************************************************
    channel=channel+1
    if (channel>dimSigma) then
       call TraceBack('index value channel > bounds')
    end if
    sigma_name(channel,1)=intTochar(channel)//': p w/o pi'
    sigma_name(channel,2)=intTochar(channel+dimSigma)//': error p w/o pi'

    do j=lbound(E,dim=1),ubound(E,dim=1)
      ! make sure that there is 1 proton present:
       if (E(j)%numberParticles(numStableMesons+1,1).eq.1) then
      ! make sure that there are no pions
          if (sum(E(j)%numberParticles(1,-1:1)).ne.0) cycle
!          !no neutrons
!          if(E(j)%numberParticles(numStableMesons+1,0).ne.0) cycle
          if (event_GetParticle(E(j),nucleon,1,1,part)) then
             sigma(channel,1)=sigma(channel,1)+part%perweight
          else
             call Traceback('Error in event_sigma: p w/o pi -> Stop!')
          end if
       end if
    end do

    !**************************************************************************
    ! Special channels: muon, but no pions
    !**************************************************************************
    channel=channel+1
    if (channel>dimSigma) then
       call TraceBack('index value channel > bounds')
    end if
    sigma_name(channel,1)=intTochar(channel)//': mu w/o pi'
    sigma_name(channel,2)=intTochar(channel+dimSigma)//': error mu w/o pi'

    do j=lbound(E,dim=1),ubound(E,dim=1)
       if (sum(E(j)%numberParticles(1,-1:1)).eq.0) then
          perweight=0.
          if (event_GetParticle(E(j),outLeptonID,outLeptonCharge,1,part)) &
          &   perweight(1)=part%perweight
          !if(associated(E(j)%particleList%first))  &
          !  &perweight(1)=E(j)%particleList%first%V%perweight
          sigma(channel,1)=sigma(channel,1)+perweight(1)
       end if
    end do

    !**************************************************************************
    !
    ! One particle sigmas + X
    ! - no sum over the different charge states
    ! - If there is another particle of the same species, but with another charge,
    !   then we count the event.
    ! - E.g. for pi^+ pi^0 the event contributes both to "pi+ + X" and "pi0 + X"
    !**************************************************************************
    idLoop_single: do i=lbound(particleIDs,dim=1),ubound(particleIDs,dim=1)
       chargeLoop_single: do charge=-2,2
          if (validCharge_ID(particleIDs(i),charge)) then
          ! Check if channel is valid
             channel=channel+1
             if (channel>dimSigma) then
                call TraceBack('index value channel > bounds')
             end if

             ! (0) Create the name of the cross section channel to be stored
             ! in sigma_name
             if (isHadron(particleIds(i))) then
                name=trim(hadron(particleIDs(i))%name)
             else
                call TraceBack('severe problem in event_sigma!')
             end if
             name = PartName(particleIds(i),charge,.false.)
             sigma_name(channel,1)=intTochar(channel)//': '// trim(name)//' +X'
             sigma_name(channel,2)=intTochar(channel+dimSigma)//': error '//&
                        & trim(name)//' +X'

             ! (1) Sum over all events
             eventLoop_single: do j=lbound(E,dim=1),ubound(E,dim=1)
                ! (2) Count contributions of any particles in the event with
                ! charge="charge" and ID="particleIDs(i):
                if ( E(j)%numberParticles(i,charge).ge.1) then
                   if (event_GetParticle(E(j),particleIDs(i),charge,1,part)) then
                      sigma(channel,1)=sigma(channel,1)+part%perweight
                   else
                      write(*,*) 'numberParticles =', &
                                 & E(j)%numberParticles(i,charge)
                      call TraceBack()
                   end if
                end if
             end do eventLoop_single
          end if
       end do chargeLoop_single
    end do idLoop_single

    !**************************************************************************
    !
    ! Special pion cross section according to Krusche
    ! - count each pi^0 separately, e.g. if there are 2 they count twice!
    !**************************************************************************

    channel=channel+1
    if (channel>dimSigma) then
       call TraceBack('index value channel > bounds')
    end if

    sigma_name(channel,1)=intTochar(channel)//': Krusche pi0'
    sigma_name(channel,2)=intTochar(channel+dimSigma)//': error Krusche pi0'

    ! (1) Sum over all events
    eventLoop_Krusche1: do j=lbound(E,dim=1),ubound(E,dim=1)
       ! (2) Search for pi0 and count each pi0 seperately:
       if ( E(j)%numberParticles(1,0).ge.1) then
          if ( event_GetParticle(E(j),particleIDs(1),0,1,part)) then
             sigma(channel,1)=sigma(channel,1)   &
                 & +part%perweight* E(j)%numberParticles(1,0)
          else
             write(*,*) 'numberParticles =', E(j)%numberParticles(1,0)
             call TraceBack()
          end if
       end if
    end do eventLoop_Krusche1


    channel=channel+1
    if (channel>dimSigma) then
       call TraceBack('index value channel > bounds')
    end if

    sigma_name(channel,1)=intTochar(channel)//': Krusche pi+'
    sigma_name(channel,2)=intTochar(channel+dimSigma)//': error Krusche pi+'

    ! (1) Sum over all events
    eventLoop_Krusche2: do j=lbound(E,dim=1),ubound(E,dim=1)
       ! (2) Search for pi+ and count each pi+ seperately:
       if ( E(j)%numberParticles(1,1).ge.1) then
          if ( event_GetParticle(E(j),particleIDs(1),1,1,part)) then
             sigma(channel,1)=sigma(channel,1)+part%perweight* &
                & E(j)%numberParticles(1,1)
          else
             write(*,*) 'numberParticles =', E(j)%numberParticles(1,1)
             call TraceBack()
          end if
       end if
    end do eventLoop_Krusche2


    channel=channel+1
    if (channel>dimSigma) then
       call TraceBack('index value channel > bounds')
    end if

    sigma_name(channel,1)=intTochar(channel)//': Krusche p'
    sigma_name(channel,2)=intTochar(channel+dimSigma)//': error Krusche p'

    ! (1) Sum over all events
    eventLoop_Krusche3: do j=lbound(E,dim=1),ubound(E,dim=1)
       ! (2) Search for p and count each p seperately:
       if ( E(j)%numberParticles(numStableMesons+1,1).ge.1) then
          if ( event_GetParticle(E(j),particleIDs(numStableMesons+1),1,1,part)) &
            & then
             sigma(channel,1)=sigma(channel,1)    &
             &  +part%perweight* E(j)%numberParticles(numStableMesons+1,1)
          else
             write(*,*) 'numberParticles =', &
                  & E(j)%numberParticles(numStableMesons+1,1)
             call TraceBack()
          end if
       end if
    end do eventLoop_Krusche3


    channel=channel+1
    if (channel>dimSigma) then
       call TraceBack('index value channel > bounds')
    end if

    sigma_name(channel,1)=intTochar(channel)//': Krusche n'
    sigma_name(channel,2)=intTochar(channel+dimSigma)//': error Krusche n'

    ! (1) Sum over all events
    eventLoop_Krusche4: do j=lbound(E,dim=1),ubound(E,dim=1)
       ! (2) Search for p and count each p seperately:
       if ( E(j)%numberParticles(numStableMesons+1,0).ge.1) then
         if (event_GetParticle(E(j),particleIDs(numStableMesons+1),0,1,part)) then
             sigma(channel,1)=sigma(channel,1) + &
               &  part%perweight* E(j)%numberParticles(numStableMesons+1,0)
         else
             write(*,*) 'numberParticles =',   &
                        & E(j)%numberParticles(numStableMesons+1,0)
             call TraceBack()
         end if
       end if
    end do eventLoop_Krusche4


    channel=channel+1
    if (channel>dimSigma) then
       call TraceBack('index value channel > bounds')
    end if

    sigma_name(channel,1)=intTochar(channel)//': Krusche pi-'
    sigma_name(channel,2)=intTochar(channel+dimSigma)//': error Krusche pi-'

    ! (1) Sum over all events
    eventLoop_Krusche5: do j=lbound(E,dim=1),ubound(E,dim=1)
       ! (2) Search for pi- and count each pi- seperately:
       if ( E(j)%numberParticles(1,-1).ge.1) then
          if ( event_GetParticle(E(j),particleIDs(1),-1,1,part)) then
             sigma(channel,1)=sigma(channel,1)  &
               & + part%perweight* E(j)%numberParticles(1,-1)
          else
             write(*,*) 'numberParticles =', E(j)%numberParticles(1,-1)
             call TraceBack()
          end if
       end if
    end do eventLoop_Krusche5


    !**************************************************************************
    ! Special channel: muon, but no nucleons
    !**************************************************************************
    channel=channel+1
    if (channel>dimSigma) then
       call TraceBack('index value channel > bounds')
    end if
    sigma_name(channel,1)=intTochar(channel)//': 0N'
    sigma_name(channel,2)=intTochar(channel+dimSigma)//': error 0N'

    do j=lbound(E,dim=1),ubound(E,dim=1)
       if (sum(E(j)%numberParticles(numStableMesons + 1,0:1)).eq.0) then
          perweight=0.
          if (event_GetParticle(E(j),outLeptonID,outLeptonCharge,1,part)) &
          &  perweight(1)=part%perweight
          sigma(channel,1)=sigma(channel,1)+perweight(1)
       end if
    end do


    !**************************************************************************
    ! Special channel: 5 or more nucleons and everything
    !**************************************************************************
    channel=channel+1
    if (channel>dimSigma) then
       call TraceBack('index value channel > bounds')
    end if
    sigma_name(channel,1)=intTochar(channel)//': >4N'
    sigma_name(channel,2)=intTochar(channel+dimSigma)//': error >4N'

    do j=lbound(E,dim=1),ubound(E,dim=1)
       if (sum(E(j)%numberParticles(numStableMesons+1,0:1)).ge.5) then
          perweight=0.
           ! read perweight from the first particle in the list
          if (associated(E(j)%particleList%first)) &
          &  perweight(1)=E(j)%particleList%first%V%perweight
          sigma(channel,1)=sigma(channel,1)+perweight(1)
       end if
    end do



    !**************************************************************************
    !
    ! Output Results to files and to input arrays:
    !
    !**************************************************************************

    if (initSigma) then
       initSigma=.false.
       !***********************************************************************
       !****o* AnaEvent/sigma.dat
       ! NAME
       ! file sigma.dat
       ! PURPOSE
       ! * First column  : Naming for different runs or runNumber
       ! * Columns 2-... : Single and multiparticle cross sections.
       !                   The header of the file names the explicit channels.
       ! NOTES
       ! Cf. documenation in AnaEvent/event_sigma for the definition of "single"
       ! and "multi" particle cross sections. Esp. note the exclusive aspects.
       !
       !***********************************************************************
       open(12,file='sigma.dat')
       do j=1,channel
          write(12,'("# ",A,"    ",A)') sigma_name(j,1),sigma_name(j,2)
       end do
       formatOut='("#",A8,'//intTochar(2*channel)//'A14)'
       write(12,formatOut) sigma_name(:channel,1)(1:4),   &
             & sigma_name(:channel,2)(1:4)

    else
       open(12,file='sigma.dat',position='append')
    end if


    if (present(sigmaOut).and.present(runNumber).and.present(newInit)) then

! For error evaluation we save the square of the change in the present run
! into sigma(:,2) :
!!$       do i=1,100
!!$          sigma(i,2)=sigma(i,2)+(sigma(i,1)-sigmaOut(i,1))**2
!!$       end do
       sigma(:,2)=sigma(:,2)+(sigma(:,1)-sigmaOut(:,1))**2
       sigmaOut=sigma
       do i=1,100
          if (runNumber.gt.1) then
             ! Statististical error of the mean value of several runs:
             sigma(i,2)=sqrt(max(((sigma(i,2)-sigma(i,1)**2/float(runNumber))/ &
             &  (float((runNumber-1)*(runNumber)))),0.))
          else
             sigma(i,2)=0.
          end if
       end do
       !avoid output since the files is large and of no immediate use
       !formatOut='(A,F8.4,'//intTochar(dimSigma+channel)//'E14.5)'
       !write(12,formatOut) ' ',sigma(1,1),sigma(2:,1)/float(runNumber), &
       ! & sigma(:channel,2)

    else
       ! NO ERROR EVALUATION, SINCE NO KNOWLEDGE ABOUT EARLIER RUNS

       !avoid output since the files is large and of no immediate use
       !formatOut='(A,'//intTochar(channel)//'E14.5)'
       !write(12,formatOut) ' ',sigma(1:channel,1)
    end if
    close(12)


  end subroutine event_sigma




  !****************************************************************************
  !****s* AnaEvent/event_dSigma_dOmega
  ! NAME
  ! subroutine event_dSigma_dOmega(E, dTheta, dPhi, string, runNumber,
  ! hists_theta, hists_phi, hists2D, initHists, makeOutputIn, sameFileNameIn)
  !
  ! PURPOSE
  ! Output for single particle production dsigma/dTheta and dsigma/dTheta/dPhi
  ! with dTheta,dPhi in units of sr=radian for all particle-IDs includede in
  ! "particleIDs" (see above)  and all charge states.
  !
  ! INPUTS
  ! * type(tAnaEvent), dimension(:) :: E -- List of events
  ! * real         :: dTheta,dPhi -- Delta(Theta) and Delta(Phi) in [Radian]
  ! * character(*) :: string -- used as prefix for all output file
  ! * integer      :: runNumber -- number of the run =1,2,3,4...
  !
  ! Optional:
  ! * type(histogram) intent(inout),
  !   dimension(lBound(particleIDs,dim=1):uBound(particleIDs,dim=1),-2:2) ::
  !   hists_theta,hists_phi
  ! * type(histogram2D), intent(inout),
  !   dimension(lBound(particleIDs,dim=1):uBound(particleIDs,dim=1),-2:2) ::
  !   hists2D
  ! * logical, intent(in)  :: initHists
  ! * logical :: makeOutputIn -- Flag to switch off output (.false.=no output)
  ! * logical, intent(in) :: sameFileNameIn  --
  !   .true. = always print to same filenames,
  !   .false. = different files for different runs
  ! These optional input variables return the evaluated histograms.
  ! In hists_theta you will find the dsigma_dtheta histograms,
  ! in hists_phi you will find the dsigma_dphi histograms and
  ! in hists2D dsigma_dOmega.
  ! They must be given all together, to make an effect on the program.
  !
  ! If initHists=.true., then the input histograms "hists_theta","hists_phi"
  ! and "hists2D" are initialized again. If inithists=.false. then the input
  ! histograms are used as results of earlier runs and the results of the
  ! present run are just added. If .true. then we assume that runNumber is
  ! the number to use as a normalization for the histogram results = We
  ! assume that the histograms contain the results of "runNumber" runs.
  !
  ! The fourth column in the histogram defines the error of the 2nd column
  ! (only if optional input is given!)
  !
  ! NOTES
  ! If there are more than one particle of the same species produced, then the
  ! event is not considered for this "single particle cross section".
  !
  ! OUTPUT
  ! * filenames 'STRING_dsigma_dTheta...*.#.dat',
  !   'STRING_dsigma_dTheta_Phi...*.#.dat'
  !   where #=1,2,3,... is given by runNumber and "STRING" by the input string.
  !****************************************************************************
  subroutine event_dSigma_dOmega(E,dTheta,dPhi,string,runNumber,hists_theta,  &
                                 & hists_phi,hists2D, &
                                 &  initHists,makeOutputIn,sameFileNameIn)
    use particleDefinition
    use rotation, only:  get_phi_theta
    use particleProperties, only: PartName, validCharge_ID
    use output, only: intTochar_pm,intToChar
    use histf90
    use hist2Df90
    use constants, only: pi

    logical, intent(in), optional :: sameFileNameIn
    logical, intent(in), optional :: makeOutputIn
    real, intent(in) :: dTheta, dphi
    character(*), intent(in) :: string
    integer, intent(in) :: runNumber
    ! Optional
    type(histogram)  ,optional, intent(inout),dimension(1:numStableParts,-2:2) &
        & :: hists_theta,hists_phi
    type(histogram2D),optional, intent(inout),dimension(1:numStableParts,-2:2) &
        & :: hists2D
    logical, intent(in),optional :: initHists
    ! Input
    type(tAnaEvent) , dimension(:) :: E ! List of events
    ! Local
    type(particle) :: part
    type(histogram) ::    dsigma_dTheta,dsigma_dPhi
    type(histogram2D) ::  dsigma_dTheta_dPhi
    integer :: charge,i,j
    character(80) :: filename_theta,fileName2D,name,filename_phi
    real :: theta, phi
    logical :: makeInitHists

    logical :: makeOutput, optionals,sameFileName
    real :: fak

    ! Check optionals:
    if (present(makeOutputIn)) then
       makeOutput=makeOutputIN
    else
       makeOutput=.true.
    end if

    if (present(sameFileNameIn)) then
       sameFileName=sameFileNameIn
    else
       sameFileName=.false.
    end if


    if (present(hists_theta).neqv.present(hists_phi) &
      & .or. present(hists2D).neqv.present(initHists)  &
      & .or. present(hists_phi).neqv.present(hists2D)) then
       call TraceBack('You have to give either no or all optional arguments!! &
                      &  Stop')
    end if
    optionals = (present(hists_theta).and.present(hists_phi) &
         .and.present(hists2D).and.present(initHists))

    if (runNumber.lt.1) then
       write(*,*) 'Error in event_dSigma_dOmega. runNumber.lt.1'
       stop
    end if

    idLoop: do i=lbound(particleIDs,dim=1),ubound(particleIDs,dim=1)
       if (.not.particleIDs_flag(i)) cycle
       chargeLoop: do charge=-2,2
          if (.not.validCharge_ID(particleIDs(i),charge)) cycle

          makeInitHists=.true.
          if (optionals) then
             ! Do not initialize histograms, but use input as basis:
             if (.not.initHists) makeInitHists=.false.
          end if

          if (makeInitHists) then
             ! Initialize histogram
             name = PartName(particleIDs(i),charge,.false.)
             if (i <= numStableMesons .and. exclusive_hadron)  &
                    & name = trim('excl'//name)
             call createHist(dsigma_dTheta, 'dSigma/dTheta('// trim(name) //&
                  & ') for single '//trim(name)// ' production' &
                  &  ,0. ,pi,dTheta)
             call createHist(dsigma_dPhi, 'dSigma/dPhi('// trim(name) //&
                  & ') for single '//trim(name)// ' production' &
                  &  ,0. ,2.*pi,dPhi)
             call createHist2D(dsigma_dTheta_dPhi,              &
                          &'dSigma/dTheta/dPhi('// trim(name) //&
                  & ') for single '//trim(name)// ' production' &
                  &  ,(/0.,0./) ,(/pi,2.*pi/),(/dTheta,dPhi/))
          else
             dsigma_dTheta     = hists_theta(i,charge)
             dsigma_dPhi       = hists_phi(i,charge)
             dsigma_dTheta_dPhi= hists2D(i,charge)
          end if

          ! instead of PartName(ID)//'_charge_'//intTochar_pm(charge) I would
          ! prefer to use PartName(ID,charge,.false) in the following filenames
          ! [KG]

          name = PartName(particleIDs(i))
          if (i <= numStableMesons .and. exclusive_hadron)   &
                  & name = trim('excl'//name)
          if (sameFileName) then
             filename2D=trim(string)//'_dSigma_dTheta_dPhi_'// trim(name)&
                  & //'_charge_'//intTochar_pm(charge) //'.dat'
             filename_theta=trim(string)//'_dSigma_dTheta_'// trim(name)&
                  & //'_charge_'//intTochar_pm(charge)//'.dat'
             filename_phi=trim(string)//'_dSigma_dPhi_'// trim(name)&
                  & //'_charge_'//intTochar_pm(charge)//'.dat'
          else
             filename2D=trim(string)//'_dSigma_dTheta_dPhi_'// trim(name)&
                  & //'_charge_'//intTochar_pm(charge) //'.' &
                  & //trim(intTochar(runNumber))//'.dat'
             filename_theta=trim(string)//'_dSigma_dTheta_'// trim(name)&
                  & //'_charge_'//intTochar_pm(charge)//'.'  &
                  & //trim(intTochar(runNumber))//'.dat'
             filename_phi=trim(string)//'_dSigma_dPhi_'// trim(name)&
                  & //'_charge_'//intTochar_pm(charge)//'.'  &
                  & //trim(intTochar(runNumber))//'.dat'
          end if

          ! Count contributions of any particles in the event with charge="charge"
          ! and ID="particleIDs(i):
          eventLoop: do j=lbound(E,dim=1),ubound(E,dim=1)
             if (exclusive_hadron) then
                if ( ((i <= numStableMesons) .and.  &
                  & ( E(j)%numberParticles(i,charge).eq.1)  &
                  & .and.( sum(E(j)%numberParticles(1:numStableMesons,:)).eq.1))&
               ! so far exclusive 1 meson in final state, no other mesons
               ! now exclusive 1 baryon in final state, plus mesons
                  &          .or.                                          &
                  &  ((i > numStableMesons) .and.  &
                  &     ( E(j)%numberParticles(i,charge).eq.1) &
                  &  .and.( sum(E(j)%numberParticles(i,:)).eq.1)) ) then

                      if ( event_GetParticle(E(j),particleIDs(i),charge,1,part))&
                      &  then
                         call get_phi_Theta(part%momentum(1:3),theta,phi)
                   !                    theta=degrees(theta)
                   !                    phi  =degrees(phi)
                   call addHist(dsigma_dTheta       , theta         , &
                                & part%perweight)
                   call addHist(dsigma_dPhi         , phi           , &
                                & part%perweight)
                   call addHist2d(dsigma_dTheta_dPhi, (/theta,phi/) , &
                                & part%perweight)
                      else
                         write(*,*) 'numberParticles =', &
                                    & E(j)%numberParticles(i,charge)
                         call Traceback()
                      end if
                end if
          else
             if ((E(j)%numberParticles(i,charge).eq.1).and.  &
                & ( sum(E(j)%numberParticles(i,:)).eq.1))&
                & then
                if ( event_GetParticle(E(j),particleIDs(i),charge,1,part)) then
                   call get_phi_Theta(part%momentum(1:3),theta,phi)
                   !                    theta=degrees(theta)
                   !                    phi  =degrees(phi)
                   call addHist(dsigma_dTheta       , theta         ,  &
                               & part%perweight)
                   call addHist(dsigma_dPhi         , phi           ,  &
                               & part%perweight)
                   call addHist2d(dsigma_dTheta_dPhi, (/theta,phi/) ,  &
                               & part%perweight)
                else
                   write(*,*) 'numberParticles =',E(j)%numberParticles(i,charge)
                   call Traceback()
                end if
             end if
           end if
          end do eventLoop

          ! Output histograms
          if (optionals) then
             ! Make error calculation, therefore subtract the input
             ! from the present histogram to get the change
             if (.not.makeInitHists) then
                call makeError_hist(hists_theta(i,charge),dSigma_dTheta)
                call makeError_hist(hists_Phi(i,charge),dSigma_dPhi)
             else
                call makeError_hist(b=dSigma_dTheta)
                call makeError_hist(b=dSigma_dPhi)
             end if
             ! Save the lists:
             hists_theta(i,charge)=dsigma_dTheta
             hists_phi(i,charge)=dsigma_dPhi
             hists2D(i,charge)=dsigma_dTheta_dPhi
             ! Set yVal(:,3) to the statistical error before doing the output :
             if (.not.makeInitHists) then
                call setError_hist(dSigma_dTheta,runNumber)
                call setError_hist(dSigma_dPhi,runNumber)
             end if
          end if

          if (makeoutput) then
             fak = 1.0
             if (.not.makeInitHists) fak=1./float(runNumber)
             ! normalize to number of runs

             call WriteHist(dsigma_dTheta,file=filename_theta,mul=fak)
             call WriteHist(dsigma_dPhi,file=filename_phi,mul=fak)
             call WriteHist2D_gnuplot(dsigma_dTheta_dPhi,file=filename2D,mul=fak)
          end if

          call RemoveHist(dsigma_dTheta)
          call RemoveHist(dsigma_dPhi)
          call RemoveHist2D(dsigma_dTheta_dPhi)

       end do chargeLoop
    end do idLoop

  end subroutine event_dSigma_dOmega


  !****************************************************************************
  !****s* AnaEvent/event_dSigma_dE
  ! NAME
  ! subroutine event_dSigma_dE(E, eMin, eMax, dE, string, runNumber, hists,
  ! initHists, makeOutputIn, sameFileNameIn, histsMulti, hists1X, hists2X)
  !
  ! PURPOSE
  ! Output for single particle production dsigma_dE for all particle IDs
  ! included in the  array "particleIDs" (see above) and all charge states.
  !
  ! INPUTS
  ! * type(tAnaEvent), dimension(:) :: E -- List of events
  ! * real :: eMin, eMax, dE -- Minimal/Maximal Ekin, Delta(Ekin) in [GeV]
  ! * character(*) :: string -- used as prefix for all output files
  ! * integer      :: runNumber -- number of the run
  !
  ! Optional:
  ! * type(histogram), intent(inout),
  !   dimension(lBound(particleIDs,dim=1):uBound(particleIDs,dim=1),-2:2) ::
  !   hists
  ! * logical, intent(in) :: initHists
  ! * logical :: makeOutputIn -- Flag to switch off output (.false.=no output)
  ! * logical, intent(in) :: sameFileNameIn  --
  !   .true. = always print to same filenames,
  !   .false. = different files for different runs
  !
  ! These optional input variables return the evaluated histograms.
  ! In hists you will find the dsigma_dE.
  !
  ! If initHists=.true., then the input histograms "hists" and "histsMulti" are
  ! initialized again.
  ! If inithists=.false., then the input histograms are used as results of
  ! earlier runs and the results of the present run are just added. If .true.
  ! then we assume that runNumber is the number to use as a normalization for
  ! the histogram results = We assume that the histograms contain the results
  ! of "runNumber" runs.
  !
  ! NOTES
  ! If there are more than one particle of the same species produced, then the
  ! event is not considered for the "single particle cross section", but for the
  ! "multi particle cross section".
  !
  ! OUTPUT
  ! * filenames 'STRING_dsigma_dEkin...*.#.dat'
  !   where #=1,2,3,... is given by runNumber
  !   and "STRING" by the input string.
  !****************************************************************************
  subroutine event_dSigma_dE(E,eMin,eMax,dE,string,runNumber,hists,initHists, &
    & makeOutputIn,sameFileNameIn,histsMulti,hists1X,hists2X)

    use particleDefinition
    use particleProperties, only: validCharge_ID, PartName
    use output, only: intTochar_pm,intToChar
    use histf90


    logical , intent(in), optional ::sameFileNameIn
    ! .true. = always print to same filenames.,
    ! .false. = different files for different runs

    logical , intent(in), optional ::makeOutputIn
    ! Flag to switch off output: .false. = no output
    real, intent(in) :: eMin, eMax, dE
    character(*), intent(in) :: string ! used as prefix for all output files
    integer, intent(in) :: runNumber ! number of the run
    type(tAnaEvent) , dimension(:) :: E ! List of events
    type(particle) :: part
    type(histogram) :: dsigma_dE,dsigma_dE_multi, dsigma_dE_1X, dsigma_dE_2X
    integer :: charge,i,j, particle_number
    character(80) :: filename,filename1X,filename2X,filenameMulti,name,name1
 !   character(100) :: fname
    real :: ekin
    type(histogram),optional, intent(inout),dimension(1:numStableParts,-2:2) ::&
        & hists
    type(histogram),optional, intent(inout),dimension(1:numStableParts,-2:2) ::&
        & hists1X
    type(histogram),optional, intent(inout),dimension(1:numStableParts,-2:2) ::&
        & hists2X
    type(histogram),optional, intent(inout),dimension(1:numStableParts,-2:2) ::&
        & histsMulti
    logical, intent(in),optional :: initHists
    logical ::  makeInitHists

    logical :: makeOutput,sameFileName
    real :: fak

    if (present(makeOutputIn)) then
       makeOutput=makeOutputIN
    else
       makeOutput=.true.
    end if

    if (present(sameFileNameIn)) then
       sameFileName=sameFileNameIn
    else
       sameFileName=.false.
    end if

    if (runNumber.lt.1) then
       call Traceback('Error in event_dSigma_dE. runNumber.lt.1')
    end if

 ! first loop over all particles (mesons and baryons) in IDTABLE
    idLoop: do i=lbound(particleIDs,dim=1),ubound(particleIDs,dim=1)
       if (.not.particleIDs_flag(i)) cycle

 ! then loop over individual charge states
       chargeLoop: do charge=-2,2
          if (.not.validCharge_ID(particleIDs(i),charge)) cycle

          ! Initialize histogram
          makeInitHists=.true.
          if (present(hists).and.present(initHists)) then
             if (.not.initHists) makeInitHists=.false.
          end if

          name = PartName(particleIds(i),charge,.false.)
          name1 = name
          if (i <= numStableMesons .and. exclusive_hadron) name1 =&
             & trim('excl'//name)
          if (makeInitHists) then
             call createHist(dsigma_dE, 'dSigma/dEkin('// trim(name1) //&
                  & ') for single '//trim(name1)// ' production' &
                  &  ,eMin ,eMax,dE)
             call createHist(dsigma_dE_1X, 'dSigma/dEkin('// trim(name) //&
                  & ') for 1'//trim(name)//  &
                  & 'and X(including the same id but different charge)&
                  & production',eMin ,eMax,dE)
             call createHist(dsigma_dE_2X, 'dSigma/dEkin('// trim(name) //&
                  & ') for 2'//trim(name)//  &
                  &'and X(including the same id but different charge)&
                  & production',eMin ,eMax,dE)
             call createHist(dsigma_dE_multi, 'dSigma/dEkin('// trim(name) //&
                  & ') for multi '//trim(name)// ' production' &
                  &  ,eMin ,eMax,dE)
          else
             dsigma_dE=hists(i,charge)
             dsigma_dE_1X=hists1X(i,charge)
             dsigma_dE_2X=hists2X(i,charge)
             dsigma_dE_multi=histsMulti(i,charge)
          end if

          ! instead of PartName(ID)//'_charge_'//intTochar_pm(charge) I would
          ! prefer to use PartName(ID,charge,.false) in the following filenames
          ! [KG]

         name = PartName(particleIDs(i))
         name1 = name
         if (i <= numStableMesons .and. exclusive_hadron) &
            &name1 = trim('excl'//name)
          if (sameFileName) then
             filename=trim(string)//'_dSigma_dEkin_'// trim(name1)//   &
             & '_charge_'//  &
             &  intTochar_pm(charge) //'.dat'
             filename1X=trim(string)//'_dSigma_dEkin_'// trim(name)//  &
             & '_charge_'//  &
             &  intTochar_pm(charge) //'_1X.dat'
             filename2X=trim(string)//'_dSigma_dEkin_'// trim(name)//  &
             & '_charge_'//  &
             &  intTochar_pm(charge) //'_2X.dat'
             filenameMulti=trim(string)//'_dSigma_dEkin_'// trim(name)// &
             & '_charge_'//  &
             &  intTochar_pm(charge) //'_MULTI.dat'
          else
             filename=trim(string)//'_dSigma_dEkin_'// trim(name1)&
                  & //'_charge_'//intTochar_pm(charge) //'.' //&
                  &  trim(intTochar(runNumber))//'.dat'
             filename1X=trim(string)//'_dSigma_dEkin_'// trim(name)&
                  & //'_charge_'//intTochar_pm(charge) //'.' //&
                  &  trim(intTochar(runNumber))// &
                  &  '_1X.dat'
             filename2X=trim(string)//'_dSigma_dEkin_'// trim(name)&
                  & //'_charge_'//intTochar_pm(charge) //'.' //&
                  & trim(intTochar(runNumber))//  &
                  &  '_2X.dat'
             filenameMulti=trim(string)//'_dSigma_dEkin_'// trim(name)&
                  & //'_charge_'//intTochar_pm(charge) //'.' //&
                  & trim(intTochar(runNumber))//  &
                  &  '_MULTI.dat'
          end if


      ! Find out whether there are any particles in the event 'j'
      ! with charge="charge" and ID="particleIDs(i)

          eventLoop: do j=lbound(E,dim=1),ubound(E,dim=1)

      ! Now exclusive X-sections for hadrons.

      ! So far, events are checked for exclusivity only for 'stable particles',
      ! as defined earlier in AnaEvent im array ParticleIDs.
      ! It could thus be that an 'exclusive event' still contains in addition
      ! a particle that is not included in ParticleIDs. These are only
      ! heavy mesons beyond the D and heavy baryons, beyond the OmegaResonance.

    ! For mesons (i <= numStableMesons) there is
    ! exactly 1 meson of given flavor 'i' and charge 'charge' and no other meson

 Exclhadr:if(exclusive_hadron .eqv. .true.) then
  Exclmes: if (i <= numStableMesons) then
    Excl:   if ( E(j)%numberParticles(i,charge) .eq. 1  &
            & .and. sum(E(j)%numberParticles(1:numStableMesons,:)) .eq. 1) then
     partexist:   if ( event_GetParticle(E(j),particleIDs(i),charge,1,part)) then
                         ekin=part%momentum(0)-part%mass
                         call addHist(dsigma_dE , ekin  ,part%perweight)

                  else  partexist
                   write(*,*) 'numberParticles =',E(j)%numberParticles(i,charge)
                         call Traceback()
                      end if   partexist
                end if   Excl
            end if Exclmes

    ! Now exclusive X-Section for baryons,
    ! There is exactly 1 baryon of given flavor 'i' and charge 'charge'
    ! and no other particle

 Exclbar: if (i > numStableMesons) then
    Exclb:     if ( ( E(j)%numberParticles(i,charge).eq.1) &
                 & .and.( sum(E(j)%numberParticles(:,:)).eq.1)) then
       partex:    if ( event_GetParticle(E(j),particleIDs(i),charge,1,part)) then
                         ekin=part%momentum(0)-part%mass
                         call addHist(dsigma_dE , ekin  ,part%perweight)

                  else  partex
                   write(*,*) 'numberParticles =',E(j)%numberParticles(i,charge)
                        call Traceback()
                      end if   partex
                end if   Exclb
            end if Exclbar
       end if Exclhadr

    ! Now X-section for exactly 1 hadron of a given kind and charge,
    ! but all sorts of other  particles with different flavor and charge

       if (exclusive_hadron .eqv. .false.) then
       Incl:  if (( E(j)%numberParticles(i,charge).eq.1)  .and. &
                  &  ( sum(E(j)%numberParticles(i,:)).eq.1)) then
       Incltrue:  if ( event_GetParticle(E(j),particleIDs(i),charge,1,part)) then
                       ekin=part%momentum(0)-part%mass
                       call addHist(dsigma_dE , ekin ,part%perweight)

                  else Incltrue
                   write(*,*) 'numberParticles =',E(j)%numberParticles(i,charge)
                       call Traceback()
                    end if  Incltrue
                 end if   Incl
              end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! event selection check
!   if (charge == -1 .and. i == 1 .and. (trim(string) == 'diff_000')) then
!    fname = 'SpecialEvent'
!    call event_dump(1,E,fname)
!    write(*,111) 'SECOND i= ',i,'j= ',j,ekin,part%momentum(0),part%perweight, &
!    & E(j)%numberParticles(i,charge),&
!    & E(j)%numberParticles(1:numStableMesons,-2), &
!    & E(j)%numberParticles(1:numStableMesons,-1), &
!    & E(j)%numberParticles(1:numStableMesons,0),&
!    & E(j)%numberParticles(1:numStableMesons,+1), &
!    & E(j)%numberParticles(1:numStableMesons,+2), &
!    & sum(E(j)%numberParticles(1,:)), &
!    & sum(E(j)%numberParticles(1:numStableMesons,:))
!    111 Format(2(A,I5),3F10.3, I5/6I5/6I5/6I5/6I5/6I5/2I5)
!    stop
!   end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


             ! Loop for 1-particle-plus-X cross sections (e.g. 1pi+
             ! and other pi0, pi- in the event)
             !Find out whether there are any particles in the event
             !with charge="charge" and ID="particleIDs(i):
             if (( E(j)%numberParticles(i,charge).eq.1)) then
                if ( event_GetParticle(E(j),particleIDs(i),charge,1,part)) then
                  ekin=part%momentum(0)-part%mass
                  call addHist(dsigma_dE_1X       , ekin        ,part%perweight)
                else
                   write(*,*) 'numberParticles =',E(j)%numberParticles(i,charge)
                   call Traceback()
                end if
             end if



             ! Loop for 2-particle-plus-X cross sections
             ! (e.g. 2pi+ and other pi0, pi- in the event)
             !Find out whether there are any particles in the event
             ! with charge="charge"  and ID="particleIDs(i):
             if (E(j)%numberParticles(i,charge).eq.2) then
                do particle_number=1,2
                   if ( event_GetParticle(E(j),particleIDs(i),charge, &
                       & particle_number,part)) then
                      ekin=part%momentum(0)-part%mass
                      call addHist(dsigma_dE_2X  , ekin   ,part%perweight)
                   else
                    write(*,*) 'numberParticles =',E(j)%numberParticles(i,charge)
                      call Traceback()
                   end if
                end do
             end if



             ! for multi-particle cross sections
             ! Find out whether there are any particles in the event
             ! with charge="charge" and ID="particleIDs(i):
             if (E(j)%numberParticles(i,charge).ge.1) then
                !write(*,'(3(A,I5))') 'i=',i, '   charge=',charge,&
                ! & 'numberParticles=', &
                ! & E(j)%numberParticles(i,charge)
                do particle_number=1,E(j)%numberParticles(i,charge)
                   if ( event_GetParticle(E(j),particleIDs(i),charge,  &
                       & particle_number,part)) then
                      ekin=part%momentum(0)-part%mass
                      !write(*,'(A,I5,2(A,g15.5))') 'particle_number=', &
                      !    & particle_number, ' &
                      ! &   Ekin=',ekin, '    perweight=',part%perweight
                      call addHist(dsigma_dE_multi , ekin  ,part%perweight)
                   else
                      write(*,*) 'numberParticles =', &
                         & E(j)%numberParticles(i,charge)
                      call Traceback()
                   end if
                end do
             end if

          end do eventLoop


          if (present(hists)) then
             if (.not.makeInitHists) then
                ! Calculate the square of the changes compared to last run:
                call makeError_hist(hists(i,charge),dsigma_dE)
                ! Save the list
                hists(i,charge)=dsigma_dE
                ! Set yVal(:,3) to the statistical error before doing the output:
                call setError_hist(dsigma_dE,runNumber)
             else
                call makeError_hist(b=dsigma_dE)
                ! Save the list
                hists(i,charge)=dsigma_dE
             end if
          end if


          if (present(hists1X)) then
             if (.not.makeInitHists) then
                ! Calculate the square of the changes compared to last run:
                call makeError_hist(hists1X(i,charge),dsigma_dE_1X)
                ! Save the list
                hists1X(i,charge)=dsigma_dE_1X
                ! Set yVal(:,3) to the statistical error before doing the output:
                call setError_hist(dsigma_dE_1X,runNumber)
             else
                call makeError_hist(b=dsigma_dE_1X)
                ! Save the list
                hists1X(i,charge)=dsigma_dE_1X
             end if
          end if

          if (present(hists2X)) then
             if (.not.makeInitHists) then
                ! Calculate the square of the changes compared to last run:
                call makeError_hist(hists2X(i,charge),dsigma_dE_2X)
                ! Save the list
                hists2X(i,charge)=dsigma_dE_2X
                ! Set yVal(:,3) to the statistical error before doing the output:
                call setError_hist(dsigma_dE_2X,runNumber)
             else
                call makeError_hist(b=dsigma_dE_2X)
                ! Save the list
                hists2X(i,charge)=dsigma_dE_2X
             end if
          end if

          if (present(histsMulti)) then
             if (.not.makeInitHists) then
                ! Calculate the square of the changes compared to last run:
                call makeError_hist(histsMulti(i,charge),dsigma_dE_multi)
                ! Save the list
                histsMulti(i,charge)=dsigma_dE_multi
                ! Set yVal(:,3) to the statistical error before doing the output:
                call setError_hist(dsigma_dE_multi,runNumber)
             else
                call makeError_hist(b=dsigma_dE_multi)
                ! Save the list
                histsMulti(i,charge)=dsigma_dE_multi
             end if
          end if

          ! Output histograms
          if (makeOutput) then
             fak = 1.0
             if (.not.makeInitHists) fak=1./float(runNumber)
             ! normalize to number of runs

             call WriteHist(dsigma_dE,file=filename,mul=fak)
             call WriteHist(dsigma_dE_1X,file=filename1X,mul=fak)
             call WriteHist(dsigma_dE_2X,file=filename2X,mul=fak)
             call WriteHist(dsigma_dE_multi,file=filenameMulti,mul=fak)


          end if

          call RemoveHist(dsigma_dE)
          call RemoveHist(dsigma_dE_1X)
          call RemoveHist(dsigma_dE_2X)
          call RemoveHist(dsigma_dE_multi)

       end do chargeLoop

    end do idLoop

  end subroutine event_dSigma_dE



  !****************************************************************************
  !****s* AnaEvent/event_dSigma_dEcostheta
  ! NAME
  ! subroutine event_dSigma_dEcostheta(E, EcosthetaMin, EcosthetaMax,
  ! dEcostheta, string, runNumber, hists, hists_c, hists_MULTI, hists_c_MULTI,
  ! initHists, makeOutputIn, sameFileNameIn)
  !
  ! PURPOSE
  ! Output for single particle production  dsigma/d(E(1-cosTheta))
  ! for all particle-IDs includede in "particleIDs" (see above) and all
  ! charge states.
  ! Theta is the angle of the outgoing particle
  ! with respect to the neutrino beam direction
  !
  ! The output is also given for dsigma/dcosTheta
  !
  ! INPUTS
  ! * type(tAnaEvent), dimension(:) :: E -- List of events
  ! * real ::
  ! * character(*), intent(in) :: string -- used as prefix for all output file
  ! * integer, intent(in) :: runNumber -- number of the run =1,2,3,4...
  !
  ! Optional:
  ! * type(histogram), intent(inout),
  !   dimension(lBound(particleIDs,dim=1):uBound(particleIDs,dim=1),-2:2) ::
  !   hists
  ! * logical, intent(in) :: initHists
  ! * logical ::makeOutputIn -- Flag to switch off output (.false.=no output)
  ! * logical , intent(in), optional ::sameFileNameIn  --
  !   .true. = always print to same filenames,
  !   .false. = different files for different runs
  ! These optional input variables return the evaluated histograms.
  !
  ! If initHists=.true., then the input histograms are initialized again.
  ! If inithists=.false. then the input histograms are used as results of
  ! earlier runs and the results of the present run are just added. If .true.
  ! then we assume that runNumber is the number to use as a normalization
  ! for the histogram results = We assume that the histograms
  ! contain the results of "runNumber" runs.
  !
  ! The fourth column in the histogram defines the error of the 2nd column
  ! (only if optional input is given!)
  !
  ! NOTES
  ! If there are more than one particle of the same species produced,
  ! then the event is not considered for
  ! this "single particle cross section".
  !
  ! OUTPUT
  ! * filenames 'STRING_dsigma_dEcostheta...*.#.dat',
  !   where #=1,2,3,... is given by runNumber
  !   and "STRING" by the input string.
  !****************************************************************************
  subroutine event_dSigma_dEcostheta(E,EcosthetaMin,EcosthetaMax,dEcostheta,  &
             & string,runNumber, hists,hists_c, hists_MULTI, hists_c_MULTI,   &
             & initHists,  makeOutputIn, sameFileNameIn)
    use particleDefinition
    use IdTable, only: isHadron
    use particleProperties, only: hadron, validCharge_ID
    use output, only: intTochar_pm,intToChar
    use histf90

    logical , intent(in), optional ::sameFileNameIn
    logical , intent(in), optional ::makeOutputIn
    real, intent(in) :: EcosthetaMin,EcosthetaMax,dEcostheta
    character(*),intent(in) :: string
    integer, intent(in) :: runNumber
    ! Optional
    type(histogram)  ,optional, intent(inout),dimension(1:numStableParts,-2:2) &
                 & :: hists
    type(histogram)  ,optional, intent(inout),dimension(1:numStableParts,-2:2) &
                 & :: hists_MULTI
    type(histogram)  ,optional, intent(inout),dimension(1:numStableParts,-2:2) &
                 & :: hists_c
    type(histogram)  ,optional, intent(inout),dimension(1:numStableParts,-2:2) &
                 & :: hists_c_MULTI
    logical, intent(in),optional :: initHists
    ! Input
    type(tAnaEvent) , intent(in), dimension(:) :: E ! List of events
    ! Local
    type(particle) :: part
    type(histogram) ::    dsigma_dEcostheta, dsigma_dcostheta,   &
       &   dsigma_dEcostheta_MULTI,&
       &   dsigma_dcostheta_MULTI
    integer :: charge,i,j, particle_number
    character(80) :: name, name1, filename_Ecostheta, filename_costheta, &
       & filename_Ecostheta_MULTI, filename_costheta_MULTI
!    character(100) :: fname
    logical :: makeInitHists

    real :: Ecostheta,costheta

    logical :: makeOutput, optionals,sameFileName

    ! Check optionals:
    if (present(makeOutputIn)) then
       makeOutput=makeOutputIN
    else
       makeOutput=.true.
    end if

    if (present(sameFileNameIn)) then
       sameFileName=sameFileNameIn
    else
       sameFileName=.false.
    end if


    if (present(hists).and.present(initHists)) then
       optionals=.true.
    else
       optionals=.false.
    end if


    if (runNumber.lt.1) then
       call TraceBack('runNumber.lt.1')
    end if

    idLoop: do i=lbound(particleIDs,dim=1),ubound(particleIDs,dim=1)
       if (.not.particleIDs_flag(i)) cycle
       chargeLoop: do charge=-2,2
          if (validCharge_ID(particleIDs(i),charge)) then
             if (isHadron(particleIds(i))) then
                name=trim(hadron(particleIDs(i))%name)
                name1 = name
                if (i <= numStableMesons .and. exclusive_hadron) name1 = &
                  & trim('excl'//name)
             else
                call TraceBack()
             end if


             makeInitHists=.true.
             if (optionals) then
                if (.not.initHists) then
                   ! Do not initialze histograms, but use input as basis.
                   ! write(*,*)'  makeInitHists=.false.'
                   makeInitHists=.false.
                end if
             end if

             if (makeInitHists) then
             ! Initialize histogram
                call createHist(dsigma_dEcostheta, 'dSigma/dEcostheta('//  &
                     & trim(name1) //&
                     & ') for single '//trim(name1)// ' production' &
                     &  ,EcosthetaMin,EcosthetaMax,dEcostheta)
                call createHist(dsigma_dcostheta, 'dSigma/dcostheta('//    &
                     & trim(name1) //&
                     & ') for single '//trim(name1)// ' production' &
                     &  ,-1.0,1.0,0.02)
                call createHist(dsigma_dEcostheta_MULTI, &
                     &'dSigma/dEcostheta('// trim(name) //&
                     & ') for multiple '//trim(name)// ' production' &
                     &  ,EcosthetaMin,EcosthetaMax,dEcostheta)
                call createHist(dsigma_dcostheta_MULTI, &
                     & 'dSigma/dcostheta('// trim(name) //&
                     & ') for multiple '//trim(name)// ' production' &
                     &  ,-1.0,1.0,0.02)
             else
                dsigma_dEcostheta=hists(i,charge)
                dsigma_dcostheta=hists_c(i,charge)
                dsigma_dEcostheta_MULTI=hists_MULTI(i,charge)
                dsigma_dcostheta_MULTI=hists_c_MULTI(i,charge)
             end if

             if (sameFileName) then
                filename_Ecostheta=trim(string)//'_dSigma_dEcostheta_'//  &
                     & trim(name1)&
                     & //'_charge_'//intTochar_pm(charge) //'.dat'
                filename_costheta=trim(string)//'_dSigma_dcostheta_'//  &
                     & trim(name1)&
                     & //'_charge_'//intTochar_pm(charge) //'.dat'
                filename_Ecostheta_MULTI=trim(string)//'_dSigma_dEcostheta_'// &
                     & trim(name)&
                     & //'_charge_'//intTochar_pm(charge) //'_MULTI.dat'
                filename_costheta_MULTI=trim(string)//'_dSigma_dcostheta_'//   &
                     & trim(name)&
                     & //'_charge_'//intTochar_pm(charge) //'_MULTI.dat'
             else
                filename_Ecostheta=trim(string)//'_dSigma_dEcostheta_'// &
                     & trim(name1)&
                     & //'_charge_'//intTochar_pm(charge)//'.'   &
                     & //trim(intTochar(runNumber))//'.dat'
                filename_costheta=trim(string)//'_dSigma_dcostheta_'// &
                     & trim(name1)&
                     & //'_charge_'//intTochar_pm(charge)//'.'  &
                     & //trim(intTochar(runNumber))//'.dat'
                filename_Ecostheta_MULTI=trim(string)//'_dSigma_dEcostheta_'// &
                     & trim(name)&
                     & //'_charge_'//intTochar_pm(charge)//'.'   &
                     & //trim(intTochar(runNumber))// &
                     & '_MULTI.dat'
                filename_costheta_MULTI=trim(string)//'_dSigma_dcostheta_'//  &
                     & trim(name)&
                     & //'_charge_'//intTochar_pm(charge)//'.'   &
                     & //trim(intTochar(runNumber))// &
                     & '_MULTI.dat'
             end if

             ! Count contributions of any particles in the event
             ! with charge="charge" and ID="particleIDs(i):

             eventLoop: do j=lbound(E,dim=1),ubound(E,dim=1)
             ! single particle events

    excl: if (exclusive_hadron) then

         sel:  if ( ((i <= numStableMesons) .and.  &
                  & ( E(j)%numberParticles(i,charge).eq.1)  &
                  & .and.( sum(E(j)%numberParticles(1:numStableMesons,:)).eq.1))&
               ! so far exclusive 1 meson in final state, no other mesons
               ! now exclusive 1 baryon in final state, plus mesons
                  &          .or.                                          &
                  &  ((i > numStableMesons) .and.  &
                  & ( E(j)%numberParticles(i,charge).eq.1) &
                  &  .and.( sum(E(j)%numberParticles(i,:)).eq.1)) ) then

         getpart:  if ( event_GetParticle(E(j),particleIDs(i),charge,1,part)) then
                        if (dot_product(part%momentum(1:3),part%momentum(1:3)) &
                          & .gt.0) then
                           costheta=part%momentum(3)/   &
                          & sqrt(dot_product(part%momentum(1:3), &
                          &  part%momentum(1:3)))
                        else
                           costheta=0.
                        end if

                      Ecostheta=part%momentum(0)*(1.-costheta)
                      call addHist(dsigma_dEcostheta, Ecostheta ,part%perweight)
                      call addHist(dsigma_dcostheta, costheta ,part%perweight)

                   else  getpart
                      write(*,*) 'numberParticles =',E(j)%numberParticles(i,charge)
                      call Traceback()
                   end if  getpart
                 end if   sel

          else excl

              if ((E(j)%numberParticles(i,charge).eq.1) .and.  &
                &  ( sum(E(j)%numberParticles(i,:)).eq.1))&
                &  then
                   if ( event_GetParticle(E(j),particleIDs(i),charge,1,part)) then
                      if (dot_product(part%momentum(1:3),part%momentum(1:3))  &
                        & .gt.0) then
                         costheta=part%momentum(3)/   &
                                  & sqrt(dot_product(part%momentum(1:3), &
                                  &  part%momentum(1:3)))
                      else
                         costheta=0.
                      end if
                      Ecostheta=part%momentum(0)*(1.-costheta)
                      call addHist(dsigma_dEcostheta, Ecostheta ,part%perweight)
                      call addHist(dsigma_dcostheta, costheta ,part%perweight)
                   else
                      write(*,*) 'numberParticles =', &
                                 & E(j)%numberParticles(i,charge)
                      call Traceback()
                   end if
              end if
          end if excl


             ! multi particle events
             if (E(j)%numberParticles(i,charge).ge.1) then
               do particle_number=1,E(j)%numberParticles(i,charge)
                   if ( event_GetParticle(E(j),particleIDs(i),charge,  &
                      & particle_number,part)) then
                      if (dot_product(part%momentum(1:3),part%momentum(1:3)) &
                      & .gt.0) then
                         costheta=part%momentum(3)/  &
                         & sqrt(dot_product(part%momentum(1:3), &
                         & part%momentum(1:3)))
                      else
                         costheta=0.
                      end if
                      Ecostheta=part%momentum(0)*(1.-costheta)
                      call addHist(dsigma_dEcostheta_MULTI,   &
                                  & Ecostheta,part%perweight)
                      call addHist(dsigma_dcostheta_MULTI,    &
                                  & costheta ,part%perweight)
                   else
                      write(*,*) 'numberParticles =',   &
                                 & E(j)%numberParticles(i,charge)
                      call Traceback()
                   end if
                end do
             end if


             end do eventLoop

             ! Output histograms
             if (optionals) then
                ! Make error calculation, therefore subtract the input
                ! from the present histogram to get the change
                if (.not.makeInitHists) then
                   call makeError_hist(hists(i,charge),dSigma_dEcostheta)
                   call makeError_hist(hists_c(i,charge),dSigma_dcostheta)
                   call makeError_hist(hists_MULTI(i,charge),  &
                        & dSigma_dEcostheta_MULTI)
                   call makeError_hist(hists_c_MULTI(i,charge),  &
                        & dSigma_dcostheta_MULTI)
                else
                   call makeError_hist(b=dSigma_dEcostheta)
                   call makeError_hist(b=dSigma_dcostheta)
                   call makeError_hist(b=dSigma_dEcostheta_MULTI)
                   call makeError_hist(b=dSigma_dcostheta_MULTI)
                end if
                ! Save the lists:
                hists(i,charge)=dsigma_dEcostheta
                hists_c(i,charge)=dsigma_dcostheta
                hists_MULTI(i,charge)=dsigma_dEcostheta_MULTI
                hists_c_MULTI(i,charge)=dsigma_dcostheta_MULTI
                ! Set yVal(:,3) to the statistical error before doing the output:
                if (.not.makeInitHists) then
                   call setError_hist(dSigma_dEcostheta,runNumber)
                   call setError_hist(dSigma_dcostheta,runNumber)
                   call setError_hist(dSigma_dEcostheta_MULTI,runNumber)
                   call setError_hist(dSigma_dcostheta_MULTI,runNumber)
                end if
             end if

             if (makeoutput) then
                open(12,file=filename_Ecostheta)
                if (makeInitHists) then
                   call WriteHist(dsigma_dEcostheta,12,0.,1.)
                else
                   ! normalize to number of runs
                   call WriteHist(dsigma_dEcostheta,12,0.,1./float(runNumber))
                end if
                close(12)

                open(12,file=filename_costheta)
                if (makeInitHists) then
                   call WriteHist(dsigma_dcostheta,12,0.,1.)
                else
                   ! normalize to number of runs
                   call WriteHist(dsigma_dcostheta,12,0.,1./float(runNumber))
                end if
                close(12)

                open(12,file=filename_Ecostheta_MULTI)
                if (makeInitHists) then
                   call WriteHist(dsigma_dEcostheta_MULTI,12,0.,1.)
                else
                   ! normalize to number of runs
                   call WriteHist(dsigma_dEcostheta_MULTI,12,0.,1./  &
                        & float(runNumber))
                end if
                close(12)

                open(12,file=filename_costheta_MULTI)
                if (makeInitHists) then
                   call WriteHist(dsigma_dcostheta_MULTI,12,0.,1.)
                else
                   ! normalize to number of runs
                   call WriteHist(dsigma_dcostheta_MULTI,12,0.,1./  &
                        & float(runNumber))
                end if
                close(12)

             end if

             call RemoveHist(dsigma_dEcostheta)
             call RemoveHist(dsigma_dcostheta)
             call RemoveHist(dsigma_dEcostheta_MULTI)
             call RemoveHist(dsigma_dcostheta_MULTI)

          end if
       end do chargeLoop
    end do idLoop

  end subroutine event_dSigma_dEcostheta





  !****************************************************************************
  !****s* AnaEvent/event_dSigma_dInvMass
  ! NAME
  ! subroutine event_dSigma_dInvMass(E, string, runNumber,
  ! dW_Npi, Wmin_Npi, Wmax_Npi, histsW_nucleon_pion,
  ! dW_mupi, Wmin_mupi, Wmax_mupi, histsW_muon_pion,
  ! dW_muN, Wmin_muN, Wmax_muN, histsW_muon_nucleon,
  ! initHists, makeOutputIn, sameFileNameIn)
  !
  ! PURPOSE
  ! Output  dsigma_d(InvariantMass) for 1 pion production process
  ! for 3 invariant masses: W(pion-nucleon'), W(muon-nucleon'), W(pion-nucleon')
  ! can be compared with the corresponding ANL and BNL data
  !
  ! INPUTS
  ! * type(tAnaEvent), dimension(:) :: E -- List of events
  ! * real ::  dInvMass --  Delta(InvMass) in [GeV]
  ! * character(*) :: string -- used as prefix for all output files
  ! * integer      :: runNumber -- number of the run
  !
  ! Optional:
  ! * type(histogram), intent(inout),
  !   dimension(lBound(particleIDs,dim=1):uBound(particleIDs,dim=1),-2:2) ::
  !   hists
  ! * logical, intent(in) :: initHists
  ! * logical :: makeOutputIn -- Flag to switch off output (.false.=no output)
  ! * logical, intent(in) :: sameFileNameIn  --
  !   .true. = always print to same filenames,
  !   .false. = different files for different runs
  !
  ! These optional input variables return the evaluated histograms.
  ! In hists you will find the dsigma_dInvMass.
  !
  ! If initHists=.true., then the input histograms "hists" is initialized again.
  ! If inithists=.false., then the input histograms are used as results of
  ! earlier runs and the results of the present run are just added. If .true.
  ! then we assume that runNumber is the number to use as a normalization for
  ! the histogram results = We assume that the histograms contain the results
  ! of "runNumber" runs.
  !
  !
  !
  !
  ! OUTPUT
  ! * filenames 'STRING_dsigma_W...*.#.dat'
  !   where #=1,2,3,... is given by runNumber
  !   and "STRING" by the input string.
  !****************************************************************************
  subroutine event_dSigma_dInvMass(E,string,runNumber, &
       & dW_Npi,Wmin_Npi,Wmax_Npi,histsW_nucleon_pion, dW_mupi,Wmin_mupi,  &
       & Wmax_mupi, &
       & histsW_muon_pion,dW_muN,Wmin_muN,Wmax_muN,histsW_muon_nucleon,    &
       & initHists,makeOutputIn, &
       & sameFileNameIn)

    use particleDefinition
    use output, only: intTochar_pm,intToChar
    use histf90
    use minkowski, only: abs4Sq!, abs4
    use initNeutrino, only: get_init_namelist

    type(tAnaEvent) , intent(in), dimension(:) :: E ! List of events
    character(*), intent(in) :: string ! used as prefix for all output files
    integer, intent(in) :: runNumber ! number of the run
    type(histogram),optional, intent(inout),dimension(-1:1) :: &
              & histsW_nucleon_pion, &
              &  histsW_muon_pion, histsW_muon_nucleon
    real, intent(in)   :: dW_Npi,Wmin_Npi,Wmax_Npi,dW_mupi,Wmin_mupi,Wmax_mupi,&
       &  dW_muN,Wmin_muN,Wmax_muN
    logical, intent(in),optional :: initHists
    logical , intent(in), optional ::makeOutputIn
    ! Flag to switch off output : .false. = no output
    logical , intent(in), optional ::sameFileNameIn
    ! .true. = always print to same filenames.,
    ! .false. = different files for different runs

    type(particle) :: part_pion, part_nucleon, part_muon
    type(histogram) :: dsigma_dW, dsigma_dW_muon_pion, dsigma_dW_muon_nucleon
    integer :: charge,j, k, nucleon_charge
    character(80) :: filenameW, filenameW_muon_pion, filenameW_muon_nucleon
    real :: W2, W2_muon_pion, W2_muon_nucleon

    logical :: makeInitHists
    logical :: makeOutput,sameFileName

    integer:: outLeptonID, outLeptonCharge

    call get_init_namelist(outLepton_ID=outLeptonID, &
             & outLepton_charge=outLeptonCharge)


    if (present(makeOutputIn)) then
       makeOutput=makeOutputIN
    else
       makeOutput=.true.
    end if

    if (present(sameFileNameIn)) then
       sameFileName=sameFileNameIn
    else
       sameFileName=.false.
    end if

    if (runNumber.lt.1) then
       write(*,*) 'Error in event_dSigma_dInvMass. runNumber.lt.1'
       stop
    end if



    pionChargeLoop: do charge=-1,1

       ! Initialize histogram
       makeInitHists=.true.

       if (present(histsW_nucleon_pion) .and. present(histsW_muon_pion) .and. &
            &  present(histsW_muon_nucleon) .and. present(initHists)) then
          if (.not.initHists) then
             ! Do not initialze histograms, but use input as basis.
             !                   write(*,*)'  makeInitHists=.false.'
             makeInitHists=.false.
          end if
       end if

       if (makeInitHists) then
          call createHist(dsigma_dW, 'dSigma/dW for 1 pion production', &
               &  Wmin_Npi, Wmax_Npi, dW_Npi)
          call createHist(dsigma_dW_muon_pion, 'dSigma/dW_muon_pion for 1 pion &
               & production', &
               &  Wmin_mupi, Wmax_mupi, dW_mupi)
          call createHist(dsigma_dW_muon_nucleon, 'dSigma/dW_muon_nucleon for  &
               & 1 pion production', &
               &  Wmin_muN, Wmax_muN, dW_muN)

          !          write(*,'(A,I5)') 'pion charge=', charge
          !          write(*,'(9(A,g10.3))') 'Wmin_Npi=',Wmin_Npi, &
          !          & '   Wmax_Npi=',Wmax_Npi,&
          !          & '  dW_Npi=',dW_Npi, '   Wmin_Npi=',Wmin_mupi, &
          !          & '   Wmax_Npi=',Wmax_mupi,  '&
          !          &  dW_Npi=',dW_mupi,'   Wmin_Npi=',Wmin_muN,    &
          !          & '   Wmax_Npi=',Wmax_muN, &
          !          & '  dW_Npi=',dW_muN

       else
          dsigma_dW=histsW_nucleon_pion(charge)
          dsigma_dW_muon_pion=histsW_muon_pion(charge)
          dsigma_dW_muon_nucleon=histsW_muon_nucleon(charge)
          ! this is the charge of pion, not nucleon
       end if

       if (sameFileName) then
          filenameW=trim(string)//'_dSigma_dW_nucleon_pion'&
               & //'_charge_'//intTochar_pm(charge) //'.dat'
          filenameW_muon_pion=trim(string)//'_dSigma_dW_muon_pion'&
               & //'_charge_'//intTochar_pm(charge) //'.dat'
          filenameW_muon_nucleon=trim(string)//'_dSigma_dW_muon_nucleon'&
               & //'_charge_'//intTochar_pm(charge) //'.dat'
       else
          filenameW=trim(string)//'_dSigma_dW_nucleon_pion'&
               & //'_charge_'//intTochar_pm(charge)//'.'  &
               & //trim(intTochar(runNumber)) //'.dat'
          filenameW_muon_pion=trim(string)//'_dSigma_dW_muon_pion'&
               & //'_charge_'//intTochar_pm(charge)//'.'  &
               & //trim(intTochar(runNumber)) //'.dat'
          filenameW_muon_nucleon=trim(string)//'_dSigma_dW_muon_nucleon'&
               & //'_charge_'//intTochar_pm(charge)//'.'  &
               & //trim(intTochar(runNumber)) //'.dat'
       end if


       !Find out whether there are pairs
       ! (nucleon-pion, muon-pion, muon-nucleon) in the event:
       eventLoop: do j=lbound(E,dim=1),ubound(E,dim=1)
          ! search for event with 1 pion  and 1 nucleon
          ! (and 1 muon automatically)
          if ((          E(j)%numberParticles(1,charge).eq.1) .and. &
               ! in numberParticles pion=1
               & sum(E(j)%numberParticles(1,:)).eq.1      .and. &
               & sum(E(j)%numberParticles(numStableMesons+1,:)).eq.1       ) then
                 ! in numberParticles nucleon=numStableMesons + 1

             !write(*,'(A,I5)') 'Event number j=',j

             ! ---------------------------------
             if ( event_GetParticle(E(j),101,charge,1,part_pion)) then
                ! ID of pion is 101

                !write(*,'(A,4g12.5,A,g12.5)') 'Found pion with 4-momentum',   &
                !      & part_pion%momentum, &
                !                         & '  and abs4(4-momentum) is',       &
                !      & abs4(part_pion%momentum)
                !write(*,'(A,g12.5)') 'pion perweight=', part_pion%perweight


                nucleon_charge=-1
                do k=0,1
                   if (  event_GetParticle(E(j),1,k,1,part_nucleon) ) then
                   ! ID of nucleon is 1
                      nucleon_charge=k
                      !write(*,'(A,4g12.5,A,g12.5)') 'Found nucleon with 4-momentum',
                      !              & part_nucleon%momentum, &
                      !              & '  and abs4(4-momentum) is',   &
                      !              & abs4(part_nucleon%momentum), &
                      !              & '  and charge is', nucleon_charge
                      !         write(*,'(A,g12.5)') 'nucleon perweight=',   &
                      !              &part_nucleon%perweight
                      exit
                   end if
                end do

                if (nucleon_charge.ge.0 .and. event_GetParticle(E(j),  &
                    & outLeptonID,outLeptonCharge, 1,part_muon) ) then

                   !write(*,'(A,4g12.5,A,g12.5)') 'Found lepton with 4-momentum', &
                   !   &  part_muon%momentum, &
                   !   & '  and abs4(4-momentum) is', abs4(part_muon%momentum)
                   !write(*,'(A,g12.5)') 'muon perweight=', part_muon%perweight

                   W2=abs4Sq(part_pion%momentum+part_nucleon%momentum)
                   W2_muon_pion=abs4Sq(part_muon%momentum+part_pion%momentum)
                   W2_muon_nucleon=abs4Sq(part_muon%momentum+part_nucleon%momentum)

                else
                   call TraceBack('Cannot get muon from the event')
                end if


             else
                call TraceBack('Cannot get pion from the event')
             end if
             ! ----------------------------------

             ! if all invariant masses positive, add to histograms
             if (W2.gt.0 .and. W2_muon_pion.gt.0 .and. W2_muon_nucleon.gt.0) then
                call addHist(dsigma_dW              , sqrt(W2)              ,  &
                            & part_muon%perweight)
                call addHist(dsigma_dW_muon_pion    , sqrt(W2_muon_pion)    ,  &
                            & part_muon%perweight)
                call addHist(dsigma_dW_muon_nucleon , sqrt(W2_muon_nucleon) ,  &
                            & part_muon%perweight)
                !write(*,'(4(A,g12.5))') 'W2=',  W2, '   W2_muon_pion=',   &
                !                        & W2_muon_pion, &
                !                        & '   W2_muon_nucleon=', W2_muon_nucleon
                !write(*,*) ''
             else
                write(*,'(4(A,g12.5))') 'One of invariant masses &
                     & is negative W2=', &
                     &  W2, '   W2_muon_pion=', W2_muon_pion, &
                     & '   W2_muon_nucleon=', W2_muon_nucleon, '!  Stop!'
                call TraceBack()
             end if

          end if
          ! end of search for event with 1 pion  and 1 nucleon
          ! (and 1 muon automatically)

       end do eventLoop


       if (present(histsW_nucleon_pion) .and. present(histsW_muon_pion)  &
          & .and. present(histsW_muon_nucleon)) then
          if (.not.makeInitHists) then
             ! Calculate the square of the changes compared to last run:
             call makeError_hist(histsW_nucleon_pion(charge),dsigma_dW)
             call makeError_hist(histsW_muon_pion(charge),dsigma_dW_muon_pion)
             call makeError_hist(histsW_muon_nucleon(charge),  &
                                & dsigma_dW_muon_nucleon)
             ! Save the list
             histsW_nucleon_pion(charge)=dsigma_dW
             histsW_muon_pion(charge)=dsigma_dW_muon_pion
             histsW_muon_nucleon(charge)=dsigma_dW_muon_nucleon
             ! Set yVal(:,3) to the statistical error before doing the output :
             call setError_hist(dsigma_dW,runNumber)
             call setError_hist(dsigma_dW_muon_pion,runNumber)
             call setError_hist(dsigma_dW_muon_nucleon,runNumber)
          else
             call makeError_hist(b=dsigma_dW)
             call makeError_hist(b=dsigma_dW_muon_pion)
             call makeError_hist(b=dsigma_dW_muon_nucleon)
             ! Save the list
             histsW_nucleon_pion(charge)=dsigma_dW
             histsW_muon_pion(charge)=dsigma_dW_muon_pion
             histsW_muon_nucleon(charge)=dsigma_dW_muon_nucleon
          end if
       end if


       ! Output histogram
       if (makeOutput) then

          open(12,file=filenameW)
          if (makeInitHists) then
             write(12,*)'# Note: charge in the name of this file is pion charge'
             call WriteHist(dsigma_dW,12,0.,1.)
          else
             ! normalize to number of runs
             call WriteHist(dsigma_dW,12,0.,1./float(runNumber))
          end if
          close(12)

          open(12,file=filenameW_muon_pion)
          if (makeInitHists) then
             call WriteHist(dsigma_dW_muon_pion,12,0.,1.)
          else
             ! normalize to number of runs
             call WriteHist(dsigma_dW_muon_pion,12,0.,1./float(runNumber))
          end if
          close(12)

          open(12,file=filenameW_muon_nucleon)
          if (makeInitHists) then
             write(12,*)'# Note! charge in the name of this file is pion charge'
             call WriteHist(dsigma_dW_muon_nucleon,12,0.,1.)
          else
             ! normalize to number of runs
             call WriteHist(dsigma_dW_muon_nucleon,12,0.,1./float(runNumber))
          end if
          close(12)

       end if


       call RemoveHist(dsigma_dW)
       call RemoveHist(dsigma_dW_muon_pion)
       call RemoveHist(dsigma_dW_muon_nucleon)

    end do pionChargeLoop

  end subroutine event_dSigma_dInvMass

                                                                                  
                                                                                        




  !****************************************************************************
  !****s* AnaEvent/event_dSigma_dLeptonVariables
  ! NAME
  ! subroutine event_dSigma_dLeptonVariables(E, eMin, eMax, dE, cost_min,
  ! cost_max, delta_cost,string,
  ! runNumber, hist_Enu, hist_E, hist_cos, hist_Q2, hist_Q2p, hist_dEdcost,
  ! initHists, makeOutputIn, sameFileNameIn, specificEvent)
  !
  ! PURPOSE
  ! Output dsi/dEkin, dsi/dcostheta for outgoing lepton and dsi/dQ^2
  !
  ! INPUTS
  ! * type(tAnaEvent), dimension(:) :: E -- List of events
  ! * real :: eMin, eMax, dE -- Minimal/Maximal Ekin, Delta(Ekin) in [GeV]
  ! * real :: cost_min, cost_max, delta_cost -- limits for lepton scatt angle
  ! * character(*) :: string -- used as prefix for all output files
  ! * integer      :: runNumber -- number of the run
  !
  ! Optional:
  ! * type(histogram), intent(inout),
  !   dimension(lBound(particleIDs,dim=1):uBound(particleIDs,dim=1),-2:2) ::
  !   hist_Enu, hist_E,  hist_cos, hist_Q2, hist_Q2p
  ! * type(histogram2D),intent(inout) :: hist_dEdcost
  ! * logical, intent(in) :: initHists
  ! * logical :: makeOutputIn -- Flag to switch off output (.false.=no output)
  ! * logical, intent(in) :: sameFileNameIn  --
  !   .true. = always print to same filenames,
  !   .false. = different files for different runs
  ! * integer, intent(in) :: specificEvent  -- specify final state
  !   (see subroutine SpecificEvent_Name)
  !
  ! These optional input variables return the evaluated histograms.
  ! In hist_E you will find the dsigma_dE(kinetic energy of the outgoing lepton).
  ! In hist_cos you will find the dsigma_dcostheta(angle of the outgoing lepton
  ! relative to incoming neutrino).
  ! In hist_Q2 you will find the dsigma_dQ^2.
  ! In hist_Q2p you will find the dsigma_dQ_p^2, Q_p^2 is calculated from kinetic
  ! energy of leading proton
  ! The Q_p^2 distribution can thus also be used to get the 0pi kinetic energy 
  ! distribution of outgoing protons
  ! In hist_dEdcost you will find the double-differential
  ! d2sigma/(dEkin dcostheta) for outgoing lepton in a 0-pion event
  !
  ! If initHists=.true., then the input histograms "hists" and "histsMulti" are
  ! initialized again.
  ! If inithists=.false., then the input histograms are used as results of
  ! earlier runs and the results of the present run are just added. If .true.
  ! then we assume that runNumber is the number to use as a normalization for
  ! the histogram results = We assume that the histograms contain the results
  ! of "runNumber" runs.
  !
  !
  ! OUTPUT
  ! * filenames 'STRING_dsigma_dEkin_lepton.#.dat'
  ! * filenames 'STRING_dsigma_dcos_lepton.#.dat'
  ! * filenames 'STRING_dsigma_dQ2_lepton.#.dat'
  ! * filenames 'STRING_dsigma_dQ2p_lepton.#.dat'
  !   where #=1,2,3,... is given by runNumber
  !   and "STRING" by the input string.
  !****************************************************************************
  subroutine event_dSigma_dLeptonVariables(E,eMin,eMax,dE,cost_min,cost_max,  &
     & delta_cost,string,runNumber,hist_Enu,hist_E, &
     & hist_cos,hist_Q2, hist_Q2p, hist_dEdcost,initHists, makeOutputIn, &
     & sameFileNameIn, specificEvent)

    use output, only: intToChar
    use neutrinoProdInfo, only: neutrinoProdInfo_Get
    use minkowski, only: abs4Sq !,abs4
    use histf90
    use hist2Df90 
    
    use particleDefinition

    type(tAnaEvent) , dimension(:) :: E ! List of events  
    type(particle) :: part
    real, intent(in) :: eMin, eMax, dE, cost_min, cost_max, delta_cost
    character(*), intent(in) :: string ! used as prefix for all output files
    integer, intent(in) :: runNumber ! number of the run
    type(histogram),optional, intent(inout) :: hist_Enu
    type(histogram),optional, intent(inout) :: hist_E
    type(histogram),optional, intent(inout) :: hist_cos
    type(histogram),optional, intent(inout) :: hist_Q2
    type(histogram2D),optional, intent(inout) :: hist_dEdcost
    
    type(histogram),optional, intent(inout) :: hist_Q2p
    
    logical, intent(in), optional :: initHists
    logical, intent(in), optional :: makeOutputIn
    logical, intent(in), optional :: sameFileNameIn
    integer, intent(in), optional :: specificEvent

    type(histogram) :: dsigma_dEnu, dsigma_dE,dsigma_dcos, dsigma_dQ2,dsigma_dQ2p
    type(histogram2D) :: dSigma_dEdCos
    integer :: j, prod_id
    character(80) :: filename_Enu, filename_E,filename_cos,filename_Q2,  &
                     & filename_Q2p,filename_Ecos
    character(13) :: name
    real :: enu,  ekin, costheta, Q2, Q2prot,perweight
    real, dimension(0:3) :: lepIn_mom, lep_mom, boson_mom,nuc_mom
    real :: nuc_mom_p,epsB,mNeut,mProt,Tp,Eprot
    integer :: Chrg_Nuc,particle_number
    logical ::  makeInitHists, pass
   

    logical :: makeOutput,sameFileName

    if (present(makeOutputIn)) then
       makeOutput=makeOutputIN
    else
       makeOutput=.true.
    end if

    if (present(sameFileNameIn)) then
       sameFileName=sameFileNameIn
    else
       sameFileName=.false.
    end if

    if (runNumber.lt.1) then
       call Traceback('runNumber.lt.1')
    end if


    if (present(specificEvent)) then
       call SpecificEvent_Name(specificEvent,name)
    else
       name=''
    end if

    ! Initialize histogram

    makeInitHists=.true.
    if (present(hist_E).and.present(initHists)) then
       if (.not.initHists) then
          ! Do not initialize histograms, but use input as basis.
          !                   write(*,*)'  makeInitHists=.false.'
          makeInitHists=.false.
       end if
    end if

    if (makeInitHists) then
       call createHist(dsigma_dEnu,  &
          & 'dSigma/dEnu( incoming neutrino ) for events'//trim(name) , &
          &  eMin ,eMax,dE)
       call createHist(dsigma_dE,    &
          & 'dSigma/dEkin( outgoing lepton ) for events'//trim(name) , &
          &  eMin ,eMax,dE)
       call createHist(dsigma_dcos,  &
          & 'dSigma/dcostheta( outgoing  lepton ) for events'//trim(name),&
          &  cost_min, cost_max, delta_cost)
       ! Q2bin and Q2max are taken numerically the same as those for energy
       call createHist(dsigma_dQ2,   &
          & 'dSigma/dQ^2( outgoing  lepton ) for events'//trim(name), &
          &  0.0, eMax, dE) 
        call createHist(dsigma_dQ2p,   &
          & 'dSigma/dQ^2( outgoing  lepton ) for events with 0pi, 1p'//trim(name), &
          &  0.0, eMax, dE)   
       call CreateHist2D(dSigma_dEdCos, &
              & 'd2sigma vs Ekin and cost (lepton) for 0pi events '//trim(name), &
              & (/0.,cost_min/),(/eMax,cost_max/),(/dE,delta_cost/))
    else
       dsigma_dEnu=hist_Enu
       dsigma_dE=hist_E
       dsigma_dcos=hist_cos
       dsigma_dQ2=hist_Q2
       dsigma_dQ2p=hist_Q2p
       dsigma_dEdcos = hist_dEdcost
    end if

    if (sameFileName) then
       filename_Enu=trim(string)//'_dSigma_dEnu_lepton_'//trim(name)//'.dat'
       filename_E=trim(string)//'_dSigma_dEkin_lepton_'//trim(name)//'.dat'
       filename_cos=trim(string)//'_dSigma_dcos_lepton_'//trim(name)//'.dat'
       filename_Q2=trim(string)//'_dSigma_dQ2_lepton_'//trim(name)//'.dat'
       filename_Q2p=trim(string)//'_dSigma_dQ2p_lepton_'//trim(name)//'.dat'
       filename_Ecos=trim(string)//'_d2Sigma_dEkindcost_lepton_'//&
             & trim(name)//'.dat'
    else
       filename_Enu=trim(string)//'_dSigma_dEnu_lepton_'//trim(name)//'.'// &
          &  trim(intTochar(runNumber))//'.dat'
       filename_E=trim(string)//'_dSigma_dEkin_lepton_'//trim(name)//'.'//  &
          &  trim(intTochar(runNumber))//'.dat'
       filename_cos=trim(string)//'_dSigma_dcos_lepton_'//trim(name)//'.'// &
          &  trim(intTochar(runNumber))//'.dat'
       filename_Q2=trim(string)//'_dSigma_dQ2_lepton_'//trim(name)//'.'// &
          &  trim(intTochar(runNumber))//'.dat'
        filename_Q2p=trim(string)//'_dSigma_dQ2p_lepton_'//trim(name)//'.'// &
          &  trim(intTochar(runNumber))//'.dat'             
       filename_Ecos=trim(string)//'_d2Sigma_dEd_cost_lepton'// &
          &  trim(name)//'.'// trim(intTochar(runNumber))//'.dat'
    end if


    eventLoop: do j=lbound(E,dim=1),ubound(E,dim=1)

       if (present(specificEvent)) then
          pass=IfPass_SpecificEvent(specificEvent,E(j))
       else
          pass=.true.
       end if


ifpass:if (pass) then
          if (.not.neutrinoProdInfo_Get(j,prod_id,perweight,lepIn_mom,lep_mom,  &
             & boson_mom,nuc_mom,Chrg_Nuc)) then
             write(*,*) j,prod_id,perweight
             call Traceback()
          end if

          enu=lep_mom(0)+boson_mom(0)
          call addHist(dsigma_dEnu, enu  ,perweight)

          ekin=lep_mom(0)-sqrt( max(0.,abs4Sq(lep_mom)) )
          call addHist(dsigma_dE, ekin  ,perweight)

          if (dot_product(lep_mom(1:3),lep_mom(1:3)).gt.0) then
             costheta=lep_mom(3)/sqrt(dot_product(lep_mom(1:3),lep_mom(1:3)))
          else
             write(*,*) 'In event_dSigma_dLeptonVariables strange |k2_prime|=', &
                &  dot_product(lep_mom(1:3),lep_mom(1:3))
             costheta=-2.
          end if

          call addHist(dsigma_dcos, costheta, perweight)


! if(specificEvent==1)  write(*,*) 'ekin= ',ekin,'costheta= ',costheta, &
!   &'perweight= ', perweight



! Now double differential cross section for 0 pion events
!
        if (specificEvent==1) then
          call AddHist2D(dSigma_dEdCos,&
               & (/ekin,costheta/),&
               & perweight)
        end if
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Now switch analysis for MINERvAQE-like with at least 1 proton with
!!  p > nuc_mom_pcut, 0pi events, Q2 calculated from proton kinematics
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
 Q2proton: if(Q2p) then         
            
            if (specificEvent==1) then
            
    
    ! now determine energy of leading proton        
               if (E(j)%numberParticles(numStableMesons+1,1) .ge.1) then
                 Eprot = 0.0
                   do particle_number=1,E(j)%numberParticles(numStableMesons+1,1)      
                     if(event_GetParticle(E(j),1,1,particle_number,part)) then 
                        If(part%momentum(0) > Eprot) Eprot = part%momentum(0)
                     end if     
                   end do   
                end if               
                                    
               epsB = 0.034        ! phenomenol. binding energy in MINERvA analysis
               mNeut = 0.939565
               mProt = 0.938272
               Tp = Eprot - mProt
               nuc_mom_p = sqrt(Eprot**2 - mProt**2) 
             
      ! Q2 calculated from outgoing proton kinematics       
                Q2prot = (mNeut - epsB)**2 - mProt**2  &
                & + 2*(mNeut - epsB)*(Tp + mProt - mNeut + epsB)
                 if(nuc_mom_p > nuc_mom_pcut) then
                   call addHist(dsigma_dQ2p, Q2prot, perweight) 
                 end if     
               end if 
           end if  Q2proton                   

!!  Now Q2 calculated from exchanged boson          
        
          if ( -abs4Sq(boson_mom) .gt. 0) then 
           ! Q2 calculated from incoming boson 4-momentum    
             Q2=-abs4Sq(boson_mom)
             call addHist(dsigma_dQ2, Q2, perweight)  
          else
             write(*,*) 'In event_dSigma_dLeptonVariables strange Q2=', &
                  & -abs4Sq(boson_mom)
             Q2=-1.
          end if
                    
     end if ifpass
       
    end do eventLoop


    hist_dEdcost = dSigma_dEdCos

    if (present(hist_Enu)) then
       if (.not.makeInitHists) then
          ! Calculate the square of the changes compared to last run:
          call makeError_hist(hist_Enu,dsigma_dEnu)
          ! Save the list
          hist_Enu=dsigma_dEnu
          ! Set yVal(:,3) to the statistical error before doing the output :
          call setError_hist(dsigma_dEnu,runNumber)
       else
          call makeError_hist(b=dsigma_dEnu)
          ! Save the list
          hist_Enu=dsigma_dEnu
       end if
    end if


    if (present(hist_E)) then
       if (.not.makeInitHists) then
          ! Calculate the square of the changes compared to last run:
          call makeError_hist(hist_E,dsigma_dE)
          ! Save the list
          hist_E=dsigma_dE
          ! Set yVal(:,3) to the statistical error before doing the output :
          call setError_hist(dsigma_dE,runNumber)
       else
          call makeError_hist(b=dsigma_dE)
          ! Save the list
          hist_E=dsigma_dE
       end if
    end if


    if (present(hist_cos)) then
       if (.not.makeInitHists) then
          ! Calculate the square of the changes compared to last run:
          call makeError_hist(hist_cos,dsigma_dcos)
          ! Save the list
          hist_cos=dsigma_dcos
          ! Set yVal(:,3) to the statistical error before doing the output :
          call setError_hist(dsigma_dcos,runNumber)
       else
          call makeError_hist(b=dsigma_dcos)
          ! Save the list
          hist_cos=dsigma_dcos
       end if
    end if

    if (present(hist_Q2)) then
       if (.not.makeInitHists) then
          ! Calculate the square of the changes compared to last run:
          call makeError_hist(hist_Q2,dsigma_dQ2)
          ! Save the list
          hist_Q2=dsigma_dQ2
          ! Set yVal(:,3) to the statistical error before doing the output :
          call setError_hist(dsigma_dQ2,runNumber)
       else
          call makeError_hist(b=dsigma_dQ2)
          ! Save the list
          hist_Q2=dsigma_dQ2
       end if
    end if   
    
    if (present(hist_Q2p)) then
       if (.not.makeInitHists) then
          ! Calculate the square of the changes compared to last run:
          call makeError_hist(hist_Q2p,dsigma_dQ2p)
          ! Save the list
          hist_Q2p=dsigma_dQ2p
          ! Set yVal(:,3) to the statistical error before doing the output :
          call setError_hist(dsigma_dQ2p,runNumber)
       else
          call makeError_hist(b=dsigma_dQ2p)
          ! Save the list
          hist_Q2p=dsigma_dQ2p
       end if
    end if


    ! Output histogram
    if (makeOutput) then
       open(12,file=filename_Enu)
       if (makeInitHists) then
          call WriteHist(dsigma_dEnu,12,0.,1.)
       else
          ! normalize to number of runs
          call WriteHist(dsigma_dEnu,12,0.,1./float(runNumber))
       end if
       close(12)

       open(12,file=filename_E)
       if (makeInitHists) then
          call WriteHist(dsigma_dE,12,0.,1.)
       else
          ! normalize to number of runs
          call WriteHist(dsigma_dE,12,0.,1./float(runNumber))
       end if
       close(12)

       open(12,file=filename_cos)
       if (makeInitHists) then
          call WriteHist(dsigma_dcos,12,0.,1.)
       else
          ! normalize to number of runs
          call WriteHist(dsigma_dcos,12,0.,1./float(runNumber))
       end if
       close(12)

    if (specificEvent==1) then
       open(12,file=filename_Ecos)
       if (makeInitHists) then
          call WriteHist2D_Gnuplot(dsigma_dEdcos,12,0.,1.)
       else
          ! normalize to number of runs
          call WriteHist2D_Gnuplot(dsigma_dEdcos,12,0.,1./float(runNumber))
       end if
       close(12)
    end if

       open(12,file=filename_Q2)
       if (makeInitHists) then
          call WriteHist(dsigma_dQ2,12,0.,1.)
       else
          ! normalize to number of runs
          call WriteHist(dsigma_dQ2,12,0.,1./float(runNumber))
       end if
       close(12)  
       
       
       if (Q2p .and. specificEvent==1) then
       open(12,file=filename_Q2p)
         if (makeInitHists) then
            call WriteHist(dsigma_dQ2p,12,0.,1.)
         else
            ! normalize to number of runs
            call WriteHist(dsigma_dQ2p,12,0.,1./float(runNumber))
         end if 
       end if
       close(12)
       
    end if


    call RemoveHist(dsigma_dE)
    call RemoveHist(dsigma_dcos)
    call RemoveHist(dsigma_dQ2) 
    call RemoveHist(dsigma_dQ2p)
    call RemoveHist2D(dsigma_dEdcos)

  end subroutine event_dSigma_dLeptonVariables


!   !*************************************************************************
!   !****s* AnaEvent/event_Print
!   ! NAME
!   ! subroutine event_Print(L)
!   !
!   ! PURPOSE
!   ! Output the event.
!   !
!   ! INPUTS
!   ! * type(tAnaEvent) :: L -- The event
!   ! * integer     :: iFile -- The number of the buffer to output into
!   !
!   !*************************************************************************
!   subroutine event_Print(L,iFile)
!
!     use particlePointerListDefinition
!
!     type(tAnaEvent),intent(in) :: L
!     integer,intent(in)      :: iFile
!     type(tParticleListNode),Pointer  :: pNode
!     character(20),parameter :: out='(I8,4E15.6)'
!
!     write(iFile,'(12I8)') L%numberParticles
!
!     pNode => L%particleList%first
!     do
!        if (.not. ASSOCIATED(pNode)) return
!        write(iFile,out) pNode%V%ID,pNode%V%perweight,pNode%V%position
!        pNode => pNode%next
!     end do
!
!   end subroutine event_Print



  !****************************************************************************
  !****s* AnaEvent/event_dump
  ! NAME
  ! subroutine event_dump(run,L,filename,onlyFreeParticles,writeNeutrinoProd)
  !
  ! PURPOSE
  ! Output a list of events.
  !
  ! INPUTS
  ! * type(tAnaEvent),dimension(:), intent(in) :: L -- The eventlist
  ! * character(100),intent(in)      :: FileName -- The file name to write to
  ! * integer,intent(in)             :: run -- run number
  ! * logical, optional              :: onlyFreeParticles -- if .true. then
  !   only "free particles" are dumped.
  !   Here "free" means that p(0)^2-m_free^2 > 0. So far this is only
  !   implemented for nucleons.
  ! * logical, optional              :: writeNeutrinoProd -- switch to write
  !   production info
  !
  !****************************************************************************
  subroutine event_dump(run,L,fileName,onlyFreeParticles,writeNeutrinoProdID)
    use particlePointerListDefinition
    use IDTAble, only: nucleon
    use neutrinoProdInfo, only: neutrinoProdInfo_Get
    use constants, only: mN

    integer,intent(in) :: run
    type(tAnaEvent),dimension(:), intent(in) :: L
    character(100),intent(in)      :: FileName
    logical, optional              :: onlyFreeParticles
    logical, optional              :: writeNeutrinoProdID
    type(tParticleListNode),Pointer  :: pNode
    character(30),parameter :: out='(4I7,8E14.6,I11,I4,E14.6)'
    character(30),parameter :: outShort='(4I7,8E14.6,I11,I4)'
    character(35),parameter :: outA='("#",A6,A9,A6,A9,8A16,A15,A18,A8)'
    integer :: i, prod_ID, Chrg_Nuc
    real :: perweight, enu
    real, dimension(0:3) :: lepIn_mom, lepton_mom, boson_mom, nuc_mom
    logical :: newline,onlyFree,printParticle,ex,NeutrinoProdID

    if (present(onlyFreeParticles)) then
       onlyFree=onlyFreeParticles
    else
       onlyFree=.false.
    end if

    if (present(writeNeutrinoProdID)) then
       NeutrinoProdID=writeNeutrinoProdID
    else
       NeutrinoProdID=.false.
    end if


    Inquire(file=filename,exist=ex)

    if (run>1 .and. ex) then
      open(47, file=filename, status='old', position='append')
    else
      open(47, file=filename, status='replace')
      write(47,outA) '1:Run','2:Event','3:ID','4:Charge','5:perweight',  &
                     & '6:position(1)', &
                     & '7:position(2)','8:position(3)',&
                     & '9:momentum(0)', '10:momentum(1)','11:momentum(2)',  &
                     & '12:momentum(3)', &
                     & '13:history', '14:production_ID', &
                     & '15:enu'
    end if

    newline=.false.
    eventLoop:  do i=lbound(L,dim=1),ubound(L,dim=1)
       ! write(47,'(12I8)') L%numberParticles
       if (newline) write(47,*)
       newline=.false.
       pNode => L(i)%particleList%first
       do
          if (.not. associated(pNode)) cycle eventLoop
          if (onlyFree.and.pNode%V%ID.eq.nucleon) then
             ! Check that nucleon may become free,
             ! i.e. that it has in vacuum a well-defined momentum,
             ! i.e. p(0)^2-m_n^2 > 0
             printParticle=(pNode%V%momentum(0)**2-mN**2.gt.0.)
          else
             printParticle=.true.
          end if

          if (NeutrinoProdID) then
             if (.not.neutrinoProdInfo_Get(i,prod_id,perweight,lepIn_mom,  &
                & lepton_mom,boson_mom,nuc_mom,Chrg_Nuc)) then
                write(*,*) 'error in getting production info in event_dump, stop'
                stop
             end if
             enu=lepton_mom(0)+boson_mom(0)
             if (printParticle) write(47,out) run,i,pNode%V%ID,pNode%V%Charge,  &
                  &  pNode%V%perweight, &
                  &  pNode%V%position,pNode%V%momentum &
                  & ,pNode%V%history,prod_id, enu
          else
             if (printParticle) write(47,outShort) run,i,pNode%V%ID,pNode%V%Charge, &
                  &  pNode%V%perweight,pNode%V%position &
                  & ,pNode%V%momentum ,pNode%V%history
          end if
          pNode => pNode%next
          !newline=.true.
       end do
    end do eventLoop

    close(47)
  end subroutine event_dump

  !****************************************************************************
  !****s* AnaEvent/event_INIT
  ! NAME
  ! subroutine event_INIT(L)
  !
  ! PURPOSE
  ! Set the multiplicity counters back to zero.
  !
  ! INPUTS
  ! * type(tAnaEvent) :: L -- The event
  !
  ! OUTPUT
  ! (none)
  !****************************************************************************
  subroutine event_INIT(L)
    use particlePointerList, only: ParticleList_INIT

    type(tAnaEvent) :: L
    integer :: charge

    do charge=lbound(L%numberParticles,dim=2), ubound(L%numberParticles,dim=2)
       L%numberParticles(:,charge) =0
    end do

    call ParticleList_INIT(L%particleList)
  end subroutine event_INIT


  !****************************************************************************
  !****s* AnaEvent/event_CLEAR
  ! NAME
  ! subroutine event_CLEAR(L)
  !
  ! PURPOSE
  ! Reset the event: Delete the particleList and clear its memory.
  ! Set the multiplicity counters back to zero.
  !
  ! INPUTS
  ! * type(tAnaEvent) :: L -- The event
  !
  ! OUTPUT
  ! (none)
  !****************************************************************************
  subroutine event_CLEAR(L)
    use particlePointerList, only: ParticleList_CLEAR

    type(tAnaEvent) :: L
    integer :: charge

    do charge=lbound(L%numberParticles,dim=2), ubound(L%numberParticles,dim=2)
       L%numberParticles(:,charge) =0
    end do

    call ParticleList_CLEAR(L%particleList)
  end subroutine event_CLEAR


  !****************************************************************************
  !****s* AnaEvent/event_add
  ! NAME
  ! subroutine event_add(L,V)
  !
  ! PURPOSE
  ! Adds the particle (which V points at) to the event.
  !
  !
  ! INPUTS
  ! * type(tAnaEvent)            :: L -- The event
  ! * type(particle), POINTER :: V -- The particle to add
  !
  ! OUTPUT
  ! (none)
  !****************************************************************************
  subroutine event_add(L, V)
    use particleDefinition
    use particlePointerList, only: ParticleList_APPEND
    use output, only: writeParticle_debug
    use IdTable, only: isMeson

    type(tAnaEvent) :: L
    type(particle), POINTER :: V
    integer :: i

    call ParticleList_APPEND(L%particleList, V)

    if (V%perweight < 1e-20) return

    ! Count particle multiplicities:
    if (.not.V%antiparticle) then
       select case (V%ID)
 !       case(particleIDs(1):)
        case (pion,eta,kaon,kaonBar,dMeson,dBar,ds_plus,ds_minus,nucleon,  &
            & lambda,sigmaResonance,Xi,OmegaResonance)
 !
 ! NOTE: the hadron names appearing in this case have to be the same as defined in
 !       the array particleIDs defined early in the module AnaEventDefinition
 !
          do i=lbound(particleIDs,dim=1),ubound(particleIDs,dim=1)
             if (particleIDs(i).eq.V%ID) exit
          end do
          L%numberParticles(i,V%charge)=L%numberParticles(i,V%charge)+1
       case default
          L%numberParticles(numStableParts+1,V%charge)=    &
             & L%numberParticles(numStableParts+1, &
             &  V%charge)+1
       end select
    else if (isMeson(V%ID)) then
       write(*,*) 'Error in event_add: We expect no ANTI-Mesons!!! &
                  & Convert these first!! Stop!!'
       call writeParticle_debug(V)
       stop
    end if
  end subroutine event_add


  !****************************************************************************
  !****s* AnaEvent/event_getParticle
  ! NAME
  ! function event_getParticle(E,ID,Charge,n,P,antiparticle) RESULT (success)
  !
  ! PURPOSE
  ! Search in the particle list of the event "E" the "n"-th particle in the
  ! list with %ID="ID". This particle "P" is then returned.
  !
  ! INPUTS
  ! * type(tAnaEvent) :: E      -- The event
  ! * integer      :: ID     -- ID of particle which shall be returned
  ! * integer      :: Charge -- charge of particle which shall be returned
  ! * integer      :: n      -- We return the n-th particle in the list
  !   with %ID=ID
  ! * logical, optional :: antiparticle -- Whether the particle should be
  !   an antiparticle, if not specified, then we search for a particle.
  !
  ! OUTPUT
  ! * type(particle) ::  P      -- n-th particle with wished ID
  ! * logical        :: success -- True if it was possible to find n-Particles
  !   with the wished ID, False otherwise
  !****************************************************************************
  function event_getParticle(E,ID,Charge,n,P,antiparticle) RESULT (success)
    use particleDefinition
    use particlePointerList, only: ParticleList_getParticle
    integer,intent(in) :: ID
    integer,intent(in) :: charge
    integer,intent(in) :: n
    type(tAnaEvent),intent(in) :: E
    type(particle), intent(out) ::  P
    logical, intent(in),optional :: antiparticle

    logical :: success

    !type(tParticleListNode), POINTER :: pNode
    !integer :: foundIDs

    ! Default return values
    success=.false.
    call SetToDefault(p)

    ! Search particle:
    if (present(antiparticle)) then
       success=ParticleList_getParticle(E%particleList,ID,charge,n,P,   &
               & antiparticle,.true.)
    else
       success=ParticleList_getParticle(E%particleList,ID,charge,n,P,   &
               & .false.,.true.)
    end if
  end function event_getParticle


  !****************************************************************************
  !****s* AnaEvent/makeError_hist
  ! NAME
  ! subroutine makeError_hist(a,b)
  !
  ! PURPOSE
  ! Set b%y(:,3)=a%y(:,3)+( b%y(:,1)- a%y(:,1))**2
  !
  ! Mainly useful for error calculations. Usually performed to the histogram
  ! after a run with a being the histogram before the run and b being the
  ! histogram after the run.
  !
  ! INPUTS
  ! * type(histogram), optional      ::  a   -- Histogram of last run
  ! * type(histogram), intent(inOut) ::  b   -- Histogram now
  !
  ! If a is missing then we assume that the histogram was 0 at the last run
  ! OUTPUT
  ! * type(histogram), intent(inOut) ::  b   -- Histogram now
  !****************************************************************************
  subroutine makeError_hist(a,b)
    use histf90

    type(histogram), intent(in) ,optional    ::  a    ! Histogram of last run
    type(histogram), intent(inOut) ::  b    ! Histogram now
    integer :: i

    do i=lbound(b%yval,dim=1),ubound(b%yval,dim=1)
       if (present(a)) then
          b%yval(i,3)=b%yval(i,3)+(b%yval(i,1)-a%yval(i,1))**2
       else
          b%yval(i,3)=b%yval(i,3)+b%yval(i,1)**2
       end if
    end do

  end subroutine makeError_hist


  !****************************************************************************
  !****s* AnaEvent/setError_hist
  ! NAME
  ! subroutine setError_hist(a,numRuns)
  !
  ! PURPOSE
  ! Assume that a%yval(:,1) contains the following : sum over all
  ! runs(xSection of run).
  ! And assume that a%yval(:,3) contains the following : sum over all
  ! runs(xSection of run**2). The number of runs is given by numRuns.
  !
  ! Sets a%yval(:,3) to the statistical error of a%yval(:,1) times the
  ! number of runs (This is done because a%yval(:,1) is the mean value times
  ! the number of runs, too.
  !
  ! INPUTS
  ! * type(histogram), intent(inout)    ::  a  --
  ! * integer,intent(in) :: numRuns --
  !
  ! OUTPUT
  ! * type(histogram), intent(inOut) ::  a --
  !****************************************************************************
  subroutine setError_hist(a,numRuns)
    use histf90

    type(histogram), intent(inout)    ::  a
    integer,intent(in) :: numRuns
    integer :: i
    real :: dummy
    real :: epsilon=0.0001 ! Allowed rounding error

    if (numRuns.le.1) then
       a%yval(:,3)=0
    else
       do i=lbound(a%yval,dim=1),ubound(a%yval,dim=1)
          ! Check that sum(x_i^2)>(sum(x_i))^2:
          if (a%yval(i,3).lt.(a%yval(i,1)**2)/float(numRuns)) then
             write(*,'(A,2E20.10)') 'Error in setError_hist',a%yval(i,3),  &
                  & (a%yval(i,1)**2)&
                  & /float(numRuns)
             ! If there is an error check whether it is due to rounding errors:
             if (abs(a%yval(i,3)).gt.epsilon) then
                if (abs(a%yval(i,3)-(a%yval(i,1)**2)/float(numRuns))  &
                    & /abs(a%yval(i,3)).gt.epsilon) then
                   call traceback()
                end if
             end if
             a%yval(i,3)=-999999.
          else
             dummy=float(numRuns)*sqrt(max(1./float(numRuns-1)/float(numRuns)  &
                   & *(a%yval(i,3)-(a%yval(i,1)**2)/float(numRuns)),0.))
             a%yval(i,3)=dummy
          end if
       end do
    end if
  end subroutine setError_hist


  !****************************************************************************
  !****s* AnaEvent/event_DumpNumbers
  ! NAME
  ! subroutine event_DumpNumbers
  !
  ! PURPOSE
  ! ...
  !
  ! INPUTS
  ! * type(tAnaEvent) :: E       ! The event
  !
  ! OUTPUT
  ! * to stdout
  !****************************************************************************
!   subroutine event_DumpNumbers(E)
!     type(tAnaEvent), intent(in) :: E
!
!     integer :: i!,j
!     !integer, dimension(-2,2) :: NN
!
!     !NN = 0
!
!     do i=1,12
!        write(*,'(i3,":",5i5," !",i5)') i,E%numberparticles(i,:),   &
!             & sum(E%numberparticles(i,:))
!     enddo
!     write(*,'("---:-------------------------------")')
!     write(*,*)
!
!   end subroutine event_DumpNumbers



  !****************************************************************************
  !****s* AnaEvent/event_pairPhotons
  ! NAME
  ! subroutine event_pairPhotons(L,filename,dsigma_dm_threeGammas,   &
  !                              & Hist_massPhotons)
  !
  ! PURPOSE
  ! Outputs number of photons in one event and a histogram
  ! which gives the mass of all possible pairs of photons versus the mass
  ! of all possible triples.
  ! Additionally a cross section dsigma/dm_{3\gammas} is generated.
  !
  ! This is meant to reproduce an analysis procedure by the TAPS group.
  !
  !
  !
  ! INPUTS
  ! * type(tAnaEvent),dimension(:) :: L -- The eventlist
  ! * character(100)            :: fileName_data,fileName_statistics
  !   -- The file name to write to
  !
  !****************************************************************************
  subroutine event_pairPhotons(L,fileName,dsigma_dm_threeGammas,Hist_massPhotons)
    use particlePointerListDefinition
    use hist2Df90
    use histf90
    use output, only: intToChar
    use idTable, only: photon
    use inputGeneral, only: num_Runs_sameEnergy, num_Energies

    type(tAnaEvent),dimension(:), intent(in) :: L
    character(*),intent(in) :: fileName
    type(histogram2D),intent(inout) :: Hist_massPhotons
    type(histogram), dimension(0:8),intent(inout) :: dsigma_dm_threeGammas
    character(100) :: fileName_data,fileName_statistics,filename_dsdm
    ! file name to write to

    type(tParticleListNode),Pointer       :: pNode
    integer :: i,indexPhoton,numPhotons
    logical :: atLeastOneParticle
    real :: default_perweight
    real, dimension(0:20,1:2) :: list_numPhotons
    real, dimension(1:20,0:3) :: momentum
    real, dimension(1:20)   :: perweight
    real :: mass_all
    real, dimension(1:3) :: mass_pairs
    integer :: index, massLoop_i, massLoop_j, massLoop_k
    integer, save :: numRuns

    if (.not.allocated(Hist_massPhotons%yVal)) then
       call CreateHist2D(Hist_massPhotons, 'Mass(2 gamma) vs. Mass(3 gamma)',  &
          & (/0.,0./), (/1.0,1.0/),(/0.03,0.03/),.true.)
       do i=lbound(dsigma_dm_threeGammas,dim=1),  &
            & ubound(dsigma_dm_threeGammas,dim=1)
          call CreateHist  (dsigma_dm_threeGammas(i), 'dsigma/dm (3 Gammas)',  &
                           & 0.,1.,0.01)
       end do
       numRuns = num_Runs_sameEnergy * num_Energies
    end if

    filename_data =       "Photons_m2_vs_m3"//fileName//".dat"
    filename_statistics = "Photons_m2_vs_m3"//fileName//".txt"
    filename_dsdm =       "Photons_dsigma_dm"//fileName//".dat"

    open(47, file=filename_data)
    open(48, file=filename_statistics)


    write(48,'(A)') '# This statistic gives the number of produced photons in&
     &each event'
    write(48,'(A)') '# Note: Events with no final state particles,   &
         & i.e. Pauli blocked events, ' &
         & //' are not considered in this statistic!!'
    write(48,'(A4,A8,A30)')  '#','Event','number of Photons'

    list_numPhotons=0

    eventLoop:  do i=lbound(L,dim=1),ubound(L,dim=1)
       ! write(47,'(12I8)') L%numberParticles

       pNode => L(i)%particleList%first

       ! (1) Collect informations about photons in event "i"
       momentum=0.
       indexPhoton=1
       numPhotons=0
       perweight=0
       atLeastOneParticle=.false.
       default_perweight=0.
       particleLoop: do
          if (.not. associated(pNode)) exit particleLoop
          atLeastOneParticle=.true.
          default_perweight=pNode%V%perweight
          if (pNode%V%ID.eq.photon) then
             if (indexPhoton.eq.1) perweight=0
             ! Initialize for the case
             ! that there is a photon
             if (indexPhoton.le.ubound(momentum,dim=1)) then
                momentum(indexPhoton,0:3)=pNode%V%momentum
                perweight(indexPhoton)=pNode%V%perweight
                indexPhoton=indexPhoton+1
                numPhotons=numPhotons+1
             else
                write(*,*) 'WARNING: More than ', ubound(momentum,dim=1),   &
                           & ' photons in event #', i
             end if
          end if
          pNode => pNode%next
       end do particleLoop
       ! (2) Create statistics
       !     (a) Protocol number of photons and their perweight
       if (atLeastOneParticle) then
          write(48,'(A4,I8,I8)') '#',i, numPhotons
          ! Number of photons
          list_numPhotons(numPhotons,1) =list_numPhotons(numPhotons,1)+1.
          ! Sum of perweights
          if (numPhotons.ne.0) then
             list_numPhotons(numPhotons,2) =list_numPhotons(numPhotons,2) &
                & +minval(perweight(1:numPhotons))
          else
             list_numPhotons(numPhotons,2) =list_numPhotons(numPhotons,2)  &
                  & +default_perweight
          end if
       end if
       ! (b) Create mass pairings
       ! Loop over all possible triples of photons
       if (numPhotons.ge.3) then
          do massLoop_i=1,numPhotons-2
             do massLoop_j=massLoop_i+1,numPhotons-1
                do massLoop_k=massLoop_j+1,numPhotons
                   !write(*,*) massLoop_i,massLoop_j,massLoop_k, &
                   ! &perweight((/massLoop_i, &
                   ! massLoop_j, massLoop_k/))
                   ! Evaluate the mass of the triple and the masses
                   ! of all possible pairs:
                   call getMasses( momentum((/massLoop_i, massLoop_j, &
                      &  massLoop_k/),0:3), &
                      & mass_all,mass_pairs)
                   ! Protocol result
                   call AddHist(dsigma_dm_threegammas(0),mass_all,   &
                      & minval(perweight &
                      & ((/massLoop_i, massLoop_j, massLoop_k/)))/numRuns)
                   if (numPhotons.le.ubound(dsigma_dm_threegammas,dim=1)) then
                      call AddHist(dsigma_dm_threegammas(numPhotons),mass_all, &
                         & minval(perweight((/massLoop_i, massLoop_j,&
                           &  massLoop_k/)))/numRuns)
                   end if
                   do index=1,3
                      call AddHist2D(hist_massPhotons,    &
                           & (/mass_all,mass_pairs(index) /),&
                           & minval(perweight((/massLoop_i, massLoop_j, &
                           & massLoop_k/)))/numRuns)
                      if (abs(minval(perweight((/massLoop_i, massLoop_j, &
                           &massLoop_k/)))) &
                           & .lt.0.0001) then
                         write(*,*) 'Error in eventLoop!!',perweight
                      end if
                   end do
                end do
             end do
          end do
       end if
    end do eventLoop

    ! Write out statisitics:
    write(48,'(A,A14,2A15)') '#','Number of gammas','Frequency','Total weight'
    do i=lbound(list_numPhotons,dim=1),ubound(list_numPhotons,dim=1)
       write(48,'(2I15,E15.6)') i, nint(list_numPhotons(i,1)),list_numPhotons(i,2)
    end do

    call WriteHist2D_Gnuplot(hist_massPhotons,47)
    !call RemoveHist2D(hist_massPhotons)
    ! Write out pairings
    close(47)
    close(48)

    ! Write dsigma/dm
    do i=lbound(dsigma_dm_threeGammas,dim=1),ubound(dsigma_dm_threeGammas,dim=1)
       select case (i)
       case (0)
          open(100, file=trim(filename_dsdm))
       case (3:)
          open(100, file=trim(filename_dsdm)//'.'//intToChar(i))
       end select
       call WriteHist(dsigma_dm_threegammas(i),100)
       !call removeHist(dsigma_dm_threegammas(i))
       close(100)
    end do

  contains

    subroutine getMasses( p, m_3, m_2 )
      use minkowski, only: abs4
      real, dimension(1:3,0:3),intent(in)  :: p
      real                    ,intent(out) :: m_3
      real, dimension (1:3)   ,intent(out) :: m_2

      m_3   =abs4(p(1,:)+p(2,:)+p(3,:))
      m_2(1)=abs4(p(1,:)+p(2,:))
      m_2(2)=abs4(p(1,:)+p(3,:))
      m_2(3)=abs4(p(2,:)+p(3,:))
    end subroutine getMasses

  end subroutine event_pairPhotons


  !****************************************************************************
  !****s* AnaEvent/event_hadronicEnergy
  ! NAME
  ! subroutine event_hadronicEnergy(L,filename,Ehad_versusNu, dSigdNu, dSigdEhad)
  !
  ! PURPOSE
  ! Outputs histograms:
  ! * Ehad_versusNu (2D) ---  transferred energy (col#1);
  !   hadronic energy(col#2) xsec(col#3)
  ! * dSigdNu  --- dsi/dnu
  ! * dSigdEhad --- dsi/dEhad
  ! * Enurestored_versusEnu --- neutrino energy(col#1);
  !   restored neutrino energy (col#2) xsec(col#3)
  ! * dSigdEnu  --- dsi/dEnu
  ! * dSigdEnurestored  --- dsi/dEnurestored
  !
  ! Additionally info "prod_id  Enu  Enu_restored  nu  Ehad  perweight" for
  ! all events is written into the file
  !
  ! This is meant to reproduce the calorimetric analysis by the MINOS group.
  !
  ! NOTES
  ! The results are averaged over num_Energies (and, of course, also over num_Runs_sameEnergy)
  !
  ! INPUTS
  ! * type(tAnaEvent),dimension(:) :: L -- The eventlist
  ! * real               :: numin, numax, nubin  --  min, max and bin for the
  !                                nu,Ehad-related  histograms
  ! * real               :: Enumin, Enumax, Enubin  --  min, max and bin for the
  !                                Enu,Enurestored-related histograms
  !
  ! OUTPUT
  ! * type(histogram),dimension(0:5) ::  Ehad_versusNu, dSigdNu, dSigdEhad
  ! * type(histogram),dimension(0:5) ::  Enurestored_versusEnu, dSigdEnu,
  ! * dSigdEnurestored
  !   where  0=all events, 1=QE, 2=Delta, 3=highRES, 4=bgr, 5=DIS
  !****************************************************************************
  subroutine event_hadronicEnergy(L, numin, numax, nubin, Ehad_versusNu,   &
                                           & dSigdNu, dSigdEhad, &
                                           & Enumin, Enumax, Enubin,   &
                                           & Enurestored_versusEnu,&
                                           & dSigdEnu, dSigdEnurestored)
    use particlePointerListDefinition
    use histf90
    use hist2Df90
    use output, only: intToChar
    use inputGeneral, only: num_Runs_sameEnergy, num_Energies
    use IdTable, only: isBaryon, isMeson
    use neutrinoProdInfo, only: neutrinoProdInfo_Get
    use output, only: WriteParticle_debug
    use constants, only: mN
    use initNeutrino, only: K2Hist, max_Hist, includeHist

    type(tAnaEvent),dimension(:), intent(in) :: L
    real, intent(in) :: numin, numax, nubin ! boundaries for the hystograms
                                            ! in the namelist    &
                                            ! "nl_calorimetric_analysis"
                                            ! in neutrinoAnalysis.f90
    real, intent(in) :: Enumin, Enumax, Enubin ! boundaries for the hystograms
                                               ! in the namelist
                                               ! "nl_calorimetric_analysis"
    type(histogram2D),dimension(0:max_Hist), intent(inout) :: Ehad_versusNu,  &
                      & Enurestored_versusEnu
    type(histogram), dimension(0:max_Hist),intent(inout) :: dSigdNu, dSigdEhad,&
       & dSigdEnu, &
       &  dSigdEnurestored
    character(100) :: filename_Ehad_versusNu, filename_dSigdNu, filename_dSigdEhad, &
                    & filename_Enurestored_versusEnu, filename_dSigdEnu, &
                    &  filename_dSigdEnurestored ! file names to write to
    type(tParticleListNode),Pointer       :: pNode
    integer :: prod_id
    real    :: Ehad, perweight
    real, dimension(0:3) :: lepIn_mom, lep_mom, boson_mom, nuc_mom
    integer :: Chrg_Nuc
    integer :: i, iHist
    integer, save :: numRuns




    if (.not.allocated(Ehad_versusNu(0)%yVal)) then
       do iHist=0, max_Hist ! 0=all 1=QE 2=Delta 3=highRES 4=1pibgr
                            ! 5=DIS 6=2p2hQE 7=2p2hDelta &
                            ! 8=2pibgr
          if (.not.includeHist(iHist)) cycle
          call CreateHist2D(Ehad_versusNu(iHist), 'Nu -- Ehad plane',  &
             & (/numin,numin/),  &
             &  (/numax,(numax+1.)/),(/nubin,nubin/),.true.)
          call CreateHist  (dSigdNu(iHist), 'dSigdNu', numin,numax,nubin)
          call CreateHist  (dSigdEhad(iHist), 'dSigdEhad', numin,(numax+1.),  &
               &  nubin)
          call CreateHist2D(Enurestored_versusEnu(iHist),&
               & 'Enu -- Enu_restored plane', &
               & (/Enumin,Enumin/),(/Enumax,Enumax/),(/Enubin,Enubin/),.true.)
          call CreateHist  (dSigdEnu(iHist), 'dSig/dEnu', Enumin,Enumax,Enubin)
          call CreateHist  (dSigdEnurestored(iHist), 'dSig/dEnu_restored',  &
               & Enumin,Enumax,Enubin)
       end do
       numRuns = num_Runs_sameEnergy * num_Energies
    end if

    open(60,file="calorim_events_analysis.dat")
    write(60,*) '# 1:prod_id 2:enu 3:enu-restored 4:nu=energy-transfer &
                &  5:Ehad 6:perweight'
    close(60)


    filename_Ehad_versusNu = "calorim_Ehad_versusNu."
    filename_dSigdNu =    "calorim_dSigdNu."
    filename_dSigdEhad =  "calorim_dSigdEhad."
    filename_Enurestored_versusEnu = "calorim_Enurec_versusEnu."
    filename_dSigdEnu = "calorim_dSigdEnu."
    filename_dSigdEnurestored ="calorim_dSigdEnurec."

    ! write info about the event in file
    open(60,file="calorim_events_analysis.dat")
    write(60,'(A)') '1:prod_id  2:Enu 3: enu_restored 4:nu  5:Ehad 6:perweight'

    eventLoop:  do i=lbound(L,dim=1),ubound(L,dim=1) ! "i" - first event

       pNode => L(i)%particleList%first

       ! (1) Collect informations about boson in event "i"
       ! event:  prod_ID (1 to 34) , perweight, outgoing lepton
       ! and intermediate boson momenta
       if (.not.neutrinoProdInfo_Get(i,prod_id,perweight,lepIn_mom,lep_mom,  &
         & boson_mom,nuc_mom,Chrg_Nuc)) then
          write(*,*) 'error in getting production info, stop'
          stop
       end if
       if (debug) write(*,'(2(A,I5),A,g14.5)') 'first_event=',i,  &
          & '      prod_ID=',prod_ID, &
          &  '     nu=', boson_mom(0)
       Ehad=0.

       ! loop over outgoing particles, summing up hadron energy
       particleLoop: do
          if (.not. associated(pNode)) exit particleLoop

          if (debug) call WriteParticle_debug(pNode%V)

          if ( isBaryon(pNode%V%ID) ) then
             ! for baryons add  energy-nucleon mass,
             ! for antybaryons energy+nucleon mass
             if (pNode%V%antiparticle) then
                Ehad=Ehad+pNode%V%momentum(0)+mN !antibaryon energy
                                                 ! + compensation for an
                                                 ! accompaning baryon
             else
                Ehad=Ehad+pNode%V%momentum(0)-mN
             end if
          else if ( isMeson(pNode%V%ID) ) then
             ! for mesons add energy
             Ehad=Ehad+pNode%V%momentum(0)
             if (debug) write(*,*) 'Ehad=', Ehad
          end if

          pNode => pNode%next
       end do particleLoop
       if (debug) write(*,*) ' ------------------------------------ '
       write(60,'(I5,5g14.5)') prod_id, boson_mom(0)+lep_mom(0),   &
          &  Ehad+lep_mom(0), boson_mom(0), &
          &  Ehad, perweight

       ! (2) Create statistics

       iHist=K2Hist(prod_id)

       call AddHist2D( Ehad_versusNu(0),  Ehad_versusNu(iHist),  &
            & (/boson_mom(0), Ehad/), perweight)
       call AddHist( dSigdNu(0), dSigdNu(iHist),  boson_mom(0),  perweight)
       call AddHist( dSigdEhad(0), dSigdEhad(iHist), Ehad,          perweight)

       call AddHist2D( Enurestored_versusEnu(0), Enurestored_versusEnu(iHist), &
                &  (/boson_mom(0)+lep_mom(0), Ehad+lep_mom(0)/), perweight)
       call AddHist( dSigdEnu(0), dSigdEnu(iHist),  boson_mom(0)+lep_mom(0),  &
                & perweight)
       call AddHist( dSigdEnurestored(0), dSigdEnurestored(iHist),   &
            & Ehad+lep_mom(0),    perweight)

    end do eventLoop
    close(60)

    ! (3) Write out statisitics:

    do iHist=0,max_Hist
       if (.not.includeHist(iHist)) cycle

       call WriteHist2D_Gnuplot(Ehad_versusNu(iHist),&
            & mul=1.0/numRuns,&
            & file=trim(filename_Ehad_versusNu)//trim(intToChar(iHist))//'.dat',&
            & dump=.false.)
       call WriteHist(dSigdNu(iHist),&
            & mul=1.0/numRuns,&
            & file=trim(filename_dSigdNu)//trim(intToChar(iHist))//'.dat',&
            & dump=.false.)
       call WriteHist(dSigdEhad(iHist),&
            & mul=1.0/numRuns,&
            & file=trim(filename_dSigdEhad)//trim(intToChar(iHist))//'.dat',&
            & dump=.false.)

       call WriteHist2D_Gnuplot(Enurestored_versusEnu(iHist),&
            & mul=1.0/numRuns,&
            & file=trim(filename_Enurestored_versusEnu)//trim(intToChar(iHist))&
            & //'.dat',&
            & dump=.false.)
       call WriteHist(dSigdEnu(iHist),&
            & mul=1.0/numRuns,&
            & file=trim(filename_dSigdEnu)//trim(intToChar(iHist))//'.dat',&
            & dump=.false.)
       call WriteHist(dSigdEnurestored(iHist), &
            & mul=1.0/numRuns,&
            & file=trim(filename_dSigdEnurestored)//trim(intToChar(iHist))//&
            & '.dat',&
            & dump=.false.)
    end do

  end subroutine event_hadronicEnergy



  subroutine SpecificEvent_Name(specificEvent,name)
    integer, intent(in) :: specificEvent
    character(13), intent(out)  :: name
    select case (specificEvent)
    case (1)
       name = 'no_pi'
    case (2)
       name = 'p_Xn_no_pi'
    case (3)
       name = 'piplus'
    case (4)
       name = 'piplus_MULTI'
    case (5)
       name = 'pi0'
    case (6)
       name = 'pi0_MULTI'
    case (7)
       name = 'piminus'
    case (8)
       name = 'piminus_MULTI'
    case (9)
       name = 'pp_no_pi'
    case (10)
       name = 'pn_no_pi'
    case (11)
       name = 'nn_no_pi'
    case (12)
       name = 'pp_Xn_no_pi'
    case (13)
       name = 'nn_Xp_no_pi'
    case (14)
       name = 'ppp_Xn_no_pi'
    case (15)
       name = 'pppp_Xn_no_pi'
    case (16)
       name = 'p_no_pi'
    case (17)
       name = 'n_no_pi'
    case (18)
       name = 'Xn_no_pi'
    case (19)
       name = 'excl_pi0'
    case (20)
       name = 'excl_piplus'
    case (21)
       name = 'excl_piminus'  
    case default
       write(*,*) 'In SpecificEvent_Name WRONG specificEvent=', specificEvent, &
                 & 'STOP'
       stop
    end select
  end subroutine SpecificEvent_Name

  function IfPass_SpecificEvent(specificEvent,event) result(r)
    logical :: r
    integer, intent(in) :: specificEvent ! predefined type of event  
    type(tAnaEvent), intent(in) :: event ! particular event under consideration
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! In the following the first index in numberParticles refers to the particleId
! given in the vector 'particleIds' which is defined in the 
! module AnaEventDefinition. These particleIds refer only to output particles
! in a detaileed analysis and are *not* identical to the GiBUU particle ids.
!integer, dimension(1:numStableParts),parameter, public :: particleIDs=(/ &
!       & pion, eta, kaon, kaonBar, DMeson, dBar,ds_plus,ds_minus,&
!       & nucleon, lambda, sigmaResonance, Xi, OmegaResonance/)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    select case (specificEvent)
    case (1) ! 0 pions
       r = ( sum(event%numberParticles(1,:)).eq.0 .and. &
            & sum(event%numberParticles(numStableMesons+1,:)) >0 )
    case (2) ! 1 proton, X neutrons, no pions
       r = ( sum(event%numberParticles(1,:)).eq.0 .and. &
            & event%numberParticles(numStableMesons+1,1).eq.1)
    case (3)  ! 1 pi+, no other pions, but possible mesons of other flavor
       r = ( sum(event%numberParticles(1,:)).eq.1 .and. &
           & event%numberParticles(1,1).eq.1)
    case (4) ! 1 or more pi+, X other pions
       r = ( event%numberParticles(1,1).ge.1 )
    case (5)  ! 1 pi0, no other pions, but possible mesons of other flavor
       r = ( sum(event%numberParticles(1,:)).eq.1 .and. &
           & event%numberParticles(1,0).eq.1)
    case (6) ! 1 or more pi0, X other pions or other mesons
       r = ( event%numberParticles(1,0).ge.1 )
    case (7) ! 1 pi-, no other pi-, but possibly mesons of other flavor
       r = ( sum(event%numberParticles(1,:)).eq.1 .and. &
           & event%numberParticles(1,-1).eq.1)
    case (8) ! 1 or more pi-, X other mesons
        r = ( event%numberParticles(1,-1).ge.1 )
    case (9) ! 2 protons, no neutrons, no pions
       r = ( sum(event%numberParticles(1,:)).eq.0 .and. &
            & event%numberParticles(numStableMesons+1,1).eq.2 .and. &
            & event%numberParticles(numStableMesons+1,0).eq.0)
    case (10) ! 1 proton 1 neutron, no pions
       r = ( sum(event%numberParticles(1,:)).eq.0 .and. &
            & event%numberParticles(numStableMesons+1,1).eq.1 .and. &
            & event%numberParticles(numStableMesons+1,0).eq.1)
    case (11) ! no protons, 2 neutrons, no pions
       r = ( sum(event%numberParticles(1,:)).eq.0 .and. &
            & event%numberParticles(numStableMesons+1,0).eq.2 .and. &
            & event%numberParticles(numStableMesons+1,1).eq.0)
    case (12) ! 2 protons, any number of neutrons, no pions
       r = ( sum(event%numberParticles(1,:)).eq.0 .and. &
            & event%numberParticles(numStableMesons+1,1).eq.2 )
    case (13) ! 2 neutrons, any number of protons, no pions
       r = ( sum(event%numberParticles(1,:)).eq.0 .and. &
            & event%numberParticles(numStableMesons+1,0).eq.2 )
    case (14) ! 3 protons, X neutrons,  no pions
       r = ( sum(event%numberParticles(1,:)).eq.0 .and. &
            &event%numberParticles(numStableMesons+1,1).eq.3 )
    case (15) ! 4 protons, X neutrons,  no pions
       r = ( sum(event%numberParticles(1,:)).eq.0 .and. &
            & event%numberParticles(numStableMesons+1,1).eq.4 )
    case (16) ! 1 proton, 0 neutrons, no pions
       r = ( sum(event%numberParticles(1,:)).eq.0 .and. &
            & event%numberParticles(numStableMesons+1,0).eq.0  .and.  &
            & event%numberParticles(numStableMesons+1,1).eq.1 )
    case (17) ! 0 protons, 1 neutron, no pions
       r = ( sum(event%numberParticles(1,:)).eq.0 .and. &
            & event%numberParticles(numStableMesons+1,1).eq.0 .and. &
            & event%numberParticles(numStableMesons+1,0).eq.1)
    case (18) ! 0 protons, X neutrons, no pions
       r = ( sum(event%numberParticles(1,:)).eq.0 .and. &
            &event%numberParticles(numStableMesons+1,1).eq.0 &
            & .and. event%numberParticles(numStableMesons+1,0) > 0 )
    case (19) ! exclusive pi0: 1 pi0 and no other mesons of any flavor
       r = ( sum(event%numberParticles(1:numStableMesons,:)).eq.1 .and. &
            & event%numberParticles(1,0).eq.1)
    case (20) ! exclusive pi+: 1 pi+ and no other mesons of any flavor
       r = ( sum(event%numberParticles(1:numStableMesons,:)).eq.1 .and. &
            & event%numberParticles(1,1).eq.1)
    case (21) ! exclusive pi-: 1 pi- and no other mesons of any flavor
       r = ( sum(event%numberParticles(1:numStableMesons,:)).eq.1 .and. &
            & event%numberParticles(1,-1).eq.1) 
    case default
       r = .false.
       write(*,*) 'WRONG specificEvent=', specificEvent
       call Traceback()
    end select
  end function IfPass_SpecificEvent



 ! This subroutine sets the switch 'exclusive_hadron' for
 ! truly exclusive 1 meson production
 ! The calling variable 'A' is read in neutrinoAnalysis

   subroutine set_Exclusive(A)
   logical, intent(in) :: A 
   
   exclusive_hadron = A
   
   end subroutine set_Exclusive
   
 ! This subroutine sets the switches for the analysis of QE-like events 
 ! with 1 muon, 0 pion, (at least) 1 proton
 ! The calling variables are read in neutrinoAnalysis  
 
   subroutine set_QElike(QEp,pcut)
   logical, intent(in) :: QEp
   real, intent(in) :: pcut
   
   Q2p = QEp
   nuc_mom_pcut = pcut 
 
   end subroutine set_QElike


end module AnaEvent
