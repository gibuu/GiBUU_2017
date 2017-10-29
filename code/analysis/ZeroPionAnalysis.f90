!******************************************************************************
!****m* /ZeroPionAnalysis
! NAME
! module ZeroPionAnalysis
!
! PURPOSE
! This module does the analysis of the output of lepton induced processes
! selecting the events with 0pions  potentially interesting for LAr detectors,
! which are sensitive to outgoing nucleons of kinetic energy > 25-30 MeV
!
! INPUTS
! No Namelist available.
!
! USES
! We use the module AnaEvent.f90
!
! NOTES
! ZeroPionAnalysis = the code of AnaEvent.f90  is just adapted for additional
! final states. It was created as separate module because there were no vacant
! colums for extra final states left in sigma.dat. The output is in files
! sigma_0pions.dat, neutrino_0pions.dat, neutrino_0pions_QE.dat and so on
!
!******************************************************************************
module ZeroPionAnalysis

  use AnaEventDefinition
  use AnaEvent

  implicit none

  private

  public :: event_sigma_0pions, event_dSigma_dE_0pions

contains



  !****************************************************************************
  !****s* ZeroPionAnalysis/event_sigma_0pions
  ! NAME
  ! subroutine event_sigma_0pions(E,sigma,newInit,runNumber,identifier)
  !
  ! PURPOSE
  ! the same as event_sigma in the module AnaEvent
  ! the channels for analysis are different
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
  !
  ! OUTPUT
  ! * filename 'sigma_0pions.dat'
  !****************************************************************************
  subroutine event_sigma_0pions(E,sigmaOut,newInit,runNumber,identifier)
    use output, only: intToChar
    use particleDefinition
    use particleProperties, only: PartName
    use initNeutrino, only: get_init_namelist

    ! Input
    real, optional :: identifier
    integer, parameter :: dimSigma=100      ! number of possible channels
    type(tAnaEvent), intent(in), dimension(:) :: E
    real, dimension(1:dimSigma,1:2),intent(inout),optional :: sigmaOut
    logical,intent(in),optional :: newInit
    integer,intent(in),optional :: runNumber
    real :: perweight

    ! Field to store the Xsections in
    real, dimension(1:dimSigma,1:2) :: sigma

    character(31), dimension(1:dimSigma,1:2) :: sigma_name
    ! Denotes the different channels, used for output
    integer :: channel,i,j,k,l,outLeptonID,outLeptoncharge
    type(particle) :: part
    character(20) :: formatOut
    logical :: initSigma=.true.

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
          sigma(channel,1)=identifier ! the first column will be running variable
       else
          sigma(channel,1)=runNumber
       end if
       sigma(channel,2)=0
    else
       sigma(:,1)=0.
       sigma(:,2)=0.
    end if


    !**************************************************************************
    ! Special channels: muon, but no pions
    !**************************************************************************
    channel=channel+1
    sigma_name(channel,1)=intTochar(channel)//': mu w/o pi'
    sigma_name(channel,2)=intTochar(channel+dimSigma)//': error mu w/o pi'

    do j=lbound(E,dim=1),ubound(E,dim=1)
       if (sum(E(j)%numberParticles(1,-1:1)).eq.0) then
          perweight=0.
          if (event_GetParticle(E(j),outLeptonID,outLeptonCharge,1,part)) &
                      &  perweight=part%perweight
          sigma(channel,1)=sigma(channel,1)+perweight
       end if
    end do


    !**************************************************************************
    ! Special channels: no pions and "k" protons and "i" neutrons
    ! and anything else
    ! = 0 pions  "k" protons  "i" neutrons  ! for LAr detector
    !**************************************************************************

    do k=0,4

    do i=0,4

       if ((k+i) .gt. 4) cycle

    channel=channel+1
    sigma_name(channel,1)=  &
           & intTochar(channel)//': 0pi '//intTochar(k)//'p '//intTochar(i)//'n'
    sigma_name(channel,2)=  &
           & intTochar(channel+dimSigma)//': error 0pi '//intTochar(k)//'p ' &
           & //intTochar(i)//'n'


    do j=lbound(E,dim=1),ubound(E,dim=1)
      if (sum(E(j)%numberParticles(1,-1:1)).eq.0) then ! no pions
          if (E(j)%numberParticles(7,1).eq.k) then ! "k" protons
              if (E(j)%numberParticles(7,0).eq.i) then ! "i" neutrons
                 perweight=0.
                 if (event_GetParticle(E(j),outLeptonID,outLeptonCharge,1,part))&
                          &  perweight=part%perweight
                 !write(*,*) 'found event number j=',j, ' with  k=',k, 'protons&
                 !       & and i=', i, '  neutrons', &
                 !       & ' muon charge is', charge, '   perweight=', perweight
                 sigma(channel,1)=sigma(channel,1)+perweight
              end if  ! neutrons
           end if ! protons
      end if ! pions
    end do
    end do ! i
    end do ! k



    !**************************************************************************
    ! Special channels: 2 pions
    !**************************************************************************

    do k=0,2
    do i=0,2
       if ((k+i) .gt. 2) cycle
    l=2-i-k

    channel=channel+1
    sigma_name(channel,1)=intTochar(channel)//': '//intTochar(k)//'pi- '   &
                          & //intTochar(i)//'pi0 '//intTochar(l)//'pi+'
    sigma_name(channel,2)=intTochar(channel+dimSigma)//': error '   &
              & //intTochar(k)//'pi- '//intTochar(i)//'pi0 '//intTochar(l)//'pi+'


        do j=lbound(E,dim=1),ubound(E,dim=1)
          if (sum(E(j)%numberParticles(1,-1:1)).eq.2) then ! 2 pions
              if (E(j)%numberParticles(1,-1).eq.k) then ! "k" pi-
                  if (E(j)%numberParticles(1,0).eq.i) then ! "i" pi0
                     if (E(j)%numberParticles(1,1).eq.l) then ! "l" pi+
                        perweight=0.
                        ! read perweight from the first particle in the list
                        if (associated(E(j)%particleList%first))  &
                                 & perweight=E(j)%particleList%first%V%perweight
                        !write(*,'(5(A,I5),A,f12.4)') 'found event number j=',j,&
                        !& '  with  k=',k, ' pi- and i=', i, ' pi0  and l=',l,  &
                        !& 'pi+ ;     nucleon charge is', charge, &
                        !&'   perweight=', perweight(1)
                        sigma(channel,1)=sigma(channel,1)+perweight
                        !write (*,*) 'channel=', channel, '  sigma=',&
                        !& sigma(channel,1)
                     end if  ! pi+
                  end if  ! pi0
               end if  ! pi-
          end if ! pions
        end do
    end do ! i
    end do ! k



    !**************************************************************************
    ! Special channels: 1 pion of a given charge  X other pions
    !**************************************************************************

    do k=-1,1

    channel=channel+1
    sigma_name(channel,1)=intTochar(channel)//': 1'//  &
                          & trim(PartName(101,k,.false.))//' X'
    sigma_name(channel,2)=intTochar(channel+dimSigma)//': error 1'//  &
                          & trim(PartName(101,k,.false.))//' X'


    do j=lbound(E,dim=1),ubound(E,dim=1)
       if (E(j)%numberParticles(1,k).eq.1) then ! "k" is the charge of the pion
          perweight=0.
          if (event_GetParticle(E(j),1,k,1,part)) perweight=part%perweight
          ! read perweight from outgoing pion
          sigma(channel,1)=sigma(channel,1)+perweight
       end if
    end do
    end do ! k


    !**************************************************************************
    ! Special channels: no pions and "k" protons and any number of  neutrons
    ! and anything else
    ! = 0 pions  "k" protons  "X" neutrons ! for LAr detector
    !**************************************************************************

    do k=0,12

    channel=channel+1
    sigma_name(channel,1)=intTochar(channel)//': 0pi '//intTochar(k)//'p X'
    sigma_name(channel,2)=intTochar(channel+dimSigma)//': error 0pi ' &
                           & //intTochar(k)//'p X'


    do j=lbound(E,dim=1),ubound(E,dim=1)
       if (sum(E(j)%numberParticles(1,-1:1)).eq.0) then ! no pions
          if (E(j)%numberParticles(7,1).eq.k) then ! "k" protons
             perweight=0.
             if (event_GetParticle(E(j),outLeptonID,outLeptonCharge,1,part))  &
                   &  perweight=part%perweight
             !write(*,*) 'found event number j=',j, '  with  k=',k, &
             ! & 'protons X neutrons', &
             ! & ' nucleon charge is', charge, '   perweight=', perweight(1)
             sigma(channel,1)=sigma(channel,1)+perweight
          end if ! protons
       end if ! pions
    end do ! j
    end do ! k



    !**************************************************************************
    !
    ! Output Results to files and to input arrays:
    !
    !**************************************************************************

    if (initSigma) then
       initSigma=.false.
       !***********************************************************************
       !****o* AnaEvent/sigma_0pions.dat
       ! NAME
       ! file sigma_0pions.dat
       ! PURPOSE
       ! * First column  : Naming for different runs or runNumber
       ! * Columns 2-... : Single and multiparticle cross sections.
       !                   The header of the file names the explicit channels.
       ! NOTES
       ! the same logic as sigma.dat (=output of AnaEvent.f90), here different
       ! channels as output
       !
       !***********************************************************************

       open(14,file='sigma_0pions.dat')
       do j=1,channel
          write(14,'("# ",A,"    ",A)') sigma_name(j,1),sigma_name(j,2)
       end do
       formatOut='("#",A8,'//intTochar(2*channel)//'A14)'
       write(14,formatOut)  sigma_name(:channel,1)(1:4),   &
                          & sigma_name(:channel,2)(1:4)
!       formatOut='(A,'//intTochar(2*channel)//'A20)'
!       write(14,formatOut)   '# ',sigma_name(:channel,1),sigma_name(:channel,2)
     else
       open(14,file='sigma_0pions.dat',position='append')
    end if


    if (present(sigmaOut).and.present(runNumber).and.present(newInit)) then

       ! For error evaluation we save the square of the change
       ! in the present run into sigma(:,2) :
       do i=1,100
          sigma(i,2)=sigma(i,2)+(sigma(i,1)-sigmaOut(i,1))**2
       end do
       sigmaOut=sigma
       do i=1,100
          if (runNumber.gt.1) then
             ! Statististical error of the mean value of several runs:
             sigma(i,2)=sqrt(max(((sigma(i,2)-sigma(i,1)**2/float(runNumber)) &
                        & /float(runNumber-1)/float(runNumber)),0.))
          else
             sigma(i,2)=0.
          end if
       end do
       !avoid output since the files is large and of no immediate use
       !formatOut='(A,F8.4,'//intTochar(dimSigma+channel)//'E14.5)'
       !write(14,formatOut) ' ',sigma(1,1),sigma(2:,1) &
       !                    & /float(runNumber),sigma(:channel,2)

    else
       ! NO ERROR EVALUATION, SINCE NO KNOWLEDGE ABOUT EARLIER RUNS

       !avoid output since the files is large and of no immediate use
       !formatOut='(A,'//intTochar(channel)//'E14.5)'
       !write(14,formatOut) ' ',sigma(1:channel,1)
    end if
    close(14)


  end subroutine event_sigma_0pions



! for ZeroPion analysis, where we are iterested only in events with 0 pions
subroutine event_dSigma_dE_0pions(E,eMin,eMax,dE,string,runNumber,hists, &
                                   & initHists,makeOutputIn,&
                                   & sameFileNameIn,histsMulti,hists1X,hists2X)
    use particleDefinition
    use IdTable, only: isHadron
    use particleProperties, only: hadron, validCharge_ID
    use output, only: intTochar_pm,intToChar
    use histf90

    logical , intent(in), optional ::sameFileNameIn
    ! .true. = always print to same filenames,
    ! .false. = different files for different runs

    logical , intent(in), optional ::makeOutputIn
    ! Flag to switch off output : .false. = no output

    real, intent(in) :: eMin, eMax, dE
    character(*), intent(in) :: string ! used as prefix for all output files
    integer, intent(in) :: runNumber ! number of the run
    type(tAnaEvent) , dimension(:) :: E ! List of events
    type(particle) :: part
    type(histogram) :: dsigma_dE,dsigma_dE_multi, dsigma_dE_1X, dsigma_dE_2X
    integer :: charge,i,j, particle_number
    character(60) :: filename,name,filenameMulti, filename1X, filename2X
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
       write(*,*) 'Error in event_dSigma_0pions_dE. runNumber.lt.1'
       stop
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Analysis for 0 pion events so far only for the nucleon, can be generalized by
! introducing loop over i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       i=numStableMesons + 1  ! =nucleon   
       
       chargeLoop: do charge=0,1
          if (validCharge_ID(particleIDs(i),charge)) then
             ! Initialize histogram

             makeInitHists=.true.
             if (present(hists).and.present(initHists)) then
                if (.not.initHists) then
                   ! Do not initialze histograms, but use input as basis.
!                   write(*,*)'  makeInitHists=.false.'
                   makeInitHists=.false.
                end if
             end if

             if (isHadron(particleIds(i))) then
                name=trim(hadron(particleIDs(i))%name)
             else
                write(*,*) 'severe problem in event_dSigma_0pions_dE! Stop!'
                stop
             end if


             if (makeInitHists) then
                call createHist(dsigma_dE, 'dSigma_0pions/dEkin  &
                     & ('// trim(name)//') for single            &
                     & '//trim(name)// ' production'             &
                     &  ,eMin ,eMax,dE)
                call createHist(dsigma_dE_1X, 'dSigma_0pions/dEkin  &
                     & ('// trim(name) //                           &
                     & ') for 1 '//trim(name)//                     &
                     & 'and X(including the same id but different charge) &
                     &  production' &
                     &  ,eMin ,eMax,dE)
                call createHist(dsigma_dE_2X, 'dSigma_0pions/dEkin( &
                     & '// trim(name) //&
                     & ') for 2 '//trim(name)// 'and X(including the same id &
                     & but different charge) production' &
                     &  ,eMin ,eMax,dE)
                call createHist(dsigma_dE_multi, 'dSigma_0pions/dEkin  &
                     & ('// trim(name) //&
                     & ') for multi '//trim(name)// ' production' &
                     &  ,eMin ,eMax,dE)
             else
                dsigma_dE=hists(i,charge)
                dsigma_dE_1X=hists1X(i,charge)
                dsigma_dE_2X=hists2X(i,charge)
                dsigma_dE_multi=histsMulti(i,charge)
             end if

             if (sameFileName) then
                filename=trim(string)//'_dSigma_dEkin_'// trim(name)&
                     & //'_charge_'//intTochar_pm(charge) //'_no_pi.dat'
                filename1X=trim(string)//'_dSigma_dEkin_'// trim(name)&
                     & //'_charge_'//intTochar_pm(charge) //'_1X_no_pi.dat'
                filename2X=trim(string)//'_dSigma_dEkin_'// trim(name)&
                     & //'_charge_'//intTochar_pm(charge) //'_2X_no_pi.dat'
                filenameMulti=trim(string)//'_dSigma_dEkin_'//trim(name)&
                     & //'_charge_'//intTochar_pm(charge) //'_MULTI_no_pi.dat'
             else
                filename=trim(string)//'_dSigma_dEkin_'// trim(name)&
                     & //'_charge_'//intTochar_pm(charge) //'.'   &
                     & //trim(intTochar(runNumber))//'_no_pi.dat'
                filename1X=trim(string)//'_dSigma_dEkin_'// trim(name)&
                     & //'_charge_'//intTochar_pm(charge) //'.'   &
                     & //trim(intTochar(runNumber))//'_1X_no_pi.dat'
                filename2X=trim(string)//'_dSigma_dEkin_'// trim(name)&
                     & //'_charge_'//intTochar_pm(charge) //'.'   &
                     & //trim(intTochar(runNumber))//'_2X_no_pi.dat'
                filenameMulti=trim(string)//'_dSigma_dEkin_'//trim(name)&
                     & //'_charge_'//intTochar_pm(charge) //'.'   &
                     & //trim(intTochar(runNumber))//'_MULTI_no_pi.dat'
             end if

       ! Loop over events
       eventLoop: do j=lbound(E,dim=1),ubound(E,dim=1)
       if (  sum(E(j)%numberParticles(1,:)).eq.0  ) then   ! no pions

               ! for single-particle cross sections
              !Find out whether there are any particles in the event
              ! with charge="charge" and ID="particleIDs(i):
               if (( E(j)%numberParticles(i,charge).eq.1).and. &
                  &( sum(E(j)%numberParticles(i,:)).eq.1)) then
                  if ( event_GetParticle(E(j),particleIDs(i),charge,1,part)) then
                      ekin=part%momentum(0)-part%mass
                      call addHist(dsigma_dE       , ekin       ,part%perweight)
                      !write(*,*) 'Single-particle: event j=',j,  &
                      !& '   charge=',charge, '  Ekin=',ekin
                   else
                      write(*,*) 'Error in event_dSigma_0pions_dE',  &
                                 &  E(j)%numberParticles(i,charge), 'Stop!'
                      stop
                   end if
                end if


                ! for 1-particle-plus-X cross sections
                if (E(j)%numberParticles(i,charge).eq.1) then
                   if ( event_GetParticle(E(j),particleIDs(i),charge,1,part)) then
                      ekin=part%momentum(0)-part%mass
                      call addHist(dsigma_dE_1X       , ekin    ,part%perweight)
                      !!write(*,*) '1-particle-plus-X: event j=',j,  &
                      ! & '   charge=',charge, '  Ekin=',ekin
                   else
                      write(*,*) 'Error in event_dSigma_0pions_dE',  &
                                 & E(j)%numberParticles(i,charge), 'Stop!'
                      stop
                   end if
                end if


                ! for 2-particle-plus-X cross sections
                if (E(j)%numberParticles(i,charge).eq.2) then
                     do particle_number=1,2
                          if ( event_GetParticle(E(j),particleIDs(i),charge,  &
                              & particle_number,part)) then
                             ekin=part%momentum(0)-part%mass
                             call addHist(dsigma_dE_2X , ekin  , &
                                & part%perweight*E(j)%numberParticles(i,charge))
                              !!write(*,*) 'Multi-particle: event j=',j,&
                              ! &'   charge=',charge, '  Ekin=',ekin
                           else
                            write(*,*) 'Error in event_dSigma_0pions_dE_multi',&
                                       & E(j)%numberParticles(i,charge), 'Stop!'
                            stop
                           end if
                     end do
                end if


                ! for multi-particle cross sections
                !Find out whether there are any particles in the event
                ! with charge="charge" and ID="particleIDs(i):
                if (E(j)%numberParticles(i,charge).ge.1) then
                     do particle_number=1,E(j)%numberParticles(i,charge)
                          if ( event_GetParticle(E(j),particleIDs(i),charge,  &
                                 & particle_number,part)) then
                             ekin=part%momentum(0)-part%mass
                             call addHist(dsigma_dE_multi , ekin  ,          &
                                & part%perweight*E(j)%numberParticles(i,charge))
                              !!write(*,*) 'Multi-particle: event j=',j, &
                              ! &'   charge=',charge, '  Ekin=',ekin
                           else
                            write(*,*) 'Error in event_dSigma_0pions_dE_multi',&
                                       & E(j)%numberParticles(i,charge), 'Stop!'
                            stop
                           end if
                     end do
                end if
       end if ! no pions
       end do eventLoop


             if (present(hists)) then
                if (.not.makeInitHists) then
                   ! Calculate the square of the changes compared to last run:
                   call makeError_hist(hists(i,charge),dsigma_dE)
                   ! Save the list
                   hists(i,charge)=dsigma_dE
                   ! Set yVal(:,3) to the statistical error
                   ! before doing the output :
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
                   ! Set yVal(:,3) to the statistical error
                   ! before doing the output :
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
                   ! Set yVal(:,3) to the statistical error
                   ! before doing the output :
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
                   ! Set yVal(:,3) to the statistical error
                   ! before doing the output :
                   call setError_hist(dsigma_dE_multi,runNumber)
                else
                   call makeError_hist(b=dsigma_dE_multi)
                   ! Save the list
                   histsMulti(i,charge)=dsigma_dE_multi
                end if
             end if

             ! Output histogram
             if (makeOutput) then
                open(14,file=filename)
                if (makeInitHists) then
                   call WriteHist(dsigma_dE,14,0.,1.)
                else
                   ! normalize to number of runs
                   call WriteHist(dsigma_dE,14,0.,1./float(runNumber))
                end if
                close(14)

                open(14,file=filename1X)
                if (makeInitHists) then
                   call WriteHist(dsigma_dE_1X,14,0.,1.)
                else
                   ! normalize to number of runs
                   call WriteHist(dsigma_dE_1X,14,0.,1./float(runNumber))
                end if
                close(14)

                open(14,file=filename2X)
                if (makeInitHists) then
                   call WriteHist(dsigma_dE_2X,14,0.,1.)
                else
                   ! normalize to number of runs
                   call WriteHist(dsigma_dE_2X,14,0.,1./float(runNumber))
                end if
                close(14)

                open(14,file=filenameMulti)
                if (makeInitHists) then
                   call WriteHist(dsigma_dE_multi,14,0.,1.)
                else
                   ! normalize to number of runs
                   call WriteHist(dsigma_dE_multi,14,0.,1./float(runNumber))
                end if
                close(14)
             end if


             call RemoveHist(dsigma_dE)
             call RemoveHist(dsigma_dE_1X)
             call RemoveHist(dsigma_dE_2X)
             call RemoveHist(dsigma_dE_multi)

          end if
       end do chargeLoop
  end subroutine event_dSigma_dE_0pions



end module ZeroPionAnalysis
