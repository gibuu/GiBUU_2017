!******************************************************************************
!****m* /Annihilation
! NAME
! module Annihilation
! PURPOSE
! Includes the routines performing the baryon-antibaryon annihilation.
!******************************************************************************
module Annihilation

  implicit none
  private

  !****************************************************************************
  !****g* Annihilation/model
  ! SOURCE
  !
  integer, save :: model=2
  ! PURPOSE
  ! Switch between the models of annihilation:
  ! * 1 -- string based model,
  ! * 2 -- statistical model
  !****************************************************************************

  !****************************************************************************
  !****g* Annihilation/position_flag
  ! SOURCE
  !
  integer, save :: position_flag=1
  ! PURPOSE
  ! Switch between the choices of position of outgoing mesons:
  ! * 1 -- at the c.m. of the baryon and antibaryon,
  ! * 2 -- at the antibaryon position
  !****************************************************************************

  public :: annihilate, annihilate_to_meson


contains

  !****************************************************************************
  !****s* Annihilation/readInput
  ! NAME
  ! subroutine readInput
  !
  ! PURPOSE
  ! Reads input in jobcard out of namelist "annihilation"
  !****************************************************************************
  subroutine readInput

    use output

    integer :: ios

    !**************************************************************************
    !****n* Annihilation/annihilation
    ! NAME
    ! NAMELIST /annihilation/
    ! PURPOSE
    ! Includes the switches:
    ! * model
    ! * position_flag
    !**************************************************************************
    NAMELIST /annihilation/ model,position_flag

    call Write_ReadingInput('annihilation',0)
    rewind(5)
    read(5,nml=annihilation,IOSTAT=ios)
    call Write_ReadingInput('annihilation',0,ios)

    write(*,*) ' Model (1 -- string-based, 2 -- statistical): ', model
    write(*,*) ' position_flag (1 --- at c.m. , 2 --- at AntiBar position): ', position_flag

    call Write_ReadingInput('annihilation',1)

  end subroutine readInput



  !****************************************************************************
  !****s* Annihilation/annihilate
  ! NAME
  ! subroutine annihilate(antibar,bar,time,finalState,flag_anni,HiEnergyType)
  ! PURPOSE
  ! Simulates an annihilation event of an antibaryon and a baryon into some finalState.
  ! INPUTS
  ! * type(particle), intent(in) ::  antibar                      ! antibaryon particle
  ! * type(particle), intent(in) ::  bar                          ! baryon particle
  ! * real, intent(in)                              :: time       ! current time
  ! OUTPUT
  ! * type(particle), intent(inout), dimension(:) :: finalState   ! Array of the final state particles
  ! * logical, intent(out) :: flag_anni                           ! .true. if annihilation successfull
  ! * integer, intent(out) :: HiEnergyType                        ! must be=-3 at the output
  ! *                                                             ! (see docu of setKinematicsHiEnergy)
  ! NOTES
  ! This subroutine has a similar structure as generateFinalState. The main difference
  ! is that the collision criterion is not checked here, i.e. the annihilation
  ! is done independently on the actual annihilation cross section.
  !****************************************************************************
    subroutine annihilate(antibar,bar,time,finalState,flag_anni,HiEnergyType)

    use particleDefinition, only: particle, sqrtS
    use particleProperties, only: hadron
    use IdTable
    use lorentzTrafo
    use twoBodyTools, only: sqrtS_free
    use RMF, only: getRMF_flag
    use propagation, only: updateVelocity, checkVelo
    use densitymodule, only: Particle4Momentum_RMF
    use mediumDefinition
    use energyCalc, only: energyCorrection
    use nBodyPhaseSpace, only: momenta_in_nBodyPS
    use random
    use Coll_BaB
    use hadronFormation, only: forceInitFormation
    use PythiaSpecFunc, only: Init_VM_Mass
    use annihilation_channels
    use finalState_Full, only: massass_nBody
    use constants, only: mPi

    type(particle), intent(in) ::  antibar                      ! antibaryon particle
    type(particle), intent(in) ::  bar                          ! baryon particle
    real, intent(in)                              :: time       ! current time
    type(particle), intent(inout), dimension(:) :: finalState   ! Array of the final state particles
    logical, intent(out) :: flag_anni                           ! .true. if annihilation successfull
    integer, intent(out) :: HiEnergyType                        ! must be=-3 at the output
                                                                ! (see docu of setKinematicsHiEnergy)

    real :: srtS_star, srtS_vacuum, srtS, mstar_antibar, mstar_bar, xx
    real, dimension(0:3) :: momentum_star, momentum_tot, momentum_antibar, momentum_bar
    real, dimension(0:3) :: momentum_final, p_cm, impuls
    real, dimension(1:3) :: betaToCV, betaToCM, position
    real, dimension(3,6)  :: pn      ! first index : particle, second index : momentum index

    type(medium) :: mediumAtColl
    real, dimension(1:3) :: betaToLRF

    logical :: successFlag
    integer :: i, k, maxId, izt, totalCharge, totalStrangeness, totalCharm, model1

    real, parameter :: branching_2pi=0.0038, branching_3pi=0.025  ! Branching ratios of
                         ! the 2pi and 3pi channels for p pbar annihilation at rest.
                         ! See Table VI in I.N. Mishustin et al, PRC 71, 035201 (2005)

    integer, save :: ncalls=0
    logical, parameter :: debugFlag=.false.

    integer, save :: nBlocked=0, nTotal=0

    ncalls=ncalls+1

    if (debugFlag) write(*,*) 'In annihilate: nBlocked, nTotal:', nBlocked, nTotal

    if (ncalls.eq.1) then
       call readInput
       call forceInitFormation
    end if

    if ( .not.antibar%antiparticle ) then
      write(*,*) ' wrong antibaryon in annihilate: ', antibar%Id, antibar%charge, antibar%antiparticle
      stop
    else if ( bar%antiparticle ) then
      write(*,*) ' wrong baryon in annihilate: ', bar%Id, bar%charge, bar%antiparticle
      stop
    end if

    ! Dummy settings: *************************
    mediumAtColl%useMedium=.false.
    betaToLRF(1:3)=0.
    ! *****************************************
    ! In the case of RMF mode momentum_star is the sum of the kinetic
    ! 4-momenta and betaToCV is the speed of the center-of-velocity frame:

    momentum_star(0:3)= antibar%momentum(0:3) + bar%momentum(0:3)
    betaToCV = lorentzCalcBeta (momentum_star, 'annihilate(1)')

    srtS_star = sqrtS((/antibar,bar/),"annihilate, srtS_star")

    if (.not.getRMF_flag() ) then

      srtS_vacuum=sqrtS_free((/antibar,bar/))
      momentum_tot=momentum_star
      betaToCM=betaToCV
      srtS=srtS_star

    else   ! RMF mode:

      mstar_antibar = sqrtS(antibar,'annihilate, mstar_antibar')
      mstar_bar = sqrtS(bar,'annihilate, mstar_bar')

      srtS_vacuum = srtS_star - mstar_antibar - mstar_bar  + antibar%mass + bar%mass

      ! Compute canonical momenta of the antibaryon and baryon:

      call Particle4Momentum_RMF(antibar,momentum_antibar)
      call Particle4Momentum_RMF(bar,momentum_bar)

      ! Compute srtS and betaToCM  with canonical momenta:

      momentum_tot(0:3)= momentum_antibar(0:3) + momentum_bar(0:3)
      srtS= momentum_tot(0)**2 - dot_product(momentum_tot(1:3),momentum_tot(1:3))
      if (srtS.le.0.) then
        if (debugFlag) then
          write(*,*) '# Space-like total 4-momentum, s**2: ', srtS
          write(*,*) '# srtS_star: ', srtS_star
          write(*,*) '# Bar. and antibar. Dirac masses: ', mstar_bar, mstar_antibar
          write(*,*) '# baryon canonical 4-momentum: ',  momentum_bar(:)
          write(*,*) '# antibaryon canonical 4-momentum: ',  momentum_antibar(:)
          write(*,*) '# baryon kinetic 4-momentum: ', bar%momentum(:)
          write(*,*) '# antibaryon kinetic 4-momentum: ', antibar%momentum(:)
          write(*,*) '# Vector field for the baryon: ',     momentum_bar(:)-bar%momentum(:)
          write(*,*) '# Vector field for the antibaryon: ', momentum_antibar(:)-antibar%momentum(:)
        end if
        flag_anni=.false.
        return
      end if

      srtS= sqrt(srtS)

      betaToCM = lorentzCalcBeta (momentum_tot, 'annihilate(2)')

    end if

!*** In some cases the statistical model is not applicable anyway
!*** and we have to rely on string model:

    totalCharge= antibar%charge + bar%charge
    totalStrangeness= hadron(bar%Id)%strangeness - hadron(antibar%Id)%strangeness
    totalCharm= hadron(bar%Id)%charm - hadron(antibar%Id)%charm

    if (abs(totalCharge).gt.1 .or. totalStrangeness.ne.0 .or. totalCharm.ne.0) then
       model1=1
    else
       model1=model
    end if

    select case (position_flag)
    case (1)
       position=(antibar%position+bar%position)/2.
    case (2)
       position=antibar%position
    case default
        write(*,*) ' Position of annihilation final state is not defined: ', position_flag
        stop
    end select

    finalState%Id=0               ! Initialize output

    select case (model1)

    case (1)  ! String-based model

        if ( srtS .gt. 4.*mPi ) then   ! Use HiEnergy package:

              ! Define CM - Momentum :
              p_cm=bar%momentum
              call lorentz(betaToCV, p_cm, 'annihilate(1)')

              if (debugFlag) write(*,*) 'Vor DoColl_BaB'
              call Init_VM_Mass(srtS_vacuum,position)
              call DoColl_BaB((/bar,antibar/),finalState,flag_anni,srtS_vacuum,p_cm,(/0.,0.,0./))
              if (debugFlag) write(*,*) 'Nach DoColl_BaB'

              call ResetPosition ! this also sets maxID

        else if ( srtS .gt. 2.*mPi ) then  ! Annihilation to 2pi or 3pi:

              izt= antibar%charge + bar%charge

              xx=rn()

              if ( xx .lt. branching_2pi ) then

                ! Generate  2pi final state:

                finalState(1:2)%ID=pion
                finalState(1:2)%mass=mPi

                if (izt.eq.0) then
                  xx = rn()
                  ! 0.82=the branching ratio of the pi^+ pi^- channel for p pbar annihilation
                  ! at rest to 2pi.
                  ! See Table VI in I.N. Mishustin et al, PRC 71, 035201 (2005)
                  if (xx.lt.0.82) then
                    finalState(1)%charge = 1
                    finalState(2)%charge = -1
                  else
                    finalState(1:2)%charge = 0
                  end if
                else if (abs(izt).eq.1) then
                  finalState(1)%charge = izt
                  finalState(2)%charge = 0
                else if (abs(izt).eq.2) then
                  finalState(1:2)%charge = izt/2
                else
                  flag_anni=.false.
                  return
                end if

              else if ( xx .lt. branching_2pi+branching_3pi ) then

                ! Generate 3pi final state:

                finalState(1:3)%ID=pion
                finalState(1:3)%mass=mPi

                if (izt.eq.0) then
                  xx = rn()
                  ! 0.72=the branching ratio of the pi^+ pi^- pi^0 direct channel for p pbar annihilation
                  ! at rest to 3pi.
                  ! See Table VI in I.N. Mishustin et al, PRC 71, 035201 (2005)
                  if (xx.lt.0.72) then
                    finalState(1)%charge = 1
                    finalState(2)%charge = -1
                    finalState(3)%charge = 0
                  else
                    finalState(1:3)%charge = 0
                  end if
                else if (abs(izt).eq.1) then
                  ! No data, thus assume equal probabilities for the different charge channels:
                  finalState(1)%charge = izt
                  xx = rn()
                  if (xx.lt.0.5) then
                    finalState(2)%charge = 1
                    finalState(3)%charge = -1
                  else
                    finalState(2:3)%charge = 0
                  end if
                else if (abs(izt).eq.2) then
                  finalState(1:2)%charge = izt/2
                  finalState(3)%charge = 0
                else if (abs(izt).eq.3) then
                  finalState(1:3)%charge = izt/3
                else
                  flag_anni=.false.
                  return
                end if

              else

                ! srtS is too small to generate other final states
                flag_anni=.false.
                return

              end if

              call ResetPosition ! this also sets maxID

              if ( srtS.lt.sum(finalState(1:maxID)%mass)+0.001 ) then
                if (debugFlag) write(*,*) ' In annihilate: event is blocked: ',&
                               & srtS,finalState(1:maxID)%Id,sum(finalState(1:maxID)%mass)
                nBlocked=nBlocked+1
                flag_anni=.false.
                return
              end if

              ! Vacuum dispersion relations are used for outgoing pions below:
              call momenta_in_nBodyPS(srtS,finalState(1:maxID)%mass,pn(1:3,1:maxID))
              do i=1,maxID
                finalState(i)%momentum(1:3)=pn(1:3,i)
                finalState(i)%momentum(0)=sqrt( finalState(i)%mass**2 &
                                        &      +dot_product(pn(1:3,i),pn(1:3,i)) )
              end do

              flag_anni=.true.

        else

              maxID=0
              flag_anni =.false.

        end if


    case (2)  ! Statistical model

        if ( srtS .gt. 4.*mPi ) then

              !maxID=3
              !finalState(1:2)%ID=omegaMeson
              !finalState(3)%ID=pion
              !finalState(1:3)%charge = 0
              !successFlag=.true.

              call choose_channel(srtS_vacuum,antibar,bar,time,finalState,successFlag)

              call ResetPosition ! this also sets maxID

              if (debugFlag) then
                 write(*,*) 'In annihilation, before massass_nBody, srts: ', srts
                 do i=1,maxID
                   write(*,*) i, finalState(i)%Id
                 end do
              end if


              if (successFlag) then
                 call massass_nBody(srtS_vacuum,betaToLRF,betaToCM,mediumAtColl,&
                                   &finalState(1:maxID),flag_anni)
              else
                 flag_anni=.false.
              end if

              if (debugFlag) then
                 write(*,*) 'In annihilation, after massass_nBody, srts: ', srts
                 do i=1,maxID
                   write(*,*) i, finalState(i)%mass, finalState(i)%momentum(1:3)
                 end do
              end if

        else

              maxID=0
              flag_anni =.false.

        end if

    case default

        write(*,*) ' Annihilation model is not defined: ', model1
        stop

    end select


    if (debugFlag) then
        if (ncalls.eq.1) open(3,file='annihilation_event.dat',status='unknown')
        write(3,*)' Call #: ', ncalls
        write(3,*)' srtS_star, srtS, srtS_vacuum: ', srtS_star, srtS, srtS_vacuum
        if (maxID > 0) then
          write(3,*)' sum final masses: ', sum(finalState(1:maxID)%mass)
          do i=1,maxID
            write(3,*) finalState(i)%Id,finalState(i)%charge,finalState(i)%mass
          end do
        end if
    end if

    if (.not.flag_anni) then
        if (debugFlag) then
          write(*,*) ' In annihilate: annihilation is impossible !!!'
          write(*,*) ' antibar%Id, antibar%charge: ', antibar%Id, antibar%charge
          write(*,*) ' bar%Id, bar%charge: ', bar%Id, bar%charge
          write(*,*) ' srtS_star, srtS, srtS_vacuum: ', srtS_star, srtS, srtS_vacuum
        end if
        return
    end if


    nTotal=nTotal+1

    if ( srtS.lt.sum(finalState(1:maxID)%mass)+0.001 ) then
         if (debugFlag) write(*,*) ' In annihilate: event is blocked: ',&
                               & srtS,finalState(1:maxID)%Id,sum(finalState(1:maxID)%mass)
         nBlocked=nBlocked+1
         flag_anni=.false.
         return
    end if


    ! Lorentz boost / energy correction:

    if ( srtS .gt. 4.*mPi ) then

        !  energyCorrection takes finalState in the CM frame
        !  and returns it in the calculational frame

        call energyCorrection(srtS, betaToLRF, betaToCM, mediumAtColl, &
                             &finalState(1:maxId), successFlag)

        if (.not.successFlag) then
          flag_anni=.false.
          if (debugFlag) write(*,*) ' Energy correction failed '
          return
        end if

    else

        do i=1,maxID
           impuls=finalState(i)%momentum
           call lorentz(-betaToCM,impuls,'annihilate(2)')
           finalState(i)%momentum=impuls
        end do

    end if

    finalState(1:maxID)%productionTime = time

    if (model1.eq.1 .and. srtS.gt.4.*mPi ) then
        finalState(1:maxID)%in_Formation = .TRUE.
    else
        finalState(1:maxID)%formationTime = time
        finalState(1:maxID)%scaleCS=1.
        finalState(1:maxID)%in_formation=.false.
    end if

    do k=0,3
      momentum_final(k)=sum(finalState(1:maxId)%momentum(k))
    end do

    ! Check energy-momentum conservation:
    if ( max( abs(momentum_tot(0)-momentum_final(0)),&
       &     abs(momentum_tot(1)-momentum_final(1)),&
       &     abs(momentum_tot(2)-momentum_final(2)),&
       &     abs(momentum_tot(3)-momentum_final(3)) ) .gt. 1.e-03 ) then
      write(*,*) ' Problem with energy-momentum conservation in annihilate'
      write(*,*) ' antibar%Id, antibar%charge: ', antibar%Id, antibar%charge
      write(*,*) ' bar%Id, bar%charge: ', bar%Id, bar%charge
      write(*,*) ' srtS_star, srtS, srtS_vacuum: ', srtS_star, srtS, srtS_vacuum
      write(*,*) ' momentum_tot - momentum_final: ', momentum_tot - momentum_final
    end if

    if ( .not.getRMF_flag() ) then
       call updateVelocity(finalState)
    else
       do i = 1,maxId
         finalstate(i)%velocity = finalState(i)%momentum(1:3)/finalState(i)%momentum(0)
         if (.not. checkVelo(finalState(i))) then
           write(*,*) ' In annihilate: checkVelo failed !!!'
           flag_anni=.false.
           return
         end if
       end do
    end if

    HiEnergyType = -3


  contains

    subroutine ResetPosition
      integer :: i
      maxId=size(finalState,dim=1)
      do i=1,size(finalState,dim=1)
         if (finalState(i)%Id.eq.0) then
            maxID=i-1
            exit
         end if
      end do
      do i=1,maxID
        finalState(i)%position=position
      end do
    end subroutine ResetPosition

  end subroutine annihilate



  !****************************************************************************
  !****s* Annihilation/annihilate_to_meson
  ! NAME
  ! subroutine annihilate_to_meson(pairIn,time,finalState,flag_anni,HiEnergyType)
  ! PURPOSE
  ! Simulates an annihilation event of an antibaryon and a baryon into a meson.
  ! INPUTS
  ! * type(particle), dimension(1:2), intent(in) ::  pairIn       ! incoming baryon and antibaryon (order does not matter)
  ! * real, intent(in)                              :: time       ! current time
  ! OUTPUT
  ! * type(particle), intent(inout), dimension(:) :: finalState   ! Array of the final state particles
  ! * logical, intent(out) :: flag_anni                           ! .true. if annihilation successfull
  ! * integer, intent(out) :: HiEnergyType                        ! must be=-3 at the output (see docu of setKinematicsHiEnergy)
  !****************************************************************************
  subroutine annihilate_to_meson(pairIn,time,finalState,flag_anni,HiEnergyType)

    use particleDefinition, only: particle
    use RMF, only: getRMF_flag
    use densitymodule, only: Particle4Momentum_RMF

    type(particle), dimension(1:2), intent(in) ::  pairIn       ! incoming baryon and antibaryon (order does not matter)
    real, intent(in)                              :: time       ! current time
    type(particle), intent(inout), dimension(:) :: finalState   ! Array of the final state particles
    logical, intent(out) :: flag_anni                           ! .true. if annihilation successfull
    integer, intent(out) :: HiEnergyType                        ! must be=-3 at the output

    real, dimension(0:3) :: momentum_tot, momentum1, momentum2
    real :: srtS
    logical, parameter :: debugFlag=.false.

    if (.not.getRMF_flag()) then
       momentum_tot=pairIN(1)%momentum+pairIN(2)%momentum
    else
       call Particle4Momentum_RMF(pairIN(1),momentum1)
       call Particle4Momentum_RMF(pairIN(2),momentum2)
       momentum_tot=momentum1+momentum2
    end if

    srtS= momentum_tot(0)**2 - dot_product(momentum_tot(1:3),momentum_tot(1:3))

    if (srtS.le.0.) then
       if (debugFlag) write(*,*) 'Space-like total 4-momentum in annihilate_to_meson, s**2: ', srtS
       flag_anni=.false.
       return
    end if

    srtS= sqrt(srtS)

    finalState(1)%mass=srts
    finalState(1)%momentum=momentum_tot

    finalState(1)%productionTime = time
    finalState(1)%formationTime = time
    finalState(1)%scaleCS=1.
    finalState(1)%in_formation=.false.

    HiEnergyType=-3
    flag_anni=.true.

  end subroutine annihilate_to_meson

end module Annihilation
