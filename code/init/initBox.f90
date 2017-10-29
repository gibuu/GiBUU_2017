!******************************************************************************
!****m* /initBox
! NAME
! module initBox
! PURPOSE
! Initializes particles within a box
!******************************************************************************
module initBox

  implicit none

  private

  !****************************************************************************
  !****g* initBox/thermalInit
  ! SOURCE
  !
  logical, save :: thermalInit = .false.
  ! PURPOSE
  ! flag how to initialize
  !****************************************************************************

  !****************************************************************************
  !****g* initBox/nDens
  ! SOURCE
  !
  real, save :: nDens = 1.0
  ! PURPOSE
  ! particle density [fm^-3]
  !****************************************************************************

  !****************************************************************************
  !****g* initBox/ChargeSelection
  ! SOURCE
  !
  integer, save :: ChargeSelection = 0
  ! PURPOSE
  ! define the type of the charge selection:
  ! * 0: only pi0
  ! * 1: 50% pi+, 50% pi-
  ! * 2: 33% for +,0,-
  !****************************************************************************

  !****************************************************************************
  !****g* initBox/pInit
  ! SOURCE
  !
  real, save :: pInit = 0.5
  ! PURPOSE
  ! initial momentum of particles [GeV/c]
  !****************************************************************************

  !****************************************************************************
  !****g* initBox/BoostZ
  ! SOURCE
  !
  real, save :: BoostZ = 0.0
  ! PURPOSE
  ! additional boost for all particles in z-direction
  !****************************************************************************

  !****************************************************************************
  !****g* initBox/Temp
  ! SOURCE
  real, dimension(1:122), save :: Temp = 0.0
  ! PURPOSE
  ! for thermal init: temperature of every meson species in GeV,
  ! if larger than 0. otherwise this species is not initialized
  !****************************************************************************

  !****************************************************************************
  !****g* initBox/Fugacity
  ! SOURCE
  real, dimension(1:122), save :: Fugacity = 1.0
  ! PURPOSE
  ! for thermal init: fugacity of every hadron species.
  !****************************************************************************

  !****************************************************************************
  !****g* initBox/correctMovingBox
  ! SOURCE
  integer, save :: correctMovingBox = 1
  ! PURPOSE
  ! switch to indicate, whether a correction of the momenta after
  ! initialization should be done to enforce vanishing 3-momenta.
  ! possibilities are:
  ! * 0 : no correction
  ! * 1 : global correction
  ! * 2 : per ensemble correction
  !****************************************************************************

  public :: initializeBox


  real,dimension(101:122), save :: densMes = 0
  real,dimension(1:2,1:61), save :: densBar = 0 ! Baryon, Antibaryon

contains

  !****************************************************************************
  !****s* initBox/initInput
  ! NAME
  ! subroutine initInput
  ! PURPOSE
  ! Reads input out of jobcard. Namelist 'initBox'.
  !****************************************************************************
  subroutine initInput
    use output, only: Write_ReadingInput
    use callstack, only: traceBack
    use constants, only: mPi

    !**************************************************************************
    !****n* initBox/Box
    ! NAME
    ! NAMELIST Box
    ! PURPOSE
    ! Includes the input parameters:
    ! * thermalInit
    ! * nDens
    ! * ChargeSelection
    ! * pInit
    ! * BoostZ
    ! * Temp
    ! * Fugacity
    ! * correctMovingBox
    !**************************************************************************
    NAMELIST /Box/ thermalInit, &
         nDens,ChargeSelection,pInit,BoostZ, &
         Temp, Fugacity, &
         correctMovingBox

    integer :: ios, ID
    character(20), dimension(0:2), parameter :: cSel = (/ &
         'only pi0        ',&
         '50% pi+, 50% pi-',&
         '33% for +,0,-   ' /)
    real :: eDens,beta(1:2),sumDens, h

    call Write_ReadingInput('Box',0)
    rewind(5)
    read(5,nml=Box,iostat=ios)

    do ID=1,122
       if ((Temp(ID)/=0.0).or.(Fugacity(ID)/=1.0)) then
          write(*,'(" ID = ",i3,"  Temp = ",f5.3,"  Fugacity = ",f5.3)') &
               ID,Temp(ID),Fugacity(ID)
       end if
    end do

    call Write_ReadingInput('Box',0,ios)

    densMes = 0
    densBar = 0

    write(*,*)
    if (thermalInit) then
       write(*,*) '##### thermal init #####'

       do ID=101,122
          if ((ID>=114).and.(ID<=121)) cycle ! no charm
          if (Temp(ID) >= 0.01) then
             densMes(ID) = integrateMeson(ID, Temp(ID)) * Fugacity(ID)
!!$             write(*,'("  ",i3,"  ",f5.3,"  ",f5.3," -> ",1P,e13.4)') &
!!$                  ID, Temp(ID), Fugacity(ID), densMes(ID)
          end if
       end do

       sumDens = sum(densMes)
       write(*,*) 'total density of mesons   : ',sum(densMes)

       do ID=1,61
          if ((ID>=56).and.(ID<=61)) cycle ! no charm
          if (Temp(ID) >= 0.01) then
             h = integrateBaryon(ID, Temp(ID))
             densBar(1,ID) = h * Fugacity(ID)
             densBar(2,ID) = h / Fugacity(ID)

!!$             write(*,'("  ",i3,"  ",f5.3,"  ",f5.3," -> ",1P,2e13.4)') &
!!$                  ID, Temp(ID), Fugacity(ID), densBar(1:2,ID)
          end if
       end do

       sumDens = sumDens + sum(densBar)
       write(*,*) 'total density of baryons  : ',sum(densBar)
       write(*,*) 'total density of particles: ',sumDens
       write(*,*)

    else
       write(*,*) '### non-thermal init ###'

       select case (ChargeSelection)
       case (0:2)
          write(*,*) 'charge selection: ',ChargeSelection,"= ",&
               cSel(ChargeSelection)
       case default
          write(*,*) 'charge selection: ',ChargeSelection,"= WRONG!"
          call traceBack('wrong input value')
       end select
       write(*,'(A,F8.3,A)') ' particle density: nDens =',nDens,' fm^(-3)'
       write(*,'(A,F8.3,A)') ' initial momentum: pInit =',pInit,' GeV/c'

       write(*,*)
       eDens = sqrt(mPi**2+pInit**2)*nDens

       ! the temperature is given by
       !   eDens/(m*nDens) == 3/betam + K1(betam)/K2(betam)
       ! As approx, we can give K1/K2 ~ betam/2
       ! This is then solved for beta:
       ! (The resulting T is off by 5-10 MeV)
       beta(1) = eDens + sqrt(eDens**2-6*(mPi*nDens)**2)
       beta(2) = eDens - sqrt(eDens**2-6*(mPi*nDens)**2)
       beta = beta/(mPi**2*nDens)
       write(*,'(A,F8.3,A)') ' energy density:   eDens =',eDens,' GeV fm^(-3)'
       write(*,'(A,F6.2,A)') ' Estimate temperature: T = ',1000/beta(2),' MeV'
       write(*,'(A,F6.2,A)') ' Estimate fugacity   : -- to be done-- '

       write(*,*)
       write(*,'(A,F8.3,A)') ' additional boost: betaZ =',BoostZ

    end if

    write(*,*) 'correct moving box: ', correctMovingBox

    call Write_ReadingInput('Box',1)

  contains
    !**************************************************************************
    real function integrateMeson(ID, T)

      use constants, only: pi
      use particleProperties, only: hadron
      use mesonWidth, only: fullWidthMeson

      integer, intent(in) :: ID
      real, intent(in) :: T

      integer :: im, nm
      real, parameter :: massMax = 3.0 ! upper boundary of integration
      real :: mass0, gamma0
      real :: mmin, mmax, m
      real :: ymin, ymax, dy, y
      real :: gamma, spectral, intfac, deg, sum
      real, parameter :: dy0 = pi/200.

      integrateMeson = 0.0

      mass0  = hadron(id)%mass
      gamma0 = hadron(id)%width
      deg = (hadron(ID)%Spin*2+1)*(hadron(ID)%isoSpinTimes2+1)

      if (gamma0 < 1e-3) then
         if (massMax > mass0) then
            integrateMeson = BoltzmannN(mass0, T)*deg
         end if
         return
      end if

      mmin = hadron(id)%minmass
      mmax = massMax
      if (mmax < mmin) return ! integral is zero

      ymax = 2*atan(2*(mmax-mass0) / gamma0)
      ymin = 2*atan(2*(mmin-mass0) / gamma0)

      nm = max(int((ymax-ymin)/dy0),1)
      dy  = (ymax-ymin)/float(nm)

      sum = 0.0
      do im=1,nm
         y = ymin+(float(im)-0.5)*dy
         m = 0.5*tan(0.5*y)*gamma0 + mass0
         m = min(max(m,mmin),mmax)
         gamma = fullWidthMeson(id, m)
         spectral = 2./pi * m**2 * gamma / ((m**2-mass0**2)**2+m**2*gamma**2)
         intfac = gamma0 / ((m-mass0)**2+gamma0**2/4.)
         sum = sum + spectral/intfac * BoltzmannN(m, T)
      end do
      integrateMeson = sum * dy * deg

    end function integrateMeson

    !**************************************************************************
    real function integrateBaryon(ID, T)

      use constants, only: pi
      use particleProperties, only: hadron
      use BaryonWidth, only: fullWidthBaryon

      integer, intent(in) :: ID
      real, intent(in) :: T

      integer :: im, nm
      real, parameter :: massMax = 3.0 ! upper boundary of integration
      real :: mass0, gamma0
      real :: mmin, mmax, m
      real :: ymin, ymax, dy, y
      real :: gamma, spectral, intfac, deg, sum
      real, parameter :: dy0 = pi/200.

      integrateBaryon = 0.0

      mass0  = hadron(id)%mass
      gamma0 = hadron(id)%width
      deg = (hadron(ID)%Spin*2+1)*(hadron(ID)%isoSpinTimes2+1)

      if (gamma0 < 1e-3) then
         if (massMax > mass0) then
            integrateBaryon = BoltzmannN(mass0, T)*deg
         end if
         return
      end if

      mmin = hadron(id)%minmass
      mmax = massMax
      if (mmax < mmin) return ! integral is zero

      ymax = 2*atan(2*(mmax-mass0) / gamma0)
      ymin = 2*atan(2*(mmin-mass0) / gamma0)

      nm = max(int((ymax-ymin)/dy0),1)
      dy  = (ymax-ymin)/float(nm)

      sum = 0.0
      do im=1,nm
         y = ymin+(float(im)-0.5)*dy
         m = 0.5*tan(0.5*y)*gamma0 + mass0
         m = min(max(m,mmin),mmax)
         gamma = fullWidthBaryon(id, m)
         spectral = 2./pi * m**2 * gamma / ((m**2-mass0**2)**2+m**2*gamma**2)
         intfac = gamma0 / ((m-mass0)**2+gamma0**2/4.)
         sum = sum + spectral/intfac * BoltzmannN(m, T)
      end do
      integrateBaryon = sum * dy * deg

    end function integrateBaryon

    !**************************************************************************
    pure real function BoltzmannN(mass, T)
      use constants, only: pi
      use besselK, only: BesselK2

      real, intent(in) :: mass, T

      real, parameter :: fak = 4*pi/(2*pi*0.197)**3

      if (T==0.0) then
         BoltzmannN = 0.0
      else
         BoltzmannN = fak * BesselK2(mass/T) * mass**2 * T
      end if

      return
    end function BoltzmannN

  end subroutine initInput

  !****************************************************************************
  !****s* initBox/initializeBox
  ! NAME
  ! subroutine initializeBox(part)
  ! PURPOSE
  ! Initialize nucleons in a box
  !****************************************************************************
  subroutine initializeBox(part)

    use particleDefinition
    use IdTable, only: pion
    use constants, only: mPi
    use densityModule, only: gridsize,get_densitySwitch
    use output, only: Write_InitStatus
    use callstack, only: traceBack
    use insertion, only: GarbageCollection
    use collisionNumbering, only: real_firstnumbering


    type(particle), dimension(:,:),intent(inOut),TARGET :: part

    integer :: dummy, iEns, iPart, nEns
    integer, dimension(0:5) :: nCharge

    call Write_InitStatus('Box',0)
    call initInput

    dummy = get_densitySwitch() ! force density module to be initialized
    nEns = size(part(:,1))

    write(*,'(A,3F9.3,A)') ' Gridsize = (',gridsize(:),') fm'
    write(*,'(A,1F9.3,A)') ' Size of Box= (2*Gridsize)^3 = ',8.*gridsize(1)*gridsize(2)*gridsize(3),' fm^3'
    write(*,*) 'Number Ensembles:              nEns =', nEns
    write(*,*) 'Number Particles per Ensemble: nPart=', size(part(1,:))

    if (thermalInit) then
!       call doThermal
       call doThermalGlobal
    else
       call doNonThermal
    end if

    call GarbageCollection(part)

    select case (correctMovingBox)
    case (0)
       ! nothing to do
    case (1)
       call doCorrectMovingBoxGlobal
    case (2)
       call doCorrectMovingBoxEnsemble
    case default
       write(*,*) 'correctMovingBox = ',correctMovingBox
       call traceback('wrong value')
    end select

    call calcSumGlobal

    call Write_InitStatus('Box',1)

  contains
    !**************************************************************************
    !****is* initializeBox/doThermal
    ! NAME
    ! subroutine doThermal
    ! PURPOSE
    ! Do the thermal init
    !**************************************************************************
    subroutine doThermal

      integer :: nPart
      integer :: ID

      nPart = 0

      write(*,*) 'Mesons:'
      do ID=101,122
         call doThermalPart(ID, .false., .true., densMes(ID), nPart)
      end do

      write(*,*)
      write(*,*) 'Baryons:'
      do ID=1,61
         call doThermalPart(ID, .false., .false., densBar(1,ID), nPart)
         call doThermalPart(ID, .true.,  .false., densBar(2,ID), nPart)
      end do

      if (nPart==0) &
           call Traceback('ERROR: density too low, no particles in box!!!')

      write(*,*)
      write(*,'(" Number of ALL per ensemble:         ",i7)') &
           nPart
      write(*,'(" Size          per ensemble:         ",i7)') &
           size(part(1,:))

    end subroutine doThermal

    !**************************************************************************
    subroutine doThermalPart(ID, anti, isMes, dens, nPart)
      use randomMaxwell, only: initMaxwell
      use AssignMassMC, only: AssignMass_1_Therm
      use mediumDefinition
      use particleProperties, only: hadron, validCharge

      integer, intent(in) :: ID
      logical, intent(in) :: anti, isMes
      real, intent(in) :: dens
      integer, intent(inOut) :: nPart

      type(particle), POINTER :: pPart
      integer :: nPart0, nPart1
      real :: mass
      real, dimension(0:3) :: momLRF = (/0,0,0,0/)

      nPart0 = nPart

      nPart1 = NINT(8.*gridsize(1)*gridsize(2)*gridsize(3)*dens)
      if (anti) then
         write(*,'(" Number of ",i4," per ensemble:         ",i7,1P,e13.4)') &
              -ID,nPart1,dens
      else
         write(*,'(" Number of ",i4," per ensemble:         ",i7,1P,e13.4)') &
              ID,nPart1,dens
      end if

      if (nPart1 < 1) return

      if (nPart0+nPart1 > size(part(1,:))) then
         call traceback('particle vector too small!')
      end if

      if (hadron(ID)%width < 1e-03) then
         call initMaxwell(hadron(ID)%mass, temp(ID))
      end if

      do iEns=1,nEns

         select case (hadron(ID)%isoSpinTimes2+1)
         case (1)
            if (hadron(ID)%strangeness /= -3) then
               call prepareChooseCharge1(nPart1) ! only neutral
            else
               if (.not.anti) then
                  call prepareChooseCharge1minus(nPart1) ! only -1
               else
                  call prepareChooseCharge1plus(nPart1) ! only +1
               end if
            end if

         case (2)
            if (hadron(ID)%charm.ne.0) then
               call traceback('charm not to be initialized!')
            end if
            if (isMes) then
               if (hadron(ID)%strangeness > 0) then
                  call prepareChooseCharge2plus(nPart1) ! 0,+1
               else
                  call prepareChooseCharge2minus(nPart1) ! -1,0
               end if
            else
               if (anti .neqv. (hadron(ID)%strangeness == 0)) then
                  call prepareChooseCharge2plus(nPart1) ! 0,+1
               else
                  call prepareChooseCharge2minus(nPart1) ! -1,0
               end if
            end if

         case (3)
            call prepareChooseCharge3(nPart1) ! -1,0,+1

         case (4)
            if (.not.anti) then
               call prepareChooseCharge4plus(nPart1) ! -1,0,+1,+2
            else
               call prepareChooseCharge4minus(nPart1) ! -2,-1,0,+1
            end if

         case default
            write(*,*) '(2I+1)=',hadron(ID)%isoSpinTimes2+1
            call traceback('Isospin not yet implemented!')

         end select


         do iPart=nPart0+1,nPart0+nPart1
            pPart => part(iEns,iPart)
            call setToDefault(pPart)
            pPart%antiparticle = anti
            call setNumber(pPart) ! give it a unique number

            pPart%event = real_firstnumbering()

            pPart%ID = ID
            if (hadron(ID)%width < 1e-03) then
               pPart%mass = hadron(ID)%mass
            else
               call AssignMass_1_Therm(ID, temp(ID), 3.0, momLRF,vacuum, mass)
               pPart%mass = mass
               call initMaxwell(mass, temp(ID)) ! re-init with new mass
            end if

            pPart%charge = chooseCharge()
            pPart%position = choosePosition()
            pPart%momentum = chooseMomentumTherm(pPart%mass)

            if (.not.validCharge(pPart)) then
               write(*,*) 'ID,anti,charge = ',pPart%ID,&
                    pPart%antiparticle, pPart%charge
               call Traceback('invalid charge')
            end if

         end do ! iPart
      end do ! iEns

      nPart = nPart0+nPart1

    end subroutine doThermalPart

    !**************************************************************************
    !****is* initializeBox/doThermalGlobal
    ! NAME
    ! subroutine doThermalGlobal
    ! PURPOSE
    ! Do the thermal init
    !**************************************************************************
    subroutine doThermalGlobal

      integer, dimension(1:nEns) :: nPart
      integer :: ID

      nPart = 0

      write(*,*)
      write(*,*) 'Mesons:'
      do ID=101,122
         call doThermalPartGlobal(ID, .false., .true., densMes(ID), nPart)
      end do

      write(*,*)
      write(*,*) 'Baryons:'
      do ID=1,61
         call doThermalPartGlobal(ID, .false., .false., densBar(1,ID), nPart)
         call doThermalPartGlobal(ID, .true.,  .false., densBar(2,ID), nPart)
      end do

      if (sum(nPart)==0) &
           call Traceback('ERROR: density too low, no particles in box!!!')

      write(*,*)
      write(*,'(" Number of ALL per ensemble:         ",f13.4)') &
           sum(nPart)/real(nEns)
      write(*,'(" Size          per ensemble:         ",f9.0)') &
           1.0*size(part(1,:))

    end subroutine doThermalGlobal

    !**************************************************************************
    subroutine doThermalPartGlobal(ID, anti, isMes, dens, nPart)
      use random, only: rn
      use randomMaxwell, only: initMaxwell
      use AssignMassMC, only: AssignMass_1_Therm
      use mediumDefinition
      use particleProperties, only: hadron, validCharge

      integer, intent(in) :: ID
      logical, intent(in) :: anti, isMes
      real, intent(in) :: dens
      integer, dimension(:), intent(inOut) :: nPart

      type(particle), POINTER :: pPart
      integer :: nPart1, iPart1
      real :: mass
      real, dimension(0:3) :: momLRF = (/0,0,0,0/)

      nPart1 = NINT(8.*gridsize(1)*gridsize(2)*gridsize(3)*dens*nEns)
      if (nPart1 < 1) return

      if (anti) then
         write(*,'("   Number of ",i4," in box:               ",i7,1P,e13.4)') &
              -ID,nPart1,dens
      else
         write(*,'("   Number of ",i4," in box:               ",i7,1P,e13.4)') &
              ID,nPart1,dens
      end if

      if (sum(nPart)+nPart1 > nEns*size(part(1,:))) then
         call traceback('particle vector too small!')
      end if

      if (hadron(ID)%width < 1e-03) then
         call initMaxwell(hadron(ID)%mass, temp(ID))
      end if

      select case (hadron(ID)%isoSpinTimes2+1)
      case (1)
         if (hadron(ID)%strangeness /= -3) then
            call prepareChooseCharge1(nPart1) ! only neutral
         else
            if (.not.anti) then
               call prepareChooseCharge1minus(nPart1) ! only -1
            else
               call prepareChooseCharge1plus(nPart1) ! only +1
            end if
         end if

      case (2)
         if (hadron(ID)%charm.ne.0) then
            call traceback('charm not to be initialized!')
         end if
         if (isMes) then
            if (hadron(ID)%strangeness > 0) then
               call prepareChooseCharge2plus(nPart1) ! 0,+1
            else
               call prepareChooseCharge2minus(nPart1) ! -1,0
            end if
         else
            if (anti .neqv. (hadron(ID)%strangeness == 0)) then
               call prepareChooseCharge2plus(nPart1) ! 0,+1
            else
               call prepareChooseCharge2minus(nPart1) ! -1,0
            end if
         end if

      case (3)
         call prepareChooseCharge3(nPart1) ! -1,0,+1

      case (4)
         if (.not.anti) then
            call prepareChooseCharge4plus(nPart1) ! -1,0,+1,+2
         else
            call prepareChooseCharge4minus(nPart1) ! -2,-1,0,+1
         end if

      case default
         write(*,*) '(2I+1)=',hadron(ID)%isoSpinTimes2+1
         call traceback('Isospin not yet implemented!')

      end select


      do iPart1=1,nPart1
         do
            iEns = 1 + int(nEns*rn())
            if (nPart(iEns) < size(part(1,:))) exit
         end do
         nPart(iEns) = nPart(iEns) + 1
         iPart = nPart(iEns)
         pPart => part(iEns,iPart)

         call setToDefault(pPart)
         pPart%antiparticle = anti
         call setNumber(pPart) ! give it a unique number

         pPart%event = real_firstnumbering()

         pPart%ID = ID
         if (hadron(ID)%width < 1e-03) then
            pPart%mass = hadron(ID)%mass
         else
            call AssignMass_1_Therm(ID, temp(ID), 3.0, momLRF,vacuum, mass)
            pPart%mass = mass
            call initMaxwell(mass, temp(ID)) ! re-init with new mass
         end if

         pPart%charge = chooseCharge()
         pPart%position = choosePosition()
         pPart%momentum = chooseMomentumTherm(pPart%mass)

         if (.not.validCharge(pPart)) then
            write(*,*) 'ID,anti,charge = ',pPart%ID,&
                 pPart%antiparticle, pPart%charge
            call Traceback('invalid charge')
         end if

      end do ! iPart

    end subroutine doThermalPartGlobal


    !**************************************************************************
    !****is* initializeBox/doNonThermal
    ! NAME
    ! subroutine doNonThermal
    ! PURPOSE
    ! Do the non-thermal init
    !**************************************************************************
    subroutine doNonThermal
      use lorentzTrafo, only: lorentz

      integer :: nTot
      real, dimension(0:3) :: pSum=0
      real, dimension(1:3) :: betaVec=0
      type(particle), POINTER :: pPart

      nTot = NINT(8.*gridsize(1)*gridsize(2)*gridsize(3)*nDens)
      write(*,*) ' Number of pions per ensemble:        ', nTot

      if (nTot > size(part(1,:))) then
         call traceback('particle vector too small!')
      end if

      do iEns=1,nEns
         pSum = 0

         select case (ChargeSelection)
         case (0)
            call prepareChooseCharge1(nTot)
         case (1)
            call prepareChooseCharge2zero(nTot)
         case (2)
            call prepareChooseCharge3(nTot)
         end select

         do iPart=1,nTot
            pPart => part(iEns,iPart)
            call setToDefault(pPart)
            call setNumber(pPart) ! give it a unique number

            pPart%event = real_firstnumbering()
            pPart%ID = pion
            pPart%mass = mPi
            pPart%charge = chooseCharge()
            pPart%position = choosePosition()
            pPart%momentum = chooseMomentum(mPi)

            pSum = pSum + pPart%momentum(0:3)

         end do

      end do


      if (BoostZ>0) then
         betaVec = (/ 0.0, 0.0, -BoostZ /)
         pSum = 0
         do iEns=1,nEns
            do iPart=1,nTot
               call lorentz(betaVec, part(iEns,iPart)%momentum )
               pSum = pSum + part(iEns,iPart)%momentum(0:3)
            end do
         end do
         pSum = pSum/(nEns*nTot)
         write(*,*) 'momentum/particle: ',pSum
      end if


    end subroutine doNonThermal

    !**************************************************************************
    !****is* initializeBox/prepareChooseCharge1
    ! NAME
    ! subroutine prepareChooseCharge1(nTot)
    ! PURPOSE
    ! prepare the random selection of the charges:
    ! (-2,-1,0,+1,+2) = ( - , - , 100%, - , - )
    !**************************************************************************
    subroutine prepareChooseCharge1(nTot)
      integer, intent(in) :: nTot

      nCharge = (/ nTot, 0, 0, nTot, 0, 0 /)

    end subroutine prepareChooseCharge1

    !**************************************************************************
    !****is* initializeBox/prepareChooseCharge1minus
    ! NAME
    ! subroutine prepareChooseCharge1minus(nTot)
    ! PURPOSE
    ! prepare the random selection of the charges:
    ! (-2,-1,0,+1,+2) = ( - , 100%, - , - , - )
    !**************************************************************************
    subroutine prepareChooseCharge1minus(nTot)
      integer, intent(in) :: nTot

      nCharge = (/ nTot, 0, nTot, 0, 0, 0 /)

    end subroutine prepareChooseCharge1minus

    !**************************************************************************
    !****is* initializeBox/prepareChooseCharge1plus
    ! NAME
    ! subroutine prepareChooseCharge1plus(nTot)
    ! PURPOSE
    ! prepare the random selection of the charges:
    ! (-2,-1,0,+1,+2) = ( - , 100%, - , - , - )
    !**************************************************************************
    subroutine prepareChooseCharge1plus(nTot)
      integer, intent(in) :: nTot

      nCharge = (/ nTot, 0, 0, 0, nTot, 0 /)

    end subroutine prepareChooseCharge1plus

    !**************************************************************************
    !****is* initializeBox/prepareChooseCharge2minus
    ! NAME
    ! subroutine prepareChooseCharge2minus(nTot)
    ! PURPOSE
    ! prepare the random selection of the charges:
    ! (-2,-1,0,+1,+2) = ( - , 50%, 50%, - , - )
    ! If the total number of particles is not even, the number of
    ! neutral particles is increased by one
    !**************************************************************************
    subroutine prepareChooseCharge2minus(nTot)
      integer, intent(in) :: nTot
      integer :: nTot2

      nTot2 = int(nTot/2)
      select case (nTot-2*nTot2)
      case (0)
         nCharge = (/ nTot, 0, nTot2, nTot2, 0, 0 /)
      case (1)
         nCharge = (/ nTot, 0, nTot2, nTot2+1, 0, 0 /)
      case default
         write(*,*) nTot,nTot2
         call traceback('fishy nCharge')
      end select

    end subroutine prepareChooseCharge2minus

    !**************************************************************************
    !****is* initializeBox/prepareChooseCharge2zero
    ! NAME
    ! subroutine prepareChooseCharge2zero(nTot)
    ! PURPOSE
    ! prepare the random selection of the charges:
    ! (-2,-1,0,+1,+2) = ( -, 50%, - , 50%, - )
    ! If the total number of particles is not even, the number of
    ! negative particles is increased by one, thus the symmetry is
    ! broken (but who cares?)
    !**************************************************************************
    subroutine prepareChooseCharge2zero(nTot)
      integer, intent(in) :: nTot
      integer :: nTot2

      nTot2 = int(nTot/2)
      select case (nTot-2*nTot2)
      case (0)
         nCharge = (/ nTot, 0, nTot2, 0, nTot2, 0/)
      case (1)
         ! this breaks the symmetry, because we have more pi- than pi+,
         ! but who cares? ;)
         nCharge = (/ nTot, 0, nTot2+1, 0, nTot2, 0 /)
      case default
         write(*,*) nTot,nTot2
         call traceback('fishy nCharge')
      end select

    end subroutine prepareChooseCharge2zero

    !**************************************************************************
    !****is* initializeBox/prepareChooseCharge2plus
    ! NAME
    ! subroutine prepareChooseCharge2plus(nTot)
    ! PURPOSE
    ! prepare the random selection of the charges:
    ! (-2,-1,0,+1,+2) = ( -, - , 50%, 50%, - )
    ! If the total number of particles is not even, the number of
    ! neutral particles is increased by one
    !**************************************************************************
    subroutine prepareChooseCharge2plus(nTot)
      integer, intent(in) :: nTot
      integer :: nTot2

      nTot2 = int(nTot/2)
      select case (nTot-2*nTot2)
      case (0)
         nCharge = (/ nTot, 0, 0, nTot2, nTot2, 0 /)
      case (1)
         nCharge = (/ nTot, 0, 0, nTot2+1, nTot2, 0 /)
      case default
         write(*,*) nTot,nTot2
         call traceback('fishy nCharge')
      end select

    end subroutine prepareChooseCharge2plus

    !**************************************************************************
    !****is* initializeBox/prepareChooseCharge3
    ! NAME
    ! subroutine prepareChooseCharge3(nTot)
    ! PURPOSE
    ! prepare the random selection of the charges:
    ! (-2,-1,0,+1,+2) = ( - , 33%, 33%, 33%, - )
    ! If the total number of particles is not dividable by 3, the
    ! number of neutral or the number of both charged states is
    ! increased by one in order to get the total number of particles,
    ! but always in charge neutral sum.
    !**************************************************************************
    subroutine prepareChooseCharge3(nTot)
      integer, intent(in) :: nTot
      integer :: nTot3
      nTot3 = int(nTot/3)
      select case (nTot-3*nTot3)
      case (0)
         nCharge = (/ nTot, 0, nTot3, nTot3, nTot3, 0 /)
      case (1)
         nCharge = (/ nTot, 0, nTot3, nTot3+1, nTot3, 0 /)
      case (2)
         nCharge = (/ nTot, 0, nTot3+1, nTot3, nTot3+1, 0 /)
      case default
         write(*,*) nTot,nTot3
         call traceback('fishy nCharge')
      end select
    end subroutine prepareChooseCharge3

    !**************************************************************************
    !****is* initializeBox/prepareChooseCharge4minus
    ! NAME
    ! subroutine prepareChooseCharge4minus(nTot)
    ! PURPOSE
    ! prepare the random selection of the charges:
    ! (-2,-1,0,+1,+2) = ( 25%, 25%, 25%, 25%, - )
    ! If the total number of particles is not a multiple of 4, the number
    ! of particles is increased such that the overall charge is not
    ! modified
    !**************************************************************************
    subroutine prepareChooseCharge4minus(nTot)
      integer, intent(in) :: nTot
      integer :: nTot2

      nTot2 = int(nTot/4)
      select case (nTot-4*nTot2)
      case (0)
         nCharge = (/ nTot, nTot2, nTot2, nTot2, nTot2, 0 /)
      case (1)
         nCharge = (/ nTot, nTot2, nTot2, nTot2+1, nTot2, 0 /)
      case (2)
         nCharge = (/ nTot, nTot2, nTot2+1, nTot2, nTot2+1, 0 /)
      case (3)
         nCharge = (/ nTot, nTot2, nTot2+1, nTot2+1, nTot2+1, 0 /)
      case default
         write(*,*) nTot,nTot2
         call traceback('fishy nCharge')
      end select

    end subroutine prepareChooseCharge4minus

    !**************************************************************************
    !****is* initializeBox/prepareChooseCharge4plus
    ! NAME
    ! subroutine prepareChooseCharge4plus(nTot)
    ! PURPOSE
    ! prepare the random selection of the charges:
    ! (-2,-1,0,+1,+2) = ( -, 50%, 50%, 50%, 50% )
    ! If the total number of particles is not a multiple of 4, the number
    ! of particles is increased such that the overall charge is not
    ! modified
    !**************************************************************************
    subroutine prepareChooseCharge4plus(nTot)
      integer, intent(in) :: nTot
      integer :: nTot2

      nTot2 = int(nTot/4)
      select case (nTot-4*nTot2)
      case (0)
         nCharge = (/ nTot, 0, nTot2, nTot2, nTot2, nTot2 /)
      case (1)
         nCharge = (/ nTot, 0, nTot2, nTot2+1, nTot2, nTot2 /)
      case (2)
         nCharge = (/ nTot, 0, nTot2+1, nTot2, nTot2+1, nTot2 /)
      case (3)
         nCharge = (/ nTot, 0, nTot2+1, nTot2+1, nTot2+1, nTot2 /)
      case default
         write(*,*) nTot,nTot2
         call traceback('fishy nCharge')
      end select

    end subroutine prepareChooseCharge4plus



    !**************************************************************************
    !****if* initializeBox/chooseCharge
    ! NAME
    ! integer function chooseCharge()
    ! PURPOSE
    ! Return a charge randomly selected according the numbers
    ! given in the array nCharge
    !**************************************************************************
    integer function chooseCharge()
      use random, only: rn
      real :: r

      r = rn() * nCharge(0)
      if ( r < nCharge(1) ) then
         chooseCharge = -2
         nCharge(1) = nCharge(1)-1
      else if ( r < sum(nCharge(1:2)) ) then
         chooseCharge = -1
         nCharge(2) = nCharge(2)-1
      else if ( r < sum(nCharge(1:3)) ) then
         chooseCharge = 0
         nCharge(3) = nCharge(3)-1
      else if ( r < sum(nCharge(1:4)) ) then
         chooseCharge = 1
         nCharge(4) = nCharge(4)-1
      else
         chooseCharge = 2
         nCharge(5) = nCharge(5)-1
      end if
      nCharge(0) = nCharge(0)-1

    end function chooseCharge

    !**************************************************************************
    !****if* initializeBox/choosePosition
    ! NAME
    ! function choosePosition()
    ! PURPOSE
    ! Set the position of the particle randomly in the box
    ! RESULT
    ! * real, dimension(1:3) : function value
    !**************************************************************************
    function choosePosition()
      use random, only: rn

      real, dimension(1:3) :: choosePosition

      choosePosition(1)=(1.-2.*rn())*gridSize(1)
      choosePosition(2)=(1.-2.*rn())*gridSize(2)
      choosePosition(3)=(1.-2.*rn())*gridSize(3)

    end function choosePosition

    !**************************************************************************
    !****if* initializeBox/chooseMomentum
    ! NAME
    ! subroutine chooseMomentum(mass)
    ! PURPOSE
    ! Return a momentum of the particle within the non-thermal init.
    ! The direction is choosen isotropically, while the absolute value
    ! is fixed.
    !**************************************************************************
    function chooseMomentum(mass)
      use random, only: rnOmega
      real, dimension(0:3) :: chooseMomentum
      real, intent(in) :: mass

      chooseMomentum(1:3) = rnOmega() * pInit
      chooseMomentum(0) = sqrt(mass**2+pInit**2)

    end function chooseMomentum

    !**************************************************************************
    !****if* initializeBox/chooseMomentumTherm
    ! NAME
    ! subroutine chooseMomentumTherm(mass)
    ! PURPOSE
    ! Set the momentum of the particle within the thermal init.
    ! The direction is choosen isotropically. The absolute value of the
    ! momentum is choosen according a relativistic Maxwell distribution.
    !
    ! The Maxwell distribution has to be initialized before!
    !**************************************************************************
    function chooseMomentumTherm(mass)
      use random, only: rnOmega
      use randomMaxwell, only: rnMaxwell

      real, dimension(0:3) :: chooseMomentumTherm
      real, intent(in) :: mass

      real :: p

      p = rnMaxwell()
      chooseMomentumTherm(1:3) = rnOmega() * p
      chooseMomentumTherm(0) = sqrt(mass**2+p**2)
    end function chooseMomentumTherm

    !**************************************************************************
    !****is* initializeBox/calcSumGlobal
    ! NAME
    ! subroutine calcSumGlobal
    ! PURPOSE
    ! calculate the sums of many observables from all particles
    !**************************************************************************
    subroutine calcSumGlobal

      use output
      use IdTable, only: isMeson, isBaryon
      use particleProperties, only: hadron

      type(particle), POINTER :: pPart
      real, dimension(0:3, 0:2) :: sumMom=0
      integer, dimension(0:3, 0:2) :: sumC=0 ! N,B,S,Q; total,meson,baryon
      integer, dimension(0:3, 0:2) :: sumC2=0

      integer :: ID,iH
      real :: Vol,mulFak
      integer, dimension(0:3) :: C = (/ 1, 0, 0, 0 /)

      Vol = 8.*gridsize(1)*gridsize(2)*gridsize(3)
      mulFak = 1.0/(nEns*Vol)

      do iEns=1,nEns
         do iPart=1,size(part(1,:))
            pPart => part(iEns,iPart)
            if (pPart%Id <  0) exit
            if (pPart%Id <= 0) cycle
            ID = pPart%Id

            if (isMeson(ID)) then
               iH = 1
            else if (isBaryon(ID)) then
               iH = 2
            else
               cycle
            end if

            C(1) = 0
            if (isBaryon(ID)) C(1) = 1
            C(2) = hadron(ID)%strangeness
            C(3) = pPart%charge

            if (pPart%antiparticle) then
               C(1:3) = -C(1:3)
            end if

            sumMom(:,0) = sumMom(:,0) + pPart%momentum(0:3)
            sumC(:,0) = sumC(:,0) + C
            sumC2(:,0) = sumC2(:,0) + C**2

            sumMom(:,iH) = sumMom(:,iH) + pPart%momentum(0:3)
            sumC(:,iH) = sumC(:,iH) + C
            sumC2(:,iH) = sumC2(:,iH) + C**2


         end do
      end do

      write(*,*)
      write(*,subchapter) 'Initial Densities'
      write(*,'(A10,4A13)') '','total','meson','baryon'
      write(*,'(A10,1P,4e13.4)') 'n =',sumC(0,:)*mulFak
      write(*,'(A10,1P,4e13.4)') 'e =',sumMom(0,:)*mulFak
      write(*,*)
!!$      write(*,'(A10,1P,4e13.4)') 'px =',sumMom(1,:)*mulFak
!!$      write(*,'(A10,1P,4e13.4)') 'py =',sumMom(2,:)*mulFak
!!$      write(*,'(A10,1P,4e13.4)') 'pz =',sumMom(3,:)*mulFak
!!$      write(*,*)
      write(*,'(A10,1P,4e13.4)') 'b =',sumC(1,:)*mulFak
      write(*,'(A10,1P,4e13.4)') 's =',sumC(2,:)*mulFak
      write(*,'(A10,1P,4e13.4)') 'q =',sumC(3,:)*mulFak
      write(*,*)
      write(*,'(A10,1P,4e13.4)') 'sqrt(b^2) =',sqrt(real(sumC2(1,:)))*mulFak
      write(*,'(A10,1P,4e13.4)') 'sqrt(s^2) =',sqrt(real(sumC2(2,:)))*mulFak
      write(*,'(A10,1P,4e13.4)') 'sqrt(q^2) =',sqrt(real(sumC2(3,:)))*mulFak
      write(*,'(79("-"))')
      write(*,*)

    end subroutine calcSumGlobal

    !**************************************************************************
    !****is* initializeBox/doCorrectMovingBoxGlobal
    ! NAME
    ! subroutine doCorrectMovingBoxGlobal
    ! PURPOSE
    ! correct momenta on a global basis
    !**************************************************************************
    subroutine doCorrectMovingBoxGlobal

      type(particle), POINTER :: pPart
      real, dimension(0:3) :: pSum1=0, pSum2=0
      integer :: nTot = 0

      do iEns=1,nEns
         do iPart=1,size(part(1,:))
            pPart => part(iEns,iPart)
            if (pPart%Id <  0) exit
            if (pPart%Id <= 0) cycle

            pSum1 = pSum1 + pPart%momentum(0:3)
            nTot = nTot+1
         end do
      end do

      pSum1 = pSum1/nTot

      do iEns=1,nEns
         do iPart=1,size(part(1,:))
            pPart => part(iEns,iPart)
            if (pPart%Id <  0) exit
            if (pPart%Id <= 0) cycle

            pPart%momentum(1:3) = pPart%momentum(1:3) - pSum1(1:3)
            pPart%momentum(0) = sqrt(pPart%mass**2 + sum(pPart%momentum(1:3)**2))

            pSum2 = pSum2 + pPart%momentum(0:3)
         end do
      end do

      pSum2 = pSum2/nTot

      write(*,'(A,1P,4e13.4)') ' before moving box correction: ',pSum1
      write(*,'(A,1P,4e13.4)') ' after  moving box correction: ',pSum2

    end subroutine doCorrectMovingBoxGlobal

    !**************************************************************************
    !****is* initializeBox/doCorrectMovingBoxEnsemble
    ! NAME
    ! subroutine doCorrectMovingBoxEnsemble
    ! PURPOSE
    ! correct momenta for every single ensemble
    !**************************************************************************
    subroutine doCorrectMovingBoxEnsemble

      type(particle), POINTER :: pPart
      real, dimension(0:3) :: pSum1=0, pSum2=0, pSum1a
      integer :: nTot = 0, nTota

      do iEns=1,nEns
         pSum1a = 0
         nTota = 0
         do iPart=1,size(part(1,:))
            pPart => part(iEns,iPart)
            if (pPart%Id <  0) exit
            if (pPart%Id <= 0) cycle

            pSum1a = pSum1a + pPart%momentum(0:3)
            nTota = nTota+1
         end do

         pSum1 = pSum1 + pSum1a
         nTot = nTot + nTota

         pSum1a = pSum1a/nTota

         do iPart=1,size(part(1,:))
            pPart => part(iEns,iPart)
            if (pPart%Id <  0) exit
            if (pPart%Id <= 0) cycle

            pPart%momentum(1:3) = pPart%momentum(1:3) - pSum1a(1:3)
            pPart%momentum(0) = sqrt(pPart%mass**2 + sum(pPart%momentum(1:3)**2))

            pSum2 = pSum2 + pPart%momentum(0:3)
         end do

      end do

      pSum1 = pSum1/nTot
      pSum2 = pSum2/nTot

      write(*,'(A,1P,4e13.4)') ' before moving box correction: ',pSum1
      write(*,'(A,1P,4e13.4)') ' after  moving box correction: ',pSum2

    end subroutine doCorrectMovingBoxEnsemble


  end subroutine initializeBox

end module initBox
