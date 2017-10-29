!******************************************************************************
!****m* /mesonWidth
! NAME
! module mesonWidth
! PURPOSE
! When this module is initialized then all information for the VACUUM width
! is once calculated for all meson resonances and then stored into the
! field gammaField, which is of type gammaFieldType. This is done by
! initWidth. Afterwards this field is used to return full and
! partial width of the meson resonances in the vacuum by the subroutine
! "partialWidthMeson, fullWidthMeson"
! USES
! MesonWidthVacuum
!******************************************************************************
module mesonWidth

  use IdTable, only: pion, nmes
  use particleProperties, only: nDecays

  implicit none
  private

  ! Parameters to store the partial widths depending on mass  in the range
  ! (minimalmass, minimalMass+maxIndex*deltaMass)
  integer, parameter :: maxIndexMass=2000
  real, parameter :: deltaMass=0.004

  ! The type used for storage:
  type tGammaField
     real :: gammatotal                ! total decay rate
     real, dimension(nDecays) :: ratio ! ratio of different  decay channels
  end type tGammaField

  ! Field which holds all the information for concerning  the vacuum width,
  ! initialized in "initWidth"
  ! First Index: ID of Resonance
  ! Second Index: Mass Index in the range (0,maxIndexMass)
  type(tGammaField), dimension(:,:), allocatable, save :: gammaField

  ! Flag to check wether this module is initialized by initWidth
  logical, save :: initFlag = .true.

  ! Switching debugging infos off and on
  logical, parameter :: debug=.false.

  ! Array to hold the integrals over the spectral functions
  real, save :: arrSpectralIntegral(pion:pion+nmes-1) = 0.0

  public :: partialWidthMeson
  public :: fullWidthMeson
  public :: decayWidthMeson
  public :: cleanUp
  public :: GetMinMaxMass
  public :: GetMaxQ
  public :: GetSpectralIntegral
  public :: calcSpectralIntegral

contains


  !****************************************************************************
  !****f* mesonWidth/partialWidthMeson
  ! NAME
  ! function partialWidthMeson (ID, mass, IDout1, IDout2, IDout3, iQ1, iQ2, iQ3)
  ! PURPOSE
  ! This function calculates the partial width (energy dependent) of all
  ! meson resonances.
  ! INPUTS
  ! * integer :: ID -- ID of resonance
  ! * real :: mass -- sqrt(p_mu p^mu) = mass of the resonance (offshell)
  ! * integer :: IDout1, IDout2 -- IDs of decay products
  !   (selecting channel of interest)
  ! * integer, OPTIONAL :: IDout3 -- ID of third decay product
  ! * integer, OPTIONAL :: iQ1, iQ2, iQ3 -- Charges of decay products
  !   (only relevant for 3-body decays)
  !****************************************************************************
  real function partialWidthMeson (ID, mass, IDout1, IDout2, IDout3, iQ1, iQ2, iQ3)
    use DecayChannels, only: Decay2bodyMeson
    use particleProperties, only: hadron, nDecays
    use idTable, only: pion, eta, photon, isMeson
    use CallStack, only: TRACEBACK

    real,   intent(in) :: mass
    integer,intent(in) :: ID
    integer,intent(in) :: IDout1, IDout2

    ! Only for three-body-channels necessary
    integer, intent(in),optional :: IDout3
    integer, intent(in),optional :: iQ1,iQ2,iQ3

    real  :: width
    integer :: down
    real    :: mass_down,weight
    integer :: i, dId, dId_wished

    if (initFlag) call initWidth

    !**************************************************************************
    ! (1) Check Input :
    if (.not.IsMeson(ID)) call TRACEBACK()
    if (.not.((IsMeson(IDout1).or.(IDout1.eq.photon)))) call TRACEBACK()
    if (.not.((IsMeson(IDout2).or.(IDout2.eq.photon)))) call TRACEBACK()
    if (present(IDout3)) then
       if (.not.((IsMeson(IDout3).or.(IDout3.eq.photon)))) call TRACEBACK()
       if (.not.(present(iQ1).and.present(iQ2).and.present(iQ3))) call TRACEBACK()
    end if


    !**************************************************************************
    down = floor((mass-hadron(ID)%minmass)/deltaMass)
    mass_down = hadron(ID)%minmass+float(down)*deltaMass
    weight = (mass-mass_down)/deltaMass ! weight for the interpolation


    if (mass <= hadron(ID)%minmass) then
       ! If mass is lower than minimalMas than gamma should be zero
       partialwidthMeson=0.

    else if (down >= 0) then
       ! Assume constant partial width at very high mass:
       if (down >= maxIndexMass) down = maxIndexMass-1

       !***********************************************************************
       width=0.

       ! (3) decide on three or two body decay
       if (.not. present(IDout3)) then  !two body
          ! Loop over all decay channels and search for demanded final state
          do i=1,nDecays
             dID = hadron(ID)%decaysID(i)
             if (dId<=0) cycle !  not 2Body
             if ((Decay2BodyMeson(dId)%ID(1)==IDout1 .and. Decay2BodyMeson(dId)%ID(2)==IDout2) .or. &
                 (Decay2BodyMeson(dId)%ID(2)==IDout1 .and. Decay2BodyMeson(dId)%ID(1)==IDout2)) &
                 width=width &
                 + gammaField(ID,down  )%ratio(i)*gammaField(ID,down  )%gammaTotal*(1.-weight) &
                 + gammaField(ID,down+1)%ratio(i)*gammaField(ID,down+1)%gammaTotal*weight
          end do

       else   !three body decay

          dId_wished = 999

          if ((IDout1.eq.pion).and.(IDout2.eq.pion).and.(IDout3.eq.pion)) then
             ! channel 2 or 3 :  pion pion pion channels
             if ((iQ1+iQ2+iQ3).eq.0) then
                if ((iQ1.ne.0).or.(iQ2.ne.0).or.(iQ3.ne.0)) then
                   dId_wished = -2  ! piPlus, piMinus, piNull channel
                else
                   dId_wished = -3  ! piNull, piNull, piNull channel
                end if
             end if
          else if ((IDout1.eq.eta).or.(IDout2.eq.eta).or.(IDout3.eq.eta)) then
             if ((IDout1+IDout2+IDout3).eq.(2*pion+eta)) then
                ! eta pion pion channel
                if ((iQ1.eq.0).and.(iQ2.eq.0).and.(iQ3.eq.0)) then
                   dId_wished = -1  ! pi0 pi0 eta
                else if (iQ1+iQ2+iQ3.eq.0) then
                   dId_wished = -4  ! pi+ pi- eta
                end if
             end if
          end if

          do i=1,nDecays
             dID = hadron(ID)%decaysID(i) ! 3Body: <0
             if (dId.ne.dId_wished) cycle
             width=width &
                  + gammaField(ID,down  )%ratio(i)*gammaField(ID,down  )%gammaTotal*(1.-weight) &
                  + gammaField(ID,down+1)%ratio(i)*gammaField(ID,down+1)%gammaTotal*weight
          end do

       end if

       partialwidthMeson=width

    else
       write(*,*) 'strange error in partialWidthMeson. STOP!', down, mass,hadron(ID)%minmass,ID
       stop

    end if

  end function partialWidthMeson



  !****************************************************************************
  !****f* mesonWidth/FullWidthMeson
  ! NAME
  ! real function FullWidthMeson(ID,mass)
  ! PURPOSE
  ! This function calculates the full width (energy dependent) of all meson
  ! resonances.
  ! INPUTS
  ! * integer :: ID -- ID of resonance
  ! * real :: mass -- sqrt(p_mu p^mu) = mass of the resonance (offshell)
  ! NOTES
  ! We declare this function to be recursive, since there may be a cycle as
  !   FullWidthMeson -> initWidth -> InitializeSpectralIntegral
  !   -> FullWidthMeson
  !****************************************************************************
  recursive real function FullWidthMeson(ID,mass)
    use IdTable, only: isMeson
    use particleProperties, only: hadron
    use CallStack, only: TRACEBACK

    integer, intent(in) :: ID
    real,    intent(in) :: mass

    integer :: down
    real    :: mass_down,weight

    if (initFlag) call initWidth

    ! (1) Check Input
    if (.not. isMeson(ID)) call TRACEBACK()

    ! (2) Return the full width
    if (mass <= hadron(ID)%minmass) then
       ! If mass is lower than minimalMass, then gamma should be zero
       FullWidthMeson=0.
    else if (mass > hadron(ID)%minmass+deltaMass*(maxIndexMass-1)) then
!       write(*,*) 'Warning in mesonResonanceWidth/widht'
!       Write(*,*) 'Mass of resonance is out of bounds. Mass=', mass
!       write(*,*) 'ID=',ID
       FullWidthMeson=gammaField(ID,maxIndexMass)%gammaTotal
    else if (mass > hadron(ID)%minmass) then
       ! Do linear interpolation between next two grid points "down" and "down+1"
       down = floor((mass-hadron(ID)%minmass)/deltaMass)
       mass_down = hadron(ID)%minmass+float(down)*deltaMass
       weight = (mass-mass_down)/deltaMass ! weight for the interpolation
       FullWidthMeson = gammaField(ID,down  )%gammaTotal * (1.-weight) &
                      + gammaField(ID,down+1)%gammaTotal * weight
    else
      write(*,*) "problem in fullWidthMeson:", ID, mass
      call TRACEBACK()
    end if

    !If(debug) Print *, "In mesonWidth", mass,FullWidthMeson

  end function fullWidthMeson


  !****************************************************************************
  !****s* mesonWidth/initWidth
  ! NAME
  ! subroutine initWidth
  ! PURPOSE
  ! Stores the vacuum width of each meson to the field "gammaField".
  ! Should be called only once.
  !****************************************************************************
  subroutine initWidth
    use IdTable, only: pion, rho, nMes
    use mesonWidthVacuum, only: vacuumWidth
    use output, only: Write_InitStatus
    use particleProperties, only: nDecays, hadron, get_rho_dilep

    integer :: massIndex, ID
    real :: mass, gammaTotal
    real, dimension(nDecays) :: ratio

    if (.not.initFlag) return ! avoid cyclic init calls
    initFlag=.false.

    call Write_InitStatus("widths of the mesons",0)

    ! allocate the field which holds the decay ratio information for each
    ! meson, depending on mass,
    ! first index: meson ID
    ! second index : mass
    allocate(gammaField(pion:pion+nMes-1,0:maxIndexMass))

    ! Initialize the gamma fields for the mesons by calling vacuumWidth
    do ID=pion,pion+nMes-1
       if (debug) write(*,*) "Resonance=", ID
       do massIndex=0,MaxIndexMass
          mass=real(massIndex)*deltaMass+hadron(ID)%minmass
          gammaTotal = vacuumWidth (mass, ID, ratio)
          gammaField(ID,MassIndex)%gammaTotal=gammaTotal
          gammaField(ID,MassIndex)%ratio=ratio
          if (debug) write(100+Id,'(6F14.9)') mass, gammaTotal
          if (.not. (ID==rho .and. get_rho_dilep()) .and. (abs(sum(ratio)-1)>1E-6) .and. (abs(sum(ratio))>1E-6)) then
             write(*,*) 'Problem in mesonWidth/initWidth'
             write(*,*) "Ratios don't add up to 1! ",gammaTotal
             write(*,*) "Ratio :  ", ratio
             write(*,*) 'ResonanceID:', ID
          end if
       end do
    end do

    call InitializeSpectralIntegral( 3.0 )
    call Write_InitStatus("widths of the mesons",1)

  end subroutine initWidth


  !****************************************************************************
  !****s* mesonWidth/cleanUp
  ! subroutine cleanUp
  ! PURPOSE
  ! Deallocate all fields
  !****************************************************************************
  subroutine cleanUp()
    if (initFlag .or. .not. allocated(gammaField)) return
    deallocate(gammaField)
  end subroutine


  !****************************************************************************
  !****f* mesonWidth/decayWidthMeson
  ! NAME
  ! function decayWidthMeson (ID, mass, ch) result(decayWidth)
  ! PURPOSE
  ! This function returns the partial out width (energy dependent) for all decay channels.
  ! INPUTS
  ! * integer :: ID -- ID of resonance
  ! * real    :: mass       -- baremass of the resonance (offshell)
  ! * integer :: ch         -- charge of the resonance
  ! OUTPUT
  ! * real, dimension(nDecays)  :: decayWidth -- widths
  ! NOTES
  ! The D* resonances are treated explicitly since they are the only resonances
  ! which have decay ratios which depend on the charge of the resonance.
  !****************************************************************************
  function decayWidthMeson (ID, mass, ch) result(decayWidth)
    use IDTable, only: dStar, dStarBar, isMeson
    use particleProperties, only: hadron, nDecays
    use decayChannels, only: decay2BodyMeson, decay3BodyMeson
    use CallStack, only: TRACEBACK

    integer,intent(in) :: ID
    real,   intent(in) :: mass
    integer,intent(in) :: ch
    real, dimension(nDecays) :: decayWidth

    integer :: i, dID, down
    real :: thr, mass_down, weight

    if (initFlag) call initWidth

    ! Check Input
    if (.not. isMeson(ID)) call TRACEBACK()

    ! Initialize to zero
    decayWidth = 0.

    if (ID==dStar .or. ID==dStarBar) then
      ! The decays of the D* mesons are treated separately here,
      ! since the D* decays are strongly isospin-dependent.
      if (ch == 0) then
        decayWidth(1) =  2E-3 * 0.381  ! D + photon
        decayWidth(2) =  2E-3 * 0.619  ! D + pion
      else
        decayWidth(1) = 96E-6 * 0.016  ! D + photon
        decayWidth(2) = 96E-6 * 0.984  ! D + pion
      end if
      return
    end if

    do i=1,nDecays
      dID = hadron(ID)%decaysID(i)
      if (dID==0) then
        cycle
      else if (dID>0) then
        thr = decay2BodyMeson(dId)%threshold
      else
        thr = decay3BodyMeson(-dId)%threshold
      end if
      if (mass < thr) then
        decayWidth(i) = 0.
      else if (mass > hadron(ID)%minmass+deltaMass*(maxIndexMass-1)) then
        decayWidth(i) = gammaField(ID,maxIndexMass)%gammaTotal * gammaField(ID,maxIndexMass)%ratio(i)
      else
        ! Do linear interpolation between next two grid points "down" and "down+1"
        down = floor((mass-hadron(ID)%minmass)/deltaMass)
        mass_down = hadron(ID)%minmass+float(down)*deltaMass
        weight = (mass-mass_down)/deltaMass ! weight for the interpolation
        decayWidth(i) = gammaField(ID,down  )%gammaTotal * gammaField(ID,down  )%ratio(i) * (1.-weight) &
                      + gammaField(ID,down+1)%gammaTotal * gammaField(ID,down+1)%ratio(i) * weight
      end if
    end do

  end function decayWidthMeson


  !****************************************************************************
  !****s* mesonWidth/GetMinMaxMass
  ! NAME
  ! subroutine GetMinMaxMass(ID,MinMass,MaxMass,InMedium)
  ! PURPOSE
  ! return values of minimal and maximal mass according the mass tabulation
  ! INPUTS
  ! * integer :: ID -- ID of particle
  ! * logical :: InMedium -- Flag to override minimal mass of vector mesons
  ! OUTPUT
  ! * real :: MinMass -- minimal mass value
  ! * real :: MaxMass -- maximal mass value
  ! NOTES
  ! This returns the minimal mass as stored as default value.
  ! For in-medium vector mesons, the minimal mass is reduced to zero.
  !
  ! The maximal mass is given by the size of the array times the bin width.
  ! All masses are restricted by an upper bound (=3 GeV).
  !****************************************************************************
  subroutine GetMinMaxMass(ID,MinMass,MaxMass,InMedium)
    use particleProperties, only: hadron
    use IdTable, only: rho,omegaMeson,phi

    integer, intent(in) :: ID
    logical, intent(in) :: InMedium
    real, intent(out) :: MinMass,MaxMass

    MinMass = hadron(ID)%minmass
    if (inMedium) then
       select case (ID)
       case (rho,omegaMeson,phi)
          MinMass = 0.01
       end select
    end if

    MaxMass = MinMass + deltaMass*(maxIndexMass-1)

    MaxMass = 3.0
!    if (ID==107) maxMass = 2.0

  end subroutine GetMinMaxMass

  !****************************************************************************
  !****s* mesonWidth/GetMaxQ
  ! NAME
  ! subroutine GetMaxQ(ID,mass0,gamma0,gammaDelta,BinM,BinMaxQ)
  ! PURPOSE
  ! Calculate the maximal values of the Q weight for bins according BinM
  ! INPUTS
  ! * integer :: ID -- ID of resonance
  ! * real :: mass0 -- pole mass
  ! * real :: gamma0 -- width at pole mass
  ! * real :: gammaDelta -- additional width to be added during calculations
  ! * real, dimension(:) :: BinM -- array with boundaries for M binning
  ! OUTPUT
  ! * real, dimension(:) :: BinMaxQ -- the maximal Q values for each bin.
  !
  ! NOTES
  ! * The size of BinMaxQ has to be at least the size of BinM minus 1.
  ! * It first calculates Q at the boundaries, then it iterates over the
  !   tabulated width values in order to take into account, that the Q
  !   value may be larger inbetween the boundaries.
  ! * if the Q value is maximal at the upper bound, we store its value as -Q.
  !****************************************************************************
  subroutine GetMaxQ(ID,mass0,gamma0,gammaDelta,BinM,BinMaxQ)
    use MassAssInfoDefinition, only: MassAssInfoQ

    integer, intent(in) :: ID
    real, intent(in) :: mass0,gamma0,gammaDelta
    real, dimension(:),intent(in)  :: BinM
    real, dimension(:),intent(out) :: BinMaxQ

    integer :: iB, iM,iM1,iM2
    real :: mass,gamma, Q

    ! 1) calculate Q at the boundaries:

    do iB=1,size(BinM)
       gamma = FullWidthMeson(ID,BinM(iB))
       BinMaxQ(iB) = MassAssInfoQ(mass0,gamma0+gammaDelta,BinM(iB),gamma+gammaDelta)
    end do

    do iB=1,size(BinM)-1
       if (BinMaxQ(iB+1).gt.BinMaxQ(iB)) BinMaxQ(iB) = -BinMaxQ(iB+1) ! sign!!
    end do

    ! 2) calculate Q between the boundaries:

    do iB=1,size(BinM)-1
       iM1 = INT((BinM(iB)  -BinM(1))/deltaMass)
       iM2 = INT((BinM(iB+1)-BinM(1))/deltaMass)

       do iM=iM1+1,iM2
          mass  = iM*deltaMass+BinM(1)
          gamma = gammaField(ID,iM)%gammaTotal
          Q = MassAssInfoQ(mass0,gamma0+gammaDelta,mass,gamma+gammaDelta)
          if (Q.gt.abs(BinMaxQ(iB))) BinMaxQ(iB) = Q
       end do
    end do

  end subroutine GetMaxQ


  !****************************************************************************
  !****is* mesonWidth/InitializeSpectralIntegral
  ! NAME
  ! subroutine InitializeSpectralIntegral(massMax)
  ! PURPOSE
  ! Calculates the integral over the spectral functions and stores the
  ! values. Prints a warning message, if the difference to unity is too large.
  !****************************************************************************
  subroutine InitializeSpectralIntegral(massMax)

    use output, only: Write_InitStatus

    real, intent(in) :: massMax ! upper boundary of integration

    integer :: id, nDev
    character(31), parameter :: III='!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*'

    call Write_InitStatus("mesonWidth/InitializeSpectralIntegral",0)

    do id=pion,pion+nmes-1
       arrSpectralIntegral(id) = calcSpectralIntegral(id, massMax)
    end do

    ! the number of exclamation marks indicates the size of derivation
    do id=pion,pion+nmes-1
       nDev = int(abs(arrSpectralIntegral(id)-1.0)/0.05)
       if (nDev>0) then
          write(*,'(" Normalization off for id=",i3," at ",f5.3," GeV: ",f7.5," ",A)') &
               id,massMax,arrSpectralIntegral(id),III(1:nDev)
       end if
    end do

    call Write_InitStatus("mesonWidth/InitializeSpectralIntegral",1)

  end subroutine InitializeSpectralIntegral

  !****************************************************************************
  !****f* mesonWidth/calcSpectralIntegral
  ! NAME
  ! real function calcSpectralIntegral(id,massMax)
  ! PURPOSE
  ! calculate the integral of the spectral function over the mass
  !****************************************************************************
  real function calcSpectralIntegral(id,massMax)

    use constants, only: pi
    use particleProperties, only: hadron

    integer, intent(in) :: id
    real, intent(in) :: massMax ! upper boundary of integration

    integer :: im, nm
    real :: mass0, gamma0
    real :: mmin, mmax, m
    real :: ymin, ymax, dy, y
    real :: sum, gamma, spectral, intfac
    real, parameter :: dy0 = pi/100.

    calcSpectralIntegral = 0.0

    mass0  = hadron(id)%mass
    gamma0 = hadron(id)%width

    if (gamma0 < 1e-3) then
       if (massMax > mass0) calcSpectralIntegral = 1.0
       return ! --> 1.0 or 0.0
    end if

    mmin = hadron(id)%minmass
    mmax = massMax
    if (mmax < mmin) return ! --> 0.0

    ymax = 2*atan(2*(mmax-mass0) / gamma0)
    ymin = 2*atan(2*(mmin-mass0) / gamma0)

    nm = max(int((ymax-ymin)/dy0),1)
    dy  = (ymax-ymin)/float(nm)

    sum = 0.0
    do im=1,nm
       y = ymin+(float(im)-0.5)*dy
       m = 0.5*tan(0.5*y)*gamma0 + mass0
       m = min(max(m,mmin),mmax)
       ! This may create a cyclic/recursive call of fullWidthMeson:
       gamma = fullWidthMeson(id, m)
       spectral = 2./pi * m**2 * gamma / ((m**2-mass0**2)**2+m**2*gamma**2)
       intfac = gamma0 / ((m-mass0)**2+gamma0**2/4.)
       sum = sum + spectral/intfac
    end do
   calcSpectralIntegral = sum * dy

  end function calcSpectralIntegral


  !****************************************************************************
  !****f* mesonWidth/GetSpectralIntegral
  ! NAME
  ! real function GetSpectralIntegral(id)
  ! PURPOSE
  ! return the value of the integral of the spectral function over the mass
  !****************************************************************************
  real function GetSpectralIntegral(id)
    integer, intent(in) :: id
    GetSpectralIntegral = arrSpectralIntegral(id)
  end function GetSpectralIntegral

end module mesonWidth
