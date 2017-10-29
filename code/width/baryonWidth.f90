!******************************************************************************
!****m* /baryonWidth
! NAME
! module baryonWidth
! PURPOSE
! When this module is initialized then all information
! for the VACUUM width is once calculated for all resonances
! and then stored into the field gammaField, which is of
! type gammaFieldType. This is done by initWidth.
! Afterwards this field is used to return full and
! partial width of the resonances in the vacuum by the
! subroutine "partialWidthBaryon, fullWidthBaryon"
! USES
! module baryonWidthVacuum
!******************************************************************************
module baryonWidth

  use IdTable, only: nbar
  use particleProperties, only: nDecays

  implicit none
  private

  ! Parameters to store the partial widths depending on mass
  ! in the range (minimalmass, minimalMass+maxIndex*deltaMass)
  integer, parameter :: maxIndexMass=2000
  real, parameter :: deltaMass=0.004

  ! The type used for storage:
  type tGammaField
     real :: gammatotal           !total decay rate
     real, dimension(nDecays) :: ratio      ! ratio of different decay channels
     real, dimension(nDecays) :: rho_AB_atPole ! Formula 2.76 in Effes PhD thesis evaluated at the polemass of the resonance
  end type tGammaField

  ! Field which holds all the information for concerning  the vacuum width, initialized in "initWidth"
  ! First Index: ID of Resonance
  ! Second Index: Mass Index in the range (0,maxIndexMass)
  type(tGammaField), dimension(:,:), allocatable, save  :: gammaField

  ! Flag to check wether this module is initialized by initWidth
  logical, save :: initFlag=.true.

  ! Switching debugging infos off and on
  logical, parameter :: debug=.false.

  ! Array to hold the integrals over the spectral functions
  real, save :: arrSpectralIntegral(1:nbar) = 0.0


  public :: partialWidthBaryon
  public :: fullWidthBaryon
  public :: decayWidthBaryon
  public :: baryonWidth_gammaN
  public :: GetMaxQ
  public :: cleanUp


  !****************************************************************************
  !****g* baryonWidth/readTable
  ! PURPOSE
  ! There is a tabulation of the widths saved in buuinput which is used to
  ! initialize ('baryonWidthVacuum.dat.bz2').
  ! If you don't want to use this pre-tabulated input,
  ! then you can set "readTable=.false". This is useful for runs on
  ! a cluster where you want to minimize input/output. Also it is necessary
  ! if the decay channels have been modified (cf. DecayChannels.dat).
  ! SOURCE
  !
  logical, save :: readTable = .true.
  !****************************************************************************


  !****************************************************************************
  !****g* baryonWidth/writeTable
  ! PURPOSE
  ! This flag determines whether we write out a new tabulation of the widths
  ! ('baryonWidthVacuum.dat.bz2').
  ! It will only have an effects if readTable == .false. or reading of the
  ! tabulation file fails for some reason.
  ! SOURCE
  !
  logical, save :: writeTable = .false.
  !****************************************************************************


contains


  !****************************************************************************
  !****s* baryonWidth/partialWidthBaryon
  ! NAME
  ! real function partialwidthBaryon (ID, mass, inWidth, mesonID, baryonID, mesonMass, baryonMass)
  ! PURPOSE
  ! This function calculates the mass-dependent partial width for a specific
  ! decay channel of a baryon resonance.
  ! INPUTS
  ! * integer :: ID                           -- id of resonance
  ! * real :: mass                            -- p_mu p^mu = mass of the resonance (offshell)
  ! * logical :: inWidth                      -- .true. => in width(only important for channels with unstable particles), .false. => out width
  ! * integer :: mesonID, baryonID            -- ID's of the decay products which one is interested in
  ! * real, optional :: mesonMass, baryonMass -- Possibility to define the masses of the incoming baryon and meson,
  !                                              needed in the case of the In-Width if one of them is off-shell.
  !                                              Otherwise not relevant.
  !****************************************************************************
  function partialwidthBaryon (ID, mass, inWidth, mesonID, baryonID, mesonMass, baryonMass) result(width)
    use DecayChannels, only: Decay2BodyBaryon
    use IdTable, only: isMeson, isBaryon
    use particleProperties, only: Hadron, nDecays, getAngularMomentum_baryon
    use baryonWidthVacuum, only: interactionRadius, srts_srt_switch
    use distributions, only: BlattWeisskopf
    use CallStack, only: TRACEBACK
    use twoBodyTools, only: pCM

    integer, intent(in) :: ID
    real, intent(in) :: mass
    logical, intent(in) ::  inWidth
    integer, intent(in) :: mesonID, baryonID
    real, intent(in), optional :: mesonMass, baryonMass

    real    :: width, partialVacuumWidth, momentum, gamma_tot, mass_down, weight
    integer :: L, massIndex, i, dId, down, up
    logical :: lDummy

    if (initFlag) call initWidth

    ! (1) Check Input
    if (.not.IsBaryon(ID)) call TRACEBACK()
    if (.not.IsBaryon(baryonID)) call TRACEBACK()
    if (.not.IsMeson(mesonID)) call TRACEBACK()

    ! Evaluate the index of the given mass:
    massIndex=NINT((mass-hadron(ID)%minmass)/deltaMass)

    if (mass<=hadron(ID)%minmass) then
       ! If mass is lower than minimalMas than gamma should be zero
       width=0.

    else if (massIndex>=0) then
       ! Assume constant partial width at very high mass:
       if (massIndex>maxIndexMass) then
          !          write(*,*) 'Warning in baryonResonanceWidth/partialWidth'
          !          Write(*,*) 'Mass of resonance is out of bounds. Mass=', mass
          !          write(*,*) 'ID=',ID
          massIndex=maxIndexMass
       end if

       ! Decide wether you need the width to scatter into the resonance or out of the resonance. And return
       ! the compressed decay ratios and the full width. Compressed means summed over all angular momenta of a
       ! specific decay product.
       width=0.

       if (InWidth) then  !in width

          down = min(maxIndexMass, floor((mass-hadron(ID)%minmass)/deltaMass))
          up   = min(maxIndexMass, down+1)
          mass_down = hadron(ID)%minmass+float(down)*deltaMass
          weight = (mass-mass_down)/deltaMass ! weight for the interpolation

          ! Loop over all decay channels and sum up all width which stem from different angular momentum but same decay products
          do i=1,nDecays
             dID = hadron(ID)%decaysID(i)
             if (dId<=0) cycle !  not 2Body

             ! Decide whether decay channel has demanded final state :
             if ((Decay2BodyBaryon(dId)%ID(1)==mesonID) .and. (Decay2BodyBaryon(dId)%ID(2)==baryonID)) then
                if (Decay2BodyBaryon(dId)%stable(1) .and. Decay2BodyBaryon(dId)%stable(2)) then
                   !STABLE FINAL STATE
                   !hier fehlen die IsoSpion-Clebsche, ansonsten Formel 2.77 aus Effenbergers Dr.-Arbeit:
                   width = width &
                        + gammaField(ID,down)%ratio(i) * gammaField(ID,down)%gammaTotal*(1.-weight) &
                        + gammaField(ID,up  )%ratio(i) * gammaField(ID,up  )%gammaTotal*weight
                else
                   ! UNSTABLE FINAL STATE
                   ! Evaluate momentum of products in resonance rest frame
                   if (present(baryonmass) .and. present(mesonMass)) then
                      momentum = pCM(mass, mesonMass, baryonMass, lDummy)
                   else
                      momentum = pCM(mass, hadron(mesonID)%mass, hadron(baryonID)%mass, lDummy)
                   end if
                   partialVacuumWidth=hadron(ID)%width*hadron(ID)%decays(i)
                   if (partialVacuumWidth<1E-6) cycle
                   ! Get Angular momentum of decay products
                   L = getAngularMomentum_baryon(i,ID)
                   !hier fehlen die IsoSpion-Clebsche, ansonsten Formel 2.77 aus Effenbergers Dr.-Arbeit:
                   if (srts_srt_switch) then
                      width = width + partialVacuumWidth * BlattWeisskopf(momentum*interactionRadius,L)**2 * momentum &
                                                         / (gammaField(ID,massIndex)%rho_AB_atPole(i) * mass**2)
                   else
                      width = width + partialVacuumWidth * BlattWeisskopf(momentum*interactionRadius,L)**2 * momentum &
                                                         / (gammaField(ID,massIndex)%rho_AB_atPole(i) * mass)
                   end if
                end if
             end if
          end do

       else   ! out width
          down = min(maxIndexMass, floor((mass-hadron(ID)%minmass)/deltaMass))
          up   = min(maxIndexMass, down+1)
          mass_down = hadron(ID)%minmass+float(down)*deltaMass
          weight = (mass-mass_down)/deltaMass ! weight for the interpolation
          gamma_tot = gammaField(ID,down)%gammaTotal*(1.-weight) &
                    + gammaField(ID,up  )%gammaTotal*    weight
          ! Loop over all decay channels and sum up all width which stem from different angular momentum but same decay products
          do i=1,nDecays
             dID = hadron(ID)%decaysID(i)
             if (dId<=0) cycle !  not 2Body
             if ((Decay2BodyBaryon(dId)%ID(1)==mesonID) .and. (Decay2BodyBaryon(dId)%ID(2)==baryonID)) then
                width = width + gamma_tot * (gammaField(ID,down)%ratio(i)*(1.-weight) + &
                                             gammaField(ID,up  )%ratio(i)*    weight)
             end if
          end do
       end if

    else
       write(*,*) 'strange error in partialWidthBaryon. STOP!', massIndex, mass, hadron(ID)%minmass, ID, inWidth, mesonID, baryonID
       call traceback()
    end if

  end function partialWidthBaryon


  !****************************************************************************
  !****f* baryonWidth/FullWidthBaryon
  ! NAME
  ! real function FullWidthBaryon(ID,mass)
  ! PURPOSE
  ! This function calculates the mass-dependent total decay width of a
  ! baryon resonance.
  ! INPUTS
  ! * integer :: ID -- id of resonance
  ! * real :: mass -- p_mu p^mu = mass of the resonance (offshell)
  ! NOTES
  ! We declare this function to be recursive, since there may be a cycle as
  !   FullWidthBaryon -> initWidth -> InitializeSpectralIntegral
  !   -> FullWidthBaryon
  !****************************************************************************
  recursive real function FullWidthBaryon (ID, mass)
    use IdTable, only: isBaryon
    use particleProperties, only: hadron
    use CallStack, only: TRACEBACK

    integer, intent(in) :: ID
    real,    intent(in) :: mass

    integer :: down
    real    :: mass_down,weight

    if (initFlag) call initWidth

    ! (1) Check Input
    if (.not. IsBaryon(ID)) call TRACEBACK()

    ! (2) Return the full width
    if (mass < hadron(ID)%minmass) then
       ! if mass is lower than minimalMass, then gamma should be zero
       FullWidthBaryon = 0.
    else if (mass > hadron(ID)%minmass+deltaMass*(maxIndexMass-1)) then
       !       write(*,*) 'Warning in baryonResonanceWidth/widht'
       !       Write(*,*) 'Mass of resonance is out of bounds. Mass=', mass
       !       write(*,*) 'ID=',ID
       FullWidthBaryon = gammaField(ID,maxIndexMass)%gammaTotal
    else
       ! Do linear interpolation between next two grid points "down" and "down+1"
       down = floor((mass-hadron(ID)%minmass)/deltaMass)
       mass_down = hadron(ID)%minmass+float(down)*deltaMass
       weight = (mass-mass_down)/deltaMass ! weight for the interpolation
       FullWidthBaryon = gammaField(ID,down  )%gammaTotal*(1.-weight) &
                       + gammaField(ID,down+1)%gammaTotal*    weight
    end if

  end function fullWidthBaryon


  !****************************************************************************
  !****f* baryonWidth/decayWidthBaryon
  ! NAME
  ! function decayWidthBaryon (ID, mass) result(decayWidth)
  ! PURPOSE
  ! This function returns the mass-dependent partial widths of all decay channels.
  ! INPUTS
  ! * integer :: ID -- id of resonance
  ! * real :: mass -- p_mu p^mu = mass of the resonance (offshell)
  ! OUTPUT
  ! * real, dimension(1:nDecays)) :: decayWidth --  widths of all decay channels
  !****************************************************************************
  function decayWidthBaryon (ID, mass) result(decayWidth)
    use IdTable, only: IsBaryon
    use particleProperties, only: hadron
    use decayChannels, only: decay2BodyBaryon
    use CallStack, only: TRACEBACK

    integer, intent(in) :: ID
    real, intent(in) :: mass
    real, dimension(1:nDecays) :: decayWidth

    integer :: i, dID, down
    real :: mass_down, weight, gamma_tot

    if (initFlag) call initWidth

    ! Check Input
    if (.not. IsBaryon(ID)) call TRACEBACK()

    ! Initialize to zero
    decayWidth = 0.

    if (mass < hadron(ID)%minmass) then
      return
    else if (mass > hadron(ID)%minmass+deltaMass*(maxIndexMass-1)) then
      gamma_tot = gammaField(ID,maxIndexMass)%gammaTotal
    else
      down = floor((mass-hadron(ID)%minmass)/deltaMass)
      mass_down = hadron(ID)%minmass+float(down)*deltaMass
      weight = (mass-mass_down)/deltaMass                          ! weight for the interpolation
      gamma_tot = gammaField(ID,down  )%gammaTotal*(1.-weight) &
                + gammaField(ID,down+1)%gammaTotal*    weight      ! interpolated total width
    end if

    do i=1,nDecays
      dID = hadron(ID)%decaysID(i)
      if (dID<=0) cycle ! 2Body: >0, 3Body: <0
      if (mass < decay2BodyBaryon(dID)%threshold) then
        ! if mass is lower than minimalMass, then gamma should be zero
        decayWidth(i) = 0.
      else if (mass > hadron(ID)%minmass+deltaMass*(maxIndexMass-1)) then
        decayWidth(i) = gamma_tot * gammaField(ID,maxIndexMass)%ratio(i)
      else
        ! Do linear interpolation between next two grid points "down" and "down+1"
        decayWidth(i) = gamma_tot * (gammaField(ID,down  )%ratio(i) * (1.-weight) + &
                                     gammaField(ID,down+1)%ratio(i) *     weight)
      end if
    end do

  end function decayWidthBaryon


  !****************************************************************************
  !****s* baryonWidth/initWidth
  ! NAME
  ! subroutine initWidth
  ! PURPOSE
  ! Stores the vacuum width of each baryon to the field "gammaField"
  ! Should be called only once.
  !****************************************************************************
  subroutine initWidth
    use output, only: Write_InitStatus, intToChar
    use inputGeneral, only: path_to_input
    use baryonWidthVacuum, only: vacuumWidth
    use particleProperties, only: nDecays, hadron, get_rho_dilep
    use IdTable, only: nbar
    use decayChannels, only: get_rhoDelta_is_sigmaDelta

    integer :: massIndex, resonanceID, i
    real :: mass, gammaTotal
    real, dimension(nDecays) :: ratio, rho_ab_atPole
    logical :: fileFailure

    call readinput

    call Write_InitStatus("widths of the baryons",0)

    ! allocate the field which holds the decay ratio information for each baryon, depending on mass,
    ! first index: baryon ID
    ! second index : mass
    allocate(gammaField(1:nbar,0:maxIndexMass))

    if (readTable .and. .not. get_rhoDelta_is_sigmaDelta() .and. .not. get_rho_dilep()) then
       call readFile(fileFailure)
       if (.not. fileFailure) then
          initFlag = .false.
          call InitializeSpectralIntegral( 3.0 )
          call Write_InitStatus("widths of the baryons",1)
          return
       end if
    end if

    write(*,*) "Tabulating the baryonic vacuum widths"

    ! Initialize the gamma fields for the baryons by calling vacuumWidth
    do resonanceID=1,nbar
       if (debug) then
          write(*,*) "Resonance=", resonanceID
       else
          write(*,'(A)',advance='no') '.'
       end if
       do massIndex=0,MaxIndexMass
          mass=real(massIndex)*deltaMass+hadron(resonanceID)%minmass
          ! minimalMass is defined in particleProperties
          gammaTotal = vacuumWidth (mass, resonanceID, ratio, rho_AB_atPole)
          gammaField(resonanceID,MassIndex)%gammaTotal=gammaTotal
          gammaField(resonanceID,MassIndex)%ratio=ratio
          gammaField(resonanceID,MassIndex)%rho_AB_atPole=rho_AB_atPole

          !****************Debugging
          if (debug) write(100+resonanceId,'(6F14.9)') mass, gammaTotal
          if (debug .and. massIndex/=0) then
             ! Just an important check for consistency:
             ! rho_AB_atPole must be independent of the mass index.
             do i=1,size(rho_AB_atPole)
                if (rho_AB_atPole(i).ne.gammaField(resonanceID,MassIndex-1)%rho_AB_atPole(i)) then
                   write(*,*) "Bug: Rho ab at pole is not constant"
                   write(*,*) i, mass,resonanceID,rho_AB_atPole(i), gammaField(resonanceID,MassIndex-1)%rho_AB_atPole(i)
                   stop
                end if
                if ((abs(sum(ratio)-1).gt.1E-6).and.(abs(sum(ratio)).gt.1E-6)) then  !Sum not equal 1 or 0
                   write(*,*) "Bug : sum of ratios not equal 1 in baryonwidth/initWidth", ratio
                   write(*,*) "sum=",sum(ratio)
                   stop
                end if
             end do
          end if
          !****************Debugging end
       end do
    end do

    if (writeTable) call writeFile()
    write(*,*)

    initFlag=.false.
    call InitializeSpectralIntegral( 3.0 )
    call Write_InitStatus("widths of the baryons",1)

  contains


      subroutine checkIOS(i,fileFailure)
        integer, intent(in) :: i
        logical, intent(out) :: fileFailure
        if (i.ne.0) then
           write(*,*) 'Error in opening input file: ',trim(path_to_input)//'/baryonWidthVacuum.dat'
           write(*,*) '  --> We must first tabulate the width'
           fileFailure=.true.
        else
           fileFailure=.false.
        end if
      end subroutine checkIOS


      subroutine writeFile()
        use baryonWidthVacuum, only: getParameters
        use bzip
        real :: a_para,b_para,c_para
        logical :: flag_para
        integer :: resonanceID, massINDEX, delta_flag
        type(bzFile) :: f
        character(len=1200) :: line

        f = bzOpenW(trim(path_to_input)//'/baryonWidthVacuum.dat.bz2')

        ! Write field boundaries to file
        call getParameters (a_para, b_para, c_para, flag_para, delta_flag)
        write(line,'(3E15.4,L8)') a_para,b_para,c_para,flag_para
        call bzWriteLine(f,line)
        write(line,'(3I12,E15.4)') 1,nbar,maxIndexMass,deltaMass
        call bzWriteLine(f,line)
        do resonanceID=1,nbar
          write(line,*) hadron(resonanceID)%minmass
          call bzWriteLine(f,line)
        end do
        do resonanceID=1,nbar
          do massIndex=0,MaxIndexMass
            write(line,'('//intToChar(1+2*nDecays)//'E16.8)') gammaField(resonanceID,MassIndex)%gammaTotal,   &
                                                              gammaField(resonanceID,MassIndex)%ratio,        &
                                                              gammaField(resonanceID,MassIndex)%rho_AB_atPole
            call bzWriteLine(f,line)
          end do
        end do
        call bzCloseW(f)
      end subroutine writeFile


      subroutine readFile(fileFailure)
        use baryonWidthVacuum, only: getParameters
        use bzip
        use output, only: Write_ReadingInput

        logical, intent(out) :: fileFailure
        logical :: flag_para_in, flag_para
        real :: a_para_in, b_para_in, c_para_in, minMass_in, deltaMass_in, a_para, b_para, c_para
        integer :: lb, ub, maxIndexMass_in, ios, ll, delta_flag
        type(bzFile) :: f
        character(len=1200) :: line

        call getParameters (a_para, b_para, c_para, flag_para, delta_flag)

        fileFailure=.false.
        call Write_ReadingInput('baryonWidthVacuum.dat.bz2',0)
        f = bzOpenR(trim(path_to_input)//"/baryonWidthVacuum.dat.bz2")

        ll = 0
        call bzReadLine(f,line,ll)
        read(line(1:ll),*,iostat=ios) a_para_in,b_para_in,c_para_in,flag_para_in
        call checkIOS(ios,fileFailure)
        if (fileFailure) then
           call bzCloseR(f)
           call Write_ReadingInput('baryonWidthVacuum.dat.bz2',1)
           return
        end if

        if ((abs(a_para_in-a_para)>0.00001) .or. (abs(b_para_in-b_para)>0.00001) .or. (abs(c_para_in-c_para)>0.00001) .or. &
            (flag_para_in.neqv.flag_para) .or. delta_flag/=1) then
           write(*,*) a_para_in,a_para
           write(*,*) b_para_in,b_para
           write(*,*) c_para_in,c_para
           write(*,*) flag_para_in,flag_para
           fileFailure=.true.
           call bzCloseR(f)
           call Write_ReadingInput('baryonWidthVacuum.dat.bz2',1)
           write(*,*) 'Error in reading input file (0): ',trim(path_to_input)//'/baryonWidthVacuum.dat.bz2'
           write(*,*) '  --> We must first tabulate the width'
           return
        end if

        ll = 0
        call bzReadLine(f,line,ll)
        read(line(1:ll),*,iostat=ios) lb,ub, maxIndexMass_in,deltaMass_in
        call checkIOS(ios,fileFailure)
        if (fileFailure) then
           call bzCloseR(f)
           call Write_ReadingInput('baryonWidthVacuum.dat.bz2',1)
           return
        end if

        if ((lb/=1) .or. (ub/=nbar) .or. (maxIndexMass_in/=maxIndexMass) .or. (abs(deltaMass-deltaMass_in)>0.00001)) then
           fileFailure=.true.
           call bzCloseR(f)
           call Write_ReadingInput('baryonWidthVacuum.dat.bz2',1)
           write(*,*) lb,ub, maxIndexMass_in,deltaMass_in
           write(*,*) 1,nbar, maxIndexMass,deltaMass
           write(*,*) 'Error in reading input file (2): ',trim(path_to_input)//'/baryonWidthVacuum.dat.bz2'
           write(*,*) '  --> We must first tabulate the width'
           return
        end if

        ll = 0
        do resonanceID=1,nbar
           call bzReadLine(f,line,ll)
           call checkIOS(ios,fileFailure)
           read(line(1:ll),*,iostat=ios) minMass_in
           if (fileFailure) then
              call bzCloseR(f)
              call Write_ReadingInput('baryonWidthVacuum.dat.bz2',1)
              return
           end if

           if ((minMass_in-hadron(resonanceID)%minmass)>0.00001) then
              fileFailure=.true.
              call bzCloseR(f)
              call Write_ReadingInput('baryonWidthVacuum.dat.bz2',1)
              write(*,*) resonanceID, minMass_in, hadron(resonanceID)%minmass
              write(*,*) 'Error in reading input file (3): ',trim(path_to_input)//'/baryonWidthVacuum.dat.bz2'
              write(*,*) '  --> We must first tabulate the width'
              return
           end if
        end do

        ll = 0
        do resonanceID=1,nbar
           do massIndex=0,MaxIndexMass
              call bzReadLine(f,line,ll)
              read(line(1:ll),*,iostat=ios) &
                   &  gammaField(resonanceID,MassIndex)%gammaTotal,gammaField(resonanceID,MassIndex)%ratio,&
                   &  gammaField(resonanceID,MassIndex)%rho_AB_atPole
              call checkIOS(ios,fileFailure)
              if (fileFailure) then
                 call bzCloseR(f)
                 call Write_ReadingInput('baryonWidthVacuum.dat.bz2',1)
                 return
              end if
           end do
        end do

        call bzCloseR(f)
        call Write_ReadingInput('baryonWidthVacuum.dat.bz2',1)

      end subroutine readFile

  end subroutine initWidth


  !****************************************************************************
  !****************************************************************************
  subroutine cleanUp()
    if (initFlag .or. .not.allocated(gammaField)) return
    deallocate(gammaField)
  end subroutine



  !****************************************************************************
  !****f* baryonWidth/BaryonWidth_gammaN
  ! NAME
  ! real function BaryonWidth_gammaN (ID, m_R, m, charge)
  ! PURPOSE
  ! This function calculates the decay width of a baryon resonance going
  ! into a nucleon and a gamma*,
  ! using the matrix elements from hadronTensor_ResProd.
  ! INPUTS
  ! * integer :: ID      -- ID of the baryon resonance
  ! * real :: m_R        -- mass of the baryon resonance
  ! * real :: m          -- invariant mass of the gamma* (dilepton)
  ! * integer :: charge  -- charge of the Delta
  ! OUTPUT
  ! Width in GeV.
  !****************************************************************************
  real function baryonWidth_gammaN (ID, m_R, m, charge)
    use constants, only: pi, electronChargeSQ
    use particleProperties, only: hadron
    use hadronTensor_ResProd, only: hadronTensor_R
    use leptonicID, only: EM
    use constants, only: mN

    integer, intent(in) :: ID
    real, intent(in) :: m_R, m
    integer, intent(in) :: charge

    real, dimension(0:3) :: p_N,p_R
    complex, dimension(0:3,0:3) :: hadronTensor
    complex :: matrixElement_Squared
    real :: p_cm,s

    ! center of mass momentum of nucleon and gamma*
    p_cm = sqrt((m_R**2-(mN+m)**2)*(m_R**2-(mN-m)**2))/(2*m_R)
    ! momenta of nucleon and Delta in CoM frame
    p_N = (/sqrt(mN**2+p_cm**2),0.,0.,p_cm/)
    p_R = (/m_R,0.,0.,0./)
    ! spin of resonance
    s = hadron(ID)%spin
    ! compute matrix element
    if (hadronTensor_R(p_N,p_R,ID,charge,EM,hadronTensor,m_R) ) then
      matrixElement_Squared = 2./(2.*s+1.) * electronChargeSQ &
                              * (-hadronTensor(0,0)+hadronTensor(1,1)+hadronTensor(2,2)+hadronTensor(3,3))
      ! 2s+1 comes from averaging over the resonance spin
      ! factor 2 cancels the average over nucleon spin (which is included in hadronTensor)
    else
      matrixElement_Squared = 0.
    end if

    ! two-body decay width, see Particle Data Book 2006 (38.17)
    BaryonWidth_gammaN = 1./(8.*pi)*REAL(matrixElement_Squared)*p_cm/m_R**2

  end function


  !****************************************************************************
  !****s* baryonWidth/readInput
  ! NAME
  ! subroutine readInput
  ! PURPOSE
  ! Reads input in jobcard out of namelist "BaryonWidth".
  !****************************************************************************
  subroutine readInput
    use output, only: Write_ReadingInput

    !**************************************************************************
    !****n*  baryonWidth/BaryonWidth
    ! NAME
    ! NAMELIST /BaryonWidth/
    ! PURPOSE
    ! Includes the input switches:
    ! * readTable
    ! * writeTable
    !**************************************************************************
    NAMELIST /baryonWidth/ readTable, writeTable

    integer :: ios

    call Write_ReadingInput('baryonWidth',0)
    rewind(5)
    read(5,nml=baryonWidth,IOSTAT=ios)
    call Write_ReadingInput('baryonWidth',0,ios)

    write(*,*) 'Using the tabulation    : ', readTable
    write(*,*) 'Write out new tabulation: ', writeTable

    call Write_ReadingInput('baryonWidth',1)

  end subroutine readInput


  !****************************************************************************
  !****s* baryonWidth/GetMaxQ
  ! NAME
  ! subroutine GetMaxQ (ID, mass0, gamma0, gammaDelta, BinM, BinMaxQ)
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
  ! * if the Q value is maximal at the upper bound, we store its value as
  !   -Q.
  !
  !****************************************************************************
  subroutine GetMaxQ (ID, mass0, gamma0, gammaDelta, BinM, BinMaxQ)
    use particleProperties, only: hadron
    use MassAssInfoDefinition, only: MassAssInfoQ

    integer, intent(in) :: ID
    real, intent(in) :: mass0,gamma0,gammaDelta
    real, dimension(:),intent(in)  :: BinM
    real, dimension(:),intent(out) :: BinMaxQ

    integer :: iB, iM,iM1,iM2
    real :: mass,gamma, Q

    ! 1) calculate Q at the boundaries:

    do iB=1,size(BinM)
       gamma = FullWidthBaryon(ID,BinM(iB))
       BinMaxQ(iB) = MassAssInfoQ(mass0,gamma0+gammaDelta,BinM(iB),gamma+gammaDelta)
    end do

!    write(*,'(A10,10f13.4)') 'BinMaxBWD',BinMaxBWD(1:size(BinBound))

    do iB=1,size(BinM)-1
       if (BinMaxQ(iB+1).gt.BinMaxQ(iB)) BinMaxQ(iB) = -BinMaxQ(iB+1) ! sign!!
    end do

    ! 2) calculate Q between the boundaries:

    do iB=1,size(BinM)-1
       iM1 = INT((BinM(iB)  -hadron(ID)%minmass)/deltaMass)
       iM2 = INT((BinM(iB+1)-hadron(ID)%minmass)/deltaMass)

       do iM=iM1+1,iM2
          mass  = iM*deltaMass+hadron(ID)%minmass
          gamma = gammaField(ID,iM)%gammaTotal
          Q = MassAssInfoQ(mass0,gamma0+gammaDelta,mass,gamma+gammaDelta)
          if (Q.gt.abs(BinMaxQ(iB))) BinMaxQ(iB) = Q
       end do
    end do

  end subroutine GetMaxQ

  !****************************************************************************
  !****is* baryonWidth/InitializeSpectralIntegral
  ! NAME
  ! subroutine InitializeSpectralIntegral(massMax)
  ! PURPOSE
  ! Calculates the integral over the spectral functions and stores the
  ! values.
  ! Prints a warning message, if the difference to unity is too large
  !****************************************************************************
  subroutine InitializeSpectralIntegral(massMax)

    use constants, only: pi
    use output, only: Write_InitStatus
    use particleProperties, only: hadron

    real, intent(in) :: massMax ! upper boundary of integration

    integer :: id
    real :: mass0, gamma0
    real :: mmin, mmax, m
    real :: ymin, ymax, dy, y
    real :: sum, gamma, spectral, intfac
    integer :: im, nm
    integer :: nDev
    real, parameter :: dy0 = pi/100.
    character(31), parameter :: III='!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*'

    call Write_InitStatus("baryonWidth/InitializeSpectralIntegral",0)

    do id=1,nbar
       mass0  = hadron(id)%mass
       gamma0 = hadron(id)%width

       if (gamma0 < 1e-3) then
          if (massMax > mass0) arrSpectralIntegral(id) = 1.0
          cycle
       end if

       mmin = hadron(id)%minmass
       mmax = massMax
       if (mmax < mmin) cycle ! integral is zero

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
          sum = sum + spectral/intfac
       end do
       arrSpectralIntegral(id) = sum * dy

    end do

    ! the number of exclamation marks indicates the size of deviation
    do id=1,nbar
       nDev = int(abs(arrSpectralIntegral(id)-1.0)/0.05)
       if (nDev>0) then
          write(*,'("  Normalization off for id=",i3," at ",f5.3," GeV: ",f7.5," ",A)') &
               id,massMax,arrSpectralIntegral(id),III(1:nDev)
       end if
    end do

    call Write_InitStatus("baryonWidth/InitializeSpectralIntegral",1)


  end subroutine InitializeSpectralIntegral



end module baryonWidth
