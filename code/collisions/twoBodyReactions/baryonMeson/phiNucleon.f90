!******************************************************************************
!****m* /phiNucleon
! NAME
! module phiNucleon
! PURPOSE
! Includes the cross sections for phi-nucleon scattering in the resonance
! regime.
!******************************************************************************
module phiNucleon
  implicit none
  private

  ! Debug-flags:
  logical,parameter :: debugFlag=.false.
  logical,parameter :: debugFlagAnti=.false.

  ! use the flux correction for the incoming particle velocities:
  logical, parameter :: fluxCorrector_flag=.true.

  public :: phiNuc

contains

  !****************************************************************************
  !****s* phiNucleon/phiNuc
  ! NAME
  ! subroutine phiNuc(srts,partIn,mediumAtColl,partOut,sigmaTot,
  ! sigmaElast,useHiEnergy,HiEnergySchwelle,plotFlag,sigmaArr)
  !
  ! PURPOSE
  ! Evaluates phi Nucleon -> anything cross sections and returns also a
  ! "preevent"
  !
  ! INPUTS
  ! * real :: srts --- sqrt(s) in the process
  ! * type(particle),dimension(1:2) :: partIn --- colliding particles
  ! * type(medium) :: mediumAtColl --- Medium at position of the collision
  !
  ! High energy matching:
  ! * logical :: useHiEnergy --- .true. if High-Energy cross sections are
  !   given by paramBarMesHE
  ! * real :: HiEnergySchwelle --- threshold sqrt(s) for paramBarMesHE,
  !   i.e. at which energy the cross sections of paramBarMesHE are used
  !
  ! Debugging:
  ! * logical, optional :: plotFlag --- Switch on plotting of the  Xsections
  !
  ! OUTPUT
  ! * real :: sigmaTot --- total Xsection
  ! * real :: sigmaElast --- elastic Xsection
  !
  ! This routine does a Monte-Carlo-decision according to the partial cross
  ! sections to decide on a final state with maximal 3 final state particles.
  ! These are returned in the vector partOut. The kinematics of these
  ! particles is only fixed in the case of a single produced resonance.
  ! Otherwise the kinematics still need to be established. The result is:
  ! * type(preEvent),dimension(1:3) :: partOut --- colliding particles
  ! * real, dimension(3), optional :: sigmaArr -- partial cross sections
  !
  ! The cross sections are based upon a parametrization by Golubeva.
  ! See routine golub_phi in parametrizationBarMes.
  ! NOTES
  ! Possible final states are :
  ! * 1-particle : baryon Resonances
  ! * 2-particle : pi N, phi N, pi pi N
  !****************************************************************************
  subroutine phiNuc(srts, partIn, mediumAtColl, partOut, sigmaTot, sigmaElast, &
       useHiEnergy, HiEnergySchwelle, plotFlag, sigmaArr)

    use idTable
    use particleDefinition
    use mediumDefinition
    use preEventDefinition, only: preEvent
    use twoBodyTools, only: velocity_correction, convertToAntiParticles, &
         pcm,searchInInput
    use RMF, only: getRMF_flag
    use callstack, only: traceBack

    real, intent(in)                              :: srts
    type(particle),dimension(1:2), intent(in)     :: partIn
    type(medium), intent(in)                      :: mediumAtColl

    logical, intent(in),optional                  :: plotFlag
    logical,intent(in)                            :: useHiEnergy
    real,intent(in)                               :: HiEnergySchwelle

    type(preEvent),dimension(1:3), intent(out) :: partOut
    real, intent(out)                          :: sigmaTot
    real, intent(out)                          :: sigmaElast
    real, dimension(3), intent(out), optional  :: sigmaArr


    ! Cross sections
    real,dimension(-1:1) :: piN  ! -> pi N, index denotes pion charge
    real :: phiN         ! -> phi N
    real :: pipiN        ! -> pi pi N

    ! Local variables
    real :: fluxCorrector
    type(particle) :: partPhi, partNucl
    logical :: antiParticleInput, failFlag

    type(preEvent) :: preEv0


    ! partial cross sections for phi N -> R
    !real , dimension(:),Allocatable :: sigmaRes

    ! Field to store the resonance masses
    !real , dimension(:),Allocatable ::massRes      ! Resonance masses
    !Allocate(sigmaRes(lBound(baryon,dim=1):uBound(baryon,dim=1)))
    !Allocate(massRes(lBound(baryon,dim=1):uBound(baryon,dim=1)))


    antiParticleINPUT=.false. ! .true. if antiparticle in the input

    ! Initialize output
    partOut = preEv0

    ! (1) Check  Input
    call searchInInput(partIn,phi,nucleon,partPhi,partNucl, failFlag)
    if (failFlag)  call traceBack('wrong input')

    if (partPhi%antiParticle)  call traceBack('meson is anti')

    if (partNucl%antiParticle) then
       ! Invert all particles in antiparticles
       partNucl%Charge        =  -partNucl%Charge
       partNucl%antiparticle  = .false.
       partPhi%Charge          =  -partPhi%Charge
       antiParticleInput=.true.
    else
       antiParticleInput=.false.
    end if

    ! Correction of the fluxfactor due to different velocities in the medium compared to the vacuum
    if ( .not.getRMF_flag() ) then
      fluxCorrector=velocity_correction(partIn)
    else
      fluxCorrector=1.
    end if

    ! (2) Evaluate the cross sections
    call evaluateXsections

    ! (2a) only fill array, no event
    if (present(sigmaArr)) then
       sigmaArr = (/ &
            phiN, sum(piN), pipiN /)
       return
    end if

    ! Cutoff to kick the case out, that the cross section is zero
    if (sigmaTot.lt.1E-12) then
       sigmatot=0.
       sigmaElast=0.
       return
    end if

    ! (3) Plot them if wished
    if (present(PlotFlag).or.debugFlag) then
       if (plotFlag.or.debugFlag)  call makeOutput
    end if

    ! (4) Define final state
    call MakeDecision

    ! (5) Check Output
    if (Sum(partOut(:)%Charge).ne.partNucl%charge+partPhi%charge) then
       write(*,*) partPhi%Charge, partNucl%Charge, &
            partOut(:)%Charge,partOut(:)%ID
       call traceback('Charge not conserved')
    end if

    ! (6) Invert particles in antiParticles if input included antiparticles
    if (antiParticleInput) then
       if (debugFlagAnti) write(*,*) partOut
       call convertToAntiParticles(partOut)
       if (debugFlagAnti) write(*,*) partOut
    end if

  contains

    !**************************************************************************
    !****s* phiNuc/evaluateXsections
    ! NAME
    ! subroutine evaluateXsections
    ! PURPOSE
    ! Evaluates phi Nucleon -> anything cross sections
    ! NOTES
    ! There are no resonance contributions to phi N scattering.
    ! The contributions are given by Golubeva (see golub_phi).
    !**************************************************************************
    subroutine evaluateXsections

      use parametrizationsBarMes, only: golub_phi
      use parBarMes_HighEnergy, only: paramBarMesHE
      use idTable, only: nucleon, phi
      use output, only: writeparticle
      use constants, only: mN, mPi
      use particleProperties, only: hadron

      real, dimension(1:4) :: sigmaGolub
      real :: pFinal, pInitial,detailedBalanceFactor

      real :: sigmaTotal_HE,sigmaElast_HE ! High energy matchin


      !########################################################################
      ! Evaluate partial cross sections
      !########################################################################

      sigmaGolub = golub_phi(srts, partPhi%mass)


      !************************************************************************
      ! phi N -> pi N
      !************************************************************************
      ! piN = cross section by Golubeva (pi^- p-> phi N)
      ! by detailed balance - vacuum cross section by resonances

      pFinal=pcm(srts,mPi,mN)
      pInitial=pcm(srts,partPhi%mass,partNucl%mass)
      if (pinitial.lt.1E-12) then
         write(*,*) 'WARNING: pInitial is zero in phiNuc', pinitial
         write(*,*) 'phi meson:'
         call writeparticle(6,0,0,partPhi)
         write(*,*) 'nucleon:'
         call writeparticle(6,0,0,partNucl)
         detailedBalanceFactor= 0.
      else
         detailedBalanceFactor= 1./3.*(pFinal/pInitial)**2
         ! given by detailed balance: factor 1/3 due to (2j+1)-Terms in cross
         ! section and different spins in initial and final state
      end if

      sigmaGolub(1)=3./2.* sigmaGolub(1)
      ! 3./2. since Golub returns cross section pi^- p -> phi N,
      ! with 3./.2 we divide by the isospin-clebsch

      if (partNucl%charge.eq.0) then
         piN = (/ 2./3., 1./3., 0. /) * sigmaGolub(1) * detailedBalanceFactor
      else if (partNucl%charge.eq.1) then
         piN = (/ 0., 1./3., 2./3. /) * sigmaGolub(1) * detailedBalanceFactor
      else
         write(*,*) 'Error in phiNuc', partNucl
      end if

      !************************************************************************
      ! phi N -> phi N
      !************************************************************************

      if (srts>mN+hadron(phi)%minmass) then
        phiN=Max(0., sigmaGolub(3))
      else
        phiN=0.
      end if

      !************************************************************************
      ! phi N -> pi pi N
      !************************************************************************

      ! Evaluate phi N -> pi pi N by taking the total inelastic crossection
      ! by Golub and subtracting then the inelastic cross sections for
      ! "phi N -> N pi"

      if (srts.gt.(2*mPi+mN)) then
         pipiN=Max(0.,sigmaGolub(4)-Sum(piN))
         ! Matching to High energy region
         if (useHiEnergy) then
            call paramBarMesHE(HiEnergySchwelle,phi,nucleon, &
                 partPhi%charge,partNucl%charge,&
                 mediumAtColl,sigmaTotal_HE,sigmaElast_HE)
            pipiN=max(pipiN,sigmaTotal_HE-Sum(piN)-phiN)
         end if
      else
         piPiN=0.
      end if


      !########################################################################
      ! evaluate elastic Xsection
      !########################################################################

      sigmaElast=phiN

      !########################################################################
      ! Do the flux correction for each channel
      !########################################################################

      if (fluxCorrector_flag) then
         ! We do this for each channel since they might show up seperately in the output if makeoutput is called
         sigmaElast=sigmaElast*fluxcorrector
         phiN=phiN*fluxcorrector
         piN=piN*fluxcorrector
         pipiN=pipiN*fluxcorrector
      end if

      !########################################################################
      ! Sum up everything for the total cross section
      !########################################################################

      sigmaTot=phiN + sum( piN ) + pipiN

    end subroutine evaluateXsections

    !**************************************************************************


    subroutine makeDecision
      use random, only: rn, ranCharge

      real :: cut!,cut2
      integer :: totalCharge!,resID

      integer,dimension (1:3) :: izmin,izmax,izout     ! needed for ranCharge
      logical :: ranChargeFlag

      integer :: pionCharge

      cut=rn()*sigmaTot ! random number for Monte-Carlo decision

      totalCharge=partPhi%Charge+partNucl%Charge

      !########################################################################
      ! (1) Two -body final states
      !########################################################################

      ! phi N production
      if (phiN.ge.cut) then
         partOut(1:2)%ID = (/ phi, nucleon /)
         partOut(1:2)%Charge = (/ 0, totalCharge /)
         return
      end if
      cut=cut-phiN

      ! piN production
      do pionCharge=-1,1
         if (piN(pionCharge).ge.cut) then
            partOut(1:2)%ID = (/ pion, nucleon /)
            partOut(1:2)%Charge = (/ pionCharge, totalCharge-pionCharge /)
            return
         end if
         cut=cut-piN(pionCharge)
      end do

      !########################################################################
      ! (2) Three body final state
      !########################################################################
      ! pi pi N production

      if (pipiN.ge.cut) then
         partOut(1:3)%ID = (/pion,pion,nucleon/)
         izmin=(/-1,-1,0/)
         izmax=(/1,1,1/)
         call rancharge(izmin,izmax,totalCharge,izout,ranChargeFlag)
         if (.not.ranChargeFlag) write(*,*) 'Error in rancharge :',izmin,izmax,totalCharge,izout,ranChargeFlag
         partOut(1:3)%Charge = izout(1:3)
         return
      end if

      ! No event was generated:
      write(*,*) 'Error in makedecision of phiNuc'
      write(*,*) srts
      write(*,*) sigmatot
      write(*,*) phiN
      write(*,*) piN
      write(*,*) pipiN
      write(*,*) cut
      stop

    end subroutine makeDecision


    !**************************************************************************
    !****s* phiNuc/makeOutput
    ! NAME
    ! subroutine makeOutput
    ! PURPOSE
    ! Writes all cross sections to file as function of srts and plab [GeV].
    ! Filenames:
    ! * 'phiN_sigTotElast.dat'        : sigmaTot, sigmaElast
    ! * 'phiN_nonStrange_nuk.dat'     : non-strange meson with nucleon in final state
    !**************************************************************************
    subroutine makeOutPut

      logical, save :: initFlag=.true.

      ! The output files
      character(30), dimension(1:2), parameter :: outputFile = (/ 'phiN_sigTotElast.dat   ', 'phiN_nonStrange_nuk.dat' /)
      real :: plab

      plab=SQRT(((srts**2-partPhi%mass**2-partNucl%mass**2)/2./partNucl%mass)**2-partPhi%mass**2)

      if (initFlag) then
         open(file=outputFile(1),UNIT=101,Status='Replace',Action='Write')
         open(file=outputFile(2),UNIT=102,Status='Replace',Action='Write')
         write(101,*) '# srts, plab, sigmaTot, sigmaElast '
         write(102,*) '# srts, plab, piN(-1:1), phiN, pipiN '
         initFlag=.false.
      else
         open(file=outputFile(1),UNIT=101,Status='old',Position='Append',Action='Write')
         open(file=outputFile(2),UNIT=102,Status='old',Position='Append',Action='Write')
      end if
      write(101,'(4F9.3)') srts, plab,sigmaTot, sigmaElast
      write(102,'(7F9.3)') srts, plab,piN, phiN, pipiN
      close(101)
      close(102)
    end subroutine makeOutPut

  end subroutine phiNuc


end module phiNucleon
