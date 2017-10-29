program testQE
  use inputGeneral
  use particleDefinition
  use nucleusDefinition
  use nucleus, only : getTarget
  use baryonPotentialModule, only: HandPotentialToDensityStatic
  use particleProperties, only: initParticleProperties
  use output, only: WriteParticle
  use minkowski, only: abs3
  use eN_eventDefinition
  use eN_event, only: eNev_SetProcess, eNev_init_BnQ, eNev_init_Target  ! ,write_electronNucleon_event
  use quasiElastic_electron, only: dSigmadcosTheta_l_dE_l_BW_eN
  use lepton2p2h, only: lepton2p2h_DoQE, lepton2p2h_DoDelta
  use pauliBlockingModule, only: checkPauli

  implicit none

  type(particle),Allocatable :: realParticles(:,:)
  type(particle) :: Part
  type(tNucleus),pointer :: targetNuc => NULL()
  integer :: lengthReal=0
  integer :: iEns, iPart, iCost, iNu
!  integer, parameter :: nEns=100
  real :: Cost, pTot, sigma, sigmaPS
  real,dimension(0:1) :: SigmaSum,SigmaSumPS,SigmaSum2p2hQE,SigmaSum2p2hDelta
  type(electronNucleon_event), save :: eNev_InitData, eNev_RunData
  logical :: flagOK
  real :: Ebeam,nu,Q2
  real, dimension(0:3) :: pf
  real :: nuc_bareMass

  integer, parameter :: iNuMin = 0
  integer :: iNuMax=100
  real, parameter :: NuDelta = 0.002

  type(particle),dimension(2) :: OutPart

!  Ebeam = 0.7
!  Q2 = 0.2

!  Ebeam = 0.262
!  Q2 = 0.107

  Ebeam = 0.354
  Q2 = 0.107
  iNuMax = nint((Ebeam-Q2/(4*Ebeam))/NuDelta)


  call readInputGeneral
  call initParticleProperties

  targetNuc => getTarget()
  lengthReal = targetNuc%mass

  if (targetNuc%ReAdjustForConstBinding) then
     write(*,*) 'we now have to readjust the density'
     call HandPotentialToDensityStatic(targetNuc)
  end if

  Allocate(realparticles(1:numEnsembles,1:lengthReal))

  realparticles(:,:)%ID = -1 ! empty vector
  call SetUpTarget(.true.)
  ! setting all fermi momenta in z-direction
  do iEns=1,numEnsembles
     do iPart=1,lengthReal
        realparticles(iEns,iPart)%momentum(3)=abs3(realparticles(iEns,iPart)%momentum)
        realparticles(iEns,iPart)%momentum(1:2)=0
     end do
  end do
  !     call WriteParticle(101,1, realparticles(1,:))


  call eNev_SetProcess(eNev_InitData, 1,1)  ! set to EM and electron


  
!  do iNu=0,30
!!     nu = iNu*0.01
!     nu = iNu*0.005
!  do iNu=0,75
!     nu = iNu*0.01
!  nu = iNu*0.002

  do iNu=iNuMin,iNuMax
     nu = iNu*NuDelta

     SigmaSum = 0.0
     SigmaSumPS = 0.0
     SigmaSum2p2hQE = 0.0
     SigmaSum2p2hDelta = 0.0

     call eNev_init_BnQ(eNev_InitData, Ebeam,nu,Q2, flagOK)
!     call write_electronNucleon_event(eNev_InitData, .false., .true.)

     do iEns = 1,numEnsembles
        do iPart=1,lengthReal
           pTot = realparticles(iEns,iPart)%momentum(3)

           !        write(*,*) '################'

           do iCost = 100,-100,-1
              cost = iCost*0.01
              eNev_RunData = eNev_InitData
              Part = realparticles(iEns,iPart)
              Part%momentum(3)=pTot*cost
              Part%momentum(2)=pTot*sqrt(1-cost**2)
              call eNev_init_Target(eNev_RunData,Part,flagOK)


              !***** QE ******

              sigma = dSigmadcosTheta_l_dE_l_BW_eN(eNev_RunData,pf,nuc_bareMass,sigmaPS)

              SigmaSum(Part%charge) = SigmaSum(Part%charge) + sigma
              SigmaSumPS(Part%charge) = SigmaSumPS(Part%charge) + sigmaPS

              !***** 2p2h -- QE ******

              call lepton2p2h_DoQE(eNev_RunData,outPart,sigma)

              flagOK = checkPauli(outPart, realparticles)
              if (.not.flagOK) then
!                 write(*,*) 'Redo event: pauli blocked! '
                 sigma = 0
              end if

              SigmaSum2p2hQE(Part%charge) = SigmaSum2p2hQE(Part%charge) + sigma

              !***** 2p2h -- Delta ******

              call lepton2p2h_DoDelta(eNev_RunData,outPart,sigma)

              flagOK = checkPauli(outPart, realparticles)
              if (.not.flagOK) then
!                 write(*,*) 'Redo event: pauli blocked! '
                 sigma = 0
              end if

              SigmaSum2p2hDelta(Part%charge) = SigmaSum2p2hDelta(Part%charge) + sigma

              !           write(*,*) iCost,sigma

              


           end do

        end do

     end do

     write(*,'(f13.5,1P,10e13.5)') nu,&
          & SigmaSum(1)/numEnsembles,SigmaSumPS(1)/numEnsembles,&
          & SigmaSum(0)/numEnsembles,SigmaSumPS(0)/numEnsembles,&
          & SigmaSum2p2hQE(1)/numEnsembles,SigmaSum2p2hQE(0)/numEnsembles,&
          & SigmaSum2p2hDelta(1)/numEnsembles,SigmaSum2p2hDelta(0)/numEnsembles
     write(102,'(f13.5,1P,10e13.5)') nu,&
          & SigmaSum(1)/numEnsembles,SigmaSumPS(1)/numEnsembles,&
          & SigmaSum(0)/numEnsembles,SigmaSumPS(0)/numEnsembles
     write(103,'(f13.5,1P,10e13.5)') nu,&
          & SigmaSum(1)/numEnsembles,SigmaSumPS(1)/numEnsembles,&
          & SigmaSum(0)/numEnsembles,SigmaSumPS(0)/numEnsembles,&
          & SigmaSum2p2hQE(1)/numEnsembles,SigmaSum2p2hQE(0)/numEnsembles,&
          & SigmaSum2p2hDelta(1)/numEnsembles,SigmaSum2p2hDelta(0)/numEnsembles

  end do

contains

  subroutine SetUpTarget(DoUpdateEnergies)
    use initNucleus_in_PS,only : initNucPhaseSpace
    use densityModule, only: updateDensity
    use coulomb, only: updateCoulomb
    use yukawa, only: updateYukawa
    use propagation, only: updateVelocity
    use energyCalc, only: updateEnergies
    
    logical, intent(in), OPTIONAL :: DoUpdateEnergies
    
    logical :: DoUE
    
    DoUE = .true.
    if (present(DoUpdateEnergies)) DoUE = DoUpdateEnergies
    
    call initNucPhaseSpace(realparticles,targetNuc)
    call updateDensity(realParticles)
    call updateCoulomb
    call updateYukawa(.true.)
    call updateVelocity(realParticles)
    if (DoUE) call updateEnergies(realParticles)
    
  end subroutine SetUpTarget


end program testQE

