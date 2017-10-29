!******************************************************************************
!****m* /JPsiNucleon
! NAME
! module JPsiNucleon
!
! PURPOSE
! Includes the cross sections for J/Psi-nucleon
! elastic scattering and J/Psi dissociation
!
! Public routines:
! * JPsiNuc
!******************************************************************************
module JPsiNucleon

  implicit none
  private

  logical,parameter :: debugFlag=.false.

  ! To decide whether we use the flux correction for the incoming particle velocities
  logical, parameter :: fluxCorrector_flag=.true.

  public :: JPsiNuc

contains

  !****************************************************************************
  !****s* JPsiNucleon/JPsiNuc
  ! NAME
  ! subroutine JPsiNuc(srts,teilchenIN,teilchenOUT,sigmaTot,sigmaElast,plotFlag)
  !
  ! PURPOSE
  ! Evaluates the cross sections for
  ! * J N -> J N,
  ! * J N -> Lambda_c Dbar,
  ! * J N -> Lambda_c D*bar ,
  ! * J N -> N D Dbar
  ! and returns also a "preevent":
  ! * type(preEvent),dimension(1:3), intent(out)               :: teilchenOut    !   outgoing particles
  !
  ! The cross sections are based on calculations of A. Sibirtsev (parameterized by A.L.).
  !
  ! INPUTS
  ! * real, intent(in)                              :: srts                  ! sqrt(s) in the process
  ! * type(particle),dimension(1:2), intent(in)     :: teilchenIn            ! colliding particles
  !
  ! Debugging:
  ! * logical, intent(in),optional                  :: plotFlag              ! Switch on plotting of the  Xsections
  !
  ! OUTPUT
  ! * real, intent(out)                             :: sigmaTot              ! total Xsection
  ! * real, intent(out)                             :: sigmaElast            ! elastic Xsection
  !
  !****************************************************************************
  subroutine JPsiNuc(srts,partIn,partOut,sigmaTot,sigmaElast,plotFlag)
  use idTable, only: nucleon,JPsi,dMeson,dBar,dStarBar,Lambda_cPlus
  use particleDefinition
  use preEventDefinition, only: preEvent
  use twoBodyTools, only: velocity_correction,p_lab,searchInInput,convertToAntiParticles
  use random, only: ranCharge
  use RMF, only: getRMF_flag

  ! INPUT:
  real, intent(in)                              :: srts                  ! sqrt(s) in the process
  type(particle),dimension(1:2), intent(in)     :: partIn            ! colliding particles
  logical, intent(in),optional                  :: plotFlag    ! Switch on plotting of the  Xsections

  ! OUTPUT:
  real, intent(out)                             :: sigmaTot   ! total Xsection
  real, intent(out)                             :: sigmaElast ! elastic Xsection
  type(preEvent),dimension(1:3), intent(out)    :: partOut !   outgoing particles

  !***** Cross sections: **********************************
  real, dimension(1:4) :: sigma
                         ! 1 --- J N -> J N (elastic),
                         ! 2 --- J N -> Lambda_c Dbar,
                         ! 3 --- J N -> Lambda_c D*bar ,
                         ! 4 --- J N -> N D Dbar
  !****************************************************************************

  ! Local variables
  real :: fluxCorrector        ! Correction of the fluxfactor due to different velocities
                               ! in the medium compared to the vacuum
  type(particle) :: jpsi_particle, partNucl
  logical :: antiParticleInput, failFlag

  ! Initialize output
  partOut(:)%ID=0                    ! ID of produced particles
  partOut(:)%charge=0                ! Charge of produced particles
  partOut(:)%antiParticle=.false.    ! Whether produced particles are particles or antiparticles
  partOut(:)%mass=0.                 ! Mass of produced particles

  ! (1) Check  Input
  call searchInInput(partIn,JPsi,nucleon,jpsi_particle,partNucl,failFlag)
  if (failFlag) then
     write(*,*) 'Wrong input in JPsiNuc', partIn%ID
     stop
  end if

  if (jpsi_particle%antiParticle) then   ! This must not happen!
    write(*,*) 'JPsi is antiparticle in "JPsiNuc"!!!',partIn%ID,partIn%antiparticle
    stop
  end if

  if (partNucl%antiParticle) then
    ! Invert an antinucleon in a nucleon:
    partNucl%Charge=-partNucl%Charge
    partNucl%antiparticle=.false.
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
  if (Sum(partOut(:)%Charge).ne.partNucl%charge) then
       write(*,*) 'No charge conservation in JPsiNuc!!! Critical error',&
       &jpsi_particle%Charge,partNucl%Charge,&
       &partOut(:)%Charge,partOut(:)%ID
       stop
  end if

  ! (6) Invert particles in antiParticles if input included antiparticles
  if (antiParticleInput) call convertToAntiParticles(partOut)

contains


  !****************************************************************************
  !****s* JPsiNuc/evaluateXsections
  ! NAME
  ! subroutine evaluateXsections
  !
  ! PURPOSE
  ! Evaluates J N -> J N, J N -> Lambda_c Dbar, J N -> Lambda_c D*bar
  ! and J N -> N D Dbar cross sections
  !
  ! NOTES
  ! There are no resonance contributions to JPsi N scattering.
  !****************************************************************************

  subroutine evaluateXsections

    use parametrizationsBarMes, only: JPsiN

    call JPsiN(srts,sigma)

    ! Flux correction for each channel:
    if (fluxCorrector_flag) sigma=sigma*fluxcorrector

    sigmaElast=sigma(1)
    sigmaTot=sum(sigma(:))

  end subroutine evaluateXsections



  !****************************************************************************
  !****s* JPsiNuc/makeDecision
  ! NAME
  ! subroutine makeDecision
  !
  ! PURPOSE
  ! Chooses randomly one of possible outgoing channels in JPsi-nucleon
  ! collision. Also the charges of outgoing particles are selected.
  !****************************************************************************
  subroutine makeDecision
      use random, only: rn

      real :: cut
      integer,dimension (1:3)  :: izmin,izmax,izout     ! needed for ranCharge
      logical :: ranChargeFlag

      cut=rn()*sigmaTot ! random number for Monte-Carlo decision

      if (sigma(1) >= cut) then                ! Elastic scattering

         partOut(1)%Id=JPsi
         partOut(2)%Id=nucleon

         partOut(1)%charge=0
         partOut(2)%charge=partNucl%charge

      else if (sigma(1)+sigma(2) >= cut) then   ! Dbar Lambda_c

         partOut(1)%Id=Dbar
         partOut(2)%Id=Lambda_cPlus

         partOut(1)%charge=partNucl%charge-1
         partOut(2)%charge=1

      else if (sigma(1)+sigma(2)+sigma(3) >= cut) then   ! D*bar Lambda_c

         partOut(1)%Id=dStarBar
         partOut(2)%Id=Lambda_cPlus

         partOut(1)%charge=partNucl%charge-1
         partOut(2)%charge=1

      else          ! D Dbar N

         partOut(1)%Id=dMeson
         partOut(2)%Id=Dbar
         partOut(3)%Id=nucleon

         izmin=(/0,-1,0/)
         izmax=(/1,0,1/)
         call rancharge(izmin,izmax,partNucl%charge,izout,ranChargeFlag)
         if (.not.ranChargeFlag) write(*,*) 'JPsiNuc/makeDecision error in rancharge :',&
                               & izmin,izmax,partNucl%charge,izout,ranChargeFlag

         partOut(1)%Charge=izout(1)
         partOut(2)%Charge=izout(2)
         partOut(3)%Charge=izout(3)

      end if

  end subroutine makeDecision



    !**************************************************************************
    !****s* JPsiNuc/makeOutput
    ! NAME
    ! subroutine makeOutput
    !
    ! PURPOSE
    ! Writes all cross sections to file as function of srts and plab [GeV]
    ! Filenames:
    ! * 'JPsiN_sigTotElast.dat'        : sigmaTot, sigmaElast
    ! * 'JPsiN_diss.dat'        : outgoing channels with JPsi dissociation
    !**************************************************************************
    subroutine makeOutPut

      logical, save :: initFlag=.true.

      ! The output files
      character(30), dimension(1:2) :: outputFile
      real :: plab


      outputFile(1)='JPsiN_sigTotElast.dat'
      outputFile(2)='JPsiN_diss.dat'

      plab=p_lab(srts,jpsi_particle%mass,partNucl%mass)

      if (initFlag) then
         open(file=outputFile(1),UNIT=101,Status='Replace',Action='Write')
         open(file=outputFile(2),UNIT=102,Status='Replace',Action='Write')
         write(101,*) '#   srts,    plab,   sigmaTot,  sigmaElast'
         write(102,*) '# Col. No.:   Quantity:'
         write(102,*) '# 1           srts (GeV)'
         write(102,*) '# 2           plab (GeV/c)'
         write(102,*) '# 3           cr. sec. in mb of J N -> J N'
         write(102,*) '# 4                             J N -> Lambda_c Dbar'
         write(102,*) '# 5                             J N -> Lambda_c D*bar'
         write(102,*) '# 6                             J N -> N D Dbar'
         initFlag=.false.
      else
         open(file=outputFile(1),UNIT=101,Status='old',Position='Append',Action='Write')
         open(file=outputFile(2),UNIT=102,Status='old',Position='Append',Action='Write')
      end if
      write(101,'(4F9.3)') srts, plab, sigmaTot, sigmaElast
      write(102,'(6F9.3)') srts, plab, sigma(1:4)
      close(101)
      close(102)
    end subroutine makeOutPut

  end subroutine JPsiNuc

end module JPsiNucleon
