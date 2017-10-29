!***************************************************************************
!****m* /CrossSectionPlotter2
! NAME
! program CrossSectionPlotter2
!
! PURPOSE
! This is a standalone main routine. It is used to plot the cross sections
! of specific channels, e.g. p p -> p p omega  or
!                            p p -> p Delta+
!
! INPUTS
! The Namelist "Plotter2" in the Jobcard and additional (usual namelists)
!
! NOTES
! particle 2 is the moving one, i.e. everything is in the rest frame of
! particle 1.
!
! The cross sections plotted here will either come from the low-energy
! resonance model or the high-energy string model or a mixture of both,
! according to the high-energy switching parameters in the namelist
! 'master_2Body', such as 'HiEnergySchwelleBarBar' etc.
!
! In contrast to the other 'CrossSectionPlotter.f90' test case, this one
! not only spits out the total and elastic XS, but also the exclusive and
! inclusive XS for producing a specific particle species (given by 'id_out'),
! as well as the 'semi-exclusive' XS with a number of pions.
!***************************************************************************
program CrossSectionPlotter2

  use version, only: PrintVersion
  use inputGeneral, only: readInputGeneral
  use particleProperties, only: initParticleProperties, hadron, partname
  use master_2Body, only: generateFinalState, HiEnergyContrib
  use particleDefinition
  use idTable, only: nucleon, pion, isBaryon, isMeson
  use preEventDefinition
  use mediumDefinition
  use mediumModule, only: mediumAt
  use energyCalc, only: energyDetermination
  use lorentzTrafo, only: lorentzCalcBeta
  use collisionTerm, only: forceDecays

  implicit none

  type(particle), dimension(1:2) :: pair               ! incoming pair of particles
  type(particle), dimension(1,1:20) :: fState          ! produced final state
  type(particle), dimension(1,1) :: realVectorDummy
  logical :: HiFlag,collisionFlag,anti1,anti2
  integer :: i,j,k,q1,q2,id1,id2,id_out,q_out = 0
  real, dimension(0:3) :: momentum_calc
  real, dimension(1:3) :: betaToCM
  type(medium) :: tMedium
  real :: srts, sigTot, rHiEnergy, sigmaTotal
  real :: sigma_out_excl(0:3)                                               ! index = number of additional pions (semi-excl)
  real :: sigma_out_incl(-1:2), sigma1(-1:2), sigma2(-1:2,-1:2), sigma3(-1:2,-1:2,-1:2)           ! index = charge
  real :: sigma_prongs(0:6)                                                 ! index = number of prongs (charged particles)
  integer :: charge(1:3)

  integer, save :: Nevents  = 1000     ! number of events
  real, save    :: srts_min = 0.0      ! minimum sqrt(s)
  real, save    :: srts_max = 4.0      ! maximum sqrt(s)
  real, save    :: dp       = 0.1      ! momentum interval

  integer :: hiEnergyType, fail, pc1, pc2, pc3
  integer :: N_tot, N_nuc, N_pi, N_prongs, N_out_tot, N_nuc_in
  integer :: N_out(-1:2)                          ! index = charge

  call PrintVersion
  call readInputGeneral
  call initParticleProperties

  call readinput

  ! make ID_OUT particle stable
  hadron(id_out)%stability = 0

  pair%Id=(/id1,id2/)
  pair%charge=(/q1,q2/)
  pair%antiparticle=(/anti1,anti2/)
  pair%perturbative=.false.
  pair(1)%event=1
  pair(2)%event=2

  pair(1)%mass=hadron(id1)%mass
  pair(2)%mass=hadron(id2)%mass

  pair(1)%position = (/ 0., 0., 0.1 /)
  pair(2)%position = (/ 0., 0., 0.0 /)

  pair(1)%momentum(1:3)=0.

  !pair(1)%momentum(0)=sqrt(pair(1)%mass**2+Dot_Product(pair(1)%momentum(1:3),pair(1)%momentum(1:3)))
  call energyDetermination(pair(1),(/0.,0.,0./))
  call energyDetermination(pair(2),(/0.,0.,0./))

  pair(1)%velocity=pair(1)%momentum(1:3)/pair(1)%momentum(0)

  tMedium = mediumAt(pair(1)%position)

  if (isMeson(id_out)) srts_min = max(srts_min, hadron(pair(1)%ID)%mass+hadron(pair(1)%ID)%mass+hadron(id_out)%minmass)

  write(*,*) '***********************'
  write(*,*) 'positions:'
  write(*,*) pair(1)%position
  write(*,*) pair(2)%position
  write(*,*) 'momenta:'
  write(*,*) pair(1)%momentum
  write(*,*) pair(2)%momentum
  write(*,*) 'sqrt(s)=',sqrts(pair(1),pair(2))
  write(*,*) 'Total momentum=',pair(1)%momentum+pair(2)%momentum
  write(*,*) 'sqrts_min = ', srts_min
  write(*,*) '***********************'

  betaToCM = (/0.,0.,0./)

  N_nuc_in = 0
  if (id1==nucleon) N_nuc_in = N_nuc_in + 1
  if (id2==nucleon) N_nuc_in = N_nuc_in + 1

  open(23,file='xs.dat')
  open(24,file='xs_inclexcl.dat')
  open(25,file='xs2.dat')
  open(26,file='xs3.dat')
  open(27,file='xs_prongs.dat')

  Do i=1,1000
     pair(2)%momentum(1:3) = (/ 0., 0., float(i)*dp /)
     pair(2)%momentum(0)   = sqrt(pair(2)%mass**2+Dot_Product(pair(2)%momentum(1:3),pair(2)%momentum(1:3)))

     momentum_calc = pair(1)%momentum(0:3)+pair(2)%momentum(0:3)

     betaToCM = lorentzCalcBeta (momentum_calc)

     call energyDetermination (pair(2), betaToCM)
     pair(2)%velocity=pair(2)%momentum(1:3)/pair(2)%momentum(0)

     srts=sqrts(pair(1),pair(2))

     if (srts < srts_min) cycle
     if (srts > srts_max) exit

     write(*,*) '--- srts = ',srts

     rHiEnergy = HiEnergyContrib(srts,pair%ID,pair%Antiparticle)

     sigmaTotal = 0.
     sigma_out_incl = 0.
     sigma_out_excl = 0.
     sigma1 = 0.
     sigma2 = 0.
     sigma3 = 0.
     sigma_prongs = 0.

evtloop: do j=1,Nevents  ! event loop to produce final states

            fail = 0

            ! (1) generate final state
            do
              call setToDefault(fState)
              fState%ID = -1
              call generateFinalState (pair, fState(1,:),1.,1,0.,collisionFlag,HiFlag,HiEnergyType,sigTot_out=sigTot)
              if (collisionFlag) then
                exit
              else
                fail = fail + 1
                if (fail == 100) then
                  print *,"fail = ",fail,"!"
                  exit evtloop
                end if
              end if
            end do

            if (fail>0) print *,"fail = ",fail

            ! (2) force decays
            call ForceDecays(fState, realVectorDummy, 0.)

            ! (3) analysis
            N_tot = 0
            N_nuc = 0
            N_pi  = 0
            N_out = 0
            N_out_tot = 0
            N_prongs = 0
            charge = 0

            do k=lbound(fState,2),ubound(fState,2)

              if (fState(1,k)%ID < 0) exit
              if (fState(1,k)%ID == 0) cycle

              N_tot = N_tot + 1
              !if (fState(1,k)%ID == ID_out .and. fState(1,k)%charge == q_out) N_out=N_out+1
              if (fState(1,k)%ID == nucleon) N_nuc=N_nuc+1
              if (fState(1,k)%ID == pion) N_pi=N_pi+1
              if (fState(1,k)%charge /= 0) N_prongs=N_prongs+1
              if (fState(1,k)%ID == ID_out) then
                N_out(fState(1,k)%charge) = N_out(fState(1,k)%charge) + 1
                N_out_tot = sum(N_out)
                if (N_out_tot<4) charge(N_out_tot) = fState(1,k)%charge
              end if

            end do

            if (N_prongs < sum(pair%charge) .and. .not. isBaryon(ID_out)) then
              print *, "Error: N_prongs = ", N_prongs
              print *, "IDs: ", fState(1,:)%ID
              print *, "charges: ", fState(1,:)%charge
              print *, "masses: ", fState(1,:)%mass
              stop
            end if

            sigmaTotal = sigmaTotal + sigTot/Nevents
            if (N_out_tot>0) then
              sigma_out_incl(:) = sigma_out_incl(:) + sigTot*N_out(:)/Nevents      ! inclusive cross section
              if (N_out_tot==1 .and. N_pi<4 .and. N_tot==N_out_tot+N_nuc+N_pi) then
                 ! (semi-)exclusive cross sections
                 if (isBaryon(ID_out)) then
                    ! baryon production
                    if (N_nuc==N_nuc_in-1) sigma_out_excl(N_pi) = sigma_out_excl(N_pi) + sigTot/Nevents
                 else
                    ! meson production
                    if (N_nuc==N_nuc_in) sigma_out_excl(N_pi) = sigma_out_excl(N_pi) + sigTot/Nevents
                 end if
              end if
            end if

            if ((N_tot == N_nuc + N_out_tot) .or. &
                (ID_out==nucleon .and. N_tot==N_out_tot)) then
               ! (semi-)exclusive particle production (1-3)
               select case(N_out_tot)
               case(1)
                  sigma1(charge(1)) = sigma1(charge(1)) + sigTot/Nevents
               case(2)
                  pc1 = min(charge(1),charge(2))
                  pc2 = max(charge(1),charge(2))
                  sigma2(pc1,pc2) = sigma2(pc1,pc2) + sigTot/Nevents
               case(3)
                  pc1 = min(charge(1),charge(2))
                  pc3 = max(charge(1),charge(2))
                  pc2 = sum(charge(1:3))-pc1-pc3
                  sigma3(pc1,pc2,pc3) = sigma3(pc1,pc2,pc3) + sigTot/Nevents
               end select
            end if

            if (N_prongs<=6) sigma_prongs(N_prongs) = sigma_prongs(N_prongs) + sigTot/Nevents

         end do evtloop


     write(23,'(9F12.4)') absmom(pair(2)), kineticEnergy(pair(2)), srts, sigmaTotal, sigma_out_excl(0:3), rHiEnergy
     flush(23)

     write(24,'(11F12.4)') absmom(pair(2)), kineticEnergy(pair(2)), srts, sigma_out_incl(-1:2), sigma1(-1:2)
     flush(24)

     write(25,'(13F12.4)') absmom(pair(2)), kineticEnergy(pair(2)), srts, &
                           sigma2(-1,-1:2), sigma2(0,0:2), sigma2(1,1:2), sigma2(2,2)
     flush(25)

     write(26,'(22F12.4)') absmom(pair(2)), kineticEnergy(pair(2)), srts, &
                           sigma3(-1,-1,-1:2), sigma3(-1,0,0:2), sigma3(-1,1,1:2), sigma3(-1,2,2), sigma3(0,0,0:2), &
                           sigma3(0,1,1:2), sigma3(0,2,2), sigma3(1,1,1:2), sigma3(1,2,2)
     flush(26)

     write(27,'(11F12.4)') absmom(pair(2)), kineticEnergy(pair(2)), srts, sigmaTotal, sigma_prongs(0:6)
     flush(27)

  End do

  close(23)
  close(24)
  close(25)
  close(26)
  close(27)

contains

  subroutine readInput
    use output, only: Write_ReadingInput

    NAMELIST /Plotter2/ id1,id2,q1,q2,anti1,anti2,id_out,q_out,Nevents,srts_max,srts_min,dp

    call Write_ReadingInput('Plotter2',0)
    rewind(5)
    read(5,nml=Plotter2)
    write(*,*) '  Id of first particle      : ', id1
    write(*,*) '  Id of second particle     : ', id2
    write(*,*) '  Charge of first particle  : ', q1
    write(*,*) '  Charge of second particle : ', q2
    write(*,*) '  antiparticle 1            : ', anti1
    write(*,*) '  antiparticle 2            : ', anti2
    write(*,*) '  Id of outgoing particle   : ', id_out
    write(*,*) '  Charge of outgoing part.  : ', q_out
    write(*,*) '  Number of events          : ', Nevents
    write(*,*) '  minimum sqrt(s)           : ', srts_min
    write(*,*) '  maximum sqrt(s)           : ', srts_max
    write(*,*) '  momentum interval dp      : ', dp
    write(*,*)
    write(*,*) '  >> ',trim(PartName(id2,q2,anti2)),' --> ',trim(PartName(id1,q1,anti1)),' <<'
    write(*,*)

    call Write_ReadingInput('Plotter2',1)

  end subroutine readInput


end program CrossSectionPlotter2
