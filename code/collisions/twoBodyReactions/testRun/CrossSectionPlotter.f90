!******************************************************************************
!****m* /CrossSectionPlotter
! NAME
! program CrossSectionPlotter
! PURPOSE
! This is a standalone main routine. It is used to plot the cross sections
! used in the code, but is restriced to the total and elastic cross section.
! An interactive web interface to this plotter is available on the GiBUU Homepage:
! https://gibuu.hepforge.org/XSection/
! INPUTS
! The Namelist "Plotter" in the Jobcard and additional (usual) namelists.
! NOTES
! Particle 1 is the (moving) projectile, particle 2 is the (resting) target.
!******************************************************************************
program CrossSectionPlotter
  use inputGeneral, only: readInputGeneral, path_to_input
  use version, only: printVersion
  use particleProperties, only: initParticleProperties, hadron, PartName, validCharge
  use master_2Body, only: XsectionMaster, HiEnergyContrib
  use particleDefinition, only: particle, sqrtS, absMom, kineticEnergy
  use mediumDefinition, only: medium
  use mediumModule, only: mediumAt
  use preEventDefinition, only: preEvent
  use energyCalc, only: energyDetermination
  use RMF, only: getRMF_flag
  use twoBodyTools, only: sqrtS_free
  use CallStack, only: system_command

  implicit none

  type(particle), dimension(1:2)  :: pair               ! incoming pair of particles
  type(preEvent), dimension(1:4)  :: finalState         ! produced final state
  real :: momLRF(0:3),betaToCF(1:3),rHiEnergy,srts,srtS_XS,srtS_vacuum,srtS_corr,mstar(2)
  real :: sigma_lo(0:7),sigma_hi(0:7)
  integer :: i,q1,q2,id1,id2
  logical :: HiFlag,anti1,anti2
  type(medium) :: tMedium
  real :: srts_max = 100.
  character(len=15) :: name(2)
  character(len=100) :: dataFileTotal, dataFileElast

  call PrintVersion

  call readInputGeneral
  call initParticleProperties

  call readinput

  betaToCF = 0.

  pair%Id = (/id1,id2/)
  pair%charge = (/q1,q2/)
  pair%antiparticle = (/anti1,anti2/)
  pair%perturbative = .false.
  pair(1)%event = 1
  pair(2)%event = 2

  pair(1)%mass = hadron(id1)%mass
  pair(2)%mass = hadron(id2)%mass

  pair(1)%position = (/0.,0.,-1./)
  pair(2)%position = (/0.,0.,0./)

  pair(1)%momentum(1:3)=(/0.,0.,0.01/)

  call energyDetermination (pair(1), betaToCF)
  call energyDetermination (pair(2), betaToCF)

  pair(1)%velocity=pair(1)%momentum(1:3)/pair(1)%momentum(0)

  if (.not. validCharge(pair(1))) then
    write(*,'(A,I3,A,I3)') "Error: particle 1 with ID ", id1, " has invalid charge", q1
    stop
  end if
  if (.not. validCharge(pair(2))) then
    write(*,'(A,I3,A,I3)') "Error: particle 2 with ID ", id2, " has invalid charge", q2
    stop
  end if

  dataFileTotal = trim(path_to_input) // "/PDG_data/rpp2012-" // trim(name(1)) // trim(name(2)) // "_total.dat"
  dataFileElast = trim(path_to_input) // "/PDG_data/rpp2012-" // trim(name(1)) // trim(name(2)) // "_elastic.dat"

  tMedium = mediumAt(pair(2)%position)

  write(*,*) '***********************'
  write(*,*) dataFileTotal
  write(*,*) dataFileElast
  write(*,*) '***********************'

  ! produce symbolic links to the data files
  call system_command ("ln -sf "//trim(dataFileTotal)//" dataTotal.dat")
  call system_command ("ln -sf "//trim(dataFileElast)//" dataElast.dat")

  open(23,file='XS.dat')

  do i=1,200000
!     pair(1)%momentum(3)= pair(1)%momentum(3) * 10.**(1./50.)
!     pair(1)%momentum(3)= pair(1)%momentum(3) * 10.**(1./100.)
     pair(1)%momentum(3)= pair(1)%momentum(3) * 10.**(1./500.)

     call energyDetermination(pair(1),betaToCF)
     pair(1)%velocity=pair(1)%momentum(1:3)/pair(1)%momentum(0)

     srtS = sqrtS(pair,"generateFinalState, srtS")
     srtS_vacuum=sqrtS_free(pair)

     if (.not. getRMF_flag()) then
        srtS_XS = srtS_vacuum ! This is the srtS value the XS should be calculated with
     else
        mstar(1) = sqrtS(pair(1),'generateFinalState, mstar(1)')
        mstar(2) = sqrtS(pair(2),'generateFinalState, mstar(2)')
        srtS_corr = srtS - mstar(1) - mstar(2) + pair(1)%mass + pair(2)%mass
        srtS_XS = srtS_corr
     end if

     if (srts > srts_max) exit

     write(*,*) 'plab = ', absmom(pair(1)), ', srts = ', srts, srts_XS

     momLRF = pair(1)%momentum+pair(2)%momentum
     rHiEnergy = HiEnergyContrib(srts,pair%ID,pair%Antiparticle)
     sigma_lo = 0.
     sigma_hi = 0.

     if (rHiEnergy<1.0) then
        call XsectionMaster (srts_XS, pair, tMedium, momLRF, finalState, sigma_lo, &
             HiFlag, .true., ForceHiEnergy=.false.)
     end if

     if (rHiEnergy>0.0) then
        call XsectionMaster (srts_XS, pair, tMedium, momLRF, finalState, sigma_hi, &
             HiFlag, .true., ForceHiEnergy=.true.)
     end if

     write(23,'(10G15.7)') absmom(pair(1)), srts, &
          (1.0-rHiEnergy)*sigma_lo(0:1)+rHiEnergy*sigma_hi(0:1), &
          sigma_lo(0:1), sigma_hi(0:1), &
          rHiEnergy, kineticEnergy(pair(1))


  end do

  close(23)

contains

  subroutine readInput

    use output, only: Write_ReadingInput

    NAMELIST /Plotter/ id1,id2,q1,q2,anti1,anti2,srts_max

    call Write_ReadingInput('Plotter',0)
    rewind(5)
    read(5,nml=Plotter)
    write(*,*) '  Id of first particle      : ',id1
    write(*,*) '  Id of second particle     : ',id2
    write(*,*) '  Charge of first particle  : ', q1
    write(*,*) '  Charge of second particle : ', q2
    write(*,*) '  antiparticle 1            : ', anti1
    write(*,*) '  antiparticle 2            : ', anti2
    write(*,*) '  srt(s)_max                : ', srts_max

    name(1) = PartName(id1,q1,anti1)
    name(2) = PartName(id2,q2,anti2)

    write(*,*)
    write(*,*) '  >> ',trim(name(2)),' --> ',trim(name(1)),' <<'
    write(*,*)

    call Write_ReadingInput('Plotter',1)

    !**************************************************************************
    ! The following generates a html-formatted table of the collision scenario.
    !**************************************************************************

    open(100,file="XS.html")
    write(100,*) '<table border="1" cellspacing="5" cellpadding="10">'
    write(100,*) '<tr>'
    write(100,*) '<th> Scattering particles </th> <th>', name(1),' </th> <th> ', name(2),' </th>'
    write(100,*) '</tr>'
    write(100,*) '<tr>'
    write(100,*) '<td> ID           </td> <td> ',   id1,' </td> <td> ',   id2,' </td>'
    write(100,*) '</tr>'
    write(100,*) '<tr>'
    write(100,*) '<td> Charge       </td> <td> ',    q1,' </td> <td> ',    q2,' </td>'
    write(100,*) '</tr>'
    write(100,*) '<tr>'
    write(100,*) '<td> Antiparticle </td> <td> ', anti1,' </td> <td> ', anti2,' </td>'
    write(100,*) '</tr>'
    write(100,*) '</table>'
    close(100)

  end subroutine readInput


end program CrossSectionPlotter
