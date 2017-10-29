program test
  use inputGeneral, only: readinputGeneral
  use version, only: printversion
  use particleProperties, only: initParticleProperties
  implicit none

  write(*,*) '**************************************************'
  call printversion
  write(*,*) '**************************************************'
  write(*,*) 'Initializing database'
  call readinputGeneral
  call initParticleProperties

!  baryon(delta)%width=0.05

!   call testFullCrossSections
!   call testExclusiveCrossSections
!   call testTwoPion
!   call testRanCharge
!   call testPionBack
!   call testpiTot
!   call testPionNuk_elas
!  call testPionNuk
  call testThreePi

contains

  !****************************************************************************

  subroutine testPionNuk_elas
    use parametrizationsBarMes

    real :: srts, sigma
    integer :: i,charge
    logical :: success

    write(*,*) '**************************************************'
    write(*,*) 'Testing PionNuk elastic background'
    write(*,*) '**************************************************'

    do charge=-1,1
       do i=0,1000
          srts=0.9+i*0.002
          sigma=piN_elastic(charge,srts,success)
          write(900+charge,*) srts,sigma,success
       end do
    end do

  end subroutine testPionNuk_elas

  !****************************************************************************

  subroutine testpiTot

    use parBarMes_HighEnergy, only: paramBarMesHE_pion
    use constants, only: mPi, mN

    integer :: i
    real :: plab,srts
    real, dimension(-1:1)  :: elastic, total

    write(*,*) '**************************************************'
    write(*,*) 'Testing PiTot'
    write(*,*) '**************************************************'

!    Do i=1,6000
!       plab=i*0.001
    do i=0,300
       plab = 0.5*exp( (log(500.)-log(0.5))*i/300. )

       srts=SQRT((SQRT(mPi**2+plab**2)+mN)**2-plab**2)
       call paramBarMesHE_pion(srts,total,elastic)
       write(701,'(12F9.2)') plab,srts,total,Elastic
!       write(701,'(12F9.2)') plab,total(-1),total(1)
    end do

  end subroutine testpiTot

  !****************************************************************************

  subroutine testPionback
    use parametrizationsBarMes
    use constants, only: mPi, mN

    integer :: i
    real :: plab,srts,elastic, chargeEx

    write(*,*) '**************************************************'
    write(*,*) 'Testing PionBack'
    write(*,*) '**************************************************'

    do i=1,600
       plab=i*0.001
       srts=SQRT((SQRT(mPi**2+plab**2)+mN-0.06)**2-plab**2)
       call pionquasibg(1,-1,srts,elastic,chargeEx)
       write(901,'(4F8.3)') plab,srts,Elastic, chargeEx
    end do

  end subroutine testPionback

  !****************************************************************************

  subroutine testPionNuk
    use IdTable
    use mediumDefinition
    use particleDefinition
    use pionNucleon
    use preEventDefinition
    use constants, only: mPi, mN

    real                           :: srts                  ! sqrt(s) in the process
    type(particle), dimension(1:2) :: partIn        ! colliding particles
    type(medium)                   :: mediumAtColl    ! Medium informations at the position of the collision

    real,dimension(0:3)       :: momLRF        ! Total Momentum in LRF
    type(preEvent),dimension(1:4) :: partOut     ! produced particles
    real                            :: sigmaTot         ! total Xsection
    real                            :: sigmaElast      ! elastic Xsecti
    integer :: i,chargePion,chargeNuk, k,j
    real :: plab, pi_Lambda, deltaS, phipiN  !,ekin
    real, save :: dens=0.
    integer, save :: numLoops=10000
    real, dimension (0:1,-1:1) :: kaonPionSigma, kaonPionLambda

    NAMELIST /initTestPionNuk/ chargePion, chargeNuk, dens,numLoops

    write(*,*) '**************************************************'
    write(*,*) 'Testing PionNuk'
    write(*,*) '**************************************************'
    write(*,*)
    write(*,*) '**Initializing testing parameter'
    write(*,*) ' Reading input ....'
    rewind(5)
    read(5,nml=initTestPionNuk)
    write(*,*) ' Set pion charge to ', chargePion,'.'
    write(*,*) ' Set nuk charge to ',  chargeNuk,'.'
    write(*,*) ' Set density to ',  dens,'.'

    mediumAtColl%useMedium=.true.
    mediumAtColl%densityProton=dens/2.
    mediumAtColl%densityNeutron=dens/2.

    partIn(1:2)%Id = (/pion, nucleon/)
    partIn(1:2)%charge = (/chargePion, chargeNuk /)
    partIn(1:2)%mass = (/mPi, mN/)

    partIn(2)%momentum = (/partIn(2)%mass,0.,0.,0./)
    partIn(2)%velocity = 0.

    open(301,file='PionNucleonXsections.dat',status='unknown')
    write(301,*) '# nukCharge=',chargeNuk,'    pionCharge=', chargePion
    do i=1,500
       plab=i*0.01
!    do i=0,300
!       plab = 0.5*exp( (log(500.)-log(0.5))*i/300. )
       partIn(1)%momentum=(/sqrt(partIn(1)%mass**2+plab**2),plab,0.,0./)
       partIn(1)%velocity=(/partIn(1)%momentum(1)/partIn(1)%momentum(0),0.,0./)

       srts=SQRT((SQRT(mPi**2+plab**2)+mN)**2-plab**2)
       momLRF=(/SQRT(mPi**2+plab**2)+mN ,plab,0.,0./)
!        ekin=SQRT(mPi**2+plab**2)-mPi
       call pionNuc(srts,partIn,mediumAtColl,momLRF,partOut,sigmaTot,sigmaElast,.true.)
!       write(301,'(5F8.3)') srts,plab,ekin,sigmaTot,sigmaElast
       write(301,'(5F8.3)') plab,sigmaTot,sigmaElast

       do k=lbound(partOut,dim=1),ubound(partOut,dim=1)
          if ( (partOut(k) %ID.eq.nucleon).and.(partOut(k)%charge.lt.0)) then
             write(*,*) partOut%charge
             write(*,*) partOut%ID
             stop 'nukcharge'
          end if
       end do
       if (plab>=2.5) exit
    end do
    close(301)


    ! Test Event Generator, check that it returns right cross sections
    pi_Lambda=0.
    deltaS=0.
    phipiN=0.

    plab=2.0
    srts=SQRT((SQRT(mPi**2+plab**2)+mN)**2-plab**2)
    momLRF=(/SQRT(mPi**2+plab**2)+mN ,plab,0.,0./)
    partIn(1)%momentum=(/sqrt(partIn(1)%mass**2+plab**2),plab,0.,0./)
    partIn(1)%velocity=(/partIn(1)%momentum(1)/partIn(1)%momentum(0),0.,0./)

    kaonPionSigma=0
    kaonPionLambda=0
    write(*,*) 'Number of Tries:', numLoops

    do i=1,NumLoops
       call pionNuc(srts,partIn,mediumAtColl,momLRF,partOut,sigmaTot,sigmaElast,.false.)
       if ((partOut(1)%ID.eq.delta).and.(partOut(2)%ID.eq.0))  deltaS=deltaS+sigmatot
       if ((partOut(1)%ID.eq.pion).and.(partOut(2)%ID.eq.nucleon).and.(partOut(3)%ID.eq.0)) &
            & pi_Lambda= pi_Lambda+sigmatot
       if ((partOut(1)%ID.eq.phi).and.(partOut(2)%ID.eq.pion).and.(partOut(3)%ID.eq.nucleon)) &
            & phipiN= phipiN+sigmatot
       if ((partOut(1)%ID.eq.kaon).and.(partOut(2)%ID.eq.pion).and.(partOut(3)%ID.eq.sigmaResonance) ) then
          kaonPionSigma(partOut(1)%charge,partOut(2)%charge)= &
               &   kaonPionSigma(partOut(1)%charge,partOut(2)%charge)+sigmatot
       end if
       if ((partOut(1)%ID.eq.kaon).and.(partOut(2)%ID.eq.pion).and.(partOut(3)%ID.eq.Lambda) ) then
          kaonPionLambda(partOut(1)%charge,partOut(2)%charge)= &
               &   kaonPionLambda(partOut(1)%charge,partOut(2)%charge)+sigmatot
       end if
    end do
    write(*,*) 'At plab= 2 GeV:'
    write(*,*) 'Simga for delta prod. :', deltas/float(NumLoops)
    write(*,*) 'Simga for pion nucleon prod. :', pi_Lambda/float(NumLoops)
    write(*,*) 'Simga for phiPiN prod. :', phipiN/float(NumLoops)
    write(*,*)
    write(*,*) 'Simga for kaon pion sigma prod. :'
    do i=0,1
       do j=-1,1
          write(*,'(2(A,I4),A,F15.9)') 'Kaon charge=',i,' Pion charge=',j,' Xsection=', kaonPionSigma(i,j)/float(NumLoops)
       end do
    end do
    write(*,*)
    write(*,*) 'Simga for kaon pion Lambda prod. :'
    do i=0,1
       do j=-1,1
          write(*,'(2(A,I4),A,F15.9)') 'Kaon charge=',i,' Pion charge=',j,' Xsection=', kaonPionLambda(i,j)/float(NumLoops)
       end do
    end do
    write(*,*) 'Incoming particles:'
    write(*,*) 'pion',    ChargePion
    write(*,*) 'nucleon', ChargeNuk
  end subroutine testPionNuk

  !****************************************************************************

  subroutine testRanCharge
    use twoBodyTools

    integer i, iztot
    integer, dimension(1:3)  :: izmin, izmax, izout
    logical :: flag

    write(*,*) '**************************************************'
    write(*,*) 'Testing RanCharge'
    write(*,*) '**************************************************'

    iztot=5
    izmin=(/ -1,1 ,0/)
    izmax=(/ 5,4,3/)

    do i=1,1000
       call rancharge(izmin,izmax,iztot,izout,flag)
       if (.not.flag) then
          write(*,*) 'flag',i
          stop
       end if
       if (sum(izout).ne.iztot) stop'iztot.ne.sum(izout)'
       write(100,*)  izout
    end do

  end subroutine

  !****************************************************************************

  subroutine testTwoPion
    use resonanceCrosssections
    use constants, only: mPi, mN

    real, dimension(-2:2) ::sigma1,sigma
    integer :: i
    real, parameter :: deltaS=0.025
    real :: mom,srts

    write(*,*) '**************************************************'
    write(*,*) 'Testing TwoPion production'
    write(*,*) '**************************************************'

    open(900,File="piMinus_Proton_2pi.dat")
    open(901,File="piPlus_Proton_2pi.dat")
    write(900,*) '# srts[GeV],mom[GeV],sigma_tot[mb],sigma_bg[mb]'
    write(901,*) '# srts[GeV],mom[GeV],sigma_tot[mb],sigma_bg[mb]'
    do i=0,100
       mom=i*deltaS
       srts=SQRT((SQRT(mPi**2+mom**2)+mN)**2-mom**2)
       call sigma_npi_n2pi_resonances(srts,-1,.true.,sigma1)
       call sigma_npi_n2pi_resonances(srts,-1,.false.,sigma)
       write(900,'(12F8.3)') srts, mom, sigma,sigma1
       call sigma_npi_n2pi_resonances(srts,1,.true.,sigma1)
       call sigma_npi_n2pi_resonances(srts,1,.false.,sigma)
       write(901,'(12F8.3)') srts, mom, sigma,sigma1
    end do
    close(900)
    close(901)
  end subroutine testTwoPion

  !****************************************************************************

!!$  subroutine testExclusiveCrossSections
!!$    type(medium) :: mediumAtCollision
!!$    integer :: i
!!$    real, parameter :: deltaS=0.025
!!$    real ,dimension(0:3) :: momLRF
!!$    real ,dimension(1:3) :: position
!!$    real :: mom,srts,mass
!!$    momLRF=0.
!!$
!!$    write(*,*) '**************************************************'
!!$    write(*,*) 'Testing EXCLUSIVE CrossSections'
!!$    write(*,*) '**************************************************'
!!$
!!$    Do i=0,100
!!$       mom=i*deltaS
!!$       srts=SQRT((SQRT(meson(pion)%mass**2+mom**2)+baryon(nucleon)%mass)**2-mom**2)
!!$
!!$
!!$!  real function barMes_R_barMes(srts,idMeson_Initial,idBaryon_Initial,idMeson_Final,idBaryon_Final, &
!!$!       &    chargeMeson_Initial,chargeBaryon_Initial,chargeMeson_Final,chargeBaryon_Final,         &
!!$!       &    background,propagated,MediumAtCollision,momentumLRF,mesonMass_initial,baryonMass_initial,position,perturbative)
!!$
!!$       position=(/0.,0.,0./)
!!$
!!$       write(500,'(8F8.3)') srts, mom, &
!!$&    barMes_R_barMes(srts,pion,nucleon,pion,nucleon,-1,1,-1,1,.false.,.true.,MediumAtCollision,momLRF,position,.false.),  &
!!$&    barMes_R_barMes(srts,pion,nucleon,pion,nucleon,-1,1,0,0,.false.,.true.,MediumAtCollision,momLRF,position,.false.),   &
!!$&    barMes_R_barMes(srts,pion,nucleon,pion,nucleon,1,1,1,1,.false.,.true.,MediumAtCollision,momLRF,position,.false.),     &
!!$&    barMes_R_barMes(srts,pion,nucleon,rho,nucleon,-1,1,0,0,.false.,.true.,MediumAtCollision,momLRF,position,.false.), &
!!$&    barMes_R_barMes(srts,pion,nucleon,eta,nucleon,-1,1,0,0,.false.,.true.,MediumAtCollision,momLRF,position,.false.), &
!!$&    barMes_R_barMes(srts,pion,nucleon,omegaMeson,nucleon,-1,1,0,0,.false.,.false.,MediumAtCollision,momLRF,position,.false.)
!!$
!!$    End do
!!$
!!$    Do i=0,100
!!$       mom=i*deltaS
!!$       srts=SQRT((SQRT(meson(pion)%mass**2+mom**2)+baryon(nucleon)%mass)**2-mom**2)
!!$       write(501,'(8F8.3)') srts, mom, &
!!$&    barMes_R_barMes(srts,pion,nucleon,pion,nucleon,-1,1,-1,1,.false.,.false.,MediumAtCollision,momLRF,position,.false.),  &
!!$&    barMes_R_barMes(srts,pion,nucleon,pion,nucleon,-1,1,0,0,.false.,.false.,MediumAtCollision,momLRF,position,.false.),   &
!!$&    barMes_R_barMes(srts,pion,nucleon,pion,nucleon,1,1,1,1,.false.,.false.,MediumAtCollision,momLRF,position,.false.) ,    &
!!$&    barMes_R_barMes(srts,pion,nucleon,rho,nucleon,-1,1,0,0,.false.,.false.,MediumAtCollision,momLRF,position,.false.), &
!!$&    barMes_R_barMes(srts,pion,nucleon,eta,nucleon,-1,1,0,0,.false.,.false.,MediumAtCollision,momLRF,position,.false.),    &
!!$&    barMes_R_barMes(srts,pion,nucleon,omegaMeson,nucleon,-1,1,0,0,.false.,.false.,MediumAtCollision,momLRF,position,.false.)
!!$
!!$    End do
!!$
!!$    ! Checking conservation of charge, charme, strangenesss
!!$    Print *, barMes_R_barMes(srts,pion,nucleon,eta,nucleon,-1,1,0,1,.false.,.true.,MediumAtCollision,momLRF,position,.false.)
!!$    Print *, barMes_R_barMes(srts,pion,nucleon,kaon,nucleon,-1,1,-1,1,.false.,.true.,MediumAtCollision,momLRF,position,.false.)
!!$    Print *, barMes_R_barMes(srts,pion,nucleon,ds_plus,nucleon,-1,1,1,0,.false.,.true.,MediumAtCollision,momLRF,position,.false.)
!!$
!!$  end subroutine testExclusiveCrossSections

  !****************************************************************************

!!$  subroutine testFullCrossSections
!!$    type(medium) :: mediumAtCollision
!!$    integer :: i
!!$    real, parameter :: deltaS=0.025
!!$    real, dimension(:),Allocatable :: sigma,sigma2,masses
!!$    real :: mom,srts,mass
!!$    real, dimension(0:3) :: momLRF
!!$
!!$    write(*,*) '**************************************************'
!!$    write(*,*) 'Testing FullCrossSections'
!!$    write(*,*) '**************************************************'
!!$
!!$    Allocate(sigma(lbound(baryon,dim=1):ubound(baryon,dim=1)))
!!$    Allocate(sigma2(lbound(baryon,dim=1):ubound(baryon,dim=1)))
!!$    Allocate(masses(lbound(baryon,dim=1):ubound(baryon,dim=1)))
!!$
!!$    ! Testing pi+ proton->X
!!$    Do i=0,100
!!$       mom=i*deltaS
!!$       srts=SQRT((SQRT(meson(pion)%mass**2+mom**2)+baryon(nucleon)%mass)**2-mom**2)
!!$       momLRF=(/sQRT(meson(pion)%mass**2+mom**2)+baryon(nucleon)%mass,mom,0.,0./)
!!$       Print *, 'hello'
!!$       call  barMes2resonance(srts,pion,nucleon,1,1,.false.,mediumAtCollision,momLRF,sigma,masses)
!!$       call  barMes2resonance(srts,pion,nucleon,1,1,.true.,mediumAtCollision,momLRF,sigma2,masses)
!!$       write(100,'(4F8.3)') srts, mom,Sum(sigma),sum(sigma2)
!!$    End do
!!$
!!$    ! Testing pi- proton->X
!!$    Do i=0,100
!!$       mom=i*deltaS
!!$       momLRF=(/sQRT(meson(pion)%mass**2+mom**2)+baryon(nucleon)%mass,mom,0.,0./)
!!$       srts=SQRT((SQRT(meson(pion)%mass**2+mom**2)+baryon(nucleon)%mass)**2-mom**2)
!!$       call  barMes2resonance(srts,pion,nucleon,-1,1,.false.,mediumAtCollision,momLRF,sigma,masses)
!!$       call  barMes2resonance(srts,pion,nucleon,-1,1,.true.,mediumAtCollision,momLRF,sigma2,masses)
!!$       write(101,'(4F8.3)') srts, mom, Sum(sigma), Sum(sigma2)
!!$       write(199,'(12F8.3)') srts, mom, sigma(1:10)
!!$    End do
!!$
!!$    ! Testing rho_0 nucleon->X
!!$    mass=0.3
!!$    Do i=0,100
!!$       mom=i*deltaS
!!$       momLRF=(/sQRT(mass**2+mom**2)+baryon(nucleon)%mass,mom,0.,0./)
!!$       srts=SQRT((SQRT(mass**2+mom**2)+baryon(nucleon)%mass)**2-mom**2)
!!$       call  barMes2resonance(srts,rho,nucleon,0,1,.false.,mediumAtCollision,momLRF,sigma,masses,mass,baryon(nucleon)%mass)
!!$       call  barMes2resonance(srts,rho,nucleon,0,1,.true.,mediumAtCollision,momLRF,sigma2,masses,mass,baryon(nucleon)%mass)
!!$       write(103,'(4F8.3)') srts, mom, Sum(sigma), Sum(sigma2)
!!$    End do
!!$    ! Testing rho_0 nucleon->X
!!$    mass=0.5
!!$    Do i=0,100
!!$       mom=i*deltaS
!!$       momLRF=(/sqrt(mass**2+mom**2)+baryon(nucleon)%mass,mom,0.,0./)
!!$       srts=SQRT((SQRT(mass**2+mom**2)+baryon(nucleon)%mass)**2-mom**2)
!!$       call  barMes2resonance(srts,rho,nucleon,0,1,.false.,mediumAtCollision,momLRF,sigma,masses,mass,baryon(nucleon)%mass)
!!$       call  barMes2resonance(srts,rho,nucleon,0,1,.true.,mediumAtCollision,momLRF,sigma2,masses,mass,baryon(nucleon)%mass)
!!$       write(104,'(4F8.3)') srts, mom, Sum(sigma), Sum(sigma2)
!!$    End do
!!$    ! Testing rho_0 nucleon->X
!!$    mass=0.7
!!$    Do i=0,100
!!$       mom=i*deltaS
!!$       momLRF=(/sQRT(mass**2+mom**2)+baryon(nucleon)%mass,mom,0.,0./)
!!$       srts=SQRT((SQRT(mass**2+mom**2)+baryon(nucleon)%mass)**2-mom**2)
!!$       call  barMes2resonance(srts,rho,nucleon,0,1,.false.,mediumAtCollision,momLRF,sigma,masses,mass,baryon(nucleon)%mass)
!!$       call  barMes2resonance(srts,rho,nucleon,0,1,.true.,mediumAtCollision,momLRF,sigma2,masses,mass,baryon(nucleon)%mass)
!!$       write(105,'(4F8.3)') srts, mom, Sum(sigma), Sum(sigma2)
!!$    End do
!!$    ! Testing rho_0 nucleon->X
!!$    mass=0.9
!!$    Do i=0,100
!!$       mom=i*deltaS
!!$       momLRF=(/sQRT(mass**2+mom**2)+baryon(nucleon)%mass,mom,0.,0./)
!!$       srts=SQRT((SQRT(mass**2+mom**2)+baryon(nucleon)%mass)**2-mom**2)
!!$       call  barMes2resonance(srts,rho,nucleon,0,1,.false.,mediumAtCollision,momLRF,sigma,masses,mass,baryon(nucleon)%mass)
!!$       call  barMes2resonance(srts,rho,nucleon,0,1,.true.,mediumAtCollision,momLRF,sigma2,masses,mass,baryon(nucleon)%mass)
!!$       write(106,'(4F8.3)') srts, mom, Sum(sigma), Sum(sigma2)
!!$    End do
!!$
!!$    ! Testing kaon nucleon->X
!!$    Do i=0,100
!!$       mom=i*deltaS
!!$       srts=SQRT((SQRT(meson(kaonBAR)%mass**2+mom**2)+baryon(nucleon)%mass)**2-mom**2)
!!$       momLRF=(/sQRT(meson(kaonBar)%mass**2+mom**2)+baryon(nucleon)%mass,mom,0.,0./)
!!$       call  barMes2resonance(srts,kaonBAR,nucleon,-1,0,.false.,mediumAtCollision,momLRF,sigma,masses)
!!$       call  barMes2resonance(srts,kaonBAR,nucleon,-1,1,.true.,mediumAtCollision,momLRF,sigma2,masses)
!!$       write(107,'(4F8.3)') srts, mom, Sum(sigma), Sum(sigma2)
!!$    End do
!!$
!!$  end subroutine testFullCrossSections

  !****************************************************************************

  subroutine TestThreePi

    use IdTable
    use mediumDefinition
    use particleDefinition
    use pionNucleon
    use preEventDefinition
    use constants, only: mPi, mN

    real                           :: srts                  ! sqrt(s) in the process
    type(particle), dimension(1:2) :: partIn        ! colliding particles
    type(medium)                   :: mediumAtColl    ! Medium informations at the position of the collision

    real,dimension(0:3)       :: momLRF        ! Total Momentum in LRF
    type(preEvent),dimension(1:4) :: partOut     ! produced particles
    real                            :: sigmaTot         ! total Xsection
    real                            :: sigmaElast      ! elastic Xsecti
    integer :: i,chargePion,chargeNuk, k,j
    real :: plab, pi_Lambda, deltaS, phipiN  !,ekin
    real, save :: dens=0.
    integer, save :: numLoops=10000
    real, dimension(0:1,-1:1) :: kaonPionSigma, kaonPionLambda
    real, dimension(1:20) :: sigmaArr

    NAMELIST /initTestPionNuk/ chargePion, chargeNuk, dens,numLoops

    write(*,*) '**************************************************'
    write(*,*) 'Testing PionNuk'
    write(*,*) '**************************************************'
    write(*,*)
    write(*,*) '**Initializing testing parameter'
    write(*,*) ' Reading input ....'
    rewind(5)
    read(5,nml=initTestPionNuk)
    write(*,*) ' Set pion charge to ', chargePion,'.'
    write(*,*) ' Set nuk charge to ',  chargeNuk,'.'
    write(*,*) ' Set density to ',  dens,'.'

    if (dens>0) then
       mediumAtColl%useMedium=.true.
       mediumAtColl%densityProton=dens/2.
       mediumAtColl%densityNeutron=dens/2.
    end if


    partIn(1:2)%Id = (/pion, nucleon/)
    partIn(1:2)%charge = (/chargePion, chargeNuk /)
    partIn(1:2)%mass = (/mPi, mN/)

    partIn(2)%momentum = (/partIn(2)%mass,0.,0.,0./)
    partIn(2)%velocity = 0.

    open(301,file='PionNucleonXsections.dat',status='unknown')
    write(301,*) '# nukCharge=',chargeNuk,'    pionCharge=', chargePion
    do i=1,500
       plab=i*0.01
!    do i=0,300
!       plab = 0.5*exp( (log(500.)-log(0.5))*i/300. )
       partIn(1)%momentum=(/sqrt(partIn(1)%mass**2+plab**2),plab,0.,0./)
       partIn(1)%velocity=(/partIn(1)%momentum(1)/partIn(1)%momentum(0),0.,0./)

       srts=SQRT((SQRT(mPi**2+plab**2)+mN)**2-plab**2)
       momLRF=(/SQRT(mPi**2+plab**2)+mN ,plab,0.,0./)
!        ekin=SQRT(mPi**2+plab**2)-mPi
       call pionNuc(srts,partIn,mediumAtColl,momLRF,partOut,sigmaTot,sigmaElast,.false.,sigmaArr)
!       write(301,'(5F8.3)') srts,plab,ekin,sigmaTot,sigmaElast
!       write(301,'(5F8.3)') plab,sigmaTot,sigmaElast

       write(301,'(50F8.3)') srts,plab,sigmaTot,sigmaElast,sigmaArr

       if (plab>=2.5) exit
    end do
    close(301)




  end subroutine TestThreePi

end program test
