program test
  use inputGeneral

  implicit none

  write(*,*) '**************************************************'
  write(*,*) 'Initializing database'
  call readinputGeneral
  call init_Database
  write(*,*) '**************************************************'

  write(*,*) '**************************************************'
  write(*,*) 'Testing Kaon Sigma -> '
  write(*,*) '**************************************************'

  call testmesonY


  contains


  subroutine testmesonY
    use IDTABLE
    use mediumDefinition
    use particleDefinition
    use particleProperties
    use mesonHyperon
    use preEventDefinition
    implicit none
    real                                          :: srts                  ! sqrt(s) in the process
    type(particle),dimension(1:2)   :: teilchenIn        ! colliding particles
    type(medium)                      :: mediumATcollision    ! Medium informations at the position of the collision

    logical                         :: plotFlag          ! Switch on plotting of the  Xsections
    real,dimension(0:3)       :: momentumLRF        ! Total Momentum in LRF
    type(preEvent),dimension(1:3) :: teilchenOut     ! colliding particles
    real                            :: sigmaTot         ! total Xsection
    real                            :: sigmaElast      ! elastic Xsecti
    integer :: i,chargeMeson,chargeNuk, k
    real :: plab,ekin
    real :: piN,piD,kaonN, pipiN,resonanz
    real :: dens
    integer :: numTries=10000


    mediumAtCollision%useMedium=.false.
    mediumAtCollision%densityProton=0.
    mediumAtCollision%densityNeutron=0.

    teilchenIN(1)%Id=pion
    teilchenIN(2)%Id=SigmaResonance
    teilchenIN(1)%charge=-1
    teilchenIN(2)%charge=1
    teilchenIN(1)%mass=meson(pion)%mass
    teilchenIN(2)%mass=baryon(sigmaResonance)%mass


    teilchenIN(1)%momentum=(/0.836625870459866,       0.196406508164487,  0.108609905565974,       0.794057440671839 /)
    teilchenIN(2)%momentum=(/ 1.84423195111350,       0.262005851924195, -0.301559962962424,        1.29413168461444 /)
    
    
    teilchenIN(1)%velocity=teilchenIn(1)%momentum(1:3)/teilchenIn(1)%momentum(0)
    teilchenIN(2)%velocity=teilchenIn(2)%momentum(1:3)/teilchenIn(2)%momentum(0)

    teilchenIn%antiparticle=.false.

    srts=sqrts(teilchenIn(1:2))

    momentumLRF=teilchenIN(1)%momentum+teilchenIN(2)%momentum

    i=0
    Do 

       call mesonY(srts,teilchenIN,mediumATcollision,momentumLRF,teilchenOUT,sigmaTot,sigmaElast,.false.)
       ! Check charge of each particle
       Do k=lbound(teilchenOut,dim=1),ubound(teilchenOut,dim=1)
          If(.not.validCharge_ID(teilchenOut(k)%ID,teilchenOut(k)%charge)) then
             Write(*,*) 'Charge of particle is not valid in finalCheck'
             Write(*,*) 'Charge=',teilchenOut(k)%charge
             Write(*,*) 'Id=',teilchenOut(k)%id
             Write(*,*) 'Antiparticle=',teilchenOut(k)%antiparticle
             write(*,*) 'Severe problem. Code stops!!!!'
             stop
          end if
       End do
       write(*,*) i
       i=i+1
    end do





  end subroutine testmesonY

end program test
