program test
  use inputGeneral

  implicit none

  write(*,*) '**************************************************'
  write(*,*) 'Initializing database'
  call readinputGeneral
  call init_Database
  write(*,*) '**************************************************'

  write(*,*) '**************************************************'
  write(*,*) 'Testing SigmaMeson Nuk -> '
  write(*,*) '**************************************************'

  call testSigmaMesonNuc


  contains


  subroutine testSigmaMesonNuc
    use IDTABLE
    use mediumDefinition
    use particleDefinition
    use particleProperties
    use sigmaNucleon
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
    real :: piN,sigmaMesonN, pipiN,resonanz
    real :: dens
    integer :: numTries=10000

    NAMELIST /initTest/ chargeMeson, chargeNuk, dens
    write (*,*)
    write (*,*) '**Initializing testing parameter'
    write(*,*)  ' Reading input ....'
    rewind(5)
    read(5,nml=initTest)
    write(*,*) ' Set charge to ', chargeMeson,'.'
    write(*,*) ' Set nuk charge to ',  chargeNuk,'.'
    write(*,*) ' Set density to ',  dens,'.'

    mediumAtCollision%useMedium=.true.
    mediumAtCollision%densityProton=dens/2.
    mediumAtCollision%densityNeutron=dens/2.

    teilchenIN(1)%Id=sigmaMeson
    teilchenIN(2)%Id=nucleon
    teilchenIN(1)%charge=chargeMeson
    teilchenIN(2)%charge=chargeNuk
    teilchenIN(1)%mass=meson(sigmaMeson)%mass
    teilchenIN(2)%mass=0.938
    teilchenIN(2)%momentum=(/teilchenIN(2)%mass,0.,0.,0./)
    teilchenIN(2)%velocity=(/0.,0.,0./)

    Open(301,file='SigmaMesonNucleonXsections.dat',status='unknown')
    write(301,*) '# nukCharge=',chargeNuk,'    mesonCharge=', chargeMeson
    Do i=1,250
       plab=i*0.01
!    do i=0,300
!       plab = 0.5*exp( (log(500.)-log(0.5))*i/300. )
       teilchenIN(1)%momentum=(/sqrt(teilchenIN(1)%mass**2+plab**2),plab,0.,0./)
       teilchenIN(1)%velocity=(/teilchenIN(1)%momentum(1)/teilchenIN(1)%momentum(0),0.,0./)

       srts=SQRT((SQRT(meson(sigmaMeson)%mass**2+plab**2)+baryon(nucleon)%mass)**2-plab**2)
       momentumLRF=(/SQRT(meson(sigmaMeson)%mass**2+plab**2)+baryon(nucleon)%mass ,plab,0.,0./)
       ekin=SQRT(meson(sigmaMeson)%mass**2+plab**2)-meson(sigmaMeson)%mass
       call sigmaNuc(srts,teilchenIN,mediumATcollision,momentumLRF,teilchenOUT,sigmaTot,sigmaElast,.true.,2.2,.true.)

       write(301,'(5F8.3)') plab,sigmaTot,sigmaElast

       Do k=lBOund(teilchenOut,dim=1),uBound(teilchenOut,dim=1)
          If ( (teilchenOut(k) %ID.eq.nucleon).and.(teilchenOut(k)%charge.lt.0)) then
             Print *, teilchenOut%charge
             Print *, teilchenOut%ID
             stop 'nukcharge'
          end if
       End do
    end do
    close(301)

    piN=0.
    sigmaMesonN=0.
    pipiN=0.
    Do i=1,numTries
       plab=1
       srts=SQRT((SQRT(meson(sigmaMeson)%mass**2+plab**2)+baryon(nucleon)%mass)**2-plab**2)
       teilchenIN(1)%momentum=(/sqrt(teilchenIN(1)%mass**2+plab**2),plab,0.,0./)
       teilchenIN(1)%velocity=(/teilchenIN(1)%momentum(1)/teilchenIN(1)%momentum(0),0.,0./)

       momentumLRF=(/SQRT(meson(sigmaMeson)%mass**2+plab**2)+baryon(nucleon)%mass ,plab,0.,0./)
       call sigmaNuc(srts,teilchenIN,mediumATcollision,momentumLRF,teilchenOUT,sigmaTot,sigmaElast,.true.,2.3,.false.)
       If ((teilchenOut(1)%ID.eq.pion).and.(teilchenOut(2)%ID.eq.pion).and.(teilchenOut(3)%ID.eq.nucleon)) then
          pipiN= pipiN+sigmatot
       else If ((teilchenOut(1)%ID.eq.sigmaMeson).and.(teilchenOut(2)%ID.eq.nucleon)) then
          sigmaMesonN=sigmaMesonN+sigmatot
       else If ((teilchenOut(1)%ID.eq.P11_1440)) then
          resonanz=resonanz+sigmatot
       else If ((teilchenOut(1)%ID.eq.pion).and.(teilchenOut(2)%ID.eq.nucleon)) then
          piN= piN+sigmatot
       end if
    end do
    write(*,*) 'Xsections at plab= 1. GeV'
    write(*,*) 'Simga for P11(1440) :', resonanz/float(numTries)
    write(*,*) 'Simga for piN. :', piN/float(numTries)
    write(*,*) 'Simga for sigmaMeson nucleon prod. :', sigmaMesonN/float(numTries)
    write(*,*) 'Simga for piPiN prod. :', pipiN/float(numTries)



    teilchenIN(1)%Id=sigmaMeson
    teilchenIN(2)%Id=nucleon
    teilchenIN(1)%charge=chargeMeson
    teilchenIN(2)%charge=-chargeNuk
    teilchenIN(2)%antiparticle=.true.
    teilchenIN(1)%mass=meson(sigmaMeson)%mass
    teilchenIN(2)%mass=0.938
    teilchenIN(2)%momentum=(/teilchenIN(2)%mass,0.,0.,0./)
    teilchenIN(2)%velocity=(/0.,0.,0./)

    Open(301,file='SigmaMesonNucleonXsections_anti.dat',status='unknown')
    write(301,*) '# nukCharge=',chargeNuk,'    mesonCharge=', chargeMeson
    Do i=1,250
       plab=i*0.01
!    do i=0,300
!       plab = 0.5*exp( (log(500.)-log(0.5))*i/300. )
       teilchenIN(1)%momentum=(/sqrt(teilchenIN(1)%mass**2+plab**2),plab,0.,0./)
       teilchenIN(1)%velocity=(/teilchenIN(1)%momentum(1)/teilchenIN(1)%momentum(0),0.,0./)

       srts=SQRT((SQRT(meson(sigmaMeson)%mass**2+plab**2)+baryon(nucleon)%mass)**2-plab**2)
       momentumLRF=(/SQRT(meson(sigmaMeson)%mass**2+plab**2)+baryon(nucleon)%mass ,plab,0.,0./)
       ekin=SQRT(meson(sigmaMeson)%mass**2+plab**2)-meson(sigmaMeson)%mass
       call sigmaNuc(srts,teilchenIN,mediumATcollision,momentumLRF,teilchenOUT,sigmaTot,sigmaElast,.true.,2.3,.false.)

       write(301,'(5F8.3)') plab,sigmaTot,sigmaElast

       Do k=lBOund(teilchenOut,dim=1),uBound(teilchenOut,dim=1)
          If ( (teilchenOut(k) %ID.eq.nucleon).and.(teilchenOut(k)%charge.gt.0)) then
             Print *, teilchenOut%charge
             Print *, teilchenOut%ID
             stop 'nukcharge'
          end if
       End do
    end do
    close(301)



  end subroutine testSigmaMesonNuc

end program test
