program test
  use inputGeneral, only: readinputGeneral
  use particleProperties, only: initParticleProperties

  implicit none

  write(*,*) '**************************************************'
  write(*,*) 'Initializing database'
  call readinputGeneral
  call initParticleProperties
  write(*,*) '**************************************************'

  write(*,*) '**************************************************'
  write(*,*) 'Testing Rho Nuk -> '
  write(*,*) '**************************************************'

  call testRhoNuc


contains


  subroutine testRhoNuc
    use IDTABLE, only: nucleon, pion, rho, P13_1900
    use mediumDefinition
    use particleDefinition
    use particleProperties, only: hadron
    use rhoNucleon, only: rhoNuc
    use preEventDefinition
    use constants, only: mN, rhoNull

    real                          :: srts               ! sqrt(s) in the process
    type(particle),dimension(1:2) :: teilchenIn         ! colliding particles
    type(medium)                  :: mediumATcollision  ! Medium informations at the position of the collision

    real,dimension(0:3)           :: momentumLRF      ! Total Momentum in LRF
    type(preEvent),dimension(1:3) :: teilchenOut      ! colliding particles
    real                          :: sigmaTot         ! total Xsection
    real                          :: sigmaElast       ! elastic Xsection
    integer :: i,j
    integer :: chargeMeson = 0, chargeNuk = 1
    real :: dens = rhoNull
    real :: plab,ekin,mass,piN,rhoN, pipiN,r_p13
    integer, parameter :: numTries=10000

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

    teilchenIN(1)%Id=rho
    teilchenIN(1)%charge=chargeMeson
    teilchenIN(1)%mass=hadron(rho)%mass
    mass = hadron(rho)%mass

    teilchenIN(2)%Id=nucleon
    teilchenIN(2)%charge=chargeNuk
    teilchenIN(2)%mass=mN
    teilchenIN(2)%momentum=(/teilchenIN(2)%mass,0.,0.,0./)
    teilchenIN(2)%velocity=(/0.,0.,0./)

    ! (1) Test momentum dependence (with all channels)
    Do i=1,300
      plab=i*0.01
      teilchenIN(1)%momentum=(/sqrt(mass**2+plab**2),plab,0.,0./)
      teilchenIN(1)%velocity=(/teilchenIN(1)%momentum(1)/teilchenIN(1)%momentum(0),0.,0./)

      srts=SQRT((SQRT(mass**2+plab**2)+mN)**2-plab**2)
      momentumLRF=(/SQRT(mass**2+plab**2)+mN,plab,0.,0./)
      call rhoNuc(srts,teilchenIN,mediumATcollision,momentumLRF,teilchenOUT,sigmaTot,sigmaElast,.true.,2.2,.true.)
    end do


    ! (2) Test mass & momentum dependence (tot. & elast. only)
    Open(301,file='RhoNucleonXsections.dat',status='unknown')
    write(301,*) '# nukCharge=',chargeNuk,'    mesonCharge=', chargeMeson
    Do i=1,100
      plab=i*0.02
      Do j=1,100
        mass=j*0.02
        teilchenIN(1)%mass=mass
        teilchenIN(1)%momentum=(/sqrt(mass**2+plab**2),plab,0.,0./)
        teilchenIN(1)%velocity=(/teilchenIN(1)%momentum(1)/teilchenIN(1)%momentum(0),0.,0./)

        srts=SQRT((SQRT(mass**2+plab**2)+mN)**2-plab**2)
        momentumLRF=(/SQRT(mass**2+plab**2)+mN,plab,0.,0./)
        ekin=SQRT(mass**2+plab**2)-mass
        call rhoNuc(srts,teilchenIN,mediumATcollision,momentumLRF,teilchenOUT,sigmaTot,sigmaElast,.true.,2.2,.false.)

        write(301,'(5F12.3)') plab,mass,srts,sigmaTot,sigmaElast
      end do
      write(301,*)
      print *,i
    end do
    close(301)

    ! (3) check special channels at fixed plab
    piN=0.
    rhoN=0.
    pipiN=0.
    Do i=1,numTries
       plab=0.2
       srts=SQRT((SQRT(hadron(rho)%mass**2+plab**2)+mN)**2-plab**2)
       teilchenIN(1)%momentum=(/sqrt(teilchenIN(1)%mass**2+plab**2),plab,0.,0./)
       teilchenIN(1)%velocity=(/teilchenIN(1)%momentum(1)/teilchenIN(1)%momentum(0),0.,0./)

       momentumLRF=(/SQRT(hadron(rho)%mass**2+plab**2)+mN ,plab,0.,0./)
       call rhoNuc(srts,teilchenIN,mediumATcollision,momentumLRF,teilchenOUT,sigmaTot,sigmaElast,.true.,2.3,.false.)
       If ((teilchenOut(1)%ID==pion).and.(teilchenOut(2)%ID==pion).and.(teilchenOut(3)%ID==nucleon)) then
          pipiN= pipiN+sigmatot
       else If ((teilchenOut(1)%ID==rho).and.(teilchenOut(2)%ID==nucleon)) then
          rhoN=rhoN+sigmatot
       else If ((teilchenOut(1)%ID==P13_1900)) then
          R_p13=R_p13+sigmatot
       else If ((teilchenOut(1)%ID==pion).and.(teilchenOut(2)%ID==nucleon)) then
          piN= piN+sigmatot
       end if
    end do
    write(*,*) 'Xsections at plab= 0.2 GeV'
    write(*,*) 'Simga for P13_1900. :', r_p13/float(numTries)
    write(*,*) 'Simga for piN. :', piN/float(numTries)
    write(*,*) 'Simga for rho nucleon prod. :', rhoN/float(numTries)
    write(*,*) 'Simga for piPiN prod. :', pipiN/float(numTries)


    ! (4) check rho-Nbar
    teilchenIN(1)%Id=rho
    teilchenIN(1)%charge=chargeMeson
    teilchenIN(1)%mass=hadron(rho)%mass

    teilchenIN(2)%Id=nucleon
    teilchenIN(2)%charge=-chargeNuk
    teilchenIN(2)%antiparticle=.true.
    teilchenIN(2)%mass=mN
    teilchenIN(2)%momentum=(/teilchenIN(2)%mass,0.,0.,0./)
    teilchenIN(2)%velocity=(/0.,0.,0./)

    Open(301,file='RhoNucleonXsections_anti.dat',status='unknown')
    write(301,*) '# nukCharge=',chargeNuk,'    mesonCharge=', chargeMeson
    Do i=1,250
       plab=i*0.01
       teilchenIN(1)%momentum=(/sqrt(teilchenIN(1)%mass**2+plab**2),plab,0.,0./)
       teilchenIN(1)%velocity=(/teilchenIN(1)%momentum(1)/teilchenIN(1)%momentum(0),0.,0./)

       srts=SQRT((SQRT(hadron(rho)%mass**2+plab**2)+mN)**2-plab**2)
       momentumLRF=(/SQRT(hadron(rho)%mass**2+plab**2)+mN ,plab,0.,0./)
       ekin=SQRT(hadron(rho)%mass**2+plab**2)-hadron(rho)%mass
       call rhoNuc(srts,teilchenIN,mediumATcollision,momentumLRF,teilchenOUT,sigmaTot,sigmaElast,.true.,2.3,.false.)

       write(301,'(5F12.3)') plab,sigmaTot,sigmaElast

    end do
    close(301)



  end subroutine testRhoNuc

end program test
