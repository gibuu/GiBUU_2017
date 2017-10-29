program test
  use particleProperties, only: initParticleProperties
  USE inputGeneral, only: readinputGeneral
  implicit none

  write(*,*) '**************************************************'
  write(*,*) 'Initializing database'
  call readinputGeneral
  call initParticleProperties
  write(*,*) '**************************************************'

  write(*,*) '**************************************************'
  write(*,*) 'Testing Omega Nuk -> '
  write(*,*) '**************************************************'

  call testOmegaNuc


contains


  subroutine testOmegaNuc
    use mediumDefinition
    use particleDefinition
    use preEventDefinition
    use IDTABLE, only: nucleon, P13_1900, pion, omegaMeson
    use particleProperties, only: hadron
    use constants, only: mN
    use omegaNucleon, only: omegaNuc

    type(particle) :: teilchenIn(1:2)        ! colliding particles
    type(medium)   :: mediumATcollision      ! Medium informations at the position of the collision
    type(preEvent) :: teilchenOut(1:3)       ! outgoing particles
    real :: srts,sigmaTot,sigmaElast,momentumLRF(0:3),plab,ekin,piN,omegaN, pipiN,r_p13,dens=0.,K_factor=1.
    integer :: i,chargeNuk=1
    integer, parameter :: numTries=10000

    NAMELIST /initTest/ chargeNuk, dens, K_factor

    write (*,*)
    write (*,*) '**Initializing testing parameter'
    write(*,*)  ' Reading input ....'
    rewind(5)
    read(5,nml=initTest)
    write(*,*) ' Nucleon charge: ',  chargeNuk
    write(*,*) ' Density: ',  dens
    write(*,*) ' K-factor: ', K_factor

    mediumAtCollision%useMedium=.true.
    mediumAtCollision%densityProton=dens/2.
    mediumAtCollision%densityNeutron=dens/2.

    !****************
    ! omega-Nucleon !
    !****************

    teilchenIN(1)%Id=omegaMeson
    teilchenIN(1)%charge=0
    teilchenIN(1)%mass=hadron(omegaMeson)%mass

    teilchenIN(2)%Id=nucleon
    teilchenIN(2)%charge=chargeNuk
    teilchenIN(2)%mass=mN
    teilchenIN(2)%momentum=(/teilchenIN(2)%mass,0.,0.,0./)
    teilchenIN(2)%velocity=0.

!     Open(301,file='OmegaNucleonXsections.dat',status='unknown')
!     write(301,*) '# nukCharge=',chargeNuk
    Do i=1,500
      plab=i*0.01
      teilchenIN(1)%momentum=(/sqrt(teilchenIN(1)%mass**2+plab**2),plab,0.,0./)
      teilchenIN(1)%velocity=(/teilchenIN(1)%momentum(1)/teilchenIN(1)%momentum(0),0.,0./)

      srts=SQRT((SQRT(teilchenIN(1)%mass**2+plab**2)+mN)**2-plab**2)
      momentumLRF=(/SQRT(teilchenIN(1)%mass**2+plab**2)+mN,plab,0.,0./)
      ekin=SQRT(teilchenIN(1)%mass**2+plab**2)-teilchenIN(1)%mass
      call omegaNuc(srts,teilchenIN,mediumATcollision,momentumLRF,teilchenOut,sigmaTot,sigmaElast,.true.,2.3,.true.,K_factor)

!       write(301,'(5F12.3)') plab,srts,ekin,sigmaTot,sigmaElast
    end do
!     close(301)

    piN=0.
    omegaN=0.
    pipiN=0.
    Do i=1,numTries
       plab=0.2
       srts=SQRT((SQRT(hadron(omegaMeson)%mass**2+plab**2)+mN)**2-plab**2)
       teilchenIN(1)%momentum=(/sqrt(teilchenIN(1)%mass**2+plab**2),plab,0.,0./)
       teilchenIN(1)%velocity=(/teilchenIN(1)%momentum(1)/teilchenIN(1)%momentum(0),0.,0./)

       momentumLRF=(/SQRT(hadron(omegaMeson)%mass**2+plab**2)+mN,plab,0.,0./)
       call omegaNuc(srts,teilchenIN,mediumATcollision,momentumLRF,teilchenOut,sigmaTot,sigmaElast,.true.,2.3,.false.,K_factor)
       If ((teilchenOut(1)%ID.eq.pion).and.(teilchenOut(2)%ID.eq.pion).and.(teilchenOut(3)%ID.eq.nucleon)) then
          pipiN= pipiN+sigmatot
       else If ((teilchenOut(1)%ID.eq.omegaMeson).and.(teilchenOut(2)%ID.eq.nucleon)) then
          omegaN=omegaN+sigmatot
       else If ((teilchenOut(1)%ID.eq.P13_1900)) then
          R_p13=R_p13+sigmatot
       else If ((teilchenOut(1)%ID.eq.pion).and.(teilchenOut(2)%ID.eq.nucleon)) then
          piN= piN+sigmatot
       end if
    end do
    write(*,*) 'Simga for P13. :', r_p13/float(numTries)
    write(*,*) 'Simga for piN. :', piN/float(numTries)
    write(*,*) 'Simga for omega nucleon prod. :', omegaN/float(numTries)
    write(*,*) 'Simga for piPiN prod. :', pipiN/float(numTries)


    !********************
    ! omega-Antinucleon !
    !********************

    teilchenIN(1)%Id=omegaMeson
    teilchenIN(2)%Id=nucleon
    teilchenIN(1)%charge=0
    teilchenIN(2)%charge=-chargeNuk
    teilchenIN(2)%antiparticle=.true.
    teilchenIN(1)%mass=hadron(omegaMeson)%mass
    teilchenIN(2)%mass=mN
    teilchenIN(2)%momentum=(/teilchenIN(2)%mass,0.,0.,0./)
    teilchenIN(2)%velocity=0.

    Open(301,file='OmegaNucleonXsections_anti.dat',status='unknown')
    write(301,*) '# nukCharge=',chargeNuk
    Do i=1,500
       plab=i*0.01
       teilchenIN(1)%momentum=(/sqrt(teilchenIN(1)%mass**2+plab**2),plab,0.,0./)
       teilchenIN(1)%velocity=(/teilchenIN(1)%momentum(1)/teilchenIN(1)%momentum(0),0.,0./)

       srts=SQRT((SQRT(hadron(omegaMeson)%mass**2+plab**2)+mN)**2-plab**2)
       momentumLRF=(/SQRT(hadron(omegaMeson)%mass**2+plab**2)+mN,plab,0.,0./)
       ekin=SQRT(hadron(omegaMeson)%mass**2+plab**2)-hadron(omegaMeson)%mass
       call omegaNuc(srts,teilchenIN,mediumATcollision,momentumLRF,teilchenOUT,sigmaTot,sigmaElast,.true.,2.3,.false.,K_factor)

       write(301,'(5F8.3)') plab,sigmaTot,sigmaElast

    end do
    close(301)

  end subroutine testOmegaNuc


end program test
