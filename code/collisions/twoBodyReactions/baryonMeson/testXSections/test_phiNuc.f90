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
  write(*,*) '**************************************************'

  write(*,*) '**************************************************'
  write(*,*) 'Testing Phi Nuk -> '
  write(*,*) '**************************************************'

  call testPhiNuc


  contains


  subroutine testPhiNuc
    use IDTABLE
    use mediumDefinition
    use particleDefinition
    use particleProperties
    use phiNucleon
    use preEventDefinition
    implicit none
    real                :: srts                  ! sqrt(s) in the process
    type(particle),dimension(1:2)   :: partIn        ! colliding particles
    type(medium)                      :: mediumAtColl    ! Medium informations at the position of the collision

    logical                         :: plotFlag          ! Switch on plotting of the  Xsections
    real,dimension(0:3)       :: momLRF        ! Total Momentum in LRF
    type(preEvent),dimension(1:3) :: partOut     ! colliding particles
    real                            :: sigmaTot         ! total Xsection
    real                            :: sigmaElast      ! elastic Xsecti
    integer :: i,j,chargeMeson,chargeNuk, k
    real :: plab,ekin,mass
    real :: piN,phiN, pipiN,r_p13
    real :: dens
    integer :: numTries=1000
    real, dimension(3) :: sigmaArr



    NAMELIST /initTest/ chargeMeson, chargeNuk, dens

    chargeMeson = 0
    chargeNuk = 0
    dens = 0.0

    write(*,*)
    write(*,*) '**Initializing testing parameter'
    write(*,*) ' Reading input ....'
    rewind(5)
    read(5,nml=initTest)
    write(*,*) ' Set charge to ', chargeMeson,'.'
    write(*,*) ' Set nuk charge to ',  chargeNuk,'.'
    write(*,*) ' Set density to ',  dens,'.'

    if (dens>0) then
       mediumAtColl%useMedium=.true.
       mediumAtColl%densityProton=dens/2.
       mediumAtColl%densityNeutron=dens/2.
    end if

    partIn(1)%Id=phi
    partIn(1)%charge=chargeMeson
    partIn(1)%mass=hadron(phi)%mass

    partIn(2)%Id=nucleon
    partIn(2)%charge=chargeNuk
    partIn(2)%mass=0.938
    partIn(2)%momentum=(/partIn(2)%mass,0.,0.,0./)
    partIn(2)%velocity=(/0.,0.,0./)

    open(301,file='PhiNucleonXsections.dat',status='unknown')
    write(301,*) '# nukCharge=',chargeNuk,'    mesonCharge=', chargeMeson
    write(301,*) '# plab,mass,srts,sigmaTot,sigmaElast,sigmaArr'

    do j=1,100
       mass=j*0.02

       do i=1,100


!          plab=i*0.02
!          srts=SQRT((SQRT(mass**2+plab**2)+hadron(nucleon)%mass)**2-plab**2)

          srts = 0.940 + i*0.02
          plab =srts**4 + mass**4 + hadron(nucleon)%mass**4 &
               - 2 * mass**2 * hadron(nucleon)%mass**2 &
               - 2 * mass**2 * srts**2 &
               - 2 * srts**2 * hadron(nucleon)%mass**2

          if (plab > 0) then
             plab = sqrt(plab)/(2*hadron(nucleon)%mass)

             partIn(1)%mass=mass
             partIn(1)%momentum=(/sqrt(mass**2+plab**2),plab,0.,0./)
             partIn(1)%velocity=(/partIn(1)%momentum(1)/partIn(1)%momentum(0),0.,0./)


             momLRF=(/SQRT(mass**2+plab**2)+hadron(nucleon)%mass ,plab,0.,0./)
             ekin=SQRT(mass**2+plab**2)-mass

          !       call phiNuc(srts,partIn,mediumAtColl,momLRF,partOut,sigmaTot,sigmaElast,.true.,2.2,.true.)
          !       write(301,'(5F12.3)') plab,mass,srts,sigmaTot,sigmaElast
             call phiNuc(srts,partIn,mediumAtColl,partOut,sigmaTot,sigmaElast,.true.,2.2,.false.,sigmaArr)

          else
             sigmaTot = 0.
             sigmaElast = 0.
             sigmaArr = 0.
          end if


          write(301,'(50F12.3)') plab,mass,srts,sigmaTot,sigmaElast,sigmaArr

       end do
       write(301,*)
       write(*,*) i
    end do
    close(301)

    piN=0.
    phiN=0.
    pipiN=0.
    do i=1,numTries
       plab=0.2
       srts=SQRT((SQRT(hadron(phi)%mass**2+plab**2)+hadron(nucleon)%mass)**2-plab**2)
       partIn(1)%momentum=(/sqrt(partIn(1)%mass**2+plab**2),plab,0.,0./)
       partIn(1)%velocity=(/partIn(1)%momentum(1)/partIn(1)%momentum(0),0.,0./)

       momLRF=(/SQRT(hadron(phi)%mass**2+plab**2)+hadron(nucleon)%mass ,plab,0.,0./)
       call phiNuc(srts,partIn,mediumAtColl,partOut,sigmaTot,sigmaElast,.true.,2.3,.false.)
       if ((partOut(1)%ID.eq.pion).and.(partOut(2)%ID.eq.pion).and.(partOut(3)%ID.eq.nucleon)) then
          pipiN= pipiN+sigmatot
       else if ((partOut(1)%ID.eq.phi).and.(partOut(2)%ID.eq.nucleon)) then
          phiN=phiN+sigmatot
!       else if ((partOut(1)%ID.eq.P13)) then
!          R_p13=R_p13+sigmatot
       else if ((partOut(1)%ID.eq.pion).and.(partOut(2)%ID.eq.nucleon)) then
          piN= piN+sigmatot
       end if
    end do
    write(*,*) 'Xsections at plab= 0.2 GeV'
    write(*,*) 'Simga for P13. :', r_p13/float(numTries)
    write(*,*) 'Simga for piN. :', piN/float(numTries)
    write(*,*) 'Simga for phi nucleon prod. :', phiN/float(numTries)
    write(*,*) 'Simga for piPiN prod. :', pipiN/float(numTries)




    partIn(1)%Id=phi
    partIn(1)%charge=chargeMeson
    partIn(1)%mass=hadron(phi)%mass

    partIn(2)%Id=nucleon
    partIn(2)%charge=-chargeNuk
    partIn(2)%antiparticle=.true.
    partIn(2)%mass=0.938
    partIn(2)%momentum=(/partIn(2)%mass,0.,0.,0./)
    partIn(2)%velocity=(/0.,0.,0./)

    open(301,file='PhiNucleonXsections_anti.dat',status='unknown')
    write(301,*) '# nukCharge=',chargeNuk,'    mesonCharge=', chargeMeson
    do i=1,250
       plab=i*0.01

       partIn(1)%momentum=(/sqrt(partIn(1)%mass**2+plab**2),plab,0.,0./)
       partIn(1)%velocity=(/partIn(1)%momentum(1)/partIn(1)%momentum(0),0.,0./)

       srts=SQRT((SQRT(hadron(phi)%mass**2+plab**2)+hadron(nucleon)%mass)**2-plab**2)
       momLRF=(/SQRT(hadron(phi)%mass**2+plab**2)+hadron(nucleon)%mass ,plab,0.,0./)
       ekin=SQRT(hadron(phi)%mass**2+plab**2)-hadron(phi)%mass
       call phiNuc(srts,partIn,mediumAtColl,partOut,sigmaTot,sigmaElast,.true.,2.3,.false.)

       write(301,'(5F8.3)') plab,sigmaTot,sigmaElast
    end do
    close(301)



  end subroutine testPhiNuc

end program test
