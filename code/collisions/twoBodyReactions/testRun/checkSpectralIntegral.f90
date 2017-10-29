program checkSpectralIntegral

  use IdTable
  use particleProperties, only: initParticleProperties, hadron
  use constants, only: mPi, mK
  use preEventDefinition
  use inputGeneral, only: readInputGeneral, path_to_input
  use version, only: printVersion
  use mesonWidth, only: fullWidthMeson
  use twoBodyTools, only: pCM

  implicit none

!  integer, save :: srts_maxIndex = 220
  real, parameter :: delta_Srts = 0.01
!  real, parameter :: srtsMin = 0.1


  type(preEvent),dimension(1:2) :: idIn
  type(preEvent),dimension(1:2) :: idOut

  integer :: iCh, iSrts
  real :: srts, srts0
  real :: valI, valC, r, dummy
  logical,dimension(1:2) :: stable
  logical :: flagOK

  integer, parameter, dimension(2,21) :: channelID = reshape((/ &
       pion, pion, & ! 1
       rho, rho, & ! 2
       pion, rho, & ! 3
       pion, eta, & ! 4
       pion, sigmaMeson, & ! 5
       pion, omegaMeson, & ! 6
       pion, etaPrime, & ! 7
       eta, eta, & ! 8
       eta, rho, & ! 9
       eta, sigmaMeson, & ! 10
       eta, omegaMeson, & ! 11
       eta, etaPrime, & ! 12
       rho, sigmaMeson, & ! 13
       rho, omegaMeson, & ! 14
       rho, etaPrime, & ! 15
       sigmaMeson, sigmaMeson, & ! 16
       sigmaMeson, omegaMeson, & ! 17
       sigmaMeson, etaPrime, & ! 18
       omegaMeson, omegaMeson, & ! 19
       omegaMeson, etaPrime, & ! 20
       etaPrime, etaPrime /), & ! 21
       (/2,21/)) ! 21

  do iCH=1,21
     write(*,*) iCH,channelID(:,iCH)
  end do
  stop


  call PrintVersion

  call readInputGeneral
  call initParticleProperties

  dummy = fullWidthMeson(103, 0.9) ! Dummy call for enforce init

  call print105

  idIn(1:2)% ID = (/kaon, kaonbar/)
  do iCh=1,21

     call getOut(iCh)
     stable(1) = (hadron(idOut(1)%id)%width < 1e-3)
     stable(2) = (hadron(idOut(2)%id)%width < 1e-3)
     srts0 = hadron(idOut(1)%id)%minmass +hadron(idOut(2)%id)%minmass

     write(*,'(A,3i5,10g12.5)') "# ",iCh-1,idOut(1)%id,idOut(2)%id
     write(101,'(A,3i5,10g12.5)') "# ",iCh-1,idOut(1)%id,idOut(2)%id
     do iSrts=1,500
        srts = iSrts*0.01

        valC = 0.0
        if (srts>hadron(idOut(1)%id)%mass +hadron(idOut(2)%id)%mass) then
           valC = pCM(srts, hadron(idOut(1)%id)%mass, hadron(idOut(2)%id)%mass, flagOK)
        end if

        if (srts > srts0) then

           if (stable(1).and.stable(2)) then
              valI = 0.0
           else if (stable(1)) then
              valI = calculate1(srts,idOut(1)%id,idOut(2)%id)
           else if (stable(2)) then
              valI = calculate1(srts,idOut(2)%id,idOut(1)%id)
           else
              valI = calculate2(srts,idOut(1)%id,idOut(2)%id)
           end if
        else
           valI = 0.0
        endif
        
        if (valC>0) then 
           r = valI / valC
        else
           r = 0.0
        end if
        

        write(101,'(8g12.5)') srts,valI,valC,r

     end do
     write(101,*)
     write(101,*)

  enddo



contains


  function calculate1(srts, id1, id2) result(integral)

    use IdTable, only: isBaryon
    use particleProperties, only: hadron
    use mesonWidth, only: fullWidthMeson
    use constants, only: pi

    real,    intent(in)  :: srts
    integer, intent(in)  :: id1,id2
    real                 :: integral

    real mass1,mass2,mu1,mu2,pfinal2,pfinal
    integer nmu2,j
    real gamma1,gamma2,gamma2Tot
    real , parameter :: dy=pi/100.
    real dy2,minmu2,maxmu2,ymax2,ymin2,y2
    real spectral2,intfac2

    integral=0.

    mass1=hadron(id1)%mass
    mass2=hadron(id2)%mass
    gamma1=hadron(id1)%width
    gamma2=hadron(id2)%width

!    write(*,*) "xx",mass1,mass2,gamma1,gamma2

    mu1 = mass1

    minmu2=hadron(id2)%minmass
    maxmu2=srts-mu1

    if (maxmu2 < minmu2) return

    ymax2=2.*atan((maxmu2-mass2) / gamma2*2.)
    ymin2=2.*atan((minmu2-mass2) / gamma2*2.)

    nmu2=max(int((ymax2-ymin2)/dy),1)
    dy2=(ymax2-ymin2)/float(nmu2)

    do j=1,nmu2 ! loop over second particle's mass
       y2=ymin2+(float(j)-0.5)*dy2
       mu2=.5*tan(y2/2.)*gamma2+mass2
       mu2=min(max(mu2,minmu2),maxmu2)

       pfinal2=(srts**2-(mu1+mu2)**2)*(srts**2-(mu1-mu2)**2)/ (4.*srts**2)
       if (pfinal2.lt.0) then
!             write(*,*)'pfinal2 lt 0',pfinal2
          pfinal=0.
       else
          pfinal=sqrt(pfinal2)
       end if

       gamma2Tot = fullWidthMeson(ID2, mu2)

       spectral2 = 2./pi * mu2**2 * gamma2tot / ((mu2**2-mass2**2)**2+mu2**2*gamma2tot**2)
       intfac2 = gamma2 / ((mu2-mass2)**2+gamma2**2/4.)

!       write(*,*) spectral2,intfac2

       integral=integral+pfinal*spectral2/intfac2*dy2
    end do


  end function calculate1


  function calculate2(srts, id1, id2) result(integral)

    use IdTable, only: isBaryon
    use particleProperties, only: hadron
    use mesonWidth, only: fullWidthMeson
    use constants, only: pi

    real,    intent(in)  :: srts
    integer, intent(in)  :: id1,id2
    real                 :: integral

    real mass1,mass2,mu1,mu2,pfinal2,pfinal
    integer nmu1,nmu2,i,j
    real gamma1,gamma2,gamma1Tot,gamma2Tot
    real , parameter :: dy=pi/100.
    real dy1,dy2,minmu1,maxmu1,minmu2,maxmu2,ymax1,ymin1,y1,y2
    real spectral1,spectral2,ymax2,ymin2,intfac1,intfac2

    integral=0.

    mass1=hadron(id1)%mass
    mass2=hadron(id2)%mass
    gamma1=hadron(id1)%width
    gamma2=hadron(id2)%width

    minmu1=hadron(id1)%minmass
    maxmu1=srts-hadron(id2)%minmass

    if (maxmu1 < minmu1) return

    ymax1=2.*atan((maxmu1-mass1) / gamma1*2.)
    ymin1=2.*atan((minmu1-mass1) / gamma1*2.)

    nmu1=max(int((ymax1-ymin1)/dy),1)
    dy1=(ymax1-ymin1)/float(nmu1)

    do i=1,nmu1 ! loop over first particle's mass
       y1=ymin1+(float(i)-0.5)*dy1
       mu1=.5*tan(y1/2.)*gamma1+mass1
       mu1=min(max(mu1,minmu1),maxmu1)

       minmu2=hadron(id2)%minmass

       maxmu2=srts-mu1

       ymax2=2.*atan((maxmu2-mass2) / gamma2*2.)
       ymin2=2.*atan((minmu2-mass2) / gamma2*2.)

       nmu2=max(int((ymax2-ymin2)/dy),1)
       dy2=(ymax2-ymin2)/float(nmu2)

       do j=1,nmu2 ! loop over second particle's mass
          y2=ymin2+(float(j)-0.5)*dy2
          mu2=.5*tan(y2/2.)*gamma2+mass2
          mu2=min(max(mu2,minmu2),maxmu2)

          pfinal2=(srts**2-(mu1+mu2)**2)*(srts**2-(mu1-mu2)**2)/ (4.*srts**2)
          if (pfinal2.lt.0) then
!             write(*,*)'pfinal2 lt 0',pfinal2
             pfinal=0.
          else
             pfinal=sqrt(pfinal2)
          end if

          ! Evaluate widht of both particles
          gamma1Tot = fullWidthMeson(ID1, mu1)
          gamma2Tot = fullWidthMeson(ID2, mu2)

          spectral1 = 2./pi * mu1**2 * gamma1Tot / ((mu1**2-  mass1**2)**2+mu1**2*gamma1Tot**2)
          intfac1 = gamma1 / ((mu1-mass1)**2+gamma1**2/4.)

          spectral2 = 2./pi * mu2**2 * gamma2tot / ((mu2**2-mass2**2)**2+mu2**2*gamma2tot**2)
          intfac2 = gamma2 / ((mu2-mass2)**2+gamma2**2/4.)

          integral=integral+pfinal*spectral1*spectral2/intfac1/intfac2*dy1*dy2
       end do
    end do

  end function calculate2


  subroutine getOut(iCh)
    integer, intent(in) :: iCh
    select case(iCh)
    case( 1)
       idOut(1:2)% ID = (/pion, pion/)
    case( 2)
       idOut(1:2)%ID = (/rho, rho/)
    case( 3)
       idOut(1:2)%ID = (/pion, rho/)
    case( 4)
       idOut(1:2)%ID = (/pion, eta/)
    case( 5)
       idOut(1:2)%ID = (/pion, sigmaMeson/)
    case( 6)
       idOut(1:2)%ID = (/pion, omegaMeson/)
    case( 7)
       idOut(1:2)%ID = (/pion, etaPrime/)
    case( 8)
       idOut(1:2)%ID = (/eta, eta/)
    case( 9)
       idOut(1:2)%ID = (/rho, eta/)
    case(10)
       idOut(1:2)%ID = (/eta, sigmaMeson/)
    case(11)
       idOut(1:2)%ID = (/eta, omegaMeson/)
    case(12)
       idOut(1:2)%ID = (/eta, etaPrime/)
    case(13)
       idOut(1:2)%ID = (/rho, sigmaMeson/)
    case(14)
       idOut(1:2)%ID = (/rho, omegaMeson/)
    case(15)
       idOut(1:2)%ID = (/rho, etaPrime/)
    case(16)
       idOut(1:2)%ID = (/sigmaMeson, sigmaMeson/)
    case(17)
       idOut(1:2)%ID = (/sigmaMeson, omegaMeson/)
    case(18)
       idOut(1:2)%ID = (/sigmaMeson, etaPrime/)
    case(19)
       idOut(1:2)%ID = (/omegaMeson, omegaMeson/)
    case(20)
       idOut(1:2)%ID = (/omegaMeson, etaPrime/)
    case(21)
       idOut(1:2)%ID = (/etaPrime, etaPrime/)
    end select
  end subroutine getOut

  subroutine Print105
    use particleProperties, only: hadron
    use mesonWidth, only: fullWidthMeson
    use constants, only: pi

    real :: mass, mass0, gammaTot,spectral,sum, integral
    integer :: iMass, id


    do id=102,107
       mass0=hadron(id)%mass
       sum = 0
       do iMass=1,500
          mass = iMass*0.01
          gammaTot = fullWidthMeson(id, mass)
          spectral = 2./pi * mass**2 * gammaTot / ((mass**2-mass0**2)**2+mass**2*gammaTot**2)
          sum = sum + spectral
          integral = integrate1(mass, id)
          write(102,*) mass,spectral,sum*0.01,integral
       end do
       write(102,*)
       write(102,*)
    end do
  end subroutine Print105

  function integrate1(srts, id2) result(integral)
    
    use IdTable, only: isBaryon
    use particleProperties, only: hadron
    use mesonWidth, only: fullWidthMeson
    use constants, only: pi

    real,    intent(in)  :: srts
    integer, intent(in)  :: id2
    real                 :: integral

    real mass2,mu2
    integer nmu2,j
    real gamma2,gamma2Tot
    real , parameter :: dy=pi/100.
    real dy2,minmu2,maxmu2,ymax2,ymin2,y2
    real spectral2,intfac2

    integral=0.

    mass2=hadron(id2)%mass
    gamma2=hadron(id2)%width

!    write(*,*) "xx",mass1,mass2,gamma1,gamma2


    minmu2=hadron(id2)%minmass
    maxmu2=srts

    if (maxmu2 < minmu2) return

    ymax2=2.*atan((maxmu2-mass2) / gamma2*2.)
    ymin2=2.*atan((minmu2-mass2) / gamma2*2.)

    nmu2=max(int((ymax2-ymin2)/dy),1)
    dy2=(ymax2-ymin2)/float(nmu2)

    do j=1,nmu2 ! loop over second particle's mass
       y2=ymin2+(float(j)-0.5)*dy2
       mu2=.5*tan(y2/2.)*gamma2+mass2
       mu2=min(max(mu2,minmu2),maxmu2)

       gamma2Tot = fullWidthMeson(ID2, mu2)

       spectral2 = 2./pi * mu2**2 * gamma2tot / ((mu2**2-mass2**2)**2+mu2**2*gamma2tot**2)
       intfac2 = gamma2 / ((mu2-mass2)**2+gamma2**2/4.)

!       write(*,*) spectral2,intfac2

       integral=integral+spectral2/intfac2*dy2
    end do


  end function integrate1

  
end program checkSpectralIntegral
