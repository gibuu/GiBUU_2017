program GlauberAbs

! This program implements simple Glauber model absporption formulae.

  use inputGeneral
  use version, only: PrintVersion
  use particleProperties, only: initParticleProperties
  use nucleusDefinition
  use nucleus, only : getTarget
  use dichteDefinition
  use densityStatic
  use output, only: Write_ReadingInput

  implicit none


  type(tNucleus),pointer :: targetNuc
  type(dichte)   :: dens

  integer :: ib, ibMax, iz
  real :: b, bMin=  0.0, bMax=14.0, db=0.05 !db=0.01
  real :: z, zMin=-14.0, zMax=14.0, dz=0.05 !dz=0.01

  real,dimension(0:5) :: hS,hIntegrand

  real, dimension(1:3) :: r
  real hh1

  real :: PartM,PartG,PartE, Part_gamma,Part_beta

  NAMELIST /Glauber/ PartM,PartG,PartE
  integer :: ios

  real :: sigma1=3.0, sigma2=6.0 ! in fm^2 = 10*mb
!  real :: sigma1=2.5, sigma2=5.0 ! in fm^2 = 10*mb

  call PrintVersion

!  call readInputGeneral
  call initParticleProperties

  targetNuc => getTarget()

  call Write_ReadingInput('Glauber',0)
  rewind(5)
  read(5,nml=Glauber,iostat=ios)
  call Write_ReadingInput('Glauber',0,ios)

  write(*,'(A,f12.3)') "Particle Mass:  ", PartM
  write(*,'(A,f12.3)') "Particle Gamma: ", PartG
  write(*,'(A,f12.3)') "Particle Energy:", PartE
  write(*,*)
  Part_gamma = PartE/PartM
  Part_beta = sqrt(1-1/Part_gamma**2)
  write(*,'(A,f12.3)') "         beta:  ", Part_beta
  write(*,'(A,f12.3)') "         gamma: ", Part_gamma

  call Write_ReadingInput('Glauber',1)

  hS = 0.0

! Everything is a impact parameter integration:

  ibMax = nint((bMax-bMin)/db)
  do ib=0,ibMax
     b = bMin+ib*db

     write(*,'("ib=",2i7.4)') ib,ibMax

     hIntegrand = 0.0
     do iz=0,nint((zMax-zMin)/dz)
        z = zMin+iz*dz

        r = (/b, 0.0, z/)
        dens = staticDensity(r,targetNuc)

        if (dens%baryon(0).lt.1e-5) then
           if (z<0) cycle
           exit
        end if

        !-------------

        hIntegrand(0)=hIntegrand(0)+dens%baryon(0)

        !-------------

        hh1 = Integrand1()
        hIntegrand(1)=hIntegrand(1)+dens%baryon(0)*exp(-sigma1*hh1)

        !-------------

        hIntegrand(2)=hIntegrand(2)+dens%baryon(0)*Integrand2()

        !-------------

        hIntegrand(3)=hIntegrand(3)+dens%baryon(0)*Integrand3()

        !-------------

        hIntegrand(4)=hIntegrand(4)+dens%baryon(0)*Integrand4()

        !-------------
        
     end do
     hIntegrand = hIntegrand * dz * b
     hS = hS + hIntegrand

  end do
  hS = hS * db * (2*3.14159d0)

  write( *,'(5f13.5,1P,30e13.5)') PartM,PartG,PartE,Part_beta,Part_gamma,hS
  write(21,'(5f13.5,1P,30e13.5)') PartM,PartG,PartE,Part_beta,Part_gamma,hS

  write( *,'(5f13.5,1P,30e13.5)') PartM,PartG,PartE,Part_beta,Part_gamma,hS/hS(0)
  write(22,'(5f13.5,1P,30e13.5)') PartM,PartG,PartE,Part_beta,Part_gamma,hS/hS(0)

contains
  !
  ! \int_z^\infty{dz' rho(z')}
  !
  real function Integrand1()
    integer :: iZA
    real :: zA
    real, dimension(1:3) :: rA
    type(dichte)   :: densA

    Integrand1 = 0.0
    do izA=iz,nint((zMax-zMin)/dz)
       zA= zMin+izA*dz
       rA = (/b, 0.0, zA/)
       densA = staticDensity(rA,targetNuc)
       if (densA%baryon(0).lt.1e-5) then
           if (zA>0) exit
        end if
       Integrand1 = Integrand1+densA%baryon(0)
    end do
    Integrand1 = Integrand1 * dz
    return
  end function Integrand1

  !
  ! \int_0^\infty{dt \Gamma/\gamma \exp[-t\Gamma/\gamma]
  !                    \exp[-\sigma_1\int_z^{z+\beta t} dz' \rho(z')] }
  !

  real function Integrand2()
    integer :: iZA, iT, iZAmax
    real :: zA
    real, dimension(1:3) :: rA
    type(dichte)   :: densA
    real :: hh,hh1
    real :: T, TMin=  0.0, TMax=50.0, dT=0.01 ! time

    real, dimension(0:10000), save :: Arr
    real :: term1, term2

! tabulate integral

    hh = 0
    iZAmax = nint((zMax-zMin)/dz)
    do izA=iz,iZAmax
       zA= zMin+izA*dz
       rA = (/b, 0.0, zA/)
       densA = staticDensity(rA,targetNuc)
       hh = hh+densA%baryon(0)*dz
       Arr(iZA) = exp(-sigma1*hh)
       if (zA>0) then
          if (densA%baryon(0).lt.1e-5) then
             iZAmax = iZA
             exit
          end if
       endif
    end do
    

    TMax = (iZAmax-iz)*dz/Part_beta

!    write(*,*) 'tMax = ',tmax,iZAmax,iz

    term1 = Arr(iZAmax) * exp(-PartG*TMax/(0.197*Part_gamma))


    term2 = 0.0
    do iT=0,nint((TMax-TMin)/dT)
       T = TMin+iT*dT
       
       hh1 = exp(-PartG*T/(0.197*Part_gamma))
       if (hh1.lt.1e-5) exit

       izA = nint((z+Part_beta*T-zMin)/dz)
       if (iZA>iZAmax) then
          write(*,*) '(iZA>iZAmax):',iZA,iZAmax
          stop
       endif

       term2 = term2 + hh1*Arr(izA)

    end do

    term2 = term2 * PartG*dT/(0.197*Part_gamma)

    Integrand2 = term1 + term2

    return
  end function Integrand2

  !
  ! \int_0^\infty{dt \Gamma/\gamma \exp[-t\Gamma/\gamma]
  !                    \exp[-\sigma_1\int_z^{z+\beta t} dz' \rho(z')] 
  !                    \exp[-\sigma_2\int_{z+\beta t}^\infty dz' \rho(z')] }
  !

  real function Integrand3()
    integer :: iZA, iT, iZAmax
    real :: zA
    real, dimension(1:3) :: rA
    type(dichte)   :: densA
    real :: hh,hh1
    real :: T, TMin=  0.0, TMax=50.0, dT=0.01 ! time

    real, dimension(0:10000), save :: Arr
    real, dimension(0:10000), save :: Arr2
    real :: term1, term2

    ! tabulate integral:
    ! * Arr(izA) : exp[-sigma1*Int_iz^izA rho(i)]
    ! * Arr2(izA): exp[-sigma2*Int_izA^izAmax rho(i)]

    hh = 0
    iZAmax = nint((zMax-zMin)/dz)
    do izA=iz,iZAmax
       zA= zMin+izA*dz
       rA = (/b, 0.0, zA/)
       densA = staticDensity(rA,targetNuc)
       hh = hh+densA%baryon(0)*dz
       Arr(iZA) = exp(-sigma1*hh)
       Arr2(izA) = hh
       if (zA>0) then
          if (densA%baryon(0).lt.1e-5) then
             iZAmax = iZA
             exit
          end if
       endif
    end do
    do izA=iz,iZAmax
       Arr2(izA) = exp(-sigma2*(Arr2(izAmax) - Arr2(izA)))
    end do

    TMax = (iZAmax-iz)*dz/Part_beta

!    write(*,*) 'tMax = ',tmax,iZAmax,iz

    term1 = Arr(iZAmax) * exp(-PartG*TMax/(0.197*Part_gamma))


    term2 = 0.0
    do iT=0,nint((TMax-TMin)/dT)
       T = TMin+iT*dT
       
       hh1 = exp(-PartG*T/(0.197*Part_gamma))
       if (hh1.lt.1e-5) exit

       izA = nint((z+Part_beta*T-zMin)/dz)
       if (iZA>iZAmax) then
          write(*,*) '(iZA>iZAmax):',iZA,iZAmax
          stop
       endif

       term2 = term2 + hh1*Arr(izA)*Arr2(izA)

    end do

    term2 = term2 * PartG*dT/(0.197*Part_gamma)

    Integrand3 = term1 + term2

    return
  end function Integrand3


  !
  ! \int_0^\infty{dt \Gamma/\gamma \exp[-t\Gamma/\gamma]
  !                    \exp[-\sigma_1\int_z^{z+\beta (t+tF/2)} dz' \rho(z')] 
  !                    \exp[-\sigma_2\int_{z+\beta (t+tF/2)}^\infty dz' \rho(z')] }
  !

  real function Integrand4()
    integer :: iZA, iT, iZAmax
    real :: zA
    real, dimension(1:3) :: rA
    type(dichte)   :: densA
    real :: hh,hh1
    real :: T, TMin=  0.0, TMax=50.0, dT=0.01 ! time

    real, dimension(0:10000), save :: Arr
    real, dimension(0:10000), save :: Arr2
    real :: term1, term2

    ! tabulate integral:
    ! * Arr(izA) : exp[-sigma1*Int_iz^izA rho(i)]
    ! * Arr2(izA): exp[-sigma2*Int_izA^izAmax rho(i)]

    hh = 0
    iZAmax = nint((zMax-zMin)/dz)
    do izA=iz,iZAmax
       zA= zMin+izA*dz
       rA = (/b, 0.0, zA/)
       densA = staticDensity(rA,targetNuc)
       hh = hh+densA%baryon(0)*dz
       Arr(iZA) = exp(-sigma1*hh)
       Arr2(izA) = hh
       if (zA>0) then
          if (densA%baryon(0).lt.1e-5) then
             iZAmax = iZA
             exit
          end if
       endif
    end do
    do izA=iz,iZAmax
       Arr2(izA) = exp(-sigma2*(Arr2(izAmax) - Arr2(izA)))
    end do

    TMax = (iZAmax-iz)*dz/Part_beta

!    write(*,*) 'tMax = ',tmax,iZAmax,iz

    term1 = Arr(iZAmax) * exp(-PartG*TMax/(0.197*Part_gamma))


    term2 = 0.0
    do iT=0,nint((TMax-TMin)/dT)
       T = TMin+iT*dT
       
       hh1 = exp(-PartG*T/(0.197*Part_gamma))
       if (hh1.lt.1e-5) exit

       izA = nint((z+Part_beta*(T+PartE/2)-zMin)/dz)
       if (iZA>iZAmax) then
          iZA=iZAmax
!          write(*,*) '(iZA>iZAmax):',iZA,iZAmax
!          stop
       endif

       term2 = term2 + hh1*Arr(izA)*Arr2(izA)

    end do

    term2 = term2 * PartG*dT/(0.197*Part_gamma)

    Integrand4 = term1 + term2

    return
  end function Integrand4


end program GlauberAbs
