program test
  use particleProperties, only: InitParticleProperties, hadron
  use PIL_rhoDiff,  only: PIL_rhoDiffractive_PUT
  use MediumDefinition, only: medium, vacuum
  use IDTable
  use particleDefinition
  use winkelVerteilung
  use constants, only: mN, mPi
  use inputGeneral, only: readInputGeneral

  implicit none

  real                            :: srts            ! sqrt(s)
  type(particle), dimension (1:2) :: pairIn          ! incoming particles a and b
  real, dimension(1:3)            :: betaToCM        ! beta for boost to CM-Frame
  type(particle), dimension (1:2) :: pairOut         ! outgoing particles c and d
  real, dimension(1:3)            :: pscatt
  integer :: i

  write(*,*)"************************************************"
  write(*,*)"Test WinkelVerteilung"
  write(*,*)"************************************************"
  write(*,*)

  write(*,*) 'Initializing databench for particles'
  call readInputGeneral
  call InitParticleProperties

  betaToCM=0.

!    call piN
!     call NN_NDelta
!    call delta_dec
!    call rho_pipi
!     call omegaN
!     call NN_NR
  call NN_BYK

contains

    !***********************************************************************************************************+

    subroutine piN
      write(*,*) ' Testing pion nucleon -> pion nucleon'

      pairOut%ID   = (/nucleon, pion/)
      pairOut%mass = (/mN, mPi/)

      pairIn%ID     = (/nucleon,pion/)
      pairIN%mass   = (/mN,mPi/)
      pairIN%charge = 0

      pairIN(1)%momentum(0:3) = (/mN,0.,0.,0./)

      pairIN(2)%momentum(0)   = sqrt(mPi**2+0.1**2)
      pairIN(2)%momentum(1:3) = (/0.,0.,0.1/)


      srts=sqrts(pairIn)


      pairIN(2)%momentum(1:3)=(/0.05,0.,0./)
      write(*,*) ' Result in fort.10 for pion in x-direction for sqrt(s)=' ,srts
      Do i=1,1000
         pscatt = winkel (pairIn, pairOut, srts, betaToCM, vacuum)
         write(10,*) pscatt
      End do

      pairIN(2)%momentum(1:3)=(/0.,0.05,0./)
      write(*,*) ' Result in fort.11 for pion in y-direction for sqrt(s)=' ,srts
      Do i=1,1000
         pscatt = winkel (pairIn, pairOut, srts, betaToCM, vacuum)
         write(11,*) pscatt
      End do

      write(*,*) ' Result in fort.12 for pion in z-direction for sqrt(s)=' ,srts
      pairIN(2)%momentum(1:3)=(/0.,0.,0.05/)
      Do i=1,1000
         pscatt = winkel (pairIn, pairOut, srts, betaToCM, vacuum)
         write(12,*) pscatt
      End do


    end subroutine piN


    !***********************************************************************************************************+

    subroutine NN_NDelta
      use constants, only: pi
      use histf90
      use lorentzTrafo, only: lorentzCalcBeta
      use dimi, only: dimiIntegrated
      use master_2body, only: setKinematics

      type(histogram) :: hist_cost, hist_theta, hist_y
      integer, parameter :: Nevt = 100000
      real :: plab,cost,sig,m,E,y
      real, dimension(1:3) :: betaToLRF = 0.
      logical :: success
      integer :: Nsucc = 0

      call CreateHist (hist_cost, 'angular distribution [cos(theta)] for N N -> N Delta', -1., 1., 0.01)
      call CreateHist (hist_theta, 'angular distribution [theta] for N N -> N Delta', 0., 180., 1.)
      call CreateHist (hist_y, 'rapidity distribution for N N -> N Delta', -2., 2., 0.02)

      write(*,*) ' Testing N N -> N Delta'

      pairOut%ID   = (/nucleon, Delta/)
      pairOut%mass = (/mN, hadron(Delta)%mass/)

      pairIn%ID     = (/nucleon, nucleon/)
      pairIN%mass   = (/mN, mN/)
      pairIN%charge = 1

      pairIN(1)%momentum(0:3) = (/mN,0.,0.,0./)

      plab = 1.98   ! 1.66 ! (Dmitriev)
      pairIN(2)%momentum(0)   = sqrt(mN**2+plab**2)
      pairIN(2)%momentum(1:3) = (/0.,0.,plab/)

      srts=sqrts(pairIn)

      betaToCM = lorentzCalcBeta(pairIn(1)%momentum(0:3)+pairIn(2)%momentum(0:3))

      sig = dimiIntegrated(srts)

      print *,"p_lab=",plab,"sqrt(s)=",srts,"sig_tot=",sig

      Do i=1,Nevt
         ! (1) produce angle in CM frame
         pscatt = winkel (pairIn, pairOut, srts, betaToCM, vacuum)  ! pscatt is in CM frame
         cost = pscatt(3)/sqrt(sum(pscatt*pscatt))
         m = hadron(Delta)%mass                                   ! for simplicity: assume constant mass
         E = sqrt(m**2 + sum(pscatt(1:3)**2))
         y = 0.5*log((E+pscatt(3))/(E-pscatt(3)))                 ! rapidity
         call AddHist (hist_cost, cost, sig/float(Nevt))
         call AddHist (hist_theta, acos(cost)*180./pi, sig/float(Nevt))
         call AddHist (hist_y, y, sig/float(Nevt))

         ! (2) produce 'full event' in lab frame
         call setKinematics (srtS, srtS, betaToLRF, betaToCM, vacuum, pairIn, pairOut, success)
         if (success) then
           pscatt = pairOut(2)%momentum(1:3)                      ! pscatt is in LAB frame
           E = pairOut(2)%momentum(0)
           cost = pscatt(3)/sqrt(sum(pscatt*pscatt))
           y = 0.5*log((E+pscatt(3))/(E-pscatt(3)))               ! rapidity
           call AddHist (hist_cost, cost, 0., sig/float(Nevt))
           call AddHist (hist_theta, acos(cost)*180./pi, 0., sig/float(Nevt))
           call AddHist (hist_y, y, 0., sig/float(Nevt))
           Nsucc = Nsucc + 1
         end if
      End do

      print *,Nevt,Nsucc

      call WriteHist(hist_cost,  file="NN_NDelta_cost.txt")
      call WriteHist(hist_theta, file="NN_NDelta_theta.txt")
      call WriteHist(hist_y,     file="NN_NDelta_y.txt")

    end subroutine NN_NDelta

    !***********************************************************************************************************

    subroutine omegaN
      ! This routine tests the angular distribution of "gamma N -> omega N".
      ! The setup is such that the results can be compared directly to the SAPHIR data,
      ! cf. J. Barth et al., Low-energy photoproduction of omega mesons, EPJ A 18 (2003), 117-127.
      use histf90
      use constants, only: pi, mN
      use twoBodyTools, only: pcm
      use random, only: rnFlat,rnExp
      use inputGeneral, only: path_To_Input
      type(histogram) :: hist_theta,hist_t
      real :: theta,t_prime,dt_max,p_i,p_f,photon_energy,E_lo,E_mid,E_hi,E_thres
      integer,parameter :: Nevt = 500000
      real,parameter :: dE = 0.025
      integer :: iParam,j,k=0
      character(len=100) :: suffix
      real, dimension(1:23) :: x,y

      write(*,*) ' Testing gamma N -> omega N'

      pairOut%ID   = (/omegaMeson,nucleon/)
      pairOut%mass = (/hadron(omegaMeson)%mass,mN/)

      pairIn%ID     = (/photon,nucleon/)
      pairIN%mass   = (/0.,mN/)
      pairIN%charge = (/0,1/)

      pairIN(2)%momentum = (/mN,0.,0.,0./)  ! nucleon at rest

      E_thres = hadron(omegaMeson)%mass/(2*mN) * (2*mN+hadron(omegaMeson)%mass)  ! threshold for omega production

      ! read total cross section data from file (for normalization)
      open(202,file=trim(path_To_Input)//"/gammaN_omegaN_saphir.dat")
      Do i=lBound(y,dim=1),uBound(y,dim=1)
         read(202,*) x(i),y(i)
      end do
      close(202)

      do j=0,59  ! loop over photon energy bins

        E_lo = 1.1 +j*dE  ! lower bound of bin

        if ((mod(j,4)==0) .or. (E_lo<1.7 .and. mod(j,2)==0) .or. E_lo<1.2) then  ! SAPHIR binning!

           k=k+1
           ! determine upper bound of bin
           if (E_lo>1.7) then
              E_hi = E_lo + 4*dE
           else if (E_lo>1.2) then
              E_hi = E_lo + 2*dE
           else
              E_hi = E_lo + dE
           end if
           E_mid = (E_hi+E_lo)/2.  ! bin center

           ! go to upper bound of bin to determine dt_max (for histogram bounds)
           pairIN(1)%momentum = (/E_hi,0.,0.,E_hi/)
           srts = sqrts(pairIn)
           p_i = pcm (srts, pairIn(1)%mass , pairIn(2)%mass)
           p_f = pcm (srts, pairOut(1)%mass, pairOut(2)%mass)
           ! via PDG (38.32):
           dt_max = 4.*p_i*p_f                                                                      !  = (p_i+p_f)**2 - (p_i-p_f)**2

           call CreateHist (hist_theta, 'angular distribution (theta) for gamma N -> omega N', 0., 180., 1.)
           call CreateHist (hist_t,     'angular distribution (t) for gamma N -> omega N'    , 0., dt_max , dt_max/20.)

           ! go to center of bin to determine p_i / p_f per bin ?!?
!           pairIN(1)%momentum = (/E_mid,0.,0.,E_mid/)
!           srts = sqrts(pairIn)
!           p_i = pcm (srts, pairIn(1)%mass , pairIn(2)%mass)
!           p_f = pcm (srts, pairOut(1)%mass, pairOut(2)%mass)

           print '(A,8F9.5)',"***",E_lo,E_mid,E_hi,dt_max,p_i,p_f,x(k),y(k)

        else
           cycle
        end if


        Do i=1,Nevt  ! loop over MC events

           photon_energy = rnFlat (max(E_thres,E_lo), E_hi)             ! choose random photon energy
           !photon_energy = rnExp (-1.,max(E_thres,E_lo), E_hi)          ! alternative: exponential distribution of photon energy ?!?

           pairIN(1)%momentum = (/photon_energy,0.,0.,photon_energy/)
           srts = sqrts(pairIn)

           ! calculate p_i and p_f per event
           p_i = pcm (srts, pairIn(1)%mass , pairIn(2)%mass)
           p_f = pcm (srts, pairOut(1)%mass, pairOut(2)%mass)

           pscatt = winkel (pairIn, pairOut, srts, betaToCM, vacuum)
           theta   = acos(pscatt(3)/sqrt(sum(pscatt*pscatt)))       ! pscatt is in CM frame !
           t_prime = 4.*p_i*p_f*sin(theta/2.)**2                    ! t_prime = abs(t - t_min) = t_min - t
           call AddHist (hist_theta, theta*180./pi, 1./sin(theta))
           call AddHist (hist_t,     t_prime      , 1.)

        End do

        iParam = getIParam()
        write(suffix,'(F6.4,A,i1,A)') E_mid,"GeV_",iParam,".dat"
        call WriteHist(hist_theta, 100, file="omegaN_theta_"//suffix,mul=y(k)/float(Nevt))
        call WriteHist(hist_t    , 100, file="omegaN_t_"//suffix,mul=y(k)/float(Nevt))

      end do

    end subroutine omegaN

    !***********************************************************************************************************+

    subroutine delta_dec
      write(*,*) ' Testing Delta -> pion Nucleon'

      pairOut%ID=(/nucleon, pion/)
      pairOut%mass=(/mN, mPi/)

      pairIn%ID= (/delta,0/)
      pairIN(1)%mass=hadron(delta)%mass
      pairIN(1)%charge = 0
      pairIN(1)%momentum(0)=hadron(delta)%mass

      srts=sqrts(pairIN)

      write(*,*) ' Result in fort.100 for mesonMomentum in z-direction for sqrt(s)=' ,srts
      pairIN(1)%momentum(1:3)=(/0.,0.,0.1/)
      Do i=1,10000
         pscatt = winkel (pairIn, pairOut, srts, betaToCM, vacuum)
         write(100,*) pscatt
      End do

      write(*,*) ' Result in fort.101 for mesonMomentum in x-direction for sqrt(s)=' ,srts
      pairIN(1)%momentum(1:3)=(/0.1,0.,0./)
      Do i=1,10000
         pscatt = winkel (pairIn, pairOut, srts, betaToCM, vacuum)
         write(101,*) pscatt
      End do

      write(*,*)
    end subroutine delta_dec

    !**************************************************************************************************************

    subroutine rho_pipi
      write(*,*) ' Testing rho -> pion pion'

      pairOut%ID=(/pion, pion/)
      pairOut%mass=(/mPi, mPi/)

      pairIn(1)%ID=rho
      pairIN(2)%ID=0
      pairIN(1)%mass=hadron(rho)%mass
      pairIN(2)%mass=0
      pairIN(1)%charge=0
      pairIN(2)%Charge=0
      pairIN(1)%momentum(1:3)=(/1.,0.,0./)
      pairIN(1)%momentum(0)=sqrt(hadron(rho)%mass**2+Dot_product(pairIN(1)%momentum(1:3),pairIN(1)%momentum(1:3)))
      pairIn(1)%number=5

      srts=sqrts(pairIN)

      write(*,*) ' Result in fort.199 no information.' ,srts
      ! call PIL_rhoDiffractive_PUT(pairIn(1)%number,.false.)
      Do i=1,100
         pscatt = winkel (pairIn, pairOut, srts, betaToCM, vacuum)
         write(199,*) pscatt
      End do

      write(*,*) ' Result in fort.200 for diffractive=.true.' ,srts
      !call PIL_rhoDiffractive_PUT(pairIn(1)%number,.true.)
      Do i=1,100
         pscatt = winkel (pairIn, pairOut, srts, betaToCM, vacuum)
         write(200,*) pscatt
      End do

      write(*,*) ' Result in fort.201 for for diffractive=.false.' ,srts
      !call PIL_rhoDiffractive_PUT(pairIn(1)%number,.false.)
      Do i=1,100
         pscatt = winkel (pairIn, pairOut, srts, betaToCM, vacuum)
         write(201,*) pscatt
      End do

    end subroutine rho_pipi

    !**************************************************************************************************************

    subroutine NN_NR
      use histf90
      use constants, only: pi, mN
      use twoBodyTools, only: pcm

      type(histogram) :: hist_theta,hist_t
      real :: theta,t_prime,dt_max,p_i,p_f
      integer, parameter :: Nevt = 1000000
      integer, parameter :: resID = P11_1440 !D13_1520
      real, parameter :: Ekin = 3.5

      write(*,*) ' Testing N N -> N R'

      pairIn%ID     = nucleon
      pairIN%mass   = mN
      pairIN%charge = 1

      pairOut%ID   = (/nucleon,resID/)
      pairOut%mass = (/mN,hadron(resID)%mass/)

      pairIN(1)%momentum = (/mN+Ekin,0.,0.,sqrt(Ekin**2+2*mN*Ekin)/)
      pairIN(2)%momentum = (/mN,0.,0.,0./)  ! nucleon at rest

      srts = sqrts(pairIn)

      ! calculate p_i and p_f
      p_i = pcm (srts, pairIn(1)%mass , pairIn(2)%mass)
      p_f = pcm (srts, pairOut(1)%mass, pairOut(2)%mass)

      dt_max = 4.*p_i*p_f                                                                      !  = (p_i+p_f)**2 - (p_i-p_f)**2

      call CreateHist (hist_theta, 'angular distribution (theta) for N N -> N R', 0., 180., 1.)
      call CreateHist (hist_t,     'angular distribution (t) for N N -> N R'    , 0., dt_max , dt_max/50.)

      Do i=1,Nevt  ! loop over MC events

         pscatt  = winkel (pairIn, pairOut, srts, betaToCM, vacuum)
         theta   = acos(pscatt(3)/sqrt(sum(pscatt*pscatt)))       ! pscatt is in CM frame !
         t_prime = 4.*p_i*p_f*sin(theta/2.)**2                    ! t_prime = abs(t - t_min) = t_min - t
         call AddHist (hist_theta, theta*180./pi, 1./sin(theta))
         call AddHist (hist_t,     t_prime      , 1.)

      End do

      call WriteHist(hist_theta, 100, file="NN_NR_theta.dat")
      call WriteHist(hist_t    , 100, file="NN_NR_t.dat")

    end subroutine NN_NR

    !**************************************************************************************************************

    subroutine NN_BYK
      use histf90
      use constants, only: pi, mN
      use nBodyPhaseSpace, only: momenta_in_3BodyPS, momenta_in_3Body_BYK

      type(histogram) :: hist_theta,hist_Ek
      real :: theta,Ek_max,pK,Ek
      integer, parameter :: Nevt = 1000000
      integer, parameter :: Bid = nucleon ! Delta
      integer, parameter :: Yid = Lambda ! SigmaResonance
      real, parameter :: Ekin = 3.5
      real, dimension(1:3) :: masses
      real, dimension(3,3) :: mom

      write(*,*) ' Testing N N -> B Y K'

      pairIn%ID     = nucleon
      pairIN%mass   = mN
      pairIN%charge = 1

      masses(1:3) = (/ hadron(Bid)%mass, hadron(Yid)%mass, hadron(Kaon)%mass /)

      pairIN(1)%momentum = (/mN+Ekin,0.,0.,sqrt(Ekin**2+2*mN*Ekin)/)
      pairIN(2)%momentum = (/mN,0.,0.,0./)  ! nucleon at rest

      srts = sqrts(pairIn)

      Ek_max = srts - sum(masses(1:2))

      call CreateHist (hist_theta, 'angular distribution (theta) for N N -> B Y K', 0., 180., 1.)
      call CreateHist (hist_Ek,    'angular distribution (t) for N N -> B Y K'    , 0., Ek_max , Ek_max/100.)

      Do i=1,Nevt  ! loop over MC events

         ! (1) phase space
         mom = momenta_in_3BodyPS (srts, masses)
         pK = sqrt(sum(mom(1:3,3)**2))
         theta   = acos(mom(3,3)/pK)       ! mom is in CM frame !
         Ek      = sqrt(masses(3)**2 + pK**2)
         call AddHist (hist_theta, theta*180./pi, y=1./sin(theta))
         call AddHist (hist_Ek,    Ek, y=1.)

         ! (2) BYK model
         mom = momenta_in_3Body_BYK (srts, pairIN(1)%momentum(1:3), masses)
         pK = sqrt(sum(mom(1:3,3)**2))
         theta   = acos(mom(3,3)/pK)       ! mom is in CM frame !
         Ek      = sqrt(masses(3)**2 + pK**2)
         call AddHist (hist_theta, theta*180./pi, y=0., y2=1./sin(theta))
         call AddHist (hist_Ek,    Ek, y=0., y2=1.)

      End do

      call WriteHist(hist_theta, 100, file="NN_BYK_theta.dat")
      call WriteHist(hist_Ek   , 100, file="NN_BYK_Ek.dat")

    end subroutine NN_BYK

end program test
