
module paramamp
  implicit none

!   Global constants
      real,parameter::hb=197.33
!      real,parameter::hb=1.
      real,parameter::fc=10000. ! Conversion factor from fm**2 to microbarns
      real,parameter::e=0.3027 ! electron charge: sqrt(4.*pi/137.)

!     pion
      real,parameter::mpi=139.6/hb
!       real,parameter::mpi=134.9764/hb
!     rho
      real,parameter::mrho=770./hb   ! mass
      real,parameter::frho=6.14     ! rhopipi coupling
!      real,parameter::frho=6.03     ! rhopipi coupling (my estimate)

!     Nucleon
      real,parameter::mn=939./hb  !mass
      real,parameter::f=1.  !NNpi coupling
      real,parameter::lam2=(1250./hb)**2  !ff for the NNpi vertex (and ResNpi)
      real,parameter::frhonn=2./mpi*mrho  ! NNrho coupling
      real,parameter::lamrho2=(1400./hb)**2  !ff for the NNrho vertex (and ResNrho)
!     (anomalous) magnetic moments (in units of nuclear magnetons)
      real,parameter::mup=1.79
      real,parameter::mun=-1.91

!     Delta(1232)
      real,parameter::mdel=1232./hb ! mass
      real,parameter::fs=2.13  ! DeltaNpi coupling
      real,parameter::fdel=0.802  ! DeltaDeltapi coupling
      real,parameter::fgam=0.12  ! DeltaNgamma coupling (q2=0)

!     N*(1440)
      real,parameter::mns=1440./hb
      real,parameter::ftil=0.477  ! N*Npi coupling
      real,parameter::gnsdel=2.07 ! N*Deltapi coupling
      real,parameter::ctil=-2.29/mpi ! N*Npipi(s wave) coupling
!     Form factors at the photon point
      real,parameter::f2ropp=-0.067/mn
!      real,parameter::f2ropp=-0.078/mn ! my estimate
      real,parameter::f1ropp=-0.287/mn**2 ! I get 0.03-0.3 depending on S(1/2)
      real,parameter::f2ropn=-0.048/mn ! my estimate using A(1/2) from PDG
      real,parameter::f1ropn=-0.02/mn**2 ! my estimate using S(1/2)=0 from NRQM

!     N*(1520)
      real,parameter::mnsp=1520./hb  ! mass
      real,parameter::gnsp0=66./hb   ! Npi partial width (55 % of 120 MeV)
      real,parameter::fnspdel=-1.061 ! one of the N*Deltapi couplings
      real,parameter::gnspdel=0.64   ! the other N*Deltapi coupling
      real,parameter::grho=5.09      ! N*Nrho coupling (paper)
!      real,parameter::grho=5.7      ! N*Nrho coupling (my estimate)
!     Form factors at the photon point
      real,parameter::g1p=0.782/mn     ! my estimate using PGD
      real,parameter::g2p=-0.410/mn**2 ! my estimate using PGD
      real,parameter::g3p=0.091/mn**2 ! from paper
      real,parameter::g1n=-0.14/mn   ! my estimate using PGD
      real,parameter::g2n=-0.06/mn**2 ! my estimate using PGD
      real,parameter::g3n=0. ! ??? find out

!     Delta*(1700)
      real,parameter::mdels=1700./hb ! mass
      real,parameter::gds0=45./hb    ! Npi partial width (15 % of 300 MeV)
      real,parameter::fdelsdel=-1.325 ! one of the Delta*Deltapi couplings
      real,parameter::gdelsdel=0.146  ! the other Delta*Deltapi coupling
      real,parameter::gprho=2.6       ! Delta*Nrho coupling
!     Form factors at the photon point
      real,parameter::gp1=-0.26/mn
      real,parameter::gp2=0.27/mn**2
      real,parameter::gp3=0.  !??????????????????

end module paramamp
