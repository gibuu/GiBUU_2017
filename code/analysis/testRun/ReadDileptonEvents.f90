program ReadDileptonEvents

  ! This is a tool to read dilepton events from a file and produce the usual histograms (mass spectrum etc).
  !
  ! It is useful for external filtering of GiBUU-generated dilepton events:
  ! (1) Make a GiBUU run and write all dilepton events to file ("WriteEvents = .true.")
  ! (2) Filter this file by some external method (Pluto, ROOT-macros, etc), setting the filter result (last column) appropriately.
  ! (3) Run this test case to produce spectra from the filtered events.
  !
  ! Most of the code here is borrowed/duplicated from DileptonAnalysis.f90.
  ! NOTE: The file name is hard-wired and must be set before compiling.

  use constants, only: pi
  use histMC
  implicit none

  character(len=100), save :: filename = "/home/nucleus/jweil/HADES/pluto_v5.34/macros/Dilepton_Events.filtered.dat"

  real(kind=8), dimension(1:2,0:3) :: pout
  real,dimension(0:3)::ptot
  real,dimension(1:3)::p1,p2
  real(kind=8) :: pw,plab,mass,ekin,y,pt,beta_gamma,p1abs,p2abs,angle
  integer :: i,iso,f
  character :: c

  type(histogramMC),save :: msigma,psigma,esigma,ysigma,ptsigma,bgsigma,oasigma!,dsigma

  real :: d!=0.02
  integer,parameter :: nbins=200
  integer,parameter :: nCh=15         ! number of channels
  real :: binsz = 0.005
  real :: projectileEnergy = 1.25

  ! adjust bin size according to available energy
  d=ProjectileEnergy/nbins
  ! create histograms
  call CreateHistMC(msigma,'Dilepton Cross Section dSigma/dM [mubarn/GeV]',0.,3.2,binsz,nCh)
  call CreateHistMC(psigma,'Dilepton Cross Section dSigma/dP_lab [mubarn/GeV]',0.,d*nbins,d,nCh)
  call CreateHistMC(esigma,'Dilepton Cross Section dSigma/dE_kin [mubarn/GeV]',0.,d*nbins,d,nCh)
  call CreateHistMC(ptsigma,'Dilepton Cross Section dSigma/dP_t [mubarn/GeV]',0.,d*nbins,d,nCh)
  call CreateHistMC(ysigma,'Dilepton Cross Section dSigma/dy [mubarn]',-10.,10.,0.1,nCh)
  call CreateHistMC(bgsigma,'Dilepton Cross Section dSigma/d(beta*gamma) [mubarn]',0.,50.,0.1,nCh)
  !call CreateHistMC(oasigma,'Dilepton Cross Section dSigma/dTheta [mubarn]',0.,180.,1.,nCh)
  !call CreateHistMC(dsigma,'Dilepton Cross Section dSigma/dRho [mubarn]',0.,0.2,0.002,nCh)
  ! set channel descriptions
  msigma%yDesc(1:4)   = (/ 'rho -> e+e-        ', 'omega -> e+e-      ', 'phi -> e+e-        ', 'omega -> pi0 e+e-  ' /)
  msigma%yDesc(5:8)   = (/ 'pi0 -> e+e- gamma  ', 'eta -> e+e- gamma  ', 'Delta -> N e+e-    ', 'eta -> e+e-        ' /)
  msigma%yDesc(9:12)  = (/ 'pn -> pn e+e-      ', 'pp -> pp e+e-      ', 'pi- n -> pi- n e+e-', 'pi- p -> pi- p e+e-' /)
  msigma%yDesc(13:15) = (/ 'pi+ n -> pi+ n e+e-', 'pi+ p -> pi+ p e+e-', 'Bethe-Heitler      ' /)
  call CopyDesc(psigma,msigma)
  call CopyDesc(esigma,msigma)
  call CopyDesc(ysigma,msigma)
  call CopyDesc(ptsigma,msigma)
  call CopyDesc(bgsigma,msigma)
  !call CopyDesc(oasigma,msigma)
  !call CopyDesc(dsigma,msigma)
  msigma%xDesc='Dilepton Mass [GeV]'
  psigma%xDesc='Momentum [GeV]'
  esigma%xDesc='Kinetic Energy [GeV]'
  ysigma%xDesc='Rapidity y'
  ptsigma%xDesc='Transverse Momentum [GeV]'
  bgsigma%xDesc='beta*gamma = p/m'
  !oasigma%xDesc='Opening Angle'
  !dsigma%xDesc='Density [fm^-3]'

  i = 0
  open(999,file=filename)
  do
    read(999,*,end=100) c,pout(1,0:3),pw,iso,f ! positron
    read(999,*,end=100) c,pout(2,0:3),pw,iso,f ! electron

    i = i + 1
    if (mod(i,10000)==0)   write (unit=*,fmt='(a)',advance='no') '.'
    if (mod(i,1000000)==0) write (unit=*,fmt='(a)')

    if (f == 0) cycle
    !!! === set up kinetic variables ===
    ptot=pout(1,:)+pout(2,:)                       ! total 4-momentum
    plab=sqrt(ptot(1)**2+ptot(2)**2+ptot(3)**2)    ! absolute 3-momentum
    mass=sqrt(ptot(0)**2-plab**2)                  ! invariant mass
    ekin=ptot(0)-mass                              ! kinetic energy
    y=0.5*log((ptot(0)+ptot(3))/(ptot(0)-ptot(3))) ! rapidity
    pt=sqrt(ptot(1)**2+ptot(2)**2)                 ! transverse momentum
    beta_gamma=plab/mass                           ! beta*gamma=p/m
    p1=pout(1,1:3)
    p2=pout(2,1:3)
    p1abs=sqrt(sum(p1*p1))
    p2abs=sqrt(sum(p2*p2))
    !angle=acos(dot_product(p1,p2)/sqrt(p1abs**2*p2abs**2))*180./pi ! opening angle in degrees
    !!! === put cross sections (averaged over beam energies) into histograms ===
    call AddHistMC(msigma,mass,iso,pw)
    call AddHistMC(psigma,plab,iso,pw)
    call AddHistMC(esigma,ekin,iso,pw)
    call AddHistMC(ysigma,y,iso,pw)
    call AddHistMC(ptsigma,pt,iso,pw)
    call AddHistMC(bgsigma,beta_gamma,iso,pw)
    !call AddHistMC(oasigma,angle,iso,pw)
    !call AddHistMC(dsigma,dens,iso,pw)
  end do
100  close(999)

   print *,""
   print *,"read ",i," dilepton pairs!"

  ! mass differential cross section
  call WriteHistMC(msigma,'DileptonMass.dat')
  !call WriteHistMC_Gauss(msigma,'DileptonMass_g.dat')
  ! momentum differential cross section
  call WriteHistMC(psigma,'DileptonPlab.dat')
  !call WriteHistMC_Gauss(psigma,'DileptonPlab_g.dat')
  ! kinetic energy differential cross section
  call WriteHistMC(esigma,'DileptonEkin.dat')
  !call WriteHistMC_Gauss(esigma,'DileptonEkin_g.dat')
  ! rapidity
  call WriteHistMC(ysigma,'DileptonY.dat')
  ! transverse momentum
  call WriteHistMC(ptsigma,'DileptonPt.dat')
  ! beta*gamma
  call WriteHistMC(bgsigma,'DileptonBetaGamma.dat')
  ! opening angle
  !call WriteHistMC(oasigma,'DileptonAngle.dat')
  ! density
  !call WriteHistMC(dsigma,'DileptonDensity.dat')

end program
