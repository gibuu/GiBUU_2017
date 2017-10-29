!*******************************************************************************************************
!****p* test_elec/
! NAME
! program test_elec
! PURPOSE
! * Basic test of hadron tensor and the vacuum cross section for e^- n -> e^- n pi
! * Results have to be compared to MAID results for virtual photon cross sections.
! * See http://www.kph.uni-mainz.de/MAID/maid2003
! NOTES
! * There are three namelists in file "job". One of them is the main one:
!   - NAMELIST /testing/ switch  ! Main namelist. Specifies kind of test.
! * If SWITCH:
!
!      =1 : Calculate dSigma/dOmega_f/dE_f/dOmega_k in the CM-Frame of the hadronic vertex. Output to fort.11.
!           Namelist Xsection must be given
!
!      =2 : Evaluate virtual photon Xsections in the CM-Frame of the hadronic vertex as function of W at fixed Q^2.
!    Output to fort.77  QSquared in namelist virtual must be given as input
!
!
!      =3 : Evaluate virtual photon Xsections in the CM-Frame of the hadronic vertex as function of theta and phi.
!    Output to fort.300+phi_k.  Namelist virtual must be given
!
!      =4 : Calculate dSigma/dOmega_f/dE_f/dOmega_k in the CM-Frame of the hadronic vertex. Only for phi_k=30°.
!           Namelist Xsection must be given
!
!      =5 : Calculate dSigma/dOmega_f/dE_f/dOmega_k in the lab-Frame. Output to fort.11.
!           Namelist Xsection must be given
!
!      =6 : Calculate dSigma/dOmega_f/dE_f by integrating over dOmega_k in the lab-Frame.
!           Namelist Xsection must be given, however theta_k must not be set
!************************************************************************************

program test_elec


  use electronPionProduction_xSection
  use degRad_conversion, only : degrees, radian
  use inputGeneral
  implicit none
  real, parameter:: me=0.000510
  real, parameter:: mn=0.938
  real :: dummy
  real :: phi_k, theta_k,theta_lf
  real:: energy_li, energy_lf
  real :: qSquared, W,thetaCM,sigmaV
  integer :: switch,j,i
  real , dimension (0:3) :: q
  real :: sigmaT,sigmaTT, sigmaL,sigmaTL,  sigmaVirtual,eps,gamma
  real :: sigma,cost
  integer  :: numCost
  integer , parameter :: numPhi=1000
  real :: dCost,dPhi
  integer :: charge_pionOut,charge_nucOut

  NAMELIST /testing/ switch
  NAMELIST /Xsection/ energy_li,energy_lf,theta_lf,theta_k,charge_pionOut,charge_nucOut
  NAMELIST /virtual/ W, Qsquared,charge_pionOut,charge_nucOut


  call readinputGeneral
  call init_database


  rewind(5)
  read(5,nml=testing)
  write(*,*) 'switch=',switch


  ! Reading out jobcards:
  select case(switch)
  case(1,4,5,6)
     rewind(5)
     read(5,nml=Xsection)

     write(*,*) 'Outgoing pion     charge:',charge_pionOut
     write(*,*) 'Outgoing nucleon  charge:',charge_nucOut
     write(*,*)
     write(*,*) 'Energy incoming:',energy_li
     write(*,*) 'Energy outgoing:',energy_lf
     write(*,*) 'Electron scattering angle (degree):', theta_lf
     write(*,*) 'Pion scattering angle theta (degree) :', theta_k
     write(*,*)
  case(2,3)
     rewind(5)
     read(5,nml=virtual)
     write(*,*) 'Outgoing pion     charge:',charge_pionOut
     write(*,*) 'Outgoing nucleon  charge:',charge_nucOut
     write(*,*)
     write(*,*) 'W:', W
     write(*,*) 'QSquared:', Qsquared
   case default
      write(*,*) 'This choice for switch is not valid'
      stop
  end select

  ! Performing job

  if(switch.eq.1) then
     write(11,*) '# phi_k ,theta_k,theta_k (CM) ,dsigma/d(...)'
     do phi_k=0,360,30
        dummy=dSigmadOmega_fdE_fdOmega_k(charge_pionOut,charge_nucOut,&
             & me,me,mn,mn,energy_li,energy_lf,theta_lf,phi_k,theta_k,&
             & .true.,thetaCM)
        sigmaV=dSigmadOmega_fdE_fdOmega_k_withVirtual(charge_pionOut,charge_nucOut,&
             & me,me,mn,mn,energy_li,energy_lf,theta_lf,phi_k,theta_k,&
             & .true.,thetaCM)
        write(*,'(5E16.5)')  phi_k,theta_k,thetaCM,dummy,sigmaV
        write(11,'(5E16.5)') phi_k ,theta_k,thetaCM,dummy,sigmaV
     end do


  else if(switch.eq.2) then
     write(*,*) '#####################################################'
     w=1.1
     write(77,'(A)') '#W, Qsquared,sigmaT,sigmaL, sigmaTT,sigmaTL'
     do i=0,15
        w=w+0.05
        call integrate_virtualPhotonXsection(charge_pionOut,charge_nucOut,sigmaT,sigmaL,sigmaTL, sigmaTT,W,QSquared)
        write(77,'(6E22.10)') W, Qsquared,sigmaT,sigmaL, sigmaTT,sigmaTL
        write(*,'(6E22.10)') W, Qsquared,sigmaT,sigmaL, sigmaTT,sigmaTL
     end do


  else if(switch.eq.3) then
     write(*,*) '#####################################################'
     do phi_k=0,360,30
        write(300+phi_k,'(A,F10.3)')'# W:', W
        write(300+phi_k,'(A,F10.3)')'# QSquared:', Qsquared
        write(300+phi_k,'(A,F10.3)')'# phi=',phi_k
        write(300+phi_k,'(A)')'# theta(pion),phi(pion), sigmaT,sigmaL, sigmaTT,sigmaTL (mcb/sr)'
        do theta_k=0,180,10
           call virtualPhotonXsection(charge_pionOut,charge_nucOut,sigmaT,sigmaL,sigmaTL, sigmaTT,W,QSquared,theta_k,phi_k,.true.)
           write(*,*)
           write(*,*)  'sigma_T=', sigmaT
           write(*,*)  'sigma_L=', sigmaL
           write(*,*)  'sigma_TL=', sigmaTL
           write(*,*)  'sigma_TT=', sigmaTT
           write(300+phi_k,'(6E22.10)')theta_k,phi_k, sigmaT,sigmaL, sigmaTT,sigmaTL
        end do
     end do

  else if(switch.eq.4) then
     rewind(5)
     read(5,nml=Xsection)

     write(*,*)
     write(*,*) 'Energy incoming:',energy_li
     write(*,*) 'Energy outgoing:',energy_lf
     write(*,*) 'Electron scattering angle (degree):', theta_lf
     write(*,*) 'Pion scattering angle theta (degree) :', theta_k
     write(*,*)
     dummy=dSigmadOmega_fdE_fdOmega_k(charge_pionOut,charge_nucOut,&
          & me,me,mn,mn,energy_li,energy_lf,theta_lf,30.,theta_k,&
          & .true.,thetaCM)
     write(*,*)  30.,thetaCM, dummy

  else if(switch.eq.5) then
     rewind(5)
     read(5,nml=Xsection)

     write(*,*)
     write(*,*) 'Energy incoming:',energy_li
     write(*,*) 'Energy outgoing:',energy_lf
     write(*,*) 'Electron scattering angle (degree):', theta_lf
     write(*,*) 'Pion scattering angle theta (degree) :', theta_k
     write(*,*)
     write(11,*) '# phi_k ,theta_k,theta_k (CM) ,dsigma/d(...)'
     do phi_k=0,360,30
        dummy=dSigmadOmega_fdE_fdOmega_k(charge_pionOut,charge_nucOut,&
             & me,me,mn,mn,energy_li,energy_lf,theta_lf,phi_k,theta_k,.false.)
        write(*,'(5E16.5)')  phi_k,theta_k,dummy
        write(11,'(5E16.5)') phi_k ,theta_k,dummy
     end do
  else if(switch.eq.6) then
     rewind(5)
     read(5,nml=Xsection)

     write(*,*)
     write(*,*) 'Energy incoming:',energy_li
     write(*,*) 'Energy outgoing:',energy_lf
     write(*,*) 'Electron scattering angle (degree):', theta_lf
     write(*,*) 'Pion scattering angle theta (degree) :', theta_k
     write(*,*)
     write(11,*) '# phi_k ,theta_k,theta_k (CM) ,dsigma/d(...)'

     numcost=1
     do
        if(numcost.gt.10000) exit
        write(*,*)
        write(*,'(A,I7)') 'number of integration points in cos(theta):', numcost
        write(*,'(A,I7)') 'number of integration points in phi       :', numphi
        write(33,'(A,I7)') 'number of integration points in cos(theta):', numcost
        write(33,'(A,I7)') 'number of integration points in phi       :', numphi

        dcost=2./float(numCost)
        dphi=360/float(numPhi)

        sigma=0.
        do j=0,numCost
           cost=-1+j*dcost
           theta_k=degrees(acos(cost))
           do i=0,numPhi
              phi_k=dPhi*i
              sigma=sigma+dSigmadOmega_fdE_fdOmega_k(charge_pionOut,&
                   & charge_nucOut,me,me,mn,mn,energy_li,energy_lf,theta_lf,&
                   & phi_k,theta_k,.false.)*dcost*radian(dphi)
           end do
        end do
        write(*,*) 'sigma=', sigma
        write(33,*) 'sigma=', sigma
        numcost=numcost*5
     end do

  end if



end program test_elec


!*******************************************************************************************************
!****s* virtualPhotonXsection/
! NAME
! subroutine  virtualPhotonXsection(charge_pionOut,charge_nucOut,sigmaT,sigmaL,sigmaTL, sigmaTT,W,QSquared,theta_k,phi_k)
! PURPOSE
! * Basic test of hadron tensor. Results have to be compared to MAID results for virtual photon cross sections.
! * See http://www.kph.uni-mainz.de/MAID/maid2003/observ.html
! * Evaluates the virtual photon nucleon -> pion nucleon cross sections dsigma_T, dsigma_L, dsigma_TT,dsigma_TL .
! * All in units of  (mcb/sr) .
! * These cross sections are calculated in the CM-Frame.
! INPUTS
! * W      ! sqrt(s)
! * QSquared ! -q^2
! * real,intent(in) :: phi_k, theta_k             ! pion scattering angles in units of degree
! * All angles in degree.
! * integer, intent(in) :: charge_pionOut,charge_nucOut ! Outgoign nucleon and pion charge
! RESULT
! * real :: sigmaT,sigmaL,sigmaTL, sigmaTT,
! * sigmaVirtual has not been considered (there is a problem since epsilon is evaluated in CM instead of lab frame
!****************************************************************************************************
subroutine  virtualPhotonXsection(charge_pionOut,charge_nucOut,sigmaT,sigmaL,sigmaTL, sigmaTT,W,QSquared,theta_k,phi_k,print)


  use twoBodyTools, only : pcm
  use idTable, only : nucleon, pion
  use particleProperties, only : baryon, meson
  use minkowski, only : pair => SP
  use degRad_conversion, only : degrees, radian
  use hadronTensor_npi
  use formFactors_A_main, only : getA
  implicit none

  real, intent(out) :: sigmaT, sigmaL, sigmaTL, sigmaTT
  real, intent(in)  :: W ! =sqrt(s)
  real, intent(in)  :: QSquared, theta_k,phi_k
  logical, intent(in) :: print
  real, dimension(1:3) :: k_unit
  real :: p_CM,kOverkGamma,Normalization
  real, dimension(0:3) :: q,pin,pout,k
  complex, dimension(1:6) :: A
  integer, intent(in) :: charge_pionOut,charge_nucOut ! Outgoign nucleon and pion charge
  ! Define unit vector in direction of outgoing pion:
  real :: dummy,mn
  real :: sol1,sol2
  if(print) then
     write(*,*) 'W', W
     write(*,*) 'QSquared', QSquared
  end if
  mn=baryon(nucleon)%mass

  p_CM=pcm(W,baryon(nucleon)%mass,meson(pion)%mass)

  ! Outgoing particles in CM-Frame:
  k(1:3)=p_CM*(/sin(radian(theta_k))*cos(radian(phi_k)),sin(radian(theta_k))*sin(radian(phi_k)),cos(radian(theta_k))/)
  k(0)=sqrt(meson(pion)%mass**2+p_cm**2)
  pout(1:3)=-k(1:3)
  pout(0)=sqrt(baryon(nucleon)%mass**2+p_cm**2)
  if(print) then
     write(*,*) 'W', W
     write(*,*) 'QSquared', QSquared
  end if

  dummy=(W**2+Qsquared-mn**2)/2.

  ! Incoming particles in CM-Frame:

  q(3)=sqrt((dummy**2+QSquared*mn**2)/(2.*dummy+mn**2-Qsquared))
  if(print) write(*,*) 'q3',q(3)
  q(3)=sqrt((QSquared+(w+mn)**2)*(QSquared+(w-mn)**2))/2/W
  if(print) write(*,*) 'q3',q(3)

  q(0)=sqrt(-QSquared+q(3)**2)
  pin(3)=-q(3)
  pin(1:2)=0
  pin(0)=sqrt(baryon(nucleon)%mass**2+pin(3)**2)
  if(abs(sqrt(pair(pin+q,pin+q))-W).gt.0.001) then
     q(0)=-sqrt(-QSquared+q(3)**2)
  end if

  if(print) then
  write(*,*)
  write(*,'(A)')'********** Incoming :'
  write(*,'(A,4F9.5)') 'pi=',pin
  write(*,'(A,4F9.5)') 'q=',q
  write(*,'(A,4F9.5)') 'Total=' , q+pin
  write(*,*)
  write(*,'(A)')'***********Outgoing :'
  write(*,'(A,4F9.5)') 'pf=',pout
  write(*,'(A,4F9.5)') 'k=',k
  write(*,'(A,4F9.5)') 'Total=' , k+pout
  write(*,*)
  write(*,'(A,4F9.5)') 's        =' , pair(pin+q,pin+q)
  write(*,'(A,4F9.5)') 't        =' , pair(k-q,k-q)
  write(*,*)
  write(*,'(A,4F9.5)') 'W=sqrt(s)=' , sqrt(pair(pin+q,pin+q))
  write(*,'(A,4F9.5)') 'Q^2      =' , -pair(q,q)
  end if


  A=getA(charge_pionOut,charge_nucOut,theta_k,w**2,QSquared)

  !H_munu(mu,nu,pin,pout,k,q,B)

  kOverkGamma=sqrt(dot_Product(k(1:3),k(1:3)))/ &
       & ((w**2-baryon(nucleon)%mass**2)/2./w)

  Normalization=(baryon(nucleon)%mass/(4.*3.14*W))**2


  if(print) write(*,*) kOVerkGamma,Normalization

  sigmaT=(H_munu(1,1,pin,pout,k,q,A)+H_munu(2,2,pin,pout,k,q,A))/2.*Normalization* kOverkGamma

  ! sigmaT is in GEV**-2=0.197fm**2=10*0.197**2 mb=10000 * 0.197**2 mcb
  sigmaT=sigmaT*0.197**2*10000.

  sigmaL=H_munu(3,3,pin,pout,k,q,A)*kOverKGamma*Normalization
  sigmaL=sigmaL*0.197**2*10000./q(0)**2*qSquared


  sigmaTL=-1./cos(radian(phi_k))     *REAl(H_munu(1,3,pin,pout,k,q,A))*kOverKGamma*Normalization
  sigmaTL=sigmaTL*0.197**2*10000./sqrt(q(0)**2/qSquared)


  sigmaTT=1./cos(2*radian(phi_k))   *(H_munu(1,1,pin,pout,k,q,A)-H_munu(2,2,pin,pout,k,q,A))/2.*kOverKGamma*Normalization
  sigmaTT=sigmaTT*0.197**2*10000.




end subroutine virtualPhotonXsection



!*******************************************************************************************************
!****s* integrate_virtualPhotonXsection
! NAME
! subroutine  integrate_virtualPhotonXsection(charge_pionOut,charge_nucOut,sigmaT,sigmaL,sigmaTL, sigmaTT,W,QSquared)
! PURPOSE
! * Basic test of hadron tensor. Results have to be compared to MAID results for virtual photon cross sections.
! * http://www.kph.uni-mainz.de/MAID/maid2007/total.html
! * Evaluates the virtual photon nucleon -> pion nucleon cross sections sigma_T, sigma_L, sigma_TT,sigma_TL .
! * All in units of  (mcb) .
! * These cross sections are calculated in the CM-Frame !!!
! INPUTS
! * W      ! sqrt(s)
! * QSquared ! -q^2
! * integer, intent(in) :: charge_pionOut,charge_nucOut ! Outgoign nucleon and pion charge
! RESULT
! * real :: sigmaT,sigmaL,sigmaTL, sigmaTT integrate over pion angle
!*******************************************************************************************************

subroutine  integrate_virtualPhotonXsection(charge_pionOut,charge_nucOut,sigmaT,sigmaL,sigmaTL, sigmaTT,W,QSquared)
  use constants, only : pi
  use degrad_conversion
  implicit none
  real, intent(out) :: sigmaT, sigmaL, sigmaTL, sigmaTT
  real, intent(in)  :: W ! =sqrt(s)
  real, intent(in)  :: QSquared
  integer :: charge_pionOut,charge_nucOut
  real ::  theta_k,phi_k
  real :: sigmaT_tot, sigmaL_tot, sigmaTL_tot, sigmaTT_tot
  integer :: i,j
  real :: dphi,dcostheta,dOmega
  integer:: numsteps_phi=100
  integer:: numsteps_theta=1000

  dphi=2*pi/float(numsteps_phi)
  dcosTheta=2./float(numsteps_theta)
  sigmaT_tot=0.
  sigmaL_tot=0.
  sigmaTL_tot=0.
  sigmaTT_tot=0.
  do i=1,numsteps_theta
     do j=1,numsteps_phi
        theta_k=degrees(acos(-1.+(float(i)-0.5)*dcosTheta))
        phi_k=degrees((float(j)-0.5)*dPhi)
        call  virtualPhotonXsection(charge_pionOut,charge_nucOut,sigmaT,sigmaL,sigmaTL, sigmaTT,W,QSquared,theta_k,phi_k,.false.)
        sigmaT_tot=sigmaT_tot+sigmaT
        sigmaL_tot=sigmaL_tot+sigmaL
        sigmaTL_tot=sigmaTL_tot+sigmaTL
        sigmaTT_tot=sigmaTT_tot+sigmaTT
     end do
  end do
  dOmega=dcosTheta*dphi
  sigmaT=sigmaT_tot*dOmega
  sigmaL=sigmaL_tot*dOmega
  sigmaTL=sigmaTL_tot*dOmega
  sigmaTT=sigmaTT_tot*dOmega



end subroutine integrate_virtualPhotonXsection
