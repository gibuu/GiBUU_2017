

!*******************************************************************************************************
!****p* test_med/
! NAME
! program test_med
! PURPOSE
! * Basic test of the in-medium X section for e^- n -> e^- n pi
! * Results have to be compared to MAID results for virtual photon cross sections in the vacuum
! * See http://www.kph.uni-mainz.de/MAID/maid2003
!************************************************************************************


program test_med
  use inputGeneral, only: readinputGeneral
  implicit none
  call init_database
  call readinputGeneral
 
!  call testKine
  call testCross

end program test_med


!*******************************************************************************************************
!****s* testCross/
!
! PURPOSE
! Calculates dSigma/dOmega_f/dE_f/dOmega_k in the nucleon rest-frame. 
!
! INPUTS
! * Namelist Xsection : energy_li,energy_lf,theta_lf,theta_k,charge_pionOut,charge_nucOut
!
! OUTPUT
! * File fort.111
!************************************************************************************
subroutine testCross

  use electronPionProduction_medium
  use particleDefinition
  use idTable
  use energyCalc
  use output
  use propagation, only: updateVelocity
  use particleProperties
  use IdTable
  use pauliBlockingModule, only : pauliBlocking

  use degRad_conversion, only : degrees, radian 

  implicit none
  real, parameter:: me=0.000510
  real, parameter:: mn=0.938
  real :: dummy
  real :: phi_k, theta_k,theta_lf
  real:: energy_li, energy_lf
  real :: qSquared, W,thetaCM,sigmaV
  integer :: pionCharge,i,j
  real , dimension (0:3) :: q
  real ,dimension(0:3) :: lf,k,pf
  type(particle) :: initNuc
  real :: sigma,cost
  integer  :: numCost
  integer , parameter :: numPhi=1000
  real :: dCost,dPhi
  integer :: charge_pionOut,charge_nucOut
  NAMELIST /Xsection/ energy_li,energy_lf,theta_lf,theta_k,charge_pionOut,charge_nucOut


  rewind(5)
  read(5,nml=Xsection)
  write(*,*)   
  write(*,*) 'Outgoing pion     charge:',charge_pionOut
  write(*,*) 'Outgoing nucleon  charge:',charge_nucOut
  write(*,*)
  write(*,*) 'Energy incoming:',energy_li
  write(*,*) 'Energy outgoing:',energy_lf
  write(*,*) 'Electron scattering angle (degree):', theta_lf
  write(*,*) 'Pion scattering angle theta (degree) :', theta_k
  write(*,*)
  write(111,*) '# phi_k ,theta_k,dsigma/d(...)'

  initNuc%position=(/1.,2.,3./)
  initNuc%mass=0.938
  initNuc%momentum(1:3)=(/0.,0.,0./)
  initNuc%ID=nucleon
  initNuc%charge=charge_nucOut+charge_pionOut
  initNuc%perturbative=.false.
  call energyDetermination(initNuc)
  call updateVelocity     (initNuc)
  write(*,*)
  write(*,*)
  write(*,*) 'pi(0)', initNuc%momentum(0)
  write(*,*)
  write(*,*)
  write(*,*)
  write(*,*)
  write(*,*) 'velocity', initNuc%velocity
  write(*,*)
  write(*,*)

  write(111,'(A)') 'phi_k ,theta_k, dSigma/dOmega_f/dE_f/dOmega_k in the nucleon rest-frame'
  do phi_k=0,360,10
     dummy=dSigmadOmega_fdE_fdOmega_k_med(initNuc,charge_pionOut,&
          & energy_li,energy_lf,theta_lf,phi_k,theta_k,q,lf,k,pf)
!     if(pauliBlocking(pf,initNuc%position,initNuc%charge-pionCharge) ) dummy=0.
     write(*,'(5E16.5)')   phi_k,theta_k,dummy
     write(111,'(5E16.5)') phi_k ,theta_k,dummy
  end do

  
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
           sigma=sigma+dSigmadOmega_fdE_fdOmega_k_med(initNuc,charge_pionOut,energy_li,&
                & energy_lf,theta_lf,phi_k,theta_k,q,lf,k,pf)*dcost*radian(dphi)
        end do
     end do
     write(*,*) 'sigma=', sigma
     write(33,*) 'sigma=', sigma
     numcost=numcost*5
  end do



end subroutine testCross




subroutine testKine
  use electronPionProduction_medium
  use particleDefinition
  use idTable
  use energyCalc
  use output
  use propagation, only: updateVelocity

  implicit none
  integer :: i,ios
  real, dimension(1:3) :: position
  real, dimension(0:3) :: pin,pf,li,lf,k,q
  real,save :: energy_lepton_in=1.2
  real,save :: energy_lepton_out=0.8
  real,save :: theta_lepton=40
  real,save :: theta_pion=30
  real,save :: phi_pion=0
  type(particle) :: initNuc
  integer :: pionCharge
  logical :: success
  integer :: charge_pionOut,charge_nucOut

  NAMELIST /init_lowElectron/ pionCharge, energy_lepton_in, energy_lepton_out,theta_lepton,phi_pion, theta_pion

  character(*), parameter :: format1 = '("   Simulation type is ",A)'

  call Write_ReadingInput("init_lowElectron",0)
  rewind(5)
  read(5,nml=init_lowElectron,IOSTAT=ios)
  call Write_ReadingInput("init_lowElectron",0,ios)

  write(11,*) 'pionCharge, energy_lepton_in, energy_lepton_out,theta_lepton,phi_pion, theta_pion'
  write(11,*) pionCharge, energy_lepton_in, energy_lepton_out,theta_lepton,phi_pion, theta_pion

  call Write_ReadingInput("init_lowElectron",1)


  initNuc%position=(/1.,2.,3./)
  initNuc%mass=0.938
  initNuc%momentum(1:3)=(/0.2,0.1,0.05/)
  initNuc%ID=nucleon
  initNuc%charge=1
  initNuc%perturbative=.false.
  call energyDetermination(initNuc)
  call updateVelocity     (initNuc)
  write(*,*)
  write(*,*)
  write(*,*) 'pi(0)', initNuc%momentum(0)
  write(*,*)
  write(*,*)
  write(*,*)
  write(*,*)
  write(*,*) 'velocity', initNuc%velocity
  write(*,*)
  write(*,*)


  write(*,*) 'call getKinematics'
  !  subroutine getKinematics(pionCharge,initNuc,energy_li,energy_lf,theta_lf,phi_k,theta_k,li,lf,q,k,pf)

  call getKinematics(pionCharge,initNuc,energy_lepton_in,energy_lepton_out,theta_lepton,phi_pion,theta_pion,li,lf,q,k,pf,success)

  if(success) then

     write(*,'(4E12.4)') li
     write(*,'(4E12.4)') lf
     write(*,*)
     write(*,'(4E12.4)') q
     write(*,'(4E12.4)') lf-li
     write(*,'(4E12.4)') pf+k-initNuc%momentum
     write(*,*)
     write(*,'(4E12.4)') k
     write(*,'(4E12.4)') pf
     write(*,*)
     write(*,'(4E12.4)') li+initNuc%momentum
     write(*,'(4E12.4)') lf+k+pf


     write(*,*)
     write(*,*)
     write(*,*) 'pf(0)', pf(0)
     write(*,*)
     write(*,*)

     initNuc%position=(/1.,2.,3./)
     initNuc%mass=0.938
     initNuc%momentum(0:3)=pf
     initNuc%ID=nucleon
     initNuc%charge=1
     initNuc%perturbative=.false.
     call energyDetermination(initNuc)
     write(*,*) 'pf(0)', initNuc%momentum(0)
  else
     write(*,*) 'Failure in establishing kinematics'
  end if


 write(*,*)  dSigmadOmega_fdE_fdOmega_k_med(initNuc,pionCharge,energy_lepton_in,energy_lepton_out,theta_lepton,phi_pion,theta_pion,q,lf,k,pf)
  if(success) then

     write(*,'(4E12.4)') li
     write(*,'(4E12.4)') lf
     write(*,*)
     write(*,'(4E12.4)') q
     write(*,'(4E12.4)') lf-li
     write(*,'(4E12.4)') pf+k-initNuc%momentum
     write(*,*)
     write(*,'(4E12.4)') k
     write(*,'(4E12.4)') pf
     write(*,*)
     write(*,'(4E12.4)') li+initNuc%momentum
     write(*,'(4E12.4)') lf+k+pf


     write(*,*)
     write(*,*)
     write(*,*) 'pf(0)', pf(0)
     write(*,*)
     write(*,*)

     initNuc%position=(/1.,2.,3./)
     initNuc%mass=0.938
     initNuc%momentum(0:3)=pf
     initNuc%ID=nucleon
     initNuc%charge=1
     initNuc%perturbative=.false.
     call energyDetermination(initNuc)
     write(*,*) 'pf(0)', initNuc%momentum(0)
  else
     write(*,*) 'Failure in establishing kinematics'
  end if


end subroutine testKine




subroutine testKai
  ! USED to find "abs()"-problem in electronPionProduction_medium in revision 2184
  use electronPionProduction_medium
  use particleDefinition
  use idTable
  use energyCalc
  use output
  use propagation, only: updateVelocity
  use particleProperties
  use IdTable
  use pauliBlockingModule, only : pauliBlocking

  use degRad_conversion, only : degrees, radian 

  implicit none
  real, parameter:: me=0.000510
  real, parameter:: mn=0.938
  real :: dummy
  real :: phi_k, theta_k,theta_lf
  real:: energy_li, energy_lf
  real :: qSquared, W,thetaCM,sigmaV
  integer :: pionCharge,i,j
  real , dimension (0:3) :: q
  real ,dimension(0:3) :: lf,k,pf
  type(particle) :: initNuc
  real :: sigma,cost
  integer  :: numCost
  integer , parameter :: numPhi=1000
  real :: dCost,dPhi
  integer :: charge_pionOut,charge_nucOut
  real ::phi_k_center,theta_k_center



!!$Input:
!!$
!!$Initial nucleon   :
!!$ * Momentum :  0.93800  0.00000  0.00000  0.00000
!!$ * Position :  1000.00000  1000.00000  1000.00000
!!$ * Charge   :        1
!!$ * Mass     :  0.93800
!!$ * ID       :       1
!!$ * perturbative      :       F
!!$
!!$Initial electron   :
!!$* electron energy incoming :  8.97157
!!$
!!$Final electron   :
!!$* electron energy outgoing :  7.07383
!!$* Theta of electron : 12.48190
!!$
!!$Final pion   :
!!$Phi of pion       :256.19640
!!$Theta of pion     : 11.64209


  initNuc%position=(/1000.,1000.,1000./)
  initNuc%mass=0.938
  initNuc%momentum(1:3)=(/0.,0.,0./)
  initNuc%ID=nucleon
  initNuc%charge=1
  initNuc%perturbative=.false.
  call energyDetermination(initNuc)
  call updateVelocity     (initNuc)
  write(*,*)
  write(*,*)
  write(*,*) 'pi(0)', initNuc%momentum(0)
  write(*,*)
  write(*,*)
  write(*,*)
  write(*,*)
  write(*,*) 'velocity', initNuc%velocity
  write(*,*)
  write(*,*)
  energy_lf=7.07383
  energy_li=8.97157
  theta_lf= 12.48190
  phi_k=256.19640
  theta_k=11.64209
  charge_pionOut=0


  write(111,'(A)') '#phi_k ,theta_k, dSigma/dOmega_f/dE_f/dOmega_k in the nucleon rest-frame'
  do phi_k=260,260,10
     dummy=dSigmadOmega_fdE_fdOmega_k_med(initNuc,charge_pionOut,&
          & energy_li,energy_lf,theta_lf,phi_k,theta_k,q,lf,k,pf)
!     if(pauliBlocking(pf,initNuc%position,initNuc%charge-pionCharge) ) dummy=0.
     write(*,'(5E16.5)')   phi_k,theta_k,dummy
     write(111,'(5E16.5)')   phi_k,theta_k,dummy
     !write(111,'(5E16.5)') phi_k,theta_k,phi_k_center,theta_k_center,dummy
  end do

  close(111)
end subroutine testKai
