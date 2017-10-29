program testQE
use inputGeneral

call  readInputGeneral
call init_Database()
call test2
!call test_fg

end program testQE



subroutine test2
  use quasiElastic_electron, only : dSigmadcosTheta_l_dE_l_BW
  use particleDefinition
  use particleProperties
  use IDTable
  use energyCalc
  use output
  use vector
  use random
  implicit none

  real, parameter :: E_lepton_in=0.7
  real:: E_lepton_out,nuc_bareMass
  real,parameter:: theta_lepton_out=32.
  real, dimension(0:3) :: lf,pf,q
  type(particle) :: targetNuc
  integer :: i,j
  real, parameter :: dE=0.005
  integer :: maxSteps
  integer, dimension (1:2) :: num
  real :: absP, dummy,sigma_to,sigma_away,sigma
  real,parameter :: pfermi=0.2

  maxSteps=int(0.3/dE)

  targetNuc%charge=1
  targetNuc%position=0.
  targetNuc%ID=nucleon
  targetNuc%momentum(1:2)=0.
  targetNuc%mass=baryon(nucleon)%mass

  open(100,file="QE_toAway.dat")
  write(100,'("#",10A20)') 'E_gamma', 'sigma_away', 'sigma_to', 'sigma_away+sigma_to'

  do i=1,maxSteps
     num=0
     sigma_to=0.
     sigma_away=0.
     E_lepton_out=0.400+dE*float(i)
     write(*,*) E_lepton_out
     
     do j=-2000,2000,3
        call energyDetermination(targetNuc)
        ! (1) Find out q-direction by getting the cross section with some target input
        targetNuc%momentum(1:3)=0.
        call energyDetermination(targetNuc)
        dummy= dSigmadcosTheta_l_dE_l_BW(targetNuc,E_lepton_in,E_lepton_out,theta_lepton_out,lf,pf,nuc_bareMass,q)

        ! (2) Define target momentum in q-direction
        !absP=0.001*float(j)
        do 
           targetNuc%momentum(1)=(-1.+2.*rn())*pfermi
           targetNuc%momentum(2)=(-1.+2.*rn())*pfermi
           targetNuc%momentum(3)=(-1.+2.*rn())*pfermi
           if(dot_product(targetNuc%momentum(1:3),targetNuc%momentum(1:3)).lt.pfermi**2) exit
        end do

        call energyDetermination(targetNuc)

        sigma=  dSigmadcosTheta_l_dE_l_BW(targetNuc,E_lepton_in,E_lepton_out,theta_lepton_out,lf,pf,nuc_bareMass)
        if(absVec(pf(1:3)).gt.pfermi) then
           if(dot_Product(targetNuc%momentum(1:3),q(1:3)).gt.0) then
              sigma_away= sigma_away+ sigma
              num(1)=num(1)+1 
           else if(dot_Product(targetNuc%momentum(1:3),q(1:3)).lt.0) then
              sigma_to= sigma_to+ sigma
              num(2)=num(2)+1 
           end if
        end if
     end do
     write(100,'(10E20.5,L8)') E_lepton_in-E_lepton_out, sigma_away&
          & , sigma_to, sigma_away+sigma_to
  end do
end subroutine test2

subroutine test1
use quasiElastic_electron, only : dSigmadcosTheta_l_dE_l_BW
use particleDefinition
use particleProperties
use IDTable
use energyCalc
use output
use vector
implicit none

real, parameter :: E_lepton_in=0.7
real:: E_lepton_out,nuc_bareMass
real,parameter:: theta_lepton_out=32.
real, dimension(0:3) :: lf,pf,q
type(particle) :: targetNuc
integer :: i,j
real, parameter :: dE=0.0001
integer :: maxSteps
real :: absP, dummy

maxSteps=int(E_lepton_in/dE)

targetNuc%charge=1
targetNuc%position=0.
targetNuc%ID=nucleon
targetNuc%momentum(1:2)=0.
targetNuc%mass=baryon(nucleon)%mass



do j=-1,1

   absP=0.15*float(j)
   write(*,*) absP
   if(absP.lt.0) then
      open(100,file='QE_mom_-'//intToChar(nint(abs(absP*1000.)))//'_MeV.dat')
   else
      open(100,file='QE_mom_+'//intToChar(nint(abs(absP*1000.)))//'_MeV.dat')
   end if
   call energyDetermination(targetNuc)

   !targetNuc%momentum(0)=sqrt(baryon(nucleon)%mass**2 &
   !     & + Dot_Product(targetNuc%momentum(1:3),targetNuc%momentum(1:3)))+ &
   !     & + scalarPotential_nucleon(targetNuc%momentum(3),targetNuc%charge,targetNuc%position)

   write(100,'(A,4F9.4)')  ' # Momentum',absP
   write(100,'(A,I5)')    ' # Charge  ',targetNuc%charge
   write(100,'(A,3F9.4)') ' # Position',targetNuc%position
   do i=1,maxSteps
      E_lepton_out=dE*float(i)

      ! (1) Find out q-direction by getting the cross section with some target input
      targetNuc%momentum(1:3)=0.
      call energyDetermination(targetNuc)
      dummy= dSigmadcosTheta_l_dE_l_BW(targetNuc,E_lepton_in,E_lepton_out,theta_lepton_out,lf,pf,nuc_bareMass,q)
           
      ! (2) Define target momentum perpendicular to q-direction
      !targetNuc%momentum(1:3)=crossProduct(q(1:3)/absVec(q(1:3)),(/1.,0.,0./))*absP


      ! (2) Define target momentum in q-direction
      targetNuc%momentum(1:3)=q(1:3)/absVec(q(1:3))*absP
      call energyDetermination(targetNuc)
      write(100,'(10E20.5,L8)') E_lepton_in-E_lepton_out, &
           &  dSigmadcosTheta_l_dE_l_BW(targetNuc,E_lepton_in,E_lepton_out,theta_lepton_out,lf,pf,nuc_bareMass)&
           & , pf,lf,(absVec(pf(1:3)).lt.0.3)
      
   end do
end do
end subroutine test1


! Makes a FG integration
subroutine test_fg
use quasiElastic_electron, only : dSigmadcosTheta_l_dE_l_BW
use particleDefinition
use particleProperties
use IDTable
use random
implicit none

real, parameter :: E_lepton_in=0.7
real:: E_lepton_out, sigma, sigmaTot,absP,nuc_bareMass
real,parameter:: theta_lepton_out=32.
real, dimension(0:3) :: lf,pf
type(particle) :: targetNuc
integer :: i,j
real, parameter :: dE=0.005
integer :: maxSteps
real, parameter :: pfermi=0.205
integer :: numPoints=5000

maxSteps=int((E_lepton_in-0.3)/dE)

targetNuc%charge=1
targetNuc%position=0.
targetNuc%ID=nucleon
targetNuc%momentum(1:2)=0.
targetNuc%mass=baryon(nucleon)%mass

sigmaTot=0.


do i=1,maxSteps

   E_lepton_out=0.3+dE*float(i)

   sigma=0.

   do j=1,numPoints
      do
         targetNuc%momentum(1)=(1.-2.*rn())*pfermi
         targetNuc%momentum(2)=(1.-2.*rn())*pfermi
         targetNuc%momentum(3)=(1.-2.*rn())*pfermi
         absP=sqrt(Dot_product(targetNuc%momentum(1:3),targetNuc%momentum(1:3)))
         if (absP.lt.pfermi) exit
      end do


      targetNuc%momentum(0)=sqrt(baryon(nucleon)%mass**2 + absP**2)+ &
           & + scapot(1,targetNuc%charge,(/0.,0.,p/),targetNuc%position)

      sigma=sigma+ dSigmadcosTheta_l_dE_l_BW(targetNuc,E_lepton_in,E_lepton_out,theta_lepton_out,lf,pf,nuc_bareMass)
   end do
   sigma=sigma/float(numPoints)
   sigmaTot=sigmaTot+sigma*dE
   write(11,*) E_lepton_in-E_lepton_out, sigma
end do
   write(*,*) 'total Xsection=',sigmaTot

 end subroutine test_fg

