program test
  ! Tests the lowPhototon analysis: in particular kruscheAnalyse


  use particleProperties
  use particleDefinition
  use initLowPhoton, only  : getEnergyGamma
  use LowPhotonAnalysis, only : kruscheAnalyse
  use vector
  use minkowski
  use inputGeneral
  use idtable
  use random
  use lorentzTrafo
  use twoBodyTools


  implicit none

  type(particle), dimension(1:100000,1:2) :: part
  type(particle) :: t
  real, dimension(0:3) :: p, ptot
  real, dimension(1:3) :: beta
  real, parameter :: pfermi=0.001
  real :: s, egamma
  integer :: index,i
  logical, parameter :: debug=.false.
  

  call init_database
  call readinputGeneral  
  egamma= getEnergyGamma()

  write(*,*) lbound(part,dim=1),ubound(part,dim=1)

  ! Initialize pions and nucleons isotropic in CM-Frame and boost them to LAB frame
  do index=lbound(part,dim=1),ubound(part,dim=1)
     pLoop: do
        do i=1,3
           p(i)=pfermi*(1.-2.*rn())
        end do
        if(absVec(p(1:3)).lt.pfermi) exit 
     end do pLoop
     p(0)=sqrt(baryon(nucleon)%mass**2+Dot_product(p(1:3),p(1:3)))
     ptot(0:3)=p(0:3)+(/Egamma,0.,0.,Egamma/)
     s=sp(pTot,ptot)
     beta=ptot(1:3)/ptot(0)
     part(index,1)%ID=pion
     part(index,1)%charge=0
     part(index,1)%mass=meson(pion)%mass
     part(index,2)%ID=nucleon
     part(index,2)%charge=0
     part(index,2)%mass=baryon(nucleon)%mass

     part(index,1)%momentum(1:3)=   pCM(sqrt(s), meson(pion)%mass, baryon(nucleon)%mass)    *rnOmega()
     part(index,1)%momentum(0)  =freeEnergy(part(index,1))


     part(index,1:2)%perweight=1.

     part(index,2)%momentum(1:3)=   -part(index,1)%momentum(1:3)
     part(index,2)%momentum(0)  =freeEnergy(part(index,2))

     
     call lorentz(-beta, part(index,1)%momentum)
     call lorentz(-beta, part(index,2)%momentum)
     if(debug) then
        write(*,'(4E14.4)') part(index,1)%momentum+part(index,2)%momentum
        write(*,'(4E14.4)') ptot
        write(*,*) s, sqrtS(part(index,1),part(index,2))**2
     end if

  end do

  ! Evaluate deltaE
  call kruscheAnalyse(part,.true.)
end program test
