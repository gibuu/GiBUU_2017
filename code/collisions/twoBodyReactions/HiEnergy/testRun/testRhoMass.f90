program testRhoMass
use PythiaSpecFunc!, only: srtFreeVMMass
use histf90
use mesonWidth, only : FullWidthMeson
use particleProperties
use idTable, only : rho
use output
implicit none
integer :: i,j
type(histogram) :: massDistri
real :: mass,pcm,mres0, gamtot,s,spectral, norm
character(3) :: zahl
!real,external:: VM_Mass

call init_Database

Do i=0,350
   zahl=intToChar(i)
   Open(101,File='Simulated.'//zahl) ! Contains results of "rhoMass"
   Open(501,File='Real.'//zahl)      ! Contains original distribution function

   !************************************************************************************************************
   ! Plot histogram of distribution of rhoMass
   s = i*0.05
   call Init_VM_Mass(s)
   call CreateHist(massDistri,"mass Distribution",0.,30.,0.01)
   do j=0,100000
      mass=VM_Mass(113)
      call AddHist(massDistri, mass,1.)
   end do
   write(101,*) '#### Simulated Distribution function for dsigma/dmass(rhoMeson)'
   write(101,*) '#### srts=', s
!   call writeHist(massDistri,i+101,0.,1./massdistri%yval(0,1))
   call writeHist(massDistri,101,0.,1./10000.)
   call RemoveHist(massDistri)

   !************************************************************************************************************
   ! Plot original distribution function, which ought to have same shape as the results from above, just scaled...
   ! Normalized to one at the rho pole mass

   ! * Evaluate normalization
   mass=meson(rho)%mass
   gamtot=FullWidthMeson(rho,mass)
   pcm=0.25*(s**2-mass**2+0.938**2)**2/s-0.938**2
   norm=mass**2*gamtot/  (gamtot**2*mass**2)*sqrt(pcm)
   ! * Evaluate distribution function
   write(501,*) '#### Real Distribution function for dsigma/dmass(rhoMeson)'
   write(501,*) '#### srts=', s
   do j=0,500
      mass=float(j)*0.01
      if(mass.gt.s-0.938) then
         spectral=0.
         cycle
      end if
      pcm=0.25*(s**2-mass**2+0.938**2)**2/s-0.938**2
      gamtot=FullWidthMeson(rho,mass)
      mres0=meson(rho)%mass
      if(pcm.gt.0) then
         spectral=mass**2*gamtot/  ((mres0**2-mass**2)**2+gamtot**2*mass**2)*sqrt(pcm)
      else 
         spectral=0.
      end if
      write(501,*) mass, spectral/norm
   end do
   close(101)
   close(501)
end do

end program testRhoMass
