program vesmes_test

use particleProperties, only: initParticleProperties
use IdTable, only: nucleon, Delta, rho, omegaMeson, phi, eta
use inputGeneral, only: readInputGeneral
use mediumDefinition
use photonXSections, only: calcXS_gammaN2VN, calcXS_gammaN2VDelta, setIParam
use output
use particleDefinition
use resprod_lepton, only: sigma_barMes_res_vac
use minkowski, only: abs4
use XS_VMD, only: vmd,gvmd
use constants, only: mN

implicit none

real:: E_gamma,srts,mom_gamma(0:3)
integer :: i,j
real:: clock_start, clock_finish ! time measurement
real :: XS_VN(0:3), XS_VNpi(0:3), XS_VD(0:3), XS_Res_VN(1:4),XS_Res_VD(1:3),XSvmd(0:4),XSgvmd(0:4)
type(medium)  :: med
type(particle) :: nuc

call readInputGeneral
call initParticleProperties

! set up proton at rest
nuc%ID = nucleon
nuc%mass = mN
nuc%momentum = (/mN, 0., 0., 0./)
nuc%charge = 1

call cpu_time(clock_start)

do i = 1,3  ! loop over different parametrizations of the high energy behaviour
   open(unit=10,file="vecmes"//inttochar(i)//".dat")
   call setIParam(i,.true.)
   do j=1,2000
      E_gamma = j*0.01
      mom_gamma = (/E_gamma, 0., 0., E_gamma/)
      srts = abs4(nuc%momentum+mom_gamma)
      call calcXS_gammaN2VN (srts, med, XS_VN, XS_VNpi)
      call calcXS_gammaN2VDelta (srts, XS_VD, med)
      XS_Res_VN(1:4) = sigma_barMes_res_vac(nuc,mom_gamma,nucleon,(/rho,omegaMeson,phi,eta/))*1000.
      XS_Res_VD(1:3) = sigma_barMes_res_vac(nuc,mom_gamma,Delta,(/rho,omegaMeson,phi/))*1000.
      call vmd(srts,0.,0.,XSvmd)
      call gvmd(srts,0.,0.,XSgvmd)
      write(10,"(22F12.6)") E_gamma,srts,XS_VN,XS_VNpi,XS_VD,XS_Res_VN(1:3),XS_Res_VD(1:3),XS_Res_VN(4),XSvmd(0)+XSgvmd(0)
   end do
  call cpu_time(clock_finish)
  close(10)
  write(*,*) "--------------------------------------------------"
  write(*,*) "CPU Time = ",time_format(clock_finish - clock_start)
  write(*,*) "--------------------------------------------------"
end do

contains

function time_format(t)
   real,intent(in)::t
   character::time_format*13
   integer:: ms,s,m,h,d
   ms=mod(int(t*10),10)
   s=mod(int(t),60)
   m=mod(int(t/60),60)
   h=mod(int(t/60/60),24)
   d=int(t/60/60/24)
   write(time_format,"(i2,a,i2.2,a,i2.2,a,i2.2,a,i1.1)") d,":",h,":",m,":",s,".",ms
end function

end program
