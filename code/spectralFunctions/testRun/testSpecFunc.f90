program test

use inputGeneral, only: readinputGeneral
use IDTABLE,only : delta,nucleon
use output, only: intToChar, realToChar4
use particleProperties, only: initParticleProperties, hadron, isNonExotic
use spectralFunc, only : specFunc
use baryonWidthMedium_tables, only : get_deltaOset
use potentialModule, only: scapot

implicit none
integer :: ID, i
integer, parameter :: charge=1
real, parameter :: rho=0.00000000000
real :: mass,pAbs,dE,E,m_min,spec
real, parameter :: m_max=10.5
integer, parameter :: numPoints_p=10
real , dimension(0:numPoints_p) :: int
real, dimension(1:3) :: position
real, dimension(0:3) :: p


call readinputGeneral
call initParticleProperties


dE=0.00025
open(12,File='norm.dat')
write(12,'(A,F9.3,A)') '# Normalization of spectral function. Upper integration boundary=',m_max,'GeV'
write(12,*) '# id; p= 0.2,0.4.,0.6,...'

position=0.

do id=nucleon,delta
   if(.not.isNonExotic(id)) exit
   open(22,File='norm'//intToChar(id)//'.dat')
   write(22,'(A,F9.3,A)') '# Normalization of spectral function. Upper integration boundary=',m_max,'GeV'
   write(22,*) '# id=', id

   m_min=hadron(id)%minmass+0.01
   if(id.eq.nucleon)    m_min=0.7

   do i=0,numPoints_p
      pAbs=float(i)*0.2
      if (get_deltaOset().and.id.eq.delta.and.pabs.gt.1.2) exit
      int(i)=0
      E=sqrt(0.1**2+pabs**2)
      energyLoop : do
         mass=sqrt(E**2-pabs**2)
         if (mass.gt.m_max) exit energyLoop
         if (get_deltaOset().and.id.eq.delta.and.mass.gt.2.5) exit energyLoop
         !write (*,'(A,F9.3,A,F9.3,A,I3,A,F9.3,F6.2,A,F15.4)') &
         !     & 'E=',E,' Mass=',mass,' ID=',id,' p=',pabs, mass/m_max,'%. Integral=',int(i)
         p=(/E,pabs,0.,0./)
         !        write(*,*) Id,charge,p,position
         spec=specFunc(ID,charge,p,position)
         write(14,*) pabs,E, spec
         int(i)=int(i)+spec*2*E*dE
         E=E+dE
      end do energyLoop
      write(14,*)
      write(22,*) pabs, int(i)
   end do

   close(22)
   write(12,'(I4,'//intToChar(numPoints_p+1)//'F9.4)')  id,int

   ! Plots: less points, only some momenta

   do i=0,6
      pAbs=float(i)*0.4
      open(14,File='specFunc_'//intToChar(id)//'.pabs_'//realToChar4(pabs*1000.)//'MeV.dat')
      int(i)=0
      E=sqrt(0.1**2+pabs**2)
      energyLoop_plot : do
         mass=sqrt(E**2-pabs**2)
         if (mass.gt.m_max) exit energyLoop_plot
         p=(/E,pabs,0.,0./)
         spec=specFunc(ID,charge,p,position)
         write(14,*) pabs,E, spec
         int(i)=int(i)+spec*2*E*dE
         E=E+4.*dE
      end do energyLoop_plot
      close(14)
   end do

end do


end program test
