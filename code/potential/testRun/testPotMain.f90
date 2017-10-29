program testPotMain
! After executing testPotMain.x < job in potential/testRun all potentials of the different particle species are plotted
! to files ("fort.*").
! fort.400+id : baryon potentials
! fort.id : baryon potentials
! Loading "plotPotentials.gnu" in gnuplot with "load plotPotentials.gnu" will show all potentials as 3-D-plots.

use particleDefinition
use potentialModule, only: potential_LRF
use IdTable
use particleProperties, only: initParticleProperties, hadron

implicit none
type(particle) :: teilchen

real :: dens,mom
integer :: index_dens
integer :: index_mom
integer :: id
integer, parameter :: densPoints=30  !number of points in density
integer, parameter :: momPoints=50   !number of points in momentum

Open(99,file='plotPotentials.gnu')
write(99,*) 'set dgrid3d',densPoints,",",momPoints
write(99,*) "set xlabel 'density[fm^-3]'"
write(99,*) "set ylabel 'momentum[GeV]'"
write(99,*) "set zlabel 'potential[GeV]'"
write(99,*) "set dgrid3d 10,30"
write(99,*) "set contour base"

call initParticleProperties

do id=pion, dSStar_minus
   write(id,*) '#',hadron(id)%name,"dens, mom, potential_LRF(teilchen,dens)"
   write(99,'(A,I3,A)') 'set title "meson #',id,'"'
   write(99,'(A,I3,A)') 'splot "fort.',id,'" w l'
   do index_dens=0,densPoints-1
      dens=float(index_dens)*0.175/float(densPoints-1)
      do index_mom=0,momPoints
         mom=float(index_mom)*0.5/float(momPoints-1)
         teilchen%Id=id
         teilchen%mass=hadron(id)%mass
         teilchen%momentum(1:3)=(/mom,0.,0./)
         teilchen%momentum(0)=FreeEnergy(teilchen)
         write(id,*) dens, mom, potential_LRF(teilchen,dens)
      end do
      write(id,*) 
   end do
   write(99,'(A)') 'pause 2'
end do

do id=nucleon, Omega_c
   write(400+id,*) '#',hadron(id)%name,"dens, mom, potential_LRF(teilchen,dens)"
   write(99,'(A,I3,A)') 'set title "baryon #',id,'"'
   write(99,'(A,I3,A)') 'splot "fort.',400+id,'" w l'
   do index_dens=0,densPoints-1
      dens=index_dens*1.*0.5/float(densPoints-1)
      do index_mom=0,momPoints
         mom=index_mom*1.0/float(momPoints-1)
         teilchen%Id=id
         teilchen%mass=hadron(id)%mass
         teilchen%momentum(1:3)=(/mom,0.,0./)
         teilchen%momentum(0)=FreeEnergy(teilchen)
         write(400+id,*) dens, mom, potential_LRF(teilchen,dens)
      end do
      write(400+id,*) 
   end do
   write(99,'(A)') 'pause 2'
end do

close(99)

end program testPotMain
