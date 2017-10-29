program aaa

  character*(50) :: BUF

!!$!  do iEn=50,-4,-1
!!$  do iEn=50,0,-1
!!$     Ebeam = 10d0**(iEn*0.05d0)
!!$     write(*,'(i5.3,f12.5)') iEn,Ebeam
!!$  end do

  do iEn=20,0,-2
     Ebeam = 10d0**(iEn*0.05d0)
!     write(*,'(i5.3,f12.5)') iEn,Ebeam
     write(BUF,'(1P,e11.5)') Ebeam

     write(12,'(A,A,A,i3.3)') 'sed s/XXX/',trim(BUF),'/ jobDISXXX > jobDIS',iEn

     write(13,'(A,i3.3,A,i3.3,A,A,i3.3,A)') &
          & '../PlotEth.x < jobDIS',iEn,&
          & '; cp Neutrino.h2D.X_Y.dat Neutrino.h2D.X_Y.run',iEn,'.dat',&
          & '; cp Neutrino.h2D.X_Y.dat.bin Neutrino.h2D.X_Y.run',iEn,'.dat.bin'

     write(14,'(A,i3.3,A,A,A,i3.3,A)') &
          & 'cp Neutrino.h2D.X_Y.run',iEn,'.dat Neutrino.h2D.X_Y.dat',&
          & '; gnuplot CreateContDat.gp',&
          & '; cp ContDat.dat ContDat.run',iEn,'.dat'


  end do

end program aaa
