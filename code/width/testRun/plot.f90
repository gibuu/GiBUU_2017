program plot
implicit none
integer :: i
open(1,File="plot.gnu",Status="Unknown")
Write(1,'(A)') 'set logscale y' 
Do i=1,65
Write(1,'(3(A,I3),A)') ' plot  [][0.01:100]"fort.',i+100,'" w p , "old/fort.',i+100,'" w l '
!Write(1,'(3(A,I3),A)') ' plot [][0.01:100] "fort.',i+100,'" w p t "ID=',i,'"' 
write(1,*) "pause(-1)"
End do

end program plot
