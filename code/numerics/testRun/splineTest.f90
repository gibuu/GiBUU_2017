program test

call splineTest


contains


  subroutine splineTest
    use spline
    implicit none
    real, dimension(1:5) :: y, y2
    real, dimension(1:5) :: x
    real :: value, res
    integer :: i
    x=(/ 1.,2.,12.,30.,50. /)
    y=(/ 1.,3.,126.,200.,-50.1 /)
    ! spline initialisieren
    call spline2(x,y,y2)
    ! schreibe Originalwerte nach fort.201
    Do i=lBound(y,dim=1),uBound(y,dim=1) 
       write(201,*) x(i), y(i)
    End do
    ! schreibe spline nach fort.200
    Do i=-100,100
       value=0.6*i
       call splint2(x,y,y2,value,res)
       write(200,*) value, res
    End do
  end subroutine splineTest

end program test
