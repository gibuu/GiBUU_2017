! This test case is used to create a 2-dimensional tabulation of nucleon-nucleon
! bremsstrahlung results from an OBE calculation by Radhey Shyam. These results
! are given by several tables of dsigma/dm for different projectile kinetic
! energies (0.25 - 3.5 GeV), one file per energy, which are read into a 2-dim.
! table and exported into a single file.
!
! We create three tables: one for pp collisions and two for pn (with and without
! pion electromagnetic form factor).


program createBremsTable
  use tabulation, only: table_2dim
  implicit none

  call createTable('pp')        ! proton-proton
  call createTable('pn')        ! proton-neutron (without pion FF)
  call createTable('pn-piFF')   ! proton-neutron (with pion FF)

contains


  subroutine createTable (str)
    use tabulation, only: printTable
    
    character(len=*), intent(in) :: str

    type(table_2dim) :: tab
    integer :: numsteps(1:2),i

    ! first component: invariant mass, second component: beam energy
    tab%min   = (/ 0.010, 0.25 /)
    tab%max   = (/ 1.274, 3.50 /)
    tab%delta = (/ 0.001, 0.25 /)

    numsteps = nint((tab%max-tab%min)/tab%delta)
    Allocate(tab%entry(0:numsteps(1),0:numsteps(2)))

    do i=1,14
      call readFile (tab, str, i)
    end do

    call printTable (tab, "brems-"//trim(str)//".bz2")
  
  end subroutine


  subroutine readFile (tab, str, i)
    type(table_2dim) :: tab
    character(len=*), intent(in) :: str
    integer, intent(in) :: i

    character(len=100) :: filename
    integer :: ios, j
    real :: m, itg

    write(filename,'(A,f4.2)') "/home/jweil/Radhey/tables/"//trim(str)//"/brems-", i*0.25
    open(88,file=filename)
    itg = 0.
    do j=0,ubound(tab%entry,dim=1)
      read(88,*,iostat=ios) m, tab%entry(j,i-1)
      if (m-(tab%min(1)+j*tab%delta(1))>1E-4) print *,"warning: ",i,j,m
      if (ios/=0) exit
      itg = itg + tab%entry(j,i-1) * tab%delta(1)
    end do
    close(88)
    print *,"read: ",i,i*tab%delta(2),j,itg
  end subroutine


end program
