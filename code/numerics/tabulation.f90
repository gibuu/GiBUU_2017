!******************************************************************************
!****m* /tabulation
! NAME
! module tabulation
!
! PURPOSE
! Implements multi-dimensional tabulation with constant grid size, single
! precision (4 byte) and linear interpolation.
! Tabulation routines for 1-4 dimensions are available.
!******************************************************************************
module tabulation

  implicit none
  private

!   type table_1dim
!     real(4) :: min,max,delta
!     real(4),dimension(:), Allocatable :: entry
!   end type table_1dim

  type table_2dim
    real(4),dimension(1:2) :: min,max,delta
    real(4),dimension(:,:), Allocatable :: entry
  end type table_2dim

  type table_3dim
    real(4),dimension(1:3) :: min,max,delta
    real(4),dimension(:,:,:), Allocatable :: entry
  end type table_3dim

  type table_4dim
    real(4),dimension(1:4) :: min,max,delta
    real(4),dimension(:,:,:,:), Allocatable :: entry
  end type table_4dim


  Interface getValue
    Module procedure getValue_2dim,getValue_3dim,getValue_4dim  !,getValue_1dim
  end interface

  Interface printTable
    Module procedure print_2dim,print_3dim,print_4dim  !,print_1dim
  end interface

  Interface readTable
    module procedure read_2dim,read_3dim,read_4dim  !,read_1dim
  end interface

  Interface cleanupTable
    module procedure cleanUp_2dim,cleanUp_3dim,cleanUp_4dim  !,cleanUp_1dim
  end interface

!   Interface initTable  ! (ambiguous interface)
!     Module procedure init_table_1dim,init_table_2dim,init_table_3dim,init_table_4dim
!   end interface

  public :: table_2dim, table_3dim, table_4dim  !,table_1dim
  public :: getValue, readTable, printTable, cleanupTable
  public :: init_table_2dim, init_table_3dim, init_table_4dim
  public :: interpolate2, interpolate4

  logical,parameter :: debugflag=.false.

contains

  ! ****************************************************************************************************************
  ! INITIALIZING
  ! ****************************************************************************************************************

!   function init_table_1dim(min,max,delta,f) result(t)
!     type(table_1dim) :: t
!     real(4) :: min,max,delta
!     integer :: numsteps
!     integer :: i
!     interface
!       real function f(x1)
!         real, intent(in) :: x1
!       end function
!     end interface
!
!     t%min=min
!     t%delta=delta
!     numsteps=Nint((max-min)/delta)
!     t%max=min+float(numsteps)*delta
!
!     call cleanUp_1dim(t)
!     Allocate(t%entry(0:numsteps))
!
!     do i=lbound(t%entry,dim=1),ubound(t%entry,dim=1)
!        t%entry(i)=f(t%min+t%delta*float(i))
!     end do
!
!   end function init_table_1dim


  function init_table_2dim(min,max,delta,f) result(t)
    type(table_2dim) :: t
    real,dimension(1:2),intent(in) :: min,max,delta
    interface
      real function f(x1,x2)
        real, intent(in) :: x1,x2
      end function
    end interface

    integer,dimension(1:2) :: numsteps
    integer :: i,j
    real, dimension(1:2) :: x

    t%min=min
    t%delta=delta
    numsteps=Nint((max-min)/t%delta)
    t%max=t%min+float(numsteps)*delta

    call cleanUp_2dim(t)
    allocate(t%entry(0:numsteps(1),0:numsteps(2)))
    call writeMemoryUsage(numSteps(1)*numSteps(2))

    do i=lbound(t%entry,dim=1),ubound(t%entry,dim=1)
       x(1)=min(1)+float(i)*delta(1)
       do j=lbound(t%entry,dim=2),ubound(t%entry,dim=2)
          x(2)=min(2)+float(j)*delta(2)
          t%entry(i,j)=f(x(1),x(2))
       end do
    end do

  end function init_table_2dim



  function init_table_3dim(min,max,delta,f) result(t)
    type(table_3dim) :: t
    real(4),dimension(1:3) :: min,max,delta
    interface
      real function f(x1,x2,x3)
        real, intent(in) :: x1,x2,x3
      end function
    end interface

    integer,dimension(1:3) :: numsteps
    integer :: i,j,k
    real, dimension(1:3) :: x

    t%min=min
    t%delta=delta
    numsteps=Nint((max-t%min)/t%delta)
    t%max=min+float(numsteps)*delta

    call cleanUp_3dim(t)
    allocate(t%entry(0:numsteps(1),0:numsteps(2),0:numsteps(3)))
    call writeMemoryUsage(numSteps(1)*numSteps(2)*numSteps(3))

    do i=lbound(t%entry,dim=1),ubound(t%entry,dim=1)
       x(1)=min(1)+float(i)*delta(1)
       do j=lbound(t%entry,dim=2),ubound(t%entry,dim=2)
          x(2)=min(2)+float(j)*delta(2)
          do k=lbound(t%entry,dim=3),ubound(t%entry,dim=3)
             x(3)=min(3)+float(k)*delta(3)
             t%entry(i,j,k)=f(x(1),x(2),x(3))
          end do
       end do
    end do

  end function init_table_3dim


  function init_table_4dim(min,max,delta,f,printStatus) result(t)
    type(table_4dim) :: t
    real(4),dimension(1:4) :: min,max,delta
    interface
      real function f(x1,x2,x3,x4)
        real, intent(in) :: x1,x2,x3,x4
      end function
    end interface
    logical, optional :: printStatus ! if .true. then print out the status of the initialization

    integer,dimension(1:4) :: numsteps
    integer :: i,j,k,m
    real, dimension(1:4) :: x

    t%min=min
    t%delta=delta
    numsteps=Nint((max-t%min)/t%delta)
    t%max=t%min+float(numsteps)*delta

    call cleanup_4dim(t)
    allocate(t%entry(0:numsteps(1),0:numsteps(2),0:numsteps(3),0:numsteps(4)))

    call writeMemoryUsage(numSteps(1)*numSteps(2)*numSteps(3)*numSteps(4))

    do i=lbound(t%entry,dim=1),ubound(t%entry,dim=1)
       x(1)=min(1)+float(i)*delta(1)
       do j=lbound(t%entry,dim=2),ubound(t%entry,dim=2)
          if (present(printStatus)) then
             if (printStatus) write(*,'(I6,A,I6)') size(t%entry,dim=2)*(i-lbound(t%entry,dim=1))+j,'/', &
                  &  size(t%entry,dim=1)*size(t%entry,dim=2)
          end if
          x(2)=min(2)+float(j)*delta(2)
          do k=lbound(t%entry,dim=3),ubound(t%entry,dim=3)
             x(3)=min(3)+float(k)*delta(3)
             do m=lbound(t%entry,dim=4),ubound(t%entry,dim=4)
                x(4)=min(4)+float(m)*delta(4)
                t%entry(i,j,k,m)=f(x(1),x(2),x(3),x(4))
             end do
          end do
       end do
    end do

  end function init_table_4dim


  ! ****************************************************************************************************************
  ! GETVALUE
  ! ****************************************************************************************************************


!   function getValue_1dim(x,t,outOfBounds) result(value)
!     real(4)         ,intent(in ) :: x
!     type(table_1dim),intent(in)  :: t
!     integer         ,intent(out) :: outOfBounds
!
!     real :: value
!     real :: weight
!     integer :: index
!
!     if(x.le.t%min) then
!        outOfBounds=-1
!        value=t%entry(0)
!     else if(x.ge.t%max) then
!        outOfBounds=1
!        value=t%entry(ubound(t%entry,dim=1))
!     else
!        outOfBounds=0
!        index=int((x-t%min)/t%delta)
!        weight=(x-t%min-float(index)*t% delta)/t% delta
!        value=(1.-weight)*t%entry(index)+weight * t%entry(index+1)
!     end if
!   end function getValue_1dim


  function getValue_2dim(x,t,outOfBounds) result(value)
    real(4),dimension(1:2),intent(in) :: x
    type(table_2dim),intent(in)  :: t
    integer         ,intent(out) :: outOfBounds
    real :: value

    if (x(1).le.t%min(1).or.x(2).le.t%min(2)) then
       outOfBounds=-1
    else if (x(1).ge.t%max(1).or.x(2).ge.t%max(2)) then
       outOfBounds=1
    else
       outOfBounds=0
    end if

    value = interpolate2(t, x)

  end function getValue_2dim


  function getValue_3dim(x,t,outOfBounds) result(value)
    real(4),dimension(1:3),intent(in) :: x
    type(table_3dim),intent(in)  :: t
    integer         ,intent(out) :: outOfBounds
    real :: value

    if (x(1).le.t%min(1).or.x(2).le.t%min(2).or.x(3).le.t%min(3)) then
       outOfBounds=-1
    else if (x(1).ge.t%max(1).or.x(2).ge.t%max(2).or.x(3).ge.t%max(3)) then
       outOfBounds=1
    else
       outOfBounds=0
    end if

    value=interpolate3(x,t%min,t%max,t%delta,t%entry)

  end function getValue_3dim


  function getValue_4dim(x,t,outOfBounds) result(value)
    real(4),dimension(1:4),intent(in) :: x
    type(table_4dim),intent(in)  :: t
    integer         ,intent(out) :: outOfBounds
    real :: value

    if (x(1).le.t%min(1).or.x(2).le.t%min(2).or.x(3).le.t%min(3).or.x(4).le.t%min(4)) then
       outOfBounds=-1
    else if (x(1).ge.t%max(1).or.x(2).ge.t%max(2).or.x(3).ge.t%max(3).or.x(4).ge.t%max(4)) then
       outOfBounds=1
    else
       outOfBounds=0
    end if

    value=interpolate4(x,t%min,t%max,t%delta,t%entry)

  end function getValue_4dim




  ! ****************************************************************************************************************
  ! PRINTING
  ! ****************************************************************************************************************


!   subroutine print_1dim(t,filename)
!     type(table_1dim),intent(in)  :: t
!     character(*) :: filename
!
!     integer :: ios
!     integer :: i
!
!     open(77,file=filename,iostat=ios)
!     if(ios.ne.0) then
!        write(*,*) 'Error opening file', filename
!     end if
!     write(77,'(1E16.8)') t%min
!     write(77,'(1E16.8)') t%delta
!     write(77,'(1I8)')     ubound(t%entry,dim=1)
!
!     do i=lbound(T%entry,dim=1),ubound(T%entry,dim=1)
!        write(77,'(1E16.8)') T%entry(i)
!     end do
!     close(77)
!   end subroutine print_1dim


  subroutine print_2dim(t,filename)
    use output, only:intToChar
    use bzip

    type(table_2dim),intent(in)  :: t
    character(*) :: filename

    integer :: i
    type(bzFile) :: f
    character(len=300) :: line

    f = bzOpenW(filename)
    write(line,'(2E16.8)') t%min
    call bzWriteLine(f,line)
    write(line,'(2E16.8)') t%delta
    call bzWriteLine(f,line)
    write(line,'(2I8)')    ubound(t%entry)
    call bzWriteLine(f,line)

    do i=lbound(T%entry,dim=1),ubound(T%entry,dim=1)
       write(line,'('//intToChar(ubound(T%entry,dim=2)+1)//'E16.8)') T%entry(i,:)
       call bzWriteLine(f,line)
    end do
    call bzCloseW(f)
  end subroutine print_2dim


  subroutine print_3dim(t,filename)
    use output, only:intToChar

    type(table_3dim),intent(in)  :: t
    character(*) :: filename

    integer :: ios,i,j

    open(77,file=filename,iostat=ios)
    if (ios.ne.0) then
       write(*,*) 'Error opening file', filename
    end if
    write(77,'(3E16.8)') t%min
    write(77,'(3E16.8)') t%delta
    write(77,'(3I8)')     ubound(t%entry)

    do i=lbound(T%entry,dim=1),ubound(T%entry,dim=1)
       do j=lbound(T%entry,dim=2),ubound(T%entry,dim=2)
          write(77,'('//intToChar(ubound(T%entry,dim=3)+1)//'E16.8)') T%entry(i,j,:)
       end do
    end do
    close(77)
  end subroutine print_3dim


  subroutine print_4dim(t,filename)
    use output, only:intToChar
    use bzip

    type(table_4dim),intent(in)  :: t
    character(*) :: filename

    integer :: i,j,k
    type(bzFile) :: f
    character(len=300) :: line

    f = bzOpenW(filename)
    write(line,'(4E16.8)') t%min
    call bzWriteLine(f,line)
    write(line,'(4E16.8)') t%delta
    call bzWriteLine(f,line)
    write(line,'(4I8)')     ubound(t%entry)
    call bzWriteLine(f,line)

    do i=lbound(T%entry,dim=1),ubound(T%entry,dim=1)
       do j=lbound(T%entry,dim=2),ubound(T%entry,dim=2)
          do k=lbound(T%entry,dim=3),ubound(T%entry,dim=3)
             write(line,'('//intToChar(ubound(T%entry,dim=4)+1)//'E16.8)') T%entry(i,j,k,:)
             call bzWriteLine(f,line)
          end do
       end do
    end do
    call bzCloseW(f)
  end subroutine print_4dim


  ! ****************************************************************************************************************
  ! READING
  ! ****************************************************************************************************************


!   logical function read_1dim(t,filename,min,max,delta)
!     type(table_1dim),intent(out)  :: t
!     character(*) :: filename
!     real(4) :: min,max,delta
!     real :: eps
!     integer :: ios
!     integer :: i
!     integer :: numsteps
!
!     t%delta=delta
!     t%min=min
!     eps=delta*1E-4
!
!     call cleanup_1dim(t)
!     Allocate(t%entry(0:INT((max-min)/delta)+1))
!     t%max=min+delta*ubound(t%entry,dim=1)
!
!     open(77,file=filename,iostat=ios)
!     if(ios.ne.0) then
!        write(*,*) 'Error opening file ', filename
!        read_1dim=.false.
!        close(77)
!        return
!     end if
!
!
!     read(77,*,iostat=ios) t%min
!     if(ios.ne.0) then
!        write(*,*) 'Error opening file ', filename
!        read_1dim=.false.
!        close(77)
!        return
!     end if
!
!     read(77,*,iostat=ios) t%delta
!     if(ios.ne.0) then
!        write(*,*) 'Error opening file ', filename
!        read_1dim=.false.
!        close(77)
!        return
!     end if
!
!
!     read(77,*,iostat=ios) numSteps
!     if(ios.ne.0) then
!        write(*,*) 'Error opening file ', filename
!        read_1dim=.false.
!        close(77)
!        return
!     end if
!
!
!     if(abs(t%min-min).gt.eps) then
!        write(*,*) 'min does not agree in file ', filename
!        read_1dim=.false.
!        write(*,*) min, t%min
!        close(77)
!        return
!     end if
!     if( abs(ubound(t%entry,dim=1)- numsteps).gt.eps) then
!        write(*,*) 'numsteps do not agree in file ', filename
!        write(*,*) numsteps, ubound(t%entry,dim=1)
!        read_1dim=.false.
!        close(77)
!        return
!     end if
!     if(abs(t%delta-delta).gt.eps) then
!        write(*,*) 'delta does not agree in file ', filename
!        write(*,*) delta, t%delta
!        read_1dim=.false.
!        close(77)
!        return
!     end if
!
!     do i=lbound(T%entry,dim=1),ubound(T%entry,dim=1)
!        read(77,*,iostat=ios) T%entry(i)
!        if(ios.ne.0) then
!           write(*,*) 'Error opening file ', filename
!           read_1dim=.false.
!           close(77)
!           return
!        end if
!     end do
!     read_1dim=.true.
!     close(77)
!
!   end function read_1dim


  logical function read_2dim(t,filename)
    use bzip

    type(table_2dim),intent(out)  :: t
    character(*) :: filename
!     real(4),dimension(1:2) :: minimum,max,delta

    integer :: ios,i,ll
    integer,dimension(1:2) :: numsteps
    type(bzFile) :: f
    character(len=300) :: line

    call cleanup_2dim(t)

    f = bzOpenR(filename)

    ll = 0
    call bzReadLine(f,line,ll)
    read(line(1:ll),*,iostat=ios) t%min
    if (ios.ne.0) then
       write(*,*) 'Error reading file ', filename
       read_2dim=.false.
       call bzCloseR(f)
       return
    end if

    ll = 0
    call bzReadLine(f,line,ll)
    read(line(1:ll),*,iostat=ios) t%delta
    if (ios.ne.0) then
       write(*,*) 'Error reading file ', filename
       read_2dim=.false.
       call bzCloseR(f)
       return
    end if

    ll = 0
    call bzReadLine(f,line,ll)
    read(line(1:ll),*,iostat=ios) numSteps
    if (ios.ne.0) then
       write(*,*) 'Error reading file ', filename
       read_2dim=.false.
       call bzCloseR(f)
       return
    end if

    t%max = t%min + t%delta * numSteps
    allocate (t%entry (0:numSteps(1),0:numSteps(2)))

    ll = 0
    do i=lbound(T%entry,dim=1),ubound(T%entry,dim=1)
       call bzReadLine(f,line,ll)
       read(line(1:ll),*,iostat=ios) T%entry(i,:)
       if (ios.ne.0) then
          write(*,*) 'Error reading file ', filename
          read_2dim=.false.
          call bzCloseR(f)
          return
       end if
    end do

    read_2dim=.true.
    call bzCloseR(f)
    call writeMemoryUsage(numSteps(1)*numSteps(2))

  end function read_2dim


  logical function read_3dim(t,filename,minimum,max,delta)
    type(table_3dim),intent(out)  :: t
    character(*) :: filename
    real(4),dimension(1:3) :: minimum,max,delta

    real :: eps
    integer :: ios,i,j
    integer,dimension(1:3) :: numsteps

    t%delta=delta
    t%min=minimum
    eps=minval(delta)*1E-4

    call cleanup_3dim(t)
    allocate(t%entry(0:INT((max(1)-minimum(1))/delta(1))+1 &
         & ,0:INT((max(2)-minimum(2))/delta(2))+1&
         & ,0:INT((max(3)-minimum(3))/delta(3))+1))

    do i=lbound(minimum,dim=1),ubound(minimum,dim=1)
       t%max(i)=minimum(i)+delta(i)*ubound(t%entry,dim=i)
    end do

    open(77,file=filename,iostat=ios)
    if (ios.ne.0) then
       write(*,*) 'Error opening file ', filename
       read_3dim=.false.
       close(77)
       return
    end if

    read(77,*,iostat=ios) t%min
    if (ios.ne.0) then
       write(*,*) 'Error opening file ', filename
       read_3dim=.false.
       close(77)
       return
    end if

    read(77,*,iostat=ios) t%delta
    if (ios.ne.0) then
       write(*,*) 'Error opening file ', filename
       read_3dim=.false.
       close(77)
       return
    end if

    read(77,*,iostat=ios) numSteps
    if (ios.ne.0) then
       write(*,*) 'Error opening file ', filename
       read_3dim=.false.
       close(77)
       return
    end if

    do i=lbound(minimum,dim=1),ubound(minimum,dim=1)
       if (abs(t%min(i)-minimum(i)).gt.eps) then
          write(*,*) 'minimum(',i,') does not agree in file ', filename
          read_3dim=.false.
          write(*,*) minimum(i), t%min(i)
          close(77)
          return
       end if
       if ( abs(ubound(t%entry,dim=i)- numsteps(i)).gt.eps) then
          write(*,*) 'numsteps(',i,') do not agree in file ', filename
          write(*,*) numsteps(i), ubound(t%entry,dim=i)
          read_3dim=.false.
          close(77)
          return
       end if
       if (abs(t%delta(i)-delta(i)).gt.eps) then
          write(*,*) 'delta(',i,') does not agree in file ', filename
          write(*,*) delta(i), t%delta(i)
          read_3dim=.false.
          close(77)
          return
       end if
    end do

    do i=lbound(T%entry,dim=1),ubound(T%entry,dim=1)
       do j=lbound(T%entry,dim=2),ubound(T%entry,dim=2)
          read(77,*,iostat=ios) T%entry(i,j,:)
          if (ios.ne.0) then
             write(*,*) 'Error opening file ', filename
             read_3dim=.false.
             close(77)
             return
          end if
       end do
    end do
    read_3dim=.true.
    close(77)
    call writeMemoryUsage(numSteps(1)*numSteps(2)*numSteps(3))

  end function read_3dim


  logical function read_4dim(t,filename,minimum,delta)
    use bzip

    type(table_4dim),intent(out)  :: t
    character(*) :: filename
    real(4),dimension(1:4) :: minimum,delta

    real :: eps
    integer :: ios,i,j,k,ll
    integer,dimension(1:4) :: numsteps
    type(bzFile) :: f
    character(len=300) :: line

    call cleanup_4dim(t)

    eps=minval(delta)*1E-4

    f = bzOpenR(filename)

    ll = 0
    call bzReadLine(f,line,ll)
    read(line(1:ll),*,iostat=ios) t%min
    if (ios.ne.0) then
       write(*,*) 'Error opening file ', filename
       read_4dim=.false.
       call bzCloseR(f)
       return
    end if

    ll = 0
    call bzReadLine(f,line,ll)
    read(line(1:ll),*,iostat=ios) t%delta
    if (ios.ne.0) then
       write(*,*) 'Error opening file ', filename
       read_4dim=.false.
       call bzCloseR(f)
       return
    end if

    ll = 0
    call bzReadLine(f,line,ll)
    read(line(1:ll),*,iostat=ios) numSteps
    if (ios.ne.0) then
       write(*,*) 'Error opening file ', filename
       read_4dim=.false.
       call bzCloseR(f)
       return
    end if

    do i=lbound(minimum,dim=1),ubound(minimum,dim=1)
       if (abs(t%min(i)-minimum(i)).gt.eps) then
          write(*,*) 'minimum(',i,') does not agree in file ', filename
          read_4dim=.false.
          write(*,*) minimum(i), t%min(i)
          call bzCloseR(f)
          return
       end if
!       if( abs(ubound(t%entry,dim=i)- numsteps(i)).gt.eps) then
!          write(*,*) 'numsteps(',i,') do not agree in file ', filename
!          write(*,*) numsteps(i), ubound(t%entry,dim=i)
!          read_4dim=.false.
!          call bzCloseR(f)
!          return
!       end if
       if (abs(t%delta(i)-delta(i)).gt.eps) then
          write(*,*) 'delta(',i,') does not agree in file ', filename
          write(*,*) delta(i), t%delta(i)
          read_4dim=.false.
          call bzCloseR(f)
          return
       end if
    end do

    allocate(t%entry(0:numsteps(1),0:numsteps(2),0:numsteps(3),0:numsteps(4)))
    t%max=t%min+float(numsteps)*delta

    ll = 0
    do i=lbound(T%entry,dim=1),ubound(T%entry,dim=1)
       do j=lbound(T%entry,dim=2),ubound(T%entry,dim=2)
          do k=lbound(T%entry,dim=3),ubound(T%entry,dim=3)
             call bzReadLine(f,line,ll)
             read(line(1:ll),*,iostat=ios) T%entry(i,j,k,:)
             if (ios.ne.0) then
                write(*,*) 'Error opening file ', filename
                read_4dim=.false.
                call bzCloseR(f)
                return
             end if
          end do
       end do
    end do
    read_4dim=.true.
    call bzCloseR(f)
    call writeMemoryUsage(numSteps(1)*numSteps(2)*numSteps(3)*numSteps(4))

  end function read_4dim


  ! ****************************************************************************************************************
  ! CLEANUP
  ! ****************************************************************************************************************

!   subroutine cleanUp_1dim(t)
!     type(table_1dim),intent(inout)  :: t
!     if(allocated(t%entry)) DeAllocate(t%entry)
!   end subroutine

  subroutine cleanUp_2dim(t)
    type(table_2dim),intent(inout)  :: t
    if (allocated(t%entry)) deallocate(t%entry)
  end subroutine

  subroutine cleanUp_3dim(t)
    type(table_3dim),intent(inout)  :: t
    if (allocated(t%entry)) deallocate(t%entry)
  end subroutine

  subroutine cleanUp_4dim(t)
    type(table_4dim),intent(inout)  :: t
    if (allocated(t%entry)) deallocate(t%entry)
  end subroutine


  ! ****************************************************************************************************************
  ! INTERPOLATION
  ! ****************************************************************************************************************

  function interpolate2(t,X_in) result(interValue)
    type(table_2dim),intent(in) :: t
    real(4), dimension(1:2),intent(in) :: X_in
    real(4) :: interValue

    real(4), dimension(1:2) :: X
    real(4),dimension (1:2,1:2) :: weights
    integer,dimension(1:2,1:2) :: index
    integer :: j,i1,i2

    if (debugflag) then
      do j=1,2
         if (x_in(j).gt.t%max(j) .or. x_in(j).lt.t%min(j)) then
            write(*,*) 'WARNING in interpolate2: Out of bounds!'
            write(*,*) j,x_in(j),t%min(j),t%max(j)
            !stop
         end if
      end do
    end if

    x = min(t%max,max(x_in,t%min))

    ! Get smallest index which is closest to x
    index(:,1) = min(lbound(t%entry)+int((X-t%min)/t%delta),ubound(t%entry)-1)
    index(:,2) = index(:,1)+1

    ! calculate weights
    weights(:,1) = 1.-(abs(X-(t%min+t%delta*float(index(:,1)-lbound(t%entry))))/t%delta)
    weights(:,2) = 1.-weights(:,1)

    ! Evaluate interpolation
    interValue=0.
    do i2=1,2
      do i1=1,2
        interValue = interValue + weights(1,i1) * weights(2,i2) * t%entry(index(1,i1),index(2,i2))
      end do
    end do

  end function interpolate2


  function interpolate3(X_in,Xmin,Xmax,delta,Table) result(interValue)
    real(4), dimension(:,:,:),intent(in) :: table
    real(4), dimension(1:3),intent(in) :: X_in,Xmin,Xmax,delta
    real(4) :: interValue

    real(4), dimension(1:3) :: X
    real(4),dimension (1:3,1:2) :: weights
    integer,dimension(1:3,1:2) :: index
    integer :: j,i1,i2,i3

    if (debugflag) then
      do j=1,3
        if (x_in(j).gt.Xmax(j) .or. x_in(j).lt.Xmin(j)) then
            write(*,*) 'WARNING in interpolate3: Out of bounds!'
            write(*,*) j,x_in(j),Xmin(j),Xmax(j)
            !stop
         end if
      end do
    end if

    x=min(Xmax,max(x_in,Xmin))

    ! Get smallest index which is closest to x
    index(:,1)=min(lbound(table)+int((X-Xmin)/delta),ubound(table)-1)
    index(:,2)=index(:,1)+1

    ! calculate weights
    weights(:,1)=1.-(abs(X-(Xmin+delta*float(index(:,1)-lbound(table))))/delta)
    weights(:,2)=1.-weights(:,1)

    ! Evaluate interpolation
    interValue=0.
      do i3=1,2
        do i2=1,2
          do i1=1,2
            interValue = interValue + weights(1,i1) * weights(2,i2) * weights(3,i3) * &
                         table(index(1,i1),index(2,i2),index(3,i3))
          end do
        end do
      end do

  end function interpolate3


  !****************************************************************************
  !****f* tabulation/interpolate4
  ! NAME
  ! function interpolate4(X,min,max,delta,Table) result(interValue)
  !
  ! PURPOSE
  ! Implements 4-dimensional linear interpolations for constant grid distances.
  ! Assume that f(x)=y where x 4-dim and y 1-dim. Given a table of function values
  ! f_{m,n,o,p}=f(x_1=min(1)+m*delta(1),x_2=min(2)+n*delta(2),x_3=...,x_4=...), this function estimates
  ! the value of f(X)  where X is not a grid-point.
  !
  ! INPUTS
  ! * real, dimension(:,:),intent(in) :: table         -- Table of function values f_{m,n,o,p}
  ! * real, dimension(1:4),intent(in) :: min,max,delta -- Minimum, maximum and step size of the underlying x-grid
  ! * real, dimension(1:4),intent(in) :: X             -- X-Value at which we want to know f(X)
  !
  ! OUTPUT
  ! * real :: interValue
  !
  ! NOTES
  ! If the point x lies in one of its dimensions outside the grid,
  ! then we take either the lowest or the highest grid point of this
  ! dimension to evaluate the interpolated value.
  !****************************************************************************
  function interpolate4(X_in,Xmin,Xmax,delta,Table) result(interValue)
    real(4), dimension(:,:,:,:),intent(in) :: table
    real(4), dimension(1:4),intent(in) :: X_in,Xmin,Xmax,delta
    real(4) :: interValue

    real(4), dimension(1:4) :: X
    real(4),dimension (1:4,1:2) :: weights
    integer,dimension(1:4,1:2) :: index
    integer :: j,i1,i2,i3,i4

    if (debugflag) then
       do j=1,4
          if (x_in(j).gt.Xmax(j) .or. x_in(j).lt.Xmin(j)) then
             write(*,*) 'WARNING in interpolate4: Out of bounds!'
             write(*,*) j,x_in(j),Xmin(j),Xmax(j)
             !stop
          end if
       end do
    end if

    x=min(Xmax,max(x_in,Xmin))

    ! Get smallest index which is closest to x
    index(:,1)=min(lbound(table)+int((X-Xmin)/delta),ubound(table)-1)
    index(:,2)=index(:,1)+1

    ! calculate weights
    weights(:,1)=1.-(abs(X-(Xmin+delta*float(index(:,1)-lbound(table))))/delta)
    weights(:,2)=1.-weights(:,1)

    ! Evaluate interpolation
    interValue=0.
    do i4=1,2
      do i3=1,2
        do i2=1,2
          do i1=1,2
            interValue = interValue + weights(1,i1) * weights(2,i2) * weights(3,i3) * weights(4,i4) * &
                         table(index(1,i1),index(2,i2),index(3,i3),index(4,i4))
          end do
        end do
      end do
    end do

  end function interpolate4


  subroutine writeMemoryUsage(entries)
    integer :: entries
    ! 8 Byte per real:
    write(*,'(A,I8,A,G12.2,A)') 'A tabulated field with ',entries,' entries of ', entries*8./1024.**2  &
         & , ' MB has been initialized!'
  end subroutine writeMemoryUsage


end module tabulation
