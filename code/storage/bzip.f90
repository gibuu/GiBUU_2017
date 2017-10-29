!******************************************************************************
!****m* /bzip
! NAME
! module bzip
! NOTES
! This module provides routines for reading and writing bz2-compressed files.
!******************************************************************************
module bzip

  use, intrinsic :: ISO_C_BINDING, only: c_int, c_char, c_ptr, c_null_ptr, c_associated, c_new_line

  implicit none
  private

  public :: bzOpenR, bzOpenW
  public :: bzReadLine, bzWriteLine
  public :: bzCloseR, bzCloseW


  !****************************************************************************
  !****t* bzip/bzFile
  ! SOURCE
  !
  type, public :: bzFile
    type(c_ptr) :: handle, bzHandle           ! file handles
    character(kind=c_char,len=100) :: fname   ! file name
    logical :: eof = .false.                  ! end of file reached?
  end type
  !
  ! NOTES
  ! This is a data type to represent bzip files, including file name and file handles.
  ! 'handle' is used by the stdio functions 'fopen' and 'close'.
  ! 'bzHandle' is used by the libbz2 functions 'BZ2_bzReadOpen', 'BZ2_bzRead', 'BZ2_bzReadClose' etc.
  !****************************************************************************


  interface

    ! Fortran interfaces for the C functions 'fopen' and 'fclose', defined in stdio.h

    type(c_ptr) function fOpen (path, mode) bind (c, NAME='fopen')
      use, intrinsic :: ISO_C_BINDING
      character(kind=c_char) :: path(*), mode(*)
    end function

    integer(c_int) function fClose (file) bind (c, NAME='fclose')
      use, intrinsic :: ISO_C_BINDING
      type(c_ptr), value :: file
    end function

    ! Fortran interfaces for the functions from libbz2, defined in bzlib.h

    type(c_ptr) function bzReadOpen (bzerror, f, verbosity, small, unused, nUnused) bind(c, NAME='BZ2_bzReadOpen')
      use, intrinsic :: ISO_C_BINDING
      integer(c_int) :: bzerror
      type(c_ptr), value :: f, unused
      integer(c_int), value :: verbosity, small, nUnused
    end function

    integer(c_int) function bzRead (bzerror, f, buf, len) bind (c, NAME='BZ2_bzRead')
      use, intrinsic :: ISO_C_BINDING
      integer(c_int) :: bzerror
      type(c_ptr), value :: f
      character(kind=c_char) :: buf(*)
      integer(c_int), value :: len
    end function

    subroutine bzReadClose (bzerror, f) bind (c, NAME='BZ2_bzReadClose')
      use, intrinsic :: ISO_C_BINDING
      integer(c_int) :: bzerror
      type(c_ptr), value :: f
    end subroutine

    type(c_ptr) function bzWriteOpen (bzerror, f, blockSize100k, verbosity, workFactor) bind (c, NAME='BZ2_bzWriteOpen')
      use, intrinsic :: ISO_C_BINDING
      integer(c_int) :: bzerror
      type(c_ptr), value :: f
      integer(c_int), value :: blockSize100k, verbosity, workFactor
    end function

    subroutine bzWrite (bzerror, f, buf, len) bind (c, NAME='BZ2_bzWrite')
      use, intrinsic :: ISO_C_BINDING
      integer(c_int) :: bzerror
      type(c_ptr), value :: f
      character(kind=c_char) :: buf(*)
      integer(c_int), value :: len
    end subroutine

    subroutine bzWriteClose (bzerror, f, abandon, nbytes_in, nbytes_out) bind (c, NAME='BZ2_bzWriteClose')
      use, intrinsic :: ISO_C_BINDING
      integer(c_int) :: bzerror, nbytes_in, nbytes_out
      type(c_ptr), value :: f
      integer(c_int), value :: abandon
    end subroutine

  end interface


contains

  !****************************************************************************
  !****f* bzip/bzOpenR
  ! NAME
  ! function bzOpenR (fname)
  ! PURPOSE
  ! Open a bz2 file for reading.
  ! INPUTS
  ! character(len=*) :: fname   --- full file name
  ! OUTPUT
  ! type(bzFile) :: bzOpenR     --- bzFile structure containing valid file handles (if successful)
  !****************************************************************************
  type(bzFile) function bzOpenR (fname)
    character(len=*), intent(in) :: fname

    character(kind=c_char,len=2) :: mode = "r"//achar(0)
    type(c_ptr) :: u = c_null_ptr
    integer(c_int) :: err

    bzOpenR%fname=trim(fname)//achar(0)
    bzOpenR%handle = fOpen(bzOpenR%fname,mode)
    if (.not. c_associated(bzOpenR%handle)) then
      write(*,*) "Error opening file ",bzOpenR%fname
      stop
    end if

    bzOpenR%bzHandle = bzReadOpen(err,bzOpenR%handle,0,0,u,0)
    if (.not. c_associated(bzOpenR%bzHandle) .or. err/=0) then
      write(*,*) "Error opening bz file ",bzOpenR%fname
      write(*,*) "Error code: ",err
      stop
   end if

  end function


  !****************************************************************************
  !****f* bzip/bzOpenW
  ! NAME
  ! function bzOpenW (fname)
  ! PURPOSE
  ! Open a bz2 file for writing.
  ! INPUTS
  ! character(len=*) :: fname --- full file name
  ! OUTPUT
  ! type(bzFile) :: bzOpenW   --- bzFile structure containing valid file handles (if successful)
  !****************************************************************************
  type(bzFile) function bzOpenW (fname)
    character(len=*), intent(in) :: fname

    character(kind=c_char,len=2) :: mode = "w"//achar(0)
    integer(c_int) :: err

    bzOpenW%fname=trim(fname)//achar(0)
    bzOpenW%handle = fOpen(bzOpenW%fname,mode)
    if (.not. c_associated(bzOpenW%handle)) then
      write(*,*) "Error opening file ",bzOpenW%fname
      stop
    end if

    bzOpenW%bzHandle = bzWriteOpen(err,bzOpenW%handle,9,0,0)
    if (.not. c_associated(bzOpenW%bzHandle) .or. err/=0) then
      write(*,*) "Error opening bz file ",bzOpenW%fname
      write(*,*) "Error code: ",err
      stop
   end if

  end function


  !****************************************************************************
  !****s* bzip/bzReadline
  ! NAME
  ! subroutine bzReadline (f, buf, length, debug)
  ! PURPOSE
  ! Read one line from a bz2 file.
  ! INPUTS
  ! * type(bzFile) :: f --- bzFile structure (which should have been obtained from a previous bzOpenR)
  ! * integer :: length --- minimal length of the line to read (excluding '\n')
  ! OUTPUT
  ! * character(kind=C_CHAR,len=*) :: buf --- buffer containing the data which has been read from the file
  ! * integer,intent(inout) :: length --- actual length of the line (excluding '\n')
  ! NOTES
  ! The input "length" must be set before calling this routine. Smaller values provide more safety and flexibility,
  ! while larger value provide better performance. If the line length is unknown, this can just be set to zero,
  ! and the line length will be determined automatically (note: this is slow). If the line length is known,
  ! the fastest method is to set "length" to the exact line length (excluding the newline character '\n').
  ! "Length" can also be set to any value between zero and the actual line length (as a trade-off between performance and flexibility).
  ! However, it should never be set to values larger than the actual line length, or you will get two lines instead of one.
  !****************************************************************************
  subroutine bzReadline (f, buf, length, debug)
    type(bzFile) :: f
    character(kind=C_CHAR,len=*) :: buf
    integer :: length
    logical, optional :: debug

    integer(c_int) :: i,r,err
    character(kind=C_CHAR) :: c

    if (length>0) then
      length = length+1   ! add '\n'
      r = bzRead(err,f%bzHandle,buf,length)
      if (err == 4) then
        f%eof = .true.
        if (present(debug) .and. debug) write(*,*) "bzRead(1): EOF"
      else if (err/=0 .or. r<1) then
        write(*,*) "bzRead(1): Error ",err," in file ",trim(f%fname)
        if (err<0) stop
      end if
      if (buf(length:length)==c_new_line .or. f%eof) then
        length=length-1
        if (present(debug) .and. debug) print "(A,i4,2A)","bzReadLine(",length,"):",buf(1:length)
        return
      end if
    end if

    do i=length+1,len(buf)
      r = bzRead(err,f%bzHandle,c,1)
      if (err == 4) then
        f%eof = .true.
        if (present(debug) .and. debug) write(*,*) "bzRead(2): EOF"
      else if (err/=0 .or. r<1) then
        write(*,*) "bzRead(2): Error ",err," in file ",trim(f%fname)
        if (err<0) stop
      end if
      buf(i:i) = c
      if (c==c_new_line .or. f%eof) then
        length = i-1
        if (present(debug) .and. debug) print "(A,i4,2A)","bzReadLine(",length,"):",buf(1:length)
        return
      end if
    end do

    if (c/=c_new_line) then
      print '(3A,I4,2A)', "Error in file ",trim(f%fname),"! bzReadLine: buffer too small (",len(buf),")! Data: ",trim(buf)
      stop
    end if

  end subroutine


  !****************************************************************************
  !****s* bzip/bzWriteline
  ! NAME
  ! subroutine bzWriteline (f, buf)
  ! PURPOSE
  ! Write one line to a bz2 file.
  ! INPUTS
  ! * type(bzFile) :: f --- bzFile structure (which should have been obtained from a previous bzOpenW)
  ! * character(kind=C_CHAR,len=*) :: buf --- buffer containing the data which should be written the file
  !****************************************************************************
  subroutine bzWriteline (f, buf)
    type(bzFile) :: f
    character(kind=C_CHAR,len=*) :: buf

    integer(c_int) :: err

    buf = trim(buf) // c_new_line
    call bzWrite(err,f%bzHandle,buf,len(trim(buf)))
    if (err/=0) then
      write(*,*) "bzWrite: Error ",err
      stop
    end if

  end subroutine


  !****************************************************************************
  !****s* bzip/bzCloseR
  ! NAME
  ! subroutine bzCloseR (f)
  ! PURPOSE
  ! Close a bz2 file, which has been opened for reading.
  ! INPUTS
  ! type(bzFile) :: f --- bzFile structure representing the file which should be closed
  !****************************************************************************
  subroutine bzCloseR (f)
    type(bzFile) :: f

    integer(c_int) :: err

    call bzReadClose(err,f%bzHandle)
    if (err/=0) then
      write(*,*) "Error closing bz file ",f%fname
      write(*,*) "Error code: ",err
      stop
   end if

    err = fClose(f%handle)
    if (err /= 0) then
      write(*,*) "Error closing file ",f%fname
      write(*,*) "Error code: ",err
      stop
    end if

  end subroutine


  !****************************************************************************
  !****s* bzip/bzCloseW
  ! NAME
  ! subroutine bzCloseW (f)
  ! PURPOSE
  ! Close a bz2 file, which has been opened for writing.
  ! INPUTS
  ! type(bzFile) :: f --- bzFile structure representing the file which should be closed
  !****************************************************************************
  subroutine bzCloseW (f)
    type(bzFile) :: f

    integer(c_int) :: err,nbytes_in,nbytes_out

    call bzWriteClose(err,f%bzHandle,0,nbytes_in,nbytes_out)
    if (err/=0) then
      write(*,*) "Error closing bz file ",f%fname
      write(*,*) "Error code: ",err
      stop
    else
      write(*,*) "bzCloseW: writing file ",f%fname
      write(*,*) nbytes_in," Bytes raw"
      write(*,*) nbytes_out," Bytes compressed"
    end if

    err = fClose(f%handle)
    if (err /= 0) then
      write(*,*) "Error closing file ",f%fname
      write(*,*) "Error code: ",err
      stop
    end if

  end subroutine


end module
