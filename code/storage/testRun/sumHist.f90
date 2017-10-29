! This is a test case for histogram summation.
!
! It will read histogram files from several directories, indicated by 'dirName' and a number,
! and perform a (weighted) summation/average of all the histograms.
!
! Typcial application: Summing histogram files from several runs.

program sum_Hist

  use histf90, only: histogram, ReadHist, WriteHist, sumHist

  implicit none

  type(histogram) :: hist, tmp_hist

  integer, save :: nFiles = 10                                ! number of file names to be processed
  character(len=100), allocatable, dimension(:) :: fileNames  ! histogram files to be read
  character(len=100), save :: pattern = ""                    ! file-name pattern (optional), given as format specifier
  character(len=100) :: dirName = "run_"                 ! directory base name
  integer, save :: nDirs = 96                            ! number of directories
  integer, save :: weighting = 0                         ! 0=equal weight, 1=weights proportional to dir number (=impact parameter?)

  integer :: ios, i, i0, j
  character(len=100) :: f
  real :: w,s

  call init

  print *,"nFiles = ", nFiles
  print *,"nDirs  = ", nDirs
  print *,"dir name: ", trim(dirName)
  print *,"weighting = ", weighting
  print *,"pattern   = ", trim(pattern)
  print *,""

  do j=1,nFiles

    ! if a pattern is given, produce a fileName from it (using the running file number)
    if (trim(pattern) /= "") write(fileNames(j),pattern) j

    if (trim(fileNames(j)) == "") cycle
    print *,"file name: ", trim(fileNames(j))

    if (weighting == 0) then
      w = 1.
    else
      w = 0.
    end if

    ! initialize summed histogram by reading first summand
    do i0=0,nDirs
      write(f,'(A,i3.3,A)') trim(dirName),i0,"/"//trim(fileNames(j))
      if (weighting == 1) w = real(i0)
      ! read histogram file
      if (ReadHist (hist, f)) then
        s = w
        write(*,'(A)',advance="NO") "_"
        exit
      else
        s = 0.
        write(*,'(A)',advance="NO") "I"
      end if
    end do

    ! sum up remaining histograms
    do i=i0+1,nDirs
      if (i>999) then
         write(f,'(A,i4.4,A)') trim(dirName),i,"/"//trim(fileNames(j))
      else
         write(f,'(A,i3.3,A)') trim(dirName),i,"/"//trim(fileNames(j))
      end if
      if (weighting == 1) w = real(i)
      ! read histogram file
      if (ReadHist (tmp_hist, f)) then
        s = s + w
        ! add histogram
        call sumHist (hist, tmp_hist, w)
        write(*,'(A)',advance="NO") "."
      else
        write(*,'(A)',advance="NO") "X"
      end if
    end do

    ! write rebinned histogram
    call WriteHist (hist, file=trim(fileNames(j)), mul=1./s)

    print "(A,f5.1)",new_line("")//"weight sum: ",s
    print *,""

  end do

contains

  subroutine init
    use version, only: printVersion
    use output, only: Write_ReadingInput

    namelist /sumHist/ nFiles, nDirs, dirName, weighting
    namelist /files/ fileNames, pattern

    call printVersion

    ! read namelist input
    call Write_ReadingInput('sumHist',0)
    rewind(5)
    read(5,nml=sumHist,iostat=ios)
    call Write_ReadingInput('sumHist',0,ios)
    call Write_ReadingInput('sumHist',1)

    ! allocate fileNames array
    allocate(fileNames(1:nFiles))
    fileNames = ""

    ! read namelist input
    call Write_ReadingInput('files',0)
    rewind(5)
    read(5,nml=files,iostat=ios)
    call Write_ReadingInput('files',0,ios)
    call Write_ReadingInput('files',1)


  end subroutine

end program
