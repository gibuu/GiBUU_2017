! This is a test case for histogram summation.
!
! It will read histogramMC_avg files from several directories, indicated by 'dirName' and a number,
! and perform a (weighted) summation/average of all the histograms.
!
! Typcial application: Summing histogram files from several heavy-ion runs (possibly with different impact parameters etc).

program sum_HistMC_avg

  use constants, only: pi
  use histMC_avg, only: histogramMC_avg, ReadHistMC_avg, WriteHistMC_avg, sumHistMC_avg

  implicit none

  type(histogramMC_avg) :: sum_hist, tmp_hist

  integer, parameter :: NumFiles = 30
  integer, parameter :: strLen = 100
  character(len=strLen), dimension(1:NumFiles) :: fileNames = ""  ! histogram files to be read
  character(len=strLen) :: dirName = "run_"                       ! directory base name
  integer, save :: nDirs = 96                                     ! number of directories
  integer, save :: weighting = 0                                  ! 0=equal weight, 1=weights proportional to dir number (=impact parameter?)
  real, save :: bmax = -1.                                        ! maximum impact parameter

  integer :: ios, i, i0, j
  character(len=2*strLen) :: f
  real :: w,s,xs

  call init

  print *,"dir name  = ", trim(dirName)
  print *,"num       = ", nDirs
  print *,"weighting = ", weighting
  print *,"bmax      = ", bmax

  if (bmax>0.) then
    ! total XS in mb
    xs = pi * bmax**2 * 10.
    print *,"total XS  = ",xs
  else
    xs = 1.
  end if

  print *,""

  if (weighting == 0) then
    ! equal weighting
    w = 1.
  else
    w = 0.
  end if

  do j=1,NumFiles

    if (trim(fileNames(j)) == "") cycle
    print *,"file name: ", trim(fileNames(j))

    ! initialize summed histogram by reading first summand
    do i0=0,nDirs
      write(f,'(A,i3.3,A)') trim(dirName),i0,"/"//trim(fileNames(j))
      if (weighting == 1) w = real(i0)
      ! read histogram file
      if (ReadHistMC_avg (sum_hist, f, w)) then
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
      if (ReadHistMC_avg (tmp_hist, f)) then
         s = s + w
         ! add histogram
         call sumHistMC_avg(sum_hist, tmp_hist, w)
         write(*,'(A)',advance="NO") "."
      else
         write(*,'(A)',advance="NO") "X"
      end if
    end do

    ! write rebinned histogram
    call WriteHistMC_avg (sum_hist, file=trim(fileNames(j)))  !, mul=xs/s)

    print "(A,f6.1)",new_line("")//"weight sum: ",s
    print *,""

  end do

contains

  subroutine init
    use version, only: printVersion
    use output, only: Write_ReadingInput

    namelist /sumHist/ fileNames, dirName, nDirs, weighting, bmax

    call printVersion

    ! read namelist input
    call Write_ReadingInput('sumHist',0)
    rewind(5)
    read(5,nml=sumHist,iostat=ios)
    call Write_ReadingInput('sumHist',0,ios)
    call Write_ReadingInput('sumHist',1)

  end subroutine

end program
