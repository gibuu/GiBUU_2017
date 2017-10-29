program main

  use histf90
  use hist2Df90

  type(histogram),save :: hist1
  type(histogram),save :: histPT_pi0,histPT_pip,histPT_pim
  type(histogram2D),save :: H2D_AllPion,H2D_AllOther
  type(histogram2D),save :: H2D_pi0,H2D_pip,H2D_pim
  



!  integer, parameter :: nDir = 10
!  integer, parameter :: nPlus = 1000
!  character*(*), dimension(nDir), parameter :: sDir = (/&
!       '00p_00/run', '00p_01/run', '00p_02/run', '00p_03/run', '00p_04/run', &
!       '00p_05/run', '00p_06/run', '00p_07/run', '00p_08/run', '00p_09/run' /)

!!$  integer, parameter :: nDir = 10
!!$  integer, parameter :: nPlus = 2000
!!$  character*(*), dimension(nDir), parameter :: sDir = (/&
!!$       'Cu_00/runBAK', 'Cu_01/runBAK', 'Cu_02/runBAK', 'Cu_03/runBAK', 'Cu_04/runBAK', &
!!$       'Cu_05/runBAK', 'Cu_06/runBAK', 'Cu_07/runBAK', 'Cu_08/runBAK', 'Cu_09/runBAK' /)

  integer, parameter :: nDir = 10
  integer, parameter :: nPlus = 3000
  character*(*), dimension(nDir), parameter :: sDir = (/&
       'Be_00/run', 'Be_01/run', 'Be_02/run', 'Be_03/run', 'Be_04/run', &
       'Be_05/run', 'Be_06/run', 'Be_07/run', 'Be_08/run', 'Be_09/run' /)


  call CreateHist(hist1,'p_Z, all',-10.,510.,5.)
  call CreateHist(histPT_pi0,'pT, pi0',0.,10.,0.2)
  call CreateHist(histPT_pip,'pT, pip',0.,10.,0.2)
  call CreateHist(histPT_pim,'pT, pim',0.,10.,0.2)


  call CreateHist2D(H2D_AllPion,  'H2D_AllPion',  (/-10.,0./), (/10.,10./), (/0.2,0.2/)) ! y_L, p_T
  call CreateHist2D(H2D_AllOther, 'H2D_AllOther', (/-10.,0./), (/10.,10./), (/0.2,0.2/)) ! y_L, p_T

  call ReadMultipleHist(nDir,sDir,hist1,100)
  call WriteHist(hist1,100 + nPlus)

  call ReadMultipleHist(nDir,sDir,histPT_pip,101)
  call WriteHist(histPT_pip,101 + nPlus)
  call ReadMultipleHist(nDir,sDir,histPT_pi0,102)
  call WriteHist(histPT_pi0,102 + nPlus)
  call ReadMultipleHist(nDir,sDir,histPT_pim,103)
  call WriteHist(histPT_pim,103 + nPlus)

  
  call ReadMultipleHist2D(nDir,sDir,H2D_AllPion,310)
  call WriteHist2D_Gnuplot(H2D_AllPion,310 + nPlus)

  call ReadMultipleHist2D(nDir,sDir,H2D_AllOther,311)
  call WriteHist2D_Gnuplot(H2D_AllOther,311 + nPlus)

contains

  subroutine ReadMultipleHist(nDir,sDir, H, iFile)
    implicit none
    integer :: nDir
    character*(*), dimension(nDir) :: sDir
    type(histogram) :: H
    integer :: iFile


    integer :: iD,iL, iBin, iBinMax
    character*(80) :: text
    character*(200) :: BUF
    real, dimension(4) :: val

    iBinMax = ubound(H%yVal,dim=1)

    do iD=1,nDir
       write(text,'(A,"/fort.",i3.0)') sDir(iD),iFile
!       text = sDir(i)//"/fort."
       write(*,'("Reading: ",A,"...")') text
       
       open(13,file=text,status='unknown')

       do iL=1,9
          read(13,'(A)') BUF
       end do



       do iBin=1,iBinMax
          read(13,*) val

          H%yVal(iBin,1) = H%yVal(iBin,1) + val(2)*H%xBin
          H%yVal(iBin,2) = H%yVal(iBin,2) + val(3)
          H%yVal(iBin,3) = H%yVal(iBin,3) + val(4)*H%xBin
       end do
          
       close (13)
    end do


    H%yVal(:,1) = H%yVal(:,1) / dble(nDir)
    H%yVal(:,3) = H%yVal(:,3) / dble(nDir)

    
  end subroutine ReadMultipleHist

  subroutine ReadMultipleHist2D(nDir,sDir, H, iFile)
    implicit none
    integer :: nDir
    character*(*), dimension(nDir) :: sDir
    type(histogram2D) :: H
    integer :: iFile

    integer :: iD,iL
    integer, dimension(2) :: iBinMax
    integer :: iBin1,iBin2

    character*(80) :: text
    character*(200) :: BUF
    real, dimension(5) :: val

    iBinMax(1) = uBound(H%yVal,dim=1)
    iBinMax(2) = uBound(H%yVal,dim=2)


    do iD=1,nDir
       write(text,'(A,"/fort.",i3.0)') sDir(iD),iFile
!       text = sDir(i)//"/fort."
       write(*,'("Reading: ",A,"...")') text
       
       open(13,file=text,status='unknown')

       do iBin1=1,iBinMax(1)
          do iBin2=1,iBinMax(2)
             read(13,*) val
             H%yVal(iBin1,iBin2,1:3) = H%yVal(iBin1,iBin2,1:3)+val(3:5)*(H%xBin(1)*H%xBin(2))
          end do
       end do

       close (13)
    end do

    H%yVal = H%yVal / dble(nDir)

  end subroutine ReadMultipleHist2D

end program main
