module Read
#define BB 7
#define SS 20
#define II 20
#define NN 1000
  integer, parameter :: qBmax = BB
  integer, parameter :: qSmax = SS
  integer, parameter :: qIx2max = II
  integer, parameter :: nBin=NN
  integer, parameter :: nBin1=300
  real, parameter :: dM = 0.01

  double precision, dimension(1:nBin1,  1: 61), save :: WarrB
  double precision, dimension(1:nBin1,101:122), save :: WarrM
  double precision, dimension(1:nBin1,-1:1,-3:3,0:3), save :: ArrTauHadron

  double precision, dimension(1:NN,-BB:BB,-SS:SS,0:II), save :: ArrTau
  double precision, dimension(1:NN,-BB:BB,-SS:SS,0:II), save :: ArrGamma
  integer, dimension(-BB:BB,-SS:SS,0:II), save :: ArrMinIM

contains

  subroutine readArrays

    integer :: ios
    integer :: iBmax, iSmax, iIx2max, iBin, iBin1
    integer :: qB, qS, qIx2, iM, iH1, iH2, iiH

    double precision, dimension(1:nBin1) :: hArr1
    double precision, dimension(1:NN) :: hArrTau, hArrGamma
    integer :: hArrMinIm

    write(*,*) "Reading the file 'Bootstrap.dat': ..."
    open(131,file='Bootstrap.dat',status='old', &
         & form='UNFORMATTED',iostat=ios,position='rewind')
    if (ios/=0) then
       write(*,*) "Error while opening the file."
       stop
    endif

    read(131) iBmax, iSmax, iIx2max, iBin, iBin1
    write(*,*) iBmax, iSmax, iIx2max, iBin, iBin1
    if (iBmax /= BB) then
       write(*,*) "iBmax /= qBmax"
       stop
    end if
    if (iSmax /= SS) then
       write(*,*) "iSmax /= qSmax"
       stop
    end if
    if (iIx2max /= II) then
       write(*,*) "iIx2max /= qIx2max"
       stop
    end if
    if (iBin /= NN) then
       write(*,*) "iBin /= nBin"
       stop
    end if
    if (iBin1 /= 300) then
       write(*,*) "iBin1 /= nBin1"
       stop
    end if

    read(131,iostat=ios) iM
    if (ios/=0) then
       write(*,*) "Error while reading the file: iM"
       stop
    end if
    read(131,iostat=ios) WarrB
    if (ios/=0) then
       write(*,*) "Error while reading the file: WarrB"
       stop
    end if
    read(131,iostat=ios) WarrM
    if (ios/=0) then
       write(*,*) "Error while reading the file: WarrM"
       stop
    end if
    read(131,iostat=ios) ArrTauHadron
    if (ios/=0) then
       write(*,*) "Error while reading the file: ArrTauHadron"
       stop
    end if
    read(131,iostat=ios) ArrTau
    if (ios/=0) then
       write(*,*) "Error while reading the file: ArrTau"
       stop
    end if
    read(131,iostat=ios) ArrGamma
    if (ios/=0) then
       write(*,*) "Error while reading the file: ArrGamma"
       stop
    end if
    read(131,iostat=ios) ArrMinIM
    if (ios/=0) then
       write(*,*) "Error while reading the file: ArrMinIM"
       stop
    end if

    close(131)

    write(*,*) "Reshuffling the entries"

    ! Since Python arrays always start with index 0, but have the
    ! ability to handel negative indizes as give from the end,
    ! we are shuffling all negative index entries to the end.
    ! This has to be done in B and in S direction seperately

    do qS=-3,3
       do qIx2=0,3
          hArr1 = ArrTauHadron(:,-1,qS,qIx2)
          ArrTauHadron(:,-1,qS,qIx2) = ArrTauHadron(:,0,qS,qIx2)
          ArrTauHadron(:, 0,qS,qIx2) = ArrTauHadron(:,1,qS,qIx2)
          ArrTauHadron(:, 1,qS,qIx2) = hArr1
       end do
    end do

!    -3 -2 -1  0  1  2  3
!
!    xx -2 -1  0  1  2  3
!     0 -2 -1 xx  1  2  3
!     0 -2 -1  3  1  2 xx
!     0 -2 xx  3  1  2 -1
!     0 -2  2  3  1 xx -1
!     0 xx  2  3  1 -2 -1
!     0  1  2  3 xx -2 -1
!
!     0  1  2  3 -3 -2 -1

    do qB=-1,1
       do qIx2=0,3
          iH1 = -3
          hArr1 = ArrTauHadron(:,qB,iH1,qIx2)
          do iiH=1,6
             iH2 = iH1+3
             if (iH2>3) iH2=iH2-7
             ArrTauHadron(:,qB,iH1,qIx2) = ArrTauHadron(:,qB,iH2,qIx2)
             iH1 = iH2
          end do
          ArrTauHadron(:,qB,iH1,qIx2) = hArr1
       end do
    end do

    do qS=-SS,SS
       do qIx2=0,II
          iH1 = -BB
          hArrTau = ArrTau(:,iH1,qS,qIx2)
          hArrGamma = ArrGamma(:,iH1,qS,qIx2)
          hArrMinIM = ArrMinIM(iH1,qS,qIx2)
          do iiH=1,2*BB
             iH2 = iH1+BB
             if (iH2>BB) iH2=iH2-(2*BB+1)
             ArrTau(:,iH1,qS,qIx2) = ArrTau(:,iH2,qS,qIx2)
             ArrGamma(:,iH1,qS,qIx2) = ArrGamma(:,iH2,qS,qIx2)
             ArrMinIM(iH1,qS,qIx2) = ArrMinIM(iH2,qS,qIx2)
             iH1 = iH2
          end do
          ArrTau(:,iH1,qS,qIx2) = hArrTau
          ArrGamma(:,iH1,qS,qIx2) = hArrGamma
          ArrMinIM(iH1,qS,qIx2) = hArrMinIM
       end do
    end do

    do qB=-BB,BB
       do qIx2=0,II
          iH1 = -SS
          hArrTau = ArrTau(:,qB,iH1,qIx2)
          hArrGamma = ArrGamma(:,qB,iH1,qIx2)
          hArrMinIM = ArrMinIM(qB,iH1,qIx2)
          do iiH=1,2*SS
             iH2 = iH1+SS
             if (iH2>SS) iH2=iH2-(2*SS+1)
             ArrTau(:,qB,iH1,qIx2) = ArrTau(:,qB,iH2,qIx2)
             ArrGamma(:,qB,iH1,qIx2) = ArrGamma(:,qB,iH2,qIx2)
             ArrMinIM(qB,iH1,qIx2) = ArrMinIM(qB,iH2,qIx2)
             iH1 = iH2
          end do
          ArrTau(:,qB,iH1,qIx2) = hArrTau
          ArrGamma(:,qB,iH1,qIx2) = hArrGamma
          ArrMinIM(qB,iH1,qIx2) = hArrMinIM
       end do
    end do

    write(*,*) "Reading the file 'Bootstrap.dat': Finished"

  end subroutine readArrays


end module Read
