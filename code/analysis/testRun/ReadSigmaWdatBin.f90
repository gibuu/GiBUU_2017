program ReadSigmaWdatBin

  use histf90

  type(histogram) :: h

  call Read("1.6")
  call Read("1.7")
  call Read("1.8")
  call Read("1.85")
  call Read("1.9")
  call Read("1.95")
  call Read("1.975")
  call Read("2.0")
  call Read("2.025")
  call Read("2.05")
  call Read("2.1")
  call Read("2.15")
  call Read("2.2")
  call Read("2.3")
  call Read("2.4")
  call Read("2.5")
  call Read("2.6")
  call Read("2.7")
  call Read("2.8")
  call Read("2.9")
  call Read("3.0")



contains

  subroutine Read(bb)

    character*(*) :: bb

    real,dimension(-1:4) :: xx
    real,dimension(1:3) :: S
    real :: mul 

    xx = 0
    S = 0
    mul = 1

    call FetchHist(h,"HiLep.JLABrho0sigmaW.dat.bin."//bb,mul=mul)
    S = SUM(H%yVal,dim=1) - H%yVal(-1,1:3) - H%yVal(0,1:3)
    xx(-1) = S(1)*H%xBin

    call FetchHist(h,"HiLep.JLABrho0sigmaW_proc000.dat.bin."//bb)
    S = SUM(H%yVal,dim=1) - H%yVal(-1,1:3) - H%yVal(0,1:3)
    xx(0) = S(1)*H%xBin

    call FetchHist(h,"HiLep.JLABrho0sigmaW_proc001.dat.bin."//bb)
    S = SUM(H%yVal,dim=1) - H%yVal(-1,1:3) - H%yVal(0,1:3)
    xx(1) = S(1)*H%xBin

    call FetchHist(h,"HiLep.JLABrho0sigmaW_proc002.dat.bin."//bb)
    S = SUM(H%yVal,dim=1) - H%yVal(-1,1:3) - H%yVal(0,1:3)
    xx(2) = S(1)*H%xBin

    call FetchHist(h,"HiLep.JLABrho0sigmaW_proc003.dat.bin."//bb)
    S = SUM(H%yVal,dim=1) - H%yVal(-1,1:3) - H%yVal(0,1:3)
    xx(3) = S(1)*H%xBin

    call FetchHist(h,"HiLep.JLABrho0sigmaW_proc004.dat.bin."//bb)
    S = SUM(H%yVal,dim=1) - H%yVal(-1,1:3) - H%yVal(0,1:3)
    xx(4) = S(1)*H%xBin

    write(*,'(A7,1P,7e13.5)') bb,mul,xx

  end subroutine Read
    

end program ReadSigmaWdatBin
