program HiLeptonCLASrho0

  use histf90
  use hist2Df90

  implicit none

  type(histogram)   :: hLep_Q2
  type(histogram)   :: H
  type(histogram2D) :: H2D
  type(histogram)   :: Htry

  ! Fe, 5GeV:

  real, dimension(-1:6), parameter :: cutsQ2 = (/0.,0.8,1.1,1.3,1.5,2.0,3.0,99.0/)
  real, dimension(0:6), parameter :: numin  = (/0.0,2.15,2.25,2.35,2.45,2.70, 0.0/)
  real, dimension(0:6), parameter :: numax  = (/0.0,4.85,4.80,4.75,4.70,4.45, 0.0/)

  ! C, 5 GeV:

!!$  real, dimension(0:6), parameter :: cutsQ2 = (/0.8,1.15,1.4,3.96,99.0,99.0,99.0/)
!!$  real, dimension(0:6), parameter :: numin  = (/0.0,2.15,2.25,2.40,0.0, 0.0, 0.0/)
!!$  real, dimension(0:6), parameter :: numax  = (/0.0,4.80,4.80,4.75,0.0, 0.0, 0.0/)

  ! Fe, 4GeV:

!!$  real, dimension(0:6), parameter :: cutsQ2 = (/0.6,0.875,1.075,2.4,99.0,99.0,99.0/)
!!$  real, dimension(0:6), parameter :: numin  = (/0.0,2.05,2.10,2.20,0.0, 0.0, 0.0/)
!!$  real, dimension(0:6), parameter :: numax  = (/0.0,3.85,3.80,3.75,0.0, 0.0, 0.0/)


  real, dimension(0:6,0:1) :: ArrQ2,Arrnu,ArrY

  integer :: i1,i2,ibin
  integer :: iMax1, iMax2

  real :: nu,Q2,Y,nLeptons,HH(2),AA(2)
  real :: add,mul


  call DoIt2

  contains

    subroutine DoIt2

      call FetchHist(H,"HiLep.Rho0_NuQ2.Int2.dat.bin",add=add,mul=mul)

      iMax1=ubound(H%yVal,dim=1)
      ArrY  = 0
      ArrQ2 = 0
      ibin = 0
      do i1=1,iMax1
         Q2 = H%xMin+(real(i1)-0.5)*H%xBin
         if (Q2.gt.cutsQ2(ibin)) then
            ibin=ibin+1
         end if
!         Y = H%yVal(i1,1)
         Y = H%yVal(i1,3)

         ArrY(ibin,:)=ArrY(ibin,:)+ (/1.0,Y /)
         ArrQ2(ibin,:)=ArrQ2(ibin,:)+ (/Y,Y*Q2 /)

      end do

      ArrQ2(0,0)=ArrQ2(0,0)+1.0 ! DUMMY

      open(121,file="HiLep.Rho0_NuQ2.Int2.Rebin.dat",status="unknown")

      do i1 = 0,6
         if (ArrY(i1,0).gt.0) then
            Q2 = ArrQ2(i1,1)/ArrQ2(i1,0)
            write(121,'(1P,9e14.3)')  Q2,cutsQ2(i1-1),cutsQ2(i1),&
                 & ArrY(i1,1)/ArrY(i1,0)*mul
         end if
      end do

      close(121)

    end subroutine DoIt2


    subroutine DoIt1


      call FetchHist(hLep_Q2,"HiLep.lep.Q2.kinematics.dat.bin")
      nLeptons = SUM(hLep_Q2%yVal(:,1),dim=1)

      call FetchHist2D(H2D,"HiLep.JLABrhoNuQ2.dat.bin")
      !  call FetchHist2D(H2D,"HiLep.JLABrhoNuQ2_proc001.dat.bin")

      call IntegrateHist2D(H2D,Htry,2)
      call WriteHist(Htry,231,mul=1.0/nLeptons)

      call AverageHist2D(H2D,Htry,2)
      call WriteHist(Htry,531)

      ArrQ2 = 0
      Arrnu = 0
      ArrY  = 0

      ibin=0

      iMax1 = uBound(H2D%yVal,dim=1)
      iMax2 = uBound(H2D%yVal,dim=2)

      do i2=1,iMax2

         AA = 0

         do i1=1,iMax1

            !        (/nu,Q2/) = H2D%xMin+((/i1,i2/)*1.0-0.5)*H2D%xBin
            HH=H2D%xMin+((/i1,i2/)*1.0-0.5)*H2D%xBin
            nu=HH(1)
            Q2=HH(2)

            !        write(*,*) ibin,nu,Q2

            if (Q2.gt.cutsQ2(ibin)) then
               ibin=ibin+1
            end if

!!!        if (nu.lt.0.7+Q2/0.938) cycle
            !        if (nu.lt.1.0+Q2/0.938) cycle

            Y = H2D%yVal(i1,i2,1)/nLeptons

            AA = AA + (/1.0,Y /)

            if (nu.lt.numin(ibin)) cycle
            if (nu.gt.numax(ibin)) cycle

            ArrQ2(ibin,:) = ArrQ2(ibin,:)+ (/Y,Y*Q2 /)
            Arrnu(ibin,:) = Arrnu(ibin,:)+ (/Y,Y*nu /)
            ArrY(ibin,:) = ArrY(ibin,:)+ (/1.0,Y /)



         end do
         !     write(*,*) 

         !     if (AA(1).gt.0) AA(2)=AA(2)/AA(1)
         write(*,'(1P,9e14.3)') Q2,AA(2)
         write(31,'(1P,9e14.3)') Q2,AA(2)

      end do
      do i1 = 0,6
         if (ArrY(i1,0).gt.0) then
            !        write(*,'(1P,9e14.3)')  ArrQ2(i1,0),ArrQ2(i1,1),ArrQ2(i1,1)/ArrQ2(i1,0),&
            !             & Arrnu(i1,0),Arrnu(i1,1),Arrnu(i1,1)/Arrnu(i1,0), &
            !             & ArrY(i1,0),ArrY(i1,1),ArrY(i1,1)/ArrY(i1,0)
            Q2 = ArrQ2(i1,1)/ArrQ2(i1,0)
            write(*,'(1P,9e14.3)')  Q2,cutsQ2(i1-1),cutsQ2(i1),&
                 & Arrnu(i1,1)/Arrnu(i1,0), &
                 & ArrY(i1,1)/ArrY(i1,0)
            write(32,'(1P,9e14.3)')  Q2,cutsQ2(i1-1),cutsQ2(i1),&
                 & Arrnu(i1,1)/Arrnu(i1,0), &
                 & ArrY(i1,1)/ArrY(i1,0)
         end if
      end do


    end subroutine DoIt1


  

end program HiLeptonCLASrho0
