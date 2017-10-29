program HiLeptonMix10Runs

  use histf90
  use hist2Df90
  use histMPf90
  use histMC
  use output

  implicit none

  type(histogram)   :: HHH,hLep_Q2,hLep_nu
  type(histogramMP) :: hMP
  type(histogram2D) :: h2D
  type(histogramMC) :: hMC
  type(histogram)   :: HHH1,HHH2,HHH3


  real, save :: nLeptons
  character(3), dimension(-1:1), parameter :: piName  = (/'pi-','pi0','pi+'/)
  character(2), dimension(0:1),  parameter :: nucName = (/'N0','N+'/)
  integer :: i
  logical :: lDum

  ! ===== Lepton Kinematics =====


  call FetchHist_10(hLep_Q2,"HiLep.lep.Q2.kinematics.dat")
  call WriteHist(hLep_Q2,101,add=1e-20,file="HiLep.lep.Q2.kinematics.dat",dump=.true.)
  call WriteHist(hLep_Q2,101,add=1e-20,DoAve=.true.,&
       &file="HiLep.lep.nuAve_Q2.kinematics.dat")

  nLeptons = SUM(hLep_Q2%yVal(:,1),dim=1)
  write(*,*) nLeptons

  call FetchHist_10(hLep_nu,"HiLep.lep.nu.kinematics.dat")
  call WriteHist(hLep_nu,101,add=1e-20,file="HiLep.lep.nu.kinematics.dat")
  call WriteHist(hLep_nu,101,add=1e-20,DoAve=.true.,&
       &file="HiLep.lep.Q2Ave_nu.kinematics.dat")

  ! to be continued...


  ! ===== nu-spectra =====

  call FetchHistMP_10(hMP,"HiLep.nu.idH.Acc.dat")
  call WriteHistMP(hMP,101,add=1e-20,H2=hLep_nu,iColumn=3,&
       &file='HiLep.nu.idH.noAcc.dat') ! noAcc
  call WriteHistMP(hMP,102,add=1e-20,H2=hLep_nu,iColumn=1,&
       &file='HiLep.nu.idH.Acc.dat') ! Acc

  call FetchHistMP_10(hMP,"HiLep.AvePT2.nu.dat")
  call WriteHistMP(hMP,101,DoAve=.true.,&
       &file='HiLep.AvePT2.nu.dat')


  ! to be continued...


  ! ===== Q2-spectra =====

  call FetchHistMP_10(hMP,"HiLep.Q2.idH.Acc.dat")
  call WriteHistMP(hMP,101,add=1e-20,H2=hLep_Q2,iColumn=3,&
       &file='HiLep.Q2.idH.noAcc.dat') ! noAcc
  call WriteHistMP(hMP,102,add=1e-20,H2=hLep_Q2,iColumn=1,&
       &file='HiLep.Q2.idH.Acc.dat') ! Acc

  call FetchHistMP_10(hMP,"HiLep.AvePT2.Q2.dat")
  call WriteHistMP(hMP,101,DoAve=.true.,&
       &file='HiLep.AvePT2.Q2.dat')


  ! to be continued...


  ! ===== zH-spectra =====

  call FetchHistMP_10(hMP,"HiLep.zH.idH.Acc.dat")
  call WriteHistMP(hMP,101,add=1e-20,mul=1./NLeptons,iColumn=3,&
       &file='HiLep.zH.idH.noAcc.dat') ! noAcc
  call WriteHistMP(hMP,102,add=1e-20,mul=1./NLeptons,iColumn=1,&
       &file='HiLep.zH.idH.Acc.dat') ! Acc

  call FetchHistMP_10(hMP,"HiLep.AvePT2.zH.dat")
  call WriteHistMP(hMP,101,DoAve=.true.,&
       &file='HiLep.AvePT2.zH.dat')


  ! to be continued...

  ! ===== pT2-spectra =====

  call FetchHistMP_10(hMP,"HiLep.pT2.idH.Acc.dat")
  call WriteHistMP(hMP,101,add=1e-20,mul=1./NLeptons,iColumn=3,&
       &file='HiLep.pT2.idH.noAcc.dat') ! noAcc
  call WriteHistMP(hMP,102,add=1e-20,mul=1./NLeptons,iColumn=1,&
       &file='HiLep.pT2.idH.Acc.dat') ! Acc

  ! =====================

  do i=-1,1
     call FetchHist2D_10(h2D,'HiLep.pT2Ave_Plane.'//piName(i)//'.dat')
     call WriteHist2D_Gnuplot(h2D,140,DoAve=.true.,maxval=0.0,&
          &file='HiLep.pT2Ave_Plane.'//piName(i)//'.dat')

     call FetchHist2D_10(h2D,'HiLep.pT2_zH.'//piName(i)//'.dat')
     call WriteHist2D_Gnuplot(h2D,140, mul=1./nLeptons, add=1e-20,&
          &file='HiLep.pT2_zH.'//piName(i)//'.dat')

     call FetchHist2D_10(h2D,'HiLep.pT2_nu.'//piName(i)//'.dat')
     call WriteHist2D_Gnuplot(h2D,140, H2=hLep_nu, add=1e-20,&
          &file='HiLep.pT2_nu.'//piName(i)//'.dat')

     call FetchHist2D_10(h2D,'HiLep.pT2_Q2.'//piName(i)//'.dat')
     call WriteHist2D_Gnuplot(h2D,140, H2=hLep_Q2, add=1e-20,&
          &file='HiLep.pT2_Q2.'//piName(i)//'.dat')
  end do


  ! =====================
  ! Binning_11:
  ! =====================

  if (TryFile("0/HiLep.dN_id.nu.zH.001.dat.bin")) then
     do i=1,3
        call FetchHistMP_10(hMP,&
             &file='HiLep.dN_id.nu.zH.'//trim(intToChar(i))//'.dat')
        call WriteHistMP(hMP,141,add=1e-20,mul=1./NLeptons,iColumn=1,&
             &file='HiLep.dN_id.nu.zH.'//trim(intToChar(i))//'.dat') ! Acc
        call FetchHistMP_10(hMP,&
             &file='HiLep.dN_id.Q2.zH.'//trim(intToChar(i))//'.dat')
        call WriteHistMP(hMP,141,add=1e-20,mul=1./NLeptons,iColumn=1,&
             &file='HiLep.dN_id.Q2.zH.'//trim(intToChar(i))//'.dat') ! Acc
        call FetchHistMP_10(hMP,&
             &file='HiLep.dN_id.pT2.zH.'//trim(intToChar(i))//'.dat')
        call WriteHistMP(hMP,141,add=1e-20,mul=1./NLeptons,iColumn=1,&
             &file='HiLep.dN_id.pT2.zH.'//trim(intToChar(i))//'.dat') ! Acc

        call FetchHistMP_10(hMP,&
             &file='HiLep.dN_id.zH.nu.'//trim(intToChar(i))//'.dat')
        call WriteHistMP(hMP,141,add=1e-20,mul=1./NLeptons,iColumn=1,&
             &file='HiLep.dN_id.zH.nu.'//trim(intToChar(i))//'.dat') ! Acc
        call FetchHistMP_10(hMP,&
             &file='HiLep.dN_id.Q2.nu.'//trim(intToChar(i))//'.dat')
        call WriteHistMP(hMP,141,add=1e-20,mul=1./NLeptons,iColumn=1,&
             &file='HiLep.dN_id.Q2.nu.'//trim(intToChar(i))//'.dat') ! Acc
        call FetchHistMP_10(hMP,&
             &file='HiLep.dN_id.pT2.nu.'//trim(intToChar(i))//'.dat')
        call WriteHistMP(hMP,141,add=1e-20,mul=1./NLeptons,iColumn=1,&
             &file='HiLep.dN_id.pT2.nu.'//trim(intToChar(i))//'.dat') ! Acc
     end do
     do i=1,2
        call FetchHistMP_10(hMP,&
             &file='HiLep.dN_id.nu.pT2.'//trim(intToChar(i))//'.dat')
        call WriteHistMP(hMP,141,add=1e-20,mul=1./NLeptons,iColumn=1,&
             &file='HiLep.dN_id.nu.pT2.'//trim(intToChar(i))//'.dat') ! Acc
        call FetchHistMP_10(hMP,&
             &file='HiLep.dN_id.zH.pT2.'//trim(intToChar(i))//'.dat')
        call WriteHistMP(hMP,141,add=1e-20,mul=1./NLeptons,iColumn=1,&
             &file='HiLep.dN_id.zH.pT2.'//trim(intToChar(i))//'.dat') ! Acc
        call FetchHistMP_10(hMP,&
             &file='HiLep.dN_id.Q2.pT2.'//trim(intToChar(i))//'.dat')
        call WriteHistMP(hMP,141,add=1e-20,mul=1./NLeptons,iColumn=1,&
             &file='HiLep.dN_id.Q2.pT2.'//trim(intToChar(i))//'.dat') ! Acc
     end do

  endif

!!$  call FetchHist2D_10(h2D,'HiLep.CollHistPT2.1.dat')
!!$  call WriteHist2D_Gnuplot(h2D,140, DoAve=.true., mul=1./nLeptons, add=1e-20,&
!!$       & file='HiLep.CollHistPT2.1.dat') ! Acc
!!$
!!$  call FetchHist2D_10(h2D,'HiLep.CollHistPT2.2.dat')
!!$  call WriteHist2D_Gnuplot(h2D,140, DoAve=.true., mul=1./nLeptons, add=1e-20,&
!!$       & file='HiLep.CollHistPT2.2.dat') ! Acc

  if (TryFile("HiLep.nleadPT.N.000.dat.bin")) then
     do i=0,3
        call FetchHistMP_10(hMP,'HiLep.nleadPT.N.'//trim(intToChar(i))//'.dat')
        call WriteHistMP(hMP, 56, add=1e-20,mul=1./NLeptons,iColumn=1,&
             & file='HiLep.nleadPT.N.'//trim(intToChar(i))//'.dat') ! noAcc
        call WriteHistMP(hMP, 56 ,mul=1./NLeptons,DoAve=.true.,&
             & file='HiLep.nleadPT.PT2.'//trim(intToChar(i))//'.dat') ! noAcc
     end do
  end if

  if (TryFile("0/HiLep.Rho0_NuQ2.dat.bin")) then
     call FetchHist2D_10(h2D,file='HiLep.Rho0_NuQ2.dat')
     call WriteHist2D_Gnuplot(h2D,140, mul=1./nLeptons, add=1e-20,&
          & file='HiLep.Rho0_NuQ2.dat',dump=.true.) ! Acc
     call WriteHist2D_Gnuplot(h2D,140, DoAve=.true., mul=1./nLeptons, add=1e-20,&
          & file='HiLep.Rho0_NuQ2.AveMis.dat',dump=.false.) ! Acc

     call IntegrateHist2D(h2D,HHH,2)
     call WriteHist(HHH,140,add=1e-20,mul=1./nLeptons,&
          & file="HiLep.Rho0_NuQ2.Int2.dat",dump=.true.)
     call AverageHist2D(h2D,HHH,2)
     call WriteHist(HHH,140,&
          & file="HiLep.Rho0_NuQ2.Ave2.dat",dump=.true.)

  endif

  do i=0,4
     if (TryFile("0/HiLep.Rho0_NuQ2_proc"//trim(intToChar(i))//".dat.bin")) then
        call FetchHist2D_10(h2D,file='HiLep.Rho0_NuQ2_proc'//trim(intToChar(i))//'.dat')
        call WriteHist2D_Gnuplot(h2D,140, mul=1./nLeptons, add=1e-20,&
             & file='HiLep.Rho0_NuQ2_proc'//trim(intToChar(i))//'.dat',dump=.true.) ! Acc
        call WriteHist2D_Gnuplot(h2D,140, DoAve=.true., mul=1./nLeptons, add=1e-20,&
             & file='HiLep.Rho0_NuQ2_proc'//trim(intToChar(i))//'.AveMis.dat',dump=.false.) ! Acc

        call IntegrateHist2D(h2D,HHH,2)
        call WriteHist(HHH,140,add=1e-20,mul=1./nLeptons,&
             & file='HiLep.Rho0_NuQ2_proc'//trim(intToChar(i))//'.Int2.dat',dump=.true.) ! Acc
        call AverageHist2D(h2D,HHH,2)
        call WriteHist(HHH,140,&
             & file='HiLep.Rho0_NuQ2_proc'//trim(intToChar(i))//'.Ave2.dat',dump=.true.) ! Acc
     end if
  end do

  if (TryFile("0/HiLep.Rho0_WQ2.dat.bin")) then
     call FetchHist2D_10(h2D,file='HiLep.Rho0_WQ2.dat')
     call WriteHist2D_Gnuplot(h2D,140, mul=1./nLeptons, add=1e-20,&
          & file='HiLep.Rho0_WQ2.dat',dump=.true.) ! Acc
     call WriteHist2D_Gnuplot(h2D,140, DoAve=.true., mul=1./nLeptons, add=1e-20,&
          & file='HiLep.Rho0_WQ2.AveMis.dat',dump=.false.) ! Acc

     call IntegrateHist2D(h2D,HHH,2)
     call WriteHist(HHH,140,add=1e-20,mul=1./nLeptons,&
          & file="HiLep.Rho0_WQ2.Int2.dat",dump=.true.)
     call AverageHist2D(h2D,HHH,2)
     call WriteHist(HHH,140,&
          & file="HiLep.Rho0_WQ2.Ave2.dat",dump=.true.)

  endif

  do i=0,4
     if (TryFile("0/HiLep.Rho0_WQ2_proc"//trim(intToChar(i))//".dat.bin")) then
        call FetchHist2D_10(h2D,file='HiLep.Rho0_WQ2_proc'//trim(intToChar(i))//'.dat')
        call WriteHist2D_Gnuplot(h2D,140, mul=1./nLeptons, add=1e-20,&
             & file='HiLep.Rho0_WQ2_proc'//trim(intToChar(i))//'.dat',dump=.true.) ! Acc
        call WriteHist2D_Gnuplot(h2D,140, DoAve=.true., mul=1./nLeptons, add=1e-20,&
             & file='HiLep.Rho0_WQ2_proc'//trim(intToChar(i))//'.AveMis.dat',dump=.false.) ! Acc

        call IntegrateHist2D(h2D,HHH,2)
        call WriteHist(HHH,140,add=1e-20,mul=1./nLeptons,&
             & file='HiLep.Rho0_WQ2_proc'//trim(intToChar(i))//'.Int2.dat',dump=.true.) ! Acc
        call AverageHist2D(h2D,HHH,2)
        call WriteHist(HHH,140,&
             & file='HiLep.Rho0_WQ2_proc'//trim(intToChar(i))//'.Ave2.dat',dump=.true.) ! Acc
     end if
  end do



  if (TryFile("0/HiLep.Rho0_NuQ2_VMD.dat.bin")) then
     call FetchHist2D_10(h2D,file='HiLep.Rho0_NuQ2_VMD.dat')
     call WriteHist2D_Gnuplot(h2D,140, mul=1./nLeptons, add=1e-20,&
          & file='HiLep.Rho0_NuQ2_VMD.dat',dump=.true.) ! Acc
  endif

  if (TryFile("0/HiLep.Rho0_MV.dat.bin")) then
     call FetchHist_10(hhh,file='HiLep.Rho0_MV.dat')
     call WriteHist(hhh,140,add=1e-20,mul=1./nLeptons,&
          & file='HiLep.Rho0_MV.dat',dump=.true.)

     call FetchHist_10(hhh,file='HiLep.Rho0_MX.dat')
     call WriteHist(hhh,140,add=1e-20,mul=1./nLeptons,&
          & file='HiLep.Rho0_MX.dat',dump=.true.)

     call FetchHist_10(hhh,file='HiLep.Rho0_DE.dat')
     call WriteHist(hhh,140,add=1e-20,mul=1./nLeptons,&
          & file='HiLep.Rho0_DE.dat',dump=.true.)
  end if

  if (TryFile("0/HiLep.Shadowing.dat.bin")) then
     call FetchHist2D_10(h2D,file='HiLep.Shadowing.dat')
     call WriteHist2D_Gnuplot(h2D,140, DoAve=.true.,MaxVal=-1.0,&
          & file='HiLep.Shadowing.dat',dump=.true.) ! Acc
  endif

  if (TryFile("0/HiLep.nuQ2.Flux.dat")) then
     call FetchHist2D_10(h2D,file='HiLep.nuQ2.Flux.dat')
     call WriteHist2D_Gnuplot(h2D,140,add=1e-20,mul=1./nLeptons,&
          & file='HiLep.nuQ2.Flux.dat',dump=.true.)
     call FetchHist2D_10(h2D,file='HiLep.nuQ2.FluxW.dat')
     call WriteHist2D_Gnuplot(h2D,140,add=1e-20,mul=1./nLeptons,&
          & file='HiLep.nuQ2.FluxW.dat',dump=.true.)
  end if
!!$  do i=0,4
!!$     if (TryFile("0/HiLep.nuQ2."//intTochar(i)//".dat")) then
!!$        call FetchHist2D_10(h2D,file='HiLep.nuQ2.'//intTochar(i)//'.dat')
!!$        call WriteHist2D_Gnuplot(h2D,140,add=1e-20,mul=1./nLeptons,&
!!$             & file='HiLep.nuQ2.'//intTochar(i)//'.dat',dump=.true.)
!!$     end if
!!$  end do

  if (TryFile("0/HiLep.Brooks.dat")) then
     call FetchHistMC_10(hMC,file='HiLep.Brooks.dat')
     call WriteHistMC(hMC,'HiLep.Brooks.dat',add=1e-20,mul=1./nLeptons)
  end if

  if (TryFile("0/HiLep.NuQ2planeXS.SYS.dat")) then
     do i=0,9
        call FetchHist2D_10(h2D,file='HiLep.NuQ2planeXS.'//trim(intToChar(i))//'.dat')
        call WriteHist2D_Gnuplot(h2D,140,add=1e-20,mul=1./nLeptons,&
             & file='HiLep.NuQ2planeXS.'//trim(intToChar(i))//'.dat',dump=.true.)
     end do
  end if

  if (TryFile("0/HiLep.CentralN_b.N0.dat")) then
     call FetchHist_10(HHH1,file='HiLep.CentralN_b.C_CN_b.dat')
     call FetchHist_10(HHH2,file='HiLep.CentralN_b.C_CN_bT.dat')
     call FetchHist_10(HHH3,file='HiLep.CentralN_b.C_CN_bZ.dat')
     
     do i=0,1
        call FetchHist2D_10(h2D,file='HiLep.CentralN_b.'//nucName(i)//'.dat')

        call WriteHist2D_Gnuplot(H2D, mul=1./nLeptons, add=1e-20,&
               & file='HiLep.CentralN_b.'//nucName(i)//'.dat',dump=.true.)
        call WriteHist2D_Gnuplot(H2D, MaxVal=0.0,H2=HHH1,&
               & file='HiLep.CentralN_b.norm.'//nucName(i)//'.dat')

        call FetchHist2D_10(h2D,file='HiLep.CentralN_bT.'//nucName(i)//'.dat')

        call WriteHist2D_Gnuplot(H2D, mul=1./nLeptons, add=1e-20,&
               & file='HiLep.CentralN_bT.'//nucName(i)//'.dat',dump=.true.)
        call WriteHist2D_Gnuplot(H2D, MaxVal=0.0,H2=HHH2,&
               & file='HiLep.CentralN_bT.norm.'//nucName(i)//'.dat')

        call FetchHist2D_10(h2D,file='HiLep.CentralN_bZ.'//nucName(i)//'.dat')

        call WriteHist2D_Gnuplot(H2D, mul=1./nLeptons, add=1e-20,&
               & file='HiLep.CentralN_bZ.'//nucName(i)//'.dat',dump=.true.)
        call WriteHist2D_Gnuplot(H2D, MaxVal=0.0,H2=HHH3,&
               & file='HiLep.CentralN_bZ.norm.'//nucName(i)//'.dat')

        call FetchHist2D_10(h2D,file='HiLep.CentralN_bZ_Ekin.'//nucName(i)//'.dat')

        call WriteHist2D_Gnuplot(H2D, mul=1./nLeptons, add=1e-20,&
               & file='HiLep.CentralN_bZ_Ekin.'//nucName(i)//'.dat',dump=.true.)
        call WriteHist2D_Gnuplot(H2D, MaxVal=0.0,H2=HHH3,&
               & file='HiLep.CentralN_bZ_Ekin.norm.'//nucName(i)//'.dat')

        call FetchHist2D_10(h2D,file='HiLep.CentralN.pTz.'//nucName(i)//'.dat')
        call WriteHist2D_Gnuplot(H2D, mul=1./nLeptons, add=1e-20,&
               & file='HiLep.CentralN.pTz.'//nucName(i)//'.dat',dump=.true.)
        
     end do
  end if

  if (TryFile("0/HiLep.EIC_zH.idH.Acc.dat")) then
     call FetchHistMP_10(hMP,"HiLep.EIC_zH.idH.Acc.dat")
     call WriteHistMP(hMP,add=1e-20,mul=1./nLeptons,iColumn=1,&
          &file='HiLep.EIC_zH.idH.Acc.dat',dump=.true.)

     call FetchHistMP_10(hMP,"HiLep.EIC_nu.idH.Acc.dat")
     call WriteHistMP(hMP,add=1e-20,H2=hLep_nu,iColumn=1,&
          &file='HiLep.EIC_nu.idH.Acc.dat',dump=.true.)

     call FetchHistMP_10(hMP,"HiLep.EIC_Q2.idH.Acc.dat")
     call WriteHistMP(hMP,add=1e-20,H2=hLep_Q2,iColumn=1,&
          &file='HiLep.EIC_Q2.idH.Acc.dat',dump=.true.)

     do i=1,5
        call FetchHistMP_10(hMP,'HiLep.EIC_zH.idH.Q2_'//Achar(i+48)//'.Acc.dat')
        call WriteHistMP(hMP,add=1e-20,mul=1./nLeptons,iColumn=1,&
             &file='HiLep.EIC_zH.idH.Q2_'//Achar(i+48)//'.Acc.dat',dump=.true.)

        call FetchHistMP_10(hMP,'HiLep.EIC_nu.idH.Q2_'//Achar(i+48)//'.Acc.dat')
        call WriteHistMP(hMP,add=1e-20,H2=hLep_nu,iColumn=1,&
             &file='HiLep.EIC_nu.idH.Q2_'//Achar(i+48)//'.Acc.dat',dump=.true.)
     end do

  end if

contains

  ! ---------------------

  subroutine FetchHist_10(H,file)
    use output

    implicit none
    type(histogram),intent(inout)       :: H
    character*(*),  intent(in)          :: file

    type(histogram) :: H1
    integer :: i
    logical :: flagOK

    if (.not.TryFile("0/"//file//".bin")) return

    call FetchHist(H,"0/"//file//".bin",flagOK=flagOK)
    if (.not.flagOK) return

    do i=1,9
       call FetchHist(H1,Achar(i+48)//"/"//file//".bin",flagOK=flagOK)
       if (.not.flagOK) cycle

       H%xExtreme(1) = min(H%xExtreme(1),H1%xExtreme(1))
       H%xExtreme(2) = max(H%xExtreme(2),H1%xExtreme(2))

       H%yVal = H%yVal + H1%yVal
    end do

  end subroutine FetchHist_10

  ! ---------------------

  subroutine FetchHistMP_10(H,file)
    use output

    implicit none
    type(histogramMP),intent(inout)       :: H
    character*(*),  intent(in)          :: file

    type(histogramMP) :: H1
    integer :: i
    logical :: flagOK

    if (.not.TryFile("0/"//file//".bin")) return

    call FetchHistMP(H,"0/"//file//".bin",flagOK=flagOK)
    if (.not.flagOK) return

    do i=1,9
       call FetchHistMP(H1,Achar(i+48)//"/"//file//".bin",flagOK=flagOK)
       if (.not.flagOK) cycle

       H%xExtreme(:,1) = min(H%xExtreme(:,1),H1%xExtreme(:,1))
       H%xExtreme(:,2) = max(H%xExtreme(:,2),H1%xExtreme(:,2))

       H%yVal = H%yVal + H1%yVal
    end do

  end subroutine FetchHistMP_10

  ! ---------------------

  subroutine FetchHistMC_10(H,file)
    use output

    implicit none
    type(histogramMC),intent(inout)       :: H
    character*(*),  intent(in)          :: file

    type(histogramMC) :: H1
    integer :: i
    logical :: flagOK

    if (.not.TryFile("0/"//file//".bin")) return

    call FetchHistMC(H,"0/"//file//".bin",flagOK=flagOK)
    if (.not.flagOK) return

    do i=1,9
       call FetchHistMC(H1,Achar(i+48)//"/"//file//".bin",flagOK=flagOK)
       if (.not.flagOK) cycle

       H%xExtreme(1) = min(H%xExtreme(1),H1%xExtreme(1))
       H%xExtreme(2) = max(H%xExtreme(2),H1%xExtreme(2))

       H%yVal = H%yVal + H1%yVal
    end do

  end subroutine FetchHistMC_10

  ! ---------------------

  subroutine FetchHist2D_10(H,file)
    use output

    implicit none
    type(histogram2D),intent(inout)       :: H
    character*(*),  intent(in)          :: file

    type(histogram2D) :: H1
    integer :: i
    logical :: flagOK

    if (.not.TryFile("0/"//file//".bin")) return

    call FetchHist2D(H,"0/"//file//".bin",flagOK=flagOK)
    if (.not.flagOK) return

    do i=1,9
       call FetchHist2D(H1,Achar(i+48)//"/"//file//".bin",flagOK=flagOK)
       if (.not.flagOK) cycle

       H%xExtreme(:,1) = min(H%xExtreme(:,1),H1%xExtreme(:,1))
       H%xExtreme(:,2) = max(H%xExtreme(:,2),H1%xExtreme(:,2))

       H%yVal = H%yVal + H1%yVal
    end do

  end subroutine FetchHist2D_10

  ! ---------------------

  logical function TryFile(file)
    implicit none
    character*(*),  intent(in)          :: file

    integer :: iVar

    open(9932,file=file,status='OLD',IOSTAT=iVar)
    close(9932)
    TryFile = (iVar.eq.0)

  end function TryFile


end program HiLeptonMix10Runs
