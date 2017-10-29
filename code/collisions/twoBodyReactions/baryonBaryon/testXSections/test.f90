program test
  use inputGeneral, only: readinputGeneral
  use particleProperties, only: InitParticleProperties
  use version, only: printVersion
  implicit none

  real, parameter :: deltaP = 0.01
  real, parameter :: deltaS = 0.01

  call printVersion
  call readinputGeneral
  call InitParticleProperties

  call test_NN_tot
  call test_NN_NNpi

!  call testBarBar_BarBar
!  call testBarBar_BarBar_chooseCharge
!   call testDimi

!   call test_NN_NNeta
!   call test_NN_NNpieta
!   call test_NN_NNrho
!   call test_NN_NNomega
!   call test_NN_NNpiomega
!   call test_np_deta

contains


!***********************************************************************
!***********************************************************************

subroutine test_NN_tot
  use particleDefinition
  use preEventDefinition
  use mediumDefinition
  use IDtable, only: nucleon
  use constants, only: mN
  use barBar_main, only: XsectionBarBar
  use master_2body, only: hiEnergyContrib

  integer :: i
  real :: srts,sigmaTot, sigmaElast
  logical :: pauliIncluded
  type(particle) :: teilchenIN(1:2)
  type(preEvent) :: teilchenOut(1:4)
  type(medium) :: mediumATcollision

  teilchenIn%ID = nucleon
  teilchenIn%mass = mN

  teilchenIn(1)%momentum = (/ mN, 0., 0., 0. /)
  teilchenIn(2)%momentum = (/ mN, 0., 0., 0. /)

  ! force initialization of module 'master_2body'
  print *,HiEnergyContrib(2.0,teilchenIn%Id,teilchenIn%antiParticle)

  write(*,*) '##########################################################'
  write(*,*) 'Testing pp -> X'
  write(*,*) '##########################################################'

  teilchenIn%charge = (/1,1/)

  do i=1,6000
    teilchenIn(2)%momentum(3) = i * deltaP  ! p_lab
    teilchenIn(2)%momentum(0) = Sqrt(teilchenIn(2)%mass**2+Dot_Product(teilchenIn(2)%momentum(1:3),teilchenIn(2)%momentum(1:3)))
    srts = sqrtS(teilchenIn)
    if (HiEnergyContrib(srts,teilchenIn%ID,teilchenIn%antiParticle)>=1.0) exit
    call XsectionBarBar(srts,teilchenIN,mediumATcollision,teilchenOUT,sigmaTot,sigmaElast,pauliIncluded,"pp")
  end do

  write(*,*) '##########################################################'
  write(*,*) 'Testing pn -> X'
  write(*,*) '##########################################################'

  teilchenIn%charge = (/1,0/)

  do i=1,6000
    teilchenIn(2)%momentum(3) = i * deltaP  ! p_lab
    teilchenIn(2)%momentum(0) = Sqrt(teilchenIn(2)%mass**2+Dot_Product(teilchenIn(2)%momentum(1:3),teilchenIn(2)%momentum(1:3)))
    srts = sqrtS(teilchenIn)
    if (HiEnergyContrib(srts,teilchenIn%ID,teilchenIn%antiParticle)>=1.0) exit
    call XsectionBarBar(srts,teilchenIN,mediumATcollision,teilchenOUT,sigmaTot,sigmaElast,pauliIncluded,"pn")
  end do

  write(*,*) '##########################################################'
  write(*,*) 'NN -> X: Subtraction'
  write(*,*) '##########################################################'

  call NN_tot_subtract_data("pp", 1)
  call NN_tot_subtract_data("np", 2)
  call NN_tot_subtract_data("pn", 2)

end subroutine test_NN_tot


subroutine NN_tot_subtract_data(filename, ch)
    use constants, only: mN
    use barBar_barBar, only: nukNuk_nukNuk

    character(len=*), intent(in) :: filename
    integer, intent(in) :: ch

    integer :: i, ios, charge(1:2)
    character(len=64) :: str
    real :: plab, plab_min, plab_max, sig_tot, sta_err(1:2), sy_err(1:2), srts, sig_el

    open(25,file="PDG_"//filename//"_total.dat")
    open(26,file=filename//"_inelast.dat")

    print *, filename

    do i=1,11
      ! discard file header
      read(25,'(A)') str
    end do

    if (ch==1) then
      charge = (/1,1/)
    else if (ch==2) then
      charge = (/1,0/)
    end if

    do
      read(25,*,iostat=ios) i, plab, plab_min, plab_max, sig_tot, sta_err(:), sy_err(:), str
      if (ios /= 0) exit
      ! print *, plab, sigma, err, str
      srts = sqrt(2*mN**2 + 2*mN*sqrt(mN**2+plab**2))
      sig_el = nukNuk_nukNuk (srts, (/mN,mN/), sum(charge))
      write(26,'(i4,8F10.5,2X,A)') i, plab, plab_min, plab_max, sig_tot-sig_el, sta_err(:), sy_err(:), str
    end do

    close(25)
    close(26)

end subroutine NN_tot_subtract_data


!***********************************************************************
!***********************************************************************


subroutine test_NN_NNpi
  use IdTable, only: nucleon
  use master_2body, only: HiEnergyContrib
  use barBar_barBarMes, only: Np_NRplus_NNPion, sstatepi

  real :: srtFree, sigmaDelta_p, sigma1_2_p, sigma3_2_p, sigmaDelta_n, sigma1_2_n, sigma3_2_n
  real, dimension(1:4) :: Sigma_NNPion, BG_Teis, BG_Buss, BG_Weil
  integer :: i, j

  real, parameter :: isofac12(1:4) = (/ 1./3., 1./3., 1./3., 2./3. /)
  real, parameter :: isofac32(1:4) = (/ 2./3., 1./3., 4./3., 10./3. /)

  write(*,*) '##########################################################'
  write(*,*) 'Testing NN -> NN pi'
  write(*,*) '##########################################################'

  Open(21,File='pp_pppi0.dat')
  Open(22,File='pn_pppi-.dat')
  Open(23,File='pn_pnpi0.dat')
  Open(24,File='pp_pnpi+.dat')

  do i=1,600
    srtFree = 2.0 + i * deltaS
    if (HiEnergyContrib(srtFree,(/nucleon,nucleon/),(/.false.,.false./))>=1.0) exit

    BG_Teis = sstatepi (srtFree, 1)                                           ! non-resonant background
    BG_Buss = sstatepi (srtFree, 2)
    BG_Weil = sstatepi (srtFree, 3)

    ! pp, full resonance contr.
    call Np_NRplus_NNPion (srtFree, .false., 1, sigmaDelta_p, sigma1_2_p, sigma3_2_p)
    Sigma_NNPion(1) = isofac12(1) * sigma1_2_p + isofac32(1) * (sigmaDelta_p+sigma3_2_p)
    Sigma_NNPion(4) = isofac12(4) * sigma1_2_p + isofac32(4) * (sigmaDelta_p+sigma3_2_p)

    ! pn, full resonance contr.
    call Np_NRplus_NNPion (srtFree, .false., 0, sigmaDelta_n, sigma1_2_n, sigma3_2_n)
    Sigma_NNPion(2) = isofac12(2) * sigma1_2_n + isofac32(2) * (sigmaDelta_n+sigma3_2_n)
    Sigma_NNPion(3) = isofac12(3) * sigma1_2_n + isofac32(3) * (sigmaDelta_n+sigma3_2_n)

    do j=1,4
      ! 1 = pp -> pp pi0
      ! 2 = pn -> pp pi-
      ! 3 = pn -> pn pi0
      ! 4 = pp -> pn pi+
      if (j==1 .or. j==4) then
        write(20+j,'(8F7.3)') srtFree, Sigma_NNPion(j), sigmaDelta_p*isofac32(j), sigma1_2_p*isofac12(j), sigma3_2_p*isofac32(j), &
                              BG_Teis(j), BG_Buss(j), BG_Weil(j)
      else
        write(20+j,'(8F7.3)') srtFree, Sigma_NNPion(j), sigmaDelta_n*isofac32(j), sigma1_2_n*isofac12(j), sigma3_2_n*isofac32(j), &
                              BG_Teis(j), BG_Buss(j), BG_Weil(j)
      end if
    end do
  end do

  close(21)
  close(22)
  close(23)
  close(24)

  !!!!! produce subtracted data for fitting background (data minus resonance contribution)

  write(*,*) '##########################################################'
  write(*,*) 'NN -> NN pi: Subtraction'
  write(*,*) '##########################################################'

  call singlePi_subtract_data("pp_pppi0", 1)
  call singlePi_subtract_data("pn_pppi-", 2)
  call singlePi_subtract_data("pn_pnpi0", 3)
  call singlePi_subtract_data("pp_pnpi+", 4)

end subroutine test_NN_NNpi


  subroutine singlePi_subtract_data(filename, ch)
    use constants, only: mN
    use barBar_barBarMes, only: Np_NRplus_NNPion

    character(len=*), intent(in) :: filename
    integer, intent(in) :: ch

    integer :: i, ios
    character(len=64) :: str
    real :: plab, sigma, err, srts, sigmaDelta, sigma1_2, sigma3_2, sigma_res

    real, parameter :: isofac12(1:4) = (/ 1./3., 1./3., 1./3., 2./3. /)
    real, parameter :: isofac32(1:4) = (/ 2./3., 1./3., 4./3., 10./3. /)

    open(25,file="LB_"//filename//".dat")
    open(26,file=filename//"_subtracted.dat")

    print *, filename

    do i=1,5
      ! discard file header
      read(25,'(A)') str
    end do

    do
      read(25,*,iostat=ios) plab, sigma, err, str
      if (ios /= 0) exit
      ! print *, plab, sigma, err, str
      srts = sqrt(2*mN**2 + 2*mN*sqrt(mN**2+plab**2))
      if (ch==1 .or. ch==4) then
        call Np_NRplus_NNPion (srts, .false., 1, sigmaDelta, Sigma1_2, Sigma3_2)       ! full resonance contr.
      else
        call Np_NRplus_NNPion (srts, .false., 0, sigmaDelta, Sigma1_2, Sigma3_2)
      end if
      sigma_res = isofac12(ch) * sigma1_2 + isofac32(ch) * (sigmaDelta+sigma3_2)
      write(26,'(4F10.5,2X,A)') plab, srts, sigma-sigma_res, err, str
    end do

    close(25)
    close(26)

  end subroutine


!***********************************************************************
!***********************************************************************


subroutine testbarBar_barBar
  ! to test the module barBar_barBar
  use particleDefinition
  use particleProperties, only: hadron
  use IdTable, only: nucleon, Delta, P11_1440
  use barBar_BarBar, only: sigmaBB
  use mediumDefinition, only: medium, vacuum
  use constants, only: mN

  type(particle), dimension(1:2) :: teilchenIN
  integer , dimension(1:2) ::idOut
  real :: srts
  integer :: i
  logical :: pauliIncluded

  write(*,*) '##########################################################'
  write(*,*) 'Testing barBar_barBar:sigmabb'
  write(*,*) '##########################################################'

  !*******************************************************************************

  ! (1) p p -> p p
  teilchenIn(1:2)%ID = nucleon
  teilchenIn(1:2)%mass = mN
  teilchenIN(1:2)%charge = 1
  idOut = (/nucleon,nucleon/)

  teilchenIn(1)%momentum(1:3)=0.
  teilchenIn(1)%momentum(0)=Sqrt(teilchenIn(1)%mass**2+Dot_Product(teilchenIn(1)%momentum(1:3),teilchenIn(1)%momentum(1:3)))
  teilchenIn(2)%momentum(1:3)=0.

  Open(100,File='barBar_Test_protProt.dat')
  write(100,*) '# srts, pLab, sigma proton proton-> proton proton'
  write(*,*) 'nuk nuk -> nuk  nuk'
  Do i=1,2000
     teilchenIn(2)%momentum(1)=i*0.0025
     teilchenIn(2)%momentum(0)=Sqrt(teilchenIn(2)%mass**2+Dot_Product(teilchenIn(2)%momentum(1:3),teilchenIn(2)%momentum(1:3)))
     srts=sqrtS(teilchenIn(1),teilchenIn(2))
     write(100,'(3F12.3)') srts, teilchenIn(2)%momentum(1),sigmaBB(teilchenIN,idOut,vacuum,srts,pauliIncluded)
     !write(*,'(3F10.3)') srts, teilchenIn(2)%momentum(1),sigmaBB(teilchenIN,idOut,vacuum,srts,pauliIncluded)
  end Do
  close(100)

  ! (2) n n -> n n
  teilchenIN(1:2)%charge = 0

  Open(100,File='barBar_Test_neutNeut.dat')
  write(100,*) '# srts, pLab, sigma neutron neutron-> neutron neutron'
  write(*,*) 'nuk nuk -> nuk  nuk'
  Do i=1,2000
     teilchenIn(2)%momentum(1)=i*0.0025
     teilchenIn(2)%momentum(0)=Sqrt(teilchenIn(2)%mass**2+Dot_Product(teilchenIn(2)%momentum(1:3),teilchenIn(2)%momentum(1:3)))
     srts=sqrtS(teilchenIn(1),teilchenIn(2))
     write(100,'(3F12.3)') srts, teilchenIn(2)%momentum(1),sigmaBB(teilchenIN,idOut,vacuum,srts,pauliIncluded)
     !write(*,'(3F10.3)') srts, teilchenIn(2)%momentum(1),sigmaBB(teilchenIN,idOut,vacuum,srts,pauliIncluded)
  end Do
  close(100)

  ! (3) n p -> n p
  teilchenIN(1:2)%charge = (/0,1/)

  Open(100,File='barBar_Test_neutProt.dat')
  write(100,*) '# srts, pLab, sigma neutron proton-> neutron proton'
  write(*,*) 'nuk nuk -> nuk  nuk'
  Do i=1,2000
     teilchenIn(2)%momentum(1)=i*0.0025
     teilchenIn(2)%momentum(0)=Sqrt(teilchenIn(2)%mass**2+Dot_Product(teilchenIn(2)%momentum(1:3),teilchenIn(2)%momentum(1:3)))
     srts=sqrtS(teilchenIn(1),teilchenIn(2))
     write(100,'(3F12.3)') srts, teilchenIn(2)%momentum(1),sigmaBB(teilchenIN,idOut,vacuum,srts,pauliIncluded)
     !write(*,'(3F10.3)') srts, teilchenIn(2)%momentum(1),sigmaBB(teilchenIN,idOut,vacuum,srts,pauliIncluded)
  end Do
  close(100)

  !*******************************************************************************

  ! (4) p p -> N D
  teilchenIN%charge = 1
  idOut = (/nucleon,delta/)

  Open(101,File='barBar_Test_nuknuk_nukdelta.dat')
  write(101,*) '# srts, pLab, sigma proton proton-> nucleon Delta'
  write(*,*) 'nuk nuk -> nuk  delta'
  Do i=1,2000
     teilchenIn(2)%momentum(1)=i*0.0025
     teilchenIn(2)%momentum(0)=Sqrt(teilchenIn(2)%mass**2+Dot_Product(teilchenIn(2)%momentum(1:3),teilchenIn(2)%momentum(1:3)))
     srts=sqrtS(teilchenIn(1),teilchenIn(2))
     write(101,'(3F10.3)') srts, teilchenIn(2)%momentum(1),sigmaBB(teilchenIN,idOut,vacuum,srts,pauliIncluded)
     !write(*,'(3F10.3)') srts, teilchenIn(2)%momentum(1),sigmaBB(teilchenIN,idOut,vacuum,srts,pauliIncluded)
  End do
  close(101)

  ! (5) p D+ -> p p
  teilchenIn(2)%ID=delta
  teilchenIn(2)%mass=hadron(delta)%mass
  teilchenIN%charge = 1
  idOut = (/nucleon,nucleon/)

  Open(101,File='barBar_Test_nukdelta_nuknuk.dat')
  write(101,*) '# srts, pLab, sigma proton delta+-> nucleon nucleon'
  write(*,*) 'nuk nuk -> nuk  delta'
  Do i=1,2000
     teilchenIn(2)%momentum(1)=i*0.0025
     teilchenIn(2)%momentum(0)=Sqrt(teilchenIn(2)%mass**2+Dot_Product(teilchenIn(2)%momentum(1:3),teilchenIn(2)%momentum(1:3)))
     srts=sqrtS(teilchenIn(1),teilchenIn(2))
     write(101,'(3F10.3)') srts, teilchenIn(2)%momentum(1),sigmaBB(teilchenIN,idOut,vacuum,srts,pauliIncluded)
     !write(*,'(3F10.3)') srts, teilchenIn(2)%momentum(1),sigmaBB(teilchenIN,idOut,vacuum,srts,pauliIncluded)
  End do
  close(101)

  ! (6) p D+ -> N D
  teilchenIn(2)%ID=delta
  teilchenIn(2)%mass=hadron(delta)%mass
  teilchenIN%charge = 1
  idOut = (/nucleon,delta/)

  Open(101,File='barBar_Test_nukdelta_nukdelta.dat')
  write(101,*) '# srts, pLab, sigma proton delta+-> nucleon delta'
  write(*,*) 'nuk nuk -> nuk  delta'
  Do i=1,2000
     teilchenIn(2)%momentum(1)=i*0.0025
     teilchenIn(2)%momentum(0)=Sqrt(teilchenIn(2)%mass**2+Dot_Product(teilchenIn(2)%momentum(1:3),teilchenIn(2)%momentum(1:3)))
     srts=sqrtS(teilchenIn(1),teilchenIn(2))
     write(101,'(3F10.3)') srts, teilchenIn(2)%momentum(1),sigmaBB(teilchenIN,idOut,vacuum,srts,pauliIncluded)
     !write(*,'(3F10.3)') srts, teilchenIn(2)%momentum(1),sigmaBB(teilchenIN,idOut,vacuum,srts,pauliIncluded)
  End do
  close(101)

  !*******************************************************************************

  ! (7) n p -> D D
  teilchenIn(2)%ID=nucleon
  teilchenIn(2)%mass=mN
  teilchenIN%charge = (/0,1/)
  idOut = (/delta,delta/)

  Open(102,File='barBar_Test_nuknuk_deltadelta.dat')
  write(102,*) '# srts, pLab, sigma proton neutron-> Delta Delta'
  write(*,*) 'nuk nuk -> delta  delta'
  Do i=1,200
     teilchenIn(2)%momentum(1)=i*0.05
     teilchenIn(2)%momentum(0)=Sqrt(teilchenIn(2)%mass**2+Dot_Product(teilchenIn(2)%momentum(1:3),teilchenIn(2)%momentum(1:3)))
     srts=sqrtS(teilchenIn(1),teilchenIn(2))
     write(102,'(3F10.3)') srts, teilchenIn(2)%momentum(1),sigmaBB(teilchenIN,idOut,vacuum,srts,pauliIncluded)
     !write(*,'(3F10.3)') srts, teilchenIn(2)%momentum(1),sigmaBB(teilchenIN,idOut,vacuum,srts,pauliIncluded)
  End do
  close(102)

  !*******************************************************************************

  ! (8) D0 D+ -> n p
  teilchenIn(1:2)%ID = delta
  teilchenIN(1:2)%charge = (/1,0/)
  teilchenIn(1:2)%mass = hadron(delta)%mass
  idOut = (/nucleon,nucleon/)

  teilchenIn(1)%momentum(0)=Sqrt(teilchenIn(1)%mass**2+Dot_Product(teilchenIn(1)%momentum(1:3),teilchenIn(1)%momentum(1:3)))

  Open(103,File='barBar_Test_deltadelta_nuknuk.dat')
  write(103,*) '# srts, pLab, sigma proton neutron <- Delta Delta'
  write(*,*) 'delta  delta -> nuk nuk'
  Do i=1,200
     teilchenIn(2)%momentum(1)=i*0.02
     teilchenIn(2)%momentum(0)=Sqrt(teilchenIn(2)%mass**2+Dot_Product(teilchenIn(2)%momentum(1:3),teilchenIn(2)%momentum(1:3)))
     srts=sqrtS(teilchenIn(1),teilchenIn(2))
     write(103,'(3F10.3)') srts, teilchenIn(2)%momentum(1),sigmaBB(teilchenIN,idOut,vacuum,srts,pauliIncluded)
     !write(*,'(3F10.3)') srts, teilchenIn(2)%momentum(1),sigmaBB(teilchenIN,idOut,vacuum,srts,pauliIncluded)
  End do
  close(103)

  !************************************************************************************************

  ! (9) n p -> N P11_1440
  teilchenIn(1:2)%ID = nucleon
  teilchenIN(1:2)%charge = (/1,0/)
  teilchenIn(1:2)%mass = mN
  idOut = (/nucleon,P11_1440/)

  teilchenIn(1)%momentum(0)=Sqrt(teilchenIn(1)%mass**2+Dot_Product(teilchenIn(1)%momentum(1:3),teilchenIn(1)%momentum(1:3)))

  Open(104,File='barBar_Test_nuknuk_nukP11.dat')
  write(104,*) '# srts, pLab, sigma proton neutron -> nucleon P11_1440'
  write(*,*) 'nuk nuk -> nuk P11_1440'
  Do i=1,200
     teilchenIn(2)%momentum(1)=i*0.02
     teilchenIn(2)%momentum(0)=Sqrt(teilchenIn(2)%mass**2+Dot_Product(teilchenIn(2)%momentum(1:3),teilchenIn(2)%momentum(1:3)))
     srts=sqrtS(teilchenIn(1),teilchenIn(2))
     write(104,'(3F10.3)') srts, teilchenIn(2)%momentum(1),sigmaBB(teilchenIN,idOut,vacuum,srts,pauliIncluded)
     !write(*,'(3F10.3)') srts, teilchenIn(2)%momentum(1),sigmaBB(teilchenIN,idOut,vacuum,srts,pauliIncluded)
  End do
  close(104)

  !************************************************************************************************

  ! (10) p P11_1440 -> n p
  teilchenIn%ID=(/nucleon,P11_1440/)
  teilchenIN%charge=(/1,0/)
  teilchenIn(1)%mass=mN
  teilchenIn(2)%mass=hadron(P11_1440)%mass

  Open(105,File='barBar_Test_nukP11_nuknuk.dat')
  write(105,*) '# srts, pLab, sigma proton neutron <- nucleon P11_1440'
  write(*,*) 'nuk nuk -> nuk P11_1440'
  Do i=1,200
     teilchenIn(2)%momentum(1)=i*0.02
     teilchenIn(2)%momentum(0)=Sqrt(teilchenIn(2)%mass**2+Dot_Product(teilchenIn(2)%momentum(1:3),teilchenIn(2)%momentum(1:3)))
     srts=sqrtS(teilchenIn(1),teilchenIn(2))

     idOut=(/nucleon,nucleon/)
     write(105,'(3F10.3)') srts, teilchenIn(2)%momentum(1),sigmaBB(teilchenIN,idOut,vacuum,srts,pauliIncluded)
     !write(*,'(3F10.3)') srts, teilchenIn(2)%momentum(1),sigmaBB(teilchenIN,idOut,vacuum,srts,pauliIncluded)

  End do
  close(104)

end subroutine testbarBar_barBar


!***********************************************************************
!***********************************************************************


subroutine testbarBar_barBar_chooseCharge
  ! to test the module barBar_barBar
  use particleDefinition
  use particleProperties, only: hadron
  use IdTable, only: nucleon, Delta
  use barBar_BarBar, only: chooseCharge
  use constants, only: mN

  type(particle), dimension(1:2) :: teilchenIN
  integer , dimension(1:2) ::idOut, chargeOut
  integer :: i, count

  write(*,*) '##########################################################'
  write(*,*) 'Testing barBar_barBar:chooseCharge'
  write(*,*) '##########################################################'

  !!! Proton Proton -> Delta nucleon
  teilchenIn%ID=nucleon
  teilchenIN%charge=1
  teilchenIn(1)%mass=mN
  teilchenIn(1)%momentum(1:3)=0.
  teilchenIn(1)%momentum(0)=Sqrt(teilchenIn(1)%mass**2+Dot_Product(teilchenIn(1)%momentum(1:3),teilchenIn(1)%momentum(1:3)))

  teilchenIn(2)%mass=mN
  teilchenIn(2)%momentum(1:3)=0.
  teilchenIn(2)%momentum(0)=Sqrt(teilchenIn(2)%mass**2+Dot_Product(teilchenIn(2)%momentum(1:3),teilchenIn(2)%momentum(1:3)))

  idOut=(/Delta, nucleon/)

  count=0
  do i=1,1000
     chargeOut = chooseCharge(teilchenIN,idOut)
!     Write (*,*) chargeOut
     if (chargeOut(1)==1) count=count+1
  End do
  write(*,*) 'Probability for 1,1:', count/1000.


  !!! Proton Delta+ -> Delta nucleon
  teilchenIn%ID=(/nucleon,delta/)
  teilchenIN%charge=1
  teilchenIn(1)%mass=mN
  teilchenIn(1)%momentum(1:3)=0.
  teilchenIn(1)%momentum(0)=Sqrt(teilchenIn(1)%mass**2+Dot_Product(teilchenIn(1)%momentum(1:3),teilchenIn(1)%momentum(1:3)))

  teilchenIn(2)%mass=hadron(delta)%mass
  teilchenIn(2)%momentum(1:3)=0.
  teilchenIn(2)%momentum(0)=Sqrt(teilchenIn(2)%mass**2+Dot_Product(teilchenIn(2)%momentum(1:3),teilchenIn(2)%momentum(1:3)))

  idOut=(/Delta, nucleon/)

  count=0
  do i=1,1000
     chargeOut = chooseCharge(teilchenIN,idOut)
!     Write (*,*) chargeOut
     if (chargeOut(1)==1) count=count+1
  End do
  write(*,*) 'Probability for 1,1:', count/1000.

end subroutine testbarBar_barBar_chooseCharge


!************************************************************************
!************************************************************************


subroutine testdimi
  ! to test dimi.f90
  use dimi, only: dimidm
  use particleDefinition
  use twobodyTools, only : pcm
  use constants, only: mN

  real, parameter :: mdel=1.232
  real, parameter :: deltaM=0.001
  integer :: i
  real :: srtS,mass, p_cm, dsdm, dsdm_integrated, sigma

  write(*,*) '##########################################################'
  write(*,*) 'Testing dimi'
  write(*,*) '##########################################################'

  Open(100,File='dimiTest_srts.dat')
  write(100,*) '# srts, ds/dm, Integral(ds/dm, dm), sigma'
  Do i=1,300
     srts = 1. + deltaS*real(i)
     call dimidm(dsdm,dsdm_Integrated,sigma,srts,mdel)
     write(100,'(4F10.3)') srts,dsdm, dsdm_Integrated, sigma
  End do
  close(100)

  Open(100,File='dimiTest_mass1.dat')
  srts=2.05
  write(100,*) '# srts=',srts
  write(100,*) 'mass, ds/dm, Integral(ds/dm, dm), sigma'
  p_cm=pcm(srts,mN,mN)
  Do i=1,500
     mass=0.9+deltaM*real(i)
     call dimidm(dsdm,dsdm_Integrated,sigma,srts,mass)
     write(100,'(5F10.3)') mass,dsdm, dsdm_Integrated, sigma,dsdm/p_cm
  End do
  close(100)

  srts=2.51
  Open(100,File='dimiTest_mass2.dat')
  write(100,*) '# srts=',srts
  write(100,*) 'mass, ds/dm, Integral(ds/dm, dm), sigma'
  p_cm=pcm(srts,mN,mN)
  Do i=1,500
     mass=0.9+deltaM*real(i)
     call dimidm(dsdm,dsdm_Integrated,sigma,srts,mass)
     write(100,'(5F10.3)') mass,dsdm, dsdm_Integrated, sigma,dsdm/p_cm
  End do
  close(100)

end subroutine testdimi


!************************************************************************
!************************************************************************


  subroutine test_NN_NNeta
    ! test for single eta production via resonances
    use IdTable
    use particleProperties, only: hadron
    use constants, only: mN
    use barBar_BarBar, only : NN_NRes
    integer :: i, resID, isoIn(1:2), chargeIn(1:2), idOut(1:2)
    real :: srts, pcm, sig, sigma1_2, sigma3_2, sigma1535

    write(*,*) '##########################################################'
    write(*,*) 'Testing p p -> p p eta'
    write(*,*) '##########################################################'

    ! Define two protons in the incoming channel
    isoIn = (/1,1/) ! Isospin of incoming nucleons times 2
    chargeIN = (/1,1/) ! Charges of incoming nucleons

    Open(201,File='pp_eta.dat')
    do i=1,400

      srts = 2.*mN + i*deltaS
      pcm = sqrt(Max(0.,srts**2/4.-mN**2))

      ! Initialize output
      sigma1_2 = 0.
      sigma3_2 = 0.
      sigma1535 = 0.

      Do resID = delta+1, F37_1950
         If(.not. hadron(resID)%usedForXsections ) cycle         ! Exclude resonances from the Xsections
         idOut = (/nucleon,resID/)
         sig = NN_NRes(srts,pCM,idOut,isoIn,chargeIN,3)
         If (hadron(resID)%IsoSpinTimes2 == 1) then ! Isospin 1/2 resonances
            sigma1_2 = sigma1_2  + sig
         else if (hadron(resID)%IsoSpinTimes2 == 3) then ! Isospin 3/2 resonances
            ! factor 1/4 because the NN_NRES gives the result summed over all final states. And Clebsch(p p->p R+)=1/4 for resonance R with I=3/2.
            sigma3_2 = sigma3_2  + sig/4.
         end if
         if (resID == S11_1535) sigma1535 = sig
      end do

      write (201,'(4F10.6)') srts, sigma1_2, sigma3_2, sigma1535

    end do
    close(201)

  end subroutine


!************************************************************************
!************************************************************************


  subroutine test_NN_NNpieta
    ! test for pi eta production via resonances
    use IdTable
    use particleProperties, only: hadron
    use constants, only: mN
    use barBar_BarBar, only : NN_ResRes

    integer :: i, isoIn(1:2), chargeIn(1:2), idOut(1:2)
    real :: srts, pcm, sigma, BR

    write(*,*) '##########################################################'
    write(*,*) 'Testing p p -> p p pi eta'
    write(*,*) '##########################################################'

    ! Define two protons in the incoming channel
    isoIn = (/1,1/) ! Isospin of incoming nucleons times 2
    chargeIN = (/1,1/) ! Charges of incoming nucleons

    Open(201,File='pp_pi_eta.dat')
    do i=1,400

      srts = 2.*mN + i*deltaS
      pcm = sqrt(Max(0.,srts**2/4.-mN**2))

      idOut = (/Delta,S11_1535/)

      BR = hadron(S11_1535)%decays(2)                       ! branching ratio into eta N

      sigma = NN_ResRes (srts,pCM,idOut,isoIn,chargeIN) * BR

      write (201,'(4F10.6)') srts, sigma

    end do
    close(201)

  end subroutine


!************************************************************************
!************************************************************************


  subroutine test_NN_NNrho
    ! test for single rho0 production via resonances
    use IdTable
    use particleProperties, only: hadron
    use constants, only: mN
    use barBar_BarBar, only : NN_NRes
    integer :: i, resID, isoIn(1:2), chargeIn(1:2), idOut(1:2)
    real :: srts, pcm, sig, sigma1_2, sigma3_2, sigR(2:31)

    write(*,*) '##########################################################'
    write(*,*) 'Testing p p -> p p rho^0'
    write(*,*) '##########################################################'

    ! Define two protons in the incoming channel
    isoIn = (/1,1/) ! Isospin of incoming nucleons times 2
    chargeIN = (/1,1/) ! Charges of incoming nucleons

    Open(202,File='pp_rho.dat')
    do i=1,400

      srts = 2.*mN + i*deltaS
      pcm = sqrt(Max(0.,srts**2/4.-mN**2))

      ! Initialize output
      sigma1_2 = 0.
      sigma3_2 = 0.
      sigR = 0.

      Do resID = delta+1, F37_1950
         If(.not. hadron(resID)%usedForXsections ) cycle         ! Exclude resonances from the Xsections
         idOut = (/nucleon,resID/)
         sig = NN_NRes(srts,pCM,idOut,isoIn,chargeIN,4)
         If (hadron(resID)%IsoSpinTimes2 == 1) then ! Isospin 1/2 resonances, isospin factor 1/3 for the N* -> N rho^0 decay
            sigR(resID) = sig * 1./3.
            sigma1_2 = sigma1_2  + sigR(resID)
         else if (hadron(resID)%IsoSpinTimes2 == 3) then ! Isospin 3/2 resonances, isospin factor 2/3 for the D* -> N rho^0 decay
            ! factor 1/4 because the NN_NRES gives the result summed over all final states. And Clebsch(p p->p R+)=1/4 for resonance R with I=3/2.
            sigR(resID) = sig * 1./4. * 2./3.
            sigma3_2 = sigma3_2  + sigR(resID)
         end if
      end do

      write (202,'(33F10.6)') srts, sigma1_2, sigma3_2, sigR

    end do
    close(202)

  end subroutine


!************************************************************************
!************************************************************************


  subroutine test_NN_NNomega
    ! test for single omega production via resonances
    use IdTable
    use particleProperties, only: hadron
    use constants, only: mN
    use barBar_BarBar, only : NN_NRes
    integer :: i, resID, isoIn(1:2), chargeIn(1:2), idOut(1:2)
    real :: srts, pcm, sig, sigma1_2, sigma3_2, sigR(2:31)

    write(*,*) '##########################################################'
    write(*,*) 'Testing p p -> p p omega'
    write(*,*) '##########################################################'

    ! Define two protons in the incoming channel
    isoIn = (/1,1/) ! Isospin of incoming nucleons times 2
    chargeIN = (/1,1/) ! Charges of incoming nucleons

    Open(203,File='pp_omega.dat')
    do i=1,400

      srts = 2.*mN + i*deltaS
      pcm = sqrt(Max(0.,srts**2/4.-mN**2))

      ! Initialize output
      sigma1_2 = 0.
      sigma3_2 = 0.
      sigR = 0.

      Do resID = delta+1, F37_1950
         If(.not. hadron(resID)%usedForXsections ) cycle         ! Exclude resonances from the Xsections
         idOut = (/nucleon,resID/)
         sig = NN_NRes(srts,pCM,idOut,isoIn,chargeIN,5)
         If (hadron(resID)%IsoSpinTimes2 == 1) then ! Isospin 1/2 resonances
            sigR(resID) = sig
            sigma1_2 = sigma1_2  + sigR(resID)
         else if (hadron(resID)%IsoSpinTimes2 == 3) then ! Isospin 3/2 resonances
            ! factor 1/4 because the NN_NRES gives the result summed over all final states. And Clebsch(p p->p R+)=1/4 for resonance R with I=3/2.
            sigR(resID) = sig/4.
            sigma3_2 = sigma3_2  + sigR(resID)
         end if
      end do

      write (203,'(33F10.6)') srts, sigma1_2, sigma3_2, sigR

    end do
    close(203)

  end subroutine


!************************************************************************
!************************************************************************


  subroutine test_NN_NNpiomega
    ! test for pi omega production via resonances
    use IdTable
    use particleProperties, only: hadron
    use constants, only: mN
    use barBar_BarBar, only : NN_ResRes

    integer :: i, isoIn(1:2), chargeIn(1:2), idOut(1:2)
    real :: srts, pcm, sigma_1900, sigma_2190, BR

    write(*,*) '##########################################################'
    write(*,*) 'Testing p p -> p p pi omega'
    write(*,*) '##########################################################'

    ! Define two protons in the incoming channel
    isoIn = (/1,1/) ! Isospin of incoming nucleons times 2
    chargeIN = (/1,1/) ! Charges of incoming nucleons

    Open(201,File='pp_pi_omega.dat')
    do i=1,400

      srts = 2.*mN + i*deltaS
      pcm = sqrt(Max(0.,srts**2/4.-mN**2))

      idOut = (/Delta,P13_1900/)
      BR = hadron(P13_1900)%decays(2)                       ! branching ratio into omega N
      sigma_1900 = NN_ResRes (srts,pCM,idOut,isoIn,chargeIN) * BR

      idOut = (/Delta,G17_2190/)
      BR = hadron(G17_2190)%decays(2)                       ! branching ratio into omega N
      sigma_2190 = NN_ResRes (srts,pCM,idOut,isoIn,chargeIN) * BR

      write (201,'(4F10.6)') srts, sigma_1900+sigma_2190, sigma_1900, sigma_2190

    end do
    close(201)

  end subroutine


!************************************************************************
!************************************************************************

  subroutine test_np_deta
    use barBar_Main, only : eta_deuteron
    use constants, only: mN
    integer :: i
    real :: srts

    Open(201,File='eta_deuteron.dat')

    do i=1,300
      srts = 2.*mN + 0.400 + i*0.001
      write (201,'(2F10.6)') srts, eta_deuteron(srts)
    end do

    close(201)

  end subroutine

end program test
