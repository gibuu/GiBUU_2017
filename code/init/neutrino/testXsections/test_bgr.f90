program test_bgr

  use leptonicID
  use constants
  use output
  use version
  use inputGeneral
  use insertion, only: GarbageCollection
  use densityModule
  use particleProperties, only : baryon, meson
  use idTable, only : nucleon, pion, Delta
  use lepton_kinematics_free
  use lepton_kinematics_free_FULL
  use lepton_xsec_free

  implicit none

  integer     :: process_ID=2 ! 1=EM 2=CC
  integer     :: out_lepton=2 ! 1=electron  2=muon 3=tau
  logical     :: incl_E1_ct_ctWpi=.true., incl_E1_ct_Epi=.true.
  logical     :: incl_ct_E1=.false., incl_E1=.false., incl_Epi=.false. ! which xsec to calculate
  logical     :: debugFlag=.false., firstTime=.true. ! in the namelist  nl_test_bgr
  logical     :: cNpCTpppF=.true. ! calculate cNp, CT, pp, pF diagrams if .true., calculate Dp if .false.
  logical     :: incl_minmaxctqPi_Q2_W2_Epi=.true., incl_minmaxt_Q2_W2=.true., incl_minmaxctnupi_Enu_sl_Epi=.true., &
       & incl_minmaxTpi_Enu_sl=.true., incl_minmaxEpi_Enu=.true.

  real        :: mN, mN1, mpi, ml
  real        :: Enu=1.2, E1=0.45, costheta=0.85 ! in the namelist  nl_test_bgr
  real        :: costhetamin,costhetamax
  integer     :: i, k
  real        :: xsec_CT, xsec_pp, xsec_pF
  character(45) :: directory='results-bgr', name='dummy', name_lepton, name_process, name_xsec



  call PrintVersion
  call readInputGeneral
  call init_DataBase  ! Masses, Name, Index, width of the decay channels


  ! program here

  mN=baryon(Nucleon)%mass
  mN1=mN
  mpi=meson(pion)%mass

  write(*,*) 'Nieves cross sections for 1-pion final state'

  call readinput_test_bgr()

  select case(out_lepton)
  case (1)
     ml=melec
     name_lepton='elec_'
  case (2)
     ml=mmuon
     name_lepton='muon_'
  case (3)
     ml=mtau
     name_lepton='tau_'
  case default
     stop 'Wrong outgoing lepton. It should be 1(electon) or 2(muon)  or 3(tau)'
  end select


  select case(process_ID)
  case (EM)
     name_process='EM_'
  case (CC)
     name_process='CC_'
  case default
     stop 'Wrong process_ID. It should be 1(EM) or 2(CC)'
  end select





  ! ****************************
  ! show kinematics
  ! ****************************

  if (incl_minmaxctqPi_Q2_W2_Epi) call show_minmaxctqPi_Q2_W2_Epi()

  if (incl_minmaxt_Q2_W2) call show_minmaxt_Q2_W2()

  if (incl_minmaxctnupi_Enu_sl_Epi) call show_minmaxctnupi_Enu_sl_Epi()

  if (incl_minmaxTpi_Enu_sl) call show_minmaxTpi_Enu_sl()

  if (incl_minmaxEpi_Enu) call show_minmaxEpi_Enu()

  ! ****************************
  ! triple differential versus costheta_Wpi
  ! ****************************

  if (incl_E1_ct_ctWpi) then

     name_xsec='E1ctctWpi'
     name=trim(directory) // '/' //'Nieves_'//trim(name_process)//trim(name_lepton)//trim(name_xsec)//'.dat'

     Open(10,File=name)
     write(10,*)'# Nieves 1-pion production  d', name_xsec, '    ', name_process, '    ', name_lepton

     ! min and max costheta for a given Enu and nu
     call minmaxcostheta_Enu_nu(Enu,(Enu-E1),(mN+mpi),mN,ml,costhetamin,costhetamax)
     write(10,'(2(A,g10.4))') '# costhetamin=', costhetamin, '    costhetamax=', costhetamax

     do i=1,8
        costheta=0.1+0.1*i
        if ( costheta<costhetamin  .or.  costheta>costhetamax ) then
           cycle
        else
    !     call Nieves_E1_ct_ctWpi(CC,melec,Enu,E1,costheta,debugflag)
           call show_HNV_free_elepton_ct_ctPi(process_ID,ml,Enu,costheta,E1,debugflag)

        end if
        write (10,*) ''
        write (10,*) ''
     end do
     close(10)

  end if




  ! ****************************
  ! triple differential versus Epi
  ! ****************************

  if (incl_E1_ct_Epi) then

     name_xsec='E1ctEpi'
     name=trim(directory) // '/' //'Nieves_'//trim(name_process)//trim(name_lepton)//trim(name_xsec)//'.dat'
     Open(10,File=name)
     write(10,*)'# Nieves 1-pion production  d', name_xsec, '    ', name_process, '    ', name_lepton

     ! min and max costheta for a given Enu and nu
     call minmaxcostheta_Enu_nu(Enu,(Enu-E1),(mN+mpi),mN,ml,costhetamin,costhetamax)
     write(10,'(2(A,g10.4))') '# costhetamin=', costhetamin, '    costhetamax=', costhetamax

     do i=1,8
        costheta=0.1+0.1*i
        if (costheta<costhetamin .or. costheta>costhetamax) then
           cycle
        else
           call Nieves_E1_ct_Epi(process_ID,ml,Enu,E1,costheta,debugflag)
        end if
        write (10,*) ''
        write (10,*) ''
     end do
     close(10)

  end if


  ! ****************************
  ! double differential     dsi/dcostheta dE1        versus E1
  ! ****************************
  if (incl_ct_E1) then

     name_xsec='ct-E1'
     write(*,'(8A12)') 'directory=',directory, '    name_process=',name_process, &
          & '     name_lepton=', name_lepton, '    name_xsec=',name_xsec
     name=trim(directory)//'/'//'Nieves_'//trim(name_process)//trim(name_lepton)//trim(name_xsec)//'.dat'
     write(*,*) 'name=',name
     Open(10,File=name)

     write(10,*)'# Nieves 1-pion production  d', name_xsec, '    ', name_process, '    ', name_lepton

     ! min and max costheta for a given Enu (-1 and 1)
     call minmaxcostheta_Enu(Enu,(mN+mpi),mN,ml,costhetamin,costhetamax)
     write(10,'(2(A,g10.4))') '# costhetamin=', costhetamin, '       costhetamax=', costhetamax

     if (costheta<costhetamin .or. costheta>costhetamax) then
     else
        write (10,*) '# ------------------------------------------------------'
        if (cNpCTpppF .eq. .false.) call Nieves_ct_E1(process_ID,ml,Enu,costheta,debugflag) ! for direct Delta, integrated over Epi
        !call Nieves_ct_E1_(CC,mmuon,Enu,costheta,debugflag) ! for direct Delta, integrated over costheta_Wpi
        if (cNpCTpppF .eq. .true.) call show_HNV_free_elepton_ct(process_ID,ml,Enu,costheta,debugflag)
        write (10,*) ''
        write (10,*) ''
        write (10,*) ''

        !call Nieves_ct_E1_(CC,mmuon,Enu,costheta,debugflag) ! for direct Delta
     end if
  end if


  ! ****************************
  ! one differential dsi/dE1   versus E1
  ! ****************************
  if (incl_E1) then
     name_xsec='E1'
     name=trim(directory) // '/' //'Nieves_'//trim(name_process)//trim(name_lepton)//trim(name_xsec)//'.dat'
     Open(10,File=name)

     write(10,*)'# Nieves 1-pion production  d', name_xsec, '    ', name_process, '    ', name_lepton

     write (10,*) '# ------------------------------------------------------'
     call Nieves_E1(process_ID,ml,Enu,debugflag)
     write (10,*) ''
     write (10,*) ''
     write (10,*) ''
     write (10,*) ''
  end if
  close(10)




  ! ****************************
  ! one differential dsi/dEpi   versus Epi
  ! ****************************
  if (incl_Epi) then
     name_xsec='Epi'
     name=trim(directory) // '/' //'Nieves_'//trim(name_process)//trim(name_lepton)//trim(name_xsec)//'.dat'
     Open(10,File=name)

     write(10,*)'# Nieves 1-pion production  d', name_xsec, '    ', name_process, '    ', name_lepton

     do i=0,3
        Enu = 0.5 + 0.5*i
        write (10,*) '# ------------------------------------------------------'
        call Nieves_Epi(CC,mmuon,Enu,debugflag)
        write (10,*) ''
        write (10,*) ''
        write (10,*) ''
        write (10,*) ''
     end do

     close(10)
  end if


contains


  subroutine readinput_test_bgr()
    use output
    implicit none
    integer :: IOS

    !***************************************************************************
    !****n* test_bgr/nl_test_bgr
    ! NAME
    ! NAMELIST /nl_test_bgr/
    ! PURPOSE
    ! This Namelist includes:
    ! * debugflag
    ! * Enu
    ! * E1
    ! * costheta
    !***************************************************************************
    NAMELIST /nl_test_bgr/ process_ID, out_lepton, debugflag, cNpCTpppF, &
         &incl_E1_ct_ctWpi, incl_E1_ct_Epi, incl_ct_E1, incl_E1, incl_Epi, &
         & Enu, E1, costheta, &
         & incl_minmaxctqPi_Q2_W2_Epi, incl_minmaxt_Q2_W2, incl_minmaxctnupi_Enu_sl_Epi, incl_minmaxTpi_Enu_sl, &
         & incl_minmaxEpi_Enu

    call Write_ReadingInput('nl_test_bgr',0)
    rewind(5)
    read(5,nml=nl_test_bgr,IOSTAT=IOS)
    call Write_ReadingInput('nl_test_bgr',0,IOS) ! if IOS=0 writes nothing
    call Write_ReadingInput('nl_test_bgr',1)

    write(*,'(3(A,g10.3))') 'Input parameters for test_bgr:    Enu=', Enu, '    E1=', E1, '    costheta=', costheta

  end subroutine readinput_test_bgr




  ! provides output for 1-fold differential xsec dsi/dEpi  versus Epi
  subroutine Nieves_Epi(process_ID,ml_out,Enu,debugflag)

    use particleProperties, only : baryon, meson
    use idTable, only : nucleon, pion, Delta
    use random, only : rn
    use lepton_kinematics_free
    use lepton_kinematics_free_FULL
    use gauss_integration


    implicit none
    integer, intent(in) :: process_ID
    real, intent(in) :: ml_out, Enu
    logical, optional, intent(in) :: debugflag


    logical :: debug
    real :: mN, mDelta, mpi
    real :: Epi, Epimin, Epimax
    real, dimension(0:3) :: k_in, p_in
    real, dimension(1:3) :: position
    integer :: n1, n, j
    real, allocatable, dimension(:) :: x, xsec_neut, xsec_prot
    real :: sig_neut, sig_prot


    if (present(debugflag)) then
       debug=debugflag
    else
       debug=.true.
    end if

    write(10,*) '# Differential x-sec: dsi/Epi:   Enu, Epi, neut, prot'
    write(10,'(1(A,g10.4))') '# Enu=', Enu

    mN=baryon(nucleon)%mass
    mDelta=baryon(Delta)%mass
    mpi=meson(pion)%mass


    call  minmaxEpi_Enu( Enu , mN , mpi , mN , ml_out , Epimin , Epimax )
    if (debug) write(10,'(2(A,g10.4))') '# Epimin=', Epimin, '   Epimax=', Epimax

    n=1
    allocate (x(64*n))
    allocate (xsec_neut(64*n))
    allocate (xsec_prot(64*n))
    call sg64r(Epimin,Epimax,n,x,n1)


    Epi_loop:        do j=1,n1
       Epi=x(j)

       k_in=(/ Enu,0.,0.,Enu /)
       p_in=(/ mN,0.,0.,0. /)
       position=(/0.,0.,0./)

       call Nieves1piN_free_Epi(process_ID, ml_out, k_in, Epi, mN, p_in, position, 0, & ! 0 is neutron charge
            & 0, 1, mpi, & ! 0 is neutron charge,  1 is pion charge
            & xsec_neut(j)) ! OUT

       call Nieves1piN_free_Epi(process_ID, ml_out, k_in, Epi, mN, p_in, position, 1, & ! 0 is proton charge
            & 1, 1, mpi, & ! 0 is proton charge,  1 is pion charge
            & xsec_prot(j)) ! OUT

    end do Epi_loop

    call rg64r(Epimin,Epimax,n,xsec_neut,sig_neut)
    call rg64r(Epimin,Epimax,n,xsec_prot,sig_prot)

    write(10,*) '# Integrated over Epi:   Enu, sig_neut, sig_prot'
    write(10,'(A,3g12.5)') '#', Enu,  sig_neut, sig_prot

    deallocate(x,xsec_neut,xsec_prot)




  end subroutine Nieves_Epi





  ! provides output for 1-fold differential xsec dsi/dE1 (after  integrating the 2-fold one)  versus E1
  subroutine Nieves_E1(process_ID,ml_out,Enu,debugflag)

    use particleProperties, only : baryon, meson
    use idTable, only : nucleon, pion, Delta
    use random, only : rn
    use gauss_integration
    use lepton_kinematics_free
    use lepton_kinematics_free_FULL




    implicit none
    integer, intent(in) :: process_ID
    real, intent(in) :: ml_out, Enu
    logical, optional, intent(in) :: debugflag



    real :: mN, mDelta, mpi
    real :: E1, numin, numax, E1min, E1max
    real :: Epimin, Epimax, Epi
    integer :: j, k ! counters
    integer :: n1, n, nn1, nn ! number of integration points
    real, dimension(0:3) :: k_in, k_out, q, p_in, p_out, W, ppi_out
    real, dimension(1:3) :: position, vecW
    real :: dsig_neut, dsig_prot ! for one differential
    logical :: debug
    real, allocatable, dimension(:) :: x, xsec_neut, xsec_prot, xx, yy1, yy2

    if (present(debugflag)) then
       debug=debugflag
    else
       debug=.true.
    end if

    write(10,*) '# Differential x-sec: dsi/E1:   Enu, E1, neut, prot'
    write(10,'(1(A,g10.4))') '# Enu=', Enu

    mN=baryon(nucleon)%mass
    mDelta=baryon(Delta)%mass
    mpi=meson(pion)%mass


    call  minmaxNu_Enu(Enu,(mN+mpi),mN,ml_out,numin,numax) ! mN+mpi is Wmin
    E1min=Enu-numax
    E1max=Enu-numin
    if (debug) write(10,'(2(A,g10.4))') '# E1min=', E1min, '  E1max=', E1max

    n=1
    allocate (x(64*n))
    allocate (xsec_neut(64*n))
    allocate (xsec_prot(64*n))
    call sg64r(E1min,E1max,n,x,n1)


    energyloop:        do j=1,n1
       E1=x(j)
       k_in=(/ Enu,0.,0.,Enu /)
       p_in=(/ mN,0.,0.,0. /)
       position=(/0.,0.,0./)


       call Nieves1piN_free_E1(process_ID, ml_out, k_in, E1, mN, p_in, position, 0, &
            & 0, 0, mpi, &
            & xsec_neut(j))

       call Nieves1piN_free_E1(process_ID, ml_out, k_in, E1, mN, p_in, position, 1, &
            & 1, 1, mpi, &
            & xsec_prot(j))

       write(10,'(4g14.5)') Enu, E1, xsec_neut, xsec_prot

    end do energyloop

    call rg64r(E1min,E1max,n,xsec_neut,dsig_neut)
    call rg64r(E1min,E1max,n,xsec_prot,dsig_prot)

    write(10,*) '# Integrated over E1:   Enu, dsig_neut, dsig_prot'
    write(10,'(A,3g12.5)') '#', Enu,  dsig_neut, dsig_prot

    deallocate(x,xsec_neut,xsec_prot)

  end subroutine Nieves_E1







  ! provides output for 2-fold differential xsec dsi/dE1dcostheta  versus E1 for given Enu, costheta , obtained by integrating over costheta_Wpi
  subroutine Nieves_ct_E1_(process_ID,ml_out,Enu,costheta,debugflag)

    use particleProperties, only : baryon, meson
    use idTable, only : nucleon, pion, Delta
    use random, only : rn
    use gauss_integration
    use lepton_kinematics_free
    use lepton_kinematics_free_FULL



    implicit none
    integer, intent(in) :: process_ID
    real, intent(in) :: ml_out, Enu, costheta
    logical, optional, intent(in) :: debugflag



    real :: mN, mDelta, mpi
    real :: E1, k1, sintheta, numin, numax, E1min, E1max
    integer :: j, k ! counters
    integer :: n1, n, nn1, nn ! number of integration points
    real, dimension(0:3) :: k_in, k_out, q, p_in, p_out, W, ppi_out
    real, dimension(1:3) :: position, vecW
    real :: dsig_neut, dsig_prot ! for one differential
    logical :: debug
    real, allocatable, dimension(:) :: x, xsec_neut, xsec_prot

    if (present(debugflag)) then
       debug=debugflag
    else
       debug=.true.
    end if

    write(10,*) '# Double differential x-sec: dsi/E1/dcostheta:   Enu, costheta, E1, neut, prot'
    write(10,*) '# obtained by integrating over costheta_Wpi'
    write(10,'(2(A,g10.4))') '# Enu=', Enu, '   costheta=', costheta

    mN=baryon(nucleon)%mass
    mDelta=baryon(Delta)%mass
    mpi=meson(pion)%mass

    sintheta=sqrt(1.-costheta**2)

    call  minmaxNu_Enu_ct(Enu,costheta,(mN+mpi),mN,ml_out,numin,numax) ! 1.08 is Wmin
    E1min=Enu-numax
    E1max=Enu-numin
    if (debug) write(10,'(2(A,g10.4))') '# E1min=', E1min, '   E1max=', E1max

    n=1
    allocate (x(64*n))
    allocate (xsec_neut(64*n))
    allocate (xsec_prot(64*n))
    call sg64r(E1min,E1max,n,x,n1)


    energyloop:        do j=1,n1
       E1=x(j)
       k1=sqrt(E1**2-ml_out**2)

       k_in=(/ Enu, Enu  , 0. , 0. /)
       p_in=(/ mN,0.,0.,0. /)
       k_out=(/ E1,k1*costheta,k1*sintheta,0./)
       q=k_in-k_out
       W=p_in+q
       vecW=W(1:3)
       position=(/0.,0.,0./)

       call Nieves1piN_free_elepton_ct_(process_ID, ml_out, k_in, k_out,  mN, p_in, position, 0, & ! 0 is neutron charge
            & 0, 0, mpi, &  ! 0 is neutron charge, 0 is pion charge
            & xsec_neut(j))

       call Nieves1piN_free_elepton_ct_(process_ID, ml_out, k_in, k_out,  mN, p_in, position, 1, & ! 1 is proton charge
            & 1, 1, mpi, &  ! 1 is proton charge, 1 is pion charge
            & xsec_prot(j))

       write(10,'(5g14.5)') Enu, costheta, E1, xsec_neut(j), xsec_prot(j)

    end do energyloop

    call rg64r(E1min,E1max,n,xsec_neut,dsig_neut)
    call rg64r(E1min,E1max,n,xsec_prot,dsig_prot)

    write(10,*) '# Integrated over E1:   Enu, costheta, dsig_neut, dsig_prot'
    write(10,'(A,5g12.5)') '#', Enu, costheta,  dsig_neut, dsig_prot

    deallocate(x,xsec_neut,xsec_prot)

  end subroutine Nieves_ct_E1_







  ! provides output for 2-fold differential xsec dsi/dE1dcostheta  versus E1 for given Enu, costheta , obtained by integrating over Epi
  subroutine Nieves_ct_E1(process_ID,ml_out,Enu,costheta,debugflag)

    use particleProperties, only : baryon, meson
    use idTable, only : nucleon, pion, Delta
    use random, only : rn
    use gauss_integration
    use lepton_kinematics_free
    use lepton_kinematics_free_FULL

    implicit none
    integer, intent(in) :: process_ID
    real, intent(in) :: ml_out, Enu, costheta
    logical, optional, intent(in) :: debugflag



    real :: mN, mDelta, mpi
    real :: E1, k1, sintheta, numin, numax, E1min, E1max
    real :: Epimin, Epimax, Epi, costheta_PiW
    integer :: j ! counter
    real, dimension(0:3) :: k_in, k_out, q, p_in, p_out, W, ppi_out
    real, dimension(1:3) :: position, vecW
    real :: xsec_ntonpiplus, xsec_ptoppiplus
    logical :: debug

    if (present(debugflag)) then
       debug=debugflag
    else
       debug=.true.
    end if

    write(10,*) '# Double differential x-sec: dsi/E1/dcostheta:   Enu, costheta, E1, neut, prot'
    write(10,'(2(A,g10.4))') '# Enu=', Enu, '   costheta=', costheta

    mN=baryon(nucleon)%mass
    mDelta=baryon(Delta)%mass
    mpi=meson(pion)%mass


    sintheta=sqrt(1.-costheta**2)


    call  minmaxNu_Enu_ct(Enu,costheta,(mN+mpi),mN,ml_out,numin,numax) ! 1.08 is Wmin
    E1min=Enu-numax
    E1max=Enu-numin
    if (debug) write(10,'(2(A,g10.4))') '# E1min=', E1min, '   E1max=', E1max

    write(10,*) '# Enu, costheta, E1,  xsec_ntonpiplus, xsec_ptoppiplus'



    E1loop:        do j=1,200
       E1=0.01*j
       if (E1.le.E1min) cycle
       if (E1.ge.E1max) exit


       k1=sqrt(E1**2-ml_out**2)

       k_in=(/ Enu, 0. , 0. , Enu /)
       p_in=(/ mN,0.,0.,0. /)
       k_out=(/ E1, -k1*sintheta, 0., k1*costheta/)
       q=k_in-k_out
       W=p_in+q
       vecW=W(1:3)
       position=(/0.,0.,0./)

       call Nieves1piN_free_elepton_ct(process_ID, ml_out, k_in, k_out,  mN, p_in, position, 0, & ! 0 is neutron charge
            & 0, 1, mpi, & ! 0 is neutron charge, 1 is pion charge
            & xsec_ntonpiplus) ! OUT


       call Nieves1piN_free_elepton_ct(process_ID, ml_out, k_in, k_out,  mN, p_in, position, 1, & ! 0 is proton charge
            & 1, 1, mpi, & ! 1 is proton charge, 0 is pion charge
            & xsec_ptoppiplus) ! OUT


       write(10,'(5g14.5)') Enu, costheta, E1, xsec_ntonpiplus, xsec_ptoppiplus

    end do E1loop

  end subroutine Nieves_ct_E1







  ! provides output for 3-fold differential xsection dsi/dE1 dcostheta dctWpi   versus costheta_Wpi  for given Enu, E1, costheta_Wpi
  subroutine Nieves_E1_ct_ctWpi(process_ID,ml_out,Enu,E1,costheta,debugflag)

    use particleProperties, only : baryon, meson
    use idTable, only : nucleon, pion, Delta
    use minkowski, only : abs4, SP
    use singlePionProductionNHVlike
    use constants, only : pi
    use rotation
    use random, only : rn
    use gauss_integration
    use lepton_kinematics_free
    use lepton_kinematics_free_FULL

    implicit none
    integer, intent(in) :: process_ID
    real, intent(in) :: ml_out, Enu, E1, costheta
    logical, intent(in), optional :: debugflag

    logical :: debug
    real :: mN, mDelta, mpi
    real :: k1, sintheta, numin, numax, E1min, E1max
    real :: Epi, ctWpi, ctWpi_min, ctWpi_max
    integer :: j, n1, n
    real, dimension(0:3) :: k_in, k_out, q, p_in, p_out, W, ppi_out
    real, dimension(1:3) :: position, vecW, vecppi_W, vecppi_out
    ! real :: thetaW, phiW  ! to be calculated
    real :: dsig_neut, dsig_prot, kinem_factor
    real, allocatable, dimension(:) :: x,y1,y2
    complex :: MaEl_prot, MaEl_neut


    if (present(debugflag)) then
       debug=debugflag
    else
       debug=.true.
    end if


    write(10,*) '# Triple differential x-sec: dsi/E1/dcostheta/dcostheta_PiW'
    write(10,*) '#Enu E1 costheta costheta_PiW Epi neut prot'


    mN=baryon(nucleon)%mass
    mDelta=baryon(Delta)%mass
    mpi=meson(pion)%mass

    !if (debug==.true.) write(10,'(3(A,g10.3))') '# mN=', mN,'    mDelta=', mDelta, '     mpi=', mpi

    k1=sqrt(E1**2-ml_out**2)
    if (debug==.true.) write(10,'(4(A,g10.3))') '# In Nieves_E1_ct_ctWpi: Enu=', Enu, &
         & '    costheta=', costheta, '    E1=', E1, '   k1=',k1
    sintheta=sqrt(1.-costheta**2)


    ! to check if we are within the limits for E1
    call  minmaxNu_Enu_ct(Enu,costheta,1.08,mN,ml_out,numin,numax) ! 1.08 is Wmin
    E1min=Enu-numax
    E1max=Enu-numin
    if (debug==.true.) write(10,'(2(A,g10.3))') '# In Nieves_E1_ct_ctWpi: To check the limits:   E1min=', E1min,'    E1max=', E1max



    if(E1<E1min .or. E1>E1max ) then
       write(*,*) '# In Nieves_E1_ct_ctWpi:  Enu=',Enu,'     costheta=',costheta
       write(*,*) '# E1=', E1, 'is out of kimatically allowed range E1min=', E1min, '   E1max=',E1max
       return
    end if


    k_in=(/ Enu,Enu,0.,0. /)
    p_in=(/ mN,0.,0.,0. /)
    k_out=(/ E1,k1*costheta,0.,-k1*sintheta/)
    q=k_in-k_out
    write(10,*) '# In Nieves_E1_ct_ctWpi:    Q2:', -SP(q,q)
    W=p_in+q
    !write(*,*) '# In Nieves_E1_ct_ctWpi:    Q2:', -SP(q,q), '    W=', W
    vecW=W(1:3)
    position=(/0.,0.,0./)


    ! generate an array of cosThetaWPi
    call minmaxCosThetaWPi_W(W,mN,mpi,ctWpi_min,ctWpi_max)
    if (debug==.true.) write(10,'(2(A,g10.3))') '# costheta_Wpi_min=', ctWpi_min,'    costheta_Wpi_max=', ctWpi_max
    n=1
    allocate (x(64*n))
    allocate (y1(64*n))
    allocate (y2(64*n))
    call sg64r(ctWpi_min,ctWpi_max,n,x,n1)


    energyPionLoop: do j=1,n1
       ctWpi=x(j)
       !! if (debug==.true.) write(10,'(A,g10.3)') '#  ppi=', ppi

       call Nieves1piN_free_elepton_ct_ctWpi(process_ID, ml_out, k_in, k_out,  mN, p_in, position, 0, & ! 0 is incoming neutron charge
            & 0, 1, mpi, ctWpi, &  ! 0 is neutron  1 is pion charges
            & Epi, p_out, ppi_out, y1(j))

       call Nieves1piN_free_elepton_ct_ctWpi(process_ID, ml_out, k_in, k_out,  mN, p_in, position, 1, & ! 0 is incoming proton charge
            & 1, 1, mpi, ctWpi, &  ! 1 is proton  1 is pion charges
            & Epi, p_out, ppi_out, y2(j))

       write(10,'(7g12.5)') Enu, E1, costheta, ctWpi, Epi, y1(j), y2(j)

    end do energyPionLoop

    ! integrate over Epi
    call rg64r(ctWpi_min,ctWpi_max,n,y1,dsig_neut)
    call rg64r(ctWpi_min,ctWpi_max,n,y2,dsig_prot)

    write(10,*) '# Integrated over ctWpi:   Enu, costheta, E1, dsig_neut, dsig_prot'
    write(10,'(A,5g12.5)') '#', Enu, costheta, E1, dsig_neut, dsig_prot

    deallocate(x,y1,y2)

  end subroutine Nieves_E1_ct_ctWpi


  ! provides output for 3-fold differential xsection dsi/dE1 dcostheta dEpi   versus Epi  for given Enu, E1, costheta
  subroutine Nieves_E1_ct_Epi(process_ID,ml_out,Enu,E1,costheta,debugflag)

    use particleProperties, only : baryon, meson
    use idTable, only : nucleon, pion, Delta
    use minkowski, only : abs4, SP
    use singlePionProductionNHVlike
    use constants, only : pi
    use rotation
    use random, only : rn
    use gauss_integration
    use lepton_kinematics_free
    use lepton_kinematics_free_FULL

    implicit none
    integer, intent(in) :: process_ID
    real, intent(in) :: ml_out, Enu, E1, costheta
    logical, intent(in), optional :: debugflag

    logical :: debug
    real :: mN, mDelta, mpi
    real :: k1, sintheta, numin, numax, E1min, E1max
    real :: Epimin, Epimax, Epi, ppi, costheta_PiW!, sintheta_PiW, phi_PiW
    integer :: j, n1, n
    real, dimension(0:3) :: k_in, k_out, q, p_in, p_out, W, ppi_out
    real, dimension(1:3) :: position, vecW, vecppi_W, vecppi_out
    ! real :: thetaW, phiW ! to be calculated
    real :: dsig_neut, dsig_prot, kinem_factor
    real, allocatable, dimension(:) :: x,y1,y2
    complex :: MaEl_prot, MaEl_neut


    if (present(debugflag)) then
       debug=debugflag
    else
       debug=.true.
    end if


    write(10,*) '# Triple differential x-sec: dsi/E1/dcostheta/dEpi'
    write(10,*) '#Enu E1 costheta Epi costheta_PiW neut prot'


    mN=baryon(nucleon)%mass
    mDelta=baryon(Delta)%mass
    mpi=meson(pion)%mass

    if (debug==.true.) write(10,'(3(A,g10.3))') '# mN=', mN,'    mDelta=', mDelta, '     mpi=', mpi

    k1=sqrt(E1**2-ml_out**2)
    if (debug==.true.) write(10,'(4(A,g10.3))') '# Enu=', Enu,'    costheta=', costheta, '    E1=', E1, '   k1=',k1
    sintheta=sqrt(1.-costheta**2)


    call  minmaxNu_Enu_ct(Enu,costheta,1.08,mN,ml_out,numin,numax) ! 1.08 is Wmin
    E1min=Enu-numax
    E1max=Enu-numin
    if (debug==.true.) write(10,'(2(A,g10.3))') '#To check the limits:   E1min=', E1min,'    E1max=', E1max



    if(E1<E1min .or. E1>E1max ) then
       write(*,*) '# In Nieves_E1_ct_Epi  Enu=',Enu,'     costheta=',costheta
       write(*,*) '# E1=', E1, 'is out of kimatically allowed range E1min=', E1min, '   E1max=',E1max
       return
    end if


    k_in=(/ Enu,Enu,0.,0. /)
    p_in=(/ mN,0.,0.,0. /)
    k_out=(/ E1,k1*costheta,0.,-k1*sintheta/)
    q=k_in-k_out
    write(10,*) '# In test_bgr:    Q2:', -SP(q,q)
    W=p_in+q
    !write(*,*) '# In test_bgr:    Q2:', -SP(q,q), '    W=', W
    vecW=W(1:3)
    position=(/0.,0.,0./)


    ! generate an array of Epi
    call minmaxEpi_W(W,mN,mpi,Epimin,Epimax)
    if (debug==.true.) write(10,'(2(A,g10.3))') '#Epi_min=', Epimin,'    Epi_max=', Epimax

    n=1
    allocate (x(64*n))
    allocate (y1(64*n))
    allocate (y2(64*n))
    call sg64r(Epimin,Epimax,n,x,n1)


    energyPionLoop: do j=1,n1
       Epi=x(j)
       ppi=sqrt(Epi**2-mpi**2)
       !! if (debug==.true.) write(10,'(A,g10.3)') '#  ppi=', ppi

       call Nieves1piN_free_elepton_ct_Epi(process_ID, ml_out, k_in, k_out,  mN, p_in, position, 0, & ! 0 is incoming neutron charge
            & 0, 1, mpi, Epi, &  ! 0 is neutron  1 is pion charges
            & costheta_PiW, p_out, ppi_out, y1(j))

       call Nieves1piN_free_elepton_ct_Epi(process_ID, ml_out, k_in, k_out,  mN, p_in, position, 1, & ! 0 is incoming proton charge
            & 1, 1, mpi, Epi, &  ! 1 is proton  1 is pion charges
            & costheta_PiW, p_out, ppi_out, y2(j))

       write(10,'(7g12.5)') Enu, E1, costheta, Epi, costheta_PiW, y1(j), y2(j)

    end do energyPionLoop

    ! integrate over Epi
    call rg64r(Epimin,Epimax,n,y1,dsig_neut)
    call rg64r(Epimin,Epimax,n,y2,dsig_prot)

    write(10,*) '# Integrated over Epi:   Enu, costheta, E1, dsig_neut, dsig_prot'
    write(10,'(A,5g12.5)') '#', Enu, costheta, E1, dsig_neut, dsig_prot

    deallocate(x,y1,y2)

  end subroutine Nieves_E1_ct_Epi







  ! calculates 3-fold differential xsection
  subroutine Nieves1piN_free_elepton_ct_Epi(process_ID, ml, k_in, k_out,  mN, p_in, position, charge_in, &
       & charge_out, pion_charge_out, mpi, Epi, &
       & costheta_PiW, p_out, ppi_out, xsec) ! OUT

    use singlePionProductionNHVlike
    use constants, only : pi
    use rotation
    use random, only : rn
    use lepton_kinematics_free
    use lepton_kinematics_free_FULL
    use gauss_integration


    implicit none
    integer, intent(in) :: process_ID !CC, NC, em
    real, dimension(0:3), intent(in) :: k_in, k_out, p_in
    integer, intent(in) :: charge_in, charge_out, pion_charge_out
    real, intent(in) :: ml, mN, mpi, Epi
    real, dimension(1:3), intent(in) :: position
    real, dimension(0:3), intent(out) :: p_out, ppi_out
    real, intent(out) :: costheta_PiW, xsec


    real, dimension(0:3) :: W
    real :: Enu, k1, ppi, sintheta_PiW, phi_PiW
    real, dimension(1:3) :: vecW, vecppi_W, vecppi_out
    complex :: MaEl
    real :: thetaW, phiW, kinem_factor!, kinem_factor2

    real, allocatable, dimension(:) :: xx, yy
    integer :: nn, nn1, k


    W=k_in-k_out+p_in
    vecW=W(1:3)
    ppi=sqrt(Epi**2-mpi**2)
    Enu=k_in(0)
    k1=sqrt( dot_product(k_out(1:3),k_out(1:3)) )
    !! if (debug) write(10,'(2(A,g10.3))') '#     Epi=', Epi, '     ppi=', ppi


    ! calculate costheta_PiW from kinematics of the off-mass-shell Delta decay
    call CosThetaPiW_W_Epi(W,Epi,mN,mpi,costheta_PiW) ! costheta between vec{W} and vec{ppi}
    sintheta_PiW=sqrt(1.-costheta_PiW**2)

    nn=1
    allocate (xx(64*nn))
    allocate (yy(64*nn))

    ! integration over phi_PiW (azimuthal angle in W-along-z frame) here , start loop
    ! is not 2pi because matrix element can depend on PhiWPi
    call sg64r(0.,2.*pi,nn,xx,nn1)
    over_PhiWPi: do k=1,nn1
       phi_PiW=xx(k)

       ! set up ppi in the frame that W is along z
       vecppi_W=(/ppi*sintheta_PiW*cos(phi_PiW), ppi*sintheta_PiW*sin(phi_PiW), ppi*costheta_PiW /)

       ! determine the angle of W with respect to z-axis
       call get_phi_Theta(vecW, thetaW, phiW)

       ! rotate pi on this angle
       call rotateYZ( thetaW, phiW, vecppi_W, vecppi_out)

       ! checked 2009-02-20 OK
       !write(*,*) ' |vecppi_W|^2=  ', DOT_PRODUCT(vecppi_W,vecppi_W), ' must equals to   |vecppi_out|^2=',  DOT_PRODUCT(vecppi_out,vecppi_out)

       ! set up the final  ppion_out and p_out
       ppi_out=(/ Epi, vecppi_out(1:3) /)
       p_out=W-ppi_out
       !!if (debug==.true.) write(10,'(A,g10.3)') '# Checking ppi_out,  abs4(ppi_out)=', abs4(ppi_out)
       !!if (debug==.true.) write(10,'(A,g10.3)') '#     |ppi_out|=', sqrt(Dot_product(ppi_out(1:3),ppi_out(1:3)))

       ! matrix element is invariant
       MaEl=Nieves_MaEl_1pi(process_ID,  k_in, k_out,  p_in, position, charge_in, &
            & p_out, charge_out,  ppi_out, pion_charge_out)


       ! absolute value of W does not depend on W-direction
       kinem_factor = k1/4./mN/Enu/(2*pi)**4/8./sqrt( dot_product(vecW,vecW) )

       ! not relevant here; used for for inttegration over Epi for costheta_Wpi distribution
       !kinem_factor2 = k1/4./mN/Enu/(2.*pi)**4*k1*ppi**2/8./abs( sqrt(dot_product(vecW,vecW))*costheta_PiW*Epi - W(0)*ppi )
       !write(*,'(3(A,g11.4))') '# MatEl=',REAL(MaEl), '    Kinem factror 1 =', kinem_factor, '    Kin factor 2=', kinem_factor2

       yy(k) = REAL(MaEl)*kinem_factor*3.8938E-28

    end do over_PhiWPi

    call rg64r(0.,2.*pi,nn,yy,xsec) ! the integration itself

    ! write(10,'(A,g14.4)') '# integrated with the result    2-bl xsec=',xsec
    deallocate(xx,yy)

  end subroutine Nieves1piN_free_elepton_ct_Epi




  ! calculates 3-fold differential xsection  dsi/dE1 dcostheta dcostheta_Wpi
  subroutine Nieves1piN_free_elepton_ct_ctWpi(process_ID, ml, k_in, k_out,  mN, p_in, position, charge_in, &
       & charge_out, pion_charge_out, mpi, costheta_PiW, &
       & Epi, p_out, ppi_out, xsec)

    use singlePionProductionNHVlike
    use constants, only : pi
    use rotation
    use random, only : rn
    use lepton_kinematics_free
    use lepton_kinematics_free_FULL



    implicit none
    integer, intent(in) :: process_ID !CC, NC, em
    real, dimension(0:3), intent(in) :: k_in, k_out, p_in
    integer, intent(in) :: charge_in, charge_out, pion_charge_out
    real, intent(in) :: ml, mN, mpi, costheta_PiW
    real, dimension(1:3), intent(in) :: position
    real, dimension(0:3), intent(out) :: p_out, ppi_out
    real, intent(out) :: Epi, xsec


    real, dimension(0:3) :: W
    real :: Enu, k1, ppi, sintheta_PiW, phi_PiW
    real, dimension(1:3) :: vecW, vecppi_W, vecppi_out
    complex :: MaEl
    real :: thetaW, phiW, kinem_factor

    ! for the case that there are two solutions for Epi
    integer :: numRoots
    real :: Epi1,Epi2


    W=k_in-k_out+p_in
    vecW=W(1:3)
    ppi=sqrt(Epi**2-mpi**2)
    Enu=k_in(0)
    k1=sqrt( dot_product(k_out(1:3),k_out(1:3)) )
    !! if (debug==.true.) write(10,'(2(A,g10.3))') '#     Epi=', Epi, '     ppi=', ppi


    ! calculate Epi from kinematics of the off-mass-shell Delta decay
    Epi=PionEnergy_W_ctWpi(W,costheta_PiW,mN,mpi,numRoots,Epi1,Epi2)
    ppi=sqrt(Epi**2-mpi**2)
    sintheta_PiW=sqrt(1.-costheta_PiW**2)
    ! generate phi_PiW randomly
    phi_PiW=pi/2. !rn()*2*pi
    ! set up ppi in the frame that W is along z
    vecppi_W=(/ppi*sintheta_PiW*cos(phi_PiW), ppi*sintheta_PiW*sin(phi_PiW), ppi*costheta_PiW /)

    ! determine the angle of W with respect to z-axis
    call get_phi_Theta(vecW, thetaW, phiW)

    ! rotate pi on this angle
    call rotateYZ( thetaW, phiW, vecppi_W, vecppi_out)

    ! set up the final  ppion_out and p_out
    ppi_out=(/ Epi, vecppi_out(1:3) /)
    p_out=W-ppi_out


    MaEl=Nieves_MaEl_1pi(process_ID,  k_in, k_out,  p_in, position, charge_in, &
         & p_out, charge_out,  ppi_out, pion_charge_out)


    kinem_factor = k1/4./mN/Enu/(2.*pi)**3*k1*ppi**2/8./abs( sqrt(dot_product(vecW,vecW))*costheta_PiW*Epi - W(0)*ppi )

    xsec = REAL(MaEl)*kinem_factor*3.8938E-28

  end subroutine Nieves1piN_free_elepton_ct_ctWpi


  ! calculates 2-fold differential xsection dsi/ dE1 dcostheta by integrating over Epi
  subroutine Nieves1piN_free_elepton_ct(process_ID, ml, k_in, k_out,  mN, p_in, position, charge_in, &
       & charge_out, pion_charge_out, mpi, &
       & xsec) ! OUT


    use lepton_kinematics_free
    use lepton_kinematics_free_FULL

    use gauss_integration


    implicit none
    integer, intent(in) :: process_ID !CC, NC, em
    real, dimension(0:3), intent(in) :: k_in, k_out, p_in
    integer, intent(in) :: charge_in, charge_out, pion_charge_out
    real, intent(in) :: ml, mN, mpi ! mpi is needed for Epi-integration-limits
    real, dimension(1:3), intent(in) :: position
    real, intent(out) :: xsec


    real, dimension(0:3) :: W, p_out, ppi_out
    real :: Epimin, Epimax, Epi, costheta_PiW
    real, allocatable, dimension(:) :: xx, yy
    integer :: nn, nn1, k

    W=k_in-k_out+p_in

    call minmaxEpi_W(W,mN,mpi,Epimin,Epimax)
    if (debugflag) write(10,'(4(A,g10.4))') '#      Epimin=', Epimin, '   Epimax=', Epimax

    ! integration over Epi
    nn=1
    allocate (xx(64*nn))
    allocate (yy(64*nn))

    call sg64r(Epimin,Epimax,nn,xx,nn1)
    integration_loop: do k=1,nn1
       Epi=xx(k)

       call Nieves1piN_free_elepton_ct_Epi(process_ID, ml, k_in, k_out,  mN, p_in, position, charge_in, &
            & charge_out, pion_charge_out, mpi, Epi, &
            & costheta_PiW, p_out, ppi_out, yy(k)) ! output



    end do integration_loop
    call rg64r(Epimin,Epimax,nn,yy,xsec) ! the integration itself

    deallocate(xx,yy)
    ! integration over Epi finished

  end subroutine Nieves1piN_free_elepton_ct








  ! calculates 2-fold differential xsection dsi/ dE1 dcostheta by integrating over costheta_WPi
  subroutine Nieves1piN_free_elepton_ct_(process_ID, ml, k_in, k_out,  mN, p_in, position, charge_in, &
       & charge_out, pion_charge_out, mpi, &
       & xsec)
    ! to be corrected --- we cannot integrate over costheta without taking into account that there are two physical solutions for Epi sometimes

    use lepton_kinematics_free
    use lepton_kinematics_free_FULL

    use gauss_integration


    implicit none
    integer, intent(in) :: process_ID !CC, NC, em
    real, dimension(0:3), intent(in) :: k_in, k_out, p_in
    integer, intent(in) :: charge_in, charge_out, pion_charge_out
    real, intent(in) :: ml, mN, mpi ! mpi is needed for Epi-integration-limits
    real, dimension(1:3), intent(in) :: position
    real, intent(out) :: xsec


    real, dimension(0:3) :: W, p_out, ppi_out
    real :: ctWpi_min, ctWpi_max, ctWpi, Epi
    real, allocatable, dimension(:) :: xx, yy
    integer :: nn, nn1, k

    W=k_in-k_out+p_in

    ! integration limits over costheta_PiW

    call minmaxCosThetaWPi_W(W,mN,mpi,ctWpi_min,ctWpi_max)

    ! integration over Epi
    nn=1
    allocate (xx(64*nn))
    allocate (yy(64*nn))

    call sg64r(ctWpi_min,ctWpi_max,nn,xx,nn1)
    integration_loop: do k=1,nn1
       ctWpi=xx(k)

       call Nieves1piN_free_elepton_ct_ctWpi(process_ID, ml, k_in, k_out,  mN, p_in, position, charge_in, &
            & charge_out, pion_charge_out, mpi, ctWpi, &
            & Epi, p_out, ppi_out, yy(k)) ! output

       !    write(10,'(A,I4.3,3(A,g14.4))') '# k=', k, '    ctWpi=',ctWpi,  '   Epi=',Epi,  '    3-pl xsec=',yy(k)

    end do integration_loop

    call rg64r(ctWpi_min,ctWpi_max,nn,yy,xsec) ! the integration itself

    ! write(10,'(A,g14.4)') '# integrated with the result    2-bl xsec=',xsec
    deallocate(xx,yy)
    ! integration over costheta_Wpi finished

  end subroutine Nieves1piN_free_elepton_ct_


  ! calculates 2-fold differential xsection dsi/ dE1 dEpi
  subroutine Nieves1piN_free_elepton_Epi(process_ID, ml, k_in, E1, Epi,  mN, p_in, position, charge_in, &
       & charge_out, pion_charge_out, mpi, &
       & xsec) ! OUT


    use lepton_kinematics_free
    use lepton_kinematics_free_FULL

    use gauss_integration
    use rotation


    implicit none
    integer, intent(in) :: process_ID !CC, NC, em
    real, dimension(0:3), intent(in) :: k_in, p_in
    integer, intent(in) :: charge_in, charge_out, pion_charge_out
    real, intent(in) :: E1, Epi, ml, mN, mpi ! mpi is needed for
    real, dimension(1:3), intent(in) :: position
    real, intent(out) :: xsec


    real, dimension(0:3) :: k_out, W, p_out, ppi_out
    real, dimension(1:3) :: veckout_k, veck
    real :: costhetamin, costhetamax, costheta, sintheta, k1, costheta_PiW
    real :: thetak, phik
    real, allocatable, dimension(:) :: xx, yy
    integer :: nn, nn1, k


    ! determine the angles of k with respect to z-axis
    call get_phi_Theta(veck, thetak, phik)
    ! the absolute value of k1
    k1=sqrt((E1+ml)*(E1-ml))
    ! zero-th komponent of k_out
    k_out(0)=E1

    call minmaxcostheta_Enu_nu(Enu,(Enu-E1),(mN+mpi),mN,ml,costhetamin,costhetamax)
    !if (debug) write(10,'(3(A,g10.4))') '# E1=',E1, '     Epimin=', Epimin, '   Epimax=', Epimax

    ! integration over costheta
    nn=1
    allocate (xx(64*nn))
    allocate (yy(64*nn))

    call sg64r(costhetamin,costhetamax,nn,xx,nn1)
    integration_loop: do k=1,nn1
       costheta=xx(k)
       sintheta=sqrt((1.-costheta)*(1.+costheta))

       ! set up veckout in the frame that k_in  is along z
       veckout_k=(/k1*sintheta, 0., k1*costheta /)

       ! rotate veckout on this angle, thus setting up spatial components of k_out
       call rotateYZ( thetak, phik, veckout_k, k_out(1:3))

       call Nieves1piN_free_elepton_ct_Epi(process_ID, ml, k_in, k_out,  mN, p_in, position, charge_in, &
            & charge_out, pion_charge_out, mpi, Epi, &
            & costheta_PiW, p_out, ppi_out, yy(k)) ! output

    end do integration_loop
    call rg64r(costhetamin,costhetamax,nn,yy,xsec) ! the integration itself

    deallocate(xx,yy)
    ! integration over costheta finished

  end subroutine Nieves1piN_free_elepton_Epi


  ! calculates 1-fold differential xsection  dsi/dEpi
  subroutine Nieves1piN_free_Epi(process_ID, ml, k_in, Epi, mN, p_in, position, charge_in, &
       & charge_out, pion_charge_out, mpi, &
       & xsec) ! OUT
    use lepton_kinematics_free
    use lepton_kinematics_free_FULL

    use gauss_integration


    implicit none
    integer, intent(in) :: process_ID !CC, NC, em
    real, dimension(0:3), intent(in) :: k_in, p_in
    integer, intent(in) :: charge_in, charge_out, pion_charge_out
    real, intent(in) :: Epi, ml, mN, mpi ! mpi is needed for
    real, dimension(1:3), intent(in) :: position
    real, intent(out) :: xsec

    real :: E1min, E1max
    real, allocatable, dimension(:) :: xx, yy
    integer :: nn, nn1, k



    call minmaxE1_Enu_Epi(Enu,Epi,(mpi+mN),ml,mN,mpi,E1min,E1max)
    write(10,'(3(A,g10.4))') '# Enu=',Enu, '   Epi=',Epi, '     E1min=', E1min, '   E1max=', E1max

    ! integration over E1
    nn=1
    allocate (xx(64*nn))
    allocate (yy(64*nn))

    call sg64r(E1min,E1max,nn,xx,nn1)
    integration_loop: do k=1,nn1
       E1=xx(k)

       call Nieves1piN_free_elepton_Epi(process_ID, ml, k_in, E1, Epi, mN, p_in, position, charge_in, &
            & charge_out, pion_charge_out, mpi, &
            & yy(k) ) ! output

    end do integration_loop
    call rg64r(E1min,E1max,nn,yy,xsec) ! the integration itself

    deallocate(xx,yy)
    ! integration over costheta finished

  end subroutine Nieves1piN_free_Epi


  ! calculates 1-fold differential xsection  dsi/dE1 by integrating over costheta
  subroutine  Nieves1piN_free_E1(process_ID, ml, k_in, E1, mN, p_in, position, charge_in, &
       & charge_out, pion_charge_out, mpi, &
       & xsec) ! OUT
    use lepton_kinematics_free
    use lepton_kinematics_free_FULL

    use gauss_integration


    implicit none
    integer, intent(in) :: process_ID !CC, NC, em
    real, dimension(0:3), intent(in) :: k_in, p_in
    integer, intent(in) :: charge_in, charge_out, pion_charge_out
    real, intent(in) :: E1, ml, mN, mpi
    real, dimension(1:3), intent(in) :: position
    real, intent(out) :: xsec

    real :: k1, ctmin, ctmax, costheta, sintheta
    real, dimension(0:3) :: k_out
    real, allocatable, dimension(:) :: xx, yy
    integer :: nn, nn1, k


    k1=sqrt(E1**2-ml**2)

    call minmaxcostheta_Enu_nu(Enu,(Enu-E1),(mN+mpi),mN,ml,ctmin,ctmax)
    !if (debug) write(10,'(4(A,g10.4))') '# Enu=',Enu, '    E1=',E1, '     ctmin=', ctmin, '   ctmax=', ctmax

    ! integration over costheta
    nn=1
    allocate (xx(64*nn))
    allocate (yy(64*nn))

    call sg64r(ctmin,ctmax,nn,xx,nn1)
    integration_loop: do k=1,nn1
       costheta=xx(k)
       sintheta=sqrt(1.-costheta**2)

       k_out=(/ E1,k1*sintheta,0,k1*costheta/)

       call Nieves1piN_free_elepton_ct(process_ID, ml, k_in, k_out,  mN, p_in, position, charge_in, &  ! 0 is neutron charge
            & charge_out, pion_charge_out, mpi, &     ! 0 is neutron charge, 0 is pion charge
            & yy(k))

    end do integration_loop

    call rg64r(ctmin,ctmax,nn,yy,xsec) ! the integration itself

    !write(10,'(5g12.5)') Enu, E1, xsec

    deallocate(xx,yy)
    ! integration over ct finished

  end subroutine Nieves1piN_free_E1

























  ! ***************************************************
  ! output for free lepton kinematics
  ! ***************************************************




  subroutine show_minmaxt_Q2_W2()

    implicit none
    real :: mN,Q2,W2,mN1,mpi,tmin,tmax

    name='show_minmaxt_Q2_W2.dat'
    Open(10,File=name)


    mN=0.939
    do k=0,4
       Q2=0.5+0.5*k
       mN1=mN
       mpi=0.14

       write(10,'(5(A,g12.5))') '# In show_minmaxt_Q2_W2    mN=', mN, '   Q2=',Q2, &
            & '    mN1=', mN1, '    mpi=',mpi
       write(10,'(5(A,g12.5))') '#  Q2 W2 tmin tmax'

       do i=1,600
          W2=(1.070+0.002*i)**2
          call minmaxt_Q2_W2(mN,Q2,W2,mN1,mpi,tmin,tmax)
          write(10,'(6(g16.5))') Q2, W2, tmin, tmax
       end do
       write(10,*) ''
       write(10,*) ''
       write(10,*) ''
    end do
    close(10)
  end subroutine show_minmaxt_Q2_W2




  subroutine show_minmaxctqPi_Q2_W2_Epi()

    implicit none
    real :: mN,Q2,W2,Epi,mN1,mpi,ctmin,ctmax

    name='show_minmaxctqPi_Q2_W2_Epi.dat'
    Open(10,File=name)


    mN=0.939
    do k=0,4
       Q2=0.5+0.5*k
       W2=1.6**2
       mN1=mN
       mpi=0.14

       write(10,'(5(A,g12.5))') '# In show_minmaxctqPi_Q2_W2_Epi mN=', mN, '   Q2=',Q2, '   W2=',W2, &
            & '    mN1=',mN1, '   mpi=',mpi
       write(10,'(5(A,g12.5))') '#  Q2 W2 Epi ctmin ctmax'

       do i=1,260
          Epi=mpi+0.02*i
          call minmaxcosthetaqPi_Q2_W2_Epi(mN,Q2,W2,Epi,mN1,mpi,ctmin,ctmax)
          write(10,'(6(g16.5))') Q2, W2, Epi, ctmin, ctmax
       end do
       write(10,*) ''
       write(10,*) ''
       write(10,*) ''
    end do
    close(10)
  end subroutine show_minmaxctqPi_Q2_W2_Epi









  subroutine show_minmaxTpi_Enu_sl()

    implicit none
    real :: ml,mN,Enu,sl,mN1,mpi,Tpimin,Tpimax, slmin,slmax

    name='show_minmaxTpi_Enu_sl.dat'
    Open(10,File=name)


    ml=0.105
    mN=0.939
    do k=0,4
       Enu=0.5+0.5*k
       mN1=mN
       mpi=0.14

       call minmaxsl_Enu(ml,mN,Enu,mN1,mpi,slmin,slmax)

       write(10,'(7(A,g12.5))') '# In show_minmaxTpi_Enu_sl    mN=', mN, '   Enu=',Enu, &
            & '    mN1=', mN1, '    mpi=',mpi, '   slmin=',slmin, '   slmax=',slmax
       write(10,'(A)') '#  Enu sl Tpimin Tpimax'

       do i=1,200
          sl=slmin+(slmax-slmin)/200.*i
          call minmaxTpi_Enu_sl(mN,Enu,sl,mN1,mpi,Tpimin,Tpimax)
          write(10,'(6(g16.5))') Enu, sl, Tpimin, Tpimax
       end do
       write(10,*) ''
       write(10,*) ''
       write(10,*) ''
    end do
    close(10)
  end subroutine show_minmaxTpi_Enu_sl




  subroutine show_minmaxctnupi_Enu_sl_Epi()

    implicit none
    real :: mN,Enu,sl,mN1,mpi,Epi,ctmin,ctmax

    real :: slmin,slmax,ml

    name='show_minmaxctnupi_Enu_sl_Epi.dat'
    Open(10,File=name)

    ml=0.105
    mN=0.939
    mN1=mN
    mpi=0.14
    Enu=1.

    call minmaxsl_Enu(ml,mN,Enu,mN1,mpi,slmin,slmax)


    do k=0,5
       !sl= (mpi+0.105)**2 + ( (Enu-mpi)**2 - (mpi+0.105)**2 )/2. ! average between min and max
       sl=slmax*(0.98+0.004*k)*(1.-0.00001)

       write(10,'(5(A,g12.5))') '# In show_minmaxctnupi_Enu_sl_Epi    mN=', mN, '   Enu=',Enu, '   sl=',sl, &
            & '    mN1=',mN1, '   mpi=',mpi
       write(10,'(5(A,g12.5))') '#  Enu sl Epi ctmin ctmax'

       do i=1,260
          Epi=mpi+0.02*i
          call minmaxctnupi_Enu_sl_Epi(mN,Enu,sl,mN1,mpi,Epi,ctmin,ctmax)
          write(10,'(6(g16.5))') Enu, sl, Epi, ctmin, ctmax
       end do
       write(10,*) ''
       write(10,*) ''
       write(10,*) ''
    end do
    close(10)
  end subroutine show_minmaxctnupi_Enu_sl_Epi




  subroutine show_minmaxEpi_Enu()

    implicit none
    real :: Enu , mN , mpi , mN1 , ml , Epimin , Epimax

    name='show_minmaxEpi_Enu.dat'
    Open(10,File=name)


    ml=0.105
    mN=0.939
    mN1=mN
    mpi=0.14



    write(10,'(7(A,g12.5))') '# In show_minmaxEpi_Enu    mN=', mN, &
         & '    mN1=', mN1, '    mpi=',mpi
    write(10,'(A)') '#  Enu Epimin Epimax'

    do i=1,200
       Enu=0.5+0.01*i
       call minmaxEpi_Enu( Enu , mN , mpi , mN1 , ml , Epimin , Epimax )
       write(10,'(3(g16.5))') Enu,  Epimin, Epimax
    end do
    close(10)
  end subroutine show_minmaxEpi_Enu









  ! **********************************************************************************
  ! output for cNp, CT, pp, pF diagrams
  !
  ! the functions used for the actual calculation of the xsec are in lepton_xsec_free.f90
  ! **********************************************************************************


  ! provides output for 2-fold differential xsec dsi/dE1dcostheta  versus E1 for given Enu, costheta
  subroutine show_HNV_free_elepton_ct(process_ID,ml_out,Enu,costheta,debugflag)

    use particleProperties, only : baryon, meson
    use idTable, only : nucleon, pion
    use lepton_kinematics_free
    use lepton_kinematics_free_FULL



    implicit none
    integer, intent(in) :: process_ID
    real, intent(in) :: ml_out, Enu, costheta
    logical, optional, intent(in) :: debugflag



    real :: mN, mpi
    real :: E1, numin, numax, E1min, E1max
    integer :: j ! counter
    real :: xsec_cDp, xsec_Np, xsec_cNp, xsec_CT,xsec_pp,xsec_pF
    logical :: debug

    if (present(debugflag)) then
       debug=debugflag
    else
       debug=.true.
    end if

    write(10,*) '# Double differential x-sec: dsi/E1/dcostheta'
    write(10,'(2(A,g10.4))') '# Enu=', Enu, '   costheta=', costheta

    mN=baryon(nucleon)%mass
    mpi=meson(pion)%mass



    call  minmaxNu_Enu_ct(Enu,costheta,(mN+mpi),mN,ml_out,numin,numax) ! 1.08 is Wmin
    E1min=Enu-numax
    E1max=Enu-numin
    if (debug) write(10,'(2(A,g10.4))') '# E1min=', E1min, ' E1max=', E1max

    write(10,*) '# Integrated over costheta_Wpi:'
    write(10,*) '# Enu costheta  E1   xsec_cDp  xsec _Np  xsec_cNp xsec_CT xsec_pp xsec_pF'



    E1loop:        do j=1,200
       E1=0.01*j
       if (E1.le.E1min) cycle
       if (E1.ge.E1max) exit


       call HNV_free_elepton_ct(process_ID,Enu,E1,costheta,mN,ml_out,mN,mpi,xsec_cDp,xsec_Np,xsec_cNp,xsec_CT,xsec_pp,xsec_pF) ! this is for p pi+ final state


       write(10,'(3g11.5,6g14.5)')  Enu, costheta,  E1, xsec_cDp, xsec_Np, xsec_cNp, xsec_CT, xsec_pp, xsec_pF

    end do E1loop

  end subroutine show_HNV_free_elepton_ct





  ! provides output for 3-fold differential xsec dsi/dE1dcosthetadcosthetaPi   versus costhetapi for given Enu, costheta, E1
  subroutine show_HNV_free_elepton_ct_ctPi(process_ID,ml_out,Enu,costheta,E1,debugflag)

    use particleProperties, only : baryon, meson
    use idTable, only : nucleon, pion
    use lepton_kinematics_free
    use lepton_kinematics_free_FULL



    implicit none
    integer, intent(in) :: process_ID
    real, intent(in) :: ml_out, Enu, costheta, E1
    logical, optional, intent(in) :: debugflag



    real :: mN, mpi
    real :: ctPi, numin, numax, E1min, E1max
    integer :: j ! counter
    real :: xsec_cDp, xsec_Np, xsec_cNp, xsec_CT,xsec_pp,xsec_pF
    logical :: debug

    if (present(debugflag)) then
       debug=debugflag
    else
       debug=.true.
    end if

    write(10,*) '# Double differential x-sec: dsi/E1/dcostheta'
    write(10,'(2(A,g10.4))') '# Enu=', Enu, ' costheta=', costheta

    mN=baryon(nucleon)%mass
    mpi=meson(pion)%mass



    call  minmaxNu_Enu_ct(Enu,costheta,(mN+mpi),mN,ml_out,numin,numax) ! 1.08 is Wmin
    E1min=Enu-numax
    E1max=Enu-numin
    !if (E1.ge.E1max .or. E1.le.E1min)  exit

    write(10,*) '# over Phi_pi:  Enu, costheta, E1, costheta_pi xsec_cDp xsec_Np  xsec_cNP  xsec_CT, xsec_pp, xsec_pF'



    ct_pi:        do j=1,200
       ctPi=-1.+0.01*j


       call HNV_free_elepton_ct_ctPi(process_ID,Enu,E1,costheta,ctPi,&
            &mN,ml_out,mN,mpi,xsec_cDp,xsec_Np,xsec_cNp,xsec_CT,xsec_pp,xsec_pF)
       ! this is for proton target, p pi+ final state for CC,  p pi0 + n pi+ for NC
       ! the only exception is directNucleon for CC, which is for the neutron target

       write(10,'(4g11.5,6g14.5)')  Enu, costheta,  E1, ctPi, xsec_cDp, xsec_Np, xsec_cNp, xsec_CT, xsec_pp, xsec_pF

    end do ct_pi

  end subroutine show_HNV_free_elepton_ct_ctPi








end program test_bgr
