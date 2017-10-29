program test_F2

  use leptonicID
  use constants
  use output
  use version
  use inputGeneral
  use particleProperties, only: initParticleProperties

  implicit none

  character(35) :: name='dummy'
  logical :: F2overXi, get_xsec, debug_F2
  integer :: process_ID, charge_in, F2from, numsteps
  real :: ml_out, Enu, costheta, Q2begin, Q2step


  ! defaul values

  F2overXi=.true.
  get_xsec=.false.
  process_ID=EM
  charge_in=1
  F2from=1 ! 1=from background xsec, 2=from full xsec, 3=Oliver's resonance
  debug_F2=.true.

  Q2begin=0.225
  Q2step=0.05
  numsteps=38

  Enu=3.
  costheta=0.6
  !

  call PrintVersion
  call readInputGeneral
  call initParticleProperties

  ! program here

  call readinput_test_F2()


  ! ****************************
  ! integral F2  over Nachtman variable xi
  ! ****************************
  if (F2overXi) then
     name='F2overXi.dat'
     Open(10,File=name)
     write(10,*)'# F2 structure function extracted from the xsec and then integrated over xi'

     if(debug_F2) then
        name='F2.dat'
        Open(11,File=name)
        write(11,*) '# F2 structure function extracted from the xsec'
        write(11,*) '#  Q2 costheta E1 xsec x xi_out F2'
     end if

     call show_integral_over_F2(process_ID,ml_out,charge_in)

  end if



  ! ****************************
  ! double differential xsection
  ! ****************************
  if (get_xsec) then

     name='xsec-for-F2.dat'
     Open(10,File=name)
     write(10,*)'# double differentila xsec dsi/ dcostheta dE1 '

     call show_xsec(process_ID,charge_in,ml_out,Enu,costheta)

  end if


contains



  subroutine readinput_test_F2()
    use output, only: Write_ReadingInput

    integer :: ios

    NAMELIST /nl_test_F2/ F2overXi, get_xsec, process_ID, charge_in, F2from
    NAMELIST /nl_xsec/ ml_out,Enu,costheta
    NAMELIST /nl_F2/ debug_F2, Q2begin,Q2step, numsteps

    call Write_ReadingInput('nl_test_F2',0)
    rewind(5)
    read(5,nml=nl_test_F2,IOSTAT=ios)
    call Write_ReadingInput('nl_test_F2',0,ios)
    write(*,*) 'process_ID=', process_ID, '   charge_in=',charge_in, '    F2from=',F2from
    call Write_ReadingInput('nl_test_F2',1)

    if (process_ID.eq.EM) ml_out=melec
    if (process_ID.eq.CC) ml_out=mmuon



    if (F2overXi) then
       write(*,*) ' Calculating integral of F2 over Nachtman variable xi'
       call Write_ReadingInput('nl_F2',0)
       rewind(5)
       read(5,nml=nl_F2,IOSTAT=ios)
       call Write_ReadingInput('nl_F2',0,ios)
       write(*,*) 'debug_F2=', debug_F2
       call Write_ReadingInput('nl_F2',1)
    end if

    if (get_xsec) then
       write(*,*) ' Calculating double differential xsec'
       call Write_ReadingInput('nl_xsec',0)
       rewind(5)
       read(5,nml=nl_xsec,IOSTAT=ios)
       call Write_ReadingInput('nl_xsec',0,ios)
       write(*,*) 'ml_out=', ml_out, '   Enu=',Enu, '   costheta=',costheta
       call Write_ReadingInput('nl_xsec',1)
    end if

  end subroutine readinput_test_F2





  subroutine show_integral_over_F2(process_ID,ml,charge_in)

    use constants, only: mN
    use StructureFunctions

    integer, intent(in) :: process_ID
    real,  intent(in)   :: ml ! mass of the outgoing lepton
    integer, intent(in) :: charge_in

    real :: Q2, ximin, ximax, F2overXi
    integer :: i

    write(10,'(A)') '# Q2 ximin ximax F2overXi'
    if(debug_F2 .and. charge_in.eq.10) write(11,*) '#  Even  indexies are protons, odd indexies are neutrons'


    Q2loop: do i=0,numsteps
       Q2=Q2begin+Q2step*i

       ximin = Nachtmann_Q2_W(Q2,2.0,mN)
       ximax = Nachtmann_Q2_W(Q2,1.1,mN)

       ! write(10,'(A,6g12.5)') '# Q2=',Q2, 'ximin=',ximin, 'ximax=', ximax


       if (charge_in.eq.10) then

          F2overXi=(intF2_to_duality(process_ID,ml,F2from,1,Q2,mN,ximin,ximax,debug_F2) &
               & + intF2_to_duality(process_ID,ml,F2from,0,Q2,mN,ximin,ximax,debug_F2))/2.

       else

          F2overXi=intF2_to_duality(process_ID,ml,F2from,charge_in,Q2,mN,ximin,ximax,debug_F2)

       end if

       write(10,'(6g12.5)') Q2, ximin, ximax, F2overXi

       write(*,*) ''
       write(*,*) ''
       write(*,*) ''

       if (debug_F2) then
          write(11,*) ''
          write(11,*) ''
       end if

    end do Q2loop

  end subroutine show_integral_over_F2







  subroutine show_xsec(process_ID,charge_in,ml_out,Enu,costheta)

    use singlePionProductionMAIDlike
    use StructureFunctions, only : MAIDlike_dE1dcostheta_res
    use lepton_kinematics_free
    use minkowski, only: abs4, abs4Sq
    use constants, only: mN

    integer, intent(in) :: process_ID
    integer, intent(in) :: charge_in
    real, intent(in) :: ml_out, Enu, costheta

    real :: xsec, E1min, E1max, knu, sintheta, E1, k1, W, Q2
    integer :: charge_out, pion_charge_out
    integer :: i
    real, dimension(0:3) :: k_in, k_out, p_in
    real, dimension(1:3) :: position

    call minmaxE1_Enu_ct(Enu,costheta,1.1,mN,ml_out,E1min,E1max)
    !E1min=0.57

    write(10,'(A,I5,4(A,g12.5))') '#  process_ID=', process_ID, '   Enu=  ', Enu, &
         & '   costheta=',costheta, '   E1min=',E1min, '   E1max=',E1max
    write(*,'(A,I5,4(A,g12.5))') '#  process_ID=', process_ID, '   Enu=  ', Enu, &
         & '   costheta=',costheta, '   E1min=',E1min, '   E1max=',E1max
    write(10,*) '#  costheta E1 xsec W2 Q2'

    ! incoming 4-vectors
    if(process_ID.eq.CC) knu=Enu
    if(process_ID.eq.EM) knu=sqrt((Enu-ml_out)*(Enu+ml_out))
    sintheta=sqrt((1.-costheta)*(1.+costheta))
    k_in=(/Enu,0.,0.,knu/)
    p_in = (/mN, 0., 0., 0./)
    position=(/0., 0., 0./)

    do i=1,199

       ! outgoing lepton   4-vector
       E1=E1max-(E1max-E1min)*i/200.
       k1=sqrt((E1-ml_out)*(E1+ml_out))
       k_out=(/E1, -k1*sintheta, 0., k1*costheta/)
       W=abs4(p_in+k_in-k_out)
       Q2=-abs4Sq(k_in-k_out)
       !write(*,*) 'W=',W
       if (W>2.) cycle


       if (F2from.eq.3) then

          xsec=MAIDlike_dE1dcostheta_res(process_ID,charge_in,charge_out,  &
               & k_in, k_out, costheta, p_in, position, .false., ml_out)

       else

          xsec=0
          do charge_out=0,1
             !charge_out=1
             if (process_ID.eq.1) then
                pion_charge_out=charge_in-charge_out
             else
                pion_charge_out=charge_in+1-charge_out

                if (pion_charge_out>1 .or. pion_charge_out<-1) cycle

             end if

             if (F2from.eq.2) then
                xsec=xsec+MAIDlike_dE1dcostheta_full(process_ID,charge_in,charge_out,pion_charge_out, &
                     & k_in,k_out,costheta,p_in,position,.false.,ml_out)

             else
                xsec=xsec+MAIDlike_dE1dcostheta_bgr(process_ID,charge_in,charge_out,pion_charge_out, &
                     & k_in,k_out,costheta,p_in,position,.false.,ml_out)

             end if
          end do

       end if

       xsec=xsec*0.38938*1.e6 ! convertion from GeV^{-2} and then to nbarn
       write(10,'(5(g14.5))') costheta, E1, xsec, W*W, Q2

       write(*,'(5(A,g12.5))') 'E1=',E1, '   costheta=',costheta, '   W2=',W**2, '   xsec=', xsec, 'Q2=', Q2

    end do

  end subroutine show_xsec





!!$! was done to check the second solution for Oliver, do not need it anymore
!!$subroutine show_lowBoundary_for2costhetaPiSolutions()
!!$
!!$use lepton_kinematics_free, only : lowBoundary_for2costhetaPiSolutions
!!$
!!$ implicit none
!!$
!!$ real :: W, W2, mN1,mpi,EWmin,numin,Q2min
!!$ integer :: i
!!$
!!$    name='low-boundary.dat'
!!$     Open(20,File=name)
!!$     write(20,*)'# Low boundaries  EW_min nu_min Q2_min for a given W2,', &
!!$     'where the second solution for the angle between the outgoing  pion and nucleon appears'
!!$     write(20,*)'# W EW_min nu_min Q2_min'
!!$
!!$
!!$
!!$ mN1=0.939
!!$ mpi=0.138
!!$ W=1.08
!!$    do i=1,1000
!!$ W=W+0.001
!!$ W2=W**2
!!$ call lowBoundary_for2costhetaPiSolutions(W2,mN1,mpi,EWmin,numin,Q2min)
!!$ write(20,'(4g12.5)')    W, EWmin, numin, Q2min
!!$ end do
!!$
!!$end subroutine show_lowBoundary_for2costhetaPiSolutions


end program test_F2
