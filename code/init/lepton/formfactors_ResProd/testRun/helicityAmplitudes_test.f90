!***************************************************************************************************
!****p* /helicityAmplitudes_test
! NAME
! program helicityAmplitudes_test
! PURPOSE
! * Checks our form factors for consistency with the underlying helicity amplitudes.
! * Therefore the MAID amplitudes are printed to files and ours are calculated using  
!   formulas E.7-E.16 of the Phd-thesis of Oliver Buss.
!***************************************************************************************************
program helicityAmplitudes_test
  ! Tests our form factors

  use particleProperties, only: initParticleProperties
  implicit none

  call initParticleProperties
  
  ! Print MAID's helicity amplitudes
  call printMAID()

  ! Print our helicity amplitudes
  call printOurs()
  call printOurs_atPhotonPoint()


contains


subroutine printOurs()
  ! Prints the helicity amplitudes at the pole position as a function of Q^2 for the following resonances:
  ! * delta,D13_1520,S11_1535,P11_1440
  ! Output to files
  ! * Helicities_proton.dat
  ! * Helicities_neutron.dat
  ! 
  ! Here we use the Formfactors as input and apply formulas E.7-E.16 of the Phd-thesis of Oliver Buss
  !

  use IdTable, only : delta,D13_1520,S11_1535,P11_1440,nucleon
  use leptonicID,only : EM
  !se helicityAmplitudes
  use formFactor_ResProd
  use particleProperties, only: hadron
  use constants, only: pi, alphaQED, mN

  integer, dimension (1:4) :: id=(/delta,D13_1520,S11_1535,P11_1440 /)
  real, dimension(1:8) :: f
  logical :: ff_set
  real, dimension(1:4) :: A12,A32, S12      ! First Index Resonance ID, Second Index : 1=MAID03, 2=MAID05
  integer :: i,targetCharge
  real :: QS, W, factor
  real, parameter :: delta_QS=0.01

  write(*,'(A)') 
  write(*,'(2A)') "Routine prints our helicity amplitudes at the pole position as a function", &
       & " of Q^2 for the following resonances:"
  write(*,'(A)') " * delta,D13_1520,S11_1535,P11_1440"
  write(*,'(A)') 
  write(*,'(A)') "Output to files"
  write(*,'(A)') " * Helicities_proton.dat"
  write(*,'(A)') " * Helicities_neutron.dat"
  write(*,'(A)') 

  nucleonCharge_loop: do targetCharge=0,1
     if(targetCharge.eq.1) then
        open(100,file="Helicities_proton.dat")
        write(100,*) "# Helicity Amplitudes for gamma Proton from MAID 03"
     else
        open(100,file="Helicities_neutron.dat")
        write(100,*) "# Helicity Amplitudes for gamma Neutron from MAID 03"
     end if
     write(100,*) "# QSquared (GeV^2), (A_1_2, A_3_2, S_1_2) for  delta,D13_1520,S11_1535,P11_1440"
     write(100,*) "# e.g.: column 2= A_1_2 (Delta), 3= A_3_2(Delta), ..., 5=A_1_2(D13_1520), ...."

     QS=0.
     QS_loop: do
        resonance_loop : do i=1,4
           W=hadron(id(i))%mass
           f=getFormfactor_Res(Qs,W,id(i),targetCharge,EM,FF_set)
           Select case(id(i))
           case(delta)
              ! S=3/2, P=+1
              factor=sqrt(pi*alphaQED/(3.*mN)*((w-mN)**2+QS)/(W**2-mN**2))
              A12(i)=factor*( &
                   &   f(1)/mN*(mN**2+mN*W+QS)/W           &
                   &  -f(2)/mN**2*(W**2-mN**2-QS)/2.      &
                   &  -f(3)/mN**2*(W**2-mN**2+QS)/2.      &
                   & )
              !
              factor=sqrt(pi*alphaQED/(mN)*((w-mN)**2+QS)/(W**2-mN**2))
              A32(i)=factor*( &
                   &   f(1)/mN*(mN+W)           &
                   &  +f(2)/mN**2*(W**2-mN**2-QS)/2.      &
                   &  +f(3)/mN**2*(W**2-mN**2+QS)/2.      &
                   & )
              !
              factor=sqrt(pi*alphaQED/(6.*mN)*((w-mN)**2+QS)/(W**2-mN**2)) &
                   &   * sqrt(((W-mN)**2+QS)*((W+mN)**2+QS))/W**2
              S12(i)=factor*( &
                   &   f(1)/mN*W           &
                   &  +f(2)/mN**2*W**2      &
                   &  +f(3)/mN**2*(W**2+mN**2+QS)/2.      &
                   & )
              case(D13_1520)
              ! S=3/2, P=+1
              factor=sqrt(pi*alphaQED/(3.*mN)*((w+mN)**2+QS)/(W**2-mN**2))
              A12(i)=factor*( &
                   &   f(1)/mN*(mN**2-mN*W+QS)/W           &
                   &  -f(2)/mN**2*(W**2-mN**2-QS)/2.      &
                   &  -f(3)/mN**2*(W**2-mN**2+QS)/2.      &
                   & )
              !
              factor=sqrt(pi*alphaQED/(mN)*((w+mN)**2+QS)/(W**2-mN**2))
              A32(i)=factor*( &
                   &   f(1)/mN*(mN-W)           &
                   &  -f(2)/mN**2*(W**2-mN**2-QS)/2.      &
                   &  -f(3)/mN**2*(W**2-mN**2+QS)/2.      &
                   & )
              !
              factor=-sqrt(pi*alphaQED/(6.*mN)*((w+mN)**2+QS)/(W**2-mN**2)) &
                   &   * sqrt(((W-mN)**2+QS)*((W+mN)**2+QS))/W**2
              S12(i)=factor*( &
                   &   f(1)/mN*W           &
                   &  +f(2)/mN**2*W**2      &
                   &  +f(3)/mN**2*(W**2+mN**2+QS)/2.      &
                   & )

           case(S11_1535)
              ! S=1/2, P=-1
              factor=sqrt(2.*pi*alphaQED/mN*((w+mN)**2+QS)/(W**2-mN**2))
              A12(i)=factor*(QS/(4.*mN**2)*f(1)+(W-mN)/2./mN*f(2))
              A32(i)=0.
              factor=sqrt(pi*alphaQED/mN*((mN-W)**2+QS)/(W**2-mN**2))*((W+Mn)**2+qs)/(4.*w*Mn)
              S12(i)=factor*((W-mN)/(2.*mN)*f(1)-f(2))
           case(P11_1440)
              ! S=1/2, P=+1
              factor=sqrt(2.*pi*alphaQED/mN*((w-mN)**2+QS)/(W**2-mN**2))
              A12(i)=factor*(QS/(4.*mN**2)*f(1)+(W+mN)/2./mN*f(2))
              A32(i)=0.
              factor=-sqrt(pi*alphaQED/mN*((mN+W)**2+QS)/(W**2-mN**2))*((W-Mn)**2+qs)/(4.*w*Mn)
              S12(i)=factor*((W+mN)/(2.*mN)*f(1)-f(2))
              case Default
              STOP 'Resonance not included! STOP!'
           end Select
        end do resonance_loop
        write(100,'(13E15.4)')  QS,( (/ A12(i),A32(i),S12(i) /), i=1,4 )
        QS=QS+Delta_QS
        if(QS.gt.2) exit QS_Loop
     end do QS_loop
     close(100)
  end do nucleonCharge_loop


end subroutine printOurs



subroutine printOurs_atPhotonPoint()
  ! Prints the helicity amplitudes at the pole position as a function of Q^2 for the following resonances:
  ! * delta,D13_1520,S11_1535,P11_1440
  ! Output to files
  ! * Helicities_proton.dat
  ! * Helicities_neutron.dat
  ! 
  ! Here we use the Formfactors as input and apply formulas E.7-E.16 of the Phd-thesis of Oliver Buss
  !

  use IdTable, only : delta,D13_1520,S11_1535,P11_1440,nucleon
  use leptonicID,only : EM
  !se helicityAmplitudes
  use formFactor_ResProd
  use constants, only: pi, alphaQED, mN

  integer, dimension (1:4) :: id=(/delta,D13_1520,S11_1535,P11_1440 /)
  real, dimension(1:8) :: f
  logical :: ff_set
  real, dimension(1:4) :: A12,A32, S12      ! First Index Resonance ID, Second Index : 1=MAID03, 2=MAID05
  integer :: i,targetCharge
  real :: QS, W, factor
  real, parameter :: delta_W=0.01

  write(*,'(A)') 
  write(*,'(2A)') "Routine prints our helicity amplitudes at the photon point (Q^2=0) as a function of W" & 
       & ,"for the following resonances:"
  write(*,'(A)') " * delta,D13_1520,S11_1535,P11_1440"
  write(*,'(A)') 
  write(*,'(A)') "Output to files"
  write(*,'(A)') " * Helicities_proton_.dat"
  write(*,'(A)') " * Helicities_neutron.dat"
  write(*,'(A)') 

  nucleonCharge_loop: do targetCharge=0,1
     if(targetCharge.eq.1) then
        open(100,file="Helicities_proton_atphotonPoint.dat")
        write(100,*) "# Helicity Amplitudes for gamma Proton from MAID 03"
     else
        open(100,file="Helicities_neutron_atphotonPoint.dat")
        write(100,*) "# Helicity Amplitudes for gamma Neutron from MAID 03"
     end if
     write(100,*) "# W (GeV), (A_1_2, A_3_2, S_1_2) for  delta,D13_1520,S11_1535,P11_1440"
     write(100,*) "# e.g.: column 2= A_1_2 (Delta), 3= A_3_2(Delta), ..., 5=A_1_2(D13_1520), ...."

     W=1.
     QS=0.
     W_loop: do
        resonance_loop : do i=1,4
           f=getFormfactor_Res(Qs,W,id(i),targetCharge,EM,FF_set)
           Select case(id(i))
           case(delta)
              ! S=3/2, P=+1
              factor=sqrt(pi*alphaQED/(3.*mN)*((w-mN)**2+QS)/(W**2-mN**2))
              A12(i)=factor*( &
                   &   f(1)/mN*(mN**2+mN*W+QS)/W           &
                   &  -f(2)/mN**2*(W**2-mN**2-QS)/2.      &
                   &  -f(3)/mN**2*(W**2-mN**2+QS)/2.      &
                   & )
              !
              factor=sqrt(pi*alphaQED/(mN)*((w-mN)**2+QS)/(W**2-mN**2))
              A32(i)=factor*( &
                   &   f(1)/mN*(mN+W)           &
                   &  +f(2)/mN**2*(W**2-mN**2-QS)/2.      &
                   &  +f(3)/mN**2*(W**2-mN**2+QS)/2.      &
                   & )
              !
              factor=sqrt(pi*alphaQED/(6.*mN)*((w-mN)**2+QS)/(W**2-mN**2)) &
                   &   * sqrt(((W-mN)**2+QS)*((W+mN)**2+QS))/W**2
              S12(i)=factor*( &
                   &   f(1)/mN*W           &
                   &  +f(2)/mN**2*W**2      &
                   &  +f(3)/mN**2*(W**2+mN**2+QS)/2.      &
                   & )
              case(D13_1520)
              ! S=3/2, P=+1
              factor=sqrt(pi*alphaQED/(3.*mN)*((w+mN)**2+QS)/(W**2-mN**2))
              A12(i)=factor*( &
                   &   f(1)/mN*(mN**2-mN*W+QS)/W           &
                   &  -f(2)/mN**2*(W**2-mN**2-QS)/2.      &
                   &  -f(3)/mN**2*(W**2-mN**2+QS)/2.      &
                   & )
              !
              factor=sqrt(pi*alphaQED/(mN)*((w+mN)**2+QS)/(W**2-mN**2))
              A32(i)=factor*( &
                   &   f(1)/mN*(mN-W)           &
                   &  -f(2)/mN**2*(W**2-mN**2-QS)/2.      &
                   &  -f(3)/mN**2*(W**2-mN**2+QS)/2.      &
                   & )
              !
              factor=-sqrt(pi*alphaQED/(6.*mN)*((w+mN)**2+QS)/(W**2-mN**2)) &
                   &   * sqrt(((W-mN)**2+QS)*((W+mN)**2+QS))/W**2
              S12(i)=factor*( &
                   &   f(1)/mN*W           &
                   &  +f(2)/mN**2*W**2      &
                   &  +f(3)/mN**2*(W**2+mN**2+QS)/2.      &
                   & )

           case(S11_1535)
              ! S=1/2, P=-1
              factor=sqrt(2.*pi*alphaQED/mN*((w+mN)**2+QS)/(W**2-mN**2))
              A12(i)=factor*(QS/(4.*mN**2)*f(1)+(W-mN)/2./mN*f(2))
              A32(i)=0.
              factor=sqrt(pi*alphaQED/mN*((mN-W)**2+QS)/(W**2-mN**2))*((W+Mn)**2+qs)/(4.*w*Mn)
              S12(i)=factor*((W-mN)/(2.*mN)*f(1)-f(2))
              
           case(P11_1440)
              ! S=1/2, P=+1
              factor=sqrt(2.*pi*alphaQED/mN*((w-mN)**2+QS)/(W**2-mN**2))
              A12(i)=factor*(QS/(4.*mN**2)*f(1)+(W+mN)/2./mN*f(2))
              A32(i)=0.
              factor=-sqrt(pi*alphaQED/mN*((mN+W)**2+QS)/(W**2-mN**2))*((W-Mn)**2+qs)/(4.*w*Mn)
              S12(i)=factor*((W+mN)/(2.*mN)*f(1)-f(2))
              case Default
              STOP 'Resonance not included! STOP!'
           end Select
        end do resonance_loop
        write(100,'(13E15.4)')  W,( (/ A12(i),A32(i),S12(i) /), i=1,4 )
        W=W+Delta_W
        if(W.gt.3) exit W_Loop
     end do W_loop
     close(100)
  end do nucleonCharge_loop


end subroutine printOurs_atPhotonPoint




subroutine printMAID()
  ! Prints MAID helicity amplitudes at the pole position as a function of Q^2 for the following resonances:
  ! * delta,D13_1520,S11_1535,P11_1440
  ! Output to files
  ! * Helicities_MAID03_proton.dat
  ! * Helicities_MAID05_proton.dat
  ! * Helicities_MAID03_neutron.dat
  ! * Helicities_MAID05_neutron.dat
  use IdTable, only : delta,D13_1520,S11_1535,P11_1440,nucleon
  use leptonicID,only : EM
  use helicityAmplitudes

  integer, dimension (1:4) :: id=(/delta,D13_1520,S11_1535,P11_1440 /)
  integer :: i
  real, dimension(1:4,1:2) :: A12maid,A32maid, S12MAID      ! First Index Resonance ID, Second Index : 1=MAID03, 2=MAID05
  real :: QS
  real, parameter :: delta_QS=0.01


  write(*,'(A)') 
  write(*,'(A)') "Routine prints MAID helicity amplitudes at the pole position as a", &
       & "function of Q^2 for the following resonances:"
  write(*,'(A)') " * delta,D13_1520,S11_1535,P11_1440"
  write(*,'(A)') 
  write(*,'(A)') "Output to files"
  write(*,'(A)') " * Helicities_MAID03_proton.dat"
  write(*,'(A)') " * Helicities_MAID05_proton.dat"
  write(*,'(A)') " * Helicities_MAID03_neutron.dat"
  write(*,'(A)') " * Helicities_MAID05_neutron.dat"
  write(*,'(A)') 


  open(101,file="Helicities_MAID03_proton.dat")
  open(102,file="Helicities_MAID05_proton.dat")
  write(101,*) "# Helicity Amplitudes for gamma Proton from MAID 03"
  write(101,*) "# QSquared (GeV^2), (A_1_2, A_3_2, S_1_2) for  delta,D13_1520,S11_1535,P11_1440"
  write(101,*) "# e.g.: column 2= A_1_2 (Delta), 3= A_3_2(Delta), ..., 5=A_1_2(D13_1520), ...."
  write(102,*) "# Helicity Amplitudes for gamma Proton from MAID 05"
  write(102,*) "# QSquared (GeV^2), (A_1_2, A_3_2, S_1_2) for  delta,D13_1520,S11_1535,P11_1440"
  write(102,*) "# e.g.: column 2= A_1_2 (Delta), 3= A_3_2(Delta), ..., 5=A_1_2(D13_1520), ...."

  QS=0.
  do 
     do i=1,4
        call get_helicityAmplitudes(1,id(i),Qs,A12maid(i,1),A32maid(i,1),S12maid(i,1),1)  !MAID 2003
        call get_helicityAmplitudes(1,id(i),Qs,A12maid(i,2),A32maid(i,2),S12maid(i,2),2)  !MAID 2005
     end do
     write(101,'(13E15.4)')  QS,( (/ A12maid(i,1),A32maid(i,1),S12maid(i,1) /), i=1,4 )
     write(102,'(13E15.4)')  QS,( (/ A12maid(i,2),A32maid(i,2),S12maid(i,2) /), i=1,4 )
     QS=QS+delta_QS
     if(QS.gt.2) exit
  End do
  close(101)
  close(102)

  open(101,file="Helicities_MAID03_neutron.dat")
  open(102,file="Helicities_MAID05_neutron.dat")
  write(101,*) "# Helicity Amplitudes for gamma Neutron from MAID 03"
  write(101,*) "# QSquared (GeV^2), (A_1_2, A_3_2, S_1_2) for  delta,D13_1520,S11_1535,P11_1440"
  write(101,*) "# e.g.: column 2= A_1_2 (Delta), 3= A_3_2(Delta), ..., 5=A_1_2(D13_1520), ...."
  write(102,*) "# Helicity Amplitudes for gamma Neutron from MAID 05"
  write(102,*) "# QSquared (GeV^2), (A_1_2, A_3_2, S_1_2) for  delta,D13_1520,S11_1535,P11_1440"
  write(102,*) "# e.g.: column 2= A_1_2 (Delta), 3= A_3_2(Delta), ..., 5=A_1_2(D13_1520), ...."

  QS=0.
  do 
     do i=1,4
        call get_helicityAmplitudes(0,id(i),Qs,A12maid(i,1),A32maid(i,1),S12maid(i,1),1)  !MAID 2003
        call get_helicityAmplitudes(0,id(i),Qs,A12maid(i,2),A32maid(i,2),S12maid(i,2),2)  !MAID 2005
     end do
     write(101,'(13E15.4)')  QS,( (/ A12maid(i,1),A32maid(i,1),S12maid(i,1) /), i=1,4 )
     write(102,'(13E15.4)')  QS,( (/ A12maid(i,2),A32maid(i,2),S12maid(i,2) /), i=1,4 )
     QS=QS+delta_QS
     if(QS.gt.2) exit
  End do
  close(101)
  close(102)

end subroutine printMAID

end program helicityAmplitudes_test
