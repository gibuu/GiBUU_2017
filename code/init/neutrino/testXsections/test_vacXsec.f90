
program testXsections
  use inputGeneral
  use neutrino_IDTable
  use leptonicID
  use particleDefinition
  use idTable, only : nucleon,delta,P11_1440,S11_1535,S11_1650,S11_2090, & 
       & D13_1520,D13_1700,D13_2080,D15_1675,G17_2190,P11_1710,P11_2100, &
       & P13_1720,P13,F15_1680,F15_2000,F17_1990,S31_1620,S31_1900,      &
       & D33_1700,D33_1940,D35_1930,D35_2350,P31,P31_1910,P33_1600,      &
       & P33_1920,F35,F35_1905,F37_1950,pion
  use particleProperties
  use constants
  use NeutrinoMatrixElement
  use baryonwidth
  use random
  use FF_QE_nucleonScattering
  use Formfactor_ResProd
  use helicityAmplitudes
  use gauss_integration
  use singlePionProductionMAIDlike
  use neutrinoXsection
  use version
  use output


  implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real, parameter :: scaling=1.        !used to scale the xsection, now 10^-38 cm^2

  character(100) :: name='dummy'

  character(0) :: QEadd=''
  !character(7) :: QEadd='_newMEl'

  character(0) :: RESadd=''
  !character(7) :: RESadd='_nMF100'

  logical :: skipdoubledifferential=.true.
  logical :: skipdifferential=.false.
  logical :: skiptotal=.false.
  logical :: onlyQE=.false.
  logical :: onlymuonCC=.false.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  logical :: FF_set
  real, parameter :: dq=0.005
  real,dimension(0:1) :: A12,A32,S12


  real :: MNucleon, MDelta, MPion

  character(8) :: resname

  integer :: k
  integer, parameter :: max_finalstateID=31

  real:: Qs, F1,F2,FA,FP
  integer :: i

  real :: fpl1,fpl3,fmi1,fmi3,f0pl,f0mi

  integer :: ios

  call PrintVersion

  call readinputGeneral
  call init_database

  NAMELIST /testnuVac/skipdoubledifferential,skipdifferential,skiptotal,onlyQE,onlymuonCC
  call Write_ReadingInput('testnuVac',0)
  rewind(5)
  read(5,nml=testnuVac,IOSTAT=ios)
  call Write_ReadingInput("testnuVac",0,ios)
  call Write_ReadingInput('testnuVac',1)
  write(*,*) 'skipdoubledifferential,skipdifferential,skiptotal,onlyQE,onlymuonCC', &
       & skipdoubledifferential,skipdifferential,skiptotal,onlyQE,onlymuonCC
  

  MNucleon=baryon(Nucleon)%mass
  MDelta=baryon(delta)%mass
  MPion=meson(pion)%mass


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if(.not.skipdifferential) then

     !QE differential cross sections
     write(*,*) '++++QE differential cross sections'

     write(*,*) '++++++++++EM electron' 
     name='QE_dsigdQs_EM_electron'//QEadd//'.dat'
     Open(10,File=name)
     write(10,*)'# QE dsig/dQs EM electron'
     call calcQEdsigdQs(EM,electron)
     close(10)

     write(*,*) '++++++++++charged current muon neutrinos' 
     name='QE_dsigdQs_CC_muon_neutrino'//QEadd//'.dat'
     Open(10,File=name)
     write(10,*)'# QE dsig/dQs CC muon neutrino'
     call calcQEdsigdQs(CC,muon)
     close(10)

     write(*,*) '++++++++++charged current electron neutrinos' 
     name='QE_dsigdQs_CC_electron_neutrino'//QEadd//'.dat'
     Open(10,File=name)
     write(10,*)'# QE dsig/dQs CC electron neutrino'
     call calcQEdsigdQs(CC,electron)
     close(10)

     write(*,*) '++++++++++charged current tau neutrinos' 
     name='QE_dsigdQs_CC_tau_neutrino'//QEadd//'.dat'
     Open(10,File=name)
     write(10,*)'# QE dsig/dQs CC tau neutrino'
     call calcQEdsigdQs(CC,taulepton)
     close(10)

     write(*,*) '++++++++++neutral current neutrinos' 
     name='QE_dsigdQs_NC_neutrino'//QEadd//'.dat'
     Open(10,File=name)
     write(10,*)'# QE dsig/dQs NC neutrino'
     call calcQEdsigdQs(NC,muon)
     close(10)

     write(*,*) '++++++++++charged current muon antineutrinos' 
     name='QE_dsigdQs_CC_muon_antineutrino'//QEadd//'.dat'
     Open(10,File=name)
     write(10,*)'# QE dsig/dQs CC muon antineutrino'
     call calcQEdsigdQs(antiCC,muon)
     close(10)

     write(*,*) '++++++++++charged current electron antineutrinos' 
     name='QE_dsigdQs_CC_electron_antineutrino'//QEadd//'.dat'
     Open(10,File=name)
     write(10,*)'# QE dsig/dQs CC electron antineutrino'
     call calcQEdsigdQs(antiCC,electron)
     close(10)

     write(*,*) '++++++++++charged current tau antineutrinos' 
     name='QE_dsigdQs_CC_tau_antineutrino'//QEadd//'.dat'
     Open(10,File=name)
     write(10,*)'# QE dsig/dQs CC tau antineutrino'
     call calcQEdsigdQs(antiCC,taulepton)
     close(10)

     write(*,*) '++++++++++neutral current antineutrinos' 
     name='QE_dsigdQs_NC_antineutrino'//QEadd//'.dat'
     Open(10,File=name)
     write(10,*)'# QE dsig/dQs NC antineutrino'
     call calcQEdsigdQs(antiNC,muon)
     close(10)


     write(*,*) '++++end of QE differential cross sections'

  end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! QE form factors

  name='QE_CC_FF'//QEadd//'.dat'
  Open(10,File=name)
  write(10,*)'# nucleon FF'
  do i=0,800
     Qs=dq*i
     call formfactors_QE(Qs,CC,0,F1,F2,FA,FP) 
     write(10,'(9F12.5)') Qs, F1,F2,FA,FP
  end do
  close(10)

  name='QE_EM_FF_proton'//QEadd//'.dat'
  Open(10,File=name)
  write(10,*)'# nucleon FF'
  do i=0,800
     Qs=dq*i
     call formfactors_QE(Qs,EM,1,F1,F2,FA,FP) 
     write(10,'(9F12.5)') Qs, F1,F2,FA,FP
  end do
  close(10)

  name='QE_EM_FF_neutron'//QEadd//'.dat'
  Open(10,File=name)
  write(10,*)'# nucleon FF'
  do i=0,800
     Qs=dq*i
     call formfactors_QE(Qs,EM,0,F1,F2,FA,FP) 
     write(10,'(9F12.5)') Qs, F1,F2,FA,FP
  end do
  close(10)

  name='QE_NC_FF_proton'//QEadd//'.dat'
  Open(10,File=name)
  write(10,*)'# nucleon FF'
  do i=0,800
     Qs=dq*i
     call formfactors_QE(Qs,NC,1,F1,F2,FA,FP) 
     write(10,'(9F12.5)') Qs, F1,F2,FA,FP
  end do
  close(10)

  name='QE_NC_FF_neutron'//QEadd//'.dat'
  Open(10,File=name)
  write(10,*)'# nucleon FF'
  do i=0,800
     Qs=dq*i
     call formfactors_QE(Qs,NC,0,F1,F2,FA,FP) 
     write(10,'(9F12.5)') Qs, F1,F2,FA,FP
  end do
  close(10)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !QE total cross sections
  write(*,*) '++++QE total cross sections'

  write(*,*) '++++++++++EM electron' 
  name='QE_sigma_EM_electron'//QEadd//'.dat'
  Open(10,File=name)
  write(10,*)'# QE sigma EM electron'
  call calcQEsigma(EM,electron)
  close(10)

  write(*,*) '++++++++++charged current muon neutrinos' 
  name='QE_sigma_CC_muon_neutrino'//QEadd//'.dat'
  Open(10,File=name)
  write(10,*)'# QE sigma CC muon neutrino'
  call calcQEsigma(CC,muon)
  close(10)

  write(*,*) '++++++++++charged current electron neutrinos' 
  name='QE_sigma_CC_electron_neutrino'//QEadd//'.dat'
  Open(10,File=name)
  write(10,*)'# QE sigma CC electron neutrino'
  call calcQEsigma(CC,electron)
  close(10)

  write(*,*) '++++++++++charged current tau neutrinos' 
  name='QE_sigma_CC_tau_neutrino'//QEadd//'.dat'
  Open(10,File=name)
  write(10,*)'# QE sigma CC tau neutrino'
  call calcQEsigma(CC,taulepton)
  close(10)

  write(*,*) '++++++++++neutral current neutrinos' 
  name='QE_sigma_NC_neutrino'//QEadd//'.dat'
  Open(10,File=name)
  write(10,*)'# QE sigma NC neutrino'
  call calcQEsigma(NC,muon)
  close(10)

  write(*,*) '++++++++++charged current muon antineutrinos' 
  name='QE_sigma_CC_muon_antineutrino'//QEadd//'.dat'
  Open(10,File=name)
  write(10,*)'# QE sigma CC muon antineutrino'
  call calcQEsigma(antiCC,muon)
  close(10)

  write(*,*) '++++++++++charged current electron antineutrinos' 
  name='QE_sigma_CC_electron_antineutrino'//QEadd//'.dat'
  Open(10,File=name)
  write(10,*)'# QE sigma CC electron antineutrino'
  call calcQEsigma(antiCC,electron)
  close(10)

  write(*,*) '++++++++++charged current tau antineutrinos' 
  name='QE_sigma_CC_tau_antineutrino'//QEadd//'.dat'
  Open(10,File=name)
  write(10,*)'# QE sigma CC tau antineutrino'
  call calcQEsigma(antiCC,taulepton)
  close(10)

  write(*,*) '++++++++++neutral current antineutrinos' 
  name='QE_sigma_NC_antineutrino'//QEadd//'.dat'
  Open(10,File=name)
  write(10,*)'# QE sigma NC antineutrino'
  call calcQEsigma(antiNC,muon)
  close(10)


  write(*,*) '++++end of QE total cross sections'

  if(onlyQE) then
     write(*,*) 'END OF QE'
     stop
  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !RES form factors and helicity amplitudes
  do k=2, max_finalstateID

     !exclude those resonances which are not implemented
     if(k.eq.S11_2090) cycle
     if(k.eq.D13_2080) cycle
     if(k.eq.G17_2190) cycle
     if(k.eq.P11_2100) cycle
     if(k.eq.P13) cycle
     if(k.eq.F15_2000) cycle
     if(k.eq.S31_1900) cycle
     if(k.eq.D33_1940) cycle
     if(k.eq.D35_1930) cycle
     if(k.eq.D35_2350) cycle
     if(k.eq.P31) cycle
     if(k.eq.F35) cycle

     resname=baryon(k)%name
     if(k.eq.2) resname='P33_1232'

     name=resname//'_CC_proton_FFres'//RESadd//'.dat'
     Open(10,File=name)
     write(10,*)'# ',resname,' FF'

     name=resname//'_CC_neutron_FFres'//RESadd//'.dat'
     Open(15,File=name)
     write(15,*)'# ',resname,' FF'

     name=resname//'_EM_proton_FFres'//RESadd//'.dat'
     Open(12,File=name)
     write(12,*)'# ',resname,' FF'

     name=resname//'_EM_neutron_FFres'//RESadd//'.dat'
     Open(13,File=name)
     write(13,*)'# ',resname,' FF'

     name=resname//'_NC_proton_FFres'//RESadd//'.dat'
     Open(14,File=name)
     write(14,*)'# ',resname,' FF'

     name=resname//'_NC_neutron_FFres'//RESadd//'.dat'
     Open(16,File=name)
     write(16,*)'# ',resname,' FF'


     name=resname//'_helicity'//RESadd//'.dat'
     Open(11,File=name)
     write(11,*)'# ',resname,' heli'

     do i=0,800
        Qs=dq*i
        write(10,'(9F12.5)') Qs, getFormfactor_Res(Qs,baryon(k)%mass,k,1,CC,FF_set)
        write(15,'(9F12.5)') Qs, getFormfactor_Res(Qs,baryon(k)%mass,k,0,CC,FF_set)
        write(12,'(9F12.5)') Qs, getFormfactor_Res(Qs,baryon(k)%mass,k,1,EM,FF_set)
        write(13,'(9F12.5)') Qs, getFormfactor_Res(Qs,baryon(k)%mass,k,0,EM,FF_set)
        write(14,'(9F12.5)') Qs, getFormfactor_Res(Qs,baryon(k)%mass,k,1,NC,FF_set)
        write(16,'(9F12.5)') Qs, getFormfactor_Res(Qs,baryon(k)%mass,k,0,NC,FF_set)
        call get_helicityAmplitudes(proton,k,Qs,A12(proton),A32(proton),S12(proton))
        call get_helicityAmplitudes(neutron,k,Qs,A12(neutron),A32(neutron),S12(neutron))
        write(11,'(7F12.5)') Qs,A12(proton),A32(proton),S12(proton),A12(neutron),A32(neutron),S12(neutron)
     end do
     close(10)
     close(11)
     close(15)
     close(12)
     close(13)
     close(14)
     close(16)


     !Rein-Sehgal helis for four lowest lying res
     if(k.eq.delta.or.k.eq.P11_1440.or.k.eq.S11_1535.or.k.eq.D13_1520) then
        
        resname=baryon(k)%name
        if(k.eq.2) resname='P33_1232'

        name=resname//'_helicity_reinsehgal_proton'//RESadd//'.dat'
        Open(11,File=name)
        write(11,*)'# ',resname,' heli proton'
        
        name=resname//'_helicity_reinsehgal_neutron'//RESadd//'.dat'
        Open(12,File=name)
        write(12,*)'# ',resname,' heli neutron'

        do i=0,800
           Qs=dq*i
           call reinsehgalhelis(Qs,proton,k,fpl1,fpl3,fmi1,fmi3,f0pl,f0mi)
           write(11,'(7F12.5)') Qs,fpl1,fpl3,fmi1,fmi3,f0pl,f0mi
           call reinsehgalhelis(Qs,neutron,k,fpl1,fpl3,fmi1,fmi3,f0pl,f0mi)
           write(12,'(7F12.5)') Qs,fpl1,fpl3,fmi1,fmi3,f0pl,f0mi
        end do
        close(11)
        close(12)
        
     end if







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     if(.not.skipdoubledifferential) then

        !RES double differential cross sections
        write(*,*) '++++', resname,' double differential cross sections'

        write(*,*) '++++++++++charged current muon neutrinos' 
        name=resname//'_dsigdcosthetadEl_CC_muon_neutrino'//RESadd//'.dat'
        Open(10,File=name)
        write(10,*)'# ',resname,'  dsig/dcosthetadEl CC muon neutrino'
        call calcRESdsigdcosthetadEl(k,CC,muon)
        close(10)

        if(.not.onlymuonCC) then
           write(*,*) '++++++++++EM electron' 
           name=resname//'_dsigdcosthetadEl_EM_electron'//RESadd//'.dat'
           Open(10,File=name)
           write(10,*)'# ',resname,'  dsig/dcosthetadEl EM electron'
           call calcRESdsigdcosthetadEl(k,EM,electron)
           close(10)

           write(*,*) '++++++++++charged current electron neutrinos' 
           name=resname//'_dsigdcosthetadEl_CC_electron_neutrino'//RESadd//'.dat'
           Open(10,File=name)
           write(10,*)'# ',resname,'  dsig/dcosthetadEl CC electron neutrino'
           call calcRESdsigdcosthetadEl(k,CC,electron)
           close(10)

           write(*,*) '++++++++++charged current tau neutrinos' 
           name=resname//'_dsigdcosthetadEl_CC_tau_neutrino'//RESadd//'.dat'
           Open(10,File=name)
           write(10,*)'# ',resname,'  dsig/dcosthetadEl CC tau neutrino'
           call calcRESdsigdcosthetadEl(k,CC,taulepton)
           close(10)

           write(*,*) '++++++++++neutral current neutrinos' 
           name=resname//'_dsigdcosthetadEl_NC_neutrino'//RESadd//'.dat'
           Open(10,File=name)
           write(10,*)'# ',resname,'  dsig/dcosthetadEl NC neutrino'
           call calcRESdsigdcosthetadEl(k,NC,muon)
           close(10)

!!$        write(*,*) '++++++++++charged current muon antineutrinos' 
!!$        name=resname//'_dsigdcosthetadEl_CC_muon_antineutrino'//RESadd//'.dat'
!!$        Open(10,File=name)
!!$        write(10,*)'# ',resname,'  dsig/dcosthetadEl CC muon antineutrino'
!!$        call calcRESdsigdcosthetadEl(k,antiCC,muon)
!!$        close(10)
!!$
!!$        write(*,*) '++++++++++charged current electron antineutrinos'
!!$        name=resname//'_dsigdcosthetadEl_CC_electron_antineutrino'//RESadd//'.dat'
!!$        Open(10,File=name)
!!$        write(10,*)'# ',resname,'  dsig/dcosthetadEl CC electron antineutrino'
!!$        call calcRESdsigdcosthetadEl(k,antiCC,electron)
!!$        close(10)
!!$
!!$        write(*,*) '++++++++++charged current tau antineutrinos' 
!!$        name=resname//'_dsigdcosthetadEl_CC_tau_antineutrino'//RESadd//'.dat'
!!$        Open(10,File=name)
!!$        write(10,*)'# ',resname,'  dsig/dcosthetadEl CC tau antineutrino'
!!$        call calcRESdsigdcosthetadEl(k,antiCC,taulepton)
!!$        close(10)
!!$
!!$        write(*,*) '++++++++++neutral current antineutrinos' 
!!$        name=resname//'_dsigdcosthetadEl_NC_antineutrino'//RESadd//'.dat'
!!$        Open(10,File=name)
!!$        write(10,*)'# ',resname,'  dsig/dcosthetadEl NC antineutrino'
!!$        call calcRESdsigdcosthetadEl(k,antiNC,muon)
!!$        close(10)

        end if
        write(*,*) '++++end of ', resname,' double differential cross sections'




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !RES double differential cross sections
        write(*,*) '++++', resname,' double differential cross sections'


        write(*,*) '++++++++++charged current muon neutrinos' 
        name=resname//'_dsigdQsdW_CC_muon_neutrino'//RESadd//'.dat'
        Open(10,File=name)
        write(10,*)'# ',resname,'  dsig/dQsdW CC muon neutrino'
        call calcRESdsigdQsdW(k,CC,muon)
        close(10)

        if(.not.onlymuonCC) then

           write(*,*) '++++++++++EM electron' 
           name=resname//'_dsigdQsdW_EM_electron'//RESadd//'.dat'
           Open(10,File=name)
           write(10,*)'# ',resname,'  dsig/dQsdW EM electron'
           call calcRESdsigdQsdW(k,EM,electron)
           close(10)


           write(*,*) '++++++++++charged current electron neutrinos' 
           name=resname//'_dsigdQsdW_CC_electron_neutrino'//RESadd//'.dat'
           Open(10,File=name)
           write(10,*)'# ',resname,'  dsig/dQsdW CC electron neutrino'
           call calcRESdsigdQsdW(k,CC,electron)
           close(10)

           write(*,*) '++++++++++charged current tau neutrinos' 
           name=resname//'_dsigdQsdW_CC_tau_neutrino'//RESadd//'.dat'
           Open(10,File=name)
           write(10,*)'# ',resname,'  dsig/dQsdW CC tau neutrino'
           call calcRESdsigdQsdW(k,CC,taulepton)
           close(10)

           write(*,*) '++++++++++neutral current neutrinos' 
           name=resname//'_dsigdQsdW_NC_neutrino'//RESadd//'.dat'
           Open(10,File=name)
           write(10,*)'# ',resname,'  dsig/dQsdW NC neutrino'
           call calcRESdsigdQsdW(k,NC,muon)
           close(10)

!!$        write(*,*) '++++++++++charged current muon antineutrinos' 
!!$        name=resname//'_dsigdQsdW_CC_muon_antineutrino'//RESadd//'.dat'
!!$        Open(10,File=name)
!!$        write(10,*)'# ',resname,'  dsig/dQsdW CC muon antineutrino'
!!$        call calcRESdsigdQsdW(k,antiCC,muon)
!!$        close(10)
!!$
!!$        write(*,*) '++++++++++charged current electron antineutrinos'
!!$        name=resname//'_dsigdQsdW_CC_electron_antineutrino'//RESadd//'.dat'
!!$        Open(10,File=name)
!!$        write(10,*)'# ',resname,'  dsig/dQsdW CC electron antineutrino'
!!$        call calcRESdsigdQsdW(k,antiCC,electron)
!!$        close(10)
!!$
!!$        write(*,*) '++++++++++charged current tau antineutrinos' 
!!$        name=resname//'_dsigdQsdW_CC_tau_antineutrino'//RESadd//'.dat'
!!$        Open(10,File=name)
!!$        write(10,*)'# ',resname,'  dsig/dQsdW CC tau antineutrino'
!!$        call calcRESdsigdQsdW(k,antiCC,taulepton)
!!$        close(10)
!!$
!!$        write(*,*) '++++++++++neutral current antineutrinos' 
!!$        name=resname//'_dsigdQsdW_NC_antineutrino'//RESadd//'.dat'
!!$        Open(10,File=name)
!!$        write(10,*)'# ',resname,'  dsig/dQsdW NC antineutrino'
!!$        call calcRESdsigdQsdW(k,antiNC,muon)
!!$        close(10)

        end if

        write(*,*) '++++end of ', resname,' double differential cross sections'

     end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     if(.not.skipdifferential) then

        !RES differential cross sections
        write(*,*) '++++', resname,' differential cross sections'



        write(*,*) '++++++++++charged current muon neutrinos' 
        name=resname//'_dsigdQs_CC_muon_neutrino'//RESadd//'.dat'
        Open(10,File=name)
        write(10,*)'# ',resname,'  dsig/dQs CC muon neutrino'
        call calcRESdsigdQs(k,CC,muon)
        close(10)

        if(.not.onlymuonCC) then

           write(*,*) '++++++++++EM electron' 
           name=resname//'_dsigdQs_EM_electron'//RESadd//'.dat'
           Open(10,File=name)
           write(10,*)'# ',resname,'  dsig/dQs EM electron'
           call calcRESdsigdQs(k,EM,electron)
           close(10)


           write(*,*) '++++++++++charged current electron neutrinos' 
           name=resname//'_dsigdQs_CC_electron_neutrino'//RESadd//'.dat'
           Open(10,File=name)
           write(10,*)'# ',resname,'  dsig/dQs CC electron neutrino'
           call calcRESdsigdQs(k,CC,electron)
           close(10)

           write(*,*) '++++++++++charged current tau neutrinos' 
           name=resname//'_dsigdQs_CC_tau_neutrino'//RESadd//'.dat'
           Open(10,File=name)
           write(10,*)'# ',resname,'  dsig/dQs CC tau neutrino'
           call calcRESdsigdQs(k,CC,taulepton)
           close(10)

           write(*,*) '++++++++++neutral current neutrinos' 
           name=resname//'_dsigdQs_NC_neutrino'//RESadd//'.dat'
           Open(10,File=name)
           write(10,*)'# ',resname,'  dsig/dQs NC neutrino'
           call calcRESdsigdQs(k,NC,muon)
           close(10)

!!$        write(*,*) '++++++++++charged current muon antineutrinos' 
!!$        name=resname//'_dsigdQs_CC_muon_antineutrino'//RESadd//'.dat'
!!$        Open(10,File=name)
!!$        write(10,*)'# ',resname,'  dsig/dQs CC muon antineutrino'
!!$        call calcRESdsigdQs(k,antiCC,muon)
!!$        close(10)
!!$
!!$        write(*,*) '++++++++++charged current electron antineutrinos' 
!!$        name=resname//'_dsigdQs_CC_electron_antineutrino'//RESadd//'.dat'
!!$        Open(10,File=name)
!!$        write(10,*)'# ',resname,'  dsig/dQs CC electron antineutrino'
!!$        call calcRESdsigdQs(k,antiCC,electron)
!!$        close(10)
!!$
!!$        write(*,*) '++++++++++charged current tau antineutrinos' 
!!$        name=resname//'_dsigdQs_CC_tau_antineutrino'//RESadd//'.dat'
!!$        Open(10,File=name)
!!$        write(10,*)'# ',resname,'  dsig/dQs CC tau antineutrino'
!!$        call calcRESdsigdQs(k,antiCC,taulepton)
!!$        close(10)
!!$
!!$        write(*,*) '++++++++++neutral current antineutrinos' 
!!$        name=resname//'_dsigdQs_NC_antineutrino'//RESadd//'.dat'
!!$        Open(10,File=name)
!!$        write(10,*)'# ',resname,'  dsig/dQs NC antineutrino'
!!$        call calcRESdsigdQs(k,antiNC,muon)
!!$        close(10)
        end if
        write(*,*) '++++end of ', resname,' differential cross sections'

     end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     if(.not.skiptotal) then

        !RES total cross sections
        write(*,*) '++++', resname,' total cross sections'
        
        write(*,*) '++++++++++charged current muon neutrinos' 
        name=resname//'_sigma_CC_muon_neutrino'//RESadd//'.dat'
        Open(10,File=name)
        write(10,*)'# ',resname,'  sigma CC muon neutrino'
        call calcRESsigma(k,CC,muon)
        close(10)

        if(.not.onlymuonCC) then
           
           write(*,*) '++++++++++EM electron' 
           name=resname//'_sigma_EM_electron'//RESadd//'.dat'
           Open(10,File=name)
           write(10,*)'# ',resname,'  sigma EM electron'
           call calcRESsigma(k,EM,electron)
           close(10)

           write(*,*) '++++++++++charged current electron neutrinos' 
           name=resname//'_sigma_CC_electron_neutrino'//RESadd//'.dat'
           Open(10,File=name)
           write(10,*)'# ',resname,'  sigma CC electron neutrino'
           call calcRESsigma(k,CC,electron)
           close(10)

           write(*,*) '++++++++++charged current tau neutrinos' 
           name=resname//'_sigma_CC_tau_neutrino'//RESadd//'.dat'
           Open(10,File=name)
           write(10,*)'# ',resname,'  sigma CC tau neutrino'
           call calcRESsigma(k,CC,taulepton)
           close(10)

           write(*,*) '++++++++++neutral current neutrinos' 
           name=resname//'_sigma_NC_neutrino'//RESadd//'.dat'
           Open(10,File=name)
           write(10,*)'# ',resname,'  sigma NC neutrino'
           call calcRESsigma(k,NC,muon)
           close(10)

!!$        write(*,*) '++++++++++charged current muon antineutrinos' 
!!$        name=resname//'_sigma_CC_muon_antineutrino'//RESadd//'.dat'
!!$        Open(10,File=name)
!!$        write(10,*)'# ',resname,'  sigma CC muon antineutrino'
!!$        call calcRESsigma(k,antiCC,muon)
!!$        close(10)
!!$
!!$        write(*,*) '++++++++++charged current electron antineutrinos' 
!!$        name=resname//'_sigma_CC_electron_antineutrino'//RESadd//'.dat'
!!$        Open(10,File=name)
!!$        write(10,*)'# ',resname,'  sigma CC electron antineutrino'
!!$        call calcRESsigma(k,antiCC,electron)
!!$        close(10)
!!$
!!$        write(*,*) '++++++++++charged current tau antineutrinos' 
!!$        name=resname//'_sigma_CC_tau_antineutrino'//RESadd//'.dat'
!!$        Open(10,File=name)
!!$        write(10,*)'# ',resname,'  sigma CC tau antineutrino'
!!$        call calcRESsigma(k,antiCC,taulepton)
!!$        close(10)
!!$
!!$        write(*,*) '++++++++++neutral current antineutrinos' 
!!$        name=resname//'_sigma_NC_antineutrino'//RESadd//'.dat'
!!$        Open(10,File=name)
!!$        write(10,*)'# ',resname,'  sigma NC antineutrino'
!!$        call calcRESsigma(k,antiNC,muon)
!!$        close(10)

        end if
        write(*,*) '++++end of ', resname,' sigma cross sections'

     end if

  end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! single pion

  write(*,*) '++++charged current muon neutrinos single-pion production' 
  name='single_pion_sigma_CC_muon_neutrino'//RESadd//'.dat'
  Open(10,File=name)
  write(10,*)'# single-pion  sigma CC muon neutrino'
  call  calcCCpionprodsigma
  close(10)
  
  if(.not.onlymuonCC) then
     write(*,*) '++++neutral current neutrinos single-pion production' 
     name='single_pion_sigma_NC_neutrino'//RESadd//'.dat'
     Open(10,File=name)
     write(10,*)'# single-pion  sigma NC neutrino'
     call  calcNCpionprodsigma
     close(10)
  end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  write(*,*) '++++++++++++++++++++++++++++++end of code'

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine calcQEsigma(process_ID,flavor_ID)
    implicit none
    integer, intent(in) :: process_ID, flavor_ID
    real :: enu, QE_thres
    real :: Qs, Qsmin, Qsmax,s,ml_in,ml_out
    real :: dsig_neut, dsig_prot
    integer :: n,n1,j
    real, dimension(3000) :: x,y1,y2
    real :: matrixelementproton, matrixelementneutron
    real, dimension(0:3) :: p_in,p_out,k_in,k_out


    n=5 !integral precision

    call setMasses(process_ID,flavor_ID,ml_in,ml_out)

    write(10,*)'# enu   |   sigma on proton   |   sigma on neutron'
    !set starting energy
    enu=0.
    QE_thres=energythreshold(ml_out,MNucleon)

    energyloop: do 

       enu=enu+0.01
       if(enu.le.QE_thres) cycle
       if(enu.gt.100.) exit  !go until 10 GeV neutrino energy

       s=MNucleon**2+2.*enu*MNucleon 

       call  maxminQs(MNucleon**2+2.*enu*MNucleon, MNucleon,ml_out, Qsmax, Qsmin)

       if(Qsmax.lt.Qsmin) then
          cycle
       end if

       call sg20r(Qsmin,Qsmax,n,x,n1)
       do j=1,n1
          Qs=x(j)

          k_in=(/ enu,0.,0.,enu /)
          p_in=(/ MNucleon,0.,0.,0. /)
          p_out(0)=(2.*MNucleon**2+Qs)/(2.*MNucleon)
          k_out(0)=p_in(0)+k_in(0)-p_out(0)
          k_out(3)=(-Qs-ml_out**2+2.*enu*k_out(0))/(2.*enu)
          k_out(1)=sqrt(k_out(0)**2-k_out(3)**2-ml_out**2)
          k_out(2)=0.
          p_out(1)=-k_out(1)
          p_out(2)=0.
          p_out(3)=enu-k_out(3)

          matrixelementneutron=nuMaEl(process_ID,nucleon,neutron, k_in, k_out, p_in, p_out,MNucleon)
          matrixelementproton=nuMaEl(process_ID,nucleon,proton, k_in, k_out, p_in, p_out,MNucleon)


          y1(j)=scaling*1.9732696817**2*10.**10/16./pi /(s-MNucleon**2)**2*matrixelementneutron

          y2(j)=scaling*1.9732696817**2*10.**10/16./pi /(s-MNucleon**2)**2*matrixelementproton

       end do

       call rg20r(Qsmin,Qsmax,n,y1,dsig_neut)
       call rg20r(Qsmin,Qsmax,n,y2,dsig_prot)

       write(10,'(10g12.5)') enu, dsig_prot, dsig_neut

    end do energyloop

  end subroutine calcQEsigma


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine calcQEdsigdQs(process_ID,flavor_ID)
    implicit none

    integer, intent(in) :: process_ID,flavor_ID
    real :: enu, Qs, Qsmin, Qsmax, QE_thres
    real :: dsig_neut, dsig_prot
    real :: ml_in,ml_out,s
    real :: matrixelementproton, matrixelementneutron
    real, dimension(0:3) :: p_in,p_out,k_in,k_out


    call setMasses(process_ID,flavor_ID,ml_in,ml_out)
    write(10,*)'# enu   |   Qs   |   dsig on proton   |   dsig on neutron'
    !set starting energy
    enu=0.
    QE_thres=energythreshold(ml_out,MNucleon)

    energyloop: do 

       !write(*,*) 'energy loop', enu, QE_thres

       enu=enu+0.25

       if(enu.le.QE_thres) cycle
       if(enu.gt.5.) exit  !go until 5 GeV neutrino energy

       s=MNucleon**2+2.*enu*MNucleon

       !set starting Qs
       Qs=0.

       qsloop: do

          Qs=Qs+0.005 

          call  maxminQs(s, MNucleon, ml_out,Qsmax, Qsmin) 

          !write(*,*) 'Qs loop', Qs, Qsmin, Qsmax

          if(Qs.le.Qsmin) cycle
          if(Qs.ge.Qsmax) exit


          k_in=(/ enu,0.,0.,enu /)
          p_in=(/ MNucleon,0.,0.,0. /)
          p_out(0)=(2.*MNucleon**2+Qs)/(2.*MNucleon)
          k_out(0)=p_in(0)+k_in(0)-p_out(0)
          k_out(3)=(-Qs-ml_out**2+2.*enu*k_out(0))/(2.*enu)
          k_out(1)=sqrt(k_out(0)**2-k_out(3)**2-ml_out**2)
          k_out(2)=0.
          p_out(1)=-k_out(1)
          p_out(2)=0.
          p_out(3)=enu-k_out(3)

          matrixelementneutron=nuMaEl(process_ID,nucleon,neutron, k_in, k_out, p_in, p_out,MNucleon)
          matrixelementproton=nuMaEl(process_ID,nucleon,proton, k_in, k_out, p_in, p_out,MNucleon)


          dsig_neut=scaling*1.9732696817**2*10.**10/16./pi /(s-MNucleon**2)**2*matrixelementneutron

          dsig_prot=scaling*1.9732696817**2*10.**10/16./pi /(s-MNucleon**2)**2*matrixelementproton

          write(10,'(10g12.5)') enu, Qs, dsig_prot, dsig_neut

       end do qsloop

    end do energyloop

  end subroutine calcQEdsigdQs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine calcRESsigma(resonance_ID,process_ID,flavor_ID)
    implicit none

    integer, intent(in) :: resonance_ID
    integer, intent(in) :: process_ID,flavor_ID
    integer :: n,n1,l,n2,j
    real :: enu, QE_thres
    real :: Wmin, Wmax
    real :: dsig_neut, dsig_prot
    real, dimension(3000) :: x,xx,y1,y2,yy1,yy2
    real :: ml_in,ml_out,Mres
    real :: W,width,breitwig,s,Qs,Qsmin,Qsmax
    real :: matrixelementproton, matrixelementneutron
    real, dimension(0:3) :: p_in,p_out,k_in,k_out

    n=3 !integral precision

    MRes=baryon(resonance_ID)%mass

    call setMasses(process_ID,flavor_ID,ml_in,ml_out)

    write(10,*)'# enu   |   sigma on proton   |   sigma on neutron'
    !set starting energy
    enu=0.
    QE_thres=energythreshold(ml_out,MNucleon+Mpion)

    energyloop: do 

       enu=enu+0.05

       if(enu.le.QE_thres) cycle
       if(enu.gt.3.5) exit  !go until 3 GeV neutrino energy

       s=MNucleon**2+2.*enu*MNucleon

       Wmin=MNucleon+Mpion
       Wmax=sqrt(s)-ml_out

       call sg20r(Wmin,Wmax,n,x,n1)
       do j=1,n1
          W=x(j)

          call maxminQs(s,W, ml_out,Qsmax, Qsmin) 

          if(Qsmax.lt.Qsmin) then
             cycle
          end if

          call sg20r(Qsmin,Qsmax,n,xx,n2)
          do l=1,n2
             Qs=xx(l)

             width=FullWidthBaryon(resonance_ID,W)
             breitwig=width/pi*W/((W**2-MRes**2)**2+W**2*width**2)  

             k_in=(/ enu,0.,0.,enu /)
             p_in=(/ MNucleon,0.,0.,0. /)
             p_out(0)=(W**2+MNucleon**2+Qs)/(2.*MNucleon)
             k_out(0)=p_in(0)+k_in(0)-p_out(0)
             k_out(3)=(-Qs-ml_out**2+2.*enu*k_out(0))/(2.*enu)
             k_out(1)=sqrt(k_out(0)**2-k_out(3)**2-ml_out**2)
             k_out(2)=0.
             p_out(1)=-k_out(1)
             p_out(2)=0.
             p_out(3)=enu-k_out(3)

             matrixelementneutron=nuMaEl(process_ID,resonance_ID,neutron, k_in, k_out, p_in, p_out,W)
             matrixelementproton=nuMaEl(process_ID,resonance_ID,proton, k_in, k_out, p_in, p_out,W)


             yy1(l)=scaling*1.9732696817**2*10.**10*W/8./pi/(s-MNucleon**2)**2*breitwig*matrixelementneutron
             yy2(l)=scaling*1.9732696817**2*10.**10*W/8./pi/(s-MNucleon**2)**2*breitwig*matrixelementproton


          end do
          call rg20r(Qsmin,Qsmax,n,yy1,dsig_neut)
          call rg20r(Qsmin,Qsmax,n,yy2,dsig_prot)

          y1(j)=dsig_neut
          y2(j)=dsig_prot
       end do

       call rg20r(Wmin,Wmax,n,y1,dsig_neut)
       call rg20r(Wmin,Wmax,n,y2,dsig_prot)

       write(10,'(10g12.5)') enu, dsig_prot, dsig_neut


    end do energyloop

  end subroutine calcRESsigma


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine calcRESdsigdQs(resonance_ID,process_ID,flavor_ID)
    implicit none

    integer, intent(in) :: resonance_ID,process_ID,flavor_ID
    integer :: n,n1,l
    real :: enu, Qs, s, Qsmin, Qsmax, QE_thres
    real :: Wmin, Wmax
    real :: dsig_neut, dsig_prot
    real, dimension(3000) :: x,y1,y2
    real :: width,breitwig,W
    real :: ml_in,ml_out,Mres
    real :: matrixelementproton, matrixelementneutron
    real, dimension(0:3) :: p_in,p_out,k_in,k_out


    n=5 !integral precision

    MRes=baryon(resonance_ID)%mass
    call setMasses(process_ID,flavor_ID,ml_in,ml_out)

    write(10,*)'# enu   |   Qs   |   dsig on proton   |   dsig on neutron'
    !set starting energy
    enu=0.
    QE_thres=energythreshold(ml_out,MNucleon+Mpion)

    energyloop: do 

       !write(*,*) 'energy loop', enu, QE_thres

       enu=enu+0.5

       if(enu.le.QE_thres) cycle
       if(enu.gt.2.) exit  !go until 2 GeV neutrino energy

       s=MNucleon**2+2.*enu*MNucleon

       !set starting Qs
       Qs=0.

       qsloop: do

          Qs=Qs+0.005

          call maxminQs(s,MNucleon+Mpion, ml_out,Qsmax, Qsmin) 
          !write(*,*) 'Qs loop', Qs, Qsmin, Qsmax
          if(Qs.le.Qsmin) cycle
          if(Qs.ge.Qsmax) exit

          !write(*,*) 'Qs loop', Qs, Qsmin, Qsmax

          Wmin=MNucleon+Mpion
          Wmax=sqrt(MNucleon**2-Qs+2.*MNucleon*(-Qs*enu/(-Qs-ml_out**2)-(ml_out**2+Qs)/(4.*enu)))

          if(Wmax.lt.Wmin) then
             write(*,*)Wmax,Wmin
             cycle
          end if


          call sg20r(Wmin,Wmax,n,x,n1)
          do l=1,n1
             W=x(l)


             width=FullWidthBaryon(resonance_ID,W)
             breitwig=width/pi*W/((W**2-MRes**2)**2+W**2*width**2)  

             k_in=(/ enu,0.,0.,enu /)
             p_in=(/ MNucleon,0.,0.,0. /)
             p_out(0)=(W**2+MNucleon**2+Qs)/(2.*MNucleon)
             k_out(0)=p_in(0)+k_in(0)-p_out(0)
             k_out(3)=(-Qs-ml_out**2+2.*enu*k_out(0))/(2.*enu)
             k_out(1)=sqrt(k_out(0)**2-k_out(3)**2-ml_out**2)
             k_out(2)=0.
             p_out(1)=-k_out(1)
             p_out(2)=0.
             p_out(3)=enu-k_out(3)

             matrixelementneutron=nuMaEl(process_ID,resonance_ID,neutron, k_in, k_out, p_in, p_out,W)
             matrixelementproton=nuMaEl(process_ID,resonance_ID,proton, k_in, k_out, p_in, p_out,W)

             y1(l)=scaling*1.9732696817**2*10.**10*W/8./pi/(s-MNucleon**2)**2*breitwig*matrixelementneutron
             y2(l)=scaling*1.9732696817**2*10.**10*W/8./pi/(s-MNucleon**2)**2*breitwig*matrixelementproton


          end do
          call rg20r(Wmin,Wmax,n,y2,dsig_prot)
          call rg20r(Wmin,Wmax,n,y1,dsig_neut)

          write(10,'(10g12.5)') enu, Qs, dsig_prot, dsig_neut

       end do qsloop


    end do energyloop

  end subroutine calcRESdsigdQs


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine calcRESdsigdQsdW(resonance_ID,process_ID,flavor_ID)
    implicit none


    integer, intent(in) :: resonance_ID,process_ID,flavor_ID
    real :: W, Wmin, Wmax, wurzel, Mres
    real :: breitwig, width
    real :: enu, Qs, s, Qsmin, Qsmax, QE_thres
    real :: dsig_neut, dsig_prot
    real :: ml_in,ml_out
    real :: matrixelementproton, matrixelementneutron
    real, dimension(0:3) :: p_in,p_out,k_in,k_out



    MRes=baryon(resonance_ID)%mass

    call setMasses(process_ID,flavor_ID,ml_in,ml_out)

    write(10,*)'# enu   |   W   |   Qs   |   dsig on proton   |   dsig on neutron'
    !set starting energy
    enu=0.
    QE_thres=energythreshold(ml_out,MNucleon+Mpion)

    energyloop: do 

       !write(*,*) 'energy loop', enu, QE_thres

       enu=enu+0.5

       if(enu.le.QE_thres) cycle
       if(enu.gt.2.) exit  !go until 2 GeV neutrino energy

       s=MNucleon**2+2.*enu*MNucleon


       Wmin=MNucleon+Mpion
       Wmax=sqrt(s)-ml_out

       !set starting W
       W=1.0

       Wloop: do

          W=W+0.005

          !write(*,*)'Wloop', W, Wmin, Wmax

          if(W.le.Wmin) cycle
          if(W.gt.Wmax) exit  

          !set starting Qs
          Qs=0.

          qsloop: do

             Qs=Qs+0.005

             call maxminQs(s, W, ml_out,Qsmax, Qsmin) 

             !write(*,*) 'Qs loop', Qs, Qsmin, Qsmax

             if(Qs.le.Qsmin) cycle
             if(Qs.ge.Qsmax) exit

             width=FullWidthBaryon(resonance_ID,W)
             breitwig=width/pi*W/((W**2-MRes**2)**2+W**2*width**2)  


             k_in=(/ enu,0.,0.,enu /)
             p_in=(/ MNucleon,0.,0.,0. /)
             p_out(0)=(W**2+MNucleon**2+Qs)/(2.*MNucleon)
             k_out(0)=p_in(0)+k_in(0)-p_out(0)
             k_out(3)=(-Qs-ml_out**2+2.*enu*k_out(0))/(2.*enu)
             k_out(1)=sqrt(k_out(0)**2-k_out(3)**2-ml_out**2)
             k_out(2)=0.
             p_out(1)=-k_out(1)
             p_out(2)=0.
             p_out(3)=enu-k_out(3)

             matrixelementneutron=nuMaEl(process_ID,resonance_ID,neutron, k_in, k_out, p_in, p_out,W)
             matrixelementproton=nuMaEl(process_ID,resonance_ID,proton, k_in, k_out, p_in, p_out,W)


             dsig_neut=scaling*1.9732696817**2*10.**10*W/8./pi/(s-MNucleon**2)**2*breitwig*matrixelementneutron
             dsig_prot=scaling*1.9732696817**2*10.**10*W/8./pi/(s-MNucleon**2)**2*breitwig*matrixelementproton


             write(10,'(10g12.5)') enu, W, Qs, dsig_prot, dsig_neut

          end do qsloop

       end do Wloop

    end do energyloop

  end subroutine calcRESdsigdQsdW


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine calcRESdsigdcosthetadEl(resonance_ID,process_ID,flavor_ID)
    use minkowski, only : sp
    implicit none

    integer, intent(in) :: resonance_ID
    integer, intent(in) :: process_ID,flavor_ID
    real :: elepton, wurzel, Mres
    real :: breitwig, width
    real :: enu, costheta, QE_thres
    real :: dsig_neut, dsig_prot
    real :: ml_in,ml_out
    real :: matrixelementproton, matrixelementneutron
    real, dimension(0:3) :: p_in,p_out,k_in,k_out
    real :: phi, W,phasespace,s,t


    MRes=baryon(resonance_ID)%mass

    call setMasses(process_ID,flavor_ID,ml_in,ml_out)

    write(10,*)'# enu   |   costheta   |   elepton   |   dsig on proton   |   dsig on neutron'
    !set starting energy
    enu=0.
    QE_thres=energythreshold(ml_out,MNucleon+Mpion)

    energyloop: do 

       !write(*,*) 'energy loop', enu, QE_thres

       enu=enu+0.5

       if(enu.le.QE_thres) cycle
       if(enu.gt.2.) exit  !go until 2 GeV neutrino energy

       !set starting costheta
       costheta=-1.

       costhetaloop: do

          costheta=costheta+0.1

          if(costheta.gt.1.) exit  

          !set starting El
          elepton=ml_out

          eleptonloop: do

             elepton=elepton+0.005

             if(elepton.ge.enu) exit

             !random generation of azimutal angle
             phi=rn()*2.*pi

             k_in=(/ enu,0.,0.,enu /)
             p_in=(/ MNucleon,0.,0.,0. /)

             !direction of outgoing lepton
             k_out(0)=elepton
             k_out(1)=sqrt(max((1.-costheta**2),0.))*cos(phi)*sqrt(elepton**2-ml_out**2)
             k_out(2)=sqrt(max((1.-costheta**2),0.))*sin(phi)*sqrt(elepton**2-ml_out**2)
             k_out(3)=costheta*sqrt(elepton**2-ml_out**2)

             !momentum of outgoing particle
             p_out=p_in+k_in-k_out

             if(SP(p_out,p_out).le.(Mnucleon+mpion)) cycle !reaction not possible

             W=sqrt(SP(p_out,p_out))

             width=FullWidthBaryon(resonance_ID,W)
             breitwig=width/pi*W/((W**2-MRes**2)**2+W**2*width**2)  

             phasespace=sqrt(elepton**2-ml_out**2)/pi/16./(SP(k_in,p_in))*  breitwig


             matrixelementneutron=nuMaEl(process_ID,resonance_ID,neutron, k_in, k_out, p_in, p_out,W)
             matrixelementproton=nuMaEl(process_ID,resonance_ID,proton, k_in, k_out, p_in, p_out,W)


             dsig_neut=scaling*1.9732696817**2*10.**10*phasespace*matrixelementneutron
             dsig_prot=scaling*1.9732696817**2*10.**10*phasespace*matrixelementproton


             write(10,'(10g12.5)') enu, costheta, elepton, dsig_prot, dsig_neut

          end do eleptonloop

       end do Costhetaloop

    end do energyloop

  end subroutine calcRESdsigdcosthetadEl


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine calcCCpionprodsigma
    implicit none

    integer :: n,n1,l,n2,j,res
    real :: enu, QE_thres,br
    real :: Wmin, Wmax
    real :: dsig_neut, dsig_prot
    real, dimension(3000) :: x,xx,y1,y2,yy1,yy2
    real :: ml_in,ml_out
    real :: W,width,breitwig,s,Qs,Qsmin,Qsmax
    real :: matrixelementproton, matrixelementneutron
    real, dimension(0:3) :: p_in,p_out,k_in,k_out,pion_momentum_out
    integer :: flavor_ID,process_ID,resID
    real, dimension(1:3) :: position
    real :: dummy1,mass_out
    integer :: charge_out,pion_charge_out
    
    integer :: numtry

    real :: sigtotalproton,sigtotalneutron,sigprot_protpipl,signeut_neutpipl,signeut_protpinu
    real :: backgroundppinu,backgroundneutpipl,backgroundneutpipl_final,backgroundppinu_final

    position(:)=1000.

    sigtotalproton=0.
    sigtotalneutron=0.
    sigprot_protpipl=0.
    signeut_neutpipl=0.
    signeut_protpinu=0.
    backgroundppinu=0.
    backgroundneutpipl=0.

    n=3 !integral precision

    process_ID=CC
    flavor_ID=muon

    call setMasses(process_ID,flavor_ID,ml_in,ml_out)

    write(10,*)'# enu   |   sigma on proton   |   sigma on neutron'
    !set starting energy
    enu=0.
    QE_thres=energythreshold(ml_out,MNucleon+Mpion)

    energyloop: do 

       enu=enu+0.1

       if(enu.le.QE_thres) cycle
       if(enu.gt.10.) exit  !go until 10 GeV neutrino energy

       sigtotalproton=0.
       sigtotalneutron=0.
       sigprot_protpipl=0.
       signeut_neutpipl=0.
       signeut_protpinu=0.
       backgroundppinu=0.
       backgroundneutpipl=0.

       do res=2, max_finalstateID

          resID=res

          s=MNucleon**2+2.*enu*MNucleon

          Wmin=MNucleon+Mpion
          Wmax=sqrt(s)-ml_out

          call sg20r(Wmin,Wmax,n,x,n1)
          do j=1,n1
             W=x(j)

             call maxminQs(s,W, ml_out,Qsmax, Qsmin) 

             if(Qsmax.lt.Qsmin) then
                cycle
             end if

             call sg20r(Qsmin,Qsmax,n,xx,n2)
             do l=1,n2
                Qs=xx(l)


                width=FullWidthBaryon(res,W)
                breitwig=width/pi*W/((W**2-baryon(res)%mass**2)**2+W**2*width**2)  

                k_in=(/ enu,0.,0.,enu /)
                p_in=(/ MNucleon,0.,0.,0. /)
                p_out(0)=(W**2+MNucleon**2+Qs)/(2.*MNucleon)
                k_out(0)=p_in(0)+k_in(0)-p_out(0)
                k_out(3)=(-Qs-ml_out**2+2.*enu*k_out(0))/(2.*enu)
                k_out(1)=sqrt(k_out(0)**2-k_out(3)**2-ml_out**2)
                k_out(2)=0.
                p_out(1)=-k_out(1)
                p_out(2)=0.
                p_out(3)=enu-k_out(3)

                matrixelementneutron=nuMaEl(process_ID,resID,neutron, k_in, k_out, p_in, p_out,W)
                matrixelementproton=nuMaEl(process_ID,resID,proton, k_in, k_out, p_in, p_out,W)

                yy1(l)=scaling*1.9732696817**2*10.**10*W/8./pi/(s-MNucleon**2)**2*breitwig*matrixelementneutron
                yy2(l)=scaling*1.9732696817**2*10.**10*W/8./pi/(s-MNucleon**2)**2*breitwig*matrixelementproton


             end do

             call rg20r(Qsmin,Qsmax,n,yy1,dsig_neut)
             call rg20r(Qsmin,Qsmax,n,yy2,dsig_prot)

             y1(j)=dsig_neut
             y2(j)=dsig_prot

          end do

          call rg20r(Wmin,Wmax,n,y1,dsig_neut)
          call rg20r(Wmin,Wmax,n,y2,dsig_prot)

          sigtotalneutron=sigtotalneutron+dsig_neut
          sigtotalproton=sigtotalproton+dsig_prot

          br=baryon(res)%decays2Body(1)

          if(baryon(res)%isoSpinTimes2.eq.3) then
             sigprot_protpipl=sigprot_protpipl+br*dsig_prot
             signeut_neutpipl=signeut_neutpipl+1./3.*br*dsig_neut
             signeut_protpinu=signeut_protpinu+2./3.*br*dsig_neut
          end if

          if(baryon(res)%isoSpinTimes2.eq.1) then 
             signeut_neutpipl=signeut_neutpipl+2./3.*br*dsig_neut
             signeut_protpinu=signeut_protpinu+1./3.*br*dsig_neut
          end if

       end do

       p_in=(/ MNucleon,0.,0.,0. /)


       !BG contribution
       backgroundneutpipl_final=0.
       backgroundppinu_final=0.
       numtry=500
       do i=1,numtry
          call Xsec_integratedSigma(process_ID,flavor_ID,p_in,position,neutron,.false.,33,dummy1, &
               &  k_in,k_out,p_out,mass_out,charge_out,pion_momentum_out,pion_charge_out,backgroundppinu,enu)
          call Xsec_integratedSigma(process_ID,flavor_ID,p_in,position,neutron,.false.,32,dummy1, &
               &  k_in,k_out,p_out,mass_out,charge_out,pion_momentum_out,pion_charge_out,backgroundneutpipl,enu)
          backgroundneutpipl_final=backgroundneutpipl/float(numtry)+backgroundneutpipl_final
          backgroundppinu_final=backgroundppinu/float(numtry)+backgroundppinu_final
       end do


       write(10,'(10g12.4)') enu,sigtotalproton,sigtotalneutron,sigprot_protpipl,signeut_neutpipl,   &
            &  signeut_protpinu,backgroundppinu_final,backgroundneutpipl_final

    end do energyloop

  end subroutine calcCCpionprodsigma


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine calcNCpionprodsigma
    implicit none

    integer :: n,n1,l,n2,j,res
    real :: enu, QE_thres,br
    real :: Wmin, Wmax
    real :: dsig_neut, dsig_prot
    real, dimension(3000) :: x,xx,y1,y2,yy1,yy2
    real :: ml_in,ml_out
    real :: W,width,breitwig,s,Qs,Qsmin,Qsmax
    real :: matrixelementproton, matrixelementneutron
    real, dimension(0:3) :: p_in,p_out,k_in,k_out,pion_momentum_out
    integer :: flavor_ID,process_ID,resID
    real, dimension(1:3) :: position
    real :: dummy1,mass_out
    integer :: charge_out,pion_charge_out

   
    real :: sigtotalproton,sigtotalneutron,sigprot_neutpipl, sigprot_protpinu, signeut_neutpinu, signeut_protpimi

    position(:)=1000.

    sigtotalproton=0.
    sigtotalneutron=0.
    sigprot_neutpipl=0.
    sigprot_protpinu=0.
    signeut_neutpinu=0.
    signeut_protpimi=0.

    n=3 !integral precision

    process_ID=NC
    flavor_ID=muon

    call setMasses(process_ID,flavor_ID,ml_in,ml_out)

    write(10,*)'# enu   |   sigma on proton   |   sigma on neutron'
    !set starting energy
    enu=0.
    QE_thres=energythreshold(ml_out,MNucleon+Mpion)

    energyloop: do 

       enu=enu+0.05

       if(enu.le.QE_thres) cycle
       if(enu.gt.2.2) exit  !go until 2 GeV neutrino energy

       sigtotalproton=0.
       sigtotalneutron=0.
       sigprot_neutpipl=0.
       sigprot_protpinu=0.
       signeut_neutpinu=0.
       signeut_protpimi=0.

       do res=2, max_finalstateID

          resID=res

          s=MNucleon**2+2.*enu*MNucleon

          Wmin=MNucleon+Mpion
          Wmax=sqrt(s)-ml_out

          call sg20r(Wmin,Wmax,n,x,n1)
          do j=1,n1
             W=x(j)

             call maxminQs(s,W, ml_out,Qsmax, Qsmin) 

             if(Qsmax.lt.Qsmin) then
                cycle
             end if

             call sg20r(Qsmin,Qsmax,n,xx,n2)
             do l=1,n2
                Qs=xx(l)


                width=FullWidthBaryon(res,W)
                breitwig=width/pi*W/((W**2-baryon(res)%mass**2)**2+W**2*width**2)  

                k_in=(/ enu,0.,0.,enu /)
                p_in=(/ MNucleon,0.,0.,0. /)
                p_out(0)=(W**2+MNucleon**2+Qs)/(2.*MNucleon)
                k_out(0)=p_in(0)+k_in(0)-p_out(0)
                k_out(3)=(-Qs-ml_out**2+2.*enu*k_out(0))/(2.*enu)
                k_out(1)=sqrt(k_out(0)**2-k_out(3)**2-ml_out**2)
                k_out(2)=0.
                p_out(1)=-k_out(1)
                p_out(2)=0.
                p_out(3)=enu-k_out(3)

                matrixelementneutron=nuMaEl(process_ID,resID,neutron, k_in, k_out, p_in, p_out,W)
                matrixelementproton=nuMaEl(process_ID,resID,proton, k_in, k_out, p_in, p_out,W)

                yy1(l)=scaling*1.9732696817**2*10.**10*W/8./pi/(s-MNucleon**2)**2*breitwig*matrixelementneutron
                yy2(l)=scaling*1.9732696817**2*10.**10*W/8./pi/(s-MNucleon**2)**2*breitwig*matrixelementproton


             end do

             call rg20r(Qsmin,Qsmax,n,yy1,dsig_neut)
             call rg20r(Qsmin,Qsmax,n,yy2,dsig_prot)

             y1(j)=dsig_neut
             y2(j)=dsig_prot

          end do

          call rg20r(Wmin,Wmax,n,y1,dsig_neut)
          call rg20r(Wmin,Wmax,n,y2,dsig_prot)

          sigtotalneutron=sigtotalneutron+dsig_neut
          sigtotalproton=sigtotalproton+dsig_prot

          br=baryon(res)%decays2Body(1)

          if(baryon(res)%isoSpinTimes2.eq.3) then
             sigprot_neutpipl=sigprot_neutpipl+1./3.*br*dsig_prot
             sigprot_protpinu=sigprot_protpinu+2./3.*br*dsig_prot
             signeut_neutpinu=signeut_neutpinu+2./3.*br*dsig_neut
             signeut_protpimi=signeut_protpimi+1./3.*br*dsig_neut
          end if

          if(baryon(res)%isoSpinTimes2.eq.1) then 
             sigprot_neutpipl=sigprot_neutpipl+2./3.*br*dsig_prot
             sigprot_protpinu=sigprot_protpinu+1./3.*br*dsig_prot
             signeut_neutpinu=signeut_neutpinu+1./3.*br*dsig_neut
             signeut_protpimi=signeut_protpimi+2./3.*br*dsig_neut
          end if

       end do

       p_in=(/ MNucleon,0.,0.,0. /)

       write(10,'(10g12.4)') enu,sigtotalproton,sigtotalneutron, &
            & sigprot_neutpipl, sigprot_protpinu, signeut_neutpinu, signeut_protpimi



    end do energyloop

  end subroutine calcNCpionprodsigma


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real function energythreshold(finallepton_mass,finalhadron_minimalmass)
    !note: finalhadron_minimalmass is the minimal mass the final hadron can have
    !e.g. for the delta it is MNucleon+Mpion

    implicit none
    real, intent(in) :: finallepton_mass, finalhadron_minimalmass

    energythreshold=((finallepton_mass+finalhadron_minimalmass)**2-MNucleon**2)/(2.*MNucleon)

  end function energythreshold


  subroutine maxminQs(s, W, ml_out, Qsmax, Qsmin) 

    IMPLICIT NONE      
    REAL, intent(in) ::  s,W,ml_out
    REAL, intent(out) :: Qsmax,Qsmin
    real :: Ecm,Ecmpr

    !real :: wurzel
    !wurzel=max(0.,(s-ml_out**2)**2-2.*(s+ml_out**2)*W**2 + W**4)  !max(0,..) since it might cause problems if W approx. Wmax
    !Qsmax=(2.*enu**2*Mnucleon-Mnucleon*ml_out**2+enu*(MNucleon**2-ml_out**2-W**2)+enu*sqrt(wurzel))/(2.*enu+MNucleon)
    !Qsmin=(2.*enu**2*Mnucleon-Mnucleon*ml_out**2-enu*(-MNucleon**2+ml_out**2+W**2)-enu*sqrt(wurzel))/(2.*enu+MNucleon)

    Ecm=(s-MNucleon**2)/2./sqrt(s)
    Ecmpr=(s-W**2+ml_out**2)/2./sqrt(s)

    Qsmax=-ml_out**2+2.*Ecm*(Ecmpr+sqrt(Ecmpr**2-ml_out**2))
    Qsmin=-ml_out**2+2.*Ecm*(Ecmpr-sqrt(Ecmpr**2-ml_out**2))

  end subroutine maxminQs


  subroutine setMasses(process_ID,flavor_ID,ml_in,ml_out)
    use constants, ONLY: mmuon, melec, mtau
    use neutrino_IDTable
    use leptonicID
    implicit none

    integer, intent(in) :: process_ID
    integer, intent(in) :: flavor_ID
    real, intent(out) :: ml_in,ml_out

    select case (flavor_ID)
    case(muon) 
       ml_in=mmuon
       ml_out=mmuon
    case(electron) 
       ml_in=melec
       ml_out=melec
    case(taulepton) 
       ml_in=mtau
       ml_out=mtau
    case default
       write(*,*) 'unknown flavor, error!'
       stop
    end select

    select case (process_ID)
    case(CC,antiCC) 
       ml_in=0.
    case(NC,antiNC) 
       ml_in=0.
       ml_out=0.
    case(EM,antiEM)

    case default      
       write(*,*) 'CC, NC or EM?? error!'
       stop
    end select

  end subroutine setMasses






  subroutine reinsehgalhelis(Qs,charge,id,fPlus1,fPlus3,fMinus1,fMinus3,f0Plus,f0Minus)
    use idtable, only : nucleon,delta,P11_1440,S11_1535,D13_1520
    use ParticleProperties
    implicit none
    real, intent(in) :: Qs
    integer, intent(in) :: charge, id
    real :: MV = 0.84
    real :: m_target,W
    integer :: n
    real :: ffcorr, gv, lambda,t,r,s, omega,root_half_omega,qq,nu,qsq
    real, intent(out) :: fPlus1,fPlus3,fMinus1,fMinus3,f0Plus,f0Minus
    

    m_target=baryon(Nucleon)%mass
    w=baryon(id)%mass


    nu=(W**2-m_target**2-Qs)/2./m_target

    qq=sqrt(nu**2+Qs)
    qsq=-Qs

    

    omega=1.05
    root_half_omega=sqrt(omega/2.)
    n=2
    if(id.eq.delta) n=0
    if(id.eq.D13_1520.or.id.eq.S11_1535) n=1


    ffcorr=sqrt(1.+Qs/(4.*M_target**2))
    gv=ffcorr**(1-2*n)*(1./(1.+Qs/MV**2))**2
  
    lambda=(m_target*qq)/w/root_half_omega
    t=root_half_omega*gv/(3.*w)
    r=sqrt(2.)*(m_target/w)*(((w+m_target)*qq)/((w+m_target)**2-qsq))*gv
    s=(-qsq/qq**2)*((3.*w*m_target+qsq-m_target**2)/(6.*m_target**2))*gv
  
    if(charge.eq.1) then 
       !EM prot

       select case (id)
          
       case (delta) 
          
         fPlus1  =  Sqrt(2.) * R
         fPlus3  =  Sqrt(6.) * R
         fMinus1 = -1 * fPlus1
         fMinus3 = -1 * fPlus3
         f0Minus =  0.
         f0Plus  =  0.

      case (P11_1440) 

         fMinus1 = -0.5*Sqrt(3.) * Lambda**2 * R
         fPlus1  =  fMinus1
         fMinus3 =  0.
         fPlus3  =  0.
         f0Minus = -0.5*Sqrt(3.) * Lambda**2 * S
         f0Plus  =  f0Minus


      case (S11_1535) 

         fMinus1 =  Sqrt(3.) * T + sqrt(3./2.) * Lambda * R
         f0Minus = -sqrt(3./2.) * Lambda * S
         fPlus1  = -1. * fMinus1
         f0Plus  = -1. * f0Minus
         fMinus3 =  0.
         fPlus3  =  0.


      case (D13_1520) 

         fMinus1 =  sqrt(3./2.) * T - Sqrt(3.) * Lambda * R
         fMinus3 =  3./Sqrt(2.) * T
         f0Minus = -Sqrt(3.) * Lambda * S
         fPlus1  =  fMinus1
         fPlus3  =  fMinus3
         f0Plus  =  f0Minus

      case default
    
         fMinus1 =  0.
         fPlus1  =  0.
         fMinus3 =  0.
         fPlus3  =  0.
         f0Minus =  0.
         f0Plus  =  0.

      end select

   else if(charge.eq.0) then

!EM neut
  select case (id)

      case (delta) 

         fPlus1  =  sqrt(2.) * R
         fPlus3  =  sqrt(6.) * R
         fMinus1 = -1. * fPlus1
         fMinus3 = -1. * fPlus3
         f0Minus =  0.
         f0Plus  =  0.


      case (S11_1535) 

         fPlus1  =  sqrt(3.) * T + 1./sqrt(6.) * Lambda * R
         f0Minus =  sqrt(3./2.) * Lambda * S
         fMinus1 = -1. * fPlus1
         f0Plus  = -1. * f0Minus
         fMinus3 =  0.
         fPlus3  =  0.



      case (D13_1520) 

         fMinus1 = -sqrt(3./2.) * T + sqrt(1./3.) * Lambda * R
         fMinus3 = -sqrt(9./2.) * T
         f0Minus =  sqrt(3.) * Lambda * S
         fPlus1  =  fMinus1
         fPlus3  =  fMinus3
         f0Plus  =  f0Minus


    
      case (P11_1440) 

         fMinus1 = sqrt(1./3.) * Lambda**2 * R
         fPlus1  = fMinus1
         fMinus3 =  0.
         fPlus3  =  0.
         f0Minus =  0.
         f0Plus  =  0.

      case default

         fMinus1 =  0.
         fPlus1  =  0.
         fMinus3 =  0.
         fPlus3  =  0.
         f0Minus =  0.
         f0Plus  =  0.


      end select

   end if



    end subroutine reinsehgalhelis











end program testXsections

