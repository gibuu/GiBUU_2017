program test

  use inputGeneral, only: readInputGeneral
  use output, only: intToChar
  use constants, only: pi
  use mediumDefinition
  use mesonWidthMedium_tables, only : get_inMediumWidth,getMaxMes,getMinMes,getMaxAbsP
  use mesonWidthMedium, only: WidthMesonMedium
  use mesonWidth, only: FullWidthMeson
  use particleProperties, only: initParticleProperties, hadron
  implicit none

  real,parameter :: rho=0.17
  real,parameter :: dp=0.02,dm=0.02
  integer :: j,k,ID,jmax,kmax
  real :: pabs,mass,Gtot,mres0,spec,x_off,dGdp,dGdm,E,Gcoll
  real :: mom(0:3)
  type(medium) :: med

  call initParticleProperties
  call readInputGeneral

  med%densityProton=rho/2
  med%densityNeutron=rho/2
  med%useMedium=.true.

  ! write out coll. width, total width and spectral function as a function of mass and momentum at density rho_0=0.17fm^(-3)
  ! also: offshell parameter and dGamma/dp
  do ID=getMinMes(),getMaxMes()
    write(*,*) "ID=",ID
    Open(57,File='Width_coll_'//trim(intToChar(ID))//'.dat')
    Open(58,File='Width_tot_'//trim(intToChar(ID))//'.dat')
    Open(59,File='Spectral_'//trim(intToChar(ID))//'.dat')
    Open(60,File='OSParam_'//trim(intToChar(ID))//'.dat')
    Open(61,File='dGammadp_'//trim(intToChar(ID))//'.dat')
    Open(62,File='dGammadm_'//trim(intToChar(ID))//'.dat')
    mres0=hadron(ID)%mass
    kmax=nint(getMaxAbsP(ID)/dp)
    jmax=nint(3.0/dm)
    write(*,*) "kmax=",kmax
    do j=1,jmax
      mass=float(j)*dm
      do k=0,kmax
        pabs=float(k)*dp
        E=sqrt(mass**2+pabs**2)
        mom=(/E,pabs,0.,0./)
        Gcoll = mass*get_inMediumWidth(ID,pabs,mass,med)
        Gtot = mass*WidthMesonMedium(ID,mass,mom,med)
        spec = 2./pi*mass**2*Gtot/((mass**2-mres0**2)**2+Gtot**2*mass**2)
        x_off=(mass-mres0)/Gtot
        if (k==0) then
          dGdp=(mass*WidthMesonMedium(ID,mass,(/E,pabs+dp,0.,0./),med)-mass*WidthMesonMedium(ID,mass,(/E,pabs,0.,0./),med))/dp
        else if (k==kmax) then
          dGdp=(mass*WidthMesonMedium(ID,mass,(/E,pabs,0.,0./),med)-mass*WidthMesonMedium(ID,mass,(/E,pabs-dp,0.,0./),med))/dp
        else
          dGdp=(mass*WidthMesonMedium(ID,mass,(/E,pabs+dp,0.,0./),med) - &
                mass*WidthMesonMedium(ID,mass,(/E,pabs-dp,0.,0./),med))/(2*dp)
        end if
        if (j==1) then
          dGdm=((mass+dm)*WidthMesonMedium(ID,mass+dm,mom,med)-mass*WidthMesonMedium(ID,mass,mom,med))/dm
        else if (j==jmax) then
          dGdm=(mass*WidthMesonMedium(ID,mass,mom,med)-(mass-dm)*WidthMesonMedium(ID,mass-dm,mom,med))/dm
        else
          dGdm=((mass+dm)*WidthMesonMedium(ID,mass+dm,mom,med)-(mass-dm)*WidthMesonMedium(ID,mass-dm,mom,med))/(2*dm)
        end if
        write(57,'(3E18.6)') mass,pabs,Gcoll
        write(58,'(3E18.6)') mass,pabs,Gtot
        write(59,'(3E18.6)') mass,pabs,spec
        write(60,'(3E18.6)') mass,pabs,x_off
        write(61,'(3E18.6)') mass,pabs,dGdp
        write(62,'(3E18.6)') mass,pabs,dGdm
      end do
      write(57,*)
      write(58,*)
      write(59,*)
      write(60,*)
      write(61,*)
      write(62,*)
    end do
    close(57)
    close(58)
    close(59)
    close(60)
    close(61)
    close(62)

    ! write out collisional width for onshell mass, as function of momentum
    Open(63,File='Width_coll_onshell_'//trim(intToChar(ID))//'.dat')
    mass=mres0
    do k=0,kmax
      pabs=float(k)*dp
      E=sqrt(mass**2+pabs**2)
      mom=(/E,pabs,0.,0./)
      Gcoll = get_inMediumWidth(ID,pabs,mass,med)
      write(63,'(2E18.6)') pabs,Gcoll
    end do
    close(63)
  end do

end program test
