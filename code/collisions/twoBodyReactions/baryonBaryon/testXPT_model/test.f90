program test
  use barBar_to_barBar_model
  use idTable, only: nucleon, delta
  use inputGeneral, only: readinputGeneral
  use particleProperties, only: initParticleProperties, hadron
  use twoBodyTools, only : p_lab
  use constants, only: mPi, mN
  implicit none

  integer, dimension (1:2) :: idIn,chargeIn,idOut,chargeOut
  integer :: i
  real :: srts,cross

  call readinputGeneral
  call initParticleProperties

!  call testDeltaMass
!  call testDeltaWinkel
!  call testDelta
!  call testISO
!  call testNuc
  
  call testDelta_nuc

contains

  subroutine testISO
    ! Tests isospin factors and charge choosing

    integer :: deltaCharge,i
    real, dimension(-1:2,0:1) :: counter
    idIn=(/nucleon,delta/)
    idOUT=(/delta,nucleon/)

    ! Proton
    write(*,*)
    write(*,*) 'PROTON'
    write(*,*)

    do  deltaCharge=-1,2
       counter=0
       chargeIn=(/1,deltaCharge/)
       do i=1,10000
          chargeOUT = ND_ND_choosecharge(idIN,idOUT,chargeIN)
          counter(chargeOUT(1),chargeOUT(2))=counter(chargeOUT(1),chargeOUT(2))+1
       end do
       counter=counter/float(10000)
       write(*,'(A,I3,A,4F10.5)') 'p Delta(',deltaCharge,') -> n + ', counter(:,0)
       write(*,'(A,I3,A,4F10.5)') 'p Delta(',deltaCharge,') -> p + ', counter(:,1)
    end do


    write(*,*) 'Neutron'
    write(*,*)

    do  deltaCharge=-1,2
       counter=0
       chargeIn=(/0,deltaCharge/)
       do i=1,10000
          chargeOUT = ND_ND_choosecharge(idIN,idOUT,chargeIN)
          counter(chargeOUT(1),chargeOUT(2))=counter(chargeOUT(1),chargeOUT(2))+1
       end do
       counter=counter/float(10000)
       write(*,'(A,I3,A,4F10.3)') 'n Delta(',deltaCharge,') -> n + ', counter(:,0)
       write(*,'(A,I3,A,4F10.3)') 'n Delta(',deltaCharge,') -> p + ', counter(:,1)
    end do
  end subroutine testISO


!   subroutine testNuc
!     idIn=(/nucleon,nucleon/)
!     idOUT=(/nucleon,nucleon/)
!     chargeIn=(/1,2/)
!     chargeOut=(/1,2/)
!     open(100,file='barbar_nuc.dat')
!     do i=0,100
!        srts=2*hadron(nucleon)%mass+float(i)*0.01
!        cross=bb_bb(srts,idIn,chargeIn,idOut,chargeOut,chargeSUM=.false.)
!        write(*,*)   srts, p_lab(srts, hadron(nucleon)%mass , hadron(nucleon)%mass), cross
!        write(100,*) srts, p_lab(srts, hadron(nucleon)%mass , hadron(nucleon)%mass), cross
!     end do
!     close(100)
!   end subroutine testNuc

  subroutine testDelta
    idIn=(/nucleon,nucleon/)
    idOUT=(/nucleon,delta/)
    chargeIn=(/1,1/)
    chargeOut=(/1,1/)
    open(100,file='barbar_delta.dat')
    do i=0,300
       srts=2*mN+float(i)*0.01
       cross = NN_ND_model (srts) !,idIn,chargeIn,idOut,chargeOut,chargeSUM=.false.)
       write(*,*)   srts, p_lab(srts, mN , mN), cross
       write(100,*) srts, p_lab(srts, mN , mN), cross
    end do
    close(100)
  end subroutine testDelta


  subroutine testDelta_Nuc
    use constants, only: mN, GeVSquared_times_mb, pi

    real, dimension(1:2) :: mass
    integer:: i,j
    real :: massout, costheta, pdel,edel
    idIn=(/nucleon,delta/)
    idOUT=(/nucleon,delta/)
    chargeIn=(/1,1/)
    chargeOut=(/1,1/)

    open(100,file='deltaNuc_elast_pole.dat')
    write(100,'(A,4I5,A)') '# Reaction:', idIn(1), chargeIn(1),idIn(2), chargeIn(2),&
         & '=> N Delta (summed over outgoing charges!!)'
    write(100,*) '#p_delta incoming , m_delta incoming , cross section'
    mass=(/mN,hadron(delta)%mass/)

    do i=0,300
       pdel=float(i)*0.005
       edel=sqrt(pdel**2+mass(2)**2)
       srts=sqrt((mN+edel)**2-pdel**2)
       cross = ND_ND_model (srts,idIn,chargeIn,idOut,chargeOut,mass,chargeSUM=.true.)
       write(*,'(5G18.4)')   pdel, srts, mass(2), cross, 0.75*0.168*0.197**3*pdel/edel*cross* GeVSquared_times_mb
       write(100,'(5G18.4)') pdel, srts, mass(2), cross, 0.75*0.168*0.197**3*pdel/edel*cross* GeVSquared_times_mb
    end do

    close(100)
    stop 

    if(.true.) then
    open(100,file='deltaNuc_elast_mass_mom.dat')
    write(100,'(A,4I5,A)') '# Reaction:', idIn(1), chargeIn(1),idIn(2), chargeIn(2),&
         & '=> N Delta (summed over outgoing charges!!)'
    write(100,*) '#p_delta incoming , m_delta incoming , cross section'
    do j=0,50
              mass=(/mN,mN+mPi+float(j)*0.005/)
       do i=0,300
          pdel=float(i)*0.005
          edel=sqrt(pdel**2+mass(2)**2)
          srts=sqrt((mN+edel)**2-pdel**2)
          cross = ND_ND_model(srts,idIn,chargeIn,idOut,chargeOut,mass,chargeSUM=.true.)
          write(*,*)   pdel, mass(2), cross
          write(100,*) pdel, mass(2), cross
       end do
       write(100,*)
    end do
    close(100)

       mass=(/mN,hadron(delta)%mass/)
       open(100,file='barbar_deltaNuc_elast.dat')
       write(100,'(A,4I5,A,4I5)') '# Reaction:', idIn(1), chargeIn(1),idIn(2), chargeIn(2),&
            & '=>',idOut(1), chargeOut(1),idOut(2), chargeOut(2)
       write(100,*) '# Mass of incoming Delta=', mass(2)
       do i=0,100
          srts=2*mN+float(i)*0.01
          cross = ND_ND_model(srts,idIn,chargeIn,idOut,chargeOut,mass,chargeSUM=.false.)
          write(*,*)   srts, p_lab(srts,mass(2) , mass(1) ) , cross
          write(100,*) srts, p_lab(srts,mass(2) , mass(1) ) , cross
       end do
       close(100)

       mass=(/mN,hadron(delta)%mass/)
       open(100,file='barbar_deltaNuc_elast_chargeSUM.dat')
       write(100,'(A,4I5,A,4I5)') '# Reaction:', idIn(1), chargeIn(1),idIn(2), chargeIn(2),&
            & '=>',idOut(1), chargeOut(1),idOut(2), chargeOut(2)
       write(100,*) '# Mass of incoming Delta=', mass(2)
       do i=0,100
          srts=2*mN+float(i)*0.01
          cross = ND_ND_model(srts,idIn,chargeIn,idOut,chargeOut,mass,chargeSUM=.true.)
          write(*,*)   srts, p_lab(srts,mass(2) , mass(1) ) , cross
          write(100,*) srts, p_lab(srts,mass(2) , mass(1) ) , cross
       end do
       close(100)


       open(100,file='barbar_deltaNuc_elast_dm.dat')
       write(100,'(A,4I5,A,4I5)') '# Reaction:', idIn(1), chargeIn(1),idIn(2), chargeIn(2),&
            & '=>',idOut(1), chargeOut(1),idOut(2), chargeOut(2)
       do j=0,30
          mass(2)=1.+float(j)*0.015
          do i=0,30
             srts=2*mN+float(i)*0.0125
             cross = ND_ND_model (srts,idIn,chargeIn,idOut,chargeOut,mass,chargeSUM=.false.)
             write(*  ,'(4G18.5)') srts, p_lab(srts,mass(2) , mass(1) ) , mass(2),cross
             write(100,'(4G18.5)') srts, p_lab(srts,mass(2) , mass(1) ) , mass(2),cross
          end do
          write(100,*) 
       end do
       close(100)
    end if

    chargeIn=(/1,1/)
    chargeOut=(/0,2/)
    mass=(/mN,hadron(delta)%mass/)
    srts=2.18
    open(100,file='barbar_deltaNuc_effe')
    write(100,'(A,4I5,A,4I5)') '# Reaction:', idIn(1), chargeIn(1),idIn(2), chargeIn(2),&
         & '=>',idOut(1), chargeOut(1),idOut(2), chargeOut(2)
    write(100,*) '# srts=', srts
    write(100,*) '# Mass of incoming Delta=', mass(2)
    write(100,*) '# mass outgoing, cosTheta, theta [rad], dsigma/dOmega/dMass'
    massOut=1.08
    do 
       massOut=massout+0.02
       do i=0,40
          costheta=-1+float(i)*2./float(40)
          cross=  MSquared_ND_ND(chargeIN,chargeOUT,srts,acos(costheta),mass(2),massout,.false., .false.,chargeSUM=.false.) &
               & /(64.*pi**2*srts**2)/GeVSquared_times_mb
          write(100,'(4G18.5)')   massout,costheta,acos(costheta), cross
       end do
       write(100,*) 
       if(massout.gt.1.4) exit
    end do
    close(100)

  end subroutine testDelta_Nuc


  subroutine testDeltaMass
    use constants, only : pi, GeVSquared_times_mb
    real :: plab,mass
    idIn=(/nucleon,nucleon/)
    idOUT=(/nucleon,delta/)
    chargeIn=(/1,0/)
    chargeOut=(/1,2/)
    plab=2.23
    open(100,file='barbar_delta_dm.dat')
    srts=sqrt((mN+sqrt(mN**2+plab**2))**2-plab**2)
    write(100,*) '# sqrt(s)=',srts
    write(100,*) '# plab   =',plab
    do i=0,100
       mass=mN+mPi+float(i)*0.01
       cross = MSquared_NN_ND(srts,0.,mass,.true., .false.)/(64.*pi**2*srts**2)/GeVSquared_times_mb
       write(*,*)   mass, cross
       write(100,*) mass, cross
    end do
    close(100)
  end subroutine testDeltaMass


  subroutine testDeltaWinkel
    use constants, only : pi, GeVSquared_times_mb
    real :: plab,costheta
    idIn=(/nucleon,nucleon/)
    idOUT=(/nucleon,delta/)
    chargeIn=(/1,0/)
    chargeOut=(/1,2/)
    plab=1.66

    open(100,file='barbar_delta_dcosTheta.dat')
    srts=sqrt((mN+sqrt(mN**2+plab**2))**2-plab**2)
    write(100,*) '# sqrt(s)=',srts
    write(100,*) '# plab   =',plab
    do i=0,100
       costheta=-1.+float(i)*0.02
       cross = MSquared_NN_ND(srts,acos(costheta),99.,.false., .true.)/(64.*pi**2*srts**2)/GeVSquared_times_mb*2.*pi
       write(100,*)   costheta,acos(costheta), cross
    end do
    close(100)
  end subroutine testDeltaWinkel


end program test
