!***********************************************************************
! Test case for BB -> BYK
!***********************************************************************
program testbarBar_barHypKaon

  use inputGeneral, only: readinputGeneral
  use particleProperties, only: initParticleProperties, hadron, PartName
  use particleDefinition
  use IdTable, only: nucleon, Delta
  implicit none

  type(particle), dimension(1:2) :: teilchenIN
  integer :: ID1 = nucleon, charge1 = 1, ID2 = Delta, charge2 = 0
  real :: sigma_BYK(1:18)

  call readinputGeneral
  call initParticleProperties

  write(*,*) '##########################################################'
  write(*,*) 'testing the module barBar_BarHypKaon.f90'
  write(*,*) '##########################################################'

  call get_InputBYK

  teilchenIn(1)%ID=ID1
  teilchenIN(1)%charge=charge1
  teilchenIn(1)%mass=hadron(ID1)%mass
  teilchenIn(1)%momentum(1:3)=0.
  teilchenIn(1)%momentum(0)=Sqrt(teilchenIn(1)%mass**2+sum(teilchenIn(1)%momentum(1:3)**2))

  teilchenIn(2)%ID=ID2
  teilchenIN(2)%charge=charge2
  teilchenIn(2)%mass=hadron(ID2)%mass
  teilchenIn(2)%momentum(1:3)=0.
  teilchenIn(2)%momentum(0)=Sqrt(teilchenIn(2)%mass**2+sum(teilchenIn(2)%momentum(1:3)**2))

  call write_XS
  call check_charge

 contains


  subroutine get_InputBYK
    integer :: ios
    NAMELIST /inputBYK/ ID1,charge1,ID2,charge2
    rewind(5)
    read(5,nml=inputBYK,iostat=ios)
  end subroutine get_InputBYK


  subroutine write_XS
    use barBar_BarHypKaon, only: barBar_barBarMeson_strange
    integer :: i
    real :: srts

    Open(100,File=trim(PartName(teilchenIn(1)))//trim(PartName(teilchenIn(2)))//'_BarHypKaon_XS.dat')
    write(100,*) '# srts, pLab, sigma BB-> BYK'

    !write the xsections (mb) for all final isospin channels
    Do i=1,600
      teilchenIn(2)%momentum(1) = i*0.025
      teilchenIn(2)%momentum(0) = sqrt(teilchenIn(2)%mass**2 + sum(teilchenIn(2)%momentum(1:3)**2))
      srts=sqrtS(teilchenIn(1),teilchenIn(2))
      sigma_BYK = barBar_barBarMeson_strange (srts, teilchenIN)
      write(100,'(18E13.4)') srts, teilchenIn(2)%momentum(1), sum(sigma_BYK), sigma_BYK(1:15)
  !     write(*,'(18E13.4)') srts, teilchenIn(2)%momentum(1), sum(sigma_BYK), sigma_BYK(1:15)
    end Do
    close(100)
  end subroutine write_XS


  subroutine check_charge
    use preEventDefinition
    use barBar_BarHypKaon, only: get_Channels_BYK
    integer :: i,TotalCharge_In,TotalCharge_Out,CheckCharge
    type(preEvent), dimension(1:3) :: PartOut

    !check charge conservation
    Open(101,File=trim(PartName(teilchenIn(1)))//trim(PartName(teilchenIn(2)))//'_BarHypKaon_chargeCons.dat')
    write(101,'(2A)') 'Incoming particle #1: ',PartName(teilchenIn(1))
    write(101,'(2A)') 'Incoming particle #2: ',PartName(teilchenIn(2))
    write(*,'(2A)') 'Incoming particle #1: ',PartName(teilchenIn(1))
    write(*,'(2A)') 'Incoming particle #2: ',PartName(teilchenIn(2))
    write(101,'(A)') 'Possible final channels (BB->BYK):'
    write(*,'(A)') 'Possible final channels (BB->BYK):'

    TotalCharge_In = Sum(teilchenIn%charge)

    do i=1,Size(sigma_BYK)
  !     write(101,'(7i4)') i,(get_Channels_BYK(i,j),j=1,3)
      PartOut(1:3) = get_Channels_BYK(i)
      TotalCharge_Out = Sum(PartOut%charge)
      CheckCharge     = TotalCharge_In - TotalCharge_Out
      if (PartOut(1)%ID==0 .or. PartOut(1)%ID==0) then
         write(101,'(A,i3,A)') 'Channel',i,':  ---'
         write(*,  '(A,i3,A)') 'Channel',i,':  ---'
      else
         write(101,'(A,i3,5A,i3)') 'Channel',i,':  ', &
              PartName (PartOut(1)%ID, PartOut(1)%charge, PartOut(1)%antiParticle), &
              PartName (PartOut(2)%ID, PartOut(2)%charge, PartOut(2)%antiParticle), &
              PartName (PartOut(3)%ID, PartOut(3)%charge, PartOut(3)%antiParticle), &
              'Charge conservation: ',CheckCharge
         write(*,'(A,i3,5A,i3)') 'Channel',i,':  ', &
              PartName (PartOut(1)%ID, PartOut(1)%charge, PartOut(1)%antiParticle), &
              PartName (PartOut(2)%ID, PartOut(2)%charge, PartOut(2)%antiParticle), &
              PartName (PartOut(3)%ID, PartOut(3)%charge, PartOut(3)%antiParticle), &
              'Charge conservation: ',CheckCharge
      endif
    end do

    close(101)

  end subroutine


end program testbarBar_barHypKaon
