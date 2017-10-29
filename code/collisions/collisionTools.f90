!******************************************************************************
!****m* /collisionTools
! NAME
! module collisionTools
! PURPOSE
! Some helper routines for collisions
!******************************************************************************
module collisionTools

  implicit none
  private

  public:: finalCheck
  public:: setEnergyCheck

  logical, save :: initFlag=.true.
  real, save :: energyCheck=0.01 ! copy of collisionTerm/energyCheck

contains



  !****************************************************************************
  !****s* collisionTools/setEnergyCheck
  ! NAME
  ! subroutine setEnergyCheck(va)
  ! PURPOSE
  ! set the internal variable energyCheck
  !****************************************************************************
  subroutine setEnergyCheck(val)
    real, intent(in) :: val
    energyCheck = val
    initFlag=.false.
  end subroutine setEnergyCheck


  !****************************************************************************
  !****f* collisionTools/finalCheck
  ! NAME
  ! function finalCheck(partIn, partOut, HiEnergyFlagge, woher) result(flag)
  ! PURPOSE
  ! Checks the final state of a collision for the conservation of all
  ! quantum numbers,  also including momentum and energy conservation.
  !
  ! For HiEnergy events we do not check charge and momentum conservation,
  ! this MUST be done separately.
  ! The reason for this is, that some particles which are produced by
  ! Pythia/Fritiof can not be propagated by BUU and therefore do not show up in
  ! the final state vector "partOut"
  ! ("unknown particles wont be propagated").
  !
  ! INPUTS
  ! * type(particle),dimension(:)  :: partIn  -- Incoming particles
  ! * type(particle),dimension(:)  :: partOut -- Outgoing particles
  ! * logical, optional            :: HiEnergyFlag --
  !   .true. if it was a HiEnergy event.
  !   if .true. then energy conservation is not checked
  !   and code does not stop if charge conservation is violated
  ! * character(*), optional      :: woher -- ...
  !
  ! OUTPUT
  ! * logical :: flag -- .true. if quantum numbers are conserved
  !****************************************************************************
  function finalCheck(partIn, partOut, HiEnergyFlag, woher) result(flag)

    use particleDefinition
    use IdTable, only: isBaryon, isHadron
    use particleProperties, only: hadron, validCharge
    use callStack, only: traceBack
    use output, only: WriteParticle

    type(particle), intent(in), dimension(:) :: partIn, partOut
    logical, intent(in), optional :: HiEnergyFlag
    character(*), intent(in), optional :: woher
    logical :: flag

    integer :: totalCharge_in, totalCharge_out, baryon_number_In, baryon_number_Out, strangeness_In, strangeness_Out, k
    real, dimension(0:3) :: momentum_In, momentum_Out
    real :: m_miss2  ! ,s_In, s_Out

    flag=.false.

    ! Initialize at first call
    if (initFlag) call Traceback('finalCheck: not initialized!')

    ! Check charge of each particle
    do k=lbound(partOut,dim=1),ubound(partOut,dim=1)
       if (.not.validCharge(partOut(k))) then
          if (present(woher)) write(*,*) 'Problem in ',trim(woher),' !!!!!!'
          write(*,*) 'Charge of particle is not valid in finalCheck'
          call WriteParticle(6,1,k,partOut(k))
          write(*,*) 'Charge=',partOut(k)%charge
          write(*,*) 'Id=',partOut(k)%id
          write(*,*) 'Antiparticle=',partOut(k)%antiparticle
          write(*,*) 'Momentum=',partOut(k)%momentum
          write(*,*)
          if (present(HiEnergyFlag)) then
             if (HiEnergyFlag) then
                write(*,*) 'It was a HiEnergy event!'
                call PYLIST(2)
             else
                write(*,*) 'It was a low-energy event!'
             end if
          else
             write(*,*) 'It was a decay-event!'
          end if
          write(*,*)
          write(*,'(79("*"))')
          call printPart(partIn, .true.)
          call printPart(partOut, .false.)
          call writeParticle(6,1,partIn)
          call writeParticle(6,2,partOut)
          call Traceback('Severe problem. Code stops!!!!')
       end if
    end do

    if (present(HiEnergyFlag)) then
       if (HiEnergyFlag) then
          flag = .true.
          return
       end if
    end if

    ! Check baryon number conservation:
    baryon_number_In = 0
    do k=lbound(partIn,dim=1),ubound(partIn,dim=1)
       if (isBaryon(partIn(k)%Id)) then
          if (.not.partIn(k)%antiparticle) then
             baryon_number_In = baryon_number_In + 1
          else
             baryon_number_In = baryon_number_In - 1
          end if
       end if
    end do
    baryon_number_Out = 0
    do k=lbound(partOut,dim=1),ubound(partOut,dim=1)
       if (isBaryon(partOut(k)%Id)) then
          if (.not.partOut(k)%antiparticle) then
             baryon_number_Out = baryon_number_Out + 1
          else
             baryon_number_Out = baryon_number_Out - 1
          end if
       end if
    end do
    if (baryon_number_In /= baryon_number_Out) then
       if (present(woher)) then
          write(*,*) 'Problems in '//trim(woher)//' : Baryon number conservation'
       else
          write(*,*) 'Problems in collisionTerm: Baryon number conservation'
       end if
       write(*,*) 'Baryon number in: ', baryon_number_In
       write(*,*) 'Baryon number out:', baryon_number_Out
       write(*,'(79("*"))')
       call printPart(partIn, .true.)
       call printPart(partOut, .false.)
       call writeParticle(6,1,partIn)
       call writeParticle(6,2,partOut)
       call Traceback('stop in finalCheck, baryon number check')
    end if

    ! Check strangeness conservation:
    strangeness_In = 0
    do k=lbound(partIn,dim=1),ubound(partIn,dim=1)
       if (.not. isHadron(partIn(k)%Id)) cycle
       if (.not. partIn(k)%antiparticle) then
          strangeness_In = strangeness_In + hadron(partIn(k)%Id)%strangeness
       else
          strangeness_In = strangeness_In - hadron(partIn(k)%Id)%strangeness
       end if
    end do
    strangeness_Out = 0
    do k=lbound(partOut,dim=1),ubound(partOut,dim=1)
       if (.not. isHadron(partOut(k)%Id)) cycle
       if (.not. partOut(k)%antiparticle) then
          strangeness_Out = strangeness_Out + hadron(partOut(k)%Id)%strangeness
       else
          strangeness_Out = strangeness_Out - hadron(partOut(k)%Id)%strangeness
       end if
    end do
    if (strangeness_In /= strangeness_Out) then
       if (present(woher)) then
          write(*,*) 'Problems in '//trim(woher)//' : Strangeness conservation'
       else
          write(*,*) 'Problems in collisionTerm: Strangeness conservation'
       end if
       write(*,*) 'Strangeness in: ', strangeness_In
       write(*,*) 'Strangeness out:', strangeness_Out
       write(*,'(79("*"))')
       call printPart(partIn, .true.)
       call printPart(partOut, .false.)
       call writeParticle(6,1,partIn)
       call writeParticle(6,2,partOut)
       !  Check conservation of srts
       do k=0,3
          momentum_IN(k) =Sum(partIn%momentum(k))
          momentum_Out(k)=Sum(partOut%momentum(k))
       end do
       m_miss2 = (momentum_IN(0)-momentum_Out(0))**2 - (momentum_IN(1)-momentum_Out(1))**2 &
            -(momentum_IN(2)-momentum_Out(2))**2 - (momentum_IN(3)-momentum_Out(3))**2
       write(*,*) ' (Missing mass)^2 :', m_miss2
       call Traceback('stop in finalCheck, strangeness check')
    end if


    !  Check total charge
    totalCharge_In  =   Sum(partIn%charge)
    totalCharge_Out =   Sum(partOut%charge)
    if (totalCharge_In /= totalCharge_Out) then
       if (present(woher)) then
          write(*,*) 'Problems in '//trim(woher)//' : Charge conservation'
       else
          write(*,*) 'Problems in collisionTerm: Charge conservation'
       end if
       write(*,*) 'Charge in: ', totalCharge_In
       write(*,*) 'Charge out:', totalCharge_Out
       write(*,'(79("*"))')
       call printPart(partIn, .true.)
       call printPart(partOut, .false.)
       call writeParticle(6,1,partIn)
       call writeParticle(6,2,partOut)
       call Traceback('stop in finalCheck, charge check')
    end if

    !       if (debug) Print *, 'CollisionTerm: SqrtS=',sqrtS(partIN(1),partIN(2))

    !  Check conservation of srts
    do k=0,3
       momentum_IN(k) =Sum(partIn%momentum(k))
       momentum_Out(k)=Sum(partOut%momentum(k))
    end do

    !s_In=momentum_IN(0)**2-Dot_Product(momentum_IN(1:3), momentum_IN(1:3))
    !s_Out=momentum_Out(0)**2-Dot_Product(momentum_Out(1:3), momentum_Out(1:3))

    do k=0,3
       if (abs(momentum_In(k)-momentum_Out(k)).gt.energyCheck) then
          if (present(woher)) then
             write(*,*) 'Problems in '//trim(woher)//' : Energy/momentum conservation'
          else
             write(*,*) 'Problems in collisionTerm: Energy/momentum conservation'
          end if
          write(*,'(A,4G13.5)')'Momentum in: ', Momentum_In
          write(*,'(A,4G13.5)')'Momentum out:', Momentum_Out
          write(*,*)
          write(*,'(79("*"))')
          call printPart(partIn, .true.)
          call printPart(partOut, .false.)
          call writeParticle(6,1,partIn)
          call writeParticle(6,2,partOut)
          call Traceback('stop in finalCheck, momentum check')
       end if
    end do

    flag=.true.

  contains

    subroutine printPart(part, in)
      use mediumModule, only: mediumAt
      use output, only: writeParticle_debug
      type(particle), intent(in), dimension(:) :: part
      logical, intent(in) :: in
      integer :: i
      do i=lbound(part,dim=1) , ubound(part,dim=1)
         if (part(i)%ID==0) cycle
         if (in) then
            write(*,*) 'Incoming Particle number #', i
         else
            write(*,*) 'Outgoing Particle number #', i
         end if
         call writeParticle_debug(part(i),mediumAt(part(i)%position))
      end do
      write(*,'(79("*"))')

    end subroutine printPart

  end function finalCheck

end module collisionTools
