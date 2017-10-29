!******************************************************************************
!****m* /energyCalc
! NAME
! module energyCalc
! PURPOSE
! Incorporates routines to evaluate the one-particle energy of a given
! particle and a routine for making a so called energy correction to a given
! final state.
!******************************************************************************
module energyCalc
  implicit none

  private

  !****************************************************************************
  !****g* energyCalc/accuracy
  ! PURPOSE
  ! Determines accuracy of energy determination in the determination of
  ! the one particle energy. Units : GeV.
  ! SOURCE
  !
  real, parameter :: accuracy=0.0005
  !****************************************************************************


  public :: energyDetermination
  public :: energyCorrection
  public :: updateEnergies

contains


  !****************************************************************************
  !****s* energyCalc/updateEnergies
  ! NAME
  ! subroutine updateEnergies(part,check_in)
  ! PURPOSE
  ! Routine updates energies of the particles
  ! INPUTS
  ! * type(particle), dimension(:) :: part  -- particle array
  ! * logical, OPTIONAL :: check_in -- flag for 'energyDetermination'
  ! OUTPUT
  ! * type(particle), dimension(:) :: part  -- particle array
  !****************************************************************************
  subroutine updateEnergies(part,check_in)
    use particleDefinition
    use output, only: DoPR
    type(particle),intent(inOut), dimension(:,:) :: part
    logical,intent(in),optional :: check_in
    integer :: i,j
    logical :: check

    check = .false.
    if (present(check_in)) check = check_in

    if (DoPr(2)) write(*,*) 'Updating energies of Particles'

    do i=1,Size(part,dim=1)
       do j=1,Size(part,dim=2)
          if (part(i,j)%ID < 0 ) exit
          if (part(i,j)%ID == 0 ) cycle
          call energyDetermination(part(i,j),check=check)
       end do
    end do

  end subroutine updateEnergies


  !****************************************************************************
  !****s* energyCalc/energyDetermination
  ! NAME
  ! subroutine energyDetermination(part,betaToCF,countMax_in,error_out,
  ! check,warn,ForbidCoulomb,skipTestE0)
  !
  ! PURPOSE
  ! This subroutine determines the one-particle energy E with a
  ! momentum-dependent scalar potential that is defined in the local
  ! rest frame.
  !
  ! NOTES
  ! * We keep p(1:3) and m constant.
  ! * To solve: p**2=E**2-p(1:3)**2=(m+s(E))**2
  ! * The scalar potential s is only defined in the LRF.
  ! * First we try to find the solution f(E)=0 via a
  !   the newton method,
  !     f(E)=m_eff**2-p**2=m_eff**2+p(1:3)**2-E**2
  !   with
  !     m_eff=m+s(E)
  !   by the iteration
  !     E(i+2)=E(i+1)-f(i+1)*(E(i+1)-E(i))/(f(i+1)-f(i)).
  ! * If the newton method fails, we try it again via a
  !   bisection method
  ! * At the end, the Coulomb potential is added.
  !
  ! INPUTS
  ! * type(particle) :: part --  Particle to be considered.
  ! * real,dimension (1:3), OPTIONAL :: betaToCF -- Velocity of
  !   "Calculation Frame" in the frame where the particle's momentum is
  !    defined.
  ! * integer, OPTIONAL :: countMax_in  -- maximal number of cycles to find
  !   a solution (default=10000)
  ! * logical, OPTIONAL :: check -- check for energy conservation
  !   (Default: .false.)
  ! * logical, OPTIONAL :: warn -- if true, a warning message is issued,
  !   if the incomming momentum is too large (otherwise, this is skipped
  !   silently). (Default: .true.)
  ! * logical, OPTIONAL :: ForbidCoulomb -- if true, the Coulomb potential
  !   is NOT added to p0 (needed by 'propagate', where the coulomb force
  !   is treated seperately) (Default: .false.)
  ! * logical, OPTIONAL :: skipTestE0 -- if true, the input value of E will
  !   not be tested to be a already a good solution
  ! OUTPUT
  ! * type(particle) :: part
  ! * real, OPTIONAL    :: error_out    -- approximate error of E = sqrt(f(e))
  !****************************************************************************
  subroutine energyDetermination(part, betaToCF, countMax_in, error_out, check, warn, ForbidCoulomb, skipTestE0)
    use IdTable, only: rho,omegaMeson,phi,isMeson
    use particleDefinition
    use minkowski, only: abs4
    use offShellPotential, only: treatParticleOffShell, get_offshell_debug
    use output, only: WriteParticle_debug!, WriteParticle
    use hist2Df90
    use callStack, only: traceback
    use dichteDefinition
    use densityModule, only: densityAt
    use coulomb, only: emfoca, getCoulombFlag
    use lorentzTrafo, only: lorentzCalcBeta, lorentz

    type(particle), intent(inOut) :: part
    real, optional, dimension (1:3), intent(in) :: betaToCF
    integer, optional, intent(in) :: countMax_in
    real, optional, intent(out)   :: error_out
    logical, optional, intent(in) :: check, warn, ForbidCoulomb, skipTestE0

    integer, save :: failures=0, numtry=0  !,hugeCounter=0
    real :: mass,psqr,E0, E_min, Cpot
    real,dimension(1:3) :: f,E  ! f=f(E), need 1:3 for iteration
    logical :: CFBoostFlag      ! Determines whether it's necessary to boost to Calculation Frame
    integer :: counter,countmax
    integer,parameter :: countmax_default=100
    real,parameter :: maxError=accuracy**2   !wished accuracy of f(E) in GeV**2
    ! --> Energy is correct within (maxError)**(1/2)
    logical,parameter :: debugflag=.false.
    real,parameter :: delta_E = 0.01
    logical :: checkE,success,outOfBounds_offshell,warn_in,DoCoulomb,skipE0
    type(dichte) :: density
!    real, dimension(0:3) :: mm
!    real, dimension(1:3) :: bb
    real :: f0 = 99.999 ! dummy

    warn_in=.true.
    if (present(warn)) warn_in=warn

!     If(AbsMom(part).gt.10000) then
!        !*********************************************************************
!        ! This is to explore some problem which seems to occur in  the
!        ! energyCorrection loops:
!        ! When no solution can be found, the regula falsi is producing a X to
!        ! scale the momentum which is way to high!!
!        ! This should be no problem to finished jobs. It just crashes the code
!        ! from time to time since we can not find a solution for 100 TEV
!        ! particles in this routine here.
!        ! Oliver
!        !*********************************************************************
!        if(warn_in) then
!           write(*,*)
!           write(*,*) 'WARNING: Particle with huge momentum in EnergyDetermination'
!           call errorMessage
!           hugeCounter=hugeCounter+1
!           write(*,*) 'It is the ',hugeCounter, 'th time'
!        end if
!        part%momentum(0)=FreeEnergy(part)
!        return
!        write(*,*)
!     end if

!!$    If(AbsMom(part).gt.1000) then
!!$       write(*,*) "energyDetermination: Particle with large momentum: ",AbsMom(part)
!!$       call WriteParticle(6,0,0,part)
!!$       call traceback("calllist",-1)
!!$    endif


    countMax=countMax_default
    if (present(countMax_in))  countMax=countMax_in

    CFBoostFlag=.false.
    if (present(betaToCF)) CFBoostFlag=.true.

    checkE=.false.
    if (present(check)) checkE=check

    DoCoulomb = .true.
    if (present(ForbidCoulomb)) DoCoulomb = .not.ForbidCoulomb

    skipE0 = .false.
    if (present(skipTestE0)) skipE0 = skipTestE0

    !----------------------------------------------

    E0=part%momentum(0)

    mass=part%mass
    psqr=Sum(part%momentum(1:3)**2)
    density = densityAt(part%position)

    if (.not.skipE0 .or. checkE) then
       f0 = m_eff_sqr(part,density)-(E0**2-psqr)
    end if
    if (.not.skipE0) then
       if (abs(f0) <= maxError) then
          if (DebugFlag) write(*,*) 'E0 already okay!'
          if (present(error_out)) error_out=sqrt(abs(f0))
          return
       end if
    end if

    if (getCoulombFlag()) then
      E_min = max(0.,sqrt(psqr)-0.2) ! just a guess (potentials can be negative)
    else
      E_min = sqrt(psqr+0.001**2)
    end if

!!$    if (DebugFlag) then
!!$       if (part%ID==101) then
!!$          write(124,'(1P,4e15.7)') part%position
!!$          write(124,'(1P,4e15.7)') part%momentum
!!$          mm = part%momentum
!!$          if (present(betaToCF)) then
!!$             write(124,'(A,1P,4e15.7)') 'betaToCF: ',betaToCF
!!$             call lorentz(betaToCF, mm)
!!$             write(124,'(1P,4e15.7)') mm
!!$          end if
!!$
!!$          if (density%baryon(0)>1E-8 .and. sum(density%baryon(1:3)**2)>1E-8) then
!!$             bb = lorentzCalcBeta(density%baryon, 'm_eff_sqr')
!!$             write(124,'(A,1P,4e15.7)') 'betaToLRF:',bb
!!$             call lorentz(bb, mm)
!!$             write(124,'(1P,4e15.7)') mm
!!$          end if
!!$
!!$          write(124,*)
!!$          write(124,*)
!!$          flush(124)
!!$       end if
!!$    end if


    if (isMeson(part%ID) &
         .and. treatParticleOffShell(part%ID,part%OffShellParameter) &
         .and. checkE &
         .and. E0>sqrt(psqr)) then

       ! special treatment for offshell mesons, since f(E)=0 may have more
       ! than one solution:
       ! find a solution to f(E)=0, restricting E to E0 +- a few MeV

       if (get_offshell_debug()) then
          call MesonHistFailures(1)
       end if

       numtry=numtry+1

       success=bisection(E0-delta_E,E0+delta_E)

       if (.not. success) then
          failures=failures+1
          if (get_offshell_debug()) then !  .and. checkE
             part%momentum(0)=E0
             call MesonHistFailures(2)
          end if
          write(*,*) "WARNING in EnergyDetermination: Impossible to determine energy for offshell meson!"
          write(*,*) "Number of Failures: ",failures,' (',float(failures)/float(numtry)*100., '%)'
          write(*,*)
          !call plot_f(part)
          part%ID = 0
          !call traceback('Error in energyDetermination',0)
       end if

    else

       ! default case

       outOFBounds_offShell=.false.
       success=newton()
       if (outOFBounds_offShell) success=vary_mass_offshell()
       if (.not. success) success=bisection()
       if (.not. success) then
          write(*,*) "WARNING in EnergyDetermination: Impossible to determine energy!"
          write(*,*)
          if (warn_in) then
             call plot_f(part)
             call traceback('Error in energyDetermination')
          end if
       end if

    end if

    if (checkE .and. abs(part%momentum(0)-E0)>delta_E) then
       write(*,*) "WARNING in EnergyDetermination: Energy not conserved!",part%momentum(0)-E0
       write(*,*) 'E0 = ', E0
       write(*,*) 'f0 = ', f0, abs(f0)/maxError
       if (present(betaToCF)) then
          write(*,'(A,1P,4e15.7)') 'betaToCF: ',betaToCF
       end if
       call WriteParticle_debug(part)

!       call plot_f(part)
       !stop
    end if


  contains

    !**************************************************************************
    !****is* energyDetermination/MesonHistFailures
    ! NAME
    ! subroutine MesonHistFailures(iTyp)
    !
    ! PURPOSE
    ! Do statistics of failures connected with mesons and off-shellness
    !**************************************************************************
    subroutine MesonHistFailures(iTyp)
      integer, intent(in) :: iTyp

      logical, save :: firstTime=.true.
      type(histogram2D),save :: hist2D_rho,hist2D_omega,hist2D_phi

      select case (iTyp)
      case (1)
         if (firstTime) then
            call createHist2D(hist2D_rho,&
                 'energy violation as function of mass and momentum',&
                 (/0.,0./),(/2.,2./),(/0.02,0.02/))
            call createHist2D(hist2D_omega,&
                 'energy violation as function of mass and momentum',&
                 (/0.,0./),(/2.,2./),(/0.02,0.02/))
            call createHist2D(hist2D_phi,&
                 'energy violation as function of mass and momentum',&
                 (/0.,0./),(/2.,2./),(/0.02,0.02/))
            firstTime=.false.
         end if

         select case (part%ID)
         case (rho)
            call AddHist2D(hist2D_rho,(/abs4(part%momentum),absMom(part)/),0.,1.)
         case (omegaMeson)
            call AddHist2D(hist2D_omega,(/abs4(part%momentum),absMom(part)/),0.,1.)
         case (phi)
            call AddHist2D(hist2D_phi,(/abs4(part%momentum),absMom(part)/),0.,1.)
         end select

      case (2)

         select case (part%ID)
         case (rho)
            call AddHist2D(hist2D_rho,(/abs4(part%momentum),absMom(part)/),1.)
         case (omegaMeson)
            call AddHist2D(hist2D_omega,(/abs4(part%momentum),absMom(part)/),1.)
         case (phi)
            call AddHist2D(hist2D_phi,(/abs4(part%momentum),absMom(part)/),1.)
         end select

         if (mod(failures,100)==0) then
            call writeHist2D_Gnuplot(hist2D_rho,44,file='Energy2D_rho.dat')
            call writeHist2D_Gnuplot(hist2D_omega,45,file='Energy2D_omega.dat')
            call writeHist2D_Gnuplot(hist2D_phi,46,file='Energy2D_phi.dat')
         end if

      end select

    end subroutine MesonHistFailures

    !**************************************************************************
    !****if* energyDetermination/vary_mass_offshell
    ! NAME
    ! logical function vary_mass_offshell()
    !
    ! PURPOSE
    ! Try maximizing p_0 such that equation is still fulfilled
    !**************************************************************************
    logical function vary_mass_offshell()
      use particleDefinition
      type(particle) :: part2
      real, parameter :: delta_pNull=0.0001
      real :: error
      part2=part
      vary_mass_offshell=.true.
      do
         part2%momentum(0)=part2%momentum(0)+delta_pNull
         error=m_eff_sqr(part2,density)-(part2%momentum(0)**2-psqr)
         if (abs(error).lt.maxError) then !  Successful!!!
            if (DebugFlag) write(*,*) 'SuccessFull with f(2)',f(2),E(2)
            if (present(error_out)) error_out=sqrt(abs(f(2)))
            return
         else
            part2%momentum(0)=part2%momentum(0)-delta_pNull
            exit
         end if
         if (part2%momentum(0).gt.10.) then
            ! FAILURE
            write(*,*) 'FAILURE in vary_mass_offshell'
            vary_mass_offshell=.false.
            exit
         end if
      end do

    end function vary_mass_offshell

    !**************************************************************************
    !****if* energyDetermination/newton
    ! NAME
    ! logical function newton()
    !
    ! PURPOSE
    ! solve f(E)=0 by Newton Method
    !**************************************************************************
    logical function newton()

      newton=.true.

      if (DebugFlag) write(*,*) '**In EnergyDetermination (Newton)'

      !Make Vacuum Ansatz to evaluate starting values
      E(1)=SQRT(psqr+mass**2)
      part%momentum(0)=E(1)

      f(1)=m_eff_sqr(part,density)-(E(1)**2-psqr)
      if (abs(f(1)).lt.maxError) then !  Successful!!!
         if (DebugFlag) write(*,*) '                 E0: ',0.0,E0
         if (DebugFlag) write(*,*) 'SuccessFull with f(1)',f(1),E(1)
         if (present(error_out)) error_out=sqrt(abs(f(1)))
         return
      end if

      E(2)=sqrt(max(f(1)+E(1)**2,E_min**2))
      part%momentum(0)=E(2)
      f(2)=m_eff_sqr(part,density)-(E(2)**2-psqr)

      if (abs(f(2)).lt.maxError) then !  Successful!!!
         if (DebugFlag) write(*,*) '                 E0: ',0.0,E0
         if (DebugFlag) write(*,*) '                 f(1)',f(1),E(1)
         if (DebugFlag) write(*,*) 'SuccessFull with f(2)',f(2),E(2)
         if (present(error_out)) error_out=sqrt(abs(f(2)))
         return
      end if

      if (debugFlag) write(*,*) 'Begin Iteration'

      ! Begin Iteration

      counter=0
      do
         if (abs(f(2)-f(1)).eq.0) then
            if (DebugFlag) then
               write(*,*) 'Error in EnergyDetermination (Newton). Derivative too small:', f(1:2),f(2)-f(1)
               call errorMessage
               write(*,*) 'Counter:' , counter
            end if
            newton=.false.
            return
         end if

         E(3)=E(2)-f(2)*(E(2)-E(1))/(f(2)-f(1))
         part%momentum(0)=E(3)
         f(3)=m_eff_sqr(part,density)-(E(3)**2-psqr)

         if (debugFlag) write(11,*) e(3),f(3)

         if (abs(f(3)).lt.maxError) then !Successful!!!
            if (DebugFlag) write(*,*) '                 E0: ',0.0,E0
            if (DebugFlag) write(*,*) 'SuccessFull with f(3)',f(3),E(3)
            if (DebugFlag) write(*,*) 'counter=',counter
            if (present(error_out)) error_out=sqrt(abs(f(3)))
            if (part%momentum(0)<0.) newton=.false.
            return
         end if

         if (counter.gt.countmax) then
            if (DebugFlag) then
               write(*,*) 'Error in Energydetermination (Newton):  counter.gt.countmax'
               call errorMessage
               call plot_f(part)
            end if
            newton=.false.
            return
         end if
         counter=counter+1
         E(1)=E(2)
         E(2)=E(3)
         f(1)=f(2)
         f(2)=f(3)
      end do

    end function newton

    !**************************************************************************
    !****if* energyDetermination/bisection
    ! NAME
    ! logical function bisection(E1,E2)
    !
    ! PURPOSE
    ! solve f(E)=0 by Bisection Method
    !**************************************************************************
    logical function bisection(E1,E2)

      real,intent(in),optional :: E1,E2 ! initial boundaries

      bisection=.true.

      if (DebugFlag) write(*,*) '**In EnergyDetermination (Bisection)'

      if (present(E1) .and. Present(E2)) then

         E(1)=max(E_min,E1)
         part%momentum(0)=E(1)
         f(1)=m_eff_sqr(part,density)-(E(1)**2-psqr)

         if (abs(f(1)).lt.maxError) then !  Successful!!!
            if (DebugFlag) write(*,*) 'SuccessFull with f(1)',f(1),E(1)
            if (present(error_out)) error_out=sqrt(abs(f(1)))
            return
         end if

         E(2)=E2
         part%momentum(0)=E(2)
         f(2)=m_eff_sqr(part,density)-(E(2)**2-psqr)

         if (abs(f(2)).lt.maxError) then !  Successful!!!
            if (DebugFlag) write(*,*) 'SuccessFull with f(2)',f(2),E(2)
            if (present(error_out)) error_out=sqrt(abs(f(2)))
            return
         end if

         ! E1 and E2 are given, but not bracketing a zero. Maybe the range
         ! was just too big?
         counter=0
         do while (f(1)*f(2)>0. .and. counter<20)
            E(3)=(E(1)+E(2))/2.  ! half interval
            part%momentum(0)=E(3)
            f(3)=m_eff_sqr(part,density)-(E(3)**2-psqr)
            if (abs(f(3)).lt.maxError) then !  Successful!!!
               if (DebugFlag) write(*,*) 'SuccessFull with f(3)',f(3),E(3)
               if (present(error_out)) error_out=sqrt(abs(f(2)))
               return
            end if
            if (f(1)*f(3)<0.) then
               f(2)=f(3)
               E(2)=E(3)
            else if (f(2)*f(3)<0.) then
               f(1)=f(3)
               E(1)=E(3)
            else if (abs(f(1))<abs(f(2))) then
               f(2)=f(3)
               E(2)=E(3)
            else if (abs(f(2))<abs(f(1))) then
               f(1)=f(3)
               E(1)=E(3)
            else
               f(mod(counter,2)+1)=f(3)
               E(mod(counter,2)+1)=E(3)
            end if
            counter=counter+1
         end do

      else

         !if (E0==0. .and. debugFlag) write(*,*) 'EnergyDetermination (Bisection): initializing E!'

         ! set starting values: E(1) and E(2)
         E(1)=E_min
         part%momentum(0)=E(1)
         f(1)=m_eff_sqr(part,density)-(E(1)**2-psqr)

         if (abs(f(1)).lt.maxError) then !  Successful!!!
            if (DebugFlag) write(*,*) 'SuccessFull with f(1)',f(1),E(1)
            if (present(error_out)) error_out=sqrt(abs(f(1)))
            return
         end if

         E(2)=E(1)+0.001
         part%momentum(0)=E(2)
         f(2)=m_eff_sqr(part,density)-(E(2)**2-psqr)

         if (abs(f(2)).lt.maxError) then !  Successful!!!
            if (DebugFlag) write(*,*) 'SuccessFull with f(2)',f(2),E(2)
            if (present(error_out)) error_out=sqrt(abs(f(2)))
            return
         end if

         ! increase E(2)
         counter=0
         do while (f(1)*f(2)>0. .and. counter<20)
            E(2)=2*E(2)-E(1)  ! double interval
            if (DebugFlag) write(*,*) 'increasing E(2): ',E(2)
            part%momentum(0)=E(2)
            f(2)=m_eff_sqr(part,density)-(E(2)**2-psqr)
            if (abs(f(2)).lt.maxError) then !  Successful!!!
               if (DebugFlag) write(*,*) 'SuccessFull with f(2)',f(2),E(2)
               if (present(error_out)) error_out=sqrt(abs(f(2)))
               return
            end if
            counter=counter+1
         end do

      end if


      if (f(1)*f(2)>0) then
         write(*,*) "Error in EnergyDetermination (Bisection): Impossible to find starting values!"
         if (present(E1) .and. Present(E2)) write(*,*) E1,E2
         call errorMessage
         bisection=.false.
         return
         !call plot_f(part)
         !Stop 'Stop in EnergyDetermination (Bisection)'
      end if

      ! Begin Iteration
      !If (DebugFlag) write(*,*) "Bisection: starting iteration!"

      counter=0
      do
         if (f(1)*f(2)>0) then
            write(*,*) 'Error in EnergyDetermination (Bisection):  Same Sign!'
            call errorMessage
            call plot_f(part)
            Stop 'Stop in EnergyDetermination (Bisection)'
         end if

         E(3)=(E(1)+E(2))/2
         part%momentum(0)=E(3)
         f(3)=m_eff_sqr(part,density)-(E(3)**2-psqr)

         if (abs(f(3))<maxError .or. E(2)-E(1)<accuracy) then !Successful!!!
            if (DebugFlag) write(*,*) 'SuccessFull with f(3)',f(3),E(3)
            if (present(error_out)) error_out=sqrt(abs(f(3)))
            return
         end if

         if (counter>40) write(*,*) E,f

         if (counter.gt.countmax) then
            write(*,*) 'Error in EnergyDetermination (Bisection):  counter.gt.countmax'
            call errorMessage
            bisection=.false.
            return
            !call plot_f(part)
            !Stop 'Stop in EnergyDetermination (Bisection)'
         end if
         counter=counter+1

         if (f(1)*f(3)<0.) then
            if (DebugFlag) write(*,'(A)',ADVANCE='NO') 'L'
            E(2)=E(3)
            f(2)=f(3)
         else if (f(2)*f(3)<0.) then
            if (DebugFlag) write(*,'(A)',ADVANCE='NO') 'U'
            E(1)=E(3)
            f(1)=f(3)
         else
            write(*,*) 'Error in EnergyDetermination (Bisection): ZERO!'
            write(*,*) 'counter = ',counter
            call errorMessage
            call plot_f(part)
            call traceback('Stop in EnergyDetermination (Bisection)')
         end if
      end do

    end function bisection

    !**************************************************************************
    !****is* energyDetermination/plot_f
    ! NAME
    ! subroutine plot_f(part_in)
    !
    ! PURPOSE
    ! This subroutine plots the function
    !   m_eff**2(part,ener)-(Energy**2-Sum(p**2)) as a
    ! function of energy. Zoomes into the region where the function is zero.
    ! Useful for debugging!!
    !**************************************************************************
    subroutine plot_f(part_in)
      use output

      type(particle),intent(in) :: part_in
      type(particle) :: part2
      real :: dE
      real :: resu,energy,me2
      real :: down,ener!,old
      integer :: i!,j,ende
      integer :: max=10000

      part2=part_in
      down=E_min
      dE=1./float(max)
      open(111,file='plot_f.dat')
      !do j=1,5 ! Loop which zooms into the region of the zero
      !ende=max
      do i=0,max ! Loop which steps through the energy
         energy=down+float(i)*dE
         part2%momentum(0)=energy
         me2 = m_eff_sqr(part2,density,ener)
         resu= me2 - (energy**2-psqr)
         !if(i.eq.0) old=resu
         write(111,'(4E15.7)') part2%momentum(0),resu,me2,ener
         !if(old*resu.lt.0) then
         !   ! Change of sign in the last dE step
         !   down=energy-dE
         !   ende=i+10 ! stop plotting after the next 10 steps
         !end if
         !old=resu
         !if(i.gt.ende) exit
      end do
      !dE=dE/float(max) ! Decrease step size during the zooming process
      !end do
      close(111)
    end subroutine plot_f

    !**************************************************************************
    !****is* energyDetermination/errorMessage
    ! NAME
    ! subroutine errorMessage
    !
    ! PURPOSE
    ! Print some info about particle
    !**************************************************************************
    subroutine errorMessage
      use output, only:  WriteParticle_debug

      write(*,'(a,3E20.7)')'SQRT( f(1:3) )       :', SQRT(abs(f))
      write(*,'(a,3E20.7)')'f(1:3)               :', f
      write(*,'(a,3E20.7)')'E(1:3)               :', E
      write(*,'(a,i4)')    'Id of particle       :', part%Id
      write(*,'(a,i4)')    'Charge of particle   :', part%charge
      write(*,'(a,E12.4)') 'Mass of particle     :', part%mass
      write(*,'(a,4E12.4)')'Momentum of particle :', part%momentum
      write(*,'(a,E12.4)') 'Abs. mom. of particle:', absMom(part)
      write(*,'(a,L2)')    'Perturbative         :', part%perturbative
      write(*,'(a,E12.4)') 'Offshell parameter   :', part%offshellParameter
      if (present(betaToCF)) then
         write(*,'(a,3E12.4)')'beta to Calculation frame :', betaToCF
      else
         write(*,*) 'optional parameter "betaToCF" not used'
      end if

      call  WriteParticle_debug(part)

    end subroutine errorMessage

    !**************************************************************************
    !****if* energyDetermination/m_eff_sqr
    ! NAME
    ! real function m_eff_sqr(partIN,density,enerOut)
    ! PURPOSE
    ! Evaluate m_eff**2=(m+S)**2 of a given particle by boosting to the
    ! LRF-Frame, where S is the scalar potential.
    !**************************************************************************
    real function m_eff_sqr(partIn, density, enerOut)
      use particleDefinition
      use potentialModule, only: potential_LRF, massDetermination
      use lorentzTrafo, only: lorentzCalcBeta, lorentz
      use offShellPotential, only: hamiltonFunc_offshell
      use minkowski, only: SP

      type(particle),intent(in):: partIn
      type(dichte),intent(in) :: density
      type(particle) :: part
      real, optional :: enerOut
      real, dimension(1:3) :: betaToLRF

      part=partIn ! local copy

      ! (1.1) Boost to calculation frame
      if (CFBoostFlag) call lorentz(betaToCF,part%momentum, 'm_eff_sqr')
      ! (1.2) Boost from calculation frame to LRF
      if (density%baryon(0)>1E-8 .and. sum(density%baryon(1:3)**2)>1E-8) then
         betaToLRF = lorentzCalcBeta(density%baryon, 'm_eff_sqr')
         call lorentz(betaToLRF,part%momentum,'m_eff_sqr(2)')
      else
         betaToLRF=0.
      end if

      ! (1.3) Evaluate scalar potential V_S:
      !       (m+V_S)**2=p(0)**2-p(1:3)**2=p**2 =>
      if (treatParticleOffShell(part%ID,part%OffShellParameter)) then
         call massDetermination(part,success=success,verbose=.false.)
         outOfBounds_offshell=.false.
         if (.not.success .or. SP(part%momentum,part%momentum)<0.) then
            ! BAD Solution !!
            part%momentum(0)=999999999.
         else
            part%momentum(0)=HamiltonFunc_offshell(part,outOfBounds_offshell,.true.)
         end if
      else
         part%momentum(0)=FreeEnergy(part)+potential_LRF(part,density,addCoulomb=.false.)
      end if

      if (DoCoulomb) then

         Cpot = emfoca(part%position,part%momentum(1:3),part%charge,part%ID)

         ! This Lorentz transformation to LRF is needed because Cpot from
         ! emfoca is defined in CF:
         Cpot = Cpot  * sqrt(1.-sum(betaToLRF(1:3)**2))

         part%momentum(0)=part%momentum(0)+Cpot
      end if

      if (present(enerOut)) enerOUT=part%momentum(0)

      m_eff_sqr=part%momentum(0)**2-Dot_product(part%momentum(1:3),part%momentum(1:3))

    end function m_eff_sqr

  end subroutine energyDetermination



  !****************************************************************************
  !****s* energyCalc/energyCorrection
  ! NAME
  ! subroutine energyCorrection(srqts, betaToLRF, betaToCM, mediumAtCollision,
  ! finalState, successFlag, potentialFailure, verbose)
  !
  ! PURPOSE
  ! Gets the finalState particles with vacuum kinematics in the CM Frame.
  ! Also the real srts of the final state is handed over to the routine.
  ! Then this routine tries to solve energy and momentum conservation by the
  ! ansatz described in
  ! |html  <a href="../../Documentation_Extra/crossSections/Xsections/node25.html"> the cross section documentation. </a>.
  ! Therefore it evaluates the scaling factor "x" which is used to scale the
  ! CM-momenta of the final states. This is done by  a Regula-Falsi-Method.
  ! We search of a zero of the function "deltaSrts" which depends on this
  ! scaling factor x.
  !
  ! IMPORTANT :
  ! Output are the final state particles in the calculation frame!
  ! INPUTS
  ! * real, intent(in) :: srtS
  ! * real, intent(in), dimension(1:3) :: betaToLRF -- velocity of LRF - NOT USED right now!!!
  ! * real, intent(in), dimension(1:3) :: betaToCM  -- velocity of CM frame in calculation frame
  ! * type(medium), intent(in) :: mediumAtCollision
  ! * type(particle), dimension(:), intent (INOUT) :: finalState -- In CM frame
  ! * logical, OPTIONAL :: verbose -- flag print messages on/off (default: on)
  ! OUTPUT
  ! * type(particle), dimension(:), intent (INOUT) :: finalState -- In calculation frame
  ! * logical, intent(out) :: successFlag
  ! * logical, optional :: potentialFailure -- .true. if it fails due to neglecting perturbative potentials
  ! NOTES
  ! If the iteration fails, successFlag is set to .false. as output.
  ! Now is used also in the RMF mode (only function deltaSrts is modified).
  !****************************************************************************
  subroutine energyCorrection(srts, betaToLRF, betaToCM, mediumAtColl, &
       & finalState, successFlag, potentialFailure, verbose)

    use mediumDefinition, only: medium
    use particleDefinition, only: particle,sqrts
    use LorentzTrafo, only: lorentz
    use random, only: rn
    use mesonPotentialModule, only: getNoPertPot_meson
    use baryonPotentialModule, only: getNoPertPot_baryon
    use RMF, only: getRMF_flag

    real, intent(in) :: srtS
    real, intent(in), dimension(1:3) :: betaToLRF, betaToCM
    type(medium), intent(in) :: mediumAtColl
    type(particle), dimension(:), intent (inout)  :: finalState
    logical, intent(out) :: successFlag
    logical, optional,intent(out) :: potentialFailure
    logical, OPTIONAL, intent(in)    :: verbose

    integer, parameter :: maxLoops=10  ! default maximal number of tries to find a solution for x

    ! local variables
    type(particle), allocatable, dimension(:):: part
    integer :: i, counterEqual
    real, dimension (-1:1) :: f, x

    logical, parameter :: debug=.false.
    real, parameter :: accuracySiter=1E-10
    logical :: verb

    successFlag=.false.
    if (present(potentialFailure)) potentialFailure=.false.

    verb = .true.
    if (present(verbose)) verb = verbose

    if (.not.getRMF_flag()) then
       if (getNoPertPot_meson().and.getNoPertPot_baryon()) then
          ! If perturbative potentials are switched off and finalState is
          ! perturbative,  then we can only initialize particles above the
          ! mass threshold.
          ! Using the following if-statement endles loops are prohibited and
          ! computing time shall be saved.
          if (finalState(1)%perturbative.and.(Sum(finalstate%mass)>srts)) then
             if (present(potentialFailure)) potentialFailure=.true.
             return
          end if
       end if
    end if

    allocate(part(1:size(finalState,dim=1)))

    ! Initialize
    x(0)=1.
    f(0)=deltaSrts( x(0) )


    if (debug) write(*,*) 'f(0)=', f(0)
    if (debug) write(*,*) 'srst=', srts

    if (abs(f(0)).lt. accuracySiter) then
       call setFinalState
       successFlag=.true.
       return
    end if

    x(1)=0.9 ! first wild guess for starting point
    f(1)=deltaSrts( x(1) )

    if (debug) write(*,*) 'f(1)=', f(1)

    if (abs(f(1)).lt. accuracySiter) then
       call setFinalState
       successFlag=.true.
       return
    end if

    if (abs(f(1)-f(0)).lt.1E-6) then
       ! Try to find new starting point for the regula-falsi-method
       counterEqual=0
       findX_loop : do
          x(1)=0.5+rn()
          f(1)=deltaSrts( x(1) )
          if (abs(f(1)-f(0)).gt.1E-6) exit findX_loop
          if (counterEqual.gt.10) then
             if (verb) then
                write(*,*) ' Problem in energyCorrection . f(1) = f(0) :'
                write(*,*) f(0:1)
                write(*,*) ' x=',x(0:1)
                write(*,*) ' Hence iteration not possible'
             end if
             successFlag=.false.
             deallocate(Part)
             return
          end if
          counterEqual=counterEqual+1
       end do findX_loop
    end if


    ! Start Iteration
    do i=1, maxLoops
       ! Store results of last loop
       x(-1:0)=x(0:1)
       f(-1:0)=f(0:1)
       ! Evaluate new guess
       x(1)=x(0)-f(0)*(x(0)-x(-1))/(f(0)-f(-1))
       if (debug) write(*,*) 'Loop', i, 'x=',x(1)
       if (x(1).lt.0) then
          x(1)=-x(1)/float(maxLoops)
          ! /float(maxLoops) because if there is a solution then probably
          ! close to 0
       end if
       f(1)=deltaSrts(x(1))
       if (debug) write(*,*) 'Loop', i, 'f(x)=',f(1)
       if (abs(f(1)).lt.accuracySiter) then
          call setFinalState
          successFlag=.true.
          return
       else if (f(1)-f(0).eq.0.) then
          if (verb) then
             write(*,*) ' WARNING : In energycorrection f(1)-f(0).eq.0',f(0:1),x(0:1)
             write(*,*) 'Hence iteration not possible'
          end if
          successFlag=.false.
          deallocate(Part)
          return
       end if
    end do


    !    If(.not.successFlag) write(*,*) 'energyCorrection:  Energy correction failed'

  contains

    !**************************************************************************
    !****is* energyCorrection/setFinalState
    ! NAME
    ! subroutine setFinalState
    ! PURPOSE
    ! boost the particles and set them into 'finalState'
    !**************************************************************************
    subroutine setFinalState
      integer :: i
      do i=1,size(part,dim=1)
         call lorentz(-betaToCM,part(i)%momentum(0:3), 'energyDetermination,setFinalState')  ! Boost to Calculation frame
      end do
      finalState = part
      if (debug) write(*,*) 'Setting finalState. Sqrts=',sqrtS(finalState(1),finalState(2))
      deallocate(Part)
    end subroutine setFinalState


    !**************************************************************************
    !****if* energyCorrection/deltaSrts
    ! NAME
    ! real function deltaSrts(x)
    ! PURPOSE
    ! Evaluates sqrt(s)-Sum(Energies in CM frame).
    ! This value is the size of the energy non-conservation.
    ! The sum depends on the scaling factor for the three momenta
    !**************************************************************************
    real function deltaSrts(x)

      use densitymodule, only: energyDeterminationRMF

      real, intent(in) :: x
      integer :: i

      part=finalState

      do i=1, size( finalState, dim=1 )
         part(i)%momentum(1:3)=x * part(i)%momentum(1:3)
      end do

      if (debug) write(*,*) 'moms vorher =' , &
           & sum(part%momentum(0)),sum(part%momentum(1)), &
           & sum(part%momentum(2)),sum(part%momentum(3))

      do i=1, size( finalState, dim=1 )
         if (debug) write(*,*) 'part ',i,'vorher :',part(i)%mass, part(i)%momentum,part(i)%perturbative
         if ( .not.getRMF_flag() ) then
            call energyDetermination(part(i),-betaToCM,skipTestE0=.true.)
         else
            call energyDeterminationRMF(part(i))
         end if
         if (debug) write(*,*) 'part ',i,'nachher :',part(i)%mass, part(i)%momentum
      end do

      deltaSrts=srts-Sum(part%momentum(0))

      if (debug) then
         write(*,*) 'real srst=', srts
         write(*,*) 'moms nachher=', &
              & sum(part%momentum(0)),sum(part%momentum(1)), &
              & sum(part%momentum(2)),sum(part%momentum(3))
         write(*,*) 'delta=',deltaSrts
         write(*,*) '-------------------------------------------------------------------------------------------------'
      end if
    end function deltaSrts


  end subroutine energyCorrection





end module energyCalc
