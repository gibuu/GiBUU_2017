!******************************************************************************
!****m* /deuterium_PL
! module deuterium_PL
!
! PURPOSE
! This module defines a pointerList to store informations of a deuterium initialization.
! We must remember which particles belong together and this information is then used in
! the module "baryonPotential".
!
! NOTES
! This whole module can not be used in:
!  * full-ensemble mode
!  * real-particle mode
!******************************************************************************
module deuterium_PL

  use particleDefinition, only: particle
  implicit none
  private


  !****************************************************************************
  !****t* deuterium_PL/type_deuteriumPL
  !
  ! SOURCE
  !
  type, public :: type_deuteriumPL
     type(particle), pointer :: part1 ! first particle
     type(particle), pointer :: part2 ! second particle
  end type type_deuteriumPL
  !
  !****************************************************************************


  !****************************************************************************
  !****t* deuterium_PL/deuterium_pointerList
  ! PURPOSE
  ! List of all pairs of particles.
  !
  ! SOURCE
  !
  type(type_deuteriumPL), allocatable, dimension(:), public :: deuterium_pointerList
  !
  !****************************************************************************


  type(particle), pointer, public :: deuterium_pertOrigin
  integer, public                 :: deuterium_pertOrigin_flag
  integer, public                 :: deuteriumPL_ensemble


  public :: get_DeuteriumPotential
  public :: deuteriumPL_clear, deuteriumPL_getSpectra, deuteriumPL_inUse, deuteriumPL_assign

contains


  !****************************************************************************
  !****f* deuterium_PL/deuteriumPL_inUse
  ! NAME
  ! logical function deuteriumPL_inUse()
  !
  ! PURPOSE
  ! Returns true if deuterium_pointerList is allocated,
  ! unless we're working in real-particle mode.
  !****************************************************************************
  logical function deuteriumPL_inUse()
    use inputGeneral, only: eventtype
    use eventtypes, only: HeavyIon, hadron
    deuteriumPL_inUse = (eventtype/=HeavyIon) .and. (eventtype/=hadron) .and. allocated(deuterium_pointerList)
  end function deuteriumPL_inUse


  !****************************************************************************
  !****s* deuterium_PL/deuteriumPL_assign
  ! NAME
  ! subroutine  deuteriumPL_assign(realP)
  !
  ! PURPOSE
  ! Assigns the nucleons in the particle vector "realP" to the pointerList
  ! "deuterium_pointerList". Checks that there are exactly 2 nucleons in each
  ! ensemble. Does not work in fullensemble mode!!!
  !
  ! INPUTS
  ! type(particle), dimension(:,:),target :: realP
  !****************************************************************************
  subroutine deuteriumPL_assign(realP)
    use IdTable, only: nucleon
    use particleDefinition
    use output, only: DoPr, writeParticle_debug
    use callStack, only: traceback
    type(particle), dimension(:,:),target :: realP
    integer :: i,j,numNucs

    if (.not. deuteriumPL_inUse()) return

    if (DoPr(1)) write(*,*) 'DEUTERIUM RUN: Assign particles to pointerlist...'

    do i=lbound(realP,dim=1),ubound(realP,dim=1)
      numNucs=0
      do j=lbound(realP,dim=2),ubound(realP,dim=2)
        if (realP(i,j)%ID /= nucleon) cycle
        numNucs=numNucs+1
        select case (numNucs)
        case (1)
          deuterium_pointerList(i)%part1 => realP(i,j)
        case (2)
          deuterium_pointerList(i)%part2 => realP(i,j)
        case default
          write(*,*) 'ERROR: More than two nucleons in deuterium ! ensemble #', i
          call writeParticle_debug (deuterium_pointerList(i)%part1)
          call writeParticle_debug (deuterium_pointerList(i)%part2)
          call writeParticle_debug (realP(i,j))
          call traceback()
        end select
      end do
      if (numNucs<2) then
        write(*,*) 'ERROR: Less than two nucleons in deuterium ! ensemble #', i, "count=", numNucs
        do j=lbound(realP,dim=2),ubound(realP,dim=2)
          if (realP(i,j)%ID > 0) call writeParticle_debug (realP(i,j))
        end do
        call traceback()
      end if
    end do

    if (DoPr(1)) write(*,*) '... done'

  end subroutine deuteriumPL_assign


  !****************************************************************************
  !****s* deuterium_PL/deuteriumPL_getSpectra
  ! NAME
  ! subroutine deuteriumPL_getSpectra(pl,file_mom,file_dist,time)
  !
  ! PURPOSE
  ! Generates output based on the pointer list "pl". I.e.:
  ! * Binding energy
  ! * Relative momentum and position spectra
  !
  ! INPUTS
  ! * type (type_deuteriumPL),dimension(:) :: pl
  ! * integer :: file_mom,file_dist  -- file identifiers for momentum/position histograms
  ! * real :: time    -- time
  !****************************************************************************
  subroutine deuteriumPL_getSpectra(pl,file_mom,file_dist,time)
    use histf90
    use vector, only: absVec
    use argonneV18, only: argonne_deuteriumPot
    use constants, only: mN
    use particleDefinition, only: freeEnergy

    type (type_deuteriumPL),dimension(:) :: pl
    type(histogram) :: momentumHist, distanceHist
    integer :: file_mom,file_dist
    ! set density to dynamic
    !real :: abs_mom,abs_dist
    real :: integral
    integer :: i
    real :: bindingEnergy_squares,bindingEnergy_error,bindingEnergy,absDist
    real :: radius,radius_squares,radius_error
    real :: rms,rms_squares,rms_error
    real, dimension (0:3) :: totalMom
    real :: time
    logical :: firstFlag=.true.

    if (allocated(deuterium_pointerlist)) then
       if (firstFlag) then
          open(99,file='Binding_deuterium.dat')
          firstFlag=.false.
       else
          open(99,file='Binding_deuterium.dat',position='append')
       end if

       integral=0.
       bindingEnergy=0.
       BindingEnergy_squares=0.
       totalMom=0.
       radius=0.
       radius_squares=0.
       rms=0.
       rms_squares=0.
       call CreateHist(momentumHist, 'momentum distribution',0.,0.5,0.002)
       call CreateHist(distanceHist, 'distance distribution',0.,15.,0.05)
       do i= lbound(pl,dim=1),ubound(pl,dim=1)
          call AddHist(momentumHist,absVec(pl(i)%part1%momentum(1:3)-pl(i)%part2%momentum(1:3)),1.)
          absDist=absVec(pl(i)%part1%position(1:3)-pl(i)%part2%position(1:3))
          call AddHist(distanceHist,absDist,1.)
          radius=radius+absDist
          radius_squares=radius_squares+absDist**2
          rms=rms+absDist**2
          rms_squares=rms_squares+absDist**4
          integral=integral+1.
          totalMom=totalMom+pl(i)%part1%momentum+pl(i)%part2%momentum
          bindingEnergy=bindingEnergy+freeEnergy(pl(i)%part1)+freeEnergy(pl(i)%part2)+argonne_deuteriumPot(absDist)&
               &- 2* mN
          bindingEnergy_squares=bindingEnergy_squares+(freeEnergy(pl(i)%part1)&
               & +freeEnergy(pl(i)%part2)+argonne_deuteriumPot(absDist)&
               &- 2* mN)**2
       end do

       call WriteHist(momentumHist,file_Mom,mul=1./integral)
       call WriteHist(distanceHist,file_Dist,mul=1./integral)

       call removeHist(momentumHist)
       call removeHist(distanceHist)
       bindingEnergy=bindingEnergy/Integral
       radius=radius/Integral
       rms=rms/integral
       if (Integral-1.gt.0.1) then
          bindingEnergy_Error=sqrt(1./Integral/(Integral-1.)*(BindingEnergy_squares-Integral*BindingEnergy**2))
          radius_Error=sqrt(1./Integral/(Integral-1.)*(radius_squares-Integral*radius**2))
          rms_error=sqrt(1./Integral/(Integral-1.)*(rms_squares-Integral*rms**2))
       else
          bindingEnergy_Error=999.
       end if
       rms=sqrt(rms)
       rms_error=1./2.*rms_error/rms
       totalMom=totalMom/Integral
       write(99,'(11G18.4)') time, bindingEnergy,bindingEnergy_Error, rms, rms_Error,radius, radius_Error, totalMom
       write(*,'(2(A,G12.3))') 'Deuterium binding energy=', bindingEnergy, ' +/- ',bindingEnergy_Error
       write(*,'(2(A,G12.3))') 'Deuterium radius        =', radius/2.      , ' +/- ',radius_Error/2.
       write(*,'(2(A,G12.3))') 'Deuterium radius (RMS)  =', rms/2.      , ' +/- ',rms_Error/2.
       close(99)
    end if

  end subroutine deuteriumPL_getSpectra



  !****************************************************************************
  !****s* deuterium_PL/deuteriumPL_clear
  ! NAME
  ! subroutine deuteriumPL_clear(pl)
  !
  ! PURPOSE
  ! Clears an entry of type type_deuteriumPL
  !
  ! INPUTS
  ! * type (type_deuteriumPL) :: pl
  !****************************************************************************
  subroutine deuteriumPL_clear(pl)
    type (type_deuteriumPL) :: pl
    nullify(pl%part1)
    nullify(pl%part2)
  end subroutine deuteriumPL_clear


  !****************************************************************************
  !****s* deuterium_PL/get_DeuteriumPotential
  ! NAME
  ! function get_DeuteriumPotential(p) result(pot)
  !
  ! PURPOSE
  ! Evaluates the potential energy of a pair containing particle p
  !
  ! INPUTS
  ! type(particle),target :: p
  !
  ! OUTPUT
  ! real :: pot -- Argonne V18 potential
  !****************************************************************************
  function get_DeuteriumPotential(p) result(pot)
    use particleDefinition
    use vector, only: absVec
    use argonneV18, only: argonne_deuteriumPot

    type(particle) :: p
    real :: pot,r
    integer :: i
    logical :: success

    pot=0.
    success=.false.

    if (.not.p%perturbative) then

       if (deuterium_pertOrigin_flag.eq.99) then
          ! 2 Body collisions
          r=absVec(deuterium_pointerList(deuteriumPL_ensemble)%part1%position(1:3)&
               & -deuterium_pointerList(deuteriumPL_ensemble)%part2%position(1:3))
          if (r.lt.0.00000001) then
             write(*,*) 'Real Event'
             write(*,'(A,3F10.4)') 'Partner ',deuterium_pointerList(deuteriumPL_ensemble)%part1%position(1:3)
             write(*,'(A,3F10.4)') 'Original',deuterium_pointerList(deuteriumPL_ensemble)%part2%position(1:3)
             write(*,*)
             write(*,*) 'Error in REAL_search (2). STOP',deuterium_pertOrigin
             stop
          end if
          pot=argonne_deuteriumPot(r)
          success=.true.
          return
       end if

       if (.not.success) then
          !       write(*,*) 'Direct search'
          directSearch: do i=lbound(deuterium_pointerList,dim=1),ubound(deuterium_pointerList,dim=1)
             if (equalParticles(deuterium_pointerList(i)%part1,p)) then
                r=absVec(deuterium_pointerList(i)%part2%position(1:3)&
                     & -p%position(1:3))
                pot=argonne_deuteriumPot(r)
                success=.true.
                exit directSearch
             else if (equalParticles(deuterium_pointerList(i)%part2,p)) then
                r=absVec(deuterium_pointerList(i)%part1%position(1:3)&
                     & -p%position(1:3))
                pot=argonne_deuteriumPot(r)
                success=.true.
                exit directSearch
             end if
          end do directSearch
       else
          return
       end if

       if (.not.success) then
          if (associated(deuterium_pertOrigin)) then
             realSearch: do i=lbound(deuterium_pointerList,dim=1),ubound(deuterium_pointerList,dim=1)
                if (equalParticles(deuterium_pointerList(i)%part1,deuterium_pertOrigin)) then
                   r=absVec(deuterium_pointerList(i)%part2%position(1:3)&
                        & -p%position(1:3))
                   if (r.lt.0.00000001) then
                      write(*,*) 'Perturbative Event'
                      write(*,'(A,3F10.4)') 'Partner ',deuterium_pointerList(i)%part2%position(1:3)
                      write(*,'(A,3F10.4)') 'Original',deuterium_pointerList(i)%part1%position(1:3)
                      write(*,'(A,3F10.4)') 'Particle',p%position
                      write(*,*)

                      write(*,*) 'Error in pert_search (1). STOP',deuterium_pertOrigin_flag
                      stop
                   end if
                   pot=argonne_deuteriumPot(r)
                   success=.true.
                   exit realSearch
                else if (equalParticles(deuterium_pointerList(i)%part2,deuterium_pertOrigin)) then
                   r=absVec(deuterium_pointerList(i)%part1%position(1:3)&
                        & -p%position(1:3))
                   if (r.lt.0.00000001) then
                      write(*,*) 'Perturbative Event'
                      write(*,'(A,3F10.4)') 'Partner ',deuterium_pointerList(i)%part1%position(1:3)
                      write(*,'(A,3F10.4)') 'Original',deuterium_pointerList(i)%part2%position(1:3)
                      write(*,'(A,3F10.4)') 'Particle',p%position
                      write(*,*)
                      write(*,*) 'Error in pert_search (2). STOP',deuterium_pertOrigin_flag
                      stop
                   end if
                   pot=argonne_deuteriumPot(r)
                   success=.true.
                   exit realSearch
                end if
             end do realSearch
          end if
       end if

       if (.not.success) then
          write(*,*) 'Could not find particle (REAL)',deuterium_pertOrigin_flag
          pot=0.
          success=.false.
          return
       end if

    else

       if (associated(deuterium_pertOrigin)) then
          pertSearch: do i=lbound(deuterium_pointerList,dim=1),ubound(deuterium_pointerList,dim=1)
             if (equalParticles(deuterium_pointerList(i)%part1,deuterium_pertOrigin)) then
                r=absVec(deuterium_pointerList(i)%part2%position(1:3)&
                     & -p%position(1:3))
                if (r.lt.0.00000001) then
                   write(*,*) 'Perturbative Event'
                   write(*,'(A,3F10.4)') 'Partner ',deuterium_pointerList(i)%part2%position(1:3)
                   write(*,'(A,3F10.4)') 'Original',deuterium_pointerList(i)%part1%position(1:3)
                   write(*,'(A,3F10.4)') 'Particle',p%position
                   write(*,*)

                   write(*,*) 'Error in pert_search (1). STOP',deuterium_pertOrigin_flag
                   stop
                end if
                pot=argonne_deuteriumPot(r)
                success=.true.
                exit pertSearch
             else if (equalParticles(deuterium_pointerList(i)%part2,deuterium_pertOrigin)) then
                r=absVec(deuterium_pointerList(i)%part1%position(1:3)&
                     & -p%position(1:3))
                if (r.lt.0.00000001) then
                   write(*,*) 'Perturbative Event'
                   write(*,'(A,3F10.4)') 'Partner ',deuterium_pointerList(i)%part1%position(1:3)
                   write(*,'(A,3F10.4)') 'Original',deuterium_pointerList(i)%part2%position(1:3)
                   write(*,'(A,3F10.4)') 'Particle',p%position
                   write(*,*)
                   write(*,*) 'Error in pert_search (2). STOP',deuterium_pertOrigin_flag
                   stop
                end if
                pot=argonne_deuteriumPot(r)
                success=.true.
                exit pertSearch
             end if
          end do pertSearch
       else
          pot=0.
       end if

    end if

  contains

      !************************************************************************
      !****s* get_DeuteriumPotential/equalParticles
      ! NAME
      ! logical function equalParticles(a,b)
      !
      ! PURPOSE
      ! Checks whether two particles are equal. Works only if all particles have different
      ! %number entry (i.e. not in fullEnsemble mode)!!!
      !
      ! INPUTS
      ! type(particle) :: a,b
      !************************************************************************
      logical function equalParticles(a,b)
        use output, only: WriteParticle_debug

        type(particle) :: a,b
        real,parameter :: eps_pos=0.5
        real,parameter :: eps_mom=0.5
        logical, parameter :: debug=.false.

        equalParticles=.false.

        if (a%number.eq.b%number) then
           if (debug) then
              write(*,*)
              equalParticles=.true.
              write(*,*) a%position-b%position
              write(*,*) a%momentum-b%momentum
              write(*,*) a%number, b%number
              return
           end if
           if (AbsVec(a%position-b%position).gt.eps_pos) then
              write(*,*) 'WARNING: PARTICLES NOT EQUAL: Positions'
              write(*,*) 'diff=',a%position-b%position
              call WriteParticle_debug(a)
              call WriteParticle_debug(b)
           end if
           if (AbsVec(a%momentum(1:3)-b%momentum(1:3)).gt.eps_mom) then
              write(*,*) 'WARNING: PARTICLES NOT EQUAL: Momenta'
              write(*,*) 'diff=',a%momentum(1:3)-b%momentum(1:3)
              call WriteParticle_debug(a)
              call WriteParticle_debug(b)
           end if
           equalParticles=.true.
           return
        end if

      end function equalParticles

  end function get_DeuteriumPotential


end module deuterium_PL
