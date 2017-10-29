program test
  use inputGeneral



  call readinputGeneral
  call init_Database
  call tester

end program test




subroutine tester
  use particleDefinition
  use ParticleProperties
  use energyCalc
  use offShellPotential
  use minkowski
  use propagation
  use mediumDefinition
  use mediumModule, only: mediumAt
  use baryonwidthmedium, only :WidthBaryonMedium
  use mesonwidthmedium, only: WidthMesonMedium
  use version

  implicit none

  type(particle),dimension(1:1,1:1) :: finalstate
  integer :: i,j,n
  type(medium) :: mediumAtPosition
  real :: dummy,width,energy_previous
  logical :: success
  real    :: error_perTimeStep

  integer,save :: charge=1
  real, save :: delta_T=0.1
  integer, save :: particle_id=2
  real, save :: mass=1.12
  real, save :: mom=0.5
  real, save :: pos=0.

  logical :: outOfBounds

  call PrintVersion

  write(200,*)
  write(101,*)
  write(102,*)
  write(4711,*)
  write(4712,*)
  write(4713,*)
  write(4714,*)


  call readInput

  call setToDefault(finalState(1,1))



  finalState%antiparticle=.false.
  finalState%perturbative=.true.
  finalState%productionTime=0.
  finalState%lastCollisionTime=0.
  finalState%formationTime=0.
  finalState%scaleCS=1.
  finalState%in_Formation=.false.
  finalstate%history=0
  finalState(1,1)%position(1)=0.
  finalState(1,1)%position(2)=0.
  finalState(1,1)%position(3)=pos
  finalstate%id=particle_id
  finalstate%charge=charge
  finalstate%mass=mass
  finalstate(1,1)%momentum(1)=0.
  finalstate(1,1)%momentum(2)=0.
  finalstate(1,1)%momentum(3)=mom

  call energyDetermination(finalstate(1,1))
  finalstate%offshellparameter= &
       & getOffShellParameter(finalstate(1,1)%ID,finalstate(1,1)%Mass,finalstate(1,1)%momentum,finalstate(1,1)%position,success)
  mediumAtPosition=mediumAt(finalstate(1,1)%position)

  if (isBaryon(particle_ID)) then
    width=WidthBaryonMedium(finalstate(1,1)%ID,finalstate(1,1)%mass,finalstate(1,1)%momentum,mediumATposition)
  else if (isMeson(particle_ID)) then
    width=WidthMesonMedium(finalstate(1,1)%ID,finalstate(1,1)%mass,finalstate(1,1)%momentum,mediumATposition)
  end if

  write(*,*) width,abs4(finalstate(1,1)%momentum)
  write(*,*) finalstate(1,1)%momentum
  write(*,*) finalstate(1,1)%offshellparameter,mediumAtPosition%density

  if(.not.success) then
     ! REJECT event
     stop 'offshellparameter'
  end if
  if(.not.treatParticleOffShell(finalstate(1,1)%ID,finalstate(1,1)%OffShellParameter)) stop 'not offshell'

  write(200,*)   finalState(1,1)%momentum(0), HamiltonFunc_offshell(finalstate(1,1),outOfBounds)

  write(*,*) 'starting propa loop',finalstate(1,1)%offshellparameter, finalstate(1,1)%momentum(0)
  energy_previous=finalstate(1,1)%momentum(0)
  error_perTimeStep=energy_previous
  propaloop: do j=1,100000000
     call propagate(finalstate,finalstate,delta_T,(/0.,0.,0./),(/0.,0.,0./),.true.)
     if(finalstate(1,1)%id.eq.0) stop 'id null'
     call energyDetermination(finalstate(1,1))
     write(200,*)   finalState(1,1)%momentum(0), HamiltonFunc_offshell(finalstate(1,1),outOfBounds)

     !if(j*delta_T.lt.3.0+delta_T/2..and.j*delta_T.gt.3.0-delta_T/2.) then
     !   call rloop(finalstate(1,1))
     !   call ploop(finalstate(1,1))
     !   call eloop(finalstate(1,1))
     !end if

     mediumAtPosition=mediumAt(finalstate(1,1)%position)
     if (isBaryon(particle_ID)) then
       width=WidthBaryonMedium(finalstate(1,1)%ID,finalstate(1,1)%mass,finalstate(1,1)%momentum,mediumATposition)
     else if (isMeson(particle_ID)) then
       width=WidthMesonMedium(finalstate(1,1)%ID,finalstate(1,1)%mass,finalstate(1,1)%momentum,mediumATposition)
     end if
     !write(*,'(10G20.6)') j*delta_T,sqrt(dot_product(finalstate(1,1)%momentum(1:3),finalstate(1,1)%momentum(1:3))), &
     !     & sqrt(dot_product(finalstate(1,1)%position,finalstate(1,1)%position)), &
     !     & mediumAtPosition%density, Width, &
     !     & finalstate(1,1)%momentum(0), finalstate(1,1)%mass
     !write(*,*)
     write(101,'(10G20.4)') finalstate(1,1)%position,finalstate(1,1)%momentum

     write(4714,'(10G20.6)') j*delta_T,sqrt(dot_product(finalstate(1,1)%momentum(1:3),finalstate(1,1)%momentum(1:3))), &
          & sqrt(dot_product(finalstate(1,1)%position,finalstate(1,1)%position)), &
          & mediumAtPosition%density, width, &
          & finalstate(1,1)%momentum(0), finalstate(1,1)%mass

     !check energy conservation:
     if(abs(energy_previous-finalstate(1,1)%momentum(0)).gt.0.001) then
        write(*,*) 'energy cons failed: ',finalstate(1,1)%momentum(0)-energy_previous
        call rloop(finalstate(1,1))
        call ploop(finalstate(1,1))
        call eloop(finalstate(1,1))
        write(*,*) finalState(1,1)%position
        write(*,*) 'Time=', delta_T*float(j),j
  !      stop
     end if

     energy_previous=finalstate(1,1)%momentum(0)

     write(102,*) error_perTimeStep- finalstate(1,1)%momentum(0)
     error_perTimeStep=finalstate(1,1)%momentum(0)

     if(absPos(finalState(1,1)).gt.7) exit propaLoop
     if(j*delta_T.gt.50) exit propaLoop
  end do propaloop


  write(*,*)'end of hamilton test'


contains

  subroutine eloop(finalstate)
    implicit none
    type(particle),dimension(1:1,1:1), intent(in) :: finalstate
    type(particle),dimension(1:1,1:1) :: dummyfinalstate

    real :: E_in
    real :: E
    call setToDefault(dummyfinalState(1,1))
    dummyfinalstate=finalstate

    E_in=finalstate(1,1)%momentum(0)

    write(*,*) 'starting E loop',dummyfinalstate(1,1)%offshellparameter,E_in
    do i=-5,5
       write(*,*) 'loop',i
       E=E_in+i*0.01
       dummyfinalstate(1,1)%momentum(0)=E
       mediumAtPosition=mediumAt(dummyfinalstate(1,1)%position)
       if (isBaryon(particle_ID)) then
         width=WidthBaryonMedium(dummyfinalstate(1,1)%ID,dummyfinalstate(1,1)%mass,dummyfinalstate(1,1)%momentum,mediumATposition)
       else if (isMeson(particle_ID)) then
         width=WidthMesonMedium(dummyfinalstate(1,1)%ID,dummyfinalstate(1,1)%mass,dummyfinalstate(1,1)%momentum,mediumATposition)
       end if
       write(4713,'(16E20.5)') E,abs4(dummyfinalstate(1,1)%momentum),HamiltonFunc_offshell(dummyfinalstate(1,1),outOfBounds),  &
            & getOffShellMass(dummyfinalstate(1,1)%ID,dummyfinalstate(1,1)%offshellparameter,dummyfinalstate(1,1)%momentum,&
            & dummyfinalstate(1,1)%mass,dummyfinalstate(1,1)%position,outOfBounds),  &
            & Width,E_in
    end do

  end  subroutine eloop


  subroutine rloop(finalstate)
    implicit none
    type(particle),dimension(1:1,1:1), intent(in) :: finalstate
    type(particle),dimension(1:1,1:1) :: dummyfinalstate

    real :: r_in
    real :: r
    call setToDefault(dummyfinalState(1,1))
    dummyfinalstate=finalstate

    r_in=finalstate(1,1)%position(3)

    write(*,*) 'starting r loop',dummyfinalstate(1,1)%offshellparameter,r_in
    do i=-5,5
       write(*,*) 'loop',i
       r=r_in+i*0.4
       dummyfinalstate(1,1)%position(3)=r
       mediumAtPosition=mediumAt(dummyfinalstate(1,1)%position)
       if (isBaryon(particle_ID)) then
         width=WidthBaryonMedium(dummyfinalstate(1,1)%ID,dummyfinalstate(1,1)%mass,dummyfinalstate(1,1)%momentum,mediumATposition)
       else if (isMeson(particle_ID)) then
         width=WidthMesonMedium(dummyfinalstate(1,1)%ID,dummyfinalstate(1,1)%mass,dummyfinalstate(1,1)%momentum,mediumATposition)
       end if
       write(4712,'(16E20.5)') r,abs4(dummyfinalstate(1,1)%momentum),HamiltonFunc_offshell(dummyfinalstate(1,1),outOfBounds), &
            & getOffShellMass(dummyfinalstate(1,1)%ID,dummyfinalstate(1,1)%offshellparameter,dummyfinalstate(1,1)%momentum,&
            & dummyfinalstate(1,1)%mass,dummyfinalstate(1,1)%position,outOfBounds),  &
            & width,r_in
    end do

  end subroutine rloop


  subroutine ploop(finalstate)
    implicit none
    type(particle),dimension(1:1,1:1), intent(in) :: finalstate
    type(particle),dimension(1:1,1:1) :: dummyfinalstate

    real :: p_in
    real :: p

    call setToDefault(dummyfinalState(1,1))
    dummyfinalstate=finalstate

    p_in=finalstate(1,1)%momentum(3)

    write(*,*) 'starting p loop',dummyfinalstate(1,1)%offshellparameter,p_in
    do i=-5,5
       write(*,*) 'loop',i
       p=p_in+i*0.01
       dummyfinalstate(1,1)%momentum(3)=p

       mediumAtPosition=mediumAt(dummyfinalstate(1,1)%position)
       if (isBaryon(particle_ID)) then
         width=WidthBaryonMedium(dummyfinalstate(1,1)%ID,dummyfinalstate(1,1)%mass,dummyfinalstate(1,1)%momentum,mediumATposition)
       else if (isMeson(particle_ID)) then
         width=WidthMesonMedium(dummyfinalstate(1,1)%ID,dummyfinalstate(1,1)%mass,dummyfinalstate(1,1)%momentum,mediumATposition)
       end if
       write(4711,'(16E20.5)') p,abs4(dummyfinalstate(1,1)%momentum),HamiltonFunc_offshell(dummyfinalstate(1,1),outOfBounds), &
            & getOffShellMass(dummyfinalstate(1,1)%ID,dummyfinalstate(1,1)%offshellparameter,dummyfinalstate(1,1)%momentum,&
            & dummyfinalstate(1,1)%mass,dummyfinalstate(1,1)%position,outOfBounds),  &
            & Width,p_in
    end do

  end subroutine ploop

  subroutine readInput
    use output

    implicit none
    integer :: ios
    NAMELIST /tester/ delta_T,particle_Id,mass,mom,charge,pos

    call Write_ReadingInput('tester',0)
    rewind(5)
    read(5,nml=tester,IOSTAT=IOS)
    call Write_ReadingInput('tester',0,IOS)

    write(*,*) 'delta_T', delta_T
    write(*,*) 'particle_Id',particle_ID
    write(*,*) 'mass', mass
    write(*,*) 'z-momentum', mom
    write(*,*) 'charge',charge
    write(*,*) 'z-pos', pos

    call Write_ReadingInput('offShellPotential',1)
  end subroutine readInput



end subroutine tester
