
program test
implicit none

call init_dataBase
call testMaster


end program


  subroutine testMaster
    use master_2Body
    use particleDefinition 
    use idTable
    use particleProperties
    use inputGeneral , only : delta_t
    use timing
    use preEventDefinition
    implicit none

    type(particle),dimension(1:2)    :: pair                   ! incoming pair of particles
    integer               :: numEnsembles ! number of ensembles
    type(particle), dimension(1:30)         :: finalState         ! produced final state
    logical        :: collisionFlag     ! true if collision takes plac
    integer :: i,j
    integer :: id1,id2, q1,q2
    logical :: HiEnergyFlag
    integer :: HiEnergyType
    real :: mom
    real , dimension(0:3) :: totMom_Vorher,totMom_VorherSum,totMom_Nachher,totMom_NachherSum
    
    write(*,*) 'Testing collide_2body'

    call readinput


    pair%Id=(/id1,id2/)
    pair%charge=(/q1,q2/)
    pair(1)%event=(/1,1/)
    pair(2)%event=(/2,2/)


    If(isMeson(id1)) then
       pair(1)%mass=meson(id1)%mass
    else if (isBaryon(id1)) then
       pair(1)%mass=baryon(id1)%mass
    else
       write(*,*) 'id1 is no particle:', id1
    end if


    If(isMeson(id2)) then
       pair(2)%mass=meson(id2)%mass
    else if (isBaryon(id2)) then
       pair(2)%mass=baryon(id2)%mass
    else
       write(*,*) 'id2 is no particle:', id2
    end if

    pair(1)%position=(/0.001,0.,0./)
    pair(2)%position=(/0.,0.,0./)

    pair(1)%momentum(1:3)=(/0.,0.,0./)
    pair(1)%momentum(0)=sqrt(pair(1)%mass**2+Dot_Product(pair(1)%momentum(1:3),pair(1)%momentum(1:3)))

    pair(2)%momentum(1:3)=(/0.,mom,0./)
    pair(2)%momentum(0)=sqrt(pair(2)%mass**2+Dot_Product(pair(2)%momentum(1:3),pair(2)%momentum(1:3)))

    pair(1)%velocity=pair(1)%momentum(1:3)/pair(1)%momentum(0)
    pair(2)%velocity=pair(2)%momentum(1:3)/pair(2)%momentum(0)

    numEnsembles=100

    totMom_NachherSum=0.
    totMom_Vorher=pair(1)%momentum+pair(2)%momentum

    Do i=1,2000
       Print *
       Print *, 'momenta'
       Print *, pair(1)%momentum
       Print *, pair(2)%momentum
       Print *, 'Wurzel(s)=',sqrts(pair(1),pair(2))
       Print *, 'Total momentum=',pair(1)%momentum+pair(2)%momentum
       print * ,'***********************'
          call setToDefault(finalState)
       call collide_2body(pair,finalState,numEnsembles,collisionFlag,0.57,HiEnergyFlag,HiEnergyType)
       If (collisionFlag) then
          totMom_VorherSum=totMom_VorherSum+totMom_Vorher
           Print *, 'Finalstate:       ############################',collisionFlag
          Print *, finalState(1:4)%ID, collisionFlag, HiEnergyFlag,HiEnergyType
          Print *, 'Charges' , finalState(1:4)%charge, finalState(1:4)%antiParticle
          totMom_Nachher=0.
          Do j=lBound(finalState,dim=1),uBound(finalState,dim=1)
             If(finalState(j)%ID.ne.0) then
                totMom_Nachher=totMom_Nachher+finalState(j)%momentum
                write(9,*) j, finalstate(j)%ID, finalState(j)%momentum
             end if
          end do
          Print *, ' Total momentum=',totMom_Nachher
          totMom_NachherSum=totMom_NachherSum+totMom_Nachher
       else
          Print *, 'No collision',collisionFlag
       end if

    End do
 
    Print*, "total Momentum before" ,totMom_VorherSum
    Print*, "total Momentum after" ,totMom_NachherSum


    call timeMeasurement


contains

  subroutine readInput
    !****s* master_2Body/readInput
    ! NAME
    ! subroutine readInput
    ! FUNCTION
    ! Reads input in jobcard out of namelist "master_2Body"
    !***

    use output
    
    NAMELIST /test/ id1,id2,q1,q2,mom

    call Write_ReadingInput('test',0)
    rewind(5)
    read(5,nml=test)
    write(*,*) '  Id of first particle',id1

    write(*,*) '  Id of second particle'    ,id2
    write(*,*) '  Charge of first particle', q1
    write(*,*) '  Charge of second particle', q2
    write(*,*) '  Mom of second particle', mom

    call Write_ReadingInput('test',1)

    
  end subroutine readInput

  end subroutine testMaster

