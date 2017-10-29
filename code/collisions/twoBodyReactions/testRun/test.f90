program test
use inputGeneral
use version

implicit none

call PrintVersion

call readInputGeneral           ! in module inputGeneral
call init_dataBase

!call testMaster
call crossMaster

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
    integer :: i
    integer :: id1,id2,q1,q2
    logical :: HiEnergyFlag
    integer :: HiEnergyType


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

    pair(2)%momentum(1:3)=(/0.,1.2,0./)
    pair(2)%momentum(0)=sqrt(pair(2)%mass**2+Dot_Product(pair(2)%momentum(1:3),pair(2)%momentum(1:3)))

    pair(1)%velocity=pair(1)%momentum(1:3)/pair(1)%momentum(0)
    pair(2)%velocity=pair(2)%momentum(1:3)/pair(2)%momentum(0)

    numEnsembles=100

    Do i=1,20
       pair(2)%momentum(1:3)=(/0.,i*0.1,0./)
       pair(2)%momentum(0)=sqrt(pair(2)%mass**2+Dot_Product(pair(2)%momentum(1:3),pair(2)%momentum(1:3)))
       Print *
       Print *, 'momenta'
       Print *, pair(1)%momentum
       Print *, pair(2)%momentum
       Print *, 'Wurzel(s)=',sqrts(pair(1),pair(2))
       Print *, 'Total momentum=',pair(1)%momentum+pair(2)%momentum
       print * ,'***********************'
       call collide_2body(pair,finalState,numEnsembles,collisionFlag,0.57,HiEnergyFlag,HiEnergyType)
       If (collisionFlag) then
          Print *, 'Finalstate:       ############################'
          Print *, finalState(1:4)%ID, collisionFlag, sqrtS(finalState(1),finalState(2),finalState(3))
          Print *, '             Charges' , finalState(1:4)%charge, finalState(1:4)%antiParticle
          Print *, ' Total momentum=',finalState(1)%momentum+finalState(2)%momentum+finalState(3)%momentum
       else
          Print *, 'No collision'
       end if
    End do


    write(*,*)
    write(*,*) 'checking speed'
    pair(2)%momentum(1:3)=(/0.,1.5,0./)
    pair(2)%momentum(0)=sqrt(pair(2)%mass**2+Dot_Product(pair(2)%momentum(1:3),pair(2)%momentum(1:3)))
    write(*,*) 'Wurzel(s)=', sqrtS(pair(1),pair(2))
    write(*,*) pair%ID
    write(*,*)


    call timeMeasurement(.true.)
    Do i=1,200
       call collide_2body(pair,finalState,numEnsembles,collisionFlag,0.57,HiEnergyFlag,HiEnergyType)
    End do

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

    NAMELIST /test/ id1,id2,q1,q2

    call Write_ReadingInput('test',0)
    rewind(5)
    read(5,nml=test)
    write(*,*) '  Id of first particle',id1

    write(*,*) '  Id of second particle'    ,id2
    write(*,*) '  Charge of first particle', q1
   write(*,*) '  Charge of second particle', q2

    call Write_ReadingInput('test',1)


  end subroutine readInput

  end subroutine testMaster








!!$!****************************************************************************************************



  subroutine crossMaster
    use master_2Body
    use particleDefinition
    use idTable
    use mediumDefinition
    use particleProperties
    use inputGeneral , only : delta_t
    use timing
    use preEventDefinition
    use energyCalc
    implicit none

    type(particle),dimension(1:2)    :: pair                   ! incoming pair of particles
    integer               :: numEnsembles ! number of ensembles
    type(preEvent), dimension(1:3)         :: finalState         ! produced final state
    logical        :: HiFlag
    integer :: i,j
    real, dimension(0:3) :: momentumLRF
    real, dimension(1:3) :: betaToCF
    TYPE(medium) :: vacuum
    real :: srts, sigmaTot, sigmaElast
    real :: sigma_F(2), sigma_R(2)

    integer :: q1,q2,id1,id2
    logical :: anti1,anti2
    real :: dens=0.168
    real :: Sum1,Sum2, rHiEnergy

    call readinput

    vacuum%useMedium=.true.
    vacuum%temperature=0.
    vacuum%densityProton=dens/2.
    vacuum%densityNeutron=dens/2.
    vacuum%density=dens

    write(*,*) 'Testing xSectionMaster'
    pair%Id=(/id1,id2/)
    pair%charge=(/q1,q2/)
    pair%antiparticle=(/anti1,anti2/)
    pair(1)%event=(/1,1/)
    pair(1)%perturbative=.true.
    pair(2)%perturbative=.false.
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


    pair(1)%position=(/1.,0.,0./)
    pair(2)%position=(/0.,0.,0./)

    pair(1)%momentum(1:3)=(/0.,0.,0./)

    !pair(1)%momentum(0)=sqrt(pair(1)%mass**2+Dot_Product(pair(1)%momentum(1:3),pair(1)%momentum(1:3)))
    call energyDetermination(pair(1),(/0.,0.,0./))
    call energyDetermination(pair(2),(/0.,0.,0./))

    pair(1)%velocity=pair(1)%momentum(1:3)/pair(1)%momentum(0)


    numEnsembles=100

    Print *, 'positions'
    Print *, pair(1)%position
    Print *, pair(2)%position
    Print *, 'momenta'
    Print *, pair(1)%momentum
    Print *, pair(2)%momentum
    Print *, 'Wurzel(s)=',sqrts(pair(1),pair(2))
    Print *, 'Total momentum=',pair(1)%momentum+pair(2)%momentum
    print * ,'***********************'

    Do i=1,2000
       pair(2)%momentum(1:3)=(/0.,float(i)*0.05,0./)
       !pair(2)%momentum(0)=sqrt(pair(2)%mass**2+Dot_Product(pair(2)%momentum(1:3),pair(2)%momentum(1:3)))

       betaToCF=(/0.,0.,0./)
       call energyDetermination(pair(2),betaToCF)
       pair(2)%velocity=pair(2)%momentum(1:3)/pair(2)%momentum(0)

       srts=sqrts(pair(1),pair(2))

       if (srts.gt.10.5) exit

       write(*,*) '--- srts = ',srts

       momentumLRF=pair(2)%momentum+pair(1)%momentum


       rHiEnergy = HiEnergyContib(srts,pair)

       sigma_R=0.0
       sigma_F=0.0
       if (rHiEnergy.lt.1.0) then
          call XsectionMaster(srts,pair,vacuum,momentumLRF,finalState,sigma_R(1),sigma_R(2),HiFlag,.true.,ForceHiEnergy=.false.)
       endif
       if (rHiEnergy.gt.0.0) then
          call XsectionMaster(srts,pair,vacuum,momentumLRF,finalState,sigma_F(1),sigma_F(2),HiFlag,.true.,ForceHiEnergy=.true.)
       endif

       write(23,'(10F12.4)') absmom(pair(2)),srts,&
            (1.0-rHiEnergy)*sigma_R+rHiEnergy*sigma_F, sigma_R, sigma_F, rHiEnergy,kineticEnergy(pair(2))


    End do


contains

  subroutine readInput
    !****s* master_2Body/readInput
    ! NAME
    ! subroutine readInput
    ! FUNCTION
    ! Reads input in jobcard out of namelist "master_2Body"
    !***

    use output
    use particleProperties
    character(40) :: name1,name2

    NAMELIST /test/ id1,id2,q1,q2,anti1,anti2,dens

    call Write_ReadingInput('test',0)
    rewind(5)
    read(5,nml=test)
    write(*,*) '  Id of first particle      : ',id1
    write(*,*) '  Id of second particle     : ',id2
    write(*,*) '  Charge of first particle  : ', q1
    write(*,*) '  Charge of second particle : ', q2
    write(*,*) '  antiparticle 1            : ', anti1
    write(*,*) '  antiparticle 2            : ', anti2
    write(*,*) '  Density=', dens

    call Write_ReadingInput('test',1)

    if(ismeson(id1)) then
         if(ismeson(id2)) then
            !    write(100,*) '<h3>', meson(id1)%name, meson(id2)%name, ' cross section','</h3><BR>'
            name1=meson(id1)%name
            name2=meson(id2)%name
         else
            !    write(100,*) '<h3>', meson(id1)%name, baryon(id2)%name, ' cross section','</h3><BR>'
            name1=meson(id1)%name
            name2=baryon(id2)%name
         end if
    else
         if(ismeson(id2)) then
            !    write(100,*) '<h3>', baryon(id1)%name, meson(id2)%name, ' cross section','</h3><BR>'
            name1=baryon(id1)%name
            name2=meson(id2)%name
         else
            !    write(100,*) '<h3>', baryon(id1)%name, baryon(id2)%name, ' cross section','</h3><BR>'
            name1=baryon(id1)%name
            name2=baryon(id2)%name
         end if
    end if
    write(100,*) '<table border="1" cellspacing="5" cellpadding="10">'
    write(100,*) '<tr>'
    write(100,*) '<th>  Scattering particles  </th> <th>',  name1,'  </th> <th> ', name2,'  </th> '
    write(100,*) '</tr>'
    write(100,*) '<tr>'
    write(100,*) '<td>  Charge  </td> <td>',  q1,'  </td> <td> ', q2,'  </td> '
    write(100,*) '</tr>'
    write(100,*) '<tr>'
    write(100,*) '<td>  Antiparticle  </td> <td>',  anti1,'  </td> <td> ', anti2,'  </td> '
    write(100,*) '</tr>'
    write(100,*) '</table>'

  end subroutine readInput



end subroutine crossMaster
