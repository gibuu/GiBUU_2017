program test_massass_omegaN

    use inputGeneral
    use particleDefinition
    use particleProperties, only: baryon
    use master_2Body,only: generateFinalState
    use IDTable, only: nucleon, omegaMeson

    implicit none

    integer :: i
    real    :: stringFactor=1.
    type(particle), dimension(1:2)  :: pair
    type(particle), dimension(1:10) :: finalState
    real    :: time=999.
    logical :: collisionFlag,pauliIsUsedforXsection
    logical :: HiEnergyFlag  ! .true. if fritiof was used
    integer :: HiEnergyType  ! 0:LowEnergy, 1:Fritiof, 2:Pythia
    real    :: sigmaTot,srts

    call init_Database
    call readinputGeneral

    call setToDefault(pair)

    pair(1)%ID=omegaMeson
    pair(1)%mass=0.15
    pair(1)%momentum(1:2)= 0.
    pair(1)%momentum(3)=0.02
    pair(1)%momentum(0) = freeEnergy(pair(1))
    pair(1)%velocity=pair(1)%momentum(1:3)/pair(1)%momentum(0)
    !pair(1)%position=999.   ! outside nucleus
    pair(1)%charge=0

    pair(2)%ID=nucleon
    pair(2)%mass=baryon(nucleon)%mass
    pair(2)%position=pair(1)%position
    pair(2)%charge=1
    pair(2)%momentum(1:3) = 0.
    pair(2)%momentum(0) = freeEnergy(pair(2))
    pair(2)%velocity    = pair(2)%momentum(1:3)/pair(2)%momentum(0)

    srts=sqrtS(pair) 
    write(*,*) 'Srts: ', srts

    do i=1,1000

       call setToDefault(finalState)
       finalState%ID=0
       call generateFinalState(pair,finalState,stringFactor,numEnsembles,collisionFlag,time,HiEnergyFlag,HiEnergyType, &
             sigTot_out=sigmaTot,pauliIncluded_out=pauliIsUsedforXsection)

       if (.not. collisionFlag) then
         write(*,*) 'In:',pair%ID
         write(*,*) 'OUT:',finalState(1:3)%ID
         write(*,*) 'sigma:',sigmaTot
         write(*,*) 'flags:',collisionFlag,HiEnergyFlag,HiEnergyType
         write(*,'(4G18.3)') absmom(pair(1)), absMom(pair(2))
         write(*,'(4G18.3)') pair(1)%mass, pair(2)%mass
         write(*,'(4G18.3)') finalstate(1)%mass, finalstate(2)%mass
         write(*,'(4G18.3)') AbsMom(finalState(1)),AbsMom(finalState(2))
         write(*,*)
       end if

       write(*,*) "i = ",i

    end do


end program
