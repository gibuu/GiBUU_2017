module KodamaCluster

  implicit none
  PRIVATE

  public :: TryKodamaClusters

contains

  subroutine TryKodamaClusters (realParticles, time)

    use particleDefinition
    use IdTable, only: isMeson, isBaryon
    use collisionCriteria, only : kodama_position, kodama_time
    use master_2Body, only: check_veryRough, checkKodama_rough
!     use output, only: WriteParticle
    use collisionNumbering, only: check_justCollided
    use constants, only: twopi, hbarc
    use inputGeneral, only : delta_T

    type(particle),intent(inOUT),dimension(:,:) :: realParticles
    real :: time

    integer iEns,iPart1,iPart2, nEns
    type(particle),dimension(2)    :: pair 
    integer :: scenario

    integer :: nGroup
    integer, parameter :: nGroupMax = 200, sizeGroupMax = 10000
    integer :: sizeGroup(nGroupMax), Group(nGroupMax,sizeGroupMax)

    integer :: iGr,iP, iGr1,iGr2, iP1,iP2, iP_Min(2)

    real :: RelMom(3), AbsRelMom, AbsRelMom_Min, lambda_Min
    real :: RelPos(3)!, AbsRelPos

    integer, parameter:: sizeClusterMax = sizeGroupMax
    integer :: sizeCluster, Cluster(sizeClusterMax) ! not the iP but the number of the iP !

    integer,parameter :: baryonBaryon_scenario=1
    integer,parameter :: baryonMeson_scenario=2
    integer,parameter :: MesonMeson_scenario=3

    integer, dimension(50) :: AnzColl

    real :: stringFactor

    AnzColl = 0


    nEns = Size(realParticles,dim=1)
    EnsLoop: do iEns = 1,nEns
       nGroup = 0
       sizeGroup = 0


       ! ===== 1) Building up the List of Groups =====


       Part1Loop: do iPart1=1,Size(realParticles,dim=2)
          if (realparticles(iEns,iPart1)%ID <= 0) cycle Part1Loop
          pair(1) = realparticles(iEns,iPart1)
          Part2Loop: do iPart2=iPart1+1,Size(realParticles,dim=2)
             if (realparticles(iEns,iPart2)%ID <= 0)            cycle Part2Loop
             pair(2) = realparticles(iEns,iPart2)

             ! --- Check for collision [Start] ---

             If (check_justCollided(pair(1),pair(2)))           cycle Part2Loop

             ! --- Set collision scenario ---

             If (isBaryon(pair(1)%ID).and.isBaryon(pair(2)%ID)) then    ! baryon-baryon-scattering
                scenario=baryonBaryon_scenario
             else if (isMeson(pair(1)%ID).and.isMeson(pair(2)%ID)) then ! meson-meson-scattering
                scenario=MesonMeson_scenario
             else
                scenario=baryonMeson_scenario
             endif

             ! --- Set cross section scaling factor ---

             stringFactor=1.
             If(pair(1)%in_Formation) stringFactor=stringFactor*pair(1)%scaleCS
             If(pair(2)%in_Formation) stringFactor=stringFactor*pair(2)%scaleCS

             ! --- Check for collision [continue] ---

             If(.not.check_veryRough(pair,stringFactor,nEns,scenario))   cycle Part2Loop
             if(.not.kodama_time(pair,delta_T))                          cycle Part2Loop
             If(.not.checkKodama_rough(pair,stringFactor,nEns,scenario)) cycle Part2Loop

             !               write(*,*) 'Possible: ',iEns,iPart1,iPart2
             !               call WriteParticle(6,iEns,iPart1,pair(1))
             !               call WriteParticle(6,iEns,iPart2,pair(2))


             ! === Search in Groups ===

             iGr1 = 0
             iGr2 = 0
             grLoop1: do iGr=1,nGroup
                do iP=1,sizeGroup(iGr)
                   if (Group(iGr,iP) == iPart1) then
                      iGr1 = iGr
                      exit grLoop1
                   endif
                enddo
             end do grLoop1
             grLoop2: do iGr=1,nGroup
                do iP=1,sizeGroup(iGr)
                   if (Group(iGr,iP) == iPart2) then
                      iGr2 = iGr
                      exit grLoop2
                   endif
                enddo
             end do grLoop2

             ! === Insert in Groups ===

             if (iGr1==0 .and. iGr2==0) then ! create new group
                !                  write(*,*) 'not found: ',iPart1,iPart2

                if (nGroup == nGroupMax) then
                   write(*,*) 'nGroup == nGroupMax'
                   stop
                endif
                nGroup = nGroup +1
                sizeGroup(nGroup) = 2
                Group(nGroup,1) = iPart1
                Group(nGroup,2) = iPart2
             else if (iGr1==0) then          ! insert iPart1 in iGr2
                if (sizeGroup(iGr2)==sizeGroupMax) then
                   write(*,*) 'sizeGroup(iGr2)==sizeGroupMax'
                   stop
                endif
                !                  write(*,*) 'inserting:',iGr2,iPart1
                sizeGroup(iGr2) = sizeGroup(iGr2)+1
                Group(iGr2,sizeGroup(iGr2)) = iPart1
             else if (iGr2==0) then          ! insert iPart2 in iGr1
                if (sizeGroup(iGr1)==sizeGroupMax) then
                   write(*,*) 'sizeGroup(iGr1)==sizeGroupMax'
                   stop
                endif
                !                  write(*,*) 'inserting:',iGr1,iPart2
                sizeGroup(iGr1) = sizeGroup(iGr1)+1
                Group(iGr1,sizeGroup(iGr1)) = iPart2
             else
                if (iGr1==iGr2) then         ! nothing to do
                   !                     write(*,*) 'both found'
                else                         ! iPart1 and iPart2 in diff groups !!!
                   if (sizeGroup(iGr1)+sizeGroup(iGr2)>sizeGroupMax) then
                      write(*,*) 'sizeGroup(iGr1)+sizeGroup(iGr2)>sizeGroupMax'
                      stop
                   endif
                   !                     write(*,*) 'adding :', iGr2,'->',iGr1
                   do iP=1,sizeGroup(iGr2)
                      group(iGr1,sizeGroup(iGr1)+iP) = group(iGr2,iP)
                   end do
                   sizeGroup(iGr1) = sizeGroup(iGr1)+sizeGroup(iGr2)

                   if (iGr2 .ne. nGroup) then
                      !                        write(*,*) 'copying :', nGroup,'->',iGr2
                      sizeGroup(iGr2) = sizeGroup(nGroup)
                      group(iGr2,:) = group(nGroup,:)
                   endif
                   nGroup = nGroup-1

                endif
             endif

          end do Part2Loop
       end do Part1Loop

       if (nGroup == 0) cycle EnsLoop

       ! ===== 2) Printing the List of Groups =====

       !         do iGr=1,nGroup
       !            write(*,'(i3," [",i3,"] : ",500(i5))') iGr,sizeGroup(iGr), group(iGr,1:sizeGroup(iGr))
       !         end do

       ! ===== 3) Build Clusters for every Group =====


       grLoop: do iGr=1,nGroup

1003      continue

          ! --- 3a) find minimal (relative) de Broglie wavelength

          absRelMom_Min = 0.0
          do iP1=1,sizeGroup(iGr)
             do iP2=iP1+1,sizeGroup(iGr)

                ! test, whether this pair is allowed to scatter at all:
                ! (if we would have saved it above, this would do the same job)

                pair(1) = realparticles(iEns,group(iGr,iP1))
                pair(2) = realparticles(iEns,group(iGr,iP2))

                if (check_justCollided(pair(1),pair(2))) cycle

                If (isBaryon(pair(1)%ID).and.isBaryon(pair(2)%ID)) then ! baryon-baryon-scattering
                   scenario=baryonBaryon_scenario
                else if (isMeson(pair(1)%ID).and.isMeson(pair(2)%ID)) then ! meson-meson-scattering
                   scenario=MesonMeson_scenario
                else
                   scenario=baryonMeson_scenario
                endif

                stringFactor=1.
                If(pair(1)%in_Formation) stringFactor=stringFactor*pair(1)%scaleCS
                If(pair(2)%in_Formation) stringFactor=stringFactor*pair(2)%scaleCS

                If(.not.check_veryRough(pair,stringFactor,nEns,scenario))   cycle
                if(.not.kodama_time(pair,delta_T))                          cycle
                If(.not.checkKodama_rough(pair,stringFactor,nEns,scenario)) cycle

                RelMom(1:3) = pair(1)%momentum(1:3)-pair(2)%momentum(1:3)
                AbsRelMom = sqrt(DOT_PRODUCT(RelMom,RelMom))

                if (AbsRelMom > AbsRelMom_Min) then
                   AbsRelMom_Min = AbsRelMom
                   iP_Min = (/iP1,iP2/)
                end if

             end do
          end do

          if (absRelMom_Min < 1e-5) cycle grLoop

          lambda_min = twopi*hbarc/AbsRelMom_Min

          !            write(*,*) 'lambda_min : ',group(iGr,iP_Min(1:2)),lambda_min

          sizeCluster = 2
          Cluster(1:2) = iP_Min

          do iP=1,sizeGroup(iGr)
             if (iP == Cluster(1)) cycle
             if (iP == Cluster(2)) cycle

             RelPos = realparticles(iEns,group(iGr,iP))%position - realparticles(iEns,group(iGr,Cluster(1)))%position
             if (DOT_PRODUCT(RelPos,RelPos).lt.lambda_Min**2) then
                !                  write(*,*) 'insert Cluster; 1'
                sizeCluster = sizeCluster+1
                Cluster(sizeCluster) = iP
                cycle
             end if

             RelPos = realparticles(iEns,group(iGr,iP))%position - realparticles(iEns,group(iGr,Cluster(2)))%position
             if (DOT_PRODUCT(RelPos,RelPos).lt.lambda_Min**2) then
                !                  write(*,*) 'insert Cluster; 2'
                sizeCluster = sizeCluster+1
                Cluster(sizeCluster) = iP
             end if

          enddo

          if (sizeCluster < 50)  AnzColl(sizeCluster) = AnzColl(sizeCluster)+1


          !            write(*,'(i3,f12.3,": Doing a ",i3,"-body collision: ",f12.3," ",100(i5))') iEns,time,sizeCluster,lambda_Min,group(iGr,Cluster(1:sizeCluster))

          ! --- Remove Cluster from Group

          group(iGr,Cluster(1:sizeCluster)) = 0
          sizeGroup(iGr) = sizeGroup(iGr)-sizeCluster

          iP2 = 0
          do iP1=1,sizeGroup(iGr)
1010         iP2 = iP2+1
             if (group(iGr,iP2)==0) goto 1010
             group(iGr,iP1)=group(iGr,iP2)
          end do


          !            write(*,'(i3," [",i3,"] : ",500(i5))') iGr,sizeGroup(iGr), group(iGr,1:sizeGroup(iGr))


          if (sizeGroup(iGr)>1) goto 1003 ! redo, if two or more particles left

       end do grLoop


       !         if (nGroup > 0) stop

    end do EnsLoop


    write(*,'(f12.3, 50(i5))') time, AnzColl(1:20)
    write(89,'(f12.3, 50(i5))') time, AnzColl


  end subroutine TryKodamaClusters



end module KodamaCluster
