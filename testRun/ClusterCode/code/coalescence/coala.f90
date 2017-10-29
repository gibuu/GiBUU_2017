!***************************************************************************
!****m* /coala
! PURPOSE
! Phase space coalescence between nucleons + Statistical Evaporation Model
! NOTES
! * Stability check of produced clusters by means of generalized
!   evaporation model improved.
! * Fermi break-up NOT included. This is done now in the SMM-part of
!   the cluster code.
! * The variables "ifrm", "ipf" and "npar" originate from the old
!   cluster code (Munich-Catania). They have now PUBLIC status, since they
!   are still needed in the HypCoala module (Hyperon-Cluster coalescence).
!   In case that Hyperon-Cluster coalescence is switched on, they are
!   now deallocated in the HypCoala module.
!***************************************************************************
module coala

use typeDefinitions, only : cluster, particle
use WriteStatus,     only : IOControl


PRIVATE

!*************************************************************************
!****g* coala/ifrm
! NAME
! Integer, Allocatable,dimension(:),save,PUBLIC :: ifrm
! PURPOSE
! indicates which particles from the original BUU ParticleVector
! are bounded/free. Important to check particle conservation.
!*************************************************************************
integer,dimension(:),ALLOCATABLE,save,PUBLIC   :: ifrm

!*************************************************************************
!****g* coala/ipf
! NAME
! Integer, Allocatable,dimension(:),save,PUBLIC :: ipf
! PURPOSE
! local variable, needed to construct variable "npar"
!*************************************************************************
integer,dimension(:),ALLOCATABLE,save,PUBLIC   :: ipf

!*************************************************************************
!****g* coala/npar
! NAME
! Integer, Allocatable,dimension(:,:),save,PUBLIC :: npar
! PURPOSE
! indicates which particles from the original BUU ParticleVector
! form the clusters. Important to check particle conservation
! and construct the phase-space coordinates of the clusters.
!*************************************************************************
integer,dimension(:,:),ALLOCATABLE,save,PUBLIC :: npar

PUBLIC :: Coalescence


contains

  !*************************************************************************
  !****s* coala/Coalescence
  ! NAME
  ! subroutine Coalescence
  !
  ! PURPOSE
  ! Phase Space Coalescence from final state of GiBUU
  ! Parameters are R_c and P_c defined in InputClusters
  !
  ! INPUTS
  ! * type(particle),dimension(:) :: ParticleVector  -- GiBUU particle vector
  ! * integer :: particles -- number of particles/event
  ! * integer :: Hyperons  -- number of hyperons per event
  ! * integer :: Mesons    -- number of mesons per event
  ! * logical :: Get_Hyp   -- Switch for hyperclusters (needed only for deallocation)
  ! * logical :: Get_Asy   -- Switch for correct N/Z-Asymmetry of clusters
  ! * Real    :: R_c, P_c  -- Coalescence parameters in coordinate and momentum space
  !
  ! OUTPUT
  ! * type(particle),dimension(:) :: FragmentVector  -- produced final state particles
  ! * integer :: Clusters -- Number of clusters per event
  ! * integer :: Global   -- Number of clusters+free particles per event
  !
  ! NOTES
  ! * In the FragmentVector first the clusters, then the free baryons and
  !   then the mesons are stored. Its dimension is equal to variable Global.
  ! * Please use always Get_Asy=true, otherwise you would get very proton-
  !   neutron-rich unrealistic clusters!!!
  !*************************************************************************
  subroutine Coalescence(ParticleVector,particles,SpectatorPart,&
       &                 Hyperons,Mesons,R_c,P_c, &
       &                 MaxGEM,FragmentVector,Clusters,Global, &
       &                 Get_Hyp,Get_Asy,Get_GEM,ALADIN_Flag)

    use GEM_Evaporation,only : StabilityGEM !evaporation model
    use CheckAsymmetry, only : GEMAsymmetry !precise asymmetry checks

    implicit none

    !----------------------------------------------------------------------
    ! Input-Output Variables
    !----------------------------------------------------------------------
    integer, intent(in)  :: particles,Hyperons,Mesons,MaxGEM,SpectatorPart
    real,    intent(in)  :: R_c,P_c
    logical, intent(in)  :: Get_Hyp,Get_Asy,Get_GEM,ALADIN_Flag
    type(particle) ,dimension(:),intent(inout) :: ParticleVector
    type(cluster)  ,dimension(:),intent(out)   :: FragmentVector
    integer,                     intent(out)   :: Clusters,global
    !----------------------------------------------------------------------
    ! local Variables
    !----------------------------------------------------------------------
    type(Cluster) :: TheFragment
    integer, dimension(particles) :: npart,label,newlabel
    real,    dimension(particles) :: fsort
    real,    dimension(1:3)       :: dist
    integer :: i1,nfnew,jj,massn,protons,status,check1
    integer :: k1,inix,j5,jjj,mnext,nachbarn,chf,ll,ii,j,i,l1,i3
    integer :: bound,free,iGEM,newFragments
    integer :: ll_check,chf_check,nhf_check
    real    :: aint,rinit,aint2,rsf
    real    :: emass,rrx,rry,rrz,pp0,ppx,ppy,ppz
    real    :: rr1,rr2,rr3,rf1,rf2,rf3,rabs
    real    :: rrk1,rrk2,rrk3,px,py,pz,pabs,aaf
    logical :: lfragm,lcc,CheckAsy
    !----------------------------------------------------------------------
    allocate(ifrm(1:particles),STAT=status)
    call IOControl(1,status,'coala','ifrm')
    allocate(ipf(1:particles),STAT=status)
    call IOControl(1,status,'coala','ipf')
    allocate(npar(1:particles,1:particles),STAT=status)
    call IOControl(1,status,'coala','npar')
    !----------------------------------------------------------------------
    nfnew    = 0
    Clusters = 0 !number of clusters in an events
    global   = 0 !number of clusters & singles in an event
    do i=1,particles
       ifrm(i) = 0 !bounded/unbounded particle (1/0)
       npart(i)= 0
       label(i)= 0
       newlabel(i) = 0
       fsort(i) = 0.
    end do
    !//////////////////////////////////////////////////////////////////////
    ! START OF PHASE-SPACE COALESCENCE+EVAPORATION ////////////////////////
    !//////////////////////////////////////////////////////////////////////
    Clusters = nfnew
    !**********************************************************************
    do i1=1,particles !loop over Particles of event ***********************
    !**********************************************************************
       if (ParticleVector(i1)%ID==999) then
          write(*,*) 'from Module Coala, routine Coalescence :'
          write(*,*) 'this should not be happened!!! i1,ID = ', &
               & i1,ParticleVector(i1)%ID
          write(*,*) '!!! TERMINATION OF PROGRAM NOW !!!'
          STOP
       endif

!--> this option (%sd) not used any more!!!
!       if (ALADIN_Flag) then
!          if (ParticleVector(i1)%SD /= 2) CYCLE
!       endif

       lfragm = .false.
       lcc    = .false.
       if (ParticleVector(i1)%id.ne.1) cycle !only nucleons!!!!
       if (ifrm(i1).eq.1) cycle ! particle with label i1 is
                                ! already bounded, cycle particle loop

       if (Clusters.ne.0) & !check if particle i1 can be bounded to existing clusters
            call fff2(i1,particles,ParticleVector,FragmentVector, &
                      Clusters,R_c,P_c,lcc,Get_Asy)
       if (lcc) cycle !particle i1 is bounded to an existing cluster,
                      !therefore cycle particle loop

       jj = 0
       !--------------------------------------------------------------------
       ! First coalescence of particles in coordinate space
       ! The Coalescence radius is R_c
       !--------------------------------------------------------------------
       do k1=1,particles
          if (ParticleVector(k1)%ID==999) then
             write(*,*) 'from Module Coala, routine Coalescence:'
             write(*,*) 'this should not be happened!!! k1,ID = ', &
                  & k1,ParticleVector(i1)%ID
             write(*,*) '!!! TERMINATION OF PROGRAM NOW !!!'
             STOP
          endif

!--> see above
!          if (ALADIN_Flag) then
!             if (ParticleVector(k1)%SD /= 2) CYCLE
!          endif

          if(ParticleVector(k1)%id.ne.1) cycle !only nucleons!!!
          if (k1==i1) cycle
          if (ifrm(k1).eq.1) cycle !only unbounded nucleons
          dist(:)=ParticleVector(i1)%position(:)-ParticleVector(k1)%position(:)
          aint= sqrt( dist(1)**2 + dist(2)**2 + dist(3)**2 )
          if (aint.gt.R_c) cycle
          jj=jj+1
          npart(jj)=k1
          fsort( npart(jj) )=0.0
       end do
       call panic (jj,particles) !check overflows, don't get panic...
       if (jj.eq.0) cycle
       !--------------------------------------------------------------------
       ! sorting of possible coalescence partners acc. their distance
       !--------------------------------------------------------------------
       inix = 0
900    continue
       if (inix.ge.jj) goto 800
       rinit  = 50.
       do j5=1,jj
          jjj   = npart(j5)
          if (fsort(jjj).eq.1.0) cycle
          dist(:) = ParticleVector(i1)%position(:) - ParticleVector(jjj)%position(:)
          aint2   = sqrt( dist(1)**2 + dist(2)**2 + dist(3)**2 )
          if (aint2.lt.rinit) then
             rinit = aint2
             mnext = jjj
          endif
       end do
       inix        = inix + 1
       label(inix) = mnext
       fsort(mnext)= 1.0
       goto 900
800    continue
       call panic (inix,particles) !check overflows, don't get panic...
       nachbarn = inix
       !-------------------------------------------------------------------
       ! Start of Coalescence in Momentum Space
       ! The Coalescence Radius is P_c
       !-------------------------------------------------------------------
       emass     = ParticleVector(i1)%mass               !efm0(i1)
       rrx       = ParticleVector(i1)%position(1)        !x(i1)
       rry       = ParticleVector(i1)%position(2)        !y(i1)
       rrz       = ParticleVector(i1)%position(3)        !z(i1)
       pp0       = ParticleVector(i1)%momentum(0)        !k0(i1)
       ppx       = ParticleVector(i1)%momentum(1)        !k0(i1)
       ppy       = ParticleVector(i1)%momentum(2)        !k0(i1)
       ppz       = ParticleVector(i1)%momentum(3)        !k0(i1)
       chf       = ParticleVector(i1)%charge             !char(i1)
       ll        = 1
       ipf(ll)   = i1
       rsf       = 0.0
       !-------------------------------------------------------------------
       search : do ii=1,nachbarn
          !*****************************
          j = label(ii)
          if (j.le.0) then
             OPEN(96,FILE='error_file')
             write (96,*) 'AUS FRAGMENT: FEHLER BEIM SORTIEREN !!'
             write (96,*) '!!! PROGRAMM-ABBRUCH !!!'
             CLOSE(96)
             STOP
          endif
          !*****************************
          rr1  = rrx
          rr2  = rry
          rr3  = rrz

          rf1  = rr1 - ParticleVector(j)%position(1)
          rf2  = rr2 - ParticleVector(j)%position(2)
          rf3  = rr3 - ParticleVector(j)%position(3)

          rabs = sqrt( rf1*rf1 + rf2*rf2 + rf3*rf3 )

          rrk1 = ppx
          rrk2 = ppy
          rrk3 = ppz

          px   = rrk1-ParticleVector(j)%momentum(1)
          py   = rrk2-ParticleVector(j)%momentum(2)
          pz   = rrk3-ParticleVector(j)%momentum(3)

          pabs = sqrt( px*px + py*py + pz*pz )

          !-----------------------------------------------------------------
          ! * Phase-Space Coalescence
          ! * an approximate cluster radius acc. Rsf=1.3*A^1/3 is defined
          ! * all particles inside this radius do not undergo
          !   momentum coalescence
          ! * For heavier clusters than 2H N/Z-ratio is already here checked
          !-----------------------------------------------------------------
          PhaseSpaceCoalescence : if ( (rabs.le.rsf).or. &
               &  ( (rabs.gt.rsf).and. &
               &    (rabs.lt.(R_c+rsf)).and.(pabs.lt.P_c) ) ) then

             lfragm  = .true.
             ifrm(j) = 1
             ll      = ll + 1
             !errmessage='aus fragment: 120-schleife'
             CALL PANIC (LL,particles)
             aaf     = float(ll)
             if (ll.le.3) then
                rsf = 2.25
             else
                rsf = 1.3*aaf**0.33333333
             endif
             if ( (R_c.eq.0.0).and.(P_c.eq.0.0) ) rsf = 0.0
             ipf(ll) = j

             emass   = emass + ParticleVector(j)%mass

             rrx = (emass-ParticleVector(j)%mass)*rr1/emass + &
                  ParticleVector(j)%mass*ParticleVector(j)%position(1)/emass
             rry = (emass-ParticleVector(j)%mass)*rr2/emass + &
                  ParticleVector(j)%mass*ParticleVector(j)%position(2)/emass
             rrz = (emass-ParticleVector(j)%mass)*rr3/emass + &
                  ParticleVector(j)%mass*ParticleVector(j)%position(3)/emass

             pp0  = ( pp0*float(ll-1) + ParticleVector(j)%momentum(0) )/float(ll)
             ppx  = ( ppx*float(ll-1) + ParticleVector(j)%momentum(1) )/float(ll)
             ppy  = ( ppy*float(ll-1) + ParticleVector(j)%momentum(2) )/float(ll)
             ppz  = ( ppz*float(ll-1) + ParticleVector(j)%momentum(3) )/float(ll)

             chf  = chf + ParticleVector(j)%charge

          endif PhaseSpaceCoalescence

       end do search
       !-------------------------------------------------------------------
       buildFragments : if (lfragm) then

          ifrm(i1)  = 1
          Clusters        = Clusters + 1
          FragmentVector(Clusters)%MassNumber = ll

          do i3=1,FragmentVector(Clusters)%MassNumber
             npar(Clusters,i3) = ipf(i3)
             !----------------------------------------------
             if ( (FragmentVector(Clusters)%MassNumber.gt.particles).or. &
                  & (npar(Clusters,i3).le.0) ) then
                OPEN(96,FILE='error_file')
                write (96,*) 'AUS FRAGMENT:FALSCH !', &
                     Clusters,FragmentVector(Clusters)%MassNumber,npar(Clusters,i3)
                write (96,*) ' !!! PROGRAMM-ABBRUCH !!!'
                CLOSE(96)
                STOP
             endif
             !---------------------------------------------
          end do
          !errmessage='aus fragment: lfragm-schleife'
          call panic (Clusters,particles)
          FragmentVector(Clusters)%position(1)  = rrx
          FragmentVector(Clusters)%position(2)  = rry
          FragmentVector(Clusters)%position(3)  = rrz
          FragmentVector(Clusters)%momentum(0)  = &
               & pp0*float( FragmentVector(Clusters)%MassNumber )
          FragmentVector(Clusters)%momentum(1)  = &
               & ppx*float( FragmentVector(Clusters)%MassNumber )
          FragmentVector(Clusters)%momentum(2)  = &
               & ppy*float( FragmentVector(Clusters)%MassNumber )
          FragmentVector(Clusters)%momentum(3)  = &
               & ppz*float( FragmentVector(Clusters)%MassNumber )
          FragmentVector(Clusters)%Mass         = emass
          FragmentVector(Clusters)%ChargeNumber = chf
          FragmentVector(Clusters)%HypNumber    = 0
          FragmentVector(Clusters)%ID           = 1
          FragmentVector(Clusters)%FreeBound    = .true. !assume it stable

          !-----------------------------------------------------------------
          ! - Check N/Z-Asymmetry for heavy Clusters
          !   Funny clusters with Z=0 or N=0 are checked later on.
          !   These funny clusters are forced to completely decay into
          !   free nucleons.
          ! - Otherwise (invalid N/Z-ratio, but N /=0 and Z/=0) the routine
          !   ConstructStableElement search for the nearest stable element
          !   with respect to the cluster with label m1, and its quantities
          !   are then re-calculated.
          !-----------------------------------------------------------------
          ll_check  = FragmentVector(Clusters)%MassNumber
          chf_check = FragmentVector(Clusters)%ChargeNumber
          nhf_check = ll_check - chf_check
          CheckIso : if (Get_Asy .and. ll_check > 2 .and. &
               &        (chf_check /= 0 .and. nhf_check /= 0) )  then
             TheFragment = FragmentVector(Clusters)
             call ConstructStableElement(Clusters,ll_check,chf_check,&
                  &                      TheFragment,ParticleVector)
             FragmentVector(Clusters) = TheFragment
          endif CheckIso
          !-----------------------------------------------------------------

       endif buildFragments
    !**********************************************************************
    end do !end of ParticlesLoop ******************************************
    !**********************************************************************

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !IN THE FOLLOWING SOME STABILITY CHECKS ARE PERFORMED
    ! (1) De-Excitation:
    !     --> application of statistical evaporation model (GEM)
    ! (2) CheckAymmetry:
    !     FRAGMENTS MUST LIE INSIDE THE (N,Z)-STABILITY REGION
    !     ACCORDING THE NNDC-TABLE (See InitGEM Module)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    bound = 0
    free  = 0
    do i=1,particles
       if (ParticleVector(i)%ID==999) then
          write(*,*) 'from coala : this should not be happened!!! i,ID = ', &
               & i,ParticleVector(i1)%ID
          write(*,*) '!!! TERMINATION NOW !!!'
          STOP
       endif
       if (ParticleVector(i)%id.gt.100) cycle !only Baryons!!!!
!       if (ALADIN_Flag) then
!          if (ParticleVector(i)%SD /= 2) CYCLE
!       endif
       if (ifrm(i).eq.1) bound = bound + 1
       if (ifrm(i).eq.0) free  = free + 1
       if (ifrm(i).lt.0.or.ifrm(i).gt.1) then
          write(*,*) 'ifrm not well determined !!!',ifrm
          write(*,*) '!!! stop the run now !!!'
          STOP
       endif
    end do

    write(*,*) 'Clusters         = ',Clusters
    write(*,*) 'Mass of Clust.   = ',(FragmentVector(i)%MassNumber,i=1,Clusters)
    write(*,*) 'Charge of Clust. = ',(FragmentVector(i)%ChargeNumber,i=1,Clusters)
    write(*,*) 'bounded Baryons  = ',bound
    write(*,*) 'free Baryons     = ',free

    !----------------------------------------------------------------------
    CkeckAsymmetry : IF (Get_Asy) THEN
       nfnew = 0
       CLusterLoop : do i=1,Clusters
          massn   = FragmentVector(i)%MassNumber
          protons = FragmentVector(i)%ChargeNumber
          call GEMAsymmetry(massn,protons,CheckAsy) !precise check on N/Z
          StabilityControl : if (.not.CheckAsy) then
             do j=1,FragmentVector(i)%MassNumber
                ifrm( npar(i,j) ) = 0 !set particles of unbounded cluster to free!
             end do
          else
             nfnew           = nfnew + 1
             newlabel(nfnew) = i
          endif StabilityControl
       END DO CLusterLoop
       Clusters = nfnew
       !--------------------------------------------------------------------
       call panic (nfnew,particles)
       !--------------------------------------------------------------------
       do i=1,Clusters   ! Sort the Clusters into the Global FragmentVector !
          FragmentVector(i)%position(:) = FragmentVector(newlabel(i))%position(:)
          FragmentVector(i)%momentum(:) = FragmentVector(newlabel(i))%momentum(:)
          FragmentVector(i)%MassNumber  = FragmentVector(newlabel(i))%MassNumber
          FragmentVector(i)%ChargeNumber= FragmentVector(newlabel(i))%ChargeNumber
          FragmentVector(i)%Mass        = FragmentVector(newlabel(i))%Mass
          FragmentVector(i)%HypNumber   = 0
          FragmentVector(i)%FreeBound   = .true.
          FragmentVector(i)%ID          = 1
          do j=1,FragmentVector(i)%MassNumber
             npar(i,j) = npar( newlabel(i), j )
          end do
       end do
    endif CkeckAsymmetry
    !----------------------------------------------------------------------
    GEMEvaporation : if (Get_GEM) then
       GEMLoop : do iGEM=1,MaxGEM !max. value of GEM steps
          call stabilityGEM(iGEM,Particles,Clusters,npar,ifrm,&
               &            FragmentVector,ParticleVector, &
               &            newFragments)
          GEMAsyCheck : if (MaxGEM > 1) then
             nfnew = 0
             do i=1,Clusters
                massn   = FragmentVector(i)%MassNumber
                protons = FragmentVector(i)%ChargeNumber
                call GEMAsymmetry(massn,protons,CheckAsy) !Again precise check on N/Z
                if (.not.CheckAsy) then
                   write(*,*) 'this should actually NOT happened!!!!!'
                   STOP
                   do j=1,FragmentVector(i)%MassNumber
                      ifrm( npar(i,j) ) = 0 !set particles of unbounded cluster to free!
                   END DO
                else
                   nfnew           = nfnew + 1
                   newlabel(nfnew) = i
                endif
             END DO
             Clusters = nfnew
             !--------------------------------------------------------------
             call panic (nfnew,particles)
             !--------------------------------------------------------------
             do i=1,Clusters   ! Sort the Clusters into the Global FragmentVector !
                FragmentVector(i)%position(:) = FragmentVector(newlabel(i))%position(:)
                FragmentVector(i)%momentum(:) = FragmentVector(newlabel(i))%momentum(:)
                FragmentVector(i)%MassNumber  = FragmentVector(newlabel(i))%MassNumber
                FragmentVector(i)%ChargeNumber= FragmentVector(newlabel(i))%ChargeNumber
                FragmentVector(i)%Mass        = FragmentVector(newlabel(i))%Mass
                FragmentVector(i)%HypNumber   = 0
                FragmentVector(i)%FreeBound   = .true.
                FragmentVector(i)%ID          = 1
                do j=1,FragmentVector(i)%MassNumber
                   npar(i,j) = npar( newlabel(i), j )
                end do
             end do
          endif GEMAsyCheck
          if (newFragments==0) EXIT !GEM loop when all clusters are stable
       end do GEMLoop
    endif GEMEvaporation
    !--------------------------------------------------------------------
    call panic (nfnew,particles)
    !//////////////////////////////////////////////////////////////////////
    ! END OF PHASE-SPACE COALESCENCE+EVAPORATION //////////////////////////
    !//////////////////////////////////////////////////////////////////////
    bound = 0
    free  = 0
    global = Clusters
    do l1=1,Particles ! ! Sort Free Baryons into the Global FragmentVector ! !
       if (ParticleVector(l1)%ID==999) then
          write(*,*) 'from coala : this should not be happened!!! l1,ID = ', &
               & l1,ParticleVector(i1)%ID
          write(*,*) '!!! TERMINATION NOW !!!'
          STOP
       endif
!       if (ALADIN_Flag) then
!          if (ParticleVector(l1)%SD /= 2) CYCLE
!       endif
       if(ParticleVector(l1)%id.gt.33.and.ParticleVector(l1)%id.lt.100) cycle
!       if (ParticleVector(l1)%id.gt.200) cycle
       if (ParticleVector(l1)%id.lt.4) then
          if (ifrm(l1)==1) bound = bound + 1
          if (ifrm(l1)==0) free  = free  + 1
       endif
       if (ifrm(l1) == 1) cycle
       global = global + 1
       FragmentVector(global)%position(:) = ParticleVector(l1)%position(:)
       FragmentVector(global)%momentum(:) = ParticleVector(l1)%momentum(:)
       FragmentVector(global)%Mass        = ParticleVector(l1)%Mass
       FragmentVector(global)%MassNumber  = 0   !
       FragmentVector(global)%HypNumber   = 0   !
       FragmentVector(global)%ChargeNumber= ParticleVector(l1)%Charge !
       FragmentVector(global)%FreeBound   = .false.    !true=cluster/false=free nucleon
       FragmentVector(global)%ID          = ParticleVector(l1)%ID
    end do

    write(*,*)
    write(*,*) '               after StabilityChecks     '
    write(*,*)
    write(*,*) 'Clusters         = ',Clusters
    write(*,*) 'Mass of Clust.   = ',(FragmentVector(i)%MassNumber,i=1,Clusters)
    write(*,*) 'Charge of Clust. = ',(FragmentVector(i)%ChargeNumber,i=1,Clusters)
    write(*,*) 'bounded Nucleons = ',bound
    write(*,*) 'free Nucleons    = ',free
    write(*,*) 'Hyperons         = ',Hyperons
    write(*,*) 'Mesons           = ',Mesons
    write(*,*)
    !------------------------------------------------------------------------
    ! check particle conservation
    !------------------------------------------------------------------------
    check1 = 0
    do i=1,global
       if (FragmentVector(i)%FreeBound) then
          do j=1,FragmentVector(i)%MassNumber
             check1 = check1 + 1
          end do
       else
          check1 = check1 + 1
       endif
    end do

!    if ((check1+Mesons+Hyperons).ne.Particles) then
    if (.not.ALADIN_Flag) then
       if (check1.ne.Particles) then
          write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
          write(*,*) '!!!! NO PARTICLE CONSERVATION !!!!!'
          write(*,*) 'check1, Particles = ',check1,Particles
          write(*,*) '!!!! PROGRAM TERMINATION !!!!!!!!!!!'
          STOP
       endif
    else
       if (check1.ne.SpectatorPart) then
          write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
          write(*,*) '!!!! NO PARTICLE CONSERVATION in Spectator fragmentation !!!!!'
          write(*,*) 'check1, Particles = ',check1,Particles
          write(*,*) '!!!! PROGRAM TERMINATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
          STOP
       endif
    endif
    !------------------------------------------------------------------------
    ! Delallocation here, if Get_Hyp = .false.
    ! If Get_Hyp = .true., these fields are needed also in Hypcoala,
    ! so they are deallocaled in HypCoala
    !------------------------------------------------------------------------
    if ( ( .not.Get_Hyp ) .or. &
       & ( Get_Hyp .and. (Clusters==0.or.Hyperons==0) ) ) then
    !------------------------------------------------------------------------
       deallocate(ifrm,STAT=status)
       call IOControl(2,status,'coala','ifrm')
       deallocate(ipf,STAT=status)
       call IOControl(2,status,'coala','ipf')
       deallocate(npar,STAT=status)
       call IOControl(2,status,'coala','npar')
    !------------------------------------------------------------------------
    endif
    !------------------------------------------------------------------------
  !**************************************************************************
  end subroutine Coalescence !***********************************************
  !**************************************************************************

  !********************************************************************
  ! just check for overflows...
  !********************************************************************
  subroutine panic(idummy1,idummy2)
    implicit none
    integer ierr
    integer, intent(in) :: idummy1,idummy2
    if ( (idummy1.lt.0).or.(idummy1.gt.idummy2) ) then
       open(96,file='error_file',IOSTAT=ierr)
       write (96,*)
       write (96,*) '******************************'
       write (96,*)
       write (96,*)
       write (96,*)
       write (96,*) ' !!! PROGRAMM STOP !!!'
       write (96,*)
       write (96,*) '******************************'
       close(96,IOSTAT=ierr)
       STOP
    endif
  end subroutine panic
  !********************************************************************

  !*************************************************************************
  !****s* coala/fff2
  ! NAME
  ! subroutine fff2
  !
  ! PURPOSE
  ! Checks coalescence between nucleons and already produced clusters
  !
  ! INPUTS
  ! * type(particle),dimension(:) :: ParticleVector  -- GiBUU particle vector
  ! * integer :: particles       -- number of particles/event
  ! * integer :: i1              -- particle with label i1
  ! * integer :: Clusters        -- Number of clusters per event
  !
  ! OUTPUT
  ! * type(particle),dimension(:) :: FragmentVector  -- produced final state particles
  ! * logical :: lcc -- true in case of nucleon+cluster-Coalescence
  !
  ! NOTES
  ! FragmentVector enters as input, in case of coalescence with a free nucleon
  ! its quantities are re-calculated.
  !*************************************************************************
  subroutine fff2(i1,particles,ParticleVector,FragmentVector, &
                  Clusters,R_c,P_c,lcc,Get_Asy)
  !********************************************************************
    implicit none
    !-----------------------------------------------------------------*
    ! Input-Output Variables
    !-----------------------------------------------------------------*
    integer,                     intent(in)    :: i1,particles,Clusters
    real,                        intent(in)    :: R_c,P_c
    type(particle), dimension(:),intent(in)    :: ParticleVector
    logical,                     intent(in)    :: Get_Asy
    type(cluster), dimension(:), intent(inout) :: FragmentVector
    logical,                     intent(inout) :: lcc
    !-----------------------------------------------------------------*
    ! Local Variables
    !-----------------------------------------------------------------*
    type(Cluster) :: TheFragment
    integer :: m1,l1,ll,chh,chnew,ll_check,chf_check,nhf_check
    real    :: rinit,r1,r2,r3,rabs,k1h,k2h,k3h,emh,rp1,rp2,rp3
    real    :: pabs,an,rs,xnew,ynew,znew,pp0,ppx,ppy,ppz,emnew
    !-----------------------------------------------------------------*
    rinit = 500. !fm
    m1    = 0

    lcc = .false.
    ifrm(i1) = 0
    !just to be sure...
    if (Clusters.eq.0) then
       return
    endif
    !choose the closest fragment
    if (Clusters.gt.1) then
       do l1=1,Clusters
          r1 = ParticleVector(i1)%position(1)-FragmentVector(l1)%position(1)
          r2 = ParticleVector(i1)%position(2)-FragmentVector(l1)%position(2)
          r3 = ParticleVector(i1)%position(3)-FragmentVector(l1)%position(3)
          rabs = sqrt(r1**2+r2**2+r3**2)
          if (rabs.lt.rinit) then
             rinit = rabs
             m1 = l1
          endif
       end do
    else
       m1 = Clusters
    endif
    !------------------------------------------------------------------*
    if (m1.eq.0) return
    !****************************************************
    k1h = FragmentVector(m1)%momentum(1)/float(FragmentVector(m1)%MassNumber)
    k2h = FragmentVector(m1)%momentum(2)/float(FragmentVector(m1)%MassNumber)
    k3h = FragmentVector(m1)%momentum(3)/float(FragmentVector(m1)%MassNumber)

    emh = FragmentVector(m1)%Mass
    chh = FragmentVector(m1)%ChargeNumber

    r1 = ParticleVector(i1)%position(1)-FragmentVector(m1)%position(1)
    r2 = ParticleVector(i1)%position(2)-FragmentVector(m1)%position(2)
    r3 = ParticleVector(i1)%position(3)-FragmentVector(m1)%position(3)
    rabs = sqrt(r1**2+r2**2+r3**2)

    rp1 = ParticleVector(i1)%momentum(1)-k1h
    rp2 = ParticleVector(i1)%momentum(2)-k2h
    rp3 = ParticleVector(i1)%momentum(3)-k3h
    pabs = sqrt(rp1**2+rp2**2+rp3**2)

    an = float( FragmentVector(m1)%MassNumber ) + 1.0
!    rs = 1.3*an**0.3333333
    if (an.le.3.) then
       rs = 2.25
    else
       rs = 1.3*an**0.33333333
    endif
    if ( (R_c.eq.0.0).and.(P_c.eq.0.0) ) rs = 0.0

    !----------- Abfrage zur Fragmentbildung --------------------------------*

!    if ( (rabs.gt.R_c).or.(pabs.gt.P_c) ) return

    if (rabs.gt.(rs+R_c)) return
    if ( (rabs.gt.rs).and.(rabs.le.(R_c+rs)) ) then
       if (pabs.gt.P_c) return
    endif

    lcc       = .true.
    ifrm(i1)  = 1
    ll        = FragmentVector(m1)%MassNumber + 1
    ipf(ll)   = i1
    call panic (ll,particles)

    xnew  = ( ParticleVector(i1)%mass*ParticleVector(i1)%position(1) &
         + emh*FragmentVector(m1)%position(1) )/ &
         ( ParticleVector(i1)%mass + emh )
    ynew  = ( ParticleVector(i1)%mass*ParticleVector(i1)%position(2) &
         + emh*FragmentVector(m1)%position(2) )/ &
         ( ParticleVector(i1)%mass + emh )
    znew  = ( ParticleVector(i1)%mass*ParticleVector(i1)%position(3) &
         + emh*FragmentVector(m1)%position(3) )/ &
         ( ParticleVector(i1)%mass + emh )

    pp0   = ParticleVector(i1)%momentum(0) + FragmentVector(m1)%momentum(0)
    ppx   = ( ParticleVector(i1)%momentum(1) + &
         k1h*float(FragmentVector(m1)%MassNumber) ) / float(ll)
    ppy   = ( ParticleVector(i1)%momentum(2) + &
         k2h*float(FragmentVector(m1)%MassNumber) ) / float(ll)
    ppz   = ( ParticleVector(i1)%momentum(3) + &
         k3h*float(FragmentVector(m1)%MassNumber) ) / float(ll)

    emnew = ParticleVector(i1)%mass + emh     !new fragment mass
    chnew = ParticleVector(i1)%charge + chh   !new fragment charge

    !--------------- neue schwerpunktskoopdimaten ---------------------------

    FragmentVector(m1)%MassNumber   = ll

    npar(m1,FragmentVector(m1)%MassNumber) = ipf(FragmentVector(m1)%MassNumber)

    FragmentVector(m1)%position(1)  = xnew
    FragmentVector(m1)%position(2)  = ynew
    FragmentVector(m1)%position(3)  = znew
    FragmentVector(m1)%momentum(0)  = pp0
    FragmentVector(m1)%momentum(1)  = ppx*float(FragmentVector(m1)%MassNumber)
    FragmentVector(m1)%momentum(2)  = ppy*float(FragmentVector(m1)%MassNumber)
    FragmentVector(m1)%momentum(3)  = ppz*float(FragmentVector(m1)%MassNumber)
    FragmentVector(m1)%Mass         = emnew
    FragmentVector(m1)%ChargeNumber = chnew
    FragmentVector(m1)%FreeBound    = .true. !assume it to be stable

    !-----------------------------------------------------------------
    ! - Check N/Z-Asymmetry for heavy Clusters
    !   Funny clusters with Z=0 or N=0 are checked later on.
    !   These funny clusters are forced to completely decay into
    !   free nucleons.
    ! - Otherwise (invalid N/Z-ratio, but N /=0 and Z/=0) the routine
    !   ConstructStableElement search for the nearest stable element
    !   with respect to the cluster with label m1, and its quantities
    !   are then re-calculated.
    !-----------------------------------------------------------------
    ll_check  = FragmentVector(m1)%MassNumber
    chf_check = FragmentVector(m1)%ChargeNumber
    nhf_check = ll_check - chf_check

    CheckIso : if (Get_Asy .and. ll_check > 2 .and. &
         &        (chf_check /= 0 .and. nhf_check /= 0) )  then
       TheFragment = FragmentVector(m1)
       call ConstructStableElement(m1,ll_check,chf_check,TheFragment,&
            &                      ParticleVector)
       FragmentVector(m1) = TheFragment
    endif CheckIso
    !-----------------------------------------------------------------

  !********************************************************************
  end subroutine fff2 !************************************************
  !********************************************************************


  !********************************************************************
  subroutine ConstructStableElement(FraLabel,ll_check,chf_check,TheFragment,&
       &                            ParticleVector)
  !********************************************************************
    use CheckAsymmetry, only : GEMAsymmetry,FindNearestStableElement
    implicit none
    !-----------------------------------------------------------------*
    ! Input-Output Variables
    !-----------------------------------------------------------------*
    type(particle), dimension(:),intent(in)    :: ParticleVector
    integer,                     intent(in)    :: ll_check,chf_check
    integer,                     intent(in)    :: FraLabel
    type(Cluster),               intent(inout) :: TheFragment
    !-----------------------------------------------------------------*
    ! Local Variables
    !-----------------------------------------------------------------*
    logical, dimension(:),allocatable :: setItFree
    integer, dimension(:),allocatable :: ipf_new
    real, dimension(1:3) :: xnew
    real, dimension(0:3) :: pnew
    integer :: status
    integer :: NewN,NewZ,Ndummy,Pdummy,i3,chf_new,ll_new,FraDim
    logical :: CheckAsy
    !-----------------------------------------------------------------*
    call GEMAsymmetry(ll_check,chf_check,CheckAsy)
    if (.not.CheckAsy) then
       call FindNearestStableElement(ll_check,chf_check,&
            & NewN,NewZ)
       Ndummy  = 0
       Pdummy  = 0
       allocate(setItFree(1:ll_check),STAT=status)
       call IOControl(1,status,'coala','setItFree')
       do i3=1,ll_check
          setItFree(i3) = .true.
          if (ParticleVector(npar(FraLabel,i3))%Charge==0 .and. Ndummy .le. NewN) then
             Ndummy        = Ndummy + 1
             setItFree(i3) = .false.
          endif
          if (ParticleVector(npar(FraLabel,i3))%Charge==1 .and. Pdummy .le. NewZ) then
             Pdummy        = Pdummy + 1
             setItFree(i3) = .false.
          endif
       end do
       chf_new = 0
       FraDim  = 0
       do i3=1,ll_check
          if (setItFree(i3)) then
             ifrm(npar(FraLabel,i3)) = 0
          else
             ifrm(npar(FraLabel,i3)) = 1
             chf_new = chf_new + ParticleVector(npar(FraLabel,i3))%Charge
             FraDim  = FraDim + 1
          endif
       end do
       ll_new  = 0
       allocate(ipf_new(1:FraDim),STAT=status)
       call IOControl(1,status,'coala','ipf_new')
       do i3=1,ll_check
          if(.not.setItFree(i3)) then
             ll_new          = ll_new + 1
             ipf_new(ll_new) = npar(FraLabel,i3)
          end if
       end do
       !build the fragment & re-calculate its PhaseSpace Coordinates
       TheFragment%MassNumber   = ll_new
       TheFragment%ChargeNumber = chf_new
       TheFragment%Mass         = float(ll_new)*0.939/0.19733
       TheFragment%HypNumber    = 0
       TheFragment%FreeBound    = .true.
       xnew(:) = 0.0
       pnew(:) = 0.0
       do i3=1,ll_new
          npar(FraLabel,i3) = ipf_new(i3)
          !----------------------------------------------
          if ( (TheFragment%MassNumber.gt.size(ParticleVector)).or. &
               & (npar(FraLabel,i3).le.0) ) then
             OPEN(96,FILE='error_file')
             write (96,*) 'AUS FRAGMENT:FALSCH !', &
                  FraLabel,TheFragment%MassNumber,npar(FraLabel,i3)
             write (96,*) ' !!! PROGRAMM-ABBRUCH !!!'
             CLOSE(96)
             STOP
          endif
          !---------------------------------------------
          xnew(:) = xnew(:) + &
               & ParticleVector(ipf_new(i3))%position(:) * &
               & ParticleVector(ipf_new(i3))%Mass
          pnew(:) = pnew(:) + &
               & ParticleVector(ipf_new(i3))%momentum(:)
       end do
       TheFragment%position(:) =  xnew(:)/TheFragment%mass
       TheFragment%momentum(:) = pnew(:)

       deallocate(setItFree,STAT=status)
       call IOControl(2,status,'coala','setItFree')
       deallocate(ipf_new,STAT=status)
       call IOControl(2,status,'coala','ipf_new')
    endif

  !*************************************************************************
  end subroutine ConstructStableElement
  !*************************************************************************





end module coala
