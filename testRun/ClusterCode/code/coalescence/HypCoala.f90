!***************************************************************************
!****m* /HypCoalescence
! PURPOSE
! Phase space coalescence between hyperons and clusters
! NOTES
! * Done with simple coalescence, no criteria on hypernuclear structure!!!
! * It has very similar structure as the coala subroutine
!***************************************************************************
module HypCoalescence

use typeDefinitions, only : Cluster, particle
use coala,           only : ifrm,ipf,npar

PRIVATE
PUBLIC :: HypClusters


contains

  !*************************************************************************
  !****s* HypCoalescence/HypClusters
  ! NAME
  ! subroutine HypClusters
  !
  ! PURPOSE
  ! Phase Space Coalescence between Clusters and Hyperons
  !
  ! INPUTS
  ! * type(particle),dimension(:) :: ParticleVector  -- GiBUU particle vector
  ! * integer :: particles-- number of particles/event
  ! * integer :: Hyperons -- number of hyperons per event
  ! * integer :: Clusters -- Number of clusters per event
  ! * Real    :: R_c, P_c -- The coalescence parameters in coordinate and momentum space
  !
  ! OUTPUT
  ! * type(particle),dimension(:) :: FragmentVector  -- produced final state particles
  ! * integer :: Global  -- Number of clusters+free particles per event
  ! * logical :: FormHyp -- only important for output, if FormHyp=false, no
  !                         output-info on hyperclusters is written
  !
  ! NOTES
  ! In the FragmentVector first the clusters, then the free baryons and
  ! then the mesons are stored. Its dimension is equal to variable Global.
  !*************************************************************************
  subroutine HypClusters(ParticleVector,particles,Hyperons, &
                         Clusters,Global,R_c,P_c,FragmentVector, &
                         FormHyp)
    use WriteStatus, only : IOControl
    implicit none

    !----------------------------------------------------------------------
    integer, intent(in) :: particles,Hyperons,Clusters
    real, intent(in) :: R_c,P_c
    type(particle) ,dimension(:),intent(in) :: ParticleVector
    type(cluster)  ,dimension(:),intent(inout) :: FragmentVector
    integer, intent(inout) :: Global
    logical :: lcc
    logical, intent(out) :: FormHyp
    integer :: i1,l1,status
    !//////////////////////////////////////////////////////////////////////
    ! START OF PHASE-SPACE COALESCENCE ////////////////////////////////////
    !//////////////////////////////////////////////////////////////////////
    if (Clusters.eq.0) return
    if (Hyperons.eq.0) return
    FormHyp = .false.
    !**********************************************************************
    do i1=1,particles !loop over Hyperons of event ************************
    !**********************************************************************
       lcc = .false.
       if ( ParticleVector(i1)%id.lt.32 .or. &
            ParticleVector(i1)%id.gt.33 ) cycle !only Hyperons
       !-------------------------------------------------------------------
       if (ifrm(i1).eq.1) cycle !Hyperon is already bounded
       call FormHypClusters(i1,ParticleVector,FragmentVector, &
                            Clusters,R_c,P_c,lcc)
       if (lcc) FormHyp = .true.
    !**********************************************************************
    end do !end of ParticlesLoop ******************************************
    !**********************************************************************
    Global = Clusters
    do l1=1,Particles ! ! Sort singles into the GlobalVector ! !
       if (ParticleVector(l1)%id.gt.200) cycle !only Baryons!!!!
       if (ifrm(l1)==1) cycle !Baryons bounded to Clusters
       global = global + 1
       FragmentVector(global)%position(1) = ParticleVector(l1)%position(1)
       FragmentVector(global)%position(2) = ParticleVector(l1)%position(2)
       FragmentVector(global)%position(3) = ParticleVector(l1)%position(3)
       FragmentVector(global)%momentum(0) = ParticleVector(l1)%momentum(0)
       FragmentVector(global)%momentum(1) = ParticleVector(l1)%momentum(1)
       FragmentVector(global)%momentum(2) = ParticleVector(l1)%momentum(2)
       FragmentVector(global)%momentum(3) = ParticleVector(l1)%momentum(3)
       FragmentVector(global)%Mass        = ParticleVector(l1)%Mass
       FragmentVector(global)%ID          = ParticleVector(l1)%ID
       FragmentVector(global)%ChargeNumber= ParticleVector(l1)%Charge
       FragmentVector(global)%MassNumber  = 0   !
       FragmentVector(global)%HypNumber   = 0
       FragmentVector(global)%ChargeNumber= ParticleVector(l1)%Charge
       FragmentVector(global)%FreeBound   = .false. !free Baryons
    end do
    !------------------------------------------------------------------------
    ! Deallocation of Global Fields needed for the Coalescence procedure
    !------------------------------------------------------------------------
    deallocate(ifrm,STAT=status)
    call IOControl(2,status,'HypClusters','ifrm')
    deallocate(ipf,STAT=status)
    call IOControl(2,status,'HypClusters','ipf')
    deallocate(npar,STAT=status)
    call IOControl(2,status,'HypClusters','npar')
  !**************************************************************************
  end subroutine HypClusters  !**********************************************
  !**************************************************************************


  !*************************************************************************
  !****s* HypCoalescence/FormHypClusters
  ! NAME
  ! subroutine FormHypClusters
  !
  ! PURPOSE
  ! Checks the criteria for coalescence between hyperons and clusters
  !
  ! INPUTS
  ! * type(particle),dimension(:) :: ParticleVector  -- GiBUU particle vector
  ! * integer :: Clusters -- Number of clusters per event
  ! * Real    :: R_c, P_c -- The coalescence parameters in coordinate and momentum space
  !
  ! OUTPUT
  ! * type(particle),dimension(:) :: FragmentVector  -- produced final state particles
  ! * logical :: lcc -- true, if hypercluster has been formed
  !
  ! NOTES
  ! More hypernuclear structure would be important here...
  ! --> SMM part of the cluster code.
  !*************************************************************************

  subroutine FormHypClusters(i1,ParticleVector,FragmentVector, &
                              Clusters,R_c,P_c,lcc)
  !global variables: Clusters, ifrm(:),ipf(:),npar(:,:)
  !********************************************************************
    implicit none

    integer, intent(in) :: i1,Clusters
    real, intent(in) :: R_c,P_c
    type(particle), dimension(:),intent(in) :: ParticleVector
    type(cluster), dimension(:),intent(inout) :: FragmentVector
    logical, intent(inout) :: lcc
    integer :: m1,l1,ll,chh,chnew,lhyp
    real :: rinit,r1,r2,r3,rabs
    real :: k1h,k2h,k3h,emh
    real :: rp1,rp2,rp3,pabs,an,rs
    real :: xnew,ynew,znew
    real :: pp0,ppx,ppy,ppz,emnew,RadiusFra
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

    if (FragmentVector(m1)%ID==0 .or. FragmentVector(m1)%MassNumber==0 ) then
       write(*,*) 'from FormHypClusters: cluster not defined!!! ID,Mass = ', &
            & FragmentVector(m1)%ID,FragmentVector(m1)%MassNumber
       write(*,*) '!!! TERMINATION NOW !!!'
       STOP
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
    rs = 1.3*an**0.3333333
    if ( (R_c.eq.0.0).and.(P_c.eq.0.0) ) rs = 0.0

    !----------- Abfrage zur Fragmentbildung --------------------------------*

    RadiusFra = 1.3*float( FragmentVector(m1)%MassNumber )**0.333333

    !Coalescence in X-space (including radius of cluster, which is assumed to be
    !in its ground state)
    if (R_c .ne. 0.) then
       if (rabs > (RadiusFra+R_c)) return
    endif
    !Coalescence in P-space
!    if ( (rabs >RadiusFra.and.rabs<(RadiusFra+R_c)) .and. (pabs.gt.P_c) ) return

!    !Old prescriptions, not used any more:
!    if ( (rabs.gt.(R_c+rs)).or.(pabs.gt.P_c) ) return
!    if (rabs.gt.(rs+R_c)) return
!    if ( (rabs.gt.rs).and.(rabs.le.R_c) ) then
!       if (pabs.gt.P_c) return
!    endif
    !*******************
    lcc       = .true.
    ifrm(i1)  = 1 !Hyperon with label i1 is bounded to a Cluster with label m1

    lhyp    = FragmentVector(m1)%HypNumber  + 1 !Hyperon content of the cluster

    ll      = FragmentVector(m1)%MassNumber + 1 !Mass Number
    ipf(ll) = i1

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

    FragmentVector(m1)%position(1) = xnew
    FragmentVector(m1)%position(2) = ynew
    FragmentVector(m1)%position(3) = znew

    FragmentVector(m1)%momentum(0) = pp0
    FragmentVector(m1)%momentum(1) = ppx*float(FragmentVector(m1)%MassNumber)
    FragmentVector(m1)%momentum(2) = ppy*float(FragmentVector(m1)%MassNumber)
    FragmentVector(m1)%momentum(3) = ppz*float(FragmentVector(m1)%MassNumber)
    FragmentVector(m1)%Mass        = emnew

    FragmentVector(m1)%HypNumber   = lhyp
    FragmentVector(m1)%ChargeNumber= chnew

    if (lhyp==1) then !single Y-Hypernuclei (Y='Lambda','Sigma')
       if (particleVector(i1)%ID==32) FragmentVector(m1)%HypType = 'L'
       if (particleVector(i1)%ID==33) FragmentVector(m1)%HypType = 'S'
    else
       FragmentVector(m1)%HypType = 'M' !multi Y-Hypernuclei
    endif

  !********************************************************************
  end subroutine FormHypClusters   !***********************************
  !********************************************************************


end module HypCoalescence
