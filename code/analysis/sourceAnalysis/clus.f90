!******************************************************************************
!****m* /clus
! NAME
! module clus
!
! PURPOSE
! Includes the programs which perform the search for clusters in coordinate
! (and momentum) space.
! NOTES
!******************************************************************************


 module clus

  implicit none
  private


  public :: DoClustering

 contains


  subroutine DoClustering(numEnsembles,numParticles,particleVector,hyperSource,sourceType,NumSources)

    use particleDefinition
    use IdTable, only: nucleon,Lambda,SigmaResonance

    integer,                        intent(in) :: numEnsembles,numParticles
    type(particle), dimension(:,:), intent(in) :: particleVector
    logical,                        intent(in) :: hyperSource
    integer, dimension(1:numEnsembles),                intent(out) :: NumSources
    integer, dimension(1:numEnsembles,1:numParticles), intent(out) :: sourceType

    real, parameter :: dcl=3., dcl2=dcl**2
    real, dimension(1:3) :: dr
    real :: dr2
    integer, dimension(1:numParticles) :: ind
    integer :: i,j,nclus,ia,n1,n2,k,n,jj
    logical :: flagFirst

    sourceType=999

    Loop_over_ensembles : do i = 1,numEnsembles

        ! find first particle in a cluster:
        flagFirst=.false.
        do j=1,numParticles

          if (particleVector(i,j)%id<=0) cycle
          if (particleVector(i,j)%antiparticle) cycle
          if (.not.hyperSource) then
              if (ParticleVector(i,j)%ID .ne. nucleon) cycle ! only nucleons
          else
              if (       ParticleVector(i,j)%ID.ne.nucleon &
                 & .and. ParticleVector(i,j)%ID.ne.Lambda &
                 & .and. ParticleVector(i,j)%ID.ne.SigmaResonance &
                 & .or.  ParticleVector(i,j)%ID.eq.SigmaResonance &
                 & .and. ParticleVector(i,j)%charge.ne.0     ) cycle  ! only nucleons, Lambda's
                                                                      ! and Sigma0's
          end if

          flagFirst=.true.
          exit

        end do

        if (.not.flagFirst) cycle

        nclus=0
        loop_over_clusters : do

           nclus=nclus+1
           sourceType(i,j)=nclus
           ia=1
           n1=1
           n2=1
           ind(1)=j

           loop_over_groups_of_particles : do
!              (each group includes particles
!               close to the particles from previous group)

              k=0
              loop_over_particles_from_group : do n=n1,n2

                 search_loop : do jj=j+1,numParticles

                    if (particleVector(i,jj)%id<=0) cycle
                    if (particleVector(i,jj)%antiparticle) cycle

                    if (.not.hyperSource) then
                       if (ParticleVector(i,jj)%ID .ne. nucleon) cycle !only nucleons
                    else
                       if (      ParticleVector(i,jj)%ID.ne.nucleon &
                         & .and. ParticleVector(i,jj)%ID.ne.Lambda &
                         & .and. ParticleVector(i,jj)%ID.ne.SigmaResonance &
                         & .or.  ParticleVector(i,jj)%ID.eq.SigmaResonance &
                         & .and. ParticleVector(i,jj)%charge.ne.0     ) cycle  ! only nucleons, Lambda's
                                                                               ! and Sigma0's
                    end if

                    if (sourceType(i,jj).ne.999) cycle

                    dr(1:3)=particleVector(i,jj)%position(1:3)-particleVector(i,ind(n))%position(1:3)
                    dr2=dot_product(dr(1:3),dr(1:3))

                    if (dr2 <= dcl2) then
                       sourceType(i,jj)=nclus
                       ia=ia+1
                       ind(ia)=jj
                       k=k+1
                    end if

                 end do search_loop

              end do loop_over_particles_from_group

              if (k.gt.0) then
                 n1=n2+1
                 n2=ia
                 cycle loop_over_groups_of_particles
              end if

              exit

           end do loop_over_groups_of_particles

           ! test whether there are particles
           ! not included into clusters yet:

           do jj=j+1,numParticles
              if (particleVector(i,jj)%id<=0) cycle
              if (particleVector(i,jj)%antiparticle) cycle

              if (.not.hyperSource) then
                  if (ParticleVector(i,jj)%ID .ne. nucleon) cycle !only nucleons
              else
                  if (      ParticleVector(i,jj)%ID.ne.nucleon &
                    & .and. ParticleVector(i,jj)%ID.ne.Lambda &
                    & .and. ParticleVector(i,jj)%ID.ne.SigmaResonance &
                    & .or.  ParticleVector(i,jj)%ID.eq.SigmaResonance &
                    & .and. ParticleVector(i,jj)%charge.ne.0     ) cycle  ! only nucleons, Lambda's
                                                                          ! and Sigma0's
              end if

              if (sourceType(i,jj)==999) then
                 j=jj
                 cycle loop_over_clusters
              end if
           end do

           exit

        end do loop_over_clusters

        NumSources(i)=nclus

    end do Loop_over_ensembles

  end subroutine DoClustering

 end module clus
