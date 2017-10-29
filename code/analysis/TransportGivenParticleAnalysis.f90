!******************************************************************************
!****m* /transportGivenParticleAnalysis
! NAME
! module transportGivenParticleAnalysis
!
! PURPOSE
! This module does the analysis of the output of eventtype
! transportGivenParticle
!******************************************************************************

module transportGivenParticleAnalysis
  implicit none
  private

  public:: transportGivenParticle_analyze

  ! plot z coordinate instead of radius:
  logical, parameter :: plot_z_coord = .true.

contains

  subroutine transportGivenParticle_analyze(Particles,timestep,id_in,freq_in)
    use particleDefinition
    use particleProperties, only: hadron
    use IDTable, only: nucleon, isMeson
    use inputGeneral, only: numEnsembles,numTimeSteps,num_Runs_SameEnergy
    use nucleusDefinition
    use nucleus, only: getTarget
    use histf90
    use hist2Df90
    use output, only: intToChar
    use minkowski, only: abs4

    type(particle), intent(in),dimension(:,:) ,target :: Particles
    integer, intent(in) :: timestep
    integer, intent(in), optional :: id_in, freq_in

    type(histogram),save,dimension(:),allocatable:: SF
    type(histogram2D),save,dimension(:),allocatable:: SF_vsRadius

    real, save :: maxmass, minmass
    integer, save :: id,frequency,Nsteps
    logical, save :: initflag = .true.
    real, save :: Ninv
    type(tnucleus),pointer :: nuc
    integer :: i,j,k

    if (initflag) then
       if (present(id_in)) then
          id = id_in
       else
          stop "no particle ID in TGP analysis"
       end if
       if (present(freq_in)) then
          frequency = freq_in
       else
          stop "no frequency in TGP analysis"
       end if
       if (id.eq.nucleon) then
          maxmass = 1.5
       else
          maxmass = 3.0
       end if
       if (isMeson(ID)) then
          minmass = 0.
       else
          minmass = hadron(ID)%minmass-0.2
       end if
       nuc => getTarget()
       Ninv = 1./float(numEnsembles * nuc%mass * num_Runs_SameEnergy)

       Nsteps = numTimeSteps/frequency
       allocate(SF(0:Nsteps))
       allocate(SF_vsRadius(0:Nsteps))
       do k=0,Nsteps
          call CreateHist(SF(k), &
               'reconstructed spectral function at timestep '//intTochar(k), &
               minmass,maxmass,0.005)
          call CreateHist2D(SF_vsRadius(k), &
               'reconstructed spectral function at timestep '//intTochar(k), &
               (/minmass,-10./),(/maxmass,10./),(/0.01,0.2/))
       end do

       initflag = .false.
    end if

    if (mod(timestep,frequency)==0) then
       k = timestep/frequency

       ! Spectral Function
       do i=lbound(Particles,dim=1),ubound(Particles,dim=1)
          do j=lbound(Particles,dim=2),ubound(Particles,dim=2)
             if (Particles(i,j)%ID.ne.id) cycle
             call AddHist(SF(k),abs4(Particles(i,j)%momentum),&
                  Particles(i,j)%perweight*Ninv)
          end do
       end do
       call writeHist(SF(k),10,file='TGP_SF_'//trim(inttochar(k))//'.dat')

       ! Spectral Function versus radius (or z coordinate)
       do i=lbound(Particles,dim=1),ubound(Particles,dim=1)
          do j=lbound(Particles,dim=2),ubound(Particles,dim=2)
             if (Particles(i,j)%ID.ne.id) cycle
             if (plot_z_coord) then
                call AddHist2D(SF_vsRadius(k), &
                     (/abs4(Particles(i,j)%momentum),particles(i,j)%position(3)/), &
                     Particles(i,j)%perweight*Ninv)
             else
                call AddHist2D(SF_vsRadius(k), &
                     (/abs4(Particles(i,j)%momentum),absPos(particles(i,j))/), &
                     Particles(i,j)%perweight*Ninv)
             end if
          end do
       end do
       call writeHist2D_Gnuplot(SF_vsRadius(k), 10, &
            file='TGP_SFvsRadius_'//trim(inttochar(k))//'.dat')
    end if

  end subroutine transportGivenParticle_analyze

end module transportGivenParticleAnalysis
