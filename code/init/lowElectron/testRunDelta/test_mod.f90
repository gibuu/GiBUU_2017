program main
  ! Test the 5-fold differential cross section by integrating over dOmega(pion). 
  ! The resulting dSigma/dOmega(electron)/dE(electron) can be compared to the 
  ! experimental results of O'Connell et al. (PRL 53(17))
  ! Writes result to file 'dsigma_dOmegadE.dat'.
  use inputGeneral

  call init_Database
  call readInputGeneral
  call test_piN

contains

  subroutine test_piN
    use constants, only :pi
    use electronNuc_to_Delta, only : eN_delta_piN
    use particleDefinition
    use particleProperties
    use Idtable, only : nucleon
    use degRad_conversion
    implicit none

    ! The O'Connell experiment: PRL 53(17)
    real,parameter :: initial_energy=0.73
    real,parameter :: thetaElec_deg=37.1

    integer, parameter :: maxSteps=120               ! Number of photon energies
    integer, parameter :: maxSteps_phi=50          ! Number of integration steps in phi(pion)
    integer, parameter :: maxSteps_dCosTheta=50   ! Number of integration steps in cosTheta(pion)

    real     :: final_energy
    real     :: phiPion_deg,thetaPion_deg
    real     :: cs,total
    integer  :: i,j,kk
    real     :: dcosTheta, dPhi

    real, dimension(0:3) :: q,lf,k,pf
    
    type(particle) :: initNuc

    initNuc%mass=baryon(nucleon)%mass
    initNuc%ID  =nucleon
    initNuc%charge  =1
    initNuc%position=0
    
    initNuc%momentum(1:3)=0.
    initNuc%momentum(0)=freeEnergy(initNuc)
    
    write(*,*) 'p=', initNuc%momentum


    ! Step sizes for dOmega(pion) integration:
    dcosTheta=2./float(maxsteps_dCosTheta)
    dPhi=2.*pi/float(maxsteps_phi)

    write(*,*) 'dPhi=',dphi
    write(*,*) 'dCos(Theta)=',dCosTheta
    write(*,*)
    open(11,file='dsigmadOmegadE.dat')

    thetaPion_deg=10
    phiPion_deg=60
    final_energy=0.3

    total=total+eN_delta_piN(1,initial_energy,final_energy,thetaElec_deg,thetaPion_deg,phiPion_deg,initNuc,q,lf,k,pf)
    stop

    energyLoop: do i=1,maxsteps-1
       final_Energy=initial_energy/float(maxsteps)*float(i)

       ! Integrate over phi and cos(theta) of the outgoing pion
       total=0.
       phiLoop: do j=1,maxsteps_phi
          phiPion_deg=degrees((float(j)-0.5)*dphi)
          cosThetaLoop: do kk=1,maxsteps_dCostheta
             thetaPion_deg=degrees(acos(-1.+(float(kk)-0.5)*dcosTheta))
             ! Sum both proton contributions:

             !  real function eN_delta_piN(n,e,ef,thf,thpi,phipi,initNuc,q,kf,ppi,pf) 
             total=total+eN_delta_piN(1,initial_energy,final_energy,thetaElec_deg,thetaPion_deg,phiPion_deg,initNuc,q,lf,k,pf)
             total=total+eN_delta_piN(2,initial_energy,final_energy,thetaElec_deg,thetaPion_deg,phiPion_deg,initNuc,q,lf,k,pf)
          end do cosThetaLoop
       end do phiLoop
       total=total*dphi*dcosTheta
       write(11,'(6E15.3)')initial_energy-final_energy,  total
    end do energyLoop
    close(11)
  end subroutine test_piN
end program main
