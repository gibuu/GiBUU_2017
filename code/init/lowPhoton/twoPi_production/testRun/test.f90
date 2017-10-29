program test
  use inputGeneral
  use particleProperties, only: initParticleProperties
  implicit none

  call initParticleProperties
  call readInputGeneral
  call testXsection

contains

  subroutine testXsection
    use gamma2Pi_Xsections
    use Idtable
    use mediumDefinition
    use output
    use particleDefinition
    use particleProperties, only: hadron, nDecays
    use resProd_lepton, only: sigma_pipi_res_vac, sigma_barmes_res_vac
    use ClebschGordan
    use photonXSections, only: calcXS_gammaN2VN
    use constants, only: mN

    real, dimension(0:3) :: sig2pi, sigRes_2Pi, sig2pi_VM
    real, dimension(1:4) :: sigVM
    real, dimension(1:3) :: sigRes_VM
    real, dimension(0:3) :: photon_mom
    integer :: i,ios, j, ID ,k
    real :: srts,photonEnergy,BR_2pi,IS
    real, save :: density=0.168
    real, dimension(1:3) :: betaToLRF,position
    type(medium) ::  mediumAtPosition
    type(particle) :: pro,neu

    NAMELIST /test/ density

    call Write_ReadingInput('test',0)

    rewind(5)
    read(5,nml=test,IOSTAT=IOS)
    call Write_ReadingInput('test',0,IOS)

    write(*,*) 'Set density      :' ,density
    call Write_ReadingInput('test',1)

    Open(11,File='gammaProton.dat')
    Open(10,File='gammaNeutron.dat')

    write(10,*) '# E_gamma,srts,sig2pi,sigRes_2Pi'
    write(10,*)
    write(10,*) '# * sig2pi(0:3) Cross sections for gamma+nucleon->nucleon+2Pi production'
    write(10,*) '# * sig2pi(1) -> nucleon piMinus piPlus'
    write(10,*) '# * sig2pi(2) -> nucleon piPlus piNull or nucleon piMinus piNull'
    write(10,*) '# * sig2pi(3) ->nucleon piNull piNull'
    write(10,*) '# * sig2pi(0) Total Xsection into nucleon+2 Pions'

    write(11,*) '# E_gamma,srts,sig2pi,sigRes_2Pi'
    write(11,*)
    write(11,*) '# * sig2pi(0:3) Cross sections for gamma+nucleon->nucleon+2Pi production'
    write(11,*) '# * sig2pi(1) -> nucleon piMinus piPlus'
    write(11,*) '# * sig2pi(2) -> nucleon piPlus piNull or nucleon piMinus piNull'
    write(11,*) '# * sig2pi(3) ->nucleon piNull piNull'
    write(11,*) '# * sig2pi(0) Total Xsection into nucleon+2 Pions'


    position=0.
    mediumAtPosition%useMedium=.false.
    mediumAtPosition%densityNeutron=density/2.
    mediumAtPosition%densityProton=density/2.

    pro%ID = nucleon
    pro%charge = 1
    pro%mass = mN
    pro%momentum(0) = mN
    neu%ID = nucleon
    neu%charge = 0
    neu%mass = mN
    neu%momentum(0) = mN

    do i=1,150
       photonEnergy = 0.3 + i*0.01
       srts=sqrt((photonEnergy+mN)**2-photonEnergy**2)
       betaToLRF=-(/photonEnergy,0.,0./)/(photonEnergy+mN)
       photon_mom = (/photonEnergy,0.,0.,photonEnergy/)

       ! calculate vector meson contribution
       call calcXS_gammaN2VN (srts, mediumAtPosition, sigVM(1:4))
       ! Subtract resonance contribution from V X channels (only important for rho)
       sigRes_VM(1:3) = sigma_barmes_res_vac(pro,photon_mom,nucleon,(/rho,omegaMeson,phi/))*1000.
       do j=1,3
          sigVM(j)=max(0.,sigVM(j)-sigRes_VM(j))
       end do
       ! Calculate vector meson contribution
       sig2Pi_VM = 0.
       do j=1,3                                ! loop over vector mesons
          ID = 101+2*j                         ! ID = 103,105,107
          BR_2pi = 0.
          do k=1,nDecays
             if (hadron(ID)%decaysID(k)==1) BR_2pi = hadron(ID)%decays(k)        ! branching ratio (V->2pi) at pole mass
          end do
          IS = float(hadron(ID)%isoSpinTimes2)/2.   ! isospin of vector meson
          ! total :
          sig2Pi_VM(0)=sig2Pi_VM(0)+sigVM(j)*BR_2pi
          ! pi+ pi- :
          sig2Pi_VM(1)=sig2Pi_VM(1)+sigVM(j)*BR_2pi*(ClebschSquared(1.,1.,IS,1.,-1.)+ClebschSquared(1.,1.,IS,-1.,1.))
          ! pi+ pi0 and pi- pi0 have no contributions from neutral vector mesons
          ! pi0 pi0 :
          sig2Pi_VM(3)=sig2Pi_VM(3)+sigVM(j)*BR_2pi*ClebschSquared(1.,1.,IS,0.,0.)
       end do

       !!! PROTON
       call gamma2pi(1,srts,sig2pi,betaToLRF,mediumAtPosition,position)
       sigRes_2Pi = sigma_pipi_res_vac(pro,photon_mom)*1000.
       write(11,'(14F9.3)') photonEnergy,srts,sig2pi,sigRes_2Pi,sig2Pi_VM
       write(*,'(14F9.3)') photonEnergy,srts,sig2pi,sigRes_2Pi,sig2Pi_VM

       !!! NEUTRON
       call gamma2pi(0,srts,sig2pi,betaToLRF,mediumAtPosition,position)
       sigRes_2Pi = sigma_pipi_res_vac(neu,photon_mom)*1000.
       write(10,'(10F9.3)') photonEnergy,srts,sig2pi,sigRes_2Pi
       write(*,'(10F9.3)') photonEnergy,srts,sig2pi,sigRes_2Pi
    end do
  end subroutine testXsection
end program test
