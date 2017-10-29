
program test_inMedRecipes
  use version, only: PrintVersion
  use output
  use inputGeneral
  use particleProperties, only: initParticleProperties
  use insertion, only: GarbageCollection
  use densityModule
  use coulomb, only: updateCoulomb
  use propagation
  use yukawa
  use particleDefinition
  use nucleusDefinition
  use initNucleus_in_PS,only : initNucPhaseSpace
  use energyCalc
  use neutrino_IDTable 
  use nucleus, only : getTarget

  implicit none  

  type(particle),Allocatable :: realParticles(:,:)

  type(tNucleus),pointer :: targetNuc      
  integer :: lengthReal=0 !maximal number of real particles per ensemble

  integer :: i,j,k,l
  logical :: raiseEnergy

  real,save :: costheta=0.848048096
  real,save :: elepton=0.5
  real,save :: delta_elepton=0.005
  real,save :: enu=.7

  call PrintVersion
  call readInputGeneral
  call initParticleProperties  

  targetNuc => getTarget()                 !set up target resting at 0
  lengthReal = targetNuc%mass ! Real particles fixed by nucleons in target

  write(*,*)'################################################################'

  raiseEnergy=.false.
  energyloop: do j=1,num_energies
     write(*,*) 'Energy loop :',j,num_Energies

     subsequentRunsLoop: Do i=1,num_runs_SameEnergy   ! loop over subsequent runs
        write(*,*) 'Run loop :',i,num_runs_SameEnergy


        Allocate(realparticles(1:numEnsembles,1:lengthReal))

        call GarbageCollection(realParticles)

        call initNucPhaseSpace(realparticles,targetNuc)
        call GarbageCollection(realParticles)

        call updateDensity(realParticles)
!        call updateCoulomb(realParticles,.true.) !Initialization of Coulomb
        call updateCoulomb !Initialization of Coulomb
        call updateYukawa(.true.)
        call updateVelocity(realParticles)
        Do l=1,Size(realParticles,dim=2)
           Do k=1,Size(realParticles,dim=1)
              If (realparticles(k,l)%ID > 0) then
                 call energyDetermination(realParticles(k,l))
              end if
           End do
        end do


        call init_event(realParticles,raiseEnergy,num_runs_sameEnergy)

        raiseEnergy=.false.

        deAllocate(realParticles)

        write(*,*)'################################################################'

     End do subsequentRunsLoop
     raiseEnergy=.true.
  End do energyLoop

  write(*,*) 'all is finished - dumdidum'


contains

  subroutine initInput
    use output, only: Write_ReadingInput
    
    integer :: ios
    
    NAMELIST /inMedQE/ costheta,elepton,delta_elepton,enu

    rewind(5)
    read(5,nml=inMedQE,IOSTAT=ios)
    call Write_ReadingInput("inMedQE",0,ios)
    
    write(*,*) 'enu= ', enu
    write(*,*) 'costheta= ', costheta
    write(*,*) 'Elepton=', elepton
    write(*,*) 'delta_Elepton=', delta_elepton
    
    call Write_ReadingInput('inMedQE',1)

  end subroutine initInput


  subroutine init_event(realParticles,raiseFlag,num_runs_sameEnergy)

    use leptonTensor
    use hadronicTensorQE, only: hadronicTensorforQE
    use particleDefinition
    use idtable, only : nucleon
    use pauliBlockingModule, only : pauliBlocking
    use constants, only: pi, mN
    use minkowski, only : contract,SP
    use neutrino_IDTable
    use leptonicID
    use lorentztrafo
    use NeutrinoMatrixElement

    integer, save :: process_ID=1
    !integer, save :: flavor_ID=2

    type(particle), dimension(:,:),intent(inOut) :: realParticles
    integer, intent(in) :: num_runs_sameEnergy
    logical, intent(in) :: raiseFlag ! if .true. then the energy etc is raised by raiseValue 

    integer :: numtry=0
    integer, save :: numberofcalls=0

    real, dimension(1:10)  ::  sigma=0.
    real, dimension(1:10) :: sigabs=0.
    real, save, dimension(1:10) :: sigabsfinal=0.
    real, dimension(1:10) :: MaEl
    real :: kin_factor

    real ::  mass_out=0.
    integer::  charge=0
    real, dimension(0:3) :: p_out,k_out
    real, dimension(0:3) :: k_in,p_in
    real, dimension(0:3) :: pi_dF,k_in_lor,k_out_lor,p_out_lor,p_in_lor
    
    real :: lorentzfactor,qz
    real :: invMass
    real :: plep
    real :: ml_out
    real :: width
    integer :: i,j
    complex,dimension(0:3,0:3) :: hadronTensor,leptonTens
    logical, save ::initflag=.true.


    logical, save :: q_in_z_dir=.true.
    if(initFlag) then
       call initInput
       initFlag=.false.
       Open(10,File='QE-inmed.dat') 
       write(10,*) '#  raiseFlagVariable, xsection'
       Close(10)
    end if

    write(*,*) 'enu= ', enu
    write(*,*) 'costheta= ', costheta
    if(raiseFlag) then
       elepton=elepton+delta_elepton
       write(*,*) 'Elepton is raised by ...', delta_elepton, ' to ', elepton
    end if

    !neglect electron mass
    ml_out=0.

    plep=sqrt(max((elepton**2-ml_out**2),0.))


    MaEl=0. 

    numberofcalls=numberofcalls+1

    if(numberofcalls.gt.num_runs_sameEnergy) then
       numberofcalls=1
       sigabsfinal=0.
    end if

    !absorption cross section
    sigabs=0.

    !loop to determine numtry (number of testteilchen)
    numtry=0
    loopOverEnsemble0: do i = lBound(realParticles,dim=1),uBound(realParticles,dim=1)
       loopVector0: do j = lBound(realParticles,dim=2),uBound(realParticles,dim=2)
          if(realParticles(i,j)%ID.ne.nucleon) cycle
          numtry=numtry+1
       end do loopVector0
    end do loopOverEnsemble0


    loopOverEnsemble: do i = lBound(realParticles,dim=1),uBound(realParticles,dim=1)
       loopVector: do j = lBound(realParticles,dim=2),uBound(realParticles,dim=2)

          if(realParticles(i,j)%ID.ne.nucleon) cycle

          sigma=0.

          p_in=realParticles(i,j)%momentum 
          charge=realParticles(i,j)%charge

          !set kinematics of leptons 
          !second option: q in z-direction
          if(q_in_z_dir) then !q in z-direction
             k_in(0)=enu 
             k_out(0)=elepton
             qz=sqrt(k_in(0)**2+k_out(0)**2- 2.*k_in(0)*k_out(0)*costheta)
             k_in(3)=enu*(enu-elepton*costheta)/qz
             k_in(1)=sqrt(enu**2-k_in(3)**2)
             k_in(2)=0.
             k_in(3)=enu*(enu-elepton*costheta)/qz
             k_out(1)=k_in(1)
             k_out(2)=0.
             k_out(3)=elepton*(enu*costheta-elepton)/qz

          else !incoming neutrino in z-direction
             k_in=(/enu,0.,0.,enu/) 
             k_out(0)=elepton
             k_out(1)=sqrt(max((1.-costheta**2),0.))*plep
             k_out(2)=0.
             k_out(3)=costheta*plep
             
             write(*,*)'q_in_z_dir=.false. leads to problems with deForest current conservation -> STOP'
             stop
             
          end if

          !momentum of outgoing particle
          p_out=p_in+k_in-k_out

          if(SP(p_out,p_out).le.0.) cycle !reaction not possible

          invMass=sqrt(SP(p_out,p_out))


          !set phasespace: common to all approaches
          !assumption: outgoing nucleon free

          width=0.001
          kin_factor=plep/pi**2/32./k_in(0)/p_in(0)*width/pi*invMass/((invMass**2-mN**2)**2 +invMass**2*width**2)

          mass_out=mN

          leptonTens=leptonicTensor(EM,k_in,k_out)


          
          !different hadronic tensors

          !1: our standard approach with current-conserving term 
          if(hadronicTensorforQE(p_in,p_out,charge,process_ID,hadronTensor,1)) & 
               & MaEl(1)=kin_factor*REAL(Contract(leptonTens,hadronTensor))
          
          !2: our standard approach WITHOUT current-conserving term
          if(hadronicTensorforQE(p_in,p_out,charge,process_ID,hadronTensor,2)) & 
               & MaEl(2)=kin_factor*REAL(Contract(leptonTens,hadronTensor))
      
          !2': deForest current no. 2 (same as ours) WITHOUT current-conserving term
          pi_dF=p_in
          pi_dF(0)=sqrt(mN**2+p_in(1)**2+p_in(2)**2+p_in(3)**2)
          if(hadronicTensorforQE(pi_dF,p_out,charge,process_ID,hadronTensor,2)) & 
               & MaEl(3)=kin_factor*REAL(Contract(leptonTens,hadronTensor))
        
          !3: deForest current no. 2 (same as ours) with current-conserving via longitudinal component (see Tinas notes)
          pi_dF=p_in
          pi_dF(0)=sqrt(mN**2+p_in(1)**2+p_in(2)**2+p_in(3)**2)
          if(hadronicTensorforQE(pi_dF,p_out,charge,process_ID,hadronTensor,3)) & 
               & MaEl(4)=kin_factor*REAL(Contract(leptonTens,hadronTensor))
          
          !4: deForest current no. 1 WITHOUT current-conserving term
          pi_dF=p_in
          pi_dF(0)=sqrt(mN**2+p_in(1)**2+p_in(2)**2+p_in(3)**2)
          if(hadronicTensorforQE(pi_dF,p_out,charge,process_ID,hadronTensor,4)) & 
               & MaEl(5)=kin_factor*REAL(Contract(leptonTens,hadronTensor))

          !5: deForest current no. 1 with current-conserving via longitudinal component (see Tinas notes)
          pi_dF=p_in
          pi_dF(0)=sqrt(mN**2+p_in(1)**2+p_in(2)**2+p_in(3)**2)
          if(hadronicTensorforQE(pi_dF,p_out,charge,process_ID,hadronTensor,5)) & 
               & MaEl(6)=kin_factor*REAL(Contract(leptonTens,hadronTensor))
        
          !6: Saclay 
          if(hadronicTensorforQE(p_in,p_out,charge,process_ID,hadronTensor,6)) & 
               & MaEl(7)=kin_factor*REAL(Contract(leptonTens,hadronTensor))
          

          !for checking: old Mathematica ansatz
          MaEl(8)=kin_factor*REAL(nuMaEl(process_ID,1, charge, k_in, k_out, p_in, p_out,mN))




          !NOW WITH THE GENERAL TENSOR 
          !Lorentztrafos
          k_in_lor=k_in
          k_out_lor=k_out
          p_in_lor=p_in
          p_out_lor=p_out
!!$          call lorentz(realParticles(i,j)%velocity,k_in_lor)
!!$          call lorentz(realParticles(i,j)%velocity,k_out_lor)
!!$          call lorentz(realParticles(i,j)%velocity,p_in_lor)
!!$          call lorentz(realParticles(i,j)%velocity,p_out_lor)

          !write(*,'(10g12.5)')'check lorentz',p_in,p_in_lor

          lorentzfactor=1.    !p_in_lor(0)/p_in(0) !*k_in_lor(0)/k_in(0)
          kin_factor=k_out_lor(0)/pi**2/32./k_in_lor(0)/p_in_lor(0)* & 
               & width/pi*invMass/((invMass**2-mN**2)**2 +invMass**2*width**2)*lorentzfactor

          leptonTens=leptonicTensor(EM,k_in_lor,k_out_lor)

          !Benhar
          pi_dF=p_in_lor
          pi_dF(0)=sqrt(mN**2+p_in_lor(1)**2+p_in_lor(2)**2+p_in_lor(3)**2)
          if(hadronicTensorforQE_genTensor(pi_dF,p_out_lor,charge,process_ID,hadronTensor)) & 
               & MaEl(9)=kin_factor*REAL(Contract(leptonTens,hadronTensor))
          
          !Ferree
          if(hadronicTensorforQE_genTensor(p_in_lor,p_out_lor,charge,process_ID,hadronTensor)) & 
               & MaEl(10)=kin_factor*REAL(Contract(leptonTens,hadronTensor))
         


    
          sigma=MaEl

          sigma=0.197**2*10.*sigma  !cross section dsig/dElep dOmega  in millibarn/GeV/sr

          if(pauliBlocking(p_out,realParticles(i,j)%position,charge,realParticles)) sigma=0.

          sigabs=sigabs+sigma/float(numtry)

       end do loopVector
    end do loopOverEnsemble

    write(*,'(11g12.5)')'sigabs',sigabs,'numberofcalls',numberofcalls

    sigabsfinal=sigabs+sigabsfinal
    write(*,'(11g12.5)') 'sigabsfinal',enu-elepton, sigabsfinal/real(numberofcalls)
    
    Open(10,File='QE-inmed.dat',position='append')
    if(numberofcalls.ne.1) backspace(10)
    write(10,'(11g12.5)') enu-elepton, sigabsfinal/real(numberofcalls)
    close(10)

  end subroutine init_event
  
!#######################################################################################################
  
  function hadronicTensorforQE_genTensor(pi,pf,targetCharge,process,matrix)   result(success)
    use  FF_QE_nucleonScattering
    use leptonicID, only : EM
    use minkowski, only          : abs4,SP,metricTensor

    real, dimension(0:3),intent(in)             :: pi, pf
    integer,intent(in)                          :: process, targetCharge
    complex, dimension(0:3,0:3),intent(inout)   :: matrix
    logical                                     :: success 
    real, dimension(1:2):: G   ! Vector or EM-Form factors
    real, dimension(1:2):: GA  ! Axial form factors
    real                        :: mi,mf
    integer                     :: mu,nu
    real, dimension(0:3)        :: q
    real :: W1,W2,GM,GE,QSquared,tau

    success=.false.
    matrix=0.  

    if(sp(pf,pf).le.0) then
       return
    end if

    call formfactors_QE(-SP(pf-pi,pf-pi),process,targetcharge,G(1),G(2),GA(1),GA(2),GE,GM) 
  
    mi=abs4(pi)
    mf=abs4(pf)

     
    q= pf-pi
    QSquared=-SP(q,q)
    tau=QSquared/4./mi**2
        
    W1=tau*GM**2
    W2=(GE**2+tau*GM**2)/(1.+tau)

    matrix=0.
    do mu=0,3
       do nu=0,3
          matrix(mu,nu)=W1*(-metricTensor(mu,nu)-q(mu)*q(nu)/QSquared)+& 
               & W2/mi**2*(pi(mu)+SP(pi,q)*q(mu)/QSquared)*(pi(nu)+SP(pi,q)*q(nu)/QSquared)
          matrix(mu,nu)=matrix(mu,nu)*4.*mi*mf
       end do
    end do

    success=.true.

  end function hadronicTensorforQE_genTensor
 
!####################################################################


end program Test_inMedRecipes
