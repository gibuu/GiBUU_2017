!******************************************************************************
!****m* /FragmentNucleons
! NAME
! module FragmentNucleons
!
! PURPOSE
! This does a fragmentation according SMM and deletes bound nucleons and
! inserts nucleons stemming from this fragmentation.
!******************************************************************************
module FragmentNucleons

  implicit none
  private

  !****************************************************************************
  !****g*
  ! SOURCE
  !
  integer, save :: BoundMethod = 1
  !
  ! PURPOSE
  ! select method, how one defines "is bound":
  ! * =1: E_kin < 0, i.e. momentum(0) < mass0
  ! * =2: rho < rho0/100
  !****************************************************************************

  !****************************************************************************
  !****g*
  ! SOURCE
  !
  logical, save :: DoRemoveBound = .true.
  !
  ! PURPOSE
  !****************************************************************************

  public :: DoAddFragmentNucleons

contains

  !****************************************************************************
  !****s* FragmentNucleons/DoAddFragmentNucleons
  ! NAME
  ! subroutine DoAddFragmentNucleons(Part,target)
  ! INPUTS
  ! ...
  ! OUTPUT
  ! ...
  ! NOTES
  ! * we assume, that the run was in parallel mode. (This is no
  !   restriction, since the real part vector is only modified in
  !   real-real collisions, and here we have no 'local/full' mode yet.)
  ! * At the moment, we assume only ONE single source. (cf.
  !   module sourceAnalysis how one should extent this.)
  ! * At the moment, RMF is not implemented yet.
  !****************************************************************************
  subroutine DoAddFragmentNucleons(Part,targetNuc)
    use particleDefinition
    use sourceTypeDefinition
    use Insertion, only: GarbageCollection
    use RMF, only: getRMF_flag
    use lorentzTrafo, only: lorentz
    use potentialModule, only: trueEnergy
    use nucleusDefinition
    use EventInfo_HiLep, only: EventInfo_HiLep_Get
    use constants, only: mN
    use CallStack, only: TRACEBACK

    type(particle), dimension(:,:), intent(inOut) :: Part
    type(tNucleus),pointer :: targetNuc

    integer :: iEns, iPart
    type(quelle) :: Source
    type(quelle) :: Source0
    real, dimension(0:3) :: SourceMom, SourceMomStar, mom
    real :: Eout_part, Eout_all
    real :: BEin,BEout
    real :: Eexc,EexcB
    type(particle) :: PartSRF
    real :: weight,nu,Q2,epsilon
    integer :: EventType


    Source0%position = 0.0


    write(*,'("B(N=",i3,",Z=",i3,") = ",f12.4,"MeV")') &
         & targetNuc%mass-targetNuc%charge,targetNuc%charge,&
         & 1000*BWformula( (/real(targetNuc%mass-targetNuc%charge),real(targetNuc%charge)/) )/targetNuc%mass
    write(*,'("B(N=",i3,",Z=",i3,") = ",f12.4,"MeV")') &
         & targetNuc%mass-targetNuc%charge-1,targetNuc%charge,&
         & 1000*BWformula( (/real(targetNuc%mass-targetNuc%charge-1),real(targetNuc%charge)/) )/(targetNuc%mass-1)
    write(*,'("B(N=",i3,",Z=",i3,") = ",f12.4,"MeV")') &
         & targetNuc%mass-targetNuc%charge,targetNuc%charge-1,&
         & 1000*BWformula( (/real(targetNuc%mass-targetNuc%charge),real(targetNuc%charge-1)/) )/(targetNuc%mass-1)


    do iEns = 1,size(Part,dim=1)
       Source = Source0
       SourceMom = 0
       SourceMomStar = 0

       ! Step 1: determine LRF of Source

       Eout_all = 0.0
       do iPart=1,size(Part,dim=2)
          if (Part(iEns,iPart)%Id <  0) exit
          if (Part(iEns,iPart)%Id == 0) cycle

          Eout_all = Eout_all + Part(iEns,iPart)%momentum(0)

          if (Part(iEns,iPart)%Id.ne.1) cycle ! no nucleon
          if (Part(iEns,iPart)%antiparticle) cycle
          if (.not.IsBound()) cycle

          Source%Size   = Source%Size+1
          Source%Charge = Source%Charge+Part(iEns,iPart)%Charge
          Source%position = Source%position + Part(iEns,iPart)%position
          SourceMomStar = SourceMomStar + Part(iEns,iPart)%momentum
       end do
       Eout_part = Eout_all - SourceMomStar(0)

       if (Source%Size==0) cycle

       Source%velocity = SourceMomStar(1:3) / SourceMomStar(0)
       Source%position = Source%position/Source%Size

       BEin = BWformula( (/float(targetNuc%mass-targetNuc%charge), float(targetNuc%charge)/) )
       BEout= BWformula( (/float(Source%Size-Source%Charge), float(Source%Charge)/) )

       BEin = BEin/targetNuc%mass
       BEout = BEout/Source%Size

       Eexc = Eout_all-Eout_part-(mN-BEout)*Source%Size &
            & - (mN-BEout)*Source%Size*Dot_product(Source%velocity,Source%velocity)/2

       if (EventInfo_HiLep_Get(iEns,-1,Weight,nu,Q2,epsilon,EventType)) then
          write(*,'(A,f12.4)') 'nu       = ', nu
       else
          write(*,*) 'ooops.'
          stop
       end if

       EexcB = nu+(mN-BEin)*targetNuc%mass-Eout_part-(mN-BEout)*Source%Size &
            & - (mN-BEout)*Source%Size*Dot_product(Source%velocity,Source%velocity)/2

       write(*,'(A,f12.4)') 'Eout_all = ', Eout_all
       write(*,'(A,f12.4)') 'Eout_part= ', Eout_part
       write(*,'(A,f12.4)') '(m-B)A   = ', (mN-BEin)*targetNuc%mass
       write(*,'(A,f12.4)') "(m-B')A' = ", (mN-BEout)*Source%Size
       write(*,'(A,f12.4)') "E_beam   = ", Eout_all-(mN-BEin)*targetNuc%mass
       write(*,'(A,2f12.4)') "E_exc    = ", Eexc,Eexc/Source%Size
       write(*,'(A,2f12.4)') "E_exc (B)= ", EexcB,EexcB/Source%Size


       ! Step 2: determine total Energy after boosting into LRF of Source
       !         (=SRF, Source Rest Frame)
       !         (only done in Skyrme mode)

       if ( .not. getRMF_flag() ) then

          do iPart=1,size(Part,dim=2)
             if (Part(iEns,iPart)%Id <  0) exit
             if (Part(iEns,iPart)%Id.ne.1) cycle ! no nucleon
             if (Part(iEns,iPart)%antiparticle) cycle
             if (.not.IsBound()) cycle

             PartSRF = Part(iEns,iPart)
             mom(0:3) = Part(iEns,iPart)%momentum(0:3)
             call lorentz(Source%velocity,mom) ! boost to SRF
             PartSRF%momentum = mom

             SourceMom(0) = SourceMom(0) + trueEnergy(PartSRF)
          end do
          SourceMom(1:3) = 0.0

       else
          call TRACEBACK('RMF not yet implemented')
       end if

       ! Step 3: Remove bound nucleons from particle vector, if desired

       if (DoRemoveBound) then
          do iPart=1,size(Part,dim=2)
             if (Part(iEns,iPart)%Id <  0) exit
             if (Part(iEns,iPart)%Id.ne.1) cycle ! no nucleon
             if (Part(iEns,iPart)%antiparticle) cycle
             if (.not.IsBound()) cycle

             Part(iEns,iPart)%Id = 0 ! delete particle
          end do

          call GarbageCollection(Part, iiEns = iEns)
       end if

       ! Step 4:

       Source%status=.true.
       Source%ExEnergy = SourceMom(0)
       Source%radEnergy= 0.0

       write(*,'(A,2i4)') 'Source:',Source%Size, Source%Charge
       write(*,'(A,4es12.3)') '  pos:  ',0.0,Source%position
       write(*,'(A,4es12.3)') '  p*:   ',SourceMomStar
       write(*,'(A,4es12.3)') '  vel:  ',&
            & 1./sqrt(1.-Dot_product(Source%velocity,Source%velocity)),Source%velocity
       write(*,'("==>",es12.3)') SourceMom(0)/Source%Size

    end do

  contains

    logical function IsBound()
      use constants, only: rhoNull
      use densitymodule, only: DensityAt
      use dichteDefinition, only: dichte

      real, parameter :: densCut = rhoNull/100
      type(dichte)    :: LocalDens

      IsBound = .false.

      select case (BoundMethod)
      case (1)
         IsBound = (Part(iEns,iPart)%momentum(0) < Part(iEns,iPart)%mass)
      case (2)
         LocalDens = DensityAt(Part(iEns,iPart)%position)
         IsBound = (LocalDens%baryon(0) > densCut)
      case default
         call TRACEBACK('wrong BoundMethod')
      end select

    end function IsBound


    real function BWformula(A)
      real, dimension(1:2), intent(in) :: A
      real, parameter :: av=15.67, a0=17.23, ac=0.714, as=23.3 !MeV
      real :: M,N,Z

      M = A(1)+A(2)
      N = A(1)
      Z = A(2)

      BWformula = av*M - a0*M**(2./3.) - ac*Z*(Z-1.)*M**(-0.33333333) - &
           & as*(N-Z)**2/M

      BWformula = abs(BWformula/1000.)

    end function BWformula

  end subroutine DoAddFragmentNucleons

end module FragmentNucleons
