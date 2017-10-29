!***************************************************************************
!****m* /gem
! PURPOSE
! Fragment De-Excitation using the Generalized Evaporation Model (GEM)
! NOTES
! * Related literature:
! * S. Furihata, Nuclear Instruments and Methods in Physics Research, 
! * B171 (2000) 251-258.
! * S.G. Mashnik et al., nucl-th/0208048
! * GEM module: 4-momenta & masses in units of [MeV]!!!
!***************************************************************************
module GEM_Evaporation

  implicit none

  PRIVATE

  real, save :: pi = 3.14159

  PUBLIC :: StabilityGEM

  contains 


    !***********************************************************************
    !****s* GEM_Evaporation/StabilityGEM
    ! NAME
    ! subroutine StabilityGEM
    ! 
    ! PURPOSE
    ! Main routine of GEM Evaporation. 
    ! NOTES
    ! * Calling init for Reading of data-tables from NNDC (mass excess, spin)
    ! * Calling of important functions: 
    ! * Calculation of Decay probabilities & Monte-Carlo decision
    ! * Calculation of Energies & Momenta of daugther & ejectile
    ! * Re-Construction of the new FragmentVector
    !***********************************************************************
    subroutine StabilityGEM(iGEM,Particles,Clusters,npar,ifrm,& 
         &                  FragmentVector,ParticleVector, & 
         &                  newFragments)
      use ClusterEnergy,   only : energy,boost
      use typeDefinitions, only : Cluster,Particle
      use InitGEM,         only : Init_GEM,Nuclide,Elements,NucPosition
      use WriteStatus,     only : IOControl

      !---------------------------------------------------------------------
      !Input-Variables
      !---------------------------------------------------------------------
      integer, intent(in) :: Particles,iGEM  !,Clusters
      !---------------------------------------------------------------------
      !Input-Output Variables
      !---------------------------------------------------------------------
      type(particle), dimension(:),  intent(inout) :: ParticleVector
      type(cluster),  dimension(:),  intent(inout) :: FragmentVector
      integer,        dimension(:,:),intent(inout) :: npar
      integer,        dimension(:),  intent(inout) :: ifrm
      integer,                       intent(inout) :: Clusters
      integer,                       intent(out)   :: newFragments
      !---------------------------------------------------------------------
      !Local Variables
      !---------------------------------------------------------------------
      integer, allocatable, dimension(:)     :: nreal,PosFra
      integer, allocatable, dimension(:,:)   :: npar_new
      integer, dimension(1:Clusters,1:2)     :: NewA,NewZ
      real,    dimension(0:3)                :: Impulsj,Impulsd
      real,    dimension(1:Clusters,1:2,0:3) :: Mom
      logical, dimension(1:Particles)        :: RepeatGEM

      real,    parameter :: NucMassFm = 4.7534587 !Nuclear mass in [1/fm]

      integer :: status,i,FMass,ZMass,AFinalj,ZFinalj,AFinald,ZFinald,k,m
      integer :: nfnew,Index
      real    :: BE,Etot,Ex,betax,betay,betaz
      real    :: xnew(1:3)
      logical :: DecayFlag
      !---------------------------------------------------------------------
      ! GEM LOOP, MaxGEM-times
      !---------------------------------------------------------------------
      RepeatGEM = .false. !indicates that an unstable cluster has passed 
                          !through the evaporation procedure (true)
      !---------------------------------------------------------------------
      ! allocate & initialize local fields
      !---------------------------------------------------------------------
      allocate(nreal(1:Clusters),STAT=status)
      call IOControl(1,status,'StabilityGEM','nreal')
      allocate(npar_new(1:Particles,1:Particles),STAT=status)
      call IOControl(1,status,'StabilityGEM','npar_new')

      do i=1,Clusters
         nreal(i) = 999
      end do
      NewA      = 0
      NewZ      = 0
      Mom       = 0.
      npar_new  = 0
      !---------------------------------------------------------------------
      ! -> check weather clusters are stable
      ! -> Hot clusters can decay according Furihata's Evaporation Model
      ! -> repeat evaporation procedure until all clusters are stable
      !---------------------------------------------------------------------
      HotCluster : do i=1,Clusters
         FMass = FragmentVector(i)%MassNumber
         ZMass = FragmentVector(i)%ChargeNumber
         if (FMass .le. 3) then !asume very light clusters stable
            nreal(i) = -1
            cycle
         endif
         if (RepeatGEM(i)) then
            nreal(i) = -1
            cycle
         endif
         allocate(PosFra(1:FMass),STAT=status)
         call IOControl(1,status,'StabilityGEM','PosFra')
         PosFra(:) = npar(i,:)
         !---------------------------------------------------------------------
         !* returns average binding energy (exp. values), total energy 
         !  and Excitation energy per nucleon
         !---------------------------------------------------------------------
         CALL ENERGY(FMass,PosFra,ParticleVector,FragmentVector(i)%momentum, & 
              &      Etot,betax,betay,betaz)
         !---------------------------------------------------------------------
         !get the TOTAL excitation energy of the cluster-->input for GEM 
         !---------------------------------------------------------------------
         Index = NucPosition((FMass-ZMass),ZMass)
         if (Index < 0) then
            write(*,*) 'module GEM-Evaporation,routine StabilityGEM:'
            write(*,*) 'wrong element index!!! ',i,FMass,ZMass,Index
            write(*,*) '!!! TERMINATION NOW !!!'
            STOP
         endif
         BE = Elements(Index)%Bind !RMF binding energy/A (MeV)
         Ex = Etot - BE            !Excitation energy (MeV)
!         write(*,222) FMass,ZMass,Etot,BE,Ex
!222      format('***',2i5,3f12.4)
!         write(*,*) 'ENERGY      = ',Etot,BE,Ex,FMass,ZMass,betax,betay,betaz
!         write(*,*)
         !---------------------------------------------------------------------
         deallocate(PosFra,STAT=status)
         call IOControl(2,status,'StabilityGEM','PosFra')
         if (Ex.le.(0.0)) then
            nreal(i) = -1 !--> fragment stable
            cycle
         endif
         nreal(i) = 2 !fragment is excited, apply evaporation procedure:
         call GetDecayChannel(FMass,ZMass,Ex,betax,betay,betaz,& !--> Input
              &               AFinalj,ZFinalj,Impulsj,&     !-->Output1
              &               AFinald,ZFinald,Impulsd,DecayFlag)  !-->Output2
         if (.not.DecayFlag) then!emission patern must exceed Q-value & pairing energy!
            nreal(i) = -1        !otherwise cannot decay, set residual to stable
         else
            nreal(i)   = 2       !sucessful decay
            NewA(i,1)  = Afinalj
            NewZ(i,1)  = ZFinalj
            Mom(i,1,:) = Impulsj(:)
            NewA(i,2)  = Afinald
            NewZ(i,2)  = ZFinald
            Mom(i,2,:) = Impulsd(:)
         endif
      end do HotCluster
      !---------------------------------------------------------------------
      ! Construct the new FragmentVector
      ! Just to remember the units of the FragmentVector:
      ! !!! Mass in [1/fm], Momenta in [GeV], positions in [fm] !!!!!!!!!!!!
      !---------------------------------------------------------------------
      nfnew = 0 !new clusters
      StableClusters : do i=1,Clusters
         if (nreal(i)==999) then
            write(*,*) 'from StabilityGEM: something wrong in StabilityGEM...'
            write(*,*) 'i,nreal()',i,nreal(i)
            STOP
         endif
         if (nreal(i)==-1) then !no changes in FragmentVector
            nfnew                 = nfnew + 1
            RepeatGEM(nfnew)      = .true.
            FragmentVector(nfnew) = FragmentVector(i) 
            do m=1,FragmentVector(nfnew)%MassNumber
               npar_new(nfnew,m)       = npar(i,m)
            end do
         endif
      end do StableClusters

      newFragments = 0

      NewClusters : do i=1,Clusters
         if (nreal(i) .ne. 2) cycle
         if (NewA(i,1)==0 .or. NewA(i,2)==0) then
            write(*,*) 'NewA=0??? something wrong in resorting loop...'
            stop
         endif
         k = 1
         Patern1 : if (NewA(i,k)==1) then !single free particle emission, no cluster
            ifrm(npar(i,NewA(i,k))) = 0
            ParticleVector(npar(i,NewA(i,k)))%Momentum(:) = Mom(i,k,:)*0.001
         else!cluster emission, select bounded particles of patern2
            nfnew = nfnew + 1
            newFragments = newFragments + 1
            RepeatGEM(nfnew) = .false.
            do m=1,NewA(i,k)
               npar_new(nfnew,m)       = npar(i,m)
            end do
            FragmentVector(nfnew)%momentum(:) = Mom(i,k,:)*0.001
            FragmentVector(nfnew)%mass        = real(NewA(i,k))*NucMassFm
            FragmentVector(nfnew)%ID          = 1
            FragmentVector(nfnew)%MassNumber  = NewA(i,k)
            FragmentVector(nfnew)%ChargeNumber= NewZ(i,k)
            FragmentVector(nfnew)%HypNumber   = 0
            FragmentVector(nfnew)%FreeBound   = .true.
            xnew(:) = 0.0
            do m=1,NewA(i,k)
               xnew(:) = xnew(:) + & 
                    & ParticleVector(npar_new(nfnew,m))%position(:) * & 
                    & ParticleVector(npar_new(nfnew,m))%Mass
            end do
            FragmentVector(nfnew)%position(:) =  xnew(:)/FragmentVector(nfnew)%mass
         endif Patern1
         k = k + 1
         Patern2 : if (NewA(i,k)==1) then !single free particle emission, no cluster
            ifrm(npar(i,(NewA(i,k)+NewA(i,k-1)))) = 0 
            ParticleVector(npar(i,(NewA(i,k)+NewA(i,k-1))))%Momentum(:) = Mom(i,k,:)*0.001
         else !cluster emission, select bounded particles of patern2
            nfnew = nfnew + 1
            newFragments = newFragments + 1
            RepeatGEM(nfnew) = .false.
            do m=1,NewA(i,k)
               npar_new(nfnew,m)       = npar(i,m+NewA(i,k-1))
            end do
            FragmentVector(nfnew)%momentum(:) = Mom(i,k,:)*0.001
            FragmentVector(nfnew)%mass        = real(NewA(i,k))*NucMassFm
            FragmentVector(nfnew)%ID          = 1
            FragmentVector(nfnew)%MassNumber  = NewA(i,k)
            FragmentVector(nfnew)%ChargeNumber= NewZ(i,k)
            FragmentVector(nfnew)%HypNumber   = 0
            FragmentVector(nfnew)%FreeBound   = .true.
            xnew(:) = 0.0
            do m=1,NewA(i,k)
               xnew(:) = xnew(:) + & 
                    & ParticleVector(npar_new(nfnew,m))%position(:) * & 
                    & ParticleVector(npar_new(nfnew,m))%Mass
            end do
            FragmentVector(nfnew)%position(:) =  xnew(:)/FragmentVector(nfnew)%mass
         endif Patern2
      end do NewClusters
      !---------------------------------------------------------------------
      do i=1,nfnew
         do m=1,FragmentVector(i)%MassNumber
            npar(i,m) = npar_new(i,m)
         end do
      end do
      !---------------------------------------------------------------------
      deallocate(nreal,STAT=status)
      call IOControl(2,status,'StabilityGEM','nreal')
      deallocate(npar_new,STAT=status)
      call IOControl(2,status,'StabilityGEM','npar_new')
      Clusters = nfnew

      write(*,*) 'from GEM-Loop: ',iGEM,Clusters,newFragments

    !***********************************************************************
    end subroutine StabilityGEM !*******************************************
    !***********************************************************************


    !***********************************************************************
    !****s* GEM_Evaporation/GetDecayChannel
    ! NAME
    ! subroutine GetDecayChannel
    ! 
    ! PURPOSE
    ! * Main routine of finding decay channel
    ! NOTES
    ! * This routine is called for each "hot" cluster in StabilityGEM
    !***********************************************************************
    subroutine GetDecayChannel(FMass,ZMass,Ex,betax,betay,betaz,& !--> Input
              &               AFinalj,ZFinalj,Impulsj,&      !-->Output1
              &               AFinald,ZFinald,Impulsd,DecayFlag)!-->Output2
      use initGEM, only     : Nuclide, Elements,NucPosition
      use WriteStatus, only : IOControl
      use random, only: rn
      !---------------------------------------------------------------------
      !Input-Variables
      !---------------------------------------------------------------------
      integer, intent(in) :: FMass,ZMass
      real,    intent(in) :: Ex,betax,betay,betaz
      !---------------------------------------------------------------------
      !Input-Output-Variables
      !---------------------------------------------------------------------
      integer,                 intent(out) :: AFinalj,ZFinalj,AFinald,ZFinald
      real,    dimension(0:3), intent(out) :: Impulsj,Impulsd
      logical,                 intent(out) :: DecayFlag
      !---------------------------------------------------------------------
      !Local-Variables
      !---------------------------------------------------------------------
      integer :: A,N,Index,Ai,Zi,k,IndexEjectile,IndexDaughter, & 
           &     Aj,Zj,Ad,Zd,Status,Finalj,Finald,FinalChannel
      real    :: MassExAi,MassExAj,SpinAj,MassExAd,ProbTot, &  !rand, & 
           &     x,chooseChannel,Normalization,Width_k
      real, allocatable, dimension(:) :: Width
      !---------------------------------------------------------------------
      DecayFlag= .true.
      A        = FMass
      N        = FMass-ZMass
      Index    = NucPosition(N,ZMass) !Element with label "Index"
      Ai       = Elements(Index)%ZahlN+Elements(Index)%ZahlZ
      Zi       = Elements(Index)%ZahlZ
      if (A.ne.Ai .and. ZMass.ne.Zi) then
         write(*,*) 'from Module GEM_Evaporation:'
         write(*,*) 'Residual cluster does not exist !!!'
         write(*,*) 'something wrong in input or in Asymmetry routine...'
         write(*,*) 'BUU Cluster: FMass,ZMass',FMass,ZMass
         write(*,*) 'Element in NNDC-Table: Index,Ai,Zi',Index,Ai,Zi
         write(*,*) '!!! TERMINATION OF PROGRAM NOW !!!'
         STOP
      endif
      MassExAi = Elements(Index)%MassEx
      !---------------------------------------------------------------------
      ! Determine probability for each decay channel
      !---------------------------------------------------------------------
      allocate(Width(1:size(Elements(Index)%Ejectile)),STAT=status)
      call IOControl(1,status,'GetDecayChannel','Width')
      Width   = 0.0
      ProbTot = 0.0
      Channels : do k=1,size(Elements(Index)%Ejectile)
         if (.not.Elements(Index)%ChannelStatus(k)) exit !invalid channel
         IndexEjectile = Elements(Index)%Ejectile(k)
         IndexDaughter = Elements(Index)%Daughter(k)
         Aj       = Elements(IndexEjectile)%ZahlN+Elements(IndexEjectile)%ZahlZ
         Zj       = Elements(IndexEjectile)%ZahlZ
         MassExAj = Elements(IndexEjectile)%MassEx
         SpinAj   = Elements(IndexEjectile)%Spin
         Ad       = Elements(IndexDaughter)%ZahlN+Elements(IndexDaughter)%ZahlZ
         Zd       = Elements(IndexDaughter)%ZahlZ
         MassExAd = Elements(IndexDaughter)%MassEx
         call DecayWidth(Ai,Zi,Aj,Zj,Ad,Zd,SpinAj,Ex,Width_k)
         Width(k)= Width_k
         ProbTot = ProbTot + Width(k)
      end do Channels
      !---------------------------------------------------------------------
      if (ProbTot==0.) then
         DecayFlag = .false.
         return
      endif
      !---------------------------------------------------------------------
      ! Decide final emission channel by monte-carlo
      !---------------------------------------------------------------------
      x             = rn()*ProbTot
      chooseChannel = 0.0
      FinalChannel  = 0
      do k=1,size(Elements(Index)%Ejectile)
         if (.not.Elements(Index)%ChannelStatus(k)) exit !invalid channel
         chooseChannel = choosechannel + Width(k)
         if (chooseChannel.gt.x) then
            FinalChannel = k
            exit
         endif
      end do
      if (FinalChannel == 0) then
         write(*,*) 'from module GEM_Evaporation, routine GetDecayChannel:'
         write(*,*) 'Final channel not found in monte-carlo decision'
         write(*,*) '!!! TERMINATION OF PROGRAM !!!'
         write(*,*) Index,Ai,Zi
         write(*,*) size(Elements(Index)%Ejectile)
         write(*,*) ProbTot
         write(*,*) Width
         STOP
      endif

      Normalization = Width(FinalChannel)
      if (Normalization < 1.e-06) then
!         write(*,*) 'this should not happened...',Normalization
         DecayFlag = .false.
         return
      endif

      !---------------------------------------------------------------------
      ! Get ejectile and daughter properties of final channel
      !---------------------------------------------------------------------
      Finalj  = Elements(Index)%Ejectile(FinalChannel)
      Finald  = Elements(Index)%Daughter(FinalChannel)
      AFinalj = Elements(Finalj)%ZahlN+Elements(Finalj)%ZahlZ
      ZFinalj = Elements(Finalj)%ZahlZ
      AFinald = Elements(Finald)%ZahlN+Elements(Finald)%ZahlZ
      ZFinald = Elements(Finald)%ZahlZ

      call Get4Momentum(Ai,Zi,Afinald,Zfinald,AFinalj,ZFinalj,& 
           &            Ex,Normalization, & 
           &            betax,betay,betaz,Impulsj,Impulsd)
      !---------------------------------------------------------------------
      deallocate(Width,STAT=status)
      call IOControl(2,status,'GetDecayChannel','Width')
    !***********************************************************************
    end subroutine GetDecayChannel !****************************************
    !***********************************************************************


    !***********************************************************************
    !****s* GEM_Evaporation/DecayWidth
    ! NAME
    ! subroutine DecayWidth
    ! 
    ! PURPOSE
    ! Main routine of calculating decay probability "Width_k"
    !***********************************************************************
    subroutine DecayWidth(Ai,Zi,Aj,Zj,Ad,Zd,SpinAj,Ex,Width_k)
      use initGEM, only : Nuclide,Elements,NucPosition
      !---------------------------------------------------------------------
      !Input-Variables
      !---------------------------------------------------------------------
      integer, intent(in) :: Ai,Zi,Aj,Zj,Ad,Zd
      real,    intent(in) :: SpinAj,Ex
      !---------------------------------------------------------------------
      !Input-Output-Variables
      !---------------------------------------------------------------------
      real, intent(out) :: Width_k
      !---------------------------------------------------------------------
      !Local-Variables & external functions
      !---------------------------------------------------------------------
      integer :: Ni,Nj,Nd,IDi,IDj,IDd
      real    :: U,EjectileMass,DaughterMass,a,delta0,Alpha,Beta,Coulomb
      real    :: Ux,T,E0,Qi,Qj,Qd,MaximalKineticEnergy,tm,tx,s,sx,g
      real    :: ConstantFactor,R1,R2,Rb,GeometricalXS,Exx
      real    :: InitialLevelDensity
      !---------------------------------------------------------------------
      width_k = 0.0
      !---------------------------------------------------------------------
      ! Copy input parameters into local variables
      !---------------------------------------------------------------------
      Ni          = Ai - Zi !NeutronNumber
      Nj          = Aj - Zj !NeutronNumber
      Nd          = Ad - Zd !NeutronNumber
      U           = Ex !Excitation energy of fragment with mass number Ai [MeV]
      EjectileMass= real(Aj)*939. !fragment mass of ejectile [MeV]
      DaughterMass= real(Ad)*939. !fragment mass of daughter [MeV]
      a           = GetLevelDensityParameter(Ad,Zd,U) !Level density parameter a [1/MeV]
!      a           = float(Ad)/8. !1/MeV
      delta0      = GetPairingCorrection(Ad,Zd) !pairing energy [MeV]
      Alpha       = GetAlphaParameter(Aj,Zj,Ad) !alpha parameter
      Beta        = GetBetaParameter(Aj,Zj,Ad,Zd) !beta parameter
      Coulomb     = GetCoulombBarrier(Aj,Zj,Ad,Zd) !Coulomb barrier [MeV]

      Ux = (2.5 + 150.0/float(Ad))
      Exx= Ux + delta0
      T  = 1./(sqrt(a/Ux) - 1.5/Ux)
      E0 = Exx - T*(log(T) - log(a)/4.0 - 1.25*log(Ux) + 2.0*sqrt(a*Ux))
      !---------------------------------------------------------------------
      ! Maximal kinetic energy E_max = Ex - Q - V
      !---------------------------------------------------------------------
      IDi = NucPosition(Ni,Zi)
      IDj = NucPosition(Nj,Zj)
      IDd = NucPosition(Nd,Zd)
      Qi  = Elements(IDi)%MassEx
      Qj  = Elements(IDj)%MassEx
      Qd  = Elements(IDd)%MassEx
!      MaximalKineticEnergy = ((Qi+U-delta0)**2+Qj**2-Qd**2)/ & 
!           &                 (2.*(Qi+U-delta0)) - Qj
      MaximalKineticEnergy = U - (Qj+Qd-Qi) - delta0
      !---------------------------------------------------------------------
      !(Excitation energy - (Q-Value) - Coulomb) must be positive:
      !---------------------------------------------------------------------
      if ((MaximalKineticEnergy-Coulomb) < 0.) then
         Width_k = 0.
         return
      endif
      !---------------------------------------------------------------------
      ! Get initial level density rho_i and start to determine width
      !---------------------------------------------------------------------
      if ( MaximalKineticEnergy - Coulomb < Exx ) then
         tm                  = (MaximalKineticEnergy-Coulomb)/T
         Width_k               = (I1(tm,tm)*T + (Beta+Coulomb)*I0(tm))/exp(E0/T)
         InitialLevelDensity = (pi/12.0)*exp((U-E0)/T)/T
      else 
         tm = (MaximalKineticEnergy-Coulomb)/T
         tx = Exx/T
         s  = 2.0*sqrt(a*(MaximalKineticEnergy-Coulomb-delta0))
         sx = 2.0*sqrt(a*(Exx-delta0))
        Width_k = I1(tm,tx)/exp(E0/T) + I3(s,sx,a)*exp(s)
        !For charged particles (Beta+Coulomb) = 0 because Beta = -Coulomb
        if (Zj == 0) then
           Width_k = Width_k + & 
       &  (Beta+Coulomb)*(I0(tx)/exp(E0/T)+I2(s,sx)*exp(s))
        endif
        InitialLevelDensity = (pi/12.0)*exp(2.*sqrt(a*(U-delta0)))/ & 
              &               (a**(1./4.)*(U-delta0)**(5./4.))
      endif
      !---------------------------------------------------------------------
      if (InitialLevelDensity < 0.00000001) then
         write(*,*) 'from GetDecayWidth: Leveldensity=0!!!'
         write(*,*) U,delta0,InitialLevelDensity
         STOP
      endif
      !---------------------------------------------------------------------
      ! Get spin factor
      !---------------------------------------------------------------------
      g = (2.0*SpinAj+1.0)*EjectileMass/(pi*pi)
      !---------------------------------------------------------------------
      ! Get geometrical cross section
      !---------------------------------------------------------------------
      if (Aj > 4) then
         R1 = float(Ad)**(1./3.)
         R2 = float(Aj)**(1./3.)
         Rb = 1.12*(R1 + R2) - 0.86*((R1+R2)/(R1*R2)) + 2.85 !fm
      else 
         R1 = float(Ad)**(1./3.)
         R2 = float(Aj)**(1./3.)
         Rb = 1.5*(R1+R2) !fm
      endif
      GeometricalXS = pi*Rb*Rb/(197.33*197.33)
      !---------------------------------------------------------------------
      ! get constant factors
      !---------------------------------------------------------------------
      ConstantFactor = g*GeometricalXS*Alpha/InitialLevelDensity
      ConstantFactor = ConstantFactor*pi/12.0
      !---------------------------------------------------------------------
      ! get total decay width
      !---------------------------------------------------------------------
      Width_k = ConstantFactor*Width_k

    !***********************************************************************
    end subroutine DecayWidth !*********************************************
    !***********************************************************************

    !***********************************************************************
    ! here the functions I0,I1,I2 and I3 needed in DecayWidth
    !***********************************************************************
    real function I0(t)
      real, intent(in) :: t
      I0 = (exp(t) - 1.0)
    end function I0
    !***********************************************************************
    real function I1(t,tx)
      real, intent(in) :: t,tx
      real :: temp
      temp = t - tx + 1.0
      temp = temp*exp(tx)
      temp = temp - (t + 1.0)
      I1   = temp
    end function I1
    !***********************************************************************
    real function I2(sa,sax)
      real, intent(in) :: sa,sax
      real :: S,Sx,p1,p2
      real, parameter :: const = 2.82843 !2*sqrt(2)

      S  = 1.0/sqrt(sa)
      Sx = 1.0/sqrt(sax)
    
      p1 = S*S*S*( 1.0 + S*S*( 1.5 + 3.75*S*S) )
      p2 = Sx*Sx*Sx*( 1.0 + Sx*Sx*( 1.5 + 3.75*Sx*Sx) )*exp(sax-sa)
    
      I2 = const*(p1-p2)
    end function I2
    !***********************************************************************
    real function I3(sa,sax,a)
      real, intent(in) :: sa,sax,a
      real :: s2,sx2,S,St2,Stx,Stx2,p1,p2
      real, parameter :: const = 1.414 !sqrt(2)

      s2 = sa*sa
      sx2 = sax*sax
      S = 1.0/sqrt(sa)
      St2 = S*S
      Stx = 1.0/sqrt(sax)
      Stx2 = Stx*Stx
    
      p1 = S *(2.0 + St2 *( 4.0 + St2 *( 13.5 + St2 *( 60.0 + St2 * 325.125 ))))
      p2 = Stx*Stx2 *( & 
           & (s2-sx2) + Stx2 *( & 
           &  (1.5*s2+0.5*sx2) + Stx2 *( & 
           &      (3.75*s2+0.25*sx2) + Stx2 *( & 
           &          (12.875*s2+0.625*sx2) + Stx2 *( & 
           &              (59.0625*s2+0.9375*sx2) + Stx2 *(324.8*s2+3.28*sx2))))))
      p2 = p2*exp(sax-sa)
      I3 = (p1-p2)/(1.414*a)

    end function I3
    !***********************************************************************


    !***********************************************************************
    !****s* GEM_Evaporation/Get4Momentum
    ! NAME
    ! subroutine Get4Momentum
    ! 
    ! PURPOSE
    ! * Calculates total kinetic energy of fragment according GEM-Probability 
    !   in the rest frame of the evaporating residual cluster
    ! * Calculates angular distribution of motion isotropically
    ! * Calculates 4-momenta of the daughter cluster
    ! * Finaly boost the 4-momenta of ejectile & Daughter clusters into the 
    !   global CMS (BUU calculational frame)
    !***********************************************************************
    subroutine Get4Momentum(Ai,Zi,Ad,Zd,A,Z,Ex,Normalization, & 
           &                betax,betay,betaz,Impulsj,Impulsd)

      use initGEM, only       : Nuclide,Elements,NucPosition
      use ClusterEnergy, only : boost
      use random, only: rn

      !---------------------------------------------------------------------
      !Input-Variables
      !---------------------------------------------------------------------
      integer, intent(in) :: Ai,Zi,Ad,Zd,A,Z
      real,    intent(in) :: Ex,betax,betay,betaz,Normalization
      !---------------------------------------------------------------------
      !Output-Variables
      !---------------------------------------------------------------------
      real, dimension(0:3), intent(out) :: Impulsj,Impulsd
      !---------------------------------------------------------------------
      !Local-Variables & external functions
      !---------------------------------------------------------------------

      integer :: i,MaxTest
      real    :: NormTest,MaxProb


      integer :: N,Index,IDi,IDj,IDd,isteps
      real    :: Qi,Qj,Qd,Factor
      real    :: U,EjectileMass,DaughterMass,alevel,Exx
      real    :: delta0,Alpha,Beta,KineticEnergy !,rand
      real    :: Ux,T,E0,MaximalKineticEnergy,g,Spin,Probability
      real    :: ConstantFactor,R1,R2,Rb,GeometricalXS,theEnergy
      real    :: InitialLevelDensity,Ekin,ModMom,CoulombBarrier
      real, dimension(1:3) :: Vector
      real, dimension(0:3) :: MomentumEj,MomentumDt
      !---------------------------------------------------------------------
      ! Get input parameters
      !---------------------------------------------------------------------
      N              = A - Z !NeutronNumber
      U              = Ex !Excitation energy
      EjectileMass   = real(A)*939. !fragment mass of ejectile
      DaughterMass   = real(Ad)*939. !fragment mass of residual
      alevel         = GetLevelDensityParameter(Ad,Zd,U) !Level density parameter a
!      alevel         = float(Ad)/8. !1/MeV
      delta0         = GetPairingCorrection(Ad,Zd) !pairing energy
      Alpha          = GetAlphaParameter(A,Z,Ad) !alpha parameter
      Beta           = GetBetaParameter(A,Z,Ad,Zd) !beta parameter
      CoulombBarrier = GetCoulombBarrier(A,Z,Ad,Zd) !Coulomb barrier

      Ux = (2.5 + 150.0/float(Ad))
      Exx= Ux + delta0
      T  = 1.0/(sqrt(alevel/Ux) - 1.5/Ux)
      E0 = Exx - T*(log(T) - log(alevel)/4.0 - 1.25*log(Ux) + & 
         & 2.0*sqrt(alevel*Ux))
      !---------------------------------------------------------------------
      ! Get initial level density rho_i
      !---------------------------------------------------------------------
      if ( U < Exx ) then
         InitialLevelDensity = (pi/12.0)*exp((U-E0)/T)/T
      else 
         InitialLevelDensity = (pi/12.0)*exp(2.*sqrt(alevel*(U-delta0)))/ & 
              & (alevel**(1./4.)*(U-delta0)**(5./4.))
      endif
      !---------------------------------------------------------------------
      ! Get spin factor
      !---------------------------------------------------------------------
      Index = NucPosition(N,Z)
      Spin  = Elements(Index)%Spin
      g     = (2.0*Spin+1.0)*EjectileMass/(pi*pi)
      !---------------------------------------------------------------------
      ! Get geometrical cross section
      !---------------------------------------------------------------------
      if (A > 4) then
         R1 = float(Ad)**(1./3.)
         R2 = float(A)**(1./3.)
         Rb = 1.12*(R1 + R2) - 0.86*((R1+R2)/(R1*R2)) + 2.85 !fm
      else 
         R1 = float(Ad)**(1./3.)
         R2 = float(A)**(1./3.)
         Rb = 1.5*(R1+R2) !fm
      endif
      GeometricalXS = pi*Rb*Rb/(197.33*197.33)
      !---------------------------------------------------------------------
      ! Maximal kinetic energy E_max = Ex - Q - V
      !---------------------------------------------------------------------
      IDi = NucPosition((Ai-Zi),Zi)
      IDj = Index
      IDd = NucPosition((Ad-Zd),Zd)
      Qi  = Elements(IDi)%MassEx
      Qj  = Elements(IDj)%MassEx
      Qd  = Elements(IDd)%MassEx
!      MaximalKineticEnergy = ((Qi+U-delta0)**2+Qj**2-Qd**2)/ & 
!           &                 (2.*(Qi+U-delta0)) - Qj - CoulombBarrier

      MaximalKineticEnergy = U - (Qj+Qd-Qi) - delta0 - CoulombBarrier
      
      if (MaximalKineticEnergy < 0) then
         write(*,*) 'from Get4Momentum: Max_Ekin<0???'
         STOP
      endif
      !---------------------------------------------------------------------
      ! get constant factors
      !---------------------------------------------------------------------
      ConstantFactor = g*GeometricalXS*Alpha/InitialLevelDensity
      ConstantFactor = ConstantFactor*pi/12.0
      theEnergy      = MaximalKineticEnergy + CoulombBarrier
      !---------------------------------------------------------------------
      ! Start monte-carlo for ejectile energy according GEM-Probability
      ! - one random number for selecting kinetic energy X between 
      !   E_min and E_max
      ! - another random number to decide whether X is accepted 
      !   (if P(X) > rand()) or rejected (P(X)<rand())
      !---------------------------------------------------------------------
      ! Find max. value X_max, at which P(X_max)=P_max
      ! P_max is the normalization for Probability
      !---------------------------------------------------------------------
      KineticEnergy = CoulombBarrier
      MaxTest       = (theEnergy - CoulombBarrier)/0.01
      NormTest      = 0.0
      MaxProb = 0.0
      do i=1,MaxTest
         Probability   = 0.
         Factor   = ConstantFactor*(KineticEnergy + Beta)
         if ( theEnergy-KineticEnergy < Exx) then
            Probability = Factor* & 
                 & ( exp((theEnergy-KineticEnergy-E0)/T)/T )
         else 
            Probability = Factor* & 
                 & ( exp(2.*sqrt(a*(theEnergy-KineticEnergy-delta0))))/ & 
                 & ( a**(1./4.)*(theEnergy-KineticEnergy-delta0)**(5./4.) )
         endif
         if (Probability > MaxProb) then
            MaxProb = Probability
         endif
         NormTest = NormTest + Probability*0.01
         kineticEnergy = kineticEnergy + 0.01
      end do
      !---------------------------------------------------------------------
      Probability = 0. !starting value
      isteps      = 0
      do while (rn() > Probability)
         if (isteps > 10000) then
            write(*,*) 'Get4Momentum: infinite loop...',& 
                 & isteps,Normalization,Probability,Ai,Ad,A
            stop
         endif
         KineticEnergy =  CoulombBarrier + rn()*MaximalKineticEnergy
         Factor   = ConstantFactor*(KineticEnergy + Beta)
         if ( theEnergy-KineticEnergy < Exx) then
            Probability = Factor* & 
                 & ( exp((theEnergy-KineticEnergy-E0)/T)/T )
         else 
            Probability = Factor* & 
                 & ( exp(2.*sqrt(a*(theEnergy-KineticEnergy-delta0))))/ & 
                 & ( a**(1./4.)*(theEnergy-KineticEnergy-delta0)**(5./4.) )
         endif
         Probability = Probability/MaxProb
         isteps = isteps + 1
      end do
      !---------------------------------------------------------------------
      ! Get 4-momenta
      !---------------------------------------------------------------------
      !Ejectile
      Ekin   = KineticEnergy
      ModMom = sqrt((Ekin+EjectileMass)*(Ekin+EjectileMass)-EjectileMass*EjectileMass)
      call IsotropicVector(Modmom,Vector)
      MomentumEj(0)   = Ekin+EjectileMass
      MomentumEj(1:3) = Vector(1:3)
      !Daughter
      MomentumDt(1:3) = -Vector(1:3)
      MomentumDt(0)   = sqrt(DaughterMass*DaughterMass+Sum(MomentumDt(1:3)*MomentumDt(1:3)))
      !and boost them back to calculational BUU frame
!      checkEj = sqrt(MomentumEj(0)**2-Sum(MomentumEj(1:3)*MomentumEj(1:3)))
!      checkDt = sqrt(MomentumDt(0)**2-Sum(MomentumDt(1:3)*MomentumDt(1:3)))
      call boost(betax,betay,betaz,Impulsj,MomentumEj)
      call boost(betax,betay,betaz,Impulsd,MomentumDt)
    !***********************************************************************
    end subroutine Get4Momentum !*******************************************
    !***********************************************************************


    !***********************************************************************
    subroutine IsotropicVector(ModMom,Vector)
    !***********************************************************************
      use random, only: rn
      real, intent(in)                  :: ModMom
      real, dimension(1:3), intent(out) :: Vector
      real                              :: CosTheta,SinTheta,Phi !,rand
      
      CosTheta  = 1.0 - 2.0*rn()
      SinTheta  = sqrt(1.0 - CosTheta*CosTheta)
      Phi       = pi*pi*rn()
      Vector(1) = ModMom*cos(Phi)*SinTheta
      Vector(2) = ModMom*sin(Phi)*SinTheta
      Vector(3) = ModMom*CosTheta
    !***********************************************************************
    end subroutine IsotropicVector !****************************************
    !***********************************************************************

    !***********************************************************************
    real function GetPairingCorrection(A,Z)
    !***********************************************************************
      integer, intent(in) :: A,Z
      integer             :: N,IN,IZ
      real                :: delta
      !---------------------------------------------------------------------
      N     = A-Z
      IN    = mod(N,2)
      IZ    = mod(Z,2)
      delta = 0.0
      if (IN==0 .and. IZ==0) delta = 11.46/sqrt(float(A))
      if (IN==1 .and. IZ==1) delta = -11.46/sqrt(float(A))
      GetPairingCorrection = delta
    !***********************************************************************
    end function GetPairingCorrection !*************************************
    !***********************************************************************

    !***********************************************************************
    real function GetLevelDensityParameter(A,Z,Ex)
    !***********************************************************************
      integer, intent(in) :: A,Z
      real,    intent(in) :: Ex
      integer             :: N
      real                :: delta,ul,ail,atl,sl,apl,temp
      !---------------------------------------------------------------------
      
      N     = A-Z
      delta = GetPairingCorrection(A,Z)
      ul    = 0.05*(Ex-delta)
      ail   = (0.1375-0.0000836*real(A))*real(A)
      if (Z < 9) then
         atl = real(A)/8.
      else
         sl  = 0.0 !Shell corrections...skip them for the moment
         if ((Z.ge.54.and.Z.le.78).or.(Z.ge.86.and.Z.le.98).or. & 
          &  (N.ge.86.and.N.le.122).or.(N.ge.130.and.N.le.150) ) then
            apl = 0.12
         else
            apl = 0.142
         endif
         atl = real(A)*(apl+0.00917*sl)
      endif
      !---------------------------------------------------------------------
      temp = (1-exp(-ul))/ul
      GetLevelDensityParameter = atl*temp+ail*(1.-temp)

    !***********************************************************************
    end function GetLevelDensityParameter !*********************************
    !***********************************************************************


    !***********************************************************************
    real function GetAlphaParameter(Aj,Zj,Ad)
    !***********************************************************************
      integer, intent(in) :: Aj,Zj,Ad
      real                :: cj
      !---------------------------------------------------------------------
      !neutrons
      if (Zj==0) then
         GetAlphaParameter = 0.76+1.93*float(Ad)**(-0.333333333)
         return
      endif
      !charged particles
      cj = 1.0
      if (Aj==1 .and. Zj==1) cj = 0.51
      if (Aj==2 .and. Zj==1) cj = 1.255
      if (Aj==3 .and. Zj==1) cj = 1.17
      if ( (Aj==4 .and. Zj==2) .or. (Aj==3.and.Zj==2) ) cj = 0.
      GetAlphaParameter = 1.+cj

    !***********************************************************************
    end function GetAlphaParameter !****************************************
    !***********************************************************************

    !***********************************************************************
    real function GetBetaParameter(Aj,Zj,Ad,Zd)
    !***********************************************************************
      integer, intent(in) :: Aj,Zj,Ad,Zd
      integer             :: Nj
      real                :: alpha,Coulomb
      !---------------------------------------------------------------------
      Nj = Aj-Zj
      !neutrons
      if (Zj==0) then
         alpha            = GetAlphaParameter(Aj,Zj,Ad)
         GetBetaParameter = (1.66*float(Ad)**(-2./3.)-0.05)/alpha
         return
      endif
      !charged particles
      Coulomb = GetCoulombBarrier(Aj,Zj,Ad,Zd)
      GetBetaParameter = -Coulomb


    !***********************************************************************
    end function GetBetaParameter !*****************************************
    !***********************************************************************


    !***********************************************************************
    real function GetCoulombBarrier(Aj,Zj,Ad,Zd)
    !***********************************************************************
      integer, intent(in) :: Zj,Zd,Aj,Ad
      real, parameter     :: ladung = 1./137.
      real                :: Rc,kj
      !---------------------------------------------------------------------
      !neutral particles (neutrons)
      if (Zj==0 .or. Zd==0) then
         GetCoulombBarrier = 0.0
         return
      endif
      !charged particles

      if (Aj == 1) Rc = 1.7*Aj**0.3333
      if (Aj > 1 .and. Aj <= 4) then
         Rc = 1.7*Aj**0.3333 + 1.2
      endif
      if (Aj > 4) then
         Rc = 1.12*(Aj**0.3333+Ad**0.33333)-0.86*(Aj**(-0.3333)+Ad**(-0.33333))
      endif

      if (Aj < 4) then
         if (Aj==1 .and. Zj==1) kj = 0.51
         if (Aj==2 .and. Zj==1) kj = 0.57
         if (Aj==3 .and. Zj==1) kj = 0.63
         if (Aj==3 .and. Zj==2) kj = 0.75
      else
         kj = 1.
      endif
      GetCoulombBarrier = ( kj*Zj*Zd*ladung/Rc )*197.33 !MeV

    !***********************************************************************
    end function GetCoulombBarrier !****************************************
    !***********************************************************************





  end module GEM_Evaporation
