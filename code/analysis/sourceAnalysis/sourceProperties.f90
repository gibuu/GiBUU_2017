!******************************************************************************
!****m* /sourceProperties
! NAME
! module sourceProperties
! PURPOSE
! Administrates the calculation of local thermodynamical quantities
! NOTES
! * Calculation of longitudinal and transversal pressure densities.
! * This module can be used in RMF- and Skyrme-mode, however, in Skyrme-mode
!   the calculation of the pressures has to be still implemented.
! * This module estimates the time for switching from dynamical BUU to
!   statistical fragmentation.
!******************************************************************************
module sourceProperties

  use constants, only: singlePrecision

  private


  !****************************************************************************
  !****g* sourceProperties/TheTensor
  ! SOURCE
  !
  real(singlePrecision), allocatable, dimension(:,:,:,:,:), SAVE :: TheTensor
  ! PURPOSE
  ! Energy-momentum tensor T^{\mu\nu}(x) at each grid point.
  ! NOTES
  ! * The structure of the TheTensor:
  ! * TheTensor(:,:,:,:,:) = T(x,y,z,0:3,0:3)
  !
  !****************************************************************************

  !****************************************************************************
  !****g* sourceProperties/PressureField
  ! SOURCE
  !
  real(singlePrecision), dimension(:,:), allocatable, SAVE :: PressureField
  ! PURPOSE
  ! Components of the pressure density field at the center of the source(s).
  ! NOTES
  ! * 1-Index: Source label.
  ! * 2-Index: diagonal components (P_xx,P_yy,P_zz) of T^{\mu\nu}(x)
  !****************************************************************************

  !****************************************************************************
  !****g* sourceProperties/Anisotropy
  ! SOURCE
  !
  real(singlePrecision), dimension(:), allocatable, SAVE :: Anisotropy
  ! PURPOSE
  ! Anisotropy factor Q(x)=2*T^{zz}(x)-(T^{xx}(x)+T^{yy}(x)) at the
  ! center of the source(s).
  !****************************************************************************

  !****************************************************************************
  !****g* sourceProperties/MeanAnisotropy
  ! SOURCE
  !
  real, SAVE :: meanAnisotropy = 0.0
  ! PURPOSE
  ! Average anisotropy factor <Q> (over the sources).
  !****************************************************************************

  !****************************************************************************
  !****g* sourceProperties/InvariantDensity
  ! SOURCE
  !
  real(singlePrecision), dimension(:), allocatable, SAVE :: InvariantDensity
  ! PURPOSE
  ! Invariant Density at the center of the source(s).
  !****************************************************************************


  public :: sourceProperties_Main,WriteSourceInfo,deallocate_sourceFields

contains


  !****************************************************************************
  !****s* sourceProperties/sourceProperties_Main
  ! NAME
  ! subroutine sourceProperties_Main(time,NumSources,TheSource,instability_Flag)
  ! PURPOSE
  ! main routine: administrates the calculation of local thermodynamical
  ! properties of the source(s) and the determination of local spinodal instabilities.
  ! INPUTS
  ! * real,                             :: time         -- actual time (fm/c)
  ! * integer                           :: NumEnsemples -- number of ensemples
  ! * type(particle), dimension(:,:),   :: realPV       -- the real particle vector
  ! * integer,                          :: NumSources   -- number of valid sources
  !****************************************************************************
  subroutine sourceProperties_Main(time,NumEnsemples,NumSources,realPV)

    use RMF, only: getRMF_flag
    use particleDefinition
    use sourceTypeDefinition
    use determineSource, only: TheSource,MaxNumSources

    implicit none
    !-----------------------------------------------------------------------
    !Input-Output variables
    !-----------------------------------------------------------------------
    type(particle), dimension(:,:),           intent(in)  :: realPV
    integer,                                  intent(in)  :: NumEnsemples
    real,                                     intent(in)  :: time
    integer, dimension(1:Size(realPV,dim=1)), intent(in)  :: NumSources
    !-----------------------------------------------------------------------
    !Local variables
    !-----------------------------------------------------------------------
    logical, SAVE :: firstCall   = .true.
    !-----------------------------------------------------------------------
    if (firstCall) then
       if ( .not. getRMF_flag() ) then
          write(*,*) '  sourceProperties/Instability_Main: '
          write(*,*) '  Invalid mode for mean field! '
          write(*,*) '  Termination of program now'
          STOP
       end if
       firstCall = .false.
    end if
    !-----------------------------------------------------------------------
    !    Determine at each time step delta_T the energy-momentum tensor
    !    on the grid.
    !    The local pressure components and the density are calculated also
    !    at the center of the source(s).
    !-----------------------------------------------------------------------
    call EnergyMomentumTensor(realPV)
    call allocate_sourceFields(MaxNumSources)
    call get_PressureField(NumEnsemples,MaxNumSources,NumSources,TheSource)
    !-----------------------------------------------------------------------
    ! Printing the time dependence of the average (over the centroids of the
    ! valid sources) of the anisotropy factor. This file gives important
    ! information on achieved degree of equilibration at that time, when
    ! statistical fragmentation is applied.
    write(*,*)
    write(*,'(A,(f12.5))') '### <Anisotropy> = ',MeanAnisotropy
    write(700,222) time,MeanAnisotropy
222 format(2f15.5)

  !****************************************************************************
  end subroutine SourceProperties_Main !**********************************
  !****************************************************************************

  !****************************************************************************
  !****s* sourceProperties/allocateFields
  ! PURPOSE
  ! Allocation of the energy-momentum Tensor at the first call.
  ! NOTES
  ! The allocation of the tensor is done only at the first call, since
  ! its size remains constant. For other fields of this module, see
  ! notes in the routines (de)allocate_sourceFields.
  !****************************************************************************
  subroutine allocateFields

    use densitymodule, only: gridPoints

    allocate(TheTensor(-gridPoints(1):gridPoints(1),&
         &    -gridPoints(2):gridPoints(2),&
         &    -gridPoints(3):gridPoints(3),&
         &     0:3,0:3))

  !****************************************************************************
  end subroutine allocateFields !*****************************************
  !****************************************************************************

  !****************************************************************************
  !****s* sourceProperties/get_localPressure
  ! PURPOSE
  ! Determination of the local Pressure density, baryon density and
  ! of local anisotropies at the center of the source(s), as well as
  ! determination of the average anisotropy ratio.
  ! INPUTS
  ! integer               :: NumEnsemples     -- Number of ensemples
  ! integer               :: MaxNumSources    -- Max. number of valid sources
  ! integer,dimension(:)  :: NumSources       -- Number of valid sources
  ! type(quelle), dimension(:,:) :: TheSource -- source(s) properties
  !
  ! USES
  ! sourceTypeDefinition, densitymodule,lorentzTrafo
  !****************************************************************************
  subroutine get_PressureField(NumEnsemples,MaxNumSources,NumSources,TheSource)

    use densitymodule, only: densityField, gridSpacing, gridPoints !,gridSize
    use lorentzTrafo, only: BoostTensor
    use sourceTypeDefinition

    implicit none
    !-----------------------------------------------------------------------
    ! Input-Output variables
    !-----------------------------------------------------------------------
    integer,                                 intent(in) :: NumEnsemples,MaxNumSources
    type(quelle), dimension(:,:),            intent(in) :: TheSource
    integer,      dimension(1:NumEnsemples), intent(in) :: NumSources
    !-----------------------------------------------------------------------
    ! Local variables
    !-----------------------------------------------------------------------
    integer                             :: i,i1,i2,i3,m
    real                                :: localRho,MeanAn,NormAn
    real,   dimension(1:MaxNumSources)  :: Dichte,NonEquil,NormS,Norm,Norm2
    real,   dimension(1:MaxNumSources,1:3) :: Druck
    real,   dimension(1:3)              :: v,LocalPos
    integer,dimension(1:3)              :: index
    real,   dimension(0:3)              :: uu,strom
    real,   dimension(0:3,0:3)          :: ee,ee0
    !-----------------------------------------------------------------------
    InvariantDensity(:)  = 0.0 !rho at the center of the source(s)
    PressureField(:,:)     = 0.0 !P           -//-
    Anisotropy(:)        = 0.0 !Q           -//-
    NormS(:)             = 0.0 !normalization factor
    MeanAnisotropy       = 0.0 !<Q> over coordinate space
    !-----------------------------------------------------------------------

    NormAn = 0.0
    MeanAn = 0.0

    Loop_over_Ensemples : do m=1,NumEnsemples

       !skip here events with non-existing projectile/target/fireball
!       if (NumSources(m) > 1 .and. NumSources(m) < 3) cycle

       Loop_over_Sources : do i=1,NumSources(m)

          if (i > MaxNumSources) then
             write(*,*) 'module sourceProperties, routine get_PressureField:'
             write(*,*) 'Wrong determination of NumSources(:) and/or MaxNumSOurces'
             write(*,*) '!!! Termination of the program !!!'
             STOP
          end if

          if (.not.TheSource(m,i)%status) cycle

          Norm(i)     = 0.0
          Norm2(i)    = 0.0
          Druck(i,:)    = 0.0
          NonEquil(i) = 0.0
          Dichte(i)   = 0.0

          LocalPos(1:3) = TheSource(m,i)%Position(1:3)
          index         = nint(LocalPos/gridSpacing)

          Loop_xSpace : do i1=index(1)-2,index(1)+2
             Loop_ySpace : do i2=index(2)-2,index(2)+2
                Loop_zSpace : do i3=index(3)-2,index(3)+2

                   if (   abs(i1) > gridPoints(1) .or. &
                        & abs(i2) > gridPoints(2) .or. &
                        & abs(i3) > gridPoints(3) ) cycle

                   strom(:) = densityField(i1,i2,i3)%baryon(:)
                   localRho = sqrt( max(0.,(strom(0)**2 -dot_product(strom(1:3), strom(1:3)))) )
                   Dichte(i) = Dichte(i) + localRho

                   !collective 4-velocity of the source(s):
                   v(1:3) = TheSource(m,i)%velocity(1:3)
                   uu(0)  = 1./sqrt(1.-v(1)**2-v(2)**2-v(3)**2)
                   uu(1:3)= uu(0)*v(1:3)
                   !performe boost of the tensor to the RF of the source(s):
                   ee(0:3,0:3) = TheTensor(i1,i2,i3,0:3,0:3)
                   call boostTensor(uu,ee,ee0)

                   !local pressure components
                   Druck(i,1) = Druck(i,1) + ee0(1,1)
                   Druck(i,2) = Druck(i,2) + ee0(2,2)
                   Druck(i,3) = Druck(i,3) + ee0(3,3)

                   if (localRho > 0.00001) then
                      NonEquil(i) = NonEquil(i) + ( 2.*ee0(3,3) - (ee0(1,1) + ee0(2,2)) )
                      Norm2(i) = Norm2(i) + 1.
                   end if

                   Norm(i) = Norm(i) + 1.

                end do Loop_zSpace
             end do Loop_ySpace
          end do Loop_xSpace

          if (Norm(i) /= 0.)  Dichte(i)   = Dichte(i) / Norm(i)
          if (Norm(i) /= 0.)  Druck(i,:)    = Druck(i,:)  / Norm(i)
          if (Norm2(i) /= 0.) NonEquil(i) = NonEquil(i) / Norm2(i)

          InvariantDensity(i) = InvariantDensity(i) + Dichte(i)
          PressureField(i,:)    = PressureField(i,:)    + Druck(i,:)
          Anisotropy(i)       = Anisotropy(i)       + NonEquil(i)

          NormS(i) = NormS(i) + 1.

          MeanAn = MeanAn + NonEquil(i)
          NormAn = NormAn + 1.

       end do Loop_over_Sources

    end do Loop_over_Ensemples

    InvariantDensity(:) = InvariantDensity(:) / max(NormS(:),0.0001)
    PressureField(:,1)    = PressureField(:,1) / max(NormS(:),0.0001)
    PressureField(:,2)    = PressureField(:,2) / max(NormS(:),0.0001)
    PressureField(:,3)    = PressureField(:,3) / max(NormS(:),0.0001)
    Anisotropy(:)       = Anisotropy(:) / max(NormS(:),0.0001)

    if (NormAn == 0.0) then
       MeanAnisotropy = 0.0
    else
       MeanAnisotropy = MeanAn / NormAn
    end if

  !****************************************************************************
     end subroutine get_PressureField !***********************************
  !****************************************************************************


  !****************************************************************************
  !****s* sourceProperties/EnergyMomentumTensor
  ! PURPOSE
  ! Determination of the energy-momentum tensor at each grid point.
  ! INPUTS
  ! * type(particle), dimension(:,:), intent(in) :: teilchen
  ! OUTPUT
  ! The Energy-Momentum-Tensor Tensor(:,:,:,:,:) T^{\mu\nu}(x)
  ! NOTES
  ! * Valid only the RMF mode!
  ! * Calculation of T^{\mu\nu} in the same system as done for the currents.
  !****************************************************************************
  subroutine EnergyMomentumTensor(teilchen)

    use particleDefinition, only: particle
    use IdTable, only: pion
    use constants, only: hbarc
    use RMF
    use densitymodule

    implicit none
    !-----------------------------------------------------------------------
    type(particle), dimension(:,:), intent(in) :: teilchen

    real,    dimension(1:3) :: rpos
    integer, dimension(1:3) :: posOrig,IndexSmall
    real,    dimension(0:3) :: ujj,strom,stromi

    integer :: i,j,m,l,small,large,Index1,Index2,Index3
    real    :: mstar,stromq,stromiq,factor,phi
    real    :: ScalarPart,VectorPart

    real,   dimension(0:3,0:3), SAVE :: Metrik=0.
    logical,                    SAVE :: start_Flag  = .true.
    !-----------------------------------------------------------------------
    if (start_Flag) then
       call allocateFields
       Metrik(0,0) = 1.0
       do i=1,3
          Metrik(i,i) = -1.0
       end do
       start_Flag = .false.
    end if
    !-----------------------------------------------------------------------
    ! initialize the tensor:
    ! 1-st 3 indexes: coordinates on the 3D-grid
    ! last 2 indexes: space-time components of the tensor
    ! TheTensor(:,:,:,:,:) = T(xi,yi,zi,0:3,0:3)
    !-----------------------------------------------------------------------
    TheTensor(:,:,:,:,:) = 0.0
    !-----------------------------------------------------------------------
    ! calculation of the kinetic contribution t_{kin}^{\mu\nu}(x)
    !-----------------------------------------------------------------------
    Loop_over_ensembles_1 : do j=1,Size(Teilchen,dim=1)
       Loop_over_particles_1 : do i=1,Size(Teilchen,dim=2)

          if ( teilchen(j,i)%id == 0 ) then
             cycle Loop_over_particles_1
          else if ( teilchen(j,i)%id < 0 ) then
             exit Loop_over_particles_1
          end if

          if ( teilchen(j,i)%id >= pion) cycle  Loop_over_particles_1 ! Only baryons are accounted for presently

          ! Modification factor for the coupling constants:
          factor = ModificationFactor(teilchen(j,i)%Id,teilchen(j,i)%antiparticle)

          ! position in large grid:
          rpos = Teilchen(j,i)%position(1:3)/gridSpacing
          posOrig=NINT(rpos)

          if (abs(posOrig(1)).gt.gridPoints(1) &
               &.or. abs(posOrig(2)).gt.gridPoints(2) &
               &.or. abs(posOrig(3)).gt.gridPoints(3) ) cycle

          ! position in small grid:
          indexSmall=NINT((rpos-posOrig)*(2.*SmallergridPoints+1.))

          ! Test for errors:
          if ((abs(indexSmall(1)).gt.smallergridpoints)&
               &.or.(abs(indexSmall(2)).gt.smallergridpoints)&
               &.or.(abs(indexSmall(3)).gt.smallergridpoints)) then
             write(*,*) 'Problem in InstabilityConditions, module FragmentationAnalysis.f90'
             write(*,*) IndexSmall, 'too big, choose different grid'
             stop
          end if

          small=1+(SmallergridPoints+indexSmall(3))                               &
               & +(SmallergridPoints+indexSmall(2))*(2*SmallerGridPoints+1)       &
               & +(SmallergridPoints+indexSmall(1))*(2*SmallerGridPoints+1)**2

          large=0

          ! Smearing particle over points in Neighborhood:
          Loop_Index1 : do Index1=posOrig(1)-numberLargePoints,posOrig(1)+numberlargePoints
             Loop_Index2 :Do Index2=posOrig(2)-numberLargePoints,posOrig(2)+numberlargePoints
                Loop_Index3 :Do Index3=posOrig(3)-numberLargePoints,posOrig(3)+numberlargePoints

                   large=large+1

                   InsideGrid : if (       abs(Index1).le.gridPoints(1) &
                        & .and. abs(Index2).le.gridPoints(2) &
                        & .and. abs(Index3).le.gridPoints(3) &
                        & .and. smearingWeights(small,large).gt.0. ) then

                      mstar = teilchen(j,i)%mass + factor*g_sigma*sigmaField(Index1,Index2,Index3)
                      ujj(:) = teilchen(j,i)%momentum(:)/mstar
                      do m=0,3
                         do l=0,3
                            TheTensor(Index1,Index2,Index3,m,l) = TheTensor(Index1,Index2,Index3,m,l)  + &
                                 &  mstar*ujj(m)*ujj(l)*smearingWeights(small,large)/ujj(0)
                         end do
                      end do


                   end if InsideGrid

                end do Loop_Index3
             end do Loop_Index2
          end do Loop_Index1

       end do Loop_over_particles_1
    end do Loop_over_ensembles_1
    !-----------------------------------------------------------------------
    ! calculation of T^{\mu\nu(x) = }t_{kin}^{\mu\nu}(x) + g^{\mu\nu}*Potential(x)
    !-----------------------------------------------------------------------
    do Index1 = -gridPoints(1),gridPoints(1)
       do Index2 = -gridPoints(2),gridPoints(2)
          do Index3 = -gridPoints(3),gridPoints(3)

             !total baryon current
             strom(:) = densityField(Index1,Index2,Index3)%baryon(:)
             stromq   = strom(0)**2 - dot_product(strom(1:3),strom(1:3))

             !total isospin current (difference between proton and neutron currents)
             stromi(:) = &
                  & densityField(Index1,Index2,Index3)%proton(:)- &
                  & densityField(Index1,Index2,Index3)%neutron(:)
             stromiq   = stromi(0)**2 - dot_product(stromi(1:3),stromi(1:3))

             phi  = sigmaField(Index1,Index2,Index3)

             ScalarPart = ( 0.5*(m_sigma*phi)**2 + g_2*phi**3/3. + g_3*phi**4/4. ) / hbarc**3
             VectorPart = 0.5*a_6*stromq + 0.5*a_7*stromiq

             do m=0,3
                do l=0,3
                   TheTensor(Index1,Index2,Index3,m,l) = TheTensor(Index1,Index2,Index3,m,l)  &
                        & + a_6*strom(m)*strom(l) + a_7*stromi(m)*stromi(l) & !still kinetic contrinution
                        & + metrik(m,l)*( ScalarPart - VectorPart ) !potential contribution
                end do
             end do

          end do
       end do
    end do

  !****************************************************************************
  end subroutine EnergyMomentumTensor !***********************************
  !****************************************************************************


  !****************************************************************************
  !****s* sourceProperties/WriteSourceInfo
  ! PURPOSE
  ! Administrates the printing of source properties.
  ! NOTES
  ! * Print information on sources as function of time
  ! * PrintCotrol = true --> prints at onset of equilibration
  !   (SourceFile.dat)
  ! * FinalFlag   = true --> prints at onset of equilibration
  !   after forced decays (SourceFile_NoResonances.dat)
  ! * The file "SourceFile.dat" or "SourceFile_NoResonances.dat" serve
  !   as input for the statistical multifragmentation code.
  !****************************************************************************
  subroutine WriteSourceInfo(time,stossParameter,isut,NumSources,hyperSource,FinalFlag)

    use sourceTypeDefinition
    use RMF, only: getRMF_flag
    use determineSource, only: TheSource,MaxNumSources
    use output, only: realTochar

    implicit none
    !-----------------------------------------------------------------------
    ! Input variables
    !-----------------------------------------------------------------------
    integer, dimension(:), intent(in) :: NumSources
    integer,               intent(in) :: isut
    real,                  intent(in) :: time,stossParameter
    logical,               intent(in) :: hyperSource
    logical, optional,     intent(in) :: FinalFlag
    !-----------------------------------------------------------------------
    ! Local variables
    !-----------------------------------------------------------------------
    integer                                     :: i,m
    type(quelle),dimension(1:MaxNumSources)     :: MeanSource
    real,        dimension(1:MaxNumSources)     :: NormSource,Pmean,P_zz,P_tr
    integer, dimension(1:Size(TheSource,dim=1)) :: NumValidSources
    !-----------------------------------------------------------------------
    ! Determine number of valid sources for each ensemple
    !-----------------------------------------------------------------------
    NumValidSources(:) = 0
    do m=1,Size(TheSource,dim=1)
       do i=1,NumSources(m)
          if (.not.TheSource(m,i)%status) cycle
          NumValidSources(m) = NumValidSources(m) + 1
       end do
    end do
    !-----------------------------------------------------------------------
    ! Print source info at the end of the simulation after the forced decays
    ! and do nothing else!
    !-----------------------------------------------------------------------
    !Open file for printing source info at the end after forced decays:
    ForcedDecays : if (present(FinalFlag)) then
       open(105,file='SourceFile_NoResonances.dat',position='Append')

       do m=1,Size(TheSource,dim=1)

          write(105,'(A,1x,i5)') 'Number of valid sources: ',NumValidSources(m)

          do i=1,NumSources(m)

          if (.not.TheSource(m,i)%status) cycle

          if (.not.hyperSource) then
             write(105,100) time,&
                  & TheSource(m,i)%Size,TheSource(m,i)%Charge, &
                  & TheSource(m,i)%position,TheSource(m,i)%velocity, &
                  & TheSource(m,i)%ExEnergy,TheSource(m,i)%radEnergy,&
                  & stossParameter,m,isut+1
          else
             write(105,102) time,&
                  & TheSource(m,i)%Size,TheSource(m,i)%Charge, &
                  & TheSource(m,i)%nLambda,TheSource(m,i)%nSigma0, &
                  & TheSource(m,i)%position,TheSource(m,i)%velocity, &
                  & TheSource(m,i)%ExEnergy,TheSource(m,i)%radEnergy,&
                  & stossParameter,m,isut+1
          end if

       end do

    end do

    close(105)

    RETURN

    end if ForcedDecays

    !-----------------------------------------------------------------------
    ! Print source info into several files separated in time.
    !-----------------------------------------------------------------------

    !Open file for printing source info at different times:
    open(103,file='Source'//realTochar(time)//'.dat',position='Append')

    Loop_over_Ensemples : do m=1,Size(TheSource,dim=1)

       write(103,'(A,1x,i5)') 'Number of valid sources: ',NumValidSources(m)

       Loop_over_Sources : do i=1,NumSources(m)

          if (.not.TheSource(m,i)%status) cycle

          if (.not.hyperSource) then
             write(103,100) time,&
                  & TheSource(m,i)%Size,TheSource(m,i)%Charge, &
                  & TheSource(m,i)%position,TheSource(m,i)%velocity, &
                  & TheSource(m,i)%ExEnergy,TheSource(m,i)%radEnergy,&
                  & stossParameter,m,isut+1
          else
             write(103,102) time,&
                  & TheSource(m,i)%Size,TheSource(m,i)%Charge, &
                  & TheSource(m,i)%nLambda,TheSource(m,i)%nSigma0, &
                  & TheSource(m,i)%position,TheSource(m,i)%velocity, &
                  & TheSource(m,i)%ExEnergy,TheSource(m,i)%radEnergy,&
                  & stossParameter,m,isut+1
          end if

       end do Loop_over_Sources

    end do Loop_over_Ensemples

    close(103)

    !-----------------------------------------------------------------------
    ! Print source info on the terminal
    !-----------------------------------------------------------------------
    NormSource(:) = 0.0
    MeanSource(:)%Size = 0
    MeanSource(:)%Charge = 0
    MeanSource(:)%nLambda = 0
    MeanSource(:)%nSigma0 = 0
    MeanSource(:)%Position(1) = 0.0
    MeanSource(:)%Position(2) = 0.0
    MeanSource(:)%Position(3) = 0.0
    MeanSource(:)%Velocity(1) = 0.0
    MeanSource(:)%Velocity(2) = 0.0
    MeanSource(:)%Velocity(3) = 0.0
    MeanSource(:)%radEnergy = 0
    MeanSource(:)%ExEnergy = 0
    Loop_over_Ensemples2 : do i=1,Size(TheSource,dim=1)
!       if (NumSources(i) > 1 .and. NumSources(i) < 3) cycle
       Loop_over_Sources2 : do m=1,NumSources(i)
          if (.not.TheSource(i,m)%Status) cycle
          NormSource(m) = NormSource(m) + 1.
          MeanSource(m)%Size = MeanSource(m)%Size  + TheSource(i,m)%Size
          MeanSource(m)%Charge = MeanSource(m)%Charge  + TheSource(i,m)%Charge
          MeanSource(m)%nLambda = MeanSource(m)%nLambda  + TheSource(i,m)%nLambda
          MeanSource(m)%nSigma0 = MeanSource(m)%nSigma0  + TheSource(i,m)%nSigma0
          MeanSource(m)%Position(:) = MeanSource(m)%Position(:)  + TheSource(i,m)%Position(:)
          MeanSource(m)%Velocity(:) = MeanSource(m)%Velocity(:)  + TheSource(i,m)%Velocity(:)
          MeanSource(m)%radEnergy = MeanSource(m)%radEnergy  + TheSource(i,m)%radEnergy
          MeanSource(m)%ExEnergy = MeanSource(m)%ExEnergy  + TheSource(i,m)%ExEnergy
       end do Loop_over_Sources2
    end do Loop_over_Ensemples2

    open(106,file='SourceEvol.dat',position='Append')

    Loop_Sources : do i=1,MaxNumSources

       if (NormSource(i)==0.) cycle

       if (getRMF_Flag()) then
          Pmean(i) = ( Sum(PressureField(I,1:3))/3. ) * 1000. !output in MeV
          P_zz(i) = PressureField(i,3)*1000. !output in MeV
          P_tr(i) = 0.5*(PressureField(i,1)+PressureField(i,2))*1000. !output in MeV
       end if

       !--------------------------------------------------------------------
       write(*,*)
       write(*,'(A,f12.4)')       '### <A>                = ',float(MeanSource(i)%Size)/NormSource(i)
       write(*,'(A,f12.4)')       '### <Z>                = ',float(MeanSource(i)%Charge)/NormSource(i)
       if (hyperSource) then
          write(*,'(A,e13.6)')       '### <N_Lambda>         = ',float(MeanSource(i)%nLambda)/NormSource(i)
          write(*,'(A,e13.6)')       '### <N_Sigma0>         = ',float(MeanSource(i)%nSigma0)/NormSource(i)
       end if
       write(*,'(A,3(f12.4))')    '### <Pos> [fm]         = ',MeanSource(i)%Position/NormSource(i)
       write(*,'(A,3(f12.4))')    '### <Vel>              = ',MeanSource(i)%Velocity/NormSource(i)
       write(*,'(A,f12.4)')       '### <Erad> [GeV/A]     = ',MeanSource(i)%radEnergy/NormSource(i)
       write(*,'(A,f12.4)')       '### <Etot> [GeV/A]     = ',MeanSource(i)%ExEnergy/NormSource(i)
       if ( getRMF_Flag() ) then
          write(*,'(A,f12.4)')    '### <rho> [fm^-3]      = ',InvariantDensity(i)
          write(*,'(A,3f12.4,A)') '### <P_mean,zz,tr>     = ',Pmean(i),P_zz(i),P_tr(i),' [MeV*fm^-3]'
          write(*,'(A,f12.4)')    '### <Q>                = ',Anisotropy(i)
       end if
       !--------------------------------------------------------------------

       if ( getRMF_Flag() ) then

          if (.not.hyperSource) then
              write(106,101) time, &
                   &   float(MeanSource(i)%Size)/NormSource(i), &
                   &   float(MeanSource(i)%Charge)/NormSource(i), &
                   &   MeanSource(i)%Position/NormSource(i), &
                   &   MeanSource(i)%Velocity/NormSource(i), &
                   &   MeanSource(i)%radEnergy/NormSource(i), &
                   &   MeanSource(i)%ExEnergy/NormSource(i), &
                   &   InvariantDensity(i),Pmean(i),P_zz(i),P_tr(i),Anisotropy(i)
          else
              write(106,101) time, &
                   &   float(MeanSource(i)%Size)/NormSource(i), &
                   &   float(MeanSource(i)%Charge)/NormSource(i), &
                   &   float(MeanSource(i)%nLambda)/NormSource(i), &
                   &   float(MeanSource(i)%nSigma0)/NormSource(i), &
                   &   MeanSource(i)%Position/NormSource(i), &
                   &   MeanSource(i)%Velocity/NormSource(i), &
                   &   MeanSource(i)%radEnergy/NormSource(i), &
                   &   MeanSource(i)%ExEnergy/NormSource(i), &
                   &   InvariantDensity(i),Pmean(i),P_zz(i),P_tr(i),Anisotropy(i)
          end if


       else

          if (.not.hyperSource) then
              write(106,101) time, &
                   &   float(MeanSource(i)%Size)/NormSource(i), &
                   &   float(MeanSource(i)%Charge)/NormSource(i), &
                   &   MeanSource(i)%Position/NormSource(i), &
                   &   MeanSource(i)%Velocity/NormSource(i), &
                   &   MeanSource(i)%radEnergy/NormSource(i), &
                   &   MeanSource(i)%ExEnergy/NormSource(i)
          else
              write(106,101) time, &
                   &   float(MeanSource(i)%Size)/NormSource(i), &
                   &   float(MeanSource(i)%Charge)/NormSource(i), &
                   &   float(MeanSource(i)%nLambda)/NormSource(i), &
                   &   float(MeanSource(i)%nSigma0)/NormSource(i), &
                   &   MeanSource(i)%Position/NormSource(i), &
                   &   MeanSource(i)%Velocity/NormSource(i), &
                   &   MeanSource(i)%radEnergy/NormSource(i), &
                   &   MeanSource(i)%ExEnergy/NormSource(i)
          end if

       end if


    end do Loop_Sources

    close(106)

101 format(30f15.7)
100 format(f6.2,2(1x,i4),8(1x,f15.7),1x,f8.3,1x,i5,1x,i4)
102 format(f6.2,2(1x,i4),2(1x,i2),8(1x,f15.7),1x,f8.3,1x,i5,1x,i4)

  !****************************************************************************
  end subroutine WriteSourceInfo !****************************************
  !****************************************************************************

  !****************************************************************************
  !****s* sourceProperties/allocate_source
  ! NAME
  ! subroutine allocate_source
  !
  ! PURPOSE
  ! allocates source's properties
  ! NOTES
  ! The allocation of the source's properties has to be done at each
  ! time when this module is called from outside, because the size
  ! of the source properties is variable in time!
  !
  !****************************************************************************
  subroutine allocate_sourceFields(MaxNumSources)
    implicit none

    integer, intent(in) :: MaxNumSources

    allocate(PressureField(1:MaxNumSources,1:3))
    allocate(InvariantDensity(1:MaxNumSources))
    allocate(Anisotropy(1:MaxNumSources))

  !****************************************************************************
  end subroutine allocate_sourceFields !**********************************
  !****************************************************************************


  !****************************************************************************
  !****s* sourceProperties/deallocate_source
  ! NAME
  ! subroutine deallocate_source
  !
  ! PURPOSE
  ! deallocates source's properties.
  ! NOTES
  ! The size of the source properties is variable in time! see also
  ! notes in routine allocate_sourceFields.
  !
  !****************************************************************************
  subroutine deallocate_sourceFields
    implicit none

    if (allocated(PressureField)) deallocate(PressureField)
    if (allocated(InvariantDensity)) deallocate(InvariantDensity)
    if (allocated(Anisotropy)) deallocate(Anisotropy)

  !****************************************************************************
  end subroutine deallocate_sourceFields !********************************
  !****************************************************************************


end module sourceProperties
