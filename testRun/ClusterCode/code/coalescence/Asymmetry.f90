!***************************************************************************
!****m* /CheckAsymmetry
! NAME
! module CheckAsymmetry
! FUNCTION
! Check isospin asymmetry of produced clusters, before starting 
! evaporation procedure. Only clusters with correct N/Z-ratio 
! (those included in the nndc.dat table) are considered in the 
! evaporation procedure!
!***************************************************************************
module CheckAsymmetry

  PRIVATE

  PUBLIC :: GEMAsymmetry,FindNearestStableElement

  contains

    !***********************************************************************
    subroutine GEMAsymmetry(A,Z,CheckAsy)
    !***********************************************************************
      use InitGEM, only : NucPosition
      implicit  none
      !---------------------------------------------------------------------
      integer, intent(in)  :: A,Z
      logical, intent(out) :: CheckAsy
      integer, parameter   :: MaxN = 130
      integer, parameter   :: MaxZ = 90
      integer              :: N,Search
      !---------------------------------------------------------------------
      ! single particles (neutrons,protons) should NOT enter here
      !---------------------------------------------------------------------
      if (A==1) then
         write(*,*) 'from Module CheckAsymmetry, routine Asymmetry:'
         write(*,*) '!!!free proton/neutron as input in Asymmetry!!!'
         write(*,*) 'wrong input in Asymmetry! A,Z = ',A,Z
         write(*,*) 'Asymmetry check should be applied only to clusters!'
         write(*,*) '!!! TERMINATION OF PROGRAM !!!'
         STOP 
      endif
      !---------------------------------------------------------------------
      ! Clusters heavier than (N,Z)=(130,90) are not defined
      ! see also NOTES in InitGEM-module, routine Init_GEM
      !---------------------------------------------------------------------
      N = A-Z
      if (N > MaxN .or. Z > MaxZ) then
         write(*,*) 'from Module CheckAsymmetry, routine Asymmetry:'
         write(*,*) 'Cluster with N = ',N,' and Z = ',Z,' not defined!'
         write(*,*) 'Increase "MaxN" and "MaxZ" here and in InitGEM'
         write(*,*) '!!! TERMINATION OF PROGRAM !!!'
         STOP
      endif
      !---------------------------------------------------------------------
      Search   = NucPosition(N,Z)
      CheckAsy = .false.
      if (Search > 0) CheckAsy = .true.

    !***********************************************************************
    end subroutine GEMAsymmetry !*******************************************
    !***********************************************************************

    !***********************************************************************
    subroutine FindNearestStableElement(A,Z,NewN,NewZ)
    !***********************************************************************
      use random, only : rn
      use InitGEM, only : NucPosition,Elements
      implicit none
      !---------------------------------------------------------------------
      integer, intent(in)  :: A,Z
      integer, intent(out) :: NewN,NewZ
      integer :: N,Position,NewPosition,dummy,dummy2,NewE,NeutronEx,ProtonEx
      integer :: MassE
!      real    :: rand
      logical :: SearchWay,find
      !---------------------------------------------------------------------
      N = A-Z !neutron number
      !---------------------------------------------------------------------
      if (A.le.0 .or. N.le.0 .or. Z.le.0) then
         write(*,*) 'from module CheckAsymmetry, routine FindNearestStableElement: '
         write(*,*) 'Something wrong in input-variables in AsymmetryCheck...'
         write(*,*) 'Input-Variables: A,N,Z = ',A,N,Z
         write(*,*) '!!! TERMINATION OF PROGRAM NOW !!!'
         STOP
      endif
      !---------------------------------------------------------------------
      Position   = NucPosition(N,Z)
      if (Position > 0) then 
         write(*,*) 'from module CheckAsymmetry, routine FindNearestStableElement: '
         write(*,*) 'Something wrong in AsymmetryCheck...N,Z,Position = ',&
              &     N,Z,Position
         write(*,*) '!!! TERMINATION OF PROGRAM NOW !!!'
         STOP
      endif
      !---------------------------------------------------------------------
      if (N<Z) then
         dummy  = Z
         dummy2 = N
         SearchWay = .true.
      endif
      if (N>Z) then
         dummy  = N
         dummy2 = Z
         SearchWay = .false.
      endif
      if (N==Z) then
!         if (rand(0)<0.5) then
         if (rn()<0.5) then
            dummy  = Z
            dummy2 = N
            SearchWay = .true.
         else
            dummy  = Z
            dummy2 = N
            SearchWay = .false.
         endif
      endif
      !---------------------------------------------------------------------
      find = .false.
      do 
         NewE = dummy - 1
         if (SearchWay) then
            NewPosition = NucPosition(N,NewE)
            MassE = N + NewE
         else
            NewPosition = NucPosition(NewE,Z)
            MassE = Z + NewE
         endif
         if (NewPosition < 0) then
            dummy = NewE
         else
            find = .true.
            exit
         endif
         if (NewE==0 .or. MassE==1) then
            exit
         endif
      end do
      
      if (.not.find) then
         find = .false.
         do 
            NewE = dummy2 - 1
            if (SearchWay) then
               NewPosition = NucPosition(NewE,Z)
               MassE = Z + NewE
            else
               NewPosition = NucPosition(N,NewE)
               MassE = N + NewE
            endif
            if (NewPosition < 0) then
               dummy2 = NewE
            else
               find = .true.
               exit
            endif
            if (NewE==0 .or. MassE==1) then
               exit
            endif
         end do
      endif

      if (.not.find) then
         write(*,*) 'from Module CheckAsymmetry, routine FindNearestStableElement: '
         write(*,*) 'Search for nearest stable element not succefull...'
         write(*,*) 'this should not be happened!'
         write(*,*) '!!! TERMINATION OF THE PROGRAM NOW !!!'
         STOP         
      endif
      !---------------------------------------------------------------------
      !define the new values for Mass and Charge
      !set the proton and/or neutron excess to free particles
      !---------------------------------------------------------------------
      NewN = Elements(NewPosition)%ZahlN
      NewZ = Elements(NewPosition)%ZahlZ
      NeutronEx = N - NewN
      ProtonEx  = Z - NewZ
      if(NeutronEx < 0 .or. ProtonEx < 0) then
         write(*,*) 'from Module CheckAsymmetry, routine FindNearestStableElement: '
         write(*,*) 'something wrong in finding the nearest element...'
         write(*,*) 'NewPosition,NewN,NewZ = ',NewPosition,NewN,NewZ
         write(*,*) 'Old N,Old Z           = ',N,Z
         write(*,*) '!!! TERMINATION OF THE PROGRAM NOW !!!'
         STOP
      endif
    !***********************************************************************
    end subroutine FindNearestStableElement
    !***********************************************************************


  !*************************************************************************
  end module CheckAsymmetry !***********************************************
  !*************************************************************************
