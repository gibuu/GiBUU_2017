!***************************************************************************
!****m* /velocityFieldModule
! NAME
! module velocityFieldModule

! FUNCTION
! calculates the collective flow field on the grid points.
! NOTES
! The dimensions of the grid must be consistent with the GiBUU-run, 
! which is analyzed in the CLusterCode/code/smm directory in terms 
! of fragment production.
!***************************************************************************

Module velocityFieldModule

  implicit none
  PRIVATE

  PUBLIC :: Get_collectiveFlow, Get_VelocityFieldAt,& 
       & Get_RadialFlowProfile !,CelRatio


  !***********************************************************************
  !****g* velocityFieldModule/gridSize
  ! SOURCE
  !
  real,    dimension(1:3), save ::  gridSize=(/40.,40.,50./)
  !
  ! PURPOSE
  ! Size of 3D-grid in coordinate space (fm).
  !
  !***********************************************************************

  !***********************************************************************
  !****g* velocityFieldModule/gridPoints
  ! SOURCE
  !
  integer, dimension(1:3), save ::  gridPoints=(/40,40,50/) 
  !
  ! PURPOSE
  ! Number of gridpoints in each space direction
  !
  !***********************************************************************

  !***********************************************************************
  !****g* velocityFieldModule/gridSpacing
  ! SOURCE
  !
  real,    dimension(1:3), save ::  gridSpacing=1.
  !
  ! PURPOSE
  ! Spacing of the grid in each dimension.
  !
  !***********************************************************************

  !***********************************************************************
  !****g* velocityFieldModule/VelocityField
  ! SOURCE
  !
  Real, allocatable, dimension(:,:,:,:), SAVE :: VelocityField
  !
  ! PURPOSE
  ! The collective velocity vector: contains the collective flow 
  ! on the gridpoints for fireball sources.
  !
  !***********************************************************************

  !***********************************************************************
  !****g* velocityFieldModule/Norm
  ! SOURCE
  !
  Real, allocatable, dimension(:,:,:), SAVE :: Norm
  !
  ! PURPOSE
  ! Number of fireball particles at each grid cell.
  !***********************************************************************

  !***********************************************************************
  !****g* velocityFieldModule/Slope
  ! SOURCE
  !
  Real, dimension(1:3), SAVE :: Slope
  !
  ! PURPOSE
  ! Slope of linear radial flow profiles along x,y,z.
  !***********************************************************************

  !***********************************************************************
  !****g* velocityFieldModule/CelRatio
  ! SOURCE
  !
  !  Real, SAVE :: CelRatio=1.
  !
  ! PURPOSE
  ! Estimated axes of prolonged ellipsoidal fireball in transverse direction.
  !***********************************************************************

  contains

!***************************************************************************
  subroutine Get_collectiveFlow(teilchen)
!***************************************************************************
    use TypeDefinitions, only : particle

    type(particle), dimension(:,:), intent(in) :: teilchen

    integer :: i,j,i1,i2,i3

    real,    dimension(1:3) :: r
    integer, dimension(1:3) :: Index

    logical, SAVE :: firstCall=.true.

    !the option of an anisotropic fireball is not used any more
    !    integer :: normx,normy,normz

    if (firstCall) then
       call Get_GridParameters
       firstCall = .false.
    endif

    Norm(:,:,:)            = 0.0
    VelocityField(:,:,:,:) = 0.0

    do i=1,Size(teilchen,dim=1)
       do j=1,Size(teilchen,dim=2)

          if (teilchen(i,j)%ID .ge.100) cycle
!          if (teilchen(i,j)%sd == 999) cycle
          
          r(1:3)=teilchen(i,j)%position(1:3)
          Index(1:3) = nint(r(1:3)/gridspacing)

          if ( (abs(Index(1)).gt.gridSize(1)) .or.(abs(Index(2)).gt.gridSize(2)) .or. & 
               & ((abs(Index(3)).gt.gridSize(3))) ) then
             write(*,*) 'smmModule/Get_collectiveVelocity:'
             write(*,*) 'Increase grid!, Index(1:3) = ',Index(1:3),' STOP!'
             STOP
          endif
                   
          velocityField(1:3,Index(1),Index(2),Index(3)) = & 
               & velocityField(1:3,Index(1),Index(2),Index(3)) + & 
               & teilchen(i,j)%momentum(1:3)/teilchen(i,j)%momentum(0)
          Norm(Index(1),Index(2),Index(3)) = Norm(Index(1),Index(2),Index(3)) + 1.

       end do
    end do

    do I1 = -gridPoints(1),gridPoints(1)
       do I2 = -gridPoints(2),gridPoints(2)
          do I3 = -gridPoints(3),gridPoints(3)

             if (Norm(i1,i2,i3)==0.) then
                velocityField(:,i1,i2,i3) = 0.0
             else
                velocityField(:,i1,i2,i3) = velocityField(:,i1,i2,i3)/Norm(i1,i2,i3)
             endif

          end do
       end do
    end do   

!!$    do i1=-gridPoints(1),0
!!$       if (Norm(i1,0,0).ne.0.) then
!!$          normx = i1
!!$          exit
!!$       endif
!!$    end do
!!$    do i1=-gridPoints(2),0
!!$       if (Norm(0,i1,0).ne.0.) then
!!$          normy = i1
!!$          exit
!!$       endif
!!$    end do
!!$    do i1=-gridPoints(3),0
!!$       if (Norm(0,0,i1).ne.0.) then
!!$          normz = i1
!!$          exit
!!$       endif
!!$    end do
!!$
!!$    CelRatio = ( float(abs(normx))+float(abs(normy)) ) / 2. / float(abs(normz))
!***************************************************************************
    end subroutine Get_collectiveFlow !*************************************
!***************************************************************************

!***************************************************************************
    subroutine Get_RadialFlowProfile(isu)
!***************************************************************************

      integer, intent(in) :: isu
      integer :: i1
      real    :: x_points,y_points,slope_x1,slope_x2
      !---------------------------------------------------------------------
      !v_x(x)
      !---------------------------------------------------------------------
!      do i1=-gridPoints(1),gridPoints(1)
!         if (velocityField(1,i1,0,0).ne.0.) write(799+isu,1000) i1,velocityField(1,i1,0,0)
!      end do
      y_points = 0.0
      x_points = 0.0
      do i1=0,gridPoints(1)
         if (velocityField(1,i1,0,0).ne.0.) then
            y_points  = y_points + velocityField(1,i1,0,0)
            x_points  = x_points + float(i1)
         endif
      end do
      slope_x1 = y_points/x_points

      y_points = 0.0
      x_points = 0.0
      do i1=-gridPoints(1),0
         if (velocityField(1,i1,0,0).ne.0.) then
            y_points  = y_points + velocityField(1,i1,0,0)
            x_points  = x_points + float(i1)
         endif
      end do
      slope_x2 = y_points/x_points

      slope(1) = ( slope_x1 + slope_x2) / 2.

      !---------------------------------------------------------------------
      !v_y(y)
      !---------------------------------------------------------------------
!      do i1=-gridPoints(2),gridPoints(2)
!         if (velocityField(2,0,i1,0).ne.0.) write(849+isu,1000) i1,velocityField(2,0,i1,0)
!      end do
      y_points = 0.0
      x_points = 0.0
      do i1=0,gridPoints(2)
         if (velocityField(2,0,i1,0).ne.0.) then
            y_points  = y_points + velocityField(2,0,i1,0)
            x_points  = x_points + float(i1)
         endif
      end do
      slope_x1 = y_points/x_points

      y_points = 0.0
      x_points = 0.0
      do i1=-gridPoints(2),0
         if (velocityField(2,0,i1,0).ne.0.) then
            y_points  = y_points + velocityField(2,0,i1,0)
            x_points  = x_points + float(i1)
         endif
      end do
      slope_x2 = y_points/x_points

      slope(2) = ( slope_x1 + slope_x2) / 2.

      !---------------------------------------------------------------------
      !v_z(z)
      !---------------------------------------------------------------------
!      do i1=-gridPoints(3),gridPoints(3)
!         if (velocityField(3,0,0,i1).ne.0.) write(899+isu,1000) i1,velocityField(3,0,0,i1)
!      end do
      y_points = 0.0
      x_points = 0.0
      do i1=0,gridPoints(3)
         if (velocityField(3,0,0,i1).ne.0.) then
            y_points  = y_points + velocityField(3,0,0,i1)
            x_points  = x_points + float(i1)
         endif
      end do
      slope_x1 = y_points/x_points

      y_points = 0.0
      x_points = 0.0
      do i1=-gridPoints(3),0
         if (velocityField(3,0,0,i1).ne.0.) then
            y_points  = y_points + velocityField(3,0,0,i1)
            x_points  = x_points + float(i1)
         endif
      end do
      slope_x2 = y_points/x_points

      slope(3) = ( slope_x1 + slope_x2) / 2.

      write(*,'(A,3f9.4)') 'slopes(x,y,z) = ',slope
! 1000  format(i4,2x,f8.4)

!***************************************************************************
    end subroutine Get_RadialFlowProfile !**********************************
!***************************************************************************


!***************************************************************************
    subroutine Get_GridParameters
!***************************************************************************

      integer :: ios
      !---------------------------------------------------------------------
      NAMELIST /inputGrid/GridSize,GridPoints

      write(*,*) '----- Reading Namelist "inputGrid": Start -----'
      rewind(5)
      read(5,nml=inputGrid,IOSTAT=ios)
      !---------------------------------------------------------------------
      if(ios.ne.0) then
         write(*,*) 'Error in  namelist "inputGrid" : This namelist is crucial. STOP!'
         stop
      end if
      gridSpacing=gridSize/float(gridPoints)

      write(*,'(A,2(F5.2,","),F5.2,A)') & 
           & '  The gridsize of the density grid is        = (', gridsize, ') fm'
      write(*,'(A,2(I5,","),I5,A)')     & 
           & '  The number of gridpoints per dimension are = (', gridPoints, ') '
      write(*,'(A,2(F5.2,","),F5.2,A)') & 
           & '  The grid spacing is                        = (', gridSpacing, ') fm'

      allocate(VelocityField(1:3,& 
           & -gridPoints(1):gridPoints(1),-gridPoints(2):gridPoints(2),&
           & -gridPoints(3):gridPoints(3)))
      allocate(Norm(-gridPoints(1):gridPoints(1),-gridPoints(2):gridPoints(2),&
           &        -gridPoints(3):gridPoints(3)))

      write(*,*) '----- Reading Namelist "inputGrid": END   -----'
      write(*,*)

!***************************************************************************
    end subroutine Get_GridParameters !*************************************
!***************************************************************************

!***************************************************************************
    subroutine Get_VelocityFieldAt(r,velo)
!***************************************************************************

      Real, dimension(1:3), intent(in)  :: r
      Real, dimension(1:3), intent(out) :: velo
      
      integer, dimension(1:3) :: ipos

!      real :: meanSlope

      ipos = nint(r/gridSpacing)
      if ( (abs(iPos(1))>gridSize(1)).or.(abs(iPos(2))>gridSize(2)).or.& 
           & (abs(iPos(3))>gridSize(3)) ) then
         write(*,*) 'upss....increase grid....'
         STOP
      endif

!      velo(1:3) = VelocityField(1:3,ipos(1),ipos(2),ipos(3))
      velo(1:3) = slope(1:3)*r(1:3)

!      velo_sp(1:3) = VelocityField_sp(1:3,ipos(1),ipos(2),ipos(3))
!      meanSlope = ( slope(1) + slope(2) ) / 2.
!      velo(1:3) = meanSlope*r
!***************************************************************************
    end subroutine Get_VelocityFieldAt !************************************
!***************************************************************************


end module velocityFieldModule
