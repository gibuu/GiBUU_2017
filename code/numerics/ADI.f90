!******************************************************************************
!****m* /ADI
! NAME
! module ADI
!
! PURPOSE
! This module contains the subroutine ADI_solve for the solution of the matrix
! equations, which appear when the gradient terms are taken into account in RMF.
! Can also be applied for the Yukawa or Coulomb potential calculations.
! Uses the Alternating Direction Implicit (ADI) iterative method, similar
! to the one from module yukawa.
! (See S.Teis, PhD thesis, and the book of R.S.Varga,
! "Matrix Iterative Analysis".)
!******************************************************************************
module ADI

  implicit none

  private

  !****************************************************************************
  !****g* ADI/debugADI
  ! SOURCE
  !
  logical, parameter :: debugADI=.false.
  ! PURPOSE
  ! Switch on some debug output.
  !****************************************************************************

  public :: ADI_solve
  public :: ADI_solve_Douglas, ADI_Coulomb, cleanUp

  real, save, dimension(:), allocatable :: alpx, betx, &  ! coefficients needed
       alpy, bety, &  ! for the tridiagonal inversion
       alpz, betz

  real, save, dimension(:,:,:), allocatable :: U1, U2, U3 ! intermediate iterations of the field U

  integer, save :: Nx, Ny, Nz

  logical, save :: flagini=.true.

contains

  !****************************************************************************
  !****s* ADI/init
  ! NAME
  ! subroutine init
  ! PURPOSE
  ! Reads input switches. Initializes fields.
  !****************************************************************************
  subroutine init
!!$    use output
!!$
!!$    integer :: ios
!!$    !*************************************************************************
!!$    !****n* ADI/ADI_input
!!$    ! NAME
!!$    ! NAMELIST /ADI_input/
!!$    ! PURPOSE
!!$    ! Includes the switches:
!!$    ! * debugADI
!!$    !*************************************************************************
!!$    NAMELIST /ADI_input/ debugADI
!!$
!!$    call Write_ReadingInput('ADI_input',0)
!!$    rewind(5)
!!$    read(5,nml=ADI_input,iostat=ios)
!!$    call Write_ReadingInput('ADI_input',0,ios)
!!$    write(*,*) 'Set debugADI to', debugADI
!!$    call Write_ReadingInput('ADI_input',1)
  end subroutine init

  subroutine cleanUp
    if (allocated(alpx)) deallocate(alpx)
    if (allocated(betx)) deallocate(betx)
    if (allocated(alpy)) deallocate(alpy)
    if (allocated(bety)) deallocate(bety)
    if (allocated(alpz)) deallocate(alpz)
    if (allocated(betz)) deallocate(betz)
    if (allocated(U1)) deallocate(U1)
    if (allocated(U2)) deallocate(U2)
    if (allocated(U3)) deallocate(U3)
  end subroutine cleanUp

  !****************************************************************************
  !****s* ADI/ADI_solve_Douglas
  ! NAME
  ! subroutine ADI_solve_Douglas(U, S, Diag, eta_x, eta_y, CoulombFlag, rconv, niter_out, error_out)
  ! PURPOSE
  ! Solve the linear system
  !   ( tildeX + tildeY + tildeZ )* U = S
  ! with respect to the field U(i,j,k), i=2,...,Nx-1, j=2,...,Ny-1, k=2,...,Nz-1,
  ! where Nx=2*gridPoints(1)+1, Ny=2*gridPoints(2)+1,  Nz=2*gridPoints(3)+1.
  ! Here tildeX, tildeY and tildeZ are the linear operators:
  !   [ tildeX * U ](i,j,k) = ( - U(i-1,j,k) + 2.*U(i,j,k) -  U(i+1,j,k) )*eta_x
  !                           + Diag(i,j,k)*U(i,j,k),
  !   [ tildeY * U ](i,j,k) = ( - U(i,j-1,k) + 2.*U(i,j,k) -  U(i,j+1,k) )*eta_y
  !                           + Diag(i,j,k)*U(i,j,k),
  !   [ tildeZ * U ](i,j,k) = - U(i,j,k-1) + 2.*U(i,j,k) -  U(i,j,k+1)
  !                           + Diag(i,j,k)*U(i,j,k).
  ! with:
  ! * S(i,j,k), i=2,...,Nx-1, j=2,...,Ny-1, k=2,...,Nz-1 --- the source field.
  ! The values of U(i,j,k) at the boundary ( i=1,Nx or j=1,Ny  or k=1,Nz ) are fixed.
  !
  ! INPUTS
  ! * real, dimension(:,:,:), intent(inout) :: U    ! -- starting value of the field U for iterations,
  ! * real, dimension(:,:,:), intent(in)    :: S    ! -- the source field,
  ! * real, dimension(:,:,:), intent(in)    :: Diag ! -- diagonal coefficient field,
  ! * real, intent(in) :: eta_x          ! =(gridSpacing(3)/gridSpacing(1))**2,
  ! * real, intent(in) :: eta_y          ! =(gridSpacing(3)/gridSpacing(2))**2
  ! * logical, optional, intent(in) :: CoulombFlag   ! please, set .true. when Coulomb field is computed
  ! * real, optional, intent(in)    :: rconv         ! convergence parameter
  ! * integer, optional, intent(out) :: niter_out    ! number of iterations
  ! * real, optional, intent(out) :: error_out       ! error
  ! OUTPUT
  ! * real, dimension(:,:,:), intent(inout) :: U    ! -- final iterated value of the field U,
  ! NOTES
  ! Uses the Douglas Iterative Method
  ! (c.f. R.S.Varga, "Matrix Iterative Analysis")
  ! The tridiagonal inversion procedure is explained, e.g. in
  ! S.E. Koonin, D.C. Meredith, "Computational Physics"
  !****************************************************************************
  subroutine ADI_solve_Douglas(U, S, Diag, eta_x, eta_y, CoulombFlag, rconv, niter_out, error_out)

    use constants, only: pi

    real, dimension(:,:,:), intent(inout) :: U
    real, dimension(:,:,:), intent(in)    :: S
    real, dimension(:,:,:), intent(in)    :: Diag
    real, intent(in)    :: eta_x
    real, intent(in)    :: eta_y

    logical, optional, intent(in) :: CoulombFlag
    real, optional, intent(in)    :: rconv
    integer, optional, intent(out) :: niter_out
    real, optional, intent(out) :: error_out

    real, save :: X_min, Y_min, Z_min, XYZ_min, XYZ_max
    real :: Diag_min, Diag_max, absS_max
    real :: alpha, beta , r !, gamma
    integer :: j, ix, iy, iz
    real :: A0x, A0y, A0z, A0, b, denom, error, lhs, error_overall
    real :: cx, cy, cz, cdia

    integer, SAVE :: numberIterations_max=20       ! Maximum number of iterations
    real,    SAVE :: relativeError_max=1.e-04      ! Maximum relative error required

    if (flagini) then
       call init
       Nx=size(U,dim=1)
       Ny=size(U,dim=2)
       Nz=size(U,dim=3)
       if (debugADI) write(*,'(A,3i5)') 'In ADI_solve_Douglas, Nx, Ny, Nz:', Nx, Ny, Nz
       allocate(alpx(1:Nx-1))
       allocate(betx(1:Nx-1))
       allocate(alpy(1:Ny-1))
       allocate(bety(1:Ny-1))
       allocate(alpz(1:Nz-1))
       allocate(betz(1:Nz-1))
       allocate(U1(1:Nx,1:Ny,1:Nz))
       allocate(U2(1:Nx,1:Ny,1:Nz))
       allocate(U3(1:Nx,1:Ny,1:Nz))
       X_min=4.*eta_x*sin(pi/(2.*(Nx-1)))**2
       Y_min=4.*eta_y*sin(pi/(2.*(Ny-1)))**2
       Z_min=4.*sin(pi/(2.*(Nz-1)))**2
       if (debugADI) write(*,*) 'In ADI_solve_Douglas, X_min, Y_min, Z_min:', X_min, Y_min, Z_min
       XYZ_min=min(X_min,Y_min,Z_min)
       XYZ_max=4.*max(eta_x,eta_y,1.)
       !       XYZ_min=0.33
       !       XYZ_max=5.7
       flagini=.false.
    end if

    if (present(CoulombFlag)) then
       numberIterations_max = 600
       relativeError_max=3.e-04
    else
       numberIterations_max = 20
       relativeError_max=1.e-04
    end if

    if (debugADI) write(*,'(A,i5,2x,e12.5)') 'In ADI_solve_Douglas, numberIterations_max,relativeError_max: ',&
         & numberIterations_max,relativeError_max

    if (present(rconv)) then
       r=rconv
    else if (present(CoulombFlag)) then
       r=0.5   ! Grid step 1 fm
       !r=0.3   ! Grid step 0.5 fm
    else

       Diag_min=1.e+06
       Diag_max=-1.e+06
       do ix = 2,Nx-1
          do iy = 2,Ny-1
             do iz = 2,Nz-1
                if (Diag(ix,iy,iz).lt.Diag_min) Diag_min=Diag(ix,iy,iz)
                if (Diag(ix,iy,iz).gt.Diag_max) Diag_max=Diag(ix,iy,iz)
             end do
          end do
       end do

       if (Diag_min.lt.0.) write(*,*) ' Warning in ADI_solve_Douglas: Diag_min < 0'

       ! Estimate lower (alpha) and upper (beta) boundaries of the eigenvalues
       ! of the matrices  tildeX, tildeY and tildeZ:

       alpha = XYZ_min + Diag_min
       beta = XYZ_max + Diag_max

       if (alpha.lt.0.) then
          write(*,*) ' In ADI_solve_Douglas: alpha < 0'
          stop
       end if

       !gamma=(beta/alpha)**(0.5/float(numberIterations_Max))

       if (debugADI) write(*,*) 'In ADI_solve_Douglas, alpha, beta:', alpha, beta

       r=sqrt(alpha*beta)

    end if

    absS_max=0.
    do ix = 2,Nx-1
       do iy = 2,Ny-1
          do iz = 2,Nz-1
             if (abs(S(ix,iy,iz)).gt.absS_max) absS_max=abs(S(ix,iy,iz))
          end do
       end do
    end do

    if (debugADI) write(*,*) ' In ADI_solve_Douglas: absS_max = ', absS_max

    ! Copying in oder to have the correct boundary conditions also for U1, U2 and U3:
    U1=U
    U2=U
    U3=U

    alpx(Nx-1)=0.
    alpy(Ny-1)=0.
    alpz(Nz-1)=0.

    j=0
    Loop_over_iterations : do

       j=j+1

       !r=alpha*gamma**(2.*float(j)-1.)

       if (debugADI) write(*,*) 'In ADI_solve_Douglas, r:', r

       A0x= 2. + r/eta_x
       cx= r - 2.*eta_x - 4.*eta_y - 4.

       do iy = 2,Ny-1
          do iz = 2,Nz-1

             ! Inversion along x-axis for fixed y and z:

             betx(Nx-1) = U1(Nx,iy,iz)

             do ix = Nx-1,2,-1

                A0= A0x + Diag(ix,iy,iz) / eta_x

                b = (  ( cx - 5.*Diag(ix,iy,iz) ) * U(ix,iy,iz) &
                     &+      eta_x * ( U(ix-1,iy,iz) + U(ix+1,iy,iz) ) &
                     &+ 2. * ( eta_y * ( U(ix,iy-1,iz) + U(ix,iy+1,iz) ) &
                     &+        U(ix,iy,iz-1) + U(ix,iy,iz+1) + S(ix,iy,iz) ) ) / eta_x

                denom = A0 - alpx(ix)
                if (denom.eq.0.) then
                   write(*,*) 'In ADI, denom_x=0.'
                   stop
                end if

                alpx(ix-1) = 1. / denom
                betx(ix-1) = ( betx(ix) + b ) / denom
             end do

             do ix = 1,Nx-2
                U1(ix+1,iy,iz) = alpx(ix) *  U1(ix,iy,iz) + betx(ix)
             end do

          end do
       end do

       A0y= 2. + r/eta_y
       cy=  r - 2. * ( eta_x + eta_y + 2. )

       do ix = 2,Nx-1
          do iz = 2,Nz-1

             ! Inversion along y-axis for fixed x and z:
             bety(Ny-1) = U2(ix,Ny,iz)

             do iy = Ny-1,2,-1

                A0 = A0y + Diag(ix,iy,iz) / eta_y

                b = ( ( cy - 4.*Diag(ix,iy,iz) ) * U(ix,iy,iz) &
                     &+ eta_x * ( U(ix-1,iy,iz) + U(ix+1,iy,iz) + U1(ix-1,iy,iz) + U1(ix+1,iy,iz) ) &
                     &+ eta_y * ( U(ix,iy-1,iz) + U(ix,iy+1,iz) ) &
                     &+    2. * ( U(ix,iy,iz-1) + U(ix,iy,iz+1) + S(ix,iy,iz) ) &
                     &- ( 2.*eta_x + Diag(ix,iy,iz) ) * U1(ix,iy,iz) ) / eta_y


                denom = A0 - alpy(iy)
                if (denom.eq.0.) then
                   write(*,*) 'In ADI_solve_Douglas, denom_y=0.'
                   stop
                end if

                alpy(iy-1) = 1. / denom
                bety(iy-1) = ( bety(iy) + b ) / denom
             end do

             do iy = 1,Ny-2
                U2(ix,iy+1,iz) = alpy(iy) *  U2(ix,iy,iz) + bety(iy)
             end do

          end do
       end do


       A0z= 2. + r
       cz= r - 2. * ( eta_x + eta_y + 1. )

       do ix = 2,Nx-1
          do iy = 2,Ny-1

             ! Inversion along z-axis for fixed x and y:
             betz(Nz-1) = U3(ix,iy,Nz)

             do iz = Nz-1,2,-1

                A0= A0z + Diag(ix,iy,iz)

                b = ( cz - 3. * Diag(ix,iy,iz) ) * U(ix,iy,iz) &
                     &+ eta_x * ( U(ix-1,iy,iz) + U(ix+1,iy,iz) + U1(ix-1,iy,iz) + U1(ix+1,iy,iz) ) &
                     &+ eta_y * ( U(ix,iy-1,iz) + U(ix,iy+1,iz) + U2(ix,iy-1,iz) + U2(ix,iy+1,iz) ) &
                     &+           U(ix,iy,iz-1) + U(ix,iy,iz+1) &
                     &- ( 2.*eta_x + Diag(ix,iy,iz) ) * U1(ix,iy,iz)  &
                     &- ( 2.*eta_y + Diag(ix,iy,iz) ) * U2(ix,iy,iz)  + 2. * S(ix,iy,iz)

                denom = A0 - alpz(iz)
                if (denom.eq.0.) then
                   write(*,*) 'In ADI_solve_Douglas, denom_z=0.'
                   stop
                end if

                alpz(iz-1) = 1. / denom
                betz(iz-1) = ( betz(iz) + b ) / denom
             end do

             do iz = 1,Nz-2
                U3(ix,iy,iz+1) = alpz(iz) *  U3(ix,iy,iz) + betz(iz)
             end do

          end do
       end do

       ! Copy final iteration to the ouput field U:
       U=U3


       ! Check  the error:
       error_overall = 0.
       cdia= 2.*(eta_x+eta_y+1.)

       do ix = 2,Nx-1
          do iy = 2,Ny-1
             do iz = 2,Nz-1

                lhs =  -eta_x * ( U(ix-1,iy,iz) + U(ix+1,iy,iz) ) &
                     & -eta_y * ( U(ix,iy-1,iz) + U(ix,iy+1,iz) ) &
                     & -U(ix,iy,iz-1) - U(ix,iy,iz+1) &
                     & + ( cdia + 3.*Diag(ix,iy,iz) ) * U(ix,iy,iz)

                error = abs(lhs - S(ix,iy,iz))
                if (error.gt.error_overall) error_overall=error

             end do
          end do
       end do

       if (absS_max.gt.1.e-10) error_overall=error_overall/absS_max

       if (debugADI) write(*,*) 'In ADI_solve_Douglas, iteration, error:', j, error_overall

       if ( error_overall .lt. relativeError_max ) then
          exit
       else if ( j .eq. numberIterations_max ) then
          write(*,*) ' Warning in ADI_solve_Douglas: exit after', numberIterations_max,' iterations'
          write(*,*) ' Accuracy reached: ', error_overall
          exit
       end if

    end do Loop_over_iterations

    if (present(niter_out)) niter_out=j
    if (present(error_out)) error_out=error_overall

  end subroutine ADI_solve_Douglas

  !****************************************************************************
  !****s* ADI/ADI_Coulomb
  ! NAME
  ! subroutine ADI_Coulomb(U, S, eta_x, eta_y, rconv, niter_out, error_out)
  ! PURPOSE
  ! Solve the linear system
  !   ( tildeX + tildeY + tildeZ )* U = S
  ! with respect to the field U(i,j,k), i=2,...,Nx-1, j=2,...,Ny-1, k=2,...,Nz-1,
  ! where Nx=2*gridPoints(1)+1, Ny=2*gridPoints(2)+1,  Nz=2*gridPoints(3)+1.
  ! Here tildeX, tildeY and tildeZ are the linear operators:
  !   [ tildeX * U ](i,j,k) = ( - U(i-1,j,k) + 2.*U(i,j,k) -  U(i+1,j,k) )*eta_x
  !   [ tildeY * U ](i,j,k) = ( - U(i,j-1,k) + 2.*U(i,j,k) -  U(i,j+1,k) )*eta_y
  !   [ tildeZ * U ](i,j,k) = - U(i,j,k-1) + 2.*U(i,j,k) -  U(i,j,k+1)
  ! with:
  ! * S(i,j,k), i=2,...,Nx-1, j=2,...,Ny-1, k=2,...,Nz-1 --- the source field.
  ! The values of U(i,j,k) at the boundary ( i=1,Nx or j=1,Ny  or k=1,Nz ) are fixed.
  ! INPUTS
  ! * real, dimension(:,:,:), intent(inout) :: U    ! -- starting value of the field U for iterations,
  ! * real, dimension(:,:,:), intent(in)    :: S    ! -- the source field,
  ! * real, intent(in) :: eta_x          ! =(gridSpacing(3)/gridSpacing(1))**2,
  ! * real, intent(in) :: eta_y          ! =(gridSpacing(3)/gridSpacing(2))**2
  ! * real, optional, intent(in)    :: rconv         ! convergence parameter
  ! * integer, optional, intent(out) :: niter_out    ! number of iterations
  ! * real, optional, intent(out) :: error_out       ! error
  ! OUTPUT
  ! * real, dimension(:,:,:), intent(inout) :: U    ! -- final iterated value of the field U,
  ! NOTES
  ! Uses the Douglas Iterative Method (c.f. R.S. Varga "Matrix Iterative Analysis")
  ! The tridiagonal inversion procedure is explained, e.g. in S.E. Koonin, D.C. Meredith,
  ! "Computational Physics". This subroutine is the same as ADI_solve_Douglas, but
  ! without diagonal coefficients.
  !****************************************************************************
  subroutine ADI_Coulomb(U, S, eta_x, eta_y, rconv, niter_out, error_out)

    real, dimension(:,:,:), intent(inout) :: U
    real, dimension(:,:,:), intent(in)    :: S
    real, intent(in)    :: eta_x
    real, intent(in)    :: eta_y

    real, optional, intent(in)    :: rconv
    integer, optional, intent(out) :: niter_out
    real, optional, intent(out) :: error_out

    real :: absS_max
    real :: r
    integer :: j, ix, iy, iz
    real :: A0x, A0y, A0z, b, denom, error, lhs, error_overall
    real :: cx, cy, cz, cdia

    integer, SAVE :: numberIterations_max=20       ! Maximum number of iterations
    real,    SAVE :: relativeError_max=1.e-04      ! Maximum relative error required

    if (flagini) then
       call init
       Nx=size(U,dim=1)
       Ny=size(U,dim=2)
       Nz=size(U,dim=3)
       if (debugADI) write(*,'(A,3i5)') 'In ADI_Coulomb, Nx, Ny, Nz:', Nx, Ny, Nz
       allocate(alpx(1:Nx-1))
       allocate(betx(1:Nx-1))
       allocate(alpy(1:Ny-1))
       allocate(bety(1:Ny-1))
       allocate(alpz(1:Nz-1))
       allocate(betz(1:Nz-1))
       allocate(U1(1:Nx,1:Ny,1:Nz))
       allocate(U2(1:Nx,1:Ny,1:Nz))
       allocate(U3(1:Nx,1:Ny,1:Nz))
       flagini=.false.
    end if

    numberIterations_max = 600
    relativeError_max=3.e-04

    if (debugADI) write(*,'(A,i5,2x,e12.5)') 'In ADI_Coulomb, numberIterations_max,relativeError_max: ',&
         & numberIterations_max,relativeError_max

    if (present(rconv)) then
       r=rconv
    else
       r=0.5   ! Grid step 1 fm
       !r=0.3   ! Grid step 0.5 fm
    end if

    absS_max=0.
    do ix = 2,Nx-1
       do iy = 2,Ny-1
          do iz = 2,Nz-1
             if (abs(S(ix,iy,iz)).gt.absS_max) absS_max=abs(S(ix,iy,iz))
          end do
       end do
    end do

    if (debugADI) write(*,*) ' In ADI_Coulomb: absS_max = ', absS_max

    ! Copying in oder to have the correct boundary conditions also for U1, U2 and U3:
    U1=U
    U2=U
    U3=U

    alpx(Nx-1)=0.
    alpy(Ny-1)=0.
    alpz(Nz-1)=0.

    A0x= 2. + r/eta_x
    A0y= 2. + r/eta_y
    A0z= 2. + r

    cx= r - 2.*eta_x - 4.*eta_y - 4.
    cy=  r - 2. * ( eta_x + eta_y + 2. )
    cz= r - 2. * ( eta_x + eta_y + 1. )
    cdia= 2.*(eta_x+eta_y+1.)

    j=0
    Loop_over_iterations : do

       j=j+1

       if (debugADI) write(*,*) 'In ADI_Coulomb, r:', r

       do iy = 2,Ny-1
          do iz = 2,Nz-1

             ! Inversion along x-axis for fixed y and z:

             betx(Nx-1) = U1(Nx,iy,iz)

             do ix = Nx-1,2,-1

                b = (  cx * U(ix,iy,iz) &
                     &+      eta_x * ( U(ix-1,iy,iz) + U(ix+1,iy,iz) ) &
                     &+ 2. * ( eta_y * ( U(ix,iy-1,iz) + U(ix,iy+1,iz) ) &
                     &        + U(ix,iy,iz-1) + U(ix,iy,iz+1) + S(ix,iy,iz) ) ) / eta_x

                denom = A0x - alpx(ix)
                if (denom.eq.0.) then
                   write(*,*) 'In ADI_Coulomb, denom_x=0.'
                   stop
                end if

                alpx(ix-1) = 1. / denom
                betx(ix-1) = ( betx(ix) + b ) / denom
             end do

             do ix = 1,Nx-2
                U1(ix+1,iy,iz) = alpx(ix) *  U1(ix,iy,iz) + betx(ix)
             end do

          end do
       end do


       do ix = 2,Nx-1
          do iz = 2,Nz-1

             ! Inversion along y-axis for fixed x and z:
             bety(Ny-1) = U2(ix,Ny,iz)

             do iy = Ny-1,2,-1

                b = (  cy * U(ix,iy,iz) &
                     &+ eta_x * ( U(ix-1,iy,iz) + U(ix+1,iy,iz) + U1(ix-1,iy,iz) + U1(ix+1,iy,iz) ) &
                     &+ eta_y * ( U(ix,iy-1,iz) + U(ix,iy+1,iz) ) &
                     &+    2. * ( U(ix,iy,iz-1) + U(ix,iy,iz+1)  &
                     &           - eta_x * U1(ix,iy,iz) + S(ix,iy,iz) ) ) / eta_y

                denom = A0y - alpy(iy)
                if (denom.eq.0.) then
                   write(*,*) 'In ADI_Coulomb, denom_y=0.'
                   stop
                end if

                alpy(iy-1) = 1. / denom
                bety(iy-1) = ( bety(iy) + b ) / denom
             end do

             do iy = 1,Ny-2
                U2(ix,iy+1,iz) = alpy(iy) *  U2(ix,iy,iz) + bety(iy)
             end do

          end do
       end do


       do ix = 2,Nx-1
          do iy = 2,Ny-1

             ! Inversion along z-axis for fixed x and y:
             betz(Nz-1) = U3(ix,iy,Nz)

             do iz = Nz-1,2,-1

                b = cz * U(ix,iy,iz) &
                     &+ eta_x * ( U(ix-1,iy,iz) + U(ix+1,iy,iz) + U1(ix-1,iy,iz) + U1(ix+1,iy,iz) ) &
                     &+ eta_y * ( U(ix,iy-1,iz) + U(ix,iy+1,iz) + U2(ix,iy-1,iz) + U2(ix,iy+1,iz) ) &
                     &+           U(ix,iy,iz-1) + U(ix,iy,iz+1) &
                     &- 2.*( eta_x * U1(ix,iy,iz)  &
                     &      +eta_y * U2(ix,iy,iz)  - S(ix,iy,iz) )

                denom = A0z - alpz(iz)
                if (denom.eq.0.) then
                   write(*,*) 'In ADI_Coulomb, denom_z=0.'
                   stop
                end if

                alpz(iz-1) = 1. / denom
                betz(iz-1) = ( betz(iz) + b ) / denom
             end do

             do iz = 1,Nz-2
                U3(ix,iy,iz+1) = alpz(iz) *  U3(ix,iy,iz) + betz(iz)
             end do

          end do
       end do

       ! Copy final iteration to the ouput field U:
       U=U3


       ! Check  the error:
       error_overall = 0.

       do ix = 2,Nx-1
          do iy = 2,Ny-1
             do iz = 2,Nz-1

                lhs =  -eta_x * ( U(ix-1,iy,iz) + U(ix+1,iy,iz) ) &
                     & -eta_y * ( U(ix,iy-1,iz) + U(ix,iy+1,iz) ) &
                     & -U(ix,iy,iz-1) - U(ix,iy,iz+1) &
                     & + cdia * U(ix,iy,iz)

                error = abs(lhs - S(ix,iy,iz))
                if (error.gt.error_overall) error_overall=error

             end do
          end do
       end do

       if (absS_max.gt.1.e-10) error_overall=error_overall/absS_max

       if (debugADI) write(*,*) 'In ADI_Coulomb, iteration, error:', j, error_overall

       if ( error_overall .lt. relativeError_max ) then
          exit
       else if ( j .eq. numberIterations_max ) then
          write(*,*) ' Warning in ADI_Coulomb: exit after', numberIterations_max,' iterations'
          write(*,*) ' Accuracy reached: ', error_overall
          exit
       end if

    end do Loop_over_iterations

    if (present(niter_out)) niter_out=j
    if (present(error_out)) error_out=error_overall

  end subroutine ADI_Coulomb

  !****************************************************************************
  !****s* ADI/ADI_solve
  ! NAME
  ! subroutine ADI_solve(U, S, Diag, eta_x, eta_y, CoulombFlag, rconv, niter_out, error_out)
  ! PURPOSE
  ! Solve the linear system
  !   ( tildeX + tildeY + tildeZ )* U = S
  ! with respect to the field U(i,j,k), i=2,...,Nx-1, j=2,...,Ny-1, k=2,...,Nz-1,
  ! where Nx=2*gridPoints(1)+1, Ny=2*gridPoints(2)+1,  Nz=2*gridPoints(3)+1.
  ! Here tildeX, tildeY and tildeZ are the linear operators:
  !   [ tildeX * U ](i,j,k) = ( - U(i-1,j,k) + 2.*U(i,j,k) -  U(i+1,j,k) )*eta_x
  !                         + Diag(i,j,k)*U(i,j,k),
  !   [ tildeY * U ](i,j,k) = ( - U(i,j-1,k) + 2.*U(i,j,k) -  U(i,j+1,k) )*eta_y
  !                         + Diag(i,j,k)*U(i,j,k),
  !   [ tildeZ * U ](i,j,k) = - U(i,j,k-1) + 2.*U(i,j,k) -  U(i,j,k+1)
  !                         + Diag(i,j,k)*U(i,j,k).
  ! with:
  ! * S(i,j,k), i=2,...,Nx-1, j=2,...,Ny-1, k=2,...,Nz-1 --- the source field.
  ! The values of U(i,j,k) at the boundary ( i=1,Nx or j=1,Ny  or k=1,Nz ) are fixed.
  ! INPUTS
  ! * real, dimension(:,:,:), intent(inout) :: U    ! -- starting value of the field U for iterations,
  ! * real, dimension(:,:,:), intent(in)    :: S    ! -- the source field,
  ! * real, dimension(:,:,:), intent(in)    :: Diag ! -- diagonal coefficient field,
  ! * real, intent(in)    :: Diag_min, Diag_max     ! -- min and max values of Diag,
  ! * real, intent(in) :: eta_x          ! =(gridSpacing(3)/gridSpacing(1))**2,
  ! * real, intent(in) :: eta_y          ! =(gridSpacing(3)/gridSpacing(2))**2
  ! * logical, optional, intent(in) :: CoulombFlag   ! please, set .true. when Coulomb field is computed
  ! * real, optional, intent(in)    :: rconv         ! convergence parameter
  ! * integer, optional, intent(out) :: niter_out    ! number of iterations
  ! * real, optional, intent(out) :: error_out       ! error
  ! OUTPUT
  ! * real, dimension(:,:,:), intent(inout) :: U    ! -- final iterated value of the field U,
  ! * real, intent(out)                     :: funMax ! =max(abs(( tildeX + tildeY + tildeZ )* U - S))
  ! NOTES
  ! Uses the Peaceman-Rachford Iterative Method (c.f. R.S. Varga "Matrix Iterative Analysis")
  ! The tridiagonal inversion procedure is explained, e.g. in S.E. Koonin, D.C. Meredith,
  ! "Computational Physics"
  !****************************************************************************
  subroutine ADI_solve(U, S, Diag, eta_x, eta_y, CoulombFlag, rconv, niter_out, error_out)

    use constants, only: pi

    real, dimension(:,:,:), intent(inout) :: U
    real, dimension(:,:,:), intent(in)    :: S
    real, dimension(:,:,:), intent(in)    :: Diag
    real, intent(in)    :: eta_x
    real, intent(in)    :: eta_y

    logical, optional, intent(in) :: CoulombFlag
    real, optional, intent(in)    :: rconv
    integer, optional, intent(out) :: niter_out
    real, optional, intent(out) :: error_out

    real, save :: X_min, Y_min, Z_min, XYZ_min, XYZ_max
    real :: alpha, beta!, gamma
    real, dimension(2) :: r
    real :: A0, b, denom, error, lhs, error_overall
    real :: Diag_min, Diag_max, absS_max

    integer, SAVE :: numberIterations_max=20       ! Maximum number of iterations
    real,    SAVE :: relativeError_max=1.e-04       ! Maximum relative error required
    integer :: i, j, ix, iy, iz

    logical, save :: flagini=.true.

    if (flagini) then
       call init
       Nx=size(U,dim=1)
       Ny=size(U,dim=2)
       Nz=size(U,dim=3)
       if (debugADI) write(*,*) 'In ADI, Nx, Ny, Nz:', Nx, Ny, Nz
       allocate(alpx(1:Nx-1))
       allocate(betx(1:Nx-1))
       allocate(alpy(1:Ny-1))
       allocate(bety(1:Ny-1))
       allocate(alpz(1:Nz-1))
       allocate(betz(1:Nz-1))
       allocate(U1(1:Nx,1:Ny,1:Nz),U2(1:Nx,1:Ny,1:Nz))
       X_min=4.*eta_x*sin(pi/(2.*(Nx-1)))**2
       Y_min=4.*eta_y*sin(pi/(2.*(Ny-1)))**2
       Z_min=4.*sin(pi/(2.*(Nz-1)))**2
       if (debugADI) write(*,*) 'In ADI, X_min, Y_min, Z_min:', X_min, Y_min, Z_min
       XYZ_min=min(X_min,Y_min,Z_min)
       XYZ_max=4.*max(eta_x,eta_y,1.)
       !       XYZ_min=0.33
       !       XYZ_max=5.7
       flagini=.false.
    end if

    if (present(CoulombFlag)) then
       numberIterations_max = 600
       relativeError_max=3.e-04
    else
       numberIterations_max = 20
       relativeError_max=1.e-04
    end if

    if (debugADI) write(*,'(A,i5,2x,e12.5)') 'In ADI_solve, numberIterations_max,relativeError_max: ',&
         & numberIterations_max,relativeError_max

    if (present(rconv)) then
       r=rconv
    else if (present(CoulombFlag)) then
       r=2.0
    else

       Diag_min=1.e+06
       Diag_max=-1.e+06
       do ix = 2,Nx-1
          do iy = 2,Ny-1
             do iz = 2,Nz-1
                if (Diag(ix,iy,iz).lt.Diag_min) Diag_min=Diag(ix,iy,iz)
                if (Diag(ix,iy,iz).gt.Diag_max) Diag_max=Diag(ix,iy,iz)
             end do
          end do
       end do

       if (Diag_min.lt.0.) write(*,*) ' Warning in ADI_solve: Diag_min < 0'

       ! Estimate lower (alpha) and upper (beta) boundaries of the eigenvalues
       ! of the matrices  tildeX, tildeY and tildeZ:

       alpha = XYZ_min + Diag_min
       beta = XYZ_max + Diag_max

       if (debugADI) write(*,*) 'In ADI_solve, alpha, beta:', alpha, beta

       if (alpha.le.0.) then
          write(*,*) 'In ADI_solve: alpha less than zero, alpha=', alpha
          stop
       end if

       r=sqrt(alpha*beta)

       ! Use the formulas for the two-step P-R method (they are, actually, derived
       ! for the two-dimensional problem, but we apply them for the three-dimensional problem,
       ! see also subroutine yuint in module yukawa):
       !    gamma=(alpha*beta)**0.25*sqrt((alpha+beta)/2.)
       !    r(1)=gamma-sqrt(gamma**2-alpha*beta)
       !    r(2)=gamma+sqrt(gamma**2-alpha*beta)

       !    gamma=(beta/alpha)**(0.5/float(numberIterations))

    end if

    absS_max=0.
    do ix = 2,Nx-1
       do iy = 2,Ny-1
          do iz = 2,Nz-1
             if (abs(S(ix,iy,iz)).gt.absS_max) absS_max=abs(S(ix,iy,iz))
          end do
       end do
    end do

    if (debugADI) write(*,*) ' In ADI_solve: absS_max = ', absS_max

    ! Copying in oder to have the correct boundary conditions also for U1 and U2:
    U1 = U
    U2 = U

    alpx(Nx-1)=0.
    alpy(Ny-1)=0.
    alpz(Nz-1)=0.

    j=0
    Loop_over_iterations : do

       j=j+1

       i=2-mod(j,2)

       !r(i)=alpha*gamma**(2.*float(j)-1.)

       if (debugADI) write(*,*) 'In ADI_solve, r:', r

       do iy = 2,Ny-1
          do iz = 2,Nz-1

             ! Inversion along x-axis for fixed y and z:

             betx(Nx-1) = U1(Nx,iy,iz)

             do ix = Nx-1,2,-1
                A0 = 2. + ( Diag(ix,iy,iz) + r(i) ) / eta_x

                b = (  ( r(i) - 2.*(eta_y+1.+Diag(ix,iy,iz)) ) * U(ix,iy,iz) &
                     &+ eta_y * ( U(ix,iy-1,iz) + U(ix,iy+1,iz) ) &
                     &+ U(ix,iy,iz-1) + U(ix,iy,iz+1) + S(ix,iy,iz) ) / eta_x

                denom = A0 - alpx(ix)
                if (denom.eq.0.) then
                   write(*,*) 'In ADI_solve, denom_x=0.'
                   stop
                end if

                alpx(ix-1) = 1. / denom
                betx(ix-1) = ( betx(ix) + b ) / denom
             end do

             do ix = 1,Nx-2
                U1(ix+1,iy,iz) = alpx(ix) *  U1(ix,iy,iz) + betx(ix)
             end do

          end do
       end do


       do ix = 2,Nx-1
          do iz = 2,Nz-1

             ! Inversion along y-axis for fixed x and z:
             bety(Ny-1) = U2(ix,Ny,iz)

             do iy = Ny-1,2,-1
                A0 = 2. + ( Diag(ix,iy,iz) + r(i) ) / eta_y

                b = (  ( r(i) - 2.*(eta_x+1.+Diag(ix,iy,iz)) ) * U1(ix,iy,iz) &
                     &+ eta_x * ( U1(ix-1,iy,iz) + U1(ix+1,iy,iz) ) &
                     &+ U1(ix,iy,iz-1) + U1(ix,iy,iz+1) + S(ix,iy,iz) ) / eta_y

                denom = A0 - alpy(iy)
                if (denom.eq.0.) then
                   write(*,*) 'In ADI_solve, denom_y=0.'
                   stop
                end if

                alpy(iy-1) = 1. / denom
                bety(iy-1) = ( bety(iy) + b ) / denom
             end do

             do iy = 1,Ny-2
                U2(ix,iy+1,iz) = alpy(iy) *  U2(ix,iy,iz) + bety(iy)
             end do

          end do
       end do


       do ix = 2,Nx-1
          do iy = 2,Ny-1

             ! Inversion along z-axis for fixed x and y:
             betz(Nz-1) = U(ix,iy,Nz)

             do iz = Nz-1,2,-1
                A0 = 2. + Diag(ix,iy,iz) + r(i)

                b = ( r(i) - 2.*(eta_x+eta_y+Diag(ix,iy,iz)) ) * U2(ix,iy,iz) &
                     &+ eta_x * ( U2(ix-1,iy,iz) + U2(ix+1,iy,iz) ) &
                     &+ eta_y * ( U2(ix,iy-1,iz) + U2(ix,iy+1,iz) ) + S(ix,iy,iz)

                denom = A0 - alpz(iz)
                if (denom.eq.0.) then
                   write(*,*) 'In ADI_solve, denom_z=0.'
                   stop
                end if

                alpz(iz-1) = 1. / denom
                betz(iz-1) = ( betz(iz) + b ) / denom
             end do

             do iz = 1,Nz-2
                U(ix,iy,iz+1) = alpz(iz) *  U(ix,iy,iz) + betz(iz)
             end do

          end do
       end do


       ! Check  the error:
       error_overall = 0.

       do ix = 2,Nx-1
          do iy = 2,Ny-1
             do iz = 2,Nz-1

                lhs =  -eta_x * ( U(ix-1,iy,iz) + U(ix+1,iy,iz) ) &
                     & -eta_y * ( U(ix,iy-1,iz) + U(ix,iy+1,iz) ) &
                     & -U(ix,iy,iz-1) - U(ix,iy,iz+1) &
                     & + ( 2.*(eta_x+eta_y+1.) + 3.*Diag(ix,iy,iz) ) * U(ix,iy,iz)

                error = abs(lhs - S(ix,iy,iz))
                if (error.gt.error_overall) error_overall=error

             end do
          end do
       end do

       if (absS_max.gt.1.e-10) error_overall=error_overall/absS_max

       if (debugADI) write(*,*) 'In ADI_solve, iteration, error:', j, error_overall

       if ( error_overall .lt. relativeError_max ) then
          exit
       else if ( j .eq. numberIterations_max ) then
          write(*,*) ' Warning in ADI_solve: exit after', numberIterations_max,' iterations'
          write(*,*) ' Accuracy reached: ', error_overall
          exit
       end if

    end do Loop_over_iterations

    if (present(niter_out)) niter_out=j
    if (present(error_out)) error_out=error_overall

  end subroutine ADI_solve

end module ADI
