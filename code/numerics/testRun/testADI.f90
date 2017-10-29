! Testing routine for the ADI routines

   program testADI

   use ADI
   use constants, only : pi

   implicit none

   real, parameter :: g_coulomb = 0.0014398  !=e^2*0.2 Gev fm=1/137*0.2 GeV fm
   integer, parameter :: nx=40, ny=40, nz=40
   real, parameter :: dx=0.5, dy=0.5, dz=0.5
   real, dimension(1:3) :: r, r0
   real, dimension(-nx:nx,-ny:ny,-nz:nz) :: source, field, diag, field_exact
   real :: eta_x, eta_y, d, diff_ini, diff_fin, rconv, error
   integer :: i, j, k, n, niter

   eta_x=(dz/dx)**2
   eta_y=(dz/dy)**2

   r0(1:3)=(/30.,0.,0./)    ! Position of the point-like charge

   open(1,file='testADI.dat',status='unknown')
   write(1,*)'# Grid parameters: '
   write(1,*)'# nx, ny, nz: ', nx, ny, nz
   write(1,*)'# dx, dy, dz: ', dx, dy, dz
   write(1,*)'# Charge position: ', r0 
   write(1,*)'# rconv:          diff_ini:     diff_fin:           Niter: error:'

   open(11,file='field1.dat',status='unknown')
   open(12,file='field2.dat',status='unknown')
   open(13,file='field3.dat',status='unknown')

   do n=1,100

      rconv=0.01*float(n)

      diag=0.
      source=0.
      field=0.
   
      diff_ini=0.
      do i=-nx,nx
        r(1)=float(i)*dx
        do j=-ny,ny
          r(2)=float(j)*dy
          do k=-nz,nz
            r(3)=float(k)*dz
            d=sqrt((r(1)-r0(1))**2+(r(2)-r0(2))**2+(r(3)-r0(3))**2)
            field_exact(i,j,k)=g_coulomb/d
            if(abs(i)==nx .or. abs(j)==ny .or. abs(k)==nz) field(i,j,k)=field_exact(i,j,k) ! Boundary conditions 
            diff_ini=diff_ini+abs(field(i,j,k)-field_exact(i,j,k))
          end do
        end do
      end do

      !call ADI_solve_Douglas(field,source,diag,eta_x,eta_y,.true.,rconv,niter, error)
      call ADI_Coulomb(field,source,eta_x,eta_y,rconv,niter, error)
      !call ADI_solve(field,source,diag,eta_x,eta_y,.true.,rconv,niter, error)

      diff_fin=0.
      do i=-nx,nx
        do j=-ny,ny
          do k=-nz,nz
            diff_fin=diff_fin+abs(field(i,j,k)-field_exact(i,j,k))
          end do
        end do
      end do

      write(1,'(3(2x,e13.6),2x,i4,2x,e13.6)') rconv, diff_ini, diff_fin, niter, error 

      r(1)=-10.
      i=nint(r(1)/dx)
      write(11,*)'# rconv: ', rconv
      write(11,*)'# x= ', r(1)
      write(11,*)'# y:    z:  (field-field_exact):'
      do j=-ny,ny
        do k=-nz,nz      
          write(11,'(3(2x,e13.6))') float(j)*dy, float(k)*dz, field(i,j,k)-field_exact(i,j,k)   
        end do
      end do

      r(1)=0.
      i=nint(r(1)/dx)
      write(12,*)'# rconv: ', rconv
      write(12,*)'# x= ', r(1)
      write(12,*)'# y:    z:  (field-field_exact):'
      do j=-ny,ny
        do k=-nz,nz      
          write(12,'(3(2x,e13.6))') float(j)*dy, float(k)*dz, field(i,j,k)-field_exact(i,j,k)   
        end do
      end do

      r(1)=10.
      i=nint(r(1)/dx)
      write(13,*)'# rconv: ', rconv
      write(13,*)'# x= ', r(1)
      write(13,*)'# y:    z:  (field-field_exact):'
      do j=-ny,ny
        do k=-nz,nz      
          write(13,'(3(2x,e13.6))') float(j)*dy, float(k)*dz, field(i,j,k)-field_exact(i,j,k)   
        end do
      end do

   end do

   end program testADI   
   


         
         



 
