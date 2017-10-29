!******************************************************************************
!****m* /shadowing
! NAME
! module shadowing
!
! PURPOSE
! This module provides the calculation of shadowing for
! the vector meson components of a (high energetic) virtual photon:
! The couplings are rescaled...
!
! (cf. Th.Falter, PhD thesis)
!
!******************************************************************************
module shadowing

  use nucleusDefinition

  implicit none
  private

  public :: AeffCalc

contains


  !****************************************************************************
  !****s* shadowing/AeffCalc
  ! NAME
  ! subroutine AeffCalc(targetNuc,position,egamma,qsq2,shadfac)
  !
  ! PURPOSE
  ! calculate shadowing factor R_V for vector meson components:
  !   e^2/g_V^2 -> e^2/g_V^2* R_V
  !
  ! INPUTS
  ! * type(tNucleus)      :: targetNuc -- nucleus to consider
  ! * real,dimension(1:3) :: position  -- spatial coordinate
  ! * real                :: egamma -- energy of photon (=nu)
  ! * real                :: qsqs2  -- transfered momentum (=Q2)
  !
  ! OUTPUT
  ! * real,dimension(4)   :: shadfac -- the shadowing factors for
  !   rho0, omega, phi, J/psi
  ! NOTES
  ! * cf. PyVP[ScaleVMD]:  Pythia stores the inverse of the couplings,
  !   therefore f_V^2/4pi -> f_V^2/4pi / shadfac
  ! * cf. Th.Falter, PhD thesis, p.88ff, p117ff
  ! * the underlying assumptions seem not to hold for small A (i.e. He):
  !   assumption (5.13) is not very important;
  !   assumption (5.11) not tested!
  !****************************************************************************
  subroutine AeffCalc(targetNuc,position,egamma,qsq2,shadfac)
    use constants, only: hbarc

    integer,parameter :: nvec=4

    type(tNucleus),      pointer     :: targetNuc
    real,dimension(1:3), intent(in)  :: position
    real,                intent(in)  :: egamma,qsq2
    real,                intent(out) :: shadfac(nvec)

    real, parameter :: qc=0.78, dz2=1e-01, rhomin=1e-08
    integer :: i,k,j
    real :: z1,z2,r2s,rho2,kvsq,kgamma,inte,xbes,r,cost!,wcm2,z3,rho,imfv
    logical,dimension(nvec) :: lmeson
    real,   dimension(nvec) :: sigv,alphav,kv,qv
    real,   dimension(nvec), parameter :: mv = (/0.770,0.782,1.019,3.097/)
    complex :: css, S
    logical, parameter :: corflag=.TRUE. ! 2-body correlation correction?

    if (targetNuc%mass .le. 2) then
       write(*,*) 'Trying to calculate shadowing for A <= 2! stop.'
       stop
    end if

    r=sqrt(position(1)**2+position(2)**2+position(3)**2)
    cost = position(3)/(r+1e-20)

    kgamma=sqrt(qsq2+egamma**2)

    do i=1,nvec
       kvsq=egamma**2-mv(i)**2
       if (kvsq.gt.0.) then
          kv(i)=sqrt(kvsq)
          qv(i)=kgamma-kv(i)
          lmeson(i)=.TRUE.
       else
          kv(i)=1.
          qv(i)=1.
          lmeson(i)=.FALSE.
       end if
    end do

    do i=1,2
       sigv(i)  = 20.8*(1+0.766/sqrt(kv(i))) ! mb
       alphav(i)= -0.766/(sqrt(kv(i))+0.766)
    end do
    sigv(3)  =12.
    alphav(3)=0.
    sigv(4)  =2.2
    alphav(4)=0.


    z1=cost*r
    do k=1,nvec
       S=cmplx(0.,0.)
       if (lmeson(k)) then
          css = sigv(k)*cmplx(1,-alphav(k))/(2*10.) ! [fm^2]

          j=0
          do
             j=j+1
             z2=z1-(float(j)-0.5)*dz2

             inte = pathint(targetNuc,r,cost,z2)

             r2s=sqrt(r**2-z1**2+z2**2)
             rho2 = NucleusStaticDens(targetNuc,r2s,0)
             if (corflag) then
                xbes=qc*abs(z2-z1)/hbarc
                rho2=rho2*(1.-sin(xbes)/xbes)
             end if

             ! This is the large-A simplification:
!!$             S=S + css*rho2*exp(cmplx(0.,qv(k)*(z2-z1)/hbarc))&
!!$                  & *exp(-css*inte)

             ! This is the expression for all A:
             S=S + css*rho2*exp(cmplx(0.,qv(k)*(z2-z1)/hbarc))&
                  & *(1.-css*inte/targetNuc%mass)**targetNuc%mass


             if (z2.lt.0.and.rho2.lt.rhomin) exit
          end do
       end if
       shadfac(k)=(abs(1.-S*dz2))**2
    end do

  end subroutine AeffCalc


  !****************************************************************************
  !****f* shadowing/pathint
  ! NAME
  ! real function pathint(targetNuc,r,cost,zprime)
  !
  ! PURPOSE
  ! calculates \int_{z'}^z\rho(\vec b,z'')dz''
  ! with z=r*cos(\theta)
  !
  ! INPUTS
  ! * type(tNucleus)  :: targetNuc -- nucleus to consider
  ! * real :: r,cost,zprime -- as given
  !
  ! OUTPUT
  ! * function value -- as given
  !****************************************************************************
  real function pathint(targetNuc,r,cost,zprime)
    type(tNucleus), pointer :: targetNuc
    real,intent(in) :: r,cost,zprime

    integer,parameter :: nr=101,ncost=51
    real,   parameter :: dr=0.1,dcost=0.02

    integer :: ir,izp,icost

    real,save :: intes(0:nr,-nr:nr,-ncost:ncost)=0.
    logical,save :: initFlag=.true.

    if (initFlag) call inipath()

    ir=min(nint(r/dr),nr)
    izp=max(min(nint(zprime/dr),nr),-nr)
    icost=nint(cost/dcost)

    pathint=intes(ir,izp,icost)


  contains
    !**************************************************************************
    !****s* pathint/inipath
    ! NAME
    ! subroutine inipath
    !
    ! PURPOSE
    ! Initialise pathint
    !
    ! NOTES
    !**************************************************************************
    subroutine inipath()
      real,parameter :: dz=1e-01,rhomin=1e-08
      real :: dza,rsqr,rho,zprime2,z2,r2,cost2,inte2
      real,dimension(3) :: x
      integer :: i,nst,j,k,i2

      do i=0,nr
         r2=float(i)*dr
         do j=-nr,nr
            zprime2=float(j)*dr
            do k=-ncost,ncost
               cost2=max(min(float(k)*dcost,1.),-1.)
               z2=cost2*r2
               intes(i,j,k)=0
               if (zprime2.lt.z2) then
                  nst=int((z2-zprime2)/dz)+1
                  dza=(z2-zprime2)/float(nst)
                  x(1)=sqrt(max(1.-cost2**2,0.))*r2
                  x(2)=0.
                  inte2=0.
                  do i2=1,nst
                     x(3)=-(float(i2)-0.5)*dza+z2
                     rsqr=sqrt(x(1)**2+x(2)**2+x(3)**2)
                     rho = NucleusStaticDens(targetNuc,rsqr,0)
                     inte2=inte2+rho
                     if ((rho.lt.rhomin.and.x(3).lt.0)) exit
                  end do
                  intes(i,j,k)=inte2*dza
               end if
            end do
         end do
      end do

      initFlag=.false.
    end subroutine inipath
  end function pathint



end module shadowing
