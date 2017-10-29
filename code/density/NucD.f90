!******************************************************************************
!****m* /NucD
! NAME
! module NucD
! PURPOSE
! This module includes routines to calculate nuclear radii for
! neutrons and protons for all nuclei.
! AUTHOR
! Horst Lenske
!******************************************************************************

module NucD
  private

  public :: DFS, nuclfit

contains
  !****************************************************************************
  !****s* NucD/DFS
  ! NAME
  ! subroutine DFS(A,Z,r,a,rho,rms,rchrg)
  ! PURPOSE
  ! Proton/Neutron Density Parameters from Skyrme Systematics
  ! Fermi Shape are used with half-density radius and diffuseness:
  ! NOTES
  ! * R_{q}=R_{0q}*A^(1/3)+R_{1q}+R_{2q}*AS
  ! * a_{q}=a_{0q}        +a_{1q}*AS
  !
  ! A=Mass Number, AS=(N-Z)/A Asymmetry, \tau=1,2 for p,n
  !
  ! Convention:
  ! * Protons: q = 1,
  ! * Neutrons: q = 2,
  ! * weighted Average: q = 3, with  e.g. R_av=(N*R_n+Z*R_p)/A
  ! etc.
  !
  ! INPUTS
  ! * integer, intent(in) :: A      ! Mass number of nucleus
  ! * integer, intent(in) :: Z      ! Charge
  !
  ! OUTPUT
  ! * real, dimension(1:3),intent(out) :: r
  !   -- radii
  ! * real, dimension(1:3),intent(out) :: a
  !   -- surface parameter
  ! * real, dimension(1:3),intent(out) :: rho
  !   -- contains the central densities (from normalization to Z/N numbers)
  ! * real, dimension(1:3),intent(out),optional :: RMS
  !   -- contains the mass   rms-radii (obtained analytically)
  ! * real, dimension(1:3),intent(out),optional :: RCHRG
  !   -- contains the charge rms-radii (obtained analytically)
  !
  ! Convention for indices:
  ! * i=1: Protons
  ! * i=2: Neutrons
  ! * i=3: weighted Average with  e.g. R_av=(N*R_n+Z*R_p)/A
  !
  !****************************************************************************
  subroutine DFS(Amass_int,Z_int,rd_real,ad_real,rhod_real,rms_real,rchrg_real)
    implicit double precision (a-h,o-z)
    integer, intent(in) :: Amass_int
    integer, intent(in) :: Z_int
    real, intent(out),dimension(1:3) :: rd_real
    real, intent(out),dimension(1:3) :: ad_real
    real, intent(out),dimension(1:3) :: rhod_real
    real, intent(out),dimension(1:3),optional :: rms_real
    real, intent(out),dimension(1:3),optional :: rchrg_real
    dimension  rd(3),ad(3),rhod(3),RMS(3),RCHRG(3)
    dimension rnp(3,2),anp(2,2),xmsq(2)
    data rnp/1.2490,-0.5401,-0.9582, 1.2131,-0.4415, 0.8931/
    data anp/0.4899,-0.1236, 0.4686,0.0741/
    data xmsq/0.7429d0,-0.113d0/

    RCHRG=0.

    aMass=aMass_int
    Z=Z_int

    pi=atan(1.d0)*4
    fpi=4*pi
    A1=Amass
    an=a1-z
    MM=int(A1)
    NZ=int(Z)
    NN=MM-NZ
    rd(3)=0.d0
    ad(3)=0.d0
    anz=z
    dnz=(an-z)/a1
    a3=a1**(1./3.)
    !* Proton/Neutron Parameters
    do  i=1,2
       rd(i) =a3*rnp(1,i)+rnp(2,i)+rnp(3,i)*dnz
       ad(i) =anp(1,i)+anp(2,i)*dnz
       !c Vol.Integral of a Fermi distrib. (analytical)
       rhod(i)=0.75*anz/(rd(I)**3*pi  *(1.+(pi*ad(i)/rd(i))**2))
       !c <r^2>        of a Fermi distrib. (analytical)
       x=(pi*ad(i)/rd(i))**2
       rms(i)=0.2d0*rd(i)**5*(1.d0+x/(0.3d0)*(1.+0.7d0*x))*  fpi*rhod(i)/anz
       !c <r^2> NUMERICALLY:
       dx=0.2d0
       RX=5.d0*RD(i)
       NX=int(RX/dx)
       DMSS=0.d0
       DRMS=0.d0
       x=0.d0
       do n=1,NX
          x=x+dx
          xsq=x*x
          dd=rhod(i)/(1.d0+exp((x-rd(i))/ad(i)))
          DMSS=DMSS+dd*xsq
          DRMS=DRMS+dd*xsq*xsq
       end do
       rms(i)=DRMS/DMSS
       !c
       !c charge radius (xmsq are the <r^2> of p/n)
       !c
       if (i.eq.1) then
          rch=rms(i)+xmsq(i)
          RCHRG(i)=sqrt(rch)
       else
          rch=xmsq(i)/6.
          RCHRG(i)=-sqrt(abs(xmsq(i))/6.)
       end if
       rhod(i)=0.75*anz/(rd(i)**3*pi   *(1.+(pi*ad(i)/rd(i))**2))
       !c
       !c Prepare averaged R_0 and diffuseness:
       !c
       rd(3)=rd(3)+anz*rd(i)
       ad(3)=ad(3)+anz*ad(i)
       RCHRG(3)=RCHRG(3)+anz*rch
       anz=an
    end do
    !c
    !c Averaged R_0 and diffuseness (do averaging):
    !c
    RCHRG(3)=sqrt(RCHRG(3)/Z)
    rd(3)=rd(3)/a1
    ad(3)=ad(3)/a1
    rhod(3)=0.75*a1/(rd(3)**3*(1.+(pi*ad(3)/rd(3))**2)*pi)
    i=3
    x=(pi*ad(i)/rd(i))**2
    rms(i)=0.2d0*rd(i)**5*(1.d0+x/(0.3d0)*(1.+0.7d0*x))*   fpi*rhod(i)/a1
    rms(3)=(Z*rms(i)+AN*rms(2))/Amass
    do i=1,3
       rms(i)=sqrt(rms(i))
    end do
    !c
    write(6,6005)NZ,NN,MM,  &
         & (rd(I),i=1,3),(rd(i)/a3,i=1,3),(ad(I),i=1,3), &
         & (RMS(i),i=1,3),(RCHRG(i),i=1,3),(rhod(I),I=1,3)

6005 FORMAT(/28X,'Density Parameters from Systematics:'/  &
         & 28x,'   Protons','  Neutrons',3x,'Average'/ &
         & 10x,'   Number        :',3I10/               &
         & 10x,'   R             :',3f10.4,'  [fm]'/    &
         & 10x,'   R(reduced)    :',3f10.4,'  [fm]'/    &
         & 10x,'   a             :',3f10.4,'  [fm]'/    &
         & 10x,'   <r^2>         :',3f10.4,'  [fm]'/    &
         & 10x,'   <r^2>(charge) :',3f10.4,'  [fm]'/    &
         & 10x,'   rho(central)  :',3f10.4,'  [A/fm**3]')

    rd_real=rd
    ad_real=ad
    rhod_real=rhod
    if (present(rms_real)) rms_real=rms
    if (present(rchrg_real)) rchrg_real=rchrg

    return
  end subroutine DFS
  !c---------------------------------------------------------------------
!   double Precision Function Rqq(I3,A,Z)
!     implicit double precision (a-h,o-z)
!     !c
!     dimension rnp(3,2),anp(2,2),xmsq(2)
!     !c
!     data rnp/1.2490,-0.5401,-0.9582,  1.2131,-0.4415, 0.8931/
!     !c
!     data anp/0.4899,-0.1236, 0.4686,0.0741/
!     data xmsq/0.7429d0,-0.113d0/
!     !c
!     !c Mass and asymmetry factors:
!     !c
!     A3=A**(1.d0/3)
!     DNZ=(A-2*Z)/A
!     !c
!     !c Proton (I3=1) and neutron(I3=2) Mass Radius
!     !c
!     Rqq=A3*rnp(1,i3)+rnp(2,i3)+rnp(3,i3)*DNZ
!     !c
!     return
!   end Function Rqq
  !c---------------------------------------------------------------------
!!$  double Precision Function aqq(I3,A,Z)
!!$    implicit double precision (a-h,o-z)
!!$    !c
!!$    dimension rnp(3,2),anp(2,2),xmsq(2)
!!$    !c
!!$    data rnp/1.2490,-0.5401,-0.9582,  1.2131,-0.4415, 0.8931/
!!$    !c
!!$    data anp/0.4899,-0.1236, 0.4686,0.0741/
!!$    data xmsq/0.7429d0,-0.113d0/
!!$    !c
!!$    !c Asymmetry factor:
!!$    !c
!!$    DNZ=(A-2*Z)/A
!!$    !c
!!$    !c Proton (I3=1) and neutron(I3=2) diffusivity
!!$    !c
!!$    aqq=anp(1,i3)+anp(2,i3)*DNZ
!!$    !c
!!$    return
!!$  end Function aqq


  !****************************************************************************
  !****s* NucD/nuclfit
  ! NAME
  ! subroutine nuclfit(ma,mz,mradi,msurf,mdens)
  ! PURPOSE
  ! Evaluates basic properties of a nucleus
  ! INPUTS
  ! * integer ma :: mass of nucleus in units of nucleon masses
  ! * integer mz :: charge of nucleus in units of e
  ! RESULT
  ! * real mradi :: nuclear radius
  ! * real msurf :: nuclear surface parameter
  ! * real mdens :: nuclear maximal density
  ! AUTHOR
  ! Horst Lenske
  !****************************************************************************
  subroutine nuclfit(ma,mz,mradi,msurf,mdens)
    implicit double precision (a-h,o-z)
    real,intent(out)::    msurf, mradi, mdens
    integer ma, mz
    COMMON/RADPN/rd(3),ad(3),rhod(3),RMS(3)
    dimension rnp(3,2),anp(2,2)
    data rnp/1.2490,-0.5401,-0.9582,1.2131,-0.4415, 0.8931/
    data anp/0.4899,-0.1236, 0.4686,0.0741/
    data xmsq/0.64d0/

    !**    index 1 = protons
    !**    index 2 = neutron
    !**    index 3 = together

    !      write(*,*)'in nuclfit ', ma,mz
    if ((ma.le.0).or.(mz.le.0)) then
       write(*,*) 'Invalid input for nuclfit.'
       write(*,*) 'Mass=',ma,'Charge=',mz
       stop
    end if

    a1 = float(ma)
    z  = float(mz)
    pi=atan(1.d0)*4.d0
    fpi=4.d0*pi
    an=a1-z
    rd(3)=0.d0
    ad(3)=0.d0
    anz=z
    dnz=(an-z)/a1
    a3=a1**(1./3.)
    do i=1,2
       rd(i) =a3*rnp(1,i)+rnp(2,i)+rnp(3,i)*dnz
       ad(i) =anp(1,i)+anp(2,i)*dnz
       rhod(i)=0.75*anz/(rd(I)**3*pi*(1.+(pi*ad(i)/rd(i))**2))
       !c <r**2>
       x=(pi*ad(i)/rd(i))**2
       rms(i)=0.2d0*rd(i)**5*(1.d0+x/(0.3d0)*(1.+0.7d0*x))*fpi*rhod(i)/anz
       !c correct to point-particle distribution (xmsq is the rms of the nucleon)
       rmsp=rms(i)-xmsq
       pcor=1
       rd(i)=pcor*rd(i)
       rms(i)=rmsp
       rhod(i)=0.75*anz/(rd(i)**3*pi*(1.+(pi*ad(i)/rd(i))**2))
       !c
       rd(3)=rd(3)+anz*rd(i)
       ad(3)=ad(3)+anz*ad(i)
       anz=an
    end do
!!$    RC=sqrt(rms(1)/0.6d0)
    rd(3)=rd(3)/a1
    ad(3)=ad(3)/a1
    rhod(3)=0.75*a1/(rd(3)**3*(1.+(pi*ad(3)/rd(3))**2)*pi)
    i=3
    x=(pi*ad(i)/rd(i))**2
    rms(i)=0.2d0*rd(i)**5*(1.d0+x/(0.3d0)*(1.+0.7d0*x))*fpi*rhod(i)/a1
    do i=1,3
       rms(i)=sqrt(rms(i))
    end do

    mradi = sngl(rd(3))
    msurf = sngl(ad(3))
    mdens = sngl(rhod(3))

    return
  end subroutine nuclfit



end module NucD
