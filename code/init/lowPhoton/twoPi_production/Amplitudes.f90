!******************************************************************************
!****m* /amplitudes_2Pi
! NAME
! module amplitudes_2Pi
! PURPOSE
! Calculates 2pi photoproduction amplitudes.
!******************************************************************************
module amplitudes_2Pi

  implicit none
  private

  type dos
     complex, dimension(2,2)::t1,t2
  end type dos

  !****************************************************************************
  !****g* amplitudes_2Pi/inMedium_delta_width
  ! SOURCE
  logical, save :: inMedium_delta_width=.false.
  ! PURPOSE
  ! turn of the in-medium-width of the delta
  !****************************************************************************

  !****************************************************************************
  !****g* amplitudes_2Pi/inMedium_delta_potential
  ! SOURCE
  logical, save :: inMedium_delta_potential=.false.
  ! PURPOSE
  ! turn of the in-medium-potential of the delta
  !****************************************************************************

  !****************************************************************************
  !****g* amplitudes_2Pi/inMedium_nucleon_potential
  ! SOURCE
  logical, save :: inMedium_nucleon_potential=.false.
  ! PURPOSE
  ! turn on the in-medium-potential of the nucleon
  !****************************************************************************

  !****************************************************************************
  !****g* amplitudes_2Pi/inMedium_pion_potential
  ! SOURCE
  logical, save :: inMedium_pion_potential=.false.
  ! PURPOSE
  ! turn on the in-medium-potential of the pion
  !****************************************************************************

  !****************************************************************************
  !****g* amplitudes_2Pi/buuPotential
  ! SOURCE
  logical, save :: buuPotential=.true.
  ! PURPOSE
  ! use buu potentials, else constants
  !****************************************************************************

  !****************************************************************************
  !****g* amplitudes_2Pi/debugFlag
  ! SOURCE
  logical, parameter :: debugFlag=.false.
  ! PURPOSE
  ! Switch for debug information
  !****************************************************************************

  logical, save :: initFlag=.true.

  public :: dos, ampli

contains

  !****************************************************************************
  !****s* amplitudes_2Pi/init
  ! NAME
  ! subroutine init
  ! PURPOSE
  ! Read out job card to initialize module parameters.
  !****************************************************************************
  subroutine init

    use output

    integer :: ios

    !**************************************************************************
    !****n* amplitudes_2Pi/Amplitudes_2Pi
    ! NAME
    ! NAMELIST /Amplitudes_2Pi/
    ! PURPOSE
    ! Includes the switches:
    ! * inMedium_delta_width
    ! * inMedium_delta_potential
    ! * inMedium_nucleon_potential
    ! * inMedium_pion_potential
    ! * buuPotential
    !**************************************************************************
    NAMELIST /Amplitudes_2Pi/ inMedium_delta_width, inMedium_delta_potential, inMedium_nucleon_potential, &
                              inMedium_pion_potential, buuPotential

    call Write_ReadingInput('Amplitudes_2Pi',0)

    rewind(5)
    read(5,nml=Amplitudes_2Pi,IOSTAT=IOS)
    call Write_ReadingInput('Amplitudes_2Pi',0,IOS)

    write(*,*) 'Set inMedium_delta_width       :' ,inMedium_delta_width
    write(*,*) 'Set inMedium_delta_potential   :' ,inMedium_delta_potential
    write(*,*) 'Set inMedium_nucleon_potential :' ,inMedium_nucleon_potential
    write(*,*) 'Set inMedium_pion_potential    :' ,inMedium_pion_potential
    write(*,*) 'BuuPotential    :' ,buuPotential
    call Write_ReadingInput('Amplitudes_2Pi',1)

  end subroutine init


  !****************************************************************************
  !****f* amplitudes_2Pi/ampli
  ! NAME
  ! function ampli(r,q,p1,p2,p4,p5,o,betaToLRF,mediumAtPosition) result(t)
  ! PURPOSE
  ! Calculates the 2pi photoproduction amplitude taken from Nacher et. al. NPA 695 (2001) 295.
  ! t is given in the CM frame, only transverse components are calculated.
  ! AUTHOR
  ! L. Alvarez-Ruso
  !****************************************************************************
  function ampli(r,q,p1,p2,p4,p5,o,betaToLRF,mediumAtPosition,position) result(t)
   use constants, only: pi, ii, hbarc, rhoNull
   use matrix_module, only: unim => unit2
   use minkowski, only: sigma1, sigma2, sigma3
   use mediumDefinition
   use paramamp  ! masses, coupling constants, etc
   use gauss_integration, only: sg20r, rg20r

   integer, intent(in) :: r  ! identifies the reaction channel:
                             !    r=1: gamma p -> pi+ pi- p
                             !    r=2: gamma p -> pi+ pi0 n
                             !    r=3: gamma p -> pi0 pi0 p
                             !    r=4: gamma n -> pi+ pi- n
                             !    r=5: gamma n -> pi- pi0 p
                             !    r=6: gamma n -> pi0 pi0 n
   real, dimension(0:3), intent(in) :: q     ! incoming photon 4-momentum
   real, dimension(0:3), intent(in) :: p1,p2 ! incoming and outgoing nucleon momentum
   real, dimension(0:3), intent(in) :: p4,p5 ! outgoing  pion momenta: p4(- or 0), p5(+ or 0)
   integer, intent(in) :: o ! refers to the order in which pions are emitted:
                            !      o=1 <-> p4,p5 (the pion with p5 is emitted first)
                            !      o=2 <-> p5,p4 (the pion with p4 is emitted first)
   real, dimension(1:3), intent(in) :: betaToLRF
   type(medium), intent(in) :: mediumAtPosition
   real, dimension(1:3), intent(in) :: position
   type(dos) :: t  ! amplitude

   !     other variables used
   type(dos) :: ta,tb,tc,td,te,tf,tg,th,ti,tj,tk,tl,tm,teo,tp,tq,tr,ts,tt,tu,tv,tw,tx
   complex, dimension(2,2) :: taux,ma1,ma2,ma3,ma4
   real, dimension(0:3) :: p4p,p5p,vaux1,vaux2,vaux3,qp!,vaux4
   real :: fo,f1,gm,ed,mud,g1,g2,g3,f1rop,f2rop,smandel
   real, parameter :: GEV_TO_INVfm = 1./hbarc

   ! Some propagators depend only at s, so it is convenient to calculate them
   ! only once for all phase space points at a fixed energy.
   real,save::q0prev=0.
   complex,save::gnf,gdf,gnsf,gnspf,gdsf
   INTEGER :: J

   if (initFlag) then
      call init
      initFlag=.false.
      if (debugFlag) then
         do j=0,100
            sMandel=((1.1+j*0.005)*GEV_TO_INVfm)**2
            write(400,*) sqrt(sMandel)/GEV_TO_INVfm, gamd(sMandel)/GEV_TO_INVfm
         end do
      end if
   end if

   if (q0prev.ne.q(0)) then
      q0prev=q(0)
      gnf=gn(p1+q)
      gdf=gd(p1+q)
      gnsf=gns(p1+q)
      gnspf=gnsp(p1+q)
      gdsf=gds(p1+q)
   end if
   !      gdsf=0.
   !      gnspf=0.

   !     in photons the ff(q2=0) for the N-N* transition are chosen
   !     depending on whether r takes place on p or n
   call photons(r,g1,g2,g3,f1rop,f2rop)

   p4p=p4-p4(0)/w(p2+p4)*(p2+p4)  ! pion momentum in the delta rest frame
   !      p4p=p4

   !      goto 20
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !     Diagram (a) [in fig 1]
   taux=c(r,'a',o)*e*(f/mpi)**2*gn(p2+p4)*ff(p5-q)* &
        &     (-p4(0)*sdotp(2.*p2+p4)/(2.*mn)+sdotp(p4))

   ta%t1=matmul(taux,sigma1)
   ta%t2=matmul(taux,sigma2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !     Diagram (b)
   taux=c(r,'b',o)*e*(f/mpi)**2*gn(p1-p5)*ff(p4-q)* &
        &     (-p5(0)*sdotp(2.*p1-p5)/(2.*mn)+sdotp(p5))

   tb%t1=matmul(sigma1,taux)
   tb%t2=matmul(sigma2,taux)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !     Diagram (c)

   taux=c(r,'c',o)*e*(f/mpi)**2*gn(p2+p4)*dpi(p5-q)*ff(p5-q)* &
        &     matmul(-p4(0)*sdotp(2.*p2+p4)/(2.*mn)+sdotp(p4), &
        &     -(p5(0)-q(0))*sdotp(p1+p2+p4)/(2.*mn)+sdotp(p5-q))

   tc%t1=taux*2.*p5(1)
   tc%t2=taux*2.*p5(2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !     Diagram (d)

   taux=c(r,'d',o)*e*(f/mpi)**2*gn(p1-p5)*dpi(p4-q)*ff(p4-q)* &
        &     matmul(-(p4(0)-q(0))*sdotp(p1+p2-p5)/(2.*mn)+sdotp(p4-q), &
        &     -p5(0)*sdotp(2.*p1-p5)/(2.*mn)+sdotp(p5))

   td%t1=taux*2*p4(1)
   td%t2=taux*2.*p4(2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !     Diagram (e)

   taux=c(r,'e',o)*e*(f/mpi)**2*gn(p2+p4)*gnf* &
        &       matmul(-p4(0)*sdotp(2.*p2+p4)/(2.*mn)+sdotp(p4), &
        &       -p5(0)*sdotp(p2+p4)/(2.*mn)+sdotp(p5))

   call photoff(r,'e',o,f1,gm)

   te%t1=matmul(taux,ii*gm*svp(q,1))/(2.*mn)
   te%t2=matmul(taux,ii*gm*svp(q,2))/(2.*mn)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !     Diagram (f)

   taux=c(r,'f',o)*e*(f/mpi)**2*gn(p2+p4)*gn(p1-p5)* &
        &       (-p4(0)*sdotp(2.*p2+p4)/(2.*mn)+sdotp(p4))

   call photoff(r,'f',o,f1,gm)

   tf%t1=matmul(matmul(taux,-2.*p5(1)*unim*f1+ii*gm*svp(q,1))/(2.*mn), &
        &       -p5(0)*sdotp(2.*p1-p5)/(2.*mn)+sdotp(p5))

   tf%t2=matmul(matmul(taux,-2.*p5(2)*unim*f1+ii*gm*svp(q,2))/(2.*mn), &
        &       -p5(0)*sdotp(2.*p1-p5)/(2.*mn)+sdotp(p5))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !     Diagram (g)

   taux=c(r,'g',o)*e*(f/mpi)**2*gn(p2-q)*gn(p1-p5)* &
        &       matmul(-p4(0)*sdotp(p1+p2-p5-q)/(2.*mn)+sdotp(p4),&
        &       -p5(0)*sdotp(2.*p1-p5)/(2.*mn)+sdotp(p5))

   call photoff(r,'g',o,f1,gm)

   tg%t1=matmul(2.*p2(1)*unim*f1+ii*gm*svp(q,1),taux)/(2.*mn)
   tg%t2=matmul(2.*p2(2)*unim*f1+ii*gm*svp(q,2),taux)/(2.*mn)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !     Diagram (h)

   taux=c(r,'h',o)*(f/mpi)*(fs/mpi)*(fgam/mpi)*gn(p2+p4)*gdf* &
        &       (-p4(0)*sdotp(2.*p2+p4)/(2.*mn)+sdotp(p4))
   !      *(p1(0)+q(0))/mdel

   vaux1=vp(p5,q)

   th%t1=matmul(taux,(-2.*ii*vaux1(1)*unim-sdotp(q)*p5(1)+sp(p5,q)*sigma1))
   th%t2=matmul(taux,(-2.*ii*vaux1(2)*unim-sdotp(q)*p5(2)+sp(p5,q)*sigma2))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !     Diagram (i)

   taux=c(r,'i',o)*e*(fs/mpi)**2*gd(p2+p4)*ff(p5-q)

   ti%t1=taux*(2.*p4p(1)*unim-ii*svp(p4p,1))
   ti%t2=taux*(2.*p4p(2)*unim-ii*svp(p4p,2))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !     Diagram (j)

   taux=c(r,'j',o)*e*(fs/mpi)**2*gd(p2+p4)*dpi(p5-q)*ff(p5-q)* &
        &       (2.*sp(p4p,p5-q)*unim-ii*sdotp(vp(p4p,p5-q)))

   tj%t1=taux*2.*p5(1)
   tj%t2=taux*2.*p5(2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !     Diagram (k)

   taux=c(r,'k',o)*(f/mpi)*(fs/mpi)*(fgam/mpi)*gn(p1-p5)*gd(p2+p4)* &
        &       (-p5(0)*sdotp(2.*p1-p5)/(2.*mn)+sdotp(p5))

   qp=q-q(0)/(p2(0)+p4(0))*(p2+p4)
   !      qp=q
   vaux1=vp(p4p,qp)

   !      tk%t1=matmul((-2.*ii*vaux1(1)*unim-sdotp(qp)*p4p(1)+ &
   !&       sp(p4p,qp)*sigma1)*(p2(0)+p4(0))/mdel,taux)

   !      tk%t2=matmul((-2.*ii*vaux1(2)*unim-sdotp(qp)*p4p(2)+ &
   !&       sp(p4p,qp)*sigma2)*(p2(0)+p4(0))/mdel,taux)

   tk%t1=matmul((-2.*ii*vaux1(1)*unim-sdotp(qp)*p4p(1)+ &
        &       sp(p4p,qp)*sigma1),taux)

   tk%t2=matmul((-2.*ii*vaux1(2)*unim-sdotp(qp)*p4p(2)+ &
        &       sp(p4p,qp)*sigma2),taux)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !     Diagram (l)

   taux=1./3.*c(r,'l',o)*(fs/mpi)*gd(p2+p4)*gnspf

   ma1=2.*sp(p4p,q)*unim-ii*sdotp(vp(p4p,q))
   ma2=2.*sp(p4p,p5)*unim-ii*sdotp(vp(p4p,p5))
   ma3=2.*sp(p5,q)*unim-ii*sdotp(vp(p5,q))
   ma4=matmul(ma2,ma3)

   fo=(g1*(1+q(0)/(2.*mn))+g2*(p1(0)+q(0)))*q(0)

   tl%t1=taux*(fnspdel*(matmul(ma1,-ii*g1/(2.*mn)*svp(q,1))- &
        &       (2.*p4p(1)*unim-ii*svp(p4p,1))*fo)+ &
        &       1./3.*gnspdel/(mpi**2)*(matmul(ma4,-ii*g1/(2.*mn)*svp(q,1))- &
        &       matmul(ma2,2.*p5(1)*unim-ii*svp(p5,1))*fo))

   tl%t2=taux*(fnspdel*(matmul(ma1,-ii*g1/(2.*mn)*svp(q,2))-&
        &       (2.*p4p(2)*unim-ii*svp(p4p,2))*fo)+&
        &       1./3.*gnspdel/(mpi**2)*(matmul(ma4,-ii*g1/(2.*mn)*svp(q,2))- &
        &       matmul(ma2,2.*p5(2)*unim-ii*svp(p5,2))*fo))

   !      tl%t1=0.157*taux*(fnspdel*(2.*p4p(1)*unim-ii*svp(p4p,1))+ &
   !&       1./3.*gnspdel/(mpi**2)*matmul(ma2,2.*p5(1)*unim-ii*svp(p5,1)))

   !      tl%t2=0.157*taux*(fnspdel*(2.*p4p(2)*unim-ii*svp(p4p,2))+ &
   !&       1./3.*gnspdel/(mpi**2)*matmul(ma2,2.*p5(2)*unim-ii*svp(p5,2)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !     Diagram (m)

   taux=c(r,'m',o)*e*(fs/mpi)**2*gnf*gd(p2+p4)*ff(p5-q)* &
        &       (2.*sp(p4p,p5)*unim-ii*sdotp(vp(p4p,p5)))

   call photoff(r,'m',o,f1,gm)

   tm%t1=matmul(taux,ii*gm*svp(q,1))/(2.*mn)
   tm%t2=matmul(taux,ii*gm*svp(q,2))/(2.*mn)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !     Diagram (o)

   taux=c(r,'o',o)*(fs/mpi)*(fdel/mpi)*(fgam/mpi)*gd(p2+p4)*gdf

   vaux1=vp(p5,q)
   vaux2=vp(p4p,q)

   teo%t1=taux*(ii*5./6.*sp(p4p,q)*p5(1)*unim-&
        &       ii*5./6.*sp(p5,q)*p4p(1)*unim-1./6.*sp(p4p,p5)*svp(q,1)- &
        &       1./6.*sdotp(p4p)*vaux1(1)+2./3.*sdotp(p5)*vaux2(1))

   teo%t2=taux*(ii*5./6.*sp(p4p,q)*p5(2)*unim- &
        &       ii*5./6.*sp(p5,q)*p4p(2)*unim-1./6.*sp(p4p,p5)*svp(q,2)- &
        &       1./6.*sdotp(p4p)*vaux1(2)+2./3.*sdotp(p5)*vaux2(2))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !     Diagram (p)

   if (s(p1-p5).lt.0.) then
      tp%t1=(0.,0.)
      tp%t2=(0.,0.)

   else

      taux=c(r,'p',o)*(fs/mpi)**2*gd(p2+p4)*gd(p1-p5)*ff(p4-q)

      p5p=p5

      vaux1=p1-2.*p5p
      vaux2=vp(p5p,q)
      vaux3=vp(p4p,q)

      call photodel(r,o,ed,mud)

      tp%t1=taux*(ed/2.*vaux1(1)/mdel*(2.*sp(p4p,p5p)*unim- &
           &         ii*sdotp(vp(p4p,p5p)))+ii*e*mud/mn*(ii*5./6.*sp(p5p,q)*p4p(1)*unim-  &
           &         ii*5./6.*sp(p4p,q)*p5p(1)*unim-1./6.*sdotp(p4p)*vaux2(1)-1./6.*sdotp(p5p)*vaux3(1)+  &
           &         2./3.*sp(p4p,p5p)*svp(q,1)))

      tp%t2=taux*(ed/2.*vaux1(2)/mdel*(2.*sp(p4p,p5p)*unim- &
           &         ii*sdotp(vp(p4p,p5p)))+ii*e*mud/mn*(ii*5./6.*sp(p5p,q)*p4p(2)*unim- &
           &         ii*5./6.*sp(p4p,q)*p5p(2)*unim-1./6.*sdotp(p4p)*vaux2(1)-1./6.*sdotp(p5p)*vaux3(1)+  &
           &         2./3.*sp(p4p,p5p)*svp(q,2)))

   end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !     Diagram (q)

   taux=c(r,'q',o)*(fs/mpi)*(gnsdel/mpi)*gd(p2+p4)*gnsf* &
        &       (2.*sp(p4p,p5)*unim-ii*sdotp(vp(p4p,p5)))

   !      tq%t1=matmul(taux,f2rop*(1.+q(0)/(2.*mn))*svp(q,1))
   !      tq%t2=matmul(taux,f2rop*(1.+q(0)/(2.*mn))*svp(q,2))

   !     J.A.
   if (r.le.3) then
      tq%t1=matmul(taux,(0.0173/mpi)*svp(q,1))
      tq%t2=matmul(taux,(0.0173/mpi)*svp(q,2))
   else
      tq%t1=matmul(taux,(-0.0112/mpi)*svp(q,1))
      tq%t2=matmul(taux,(-0.0112/mpi)*svp(q,2))
   end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !     Diagram (r)

   !      taux=ctil*c(r,'r',o)*gnsf*f2rop*(1.+q(0)/(2.*mn))

   !     J.A.
   if (r.le.3) then
      taux=ctil*c(r,'r',o)*gnsf*(0.0173/mpi)
   else
      taux=ctil*c(r,'r',o)*gnsf*(-0.0112/mpi)
   end if

   tr%t1=taux*svp(q,1)
   tr%t2=taux*svp(q,2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !     Diagram (s)

   taux=c(r,'s',o)*e*(ftil/mpi)**2*gns(p2+p4)*ff(p5-q)*sdotp(p4p)

   ts%t1=matmul(taux,sigma1)
   ts%t2=matmul(taux,sigma2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !     Diagram (t)

   if (s(p1-p5).lt.0.) then
      tt%t1=(0.,0.)
      tt%t2=(0.,0.)
   else
      p5p=p5-p5(0)/w(p1-p5)*(p1-p5)
      taux=c(r,'t',o)*e*(ftil/mpi)**2*gns(p1-p5)*ff(p4-q)*sdotp(p5p)

      tt%t1=matmul(sigma1,taux)
      tt%t2=matmul(sigma2,taux)
   end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !     Diagram (u)

   taux=1./3.*c(r,'u',o)*(fs/mpi)*gd(p2+p4)*gdsf

   ma1=2.*sp(p4p,q)*unim-ii*sdotp(vp(p4p,q))
   ma2=2.*sp(p4p,p5)*unim-ii*sdotp(vp(p4p,p5))
   ma3=2.*sp(p5,q)*unim-ii*sdotp(vp(p5,q))
   ma4=matmul(ma2,ma3)

   fo=(gp1*(1.+q(0)/(2.*mn))+gp2*(p1(0)+q(0)))*q(0)

   tu%t1=taux*(fdelsdel*(matmul(ma1,-ii*gp1/(2.*mn)*svp(q,1))- &
        &       (2.*p4p(1)*unim-ii*svp(p4p,1))*fo)+&
        &       1./3.*gdelsdel/(mpi**2)*(matmul(ma4,-ii*gp1/(2.*mn)*svp(q,1))-&
        &       matmul(ma2,2.*p5(1)*unim-ii*svp(p5,1))*fo))

   tu%t2=taux*(fdelsdel*(matmul(ma1,-ii*gp1/(2.*mn)*svp(q,2))- &
        &       (2.*p4p(2)*unim-ii*svp(p4p,2))*fo)+ &
        &       1./3.*gdelsdel/(mpi**2)*(matmul(ma4,-ii*gp1/(2.*mn)*svp(q,2))-&
        &       matmul(ma2,2.*p5(2)*unim-ii*svp(p5,2))*fo))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !     Diagram (v)

   taux=1./3.*c(r,'v',o)*grho*frho*drho(p4+p5)*ffrho(p4+p5)*gnspf

   fo=(g1*(1.+q(0)/(2.*mn))+g2*(p1(0)+q(0)))*q(0)
   vaux1=p5-p4

   tv%t1=-taux*(matmul(2.*sp(vaux1,q)*unim-ii*sdotp(vp(vaux1,q)),&
        &       -ii*g1/(2.*mn)*svp(q,1))-(2.*vaux1(1)*unim-ii*svp(vaux1,1))*fo)

   tv%t2=-taux*(matmul(2.*sp(vaux1,q)*unim-ii*sdotp(vp(vaux1,q)),&
        &       -ii*g1/(2.*mn)*svp(q,2))-(2.*vaux1(2)*unim-ii*svp(vaux1,2))*fo)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !     Diagram (w)

   taux=-c(r,'w',o)*e*sqrt(2.)*frho*frhonn/mrho*drho(p4+p5)*ffrho(p4+p5)
   vaux1=p5-p4

   tw%t1=taux*svp(vaux1,1)
   tw%t2=taux*svp(vaux1,2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !     Diagram (x)

   taux=1./3.*c(r,'x',o)*gprho*frho*drho(p4+p5)*ffrho(p4+p5)*gdsf

   fo=(gp1*(1.+q(0)/(2.*mn))+gp2*(p1(0)+q(0)))*q(0)
   vaux1=p5-p4

   tx%t1=-taux*(matmul(2.*sp(vaux1,q)*unim-ii*sdotp(vp(vaux1,q)),&
        &       -ii*gp1/(2.*mn)*svp(q,1))-(2.*vaux1(1)*unim-ii*svp(vaux1,1))*fo)

   tx%t2=-taux*(matmul(2.*sp(vaux1,q)*unim-ii*sdotp(vp(vaux1,q)),&
        &       -ii*gp1/(2.*mn)*svp(q,2))-(2.*vaux1(2)*unim-ii*svp(vaux1,2))*fo)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   t%t1=ta%t1+tb%t1+tc%t1+td%t1+te%t1+tf%t1+tg%t1+th%t1+ti%t1+tj%t1+tk%t1+  &
        &       tl%t1+tm%t1+teo%t1+tp%t1+tq%t1+tr%t1+ts%t1+tt%t1+tu%t1+ &
        &       tv%t1+tw%t1+tx%t1

   t%t2=ta%t2+tb%t2+tc%t2+td%t2+te%t2+tf%t2+tg%t2+th%t2+ti%t2+tj%t2+tk%t2+ &
        &       tl%t2+tm%t2+teo%t2+tp%t2+tq%t2+tr%t2+ts%t2+tt%t2+tu%t2+  &
        &       tv%t2+tw%t2+tx%t2

   !10     t%t1=th%t1+ti%t1+tj%t1+tk%t1+tm%t1+teo%t1+tp%t1
   !       t%t2=th%t2+ti%t2+tj%t2+tk%t2+tm%t2+teo%t2+tp%t2

   !10     t%t1=ta%t1+tb%t1+tc%t1+td%t1+te%t1+tf%t1+tg%t1
   !       t%t2=ta%t2+tb%t2+tc%t2+td%t2+te%t2+tf%t2+tg%t2

   !10      t%t1=tq%t1+tr%t1+ts%t1+tt%t1
   !        t%t2=tq%t2+tr%t2+ts%t2+tt%t2

   !10     t%t1=tk%t1-tr%t1
   !       t%t2=tk%t2-tr%t2

   !10      t=tk

   !10      t%t1=ta%t1+tb%t1+tc%t1+td%t1+te%t1+tf%t1+tg%t1
   !      +th%t1+ti%t1+tj%t1+ &
   !&       tk%t1+ tl%t1+ tm%t1+teo%t1+tp%t1+ tq%t1+ &
   !&       tr%t1+ &
   !      ts%t1+tt%t1+ &
   !&      tu%t1+ tv%t1+ tw%t1+ tx%t1

   !      t%t2=ta%t2+tb%t2+tc%t2+td%t2+te%t2+tf%t2+tg%t2
   !      +th%t2+ti%t2+tj%t2+ &
   !&       tk%t2+ tl%t2+ tm%t2+teo%t2+tp%t2+tq%t2+ &
   !&       tr%t2+ &
   !      ts%t2+tt%t2+ &
   !&      tu%t2+tv%t2+tw%t2+tx%t2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 contains

   real function s(p)
     real, dimension(0:3), intent(in):: p
     s=p(0)**2-p(1)**2-p(2)**2-p(3)**2
   end function s

   real function w(p)
     !     Calculates the inv. mass that corresponds to a given 4-momentum
     real, dimension(0:3), intent(in):: p
     w=sqrt(s(p))
   end function w

   real function sp(p1,p2)
     !     Scalar product of two 3-momenta
     real, dimension(0:3), intent(in):: p1,p2
     sp=p1(1)*p2(1)+p1(2)*p2(2)+p1(3)*p2(3)
   end function sp

   function vp(p1,p2)
     !     Vector product of two 3-momenta
     real, dimension(0:3), intent(in):: p1,p2
     real, dimension(0:3)::vp
     vp(0)=0.
     vp(1)=p1(2)*p2(3)-p1(3)*p2(2)
     vp(2)=p1(3)*p2(1)-p1(1)*p2(3)
     vp(3)=p1(1)*p2(2)-p1(2)*p2(1)
   end function vp

   function sdotp(p)
     !     Calculates (sigma.p)
     real, dimension(0:3), intent(in):: p
     complex, dimension(2,2)::sdotp
     sdotp=sigma1*p(1)+sigma2*p(2)+sigma3*p(3)
   end function sdotp

   function svp(p,i)
     !     i-component of the vector product [sigma, p]
     real, dimension(0:3), intent(in):: p
     integer, intent(in)::i
     complex, dimension(2,2)::svp
     select case (i)
     case (1)
        svp=sigma2*p(3)-sigma3*p(2)
     case (2)
        svp=sigma3*p(1)-sigma1*p(3)
     case (3)
        svp=sigma1*p(2)-sigma2*p(1)
     end select
   end function svp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !     Propagators

   complex function dpi(p)
     !     pion propagator
     use lorentzTrafo, only: lorentz
     use potentialModule, only: potential_LRF
     use IdTable, only: pion
     use particleDefinition

     real, dimension(0:3), intent(in)::p
     real, dimension(0:3) ::pboosted
     real :: selfEnergy
     type(particle) ::teilchen

     if (inMedium_pion_potential) then
        pBoosted=p
        ! Boost to LRF
        call Lorentz(betaToLRF,pBoosted,'Amplitudes')
        ! Include pion potential :
        call setToDefault(teilchen)
        teilchen%ID=pion
        teilchen%position=position
        teilchen%momentum=pBoosted/GEV_to_INVfm
        teilchen%antiparticle=.false.
        teilchen%perturbative=.true.
        selfEnergy=2.*sqrt(mpi**2+Dot_Product(pBoosted(1:3),pBoosted(1:3)))*potential_LRF(teilchen)*GEV_TO_INVfm ! Converting from GEV to fm^-1
        ! Evaluate scalar potential :
        dpi=1./(s(p)-mpi**2-selfenergy)
     else
        dpi=1./(s(p)-mpi**2)
     end if
   end function dpi

   complex function drho(p)
     !     rho propagator
     real, dimension(0:3), intent(in)::p
     drho=1./(s(p)-mrho**2+ii*mrho*gamrho(s(p)))
     !        drho=1./(s(p)-mrho**2+ii*sqrt(s(p))*gamrho(s(p)))
   end function drho


   complex function gn(p)
     !     nucleon propagator
     use lorentzTrafo, only: lorentz
     use potentialModule, only: potential_LRF
     use IdTable, only: nucleon
     use particleDefinition

     real, dimension(0:3), intent(in)::p
     real, dimension(0:3) ::pBoosted
     real::ene,selfenergy,sqs
     type(particle) ::teilchen

     if (inMedium_nucleon_potential) then
        ! Boost to LRF
        pBoosted=p
        sqs=sqrt(p(0)**2-p(1)**2-p(2)**2-p(3)**2)
        call Lorentz(betaToLRF,pBoosted,'Amplitudes')
        ! Include delta potential :
        call setToDefault(teilchen)
        teilchen%ID=nucleon
        teilchen%position=position
        teilchen%momentum=pBoosted/GEV_to_INVfm
        teilchen%antiparticle=.false.
        teilchen%perturbative=.true.
        if (buuPotential) then
           selfenergy=sqrt(mn**2+dot_product(pBoosted(1:3),pBoosted(1:3)))/mn*potential_LRF(teilchen)*GEV_TO_INVfm ! Converting from GEV to fm^-1
        else
           selfenergy=-0.05*Gev_To_INVFM*(mediumAtPosition%DensityNeutron+mediumAtPosition%DensityProton)/rhoNull
           if (debugFlag) write(*,*) 'self=', selfenergy
        end if
        gn=1./(sqs-(mn+selfenergy))
        if (debugFlag) then
           write(*,*) 'nucleon selfenergy',selfenergy,potential_LRF(teilchen)
        end if
     else
        ene=sqrt(mn**2+p(1)**2+p(2)**2+p(3)**2)
        gn=1./(p(0)-ene)
        !        *mn/ene
     end if
   end function gn

   complex function gd(p)
     !     delta propagator
     use lorentzTrafo, only: lorentz
     use baryonWidthMedium, only:WidthBaryonMedium
     use potentialModule, only: potential_LRF
     use IdTable, only: delta
     use particleDefinition

     real, dimension(0:3), intent(in)::p
     real, dimension(0:3)::pBoosted
     real::s,sqs,sqsGEV!,ene
     real :: gammaMed
     type(particle)  :: teilchen
     real :: selfenergy

     if (inMedium_delta_width) then
        ! In Medium delta propagator : Include delta Width
        ! Boost to LRF

        pBoosted=p
        call Lorentz(betaToLRF,pBoosted,'Amplitudes')
        s=p(0)**2-p(1)**2-p(2)**2-p(3)**2
        sqs=sqrt(s)
        if (inMedium_delta_potential) then
           ! Include also delta potential :
           call setToDefault(teilchen)
           teilchen%ID=delta
           teilchen%position=position
           teilchen%momentum=pBoosted/GEV_To_INVfm
           teilchen%antiparticle=.false.
           teilchen%perturbative=.true.
           if (buuPotential) then
              selfenergy=sqrt(mdel**2+dot_product(pBoosted(1:3),pBoosted(1:3)))/mdel*potential_LRF(teilchen)*GEV_TO_INVfm ! Converting from GEV to fm^-1
           else
              selfenergy=-0.03*Gev_To_INVFM*(mediumAtPosition%DensityNeutron+mediumAtPosition%DensityProton)/rhoNull
              if (debugFlag) write(*,*) 'self=', selfenergy,selfEnergy/GeV_TO_InvFm
           end if
        else
           selfenergy=0.
        end if
        pBoosted=pBoosted/GEV_to_INVfm
        sqsGEV=sqs/GEV_to_INVfm
        gammaMed=WidthBaryonMedium(delta,sqsGEV,pBoosted,mediumATposition)*GEV_TO_INVfm ! Converting from GEV to fm^-1
        gd=1./(sqs-(mdel+selfEnergy)+ii/2.*gammaMed)
        if (debugFlag) then
           write(*,*) 'delta selfenergy',selfenergy
           write(*,*) 'delta width',gammamed
        end if

     else if (.not.inMedium_delta_width.and.inMedium_delta_potential) then
        ! Boost to LRF
        pBoosted=p
        s=p(0)**2-p(1)**2-p(2)**2-p(3)**2
        sqs=sqrt(s)
        call Lorentz(betaToLRF,pBoosted,'Amplitudes')
        ! Include delta potential :
        call setToDefault(teilchen)
        teilchen%ID=delta
        teilchen%position=position
        teilchen%momentum=pBoosted/GEV_to_INVfm
        teilchen%antiparticle=.false.
        teilchen%perturbative=.true.
        if (buuPotential) then
           selfenergy=sqrt(mdel**2+dot_product(pBoosted(1:3),pBoosted(1:3)))/mdel*potential_LRF(teilchen)*GEV_TO_INVfm ! Converting from GEV to fm^-1
        else
           selfenergy=-0.03*Gev_To_INVFM*(mediumAtPosition%DensityNeutron+mediumAtPosition%DensityProton)/rhoNull
        end if

        if (debugFlag) then
           write(*,*) 'delta selfenergy',selfenergy/GeV_TO_InvFm,selfEnergy
        end if

        gd=1./(sqs-(mdel+selfenergy)+ii/2.*gamd(s))
     else
        s=p(0)**2-p(1)**2-p(2)**2-p(3)**2
        sqs=sqrt(s)
        !       Original vacuum delta propagator :
        gd=1./(sqs-mdel+ii/2.*gamd(s))
     end if
     !        *mdel/ene

   end function gd

   complex function gns(p)
     !     N*(1440) propagator

     real, dimension(0:3), intent(in)::p
     real::s,sqs!,ene

     s=p(0)**2-p(1)**2-p(2)**2-p(3)**2
     if (s.lt.0.) then
        gns=(0.,0.)
     else
        sqs=sqrt(s)
        !          ene=sqrt(mns**2+p(1)**2+p(2)**2+p(3)**2)
        gns=1./(sqs-mns+ii/2.*gamns(s))
        !          *mns/ene
     end if
   end function gns

   complex function gnsp(p)
     !     N*(1520) propagator

     real, dimension(0:3), intent(in)::p
     real::s,sqs!,ene

     s=p(0)**2-p(1)**2-p(2)**2-p(3)**2
     sqs=sqrt(s)

     !        ene=sqrt(mnsp**2+p(1)**2+p(2)**2+p(3)**2)

     gnsp=1./(sqs-mnsp+ii/2.*gamnsp(s))
     !        *mnsp/ene

   end function gnsp

   complex function gds(p)
     !     Delta(1700) propagator

     real, dimension(0:3), intent(in)::p
     real::s,sqs

     s=p(0)**2-p(1)**2-p(2)**2-p(3)**2
     sqs=sqrt(s)
     gds=1./(sqs-mdels+ii/2.*gamds(s))

   end function gds

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !     Decay widths

   real function gamrho(s)
     !     rho -> pi pi

     real,intent(in)::s
     real::pcm

     if (s.le.(4.*mpi**2)) then
        gamrho=0.
     else
        pcm=sqrt(s-4.*mpi**2)/2.
        gamrho=1./(6.*pi)*frho**2*pcm**3/s
     end if

   end function gamrho

   real function gamd(s)
     !     Delta -> N pi

     real,intent(in)::s

     if (s.le.(mn+mpi)**2) then
        gamd=0.
     else
        gamd=1./(6.*pi)*(fs/mpi)**2*mn*qcm(s)**3/sqrt(s)
     end if

   end function gamd

   real function gamns(s)
     !     Full width of the N*(1440) from three contributions:
     !     N* -> N pi, N*-> N Delta, N*-> N (pi pi)s-wave

     real,intent(in)::s
     real::gam1,gam2,gam3,pmin,pmax,p,om,sd,resu,xmin,xmax,x
     real,dimension(:),allocatable::absi,orde
     integer::n,ns,i

     !     N* -> N pi
     if (s.le.(mn+mpi)**2) then
        gam1=0.
     else
        gam1=3./(2.*pi)*(ftil/mpi)**2*mn*qcm(s)**3/sqrt(s)
     end if

     !     N*-> N Delta
     if (s.le.(mn+2.*mpi)**2) then
        gam2=0.
     else
        pmin=0.
        pmax=sqrt((s-mn**2-2.*mn*mpi)**2/4./s-mpi**2)
        n=1  ! precision of the gauss integrals
        allocate (absi(20*n))
        allocate (orde(20*n))
        call sg20r(pmin,pmax,n,absi,ns)
        do i=1,ns
           p=absi(i)
           om=sqrt(p**2+mpi**2)
           sd=s-2.*sqrt(s)*om+mpi**2
           orde(i)=p**4/om/((sqrt(sd)-mdel)**2+gamd(sd)**2/4.)*gamd(sd)
        end do
        deallocate (absi)
        call rg20r(pmin,pmax,n,orde,resu)
        deallocate (orde)
        gam2=1./3./pi**2*(gnsdel/mpi)**2*resu
     end if

     !       N*-> N (pi pi)s-wave
     if (s.le.(mn+2.*mpi)**2) then
        gam3=0.
     else
        xmin=4.*mpi**2
        xmax=(sqrt(s)-mn)**2
        n=2  ! precision of the gauss integrals
        allocate (absi(20*n))
        allocate (orde(20*n))
        call sg20r(xmin,xmax,n,absi,ns)
        do i=1,ns
           x=absi(i)
           orde(i)=sqrt(lam(s,x,mn**2)*lam(x,mpi**2,mpi**2))/x
        end do
        deallocate (absi)
        call rg20r(xmin,xmax,n,orde,resu)
        deallocate (orde)
        gam3=3.*ctil**2/2.**5/pi**3*mn/s*resu
     end if
     gamns=gam1+gam2+gam3
   end function gamns

   real function gamnsp(s)
     !     Full width of the N*(1520) from the main three contributions
     !     N pi (d wave), Delta pi, N rho

     real,intent(in)::s
     real::gam1,gam2,gam3,k,as,ad,mimin,mimax,mi,&
          &         srho,resu,ffrho,e4,e4min,e4max,e5,e5min,e5max,det,drho2,alph!,enmin,enmax,en,p,erho
     real,dimension(:),allocatable::absi1,orde1,absi2,orde2,absi,orde
     integer::n,ns,i,n1,ns1,n2,ns2,j

     !     N*(1520) -> N pi
     if (s.le.(mn+mpi)**2) then
        gam1=0.
     else
        gam1=gnsp0*(qcm(s)/qcm(mnsp**2))**5  ! gnsp0 <-> 55% of the total width
     end if

     !     N*(1520) -> Delta pi
     if (s.le.(mn+2.*mpi)**2) then
        gam2=0.
     else
        mimin=mn+mpi
        mimax=sqrt(s)-mpi
        n=1  ! precision of the gauss integrals
        allocate (absi(20*n))
        allocate (orde(20*n))
        call sg20r(mimin,mimax,n,absi,ns)
        do i=1,ns
           mi=absi(i)
           k=sqrt(lam(s,mi**2,mpi**2))/2./sqrt(s)
           ad=sqrt(4.*pi)/3.*gnspdel*(k/mpi)**2
           as=-sqrt(4.*pi)*fnspdel-ad
           orde(i)=gamd(mi**2)/((mi-mdel)**2+gamd(mi**2)**2/4.)* &
                &             mi*k*(as**2+ad**2)
        end do
        deallocate (absi)
        call rg20r(mimin,mimax,n,orde,resu)
        deallocate (orde)
        gam2=resu/sqrt(s)/(2*pi)**3
     end if

     !     N*(1520) -> N rho
     if (s.le.(mn+2.*mpi)**2) then
        gam3=0.
     else
        n1=1  ! precision of the gauss integrals
        n2=1  ! precision of the gauss integrals
        allocate (absi1(20*n1))
        allocate (orde1(20*n1))
        allocate (absi2(20*n2))
        allocate (orde2(20*n2))

        e4min=mpi
        e4max=(s-mn**2-2.*mn*mpi)/2./sqrt(s)

        call sg20r(e4min,e4max,n1,absi1,ns1)
        do i=1,ns1
           e4=absi1(i)
           alph=s+mpi**2-2.*sqrt(s)*e4
           det=sqrt((e4**2-mpi**2)*(alph-(mn+mpi)**2)*(alph-(mn-mpi)**2))
           e5min=((sqrt(s)-e4)*(alph-mn**2+mpi**2)-det)/2./alph
           e5max=((sqrt(s)-e4)*(alph-mn**2+mpi**2)+det)/2./alph

           call sg20r(e5min,e5max,n2,absi2,ns2)
           do j=1,ns2
              e5=absi2(j)
              srho=2.*sqrt(s)*(e4+e5)+mn**2-s
              if (srho.lt.(4.*mpi**2)) then
              end if

              !              ffrho=(lamrho2-mrho**2)/(lamrho2-srho)
              ffrho=1.
              drho2=1./((srho-mrho**2)**2+(mrho*gamrho(srho))**2)
              !              drho2=1./((srho-mrho**2)**2+(sqrt(srho)*gamrho(srho))**2)
              orde2(j)=drho2*ffrho**2* &
                   &               ((e4-e5)**2+2.*sqrt(s)*(e4+e5)-s-4.*mpi**2+mn**2)
           end do
           call rg20r(e5min,e5max,n2,orde2,orde1(i))
        end do
        call rg20r(e4min,e4max,n1,orde1,resu)
        gam3=mn/2./(2.*pi)**3*grho**2*frho**2*resu*(mnsp/sqrt(s))

        deallocate (absi1)
        deallocate (absi2)
        deallocate (orde1)
        deallocate (orde2)
     end if
     gamnsp=gam1+gam2+gam3
   end function gamnsp

   real function gamds(s)
     !     Full width of the Delta*(1700) from the main three contributions
     !     N pi (d wave), Delta pi (s and d waves), N rho (s-wave)

     real,intent(in)::s
     real::gam1,gam2,gam3,k,as,ad,mimin,mimax,mi,srho,resu, &
          &         e4,e4min,e4max,e5,e5min,e5max,det,drho2,alph,ffrho!,enmin,enmax,en,p,erho
     real,dimension(:),allocatable::absi,orde,absi1,orde1,absi2,orde2
     integer::n,ns,i,n1,ns1,n2,ns2,j

     !     Delta(1700) -> N pi
     if (s.le.(mn+mpi)**2) then
        gam1=0.
     else
        gam1=gds0*(qcm(s)/qcm(mdels**2))**5  ! gds0 <-> 15% of the total width
     end if

     !     Delta(1700) -> Delta pi
     !     Like the N*(1520)->Delta pi decay but with
     !     f(g)nspdel repalced by f(g)deldel and a factor of 15/8
     if (s.le.(mn+2.*mpi)**2) then
        gam2=0.
     else
        mimin=mn+mpi
        mimax=sqrt(s)-mpi
        n=1  ! precision of the gauss integrals
        allocate (absi(20*n))
        allocate (orde(20*n))
        call sg20r(mimin,mimax,n,absi,ns)
        do i=1,ns
           mi=absi(i)
           k=sqrt(lam(s,mi**2,mpi**2))/2./sqrt(s)
           ad=sqrt(4.*pi)/3.*gdelsdel*(k/mpi)**2
           as=-sqrt(4.*pi)*fdelsdel-ad
           orde(i)=gamd(mi**2)/((mi-mdel)**2+gamd(mi**2)**2/4.)* &
                &           mi*k*(as**2+ad**2)
        end do
        deallocate (absi)
        call rg20r(mimin,mimax,n,orde,resu)
        deallocate (orde)
        gam2=resu/sqrt(s)/(2*pi)**3*(15./8.)
     end if

     !     Delta(1700) -> N rho
     !     Like the N*(1520) but grho replaced by gprho and a factor 1/3
     if (s.le.(mn+2.*mpi)**2) then
        gam3=0.
     else
        n1=1  ! precision of the gauss integrals
        n2=1  ! precision of the gauss integrals
        allocate (absi1(20*n1))
        allocate (orde1(20*n1))
        allocate (absi2(20*n2))
        allocate (orde2(20*n2))

        e4min=mpi
        e4max=(s-mn**2-2.*mn*mpi)/2./sqrt(s)

        call sg20r(e4min,e4max,n1,absi1,ns1)
        do i=1,ns1
           e4=absi1(i)
           alph=s+mpi**2-2.*sqrt(s)*e4
           det=sqrt((e4**2-mpi**2)*(alph-(mn+mpi)**2)*(alph-(mn-mpi)**2))
           e5min=((sqrt(s)-e4)*(alph-mn**2+mpi**2)-det)/2./alph
           e5max=((sqrt(s)-e4)*(alph-mn**2+mpi**2)+det)/2./alph

           call sg20r(e5min,e5max,n2,absi2,ns2)
           do j=1,ns2
              e5=absi2(j)
              srho=2.*sqrt(s)*(e4+e5)+mn**2-s
              if (srho.lt.(4.*mpi**2)) then
                 write(*,*) sqrt(srho)*hb
              end if

              !            ffrho=(lamrho2-mrho**2)/(lamrho2-srho)
              ffrho=1.
              drho2=1./((srho-mrho**2)**2+(mrho*gamrho(srho))**2)
              !            drho2=1./((srho-mrho**2)**2+(sqrt(srho)*gamrho(srho))**2)
              orde2(j)=drho2*ffrho**2* &
                   &               ((e4-e5)**2+2.*sqrt(s)*(e4+e5)-s-4.*mpi**2+mn**2)
           end do
           call rg20r(e5min,e5max,n2,orde2,orde1(i))
        end do
        call rg20r(e4min,e4max,n1,orde1,resu)
        gam3=mn/2./(2.*pi)**3*gprho**2*frho**2*resu/3.*(mdels/sqrt(s))

        deallocate (absi1)
        deallocate (absi2)
        deallocate (orde1)
        deallocate (orde2)
     end if
     gamds=gam1+gam2+gam3
   end function gamds

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !     Form factors

   real function ff(p)
     real,dimension(0:3),intent(in)::p
     real::p2
     p2=p(0)**2-p(1)**2-p(2)**2-p(3)**2
     ff=(lam2-mpi**2)/(lam2-p2)
   end function ff

   real function ffrho(p)
     !     Form factor for the RNrho vertex
     real,dimension(0:3),intent(in)::p
     ffrho=(lamrho2-mrho**2)/(lamrho2-s(p))
   end function ffrho


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real function qcm(s)
     !     Returns the 3-momentum of the pion formed after the decay of a
     !     resonance (R->N pi) of inv. mass s in the rest frame of the resonance
     real,intent(in)::s
     qcm=sqrt((s-mpi**2-mn**2)**2-4.*mpi**2*mn**2)/2./sqrt(s)
   end function qcm

   real function lam(x,y,z)
     real,intent(in)::x,y,z
     lam=x**2+y**2+z**2-2.*(x*y+x*z+y*z)
   end function lam

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   complex function c(r,diag,o)
     !     Determines the numerical coeff. for each reaction
     !     according to Table 3 of Nacher et al.
     integer,intent(in)::r  ! reaction
     integer,intent(in)::o  ! order in which pions are emited
     character (len=1),intent(in):: diag  !diagram level in appendix 4

     select case (r)
     case (1)  ! gamma p -> p pi+ pi-
        select case (diag)
        case ('a')
           if (o.eq.1) then
              c=2.*ii
           else if (o.eq.2) then
              c=0.
           end if
        case ('b')
           if (o.eq.1) then
              c=-2.*ii
           else if (o.eq.2) then
              c=0.
           end if
        case ('c')
           if (o.eq.1) then
              c=2.*ii
           else if (o.eq.2) then
              c=0.
           end if
        case ('d')
           if (o.eq.1) then
              c=-2.*ii
           else if (o.eq.2) then
              c=0.
           end if
        case ('e')
           if (o.eq.1) then
              c=2.*ii
           else if (o.eq.2) then
              c=0.
           end if
        case ('f')
           if (o.eq.1) then
              c=2.*ii
           else if (o.eq.2) then
              c=0.
           end if
        case ('g')
           if (o.eq.1) then
              c=2.*ii
           else if (o.eq.2) then
              c=0.
           end if
        case ('h')
           if (o.eq.1) then
              c=2.*ii/9.
           else if (o.eq.2) then
              c=0.
           end if
        case ('i')
           if (o.eq.1) then
              c=ii/9.
           else if (o.eq.2) then
              c=-ii/3.
           end if
        case ('j')
           if (o.eq.1) then
              c=ii/9.
           else if (o.eq.2) then
              c=-ii/3.
           end if
        case ('k')
           if (o.eq.1) then
              c=-2.*ii/9.
           else if (o.eq.2) then
              c=0.
           end if
        case ('l')
           if (o.eq.1) then
              c=ii/3.
           else if (o.eq.2) then
              c=ii
           end if
        case ('m')
           if (o.eq.1) then
              c=ii/9.
           else if (o.eq.2) then
              c=ii/3.
           end if
        case ('o')
           if (o.eq.1) then
              c=-2./3.
           else if (o.eq.2) then
              c=1.
           end if
        case ('p')
           if (o.eq.1) then
              c=0.
           else if (o.eq.2) then
              c=ii/3.
           end if
        case ('q')
           if (o.eq.1) then
              c=-1./9.
           else if (o.eq.2) then
              c=-1./3.
           end if
        case ('r')
           if (o.eq.1) then
              c=2.
           else if (o.eq.2) then
              c=0.
           end if
        case ('s')
           if (o.eq.1) then
              c=2.*ii
           else if (o.eq.2) then
              c=0.
           end if
        case ('t')
           if (o.eq.1) then
              c=-2.*ii
           else if (o.eq.2) then
              c=0.
           end if
        case ('u')
           if (o.eq.1) then
              c=ii*sqrt(2./3.)
           else if (o.eq.2) then
              c=-ii*sqrt(3./2.)
           end if
        case ('v')
           if (o.eq.1) then
              c=-ii
           else if (o.eq.2) then
              c=0.
           end if
        case ('w')
           if (o.eq.1) then
              c=0.
           else if (o.eq.2) then
              c=0.
           end if
        case ('x')
           if (o.eq.1) then
              c=-ii*sqrt(2./3.)
           else if (o.eq.2) then
              c=0.
           end if
        end select

     case (2)  ! gamma p -> n pi+ pi0
        select case (diag)
        case ('a')
           if (o.eq.1) then
              c=-sqrt(2.)*ii
           else if (o.eq.2) then
              c=0.
           end if
        case ('b')
           if (o.eq.1) then
              c=0.
           else if (o.eq.2) then
              c=sqrt(2.)*ii
           end if
        case ('c')
           if (o.eq.1) then
              c=-sqrt(2.)*ii
           else if (o.eq.2) then
              c=0.
           end if
        case ('d')
           if (o.eq.1) then
              c=0.
           else if (o.eq.2) then
              c=sqrt(2.)*ii
           end if
        case ('e')
           if (o.eq.1) then
              c=-sqrt(2.)*ii
           else if (o.eq.2) then
              c=sqrt(2.)*ii
           end if
        case ('f')
           if (o.eq.1) then
              c=-sqrt(2.)*ii
           else if (o.eq.2) then
              c=sqrt(2.)*ii
           end if
        case ('g')
           if (o.eq.1) then
              c=-sqrt(2.)*ii
           else if (o.eq.2) then
              c=sqrt(2.)*ii
           end if
        case ('h')
           if (o.eq.1) then
              c=-ii*sqrt(2.)/9.
           else if (o.eq.2) then
              c=-2.*ii*sqrt(2.)/9.
           end if
        case ('i')
           if (o.eq.1) then
              c=ii*sqrt(2.)/9.
           else if (o.eq.2) then
              c=0.
           end if
        case ('j')
           if (o.eq.1) then
              c=ii*sqrt(2.)/9.
           else if (o.eq.2) then
              c=0.
           end if
        case ('k')
           if (o.eq.1) then
              c=-2.*ii*sqrt(2.)/9.
           else if (o.eq.2) then
              c=ii*sqrt(2.)/9.
           end if
        case ('l')
           if (o.eq.1) then
              c=ii*sqrt(2.)/3.
           else if (o.eq.2) then
              c=-ii*sqrt(2.)/3.
           end if
        case ('m')
           if (o.eq.1) then
              c=ii*sqrt(2.)/9.
           else if (o.eq.2) then
              c=-ii*sqrt(2.)/9.
           end if
        case ('o')
           if (o.eq.1) then
              c=-2.*sqrt(2.)/3.
           else if (o.eq.2) then
              c=sqrt(2.)/3.
           end if
        case ('p')
           if (o.eq.1) then
              c=0.
           else if (o.eq.2) then
              c=-ii*sqrt(2.)/9.
           end if
        case ('q')
           if (o.eq.1) then
              c=-sqrt(2.)/9.
           else if (o.eq.2) then
              c=sqrt(2.)/9.
           end if
        case ('r')
           if (o.eq.1) then
              !                  c=0.    !el correcto
              c=2.
           else if (o.eq.2) then
              c=0.
           end if
        case ('s')
           if (o.eq.1) then
              c=-sqrt(2.)*ii
           else if (o.eq.2) then
              c=0.
           end if
        case ('t')
           if (o.eq.1) then
              c=0.
           else if (o.eq.2) then
              c=sqrt(2.)*ii
           end if
        case ('u')
           if (o.eq.1) then
              c=2.*ii/sqrt(3.)
           else if (o.eq.2) then
              c=-ii/2./sqrt(3.)
           end if
        case ('v')
           if (o.eq.1) then
              c=ii*sqrt(2.)
           else if (o.eq.2) then
              c=0.
           end if
        case ('w')
           if (o.eq.1) then
              c=1.
           else if (o.eq.2) then
              c=0.
           end if
        case ('x')
           if (o.eq.1) then
              c=ii*sqrt(1./3.)
           else if (o.eq.2) then
              c=0.
           end if
        end select

     case (3)  ! gamma p -> p pi0 pi0
        select case (diag)
        case ('a')
           if (o.eq.1) then
              c=0.
           else if (o.eq.2) then
              c=0.
           end if
        case ('b')
           if (o.eq.1) then
              c=0.
           else if (o.eq.2) then
              c=0.
           end if
        case ('c')
           if (o.eq.1) then
              c=0.
           else if (o.eq.2) then
              c=0.
           end if
        case ('d')
           if (o.eq.1) then
              c=0.
           else if (o.eq.2) then
              c=0.
           end if
        case ('e')
           if (o.eq.1) then
              c=ii
           else if (o.eq.2) then
              c=ii
           end if
        case ('f')
           if (o.eq.1) then
              c=ii
           else if (o.eq.2) then
              c=ii
           end if
        case ('g')
           if (o.eq.1) then
              c=ii
           else if (o.eq.2) then
              c=ii
           end if
        case ('h')
           if (o.eq.1) then
              c=-2./9.*ii
           else if (o.eq.2) then
              c=-2./9.*ii
           end if
        case ('i')
           if (o.eq.1) then
              c=0.
           else if (o.eq.2) then
              c=0.
           end if
        case ('j')
           if (o.eq.1) then
              c=0.
           else if (o.eq.2) then
              c=0.
           end if
        case ('k')
           if (o.eq.1) then
              c=-2.*ii/9.
           else if (o.eq.2) then
              c=-2.*ii/9.
           end if
        case ('l')
           if (o.eq.1) then
              c=ii*2./3.
           else if (o.eq.2) then
              c=ii*2./3.
           end if
        case ('m')
           if (o.eq.1) then
              c=ii*2./9.
           else if (o.eq.2) then
              c=ii*2./9.
           end if
        case ('o')
           if (o.eq.1) then
              c=-1./3.
           else if (o.eq.2) then
              c=-1./3.
           end if
        case ('p')
           if (o.eq.1) then
              c=2.*ii/9.
           else if (o.eq.2) then
              c=2.*ii/9.
           end if
        case ('q')
           if (o.eq.1) then
              c=-2./9.
           else if (o.eq.2) then
              c=-2./9.
           end if
        case ('r')
           if (o.eq.1) then
              c=2.
           else if (o.eq.2) then
              c=0.
           end if
        case ('s')
           if (o.eq.1) then
              c=0.
           else if (o.eq.2) then
              c=0.
           end if
        case ('t')
           if (o.eq.1) then
              c=0.
           else if (o.eq.2) then
              c=0.
           end if
        case ('u')
           if (o.eq.1) then
              c=ii/sqrt(6.)
           else if (o.eq.2) then
              c=ii/sqrt(6.)
           end if
        case ('v')
           if (o.eq.1) then
              c=0.
           else if (o.eq.2) then
              c=0.
           end if
        case ('w')
           if (o.eq.1) then
              c=0.
           else if (o.eq.2) then
              c=0.
           end if
        case ('x')
           if (o.eq.1) then
              c=0.
           else if (o.eq.2) then
              c=0.
           end if
        end select

     case (4)  ! gamma n -> n pi+ pi-
        select case (diag)
        case ('a')
           if (o.eq.1) then
              c=0.
           else if (o.eq.2) then
              c=-2.*ii
           end if
        case ('b')
           if (o.eq.1) then
              c=0.
           else if (o.eq.2) then
              c=2.*ii
           end if
        case ('c')
           if (o.eq.1) then
              c=0.
           else if (o.eq.2) then
              c=-2.*ii
           end if
        case ('d')
           if (o.eq.1) then
              c=0.
           else if (o.eq.2) then
              c=2.*ii
           end if
        case ('e')
           if (o.eq.1) then
              c=0.
           else if (o.eq.2) then
              c=2.*ii
           end if
        case ('f')
           if (o.eq.1) then
              c=0.
           else if (o.eq.2) then
              c=2.*ii
           end if
        case ('g')
           if (o.eq.1) then
              c=0.
           else if (o.eq.2) then
              c=2.*ii
           end if
        case ('h')
           if (o.eq.1) then
              c=0.
           else if (o.eq.2) then
              c=-2.*ii/9.
           end if
        case ('i')
           if (o.eq.1) then
              c=ii/3.
           else if (o.eq.2) then
              c=-ii/9.
           end if
        case ('j')
           if (o.eq.1) then
              c=ii/3.
           else if (o.eq.2) then
              c=-ii/9.
           end if
        case ('k')
           if (o.eq.1) then
              c=0.
           else if (o.eq.2) then
              c=2.*ii/9.
           end if
        case ('l')
           if (o.eq.1) then
              c=ii
           else if (o.eq.2) then
              c=ii/3.
           end if
        case ('m')
           if (o.eq.1) then
              c=ii/3.
           else if (o.eq.2) then
              c=ii/9.
           end if
        case ('o')
           if (o.eq.1) then
              c=-1.
           else if (o.eq.2) then
              c=2./3.
           end if
        case ('p')
           if (o.eq.1) then
              c=ii/3.
           else if (o.eq.2) then
              c=ii/9.
           end if
        case ('q')
           if (o.eq.1) then
              c=-1./3.
           else if (o.eq.2) then
              c=-1./9.
           end if
        case ('r')
           if (o.eq.1) then
              c=2.
           else if (o.eq.2) then
              c=0.
           end if
        case ('s')
           if (o.eq.1) then
              c=0.
           else if (o.eq.2) then
              c=-2.*ii
           end if
        case ('t')
           if (o.eq.1) then
              c=0.
           else if (o.eq.2) then
              c=2.*ii
           end if
        case ('u')
           if (o.eq.1) then
              c=ii*sqrt(3./2.)
           else if (o.eq.2) then
              c=-ii*sqrt(2./3.)
           end if
        case ('v')
           if (o.eq.1) then
              c=ii
           else if (o.eq.2) then
              c=0.
           end if
        case ('w')
           if (o.eq.1) then
              c=0.
           else if (o.eq.2) then
              c=0.
           end if
        case ('x')
           if (o.eq.1) then
              c=-ii*sqrt(2./3.)
           else if (o.eq.2) then
              c=0.
           end if
        end select

     case (5)  ! gamma n -> p pi0 pi-
        select case (diag)
        case ('a')
           if (o.eq.1) then
              c=0.
           else if (o.eq.2) then
              c=-sqrt(2.)*ii
           end if
        case ('b')
           if (o.eq.1) then
              c=sqrt(2.)*ii
           else if (o.eq.2) then
              c=0.
           end if
        case ('c')
           if (o.eq.1) then
              c=0.
           else if (o.eq.2) then
              c=-sqrt(2.)*ii
           end if
        case ('d')
           if (o.eq.1) then
              c=sqrt(2.)*ii
           else if (o.eq.2) then
              c=0.
           end if
        case ('e')
           if (o.eq.1) then
              c=-sqrt(2.)*ii
           else if (o.eq.2) then
              c=sqrt(2.)*ii
           end if
        case ('f')
           if (o.eq.1) then
              c=-sqrt(2.)*ii
           else if (o.eq.2) then
              c=sqrt(2.)*ii
           end if
        case ('g')
           if (o.eq.1) then
              c=-sqrt(2.)*ii
           else if (o.eq.2) then
              c=sqrt(2.)*ii
           end if
        case ('h')
           if (o.eq.1) then
              c=-2.*ii*sqrt(2.)/9.
           else if (o.eq.2) then
              c=-ii*sqrt(2.)/9.
           end if
        case ('i')
           if (o.eq.1) then
              c=0.
           else if (o.eq.2) then
              c=ii*sqrt(2.)/9.
           end if
        case ('j')
           if (o.eq.1) then
              c=0.
           else if (o.eq.2) then
              c=ii*sqrt(2.)/9.
           end if
        case ('k')
           if (o.eq.1) then
              c=ii*sqrt(2.)/9.
           else if (o.eq.2) then
              c=-2.*ii*sqrt(2.)/9.
           end if
        case ('l')
           if (o.eq.1) then
              c=ii*sqrt(2.)/3.
           else if (o.eq.2) then
              c=-ii*sqrt(2.)/3.
           end if
        case ('m')
           if (o.eq.1) then
              c=ii*sqrt(2.)/9.
           else if (o.eq.2) then
              c=-ii*sqrt(2.)/9.
           end if
        case ('o')
           if (o.eq.1) then
              c=sqrt(2.)/3.
           else if (o.eq.2) then
              c=-2.*sqrt(2.)/3.
           end if
        case ('p')
           if (o.eq.1) then
              c=0.
           else if (o.eq.2) then
              c=-ii*sqrt(2.)/9.
           end if
        case ('q')
           if (o.eq.1) then
              c=-sqrt(2.)/9.
           else if (o.eq.2) then
              c=sqrt(2.)/9.
           end if
        case ('r')
           if (o.eq.1) then
              !                  c=0.    !el correcto
              c=2.
           else if (o.eq.2) then
              c=0.
           end if
        case ('s')
           if (o.eq.1) then
              c=0.
           else if (o.eq.2) then
              c=-sqrt(2.)*ii
           end if
        case ('t')
           if (o.eq.1) then
              c=sqrt(2.)*ii
           else if (o.eq.2) then
              c=0.
           end if
        case ('u')
           if (o.eq.1) then
              c=-ii/2./sqrt(3.)
           else if (o.eq.2) then
              c=ii*2./sqrt(3.)
           end if
        case ('v')
           if (o.eq.1) then
              c=ii*sqrt(2.)
           else if (o.eq.2) then
              c=0.
           end if
        case ('w')
           if (o.eq.1) then
              c=-1.
           else if (o.eq.2) then
              c=0.
           end if
        case ('x')
           if (o.eq.1) then
              c=ii*sqrt(1./3.)
           else if (o.eq.2) then
              c=0.
           end if
        end select

     case (6)  ! gamma n -> n pi0 pi0
        select case (diag)
        case ('a')
           if (o.eq.1) then
              c=0.
           else if (o.eq.2) then
              c=0.
           end if
        case ('b')
           if (o.eq.1) then
              c=0.
           else if (o.eq.2) then
              c=0.
           end if
        case ('c')
           if (o.eq.1) then
              c=0.
           else if (o.eq.2) then
              c=0.
           end if
        case ('d')
           if (o.eq.1) then
              c=0.
           else if (o.eq.2) then
              c=0.
           end if
        case ('e')
           if (o.eq.1) then
              c=ii
           else if (o.eq.2) then
              c=ii
           end if
        case ('f')
           if (o.eq.1) then
              c=ii
           else if (o.eq.2) then
              c=ii
           end if
        case ('g')
           if (o.eq.1) then
              c=ii
           else if (o.eq.2) then
              c=ii
           end if
        case ('h')
           if (o.eq.1) then
              c=2./9.*ii
           else if (o.eq.2) then
              c=2./9.*ii
           end if
        case ('i')
           if (o.eq.1) then
              c=0.
           else if (o.eq.2) then
              c=0.
           end if
        case ('j')
           if (o.eq.1) then
              c=0.
           else if (o.eq.2) then
              c=0.
           end if
        case ('k')
           if (o.eq.1) then
              c=2.*ii/9.
           else if (o.eq.2) then
              c=2.*ii/9.
           end if
        case ('l')
           if (o.eq.1) then
              c=ii*2./3.
           else if (o.eq.2) then
              c=ii*2./3.
           end if
        case ('m')
           if (o.eq.1) then
              c=ii*2./9.
           else if (o.eq.2) then
              c=ii*2./9.
           end if
        case ('o')
           if (o.eq.1) then
              c=1./3.
           else if (o.eq.2) then
              c=1./3.
           end if
        case ('p')
           if (o.eq.1) then
              c=0.
           else if (o.eq.2) then
              c=0.
           end if
        case ('q')
           if (o.eq.1) then
              c=-2./9.
           else if (o.eq.2) then
              c=-2./9.
           end if
        case ('r')
           if (o.eq.1) then
              c=2.
           else if (o.eq.2) then
              c=0.
           end if
        case ('s')
           if (o.eq.1) then
              c=0.
           else if (o.eq.2) then
              c=0.
           end if
        case ('t')
           if (o.eq.1) then
              c=0.
           else if (o.eq.2) then
              c=0.
           end if
        case ('u')
           if (o.eq.1) then
              c=-ii/sqrt(6.)
           else if (o.eq.2) then
              c=-ii/sqrt(6.)
           end if
        case ('v')
           if (o.eq.1) then
              c=0.
           else if (o.eq.2) then
              c=0.
           end if
        case ('w')
           if (o.eq.1) then
              c=0.
           else if (o.eq.2) then
              c=0.
           end if
        case ('x')
           if (o.eq.1) then
              c=0.
           else if (o.eq.2) then
              c=0.
           end if
        end select
     end select

   end function c

   subroutine photoff(r,diag,o,f1,gm)
     !     determines the nucleon ff at the photon point deppending on the
     !     reaction and diagram
     integer,intent(in)::r  ! reaction
     integer,intent(in)::o  ! order in which pions are emited
     character (len=1),intent(in):: diag  !diagram label in appendix 4
     real,intent(out)::f1,gm

     if ((diag.eq.'e').or.(diag.eq.'m')) then
        if ((r.eq.1).or.(r.eq.2).or.(r.eq.3)) then
           !         on proton
           f1=1.
           gm=mup
        else
           !         on neutron
           f1=0.
           gm=mun
        end if
     else if (diag.eq.'g') then
        if ((r.eq.1).or.(r.eq.3).or.(r.eq.5)) then
           !         on proton
           f1=1.
           gm=mup
        else
           !         on neutron
           f1=0.
           gm=mun
        end if
     else if (diag.eq.'f') then
        select case (r)
        case (1)
           !         on neutron for o=1 (o=2) is not allowed and c=0
           f1=0.
           gm=mun
        case (2)
           if (o.eq.1) then
              f1=0.
              gm=mun
           else if (o.eq.2) then
              f1=1.
              gm=mup
           end if
        case (3)
           !           always proton
           f1=1.
           gm=mup
        case (4)
           !           on proton for o=2 (o=1) is not allowed and c=0
           f1=1.
           gm=mup
        case (5)
           if (o.eq.1) then
              f1=0.
              gm=mun
           else if (o.eq.2) then
              f1=1.
              gm=mup
           end if
        case (6)
           !           always neutron
           f1=0.
           gm=mun
        end select
     end if
   end subroutine photoff

   subroutine photodel(r,o,ed,mud)
     !     Delta's charge and magnetic moment (for diag. p)
     integer,intent(in)::r  ! reaction
     integer,intent(in)::o  ! order in which pions are emitted
     real,intent(out)::ed,mud

     select case (r)
     case (1)
        if (o.eq.1) then
           ed=0.
        else
           ed=2.
        end if
     case (2)
        if (o.eq.1) then
           ed=0.
        else
           ed=1.
        end if
     case (3)
        ed=1.
     case (4)
        if (o.eq.1) then
           ed=-1.
        else
           ed=1.
        end if
     case (5)
        if (o.eq.1) then
           ed=0.
        else
           ed=1.
        end if
     case (6)
        ed=0.
     end select

     if (ed.gt.1.) then
        mud=1.62*mup  ! for the Delta++ the exp. value is taken
     else
        mud=ed*mup
     end if
     ed=ed*e

   end subroutine photodel

   subroutine photons(r,g1,g2,g3,f1rop,f2rop)
     !     Fixes the values of the N-N* transitions ff(q2=0)
     !     for N*(1520) and  N*(1440) depending on the target (p or n)
     integer,intent(in)::r
     real,intent(out)::g1,g2,g3,f1rop,f2rop

     if (r.le.3) then
        !       on proton
        g1=g1p
        g2=g2p
        g3=g3p
        f1rop=f1ropp
        f2rop=f2ropp
     else
        !       on neutron
        g1=g1n
        g2=g2n
        g3=g3n
        f1rop=f1ropn
        f2rop=f2ropn
     end if

   end subroutine photons

 end function ampli


end module amplitudes_2Pi
