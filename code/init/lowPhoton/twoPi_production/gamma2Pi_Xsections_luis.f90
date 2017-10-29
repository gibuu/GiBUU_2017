!******************************************************************************
!****m* /gamma2Pi_Xsections_luis
! NAME
! module gamma2Pi_Xsections_luis
! PURPOSE
! Includes the total cross sections for gamma Nucleon -> Nucleon Pion Pion processes according to Tejedor.
!******************************************************************************
module gamma2Pi_Xsections_luis

  implicit none
  private

  logical, parameter :: debugflag=.false.

  public :: minva

  integer,save :: precision=3       ! Precision of integrals, higher value means better precision
  logical,save :: initFlag=.true.
  logical,save :: debug=.false.

contains


  subroutine init

    use output

    integer :: ios

    NAMELIST /gamma_2Pi_Xsection_luis/ precision

    call Write_ReadingInput('gamma_2Pi_Xsection_luis',0)

    rewind(5)
    read(5,nml=gamma_2Pi_Xsection_luis,IOSTAT=IOS)
    call Write_ReadingInput('gamma_2Pi_Xsection_luis',0,IOS)
    write(*,*) 'Precision       :' , precision
    call Write_ReadingInput('gamma_2Pi_Xsections_luis',1)

  end subroutine init


  !****************************************************************************
  !****f* gamma2Pi_Xsections_luis/minva
  ! NAME
  ! function minva(photonEnergy,channel,betaToLRF,mediumAtPosition,position) Result(totalXsection)
  ! INPUTS
  ! * real, intent(in) :: photonEnergy [GEV]
  ! * integer, intent(in) :: channel
  !   * 1: gamma p -> pi+ pi- p
  !   * 2: gamma p -> pi+ pi0 n
  !   * 3: gamma p -> pi0 pi0 p
  !   * 4: gamma n -> pi+ pi- n
  !   * 5: gamma n -> pi- pi0 p
  !   * 6: gamma n -> pi0 pi0 n
  ! * type(medium) ,intent(in) :: mediumAtPosition
  ! * real, dimension(1:3) ,intent(in) :: betaToLRF -- velocity of LRF
  ! * real, dimension(1:3) ,intent(in) :: position
  ! OUTPUT
  ! real :: totalXsection -- in [mb]
  !
  ! PURPOSE
  ! * Photoproduction of 2pi pairs according to Tejedor, Gomez et al.
  ! * Integration over 3-dim phase space of outgoing nucleon and pion pair at fixed photon energy
  !****************************************************************************
  function minva(photonEnergy,channel,betaToLRF,mediumAtPosition,position) Result(totalXsection)
    use constants, only: pi
    use amplitudes_2Pi, only: dos, ampli
    use paramamp  ! masses, coupling constants, etc
    use mediumDefinition
    use gauss_integration, only: sg20r, rg20r

    real, intent(in) :: photonEnergy
    integer, intent(in) :: channel
    type(medium) ,intent(in) :: mediumAtPosition
    real, dimension(1:3) ,intent(in) :: betaToLRF
    real, dimension(1:3) ,intent(in) :: position


    type (dos)::t  ! amplitude
    type (dos)::td,ti

    real, dimension(0:3)::q  ! photon 4-momentum
    real, dimension(0:3)::p1 ! incoming nucleon momentum
    real, dimension(0:3)::p2 ! outgoing nucleon momentum
    real, dimension(0:3)::p4 ! outgoing pi(- or 0) momentum
    real, dimension(0:3)::p5 ! outgoing pi(+ or 0) momentum
    integer:: r              ! identifies the reaction channel:
    !                                r=1: gamma p -> pi+ pi- p
    !                                r=2: gamma p -> pi+ pi0 n
    !                                r=3: gamma p -> pi0 pi0 p
    !                                r=4: gamma n -> pi+ pi- n
    !                                r=5: gamma n -> pi- pi0 p
    !                                r=6: gamma n -> pi0 pi0 n


    integer::j,k,l,n1,n2,ns2,n3,ns3,n4,ns4!,i,ns1
    real::sb,minv,minv_min,minv_max,step,s,smax, &
         &       cos45,sin45,cos4,sin4,phi45,ampl2,resu,csec,lam,wlab!,alph,det,wlab_max,wlab_min
    real::e4,e4min,e4max,e5,p4m,p5m!,e5min,e5max
    real,dimension(2,2)::tsq
    real,dimension(0:3,0:3)::rot
    real, dimension(0:3)::p5p
    real,dimension(:),allocatable::absi2,orde2, &
         &       absi3,orde3,absi4,orde4!,absi1,orde1
    real,parameter::mpi0=134.9764/hb  ! pi0 mass
    !      real,parameter::mpi0=139.6/hb
    !integer :: stepper
    real :: totalXsection
    !      If(debugFlag)  write(*,*)'In minva:',photonEnergy,channel,betaToLRF, mediumAtPosition,position


    if (initFlag) then
       call init
       initFlag=.false.
    end if

    totalXsection=0.
    ! FOR TESTING :::
    ! if(channel.ne.3) return


    !      open(unit=1,file='minv_pi0pi0_460_pi00.dat',status='unknown')

    r=channel   ! choosing the reaction channel
    !     Symmetry factor
    if ((r.eq.3).or.(r.eq.6)) then
       sb=0.5
    else
       sb=1.
    end if

    !     definitions for the integrals
    n1=20*precision  ! precision of the gauss integrals
    n2=precision
    n3=precision
    n4=precision
    !allocate (absi1(20*n1))
    !allocate (orde1(20*n1))
    allocate (absi2(20*n2))
    allocate (orde2(20*n2))
    allocate (absi3(20*n3))
    allocate (orde3(20*n3))
    allocate (absi4(20*n4))
    allocate (orde4(20*n4))
    call sg20r(0.,2.*pi,n3,absi3,ns3)
    call sg20r(-1.,1.,n4,absi4,ns4)


    !     Photon energy interval (in Lab)
    wlab=photonEnergy/hb*1000.

    smax=mn*(mn+2.*wlab)

    minv_min=2.*mpi0
    minv_max=sqrt(smax)-mn
    step=(minv_max-minv_min)/n1

    !      do minv=minv_min,minv_max,step

    !      do stepper=0,20
    minv=minv_min
    do while (minv<=minv_max)

       s=mn*(mn+2.*wlab)

       !       photon 4-momentum
       q(0)=wlab*mn/sqrt(s)
       q(1)=0.
       q(2)=0.
       q(3)=q(0)
       !       proton 4-momentum
       p1(0)=sqrt(mn**2+q(3)**2)
       p1(1)=0.
       p1(2)=0.
       p1(3)=-q(3)

       !       Phase space integrals
       !       Limits in e4
       lam=minv**4+mn**4+s**2-2.*(s*minv**2+mn**2*s+mn**2*minv**2)
       e4min=max(mpi0,(minv*(minv**2-mn**2+s)- &
            &         sqrt(max(0.,(minv**2-4.*mpi0**2)*lam)))/(4.*minv*sqrt(s)))
       e4max=min((s-mn**2-2.*mn*mpi0)/2./sqrt(s),(minv*(minv**2-mn**2+s)+ &
            &         sqrt(max(0.,(minv**2-4.*mpi0**2)*lam)))/(4.*minv*sqrt(s)))

       call sg20r(e4min,e4max,n2,absi2,ns2)
       do j=1,ns2
          e4=absi2(j)
          p4m=sqrt(e4**2-mpi0**2)

          e5=(minv**2+s-mn**2-2.*sqrt(s)*e4)/2./sqrt(s)
          p5m=sqrt(e5**2-mpi0**2)

          cos45=(s-mn**2+2.*mpi0**2-2.*sqrt(s)*(e4+e5)+2.*e4*e5)/(2.*p4m*p5m)
          if (abs(cos45).gt.1.) then
             if (debug) write(*,*) 'la hemos cagao cos45=',cos45
             sin45=0.
             orde2(j)=0.
             cycle
          else
             sin45=sqrt(1.-cos45**2)
          end if
          do k=1,ns3
             phi45=absi3(k)

             do l=1,ns4
                cos4=absi4(l)
                sin4=sqrt(1.-cos4**2)

                !             4-momentum of one of the outgoing pions
                p4(0)=e4
                p4(1)=p4m*sin4
                p4(2)=0.
                p4(3)=p4m*cos4

                rot=reshape((/1.,0.,0.,0.,0.,cos4,0.,-sin4, &
                     &               0.,0.,1.,0.,0.,sin4,0.,cos4/),(/4,4/))

                !             4-momentum of the other pion in the rotated frame (with z along p4)
                p5p(0)=e5
                p5p(1)=p5m*sin45*cos(phi45)
                p5p(2)=p5m*sin45*sin(phi45)
                p5p(3)=p5m*cos45

                p5=matmul(rot,p5p)

                !             the nucleon 4-momentum
                p2=q+p1-p4-p5

                td = ampli(r,q,p1,p2,p4,p5,1,betaToLRF,mediumAtPosition,position)
                ti = ampli(r,q,p1,p2,p5,p4,2,betaToLRF,mediumAtPosition,position)

                t%t1=td%t1+ti%t1
                t%t2=td%t2+ti%t2

                !              If(debugFlag)  write(*,*)'In minva ampli:',t
                !              If(debugFlag)  write(*,*)'In minva ampli td:',td
                !              If(debugFlag)  write(*,*)'In minva ampli ti:',ti
                !             ampl2: amplitude**2, summed and averaged over spins

                tsq=matmul(t%t1,transpose(conjg(t%t1)))+ &
                     &               matmul(t%t2,transpose(conjg(t%t2)))

                ampl2=1./4.*(tsq(1,1)+tsq(2,2))

                orde4(l)=ampl2*(minv/sqrt(s))
             end do
             call rg20r(-1.,1.,n4,orde4,orde3(k))
             !            call sumUp(-1.,1.,orde4,orde3(k))
          end do
          call rg20r(0.,2.*pi,n3,orde3,orde2(j))
          !          call sumUp(0.,2.*pi,orde3,orde2(j))
       end do
       call rg20r(e4min,e4max,n2,orde2,resu)
       !        call sumUp(e4min,e4max,orde2,resu)

       csec=mn**2/(s-mn**2)*sb/4./(2.*pi)**4*resu
       if (debug) write(111,*)minv*hb,csec/hb*fc*1000.
       !        If(debugFlag)  write(*,*)'In minva:',minv*hb,csec,csec/hb*fc*1000.,resu, s, mn, sb,pi

       totalXsection=totalXsection+csec!*minv*hb
       !      *(197.33)**2

       minv=minv+step

    end do
    totalXsection=step*totalXsection*fc

  end function minva

!!$  subroutine sumUp(a,b,array,result)
!!$    implicit none
!!$    real, intent(in)  :: a,b
!!$    real, intent(in) , dimension(:) :: array
!!$    real, intent(out) :: result
!!$    result=(b-a)/float(size(array))*sum(array)
!!$  end subroutine sumUp

end module gamma2Pi_Xsections_luis
