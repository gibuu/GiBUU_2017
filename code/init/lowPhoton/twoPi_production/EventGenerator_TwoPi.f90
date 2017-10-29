!******************************************************************************
!****m* /EventGenerator_TwoPi
! NAME
! module EventGenerator_TwoPi
! PURPOSE
! Includes the event generator for gamma Nucleon -> Nucleon Pion Pion processes.
! AUTHOR
! Pascal Muehlich
!******************************************************************************
module EventGenerator_TwoPi

  implicit none
  private

  public :: EventGenerator

contains

  !****************************************************************************
  !****s* EventGenerator_TwoPi/EventGenerator
  ! NAME
  ! subroutine EventGenerator(r,pi1,pi2,po1,po2,po3)
  ! INPUTS
  ! * r   = reaction channel (see Amplitudes.f90)
  ! * pi1 = incomming photon 4-momentum in nuclear rest frame
  ! * pi2 = incomming nucleon 4-momentum "    "      "     "
  ! OUTPUT
  ! * po1 = outgoing nucleon 4-momentum in nuclear rest frame
  ! * po2 =     "    pion1
  ! * po3 =     "    pion2
  !****************************************************************************
  subroutine EventGenerator(r,pi1,pi2,po1,po2,po3,mediumAtPosition,position,noBoostToLab)
    use lorentzTrafo, only: lorentz
    use random, only: rn
    use amplitudes_2Pi, only: ampli, dos
    use mediumDefinition
    use minkowski, only: abs4,sp
    use constants, only: pi, hbarc, mPi, mN
    use rotation, only: rotateTo

    integer, intent(in) :: r
    real,dimension(0:3) :: pi1,pi2
    real,dimension(0:3),intent(out) :: po1,po2,po3
    type(medium), intent(in) :: mediumAtPosition
    real, dimension(1:3), intent(in) :: position
    logical, intent(in) :: noBoostToLab

    real,dimension(1:3) :: beta,betaDummy
    integer :: niter, numTries
    real::w4,w5,w4min,w4max,w5min,w5max,s,srts!,alph,det,egamma
    real::cost45,sint45,p4m,p5m,phi45,cost4,sint4,phir,w5minus,w5plus
    real::m4,m5,s25,s24p,s24m!,m3
    real,dimension(0:3)::k,p1,p2,p4,p5,p5p,ptot,p4r,p5r!,p2r
    real,dimension(0:3,0:3)::rot
    type(dos) :: t, td, ti
    real :: tsq(2,2),ampl2,proba,a2max,eph
    logical::flag,flag2
    integer,parameter::nitermax=1500
    logical,parameter :: debug=.false.
    logical, parameter :: pascalBug=.false. ! Never turn this to true, except for testing!!!

    real, parameter :: mn_fm=mN/hbarc, mpi_fm=mPi/hbarc

    if (debug) then
       write(*,*) r,pi1,pi2
    end if

    pi1=pi1/hbarc
    pi2=pi2/hbarc
    ptot=pi1+pi2
    if (pascalbug) then
       s=ptot(0)**2-ptot(1)**1-ptot(2)**2-ptot(3)**2
       srts=sqrt(s)
    else
       s=sp(ptot,ptot)
       srts=abs4(ptot)
    end if
    beta(1:3)=ptot(1:3)/ptot(0)
    betaDummy=-beta
    eph=(s-mn_fm**2)/(2.*mn_fm)
    call amax(r,eph,a2max)
    !boost to cm system
    k=pi1
    call lorentz(beta(1:3),k(0:3))
    p1(1:3)=-k(1:3)
    p1(0)=sqrt(mn_fm**2+p1(1)**2+p1(2)**2+p1(3)**2)
    !masses of outgoing particles
    m4=mpi_fm
    m5=mpi_fm
    !m3=mn_fm
    !maximum phase-space
    flag=.true.
    niter=0
    do while(flag)
       niter=niter+1
       !start monte carlo
       w4min=m4
       w4max=(s+m4**2-(mn_fm+m5)**2)/2./srts
       w5min=m5
       w5max=(s+m5**2-(mn_fm+m4)**2)/2./srts
       flag2=.true.
       !choose energies
       numTries=1
       do while(flag2)
          w4=rn()*(w4max-w4min)+w4min
          w5=rn()*(w5max-w5min)+w5min
          s25=s+m4**2-2.*srts*w4
          if (s25.eq.0.or.kallen(s25,s,m4**2).lt.0.or.kallen(s25,s,m5**2).lt.0) then
             write(*,'(A,I10,3F12.5)') "In event generator for 2Pi, try number: ", &
                  & numTries,s25,kallen(s25,s,m4**2),kallen(s25,s,m5**2)
             numTries=numTries+1
             write(*,*) w4min, w4max
             write(*,*) w5min, w5max
             cycle
          end if
          s24p=mn_fm**2+m4**2-1./2./s25*((s25-s+m4**2)*(s25+mn_fm**2-m5**2)- sqrt(kallen(s25,s,m4**2)) &
               & *sqrt(kallen(s25,mn_fm**2,m5**2)))
          s24m=mn_fm**2+m4**2-1./2./s25*((s25-s+m4**2)*(s25+mn_fm**2-m5**2)+ sqrt(kallen(s25,s,m4**2)) &
               & *sqrt(kallen(s25,mn_fm**2,m5**2)))
          w5minus=(s+m5**2-s24p)/2./srts
          w5plus=(s+m5**2-s24m)/2./srts
          if (w5>=w5minus.and.w5<=w5plus)flag2=.false.
       end do
       p4m=sqrt(w4**2-m4**2)
       p5m=sqrt(w5**2-m5**2)
       cost45=(s-mn_fm**2+m4**2+m5**2-2.*srts*(w4+w5)+2.*w4*w5)  /(2.*p4m*p5m)
       if (abs(cost45)>1.0) then
          write(*,*) 'problems cost45',cost45
          cost45=sign(1.,cost45)
       end if
       sint45=sqrt(1.-cost45**2)
       !choose directions
       phi45=rn()*2.*pi
       cost4=max(min(2.*(rn()-0.5),1.),-1.)
       sint4=sqrt(1.-cost4**2)
       p4r(0)=w4
       p4r(1)=p4m*sint4
       p4r(2)=0.
       p4r(3)=p4m*cost4
       rot=reshape((/1.,0.,0.,0.,0.,cost4,0.,-sint4,  0.,0.,1.,0.,0.,sint4,0.,cost4/),(/4,4/))
       p5p(0)=w5
       p5p(1)=p5m*sint45*cos(phi45)
       p5p(2)=p5m*sint45*sin(phi45)
       p5p(3)=p5m*cost45
       p5r=matmul(rot,p5p)
       !random rotation with respect to z-axis
       phir=rn()*2.*pi
       rot=reshape((/1.,0.,0.,0.,  0.,cos(phir),sin(phir),0., 0.,-sin(phir),cos(phir),0., 0.,0.,0.,1./),(/4,4/))
       p4=matmul(rot,p4r)
       p5=matmul(rot,p5r)
       ! rotation to incoming photon direction
       p4(1:3) = rotateTo (k(1:3),p4(1:3))
       ! (now vec(k)*vec(p4)=|p4|*|k|*cost4)
       p5(1:3) = rotateTo (k(1:3),p5(1:3))
       p2=k+p1-p4-p5
       !amplitudes
       td = ampli(r,k,p1,p2,p4,p5,1,betaDummy,mediumAtPosition,position)
       ti = ampli(r,k,p1,p2,p5,p4,2,betaDummy,mediumAtPosition,position)
       t%t1=td%t1+ti%t1
       t%t2=td%t2+ti%t2
       tsq=matmul(t%t1,transpose(conjg(t%t1)))+   matmul(t%t2,transpose(conjg(t%t2)))
       ampl2=1./4.*(tsq(1,1)+tsq(2,2))
       proba=ampl2/a2max
       if (proba>1.)write(*,*)'problems proba Event',proba,a2max,ampl2
       if (rn()<proba)flag=.false.
       if (niter==nitermax) then
          write(*,*) 'problems EventGenerator nitermax' ,niter,proba,a2max,srts,p4m,p5m
          flag=.false.
       end if
    end do
    po2=p4
    po3=p5

    if (noBoostToLab) then
       ! We stay in the CM frame
       po1=(/srts,0.,0.,0./)-po2-po3
    else
       ! Boost to lab frame
       call lorentz(betaDummy, po2(0:3))
       call lorentz(betaDummy, po3(0:3))
       po1=ptot-po2-po3
    end if

    !conversion to GeV for BUU
    po1=po1*hbarc
    po2=po2*hbarc
    po3=po3*hbarc

  end subroutine EventGenerator

  !****************************************************************************

  subroutine amax(r,e,a2max)
    use inputGeneral, only: path_to_Input

    integer::r,i!,j
    real::e
    real,intent(out):: a2max
    integer,parameter::n=32
    real,dimension(1:n)::ephot
    real,dimension(1:6,1:n)::maxa
    logical::flag
    save::ephot,maxa
    logical,save :: initFlagAmax=.true.

    if (initFlagAmax) then
       open(10,file=trim(path_to_Input)//'/photo_twoPi/gam2pi_max.dat',ACTION='read')
       do i=1,n
          read(10,*)ephot(i),maxa(1:6,i)
       end do
       close(10)
       initFlagAmax=.false.
    end if

    i=0
    flag=.true.
    if (e>4.97) then
       write(*,*) 'photon energy too large for EventGenerator',   e,ephot(n)
    end if
    do while(flag)
       i=i+1
       if (ephot(i)>=e.or.i==n) then
          flag=.false.
       end if
    end do
    if (i>=2) then
       a2max=1.5*max(maxa(r,i),maxa(r,i-1))
    else if (i==1) then
       a2max=1.5*maxa(r,i)
    else
       write(*,*) 'Error in amax of EventGenerator_TwoPi.f90. i=', i
       stop
    end if

  end subroutine amax

  !****************************************************************************

  real function kallen(x,y,z)
    real::x,y,z
    kallen=(x-y-z)**2-4.*y*z
  end function kallen

end module EventGenerator_TwoPi
