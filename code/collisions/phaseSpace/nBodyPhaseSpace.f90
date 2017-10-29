!******************************************************************************
!****m*  /nBodyPhaseSpace
! NAME
! module nBodyPhaseSpace
! NOTES
! Includes the routine integrate_nBodyPS needed for calculation
! of the n-body phase space volume (n=2,...,6) and also the
! routines momenta_in_nBodyPS, momenta_in_2BodyPS, ..., momenta_in_6BodyPS
! for the simulation of the particle momenta according to the n-body phase
! space.
!******************************************************************************
module nBodyPhaseSpace

  implicit none
  private

  ! Global variables:
  real :: srtsG, qG
  real, dimension(1:6) :: massG

  logical, parameter :: debugFlag=.false.

  public :: momenta_in_nBodyPS, momenta_in_3BodyPS, momenta_in_4BodyPS
!   Public :: momenta_in_2BodyPS, momenta_in_5BodyPS, momenta_in_6BodyPS
  public :: momenta_in_3Body_BYK

  public :: integrate_nBodyPS

contains


  !****************************************************************************
  !****s*  nBodyPhaseSpace/momenta_in_nBodyPS
  ! NAME
  ! subroutine momenta_in_nBodyPS (srts, mass, pn)
  ! PURPOSE
  ! This subroutine determines momenta of outgoing particles according to
  ! n-body phase space (n=2,...,6). Vacuum prescription.
  ! INPUTS
  ! * real :: srts               -- sqrt(s)
  ! * real, dimension(:) :: mass -- masses of particles
  ! OUTPUT
  ! * real, dimension(:,:) :: pn -- three-vectors of momenta
  !   (first index: momentum index, second index: particle)
  ! NOTES
  ! Subroutines momenta_in_2BodyPS, ... , momenta_in_6BodyPS
  ! can be also used directly instead.
  !****************************************************************************
  subroutine momenta_in_nBodyPS (srts, mass, pn)

    real, intent(in) :: srts
    real, intent(in), dimension(:) :: mass
    real, intent(out), dimension(:,:) :: pn

    integer :: n

    n = size(mass,dim=1)

    if (n < 2 .or. n > 6) then
       write(*,*) ' In momenta_in_nBodyPS: size of mass array is wrong: ', n
       stop
    else if (size(pn,dim=2)/=n .or. size(pn,dim=1)/=3) then
       write(*,*) ' In momenta_in_nBodyPS: size of pn array is wrong: ', size(pn,dim=1), size(pn,dim=2)
       stop
    end if

    select case (n)
    case (2)
       if (debugFlag) write(*,*) ' Calling momenta_in_2BodyPS ...'
       pn(1:3,1:2) = momenta_in_2BodyPS (srts, mass(1:2))
       if (debugFlag) write(*,*) ' Done'

    case (3)
       if (debugFlag) write(*,*) ' Calling momenta_in_3BodyPS ...'
       pn(1:3,1:3) = momenta_in_3BodyPS (srts, mass(1:3))
       if (debugFlag) write(*,*) ' Done'

    case (4)
       if (debugFlag) write(*,*) ' Calling momenta_in_4BodyPS ...'
       pn(1:3,1:4) = momenta_in_4BodyPS (srts, mass(1:4))
       if (debugFlag) write(*,*) ' Done'

    case (5)
       if (debugFlag) write(*,*) ' Calling momenta_in_5BodyPS ...'
       pn(1:3,1:5) = momenta_in_5BodyPS (srts, mass(1:5))
       if (debugFlag) write(*,*) ' Done'

    case (6)
       if (debugFlag) write(*,*) ' Calling momenta_in_6BodyPS ...'
       pn(1:3,1:6) = momenta_in_6BodyPS (srts, mass(1:6))
       if (debugFlag) write(*,*) ' Done'

    end select

  end subroutine momenta_in_nBodyPS



  !****************************************************************************
  !****f*  nBodyPhaseSpace/momenta_in_2BodyPS
  ! NAME
  ! function momenta_in_2BodyPS (srts, mass) result (p2)
  ! PURPOSE
  ! This subroutine determines momenta of outgoing particles according to
  ! 2-body phase space. Vacuum prescription.
  ! INPUTS
  ! * real :: srts               -- sqrt(s)
  ! * real, dimension(2) :: mass -- masses of particles
  ! OUTPUT
  ! * real, dimension(3,2) :: p2 -- 2 three-vectors of momenta
  !   (first index: momentum index, second index: particle)
  ! NOTES
  !****************************************************************************
  function momenta_in_2BodyPS (srts, mass) result (p2)

    use random, only: rnOmega
    use twoBodyTools, only: pCM

    real, intent(in) :: srts
    real, intent(in), dimension(2) :: mass
    real, dimension(3,2) :: p2

    if (srts<=sum(mass)) then
      write(*,'(A,2F8.4)') 'problem in momenta_in_2BodyPS: srts too small', srts, sum(mass)
    else
      p2(1:3,1) = pCM(srts,mass(1),mass(2)) * rnOmega()
      p2(1:3,2) = -p2(1:3,1)
    end if

  end function momenta_in_2BodyPS


  !****************************************************************************
  !****f* nBodyPhaseSpace/momenta_in_3BodyPS
  ! NAME
  ! function momenta_in_3BodyPS (srts, mass) result (p3)
  ! PURPOSE
  ! This subroutine determines momenta of outgoing particles according to
  ! 3-body phase space. Vacuum prescription!
  ! INPUTS
  ! * real               :: srts  -- sqrt(s) in the reaction
  ! * real, dimension(3) :: mass  -- Masses of the three particles
  ! OUTPUT
  ! * real,dimension(3,3)  :: p3  -- 3 three-vectors of momenta
  !   (first index: momentum index, second index: particle)
  ! NOTES
  ! Original name in old code: "bopsMom"
  !****************************************************************************
  function momenta_in_3BodyPS (srts, mass) result (p3)

    use rotation, only: rotateYZ
    use random, only: rn
    use constants, only: pi

    real, intent(in) :: srts
    real, intent(in), dimension(3) :: mass
    real, dimension(3,3) :: p3

    real :: m12min,m12max,m23min,m23max,m12,m23,m23ma,m23mi
    real :: e1,e2,e3,pp3(3),cost3,phi,sint3,cost,s,mass2(3)
    integer :: k
    logical :: flag

    if (srts<=sum(mass)) then
       write(*,'(A,4F8.4)')'problem in momenta_in_3BodyPS: srts too small', srts, sum(mass)
       return
    end if

    s = srts**2
    mass2 = mass**2

    m12min=(mass(1)+mass(2))**2
    m12max=(srts-mass(3))**2
    m23min=(mass(2)+mass(3))**2
    m23max=(srts-mass(1))**2

    flag=.true.

    ! Roll the dice for m12 and m23.
    ! We assume that the phase space is equally distributed in the
    ! (m12, m23 ) plane.
    ! The borders of this phase space are evaluated by the routine "ps3bo".

    do while (flag)
       m12=m12min+rn()*(m12max-m12min)
       m23=m23min+rn()*(m23max-m23min)

       call ps3bo(mass(1),mass(2),mass(3),srts,m12,m23ma,m23mi)

       if ((m23>m23mi).and.(m23<m23ma)) flag = .false.
    end do

    ! Define third momentum vector to lie on the z-Axis
    e3=(s+mass2(3)-m12)/(2*srts)
    pp3(3)=sqrt(max(e3**2-mass2(3),0.))
    p3(:,3)=(/0., 0., pp3(3)/)

    e1=(s+mass2(1)-m23)/(2*srts)
    pp3(1)=sqrt(max(e1**2-mass2(1),0.))

    e2=srts-e1-e3
    pp3(2)=sqrt(max(e2**2-mass2(2),0.))

    !*angle between pp3(1) and pp3(3)
    cost3=(-1.)*(pp3(1)**2+pp3(3)**2-pp3(2)**2)/(2.*pp3(1)*pp3(3))

    if (abs(cost3)>1) then
       write(*,*) 'Problem in momenta_in_3BodyPS: |cos|>1',cost3
       stop
       ! Old bugfix:
       ! cost3=sign(1.,cost3)
    end if

    ! Choose the first momentum vector relative to the third one
    phi=rn()*2.*pi
    sint3=sqrt(max(1.-cost3**2,0.))

    p3(:,1)=(/sint3*cos(phi), sint3*sin(phi), cost3/)*pp3(1)

    ! The second momentum vector is now defined by momentum conservation:
    p3(:,2)=-p3(:,1)-p3(:,3)

    !...check:
    if (abs(p3(1,2)**2+p3(2,2)**2+p3(3,2)**2-pp3(2)**2)>1e-05) then
       write(*,*) 'Problems in momenta_in_3BodyPS:', &
            p3(1,2),p3(2,2),p3(3,2),pp3(2),&
            p3(1,2)**2+p3(2,2)**2+p3(3,2)**2-pp3(2)**2
       write(*,*) 'srt(s)=',srts
       stop
    end if

    ! Since we have chosen arbitrarily the z-axis as direction of the
    ! third vector, we now need to do
    ! a random rotation of the whole coordinate system. By doing this
    ! we ensure that we distribute the phase space
    ! in an isotropic way.

    cost=2*rn()-1.0
    phi=rn()*2.*pi

    do k=1,3
      p3(:,k) = rotateYZ (0.0, phi, p3(:,k), cosTheta=cost)
    end do

  contains

    !**************************************************************************
    !****is*  momenta_in_3BodyPS/ps3bo
    ! NAME
    ! subroutine ps3bo(m1,m2,m3,srts,m12,m23max,m23min)
    ! PURPOSE
    ! This subroutine determines the boundaries m23max,m23min in a Dalitz plot
    ! for given srts and m12 (m12=m12**2)
    ! notation from Particle Data Book p.176
    ! INPUTS
    ! * real :: m1,m2,m3 -- masses of the three particles
    ! * real :: srts     -- sqrt(s) in the problem
    ! * real :: m12      -- Mass of particles one and two
    ! OUTPUT
    ! * real :: m23max,m23min
    !**************************************************************************
    subroutine ps3bo(m1,m2,m3,srts,m12,m23max,m23min)

      real, intent(in) ::  m1,m2,m3,srts,m12
      real, intent(out) :: m23max,m23min

      real :: e2star,e3star

      e2star=(m12-m1**2+m2**2)/(2*sqrt(m12))
      e3star=(srts**2-m12-m3**2)/(2*sqrt(m12))
      m23max=(e2star+e3star)**2-(sqrt(max(e2star**2-m2**2,0.))-sqrt(max(e3star**2-m3**2,0.)))**2
      m23min=(e2star+e3star)**2-(sqrt(max(e2star**2-m2**2,0.))+sqrt(max(e3star**2-m3**2,0.)))**2

    end subroutine ps3bo

  end function momenta_in_3BodyPS



  !****************************************************************************
  !****f* nBodyPhaseSpace/momenta_in_3Body_BYK
  ! NAME
  ! function momenta_in_3Body_BYK (srts, pcm_in, mass) result (p3)
  ! PURPOSE
  ! Determine momenta of outgoing Baryon (1), Hyperon (2) and Kaon (3)
  ! in BB->BYK collisions using Kaon angular and momentum disitributions from:
  ! Yu-Ming Zheng et al., PRC 69 (2004) 034907.
  ! INPUTS
  ! * real :: srts         -- sqrt(s) in the reaction
  ! * real :: pcm_in(1:3)  -- CM 3-momentum of one of the incoming particles
  ! * real :: mass(1:3)    -- Masses of the three particles
  ! OUTPUT
  ! * real,dimension(3,3)  :: p3  -- 3 three-vectors of momenta
  !   (first index: momentum index, second index: particle)
  ! NOTES
  ! Original name in old code: "bykmom"
  ! TODO
  ! * 1) Use this routine not only for "N Y K", but also for "Delta Y K".
  ! * 2) Take care of potentials.
  !****************************************************************************
  function momenta_in_3Body_BYK (srts, pcm_in, mass) result (p3)
    use constants, only: pi
    use random, only: rn, rnOmega
    use twoBodyTools, only: pCM
    use minkowski, only: abs4
    use rotation, only: rotateTo
    use lorentzTrafo, only: lorentz

    real, intent(in) :: srts, pcm_in(1:3), mass(1:3)
    real, dimension(3,3) :: p3

    real :: betaby(1:3),x1,x2,pkaon,ct,st,phi,pbyabs
    real, dimension(0:4) :: pK, pby, pbar, phyp
    real, parameter :: a = 1.2

    ! Test for threshold
    if (srts <= sum(mass)) then
      write(*,*) 'problems in momenta_in_3Body_BYK: srts too small', srts, mass(1:3)
      p3(:,:) = 0.
      return
    end if

    ! determine abs. value of kaon momentum in BYK c.m. frame:
    do
      x1 = rn()
      x2 = rn()*0.035
      if (x2 <= x1**3*(1.-x1)**2) exit
    end do
    pkaon = x1* pCM(srts,mass(1)+mass(2),mass(3))

    ! determine kaon momentum components in BYK c.m. frame:
    do
      ct = 1.-2.*rn()
      x2 = (1.+a)*rn()
      if (x2 <= 1.+a*ct**2) exit
    end do
    st = sqrt(1.-ct**2)
    phi = 2.*pi*rn()

    pK(0) = sqrt(pkaon**2+mass(3)**2)
    pK(1:3) = pkaon * (/st*cos(phi),st*sin(phi),ct/)

    ! rotate according to incoming particle
    pK(1:3) = rotateTo (pcm_in(1:3), pK(1:3))

    ! Energy, momentum components and c.m. velocity of the BY pair:
    pby(0) = srts - pK(0)
    pby(1:3) = -pK(1:3)
    betaby(1:3) = pby(1:3)/pby(0)

    ! B and Y momenta in BY c.m. frame:
    pbyabs = pCM (abs4(pby), mass(1), mass(2))
    ! Isotropic distribution in BY c.m. frame:
    pbar(1:3) = rnOmega() * pbyabs
    pbar(0) = sqrt(pbyabs**2+mass(1)**2)
    phyp(1:3) = -pbar(1:3)
    phyp(0) = sqrt(pbyabs**2+mass(2)**2)

    ! Boost to the BYK c.m. frame:
    call lorentz (-betaby(1:3), pbar(:))
    call lorentz (-betaby(1:3), phyp(:))

    ! define output
    p3(1:3,1) = pbar(1:3)
    p3(1:3,2) = phyp(1:3)
    p3(1:3,3) = pK  (1:3)

  end function momenta_in_3Body_BYK



  !****************************************************************************
  !****f* nBodyPhaseSpace/momenta_in_4BodyPS
  ! NAME
  ! function momenta_in_4BodyPS (srts, mass) result (p4)
  ! PURPOSE
  ! This subroutine determines momenta of outgoing particles according to
  ! 4-body phase space. Vacuum prescription
  ! INPUTS
  ! * real :: srts               -- sqrt(s)
  ! * real, dimension(4) :: mass -- masses of the 4 particles
  ! OUTPUT
  ! * real, dimension(3,4) :: p4 -- 4 three-vectors of momenta
  !   (first index: momentum index, second index: particle)
  ! NOTES
  ! Original name in old code: "bopsMom4"
  !****************************************************************************
  function momenta_in_4BodyPS (srts, mass) result (p4)

    use random, only: rn, rnOmega
    use threeBodyPhaseSpace, only: Integrate_3BodyPS
    use lorentzTrafo, only: lorentz

    real, intent(in) :: srts
    real, intent(in), dimension(4) :: mass
    real, dimension(3,4) :: p4

    real :: s,q2max,q2min,pfinalm,pfinal,q2,q,ps3,ps3max,p3(3,3),e123,beta(3),ptot
    real, dimension(0:3) :: FourVector
    integer :: ni,j,k
    integer, parameter :: nimax=500
    logical :: flag
    real, dimension(4) :: mass2

    if (srts<=sum(mass)) then
       write(*,'(A,2F8.4)') 'problem in momenta_in_4BodyPS: srts too small', srts, sum(mass)
       return
    end if

    q2max=(srts-mass(4))**2
    q2min=(mass(1)+mass(2)+mass(3))**2
    s = srts**2
    mass2 = mass**2

    pfinalm=sqrt(max(0., (s+q2min-mass2(4))**2/(4.*s) - q2min ))

    ps3max = Integrate_3bodyPS (sqrt(q2max), mass(1), mass(2), mass(3))

    flag=.true.
    ni=0
    if (q2max>q2min) then
       do while(flag)
          ni=ni+1

          q2=rn()*(q2max-q2min)+q2min
          pfinal=sqrt(max(0., (s+q2-mass2(4))**2/(4.*s) - q2 ))

          q=sqrt(q2)
          ps3 = Integrate_3bodyPS (q, mass(1), mass(2), mass(3))

          if (pfinal*ps3 > rn()*pfinalm*ps3max) flag = .false.

          if (ni>nimax) then
             write(*,*) 'problems in momenta_in_4BodyPS ni.gt.nimax',srts, sum(mass)
             q2=0.5*(q2max-q2min)+q2min
             q=sqrt(q2)
             pfinal=sqrt(max(0., (s+q2-mass2(4))**2/(4.*s) - q2 ))
             flag=.false.
          end if
       end do
    else
       q2=q2min
       q=sqrt(q2)
       pfinal=0.
    end if

    !... select momentum of particle 4:
    p4(1:3,4)=pfinal*rnOmega()

    !... determine momenta of particles 1-3 in their rest frame:
    p3 = momenta_in_3bodyPS (q, mass(1:3))

    !... boost momenta to overall rest frame:
    e123=sqrt(q2+pfinal**2)
    beta(1:3) = p4(1:3,4)/e123

    do k=1,3
       fourVector(1:3)= p3(1:3,k)
       fourVector(0)  = sqrt(mass(k)**2+Dot_Product(fourVector(1:3),fourVector(1:3)))
       call lorentz(beta,fourVector)
       p4(1:3,k) = fourVector(1:3)
    end do

    !... checks:
    do k=1,3
       ptot = sum(p4(k,1:4))
       if (abs(ptot)>1e-05) then
          write(*,*) 'Problems in momenta_in_4BodyPS: momentum conservation',ptot
          stop
       end if
    end do

    ptot=0
    do j=1,4
       ptot=ptot+sqrt(mass(j)**2+p4(1,j)**2+p4(2,j)**2+p4(3,j)**2)
    end do
    if (abs(ptot-srts)>1e-05) then
       write(*,*) 'Problems in momenta_in_4BodyPS: energy conservation',ptot,srts
       stop
    end if

  end function momenta_in_4BodyPS



  !****************************************************************************
  !****f* nBodyPhaseSpace/momenta_in_5BodyPS
  ! NAME
  ! function momenta_in_5BodyPS (srts, mass) result (p5)
  ! PURPOSE
  ! This subroutine determines momenta of outgoing particles according to
  ! 5-body phase space. Vacuum prescription.
  ! INPUTS
  ! * real :: srts               -- sqrt(s)
  ! * real, dimension(5) :: mass -- masses of particles
  ! OUTPUT
  ! * real, dimension(3,5) :: p5 -- 5 three-vectors of momenta
  !   (first index: momentum index, second index: particle)
  ! NOTES
  ! Uses the recursion formula for the phase space (see PDG-2002, p.269)
  ! in a similar way as in subroutine momenta_in_4BodyPS.
  !****************************************************************************
  function momenta_in_5BodyPS (srts, mass) result (p5)

    use random, only: rn, rnOmega
    use lorentzTrafo, only: lorentz

    real, intent(in) :: srts
    real, intent(in), dimension(5) :: mass
    real, dimension(3,5) :: p5

    integer, parameter :: nimax=500  ! maximum number of iterations
    real :: q2max,q2min,s,pfinalm,q2,q,pfinal,ps4,ps4max,e1234,ptot
    integer :: ni,k
    logical :: flag
    real, dimension(5) :: mass2
    real, dimension(1:3) :: beta
    real, dimension(0:3) :: fourVector

    if (srts<=sum(mass)) then
       write(*,'(A,2F8.4)') 'problem in momenta_in_5BodyPS: srts too small', srts, sum(mass)
       return
    end if

    q2max=(srts-mass(5))**2
    q2min=(mass(1)+mass(2)+mass(3)+mass(4))**2
    s = srts**2
    mass2 = mass**2

    pfinalm=sqrt(max(0., (s+q2min-mass2(5))**2/(4.*s) - q2min ))

    ps4max = integrate_nBodyPS (sqrt(q2max), mass(1:4))

    flag=.true.
    ni=0
    if (q2max>q2min) then
       do while(flag)
          ni=ni+1

          q2=rn()*(q2max-q2min)+q2min
          pfinal=sqrt(max(0., (s+q2-mass2(5))**2/(4.*s) - q2 ))

          q=sqrt(q2)
          ps4 = integrate_nBodyPS (q, mass(1:4))

          if (pfinal*ps4 > rn()*pfinalm*ps4max) flag = .false.

          if (ni>nimax) then
             write(*,*) 'problems in momenta_in_5BodyPS ni.gt.nimax',srts, sum(mass)
             q2=0.5*(q2max-q2min)+q2min
             q=sqrt(q2)
             pfinal=sqrt(max(0., (s+q2-mass2(5))**2/(4.*s) - q2 ))
             flag=.false.
          end if
       end do
    else
       q2=q2min
       q=sqrt(q2)
       pfinal=0.
    end if

    !... select momentum of particle 5:
    p5(1:3,5)=pfinal*rnOmega()

    !... determine momenta of particles 1-4 in their rest frame:
    p5(1:3,1:4) = momenta_in_4BodyPS (q, mass(1:4))

    !... boost momenta to overall rest frame:
    e1234=sqrt(q2+pfinal**2)
    beta(1:3) = p5(1:3,5)/e1234

    do k=1,4
       fourVector(1:3)= p5(1:3,k)
       fourVector(0)  = sqrt(mass(k)**2+Dot_Product(fourVector(1:3),fourVector(1:3)))
       call lorentz(beta,fourVector)
       p5(1:3,k) = fourVector(1:3)
    end do

    !... checks:
    do k=1,3
       ptot = sum(p5(k,1:5))
       if (abs(ptot)>1.e-05) then
          write(*,*) 'Problems in momenta_in_5BodyPS: momentum conservation', ptot
          stop
       end if
    end do

    ptot=0.
    do k=1,5
       ptot=ptot+sqrt(mass(k)**2+p5(1,k)**2+p5(2,k)**2+p5(3,k)**2)
    end do
    if (abs(ptot-srts)>1.e-05) then
       write(*,*) 'Problems in momenta_in_5BodyPS: energy conservation', ptot,srts
       stop
    end if

    if (debugFlag) then
       write(*,*) 'In momenta_in_5BodyPS, srts: ', srts
       do k=1,5
          write(*,*) k, mass(k), p5(1:3,k)
       end do
    end if

  end function momenta_in_5BodyPS


  !****************************************************************************
  !****f* nBodyPhaseSpace/momenta_in_6BodyPS
  ! NAME
  ! function momenta_in_6BodyPS (srts, mass) result (p6)
  ! PURPOSE
  ! This subroutine determines momenta of outgoing particles according to
  ! 6-body phase space. Vacuum prescription.
  ! INPUTS
  ! * real :: srts               -- sqrt(s)
  ! * real, dimension(6) :: mass -- masses of particles
  ! OUTPUT
  ! * real, dimension(3,6) :: p6 -- 6 three-vectors of momenta
  !   (first index: momentum index, second index: particle)
  ! NOTES
  ! Uses the recursion formula for the phase space (see PDG-2002, p.269)
  ! in a similar way as in subroutine momenta_in_4BodyPS.
  !****************************************************************************
  function momenta_in_6BodyPS (srts, mass) result (p6)

    use random, only: rn, rnOmega
    use lorentzTrafo, only: lorentz

    real, intent(in) :: srts
    real, intent(in), dimension(6) :: mass
    real, dimension(3,6) :: p6

    integer, parameter :: nimax=500  ! maximum number of iterations
    real :: q2max,q2min,s,pfinalm,q2,q,pfinal,ps5,ps5max,e12345,ptot
    integer :: ni,k
    logical :: flag
    real, dimension(6) :: mass2
    real, dimension(1:3) :: beta
    real, dimension(0:3) :: fourVector

    if (srts<=sum(mass)) then
       write(*,'(A,2F8.4)') 'problem in momenta_in_6BodyPS: srts too small', srts, sum(mass)
       return
    end if

    q2max=(srts-mass(6))**2
    q2min=(mass(1)+mass(2)+mass(3)+mass(4)+mass(5))**2
    s = srts**2
    mass2 = mass**2

    pfinalm=sqrt(max(0., (s+q2min-mass2(6))**2/(4.*s) - q2min ))

    ps5max = integrate_nBodyPS (sqrt(q2max), mass(1:5))

    flag=.true.
    ni=0
    if (q2max>q2min) then
       do while(flag)
          ni=ni+1

          q2=rn()*(q2max-q2min)+q2min
          pfinal=sqrt(max(0., (s+q2-mass2(6))**2/(4.*s) - q2 ))

          q=sqrt(q2)
          ps5 = integrate_nBodyPS (q, mass(1:5))

          if (pfinal*ps5 > rn()*pfinalm*ps5max) flag = .false.

          if (ni>nimax) then
             write(*,*) 'problems in momenta_in_6BodyPS ni.gt.nimax',srts, sum(mass)
             q2=0.5*(q2max-q2min)+q2min
             q=sqrt(q2)
             pfinal=sqrt(max(0., (s+q2-mass2(6))**2/(4.*s) - q2 ))
             flag=.false.
          end if
       end do
    else
       q2=q2min
       q=sqrt(q2)
       pfinal=0.
    end if

    !... select momentum of particle 6:
    p6(1:3,6)=pfinal*rnOmega()

    !... determine momenta of particles 1-5 in their rest frame:
    p6(1:3,1:5) = momenta_in_5BodyPS (q, mass(1:5))

    !... boost momenta to overall rest frame:
    e12345=sqrt(q2+pfinal**2)
    beta(1:3) = p6(1:3,6)/e12345

    do k=1,5
       fourVector(1:3)= p6(1:3,k)
       fourVector(0)  = sqrt(mass(k)**2+Dot_Product(fourVector(1:3),fourVector(1:3)))
       call lorentz(beta,fourVector)
       p6(1:3,k) = fourVector(1:3)
    end do

    !... checks:
    do k=1,3
       ptot = sum(p6(k,1:6))
       if (abs(ptot)>1.e-05) then
          write(*,*) 'Problems in momenta_in_6BodyPS: momentum conservation', ptot
          stop
       end if
    end do

    ptot=0.
    do k=1,6
       ptot=ptot+sqrt(mass(k)**2+p6(1,k)**2+p6(2,k)**2+p6(3,k)**2)
    end do
    if (abs(ptot-srts)>1.e-05) then
       write(*,*) 'Problems in momenta_in_6BodyPS: energy conservation', ptot,srts
       stop
    end if

  end function momenta_in_6BodyPS


  !****************************************************************************
  !****f* nBodyPhaseSpace/integrate_nBodyPS
  ! NAME
  ! function integrate_nBodyPS (srts, mass) result (phi)
  ! PURPOSE
  ! Determine the phase space volume for n=2,...,6 outgoing particles
  ! (see definition in PDG-booklet-2002 p. 269)
  ! INPUT:
  ! * real, intent(in) :: srts                ! c.m. energy (GeV),
  ! * real, intent(in), dimension(:) ::  mass ! masses of outgoing particles (GeV)
  ! * OUTPUT:
  ! * real, intent(out) :: phi                ! phase space volume (GeV**(2*n-4))
  ! NOTES:
  ! n=2,...,20 for the case when flagKopylov=.true.
  !****************************************************************************
  function integrate_nBodyPS (srts, mass) result (phi)

    use constants, only: pi
    use threeBodyPhaseSpace, only: Integrate_3BodyPS
    use twoBodyTools, only: pCM
    use um1, only: SN

    real, intent(in) :: srts
    real, intent(in), dimension(:) ::  mass
    real :: phi

    logical, parameter :: flagKopylov=.true. ! If .true. --- use analytical formulae
    ! of G.I. Kopylov, NPA 36 (1962) 425
    ! (communicated by Igor Pschenichnov)

    real :: ps, q2min, q2max, thresh
    integer :: n
    real, dimension(1:21) :: mass_wk

    n = size(mass,dim=1)

    thresh = sum(mass(:))

    if (srts<=thresh) then
       phi = 0.
       return
    end if

    srtsG = srts
    massG(1:n) = mass(1:n)

    if (n==2) then
       phi = pi*pcm(srts,mass(1),mass(2))/(2.*pi)**6/srts
       return
    end if

    if (flagKopylov .and. n>=2 .and. n<=20) then
       mass_wk(:)=0.
       mass_wk(1:n)=mass(1:n)
       phi=SN(srts,n,mass_wk)/(2.*pi)**(3*n)
       return
    end if

    select case (n)
    case (3)
       ps = Integrate_3BodyPS (srts, mass(1), mass(2), mass(3))
       phi = ps/(2.*pi)**7/16.
    case (4)
       q2min = (mass(3)+mass(4))**2
       q2max = (srts-mass(1)-mass(2))**2
       phi = sint(fun1,q2min,q2max,100)*pi/(2.*pi)**10/16.
    case (5)
       q2min = (mass(3)+mass(4)+mass(5))**2
       q2max = (srts-mass(1)-mass(2))**2
       phi = sint(fun2,q2min,q2max,50)/16.**2/(2.*pi)**11
    case (6)
       q2min = (mass(3)+mass(4)+mass(5)+mass(6))**2
       q2max = (srts-mass(1)-mass(2))**2
       phi = sint(fun3,q2min,q2max,50)
    case default
       write(*,*) 'In integrate_nBodyPS n=', n
       write(*,*) 'not implemented'
       stop
    end select

  end function integrate_nBodyPS



  real function fun1(q2)
    use threeBodyPhaseSpace, only: Integrate_3BodyPS
    use twoBodyTools, only: pCM
    real, intent(in) :: q2
    real :: q

    q = sqrt(q2)
    fun1 = pCM(q,massG(3),massG(4)) / q * Integrate_3BodyPS(srtsG,q,massG(1),massG(2))
  end function fun1




  real function fun2(q2)
    use threeBodyPhaseSpace, only: Integrate_3BodyPS
    real, intent(in) :: q2
    real :: q

    q = sqrt(q2)
    fun2 = Integrate_3BodyPS(q,massG(3),massG(4),massG(5)) * Integrate_3BodyPS(srtsG,q,massG(1),massG(2))
  end function fun2




  real function fun3(q2)
    use constants, only: pi
    use threeBodyPhaseSpace, only: Integrate_3BodyPS

    real, intent(in) :: q2
    real :: qp2min,qp2max,phi4

    qG = sqrt(q2)
    qp2min = (massG(4)+massG(5)+massG(6))**2
    qp2max = (qG-massG(3))**2

    phi4 = sint(fun4,qp2min,qp2max,100)*pi/qG/16./(2.*pi)**10
    fun3 = phi4/16./(2.*pi)**4 * Integrate_3BodyPS(srtsG,qG,massG(1),massG(2))
  end function fun3



  real function fun4(qp2)
    use threeBodyPhaseSpace, only: Integrate_3BodyPS
    use twoBodyTools, only: pCM

    real, intent(in) :: qp2
    real :: qp

    qp = sqrt(qp2)
    fun4 = pcm(qG,qp,massG(3)) * Integrate_3BodyPS(qp,massG(4),massG(5),massG(6))
  end function fun4





  REAL FUNCTION SINT(F,E0,EN,N)

    ! INTEGRATE REAL*4 FUNCTION F(E) OVER REAL*4 E,
    ! FROM E0 TO EN (METOD TRAPECIJ)
    ! N --- NUMBER OF BINS,

    real, intent(in) :: E0, EN   ! Integration limits
    integer, intent(in) :: N     ! Number of bins

    integer :: I
    real :: DE, E

    Interface
       real function F(x)
         real, intent(in) :: x
       end function F
    end Interface

    SINT = 0.
    DE = (EN-E0)/N

    do I=1,N-1
       E = E0 + DE*I
       SINT = SINT + F(E)
    end do

    SINT = SINT*DE + DE/2. * (F(E0)+F(EN))

  END FUNCTION SINT


end module nBodyPhaseSpace
