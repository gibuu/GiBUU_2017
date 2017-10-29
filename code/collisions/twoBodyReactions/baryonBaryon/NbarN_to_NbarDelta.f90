!******************************************************************************
!****m* /NbarN_to_NbarDelta
! NAME
! module NbarN_to_NbarDelta
! PURPOSE
! Includes the Nbar N <-> Nbar Delta cross sections according to the one-pion
! exchange model.
! NOTES:
! This module has a similar structure as the module dimi.
!******************************************************************************
module NbarN_to_NbarDelta
  implicit none
  private


  !****************************************************************************
  !****n* NbarN_to_NbarDelta/initNbarN_to_NbarDelta
  ! NAME
  ! NAMELIST initNbarN_to_NbarDelta
  ! Includes the input variables:
  ! * delta_mass
  ! * maxPoints_mass
  ! * delta_srts
  ! * maxPoints_srts
  !****************************************************************************

  !****************************************************************************
  !****g* NbarN_to_NbarDelta/delta_mass
  ! SOURCE
  !
  real, save :: delta_mass=0.01
  !
  ! PURPOSE
  ! * grid step on a delta mass (GeV)
  !****************************************************************************

  !****************************************************************************
  !****g* NbarN_to_NbarDelta/maxPoints_mass
  ! SOURCE
  !
  integer, save :: maxPoints_mass=150
  !
  ! PURPOSE
  ! * number of the grid points on the delta mass
  !****************************************************************************

  !****************************************************************************
  !****g* NbarN_to_NbarDelta/delta_srts
  ! SOURCE
  !
  real, save :: delta_srts=0.01
  !
  ! PURPOSE
  ! * grid step on an invariant energy (GeV)
  !****************************************************************************

  !****************************************************************************
  !****g* NbarN_to_NbarDelta/maxPoints_srts
  ! SOURCE
  !
  integer, save :: maxPoints_srts=100
  !
  ! PURPOSE
  ! * number of the grid points on the invariant energy
  !****************************************************************************


  real, parameter     :: minimal_mass=1.076
  real, parameter     :: minimal_srts=2.014
  logical, save       :: initFlag=.true.


  public :: massNbarN_NbarDelta
  public :: NbarN_to_NbarDelta_Integrated
  public :: NbarDelta_to_NbarN
  public :: NbarN_to_NbarDelta_dm
  public :: calculate
  public :: calculate1
  public :: mNbarN_to_NbarD2


contains


  !****************************************************************************
  !****s* NbarN_to_NbarDelta/init
  ! NAME
  ! subroutine init
  ! PURPOSE
  ! Reads data out of namelist 'initNbarN_to_NbarDelta'
  !****************************************************************************
  subroutine init

    use output, only: Write_ReadingInput

    integer :: ios
    NAMELIST /initNbarN_to_NbarDelta/  delta_mass, maxPoints_mass, &
         & delta_srts, maxPoints_srts


    call Write_ReadingInput('initNbarN_to_NbarDelta',0)
    rewind(5)
    read(5,nml=initNbarN_to_NbarDelta,iostat=ios)
    call Write_ReadingInput('initNbarN_to_NbarDelta',0,ios)

    write(*,*) ' Set delta_mass to ', delta_mass,'.'
    write(*,*) ' Set maxPoints_mass to ', maxPoints_mass,'.'
    write(*,*) ' Set delta_srts to ', delta_srts,'.'
    write(*,*) ' Set maxPoints_srts to ', maxPoints_srts,'.'

    call Write_ReadingInput('initNbarN_to_NbarDelta',1)

  end subroutine init


  !****************************************************************************
  !****f* NbarN_to_NbarDelta/massNbarN_NbarDelta
  ! NAME
  ! function massNbarN_NbarDelta(srts) result(massDelta)
  ! PURPOSE
  ! Evaluates masses in Nbar N -> Nbar Delta process using NbarN_to_NbarDelta_dm
  ! This is done by a Monte-Carlo process utilizing dsigma/dm.
  ! NOTES
  !****************************************************************************
  function massNbarN_NbarDelta(srts) result(massDelta)

    use constants, only: mN
    use random, only: rn

    real, intent(in) :: srts
    real :: massDelta

    integer :: indexMass
    real      :: sum_dsdm
    real :: x
    real :: dsdm, dsdm_Integrated,sigma

    if (initFlag) then
       call init
       initFlag=.false.
    end if

    indexMass=0
    sum_dsdm=0.
    x=rn()
    do
       massDelta=minimal_Mass+float(indexMass)*delta_Mass
       call NbarN_to_NbarDelta_dm(dsdm,dsdm_Integrated,sigma,srts,massDelta)
       sum_dsdm=sum_dsdm+dsdm*delta_Mass
       if (sum_dsdm.ge.x*dsdm_Integrated) exit
       if (massDelta.gt.6.) then
          write(*,*) 'Problem in massNbarN_NbarDelta. massDelta > 6 GeV - Monte Carlo fails'
          write(*,*) 'mass of Delta=', massDelta,srts
          write(*,*) 'sqrt(s)=', srts
          write(*,*) sum_dsdm, dsdm_Integrated,x
          stop
       end if
       indexMass=indexMass+1
    end do

    if (indexMass==0) then
       massDelta= minimal_mass + rn()*(srts - minimal_mass - mN)
    else
       massDelta= massDelta + (rn()-0.5)*delta_Mass
       if (massDelta > srts-mN) then
          massDelta= minimal_mass + (float(indexMass)-0.5)*delta_Mass &
               + rn()*(srts-mN-minimal_mass &
               -(float(indexMass)-0.5)*delta_Mass)
          if (massDelta > srts-mN) then
             write(*,*) 'problems in massNbarN_NbarDelta:',srts,massDelta,indexMass
             stop
          end if
       end if
    end if

  end function massNbarN_NbarDelta


  !****************************************************************************
  !****s* NbarN_to_NbarDelta/NbarN_to_NbarDelta_Integrated
  ! NAME
  ! subroutine NbarN_to_NbarDelta_Integrated(dsdm_Integrated,srts)
  ! PURPOSE
  ! Returns the cross section for pbar n --> nbar delta^- scattering,
  ! integrated over the Delta mass.
  !****************************************************************************
  subroutine NbarN_to_NbarDelta_Integrated(dsdm_Integrated,srts)
    real, intent(in) ::  srts
    real, intent(out) :: dsdm_Integrated
    real, parameter :: dummy1=1.2
    real :: dummy2, dummy3
    call NbarN_to_NbarDelta_dm(dummy2,dsdm_Integrated, dummy3, srts, dummy1)
  end subroutine NbarN_to_NbarDelta_Integrated


  !****************************************************************************
  !****s* NbarN_to_NbarDelta/NbarDelta_to_NbarN
  ! NAME
  ! subroutine NbarDelta_to_NbarN(sigma,srts,mdel)
  ! PURPOSE
  ! Returns the cross section for nbar delta^- --> pbar n scattering.
  !****************************************************************************
  subroutine NbarDelta_to_NbarN(sigma,srts,mdel)
    real, intent(in) ::  srts, mdel
    real, intent(out) :: sigma
    real :: dummy1, dummy2
    call NbarN_to_NbarDelta_dm(dummy1,dummy2,sigma,srts,mdel)
  end subroutine NbarDelta_to_NbarN


  !****************************************************************************
  !****s* NbarN_to_NbarDelta/NbarN_to_NbarDelta_dm
  ! NAME
  ! subroutine NbarN_to_NbarDelta_dm(dsdm,dsdm_Integrated,sigma,srts,mdel)
  ! PURPOSE
  ! this subroutine calculates the cross sections for pbar n <--> nbar delta^-
  ! Stores cross sections at first call to a field and then it only returns
  ! those field values.
  ! INPUTS
  ! * real :: srts --- sqrt(s) in the process
  ! * real :: mdel --- mass of delta
  ! RESULT
  ! * real :: sigma --- total cross section (for reaction with incoming delta)
  ! * real :: dsdm --- mass differential cross section as function of mdel and
  !   srts (for reaction with outgoing delta)
  ! * real :: dsdm_Integrated ---mass differential cross section integrated
  !   over mass as function of srts (for reaction with outgoing delta)
  ! NOTES
  ! Actual cross section is computed by subroutine calculate below.
  !****************************************************************************
  subroutine NbarN_to_NbarDelta_dm(dsdm,dsdm_Integrated,sigma,srts,mdel)

    use output, only: Write_InitStatus

    real, intent(out) :: dsdm
    real, intent(out) :: dsdm_Integrated
    real, intent(out) :: sigma
    real, intent(in) :: srts
    real, intent(in) :: mdel

    real, allocatable, dimension(:,:), save  :: dsdm_save
    real, allocatable, dimension(:),   save  :: dsdm_integrated_save
    real, allocatable, dimension(:,:), save  :: sigma_save

    integer :: srts_Index, mass_Index
    real :: mass, wurzelS, delta_mass_new

    dsdm=0.
    dsdm_Integrated=0.
    sigma=0.

    if (srts.lt.minimal_srts .or. mdel.lt.minimal_mass) return

    if (initFlag) then

       call init

       call Write_InitStatus('Nbar N <--> Nbar Delta Xsections',0)

       allocate(dsdm_save(0:maxPoints_srts,0:maxPoints_mass))
       allocate(dsdm_integrated_save(0:maxPoints_srts))
       allocate(sigma_save(0:maxPoints_srts,0:maxPoints_mass))

       do srts_Index=0,maxPoints_srts
          wurzelS  =minimal_srts   +  Real(srts_Index)   * delta_srts
          do mass_Index=0,maxPoints_mass
             mass    =  minimal_Mass +  Real(mass_Index)* delta_mass
             call calculate(dsdm,sigma,wurzelS,mass)
             dsdm_save(srts_Index, mass_Index)=dsdm
             sigma_save(srts_Index, mass_Index)=sigma
          end do
          ! Integrate dsDm_save over the mass index
          dsdm_integrated_save(srts_Index) = Sum(dsdm_save(srts_Index, : )) * delta_mass
       end do

       initFlag=.false.
       call Write_InitStatus('Nbar N <--> Nbar Delta Xsections',1)

    end if

    srts_Index=NINT((srts-minimal_srts)/delta_srts)
    mass_Index=NINT((mdel-minimal_mass)/delta_mass)

    if ((srts_Index.gt.maxPoints_srts).or.(mass_Index.gt.maxPoints_mass)) then
       ! no fast possibility to get cross section
       ! give warning to tell user that system slows down
       write(*,*) 'Warning ! NbarN_to_NbarDelta_dm is out of bounds.'
       write(*,*) 'Need to calculate it. System therefore slow. If this happens too often, then change bounds.'
       write(*,*) 'SquareRoot(s)= ',srts, srts_Index, maxPoints_srts
       write(*,*) 'Mass of delta= ',mdel, mass_Index, maxPoints_mass
       dsdm_integrated=0.

       if (mass_Index.gt.maxPoints_mass) then
          ! Increase mass step size
          delta_mass_new=delta_mass*10.
          mass_Index=Max(   0,   NINT((mdel-minimal_mass)/delta_mass_new)  )
          if (mass_Index.gt.maxPoints_mass) then
             write(*,*) 'Critical error.stop!'
          end if
          write(*,*) 'New indizes for mass= ',mass_Index, maxPoints_mass
       else
          delta_mass_new=delta_Mass
       end if

       do mass_Index=1,maxPoints_mass
          mass    =  minimal_Mass +  Real(mass_Index)* delta_mass_new
          call calculate(dsdm,sigma,srts,mass)
          dsdm_integrated = dsdm_integrated+dsdm * delta_mass_new
       end do

       call calculate(dsdm,sigma,srts,mdel)
    else
       dsdm=dsdm_save(srts_Index, mass_Index)
       dsdm_integrated=dsdm_integrated_save(srts_Index)
       sigma=sigma_save(srts_Index, mass_Index)
    end if

  end subroutine NbarN_to_NbarDelta_dm


  !****************************************************************************
  !****s* NbarN_to_NbarDelta/calculate
  ! NAME
  ! subroutine calculate(dsdm,sigma,srts,mdel)
  ! PURPOSE
  ! Calculate the cross sections for pbar n <--> nbar delta^-
  ! INPUTS
  ! * real :: srts --- sqrt(s) in the process
  ! * real :: mdel ---  mass of delta
  ! RESULT
  ! * real :: sigma --- total cross section (for reaction with incoming delta)
  ! * real :: dsdm --- mass differential cross section as function of mdel and
  !   srts (for reaction with outgoing delta)
  !****************************************************************************
  subroutine calculate(dsdm,sigma,srts,mdel)

    use constants, only: pi, mN, mPi
    use particleProperties, only: hadron
    use idTable, only: delta
    use baryonWidth, only: fullWidthBaryon

    real,intent(out) ::  dsdm,sigma
    real,intent(in) :: mdel,srts

    real ::t,qi2,qf2,g,f,dcost,cost!,q,q0
    real :: E1,E2,E4,melem2,dsdodm,dsdo,fac,dmass!,E3
    integer :: nt,i
    logical, parameter :: debug = .false.

    dmass = hadron(Delta)%mass

    if (srts-2.*mn-mpi.lt.0..or.mdel.gt.srts-mn) then
       dsdm = 0.0
       sigma = 0.0
       return
    end if

    qi2=srts**2/4.-mn**2
    qf2=max((srts**2-(mdel+mn)**2)*(srts**2-(mdel-mn)**2)/(4.*srts**2),0.)

    E1 = sqrt(mn**2 + qi2)
    E2 = E1
    !       E3 = sqrt(mn**2+qf2)
    E4 = sqrt(mdel**2+qf2)

    if (mdel.lt.mn+mpi) then
       f=0.
    else
       g=FullWidthBaryon(delta,mdel)
       !        q=sqrt(max((mdel**2+mpi**2-mn**2)**2/(4.*mdel**2)-mpi**2,0.))
       !        q0=sqrt(max((dmass**2+mpi**2-mn**2)**2/(4.*dmass**2)-mpi**2,0.))
       !        g=0.118*(q/q0)**3*(0.2**2+q0**2)/(0.2**2+q**2)
       f=1./pi*mdel*g/((mdel**2-dmass**2)**2+(mdel*g)**2)
    end if

    if (qf2.gt.0.) then
       dcost = min(1.,0.1*mpi**2/(2.*sqrt(qi2*qf2)))
       nt = nint(1./dcost)
       dcost = 1./Real(nt)
    else
       nt = 1
       dcost = 1.
    end if

    if (debug) then
       write(*,*) 'In NbarN_to_NbarDelta/calculate srts, mdel, dcost:'
       write(*,*) srts, mdel, dcost
       write(*,*) 'qi2, qf2:', qi2, qf2
    end if

    dsdm=0.
    sigma=0.
    do i=-nt,nt

       cost=dcost*Real(i)

       t= mdel**2 + mn**2 - 2.*E4*E2 + 2.*cost*sqrt(qi2*qf2)

       !****   Relativistic matrix element^2 :
       melem2 = mNbarN_to_NbarD2(t,mn,mn,mn,mdel)

       !****   pbar n -> nbar delta^- :
       dsdodm= 1./(32.*pi*srts**2)*sqrt(qf2)*melem2*2.*mdel*f*0.389

       !****   nbar delta^- -> pbar n, factor 1./2. bec. of spin:
       dsdo= 1./(32.*pi*srts**2)*sqrt(qi2)*melem2*0.389/2.

       fac=1.
       if (abs(i).eq.nt) fac=0.5
       dsdm=dsdm+dsdodm/Real(nt)*fac
       sigma=sigma+dsdo/Real(nt)*fac

    end do

  end subroutine calculate


  !****************************************************************************
  !****s* NbarN_to_NbarDelta/calculate1
  ! NAME
  ! subroutine calculate1(s,t,dsigdt,mass_min,mass_max)
  ! PURPOSE
  ! Compute d sigma(pbar_1 n_2 -> nbar_3 delta_4^-) / d t (mb/GeV^2)
  ! integrated over mass of delta.
  ! INPUTS
  ! * real :: s --- s=(p1+p2)^2, Mandelstam variable (GeV^2)
  ! * real :: t --- t=(p1-p3)^2=(p4-p2)^2, Mandelstam variable (GeV^2)
  ! RESULT
  ! * real :: dsigdt ---  d sigma(pbar_1 n_2 -> nbar_3 delta_4^-)
  ! * real :: mass_min --- delta mass integration limits (GeV)
  ! * real :: mass_max --- delta mass integration limits (GeV)
  !****************************************************************************
  subroutine calculate1(s,t,dsigdt,mass_min,mass_max)

    use constants, only: pi, mN, mPi
    use particleProperties, only: hadron
    use idTable, only: delta
    use baryonWidth, only: fullWidthBaryon

    real, intent(in)  :: s,t
    real, intent(out) :: dsigdt,mass_min,mass_max

    real :: mdel,qi2,qi,qf2,qf,g,f,t0,A,B,C,fac!,q,q0
    real :: E1,melem2,dm,dmass!,E2,E4
    integer :: nmass, i
    logical, parameter :: debug = .false.

    dmass = hadron(Delta)%mass

    if (s-(2.*mn+mpi)**2.lt.0.) then
       dsigdt = 0.0
       return
    end if

    qi2=s/4.-mn**2
    qi=sqrt(max(qi2,0.))
    E1 = sqrt(mn**2+qi2)
    !         E2 = E1

    t0 = 2.*mn**2 - 2.*E1*mn
    A = 2.*mn**2 - t
    B = A*qi/(2.*mn**2)
    C = (A**2 - 4.*E1**2*mn**2 + A**2*qi2/mn**2)/(4.*mn**2)
    if (C.lt.0.) then
       write(*,*) ' In NbarN_to_NbarDelta/calculate1: C.lt.0., this must not happen !!!'
       stop
    end if
    if (t.ge.t0) then
       qf = B - sqrt(C)
    else
       qf = -B + sqrt(C)
    end if

    mass_max = (sqrt(s)-sqrt(mn**2+qf**2))**2 - qf**2
    if (mass_max.lt.0.) then
       write(*,*) ' In NbarN_to_NbarDelta/calculate1: mass_max.lt.0., this must not happen !!!'
       stop
    end if
    mass_max = sqrt(mass_max)
    mass_min = mn + mpi

    dm=0.001
    nmass=max(1,nint((mass_max-mass_min)/dm))
    dm = (mass_max-mass_min)/Real(nmass)

    if (debug) then
       write(*,*) 'In NbarN_to_NbarDelta/calculate1 srts, t, mass_min, mass_max, dm:'
       write(*,*) sqrt(s), t, mass_min, mass_max, dm
    end if

    dsigdt = 0.

    do i = 0,nmass

       mdel = mass_min + dm*Real(i)

       qf2 = (s + mdel**2 - mn**2)**2/(4.*s) - mdel**2

       !            E4 = sqrt(mdel**2+qf2)

       if (mdel.lt.mn+mpi) then
          f=0.
       else
          g=FullWidthBaryon(delta,mdel)
          !              q=sqrt(max((mdel**2+mpi**2-mn**2)**2/(4.*mdel**2)-mpi**2,0.))
          !              q0=sqrt(max((dmass**2+mpi**2-mn**2)**2/(4.*dmass**2)-mpi**2,0.))
          !              g=0.118*(q/q0)**3*(0.2**2+q0**2)/(0.2**2+q**2)
          f=1./pi*mdel*g/((mdel**2-dmass**2)**2+(mdel*g)**2)
       end if

       !****      Relativistic matrix element^2 :
       melem2 = mNbarN_to_NbarD2(t,mn,mn,mn,mdel)

       fac=1.
       if (i.eq.0 .or. i.eq.nmass) fac=0.5
       dsigdt = dsigdt + melem2/(64.*pi*s)*2.*mdel*f*0.389*dm*fac

    end do

    dsigdt = dsigdt/qi2

  end subroutine calculate1


  !****************************************************************************
  !****s* NbarN_to_NbarDelta/mNbarN_to_NbarD2
  ! NAME
  ! real function mNbarN_to_NbarD2 (t, m_1, m_2, m_3, m_4)
  ! PURPOSE
  ! Matrix element squared and averaged over initial and summed over
  ! final spins for pbar_1 n_2 -> nbar_3 delta_4^-
  ! within OPEM with effective nucleon and delta masses.
  ! INPUTS
  ! * real :: t --- t=(p1-p3)^2, Mandelstam variable (GeV^2)
  ! * real :: m_1,m_2,m_3,m_4  -- effective masses of particles (GeV)
  !****************************************************************************
  real function mNbarN_to_NbarD2 (t, m_1, m_2, m_3, m_4)

    use constants, only: mPi
    use particleProperties, only: hadron
    use idTable, only: Delta

    real, intent(in) :: t,m_1,m_2,m_3,m_4
    real :: dmass,qr2t,q2t,zt,ft
    ! Coupling constants from V. Dmitriev et al, NPA 459, 503 (1986):
    real, parameter :: fps=2.202d0,fp=1.008d0
    real, parameter :: lambda= 0.50d0  ! readjusted by A.L.
    real, parameter :: kappa =0.2d0

    logical, parameter :: debugflag = .false.

    dmass = hadron(Delta)%mass

    ft=(lambda**2-mpi**2)/(lambda**2-t)       ! From V. Dmitriev et al
    !ft=0.72/(1.+(-t+mpi**2)/(4.73*mpi**2)) + 0.28   ! From R. Armenteros & B. French at low energies

    qr2t=(dmass**2+m_2**2-t)**2/(4.d0*dmass**2)-m_2**2
    q2t=(m_4**2+m_2**2-t)**2/(4.d0*m_4**2)-m_2**2
    zt=(qr2t+kappa**2)/(q2t+kappa**2)

    if (debugflag) write(*,*) 'ft,qr2t,q2t,zt:', ft,qr2t,q2t,zt


    mNbarN_to_NbarD2 = 2. * (fp*fps*ft**2)**2*zt/mpi**4 / (t-mpi**2)**2 &
         &* (m_1+m_3)**2 * (t-(m_1-m_3)**2) &
         &* (t-(m_4-m_2)**2) * ((m_4+m_2)**2-t)**2/6./m_4**2


    if (debugflag) write(*,*) ' mNbarN_to_NbarD2:', mNbarN_to_NbarD2

  end function mNbarN_to_NbarD2


end module NbarN_to_NbarDelta
