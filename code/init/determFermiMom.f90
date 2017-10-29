!******************************************************************************
!****m* /determine_fmbybE
! NAME
! module determine_fmbybE
!
! PURPOSE
! determine the fermi momentum
!
! "fmbybE" stands for "fermi momentum by binding energy"
!
! NOTES
! [KG]: I cut this out of "initNucPhaseSpase.f90" and brought it to the
! code format we all agreed on.
! I do not guarantee full functionality and meaningful documentation.
!******************************************************************************

module determine_fmbybE
  use particleDefinition, only: particle
  implicit none

  private

  public :: determine_fermimomentum,determine_fermiNucDLDA


  !****************************************************************************
  !****g* determine_fmbybE/teilchen_fermi
  ! SOURCE
  !
  type(particle), save, public :: teilchen_fermi
  ! PURPOSE
  ! ???
  !****************************************************************************

  real, save :: baryondensity

contains

  !****************************************************************************
  !****s* determine_fmbybE/determine_fermimomentum
  ! NAME
  ! subroutine determine_fermimomentum(baryondens,fermimomentum)
  ! PURPOSE
  ! ...
  ! INPUTS
  ! * real :: baryondens -- ???
  ! OUTPUT
  ! * real :: fermimomentum -- ???
  !
  ! NOTES
  ! * one should better define this as function. (?)
  ! * one should add documentation. (!)
  !****************************************************************************
  subroutine determine_fermimomentum(baryondens,fermimomentum)

    use constants, only: mN
    use particleDefinition
    use cern_lib, only: dzerox

    real, intent(out) :: fermimomentum
    real, intent(in) :: baryondens

    logical, save::initflag=.true.
    if (initflag) then

       call setToDefault(teilchen_fermi)
       teilchen_fermi%mass=mN
       ! finalstate%charge =charge_out
       ! teilchen_fermi%momentum=p_out
       teilchen_fermi%position=(/0.,0.,0./)
       teilchen_fermi%antiparticle=.false.
       teilchen_fermi%id=1
       teilchen_fermi%perturbative=.false.
       ! finalState%productionTime=0.
       ! finalState%lastCollisionTime=0.
       ! finalState%formationTime=0.
       ! finalState%scaleCS=1.
       ! finalState%in_Formation=.false.
       write(*,*) 'determin_fermimomentum initialized'
       initflag=.false.
    end if

    baryondensity=baryonDens
    fermimomentum=DZEROX(0d0,1d0,1d-9,1000,func,1)

    if (fermimomentum.lt. 0.01) then
       write(*,*) 'fermimomentum=', fermimomentum, ' rho=', baryondens
    end if

  end subroutine determine_fermimomentum

  !****************************************************************************
  !****s* determine_fmbybE/determine_fermiNucDLDA
  ! NAME
  ! subroutine determine_fermiNucDLDA(baryonplace,rhocent,baryonpFermi)
  ! PURPOSE
  ! ...
  ! INPUTS
  ! * real, dimension(1:3) :: baryonplace -- ???
  ! * real                 :: rhocent -- ???
  ! OUTPUT
  ! * real :: baryonpFermi -- ???
  !
  ! NOTES
  ! * one should better define this as function. (?)
  ! * one should add documentation. (!)
  !****************************************************************************
  subroutine determine_fermiNucDLDA(baryonplace,rhocent,baryonpFermi)
    use constants, only: hbarc, pi, mN
    use NucDLDA, only: getEParticleLaplace
    use densitymodule, only: gridSpacing
    use baryonPotentialModule, only: rhoLaplace

    real, intent(out):: baryonpFermi
    real, intent(in) :: rhocent
    real, dimension(1:3),intent(in) :: baryonplace
    real, dimension(1:4):: Etemp
    real :: rhoptemp, pFermi2temp, difftemp
    integer, save :: notBound=0
    logical, save :: firsttime=.true.

    rhoptemp=rhoLaplace(baryonPlace,gridSpacing)
    pFermi2temp=(1.5*pi**2*rhocent)**(2./3.)*hbarc**2
    call getEParticleLaplace(Etemp,rhocent,rhoptemp,pFermi2temp)
    baryonpFermi=0.
    if (Etemp(1).ge.0.) then
       difftemp=pFermi2temp-(2.*mN)/(1000.)*Etemp(1)
       if (difftemp.ge.0.) then
          baryonpFermi=SQRT(difftemp)
          if (firsttime) then
             open(97,file='NotBound.dat')
             firsttime=.false.
          else
             open(97,file='NotBound.dat',position='Append')
          end if
          notBound=notBound+1
          write(97,fmt='(I5,1X,F10.4,1X,F10.4,1X,F10.4)') notBound,pFermi2temp/(2.*mN)*1000. &
               &, difftemp/(2.*mN)*1000.,Etemp(1)
          close(97)
       else
          baryonpFermi=0.
       end if
    else
       baryonpFermi=SQRT(pFermi2temp)
    end if

  end subroutine determine_fermiNucDLDA


  !****************************************************************************
  !****if* determine_fmbybE/func
  ! NAME
  ! double precision function func
  ! PURPOSE
  ! function used to find a zero in determine_fermimomentum
  !****************************************************************************
  double precision function func(p)
    use baryonPotentialModule, only: BaryonPotential
    use mediumDefinition
    real, intent(in):: p
    type(medium)    :: med
    teilchen_fermi%momentum(1:3)=(/p,0.,0./)
    med%temperature    =0.
    med%useMedium      =.true.
    med%density        = baryondensity
    med%densityProton  = baryondensity/2
    med%densityNeutron = baryondensity/2
    func=sqrt((teilchen_fermi%mass+BaryonPotential(teilchen_fermi,med,.true.))**2 +p**2)-teilchen_fermi%mass+0.016
    !non relativistic energy  -Binding energy
  end function func

end module determine_fmbybE
