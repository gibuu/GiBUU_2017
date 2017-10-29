!******************************************************************************
!****m* /nucleusDefinition
! NAME
! module nucleusDefinition
!
! PURPOSE
! Define structure "tNucleus".
!
! NOTES
! Because of internal module dependencies we do not include any routines
! working with/on the defined objects in this module but rather give them
! in the module "nucleus".
!******************************************************************************
module nucleusDefinition

  implicit none
  private

  integer, parameter :: maxIndex_def= 2000
  real   , parameter ::       dx_def=0.01

  !****************************************************************************
  !****t* nucleusDefinition/tNucleus
  ! SOURCE
  !
  type, public :: tNucleus
     real      :: radius = 0.
     real      :: surface = 0.
     real      :: density = 0.
     real,dimension(1:3) :: position = 0.  ! position in calculation frame
     real,dimension(1:3) :: velocity = 0.  ! velocity in calculation frame
     integer   :: mass = 0                 ! in units of nucleon masses
     integer   :: charge = 0
     logical   :: fermiMotion=.true. ! switch to turn on/off FermiMotion

     logical   :: DoInit =.true.     ! if true, then first initialize !

     integer   :: densitySwitch_static=3
     real      :: fermiMomentum_input=0.225

     integer    :: MaxIndex = maxIndex_def
     real       :: dx       = dx_def
     real,dimension(0:maxIndex_def,2) :: densTab ! 1: p, 2: n

     real       :: MaxDist = 0.   ! Distance-CutOff: dens < 1e-6
     real       :: MaxDens(2)=0.  ! maximum of densTab

     real       :: chemPot = 0.   ! 'chemical' Potential for LDA (relative
                                  ! units of E0)
     logical    :: ReAdjustForConstBinding = .false.
     real       :: ConstBinding = 0.
     real       :: facN = 0., facP = 0. ! scaling factors in Readjust

  end type tNucleus
  !
  ! PURPOSE
  ! This type stores informations about nuclei, namely target and projectile.
  !
  ! NOTES
  ! Following parameters are valid for all density parametrisations:
  ! * mass
  ! * charge
  !
  ! parameters for Woods-Saxon:
  ! * radius
  ! * surface
  ! * density
  !****************************************************************************


  public :: WriteNucleus, WriteNucleusStaticDens
  public :: NucleusStaticDens, NucleusAverageDensity


contains


  !****************************************************************************
  !****s* nucleusDefinition/WriteNucleus
  ! NAME
  ! subroutine WriteNucleus(nuc)
  ! PURPOSE
  ! write the main parameters to stdout
  ! INPUTS
  ! * type(tNucleus), pointer :: Nuc -- nucleus to consider
  ! OUTPUT
  ! written to stdout
  !****************************************************************************
  subroutine WriteNucleus (nuc)
    type(tNucleus), pointer :: Nuc
    write(*,*) '::: Mass of nucleus    =',Nuc%mass
    write(*,*) '::: Charge of nucleus  =',Nuc%charge
    write(*,*) '::: Radius of nucleus  =',Nuc%radius
    write(*,*) '::: Surface of nucleus =',Nuc%surface
    write(*,*) '::: Central density    =',Nuc%density
  end subroutine WriteNucleus


  !****************************************************************************
  !****s* nucleusDefinition/WriteNucleusStaticDens
  ! NAME
  ! subroutine WriteNucleusStaticDens(name,nuc)
  ! PURPOSE
  ! write the static density to file
  !****************************************************************************
  subroutine WriteNucleusStaticDens (name, nuc)
    character*(*), intent(in) :: name
    type(tNucleus),  pointer :: Nuc

    integer :: i

    open(13,file=name,status='unknown')
    write(13,'(A)') '# r [fm], rho_B, rho_p, rho_n'

    do i = 0, min(nuc%MaxIndex,nint((nuc%radius+10*nuc%surface)/nuc%dx))
      write(13,'(F6.2,3ES12.4)') i*nuc%dx, nuc%densTab(i,1)+nuc%densTab(i,2), nuc%densTab(i,1), nuc%densTab(i,2)
    end do
    close(13)

  end subroutine WriteNucleusStaticDens


  !****************************************************************************
  !****f* nucleusDefinition/NucleusStaticDens
  ! NAME
  ! real function NucleusStaticDens(nuc,r,iType)
  ! PURPOSE
  ! return value of the tabulated nuclear density
  ! INPUTS
  ! * type(tNucleus) :: Nuc   -- nucleus to consider
  ! * real           :: r     -- radius [fm]
  ! * integer        :: iType -- 0: total, 1: proton, 2: neutron
  ! OUTPUT
  ! function value: density [units????]
  !****************************************************************************
  real function NucleusStaticDens (nuc, r, iType)
    type(tNucleus), pointer :: Nuc
    real,           intent(IN) :: r
    integer,        intent(IN) :: iType

    integer :: i

    i = nint(r/nuc%dx)

    select case (iType)
    case (0)
        NucleusStaticDens=nuc%densTab(i,1)+nuc%densTab(i,2)
    case (1,2)
       NucleusStaticDens=nuc%densTab(i,iType)
    case default
       write(*,*) 'ERROR: NucleusStaticDens; iType =', iType
       stop
    end select

  end function NucleusStaticDens


  !****************************************************************************
  !****s* nucleusDefinition/NucleusAverageDensity
  ! NAME
  ! subroutine printAverageDensity(nuc)
  ! PURPOSE
  ! Prints average density of nucleus
  !
  ! INPUTS
  ! * type(tNucleus) :: nuc
  !
  ! NOTES
  ! <rho> = int(rho   rho r**2 dr)/int(rho r**2 dr)
  !****************************************************************************
  subroutine NucleusAverageDensity (nuc)
    type(tNucleus),pointer :: nuc

    integer :: i
    real :: x, rho_average,rho_integral

    rho_average = 0.
    rho_integral= 0.

    do i=0,nuc%MaxIndex
       x = i*nuc%dx
       rho_average =rho_average +(nuc%densTab(i,1)+nuc%densTab(i,2))**2 *  x**2
       rho_integral=rho_integral+(nuc%densTab(i,1)+nuc%densTab(i,2))    *  x**2
    end do

    if (rho_integral/=0.) rho_average= rho_average/rho_integral

    write(*,'(A,F9.5)') ' ::: Average density of nucleus: ', rho_average

  end subroutine NucleusAverageDensity



end module nucleusDefinition
