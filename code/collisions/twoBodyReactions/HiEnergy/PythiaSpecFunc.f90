!******************************************************************************
!****m* /PythiaSpecFunc
! NAME
! module PythiaSpecFunc
! PURPOSE
! Proper spectral functions for broad resonances in Pythia & Fritiof.
! We currently treat:
!  * vector mesons (rho, omega, phi)
!  * Delta(1232)
!******************************************************************************
module PythiaSpecFunc

  use mediumDefinition

  implicit none

  private

!   !*************************************************************************
!   !****g* PythiaSpecFunc/srtFreeVMMass
!   ! SOURCE
!   !
!   real, save :: srtFreeVMMass
!   ! NOTES
!   ! Needed by VM_Mass, so one needs to set this parameter before
!   ! PYTHIA/FRITIOF is called.
!   !*************************************************************************


  !****************************************************************************
  !****g* PythiaSpecFunc/mediumAtCollision
  ! SOURCE
  !
  type(medium), save :: mediumAtColl
  ! NOTES
  ! Needed by VM_Mass, so one needs to set this parameter before
  ! PYTHIA/FRITIOF is called.
  !****************************************************************************

  public :: Init_VM_Mass
  public :: VM_Mass
  public :: Delta_Mass
  public :: Delta_Mass_Full

  real, dimension(0:3), parameter :: momLRF = (/0.,0.,0.,0./)

contains

  !****************************************************************************
  !****s* PythiaSpecFunc/Init_VM_Mass
  ! NAME
  ! subroutine Init_VM_Mass(srts,position)
  ! PURPOSE
  ! Initialize the parameters srtFreeVMMass and mediumAtCollision before
  ! calling VM_Mass.
  ! INPUTS
  ! * real          :: srts          -- sqrt(s) of collision
  ! * real,OPTIONAL :: position(1:3) -- place of collision
  !****************************************************************************
  subroutine Init_VM_Mass(srts,position)
    use mediumModule, only: mediumAt
    real, intent(in) :: srts
    real, intent(in), optional :: position(1:3)

!     srtFreeVMMass = srts

    if (present(position)) then
       mediumAtColl = mediumAt(position)
    else
       mediumAtColl = mediumAt( (/1000.0, 0.0, 0.0/) ) ! somewhere far off
    end if

  end subroutine Init_VM_Mass

  !****************************************************************************
  !****f* PythiaSpecFunc/VM_Mass
  ! NAME
  ! real function VM_Mass(kfa,nReject)
  ! PURPOSE
  ! Select mass of vector mesons for Pythia/Fritiof (function pymass/ulmass)
  ! according to a relativistic Breit-Wigner with mass-dependent width.
  !
  ! Used for rho, omega and phi mesons.
  ! INPUTS
  ! * integer :: kfa  -- KF code of the particle
  ! OUTPUT
  ! * function value returns the mass
  ! * integer, OPTIONAL :: nReject -- number of steps in rejection method
  ! NOTES
  ! * Init_VM_Mass needs to be called before Pythia/Fritriof is called in
  !   order to set the medium at the position
  !   (THIS IS IMPORTANT!)
  ! * The mass selection is done in principle like in massass or assmass,
  !   but here we have improved the MC selection a little bit by dividing
  !   the mass region into slices.
  !****************************************************************************
  real function VM_Mass(kfa,nReject)
    use MassAssInfoDefinition, only: Get_UseMassAssInfo
    use AssignMassMC, only: AssignMass_1
    use IDTable, only: rho,omegaMeson,phi
    use CALLSTACK, only: TRACEBACK
    implicit none

    integer,intent(in):: kfa
    integer,intent(out),OPTIONAL :: nReject

    integer :: ID, nRej
    real :: maxMass

    if (Get_UseMassAssInfo()) then
       ! decode KF code into GiBUU ID
       select case (kfa)
       case (113,213)
          ID=rho
       case (223)
          ID=omegaMeson
       case (333)
          ID=phi
       case default
          write(*,*) "KFA = ",kfa
          call TRACEBACK()
       end select

!       maxMass = srtFreeVMMass
       maxMass = 3.0
       call AssignMass_1(ID,maxMass,momLRF,mediumAtColl,VM_Mass,nRej)

    else
       VM_Mass = VM_Mass_Full(kfa,nRej)
    end if

    if (present(nReject)) nReject = nRej

  end function VM_Mass

  !****************************************************************************
  !****f* PythiaSpecFunc/Delta_Mass
  ! NAME
  ! real function Delta_Mass(nReject)
  ! PURPOSE
  ! Select mass of vector mesons for Pythia/Fritiof (function pymass/ulmass)
  ! according to a relativistic Breit-Wigner with mass-dependent width.
  !
  ! Used for Delta Resonances.
  ! OUTPUT
  ! * function value returns the mass
  ! * integer, OPTIONAL :: nReject -- number of steps in rejection method
  ! NOTES
  ! * This is important for dilepton simulations,
  !   where high-mass Deltas give a large contribution to the dilepton spectrum.
  ! * Treatment as in VM_Mass, but with slightly changed mass slices.
  !   (weight bwd is not always an increasing function with M!)
  !****************************************************************************
  real function Delta_Mass(nReject)
   use MassAssInfoDefinition, only: Get_UseMassAssInfo
    use AssignMassMC, only: AssignMass_1
    implicit none

    integer,intent(out),OPTIONAL :: nReject

    integer :: nRej
    real :: maxMass

    if (Get_UseMassAssInfo()) then

!       maxMass = srtFreeVMMass
       maxMass = 3.0
       call AssignMass_1(2,maxMass,momLRF,mediumAtColl,Delta_Mass,nRej)

    else
       Delta_Mass = Delta_Mass_Full(nRej)
    end if

    if (present(nReject)) nReject = nRej
  end function Delta_Mass

  !****************************************************************************
  !****f* PythiaSpecFunc/VM_Mass_Full
  ! NAME
  ! real function VM_Mass_Full(kfa,nReject)
  ! PURPOSE
  ! Select mass of vector mesons for Pythia/Fritiof (function pymass/ulmass)
  ! according to a relativistic Breit-Wigner with mass-dependent width.
  !
  ! Used for rho, omega and phi mesons.
  !
  ! This is the internal routine, which is called, if UseMassAssInfo is set false
  ! INPUTS
  ! * integer :: kfa  -- KF code of the particle
  ! OUTPUT
  ! * function value returns the mass
  ! * integer, OPTIONAL :: nReject -- number of steps in rejection method
  ! NOTES
  ! * Init_VM_Mass needs to be called before Pythia/Fritriof is called in
  !   order to set the medium at the position
  !   (THIS IS IMPORTANT!)
  ! * The mass selection is done in principle like in massass or assmass,
  !   but here we have improved the MC selection a little bit by dividing
  !   the mass region into slices.
  !****************************************************************************
  real function VM_Mass_Full(kfa,nReject)
    use particleProperties, only: hadron
    use IDTable, only: rho,omegaMeson,phi
    use random
    use mesonWidthMedium, only: WidthMesonMedium,get_MediumSwitchMesons
    use mesonPotentialModule, only: vecMes_massShift
    use output
    use monteCarlo, only: MonteCarloChoose
    use CALLSTACK, only: TRACEBACK

    implicit none

    integer,intent(in):: kfa
    integer,intent(out),OPTIONAL :: nReject


    real:: m,minmass,maxmass,bwd,m2,ymin,ymax,y,mres0,gamres0,gamtot,spectral,intfac,maxbwd
    real:: gamtot_max,spectral_max,intfac_max
    real:: mres02,gamres024
    integer:: ncount,vm
    integer, parameter:: ncountm = 100
    logical:: lcount

    integer, parameter :: nB = 7
    real, dimension(nB) :: mbound,bweight,bmaxbwd, mbound2
    integer :: i,iB
    real :: bweightTot

    ! decode KF code into GiBUU ID
    select case (kfa)
    case (113,213)
       vm=rho
    case (223)
       vm=omegaMeson
    case (333)
       vm=phi
    case default
       write(*,*) "KFA = ",kfa
       call TRACEBACK()
    end select

    if (present(nReject)) nReject=0

    if (get_MediumSwitchMesons() .and. mediumAtColl%useMedium) then
       minMass = - vecMes_massShift(vm,mediumAtColl%density) + 0.01
    else
       minMass = hadron(vm)%minmass
    end if

    ! use a constant maxmass, to avoid changing the relative strength of different channels
    maxmass = 3.0

    mres0    = hadron(vm)%mass
    gamres0  = WidthMesonMedium(vm,mres0,momLRF,mediumAtColl)

    mres02   = mres0**2
    gamres024= gamres0**2/4

    ! Divide the mass region into slices, calculate values at the borders:

    mbound = (/minMass,mres0-2*gamres0,mres0-gamres0,mres0,&
         & mres0+gamres0,mres0+2*gamres0,maxMass/)

    do i=1,nB
       if (mbound(i).lt.minMass) mbound(i)=minMass
       if (mbound(i).gt.maxMass) mbound(i)=maxMass
    end do

    mbound2 = mbound**2

    do i=1,nB
       bweight(i) = atan2(2*(mbound(i)-mres0),gamres0)

       gamtot_max   = WidthMesonMedium(vm,mbound(i),momLRF,mediumAtColl)
       spectral_max = mbound2(i)*gamtot_max*gamres0/((mres02-mbound2(i))**2 + gamtot_max**2*mbound2(i))
       intfac_max   = gamres024/((mbound(i)-mres0)**2 + gamres024)
       bmaxbwd(i)   = spectral_max/intfac_max * 1.05
    end do

    ! Calculate the integral of the Cauchy distribution and their
    ! modification:

    do i=1,nB-1
       bweight(i)=(bweight(i+1)-bweight(i))
       bmaxbwd(i) = max(bmaxbwd(i),bmaxbwd(i+1))
       bweight(i) = bweight(i) * bmaxbwd(i) ! modify the Cauchy weights
    end do

    ! Now do the actual MC:

    ncount=0
    lcount=.true.
    do while(lcount)
       ncount=ncount+1

       ! STEP 1: Select the mass region:
       ! (it is important to have this step in the overall loop, since
       ! this cures inaccuracies in the calculation of the weights)

       iB=MonteCarloChoose(bweight(1:nB-1), bweightTot)
       ymax=2.*atan2(2*(mbound(iB+1)-mres0),gamres0)
       ymin=2.*atan2(2*(mbound(iB)-mres0),gamres0)
       maxbwd = bmaxbwd(iB)

       ! STEP 2: generate random value according Cauchy with constant width

       y=ymin+rn()*(ymax-ymin)
       m=.5*tan(y/2.)*gamres0+mres0
       m=min(max(m,minmass),maxmass)
       m2=m**2

       ! STEP 3: Do the rejection

       gamtot=WidthMesonMedium(vm,m,momLRF,mediumAtColl)
       spectral=m2*gamtot*gamres0/((mres02-m2)**2 + gamtot**2*m2)
       intfac=gamres024/((m-mres0)**2 + gamres024)
       bwd=spectral/intfac/maxbwd

       if (bwd>1) then
          if (DoPr(1)) then
             write(*,'(A,i3,A,2F12.3,A,2G12.3)') 'Problems in VM_Mass: (bwd>1)  vm=',vm, &
                  & ' mass=',mbound(iB),mbound(iB+1),' bwd:',bwd,maxbwd
             write(*,'(A,1P,10g12.3)')'If this problem occurs frequently, adjust slicing !!!'
             write(*,*)
          end if
       end if

       if (rn()<=bwd) lcount =.false.

       if (ncount>=ncountm) then
          if (DoPr(1)) then
             write(*,'(A,i3,A,2G12.3)') 'Problems in VM_Mass: (ncount)  vm=',vm, &
                  & ' mass=',mbound(iB),mbound(iB+1)
          end if
          m=mres0
          lcount=.false.
       end if
    end do

    if (present(nReject)) nReject=ncount

    VM_Mass_Full = m

  end function VM_Mass_Full


  !****************************************************************************
  !****f* PythiaSpecFunc/Delta_Mass_Full
  ! NAME
  ! real function Delta_Mass_Full(nReject)
  ! PURPOSE
  ! Select mass of vector mesons for Pythia/Fritiof (function pymass/ulmass)
  ! according to a relativistic Breit-Wigner with mass-dependent width.
  !
  ! Used for Delta Resonances.
  !
  ! This is the internal routine, which is called, if UseMassAssInfo is set false
  ! OUTPUT
  ! * function value returns the mass
  ! * integer, OPTIONAL :: nReject -- number of steps in rejection method
  ! NOTES
  ! * This is important for dilepton simulations,
  !   where high-mass Deltas give a large contribution to the dilepton spectrum.
  ! * Treatment as in VM_Mass, but with slightly changed mass slices.
  !   (weight bwd is not always an increasing function with M!)
  !****************************************************************************
  real function Delta_Mass_Full(nReject)
    use particleProperties, only: hadron
    use IDTable, only: Delta
    use random
    use baryonWidth, only: FullWidthBaryon
    use output
    use monteCarlo, only: MonteCarloChoose

    implicit none

    integer,intent(out),OPTIONAL :: nReject

    real:: m,minmass,maxmass,bwd,m2,ymin,ymax,y,mres0,gamres0,gamtot,spectral,intfac,maxbwd
    real:: gamtot_max,spectral_max,intfac_max
    real:: mres02,gamres024
    integer:: ncount
    integer, parameter:: ncountm = 100
    logical:: lcount

    integer, parameter :: nB = 8
    real, dimension(nB) :: mbound,bweight,bmaxbwd, mbound2
    integer :: i,iB
    real :: bweightTot

    if (present(nReject)) nReject=0

    minMass = hadron(Delta)%minmass
    maxmass = 3.0

    mres0=hadron(Delta)%mass
    gamres0=FullWidthBaryon(Delta,mres0)

    mres02   = mres0**2
    gamres024= gamres0**2/4

    ! Divide the mass region into slices, calculate values at the borders:

    mbound = (/minMass,mres0-2*gamres0,mres0-gamres0,mres0-0.25*gamres0,mres0,&
         & mres0+gamres0,mres0+2*gamres0,maxMass/)

    do i=1,nB
       if (mbound(i).lt.minMass) mbound(i)=minMass
       if (mbound(i).gt.maxMass) mbound(i)=maxMass
    end do

    mbound2 = mbound**2

    do i=1,nB
       bweight(i) = atan2(2*(mbound(i)-mres0),gamres0)

       gamtot_max   = FullWidthBaryon(Delta,mbound(i))
       spectral_max = mbound2(i)*gamtot_max*gamres0/((mres02-mbound2(i))**2 + gamtot_max**2*mbound2(i))
       intfac_max   = gamres024/((mbound(i)-mres0)**2 + gamres024)
       bmaxbwd(i)   = spectral_max/intfac_max * 1.05
    end do

    ! Calculate the integral of the Cauchy distribution and their
    ! modification:

    do i=1,nB-1
       bweight(i)=(bweight(i+1)-bweight(i))
       bmaxbwd(i) = max(bmaxbwd(i),bmaxbwd(i+1))
       bweight(i) = bweight(i) * bmaxbwd(i) ! modify the Cauchy weights
    end do

    ! Now do the actual MC:

    ncount=0
    lcount=.true.
    do while(lcount)
       ncount=ncount+1

       ! STEP 1: Select the mass region:
       ! (it is important to have this step in the overall loop, since
       ! this cures inaccuracies in the calculation of the weights)

       iB=MonteCarloChoose(bweight(1:nB-1), bweightTot)
       ymax=2.*atan2(2*(mbound(iB+1)-mres0),gamres0)
       ymin=2.*atan2(2*(mbound(iB)-mres0),gamres0)
       maxbwd = bmaxbwd(iB)

       ! STEP 2: generate random value according Cauchy with constant width

       y=ymin+rn()*(ymax-ymin)
       m=.5*tan(y/2.)*gamres0+mres0
       m=min(max(m,minmass),maxmass)
       m2=m**2

       ! STEP 3: Do the rejection

       gamtot=FullWidthBaryon(Delta,m)
       spectral=m2*gamtot*gamres0/((mres02-m2)**2 + gamtot**2*m2)
       intfac=gamres024/((m-mres0)**2 + gamres024)
       bwd=spectral/intfac/maxbwd

       if (bwd>1) then
          if (DoPr(1)) then
             write(*,'(A,2F12.3,A,2G12.3)') 'Problems in Delta_Mass: (bwd>1)  mass=',&
                  & mbound(iB),mbound(iB+1),' bwd:',bwd,maxbwd
             write(*,'(A,1P,10g12.3)')'If this problem occurs frequently, adjust slicing !!!'
             write(*,*)
          end if
       end if

       if (rn()<=bwd) lcount =.false.

       if (ncount>=ncountm) then
          if (DoPr(1)) then
             write(*,'(A,2G12.3)') 'Problems in Delta_Mass: (ncount)  mass=',minmass,maxmass
          end if
          m=mres0
          lcount=.false.
       end if
    end do

    if (present(nReject)) nReject=ncount
    Delta_Mass_Full = m

  end function Delta_Mass_Full



end module PythiaSpecFunc
