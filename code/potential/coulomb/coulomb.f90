!******************************************************************************
!****m* /Coulomb
! NAME
! module Coulomb
! PURPOSE
! Includes all routines and global variables which are necessary for
! the evaluation of a Coulomb potential out of a charge configuration.
! NOTES
! * The coulomb routines use now the same spatial grid as that
!   extracted in the module "densityModule".
! * New ADI routines are used to solve the Poisson equation on the grid.
!******************************************************************************
module coulomb

  use constants, only: singlePrecision

  implicit none

  private

  !****************************************************************************
  !****g* Coulomb/coulombField
  ! SOURCE
  real(singlePrecision), save, allocatable, dimension(:,:,:,:) :: coulombField
  ! PURPOSE
  ! The Coulomb field A^{\mu}(\vec{r}) = (\phi(\vec{r}),\vec{A}(\vec{r}))
  ! Indizes:
  ! * 1-3: components on the 3D-Grid: x,y,z--coordinate
  ! * 4: Lorentz-components: 0th entry is the static coulomb field,
  !   while entries 1:3 represent the magnet field.
  !   (The magnetic filed is not yet implemented !)
  !****************************************************************************

  !****************************************************************************
  !****g* Coulomb/coulombFlag
  ! SOURCE
  logical, save :: coulombFlag = .false.
  ! PURPOSE
  ! Switch to turn on/off the Coulomb potential.
  ! If turned on, also 'symmetriePotFlag' (in namelist 'baryonPotential')
  ! needs to be turned on.
  !****************************************************************************

  !****************************************************************************
  !****g* Coulomb/magnetFieldFlag
  ! SOURCE
  logical, save :: magnetFieldFlag=.false.
  ! PURPOSE
  ! Switch to turn on/off elm. vector potential.
  ! NOTES
  ! The vector potential is not yet fully implemented! Please do not use
  ! this option.
  !****************************************************************************

  !****************************************************************************
  !****g* Coulomb/cutMomentumPotential
  ! SOURCE
  real, save :: cutMomentumPotential = 0.025 ! in GeV
  ! PURPOSE
  ! If larger than 0, the coulomb potential is set to zero for all particles
  ! with momentum larger than minmass**2/(2*cutMomentumPotential)
  ! The cut-off is smeared out, if cutMomentumWidth>0
  !
  ! This cut is needed in non-RMF mode since
  !   m_eff^2 = (sqrt(m^2+p^2)+U_C)^2-p^2 can
  ! become smaller than zero for
  !   p > U_C/2 - m^2/(2*U_C).
  ! In this case we have a 'tachyon'.
  !
  ! NOTES
  ! * for RMF mode you do not need this cut
  ! * This value should correspond to the maximum vale of the Coulomb
  !   potential. Therefore you should readjust this for every nucleus.
  ! * For the pion, we take the mass (0.138) instead of minmass, since here
  !   minmass is zero!
  !****************************************************************************

  !****************************************************************************
  !****g* Coulomb/cutMomentumWidth
  ! SOURCE
  real, save :: cutMomentumWidth = 0.100 ! in GeV
  ! PURPOSE
  ! If cutMomentumPotential>0 and cutMomentumWidth>0, then the cut-off is
  ! smeared by some linear interpolation:
  !   (At-Width) = 1.0 ... (At+Width) = 0.0
  ! with At = minmass**2/(2*cutMomentumPotential)
  ! and the width given here in GeV.
  !
  ! This width is introduced in order to destroy numerical problems due to
  ! sharp edges.
  !****************************************************************************


  !****************************************************************************
  !****g* Coulomb/K_ADI
  ! SOURCE
  real, save, allocatable, dimension(:,:,:) :: K_ADI
  ! PURPOSE
  ! auxiliary arrays needed for ADI routines (see also ADI.f90)
  !****************************************************************************

  !****************************************************************************
  !****g* Coulomb/U_ADI
  ! SOURCE
  real, save, allocatable, dimension(:,:,:) :: U_ADI
  ! PURPOSE
  ! auxiliary arrays needed for ADI routines (see also ADI.f90)
  !****************************************************************************

  !****************************************************************************
  !****g* Coulomb/chatot
  ! SOURCE
  real, SAVE :: chatot
  ! PURPOSE
  ! Total charge density (needed for calculating coulomb outside of grid)
  !****************************************************************************

  !****************************************************************************
  !****g* Coulomb/rchar
  ! SOURCE
  real, dimension(1:3), SAVE :: rchar
  ! PURPOSE
  ! Position of total charge (needed for calculating coulomb outside of grid)
  !****************************************************************************

  logical, save :: initreadingFlag = .true.

  real, parameter :: g_coulomb = 0.0014398  !=e^2*0.2 Gev fm=1/137*0.2 GeV fm

  integer, save :: kmax = 0 ! if MagnetFieldFLAG=true: =3


  public :: emfoca
!  public :: coulomb_ADI
  public :: cleanUp
  public :: updateCoulomb
  public :: getCoulombFlag


contains


  !****************************************************************************
  !****is* Coulomb/getCoulombFlag
  ! NAME
  ! logical function getCoulombFlag()
  ! PURPOSE
  ! Return the value of 'coulombFlag'.
  !****************************************************************************
  logical function getCoulombFlag()
    if (initreadingFlag) call initCoulombReading
    getCoulombFlag = coulombFlag
  end function getCoulombFlag


  !****************************************************************************
  !****is* Coulomb/initCoulombReading
  ! NAME
  ! subroutine initCoulombReading
  ! PURPOSE
  ! initialize the module, read namelist 'coulomb'
  !****************************************************************************
  subroutine initCoulombReading
    use output, only: Write_ReadingInput
    use nucleusDefinition, only: tNucleus
    use nucleus, only: getTarget
    use inputGeneral, only: eventType
    use eventTypes, only: elementary, Box
    use densityModule, only: getGridPoints

    integer :: ios, gridPoints(1:3)
    type(tNucleus),pointer::target

    !**************************************************************************
    !****n* Coulomb/coulomb
    ! NAME
    ! NAMELIST /coulomb/
    ! PURPOSE
    ! Includes the switches:
    ! * coulombFlag
    ! * magnetFieldFlag
    ! * cutMomentumPotential
    ! * cutMomentumWidth
    !**************************************************************************
    NAMELIST /coulomb/ coulombFlag, magnetFieldFlag, cutMomentumPotential, &
         cutMomentumWidth

    initreadingFlag=.false.

    call Write_ReadingInput('coulomb',0)
    rewind(5)
    read(5,nml=coulomb,iostat=ios)
    call Write_ReadingInput('coulomb',0,ios)

    select case (eventType)
    case (elementary, Box)
       if (coulombFlag) then
          coulombFlag = .false.
          write(*,*) 'elementary events/box => Disabling coulomb flag!'
       end if
    case default
       target => getTarget()
       if (target%mass<3 .and. coulombFlag) then
          coulombFlag = .false.
          write(*,*) 'Hydrogen target => Disabling coulomb flag!'
       end if
    end select

    write(*,*) 'coulomb flag =     ',coulombFlag
    if (coulombFlag) then
       write(*,*) 'magnetfield flag = ',magnetFieldFlag

       if (cutMomentumPotential > 0) then
          write(*,*) 'momentum cut = m^2/(2*', cutMomentumPotential,') +- ',&
               cutMomentumWidth
       else
          write(*,*) 'momentum cut : no'
       end if
    end if


    call Write_ReadingInput('coulomb',1)

    if (.not.coulombFlag) RETURN !exit

    if (magnetFieldFlag) then
       magnetFieldFlag = .false.
       write(*,*)
       write(*,*) '!!! WARNING: Magnet fields not yet implemented!'
       write(*,*) '!!! Thus ignoring the input value and using'
       write(*,*) '!!! magnetFieldFLAG = ',magnetFieldFLAG
       write(*,*)
    end if

    if (magnetFieldFlag) kmax = 3 ! see above !!!!!

    ! allocate fields
    gridPoints = getGridPoints()

    allocate(coulombField(-gridPoints(1):gridPoints(1), &
                          -gridPoints(2):gridPoints(2), &
                          -gridPoints(3):gridPoints(3), 0:kmax))

    allocate(K_ADI(-gridPoints(1):gridPoints(1), &
                   -gridPoints(2):gridPoints(2), &
                   -gridPoints(3):gridPoints(3)))

    allocate(U_ADI(-gridPoints(1):gridPoints(1), &
                   -gridPoints(2):gridPoints(2), &
                   -gridPoints(3):gridPoints(3)))

  end subroutine initCoulombReading


  !****************************************************************************
  !****s* Coulomb/cleanUp
  ! subroutine cleanUp
  ! PURPOSE
  ! Deallocate all fields
  !****************************************************************************
  subroutine cleanUp
    use ADI, only: cleanupADI => cleanup

    if (.not.coulombFlag) return
    if (allocated(coulombField)) deallocate(coulombField)
    if (allocated(K_ADI)) deallocate(K_ADI,U_ADI)
    call cleanupADI
  end subroutine cleanUp


  !****************************************************************************
  !****f* Coulomb/emfoca
  ! NAME
  ! real function emfoca(position, momentum, chargeIn, ID, emForce)
  ! PURPOSE
  ! This routine calculates the electromagnetic forces acting on a particle
  ! (baryon or meson).
  ! INPUTS
  ! * real, dimension(1:3) :: position -- space coordinates of particle
  ! * real, dimension(1:3) :: momentum -- momentum of particle
  ! * integer              :: chargeIn -- charge of particle
  ! * integer, OPTIONAL    :: ID -- ID of particle
  ! RESULT
  ! * real                           :: cpot -- electromagn. potential (in GeV)
  ! * real, dimension(1:3), OPTIONAL :: emForce-- electromagn. force (in GeV/fm)
  !****************************************************************************
  real function emfoca(position, momentum, chargeIn, ID, emForce) result (cpot)

    use IdTable, only: isHadron
    use ParticleProperties, only: hadron
    use densityModule, only: GetGridIndex, gridSpacing
    use CallStack, only: Traceback

    real,    dimension(1:3) ,intent(in)            :: position, momentum
    integer,                 intent(in)            :: chargeIn, ID
    real,    dimension(1:3) ,intent(out), optional :: emForce

    real    :: charge, relr
    integer :: ipkt, ipktmax
    real    :: funlin(1:8), deriv(1:7), coord(1:7,1:3)
    integer, dimension(1:3) :: ic,io,iu
    real, dimension(1:3) :: lin,rel

    real :: pAbs = 0, facP, CutMomentumAt


    if (initreadingFlag) call initCoulombReading

    cpot=0.
    if (present(emForce)) emForce=0.

    if (.not.CoulombFlag) return
    if (chargeIn.eq.0) return

    if (CutMomentumPotential > 0) then ! momentum cut is active
       if (isHadron(ID)) then
          if (ID==101) then
             CutMomentumAt = hadron(ID)%mass**2/(2*CutMomentumPotential)
          else
             CutMomentumAt = hadron(ID)%minmass**2/(2*CutMomentumPotential)
          end if

          if (CutMomentumAt < CutMomentumWidth) then
             write(*,*) 'ID = ', ID
             write(*,*) 'CutMomentumAt    = ',CutMomentumAt
             write(*,*) 'CutMomentumWidth = ',CutMomentumWidth
             call Traceback('CutMomentumAt < CutMomentumWidth')
          end if

          pAbs = sqrt(Dot_product(momentum,momentum))
          if (pAbs > CutMomentumAt+CutMomentumWidth) return ! potential is zero
       else
          CutMomentumAt = -99.9
       end if
    end if

    charge = float(chargeIn)

    ipktmax = 1
    if (present(emForce)) ipktmax=7

    if (GetGridIndex(position,ic,-2)) then
       ! *---------------------------------------------------------------------
       ! inside the grid

       coord(1,1:3) = position

       if (present(emForce)) then
          do ipkt=2,ipktmax
             coord(ipkt,1:3) = position
          end do
          coord(2,1) = coord(2,1) + gridSpacing(1)
          coord(3,1) = coord(3,1) - gridSpacing(1)
          coord(4,2) = coord(4,2) + gridSpacing(2)
          coord(5,2) = coord(5,2) - gridSpacing(2)
          coord(6,3) = coord(6,3) + gridSpacing(3)
          coord(7,3) = coord(7,3) - gridSpacing(3)
       end if

       do ipkt = 1,ipktmax

          iu = int(coord(ipkt,1:3)/gridSpacing)
          io = iu+1

          funlin(1) = coulombField(iu(1),iu(2),iu(3),0)
          funlin(2) = coulombField(io(1),iu(2),iu(3),0)
          funlin(3) = coulombField(io(1),io(2),iu(3),0)
          funlin(4) = coulombField(iu(1),io(2),iu(3),0)
          funlin(5) = coulombField(iu(1),iu(2),io(3),0)
          funlin(6) = coulombField(io(1),iu(2),io(3),0)
          funlin(7) = coulombField(io(1),io(2),io(3),0)
          funlin(8) = coulombField(iu(1),io(2),io(3),0)


          lin = (coord(ipkt,1:3)-float(iu)*gridSpacing)/gridSpacing

          deriv(ipkt) =  &
               &  (1.-lin(1))*(1.-lin(2))*(1-lin(3))*funlin(1) + &
               &      lin(1) *(1.-lin(2))*(1-lin(3))*funlin(2) + &
               &      lin(1) *    lin(2) *(1-lin(3))*funlin(3) + &
               &  (1.-lin(1))*    lin(2) *(1-lin(3))*funlin(4) + &
               &  (1.-lin(1))*(1.-lin(2))*   lin(3) *funlin(5) + &
               &      lin(1) *(1.-lin(2))*   lin(3) *funlin(6) + &
               &      lin(1) *    lin(2) *   lin(3) *funlin(7) + &
               &  (1.-lin(1))*    lin(2) *   lin(3) *funlin(8)

       end do

       cpot = charge*deriv(1)

       if (present(emForce)) then
          emForce(1) = deriv(2)-deriv(3)
          emForce(2) = deriv(4)-deriv(5)
          emForce(3) = deriv(6)-deriv(7)

          emForce = -emForce *charge/(2*gridSpacing)
       end if

       ! *---------------------------------------------------------------------
    else
       ! *---------------------------------------------------------------------
       ! outside the grid
       rel = position-rchar
       relr= sum(rel**2)

       if (relr .lt. 1.0) then
          !          write(*,*)'problems in emfoca; relr**2 = ', relr
          !          write(*,*)'particle coordinates :', position
          !          write(*,*)'center of charge :    ', rchar
          relr = 1.0
          !          stop
       end if

       relr   = sqrt(relr)
       cpot   = charge * g_coulomb * chatot/relr

       if (present(emForce)) then
          emForce = cpot * rel/relr**2
       end if
       ! *---------------------------------------------------------------------
    end if

    if (CutMomentumPotential > 0) then ! momentum cut is active
       if (CutMomentumAt > 0) then ! some true value is given
          if (pAbs < CutMomentumAt-CutMomentumWidth) return ! potential is one

          if (CutMomentumWidth>0) then
             facP = (1.-(pAbs-CutMomentumAt)/CutMomentumWidth)/2

             cpot = cpot * facP
             if (present(emForce)) emForce = emForce * facP
          end if
       end if
    end if

  end function emfoca



  !****************************************************************************
  !****s* Coulomb/coulomb_ADI
  ! NAME
  ! subroutine coulomb_ADI
  ! PURPOSE
  ! This routine calculates coulomb field on the grid using ADI methods
  !****************************************************************************
  subroutine coulomb_ADI
    use constants, only: pi
    use ADI
    use densityModule, only: get_densitySwitch, densityField, gridSpacing
    use output, only: Write_InitStatus

    integer :: k
    real,    SAVE :: eta_x, eta_y
    logical, SAVE :: init=.true.

    if (get_densitySwitch()==0) then ! --- no density:
       coulombField   = 0.
       init=.false.
       return
    end if

    if (init) then
       call Write_InitStatus("Coulomb_ADI",0)

       coulombField   = 0.
       call bounds
       eta_x = (gridSpacing(3)/gridSpacing(1))**2
       eta_y = (gridSpacing(3)/gridSpacing(2))**2
       init=.false.
    end if

    do k = 0,kmax ! cf. magnetFieldFlag

       U_ADI = coulombField(:,:,:,k)
       K_ADI = densityField(:,:,:)%proton(k) * gridSpacing(3)**2*4*pi*g_coulomb

       call ADI_Coulomb(U_ADI, K_ADI, eta_x, eta_y)

       coulombField(:,:,:,k) = U_ADI

    end do
    call bounds

    if (init) then
       call Write_InitStatus("Coulomb_ADI",1)
       init=.false.
    end if

  end subroutine coulomb_ADI




  !****************************************************************************
  !****s* Coulomb/Get_TotalCharge
  ! NAME
  ! subroutine Get_TotalCharge
  ! PURPOSE
  ! This routine calculates the total charge Q & its center R_Q.
  !****************************************************************************
  subroutine Get_TotalCharge

    use dichteDefinition
    use densityModule, only: DensityAt, gridSpacing, gridPoints

    real    :: dV
    integer :: izc,iyc,ixc
    real, dimension(1:3) :: pos
    type(dichte) :: dens
    logical, save :: verbose=.false.

    chatot = 0.0
    rchar = 0.0

    dV = gridSpacing(1)*gridSpacing(2)*gridSpacing(3)


       do izc = -gridPoints(3),gridPoints(3)
          pos(3) = float(izc)*gridSpacing(3)
          do iyc = -gridPoints(2),gridPoints(2)
             pos(2)= float(iyc)*gridSpacing(2)
             do ixc = -gridPoints(1),gridPoints(1)
                pos(1) = float(ixc)*gridSpacing(1)

                dens = DensityAt(pos)
                rchar  = rchar  + pos*dens%proton(0)
                chatot = chatot +     dens%proton(0)

             end do
          end do
       end do


!!$    if(chatot.lt. 1e-3) then
!!$       write(*,*)'error in chargeDensity total charge = ', chatot
!!$       stop
!!$    end if

    if (abs(chatot).gt.0.001) then
       rchar = rchar/chatot
    else
       rchar = 0.
    end if

    chatot = chatot * dV

    if (verbose) then
       write(*,*) 'testing Get_TotalCharge :'
       write(*,*) 'total charge = ',chatot
       write(*,*) 'center of charge: ',rchar(1), rchar(2), rchar(3)
       write(*,*) 'am ende von Get_TotalCharge'
    end if

  end subroutine Get_TotalCharge



  !****************************************************************************
  ! **     berechung der bounds fuer die poisson-gleichung
  ! **     es wird angenommen, dass die dichte- bzw. die stromverteilung auf
  ! **     einem dreidimensionalen gitter gegeben ist.
  ! Routine taken and modified from old coulombroutines
  ! (see S. Theis, PhD-Thesis)
  ! ************************************************************************
  subroutine bounds

    use dichteDefinition
    use densityModule, only: densityField, gridSpacing, gridPoints

    real,save :: mhat0(0:3)
    real,save :: mhat1(0:3,1:3)
    real,save :: mhat2(0:3,1:3,1:3)
    real,save :: mhat3(0:3,1:3,1:3,1:3)

    real :: xpri(1:3), rcell, rabs, term1, term2, term3, term4
    integer :: i, j, k, l, m, n, o

    real, save :: dVol=0.
    logical, save :: lflag=.true.

    if (lflag) then
       dVol = gridSpacing(1)*gridSpacing(2)*gridSpacing(3)
       lflag = .false.
    end if

    ! **  felder auf null setzen

    mhat0 = 0.0
    mhat1 = 0.0
    mhat2 = 0.0
    mhat3 = 0.0


    ! **  calculate the tensors mhat
    do l = 0, kmax
       do i = -gridPoints(3),gridPoints(3)
          xpri(3) = float(i)*gridSpacing(3)
          do j = -gridPoints(2), gridPoints(2)
             xpri(2) = float(j)*gridSpacing(2)
             do k = -gridPoints(1), gridPoints(1)
                xpri(1) = float(k)*gridSpacing(1)

                rabs    = sqrt(sum(xpri**2))

                ! ***   evaluate mhats explicitely
                rcell  = densityField(k,j,i)%proton(l)

                mhat0(l) = mhat0(l) + rcell

                do m = 1, 3
                   mhat1(l,m)    = mhat1(l,m)    + rcell*xpri(m)
                   do n = 1, 3

                      if (n.eq.m) then
                         mhat2(l,m,n)=mhat2(l,m,n)+rcell*(3.0*xpri(m)*xpri(n)-rabs**2)
                      else
                         mhat2(l,m,n)=mhat2(l,m,n)+rcell*(3.0*xpri(m)*xpri(n))
                      end if


                      do o = 1, 3
                         if (n.eq.m) then
                            mhat3(l,m,n,o)= mhat3(l,m,n,o)+              &
                                 &  rcell*(5.0*xpri(m)*xpri(n)*xpri(o)      &
                                 &  -3.0*rabs**2*xpri(o))
                         else
                            mhat3(l,m,n,o)= mhat3(l,m,n,o)+rcell*(5.0*xpri(m)*xpri(n)*xpri(o))
                         end if

                      end do
                   end do
                end do
             end do
          end do
       end do
    end do

    mhat0 = mhat0*dVol
    mhat1 = mhat1*dVol
    mhat2 = mhat2*dVol
    mhat3 = mhat3*dVol

    ! ******** set potential to zero and calculate pots at bounds **********

    do l = 0, kmax

       do i = -gridPoints(3), gridPoints(3)
          do j = -gridPoints(2), gridPoints(2)
             do k = -gridPoints(1), gridPoints(1)

                if (  abs(i).eq.gridPoints(3) .or.  &
                     abs(j).eq.gridPoints(2) .or.  &
                     abs(k).eq.gridPoints(1)) then

                   xpri = float((/k,j,i/))*gridSpacing

                   rabs = sqrt(sum(xpri**2))
                   if (rabs.le.(0.000001)) then
                      write(*,*) 'mistake while calc. of rabs'
                      stop
                   end if


                   term1 = mhat0(l)/rabs

                   term2 = 0.0
                   do m = 1,3
                      term2 = term2 + mhat1(l,m)*xpri(m)
                   end do
                   term2 = term2/(rabs**3)

                   term3 = 0.0
                   do m = 1,3
                      do n = 1,3
                         term3 = term3 + mhat2(l,m,n)*xpri(m)*xpri(n)
                      end do
                   end do

                   term3 = term3/(2.0*rabs**5)

                   term4 = 0.0
                   do m = 1,3
                      do n = 1,3
                         do o = 1,3
                            term4 = term4 + mhat3(l,m,n,o)*xpri(m)*xpri(n)*xpri(o)
                         end do
                      end do
                   end do
                   term4 = term4/(2.0*rabs**7)

                   coulombField(k,j,i,l) = (term1+term2+term3+term4)*g_coulomb

                end if

             end do
          end do
       end do
    end do


  end subroutine bounds

  !****************************************************************************
  !****s* Coulomb/updateCoulomb
  ! NAME
  ! subroutine updateCoulomb
  ! PURPOSE
  ! Updates the coulomb field values. Necessary to call it in each time step
  ! INPUTS
  ! none
  !****************************************************************************
  subroutine updateCoulomb
    use output, only: DoPR
    use densityModule, only: get_densitySwitch
    logical, save :: SkipThis = .false.

    if (initreadingFlag) call initCoulombReading

    if (.not.CoulombFlag) return
    if (SkipThis) return

    if (DoPr(2)) write(*,*) 'Updating Coulomb'

    !determine coulomb field on the grid
    call coulomb_ADI

    ! determine total charge & center of total charge
    ! (needed in calculating the Coulomb force acting an particles
    ! outside of the grid.
    call Get_TotalCharge

    if (get_densitySwitch()==2) SkipThis = .true. ! never do an update again


  end subroutine updateCoulomb



end module coulomb
