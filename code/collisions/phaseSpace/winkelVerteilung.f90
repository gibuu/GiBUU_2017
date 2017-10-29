!******************************************************************************
!****m* /winkelVerteilung
! NAME
! module winkelVerteilung
! PURPOSE
! Incorporates the functions used for angular distributions.
!******************************************************************************
module winkelVerteilung

  implicit none
  private

  !****************************************************************************
  !****g* winkelVerteilung/deltaPWave
  ! SOURCE
  !
  logical, save :: deltaPWave=.true.
  ! PURPOSE
  ! Switch for P-Wave decay of delta in pion nucleon
  ! Only relevant for deltas which are produced in pion-nucleon collisions.
  ! -> see also master_2body
  ! Values:
  ! * .false.= isotropic in CM-Frame
  ! * .true. = 1+3*cos(theta)**2 in CM Frame (theta is angle of producing
  !   pion to outgoing pion)
  !****************************************************************************

  !****************************************************************************
  !****g* winkelVerteilung/rho_pipi_nonIsotropic
  ! SOURCE
  !
  logical, save :: rho_pipi_nonIsotropic=.true.
  ! PURPOSE
  ! Switch for non-isotropic rho -> pi pi decay:
  ! * .false.= isotropic in CM-Frame
  ! * .true. = non-isotropic
  !****************************************************************************

  !****************************************************************************
  !****g* winkelVerteilung/pionNucleon_backward
  ! SOURCE
  !
  logical, save :: pionNucleon_backward=.true.
  ! PURPOSE
  ! Switch for backward peaked pion nucleon cross section:
  ! * .true.= use backward peaked distribution
  ! * .false.= isotropic
  !****************************************************************************

  !****************************************************************************
  !****g* winkelVerteilung/pionNucleon_backward_exponent
  ! SOURCE
  !
  real, save :: pionNucleon_backward_exponent=26.5
  ! PURPOSE
  ! Exponent for backward peaked pion nucleon cross section.
  ! Distribution=(coeff-cos(theta))**exponent*(pole-sqrt(s)/pole)
  ! Only used if pionNucleon_backward=.true. .
  !****************************************************************************

  !****************************************************************************
  !****g* winkelVerteilung/pionNucleon_backward_coeff
  ! SOURCE
  !
  real, save :: pionNucleon_backward_coeff=1.9
  ! PURPOSE
  ! Exponent for backward peaked pion nucleon cross section.
  ! Distribution=(coeff-cos(theta))**exponent*(pole-sqrt(s)/pole)
  ! Only used if pionNucleon_backward=.true. .
  !****************************************************************************

  !****************************************************************************
  !****g* winkelVerteilung/NNisotropic
  ! SOURCE
  !
  logical, save :: NNisotropic=.false.
  ! PURPOSE
  ! if .true.: set isotropic nucleon-nucleon elastic cross section
  !****************************************************************************

  !****************************************************************************
  !****g* winkelVerteilung/iParam_gammaNVN
  ! SOURCE
  !
  integer, save :: iParam_gammaNVN = 3
  ! PURPOSE
  ! for gamma N -> V N events, this parameter is given to the routine
  ! vecmesa and selects there, how dsigma/dt is calculated.
  ! Only if iParam_gammaNVN >= 0 the default value of that routine
  ! is overwritten.
  !
  ! Possible values:
  ! * 0: 'old' parametrisation for gammaN->VN (cf. Effenberger PhD):
  !   dsigma/dt ~ exp(Bt). Slope paremeter B according ABBHHM collab, PR 175, 1669 (1968).
  ! * 1: Pythia parametrisation:
  !   Slope parameter B=2*b_p+2*b_V+4*s**eps-4.2
  ! * 2: 'Donnachie, Landshoff'
  !   Select t according dsig/dt as given by VecMesWinkel/dsigdt, not
  !   by a given slope parameter
  ! * 3: as 1, but for rho and W<~6GeV slope parameter adjusted
  !   according CLAS experimental data [Morrow et al, EPJ A39, 5 (2009)]
  ! * 4: Muehlich PhD, Appendix E
  ! * 5: Rho0 Toy Init
  ! * 6: Rho0 Toy Init: Fit to PYTHIA-VMD
  ! * 7: Flat (not exp.)
  !
  ! cf. VecMesWinkel/vecmesa for a detailed description.
  !****************************************************************************

  !****************************************************************************
  !****g* winkelVerteilung/NN_NR_noniso
  ! SOURCE
  !
  logical, save :: NN_NR_noniso = .false.
  ! PURPOSE
  ! If .true., use non-isotropic angular distr. for NN -> NR,
  ! according to dsigma/dt = b/t**a.
  !****************************************************************************

  !****************************************************************************
  !****g* winkelVerteilung/debug
  ! SOURCE
  logical, parameter:: debug=.false.
  ! PURPOSE
  ! Switch for debug information
  !****************************************************************************

  logical, save :: initFlag = .true.

  public :: winkel, getIParam

contains

  !****************************************************************************
  !****is* winkelVerteilung/readInput
  ! NAME
  ! subroutine readInput
  ! PURPOSE
  ! Reads input in jobcard from namelist "angular_distribution"
  !****************************************************************************
  subroutine readInput
    use output, only: Write_ReadingInput

    !**************************************************************************
    !****n* winkelVerteilung/angular_distribution
    ! NAME
    ! NAMELIST angular_distribution
    ! PURPOSE
    ! Includes the switches:
    ! * deltaPWave
    ! * pionNucleon_backward
    ! * pionNucleon_backward_exponent
    ! * pionNucleon_backward_coeff
    ! * rho_pipi_nonIsotropic
    ! * NNisotropic
    ! * iParam_gammaNVN
    ! * NN_NR_noniso
    !**************************************************************************
    NAMELIST /angular_distribution/ deltaPWave, rho_pipi_nonIsotropic, NNisotropic, iParam_gammaNVN, &
                                    pionNucleon_Backward, pionNucleon_backward_exponent, pionNucleon_backward_coeff, &
                                    NN_NR_noniso

    integer :: IOS

    ! Read special angular distribution switches:
    call Write_ReadingInput('angular_distribution',0)
    rewind(5)
    read(5,nml=angular_distribution,IOSTAT=IOS)
    call Write_ReadingInput('angular_distribution',0,IOS)

    write(*,*) 'Switch for Delta P-Wave decay:',deltaPWave
    write(*,*) 'Switch to have low energy N pi-> N pi backward peaked:',pionNucleon_Backward
    if (pionNucleon_Backward) &
       write(*,'(A,F7.3,A,F7.3,A)') ' Distribution for N pi-> N pi = (',pionNucleon_Backward_coeff, &
                                    '- cos(theta))^(',pionNucleon_Backward_exponent,'(pole-sqrt(s))/pole)'
    write(*,*) 'rho-> pi pi: rho_pipi_nonIsotropic=',rho_pipi_nonIsotropic
    write(*,*) 'NN -> NN: NNisotropic=',NNisotropic

    write(*,*) 'gammaN -> VN: Slope-Param = ',iParam_gammaNVN
    write(*,*) 'NN -> NR: non-isotropic? ', NN_NR_noniso

    call Write_ReadingInput('angular_distribution',1)

    initFlag = .false.
  end subroutine readInput


  integer function getIParam()
    if (initFlag) call readInput
    getIParam = iParam_gammaNVN
  end function


  !****************************************************************************
  !****f*  winkelVerteilung/winkel
  ! NAME
  ! function winkel (teilchenIn, teilchenOut, srts, betaToCM, mediumAtCollision, successflag) result (pscatt)
  !
  ! PURPOSE
  ! This subroutine determines the scattering angle for two particles in
  ! the final state.
  ! For 2->2 and 1->2 processes.
  !
  ! INPUTS
  ! * type(particle) :: teilchenIn(1:2)     -- Initial state particles,
  !   If teilchenIn(2)%ID=0 then it's just a decay of teilchenIn(1)
  ! * type(particle) :: teilchenOut(1:2)    -- Final state particles
  !   (only their Id's,antiflags,charges and masses are needed here)
  ! * real           :: srts                -- SQRT(s) in vacuum
  ! * real, dimension(1:3) :: betaToCM      -- Center of mass velocity in
  !   calc. frame
  ! * type(medium)   :: MediumAtCollision   -- medium information
  !
  ! OUTPUT
  ! * real, dimension(1:3):: pscatt         -- unit vector in the direction
  !   of the outgoing particle teilchenOut(1) in the CM frame
  ! * logical             :: successflag    -- .true. is routine was successful
  !****************************************************************************
  function winkel (partIn, partOut, srts, betaToCM, mediumAtColl, successflag) result (pscatt)

    use mediumDefinition
    use particleDefinition
    use IdTable, only: pion, rho, omegaMeson, phi, JPsi, Nucleon, Delta, photon, isBaryon
    use Random, only: rnOmega

    type(particle),dimension(1:2),intent(in) :: partIn,partOut
    real, intent(in) :: srts
    real, dimension(1:3), intent(in) :: betaToCM
    type(medium), intent(in) :: mediumAtColl
    logical, intent(out),optional :: successflag
    real, dimension(1:3) :: pscatt

    integer :: id1,id2,id3,id4,iii1,iii2,ooo1,ooo2
    real :: m3,m4
    logical :: successVN,successND

    if (initFlag) call readInput

    if (present(successflag)) successflag = .true.

    id1=partIn(1)%ID
    id2=partIn(2)%ID
    id3=partOut(1)%Id
    id4=partOut(2)%Id
    m3=partOut(1)%mass
    m4=partOut(2)%mass

    if (id2 == 0) then
       ! ******************PARTICLE DECAYS*************************
       if (id1==Delta .and. id3==pion .and. deltaPWave) then
          ! Delta -> pi N decay in a P-Wave mode if delta was produced in pion-nucleon->Delta scattering
          if (debug) write(*,*) ' Delta -> pi N '
          pscatt = pscatt_Delta_PWave (partIn(1), srts)
       else if ((id1==rho) .and. (id3==pion.and.id4==pion) .and. rho_pipi_nonIsotropic) then
          ! rho-> pi pi decay
          if (debug) write(*,*) ' rho -> pi pi '
          pscatt = pscatt_Rho_pipi_Distribution (partIn(1))
       else
          !  Isotropic decay:
          if (debug) write(*,*) ' Isotropic decay: ',id1,' -> ',id3,id4
          pscatt = rnOmega()
       end if

    else
       ! ******************2-BODY REACTIONS*************************

       ! please note: mesons have always larger ID than baryons!
       iii1 = min(id1,id2)
       iii2 = max(id1,id2)
       ooo1 = min(id3,id4)
       ooo2 = max(id3,id4)

       if ((iii2==pion .and. iii1==nucleon) .and. &
           (ooo2==pion .and. ooo1==nucleon) .and. &
           pionNucleon_Backward .and. (srts<1.35)) then

          ! pi N -> pi N                 with sqrt(s).lt.1.35
          if (debug) write(*,*) ' pi N -> pi N with srts<1.35'
          pscatt = pscatt_piN_piN_backward (partIn, srts)
          ! assuming that pscatt is later attributed to id3:
          if (id4==pion) pscatt=-pscatt

       else if (((ooo1==nucleon).or.(ooo1==delta)) .and. &
                ((ooo2==rho).or.(ooo2==omegaMeson).or.(ooo2==phi).or.(ooo2==JPsi)) .and. &
                (iii1==nucleon) .and. (iii2==ooo2) ) then

          ! V N -> V N
          ! V N -> V Delta
          if (debug) write(*,*) ' V N -> V X'
          pscatt = pscatt_VN_VX (partIn, id3, id4, m3, m4, srts, mediumAtColl, successVN)
          if (present(successflag)) successflag=successVN

       else if ((iii1==nucleon) .and. (iii2==photon) .and. (ooo1==nucleon) .and. &
                ((ooo2==rho) .or. (ooo2==omegaMeson) .or. (ooo2==phi))) then

          ! gamma N -> V N
          pscatt = pscatt_gammaN_VN (partIn, id3, id4, m3, m4, srts, mediumAtColl)

       else if (isBaryon(iii1) .and. isBaryon(iii2)) then

          if (debug) write(*,*) ' Bar Bar -> Bar Bar'
          pscatt = pscatt_BarBar (partIn, partOut, srts, betaToCM, successND)
          if (present(successflag)) successflag=successND

       else

          !  Isotropic scattering
          if (debug) write(*,*) ' Isotropic scattering: ',id1,id2,' -> ',id3,id4
          pscatt = rnOmega()

       end if

    end if

    if (debug) write(*,*) 'in winkel:',pscatt,dot_product(pscatt,pscatt)

  end function winkel




  !############################################################################
  !############################################################################
  !   SPECIAL DISTRIBUTIONS
  !############################################################################
  !############################################################################



  !****************************************************************************
  !****if* winkelVerteilung/pscatt_Delta_PWave
  ! NAME
  ! function pscatt_Delta_PWave (teilchen, srts) result (pScatt)
  ! PURPOSE
  ! This routine determines the outgoing unit vector in the cm-frame for
  ! a meson out of p-wave Delta resonance decay.
  ! The angle should be distributed according to
  ! f(theta)=(3*cos^2(theta)+1)*(pionNucleon_Backward_coeff-cost)**(pionNucleon_Backward_exponent*(pole-srts)/pole)
  ! INPUTS
  ! * real           :: srts
  ! * type(particle) :: teilchen  --- The Delta Particle
  ! OUTPUT
  ! * real, dimension(1:3) :: pScatt --- Scattering vector in cm-frame
  !****************************************************************************
  function pscatt_Delta_PWave (teilchen, srts) result (pScatt)
    use particleDefinition
    use PIL_mesonMom, only: PIL_mesonMom_GET
    use random, only: rn, rnOmega
    use rotation, only: rotateTo
    use constants, only: pi, mN

    type(particle),intent(in)         :: teilchen                 ! The Delta Particle
    real, intent(in)                  :: srts
    real, dimension(1:3) :: pScatt

    real, dimension(1:3) :: mesonMomentum
    real :: cost, x, phi, sint
    real, parameter :: pole=1.232
    logical :: infoAvailable

    if (initFlag) call readInput

    infoAvailable = PIL_mesonMom_GET (teilchen%number, mesonMomentum)

    if ((.not. infoAvailable) .or. dot_product(mesonMomentum,mesonMomentum)<1E-15) then

       ! MesonMomentum is zero -> isotropic decay
       pScatt = rnOmega()

    else

       if (pionNucleon_backward) then
          ! (1a) Choose cos(theta) according to
          !  f(cost)=(3*cost**2+1)*(coeff-cost)**(exponent*(pole-srts)/pole)
          do
             cost=1-(rn()*2.)
             x=rn()*(1+pionNucleon_Backward_coeff)**(pionNucleon_Backward_exponent*(pole-mN)/pole)
             if (x<(3.*cost**2+1.)/4.*(pionNucleon_Backward_coeff-cost)**(pionNucleon_Backward_exponent*(pole-srts)/pole)) exit
          end do

       else
          ! (1b) Choose cos(theta) according to f(cost)=(3*cost**2+1)
          do
             cost=1-(rn()*2.)
             x=rn()
             if (x<(3.*cost**2+1.)/4.) exit
          end do

       end if

       if (debug) write(2,*) mesonMomentum,cost

       sint=sqrt(max(1.-cost**2,0.))
       phi=rn()*2.*pi

       pScatt(1)=cos(phi)*sint
       pScatt(2)=sin(phi)*sint
       pScatt(3)=cost

       ! (2) Rotate z-Axis on the momentum direction of the meson which produced the delta
       pScatt = rotateTo (mesonMomentum, pScatt)
    end if

  end function pscatt_Delta_PWave



  !****************************************************************************
  !****if* winkelVerteilung/pscatt_VN_VX
  ! NAME
  ! function pscatt_VN_VX (teilchenIn, id3, id4, mass3, mass4, srts, mediumAtCollision, successflag) result (pscatt)
  ! PURPOSE
  ! This routine determines a random unit vector for an outgoing vector
  ! meson in a V N-> V X colllision. (X stands for N or Delta)
  !
  ! INPUTS
  ! * type(particle) :: teilchenIn(1:2)  --- incoming particles
  ! * real,          :: srts             --- sqrt(s)
  ! * integer        :: id3, id4         --- ID of outgoing particles
  ! * real           :: mass3, mass4     --- mass of outgoing particle
  !
  ! OUTPUT
  ! * real, dimension(1:3) :: pscatt --- unit vector in outgoing baryon
  !   direction in CM frame
  !
  ! NOTES
  ! The given vector is the vector of particle given by id3, mass3.
  ! This is normally a baryon (we are normally treating B+M -> B+M
  ! collisions).
  !****************************************************************************

  function pscatt_VN_VX (partIn, id3, id4, mass3, mass4, srts, mediumAtColl, successflag) result (pscatt)

    use particleDefinition
    use VecMesWinkel, only: vecmesa, vecdelta
    use IDtable, only: nucleon, isMeson
    use rotation, only: rotateTo
    use constants, only: pi
    use random, only: rn
    use mediumDefinition
    use mesonPotentialModule, only: vecMes_massShift

    type(particle), dimension(1:2), intent(in) :: partIn
    integer, intent(in) :: id3,id4
    real, intent(in) :: mass3, mass4, srts
    logical, intent(out) :: successflag
    type(medium), intent(in) :: mediumAtColl
    real, dimension(1:3) :: pscatt

    real :: cost, sint, phi, massMes, massBar
    real :: massInMes2 ! mass of incoming meson mass squared or -Q^2
    integer :: idMes,idBar

!    write(*,*) '........pscatt_VN_VX'
!    call WriteParticle(6,1,1, teilchenIn(1))
!    call WriteParticle(6,1,2, teilchenIn(2))

    if (isMeson(partIn(1)%ID)) then
       massInMes2 = partIn(1)%mass
    else
       massInMes2 = partIn(2)%mass
    end if
    ! if %mass is really the mass, then we have to square it.
    ! if %mass < 0, then the really stored value was -Q^2
    if (massInMes2 > 0.) massInMes2 = massInMes2**2

    if (id3 < id4) then !
       idBar = id3
       idMes = id4
       massBar = mass3
       massMes = mass4
    else
       idBar = id4
       idMes = id3
       massBar = mass4
       massMes = mass3
    end if

    ! take care of mass shift
    massMes = massMes + vecMes_massShift(idMes,mediumAtColl%density)

    if (idBar==nucleon) then
       call vecmesa(srts,idMes,massMes,massInMes2,cost,iParam_gammaNVN,successflag=successflag)
    else
       call vecdelta(srts,idMes,massBar,massMes,massInMes2,cost,successflag=successflag)
    end if
    if (.not.successflag) return

    sint=sqrt(max(1.-cost**2,0.))
    phi=rn()*2.*pi
    pscatt(1)=cos(phi)*sint
    pscatt(2)=sin(phi)*sint
    pscatt(3)=cost

    ! Rotate according to incoming meson
    if (isMeson(partIn(1)%ID)) then
      pscatt = rotateTo (partIn(1)%momentum(1:3), pscatt)
    else
      pscatt = rotateTo (partIn(2)%momentum(1:3), pscatt)
    end if

    if (.not. isMeson(id3)) pscatt = -pscatt

!    write(*,'(0P,3f11.3)') pscatt

  end function pscatt_VN_VX



  !****************************************************************************
  !****if*  winkelVerteilung/pscatt_gammaN_VN
  ! NAME
  ! function pscatt_gammaN_VN (teilchenIn, id3, id4, mass3, mass4, srts, mediumAtCollision) result (pscatt)
  ! PURPOSE
  ! This routine determines a random unit vector for an outgoing vector
  ! meson in a gamma N -> V N colllision.
  ! (V is a vector meson: rho, omega, phi, J/psi)
  !
  ! INPUTS
  ! * type(particle) :: teilchenIn(1:2) --- incoming photon and nucleon
  ! * real,          :: srts            --- sqrt(s)
  ! * integer        :: id3, id4        --- ID of outgoing particles
  ! * real           :: mass3, mass4    --- mass of outgoing particle
  !
  ! OUTPUT
  ! * real, dimension(1:3) :: pscatt --- unit vector in outgoing baryon
  !   direction in CM frame
  !
  ! NOTES
  ! The given vector is the vector of particle given by id3, mass3.
  !****************************************************************************
  function pscatt_gammaN_VN (partIn, id3, id4, mass3, mass4, srts, mediumAtColl) result (pscatt)
    use particleDefinition, only: particle
    use constants, only: pi
    use random, only: rn
    use IdTable, only: photon, isMeson
    use rotation, only: rotateTo
    use VecMesWinkel, only: vecmesa
    use mediumDefinition
    use mesonPotentialModule, only: vecMes_massShift
    use minkowski, only: abs4Sq

    type(particle), dimension(1:2), intent(in) :: partIn
    integer, intent(in) :: id3,id4
    real, intent(in) :: mass3, mass4, srts
    type(medium), intent(in) :: mediumAtColl
    real, dimension(1:3) :: pscatt

    real :: cost, sint, phi, massMes
    integer :: idMes
    real, dimension(0:3) :: photonMom

    if (id3 < id4) then
       idMes = id4
       massMes = mass4
    else
       idMes = id3
       massMes = mass3
    end if

    ! take care of mass shift
    massMes = massMes + vecMes_massShift(idMes,mediumAtColl%density)

    ! get cos(Theta) and phi of photon
    if (partIn(1)%ID==photon) then
      photonMom = partIn(1)%momentum
    else
      photonMom = partIn(2)%momentum
    end if

    ! get scattering angle
    call vecmesa(srts,idMes,massMes,abs4Sq(photonMom),cost,iParam_gammaNVN)

    sint=sqrt(max(1.-cost**2,0.))
    phi=rn()*2.*pi
    pscatt(1)=cos(phi)*sint
    pscatt(2)=sin(phi)*sint
    pscatt(3)=cost

    ! Rotate according to incoming photon (?)
    pscatt = rotateTo (photonMom(1:3), pscatt)

    if (.not. isMeson(id3)) pscatt = -pscatt

  end function pscatt_gammaN_VN



  !****************************************************************************
  !****if* winkelVerteilung/pscatt_piN_piN_backward
  ! NAME
  ! function pscatt_piN_piN_backward (teilchenIn, srts) result (pscatt)
  ! PURPOSE
  ! Determines a random unit vector for an outgoing pion in a piN->piN
  ! colllision.
  !
  ! Backward peaked pion nucleon cross section (e.g. nucl-ex 0403040),
  ! valid for sqrt(s)<1.35 GeV
  !
  ! INPUTS
  ! * type(particle) :: teilchenIn(1:2) --- incoming pion and nucleon
  ! * real           :: srts            --- sqrt(s)
  ! OUTPUT
  ! * real, dimension(1:3) :: pscatt --- unit vector in outgoing pion
  ! direction in CM frame
  !****************************************************************************
  function pscatt_piN_piN_backward (partIn, srts) result (pscatt)
    use rotation, only: rotateTo
    use random, only: rn
    use constants, only: pi, mN
    use idTable, only: pion
    use particleDefinition

    type(particle), dimension(1:2), intent(in) :: partIn
    real, intent(in) :: srts
    real, dimension(1:3) :: pscatt

    real :: cost, sint, phi, x

    do
       cost=1-(rn()*2.)
       x=rn()*(1+pionNucleon_Backward_coeff)**(pionNucleon_Backward_exponent*(1.235-mN)/1.235)
       if (x<(pionNucleon_Backward_coeff-cost)**(pionNucleon_Backward_exponent*(1.235-srts)/1.235)) exit
    end do
    sint=sqrt(max(1.-cost**2,0.))
    phi=rn()*2.*pi
    pscatt(1)=cos(phi)*sint
    pscatt(2)=sin(phi)*sint
    pscatt(3)=cost
    if (debug) then
       write(77,*) cost
       write(34,*) srts, partIn(1)%Id, partIn(2)%id
    end if

    ! Rotate according to incoming pion
    if (partIn(1)%ID == pion) then
      pscatt = rotateTo (partIn(1)%momentum(1:3), pscatt)
    else
      pscatt = rotateTo (partIn(2)%momentum(1:3), pscatt)
    end if

  end function pscatt_piN_piN_backward



  !****************************************************************************
  !****f* winkelVerteilung/pscatt_BarBar
  ! NAME
  ! function pscatt_BarBar (teilchenIn, teilchenOut, srts, betaToCM, successflag) result(pscatt)
  !
  ! PURPOSE
  ! This subroutine determines the scattering angle for baryon-baryon scattering processes.
  ! Especially for NN<-> NN, NN<->NDelta, Nbar N -> LambdaBar Lambda ...
  !
  ! INPUTS
  ! * type(particle) :: teilchenIn(1:2)   --- Initial state particles
  ! * type(particle) :: teilchenOut(1:2)  --- Final state particles (only their Id's,antiflags,charges and masses are needed here)
  ! * real :: srts                        --- SQRT(s) in vacuum
  ! * real, dimension(1:3) :: betaToCM    --- Center of mass velocity in calc. frame
  !
  ! OUTPUT
  ! * real, dimension(1:3) :: pscatt --- unit vector pscatt in the direction of the outgoing teilchen3 in the CM frame
  !****************************************************************************
  function pscatt_BarBar (partIn, partOut, srts, betaToCM, successflag) result(pscatt)

    use particleDefinition
    use IdTable, only: Nucleon, Delta, Lambda, SigmaResonance, nRes
    use constants, only: pi
    use random, only: rn, rnOmega, rnExp
    use lorentzTrafo, only: lorentz
    use dimi, only: mnnnd2
    use NbarN_to_NbarDelta, only: mNbarN_to_NbarD2
    use RMF, only: getRMF_flag
    use baryonWidthMedium, only: get_MediumSwitch_coll
    use barBar_barBar, only: get_icugnon
    use winkel_tools, only: Cugnon_bnp, Cugnon_bpp, tSlope_EL_pbarp, tSlope_CEX_pbarp, dsigdt_LE, dsigdt_Regge
    use constants, only: mN

    type(particle), dimension(1:2), intent(in) :: partIn, partOut
    real, intent(in) :: srts
    real, dimension(1:3), intent(in) :: betaToCM
    logical, intent(out) :: successflag
    real, dimension(1:3) :: pscatt

    integer :: id1,id2,id3,id4,i
    real :: m1,m2,m3,m4,cost
    real, dimension(0:3) :: pcm        ! CM momentum of teilchenIn(1)
    real :: bb,a,ta,t,c1,c1_min,c1_max,asrt,as,t1,t2,s1,c2,s2,ct1,st1,ct2,st2,ss
    real :: s,plab,prcm,mdel,pi2,pf2,dsdm,u,rntheta,forback,sum
    integer, parameter :: ntheta=100
    real, dimension(-ntheta:ntheta) :: dsdodm
    logical :: flag
    logical, parameter :: debugflag = .false.

    successflag=.true.

    id1=partIn(1)%Id
    id2=partIn(2)%Id
    id3=partOut(1)%Id
    id4=partOut(2)%Id

    m1=partIn(1)%mass
    m2=partIn(2)%mass
    m3=partOut(1)%mass
    m4=partOut(2)%mass

    pcm = partIn(1)%momentum
    call lorentz(betaToCM,pcm,'winkelVerteilung') ! boost from Lab to CM
    prcm = sqrt(dot_product(pcm(1:3),pcm(1:3)))

    if (prcm==0.) then
       write(*,*) ' In winkel: c.m. momentum of incoming particles = 0'
       write(*,*) id1,m1,partIn(1)%momentum, &
                  id2,m2,partIn(2)%momentum
       stop
    end if

    pcm(1:3) = pcm(1:3)/prcm  ! Renormalize pcm to 1

    if (getRMF_flag()) then ! Redefine prcm using vacuum kinematics
        prcm = ( srts**2 + m1**2 - m2**2 )**2/(4.*srts**2) - m1**2
        prcm = sqrt( max(0.,prcm) )
    end if

    if (id1==Nucleon .and. id2==Nucleon .and. id3==Nucleon .and. id4==Nucleon) then

       ! N N -> N N

       if (NNisotropic .and. .not.(partIn(1)%antiParticle .or. partIn(2)%antiParticle)) then
         ! isotropic N N -> N N scattering
         pscatt = rnOmega()
         return
       end if

       s = (srts-m1-m2+2.*mN)**2
       plab=sqrt((s/(2.*mN))**2-s)
       ! write(*,*)'in winkel: s=', s
       ! write(*,*)'in winkel: plab=', plab

       if (partIn(1)%antiParticle.eqv.partIn(2)%antiParticle) then

          ! proper NN -> NN or Nbar Nbar -> Nbar Nbar

          select case (get_icugnon())
          case (1)
             !*****  New scattering prescription:
             if (partIn(1)%charge+partIn(2)%charge==1) then
                ! pn
                bb = max(Cugnon_bnp(plab),1.e-06)
                if (plab < 0.8) then
                   a = 1.
                else
                   a = 0.64/plab**2
                end if
             else
                ! pp or nn
                bb = max(Cugnon_bpp(plab),1.e-06)
                a = 1.
             end if
             ta = -4.*prcm**2
             if (rn() <= 1./(1.+a)) then
                t = rnExp (bb, 0., ta)
             else
                t = ta - rnExp (bb, 0., ta)
             end if
             c1 = 1. - 2.*t/ta
          case (0)
             !****   Old scatt. prescription: *****************
             asrt=srts-m1-m2
             as=(3.65*asrt)**6
             a=max(6.0*as/(1.0+as),1e-06)
             ta=-2.0*prcm**2
             !*dsdo= c*exp(a*t)=P(t)
             !*t(x) follows from P(t) dt = P(x) dx
             !*and t(-1)=2*ta,t(1)=0
             t1 = rnExp (a, 0., 2.*ta)
             c1=1.0-t1/ta
          end select

       else

          ! Nbar N -> Nbar N

          if (partOut(1)%charge==partIn(1)%charge .or. partOut(1)%charge==partIn(2)%charge) then
             bb = tSlope_EL_pbarp(srts-m1-m2+2.*mN)   ! Elastic scattering
          else
             bb = tSlope_CEX_pbarp(srts-m1-m2+2.*mN)  ! Charge exchange
          end if

          ta = -4.*prcm**2
          t = rnExp (bb , 0., ta)
          c1 = 1. - 2.*t/ta

          if (partOut(1)%antiParticle.neqv.partIn(1)%antiParticle) c1=-c1

       end if

       if (abs(c1) > 1.) c1=2.0*rn()-1.0

    else if (id1+id2+id3+id4==5) then

       ! N N <-> N Delta

       if (id1+id2==2) then
          ! N N -> N Delta
          if (Id3==Delta) then
             mdel=m3
          else
             mdel=m4
          end if
       else
          ! N Delta -> N N
          if (id1==Delta) then
             mdel=m1
          else
             mdel=m2
          end if
       end if

       s=srts**2

       pi2=(s - mN**2 + mdel**2)**2/(4.*s) - mdel**2
       pf2=s/4. - mN**2
       if (pi2*pf2<0.) then
          write(*,*) 'winkelverteilung: not possible with this kinematics'
          write(*,*) 'srts,pi2,pf2', srts,pi2,pf2
          write(*,*) 'incoming particles (masses):',id1,id2,m1,m2
          write(*,*) 'outgoing particles (masses):',id3,id4,m3,m4
          if (get_MediumSwitch_coll()) then
             successflag=.false.
             return
          else
             write(*,*) 'if get_MediumSwitch_coll=.false. -> serious error -> stop'
             stop
          end if
       end if

       if (partIn(1)%antiParticle.eqv.partIn(2)%antiParticle) then

         ! proper N N <-> N Delta or Nbar Nbar <-> Nbar Deltabar

          dsdm=0.
          do i=0,ntheta-1
             cost=(float(i)+0.5)/float(ntheta)
             t = mdel**2 + mN**2 - 2.*sqrt(mdel**2+pi2)*sqrt(mN**2+pf2) + 2.*cost*sqrt(pi2*pf2)
             u =        2.*mN**2 - 2.*sqrt(mN**2+pf2)  *sqrt(mN**2+pi2) - 2.*cost*sqrt(pi2*pf2)
             if (debugflag) write(*,*) 'before mnnnd2:', t,u,mdel
             !****   Relativistic matrix element^2:
             dsdodm(i) = mnnnd2 (t,u,mdel)
             if (debugflag) write(*,*) 'dsdodm=',dsdodm(i)
             dsdm = dsdm + dsdodm(i)/float(ntheta)
          end do

          rntheta= rn()
          forback= 2.0*(nint(rn())-0.5)

          sum= 0.
          flag= .true.
          i= -1
          do while(flag)
             i= i + 1
             if (i>ntheta-1) then
                write(*,*) 'problem in pscatt_BarBar', i
                write(*,*) rntheta, dsdodm(i-1), sum , dsdm
                stop
             end if
             sum = sum + dsdodm(i)/float(ntheta)
             if (sum>=rntheta*dsdm) then
               flag = .false.
               c1 = forback * (float(i)+0.5)/float(ntheta)
             end if
          end do

       else

          ! Nbar N <-> Nbar Delta or Nbar N <-> Deltabar N

          dsdm=0.
          do i=-ntheta,ntheta-1
             cost=(float(i)+0.5)/float(ntheta)
             t = mdel**2 + mN**2 - 2.*sqrt(mdel**2+pi2)*sqrt(mN**2+pf2) + 2.*cost*sqrt(pi2*pf2)
             if (debugflag) write(*,*) 'before mNbarN_to_NbarD2:', t,mdel
             !****   Relativistic matrix element^2:
             dsdodm(i) = mNbarN_to_NbarD2 (t,mN,mN,mN,mdel)
             if (debugflag) write(*,*) 'dsdodm=', dsdodm(i)
             dsdm= dsdm + dsdodm(i)/float(ntheta)
          end do

          rntheta= rn()

          sum= 0.
          flag= .true.
          i= -ntheta-1
          do while(flag)
             i= i + 1
             if (i>ntheta-1) then
                write(*,*) 'problem in pscatt_BarBar', i
                write(*,*) rntheta, dsdodm(i-1), sum , dsdm
                stop
             end if
             sum= sum + dsdodm(i)/float(ntheta)
             if (sum>=rntheta*dsdm) then
                flag= .false.
                c1= (float(i)+0.5) /float(ntheta)
             end if
          end do

          if (partOut(1)%antiParticle.neqv.partIn(1)%antiParticle) c1=-c1

       end if

       c1_min=c1-0.5/float(ntheta)
       c1_max=c1+0.5/float(ntheta)
       c1= c1_min + rn()*(c1_max-c1_min)
       if (c1>1.) then
         c1= c1_min + rn()*(1.-c1_min)
         if (c1>1.) then
           write(*,*) 'problem in pscatt_BarBar: c1=', c1
           stop
         end if
       else if (c1<-1.) then
         c1= c1_max + rn()*(-1.-c1_max)
         if (c1<-1.) then
           write(*,*) 'problem in pscatt_BarBar: c1=', c1
           stop
         end if
       end if

    else if (id1==Nucleon .and. id2==Nucleon) then

       if (id3==Lambda .and. id4==Lambda) then

          ! Nbar N -> LambdaBar Lambda
          if (srts <= 2.37) then
             c1 = dsigdt_LE (srts, prcm)
          else
             c1 = dsigdt_Regge (srts, prcm, 1)
          end if

       else if ((id3==Lambda .and. id4==SigmaResonance) .or. (id3==SigmaResonance.and.id4==Lambda)) then

          ! Nbar N -> LambdaBar Sigma or Nbar N -> Lambda SigmaBar
          c1 = dsigdt_Regge (srts, prcm, 2)

       else if (NN_NR_noniso .and. ((id3==Nucleon .and. id4<=nres+1) .or. (id4==Nucleon .and. id3<=nRes+1))) then

          ! N N -> N R
          if (id3>nucleon) then
            c1 = dsigdt_NRes (srts, m3)
          else
            c1 = dsigdt_NRes (srts, m4)
          end if

       else

          pscatt = rnOmega()  ! isotropic
          return

       end if

       if (partOut(1)%antiParticle.neqv.partIn(1)%antiParticle) c1=-c1

    else

       ! all other baryon-baryon collisions: isotropic scattering
       pscatt = rnOmega()
       return

    end if

    ! write(*,*)'in winkel: c1=',c1

    t1=2.0*pi*rn()

    ! write(*,*)'in winkel: t1=',t1

    !**************************************************************************
    !*     com: set the new momentum coordinates
    !**************************************************************************
    if (pcm(1)==0. .and. pcm(2)==0.) then
      t2=0.
    else
      t2=atan2(pcm(2),pcm(1))
    end if
    ! write(*,*)'in winkel: t2=', t2
    if (abs(c1) > 1.) write(*,*) 'winkel ,1 ', c1
    s1=sqrt(max(1.0-c1**2,0.))
    c2=pcm(3)
    if (abs(c2) > 1.) write(*,*) 'winkel ,2 ', c2
    s2=sqrt(max(1.0-c2**2,0.))
    ct1=cos(t1)
    st1=sin(t1)
    ct2=cos(t2)
    st2=sin(t2)
    ss=c2*s1*ct1+s2*c1
    ! write(*,*)'in winkel: s1,c2,s2,ct1,st1,ct2,st2,ss=',s1,c2,s2,ct1,st1,ct2,st2,ss
    pscatt(1) = ss*ct2 - s1*st1*st2
    pscatt(2) = ss*st2 + s1*st1*ct2
    pscatt(3) = c1*c2  - s1*s2*ct1

  end function pscatt_BarBar


  !****************************************************************************
  !****if* winkelVerteilung/dsigdt_NRes
  ! NAME
  ! function dsigdt_NRes (srts, mRes) result (cost)
  ! PURPOSE
  ! This routine implements a non-isotropic angular distribution
  ! for resonance production (N N -> N R), according to dsigma/dt = b/t**a.
  !
  ! Here b is a normalization constant, and the parameter a is taken from a fit
  ! to HADES data by A. Dybczak (private communication).
  !
  ! Reference: G. Agakishiev et al., Eur.Phys.J. A50 (2014) 82,
  !            http://inspirehep.net/record/1285514
  !
  ! INPUTS
  ! * real, intent(in) :: srts --- energy in CM frame, sqrt(s)
  ! * real, intent(in) :: mRes --- mass of outgoing resonance
  !
  ! OUTPUT
  ! * real :: cost --- cos(theta_cm), cosine of polar angle in CM frame
  !****************************************************************************
  function dsigdt_NRes (srts, mRes) result (cost)
    use constants, only: mN
    use twoBodyTools, only: pCM
    use random, only: rn, rnPower

    real, intent(in) :: srts, mRes
    real :: cost

    real :: pi, pf, tmin, tmax, a, t
    real, parameter :: p(0:3) = (/ 1.46434, 5.80311, -6.89358, 1.94302 /)

    pi = pCM (srts, mN, mN)
    pf = pCM (srts, mN, mRes)
    tmin = ((mRes**2-mN**2)/(2.*srts))**2 - (pi-pf)**2
    tmax = ((mRes**2-mN**2)/(2.*srts))**2 - (pi+pf)**2

    ! fit by Adrian
    a = p(0) + p(1)*mRes**1 + p(2)*mRes**2 + p(3)*mRes**3

    t = rnPower (-a, tmin, tmax)

    cost = (2*t-tmin-tmax) / (tmin-tmax)

    if (rn()>0.5) cost = -cost   ! symmetrize

  end function dsigdt_NRes


  !****************************************************************************
  !****if* winkelVerteilung/pscatt_Rho_pipi_Distribution
  ! NAME
  ! function pscatt_Rho_pipi_Distribution (Part) result (pscatt)
  ! PURPOSE
  ! This routine determines the outgoing unit vector in the cm-frame for
  ! the pi pi final-state
  ! out of decay of a diffractive rho0.
  !
  ! The angle should be distributed according to:
  !   sin(theta)**2+2*eps*R*cos(theta)**2
  ! in the rest frame of ther rho0, where the z-axis is defined by the
  ! negative recoil momentum of the nucleon (helicity frame).
  !
  ! INPUTS
  ! * type(particle) :: Part --- incoming rho meson
  !
  ! OUTPUT
  ! * real, dimension(1:3) :: pscatt --- Scattering vector in cm-frame
  ! NOTES
  ! If initialized via HiLepton:
  ! * the momentum of Part (=rho0) is given in the system, where the
  !   photon was going in z-direction and the target was moving according
  !   fermi motion.
  ! * betaToCM = Part%momentum/Part%momentum(0)
  !****************************************************************************
  function pscatt_Rho_pipi_Distribution (Part) result (pscatt)
    use particleDefinition
    use PIL_rhoDiff, only: PIL_rhoDiffractive_GET
!     use output, only: WriteParticle
    use random, only: rn, rnOmega
    use rotation, only: rotateTo
    use constants, only: pi
    use lorentzTrafo, only: lorentz, lorentzCalcBeta

    type(particle), intent(in) :: Part
    real, dimension(1:3) :: pscatt

    logical :: isDiffractive, r
    !integer :: iEns,iNuc
    !integer :: EventType
    !real :: nu,epsilon,Q2,Weight

    real :: epsR!, dummy
    real, dimension(0:3) :: recoil
    real, dimension(3) :: beta

    real :: cost, sint, phi, x0

    r = PIL_rhoDiffractive_GET(Part%number,isDiffractive, epsR, recoil)

    if (isDiffractive) then
       if (Debug) write(*,*) 'rho: isDiffractive=.true.'

!!$       write(*,*) 'rho: isDiffractive=.true.'
!!$       write(*,*) recoil(1:3)
!!$       write(*,*) 'epsR = ',epsR
!!$       call writeParticle(6)
!!$       call writeParticle(6,0,0,Part )

       !***********************************************************************

       ! (1) Choose cos(theta) according to f(cost)=(sint**2 + 2*epsR*cost**2)

!!$       epsR=0.0 !!!! DUMMY !!!!
!!$       epsR=1.27 !!!! DUMMY !!!!

       x0 = max(1.0,2*epsR)
       do
          cost = 1.0-(rn()*2.0)
          if (rn()*x0 .lt. 1.0+(2.0*epsR-1.0)*cost**2) exit
       end do

       sint=sqrt(max(1.-cost**2,0.))
       phi=rn()*2.*pi

       pScatt(1)=cos(phi)*sint
       pScatt(2)=sin(phi)*sint
       pScatt(3)=cost

!!$       write(*,*) pscatt

!!$       pscatt=(/0.,0.,1./)   !*** DUMMY

       ! (2) Rotate according to (-recoil)

       beta = lorentzCalcBeta (Part%momentum, "pscatt_Rho_pipi")
       call lorentz(beta,recoil,"pscatt_Rho_pipi")

       pScatt = rotateTo (-recoil(1:3), pScatt)

!!$       write(*,*) '************'
!!$       write(*,*) theta_Rot,phi_Rot
!!$       write(*,*) recoil(1:3)
!!$       write(*,*) pscatt

       !***********************************************************************

!!$       stop

    else
       if (Debug) write(*,*) 'rho: isDiffractive=.false. or noInfo'
       pscatt = rnOmega()
    end if

  end function pscatt_Rho_pipi_Distribution



end module winkelVerteilung
