!******************************************************************************
!****m* /barBar_BarBar
! NOTES
! This module includes the calculation of all
! baryon baryon -> baryon baryon cross sections.
!******************************************************************************
module barBar_BarBar
  use idTable, only: Delta, F37_1950
  implicit none
  private

  public :: sigmaBB, chooseCharge, NN_NRes, NN_ResRes, nukNuk_nukNuk, delta2Body_inMedium_treatment, get_icugnon

  public :: iifac0, iifac1


  !****************************************************************************
  !****g* barBar_BarBar/mat_NR
  ! PURPOSE
  ! Squared matrix elements M**2/16pi for N N -> N R.
  ! See http://arxiv.org/abs/1203.3557.
  ! SOURCE
  !
  real, dimension(Delta:F37_1950), save :: mat_NR =  &
       (/14., & ! Delta   =02
         70., & ! P11_1440=03
          8., & ! S11_1535=04, eta/rho
          4., & ! S11_1650=05, rho
          0., & ! S11_2090=06
          4., & ! D13_1520=07, rho
          0., & ! D13_1700=08
          0., & ! D13_2080=09
         17., & ! D15_1675=10
          0., & ! G17_2190=11, omega/rho
          0., & ! P11_1710=12
          0., & ! P11_2100=13
          4., & ! P13_1720=14, rho
          0., & ! P13_1900=15, omega/rho
          4., & ! F15_1680=16, rho
          0., & ! F15_2000=17
          0., & ! F17_1990=18
          7., & ! S31_1620=19, rho
          0., & ! S31_1900=20
          7., & ! D33_1700=21, rho
          0., & ! D33_1940=22
          0., & ! D35_1930=23
          0., & ! D35_2350=24
          0., & ! P31_1750=25
         14., & ! P31_1910=26
         14., & ! P33_1600=27
          0., & ! P33_1920=28
          0., & ! F35_1750=29
          7., & ! F35_1905=30, rho
         14. /) ! F37_1950=31
  !****************************************************************************


  !****************************************************************************
  !****g* barBar_BarBar/mat_DR
  ! PURPOSE
  ! Squared matrix elements M**2/16pi for N N -> Delta R.
  ! See http://arxiv.org/abs/1203.3557.
  ! SOURCE
  !
  real, dimension(Delta:F37_1950), save :: mat_DR = (/210., & ! Delta   =02
                                                       0., & ! P11_1440=03
                                                       60., & ! S11_1535=04, eta/rho
                                                       12., & ! S11_1650=05, rho
                                                       0., & ! S11_2090=06
                                                       12., & ! D13_1520=07, rho
                                                       0., & ! D13_1700=08
                                                       0., & ! D13_2080=09
                                                       0., & ! D15_1675=10
                                                       0., & ! G17_2190=11, omega/rho
                                                       0., & ! P11_1710=12
                                                       0., & ! P11_2100=13
                                                       12., & ! P13_1720=14, rho
                                                       0., & ! P13_1900=15, omega/rho
                                                       12., & ! F15_1680=16, rho
                                                       0., & ! F15_2000=17
                                                       0., & ! F17_1990=18
                                                       21., & ! S31_1620=19, rho
                                                       0., & ! S31_1900=20
                                                       21., & ! D33_1700=21, rho
                                                       0., & ! D33_1940=22
                                                       0., & ! D35_1930=23
                                                       0., & ! D35_2350=24
                                                       0., & ! P31_1750=25
                                                       0., & ! P31_1910=26
                                                       0., & ! P33_1600=27
                                                       0., & ! P33_1920=28
                                                       0., & ! F35_1750=29
                                                       21., & ! F35_1905=30, rho
                                                       0. /) ! F37_1950=31
  !****************************************************************************


  !****************************************************************************
  !****g* barBar_BarBar/icugnon
  ! SOURCE
  !
  integer, save :: icugnon=1
  ! PURPOSE
  ! Switch for nucleon nucleon -> nucleon nucleon cross sections:
  ! * 0=old parametrization
  ! * 1=new parametrization (Alexej Larionov, Cugnon)
  !****************************************************************************


  !****************************************************************************
  !****g* barBar_BarBar/use_ND_ND_model
  ! SOURCE
  !
  logical, save :: use_ND_ND_model=.false.
  ! PURPOSE
  ! Switch for delta nucleon -> delta nucleon cross sections:
  ! * false=old parametrization
  ! * true =one pion exchange model (Effenberger, Buss)
  !****************************************************************************


  !****************************************************************************
  !****g* barBar_BarBar/new_NR_NR
  ! SOURCE
  !
  logical, save :: new_NR_NR=.true.
  ! PURPOSE
  ! * .false.= Switch off the NR-> NR improvement (improvement= better NN<->NN fit is being used)
  ! * only for debugging or comparing
  !****************************************************************************


  !****************************************************************************
  !****g* barBar_BarBar/NR_NR_massSHIFT
  ! SOURCE
  !
  logical, save :: NR_NR_massSHIFT=.false.
  ! PURPOSE
  ! * .true.= Shift the srts in NR-> NR elastic collisions.
  !****************************************************************************


  !****************************************************************************
  !****g* barBar_BarBar/oldOset_treatment
  ! SOURCE
  !
  logical, save :: oldOset_treatment=.false.
  ! PURPOSE
  ! * .true.= Use the old treatment for the Oset Delta width: Put everything into 3-body.
  ! * only for debugging or comparing
  !****************************************************************************


  !****************************************************************************
  !****g* barBar_BarBar/etafac
  ! SOURCE
  !
  real, save :: etafac = 6.5
  ! PURPOSE
  ! Parameter for enhancement of p n -> N*(1535) N, relative to  p p -> N*(1535) N,
  ! in order to enhance eta production in pn collisions.
  ! See Calen et al., PRC 58 (1998) 2667.
  !****************************************************************************


  !****************************************************************************
  !****g* barBar_BarBar/rhofac
  ! SOURCE
  !
  real, save :: rhofac = 1.
  ! PURPOSE
  ! Parameter for enhancement of p n -> N*(1520) N, relative to  p p -> N*(1520) N,
  ! in order to enhance rho production in p n collisions.
  !****************************************************************************


  !****************************************************************************
  !****g* barBar_BarBar/neufac
  ! SOURCE
  !
  real, save :: neufac = 1.
  ! PURPOSE
  ! Parameter for enhancement of p n -> N R, relative to  p p -> N R, affecting all resonances.
  !****************************************************************************


  !****************************************************************************
  !****g* barBar_BarBar/neufac_roper
  ! SOURCE
  !
  real, save :: neufac_roper = 2.
  ! PURPOSE
  ! Parameter for enhancement of p n -> N N*(1440), relative to p p -> N N*(1440).
  ! See http://arxiv.org/abs/1203.3557.
  !****************************************************************************


  !****************************************************************************
  !****g* barBar_BarBar/deltaN_densityDependence
  ! SOURCE
  !
  logical, save :: deltaN_densityDependence = .false.
  ! PURPOSE
  ! Switch to turn on the density dependence of the NN <-> N Delta cross section.
  ! We use the density dependence from: Song/Ko, arXiv:1403.7363. The strength
  ! of the in-medium modification is controlled by the parameter alpha.
  !****************************************************************************

  !****************************************************************************
  !****g* barBar_BarBar/alpha
  ! SOURCE
  !
  real, save :: alpha = 1.2
  ! PURPOSE
  ! Parameter which controls the density dependence of the NN <-> N Delta cross section.
  ! We use the density dependence from: Song/Ko, arXiv:1403.7363.
  ! See also deltaN_densityDependence.
  !****************************************************************************

  real, parameter :: eps = 0.00001     ! Some epsilon
  logical, save   :: initFlag=.true.

contains


  !****************************************************************************
  !****s* barBar_BarBar/readInput
  ! NAME
  ! subroutine readInput
  ! PURPOSE
  ! Reads input in jobcard out of namelist "barBar_barBar".
  !****************************************************************************
  subroutine readInput

    use output
    integer :: ios, i

    !**************************************************************************
    !****n* barBar_BarBar/barBar_barBar
    ! NAME
    ! NAMELIST barBar_barBar
    ! PURPOSE
    ! Includes the switches:
    ! * mat_NR
    ! * mat_DR
    ! * icugnon
    ! * use_ND_ND_model
    ! * new_NR_NR
    ! * NR_NR_massSHIFT
    ! * oldOset_treatment
    ! * etafac
    ! * rhofac
    ! * neufac
    ! * neufac_roper
    ! * deltaN_densityDependence
    ! * alpha
    !**************************************************************************
    NAMELIST /barBar_barBar/ mat_NR, mat_DR, icugnon, use_ND_ND_model, new_NR_NR, NR_NR_massSHIFT, oldOset_treatment, &
                             etafac, rhofac, neufac, neufac_roper, deltaN_densityDependence, alpha

    call Write_ReadingInput('barBar_barBar',0)
    rewind(5)
    read(5,nml=barBar_barBar,iostat=ios)
    call Write_ReadingInput('barBar_barBar',0,ios)

    write(*,*) 'mat_NR, mat_DR:'
    do i=Delta,F37_1950
       write(*,'(i7,2f12.4)') i,mat_NR(i),mat_DR(i)
    end do

    write(*,*) 'etafac: ',etafac
    write(*,*) 'rhofac: ',rhofac
    write(*,*) 'neufac: ',neufac
    write(*,*) 'neufac_roper: ',neufac_roper
    write(*,*) '  Switch for NN -> NN=',icugnon
    write(*,*) '  USE ND -> ND model ?',use_ND_ND_model
    write(*,*) ' deltaN_densityDependence: ', deltaN_densityDependence
    write(*,*) ' alpha: ', alpha
    if (.not.new_NR_NR)    write(*,*) '  WARNING: Old NN -> NN fit is being used for NR ->NR!!'
    if ( NR_NR_massSHIFT ) write(*,*) '  INFO   : We use the mass shift in NR-> NR collisions'
    if ( oldOset_treatment) write(*,*) '  WARNING: We use the old treatment of Delta scattering          &
         &      in case that the in medium Delta width is being used! Everything goes to 3-body!'
    call Write_ReadingInput('barBar_barBar',1)

    initFlag = .false.

  end subroutine readInput


  integer function get_icugnon()
    if (initFlag) call readInput
    get_icugnon = icugnon
  end function


  !****************************************************************************
  !****f* barBar_BarBar/sigmaBB
  ! NAME
  ! real function sigmaBB(teilchenIN,idOut,mediumAt,srtfree)
  ! NOTES
  ! * This function calculates all baryon-baryon -> baryon baryon cross sections.
  ! * The cross sections are given in mB. Note that we sum over
  !   the channels with different charges of the final state particles.
  ! INPUTS
  ! * type(particle), dimension (1:2), intent(in):: teilchenIN  -- incoming particles
  ! * integer, dimension (1:2), intent(in)       :: idOut       -- Id's of outgoing particles
  ! * real, intent(in)                           :: srtfree     -- SQRT(s)
  ! * type(medium),intent(in)                    :: mediumAt    -- mediuma at the collision point
  ! OUTPUT
  ! * logical, intent(out)                       :: pauliIncluded  -- true = cross section includes Pauli blocking
  ! NOTES
  ! Included Xsections:
  !   NN -> NN
  !   NN <-> Delta Delta
  !   NR <-> NR'
  !   NN <-> NR
  ! where N is a nucleon, R and R' are  non-strange non-charmed baryons.
  !****************************************************************************
  real function sigmaBB (partIn, idOut, mediumAt, srtfree, pauliIncluded)
    use mediumDefinition
    use particleProperties, only: hadron
    use particleDefinition
    use IdTable, only: nucleon, delta, isBaryon
    use twoBodytools, only: pCM, get_PInitial
    use baryonWidthMedium, only: get_MediumSwitch_Delta
    use RMF, only: getRMF_flag
    use barBar_to_barBar_model, only: ND_ND_model
    use deltaWidth, only:deltaOset_ND_ND, deltaOset_ND_NN
    use constants, only: mN, mPi

    type(particle), dimension(1:2), intent(in)  :: partIn
    integer,        dimension(1:2), intent(in)  :: idOut
    type(medium),                   intent(in)  :: mediumAt
    real,                           intent(in)  :: srtfree
    logical,                        intent(out) :: pauliIncluded

    integer, dimension(1:2) :: isoIn, idIn, chargeIn  ! isospins, IDs and charges of incoming particles
    real,    dimension(1:2) :: massIn                 ! masses of incoming particles
    real    :: massDelta, pInitial, srtFree_withoutMass
    integer :: resID, resIDOut

    if (initFlag) call readInput

    !##########################################################################
    ! (1) Initialize variables
    !##########################################################################
    pauliIncluded = .false.
    sigmabb = 0.

    idIn(1:2)     = partIn(1:2)%ID
    chargeIn(1:2) = partIn(1:2)%charge
    massIn(1:2)   = partIn(1:2)%mass
    isoIn(1:2)    = hadron(idIn(1:2))%isoSpinTimes2

    ! Check that all particles are baryons
    if (.not.(isBaryon(idIn(1)).and.isBaryon(idIn(2)).and.isBaryon(idOut(1)).and.isBaryon(idOut(2)))) then
       write(*,*) 'sigmaBB: invalid IDs', idIn, idOut
       stop
    end if

    !##########################################################################
    ! (2) Define initial momenta of incoming particles
    !##########################################################################
    if (.not. getRMF_flag()) then
      pInitial = get_pInitial (partIn, 0)
    else
      ! Here we take a vacuum expression, since srtfree has been already
      ! corrected for the in-medium thresholds (see generateFinalState):
      pInitial = pCM (srtfree, massIn(1), massIn(2))
    end if

    !##########################################################################
    ! (3) Evaluate the Xsections. We first distinguish the processes according to the incoming particles ...
    !##########################################################################

    if ((idIN(1)==nucleon) .and. (idIN(2)==nucleon)) then

       !***********************************************************************
       ! N N in incoming channel
       !***********************************************************************
       if ((idOut(1)==nucleon) .and. (IdOut(2)==nucleon)) then
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ! elastic scattering N N <-> N N
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          sigmabb = nukNuk_nukNuk (srtfree, massIn, sum(chargeIn))
       else if (srtfree > 2*mN + mPi) then
          ! Srts is greater than the treshold for resonance production and the final state is not NN
          if (((idOut(1).eq.nucleon).and.(idOut(2).eq.delta)).or.((idOut(2).eq.nucleon).and.(idOut(1).eq.delta))) then
             !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             ! N N -> N Delta
             !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             sigmabb = NN_NDelta (srtFree,pInitial,idOut,isoIn,chargeIn,mediumAt)
          else if (idOut(1).eq.nucleon.or.idOut(2).eq.nucleon) then
             !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             !  N N -> N Res
             !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             sigmabb = NN_NRes (srtFree,pInitial,idOut,isoIn,chargeIn)
          else
             !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             !  N N -> Res Res
             !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             sigmabb = NN_ResRes (srtFree,pInitial,idOut,isoIn,chargeIn)
          end if
       end if

    else if ((((idIN(1).eq.nucleon).and.(idIN(2).eq. delta)) &
         &  .or.((idIN(2).eq.nucleon).and.(idIN(1).eq. delta))) &
         & .and.((idOut(1).eq.nucleon).and.(idOut(2).eq.nucleon))) then
       !***********************************************************************
       ! N Delta -> N N
       !***********************************************************************
       if (get_MediumSwitch_Delta()) then
          ! Delta has a dressed width : See baryonWidth_Medium!
          if (delta2Body_inMedium_treatment()) then
             if ( (sum(chargeIN).gt.2) .or. ( sum(chargeIN).lt.0) ) then
                ! No charge conservation possible!!!
                sigmabb=0.
             else
                sigmabb=deltaOset_ND_NN(partIn)
                pauliIncluded=.true.
             end if
          else
             ! We consider all scattering processes of the delta as decays of the delta into a 3 body channel.
             sigmabb=0.
          end if
      else
          if (idIN(1).eq.delta) then
             massDelta=partIn(1)%mass
          else
             massDelta=partIn(2)%mass
          end if
          sigmabb = NDelta_NN (srtFree,pInitial,idOut,isoIn,chargeIn,mediumAt,massDelta)
       end if
    else if (idIN(1).eq.nucleon.or.IdIN(2).eq.nucleon) then
       !***********************************************************************
       ! N R -> X
       !***********************************************************************

       ! Define the Id of the resonance
       if (IdIn(1)==nucleon) then
          resID=idIn(2)
       else
          resID=idIn(1)
       end if

       if (IdOut(1)==nucleon .and. IdOut(2)==nucleon) then
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ! N Res -> N N
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          sigmabb = NRes_NN (srtFree,pInitial,idOut,isoIn,chargeIn,resID)
       else if (IdOUT(1).eq.nucleon.or.IdOut(2).eq.nucleon) then
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ! N Res -> N Res'
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          if (IdOut(1).eq.nucleon) then
             resIDOut=idOut(2)
          else
             resIDOut=idOut(1)
          end if
          if (resId==Delta .and. resIDOut==Delta .and. get_MediumSwitch_Delta()) then
             ! ND -> ND where Delta has in medium changes
             if (delta2Body_inMedium_treatment()) then
                sigmabb=deltaOset_ND_ND(partIn)
                pauliIncluded=.true.
             else
                ! We consider all scattering processes of the delta as decays of the delta into a 3 body channel.
                sigmabb=0.
             end if
          else if (resId==Delta .and. resIDOut==Delta .and. use_ND_ND_model) then
             sigmabb = ND_ND_model (srtFree,idIn,chargeIn,idOut,(/999,999/),massIn,.true.)
          else
             srtFree_withoutMass=srtFree-partIn(1)%mass-partIn(2)%mass
             sigmabb = NR_NR(resID,resIDout,srtFree,pinitial,idOut, isoIn, chargeIn, srtFree_withoutMass)
          end if
       else
          write(*,*) 'N Res -> Res Res not yet implemented:',IDIN,'->',IdOUT
          stop
       end if
    else if (idIn(1)==Delta .or. idIn(2)==Delta) then
       !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       ! Delta Res -> X
       !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       if (idIn(1)==Delta) then
          resID=idIn(2)
       else
          resID=idIn(1)
       end if
       if (IdOut(1)==nucleon .and. IdOut(2)==nucleon) then
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ! Delta Res -> N N
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          sigmabb = DRes_NN (srtFree, pInitial, idOut, isoIn, chargeIn, resID)
       else
          write(*,*) 'Delta R -> N R or R R not yet implemented'
          stop
       end if
    else
       write(*,*) 'Res Res collisions not yet implemented. Incoming particles:',idIn
       stop
    end if

  end function sigmabb


  !****************************************************************************
  !****f* barBar_BarBar/m2d16pi
  ! NAME
  ! real function m2d16pi (resID, totalCharge)
  ! PURPOSE
  ! Matrix Elements for NN -> NR, see:
  !  * Teis et al., Pion-production in heavy-ion collisions at SIS energies, Z. Phys. A 356 (1997) 421-435
  !  * Weil et al., http://arxiv.org/abs/1203.3557.
  ! INPUTS
  ! * integer, intent(in) :: resID         ---   resonance ID
  ! * integer, intent(in) :: totalCharge   ---   total charge of the reaction
  ! OUTPUT
  ! * matrixelement**2 / 16 pi
  !****************************************************************************
  real function m2d16pi (resID, totalCharge)
    use IdTable, only: S11_1535, D13_1520, P11_1440
    integer, intent(in) :: resID, totalCharge

    m2d16pi = mat_NR(resID)

    if (totalCharge == 1) then
      ! modification factors for pn collisions
      select case (resID)
      case (S11_1535)
        m2d16pi = m2d16pi * max(etafac,neufac)        ! p n -> N*(1535) N
      case (D13_1520)
        m2d16pi = m2d16pi * max(rhofac,neufac)        ! p n -> N*(1520) N
      case (P11_1440)
        m2d16pi = m2d16pi * max(neufac_roper,neufac)  ! p n -> N*(1440) N
      case default
        m2d16pi = m2d16pi * neufac                    ! p n -> N R (other resonances)
      end select
    end if

  end function m2d16pi


  !****************************************************************************
  !****f* barBar_BarBar/NR_NR
  ! NAME
  ! real function NR_NR(resIN,resOut,srtFree,pinitial,idOut, isoIn, chargeIn)
  ! PURPOSE
  ! Evaluates NR -> NR scattering, where the resonances do not need to be identical.
  ! Assuming that it goes through a NN intermediate state.
  ! INPUTS
  ! * real,intent(in) :: pInitial, srtfree ! Initial momentum and sqrt(s)
  ! * integer , intent(in) :: resIn, resOut  ! Incoming and outgoing resonance ID
  ! * integer, intent(in), dimension(1:2) :: idOut,isoIn, chargeIn     ! outgoing particles, incoming isospin, incoming charges
  !****************************************************************************
  real function NR_NR (resIn, resOut, srtfree, pinitial, idOut, isoIn, chargeIn, srtFree_withoutMass)

    use twoBodyPhaseSpace, only: Integrate_2bodyPS_resonance
    use particleProperties, only: hadron
    use constants, only: mN

    integer , intent(in) :: resIn, resOut
    real,intent(in) :: srtfree, pinitial
    integer, intent(in), dimension(1:2) :: idOut, isoIn, chargeIn
    real, optional, intent(in) :: srtFree_withoutMass

    real :: csrt, mat, multipli, sigma_nuknuk, srts
    real, dimension(1:5) :: ps

    if (resOUT.eq.resIN) then
       !   same resonance in incoming and outgoing channel
       if (pinitial.gt.eps) then
          if (new_NR_NR) then
             if (.not.present(srtFree_withoutMass)) then
                write(*,*) 'srtFree_withoutMass must be given!!!'
                stop 'error in barBar_barBar/NR_NR'
             end if
             if (srtFree_withoutMass.lt.0) then
                write(*,*) 'Warning: srtFree_withoutMass.lt.0', srtFree_withoutMass, ' in barBar_barBar/NR_NR'
             end if
             if (NR_NR_massSHIFT) then
                ! Shift the srts to nucleon nucleon case:
                srts = srtFree_withoutMass + 2*mN
             else
                srts = srtFree
             end if
             sigma_nuknuk = (nukNuk_nukNuk(srts,(/mN,mN/),1) + nukNuk_nukNuk(srts,(/mN,mN/),2))/2.
             ps = Integrate_2bodyPS_resonance (resOut, srtFree, mN, 0.)
             nr_nr = sigma_nuknuk * ps(1) / pinitial * iifac0(idOut,isoIn,chargeIn)
          else
             csrt=srtfree-2.*mN
             ps = Integrate_2bodyPS_resonance (resOut, srtFree, mN, 0.)
             nr_nr = (35.0/(1.+csrt*100.0)+20.0)*ps(1)/pinitial * iifac0(idOut,isoIn,chargeIn)
          end if
       else
          nr_nr=0.
       end if
    else
       mat = m2d16pi(resIn,sum(chargeIn)) + m2d16pi(resOut,sum(chargeIn))
       multipli = 2./(2*hadron(resIN)%spin+1.)
       ps = Integrate_2bodyPS_resonance (resOut, srtFree, mN, 0.)
       if (pinitial>eps .and. srtfree>eps) then
          nr_NR = mat/(2.*pinitial*srtfree**2) * ps(1) * multipli * iifac0(idOut,isoIn,chargeIN)
       else
          nr_nr=0.
       end if
    end if
  end function NR_NR


  !****************************************************************************
  !****f* barBar_BarBar/nukNuk_nukNuk
  ! NAME
  ! real function nukNuk_nukNuk (srtfree, massIn, totalcharge)
  ! PURPOSE
  ! Calculates NN -> N N scattering cross section
  ! INPUTS
  ! * real, intent(in) :: srtfree
  ! * real, dimension(1:2), intent(in) :: massIn            ! masses of incoming particles
  ! * integer, dimension (1:2), intent(in):: totalcharge    ! total charge of incoming particles
  ! NOTES
  ! The parameter icugnon can be used to swith between two different modes of the cross section.
  ! References:
  ! * Cugnon old: J. Cugnon, T. Mizutani and J. Vandermeulen, Nucl. Phys. A 352, 505 (1981).
  ! * Cugnon new: J. Cugnon, J. Vandermeulen, D. L’hote, Nuclear Instrum. Methods B 111 (1996) 215–220.
  ! * PDG: L. Montanet et al., Phys. Rev. D 50 (1994) 1173–1814
  !****************************************************************************
  real function nukNuk_nukNuk (srtfree, massIn, totalcharge)

    use constants, only: mN

    real, intent(in) :: srtfree
    real, dimension(1:2), intent(in) :: massIn  ! masses of incoming particles
    integer, intent(in) :: totalcharge          ! total charge of incoming particles

    real :: plab, srtf

    if (srtfree>1.338) then
       if (icugnon==0) then
          ! old Cugnon parametrization as used by Effenberger
          nukNuk_nukNuk = 35.0/(1.+(srtfree-massIn(1)-massIn(2))*100.0) + 20.0
       else if (icugnon==1) then
          !*  new elastic cross sections:
          srtf = (srtfree-massIn(1)-massIn(2)+2.*mN)**2
          plab = sqrt((srtf/(2.*mN))**2-srtf)
          if (totalcharge==1) then ! pn channel
             if (plab<=0.525) then
                ! low-energy part: (A.L.)
                nukNuk_nukNuk = 17.0546*mN/(srtf-4.*mN**2) - 6.82506
             else if (plab<=0.8) then
                ! Cugnon
                nukNuk_nukNuk = 33. + 196.*abs(plab-0.95)**2.5
             else if (plab<=2.0) then
                ! Cugnon
                nukNuk_nukNuk=31./sqrt(plab)
             else if (plab<=2.776) then
                ! Cugnon
                nukNuk_nukNuk=77./(plab+1.5)
             else
                ! PDG
                nukNuk_nukNuk = 11.9 + 26.9*plab**(-1.21) + 0.169*log(plab)**2 - 1.85*log(plab)
             end if
          else ! pp, nn
             if (plab<=0.435) then
                ! low-energy part: (A.L.)
                nukNuk_nukNuk = 5.11775*mN/(srtf-4.*mN**2) + 1.67284
             else if (plab<=0.8) then
                ! Cugnon
                nukNuk_nukNuk = 23.5 + 1000.*(plab-0.7)**4
             else if (plab<=2.0) then
                ! Cugnon
                nukNuk_nukNuk = 1250./(plab+50.) - 4.*(plab-1.3)**2
             else if (plab<=2.776) then
                ! Cugnon
                nukNuk_nukNuk = 77./(plab+1.5)
             else
                ! PDG
                nukNuk_nukNuk = 11.9 + 26.9*plab**(-1.21) + 0.169*log(plab)**2 - 1.85*log(plab)
             end if
          end if
       else
          write(*,*) 'value of icugnon not defined!',icugnon
          stop
       end if
    end if
  end function nukNuk_nukNuk


  !****************************************************************************
  !****f* barBar_BarBar/NN_NRes
  ! NAME
  ! real function NN_NRes (srtFree,pInitial,idOut,isoIn,chargeIN,integerInput)
  ! PURPOSE
  ! Calculates NN -> N Res scattering cross section
  ! INPUTS
  ! * real, intent(in) :: srtFree
  ! * real, intent(in) :: pInitial                    ! CM impulse of incoming particles
  ! * integer, intent(in), dimension(1:2) :: idOut    ! IDs of outgoing particles
  ! * integer, intent(in), dimension(1:2) :: isoIN    ! isospin of incoming particles (times 2!!)
  ! * integer, intent(in), dimension(1:2) :: chargeIN ! charges of incoming particles
  ! * integer, intent(in), optional :: integerInput
  ! NOTES
  ! integerInput defines the choice of the regarded final state of the resonance:
  ! * integerInput=1     ! All final states allowed
  ! * integerInput=2     ! final state of Res is pion N
  ! * integerInput=3     ! final state of Res is eta N
  ! * integerInput=4     ! final state of Res is rho N
  ! * integerInput=5     ! final state of Res is omega N
  !
  ! Xsection is summed over all possible final state charge configurations.
  !
  !****************************************************************************
  real function NN_NRes (srtFree, pInitial, idOut, isoIn, chargeIn, integerInput)

    use twoBodyPhaseSpace, only: Integrate_2bodyPS_resonance
    use IDTable
    use constants, only: mN

    real, intent(in) :: srtFree, pInitial
    integer, intent(in), dimension(1:2) :: idOut, isoIn, chargeIn
    integer, intent(in), optional :: integerInput

    integer :: resonanceID, channel
    real,dimension(1:5) ::  ps

    if (initFlag) call readInput

    ! Check input
    if (pInitial < eps) then
       write(*,*) 'Warning in NN_NRes. pInitial=0:, ', pInitial
       NN_NREs=0.
       return
    end if

    if (present(integerInput)) then
      channel=integerInput
      if (channel.gt.5.or.channel.lt.1) then
         write(*,*) 'this integerInput is not valid for NN_NRES:' , integerInput
         stop
      end if
    else
      channel=1
    end if

    if (IdOut(1)==nucleon .and. idOut(2)<=F37_1950 .and. idOut(2)>=P11_1440) then
      resonanceID = idOut(2)
    else if (IdOut(2)==nucleon .and. idOut(1)<=F37_1950 .and. idOut(1)>=P11_1440) then
      resonanceId = idOut(1)
    else
      write(*,*) 'Error in NN_NRES. Idout is not conform with expected input:', idOut
      write(*,*) 'stop program'
      stop
    end if

    ps = Integrate_2bodyPS_resonance (resonanceID, srtFree, mN, 0.)

    NN_NRes = ps(channel)/(pinitial*srtfree**2) * m2d16pi(resonanceID,sum(chargeIn)) * iifac0(idOut,isoIn,chargeIn)

  end function NN_NRes


  !****************************************************************************
  !****f* barBar_BarBar/NRes_NN
  ! NAME
  ! real function NRes_NN (srtFree,pInitial,idOut,isoIn,chargeIN,resonanceIN)
  ! PURPOSE
  ! Calculates NRes -> NN scattering cross section
  ! INPUTS
  ! * real, intent(in) :: srtFree
  ! * real, intent(in) :: pInitial                    ! CM impulse of incoming particles
  ! * integer, intent(in), dimension(1:2) :: idOut    ! IDs of outgoing particles
  ! * integer, intent(in), dimension(1:2) :: isoIN    ! isospin of incoming particles (times 2!!)
  ! * integer, intent(in), dimension(1:2) :: chargeIN ! charges of incoming particles
  ! * integer, intent(in), optional :: resonanceIN    ! ID of incoming resonance
  ! NOTES
  ! Xsection is summed over all possible final state charge configurations.
  !****************************************************************************
  real function NRes_NN (srtFree, pInitial, idOut, isoIn, chargeIn, resonanceIn)

    use particleProperties, only: hadron
    use IDTable
    use constants, only: mN

    real, intent(in) :: srtFree, pInitial
    integer, intent(in), dimension(1:2) :: idOut, isoIn, chargeIN
    integer, intent(in) :: resonanceIn

    real :: pfinal, multipli

    if (initFlag) call readInput

    ! Check input
    if (pInitial < eps) then
       write(*,*) 'Warning in NRes_NN. pInitial=0:, ', pInitial
       NRes_NN=0.
       return
    end if

    if (.not.(resonanceIn<=F37_1950 .and. resonanceIn>=P11_1440)) then
      write(*,*) 'Error in NRes_NN. ResonanceIN is not conform with expected input:', resonanceIn
      stop
    end if

    pfinal = sqrt(max(srtfree**2/4.-mN**2,0.))
    multipli = 2./(2*hadron(resonanceIn)%spin+1.)
    ! Divide out the multiplicities due to detailed balance
    if (srtfree>eps) then
      NRes_NN = m2d16pi(resonanceIn,sum(chargeIn)) * pfinal/(pinitial*srtfree**2) * multipli * iifac0(idOut,isoIn,chargeIn)
    else
      NRes_NN=0.
    end if

   end function NRes_NN


  !****************************************************************************
  !****f* barBar_BarBar/NN_NDelta
  ! NAME
  ! real function NN_NDelta (srtFree, pInitial, idOut, isoIn, chargeIn, mediumAt)
  ! PURPOSE
  ! Calculates NN -> N Delta scattering cross section
  ! INPUTS
  ! * real, intent(in) :: srtFree
  ! * real, intent(in) :: pInitial ! CM impulse of incoming particles
  ! * integer, intent(in), dimension(1:2) :: idOut
  ! * integer, intent(in), dimension(1:2) :: isoIN ! isospin of incoming particles (times 2!!)
  ! * integer, intent(in), dimension(1:2) :: chargeIN ! charges of incoming particles
  ! * type(medium), intent(in) :: mediumAt
  ! NOTES
  ! Xsection is summed over all possible final state charge configurations.
  !****************************************************************************
  real function NN_NDelta (srtFree, pInitial, idOut, isoIn, chargeIn, mediumAt)

    use mediumDefinition
    use dimi, only: dimiIntegrated
    use constants, only: rhoNull

    real, intent(in) :: srtFree
    real, intent(in) :: pInitial ! CM impulse of incoming particles
    integer, intent(in), dimension(1:2) :: idOut
    integer, intent(in), dimension(1:2) :: isoIN ! isospin of incoming particles (times 2!!)
    integer, intent(in), dimension(1:2) :: chargeIN ! charges of incoming particles
    type(medium), intent(in) :: mediumAt

    ! factor 4/3 because dimixsec contains the xsection for p p->n d++
    NN_NDelta = dimiIntegrated(srtfree) * iifac0(idOut,isoIn,ChargeIn) * 4./3./pinitial

    if (deltaN_densityDependence .and. mediumAt%useMedium) then
      NN_NDelta = NN_NDelta * exp(-alpha*(mediumAt%densityProton+mediumAt%densityNeutron)/rhoNull)
    end if

  end function NN_NDelta


  !****************************************************************************
  !****f* barBar_BarBar/NDelta_NN
  ! NAME
  ! real function NDelta_NN (srtFree, pInitial, idOut, isoIn, chargeIn, mediumAt, massDelta)
  ! PURPOSE
  ! Calculates N Delta -> N N scattering cross section
  ! INPUTS
  ! * real, intent(in) :: srtFree
  ! * real, intent(in) :: pInitial ! CM impulse of incoming particles
  ! * integer, intent(in), dimension(1:2) :: idOut
  ! * integer, intent(in), dimension(1:2) :: isoIN ! isospin of incoming particles (times 2!!)
  ! * integer, intent(in), dimension(1:2) :: chargeIN ! charges of incoming particles
  ! * type(medium), intent(in) :: mediumAt
  ! * real, intent(in) :: massDelta ! mass of the Delta resonance
  ! NOTES
  ! Xsection is summed over all possible final state charge configurations.
  !****************************************************************************
  real function NDelta_NN (srtFree, pInitial, idOut, isoIn, chargeIn, mediumAt, massDelta)

    use mediumDefinition
    use dimi, only: dimiSigma
    use constants, only: rhoNull

    real, intent(in) :: srtFree
    real, intent(in) :: pInitial ! CM impulse of incoming particles
    integer, intent(in), dimension(1:2) :: idOut
    integer, intent(in), dimension(1:2) :: isoIN ! isospin of incoming particles (times 2!!)
    integer, intent(in), dimension(1:2) :: chargeIN ! charges of incoming particles
    type(medium), intent(in) :: mediumAt
    real, intent(in)  :: massDelta ! mass of the Delta resonance

    ! factor 8./3. because dimiSigma contains n d++ -> p p
    NDelta_NN = dimiSigma(srtfree,massDelta) * iifac0(idOut,isoIn,ChargeIn) * 8./3./pinitial

    if (deltaN_densityDependence .and. mediumAt%useMedium) then
      NDelta_NN = NDelta_NN * exp(-alpha*(mediumAt%densityProton+mediumAt%densityNeutron)/rhoNull)
    end if

  end function NDelta_NN


  !****************************************************************************
  !****f* barBar_BarBar/NN_ResRes
  ! NAME
  ! real function NN_ResRes (srtFree,pInitial,idOut,isoIn,chargeIN)
  ! PURPOSE
  ! Calculates NN -> Res Res scattering cross section, assuming constant matrix elements.
  ! INPUTS
  ! * real, intent(in) :: srtFree
  ! * real, intent(in) :: pInitial                    ! CM impulse of incoming particles
  ! * integer, intent(in), dimension(1:2) :: idOut    ! IDs of outgoing particles
  ! * integer, intent(in), dimension(1:2) :: isoIN    ! isospin of incoming particles (times 2!!)
  ! * integer, intent(in), dimension(1:2) :: chargeIN ! charges of incoming particles
  ! NOTES
  ! The cross section is summed over all possible final state charge configurations.
  ! Currently two channels are implemented: Delta-Delta and Delta-N*(1535).
  ! For Delta-Delta see Effenberger, PhD, appendix A.1.2.
  ! The Delta-N*(1535) channel was fitted to roughly match Pythia's pi-eta production XS
  ! around sqrts ~ 3 GeV, in order to describe HADES dilepton data for p+p@3.5GeV.
  !****************************************************************************
  real function NN_ResRes (srtFree, pInitial, idOut, isoIn, chargeIN)

    use twoBodyPhaseSpace, only: nnRR
    use IDTable, only: Delta
    use constants, only: mN, mPi

    real, intent(in) :: srtFree
    real, intent(in) :: pInitial ! CM impulse of incoming particles
    integer, intent(in), dimension(1:2) :: idOut
    integer, intent(in), dimension(1:2) :: isoIN ! isospin of incoming particles (times 2!!)
    integer, intent(in), dimension(1:2) :: chargeIN ! charges of incoming particles

    real :: ME, integral

    NN_ResRes = 0.

    if (srtfree < 2*mN + 2*mPi) return
    if (idOut(1)/=Delta) return

    ME = mat_DR(idOut(2))

    if (ME>0.) then
      integral = nnRR(srtFree,idOut)     ! Integrate the Res Res Phasespace
      NN_ResRes = integral/pinitial/srtfree**2 * ME * iifac0(idOut,isoIn,ChargeIn)
    end if

  end function NN_ResRes


  !****************************************************************************
  !****f* barBar_BarBar/DRes_NN
  ! NAME
  ! real function DRes_NN (srtFree, pInitial, idOut, isoIn, chargeIN, resonanceIN)
  ! PURPOSE
  ! Calculates Delta Res -> N N scattering cross section.
  ! INPUTS
  ! * real, intent(in) :: srtFree
  ! * real, intent(in) :: pInitial                    ! CM impulse of incoming particles
  ! * integer, intent(in), dimension(1:2) :: idOut    ! IDs of outgoing particles
  ! * integer, intent(in), dimension(1:2) :: isoIN    ! isospin of incoming particles (times 2!!)
  ! * integer, intent(in), dimension(1:2) :: chargeIN ! charges of incoming particles
  ! * integer, intent(in), optional :: resonanceIN    ! ID of incoming resonance
  ! NOTES
  ! Xsection is summed over all possible final state charge configurations.
  !****************************************************************************
  real function DRes_NN (srtFree, pInitial, idOut, isoIn, chargeIn, resonanceIn)

    use particleProperties, only: hadron
    use IDTable
    use constants, only: mN
    use mediumDefinition

    real, intent(in) :: srtFree, pInitial
    integer, intent(in), dimension(1:2) :: idOut, isoIn, chargeIN
    integer, intent(in) :: resonanceIn

    real :: pfinal, multipli

    if (initFlag) call readInput

    ! Check input
    if (pInitial < eps) then
       write(*,*) 'Warning in DRes_NN. pInitial=0:, ', pInitial
       DRes_NN=0.
       return
    end if

    if (.not.(resonanceIn<=F37_1950 .and. resonanceIn>=Delta)) then
      write(*,*) 'Error in DRes_NN. ResonanceIN is not conform with expected input:', resonanceIn
      stop
    end if

    pfinal = sqrt(max(srtfree**2/4.-mN**2,0.))
    multipli = 1./(2*hadron(resonanceIn)%spin+1.)   ! Divide out the multiplicities due to detailed balance
    if (srtfree>eps) then
      DRes_NN = mat_DR(resonanceIn) * pfinal/(pinitial*srtfree**2) * multipli * iifac0(idOut,isoIn,chargeIn)
    else
      DRes_NN = 0.
    end if

  end function DRes_NN


  !****************************************************************************
  !****f* barBar_BarBar/iifac0
  ! NAME
  ! real function iifac0 (IDout,is,iz)
  ! NOTES
  ! This function calculates the isospin-factors given in table 3.4 of
  ! my diploma thesis (M.E.) + factor for identical particles in the
  ! final state.
  ! Sum over final state charges!
  ! INPUTS
  ! * IDout  : IDs of final state
  ! * is     : isospins of initial state times 2
  ! * iz     : charges of initial state particles
  !****************************************************************************
  real function iifac0 (IDout,is,iz) result(iifac)

    use IdTable, only: isBaryon
    use particleProperties, only: hadron

    integer, intent(in) :: IDout(1:2),is(1:2),iz(1:2)

    integer :: isOut(1:2), izT, isT, isOutT

    izT = sum(iz)
    isT = sum(is)

    if (isBaryon(IDout(1))) then
       isOut(1) = hadron(IDout(1))%isoSpinTimes2
    else
       write(*,*) 'Error in iifac0, IDout(1) is no baryon:', IDout(1)
       stop
    end if
    if (isBaryon(IDout(2))) then
       isOut(2) = hadron(IDout(2))%isoSpinTimes2
    else
       write(*,*) 'Error in iifac0, IDout(2) is no baryon:', IDout(2)
       stop
    end if
    isOutT = sum(isOut)


    if (isT==2) then ! isoSpin 1/2 x isoSpin 1/2
      if (isOutT==2) then
        iifac=1.
      else if (isOutT==4) then
        if (izt==1) then
          iifac=0.5
        else
          iifac=1.
        end if
      else if (isOutT==6) then
        iifac=1.
      else
        write(*,*) 'iifac0 isospin factor not implemented 1', isT,isOutT
        stop
      end if
    else if (isT==4) then ! isoSpin 3/2 x isoSpin 1/2
      if (isOutT==2) then
        if ((izt==3).or.(izt==-1)) then
          iifac=0.
        else if ((izt==2).or.(izt==0)) then
          if (iz(1)==iz(2)) then
            iifac=0.25
          else
            iifac=0.75
          end if
        else
          iifac=0.5
        end if
      else if (isOutT==4) then
        iifac=1.
      else
        write(*,*) 'ERROR in iifac0/barBar_barBar'
        write(*,*) 'iifac0 isospin factor not implemented 4', isT,isOutT
        write(*,*) IDout,is,iz,izt
        stop
      end if
    else if (isT==6) then ! isoSpin 3/2 x isoSpin 3/2
      if (isOutT==2) then
        if ((izt>=3).or.(izt<=-1)) then
          iifac=0.
        else if ((izt==2).or.(izt==0)) then
          if (iz(1)==iz(2)) then
            iifac=0.4
          else
            iifac=0.3
          end if
        else if (abs(iz(1)-iz(2))==1) then
          iifac=0.3
        else
          iifac=0.7
        end if
      else
        write(*,*) 'iifac0 isospin factor not implemented 2', isT,isOutT
        stop
      end if
    else
      write(*,*) 'iifac0 isospin factor not implemented 3', isT,isOutT
      stop
    end if

    ! factor for identical particles in final state
    if (IDout(1)==IDout(2)) iifac=iifac*0.5

  end function iifac0


  !****************************************************************************
  !****f* barBar_BarBar/iifac1
  ! NAME
  ! real function iifac1 (IDout,is,iz,iz3v,iz4v)
  ! NOTES
  ! This function calculates the isospin-factors given by eq. (A.5) in
  ! Effenberger's PhD thesis + factor for identical particles in the
  ! final state.
  ! No sum over final state charges!
  ! INPUTS
  ! * IDout  : IDs of final state
  ! * is     : isospins of initial state times 2
  ! * iz     : charges of initial state particles
  ! * iz3v,iz4v  : charges of final state particles
  !****************************************************************************
  real function iifac1 (IDout,is,iz,iz3v,iz4v) result(iifac)

    use IdTable, only: isBaryon
    use particleProperties, only: hadron

    integer, intent(in) :: IDout(1:2),is(1:2),iz(1:2)
    integer, intent(in) :: iz3v, iz4v  ! Charges of final state particles

    integer :: isOut(1:2), izT, isT, isOutT

    izT = sum(iz)
    isT = sum(is)

    if (isBaryon(IDout(1))) then
       isOut(1) = hadron(IDout(1))%isoSpinTimes2
    else
       write(*,*) 'Error in iifac1, IDout(1) is no baryon:', IDout(1)
       stop
    end if
    if (isBaryon(IDout(2))) then
       isOut(2) = hadron(IDout(2))%isoSpinTimes2
    else
       write(*,*) 'Error in iifac1, IDout(2) is no baryon:', IDout(2)
       stop
    end if
    isOutT = sum (isOut)


    if (isT==2) then   ! isoSpin 1/2 x isoSpin 1/2
      if (isOutT==2) then
        if (iz3v+iz4v==1) then
          iifac=0.5
        else
          iifac=1.
        end if
      else if (isOutT==4) then
        if (izt==1) then
          iifac=0.25
        else if (iz3v==iz4v) then
          iifac=0.25
        else
          iifac=0.75
        end if
      else if (isOutT==6) then
        if (izt==1) then
          if (abs(iz3v-iz4v)==1) then
            iifac=0.15
          else
            iifac=0.35
          end if
        else if (iz3v==iz4v) then
          iifac=0.4
        else
          iifac=0.3
        end if
      else
        write(*,*) 'iifac1 isospin factor not implemented 5', isT,isOutT
        stop
      end if
    else if (isT==4) then  ! isoSpin 3/2 x isoSpin 1/2
      if (isOutT==2) then
        if ((izt==3).or.(izt==-1)) then
          iifac=0.
        else if ((izt==2).or.(izt==0)) then
          if (iz(1)==iz(2)) then
            iifac=0.25
          else
            iifac=0.75
          end if
        else
          iifac=0.25
        end if
      else if (isOutT==4) then
        if (((iz3v+iz4v)==3).or.(iz3v+iz4v)==-1) then
          iifac=1.
        else if (iz3v==iz4v) then
          if (iz(1)==iz(2)) then
            iifac=0.625
          else
            iifac=0.375
          end if
        else if (iz3v+iz4v==1) then
          iifac=0.5
        else if (iz(1)==iz(2)) then
          iifac=0.375
        else
          iifac=0.625
        end if
      else
        write(*,*) 'ERROR in iifac1/barBar_barBar'
        write(*,*) 'iifac1 isospin factor not implemented 4', isT,isOutT
        write(*,*) IDout,is,iz,izt
        stop
      end if
    else if (isT==6) then  ! isoSpin 3/2 x isoSpin 3/2
      if (isOutT==2) then
        if ((izt>=3).or.(izt<=-1)) then
          iifac=0.
        else if ((izt==2).or.(izt==0)) then
          if (iz(1)==iz(2)) then
            iifac=0.4
          else
            iifac=0.3
          end if
        else if (abs(iz(1)-iz(2))==1) then
          iifac=0.15
        else
          iifac=0.35
        end if
      else
        write(*,*) 'iifac1 isospin factor not implemented 6', isT,isOutT
        stop
      end if
    else
      write(*,*) 'iifac1 isospin factor not implemented', isT,isOutT
      stop
    end if

    ! factor for identical particles in final state
    ! if charges are different, factor is cancelled by summation
    ! (e.g. p p -> Delta0 Delta++  +  p p -> Delta++ Delta0)
    if (IDout(1)==IDout(2)) iifac=iifac*0.5

  end function iifac1


  !****************************************************************************
  !****f* barBar_BarBar/chooseCharge
  ! NAME
  ! function chooseCharge (teilchenIN, idOut) result(chargeOut)
  ! PURPOSE
  ! Two baryons with defined charges and ID's are defined as incoming particles (teilchenIN),
  ! and we also pass the outgoing particles to this routine (idOUT). This subroutine then
  ! makes a Monte-Carlo decision for the charges of the outgoing particles, using
  ! the isospin factors provided by iifac or the NDelta -> NDelta model.
  !****************************************************************************
  function chooseCharge (partIn, idOut) result(chargeOut)

    use random, only: rn
    use barbar_to_barbar_model, only: ND_ND_chooseCharge
    use particleProperties, only: hadron
    use particleDefinition
    use idTable, only: nucleon,delta,isBaryon

    type(particle), dimension(1:2), intent(in) :: partIn
    integer,        dimension(1:2), intent(in) :: idOut
    integer,        dimension(1:2)             :: chargeOut

    integer                  :: totalCharge
    integer, dimension (1:2) :: idIn        !Id's of incoming particles
    integer, dimension (1:2) :: chargeIn    !Charges of incoming particles
    integer , dimension(1:2) :: isoIN ! IsoSpin of incoming paritcles
    integer                  :: charge1, charge2
    integer                  :: charge1_Max, charge2_Max
    integer                  :: charge1_Min, charge2_Min
    real                     :: x, wahrscheinlichkeit
    integer :: minCharge, maxCharge

    !##########################################################################
    ! (1) Set important variables
    !##########################################################################
    if (initFlag) call readInput

    idIN(1:2)=partIn(1:2)%ID
    chargeIN(1:2)=partIn(1:2)%Charge
    totalCharge=Sum(chargeIn)
    ! Set isoSpin of incoming particles and checkInput
    if (isBaryon(idIn(1))) then
       isoIn(1)=hadron(idIn(1))%isoSpinTimes2
    else
       write(*,*) 'Error in barBar_barBar/chooseCharge. IdIn(1) is no baryon:', idIn(1)
       stop
    end if
    if (isBaryon(idIn(2))) then
       isoIn(2)=hadron(idIn(2))%isoSpinTimes2
    else
       write(*,*) 'Error in barBar_barBar/chooseCharge. IdIn(2) is no baryon:', idIn(2)
       stop
    end if
    ! Check that final state particles are baryons as well
    if (.not.isBaryon(idOut(1)).or..not.isBaryon(idOut(2))) then
       write(*,*) 'Outgoing particles in barBar_barBar/chooseCharge are no baryons:', idOut
       stop
    end if

    if ((Hadron(idOut(1))%strangeness.ne.0).or.(Hadron(idOut(2))%strangeness.ne.0)) then
       write(*,*) 'Outgoing particles in barBar_barBar/chooseCharge are strange baryons:', idOut
       stop
    end if

    if ((Hadron(idOut(1))%charm.ne.0).or.(Hadron(idOut(2))%charm.ne.0)) then
       write(*,*) 'Outgoing particles in barBar_barBar/chooseCharge are charmed baryons:', idOut
       stop
    end if

   !###########################################################################
   ! Choose charge of final state
   !###########################################################################
    if ((sum(idIn)==nucleon+delta) .and. (sum(idout)==nucleon+delta) .and. use_ND_ND_model) then
       chargeOut = ND_ND_chooseCharge (idIn,idOut,chargeIn)
    else
       wahrscheinlichkeit=0.
       ! x= rn() * (iifac summed over all possible final state charge configurations)
       x = rn() * iifac0(idOut,isoIn,chargeIN)

       ! Charges which are possible for the first final state particle
       charge1_min=int(float(1-hadron(idOut(1))%isoSpinTimes2)/2.)
       charge1_max=int(float(1+hadron(idOut(1))%isoSpinTimes2)/2.)

       ! Charges which are possible for the second final state particle
       charge2_min=int(float(1-hadron(idOut(2))%isoSpinTimes2)/2.)
       charge2_max=int(float(1+hadron(idOut(2))%isoSpinTimes2)/2.)

       MinCharge=max(totalCharge-charge2_max,charge1_min)    ! Minimal charge which the first final state particle can have with some given totalcharge
       MaxCharge=min(totalCharge-charge2_min,charge1_max)    ! Minimal charge which the first final state particle can have with some given totalcharge

       if (maxCharge.lt.minCharge) then
          write(*,*) ' Error in barBar_barBar/chooseCharge'
          write(*,*) ' Could not find a charge configuration: maxCharge<minCharge: '
          write(*,*) maxCharge,  mincharge
          write(*,*) idOut, isoIn, chargeIn,partIn%ID, partIn%charge
          stop
       end if

       do charge1=Mincharge,Maxcharge
          charge2=totalCharge-charge1
          wahrscheinlichkeit = wahrscheinlichkeit + iifac1(idOut,isoIn,chargeIN,charge1,charge2)
          if (wahrscheinlichkeit.ge.x) exit
          if (charge1.eq.MaxCharge) then
             write(*,*) ' Error in barBar_barBar/chooseCharge'
             write(*,*) ' Could not find a charge configuration:', idOut, isoIn, chargeIn,partIn%ID, partIn%charge
             stop
          end if
       end do
       !#######################################################################
       ! Set output
       !#######################################################################
       chargeOut(1:2) = (/charge1,charge2/)
    end if

  end function chooseCharge


  !****************************************************************************
  !****s* barBar_BarBar/delta2Body_inMedium_treatment
  ! NAME
  ! logical function delta2Body_inMedium_treatment()
  ! PURPOSE
  ! This function decides whether the condition is fulfilled such that the in-medium width of
  ! the delta can be partially implemented via two body processes.
  !****************************************************************************
  logical function delta2Body_inMedium_treatment()
    use deltaWidth, only: osetDelta_used
    use baryonWidthMedium, only: get_mediumSwitch_coll,get_MediumSwitch_Delta

    if (initFlag) call readInput

    if (oldOset_treatment) then
       ! No 2-body treatment for the Delta->Everything 3-body
       delta2Body_inMedium_treatment=.false.
    else
       delta2Body_inMedium_treatment=(get_MediumSwitch_Delta().and.(osetDelta_used().or.get_mediumSwitch_coll()))
    end if
  end function delta2Body_inMedium_treatment


end module barBar_BarBar
