!******************************************************************************
!****m*  /XsectionRatios
! NAME
! module XsectionRatios
! NOTES
! Computes and stores the ratios of the in-medium
! and vacuum BB -> BB + mesons cross sections
! (see M. Wagner et al., PRC 71, 034910 (2005).
!******************************************************************************
module XsectionRatios

  implicit none
  private

  public :: accept_event, getRatio, getShift0, getSigmaTotal, getSigmaScreened

  ! Global variables:

  integer, parameter :: nMshift = 20         ! number of the mass shift bins
  real, parameter :: MshiftBin = 0.04        ! mass shift bin (GeV)
  integer, parameter :: nSrtsStar = 100      ! number of the SrtsStar bins
  real, parameter :: SrtsStarBin = 0.1       ! SrtsStar bin (GeV)
  integer, parameter :: nCh=35               ! number of reaction channels
  real, save, allocatable, dimension(:,:,:) :: rat    ! array of cross section ratios
  real, save, allocatable, dimension(:,:) :: rat_max  ! array of maximum cross section ratios
  !(for fixed mass shift and srts)

  real, save, allocatable, dimension(:) :: rat_abs_max  ! array of absolute maximum cross section ratios
  !(for fixed mass shift only)

  real, save, allocatable, dimension(:,:) :: SigmaTotal ! total in-medium pp cross section (mbarn)


  !****************************************************************************
  !****g* XsectionRatios/flagScreen
  ! SOURCE
  !
  logical, save :: flagScreen = .false.
  ! PURPOSE
  ! If .true. -- The in-medium screening is applied to the input cross section.
  ! If .false. -- No cross section modification .
  !****************************************************************************

  !****************************************************************************
  !****g* XsectionRatios/flagInMedium
  ! SOURCE
  !
  logical, save :: flagInMedium = .false.
  ! PURPOSE
  ! If .true. -- In-medium ratios are used to decide whether an event is accepted or not.
  ! If .false. -- The event is always accepted
  !****************************************************************************


  !****************************************************************************
  !****g* XsectionRatios/flagTabulate
  ! SOURCE
  !
  logical, save :: flagTabulate = .false.
  ! PURPOSE
  ! If .true. -- in-medium ratios are tabulated.
  ! If .false. -- in-medium ratios are read-in.
  ! (This flag is important only if flagInMedium = .true.)
  !****************************************************************************


  !****************************************************************************
  !****g* XsectionRatios/shift0
  ! SOURCE
  !
  real, save :: shift0 = 0.
  ! PURPOSE
  ! Mass shift m-m^* (GeV) for using in elementary particle collision mode.
  !****************************************************************************


  logical, save :: iniFlag = .true.

contains

  !****************************************************************************
  !****f* XsectionRatios/getSigmaScreened
  ! NAME
  ! real function getSigmaScreened(pair,sigma)
  ! PURPOSE
  ! Returns the in-medium screened cross section (mbarn).
  ! See P. Daniewlewicz, NPA 673, 375 (2000); Acta. Phys. Pol. B 33, 45 (2002)
  ! INPUT:
  ! * type(particle), dimension(1:2), intent(in) :: pair         ! incoming particles
  ! * real, intent(in) :: sigma                                  ! not screened cross section
  !****************************************************************************
  real function getSigmaScreened(pair,sigma)

    use particleDefinition
    use IdTable, only: isMeson
    use densitymodule, only: densityAt
    use dichtedefinition
    use lorentzTrafo

    type(particle), dimension(1:2), intent(in) :: pair         ! incoming particles
    real, intent(in) :: sigma                                  ! not screened cross section

    real, parameter :: y = 1.2    ! coefficient at rho**(-2/3) in the maximal cross section
    real, dimension(1:3)  :: position, beta
    real, dimension(0:3)  :: j_baryon, momentumStar
    type(dichte) :: density
    real :: sigma_0

    if ( iniFlag ) then
       call init
       iniFlag = .false.
    end if

    ! Currently only baryon-baryon collisions are included:
    if (isMeson(pair(1)%ID) .or. isMeson(pair(2)%ID)) then
       getSigmaScreened = sigma
       return
    end if

    if ( .not.flagScreen ) then
       getSigmaScreened = sigma
       return
    end if

    position = (pair(1)%position+pair(2)%position)/2.
    density = densityAt(position)
    j_baryon = density%baryon
    momentumStar(0:3) = pair(1)%momentum(0:3) + pair(2)%momentum(0:3)
    beta = lorentzCalcBeta (momentumStar, 'XsectionRatios/getSigmaScreened')
    call lorentz(beta, j_baryon,'XsectionRatios/getSigmaScreened')

    if ( j_baryon(0) > 1.e-03 ) then
       sigma_0 = 10. * y / j_baryon(0)**0.666667
       getSigmaScreened = sigma_0 * tanh( sigma / sigma_0 )
    else
       getSigmaScreened = sigma
    end if

  end function getSigmaScreened


  !****************************************************************************
  !****f* XsectionRatios/getSigmaTotal
  ! NAME
  ! real function getSigmaTotal(pair)
  ! PURPOSE
  ! Returns the total in-medium cross section (mbarn).
  ! INPUT:
  ! * type(particle), dimension(1:2), intent(in) :: pair         ! incoming particles
  ! NOTES:
  ! The total in-medium pp cross section is used for all types of colliding baryons/antibaryons.
  ! Returns zero if there is at least one incoming meson.
  !****************************************************************************
  real function getSigmaTotal(pair)

    use particleDefinition
    use IdTable, only: isMeson

    type(particle), dimension(1:2), intent(in) :: pair         ! incoming particles

    real :: sqrtsStar, threshold, Q, shift
    real, dimension(0:3) :: momentumStar
    real, dimension(1:2) :: mstar
    integer :: i, iQ, ishift

    if ( iniFlag ) then
       call init
       iniFlag = .false.
    end if

    if (isMeson(pair(1)%ID) .or. isMeson(pair(2)%ID)) then
       ! Use asympotic high energy value of the pi^+ p total cross section:
       getSigmaTotal = 20.
       return
    end if

    if ( .not.flagInMedium ) then
       ! Use asympotic high energy value of the pp total cross section:
       getSigmaTotal = 40.
       return
    end if

    ! Determine the sqrtsStar bin number:

    momentumStar(0:3) = pair(1)%momentum(0:3) + pair(2)%momentum(0:3)
    sqrtsStar = momentumStar(0)**2 - Dot_Product(momentumStar(1:3),momentumStar(1:3))
    sqrtsStar = sqrt( max(0.,sqrtsStar) )
    do i=1,2
       mstar(i) = pair(i)%momentum(0)**2 - dot_product( pair(i)%momentum(1:3),&
            & pair(i)%momentum(1:3) )
       mstar(i) = sqrt( max(0.,mstar(i)) )
    end do
    threshold = mstar(1) + mstar(2)
    Q = sqrtsStar - threshold
    iQ = nint(Q/SrtsStarBin)
    iQ = max(1,min(iQ,nSrtsStar))

    ! Determine the mass shift bin number:

    shift = pair(1)%mass - mstar(1)
    ishift = nint(shift/MshiftBin)
    ishift = max(0,min(ishift,nMshift))

    ! Use total in-medium pp cross section:

    getSigmaTotal = SigmaTotal(ishift,iQ)

  end function getSigmaTotal

  !****************************************************************************
  !****************************************************************************
  real function getShift0()
    if ( iniFlag ) then
       call init
       iniFlag = .false.
    end if
    getShift0 = shift0
  end function getShift0


  !****************************************************************************
  !****f* XsectionRatios/getRatio
  ! NAME
  ! real function getRatio(pair)
  ! PURPOSE
  ! Computes the ratio of the in-medium cross section to the vacuum cross section.
  ! INPUT:
  ! * type(particle), dimension(1:2), intent(in) :: pair         ! incoming particles
  !****************************************************************************
  real function getRatio(pair)
    use particleDefinition
    use IdTable, only: isMeson
    use particleProperties, only: hadron
    use inputGeneral, only: eventtype
    use eventtypes, only: elementary

    type(particle), dimension(1:2), intent(in) :: pair         ! incoming particles

    real :: sqrtsStar, threshold, Q, shift
    real, dimension(0:3) :: momentumStar
    real, dimension(1:2) :: mstar
    integer :: i, iQ, ishift

    if ( iniFlag ) then
       call init
       iniFlag = .false.
    end if

    if ( .not.flagInMedium ) then
       getRatio = 1.
       return
    end if

    if ( isMeson(pair(1)%Id) .or. isMeson(pair(2)%Id) ) then

       ! No medium modifications implemented for meson-baryon and meson-meson collisions:
       getRatio = 1.
       return

    else if ( pair(1)%antiparticle .and. pair(2)%antiparticle ) then

       ! No medium modifications implemented for antibaryon-antibaryon collisions:
       getRatio = 1.
       return

    else if (      hadron(pair(1)%Id)%strangeness.ne.0 &
         &  .or. hadron(pair(1)%Id)%charm.ne.0 &
         &  .or. hadron(pair(2)%Id)%strangeness.ne.0 &
         &  .or. hadron(pair(2)%Id)%charm.ne.0  ) then

       !  No medium modifications implemented for collisions involving
       !  a strange or a charmed particle:
       getRatio = 1.
       return

    end if

    ! Incoming baryons are nonstrange and noncharmed

    if (pair(1)%antiparticle.neqv.pair(2)%antiparticle) then
       ! Baryon-antibaryon collision:
       getRatio = 1.
       return
    end if

    ! Determine the sqrtsStar bin number:

    momentumStar(0:3) = pair(1)%momentum(0:3) + pair(2)%momentum(0:3)
    sqrtsStar = momentumStar(0)**2 - Dot_Product(momentumStar(1:3),momentumStar(1:3))
    sqrtsStar = sqrt( max(0.,sqrtsStar) )
    do i=1,2
       mstar(i) = pair(i)%momentum(0)**2 - dot_product( pair(i)%momentum(1:3),&
            & pair(i)%momentum(1:3) )
       mstar(i) = sqrt( max(0.,mstar(i)) )
    end do
    threshold = mstar(1) + mstar(2)
    Q = sqrtsStar - threshold
    iQ = nint(Q/SrtsStarBin)
    iQ = max(1,min(iQ,nSrtsStar))

    ! Determine the mass shift bin number:

    if ( .not. eventtype==elementary ) then

       shift = pair(1)%mass - mstar(1)

    else !***** Test in elementary collisions:

       shift = shift0

    end if

    ishift = nint(shift/MshiftBin)
    ishift = max(0,min(ishift,nMshift))
    if ( ishift == 0 ) then
       ! Low density, no cross section modifications:
       getRatio = 1.
       return
    end if

    ! Use the ratio of elastic cross sections:
    ! getRatio = rat(ishift,iQ,35)

    ! Use maximum ratio:
    !getRatio = rat_max(ishift,iQ)

!    getRatio = rat_abs_max(ishift)
   getRatio=1.

  end function getRatio


  !****************************************************************************
  !****f* XsectionRatios/accept_event
  ! NAME
  ! logical function accept_event(pair,finalState)
  ! PURPOSE
  ! Accepts or rejects the event randomly, using the in-medium cross section ratios
  ! * returns .true.  --- if the event is accepted
  ! * returns .false. --- if the event is rejected
  ! INPUT:
  ! * type(particle), dimension(1:2), intent(in) :: pair         ! incoming particles
  ! * type(particle), dimension(:), intent(in) :: finalState     ! outgoing particles
  ! * real, optional :: sqrts_ini      ! initial sqrts of colliding pair
  !                                    ! (before modifications by a 3-body collision)
  !****************************************************************************
  logical function accept_event(pair,finalState,sqrts_ini)

    use IdTable
    use particleDefinition
    use particleProperties, only: hadron
    use random
    use inputGeneral, only: eventtype
    use eventtypes, only: elementary

    type(particle), dimension(1:2), intent(in) :: pair         ! incoming particles
    type(particle), dimension(:), intent(in) :: finalState     ! outgoing particles
    real, optional :: sqrts_ini      ! initial sqrts of colliding pair
    ! (before modifications by a 3-body collision)

    real :: sqrtsStar, threshold, Q, Q_ini, shift
    real :: ratio_ini, ratio1, factor
    real, dimension(0:3) :: momentumStar
    real, dimension(1:2) :: mstar

    integer :: nFinal, nBaryon, nMeson, nPion, nStrangeMeson, k, iQ, iQ_ini, ishift, i

    integer, parameter :: nBaryonFinal=2, nMesonFinal=20

    integer, dimension(1:nBaryonFinal) :: idBaryon    ! Id's of outgoing baryons
    integer, dimension(1:nMesonFinal) :: idMeson     ! Id's of outgoing mesons

    real, dimension(1:nBaryonFinal+nMesonFinal) :: mstar_final

    if ( iniFlag ) then
       call init
       iniFlag = .false.
    end if

    if ( .not.flagInMedium ) then
       accept_event = .true.
       return
    end if

    ! Initialize outgoing channel:
    k = 0

    if ( isMeson(pair(1)%Id) .or. isMeson(pair(2)%Id) ) then

       ! No medium modifications implemented for meson-baryon and meson-meson collisions:
       accept_event = .true.
       return

    else if ( pair(1)%antiparticle .and. pair(2)%antiparticle ) then

       ! No medium modifications implemented for antibaryon-antibaryon collisions:
       accept_event = .true.
       return

    else if (      hadron(pair(1)%Id)%strangeness.ne.0 &
         &  .or. hadron(pair(1)%Id)%charm.ne.0 &
         &  .or. hadron(pair(2)%Id)%strangeness.ne.0 &
         &  .or. hadron(pair(2)%Id)%charm.ne.0  ) then

       !  No medium modifications implemented for collisions involving
       !  a strange or a charmed particle:
       accept_event = .true.
       return

    end if

    ! Incoming baryons are nonstrange and noncharmed

    ! Analyse the outgoing particles:

    nFinal = 0
    nBaryon = 0
    nMeson = 0
    nPion = 0
    nStrangeMeson = 0
    idBaryon = 0
    idMeson = 0
    do i=1,size(finalState,dim=1)

       if (finalState(i)%id <= 0) exit

       nFinal = nFinal + 1

       if (nFinal > nBaryonFinal + nMesonFinal) then
          !  No medium modifications implemented:
          accept_event = .true.
          return
       end if

       if ( isBaryon(finalState(i)%id) ) then

          if ( finalState(i)%antiparticle ) then
             !  No medium modifications implemented for baryon-antibaryon pair production:
             accept_event = .true.
             return
          end if

          nBaryon = nBaryon + 1
          if ( nBaryon <= nBaryonFinal ) then
             idBaryon(nBaryon) = finalState(i)%Id
          else
             !  No medium modifications implemented for more than 2 outgoing baryons
             !  (this also must be baryon-antibaryon pair production):
             accept_event = .true.
             return
          end if

       else if ( isMeson(finalState(i)%id) ) then

          nMeson = nMeson + 1
          if ( nMeson <= nMesonFinal ) then
             idMeson(nMeson) = finalState(i)%Id
          else
             !  No medium modifications implemented:
             accept_event = .true.
             return
          end if

          if ( finalState(i)%Id .eq. pion ) nPion = nPion + 1

          if ( hadron(finalState(i)%Id)%strangeness .ne. 0 ) nStrangeMeson = nStrangeMeson + 1

       end if

    end do


    if ((pair(1)%antiparticle.neqv.pair(2)%antiparticle) .and. nBaryon.eq.0) then

       ! Baryon-antibaryon annihilation:

       if ( rn() .lt.  Ratio_BaB(pair,finalState) ) then
          accept_event = .true.
       else
          accept_event = .false.
       end if

       return

    end if


    if (nBaryon.eq.2) then
       if (        hadron(idBaryon(1))%charm .ne. 0 &
            & .or.   hadron(idBaryon(2))%charm .ne. 0 &
            & .or.   hadron(idBaryon(1))%strangeness &
            & + hadron(idBaryon(2))%strangeness.le.-2  ) then
          ! No medium modifications for outgoing charmed baryons or
          ! outgoing pair of strange baryons, or a baryon with S=-2
          accept_event = .true.
          return
       end if
    else
       accept_event = .true.
       return
    end if

    ! Choose outgoing channel:

    if ( hadron(idBaryon(1))%strangeness == 0 .and. &
         &hadron(idBaryon(2))%strangeness == 0 ) then

       ! 2 outgoing nonstrange baryons:
       select case (nMeson)
       case (0)  ! B B -> B B :

          if (      (idBaryon(1) == pair(1)%Id .and. idBaryon(2) == pair(2)%Id) &
               & .or. (idBaryon(2) == pair(1)%Id .and. idBaryon(1) == pair(2)%Id) ) then
             k = 35   ! Elastic scattering
          else
             k = 31   ! Resonance production/absorption
          end if

       case (1) ! B B -> B B M :

          if ( idMeson(1) == pion ) then

             if ( max(idBaryon(1),idBaryon(2)) == nucleon ) then
                k = 1    ! N + N + pi production
             else
                k = 34   ! B + B + pi production  (at least one of final baryons is not a nucleon)
             end if
          else
             k = 5    ! heavy meson production

          end if

       case (2) ! B B -> B B M M :

          if ( nPion == 2 ) then
             k = 2    ! 2 pi production
          else if ( nPion == 1 ) then
             k = 6    ! heavy meson + pi production
          else if ( nStrangeMeson == 0 ) then
             k = 9    ! 2 heavy nonstrange meson production
          else if ( max(idMeson(1),idMeson(2)).le.kaonBar ) then
             k = 25   ! K Kbar production
          else
             k = 33   ! Kstar Kbar or K KstarBar production
          end if

       case (3) ! B B -> B B M M M:

          if ( nPion == 3 ) then
             k = 3    ! 3 pi production
          else if ( nPion == 2 ) then
             k = 7    ! pi pi + heavy nonstrange meson production
          else if ( nPion == 1 ) then
             if ( nStrangeMeson == 0 ) then
                k = 10   ! pi + 2 heavy nonstrange meson production
             else
                k = 26   ! pi + 2 strange meson production
             end if
          else
             if ( nStrangeMeson == 0 ) then
                k = 12   ! 3 heavy nonstrange meson production
             else
                k = 28   ! heavy nonstrange meson + 2 strange meson production
             end if
          end if

       case (4) ! B B -> B B M M M M:

          if ( nPion == 4 ) then
             k = 4    ! 4 pi production
          else if ( nPion == 3 ) then
             k = 8    ! 3 pi + heavy nonstrange meson production
          else if ( nPion == 2 ) then
             if ( nStrangeMeson == 0 ) then
                k = 11 ! 2 pi + 2 heavy nonstrange meson production
             else
                k = 27 ! 2 pi + 2 strange meson production
             end if
          else if ( nPion == 1 ) then
             if ( nStrangeMeson == 0 ) then
                k = 13 ! pi + 3 heavy nonstrange meson production
             else
                k = 29 ! pi + heavy nonstrange meson + 2 strange meson productio
             end if
          else
             if ( nStrangeMeson == 0 ) then
                k = 14 ! 4 heavy nonstrange meson production
             else if ( nStrangeMeson == 2 ) then
                k = 30 ! 2 strange meson + 2 heavy nonstrange meson production
             else
                ! No medium modifications for 4 outgoing strange mesons yet:
                accept_event = .true.
                return
             end if
          end if

       case default ! B B -> B B M M M M M .... (more than 4 outgoing mesons)

          k = 35

       end select

    else

       ! outgoing nonstrange baryon + hyperon:

       if ( nMeson == 0 ) then
          write(*,*) ' In accept_event: wrong reaction B B -> B Y', &
               &pair(:)%Id, finalState(1:nFinal)%Id
          stop
       else if ( nMeson == 1 ) then
          if ( nStrangeMeson == 0 ) then
             write(*,*) ' In accept_event: wrong reaction B B -> B Y + nonstr meson', &
                  &pair(:)%Id, finalState(1:nFinal)%Id
             stop
          else
             if ( idMeson(1) == kaon ) then
                k = 15     ! kaon production
             else
                k = 32     ! Kstar production
             end if
          end if
       else if ( nMeson == 2 ) then
          if ( nPion == 2 ) then
             write(*,*) ' In accept_event: wrong reaction B B -> B Y pi pi', &
                  &pair(:)%Id, finalState(1:nFinal)%Id
             stop
          else if ( nPion == 1 ) then
             k = 16     ! strange meson + pion production
          else
             k = 19     ! strange meson + heavy nonstrange meson production
          end if
       else if ( nMeson == 3 ) then
          if ( nPion == 2 ) then
             k = 17     ! strange meson + 2 pi production
          else if ( nPion == 1 ) then
             k = 20     ! strange meson + pi + heavy nonstrange meson production
          else
             if ( nStrangeMeson == 1 ) then
                k = 22   ! strange meson + 2 heavy nonstrange meson production
             else if ( nStrangeMeson == 3 ) then
                ! No medium modifications for B B -> B Y K K Kbar yet:
                accept_event = .true.
                return
             end if
          end if
       else if ( nMeson == 4 ) then
          if ( nPion == 3 ) then
             k = 18     ! strange meson + 3 pi production
          else if ( nPion == 2 ) then
             k = 21     ! strange meson + heavy nonstrange meson + 2 pi production
          else if ( nPion == 1 ) then
             if ( nStrangeMeson == 1 ) then
                k = 23   ! strange meson + 2 heavy nonstrange mesons + pi production
             else if ( nStrangeMeson == 3 ) then
                ! No medium modifications for B B -> B Y K K Kbar pi yet:
                accept_event = .true.
                return
             end if
          else
             if ( nStrangeMeson == 1 ) then
                k = 24   ! strange meson + 3 heavy nonstrange meson production
             else if ( nStrangeMeson == 3 ) then
                ! No medium modifications for B B -> B Y K K Kbar + heavy nonstrange meson yet:
                accept_event = .true.
                return
             end if
          end if
       else
          ! B B -> B Y M M M M M .... (more than 4 outgoing mesons)
          k = 35
       end if
    end if


    if ( k .le. 0 ) then

       write(*,*) ' In accept_event: reaction channel is not set', k, &
            &pair(:)%Id, finalState(1:nFinal)%Id
       stop

    else if ( k .gt. nCh ) then

       write(*,*) ' In accept_event: wrong reaction channel', k, &
            &pair(:)%Id, finalState(1:nFinal)%Id
       stop

    end if

    ! Determine the sqrtsStar bin numbers

    momentumStar(0:3) = pair(1)%momentum(0:3) + pair(2)%momentum(0:3)
    sqrtsStar = momentumStar(0)**2 - Dot_Product(momentumStar(1:3),momentumStar(1:3))
    sqrtsStar = sqrt( max(0.,sqrtsStar) )

    do i = 1,2
       mstar(i) = pair(i)%momentum(0)**2 - dot_product( pair(i)%momentum(1:3),&
            & pair(i)%momentum(1:3) )
       mstar(i) = sqrt( max(0.,mstar(i)) )
    end do
    threshold = mstar(1) + mstar(2)

    if (present(sqrts_ini)) then
       Q_ini = sqrts_ini - threshold
    else
       Q_ini = sqrtsStar - threshold
    end if
    iQ_ini = nint(Q_ini/SrtsStarBin)
    iQ_ini = max(1,min(iQ_ini,nSrtsStar))

    Q = sqrtsStar - threshold
    iQ = nint(Q/SrtsStarBin)
    iQ = max(1,min(iQ,nSrtsStar))

    ! Determine the mass shift bin number:

    if ( .not. eventtype==elementary ) then

       shift = pair(1)%mass - mstar(1)

    else !***** Test in elementary collisions:

       shift = shift0

    end if

    ishift = nint(shift/MshiftBin)
    ishift = max(0,min(ishift,nMshift))
    if ( ishift == 0 ) then
       ! Low density, no cross section modifications:
       accept_event = .true.
       return
    end if


    ! Accept the event with a probability equal to the factor
    ! sigma_med / sigma_vac / rat_max :

    if (pair(1)%antiparticle.neqv.pair(2)%antiparticle) then
       ratio_ini = 1.
    else
       !ratio_ini = rat_max(ishift,iQ_ini)
       ! ratio_ini = rat_abs_max(ishift)
       ratio_ini = 1.
    end if

    ! write(*,*) ' In accept_event: ratio_ini=', ratio_ini

    if ( ratio_ini <= 0 ) then
       write(*,*) ' In accept_event: ratio_ini=', ratio_ini, ishift, iQ_ini
       stop
    end if

!    ratio1 = rat(ishift,iQ,k)

    do i = 1,nFinal
       mstar_final(i)=finalState(i)%momentum(0)**2 &
                      & - dot_product( finalState(i)%momentum(1:3),&
                      &                finalState(i)%momentum(1:3) )
       mstar_final(i) = sqrt( max(0.,mstar_final(i)) )
    end do

    if (2.le.nFinal .and. nFinal.le.6) then
        ratio1=Ratio(sqrtsStar,(/mstar(1:2),mstar_final(1:nFinal)/),&
                    &(/pair(1:2)%mass,finalState(1:nFinal)%mass/))
    else
        accept_event = .true.
        return
    end if

    factor = ratio1 / ratio_ini

    if ( factor > 1 ) then
       write(*,*) ' Warning: in accept_event factor=', factor, ' k=', k
       write(*,*) pair(:)%Id, finalState(1:nFinal)%Id
    end if

    if ( rn() .le. factor ) then

       accept_event = .true.

    else

       accept_event = .false.

    end if

  end function accept_event


  !****************************************************************************
  !****s* XsectionRatios/init
  ! NAME
  ! subroutine init
  ! PURPOSE
  ! Reads in input switches from namelist XsectionRatios_input.
  !****************************************************************************
  subroutine init

    use output
    use inputGeneral, only: path_To_Input

    integer :: ios, i, j, k
    real :: shift, Q
    character(8) :: dummy1, dummy2

    !**************************************************************************
    !****n* XsectionRatios/XsectionRatios_input
    ! NAME
    ! NAMELIST /XsectionRatios_input/
    ! PURPOSE
    ! Includes the switches:
    ! * flagScreen
    ! * flagInMedium
    ! * flagTabulate
    ! * shift0
    !**************************************************************************
    NAMELIST /XsectionRatios_input/ flagScreen, flagInMedium, flagTabulate, shift0

    call Write_ReadingInput('XsectionRatios_input',0)
    rewind(5)
    read(5,nml=XsectionRatios_input,iostat=ios)
    call Write_ReadingInput('XsectionRatios_input',0,ios)
    write(*,*) ' Set flagScreen to', flagScreen , '.'
    write(*,*) ' Set flagInMedium to', flagInMedium ,'.'
    write(*,*) ' Set flagTabulate to', flagTabulate ,'.'
    write(*,*) ' Set shift0 to', shift0 ,'.'
    call Write_ReadingInput('XsectionRatios_input',1)

    if ( flagScreen .and. flagInMedium ) then
       write(*,*) ' In-medium screening and in-medium reduction due to Dirac masses '
       write(*,*) ' can not be done simultaneously !!!'
       stop
    end if

    if (flagInMedium) then

       allocate( rat(1:nMshift,1:nSrtsStar,1:nCh) )

       allocate( rat_max(1:nMshift,1:nSrtsStar) )

       allocate( rat_abs_max(1:nMshift) )

       allocate( SigmaTotal(0:nMshift,1:nSrtsStar) )

       if (flagTabulate) then

          write(*,*) ' Cross section ratio tabulation starts ...'
          call TabulateRatio
          write(*,*) ' Cross section ratio tabulation is finished'
          open(1,file='XsectionRatios.dat',status='unknown')
          do i = 1,nMshift
             shift = MshiftBin * float(i)
             write(1,*)'# shift:', shift
             do j = 1,nSrtsStar
                Q = SrtsStarBin * float(j)
                write(1,'(50(e13.6,1x))') Q, rat(i,j,:)
             end do
          end do
          close(1)

       else

          write(*,*) ' Cross section ratio read-in starts ...'
          ios=0
          open(1,file=trim(path_to_input)//'/XsectionRatios.dat',status='old',iostat=ios)
          do i = 1,nMshift
             read(1,*) dummy1, dummy2, shift
             do j = 1,nSrtsStar
                read(1,*) Q, rat(i,j,:)
             end do
          end do
          close(1)

          if (ios.eq.0) then
             write(*,*) ' Cross section ratio read-in is successfully finished'
          else
             write(*,*) ' Error during cross section ratio read-in, ios= ', ios
             stop
          end if

          !            open(1,file='XsectionRatios_chk.dat',status='unknown')
          !            do i = 1,nMshift
          !               shift = MshiftBin * float(i)
          !               write(1,*)'# shift:', shift
          !               do j = 1,nSrtsStar
          !                 Q = SrtsStarBin * float(j)
          !                 write(1,'(50(e13.6,1x))') Q, rat(i,j,:)
          !               end do
          !            end do
          !            close(1)

       end if

       ! Compute maximum ratios:

       do i = 1,nMshift
          rat_abs_max(i) = 0.
          do j = 1,nSrtsStar
             rat_max(i,j) = 0.
             do k = 1,nCh
                if ( rat(i,j,k) >  rat_max(i,j) ) rat_max(i,j) = rat(i,j,k)
                if ( rat(i,j,k) >  rat_abs_max(i) ) rat_abs_max(i) = rat(i,j,k)
             end do
          end do
       end do

       !         open(1,file='rat_max.dat',status='unknown')
       !         do i = 1,nMshift
       !           shift = MshiftBin * float(i)
       !           write(1,*)'# shift:', shift
       !           do j = 1,nSrtsStar
       !             Q = SrtsStarBin * float(j)
       !             write(1,'(3(e13.6,1x))') Q, rat_max(i,j), rat_abs_max(i)
       !           end do
       !         end do
       !         close(1)

       ! Read-in the total in-medium pp cross section:
       write(*,*) ' Total in-medium pp cross section read-in starts ...'
       ios=0
       open(1,file=trim(path_to_input)//'/XsectionTotal.dat',status='old',iostat=ios)
       do i = 0,nMshift
          read(1,*) dummy1, dummy2, shift
          do j = 1,nSrtsStar
             read(1,*) Q, SigmaTotal(i,j)
          end do
       end do
       close(1)
       write(*,*) ' Total in-medium pp cross section read-in is successfully finished'

       !         open(1,file='XsectionTotal_chk.dat',status='unknown')
       !         do i = 0,nMshift
       !           shift = MshiftBin * float(i)
       !           write(1,*)'# shift:', shift
       !           do j = 1,nSrtsStar
       !             Q = SrtsStarBin * float(j)
       !             write(1,'(2(e13.6,1x))') Q, SigmaTotal(i,j)
       !           end do
       !         end do
       !         close(1)

       !         open(1,file='ratio_chk.dat',status='unknown')
       !         do i = 10,10
       !           shift = MshiftBin * float(i)
       !           write(1,*)'# shift:', shift
       !           do j = 1,nSrtsStar
       !             Q = SrtsStarBin * float(j)
       !             write(1,'(2(e13.6,1x))') Q, SigmaTotal(i,j)/SigmaTotal(0,j)
       !           end do
       !         end do
       !         close(1)

    end if

  end subroutine init


  !****************************************************************************
  !****s* XsectionRatios/TabulateRatio
  ! NAME
  ! subroutine TabulateRatio
  ! PURPOSE
  ! Tabulates the ratio sigma_med / sigma_vac as a function
  ! of srtsStar and mass shift (m-mStar).
  !****************************************************************************
  subroutine TabulateRatio

    use IdTable, only: rho, kaonStar, Delta, Lambda
    use particleProperties, only: hadron
    use constants, only: mN, mPi, mK

    integer :: i, j, k
    real :: shift, Q, SrtsStar
    real, dimension(1:8) :: mStar, mass

    Loop_over_massShift : do i = 1,nMshift

       shift = MshiftBin * float(i)

       Loop_over_SrtsStar : do j = 1,nSrtsStar

          Q = SrtsStarBin * float(j)
          SrtsStar = Q + 2. * ( mN - shift )
          write(*,*) ' i:', i, 'j:', j

          Loop_over_channels : do k = 1,nCh

             if ( k == 1 ) then !***** N_1 N_2 -> N_3 N_4 pi_5
                mass(1:4) = mN
                mStar(1:4) = mass(1:4) - shift
                mass(5) = mPi
                mStar(5) = mass(5)
                rat(i,j,k) = Ratio(SrtsStar,mStar(1:5),mass(1:5))
             else if ( k == 2 ) then !***** N_1 N_2 -> N_3 N_4 pi_5 pi_6
                mass(1:4) = mN
                mStar(1:4) = mass(1:4) - shift
                mass(5:6) = mPi
                mStar(5:6) = mass(5:6)
                rat(i,j,k) = Ratio(SrtsStar,mStar(1:6),mass(1:6))
             else if ( k == 3 ) then !***** N_1 N_2 -> N_3 N_4 pi_5 pi_6 pi_7
                mass(1:4) = mN
                mStar(1:4) = mass(1:4) - shift
                mass(5:7) = mPi
                mStar(5:7) = mass(5:7)
                rat(i,j,k) = Ratio(SrtsStar,mStar(1:7),mass(1:7))
             else if ( k == 4 ) then !***** N_1 N_2 -> N_3 N_4 pi_5 pi_6 pi_7 pi_8
                mass(1:4) = mN
                mStar(1:4) = mass(1:4) - shift
                mass(5:8) = mPi
                mStar(5:8) = mass(5:8)
                rat(i,j,k) = Ratio(SrtsStar,mStar(1:8),mass(1:8))
             else if ( k == 5 ) then !***** N_1 N_2 -> N_3 N_4 rho_5
                mass(1:4) = mN
                mStar(1:4) = mass(1:4) - shift
                mass(5) = hadron(rho)%mass
                mStar(5) = mass(5)
                rat(i,j,k) = Ratio(SrtsStar,mStar(1:5),mass(1:5))
             else if ( k == 6 ) then !***** N_1 N_2 -> N_3 N_4 rho_5 pi_6
                mass(1:4) = mN
                mStar(1:4) = mass(1:4) - shift
                mass(5) = hadron(rho)%mass
                mass(6) = mPi
                mStar(5:6) = mass(5:6)
                rat(i,j,k) = Ratio(SrtsStar,mStar(1:6),mass(1:6))
             else if ( k == 7 ) then !***** N_1 N_2 -> N_3 N_4 rho_5 pi_6 pi_7
                mass(1:4) = mN
                mStar(1:4) = mass(1:4) - shift
                mass(5) = hadron(rho)%mass
                mass(6:7) = mPi
                mStar(5:7) = mass(5:7)
                rat(i,j,k) = Ratio(SrtsStar,mStar(1:7),mass(1:7))
             else if ( k == 8 ) then !***** N_1 N_2 -> N_3 N_4 rho_5 pi_6 pi_7 pi_8
                mass(1:4) = mN
                mStar(1:4) = mass(1:4) - shift
                mass(5) = hadron(rho)%mass
                mass(6:8) = mPi
                mStar(5:8) = mass(5:8)
                rat(i,j,k) = Ratio(SrtsStar,mStar(1:8),mass(1:8))
             else if ( k == 9 ) then !***** N_1 N_2 -> N_3 N_4 rho_5 rho_6
                mass(1:4) = mN
                mStar(1:4) = mass(1:4) - shift
                mass(5:6) = hadron(rho)%mass
                mStar(5:6) = mass(5:6)
                rat(i,j,k) = Ratio(SrtsStar,mStar(1:6),mass(1:6))
             else if ( k == 10 ) then !***** N_1 N_2 -> N_3 N_4 rho_5 rho_6 pi_7
                mass(1:4) = mN
                mStar(1:4) = mass(1:4) - shift
                mass(5:6) = hadron(rho)%mass
                mass(7) = mPi
                mStar(5:7) = mass(5:7)
                rat(i,j,k) = Ratio(SrtsStar,mStar(1:7),mass(1:7))
             else if ( k == 11 ) then !***** N_1 N_2 -> N_3 N_4 rho_5 rho_6 pi_7 pi_8
                mass(1:4) = mN
                mStar(1:4) = mass(1:4) - shift
                mass(5:6) = hadron(rho)%mass
                mass(7:8) = mPi
                mStar(5:8) = mass(5:8)
                rat(i,j,k) = Ratio(SrtsStar,mStar(1:8),mass(1:8))
             else if ( k == 12 ) then !***** N_1 N_2 -> N_3 N_4 rho_5 rho_6 rho_7
                mass(1:4) = mN
                mStar(1:4) = mass(1:4) - shift
                mass(5:7) = hadron(rho)%mass
                mStar(5:7) = mass(5:7)
                rat(i,j,k) = Ratio(SrtsStar,mStar(1:7),mass(1:7))
             else if ( k == 13 ) then !***** N_1 N_2 -> N_3 N_4 rho_5 rho_6 rho_7 pi_8
                mass(1:4) = mN
                mStar(1:4) = mass(1:4) - shift
                mass(5:7) = hadron(rho)%mass
                mass(8) = mPi
                mStar(5:8) = mass(5:8)
                rat(i,j,k) = Ratio(SrtsStar,mStar(1:8),mass(1:8))
             else if ( k == 14 ) then !***** N_1 N_2 -> N_3 N_4 rho_5 rho_6 rho_7 rho_8
                mass(1:4) = mN
                mStar(1:4) = mass(1:4) - shift
                mass(5:8) = hadron(rho)%mass
                mStar(5:8) = mass(5:8)
                rat(i,j,k) = Ratio(SrtsStar,mStar(1:8),mass(1:8))
             else if ( k == 15 ) then !***** N_1 N_2 -> N_3 Y_4 K_5
                mass(1:3) = mN
                mass(4) = hadron(Lambda)%mass
                mStar(1:4) = mass(1:4) - shift
                mass(5) = mK
                mStar(5) = mass(5)
                rat(i,j,k) = Ratio(SrtsStar,mStar(1:5),mass(1:5))
             else if ( k == 16 ) then !***** N_1 N_2 -> N_3 Y_4 K_5 pi_6
                mass(1:3) = mN
                mass(4) = hadron(Lambda)%mass
                mStar(1:4) = mass(1:4) - shift
                mass(5) = mK
                mass(6) = mPi
                mStar(5:6) = mass(5:6)
                rat(i,j,k) = Ratio(SrtsStar,mStar(1:6),mass(1:6))
             else if ( k == 17 ) then !***** N_1 N_2 -> N_3 Y_4 K_5 pi_6 pi_7
                mass(1:3) = mN
                mass(4) = hadron(Lambda)%mass
                mStar(1:4) = mass(1:4) - shift
                mass(5) = mK
                mass(6:7) = mPi
                mStar(5:7) = mass(5:7)
                rat(i,j,k) = Ratio(SrtsStar,mStar(1:7),mass(1:7))
             else if ( k == 18 ) then !***** N_1 N_2 -> N_3 Y_4 K_5 pi_6 pi_7 pi_8
                mass(1:3) = mN
                mass(4) = hadron(Lambda)%mass
                mStar(1:4) = mass(1:4) - shift
                mass(5) = mK
                mass(6:8) = mPi
                mStar(5:8) = mass(5:8)
                rat(i,j,k) = Ratio(SrtsStar,mStar(1:8),mass(1:8))
             else if ( k == 19 ) then !***** N_1 N_2 -> N_3 Y_4 K_5 rho_6
                mass(1:3) = mN
                mass(4) = hadron(Lambda)%mass
                mStar(1:4) = mass(1:4) - shift
                mass(5) = mK
                mass(6) = hadron(rho)%mass
                mStar(5:6) = mass(5:6)
                rat(i,j,k) = Ratio(SrtsStar,mStar(1:6),mass(1:6))
             else if ( k == 20 ) then !***** N_1 N_2 -> N_3 Y_4 K_5 rho_6 pi_7
                mass(1:3) = mN
                mass(4) = hadron(Lambda)%mass
                mStar(1:4) = mass(1:4) - shift
                mass(5) = mK
                mass(6) = hadron(rho)%mass
                mass(7) = mPi
                mStar(5:7) = mass(5:7)
                rat(i,j,k) = Ratio(SrtsStar,mStar(1:7),mass(1:7))
             else if ( k == 21 ) then !***** N_1 N_2 -> N_3 Y_4 K_5 rho_6 pi_7 pi_8
                mass(1:3) = mN
                mass(4) = hadron(Lambda)%mass
                mStar(1:4) = mass(1:4) - shift
                mass(5) = mK
                mass(6) = hadron(rho)%mass
                mass(7:8) = mPi
                mStar(5:8) = mass(5:8)
                rat(i,j,k) = Ratio(SrtsStar,mStar(1:8),mass(1:8))
             else if ( k == 22 ) then !***** N_1 N_2 -> N_3 Y_4 K_5 rho_6 rho_7
                mass(1:3) = mN
                mass(4) = hadron(Lambda)%mass
                mStar(1:4) = mass(1:4) - shift
                mass(5) = mK
                mass(6:7) = hadron(rho)%mass
                mStar(5:7) = mass(5:7)
                rat(i,j,k) = Ratio(SrtsStar,mStar(1:7),mass(1:7))
             else if ( k == 23 ) then !***** N_1 N_2 -> N_3 Y_4 K_5 rho_6 rho_7 pi_8
                mass(1:3) = mN
                mass(4) = hadron(Lambda)%mass
                mStar(1:4) = mass(1:4) - shift
                mass(5) = mK
                mass(6:7) = hadron(rho)%mass
                mass(8) = mPi
                mStar(5:8) = mass(5:8)
                rat(i,j,k) = Ratio(SrtsStar,mStar(1:8),mass(1:8))
             else if ( k == 24 ) then !***** N_1 N_2 -> N_3 Y_4 K_5 rho_6 rho_7 rho_8
                mass(1:3) = mN
                mass(4) = hadron(Lambda)%mass
                mStar(1:4) = mass(1:4) - shift
                mass(5) = mK
                mass(6:8) = hadron(rho)%mass
                mStar(5:8) = mass(5:8)
                rat(i,j,k) = Ratio(SrtsStar,mStar(1:8),mass(1:8))
             else if ( k == 25 ) then !***** N_1 N_2 -> N_3 N_4 K_5 Kbar_6
                mass(1:4) = mN
                mStar(1:4) = mass(1:4) - shift
                mass(5:6) = mK
                mStar(5:6) = mass(5:6)
                rat(i,j,k) = Ratio(SrtsStar,mStar(1:6),mass(1:6))
             else if ( k == 26 ) then !***** N_1 N_2 -> N_3 N_4 K_5 Kbar_6 pi_7
                mass(1:4) = mN
                mStar(1:4) = mass(1:4) - shift
                mass(5:6) = mK
                mass(7) = mPi
                mStar(5:7) = mass(5:7)
                rat(i,j,k) = Ratio(SrtsStar,mStar(1:7),mass(1:7))
             else if ( k == 27 ) then !***** N_1 N_2 -> N_3 N_4 K_5 Kbar_6 pi_7 pi_8
                mass(1:4) = mN
                mStar(1:4) = mass(1:4) - shift
                mass(5:6) = mK
                mass(7:8) = mPi
                mStar(5:8) = mass(5:8)
                rat(i,j,k) = Ratio(SrtsStar,mStar(1:8),mass(1:8))
             else if ( k == 28 ) then !***** N_1 N_2 -> N_3 N_4 K_5 Kbar_6 rho_7
                mass(1:4) = mN
                mStar(1:4) = mass(1:4) - shift
                mass(5:6) = mK
                mass(7) = hadron(rho)%mass
                mStar(5:7) = mass(5:7)
                rat(i,j,k) = Ratio(SrtsStar,mStar(1:7),mass(1:7))
             else if ( k == 29 ) then !***** N_1 N_2 -> N_3 N_4 K_5 Kbar_6 rho_7 pi_8
                mass(1:4) = mN
                mStar(1:4) = mass(1:4) - shift
                mass(5:6) = mK
                mass(7) = hadron(rho)%mass
                mass(8) = mPi
                mStar(5:8) = mass(5:8)
                rat(i,j,k) = Ratio(SrtsStar,mStar(1:8),mass(1:8))
             else if ( k == 30 ) then !***** N_1 N_2 -> N_3 N_4 K_5 Kbar_6 rho_7 rho_8
                mass(1:4) = mN
                mStar(1:4) = mass(1:4) - shift
                mass(5:6) = mK
                mass(7:8) = hadron(rho)%mass
                mStar(5:8) = mass(5:8)
                rat(i,j,k) = Ratio(SrtsStar,mStar(1:8),mass(1:8))
             else if ( k == 31 ) then !***** N_1 N_2 -> N_3 Delta_4
                mass(1:3) = mN
                mass(4) = hadron(delta)%mass
                mStar(1:4) = mass(1:4) - shift
                rat(i,j,k) = Ratio(SrtsStar,mStar(1:4),mass(1:4))
             else if ( k == 32 ) then !***** N_1 N_2 -> N_3 Y_4 Kstar_5
                mass(1:3) = mN
                mass(4) = hadron(Lambda)%mass
                mStar(1:4) = mass(1:4) - shift
                mass(5) = hadron(kaonStar)%mass
                mStar(5) = mass(5)
                rat(i,j,k) = Ratio(SrtsStar,mStar(1:5),mass(1:5))
             else if ( k == 33 ) then !***** N_1 N_2 -> N_3 N_4 Kstar_5 Kbar_6
                mass(1:4) = mN
                mStar(1:4) = mass(1:4) - shift
                mass(5) = hadron(kaonStar)%mass
                mass(6) = mK
                mStar(5:6) = mass(5:6)
                rat(i,j,k) = Ratio(SrtsStar,mStar(1:6),mass(1:6))
             else if ( k == 34 ) then !***** N_1 N_2 -> N_3 Delta_4 pi_5
                mass(1:3) = mN
                mass(4) = hadron(delta)%mass
                mStar(1:4) = mass(1:4) - shift
                mass(5) = mPi
                mstar(5) = mass(5)
                rat(i,j,k) = Ratio(SrtsStar,mStar(1:5),mass(1:5))
             else if ( k == 35 ) then !***** N_1 N_2 -> N_3 N_4
                mass(1:4) = mN
                mStar(1:4) = mass(1:4) - shift
                rat(i,j,k) = Ratio(SrtsStar,mStar(1:4),mass(1:4))
             end if

          end do Loop_over_channels

       end do Loop_over_SrtsStar

    end do Loop_over_massShift

  end subroutine TabulateRatio



  !****************************************************************************
  !****f* XsectionRatios/Ratio
  ! NAME
  ! real function Ratio(srtsStar,mStar,mass)
  ! PURPOSE
  ! Compute the ratio of the in-medium and vacuum
  ! B_1 + B_2 -> B_3 + B_4 + M_5 + M_6 + M_7 + M_8  cross sections
  ! INPUT:
  ! * real, intent(in) :: srtsStar                  ! c.m. energy (GeV),
  ! * real, intent(in), dimension(:) :: mStar, mass ! Dirac and vacuum masses of all particles (GeV)
  ! OUTPUT:
  ! * real :: Ratio                                 ! ratio of the in-medium and vacuum cross sections
  ! NOTES:
  ! The particles 1 and 2 are incoming baryons, the particles 3 and 4 are outgoing baryons,
  ! the particles M_5, ..., M_8 are outgoing mesons.
  ! The routine will compute one of the ratios: B_1 + B_2 -> B_3 + B_4,
  ! B_1 + B_2 -> B_3 + B_4 + M_5, ..., B_1 + B_2 -> B_3 + B_4 + M_5 + ...+ M_8
  ! depending on the the size n of the arrays mStar and mass, which must be in the interval [4:8].
  !****************************************************************************
  real function Ratio(srtsStar,mStar,mass)

    use nBodyPhaseSpace
    use twoBodyTools, only: pCM

    real, intent(in) :: srtsStar                  ! c.m. energy (GeV),
    real, intent(in), dimension(:) :: mStar, mass ! Dirac and vacuum masses of all particles (GeV)

    integer :: i,n
    real :: srts, Ifac, IfacStar, phaseSpace, phaseSpaceStar

    n = size(mStar,dim=1)

    if ( n < 4 .or. n > 8 ) then
       write(*,*) ' In XsectionRatios/Ratio: wrong size of input array, n= ', n
       stop
    end if

    if ( n .ne. size(mass,dim=1) ) then
       write(*,*) ' In XsectionRatios/Ratio: different sizes of input arrays ', &
            & n, size(mass,dim=1)
       stop
    end if

    if ( srtsStar <= sum(mStar(1:2)) ) then
       write(*,*) ' In XsectionRatios/Ratio: srtsStar, mstar(1:2) : ', srtsStar, mstar(1:2)
       write(*,*) ' Sum of incoming Dirac masses', sum(mStar(1:2))
       stop
    end if

    if ( srtsStar <= sum(mStar(3:n)) + 0.001 ) then
       Ratio = 0.
       return
    end if

    ! c.m. energy in vacuum:
    srts = srtsStar - sum(mStar(3:n)) + sum(mass(3:n))

    ! Vacuum and in-medium flux factors:
    Ifac     = pcm(srts,mass(1),mass(2)) * srts
    IfacStar = pcm(srtsStar,mStar(1),mStar(2)) * srtsStar

    ! Vacuum and in-medium phase spaces:
    phaseSpace     = integrate_nBodyPS(srts,mass(3:n))
    phaseSpaceStar = integrate_nBodyPS(srtsStar,mStar(3:n))

    ! Cross section ratio = sigma_med / sigma_vac:
    if (IfacStar.gt.0. .and. phaseSpace.gt.0.) then
        Ratio = mStar(1)*mStar(2)/(mass(1)*mass(2)) &
              & *Ifac/IfacStar * phaseSpaceStar/phaseSpace
        do i = 3,n
           if (mass(i).gt.0.) Ratio = Ratio * mStar(i)/mass(i)
        end do
    end if

  end function Ratio


  !****************************************************************************
  !****f* XsectionRatios/Ratio_BaB
  ! NAME
  ! real function Ratio_BaB(pair,finalState)
  ! PURPOSE
  ! Compute the ratio of the in-medium and vacuum cross sections
  ! baryon + antibaryon -> mesons
  ! INPUT:
  ! * type(particle), dimension(1:2), intent(in) :: pair         ! incoming particles
  ! * type(particle), dimension(:), intent(in) :: finalState     ! outgoing particles
  ! OUTPUT:
  ! * real :: Ratio_BaB            ! ratio of the in-medium and vacuum cross sections
  ! NOTES:
  ! Maximum 6 mesons can be treated (otherwice returns Ratio_BaB=1.).
  !****************************************************************************
  real function Ratio_BaB(pair,finalState)

    use particleDefinition, only: particle, sqrtS
    use nBodyPhaseSpace
    use RMF, only: getRMF_flag
    use densitymodule, only: Particle4Momentum_RMF
    use twoBodyTools, only: sqrtS_Free,pCM

    type(particle), dimension(1:2), intent(in) :: pair         ! incoming particles
    type(particle), dimension(:), intent(in) :: finalState     ! outgoing particles

    real :: srtS_star, srtS_vacuum, srtS    ! c.m. energies (GeV),
    real :: mstar1, mstar2                   ! Dirac masses of incoming particles (GeV)
    real, dimension(0:3) :: momentum1, momentum2, momentum_tot

    integer :: n, i
    real :: Ifac_vacuum, Ifac, phaseSpace_vacuum, phaseSpace

    logical, parameter :: debugFlag = .false.

    ! Determine how many particles are in final state:
    n=0
    do i=1,size(finalState,dim=1)
       if (finalState(i)%Id.eq.0) then
          cycle
       else if (finalState(i)%Id.lt.0) then
          exit
       else
          n=n+1
       end if
    end do

    if ( n < 2 .or. n > 6 ) then
       Ratio_BaB=1.
       return
    end if

    srtS_star = sqrtS(pair,"Ratio_BaB, srtS_star")

    mstar1 = sqrtS(pair(1),'Ratio_BaB, mstar(1)')
    mstar2 = sqrtS(pair(2),'Ratio_BaB, mstar(2)')

    if ( .not.getRMF_flag() ) then
       srtS_vacuum=sqrtS_free(pair)
       srtS=srtS_star
    else
       srtS_vacuum  = srtS_star - mstar1 - mstar2 + pair(1)%mass + pair(2)%mass
       call Particle4Momentum_RMF(pair(1),momentum1)
       call Particle4Momentum_RMF(pair(2),momentum2)
       ! Compute srtS with canonical momenta:
       momentum_tot(0:3) = momentum1(0:3) + momentum2(0:3)
       srtS = momentum_tot(0)**2 - dot_product(momentum_tot(1:3),momentum_tot(1:3))
       srtS = sqrt(max(0.,srtS))
    end if

    if (debugFlag) then
       write(*,*) ' In Ratio_BaB: incoming particles: ', pair%Id,pair%antiparticle,pair%charge
       write(*,*) ' In Ratio_BaB: srtS_star, srtS, srtS_vacuum: ', srtS_star, srtS, srtS_vacuum
       write(*,*) ' In Ratio_BaB: outgoing particles: ', finalState(1:n)%Id
       write(*,*) ' In Ratio_BaB: sum of outgoing masses: ', sum(finalState(1:n)%mass)
    end if


    if ( srtS <= sum(finalState(1:n)%mass) + 0.001 ) then
       Ratio_BaB = 0.
       return
    end if

    ! Vacuum and in-medium flux factors:
    Ifac_vacuum = pcm(srtS_vacuum,pair(1)%mass,pair(2)%mass) * srtS_vacuum
    Ifac        = pcm(srtS_star,mstar1,mstar2) * srtS_star

    ! Vacuum and in-medium phase spaces:
    phaseSpace_vacuum = integrate_nBodyPS (srtS_vacuum, finalState(1:n)%mass)
    phaseSpace        = integrate_nBodyPS (srtS, finalState(1:n)%mass)

    ! Cross section ratio = sigma_med / sigma_vac:
    Ratio_BaB = mstar1*mstar2/(pair(1)%mass*pair(2)%mass) &
                * Ifac_vacuum/Ifac * phaseSpace/phaseSpace_vacuum

    if (debugFlag) then
       write(*,*) ' In Ratio_BaB: mstar1, mstar2: ', mstar1, mstar2
       write(*,*) ' In Ratio_BaB: Ifac_vacuum, Ifac: ',Ifac_vacuum, Ifac
       write(*,*) ' In Ratio_BaB: phaseSpace_vacuum, phaseSpace: ', phaseSpace_vacuum, phaseSpace
       write(*,*) ' In Ratio_BaB: Ratio_BaB', Ratio_BaB
    end if

  end function Ratio_BaB


end module XsectionRatios
