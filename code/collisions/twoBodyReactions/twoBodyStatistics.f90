!******************************************************************************
!****m* /twoBodyStatistics
! NAME
! module twoBodyStatistics
! PURPOSE
! This module contains routines for collecting 2-body statistics,
! e.g. sqrt(s) distributions and collision rates.
!******************************************************************************
module twoBodyStatistics

  implicit none
  private


  !****************************************************************************
  !****g* twoBodyStatistics/flag_sqrts
  ! PURPOSE
  ! If .true., then the calculation and output of the sqrts distributions
  ! from subroutine sqrts_distribution will be done
  ! SOURCE
  !
  logical, save :: flag_sqrts = .false.
  !****************************************************************************


  !****************************************************************************
  !****g* twoBodyStatistics/flag_rate
  ! PURPOSE
  ! If .true., then the calculation and output of the collision rates
  ! from subroutine rate will be done
  ! SOURCE
  !
  logical, save :: flag_rate = .false.
  !****************************************************************************

  !****************************************************************************
  !****g* twoBodyStatistics/flag_varirate
  ! PURPOSE
  ! If .true., then the calculation and output of the collision rates
  ! from subroutine varirate will be done
  ! SOURCE
  !
  logical, save :: flag_varirate = .false.
  !****************************************************************************

  !****************************************************************************
  !****g* twoBodyStatistics/varirate_chargeZero
  ! PURPOSE
  ! If .true., then all charge states are combined together
  ! SOURCE
  !
  logical, save :: varirate_chargeZero = .true.
  !****************************************************************************

  !****************************************************************************
  !****g* twoBodyStatistics/varirate_filterPhi
  ! PURPOSE
  ! If .true., then only channels involving a phi meson are reported
  ! SOURCE
  !
  logical, save :: varirate_filterPhi = .false.
  !****************************************************************************


  !****************************************************************************
  !****g* twoBodyStatistics/varirate_size
  ! PURPOSE
  ! size of array to hold all rates
  ! SOURCE
  !
  integer, save :: varirate_size = 100
  !****************************************************************************

  !****************************************************************************
  !****g* twoBodyStatistics/sqrts_mode
  ! PURPOSE
  ! This flag determines the way how sqrt(s) is calculated (if flag_sqrts = .true.).
  ! 1 = use vacuum sqrt(s)
  ! 2 = use in-medium, i.e. full sqrts
  ! SOURCE
  !
  integer, save ::  sqrts_mode = 1
  !****************************************************************************


  logical, save :: initFlag=.true.


  public :: sqrts_distribution, rate

contains


  !****************************************************************************
  !****s* twoBodyStatistics/init
  ! NAME
  ! subroutine init
  ! PURPOSE
  ! Reads in namelist "ColStat"
  ! INPUTS
  ! * (none)
  ! OUTPUT
  ! * Initializes global module variables
  !****************************************************************************
  subroutine init
    use output, only: Write_ReadingInput

    integer :: ios

    !**************************************************************************
    !****n* twoBodyStatistics/ColStat
    ! NAME
    ! NAMELIST ColStat
    ! PURPOSE
    ! Namelist which includes the input switches:
    ! * flag_sqrts
    ! * flag_rate
    ! * flag_varirate
    ! * sqrts_mode
    ! * varirate_chargeZero
    ! * varirate_size
    ! * varirate_filterPhi
    !**************************************************************************
    NAMELIST /ColStat/ flag_sqrts, flag_rate, flag_varirate, sqrts_mode, &
         varirate_chargeZero, varirate_size, varirate_filterPhi

    call Write_ReadingInput('ColStat',0)
    rewind(5)
    read(5,nml=ColStat,iostat=ios)
    call Write_ReadingInput('ColStat',0,ios)

    write(*,*) ' flag_sqrts = ', flag_sqrts
    write(*,*) ' flag_rate  = ', flag_rate
    write(*,*) ' flag_varirate  = ', flag_varirate
    if (flag_sqrts) &
      write(*,*) ' sqrts_mode  =  ', sqrts_mode

    if (flag_varirate) then
       write(*,*) ' chargeZero : ', varirate_chargeZero
       write(*,*) ' size       : ', varirate_size
       write(*,*) ' fiter phi  : ', varirate_filterPhi
    end if

    call Write_ReadingInput('ColStat',1)

    initFlag = .false.
  end subroutine init



  !****************************************************************************
  !****s* twoBodyStatistics/sqrts_distribution
  ! NAME
  ! subroutine sqrts_distribution (pair, itype, flag)
  ! PURPOSE
  ! Computes the distribution of colliding pairs of particles over sqrt(s).
  ! INPUTS
  ! * type(particle), dimension(1:2), intent(in) :: pair -- incoming pair of particles
  ! * integer, intent(in) :: itype
  ! * logical, optional, intent(in) :: flag  -- .true. do output
  !
  ! Possible values for itype:
  ! * 1 -- separate collision,
  ! * 2 -- collision in presence of a particle nearby before sqrt(s) correction,
  ! * 3 -- collision in presence of a particle nearby after sqrt(s) correction,
  ! * 4 -- elastic binary collisions
  ! * 5 -- inelastic binary collisions
  ! OUTPUT
  ! * If flag==.true., output is written to the files "dNdsqrts_BB.dat",
  ! "dNdsqrts_Bm.dat" and "dNdsqrts_mm.dat".
  !****************************************************************************
  subroutine sqrts_distribution (pair, itype, flag)

    use particleDefinition
    use IdTable
    use inputGeneral, only: numEnsembles, num_Runs_SameEnergy
    use twoBodyTools, only: sqrtS_free
    use densitymodule, only: Particle4Momentum_RMF
    use RMF, only: getRMF_flag
    use minkowski, only: abs4
    use histMC

    type(particle), dimension(1:2), intent(in) :: pair
    integer, intent(in) :: itype
    logical, optional, intent(in) :: flag

    real, parameter :: sqrts_min=0.2, dsqrts=0.02         ! minimum value of sqrt(s) and sqrt(s)-bin size (GeV)
    integer, parameter :: Nsqrts=400                      ! number of sqrt(s)-bins
    type(histogramMC), save :: dNdsqrts_BB, dNdsqrts_Bm, dNdsqrts_mm   ! sqrt(s)-distribution for BB/Bm
    real :: sqrtsStar, srtS, mstar1, mstar2
    real, dimension(0:3) :: momentum1, momentum2
    logical, save :: first = .true.
    real, save :: w

    if (initFlag) call init
    if (.not.flag_sqrts) return

    if (first) then
      call CreateHistMC (dNdsqrts_BB, "sqrt(s) distribution of baryon-baryon collisions", &
                         sqrts_min, sqrts_min+Nsqrts*dsqrts, dsqrts, 5)
      dNdsqrts_BB%xDesc = "sqrt(s) [GeV]"
      dNdsqrts_BB%yDesc(1:5) = (/ "dNdsqrts-sep.        ", "dNdsqrts-before corr.", &
                                  "dNdsqrts-after corr. ", "binary elastic       ", "binary inelastic     " /)
      call CreateHistMC (dNdsqrts_Bm, "sqrt(s) distribution of baryon-meson collisions", &
                         sqrts_min, sqrts_min+Nsqrts*dsqrts, dsqrts, 9)
      dNdsqrts_Bm%xDesc = "sqrt(s) [GeV]"
      dNdsqrts_Bm%yDesc(1:9) = (/ "pi-N      ", "eta-N     ", "rho-N     ", "omega-N   ", "Kbar-N    ", &
                                  "KstarBar-N", "Kbar-Delta", "pi-Y      ", "Y-N       " /)
      call CreateHistMC (dNdsqrts_mm, "sqrt(s) distribution of meson-meson collisions", &
                         sqrts_min, sqrts_min+Nsqrts*dsqrts, dsqrts, 5)
      dNdsqrts_mm%xDesc = "sqrt(s) [GeV]"
      dNdsqrts_mm%yDesc(1:5) = (/ "pi-pi ", "pi-K  ", "pi-rho", "K-K   ",  "other " /)
      first = .false.
      w = 1./float(numEnsembles*num_Runs_SameEnergy)  ! weight
    end if

    if (present(flag)) then

      if (flag) then
         !*********************************************************************
         !****o* twoBodyStatistics/dNdsqrts_BB.dat
         ! NAME
         ! file dNdsqrts_BB.dat
         ! PURPOSE
         ! Contains the sqrt(s) distribution of baryon-baryon collisions.
         ! Enabled by flag_sqrts.
         !*********************************************************************
         call WriteHistMC(dNdsqrts_BB,'dNdsqrts_BB.dat')
         !*********************************************************************
         !****o* twoBodyStatistics/dNdsqrts_Bm.dat
         ! NAME
         ! file dNdsqrts_Bm.dat
         ! PURPOSE
         ! Contains the sqrt(s) distribution of baryon-meson collisions.
         ! Enabled by flag_sqrts.
         !*********************************************************************
         call WriteHistMC(dNdsqrts_Bm,'dNdsqrts_Bm.dat')
         !*********************************************************************
         !****o* twoBodyStatistics/dNdsqrts_mm.dat
         ! NAME
         ! file dNdsqrts_mm.dat
         ! PURPOSE
         ! Contains the sqrt(s) distribution of meson-meson collisions.
         ! Enabled by flag_sqrts.
         !*********************************************************************
         call WriteHistMC(dNdsqrts_mm,'dNdsqrts_mm.dat')
         return
      end if

    end if

    select case (sqrts_mode)
    case (1)
      ! Vacuum sqrts:
      if (.not. getRMF_flag()) then
         srtS=sqrtS_free(pair)
      else
         sqrtsStar = abs4(pair(1)%momentum + pair(2)%momentum)
         mstar1 = abs4(pair(1)%momentum)
         mstar2 = abs4(pair(2)%momentum)
         srtS = sqrtsStar - (mstar1-pair(1)%mass) - (mstar2-pair(2)%mass)
      end if
    case (2)
      ! Full sqrts:
      if (.not. getRMF_flag()) then
         momentum1(:)=pair(1)%momentum
         momentum2(:)=pair(2)%momentum
      else
         call Particle4Momentum_RMF(pair(1),momentum1)
         call Particle4Momentum_RMF(pair(2),momentum2)
      end if
      srtS= abs4(momentum1 + momentum2)
    end select


    if (isBaryon(pair(1)%Id) .and. isBaryon(pair(2)%Id)) then      ! Baryon-baryon collisions
      call AddHistMC (dNdsqrts_BB, srtS, itype, w)
    end if

    if (itype>=4) then  ! binary

      if (isBaryon(pair(1)%Id).neqv.isBaryon(pair(2)%Id)) then  ! Baryon-meson collisions

         if (min(pair(1)%Id,pair(2)%Id)==nucleon) then
            select case (max(pair(1)%Id,pair(2)%Id))
            case (pion)
               call AddHistMC (dNdsqrts_Bm, srtS, 1, w)
            case (eta)
               call AddHistMC (dNdsqrts_Bm, srtS, 2, w)
            case (rho)
               call AddHistMC (dNdsqrts_Bm, srtS, 3, w)
            case (omegaMeson)
               call AddHistMC (dNdsqrts_Bm, srtS, 4, w)
            case (kaonBar)
               call AddHistMC (dNdsqrts_Bm, srtS, 5, w)
            case (kaonStarBar)
               call AddHistMC (dNdsqrts_Bm, srtS, 6, w)
            end select
         else if (min(pair(1)%Id,pair(2)%Id)==Delta .and. &
              max(pair(1)%Id,pair(2)%Id)==kaonBar) then
            call AddHistMC (dNdsqrts_Bm, srtS, 7, w)
         else if ((min(pair(1)%Id,pair(2)%Id)==Lambda .or. &
              min(pair(1)%Id,pair(2)%Id)==SigmaResonance) .and. &
              (max(pair(1)%Id,pair(2)%Id)==pion)) then
            call AddHistMC (dNdsqrts_Bm, srtS, 8, w)
         end if

      else if (isBaryon(pair(1)%Id).and.isBaryon(pair(2)%Id)) then  ! Baryon-baryon collisions

         if (min(pair(1)%Id,pair(2)%Id)==nucleon .and. &
              (max(pair(1)%Id,pair(2)%Id)==Lambda .or. &
              max(pair(1)%Id,pair(2)%Id)==SigmaResonance)) &
              call AddHistMC (dNdsqrts_Bm, srtS, 9, w)

      else if (isMeson(pair(1)%Id) .and. isMeson(pair(2)%Id)) then  ! Meson-meson collisions

        if (all(pair(1:2)%id==pion)) then
          call AddHistMC (dNdsqrts_mm, srtS, 1, w)
        else if (min(pair(1)%id,pair(2)%id)==pion .and. &
                 (max(pair(1)%id,pair(2)%id)==kaon .or. max(pair(1)%id,pair(2)%id)==kaonBar)) then
          call AddHistMC (dNdsqrts_mm, srtS, 2, w)
        else if (min(pair(1)%id,pair(2)%id)==pion .and. max(pair(1)%id,pair(2)%id)==rho) then
          call AddHistMC (dNdsqrts_mm, srtS, 3, w)
        else if ((pair(1)%id==kaon .or. pair(1)%id==kaonBar) .and. &
                 (pair(2)%id==kaon .or. pair(2)%id==kaonBar)) then
          call AddHistMC (dNdsqrts_mm, srtS, 4, w)
        else
          call AddHistMC (dNdsqrts_mm, srtS, 5, w)
        end if

      end if

    end if  ! itype if

  end subroutine sqrts_distribution


  !****************************************************************************
  !****s* twoBodyStatistics/rate
  ! NAME
  ! subroutine rate(teilchenIn,teilchenOut,time,flag)
  ! PURPOSE
  ! Computes the collision rates for the different types of collisions.
  ! INPUTS
  ! * type(particle), intent(in), dimension(:)  ::  teilchenIn -- incoming particles
  ! * type(particle), intent(in), dimension(:)  ::  teilchenOut -- outgoing particles
  ! * real, intent(in) :: time -- current time step
  ! * logical, optional, intent(in) :: flag  -- .true. do output of the accumulated collision numbers
  !
  ! OUTPUT
  ! * file "rate.dat"
  !****************************************************************************
  subroutine rate(partIn,partOut,time,flag)

  use particleDefinition
  use IdTable

  type(particle), intent(in), dimension(:)  :: partIn
  type(particle), intent(in), dimension(:)  :: partOut
  real, intent(in) :: time
  logical, optional, intent(in) :: flag

  integer, dimension(-121:121) :: num_in, num_out
  integer :: num_in_sum, num_out_sum,i,id,dLambda
  real, save, dimension(0:1,1:150) :: collision_number=0.
  real, save :: time_prev=-100.
  logical, save :: firstCall=.true.
  logical :: flagLambdaProd, flagLambdaAbs

  if (initFlag) call init
  if (flag_varirate) then
     call varirate(partIn, partOut, time)
     return
  end if
  if (.not.flag_rate) return

  if (firstCall) then
     open(file='rate.dat',UNIT=41,Status='Replace',Action='Write')
     write(41,*)'# time:  rates:'
     close(41)
     open(file='rate_accumulated_time.dat',UNIT=41,Status='Replace',Action='Write')
     write(41,*)'# time:  rates_accumulated:'
     close(41)
     firstCall=.false.
  end if

  if (present(flag)) then
     if (flag) then
        open(file='rate.dat',UNIT=41,Status='old',Position='Append',Action='Write')
        write(41,'(151(1x,e13.6))') time_prev,collision_number(1,1:150)
        close(41)
        open(file='rate_accumulated_time.dat',UNIT=41,Status='old',Position='Append',Action='Write')
        write(41,'(151(1x,e13.6))') time_prev,collision_number(0,1:150)
        close(41)
        open(file='rate_accumulated.dat',UNIT=41,Status='Replace',Action='Write')
        write(41,'(150(1x,e13.6))') collision_number(0,1:150)
        close(41)
        return
     end if
  end if

  if (time.ne.time_prev) then
     if (time_prev.gt.0.) then
        open(file='rate.dat',UNIT=41,Status='old',Position='Append',Action='Write')
        write(41,'(151(1x,e13.6))')  time_prev,collision_number(1,1:150)
        close(41)
        open(file='rate_accumulated_time.dat',UNIT=41,Status='old',Position='Append',Action='Write')
        write(41,'(151(1x,e13.6))') time_prev,collision_number(0,1:150)
        close(41)
     end if
     time_prev=time
     collision_number(1,:)=0.
  end if

  num_in=0
  do i=1,size(partIn,dim=1)
     id=partIn(i)%id
     if (partIn(i)%antiparticle) id=-id
     if (abs(id).le.121) num_in(id)=num_in(id)+1
  end do
  num_in_sum=sum(num_in(:))

  num_out=0
  do i = 1,size(partOut,dim=1)
     if (partOut(i)%ID <= 0) cycle
     id=partOut(i)%ID
     if (partOut(i)%antiparticle) id=-id
     if (abs(id).le.121) num_out(id)=num_out(id)+1
  end do
  num_out_sum=sum(num_out(:))

  flagLambdaProd=.false.
  flagLambdaAbs=.false.

  if ( sum(num_in(nucleon:F37_1950)).eq.1 .and. &
    & sum(num_in(-F37_1950:-nucleon)).eq.1 ) then         !********* Nbar(Rbar) N(R) collision ********

      if (num_out(kaon)+num_out(kaonStar).ge.1) &
         & collision_number(:,1)=collision_number(:,1)+1.        ! K(K^*) production

      if (num_out(kaonBar)+num_out(kaonStarBar).ge.1) &
         & collision_number(:,2)=collision_number(:,2)+1.        ! Kbar(K^*Bar) production

      if ( num_out(Lambda).ge.1 ) then
          collision_number(:,3)=collision_number(:,3)+num_out(Lambda)         ! Lambda production
          flagLambdaProd=.true.
      end if

      if ( num_out(SigmaResonance).ge.1 ) &
         & collision_number(:,4)=collision_number(:,4)+1.        ! Sigma production

      if ( num_out(Xi).ge.1 ) &
         & collision_number(:,5)=collision_number(:,5)+1.        ! Xi production

  else if ( sum(num_in(Lambda:sigma_1915)).eq.1 .and. &
         & sum(num_in(-F37_1950:-nucleon)).eq.1 ) then         !*********  Y(Y^*)  Nbar(Rbar) collision ********

      dLambda=num_out(Lambda)-num_in(Lambda)

      if ( dLambda.gt.0 ) then
         collision_number(:,6)=collision_number(:,6)+dLambda  !  Lambda production
         flagLambdaProd=.true.
      else
         collision_number(:,7)=collision_number(:,7)- dLambda !  Lambda absorption
         flagLambdaAbs=.true.
      end if

  else if ( sum(num_in(Lambda:sigma_1915)).eq.1 .and. &
         & sum(num_in(-sigma_1915:-Lambda)).eq.1 ) then         !*********  Y(Y^*)  Ybar(Ybar^*) collision ********

      dLambda=num_out(Lambda)-num_in(Lambda)

      if ( dLambda.gt.0 ) then
         collision_number(:,8)=collision_number(:,8)+dLambda  !  Lambda production
         flagLambdaProd=.true.
      else
         collision_number(:,9)=collision_number(:,9)- dLambda !  Lambda absorption
         flagLambdaAbs=.true.
      end if


  else if ( num_in(kaonBar)+num_in(kaonStarBar).eq.1 .and. &
         & sum(num_in(nucleon:F37_1950)).eq.1 ) then                         !********** Kbar(KstarBar) N(R) collision *******

      if ( num_out(Lambda).ge.1 ) then
          collision_number(:,11)=collision_number(:,11)+num_out(Lambda)       !  Lambda production
          flagLambdaProd=.true.
      end if

      if ( num_out(SigmaResonance).eq.1 ) &
         & collision_number(:,12)=collision_number(:,12)+1.       ! Sigma production

      if (sum(num_out(Sigma_1385:Sigma_1915))==1 .and. num_out(pion)==1) &
        collision_number(:,13)=collision_number(:,13)+1.       ! Y^* pi production

      if (sum(num_out(Sigma_1385:Sigma_1915))==1 .and. num_out_sum==1) &
        collision_number(:,14)=collision_number(:,14)+1.       ! Y^* production

      if ( num_out(Xi).eq.1 ) &
         & collision_number(:,15)=collision_number(:,15)+1.       !  Xi production

  else if ( num_in(kaon)+num_in(kaonStar).eq.1 .and. &
         & sum(num_in(nucleon:F37_1950)).eq.1 ) then                         !********** K(Kstar) N(R) collision *******

      if ( num_out(Lambda).ge.1 ) then
          collision_number(:,20)=collision_number(:,20)+num_out(Lambda)       ! Lambda production
          flagLambdaProd=.true.
      end if

  else if ( sum(num_in(kaon:kaonStarBar)).eq.1 .and. &
         & sum(num_in(Lambda:sigma_1915)).eq.1 ) then                         !********** K/Kbar(Kstar/KstarBar) Y(Y^*) collision *******

      if ( num_out(Xi).eq.1 ) &
         & collision_number(:,16)=collision_number(:,16)+1.       !  Xi production

      dLambda=num_out(Lambda)-num_in(Lambda)

      if ( dLambda.gt.0 ) then
          collision_number(:,18)=collision_number(:,18)+dLambda       !  Lambda production
          flagLambdaProd=.true.
      else
          collision_number(:,24)=collision_number(:,24)-dLambda       !  Lambda absorption
          flagLambdaAbs=.true.
      end if

  else if ( num_in(kaon)+num_in(kaonStar).eq.1 .and. &
         & sum(num_in(Xi:XiStar)).eq.1 ) then                                             !********** K(K^*) Xi(Xi^*) collision *******


      if ( num_out(kaonBar)+num_out(kaonStarBar).eq.1 ) &
         & collision_number(:,17)=collision_number(:,17)+1.       !  Kbar(K^*bar) production

      if ( num_out(Lambda).ge.1 ) then
          collision_number(:,19)=collision_number(:,19)+num_out(Lambda)       ! Lambda production
          flagLambdaProd=.true.
      end if

  else if (num_in(pion).eq.1 .and. sum(num_in(nucleon:F37_1950)).eq.1) then   !********** Pion N(R) collision *******

      if ( num_out(Lambda).eq.1 .and. num_out(kaon)+num_out(kaonStar).eq.1) then
          collision_number(:,21)=collision_number(:,21)+1.       !  Lambda K(K^*) production
          flagLambdaProd=.true.
      end if

      if ( num_out(SigmaResonance).eq.1 .and. num_out(kaon)+num_out(kaonStar).eq.1) &
         & collision_number(:,22)=collision_number(:,22)+1.     !  Sigma K(K^*) production

      if ( num_out(kaon)+num_out(kaonStar).eq.1 .and.  num_out(kaonBar)+num_out(kaonStarBar).eq.1 ) &
         & collision_number(:,23)=collision_number(:,23)+1.     !  K(K^*) Kbar(K^*bar) production

  else if (sum(num_in(pion:etaPrime)).eq.1 .and. sum(num_in(Xi:XiStar)).eq.1) then   !********** nonstr. meson-Xi(Xi^*) collision *******

      if ( sum(num_out(Lambda:sigma_1915)).eq.1 .and. (num_out(kaonBar).eq.1.or.num_out(kaonStarBar).eq.1) ) &
        & collision_number(:,27)=collision_number(:,27)+1.     !  Y(Y^*) Kbar(K^*bar) production

      if ( num_out(Lambda).ge.1 ) then
          collision_number(:,28)=collision_number(:,28)+num_out(Lambda)     ! Lambda production
          flagLambdaProd=.true.
      end if

      if ( num_out(XiStar).eq.1 ) &
        & collision_number(:,81)=collision_number(:,81)+1.     ! Xi^* production

  else if ( num_in(eta)+num_in(etaPrime).eq.1 .and. &
         & sum(num_in(nucleon:F37_1950)).eq.1) then                    !********** Eta(Eta') N(R) collision *******

       if (num_out(Lambda).eq.1 .and. num_out(kaon)+num_out(kaonStar).eq.1) then
          collision_number(:,31)=collision_number(:,31)+1.     !  Lambda K(K^*) production
          flagLambdaProd=.true.
       end if

       if (num_out(SigmaResonance).eq.1 .and. num_out(kaon)+num_out(kaonStar).eq.1) &
         & collision_number(:,32)=collision_number(:,32)+1.     !  Sigma K(K^*) production

       if ( num_out(kaon)+num_out(kaonStar).eq.1 .and.  num_out(kaonBar)+num_out(kaonStarBar).eq.1 ) &
         & collision_number(:,33)=collision_number(:,33)+1.     ! K(K^*) Kbar(K^*bar)  production

  else if (num_in(rho).eq.1 .and. sum(num_in(nucleon:F37_1950)).eq.1) then   !********** Rho N(R) collision *******

       if (num_out(Lambda).eq.1 .and. num_out(kaon)+num_out(kaonStar).eq.1) then
          collision_number(:,41)=collision_number(:,41)+1.     !  Lambda K(K^*) production
          flagLambdaProd=.true.
       end if

       if (num_out(SigmaResonance).eq.1 .and. num_out(kaon)+num_out(kaonStar).eq.1) &
         & collision_number(:,42)=collision_number(:,42)+1.     !  Sigma K(K^*) production

       if ( num_out(kaon)+num_out(kaonStar).eq.1 .and.  num_out(kaonBar)+num_out(kaonStarBar).eq.1 ) &
         & collision_number(:,43)=collision_number(:,43)+1.     !  K(K^*) Kbar(K^*bar) production

  else if (num_in(sigmaMeson).eq.1 .and. sum(num_in(nucleon:F37_1950)).eq.1) then   !********** Sigma N(R) collision *******

       if (num_out(Lambda).eq.1 .and. num_out(kaon)+num_out(kaonStar).eq.1) then
          collision_number(:,44)=collision_number(:,44)+1.     !  Lambda K(K^*) production
          flagLambdaProd=.true.
       end if

  else if (num_in(omegaMeson).eq.1 .and. sum(num_in(nucleon:F37_1950)).eq.1) then   !********** Omega N(R) collision *******

       if (num_out(Lambda).eq.1 .and. num_out(kaon)+num_out(kaonStar).eq.1) then
          collision_number(:,51)=collision_number(:,51)+1.     !  Lambda K(K^*) production
          flagLambdaProd=.true.
       end if

       if (num_out(SigmaResonance).eq.1 .and. num_out(kaon)+num_out(kaonStar).eq.1) &
         & collision_number(:,52)=collision_number(:,52)+1.     !  Sigma K(K^*) production

       if ( num_out(kaon)+num_out(kaonStar).eq.1 .and.  num_out(kaonBar)+num_out(kaonStarBar).eq.1 ) &
         & collision_number(:,53)=collision_number(:,53)+1.     !  K(K^*) Kbar(K^*bar) production

  else if (num_in(Lambda).eq.1 .and. sum(num_in(nucleon:F37_1950)).eq.1) then   !********** Lambda N(R) collision *******

      if ( num_out(SigmaResonance).eq.1 .and. num_out(Lambda).eq.0 ) then
          collision_number(:,61)=collision_number(:,61)+1.     !  Sigma production & Lambda absorption
          flagLambdaAbs=.true.
      end if

      if (sum(num_out(Sigma_1385:sigma_1915))==1 .and. num_out(Lambda)==0) then
          collision_number(:,78)=collision_number(:,78)+1.    ! Y^* production & Lambda absorption
          flagLambdaAbs=.true.
      end if

      if ( num_out(Xi).eq.1 .and. num_out(Lambda).eq.0 ) then
          collision_number(:,79)=collision_number(:,79)+1.    ! Xi production & Lambda absorption
          flagLambdaAbs=.true.
      end if

      if ( num_out(kaonBar)+num_out(kaonStarBar).eq.1 .and. num_out(Lambda).eq.0 ) then
          collision_number(:,82)=collision_number(:,82)+1.    ! Kbar(Kbar^*) production & Lambda absorption
          flagLambdaAbs=.true.
      end if


  else if (num_in(Lambda).eq.1 .and. sum(num_in(pion:etaPrime)).eq.1) then !***** Lambda - nonstr. meson collision *******

      if ( num_out(kaonBar)+num_out(kaonStarBar).eq.1 .and. num_out(Lambda).eq.0 ) then
          collision_number(:,62)=collision_number(:,62)+1.    ! Kbar(K^*bar) production & Lambda absorption
          flagLambdaAbs=.true.
      end if

      if (sum(num_out(Sigma_1385:sigma_1915))==1) then
          collision_number(:,63)=collision_number(:,63)+1.    ! Y^* production & Lambda absorption
          flagLambdaAbs=.true.
      end if

      if ( num_out(SigmaResonance).eq.1 ) then
          collision_number(:,54)=collision_number(:,54)+1.    ! Sigma production & Lambda absorption
          flagLambdaAbs=.true.
      end if

      if ( num_out(Xi)+num_out(XiStar).eq.1 ) then
          collision_number(:,57)=collision_number(:,57)+1.    ! Xi(Xi^*) production & Lambda absorption
          flagLambdaAbs=.true.
      end if


  else if (num_in(SigmaResonance).eq.1 .and. sum(num_in(nucleon:F37_1950)).eq.1) then   !********** Sigma N(R) collision *******

      if (num_out(Lambda).ge.1) then
          collision_number(:,65)=collision_number(:,65)+num_out(Lambda)     !  Lambda production
          flagLambdaProd=.true.
      end if

      if (sum(num_out(Sigma_1385:sigma_1915))==1) &
        collision_number(:,77)=collision_number(:,77)+1.    ! Y^* production

      if (num_out(Xi).eq.1) &
        & collision_number(:,80)=collision_number(:,80)+1.    ! Xi production


  else if (num_in(SigmaResonance).eq.1 .and. sum(num_in(pion:etaPrime)).eq.1) then !**** Sigma - nonstr. meson collision *******

      if (num_out(kaonBar)+num_out(kaonStarBar).eq.1) &
        &   collision_number(:,66)=collision_number(:,66)+1.    ! Kbar(K^*bar) production

      if (sum(num_out(Sigma_1385:sigma_1915))==1) &
        collision_number(:,67)=collision_number(:,67)+1.    ! Y^* production

      if (num_out(Lambda).ge.1) then
          collision_number(:,64)=collision_number(:,64)+num_out(Lambda)    ! Lambda production
          flagLambdaProd=.true.
      end if

  else if (sum(num_in(Sigma_1385:sigma_1915))==1 .and. sum(num_in(nucleon:F37_1950))==1) then !********** Y^* N(R) collision *******

      if (num_out(Lambda).ge.1) then
          collision_number(:,68)=collision_number(:,68)+num_out(Lambda)    ! Lambda production
          flagLambdaProd=.true.
      end if

      if (num_out(SigmaResonance).ge.1) &
        & collision_number(:,69)=collision_number(:,69)+num_out(SigmaResonance)    ! Sigma production

      if (num_out(Xi).eq.1) &
        & collision_number(:,70)=collision_number(:,70)+1.    ! Xi production

  else if (num_in(Xi).eq.1 .and. sum(num_in(nucleon:F37_1950)).eq.1) then   !********** Xi N(R) collision *******

      if (num_out(Lambda).eq.2) then
          collision_number(:,71)=collision_number(:,71)+2.     !  Lambda Lambda production
          flagLambdaProd=.true.
      else if (num_out(Lambda).eq.1 .and. num_out(SigmaResonance).eq.1) then
          collision_number(:,72)=collision_number(:,72)+1.     !  Lambda Sigma production
          flagLambdaProd=.true.
      else if (num_out(Lambda).eq.1 .and. num_out(kaonBar)+num_out(kaonStarBar).eq.1) then
          collision_number(:,55)=collision_number(:,55)+1.     !  Lambda Kbar(Kbar^*) production
          flagLambdaProd=.true.
      end if

  else if (num_in(XiStar).eq.1 .and. sum(num_in(nucleon:F37_1950)).eq.1) then   !********** Xi^* N(R) collision *******

      if ( num_out(Lambda).ge.1 ) then
          collision_number(:,56)=collision_number(:,56)+num_out(Lambda)       !  Lambda production
          flagLambdaProd=.true.
      end if

  else if (sum(num_in(Sigma_1385:sigma_1915))==1 .and. sum(num_in(pion:etaPrime))==1) then !********** Y^* -nonstr-meson collision *******

      if (num_out(kaonBar)+num_out(kaonStarBar).ge.1) &
        & collision_number(:,73)=collision_number(:,73)+1.    ! Kbar(K^*bar) production

      if (num_out(Lambda).eq.1) then
          collision_number(:,74)=collision_number(:,74)+1.    ! Lambda production
          flagLambdaProd=.true.
      end if

      if (num_out(SigmaResonance).eq.1) &
        & collision_number(:,75)=collision_number(:,75)+1.    ! Sigma production

      if (num_out(Xi).eq.1) &
        & collision_number(:,76)=collision_number(:,76)+1.    ! Xi production

  else if (sum(num_in(Sigma_1385:sigma_1915))==1 .and. num_in_sum==1) then !********** Y^* decay ************

      if ((num_out(kaonBar).eq.1.or.num_out(kaonStarBar).eq.1) .and. num_out(nucleon).eq.1) then
          collision_number(:,101)=collision_number(:,101)+1.     !  Kbar(KstarBar) N production
      else if (num_out(Lambda).eq.1 .and. (num_out(pion).eq.1.or.num_out(eta).eq.1)) then
          collision_number(:,102)=collision_number(:,102)+1.     !  pi(eta) Lambda production
          flagLambdaProd=.true.
      else if (num_out(SigmaResonance).eq.1 .and. num_out(pion).eq.1) then
          collision_number(:,103)=collision_number(:,103)+1.     !  pi Sigma production
      end if

  else if (sum(num_in(P11_1440:F37_1950)).eq.1 .and. num_in_sum.eq.1) then !********** R decay ************

      if (num_out(kaon).eq.1 .and. num_out(Lambda).eq.1) then
          collision_number(:,111)=collision_number(:,111)+1.
          flagLambdaProd=.true.
      end if

  else if (num_in(phi).eq.1 .and. num_in_sum.eq.1) then !********** phi decay

      if (num_out(kaonBar).eq.1 .and. num_out(kaon).eq.1) then
          collision_number(:,121)=collision_number(:,121)+1.     ! K Kbar production
      end if

  else if (num_in(XiStar).eq.1 .and. num_in_sum.eq.1) then !********** XiStar decay

      if (num_out(Xi).eq.1 .and. num_out(pion).eq.1) &
        &  collision_number(:,122)=collision_number(:,122)+1.    !  Xi pion production

  else if (sum(num_in(pion:etaPrime)).eq.2) then !********** nonstr. meson- nonstr. meson collision  **************

      if (sum(num_out(kaon:kaonStarBar)).eq.2) &
         & collision_number(:,123)=collision_number(:,123)+1.        ! K(K^*) Kbar(K^*bar) production

  else if (  num_in(kaon)+num_in(kaonStar).eq.1 .and. &
        &  num_in(kaonBar)+num_in(kaonStarBar).eq.1 ) then  !**** K(K^*) - Kbar(K^*bar) collision  **************

      if (sum(num_out(pion:etaPrime)).eq.2) &
        & collision_number(:,124)=collision_number(:,124)+1.  !  nonstr. meson- nonstr. meson production

  else if ( sum(num_in(nucleon:F37_1950)).eq.2 ) then   !********** N(R)-N(R) collision *******

      if ( num_out(Lambda).ge.1 ) then
         collision_number(:,141)=collision_number(:,141)+num_out(Lambda)  !  Lambda production
         flagLambdaProd=.true.
      end if

      if ( num_out(SigmaResonance).ge.1 ) &
       &  collision_number(:,142)=collision_number(:,142)+1.  !  Sigma production

      if ( num_out(kaon)+num_out(kaonStar).ge.1 ) &
       &  collision_number(:,143)=collision_number(:,143)+1.  !  K(K^*) production

      if ( num_out(kaonBar)+num_out(kaonStarBar).ge.1 ) &
       &  collision_number(:,144)=collision_number(:,144)+1.  !  Kbar(K^*bar) production

      if ( num_out(Xi).ge.1 ) &
       &  collision_number(:,145)=collision_number(:,145)+1.  !  Xi production

  else if ( sum(num_in(Lambda:sigma_1915)).eq.2 ) then  !**********  Y(Y^*)-Y(Y^*) collision *******

      dLambda=num_out(Lambda)-num_in(Lambda)

      if ( dLambda.gt.0 ) then
         collision_number(:,146)=collision_number(:,146)+dLambda  !  Lambda production
         flagLambdaProd=.true.
      else
         collision_number(:,147)=collision_number(:,147)- dLambda !  Lambda absorption
         flagLambdaAbs=.true.
      end if

  else if ( sum(num_in(Lambda:sigma_1915)).eq.1 .and. sum(num_in(Xi:XiStar)).eq.1 ) then  !**********  Xi(Xi^*)-Y(Y^*) collision *******

      dLambda=num_out(Lambda)-num_in(Lambda)

      if ( dLambda.gt.0 ) then
         collision_number(:,148)=collision_number(:,148)+dLambda  !  Lambda production
         flagLambdaProd=.true.
      else
         collision_number(:,149)=collision_number(:,149)- dLambda !  Lambda absorption
         flagLambdaAbs=.true.
      end if



  end if

  if ( num_out(Lambda).gt.num_in(Lambda) .and. .not.flagLambdaProd ) then
     write(*,*) ' Missed Lambda production: '
     write(*,*) ' Incoming Ids : '
     do i=1,size(partIn,dim=1)
        id=partIn(i)%id
        if (partIn(i)%antiparticle) id=-id
        write(*,*) id
     end do
     write(*,*) ' Outgoing Ids : '
     do i = 1,size(partOut,dim=1)
        if (partOut(i)%ID <= 0) cycle
        id=partOut(i)%ID
        if (partOut(i)%antiparticle) id=-id
        write(*,*) id
     end do
  else if ( num_out(Lambda).lt.num_in(Lambda) .and. .not.flagLambdaAbs ) then
     write(*,*) ' Missed Lambda absorption: '
     do i=1,size(partIn,dim=1)
        id=partIn(i)%id
        if (partIn(i)%antiparticle) id=-id
        write(*,*) id
     end do
     write(*,*) ' Outgoing Ids : '
     do i = 1,size(partOut,dim=1)
        if (partOut(i)%ID <= 0) cycle
        id=partOut(i)%ID
        if (partOut(i)%antiparticle) id=-id
        write(*,*) id
     end do
  end if

  end subroutine rate

  !****************************************************************************
  !****is* twoBodyStatistics/filterPhi
  ! NAME
  ! logical function filterPhi()
  ! PURPOSE
  ! apply filter on precence of phi meson
  !****************************************************************************
  logical function filterPhi(partIn,partOut)
    use particleDefinition

    type(particle), intent(in), dimension(:)  :: partIn
    type(particle), intent(in), dimension(:)  :: partOut

    integer :: i

    filterPhi = .true.

    do i=1,size(partIn)
       if (partIn(i)%ID == 110) return  ! eigentlich 107
    end do
    do i=1,size(partOut)
       if (partOut(i)%ID == 110) return  ! eigentlich 107
    end do

    filterPhi = .false.
    return

  end function filterPhi

  !****************************************************************************
  !****s* twoBodyStatistics/varirate
  ! NAME
  ! subroutine varirate(partIn,partOut,time)
  ! PURPOSE
  ! Computes the collision rates for alltypes of collisions.
  !
  ! Here the list of input and output channels are generated on the fly,
  ! thus this routine calculates the "variate-rate", or short: "varirate".
  !
  ! INPUTS
  ! * type(particle), intent(in), dimension(:) :: partIn -- incoming particles
  ! * type(particle), intent(in), dimension(:) :: partOut -- outgoing particles
  ! * real, intent(in) :: time -- current time step
  !
  ! OUTPUT
  ! * file "VariRate.particles.dat", "VariRate.rates.dat"
  !****************************************************************************
  subroutine varirate(partIn,partOut,time)

    use particleDefinition
    use IdTable
    use preEventDefinition
    use PreEvListDefinition
    use CallStack, only: Traceback
    use PreEvList, only: CreateSortedPreEvent, ComparePreEvent, &
         PreEvList_INIT, PreEvList_INSERT, PreEvList_Print, &
         PreEvList_GET, PreEvList_PrintEntry
    use twoBodyTools, only: LogicMatrix

    type(particle), intent(in), dimension(:)  :: partIn
    type(particle), intent(in), dimension(:)  :: partOut
    real, intent(in) :: time

    real, save :: time_prev=-100.
    logical, save :: firstCall=.true.

    type(tPreEvList), save :: ListIn, ListOut
    type(tPreEvListEntry) :: ePreEvIn, ePreEvOut, Vin, Vout
    integer :: iPosIn, iPosOut
    integer :: i,j,h,lev, scenario

    integer, dimension(:,:), allocatable, save :: EntryArr
    integer, save :: nEntryArr

    character*(*),dimension(0:3),parameter :: sS = (/ "  ","BB","BM","MM" /)

    logical :: accept

    if (firstCall) then
       call PreEvList_INIT(ListIn)
       call PreEvList_INIT(ListOut)

       nEntryArr = varirate_size
       allocate(EntryArr(nEntryArr,nEntryArr))
       EntryArr = 0

       firstCall = .false.
    end if

    if (time.ne.time_prev) then ! === do output ===

       open(313,file='VariRate.particles.dat')
       write(313,'("# time =",f13.5)') time_prev
       call PreEvList_Print(313, ListIn, 1.0, doSort=.false.)
       call PreEvList_Print(313, ListOut, 1.0, doSort=.false.)
       close(313)

       open(313,file='VariRate.rates.dat')
       write(313,'("# time =",f13.5)') time_prev

       do i=1,ListIn%nEntries
          do j=i,ListOut%nEntries

             if (EntryArr(i,j) > 0) then

                if (.not.PreEvList_GET(ListIn, Vin, i)) call Traceback()
                if (.not.PreEvList_GET(ListOut, Vout, j)) call Traceback()

                h = EntryArr(i,j)+EntryArr(j,i)
                lev = min(9, int(abs(EntryArr(i,j)-EntryArr(j,i))*10.0/h))
                scenario = 0

                if (ComparePreEvent(Vin%PreE,Vout%PreE) < 0) then
                   if (size(Vin%PreE)==2) &
                        scenario = LogicMatrix(isMeson(Vin%PreE(1)%ID),isMeson(Vin%PreE(2)%ID))

                   write(313,'("(",A,") (",i1,") [",i1,"<-->",i1,"] ")',advance='no') &
                        sS(scenario), lev, size(Vin%PreE),size(Vout%PreE)
                   call PreEvList_PrintEntry(313,Vin,1.0,doWeight=.false.,noAdvance=.true.)
                   write(313,'(" --> ")',advance='no')
                   call PreEvList_PrintEntry(313,Vout,1.0,doWeight=.false.,noAdvance=.true.)
                   write(313,'(" : ",2i10," # ",i10,f9.5)') EntryArr(i,j),EntryArr(j,i), &
                        h, (EntryArr(i,j)-EntryArr(j,i)) * 1.0/h

                else if (ComparePreEvent(Vout%PreE,Vin%PreE) < 0) then
                   if (size(Vout%PreE)==2) &
                        scenario = LogicMatrix(isMeson(Vout%PreE(1)%ID),isMeson(Vout%PreE(2)%ID))
                   write(313,'("(",A,") (",i1,") [",i1,"<-->",i1,"] ")',advance='no') &
                        sS(scenario), lev, size(Vout%PreE),size(Vin%PreE)
                   call PreEvList_PrintEntry(313,Vout,1.0,doWeight=.false.,noAdvance=.true.)
                   write(313,'(" --> ")',advance='no')
                   call PreEvList_PrintEntry(313,Vin,1.0,doWeight=.false.,noAdvance=.true.)
                   write(313,'(" : ",2i10," # ",i10,f9.5)') EntryArr(j,i),EntryArr(i,j), &
                        h, (EntryArr(j,i)-EntryArr(i,j)) * 1.0/h
                else
                   if (size(Vin%PreE)==2) &
                        scenario = LogicMatrix(isMeson(Vin%PreE(1)%ID),isMeson(Vin%PreE(2)%ID))
                   write(313,'("(",A,") (",i1,") [",i1,"====",i1,"] ")',advance='no') &
                        sS(scenario), lev, size(Vin%PreE),size(Vout%PreE)
                   call PreEvList_PrintEntry(313,Vin,1.0,doWeight=.false.,noAdvance=.true.)
                   write(313,'(" --> ")',advance='no')
                   call PreEvList_PrintEntry(313,Vout,1.0,doWeight=.false.,noAdvance=.true.)
                   write(313,'(" : ",2i10," # ",i10,f9.5)') EntryArr(i,j),EntryArr(j,i), &
                        h, (EntryArr(i,j)-EntryArr(j,i)) * 1.0/h

                end if
             end if

          end do
       end do
       close(313)

       time_prev = time
    end if

    ! ===== Now we are doing the actual work: the insertion ======

    if (varirate_filterPhi) then
       if (.not.filterPhi(partIn,partOut)) return
    end if

    if (.not.CreateSortedPreEvent(partIn,ePreEvIn%preE,chargeZero=varirate_chargeZero)) then
       call Traceback('Problems in CollHist_UpdateHist: preEvIn. STOP!')
    end if
    if (.not.CreateSortedPreEvent(partOut,ePreEvOut%preE,chargeZero=varirate_chargeZero)) then
       call Traceback('Problems in CollHist_UpdateHist: preEvOut. STOP!')
    end if

!    if ((size(ePreEvIn%preE) <= 2).and.(size(ePreEvOut%preE) <= 2)) return

    ! In the following, we insert both the in and the out channel in both lists.
    ! This is done to keep both lists simultanous. This is some overkill, but
    ! offers two different counters (->weight) for the in- and out-channels.

    ePreEvIn%weight = 0.0
    call PreEvList_INSERT(ListOut, ePreEvIn)
    ePreEvIn%weight = 1.0
    call PreEvList_INSERT(ListIn, ePreEvIn, iPosIn)

    ePreEvOut%weight = 0.0
    call PreEvList_INSERT(ListIn, ePreEvOut)
    ePreEvOut%weight = 1.0
    call PreEvList_INSERT(ListOut, ePreEvOut, iPosOut)

    if ((iPosIn > nEntryArr).or.(iPosOut > nEntryArr)) then
       write(*,*) 'more than ',nEntryArr,' input/output channels'
       call Traceback()
    end if

    ! === Now the cross array is filled: ===

    EntryArr(iPosIn,iPosOut) = EntryArr(iPosIn,iPosOut) + 1

  end subroutine varirate

end module twoBodyStatistics
