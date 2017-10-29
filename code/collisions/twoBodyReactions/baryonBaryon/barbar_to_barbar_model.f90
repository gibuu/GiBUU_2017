!******************************************************************************
!****m* /barBar_to_barBar_model
! NAME
! module barBar_to_barBar_model
! PURPOSE
! A model for nucleon delta scattering, based on one pion exchange.
!******************************************************************************
module barBar_to_barBar_model
  implicit none
  private
  ! All units in GeV

  ! Couplings according to Pascalutsa and Vanderhaeghen hep-ph/0511261v1
  real, parameter :: f_pi= 0.0924
  real, parameter :: g_A = 1.267
  !real, parameter :: h_A_deltaDelta = 9./5.*g_A
  real, parameter :: h_A_deltaNuc = 3./1.414213562*g_A  ! 3./sqrt(2.)*g_A     this change was required by Sun Studio 12.0


  ! Couplings according to Dmitriev/Sushkov, NPA 459 (1986)
  real, parameter :: f_pi_dimi= 1.008
  real, parameter :: f_pi_star_dimi = 2.202
  real, parameter :: f_pi_deltadelta_doenges=4./5.*f_pi_dimi


  integer, parameter :: dmitriev = 1
  integer, parameter :: pascalutsa = 2

  !****************************************************************************
  !****g* barBar_to_barBar_model/couplings_switch
  ! SOURCE
  !
  integer, save :: couplings_switch = 2
  ! PURPOSE
  ! Possible values:
  ! * 1 = use couplings according to Dmitriev
  ! * 2 = use couplings according to Pascalutsa (default)
  !****************************************************************************


  !****************************************************************************
  !****g* barBar_to_barBar_model/lambda_cutoff
  ! SOURCE
  !
  real   , save :: lambda_cutoff = 0.6
  ! PURPOSE
  ! Cutoff parameter in the form factor for ND->ND
  ! Possible values:
  ! * 0.6 (Dmitriev, default)
  ! * 1.2 (Doenges)
  !****************************************************************************


  logical, save :: firstTime=.true.

  public :: NN_ND_model, ND_ND_model, ND_ND_chooseCharge
  public :: MSquared_NN_ND, MSquared_ND_ND

contains


  !****************************************************************************
  !****s* barBar_to_barBar_model/readInput
  ! NAME
  ! subroutine readInput
  ! PURPOSE
  ! Reads input in jobcard out of namelist "BarBar_to_barBar_model".
  !****************************************************************************
  subroutine readInput
    use output
    integer :: ios

    !**************************************************************************
    !****n* barBar_to_barBar_model/BarBar_to_barBar_model
    ! NAME
    ! NAMELIST /BarBar_to_barBar_model/
    ! PURPOSE
    ! This namelist includes the following switches:
    ! * couplings_switch
    ! * lambda_cutoff
    !**************************************************************************
    NAMELIST /barBar_to_barBar_model/ couplings_switch,lambda_cutoff

    call Write_ReadingInput('barBar_to_barBar_model',0)
    rewind(5)
    read(5,nml=barBar_to_barBar_model,IOSTAT=IOS)
    call Write_ReadingInput('barBar_to_barBar_model',0,ios)

    write(*,*) '  Coupling Switch',couplings_switch
    write(*,*) '  Lambda         ',lambda_cutoff
    call Write_ReadingInput('barBar_to_barBar_model',1)

    firstTime = .false.
  end subroutine readInput


  !****************************************************************************
  !****f* barBar_to_barBar_model/NN_ND_model
  ! NAME
  ! real function NN_ND_model (srts)
  ! PURPOSE
  ! Evaluates cross section for NN -> NDelta (no isospin factors included: gives pp -> n \Delta^++ cross section) in mB
  ! INPUTS
  ! * real, intent (in) :: srts
  !****************************************************************************
  real function NN_ND_model (srts)
    use constants, only: pi, GeVSquared_times_mb, mN
    real, intent (in) :: srts

    if (firstTime) call readinput

    NN_ND_model = 0.
    if (srts<2.*mN) return

    NN_ND_model = 1./(64.*pi**2*srts**2) * MSquared_NN_ND (srts,999.,999.,.true.,.true.) / GeVSquared_times_mb

  end function NN_ND_model


  !****************************************************************************
  !****f* barBar_to_barBar_model/ND_ND_model
  ! NAME
  ! real function ND_ND_model(srts,idIn,chargeIn,idOut,chargeOut,massIn,chargeSum)
  ! PURPOSE
  ! Evaluates cross section for N Delta -> N Delta in mB.
  ! INPUTS
  ! * integer, dimension (1:2), intent(in) :: idIn,chargeIn,idOut,chargeOut  -- ids and charges of in- and outgoing particles
  ! * real, dimension (1:2), intent(in)    :: massIn                         -- masses of incoming particles
  ! * logical, intent(in)                  :: chargeSum                      -- .true.= sum over the final state charges
  !****************************************************************************
  real function ND_ND_model (srts,idIn,chargeIn,idOut,chargeOut,massIn,chargeSum)
    use idTable, only: Delta
    use constants, only: pi, GeVSquared_times_mb, mN

    integer, dimension (1:2) ,intent (in) :: idIn,chargeIn,idOut,chargeOut
    real, dimension (1:2) ,intent (in) :: massIn
    real, intent (in) :: srts
    real :: deltaMass_initial
    integer, dimension (1:2) :: chargeSorted_In,chargeSorted_OUT
    logical, intent(in) :: chargeSum

    if (firstTime) call readinput

    ND_ND_model = 0.
    if (srts<2.*mN) return

    if (idIn(1)==delta) then
      chargeSorted_In = (/chargeIn(2),chargeIn(1)/)
      deltaMass_initial = massin(1)
    else
       chargeSorted_In = (/chargeIn(1),chargeIn(2)/)
       deltaMass_initial = massin(2)
    end if

    if (idOUT(1)==delta) then
      chargeSorted_Out = (/chargeOut(2),chargeOut(1)/)
    else
      chargeSorted_Out = (/chargeOut(1),chargeOut(2)/)
    end if

    ND_ND_model = 1./(64.*pi**2*srts**2)/GeVSquared_times_mb &
                  * MSquared_ND_ND (chargeSorted_IN,chargeSorted_OUT,srts,999.,deltaMass_initial,999.,.true.,.true.,chargeSum)

  end function ND_ND_model


  !****************************************************************************
  !****f* barBar_to_barBar_model/ND_ND_chooseCharge
  ! NAME
  ! function ND_ND_chooseCharge (idIn,idOut,chargeIn) result(chargeOut)
  ! PURPOSE
  ! Does a Monte Carlo decision for the final state charges of N Delta -> N Delta based on isospin factors.
  ! INPUTS
  ! * integer, dimension (1:2) ,intent (in)       :: idIn,chargeIn,idOut -- ids of in- and outgoing particles, charges of incoming ones
  ! OUTPUT
  ! * integer, dimension(1:2) :: chargeOut
  ! NOTES
  !****************************************************************************
  function ND_ND_chooseCharge (idIn,idOut,chargeIn) result(chargeOut)
    use idtable, only: delta
    use random, only: rn

    integer, dimension(1:2) :: chargeOut
    integer, dimension(1:2), intent(in ) :: idIn,idOut,chargeIn
    integer, dimension(1:2) :: chargeSorted_IN
    real , dimension(0:1,-1:2) :: wahrscheinlichkeit
    integer :: nuccharge,deltaCharge
    real    :: totalWahrscheinlichkeit,x

    if (firstTime) call readinput

    if (idIn(1).eq.delta) then
       chargeSorted_In=(/chargeIn(2),chargeIn(1)/)
    else
       chargeSorted_In=(/chargeIn(1),chargeIn(2)/)
    end if

    wahrscheinlichkeit=0.
    totalWahrscheinlichkeit=0.
    do deltaCharge=-1,2
      do nucCharge=0,1
        wahrscheinlichkeit(nucCharge,deltaCharge)=isospinFactor_ND_ND(chargeSorted_In,(/nucCharge,deltaCharge/))
        totalWahrscheinlichkeit=totalWahrscheinlichkeit+wahrscheinlichkeit(nucCharge,deltaCharge)
      end do
    end do

    x=rn()*totalWahrscheinlichkeit
    deltaLOOP: do deltaCharge=-1,2
      do nucCharge=0,1
        if (wahrscheinlichkeit(nucCharge,deltaCharge).ge.x) then
          if (idOUT(1).eq.delta) then
            chargeOUT=(/deltaCharge,nucCharge/)
            exit deltaLoop
          else
            chargeOUT=(/nucCharge,deltaCharge/)
            exit deltaLoop
          end if
        end if
        x=x-wahrscheinlichkeit(nucCharge,deltaCharge)
      end do
    end do deltaLOOP

    ! Final Check
    if (sum(chargeIn)/=sum(chargeOut)) then
       write(*,*) 'Error in ND_ND_chooseCharge: Charge not conserved!!!'
       write(*,*) 'IN:', chargeIN
       write(*,*) 'OUT:', chargeOUT
       STOP 'barbar_to_barbar_model.f90'
    end if
  end function ND_ND_chooseCharge


  !****************************************************************************
  !****f* barBar_to_barBar_model/MSquared_NN_ND
  ! NAME
  ! real function MSquared_NN_ND (srts,theta_in,mass_in,integrated_theta,integrated_mass)
  ! NOTES
  ! This function  evaluates :
  ! Matrix Element^2 * p_cm(final) / p_cm(initial) * SpecFunc_Delta
  ! for NN-> NDelta according to one-pion exchange. See Dmitriev NPA 459.
  ! All variables in CM System!!
  ! INPUTS
  ! * real :: srts
  ! * real :: theta_in, mass_in -- angle in CM system in radians, mass of Delta
  ! * logical :: integrated_theta, integrated_mass -- .true.= give integrated cross section
  ! NOTES
  ! Up to now, it only works properly for pp -> n Delta^++
  !****************************************************************************
  real function MSquared_NN_ND (srts,theta_in,mass_in,integrated_theta,integrated_mass)
    use constants, only: pi, mPi, mN
    !use quadpack, only : qag

    real   , intent(in) :: srts
    real   , intent(in) :: theta_in, mass_in      ! angle in CM system in radians
    logical, intent(in) :: integrated_theta, integrated_mass ! .true.= give integrated cross section

    integer,parameter :: numsteps_theta=1000
    integer,parameter :: numsteps_mass=100

    real :: cosTheta, dcosTheta,theta,mass, dmass
    integer :: i,j
    !integer :: ier, neval
    !real :: abserr, integral

    if (firstTime) call readinput

    !Check normalization:
    !call qag (norm_integrand, baryon(nucleon)%mass+meson(pion)%mass,6.5,0.01, 0.01, 3, integral, abserr, neval, ier )
    !write(*,*) 'Normalization of delta:',ier, integral
    !stop

    MSquared_NN_ND=0.
    dcosTheta=2./float(numsteps_theta)
    dmass=(srts-2*mN-mPi)/float(numsteps_mass)

    ! Integrate over mass
    massLoop: do j=1,numsteps_mass
       if (.not.integrated_mass) then
          mass=mass_in
       else
          mass=mN+mPi+float(j)*dmass
       end if
       ! Integrate over cosTheta
       thetaLoop: do i=1,numSteps_theta
          if (.not.integrated_theta) then
             theta=theta_in
          else
             cosTheta=-1.+float(i)*dcosTheta
             theta=acos(cosTheta)
          end if
          MSquared_NN_ND = MSquared_NN_ND + ms_NN(theta,srts,mass)
          ! If no theta integration shall be performed, then we exit the integration loop after the first cycle
          if (.not.integrated_theta) exit thetaLoop
       end do thetaLoop
       ! If no mass integration shall be performed, then we exit the integration loop after the first cycle
       if (.not.integrated_mass) exit massLoop
    end do massLoop

    if (integrated_mass)  MSquared_NN_ND=MSquared_NN_ND*dmass
    if (integrated_theta) MSquared_NN_ND=MSquared_NN_ND*2.*pi*dcosTheta

  contains

    real function ms_NN(theta,srts,mass)
!       use particleProperties, only: hadron
      use matrix_module, only: matrixMult,unit4,trace
      use minkowski, only: SP, slashed,metricTensor
      use twoBodyTools, only: pcm
      use IDTABLE, only: delta
      use constants, only: mN, mPi
      use spinProjector, only: spin32proj

      real, intent(in) :: theta, srts,mass

      real :: couplings!,t,m1,mD
      real, dimension (0:3) :: a,b, p,q,pp,qp
      complex, dimension(0:3,0:3) :: mN4,lambda
      real :: pcm_initial,pcm_final,direct_traces,inter_trace!,direct_perhand,direct_dimi
      complex :: direct, crossed, interference
      integer:: mu,nu
      logical, parameter :: debug_ms=.false.

!       mD = hadron(delta)%mass

      if (couplings_switch.eq.Pascalutsa) then
         couplings=(h_A_deltaNuc/(2.*f_pi))**2*(g_A/(2.*f_pi))**2
      else if (couplings_switch.eq.Dmitriev) then
         couplings=(f_pi_star_dimi/mPi)**2    *   (f_pi_dimi/mPi)**2
      else
         stop 'strange couplings in barbar_to_barbar_model'
      end if

      mN4=unit4*mN
      !call printMatrix(mn4)
      !stop

      if (srts.le.mN+mass.or.srts.le.2.*mn) then
         ms_NN=0.
         return
      end if

      !************************************************************************
      ! Setting up kinematics in CM-FRAME
      ! Initial particles
      pcm_initial= pCM(srts, mN , mN )
      pcm_final  = pCM(srts, mN , mass )


      p(1:3)=(/0.,0.,pcm_initial/)
      q(1:3)=-p(1:3)

      ! Final ones
      pp(1:3)=(/sin(theta)*pcm_final,0.,cos(theta)*pcm_final/)
      qp(1:3)=-pp(1:3)

      ! Initial Nucleon energies:
      p(0)  =sqrt( mN**2  + Dot_product (p(1:3) ,p(1:3)  ) )
      q(0)  =sqrt( mN**2  + Dot_product (q(1:3) ,q(1:3)  ) )

      ! Final nucleon energy:
      pp(0) =sqrt( mN**2  + Dot_product (pp(1:3),pp(1:3) ) )
      ! Final delta energy :
      qp(0) =sqrt( mass**2+ Dot_product (qp(1:3),qp(1:3) ) )

      ! CHECKS:
      if (abs(p(0)+q(0)-pp(0)-qp(0)).gt.1E-5) then
         write(*,*) 'Energies not conserved',p(0)+q(0),pp(0)+qp(0)
         stop
      end if

      ! Exchanged pions
      a=p-pp
      b=p-qp
      !************************************************************************


      if (debug_ms) then
         ! Evaluate the trace numerically. Since "direct_byHand" is faster, this is only needed for debugging
         lambda=0
         do mu=0,3
            do nu=0,3
               lambda=lambda+a(nu)*a(mu)* spin32proj(delta,mu,nu,qp)*metricTensor(mu,mu)*metricTensor(nu,nu)
            end do
         end do
         direct_traces=trace( &
              &        MatrixMult( (slashed(pp)+mN4),slashed(a),(slashed(p)-mN4),slashed(a)) &
              &       ) &
              & *trace( &
              &        MatrixMult( lambda ,(slashed(q)+mN4)) &
              &       )


         direct=1./((mpi**2-SP(a,a))**2)*(  &
              &  trace( &
              &        MatrixMult( (slashed(pp)+mN4),slashed(a),(slashed(p)-mN4),slashed(a)) &
              &       ) &
              & *trace( &
              &        MatrixMult( lambda ,(slashed(q)+mN4)) &
              &       ) &
              & ) *cutoff(SP(a,a))**4
      end if

      direct=direct_byHand(SP(a,a),mass)/((mpi**2-SP(a,a))**2)*cutoff(SP(a,a))**4


      ! Crossed: replace a-> b and q <-> p
      if (debug_ms) then
         ! Evaluate the trace numerically. Since "direct_byHand" is faster, this is only needed for debugging

         lambda=0
         do mu=0,3
            do nu=0,3
               lambda(0:3,0:3)=lambda+b(nu)*b(mu)* spin32proj(delta,mu,nu,qp)*metricTensor(mu,mu)*metricTensor(nu,nu)
            end do
         end do

         crossed=1./((mpi**2-SP(b,b))**2)*(  &
              &  trace( &
              &        MatrixMult( (slashed(pp)+mN4),slashed(b),(slashed(q)-mN4),slashed(b)) &
              &       ) &
              & *trace( &
              &        MatrixMult( lambda ,(slashed(p)+mN4)) &
              &       ) &
              & ) *(cutoff(SP(b,b))**4)
      end if

      crossed=direct_byHand(SP(b,b),mass)/((mpi**2-SP(b,b))**2)*cutoff(SP(b,b))**4

      ! CHECKS:
      if (real(crossed).lt.0) then
         write(*,*) 'Crossed term less than zero!!', crossed
         stop
      else if (abs(aimag(crossed)).gt.1E-2) then
         write(*,*) 'Crossed term is imaginary!!', crossed
         stop
      else if (real(direct).lt.0) then
         write(*,*) 'Direct term less than zero!!', direct
         stop
      else if (abs(aimag(direct)).gt.1E-2) then
         write(*,*) 'Direct term is imaginary!!', crossed
         stop
      end if


      ! Interference contribution:
      lambda=0.
      do mu=0,3
         do nu=0,3
            lambda=lambda+spin32proj(delta,nu,mu,qp)*(b(nu)*a(mu))*metricTensor(mu,mu)*metricTensor(nu,nu)
         end do
      end do

      if (debug_ms) then
         inter_trace= trace( MatrixMult( (slashed(pp)-mN4),slashed(a),(slashed(p)+mN4),lambda,(slashed(q)+mN4),slashed(b))  )
      end if

      interference=1./(mpi**2-SP(b,b))/(mpi**2-SP(a,a))  &
           &  *trace( &
           &        MatrixMult( (slashed(pp)-mN4),slashed(a),(slashed(p)+mN4),lambda,(slashed(q)+mN4),slashed(b)) &
           &        ) &
           &  *cutoff(SP(a,a))**2  *  cutoff(SP(b,b))**2

      ! Factor 2 in front of crossed and direct terms due to Isospin
      ! Factor -4 in front of interference is -4=-2*2 where -2 is due to isospin!!
      ms_NN = 2.*(direct*z(mass,SP(a,a))+crossed*z(mass,SP(b,b)))-4.*real(interference)*sqrt(z(mass,SP(a,a))*z(mass,SP(b,b)))
      ms_NN = 1./4.*ms_NN*couplings*spectralDelta(mass)*pcm_final/pcm_initial*2.*mass

      if (debug_ms) then
            write(*,'(3G20.5)') interDimi(mass,SP(a,a),SP(b,b)), 2*inter_trace, 2*inter_trace/interDimi(mass,SP(a,a),SP(b,b))
            write(*,'(4G20.5)') direct_byHand(SP(a,a),mass), real(direct_traces), direct_byHand(SP(a,a),mass)/real(direct_traces),&
                 & direct_byHand_simple(SP(a,a),mass)
      end if


    end function ms_NN


    real function cutoff(t)
      ! Monopole form factor according to Dmitriev
      use constants, only: mPi
      real, intent(in) :: t                         ! Mandelstam t for the exchanged pion
      real, parameter :: lambdaSquared=0.63**2      ! = Dmitriev
      cutoff = (lambdaSquared-mPi**2) / (lambdaSquared-t)
    end function cutoff


    real function interDimi(mdel,t,u)
      ! Interference contribution to the matrix element squared according to Effenberger
      use constants, only: mN
      real, intent(in) :: mdel, u , t
      interDimi=1./(2.d0*mdel**2)*( &
           (t*u+(mdel**2-mn**2)*(t+u)-mdel**4+mn**4)* &
           (t*u+mn*(mn+mdel)*(mdel**2-mn**2)) - &
           1.d0/3.d0*(t*u-(mdel+mn)**2*(t+u)+(mn+mdel)**4)* &
           (t*u-mn*(mdel-mn)*(mdel**2-mn**2))      )*8*mn**2
    end function InterDimi


    real function direct_byHand(t,mass)
      ! Direct contribution to the matrix element squared: |M_D|^2
      use constants, only: mN
      real, intent(in) :: t, mass
      direct_byHand=-8.*mn**2*t*(-8./6.*((mn+mass)**2-t)  &
           &     *  (2.*t*(mn**2+mass**2)-t**2-(mn**2-mass**2)**2)/4./mass**2)
    end function direct_byHand


    real function direct_byHand_simple(t,mass)
      ! Direct contribution to the matrix element squared: |M_D|^2
      ! Same as direct by hand, but simplified algebra
      use constants, only: mN
      real, intent(in) :: t, mass
      direct_byHand_simple=-8*mn**2*t /(3* mass**2)* ( (mn+mass)**2-t )**2  * ( (mn-mass)**2-t )
    end function direct_byHand_simple

  end function MSquared_NN_ND



  real function spectralDelta(mass)
    ! Delta spectral function
    use constants, only: pi
    use particleProperties, only: hadron
    use idtable, only: delta
    use baryonWidth, only: fullWidthBaryon

    real, intent(in) :: mass
    real :: m0, gamma
    m0   =  hadron(delta)%mass

    gamma=  FullWidthBaryon(delta,mass) ! Use GiBUU width
    !gamma=  deltaWidth(mass)           ! Use Dmitriev width

    spectralDelta=1/pi* mass*gamma/((m0**2-mass**2)**2+mass**2*gamma**2)
  end function spectralDelta


!   real function deltaWidth(mass)
!     ! Delta width according to Dmitriev
!     use twoBodyTools, only : pcm
!     use particleProperties, only: hadron
!     use idtable, only : delta
!     use constants, only: mN, mPi
!     implicit none
!     real, parameter :: beta=0.2
!     real, intent(in) :: mass
!     real :: md,pcm_mass,pcm_pole!, gamma
!
!     md=hadron(delta)%mass
!     pcm_mass=pcm(mass,mN, mPi)
!     pcm_pole=pcm(md  ,mN, mPi)
!
!     ! Alexei paper:
!     !    deltaWidth=0.118*md/mass*(beta**2+pcm_pole**2)/(beta**2+pcm_mass**2)* &
!     !         & pcm_mass**3/pcm_pole**3
!     deltaWidth=0.120*(pcm_mass/pcm_pole)**3*z(mass,mPi)
!
!   end function deltaWidth


  real function z(mass,pionMass)
    ! (N Pi Delta) vertex form factor according to Dmitriev
    use particleProperties, only: hadron
    use Idtable, only: delta
    use TwoBodyTools, only:pcm
    use constants, only: mN

    real :: mass     ! Mass of delta
    real :: pionmass     ! Mass of pion

    real, parameter :: kappaSquared=0.2**2      ! = Dmitriev
    real :: pcm_pole, pcm_mass

    pcm_pole=pcm(hadron(delta)%mass, mN , pionMass)
    pcm_mass=pcm(mass              , mN , pionMass)

    z=(pcm_pole**2+kappaSquared)/(pcm_mass**2+kappaSquared)
  end function z



!   real function norm_integrand(mass)
!     ! Integrand for spectral function normalization
!     implicit none
!     real, intent(in) :: mass
!     !real :: m0, gamma
!     norm_integrand=spectralDelta(mass)*2.*mass
!   end function norm_integrand



  !****************************************************************************
  !****f* barBar_to_barBar_model/MSquared_ND_ND
  ! NAME
  ! real function MSquared_ND_ND(chargeIN,chargeOUT,srts,theta_in,mass_initial_in,mass_in,integrated_theta,integrated_mass,chargeSum)
  ! NOTES
  ! This function  evaluates :
  !   Matrix Element^2 * p_cm(final) / p_cm(initial) * SpecFunc_Delta
  ! for NDelta-> NDelta according to one-pion exchange. See Effenberger Diplom!!
  ! All variables in CM System!!
  ! INPUTS
  ! * integer, dimension (1:2) :: chargeIn,chargeOUT -- Charges of incoming and outgoing particles:
  !   first index nucleon, second:delta
  ! * real :: srts
  ! * real :: theta_in -- angle in CM system in radians
  ! * real ::  mass_in , mass_initial_in --  mass of final and initial Delta
  ! * logical :: integrated_theta, integrated_mass -- .true.= give integrated cross section
  ! * logical :: chargeSum -- .true.= sum over the final state charges
  !****************************************************************************
  real function MSquared_ND_ND(chargeIN,chargeOUT,srts,theta_in,mass_initial_in,mass_in,integrated_theta,integrated_mass,chargeSum)
    use constants, only: pi, mPi, mN

    integer, dimension (1:2) ,intent (in) :: chargeIn,chargeOUT
    real   , intent(in) :: srts
    real   , intent(in) :: theta_in, mass_in , mass_initial_in     ! angle in CM system in radians
    logical, intent(in) :: integrated_theta, integrated_mass       ! .true.= give integrated cross section
    logical, intent(in) :: chargeSum

    integer,parameter :: numsteps_theta=50
    integer,parameter :: numsteps_mass=50

    real :: cosTheta, dcosTheta,theta,mass, dmass
    integer :: i,j
    integer :: nucCharge, deltaCharge
    real    :: isoTotal

    if (firstTime) call readinput

    MSquared_ND_ND=0.
    dcosTheta=2./float(numsteps_theta)
    dmass=(srts-2*mN-mPi)/float(numsteps_mass)

    ! Integrate over mass
    massLoop: do j=1,numsteps_mass
       if (.not.integrated_mass) then
          mass=mass_in
       else
          mass=mN+mPi+float(j)*dmass
       end if
       ! Integrate over cosTheta
       thetaLoop: do i=1,numSteps_theta
          if (.not.integrated_theta) then
             theta=theta_in
          else
             cosTheta=-1.+float(i)*dcosTheta
             theta=acos(cosTheta)
          end if
          MSquared_ND_ND = MSquared_ND_ND + ms_ND(theta,srts,mass,mass_initial_in)
          ! If no theta integration shall be performed, then we exit the integration loop after the first cycle
          if (.not.integrated_theta) exit thetaLoop
       end do thetaLoop
       ! If no mass integration shall be performed, then we exit the integration loop after the first cycle
       if (.not.integrated_mass) exit massLoop
    end do massLoop

    if (integrated_mass)  MSquared_ND_ND=MSquared_ND_ND*dmass
    if (integrated_theta) MSquared_ND_ND=MSquared_ND_ND*2.*pi*dcosTheta

    if (chargeSum) then
       isoTotal=0.
       do nucCharge=0,1
          do deltaCharge=-1,2
             isoTotal=isoTotal+isospinFactor_ND_ND(chargeIn,(/nucCharge,deltaCharge/))
          end do
       end do
       MSquared_ND_ND=MSquared_ND_ND*isoTotal
    else
       MSquared_ND_ND=MSquared_ND_ND*isospinFactor_ND_ND(chargeIn,chargeOut)
    end if


  contains

    real function ms_ND(theta,srts,mass,massIni)
!       use particleProperties, only: hadron
!      use matrix_module, only : matrixMult,trace !, unit4
      use minkowski, only: SP !, slashed,gamma5,metricTensor
      use twoBodyTools, only: pcm
!      use IDTABLE, only : delta
      use constants, only: mN !,mPi
      use quadpack
      use spinProjector

      real, intent(in) :: theta, srts,mass,massIni

!       real :: md, couplings!,t,m1
      real, dimension (0:3) :: a,p,q,pp,qp  ! ,b
!       complex, dimension(0:3,0:3) :: mN4!,lambda_1,lambda_2
      real :: pcm_initial,pcm_final
!      complex :: traceResult
!      integer:: mu,nu
!      logical, parameter :: debug_ms=.false.

!       mD = hadron(delta)%mass

!       if(couplings_switch.eq.Pascalutsa) then
!          couplings=(h_A_deltaNuc/(2.*f_pi))**2*(g_A/(2.*f_pi))**2
!       else if(couplings_switch.eq.Dmitriev) then
!          couplings=(f_pi_deltadelta_doenges/mPi)**2    *   (f_pi_dimi/mPi)**2
!       else
!          stop'strange couplings in barbar_to_barbar_model'
!       end if

!       mN4=unit4*mN

      if (srts<=mN+max(mass,massIni)+epsilon(srts)) then
        ms_ND = 0.
        return
      end if

      !************************************************************************
      ! Setting up kinematics in CM-FRAME
      ! Initial particles
      pcm_initial= pCM(srts, mN , massIni )
      pcm_final  = pCM(srts, mN , mass )


      p(1:3)=(/0.,0.,pcm_initial/)
      q(1:3)=-p(1:3)

      ! Final ones
      pp(1:3)=(/sin(theta)*pcm_final,0.,cos(theta)*pcm_final/)
      qp(1:3)=-pp(1:3)

      ! Initial nucleon energy:
      p(0)  =sqrt( mN**2  + Dot_product (p(1:3) ,p(1:3)  ) )

      ! Initial delta energy:
      q(0)  =sqrt( massIni**2  + Dot_product (q(1:3) ,q(1:3)  ) )

      ! Final nucleon energy:
      pp(0) =sqrt( mN**2  + Dot_product (pp(1:3),pp(1:3) ) )
      ! Final delta energy :
      qp(0) =sqrt( mass**2+ Dot_product (qp(1:3),qp(1:3) ) )

      ! CHECKS:
      if (abs(p(0)+q(0)-pp(0)-qp(0)).gt.1E-5) then
         write(*,*) 'Energies not conserved',p(0)+q(0),pp(0)+qp(0)
         stop
      end if

      ! Exchanged pions
      a=p-pp
!       b=p-qp

      ms_ND = effenberger_MatrixElement(SP(a,a),massIni,mass)*spectralDelta(mass)*pcm_final/pcm_initial*2.*mass

!       if(debug_ms) then
!          ! Explicit calculation
!          ! With the following lines I tested that the result below gives
!          ! the same outcome as the line ms=... above. So the line ms=...
!          ! above is right!!!
!
!          traceResult=0.
!          do mu=0,3
!             do nu=0,3
!                traceResult=traceResult+ trace(MatrixMult(spin32proj(delta,mu,nu,qp),gamma5,&
!                     &                         slashed(a),spin32proj(delta,nu,mu,q),slashed(a),-gamma5)) &
!                     &                     *metricTensor(mu,mu)*metricTensor(nu,nu)
!             end do
!          end do
!          traceResult= traceResult &
!               &   * trace(MatrixMult(slashed(pp)-mn4,slashed(a),slashed(p)+mn4,slashed(a))) &
!               &   * (f_pi_deltadelta_doenges/mPi)**2    *   (f_pi_dimi/mPi)**2/8.    &
!               &   * cutoff(SP(a,a))**4/(SP(a,a)-mpi**2)**2
!
!          ms_ND = 1./8.*spectralDelta(mass)*pcm_final/pcm_initial*2.*mass*traceResult
!
!          ! Print test output:
!          write(*,'("(",2G18.4,")",7G18.4)') traceResult, &
!               & effenberger_MatrixElement(SP(a,a),massIni,mass),effenberger_MatrixElement(SP(a,a),massIni,mass)/real(traceresult)
!       end if

    end function ms_ND


    real function cutoff(t)
      ! Monopole form factor
      use constants, only: mPi
      real,intent(in) :: t     ! Mandelstam t for the exchanged pion
      cutoff = (lambda_cutoff**2-mPi**2) / (lambda_cutoff**2-t)
    end function cutoff


    real function effenberger_matrixElement(t,mu_i,mu_f)
      ! See Effenberger Diplomarbeit eq. (3.45)
      use constants, only: mPi,mN

      real,intent(in) :: mu_i,mu_f ! initial and final delta mass
      real,intent(in) :: t         ! Mandelstam "t"
      real :: couplings

      if (couplings_switch.eq.Pascalutsa) then
         couplings=(h_A_deltaNuc/(2.*f_pi))**2*(g_A/(2.*f_pi))**2
      else if (couplings_switch.eq.Dmitriev) then
         couplings=(f_pi_deltadelta_doenges/mPi)**2    *   (f_pi_dimi/mPi)**2
      else
         stop 'strange couplings in barbar_to_barbar_model'
      end if

      effenberger_matrixElement=1./8.*cutoff(t)**4/(t-mPi**2)**2 &
           & * (16.*(mu_i+mu_f)**2 * mN**2* t)/(9.*mu_i**2*mu_f**2) &
           & * (-mu_i**2+2.*mu_i*mu_f-mu_f**2+t)&
           & * (mu_i**4-2.*mu_i**3*mu_f+12.*mu_i**2*mu_f**2-2.*mu_i*mu_f**3&
           &    +mu_f**4-2.*mu_i**2 * t +2.*mu_i *mu_f* t-2.*mu_f**2* t+t**2) * couplings

    end function effenberger_matrixElement


  end function MSquared_ND_ND


  !****************************************************************************
  !****f* barBar_to_barBar_model/isospinFactor_ND_ND
  ! NAME
  ! real function isospinFactor_ND_ND(in_input, out_input )
  ! NOTES
  ! This function  evaluates isospin factors for ND-> ND scattering, cf. Table 3.6 of Effenberger's Diploma thesis
  ! or table (A.4) in Oliver's PhD.
  ! INPUTS
  ! Charge of incoming and outgoing particles (first index: nucleon, 2nd: Delta) :
  ! * integer, dimension(1:2), intent(in) :: in_input, out_input
  !****************************************************************************
   real function isospinFactor_ND_ND(in_input, out_input )

      ! charge of incoming and outgoing particles (first index: nucleon, 2nd: Delta) :
      integer, dimension(1:2), intent(in) :: in_input, out_input
      integer, dimension(1:2) :: in, out
      isospinFactor_ND_ND=0.
      if (in_input(1).eq.0) then
         ! Rotate all isospins such that we have an incoming proton:
         out(1)=1-out_input(1)
         out(2)=1-out_input(2)
         in(1) =1
         in(2) =1-in_input (2)
      else
         in=in_input
         out=out_input
      end if

      ! Isospin factors for incoming proton case:
      select case (in(2))
      case (2)
         if ((out(1).eq.1).and.out(2).eq.2) isospinFactor_ND_ND=9./4.
      case (1)
         if ((out(1).eq.0).and.out(2).eq.2) then
            isospinFactor_ND_ND=3
         else if ((out(1).eq.1).and.out(2).eq.1) then
            isospinFactor_ND_ND=1./4.
         end if
      case (0)
         if ((out(1).eq.1).and.out(2).eq.0) then
            isospinFactor_ND_ND=1./4.
         else if ((out(1).eq.0).and.out(2).eq.1) then
           isospinFactor_ND_ND=4
         end if
     case (-1)
         if ((out(1).eq.1).and.out(2).eq.-1) then
            isospinFactor_ND_ND=9./4.
         else if ((out(1).eq.0).and.out(2).eq.0) then
            isospinFactor_ND_ND=3.
         end if
      end select

    end function isospinFactor_ND_ND


end module barBar_to_barBar_model
