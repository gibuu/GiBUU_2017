!******************************************************************************
!****m* /comptonScattering
! NAME
! module comptonScattering
!
! PURPOSE
! This module implements Compton scattering cross sections
!
! INPUTS
! The Namelist "comptonScattering_nl" in the Jobcard.
!******************************************************************************


module comptonScattering

  implicit none
  private

  !****************************************************************************
  !****g* comptonScattering/debug
  ! SOURCE
  !
  logical, parameter :: debug=.false.
  ! PURPOSE
  ! Switch on debugging output
  !****************************************************************************


  logical, save :: initFlag=.true.

  public ::    compton_crossSection_dOmega


  integer, parameter :: lab=1
  integer, parameter :: cm=2

contains

!!$  !*******************************************************************
!!$  !****s* comptonScattering/readInput
!!$  ! NAME
!!$  ! subroutine readInput
!!$  ! PURPOSE
!!$  ! Reads input out of jobcard. Namelist 'comptonScattering_nl'.
!!$  !*******************************************************************
  subroutine readInput
!!$    use output
!!$
!!$    integer :: ios
!!$
!!$    !***************************************************************************
!!$    !****n* comptonScattering/comptonScattering_nl
!!$    ! NAME
!!$    ! NAMELIST /comptonScattering_nl/
!!$    ! PURPOSE
!!$    ! This Namelist for module "comptonScattering" includes:
!!$    ! * debug
!!$    !***************************************************************************
!!$    NAMELIST /comptonScattering_nl/ debug
!!$
!!$
!!$    call Write_ReadingInput('comptonScattering_nl',0)
!!$
!!$    rewind(5)
!!$    read(5,nml=comptonScattering_nl,IOSTAT=ios)
!!$    call Write_ReadingInput('comptonScattering_nl',0,ios)
!!$
!!$    write(*,'(A20,L4)') '  Debugging        ?',debug
!!$
!!$    call Write_ReadingInput('comptonScattering_nl',1)

  end subroutine readInput



  !****************************************************************************
  !****f* comptonScattering/compton_crossSection_dOmega
  ! NAME
  ! function compton_crossSection_dOmega(e_gamma,theta_gamma,charge_Nuc) Result (sigma)
  ! PURPOSE
  ! Evaluates Compton cross section assuming a free nucleon at rest, and a photon momentum along the z-axis.
  ! INPUTS
  ! * real                   , intent(in) :: e_gamma     ! Photon energy
  ! * real                   , intent(in) :: theta_gamma ! Photon scattering angle in degrees
  ! * integer                , intent(in) :: charge_nuc  ! Charge of Nucleon
  ! * integer                , intent(in) :: frame       ! 1=theta in lab, 2=theta in CM frame
  ! OUTPUT
  ! * real :: dsigma/dOmega    -- cross section dsigma_dOmega in mb/sr
  ! Note that the frame in which "dOmega" is calculated, is defined by the "frame" variable. E.g.: if
  ! frame=1 then dOmega=dOmega_lab or if frame=2 then dOmega=dOmega_CM
  !****************************************************************************
  function compton_crossSection_dOmega(e_gamma,theta_gamma,charge_Nuc,frame) Result (sigma)
    use degRad_conversion, only: radian

    real :: sigma

    integer                , intent(in) :: charge_Nuc  ! Charge of Nucleon
    real                   , intent(in) :: e_gamma     ! Photon energy in lab frame
    real                   , intent(in) :: theta_gamma ! Photon scattering angle
    integer                , intent(in) :: frame       ! 1=theta defined in lab, 2=theta defined in CM frame

    real, dimension(0:3)                :: p_in,p_out     ! Incoming and outgoing momentum of nucleon
    real, dimension(0:3)                :: q_in,q_out     ! Incoming and outgoing momentum of photon

    real    :: theta
    logical :: success

    if (initFlag) then
       call readinput()
       initFlag=.false.
    end if

    ! (1) Establish kinematics
    theta=radian(theta_gamma)
    call getKinematics(p_in,p_out,q_in,q_out,e_gamma,theta,frame)

    if (maxVal(abs(p_in-p_out+q_in-q_out)).gt.0.00001) then
       write(*,*) 'No energy conservation! STOP'
       STOP
    end if

    ! (2) Evaluate cross section
    sigma=getSigma(p_in,p_out,q_in,q_out,charge_Nuc,success)

  end function compton_crossSection_dOmega


  !****************************************************************************
  !****is* comptonScattering/getKinematics
  ! NAME
  ! subroutine getKinematics(p_in,p_out,q_in,q_out,e_gamma,theta_gamma,frame)
  ! PURPOSE
  ! Evaluates kinematics for compton scattering
  ! INPUTS
  ! * real                   , intent(in) :: e_gamma     ! Photon energy
  ! * real                   , intent(in) :: theta_gamma ! Photon scattering angle in radians
  ! * integer                , intent(in) :: frame       ! 1=theta in lab, 2=theta in CM frame
  ! OUTPUT
  ! * real, dimension(0:3), intent(out) :: p_in,p_out,q_in,q_out
  ! NOTE
  ! *  Same notations as in compton_crossSection_dOmega
  !****************************************************************************
  subroutine getKinematics(p_in,p_out,q_in,q_out,e_gamma,theta_gamma,frame)
    use twobodyTools, only: pCM
    use minkowski, only: abs4,SP
    use vector, only: theta_in
    use constants, only: mN

    real                   , intent(in) :: e_gamma     ! Photon energy in lab frame
    real                   , intent(in) :: theta_gamma ! Photon scattering angle in radians
    integer                , intent(in) :: frame       ! 1=theta defined in lab, 2=theta defined in CM frame
    real, dimension(0:3), intent(out) :: p_in,p_out,q_in,q_out

    real    :: p_cm,q_out_abs,W

    select case (frame)
    case (lab)
       ! Incoming Momenta
       p_in=(/mN     ,0.,0.,0.     /)
       q_in=(/e_gamma,0.,0.,e_gamma/)
       ! Define Outgoing Momenta
       q_out_abs=e_gamma*mN/(e_gamma*(1.-cos(theta_gamma))+mN)
       q_out=(/1.,sin(theta_gamma),0.,cos(theta_gamma)/)*q_out_abs
    case (cm)
       ! Incoming Momenta
       W=sqrt(mN**2+2.*mN*e_gamma) ! CM energy
       p_cm= pCM(W, mN, 0.)
       p_in=(/sqrt(mN**2+p_cm**2),0.,0.,-p_cm/)
       q_in=(/p_cm               ,0.,0., p_cm/)
       ! Define Outgoing Momenta
       q_out=(/1.,sin(theta_gamma),0.,cos(theta_gamma)/)*p_cm
    case default
       write(*,*) 'Wrong frame in  compton_crossSection! frame=', frame
       stop
    end select

    p_out=p_in+q_in-q_out

    ! Check Kinematics
    if (abs(abs4(p_out)-mN).gt.0.001) then
       write(*,*) 'Error in compton_crossSection!! Wrong mass of outgoing nucleon!!',abs4(p_out)-mN
       write(*,'(3F15.4)') abs4(P_out), SP(q_out,q_out), theta_in(q_in(1:3),q_out(1:3),1)
       stop
    else if (p_out(0).lt.0) then
       write(*,*) 'Error in compton_crossSection!! Wrong energy of outgoing nucleon!!',p_out
       stop
    end if
    if (debug) then
       write(*,'(A,4F15.4)') 'in: ',p_in+q_in
       write(*,'(A,4F15.4)') 'out:',p_out+q_out
       write(*,'(3F15.4)') abs4(P_out), SP(q_out,q_out), theta_in(q_in(1:3),q_out(1:3),1)
       write(*,'(4F15.4)') P_out
       write(*,'(4F15.4)') q_out
    end if

  end subroutine getKinematics


  !****************************************************************************
  !****if* comptonScattering/getSigma
  ! NAME
  ! function getSigma(p_in,p_out,q_in,q_out,theta_gamma,success) result(sigma)
  ! PURPOSE
  ! Evaluates cross section for compton scattering
  ! INPUTS
  ! * real, dimension(0:3), intent(in) :: p_in,p_out,q_in,q_out
  ! OUTPUT
  ! *  logical, intent(out)             :: success
  ! *  real                             :: sigma ! dsigma/dOmega (mb/sr)

  ! NOTE
  ! *  Same notations as in compton_crossSection_dOmega
  !****************************************************************************
  function getSigma(p_in,p_out,q_in,q_out,charge_Nuc,success) result(sigma)
    use vector, only: absVec
    use constants, only: GeVSquared_times_mb,pi

    real, dimension(0:3), intent(in)  :: p_in,p_out,q_in,q_out
    integer             , intent(in)  :: charge_Nuc
    logical             , intent(out) :: success

    real                              :: sigma
    real                              :: matrixElementSquared

    success=.false.
    sigma=0.


    ! Evaluate sigma

    matrixElementSquared=getMatrixElement()

    ! dsigma/dOmega(q_out)=...
    sigma=matrixElementSquared/(64.*pi**2)*absVec(p_out(1:3)) &
                          & /absVec(p_in(0)*q_in(1:3)-p_in(1:3)*q_in(0))/(p_out(0)+q_out(0))

    ! Convert Units from GeV**-2/sr to mb/sr:
    sigma=sigma/GeVSquared_times_mb

  contains

      real function getMatrixElement()
        use matrix_module, only: matrixMult, trace,unit4,printMatrix
        use minkowski, only: tilde, slashed, metricTensor
        use IdTable, only: Delta_ID => Delta
        use vertices, only: vertex_gammaN_to_spin3_2
        use propagators, only: propagator_3_2_vac
        use constants, only: electronChargeSQ, mN

        complex,dimension(0:3,0:3) :: matrix,proj_out,proj_in

        integer :: alpha,beta,gamma,delta,mu,nu,rho,sigma

        integer :: resID,resCharge

        complex,dimension(0:3,0:3,0:3,0:3) :: vertex_out, vertex_out_tilde ! Outgoing vertices
        complex,dimension(0:3,0:3,0:3,0:3) :: vertex_in , vertex_in_tilde  ! Incoming vertices
        complex,dimension(0:3,0:3,0:3,0:3) :: prop_low  , prop_low_tilde   ! Propagators
        real, dimension(0:3)               :: p_res                        ! Resonance momentum

        proj_out=slashed(p_out)+unit4*mN
        proj_in =slashed(p_in) +unit4*mN

        resID     = Delta_ID
        resCharge = charge_Nuc
        p_res     = p_in+q_in

        !p_out_inv=(/p_out(0),-p_out(1:3)/)
        !q_out_inv=(/q_out(0),-q_out(1:3)/)
        !vertex_out=vertex_gammaN_to_spin3_2(q_out_inv ,p_out_inv ,resID, resCharge)
        vertex_out=vertex_gammaN_to_spin3_2(q_out ,p_out ,resID, resCharge)
        vertex_in =vertex_gammaN_to_spin3_2(q_in  ,p_in  ,resID, resCharge)
        prop_low  =propagator_3_2_vac(resID,p_res)


        do alpha=0,3
           do beta=0,3
              vertex_out_tilde(alpha,beta,:,:) =tilde(vertex_out(alpha,beta,:,:))
              vertex_in_tilde (alpha,beta,:,:) =tilde( vertex_in(alpha,beta,:,:))
              ! Propagators:
              prop_low        (alpha,beta,:,:) =prop_low        (alpha,beta,:,:)* metricTensor(alpha,alpha)* metricTensor(beta,beta)
              prop_low_tilde  (alpha,beta,:,:) =tilde(  prop_low(alpha,beta,:,:))
           end do
        end do


        matrix=0.
        do alpha=0,3
           do beta=0,3
              do gamma=0,3
                 do delta=0,3
                    !if(delta.ne.alpha) cycle
                    do mu=0,3
                       do nu=0,3
                          do rho=0,3
                             !if(rho.ne.mu) cycle
                             do sigma=0,3
                                matrix=matrix+MatrixMult(proj_out,vertex_out_tilde(beta,alpha,:,:),prop_low(beta,nu,:,:),&
                                     &                   vertex_in(nu,mu,:,:),proj_in,vertex_in_tilde(sigma,rho,:,:), &
                                     &                   prop_low_tilde(gamma,sigma,:,:),vertex_out(gamma,delta,:,:) ) &
                                     &                   *metricTensor(alpha,delta)*metricTensor(mu,rho)
                             end do
                          end do
                       end do
                    end do
                 end do
              end do
           end do
        end do


        ! Attention: Factor (1/(2mN))^2 is omitted
        getMatrixElement=1./4.*trace(matrix)*electronChargeSQ**2

        if (getMatrixElement.lt.0) then
           write(*,*) '|M|^2 less than zero!', getMatrixElement
           call printMatrix(matrix)
           stop 'comptonScattering.f90'
        end if

      end function getMatrixElement


  end function getSigma


end module comptonScattering
