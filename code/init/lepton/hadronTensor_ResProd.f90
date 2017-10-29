!******************************************************************************
!****m* /hadronTensor_ResProd
! NAME
! module hadronTensor_ResProd
!
! PURPOSE
! This module defines the hadronic tensors for resonance production
! by EM- and weak currents
!
!******************************************************************************
module hadronTensor_ResProd

  implicit none
  private

  public:: hadronTensor_R

  !****************************************************************************
  !****g* hadronTensor_ResProd/debugFlag
  ! SOURCE
  logical, parameter :: debugFlag=.false.
  ! PURPOSE
  !
  !****************************************************************************

  !****************************************************************************
  !****g* hadronTensor_ResProd/speedup
  ! SOURCE
  logical,save      :: speedup=.true.
  ! PURPOSE
  !
  !****************************************************************************

  logical,save      :: initFlag=.true.

contains

  subroutine initInput
    use output

    integer :: ios
    !**************************************************************************
    !****n* hadronTensor_ResProd/HadronTensor_ResProd
    ! NAME
    ! NAMELIST /HadronTensor_ResProd/
    ! PURPOSE
    ! Namelist for module hadronTensor_ResProd includes:
    ! * speedup
    !**************************************************************************
    NAMELIST /HadronTensor_ResProd/ speedup

    call Write_ReadingInput('hadronTensor_ResProd',0)
    rewind(5)
    read(5,nml=hadronTensor_ResProd,IOSTAT=ios)
    call Write_ReadingInput("hadronTensor_ResProd",0,ios)

    write(*,*) 'speedup?  ', speedup

    call Write_ReadingInput('hadronTensor_ResProd',1)

  end subroutine initInput

  !****************************************************************************
  !****f* hadronTensor_ResProd/hadronTensor_R
  ! NAME
  ! function hadronTensor_R(pi,pf,resID,targetCharge,process,matrix)   &
  !  result(success)
  !
  ! PURPOSE
  ! This function returns the hadronic tensor for resonance production of the 
  ! resonance "resID" by a current type defined by "process" (e.g. Gamma* N -> R).
  !
  ! INPUTS
  ! * real, dimension(0:3),intent(in) :: pi            -- 4-momentum of incoming nucleon
  ! * real, dimension(0:3),intent(in) :: pf            -- 4-momentum of outgoing resonance
  ! * integer,intent(in)              :: resID         -- ID of resonance
  ! * integer,intent(in)              :: process       -- type of the process (EM, CC, NC)
  ! * integer,intent(in)              :: targetCharge  -- Charge of target nucleon
  !
  ! OUTPUT
  ! * complex, dimension(0:3,0:3) :: matrix -- H^mu nu
  ! * logical, intent(out)        :: success  -- true if resonance has coupling strength to the process;
  ! *                                            false if there is no coupling, therefore H^mu nu=0
  ! NOTES
  ! * We use the Peskin notation where the spin sum(u(p) ubar(p))=slashed(p)+m    !
  !****************************************************************************
  function hadronTensor_R(pi,pf,resID,targetCharge,process,matrix,bare_mass) &
                         &  result(success)
    use particleProperties, only: hadron
    use formFactor_ResProd, only: getFormfactor_Res
    use minkowski, only: SP
    use leptonicID

    real, dimension(0:3)       , intent(in)      :: pi, pf
    integer                    , intent(in)      :: resID, process, targetCharge
    complex, dimension(0:3,0:3), intent(inout)   :: matrix
    real                       , intent(in)      :: bare_mass


    logical              :: success
    integer              :: processID
    real, dimension(1:8) :: formfactor
    integer              :: parity,spinTimes2

    success=.false.
    matrix=0.

    if (initFlag) then
       call initInput
       initFlag=.false.
    end if


    if ( debugflag ) then
       write(*,*) -SP(pf-pi,pf-pi)
       write(*,*) sp(pf,pf)
    end if
    if (sp(pf,pf).le.0) then
       return
    end if


    processID=process

    if (process.eq.99) processID=CC !workaround needed for background contribution

    formfactor=getFormfactor_Res(-SP(pf-pi,pf-pi),bare_mass,resID,   &
                                & targetCharge,processID,success)
    if (.not.success) then
       matrix=0.
       return
    end if

    if (process.eq.99) processID=EM !workaround needed for background contribution


    parity=(-1)**(1+hadron(resID)%AngularMomentum)
    spinTimes2=NINT(hadron(resID)%spin*2.)

    if ( debugflag ) then
       write(*,'(8E10.2)') formfactor
       write(*,*) '2.*Spin=',spinTimes2
       write(*,*) 'parity=',parity
    end if

    select case (spinTimes2)
    case (1)
       ! e.g. S11,P11
       matrix=hadronTensor_1_2(pi,pf,parity,formfactor(1:2),formfactor(5:6),  &
                              & resID,processID)
    case (3)
       ! e.g. D13,P33
       matrix=hadronTensor_3_2(pi,pf,parity,formfactor,resID,processID)
    case (5,7)
       ! As an approximation, we use the spin 3/2 structure:
       matrix=hadronTensor_3_2(pi,pf,parity,formfactor,resID,processID)
    case default
       call writeError(resID)
       return
    end select

    success=.true.
  contains
    !**************************************************************************
    !****is* hadronTensor_R/writeError
    ! PURPOSE
    ! Error Message
    !**************************************************************************
    subroutine writeError(resID)
      use particleProperties, only: hadron

      integer,intent(in) :: resID
      write(*,*) 'This case is not yet implemented!'
      write(*,*) 'Isospin =',hadron(resID)%isospinTimes2 /2.
      write(*,*) 'Spin    =',hadron(resID)%spin
    end subroutine writeError

  end function hadronTensor_R


  !****************************************************************************
  !****************************************************************************
  !
  !  SSSSSSS    PPPPPPP    II   NN    NN                11          /     22222
  !  SSS        P     P    II   NNNN  NN    =====      111         /          22
  !  SSSSSSS    PPPPPPP    II   NN NN NN    =====     1111        /         22
  !     SSSS    P          II   NN  NNNN                11       /        22
  !  SSSSSSS    P          II   NN    NN               111      /         222222
  !
  !****************************************************************************
  !****************************************************************************


  !****************************************************************************
  !****f* hadronTensor_ResProd/hadronTensor_1_2
  ! NAME
  ! function hadronTensor_1_2(pi,pf,parity,G,GA,resID) result(matrix)
  !
  ! PURPOSE
  ! This function returns the hadronic tensor for spin=1/2 resonances.
  !
  ! INPUTS
  ! *  real, dimension(1:2),intent(in) :: F             -- Vector or EM-Form factors
  ! *  real, dimension(1:2),intent(in) :: FA            -- Axial form factors
  ! *  real, dimension(0:3),intent(in) :: pi            -- 4-momentum of incoming nucleon
  ! *  real, dimension(0:3),intent(in) :: pf            -- 4-momentum of outgoing resonance
  ! *  integer,intent(in)              :: resID         -- ID of resonance
  ! *  integer,intent(in)              :: parity        -- Parity of resonance
  ! *  integer,intent(in)              :: process       -- type of the process (EM, CC, NC)
  !
  ! OUTPUT
  ! * complex, dimension(0:3,0:3) :: matrix -- H^mu nu
  !****************************************************************************
  function hadronTensor_1_2(pi,pf,parity,G,GA,resID,process) result(matrix)
    use leptonicID, only: EM
    use minkowski, only: abs4, gamma0, slashed
    use matrix_module, only: MatrixMult, dagger, trace, unit4

    integer, intent(in)            :: parity,process
    integer, intent(in)            :: resID
    real, dimension(0:3),intent(in):: pi,pf
    real, dimension(1:2),intent(in):: G   ! Vector or EM-Form factors
    real, dimension(1:2),intent(in):: GA  ! Axial form factors
    complex, dimension(0:3,0:3)    :: matrix

    complex, dimension(0:3,0:3) :: projector_in,projector_out,j_mu,j_nu, &
                                 & j_dagger_nu
    real                        :: mi,mf
    integer                     :: mu,nu
!    logical, parameter          :: speedup=.true.
    logical                     :: do_Once


    do_Once=.true.

    if ( debugflag ) then
       write(*,*) 'Spin 1/2'
    end if

    mi=abs4(pi)
    mf=abs4(pf)

    projector_out=slashed(pf)+mf*unit4

    if (speedup.and.(parity.eq.1)) then
       projector_in=slashed(pi)-mi*unit4
    else
       projector_in=slashed(pi)+mi*unit4
    end if

    matrix=0.
    do mu=0,3
       if (speedup) j_mu=j_neg(mu)
       do nu=0,3
          select case (parity)
          case (1)
             if (speedup) then
                j_dagger_nu=j_dagger(nu)
             else
                j_mu=j_pos(mu)
                j_nu=j_pos(nu)
             end if
          case (-1)
             if (speedup) then
                j_dagger_nu=j_dagger(nu)
             else
                j_mu=j_neg(mu)
                j_nu=j_neg(nu)
             end if
          case default
             call parityError(parity,'11')
          end select

          if (speedup) then
             matrix(mu,nu)=1./2.* &
                  & Trace(MatrixMult(projector_out,j_mu,projector_in,j_dagger_nu))
          else
             matrix(mu,nu)=1./2.* &
                  & Trace(MatrixMult(projector_out,j_mu,projector_in,gamma0,  & 
                  & dagger(j_nu),gamma0))
          end if

       end do
    end do

  contains

    !**************************************************************************
    !****if* hadronTensor_1_2/j_pos
    ! NAME
    ! function j_pos(mu) result(matrix)
    !
    ! PURPOSE
    ! This function the hadronic flux "J^mu" for resonances of positive parity
    !
    ! INPUTS
    ! integer, intent(in) :: mu
    !
    ! OUTPUT
    ! * complex, dimension(0:3,0:3) :: matrix
    !**************************************************************************
    function j_pos(mu) result(matrix)
      use minkowski, only: gamma5

      integer, intent(in) :: mu
      complex, dimension(0:3,0:3) :: matrix

      matrix=MatMul(j_neg(mu),gamma5)

    end function j_pos

    !**************************************************************************
    !****if* hadronTensor_1_2/j_neg
    ! NAME
    ! function j_neg(mu) result(matrix)
    !
    ! PURPOSE
    ! This function the hadronic flux "J^mu" for resonances of negative parity
    !
    ! INPUTS
    ! integer, intent(in) :: mu
    !
    ! OUTPUT
    ! * complex, dimension(0:3,0:3) :: matrix
    !**************************************************************************
    function j_neg(mu) result(matrix)
      use minkowski, only: gamma, gamma5, SP, slashed, sigma4, metricTensor
      use matrix_module, only: unit4
      use constants, only: mN, ii

      integer, intent(in) :: mu
      complex, dimension(0:3,0:3) :: matrix,sigma4_q
      real, dimension(0:3),save        :: q
      integer :: rho,alpha
      real ,save   :: mass_mu, QSquared,G1_mass_mu2,G2_mass_mu
      complex, save,dimension(0:3,0:3) :: slashed_q


      if (do_once.or.(.not.speedup)) then
         q= pf-pi
         ! This definition differs in comparison to Lalakulich, 
         ! is however consistent with Alvarez-Ruso
         mass_mu=2*mN
         QSquared=-SP(q,q)
         if (speedup) then
            slashed_q=slashed(q)
            G1_mass_mu2=G(1)/mass_mu**2
            G2_mass_mu=G(2)/mass_mu
         end if
         do_once=.false.
      end if

      sigma4_q=0
      do rho=0,3
         do alpha=0,3
            if (alpha.ne.rho) cycle
            sigma4_q=sigma4_q+sigma4(mu,rho)*q(alpha)*metricTensor(rho,alpha)
         end do
      end do

      if (speedup) then
         matrix= MatMul(  G1_mass_mu2 * (QSquared*gamma(:,:,mu)+q(mu)*slashed_q) &
                       & +  G2_mass_mu  * ii* sigma4_q ,gamma5)
      else
         matrix= MatMul(  G(1)/mass_mu**2 *    &
              & (QSquared*gamma(:,:,mu)+q(mu)*slashed(q)) &
              & +  G(2)/mass_mu              * ii* sigma4_q         ,gamma5)
      end if

      if (process.ne.EM) then
         matrix= matrix + GA(1)                      * gamma(:,:,mu) &
              & + GA(2)/mN * q(mu)*unit4
      end if
    end function j_neg



    !**************************************************************************
    !****if* hadronTensor_1_2/j_dagger
    ! NAME
    ! function j_neg_dagger(mu) result(matrix)
    !
    ! PURPOSE
    ! This function returns "gamma_0 (J^mu)^dagger gamma_0" for resonances of 
    ! negative parity
    !
    ! INPUTS
    ! integer, intent(in) :: mu
    !
    ! OUTPUT
    ! * complex, dimension(0:3,0:3) :: matrix
    !**************************************************************************
    function j_dagger(mu) result(matrix)
      use minkowski, only: gamma, gamma5, SP, slashed, sigma4, metricTensor
      use matrix_module, only: unit4
      use constants, only: mN, ii

      integer, intent(in) :: mu
      complex, dimension(0:3,0:3) :: matrix,sigma4_q
      real, dimension(0:3)        :: q
      integer :: rho,alpha
      real    :: mass_mu, QSquared

      q= pf-pi
      mass_mu=2.*mN
      QSquared=-SP(q,q)

      sigma4_q=0
      do rho=0,3
         do alpha=0,3
            if (alpha.ne.rho) cycle
            sigma4_q=sigma4_q+sigma4(mu,rho)*q(alpha)*metricTensor(rho,alpha)
         end do
      end do

      matrix= MatMul(gamma5, - G(1)/mass_mu**2 *   &
           & (QSquared*gamma(:,:,mu)+q(mu)*slashed(q)) &
           & +  G(2)/mass_mu              * ii* sigma4_q        )

      if (process.ne.EM) then
         matrix= matrix + GA(1)                      * gamma(:,:,mu) &
              & + GA(2)/mN * q(mu)*unit4
      end if
    end function j_dagger





  end function hadronTensor_1_2


  !****************************************************************************
  !****************************************************************************
  !
  !  SSSSSSS    PPPPPPP    II   NN    NN            333333          /     22222
  !  SSS        P     P    II   NNNN  NN    =====       33         /          22
  !  SSSSSSS    PPPPPPP    II   NN NN NN    =====   333333        /         22
  !     SSSS    P          II   NN  NNNN                33       /        22
  !  SSSSSSS    P          II   NN    NN            333333      /         222222
  !
  !****************************************************************************
  !****************************************************************************


  !****************************************************************************
  !****f* hadronTensor_ResProd/hadronTensor_3_2
  ! NAME
  ! function   hadronTensor_3_2(pi,pf,parity,formfactor,resID) result(matrix)
  !
  ! PURPOSE
  ! This function returns the hadronic tensor for spin=3/2 resonances.
  !
  ! INPUTS
  !
  ! *  real, dimension(0:3),intent(in) :: pi            -- 4-momentum of incoming nucleon
  ! *  real, dimension(0:3),intent(in) :: pf            -- 4-momentum of outgoing resonance
  ! *  integer,intent(in)              :: parity        -- Parity of resonance
  ! *  real, dimension(1:8),intent(in) :: FormFactor    -- 1-4: Vector or EM, 5-8: Axial form factors
  ! *  integer,intent(in)              :: resID         -- ID of resonance
  ! *  integer,intent(in)              :: process       -- type of the process (EM, CC, NC)
  ! OUTPUT
  ! * complex, dimension(0:3,0:3) :: matrix -- H^mu nu
  !****************************************************************************
  function   hadronTensor_3_2(pi,pf,parity,F,resID,process) result(matrix)
    use leptonicID, only: EM
    use spinProjector, only: spin32proj,spin32proj_tensor
    use minkowski, only: abs4, gamma0, slashed, metricTensor
    use matrix_module, only: unit4, MatrixMult, dagger,trace,printMatrix

    integer, intent(in)            :: parity,process
    real, dimension(0:3),intent(in):: pi,pf
    integer,intent(in)             :: resID         ! ID of resonance
    real, dimension(1:8),intent(in):: F   ! 1-4: Vector, 5-8: Axial form factors
    complex, dimension(0:3,0:3)    :: matrix

    complex, dimension(0:3,0:3) :: projector_in,projector_out,j_mu,j_nu,   &
                                 & matrix_toTrace,j_nu_dagger
    integer                     :: mu,nu,alpha,beta,delta,kappa

    real                        :: mi
    logical      :: do_Once
    logical      :: do_Once2

    complex, dimension(0:3,0:3,0:3,0:3)    :: spin_proj_3_2


    do_Once=.true.
    do_Once2=.true.

    mi=abs4(pi)
    if ( debugflag ) then
       write(*,*) 'Spin 3/2'
       write(*,*) 'Parity=',parity
       write(*,*) 'resID=',resID
       write(*,*) 'mi=',mi
    end if

    if (speedup.and.(parity.eq.1)) then
       projector_in=slashed(pi)-mi*unit4
    else
       projector_in=slashed(pi)+mi*unit4
    end if

    if (speedup) spin_proj_3_2=spin32proj_tensor(resID,pf)

    matrix=0.
    if (debugFlag) call printMatrix(matrix_toTrace)
    do mu=0,3
       do nu=0,3
          Matrix_toTrace=0.
          alphaSum : do alpha=0,3
             betaSum : do beta=0,3
                if (speedup) then
                   projector_out=spin_proj_3_2(alpha,beta,:,:)
                else
                   projector_out=spin32proj(resID,alpha,beta,pf)
                end if
                kappaSum : do kappa=0,3
                   if (kappa.ne.alpha) cycle
                   if (speedup) j_nu_dagger=j_dagger(kappa,nu)
                   deltaSum : do delta=0,3
                      if (delta.ne.beta) cycle
                      if (speedup) then
                         j_mu=j_neg(delta,mu)
                         matrix_toTrace=matrix_toTrace+     &
                         & metricTensor(beta,delta)*metricTensor(alpha,kappa)* &
                         & MatrixMult(projector_out,j_mu,projector_in,j_nu_dagger)
                      else
                         select case (parity)
                         case (1)
                            j_mu=j_pos(delta,mu)
                            j_nu=j_pos(kappa,nu)
                         case (-1)
                            j_mu=j_neg(delta,mu)
                            j_nu=j_neg(kappa,nu)
                         case default
                            call parityError(parity,'11')
                         end select
                         matrix_toTrace=matrix_toTrace+    & 
                           &metricTensor(beta,delta)*metricTensor(alpha,kappa)*&
                         & MatrixMult(projector_out,j_mu,projector_in,gamma0,  &
                         & dagger(j_nu),gamma0)
                      end if
                   end do deltaSum
                end do kappaSum
             end do betaSum
          end do alphaSum
          matrix(mu,nu)=1./2.*Trace(matrix_toTrace)
       end do
    end do


  contains

    !**************************************************************************
    !****if* hadronTensor_3_2/j_pos
    ! NAME
    ! function j_pos(lambda,nu) result(matrix)
    !
    ! PURPOSE
    ! * This function evaluates the hadronic flux "J^lambda nu" for resonances 
    ! * of positive parity
    ! * Same as negative parity case times gamma(5).
    !
    ! INPUTS
    ! * integer, intent(in) :: lambda,nu  -- lambda is the lorentz index
    !   which is contracted with the 3/2 spinor
    !
    ! OUTPUT
    ! * complex, dimension(0:3,0:3) :: matrix
    !**************************************************************************
    function j_pos(lambda,nu) result(matrix)
      use minkowski, only: gamma5

      integer, intent(in)         :: lambda,nu
      complex, dimension(0:3,0:3) :: matrix

      matrix=MatMul(j_neg(lambda,nu),gamma5)

    end function j_pos

    !**************************************************************************
    !****if* hadronTensor_3_2/j_neg
    ! NAME
    ! function j_neg(alpha,nu) result(matrix)
    !
    ! PURPOSE
    ! * This function evaluates the hadronic flux "J^lambda nu" for resonances 
    !   of negative parity
    ! * See Phys Rev D74, 014009 (2006), eq. 4.3
    !
    ! INPUTS
    ! * integer, intent(in) :: alpha,nu -- lambda is the lorentz index 
    !   which is contracted with the 3/2 spinor
    !
    ! OUTPUT
    ! * complex, dimension(0:3,0:3) :: matrix
    !**************************************************************************
    function j_neg(lambda,nu) result(matrix)
      use minkowski, only: gamma, gamma5, SP, slashed, metricTensor
      use matrix_module, only: unit4
      use constants, only: mN

      integer, intent(in) ::lambda,nu
      complex, dimension(0:3,0:3) :: matrix
      real, dimension(0:3),save        :: q
      !real ,save                       :: QSquared
      complex, dimension(0:3,0:3),save :: matrix_1,matrix_2!,matrix_3

      real ,save :: F1_mN,F2_mN2,F3_mN2,F5_mN,F6_mN2,F8_mN2



      if (do_once2.or.(.not.speedup)) then
         q= pf-pi
         !QSquared=-SP(q,q)
         !mf=abs4(pf)
         if (speedup) then
            matrix_1= F(1)/mN * slashed(q) + (F(2)/mN**2 * SP(pf,q)   &
                      & + F(3)/mN**2 * SP(pi,q) + F(4))* unit4
            F1_mN=F(1)/mN
            F2_mN2=F(2)/mN**2
            F3_mN2=F(3)/mN**2
            if (process.ne.EM) then
               F5_mN=F(5)/mN
               F6_mN2=F(6)/mN**2
               F8_mN2=F(8)/mN**2
               matrix_2= F(5)/mN * MatMul(slashed(q),gamma5)          &
                         & + F(6)/mN**2* SP(pf,q) * gamma5
            end if
         end if
         do_Once2=.false.
      end if

      if (speedup) then
         matrix=  metricTensor(lambda,nu)* matrix_1 &
              & - q(lambda)*                ( F1_mN * gamma(:,:,nu)   &
              & + (F2_mN2 * pf(nu)  + F3_mN2 * pi(nu)       )* unit4 )
      else
         matrix=  metricTensor(lambda,nu)* ( F(1)/mN * slashed(q)     &
              & + (F(2)/mN**2 * SP(pf,q)&
              & + F(3)/mN**2 * SP(pi,q) + F(4))* unit4 ) &
              & - q(lambda)*               ( F(1)/mN * gamma(:,:,nu)  + (F(2)/mN**2 * pf(nu)  &
              & + F(3)/mN**2 * pi(nu)         )* unit4 )
      end if

     if (process.ne.EM) then
        if (speedup) then
           matrix=matrix &
                & + metricTensor(lambda,nu)* matrix_2  &
                & - q(lambda)* ( F5_mN * MatMul(gamma(:,:,nu) ,gamma5) &
                & + F6_mN2 * pf(nu)  * gamma5  )   &
                & + (metricTensor(lambda,nu)* F(7)+q(lambda)*q(nu)*F8_mN2)*gamma5
         else
            matrix=matrix &
            &  + metricTensor(lambda,nu)   &
            & * ( F(5)/mN * MatMul(slashed(q),gamma5) + F(6)/mN**2* SP(pf,q) &
            &  * gamma5 ) &
            & - q(lambda)* ( F(5)/mN * MatMul(gamma(:,:,nu) ,gamma5)   &
            & + F(6)/mN**2 * pf(nu) * gamma5 )&
            & + (metricTensor(lambda,nu)* F(7)+q(lambda)*q(nu)*F(8)/mN**2)*gamma5
         end if
     end if
    end function j_neg


    !**************************************************************************
    !****if* hadronTensor_3_2/j_dagger
    ! NAME
    ! function j_dagger(alpha,nu) result(matrix)
    !
    ! PURPOSE
    ! * This function evaluates "gamma_o (J^lambda nu)^dagger gamm_o" 
    !   for resonances of negative parity
    ! * See Phys Rev D74, 014009 (2006), eq. 4.3
    !
    ! INPUTS
    ! * integer, intent(in) :: alpha,nu -- lambda is the lorentz index
    !   which is contracted with the 3/2 spinor
    !
    ! OUTPUT
    ! * complex, dimension(0:3,0:3) :: matrix
    !**************************************************************************
    function j_dagger(lambda,nu) result(matrix)
      use minkowski, only: gamma, gamma5, SP, slashed, metricTensor
      use matrix_module, only: unit4
      use constants, only: mN

      integer, intent(in) ::lambda,nu
      complex, dimension(0:3,0:3) :: matrix
      real, dimension(0:3),save        :: q
      !real ,save                       :: QSquared
      complex, dimension(0:3,0:3),save :: matrix_1,matrix_2

      real ,save :: F1_mN,F2_mN2,F3_mN2,F5_mN,F6_mN2,F8_mN2


      if (do_once.or.(.not.speedup)) then
         q= pf-pi
         !QSquared=-SP(q,q)
         !mf=abs4(pf)
         if (speedup) then
            matrix_1= F(1)/mN * slashed(q) + (F(2)/mN**2 * SP(pf,q)  &
                      & + F(3)/mN**2 * SP(pi,q) + F(4))* unit4
            F1_mN=F(1)/mN
            F2_mN2=F(2)/mN**2
            F3_mN2=F(3)/mN**2
            if (process.ne.EM) then
               F5_mN=F(5)/mN
               F6_mN2=F(6)/mN**2
               F8_mN2=F(8)/mN**2
               matrix_2= -F(5)/mN * MatMul(gamma5,slashed(q))  &
                       &  - F(6)/mN**2* SP(pf,q) * gamma5
            end if
         end if
         do_Once=.false.
      end if

      if (speedup) then
         matrix=  metricTensor(lambda,nu)* matrix_1 &
              & - q(lambda)*  &
              & ( F1_mN * gamma(:,:,nu)  &
              & + (F2_mN2 * pf(nu)  + F3_mN2 * pi(nu) )* unit4 )
      else
         write(*,*) 'only used with speedup'
         stop
      end if

     if (process.ne.EM) then
           matrix=matrix &
                & + metricTensor(lambda,nu)* matrix_2  &
                & + q(lambda)*( F5_mN * MatMul(gamma5,gamma(:,:,nu)) &
                & + F6_mN2 * pf(nu)  * gamma5    )   &
                & - (metricTensor(lambda,nu)* F(7)+q(lambda)*q(nu)*F8_mN2)*gamma5
     end if
   end function j_dagger

  end function hadronTensor_3_2












  !****************************************************************************
  !****is* hadronTensor_ResProd/parityError
  ! PURPOSE
  ! Error Message
  !****************************************************************************
  subroutine parityError(parity,who)
    integer :: parity
    character(*) :: who

    write(*,*) 'Strange parity in ', who
    write(*,*) 'Parity=', parity
    write(*,*) 'Stop'
    stop
  end subroutine parityError


end module hadronTensor_ResProd
