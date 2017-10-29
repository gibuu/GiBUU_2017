
!***************************************************************************
!****m* /hadronicTensorQE
! NAME
! module hadronicTensorQE
!
! PURPOSE
! This module defines the hadronic tensors for QE EM- and weak scattering
! in particular to test different descriptions of in-medium recipes.
! For actual runs, matrixElementQE is much faster (and includes RPA correlations)
!
!**************************************************************************
module hadronicTensorQE
  implicit none
  private

  public:: hadronicTensorforQE

  logical,save      :: debugflag=.false.
  logical,save      :: initFlag=.true.
  logical,save      :: useExtraterm=.true.
  integer,save      :: whichcurrent=1

contains

  subroutine initInput
    use output, only: Write_ReadingInput

    integer :: ios

    NAMELIST /hadronicTensorQE/ debugFlag,whichcurrent

    rewind(5)
    read(5,nml=hadronicTensorQE,IOSTAT=ios)
    call Write_ReadingInput("hadronicTensorQE",0,ios)

    write(*,*)    'debugging?', debugFlag
 
    call Write_ReadingInput('hadronicTensorQE',1)

  end subroutine initInput

  !********************************************************************************************
  !****f* hadronicTensorQE/hadronicTensorforQE
  ! NAME
  !  function hadronicTensorforQE(pi,pf,resID,targetCharge,process)  result(success)
  !
  ! PURPOSE
  ! This function returns the hadronic tensor for QE scattering
  ! by a current type defined by "process" (e.g. Gamma* N -> N').
  ! 
  ! INPUTS
  ! *  real, dimension(0:3),intent(in) :: pi            -- 4-momentum of incoming nucleon
  ! *  real, dimension(0:3),intent(in) :: pf            -- 4-momentum of outgoing nucleon
  ! *  integer,intent(in)              :: process       -- type of the process (EM, CC, NC)
  ! *  integer,intent(in)              :: targetCharge  -- Charge of target nucleon
  !
  ! OUTPUT
  ! * complex, dimension(0:3,0:3) :: matrix   -- H^mu nu
  ! * logical, intent(out)        :: success  -- true if resonance has coupling strength to the process;
  ! *                                            false if there is no coupling, therefore H^mu nu=0
  ! NOTES
  ! * We use the Peskin notation where the spin sum (u(p) ubar(p))=slashed(p)+m     !!
  !**********************************************************************************************
  function hadronicTensorforQE(pi,pf,targetCharge,process,matrix,whichcurrentIN)   result(success)
    use  FF_QE_nucleonScattering
    use minkowski, only: SP,abs4
    use leptonicID

    real, dimension(0:3),intent(in)             :: pi, pf
    integer,intent(in)                          :: process, targetCharge
    complex, dimension(0:3,0:3),intent(inout)   :: matrix
    logical                                     :: success 
    integer, intent(in),optional                :: whichcurrentIn

    real :: F1,F2   ! Vector Form factors
    real :: FA,FP   ! Axial Form factors
    real :: GE,GM  ! Sachs form factors

    success=.false.
    matrix=0.  

    if(initFlag) then
       call initInput
       initFlag=.false.
    end if

    if(present(whichcurrentIN)) then
       whichcurrent=whichcurrentIn
    end if

    If( debugflag ) then
       write(*,*) -SP(pf-pi,pf-pi)
       write(*,*) sp(pf,pf)
    end if
    if(sp(pf,pf).le.0) then
       return
    end if

    if(targetcharge.eq.proton.and.process.eq.CC) return
    if(targetcharge.eq.neutron.and.process.eq.antiCC) return
    
    call formfactors_QE(-SP(pf-pi,pf-pi),process,targetcharge,F1,F2,FA,FP,GE,GM) 

    matrix=hadronicTensor(pi,pf,F1,F2,FA,FP,GE,GM,process)

    success=.true.

  end function hadronicTensorforQE

  !***************************************************************************
  !****f* hadronicTensorQE/hadronicTensor
  ! NAME
  ! function hadronicTensor(pi,pf,G,GA,process) result(matrix)
  !
  ! PURPOSE
  ! This function returns the hadronic tensor for nucleons. 
  ! 
  ! INPUTS
  ! *  real, dimension(1:2),intent(in) :: F             -- Vector or EM-Form factors
  ! *  real, dimension(1:2),intent(in) :: FA            -- Axial form factors
  ! *  real, dimension(0:3),intent(in) :: pi            -- 4-momentum of incoming nucleon
  ! *  real, dimension(0:3),intent(in) :: pf            -- 4-momentum of outgoing resonance
  ! *  integer,intent(in)              :: process       -- type of the process (EM, CC, NC)
  !
  ! OUTPUT
  ! * complex, dimension(0:3,0:3) :: matrix -- H^mu nu
  !**************************************************************************
  function hadronicTensor(pi,pf,F1,F2,FA,FP,GE,GM,process) result(matrix)
    use leptonicID, only : EM
    use minkowski, only          : abs4, gamma0, slashed
    use matrix_module, only      : MatrixMult, dagger, trace, unit4

    integer, intent(in)            :: process
    real, dimension(0:3),intent(in):: pi,pf
    real,intent(in)                :: F1,F2   ! Vector or EM-Form factors
    real,intent(in)                :: FA,FP  ! Axial form factors  
    real,intent(in)                :: GE,GM  ! Sachs form factors
    complex, dimension(0:3,0:3)    :: matrix

    complex, dimension(0:3,0:3) :: projector_in,projector_out,j_mu,j_nu
    real                        :: mi,mf
    integer                     :: mu,nu
    real, dimension(0:3),save        :: q

    mi=abs4(pi)
    mf=abs4(pf)

    q=pf-pi

    projector_out=slashed(pf)+mf*unit4
    projector_in=slashed(pi)+mi*unit4

    matrix=0.
    do mu=0,3
       do nu=0,3
          j_mu=j(mu)
          j_nu=j(nu)

          !current conserving a la deForest:
          !note that we assume here, that q is in z direction
          if(whichcurrent.eq.3.or.whichcurrent.eq.5) then
             if(mu.eq.3) j_mu=q(0)*j(0)/sqrt(Dot_product(q(1:3),q(1:3)))
             if(nu.eq.3) j_nu=q(0)*j(0)/sqrt(Dot_product(q(1:3),q(1:3)))
          end if

          matrix(mu,nu)=1./2.*Trace(MatrixMult(projector_out,j_mu,projector_in,gamma0,dagger(j_nu),gamma0))
       end do
    end do

  contains

    function j(mu) result(matrix)
      use minkowski, only          : gamma, gamma5, SP, slashed, sigma4, metricTensor,ii
      use matrix_module, only      : unit4
      use constants, only          : mN

      integer, intent(in) :: mu
      complex, dimension(0:3,0:3) :: matrix,sigma4_q
      real, dimension(0:3),save        :: k
      integer :: rho,alpha
      real ,save   :: mass_mu, QSquared,qM
      
            
      k=pf+pi
      mass_mu=2.*mN
      QSquared=-SP(q,q)

      sigma4_q=0.
      do rho=0,3
         do alpha=0,3
            if(alpha.ne.rho) cycle
            sigma4_q=sigma4_q+sigma4(mu,rho)*q(alpha)*metricTensor(rho,alpha)
         end do
      end do

      
      select case (whichcurrent)

      case (1) !our standard approach with current-conserving term 
         matrix= F1 *(gamma(:,:,mu)+q(mu)*slashed(q)/QSquared) + F2/mass_mu * ii* sigma4_q
 

      case (2,3) !2: our standard approach WITHOUT current-conserving term
                 !OR
                 !3: deForest current no. 2 with current-conserving via longitudinal component (see Tinas notes)
         matrix= F1 *gamma(:,:,mu) + F2/mass_mu * ii* sigma4_q


      case (4,5) !4: different current (see deForest, his current no. 1), no current-conserving
                 !OR
                 !5: deForest current no. 1 with current-conserving via longitudinal component (see Tinas notes)
         matrix= (F1+F2)*gamma(:,:,mu) - F2/mass_mu * k(mu)*unit4
         
      case (6) !Saclay ansatz 
         qM=Qsquared/(mi+mf)**2
         matrix= GM*gamma(:,:,mu) + k(mu)/(mi+mf)*(GE-GM)/(1+qM)*unit4+(mi-mf)/QSquared*q(mu)*(GE+qM*GM)/(1+qM)*unit4
         
      case default
         write(*,*) 'problem with whichcurrent in hadronTensorQE -> STOP'
         stop
      end select


      !axial part is the same for all
      if (process.ne.EM) then
         matrix= matrix + MatMul(FA * gamma(:,:,mu) + FP/mN * q(mu)*unit4,gamma5)
      end if
      
    end function j


  end function hadronicTensor



end module hadronicTensorQE
