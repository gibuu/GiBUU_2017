!******************************************************************************
!****m* /vertices
! NAME
! module vertices
! PURPOSE
! Provides vertex functions
!******************************************************************************
module vertices
  implicit none
  private

  !****************************************************************************
  !****g* vertices/debug
  ! SOURCE
  logical, parameter      :: debug=.false.
  ! PURPOSE
  ! switch on/off debugging output
  !****************************************************************************

  logical,save      :: initFlag=.true.


  public ::  vertex_gammaN_to_spin3_2


contains


  subroutine readInput
!!$    use output, only: Write_ReadingInput
!!$
!!$    integer :: ios
!!$
!!$    !***************************************************************************
!!$    !****n* vertices/vertices_nl
!!$    ! NAME
!!$    ! NAMELIST /vertices_nl/
!!$    ! PURPOSE
!!$    ! This Namelist includes:
!!$    ! * debug
!!$    !***************************************************************************
!!$    NAMELIST /vertices_nl/ debug
!!$
!!$    call Write_ReadingInput('vertices_nl',0)
!!$
!!$    rewind(5)
!!$    read(5,nml=vertices_nl,IOSTAT=ios)
!!$    call Write_ReadingInput('vertices_nl',0,ios)
!!$
!!$    write(*,'(A20,L4)')    '  Debug=            ',Debug
!!$
!!$    call Write_ReadingInput('vertices_nl',1)

  end subroutine readInput




  !****************************************************************************
  !****f* vertices/vertex_gammaN_to_spin3_2_negParity
  ! NAME
  ! function vertex_gammaN_to_spin3_2(q,p,resID,resCharge) result(vertex)
  !
  ! PURPOSE
  ! * Provides the vertex for gamma N -> spin 3/2 resonance
  !
  ! INPUTS
  ! * real, dimension(0:3), intent(in) :: q,p        -- photon and nucleon momenta
  ! * integer             , intent(in) :: resID      -- ID of the resonance
  ! * integer             , intent(in) :: resCharge  -- Charge of the resonance
  !
  ! OUTPUT
  ! * complex,dimension(0:3,0:3,0:3,0:3)  :: vertex
  ! * First two indices =mu,nu,  i.e.   (Gamma^{mu,nu} )_ij=vertex(mu,nu,i,j)
  !****************************************************************************
  function vertex_gammaN_to_spin3_2(q,p,resID,resCharge) result(vertex)
    use leptonicID, only: EM
    use formFactor_ResProd, only: getFormfactor_Res
    use minkowski, only: SP, gamma, slashed, metricTensor,abs4
    use matrix_module, only: unit4
    use particleProperties, only: hadron
    use constants, only: mN

    real, dimension(0:3), intent(in)    :: q,p
    integer             , intent(in)    :: resID
    integer             , intent(in)    :: resCharge

    complex,dimension(0:3,0:3,0:3,0:3)  :: vertex

    logical              :: success
    real, dimension(1:8) :: formfactor
    real                 :: Qsquared
    real, dimension(3:6) :: C
    real, dimension(0:3) :: pi,pf
    integer              :: resParity, lambda, nu

    complex, dimension(0:3,0:3) :: term,helper


    if (initflag) then
       call readInput()
       initFlag=.false.
    end if

    QSquared=-SP(q,q)
    resParity=(-1)**(1+hadron(resID)%AngularMomentum)
    if (debug) write(*,'(A,I4,A,I4,A,G15.4)') 'Resonance ID: ',resID,'=> parity=',resParity,'QSquared=',qsquared

    ! Retrieve form factors
    pf=p+q
    pi=p

    formfactor=getFormfactor_Res(QSquared,abs4(pf),resID,resCharge,EM,success)
    if (.not.success) then
       vertex=0.
       return
    end if


    c(3:6)=formFactor(1:4)
    if (debug) write(*,'(A,4G15.4)') 'C=',C




    term=C(3)/mN * slashed(q) + (C(4)/mN**2 * SP(pf,q) + C(5)/mN**2 * SP(pi,q) + C(6))* unit4
    do lambda=0,3
       do nu=0,3
          helper=  metricTensor(lambda,nu)*  term &
               & - q(lambda)* ( C(3)/mN * gamma(:,:,nu)  + (C(4)/mN**2 * pf(nu) + C(5)/mN**2 * pi(nu)         )* unit4 )
          select case (resParity)
          case (-1)
             vertex(lambda,nu,:,:)  =  helper
          case (1)
             ! Multiply by gamma(5)
             vertex(lambda,nu,:,:)  =  MatMul(helper,gamma(:,:,5))
          case default
             stop 'strange parity in vertex_gammaN_to_spin3_2'
          end select
       end do
    end do

  end function vertex_gammaN_to_spin3_2


end module vertices
