!******************************************************************************
!****m* /NievesHadronTensor
! NAME
! module NievesHadronTensor
! PURPOSE
! This module implements the hadronic tensor for 1-pion production.
!******************************************************************************
module NievesHadronTensor

  implicit none
  private

  !****************************************************************************
  !****g* NievesHadronTensor/DeltaPole
  ! SOURCE
  logical, save :: DeltaPole=.true.
  ! PURPOSE
  !****************************************************************************

  !****************************************************************************
  !****g* NievesHadronTensor/crossedDelta
  ! SOURCE
  logical, save :: crossedDelta=.true.
  ! PURPOSE
  !****************************************************************************

  !****************************************************************************
  !****g* NievesHadronTensor/nucleonPole
  ! SOURCE
  logical, save :: nucleonPole=.true.
  ! PURPOSE
  !****************************************************************************

  !****************************************************************************
  !****g* NievesHadronTensor/crossedNucleonPole
  ! SOURCE
  logical, save :: crossedNucleonPole=.true.
  ! PURPOSE
  !****************************************************************************

  !****************************************************************************
  !****g* NievesHadronTensor/contactTerm
  ! SOURCE
  logical, save :: contactTerm=.true.
  ! PURPOSE
  !****************************************************************************

  !****************************************************************************
  !****g* NievesHadronTensor/pionPole
  ! SOURCE
  logical, save :: pionPole=.true.
  ! PURPOSE
  !****************************************************************************

  !****************************************************************************
  !****g* NievesHadronTensor/pionInFlight
  ! SOURCE
  logical, save :: pionInFlight=.true.
  ! PURPOSE
  !****************************************************************************

  !****************************************************************************
  !****g* NievesHadronTensor/debug_HNV_hadrTens
  ! SOURCE
  integer, parameter :: debug_HNV_hadrTens=1
  ! PURPOSE
  !****************************************************************************

  logical, save :: initflag=.true.

  ! matrix of Clebsch-Gordon coefficients for 1-pi production for 7 diagrams
  real, dimension(1:2,1:7,-1:1,0:1), save  :: CG    ! 1:2 process_ID,  1:7 diagram, -1:1 pion charge, 0:1 nucleon charge,



  public :: NievesHadronTensor_1pi


contains

  subroutine initInput_HNV_hadrTens()
    use output

    integer :: ios

    !**************************************************************************
    !****n* NievesHadronTensor/nl_NievesHadronTensor
    ! NAME
    ! NAMELIST nl_NievesHadronTensor
    ! PURPOSE
    ! Includes the switches:
    ! * DeltaPole
    ! * crossedDelta
    ! * nucleonPole
    ! * crossedNucleonPole
    ! * contactTerm
    ! * pionPole
    ! * pionInFlight
    !**************************************************************************
    NAMELIST /nl_NievesHadronTensor/ DeltaPole, crossedDelta, nucleonPole, crossedNucleonPole, &
         contactTerm, pionPole, pionInFlight

    call Write_ReadingInput('NievesHadronTensor',0)
    rewind(5)
    read(5,nml=nl_NievesHadronTensor,IOSTAT=ios)
    call Write_ReadingInput("NievesHadronTensor",0,ios)

    call Write_ReadingInput('NievesHadronTensor',1)


    ! this is matrix of Clebsch-Gordon coefficient for 1-pi productio in EM and CC reactions
    ! indices: process_ID (1=em, 2-cc), digram_number, pion charge(-1:1), nucleon_charge(0:1)
    ! process_ID  :  1=EM  2=CC
    ! diagrams    :  1=Dp  2=cDp  3=Np  4=cNp  5=CT  6=pp  7=pF

    !                       pi-n        pi0n          pi+n           pi-p           pi0p      pi+p

    CG(1,1,:,:)= reshape( (/ 0.      , sqrt(2./3.), -sqrt(1./3.), sqrt(1./3.),  sqrt(2./3.), 0.       /) , (/3,2/) )
    CG(2,1,:,:)= reshape( (/ sqrt(3.), sqrt(2./3.),  sqrt(1./3.), sqrt(1./3.), -sqrt(2./3.), sqrt(3.) /) , (/3,2/) )

    CG(1,2,:,:)= reshape( (/ 0.         ,  sqrt(2./3.),  sqrt(1./3.), -sqrt(1./3.),  sqrt(2./3.), 0.          /) , (/3,2/) )
    CG(2,2,:,:)= reshape( (/ sqrt(1./3.), -sqrt(2./3.),  sqrt(3.),     sqrt(3.),     sqrt(2./3.), sqrt(1./3.) /) , (/3,2/) )

    CG(1,3,:,:)= reshape( (/ 0.,  -1.,  sqrt(2.),  sqrt(2.),  1.,  0. /) , (/3,2/) )
    CG(2,3,:,:)= reshape( (/ 0.,  -1.,  sqrt(2.),  sqrt(2.),  1.,  0. /) , (/3,2/) )

    CG(1,4,:,:)= reshape( (/ 0.,       -1.,  sqrt(2.),  sqrt(2.),  1.,  0.       /) , (/3,2/) )
    CG(2,4,:,:)= reshape( (/ sqrt(2.),  1.,  0.,        0.,       -1.,  sqrt(2.) /) , (/3,2/) )

    CG(1,5,:,:)= reshape( (/ 0.,  0.,       -1.,  1.,  0.,       0. /) , (/3,2/) )
    CG(2,5,:,:)= reshape( (/ 1.,  sqrt(2.), -1., -1., -sqrt(2.), 1. /) , (/3,2/) )

    CG(1,6,:,:)= reshape( (/ 0.,  0.,       0.,  0.,  0.,       0. /) , (/3,2/) )
    CG(2,6,:,:)= reshape( (/ 1.,  sqrt(2.), -1., -1., -sqrt(2.), 1. /) , (/3,2/) )

    CG(1,7,:,:)= reshape( (/ 0.,  0.,       -1.,  1.,  0.,       0. /) , (/3,2/) )
    CG(2,7,:,:)= reshape( (/ 1.,  sqrt(2.), -1., -1., -sqrt(2.), 1. /) , (/3,2/) )

    initFlag = .false.

  end subroutine initInput_HNV_hadrTens


  !****************************************************************************
  !****f* NievesHadronTensor/NievesHadronTensor_1pi
  ! NAME
  !  function   NievesHadronTensor_1pi(process_ID,q,p_in, position, charge_in, p_out,
  !                                                charge_out, pion_charge_out, matrix) result(success)
  !
  ! PURPOSE
  ! This function returns the hadronic tensor for 1-pion production (that is:  lepton A -->  lepton  nucleon pion  A-1)
  ! by a current type defined by "process" (e.g. CC, NC, em).
  !
  ! INPUTS
  ! *  integer,intent(in)               :: process       -- type of the process (EM, CC, NC)
  ! *  real, dimension(0:3), intent(in) :: q             -- momentum transferred
  ! *  real, dimension(0:3), intent(in) :: p_in           -- 4-momentum of incoming nucleon
  ! *  real, dimension(1:3), intent(in) :: position       -- position --- is used for in-medium propagator
  ! *  real,intent(in) :: charge_in      -- charge of incoming nucleon
  ! *  real, dimension(0:3),intent(in)  :: p_out          -- 4-momentum of outgoing nucleon
  ! *  real,intent(in) :: charge_out     -- charge of outgoing nucleon
  ! *  real,intent(in) :: pion_charge_out     -- charge of outgoing pion
  !
  !
  ! OUTPUT
  ! * complex, dimension(0:3,0:3) :: matrix -- H^mu nu
  ! * logical, intent(out)        :: success  -- true if there is a coupling strength to the process;
  ! *                                            false if there is no coupling, therefore H^mu nu=0
  ! NOTES
  ! * We use the Peskin notation where the spin sum(u(p) ubar(p))=slashed(p)+m     !!
  !****************************************************************************
  function   NievesHadronTensor_1pi(process_ID,q,p_in, position, charge_in, &
       p_out,charge_out,ppi_out, pion_charge_out, matrix) result(success)

    use minkowski
    use matrix_module, only: unit4, MatrixMult, Trace
    use formFactor_ResProd, only: getFormfactor_Res
    use FF_QE_nucleonScattering, only: formfactors_QE
    use spinProjector
    use idtable, only: Delta, rho
    use particleProperties, only: hadron
    use leptonicID, only: CC, EM, antiCC

    integer, intent(in) :: process_ID
    real, dimension(0:3), intent(in) :: q, p_in, p_out, ppi_out
    real, dimension(1:3), intent(in) :: position !--- is used for in-medium propagator
    integer, intent(in) :: charge_in, charge_out, pion_charge_out
    complex, dimension(0:3,0:3), intent(out) :: matrix
    logical :: success


    real, dimension(1:7,-1:1,0:1)   :: ClebGor
    real :: m_in, m_out
    real, save :: Q2
    real, dimension(1:8),save :: ff_Delta
    real, dimension(1:4),save :: ff_nucleon

    real :: FV_CT, F_rho, F_PF
    integer :: mu, nu
    complex, dimension(0:3,1:4,1:4) :: j

    logical ::  success_ff

    complex, dimension(1:4,1:4) :: projector_in,projector_out

    if (initFlag) call initInput_HNV_hadrTens()

    success = .false.
    ClebGor=CG(abs(process_ID),:,:,:)

    j=0
    Q2=-SP(q,q)

    if (debug_HNV_hadrTens.ge.1) write(*,'(4(A,I5))') '# process_ID=', process_ID, &
         '    charge_in=', charge_in, '    charge_out=', charge_out, '    pion_charge_out=',pion_charge_out


    ! currents for direct  and crossed   Delta production
    if (.not.DeltaPole .and. .not.crossedDelta) then
    else
       ! direct Delta
       if (.not.DeltaPole) then
       else

          ff_Delta=getFormfactor_Res(Q2,hadron(Delta)%mass,Delta,charge_in,process_ID,success_ff)  ! hadron(Delta)%mass is bare mass, OK
          if (.not.success_ff) then
             write(*,*) 'in NievesHadronTensor_1pi: formfactors for Delta are not available for Q2=', Q2
             return
          end if
          if (debug_HNV_hadrTens.ge.1) write(*,'(A,g10.3,A,8g12.5)') '# Dp: Delta form factors at Q2=', Q2, '   are :', ff_Delta

          j = j + j_DeltaPole(charge_out,pion_charge_out,p_in,q,ppi_out)

       end if

       ! crossed Delta
       if (.not.crossedDelta) then
          !if (debug_HNV_hadrTens.ge.4) write(*,*) 'no crossed-Delta diagram is included'
       else
          ! for the crossed Delta pole the WNN vertex is related to the outgoing and virtual nucleons
          ! important for EM interactions
          ff_Delta=getFormfactor_Res(Q2,hadron(Delta)%mass,Delta,charge_out,process_ID,success_ff)
          ! 2009-07-15 for debbuging purposes only !!!
          if (.not.success_ff) then
             write(*,*) 'in NievesHadronTensor_1pi: formfactors for Delta are not available for Q2=', Q2
             return
          end if
          if (debug_HNV_hadrTens.ge.1) write(*,'(A,g10.3,A,8g12.5)') '# cDp: Delta form factors at Q2=', Q2, '   are :', ff_Delta

          j = j + j_crossedDelta(charge_out,pion_charge_out,p_out,(-q),ppi_out)
       end if
    end if





    ! currents for direct  and crossed   nucleon poles
    if (.not.nucleonPole .and. .not.crossedNucleonPole) then
    else

       ! nucleon pole
       if (.not.nucleonPole) then
       else
          call formfactors_QE(Q2,process_ID,charge_in,ff_nucleon(1),ff_nucleon(2),ff_nucleon(3),ff_nucleon(4))
          if (debug_HNV_hadrTens.ge.1) write(*,'(A,g10.3,A,4g12.5)') '# Np: Nucleon form factors at Q2=', Q2, '   are :', ff_nucleon
          j = j + j_directNucleonPole(charge_out,pion_charge_out)
       end if

       ! crossed nucleon pole
       if (.not.crossedNucleonPole) then
       else
          ! for the crosse nucleon pole the WNN vertex is related to the OUTgoing and virtual nucleons
          ! important for EM interactions
          call formfactors_QE(Q2,process_ID,charge_out,ff_nucleon(1),ff_nucleon(2),ff_nucleon(3),ff_nucleon(4))
          if (debug_HNV_hadrTens.ge.1) write(*,'(A,g10.3,A,4g12.5)') '# cNp: Nucleon form factors at Q2=', Q2,'   are :',ff_nucleon
          j = j + j_crossedNucleonPole(charge_out,pion_charge_out)
       end if

    end if


    ! contact term
    if (.not.contactTerm) then
    else
       call formfactors_QE(Q2,CC,charge_in,ff_nucleon(1),ff_nucleon(2),ff_nucleon(3),ff_nucleon(4))
       FV_CT=ff_nucleon(1)
       if (debug_HNV_hadrTens.ge.2) write(*,'(3(A,4g12.5))') 'Contact term:  abs4Sq(q-ppi)=',abs4Sq(q-ppi_out),  '    F_rho=',F_rho

       select case (process_ID)
       case (CC,antiCC)
          if (debug_HNV_hadrTens.ge.1) write(*,'(3(A,4g12.5))') 'Contact term:  Q2=',Q2, '    FV_CT=',FV_CT
          F_rho=1./( 1. - abs4Sq(q-ppi_out)/hadron(rho)%mass**2 )
       case (EM)
          F_rho=0.    ! because axial part is zero
       case default
          write(*,*) 'STOP.  Wrong process ID=', process_ID
          stop
       end select
       j = j + j_contactTerm(charge_out,pion_charge_out)
    end if


    ! pion pole
    if (.not.pionPole) then
    else
       select case (process_ID)
       case (CC,antiCC)
          F_rho=1./( 1. - abs4Sq(q-ppi_out)/hadron(rho)%mass**2 )
          if (debug_HNV_hadrTens.ge.1) write(*,'(3(A,4g12.5))') 'Pion pole:  abs4Sq(q-ppi)=',abs4Sq(q-ppi_out),  '    F_rho=',F_rho
          j = j + j_pionPole(charge_out,pion_charge_out)
       case (EM)
          ! j_pionPole=0 ! only weak boson can be converted to pion
       case default
          write(*,*) 'STOP.  Wrong process ID=', process_ID
          stop
       end select
    end if



    ! pion in flight
    if (.not.pionInFlight) then
    else
       select case (process_ID)
       case (CC,antiCC,EM)
          call formfactors_QE(Q2,CC,charge_in,ff_nucleon(1),ff_nucleon(2),ff_nucleon(3),ff_nucleon(4))
          F_PF=ff_nucleon(1)
          if (debug_HNV_hadrTens.ge.1) write(*,'(3(A,4g12.5))') 'Pion in flight:  Q2=',Q2, '    F_PF=',F_PF
          j = j + j_pionInFlight(charge_out,pion_charge_out)
       case default
          write(*,*) 'STOP.  Wrong process ID=', process_ID
          stop
       end select

    end if




    ! --------------------------------------------------------
    ! The hadronic tensor itself
    !
    m_in=abs4(p_in)
    m_out=abs4(p_out)
    if (debug_HNV_hadrTens.ge.2) write(*,'(2(A,g11.4))') 'In  NievesHadronTensor_1pi:   m_in=', m_in, '     m_out=', m_out
    projector_in=slashed(p_in)+m_in*unit4
    projector_out=slashed(p_out)+m_out*unit4

    do mu=0,3
       do nu=0,3
          matrix(mu,nu)=Trace(  MatrixMult(projector_out,(j(mu,:,:)),projector_in,tilde(j(nu,:,:)))  )/2.
       end do
    end do
    success=.true.
    ! --------------------------------------------------------------

  contains



    ! direct Delta
    function j_DeltaPole(charge_out,pion_charge_out,p_in,q,ppi_out) result (matrix)

      use spectralFunc, only: propagator_nenner
      use minkowski, only: SP, metricTensor
      use leptonicID

      integer, intent(in) :: charge_out, pion_charge_out
      real, dimension(0:3), intent(in) :: p_in, q, ppi_out

      complex, dimension(0:3,1:4,1:4) :: matrix, j, temp
      integer :: al, be
      real :: coeff_final
      complex :: propag
      real, dimension(0:3) :: W

      matrix=0

      coeff_final=ClebGor(1,pion_charge_out,charge_out)
      if (coeff_final.eq.0) return

      W=p_in+q
      j=0
      do mu=0,3
         do al=0,3
            temp=0
            do be=0,3
               ! which way is faster?
               !j(mu,:,:)=j(mu,:,:) + vertex_DeltaNpi(al,ppi_out) * MatrixMult(spin32proj(2,al,be,W), vertex_WNDelta(be,mu,p_in,q)  ) &
               !& *metricTensor(be,be)*metricTensor(al,al)
               temp(al,:,:) = temp(al,:,:) &
                    + metricTensor(be,be)*MatrixMult( spin32proj(2,al,be,W,.false.), vertex_WNDelta(be,mu,p_in,q) )
            end do
            j(mu,:,:)=j(mu,:,:) + metricTensor(al,al)*vertex_DeltaNpi(al,ppi_out)*temp(al,:,:)
         end do
      end do

      propag=propagator_nenner(2,(charge_out+pion_charge_out),W,position)
      matrix=j*coeff_final/propag

      if (debug_HNV_hadrTens.ge.2) write(*,'(6(A,g12.5))') 'In NievesHadronTensor_1pi, j_DeltaPole:    sqrt(W2)=', &
           & sqrt(SP(W,W)), '  coeff_final=', coeff_final, &
           & '  propagator_nenner=', REAL(propag), ' + i ', AIMAG(propag)
    end function j_DeltaPole



    ! crossed Delta
    function j_crossedDelta(charge_out,pion_charge_out,momen_in,q_in,ppi_out) result (matrix)

      use spectralFunc, only: propagator_nenner
      use minkowski, only: metricTensor !, abs4sq
      use leptonicID

      integer, intent(in) :: charge_out, pion_charge_out
      real, dimension(0:3), intent(in) :: momen_in, q_in, ppi_out
      complex, dimension(0:3,1:4,1:4) :: matrix, j, temp

      integer :: al, be
      real :: coeff_final
      complex :: propag
      real, dimension(0:3) :: pDelta

      matrix=0

      coeff_final=ClebGor(2,pion_charge_out,charge_out)
      if (coeff_final.eq.0) return

      !if (debug_HNV_hadrTens.ge.2) write(*,*) 'Start j_crossedDelta'
      ! 2010-07-30 this check was good for free nucleon, BUT is not necessarily valid for reactions on nuclei
      !if( (abs4Sq(p_out-q)) > ((mN-mPi)**2) ) &
      !    & write(*,'(4(A,g12.5))') 'Error:  In j_crossedDelta  (p_out-q)^2=', (abs4Sq(pDelta)), &
      !    & '   should be equal to (p-ppi)^2=', (abs4Sq(p_in-ppi_out)), &
      !    & '   and should be less than (mN-mpi)^2=', ((mN-mPi)**2)

      pDelta=momen_in+q_in
      j=0
      do mu=0,3
         do be=0,3
            temp=0
            do al=0,3
               ! .false. in Delta--ProjectionOperator  means that M_Delta and not sqrt(p_mu*p_mu) is used
               !  this is  the only possibility in u-channel

               ! j(mu,:,:)=j(mu,:,:) + &
               ! & gg(al)*MatMul( vertex_WNDelta_DirConj(mu,al,p_out,(-q)), spin32proj(2,al,be,(p_out-q),.false.) ) &
               ! & * vertex_DeltaNpi(be)*gg(be)
               temp(be,:,:)=metricTensor(al,al) &
                    & * MatMul( vertex_WNDelta_DirConj(mu,al,momen_in,q_in), spin32proj(2,al,be,pDelta,.false.) )
            end do
            j(mu,:,:)=j(mu,:,:) + metricTensor(be,be)*temp(be,:,:)*vertex_DeltaNpi(be,ppi_out)
         end do
      end do

      propag=propagator_nenner(2,(charge_in-pion_charge_out),pDelta,position)
      matrix=j*coeff_final/propag

      if (debug_HNV_hadrTens.ge.2) write(*,'(6(A,g12.5))') 'In NievesHadronTensor_1pi, j_crossedDelta:    (p_out-q)^2=', &
           & SP(pDelta,pDelta), '  coeff_final=', coeff_final, &
           & '  propagator_nenner=', REAL(propag), ' + i ', AIMAG(propag)

    end function j_crossedDelta






    ! nucleon pole
    function j_directNucleonPole(charge_out,pion_charge_out) result (matrix)

      use spectralFunc, only: propagator_nenner

      use minkowski, only: slashed !, abs4Sq
      use matrix_module, only: MatrixMult
      use leptonicID

      integer, intent(in) :: charge_out, pion_charge_out

      complex, dimension(0:3,1:4,1:4) :: matrix,j

      real :: coeff_final
      complex, dimension(1:4,1:4) :: spin12proj
      real, dimension(0:3) :: W
      complex :: propag

      matrix=0

      coeff_final=ClebGor(3,pion_charge_out,charge_out)
      if (coeff_final.eq.0) return

      W=p_in+q
      spin12proj=slashed(W)+abs4(W)*unit4 ! projector of the virtual nucleon
      !spin12proj=slashed(W)+mN*unit4

      j=0
      do mu=0,3
         j(mu,:,:)=j(mu,:,:) + MatrixMult(  vertex_NNpi(ppi_out), spin12proj, vertex_WNN(mu,p_in,q)  )
      end do

      propag=propagator_nenner(1,(charge_out+pion_charge_out),W,position)
      matrix=-j*coeff_final/propag ! "-" according to HNV

      if (debug_HNV_hadrTens.ge.2) write(*,'(6(A,g12.5))') 'In NievesHadronTensor_1pi, j_directNucleonPole:    sqrt(W2)=', &
           & sqrt(SP(W,W)), '  coeff_final=', coeff_final, &
           & '  propagator_nenner=', REAL(propag), ' + i ', AIMAG(propag)

    end function j_directNucleonPole





    ! crossed nucleon pole
    function j_crossedNucleonPole(charge_out,pion_charge_out) result (matrix)

      use spectralFunc, only: propagator_nenner

      use minkowski, only: slashed, abs4Sq
      use matrix_module, only: MatrixMult
      use constants, only: mN, mPi

      integer, intent(in) :: charge_out, pion_charge_out
      complex, dimension(0:3,1:4,1:4) :: matrix, j

      real :: coeff_final
      complex, dimension(1:4,1:4) :: spin12proj
      complex :: propag

      !if (debug_HNV_hadrTens.ge.2) write(*,*) 'Start j_crossedNucleonPole'


      matrix=0

      coeff_final=ClebGor(4,pion_charge_out,charge_out)
      if (coeff_final.eq.0) return


      if ( (abs4Sq(p_out-q)) > ((mN-mPi)**2) ) &
           & write(*,'(4(A,g12.5))') 'Error:  In j_crossedDelta  (p_out-q)^2=', (abs4Sq(p_out-q)), &
           & '   should be equal to (p-ppi)^2=', (abs4Sq(p_in-ppi_out)), &
           & '   and should be less than (mN-mpi)^2=', ((mN-mPi)**2)


      spin12proj=slashed(p_out-q)+mN*unit4 ! projector of the virtual nucleon

      j=0
      do mu=0,3
         j(mu,:,:)=j(mu,:,:) + MatrixMult(  vertex_WNN(mu,(p_out-q),q), spin12proj, vertex_NNpi(ppi_out)  )
      end do

      propag=propagator_nenner(1,(1-charge_out),(p_out-q),position)
      matrix=-j*coeff_final/propag  ! "-" according to HNV

      if (debug_HNV_hadrTens.ge.2) write(*,'(6(A,g12.5))') 'In NievesHadronTensor_1pi, j_crossedNucleonPole:    (p_out-q)^2=', &
           & SP((p_out-q),(p_out-q)), '  coeff_final=', coeff_final, &
           & '  propagator_nenner=', REAL(propag), ' + i ', AIMAG(propag)

    end function j_crossedNucleonPole




    ! contact term
    function j_contactTerm(charge_out,pion_charge_out) result (matrix)

      use constants, only: g_A, f_pi, ii
      use minkowski, only: gamma

      integer, intent(in) :: charge_out, pion_charge_out

      complex, dimension(0:3,1:4,1:4) :: matrix, j
      real :: coeff_final

      matrix=0

      coeff_final=ClebGor(5,pion_charge_out,charge_out)
      if (coeff_final.eq.0) return

      j=0
      do mu=0,3
         j(mu,:,:) = g_A*gamma(:,:,mu+8)*FV_CT - F_rho*gamma(:,:,mu)
      end do

      matrix= -ii*j*coeff_final/sqrt(2.)/f_pi  ! "-" according to HNV

      if (debug_HNV_hadrTens.ge.2) write(*,'(3(A,g12.5))') 'In NievesHadronTensor_1pi, j_contactTerm:    FV_CT=', FV_CT,  &
           &  '   F_rho=', F_rho, '    coeff_final=', coeff_final

    end function j_contactTerm




    ! pion pole
    function j_pionPole(charge_out,pion_charge_out) result (matrix)

      use constants, only: f_pi
      use minkowski, only: slashed, SP
      use constants, only: mPi, ii

      integer, intent(in) :: charge_out, pion_charge_out

      complex, dimension(0:3,1:4,1:4) :: matrix, j
      real :: coeff_final

      matrix=0

      coeff_final=ClebGor(6,pion_charge_out,charge_out)
      if (coeff_final.eq.0) return

      j=0
      do mu=0,3
         j(mu,:,:) = F_rho*q(mu)*slashed(q)/( SP(q,q) - mPi**2 )
      end do

      matrix=  -ii*j*coeff_final/sqrt(2.)/f_pi    ! "-" according to HNV

      if (debug_HNV_hadrTens.ge.2) write(*,'(3(A,g12.5))') 'In NievesHadronTensor_1pi, j_pionPole:  F_rho=', F_rho, &
           &  '    coeff_final=', coeff_final

    end function j_pionPole






    ! pion in flight
    function j_pionInFlight(charge_out,pion_charge_out) result (matrix)

      use constants, only: g_A, f_pi, mN, mPi, ii
      use minkowski, only: SP, gamma5

      integer, intent(in) :: charge_out, pion_charge_out

      complex, dimension(0:3,1:4,1:4) :: matrix, j
      real :: coeff_final, propag

      if ( (abs4Sq(p_out-q)) > ((mN-mPi)**2) ) &
           & write(*,'(4(A,g12.5))') 'Error:  In j_pionInFlight  (p_out-q)^2=', (abs4Sq(p_out-q)), &
           & '   should be equal to (p-ppi)^2=', (abs4Sq(p_in-ppi_out)), &
           & '   and should be less than (mN-mpi)^2=', ((mN-mPi)**2)

      matrix=0

      coeff_final=ClebGor(7,pion_charge_out,charge_out)
      if (coeff_final.eq.0) return

      j=0
      do mu=0,3
         j(mu,:,:) = F_PF*( 2.*ppi_out(mu) - q(mu) )*gamma5
      end do

      propag=SP(ppi_out-q,ppi_out-q) - mpi**2

      if (debug_HNV_hadrTens.ge.2) write(*,'(6(A,g12.5))') 'In NievesHadronTensor_1pi, j_pionInFlight:    (ppi_out-q)^2=', &
           & SP((ppi_out-q),(ppi_out-q)), '     F_PF=', F_PF, '     coeff_final=', coeff_final

      matrix=  -ii*j*2. *mN /propag *coeff_final*g_A/sqrt(2.)/f_pi  ! "-" according to HNV

    end function j_pionInFlight



    !   --------------------------------------------------------------------
    !  Delta -->   N    pion
    complex function vertex_DeltaNpi(a,ppi_out)

      use minkowski, only: abs4
      use constants, only: mPi, ii

      integer, intent(in) :: a
      real, dimension(0:3), intent(in) :: ppi_out

      real :: f_DeltaNpi


      f_DeltaNpi=2.14

      vertex_DeltaNpi=ii*ppi_out(a)*f_DeltaNpi/mPi

      if (debug_HNV_hadrTens.ge.4) write(*,'(2(A,g12.5),A,4g12.5)') 'In NievesHadtonTensor_1pi: vertex_DeltaNpi: f/mpi=', &
           & (f_DeltaNpi/mPi),  'ppi_out^2=', abs4(ppi_out),   '    ppi_out= ', ppi_out

    end function vertex_DeltaNpi
    !   -------------------------------------------------------------------






    !   -------------------------------------------------------------------
    !     W  N  -->   Delta
    function vertex_WNDelta(a,nu,momen_in,q_in) result(matrix)

      use minkowski, only: metricTensor, slashed, slashed5, gamma, gamma5, SP
      use matrix_module, only: unit4
      use leptonicID
      use constants, only: mN

      integer, intent(in) :: a,nu
      real, dimension(0:3), intent(in) :: momen_in,q_in ! momenta of the incoming nucleon and incoming boson
      complex, dimension(1:4,1:4) :: matrix
      real, dimension(0:3) :: W_in


      complex, dimension(1:4,1:4), save :: matrix1, matrix2
      real, save :: F1_mN, F2_mN2, F3_mN2, F5_mN, F6_mN2, F8_mN2

      matrix=0
      W_in=momen_in+q_in

      F1_mN=ff_Delta(1)/mN     ! c3V
      F2_mN2=ff_Delta(2)/mN**2 ! c4V
      F3_mN2=ff_Delta(3)/mN**2 ! c5V
      matrix1= F1_mN * slashed5(q_in) + (F2_mN2 * SP(W_in,q_in) + F3_mN2 * SP(momen_in,q_in))* gamma5

      if (process_ID.ne.EM) then
         F5_mN=ff_Delta(5)/mN     !c3A
         F6_mN2=ff_Delta(6)/mN**2 !c4A
         F8_mN2=ff_Delta(8)/mN**2 !c6A

         matrix2= F5_mN * slashed(q_in) + (  F6_mN2* SP(W_in,q_in) + ff_Delta(7) )  * unit4
      end if


      ! vector part
      matrix=matrix + metricTensor(a,nu)* matrix1 &
           & - q_in(a)*  ( F1_mN * gamma(:,:,(nu+8))  + (F2_mN2 * W_in(nu)  + F3_mN2 * momen_in(nu) )* gamma5  )


      ! axial part
      if (process_ID.ne.EM) then
         matrix=matrix + metricTensor(a,nu)* matrix2  &
              &   - q_in(a)*  ( F5_mN * gamma(:,:,nu) + F6_mN2 * W_in(nu)  * unit4    )    &
              & + q_in(a)*q_in(nu)*F8_mN2 * unit4

      end if
    end function vertex_WNDelta
    !   -------------------------------------------------------------------



    !   -------------------------------------------------------------------
    !     W  N  -->   Delta conjugated
    function vertex_WNDelta_DirConj(a,b,momen_in,q_in) result(matrix) ! gamma_0*dagger(vertex_WNDelta)*gamma_0

      use minkowski, only: metricTensor, slashed, slashed5, gamma, gamma5, SP
      use matrix_module, only: unit4
      use leptonicID
      use constants, only: mN

      integer, intent(in) :: a,b
      real, dimension(0:3), intent(in) :: momen_in,q_in
      complex, dimension(1:4,1:4) :: matrix

      complex, dimension(1:4,1:4), save :: matrix1, matrix2
      real, save :: F1_mN, F2_mN2, F3_mN2, F5_mN, F6_mN2, F8_mN2
      real, dimension(0:3) :: W_in

      ! minus in front of the terms with gamma5 (as compared to vertex_WNDelta)
      matrix=0
      W_in=momen_in+q_in

      F1_mN=ff_Delta(1)/mN     ! c3V
      F2_mN2=ff_Delta(2)/mN**2 ! c4V
      F3_mN2=ff_Delta(3)/mN**2 ! c5V
      matrix1= F1_mN * slashed5(q_in)   - (F2_mN2 * SP(W_in,q_in) + F3_mN2 * SP(momen_in,q_in))* gamma5 ! minus here

      if (process_ID.ne.EM) then
         F5_mN=ff_Delta(5)/mN     !c3A
         F6_mN2=ff_Delta(6)/mN**2 !c4A
         F8_mN2=ff_Delta(8)/mN**2 !c6A

         matrix2= F5_mN * slashed(q_in) + (  F6_mN2* SP(W_in,q_in) + ff_Delta(7) )  * unit4
      end if

      if (debug_HNV_hadrTens.ge.3) write(*,'(7(A,g12.5))') 'In  vertex_WNDelta_DirConj  F1_mN=', F1_mN, &
           & '    F2_mN2=', F2_mN2,  '    F3_mN2=', F3_mN2,  &
           & '    F5_mN=', F5_mN,  '    F6_mN2=', F6_mN2, '    F8_mN2=', F8_mN2


      if (debug_HNV_hadrTens.ge.3) write(*,*) 'matrix1=',matrix1

      if (debug_HNV_hadrTens.ge.3) write(*,*) 'matrix2=',matrix2

      ! vector part
      matrix=matrix + metricTensor(a,b)* matrix1 &
           & - q_in(a)*  ( F1_mN * gamma(:,:,(b+8))    - (F2_mN2 * W_in(b)  + F3_mN2 * momen_in(b) )* gamma5  )   ! minus here


      ! axial part
      if (process_ID.ne.EM) then
         matrix=matrix &
              & + metricTensor(a,b)* matrix2  &
              &     - q_in(a)*( F5_mN * gamma(:,:,b) + F6_mN2 * W_in(b)  * unit4    )    &
              & + q_in(a)*q_in(b)*F8_mN2 * unit4

      end if
    end function vertex_WNDelta_DirConj
    !   -------------------------------------------------------------------







    !   -------------------------------------------------------------------
    ! N'' -->  N'   pion
    function vertex_NNpi(ppi) result(matrix)

      use constants, only: f_pi, g_A, ii
      use minkowski, only: slashed5

      real, dimension(0:3), intent(in)  :: ppi
      complex, dimension(1:4,1:4)  :: matrix

      if (debug_HNV_hadrTens.ge.3) write(*,'(A,4g12.5)') 'In NievesHadtonTensor_1pi: vertex_NNpi:  ppi_out= ', ppi

      matrix=ii*g_A/2./f_pi*slashed5(ppi)

    end function vertex_NNpi
    !   -------------------------------------------------------------------




    !   -------------------------------------------------------------------
    ! W N --> N''
    function vertex_WNN(mu,momen_in,q_in) result(matrix)

      use matrix_module, only: unit4
      use minkowski, only: gamma,gamma5
      use leptonicID
      use constants, only: mN

      integer, intent(in) :: mu
      real, dimension(0:3), intent(in) :: momen_in,q_in ! momenta of the incoming nucleon, incoming boson
      complex, dimension(1:4,1:4)  :: matrix

      matrix=0

      matrix= ( ff_nucleon(1) +ff_nucleon(2) )*gamma(:,:,mu) - ff_nucleon(2)/2/mN*( momen_in(mu) + q_in(mu) + momen_in(mu) )*unit4

      if (process_ID.ne.EM) then

         matrix= matrix + ff_nucleon(3)*gamma(:,:,mu+8) + ff_nucleon(4)/mN * q_in(mu)*gamma5 ! gA for proton is negative
      end if

      if (debug_HNV_hadrTens.ge.3) write(*,'(A,g12.5,A,4g12.5)') 'In NievesHadtonTensor_1pi: vertex_WNN:  Q2=', &
           & Q2, '   ff_nucleon= ', ff_nucleon

    end function vertex_WNN
    !   -------------------------------------------------------------------

  end function   NievesHadronTensor_1pi



end module NievesHadronTensor
