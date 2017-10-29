!******************************************************************************
!****m* /hadronTensor_QE
! NAME
! module hadronTensor_QE
! PURPOSE
! * Evaluates the hadron tensor for gamma N -> N'
! * For details see the notes about this in the PhD of Oliver Buss
!******************************************************************************

module hadronTensor_QE
  implicit none

  private

  public :: H_munu_QE

  logical, save :: initFlag=.true.
  logical, save :: debug=.false.


  !****************************************************************************
  !****g* hadronTensor_QE/extraTerm
  ! SOURCE
  logical, save :: extraTerm=.true.
  !
  ! PURPOSE
  ! Switch to take out the extra term which preserves charge conservation in
  ! case of a momentum dependend mass.
  !****************************************************************************


contains

  subroutine initInput
    use output

    integer :: ios

    NAMELIST /hadronTensorQE/ extraTerm,debug

    call Write_ReadingInput('hadronTensorQE',0)

    rewind(5)
    read(5,nml=hadronTensorQE,IOSTAT=ios)
    call Write_ReadingInput("hadronTensorQE",0,ios)

    write(*,*) 'Charge conservation term included?', extraTerm
    write(*,*) 'Debugging?', debug

    call Write_ReadingInput('hadronTensorQE',1)

  end subroutine initInput


  !****************************************************************************
  !****f* hadronTensor_QE/H_munu_QE
  ! NAME
  ! real function H_munu_QE(mu,nu,pin,pout,nucleonCharge)
  !
  ! PURPOSE
  ! calculate the hadronic tensor for QE scattering
  ! INPUTS
  ! * integer :: nu, mu                -- Lorentz indices
  ! * real, dimension(0:3) :: pin,pout -- incoming/outgoing nucleon momentum
  ! * integer :: nucleonCharge         -- charge of nucleon
  ! NOTES
  ! This routine is highly inefficient, since it recalculates many things,
  ! even if only the lorentz indizes have changed sinc the last call.
  !****************************************************************************
  real function H_munu_QE(mu,nu,pin,pout,nucleonCharge)
    use minkowski, only: metricTensor,SP
    use FF_QE_nucleonScattering, only: formfactors_QE
    use constants, only: mN, electronChargeSQ
    use leptonicID, only: EM

    integer,intent(in) :: nu, mu
    real, intent(in),dimension(0:3) :: pin,pout
    integer, intent(in) :: nucleonCharge

    real, dimension(0:3) :: q        ! incoming virtual photon momentum
    real :: F1, F2                   ! Form factors
    real:: alpha
    real, dimension(0:3)::  beta
    real :: mi, mf ! final and initial masses of the nucleons

    if (initFlag) then
       call initInput()
       initFlag=.false.
    end if

    q=pout-pin

    call formfactors_QE(-SP(q,q),EM,nucleonCharge,F1,F2)

    if ((max(mu,nu).gt.3).or.(min(mu,nu).lt.0)) then
       write(*,*) 'Error in hadronic tensor for QE: indices not well defined'
       write(*,*) mu,nu
       stop
    end if

    if (Dot_Product(pin+q-pout,pin+q-pout).gt.0.0001) then
       write(*,*) 'Error in hadronic tensor for QE: momentum not conserved'
       write(*,*) pin
       write(*,*) '+', q,  '=>'
       write(*,*) pout
       stop
    end if

    mi=sqrt(SP(pin,pin))
    mf=sqrt(SP(pout,pout))
    if (mi.lt.0.000001) then
       write(*,*) 'mi is 0 in  hadronTensor_QE', mi
    end if
    if (mf.lt.0.000001) then
       write(*,*) 'mF is 0 in  hadronTensor_QE', mf
    end if


    alpha = F1 + F2*(mi+mf)/(2*mN)
    if (extraTerm) then
       beta = F1*(mf-mi)/Sp(q,q)*q - F2*(pin+pout)/(2*mN)
    else
       beta = -F2*(pin+pout)/(2*mN)
    end if
    if (debug) then
       if (nucleonCharge.eq.1) then
          write(101,'(2I3,4F14.7)') mu,nu,F1* (mi-mf)/Sp(q,q)* q /(-1./(2.*mN) *F2* (pin+pout))
          write(201,'(3F14.7)') mi-mf,mi,mf
       else
          write(102,'(2I3,4F14.7)') mu,nu,F1* (mi-mf)/Sp(q,q)* q /(-1./(2.*mN) *F2* (pin+pout))
          write(202,'(3F14.7)') mi-mf,mi,mf
       end if
    end if

    H_munu_QE= ( &
         & alpha**2*(pout(mu)*pin(nu) + pout(nu)*pin(mu) &
         & - metricTensor(mu,nu)*SP(pin,pout) + mi*mf*metricTensor(mu,nu))  &
         & + beta(mu)*beta(nu)* (SP(pin,pout)+mi*mf)  &
         & + alpha   *beta(nu)* (mf*pin(mu)  +mi*pout(mu)) &
         & + beta(mu)*alpha   * (mf*pin(nu)  +mi*pout(nu)) &
         & ) / (2*mi*mf) * electronChargeSQ

  end function H_munu_QE
end module hadronTensor_QE
