!******************************************************************************
!****m* /FF_Delta_production
! NAME
! module FF_Delta_production
!
! PURPOSE
! Provides the electro-weak form factors for the process lepton nucleon -> lepton' Delta.
!******************************************************************************


module FF_Delta_production

  implicit none
  private

  logical, save :: initFlag=.true.

  !****************************************************************************
  !****g* FF_Delta_production/FF_Delta
  ! SOURCE
  !
  integer, save :: FF_Delta=1   ! 0: FF taken from MAID
                                ! (see notes of Luis Alvarez-Ruso)
                                ! 1: FF taken from Paschos
                                ! (hep-ph/0501109 v2)
  ! PURPOSE
  ! This switch decides whether the Paschos form factors (FF_Delta=1) or the
  ! Maid form factors (FF_Delta=0) are used. Default is FF_Delta=1.
  !****************************************************************************


  public :: formfactors_Delta


contains

  !****************************************************************************
  !****n* FF_Delta_production/input_FF_Delta
  ! NAME
  ! NAMELIST /input_FF_Delta
  ! PURPOSE
  ! Includes the input switch:
  ! * FF_Delta
  !****************************************************************************

  !****************************************************************************
  !****s* FF_Delta_production/readInputDeltaFormFactors
  ! NAME
  ! subroutine readInputDeltaFormFactors
  !
  ! PURPOSE
  ! decide, which set of form factors is used by reading in FF_Delta.
  !
  !****************************************************************************

  subroutine readInputDeltaFormFactors
    use output

    integer :: ios

    NAMELIST /input_FF_Delta/ FF_Delta
    call Write_ReadingInput('input_FF_Delta',0)
    rewind(5)
    read(5,nml=input_FF_Delta,IOSTAT=ios)
    call Write_ReadingInput("input_FF_Delta",0,ios)
    call Write_ReadingInput('input_FF_Delta',1)

    if (FF_Delta.eq.0) then
       write(*,*) 'Delta FF are taken from MAID'
    else if (FF_Delta.eq.1) then
       write(*,*) 'Delta FF are taken from Paschos'
    else
       write(*,*) 'strange FF_Delta -> STOP', FF_Delta
       stop
    end if

  end subroutine readInputDeltaFormFactors


  !****************************************************************************
  !****s* FF_Delta_production/formfactors_Delta
  ! NAME
  ! subroutine formfactors_Delta(Qs,W,processID,c3v,c4v,c5v,c6v,c3a,c4a,c5a,c6a)
  !
  ! PURPOSE
  ! * Provides the electro-weak form factors for the process lepton nucleon -> lepton' Delta.
  !
  ! INPUTS
  ! *  real, intent(in) ::    Qs ! =Q^2 : virtuality of gauge boson in units of GeV**2 = (-1)*Mandelstam-t
  ! *  real, intent(in) ::    W  ! invariant mass
  ! *  integer, intent(in) :: processID
  !
  ! "processID" specifies the reaction type: EM, CC, NC
  !
  ! OUTPUT
  ! *  real, intent(out) ::c3v,c4v,c5v,c6v             ! Vector Form factors
  ! *  real, intent(out),optional :: c3a,c4a,c5a,c6a   ! Axial Form factors
  !                                                     (not necessary as input if processID=EM)
  !
  !****************************************************************************
  subroutine formfactors_Delta(Qs,W,processID,c3v,c4v,c5v,c6v,c3a,c4a,c5a,c6a)
    use constants, only: sinsthweinbg
    use leptonicID

    real, intent(in) :: Qs ! =Q^2 : virtuality of gauge boson
    real, intent(in) :: W  ! invariant mass
    integer, intent(in) :: processID ! Specifies the reaction type

    real, intent(out) :: c3v,c4v,c5v,c6v   ! Vector Form factors
    real, intent(out),optional :: c3a,c4a,c5a,c6a   ! Axial Form factors

    !*** Read input:
    if (initFlag) then
       call readInputDeltaFormFactors
       initFlag=.false.
    end if

    !*** Check input:

    if (Qs.lt.0) then
       write(*,*) 'Qs less than zero in formfactors_Delta. STOP!', Qs
       stop
    end if


    if (FF_Delta.eq.0) then
    else if (FF_Delta.eq.1) then
    else
       write(*,*) 'nonsense input for FF_Delta in formfactors_Delta. STOP!', FF_Delta
       stop
    end if



    !*** set form factors

    c6v=0. ! Conserved vector current hypothesis



    select case (processID)

    case (CC,antiCC)
       if (.not.(present(c3a).and.present(c4a).and.present(c5a).and.present(c6a))) then
          write(*,*) 'Error in formfactors_Delta! one of the axial ff not defined'    &
               &  ,present(c3a),present(c4a),present(c5a),present(c6a),processID
          stop
       end if

       if (FF_Delta.eq.0) call vector_FF_Maid(Qs,W,c3v,c4v,c5v)
       if (FF_Delta.eq.1) call vector_FF_Paschos(Qs,W,c3v,c4v,c5v)

       call axial_FF_Paschos(Qs,W,c3a,c4a,c5a,c6a)
       !the axial form factors need to be refitted when used with Maid vector form factors!!!
       !this has not yet been done

    case (NC,antiNC)
       if (.not.(present(c3a).and.present(c4a).and.present(c5a).and.present(c6a))) then
          write(*,*) 'Error in formfactors_Delta! one of the axial ff not defined'    &
               &   ,present(c3a),present(c4a),present(c5a),present(c6a),processID
          stop
       end if

       if (FF_Delta.eq.0) call vector_FF_Maid(Qs,W,c3v,c4v,c5v)
       if (FF_Delta.eq.1) call vector_FF_Paschos(Qs,W,c3v,c4v,c5v)

       c3v=c3v*(1.-2.*sinsthweinbg)
       c4v=c4v*(1.-2.*sinsthweinbg)
       c5v=c5v*(1.-2.*sinsthweinbg)

       call axial_FF_Paschos(Qs,W,c3a,c4a,c5a,c6a)
       !the axial form factors need to be refitted when used with Maid vector form factors!!!
       !this has not yet been done

    case (EM,antiEM)
       if (present(c3a)) c3a=0.
       if (present(c4a)) c4a=0.
       if (present(c5a)) c5a=0.
       if (present(c6a)) c6a=0.

       if (FF_Delta.eq.0) call vector_FF_Maid(Qs,W,c3v,c4v,c5v)
       if (FF_Delta.eq.1) call vector_FF_Paschos(Qs,W,c3v,c4v,c5v)


    case default
       write(*,*) 'Error in formfactors_Delta! Invalid process ID:', processID
    end select


  end subroutine formfactors_Delta




  subroutine axial_FF_Paschos(Qs,W,c3a,c4a,c5a,c6a)
    use constants, only: mN, mPi

    real, intent(in) :: Qs, W
    real, intent(out) :: c3a,c4a,c5a,c6a

    real :: MA

    !formfactors as in Paschos, hep-ph/0501109 v2

    MA=1.05

    c5a=1.2/((1.+Qs/(3.*MA**2))*(1.+Qs/MA**2)**2)
    c4a=-c5a/4.
    c3a=0.
    c6a=c5a*mN**2/(mPi**2+Qs)

  end subroutine axial_FF_Paschos


  subroutine vector_FF_Paschos(Qs,W,c3v,c4v,c5v)
    use constants, only: mN

    real, intent(in) :: Qs, W
    real, intent(out) :: c3v,c4v,c5v

    real :: MV

    !formfactors as in Paschos, hep-ph/0501109 v2

    MV=0.84

    c3v=1.95/((1.+Qs/(4.*MV**2))*(1.+Qs/MV**2)**2)
    c5v=0.
    c4v=-mN/W*c3v

  end subroutine vector_FF_Paschos



  subroutine vector_FF_Maid(Qs,W,c3v,c4v,c5v)

    !     This routine calculates the vector N-Delta ff for a given Qs
    !     from the helcicity amplitudes of Tiator et al. (nucl-th/0310041)
    !     input in GeV(^2) output dimensionless

    !     subroutine written by Luis Alvarez-Ruso

    use constants, only: alphaQED, pi

    real,intent(in)::Qs,w   ! Qs and Delta inv. mass (in GeV)
    real,intent(out)::c3v,c4v,c5v

    real,parameter::mdel=1.232       !  delta mass
    real,parameter::mn=.9382723   ! taken from MAID

    real::kgcm0,egcm,qcm
    real::a12,a32,s12         ! helicity amplitudes
    real::r,kr,mpw2,mmw2,q2,Q2G

    Q2G=Qs
    q2=-Qs

    kgcm0 = (mdel*mdel-mn*mn)/2./mdel
    egcm = (mdel*mdel-Q2G-mn*mn)/2./mdel

    if (egcm**2+Q2G.lt.0.) then
       write(*,*) 'egcm**2+Q2G.lt.0. in Maid form factors for the Delta -> STOP', egcm**2+Q2G
       stop
    end if

    qcm=sqrt(egcm**2+Q2G)

    kr=(w**2-mn**2)/2./w
    mpw2=(mn+w)**2
    mmw2=(mn-w)**2

    CALL HP33(Q2G,qcm,kgcm0,a12,a32,s12)


    a12=a12*1.e-3
    a32=a32*1.e-3
    s12=s12*1.e-3
    !       M1 limit
    !        s12=0.
    !        a32=sqrt(3.)*a12

    if (((mpw2-q2)*(mmw2-q2)).lt.0.or.(6.*kr/pi/alphaQED*mn*w/(mmw2-q2)).lt.0.) then
       !write(*,*) 'sqrt problems in Maid form factors for the Delta'
       c3v=0.
       c4v=0.
       c5v=0.
       return
    end if

    r=sqrt(6.*kr/pi/alphaQED*mn*w/(mmw2-q2))

    c3v=r*mn*w/(mpw2-q2)*(a12+a32/sqrt(3.))

    c4v=2.*r*mn**2/(mpw2-q2)/(mmw2-q2)*(mn*w*a12 &
         & -(w**2+mn**2-mn*w-q2)*a32/sqrt(3.) &
         & - sqrt(2.)*w**2*(w**2-mn**2-q2)/sqrt((mpw2-q2)*(mmw2-q2))*s12)

    c5v=2.*r*(w*mn)**2/(mpw2-q2)/(mmw2-q2)*(-a12+a32/sqrt(3.) &
         & +sqrt(2.)*(w**2-mn**2+q2)/sqrt((mpw2-q2)*(mmw2-q2))*s12)

    ! to be consistent with the adopted signs
    c3v=-c3v
    c4v=-c4v
    c5v=-c5v

  contains

    subroutine HP33(Q2G,qcm,qcm0,A1,A3,S1)  ! taken from SAID

      real,intent(in)::Q2G,qcm,qcm0
      real,intent(out)::A1,A3,S1
      real::Fq,A10,A30,S10

      Fq=exp(-0.21*Q2G)/(1+Q2G/0.71)**2*(qcm/qcm0)

      ! old         A10=-140.250
      ! old          A30=-265.437
      ! old          S10=  27.577

      !MAID 2005
      A10=-137.445
      A30=-260.128
      S10=  27.577

      A1= A10*(1.+ 0.0214486*Q2G)*Fq
      A3= A30*(1.- 0.0065431*Q2G)*Fq
      S1= S10*(1.+ 0.0166667*Q2G**3)*Fq

    end subroutine HP33

  end subroutine vector_FF_Maid

end module FF_Delta_production
