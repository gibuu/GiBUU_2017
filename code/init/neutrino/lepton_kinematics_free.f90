!******************************************************************************
!****m* /lepton_kinematics_free
! NAME
! module lepton_kinematics_free
!
! PURPOSE
! calculates some kinematics for 1-pion production processes
! This module is used by modules neutrinoXsection, singleParticleProductionHNVlike,
! and by the program test_bgr (test_bgr.f90)
!******************************************************************************
module lepton_kinematics_free

implicit none
private

public ::   minmaxE1_costheta,  &    ! minmaxE1_Q2
          & minmaxEpi_W,  &          ! Epi
          & CosThetaPiW_W_Epi, &     ! theta_PiW
          &  lowBoundary_for2costhetaPiSolutions



contains


real function Wmin_IP(IP,mN)
!use particleProperties, only : baryon
  integer, intent(in)  :: IP ! ID of the process, 1=QE, 2-Delta, 3-31 highRES, 34=DIS
  real, intent(in)     :: mN ! mass of the incoming nucleon

  select case (IP)
      case (1)
          Wmin_IP=mN*0.95  ! QE
      case (2:31)
          !Wmin_IP=baryon(IP)%minmass    !  minmass of the Delta
          Wmin_IP=mN+0.14
      case (32:33)
          !Wmin_IP=baryon(2)%minmass     ! minmass of the corresponding resonance
          Wmin_IP=mN+0.14
      case (34)
          Wmin_IP=1.6
      case default
          write(*,*) 'Wrong IP. Minimal invariant mass cannot be determined. STOP.   IP=', IP
          stop
    end select

end function Wmin_IP



! real function Wmax_IP(IP,Enu,mN,ml)
!   integer, intent(in)  :: IP ! ID of the process, 1=QE, 2-Delta, 3-31 highRES, 34=DIS
!   real, intent(in)     :: Enu, mN, ml ! mass of the incoming nucleon
!
!   real  :: Wmax_Enu
!
!   Wmax_Enu=sqrt(mN*(mN+2*Enu))-ml
!
!   select case (IP)
!       case (1)
!           Wmax_IP=min( Wmax_Enu , mN*1.05 )     ! QE
!       case (2:33)
!           Wmax_IP=min( Wmax_Enu , 2.0 )
!       case (34)
!           Wmax_IP=Wmax_Enu
!       case default
!           write(*,*) 'Wrong IP. Minimal invariant mass cannot be determined. STOP.   IP=', IP
!           stop
!     end select
!
! end function Wmax_IP



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  min and  max energies of the outgoing lepton depending on Enu and Q2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine minmaxE1_Q2(IP,Enu,mN,ml,Q2,E1min,E1max, W_min, W_max, success)
!   integer, intent(in)           :: IP ! ID of the process, 1=QE, 2-Delta, 3-31 highRES, 34=DIS
!   real,    intent(in)           :: Enu, mN, ml, Q2
!   real, optional, intent(in)        :: W_min, W_max
!   real, intent(out)                 :: E1min, E1max
!   logical, optional, intent(out)    :: success
!
!   real    :: W2max_analyt, Wmin, Wmax, numin, numax, ss
!   logical :: success1
!
!   success1=.false.
!
! ! minimal invariant mass to be reached
!   if (present(W_min)) then
!     Wmin=W_min
!   else
!     Wmin=Wmin_IP(IP,mN)
!   end if
!
!
! ss=mN*(2.*Enu+mN)
! W2max_analyt=(ss-Q2-mN**2-ml**2)*(ss*Q2+mN**2*ml**2)/(ss-mN**2)/(Q2+ml**2)
!
!   if (present(W_max)) then
!     ! if W_max is given, choose min from kinematically allowed and given
!     Wmax= min( W_max, sqrt(W2max_analyt) )
!   else
!     ! if W_max is not given, choose min from kinematically allowed and  default for the process
!     Wmax=min( Wmax_IP(IP,Enu,mN,ml),  sqrt(W2max_analyt) )
!   end if
!
! numax=(Wmax**2-mN**2+Q2)/2/mN
! E1min=Enu-numax
!
! E1min=max(ml*1.001,0.8*E1min)
!
! numin=(Wmin**2-mN**2+Q2)/2/mN
! E1max=Enu-numin
! E1max=min(Enu,1.2*E1max)
!
! if (E1min.ge.ml .and. E1max.ge.E1min) success1=.true.
!
! if (present(success)) success=success1
!
! end subroutine minmaxE1_Q2




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  min and  max energies of the outgoing lepton depending on Enu and costheta_of_this_lepton
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine minmaxE1_costheta(IP,Enu,mN,ml,costheta,E1min,E1max, W_min, W_max, success)
  integer, intent(in)               :: IP ! ID of the process, 1=QE, 2-Delta, 3-31 highRES, 34=DIS
  real,         intent(in)          :: Enu, mN, ml, costheta
  real, optional, intent(in)                 :: W_min, W_max
  real, intent(out)                 :: E1min, E1max
  logical, optional, intent(out)    :: success

  real    :: Wmin, a1, a2, b, a3, numin
  logical :: success1

  success1=.false.


! minimal invariant mass to be reached
  if (present(W_min)) then
    Wmin=W_min
  else
    Wmin=Wmin_IP(IP,mN)
  end if



  a2=Wmin**2 -mN**2 -ml**2
  a1=a2+2*Enu*mN
  b=Wmin**2 -2*Enu*mN -mN**2 +ml**2
  a3=b**2 - 4*Wmin**2*ml**2 - 4*ml**2*Enu**2*(1.-costheta**2)

  numin= ( Enu*( 2*Enu**2*(1.-costheta**2) + a1) + mN*a2 - Enu*costheta*sqrt(a3) ) &
                &  /2/(mN+Enu-Enu*costheta)/(mN+Enu+Enu*costheta)
  E1max=Enu-numin
  E1max=min(Enu,1.2*E1max)


  !(  E*(2*pow(E,2)*(1.-pow(costheta,2))+a1) + p0*a2
  !       -E*costheta*sqrt(
  !          pow(b,2) - 4*pow((W_min_E(E)),2)*pow(ml,2) -4*pow(ml,2)*pow(E,2)*(1.-pow(costheta,2))  )                    )/2/(p0+E-E*costheta)/(p0+E+E*costheta);

  E1min=ml

  if (E1max>E1min) success1=.true.

  if (present(success)) success=success1

end subroutine minmaxE1_costheta





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  angle between pion and resonance
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine CosThetaPiW_W_Epi(W,Epi,mN1,mpi,ct,success)
  use minkowski, only: abs4Sq
  use random, only: rn

  real, dimension(0:3), intent(in)  :: W
  real,         intent(in)          :: Epi, mN1, mpi
  real, intent(out)                 :: ct
  logical, optional, intent(out)                 :: success

    real    :: pW, ppi
    logical :: success1

  success1=.false.

  pW=sqrt( Dot_Product(W(1:3),W(1:3)) )
  ppi=sqrt(Epi**2-mpi**2)

  if ( abs(pW).eq.0 ) then
            write(*,*) 'Decaying particle is in its rest frame. Decay is isotropic.'
            write(*,*) 'CosTheta is generated randomly'
            ct = rn()
            success1=.true.
  else
    ct = ( 2*W(0)*Epi - abs4Sq(W) + mN1**2 - mpi**2 )/2./pW/ppi
  end if

! checking that cos is within -1 .. 1
  if (ct>1.) then
                !write (10,'(4(A,g10.5))') '# Strange CosThetaPi_W_Epi = ', ct, ' for   Epi=',Epi, '   EW=',(W(0)), '   pW=',pW
                !write (10,*) '# It means NO solution'
                ct=1.
             else if (ct<-1.) then
                !write (10,'(4(A,g10.5))') '# Strange CosThetaPi_W_Epi = ', ct, ' for   Epi=',Epi, '   EW=',(W(0)), '   pW=',pW
                !write (10,*) '# It means NO solution'
                ct=-1.
            else
                success1=.true.
  end if

  if (present(success)) success=success1

end subroutine CosThetaPiW_W_Epi



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!   pion energy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine minmaxEpi_W(W,mN1,mpi,Epi_min,Epi_max,success)
! for free nucleon kinematis
use minkowski, only : abs4Sq

    real, dimension(0:3),   intent(in)  :: W
    real,           intent(in)  :: mN1,mpi
    real,           intent(out) :: Epi_min,Epi_max
    logical, optional, intent(out)                 :: success


    real    :: E_W, p_W, W2
    real    :: a, b, c, dd!, dd1
    logical :: success1

 success1=.false.

 E_W=W(0)
 p_W=sqrt( Dot_Product(W(1:3),W(1:3)) )
 W2=abs4Sq(W)

  a = 4.*W2
  b = -4.*E_W*( W2 -mN1**2  + mpi**2 )
  c = (W2-mN1**2+mpi**2)**2 + 4.*p_W**2*mpi**2
  dd=b*b/4.-a*c
  ! dd >0 for W2>(mN1+mpi)^2
  if (dd>0) then
        dd=sqrt(dd)
        Epi_max=( -b/2.+dd )/a
        Epi_min = max(  ( -b/2.-dd )/a , mpi  )
        if (Epi_min>0  .and.  Epi_max>Epi_min) success1=.true.
  else
    !write(*,'(5(A,g10.4),A)') '# In minmaxEpi_W  for Epi_max:  EW=',E_W, '    pW=',p_W, '   W=',(sqrt(W2)), &
    !      &'   mN1=', mN1, '   mpi=',mpi, '    Dis<0 which correspond to W<mN1+mpi'
    Epi_max=mpi
    Epi_min=mpi
  end if

 if (present(success))  success=success1

end subroutine minmaxEpi_W








!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  others
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine lowBoundary_for2costhetaPiSolutions(W2,mN1,mpi,EWmin,numin,Q2min)
    ! the second solution appears for the angle between the outgoing nucleon and pion

    real, intent(in)        :: W2, mN1, mpi ! W2 is invariant mass squared, mN1 and mpi masses of the outgoing nucleon and pion
    real, intent(out)   ::  EWmin, numin, Q2min


    if (sqrt(W2).le.(mN1+mpi)) then
            write(*,*) 'In EW_boundary_for2costhetaPisolutions:  W=', sqrt(W2), 'You are below 1-pion production thereshold. '
            !stop
    else
    EWmin=(W2-mN1**2+mpi**2)/2./mpi  ! actually this correspond EW > Epi_in_CM_frame * W / mpi
                                     ! this also correspond to the minimal costhetaPiW=0, as expected from relativistic kinematics
    end if
    numin=EWmin-mN1                     ! strictly speaking, here there should be mN instead of mN1
    Q2min=mN1**2+2.*mN1*numin-W2        ! strictly speaking, here there should be mN instead of mN1

end subroutine lowBoundary_for2costhetaPiSolutions




end module lepton_kinematics_free
