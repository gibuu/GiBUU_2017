!***************************************************************************
!****m* /lepton_kinematics_free_FULL
! NAME
! module lepton_kinematics_free_FULL
!
! PURPOSE
! calculates some kinematics for 1-pion production processes
! This module is used by the program test_bgr (test_bgr.f90)
!
! In this sense, this file is not a part of the "mainsteam" GiBUU code
! However, as far as it /testXsection directory,  the test_bgr.f90 does not compile ...
! So to work with it (that is to run test_bgr.f90)  move this file to the    .../init/neutrino/   directory
!
!***************************************************************************
module lepton_kinematics_free_FULL

  public :: minmaxcostheta_Enu, minmaxcostheta_Enu_nu,  &  ! costheta
       & minmaxNu_Enu, minmaxNu_Enu_ct, &                  !nu
       & minmaxE1_Enu_Epi, &                               ! E1
       & minmaxCosThetaWPi_W, &                            ! costhetaPiW
       & minmaxEpi_Enu, minmaxEpi_Enu_ct, PionEnergy_W_ctWpi, & ! Epi
       & minmaxt_Q2_W2, minmaxcosthetaqPi_Q2_W2_Epi, &
       & minmaxsl_Enu, minmaxTpi_Enu_sl, minmaxctnupi_Enu_sl_Epi




contains



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  cos of the outgoing lepton angle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine minmaxcostheta_Enu(Enu,Wmin,mN,ml,costhetamin,costhetamax)

    real, intent(in) :: Enu,Wmin,mN,ml
    real, intent(out) :: costhetamin,costhetamax

    costhetamin=-1.
    costhetamax=1.

  end subroutine minmaxcostheta_Enu




  subroutine minmaxcostheta_Enu_nu(Enu,nu,Wmin,mN,ml,costhetamin,costhetamax)

    real, intent(in) :: Enu,nu,Wmin,mN,ml
    real, intent(out)    :: costhetamin,costhetamax

    costhetamax=1.

    costhetamin=( 2.*Enu**2 + Wmin**2 - mN**2 - ml**2 - 2.*(mN+Enu)*nu )/2./Enu/( sqrt((Enu-nu)**2 - ml**2) )

    if (costhetamin>1.) then
       costhetamin=1.
       write(*,'(A,g10.3,A)') 'In minmaxcostheta_Enu_nu:  costhetamin=', costhetamin, ' > 1.,   is set to 1.'
    else
       if (costhetamin<-1.) then
          costhetamin=-1.
          write(*,'(A,g10.3,A)') 'In minmaxcostheta_Enu_nu:  costhetamin=', costhetamin, ' < -1.,   is set to -1.'
       end if
    end if

  end subroutine minmaxcostheta_Enu_nu



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  outgoing lepton energy E1  and  energy transferred  nu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine minmaxNu_Enu_ct(Enu,costheta,Wmin,mN,ml,numin,numax)

    implicit none
    real, intent(in) :: Enu,costheta,Wmin,mN,ml
    real, intent(out)    :: numin,numax

    real    :: kappa, a, b, c, Dis

    if (Enu<ml) then
       write(*,*) 'In subroutine minmaxNu_Enu_ct(Enu,costheta,Wmin,mN,ml,numin,numax)'
       write(*,'(2(A,g11.5),A)') 'Enu=' , Enu, 'is less then    ml= ' , ml, '    numax is set to  ml, reaction not possible'
       numin=ml
       numax=ml
    else

       kappa = -Wmin**2 +mN**2 +ml**2 -2.*Enu**2
       a= 4.*( mN**2 + 2*mN*Enu + Enu**2*(1.-costheta**2) )
       b= 4.*( kappa*(mN+Enu) + 2.*Enu**3*costheta**2 )
       c= kappa**2 -4.*Enu**2*(Enu**2-ml**2)*costheta**2

       Dis=b**2/4.-a*c
       if (Dis<0) then
          write(*,*) '# In subroutine minmaxNu_Enu_ct(Enu,costheta,Wmin,mN,ml,numin,numax)'
          write(*,*) '# Dis=', Dis, ' <0.   nu_min  cannot be determined'
          write(*,*) '# Enu=', Enu, '   costheta=', costheta, '   Wmin=', Wmin
          write(*,*) '        mN=', mN, '   ml=', ml
          numin=0
          return
       else
          if (costheta.ge.0) numin=( -b/2. - sqrt(Dis) )/a
          if (costheta.le.0) numin=( -b/2. + sqrt(Dis) )/a
       end if

       if (numin<ml) numin=ml

       numax=Enu-ml
       !write(10,'(4(A,g10.4))') '# In minmaxNu_Enu_ct   Enu=', Enu, '    costheta=', costheta, &
       !            & '    Wmin=', Wmin, '    numin=', numin


       if (numin>numax) then
          write(*,'(2(A,g11.5))') 'That is really interesting how you did manage to get numin=', &
               & numin, '   grosser als  numax=', numax, '    return'
          return
       else
       endif

    end if

  end subroutine minmaxNu_Enu_ct




  subroutine minmaxE1_Enu_ct(Enu,costheta,Wmin,mN,ml,E1min,E1max)

    implicit none
    real, intent(in) :: Enu,costheta,Wmin,mN,ml
    real, intent(out)    :: E1min,E1max
    real        :: numin, numax

    call minmaxNu_Enu_ct(Enu,costheta,Wmin,mN,ml,numin,numax)
    E1min=Enu-numax
    E1max=Enu-numin

  end subroutine minmaxE1_Enu_ct




  subroutine minmaxE1_Enu_Epi(Enu,Epi,Wmin,ml,mN,mpi,E1_min,E1_max)

    implicit none
    real, intent(in) :: Enu, Epi, Wmin, ml, mN, mpi
    real, intent(out) :: E1_min, E1_max

    real :: numin,numax

    E1_min=ml

    call minmaxNu_Enu(Enu,Wmin,mN,ml,numin,numax)

    !   ???????????????????
    E1_max = (Enu-numin)*0.7


  end subroutine minmaxE1_Enu_Epi



  subroutine minmaxNu_Enu(Enu,Wmin,mN,ml,numin,numax)

    implicit none
    real, intent(in) :: Enu,Wmin,mN,ml
    real, intent(out) :: numin,numax

    real :: kappa, a, b, c, Dis, numin1,mpi

    if (Enu<ml) then
       write(*,*) 'In subroutine minmaxNu_Enu_ct(Enu,costheta,Wmin,mN,ml,numin,numax)'
       write(*,'(2(A,g11.5),A)') 'Enu=' , Enu, 'is less then    ml= ' , ml, '    numax is set to  ml'
       numax=ml
       numin=ml
    else

       ! ------------------
       !  numin
       ! ------------------
       kappa = -Wmin**2 +mN**2 +ml**2 -2.*Enu**2
       a= 4.*( mN**2 + 2*mN*Enu)
       b= 4.*( kappa*(mN+Enu) + 2.*Enu**3)
       c= kappa**2 -4.*Enu**2*(Enu**2-ml**2)

       Dis=b**2/4.-a*c
       if (Dis<0) then
          write(*,*) '# In subroutine minmaxNu_Enu_ct(Enu,costheta,Wmin,mN,ml,numin,numax)'
          write(*,*) '# Dis=', Dis, ' <0.   nu_min  cannot be determined'
          write(*,*) '# Enu=', Enu,  '   Wmin=', Wmin
          write(*,*) '        mN=', mN, '   ml=', ml
          numin=0
          return
       else
          numin=( -b/2. - sqrt(Dis) )/a
          write(10,'(3(A,g10.4))') ' In minmaxNu_Enu:  Enu=',Enu, '   numin=', numin, '   ml=',ml

          ! ?? 10 Oct 2008
          ! from s1
          mpi=0.139
          numin1 = Enu - (  ( sqrt(mN**2+2.*mN*Enu)-mN )**2 - mpi**2 - ml**2   )/2./mpi
          write(10,'(3(A,g10.4))') 'This corresponds to E1max=', (Enu-numin), 'to be compared with ', numin1
       end if

       ! ------------------
       !  numax
       ! ------------------
       numax=Enu-ml
       !write(10,'(4(A,g10.4))') '# In minmaxNu_Enu_ct   Enu=', Enu, '    costheta=', costheta, &
       ! & '    Wmin=', Wmin, '    numin=', numin
    end if
  end subroutine minmaxNu_Enu








  real function PionEnergy_W_ctWpi(W,ctWpi,mN1,mpi,NumRoots,PionEnergy1,PionEnergy2)
    ! two solutions are physical for some kinematic (for the significant Lorentz boosts from lab to CM frame )

    use minkowski, only: abs4Sq

    implicit none

    real, dimension(0:3), intent(in)    :: W
    real, intent(in)            :: ctWpi, mN1, mpi
    integer, intent(out)            :: NumRoots         ! number of solutions for pion energy
    real, intent(out)       :: PionEnergy1, PionEnergy2

    real    :: E_W, p_W, W2, a,b,c,dd, cosThetaWPi_min, cosThetaWPi_max

    E_W=W(0)
    p_W=sqrt( Dot_Product(W(1:3),W(1:3)) )
    W2=abs4Sq(W)
    !! write(10,('(A,g10.4)')) ' sqrt(W2)= ', sqrt(W2)
    NumRoots=1
    PionEnergy1=0.
    PionEnergy1=0.


    a=4.*(E_W**2 - p_W**2*ctWpi**2)
    b=-4.*E_W*( W2 -mN1**2 +mpi**2 )
    c=(W2-mN1**2+mpi**2)**2 +4.*mpi**2 *p_W**2 *ctWpi**2
    dd=b*b/4.-a*c
    !   dd=4.*p_W**2*ctWpi**2*(  ( W2 -mN1**2 +mpi**2 )**2 - 4.*(E_W**2-p_W**2*ctWpi**2)*mpi**2  )

    if (abs(dd)<1.e-10) then
       PionEnergy1= -b/2./a
    else
       if (dd>0) then
          ! sign(A,B)  If B > 0 then the result is ABS(A), else it is -ABS(A)
          PionEnergy1=( -b/2.+sign(sqrt(dd),ctWpi) )/a

          ! two solutions are possible for the kinematics that costheta_Wpi_min > 0
          call minmaxCosThetaWPi_W(W,mN1,mpi,cosThetaWPi_min,cosThetaWPi_max)
          if (cosThetaWPi_min > 0) then
             NumRoots=2
             PionEnergy1=( -b/2. - sqrt(dd) )/a
             PionEnergy2=( -b/2. + sqrt(dd) )/a
          end if
       else
          !write(10,'(6(A,g10.5),A)') '# In PionEnergy_W_ctWpi:    EW=',E_W, '    pW=',p_W, '   W=', (sqrt(W2)), &
          !                       &'    mN1=', mN1, '   mpi=',mpi, '    dd=',dd, '    <0 which correspond to W<mN1+mpi'
          PionEnergy1=0.
       end if
    end if
    PionEnergy_W_ctWpi=PionEnergy1

  end function PionEnergy_W_ctWpi









  subroutine minmaxCosThetaWPi_W(W,mN1,mpi,cosThetaWPi_min,cosThetaWPi_max)
    ! free nucleon kinematics
    use minkowski, only : abs4Sq, SP

    implicit none
    real, dimension(0:3), intent(in) :: W
    real, intent(in) :: mN1, mpi
    real, intent(out) :: cosThetaWPi_min,cosThetaWPi_max

    real :: cos2ThetaWPi_min, pW2, W2

    W2=SP(W,W)
    pW2=dot_product(W(1:3),W(1:3))

    ! forward scattering is always possible
    cosThetaWPi_max=1.

    ! the same as case 3 in Byckling--Kajantie
    cos2ThetaWPi_min= ( 4.*mpi**2*W(0)**2  -  (W2 -mN1**2 +mpi**2)**2 )/ ( 4.*mpi**2*pW2 )

    if ((cos2ThetaWPi_min.ge.0.) .and. (cos2ThetaWPi_min.le.1.)) then
       cosThetaWPi_min = sqrt(cos2ThetaWPi_min)
    else
       cosThetaWPi_min=-1.
    endif

  end subroutine minmaxCosThetaWPi_W









  subroutine minmaxEpi_Enu(Enu,mN,mpi,mN1,ml,Epimin,Epimax)

    implicit none
    real, intent(in) :: Enu,mN,mpi,mN1,ml
    real, intent(out)    :: Epimin,Epimax

    real    :: kappa, a, b, c, Dis, slmin

    if (Enu<0.25) then
       write(*,*) 'In subroutine minmaxEpi_Enu(Enu,mN,ml,Epimin,Epimax)'
       write(*,'((A,g11.5),A)') 'Enu=' , Enu, 'is less then  1-pion production threshold = 0.25,    Epimax is set to  mpi'
       Epimax=mpi
       Epimin=mpi
    else

       ! ------------------
       !  Epimax   from the requirement to reach sl_min
       ! ------------------
       slmin=(mN1+ml)**2
       kappa = -slmin +mN**2 +mpi**2 -2.*Enu**2
       a= 4.*( mN**2 + 2*mN*Enu)
       b= 4.*( kappa*(mN+Enu) + 2.*Enu**3)
       c= kappa**2 -4.*Enu**2*(Enu**2-mpi**2)

       Dis=b**2/4.-a*c
       if (Dis<0) then
          write(*,'(A,2(A,g12.5))') '# In minmaxEpi_Enu(Enu,mN,mpi,mN1,ml,Epimin,Epimax)    Dis=', Dis, &
               & ' <0.   Epi_max  cannot be determined'
          write(*,'(8(A,g12.5))') '# Enu=', Enu,  '   slmin=', slmin, &
               & '        mN=', mN, '   ml=', ml, '   mpi=',mpi, '    mN1=',mN1
          Epimax=mpi
          return
       else
          Epimax=( -b/2. + sqrt(Dis) )/a
          !   write(*,'(3(A,g10.4))') ' In minmaxEpi_Enu(Enu,mN,mpi,mN1,ml,Epimin,Epimax):  Enu=',Enu, '   Epimax=', Epimax, '   ml=',ml
       end if

       ! ------------------
       !  Epimin
       ! ------------------
       Epimin=mpi
    end if
  end subroutine minmaxEpi_Enu



  subroutine minmaxEpi_Enu_ct(Enu,costheta_nupi,mN,ml,mN1,mpi,Epimin,Epimax)
    ! costheta_nupi is angle between neutrino and pion
    implicit none
    real, intent(in) :: Enu,costheta_nupi,mN,ml,mN1,mpi
    real, intent(out) :: Epimin,Epimax

    real :: kappa, a, b, c, Dis, slmin

    if (Enu<mpi) then
       write(*,*) 'In subroutine minmaxEpi_Enu_ct(Enu,costheta_nupi,mN,ml,mN1,mpi,Epimin,Epimax)'
       write(*,'(2(A,g11.5),A)') 'Enu=' , Enu, 'is less then    mpi= ' , mpi, '    Epimax is set to  mpi, reaction not possible'
       Epimax=mpi
    else
       slmin=(mN1+ml)**2
       kappa = -slmin +mN**2 +mpi**2 -2.*Enu**2
       a= 4.*( mN**2 + 2*mN*Enu + Enu**2*(1.-costheta_nupi**2) )
       b= 4.*( kappa*(mN+Enu) + 2.*Enu**3*costheta_nupi**2 )
       c= kappa**2 -4.*Enu**2*(Enu**2-mpi**2)*costheta_nupi**2

       Dis=b**2/4.-a*c
       if (Dis<0) then
          write(*,'(A,g12.5,A)') '# In subroutine minmaxEpi_Enu_ct(Enu,costheta_nupi,mN,ml,mN1,mpi,Epimin,Epimax)  Dis=', Dis, &
               & ' <0.   Epi_max  cannot be determined'
          write(*,'(8(g12.5,A))')  '  Enu=', Enu, '   costheta_nupi=', costheta_nupi, '   slmin=', slmin, &
               '        mN=', mN, '   ml=', ml,  '   mpi=',mpi, '    mN1=',mN1
          Epimax=mpi
          return
       else
          if (costheta_nupi.ge.0) Epimax=( -b/2. - sqrt(Dis) )/a
          if (costheta_nupi.le.0) Epimax=( -b/2. + sqrt(Dis) )/a
       end if

       if (Epimax<mpi) Epimax=mpi

       Epimin=mpi

       if (Epimin>Epimax) then
          write(*,'(2(A,g11.5))') 'That is really interesting how you did manage to get Epimax=', &
               & Epimax, '   bigger als  Epimin=', Epimin, '    return'
          return
       else
       endif

    end if
  end subroutine minmaxEpi_Enu_ct










  subroutine minmaxt_Q2_W2(mN,Q2,W2,mN1,mpi,tmin,tmax)

    implicit none
    real, intent(in) :: mN, Q2, W2, mN1, mpi
    real, intent(out) :: tmin,tmax

    real :: a,b,c,dd

    a = W2
    b = W2**2 - W2*(mN1**2-Q2+mN**2+mpi**2) + (Q2+mN**2)*(mN1**2-mpi**2)
    c = W2*(mN1**2-mN**2) + mN1**2*Q2*(Q2-mN1**2+mN**2+mpi**2) + mN**2*mpi**2*(mpi**2+mN**2-mN1**2+Q2)
    dd = b**2/4.-a*c

    if (dd>1.e-8) then
       dd=sqrt(dd)
    else
       if (abs(dd)<1.e-8) then
          dd=0
       else
          write(*,*) 'In minmaxt_Q2_W2   the dis<0, which correspond to Q2 and/or W2 out of', &
               & 'kinematically allowed region. ', &
               & '    Current values Q2=', Q2, '    W2=', W2, '    return'
          tmin=0
          tmax=0
          return
       end if
    end if

    tmin=(-b/2.-dd)/a
    tmax=(-b/2.+dd)/a

    write(*,'(8(A,g12.5))') 'In minmaxt_Q2_W2: Q2=', Q2, '   W2=',W2, &
         & '   a= ', a, '    b=',b, '    dd=',dd, '   tmin=',tmin, '   tmax=',tmax

  end subroutine minmaxt_Q2_W2




  ! from minmax t
  subroutine minmaxcosthetaqPi_Q2_W2_Epi(mN,Q2,W2,Epi,mN1,mpi,ctmin,ctmax)

    implicit none
    real, intent(in) :: mN, Q2, W2, Epi, mN1, mpi
    real, intent(out) :: ctmin,ctmax

    real :: nu, qz, tmin,tmax

    nu=(W2+Q2-mN**2)/2./mN
    qz=sqrt(Q2+nu**2)

    call minmaxt_Q2_W2(mN,Q2,W2,mN1,mpi,tmin,tmax)

    ctmin = ( -Q2 + 2.*nu*Epi+mpi**2 - tmax ) / ( 2.*qz*sqrt(Epi**2-mpi**2) )
    ctmax = ( -Q2 + 2.*nu*Epi+mpi**2 - tmin ) / ( 2.*qz*sqrt(Epi**2-mpi**2) )

  end subroutine minmaxcosthetaqPi_Q2_W2_Epi





  ! from minmax t
  ! to do : to be compared with minmaxEpi(W(0:3))
  subroutine minmaxEpi_Q2_W2(mN,Q2,W2,mN1,mpi,Epimin,Epimax)

    implicit none
    real, intent(in) :: mN, Q2, W2, mN1, mpi
    real, intent(out) :: Epimin, Epimax

    real :: tmin,tmax

    call minmaxt_Q2_W2(mN,Q2,W2,mN1,mpi,tmin,tmax)

    Epimin=mpi
    Epimax=W2/4.   ! 2009-07-28 strange formula

  end subroutine minmaxEpi_Q2_W2




  subroutine minmaxsl_Enu(ml,mN,Enu,mN1,mpi,slmin,slmax)
    !sl=(p1N+k1)^2 - complete analoge of W2
    implicit none
    real, intent(in) :: ml, mN, Enu, mN1, mpi
    real, intent(out) :: slmin,slmax

    slmin = (ml+mN1)**2
    slmax = (sqrt(2.*mN*Enu+mN**2) -mpi)**2

  end subroutine minmaxsl_Enu




  subroutine minmaxTpi_Enu_sl(mN,Enu,sl,mN1,mpi,Tpimin,Tpimax)
    ! Tpi = mpi**2-tpi = mpi**2  -  (k-ppi)_mu \cdot (k-ppi)^mu

    implicit none

    real, intent(in) :: mN, Enu,sl,mN1,mpi ! sl=(p1N+k1)^2 - complete analoge of W2
    real, intent(out) :: tpimin,tpimax

    !
    real :: s, stilde, dis

    stilde=2.*mN*Enu
    s=stilde+mN**2

    dis = (sl-s-mpi**2)**2 - 4.*s*mpi**2

    if (dis.ge.1.e-8) then
       dis=sqrt(dis)
    else
       if (abs(dis).le.1.e-8) then
          dis=0.
       else
          write(*,'(6(A,g12.5))') 'In minmaxTpi_Enu_sl dis<0, which correspond to sl=(p1N+k1)^2<sl_max    mN=', &
               & mN, '   Enu=', Enu, '    sl=', sl, '   mN1=', mN1, '    mpi=', mpi, 'return'
          dis=9999999.
          return
       end if
    end if

    Tpimin = stilde/2./s*( -(sl-s-mpi**2) - dis )
    Tpimax = stilde/2./s*( -(sl-s-mpi**2) + dis )

  end subroutine minmaxtpi_Enu_sl



  ! from minmaxTpi
  subroutine minmaxctnupi_Enu_sl_Epi(mN,Enu,sl,mN1,mpi,Epi,ctmin,ctmax)
    implicit none
    real, intent(in) :: mN,Enu,sl,mN1,mpi,Epi
    real, intent(out) :: ctmin,ctmax

    real :: Tpimin,Tpimax

    call minmaxTpi_Enu_sl(mN,Enu,sl,mN1,mpi,Tpimin,Tpimax)

    ctmin = (Enu*Epi - Tpimax/2.) / Enu/ sqrt(Epi**2-mpi**2)
    ctmax = (Enu*Epi - Tpimin/2.) / Enu/ sqrt(Epi**2-mpi**2)

  end subroutine minmaxctnupi_Enu_sl_Epi



end module lepton_kinematics_free_FULL
