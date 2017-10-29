!***************************************************************************
!****m* /lepton_xsec_free
! NAME
! module lepton_xsec_free
!
! PURPOSE
! calculates the contribution of various diagrams (Dp, cDp, Np, cNp, CT, pp, pF)
! to the cross section of 1-pion production  on a nucleon
! (  see Hernandez Nieves Valverde PRC 76  )
! This module is used by the program test_bgr (test_bgr.f90)
!
! In this sense, this file is not a part of the "mainsteam" GiBUU code
! However, as far as it /testXsection directory,  the test_bgr.f90 does not compile ...
! So to work with it (that is to run test_bgr.f90)  move this file to the    .../init/neutrino/   directory
!
!***************************************************************************

module lepton_xsec_free

  public :: HNV_free_elepton_ct_ctPi_phiPi,  HNV_free_elepton_ct_ctPi, HNV_free_elepton_ct, HNV_free_ctPi


contains


  ! everything is for p pi+ final state
  ! 4-fold differential xsec dsi / dE1 dcostheta dcosthetaPi dphiPi
  subroutine HNV_free_elepton_ct_ctPi_phiPi(process_ID,Enu,E1,ct,ctPi,phiPi,mN,ml,mN1,mpi, &
       &  xsec_cDp,xsec_Np,xsec_cNp,xsec_CT,xsec_pp,xsec_pF)

    use lepton_kinematics_free_FULL, only : PionEnergy_W_ctWpi
    use rotation
    use minkowski, only : SP, abs4Sq
    use idtable, only : rho, nucleon, Delta
    use particleProperties, only : meson, baryon
    use leptonicID, only : CC
    use FF_QE_nucleonScattering, only : formfactors_QE
    use formFactor_ResProd, only : getFormfactor_Res
    use constants, only : g_A, f_pi, GF, pi, coscab, alphaQED
    use electronPionProduction_kine, only : get_k_abs_improved


    use particleDefinition


    implicit none

    integer, intent(in) :: process_ID
    real, intent(in) :: Enu, E1, ct, ctPi, phiPi, mN, ml, mN1, mpi
    real, intent(out) :: xsec_cDp, xsec_Np, xsec_cNp, xsec_CT, xsec_pp, xsec_pF


    real :: pk, p1k, pq, p1q, p1p, Q2
    real :: Epi, st, k1, ppi, stPi, cphiPi, sphiPi
    real, dimension(1:2) :: Epir ! two pion energies in case of two solutions

    real :: ct_Wpi
    real, dimension(1:3) :: unit_ppi_out, unit_W

    real, dimension(0:3) :: k_in, k_out, q, p, p_out, W, ppi_out

    real :: mrho2, t, vecW ! t is the momentum of the rho-meson squared

    real :: f1V, f2V, f1A, C5A, F_rho, FV_CT, F_PF,     M2cDp, M2Np, M2cNp, M2CT, M2pp, M2pF, kinemat_factor
    real :: MR

    real, dimension(1:4) :: ff_nucleon
    real, dimension(1:8) :: ff_Delta
    logical :: success_ff


    ! for comparison with kinematics of Oliver
    logical :: debug=.false.
    type(particle) :: nucleon_in
    logical :: success
    integer :: NumRoots, r, r1
    real :: p_pion_abs2
    real, dimension(1:2,0:3) :: p_pion_out_root
    real, dimension(1:3) :: betaTOCM


    st=sqrt((1.-ct)*(1.+ct))
    k1=sqrt((E1-ml)*(E1+ml))

    k_in=(/ Enu, 0., 0., Enu /)
    p=(/ mN, 0., 0., 0. /)
    k_out=(/ E1, k1*st, 0., k1*ct /)
    q=k_in-k_out
    W=p+q
    vecW = sqrt(dot_product(W(1:3),W(1:3)))

    !determine Epi
    ! 1)  unit vector in the direction W
    unit_W=W(1:3)/vecW
    ! 2) unit_vector in the direction ppi_out
    stPi=sqrt((1.-ctPi)*(1.+ctPi))
    cphiPi=cos(phiPi)
    sphiPi=sin(phiPi)
    unit_ppi_out=(/ stPi*cphiPi, stPi*sphiPi, ctPi/)
    ! 3) cos between W and ppi
    ct_Wpi=dot_product(unit_W(1:3),unit_ppi_out(1:3))
    if (debug) write(*,'((A,3g12.5))') 'cos between W and ppi     =', ct_Wpi
    ! 4) Epi
    Epi=PionEnergy_W_ctWpi(W,ct_Wpi,mN1,mpi,numRoots,Epir(1),Epir(2))
    ! end determine Epi

    if (debug) then
       write (*,*)  ' '
       write (*,'(A,I5)')  ' Olga: number of roots =', numRoots
       write (*,*)  ' '
       write (*,'(A,4(g12.5))')  'W    = ',W
       write (*,'(A,4(g15.5),A,g12.5)')  'ppi  = ',ppi_out,'     ppi_out^2=',SP(ppi_out,ppi_out)
       write (*,'(A,4(g15.5),A,g12.5)')  'p1   = ',p_out,'     p1^2 = ',SP(p_out,p_out)
    end if


    ! nullifying cross sections before summing up over the two (possibly two) pion energies
    xsec_cDp=0.
    xsec_Np=0.
    xsec_cNp=0.
    xsec_CT=0.
    xsec_pp=0.
    xsec_pF=0.

    do r=1,numRoots

       ppi=sqrt((Epir(r)-mpi)*(Epir(r)+mpi))
       ppi_out=(/ Epir(r), ppi*stPi*cphiPi, ppi*stPi*sphiPi, ppi*ctPi/)
       p_out=W-ppi_out


       if (Epir(r)<mpi) then
          xsec_cDp=0.
          xsec_Np=0.
          xsec_cNp=0.
          xsec_CT=0.
          xsec_pp=0.
          xsec_pF=0.
          cycle
       else


          ! scalar products
          Q2=-abs4Sq(q)
          pk=mN*Enu
          p1k=SP(p_out,k_in)
          pq=mN*(Enu-E1)
          p1q=SP(p_out,q)
          p1p=SP(p,p_out)


          !  form factors for different diagrams
          if (process_ID.eq.2) then  !CC
             mrho2=meson(rho)%mass**2
             t=abs4Sq(q-ppi_out)     ! the momentum of the rho-meson squared
             F_rho=1./(1.-t/mrho2)
          else if (process_ID.eq.1)then  !EM
             F_rho=0
          else
             stop 'Wrong process number. It should be either 1 (EM)  or 2 (CC)'
          end if

          call formfactors_QE(Q2, CC, 1, ff_nucleon(1), ff_nucleon(2), ff_nucleon(3), ff_nucleon(4)) ! CC is process-ID,  1 is charge-in
          FV_CT=ff_nucleon(1)
          F_PF=ff_nucleon(1)

          ! should be corrected, because for  crossed nucleon pole the "charge-in" is the charge of the OUTGOING nucleon,
          ! for CC we have p target and p pi+ final state, which is OK
          ! for EM we have p target , but outgoing can be both p and n
          call formfactors_QE(Q2, process_ID, 1, ff_nucleon(1), ff_nucleon(2), ff_nucleon(3), ff_nucleon(4)) !   1 is charge-in
          f1V=ff_nucleon(1)
          f1A=-ff_nucleon(3)         ! in FF_QE_nucleonScattering.f90   f1A is defined as negative
          f2V=ff_nucleon(2)/2./mN

          MR=baryon(Delta)%mass
          ff_Delta=getFormfactor_Res(Q2,MR,Delta,1,CC,success_ff) ! 1=charge-in=target_charge
          C5A= ff_Delta(7)

          ! crossed Delta pole, contribution from C5A only

          ! the first multiplier is Delta-N-pi coupling
          M2cDp=  (2.14/0.14)**2 *(16*C5A**2*(4*ml**4*mN**8 - 4*ml**2*mN**10 - 4*ml**2*mN**9*MR - 2*ml**4*mN**6*MR**2 +   &
               &    6*ml**2*mN**8*MR**2 + 12*ml**2*mN**7*MR**3 - 7*ml**4*mN**4*MR**4 -                           &
               &    2*ml**4*mN**3*MR**5 - 16*ml**2*mN**5*MR**5 + 3*ml**4*mN**2*MR**6 -                           &
               &    10*ml**2*mN**4*MR**6 + 16*ml**2*mN**8*p1k + 8*ml**2*mN**7*MR*p1k -                           &
               &    20*ml**2*mN**6*MR**2*p1k - 12*ml**2*mN**5*MR**3*p1k - 4*ml**2*mN**4*MR**4*p1k -              &
               &    4*ml**2*mN**3*MR**5*p1k + 16*mN**8*p1k**2 + 16*mN**7*MR*p1k**2 -                             &
               &    32*mN**6*MR**2*p1k**2 - 24*mN**5*MR**3*p1k**2 + 20*mN**4*MR**4*p1k**2 -                      &
               &    12*mN**2*MR**6*p1k**2 - 4*ml**4*mN**6*p1p + 4*ml**2*mN**8*p1p +                              &
               &    8*ml**4*mN**5*MR*p1p + 4*ml**2*mN**7*MR*p1p - 8*ml**4*mN**4*MR**2*p1p -                      &
               &    4*ml**2*mN**6*MR**2*p1p - 24*ml**4*mN**3*MR**3*p1p - 12*ml**2*mN**5*MR**3*p1p +              &
               &    7*ml**4*mN**2*MR**4*p1p - 12*ml**2*mN**4*MR**4*p1p + 14*ml**4*mN*MR**5*p1p +                 &
               &    3*ml**4*MR**6*p1p + 4*ml**2*mN**2*MR**6*p1p - 16*ml**2*mN**6*p1k*p1p +                       &
               &    16*ml**2*mN**5*MR*p1k*p1p - 8*ml**2*mN**4*MR**2*p1k*p1p -                                    &
               &    76*ml**2*mN**3*MR**3*p1k*p1p - 8*ml**2*mN**2*MR**4*p1k*p1p +                                 &
               &    52*ml**2*mN*MR**5*p1k*p1p + 24*ml**2*MR**6*p1k*p1p - 16*mN**6*p1k**2*p1p +                   &
               &    16*mN**4*MR**2*p1k**2*p1p - 56*mN**3*MR**3*p1k**2*p1p -                                      &
               &    44*mN**2*MR**4*p1k**2*p1p + 48*mN*MR**5*p1k**2*p1p + 36*MR**6*p1k**2*p1p -                   &
               &    4*ml**4*mN**4*p1p**2 + 4*ml**2*mN**6*p1p**2 - 16*ml**4*mN**3*MR*p1p**2 +                     &
               &    4*ml**2*mN**5*MR*p1p**2 + 6*ml**4*mN**2*MR**2*p1p**2 -                                       &
               &    10*ml**2*mN**4*MR**2*p1p**2 + 24*ml**4*mN*MR**3*p1p**2 -                                     &
               &    12*ml**2*mN**3*MR**3*p1p**2 + 6*ml**4*MR**4*p1p**2 + 8*ml**2*mN**2*MR**4*p1p**2 +            &
               &    16*ml**2*mN*MR**5*p1p**2 + 6*ml**2*MR**6*p1p**2 - 16*ml**2*mN**4*p1k*p1p**2 -                &
               &    56*ml**2*mN**3*MR*p1k*p1p**2 + 12*ml**2*mN**2*MR**2*p1k*p1p**2 +                             &
               &    88*ml**2*mN*MR**3*p1k*p1p**2 + 36*ml**2*MR**4*p1k*p1p**2 -                                   &
               &    16*mN**4*p1k**2*p1p**2 - 48*mN**3*MR*p1k**2*p1p**2 + 80*mN*MR**3*p1k**2*p1p**2 +             &
               &    48*MR**4*p1k**2*p1p**2 + 4*ml**4*mN**2*p1p**3 - 4*ml**2*mN**4*p1p**3 +                       &
               &    8*ml**4*mN*MR*p1p**3 - 4*ml**2*mN**3*MR*p1p**3 + 4*ml**4*MR**2*p1p**3 +                      &
               &    8*ml**2*mN**2*MR**2*p1p**3 + 12*ml**2*mN*MR**3*p1p**3 + 4*ml**2*MR**4*p1p**3 +               &
               &    16*ml**2*mN**2*p1k*p1p**3 + 32*ml**2*mN*MR*p1k*p1p**3 +                                      &
               &    16*ml**2*MR**2*p1k*p1p**3 + 16*mN**2*p1k**2*p1p**3 + 32*mN*MR*p1k**2*p1p**3 +                &
               &    16*MR**2*p1k**2*p1p**3 - 24*ml**4*mN**6*p1q + 24*ml**2*mN**8*p1q -                           &
               &    4*ml**4*mN**5*MR*p1q + 28*ml**2*mN**7*MR*p1q + 12*ml**4*mN**4*MR**2*p1q -                    &
               &    24*ml**2*mN**6*MR**2*p1q + 6*ml**4*mN**3*MR**3*p1q - 58*ml**2*mN**5*MR**3*p1q +              &
               &    6*ml**4*mN**2*MR**4*p1q - 6*ml**2*mN**4*MR**4*p1q - 4*ml**4*mN*MR**5*p1q +                   &
               &    34*ml**2*mN**3*MR**5*p1q + 10*ml**2*mN**2*MR**6*p1q - 96*ml**2*mN**6*p1k*p1q -               &
               &    16*mN**8*p1k*p1q - 64*ml**2*mN**5*MR*p1k*p1q - 16*mN**7*MR*p1k*p1q +                         &
               &    96*ml**2*mN**4*MR**2*p1k*p1q + 32*mN**6*MR**2*p1k*p1q +                                      &
               &    80*ml**2*mN**3*MR**3*p1k*p1q + 24*mN**5*MR**3*p1k*p1q -                                      &
               &    24*ml**2*mN**2*MR**4*p1k*p1q - 20*mN**4*MR**4*p1k*p1q -                                      &
               &    32*ml**2*mN*MR**5*p1k*p1q + 12*mN**2*MR**6*p1k*p1q - 96*mN**6*p1k**2*p1q -                   &
               &    112*mN**5*MR*p1k**2*p1q + 144*mN**4*MR**2*p1k**2*p1q +                                       &
               &    136*mN**3*MR**3*p1k**2*p1q - 72*mN**2*MR**4*p1k**2*p1q - 48*mN*MR**5*p1k**2*p1q +            &
               &    16*ml**4*mN**4*p1p*p1q - 16*ml**2*mN**6*p1p*p1q - 36*ml**4*mN**3*MR*p1p*p1q -                &
               &    20*ml**2*mN**5*MR*p1p*p1q + 8*ml**4*mN**2*MR**2*p1p*p1q +                                    &
               &    22*ml**2*mN**4*MR**2*p1p*p1q + 48*ml**4*mN*MR**3*p1p*p1q +                                   &
               &    58*ml**2*mN**3*MR**3*p1p*p1q + 12*ml**4*MR**4*p1p*p1q +                                      &
               &    12*ml**2*mN**2*MR**4*p1p*p1q - 26*ml**2*mN*MR**5*p1p*p1q -                                   &
               &    6*ml**2*MR**6*p1p*p1q + 64*ml**2*mN**4*p1k*p1p*p1q + 16*mN**6*p1k*p1p*p1q -                  &
               &    80*ml**2*mN**3*MR*p1k*p1p*p1q - 40*ml**2*mN**2*MR**2*p1k*p1p*p1q -                           &
               &    16*mN**4*MR**2*p1k*p1p*p1q + 152*ml**2*mN*MR**3*p1k*p1p*p1q +                                &
               &    56*mN**3*MR**3*p1k*p1p*p1q + 96*ml**2*MR**4*p1k*p1p*p1q +                                    &
               &    44*mN**2*MR**4*p1k*p1p*p1q - 48*mN*MR**5*p1k*p1p*p1q - 36*MR**6*p1k*p1p*p1q +                &
               &    64*mN**4*p1k**2*p1p*p1q - 16*mN**3*MR*p1k**2*p1p*p1q -                                       &
               &    112*mN**2*MR**2*p1k**2*p1p*p1q + 112*mN*MR**3*p1k**2*p1p*p1q +                               &
               &    144*MR**4*p1k**2*p1p*p1q + 8*ml**4*mN**2*p1p**2*p1q - 8*ml**2*mN**4*p1p**2*p1q +             &
               &    40*ml**4*mN*MR*p1p**2*p1q + 12*ml**4*MR**2*p1p**2*p1q +                                      &
               &    10*ml**2*mN**2*MR**2*p1p**2*p1q - 2*ml**2*MR**4*p1p**2*p1q +                                 &
               &    32*ml**2*mN**2*p1k*p1p**2*p1q + 16*mN**4*p1k*p1p**2*p1q +                                    &
               &    144*ml**2*mN*MR*p1k*p1p**2*p1q + 48*mN**3*MR*p1k*p1p**2*p1q +                                &
               &    72*ml**2*MR**2*p1k*p1p**2*p1q - 80*mN*MR**3*p1k*p1p**2*p1q -                                 &
               &    48*MR**4*p1k*p1p**2*p1q + 32*mN**2*p1k**2*p1p**2*p1q +                                       &
               &    128*mN*MR*p1k**2*p1p**2*p1q + 96*MR**2*p1k**2*p1p**2*p1q -                                   &
               &    8*ml**2*mN*MR*p1p**3*p1q - 8*ml**2*MR**2*p1p**3*p1q - 16*mN**2*p1k*p1p**3*p1q -              &
               &    32*mN*MR*p1k*p1p**3*p1q - 16*MR**2*p1k*p1p**3*p1q + 48*ml**4*mN**4*p1q**2 -                  &
               &    48*ml**2*mN**6*p1q**2 + 16*ml**4*mN**3*MR*p1q**2 - 64*ml**2*mN**5*MR*p1q**2 -                &
               &    16*ml**4*mN**2*MR**2*p1q**2 + 16*ml**2*mN**4*MR**2*p1q**2 -                                  &
               &    16*ml**4*mN*MR**3*p1q**2 + 72*ml**2*mN**3*MR**3*p1q**2 +                                     &
               &    28*ml**2*mN**2*MR**4*p1q**2 + 192*ml**2*mN**4*p1k*p1q**2 + 96*mN**6*p1k*p1q**2 +             &
               &    160*ml**2*mN**3*MR*p1k*p1q**2 + 112*mN**5*MR*p1k*p1q**2 -                                    &
               &    112*ml**2*mN**2*MR**2*p1k*p1q**2 - 144*mN**4*MR**2*p1k*p1q**2 -                              &
               &    128*ml**2*mN*MR**3*p1k*p1q**2 - 136*mN**3*MR**3*p1k*p1q**2 +                                 &
               &    72*mN**2*MR**4*p1k*p1q**2 + 48*mN*MR**5*p1k*p1q**2 + 192*mN**4*p1k**2*p1q**2 +               &
               &    256*mN**3*MR*p1k**2*p1q**2 - 160*mN**2*MR**2*p1k**2*p1q**2 -                                 &
               &    192*mN*MR**3*p1k**2*p1q**2 - 16*ml**4*mN**2*p1p*p1q**2 +                                     &
               &    16*ml**2*mN**4*p1p*p1q**2 + 40*ml**4*mN*MR*p1p*p1q**2 +                                      &
               &    40*ml**2*mN**3*MR*p1p*p1q**2 + 12*ml**4*MR**2*p1p*p1q**2 -                                   &
               &    12*ml**2*mN**2*MR**2*p1p*p1q**2 - 68*ml**2*mN*MR**3*p1p*p1q**2 -                             &
               &    24*ml**2*MR**4*p1p*p1q**2 - 64*ml**2*mN**2*p1k*p1p*p1q**2 -                                  &
               &    64*mN**4*p1k*p1p*p1q**2 + 96*ml**2*mN*MR*p1k*p1p*p1q**2 +                                    &
               &    16*mN**3*MR*p1k*p1p*p1q**2 + 96*ml**2*MR**2*p1k*p1p*p1q**2 +                                 &
               &    112*mN**2*MR**2*p1k*p1p*p1q**2 - 112*mN*MR**3*p1k*p1p*p1q**2 -                               &
               &    144*MR**4*p1k*p1p*p1q**2 - 64*mN**2*p1k**2*p1p*p1q**2 +                                      &
               &    32*mN*MR*p1k**2*p1p*p1q**2 + 144*MR**2*p1k**2*p1p*p1q**2 -                                   &
               &    32*ml**2*mN*MR*p1p**2*p1q**2 - 28*ml**2*MR**2*p1p**2*p1q**2 -                                &
               &    32*mN**2*p1k*p1p**2*p1q**2 - 128*mN*MR*p1k*p1p**2*p1q**2 -                                   &
               &    96*MR**2*p1k*p1p**2*p1q**2 - 32*ml**4*mN**2*p1q**3 + 32*ml**2*mN**4*p1q**3 -                 &
               &    16*ml**4*mN*MR*p1q**3 + 48*ml**2*mN**3*MR*p1q**3 + 16*ml**2*mN**2*MR**2*p1q**3 -             &
               &    128*ml**2*mN**2*p1k*p1q**3 - 192*mN**4*p1k*p1q**3 - 128*ml**2*mN*MR*p1k*p1q**3 -             &
               &    256*mN**3*MR*p1k*p1q**3 + 160*mN**2*MR**2*p1k*p1q**3 + 192*mN*MR**3*p1k*p1q**3 -             &
               &    128*mN**2*p1k**2*p1q**3 - 192*mN*MR*p1k**2*p1q**3 - 32*ml**2*mN*MR*p1p*p1q**3 -              &
               &    24*ml**2*MR**2*p1p*p1q**3 + 64*mN**2*p1k*p1p*p1q**3 - 32*mN*MR*p1k*p1p*p1q**3 -              &
               &    144*MR**2*p1k*p1p*p1q**3 + 128*mN**2*p1k*p1q**4 + 192*mN*MR*p1k*p1q**4 -                     &
               &    8*ml**2*mN**7*MR*pk + 16*ml**2*mN**6*MR**2*pk + 28*ml**2*mN**5*MR**3*pk -                    &
               &    20*ml**2*mN**3*MR**5*pk - 16*ml**2*mN**2*MR**6*pk - 16*mN**7*MR*p1k*pk +                     &
               &    24*mN**6*MR**2*p1k*pk + 80*mN**5*MR**3*p1k*pk - 64*mN**3*MR**5*p1k*pk -                      &
               &    24*mN**2*MR**6*p1k*pk + 16*ml**2*mN**5*MR*p1p*pk - 12*ml**2*mN**3*MR**3*p1p*pk -             &
               &    16*ml**2*mN**2*MR**4*p1p*pk - 28*ml**2*mN*MR**5*p1p*pk - 8*ml**2*MR**6*p1p*pk +              &
               &    32*mN**5*MR*p1k*p1p*pk + 16*mN**4*MR**2*p1k*p1p*pk - 48*mN**3*MR**3*p1k*p1p*pk -             &
               &    40*mN**2*MR**4*p1k*p1p*pk - 32*mN*MR**5*p1k*p1p*pk - 24*MR**6*p1k*p1p*pk -                   &
               &    8*ml**2*mN**3*MR*p1p**2*pk - 16*ml**2*mN**2*MR**2*p1p**2*pk -                                &
               &    16*ml**2*mN*MR**3*p1p**2*pk - 8*ml**2*MR**4*p1p**2*pk -                                      &
               &    16*mN**3*MR*p1k*p1p**2*pk - 40*mN**2*MR**2*p1k*p1p**2*pk -                                   &
               &    32*mN*MR**3*p1k*p1p**2*pk - 8*MR**4*p1k*p1p**2*pk + 48*ml**2*mN**5*MR*p1q*pk +               &
               &    8*mN**7*MR*p1q*pk - 60*ml**2*mN**4*MR**2*p1q*pk - 12*mN**6*MR**2*p1q*pk -                    &
               &    108*ml**2*mN**3*MR**3*p1q*pk - 40*mN**5*MR**3*p1q*pk -                                       &
               &    8*ml**2*mN**2*MR**4*p1q*pk + 36*ml**2*mN*MR**5*p1q*pk + 32*mN**3*MR**5*p1q*pk +              &
               &    4*ml**2*MR**6*p1q*pk + 12*mN**2*MR**6*p1q*pk + 96*mN**5*MR*p1k*p1q*pk -                      &
               &    72*mN**4*MR**2*p1k*p1q*pk - 312*mN**3*MR**3*p1k*p1q*pk -                                     &
               &    48*mN**2*MR**4*p1k*p1q*pk + 120*mN*MR**5*p1k*p1q*pk + 24*MR**6*p1k*p1q*pk -                  &
               &    64*ml**2*mN**3*MR*p1p*p1q*pk - 16*mN**5*MR*p1p*p1q*pk -                                      &
               &    12*ml**2*mN**2*MR**2*p1p*p1q*pk - 8*mN**4*MR**2*p1p*p1q*pk +                                 &
               &    24*mN**3*MR**3*p1p*p1q*pk - 4*ml**2*MR**4*p1p*p1q*pk +                                       &
               &    20*mN**2*MR**4*p1p*p1q*pk + 16*mN*MR**5*p1p*p1q*pk + 12*MR**6*p1p*p1q*pk -                   &
               &    128*mN**3*MR*p1k*p1p*p1q*pk - 88*mN**2*MR**2*p1k*p1p*p1q*pk +                                &
               &    48*mN*MR**3*p1k*p1p*p1q*pk + 8*MR**4*p1k*p1p*p1q*pk +                                        &
               &    16*ml**2*mN*MR*p1p**2*p1q*pk + 8*mN**3*MR*p1p**2*p1q*pk +                                    &
               &    8*ml**2*MR**2*p1p**2*p1q*pk + 20*mN**2*MR**2*p1p**2*p1q*pk +                                 &
               &    16*mN*MR**3*p1p**2*p1q*pk + 4*MR**4*p1p**2*p1q*pk + 32*mN*MR*p1k*p1p**2*p1q*pk +             &
               &    32*MR**2*p1k*p1p**2*p1q*pk - 96*ml**2*mN**3*MR*p1q**2*pk -                                   &
               &    48*mN**5*MR*p1q**2*pk + 48*ml**2*mN**2*MR**2*p1q**2*pk +                                     &
               &    36*mN**4*MR**2*p1q**2*pk + 104*ml**2*mN*MR**3*p1q**2*pk +                                    &
               &    156*mN**3*MR**3*p1q**2*pk + 16*ml**2*MR**4*p1q**2*pk + 24*mN**2*MR**4*p1q**2*pk -            &
               &    60*mN*MR**5*p1q**2*pk - 12*MR**6*p1q**2*pk - 192*mN**3*MR*p1k*p1q**2*pk +                    &
               &    304*mN*MR**3*p1k*p1q**2*pk + 96*MR**4*p1k*p1q**2*pk +                                        &
               &    64*ml**2*mN*MR*p1p*p1q**2*pk + 64*mN**3*MR*p1p*p1q**2*pk +                                   &
               &    24*ml**2*MR**2*p1p*p1q**2*pk + 44*mN**2*MR**2*p1p*p1q**2*pk -                                &
               &    24*mN*MR**3*p1p*p1q**2*pk - 4*MR**4*p1p*p1q**2*pk + 128*mN*MR*p1k*p1p*p1q**2*pk +            &
               &    112*MR**2*p1k*p1p*p1q**2*pk - 16*mN*MR*p1p**2*p1q**2*pk -                                    &
               &    16*MR**2*p1p**2*p1q**2*pk + 64*ml**2*mN*MR*p1q**3*pk + 96*mN**3*MR*p1q**3*pk +               &
               &    16*ml**2*MR**2*p1q**3*pk - 152*mN*MR**3*p1q**3*pk - 48*MR**4*p1q**3*pk +                     &
               &    128*mN*MR*p1k*p1q**3*pk + 96*MR**2*p1k*p1q**3*pk - 64*mN*MR*p1p*p1q**3*pk -                  &
               &    56*MR**2*p1p*p1q**3*pk - 64*mN*MR*p1q**4*pk - 48*MR**2*p1q**4*pk -                           &
               &    8*mN**5*MR**3*pk**2 + 4*mN**4*MR**4*pk**2 + 32*mN**3*MR**5*pk**2 +                           &
               &    20*mN**2*MR**6*pk**2 + 8*mN**3*MR**3*p1p*pk**2 + 20*mN**2*MR**4*p1p*pk**2 +                  &
               &    16*mN*MR**5*p1p*pk**2 + 4*MR**6*p1p*pk**2 + 32*mN**3*MR**3*p1q*pk**2 -                       &
               &    40*mN*MR**5*p1q*pk**2 - 8*MR**6*p1q*pk**2 - 16*mN*MR**3*p1p*p1q*pk**2 -                      &
               &    16*MR**4*p1p*p1q*pk**2 - 32*mN*MR**3*p1q**2*pk**2 - 16*MR**4*p1q**2*pk**2 -                  &
               &    4*ml**4*mN**5*MR*pq - 4*ml**2*mN**7*MR*pq + 18*ml**4*mN**4*MR**2*pq -                        &
               &    14*ml**2*mN**6*MR**2*pq + 18*ml**4*mN**3*MR**3*pq + 2*ml**2*mN**5*MR**3*pq -                 &
               &    16*ml**4*mN**2*MR**4*pq + 24*ml**2*mN**4*MR**4*pq - 12*ml**4*mN*MR**5*pq +                   &
               &    10*ml**2*mN**3*MR**5*pq - 2*ml**2*mN**2*MR**6*pq + 8*mN**7*MR*p1k*pq +                       &
               &    48*ml**2*mN**4*MR**2*p1k*pq - 12*mN**6*MR**2*p1k*pq +                                        &
               &    52*ml**2*mN**3*MR**3*p1k*pq - 40*mN**5*MR**3*p1k*pq -                                        &
               &    28*ml**2*mN**2*MR**4*p1k*pq - 44*ml**2*mN*MR**5*p1k*pq + 32*mN**3*MR**5*p1k*pq -             &
               &    12*ml**2*MR**6*p1k*pq + 12*mN**2*MR**6*p1k*pq + 16*mN**5*MR*p1k**2*pq +                      &
               &    24*mN**4*MR**2*p1k**2*pq + 32*mN**3*MR**3*p1k**2*pq + 8*mN**2*MR**4*p1k**2*pq -              &
               &    40*mN*MR**5*p1k**2*pq - 24*MR**6*p1k**2*pq + 16*ml**4*mN**4*p1p*pq -                         &
               &    16*ml**2*mN**6*p1p*pq + 28*ml**4*mN**3*MR*p1p*pq - 12*ml**2*mN**5*MR*p1p*pq -                &
               &    26*ml**4*mN**2*MR**2*p1p*pq + 38*ml**2*mN**4*MR**2*p1p*pq -                                  &
               &    44*ml**4*mN*MR**3*p1p*pq + 26*ml**2*mN**3*MR**3*p1p*pq - 6*ml**4*MR**4*p1p*pq -              &
               &    24*ml**2*mN**2*MR**4*p1p*pq - 18*ml**2*mN*MR**5*p1p*pq - 2*ml**2*MR**6*p1p*pq +              &
               &    64*ml**2*mN**4*p1k*p1p*pq + 96*ml**2*mN**3*MR*p1k*p1p*pq -                                   &
               &    16*mN**5*MR*p1k*p1p*pq - 80*ml**2*mN**2*MR**2*p1k*p1p*pq -                                   &
               &    8*mN**4*MR**2*p1k*p1p*pq - 160*ml**2*mN*MR**3*p1k*p1p*pq +                                   &
               &    24*mN**3*MR**3*p1k*p1p*pq - 48*ml**2*MR**4*p1k*p1p*pq +                                      &
               &    20*mN**2*MR**4*p1k*p1p*pq + 16*mN*MR**5*p1k*p1p*pq + 12*MR**6*p1k*p1p*pq +                   &
               &    64*mN**4*p1k**2*p1p*pq + 80*mN**3*MR*p1k**2*p1p*pq -                                         &
               &    56*mN**2*MR**2*p1k**2*p1p*pq - 144*mN*MR**3*p1k**2*p1p*pq -                                  &
               &    72*MR**4*p1k**2*p1p*pq - 16*ml**4*mN**2*p1p**2*pq + 16*ml**2*mN**4*p1p**2*pq -               &
               &    24*ml**4*mN*MR*p1p**2*pq + 16*ml**2*mN**3*MR*p1p**2*pq -                                     &
               &    8*ml**4*MR**2*p1p**2*pq - 24*ml**2*mN**2*MR**2*p1p**2*pq -                                   &
               &    28*ml**2*mN*MR**3*p1p**2*pq - 4*ml**2*MR**4*p1p**2*pq -                                      &
               &    64*ml**2*mN**2*p1k*p1p**2*pq - 96*ml**2*mN*MR*p1k*p1p**2*pq +                                &
               &    8*mN**3*MR*p1k*p1p**2*pq - 32*ml**2*MR**2*p1k*p1p**2*pq +                                    &
               &    20*mN**2*MR**2*p1k*p1p**2*pq + 16*mN*MR**3*p1k*p1p**2*pq +                                   &
               &    4*MR**4*p1k*p1p**2*pq - 64*mN**2*p1k**2*p1p**2*pq - 96*mN*MR*p1k**2*p1p**2*pq -              &
               &    32*MR**2*p1k**2*p1p**2*pq + 8*ml**4*mN**4*p1q*pq - 8*ml**2*mN**6*p1q*pq +                    &
               &    20*ml**4*mN**3*MR*p1q*pq + 12*ml**2*mN**5*MR*p1q*pq -                                        &
               &    50*ml**4*mN**2*MR**2*p1q*pq + 50*ml**2*mN**4*MR**2*p1q*pq -                                  &
               &    36*ml**4*mN*MR**3*p1q*pq - 8*ml**2*mN**3*MR**3*p1q*pq + 8*ml**4*MR**4*p1q*pq -               &
               &    34*ml**2*mN**2*MR**4*p1q*pq + 4*ml**2*mN*MR**5*p1q*pq + 4*ml**2*MR**6*p1q*pq +               &
               &    32*ml**2*mN**4*p1k*p1q*pq + 16*ml**2*mN**3*MR*p1k*p1q*pq -                                   &
               &    64*mN**5*MR*p1k*p1q*pq - 128*ml**2*mN**2*MR**2*p1k*p1q*pq +                                  &
               &    12*mN**4*MR**2*p1k*p1q*pq - 104*ml**2*mN*MR**3*p1k*p1q*pq +                                  &
               &    124*mN**3*MR**3*p1k*p1q*pq - 16*ml**2*MR**4*p1k*p1q*pq +                                     &
               &    16*mN**2*MR**4*p1k*p1q*pq - 20*mN*MR**5*p1k*p1q*pq + 12*MR**6*p1k*p1q*pq +                   &
               &    32*mN**4*p1k**2*p1q*pq - 48*mN**3*MR*p1k**2*p1q*pq -                                         &
               &    56*mN**2*MR**2*p1k**2*p1q*pq - 64*mN*MR**3*p1k**2*p1q*pq -                                   &
               &    64*MR**4*p1k**2*p1q*pq - 48*ml**4*mN**2*p1p*p1q*pq + 48*ml**2*mN**4*p1p*p1q*pq -             &
               &    72*ml**4*mN*MR*p1p*p1q*pq + 24*ml**2*mN**3*MR*p1p*p1q*pq +                                   &
               &    4*ml**4*MR**2*p1p*p1q*pq - 58*ml**2*mN**2*MR**2*p1p*p1q*pq +                                 &
               &    26*ml**2*MR**4*p1p*p1q*pq - 192*ml**2*mN**2*p1k*p1p*p1q*pq -                                 &
               &    64*mN**4*p1k*p1p*p1q*pq - 256*ml**2*mN*MR*p1k*p1p*p1q*pq -                                   &
               &    16*mN**3*MR*p1k*p1p*p1q*pq - 32*ml**2*MR**2*p1k*p1p*p1q*pq +                                 &
               &    100*mN**2*MR**2*p1k*p1p*p1q*pq + 120*mN*MR**3*p1k*p1p*p1q*pq +                               &
               &    68*MR**4*p1k*p1p*p1q*pq - 192*mN**2*p1k**2*p1p*p1q*pq -                                      &
               &    224*mN*MR*p1k**2*p1p*p1q*pq - 80*MR**2*p1k**2*p1p*p1q*pq +                                   &
               &    8*ml**4*p1p**2*p1q*pq - 8*ml**2*mN**2*p1p**2*p1q*pq +                                        &
               &    16*ml**2*mN*MR*p1p**2*p1q*pq + 28*ml**2*MR**2*p1p**2*p1q*pq +                                &
               &    32*ml**2*p1k*p1p**2*p1q*pq + 64*mN**2*p1k*p1p**2*p1q*pq +                                    &
               &    80*mN*MR*p1k*p1p**2*p1q*pq + 16*MR**2*p1k*p1p**2*p1q*pq +                                    &
               &    32*p1k**2*p1p**2*p1q*pq - 32*ml**4*mN**2*p1q**2*pq + 32*ml**2*mN**4*p1q**2*pq -              &
               &    24*ml**4*mN*MR*p1q**2*pq - 8*ml**2*mN**3*MR*p1q**2*pq +                                      &
               &    32*ml**4*MR**2*p1q**2*pq - 56*ml**2*mN**2*MR**2*p1q**2*pq +                                  &
               &    8*ml**2*mN*MR**3*p1q**2*pq + 16*ml**2*MR**4*p1q**2*pq -                                      &
               &    128*ml**2*mN**2*p1k*p1q**2*pq - 32*mN**4*p1k*p1q**2*pq -                                     &
               &    32*ml**2*mN*MR*p1k*p1q**2*pq + 144*mN**3*MR*p1k*p1q**2*pq +                                  &
               &    80*ml**2*MR**2*p1k*p1q**2*pq + 56*mN**2*MR**2*p1k*p1q**2*pq -                                &
               &    88*mN*MR**3*p1k*p1q**2*pq + 16*MR**4*p1k*p1q**2*pq - 128*mN**2*p1k**2*p1q**2*pq +            &
               &    32*mN*MR*p1k**2*p1q**2*pq + 32*MR**2*p1k**2*p1q**2*pq + 32*ml**4*p1p*p1q**2*pq -             &
               &    32*ml**2*mN**2*p1p*p1q**2*pq + 32*ml**2*mN*MR*p1p*p1q**2*pq +                                &
               &    60*ml**2*MR**2*p1p*p1q**2*pq + 128*ml**2*p1k*p1p*p1q**2*pq +                                 &
               &    192*mN**2*p1k*p1p*p1q**2*pq + 160*mN*MR*p1k*p1p*p1q**2*pq +                                  &
               &    24*MR**2*p1k*p1p*p1q**2*pq + 128*p1k**2*p1p*p1q**2*pq - 32*p1k*p1p**2*p1q**2*pq +            &
               &    32*ml**4*p1q**3*pq - 32*ml**2*mN**2*p1q**3*pq + 16*ml**2*MR**2*p1q**3*pq +                   &
               &    128*ml**2*p1k*p1q**3*pq + 128*mN**2*p1k*p1q**3*pq - 96*mN*MR*p1k*p1q**3*pq -                 &
               &    80*MR**2*p1k*p1q**3*pq + 128*p1k**2*p1q**3*pq - 128*p1k*p1p*p1q**3*pq -                      &
               &    128*p1k*p1q**4*pq - 16*ml**2*mN**5*MR*pk*pq - 8*ml**2*mN**4*MR**2*pk*pq +                    &
               &    16*ml**2*mN**3*MR**3*pk*pq + 8*mN**5*MR**3*pk*pq + 32*ml**2*mN**2*MR**4*pk*pq -              &
               &    4*mN**4*MR**4*pk*pq + 24*ml**2*mN*MR**5*pk*pq - 32*mN**3*MR**5*pk*pq -                       &
               &    20*mN**2*MR**6*pk*pq - 32*mN**5*MR*p1k*pk*pq - 32*mN**4*MR**2*p1k*pk*pq +                    &
               &    56*mN**3*MR**3*p1k*pk*pq + 72*mN**2*MR**4*p1k*pk*pq + 24*mN*MR**5*p1k*pk*pq +                &
               &    8*MR**6*p1k*pk*pq + 16*ml**2*mN**3*MR*p1p*pk*pq +                                            &
               &    40*ml**2*mN**2*MR**2*p1p*pk*pq + 32*ml**2*mN*MR**3*p1p*pk*pq -                               &
               &    8*mN**3*MR**3*p1p*pk*pq + 8*ml**2*MR**4*p1p*pk*pq - 20*mN**2*MR**4*p1p*pk*pq -               &
               &    16*mN*MR**5*p1p*pk*pq - 4*MR**6*p1p*pk*pq + 32*mN**3*MR*p1k*p1p*pk*pq +                      &
               &    96*mN**2*MR**2*p1k*p1p*pk*pq + 64*mN*MR**3*p1k*p1p*pk*pq +                                   &
               &    64*ml**2*mN**3*MR*p1q*pk*pq + 16*mN**5*MR*p1q*pk*pq +                                        &
               &    44*ml**2*mN**2*MR**2*p1q*pk*pq + 16*mN**4*MR**2*p1q*pk*pq -                                  &
               &    8*ml**2*mN*MR**3*p1q*pk*pq - 60*mN**3*MR**3*p1q*pk*pq -                                      &
               &    28*ml**2*MR**4*p1q*pk*pq - 36*mN**2*MR**4*p1q*pk*pq + 28*mN*MR**5*p1q*pk*pq +                &
               &    4*MR**6*p1q*pk*pq + 128*mN**3*MR*p1k*p1q*pk*pq + 152*mN**2*MR**2*p1k*p1q*pk*pq -             &
               &    64*mN*MR**3*p1k*p1q*pk*pq - 72*MR**4*p1k*p1q*pk*pq -                                         &
               &    32*ml**2*mN*MR*p1p*p1q*pk*pq - 16*mN**3*MR*p1p*p1q*pk*pq -                                   &
               &    32*ml**2*MR**2*p1p*p1q*pk*pq - 48*mN**2*MR**2*p1p*p1q*pk*pq -                                &
               &    16*mN*MR**3*p1p*p1q*pk*pq + 16*MR**4*p1p*p1q*pk*pq - 64*mN*MR*p1k*p1p*p1q*pk*pq -            &
               &    96*MR**2*p1k*p1p*p1q*pk*pq - 64*ml**2*mN*MR*p1q**2*pk*pq -                                   &
               &    64*mN**3*MR*p1q**2*pk*pq - 56*ml**2*MR**2*p1q**2*pk*pq -                                     &
               &    76*mN**2*MR**2*p1q**2*pk*pq + 64*mN*MR**3*p1q**2*pk*pq + 52*MR**4*p1q**2*pk*pq -             &
               &    128*mN*MR*p1k*p1q**2*pk*pq - 176*MR**2*p1k*p1q**2*pk*pq +                                    &
               &    32*mN*MR*p1p*p1q**2*pk*pq + 48*MR**2*p1p*p1q**2*pk*pq + 64*mN*MR*p1q**3*pk*pq +              &
               &    88*MR**2*p1q**3*pk*pq - 8*mN**3*MR**3*pk**2*pq - 24*mN**2*MR**4*pk**2*pq -                   &
               &    16*mN*MR**5*pk**2*pq + 16*mN*MR**3*p1q*pk**2*pq + 24*MR**4*p1q*pk**2*pq -                    &
               &    12*ml**4*mN**4*pq**2 + 12*ml**2*mN**6*pq**2 - 12*ml**4*mN**3*MR*pq**2 +                      &
               &    8*ml**2*mN**5*MR*pq**2 + 20*ml**4*mN**2*MR**2*pq**2 -                                        &
               &    24*ml**2*mN**4*MR**2*pq**2 + 20*ml**4*mN*MR**3*pq**2 -                                       &
               &    16*ml**2*mN**3*MR**3*pq**2 + 8*ml**2*mN**2*MR**4*pq**2 + 4*ml**2*mN*MR**5*pq**2 -            &
               &    48*ml**2*mN**4*p1k*pq**2 - 40*ml**2*mN**3*MR*p1k*pq**2 + 16*mN**5*MR*p1k*pq**2 +             &
               &    68*ml**2*mN**2*MR**2*p1k*pq**2 + 16*mN**4*MR**2*p1k*pq**2 +                                  &
               &    72*ml**2*mN*MR**3*p1k*pq**2 - 28*mN**3*MR**3*p1k*pq**2 +                                     &
               &    12*ml**2*MR**4*p1k*pq**2 - 36*mN**2*MR**4*p1k*pq**2 - 12*mN*MR**5*p1k*pq**2 -                &
               &    4*MR**6*p1k*pq**2 - 48*mN**4*p1k**2*pq**2 - 32*mN**3*MR*p1k**2*pq**2 +                       &
               &    56*mN**2*MR**2*p1k**2*pq**2 + 64*mN*MR**3*p1k**2*pq**2 + 24*MR**4*p1k**2*pq**2 +             &
               &    20*ml**4*mN**2*p1p*pq**2 - 20*ml**2*mN**4*p1p*pq**2 + 24*ml**4*mN*MR*p1p*pq**2 -             &
               &    20*ml**2*mN**3*MR*p1p*pq**2 + 4*ml**4*MR**2*p1p*pq**2 +                                      &
               &    20*ml**2*mN**2*MR**2*p1p*pq**2 + 20*ml**2*mN*MR**3*p1p*pq**2 +                               &
               &    80*ml**2*mN**2*p1k*p1p*pq**2 + 96*ml**2*mN*MR*p1k*p1p*pq**2 -                                &
               &    16*mN**3*MR*p1k*p1p*pq**2 + 16*ml**2*MR**2*p1k*p1p*pq**2 -                                   &
               &    48*mN**2*MR**2*p1k*p1p*pq**2 - 32*mN*MR**3*p1k*p1p*pq**2 +                                   &
               &    80*mN**2*p1k**2*p1p*pq**2 + 96*mN*MR*p1k**2*p1p*pq**2 +                                      &
               &    16*MR**2*p1k**2*p1p*pq**2 + 40*ml**4*mN**2*p1q*pq**2 - 40*ml**2*mN**4*p1q*pq**2 +            &
               &    32*ml**4*mN*MR*p1q*pq**2 - 24*ml**2*mN**3*MR*p1q*pq**2 -                                     &
               &    16*ml**4*MR**2*p1q*pq**2 + 32*ml**2*mN**2*MR**2*p1q*pq**2 +                                  &
               &    4*ml**2*mN*MR**3*p1q*pq**2 - 8*ml**2*MR**4*p1q*pq**2 +                                       &
               &    160*ml**2*mN**2*p1k*p1q*pq**2 + 48*mN**4*p1k*p1q*pq**2 +                                     &
               &    112*ml**2*mN*MR*p1k*p1q*pq**2 - 32*mN**3*MR*p1k*p1q*pq**2 -                                  &
               &    40*ml**2*MR**2*p1k*p1q*pq**2 - 132*mN**2*MR**2*p1k*p1q*pq**2 -                               &
               &    32*mN*MR**3*p1k*p1q*pq**2 + 12*MR**4*p1k*p1q*pq**2 + 160*mN**2*p1k**2*p1q*pq**2 +            &
               &    96*mN*MR*p1k**2*p1q*pq**2 - 16*MR**2*p1k**2*p1q*pq**2 - 16*ml**4*p1p*p1q*pq**2 +             &
               &    16*ml**2*mN**2*p1p*p1q*pq**2 - 8*ml**2*mN*MR*p1p*p1q*pq**2 -                                 &
               &    24*ml**2*MR**2*p1p*p1q*pq**2 - 64*ml**2*p1k*p1p*p1q*pq**2 -                                  &
               &    80*mN**2*p1k*p1p*p1q*pq**2 - 64*mN*MR*p1k*p1p*p1q*pq**2 +                                    &
               &    32*MR**2*p1k*p1p*p1q*pq**2 - 64*p1k**2*p1p*p1q*pq**2 - 32*ml**4*p1q**2*pq**2 +               &
               &    32*ml**2*mN**2*p1q**2*pq**2 - 16*ml**2*MR**2*p1q**2*pq**2 -                                  &
               &    128*ml**2*p1k*p1q**2*pq**2 - 160*mN**2*p1k*p1q**2*pq**2 -                                    &
               &    32*mN*MR*p1k*p1q**2*pq**2 + 104*MR**2*p1k*p1q**2*pq**2 -                                     &
               &    128*p1k**2*p1q**2*pq**2 + 64*p1k*p1p*p1q**2*pq**2 + 128*p1k*p1q**3*pq**2 -                   &
               &    8*ml**2*mN**3*MR*pk*pq**2 - 24*ml**2*mN**2*MR**2*pk*pq**2 -                                  &
               &    16*ml**2*mN*MR**3*pk*pq**2 + 8*mN**3*MR**3*pk*pq**2 + 24*mN**2*MR**4*pk*pq**2 +              &
               &    16*mN*MR**5*pk*pq**2 - 16*mN**3*MR*p1k*pk*pq**2 - 56*mN**2*MR**2*p1k*pk*pq**2 -              &
               &    32*mN*MR**3*p1k*pk*pq**2 + 8*MR**4*p1k*pk*pq**2 + 16*ml**2*mN*MR*p1q*pk*pq**2 +              &
               &    8*mN**3*MR*p1q*pk*pq**2 + 24*ml**2*MR**2*p1q*pk*pq**2 +                                      &
               &    28*mN**2*MR**2*p1q*pk*pq**2 - 28*MR**4*p1q*pk*pq**2 + 32*mN*MR*p1k*p1q*pk*pq**2 +            &
               &    64*MR**2*p1k*p1q*pk*pq**2 - 16*mN*MR*p1q**2*pk*pq**2 - 32*MR**2*p1q**2*pk*pq**2 -            &
               &    8*ml**4*mN**2*pq**3 + 8*ml**2*mN**4*pq**3 - 8*ml**4*mN*MR*pq**3 +                            &
               &    8*ml**2*mN**3*MR*pq**3 - 4*ml**2*mN**2*MR**2*pq**3 - 4*ml**2*mN*MR**3*pq**3 -                &
               &    32*ml**2*mN**2*p1k*pq**3 - 32*ml**2*mN*MR*p1k*pq**3 + 8*mN**3*MR*p1k*pq**3 +                 &
               &    28*mN**2*MR**2*p1k*pq**3 + 16*mN*MR**3*p1k*pq**3 - 4*MR**4*p1k*pq**3 -                       &
               &    32*mN**2*p1k**2*pq**3 - 32*mN*MR*p1k**2*pq**3 + 8*ml**4*p1q*pq**3 -                          &
               &    8*ml**2*mN**2*p1q*pq**3 + 4*ml**2*MR**2*p1q*pq**3 + 32*ml**2*p1k*p1q*pq**3 +                 &
               &    32*mN**2*p1k*p1q*pq**3 + 16*mN*MR*p1k*p1q*pq**3 - 32*MR**2*p1k*p1q*pq**3 +                   &
               &    32*p1k**2*p1q*pq**3 - 32*p1k*p1q**2*pq**3 - 12*ml**4*mN**6*Q2 +                              &
               &    16*ml**2*mN**8*Q2 - 4*mN**10*Q2 + 12*ml**2*mN**7*MR*Q2 - 4*mN**9*MR*Q2 +                     &
               &    7*ml**4*mN**4*MR**2*Q2 - 15*ml**2*mN**6*MR**2*Q2 + 6*mN**8*MR**2*Q2 -                        &
               &    2*ml**4*mN**3*MR**3*Q2 - 24*ml**2*mN**5*MR**3*Q2 + 12*mN**7*MR**3*Q2 +                       &
               &    2*ml**4*mN**2*MR**4*Q2 - 6*ml**2*mN**4*MR**4*Q2 + 8*ml**2*mN**3*MR**5*Q2 -                   &
               &    16*mN**5*MR**5*Q2 + 5*ml**2*mN**2*MR**6*Q2 - 10*mN**4*MR**6*Q2 -                             &
               &    48*ml**2*mN**6*p1k*Q2 - 24*ml**2*mN**5*MR*p1k*Q2 + 52*ml**2*mN**4*MR**2*p1k*Q2 +             &
               &    20*ml**2*mN**3*MR**3*p1k*Q2 - 16*ml**2*mN**2*MR**4*p1k*Q2 -                                  &
               &    8*ml**2*mN*MR**5*p1k*Q2 - 48*mN**6*p1k**2*Q2 - 48*mN**5*MR*p1k**2*Q2 +                       &
               &    76*mN**4*MR**2*p1k**2*Q2 + 48*mN**3*MR**3*p1k**2*Q2 - 40*mN**2*MR**4*p1k**2*Q2 -             &
               &    16*mN*MR**5*p1k**2*Q2 + 12*ml**4*mN**4*p1p*Q2 - 16*ml**2*mN**6*p1p*Q2 +                      &
               &    4*mN**8*p1p*Q2 - 24*ml**4*mN**3*MR*p1p*Q2 + 4*ml**2*mN**5*MR*p1p*Q2 +                        &
               &    4*mN**7*MR*p1p*Q2 - 5*ml**4*mN**2*MR**2*p1p*Q2 + 15*ml**2*mN**4*MR**2*p1p*Q2 -               &
               &    4*mN**6*MR**2*p1p*Q2 + 30*ml**4*mN*MR**3*p1p*Q2 - 8*ml**2*mN**3*MR**3*p1p*Q2 -               &
               &    12*mN**5*MR**3*p1p*Q2 + 10*ml**4*MR**4*p1p*Q2 - 2*ml**2*mN**2*MR**4*p1p*Q2 -                 &
               &    12*mN**4*MR**4*p1p*Q2 + 8*ml**2*mN*MR**5*p1p*Q2 + 3*ml**2*MR**6*p1p*Q2 +                     &
               &    4*mN**2*MR**6*p1p*Q2 + 48*ml**2*mN**4*p1k*p1p*Q2 - 64*ml**2*mN**3*MR*p1k*p1p*Q2 -            &
               &    56*ml**2*mN**2*MR**2*p1k*p1p*Q2 + 100*ml**2*mN*MR**3*p1k*p1p*Q2 +                            &
               &    64*ml**2*MR**4*p1k*p1p*Q2 + 48*mN**4*p1k**2*p1p*Q2 - 32*mN**3*MR*p1k**2*p1p*Q2 -             &
               &    92*mN**2*MR**2*p1k**2*p1p*Q2 + 80*mN*MR**3*p1k**2*p1p*Q2 +                                   &
               &    88*MR**4*p1k**2*p1p*Q2 - 4*ml**4*mN**2*p1p**2*Q2 + 4*mN**6*p1p**2*Q2 +                       &
               &    24*ml**4*mN*MR*p1p**2*Q2 - 20*ml**2*mN**3*MR*p1p**2*Q2 + 4*mN**5*MR*p1p**2*Q2 +              &
               &    14*ml**4*MR**2*p1p**2*Q2 - 8*ml**2*mN**2*MR**2*p1p**2*Q2 -                                   &
               &    10*mN**4*MR**2*p1p**2*Q2 + 32*ml**2*mN*MR**3*p1p**2*Q2 -                                     &
               &    12*mN**3*MR**3*p1p**2*Q2 + 16*ml**2*MR**4*p1p**2*Q2 + 8*mN**2*MR**4*p1p**2*Q2 +              &
               &    16*mN*MR**5*p1p**2*Q2 + 6*MR**6*p1p**2*Q2 - 16*ml**2*mN**2*p1k*p1p**2*Q2 +                   &
               &    88*ml**2*mN*MR*p1k*p1p**2*Q2 + 68*ml**2*MR**2*p1k*p1p**2*Q2 -                                &
               &    16*mN**2*p1k**2*p1p**2*Q2 + 80*mN*MR*p1k**2*p1p**2*Q2 +                                      &
               &    80*MR**2*p1k**2*p1p**2*Q2 + 4*ml**4*p1p**3*Q2 - 4*mN**4*p1p**3*Q2 +                          &
               &    4*ml**2*mN*MR*p1p**3*Q2 - 4*mN**3*MR*p1p**3*Q2 + 8*ml**2*MR**2*p1p**3*Q2 +                   &
               &    8*mN**2*MR**2*p1p**3*Q2 + 12*mN*MR**3*p1p**3*Q2 + 4*MR**4*p1p**3*Q2 +                        &
               &    16*ml**2*p1k*p1p**3*Q2 + 16*p1k**2*p1p**3*Q2 + 48*ml**4*mN**4*p1q*Q2 -                       &
               &    72*ml**2*mN**6*p1q*Q2 + 24*mN**8*p1q*Q2 + 8*ml**4*mN**3*MR*p1q*Q2 -                          &
               &    60*ml**2*mN**5*MR*p1q*Q2 + 28*mN**7*MR*p1q*Q2 - 18*ml**4*mN**2*MR**2*p1q*Q2 +                &
               &    30*ml**2*mN**4*MR**2*p1q*Q2 - 24*mN**6*MR**2*p1q*Q2 - 8*ml**4*mN*MR**3*p1q*Q2 +              &
               &    70*ml**2*mN**3*MR**3*p1q*Q2 - 58*mN**5*MR**3*p1q*Q2 +                                        &
               &    28*ml**2*mN**2*MR**4*p1q*Q2 - 6*mN**4*MR**4*p1q*Q2 - 4*ml**2*mN*MR**5*p1q*Q2 +               &
               &    34*mN**3*MR**5*p1q*Q2 + 10*mN**2*MR**6*p1q*Q2 + 192*ml**2*mN**4*p1k*p1q*Q2 +                 &
               &    48*mN**6*p1k*p1q*Q2 + 128*ml**2*mN**3*MR*p1k*p1q*Q2 + 48*mN**5*MR*p1k*p1q*Q2 -               &
               &    120*ml**2*mN**2*MR**2*p1k*p1q*Q2 - 76*mN**4*MR**2*p1k*p1q*Q2 -                               &
               &    96*ml**2*mN*MR**3*p1k*p1q*Q2 - 48*mN**3*MR**3*p1k*p1q*Q2 +                                   &
               &    40*mN**2*MR**4*p1k*p1q*Q2 + 16*mN*MR**5*p1k*p1q*Q2 + 192*mN**4*p1k**2*p1q*Q2 +               &
               &    224*mN**3*MR*p1k**2*p1q*Q2 - 168*mN**2*MR**2*p1k**2*p1q*Q2 -                                 &
               &    160*mN*MR**3*p1k**2*p1q*Q2 - 32*ml**4*mN**2*p1p*p1q*Q2 +                                     &
               &    48*ml**2*mN**4*p1p*p1q*Q2 - 16*mN**6*p1p*p1q*Q2 + 52*ml**4*mN*MR*p1p*p1q*Q2 -                &
               &    8*ml**2*mN**3*MR*p1p*p1q*Q2 - 20*mN**5*MR*p1p*p1q*Q2 +                                       &
               &    28*ml**4*MR**2*p1p*p1q*Q2 - 32*ml**2*mN**2*MR**2*p1p*p1q*Q2 +                                &
               &    22*mN**4*MR**2*p1p*p1q*Q2 - 6*ml**2*mN*MR**3*p1p*p1q*Q2 +                                    &
               &    58*mN**3*MR**3*p1p*p1q*Q2 + 12*mN**2*MR**4*p1p*p1q*Q2 - 26*mN*MR**5*p1p*p1q*Q2 -             &
               &    6*MR**6*p1p*p1q*Q2 - 128*ml**2*mN**2*p1k*p1p*p1q*Q2 - 48*mN**4*p1k*p1p*p1q*Q2 +              &
               &    144*ml**2*mN*MR*p1k*p1p*p1q*Q2 + 32*mN**3*MR*p1k*p1p*p1q*Q2 +                                &
               &    160*ml**2*MR**2*p1k*p1p*p1q*Q2 + 92*mN**2*MR**2*p1k*p1p*p1q*Q2 -                             &
               &    80*mN*MR**3*p1k*p1p*p1q*Q2 - 88*MR**4*p1k*p1p*p1q*Q2 -                                       &
               &    128*mN**2*p1k**2*p1p*p1q*Q2 + 80*mN*MR*p1k**2*p1p*p1q*Q2 +                                   &
               &    208*MR**2*p1k**2*p1p*p1q*Q2 + 16*ml**4*p1p**2*p1q*Q2 -                                       &
               &    8*ml**2*mN**2*p1p**2*p1q*Q2 - 8*mN**4*p1p**2*p1q*Q2 +                                        &
               &    8*ml**2*mN*MR*p1p**2*p1q*Q2 + 6*ml**2*MR**2*p1p**2*p1q*Q2 +                                  &
               &    10*mN**2*MR**2*p1p**2*p1q*Q2 - 2*MR**4*p1p**2*p1q*Q2 +                                       &
               &    64*ml**2*p1k*p1p**2*p1q*Q2 + 16*mN**2*p1k*p1p**2*p1q*Q2 -                                    &
               &    80*mN*MR*p1k*p1p**2*p1q*Q2 - 80*MR**2*p1k*p1p**2*p1q*Q2 +                                    &
               &    64*p1k**2*p1p**2*p1q*Q2 - 8*mN*MR*p1p**3*p1q*Q2 - 8*MR**2*p1p**3*p1q*Q2 -                    &
               &    16*p1k*p1p**3*p1q*Q2 - 48*ml**4*mN**2*p1q**2*Q2 + 96*ml**2*mN**4*p1q**2*Q2 -                 &
               &    48*mN**6*p1q**2*Q2 - 16*ml**4*mN*MR*p1q**2*Q2 + 80*ml**2*mN**3*MR*p1q**2*Q2 -                &
               &    64*mN**5*MR*p1q**2*Q2 + 8*ml**2*mN**2*MR**2*p1q**2*Q2 +                                      &
               &    16*mN**4*MR**2*p1q**2*Q2 - 16*ml**2*mN*MR**3*p1q**2*Q2 +                                     &
               &    72*mN**3*MR**3*p1q**2*Q2 + 28*mN**2*MR**4*p1q**2*Q2 -                                        &
               &    192*ml**2*mN**2*p1k*p1q**2*Q2 - 192*mN**4*p1k*p1q**2*Q2 -                                    &
               &    160*ml**2*mN*MR*p1k*p1q**2*Q2 - 224*mN**3*MR*p1k*p1q**2*Q2 +                                 &
               &    168*mN**2*MR**2*p1k*p1q**2*Q2 + 160*mN*MR**3*p1k*p1q**2*Q2 -                                 &
               &    192*mN**2*p1k**2*p1q**2*Q2 - 256*mN*MR*p1k**2*p1q**2*Q2 +                                    &
               &    16*ml**4*p1p*p1q**2*Q2 - 32*ml**2*mN**2*p1p*p1q**2*Q2 + 16*mN**4*p1p*p1q**2*Q2 -             &
               &    8*ml**2*mN*MR*p1p*p1q**2*Q2 + 40*mN**3*MR*p1p*p1q**2*Q2 -                                    &
               &    12*ml**2*MR**2*p1p*p1q**2*Q2 - 12*mN**2*MR**2*p1p*p1q**2*Q2 -                                &
               &    68*mN*MR**3*p1p*p1q**2*Q2 - 24*MR**4*p1p*p1q**2*Q2 + 64*ml**2*p1k*p1p*p1q**2*Q2 +            &
               &    128*mN**2*p1k*p1p*p1q**2*Q2 - 80*mN*MR*p1k*p1p*p1q**2*Q2 -                                   &
               &    208*MR**2*p1k*p1p*p1q**2*Q2 + 64*p1k**2*p1p*p1q**2*Q2 -                                      &
               &    32*mN*MR*p1p**2*p1q**2*Q2 - 28*MR**2*p1p**2*p1q**2*Q2 - 64*p1k*p1p**2*p1q**2*Q2 -            &
               &    32*ml**2*mN**2*p1q**3*Q2 + 32*mN**4*p1q**3*Q2 - 16*ml**2*mN*MR*p1q**3*Q2 +                   &
               &    48*mN**3*MR*p1q**3*Q2 + 16*mN**2*MR**2*p1q**3*Q2 + 192*mN**2*p1k*p1q**3*Q2 +                 &
               &    256*mN*MR*p1k*p1q**3*Q2 - 32*mN*MR*p1p*p1q**3*Q2 - 24*MR**2*p1p*p1q**3*Q2 -                  &
               &    64*p1k*p1p*p1q**3*Q2 + 24*ml**2*mN**5*MR*pk*Q2 - 32*ml**2*mN**4*MR**2*pk*Q2 -                &
               &    44*ml**2*mN**3*MR**3*pk*Q2 + 8*ml**2*mN*MR**5*pk*Q2 + 48*mN**5*MR*p1k*pk*Q2 -                &
               &    40*mN**4*MR**2*p1k*pk*Q2 - 136*mN**3*MR**3*p1k*pk*Q2 - 16*mN**2*MR**4*p1k*pk*Q2 +            &
               &    40*mN*MR**5*p1k*pk*Q2 + 8*MR**6*p1k*pk*Q2 - 32*ml**2*mN**3*MR*p1p*pk*Q2 +                    &
               &    8*ml**2*mN**2*MR**2*p1p*pk*Q2 - 12*ml**2*mN*MR**3*p1p*pk*Q2 -                                &
               &    16*ml**2*MR**4*p1p*pk*Q2 - 64*mN**3*MR*p1k*p1p*pk*Q2 -                                       &
               &    16*mN**2*MR**2*p1k*p1p*pk*Q2 - 24*MR**4*p1k*p1p*pk*Q2 +                                      &
               &    8*ml**2*mN*MR*p1p**2*pk*Q2 - 8*ml**2*MR**2*p1p**2*pk*Q2 +                                    &
               &    16*mN*MR*p1k*p1p**2*pk*Q2 - 8*MR**2*p1k*p1p**2*pk*Q2 -                                       &
               &    96*ml**2*mN**3*MR*p1q*pk*Q2 - 24*mN**5*MR*p1q*pk*Q2 +                                        &
               &    56*ml**2*mN**2*MR**2*p1q*pk*Q2 + 20*mN**4*MR**2*p1q*pk*Q2 +                                  &
               &    84*ml**2*mN*MR**3*p1q*pk*Q2 + 68*mN**3*MR**3*p1q*pk*Q2 +                                     &
               &    8*ml**2*MR**4*p1q*pk*Q2 + 8*mN**2*MR**4*p1q*pk*Q2 - 20*mN*MR**5*p1q*pk*Q2 -                  &
               &    4*MR**6*p1q*pk*Q2 - 192*mN**3*MR*p1k*p1q*pk*Q2 + 16*mN**2*MR**2*p1k*p1q*pk*Q2 +              &
               &    264*mN*MR**3*p1k*p1q*pk*Q2 + 80*MR**4*p1k*p1q*pk*Q2 +                                        &
               &    64*ml**2*mN*MR*p1p*p1q*pk*Q2 + 32*mN**3*MR*p1p*p1q*pk*Q2 -                                   &
               &    4*ml**2*MR**2*p1p*p1q*pk*Q2 + 8*mN**2*MR**2*p1p*p1q*pk*Q2 +                                  &
               &    12*MR**4*p1p*p1q*pk*Q2 + 128*mN*MR*p1k*p1p*p1q*pk*Q2 +                                       &
               &    56*MR**2*p1k*p1p*p1q*pk*Q2 - 8*mN*MR*p1p**2*p1q*pk*Q2 +                                      &
               &    4*MR**2*p1p**2*p1q*pk*Q2 + 96*ml**2*mN*MR*p1q**2*pk*Q2 +                                     &
               &    96*mN**3*MR*p1q**2*pk*Q2 + 16*ml**2*MR**2*p1q**2*pk*Q2 -                                     &
               &    8*mN**2*MR**2*p1q**2*pk*Q2 - 132*mN*MR**3*p1q**2*pk*Q2 - 40*MR**4*p1q**2*pk*Q2 +             &
               &    192*mN*MR*p1k*p1q**2*pk*Q2 + 128*MR**2*p1k*p1q**2*pk*Q2 -                                    &
               &    64*mN*MR*p1p*p1q**2*pk*Q2 - 28*MR**2*p1p*p1q**2*pk*Q2 - 96*mN*MR*p1q**3*pk*Q2 -              &
               &    64*MR**2*p1q**3*pk*Q2 + 16*mN**3*MR**3*pk**2*Q2 - 4*mN**2*MR**4*pk**2*Q2 -                   &
               &    8*mN*MR**5*pk**2*Q2 - 8*mN*MR**3*p1p*pk**2*Q2 + 4*MR**4*p1p*pk**2*Q2 -                       &
               &    32*mN*MR**3*p1q*pk**2*Q2 - 8*MR**4*p1q*pk**2*Q2 + 16*ml**4*mN**3*MR*pq*Q2 -                  &
               &    4*ml**2*mN**5*MR*pq*Q2 - 4*mN**7*MR*pq*Q2 - 16*ml**4*mN**2*MR**2*pq*Q2 +                     &
               &    32*ml**2*mN**4*MR**2*pq*Q2 - 14*mN**6*MR**2*pq*Q2 - 24*ml**4*mN*MR**3*pq*Q2 +                &
               &    22*ml**2*mN**3*MR**3*pq*Q2 + 2*mN**5*MR**3*pq*Q2 - 20*ml**2*mN**2*MR**4*pq*Q2 +              &
               &    24*mN**4*MR**4*pq*Q2 - 12*ml**2*mN*MR**5*pq*Q2 + 10*mN**3*MR**5*pq*Q2 -                      &
               &    2*mN**2*MR**6*pq*Q2 + 32*ml**2*mN**3*MR*p1k*pq*Q2 - 24*mN**5*MR*p1k*pq*Q2 -                  &
               &    28*ml**2*mN**2*MR**2*p1k*pq*Q2 + 20*mN**4*MR**2*p1k*pq*Q2 -                                  &
               &    76*ml**2*mN*MR**3*p1k*pq*Q2 + 68*mN**3*MR**3*p1k*pq*Q2 -                                     &
               &    24*ml**2*MR**4*p1k*pq*Q2 + 8*mN**2*MR**4*p1k*pq*Q2 - 20*mN*MR**5*p1k*pq*Q2 -                 &
               &    4*MR**6*p1k*pq*Q2 + 8*mN**2*MR**2*p1k**2*pq*Q2 - 56*mN*MR**3*p1k**2*pq*Q2 -                  &
               &    48*MR**4*p1k**2*pq*Q2 - 8*ml**4*mN**2*p1p*pq*Q2 + 24*ml**2*mN**4*p1p*pq*Q2 -                 &
               &    16*mN**6*p1p*pq*Q2 - 44*ml**4*mN*MR*p1p*pq*Q2 + 48*ml**2*mN**3*MR*p1p*pq*Q2 -                &
               &    12*mN**5*MR*p1p*pq*Q2 - 14*ml**4*MR**2*p1p*pq*Q2 -                                           &
               &    24*ml**2*mN**2*MR**2*p1p*pq*Q2 + 38*mN**4*MR**2*p1p*pq*Q2 -                                  &
               &    54*ml**2*mN*MR**3*p1p*pq*Q2 + 26*mN**3*MR**3*p1p*pq*Q2 -                                     &
               &    8*ml**2*MR**4*p1p*pq*Q2 - 24*mN**2*MR**4*p1p*pq*Q2 - 18*mN*MR**5*p1p*pq*Q2 -                 &
               &    2*MR**6*p1p*pq*Q2 - 32*ml**2*mN**2*p1k*p1p*pq*Q2 -                                           &
               &    160*ml**2*mN*MR*p1k*p1p*pq*Q2 + 32*mN**3*MR*p1k*p1p*pq*Q2 -                                  &
               &    80*ml**2*MR**2*p1k*p1p*pq*Q2 + 8*mN**2*MR**2*p1k*p1p*pq*Q2 +                                 &
               &    12*MR**4*p1k*p1p*pq*Q2 - 32*mN**2*p1k**2*p1p*pq*Q2 - 144*mN*MR*p1k**2*p1p*pq*Q2 -            &
               &    104*MR**2*p1k**2*p1p*pq*Q2 - 8*ml**4*p1p**2*pq*Q2 - 8*ml**2*mN**2*p1p**2*pq*Q2 +             &
               &    16*mN**4*p1p**2*pq*Q2 - 16*ml**2*mN*MR*p1p**2*pq*Q2 + 16*mN**3*MR*p1p**2*pq*Q2 -             &
               &    12*ml**2*MR**2*p1p**2*pq*Q2 - 24*mN**2*MR**2*p1p**2*pq*Q2 -                                  &
               &    28*mN*MR**3*p1p**2*pq*Q2 - 4*MR**4*p1p**2*pq*Q2 - 32*ml**2*p1k*p1p**2*pq*Q2 -                &
               &    8*mN*MR*p1k*p1p**2*pq*Q2 + 4*MR**2*p1k*p1p**2*pq*Q2 - 32*p1k**2*p1p**2*pq*Q2 -               &
               &    16*ml**4*mN**2*p1q*pq*Q2 + 24*ml**2*mN**4*p1q*pq*Q2 - 8*mN**6*p1q*pq*Q2 -                    &
               &    36*ml**4*mN*MR*p1q*pq*Q2 + 24*ml**2*mN**3*MR*p1q*pq*Q2 + 12*mN**5*MR*p1q*pq*Q2 +             &
               &    16*ml**4*MR**2*p1q*pq*Q2 - 82*ml**2*mN**2*MR**2*p1q*pq*Q2 +                                  &
               &    50*mN**4*MR**2*p1q*pq*Q2 - 32*ml**2*mN*MR**3*p1q*pq*Q2 -                                     &
               &    8*mN**3*MR**3*p1q*pq*Q2 + 16*ml**2*MR**4*p1q*pq*Q2 - 34*mN**2*MR**4*p1q*pq*Q2 +              &
               &    4*mN*MR**5*p1q*pq*Q2 + 4*MR**6*p1q*pq*Q2 - 64*ml**2*mN**2*p1k*p1q*pq*Q2 -                    &
               &    80*ml**2*mN*MR*p1k*p1q*pq*Q2 + 96*mN**3*MR*p1k*p1q*pq*Q2 +                                   &
               &    16*ml**2*MR**2*p1k*p1q*pq*Q2 - 16*mN**2*MR**2*p1k*p1q*pq*Q2 -                                &
               &    76*mN*MR**3*p1k*p1q*pq*Q2 + 8*MR**4*p1k*p1q*pq*Q2 - 64*mN**2*p1k**2*p1q*pq*Q2 -              &
               &    16*mN*MR*p1k**2*p1q*pq*Q2 - 32*MR**2*p1k**2*p1q*pq*Q2 -                                      &
               &    48*ml**2*mN**2*p1p*p1q*pq*Q2 + 48*mN**4*p1p*p1q*pq*Q2 -                                      &
               &    40*ml**2*mN*MR*p1p*p1q*pq*Q2 + 24*mN**3*MR*p1p*p1q*pq*Q2 +                                   &
               &    34*ml**2*MR**2*p1p*p1q*pq*Q2 - 58*mN**2*MR**2*p1p*p1q*pq*Q2 +                                &
               &    26*MR**4*p1p*p1q*pq*Q2 + 32*mN**2*p1k*p1p*p1q*pq*Q2 +                                        &
               &    80*mN*MR*p1k*p1p*p1q*pq*Q2 + 76*MR**2*p1k*p1p*p1q*pq*Q2 +                                    &
               &    8*ml**2*p1p**2*p1q*pq*Q2 - 8*mN**2*p1p**2*p1q*pq*Q2 + 16*mN*MR*p1p**2*p1q*pq*Q2 +            &
               &    28*MR**2*p1p**2*p1q*pq*Q2 + 32*p1k*p1p**2*p1q*pq*Q2 + 32*ml**4*p1q**2*pq*Q2 -                &
               &    64*ml**2*mN**2*p1q**2*pq*Q2 + 32*mN**4*p1q**2*pq*Q2 -                                        &
               &    24*ml**2*mN*MR*p1q**2*pq*Q2 - 8*mN**3*MR*p1q**2*pq*Q2 +                                      &
               &    48*ml**2*MR**2*p1q**2*pq*Q2 - 56*mN**2*MR**2*p1q**2*pq*Q2 +                                  &
               &    8*mN*MR**3*p1q**2*pq*Q2 + 16*MR**4*p1q**2*pq*Q2 + 128*ml**2*p1k*p1q**2*pq*Q2 +               &
               &    64*mN**2*p1k*p1q**2*pq*Q2 - 80*mN*MR*p1k*p1q**2*pq*Q2 -                                      &
               &    32*MR**2*p1k*p1q**2*pq*Q2 + 128*p1k**2*p1q**2*pq*Q2 + 32*ml**2*p1p*p1q**2*pq*Q2 -            &
               &    32*mN**2*p1p*p1q**2*pq*Q2 + 32*mN*MR*p1p*p1q**2*pq*Q2 +                                      &
               &    60*MR**2*p1p*p1q**2*pq*Q2 + 32*ml**2*p1q**3*pq*Q2 - 32*mN**2*p1q**3*pq*Q2 +                  &
               &    16*MR**2*p1q**3*pq*Q2 - 128*p1k*p1q**3*pq*Q2 + 32*ml**2*mN**3*MR*pk*pq*Q2 +                  &
               &    8*ml**2*mN**2*MR**2*pk*pq*Q2 + 8*ml**2*mN*MR**3*pk*pq*Q2 -                                   &
               &    16*mN**3*MR**3*pk*pq*Q2 + 4*mN**2*MR**4*pk*pq*Q2 + 8*mN*MR**5*pk*pq*Q2 +                     &
               &    64*mN**3*MR*p1k*pk*pq*Q2 + 48*mN**2*MR**2*p1k*pk*pq*Q2 -                                     &
               &    8*mN*MR**3*p1k*pk*pq*Q2 - 8*MR**4*p1k*pk*pq*Q2 - 16*ml**2*mN*MR*p1p*pk*pq*Q2 +               &
               &    8*ml**2*MR**2*p1p*pk*pq*Q2 + 8*mN*MR**3*p1p*pk*pq*Q2 - 4*MR**4*p1p*pk*pq*Q2 -                &
               &    32*mN*MR*p1k*p1p*pk*pq*Q2 - 64*ml**2*mN*MR*p1q*pk*pq*Q2 -                                    &
               &    32*mN**3*MR*p1q*pk*pq*Q2 - 28*ml**2*MR**2*p1q*pk*pq*Q2 -                                     &
               &    24*mN**2*MR**2*p1q*pk*pq*Q2 + 36*mN*MR**3*p1q*pk*pq*Q2 + 12*MR**4*p1q*pk*pq*Q2 -             &
               &    128*mN*MR*p1k*p1q*pk*pq*Q2 - 120*MR**2*p1k*p1q*pk*pq*Q2 +                                    &
               &    16*mN*MR*p1p*p1q*pk*pq*Q2 + 64*mN*MR*p1q**2*pk*pq*Q2 + 60*MR**2*p1q**2*pk*pq*Q2 +            &
               &    8*mN*MR**3*pk**2*pq*Q2 + 12*ml**4*mN**2*pq**2*Q2 - 24*ml**2*mN**4*pq**2*Q2 +                 &
               &    12*mN**6*pq**2*Q2 + 20*ml**4*mN*MR*pq**2*Q2 - 28*ml**2*mN**3*MR*pq**2*Q2 +                   &
               &    8*mN**5*MR*pq**2*Q2 + 24*ml**2*mN**2*MR**2*pq**2*Q2 - 24*mN**4*MR**2*pq**2*Q2 +              &
               &    24*ml**2*mN*MR**3*pq**2*Q2 - 16*mN**3*MR**3*pq**2*Q2 + 8*mN**2*MR**4*pq**2*Q2 +              &
               &    4*mN*MR**5*pq**2*Q2 + 48*ml**2*mN**2*p1k*pq**2*Q2 + 72*ml**2*mN*MR*p1k*pq**2*Q2 -            &
               &    32*mN**3*MR*p1k*pq**2*Q2 + 12*ml**2*MR**2*p1k*pq**2*Q2 -                                     &
               &    24*mN**2*MR**2*p1k*pq**2*Q2 + 4*mN*MR**3*p1k*pq**2*Q2 + 4*MR**4*p1k*pq**2*Q2 +               &
               &    48*mN**2*p1k**2*pq**2*Q2 + 64*mN*MR*p1k**2*pq**2*Q2 + 24*MR**2*p1k**2*pq**2*Q2 +             &
               &    4*ml**4*p1p*pq**2*Q2 + 16*ml**2*mN**2*p1p*pq**2*Q2 - 20*mN**4*p1p*pq**2*Q2 +                 &
               &    20*ml**2*mN*MR*p1p*pq**2*Q2 - 20*mN**3*MR*p1p*pq**2*Q2 +                                     &
               &    4*ml**2*MR**2*p1p*pq**2*Q2 + 20*mN**2*MR**2*p1p*pq**2*Q2 +                                   &
               &    20*mN*MR**3*p1p*pq**2*Q2 + 16*ml**2*p1k*p1p*pq**2*Q2 +                                       &
               &    16*mN*MR*p1k*p1p*pq**2*Q2 + 16*p1k**2*p1p*pq**2*Q2 - 16*ml**4*p1q*pq**2*Q2 +                 &
               &    56*ml**2*mN**2*p1q*pq**2*Q2 - 40*mN**4*p1q*pq**2*Q2 +                                        &
               &    32*ml**2*mN*MR*p1q*pq**2*Q2 - 24*mN**3*MR*p1q*pq**2*Q2 -                                     &
               &    24*ml**2*MR**2*p1q*pq**2*Q2 + 32*mN**2*MR**2*p1q*pq**2*Q2 +                                  &
               &    4*mN*MR**3*p1q*pq**2*Q2 - 8*MR**4*p1q*pq**2*Q2 - 64*ml**2*p1k*p1q*pq**2*Q2 -                 &
               &    48*mN**2*p1k*p1q*pq**2*Q2 + 36*MR**2*p1k*p1q*pq**2*Q2 - 64*p1k**2*p1q*pq**2*Q2 -             &
               &    16*ml**2*p1p*p1q*pq**2*Q2 + 16*mN**2*p1p*p1q*pq**2*Q2 -                                      &
               &    8*mN*MR*p1p*p1q*pq**2*Q2 - 24*MR**2*p1p*p1q*pq**2*Q2 - 16*p1k*p1p*p1q*pq**2*Q2 -             &
               &    32*ml**2*p1q**2*pq**2*Q2 + 32*mN**2*p1q**2*pq**2*Q2 - 16*MR**2*p1q**2*pq**2*Q2 +             &
               &    64*p1k*p1q**2*pq**2*Q2 + 8*ml**2*mN*MR*pk*pq**2*Q2 - 8*mN*MR**3*pk*pq**2*Q2 +                &
               &    16*mN*MR*p1k*pk*pq**2*Q2 + 8*MR**2*p1k*pk*pq**2*Q2 - 8*mN*MR*p1q*pk*pq**2*Q2 -               &
               &    4*MR**2*p1q*pk*pq**2*Q2 - 8*ml**2*mN**2*pq**3*Q2 + 8*mN**4*pq**3*Q2 -                        &
               &    8*ml**2*mN*MR*pq**3*Q2 + 8*mN**3*MR*pq**3*Q2 - 4*mN**2*MR**2*pq**3*Q2 -                      &
               &    4*mN*MR**3*pq**3*Q2 - 8*mN*MR*p1k*pq**3*Q2 - 4*MR**2*p1k*pq**3*Q2 +                          &
               &    8*ml**2*p1q*pq**3*Q2 - 8*mN**2*p1q*pq**3*Q2 + 4*MR**2*p1q*pq**3*Q2 +                         &
               &    12*ml**4*mN**4*Q2**2 - 24*ml**2*mN**6*Q2**2 + 12*mN**8*Q2**2 -                               &
               &    12*ml**2*mN**5*MR*Q2**2 + 12*mN**7*MR*Q2**2 - 5*ml**4*mN**2*MR**2*Q2**2 +                    &
               &    12*ml**2*mN**4*MR**2*Q2**2 - 13*mN**6*MR**2*Q2**2 + 12*ml**2*mN**3*MR**3*Q2**2 -             &
               &    24*mN**5*MR**3*Q2**2 + 6*ml**2*mN**2*MR**4*Q2**2 + mN**4*MR**4*Q2**2 +                       &
               &    10*mN**3*MR**5*Q2**2 + 2*mN**2*MR**6*Q2**2 + 48*ml**2*mN**4*p1k*Q2**2 +                      &
               &    24*ml**2*mN**3*MR*p1k*Q2**2 - 32*ml**2*mN**2*MR**2*p1k*Q2**2 -                               &
               &    16*ml**2*mN*MR**3*p1k*Q2**2 + 48*mN**4*p1k**2*Q2**2 + 48*mN**3*MR*p1k**2*Q2**2 -             &
               &    44*mN**2*MR**2*p1k**2*Q2**2 - 32*mN*MR**3*p1k**2*Q2**2 -                                     &
               &    12*ml**4*mN**2*p1p*Q2**2 + 24*ml**2*mN**4*p1p*Q2**2 - 12*mN**6*p1p*Q2**2 +                   &
               &    16*ml**4*mN*MR*p1p*Q2**2 - 20*ml**2*mN**3*MR*p1p*Q2**2 - 4*mN**5*MR*p1p*Q2**2 +              &
               &    11*ml**4*MR**2*p1p*Q2**2 - 22*ml**2*mN**2*MR**2*p1p*Q2**2 +                                  &
               &    23*mN**4*MR**2*p1p*Q2**2 + 20*ml**2*mN*MR**3*p1p*Q2**2 +                                     &
               &    16*mN**3*MR**3*p1p*Q2**2 + 10*ml**2*MR**4*p1p*Q2**2 - 9*mN**2*MR**4*p1p*Q2**2 -              &
               &    6*mN*MR**5*p1p*Q2**2 - 48*ml**2*mN**2*p1k*p1p*Q2**2 +                                        &
               &    48*ml**2*mN*MR*p1k*p1p*Q2**2 + 56*ml**2*MR**2*p1k*p1p*Q2**2 -                                &
               &    48*mN**2*p1k**2*p1p*Q2**2 + 32*mN*MR*p1k**2*p1p*Q2**2 +                                      &
               &    68*MR**2*p1k**2*p1p*Q2**2 + 8*ml**4*p1p**2*Q2**2 - 12*ml**2*mN**2*p1p**2*Q2**2 +             &
               &    4*mN**4*p1p**2*Q2**2 + 16*ml**2*mN*MR*p1p**2*Q2**2 - 4*mN**3*MR*p1p**2*Q2**2 +               &
               &    18*ml**2*MR**2*p1p**2*Q2**2 - 14*mN**2*MR**2*p1p**2*Q2**2 +                                  &
               &    8*mN*MR**3*p1p**2*Q2**2 + 10*MR**4*p1p**2*Q2**2 + 32*ml**2*p1k*p1p**2*Q2**2 +                &
               &    32*p1k**2*p1p**2*Q2**2 + 4*ml**2*p1p**3*Q2**2 - 4*mN**2*p1p**3*Q2**2 -                       &
               &    4*mN*MR*p1p**3*Q2**2 + 4*MR**2*p1p**3*Q2**2 - 24*ml**4*mN**2*p1q*Q2**2 +                     &
               &    72*ml**2*mN**4*p1q*Q2**2 - 48*mN**6*p1q*Q2**2 - 4*ml**4*mN*MR*p1q*Q2**2 +                    &
               &    36*ml**2*mN**3*MR*p1q*Q2**2 - 56*mN**5*MR*p1q*Q2**2 -                                        &
               &    6*ml**2*mN**2*MR**2*p1q*Q2**2 + 18*mN**4*MR**2*p1q*Q2**2 -                                   &
               &    8*ml**2*mN*MR**3*p1q*Q2**2 + 64*mN**3*MR**3*p1q*Q2**2 +                                      &
               &    22*mN**2*MR**4*p1q*Q2**2 - 96*ml**2*mN**2*p1k*p1q*Q2**2 -                                    &
               &    48*mN**4*p1k*p1q*Q2**2 - 64*ml**2*mN*MR*p1k*p1q*Q2**2 -                                      &
               &    48*mN**3*MR*p1k*p1q*Q2**2 + 44*mN**2*MR**2*p1k*p1q*Q2**2 +                                   &
               &    32*mN*MR**3*p1k*p1q*Q2**2 - 96*mN**2*p1k**2*p1q*Q2**2 -                                      &
               &    112*mN*MR*p1k**2*p1q*Q2**2 + 16*ml**4*p1p*p1q*Q2**2 -                                        &
               &    48*ml**2*mN**2*p1p*p1q*Q2**2 + 32*mN**4*p1p*p1q*Q2**2 +                                      &
               &    28*ml**2*mN*MR*p1p*p1q*Q2**2 + 28*mN**3*MR*p1p*p1q*Q2**2 +                                   &
               &    22*ml**2*MR**2*p1p*p1q*Q2**2 - 40*mN**2*MR**2*p1p*p1q*Q2**2 -                                &
               &    54*mN*MR**3*p1p*p1q*Q2**2 - 12*MR**4*p1p*p1q*Q2**2 + 64*ml**2*p1k*p1p*p1q*Q2**2 +            &
               &    48*mN**2*p1k*p1p*p1q*Q2**2 - 32*mN*MR*p1k*p1p*p1q*Q2**2 -                                    &
               &    68*MR**2*p1k*p1p*p1q*Q2**2 + 64*p1k**2*p1p*p1q*Q2**2 +                                       &
               &    16*ml**2*p1p**2*p1q*Q2**2 - 16*mN**2*p1p**2*p1q*Q2**2 -                                      &
               &    32*mN*MR*p1p**2*p1q*Q2**2 - 6*MR**2*p1p**2*p1q*Q2**2 - 32*p1k*p1p**2*p1q*Q2**2 -             &
               &    48*ml**2*mN**2*p1q**2*Q2**2 + 48*mN**4*p1q**2*Q2**2 -                                        &
               &    16*ml**2*mN*MR*p1q**2*Q2**2 + 64*mN**3*MR*p1q**2*Q2**2 +                                     &
               &    24*mN**2*MR**2*p1q**2*Q2**2 + 96*mN**2*p1k*p1q**2*Q2**2 +                                    &
               &    112*mN*MR*p1k*p1q**2*Q2**2 + 16*ml**2*p1p*p1q**2*Q2**2 -                                     &
               &    16*mN**2*p1p*p1q**2*Q2**2 - 48*mN*MR*p1p*p1q**2*Q2**2 -                                      &
               &    24*MR**2*p1p*p1q**2*Q2**2 - 64*p1k*p1p*p1q**2*Q2**2 -                                        &
               &    24*ml**2*mN**3*MR*pk*Q2**2 + 16*ml**2*mN**2*MR**2*pk*Q2**2 +                                 &
               &    16*ml**2*mN*MR**3*pk*Q2**2 - 48*mN**3*MR*p1k*pk*Q2**2 +                                      &
               &    8*mN**2*MR**2*p1k*pk*Q2**2 + 56*mN*MR**3*p1k*pk*Q2**2 + 16*MR**4*p1k*pk*Q2**2 +              &
               &    16*ml**2*mN*MR*p1p*pk*Q2**2 - 8*ml**2*MR**2*p1p*pk*Q2**2 +                                   &
               &    32*mN*MR*p1k*p1p*pk*Q2**2 + 48*ml**2*mN*MR*p1q*pk*Q2**2 +                                    &
               &    24*mN**3*MR*p1q*pk*Q2**2 + 4*ml**2*MR**2*p1q*pk*Q2**2 -                                      &
               &    4*mN**2*MR**2*p1q*pk*Q2**2 - 28*mN*MR**3*p1q*pk*Q2**2 - 8*MR**4*p1q*pk*Q2**2 +               &
               &    96*mN*MR*p1k*p1q*pk*Q2**2 + 56*MR**2*p1k*p1q*pk*Q2**2 -                                      &
               &    16*mN*MR*p1p*p1q*pk*Q2**2 - 48*mN*MR*p1q**2*pk*Q2**2 - 28*MR**2*p1q**2*pk*Q2**2 -            &
               &    8*mN*MR**3*pk**2*Q2**2 - 12*ml**4*mN*MR*pq*Q2**2 + 20*ml**2*mN**3*MR*pq*Q2**2 -              &
               &    18*ml**2*mN**2*MR**2*pq*Q2**2 + 14*mN**4*MR**2*pq*Q2**2 -                                    &
               &    24*ml**2*mN*MR**3*pq*Q2**2 + 4*mN**3*MR**3*pq*Q2**2 - 4*mN**2*MR**4*pq*Q2**2 -               &
               &    32*ml**2*mN*MR*p1k*pq*Q2**2 + 24*mN**3*MR*p1k*pq*Q2**2 -                                     &
               &    12*ml**2*MR**2*p1k*pq*Q2**2 - 4*mN**2*MR**2*p1k*pq*Q2**2 -                                   &
               &    28*mN*MR**3*p1k*pq*Q2**2 - 8*MR**4*p1k*pq*Q2**2 - 16*mN*MR*p1k**2*pq*Q2**2 -                 &
               &    24*MR**2*p1k**2*pq*Q2**2 - 8*ml**4*p1p*pq*Q2**2 + 8*mN**4*p1p*pq*Q2**2 -                     &
               &    36*ml**2*mN*MR*p1p*pq*Q2**2 + 20*mN**3*MR*p1p*pq*Q2**2 -                                     &
               &    14*ml**2*MR**2*p1p*pq*Q2**2 + 2*mN**2*MR**2*p1p*pq*Q2**2 -                                   &
               &    10*mN*MR**3*p1p*pq*Q2**2 - 2*MR**4*p1p*pq*Q2**2 - 32*ml**2*p1k*p1p*pq*Q2**2 -                &
               &    16*mN*MR*p1k*p1p*pq*Q2**2 - 32*p1k**2*p1p*pq*Q2**2 - 8*ml**2*p1p**2*pq*Q2**2 +               &
               &    8*mN**2*p1p**2*pq*Q2**2 + 8*mN*MR*p1p**2*pq*Q2**2 - 4*MR**2*p1p**2*pq*Q2**2 +                &
               &    8*ml**4*p1q*pq*Q2**2 - 24*ml**2*mN**2*p1q*pq*Q2**2 + 16*mN**4*p1q*pq*Q2**2 -                 &
               &    36*ml**2*mN*MR*p1q*pq*Q2**2 + 4*mN**3*MR*p1q*pq*Q2**2 +                                      &
               &    20*ml**2*MR**2*p1q*pq*Q2**2 - 32*mN**2*MR**2*p1q*pq*Q2**2 +                                  &
               &    4*mN*MR**3*p1q*pq*Q2**2 + 8*MR**4*p1q*pq*Q2**2 + 32*ml**2*p1k*p1q*pq*Q2**2 -                 &
               &    32*mN*MR*p1k*p1q*pq*Q2**2 - 4*MR**2*p1k*p1q*pq*Q2**2 + 32*p1k**2*p1q*pq*Q2**2 +              &
               &    32*mN*MR*p1p*p1q*pq*Q2**2 + 30*MR**2*p1p*p1q*pq*Q2**2 + 32*p1k*p1p*p1q*pq*Q2**2 +            &
               &    32*ml**2*p1q**2*pq*Q2**2 - 32*mN**2*p1q**2*pq*Q2**2 + 16*MR**2*p1q**2*pq*Q2**2 -             &
               &    32*p1k*p1q**2*pq*Q2**2 - 16*ml**2*mN*MR*pk*pq*Q2**2 + 8*mN*MR**3*pk*pq*Q2**2 -               &
               &    32*mN*MR*p1k*pk*pq*Q2**2 - 16*MR**2*p1k*pk*pq*Q2**2 + 16*mN*MR*p1q*pk*pq*Q2**2 +             &
               &    8*MR**2*p1q*pk*pq*Q2**2 + 12*ml**2*mN**2*pq**2*Q2**2 - 12*mN**4*pq**2*Q2**2 +                &
               &    20*ml**2*mN*MR*pq**2*Q2**2 - 16*mN**3*MR*pq**2*Q2**2 +                                       &
               &    4*mN**2*MR**2*pq**2*Q2**2 + 4*mN*MR**3*pq**2*Q2**2 + 16*mN*MR*p1k*pq**2*Q2**2 +              &
               &    8*MR**2*p1k*pq**2*Q2**2 + 4*ml**2*p1p*pq**2*Q2**2 - 4*mN**2*p1p*pq**2*Q2**2 -                &
               &    4*mN*MR*p1p*pq**2*Q2**2 - 16*ml**2*p1q*pq**2*Q2**2 + 16*mN**2*p1q*pq**2*Q2**2 -              &
               &    8*MR**2*p1q*pq**2*Q2**2 - 4*ml**4*mN**2*Q2**3 + 16*ml**2*mN**4*Q2**3 -                       &
               &    12*mN**6*Q2**3 + 4*ml**2*mN**3*MR*Q2**3 - 12*mN**5*MR*Q2**3 -                                &
               &    3*ml**2*mN**2*MR**2*Q2**3 + 5*mN**4*MR**2*Q2**3 + 14*mN**3*MR**3*Q2**3 +                     &
               &    4*mN**2*MR**4*Q2**3 - 16*ml**2*mN**2*p1k*Q2**3 - 8*ml**2*mN*MR*p1k*Q2**3 -                   &
               &    16*mN**2*p1k**2*Q2**3 - 16*mN*MR*p1k**2*Q2**3 + 4*ml**4*p1p*Q2**3 -                          &
               &    16*ml**2*mN**2*p1p*Q2**3 + 12*mN**4*p1p*Q2**3 + 12*ml**2*mN*MR*p1p*Q2**3 +                   &
               &    4*mN**3*MR*p1p*Q2**3 + 11*ml**2*MR**2*p1p*Q2**3 - 17*mN**2*MR**2*p1p*Q2**3 -                 &
               &    10*mN*MR**3*p1p*Q2**3 + 16*ml**2*p1k*p1p*Q2**3 + 16*p1k**2*p1p*Q2**3 +                       &
               &    8*ml**2*p1p**2*Q2**3 - 8*mN**2*p1p**2*Q2**3 - 8*mN*MR*p1p**2*Q2**3 +                         &
               &    4*MR**2*p1p**2*Q2**3 - 24*ml**2*mN**2*p1q*Q2**3 + 24*mN**4*p1q*Q2**3 -                       &
               &    4*ml**2*mN*MR*p1q*Q2**3 + 28*mN**3*MR*p1q*Q2**3 + 12*mN**2*MR**2*p1q*Q2**3 +                 &
               &    16*mN**2*p1k*p1q*Q2**3 + 16*mN*MR*p1k*p1q*Q2**3 + 16*ml**2*p1p*p1q*Q2**3 -                   &
               &    16*mN**2*p1p*p1q*Q2**3 - 24*mN*MR*p1p*p1q*Q2**3 - 6*MR**2*p1p*p1q*Q2**3 -                    &
               &    16*p1k*p1p*p1q*Q2**3 + 8*ml**2*mN*MR*pk*Q2**3 + 16*mN*MR*p1k*pk*Q2**3 +                      &
               &    8*MR**2*p1k*pk*Q2**3 - 8*mN*MR*p1q*pk*Q2**3 - 4*MR**2*p1q*pk*Q2**3 -                         &
               &    12*ml**2*mN*MR*pq*Q2**3 + 4*mN**3*MR*pq*Q2**3 - 2*mN**2*MR**2*pq*Q2**3 -                     &
               &    8*mN*MR*p1k*pq*Q2**3 - 4*MR**2*p1k*pq*Q2**3 - 8*ml**2*p1p*pq*Q2**3 +                         &
               &    8*mN**2*p1p*pq*Q2**3 + 8*mN*MR*p1p*pq*Q2**3 + 8*ml**2*p1q*pq*Q2**3 -                         &
               &    8*mN**2*p1q*pq*Q2**3 + 4*MR**2*p1q*pq*Q2**3 - 4*ml**2*mN**2*Q2**4 +                          &
               &    4*mN**4*Q2**4 + 4*mN**3*MR*Q2**4 + 2*mN**2*MR**2*Q2**4 + 4*ml**2*p1p*Q2**4 -                 &
               &    4*mN**2*p1p*Q2**4 - 4*mN*MR*p1p*Q2**4))/(9.*MR**4)/( abs4Sq(p_out-q) - MR**2 )**2





          !nucleon pole for neutron target
          M2Np = 16*(16*f1V*f2V*ml**4*mN**5 - 16*f1A**2*ml**2*mN**6 + 16*f1V**2*ml**2*mN**6 +                    &
               &  16*f2V**2*ml**4*mN**6 - 16*f1V*f2V*ml**4*mN**3*p1p + 16*f1A**2*ml**2*mN**4*p1p -                 &
               &  16*f1V**2*ml**2*mN**4*p1p - 16*f2V**2*ml**4*mN**4*p1p - 16*f1V*f2V*ml**4*mN**3*p1q +             &
               &  16*f1A**2*ml**2*mN**4*p1q - 16*f1V**2*ml**2*mN**4*p1q - 16*f2V**2*ml**4*mN**4*p1q +              &
               &  32*f1A**2*ml**2*mN**4*pk + 32*f1V**2*ml**2*mN**4*pk - 32*f1A**2*ml**2*mN**2*p1p*pk -             &
               &  32*f1V**2*ml**2*mN**2*p1p*pk - 32*f1A**2*ml**2*mN**2*p1q*pk - 32*f1V**2*ml**2*mN**2*p1q*pk -     &
               &  64*f1A**2*mN**4*pk**2 - 64*f1V**2*mN**4*pk**2 + 64*f1A**2*mN**2*p1p*pk**2 +                      &
               &  64*f1V**2*mN**2*p1p*pk**2 + 64*f1A**2*mN**2*p1q*pk**2 + 64*f1V**2*mN**2*p1q*pk**2 +              &
               &  16*f1V*f2V*ml**4*mN**3*pq - 48*f1A**2*ml**2*mN**4*pq + 32*f1A*f1V*ml**2*mN**4*pq +               &
               &  16*f1V**2*ml**2*mN**4*pq + 16*f2V**2*ml**4*mN**4*pq + 64*f1A*f2V*ml**2*mN**5*pq +                &
               &  32*f1A**2*ml**2*mN**2*p1p*pq - 32*f1A*f1V*ml**2*mN**2*p1p*pq - 64*f1A*f2V*ml**2*mN**3*p1p*pq +   &
               &  32*f1A**2*ml**2*mN**2*p1q*pq - 32*f1A*f1V*ml**2*mN**2*p1q*pq - 64*f1A*f2V*ml**2*mN**3*p1q*pq +   &
               &  32*f1A**2*ml**2*mN**2*pk*pq + 32*f1V**2*ml**2*mN**2*pk*pq + 32*f1V*f2V*ml**2*mN**3*pk*pq +       &
               &  64*f1A**2*mN**4*pk*pq + 64*f1V**2*mN**4*pk*pq + 64*f2V**2*ml**2*mN**4*pk*pq -                    &
               &  32*f1V*f2V*ml**2*mN*p1p*pk*pq - 64*f1A**2*mN**2*p1p*pk*pq - 64*f1V**2*mN**2*p1p*pk*pq -          &
               &  64*f2V**2*ml**2*mN**2*p1p*pk*pq - 32*f1V*f2V*ml**2*mN*p1q*pk*pq - 64*f1A**2*mN**2*p1q*pk*pq -    &
               &  64*f1V**2*mN**2*p1q*pk*pq - 64*f2V**2*ml**2*mN**2*p1q*pk*pq - 64*f1A**2*mN**2*pk**2*pq -         &
               &  64*f1V**2*mN**2*pk**2*pq - 24*f1A**2*ml**2*mN**2*pq**2 + 32*f1A*f1V*ml**2*mN**2*pq**2 -          &
               &  8*f1V**2*ml**2*mN**2*pq**2 - 4*f2V**2*ml**4*mN**2*pq**2 + 96*f1A*f2V*ml**2*mN**3*pq**2 -         &
               &  32*f1V*f2V*ml**2*mN**3*pq**2 - 32*f2V**2*ml**2*mN**4*pq**2 + 16*f1V*f2V*ml**2*mN*p1k*pq**2 -     &
               &  4*f2V**2*ml**4*p1p*pq**2 - 32*f1A*f2V*ml**2*mN*p1p*pq**2 + 32*f1V*f2V*ml**2*mN*p1p*pq**2 +       &
               &  32*f2V**2*ml**2*mN**2*p1p*pq**2 - 48*f1A*f2V*ml**2*mN*p1q*pq**2 +                                &
               &  16*f1V*f2V*ml**2*mN*p1q*pq**2 + 32*f2V**2*ml**2*mN**2*p1q*pq**2 +                                &
               &  16*f1V*f2V*ml**2*mN*pk*pq**2 + 64*f1A**2*mN**2*pk*pq**2 + 64*f1V**2*mN**2*pk*pq**2 +             &
               &  64*f2V**2*ml**2*mN**2*pk*pq**2 + 32*f1A**2*p1k*pk*pq**2 + 32*f1V**2*p1k*pk*pq**2 -               &
               &  16*f1A**2*p1q*pk*pq**2 - 32*f1A*f1V*p1q*pk*pq**2 - 16*f1V**2*p1q*pk*pq**2 -                      &
               &  16*f2V**2*ml**2*p1q*pk*pq**2 + 16*f1A*f2V*ml**2*mN*pq**3 - 16*f1V*f2V*ml**2*mN*pq**3 -           &
               &  32*f2V**2*ml**2*mN**2*pq**3 - 16*f1A**2*p1k*pq**3 + 32*f1A*f1V*p1k*pq**3 -                       &
               &  16*f1V**2*p1k*pq**3 - 16*f2V**2*ml**2*p1k*pq**3 + 16*f2V**2*ml**2*p1q*pq**3 -                    &
               &  8*f1V*f2V*ml**4*mN**3*Q2 + 16*f1A**2*ml**2*mN**4*Q2 - 16*f1V**2*ml**2*mN**4*Q2 -                 &
               &  12*f2V**2*ml**4*mN**4*Q2 - 16*f1V*f2V*ml**2*mN**5*Q2 - 16*f1A**2*mN**6*Q2 +                      &
               &  16*f1V**2*mN**6*Q2 - 8*f1A**2*ml**2*mN**2*p1p*Q2 + 8*f1V**2*ml**2*mN**2*p1p*Q2 +                 &
               &  4*f2V**2*ml**4*mN**2*p1p*Q2 + 16*f1V*f2V*ml**2*mN**3*p1p*Q2 + 16*f1A**2*mN**4*p1p*Q2 -           &
               &  16*f1V**2*mN**4*p1p*Q2 - 8*f1A**2*ml**2*mN**2*p1q*Q2 + 8*f1V**2*ml**2*mN**2*p1q*Q2 +             &
               &  4*f2V**2*ml**4*mN**2*p1q*Q2 + 16*f1V*f2V*ml**2*mN**3*p1q*Q2 + 16*f1A**2*mN**4*p1q*Q2 -           &
               &  16*f1V**2*mN**4*p1q*Q2 - 16*f1A**2*ml**2*mN**2*pk*Q2 - 16*f1V**2*ml**2*mN**2*pk*Q2 -             &
               &  16*f1V*f2V*ml**2*mN**3*pk*Q2 - 64*f1A*f1V*mN**4*pk*Q2 - 128*f1A*f2V*mN**5*pk*Q2 +                &
               &  16*f1V*f2V*ml**2*mN*p1p*pk*Q2 + 64*f1A*f1V*mN**2*p1p*pk*Q2 + 128*f1A*f2V*mN**3*p1p*pk*Q2 +       &
               &  16*f1V*f2V*ml**2*mN*p1q*pk*Q2 + 64*f1A*f1V*mN**2*p1q*pk*Q2 + 128*f1A*f2V*mN**3*p1q*pk*Q2 +       &
               &  32*f1A**2*mN**2*pk**2*Q2 + 32*f1V**2*mN**2*pk**2*Q2 - 64*f2V**2*mN**4*pk**2*Q2 +                 &
               &  64*f2V**2*mN**2*p1p*pk**2*Q2 + 64*f2V**2*mN**2*p1q*pk**2*Q2 + 16*f1A**2*ml**2*mN**2*pq*Q2 -      &
               &  16*f1A*f1V*ml**2*mN**2*pq*Q2 - 48*f1A*f2V*ml**2*mN**3*pq*Q2 - 48*f1A**2*mN**4*pq*Q2 +            &
               &  32*f1A*f1V*mN**4*pq*Q2 + 16*f1V**2*mN**4*pq*Q2 + 64*f1A*f2V*mN**5*pq*Q2 -                        &
               &  16*f1V*f2V*ml**2*mN*p1k*pq*Q2 + 4*f2V**2*ml**4*p1p*pq*Q2 + 16*f1A*f2V*ml**2*mN*p1p*pq*Q2 -       &
               &  16*f1V*f2V*ml**2*mN*p1p*pq*Q2 + 32*f1A**2*mN**2*p1p*pq*Q2 - 32*f1A*f1V*mN**2*p1p*pq*Q2 -         &
               &  64*f1A*f2V*mN**3*p1p*pq*Q2 + 32*f1A*f2V*ml**2*mN*p1q*pq*Q2 + 32*f1A**2*mN**2*p1q*pq*Q2 -         &
               &  32*f1A*f1V*mN**2*p1q*pq*Q2 - 64*f1A*f2V*mN**3*p1q*pq*Q2 - 16*f1V*f2V*ml**2*mN*pk*pq*Q2 -         &
               &  32*f1A**2*mN**2*pk*pq*Q2 - 64*f1A*f1V*mN**2*pk*pq*Q2 - 32*f1V**2*mN**2*pk*pq*Q2 -                &
               &  32*f2V**2*ml**2*mN**2*pk*pq*Q2 - 192*f1A*f2V*mN**3*pk*pq*Q2 + 64*f2V**2*mN**4*pk*pq*Q2 -         &
               &  32*f1A**2*p1k*pk*pq*Q2 - 32*f1V**2*p1k*pk*pq*Q2 + 64*f1A*f2V*mN*p1p*pk*pq*Q2 -                   &
               &  64*f2V**2*mN**2*p1p*pk*pq*Q2 + 16*f1A**2*p1q*pk*pq*Q2 + 32*f1A*f1V*p1q*pk*pq*Q2 +                &
               &  16*f1V**2*p1q*pk*pq*Q2 + 16*f2V**2*ml**2*p1q*pk*pq*Q2 + 64*f1A*f2V*mN*p1q*pk*pq*Q2 -             &
               &  64*f2V**2*mN**2*p1q*pk*pq*Q2 - 64*f2V**2*mN**2*pk**2*pq*Q2 - 16*f1A*f2V*ml**2*mN*pq**2*Q2 +      &
               &  16*f1V*f2V*ml**2*mN*pq**2*Q2 - 24*f1A**2*mN**2*pq**2*Q2 + 32*f1A*f1V*mN**2*pq**2*Q2 -            &
               &  8*f1V**2*mN**2*pq**2*Q2 + 20*f2V**2*ml**2*mN**2*pq**2*Q2 + 96*f1A*f2V*mN**3*pq**2*Q2 -           &
               &  32*f1V*f2V*mN**3*pq**2*Q2 - 32*f2V**2*mN**4*pq**2*Q2 + 16*f1A**2*p1k*pq**2*Q2 -                  &
               &  32*f1A*f1V*p1k*pq**2*Q2 + 16*f1V**2*p1k*pq**2*Q2 + 16*f2V**2*ml**2*p1k*pq**2*Q2 +                &
               &  32*f1A*f2V*mN*p1k*pq**2*Q2 - 4*f2V**2*ml**2*p1p*pq**2*Q2 - 32*f1A*f2V*mN*p1p*pq**2*Q2 +          &
               &  32*f1V*f2V*mN*p1p*pq**2*Q2 + 32*f2V**2*mN**2*p1p*pq**2*Q2 - 16*f2V**2*ml**2*p1q*pq**2*Q2 -       &
               &  48*f1A*f2V*mN*p1q*pq**2*Q2 + 16*f1V*f2V*mN*p1q*pq**2*Q2 + 32*f2V**2*mN**2*p1q*pq**2*Q2 -         &
               &  32*f1A*f2V*mN*pk*pq**2*Q2 + 64*f2V**2*mN**2*pk*pq**2*Q2 + 32*f2V**2*p1k*pk*pq**2*Q2 -            &
               &  16*f2V**2*p1q*pk*pq**2*Q2 + 16*f1A*f2V*mN*pq**3*Q2 - 16*f1V*f2V*mN*pq**3*Q2 -                    &
               &  32*f2V**2*mN**2*pq**3*Q2 - 16*f2V**2*p1k*pq**3*Q2 + 16*f2V**2*p1q*pq**3*Q2 -                     &
               &  2*f1A**2*ml**2*mN**2*Q2**2 + 2*f1V**2*ml**2*mN**2*Q2**2 + f2V**2*ml**4*mN**2*Q2**2 +             &
               &  8*f1V*f2V*ml**2*mN**3*Q2**2 + 16*f1A**2*mN**4*Q2**2 - 16*f1V**2*mN**4*Q2**2 +                    &
               &  4*f2V**2*ml**2*mN**4*Q2**2 - 32*f1V*f2V*mN**5*Q2**2 - 16*f2V**2*mN**6*Q2**2 +                    &
               &  4*f1V*f2V*ml**2*mN*p1k*Q2**2 - f2V**2*ml**4*p1p*Q2**2 - 8*f1A**2*mN**2*p1p*Q2**2 +               &
               &  8*f1V**2*mN**2*p1p*Q2**2 - 4*f2V**2*ml**2*mN**2*p1p*Q2**2 + 32*f1V*f2V*mN**3*p1p*Q2**2 +         &
               &  16*f2V**2*mN**4*p1p*Q2**2 - 4*f1A*f2V*ml**2*mN*p1q*Q2**2 - 4*f1V*f2V*ml**2*mN*p1q*Q2**2 -        &
               &  8*f1A**2*mN**2*p1q*Q2**2 + 8*f1V**2*mN**2*p1q*Q2**2 - 4*f2V**2*ml**2*mN**2*p1q*Q2**2 +           &
               &  32*f1V*f2V*mN**3*p1q*Q2**2 + 16*f2V**2*mN**4*p1q*Q2**2 + 4*f1V*f2V*ml**2*mN*pk*Q2**2 +           &
               &  32*f1A*f1V*mN**2*pk*Q2**2 + 96*f1A*f2V*mN**3*pk*Q2**2 + 8*f1A**2*p1k*pk*Q2**2 +                  &
               &  8*f1V**2*p1k*pk*Q2**2 - 32*f1A*f2V*mN*p1p*pk*Q2**2 - 4*f1A**2*p1q*pk*Q2**2 -                     &
               &  8*f1A*f1V*p1q*pk*Q2**2 - 4*f1V**2*p1q*pk*Q2**2 - 4*f2V**2*ml**2*p1q*pk*Q2**2 -                   &
               &  32*f1A*f2V*mN*p1q*pk*Q2**2 + 32*f2V**2*mN**2*pk**2*Q2**2 + 4*f1A*f2V*ml**2*mN*pq*Q2**2 -         &
               &  4*f1V*f2V*ml**2*mN*pq*Q2**2 + 16*f1A**2*mN**2*pq*Q2**2 - 16*f1A*f1V*mN**2*pq*Q2**2 -             &
               &  48*f1A*f2V*mN**3*pq*Q2**2 - 16*f1V*f2V*mN**3*pq*Q2**2 - 16*f2V**2*mN**4*pq*Q2**2 -               &
               &  4*f1A**2*p1k*pq*Q2**2 + 8*f1A*f1V*p1k*pq*Q2**2 - 4*f1V**2*p1k*pq*Q2**2 -                         &
               &  4*f2V**2*ml**2*p1k*pq*Q2**2 - 32*f1A*f2V*mN*p1k*pq*Q2**2 + 4*f2V**2*ml**2*p1p*pq*Q2**2 +         &
               &  16*f1A*f2V*mN*p1p*pq*Q2**2 - 16*f1V*f2V*mN*p1p*pq*Q2**2 + 4*f2V**2*ml**2*p1q*pq*Q2**2 +          &
               &  32*f1A*f2V*mN*p1q*pq*Q2**2 + 32*f1A*f2V*mN*pk*pq*Q2**2 - 32*f2V**2*mN**2*pk*pq*Q2**2 -           &
               &  32*f2V**2*p1k*pk*pq*Q2**2 + 16*f2V**2*p1q*pk*pq*Q2**2 - 16*f1A*f2V*mN*pq**2*Q2**2 +              &
               &  16*f1V*f2V*mN*pq**2*Q2**2 + 24*f2V**2*mN**2*pq**2*Q2**2 + 16*f2V**2*p1k*pq**2*Q2**2 -            &
               &  16*f2V**2*p1q*pq**2*Q2**2 - 2*f1A**2*mN**2*Q2**3 + 2*f1V**2*mN**2*Q2**3 -                        &
               &  f2V**2*ml**2*mN**2*Q2**3 + 16*f1V*f2V*mN**3*Q2**3 + 16*f2V**2*mN**4*Q2**3 +                      &
               &  8*f1A*f2V*mN*p1k*Q2**3 - f2V**2*ml**2*p1p*Q2**3 - 8*f2V**2*mN**2*p1p*Q2**3 -                     &
               &  4*f1A*f2V*mN*p1q*Q2**3 - 4*f1V*f2V*mN*p1q*Q2**3 - 8*f2V**2*mN**2*p1q*Q2**3 -                     &
               &  8*f1A*f2V*mN*pk*Q2**3 + 8*f2V**2*p1k*pk*Q2**3 - 4*f2V**2*p1q*pk*Q2**3 + 4*f1A*f2V*mN*pq*Q2**3 -  &
               &  4*f1V*f2V*mN*pq*Q2**3 - 4*f2V**2*p1k*pq*Q2**3 + 4*f2V**2*p1q*pq*Q2**3 -                          &
               &  2*f2V**2*mN**2*Q2**4)* g_A**2/4./f_pi**2 / ( abs4Sq(p+q) - mN**2 )**2


          ! crossed nucleon pole
          M2cNp = 16*(16*f1V*f2V*ml**4*mN**5 - 16*f1A**2*ml**2*mN**6 + 16*f1V**2*ml**2*mN**6 +                   &
               &    16*f2V**2*ml**4*mN**6 - 32*f1A**2*ml**2*mN**4*p1k - 32*f1V**2*ml**2*mN**4*p1k -                  &
               &    64*f1A**2*mN**4*p1k**2 - 64*f1V**2*mN**4*p1k**2 - 16*f1V*f2V*ml**4*mN**3*p1p +                   &
               &    16*f1A**2*ml**2*mN**4*p1p - 16*f1V**2*ml**2*mN**4*p1p - 16*f2V**2*ml**4*mN**4*p1p +              &
               &    32*f1A**2*ml**2*mN**2*p1k*p1p + 32*f1V**2*ml**2*mN**2*p1k*p1p + 64*f1A**2*mN**2*p1k**2*p1p +     &
               &    64*f1V**2*mN**2*p1k**2*p1p - 16*f1V*f2V*ml**4*mN**3*p1q + 48*f1A**2*ml**2*mN**4*p1q +            &
               &    32*f1A*f1V*ml**2*mN**4*p1q - 16*f1V**2*ml**2*mN**4*p1q - 16*f2V**2*ml**4*mN**4*p1q +             &
               &    64*f1A*f2V*ml**2*mN**5*p1q + 32*f1A**2*ml**2*mN**2*p1k*p1q + 32*f1V**2*ml**2*mN**2*p1k*p1q +     &
               &    32*f1V*f2V*ml**2*mN**3*p1k*p1q + 64*f1A**2*mN**4*p1k*p1q + 64*f1V**2*mN**4*p1k*p1q +             &
               &    64*f2V**2*ml**2*mN**4*p1k*p1q + 64*f1A**2*mN**2*p1k**2*p1q + 64*f1V**2*mN**2*p1k**2*p1q -        &
               &    32*f1A**2*ml**2*mN**2*p1p*p1q - 32*f1A*f1V*ml**2*mN**2*p1p*p1q -                                 &
               &    64*f1A*f2V*ml**2*mN**3*p1p*p1q - 32*f1V*f2V*ml**2*mN*p1k*p1p*p1q -                               &
               &    64*f1A**2*mN**2*p1k*p1p*p1q - 64*f1V**2*mN**2*p1k*p1p*p1q -                                      &
               &    64*f2V**2*ml**2*mN**2*p1k*p1p*p1q - 24*f1A**2*ml**2*mN**2*p1q**2 -                               &
               &    32*f1A*f1V*ml**2*mN**2*p1q**2 - 8*f1V**2*ml**2*mN**2*p1q**2 - 4*f2V**2*ml**4*mN**2*p1q**2 -      &
               &    96*f1A*f2V*ml**2*mN**3*p1q**2 - 32*f1V*f2V*ml**2*mN**3*p1q**2 - 32*f2V**2*ml**2*mN**4*p1q**2 -   &
               &    16*f1V*f2V*ml**2*mN*p1k*p1q**2 - 64*f1A**2*mN**2*p1k*p1q**2 - 64*f1V**2*mN**2*p1k*p1q**2 -       &
               &    64*f2V**2*ml**2*mN**2*p1k*p1q**2 - 4*f2V**2*ml**4*p1p*p1q**2 +                                   &
               &    32*f1A*f2V*ml**2*mN*p1p*p1q**2 + 32*f1V*f2V*ml**2*mN*p1p*p1q**2 +                                &
               &    32*f2V**2*ml**2*mN**2*p1p*p1q**2 + 16*f1A*f2V*ml**2*mN*p1q**3 + 16*f1V*f2V*ml**2*mN*p1q**3 +     &
               &    32*f2V**2*ml**2*mN**2*p1q**3 - 16*f1V*f2V*ml**2*mN*p1q**2*pk + 32*f1A**2*p1k*p1q**2*pk +         &
               &    32*f1V**2*p1k*p1q**2*pk - 16*f1A**2*p1q**3*pk - 32*f1A*f1V*p1q**3*pk - 16*f1V**2*p1q**3*pk -     &
               &    16*f2V**2*ml**2*p1q**3*pk + 16*f1V*f2V*ml**4*mN**3*pq - 16*f1A**2*ml**2*mN**4*pq +               &
               &    16*f1V**2*ml**2*mN**4*pq + 16*f2V**2*ml**4*mN**4*pq - 32*f1A**2*ml**2*mN**2*p1k*pq -             &
               &    32*f1V**2*ml**2*mN**2*p1k*pq - 64*f1A**2*mN**2*p1k**2*pq - 64*f1V**2*mN**2*p1k**2*pq +           &
               &    32*f1A**2*ml**2*mN**2*p1q*pq + 32*f1A*f1V*ml**2*mN**2*p1q*pq + 64*f1A*f2V*ml**2*mN**3*p1q*pq +   &
               &    32*f1V*f2V*ml**2*mN*p1k*p1q*pq + 64*f1A**2*mN**2*p1k*p1q*pq + 64*f1V**2*mN**2*p1k*p1q*pq +       &
               &    64*f2V**2*ml**2*mN**2*p1k*p1q*pq - 48*f1A*f2V*ml**2*mN*p1q**2*pq -                               &
               &    16*f1V*f2V*ml**2*mN*p1q**2*pq - 32*f2V**2*ml**2*mN**2*p1q**2*pq - 16*f1A**2*p1k*p1q**2*pq +      &
               &    32*f1A*f1V*p1k*p1q**2*pq - 16*f1V**2*p1k*p1q**2*pq - 16*f2V**2*ml**2*p1k*p1q**2*pq +             &
               &    16*f2V**2*ml**2*p1q**3*pq - 8*f1V*f2V*ml**4*mN**3*Q2 + 16*f1A**2*ml**2*mN**4*Q2 -                &
               &    16*f1V**2*ml**2*mN**4*Q2 - 12*f2V**2*ml**4*mN**4*Q2 - 16*f1V*f2V*ml**2*mN**5*Q2 -                &
               &    16*f1A**2*mN**6*Q2 + 16*f1V**2*mN**6*Q2 + 16*f1A**2*ml**2*mN**2*p1k*Q2 +                         &
               &    16*f1V**2*ml**2*mN**2*p1k*Q2 + 16*f1V*f2V*ml**2*mN**3*p1k*Q2 - 64*f1A*f1V*mN**4*p1k*Q2 -         &
               &    128*f1A*f2V*mN**5*p1k*Q2 + 32*f1A**2*mN**2*p1k**2*Q2 + 32*f1V**2*mN**2*p1k**2*Q2 -               &
               &    64*f2V**2*mN**4*p1k**2*Q2 - 8*f1A**2*ml**2*mN**2*p1p*Q2 + 8*f1V**2*ml**2*mN**2*p1p*Q2 +          &
               &    4*f2V**2*ml**4*mN**2*p1p*Q2 + 16*f1V*f2V*ml**2*mN**3*p1p*Q2 + 16*f1A**2*mN**4*p1p*Q2 -           &
               &    16*f1V**2*mN**4*p1p*Q2 - 16*f1V*f2V*ml**2*mN*p1k*p1p*Q2 + 64*f1A*f1V*mN**2*p1k*p1p*Q2 +          &
               &    128*f1A*f2V*mN**3*p1k*p1p*Q2 + 64*f2V**2*mN**2*p1k**2*p1p*Q2 - 16*f1A**2*ml**2*mN**2*p1q*Q2 -    &
               &    16*f1A*f1V*ml**2*mN**2*p1q*Q2 - 48*f1A*f2V*ml**2*mN**3*p1q*Q2 + 48*f1A**2*mN**4*p1q*Q2 +         &
               &    32*f1A*f1V*mN**4*p1q*Q2 - 16*f1V**2*mN**4*p1q*Q2 + 64*f1A*f2V*mN**5*p1q*Q2 -                     &
               &    16*f1V*f2V*ml**2*mN*p1k*p1q*Q2 - 32*f1A**2*mN**2*p1k*p1q*Q2 + 64*f1A*f1V*mN**2*p1k*p1q*Q2 -      &
               &    32*f1V**2*mN**2*p1k*p1q*Q2 - 32*f2V**2*ml**2*mN**2*p1k*p1q*Q2 + 192*f1A*f2V*mN**3*p1k*p1q*Q2 +   &
               &    64*f2V**2*mN**4*p1k*p1q*Q2 + 64*f2V**2*mN**2*p1k**2*p1q*Q2 - 4*f2V**2*ml**4*p1p*p1q*Q2 +         &
               &    16*f1A*f2V*ml**2*mN*p1p*p1q*Q2 + 16*f1V*f2V*ml**2*mN*p1p*p1q*Q2 - 32*f1A**2*mN**2*p1p*p1q*Q2 -   &
               &    32*f1A*f1V*mN**2*p1p*p1q*Q2 - 64*f1A*f2V*mN**3*p1p*p1q*Q2 - 64*f1A*f2V*mN*p1k*p1p*p1q*Q2 -       &
               &    64*f2V**2*mN**2*p1k*p1p*p1q*Q2 + 16*f1A*f2V*ml**2*mN*p1q**2*Q2 +                                 &
               &    16*f1V*f2V*ml**2*mN*p1q**2*Q2 - 24*f1A**2*mN**2*p1q**2*Q2 - 32*f1A*f1V*mN**2*p1q**2*Q2 -         &
               &    8*f1V**2*mN**2*p1q**2*Q2 + 20*f2V**2*ml**2*mN**2*p1q**2*Q2 - 96*f1A*f2V*mN**3*p1q**2*Q2 -        &
               &    32*f1V*f2V*mN**3*p1q**2*Q2 - 32*f2V**2*mN**4*p1q**2*Q2 - 32*f1A*f2V*mN*p1k*p1q**2*Q2 -           &
               &    64*f2V**2*mN**2*p1k*p1q**2*Q2 - 4*f2V**2*ml**2*p1p*p1q**2*Q2 + 32*f1A*f2V*mN*p1p*p1q**2*Q2 +     &
               &    32*f1V*f2V*mN*p1p*p1q**2*Q2 + 32*f2V**2*mN**2*p1p*p1q**2*Q2 + 16*f1A*f2V*mN*p1q**3*Q2 +          &
               &    16*f1V*f2V*mN*p1q**3*Q2 + 32*f2V**2*mN**2*p1q**3*Q2 - 16*f1V*f2V*ml**2*mN*p1q*pk*Q2 +            &
               &    32*f1A**2*p1k*p1q*pk*Q2 + 32*f1V**2*p1k*p1q*pk*Q2 - 16*f1A**2*p1q**2*pk*Q2 -                     &
               &    32*f1A*f1V*p1q**2*pk*Q2 - 16*f1V**2*p1q**2*pk*Q2 - 16*f2V**2*ml**2*p1q**2*pk*Q2 +                &
               &    32*f1A*f2V*mN*p1q**2*pk*Q2 + 32*f2V**2*p1k*p1q**2*pk*Q2 - 16*f2V**2*p1q**3*pk*Q2 +               &
               &    8*f1A**2*ml**2*mN**2*pq*Q2 - 8*f1V**2*ml**2*mN**2*pq*Q2 - 4*f2V**2*ml**4*mN**2*pq*Q2 -           &
               &    16*f1V*f2V*ml**2*mN**3*pq*Q2 - 16*f1A**2*mN**4*pq*Q2 + 16*f1V**2*mN**4*pq*Q2 +                   &
               &    16*f1V*f2V*ml**2*mN*p1k*pq*Q2 - 64*f1A*f1V*mN**2*p1k*pq*Q2 - 128*f1A*f2V*mN**3*p1k*pq*Q2 -       &
               &    64*f2V**2*mN**2*p1k**2*pq*Q2 - 32*f1A*f2V*ml**2*mN*p1q*pq*Q2 + 32*f1A**2*mN**2*p1q*pq*Q2 +       &
               &    32*f1A*f1V*mN**2*p1q*pq*Q2 + 64*f1A*f2V*mN**3*p1q*pq*Q2 - 16*f1A**2*p1k*p1q*pq*Q2 +              &
               &    32*f1A*f1V*p1k*p1q*pq*Q2 - 16*f1V**2*p1k*p1q*pq*Q2 - 16*f2V**2*ml**2*p1k*p1q*pq*Q2 +             &
               &    64*f1A*f2V*mN*p1k*p1q*pq*Q2 + 64*f2V**2*mN**2*p1k*p1q*pq*Q2 + 16*f2V**2*ml**2*p1q**2*pq*Q2 -     &
               &    48*f1A*f2V*mN*p1q**2*pq*Q2 - 16*f1V*f2V*mN*p1q**2*pq*Q2 - 32*f2V**2*mN**2*p1q**2*pq*Q2 -         &
               &    16*f2V**2*p1k*p1q**2*pq*Q2 + 16*f2V**2*p1q**3*pq*Q2 - 2*f1A**2*ml**2*mN**2*Q2**2 +               &
               &    2*f1V**2*ml**2*mN**2*Q2**2 + f2V**2*ml**4*mN**2*Q2**2 + 8*f1V*f2V*ml**2*mN**3*Q2**2 +            &
               &    16*f1A**2*mN**4*Q2**2 - 16*f1V**2*mN**4*Q2**2 + 4*f2V**2*ml**2*mN**4*Q2**2 -                     &
               &    32*f1V*f2V*mN**5*Q2**2 - 16*f2V**2*mN**6*Q2**2 - 4*f1V*f2V*ml**2*mN*p1k*Q2**2 +                  &
               &    32*f1A*f1V*mN**2*p1k*Q2**2 + 96*f1A*f2V*mN**3*p1k*Q2**2 + 32*f2V**2*mN**2*p1k**2*Q2**2 -         &
               &    f2V**2*ml**4*p1p*Q2**2 - 8*f1A**2*mN**2*p1p*Q2**2 + 8*f1V**2*mN**2*p1p*Q2**2 -                   &
               &    4*f2V**2*ml**2*mN**2*p1p*Q2**2 + 32*f1V*f2V*mN**3*p1p*Q2**2 + 16*f2V**2*mN**4*p1p*Q2**2 -        &
               &    32*f1A*f2V*mN*p1k*p1p*Q2**2 + 4*f1A*f2V*ml**2*mN*p1q*Q2**2 + 4*f1V*f2V*ml**2*mN*p1q*Q2**2 -      &
               &    16*f1A**2*mN**2*p1q*Q2**2 - 16*f1A*f1V*mN**2*p1q*Q2**2 - 48*f1A*f2V*mN**3*p1q*Q2**2 +            &
               &    16*f1V*f2V*mN**3*p1q*Q2**2 + 16*f2V**2*mN**4*p1q*Q2**2 - 32*f1A*f2V*mN*p1k*p1q*Q2**2 -           &
               &    32*f2V**2*mN**2*p1k*p1q*Q2**2 - 4*f2V**2*ml**2*p1p*p1q*Q2**2 + 16*f1A*f2V*mN*p1p*p1q*Q2**2 +     &
               &    16*f1V*f2V*mN*p1p*p1q*Q2**2 + 16*f1A*f2V*mN*p1q**2*Q2**2 + 16*f1V*f2V*mN*p1q**2*Q2**2 +          &
               &    24*f2V**2*mN**2*p1q**2*Q2**2 - 4*f1V*f2V*ml**2*mN*pk*Q2**2 + 8*f1A**2*p1k*pk*Q2**2 +             &
               &    8*f1V**2*p1k*pk*Q2**2 - 4*f1A**2*p1q*pk*Q2**2 - 8*f1A*f1V*p1q*pk*Q2**2 -                         &
               &    4*f1V**2*p1q*pk*Q2**2 - 4*f2V**2*ml**2*p1q*pk*Q2**2 + 32*f1A*f2V*mN*p1q*pk*Q2**2 +               &
               &    32*f2V**2*p1k*p1q*pk*Q2**2 - 16*f2V**2*p1q**2*pk*Q2**2 - 4*f1A*f2V*ml**2*mN*pq*Q2**2 +           &
               &    4*f1V*f2V*ml**2*mN*pq*Q2**2 + 8*f1A**2*mN**2*pq*Q2**2 - 8*f1V**2*mN**2*pq*Q2**2 +                &
               &    4*f2V**2*ml**2*mN**2*pq*Q2**2 - 32*f1V*f2V*mN**3*pq*Q2**2 - 16*f2V**2*mN**4*pq*Q2**2 -           &
               &    4*f1A**2*p1k*pq*Q2**2 + 8*f1A*f1V*p1k*pq*Q2**2 - 4*f1V**2*p1k*pq*Q2**2 -                         &
               &    4*f2V**2*ml**2*p1k*pq*Q2**2 + 32*f1A*f2V*mN*p1k*pq*Q2**2 + 4*f2V**2*ml**2*p1q*pq*Q2**2 -         &
               &    32*f1A*f2V*mN*p1q*pq*Q2**2 - 16*f2V**2*p1k*p1q*pq*Q2**2 + 16*f2V**2*p1q**2*pq*Q2**2 -            &
               &    2*f1A**2*mN**2*Q2**3 + 2*f1V**2*mN**2*Q2**3 - f2V**2*ml**2*mN**2*Q2**3 +                         &
               &    16*f1V*f2V*mN**3*Q2**3 + 16*f2V**2*mN**4*Q2**3 - 8*f1A*f2V*mN*p1k*Q2**3 -                        &
               &    f2V**2*ml**2*p1p*Q2**3 - 8*f2V**2*mN**2*p1p*Q2**3 + 4*f1A*f2V*mN*p1q*Q2**3 +                     &
               &    4*f1V*f2V*mN*p1q*Q2**3 + 8*f1A*f2V*mN*pk*Q2**3 + 8*f2V**2*p1k*pk*Q2**3 -                         &
               &    4*f2V**2*p1q*pk*Q2**3 - 4*f1A*f2V*mN*pq*Q2**3 + 4*f1V*f2V*mN*pq*Q2**3 +                          &
               &    8*f2V**2*mN**2*pq*Q2**3 - 4*f2V**2*p1k*pq*Q2**3 + 4*f2V**2*p1q*pq*Q2**3 -                        &
               &    2*f2V**2*mN**2*Q2**4) * g_A**2/4./f_pi**2 / ( abs4Sq(p_out-q) - mN**2 )**2



          M2CT =     16.*( F_rho**2+g_A**2*FV_CT**2)/f_pi**2*( 4.*pk*p1k - 2.*p1k*pq - 2.*pk*p1q  ) &
               & + 16.*(-F_rho**2+g_A**2*FV_CT**2)/f_pi**2*(Q2+ml**2)*mN*mN1 &
               & + 64.*g_A*FV_CT*F_rho/f_pi**2 * ( p1k*pq - pk*p1q )


          M2pp=8.*F_rho**2*ml**2*(Q2+ml**2)/f_pi**2/(Q2+mpi**2)**2 *( 2.*pq*p1q + Q2*p1p - Q2*mN*mN1  )


          M2pF=32.*F_PF**2*g_A**2*(mN+mN1)**2*(mN*mN1-p1p)/f_pi**2/(mN**2+mN1**2-mpi**2-2.*p1p)**2*( &
               & Q2*(mN**2-2.*p1p+mN1**2) - 4.*(pk-p1k)**2 + 4.*( pk*pq - p1k*pq - pk*p1q + p1k*p1q) &
               & + ml**2*( mN**2 +mN1**2 -Q2/4. +2.*(pk-p1k-p1p) ) -ml**4/4. )



          kinemat_factor=k1*ppi**2/(16.*mN*Enu)/(2.*pi)**4/abs(2.*vecW*ct_Wpi*Epi - 2.*W(0)*ppi);

          if (process_ID.eq.2) then  !CC
             M2cDp=M2cDp/3.   ! Clebsch-Goedon for the p pi+   final state  for nu-p  reaction
             M2Np=M2Np*3.     ! Clebsch-Goedon for the sum of (n pi+   + p pi^0) final states for nu-n reaction
             M2cNp=M2cNp*2.   ! Clebsch-Goedon for p pi+ final state for nu-p  reaction
             ! CT, pp, pF coefficients for p pi+ final states are 1.
             kinemat_factor=kinemat_factor*GF**2*coscab**2/2*3.8938e-28     ! 1/GeV^2   to  cm^2;     for  CC
          else if (process_ID.eq.1)then  !EM
             !M2cDp=M2cDp*1   ! Clebsch-Goedon for the sum of (n pi+   + p pi^0) final states  for e-p  reaction
             M2Np=M2Np*3.     ! Clebsch-Goedon for the sum of (n pi+   + p pi^0) final states for e-p reaction
             M2cNp=M2cNp*3.   ! Clebsch-Goedon for the sum of (n pi+   + p pi^0)  for e-p  reaction
             M2CT=M2CT*3.
             !M2pp=0
             M2pF=M2pF*3.
             kinemat_factor=kinemat_factor *(2*pi)**2*alphaQED**2/Q2**2 *3.8938e-28     ! 1/GeV^2   to  cm^2;     for  CC
          else
             stop 'Wrong process number. It should be either 1 (EM)  or 2 (CC)'
          end if



          !  write(10,'(3(A,g12.5))') 'In HNV_free_elepton_ct_ctPi_phiPi  kinemat_factor= ', kinemat_factor, ' k1= ',k1, '  ppi=',ppi

          if ( (M2cDp*kinemat_factor) .le. 0) xsec_cDp=0.
          if ( (M2Np*kinemat_factor) .le. 0) xsec_Np=0.
          if ( (M2cNp*kinemat_factor) .le. 0) xsec_cNp=0.
          if ( (M2CT*kinemat_factor) .le. 0) xsec_CT=0.
          if ( (M2pp*kinemat_factor) .le. 0) xsec_pp=0.
          if ( (M2pF*kinemat_factor) .le. 0) xsec_pF=0.


          xsec_cDp = xsec_cDp + M2cDp/2.*kinemat_factor ! /2. is averaging over the initial nucleon spins

          xsec_Np = xsec_Np + M2Np/2.*kinemat_factor

          xsec_cNp = xsec_cNp + M2cNp/2.*kinemat_factor

          xsec_CT= xsec_CT + M2CT/2.*kinemat_factor

          xsec_pp= xsec_pp + M2pp/2.*kinemat_factor

          xsec_pF= xsec_pF + M2pF/2.*kinemat_factor


          !  write(10,*) 'In HNV_free_elepton_ct_ctPi_phiPi  xsec_CT= ', xsec_CT, '   xsec_pp= ', xsec_pp, '   xsec_pF= ', xsec_pF

       end if  !  ( Epir(r) < Epi)
    end do




    ! ********************************************
    !this is how Epi is determined in Oliver's procedure
    !
    if(debug) then
       write (*,*)  ' '
       write (*,'(4(A,g12.5))')  'E1= ', E1, '   cos(theta_Pi)=',ctPi, '   sin(theta_Pi)=',stPi, '   phi_pi=',phiPi
       write (*,*)  ' '
       write (*,'(A,4(g12.5))')  'k_in = ',k_in
       write (*,'(A,4(g12.5))')  'k_out= ',k_out
       write (*,'(A,4(g12.5))')  'q    = ',q
       write (*,'(A,4(g12.5))')  'p    = ',p
       write (*,'(A,4(g12.5))')  'W    = ',W

       write(*,'(2(A,g12.5))') 'check  vecW*ct_Wpi*ppi=', (vecW*ct_Wpi*ppi), &
            & '  should be equal to  ppi_i*W^i=' , (  Dot_product(ppi_out(1:3),W(1:3))  )



       betaTOCM(1:3)=q(1:3)/(q(0)+mN)

       ! ******************************
       ! set up incoming nucleon
       ! ******************************
       call setToDefault(nucleon_in)
       nucleon_in%position=(/ 1000., 1000., 1000. /)
       nucleon_in%momentum=p
       nucleon_in%mass=baryon(nucleon)%mass
       !nucleon_in%momentum(0)=freeEnergy(nucleon_in)! E=sqrt(p(1:3)^2+m_0^2)
       nucleon_in%ID=nucleon
       nucleon_in%charge =1
       nucleon_in%antiparticle=.false.
       nucleon_in%perturbative=.true.

       call get_k_abs_improved(NumRoots,unit_ppi_out(1),unit_ppi_out(2),unit_ppi_out(3),q,&
            & nucleon_in,1,1,success,p_pion_out_root,betaTOCM) ! 1 und 1 are charges of pion and outgoing nucleon

       if (success) then
          do r1=1,numRoots

             p_pion_abs2 = Dot_Product(p_pion_out_root(r1,1:3),p_pion_out_root(r1,1:3))

             write (*,*)  ' '
             write (*,'(2(A,I6))')  ' Oliver: number of roots: ',  numRoots, '     r1=',r1
             write (*,*)  ' '
             write (*,'(A,4(g12.5))')  'W    = ',W
             write (*,'(A,4(g12.5),A,g12.5)')  'ppi  = ',p_pion_out_root(r,:),&
                  &'     ppi_out^2=',SP(p_pion_out_root(r,:),p_pion_out_root(r,:))
             write (*,'(A,4(g12.5),A,g12.5)')  'p1   = ',(W-p_pion_out_root(r,:)),&
                  &   '     p1^2 = ',SP(W-p_pion_out_root(r,:),W-p_pion_out_root(r,:))
             write (*,*)  ' '
             write (*,*)  ' '
             write (*,*)  ' '

          end do
       end if  ! if success

    end if  ! if debug

    ! ********************************************



  end subroutine HNV_free_elepton_ct_ctPi_phiPi












  ! 3-fold differential xsec dsi / dE1 dcostheta dcosthetaPi
  subroutine HNV_free_elepton_ct_ctPi(process_ID,Enu,E1,ct,ctPi,mN,ml,mN1,mpi,xsec_cDp,xsec_Np,xsec_cNp,xsec_CT,xsec_pp,xsec_pF)

    use gauss_integration
    use constants, only : pi

    implicit none

    integer, intent(in) :: process_ID
    real, intent(in) :: Enu, E1, ct, ctPi, mN, ml, mN1, mpi
    real, intent(out) :: xsec_cDp, xsec_Np, xsec_cNp, xsec_CT,xsec_pp,xsec_pF



    integer :: intprec=1, n2, l
    real, dimension(:), allocatable :: yy_cDp, yy_Np, yy_cNp, yy_CT, yy_pp,yy_pF , xx
    real :: phiPi


    allocate (xx(20*intprec))
    allocate (yy_cDp(20*intprec))
    allocate (yy_Np(20*intprec))
    allocate (yy_cNp(20*intprec))
    allocate (yy_CT(20*intprec))
    allocate (yy_pp(20*intprec))
    allocate (yy_pF(20*intprec))


    call sg20r(0.,2.*pi,intprec,xx,n2) ! n2 - OUT

    do l=1,n2
       phiPi=xx(l)
       call  HNV_free_elepton_ct_ctPi_phiPi(process_ID,Enu,E1,ct,ctPi,phiPi,mN,ml,mN1,mpi, &
            &  yy_cDp(l),yy_Np(l),yy_cNp(l),yy_CT(l),yy_pp(l),yy_pF(l))

       ! write(*,'5(A,g14.5)') 'yy_Np(l)=', yy_Np(l), 'yy_cNp(l)=', yy_cNp(l), 'yy_CT(l)=', yy_CT(l), 'yy_pp(l)=', yy_pp(l),   'yy_pF(l)=', yy_pF(l),
    end do

    call rg20r(0.,2.*pi,intprec,yy_cDp,xsec_cDp)
    call rg20r(0.,2.*pi,intprec,yy_Np,xsec_Np)
    call rg20r(0.,2.*pi,intprec,yy_cNp,xsec_cNp)
    call rg20r(0.,2.*pi,intprec,yy_CT,xsec_CT)
    call rg20r(0.,2.*pi,intprec,yy_pp,xsec_pp)
    call rg20r(0.,2.*pi,intprec,yy_pF,xsec_pF)

    deallocate(xx, yy_cDp, yy_Np, yy_cNp, yy_CT, yy_pp, yy_pF)


  end subroutine HNV_free_elepton_ct_ctPi






  ! 2-fold differential xsec dsi / dE1 dcostheta
  subroutine HNV_free_elepton_ct(process_ID,Enu,E1,ct,mN,ml,mN1,mpi, xsec_cDp,xsec_Np,xsec_cNp,xsec_CT,xsec_pp,xsec_pF)

    use gauss_integration
    use constants, only : pi

    implicit none
    integer, intent(in) :: process_ID
    real, intent(in) :: Enu, E1, ct, mN, ml, mN1, mpi
    real, intent(out) :: xsec_cDp, xsec_Np, xsec_cNp, xsec_CT, xsec_pp, xsec_pF



    integer :: intprec=1, n2, l, n1, k
    real, dimension(:), allocatable :: yy_cDp, yy_Np, yy_cNp, yy_CT, yy_pp, yy_pF , xx,   &
         &  y_cDp, y_Np, y_cNp, y_CT, y_pp, y_pF , x
    real :: phiPi, ctPi, ctPi_min, ctPi_max


    allocate (x(20*intprec))
    allocate (y_cDp(20*intprec))
    allocate (y_Np(20*intprec))
    allocate (y_cNp(20*intprec))
    allocate (y_CT(20*intprec))
    allocate (y_pp(20*intprec))
    allocate (y_pF(20*intprec))

    allocate (xx(20*intprec))
    allocate (yy_cDp(20*intprec))
    allocate (yy_Np(20*intprec))
    allocate (yy_cNp(20*intprec))
    allocate (yy_CT(20*intprec))
    allocate (yy_pp(20*intprec))
    allocate (yy_pF(20*intprec))



    ctPi_min=-1.
    ctPi_max=1.
    call sg20r(ctPi_min,ctPi_max,intprec,x,n1)
    do k=1,n1
       ctPi=x(k)

       call sg20r(0.,2.*pi,intprec,xx,n2) ! n2 - OUT
       do l=1,n2
          phiPi=xx(l)
          call  HNV_free_elepton_ct_ctPi_phiPi(process_ID,Enu,E1,ct,ctPi,phiPi,mN,ml,mN1,mpi, &
               &  yy_cDp(l),yy_Np(l),yy_cNp(l),yy_CT(l),yy_pp(l),yy_pF(l))

          ! write(*,'3(A,g14.5)') 'yy_CT(l)=', yy_CT(l), 'yy_pp(l)=', yy_pp(l),   'yy_pF(l)=', yy_pF(l),
       end do

       call rg20r(0.,2.*pi,intprec,yy_cDp,y_cDp(k))
       call rg20r(0.,2.*pi,intprec,yy_Np,y_Np(k))
       call rg20r(0.,2.*pi,intprec,yy_cNp,y_cNp(k))
       call rg20r(0.,2.*pi,intprec,yy_CT,y_CT(k))
       call rg20r(0.,2.*pi,intprec,yy_pp,y_pp(k))
       call rg20r(0.,2.*pi,intprec,yy_pF,y_pF(k))
    end do
    call rg20r(ctPi_min,ctPi_max,intprec,y_cDP,xsec_cDp)
    call rg20r(ctPi_min,ctPi_max,intprec,y_NP,xsec_Np)
    call rg20r(ctPi_min,ctPi_max,intprec,y_cNP,xsec_cNp)
    call rg20r(ctPi_min,ctPi_max,intprec,y_CT,xsec_CT)
    call rg20r(ctPi_min,ctPi_max,intprec,y_pp,xsec_pp)
    call rg20r(ctPi_min,ctPi_max,intprec,y_pF,xsec_pF)

    deallocate(x,xx,y_cDp,y_Np,y_cNp,y_CT,y_pp,y_pF,  yy_cDp,yy_Np,yy_cNp,yy_CT,yy_pp,yy_pF)

  end subroutine HNV_free_elepton_ct




  ! 1-fold differential xsec dsi / dE_pi
  subroutine HNV_free_ctPi(process_ID,Enu,ctPi,mN,ml,mN1,mpi, xsec_cDp,xsec_Np,xsec_cNp,xsec_CT,xsec_pp,xsec_pF)

    use gauss_integration
    use constants, only : pi
    use lepton_kinematics_free_FULL, only : minmaxE1_Enu_ct

    implicit none

    integer, intent(in) :: process_ID
    real, intent(in) :: Enu, ctPi, mN, ml, mN1, mpi
    real, intent(out) :: xsec_cDp,xsec_Np, xsec_cNp,xsec_CT,xsec_pp,xsec_pF



    integer :: intprec=1, n3, j, n2, l, n1, k
    real, dimension(:), allocatable :: yyy_cDp, yyy_Np,yyy_cNp, yyy_CT, yyy_pp, yyy_pF , xxx, &
         & yy_cDp, yy_Np, yy_cNp, yy_CT, yy_pp, yy_pF ,  xx, &
         &  y_cDp, y_Np, y_cNp, y_CT, y_pp, y_pF , x
    real :: phiPi, ct_min, ct_max, ct, E1_min, E1_max, E1


    allocate (x(20*intprec))
    allocate (y_cDp(20*intprec))
    allocate (y_Np(20*intprec))
    allocate (y_cNp(20*intprec))
    allocate (y_CT(20*intprec))
    allocate (y_pp(20*intprec))
    allocate (y_pF(20*intprec))

    allocate (xx(20*intprec))
    allocate (yy_cDp(20*intprec))
    allocate (yy_Np(20*intprec))
    allocate (yy_cNp(20*intprec))
    allocate (yy_CT(20*intprec))
    allocate (yy_pp(20*intprec))
    allocate (yy_pF(20*intprec))

    allocate (xxx(20*intprec))
    allocate (yyy_cDp(20*intprec))
    allocate (yyy_Np(20*intprec))
    allocate (yyy_cNp(20*intprec))
    allocate (yyy_CT(20*intprec))
    allocate (yyy_pp(20*intprec))
    allocate (yyy_pF(20*intprec))



    ct_min=-1.
    ct_max=1.
    call sg20r(ct_min,ct_max,intprec,x,n1)
    do k=1,n1
       ct=x(k)

       call minmaxE1_Enu_ct(Enu,ct,1.1,mN,ml,E1_min,E1_max)
       call sg20r(E1_min,E1_max,intprec,xx,n2) ! n2 - OUT
       do l=1,n2
          E1=xx(l)


          call sg20r(0.,2.*pi,intprec,xxx,n3) ! n3 - OUT
          do j=1,n3
             phiPi=xxx(j)


             call  HNV_free_elepton_ct_ctPi_phiPi(process_ID,Enu,E1,ct,ctPi,phiPi,mN,ml,mN1,mpi,&
                  & yyy_cDp(j),yyy_Np(j),yyy_cNp(j),yyy_CT(j),yyy_pp(j),yyy_pF(j))

          end do
          call rg20r(0.,2.*pi,intprec,yyy_cDp,yy_cDp(l))
          call rg20r(0.,2.*pi,intprec,yyy_Np,yy_Np(l))
          call rg20r(0.,2.*pi,intprec,yyy_cNp,yy_cNp(l))
          call rg20r(0.,2.*pi,intprec,yyy_CT,yy_CT(l))
          call rg20r(0.,2.*pi,intprec,yyy_pp,yy_pp(l))
          call rg20r(0.,2.*pi,intprec,yyy_pF,yy_pF(l))


       end do
       call rg20r(E1_min,E1_max,intprec,yy_cDp,y_cDp(k))
       call rg20r(E1_min,E1_max,intprec,yy_Np,y_Np(k))
       call rg20r(E1_min,E1_max,intprec,yy_cNp,y_cNp(k))
       call rg20r(E1_min,E1_max,intprec,yy_CT,y_CT(k))
       call rg20r(E1_min,E1_max,intprec,yy_pp,y_pp(k))
       call rg20r(E1_min,E1_max,intprec,yy_pF,y_pF(k))

    end do
    call rg20r(ct_min,ct_max,intprec,y_cDp,xsec_cDp)
    call rg20r(ct_min,ct_max,intprec,y_Np,xsec_Np)
    call rg20r(ct_min,ct_max,intprec,y_cNp,xsec_cNp)
    call rg20r(ct_min,ct_max,intprec,y_CT,xsec_CT)
    call rg20r(ct_min,ct_max,intprec,y_pp,xsec_pp)
    call rg20r(ct_min,ct_max,intprec,y_pF,xsec_pF)

    deallocate(x,xx,xxx, y_cDp, y_Np,y_cNp, y_CT,y_pp,y_pF, yy_cDp, yy_Np,yy_cNp, yy_CT,yy_pp,yy_pF)
    deallocate(yyy_cDp, yyy_Np,yyy_cNp, yyy_CT,yyy_pp,yyy_pF)

  end subroutine HNV_free_ctPi


end module lepton_xsec_free
