module StructureFunctions

  implicit none
  PRIVATE

  public :: R_world_average
  public :: F2_from_xsec
  public :: Nachtmann_Q2_W
  public :: intF2_to_duality
  public :: MAIDlike_dE1dcostheta_res


contains


  !*************************************************************************
  !****f* StructureFunctions/R_world_average
  ! NAME
  ! function R_world_average(x,Q2)
  !
  ! PURPOSE
  ! This function return the world-average value of R=FL/2xF1=F2/2xF1(1+4mN^2x^2/Q2)-1
  ! as determined in Whitlow et al Phys. Lett B 250(1990) 193 for 
  ! electron--proton scattering
  !
  ! INPUTS
  ! * real   :: x    ! Bjorken variable
  ! * real   :: Q2   ! squared momentum transfered
  ! 
  ! OUTPUT
  ! * R_world_average
  !
  ! NOTES
  ! nearly the same function (BUT different at Q2<0.35) is in 
  ! /code/init/lowElectron/ParamEP.f90 :
  ! subroutine  CalcParamEP_R1990(W,Q2, R)
  !
  !*************************************************************************
  real function R_world_average(x,Q2)

    real, intent(in) :: x, Q2
    real :: R
    R=0.
!    write(*,'(5(A,g12.5))') 'Enter R_world_average:  x=',x,&
!         & '    Q2=',Q2, '   R is set to zero  R=',R,'  log(Q2/0.04)=', log(Q2/0.04)
    if (Q2<0.35) then
       R=Q2*0.666
       !write(*,*) 'R_world_average is not fitted for Q2<0.3 GeV^2, set to 0.666*Q2'
    else
       R =  0.0635/log(Q2/0.04) * (1.+12.*Q2/(1.+Q2)*0.125**2/(0.125**2+x**2)) +0.5747/Q2 -0.3534/(Q2**2+0.09)
    end if
    !write(*,'(3(A,g12.5))') 'In R_world_average:  x=',x, '    Q2=',Q2, '    R=',R
    R_world_average=R
  end function R_world_average



  real function Nachtmann_Q2_x(Q2,x,mN)
    real, intent(in) :: Q2,x,mN
    Nachtmann_Q2_x = 2.*x/(1.+sqrt(1.+(2.*mN*x)**2/Q2));
  end function Nachtmann_Q2_x



  real function Nachtmann_Q2_W(Q2,W,mN)
    real, intent(in)  :: Q2,W,mN
    real              :: x,nu,xi
    nu=(W**2+Q2-mN**2)/2./mN
    x=Q2/2./mN/nu
    xi=Nachtmann_Q2_x(Q2,x,mN)
!    write(*,'(6(A,g12.5))') '# In Nachtmann_Q2_W(Q2,W,mN):  Q2=',Q2,&
!         & '   W=',W, '   mN=',mN, '   nu=',nu, '   x=',x, '   xi=',xi
    Nachtmann_Q2_W=xi
  end function Nachtmann_Q2_W



  real function Bjorken_Q2_Na (Q2, Na, mN)
    real, intent(in) :: Q2,Na,mN
    Bjorken_Q2_Na = Na/(1.-(mN*Na)**2/Q2);
  end function Bjorken_Q2_Na



  !*************************************************************************
  !****s* StructureFunctions/F2_from_xsec
  ! NAME
  ! subroutine F2_from_xsec(process_ID,ml,Enu,E1,costheta,xsec, mN, F2, x, Q2, xi)
  !
  ! PURPOSE
  ! This subroutine calculates structure function F2 as well as (Bjorken) x, Q2, xi  
  ! from the input differential cross section dsi/dE1 dcostheta
  ! and using the world average of " R=F2/(2xF_1)*(1+4mN^2x^2/Q2)-1 "
  !
  ! INPUTS
  ! * integer :: process_ID   ! ID of the process  (EM,CC)
  ! * real    :: ml           ! mass of the outgoing lepton, 
  ! * real    :: Enu          ! energy of the incoming neutrino or electron
  ! * real    :: E1           ! energy of the outgoing lepton
  ! * real    :: costheta     ! cos(angle between incoming and outgoing leptons)
  ! * real    :: mN           ! mass of the target nucleon or nucleus
  ! * real    :: xsec         ! double-differential dsi/dE1/dcostheta cross section  in 1/GeV^3
  !
  ! 
  ! OUTPUT
  ! * real    :: F2           ! structure function F2
  ! * real    :: x            ! Bjorken variable x
  ! * real    :: Q2           ! momentuim transferred squared
  ! * real    :: xi           ! Nachtman variable x
  !*************************************************************************
  subroutine F2_from_xsec(process_ID,ml,Enu,E1,costheta,xsec, mN, F2, x, Q2, xi)

    use leptonicID
    use constants, only : alphaQED, GF, pi

    integer, intent(in) :: process_ID  ! ID of the process 
    real, intent(in)    :: ml, Enu, E1, costheta 
    real, intent(in)    :: mN          ! mass if the target nucleon or nucleus
    real, intent(in)    :: xsec        ! xsec is double-differential dsi/dE1/dcostheta cross section
    real, intent(out)   :: F2, x, Q2, xi


    real :: coeff, cos2t12, sin2t12, R

    Q2=-ml**2 +2.*Enu*E1 - 2.*Enu*sqrt(E1**2-ml**2)*costheta
    x=Q2/2./mN/(Enu-E1)
    xi=Nachtmann_Q2_x(Q2,x,mN)

    sin2t12=(1.-costheta)/2.
    cos2t12=(1.+costheta)/2.

    if( process_ID.eq.EM .or. process_ID.eq.antiEM) then
       coeff=(2.*alphaQED*E1/Q2)**2    !electromagnetic process
    else
       if ( process_ID.eq.CC .or. process_ID.eq.antiCC .or. process_ID.eq.NC .or. process_ID.eq.antiNC ) then
          coeff=2.*(GF*E1/2./pi)**2 ! weak process
       else
          write(*,'(A,6A6,A)') 'Process with ID=',process_ID, ' is not defined; return'
          stop
       end if
    endif

    !write(*,'(7(A,g12.5))') 'F2_from_xsec:   xsec=',xsec, '   coeff=',coeff, '     cos2t12=',cos2t12, '   sin2t12=',sin2t12, &
    !& '    Q2=',Q2, '   x=',x, '    xi=',xi, '    mN=',mN

    R=R_world_average(x,Q2)

    !write(*,'(2(A,g12.5))') '     R=',R

    F2= (Enu-E1) * xsec/2./pi/coeff/(  cos2t12 + 2.*sin2t12*(Q2+(2.*mN*x)**2)/(2.*mN*x)**2/(1.+ R ) )

  end subroutine F2_from_xsec



  !*************************************************************************
  !****f* StructureFunctions/intF2_to_duality
  ! NAME
  ! function intF2_to_duality(process_ID,Q2,mN,ximin,ximax)
  !
  ! PURPOSE
  ! This function integrates the structure function F2 over the range of Nachtman variable 
  ! from ximin to ximax
  ! 
  !
  ! INPUTS
  ! * integer    :: process_ID ! EM,CC
  ! * real       :: ml, Q2, mN, ximin, ximax
  ! * integer    :: F2from     ! 1=from background xsec, 2=from full xsec
  ! 
  ! OUTPUT
  ! intF2_to_duality           ! integral of F2 over xi
  !*************************************************************************
  real function intF2_to_duality(process_ID,ml,F2from,charge_in,Q2,mN,ximin,ximax,debug_F2)

    use constants, only : mmuon
    use singlePionProductionMAIDlike
    !use neutrinoXsection
    use gauss_integration
    use minkowski, only : abs4
    use leptonicID

    integer, intent(in) :: process_ID
    integer, intent(in) :: charge_in
    integer, intent(in) :: F2from                ! 1=from background xsec, 2=from full xsec
    real, intent(in)    :: Q2, mN, ximin, ximax
    real, intent(in)    :: ml                    ! mass of the outgoing lepton
    logical, intent(in), optional :: debug_F2

    ! for MAIDlike_singlePi xsec
    integer                     :: charge_out,pion_charge_out
    real                        :: costheta, sintheta, E1, k1, Enu, knu
    real, dimension(0:3)        :: k_in, k_out, p_in
    real, dimension(1:3)        :: position
    logical                     :: nuclear_phasespace

    ! for integration
    real, allocatable, dimension(:) :: xx,yy
    integer :: nn,k, nn1

    !
    real :: nu, xsec, F2, F2int
    real :: x_out, x, xi, xi_out, Q2_out


    ! construct position for MAIDlike_singlePi function, which calculates 1-pion background
    position=(/0.,0.,0./)

    ! organize a loop either in costheta or in E1 for a given Q2 and a range of xi
    Enu=3.
    if (Q2>3.) Enu=5.
    if (Q2>6.) Enu=9.
    if(process_ID.eq.CC) knu=enu 
    if(process_ID.eq.EM) knu=sqrt((Enu-ml)*(Enu+ml))

    k_in=(/Enu,0.,0.,knu/)

    ! begin to integrate F2 over xi
    nn=1
    allocate (xx(64*nn))
    allocate (yy(64*nn))

    nuclear_phasespace=.false.
    p_in = (/mN, 0., 0., 0./)

    write(*,*) ''
    !write(*,'(2(A,g12.5))') '#   Q2=',Q2


    call sg64r(ximin,ximax,nn,xx,nn1)
    integration_loop: do k=1,nn1 
       xi=xx(k)
       x=Bjorken_Q2_Na(Q2, xi, mN)


       ! construct incoming vectors for the x-sec
       nu=Q2/2./mN/x
       E1=Enu-nu

       !write(*,'(1(A,g12.5))') '#   ml=',ml
       !write(*,'(3(A,g12.5))') '#   x=',x,'   xi=',xi,  '    nu=',nu, '    E1=',E1


       k1=sqrt((E1+ml)*(E1-ml))
       costheta = (  -( Q2 + ml**2 ) + 2.*Enu*E1 )/2./Enu/k1
       sintheta=sqrt((1.-costheta)*(1.+costheta))
       k_out = (/E1, k1*sintheta, 0., k1*costheta/)

       !write(*,'((A,g12.5))') '#       costheta=',costheta


       !write(*,'(4(A,I5))') '# brgF2_to_duality:  process_ID=', process_ID, '   charge_in=',charge_in
       !write(*,'(A,4g12.5)') '#   k_in=',k_in
       !write(*,'(A,4g12.5)') '#   k_out=',k_out
       !write(*,'(A,4g12.5)') '#   p_in=',p_in
       !write(*,'(A,3g12.5)') '#   position=',position

       ! summing over outgoing nucleon charges
       xsec=0

       if (F2from.eq.3) then ! choice if we sum up over outgoing nucleon charges
          ! charge_out == resonance charge
          if (process_ID .eq. CC ) charge_out=charge_in+1
          if (process_ID .eq. EM ) charge_out=charge_in
          xsec=xsec+MAIDlike_dE1dcostheta_res(process_ID,charge_in,charge_out,k_in,k_out,costheta,  & 
               & p_in,position,nuclear_phasespace,ml)
          write(*,'(A,3g12.5)') '#   charge_out=',(charge_in+1), '   xsec=',xsec

       else ! choice if we sum up over outgoing nucleon charges


          do charge_out=0,1
             if (process_ID.eq.EM) then 
                pion_charge_out=charge_in-charge_out 
             else
                pion_charge_out=charge_in+1-charge_out 
                if (pion_charge_out>1 .or. pion_charge_out<-1) cycle
             end if

             ! calculate bgr--double-diff-xsec  from MAIDlike__dE1dcostheta
             if (F2from.eq.1) then
                xsec=xsec+MAIDlike_dE1dcostheta_bgr(process_ID,charge_in,charge_out,pion_charge_out,k_in,k_out,costheta,  &
                     & p_in,position,nuclear_phasespace,ml)
             end if

             if (F2from.eq.2) then
                xsec=xsec+MAIDlike_dE1dcostheta_full(process_ID,charge_in,charge_out,pion_charge_out,k_in,k_out,costheta,  &
                     & p_in,position,nuclear_phasespace,ml)
                !write(*,'(A,3g12.5)') '#   xsec=',xsec
             end if


             ! if (F2from.eq.3 .and. process_ID.eq.EM )then
             ! xsec=xsec+MAIDlike_dE1dcostheta_res(process_ID,charge_in,charge_out,k_in,k_out,costheta,  & ! charge_out == resonance charge for CC, nucleon_charge for EM
             ! & p_in,position,nuclear_phasespace,ml)
             ! write(*,'(A,3g12.5)') '#   charge_out=',charge_out, '   xsec=',xsec
             ! end if
          end do

       end if ! choice if we sum up over outgoing nucleon charges


       !write(*,'(A,4g12.5,A,g12.5)') '#   p_out=',p_out,  ' check p_out^2=', abs4(p_out)
       !write(*,'(A,4g12.5,A,g12.5)') '#   ppi_out=',pion_momentum_out  ,  ' check ppi_out^2=', abs4(pion_momentum_out)
       !write(*,'(A,g12.5)') '# MAID-like: xsec=', xsec  

       ! calculate F2 from xsec
       call F2_from_xsec(process_ID,ml,Enu,E1,costheta,xsec, mN, F2, x_out, Q2_out, xi_out)
       write(*,'(7(A,g12.5))') '# Q2=',Q2, '     F2=',F2,  '     x=',x, '    x_out=',x_out, '    xi_in=',xi,  '    xi_out=',xi_out

       if ( present(debug_F2)) then
          if (debug_F2) then
             write(11,'(7g14.5)')  Q2, costheta, E1, xsec,    x, xi_out, F2
          end if
       end if

       ! minimal  check if calculations are self-consistent 
       if (abs(x-x_out)>1.e-4 ) then
          write(*,'(3(A,g14.6))') 'Insufficient accuracy, x=',x, '   and x_out=', x_out, '    differ too much. Stop'
          stop
       else
          yy(k)=F2
       end if

    end do integration_loop

    ! the integration itself
    call rg64r(ximin,ximax,nn,yy,F2int)

    write(11,*) ''
    write(11,*) ''

    intF2_to_duality=F2int

  end function intF2_to_duality















  !*************************************************************************
  !****f* neutrinoXsection/MAIDlike_dE1dcostheta_res
  ! NAME
  ! function MAIDlike_dE1dcostheta_res(process_ID,charge_in,charge_out,k_in,k_out,costheta,  &
  !     & p_in,p_out,pion_momentum_out,position,nuclear_phasespace,ml_out)
  !
  ! PURPOSE
  ! This function calculates double differential xsec by gauss integrating MAIDlike_singlePi over phiPi and thetaPi
  ! shows the sum of the cross sections with different charges of outgoing pions
  ! USED for calculating structure function F2 by  brgF2_to_duality
  ! actually the resonance contribution is not directly related to MAID, but I didn' found a better name yet
  !
  ! INPUTS
  ! the same as MAIDlike_singlePi  
  ! 
  ! OUTPUT  
  ! the same as MAIDlike_singlePi,  units 1/GeV^3
  !*************************************************************************
  real function MAIDlike_dE1dcostheta_res(process_ID,charge_in,charge_out,k_in,k_out,costheta,  &
       & p_in,position,nuclear_phasespace,ml_out)

    use electronPionProd_medium_EN, only : dSdO_fdE_fdO_k_med_eN
    use eN_eventDefinition, only :  electronNucleon_event
    use eN_event, only           :  init_electronNucleon_event
    use particleDefinition
    use eN_event, only : nuclearFluxFactor_correction
    use random,only : rn, rnCos
    use constants, only: pi, mN
    use degRad_conversion, only : degrees
    use resProd_lepton, only : dSdO_fdE_fdO_k_med_res_EN
    use leptonicID
    use gauss_integration
    use neutrinoXsection, only: Xsec_dSigmadCosThetadElepton,  XsecdCosthetadElepton
    use idTable

    real, dimension(0:3), intent(in) :: k_in, p_in, k_out
    real, dimension(1:3), intent(in) :: position
    real,                 intent(in) :: costheta
    integer,              intent(in) :: charge_in
    integer,              intent(in) :: charge_out
    real,                 intent(in) :: ml_out

    integer,              intent(in) :: process_ID
    logical,              intent(in) :: nuclear_phasespace

    type(electronNucleon_event) :: eNev0, eN
    real ,dimension(0:3) :: kf,pf,ppif 
    real ,dimension(-1:1) :: sigmaRes3
    real :: phiPi, cosThetaPi, thetaPi
    real :: deltaMonteCarlo, sigmaRes, sigmaPion, massf
    type(particle) :: targetNuc
    logical :: success
    real :: sig
    integer :: flavor_ID

    ! for integration
    real, allocatable, dimension(:) :: xx,yy,x,y
    integer                         :: nn, l,nn1,  m,nn2, res, chargef 
    logical                         :: Oliver=.false.! I decides to calculate the xsec for resoanances 
    ! with the help of the Tina's function for both EM and CC processed


    sig=0.

    ! in this code resonace contributions are calculated via different methods for EM and CC processed
    ! now this is effectively erased

    ! Construct particle which is free and has the 3-momentum of the considered nucleon
    call setToDefault(targetNuc)
    targetNuc%momentum=p_in
    targetNuc%charge=charge_in
    targetNuc%id=nucleon
    targetNuc%mass=mN        ! m=m_0
    targetNuc%momentum(0)=freeEnergy(targetNuc)! E=sqrt(p(1:3)^2+m_0^2)
    targetNuc%position(:)=1000.                ! Far away the nucleus
    targetNuc%perturbative=.true.

    ! since we know ml_out, we "recalculate" the flavor:
    if (ml_out.gt.1.0) then
       flavor_ID = 3
    else if (ml_out.gt.0.01) then
       flavor_ID = 2
    else
       flavor_ID = 1
    end if


    call eNev_SetProcess(eNev0, process_ID,flavor_ID)
    call eNev_init_nuStep1(eNev0,targetNuc)

    if (Oliver) then
       select case (process_ID)

          ! xsec of all resonances for el-m electron scattering 
       case(EM) 

          ! Check for Monte-Carlo-Integrations:
          deltaMonteCarlo=2.*pi   !2pi from lepton_phi

          nn=1
          allocate (xx(20*nn))
          allocate (yy(20*nn))

          allocate (x(20*nn))
          allocate (y(20*nn))


          ! loop over costhetaPi    
          call sg20r(-1.,1.,nn,xx,nn1)
          thetaPi_loop: do l=1,nn1 
             cosThetaPi=xx(l)
             thetaPi=degrees(acos(cosThetaPi))

             ! loop over PhiPi
             call sg20r(0.,2.*pi,nn,x,nn2)
             phiPi_loop: do m=1,nn2 
                phiPi=degrees(x(m))

                call init_electronNucleon_event(en, k_in,k_out,targetNuc)
                ! called to determine the kinematics:  ppif and pf
                sigmaPion=dSdO_fdE_fdO_k_med_eN(eN, 0, phiPi, thetaPi, ppif, pf, process_ID,pionNucleonSystem=2)   
                ! Resonance contribution:
                sigmaRes3=dSdO_fdE_fdO_k_med_res_EN(eN,ppif,pf,process_ID,pionNucleonSystem=2) 

                sigmaRes=sigmaRes3(-1)+sigmaRes3(0)+sigmaRes3(1)


                if(sigmaRes.lt.1E-20) then
                   sigmaRes=0
                   y(m)=0 
                   !write(*,'(3(A,g12.5))') '# In MAIDlike_dE1dcostheta  thetaPi=',thetaPi, '   phiPi=',phiPi, '   no solution for kinematics found'
                   cycle ! no solution to kinematical problem was found
                end if

                sigmaRes=sigmaRes/0.38938  ! convertion from mbarn to GeV^{-2}

                y(m)=sigmaRes

                !write(*,'(3(A,g12.5))') '# In MAIDlike_dE1dcostheta  thetaPi=',thetaPi, '   phiPi=',phiPi, '   full xsec=',sigmaPion

             end do phiPi_loop
             ! the integration over phiPi itself
             call rg20r(0.,2.*pi,nn,y,yy(l))
             !write(*,'(A,I5,2(A,g12.5))') '# In MAIDlike_dE1dcostheta  l=',l, '   cosThetaPi=',cosThetaPi, '   yy(l)=',yy(l)

          end do thetaPi_loop
          ! the integration over costhetaPi itself 
          call rg20r(-1.,1.,nn,yy,sig)

          deallocate(x,xx,y,yy)
          !write(*,'(2(A,g12.5))') '# In MAIDlike_dE1dcostheta sig=',sig
          sig=sig*deltaMonteCarlo

          ! xsec of all resonances for CC neutrino scattering
       case(CC) 
          ! summing over all resonances
          particle_ID_loop: do res=delta, F37_1950

             if(res.eq.S11_2090) cycle
             if(res.eq.D13_2080) cycle
             if(res.eq.G17_2190) cycle
             if(res.eq.P11_2100) cycle
             if(res.eq.P13_1900) cycle
             if(res.eq.F15_2000) cycle
             if(res.eq.S31_1900) cycle
             if(res.eq.D33_1940) cycle
             if(res.eq.D35_1930) cycle
             if(res.eq.D35_2350) cycle
             if(res.eq.P31_1750) cycle
             if(res.eq.F35_1750) cycle


!!$             call call eNev_SetProcess(eNev0, process_ID,flavor_ID)
!!$             call eNev_init_nuStep1(eNev,realParticles(i,j))
!!$             call eNev_init_nuStep2(eNev,enu,IP,flagOK)
!!$             if (.not.flagOK) cycle
!!$             call eNev_init_nuStep3a(eNev,eleptonint,costhetaint,flagOK)
!!$             if (.not.flagOK) cycle
!!$             call XsecdCosthetadElepton(eNev,IP,OutPart,yy(l))

!!$             call XsecdCosthetadElepton(process_ID,k_in,p_in,position,  &
!!$                  & charge_in,k_out(0),costheta,res,kf,pf,massf,chargef,ppif,1,sigmaRes) ! 1=pion-charge - is not used for resonances

             write(*,*) '#                  res_ID=', res, '      sigmaRes=',sigmaRes
             sig=sig+sigmaRes
          end do particle_ID_loop

          sig=sig/0.38938*1e-11 ! from 10^{-38} cm^2   to GeV^-2


       case DEFAULT  
          sig=0

       end select

    else ! === if(.not.Oliver)...

       res_loop: do res=delta, F37_1950

          if(res.eq.S11_2090) cycle
          if(res.eq.D13_2080) cycle
          if(res.eq.G17_2190) cycle
          if(res.eq.P11_2100) cycle
          if(res.eq.P13_1900) cycle
          if(res.eq.F15_2000) cycle
          if(res.eq.S31_1900) cycle
          if(res.eq.D33_1940) cycle
          if(res.eq.D35_1930) cycle
          if(res.eq.D35_2350) cycle
          if(res.eq.P31_1750) cycle
          if(res.eq.F35_1750) cycle

!!$          call call eNev_SetProcess(eNev0, process_ID,flavor_ID)
!!$          call eNev_init_nuStep1(eNev,realParticles(i,j))
!!$          call eNev_init_nuStep2(eNev,enu,IP,flagOK)
!!$          if (.not.flagOK) cycle
!!$          call eNev_init_nuStep3a(eNev,eleptonint,costhetaint,flagOK)
!!$          if (.not.flagOK) cycle
!!$          call XsecdCosthetadElepton(eNev,IP,OutPart,yy(l))


!!$          call XsecdCosthetadElepton(process_ID,k_in,p_in,position,  &
!!$               & charge_in,k_out(0),costheta,res,kf,pf,massf,chargef,ppif,1,sigmaRes) ! 1=pion-charge - is not used for resonances

          write(*,*) '#                  res_ID=', res, '      sigmaRes=',sigmaRes
          sig=sig+sigmaRes
       end do res_loop

       if (process_ID.eq.CC) sig=sig/0.38938*1e-11  ! from 10^{-38} cm^2   to GeV^-2
       if (process_ID.eq.EM) sig=sig/0.38938*1e-6   ! from 10^{-33} cm^2 (nanobarn)  to GeV^-2

    end if ! end if not Oliver

    MAIDlike_dE1dcostheta_res=sig

  end function MAIDlike_dE1dcostheta_res


end module StructureFunctions
