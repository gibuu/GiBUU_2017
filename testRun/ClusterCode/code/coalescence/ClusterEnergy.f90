!***************************************************************************
!****m* /ClusterEnergy
! PURPOSE
! routines related to calculation of cluster energy
! NOTES
! * The cluster energies are calculated in an event-by-event basis, which
!   makes an accurate energy determination more difficult...The results
!   also depend on the width of the gaussian's.
! * For this purpose the SMM-part of the cluster code evaluates the
!   excitation energies better (Full Ensemple is used in the SMM part).
!
!***************************************************************************
module ClusterEnergy

  PRIVATE

  logical, save :: flag_init = .true.
  real,    save :: vkf=3.6,skf=8.5,bs=50.,cs=-6., &
       &           ms=550./197.33,ek=197.33,gr=1.4, &
       &           pi=3.14159,wr=0.,varx=0.,ecoul, &
       &           Sk1,Sk2,Sk3,Sk4,Sk5 ! Parameters for Non-Relativistic Skyrme

  PUBLIC :: energy,boost


  contains


  !*************************************************************************
  !****s* ClusterEnergy/energy
  ! NAME
  ! subroutine energy
  !
  ! PURPOSE
  ! Calculation of the total abd binding energy of the produced clusters.
  !
  ! NOTES
  ! * The calculation of the total energy is done within mean-field
  !   prescriptions. In this respect nuclear structure is not properly
  !   taken into account!
  ! * Two options for total energy are used:
  ! * Relativistic     : RMF-Walecka
  ! * Non-Relativistic : Skyrme-soft, MDI-part is not yet implemented!
  !*************************************************************************
    subroutine energy(FMass,PosFra,ParticleVector,MomFra,Etot, &
         &            bxout,byout,bzout)
      use typeDefinitions,  only : particle,Cluster
      use InputCoalescence, only : Get_model
      implicit none
  !*************************************************************************
      type(particle), dimension(:), intent(in) :: ParticleVector
      integer, dimension(:), intent(in)        :: PosFra
      integer, intent(in)                      :: FMass
      real, dimension(0:3), intent(in)         :: MomFra
      real, intent(out)                        :: Etot
      real, intent(out),optional               :: bxout,byout,bzout
      integer :: i,ii,j,jj,jf
      real    :: efmf,ustr,rsq,rsq0,dvs,betax,betay,betaz
      real    :: sd,s1,s2,s3,s,sigma
      real, dimension(1:size(ParticleVector,dim=1)) :: rhosf, &
           & curf0e,curf1e,curf2e,curf3e,coul0,coul1,coul2,coul3
      real, dimension(0:3) :: true4mom
      real                 :: true4momRF,FieldEnergy,tau,jcoul
      real, dimension(1:3) :: xi,xj
      real, dimension(0:3) :: uj,stromi,pin,pout,acoul
      real :: ScalarField,rhoScalar,MomPart,rhobar
      !---------------------------------------------------------------------
      if ( Get_Model.lt.1 .or. Get_Model.gt.2 ) then
         write(*,*) 'from Module ClusterEnergy, subroutine Energy : '
         write(*,*) &
      & 'This Option for Potential is not defined!! Get_Model = ',Get_Model
         write(*,*) '!!! TERMINATION OF PROGRAM NOW !!!'
         STOP
      endif
      !---------------------------------------------------------------------
      if (flag_init) then
         call constants
         flag_init=.false.
      endif
      call GaussianWidth(FMass)
      !---------------------------------------------------------------------
      !calculate the fields at position of the considered cluster with label IN
      !by summing up only the particles which belong to the cluster IN
      !---------------------------------------------------------------------
      betax  = MomFra(1)/MomFra(0)
      betay  = MomFra(2)/MomFra(0)
      betaz  = MomFra(3)/MomFra(0)
      if (present(bxout)) bxout = betax
      if (present(byout)) byout = betay
      if (present(bzout)) bzout = betaz
      Fields : do i=1,FMass
         ii=PosFra(i)
         sd=0.
         curf0e(i) = 0.
         curf1e(i) = 0.
         curf2e(i) = 0.
         curf3e(i) = 0.
         coul0(i)  = 0.
         coul1(i)  = 0.
         coul2(i)  = 0.
         coul3(i)  = 0.
         xi(:) = ParticleVector(ii)%position(:)
         LoopPart : do j=1,FMass
            jf=PosFra(j)
            if (jf.eq.ii) cycle
            xj(:) = ParticleVector(jf)%position(:)
            pin(:)= ParticleVector(jf)%momentum(:)
            call boost(-betax,-betay,-betaz,pout,pin)
            uj(:) = pout(:)/(ParticleVector(jf)%mass*0.19733)
            !uj(:) = ParticleVector(jf)%momentum(:)/(ParticleVector(jf)%mass*0.19733)
            s1    = ( xi(1) - xj(1)  )*uj(1)
            s2    = ( xi(2) - xj(2)  )*uj(2)
            s3    = ( xi(3) - xj(3)  )*uj(3)
            s     = -s1-s2-s3
            rsq   = -( xi(1)-xj(1) )**2 - ( xi(2)-xj(2) )**2 - ( xi(3)-xj(3) )**2 - s**2

            rsq0 = sqrt(-rsq)
            if (rsq0.lt.0.9) rsq0 = 0.9        !  cut-off, protonradius

            rsq = rsq/varx
            dvs = exp(rsq)*wr

            sd = sd + dvs
            curf0e(i) = curf0e(i) + dvs*uj(0)  !currents [fm**-3]
            curf1e(i) = curf1e(i) + dvs*uj(1)
            curf2e(i) = curf2e(i) + dvs*uj(2)
            curf3e(i) = curf3e(i) + dvs*uj(3)
            if (ParticleVector(jf)%Charge.eq.0) cycle
            coul0(i)  = coul0(i)  + uj(0)/rsq0
            coul1(i)  = coul1(i)  + uj(1)/rsq0
            coul2(i)  = coul2(i)  + uj(2)/rsq0
            coul3(i)  = coul3(i)  + uj(3)/rsq0
         end do LoopPart
         rhosf(i) = sd  !fm**-3
      end do Fields

      !---------------------------------------------------------------------
      true4mom(:) = 0.0
      !---------------------------------------------------------------------
      if (Get_Model==1) MomPart = 0.
      LoopPart2 : do jj=1,FMass
         i      = PosFra(jj)
!         efmf   = ( ParticleVector(i)%mass - skf*sigma )
!         uj(:)  = ParticleVector(i)%momentum(:)/(ParticleVector(i)%mass*0.19733)
         pin(:) = ParticleVector(i)%momentum(:)
         call boost(-betax,-betay,-betaz,pout,pin)
         uj(:)  = pout(:)/(ParticleVector(i)%mass*0.19733)
         stromi(0) = curf0e(jj)
         stromi(1) = curf1e(jj)
         stromi(2) = curf2e(jj)
         stromi(3) = curf3e(jj)
         acoul(0)  = coul0(jj)
         acoul(1)  = coul1(jj)
         acoul(2)  = coul2(jj)
         acoul(3)  = coul3(jj)
         ustr = uj(0)*stromi(0)
         do ii=1,3
            ustr = ustr - uj(ii)*stromi(ii)
         end do
         MF_Model : if (Get_Model==2) then
            ScalarField = fshift(stromi(0))
            efmf        = ParticleVector(i)%mass - ScalarField/0.19733
            sigma       = ScalarField/0.19733/skf
            rhoScalar   = (ms**2*sigma + bs*sigma**2 + cs*sigma**3)/skf
            rhosf(jj) = rhoScalar
            if (stromi(0).gt.1.e-4) then
               FieldEnergy = ( 0.5*vkf*(stromi(0)**2+stromi(1)**2+stromi(2)**2+stromi(3)**2) + &
                    &   0.5*ms**2*sigma**2 + bs*sigma**3/3. + cs*sigma**4/4.) / stromi(0)
               if (ParticleVector(i)%Charge==0) tau = -1.
               if (ParticleVector(i)%Charge==1) tau = 1.
               jcoul = stromi(0)*acoul(0) + &
                    & stromi(1)*acoul(1) +stromi(2)*acoul(2) + stromi(3)*acoul(3)
               FieldEnergy = FieldEnergy + 0.5*ecoul*(1.+tau)/2.*jcoul / stromi(0)
            else
               FieldEnergy = 0.0
            endif
         else if (Get_Model==1) then
            if (stromi(0).gt.1.e-4) then
               rhobar = stromi(0)
               FieldEnergy = ( Sk1*rhobar**2 + Sk2*rhobar**(Sk3+1.) ) / rhobar
!!$               call betaToLRF(stromi,bx_LRF,by_LRF,bz_LRF)
!!$               pin(:) = pout(:)/0.19733
!!$               call boost(-bx_LRF,-by_LRF,-bz_LRF,pi1,pin)
!!$               MomPart = 0.
!!$               do ii=1,FMass
!!$                  stromi(0) = curf0e(ii)
!!$                  stromi(1) = curf1e(ii)
!!$                  stromi(2) = curf2e(ii)
!!$                  stromi(3) = curf3e(ii)
!!$                  pin(:) = ParticleVector(PosFra(ii))%momentum(:)/0.19733
!!$                  call boost(-betax,-betay,-betaz,pout,pin)
!!$                  call betaToLRF(stromi,bx_LRF,by_LRF,bz_LRF)
!!$                  pin(:) = pout(:)
!!$                  call boost(-bx_LRF,-by_LRF,-bz_LRF,pi2,pin)
!!$                  pij     = (pi1(1)-pi2(1))**2+(pi1(2)-pi2(2))**2+(pi1(3)-pi2(3))**2
!!$                  MomPart = MomPart + Sk4/(1.+ pij/Sk5**2)
!!$               end do
!!$!               FieldEnergy = FieldEnergy + MomPart/rhobar
            else
               FieldEnergy = 0.0
               MomPart     = 0.0
            endif
            vkf = 0.
         else
            write(*,*) 'this option for Mean-Field Model is not defined!'
            write(*,*) 'Get_Model = ',Get_Model
            write(*,*) '!!! TERMINATION OF PROGRAMM NOW!!!'
            STOP
         endif MF_Model

         true4mom(0) = true4mom(0) + pout(0)/0.19733 + FieldEnergy
         true4mom(1:3) = true4mom(1:3) + pout(1:3)/0.19733 + vkf*stromi(1:3) + &
              & ecoul*acoul(1:3)

      end do LoopPart2

      !--------------------------------------------------------------------------------------
      if (Get_Model==2) then
         true4momRF = sqrt(true4mom(0)**2-true4mom(1)**2-true4mom(2)**2-true4mom(3)**2)
      else
         true4momRF = true4mom(0) + ecoul*acoul(0)
      endif
      !--------------------------------------------------------------------------------------
      Etot = (true4momRF-0.939*float(FMass)/0.19733)*ek/float(FMass) !Etot per nucleon!!!!!!!
      !--------------------------------------------------------------------------------------

    !*************************************************************************
    end subroutine energy ! of energie  **************************************
    !*************************************************************************

    !*************************************************************************
    real function fshift(rho)
    !*************************************************************************
    !****f* ClusterEnergy/fshift
    ! NAME
    ! real function fshift(rho)
    ! PURPOSE
    ! Fit of the nucleon mass shift m - m^* for the various RMF parameter sets
    ! INPUT:
    ! * real, intent(in) :: rho ! -- baryon density (fm**-3)
    ! OUTPUT:
    ! * real :: fshift ! -- m - m^* (GeV)
    !***

      implicit none

      real, intent(in) :: rho ! -- baryon density (fm**-3)
      real, parameter  :: rhoNull = 0.145
      real, parameter  :: m_nucleon = 0.938
      real :: rhorat

      rhorat = rho/rhoNull

      !**** NL2 from A. Lang et al:
      fshift = m_nucleon*(1.-1./(1.+0.201221*rhorat**0.888537))
    !*************************************************************************
    end function fshift !*****************************************************
    !*************************************************************************

  !****************************************************************************************
  !****s* ClusterEnergy/constants
  ! NAME
  ! subroutine constants
  !
  ! PURPOSE
  ! Defines some global constants for mean-field potential
  !
  ! NOTES
  ! needed for the determination of the total energy
  !****************************************************************************************
    subroutine constants
    !*************************************************************************
      real :: rhoSat

      vkf = 3.60713 !fm**2 ! parameter fuer nl2 von promotion V. Koch
      skf = 8.5

      bs = 50.57
      cs = -6.26
      ms = 2.79
      !--------------------------------------------------------------------------------------
      pi = 4.*atan(1.0)
      ek = 197.33
      !--------------------------------------------------------------------------------------
      ! Parameters for SM (Effenberger Thesi's, page 16, SM-Parametrization)
      !--------------------------------------------------------------------------------------
      rhoSat = 0.168
      Sk1    = -108.6/197.33 !MeV-->fm**-1
      Sk1    = Sk1/2./rhoSat
      Sk2    = 136.8/197.33  !MeV-->fm**-1
      Sk3    = 1.26
      Sk2    = Sk2/(Sk3+1.)/(rhoSat**Sk3)
      Sk4    = -63.6/197.33
      Sk4    = Sk4/rhoSat
      Sk5    = 2.13 !fm**-1

      ecoul = 1./137. !elementarladung
    !*************************************************************************
    end subroutine constants
    !*************************************************************************

  !*************************************************************************
  !****s* ClusterEnergy/GaussianWidth
  ! NAME
  ! subroutine GaussianWidth
  !
  ! PURPOSE
  ! Fixes the width of the Gaussians acc. the cluster radius
  !
  ! NOTES
  ! * In an event-by-event basis one increases the gaussian's width,
  !   otherwise the potential contribution of the total energy tends to zero!!!
  ! * needed for the determination of the total energy
  !*************************************************************************
    subroutine GaussianWidth(FMass)
  !*************************************************************************
      implicit none
      integer, intent(in) :: FMass
      !---------------------------------------------------------------------
      gr   = 1.3*float(FMass)**0.33333333+1.0
      wr   = 1./(pi*gr*gr)**1.5
      varx = gr**2
  !*************************************************************************
    end subroutine GaussianWidth !******************************************
  !*************************************************************************


  !*************************************************************************
  !****s* ClusterEnergy/boost
  ! NAME
  ! subroutine boost
  !
  ! PURPOSE
  ! Lorentz-Trafo: frames L and L*, L* moves with velocity vx,vy,vz in L
  !*************************************************************************
    subroutine boost(vx,vy,vz,pout,pin)
      implicit none

      real, intent(in) :: vx,vy,vz
      real, dimension(0:3), intent(in) :: pin
      real, dimension(0:3), intent(out) :: pout

      real :: e,px,py,pz,es,pxs,pys,pzs
      real :: v,gamma,vps,effm

      es  = pin(0)
      pxs = pin(1)
      pys = pin(2)
      pzs = pin(3)

      v = sqrt(vx**2 + vy**2 + vz**2)
      gamma = 1./sqrt(1. - v**2)
      vps = vx*pxs + vy*pys + vz*pzs

      if(1..le.v) then
         write(*,*) 'module ClusterEnergy, routine boost:'
         write(*,*) 'velocity > = 1'
         write(*,*) ' termination'
         STOP
      endif
      if (v.lt.0.000001) then
         effm = 0.
      else
         effm = (gamma - 1.)/v**2*vps + gamma*es
      endif

      e  = gamma*(es  + vps)
      px = pxs + vx*effm
      py = pys + vy*effm
      pz = pzs + vz*effm

      pout(0) = e
      pout(1) = px
      pout(2) = py
      pout(3) = pz

    end subroutine boost

  end module ClusterEnergy


!!$    !*************************************************************************
!!$    subroutine BetaToLRF(strom,b1,b2,b3)
!!$      implicit none
!!$      !***********************************************************************
!!$      real, dimension(0:3),intent(in) :: strom
!!$      real, intent(out) :: b1,b2,b3
!!$      !-----------------------------------------------------------------------
!!$      if (strom(0).gt.1.e-4) then
!!$         b1 = strom(1)/strom(0)
!!$         b2 = strom(2)/strom(0)
!!$         b3 = strom(3)/strom(0)
!!$      else
!!$         b1 = 0.
!!$         b2 = 0.
!!$         b3 = 0.
!!$      endif
!!$    !*************************************************************************
!!$    end subroutine BetaToLRF
!!$    !*************************************************************************
