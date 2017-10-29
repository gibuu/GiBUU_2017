!***************************************************************************
!****m* /fragmentationHyp
! NAME
! module fragmemtationHyp
!
! FUNCTION
! This module contains the main routine of coalescence between SMM-clusters 
! and hyperons.
! NOTES
! * This routine has a similar structure as the module HypCoala.f90, however, 
!   due to different field allocations between smm and coalescence, we use 
!   two different modules for hypercluster formation within coalescence and 
!   smm models, respectively.
! * Presently only for spectator fragmentation (ALADIN_Flag==true) !!!!!
!   This is because SMM code does not apply well for expanding fireballs.
!***************************************************************************
module fragmentationHyp

  use typeDefinitions, only : cluster,particle,quelle
  implicit none

  PRIVATE

  PUBLIC :: MultiFragmentationHyp

contains 

  !*************************************************************************
  !****s* /MultiFragmentationHyp
  ! NAME
  ! module MultiFragmentationHyp
  !
  ! FUNCTION
  ! The main routine of hyperfragment formation.
  !*************************************************************************
  subroutine MultiFragmentationHyp(NumEns,TheSource,FV,PV)
    !-----------------------------------------------------------------------
    ! Input-Output variables
    !-----------------------------------------------------------------------
    integer,                        intent(in)    :: NumEns
    type(quelle),  dimension(:,:),  intent(in)    :: TheSource
    type(particle), dimension(:,:), intent(inOut) :: PV
    type(cluster),  dimension(:,:), intent(inOut) :: FV
    !-----------------------------------------------------------------------
    ! Local variables
    !-----------------------------------------------------------------------
    type(cluster)                               :: fragment
    type(particle), dimension(1:Size(PV,dim=2)) :: realParticle
    integer                                     :: i,j,NumPart
    type(quelle)                                :: SourceSpect
    integer :: js,jj
    !-----------------------------------------------------------------------
    js=0
    do j=1,Size(TheSource,dim=2)
       if ( .not. TheSource(NumEns,j)%status ) cycle
       !       if (TheSource(NumEns,j)%Type /= 2) CYCLE
       if (TheSource(NumEns,j)%Type == 3) CYCLE
       js = j !--> eigentlich unnoetig...

       !--> eigentlich unnoetig...
       if (js==0) RETURN !no valid source (important if ALADIN_Flag==true)
       !<-- eigentlich unnoetig...


       !Applying HypCoalescence to spectators (projectile and target)
       SourceSpect = TheSource(NumEns,js)
       !check that it is realy projectile:
       !    if (SourceProj%velocity(3)<0.) then
       !       write(*,*) 'Module/routine MultifragmentationHyp:'
       !       write(*,'(A,2i5,2f8.4)') 'Wrong projectile selected...',&
       !            & NumEns,js,SourceProj%position(1),SourceProj%velocity(3)
       !       write(*,*) 'STOP'
       !       STOP
       !    endif

       !coalescence between SMM-fragments & Hyperons
       NumPart     = Size(PV,dim=2)
       SMMevents : do i=1,size(FV,dim=1)
          SMMclusters : do jj=1,size(FV,dim=2)
             fragment    = FV(i,jj)
             realParticle(:) = PV(NumEns,:)
             call FormHypfragments(NumPart,SourceSpect,fragment,realParticle)
             FV(i,jj) = fragment
          end do SMMclusters
       end do SMMevents

       !Separate indirect coalescence between SMM-fragments & pions
       !            pion+cluster->hypernuclei+kaon
       SMMevents2 : do i=1,size(FV,dim=1)
          SMMclusters2 : do jj=1,size(FV,dim=2)
             fragment    = FV(i,jj)
             realParticle(:) = PV(NumEns,:)
             call ToyModel(NumPart,fragment,realParticle)
             FV(i,jj)     = fragment
             PV(NumEns,:) = realParticle(:)
          end do SMMclusters2
       end do SMMevents2
       
    end do


  !*************************************************************************
  end subroutine MultiFragmentationHyp !************************************
  !*************************************************************************

  !*************************************************************************
  !****s* /FormHypFragments
  ! NAME
  ! module FormHypFragments
  
  ! FUNCTION
  ! Coalescence between SMM-clusters & hyperons (Lambdas,Sigmas).
  !*************************************************************************
  subroutine FormHypFragments(NumPart,SourceSpect,fragment,realParticle)

    use inputSMM, only : RadiusHypC, RadiusHypP, HypCoalaMethod

    !-----------------------------------------------------------------------
    ! Input-Output variables
    !-----------------------------------------------------------------------
    integer,                             intent(in)    :: NumPart
    type(quelle),                        intent(in)    :: SourceSpect
    type(cluster),                       intent(inOut) :: fragment
    type(particle), dimension(1:NumPart),intent(inOut) :: realParticle
    !-----------------------------------------------------------------------
    ! Local variables
    !-----------------------------------------------------------------------
    logical,dimension(1:NumPart) :: HypBound
    integer, dimension(1:10) :: hypNumber
    integer :: j,i,HypFound,mhyp,Lambda,Sigma
    real :: rabs,MassNew,pdiff,SourceRadius,FraRadius,hypDist,Ef
    real, dimension(1:3) :: rdist,PosNew,MomPart,MomFra,Pcoala
    real, dimension(0:3) :: MomNew
    logical, save :: initFlag=.true.
    real, save :: RadiusHypV=0.3
    !-----------------------------------------------------------------------
    if (initFlag) then
       Ef = sqrt(RadiusHypP**2 + 0.938**2)
       RadiusHypV = RadiusHypP / Ef
       initFlag = .false.
    endif

    HypBound(:) = .false.
    HypFound    = 0 !number of bound hyperons to the actual cluster

    !Radius of the residual nuclei:
    SourceRadius = 1.2*float(SourceSpect%Size)**0.3333333333+0.5

    !Radius of a fragment + coalescence radius (default coalescence):
    FraRadius    = 1.2*float(fragment%MassNumber)**0.3333333333+RadiusHypC

    !-----------------------------------------------------------------------
    Loop_over_Hyperons : do j=1,NumPart !-----------------------------------
    !-----------------------------------------------------------------------
       if (realParticle(j)%ID < 32 .or. realParticle(j)%ID > 33) cycle
       if (hypBound(j)) cycle !hyperon already bound

       rdist(1:3)   = realParticle(j)%position(1:3) - SourceSpect%position(1:3)
       hypDist      = sqrt( dot_product(rdist,rdist) )
       MomPart(1:3) = realParticle(j)%momentum(1:3) / realParticle(j)%momentum(0)
!       MomFra(1:3)  = fragment%momentum(1:3)/fragment%Mass
       MomFra(1:3)  = fragment%momentum(1:3)/fragment%momentum(0)
       Pcoala(:)    = MomPart(:)-MomFra(:)
       pdiff        = sqrt( dot_product(Pcoala,Pcoala) )
       rdist(:)     = realParticle(j)%position(:)-fragment%position(:)
       rabs         = sqrt( dot_product(rdist,rdist) )

       !--------------------------------------------------------------------
       ! Coalescence between SMM-clusters and strange Baryons (Lambda & Sigma)
       ! Three different coalescence selections: 
       ! (1) in coordinate & momentum space --> HypCoalaMethod = 1 (DEFAULT)
       ! (2) only in coordinate space       --> HypCoalaMethod = 2
       ! (3) only in momentum space (for bound hyperons) --> HypCoalaMethod = 3
       !--------------------------------------------------------------------
       Select Case(HypCoalaMethod)

       Case(1) 
          !R- & P-coalescence useful in central HIC (fireballs):
          if ( ( rabs < abs(FraRadius-RadiusHypC) ) .or. & 
             & ( rabs > abs(FraRadius-RadiusHypC)  .and. & 
             &   rabs < FraRadius .and. pdiff < RadiusHypV) ) then
             mhyp=j
          else
             mhyp = 0
          end if

       Case(2)
          !very rough estimation (only R-coalescence); use only for checks:
          if (rabs > FraRadius) then
             mhyp=0 
          else
             mhyp = j
          endif

       Case(3)
          !for p+X no mean of coalescence in R, since fragments are in any case very slow.
          if (pdiff < RadiusHypV) then
             mhyp = j
          else
             mhyp = 0
          endif

          Case(4)
             ! for X+X (HIC) in spectator fragmentation:
             ! include hyperons only inside spectator source, then p-coalescence:
             if ( hypDist > (SourceRadius+RadiusHypC) ) then
                mhyp = 0
             else
                if (pdiff < RadiusHypV) then
                   mhyp = j
                else
                   mhyp = 0
                endif
             endif

       Case default
          write(*,*) 'Module MultiFragmentationHyp / routine FormHypFragments:'
          write(*,*) & 
               & 'Wrong input for coalescence method between SMM-clusters & Hyperons!'
          write(*,*) 'HypCoalaMethod = ',HypCoalaMethod
          write(*,*) 'STOP'
          STOP
             
       end Select
       !--------------------------------------------------------------------
       if (mhyp .eq. 0) then
          cycle
       else
          HypFound            = HypFound + 1
          if(HypFound > 2) then
             write(*,*) 'module Multifragm.-Hyp:'
             write(*,*) 'not allow for hypers with more than 2 lambdas! HypFound = ',HypFound,mhyp
             HypFound = HypFound - 1
             EXIT Loop_over_Hyperons
          endif
          hypBound(mhyp)      = .true.
          hypNumber(HypFound) = mhyp
       endif

    !-----------------------------------------------------------------------
    end do Loop_over_Hyperons !---------------------------------------------
    !-----------------------------------------------------------------------

    if (HypFound==0) RETURN !no hypfragment found

    !update elements of the type(cluster) "fragment%..."
    !We do not update Mass- and Charge-Numbers
!    fragment%MassNumber  = fragment%MassNumber + HypFound
    PosNew(:) = fragment%position(:)*fragment%mass ![fm*GeV] !!!
    MassNew   = fragment%mass     ![GeV]!!!
    MomNew(:) = fragment%momentum(:) 
    Lambda    = 0
    Sigma     = 0
    do i=1,HypFound
       j=hypNumber(i)
       PosNew(:) = PosNew(:) + realParticle(j)%position(:)*realParticle(j)%mass
       MassNew   = MassNew   + realParticle(j)%mass
       MomNew(:) = MomNew(:) + realParticle(j)%momentum(:)
       if (realParticle(j)%ID==32) Lambda = Lambda + 1
       if (realParticle(j)%ID==33) Sigma  = Sigma  + 1
    end do

    fragment%position(:) = PosNew(:)/MassNew
    fragment%Mass        = MassNew
    fragment%momentum(:) = MomNew(:)
    fragment%HypNumber   = HypFound

    !store index of bound hyperons for detailed analysis:
    fragment%HypContent(1:HypFound)  = hypNumber(1:HypFound)
    
    !Determine hyperon content of SMM-clusters:

    !Lambda-Hypernuclei
    if (Sigma==0) then
       if (Lambda==1) fragment%HypType = 'L' !1-Lambda-Hypernuclei
       if (Lambda==2) fragment%HypType = 'LL' !2-Lambda-Hypernuclei
       if (Lambda > 2) then
          write(*,*) 'Module fragmentationHyp / routine FormHypFragments:'
          write(*,*) 'Warning: cluster with more than 2 Lambdas? '
          fragment%HypType = 'ML' !M-Lambdas-hyperon
       endif
    endif

    !Sigma-Hypernuclei
    if (Lambda==0) then
       if (Sigma==1) fragment%HypType = 'S' !1-Sigma-Hypernuclei
       if (Sigma==2) fragment%HypType = 'SS' !2-Sigma-Hypernuclei
       if (Sigma > 2) then
          write(*,*) 'Module fragmentationHyp / routine FormHypFragments:'
          write(*,*) 'Warning: cluster with more than 2 Sigmas? '
          fragment%HypType = 'MS' !M-Lambdas-hyperon
       endif
    endif

    !Double LambdaSigma-Hypernuclei:
    if (Lambda /= 0 .and. Sigma /= 0) then
       fragment%HypType = 'LS' !Lambda-Hypernuclei
    endif

  !*************************************************************************
  end subroutine FormHypFragments !****************************************
  !*************************************************************************

  !*************************************************************************
  !****s* /ToyModel
  ! NAME
  ! module ToyModel
  
  ! FUNCTION
  ! Very simple toy model for separating hypernuclei oroginating from 
  ! direct pion+nucleon scattering inside spectators.
  ! 
  !*************************************************************************
  subroutine ToyModel(NumPart,fragment,realParticle)

    use inputSMM, only : RadiusHypP

    !-----------------------------------------------------------------------
    ! Input-Output variables
    !-----------------------------------------------------------------------
    integer,                             intent(in)    :: NumPart
    type(particle), dimension(1:NumPart),intent(in)    :: realParticle
    type(cluster),                       intent(inOut) :: fragment
    !-----------------------------------------------------------------------
    ! Local variables
    !-----------------------------------------------------------------------
    logical,dimension(1:NumPart) :: PionBound
    integer, dimension(1:100) :: pionNumber
    integer :: j,PionFound,mpion,pCount
    real :: rabs,pdiff,SourceRadius,SRT,SIG,ProbHyp
    real, dimension(1:3) :: rdist,VV
    real, dimension(0:3) :: MomPion, MomF, Pdummy, pOut
    !-----------------------------------------------------------------------
    PionBound(:) = .false.
    PionFound    = 0 !number of bounded piona to the actual cluster

    SourceRadius = 1.2*float(fragment%MassNumber)**0.3333333333+2.52

    ProbHyp = 0.0
    pCount  = 0
    !-----------------------------------------------------------------------
    Loop_over_Pions : do j=1,NumPart !--------------------------------------
    !-----------------------------------------------------------------------
       if (realParticle(j)%ID /= 101) cycle
       if (pionBound(j)) cycle !pion already bounded

       MomPion(0:3) = realParticle(j)%momentum(0:3)
       MomF(0:3)  = fragment%momentum(0:3)/fragment%Mass
       Pdummy(:)  = MomPion(:)+MomF(:)
       pdiff      = Pdummy(0)**2 - dot_product(Pdummy(1:3),Pdummy(1:3))
       if (pdiff < 0.00001) then
          write(*,*) 'Module MultifragmentationHyp / routine ToyModel:'
          write(*,*) 'Error in calculating total kin. energy!!!'
          write(*,*) 'Pion momentum   : ',MomPion
          write(*,*) 'Fragm. momentum : ',MomF
          write(*,*) 'Pdiff           : ',Pdiff
          STOP
       else
          SRT = sqrt(pdiff)
       endif

       rdist(:) = realParticle(j)%position(:)-fragment%position(:)
       rabs     = sqrt( dot_product(rdist,rdist) )

       !--------------------------------------------------------------------
       ! 2-body phase space
       !--------------------------------------------------------------------
       call phSpace2(0.494, 1.116, srt**2, pOut) !2-body phase space in cms
       VV(1:3) = Pdummy(1:3) / Pdummy(0) !boost parameters
       call boost(VV,pOut) !Lambda momentum in calculational frame

       Pdummy(1:3)  = pOut(1:3)-MomF(1:3)
       pdiff = sqrt( dot_product(Pdummy(1:3),Pdummy(1:3)) )

       !--------------------------------------------------------------------
       if ( ( rabs < abs(SourceRadius-2.52) ) .or. & 
            & ( rabs > abs(SourceRadius-2.52)  .and. & 
            &   rabs < SourceRadius .and. pdiff < RadiusHypP) ) then
          mpion=0 
       else
          mpion = j
       endif
       !--------------------------------------------------------------------
       if (mpion .eq. 0) then
          cycle
       else
          call HuangLam(SRT,SIG)
          ProbHyp = ProbHyp + SIG / 20. !assume sig_tot=20 mb 
                                        !(valid only for high energies!!!)
          pCount  = pCount + 1
          pionBound(mpion)      = .true.
          PionFound            = PionFound + 1
          if(PionFound > 100) then
             write(*,*) 'module Multifragm.-Pion:'
             write(*,*) 'dimension overflow: PionFound = ',PionFound,mpion
             write(*,*) 'Termination of the program'
             STOP
          endif
          pionNumber(PionFound) = mpion
       endif

    !-----------------------------------------------------------------------
    end do Loop_over_Pions !------------------------------------------------
    !-----------------------------------------------------------------------

    if (PionFound==0) RETURN !no pionicfragment found

    ProbHyp = ProbHyp / float(pCount)

    if (ProbHyp > 1.) then
       write(*,*) 'ToyModel: wrong determination of P...',ProbHyp
       STOP 'STOP'
    end if

    if (ProbHyp < 0.0000000001) ProbHyp = 0.0

    !update element of the type(cluster): "fragment%PionNumber"
    fragment%Pionic  = ProbHyp !PionFound

  contains 

    !******************************************************************
    ! real sig = K+ cross section im mb (with Lambda production)
    SUBROUTINE huanglam(SRT,SIG)
      real, intent(in) :: srt
      real,intent(out):: sig
      real s0
      if(srt.lt.1.612) then
         sig=0.
         return
      end if
      S0=1.612
      SIG=0.007665*(SRT-S0)**0.1341/((SRT-1.72)**2+0.007826)
    END SUBROUTINE huanglam

!C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 
!C 
!C      2 particles randomly distributed in phase space           
!C 
    subroutine phspace2(pmas1,pmas2,s,pOut) 
!C 
!C///////////////////////////////////////////////////////////////////// 

      use random, only : rn

      real, intent(in) :: pmas1, pmas2, s
      real, dimension(0:3), intent(out) :: pOut

      REAL,  parameter :: pi=3.1415927 

!      real :: rand

      real :: e,px,py,pz,pr,c1,phis,ww,s1,cp,sp
 
      pr = s - (pmas1 + pmas2)**2 
      pr = pr*(s - (pmas1 - pmas2)**2) 
      pr = sqrt(pr/(4.*s)) 
      C1 = 1. - 2.*RN() 
      PHIS = 2.*pi*RN() 
      WW = 1. - C1*C1 
      IF(WW.LT.0.)WW = 0. 
      S1 = SQRT(WW) 
      CP = COS(PHIS) 
      SP = SIN(PHIS) 
      PZ = PR*C1 
      PX = PR*S1*CP 
      PY = PR*S1*SP 
      E  = SQRT(PMAS2**2 + PR**2) 

      pOut(0) = E
      pOut(1) = PX
      pOut(2) = PY
      pOut(3) = PZ

    end subroutine phspace2

    !C********************************************************************************
    !C******    frames L and L*, L* moves with velocity vx,vy,vz in L           ******
    !C******    input:         es,pxs,pys,pzs in L*                             ******
    !C******    output:        e, px, py, pz  in L                              ******
    !C******************************************************************************** 
    subroutine boost(VV,PP)

      real, dimension(1:3), intent(in) :: VV
      real, dimension(0:3), intent(inOut) :: PP

      real :: vx,vy,vz,e,px,py,pz,es,pxs,pys,pzs
      real :: gamma, v, vps, effm
      
      vx = VV(1)
      vy = VV(2)
      vz = VV(3)

      es  = PP(0)
      pxs = PP(1)
      pys = PP(2)
      pzs = PP(3)
            
      v = sqrt(vx**2 + vy**2 + vz**2)
      gamma = 1./sqrt(1. - v**2)
      vps = vx*pxs + vy*pys + vz*pzs
      
      if(1..le.v) then
         write (*,*) 'from cdcoll: velocity > = 1'
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

      PP(0) = e
      PP(1) = px
      PP(2) = py
      PP(3) = pz
      
    end subroutine boost



  !*************************************************************************
  end subroutine ToyModel !*************************************************
  !*************************************************************************


end module fragmentationHyp
