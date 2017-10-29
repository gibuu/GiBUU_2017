!***************************************************************************
!****m* /CoalescenceModule
! NAME
! program CoalescenceModule
!
! PURPOSE
! This is the main program for fragment formation according 
! phase-space coalescence
!***************************************************************************
module CoalescenceModule

  PRIVATE

  PUBLIC :: main_coalescence

  contains

    !*************************************************************************
    subroutine main_coalescence

      use typeDefinitions,  only : cluster,particle !TypeDefinitions
      use InputGeneral,     only : SubEvents,NumEnsemples,RealParticles,& 
           &                       pathToBUUInput,BUU_DataFile,& 
           &                       Get_Hyp, & 
           &                       ALADIN_Flag
      use InputCoalescence
      use coala,           only : Coalescence  !Coalescence between Nucleons
      use HypCoalescence,  only : HypClusters  !Coalescence between Clusters & Hyperons
      use ClusterAnalysis, only : MainAnalysis !Analysis tools (Momentum distr.)
      use Output                               !Output of FragmentVector
      use InitGEM, only : Init_GEM,SetDecays,End_GEM     !Init- and End-routines for GEM
      use WriteStatus, only : IOControl,WriteHypInfo 

      implicit none

      !*************************************************************************
      !****g* main/ParticleVector
      ! SOURCE
      !
      type(particle),allocatable,dimension(:),save :: ParticleVector
      !
      ! PURPOSE
      ! The particle vector from GiBUU
      !
      !*************************************************************************

      !*************************************************************************
      !****g* main/FragmentVector
      ! SOURCE
      !
      type(cluster),allocatable,dimension(:),save :: FragmentVector
      !
      ! PURPOSE
      ! The vector of produced clusters
      ! NOTES
      ! Also free particles (baryons and mesons) are stored into this vector
      !
      !*************************************************************************

      !*************************************************************************
      !****g* main/impactParameter
      ! SOURCE
      !
      Real, SAVE :: impactParameter=1.
      !
      ! PURPOSE
      ! impact parameter for each subsequent run.
      ! NOTES
      ! Important only in ClusterAnalysis module and in case that impact parameter 
      ! has been randomly selected in the BUU runs.
      !
      !*************************************************************************

      integer :: Clusters,global
!      integer :: ios,ior,ioe,idp,iso,subev,status,time,it,SourceIndex
      integer :: ios,ior,ioe,idp,iso,subev,status !,SourceIndex
      integer :: ActualEvent,i,j,NVALS,SpectatorPart
      integer :: nucleons,Hyperons,Mesons
      integer :: Lambdas,Sigmas,Pions(-1:1),Kaonsp,Kaons0
      integer :: GlobalID, collHis
      real    :: prodTime,lastCollTime
      real    :: mass,x,y,z,p0,px,py,pz !,rrr,rand
      logical :: FormHyp
      !------------------------------------------------------------------------*
1004  format(78('#'),/,10('#'),' SubEvents    = ',i5,/, & 
           &           10('#'),' NumEnsemples = ',i5,/, & 
           &           10('#'),' Impact Param = ',f7.2,/,78('#'))

      !-------------------------------------------------------------------------*
      !READ PARAMETERS FROM INPUT & INITIALIZE PROPERTIES OF ELEMENTS ----------*
      !-------------------------------------------------------------------------*
      call Get_InputCoalescence !general Input, contains NAMELIST for jobCard
      call Init_GEM     !Initialize nuclide elements with their properties, 
      !e.g. Spin, MassExcess,BindingEnergy.
      !Needed also for asymmetry check of produced clusters
      if (Get_GEM) then 
         call SetDecays !initialize decay channels when evaporation is switched on
      endif
      !-------------------------------------------------------------------------*
      !READ BUU-DATA FILE AND DETERMINES THE DIMENSION OF THE 
      !PARTICLE- AND FRAGMENT-VECTOR FIELDS 
      !-------------------------------------------------------------------------*

      open(Unit=1, File=trim(PathToBUUInput)// BUU_DataFile, & 
           & Status='old', Action='read',Iostat=ios)

      Open_BUU : if (ios /= 0) then
         write(*,*) 'from InputBUU: BUU_DataFile Open failed: ios = ',ios
         write(*,*) 'from InputBUU: !!! Termination of program NOW !!!'
         STOP
      else Open_BUU
         read(1,*)
         read(1,*)
         SUBRUNS : do i=1,SubEvents
            ENSEMPLES : do j=1,NumEnsemples
               allocate(ParticleVector(RealParticles),STAT=status)
               call IOControl(1,status,'main','ParticleVector')
               nvals       = 0
               SpectatorPart = 0
               nucleons    = 0
               Hyperons    = 0
               Lambdas     = 0
               Sigmas      = 0
               Mesons      = 0
               Pions(-1:1) = 0
               Kaonsp      = 0
               Kaons0      = 0

!               it   = time()
!               rrr  = rand(it)

               Particles : do 
!                     read(1,*,iostat=ior) idp,SourceIndex,iso,mass,x,y,z,px,py,pz,ActualEvent,subev

                  read(200,*,iostat=ior) idp
                  backspace(Unit=200)

                  if (idp < 1000) then
                     read(200,*,iostat=ior) idp,iso,mass,x,y,z,px,py,pz,ActualEvent,subev
                     GlobalID     = 0
                     prodTime     = -999.
                     lastCollTime = -999.
                     collHis      = -999
                  else
                     read(200,*,iostat=ior) GlobalID,idp,prodTime,lastCollTime,collHis, & 
                          & iso,mass,x,y,z,px,py,pz,ActualEvent,subev
                  end if

                  call ReadStatus(ior,i,j,nvals)
                  if (ActualEvent /= j) then
                     !              if (subev /= i) then
                     backspace(Unit=1) !move one line back!!!!
                     exit
                  endif
                  if (ior<0) exit !end of file, exit in any case...
                  if(idp.gt.33.and.idp.lt.100) cycle
                  nvals = nvals + 1
                  p0    = sqrt(mass**2+px**2+py**2+pz**2)
                  ParticleVector(nvals)%number      = globalID
                  ParticleVector(nvals)%bornTime    = prodTime
                  ParticleVector(nvals)%lastCollTime= lastCollTime
                  ParticleVector(nvals)%collHis     = collHis
                  ParticleVector(nvals)%id          = idp
                  ParticleVector(nvals)%event       = subev
                  ParticleVector(nvals)%ensemple    = ActualEvent
                  ParticleVector(nvals)%Charge      = iso
                  ParticleVector(nvals)%mass        = mass/0.19733 ![1/fm]
                  ParticleVector(nvals)%position(1) = x ![fm]
                  ParticleVector(nvals)%position(2) = y ![fm]
                  ParticleVector(nvals)%position(3) = z ![fm]
                  ParticleVector(nvals)%momentum(0) = p0 ![GeV]
                  ParticleVector(nvals)%momentum(1) = px ![GeV]
                  ParticleVector(nvals)%momentum(2) = py ![GeV]
                  ParticleVector(nvals)%momentum(3) = pz ![GeV]
                  
                  if (ALADIN_Flag) then
                     SpectatorPart = SpectatorPart + 1
                  endif

                  if (idp.lt.3) nucleons = nucleons + 1
                  if ( idp==32 .or. idp==33 ) then
                     Hyperons = Hyperons + 1
                     if (iso==0) then
                        Lambdas = Lambdas + 1 !\Lambda^{0}
                     else
                        Sigmas = Sigmas + 1   !\Sigma^{+,-,0}
                     endif
                  endif
                  if (idp.gt.100) then
                     Mesons = Mesons + 1
                     if (idp==101) Pions(iso) = Pions(iso) + 1 !\pi^{+,-,0}
                     if (idp==110) then
                        if (iso==0) Kaons0 = Kaons0 + 1 !K^{0}
                        if (iso==1) Kaonsp = Kaonsp + 1 !K^{+}
                     endif
                  endif
               end do Particles

               if(nvals.gt.RealParticles) then
                  write(*,*) 'Dim. Overflow!!! STOP!!!, nvals,RealParticles = ',nvals,RealParticles
                  STOP
               endif

               if (nvals .le. 0) then
                  write(*,*) 'Empty particle vector???...nvals,i = ',nvals,i
                  STOP
               endif

               allocate(FragmentVector(nvals),STAT=status)
               call IOControl(1,status,'main','FragmentVector')
               !----------------------------------------------------------------
               ! Performe Coalescence+Evaporation between Nucleons
               !----------------------------------------------------------------
               write(*,1004) i,j,ImpactParameter
               write(*,*)
               write(*,'(A,i5)')  'ParticleVevtor            = ',nvals
               write(*,'(A,i5)')  'Nucleons                  = ',nucleons
               if (ALADIN_Flag) & 
                    & write(*,'(A,i5)')  'Projectile nucleons= ',SpectatorPart
               write(*,'(A,3i5)') 'Hyperons,Lambdas,Sigmas   = ',Hyperons,Lambdas,Sigmas
               write(*,'(A,i5)')  'Mesons (Sum)              = ',Mesons
               write(*,'(A,3i5)') 'Non-Strange (\pi^{-,+,0}) = ',Pions(-1:1)
               write(*,'(A,2i5)') 'with strangeness: K^{0,p} = ',Kaons0,Kaonsp
               call Coalescence(ParticleVector,nvals,SpectatorPart,Hyperons,Mesons,R_c,P_c,& 
                    &           MaxGEM,FragmentVector, &
                    &           Clusters,global,Get_Hyp,Get_Asy,Get_GEM,ALADIN_Flag)
               !----------------------------------------------------------------
               ! Performe Coalescence between Clusters & Hyperons
               !----------------------------------------------------------------
               if (Get_Hyp.and.Clusters.ne.0.and.Hyperons.ne.0) then 
                  call HypClusters(ParticleVector,nvals,Hyperons,Clusters,Global, & 
                       &           R_c,P_c,FragmentVector,FormHyp)
               endif
               call WriteHypInfo(i,j,SubEvents,NumEnsemples,Clusters, & 
                    & Hyperons,Mesons,Lambdas,Sigmas,Pions,Kaonsp,Kaons0, & 
                    & FragmentVector)
               !----------------------------------------------------------------
               !performe analysis (charge-distributions, flows)
               !----------------------------------------------------------------
               if (Get_Analysis) &
                  call MainAnalysis(i,j,FragmentVector,ParticleVector,Get_Flow,Get_Zdist, &
                                    SubEvents,NumEnsemples)
               !----------------------------------------------------------------
               if (Get_Output) call WriteGlobalVector(i,FragmentVector)
               !----------------------------------------------------------------
               deallocate(ParticleVector,STAT=status)
               call IOControl(2,status,'main','ParticleVector')
               deallocate(FragmentVector,STAT=status)
               call IOControl(2,status,'main','FragmentVector')
               !----------------------------------------------------------------
            end do ENSEMPLES
         end do SUBRUNS
      endif Open_BUU
      !-------------------------------------------------------------------------*
      if (Get_GEM) call End_GEM !deallocates fields related to GEM-module
      !-------------------------------------------------------------------------*
      close(unit=1,status='keep',iostat=ioe)
      if(ioe /= 0) then
         write(*,*) 'from InputBUU: cannot close BUU_DataFile, IOSTAT = ',ioe
         write(*,*) 'from InputBUU: !!! Termination of program NOW !!!'
         STOP
      endif
      !-------------------------------------------------------------------------*
    !****************************************************************************
    end subroutine main_coalescence !********************************************
    !****************************************************************************

    !****************************************************************************
    subroutine ReadStatus(ior,i,j,nvals)
    !****************************************************************************
      implicit none
      integer, intent(in) :: ior,i,j,nvals

      if (ior > 0) then
         write(*,*) 'from ReadStatus : debug during reading of BUU file: IOSTAT = ',ior, &
              i,j,nvals
         write(*,*) '!!! Termination of program !!!'
         STOP
      endif
    !****************************************************************************
    end subroutine ReadStatus !**************************************************
    !****************************************************************************

  end module CoalescenceModule
