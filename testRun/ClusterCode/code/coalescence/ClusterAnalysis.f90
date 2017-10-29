!***************************************************************************
!****m* /ClusterAnalysis
! PURPOSE
! This module contains subroutines for calculation of Flows, 
! distributions, yields, etc.
!***************************************************************************
module ClusterAnalysis

  use InputForAnalysis

  PRIVATE
  !Flows, Rapidities, Charge and Mass distributions, Transverse Momentum Spectra, etc...

  integer, parameter :: nrap = 60, npt = 60, ip = 7
  real, parameter    :: dy = 0.1, yMax=6., dpt = 0.025, pi=3.14159
  integer, parameter :: zval = 200
  real, parameter    :: dz = 1.

  logical, save :: Init_Flow  =.true.
  logical, save :: Init_ZDist =.true.
  logical, save :: InitFlag   =.true.
  !rapidity bin (same in x-,y-,z-direction)
  real, save,dimension(-nrap:nrap)        :: ybin
  !Mean transverse Flow as function of rapidity in beam direction
  real, save,dimension(-nrap:nrap,ip)     :: Flowx
  !rapidity distributions (x,y,z) for clusters selected in Z 
  !and for free protons,d,t,4He,6He,6Li
  real, save,dimension(-nrap:nrap,ip,1:3) :: dndy,dndyfra
  !rapidity distributions (in beam direction) for all nucleons (from ParticleVector)
  real, save,dimension(-nrap:nrap)        :: dndyNuc
  !rapidity distributions (in beam direction) for hyperons,pions and kaons
  real, save,dimension(-nrap:nrap,2)      :: dndyHyp,dndyPion,dndyKaon
  !rapidity distributions (in beam direction) for hyperclusters
  real, save,dimension(-nrap:nrap,ip+2)     :: dndyHypHI
  !transverse momentum bin
  real, save,dimension(npt)         :: ptbin
  !transverse momentum spectra for different particles (according Z)
  real, save,dimension(npt,ip)      :: dndpt
  !charge and mass distributions
  real, save,dimension(zval)        :: zbin,Zdistrb,Zdistr,Zdistrb_Spect,Zdistr_Spect,Adistr,EkinZ
  real, dimension(1:13,zval), SAVE :: Adistr2 !distributions in A for different isotopes

  PUBLIC :: MainAnalysis

contains 

  !***************************************************************************
  !****s* ClusterAnalysis/MainAnalysis
  ! NAME
  ! subroutine MainAnalysis
  !
  ! PURPOSE
  ! Main Routine for the analysis
  !***************************************************************************
  subroutine MainAnalysis(i,j,FragmentVector,ParticleVector,Get_Flow,Get_Zdist, &
                          SubEvents,NumEnsemples)

    use typeDefinitions, only : cluster,particle

    implicit none

    integer, intent(in) :: i,j,SubEvents,NumEnsemples
    type(cluster), dimension(:),intent(in) :: FragmentVector
    type(particle), dimension(:),intent(in) :: ParticleVector
    logical, intent(in) :: Get_Flow,Get_Zdist

    if (initFlag) then
       call GetAnalysisParameters
       initFlag = .false.
    endif

    if (Get_Flow) then
       if (Init_Flow) then 
          call FlowInit
          Init_Flow=.false.
       endif
       call Flows(FragmentVector,ParticleVector)
       if ( i==SubEvents .and. j==NumEnsemples ) then 
          call FlowOutput(SubEvents,NumEnsemples)
       endif
    endif

    if (Get_Zdist) then
       if (Init_ZDist) then 
          call ZdistInit
          Init_ZDist=.false.
       endif
       call ChargeDistribution(FragmentVector)
       if ( i==SubEvents .and. j==NumEnsemples ) then 
          call ZdistOutput(SubEvents,NumEnsemples)
       endif
    endif
!***************************************************************************
  end subroutine MainAnalysis !*********************************************
!***************************************************************************

  !*************************************************************************
  !****s* ClusterAnalysis/Flows
  ! NAME
  ! subroutine Flows
  !
  ! PURPOSE
  ! Calculation of In-Plane Flow, transverse spectra, stopping of various particles
  !*************************************************************************
  subroutine Flows(FragmentVector,ParticleVector)
    !-----------------------------------------------------------------------
    use typeDefinitions, only : cluster,particle
    implicit none

    type(cluster), dimension(:),intent(in) :: FragmentVector
    type(particle), dimension(:),intent(in) :: ParticleVector
    integer :: m,ibin,ibinx,ibiny,ibinz,ibin2,charge,l,ich,mass,HypCon
    real :: k0f,k1f,k2f,k3f,yb,ybx,yby,ybz,pt
    logical :: PartType
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FragmentLoop : do m=1,size(FragmentVector,dim=1)
       if (FragmentVector(m)%ID==0) cycle
       if (FragmentVector(m)%ID .gt. 40) then
          cycle
       endif
       mass     = FragmentVector(m)%MassNumber
       charge   = FragmentVector(m)%ChargeNumber
       k0f      = FragmentVector(m)%momentum(0)
       k1f      = FragmentVector(m)%momentum(1)
       k2f      = FragmentVector(m)%momentum(2)
       k3f      = FragmentVector(m)%momentum(3)
       HypCon   = FragmentVector(m)%HypNumber
       PartType = FragmentVector(m)%FreeBound

!       yb = 0.5*log( (k0f+k3f)/(k0f-k3f) )/yproj

       ybz = 0.5*log( (k0f+k3f)/(k0f-k3f) ) !/yproj
       ybx = 0.5*log( (k0f+k1f)/(k0f-k1f) )  !/yproj
       yby = 0.5*log( (k0f+k2f)/(k0f-k2f) )  !/yproj

       if ( abs(ybz) > abs(yMax) ) then
          if (ybz < 0.) ybz = -yMax
          if (ybz > 0.) ybz = yMax
       endif

       if ( abs(ybx) > abs(yMax) ) then
          if (ybx < 0.) ybx = -yMax
          if (ybx > 0.) ybx = yMax
       endif

       if ( abs(yby) > abs(yMax) ) then
          if (yby < 0.) yby = -yMax
          if (yby > 0.) yby = yMax
       endif

       pt = sqrt(k1f*k1f+k2f*k2f)

!       ibinx = int( ybx/dy + 21. )
!       ibinx = min( nrap, max(1,ibinx) )
!       ibiny = int( yby/dy + 21. )
!       ibiny = min( nrap, max(1,ibiny) )
!       ibinz = int( ybz/dy + 21. )
!       ibinz = min( nrap, max(1,ibinz) )
       if (ybz < 0.) then
          ibinz = int( ybz/dy-1.)
       else
          ibinz = int( ybz/dy)
       endif
       if (ybx < 0.) then
          ibinx = int( ybx/dy-1.)
       else
          ibinx = int( ybx/dy)
       endif
       if (yby < 0.) then
          ibiny = int( yby/dy-1.)
       else
          ibiny = int( yby/dy)
       endif

       if (abs(ibinx) > nrap) then
          write(*,*) 'ClusterAnalysis/Flows:'
          write(*,*) 'Dimension overflow x!!!',ibinx,ybx
          STOP
       endif

       if (abs(ibiny) > nrap) then
          write(*,*) 'ClusterAnalysis/Flows:'
          write(*,*) 'Dimension overflow y!!!',ibiny,yby
          STOP
       endif

       if (abs(ibinz) > nrap) then
          write(*,*) 'ClusterAnalysis/Flows:'
          write(*,*) 'Dimension overflow z!!!',ibinz,ybz
          STOP
       endif

!       ibinz = min( nrap, max(1,ibinz) )
!       ibinx = min( nrap, max(1,ibinz) )
!       ibinz = min( nrap, max(1,ibinz) )

       ibin2 = int( pt/dpt + 1. )
       ibin2 = min( npt, max(1,ibin2) )

       if (ibin2.lt.0.or.ibin2.gt.npt) then
          write(*,*) 'something wrong in Flows...ptbin',m,pt,ibin2
          stop
       endif
       !Transverse Flow & Rapidity Distribution
       !------------------------------------------ all particles
       Flowx(ibinz,1) = Flowx(ibinz,1) + k1f
       dndy(ibinx,1,1) = dndy(ibinx,1,1) + 1.
       dndy(ibiny,1,2) = dndy(ibiny,1,2) + 1.
       dndy(ibinz,1,3) = dndy(ibinz,1,3) + 1.
       !------------------------------------------ charged particles (Z=>1)
       if (charge.ge.1) then
          Flowx(ibinz,2) = Flowx(ibinz,2) + k1f
          dndy(ibinx,2,1) = dndy(ibinx,2,1) + 1.
          dndy(ibiny,2,2) = dndy(ibiny,2,2) + 1.
          dndy(ibinz,2,3) = dndy(ibinz,2,3) + 1.
       endif
       !------------------------------------------ particles with Z=1,2,3,...
       ich = 1                                    !depends on ip
       do l=3,ip
          if (charge==ich.and.PartType) then
             Flowx(ibinz,l) = Flowx(ibinz,l) + k1f
             dndy(ibinx,l,1) = dndy(ibinx,l,1) + 1.
             dndy(ibiny,l,2) = dndy(ibiny,l,2) + 1.
             dndy(ibinz,l,3) = dndy(ibinz,l,3) + 1.
          endif
          ich = ich + 1
       end do
       !------------------------------------------ free protons
       if (charge.eq.1 .and. .not.PartType) then
          dndyfra(ibinx,1,1) = dndyfra(ibinx,1,1) + 1.
          dndyfra(ibiny,1,2) = dndyfra(ibiny,1,2) + 1.
          dndyfra(ibinz,1,3) = dndyfra(ibinz,1,3) + 1.
       endif
       !------------------------------------------ deuterons
       if (charge.eq.1 .and. mass.eq.2) then   
          dndyfra(ibinx,2,1) = dndyfra(ibinx,2,1) + 1.
          dndyfra(ibiny,2,2) = dndyfra(ibiny,2,2) + 1.
          dndyfra(ibinz,2,3) = dndyfra(ibinz,2,3) + 1.
       endif
       !------------------------------------------ He-isotopes
       if (charge.eq.2) then
          if (mass.eq.3) then
             dndyfra(ibinx,3,1) = dndyfra(ibinx,3,1) + 1. !He3
             dndyfra(ibiny,3,2) = dndyfra(ibiny,3,2) + 1. !He3
             dndyfra(ibinz,3,3) = dndyfra(ibinz,3,3) + 1. !He3
          endif
          if (mass.eq.4) then
             dndyfra(ibinx,4,1) = dndyfra(ibinx,4,1) + 1. !He4
             dndyfra(ibiny,4,2) = dndyfra(ibiny,4,2) + 1. !He4
             dndyfra(ibinz,4,3) = dndyfra(ibinz,4,3) + 1. !He4
          endif
          if (mass.eq.5) then
             dndyfra(ibinx,5,1) = dndyfra(ibinx,5,1) + 1. !He5
             dndyfra(ibiny,5,2) = dndyfra(ibiny,5,2) + 1. !He5
             dndyfra(ibinz,5,3) = dndyfra(ibinz,5,3) + 1. !He5
          endif
          if (mass.eq.6) then
             dndyfra(ibinx,6,1) = dndyfra(ibinx,6,1) + 1. !He6
             dndyfra(ibiny,6,2) = dndyfra(ibiny,6,2) + 1. !He6
             dndyfra(ibinz,6,3) = dndyfra(ibinz,6,3) + 1. !He6
          endif
       endif
       !------------------------------------------ Li-isotopes
       if (charge.eq.3) then
          if (mass.eq.6) then
             dndyfra(ibinx,7,1) = dndyfra(ibinx,7,1) + 1. !Li6
             dndyfra(ibiny,7,2) = dndyfra(ibiny,7,2) + 1. !Li6
             dndyfra(ibinz,7,3) = dndyfra(ibinz,7,3) + 1. !Li6
          endif
       endif
       !--------------------------------------------------------------------
       !Transverse Momentum Spectra
       !--------------------------------------------------------------------
       dndpt(ibin2,1) = dndpt(ibin2,1) + 1./pt
       !------------------------------------------ charged particles (Z=>1)
       if (charge.ge.1) then
          dndpt(ibin2,2) = dndpt(ibin2,2) + 1./pt
       endif
       !------------------------------------------ particles with Z=1,2,3,...
       ich = 1                                    !depends on ip
       do l=3,ip
          if (charge==ich.and.PartType) then
             dndpt(ibin2,l) = dndpt(ibin2,l) + 1./pt
          endif
          ich = ich + 1
       end do
       !------------------------------------------ HyperFragments
       if (abs(ybz).le.2.) then

          if (charge==2) then
             if (mass.eq.2.and.HypCon.eq.1.and.& 
                  & FragmentVector(m)%HypType=='L') &  !He3-1Lambda
                  & dndyHypHI(ibinz,1) = dndyHypHI(ibinz,1) + 1.
             if (mass.eq.2.and.HypCon.eq.1.and.& 
                  & FragmentVector(m)%HypType=='S') &  !He3-1Sigma
                  & dndyHypHI(ibinz,2) = dndyHypHI(ibinz,2) + 1.

             if (mass.eq.3.and.HypCon.eq.1.and.& 
                  & FragmentVector(m)%HypType=='L') &  !He4-1Lambda
                  & dndyHypHI(ibinz,3) = dndyHypHI(ibinz,3) + 1.
             if (mass.eq.3.and.HypCon.eq.1.and.& 
                  & FragmentVector(m)%HypType=='S') &  !He4-1Sigma
                  & dndyHypHI(ibinz,4) = dndyHypHI(ibinz,4) + 1.

             if (mass.eq.4.and.HypCon.eq.1.and.& 
                  & FragmentVector(m)%HypType=='L') &  !He5-1Lambda
                  & dndyHypHI(ibinz,5) = dndyHypHI(ibinz,5) + 1.
             if (mass.eq.4.and.HypCon.eq.1.and.& 
                  & FragmentVector(m)%HypType=='S') &  !He5-1Sigma
                  & dndyHypHI(ibinz,6) = dndyHypHI(ibinz,6) + 1.

             if (mass.eq.5.and.HypCon.eq.1.and.& 
                  & FragmentVector(m)%HypType=='L') &  !He6-1Lambda
                  & dndyHypHI(ibinz,7) = dndyHypHI(ibinz,7) + 1.
             if (mass.eq.5.and.HypCon.eq.1.and.& 
                  & FragmentVector(m)%HypType=='S') &  !He6-1Sigma
                  & dndyHypHI(ibinz,8) = dndyHypHI(ibinz,8) + 1.

          endif

       endif

    end do FragmentLoop
    !-----------------------------------------------------------------------
    ! calculate again flows from ParticleVector
    !-----------------------------------------------------------------------
    do m=1,size(ParticleVector,dim=1)
       k0f      = ParticleVector(m)%momentum(0)
       k1f      = ParticleVector(m)%momentum(1)
       k2f      = ParticleVector(m)%momentum(2)
       k3f      = ParticleVector(m)%momentum(3)

!       yb = 0.5*log( (k0f+k3f)/(k0f-k3f) )/yproj

       yb = 0.5*log( (k0f+k3f)/(k0f-k3f) )
       
       if (abs(yb) > abs(yMax) ) cycle

       if (yb < 0.0) then
          ibin = int( yb/dy -1.)
       else
          ibin = int( yb/dy)
       endif
       
!       ibin = min( nrap, max(1,ibin) )       

       if (abs(yb).le.2.) then
          if(ParticleVector(m)%ID.lt.2) dndyNuc(ibin) = dndyNuc(ibin) + 1.
!          if(ParticleVector(m)%ID==32 .or. ParticleVector(m)%ID==33) then
!             dndyHyp(ibin) = dndyHyp(ibin) + 1.
!          endif
          if(ParticleVector(m)%ID==32) dndyHyp(ibin,1) = dndyHyp(ibin,1) + 1. !Lambda
          if(ParticleVector(m)%ID==33) dndyHyp(ibin,2) = dndyHyp(ibin,2) + 1. !Sigma
          if(ParticleVector(m)%ID==101) then !Pions 
             if (ParticleVector(m)%Charge==-1) dndyPion(ibin,1) = dndyPion(ibin,1) + 1. !pi-
             if (ParticleVector(m)%Charge==1)  dndyPion(ibin,2) = dndyPion(ibin,2) + 1. !pi+
          endif
          if(ParticleVector(m)%ID==110) then !Kaons
             if (ParticleVector(m)%Charge==1) dndyKaon(ibin,1) = dndyKaon(ibin,1) + 1. !K+
             if (ParticleVector(m)%Charge==0) dndyKaon(ibin,2) = dndyKaon(ibin,2) + 1. !K0
          endif
       endif
    end do

!***************************************************************************
  end subroutine Flows !****************************************************
!***************************************************************************


  !***************************************************************************
  !****s* ClusterAnalysis/FlowInit
  ! NAME
  ! subroutine FlowInit
  !
  ! PURPOSE
  ! Initialization of Flow-Variables
  !***************************************************************************
  subroutine FlowInit

    implicit none

    integer :: iloc1,iloc2
    real :: ystep

    do iloc1 = -nrap,nrap
       do iloc2 = 1,ip
          Flowx(iloc1,iloc2)       = 0.0
          dndy(iloc1,iloc2,1:3)    = 0.0
          dndyfra(iloc1,iloc2,1:3) = 0.0
          dndyHypHI(iloc1,iloc2)   = 0.0
       end do
       dndyNuc(iloc1)      = 0.0
       dndyHyp(iloc1,1:2)  = 0.0
       dndyPion(iloc1,1:2) = 0.0
       dndyKaon(iloc1,1:2) = 0.0
    end do
    ystep = -yMax
    do iloc1=-nrap,nrap
       ybin(iloc1) = ystep
       ystep = ystep + dy
    end do
    do iloc1=-nrap,nrap
       ybin(iloc1) = ybin(iloc1) + dy/2.
    end do

    do iloc1 = 1,npt
       do iloc2 = 1,ip
          dndpt(iloc1,iloc2) = 0.0
       end do
    end do
    ystep = dpt/2.
    do iloc1=1,npt
       ptbin(iloc1) = float(iloc1)*dpt-ystep
    end do
!***************************************************************************
  end subroutine FlowInit !*************************************************
!***************************************************************************

  !*************************************************************************
  !****s* ClusterAnalysis/FlowOutput
  ! NAME
  ! subroutine FlowOutput
  !
  ! PURPOSE
  ! Output of Flow observables
  !*************************************************************************
  subroutine FlowOutput(SubEvents,NumEnsemples)
    use InputCoalescence, only : Out1,Out2,Out5,Out6,Out7
    implicit none

    integer, intent(in) :: SubEvents,NumEnsemples
    integer m,l

    open(Unit=10, File=Out1)
    open(Unit=11, File=Out2)
    open(Unit=12, File=Out5)
    open(Unit=13, File=Out6)
    open(Unit=14, File=Out7)

    do m=-nrap,nrap
       do l=1,ip
          if (dndy(m,l,3).gt.0.0001)  then
             Flowx(m,l) = Flowx(m,l)/dndy(m,l,3)     
          else
             Flowx(m,l) = 0.0
          endif
          dndy(m,l,1:3)    = dndy(m,l,1:3)/float(SubEvents*NumEnsemples)/dy
          dndyHypHI(m,l)   = dndyHypHI(m,l)/float(SubEvents*NumEnsemples)/dy
          dndyFra(m,l,1:3) = dndyFra(m,l,1:3)/float(SubEvents*NumEnsemples)/dy
       end do
       dndyNuc(m)    = dndyNuc(m)/float(SubEvents*NumEnsemples)/dy
       dndyHyp(m,:)  = dndyHyp(m,:)/float(SubEvents*NumEnsemples)/dy
       dndyKaon(m,:) = dndyKaon(m,:)/float(SubEvents*NumEnsemples)/dy
       dndyPion(m,:) = dndyPion(m,:)/float(SubEvents*NumEnsemples)/dy

       write(10,1000) ybin(m),& 
            &        (dndy(m,l,3),l=1,ip),&
            &        (dndy(m,l,1),l=1,ip),& 
            &        (dndy(m,l,2),l=1,ip)
       write(12,1000) ybin(m),& 
            &        (dndyFra(m,l,3),l=1,ip),&
            &        (dndyFra(m,l,1),l=1,ip),& 
            &        (dndyFra(m,l,2),l=1,ip)
       write(13,1000) ybin(m),& 
            & dndyNuc(m), & 
            & dndyHyp(m,1:2), & 
            & dndyKaon(m,1:2), & 
            & dndyPion(m,1:2)
       write(14,1000) ybin(m),(dndyHypHI(m,l),l=1,ip+2)
    end do
    do m=1,npt
       do l=1,ip
          dndpt(m,l) = dndpt(m,l)/dpt/float(SubEvents*NumEnsemples)
       end do
       write(11,1001) ptbin(m),(dndpt(m,l),l=1,ip)
    end do

    close(Unit=10)
    close(Unit=11)
    close(Unit=12)
    close(Unit=13)
    close(Unit=14)

1000 format(50(1x,e13.6))
1001 format(f8.4,1x,6(f12.4,1x))
!***************************************************************************
  end subroutine FlowOutput !***********************************************
!***************************************************************************

  !*************************************************************************
  !****s* ClusterAnalysis/ChargeDistribution
  ! NAME
  ! subroutine ChargeDistribution
  !
  ! PURPOSE
  ! Calculation of Charge- and Massdistributions of produced clusters
  !*************************************************************************
  subroutine ChargeDistribution(FragmentVector)

    use typeDefinitions, only : cluster
    implicit none

    type(cluster), dimension(:),intent(in) :: FragmentVector
    integer :: m,ibin
    real :: aa,zz,k0f,k3f,yb,kinEn
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ParticleLoop : do m=1,size(FragmentVector,dim=1)
       if (FragmentVector(m)%ID==0) cycle
       if (FragmentVector(m)%ID > 2) cycle !only Nucleons here!!!
!       if (.not.FragmentVector(m)%FreeBound) cycle
!       mass     = FragmentVector(m)%MassNumber
       zz = float(FragmentVector(m)%ChargeNumber)
       aa = float(FragmentVector(m)%MassNumber)
       kinEn = FragmentVector(m)%momentum(0) - FragmentVector(m)%Mass*0.19733 

       ibin = int( aa/dz)
       ibin = min( zval, max(1,ibin) )
       Adistr(ibin) = Adistr(ibin) + 1.

       if (zz==20) Adistr2(1,ibin) = Adistr2(1,ibin) + 1.
       if (zz==25) Adistr2(2,ibin) = Adistr2(2,ibin) + 1.
       if (zz==30) Adistr2(3,ibin) = Adistr2(3,ibin) + 1.
       if (zz==35) Adistr2(4,ibin) = Adistr2(4,ibin) + 1.
       if (zz==40) Adistr2(5,ibin) = Adistr2(5,ibin) + 1.
       if (zz==45) Adistr2(6,ibin) = Adistr2(6,ibin) + 1.
       if (zz==50) Adistr2(7,ibin) = Adistr2(7,ibin) + 1.
       if (zz==55) Adistr2(8,ibin) = Adistr2(8,ibin) + 1.
       if (zz==60) Adistr2(9,ibin) = Adistr2(9,ibin) + 1.
       if (zz==65) Adistr2(10,ibin) = Adistr2(10,ibin) + 1.
       if (zz==70) Adistr2(11,ibin) = Adistr2(11,ibin) + 1.
       if (zz==75) Adistr2(12,ibin) = Adistr2(12,ibin) + 1.
       if (zz==80) Adistr2(13,ibin) = Adistr2(13,ibin) + 1.

       if (zz.eq.0.) cycle

       ibin = int( zz/dz)
       ibin = min( zval, max(1,ibin) )

       if (ibin.lt.0.or.ibin.gt.zval) then
          write(*,*) 'vai a fare bagno...zval',m,zz,ibin
          stop
       endif
       Zdistrb(ibin) = Zdistrb(ibin) + 1.
       Zdistr(ibin) = Zdistr(ibin) + 1.

       EkinZ(ibin) = EkinZ(ibin) + KinEn
       !--------------------------------------------------------------------
       ! Charge distribution of spectator fragments
       !--------------------------------------------------------------------
       if (SpectCut.gt.0.) then
          k0f      = FragmentVector(m)%momentum(0)
          k3f      = FragmentVector(m)%momentum(3)
          yb       = 0.5*log( (k0f+k3f)/(k0f-k3f) )/yproj
          if (abs(yb)>SpectCut) then
             Zdistrb_Spect(ibin) = Zdistrb_Spect(ibin) + 1.
             Zdistr_Spect(ibin) = Zdistr_Spect(ibin) + 1.
          endif
       endif
       !--------------------------------------------------------------------
    end do ParticleLoop
!***************************************************************************
  end subroutine ChargeDistribution !***************************************
!***************************************************************************


  !*************************************************************************
  !****s* ClusterAnalysis/ZDistInit
  ! NAME
  ! subroutine ZDistInit
  !
  ! PURPOSE
  ! Initialization of Charge- and Massdistribution variables
  !*************************************************************************
  subroutine ZdistInit

    implicit none

    integer :: i

    do i = 1,zval
       Zdistrb(i)=0.
       Zdistrb_Spect(i) = 0.0
       Zdistr(i)=0.
       Zdistr_Spect(i) = 0.0
       Adistr(i)=0.
       Adistr2(1:13,i) = 0.0
       EkinZ(i) = 0.
    end do
    do i=1,zval
       zbin(i) = float(i)*dz
    end do
!***************************************************************************
  end subroutine ZdistInit !************************************************
!***************************************************************************

  !*************************************************************************
  !****s* ClusterAnalysis/ZDistOutput
  ! NAME
  ! subroutine ZDistOutput
  !
  ! PURPOSE
  ! Output of Charge- and Mass distributions
  !*************************************************************************
  subroutine ZdistOutput(SubEvents,NumEnsemples)

    use InputCoalescence, only : Out3
    implicit none

    integer, intent(in) :: SubEvents,NumEnsemples
    integer :: m,l,iout
    real    :: Ekinm

    open(Unit=12, File=Out3, IOStat=iout)

    do m=1,zval
       Zdistrb(m) = Zdistrb(m)/float(SubEvents*NumEnsemples)/dz
       Zdistrb_Spect(m) = Zdistrb_Spect(m)/float(SubEvents*NumEnsemples)/dz
       Zdistr(m) = Zdistr(m)/float(SubEvents*NumEnsemples)/dz
       Zdistr_Spect(m) = Zdistr_Spect(m)/float(SubEvents*NumEnsemples)/dz
       Adistr(m) = Adistr(m)/float(SubEvents*NumEnsemples)/dz
       EkinZ(m) = EkinZ(m)/float(SubEvents*NumEnsemples)/dz

       if (Zdistr(m).ne.0.) then
          Ekinm = EkinZ(m)/Zdistr(m)
       else
          Ekinm = 0.
       endif

       write(12,1002) zbin(m),Zdistrb(m),Zdistr(m),Ekinm,Adistr(m), & 
            & (Adistr2(l,m),l=1,13)
    end do

    close(Unit=12)

1002 format(f8.4,4(1x,e12.4),1x,13(1x,e12.4))
!***************************************************************************
  end subroutine ZdistOutput !**********************************************
!***************************************************************************


end module ClusterAnalysis
