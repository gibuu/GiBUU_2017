! Program for testing purposes while improving the mass selection 
! as in rhomass.f90
!
! you may run this via './testmass.x < jobDUMMY'


program testMass
  use inputGeneral
  use mediumDefinition
  use histf90
  use version
  use hadronFormation
  use particleProperties

  implicit none

  integer :: i,n,KF
  real :: m
  type(medium), save :: mediumAtCollision
  type(histogram) :: hM1,hM2,hhM1,hhM2,hBWD,hBWD0,hReject1,hReject2
  real :: nReject1, nReject2

  call PrintVersion
  call readInputGeneral


  call init_Database

  call CreateHist(hM1,"mass (old)",0.,3.,0.01)
  call CreateHist(hM2,"mass (new)",0.,3.,0.01)
  call CreateHist(hhM1,"init mass (old)",0.,3.,0.01)
  call CreateHist(hhM2,"init mass (new)",0.,3.,0.01)
  call CreateHist(hBWD,"BWD",0.,3.,0.01)
  call CreateHist(hBWD0,"BWD0",0.,3.,0.01)
  call CreateHist(hReject1,"calls in rejection (old)",-0.5,101.0,1.0)
  call CreateHist(hReject2,"calls in rejection (new)",-0.5,101.0,1.0)

!  KF = 113
!  KF = 223
!  KF = 333
  KF = -2

  nReject1=0
  nReject2=0
  n=10000000
  do i=1,n
     m = VM_Mass_ORIG(KF)
     call AddHist(hM1,m,1.0)
     m = VM_Mass_NEW(KF)
     call AddHist(hM2,m,1.0)
  end do

  call WriteHist(hM1,121,add=1e-20,mul=1.0/n,file="mass1.dat")
  call WriteHist(hM2,121,add=1e-20,mul=1.0/n,file="mass2.dat")

  call WriteHist(hhM1,121,add=1e-20,mul=1.0/(10.*n),file="mass1a.dat")
  call WriteHist(hhM2,121,add=1e-20,mul=1.0/(10.*n),file="mass2a.dat")

  call WriteHist(hBWD,121,add=1e-20,DoAve=.true.,file="bwd.dat")
  call WriteHist(hBWD0,121,add=1e-20,DoAve=.true.,file="bwd0.dat")

  call WriteHist(hReject1,121,add=1e-20,mul=1.0/n,file="Reject1.dat")
  call WriteHist(hReject2,121,add=1e-20,mul=1.0/n,file="Reject2.dat")

  write(*,*) "nReject: ",nReject1/n,nReject2/n
  write(112,*) "nReject: ",nReject1/n,nReject2/n
  
contains
  
  real function VM_Mass_ORIG(kfa)
    use particleProperties, only : meson,baryon
    use IDTable, only : rho,omegaMeson,phi,Delta
    use random
    use mesonWidthMedium, only: WidthMesonMedium_Const,get_MediumSwitchMesons
    use mesonPotentialModule, only: vecMes_massShift
    use baryonWidth, only: FullWidthBaryon
    use output
    integer,intent(in):: kfa
    real:: m,minmass,maxmass,bwd,m2,ymin,ymax,y,mres0,gamres0,gamtot,spectral,intfac,maxbwd
    real:: gamtot_max,spectral_max,intfac_max
    integer:: ncount,vm
    integer, parameter:: ncountm = 100
    logical:: lcount

    ! decode KF code into GiBUU ID
    select case (kfa)
    case(113,213)
       vm=rho
    case (223)
       vm=omegaMeson
    case(333)
       vm=phi
    case(-2)
       vm= -2 ! == DELTA
    case default
       write (*,*) "Error in VM_Mass_ORIG: KFA = ",kfa
       stop
    end select
    
    if (vm>0) then
       minMass = meson(vm)%minmass
       
       mres0=meson(vm)%mass
       gamres0=WidthMesonMedium_Const(vm,mres0,mediumAtCollision)
    else
       minMass = baryon(Delta)%minmass
       
       mres0=baryon(Delta)%mass
       gamres0=FullWidthBaryon(Delta,mres0)
    endif
    
    maxmass = 3.0
    
    ymax=2.*atan((maxmass-mres0)/gamres0*2.)
    ymin=2.*atan((minmass-mres0)/gamres0*2.)

    ! assume increasing width with increasing mass:
    if (vm>0) then
       gamtot_max   = WidthMesonMedium_Const(vm,maxmass,mediumAtCollision)
    else
       gamtot_max   = FullWidthBaryon(Delta,maxmass)
    endif
    spectral_max = maxmass**2*gamtot_max*gamres0/((mres0**2-maxmass**2)**2 + gamtot_max**2*maxmass**2)
    intfac_max   = gamres0**2/4./((maxmass-mres0)**2 + gamres0**2/4.)
    maxbwd = spectral_max/intfac_max * 1.05


    do ncount=1,10
       y=ymin+rn()*(ymax-ymin)
       m=.5*tan(y/2.)*gamres0+mres0
       m=min(max(m,minmass),maxmass)
       m2=m**2

       ! for test purposes:
       if (vm>0) then
          gamtot=WidthMesonMedium_Const(vm,maxmass,mediumAtCollision)
       else
          gamtot=FullWidthBaryon(Delta,maxmass)
       end if
       spectral=m2*gamtot*gamres0/((mres0**2-m2)**2 + gamtot**2*m2)
       intfac=gamres0**2/4./((m-mres0)**2 + gamres0**2/4.)
       bwd=spectral/intfac/maxbwd
!       call AddHist(hBWD0,m,1.0,bwd)


       if (vm>0) then
          gamtot=WidthMesonMedium_Const(vm,m,mediumAtCollision)
       else
          gamtot=FullWidthBaryon(Delta,m)
       end if
       spectral=m2*gamtot*gamres0/((mres0**2-m2)**2 + gamtot**2*m2)
       intfac=gamres0**2/4./((m-mres0)**2 + gamres0**2/4.)
       bwd=spectral/intfac/maxbwd
       call AddHist(hhM1,m,1.0)
       call AddHist(hBWD0,m,1.0,bwd)
!       call AddHist(hBWD,m,1.0,bwd)
    end do

    ncount=0
    lcount=.true.
    do while(lcount)
       ncount=ncount+1
       y=ymin+rn()*(ymax-ymin)
       m=.5*tan(y/2.)*gamres0+mres0
       m=min(max(m,minmass),maxmass)
       m2=m**2

       if (vm>0) then
          gamtot=WidthMesonMedium_Const(vm,m,mediumAtCollision)
       else
          gamtot=FullWidthBaryon(Delta,m)
       end if

       spectral=m2*gamtot*gamres0/((mres0**2-m2)**2 + gamtot**2*m2)
       intfac=gamres0**2/4./((m-mres0)**2 + gamres0**2/4.)
       bwd=spectral/intfac/maxbwd

       if(bwd>1) then
          if (DoPr(1)) then
             write(*,'(A,i3,A,2F12.3,A,2G12.3)') 'Problems in VM_Mass_ORIG: (bwd>1)  vm=',vm, &
                  & ' mass=',minmass,maxmass,' bwd:',bwd,maxbwd
             write(*,'(A,1P,10g12.3)')'If this problem occurs frequently, adjust "maxbwd" to: ',maxbwd*bwd
          end if
       end if

       if(rn()<=bwd) lcount =.false.

       if (.not.lcount) then
          nReject1=nReject1+ncount
          call AddHist(hReject1,float(ncount),1.0)
       end if

       if(ncount>=ncountm) then
          if (DoPr(1)) then
             write(*,'(A,i3,A,2G12.3)') 'Problems in VM_Mass_ORIG: (ncount)  vm=',vm, &
                  & ' mass=',minmass,maxmass
          endif
          m=mres0
          lcount=.false.
       end if
    end do
    if (.not.(m>=minmass .and. m<=maxmass)) then 
        write(*,'(A,i3,A,2G12.3)') 'Problems in VM_Mass_ORIG: (STOP)  vm=',vm, &
                  & ' mass=',minmass,maxmass
       stop
    end if
    VM_Mass_ORIG = m
  end function VM_Mass_ORIG


real function VM_Mass_NEW(kfa)
    use particleProperties, only : meson,baryon
    use IDTable, only : rho,omegaMeson,phi,Delta
    use random
    use mesonWidthMedium, only: WidthMesonMedium_Const,get_MediumSwitchMesons
    use mesonPotentialModule, only: vecMes_massShift
    use baryonWidth, only: FullWidthBaryon
    use output
    use monteCarlo, only: MonteCarloChoose
    integer,intent(in):: kfa
    real:: m,minmass,maxmass,bwd,m2,ymin,ymax,y,mres0,gamres0,gamtot,spectral,intfac,maxbwd
    real:: mres02,gamres024
    real:: gamtot_max,spectral_max,intfac_max
    integer:: ncount,vm
    integer, parameter:: ncountm = 100
    logical:: lcount

!    integer, parameter :: nB = 7
    integer, parameter :: nB = 8
!    integer, parameter :: nB = 10
    real, dimension(nB) :: mbound,bweight,bmaxbwd,mbound2
    integer :: i,iB
    real :: bweightTot


    ! decode KF code into GiBUU ID
    select case (kfa)
    case(113,213)
       vm=rho
    case (223)
       vm=omegaMeson
    case(333)
       vm=phi
    case(-2)
       vm= -2 ! == DELTA
    case default
       write (*,*) "Error in VM_Mass_NEW: KFA = ",kfa
       stop
    end select

!!!!! ATTENTION !!!
!!!!! In the original routine, minmass is somewhat more complicated !!!

    if (vm>0) then
       minMass = meson(vm)%minmass

       mres0=meson(vm)%mass
       gamres0=WidthMesonMedium_Const(vm,mres0,mediumAtCollision)
    else
       minMass = baryon(Delta)%minmass

       mres0=baryon(Delta)%mass
       gamres0=FullWidthBaryon(Delta,mres0)
    endif

    mres02   = mres0**2
    gamres024= gamres0**2/4

    
    maxmass = 3.0

    ! Divide the mass region into slices:

!    mbound = (/minMass,mres0-2*gamres0,mres0-gamres0,mres0,mres0+gamres0,mres0+2*gamres0,maxMass/)

    mbound = (/minMass,mres0-2*gamres0,mres0-gamres0,mres0-0.25*gamres0,mres0,&
         & mres0+gamres0,mres0+2*gamres0,maxMass/)

!    mbound = (/minMass,mres0-3*gamres0,mres0-2*gamres0,mres0-gamres0,mres0-0.25*gamres0,mres0,&
!         & mres0+gamres0,mres0+2*gamres0,mres0+3*gamres0,maxMass/)

    mbound2 = mbound**2

    do i=1,nB
       if (mbound(i).lt.minMass) mbound(i)=minMass
       if (mbound(i).gt.maxMass) mbound(i)=maxMass
       bweight(i) = atan2(2*(mbound(i)-mres0),gamres0)
    end do

    ! Calculate the integral of the Cauchy distribution:

    do i=1,nB-1
       bweight(i)=(bweight(i+1)-bweight(i))
    enddo
    bweight(nB) = 0
    
    ! Calculate the maximal weight for every slice

    do i=1,nB
       if (vm>0) then
          gamtot_max   = WidthMesonMedium_Const(vm,mbound(i),mediumAtCollision)
       else
          gamtot_max   = FullWidthBaryon(Delta,mbound(i))
       endif
       spectral_max = mbound2(i)*gamtot_max*gamres0/((mres02-mbound2(i))**2 + gamtot_max**2*mbound2(i))
       intfac_max   = gamres024/((mbound(i)-mres0)**2 + gamres024)
       bmaxbwd(i)   = spectral_max/intfac_max * 1.05
    end do

    ! maxbwd may have a dip around the pole mass, therefore we take
    ! the maximum value of the two boundaries

    do i=1,nB-1
       bmaxbwd(i) = max(bmaxbwd(i),bmaxbwd(i+1)) 
       bweight(i) = bweight(i) * bmaxbwd(i) ! modify the Cauchy weights
    end do

    ! Now do the actual MC:

    ncount=0
    lcount=.true.
    do while(lcount)
       ncount=ncount+1

       ! STEP 1: Select the mass region:
       ! (it is important to have this step in the overall loop, since
       ! this cures inaccuracies in the calculation of the weights)

       iB=MonteCarloChoose(bweight(1:nB-1), bweightTot)
       ymax=2.*atan2(2*(mbound(iB+1)-mres0),gamres0)
       ymin=2.*atan2(2*(mbound(iB)-mres0),gamres0)
       maxbwd = bmaxbwd(iB)
       
!!$       do i=1,nB
!!$          write(*,*) i,mbound(i),bweight(i)/bweightTot,bmaxbwd(i)
!!$       end do
!!$       stop
       
       ! only for writing out:
       do i=1,10
          y=ymin+rn()*(ymax-ymin)
          m=.5*tan(y/2.)*gamres0+mres0
          m=min(max(m,minmass),maxmass)
          m2=m**2
          
          if (vm>0) then
             gamtot=WidthMesonMedium_Const(vm,m,mediumAtCollision)
          else
             gamtot=FullWidthBaryon(Delta,m)
          end if
          spectral=m2*gamtot*gamres0/((mres02-m2)**2 + gamtot**2*m2)
          intfac=gamres024/((m-mres0)**2 + gamres024)
          bwd=spectral/intfac/maxbwd
          call AddHist(hhM2,m,1.0)
          call AddHist(hBWD,m,1.0,bwd)
       end do

       ! STEP 2: generate random value according Cauchy with constant width 

       y=ymin+rn()*(ymax-ymin)
       m=.5*tan(y/2.)*gamres0+mres0
       m=min(max(m,minmass),maxmass)
       m2=m**2

       ! STEP 3: Do the rejection

       if (vm>0) then
          gamtot=WidthMesonMedium_Const(vm,m,mediumAtCollision)
       else
          gamtot=FullWidthBaryon(Delta,m)
       end if
       spectral=m2*gamtot*gamres0/((mres02-m2)**2 + gamtot**2*m2)
       intfac=gamres024/((m-mres0)**2 + gamres024)
       bwd=spectral/intfac/maxbwd

       if(bwd>1) then
          if (DoPr(1)) then
             write(*,'(A,i3,A,2F12.3,A,2G12.3)') 'Problems in VM_Mass_NEW: (bwd>1)  vm=',vm, &
                  & ' mass=',mbound(iB),mbound(iB+1),' bwd:',bwd,maxbwd
             write(*,'(A,1P,10g12.3)')'If this problem occurs frequently, adjust "maxbwd" to: ',maxbwd*bwd
          end if
       end if

       if(rn()<=bwd) lcount =.false.

       if (.not.lcount) then
          nReject2=nReject2+ncount
          call AddHist(hReject2,float(ncount),1.0)
       end if

       if(ncount>=ncountm) then
          if (DoPr(1)) then
             write(*,'(A,i3,A,2G12.3)') 'Problems in VM_Mass_NEW: (ncount)  vm=',vm, &
                  & ' mass=',mbound(iB),mbound(iB+1)
          endif
          m=mres0
          lcount=.false.
       end if
    end do
    if (.not.(m>=minmass .and. m<=maxmass)) then 
        write(*,'(A,i3,A,2G12.3)') 'Problems in VM_Mass_NEW: (STOP)  vm=',vm, &
                  & ' mass=',minmass,maxmass
       stop
    end if
    VM_Mass_NEW = m
  end function VM_Mass_NEW
end program testMass
