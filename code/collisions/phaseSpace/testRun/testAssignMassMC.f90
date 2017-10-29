program testAssignMassMC

  use inputGeneral
  use mediumDefinition
  use histf90
  use version
  use hadronFormation
  use particleProperties
  use PythiaSpecFunc
  use mediumModule, only: mediumAt
  use output
  use AssignMassMC
  use constants, only: rhoNull

  implicit none

  integer :: iMC,iID
!  integer,parameter :: nMC = 10000000
  integer,parameter :: nMC = 1000000
!  integer,parameter :: nMC = 100

  real :: m
  type(medium), save :: mediumAtPos
  type(histogram) :: hM1,hM2,hM2b,hReject1,hReject2
  integer :: nReject
  real :: srts

  real, dimension(0:3), save :: momLRF = (/0,0,0,0/)

  integer,dimension(1:5), parameter :: IDarr = (/103,105,107,112,120/)
  integer,dimension(1:5), parameter :: KFarr = (/113,223,333, -1, -1/)

  real :: rDens

  call PrintVersion
  call readInputGeneral

  call initParticleProperties

  call CreateHist(hM1,"mass (old)",0.,5.,0.01)
  call CreateHist(hM2,"mass (new)",0.,5.,0.01)
  call CreateHist(hM2b,"massB (new)",0.,5.,0.01)
  call CreateHist(hReject1,"calls in rejection (old)",-0.5,101.0,1.0)
  call CreateHist(hReject2,"calls in rejection (new)",-0.5,101.0,1.0)

  mediumAtPos = mediumAt( (/1000.0, 0.0, 0.0/) ) ! somewhere far off

  rDens = 0.17 ! DUMMY

  mediumAtPos%density = rhoNull*rDens
  mediumAtPos%densityProton =  0.5 * rhoNull*rDens
  mediumAtPos%densityNeutron = 0.5 * rhoNull*rDens
  mediumAtPos%useMedium=.true.

  srts = 10.0 ! some dummy value
  srts = 3.0 ! some dummy value

  call Init_VM_Mass(srts)

  call checkDelta
!  call testAssign2
!  call illustrateQviolation_meson
!  call checkQviolation_meson
!  call checkMeson
!  call checkBaryon

contains
  !*************************************************************************
  subroutine testAssign2

    integer, dimension(2) :: ID2
    real, dimension(2) :: mass2

!     integer :: iS


    ID2 = (/2,103/)
!    ID2 = (/2,101/)
    srts= 10.0
!    srts= 1.5
    srts= 1.0
    iID = 0

!!$    call AssignMass_2(ID2, srts,momLRF,mediumAtPos, mass2, nReject)
!!$    write(*,*) 'ID (new) =',ID2
!!$
!!$    do iS=200,120,-2
!!$       srts=iS*0.01
!!$       call AssignMass_2(ID2, srts,momLRF,mediumAtPos, mass2, nReject)
!!$       call timeMeasurement(ForceReset=.true.)
!!$       do iMC=1,nMC
!!$          call AssignMass_2(ID2, srts,momLRF,mediumAtPos, mass2, nReject)
!!$       end do
!!$       write(*,*) 'sqrts=',srts
!!$       call timeMeasurement()
!!$    end do
!!$
!!$    stop

    call AssignMass_2(ID2, srts,momLRF,mediumAtPos, mass2, nReject)
    write(*,*) 'ID (new) =',ID2
    call timeMeasurement(ForceReset=.true.)
       
    do iMC=1,nMC
       call AssignMass_2(ID2, srts,momLRF,mediumAtPos, mass2, nReject)
       call AddHist(hM2,mass2(1),1.0)
       call AddHist(hM2b,mass2(2),1.0)
       call AddHist(hReject2,real(nReject),1.0)
    end do
    call timeMeasurement()
       
    call WriteHist(hM2,121,add=1e-20,mul=1.0/nMC,file="mass2a."//Achar(iID+48)//".dat")
    call WriteHist(hM2b,121,add=1e-20,mul=1.0/nMC,file="mass2b."//Achar(iID+48)//".dat")
    call WriteHist(hReject2,121,add=1e-20,mul=1.0/nMC,file="reject2."//Achar(iID+48)//".dat")
       
    call ClearHist(hM2)
    call ClearHist(hM2b)
    call ClearHist(hReject2)

  end subroutine testAssign2
  
  !*************************************************************************
  subroutine illustrateQviolation_meson
    use constants, only: rhoNull
    use random
    use MassAssInfoDefinition
    use monteCarlo, only: MonteCarloChoose
    use mesonWidthMedium, only : GetMassAssInfo_Meson,WidthMesonMedium
    use baryonWidthMedium, only : GetMassAssInfo_Baryon,WidthBaryonMedium

    implicit none

    integer,parameter :: nMCdens=10
    integer :: iMCdens
    real :: rDens
    type(tMassAssInfo),save :: MAI,MAI1,MAI2
    integer :: iTry
    integer,parameter :: nTry = 3000

    integer :: iB
    real :: gamma, bwd, mass
    integer :: ID

    iiD = 3
    ID = IDarr(iID)

    do iMCdens=0,nMCdens
       !       rDens = rn()
       rDens = real(iMCdens)/nMCdens+0.1*rn()

       write(*,*) '=== rDens = ',rDens,iMCdens

       mediumAtPos%density = rhoNull*rDens
       mediumAtPos%densityProton =  0.5 * rhoNull*rDens
       mediumAtPos%densityNeutron = 0.5 * rhoNull*rDens
       mediumAtPos%useMedium=.true.

       call GetMassAssInfo_Meson(MAI1,ID,momLRF,mediumAtPos,forceMix=0.01)
       call GetMassAssInfo_Meson(MAI2,ID,momLRF,mediumAtPos,forceMix=0.99)
       call GetMassAssInfo_Meson(MAI, ID,momLRF,mediumAtPos)

       do iB=1,size(MAI%W)
          if (MAI%W(iB)==0) cycle
          
          write(1200,*) MAI%M(IB),MAI%Q(iB)
          write(1200,*) MAI%M(IB+1),MAI%Q(iB)

          write(1201,*) MAI1%M(IB),MAI1%Q(iB)
          write(1201,*) MAI1%M(IB+1),MAI1%Q(iB)

          write(1202,*) MAI2%M(IB),MAI2%Q(iB)
          write(1202,*) MAI2%M(IB+1),MAI2%Q(iB)

       end do

       do iTry =1,nTry
          
          mass = iTry*3.0/nTry
          gamma = WidthMesonMedium(ID,mass,momLRF,mediumAtPos)
          bwd = MassAssInfoQ(MAI%Mass0,MAI%Gamma0,mass,gamma)

          write(1300,*) mass,bwd
          
       end do
       
       write(1200,*) 
       write(1200,*) 
       write(1201,*) 
       write(1201,*) 
       write(1202,*) 
       write(1202,*) 
       write(1300,*) 
       write(1300,*) 
       

    end do
  end subroutine illustrateQviolation_meson
  !*************************************************************************

  !*************************************************************************
  subroutine checkQviolation_meson
    use constants, only: rhoNull
    use random
    use MassAssInfoDefinition
    use monteCarlo, only: MonteCarloChoose
    use mesonWidthMedium, only : GetMassAssInfo_Meson,WidthMesonMedium
    use baryonWidthMedium, only : GetMassAssInfo_Baryon,WidthBaryonMedium

    implicit none

    integer,parameter :: nMCdens= 1000
    integer,parameter :: nTry   = 10000
    integer :: iMCdens
    integer :: iTry
    real :: rDens
    type(tMassAssInfo),save :: MassAssInfo

    integer :: iB
    real :: BinWeightTot
    real :: y, ymin, ymax, gamma, bwd, maxbwd, mass
    integer :: ID

    iiD = 3
    ID = IDarr(iID)

    do iMCdens=1,nMCdens
       rDens = rn()

       mediumAtPos%density = rhoNull*rDens
       mediumAtPos%densityProton =  0.5 * rhoNull*rDens
       mediumAtPos%densityNeutron = 0.5 * rhoNull*rDens
       mediumAtPos%useMedium=.true.

       call GetMassAssInfo_Meson(MassAssInfo,ID,momLRF,mediumAtPos)

       MassAssInfo%W = MassAssInfo%W*MassAssInfo%Q

       do iTry =1,nTry
          
          iB = MonteCarloChoose(MassAssInfo%W,BinWeightTot)

          ymin = MassAssInfo%Y(iB)
          ymax = MassAssInfo%Y(iB+1)
          maxbwd = MassAssInfo%Q(iB)

          ! STEP 2: generate random value according Cauchy with constant width 
       
          y = ymin + rn()*(ymax-ymin)
          mass = 0.5*tan(0.5*y)*MassAssInfo%Gamma0 + MassAssInfo%Mass0
       
          ! STEP 3: Do the rejection

          gamma = WidthMesonMedium(ID,mass,momLRF,mediumAtPos)
          
          bwd = MassAssInfoQ(MassAssInfo%Mass0,MassAssInfo%Gamma0,mass,gamma)

          if (bwd>maxbwd) then 
             write(131,'(20f13.6)') rDens,mass,bwd/maxbwd
!             write(*,'("viol ",20f13.6)') rDens,mass,bwd/maxbwd
          end if
       end do

    end do

  end subroutine checkQviolation_meson
  !*************************************************************************

  !*************************************************************************
  subroutine checkMeson

    implicit none

!    do iID=1,3
    do iID=1,0
       
       m = VM_Mass(KFarr(iID),nReject)
       
       write(*,*) 'ID (old) =',IDarr(iID)
       call timeMeasurement(ForceReset=.true.)
       
       do iMC=1,nMC
          m = VM_Mass(KFarr(iID),nReject)
          call AddHist(hM1,m,1.0)
          call AddHist(hReject1,real(nReject),1.0)
       end do
       call timeMeasurement()
       
       call WriteHist(hM1,121,add=1e-20,mul=1.0/nMC,file="mass1."//Achar(iID+48)//".dat")
       call WriteHist(hReject1,121,add=1e-20,mul=1.0/nMC,file="reject1."//Achar(iID+48)//".dat")
       
       call ClearHist(hM1)
       call ClearHist(hReject1)
       
    end do
    
    do iID=1,5
       
       call AssignMass_1(IDarr(iID), srts,momLRF,mediumAtPos, m, nReject)
       write(*,*) 'ID (new) =',IDarr(iID)
       call timeMeasurement(ForceReset=.true.)
       
       do iMC=1,nMC
          call AssignMass_1(IDarr(iID), srts,momLRF,mediumAtPos, m, nReject)
          call AddHist(hM2,m,1.0)
          call AddHist(hReject2,real(nReject),1.0)
       end do
       call timeMeasurement()
       
       call WriteHist(hM2,121,add=1e-20,mul=1.0/nMC,file="mass2."//Achar(iID+48)//".dat")
       call WriteHist(hReject2,121,add=1e-20,mul=1.0/nMC,file="reject2."//Achar(iID+48)//".dat")
       
       call ClearHist(hM2)
       call ClearHist(hReject2)
       
    end do

  end subroutine checkMeson
  !*************************************************************************

  !*************************************************************************
  subroutine checkBaryon
    use output

    implicit none
    
    iID=2
    m = Delta_Mass(nReject)
       
    write(*,*) 'ID (old) =',iID
    call timeMeasurement(ForceReset=.true.)
       
    do iMC=1,0
       m = Delta_Mass(nReject)
       call AddHist(hM1,m,1.0)
       call AddHist(hReject1,real(nReject),1.0)
    end do
    call timeMeasurement()

    call WriteHist(hM1,121,add=1e-20,mul=1.0/nMC,file="mass1."//intToChar(iID)//".dat")
    call WriteHist(hReject1,121,add=1e-20,mul=1.0/nMC,file="reject1."//intToChar(iID)//".dat")
    
    call ClearHist(hM1)
    call ClearHist(hReject1)


    do iID=1,61
       
       call AssignMass_1(iID, srts,momLRF,mediumAtPos, m, nReject)
       write(*,*) 'ID (new) =',iID
       call timeMeasurement(ForceReset=.true.)
       
       do iMC=1,nMC
          call AssignMass_1(iID, srts,momLRF,mediumAtPos, m, nReject)
          call AddHist(hM2,m,1.0)
          call AddHist(hReject2,real(nReject),1.0)
       end do
       call timeMeasurement()

       call WriteHist(hM2,121,add=1e-20,mul=1.0/nMC,file="mass2."//intToChar(iID)//".dat")
       call WriteHist(hReject2,121,add=1e-20,mul=1.0/nMC,file="reject2."//intToChar(iID)//".dat")
       
       call ClearHist(hM2)
       call ClearHist(hReject2)

    end do

  end subroutine checkBaryon
  !*************************************************************************

  !*************************************************************************
  subroutine checkDelta
    use output

    implicit none

    momLRF = (/0.,0.,0.,0.1/)

    rDens = 0.0 ! DUMMY
    mediumAtPos%density = rhoNull*rDens
    mediumAtPos%densityProton =  0.5 * rhoNull*rDens
    mediumAtPos%densityNeutron = 0.5 * rhoNull*rDens
    mediumAtPos%useMedium=.true.

    mediumAtPos%useMedium=.false.
    
    
    iID=2
    m = Delta_Mass(nReject)
    call AssignMass_1(iID, srts,momLRF,mediumAtPos, m, nReject)

    write(*,*) '(OLD)'
    call timeMeasurement(ForceReset=.true.)
    do iMC=1,nMC
       m = Delta_Mass_Full(nReject)
       m = Delta_Mass(nReject)
       call AddHist(hM1,m,1.0)
       call AddHist(hReject1,real(nReject),1.0)
    end do
    call timeMeasurement()

    call WriteHist(hM1,121,add=1e-20,mul=1.0/nMC,file="mass1."//intToChar(iID)//".dat")
    call WriteHist(hReject1,121,add=1e-20,mul=1.0/nMC,file="reject1."//intToChar(iID)//".dat")
    
    call ClearHist(hM1)
    call ClearHist(hReject1)

    write(*,*) '(NEW)'

    call timeMeasurement(ForceReset=.true.)
       
    do iMC=1,nMC
       call AssignMass_1(iID, srts,momLRF,mediumAtPos, m, nReject)
       call AddHist(hM2,m,1.0)
       call AddHist(hReject2,real(nReject),1.0)
    end do
    call timeMeasurement()
    
    call WriteHist(hM2,121,add=1e-20,mul=1.0/nMC,file="mass2."//intToChar(iID)//".dat")
    call WriteHist(hReject2,121,add=1e-20,mul=1.0/nMC,file="reject2."//intToChar(iID)//".dat")
    
    call ClearHist(hM2)
    call ClearHist(hReject2)

  end subroutine checkDelta
  !*************************************************************************


end program testAssignMassMC
