!******************************************************************************
!****m* /AssignMassMC
! NAME
! module AssignMassMC
! PURPOSE
! Assign masses for 1-,2- or 3-body final state according spectral functions
!******************************************************************************
module AssignMassMC
  use mediumDefinition
  implicit none

  private

  public :: AssignMass_1
  public :: AssignMass_2
  public :: AssignMass_1_Therm

contains
  !****************************************************************************
  !****s* AssignMassMC/AssignMass_1
  ! NAME
  ! subroutine AssignMass_1(ID, srts,momLRF,mediumAtPos, mass, nReject)
  !
  ! PURPOSE
  ! return a mass for a particle, which is distributed according the
  ! actual parametrizations of the spectral function, i.e. according
  !   A(m)
  !
  ! For stable particles, it returns the pole mass, while for particles
  ! with non-vanishing width the returned value is distributed according
  ! a relativistiv Breit-Wigner distribution with a mass dependend width.
  !
  ! INPUTS
  ! * integer :: ID -- ID of the particle
  ! * real :: srts -- maximal energy available
  ! * real, dimension(0:3) :: momLRF -- momentum of resonance in LRF
  ! * type(medium) :: mediumAtPos -- medium at position
  !
  ! OUTPUT
  ! * real :: mass -- the selected mass value
  ! * integer, OPTIONAL:: nReject -- number of steps in rejection method
  !
  ! NOTES
  ! * mass shifts due to potentials (vector mesons) are not considered yet
  !****************************************************************************
  subroutine AssignMass_1(ID, srts,momLRF,mediumAtPos, mass, nReject)
    use MassAssInfoDefinition
    use IdTable, only: isMeson, isBaryon
    use CALLSTACK, only: TRACEBACK
    use mesonWidthMedium, only: GetMassAssInfo_Meson, WidthMesonMedium
    use baryonWidthMedium, only: GetMassAssInfo_Baryon,WidthBaryonMedium
    use random
    use monteCarlo, only: MonteCarloChoose
    use output

    integer, intent(in) :: ID
    real, intent(in) :: srts
    real,intent(in),dimension(0:3)  :: momLRF
    type(medium),intent(in) :: mediumAtPos
    real, intent(out) :: mass
    integer,intent(out),OPTIONAL :: nReject

    type(tMassAssInfo),save :: MAI
    integer :: iTry
    integer,parameter :: nTry = 100

    integer :: iB,nB
    real :: BinWeightTot
    real :: y, ymin, ymax, gamma, Q, Qmax

    character(*), parameter :: form1 =  &
         & '("Problem in AssignMass_1: (nTry)    ID=",i3," mass= (",2G12.4,")")'
    character(*), parameter :: form2 =  &
         & '("Message in AssignMass_1: (Q-viol)  ID=",i3," mass= (",2G12.4,"): ",G12.5," Q=",G12.5)'

    mass = 0.0
    if (present(nReject)) nReject=0

!    call writeMedium(mediumAtPos)



    if (isMeson(ID)) then
       call GetMassAssInfo_Meson(MAI,ID,momLRF,mediumAtPos)
       gamma = WidthMesonMedium(ID,srts,momLRF,mediumAtPos)
    else if (isBaryon(ID)) then
       call GetMassAssInfo_Baryon(MAI,ID,momLRF,mediumAtPos)
       gamma = WidthBaryonMedium(ID,srts,momLRF,mediumAtPos)
    else
       call TRACEBACK("ID not allowed")
    end if

    call SetUpperBin(MAI,srts,gamma,nB)

    if (MAI%IsStable) then
       mass = MAI%Mass0
       return ! ==> solution found
    end if

    MAI%W = MAI%W*MAI%Q


    iTry = 0
    do
       iTry = iTry+1

       ! STEP 1: Select the mass region:
       ! (it is important to have this step in the overall loop, since
       ! this cures inaccuracies in the calculation of the weights)
       iB = MonteCarloChoose(MAI%W,BinWeightTot)

       ymin = MAI%Y(iB)
       ymax = MAI%Y(iB+1)
       Qmax = MAI%Q(iB)

       ! STEP 2: generate random value according Cauchy with constant width

       y = ymin + rn()*(ymax-ymin)
       mass = 0.5*tan(0.5*y)*MAI%Gamma0 + MAI%Mass0

       ! STEP 3: Do the rejection

       if (isMeson(ID)) then
          gamma = WidthMesonMedium(ID,mass,momLRF,mediumAtPos)
       else
          gamma = WidthBaryonMedium(ID,mass,momLRF,mediumAtPos)
       end if

       Q = MassAssInfoQ(MAI%Mass0,MAI%Gamma0,mass,gamma)

       if (Q>1.05*Qmax) then
!       if (Q>1.01*Qmax) then
!       if (Q>1.001*Qmax) then
          write(*,form2) ID,MAI%M(iB:IB+1),mass, Q/Qmax
       end if

!!$       if (Q>1.15*Qmax) then
!!$          write(*,*) '----1.10----'
!!$          write(*,*) momLRF
!!$          write(*,*) iB
!!$          write(*,*) MAI%Mass0,MAI%Gamma0
!!$          write(*,*) mass,gamma
!!$          write(*,*) Q,Qmax
!!$          write(*,'(A10,20f13.4)') 'BinM',MAI%M
!!$          write(*,'(A10,7(" "),20f13.4)') 'BinW',MAI%W(1:size(MAI%Q)-1)/sum(MAI%W(1:size(MAI%Q)-1))
!!$          write(*,'(A10,7(" "),20f13.4)') 'BinQ',MAI%Q(1:size(MAI%Q)-1)
!!$          stop
!!$       end if

       if (present(nReject)) nReject=iTry

       if (rn()*Qmax <= Q) return ! ==> solution found

       if (iTry >= nTry) then
          if (DoPr(1)) write(*,form1) ID, MAI%M(iB),MAI%M(iB+1)
          mass = MAI%Mass0
          return ! ==> solution found
       end if

    end do


  end subroutine AssignMass_1

  !****************************************************************************
  !****s* AssignMassMC/AssignMass_2
  ! NAME
  ! subroutine AssignMass_2(ID, srts,momLRF,mediumAtPos, mass, nReject)
  !
  ! PURPOSE
  ! return masses for the two outgoing particles, which are selected
  ! according
  !   A(m1)*A(m2)*p_{cm}
  !
  ! INPUTS
  ! * integer, dimension(2) :: ID -- ID of the particles
  ! * real :: srts -- maximal energy available
  ! * real, dimension(0:3) :: momLRF -- momentum of resonance in LRF
  ! * type(medium) :: mediumAtPos -- medium at position
  !
  ! OUTPUT
  ! * real, dimension(2) :: mass -- the selected mass values
  ! * integer, OPTIONAL :: nReject -- number of steps in rejection method
  !
  ! NOTES
  ! * momLRF is not respected correctly
  ! * angular dependencies are not considered
  ! * mass shifts due to potentials (vector mesons) are not considered yet
  !****************************************************************************
  subroutine AssignMass_2(ID, srts,momLRF,mediumAtPos, mass, nReject)
    use MassAssInfoDefinition
    use IdTable, only: isMeson, isBaryon
    use CALLSTACK, only: TRACEBACK
    use mesonWidthMedium, only: GetMassAssInfo_Meson, WidthMesonMedium
    use baryonWidthMedium, only: GetMassAssInfo_Baryon,WidthBaryonMedium
    use random
    use monteCarlo, only: MonteCarloChoose2Dim
    use output
    use twoBodyTools, only: pCM_sqr

    integer, intent(in),dimension(2) :: ID
    real, intent(in) :: srts
    real,intent(in),dimension(0:3)  :: momLRF
    type(medium),intent(in) :: mediumAtPos
    real, intent(out),dimension(2) :: mass
    integer,intent(out),OPTIONAL :: nReject

    integer :: iTry
    integer,parameter :: nTry = 10000
    integer :: i,i1,i2
    type(tMassAssInfo),dimension(2), save :: MAI
    integer, dimension(2) :: nBin,nnBin,iB
    real,dimension(:,:),allocatable :: Warr,Qarr
    real :: p_cd
    real :: BinWeightTot
    real :: y, ymin, ymax, gamma, Qmax, Q(2),maxmass(2)



    mass = 0.0
    if (present(nReject)) nReject=0

    do i=1,2

       if (isMeson(ID(i))) then
          call GetMassAssInfo_Meson(MAI(i),ID(i),momLRF,mediumAtPos)
       else if (isBaryon(ID(i))) then
          call GetMassAssInfo_Baryon(MAI(i),ID(i),momLRF,mediumAtPos)
       else
          write(*,*) i,ID(i)
          call TRACEBACK("ID not allowed")
       end if
       nBin(i) = size(MAI(i)%W)
       maxmass(3-i) = srts - MAI(i)%M(1)

    end do
    ! we need a second loop !!!
    do i=1,2

       if (isMeson(ID(i))) then
          gamma = WidthMesonMedium(ID(i),maxmass(i),momLRF,mediumAtPos)
       else
          gamma = WidthBaryonMedium(ID(i),maxmass(i),momLRF,mediumAtPos)
       end if

       call SetUpperBin(MAI(i),maxmass(i),gamma,nBin(i))
    end do

    allocate(Warr(nBin(1),nBin(2)))
    allocate(Qarr(nBin(1),nBin(2)))
    Warr = 0.0
    Qarr = 0.0


    MAI(1)%W = MAI(1)%W*MAI(1)%Q
    MAI(2)%W = MAI(2)%W*MAI(2)%Q

    nnBin = 0

    do i1=1,nBin(1)
       do i2=1,nBin(2)

          p_cd = pCM_sqr(srts**2,MAI(1)%M(i1)**2,MAI(2)%M(i2)**2)
          if (MAI(1)%M(i1)+MAI(2)%M(i2) .gt. srts) p_cd = 0.0
          if (p_cd.gt.0.0) then
             p_cd = sqrt(p_cd)
             Warr(i1,i2)=MAI(1)%W(i1)*MAI(2)%W(i2) * p_cd
             Qarr(i1,i2)=MAI(1)%Q(i1)*MAI(2)%Q(i2) * p_cd
             if (i1.gt.nnBin(1)) nnBin(1) = i1
             if (i2.gt.nnBin(2)) nnBin(2) = i2
          end if

       end do
    end do

    if (nnBin(1).eq.0) then
       write(*,*) ID, srts
       call TRACEBACK('all bins closed')
    end if

!!$    write(*,*) nBin
!!$    write(*,*) nnBin
!!$    do i1=1,nnBin(1)
!!$       write(*,'(1P,20e13.5)') Warr(i1,1:nnBin(2))
!!$    end do
!!$    stop

    iTry = 0
    do
       iTry = iTry+1

!!$       write(*,*)
!!$       write(*,*) '==='

       if (iTry >= nTry) then
           call TRACEBACK('iTry too large')
       end if

       ! STEP 1: Select the mass region:
       ! (it is important to have this step in the overall loop, since
       ! this cures inaccuracies in the calculation of the weights)
       iB = MonteCarloChoose2Dim(Warr(1:nnBin(1),1:nnBin(2)),BinWeightTot)

!!$       write(*,*) iB,BinWeightTot

       Qmax = Qarr(iB(1),iB(2))

       ! STEP 2: generate random value according Cauchy with constant width

       do i=1,2
          if (MAI(i)%IsStable) then
             mass(i) = MAI(i)%Mass0
             Q(i) = 1.0
          else
             ymin = MAI(i)%Y(iB(i))
             ymax = MAI(i)%Y(iB(i)+1)
             y = ymin + rn()*(ymax-ymin)
             mass(i) = 0.5*tan(0.5*y)*MAI(i)%Gamma0 + MAI(i)%Mass0
             if (isMeson(ID(i))) then
                gamma = WidthMesonMedium(ID(i),mass(i),momLRF,mediumAtPos)
             else
                gamma = WidthBaryonMedium(ID(i),mass(i),momLRF,mediumAtPos)
             end if
             Q(i) = MassAssInfoQ(MAI(i)%Mass0,MAI(i)%Gamma0,mass(i),gamma)
          end if

       end do

!!$       write(*,*) '1:',mass(1),Q(1)
!!$       write(*,*) '2:',mass(2),Q(2)

       p_cd = pCM_sqr(srts**2,mass(1)**2,mass(2)**2)
       if (mass(1)+mass(2) .gt. srts) p_cd = 0.0
       if (p_cd.le.0.0) cycle
       p_cd = sqrt(p_cd)
!!$       write(*,*) 'p_cd=',p_cd

       if (present(nReject)) nReject=iTry

       if (rn()*Qmax <= Q(1)*Q(2)*p_cd) exit ! ==> solution found

    end do


    deallocate(Warr)
    deallocate(Qarr)

  end subroutine AssignMass_2

  !****************************************************************************
  !****s* AssignMassMC/AssignMass_1_Therm
  ! NAME
  ! subroutine AssignMass_1_Therm(ID, T, srts,momLRF,mediumAtPos, mass, nReject)
  !
  ! PURPOSE
  ! return a mass for a particle, which is distributed according the
  ! actual parametrizations of the spectral function multiplied with a
  ! Boltzmann-Factor for the temperature T, i.e. according
  !   A(m) * m^2 * K_2(m/T)
  !
  ! For stable particles, it returns the pole mass, while for particles
  ! with non-vanishing width the returned value is distributed according
  ! a relativistiv Breit-Wigner distribution with a mass dependend width.
  !
  ! INPUTS
  ! * integer :: ID -- ID of the particle
  ! * real :: T -- Teemperature
  ! * real :: srts -- maximal energy available
  ! * real, dimension(0:3) :: momLRF -- momentum of resonance in LRF
  ! * type(medium) :: mediumAtPos -- medium at position
  !
  ! OUTPUT
  ! * real :: mass -- the selected mass value
  ! * integer, OPTIONAL:: nReject -- number of steps in rejection method
  !
  ! NOTES
  ! * mass shifts due to potentials (vector mesons) are not considered yet
  !****************************************************************************
  subroutine AssignMass_1_Therm(ID, T,srts,momLRF,mediumAtPos, mass, nReject)

    use particleProperties, only: hadron
    use besselK, only: BesselK2
    use random

    integer, intent(in) :: ID
    real, intent(in) :: T
    real, intent(in) :: srts
    real,intent(in),dimension(0:3)  :: momLRF
    type(medium),intent(in) :: mediumAtPos
    real, intent(out) :: mass
    integer,intent(out),OPTIONAL :: nReject

    integer :: nRejectTMP,nReject2,iTry
    real :: fak,fak0,minMass

    minMass = hadron(ID)%minMass
    fak0 = (minmass**2*BesselK2(minMass/T))

    iTry=0
    do
       iTry = iTry+1

       call AssignMass_1(ID, srts,momLRF,mediumAtPos, mass, nRejectTMP)

       nReject2 = nReject2+nRejectTMP
       !       if (present(nReject)) nReject=iTry
       if (present(nReject)) nReject=nReject2

       fak = mass**2*BesselK2(mass/T)

       if (rn()*fak0 <= fak) return ! ==> solution found
    end do


  end subroutine AssignMass_1_Therm

end module AssignMassMC
