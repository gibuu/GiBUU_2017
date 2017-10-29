!******************************************************************************
!****m* /InvMassesAnalysis
! NAME
! module InvMassesAnalysis
!
! PURPOSE
! Provides routines to produce histograms of the invariant mass distributions
! of pi+N,pi+pi,N+N pairs in the final particle vector.
!
! Used by HiPionAnalysis and HiLeptonAnalysis
!******************************************************************************
module InvMassesAnalysis
  use histf90
  use histMC
  implicit none

  private

  type(histogram), dimension(5,5), save :: hInvMasses
  type(histogramMC), dimension(1:3,4:5), save :: hInvMassesExcl
  type(histogram), save :: HnPions

  public :: InvMasses_INIT
  public :: InvMasses_Write
  public :: InvMasses_FillEvent

contains


  !****************************************************************************
  !****s* InvMassesAnalysis/InvMasses_INIT
  ! NAME
  ! subroutine InvMasses_INIT
  !
  ! PURPOSE
  ! Initialize the histograms to be used in this module
  !****************************************************************************
  subroutine InvMasses_INIT
    use IdTable, only: nres, Nucleon, F37_1950
    use particleProperties, only: hadron

    integer :: i,j,k
    character(3),dimension(5) :: NN = (/"pi-","pi0","pi+","n  ","p  "/)

    do i=1,5
       do j=i,5
          call CreateHist(hInvMasses(i,j),"InvMasses "//NN(i)//" "//NN(j),0.0,2.5,0.01)
       end do
    end do

    do i=1,3  ! pions
      do j=4,5  ! nucleons
        call CreateHistMC(hInvMassesExcl(i,j),"InvMasses (excl) "//NN(i)//" "//NN(j),0.0,2.5,0.01,2*(nres+1))
        hInvMassesExcl(i,j)%xDesc='pi-N invariant mass [GeV]'
        do k=Nucleon,F37_1950
          hInvMassesExcl(i,j)%yDesc(k)          = trim(hadron(k)%name)
          hInvMassesExcl(i,j)%yDesc(F37_1950+k) = trim(hadron(k)%name) // " (BG)"
        end do
      end do
    end do

    call CreateHist(HnPions,"number of pions",-0.5,21.0,1.0)

  end subroutine InvMasses_INIT


  !****************************************************************************
  !****s* InvMassesAnalysis/InvMasses_Write
  ! NAME
  ! subroutine InvMasses_Write (mul)
  !
  ! PURPOSE
  ! Write out the histograms
  !
  ! INPUTS
  ! * real, intent(in) :: mul -- multiplicative factor for writing histograms
  !
  !****************************************************************************
  subroutine InvMasses_Write (mul)
    use output, only: writeFileDocu, intToChar

    real, intent(in) :: mul

    integer :: i,j
    logical, save :: first = .true.

    do i=1,5
       do j=i,5
          call WriteHist(hInvMasses(i,j),140,mul=mul,file='InvMasses.'//trim(intToChar(10*i+j))//'.dat')
       end do
    end do

    do i=1,3    ! pions
      do j=4,5  ! nucleons
        call WriteHistMC (hInvMassesExcl(i,j),'InvMasses_excl.'//trim(intToChar(10*i+j))//'.dat',mul=mul)
      end do
    end do

    call WriteHist(HnPions,140,mul=mul,file='InvMasses.nPionsEvents.dat')

    if (first) then
       call writeFileDocu('nPionsEvents.dat','histogram with number of pions per event')
       call writeFileDocu('InvMasses.0nn.dat (nn=10x1..5 + 1..5)','invariant Mass of particle pair; 1..5=pi-,pi0,pi+,n,p')
       first = .false.
    end if

  end subroutine InvMasses_Write


  !****************************************************************************
  !****if* InvMassesAnalysis/partClass
  ! NAME
  ! integer function partClass(ID,IZ,anti)
  ! PURPOSE
  ! Put a particle in the classes: "pi-","pi0","pi+","n  ","p  "
  !****************************************************************************
  integer function partClass(ID,IZ,anti)
    integer,intent(in) :: ID,iZ
    logical, intent(in) :: anti

    partClass = 0
    if (anti) return

    select case (ID)
    case (101)
       partclass = 2+IZ ! 1,2,3
    case (1)
       partclass = 4+IZ ! 4,5
    end select

  end function partClass


  !****************************************************************************
  !****s* InvMassesAnalysis/InvMasses_FillEvent
  ! NAME
  ! subroutine InvMasses_FillEvent(L)
  ! PURPOSE
  ! Given one event, this routine calculates all possible invariant
  ! masses and fills the corresponding histograms.
  !
  ! INPUTS
  ! * type(tAnaEvent) :: L -- the event to analyse
  !
  ! NOTES
  ! * At the moment this routine does not respect the acceptance of any
  !   detector. We would have to find a way to give these weights as
  !   an additional argument to the routine. (A direct call to some
  !   routine, as done before directly in the module HiLeptonAnalysis,
  !   is not appropriate.)
  !
  !   It could e.g. be done as for 'FindRho0':
  !     allocate(Probs(EventArr(i,j)%particleList%nEntries))
  !     pNode => EventArr(i,j)%particleList%first
  !     Probs = 0
  !     ii = 1
  !     do while (ASSOCIATED(pNode))
  !       call checkCuts(pNode%V,nu,Q2,Ebeam,phi_Lepton,AccFlag,AccWeight)
  !       Probs(ii) = AccWeight
  !       pNode=>pNode%next
  !       ii=ii+1
  !     end do
  !****************************************************************************
  subroutine InvMasses_FillEvent(L)
    use particleDefinition, only: sqrts
    use particlePointerListDefinition
    use AnaEventDefinition
    use IdTable, only: nres, isBaryonResonance
    use history, only: history_getParents

    type(tAnaEvent), intent(in) :: L

    type(tParticleListNode), POINTER :: pNode1,pNode2
    integer :: iH1, iH2, nPion, nnPion, ch, parents1(3),parents2(3)
    real :: M2, Weight, Weight2
!    real :: accW1,AccW2
!    type(particle) :: Part

    nPion = 0
    nnPion = L%numberparticles(1,-1)+L%numberparticles(1,0)+L%numberparticles(1,1)

    pNode1 => L%particleList%first

    do
       if (.not. associated(pNode1)) exit

       Weight = pNode1%V%perWeight ! should always be the same ;)

       iH1 = PartClass(pNode1%V%ID,pNode1%V%charge,pNode1%V%antiparticle)
       if (iH1>0) then

          if (iH1 .le. 3) nPion=nPion+1

!          AccW1 = 1.0
!          Part = pNode1%V
!          call checkCuts(Part,nu,Q2,Ebeam,phi_Lepton,acceptFlag,AccW1)

          pNode2=>pNode1%next
          do
             if (.not. associated(pNode2)) exit

             iH2 = PartClass(pNode2%V%ID,pNode2%V%charge,pNode2%V%antiparticle)
             if (iH2>0) then

!                AccW2 = 1.0
!                Part = pNode2%V
!                call checkCuts(Part,nu,Q2,Ebeam,phi_Lepton,acceptFlag,AccW2)

                M2 = sqrtS(pNode1%V,pNode2%V)


                !--- Acceptance
!                Weight2 = Weight * AccW1*AccW2

                !--- Exlusivity
                Weight2 = Weight
                if (nnPion.ne.2) Weight2 = 0.0

                call AddHist(hInvMasses(min(iH1,iH2),max(iH1,iH2)),M2, Weight,Weight2)

                if (nnPion==1 .and. min(iH1,iH2)<4) then
                  ! exclusive pi-N spectra (including resonance contributions)
                  parents1 = history_getParents (pNode1%V%history)
                  parents2 = history_getParents (pNode2%V%history)

                  ! check if any of the two comes from a resonance decay
                  if (isBaryonResonance(parents1(1)) .and. parents1(1)==parents2(1)) then
                    ! both come from the same resonance
                    ch = parents1(1)
                    if (ch>nres+1) ch = 1
                  else if (isBaryonResonance(parents1(1))) then
                    ! only one of them comes a resonance
                    ch = parents1(1) + nres + 1
                    if (ch>2*(nres+1)) ch = nres + 2
                  else if (isBaryonResonance(parents2(1))) then
                    ch = parents2(1) + nres + 1
                    if (ch>2*(nres+1)) ch = nres + 2
                  else
                    ch = 1
                  end if

                  call AddHistMC (hInvMassesExcl(min(iH1,iH2),max(iH1,iH2)), M2, ch ,Weight)
                end if

             end if

             pNode2=>pNode2%next
          end do
       end if

       pNode1=>pNode1%next
    end do

    call AddHist(HnPions,real(nPion),Weight)

  end subroutine InvMasses_FillEvent



end module InvMassesAnalysis
