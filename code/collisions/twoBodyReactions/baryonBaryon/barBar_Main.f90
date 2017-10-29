!******************************************************************************
!****m* /barBar_Main
! NAME
! module barBar_Main
! PURPOSE
! This module administrates all baryon-baryon cross sections.
!******************************************************************************
module barBar_Main

  implicit none
  private

  logical, parameter :: debug = .false.

  public :: XsectionBarBar, eta_deuteron

contains


  !****************************************************************************
  !****s* barBar_Main/XsectionBarBar
  ! NAME
  ! subroutine XsectionBarBar (srts, teilchenIN, mediumATcollision, teilchenOUT, sigmaTot, sigmaElast, pauliIncluded, plotFlag)
  ! PURPOSE
  ! This routine is the main routine for baryon-baryon scattering and its cross sections.
  ! Determines total and elastic cross section and makes a Monte-Carlo decision for a special
  ! reaction channnel. This leads to a definition of ID and charge of teilchenOut, which is the final state vector.
  ! INPUTS
  ! * real, intent(in)                            :: srts               -- sqrt(s) in the process
  ! * type(particle), dimension(1:2), intent(in)  :: teilchenIn         -- colliding particles
  ! * type(medium), intent(in)                    :: mediumATcollision  -- Medium information at the position of the collision
  ! * character(len=*), intent(in), optional      :: plotFlag           -- Switch on plotting of the Xsections
  ! OUTPUT
  ! * type(particle), dimension(1:3), intent(out) :: teilchenOut        -- outgoing particles
  ! * real, intent(out)                           :: sigmaTot           -- total Xsection
  ! * real, intent(out)                           :: sigmaElast         -- elastic Xsection
  ! * logical, intent(out)                        :: pauliIncluded      -- true = cross section includes Pauli blocking
  ! NOTES
  ! plotFlag=.true. causes to make output to the files:
  ! * 'BaryonBaryon_Crosssection_1.dat'
  ! * 'BaryonBaryon_Crosssection_2.dat'
  ! The content is explained in the files.
  !****************************************************************************
  subroutine XsectionBarBar (srts, partIn, mediumAtColl, partOut, sigmaTot, sigmaElast, pauliIncluded, plotFlag)
      use mediumDefinition
      use particleDefinition
      use particleProperties, only: hadron
      use IDTable
      use random, only: rn
      use barBar_barBar, only: sigmaBB, chooseCharge
      use preEventDefinition, only: preEvent
      use barBar_barbarMes, only:  NN_NNpi_direct, sig_NNV, chooseCharge_NNpiomega
      use barBar_barHypKaon, only: barBar_barBarMeson_strange, get_Channels_BYK
      use hypNuc_hypNuc, only: hypNuc_hypNuc_Main, get_Channels_YN

      real, intent(in)                           :: srts               ! sqrt(s) in the process
      type(particle),dimension(1:2), intent(in)  :: partIn         ! colliding particles
      type(medium), intent(in)                   :: mediumAtColl  ! Medium informations at the position of the collision
      character(len=*), intent(in), optional     :: plotFlag           ! Switch on plotting of the Xsections

      type(preEvent),dimension(1:4), intent(out) :: partOut        ! outgoing particles
      real, intent(out)                          :: sigmaTot           ! total Xsection
      real, intent(out)                          :: sigmaElast         ! elastic Xsection
      logical, intent(out)                       :: pauliIncluded      ! true = cross section includes Pauli blocking

      ! Xsections of all possible channels:
      real, dimension(nucleon:F37_1950) :: sigmaNukBar      ! Cross section for producing nucleon+baryon.
                                                            ! Here the index is the ID of the produced baryon. The cross section
                                                            ! is summed over all possible charge states of the final state.
      real, dimension(Delta:F37_1950)   :: sigmaDeltaBar    ! Cross section for producing Delta+baryon,
                                                            ! summed over all possible charge states of the deltas.
      real, dimension(1:18)             :: Sigma_BYK        ! Cross section for Bar Bar -> Bar Hyperon Kaon
      real, dimension(1:4)              :: Sigma_NNPion     ! Cross section for "N N Pion" production as a background term
      real, dimension(1:3)              :: Sigma_NNV        ! Cross section for vector meson production ("N N omega/phi") as a background term
      real                              :: sigma_np_deta    ! cross section for n+p -> d+eta
      real, dimension(1:4)              :: Sigma_YN         ! Cross section for Hyperon N -> Hyperon N

      ! Local variables:
      integer :: baryonID, channelID, strange
      real :: x, wahrscheinlichkeit

      pauliIncluded=.false.

      ! Initialize output: produced particles
      partOut(:)%ID = 0
      partOut(:)%charge = 0
      partOut(:)%antiParticle = .false.
      partOut(:)%mass = 0.

      !************************************************************************
      ! Initializing the cross sections:
      !************************************************************************

      sigmaTot=0.
      sigmaElast=0.
      sigmaNukBar=0.
      sigmaDeltaBar=0.
      Sigma_BYK=0.
      Sigma_NNPion=0.
      Sigma_NNV=0.
      sigma_np_deta=0.
      Sigma_YN=0.

      if ( (hadron(partIn(1)%ID)%charm.ne.0).or.(hadron(partIn(2)%ID)%charm.ne.0) ) then
         ! At least one charmed baryon as scattering partner!
         ! Here are no cross sections implemented for this incoming particles
         return

      else if ( (hadron(partIn(1)%ID)%strangeness.ne.0).and.(hadron(partIn(2)%ID)%strangeness.ne.0) ) then
         ! Two strange baryons as scattering partner!
         ! Here are no cross sections implemented for this incoming particles
         return

      else if ( (hadron(partIn(1)%ID)%strangeness.ne.0).or.(hadron(partIn(2)%ID)%strangeness.ne.0) ) then
         ! Hyperon + Nucleon scattering
         strange = 1

         ! Lambda/Sigma/Xi + Nucleon scattering
         Sigma_YN = hypNuc_hypNuc_Main (srts, partIn)

         sigmaTot = sum(Sigma_YN)
         sigmaElast = Sigma_YN(1)

      else ! All baryons are non-strange and non-charmed

         strange = 0

         if (partIn(1)%Id==nucleon .or. partIn(2)%Id==nucleon) then
            ! N N, N Delta or N Res as initial particles
            do baryonId=nucleon,F37_1950
               ! N+N -> N+baryon
               if (hadron(baryonId)%propagated) then
                 sigmaNukBar(baryonID) = sigmaBB(partIn,(/nucleon,baryonID/),mediumAtColl,srts,pauliIncluded)
                 if (debug) &
                   & write(*,*) 'In XsectionBarBar: baryonId,sigmaNukBar:',&
                   & baryonId,sigmaNukBar(baryonID)
               else
                 sigmaNukBar(baryonID)=0.
               end if
            end do

            if (partIn(1)%Id==nucleon .and. partIn(2)%Id==nucleon) then
               do baryonId=Delta,F37_1950
                 ! N+N -> Delta+baryon
                 if (hadron(baryonId)%propagated) then
                   sigmaDeltaBar(baryonID) = sigmaBB (partIn,(/delta,baryonID/),mediumAtColl,srts,pauliIncluded)
                 else
                   sigmaDeltaBar(baryonID) = 0.
                 end if
               end do
               ! The background Xsections for NN-> NNPion
               Sigma_NNPion = NN_NNpi_direct (srts)
               ! Field Sigma_NNPION
               ! 1: p p -> p p pi^0 ; n n -> n n pi^0
               ! 2: p n -> n n pi^+ ; p n -> p p pi^-
               ! 3: p n -> p n pi^0 ;
               ! 4: p p -> p n pi^+ ; n n -> p n pi^-
               if (sum(partIn%charge) == 1) then
                  !np->X
                  Sigma_NNPION(1)=0
                  Sigma_NNPION(4)=0
                  Sigma_NNPION(2)=2.*Sigma_NNPION(2)
                  ! Factor 2 since we have two final state channels with the
                  ! same Xsection: p n -> n n pi^+ ; p n -> p p pi^-

                  sigma_np_deta = eta_deuteron(srts)
               else
                  !pp->X, nn->X
                  Sigma_NNPION(2:3)=0.
               end if
              ! Evaluate background cross sections for NN -> NN omega/phi (currently no res. with omega decay channel are produced in NN coll.)
              sigma_NNV(1:3) = sig_NNV (srts, sum(partIn%charge))
            end if
         else if (partIn(1)%Id==Delta .or. partIn(2)%Id==Delta) then
            ! Delta + Res in the initial channel
            sigmaNukBar=0.
            sigmaDeltaBar=0.
            ! Delta + Res -> N + N
            sigmaNukBar(nucleon) = sigmaBB (partIn,(/nucleon,nucleon/),mediumAtColl,srts,pauliIncluded)
         else
            ! There are no cross sections implemented for these incoming particles!
            sigmaElast=0.
            sigmaTot=0.
            return
         end if

         Sigma_BYK = barBar_barBarMeson_strange (srts, partIn)  ! BB -> BYK

         !*********************************************************************
         ! Set Sigmatot by summing over all channels
         !*********************************************************************

         sigmaTot = Sum(sigmaNukBar) + sum(sigmaDeltaBar) + SUM(Sigma_NNPion) + sigma_np_deta + sum(Sigma_NNV) + Sum(Sigma_BYK)

         if (debug) write(*,*) 'In XsectionBarBar: sigmatot:', sigmaTot

         if (present(PlotFlag)) call PlotSigmas

         !*********************************************************************
         ! Set sigmaElast. Only elastic BN -> BN cross sections are included
         !*********************************************************************

         if ((partIn(1)%Id==nucleon) .or. (partIn(2)%Id==nucleon)) then
            sigmaElast=sigmaBB(partIn,partIn(1:2)%ID,mediumAtColl,srts,pauliIncluded)
         else
            sigmaElast=0.
         end if

         if (debug) write(*,*) 'In XsectionBarBar: sigmaElast:', sigmaElast

      end if
      !************************************************************************
      ! Choose channel by Monte-Carlo
      !************************************************************************

      ! Check that there are open channels
      if (sigmaTot < 0.0001) then
         sigmaTot=0.
         return
      end if

      do
         x = rn()*sigmaTot
         if (x>1E-20) exit
      end do
      wahrscheinlichkeit = 0.

      select case (strange)

      case (0) !strange=0: baryon+baryon -> baryon+baryon

      ! nucleon baryon production
      do baryonId = lbound(sigmaNukBar,dim=1), ubound(sigmaNukBar,dim=1)
         wahrscheinlichkeit=wahrscheinlichkeit+sigmaNukBar(baryonID)
         if (wahrscheinlichkeit >= x) then
            partOut(1:2)%ID = (/nucleon,baryonID/)
            partOut(1:2)%charge = chooseCharge (partIn, partOut(1:2)%ID)
            return
         end if
      end do
      ! delta baryon production
      do baryonId = lbound(sigmaDeltaBar,dim=1), ubound(sigmaDeltaBar,dim=1)
         wahrscheinlichkeit = wahrscheinlichkeit + sigmaDeltaBar (baryonID)
         if (wahrscheinlichkeit >= x) then
            partOut(1:2)%ID = (/delta,baryonID/)
            partOut(1:2)%charge = chooseCharge (partIn, partOut(1:2)%ID)
            return
         end if
      end do
      ! Nucleon nucleon meson production
      ! NN -> NN Pion
      do channelId=1,4
         wahrscheinlichkeit=wahrscheinlichkeit+Sigma_NNPion(channelID)
         if (wahrscheinlichkeit.ge.x) then

           if (channelID.eq.1) then
               ! 1: p p -> p p pi^0 ; n n -> n n pi^0
               partOut(1:3)%ID=(/nucleon,nucleon,pion/)
               partOut(1:2)%charge=partIn%charge
               partOut(3)%charge=0
            else if (channelID.eq.2) then
               ! 2: p n -> n n pi^+ ; p n -> p p pi^-
               partOut(1:3)%ID=(/nucleon,nucleon,pion/)
               ! Monte Carlo for the two channels
               if (rn().gt.0.5) then
                  ! p n -> p p pi^-
                  partOut(1:2)%charge=(/1,1/)
                  partOut(3)%charge=-1
               else
                  !  p n -> n n pi^+
                  partOut(1:2)%charge=(/0,0/)
                  partOut(3)%charge=1
               end if
            else if (channelID.eq.3) then
               ! 3: p n -> p n pi^0 ;
               partOut(1:3)%ID=(/nucleon,nucleon,pion/)
               partOut(1:2)%charge=(/1,0/)
               partOut(3)%charge=0
            else
               ! 4: p p -> p n pi^+ ; n n -> p n pi^-
               partOut(1:3)%ID=(/nucleon,nucleon,pion/)
               partOut(1:2)%charge=(/1,0/)
               partOut(3)%charge=Sum(partIn%charge)-Sum(partOut(1:2)%charge)
               if (abs(partOut(3)%charge).ne.1) then
                  write(*,*) 'error in barbar_Main', partIn%ID, partIn%charge,'####', partOut%ID, &
                       & partOut%charge,'####', channelID,Sigma_NNPION
               end if
            end if

            return
         end if
      end do


      ! NN -> NN omega
      wahrscheinlichkeit = wahrscheinlichkeit + Sigma_NNV(1)
      if (wahrscheinlichkeit>x) then
        partOut(1:3)%ID = (/nucleon,nucleon,omegaMeson/)
        partOut(1:2)%charge = partIn%charge
        partOut(3)%charge = 0
        return
      end if
      ! NN -> NN phi
      wahrscheinlichkeit = wahrscheinlichkeit + Sigma_NNV(2)
      if (wahrscheinlichkeit>x) then
        partOut(1:3)%ID = (/nucleon,nucleon,phi/)
        partOut(1:2)%charge = partIn%charge
        partOut(3)%charge = 0
        return
      end if
      ! NN -> NN pi omega
      wahrscheinlichkeit = wahrscheinlichkeit + Sigma_NNV(3)
      if (wahrscheinlichkeit>x) then
        partOut(1:4)%ID = (/nucleon,nucleon,pion,omegaMeson/)
        partOut(1:3)%charge = chooseCharge_NNpiomega (partIn%charge)
        partOut(4)%charge = 0
        return
      end if
      ! BB -> B Hyperon Kaon with B=N,Delta; Hyperon=Lambda,Sigma; Kaon=K^{+,0}
      do channelID=1,Size(Sigma_BYK)
        wahrscheinlichkeit=wahrscheinlichkeit+Sigma_BYK(channelID)
        if (wahrscheinlichkeit>x) then
            partOut(1:3) = get_Channels_BYK(channelID)
            return
        end if
      end do
      ! n + p -> d + eta
      wahrscheinlichkeit = wahrscheinlichkeit + sigma_np_deta
      if (wahrscheinlichkeit>x) then
         partOut(1:3)%ID = (/ nucleon, nucleon, eta /)    ! assume "n p eta" instead of "d eta" for now
         partOut(1:3)%charge = (/ 1, 0 , 0 /)
         return
      end if

      case (1) !strange=1: hyperon+baryon -> hyperon+baryon

         do channelID=1,Size(Sigma_YN)
            wahrscheinlichkeit=wahrscheinlichkeit+Sigma_YN(channelID)
            if (wahrscheinlichkeit>x) then
               partOut(1:2) = get_Channels_YN (channelID)
               return
            end if
         end do

      end select

      write(*,*) 'Problem in XsectionBarBar. We did not find a channel in the Monte-Carlo-Decision.x= ', x, &
                 'wahrscheinlichkeit=', wahrscheinlichkeit
      stop

  contains

    subroutine PlotSigmas
      use twoBodyTools, only: p_lab
      character(len=100) :: plotFlag_old = ""
      real :: plab

      if (trim(plotFlag) /= trim(plotFlag_old)) then

        open(222,File=trim(plotFlag)//'_NR_DR.dat')
        write(222,'(A,2I3)') '# IDs:', partIn%ID
        write(222,'(A,2I3)') '# charges:', partIn%Charge
        write(222,'(A)')     '# Here we plot:'
        write(222,'(A)')     '# 1) srts, 2) p_lab, 3) tot, 4) sum(NN*), 6) sum(ND*), 6-36) NR, 37) sum(DN*), 38) sum(DD*) 39-68) DR'

        open(223,File=trim(plotFlag)//'_tot_BG.dat')
        write(223,'(A,2I3)') '# IDs:', partIn%ID
        write(223,'(A,2I3)') '# charges:', partIn%Charge
        write(223,'(A)')     '# Here we plot:'
        write(223,'(A)')     '# 1) srts, 2) p_lab, 3-26) cross sections (tot,NNpi,NNomega,BYK)'

        plotFlag_old = trim(plotFlag)

      else

        open(222,File=trim(plotFlag)//'_NR_DR.dat',position='append')
        open(223,File=trim(plotFlag)//'_tot_BG.dat',position='append')

      end if

      plab = p_lab(srts, partIn(1)%mass, partIn(2)%mass)
      write(222,'(68ES13.4E3)') srts, plab, sigmaTot, &
                                sum(sigmaNukBar(3:18)), sum(sigmaNukBar(19:31)), sigmaNukBar(nucleon:F37_1950), &
                                sum(sigmaDeltaBar(3:18)), sum(sigmaDeltaBar(19:31)), sigmaDeltaBar(Delta:F37_1950)
      write(223,'(26ES13.4E3)') srts, plab, sigmaTot, Sigma_NNPion, Sigma_NNV, sum(Sigma_BYK), Sigma_BYK

      close(222)
      close(223)

    end subroutine PlotSigmas

  end subroutine XsectionBarBar


  !****************************************************************************
  !****f* barBar_Main/eta_deuteron
  ! NAME
  ! real function eta_deuteron (srts)
  ! PURPOSE
  ! This routine returns the cross section for "n p -> d eta" in mb, as a function of sqrt(s).
  ! We use a tabulated spline fit of the data points, which is being read from an input file ("eta_deuteron_spline.txt").
  ! The experimental data is taken from: H. Calen et al., PRL 79 (1997) 2642.
  ! INPUTS
  ! * real, intent(in) :: srts   -- sqrt(s) in the process
  ! OUTPUT
  ! * returns the cross section in mb
  !****************************************************************************
  real function eta_deuteron (srts)
    use constants, only: mN
    use particleProperties, only: hadron
    use IdTable, only: eta
    use spline, only: Bsplint2

    real, intent(in) :: srts
    logical, save :: first = .true.
    real, save :: sigma_field(100,1:2)
    real :: Q

    if (first) call readData

    Q = (srts - 2*mN - hadron(eta)%mass)*1.0E3  ! excess energy in MeV

    if (Q < 0.) then
      eta_deuteron = 0.
    else if (Q > sigma_field(100,1)) then
      eta_deuteron = sigma_field(100,2)
    else
      eta_deuteron = Bsplint2(sigma_field(1:100,1),sigma_field(1:100,2),Q)
    end if

    eta_deuteron = eta_deuteron * 1.0E-3  ! convert to mb

  contains

    subroutine readData
      use inputGeneral, only: path_to_input

      integer :: ios, i
      character(200) :: filename

      filename = trim(path_to_input) // "/eta_deuteron_spline.txt"

      open(100,file=trim(filename),status='old',ioStat=ios)

      if (ios/=0) then
        write(*,*) 'ERROR in eta_deuteron'
        write(*,*) 'File', filename, " is not available"
        stop
      end if

      ! read data from file
      ! first column: Q in MeV
      ! second columb: cross section in microbarn
      do i=1,100
        read(100,*) sigma_field(i,1), sigma_field(i,2)
!         print *,i,sigma_field(i,1:2)
      end do

      close(100)
      first = .false.
    end subroutine

  end function eta_deuteron


end module barBar_Main
