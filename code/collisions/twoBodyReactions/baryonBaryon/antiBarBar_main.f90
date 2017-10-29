!******************************************************************************
!****m* /antiBarBar_main
! NAME
! module antiBarBar_main
! PURPOSE
! Contains the subroutine XsectionAntiBarBar, which computes various cross
! sections for antibaryon-baryon collisions in the low-energy mode,
! by using the module 'barAntiBar'.
!******************************************************************************
module antiBarBar_main

  implicit none

  private
  logical, parameter :: debug=.false.

  public :: XsectionAntiBarBar

contains


  !****************************************************************************
  !****s* antiBarBar_main/XsectionAntiBarBar
  ! NAME
  ! subroutine XsectionAntiBarBar(srts,teilchenIN,mediumATcollision,
  ! teilchenOUT,sigmaTot,sigmaElast,sigmaCEX,sigmaAnni,pauliIncluded,plotFlag)
  !
  ! PURPOSE
  ! The main routine for antibaryon-baryon elastic and inelastic scattering.
  ! Determines total, elastic, and annihilation cross section
  ! and makes a Monte-Carlo-Decision for a special reaction channel.
  ! For non-annihilation final states, ID's and charges of teilchenOut are
  ! determined, which is the finalState.
  ! In the case if the annihilation channel is selected, Id's are set to pion
  ! just to signal that annihilation is chosen (the real setting has to be
  ! done by another routine in this case).
  !
  ! INPUTS
  ! * real                          :: srts               -- sqrt(s) in the process
  ! * type(particle),dimension(1:2) :: teilchenIn         -- colliding particles
  ! * type(medium)                  :: mediumATcollision  -- medium information at collision point
  ! * logical, optional             :: plotFlag           -- Switch on plotting of the  Xsections
  !
  ! OUTPUT
  ! * type(preEvent),dimension(1:3) :: teilchenOut -- outgoing particles
  ! * real                          :: sigmaTot    -- total Xsection
  ! * real                          :: sigmaElast  -- elastic Xsection
  ! * real                          :: sigmaCEX    -- charge exchange Xsection
  ! * real                          :: sigmaAnni   -- annihilation Xsection into mesons
  ! * logical                       :: pauliInclud -- true =
  !   cross section includes Pauli blocking
  !
  ! NOTES
  ! plotFlag=.true. causes to do output to the files:
  ! * 'AntiBarBar_Tot_Elast.dat',
  ! * 'AntiBarBar_Part.dat'.
  ! The content is explained in the files.
  !****************************************************************************
  subroutine XsectionAntiBarBar(srts,partIn,mediumAtColl,partOut,&
       sigmaTot,sigmaElast,sigmaCEX,sigmaAnni,pauliIncluded,plotFlag)

    use mediumDefinition
    use particleDefinition
    use IDTable
    use random
    use preEventDefinition, only: preEvent
    use barAntiBar, only: sigmaBarAntiBar
    use NbarN_to_NbarDelta, only: NbarN_to_NbarDelta_Integrated, NbarDelta_to_NbarN
    use twoBodytools, only: get_PInitial
    use rmf, only: getRMF_flag
    use twoBodyTools, only: p_lab

    real, intent(in)                             :: srts
    type(particle),dimension(1:2), intent(in)    :: partIn
    type(medium),                  intent(in)    :: mediumAtColl
    logical, intent(in),optional                 :: plotFlag

    type(preEvent),dimension(1:3), intent(out)   :: partOut
    real, intent(out)                            :: sigmaTot
    real, intent(out)                            :: sigmaElast
    real, intent(out)                            :: sigmaCEX
    real, intent(out)                            :: sigmaAnni
    logical, intent(out)                         :: pauliIncluded

    !**************************************************************************
    ! Xsections of all possible outgoing channels summed over charge states
    ! of outgoing particles:
    !**************************************************************************
    real :: sigmaAntiNucNuc     ! AntiNucleon + Nucleon final state
    real :: sigmaAntiNucDelta   ! AntiNucleon + Delta final state
    real :: sigmaNucAntiDelta   ! Nucleon + AntiDelta final state


    ! Production cross section pbar+p -> Nbar+N+mesons, hyperon pbar+p -> Y/Ybar+X
    ! and pbar+p -> J/Psi production cross sections  (only for plotting, not for simulations):
    real :: sigmaPROD, sigmaY, sigmaJPsi

    !exclusive production Xsection pbar + p --> Lambda LambdaBar, Lambda Sigma0Bar+c.c. (uncharged),
    !Xi+XiBar and Omega+OmegaBar:
    real :: sigmaYYbar, sigmaYSbar, sigmaXiXiBar, sigmaOmegaBar

    ! Working variables:
    real :: sigma0, pInitial, s, dummy, mDelta, cut, cut2, x!,sigppbarTot, sigppbarEl
    integer :: totCharge, totId, izDelta, nAnti
    logical :: AntiDeltaFlag

    if ( partIn(1)%antiParticle .eqv. partIn(2)%antiParticle ) then
       write(*,*) 'XsectionAntiBarBar is not applicable for bar+bar or antibar+antibar collisions !!!'
       write(*,*) partIn%Id, partIn%antiParticle
       stop
    end if

    pauliIncluded=.false.

    ! Initialize output:
    partOut(:)%ID=0
    partOut(:)%charge=0
    partOut(:)%antiParticle=.false.
    partOut(:)%mass=0.
    sigmaTot=0.
    sigmaElast=0.
    sigmaCEX=0.
    sigmaAnni=0.

    if ( .not.getRMF_flag() ) then
       pInitial=get_pInitial(partIn,0)
    else
       ! Here we take a vacuum expression, since srts has been already
       ! corrected for the in-medium thresholds (see  generateFinalState):
       s= srts**2
       pInitial = (s-partIn(1)%mass**2+partIn(2)%mass**2)**2/(4.*s) - partIn(2)%mass**2
       pInitial = sqrt( max(0.,pInitial) )
    end if

    ! Initializing the cross sections:
    sigmaAntiNucNuc=0.
    sigmaAntiNucDelta=0.
    sigmaNucAntiDelta=0.

    totCharge= sum(partIn(1:2)%charge)
    totId= sum(partIn(1:2)%ID)

    call sigmaBarAntiBar(srts,partIn,mediumAtColl, &
         dummy, sigmaElast, sigmaCEX, sigmaAnni, &
         sigmaYYbar, sigmaYSbar, sigmaXiXiBar, sigmaOmegaBar, &
         sigmaJPsi,sigmaPROD,sigmaY)
    ! Note: sigTotal given by sigmaBarAntiBar is not used here, since the total
    ! cross section sigmaTot will be computed by the direct summation over
    ! included channels later on

    if ( totId.eq.2 ) then    ! nucleon+antinucleon

       sigmaAntiNucNuc= sigmaElast + sigmaCEX

       call NbarN_to_NbarDelta_Integrated(sigma0,srts)
       if ( totCharge.eq.0 ) then
          sigmaAntiNucDelta= 2./3.*sigma0/pInitial
       else
          sigmaAntiNucDelta= 4./3.*sigma0/pInitial
       end if
       sigmaNucAntiDelta= sigmaAntiNucDelta

    else if ( totId.eq.3 ) then     ! antinucleon+delta or nucleon+antidelta

       if ( partIn(1)%Id.eq.delta ) then
          izDelta=  partIn(1)%charge
          mDelta= partIn(1)%mass
          AntiDeltaFlag= partIn(1)%antiParticle
       else
          izDelta=  partIn(2)%charge
          mDelta= partIn(2)%mass
          AntiDeltaFlag= partIn(2)%antiParticle
       end if

       if (AntiDeltaFlag) then
          sigmaAntiNucDelta= 0.      ! Rough assumption
          sigmaNucAntiDelta= sigmaElast + sigmaCEX
       else
          sigmaNucAntiDelta= 0.      ! Rough assumption
          sigmaAntiNucDelta= sigmaElast + sigmaCEX
       end if

       if ( abs(totCharge) < 2 ) then
          call NbarDelta_to_NbarN(sigma0,srts,mDelta)
          if (      (      AntiDeltaFlag .and. ( izDelta.eq.-2 .or. izDelta.eq.1 ) ) &
               &.or. ( .not.AntiDeltaFlag .and. ( izDelta.eq.-1 .or. izDelta.eq.2 ) ) ) then
             sigmaAntiNucNuc= sigma0/pInitial
          else if ( totCharge.eq.0 ) then
             sigmaAntiNucNuc= 2./3.*sigma0/pInitial
          else
             sigmaAntiNucNuc= sigma0/3./pInitial
          end if
       else
          sigmaAntiNucNuc= 0.
       end if

    end if

    ! Test:
    !sigmaAnni=0.
    !sigmaAntiNucNuc=0.

    if (totId <= 3) then
       sigmaTot= sigmaAnni + sigmaAntiNucNuc + sigmaAntiNucDelta + sigmaNucAntiDelta + &
            &    sigmaYYbar + sigmaYSbar + sigmaXiXiBar + sigmaOmegaBar
    else
       sigmaTot= sigmaAnni + sigmaElast
    end if

    if (present(PlotFlag)) then
       if (plotFlag)  call makeOutput
    end if

    ! Choose channel by Monte-Carlo ***************************:

    ! Check that there are open channels:
    if (sigmaTot.lt.0.0001) return

    cut=rn()*sigmaTot

    if ( sigmaAnni >= cut ) then                         ! annihilation
       ! This setting is not the final one: it is done only to signal that
       ! annihilation happen:
       partOut(:)%ID=pion
       partOut(:)%charge=0
       return
    end if

    ! Nonannihilation channels:

    partOut(1:2)%antiParticle=partIn(1:2)%antiparticle
    if (totId > 3) then
       partOut(1:2)%ID=partIn(1:2)%Id
       partOut(1:2)%charge=partIn(1:2)%charge
       return
    end if

    if ( partIn(1)%antiparticle ) then
       nAnti=1
    else
       nAnti=2
    end if

    if ( sigmaAnni + sigmaAntiNucNuc >= cut ) then  ! antinucleon+nucleon production

       partOut(1:2)%ID=(/nucleon,nucleon/)

       if ( partIn(1)%Id.eq.nucleon .and. partIn(2)%Id.eq.nucleon ) then

          ! antinucleon+nucleon collision:

          if ( totCharge == 0 ) then
             cut2=rn()*sigmaAntiNucNuc
             if ( sigmaElast >= cut2 ) then
                partOut(nAnti)%charge=partIn(nAnti)%charge
             else
                partOut(nAnti)%charge= -1 - partIn(nAnti)%charge
             end if
          else
             partOut(nAnti)%charge=partIn(nAnti)%charge
          end if

       else

          ! antinucleon+Delta or antiDelta+nucleon collision:

          if ( totCharge == 0 ) then
             x=rn()
             if ( x < 0.5 ) then
                partOut(nAnti)%charge= -1
             else
                partOut(nAnti)%charge= 0
             end if
          else if ( totCharge == -1 ) then
             partOut(nAnti)%charge= -1
          else if ( totCharge == 1 ) then
             partOut(nAnti)%charge= 0
          else
             write(*,*) 'In XsectionAntiBarBar 1: wrong outgoing channel, totCharge=', totCharge
             stop
          end if

       end if

       partOut(3-nAnti)%charge= totCharge - partOut(nAnti)%charge

       return

    else if ( sigmaAnni + sigmaAntiNucNuc + sigmaAntiNucDelta >= cut ) then  ! antinucleon+Delta production

       partOut(nAnti)%Id=nucleon
       partOut(3-nAnti)%Id=delta

       if ( partIn(1)%Id.eq.nucleon .and. partIn(2)%Id.eq.nucleon ) then

          ! antinucleon+nucleon collision:

          x=rn()

          if ( totCharge == 0 ) then
             if ( x < 0.5 ) then
                partOut(nAnti)%charge= -1
             else
                partOut(nAnti)%charge= 0
             end if
          else if ( totCharge == -1 ) then
             if ( x < 0.25 ) then
                partOut(nAnti)%charge= -1
             else
                partOut(nAnti)%charge= 0
             end if
          else if ( totCharge == 1 ) then
             if ( x < 0.25 ) then
                partOut(nAnti)%charge= 0
             else
                partOut(nAnti)%charge= -1
             end if
          end if

       else if ( partIn(nAnti)%Id.eq.nucleon ) then

          ! antinucleon+delta collision:

          if ( abs(totCharge).eq.2 ) then
             partOut(nAnti)%charge= partIn(nAnti)%charge
          else
             cut2=rn()*sigmaAntiNucDelta
             if ( sigmaElast >= cut2 ) then
                partOut(nAnti)%charge= partIn(nAnti)%charge
             else
                partOut(nAnti)%charge= -1 - partIn(nAnti)%charge
             end if
          end if


       else

          write(*,*) 'In XsectionAntiBarBar 2: wrong outgoing channel'
          write(*,*) 'IdIn1, IdIn2, antiIn1, antiIn2:', partIn%Id, partIn%antiParticle
          write(*,*) 'IdOut1, IdOut2, antiOut1, antiOut2:', partOut%Id, partOut%antiParticle
          stop

       end if

       partOut(3-nAnti)%charge= totCharge - partOut(nAnti)%charge

       return

    else if ( sigmaAnni + sigmaAntiNucNuc + sigmaAntiNucDelta + sigmaNucAntiDelta >= cut ) then  ! nucleon+antiDelta production

       partOut(nAnti)%Id=delta
       partOut(3-nAnti)%Id=nucleon

       if ( partIn(1)%Id.eq.nucleon .and. partIn(2)%Id.eq.nucleon ) then

          ! antinucleon+nucleon collision:

          x=rn()

          if ( totCharge == 0 ) then
             if ( x < 0.5 ) then
                partOut(nAnti)%charge= -1
             else
                partOut(nAnti)%charge= 0
             end if
          else if ( totCharge == -1 ) then
             if ( x < 0.25 ) then
                partOut(nAnti)%charge= -1
             else
                partOut(nAnti)%charge= -2
             end if
          else if ( totCharge == 1 ) then
             if ( x < 0.25 ) then
                partOut(nAnti)%charge= 0
             else
                partOut(nAnti)%charge= 1
             end if
          end if

          partOut(3-nAnti)%charge= totCharge - partOut(nAnti)%charge

          return

       else if ( partIn(nAnti)%Id.eq.delta ) then

          ! antiDelta+nucleon collision:

          if ( abs(totCharge).eq.2 ) then
             partOut(3-nAnti)%charge= partIn(3-nAnti)%charge
          else
             cut2=rn()*sigmaNucAntiDelta
             if ( sigmaElast >= cut2 ) then
                partOut(3-nAnti)%charge= partIn(3-nAnti)%charge
             else
                partOut(3-nAnti)%charge= 1 - partIn(3-nAnti)%charge
             end if
          end if

          partOut(nAnti)%charge= totCharge - partOut(3-nAnti)%charge

          return

       else

          write(*,*) 'In XsectionAntiBarBar 3: wrong outgoing channel'
          write(*,*) 'IdIn1, IdIn2, antiIn1, antiIn2:', partIn%Id, partIn%antiParticle
          write(*,*) 'IdOut1, IdOut2, antiOut1, antiOut2:', partOut%Id, partOut%antiParticle

          stop

       end if

       ! Lambda AntiLambda production
    else if ( sigmaAnni + sigmaAntiNucNuc + sigmaAntiNucDelta + sigmaNucAntiDelta + sigmaYYbar >= cut ) then

       if (totCharge /= 0) then
          write(*,*) 'In XsectionAntiBarBar: wrong charge of initial state'
          write(*,*) 'for Lambda+LambdaBar production : ', totCharge
          stop
       end if

       partOut(1:2)%ID=Lambda
       partOut(1:2)%charge=0

       return

    else if ( sigmaAnni + sigmaAntiNucNuc + sigmaAntiNucDelta + sigmaNucAntiDelta + &
         &   sigmaYYbar + sigmaYSbar >= cut ) then  ! Lambda SigmaBar / LambdaBar Sigma production

       if (abs(totCharge) > 1) then
          write(*,*) 'In XsectionAntiBarBar: wrong charge of initial state'
          write(*,*) 'for  Lambda SigmaBar / LambdaBar Sigma production : ', totCharge
          stop
       end if

       x = rn()
       if (x < 0.5) then
          partOut(nAnti)%ID=Lambda  !LambdaBar Sigma
          partOut(nAnti)%charge=0
          partOut(3-nAnti)%ID=SigmaResonance
          partOut(3-nAnti)%charge=totCharge
       else
          partOut(nAnti)%ID=SigmaResonance   !Lambda SigmaBar
          partOut(nAnti)%charge=totCharge
          partOut(3-nAnti)%ID=Lambda
          partOut(3-nAnti)%charge=0
       end if

       return

    else if ( sigmaAnni + sigmaAntiNucNuc + sigmaAntiNucDelta + sigmaNucAntiDelta + &
         &   sigmaYYbar + sigmaYSbar  + sigmaXiXiBar >= cut ) then  ! Xi+XiBar production

       partOut(1:2)%ID=Xi
       if ( totCharge == 0 ) then
          x = rn()
          if (x < 0.5) then                ! Xi^- + XiBar^+
             partOut(nAnti)%charge=1
          else                             ! Xi^0 + XiBar^0
             partOut(nAnti)%charge=0
          end if
       else if ( totCharge == -1) then     ! Xi^- + XiBar^0
          partOut(nAnti)%charge=0
       else if ( totCharge == 1) then      ! Xi^0 + XiBar^+
          partOut(nAnti)%charge=1
       else
          write(*,*) 'In XsectionAntiBarBar: wrong charge of initial state'
          write(*,*) 'for Xi+XiBar production : ', totCharge
          stop
       end if

       partOut(3-nAnti)%charge= totCharge - partOut(nAnti)%charge
       return

    else if ( sigmaAnni + sigmaAntiNucNuc + sigmaAntiNucDelta + sigmaNucAntiDelta + &
         &   sigmaYYbar + sigmaYSbar  + sigmaXiXiBar + sigmaOmegaBar >= cut ) then  ! Omega+OmegaBar production

       partOut(1:2)%ID=OmegaResonance
       if ( totCharge == 0 ) then
          partOut(nAnti)%charge=1
       else
          write(*,*) 'In XsectionAntiBarBar: wrong charge of initial state'
          write(*,*) 'for Omega+OmegaBar production : ', totCharge
          stop
       end if

       partOut(3-nAnti)%charge= totCharge - partOut(nAnti)%charge
       return


    else

       write(*,*) 'Problem:  final state can not be chosen in XsectionAntiBarBar'
       write(*,*) 'srts: ', srts
       write(*,*) 'IdIn1, IdIn2, antiIn1, antiIn2:', partIn%Id, partIn%antiParticle
       write(*,*) 'sigmaAnni: ', sigmaAnni
       write(*,*) 'sigmaAntiNucNuc: ', sigmaAntiNucNuc
       write(*,*) 'sigmaAntiNucDelta: ', sigmaAntiNucDelta
       write(*,*) 'sigmaNucAntiDelta: ', sigmaNucAntiDelta
       write(*,*) 'sigmaTot: ', sigmaTot
       write(*,*) 'cut: ', cut

       ! enforce elastic scattering:
       partOut(1:2)%ID=partIn(1:2)%Id
       partOut(1:2)%charge=partIn(1:2)%charge
       partOut(1:2)%antiParticle=partIn(1:2)%antiParticle

    end if

  contains
    !**************************************************************************
    !****s* XsectionAntiBarBar/makeOutput
    ! NAME
    ! subroutine makeOutput
    !
    ! PURPOSE
    ! Writes all cross sections to file as function of srts and plab [GeV].
    !
    ! Filename:
    ! * 'AntiBarBar_Tot_Elast.dat'  : sigmaTot, sigmaElast
    ! * 'AntiBarBar_Part.dat'  : sigmaAnni, sigmaAntiNucNuc, sigmaAntiNucDelta, sigmaNucAntiDelta, sigmaYYbar
    ! NOTES
    ! * In the present implementation sigmaTot= sigmaAnni + sigmaAntiNucNuc + sigmaAntiNucDelta
    !   + sigmaNucAntiDelta + sigmaYYbar + sigmaYSbar + sigmaXiXiBar + sigmaOmegaBar
    ! * sigmaElast does not include charge exchange.
    !**************************************************************************
    subroutine makeOutput

      logical, save :: initFlag=.true.

      ! The output files
      character(30), dimension(1:2) :: outputFile
      real :: plab

      outputFile(1)='AntiBarBar_Tot_Elast.dat'
      outputFile(2)='AntiBarBar_Part.dat'

      plab=p_lab(srts,partIn(1)%mass,partIn(2)%mass)

      if (initFlag) then
         open(file=outputFile(1),UNIT=101,Status='Replace',Action='Write')
         open(file=outputFile(2),UNIT=102,Status='Replace',Action='Write')
         write(101,*) '#   srts,    plab,   sigmaTot,  sigmaElast'
         write(102,*) '# column:   quantity:'
         write(102,*) '# 1    srts'
         write(102,*) '# 2    plab'
         write(102,*) '# 3    sigmaElast'
         write(102,*) '# 4    sigmaAnni'
         write(102,*) '# 5    sigmaAntiNucNuc'
         write(102,*) '# 6    sigmaPROD (not used in simulations)'
         write(102,*) '# 7    sigmaY (not used in simulations)'
         write(102,*) '# 8    sigmaAntiNucDelta'
         write(102,*) '# 9    sigmaNucAntiDelta'
         write(102,*) '# 10   sigmaYYbar'
         write(102,*) '# 11   sigmaYSbar'
         write(102,*) '# 12   sigmaXiXiBar'
         write(102,*) '# 12   sigmaOmegaOmegaBar'
         write(102,*) '# 13   sigmaJPsi (used in high-energy mode only)'
         initFlag=.false.
      else
         open(file=outputFile(1),UNIT=101,Status='old',Position='Append',Action='Write')
         open(file=outputFile(2),UNIT=102,Status='old',Position='Append',Action='Write')
      end if

      write(101,'(4(e14.7,1x))') srts, plab, sigmaTot, sigmaElast
      write(102,'(14(e14.7,1x))') srts, plab, sigmaElast, sigmaAnni, sigmaAntiNucNuc,&
           & sigmaPROD, sigmaY, sigmaAntiNucDelta, sigmaNucAntiDelta,&
           & sigmaYYbar, sigmaYSbar, sigmaXiXiBar, sigmaOmegaBar, sigmaJPsi

    end subroutine makeOutput

  end subroutine XsectionAntiBarBar

end module antiBarBar_main
