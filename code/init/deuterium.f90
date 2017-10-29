!******************************************************************************
!****m* /deuterium
! NAME
! module deuterium
!
! PURPOSE
! Includes subroutines to model the momentum distribution of nucleons in
! a Deuterium.
!******************************************************************************
module deuterium

  implicit none
  private

  public :: initDeuterium


  !****************************************************************************
  !****g* deuterium/waveFunction_switch
  ! SOURCE
  !
  integer, save :: waveFunction_switch = 1
  !
  ! PURPOSE
  ! Possible values are:
  ! * 0 -- No Wave functions! Pointlike Deuterium
  ! * 1 -- Wave functions according to Bonn potential
  ! * 2 -- Wave functions according to Argonne V18
  !****************************************************************************

  integer, parameter :: bonn = 1
  integer, parameter :: argonne = 2


  !****************************************************************************
  !****g* deuterium/iParam
  ! SOURCE
  !
  integer, save :: iParam = 1
  !
  ! PURPOSE
  ! Choose parameterization of momentum distribution when using the Bonn potential.
  ! Possible values:
  ! * 1 --   Full Bonn  (MaH87)
  ! * 2 --   OBEPQ      (MaH87)
  ! * 3 --   OBEPQ-A    (Mac89)
  ! * 4 --   OBEPQ-B    (Mac89)
  ! * 5 --   OBEPQ-!    (Mac89)
  ! * 6 --   OBEPR      (MaH87)  self-made
  ! * 7 --   Paris
  ! References:
  ! MaH87: R. Machleidt et al. Phys. Rep. 149, 1 (1987)
  ! Mac89: R. Machleidt, Advances in Nucl. Phys. Vol 19
  !****************************************************************************

  integer, parameter, dimension(7) :: ideuteron = (/2,10,11,12,13,20,40/)


  !****************************************************************************
  !****g* deuterium/pMax
  ! SOURCE
  !
  real,save :: pMax = 0.5
  !
  ! PURPOSE
  ! Cut-off parameter for Fermi momentum
  !****************************************************************************


  !****************************************************************************
  !****g* deuterium/scaleMomentum
  ! SOURCE
  !
  real, save :: scaleMomentum = 1.0
  !
  ! PURPOSE
  ! The selected momentum is multiplied by this factor afterwards, i.e.
  ! some rescaling is done
  !****************************************************************************


  logical, parameter :: argonne_no_pWave=.true.

  logical, save :: inputFlag = .true.


contains


  subroutine readInput
    use output, only: Write_ReadingInput

    integer :: ios

    !**************************************************************************
    !****n* deuterium/deuteriumFermi
    ! NAME
    ! Namelist /deuteriumFermi/
    ! PURPOSE
    ! Includes the parameters:
    ! * waveFunction_switch
    ! * iParam
    ! * pMax
    ! * scaleMomentum
    ! which are used for the distribution of nucleons in deuterium.
    !**************************************************************************
    NAMELIST /deuteriumFermi/ waveFunction_switch, iParam, pMax, scaleMomentum

    call Write_ReadingInput('deuteriumFermi',0)
    rewind(5)
    read(5,nml=deuteriumFermi,iostat=ios)
    call Write_ReadingInput('deuteriumFermi',0,ios)
    write(*,*)
    write(*,*) '  Cut-off parameter for Fermi momentum pMax=',pMax
    select case (waveFunction_switch)
    case (BONN)
       write(*,*) '  Wave functions according to Bonn potential'
       write(*,*) '    * Parameterization of momentum distribution=',iParam
    case (ARGONNE)
       write(*,*) '  Wave functions according to Argonne V18 potential'
    case (0)
       write(*,*) '  No wave functions: Point-like deuterium without Fermi motion'
    case default
       write(*,*) '  Wrong waveFunction_switch! STOP!', waveFunction_switch
       stop
    end select

    write(*,*) '  scaleMomentum = ',scaleMomentum

    call Write_ReadingInput('deuteriumFermi',1)

    inputFlag = .false.

  end subroutine readInput


  !****************************************************************************
  !****s* deuterium/initDeuterium
  ! NAME
  ! subroutine initDeuterium(teilchen,nuc,eventNummer,keepNumber)
  ! PURPOSE
  ! Represents nucleus 'nuc' in phase space by testparticles which
  ! are stored in vector 'teilchen'.
  ! The ordering in the vector teilchen is choosen to be random.
  ! The nucleus must be Deuterium.
  !
  ! NOTES
  ! * momentum distribution implemented by Thomas Falter,
  !   spatial distribution added later
  ! * The first nucleon is initialized at the origin, while the
  !   second nucleon gets its spatial coordinate according
  !   the squared wave function.
  ! * The averaged deuterium radius is <r_d>~=2.0 fm.
  !   The variable r used here is the distance between the two
  !   nucleons. Therefore: <r_d> = 1/2 * sqrt(Int_0^Infty dr r^2 psi^2)
  !
  ! INPUTS
  ! * type(nucleus)                 :: nuc
  ! * type(particle),dimension(:,:) :: teilchen
  ! * integer                       :: eventNummer
  !   -- eventNummer is given as "%event" to any initialized nucleon.
  ! * logical                       :: keepNumber
  !   -- flag whether each testparticle gets unique number or not
  ! OUTPUT
  ! * type(particle),dimension(:,:) :: teilchen
  !****************************************************************************
  subroutine initDeuterium (teilchen, nuc, eventNummer, keepNumber, z_shift)
    use particleDefinition
    use IdTable, only: nucleon
    use nucleusDefinition
    use random, only: rn, rnOmega
    use output, only: Write_InitStatus
    use constants, only: mN
    use deuterium_PL, only: deuterium_pointerList, deuteriumPL_clear, deuteriumPL_getSpectra

    type(tNucleus),pointer :: nuc
    type(particle),target,dimension(:,:),intent(inout) :: teilchen
    integer, intent(in) :: eventNummer
    logical, intent(in) :: keepNumber
    real, intent(in) :: z_shift

    integer :: producedProtons !Counts number of produced protons
    integer :: i,k,index,offset, chargeSave,numberSave

    ! Create the deuterium_pointerList which stores information on the
    if (allocated(deuterium_pointerList)) deallocate(deuterium_pointerList)
    allocate(deuterium_pointerList(lbound(teilchen,dim=1):uBound(teilchen,dim=1)))
    do k=lbound(teilchen,dim=1),ubound(teilchen,dim=1)
       call deuteriumPL_clear(deuterium_pointerList(k))
    end do
    call Write_InitStatus("Deuterium",0)

    if (inputFlag) call readInput

    if (nuc%mass/=2 .or. nuc%charge/=1) then
       write(*,*) 'The input nucleus is no deuterium: ',nuc%mass,nuc%charge
       stop
    end if

    if (keepNumber) then ! e.g. fullensemble

       producedProtons=0
       do k=1,nuc%mass   !Loop over all particles in nucleus
          do i=1,size(teilchen,dim=1)  !Loop over all ensembles
             offset=0
             do !Search for empty space in particle vector
                index=k+offset
                if (index>size(teilchen,dim=2)) then
                   write(*,*) 'Real particle vector too small. Stop in initDeuterium. 1'
                   stop
                else if (teilchen(i,index)%ID > 0) then
                   offset=offset+1
                else
                   exit
                end if
             end do
             if (index>2) write(*,*) 'WARNING: index > 2',index
             call setToDefault(teilchen(i,index)) !set teilchen to its default values
             Teilchen(i,index)%event=eventNummer
             Teilchen(i,index)%ID=Nucleon
             Teilchen(i,index)%antiparticle=.false.
             Teilchen(i,index)%perturbative=.false.
             Teilchen(i,index)%productionTime=0.
             Teilchen(i,index)%mass=mN
             if (i==1) then
                Teilchen(i,index)%charge = chooseCharge(k)
                call setNumber(teilchen(i,index)) ! give each particle a unique number
                chargeSave=teilchen(i,index)%charge
                numberSave=teilchen(i,index)%number
             else
                Teilchen(i,index)%charge= chargeSave
                Teilchen(i,index)%number= numberSave
             end if
             if (k==1) then
                Teilchen(i,index)%position=0.
                deuterium_pointerList(i)%part1 =>Teilchen(i,index)
             else
                deuterium_pointerList(i)%part2 =>Teilchen(i,index)
                Teilchen(i,index)%position = choosePosition() ! argument of wave function is relative distance
                ! The first particle was first chosen at the origin. Now we shift the first and second particle by
                ! half of their relative distance, such that the center of mass is situated at the origin
                Teilchen(i,index)%position=0.5*Teilchen(i,index)%position
                deuterium_pointerList(i)%part1%position=-Teilchen(i,index)%position
                ! Please note: while the relative position is defined as r=r1-r2,
                ! the relative momentum is p=1/2(p1-p2) [The reduced mass is mu=m/2].
                ! Therefore the argument of the wave functionis the proton momentum,
                ! not half of it.
                Teilchen(i,index)%momentum(0:3) = chooseMomentum()
                deuterium_pointerList(i)%part1%momentum(1:3) = -Teilchen(i,index)%momentum(1:3)
                deuterium_pointerList(i)%part1%momentum(0)   =  Teilchen(i,index)%momentum(0)
                call boostIt (nuc, i)
                ! Set velocities
                Teilchen(i,index)%velocity(1:3)=Teilchen(i,index)%momentum(1:3)/Teilchen(i,index)%momentum(0)
                deuterium_pointerList(i)%part1%velocity(1:3) = &
                      deuterium_pointerList(i)%part1%momentum(1:3)/deuterium_pointerList(i)%part1%momentum(0)
             end if
         end do

       end do
       if (producedProtons/=nuc%charge) then
          write(*,*) 'Problem in initDeuterium 1', producedProtons, nuc%charge
       end if


    else
       do i=1,size(teilchen,dim=1)  !Loop over all ensembles
          producedProtons=0
          offset=0
          do k=1,nuc%mass   !Loop over all particles in nucleus
             do ! Search for empty space in particle vector
                index=k+offset
                if (teilchen(i,index)%ID > 0) then
                   offset=offset+1
                else if (index>size(teilchen,dim=2)) then
                   write(*,*) 'Real particle vector too small. Stop in initDeuterium. 2'
                   stop
                else
                   exit
                end if
             end do
             if (index>2) write(*,*) 'WARNING: index > 2',index
             call settoDefault(teilchen(i,index)) !set teilchen to its default values
             call setNumber(teilchen(i,index)) ! give each particle a unique number
             Teilchen(i,index)%charge = chooseCharge(k)
             Teilchen(i,index)%event=eventNummer
             Teilchen(i,index)%ID=Nucleon
             Teilchen(i,index)%antiparticle=.false.
             Teilchen(i,index)%perturbative=.false.
             Teilchen(i,index)%in_Formation=.false.
             Teilchen(i,index)%productionTime=0.
             Teilchen(i,index)%mass=mN
             if (k==1) then
                Teilchen(i,index)%position=0.
                Teilchen(i,index)%momentum=0.
                deuterium_pointerList(i)%part1 =>Teilchen(i,index)
             else
                Teilchen(i,index)%position = choosePosition()  ! argument of wave function is relative distance
                ! The first particle was first chosen at the origin. Now we shift the
                ! first and second particle by  half of their relative distance, such
                ! that the center of mass is situated at the origin:
                Teilchen(i,index)%position=0.5*Teilchen(i,index)%position
                deuterium_pointerList(i)%part1%position=-Teilchen(i,index)%position

                ! Please note: while the relative position is defined as r=r1-r2,
                ! the relative momentum is p=1/2(p1-p2) [The reduced mass is mu=m/2].
                ! Therefore the argument of the wave function is the proton momentum,
                ! not half of it.
                Teilchen(i,index)%momentum(0:3) = chooseMomentum()
                deuterium_pointerList(i)%part1%momentum(1:3) = -Teilchen(i,index)%momentum(1:3)
                deuterium_pointerList(i)%part1%momentum(0)   =  Teilchen(i,index)%momentum(0)
                deuterium_pointerList(i)%part2 =>Teilchen(i,index)
                call boostIt (nuc, i)
                ! Set velocities
                Teilchen(i,index)%velocity(1:3)=Teilchen(i,index)%momentum(1:3)/Teilchen(i,index)%momentum(0)
                deuterium_pointerList(i)%part1%velocity(1:3) = &
                      deuterium_pointerList(i)%part1%momentum(1:3)/deuterium_pointerList(i)%part1%momentum(0)
             end if
          end do
          if (producedProtons/=nuc%charge) then
             write(*,*) 'Problem in initDeuterium 2', producedProtons, nuc%charge
             stop
          end if
       end do
    end if

    open(100,file="RelMom_deuteriumInit.dat")
    open(101,file="RelPos_deuteriumInit.dat")
    call deuteriumPL_getSpectra(deuterium_pointerList,100,101,0.)
    close(100)
    close(101)
    call Write_InitStatus("Deuterium",1)

  contains

    integer function chooseCharge (k)
      integer, intent(in) :: k
      real :: probabilityProton  !probability to produce a proton
      !Choose randomly the charge of the nucleon in the deuterium target
      probabilityProton=float(nuc%charge-producedProtons)/float(nuc%mass-k+1)
      if (rn()<=probabilityProton) then
        chooseCharge=1
        producedProtons=producedProtons+1
      else
        chooseCharge=0
      end if
    end function chooseCharge

    function chooseMomentum () result (momentum)
      use constants, only: mN
      real :: p
      real,dimension(0:3) :: momentum

      if (nuc%fermiMotion) then
        p = deuteriumFermi() * scaleMomentum
      else
        p = 0.
      end if
      momentum(1:3)=p*rnOmega()
      ! assume vacuum dispersion relation:
      momentum(0)=Sqrt(dot_product(momentum(1:3),momentum(1:3))+mN**2)
    end function chooseMomentum

    subroutine boostIt (nucl, i)
        use lorentzTrafo, only: lorentz
        use inputGeneral, only: eventtype
        use eventtypes, only: HeavyIon

        type(tnucleus), intent(in) :: nucl
        integer :: i

        !Do nothing if nucleus rests in calculation frame :
        if (sum(abs(nucl%velocity))<0.000001) then
           ! Write(*,*) 'No boost necessary! Velocity of nucleus:',nucl%velocity
           ! Shift to actual position of center of mass of nucleus
           deuterium_pointerList(i)%part1%position(1:3)=deuterium_pointerList(i)%part1%position(1:3)+nucl%position
           deuterium_pointerList(i)%part2%position(1:3)=deuterium_pointerList(i)%part2%position(1:3)+nucl%position
        else

           if (abs(nucl%velocity(1))>0.1 .or. Abs(nucl%velocity(2))>0.1) then
              write(*,*) 'Error in deuterium::boostIt'
              write(*,*) 'Velocity of nucleus in transverse direction is not negligible'
              write(*,*) 'Velocity in x,y:' , nucl%velocity(1:2)
              write(*,*) 'Nucleus: A=', nucl%mass,'Z=', nucl%charge
              write(*,*) 'critical error. Stop program!'
              stop
           end if

           ! boost members of nucleus from restframe to calculation frame,
           ! which is moving with -beta seen from the rest frame:
           call lorentz(-nucl%velocity,deuterium_pointerList(i)%part1%momentum)
           call lorentz(-nucl%velocity,deuterium_pointerList(i)%part2%momentum)

           ! Length contraction of nucleus in z-direction(linear scaling)
           ! Neglect length contraction in transverse direction
           deuterium_pointerList(i)%part1%position(3)=deuterium_pointerList(i)%part1%position(3)*sqrt(1-nucl%velocity(3)**2)
           deuterium_pointerList(i)%part2%position(3)=deuterium_pointerList(i)%part2%position(3)*sqrt(1-nucl%velocity(3)**2)

           ! Shift to actual position of center of mass of nucleus
           deuterium_pointerList(i)%part1%position(1:3)=deuterium_pointerList(i)%part1%position(1:3)+nucl%position
           deuterium_pointerList(i)%part2%position(1:3)=deuterium_pointerList(i)%part2%position(1:3)+nucl%position
        end if

        ! shift particles along z-axis to avoid overlapping between projectile & target:
        if (eventType==HeavyIon) then
          if (eventNummer==1) then
            deuterium_pointerList(i)%part1%position(3)=deuterium_pointerList(i)%part1%position(3)+z_shift
            deuterium_pointerList(i)%part2%position(3)=deuterium_pointerList(i)%part2%position(3)+z_shift
          else
            deuterium_pointerList(i)%part1%position(3)=deuterium_pointerList(i)%part1%position(3)-z_shift
            deuterium_pointerList(i)%part2%position(3)=deuterium_pointerList(i)%part2%position(3)-z_shift
          end if
        end if

    end subroutine boostIt

    function choosePosition() result (r)
      use argonneV18, only: argonne_WF_rSpace

      real, dimension(1:3) :: r

      real, parameter :: maxDist = 20.0 ! very loosely bound system !
      real :: rAbs, Psi2, ws, wd

      if ((WaveFunction_switch==bonn .and. iParam==0) .or. WaveFunction_switch==0) then
         r=0.
         return
      end if

      do
         rAbs=rn()*maxDist

         if (WaveFunction_switch==bonn) then
            call wvfct_r_deuteron(rAbs,ws,wd,ideuteron(iParam))
            psi2 = (ws**2+wd**2)
         else if (WaveFunction_switch==argonne) then
            psi2= argonne_WF_rSpace(rAbs,argonne_no_pWave)
         else
            write(*,*) 'Wrong waveFunction_switch! STOP',waveFunction_switch
         end if
         if (psi2>1.0) write(*,*) 'Warning in deuterium/choosePosition: psi2 =',psi2,'> 1.0'
         if (rn()*1.0<psi2) exit ! the dummy 1.0 should be okay
      end do

      r=rabs*rnOmega()

    end function choosePosition

  end subroutine initDeuterium



  function deuteriumFermi() result(p)
    use argonneV18, only: argonne_max_kSpace, argonne_WF_kSpace
    use random, only: rn

    real :: p

    real,parameter :: a=0.016,b=0.008  !Parameters of comparing function f(x)
    real,save :: Area
    real :: xRan,ARan,yRan,yMax,psi2
    double precision :: ws,wd,xTest
    logical , save :: initFlag =.true.

    select case (waveFunction_switch)
    case (0)
       p=0.
       return
    case (BONN)
       if (iParam<1 .or. iParam>7) then
          write(*,*) 'Error in deuteriumFermi: iParam=',iParam
          stop
       else
          if (initFlag) then
             Area=(a-f(pMax*1000.))/b !convert to MeV
             initFlag=.false.
          end if

          MonteCarlo_bonn : do
             ARan=Area*rn()
             xRan=PofA(ARan)
             xTest=xRan
             call wvfct_p_deuteron(xTest,ws,wd,ideuteron(iParam))
             yMax=xRan**2*real(ws**2+wd**2)
             yRan=f(xRan)*rn()
             if (yRan.le.yMax) exit MonteCarlo_bonn
          end do MonteCarlo_bonn
          p=xRan/1000.!GeV
          return
       end if
    case (argonne)
       MonteCarlo_argonne : do
          p=rn()*pmax
          psi2= argonne_WF_kSpace(p,argonne_no_pWave)*p**2
          ! Assume that "argonne_max_kSpace*1.6" is maximum of distribution (1.6 is just a dummy)
          if (argonne_max_kSpace*1.6.lt.psi2) then
             write(*,*) 'WARNING in deuterium/deuteriumFermi',argonne_max_kSpace,psi2
          end if
          if (rn()*argonne_max_kSpace*1.6.lt.psi2) exit MonteCarlo_argonne
       end do MonteCarlo_argonne
       return
    case default
       write(*,*) 'Wrong waveFunction_switch! STOP',waveFunction_switch
    end select

  contains

    real function f(x)  ! comparing function for Monte Carlo
      real,intent(in) :: x
      f=a*exp(-b*x)
    end function f

    real function PofA(y)
      real,intent(in) :: y
      PofA=-log(1-y*b/a)/b
    end function PofA

  end function deuteriumFermi


  subroutine wvfct_p_deuteron(p,ws,wd,ideuteron)
    !**************************************************************************
    !     subroutine for calculating the deuteron s- and d-waves
    !     in momentum space using a Yukawa parametrization.
    !     Before calling this subroutine one has to initialize the
    !     appropriate parameter values by calling the subroutine cdmj(ideuteron)
    !     where "ideuteron" characterizes the potential model for the deuteron:
    !
    !      ideuteron=2:   Full Bonn  (MaH87)
    !      ideuteron=10:  OBEPQ      (MaH87)
    !      ideuteron=11:  OBEPQ-A    (Mac89)
    !      ideuteron=12:  OBEPQ-B    (Mac89)
    !      ideuteron=13:  OBEPQ-!    (Mac89)
    !      ideuteron=20:  OBEPR      (MaH87)  self-made
    !      ideuteron=40:  Paris
    !      MaH87: R. Machleidt et al. Phys. Rep. 149, 1 (1987)
    !      Mac89: R. Machleidt, Advances in Nucl. Phys. Vol 19
    !
    !      ws: S-wave of the deuteron at momentum p
    !      wd: D-wave of the deuteron at momentum p
    !      in Yukawa parametrization Mac89, equation D.26
    !      dimension of p [MeV], of wave function [MeV^(-3/2)]
    !
    !      normalization:
    !      integral from 0 -> infty dp*p**2*(ws**2+ wd**2) = 1
    !      consider  Mac89, equation D.12
    !
    !**************************************************************************

    double precision,intent(in) :: p
    double precision,intent(out) :: ws,wd
    integer,intent(in) :: ideuteron

    logical,save :: initFlag=.true.

    !used in cdmj
    double precision,dimension(20),save :: c,d,xm2
    integer,save :: jmax

    double precision :: PI,factor,x
    integer :: j

    if (initFlag) then
       call cdmj
       initFlag=.false.
    end if

    PI=4D0*DATAN(1D0)
    factor= dsqrt(2d0/pi)

    ws=0.d0
    wd=0.d0
    do j=1,jmax
       x=1d0/(xm2(j)+p**2)
       ws=ws+x*c(j)
       wd=wd+x*d(j)
    end do

    ws=ws*factor
    wd=wd*factor

  contains

    subroutine cdmj

      double precision,parameter :: HQC=197.3286D0
      double precision :: alpha,xm0,Aconst,Bconst,Cconst

      if (IDEUTERON.EQ.2) then !     full model
         JMAX=11
         ALPHA=.231609D0
         C( 1)= .90457337D+0
         C( 2)=-.35058661D+0
         C( 3)=-.17635927D+0
         C( 4)=-.10418261D+2
         C( 5)= .45089439D+2
         C( 6)=-.14861947D+3
         C( 7)= .31779642D+3
         C( 8)=-.37496518D+3
         C( 9)= .22560032D+3
         C(10)=-.54858290D+2
         D( 1)= .24133026D-1
         D( 2)=-.64430531D+0
         D( 3)= .51093352D+0
         D( 4)=-.54419065D+1
         D( 5)= .15872034D+2
         D( 6)=-.14742981D+2
         D( 7)= .44956539D+1
         D( 8)=-.71152863D-1
      else if (IDEUTERON.EQ.10) THEN   ! OBEPQ [MaH87]
         JMAX=11
         ALPHA=  0.231609d0
         C( 1)=  0.88628672d0
         C( 2)= -0.27591814d0
         C( 3)= -0.11610727d0
         C( 4)= -0.12975243e+02
         C( 5)=  0.77490155e+02
         C( 6)= -0.27298039e+03
         C( 7)=  0.53402693e+03
         C( 8)= -0.56328069e+03
         C( 9)=  0.30214616e+03
         C(10)= -0.64920925e+02
         D( 1)=  0.23237078E-01
         D( 2)= -0.52115578d0
         D( 3)= -0.57197401d0
         D( 4)=  0.27570246e+01
         D( 5)= -0.26157324e+02
         D( 6)=  0.84419883e+02
         D( 7)= -0.98308997e+02
         D( 8)=  0.38498490e+02

      else if (IDEUTERON.EQ.11) THEN   ! OBEPQA [Mac89]
         JMAX=11
         ALPHA=  0.231607d0
         C( 1)=  0.88681402e+00
         C( 2)= -0.27176295e+00
         C( 3)= -0.38234310e+00
         C( 4)= -0.97399200e+01
         C( 5)=  0.57873078e+02
         C( 6)= -0.21112738e+03
         C( 7)=  0.42789416e+03
         C( 8)= -0.46272723e+03
         C( 9)=  0.25255966e+03
         C(10)= -0.54964903e+02
         D( 1)=  0.23345605e-01
         D( 2)= -0.57467557e+00
         D( 3)=  0.92159360e+00
         D( 4)= -0.10072048e+02
         D( 5)=  0.21821344e+02
         D( 6)= -0.34389664e+01
         D( 7)= -0.20707396e+02
         D( 8)=  0.12048237e+02
      else if (IDEUTERON.EQ.12) THEN   ! OBEPQB [Mac89]
         JMAX=11
         ALPHA=  0.231607d0
         C( 1)=  0.88611410E+00
         C( 2)= -0.24885006E+00
         C( 3)= -0.88346659E+00
         C( 4)= -0.46847106E+01
         C( 5)=  0.34755263E+02
         C( 6)= -0.16379524E+03
         C( 7)=  0.38880024E+03
         C( 8)= -0.46566577E+03
         C( 9)=  0.27495507E+03
         C(10)= -0.64119028E+02
         D( 1)=  0.23437728E-01
         D( 2)= -0.54665750E+00
         D( 3)=  0.51669408E+00
         D( 4)= -0.73905273E+01
         D( 5)=  0.16323355E+02
         D( 6)= -0.34932110E+01
         D( 7)= -0.12845278E+02
         D( 8)=  0.74194734E+01
      else if (IDEUTERON.EQ.13) THEN   ! OBEPQC [Mac89]
         JMAX=11
         ALPHA=  0.231607d0
         C( 1)=  0.88507948E+00
         C( 2)= -0.24105451E+00
         C( 3)= -0.10338683E+01
         C( 4)= -0.29885428E+01
         C( 5)=  0.25258598E+02
         C( 6)= -0.13992344E+03
         C( 7)=  0.36051215E+03
         C( 8)= -0.45277411E+03
         C( 9)=  0.27676633E+03
         C(10)= -0.66461680E+02
         D( 1)=  0.23550301E-01
         D( 2)= -0.52404123E+00
         D( 3)=  0.15311637E+00
         D( 4)= -0.50123809E+01
         D( 5)=  0.11340227E+02
         D( 6)= -0.23474968E+01
         D( 7)= -0.81817727E+01
         D( 8)=  0.45534069E+01
      else if (IDEUTERON.EQ.20) THEN   ! OBEPR [MaH87]
         JMAX=11
         ALPHA=  0.23160176E+00
         C( 1)=  0.88610147E+00
         C( 2)= -0.24390698E+00
         C( 3)= -0.73362772E+00
         C( 4)= -0.63246271E+01
         C( 5)=  0.52336822E+02
         C( 6)= -0.24940260E+03
         C( 7)=  0.59575224E+03
         C( 8)= -0.75469583E+03
         C( 9)=  0.52716666E+03
         C(10)= -0.19708877E+03
         D( 1)=  0.23018785E-01
         D( 2)= -0.55181025E+00
         D( 3)=  0.17221384E+01
         D( 4)= -0.31920298E+02
         D( 5)=  0.22267538E+03
         D( 6)= -0.88888448E+03
         D( 7)=  0.21410284E+04
         D( 8)= -0.30936406E+04

      else if (IDEUTERON.eq.40) THEN  ! original PARIS
         JMAX=13
         ALPHA=  0.23162461D+00   ! T=0-nucleon mass
         C( 1)=  0.88688076D+00
         C( 2)= -0.34717093D+00
         C( 3)= -0.30502380D+01
         C( 4)=  0.56207766D+02
         C( 5)= -0.74957334D+03
         C( 6)=  0.53365279D+04
         C( 7)= -0.22706863D+05
         C( 8)=  0.60434469D+05
         C( 9)= -0.10292058D+06
         C(10)=  0.11223357D+06
         C(11)= -0.75925226D+05
         C(12)=  0.29059715D+05
         D( 1)=  0.23135193D-01
         D( 2)= -0.85604572D+00
         D( 3)=  0.56068193D+01
         D( 4)= -0.69462922D+02
         D( 5)=  0.41631118D+03
         D( 6)= -0.12546621D+04
         D( 7)=  0.12387830D+04
         D( 8)=  0.33739172D+04
         D( 9)= -0.13041151D+05
         D(10)=  0.19512524D+05

      else
         write(*,*) ' ERROR CDMJ/IDEUTERON'
         STOP
      end if

      XM0=.9D0
      if (ideuteron.eq.40)xm0=1d0
      do J=1,JMAX
         XM2(J)=(ALPHA+(J-1D0)*XM0)**2
      end do

      C(JMAX)=0D0
      do J=1,JMAX-1
         C(JMAX)=C(JMAX)-C(J)
      end do

      ACONST=0D0
      BCONST=0D0
      CCONST=0D0
      do J=1,JMAX-3
         ACONST=ACONST+D(J)/XM2(J)
         BCONST=BCONST+D(J)
         CCONST=CCONST+D(J)*XM2(J)
      end do
      D(jmax-2)=XM2(jmax-2)/(XM2(jmax)-XM2(jmax-2))/(XM2(jmax-1)-XM2(jmax-2))*(-XM2(jmax-1)*XM2(jmax)*ACONST &
           & +(XM2(jmax-1)+XM2(jmax))*BCONST-CCONST)
      D(jmax-1)=XM2(jmax-1)/(XM2(jmax-2)-XM2(jmax-1))/(XM2(jmax)-XM2(jmax-1))*(-XM2(jmax-2)*XM2(jmax)*ACONST &
           & +(XM2(jmax-2)+XM2(jmax))*BCONST-CCONST)
      D(jmax)=XM2(jmax)/(XM2(jmax-1)-XM2(jmax))/(XM2(jmax-2)-XM2(jmax))*(-XM2(jmax-2)*XM2(jmax-1)*ACONST &
           & +(XM2(jmax-2)+XM2(jmax-1))*BCONST-CCONST)
      do J=1,JMAX
         C(J)=C(J)*DSQRT(HQC)
         D(J)=D(J)*DSQRT(HQC)
         XM2(J)=XM2(J)*HQC**2
      end do
    END subroutine cdmj
  end subroutine wvfct_p_deuteron


  subroutine wvfct_r_deuteron(r,ws,wd,ideuteron)
    !**************************************************************************
    !     subroutine for calculating the deuteron s- and d-waves
    !     in r space using a Yukawa parametrization.
    !
    !     please note: r is the relative distance of proton and neutron
    !
    !      dimension of r [fm], of wave function [fm^(-1/2)]
    !
    !      normalization:
    !      integral from 0 -> infty dr*(ws**2+ wd**2) = 1
    !**************************************************************************

    real,intent(in) :: r
    real,intent(out) :: ws,wd
    integer,intent(in) :: ideuteron

    logical,save :: initFlag=.true.
    double precision,parameter :: HQC=197.3286D0
    !used in cdmj
    double precision,dimension(20),save :: c,d,xm2
    integer,save :: jmax

    double precision :: x,mjr
    integer :: j

    if (initFlag) then
       call cdmj
       initFlag=.false.
    end if

    ws=0.d0
    wd=0.d0
    do j=1,jmax
       mjr = sqrt(xm2(j))*r/HQC
       x = exp(-mjr)
       ws = ws+x*c(j)
       wd = wd+x*d(j)*(1.+3./mjr+3./mjr**2)
    end do

    ws = ws/sqrt(HQC)
    wd = wd/sqrt(HQC)

  contains

    subroutine cdmj
      ! this routine is the same as above (not nice, but convenient)

      double precision :: alpha,xm0,Aconst,Bconst,Cconst

      if (IDEUTERON.EQ.2) then !     full model
         JMAX=11
         ALPHA=.231609D0
         C( 1)= .90457337D+0
         C( 2)=-.35058661D+0
         C( 3)=-.17635927D+0
         C( 4)=-.10418261D+2
         C( 5)= .45089439D+2
         C( 6)=-.14861947D+3
         C( 7)= .31779642D+3
         C( 8)=-.37496518D+3
         C( 9)= .22560032D+3
         C(10)=-.54858290D+2
         D( 1)= .24133026D-1
         D( 2)=-.64430531D+0
         D( 3)= .51093352D+0
         D( 4)=-.54419065D+1
         D( 5)= .15872034D+2
         D( 6)=-.14742981D+2
         D( 7)= .44956539D+1
         D( 8)=-.71152863D-1
      else if (IDEUTERON.EQ.10) THEN   ! OBEPQ [MaH87]
         JMAX=11
         ALPHA=  0.231609d0
         C( 1)=  0.88628672d0
         C( 2)= -0.27591814d0
         C( 3)= -0.11610727d0
         C( 4)= -0.12975243e+02
         C( 5)=  0.77490155e+02
         C( 6)= -0.27298039e+03
         C( 7)=  0.53402693e+03
         C( 8)= -0.56328069e+03
         C( 9)=  0.30214616e+03
         C(10)= -0.64920925e+02
         D( 1)=  0.23237078E-01
         D( 2)= -0.52115578d0
         D( 3)= -0.57197401d0
         D( 4)=  0.27570246e+01
         D( 5)= -0.26157324e+02
         D( 6)=  0.84419883e+02
         D( 7)= -0.98308997e+02
         D( 8)=  0.38498490e+02

      else if (IDEUTERON.EQ.11) THEN   ! OBEPQA [Mac89]
         JMAX=11
         ALPHA=  0.231607d0
         C( 1)=  0.88681402e+00
         C( 2)= -0.27176295e+00
         C( 3)= -0.38234310e+00
         C( 4)= -0.97399200e+01
         C( 5)=  0.57873078e+02
         C( 6)= -0.21112738e+03
         C( 7)=  0.42789416e+03
         C( 8)= -0.46272723e+03
         C( 9)=  0.25255966e+03
         C(10)= -0.54964903e+02
         D( 1)=  0.23345605e-01
         D( 2)= -0.57467557e+00
         D( 3)=  0.92159360e+00
         D( 4)= -0.10072048e+02
         D( 5)=  0.21821344e+02
         D( 6)= -0.34389664e+01
         D( 7)= -0.20707396e+02
         D( 8)=  0.12048237e+02
      else if (IDEUTERON.EQ.12) THEN   ! OBEPQB [Mac89]
         JMAX=11
         ALPHA=  0.231607d0
         C( 1)=  0.88611410E+00
         C( 2)= -0.24885006E+00
         C( 3)= -0.88346659E+00
         C( 4)= -0.46847106E+01
         C( 5)=  0.34755263E+02
         C( 6)= -0.16379524E+03
         C( 7)=  0.38880024E+03
         C( 8)= -0.46566577E+03
         C( 9)=  0.27495507E+03
         C(10)= -0.64119028E+02
         D( 1)=  0.23437728E-01
         D( 2)= -0.54665750E+00
         D( 3)=  0.51669408E+00
         D( 4)= -0.73905273E+01
         D( 5)=  0.16323355E+02
         D( 6)= -0.34932110E+01
         D( 7)= -0.12845278E+02
         D( 8)=  0.74194734E+01
      else if (IDEUTERON.EQ.13) THEN   ! OBEPQC [Mac89]
         JMAX=11
         ALPHA=  0.231607d0
         C( 1)=  0.88507948E+00
         C( 2)= -0.24105451E+00
         C( 3)= -0.10338683E+01
         C( 4)= -0.29885428E+01
         C( 5)=  0.25258598E+02
         C( 6)= -0.13992344E+03
         C( 7)=  0.36051215E+03
         C( 8)= -0.45277411E+03
         C( 9)=  0.27676633E+03
         C(10)= -0.66461680E+02
         D( 1)=  0.23550301E-01
         D( 2)= -0.52404123E+00
         D( 3)=  0.15311637E+00
         D( 4)= -0.50123809E+01
         D( 5)=  0.11340227E+02
         D( 6)= -0.23474968E+01
         D( 7)= -0.81817727E+01
         D( 8)=  0.45534069E+01
      else if (IDEUTERON.EQ.20) THEN   ! OBEPR [MaH87]
         JMAX=11
         ALPHA=  0.23160176E+00
         C( 1)=  0.88610147E+00
         C( 2)= -0.24390698E+00
         C( 3)= -0.73362772E+00
         C( 4)= -0.63246271E+01
         C( 5)=  0.52336822E+02
         C( 6)= -0.24940260E+03
         C( 7)=  0.59575224E+03
         C( 8)= -0.75469583E+03
         C( 9)=  0.52716666E+03
         C(10)= -0.19708877E+03
         D( 1)=  0.23018785E-01
         D( 2)= -0.55181025E+00
         D( 3)=  0.17221384E+01
         D( 4)= -0.31920298E+02
         D( 5)=  0.22267538E+03
         D( 6)= -0.88888448E+03
         D( 7)=  0.21410284E+04
         D( 8)= -0.30936406E+04

      else if (IDEUTERON.eq.40) THEN  ! original PARIS
         JMAX=13
         ALPHA=  0.23162461D+00   ! T=0-nucleon mass
         C( 1)=  0.88688076D+00
         C( 2)= -0.34717093D+00
         C( 3)= -0.30502380D+01
         C( 4)=  0.56207766D+02
         C( 5)= -0.74957334D+03
         C( 6)=  0.53365279D+04
         C( 7)= -0.22706863D+05
         C( 8)=  0.60434469D+05
         C( 9)= -0.10292058D+06
         C(10)=  0.11223357D+06
         C(11)= -0.75925226D+05
         C(12)=  0.29059715D+05
         D( 1)=  0.23135193D-01
         D( 2)= -0.85604572D+00
         D( 3)=  0.56068193D+01
         D( 4)= -0.69462922D+02
         D( 5)=  0.41631118D+03
         D( 6)= -0.12546621D+04
         D( 7)=  0.12387830D+04
         D( 8)=  0.33739172D+04
         D( 9)= -0.13041151D+05
         D(10)=  0.19512524D+05

      else
         write(*,*) ' ERROR CDMJ/IDEUTERON'
         STOP
      end if

      XM0=.9D0
      if (ideuteron.eq.40)xm0=1d0
      do J=1,JMAX
         XM2(J)=(ALPHA+(J-1D0)*XM0)**2
      end do

      C(JMAX)=0D0
      do J=1,JMAX-1
         C(JMAX)=C(JMAX)-C(J)
      end do

      ACONST=0D0
      BCONST=0D0
      CCONST=0D0
      do J=1,JMAX-3
         ACONST=ACONST+D(J)/XM2(J)
         BCONST=BCONST+D(J)
         CCONST=CCONST+D(J)*XM2(J)
      end do
      D(jmax-2)=XM2(jmax-2)/(XM2(jmax)-XM2(jmax-2))/(XM2(jmax-1)-XM2(jmax-2))*(-XM2(jmax-1)*XM2(jmax)*ACONST &
           & +(XM2(jmax-1)+XM2(jmax))*BCONST-CCONST)
      D(jmax-1)=XM2(jmax-1)/(XM2(jmax-2)-XM2(jmax-1))/(XM2(jmax)-XM2(jmax-1))*(-XM2(jmax-2)*XM2(jmax)*ACONST &
           & +(XM2(jmax-2)+XM2(jmax))*BCONST-CCONST)
      D(jmax)=XM2(jmax)/(XM2(jmax-1)-XM2(jmax))/(XM2(jmax-2)-XM2(jmax))*(-XM2(jmax-2)*XM2(jmax-1)*ACONST &
           & +(XM2(jmax-2)+XM2(jmax-1))*BCONST-CCONST)
      do J=1,JMAX
         C(J)=C(J)*DSQRT(HQC)
         D(J)=D(J)*DSQRT(HQC)
         XM2(J)=XM2(J)*HQC**2
      end do
    END subroutine cdmj
  end subroutine wvfct_r_deuteron




end module deuterium
