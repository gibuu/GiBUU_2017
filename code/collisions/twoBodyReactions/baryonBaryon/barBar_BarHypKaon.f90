!******************************************************************************
!****m* /barBar_BarHypKaon
! NAME
! module barBar_BarHypKaon
! PURPOSE
! Administrates the calculation of cross sections for strangeness production
! from BB->BYK collisions (with B=N,Delta and Y=Lambda,Sigma,Lambda*,Sigma*).
! The properties of the final-state channels are calculated as well.
! NOTES
! * All isospin channels are included.
! * Ref.: K. Tsushima et al., PRC59 (1999) 369.
! * The original Tsushima model was extended and updated in arXiv:1404.7011.
! * Useful for reactions near and above threshold (up to sqrt(s)=3.6 GeV), not
!   appropriate for sqrt(s) far below threshold (cross section does not fit
!   well existing data & perturbative treatment of kaon production necessary).
!******************************************************************************
module barBar_BarHypKaon

  use preEventDefinition
  use CallStack, only: Traceback

  implicit none
  private

  integer, parameter :: Nch = 30   ! number of primary channels

  !****************************************************************************
  !****g* barBar_BarHypKaon/enable
  ! SOURCE
  !
  logical, save :: enable = .true.
  ! PURPOSE
  ! Enable the production of BB -> B Hyperon Kaon channels.
  ! B=Nucleon^{0,1},Delta^{-,0.+,++}; Hyperon=Lambda^{0},Sigma^{0,-,+}; Kaon=K^{+,0}
  !****************************************************************************

  !****************************************************************************
  !****g* barBar_BarHypKaon/parameter_set
  ! SOURCE
  integer, save :: parameter_set = 2
  ! PURPOSE
  ! Select a particular parameter set for BB->BYK collisions.
  ! Possible values:
  ! * 1 = original Tsushima model: Tsushima et al., PRC59 (1999) 369
  ! * 2 = extended/adjusted model, fitted to HADES data:
  !       Agakishiev et al., arXiv:1404.7011
  ! * 3 = custom parameters based on Tsushima values (as given by the
  !       array 'a' in the jobcard; those values not given in the jobcard
  !       are adopted from Tsushima, i.e. parameter set 1)
  ! * 4 = custom parameters based on HADES values (as given by the
  !       array 'a' in the jobcard; those values not given in the jobcard
  !       are adopted from HADES, i.e. parameter set 2)
  !****************************************************************************

  !****************************************************************************
  !****g* barBar_BarHypKaon/a
  ! SOURCE
  real, dimension(1:Nch), save :: a = -1.
  ! PURPOSE
  ! This array contains the "a" parameters (in microbarn) for the
  ! 30 primary channels, see:
  ! * Tsushima et al., PRC59 (1999) 369, table III
  ! * Agakishiev et al., arXiv:1404.7011, chapter IV
  ! Note: The values given in the jobcard are only used for
  ! parameter_set = 3 and 4.
  !****************************************************************************


  !****************************************************************************
  ! The table with the threshold and cross section parameters for the primary
  ! channels of kaon production from BB collisions:
  !      * 1-28 are from the original Tsushima model,
  !      * 29-30 are an extension to DeltaY*K channels (by K. Lapidus).
  !      * a(2,4,5,6,8,9,29,30) have been modified in arXiv:1404.7011.
  !****************************************************************************
  real, parameter :: stresh(1:Nch) = (/ &   ! threshold in free space for each channel (GeV**2)
                6.504, 6.504, 6.904, 6.904, 6.904, 6.904, 6.904,  &
                8.085, 8.531, 6.504, 6.904, 8.085, 8.085, 8.085,  &
                8.531, 8.531, 8.531, 8.531, 8.531, 8.531, 8.531,  &
                8.085, 8.085, 8.085, 8.531, 8.531, 8.531, 8.531,  &
                9.356, 8.889 /)
  real, parameter :: a_Tsushima(1:Nch) = (/ &
                1.879, 2.812, 5.321, 7.079, 6.310, 11.02, 1.466,  &
                6.166, 10.00, 8.337, 52.72, 2.704, 0.312, 2.917,  &
                10.33, 2.128, 10.57, 10.30, 1.112, 10.62, 0.647,  &
                1.054, 0.881, 0.291, 3.532, 7.047, 2.931, 5.861,  &
                0.000, 0.000 /)
  real, parameter :: a_HADES(1:Nch) = (/ &
                1.879, 1.406, 5.321, 4.955, 3.155, 5.510, 1.466,  &
                2.466, 7.000, 8.337, 52.72, 2.704, 0.312, 2.917,  &
                10.33, 2.128, 10.57, 10.30, 1.112, 10.62, 0.647,  &
                1.054, 0.881, 0.291, 3.532, 7.047, 2.931, 5.861,  &
                8.500, 3.100 /)
  real, parameter :: b(1:Nch) = (/ &
                2.176, 2.121, 2.753, 2.760, 2.773, 2.782, 2.743,  &
                2.842, 2.874, 2.227, 2.799, 2.303, 2.110, 2.350,  &
                2.743, 2.843, 2.757, 2.748, 2.846, 2.759, 2.830,  &
                2.149, 2.150, 2.148, 2.953, 2.952, 2.952, 2.952,  &
                2.842, 2.874 /)
  real, parameter :: c(1:Nch) = (/ &
                5.264, 4.893, 8.510, 8.164, 7.820, 7.674, 3.271,  &
                1.960, 2.543, 2.511, 6.303, 5.551, 2.165, 6.557,  &
                8.915, 5.986, 10.11, 9.321, 5.943, 10.20, 3.862,  &
                7.969, 7.977, 7.934, 12.06, 12.05, 12.03, 12.04,  &
                1.960, 2.543 /)


  ! ID's & charges of final states.
  ! 1-index: Number of isospin channels
  ! 2-index: number of particle in 3-body final state (Baryon,Hyperon & Kaon)
  Integer, dimension(1:18,1:3), save :: IdsOut, ChargesOut

  logical, save :: first = .true.

  public :: barBar_barBarMeson_strange, get_Channels_BYK

 contains

  subroutine readInput

    use output, only: Write_ReadingInput

    integer :: ios, i

    !**************************************************************************
    !****n* barBar_BarHypKaon/BB_BYK
    ! NAME
    ! NAMELIST BB_BYK
    ! PURPOSE
    ! Includes the switches:
    ! * enable
    ! * parameter_set
    ! * a
    !**************************************************************************
    NAMELIST /BB_BYK/ enable, parameter_set, a

    call Write_ReadingInput('BB_BYK',0)
    rewind(5)
    read(5,nml=BB_BYK,iostat=ios)
    call Write_ReadingInput('BB_BYK',0,ios)
    write(*,*) "enable:        ", enable
    write(*,*) "parameter set: ", parameter_set
    write(*,*) "modified 'a' parameters (channel nr, value used, default value, ratio):"
    select case (parameter_set)
    case (1)
      a = a_Tsushima
    case (2)
      a = a_HADES
    case (3)
      do i=1,Nch
        if (a(i)<0.) then
          a(i) = a_Tsushima(i)
        else
          write(*,'(i3,3f8.3)') i, a(i), a_Tsushima(i), a(i)/a_Tsushima(i)
        end if
      end do
    case (4)
      do i=1,Nch
        if (a(i)<0.) then
          a(i) = a_HADES(i)
        else
          write(*,'(i3,3f8.3)') i, a(i), a_HADES(i), a(i)/a_HADES(i)
        end if
      end do
    case default
      write(*,*) "BB_BYK bad parameter set!"
      stop
    end select
    call Write_ReadingInput('BB_BYK',1)

    first = .false.

  end subroutine readInput


   !***************************************************************************
   !****f* barBar_BarHypKaon/barBar_barBarMeson_strange
   ! NAME
   ! function barBar_barBarMeson_strange (srts, teilchenIN) result (sigma_BYK)
   ! PURPOSE
   ! XSections for strangeness production from BB->BYK collisions.
   ! INPUTS
   ! * real :: srts -- sqrt(s) of the incoming channel
   ! * type(particle),dimension(1:2) :: teilchenIN -- the incoming particles
   ! OUTPUT
   ! * real, dimension(1:18) :: sigma_BYK -- cross section for each channel (mb)
   ! NOTES
   ! * This routine calculates also the properties of the finalstates,
   !   needed for the construction of the Type(PreEvent) variable.
   ! * B stands for nucleon or a \Delta or a N* resonance. Other resonances
   !   not included. N* resonances are trated as nucleons.
   ! * In RMF mode sqrt(s*) already contains the self energy of the
   !   incoming channel. See its definition in generateFinalState in
   !   master_2Body.f90-module variable "srtS_XS".
   !***************************************************************************
   function barBar_barBarMeson_strange (srts, partIn) result (sigma_BYK)
     use particleDefinition
     use IdTable, only: Kaon, P11_1440
     ! Input-Output variables
     real,                           intent(in) :: srts
     type(particle), dimension(1:2), intent(in) :: partIn
     real,           dimension(1:18)            :: sigma_BYK
     ! Local variables
     real,    dimension(1:10) :: fac,fac0
     integer, dimension(1:10) :: iexit,isohf,isobf,iexit0,isohf0,isobf0
     integer :: iso1,iso2,nexit,nexit0,i,j
     real    :: srts2
     !--------------------------------------------------------------------------
     ! (1) Initialization
     !--------------------------------------------------------------------------
     nexit           = 0
     nexit0          = 0
     sigma_BYK(:)    = 0.    ! cross sections
     IdsOut(:,:)     = 0     ! IDs of 3-body final state
     ChargesOut(:,:) = 9999  ! charges of 3-body final state
     if (first) call readInput
     if (.not. enable) return
     !--------------------------------------------------------------------------
     ! (2) Only for Nucleons, Deltas, P11_1440's as incoming particles.
     !--------------------------------------------------------------------------
     if ((partIn(1)%ID > P11_1440).or.(partIn(2)%ID > P11_1440)) return
     !--------------------------------------------------------------------------
     ! (3) Define charge states of incoming particles.
     !     Convert variables related to charge states from GiBUU to RBUU scheme.
     !--------------------------------------------------------------------------
     iso1 = get_Charge (1, partIn(1)%ID, partIn(1)%charge)
     iso2 = get_Charge (1, partIn(2)%ID, partIn(2)%charge)
     !--------------------------------------------------------------------------
     ! (4) get number of final states & their properties
     !--------------------------------------------------------------------------
     call get_ChannelParameters (iso1, iso2, &
                                 nexit, iexit, fac, isohf, isobf, &
                                 nexit0, iexit0, fac0, isohf0, isobf0)
     !--------------------------------------------------------------------------
     ! (5) Calculate the cross section (in mb) for each final state.
     !     There are 30 primary channels: 1-28 are the original ones of
     !     Tsushima et al, 29-30 are NN->DeltaY*K channels added by K. Lapidus:
     !      * channel 29 = Delta-  Lambda*(1405) K+   (in analogy to channel 8)
     !      * channel 30 = Delta++ Sigma-*(1385) K+   (in analogy to channel 9)
     !     The remaining channels are determined through the primary ones using
     !     isospin-factors.
     !--------------------------------------------------------------------------

     ! the parametrizations are only valid up to sqrt(s) = 3.6 GeV
     srts2 = min(srts,3.6)**2

     ! BB-->BYK^{+}
     do i=1,nexit
        j=iexit(i)
        if (j==0) cycle  ! not yet determined channels
        sigma_BYK(i)    = fac(i) * xs_param(srts2,j)  ! Isospin factor*parametrization (mb)
        IdsOut(i,1)     = get_ID(1,isobf(i))
        IdsOut(i,2)     = get_ID(2,isohf(i))
        IdsOut(i,3)     = Kaon
        ChargesOut(i,1) = get_Charge(2,IdsOut(i,1),isobf(i))
        ChargesOut(i,2) = get_Charge(2,IdsOut(i,2),isohf(i))
        ChargesOut(i,3) = get_Charge(2,IdsOut(i,3),1)
     end do

     ! BB-->BYK^{0}
     do i=1,nexit0
        j=iexit0(i)
        if (j==0) cycle  ! not yet determined channels
        sigma_BYK(nexit+i)    = fac0(i) * xs_param(srts2,j)  ! Isospin factor*parametrization (mb)
        IdsOut(nexit+i,1)     = get_ID(1,isobf0(i))
        if (j>28) then
          IdsOut(nexit+i,2)   = get_ID(3,isohf0(i))
        else
        IdsOut(nexit+i,2)     = get_ID(2,isohf0(i))
        end if
        IdsOut(nexit+i,3)     = Kaon
        ChargesOut(nexit+i,1) = get_Charge(2,IdsOut(nexit+i,1),isobf0(i))
        ChargesOut(nexit+i,2) = get_Charge(2,IdsOut(nexit+i,2),isohf0(i))
        ChargesOut(nexit+i,3) = get_Charge(2,IdsOut(nexit+i,3),0)
     end do

   end function barBar_barBarMeson_strange



    !**************************************************************************
    !****f* barBar_BarHypKaon/xs_param
    ! NAME
    ! real function xs_param (s, j)
    ! PURPOSE
    ! Calculate the parametrized cross section for the primary channels
    ! (without isospin factor).
    ! INPUTS
    ! * real, intent(in) :: s      --- Mandelstam s in GeV**2
    ! * integer, intent(in) :: j   --- channel number
    ! OUTPUT
    ! * Cross section in mb.
    !**************************************************************************
    real function xs_param (s, j)
      real, intent(in) :: s
      integer, intent(in) :: j

      real :: x

      ! RMF: the incoming self energy is already included in srts!
      x = s/stresh(j)
      if (x<1.) then
        xs_param = 0.
      else
        xs_param = a(j) * (x-1.0)**b(j) / x**c(j)    ! in mb
      end if

    end function xs_param


    !**************************************************************************
    !****f* barBar_BarHypKaon/get_Channels_BYK
    ! NAME
    ! function get_Channels_BYK (InChan) result (parts)
    ! PURPOSE
    ! Defines the elements of the type(PreEvent) function
    ! INPUTS
    ! integer, intent(in) :: InChan -- Number of final channel
    ! OUTPUT
    ! type(preEvent), dimension(1:3) :: parts
    !**************************************************************************
    function get_Channels_BYK (InChan) result (parts)

      integer, intent(in) :: InChan
      type(preEvent), dimension(1:3) :: parts

      parts(1:3)%ID     = IdsOut (InChan,1:3)
      parts(1:3)%charge = ChargesOut (InChan,1:3)

    end function get_Channels_BYK


    !**************************************************************************
    !****f* barBar_BarHypKaon/Get_Charge
    ! NAME
    ! function Get_Charge
    ! PURPOSE
    ! Convert the variables (between GiBUU & Munich-RBUU) related to the isospin states.
    ! INPUTS
    ! integer :: direction : 1--GiBUU-->RBUU / 2--RBUU-->GiBUU.
    ! integer :: species : particle ID.
    ! integer :: iso_in  : charge of incoming particle.
    ! NOTES
    ! Only for nucleons, Delta's, N*'s, hyperons (Lambda,Sigma) and kaons (K^{+,0}).
    !**************************************************************************
    Function Get_Charge(direction,species,iso_in) result (iso_out)

      use IdTable, only: Nucleon,Delta,P11_1440,Lambda,SigmaResonance,Kaon,Lambda_1405,Sigma_1385

      integer, intent(in)  :: species   ! particle ID
      integer, intent(in)  :: iso_in    ! particle charge (input)
      integer, intent(in)  :: direction ! 1=GiBUU-->RBUU / 2=RBUU-->GiBUU

      integer :: conv,iso_out

      if (direction==1) then
         conv = 1
      else if (direction==2) then
         conv = -1
      else
         write(*,*) 'Module barBar_BarHypKaon / routine Get_Charge:'
         write(*,*) 'wrong input for direction: Direction = ',Direction
         write(*,*) 'STOP'
         STOP
      end if

      select case (species)
      case (nucleon,Lambda,Kaon,lambda_1405)
         iso_out = iso_in
      case (delta)
         iso_out = iso_in + conv*11
      case (P11_1440)
         iso_out = iso_in + conv*20
      case (SigmaResonance, Sigma_1385)
         iso_out = iso_in + conv*2
      case default
         write(*,*) 'Module barBar_BarHypKaon / routine Get_Charge:'
         write(*,*) 'Invalid input for particle ID: species = ',species
         write(*,*) 'STOP'
         STOP
      end select

    end Function Get_Charge


    !**************************************************************************
    !****f* barBar_BarHypKaon/Get_ID
    ! NAME
    ! Integer Function Get_ID (species, iso)
    ! PURPOSE
    ! Convert identification number (ID) between GiBUU and Munich-RBUU.
    ! INPUTS
    ! integer :: species : particle ID.
    ! integer :: iso_in  : charge of incoming particle.
    ! NOTES
    ! Only for nucleons, Delta's, N*'s, hyperons (Lambda,Sigma) and kaons (K^{+,0}).
    !**************************************************************************
    Integer Function Get_ID (species, iso)

      use IdTable, only: Nucleon,Delta,Lambda,SigmaResonance,Lambda_1405,Sigma_1385

      integer, intent(in) :: species,iso

      select case (species)
      case (1) ! Nucleon, Delta
         if (iso<10) then
            Get_ID = nucleon
         else
            Get_ID = Delta
         end if
      case (2) ! Hyperons: Lambda, Sigma
         if (iso==0) then
            Get_ID = Lambda
         else
            Get_ID = SigmaResonance
         end if
      case (3) ! Hyperon resonances: Lambda*(1405), Sigma*(1385)
         if (iso==0) then
            Get_ID = lambda_1405
         else
            Get_ID = Sigma_1385
         end if
      case default
         write(*,*) 'Module barBar_BarHypKaon / function get_ID:'
         write(*,*) 'wrong input for particle ID: species = ',species
         write(*,*) 'STOP'
         STOP
      end select

    end Function Get_ID



    !**************************************************************************
    !****s* barBar_BarHypKaon/get_ChannelParameters
    ! NAME
    ! subroutine get_ChannelParameters
    ! PURPOSE
    ! Returns the number of final states & their properties.
    ! INPUTS
    ! * integer :: iso1 -- charge of incoming particle Nr 1
    ! * integer :: iso2 -- charge of incoming particle Nr 2
    ! * integer :: isot -- total charge of incoming channel
    ! OUTPUT
    ! * integer :: nexit  -- Number of final states for BB-->BYK^{+}
    ! * integer :: iexit  -- Relation to table (see above)
    ! * integer :: fac    -- isospin factor for each final state in BYK^{+}
    ! * integer :: isohf  -- isospin-state of hyperon Y in BYK^{+}
    ! * integer :: isobf  -- isospin-state of the baryon B in BYK^{+}
    ! * integer :: nexit0 -- Number of final states for BB-->BYK^{0}
    ! * integer :: iexit0 -- Relation to table (see above)
    ! * integer :: fac0   -- isospin factor for each final state in BYK^{0}
    ! * integer :: isohf0 -- isospin-state of hyperon Y in BYK^{0}
    ! * integer :: isobf0 -- isospin-state of the baryon B in BYK^{0}
    ! NOTES
    ! * Notation:
    !   B : p,n,\DELTA^{-},\DELTA^{0},\DELTA^{+},\DELTA^{++}
    !   Y : \Lambda^{0},\Sigma^{-},\Sigma^{0},\Sigma^{+}
    !   K : Kaon^{+},Kaon^{0}
    ! * This is a partial copy of the Munich-RBUU code. For this reason we
    !   use here - and only here - local variables for the different isospin-states
    !   of the considered particle species, which may differ from those used
    !   in the GiBUU code. Please, DO NOT USE THE FOLLOWING PRESCRIPTIONS OUTSIDE
    !   OF THIS ROUTINE!!!
    ! * Relation between local variables and baryon/meson isospin states,
    !   used in this routine:
    ! * PROTONS     = 1
    ! * NEUTRONS    = 0
    ! * \DELTA^{-}  = 10
    ! * \DELTA^{0}  = 11
    ! * \DELTA^{+}  = 12
    ! * \DELTA^{++} = 13
    ! * N*^{0}      = 20
    ! * N*^{+}      = 21
    ! * \Lambda^{0} = 0
    ! * \Sigma^{-}  = 1
    ! * \Sigma^{0}  = 2
    ! * \Sigma^{+}  = 3
    ! * Kaon^{+}    = 1
    ! * Kaon^{0}    = 0
    !**************************************************************************
    subroutine get_ChannelParameters (iso1,iso2, &                         !<-- : Input
                                      nexit,iexit,fac,isohf,isobf, &       !<-- : Output (BB->BYK^{+})
                                      nexit0,iexit0,fac0,isohf0,isobf0)    !<-- : Output (BB->BYK^{0})

      ! Input-Output variables
      integer,                 intent(in)  :: iso1,iso2
      integer,                  intent(out) :: nexit,nexit0        ! number of K+ and K0 channels (each never larger than 10)
      integer, dimension(1:10), intent(out) :: iexit,isohf,isobf,iexit0,isohf0,isobf0
      real,    dimension(1:10), intent(out) :: fac,fac0
      ! Local variables
      integer :: isob1,isob2,isot
      ! treat N* like nucleons
      isob1 = mod(iso1,20)
      isob2 = mod(iso2,20)
      isot = isob1 + isob2  ! total charge of incoming particles
      !----------------------------------------------------------------------------
      ! nn --> BYK^{+} & BYK^{0} channels
      !----------------------------------------------------------------------------
      if (isot==0) then      ! nn
         nexit      = 7
         iexit(1:7) = (/  4,  8,   9,     9, 29,  30,    30 /)
         fac  (1:7) = (/ 1., 1., 0.5, 1./3., 1., 0.5, 1./3. /)
         isohf(1:7) = (/  1,  0,   2,     1,  0,   2,     1 /)
         isobf(1:7) = (/  0, 10,  10,    11, 10,  10,    11 /)

         nexit0      = 9
         iexit0(1:9) = (/  1,  3,  7,     8,  9,     9,    29, 30,    30 /)
         fac0  (1:9) = (/ 1., 1., 1., 1./3., 1., 1./6., 1./3., 1., 1./6. /)
         isohf0(1:9) = (/  0,  2,  1,     0,  3,     2,     0,  3,     2 /)
         isobf0(1:9) = (/  0,  0,  1,    11, 10,    11,    11, 10,    11 /)
         return
      end if
      !----------------------------------------------------------------------------
      ! pn --> BYK^{+} & BYK^{0} channels
      !----------------------------------------------------------------------------
      if (isot==1) then      ! pn
         nexit      = 9
         iexit(1:9) = (/  2,  5,  6,     8,     9,     9,    29,    30,    30 /)
         fac  (1:9) = (/ 1., 1., 1., 1./3., 1./6., 1./3., 1./3., 1./6., 1./3. /)
         isohf(1:9) = (/  0,  2,  1,     0,     2,     1,     0,     2,     1 /)
         isobf(1:9) = (/  0,  0,  1,    11,    11,    12,    11,    11,    12 /)

         nexit0      = 9
         iexit0(1:9) = (/  2,  5,  6,     8,     9,     9,    29,    30,    30 /)
         fac0  (1:9) = (/ 1., 1., 1., 1./3., 1./3., 1./6., 1./3., 1./3., 1./6. /)
         isohf0(1:9) = (/  0,  2,  3,     0,     3,     2,     0,     3,     2 /)
         isobf0(1:9) = (/  1,  1,  0,    12,    11,    12,    12,    11,    12 /)
         return
      end if
      !----------------------------------------------------------------------------
      ! pp --> BYK^{+} & BYK^{0} channels
      !----------------------------------------------------------------------------
      if (isot==2) then      ! pp
         nexit      = 10
         iexit(1:10) = (/  1,  3,  7,  0,     8,     9,  9,    29,    30, 30 /)
         fac  (1:10) = (/ 1., 1., 1., 0., 1./3., 1./6., 1., 1./3., 1./6., 1. /)
         isohf(1:10) = (/  0,  2,  3,  0,     0,     2,  1,     0,     2,  1 /)
         isobf(1:10) = (/  1,  1,  0,  0,    12,    12, 13,    12,    12, 13 /)

         nexit0      = 7
         iexit0(1:7) = (/  4,  8,   9,     9, 29,  30,    30 /)
         fac0  (1:7) = (/ 1., 1., 0.5, 1./3., 1., 0.5, 1./3. /)
         isohf0(1:7) = (/  3,  0,   2,     3,  0,   2,     3 /)
         isobf0(1:7) = (/  1, 13,  13,    12, 13,  13,    12 /)
         return
      end if
      !----------------------------------------------------------------------------
      ! n\Delta^{-} --> BYK^{+} & BYK^{0} channels
      !----------------------------------------------------------------------------
      if (isot==10) then     ! nD-
         nexit    = 1
         iexit(1) = 15
         fac  (1) = 1.
         isohf(1) = 1
         isobf(1) = 10

         nexit0      = 2
         iexit0(1:2) = (/ 12, 18 /)
         fac0  (1:2) = (/ 1., 1. /)
         isohf0(1:2) = (/  0,  2 /)
         isobf0(1:2) = (/ 10, 10 /)
         return
      end if
      !----------------------------------------------------------------------------
      ! n\Delta^{0},p\Delta^{-} --> BYK^{+} & BYK^{0} channels
      !----------------------------------------------------------------------------
      if (isot==11) then
         if (isob1==0 .or. isob2==0) then    ! nD0
            nexit      = 4
            iexit(1:4) = (/    11, 17,   13,   19 /)
            fac  (1:4) = (/ 1./3., 1., 0.75, 0.75 /)
            isohf(1:4) = (/     1,  1,    0,    2 /)
            isobf(1:4) = (/     0, 11,   10,   10 /)

            nexit0      = 6
            iexit0(1:6) = (/    10,    11, 14,   16, 20, 21 /)
            fac0  (1:6) = (/ 1./3., 1./6., 1., 0.75, 1., 1. /)
            isohf0(1:6) = (/     0,     2,  0,    3,  2,  1 /)
            isobf0(1:6) = (/     0,     0, 11,   10, 11, 12 /)
         else                                ! pD-
            nexit      = 4
            iexit(1:4) = (/ 11,   16, 12, 18 /)
            fac  (1:4) = (/ 1., 0.75, 1., 1. /)
            isohf(1:4) = (/  1,    1,  0,  2 /)
            isobf(1:4) = (/  0,   11, 10, 10 /)

            nexit0      = 6
            iexit0(1:6) = (/ 10,  11,   13, 15,   19,   21 /)
            fac0  (1:6) = (/ 1., 0.5, 0.75, 1., 0.75, 0.75 /)
            isohf0(1:6) = (/  0,   2,    0,  3,    2,    2 /)
            isobf0(1:6) = (/  0,   0,   11, 10,   11,   11 /)
         end if
         return
      end if
      !----------------------------------------------------------------------------
      ! n\Delta^{+},p\Delta^{0} --> BYK^{+} & BYK^{0} channels
      !----------------------------------------------------------------------------
      if (isot==12) then
         if (isob1==0 .or. isob2==0) then    ! nD+
            nexit      = 7
            iexit(1:7) = (/    10,    11, 19, 13,  0,    11, 17 /)
            fac  (1:7) = (/ 1./3., 1./6., 1., 1., 0., 1./3., 1. /)
            isohf(1:7) = (/     0,     2,  2,  0,  0,     1,  1 /)
            isobf(1:7) = (/     0,     0, 11, 11,  0,     1, 12 /)

            nexit0      = 7
            iexit0(1:7) = (/    10,    11,    11, 14, 16, 20,   21 /)
            fac0  (1:7) = (/ 1./3., 1./3., 1./6., 1., 1., 1., 0.75 /)
            isohf0(1:7) = (/     0,     3,     2,  0,  3,  2,    1 /)
            isobf0(1:7) = (/     1,     0,     1, 12, 11, 12,   13 /)
         else                                ! pD0
            nexit      = 7
            iexit(1:7) = (/    10,    11, 20, 14,   21,    11, 16 /)
            fac  (1:7) = (/ 1./3., 1./6., 1., 1., 0.75, 1./3., 1. /)
            isohf(1:7) = (/     0,     2,  2,  0,    3,     1,  1 /)
            isobf(1:7) = (/     0,     0, 11, 11,   10,     1, 12 /)

            nexit0      = 6
            iexit0(1:6) = (/    10,    11,    11, 13, 17, 19 /)
            fac0  (1:6) = (/ 1./3., 1./3., 1./6., 1., 1., 1. /)
            isohf0(1:6) = (/     0,     3,     2,  0,  3,  2 /)
            isobf0(1:6) = (/     1,     0,     1, 12, 11, 12 /)
         end if
         return
      end if
      !----------------------------------------------------------------------------
      ! n\Delta^{++},p\Delta^{+} --> BYK^{+} & BYK^{0} channels
      !----------------------------------------------------------------------------
      if (isot==13) then
         if (isob1==0 .or. isob2==0) then    ! nD++
            nexit      = 7
            iexit(1:7) = (/ 10,  11,  0,  0,   13,   19, 15 /)
            fac  (1:7) = (/ 1., 0.5, 0., 0., 0.75, 0.75, 1. /)
            isohf(1:7) = (/  0,   2,  0,  0,    0,    2,  1 /)
            isobf(1:7) = (/  1,   1,  0,  0,   12,   12, 13 /)

            nexit0      = 4
            iexit0(1:4) = (/ 11, 12,   16, 18 /)
            fac0  (1:4) = (/ 1., 1., 0.75, 1. /)
            isohf0(1:4) = (/  3,  0,    3,  2 /)
            isobf0(1:4) = (/  1, 13,   12, 13 /)
         else                                ! pD+
            nexit      = 7
            iexit(1:7) = (/    10,    11,  0, 21, 14, 20,   16 /)
            fac  (1:7) = (/ 1./3., 1./6., 0., 1., 1., 1., 0.75 /)
            isohf(1:7) = (/     0,     2,  0,  3,  0,  2,    1 /)
            isobf(1:7) = (/     1,     1,  0, 11, 12, 12,   13 /)

            nexit0      = 4
            iexit0(1:4) = (/    11,   13, 17,   19 /)
            fac0  (1:4) = (/ 1./3., 0.75, 1., 0.75 /)
            isohf0(1:4) = (/     3,    0,  3,    2 /)
            isobf0(1:4) = (/     1,   13, 12,   13 /)
         end if
         return
      end if
      !----------------------------------------------------------------------------
      ! p\Delta^{++} --> BYK^{+} & BYK^{0} channels
      !----------------------------------------------------------------------------
      if (isot==14) then                ! pD++
         nexit      = 4
         iexit(1:4) = (/  0,   21, 12, 18 /)
         fac  (1:4) = (/ 0., 0.75, 1., 1. /)
         isohf(1:4) = (/  0,    3,  0,  2 /)
         isobf(1:4) = (/  0,    1, 13, 13 /)

         nexit0    = 1
         iexit0(1) = 15
         fac0  (1) = 1.
         isohf0(1) = 3
         isobf0(1) = 13
         return
      end if
      !----------------------------------------------------------------------------
      ! \Delta^{-}\Delta^{-} --> no final channel
      !----------------------------------------------------------------------------
      if (isot==20) then                ! D-D-
         nexit  = 0
         nexit0 = 0
         return
      end if
      !----------------------------------------------------------------------------
      ! \Delta^{-}\Delta^{0} --> BYK^{+} & BYK^{0} channels
      !----------------------------------------------------------------------------
      if (isot==21) then                ! D-D0
         nexit    = 1
         iexit(1) = 26
         fac  (1) = 1.0
         isohf(1) = 1
         isobf(1) = 10

         nexit0      = 2
         iexit0(1:2) = (/ 22,  26 /)
         fac0  (1:2) = (/ 1., 0.5 /)
         isohf0(1:2) = (/  0,   2 /)
         isobf0(1:2) = (/ 10,  10 /)
         return
      end if
      !----------------------------------------------------------------------------
      ! \Delta^{-}\Delta^{+},\Delta^{0}\Delta^{0} --> BYK^{+} & BYK^{0} channels
      !----------------------------------------------------------------------------
      if (isot==22) then
         if (isob1==10 .or. isob2==10) then   ! D-D+
            nexit      = 3
            iexit(1:3) = (/ 28,  0, 25 /)
            fac  (1:3) = (/ 1., 0., 1. /)
            isohf(1:3) = (/  1,  0,  2 /)
            isobf(1:3) = (/ 11,  0, 10 /)

            nexit0      = 3
            iexit0(1:3) = (/ 23, 25, 27 /)
            fac0  (1:3) = (/ 1., 1., 1. /)
            isohf0(1:3) = (/  0,  3,  2 /)
            isobf0(1:3) = (/ 11, 10, 11 /)
         else                                 ! D0D0
            nexit      = 3
            iexit(1:3) = (/    26,    22,    26 /)
            fac  (1:3) = (/ 1./9., 1./3., 1./6. /)
            isohf(1:3) = (/     1,     0,     2 /)
            isobf(1:3) = (/    11,    10,    10 /)

            nexit0      = 3
            iexit0(1:3) = (/    22,    26,     26 /)
            fac0  (1:3) = (/ 1./9., 1./3., 1./18. /)
            isohf0(1:3) = (/     0,     3,      2 /)
            isobf0(1:3) = (/    11,    10,     11 /)
         end if
         return
      end if
      !----------------------------------------------------------------------------
      ! \Delta^{-}\Delta^{++},\Delta^{0}\Delta^{+} --> BYK^{+} & BYK^{0} channels
      !----------------------------------------------------------------------------
      if (isot==23) then
         if (isob1==10 .or. isob2==10) then   ! D-D++
            nexit      = 3
            iexit(1:3) = (/  0, 25, 25 /)
            fac(1:3)   = (/ 1., 1., 1. /)
            isohf(1:3) = (/  0,  1,  2 /)
            isobf(1:3) = (/  0, 12, 11 /)

            nexit0      = 2
            iexit0(1:2) = (/ 25, 25 /)
            fac0  (1:2) = (/ 1., 1. /)
            isohf0(1:2) = (/  2,  3 /)
            isobf0(1:2) = (/ 12, 11 /)
         else                                 ! D0D+
            nexit      = 4
            iexit(1:4) = (/    28, 24,  0,    27 /)
            fac  (1:4) = (/ 1./6., 1., 0., 2./3. /)
            isohf(1:4) = (/     2,  0,  0,     1 /)
            isobf(1:4) = (/    11, 11,  0,    12 /)

            nexit0      = 3
            iexit0(1:3) = (/ 24,    27,    28 /)
            fac0  (1:3) = (/ 1., 2./3., 1./6. /)
            isohf0(1:3) = (/  0,     3,     2 /)
            isobf0(1:3) = (/ 12,    11,    12 /)
         end if
         return
      end if
      !----------------------------------------------------------------------------
      ! \Delta^{0}\Delta^{++},\Delta^{+}\Delta^{+} --> BYK^{+} & BYK^{0} channels
      !----------------------------------------------------------------------------
      if (isot==24) then
         if (isob1==11 .or. isob2==11) then   ! D0D++
            nexit      = 4
            iexit(1:4) = (/  0, 23, 27, 25 /)
            fac  (1:4) = (/ 0., 1., 1., 1. /)
            isohf(1:4) = (/  0,  0,  2,  1 /)
            isobf(1:4) = (/  0, 12, 12, 13 /)

            nexit0      = 2
            iexit0(1:2) = (/ 25, 28 /)
            fac0  (1:2) = (/ 1., 1. /)
            isohf0(1:2) = (/  2,  3 /)
            isobf0(1:2) = (/ 13, 12 /)
         else                                 ! D+D+
            nexit      = 4
            iexit(1:4) = (/  0,    22,     26,    26 /)
            fac  (1:4) = (/ 0., 1./9., 1./18., 1./3. /)
            isohf(1:4) = (/  0,     0,      2,     1 /)
            isobf(1:4) = (/  0,    12,     12,    13 /)

            nexit0      = 3
            iexit0(1:3) = (/    22,    26,    26 /)
            fac0  (1:3) = (/ 1./3., 1./6., 1./9. /)
            isohf0(1:3) = (/     0,     2,     3 /)
            isobf0(1:3) = (/    13,    13,    12 /)
         end if
         return
      end if
      !----------------------------------------------------------------------------
      ! \Delta^{+}\Delta^{++} --> BYK^{+} & BYK^{0} channels
      !----------------------------------------------------------------------------
      if (isot==25) then     ! D+D++
         nexit      = 3
         iexit(1:3) = (/  0, 22,  26 /)
         fac  (1:3) = (/ 0., 1., 0.5 /)
         isohf(1:3) = (/  0,  0,   2 /)
         isobf(1:3) = (/  0, 13,  13 /)

         nexit0    = 1
         iexit0(1) = 26
         fac0  (1) = 1.
         isohf0(1) = 3
         isobf0(1) = 13
         return
      end if
      !----------------------------------------------------------------------------
      ! \Delta^{++}\Delta^{++} --> no valid final channel
      !----------------------------------------------------------------------------
      if (isot==26) then                ! D++D++
         nexit  = 0
         nexit0 = 0
         return
      end if
      !----------------------------------------------------------------------------
      ! Final check: if everything is fine, the simulation should not arrive here...
      !----------------------------------------------------------------------------
      write(*,*) 'No final channel found??? something is going wrong...'
      write(*,*) 'Charges of initial channel: iso1,iso2,iso1,iso2,isot = ', &
           & iso1,iso2,isob1,isob2,isot
      call Traceback()

    end subroutine get_ChannelParameters



end module barBar_BarHypKaon
