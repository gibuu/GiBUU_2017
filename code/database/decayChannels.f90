!******************************************************************************
!****m* /decayChannels
! NAME
! Module decayChannels
! PURPOSE
! Defines possible Decay channels for for baryons and mesons.
! After initialization all the information can be found in the
! file "./DecayChannels.tex".
!******************************************************************************
Module DecayChannels
  implicit none
  private


  !****************************************************************************
  !****t*  decayChannels/tDecay3body
  ! SOURCE
  !
  type, public :: tDecay3body
     integer, dimension(1:3) :: id = 0          ! Id's of the final state particles
     integer, dimension(1:3) :: charge = 0      ! charges of the final state particles
     real                    :: threshold = 0.  ! decay threshold in GeV
     logical, dimension(1:3) :: isAnti = .false.! flag antiparticle
  end type tDecay3body
  !****************************************************************************


  !****************************************************************************
  !****t*  decayChannels/tDecay2body
  ! SOURCE
  !
  type, public :: tDecay2body
     integer,dimension(1:2) :: id = 0              ! Id's of final state particles
     logical,dimension(1:2) :: stable = .true.     ! whether final state particle is stable
     integer                :: angularMomentum = 0 ! Angular momentum of final state
     real                   :: threshold = 0.      ! decay threshold in GeV
  end type tDecay2body
  !****************************************************************************


  integer, parameter, public :: nDecay2bodyMeson  = 18
  integer, parameter, public :: nDecay2bodyBaryon = 45
  integer, parameter, public :: nDecay3bodyMeson  =  4
  integer, parameter, public :: nDecay3bodyBaryon =  0


  !****************************************************************************
  !****g* decayChannels/Decay2bodyMeson
  ! SOURCE
  !
  type(tDecay2body), dimension(0:nDecay2bodyMeson), save, public :: Decay2bodyMeson
  ! PURPOSE
  ! Stores the information for all possible 2-body decay channels of all mesons.
  !****************************************************************************


  !****************************************************************************
  !****g* decayChannels/Decay2bodyBaryon
  ! SOURCE
  !
  type(tDecay2body), dimension(0:nDecay2bodyBaryon), save, public :: Decay2bodyBaryon
  ! PURPOSE
  ! Stores the information for all possible 2-body decay channels of all baryons.
  !****************************************************************************


  !****************************************************************************
  !****g* decayChannels/Decay3bodyMeson
  ! SOURCE
  !
  type(tDecay3body), dimension(0:nDecay3bodyMeson), save, public :: Decay3bodyMeson
  ! PURPOSE
  ! Stores the information for all possible 3-body decay channels of all mesons.
  !****************************************************************************


  !****************************************************************************
  !****g* decayChannels/Decay3bodyBaryon
  ! SOURCE
  !
  type(tDecay3body), dimension(0:nDecay3bodyBaryon), save, public :: Decay3bodyBaryon
  ! PURPOSE
  ! Stores the information for all possible 3-body decay channels of all baryons.
  !****************************************************************************


  !****************************************************************************
  !****g* decayChannels/rhoDelta_is_sigmaDelta
  ! SOURCE
  !
  logical, save :: rhoDelta_is_sigmaDelta = .false.
  ! PURPOSE
  ! If true, the rho-Delta decay channel will be replaced by sigma-Delta.
  ! For discussion, see e.g. Effenberger PhD, chapter 6.3.2.
  !****************************************************************************


  public :: InitDecayChannels, Print_DecayChannels, get_rhoDelta_is_sigmaDelta


  logical, save :: first = .true.


contains


  logical function get_rhoDelta_is_sigmaDelta ()
    if (first) call readInput
    get_rhoDelta_is_sigmaDelta = rhoDelta_is_sigmaDelta
  end function


  !****************************************************************************
  !****s*  decayChannels/InitDecayChannels
  ! NAME
  ! subroutine InitDecayChannels
  ! PURPOSE
  ! Define the decay channels for Baryons and Mesons which are used in BUU.
  !****************************************************************************
  subroutine InitDecayChannels
    use idtable

    if (first) call readInput

    ! define the decay channels ...
    call baryon_twoBody
    call meson_twoBody
    call meson_threeBody

  contains

    !**************************************************************************
    !****s* InitDecayChannels/baryon_TwoBody
    ! NAME
    ! subroutine baryon_TwoBody
    ! PURPOSE
    ! Defines possible 2-body decay channels for for the baryons.
    ! We set also flags, which decide whether the final
    ! state particles shall be treated as a stable particles.
    ! And the angular momentum of the decay is defined.
    !**************************************************************************
    subroutine baryon_TwoBody

      ! with s=0
      ! stable final states:
      Decay2bodyBaryon(1)%ID=(/pion,nucleon/)
      Decay2bodyBaryon(2)%ID=(/eta,nucleon/)
      Decay2bodyBaryon(3)%ID=(/omegaMeson,nucleon/)
      Decay2bodyBaryon(4)%ID=(/kaon,lambda/)

      ! pion Delta channels :
      Decay2bodyBaryon(5:8)%ID(1)=pion
      Decay2bodyBaryon(5:8)%ID(2)=delta
      Decay2bodyBaryon(5:8)%stable(2)=.false.
      Decay2bodyBaryon(5:8)%AngularMomentum=(/0,1,2,3/)

      ! rho Nucleon channels :
      Decay2bodyBaryon(9:12)%ID(1)=rho
      Decay2bodyBaryon(9:12)%ID(2)=nucleon
      Decay2bodyBaryon(9:12)%stable(1)=.false.
      Decay2bodyBaryon(9:12)%AngularMomentum=(/0,1,2,3/)

      Decay2bodyBaryon(13)%ID=(/sigmaMeson,nucleon/)
      Decay2bodyBaryon(13)%stable(1)=.false.

      Decay2bodyBaryon(14)%ID=(/pion,P11_1440/)
      Decay2bodyBaryon(14)%stable(2)=.false.

      ! rho Delta channels
      if (rhoDelta_is_sigmaDelta) then
        Decay2bodyBaryon(15:18)%ID(1) = sigmaMeson
        Decay2bodyBaryon(15:18)%ID(2) = Delta
      else
        Decay2bodyBaryon(15:18)%ID(1) = rho
        Decay2bodyBaryon(15:18)%ID(2) = Delta
      end if
      Decay2bodyBaryon(15:18)%stable(1) = .false.
      Decay2bodyBaryon(15:18)%stable(2) = .false.
      Decay2bodyBaryon(15:18)%AngularMomentum = (/0,1,2,3/)

      ! with S=-1:

      Decay2bodyBaryon(19)%ID=(/pion,lambda/)
      Decay2bodyBaryon(20)%ID=(/kaonBar,nucleon/)
      Decay2bodyBaryon(21)%ID=(/pion,SigmaResonance/)

      Decay2bodyBaryon(22:26)%ID(1)=pion
      Decay2bodyBaryon(22:26)%ID(2)=Sigma_1385
      Decay2bodyBaryon(22:26)%stable(2)=.false.
      Decay2bodyBaryon(22:26)%AngularMomentum = (/0,1,2,3,4/)

      Decay2bodyBaryon(27)%ID=(/eta,lambda/)
      Decay2bodyBaryon(28)%ID=(/eta,SigmaResonance/)

      Decay2bodyBaryon(29)%ID=(/kaonBar,delta/)
      Decay2bodyBaryon(29)%stable(2)=.false.

      Decay2bodyBaryon(30:34)%ID(1)=kaonStarBar
      Decay2bodyBaryon(30:34)%ID(2)=nucleon
      Decay2bodyBaryon(30:34)%stable(1)=.false.
      Decay2bodyBaryon(30:34)%AngularMomentum = (/0,1,2,3,4/)

      Decay2bodyBaryon(35)%ID=(/pion,lambda_1405/)
      Decay2bodyBaryon(35)%stable(2)=.false.

      Decay2bodyBaryon(36:40)%ID(1)=pion
      Decay2bodyBaryon(36:40)%ID(2)=lambda_1520
      Decay2bodyBaryon(36:40)%stable(2)=.false.
      Decay2bodyBaryon(36:40)%AngularMomentum = (/0,1,2,3,4/)

      Decay2bodyBaryon(41)%ID=(/pion,lambda_1820/)
      Decay2bodyBaryon(41)%stable(2)=.false.

      Decay2bodyBaryon(42)%ID=(/omegaMeson,lambda/)

      ! with S=-2:
      Decay2bodyBaryon(43)%ID=(/pion,xi/)

      ! with C=1:
      Decay2bodyBaryon(44)%ID=(/pion,lambda_CPlus/)
      Decay2bodyBaryon(45)%ID=(/pion,Xi_C/)

    end subroutine baryon_TwoBody

    !**************************************************************************
    !****s* InitDecayChannels/meson_TwoBody
    ! NAME
    ! subroutine meson_TwoBody
    ! PURPOSE
    ! Defines possible 2-body decay channels for for the mesons.
    ! We set also flags, which decide whether the final
    ! state particles shall be treated as a stable particles.
    !**************************************************************************
    subroutine meson_TwoBody

      !Define decay channels
      Decay2bodyMeson( 1)%ID(1:2) = (/pion,pion/)

      Decay2bodyMeson( 2)%ID(1:2) = (/pion,rho/)
      Decay2bodyMeson( 2)%stable(2) = .false.

      Decay2bodyMeson( 3)%ID(1:2) = (/kaon,kaonBar/)
      Decay2bodyMeson( 4)%ID(1:2) = (/kaon,pion/)

      Decay2bodyMeson( 5)%ID(1:2) = (/rho,photon/)
      Decay2bodyMeson( 5)%stable(1) = .false.

      Decay2bodyMeson( 6)%ID(1:2) = (/pion,photon/)
      Decay2bodyMeson( 7)%ID(1:2) = (/photon,photon/)
      Decay2bodyMeson( 8)%ID(1:2) = (/photon,dS_Plus/)
      Decay2bodyMeson( 9)%ID(1:2) = (/photon,dS_Minus/)
      Decay2bodyMeson(10)%ID(1:2) = (/pion,dS_Plus/)
      Decay2bodyMeson(11)%ID(1:2) = (/pion,dS_Minus/)
      Decay2bodyMeson(12)%ID(1:2) = (/kaonBar,pion/)
      Decay2bodyMeson(13)%ID(1:2) = (/dMeson,photon/)
      Decay2bodyMeson(14)%ID(1:2) = (/dBar,photon/)
      Decay2bodyMeson(15)%ID(1:2) = (/pion,dMeson/)
      Decay2bodyMeson(16)%ID(1:2) = (/pion,dBar/)

      Decay2bodyMeson(17)%ID(1:2) = (/pion,sigmaMeson/)
      Decay2bodyMeson(17)%stable(2) = .false.

      Decay2bodyMeson(18)%ID(1:2) = (/eta,sigmaMeson/)
      Decay2bodyMeson(18)%stable(2) = .false.

    end subroutine meson_twoBody

    !**************************************************************************
    !****s* InitDecayChannels/meson_ThreeBody
    ! NAME
    ! subroutine meson_ThreeBody
    ! PURPOSE
    ! Defines possible 3-body decay channels for for the mesons.
    ! Here also the charges of the final state particles are defined
    !*********************************************************************^
    subroutine meson_ThreeBody

      Decay3bodyMeson(1)%Id(1:3)=(/pion,pion,eta/)
      Decay3bodyMeson(1)%charge(1:3)=0

      Decay3bodyMeson(2)%Id(1:3)=pion
      Decay3bodyMeson(2)%charge(1:3)=(/0,-1,1/)

      Decay3bodyMeson(3)%Id(1:3)=pion
      Decay3bodyMeson(3)%charge(1:3)=0

      Decay3bodyMeson(4)%Id(1:3)=(/pion,pion,eta/)
      Decay3bodyMeson(4)%charge(1:3)=(/1,-1,0/)

    end subroutine meson_threeBody


  end subroutine InitDecayChannels


  !****************************************************************************
  !****s* decayChannels/readInput
  ! NAME
  ! subroutine readInput
  ! PURPOSE
  ! Reads input in jobcard out of namelist "DecayChannels".
  !****************************************************************************
  subroutine readInput
    use output, only: Write_ReadingInput

    !**************************************************************************
    !****n*  decayChannels/DecayChannels
    ! NAME
    ! NAMELIST /DecayChannels/
    ! PURPOSE
    ! Includes the input switches:
    ! * rhoDelta_is_sigmaDelta
    !**************************************************************************
    NAMELIST /decayChannels/ rhoDelta_is_sigmaDelta

    integer :: ios

    call Write_ReadingInput('decayChannels',0)
    rewind(5)
    read(5,nml=decayChannels,IOSTAT=ios)
    call Write_ReadingInput('decayChannels',0,ios)

    write(*,*) 'rhoDelta_is_sigmaDelta : ', rhoDelta_is_sigmaDelta

    call Write_ReadingInput('decayChannels',1)

    first = .false.

  end subroutine readInput


  !****************************************************************************
  !****s* decayChannels/Print_DecayChannels
  ! NAME
  ! subroutine Print_DecayChannels
  ! PURPOSE
  ! Prints all information about decay channels to file 'GiBUU_database_decayChannels.txt'.
  !
  ! NOTES
  ! This is preliminary: Output has to be refined.
  !****************************************************************************
  subroutine Print_DecayChannels
    integer :: i

    open(10,file='GiBUU_database_decayChannels.txt',status='unknown')

    write(10,*) '2-body decay channels of the hadrons'
    write(10,'(A)') ' i & final State 1 & final State 2 & 1 is stable & 2 is stable & angular momentum & threshold'

    do i=1,nDecay2bodyMeson
      write(10,'(3(I3," & "),2(L6," & "),I3," & ",f5.3)') i, Decay2bodyMeson(i)%id, Decay2bodyMeson(i)%stable, &
                                                          Decay2bodyMeson(i)%angularMomentum, Decay2bodyMeson(i)%threshold
    end do

    write(10,*)

    do i=1,nDecay2bodyBaryon
      write(10,'(3(I3," & "),2(L6," & "),I3," & ",f5.3)') i, Decay2bodyBaryon(i)%id, Decay2bodyBaryon(i)%stable, &
                                                          Decay2bodyBaryon(i)%angularMomentum, Decay2bodyBaryon(i)%threshold
    end do

    write(10,*)
    write(10,*) '3-body decay channels of the hadrons'
    write(10,'(A)') ' i & final State 1 & final State 2 & final State 3 & charge 1 & charge 2 & charge 3 & threshold'

    do i=1,nDecay3bodyMeson
      write(10,'(6(I3," & "),I3," & ",f5.3)') i,Decay3bodyMeson(i)%id, Decay3bodyMeson(i)%charge, Decay3bodyMeson(i)%threshold
    end do

    write(10,*)

    do i=1,nDecay3bodyBaryon
      write(10,'(6(I3," & "),I3," & ",f5.3)') i,Decay3bodyBaryon(i)%id, Decay3bodyBaryon(i)%charge, Decay3bodyBaryon(i)%threshold
    end do

    close(10)

  end subroutine Print_DecayChannels


end module DecayChannels
