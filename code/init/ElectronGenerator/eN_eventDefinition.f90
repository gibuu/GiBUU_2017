!******************************************************************************
!****m* /eN_eventDefinition
! NAME
! module eN_eventDefinition
!
! PURPOSE
! This module includes a type definition which deals with kinematics of
! lepton-nucleon scattering.
!
! NOTES
! * historically, this type was developed to hold informations in electron
!   induced reactions. Since we enhanced the usage also to neutrino induced
!   reactions, it should be renamed into "type leptonNucleon_event".
!
!******************************************************************************
module eN_eventDefinition

  use particleDefinition
  implicit none
  private


  !****************************************************************************
  !****t* eN_eventDefinition/electronNucleon_event
  ! NAME
  ! type electronNucleon_event
  ! PURPOSE
  ! This is a type to define an electron-nucleon scattering event
  !
  ! SOURCE
  !
  Type electronNucleon_event
     type(particle)        :: lepton_in  ! incoming lepton (or neutrino)
     type(particle)        :: lepton_out ! outgoing lepton
     type(particle)        :: boson      ! exchanged boson

     type(particle)        :: nucleon       ! target nucleon
     type(particle)        :: nucleon_free  ! free target nucleon  (potentials removed)

     type(particle)        :: nucleon2      ! second target nucleon (2p2h processes)

     real                  :: QSquared = 0. ! Q^2=-q^2
     real                  :: W = 0.        ! W at the hadronic vertex
     real                  :: W_free = 0.   ! Free W   (potentials removed, momentum kept)
     real                  :: W_rec = 0.    ! Free W for nucleon at rest
     real, dimension(0:3)  :: pcm = 0.      ! Lorentz-Trafo into cm-frame
     real, dimension(3)    :: betacm = 0.   ! Lorentz-Trafo into cm-frame
     real                  :: phiLepton = 0.! additional angle of lepton plane

     integer               :: idProcess = 0 ! = (+-)EM,(+-)CC,(+-)NC
     integer               :: idFamily  = 0 ! Abbrev. for: e,mu,tau

  end Type electronNucleon_event
  !
  ! NOTES
  ! 1) The integers 'idProcess' and 'idFamily' keep redundant information,
  ! but are useful as abbreviations in if-statements. (cf.module leptonicID)
  !
  !****************************************************************************


  !****************************************************************************
  !****t* eN_eventDefinition/Constants_Frame
  ! SOURCE
  !
  integer, parameter :: doNOT=0
  integer, parameter :: CM=1
  integer, parameter :: CALC=2
  integer, parameter :: THRE=3
  integer, parameter :: NucleonRest=4
  integer, parameter :: THRE2=5
  !****************************************************************************


  public :: CM, CALC, doNOT, THRE, NucleonRest, THRE2
  public :: electronNucleon_event
  public :: write_electronNucleon_event
  public :: setVacuum

contains

  !****************************************************************************
  !****s* eN_eventDefinition/write_electronNucleon_event
  ! NAME
  ! subroutine write_electronNucleon_event(e,DoDebug,DoShort)
  !
  ! PURPOSE
  ! Dump the given instance to stdout
  !****************************************************************************
  subroutine write_electronNucleon_event(e,DoDebug_,DoShort_)
    use output, only: line, writeParticle_debug,writeParticle
    use minkowski, only: abs4
    use particleProperties, only: PartName

    type(electronNucleon_event), intent(in) :: e
    logical, intent(in), optional :: DoDebug_, DoShort_
    character(50) :: format4,format1  ! ,format3
    character*(*), dimension(3), parameter :: sProcess = (/"EM","CC","NC"/)
    character*(*), dimension(3), parameter :: sFamily  = (/"e  ","mu ","tau"/)

    logical :: DoDebug, DoShort

    DoShort  = .false.
    if (present(DoShort_)) DoShort=DoShort_

    DoDebug  = .false.
    if (present(DoDebug_)) DoDebug=DoDebug_


    if (DoDebug) then
       format4='(A,1P,4E15.4)'
!        format3='(A,1P,15(" "),3E15.4)'
    else
       format4='(A,1P,"                   ",4E11.3)'
!        format3='(A,1P,"                   ",11(" "),3E15.4)'
    end if

    format1='(A,1P,E15.4)'
!     format3='(A,1P,4E15.4)'



    write(*,*)
    write(*,'(A)') line
    write(*,'(A)') line
    write(*,'(A)') '**** Lepton-Nucleon Event : ****'
    write(*,'(A)') line
    write(*,'(A)') line
    select case (e%idProcess)
    case (1:3)
       write(*,'(A,A)') '...Process: ',sProcess(e%idProcess)
    case (-3:-1)
       write(*,'(A,A)') '...Process: anti-',sProcess(e%idProcess)
    case default
       write(*,'(A,A,i5)') '...Process: ','***unknown*** ',e%idProcess
    end select
    select case (e%idFamily)
    case (1:3)
       write(*,'(A,A)') '...Family:  ',sFamily(e%idFamily)
    case default
       write(*,'(A,A,i5)') '...Family:  ','***unknown*** ',e%idFamily
    end select

    write(*,'(A)') line
    if (DoDebug) then
       write(*,format4) '* Electron incoming     ', e%lepton_in%momentum
       write(*,format4) '* Electron outgoing     ', e%lepton_out%momentum
       write(*,format4) '* Photon momentum       ', e%boson%momentum
    else
       write(*,'(A,A)') '* Incoming Lepton: ',PartName(e%Lepton_in)
       call WriteParticle(6,1,1,e%lepton_in)
       write(*,'(A,A)') '* Outgoing Lepton: ',PartName(e%Lepton_out)
       call WriteParticle(6,1,1,e%lepton_out)
       write(*,'(A,A)') '* Exchanged Boson: ',PartName(e%boson)
       call WriteParticle(6,1,1,e%boson)
    end if
    write(*,*)
    write(*,'(A,F15.5)') '* Nucleon:           m_eff =',abs4(e%nucleon%momentum)
    if (DoDebug) then
       call writeParticle_Debug(e%nucleon)
    else
       call WriteParticle(6,1,1,e%nucleon)
    end if

    if (.not.DoShort) then
       write(*,*)
       write(*,'(A,F15.5)') '* Nucleon free:      m_eff =',abs4(e%nucleon_free%momentum)
       if (DoDebug) then
          call writeParticle_Debug(e%nucleon_free)
       else
          call WriteParticle(6,1,1,e%nucleon_free)
       end if
    end if

    if (e%nucleon2%ID.gt.0) then
       write(*,*)
       write(*,'(A,F15.5)') '* Nucleon 2:         m_eff =',abs4(e%nucleon2%momentum)
       if (DoDebug) then
          call writeParticle_Debug(e%nucleon2)
       else
          call WriteParticle(6,1,1,e%nucleon2)
       end if
    end if

    write(*,*)
    write(*,format1) '* QSquared    =', e%QSquared
    write(*,format1) '* W           =' ,e%W

    if (.not.DoShort) then
       write(*,format1) '* W_free      =' , e%W_free
       write(*,*)
       write(*,format4) '* Boost: pcm            ',e%pcm
       write(*,format4) '* Boost: betacm         ',e%betacm
       write(*,format1) '* Boost: phi_lepton     ',e%phiLepton
    end if

    write(*,'(A)') line
    write(*,'(A)') line
    write(*,*)

  end subroutine  write_electronNucleon_event

  !****************************************************************************
  !****f* eN_eventDefinition/setVacuum
  ! NAME
  ! function setVacuum(eN) result(eN_vac)
  !
  ! PURPOSE
  ! Return the electron-nucleon event transformed to vacuum kinematics
  !****************************************************************************
  function setVacuum(eN) result(eN_vac)
    type(electronNucleon_event),intent(in) :: en
    type(electronNucleon_event)            :: en_Vac

    eN_Vac=eN
    en_Vac%nucleon      =  en_Vac%nucleon_free
    en_Vac%W            =  en_Vac%W_free
    return
  end function setVacuum

end module eN_eventDefinition
