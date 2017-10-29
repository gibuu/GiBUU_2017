!******************************************************************************
!****m* /expNeutrinofluxes
! NAME
! module expNeutrinofluxes
!
! PURPOSE
! This module provides specific experimental neutrino fluxes and it selects the
! neutrino energy according to the experimental flux.
! For MiniBooNE and K2K, it also extracts the reconstructed neutrino energy and
! Qs as it is done in the experiment.
!******************************************************************************

module expNeutrinofluxes

  implicit none
  private


  public ::  MiniBooNEenergy, MiniBooNEenergyBARNU
  public ::  CCQE_recQs, CCQE_recQs_Delta, CCQE_recEnergy, CCQE_recEnergy_Delta
  public ::  ANLenergy,BNLenergy
  public ::  K2Kenergy, K2K_recEnergy, K2K_recQs
  public ::  NOVAenergy_FD, NOVAenergy_ND, T2Kenergy_ND,T2Konaxisenergy
  public ::  MINOSenergyNU_fluxNU, MINOSenergyBARNU_fluxNU, &
             & MINOSenergyNU_fluxBARNU, MINOSenergyBARNU_fluxBARNU
  public ::  uniformFlux
  public ::  MINERVAenergyNU
  public ::  MINERVAenergyBARNU
  public ::  DUNEenergyNU
  public ::  DUNEenergyBARNU
  public ::  LBNOenergyNU
  public ::  NOMADenergyNU
  public ::  BNBenergyNUe, BNBenergyNUebar,  BNBenergyNUmu, BNBenergyNUmubar
  public ::  MINERVAenergy

  logical, save :: firsttime=.true.


  ! used to initialize the module via namelist nl_neutrino_energyFlux
  logical, save :: initFlag=.true.

  !variables used for CCQE energy and Q2 reconstruction:

  !****************************************************************************
  !****g* expNeutrinofluxes/Eb
  ! SOURCE
  !
  real, save    :: Eb=0.034
  ! PURPOSE
  ! contant binding energy used for energy and Q2 reconstruction based on
  ! QE scattering kinematics
  !
  !****************************************************************************

  !****************************************************************************
  !****g* expNeutrinofluxes/Eflux_min
  ! SOURCE
  !
  real, save    :: Eflux_min=0.2
  ! PURPOSE
  ! minimum energy for uniform flux distribution
  !
  ! minimum and maximum energies for the uniform neutrino flux (nuExp=10
  ! in the namelist neutrino_induced)
  ! can be changed in the namelist nl_neutrino_energyFlux
  !****************************************************************************

  !****************************************************************************
  !****g* expNeutrinofluxes/Eflux_max
  ! SOURCE
  !
  real, save    :: Eflux_max=2.5
  !
  ! PURPOSE
  ! maximum energy for uniform flux distribution
  !
  ! minimum and maximum energies for the uniform neutrino flux (nuExp=10
  ! in the namelist neutrino_induced)
  ! can be changed in the namelist nl_neutrino_energyFlux
  !****************************************************************************



contains

  !****************************************************************************
  !****s* expNeutrinofluxes/readinput
  ! NAME
  ! subroutine readinput
  ! INPUTS
  ! NONE
  ! OUTPUT
  ! NONE
  ! PURPOSE
  ! This subroutine reads out the namelist "nl_neutrino_energyFlux".
  ! Only called once to initialize the module.
  !****************************************************************************
  subroutine readInput
    use output, only: Write_ReadingInput

    integer :: ios

    !**************************************************************************
    !****n* expNeutrinofluxes/nl_neutrino_energyFlux
    ! NAME
    ! NAMELIST nl_neutrino_energyFlux
    ! PURPOSE
    ! This Namelist includes:
    ! * Eb
    ! * Eflux_min
    ! * Eflux_max
    !**************************************************************************
    NAMELIST /nl_neutrino_energyFlux/ Eb, Eflux_min, Eflux_max
    call Write_ReadingInput('nl_neutrino_energyFlux',0)
    rewind(5)
    read(5,nml=nl_neutrino_energyFlux,IOSTAT=ios)
    call Write_ReadingInput("nl_neutrino_energyFlux",0,ios)

    write(*,*) 'Only valid for uniform flux distribution: Eflux_min is', Eflux_min
    write(*,*) 'Only valid for uniform flux distribution: Eflux_max is', Eflux_max
    write(*,*) 'In CCQE energy reconstruction,  Eb is set to', Eb

    write(*,*) 'In CCQE energy reconstruction,  Eb is set to', Eb
    call Write_ReadingInput('nl_neutrino_energyFlux',1)
  end subroutine readInput


  !****************************************************************************
  !****f* expNeutrinofluxes/MiniBooNEenergy
  ! NAME
  ! real function MiniBooNEenergy()
  !
  ! PURPOSE
  ! This function gives back the neutrino energy for the MiniBooNE experiment.
  ! It determines the energy randomly weighted with the flux.
  ! Flux is taken from http://www-boone.fnal.gov/for_physicists/data_release/flux/
  !  and normalized to 1
  ! paper for reference
  ! A. A. Aguilar-Arevalo et al., "The Neutrino Flux Prediction at MiniBooNE"
  ! Phys. Rev. D. 79, 072002 (2009)
  !****************************************************************************
  real function MiniBooNEenergy()
    use random, only: rn
    use inputGeneral, only: path_To_Input

    real :: v,w,y,x
    !real,parameter:: enumax=2.525
    !real,parameter:: enumin=0.025
    !real,parameter:: ymax=0.885
    character(100) ::fileName
    integer :: status
    real, dimension(145),save :: enu, flux
    integer :: j, j0

    real :: enumax, enumin, ymax, z

    if (firsttime) then
       j=1
       fileName=trim(path_to_Input)//'/neutrino/MiniBooNE-flux.dat'
       open(13,file=filename ,status='old',action='read',iostat=status)
       if (status==0) then
          do
             read(13,*,iostat=status) enu(j),flux(j)
             if (status/=0) exit
             j=j+1
          end do
          if (status>0) then
             write(*,*) 'error reading file'
             stop
          else
             write(*,*) 'file read sucessful'
          end if
       else
          write(*,*) 'problems with file'
       end if
       close(13)
       firsttime=.false.
    end if

    do
       z=rn() ! if z<0.04180791 generate flux above 2 GeV, otherwise below 2 GeV
       if (z>0.04180791) then
          enumin=0.025
          enumax=2.0
          ymax=0.885
          j0=1
       else
          enumin=2.0
          enumax=3.975
          ymax=0.037
          j0=40
       end if
       v=rn()
       w=rn()
       x=enumin+v*(enumax-enumin)
       j=j0
       do
          if (x.lt.enu(j)) exit
          j=j+1
       end do
       y = flux(j-1) + (x - enu(j-1))*(flux(j)-flux(j-1))/(enu(j)-enu(j-1))
       if (w.lt.y/ymax) exit
    end do
    MiniBooNEenergy=x

  end function MiniBooNEenergy



  !****************************************************************************
  !****f* expNeutrinofluxes/MiniBooNEenergyBARNU
  ! NAME
  ! real function MiniBooNEenergyBARNU()
  !
  ! PURPOSE
  ! This function gives back the antineutrino energy for the MiniBooNE
  ! experiment in antineutrino mode (=negartive polarity).
  ! It determines the energy randomly weighted with the flux.
  ! Flux is taken from  http://www-boone.fnal.gov/for_physicists/data_release/flux/
  ! paper for reference
  ! A. A. Aguilar-Arevalo et al., "The Neutrino Flux Prediction at MiniBooNE"
  ! Phys. Rev. D. 79, 072002 (2009)
  !****************************************************************************
  real function MiniBooNEenergyBARNU()
    use random, only: rn
    use inputGeneral, only: path_To_Input

    real :: v,w,y,x
    character(100) ::fileName
    integer :: status
    real, dimension(104),save :: enu, flux
    integer :: j, j0

    real :: enumax, enumin, ymax, z

    if (firsttime) then
       j=1
       fileName=trim(path_to_Input)//'/neutrino/MiniBooNE-flux-barnu.dat'
       open(13,file=filename ,status='old',action='read',iostat=status)
       if (status==0) then
          do
             read(13,*,iostat=status) enu(j),flux(j)
             if (status/=0) exit
             j=j+1
          end do
          if (status>0) then
             write(*,*) 'error reading file'
             stop
          else
             write(*,*) 'file read sucessful'
          end if
       else
          write(*,*) 'problems with file'
       end if
       close(13)
       firsttime=.false.
    end if

    do
       z=rn() ! if z< 1.5*1.731/(1.5*1.731 + 2.5*0.219)=0.82585878
              ! generate flux below 1.525 GeV, otherwise above
       if (z<0.82585878) then
          enumin=0.025
          enumax=1.525
          ymax=1.731
          j0=1
       else
          enumin=1.525
          enumax=4.025
          ymax=0.219
          j0=31
       end if
       v=rn()
       w=rn()
       x=enumin+v*(enumax-enumin)
       j=j0
       do
          if (x.lt.enu(j)) exit
          j=j+1
       end do
       y = flux(j-1) + (x - enu(j-1))*(flux(j)-flux(j-1))/(enu(j)-enu(j-1))
       if (w.lt.y/ymax) exit
    end do
    MiniBooNEenergyBARNU=x

  end function MiniBooNEenergyBARNU



  !****************************************************************************
  !****f* expNeutrinofluxes/CCQE_recQs
  ! NAME
  ! real function CCQE_recQs(k_out)
  !
  ! PURPOSE
  ! This function gives back the reconstruced Qs.
  ! The reconstruction is done as in the experiment, neglecting Fermi motion
  ! (see arXiv:0706.0926v1, eq.(3)), where k_out is the "real" outgoing
  ! lepton momentum.
  !
  !****************************************************************************
  real function CCQE_recQs(k_out)
    use minkowski, only: SP
    use esample, only: Enu_lower_cut, Enu_upper_cut, energylimit_for_QSrec

    real, intent(in), dimension (0:3) :: k_out
    real :: mfsq,Qs_min,Qs_max


 !!!!
 ! Reconstruct Qs only if recEnergy within flux cuts
 !!!!
   if (energylimit_for_Qsrec) then
    if (CCQE_recEnergy(k_out) <  Enu_lower_cut   &
       & .or. CCQE_recEnergy(k_out) >  Enu_upper_cut) then
       CCQE_recQs=30.                        ! 30 = arbitrary, unphysically high value
    end if
   end if

    mfsq=max(SP(k_out,k_out),0.)
    CCQE_recQs=-mfsq+2.*CCQE_recEnergy(k_out)*(k_out(0)-k_out(3))

    Qs_max=-mfsq+2.*CCQE_recEnergy(k_out)*(k_out(0)+sqrt(k_out(0)**2-mfsq))
    Qs_min=-mfsq+2.*CCQE_recEnergy(k_out)*(k_out(0)-sqrt(k_out(0)**2-mfsq))

    if (CCQE_recQs.lt.Qs_min.or.CCQE_recQs.gt.Qs_max) then
       write(*,*) 'CCQE: in QE kinematics reconstructed Qs out of bounds:',&
                  & CCQE_recQs,Qs_min,Qs_max
       CCQE_recQs=30.  !random number, but bigger than usual Qs
    end if

  end function CCQE_recQs



!! This logical formula for Q2 reconstruction assuming Delta kinematics
!! not necessarily used by MiniBooNE
!! This formula was written because it was needed in neutrino analysis routines

  real function CCQE_recQs_Delta(k_out)
    use minkowski, only: SP

    real, intent(in), dimension (0:3) :: k_out
    real :: mfsq,Qs_min,Qs_max

    mfsq=max(SP(k_out,k_out),0.)
    CCQE_recQs_Delta=-mfsq+2.*CCQE_recEnergy_Delta(k_out)*(k_out(0)-k_out(3))

    Qs_max=-mfsq+2.*CCQE_recEnergy_Delta(k_out)*(k_out(0)+sqrt(k_out(0)**2-mfsq))
    Qs_min=-mfsq+2.*CCQE_recEnergy_Delta(k_out)*(k_out(0)-sqrt(k_out(0)**2-mfsq))

    if (CCQE_recQs_Delta.lt.Qs_min.or.CCQE_recQs_Delta.gt.Qs_max) then
       write(*,*) 'CCQE: in Delta kinematics reconstructed Qs out of bounds:', &
                  & CCQE_recQs_Delta,Qs_min,Qs_max
       CCQE_recQs_Delta=5.  !random number, but bigger than usual Qs
    end if

  end function CCQE_recQs_Delta



  !****************************************************************************
  !****f* expNeutrinofluxes/CCQE_recEnergy
  ! NAME
  ! real function CCQE_recEnergy(k_out)
  !
  ! PURPOSE
  ! This function gives back the reconstruced neutrino energy.
  ! The reconstruction is done as in the experiment, neglecting Fermi motion
  ! (see arXiv:0706.0926v1, eq.(4)), where k_out is the "real" outgoing
  ! lepton momentum.
  !
  !****************************************************************************
  real function CCQE_recEnergy(k_out)
    use minkowski, only: SP
    use constants, only: mN

    real, intent(in), dimension (0:3) :: k_out
    real :: mfsq

    if (initFlag) then
       call readInput
       initFlag=.false.
    end if

    mfsq=max(SP(k_out,k_out),0.)
    CCQE_recEnergy=(2.*(MN-Eb)*k_out(0) - &
                   & (Eb**2-2.*MN*Eb+mfsq))/(2.*(MN-Eb-k_out(0)+k_out(3)))

  end function CCQE_recEnergy

  !****************************************************************************
  !****f* expNeutrinofluxes/CCQE_recEnergy_Delta
  ! NAME
  ! real function CCQE_recEnergy_Delta(k_out)
  !
  ! PURPOSE
  ! This function gives back the reconstruced neutrino energy for pions.
  ! The reconstruction is done as in the experiment, neglecting Fermi motion
  ! and binding (see PRL 103, 081801 (2009) eq.(1)), where k_out is
  ! the "real" outgoing lepton momentum.
  !
  !****************************************************************************
  real function CCQE_recEnergy_Delta(k_out)
    use minkowski, only: SP
    use idtable, only: delta
    use ParticleProperties, only: hadron
    use constants, only: mN

    real, intent(in), dimension (0:3) :: k_out
    real :: mfsq,MD

    if (initFlag) then
       call readInput
       initFlag=.false.
    end if

    MD=hadron(delta)%mass
    mfsq=max(SP(k_out,k_out),0.)
    CCQE_recEnergy_Delta=(2.*MN*k_out(0)  &
                         & + MD**2 -MN**2 - mfsq)/(2.*(MN-k_out(0)+k_out(3)))

  end function CCQE_recEnergy_Delta





  !****************************************************************************
  !****f* expNeutrinofluxes/ANLenergy
  ! NAME
  ! real function ANLenergy()
  !
  ! PURPOSE
  ! This function gives back the neutrino energy for the ANL experiment
  ! for QE events
  ! Flux is taken from PRD 16, 3103 (1977), Fig. 7.
  !
  !****************************************************************************
  real function ANLenergy()
    use random, only: rn
    use inputGeneral, only: path_To_Input

    real :: v,w,y,x
    real,parameter:: enumax=5.98
    real,parameter:: enumin=0.125
    real,parameter:: ymax=3.2
    character(100) ::fileName
    integer :: status
    real, dimension(100),save :: enu, flux
    integer :: j

    if (firsttime) then
       j=1
       fileName=trim(path_to_Input)//'/neutrino/anlflux.dat'
       open(13,file=filename ,status='old',action='read',iostat=status)
       if (status==0) then
          do
             read(13,*,iostat=status) enu(j),flux(j)
             if (status/=0) exit
             j=j+1
          end do
          if (status>0) then
             write(*,*) 'error reading file'
             stop
          else
             write(*,*) 'file read sucessful'
          end if
       else
          write(*,*) 'problems with file'
       end if
       close(13)
       firsttime=.false.
    end if

    do
       v=rn()
       w=rn()
       x=enumin+v*(enumax-enumin)
       j=1
       do
          if (x.lt.enu(j)) exit
          j=j+1
       end do
       y = flux(j-1) + (x - enu(j-1))*(flux(j)-flux(j-1))/(enu(j)-enu(j-1))
       if (w.lt.y/ymax) exit
    end do
    ANLenergy=x

  end function ANLenergy



  !****************************************************************************
  !****f* expNeutrinofluxes/BNLenergy
  ! NAME
  ! real function BNLenergy()
  !
  ! PURPOSE
  ! This function gives back the neutrino energy for the BNL experiment.
  ! Flux is taken from K. Furuno, NUINT02 proceedings, available at
  ! http://www.ps.uci.edu/~nuint/proceedings/furuno.pdf
  ! or see Baker et al Phys Rev D23 (1981) 2499, fig.7
  !
  ! NOTES
  ! enumin is for the whole flux and is good for calculating event histograms;
  ! for calculating absolute cross sections BNL used enumin=0.5, which
  ! should be used here
  !****************************************************************************
  real function BNLenergy()
    use random, only: rn
    use inputGeneral, only: path_To_Input

    real :: v,w,y,x
    real,parameter:: enumax=5.492
    real,parameter:: enumin=0.343 ! 0.5
    real,parameter:: ymax=181.
    character(100) ::fileName
    integer :: status
    real, dimension(100),save :: enu, flux
    integer :: j

    if (firsttime) then
       j=1
       fileName=trim(path_to_Input)//'/neutrino/bnlflux.dat'
       open(13,file=filename ,status='old',action='read',iostat=status)
       if (status==0) then
          do
             read(13,*,iostat=status) enu(j),flux(j)
             if (status/=0) exit
             j=j+1
          end do
          if (status>0) then
             write(*,*) 'error reading file'
             stop
          else
             write(*,*) 'file read sucessful'
          end if
       else
          write(*,*) 'problems with file'
       end if
       close(13)
       firsttime=.false.
    end if

    do
       v=rn()
       w=rn()
       x=enumin+v*(enumax-enumin)
       j=1
       do
          if (x.lt.enu(j)) exit
          j=j+1
       end do
       y = flux(j-1) + (x - enu(j-1))*(flux(j)-flux(j-1))/(enu(j)-enu(j-1))
       if (w.lt.y/ymax) exit
    end do
    BNLenergy=x

  end function BNLenergy





  !****************************************************************************
  !****f* expNeutrinofluxes/K2Kenergy
  ! NAME
  ! real function K2Kenergy()
  !
  ! PURPOSE
  ! This function gives back the neutrino energy for the K2K experiment.
  ! Flux is taken from PLB 619 (2005), Fig. 1
  !
  !****************************************************************************
  real function K2Kenergy()
    use random, only: rn
    use inputGeneral, only: path_To_Input

    real :: v,w,y,x
    real,parameter:: enumax=3.9470
    real,parameter:: enumin=0.04525
    real,parameter:: ymax=12.4
    character(100) ::fileName
    integer :: status
    real, dimension(100),save :: enu, flux
    integer :: j

    if (firsttime) then
       j=1
       fileName=trim(path_to_Input)//'/neutrino/k2kflux.dat'
       open(13,file=filename ,status='old',action='read',iostat=status)
       if (status==0) then
          do
             read(13,*,iostat=status) enu(j),flux(j)
             if (status/=0) exit
             j=j+1
          end do
          if (status>0) then
             write(*,*) 'error reading file'
             stop
          else
             write(*,*) 'file read sucessful'
          end if
       else
          write(*,*) 'problems with file'
       end if
       close(13)
       firsttime=.false.
    end if

    do
       v=rn()
       w=rn()
       x=enumin+v*(enumax-enumin)
       j=1
       do
          if (x.lt.enu(j)) exit
          j=j+1
       end do
       y = flux(j-1) + (x - enu(j-1))*(flux(j)-flux(j-1))/(enu(j)-enu(j-1))
       if (w.lt.y/ymax) exit
    end do
    K2Kenergy=x

  end function K2Kenergy

  !****************************************************************************
  !****f* expNeutrinofluxes/K2K_recEnergy
  ! NAME
  ! real function K2K_recEnergy(k_out)
  !
  ! PURPOSE
  ! This function gives back the reconstruced neutrino energy.
  ! The reconstruction is done as in the experiment, neglecting Fermi motion
  ! (see PRL 90, 041801 (2003), eq.(1)), where k_out is the "real" outgoing
  ! lepton momentum.
  !
  ! another formula is used in recent CC-pi0/QE  measurements, see ArXiV 1012.1794
  ! it will be used if an optional parameter W is given (default K2K W=1.483 GeV)
  ! for W=MN the formular coincides with the old one
  !****************************************************************************
  real function K2K_recEnergy(k_out,W)
    use minkowski, only: SP
    use constants, only: mN

    real, intent(in), dimension (0:3) :: k_out
    real :: mfsq
    real, intent(in), optional :: W !

    mfsq=max(SP(k_out,k_out),0.)

    if (present(W)) then
    K2K_recEnergy=( W*W-mfsq +2.*k_out(0)*MN -MN*MN )/2./(MN-k_out(0)+k_out(3))
    else
    K2K_recEnergy=(MN*k_out(0) - mfsq/2.)/(MN-k_out(0)+k_out(3))
    end if

  end function K2K_recEnergy


!! this is logical formular for Q2 reconstruction assuming QE kinematics
!! not necessarily used by K2K
!! This formula was written because it was needed in neutrino analysis routines
  real function K2K_recQs(k_out,W)
    use minkowski, only: SP

    real, intent(in), dimension (0:3) :: k_out
    real :: mfsq
    real, intent(in), optional :: W !

    mfsq=max(SP(k_out,k_out),0.)

    if (present(W)) then
    K2K_recQs=-mfsq+2.*K2K_recEnergy(k_out,W)*(k_out(0)-k_out(3))
    else
    K2K_recQs=-mfsq+2.*K2K_recEnergy(k_out)*(k_out(0)-k_out(3))
    end if

  end function K2K_recQs

  !****************************************************************************
  !****f* expNeutrinofluxes/T2Kenergy_ND
  ! NAME
  ! real function T2Kenergy_ND(Flavor_ID,Process_ID)
  !
  ! PURPOSE
  ! This function gives back the neutrino energy for the T2K ND280  experiment.
  ! Flux is 2.5 degrees off-axis flux for the ND280 detector
  ! implemented is ND280_horn_250kA taken from http://t2k-experiment.org/results/
  !****************************************************************************
 real function T2Kenergy_ND(Flavor_ID,Process_ID)
    use esample, only: read_fluxfile, eneut

    integer,parameter :: NDIM = 250          !maximal dimension of fluxfile
    real, dimension (NDIM), save :: enu,flux
    real, dimension (0:NDIM), save :: sumflux
    character(100) :: fluxfilename
    integer, save :: jmax
    integer :: nswitch,flavor_id,process_id

!   Now reading of flux file from buuinput/neutrinos
!   First, set filename for desired neutrino flavor

     nswitch = sign(1,Process_ID)

     select case (nswitch)

      case (+1)
         select case (Flavor_ID)

           case (1)
           fluxfileName= 'T2K_ND280_250kA-nue.dat'   ! electron

           case (2)
           fluxfileName= 'T2K_ND280_250kA-numu.dat'  ! muon

           case default
           write(*,*) 'flavor and process IDs not compatible:1'
           stop

         end select

      case (-1)
          write(*,*) 'antineutrino fluxes for T2K ND not yet implemented'
          stop

         select case (Flavor_ID)

           case (1)
           fluxfileName= 'T2K_ND280_250kA-anue.dat'    ! anti-electron

           case (2)
           fluxfileName= 'T2K_ND280_250kA-anumu.dat'   ! anti-muon

           case default
           write(*,*) 'flavor and process IDs not compatible:2'
           stop

         end select
      case default
           write(*,*) 'flavor and process IDs not compatible:3'
     end select


    if (firsttime) then
        call read_fluxfile(NDIM,fluxfilename,jmax,enu,flux,sumflux)
        firsttime=.false.
    end if

!   Now call of sampling routine

     T2Kenergy_ND = eneut(NDIM,jmax,sumflux,enu)

 end function T2Kenergy_ND



!  !*************************************************************************
!  !****f* expNeutrinofluxes/T2K_ND_numu_energy
!  ! NAME
!  ! real function T2K_ND_numu_energy()
!  !
!  ! PURPOSE
!  ! This function gives back the neutrino energy for the T2K ND280  experiment.
!  ! Flux is 2.5 degrees off-axis flux for the ND280 detector
!  ! implemented is ND280_horn_250kA taken from http://t2k-experiment.org/results/
!  !*************************************************************************
!   real function T2K_ND_numu_energy()
!    use random, only: rn
!    use inputGeneral, only : path_To_Input
!
!    real :: v,w,y,x
!    real,parameter :: enumax=22.5
!    real,parameter :: enumin=0.0
!    real,parameter :: ymax=1.52e12
!    character(100) :: fileName
!    integer :: status
!  ! dimension= max number of energy values in flux file
!    real, dimension(221),save :: enu, flux
!    integer :: j,jmax
!
!
!     if(firsttime) then
!
!       j=1
!       fileName=trim(path_to_Input)//'/neutrino/T2K_ND280_250kA-numu.dat'
!       open(13,file=filename ,status='old',action='read',iostat=status)
!       if(status==0) then
!          do
!             read(13,*,iostat=status) enu(j),flux(j)
!             if(status/=0) exit
!             j=j+1
!          end do
!          jmax=j-1
!          if(status>0) then
!             write(*,*)'error reading file'
!             stop
!          else
!             write(*,*)'file read sucessful, jmax=',jmax
!          end if
!       else
!          write(*,*)'problems with file'
!       end if
!       close(13)
!
!       firsttime=.false.
!    end if
!
!    do
!       v=rn()
!       w=rn()
!       x=enumin+v*(enumax-enumin)
!
!
!       j=2 ! avoid interpolation problems!
!       do
!          if(x.lt.enu(j)) exit
!          j=j+1
!       end do
!       y = flux(j-1) + (x - enu(j-1))*(flux(j)-flux(j-1))/(enu(j)-enu(j-1))
!
!       if(w.lt.y/ymax) exit
!    end do
!    T2K_ND_numu_energy=x
!
!  end function T2K_ND_numu_energy
!
!


  !****************************************************************************
  !****f* expNeutrinofluxes/MINOSenergyNU_fluxNU
  ! NAME
  ! real function MINOSenergyNU_fluxNU()
  !
  ! PURPOSE
  ! This function returns the muon-neutrino energy for the MINOS neutrino experiment
  ! (NUMI low-energy flux) in neutrino mode.
  ! Flux sent to us by Minerva team (Steve Dytman)
  !
  !****************************************************************************
  real function MINOSenergyNU_fluxNU()
    use random, only: rn
    use inputGeneral, only: path_To_Input

    real :: v,w,y,x
    character(100) :: fileName
    integer :: status
    real, dimension(302),save :: enu, flux
    integer :: j, j0
    real :: enumax, enumin, ymax, z

     if (firsttime) then
       j=1
       fileName=trim(path_to_Input)//'/neutrino/MINOS-numu-numode-Minerva.dat'
       open(13,file=filename,status='old',action='read',iostat=status)
       if (status==0) then
          do
             read(13,*,iostat=status) enu(j),flux(j)
             if (status/=0) exit
             j=j+1
          end do
          if (status>0) then
             write(*,*) 'error reading file'
             stop
          else
             write(*,*) 'file read sucessful'
          end if
       else
          write(*,*) 'problems with file'
       end if
       close(13)
       firsttime=.false.
    end if

    do
       z=rn() ! if z< (0.3056*10)/(0.3056*10 + 0.00954*50 )=0.864987
              ! generate flux below 10.2 GeV
       if (z<0.864987) then
          enumin=0.2
          enumax=10.2
          ymax=0.3056
          j0=1
       else
          enumin=10.2
          enumax=60.2
          ymax=0.00954
          j0=51
       end if
       v=rn()
       w=rn()
       x=enumin+v*(enumax-enumin)
       j=j0
       do
          if (x.lt.enu(j)) exit
          j=j+1
       end do
       y = flux(j-1) + (x - enu(j-1))*(flux(j)-flux(j-1))/(enu(j)-enu(j-1))
       if (w.lt.y/ymax) exit
    end do
    MINOSenergyNU_fluxNU=x

  end function MINOSenergyNU_fluxNU




  !****************************************************************************
  !****f* expNeutrinofluxes/MINOSenergyBARNU_fluxNU
  ! NAME
  ! real function MINOSenergyBARNU_fluxNU()
  !
  ! PURPOSE
  ! This function returns the muon-antineutrino energy for the MINOS neutrino experiment
  ! (NUMI low-energy flux) in neutrino mode.
  ! Flux sent to us by LAr team (Ornella Palamara)
  !
  !****************************************************************************
  real function MINOSenergyBARNU_fluxNU()
    use random, only: rn
    use inputGeneral, only: path_To_Input

    real :: v,w,y,x
    character(100) :: fileName
    integer :: status
    real, dimension(289),save :: enu, flux
    integer :: j, j0
    real :: enumax, enumin, ymax, z

     if (firsttime) then
       j=1
       fileName=trim(path_to_Input)//'/neutrino/MINOS-barnumu-numode.dat'
       open(13,file=filename ,status='old',action='read',iostat=status)
       if (status==0) then
          do
             read(13,*,iostat=status) enu(j),flux(j)
             if (status/=0) exit
             j=j+1
          end do
          if (status>0) then
             write(*,*) 'error reading file'
             stop
          else
             write(*,*) 'file read sucessful'
          end if
       else
          write(*,*) 'problems with file'
       end if
       close(13)
       firsttime=.false.
    end if

    do
       z=rn() ! if z< (0.0395*7)/(0.0395*7 + 0.0116*53)=0.310221
              !  generate flux below 7.125 GeV
       if (z<0.310221) then
          enumin=0.125
          enumax=7.125
          ymax=0.0395
          j0=1
       else
          enumin=7.125
          enumax=60.125
          ymax=0.0116
          j0=29
       end if
       v=rn()
       w=rn()
       x=enumin+v*(enumax-enumin)
       j=j0
       do
          if (x.lt.enu(j)) exit
          j=j+1
       end do
       y = flux(j-1) + (x - enu(j-1))*(flux(j)-flux(j-1))/(enu(j)-enu(j-1))
       if (w.lt.y/ymax) exit
    end do
    MINOSenergyBARNU_fluxNU=x

  end function MINOSenergyBARNU_fluxNU





  !****************************************************************************
  !****f* expNeutrinofluxes/MINOSenergyNU_fluxBARNU
  ! NAME
  ! real function MINOSenergyNU_fluxBARNU()
  !
  ! PURPOSE
  ! This function gives back the muon-neutrino energy for
  ! NUMI low-energy flux in antineutrino mode.
  ! Flux sent to us by LAr team (Ornella Palamara)
  !
  !****************************************************************************
  real function MINOSenergyNU_fluxBARNU()
    use random, only: rn
    use inputGeneral, only: path_To_Input

    real :: v,w,y,x
    character(100) :: fileName
    integer :: status
    real, dimension(374),save :: enu, flux
    integer :: j, j0
    real :: enumax, enumin, ymax, z

     if (firsttime) then
       j=1
       fileName=trim(path_to_Input)//'/neutrino/MINOS-numu-barnumode.dat'
       open(13,file=filename ,status='old',action='read',iostat=status)
       if (status==0) then
          do
             read(13,*,iostat=status) enu(j),flux(j)
             if (status/=0) exit
             j=j+1
          end do
          if (status>0) then
             write(*,*) 'error reading file'
             stop
          else
             write(*,*) 'file read sucessful'
          end if
       else
          write(*,*) 'problems with file'
       end if
       close(13)
       firsttime=.false.
    end if

    do
       z=rn() ! if z< (0.03197*9)/(0.03197*9 + 0.008615*71)=0.319915
              !  generate flux below 9.125 GeV
       if (z<0.319915) then
          enumin=0.125
          enumax=9.125
          ymax=0.03197
          j0=1
       else
          enumin=9.125
          enumax=80.125
          ymax=0.008615
          j0=37
       end if
       v=rn()
       w=rn()
       x=enumin+v*(enumax-enumin)
       j=j0
       do
          if (x.lt.enu(j)) exit
          j=j+1
       end do
       y = flux(j-1) + (x - enu(j-1))*(flux(j)-flux(j-1))/(enu(j)-enu(j-1))
       if (w.lt.y/ymax) exit
    end do
    MINOSenergyNU_fluxBARNU=x

  end function MINOSenergyNU_fluxBARNU



  !****************************************************************************
  !****f* expNeutrinofluxes/MINOSenergyBARNU_fluxBARNU
  ! NAME
  ! real function MINOSenergyBARNU_fluxBARNU()
  !
  ! PURPOSE
  ! This function gives back the muon-antineutrino energy
  ! NUMI low-energy flux in antineutrino mode.
  ! Flux sent to us by LAr team (Ornella Palamara)
  !
  !****************************************************************************
  real function MINOSenergyBARNU_fluxBARNU()
    use random, only: rn
    use inputGeneral, only: path_To_Input

    real :: v,w,y,x
    !real,parameter :: enumax=37.5
    !real,parameter :: enumin=0.5
    !real,parameter :: ymax=3.09
    character(100) :: fileName
    integer :: status
    real, dimension(226),save :: enu, flux
    integer :: j, j0
    real :: enumax, enumin, ymax, z

     if (firsttime) then
       j=1
       fileName=trim(path_to_Input)//'/neutrino/MINOS-barnumu-barnumode.dat'
       open(13,file=filename ,status='old',action='read',iostat=status)
       if (status==0) then
          do
             read(13,*,iostat=status) enu(j),flux(j)
             if (status/=0) exit
             j=j+1
          end do
          if (status>0) then
             write(*,*) 'error reading file'
             stop
          else
             write(*,*) 'file read sucessful'
          end if
       else
          write(*,*) 'problems with file'
       end if
       close(13)
       firsttime=.false.
    end if

    do
       z=rn() ! if z< (0.09232*10)/(0.09232*10 + 0.0116*30)=0.72624292
              !  generate flux below 10.125 GeV
       if (z<0.72624292) then
          enumin=0.125
          enumax=10.125
          ymax=0.09232
          j0=1
       else
          enumin=10.125
          enumax=40.125
          ymax=0.0115
          j0=41
       end if
       v=rn()
       w=rn()
       x=enumin+v*(enumax-enumin)
       j=j0
       do
          if (x.lt.enu(j)) exit
          j=j+1
       end do
       y = flux(j-1) + (x - enu(j-1))*(flux(j)-flux(j-1))/(enu(j)-enu(j-1))
       if (w.lt.y/ymax) exit
    end do
    MINOSenergyBARNU_fluxBARNU=x

  end function MINOSenergyBARNU_fluxBARNU



  !****************************************************************************
  !****f* expNeutrinofluxes/uniformFlux
  ! NAME
  ! real function uniformFlux()
  !
  ! PURPOSE
  ! generated uniform flux from Eflux_min to Eflux_max
  ! (see namelist nl_neturino_energyFlux)
  !
  !****************************************************************************
  real function uniformFlux()
    use random
    if (initFlag) then
       call readInput
       initFlag=.false.
    end if
    uniformFlux = Eflux_min + (Eflux_max-Eflux_min)*rn()
  end  function uniformFlux




!******************************************************************************
 !****f* expNeutrinofluxes/MINERVAenergyNU
 ! NAME
 ! real function MINERVAenergyNU()
 !
 ! PURPOSE
 ! This function returns the sampled neutrino energy for the MINERvA experiment
 ! in neutrino mode.
 ! Flux is obtained from B. Tice, June 2013, now outdated
 !*****************************************************************************
 real function MINERVAenergyNU()
    use esample, only: read_fluxfile, eneut

    integer,parameter :: NDIM = 350          !maximal dimension of fluxfile
    real, dimension (NDIM), save :: enu,flux
    real, dimension (0:NDIM), save :: sumflux
    character(100) :: fluxfilename
    integer, save :: jmax

!   Now reading of flux file from buuinput/neutrinos

    fluxfilename = 'Minerva_neutrino.dat'

    if (firsttime) then
        call read_fluxfile(NDIM,fluxfilename,jmax,enu,flux,sumflux)
        firsttime=.false.
    end if

!   Now call of sampling routine

     MINERVAenergyNU = eneut(NDIM,jmax,sumflux,enu)

 end function MINERVAenergyNU



 !*****************************************************************************
 !****f* expNeutrinofluxes/MINERVAenergyBARNU
 ! NAME
 ! real function MINERVAenergyBARNU()
 !
 ! PURPOSE
 ! This function returns the sampled antineutrino energy for the MINERvA experiment
 ! in antineutrino mode.
 ! Flux is obtained from B. Tice, June 2013, now outdated
 !*****************************************************************************
 real function MINERVAenergyBARNU()
    use esample, only: read_fluxfile, eneut

    integer,parameter :: NDIM = 350          !maximal dimension of fluxfile
    real, dimension (NDIM), save :: enu,flux
    real, dimension (0:NDIM), save :: sumflux
    character(100) :: fluxfilename
    integer, save :: jmax

!   Now reading of flux file from buuinput/neutrinos

    fluxfilename = 'Minerva_antineutrino.dat'

    if (firsttime) then
        call read_fluxfile(NDIM,fluxfilename,jmax,enu,flux,sumflux)
        firsttime=.false.
    end if

!   Now call of sampling routine

     MINERVAenergyBARNU = eneut(NDIM,jmax,sumflux,enu)

 end function MINERVAenergyBARNU



  !****************************************************************************
  !****f* expNeutrinofluxes/DUNEenergyNU
  ! NAME
  ! real function DUNEenergyNU()
  !
  ! PURPOSE
  ! This function returns the sampled neutrino energy for the DUNE  experiment
  ! in neutrino mode.
  ! Flux is obtained from  
  ! http://home.fnal.gov/~ljf26/DUNE2015CDRFluxes/
  !****************************************************************************
  real function DUNEenergyNU()
    use esample, only: read_fluxfile, eneut

    integer,parameter :: NDIM = 250         !maximal dimension of fluxfile
    real, dimension (NDIM), save :: enu,flux
    real, dimension (0:NDIM), save :: sumflux
    character(100) :: fluxfilename
    integer, save :: jmax


!   Now reading of flux file from buuinput/neutrinos

    fluxfilename = 'DUNE-nu-NDFlux.dat'

    if (firsttime) then
        call read_fluxfile(NDIM,fluxfilename,jmax,enu,flux,sumflux)
        firsttime=.false.
    end if

!   Now call of sampling routine

     DUNEenergyNU = eneut(NDIM,jmax,sumflux,enu)

 end function DUNEenergyNU


 !*****************************************************************************
 !****f* expNeutrinofluxes/DUNEenergyBARNU
 ! NAME
 ! real function DUNEenergyBARNU()
 !
 ! PURPOSE
 ! This function returns the sampled antineutrino energy for the DUNE experiment
 ! in antineutrino mode.
 ! Flux is obtained from http://home.fnal.gov/~ljf26/DUNE2015CDRFluxes/
 !*****************************************************************************
 real function DUNEenergyBARNU()
    use esample, only: read_fluxfile, eneut

    integer,parameter :: NDIM = 250          !maximal dimension of fluxfile
    real, dimension (NDIM), save :: enu,flux
    real, dimension (0:NDIM), save :: sumflux
    character(100) :: fluxfilename
    integer, save :: jmax



!   Now reading of flux file from buuinput/neutrinos

    fluxfilename = 'DUNE-nubar-NDFlux.dat'

    if (firsttime) then
        call read_fluxfile(NDIM,fluxfilename,jmax,enu,flux,sumflux)
        firsttime=.false.
    end if

!   Now call of sampling routine

     DUNEenergyBARNU = eneut(NDIM,jmax,sumflux,enu)

 end function DUNEenergyBARNU

 !*****************************************************************************
 !****f* expNeutrinofluxes/LBNOenergyNU
 ! NAME
 ! real function LBNOenergyNU()
 !
 ! PURPOSE
 ! This function returns the sampled neutrino energy for the LBNO  experiment
 ! in antineutrino mode.
 ! Flux is obtained from Dario Autiero, May 2014
 !*****************************************************************************
 real function LBNOenergyNU()
    use esample, only: read_fluxfile, eneut

    integer,parameter :: NDIM = 350          !maximal dimension of fluxfile
    real, dimension (NDIM), save :: enu,flux
    real, dimension (0:NDIM), save :: sumflux
    character(100) :: fluxfilename
    integer, save :: jmax



!   Now reading of flux file from buuinput/neutrinos

    fluxfilename = 'LBNO-nu_mu-nu_mu-mode.dat'

    if (firsttime) then
        call read_fluxfile(NDIM,fluxfilename,jmax,enu,flux,sumflux)
        firsttime=.false.
    end if

!   Now call of sampling routine

     LBNOenergyNU = eneut(NDIM,jmax,sumflux,enu)

 end function LBNOenergyNU

 !*****************************************************************************
 !****f* expNeutrinofluxes/NOMADenergyNU
 ! NAME
 ! real function NOMADenergyNU()
 !
 ! PURPOSE
 ! This function returns the sampled neutrino energy for the NOMAD  experiment
 ! in antineutrino mode.
 ! Flux is obtained from Roberto Petti, Jan.  2014
 !*****************************************************************************
 real function NOMADenergyNU()
    use esample, only: read_fluxfile, eneut

    integer,parameter :: NDIM = 610          !maximal dimension of fluxfile
    real, dimension (NDIM), save :: enu,flux
    real, dimension (0:NDIM), save :: sumflux
    character(100) :: fluxfilename
    integer, save :: jmax



!   Now reading of flux file from buuinput/neutrinos

    fluxfilename = 'NOMAD-nu_mu-nu_mu-mode.dat'

    if (firsttime) then
        call read_fluxfile(NDIM,fluxfilename,jmax,enu,flux,sumflux)
        firsttime=.false.
    end if

!   Now call of sampling routine

     NOMADenergyNU = eneut(NDIM,jmax,sumflux,enu)

 end function NOMADenergyNU

 !*****************************************************************************
 !****f* expNeutrinofluxes/BNBenergyNUe
 ! NAME
 ! real function BNBenergyNUe()
 !
 ! PURPOSE
 ! This function returns the sampled neutrino energy for the BNB experiment
 ! in electron-neutrino mode.
 ! Flux is obtained from Anne Schukraft, December 2015
 !*****************************************************************************
 real function BNBenergyNUe()
    use esample, only: read_fluxfile, eneut

    integer,parameter :: NDIM = 210          !maximal dimension of fluxfile
    real, dimension (NDIM), save :: enu,flux
    real, dimension (0:NDIM), save :: sumflux
    character(100) :: fluxfilename
    integer, save :: jmax



!   Now reading of flux file from buuinput/neutrinos

    fluxfilename = 'BNB_nue.dat'

    if (firsttime) then
        call read_fluxfile(NDIM,fluxfilename,jmax,enu,flux,sumflux)
        firsttime=.false.
    end if

!   Now call of sampling routine

     BNBenergyNUe = eneut(NDIM,jmax,sumflux,enu)

 end function BNBenergyNUe

 !*****************************************************************************
 !****f* expNeutrinofluxes/BNBenergyNUebar
 ! NAME
 ! real function BNBenergyNUebar()
 !
 ! PURPOSE
 ! This function returns the sampled neutrino energy for the BNB experiment
 ! in electron-antineutrino mode.
 ! Flux is obtained from Anne Schukraft, December 2015
 !*****************************************************************************
 real function BNBenergyNUebar()
    use esample, only: read_fluxfile, eneut

    integer,parameter :: NDIM = 210          !maximal dimension of fluxfile
    real, dimension (NDIM), save :: enu,flux
    real, dimension (0:NDIM), save :: sumflux
    character(100) :: fluxfilename
    integer, save :: jmax



!   Now reading of flux file from buuinput/neutrinos

    fluxfilename = 'BNB_nuebar.dat'

    if (firsttime) then
        call read_fluxfile(NDIM,fluxfilename,jmax,enu,flux,sumflux)
        firsttime=.false.
    end if

!   Now call of sampling routine

     BNBenergyNUebar = eneut(NDIM,jmax,sumflux,enu)

 end function BNBenergyNUebar

  !****************************************************************************
 !****f* expNeutrinofluxes/BNBenergyNUmu
 ! NAME
 ! real function BNBenergyNUmu()
 !
 ! PURPOSE
 ! This function returns the sampled neutrino energy for the BNB experiment
 ! in muon-neutrino mode.
 ! Flux is obtained from Anne Schukraft, December 2015
 !*****************************************************************************
 real function BNBenergyNUmu()
    use esample, only: read_fluxfile, eneut

    integer,parameter :: NDIM = 210          !maximal dimension of fluxfile
    real, dimension (NDIM), save :: enu,flux
    real, dimension (0:NDIM), save :: sumflux
    character(100) :: fluxfilename
    integer, save :: jmax



!   Now reading of flux file from buuinput/neutrinos

    fluxfilename = 'BNB_numu.dat'

    if (firsttime) then
        call read_fluxfile(NDIM,fluxfilename,jmax,enu,flux,sumflux)
        firsttime=.false.
    end if

!   Now call of sampling routine

     BNBenergyNUmu = eneut(NDIM,jmax,sumflux,enu)

 end function BNBenergyNUmu

 !*****************************************************************************
 !****f* expNeutrinofluxes/BNBenergyNUmubar
 ! NAME
 ! real function BNBenergyNUmubar()
 !
 ! PURPOSE
 ! This function returns the sampled neutrino energy for the BNB experiment
 ! in muon-antineutrino mode.
 ! Flux is obtained from Anne Schukraft, December 2015
 !*****************************************************************************
 real function BNBenergyNUmubar()
    use esample, only: read_fluxfile, eneut

    integer,parameter :: NDIM = 210          !maximal dimension of fluxfile
    real, dimension (NDIM), save :: enu,flux
    real, dimension (0:NDIM), save :: sumflux
    character(100) :: fluxfilename
    integer, save :: jmax



!   Now reading of flux file from buuinput/neutrinos

    fluxfilename = 'BNB_numubar.dat'

    if (firsttime) then
        call read_fluxfile(NDIM,fluxfilename,jmax,enu,flux,sumflux)
        firsttime=.false.
    end if

!   Now call of sampling routine

     BNBenergyNUmubar = eneut(NDIM,jmax,sumflux,enu)

 end function BNBenergyNUmubar



!******************************************************************************
 !****f* expNeutrinofluxes/NOVAenergy_ND
 ! NAME
 ! real function NOVAenergy_ND(flavor_id,process_id)
 !
 ! PURPOSE
 ! This function returns the sampled neutrino energy for the NOvA experiment
 ! at the Near Detector, using the FHC files for neutrinos and RHC files for
 ! antineutrinos
 ! Flux is obtained from Jonathan Paley, March 2016
 !*****************************************************************************
 real function NOVAenergy_ND(Flavor_ID,Process_ID)
    use esample, only: read_fluxfile, eneut

    integer,parameter :: NDIM = 601          !maximal dimension of fluxfile
    real, dimension (NDIM), save :: enu,flux
    real, dimension (0:NDIM), save :: sumflux
    character(100) :: fluxfilename
    integer, save :: jmax
    integer :: nswitch,flavor_id,process_id



!   Now reading of flux file from buuinput/neutrinos
!   First, set filename for desired neutrino flavor

     nswitch = sign(1,Process_ID)

     select case (nswitch)

      case (+1)
         select case (Flavor_ID)

           case (1)
           fluxfileName= 'NOvA-ND-FHC-nue.dat'   ! electron

           case (2)
           fluxfileName= 'NOvA-ND-FHC-numu.dat'  ! muon

           case default
           write(*,*) 'flavor and process IDs not compatible:1'
           stop

         end select

      case (-1)

         select case (Flavor_ID)

           case (1)
           fluxfileName= 'NOvA-ND-RHC-anue.dat'    ! anti-electron

           case (2)
           fluxfileName= 'NOvA-ND-RHC-anumu.dat'   ! anti-muon

           case default
           write(*,*) 'flavor and process IDs not compatible:2'
           stop

         end select
      case default
           write(*,*) 'flavor and process IDs not compatible:3'
     end select


    if (firsttime) then
        call read_fluxfile(NDIM,fluxfilename,jmax,enu,flux,sumflux)
        firsttime=.false.
    end if

!   Now call of sampling routine

      NOVAenergy_ND = eneut(NDIM,jmax,sumflux,enu)

 end function NOVAenergy_ND

 !*****************************************************************************
 !****f* expNeutrinofluxes/NOVAenergyNU_FD
 ! NAME
 ! real function NOVAenergyNU_FD(flavor_id,process_id)
 !
 ! PURPOSE
 ! This function returns the sampled neutrino energy for the NOvA experiment
 ! at the Far Detector, using the FHC files for neutrinos and RHC files for
 ! antineutrinos
 ! Flux is obtained from Jonathan Paley, March 2016
 !*****************************************************************************
 real function NOVAenergy_FD(Flavor_ID,Process_ID)
    use esample, only: read_fluxfile, eneut

    integer,parameter :: NDIM = 601          !maximal dimension of fluxfile
    real, dimension (NDIM), save :: enu,flux
    real, dimension (0:NDIM), save :: sumflux
    character(100) :: fluxfilename
    integer, save :: jmax
    integer :: nswitch,flavor_id,process_id



!   Now reading of flux file from buuinput/neutrinos
!   First, set filename for desired neutrino flavor

     nswitch = sign(1,Process_ID)

     select case (nswitch)

      case (+1)
         select case (Flavor_ID)

           case (1)
           fluxfileName= 'NOvA-FD-FHC-nue.dat'   ! electron

           case (2)
           fluxfileName= 'NOvA-FD-FHC-numu.dat'  ! muon

           case default
           write(*,*) 'flavor and process IDs not compatible:1'
           stop

         end select

      case (-1)

         select case (Flavor_ID)

           case (1)
           fluxfileName= 'NOvA-FD-RHC-anue.dat'    ! anti-electron

           case (2)
           fluxfileName= 'NOvA-FD-RHC-anumu.dat'   ! anti-muon

           case default
           write(*,*) 'flavor and process IDs not compatible:2'
           stop

         end select
      case default
           write(*,*) 'flavor and process IDs not compatible:3'
     end select


    if (firsttime) then
        call read_fluxfile(NDIM,fluxfilename,jmax,enu,flux,sumflux)
        firsttime=.false.
    end if

!   Now call of sampling routine

      NOVAenergy_FD = eneut(NDIM,jmax,sumflux,enu)

 end function NOVAenergy_FD

 !*****************************************************************************
 !****f* expNeutrinofluxes/MINERvAenergy
 ! NAME
 ! real function MINERvAenergy(flavor_id,process_id)
 !
 ! PURPOSE
 ! This function returns the sampled neutrino energy for the MINERvA experiment
 ! at the Far Detector,
 ! flux is obtained from  Phys.Rev. D94 (2016) no.9, 092005,
 ! Addendum: Phys.Rev. D95 (2017) no.3, 039903
 !*****************************************************************************
 real function MINERvAenergy(Flavor_ID,Process_ID)
    use esample, only: read_fluxfile, eneut

    integer,parameter :: NDIM = 200          !maximal dimension of fluxfile
    real, dimension (NDIM), save :: enu,flux
    real, dimension (0:NDIM), save :: sumflux
    character(100) :: fluxfilename
    integer, save :: jmax
    integer :: nswitch,flavor_id,process_id



!   Now reading of flux file from buuinput/neutrinos
!   First, set filename for desired neutrino flavor

     nswitch = sign(1,Process_ID)

     select case (nswitch)

      case (+1)
         select case (Flavor_ID)

           case (1)
           fluxfileName= 'Minerva_enuFlux.dat'   ! electron

           case (2)
           fluxfileName= 'Minerva_munuFlux.dat'  ! muon

           case default
           write(*,*) 'flavor and process IDs not compatible:1'
           stop

         end select

      case (-1)

         select case (Flavor_ID)

           case (1)
           fluxfileName= 'Minerva_antienuFlux.dat'    ! anti-electron

           case (2)
           fluxfileName= 'Minerva_antimunuFlux.dat'   ! anti-muon

           case default
           write(*,*) 'flavor and process IDs not compatible:2'
           stop

         end select
      case default
           write(*,*) 'flavor and process IDs not compatible:3'
     end select


    if (firsttime) then
        call read_fluxfile(NDIM,fluxfilename,jmax,enu,flux,sumflux)
        firsttime=.false.
    end if

!   Now call of sampling routine

      MINERvAenergy = eneut(NDIM,jmax,sumflux,enu)

 end function MINERvAenergy


 !*****************************************************************************
 !****f* expNeutrinofluxes/T2Konaxisenergy
 ! NAME
 ! real function T2Konaxisenergy(flavor_id,process_id)
 !
 ! PURPOSE
 ! This function returns the sampled neutrino energy for the T2K on axis
 ! experiment
 !*****************************************************************************
 real function T2Konaxisenergy(Flavor_ID,Process_ID)
    use esample, only: read_fluxfile, eneut

    integer,parameter :: NDIM = 50          !maximal dimension of fluxfile
    real, dimension (NDIM), save :: enu,flux
    real, dimension (0:NDIM), save :: sumflux
    character(100) :: fluxfilename
    integer, save :: jmax
    integer :: nswitch,flavor_id,process_id



!   Now reading of flux file from buuinput/neutrinos
!   First, set filename for desired neutrino flavor

     nswitch = sign(1,Process_ID)

     select case (nswitch)

      case (+1)
         select case (Flavor_ID)

           case (1)
           write(*,*) 'T2k on axis flux for nue not implemented'
           stop
           fluxfileName= 'T2Kon_nueFlux.dat'   ! electron


           case (2)
           fluxfileName= 'T2Kon_numuFlux.dat'  ! muon

           case default
           write(*,*) 'flavor and process IDs not compatible:1'
           stop

         end select

      case (-1)

         select case (Flavor_ID)

           case (1)
           write(*,*) 'T2k on axis flux for antinue not implemented'
           stop
           fluxfileName= 'T2Kon_antienuFlux.dat'    ! anti-electron


           case (2)
           write(*,*) 'T2k on axis flux for antinumu not implemented'
           stop
           fluxfileName= 'T2Kon_antimunuFlux.dat'   ! anti-muon

           case default
           write(*,*) 'flavor and process IDs not compatible:2'
           stop

         end select
      case default
           write(*,*) 'flavor and process IDs not compatible:3'
     end select


    if (firsttime) then
        call read_fluxfile(NDIM,fluxfilename,jmax,enu,flux,sumflux)
        firsttime=.false.
    end if

!   Now call of sampling routine

      T2Konaxisenergy = eneut(NDIM,jmax,sumflux,enu)

 end function T2Konaxisenergy


 end module expNeutrinofluxes
