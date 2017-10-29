!******************************************************************************
!****m* /esample
! NAME
! module esample
!
! PURPOSE
! This module contains 2 files to read in the file with flux data for a given
! neutrino exeriment and for the sampling of neutrino energies from that flux
!
!******************************************************************************

module esample
  implicit none
  private

  public :: read_fluxfile
  public :: eneut

  !****************************************************************************
  !****g* esample/Enu_upper_cut
  ! SOURCE
  real, save, public :: Enu_upper_cut=200.
  ! PURPOSE
  ! cut events with neutrino energy above Enu_upper_cut;
  ! for ANL experiment, for example, Enu_upper_cut=1.5 for ppi0 and npi+ final
  ! state, but 5.98 for ppi+
  !****************************************************************************

  !****************************************************************************
  !****g* esample/Enu_lower_cut
  ! SOURCE
  real, save, public :: Enu_lower_cut=0.
  ! PURPOSE
  ! cut events with neutrino energy below Enu_lower_cut;
  ! for ANL experiment, for example, Enu_lower_cut=0.5 for some analysis of
  ! ppi+
  !****************************************************************************

  !****************************************************************************
  !****g* esample/energylimit_for_Qsrec
  ! SOURCE
  logical, save, public :: energylimit_for_Qsrec = .False.
  ! PURPOSE
  ! switch for using the energylimits Enu_upper_cut and Enu_lower_cut
  ! for the Q^2 reconstruction;
  ! values: .true. or .false.  (default: .false.)
  !****************************************************************************

contains

  !****************************************************************************
  !****f* esample/read_fluxfile
  ! NAME
  ! subroutine read_fluxfile(NDIM,fluxfilename,n_E,enu,flux,sumflux)
  !
  ! PURPOSE
  ! Reads the input file to obtain the energies and corresponding flux. The
  ! fluxfile must contain line by line: energy (middle of bin) and flux.
  ! The file can contain comments in the first few lines, starting with #
  ! the subroutine returns the number of elements in the fluxfile, as well as
  ! the cumulative flux distribution (sum of fluxes).
  !
  ! INPUTS
  ! * fluxfilename: name of input flux file, with energy and flux pairwise
  !   in different lines.
  ! * NDIM: maximal dimension of flux file
  !
  ! OUTPUT
  ! * n_E: number of records in data file (not counting comments).
  ! * enu: vector of neutrino energies.
  ! * flux: vector of flux values.
  ! * sumflux: vector with partial sums of flux (cumulative flux).
  !****************************************************************************
  subroutine read_fluxfile(NDIM,fluxfilename,n_E,enu,flux,sumflux)

    use inputGeneral, only: path_To_Input
    use output, only: write_ReadingInput
    use callstack, only: traceback
    implicit none

    integer :: j = 1, NDIM
    integer, intent(out) :: n_E
    integer :: ios
    integer :: ierr

    real :: sumf,en,fl
    real, dimension (:), intent(out) :: enu,flux
    real, dimension (0:NDIM), intent(out) :: sumflux
    character(len=*),intent(in) :: fluxfileName
    character(100) :: fileName,line
    logical,parameter :: verbose = .true.


    NAMELIST /nl_fluxcuts/  Enu_lower_cut, Enu_upper_cut, energylimit_for_Qsrec
    ! upper and lower energy cuts for flux and for Q^2 reconstruction

    call Write_ReadingInput('nl_fluxcuts',0)
    rewind(5)
    read(5,nml=nl_fluxcuts,IOSTAT=ios)
    call Write_ReadingInput("nl_fluxcuts",0,ios)
    write(*,*) "Enu_true_lower_cut = ",Enu_lower_cut
    write(*,*) "Enu_true_upper_cut =",Enu_upper_cut
    write(*,*) "E-Limit for Qs-Reconstruction =",energylimit_for_Qsrec

    !   Now reading of flux file from buuinput/neutrinos
    fileName=trim(path_to_Input)//'/neutrino/'//trim(fluxfilename)
    call Write_ReadingInput(trim(fileName),0)

    open(13,file=filename,status='old',action='read',iostat=ierr)
    if (ierr/=0) call traceback('ERROR: can not open file')

    ReadComments: do
       read(13,'(A)') line
       if (line(1:1) /= "#") exit ReadComments
       if (verbose) write(*,*) line
    end do ReadComments
    backspace(13)

    do
       read(13,*,iostat=ierr) en,fl
       if (ierr/=0) exit
       if (verbose) write(*,*) en,"    ",fl
       enu(j) = en
       flux(j) = fl
       j = j + 1
    end do

   if (ierr>0) call traceback('ERROR: reading flux file')

    close(13)
    call Write_ReadingInput(trim(fileName),1)

    !=== flux reading finished

    n_E = j-1  ! number of flux values in flux file

    !=== Now computation of cumulative flux distribution function

    sumflux = 0

    do j = 1,n_E
       if (enu(j) < Enu_lower_cut .or. enu(j) > Enu_upper_cut)  flux(j) = 0.0
       sumflux(j) = sumflux(j-1) + flux(j)
   !    write (*,*) sumflux(j)
    end do

    !=== now sumflux = normalized cumulative flux
    sumf = sum(flux)
    if (abs(sumflux(n_E) - sumf) >= 0.001*sumf) then
      write(*,*) 'problem with sumflux'
      stop
    end if

    write(*,*) n_E,sumflux(n_E), sumf

    sumflux = sumflux/sumf

    return

  end subroutine read_fluxfile


  !****************************************************************************
  !****f* esample/eneut
  ! NAME
  ! function eneut(NDIM,n_E,sumflux,enu)
  ! PURPOSE
  ! returns one energy value by sampling the flux distribution using discrete
  ! cumulative sampling
  !
  ! INPUTS
  ! * n_E: number of elements in input flux file and in sumflux and enu.
  ! * sumflux: vector of cumulative flux distributions.
  ! * enu: vector of neutrino energies.
  ! RESULT
  ! * eneut: one sampled neutrino energy
  !****************************************************************************
  real function eneut(NDIM,n_E,sumflux,enu)

    use random, only: rn
    implicit none

    integer, intent(in) :: n_E,NDIM
    integer :: j,l
    real, dimension (:), intent(in) :: enu
    real, dimension (0:NDIM), intent(in) :: sumflux
    real :: v,bin

    !   Now sampling by inversion of cumulative discrete cumulative flux

    v=rn()               !random number between 0 and 1

    ! linear search:
    do j = 1, n_E
       if (sumflux(j) >= v) exit
    end do
    l = j

    if (l==1) then
       bin = enu(l+1) - enu(l)        !bin = full energy bin width
    else
       bin = enu(l) - enu(l-1)
    end if

    eneut = enu(l) - 0.5*bin     &
         + bin* (v - sumflux(l-1))/(sumflux(l) - sumflux(l-1))
    ! Test printout
    !   write (*,*) "j=  ",j,"   eneut=",eneut
    ! End of testprintout
  end function eneut

end module esample
