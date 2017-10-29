program ThermalModel

  use, intrinsic :: iso_c_binding
  use inputGeneral
  use particleProperties
  use IdTable, only: pion, nmes
  use CallStack, only: Traceback

  implicit none
  
  real, parameter :: massMax = 5.0 ! upper boundary of integration
  real :: arrIntMes(0:2, pion:pion+nmes-1) = 0.0
  logical :: arrFlagMes(pion:pion+nmes-1) = .true.
  real :: arrNpiMes(pion:pion+nmes-1) = 0.0

  ! final values (after root finding):
  real :: T = 0.200  ! temperature
  real :: muPi = 0.0 ! pion chemical potential
  real :: muS = 0.0  ! strangeness chemical potential

  ! input values
  real :: eps = 0.519 ! energy density
  real :: rhoPi = 0.1 ! pion density
  real :: rhoS = 0.0  ! strangeness density

  real :: densPi = 0.0 ! pion density (fm^-3)
  real :: densRho= 0.0 ! rho density (fm^-3)
  
  integer :: modus = 0 ! 0: do not run solver, 1: run solver2, 2: run solver3
  

  
  arrNpiMes(101:122) = (/1.0, 3.0, 2.0, 2.0, 2.98, 5.0, 2.16, 0.0, 0.0, &
       & 1.0, 1.0, 2.0, 2.0 , &
       & 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.4/)

  call readInputGeneral
  call initParticleProperties
  call readInput
  call calcNpi
  write(*,*) '=== arrNpiMes: ==='
  write(*,'(50f7.3)') arrNpiMes
  write(*,*) '=================='
  write(*,*)
  
  select case (modus)
  case (0)
  case (1)
     arrFlagMes = .false.
     arrFlagMes(101) = .true.
     arrFlagMes(103) = .true.
     call RunSolver2
     write(*,*) 'eps,rhoPi,rhoS = ',eps,rhoPi,rhoS
  case (2)
     call RunSolver3
     write(*,*) 'eps,rhoPi,rhoS = ',eps,rhoPi,rhoS
  case (3)
     arrFlagMes = .false.
     arrFlagMes(101) = .true.
     arrFlagMes(103) = .true.
     write(*,*) 'start: T,muPi,muS = ',T,muPi,muS
     call RunSolverMult2
     
  end select

  write(*,*) 'T,muPi,muS = ',T,muPi,muS

  call PrintFinalList()
  
contains

  !************************************************************************
  !****s* ThermalModel/readInput
  ! NAME
  ! subroutine readInput
  ! PURPOSE
  ! Reads input in jobcard out of namelist "Thermalmodel". 
  !************************************************************************
  subroutine readInput
    use output, only: Write_ReadingInput

    NAMELIST /ThermalModel/ eps,rhoPi,rhoS, &
         & startT,startMuPi,startMuS, &
         & finalT,finalMuPi,finalMuS, &
         & densPi, densRho, &
         & modus

    real :: startT = 0.1, startMuPi = 0.1, startMuS = 0.1
    real :: finalT = 0.2, finalMuPi = 0.0, finalMuS = 0.0 
    
    integer :: ios

    call Write_ReadingInput('ThermalModel',0)
    rewind(5)
    read(5,nml=ThermalModel,IOSTAT=ios)
    call Write_ReadingInput('ThermalModel',0,ios)

    write(*,*) 'modus : ', modus
    select case(modus)
    case (0)
       T = finalT
       muPi = finalMuPi
       muS = finalMuS
    case (1)
       T = startT
       muPi = startMuPi
       muS = 0.0
    case (2)
       T = startT
       muPi = startMuPi
       muS = startMuS
    case (3)
       T = startT
       muPi = startMuPi
       muS = 0.0

       write(*,*) 'n_pi  = ',densPi
       write(*,*) 'n_rho = ',densRho
       
       
    case DEFAULT
       call TRACEBACK('wrong modus')
    end select
       

    call Write_ReadingInput('ThermalModel',1)
    
  end subroutine readInput

  !************************************************************************
  !****f* ThermalModel/funcToMinimize2
  ! NAME
  ! integer(c_int) function funcToMinimize2(x, params, f)
  ! PURPOSE
  ! The function to minimize with solver2
  !************************************************************************
  function funcToMinimize2(x, params, f) bind(c)
    use fgsl, only: fgsl_success
    use fgsl, only: fgsl_int, fgsl_double, fgsl_vector
    use fgsl, only: fgsl_obj_c_ptr, fgsl_vector_align
    use particleProperties, only: hadron
    
    type(c_ptr), value :: x, params, f
    integer(c_int) :: funcToMinimize2

    type(fgsl_vector) :: fx, ff
    real(fgsl_double), pointer :: par(:), xv(:), yv(:)
    integer(fgsl_int) :: status
    real, dimension(0:3) :: total
    integer :: id
    
    call fgsl_obj_c_ptr(fx, x)
    call fgsl_obj_c_ptr(ff, f)
    call c_f_pointer(params, par, (/ 2 /))
    status = fgsl_vector_align(xv, fx)
    status = fgsl_vector_align(yv, ff)

    T = xv(1)
    muPi = xv(2)

    call IntegrateMesons() ! T,muPi,muS --> eps,rhoPi,rhoS
    total = 0
    do id=pion,pion+nmes-1
       total(1:2) = total(1:2) + arrIntMes(1:2,id)
       total(0)  = total(0)  + arrIntMes(1,id) * arrNpiMes(id) ! yes, it's 1
       total(3)  = total(3)  + arrIntMes(1,id) * hadron(ID)%strangeness
    end do
    eps = total(2)
    rhoPi = total(0)
    rhoS = total(3)
    
    yv(1) = eps   - par(1)
    yv(2) = rhoPi - par(2)
!    write(*,*) '>>',xv(1:2)
!    write(*,*) '::',yv(1:2)
    
    funcToMinimize2 = fgsl_success
  end function funcToMinimize2

  !************************************************************************
  !****f* ThermalModel/funcToMinimize3
  ! NAME
  ! integer(c_int) function funcToMinimize3(x, params, f)
  ! PURPOSE
  ! The function to minimize with solver3
  !************************************************************************
  function funcToMinimize3(x, params, f) bind(c)
    use fgsl, only: fgsl_success
    use fgsl, only: fgsl_int, fgsl_double, fgsl_vector
    use fgsl, only: fgsl_obj_c_ptr, fgsl_vector_align
    use particleProperties, only: hadron
    
    type(c_ptr), value :: x, params, f
    integer(c_int) :: funcToMinimize3

    type(fgsl_vector) :: fx, ff
    real(fgsl_double), pointer :: par(:), xv(:), yv(:)
    integer(fgsl_int) :: status
    real, dimension(0:3) :: total
    integer :: id
    
    call fgsl_obj_c_ptr(fx, x)
    call fgsl_obj_c_ptr(ff, f)
    call c_f_pointer(params, par, (/ 3 /))
    status = fgsl_vector_align(xv, fx)
    status = fgsl_vector_align(yv, ff)

    T = xv(1)
    muPi = xv(2)
    muS = xv(3)

    call IntegrateMesons() ! T,muPi,muS --> eps,rhoPi,rhoS
    total = 0
    do id=pion,pion+nmes-1
       total(1:2) = total(1:2) + arrIntMes(1:2,id)
       total(0)  = total(0)  + arrIntMes(1,id) * arrNpiMes(id) ! yes, it's 1
       total(3)  = total(3)  + arrIntMes(1,id) * hadron(ID)%strangeness
    end do
    eps = total(2)
    rhoPi = total(0)
    rhoS = total(3)
    
    yv(1) = eps   - par(1)
    yv(2) = rhoPi - par(2)
    yv(3) = rhoS  - par(3)
!    write(*,*) '>>',xv(1:3)
!    write(*,*) '::',yv(1:3)
    
    funcToMinimize3 = fgsl_success
  end function funcToMinimize3

  !************************************************************************
  !****f* ThermalModel/funcToMinimizeMult2
  ! NAME
  ! integer(c_int) function funcToMinimizeMult2(x, params, f)
  ! PURPOSE
  ! The function to minimize with solverMult2
  !************************************************************************
  function funcToMinimizeMult2(x, params, f) bind(c)
    use fgsl, only: fgsl_success
    use fgsl, only: fgsl_int, fgsl_double, fgsl_vector
    use fgsl, only: fgsl_obj_c_ptr, fgsl_vector_align
    
    type(c_ptr), value :: x, params, f
    integer(c_int) :: funcToMinimizeMult2

    type(fgsl_vector) :: fx, ff
    real(fgsl_double), pointer :: par(:), xv(:), yv(:)
    integer(fgsl_int) :: status
    
    call fgsl_obj_c_ptr(fx, x)
    call fgsl_obj_c_ptr(ff, f)
    call c_f_pointer(params, par, (/ 2 /))
    status = fgsl_vector_align(xv, fx)
    status = fgsl_vector_align(yv, ff)

    T = xv(1)
    muPi = xv(2)

    call IntegrateMesons() ! T,muPi,muS --> eps,rhoPi,rhoS
    
    yv(1) = arrIntMes(1,101) - par(1)
    yv(2) = arrIntMes(1,103) - par(2)
!!$    write(*,*) '>>',xv(1:2)
!!$    write(*,*) '::',yv(1:2)
!!$    write(*,*) '::::',arrIntMes(1,101),par(1)
!!$    write(*,*) '::::',arrIntMes(1,103),par(2)
    
    funcToMinimizeMult2 = fgsl_success
  end function funcToMinimizeMult2

  !************************************************************************
  !****s* ThermalModel/runSolver2
  ! NAME
  ! subroutine runSolver2()
  ! PURPOSE
  ! The routine to run a 2 dimensional root finding algorithm: given
  ! (eps,rhoPi), it finds (T,muPi)
  !************************************************************************
  subroutine runSolver2()

    use fgsl ! too much stuff for 'only'

    integer(fgsl_size_t),parameter :: nParam = 2
    type(fgsl_multiroot_fsolver) :: solver
    type(fgsl_multiroot_function) :: func
    real(fgsl_double), target :: param(nParam), xv(nParam)
    type(fgsl_vector) :: xvec, fvec
    real(fgsl_double), pointer :: fptr(:), xptr(:)
    type(c_ptr) :: pParam
    integer(fgsl_int) :: status
    integer, parameter :: nIt = 100
    integer :: iIt
    

    solver = fgsl_multiroot_fsolver_alloc(fgsl_multiroot_fsolver_hybrids,nParam)

    param = (/ eps, rhoPi /)
    pParam = c_loc(param)
    func = fgsl_multiroot_function_init(funcToMinimize2,nParam,pParam)

    xv = (/ T, muPi /) ! start values

    ! create vector of values from array:
    xvec = fgsl_vector_init(1.0_fgsl_double)
    status = fgsl_vector_align(xv,nParam,xvec,nParam,0_fgsl_size_t,1_fgsl_size_t)
    
    ! setup solver:
    status = fgsl_multiroot_fsolver_set(solver, func, xvec)

    fvec = fgsl_multiroot_fsolver_f(solver)
    status = fgsl_vector_align(fptr, fvec)
    
    call fgsl_vector_free(xvec)
    xvec = fgsl_multiroot_fsolver_root(solver)
    status = fgsl_vector_align(xptr, xvec)

    iIt = 0
    do
       iIt = iIt+1
       status = fgsl_multiroot_fsolver_iterate(solver);
       if (status /= fgsl_success .or. iIt > nIt) exit

       status = fgsl_multiroot_test_residual(fvec, 1.0d-10)
       if (status == fgsl_success) exit

       write(*,*) iIt,xptr(1:2),fptr(1:2)
    end do

    call WriteStatus(status)
    
  end subroutine runSolver2


  !************************************************************************
  !****s* ThermalModel/runSolver3
  ! NAME
  ! subroutine runSolver3()
  ! PURPOSE
  ! The routine to run a 3 dimensional root finding algorithm: given
  ! (eps,rhoPi,rhoS), it finds (T,muPi,rhoS)
  !************************************************************************
  subroutine runSolver3()

    use fgsl ! too much stuff for 'only'

    integer(fgsl_size_t),parameter :: nParam = 3
    type(fgsl_multiroot_fsolver) :: solver
    type(fgsl_multiroot_function) :: func
    real(fgsl_double), target :: param(nParam), xv(nParam)
    type(fgsl_vector) :: xvec, fvec
    real(fgsl_double), pointer :: fptr(:), xptr(:)
    type(c_ptr) :: pParam
    integer(fgsl_int) :: status
    integer, parameter :: nIt = 100
    integer :: iIt
    
    solver = fgsl_multiroot_fsolver_alloc(fgsl_multiroot_fsolver_hybrids,nParam)

    param = (/ eps, rhoPi, rhoS /)
    pParam = c_loc(param)
    func = fgsl_multiroot_function_init(funcToMinimize3,nParam,pParam)

    xv = (/ T, muPi, muS /) ! start values

    ! create vector of values from array:
    xvec = fgsl_vector_init(1.0_fgsl_double)
    status = fgsl_vector_align(xv,nParam,xvec,nParam,0_fgsl_size_t,1_fgsl_size_t)
    
    ! setup solver:
    status = fgsl_multiroot_fsolver_set(solver, func, xvec)

    fvec = fgsl_multiroot_fsolver_f(solver)
    status = fgsl_vector_align(fptr, fvec)
    
    call fgsl_vector_free(xvec)
    xvec = fgsl_multiroot_fsolver_root(solver)
    status = fgsl_vector_align(xptr, xvec)

    iIt = 0
    do
       iIt = iIt+1
       status = fgsl_multiroot_fsolver_iterate(solver);
       if (status /= fgsl_success .or. iIt > nIt) exit

       status = fgsl_multiroot_test_residual(fvec, 1.0d-10)
       if (status == fgsl_success) exit

       write(*,*) iIt,xptr(1:3),fptr(1:3)
    end do

    call WriteStatus(status)
    
  end subroutine runSolver3

  !************************************************************************
  !****s* ThermalModel/runSolverMult2
  ! NAME
  ! subroutine runSolverMult2()
  ! PURPOSE
  ! The routine to run a 2 dimensional root finding algorithm: given
  ! (densPi,densRho), it finds (T,muPi)
  !************************************************************************
  subroutine runSolverMult2()

    use fgsl ! too much stuff for 'only'

    integer(fgsl_size_t),parameter :: nParam = 2
    type(fgsl_multiroot_fsolver) :: solver
    type(fgsl_multiroot_function) :: func
    real(fgsl_double), target :: param(nParam), xv(nParam)
    type(fgsl_vector) :: xvec, fvec
    real(fgsl_double), pointer :: fptr(:), xptr(:)
    type(c_ptr) :: pParam
    integer(fgsl_int) :: status
    integer, parameter :: nIt = 100
    integer :: iIt
    

    solver = fgsl_multiroot_fsolver_alloc(fgsl_multiroot_fsolver_hybrids,nParam)

    param = (/ densPi, densRho /)
    pParam = c_loc(param)
    func = fgsl_multiroot_function_init(funcToMinimizeMult2,nParam,pParam)

    xv = (/ T, muPi /) ! start values

    ! create vector of values from array:
    xvec = fgsl_vector_init(1.0_fgsl_double)
    status = fgsl_vector_align(xv,nParam,xvec,nParam,0_fgsl_size_t,1_fgsl_size_t)
    
    ! setup solver:
    status = fgsl_multiroot_fsolver_set(solver, func, xvec)

    fvec = fgsl_multiroot_fsolver_f(solver)
    status = fgsl_vector_align(fptr, fvec)
    
    call fgsl_vector_free(xvec)
    xvec = fgsl_multiroot_fsolver_root(solver)
    status = fgsl_vector_align(xptr, xvec)

    iIt = 0
    do
       iIt = iIt+1
       status = fgsl_multiroot_fsolver_iterate(solver);
       if (status /= fgsl_success .or. iIt > nIt) exit

       status = fgsl_multiroot_test_residual(fvec, 1.0d-10)
       if (status == fgsl_success) exit

       write(*,*) iIt,xptr(1:2),fptr(1:2)
    end do
    
    call WriteStatus(status)
    
  end subroutine runSolverMult2

  !************************************************************************
  !****s* ThermalModel/WriteStatus
  ! NAME
  ! subroutine WriteStatus(status)
  ! PURPOSE
  ! Write a gsl status string to stdout
  !************************************************************************
  subroutine WriteStatus(status)
    use fgsl, only: fgsl_int, fgsl_char, fgsl_strmax
    use fgsl, only: fgsl_strerror
    
    integer(fgsl_int),intent(in) :: status
    character(kind=fgsl_char, len=fgsl_strmax) :: message
    message = fgsl_strerror(status)
    write(*,*)
    write(*,*) 'Status: ',message
    write(*,*)
    
  end subroutine WriteStatus

  !************************************************************************
  !****s* ThermalModel/PrintFinalList
  ! NAME
  ! subroutine PrintFinalList()
  ! PURPOSE
  ! Write a list of densities and energy densities of all particles
  !************************************************************************
  subroutine PrintFinalList()

    integer :: id
    real, dimension(0:2) :: total
    
    call IntegrateMesons()

    write(*,'(A6,4A13)') 'id:','norm:','n:','e:',''
    do id=pion,pion+nmes-1
       write(*,'(i6,1P,3e13.4)') id, arrIntMes(0:2,id)
       total(1:2) = total(1:2) + arrIntMes(1:2,id)
       total(0)  = total(0)  + arrIntMes(1,id) * arrNpiMes(id) ! yes, it's 1
    end do
    
    write(*,*)
    write(*,'(A6,4A13)') '','','n:','e:','n_pi:'
    write(*,'(A6,A13,1P,3e13.4)') 'total', '', total(1:2),total(0)
    
    total(1:2) = arrIntMes(1:2,101) + arrIntMes(1:2,103)
    total(0)  = arrIntMes(1,101) * arrNpiMes(101) + arrIntMes(1,103) * arrNpiMes(103)
    write(*,*)
    write(*,'(A6,A13,1P,3e13.4)') 'pi+rho', '', total(1:2),total(0)
    
    write(*,*)
    write(*,*) 'multiplicity ratio: pi/rho = ',arrIntMes(1,101)/arrIntMes(1,103)
    
  end subroutine PrintFinalList

  !************************************************************************
  !****f* ThermalModel/BoltzmannN
  ! NAME
  ! real function BoltzmannN(mass,nPi,nS)
  ! PURPOSE
  ! uses the global values (T,muPi,muS) in order to calculate the particle
  ! density for given parameters
  ! INPUTS
  ! * real :: mass
  ! * real :: nPi
  ! * real :: nS
  ! OUTPUT
  ! * function value
  !************************************************************************
  real function BoltzmannN(mass,nPi,nS)

    use constants, only: pi
    use fgsl, only: kn => fgsl_sf_bessel_kcn
    
    real, intent(in) :: mass,nPi,nS
    ! global: T, muPi, muS

    real :: betamu, betam, res
    real, parameter :: fak = 4*pi/(2*pi*0.197)**3

    betamu = (nS*muS + nPi*muPi)/T
    betam  = mass/T

    res = exp(betamu) * kn(2,betam)
    BoltzmannN = fak*res * mass**2 * T

    return
  end function BoltzmannN

  !************************************************************************
  !****f* ThermalModel/BoltzmannE
  ! NAME
  ! real function BoltzmannE(mass,nPi,nS)
  ! PURPOSE
  ! uses the global values (T,muPi,muS) in order to calculate the energy
  ! density for given parameters
  !************************************************************************
  real function BoltzmannE(mass,nPi,nS)

    use constants, only: pi
    use fgsl, only: kn => fgsl_sf_bessel_kcn
    
    real, intent(in) :: mass,nPi,nS
    ! global: T, muPi, muS

    real :: betamu, betam, res
    real, parameter :: fak = 4*pi/(2*pi*0.197)**3

    betamu = (nS*muS + nPi*muPi)/T
    betam  = mass/T

    res = exp(betamu) * (kn(2,betam)+betam/3.*kn(1,betam))
    BoltzmannE = fak * res * 3 * mass**2 * T**2
    return
  end function BoltzmannE

  !************************************************************************
  !****s* ThermalModel/CalcNpi
  ! NAME
  ! subroutine CalcNpi
  ! PURPOSE
  ! fill the array holding the pion equivalent numbers
  !************************************************************************
  subroutine CalcNpi()
    
    use constants, only: pi
    use particleProperties, only: hadron, nDecays, isCharmed
    use mesonWidth, only: fullWidthMeson,decayWidthMeson
    
    integer :: id, im, nm, iid
    real :: mass0, gamma0
    real :: mmin, mmax, m
    real :: ymin, ymax, dy, y
    real :: gamma, spectral, intfac
    real, dimension(0:2) :: ss
    real, parameter :: dy0 = pi/100.
    real, dimension(nDecays)  :: decayWidth
    real :: sDecayWidth, nnpi
    
    arrNpiMes = 0
    
    do id=pion,pion+nmes-1
       if (.not.arrFlagMes(id)) cycle ! ignore this particle

!       write(*,*) '======== ID = ',ID
!       write(*,*) '#',hadron(ID)%decaysID
!       write(*,*)

       if (id==108) cycle ! skip eta_c
       if (isCharmed(id)) cycle
       
       mass0  = hadron(id)%mass
       gamma0 = hadron(id)%width
       
       if (gamma0 < 1e-10) then
          if (massMax > mass0) then
             select case(ID)
             case (pion)
                arrNpiMes(ID) = 1.0
             case (110,111) ! K, barK
                arrNpiMes(ID) = 1.0
             end select

          end if
          cycle
       end if

       mmin = hadron(id)%minmass
       mmax = massMax
       if (mmax < mmin) cycle ! integral is zero

       ymax = 2*atan(2*(mmax-mass0) / gamma0)
       ymin = 2*atan(2*(mmin-mass0) / gamma0)

       nm = max(int((ymax-ymin)/dy0),1)
       dy  = (ymax-ymin)/float(nm)

       ss = 0.0
       do im=1,nm
          y = ymin+(float(im)-0.5)*dy
          m = 0.5*tan(0.5*y)*gamma0 + mass0
          m = min(max(m,mmin),mmax)
          gamma = fullWidthMeson(id, m)
          spectral = 2./pi * m**2 * gamma / ((m**2-mass0**2)**2+m**2*gamma**2)
          intfac = gamma0 / ((m-mass0)**2+gamma0**2/4.)

          decayWidth = decayWidthMeson(id,m,0)
          sDecayWidth = sum(decayWidth)
!          write(*,*) m,decayWidth/sDecayWidth

          nnpi = 0
          do iid=1,nDecays
             if (decayWidth(iid) > 1e-10) then
                select case(hadron(ID)%decaysID(iid))
                case (7)
                   ! no pion
                case (6)
                   nnpi = nnpi + 1 * decayWidth(iid)
                case (1,5)
                   nnpi = nnpi + 2 * decayWidth(iid)
                case (2,-2,-3)
                   nnpi = nnpi + 3 * decayWidth(iid)
                case (-1,-4)
                   nnpi = nnpi + (2+arrNpiMes(102)) * decayWidth(iid)
                case (3) ! K barK as 2 pions
                   nnpi = nnpi + 2 * decayWidth(iid)
                case (4,12) ! pi K  as 2 pions
                   nnpi = nnpi + 2 * decayWidth(iid)
                case default
                   write (*,*) 'unknown:',hadron(ID)%decaysID(iid)
                   stop
                end select
             end if
          end do
          
          ss(0) = ss(0) + spectral/intfac
          ss(1) = ss(1) + spectral/intfac * sDecayWidth
          ss(2) = ss(2) + spectral/intfac * nnpi
       end do
       arrNpiMes(ID) = ss(2)/ss(1)

!       write(*,*) '::::',ss,ss(2)/ss(1)
!       write(*,*)
!       write(*,*)
       
    end do

  end subroutine CalcNpi

  !************************************************************************
  !****s* ThermalModel/IntegrateMesons
  ! NAME
  ! subroutine IntegrateMesons()
  ! PURPOSE
  ! uses the global values (T,muPi,muS) in order to calculate the particle
  ! and energy densities of all particles.
  !************************************************************************
  subroutine IntegrateMesons()

    use constants, only: pi
    use particleProperties, only: hadron
    use mesonWidth, only: fullWidthMeson
    
    integer :: id, im, nm
    real :: mass0, gamma0, nS, nPi
    real :: mmin, mmax, m
    real :: ymin, ymax, dy, y
    real :: gamma, spectral, intfac, deg
    real, dimension(0:2) :: sum
    real, parameter :: dy0 = pi/200.

    arrIntMes = 0.0
    
    do id=pion,pion+nmes-1
       if (.not.arrFlagMes(id)) cycle ! ignore this particle

       mass0  = hadron(id)%mass
       gamma0 = hadron(id)%width
       nS = hadron(ID)%strangeness
       nPi = arrNpiMes(id)
       deg = (hadron(ID)%Spin*2+1)*(hadron(ID)%isoSpinTimes2+1)
       
       if (gamma0 < 1e-3) then
          if (massMax > mass0) then
             arrIntMes(0,id) = 1.0
             arrIntMes(1,id) = BoltzmannN(mass0, nPi, nS)
             arrIntMes(2,id) = BoltzmannE(mass0, nPi, nS)

             arrIntMes(:,id) = arrIntMes(:,id) * deg
          end if
          cycle
       end if

       mmin = hadron(id)%minmass
       mmax = massMax
       if (mmax < mmin) cycle ! integral is zero

       ymax = 2*atan(2*(mmax-mass0) / gamma0)
       ymin = 2*atan(2*(mmin-mass0) / gamma0)

       nm = max(int((ymax-ymin)/dy0),1)
       dy  = (ymax-ymin)/float(nm)

       sum = 0.0
       do im=1,nm
          y = ymin+(float(im)-0.5)*dy
          m = 0.5*tan(0.5*y)*gamma0 + mass0
          m = min(max(m,mmin),mmax)
          gamma = fullWidthMeson(id, m)
          spectral = 2./pi * m**2 * gamma / ((m**2-mass0**2)**2+m**2*gamma**2)
          intfac = gamma0 / ((m-mass0)**2+gamma0**2/4.)
          sum(0) = sum(0) + spectral/intfac
          sum(1) = sum(1) + spectral/intfac * BoltzmannN(m, nPi, nS)
          sum(2) = sum(2) + spectral/intfac * BoltzmannE(m, nPi, nS)
       end do
       arrIntMes(:,id) = sum * dy * deg 
       
    end do

    
  end subroutine IntegrateMesons

  


end program ThermalModel
