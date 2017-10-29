program ThermalModelFit

  use, intrinsic :: iso_c_binding
  use fgsl
  use inputGeneral
  use particleProperties
  use IdTable, only: pion, nmes
  use CallStack, only: Traceback

  implicit none

  type tData 
     integer(fgsl_size_t) :: n
     integer(fgsl_int), allocatable :: id(:)
     real(fgsl_double), allocatable :: y(:)
     real(fgsl_double), allocatable :: sigma(:)
  end type tData

  type(tData), target :: fitData
  
  real, parameter :: massMax = 5.0 ! upper boundary of integration
  real :: arrIntMes(0:4, pion:pion+nmes-1) = 0.0
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

  integer :: modus = 0 ! 0: do not run fitter, 1: run fitter2, 2: run solver3
  

  
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
      call RunFitter2
     write(*,*) 'eps,rhoPi,rhoS = ',eps,rhoPi,rhoS
  case (2)
!     call RunFitter3
     write(*,*) 'eps,rhoPi,rhoS = ',eps,rhoPi,rhoS
     
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
         & dens, &
         & modus

    real :: startT = 0.1, startMuPi = 0.1, startMuS = 0.1
    real :: finalT = 0.2, finalMuPi = 0.0, finalMuS = 0.0

    real :: dens(pion:pion+nmes-1) = -1.0

    integer :: ios
    integer :: id,iData,nData

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
       muS = 0.0d0
!!$    case (2)
!!$       T = startT
!!$       muPi = startMuPi
!!$       muS = 0.0


    case DEFAULT
       call TRACEBACK('wrong modus')
    end select

    nData = 0
    do id=pion,pion+nmes-1
       if (dens(id) > 0.0) nData = nData+1
    end do

    if (nData < 2) then
       write(*,*) 'nData = ',nData
       call TRACEBACK('not enough data')
    end if

    allocate( fitData%id(nData), fitData%y(nData), fitData%sigma(nData) )
    fitData%n = nData
    
    iData = 0
    do id=pion,pion+nmes-1
       if (dens(id) > 0.0) then
          iData = iData + 1
          fitData%id(iData) = id
          fitData%y(iData) = dens(id)
          fitData%sigma(iData) = 1.0 ! dummy

          write(*,'(A,i3,A,1P,e12.3)') 'dens(',id,') = ',dens(id)
       end if
    end do
    
    call Write_ReadingInput('ThermalModel',1)
    
  end subroutine readInput

  !************************************************************************
  !****f* ThermalModelFit/func_f
  ! NAME
  ! integer(c_int) function func_f(x, cdata, f)
  ! PURPOSE
  ! The function to minimize
  ! INPUTS
  ! * type(c_ptr) :: x --- C-pointer to parameters
  ! * type(c_ptr) :: cdata --- C-pointer to data points
  ! * type(c_ptr) :: f --- C-pointer to calculated function values
  !************************************************************************
  integer(c_int) function func_f(x, cdata, f) bind(c)

    type(c_ptr), value :: x, cdata, f

    type(fgsl_vector) :: f_x, f_f
    type(tData), pointer :: f_data
    integer(fgsl_size_t) :: i
    integer(fgsl_int) :: status
    real(fgsl_double), pointer :: p_x(:), p_f(:)

    call fgsl_obj_c_ptr(f_x, x)
    call fgsl_obj_c_ptr(f_f, f)
    status = fgsl_vector_align(p_x, f_x)
    status = fgsl_vector_align(p_f, f_f)
    call c_f_pointer(cdata, f_data)

    T = p_x(1)
    muPi = p_x(2)
    !    muS = p_x(3)
    muS = 0.0d0
    
    call IntegrateMesons() ! T,muPi,muS --> eps,rhoPi,rhoS

    do i=1,f_data%n
       p_f(i) = (arrIntMes(1,f_data%id(i)) - f_data%y(i))/f_data%sigma(i)
    end do
    func_f = fgsl_success
  end function func_f

  !************************************************************************
  !****f* ThermalModelFit/func_df
  ! NAME
  ! integer(c_int) function func_df(x, cdata, j)
  ! PURPOSE
  ! The Jacobian of the function to minimize
  ! INPUTS
  ! * type(c_ptr) :: x --- C-pointer to parameters
  ! * type(c_ptr) :: cdata --- C-pointer to data points
  ! * type(c_ptr) :: j --- C-pointer to calculated jacobian values
  !************************************************************************
  integer(c_int) function func_df(x, cdata, j) bind(c)

    type(c_ptr), value :: x, cdata, j

    type(fgsl_vector) :: f_x
    type(fgsl_matrix) :: f_j
    type(tData), pointer :: f_data
    integer(fgsl_size_t) :: i
    integer(fgsl_int) :: status
    real(fgsl_double), pointer :: p_x(:), p_j(:,:)

    call fgsl_obj_c_ptr(f_x, x)
    call fgsl_obj_c_ptr(f_j, j)
    status = fgsl_vector_align(p_x, f_x)
    status = fgsl_matrix_align(p_j, f_j)
    call c_f_pointer(cdata, f_data)

    T = p_x(1)
    muPi = p_x(2)
    !    muS = p_x(3)
    muS = 0.0d0
    
    call IntegrateMesons() ! T,muPi,muS --> eps,rhoPi,rhoS

    do i=1,f_data%n
       p_j(1,i) = arrIntMes(2,f_data%id(i))/f_data%sigma(i) ! dN/dT
       p_j(2,i) = arrIntMes(3,f_data%id(i))/f_data%sigma(i) ! dN/dmuPi
!       p_j(3,i) = arrIntMes(4,f_data%id(i))/f_data%sigma(i) ! dN/dmuS
    end do
    func_df = fgsl_success
  end function func_df

  !************************************************************************
  !****f* ThermalModelFit/func_fdf
  ! NAME
  ! integer(c_int) function func_fdf(x, cdata, f, j)
  ! PURPOSE
  ! The function to minimize
  ! INPUTS
  ! * type(c_ptr) :: x --- C-pointer to parameters
  ! * type(c_ptr) :: cdata --- C-pointer to data points
  ! * type(c_ptr) :: f --- C-pointer to calculated function values
  ! * type(c_ptr) :: j --- C-pointer to calculated jacobian values
  !************************************************************************
  integer(c_int) function func_fdf(x, cdata, f, j) bind(c)

    type(c_ptr), value :: x, cdata, f, j
    
    integer(c_int) :: status
    status = func_f(x, cdata, f)
    status = func_df(x, cdata, j)
    func_fdf = fgsl_success
  end function func_fdf

  
  !************************************************************************
  !****s* ThermalModel/runFitter2
  ! NAME
  ! subroutine runFitter2()
  ! PURPOSE
  ! The routine to run a 2 dimensional fitting algorithm: given some
  ! dens(id) values, it finds (T,muPi)
  !************************************************************************
  subroutine runFitter2()

    type(c_ptr) :: pFitData
    integer(fgsl_size_t),parameter :: nParam = 2

    type(fgsl_multifit_function_fdf) :: nlfit_fdf
    type(fgsl_multifit_fdfsolver) :: nlfit_slv

    real(fgsl_double), target :: v_params(nParam)
    real(fgsl_double), target :: v_cov(nParam,nParam)
    real(fgsl_double), pointer :: v_fun(:), v_parptr(:)
    type(fgsl_vector) :: params
    type(fgsl_matrix) :: cov, jac
    real(fgsl_double) :: chi, c, dof

    integer :: iIt, i
    integer, parameter :: nIt = 50
    real(fgsl_double), parameter :: eps6 = 1.0d-10
    integer(fgsl_int) :: status
    
    pFitData = c_loc(fitData)

    nlfit_fdf = fgsl_multifit_function_fdf_init(func_f, func_df, func_fdf, &
         & fitData%n, nParam, pFitData)
    nlfit_slv = fgsl_multifit_fdfsolver_alloc(fgsl_multifit_fdfsolver_lmsder, &
         & fitData%n, nParam)

    v_params = (/ T, muPi /) ! start values

    ! create vector of values from array:
    params = fgsl_vector_init(1.0_fgsl_double)
    status = fgsl_vector_align(v_params,nParam,params,nParam,&
         0_fgsl_size_t,1_fgsl_size_t)

    ! setup fitter:
    status = fgsl_multifit_fdfsolver_set(nlfit_slv, nlfit_fdf, params)

    ! alignment, only needed once:
    status = fgsl_vector_align(v_fun, &
         fgsl_multifit_fdfsolver_f(nlfit_slv))
    status = fgsl_vector_align(v_parptr, &
         fgsl_multifit_fdfsolver_position(nlfit_slv))
    
    ! storage for cov within Fortran:
    cov = fgsl_matrix_init(type = 1.0_fgsl_double)
    status = fgsl_matrix_align(v_cov,nParam,nParam,nParam,cov)
    
    iIt = 0
    do
       iIt = iIt + 1
       status = fgsl_multifit_fdfsolver_iterate(nlfit_slv)
       if (status /= fgsl_success .or. iIt > nIt) exit

       write(*,*) iIt,v_parptr(:),sqrt(dot_product(v_fun,v_fun))
       
       status = fgsl_multifit_test_delta( &
            fgsl_multifit_fdfsolver_dx(nlfit_slv), &
            fgsl_multifit_fdfsolver_position(nlfit_slv), &
            eps6, eps6)
       if (status == fgsl_success) exit
    end do


    call WriteStatus(status)    

    jac = fgsl_multifit_fdfsolver_jac(nlfit_slv)
    status = fgsl_multifit_covar(jac, 0.0_fgsl_double, cov)
    chi = sqrt(dot_product(v_fun,v_fun))

    !  c only relevant for error estimate if chi*chi/dof exceeds 1
    dof = fitData%n - nParam
    if (fitData%n > nParam) then
       c = max(1.0d0, chi/sqrt(dof))
    else
       c = 1.0d0
    end if

    write(*,*) 'nIt = ',iIt
    write(6, '(''    T      = '',1PE15.8,'' +/- '',1PE12.5)') v_parptr(1), &
         c * v_cov(1,1)
    write(6, '(''    muPi   = '',1PE15.8,'' +/- '',1PE12.5)') v_parptr(2), &
         c * v_cov(2,2)

    write(6, '(''    chisq:                      '',1PE10.3)') chi*chi
    if (fitData%n > nParam) then
       write(6, '(''    chisq / degrees_of_freedom: '',1PE10.3)') chi*chi/dof
    end if
    write(*,*)
    
    call fgsl_vector_free(params)
    call fgsl_multifit_fdfsolver_free(nlfit_slv)
    call fgsl_multifit_function_fdf_free(nlfit_fdf)
    
  end subroutine runFitter2


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
  ! real function BoltzmannN(mass,nPi,nS, df)
  ! PURPOSE
  ! uses the global values (T,muPi,muS) in order to calculate the particle
  ! density for given parameters
  ! INPUTS
  ! * real :: mass
  ! * real :: nPi
  ! * real :: nS
  ! OUTPUT
  ! * function value
  ! * real, dimension(1:3), OPTIONAL:: holds dN/dT, dN/dmuPi, dN/dmuS
  !************************************************************************
  real function BoltzmannN(mass,nPi,nS, df)

    use constants, only: pi
    use fgsl, only: kn => fgsl_sf_bessel_kcn
    
    real, intent(in) :: mass,nPi,nS
    real, dimension(1:3), intent(out), optional :: df
    ! global: T, muPi, muS

    real :: betamu, betam, res
    real, parameter :: fak = 4*pi/(2*pi*0.197)**3

    betamu = (nS*muS + nPi*muPi)/T
    betam  = mass/T

    res = exp(betamu) * kn(2,betam)
    BoltzmannN = fak*res * mass**2 * T

    if (present(df)) then
       df(1) = (3.0 - betamu + betam*kn(1,betam)/kn(2,betam)) / T
       df(2) = nPi / T
       df(3) = nS / T

       df = df * BoltzmannN
    end if
    
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
    real, dimension(0:4) :: sum
    real, parameter :: dy0 = pi/200.

    real, dimension(1:3) :: partial
    
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
             arrIntMes(1,id) = BoltzmannN(mass0, nPi, nS, partial)

             arrIntMes(2:4,id) = partial(1:3)
             
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
          sum(1) = sum(1) + spectral/intfac * BoltzmannN(m, nPi, nS, partial)

          sum(2:4) = sum(2:4) + spectral/intfac * partial(1:3)
       end do
       arrIntMes(:,id) = sum * dy * deg 
       
    end do

    
  end subroutine IntegrateMesons

  


end program ThermalModelFit
