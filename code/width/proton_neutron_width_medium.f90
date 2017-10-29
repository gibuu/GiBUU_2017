!******************************************************************************
!****m* /pn_medium_width
! NAME
! module pn_medium_width
! PURPOSE
! Gives the in-medium width of a proton or a neutron.
! AUTHOR
! David F. Kalok
!******************************************************************************
module pn_medium_width
  use inputGeneral, only:  path_to_Input
  use constants

  implicit none
  private

  !****************************************************************************
  !****n* pn_medium_width/pn_medium
  ! NAME
  ! NAMELIST /pn_medium/
  ! PURPOSE
  ! Includes the input switches:
  ! * density_dependent
  ! * pn_medium_switch
  ! * form_factor
  !****************************************************************************

  !****************************************************************************
  !****ig*  pn_medium_width/debug
  ! SOURCE
  !
  logical, parameter :: debug=.false.
  ! PURPOSE
  ! If .true. debug infos are produced, if .false. not..
  !****************************************************************************

  !****************************************************************************
  !****ig*  pn_medium_width/density_dependent
  ! SOURCE
  !
  logical, save :: density_dependent=.false.
  ! PURPOSE
  ! the density of the spectral function
  !****************************************************************************

  !****************************************************************************
  !****ig*  pn_medium_width/form_factor
  ! SOURCE
  !
  logical, save :: form_factor=.true.
  ! PURPOSE
  ! If .true. the form factor for the width is used
  !****************************************************************************

  !****************************************************************************
  !****ig*  pn_medium_width/pn_medium_switch
  ! SOURCE
  !
  logical, save :: pn_medium_switch=.true.
  ! PURPOSE
  ! If .true. medium_modifications will be used
  !****************************************************************************

  logical , save :: initFlag=.true.


  public :: proton_width_medium

contains


  !****************************************************************************
  !****s* pn_medium_width/readinput
  ! NAME
  ! subroutine readinput
  ! PURPOSE
  ! This subroutine reads out the jobcard "pn_medium" which may contain:
  ! * density_dependent
  ! * form_factor
  ! * pn_medium_switch
  ! Only called once to initialize the module.
  !****************************************************************************
  subroutine readinput

    use output
    integer :: ios

    NAMELIST /pn_medium/ pn_medium_switch, density_dependent, form_factor

    call Write_ReadingInput("pn_medium",0)
    rewind(5)
    read(5,nml=pn_medium,IOSTAT=IOS)
    call Write_ReadingInput("pn_medium",0,IOS)

    write(*,*) 'Set pn_medium_switch       :' ,pn_medium_switch
    call Write_ReadingInput('pn_medium',1)

  end subroutine readinput

  !****************************************************************************
  !****f* pn_medium_width/proton_width_medium
  ! NAME
  ! real function proton_width_medium(momLRF,pdens,ndens,mu_extern)
  ! PURPOSE
  ! function for calculation of  the proton width in a medium
  ! INPUTS
  ! * real, dimension(0:3) :: momLRF -- momentum of resonance in LRF
  ! * real :: pdens -- proton density [fm^{-3}]
  ! * real :: ndens -- neutron density [fm^{-3}]
  ! * real, OPTIONAL :: mu_extern -- needed for density dependent width
  ! AUTHOR
  ! David F. Kalok
  !****************************************************************************
  real function proton_width_medium(momLRF,pdens,ndens,mu_extern)
!    use init_pn_medium
    !   use spline2d
    use constants, only: hbarc, mN
    use baryonPotentialModule, only: BaryonPotential
    use particleDefinition
    use mediumdefinition

    real,intent(in),dimension(0:3) :: momLRF
    real, intent(in) :: pdens
    real, intent(in) :: ndens
    real, optional :: mu_extern
    real:: rhoBaryon
    type(particle), save :: teilchen_fermi
    real:: omega_F,k_fermi
    type(medium) :: med


    if (initflag) then
       write(*,*) 'proton_width_medium initialized'
       call readinput

       call setToDefault(teilchen_fermi)
       teilchen_fermi%mass=mN
       ! finalstate%charge =charge_out
       ! teilchen_fermi%momentum=p_out
       teilchen_fermi%position=(/0.,0.,0./)
       teilchen_fermi%antiparticle=.false.
       teilchen_fermi%id=1
       teilchen_fermi%perturbative=.false.
       ! finalState%productionTime=0.
       ! finalState%lastCollisionTime=0.
       ! finalState%formationTime=0.
       ! finalState%scaleCS=1.
       ! finalState%in_Formation=.false.

       initflag=.false.
    end if

    if (present(mu_extern)) then
       if (density_dependent) then
          proton_width_medium=proton_width_mdd(momLRF,pdens,ndens,mu_extern)
       else
          proton_width_medium=proton_width_mdd(momLRF,pdens,ndens,mu_extern)
       end if
    else

       rhoBaryon=pdens+ndens
       k_fermi=(3./2.* pi**2.*rhoBaryon)**(1./3.)*hbarc

       med%temperature    =0.
       med%useMedium      =.true.
       med%density        = pdens+ndens
       med%densityProton  = pdens
       med%densityNeutron = ndens

       teilchen_fermi%momentum(0)=sqrt((teilchen_fermi%mass+&
            & BaryonPotential(teilchen_fermi,med,.true.))**2 +k_fermi**2)
       omega_F=teilchen_fermi%momentum(0)
       proton_width_medium=proton_width_mdd(momLRF,pdens,ndens,omega_F)
    end if

  end function proton_width_medium

  !****************************************************************************
  !****f* pn_medium_width/proton_width_mdd
  ! NAME
  ! real function proton_width_mdd(momLRF,pdens,ndens,mu_extern)
  ! PURPOSE
  ! function for calculation of  the proton width in a medium.
  ! Reads out data file and uses a spline with it.
  ! INPUTS
  ! * real, dimension(0:3) :: momLRF -- momentum of resonance in LRF
  ! * real :: pdens -- proton density [fm^{-3}]
  ! * real :: ndens -- neutron density [fm^{-3}]
  ! * real, OPTIONAL :: mu_external -- needed for density dependent width
  ! NOTES
  ! It uses the data input files /baryon/0*.34
  ! right now only for z/a =0.5
  !****************************************************************************
  real function proton_width_mdd(momLRF,pdens,ndens,mu_extern)

    real,intent(in),dimension(0:3) :: momLRF ! momentum of resonance in LRF
    real, intent(in) :: pdens
    real, intent(in) :: ndens
    real, optional :: mu_extern
    integer, parameter :: Dateiende=-1

    !   real, dimension(1:200),save :: pLabField, sigmaField, derivativeField

    logical, save :: initFlag=.true.
    integer :: ios
    real,save :: maximalEnergy,maximalRho,maxMomentum!,maximalWidth
    real, save ::  minEnergy,minimalRho,minMomentum!,minWidth
    integer :: i,j,k,l,n_count,n_rho,n_momentum,n_energy
    integer, save :: i_max
    REAL :: temp1, temp2, temp3, temp4, tmp!,temp5,temp6
    REAL :: denstmp,omegatmp, momentumtmp,gammatmp
    !REAL, save :: mass_baryon
    REAl, save :: drho, domega, dp
    REAL, ALLOCATABLE,  DIMENSION(:) :: proton_omega, proton_momentum, proton_gamma, proton_spectral
    REAL, DIMENSION(32), save :: rho_field
    REAL, DIMENSION(121), save ::momentum_field, omega_field
    REAL, DIMENSION(32,121,121),save :: gamma_matrix  !,  gam_deriv_mat
    logical, save :: allocat=.true.
    logical, save :: fill_field=.true.
    logical :: lexist
    !real :: p_fermi_extern, p_fermi_intern

1002 format(4(e12.4,1x))

    init:  if (initFlag) then

       write(*,*) 'density dependent width for protons and neutrons  is used'
       call readinput
       n_count=1
       if (debug) then
          write(*,*) 'Initializing proton_width_medium_density_dependent'
       end if
       densityloop: do l=1,32
          ! write(*,*) 'n_count= ', n_count
          rho_field(n_count)=float(n_count)*0.01
          open(100,file=trim(path_to_Input)//'/baryon/0'//nummern(n_count)//'.34',status='old')
          i=1
          INQUIRE(FILE=trim(path_to_Input)//'/baryon/0'//nummern(n_count)//'.34',EXIST=lexist)
          if (.NOT. lexist) then
             write(*,*) 'Error in subroutine proton_width_medium_density_dependent, no file: ', &
                  & trim(path_to_Input)//'/baryon/0'//nummern(n_count)//'.34', ' exists!!!'
             stop
          end if
          allo: if (allocat) then
             length: do

                read(100,1002,IOSTAT=ios)  temp1, temp2, temp3, temp4 !, temp5,  temp6 !temp

                if (ios /= 0)  EXIT
                i=i+1

             end do length

             !reading datafile input
             REWIND(100)

             i_max= i-1
             !write(*,*) 'i_max', i_max
             allocate( proton_omega(i_max), proton_momentum(i_max), proton_gamma(i_max), proton_spectral(i_max) )
             !& proton_rho(i_max), proton_za(i_max) )

             allocat=.false.
          end if allo

          do i=lbound(proton_omega,dim=1),ubound(proton_omega,dim=1)
             read(100,1002,IOSTAT=ios) proton_omega(i), proton_momentum(i), &
                  & proton_gamma(i), proton_spectral(i) !, temp5, proton_za(i) !omega1, p1, Gamma, Spectral, rho, z/a
             !proton_rho(i)=temp5 !input rho is in fm^{-3}
             !      write(*,*) proton_rho(i)

             if (IOS /= 0) exit

          end do
          close(100)   !closing proton datafile

          fill:  if (fill_field) then
             tmp=proton_momentum(lbound(proton_momentum,dim=1))
             temp1=proton_momentum(lbound(proton_momentum,dim=1))
             j=lbound(momentum_field,dim=1)
             momentum_field(lbound(momentum_field,dim=1))=tmp

             do i=(lbound(proton_momentum,dim=1)+1),ubound(proton_momentum,dim=1)

                if ( (tmp) .NE. (proton_momentum(i)) ) then
                   if ( (temp1) .EQ. (proton_momentum(i))  ) EXIT
                   j=j+1
                   momentum_field(j)=proton_momentum(i)
                end if
                tmp=proton_momentum(i)
             end do

             tmp=proton_omega(lbound(proton_omega,dim=1))
             temp1=proton_omega(lbound(proton_omega,dim=1))
             j=lbound(omega_field,dim=1)
             omega_field(lbound(omega_field,dim=1))=tmp

             do i=(lbound(proton_omega,dim=1)+1),ubound(proton_omega,dim=1)

                if ( (tmp) .NE. (proton_omega(i)) ) then
                   if ( (temp1) .EQ. (proton_omega(i))  ) EXIT
                   j=j+1
                   omega_field(j)=proton_omega(i)
                end if
                tmp=proton_omega(i)
             end do

          end if fill

          do j=lbound(gamma_matrix,dim=2),ubound(gamma_matrix,dim=2)

             do k=lbound(gamma_matrix,dim=3),ubound(gamma_matrix,dim=3)

                do i=lbound(proton_omega,dim=1),ubound(proton_omega,dim=1)
                   if ( ( proton_omega(i) .EQ. omega_field(k)  )  .AND.  (proton_momentum(i) .EQ. momentum_field(j)) ) EXIT
                end do
                if (i .gt. ubound(proton_gamma,dim=1)) then
                   i=ubound(proton_gamma,dim=1)
                end if

                gamma_matrix(n_count,j,k)=proton_gamma(i)
                ! write(*,*) 'gamma_matrix(',n_count,',', j, ',', k ,')' ,momentum_field(j), omega_field(k)  , gamma_matrix(n_count,j,k)
             end do
          end do


          n_count=n_count +1
       end do densityloop
       minimalRho=0.01
       maximalRho=0.32

       minMomentum=momentum_field(lbound(momentum_field,dim=1) )
       maxMomentum=momentum_field( ubound(momentum_field,dim=1) )
       minEnergy=omega_field( lbound(omega_field,dim=1) )
       maximalEnergy=omega_field( ubound(momentum_field,dim=1) )

       deallocate(proton_omega, proton_momentum, proton_gamma, proton_spectral)

       drho=rho_field(2)-rho_field(1)
       domega=omega_field(2)-omega_field(1)
       dp=momentum_field(2)-momentum_field(1)

       initFlag=.false.

       write(*,*) 'min Momentum ', minMomentum, ' maxMomentum ' , maxMomentum, 'dp' , dp
       write(*,*) 'min Energy ', minEnergy, ' maxEnergy ' , maximalEnergy, 'domega' , domega
       write(*,*) 'min rho ', minimalRho, ' maxRho ' , maximalRho, 'drho' , drho
       write(*,*) 'denstity dependent width initialization completed!'

    end if init

    denstmp=pdens+ndens
    ! write(*,*) 'denstmp=' , denstmp
    !p_fermi_extern= hbarc*(1.5*pi**2*denstmp)**(1./3.)

    momentumtmp=sqrt(Dot_product(momLRF(1:3),momLRF(1:3) ) )





    if (denstmp .lt.minimalRho) then
       !  write(*,*) 'input density=',denstmp ,   ' is smallerr than maximalRho, using maximalRho=', minimalrho
       denstmp = minimalRho
    end if
    if (denstmp .gt.maximalRho) then
       !  write(*,*) 'input density=',denstmp ,   ' is higher than maximalRho, using maximalRho=', maximalrho
       denstmp = maximalRho
    end if

    n_rho=nint((denstmp-minimalRho)/drho)+1
    !write(*,*) 'n_rho=', n_rho
    !write(*,*)  'rho_field(n_rho)',rho_field(n_rho)
    if (present(mu_extern)) then
       omegatmp=momLRF(0)-mu_extern+omega_f_width_field(0.5*rho_field(n_rho) )
    else
       omegatmp=momLRF(0)-mN !+omega_f_width_field(0.5*rho_field(n_rho) )
    end if
    !p_fermi_intern=hbarc*(1.5*pi**2**rho_field(n_rho))**(1./3.)
    momentumtmp=momentumtmp !-p_fermi_extern+p_fermi_intern

    if (momentumtmp .lt. minMomentum ) then
       !  write(*,*) 'input Momentum', momentumtmp  ,'  is lower than minimal Momentum, using minMomentum' , minMomentum
       momentumtmp=minMomentum
    end if
    if (momentumtmp .gt. maxMomentum) then
       !  write(*,*) 'input Momentum', maxMomentum  ,'  is higher than maximal Momentum, using maxMomentum' , maxMomentum
       momentumtmp=maxMomentum
    end if

    if (omegatmp .lt. minEnergy ) then
       !  write(*,*) 'input energy', omegatmp  ,'  is lower than minimal Energy, using minEnergy' , minEnergy
       omegatmp=minEnergy
    end if
    if (omegatmp .gt. maximalEnergy ) then
       !  write(*,*) 'input Energy', omegatmp  ,'  is higher than maximal Energy, using maximalEnergy' , maximalEnergy
       if (form_Factor) then
          n_momentum=nint( (momentumtmp-minMomentum)/dp)+1
          n_energy=nint((maximalEnergy-minEnergy)/domega)+1
          gammatmp =gamma_matrix(n_rho,n_momentum,n_energy)
          proton_width_mdd=gammatmp*exp( -( (omegatmp-maximalEnergy)/0.025 )**2)+0.002
          return
       else
          omegatmp=maximalEnergy
       end if
    end if


    n_momentum=nint( (momentumtmp-minMomentum)/dp)+1
    !write(*,*) 'n_momentum=' ,n_momentum
    n_energy=nint((omegatmp-minEnergy)/domega)+1
    ! write(*,*) 'n_energy=', n_energy
    proton_width_mdd=gamma_matrix(n_rho,n_momentum,n_energy)
    !call splint2d(momentum_field,omega_field,gamma_matrix(n_rho,:,:),gam_deriv_mat(n_rho,:,:),momentumtmp,omegatmp,proton_width_medium_density_dependent)
    !write(*,*)  gamma_matrix(16,20,60)




  end function proton_width_mdd

 !*****************************************************************************
 ! generate strings from numbers (e.g. 2 -> "02"):
 !*****************************************************************************
 character (len=2) function nummern(iteration)

    integer :: iteration,n1,n0

    n1=int(iteration/10)
    n0=mod(iteration,10)

    nummern=achar(48+n1)//achar(48+n0)

  end function nummern

  !****************************************************************************
  !****************************************************************************
  real function omega_f_width_field(dens)
    use mean_field_baryon
    use cl_splines

    real, intent(in) :: dens
    real, save :: za=0.5
    real, Dimension(64) ,save:: dens_field, omega_f_field
    integer :: i
    logical, save :: initflag=.true.
    type(tspline),save :: spline
    logical :: successFlag
    integer :: error

    if (initflag) then
       do i=1,64
          dens_field(i)=float(i)*0.005
          omega_f_field(i)=omega_f_width(dens_field(i),za)
       end do
       spline=cl_initSpline(dens_field,omega_f_field)
       initflag=.false.
    end if

    omega_f_width_field=cl_spline(spline,dens,successFlag,error)
    if (.not.successFlag) call cl_error(error,'omega_f_width_field',dens)

  end function omega_f_width_field


end module pn_medium_width
