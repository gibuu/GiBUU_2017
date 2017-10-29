!******************************************************************************
!****m* /deltaWidth
! NAME
! module deltaWidth
! PURPOSE
! Implements the routines for the medium width of the Delta-Resonance P33(1232).
!******************************************************************************
module deltaWidth

  implicit none
  private

  !****************************************************************************
  !****n* deltaWidth/deltawidth
  ! NAME
  ! NAMELIST deltaWidth
  ! PURPOSE
  ! Includes the switch:
  ! * deltaSwitch
  !****************************************************************************

  !****************************************************************************
  !****g* deltaWidth/deltaSwitch
  ! PURPOSE
  ! Switch for different prescriptions for the delta width.
  ! SOURCE
  !
  integer, save :: deltaSwitch=3
  ! NOTES
  ! * 1 = use Oset self energies+BUU input
  ! * 2 = Spreading potential
  ! * 3 = use Oset self energies
  ! * 4 = density dependent with BUU input
  !****************************************************************************


  logical, parameter :: debug=.false.

  integer, parameter :: nm=100  ! Number of mass bins
  integer, parameter :: np=30  ! Number of momentum bins
  integer, parameter :: nrho=40   ! Number of rho bins

  real, parameter :: massmin=1.071      ! Minimal Mass
  real, parameter :: rhomin=0.0      ! Minimal rho
  real, parameter :: pmin=.0          ! Minimal momentum

  real, parameter :: drho=0.1
  real, parameter :: dp=0.05
  real, parameter :: dm=0.01

!!$  ! To use Effenberger's input
!!$  integer, parameter :: nm=51  ! Number of mass bins
!!$  integer, parameter :: np=19  ! Number of momentum bins
!!$  integer, parameter :: nrho=19   ! Number of rho bins
!!$
!!$  real, parameter :: massmin=1.081      ! Minimal Mass
!!$  real, parameter :: rhomin=0.0168      ! Minimal rho
!!$  real, parameter :: pmin=.025          ! Minimal momentum
!!$
!!$  real, parameter :: drho=0.2
!!$  real, parameter :: dp=0.05
!!$  real, parameter :: dm=0.01


  !****************************************************************************
  !****g* deltaWidth/delSpread
  ! PURPOSE
  ! Spreading potential parameter.
  ! SOURCE
  !
  real, parameter :: delspread=0.04
  !****************************************************************************


  real, save ::   medwidth(0:nm,0:np,0:nrho)  ! Delta width as function of mass, momentum and density
  real, save ::   g1pi(0:nm,0:np,0:nrho)      ! partial width (Delta-> N Pi) width as function of mass, momentum and density

  real, save, dimension(:,:,:), ALLOCATABLE  ::   gcoll1, gcoll2   ! collisional width
  ! also dimension(0:nm,0:np,0:nrho), if allocated

  logical,save  :: initFlag=.true.

  public :: delta_nucleonPion
  public :: delta_fullWidth
  public :: deloset
  public :: osetDelta_used
  public :: deltaOset_ND_NN
  public :: deltaOset_ND_ND
  public :: GetRhoValue
  public :: GetRhoBin
  public :: GetMaxQ_Delta

contains

  !****************************************************************************
  !****g* deltaWidth/osetDelta_used
  ! NAME
  ! logical function osetDelta_used
  ! PURPOSE
  ! return true, if deltaSwitch==3
  !****************************************************************************
  logical function osetDelta_used()
    if (deltaSwitch.eq.3) then
       osetDelta_used=.true.
    else
       osetDelta_used=.false.
    end if
  end function osetDelta_used

  !****************************************************************************
  !****f* deltaWidth/delta_nucleonPion
  ! NAME
  ! real function delta_nucleonPion(mass,momentum,density)
  ! PURPOSE
  ! Returns partial width Delta->N Pi as a function of mass, momentum and
  ! density.
  ! INPUTS
  ! * real :: mass -- mass of Delta
  ! * real :: density -- density at position
  ! * real, dimension(1:3) :: momentum -- momentum in LRF (only abs value used)
  !****************************************************************************
  real function delta_nucleonPion(mass,momentum,density)
    use constants, only: rhoNull

    real, intent(in) :: mass
    real, intent(in) :: density
    real, intent(in),dimension(1:3)  :: momentum

    integer :: iMass, iMom, iDens

    if (initFlag) then
       call initDelWidth
       initFlag=.false.
    end if

    iMass =Max(Min(NINT((mass-massMin)/dm),nm),0)
    iMom  =Max(Min(NINT((sqrt(Dot_Product(momentum, momentum))-pMin)/dp),np),0)
    iDens =Max(Min(NINT((density-rhoMin)/(drho*rhoNull)),nrho),0)

    delta_nucleonPion=g1Pi(iMass,iMom,iDens)
    if (debug)  write(*,*) iMass, iMom, iDens,delta_nucleonPion
    if (debug)  write(*,*) mass, Dot_Product(momentum, momentum), density

  end function delta_nucleonPion

  !****************************************************************************
  !****f* deltaWidth/delta_fullWidth
  ! NAME
  ! real function delta_fullWidth(mass,momentum,density)
  ! PURPOSE
  ! Returns Delta width as a function of mass, momentum and density.
  ! INPUTS
  ! * real :: mass -- mass of Delta
  ! * real :: density -- density at position
  ! * real, dimension(1:3) :: momentum -- momentum in LRF (only abs value used)
  !****************************************************************************
  real function delta_fullWidth(mass,momentum,density)
    use constants, only: rhoNull

    real, intent(in) :: mass
    real, intent(in) :: density
    real, intent(in),dimension(1:3)  :: momentum

    integer :: iMass, iMom, iDens

    if (initFlag) then
       call initDelWidth
       initFlag=.false.
    end if

    iMass =Max(Min(NINT((mass-massMin)/dm),nm),0)
    iMom  =Max(Min(NINT((sqrt(Dot_Product(momentum, momentum))-pMin)/dp),np),0)
    iDens =Max(Min(NINT((density-rhoMin)/(drho*rhoNull)),nrho),0)

    delta_fullWidth=medWidth(iMass,iMom,iDens)
    if (debug)  write(*,*) iMass, iMom, iDens,delta_fullWidth
  end function delta_fullWidth

  !****************************************************************************
  !****is* deltaWidth/InitDelWidth
  ! NAME
  ! subroutine InitDelwidth
  ! PURPOSE
  ! Evaluates Delta width and stores it as a function of mass, momentum and
  ! density into the field medwidth.
  !****************************************************************************
  subroutine InitDelwidth
    use inputGeneral, only: path_to_input
    use constants, only: rhoNull
    use output

    integer :: i,j,k
    real      :: rho
    real      ::  mass,imsig2,imsig3,imsigq

    integer :: ios ! checks file behavior
    real , parameter :: betaD=0.
    real , parameter :: alphaD=0.

    real :: xColl1, xColl2

    logical,parameter :: dobetter=.true.
    ! If this switch is .true., then we do not use Effenbergers input since
    ! it does not include the case rho=0.
    ! and has just a small grid in rho.
    ! If false, then we try to use Effenbergers file.
    ! Set these variables before setting dobetter=.false.:
    ! massMin=1.0810
    ! rhoMin =0.0168
    ! pMin   =0.0250
    ! dp     =0.05
    ! drho   =0.2
    ! dm     =0.01
    ! nm     =51
    ! np     =19
    ! nrho   =19
    ! So you get the same grid with which Effenbergers Input file was
    ! generated!


    NAMELIST /deltaWidth/ deltaSwitch

    call Write_ReadingInput('deltaWidth',0)
    rewind(5)
    read(5,nml=deltaWidth,IOSTAT=IOS)
    call Write_ReadingInput('deltaWidth',0,IOS)

    write(*,*) 'Set deltaSwitch      :' ,deltaSwitch
    call Write_ReadingInput('deltaWidth',1)

    ios=0
    call Write_ReadingInput('delwidth_full.dat',0)
    open(13,file=trim(path_to_input)//'/delwidth_full.dat', status='old',iostat=ios)
    if (ios.ne.0.or.dobetter) then
       ! ios.ne.0: Effenberger's file not available or some other error in reading it.
       ! Calculate the in-medium pion-nucleon partial width:
       call get_g1Pi()
       ! Calculate the collision broadening due to N Delta collisions:
       if (deltaswitch.eq.1.or.deltaswitch.eq.4) then
          !call evalualuateGColl()
          write(*,*) 'Problem in delWidth: Not yet implemented.'
          write(*,*) 'This can only be implemented by reading out Effenbergers input file. Therefore:'
          write(*,*) 'Set dobetter=.false. and take care of the array bounds. See comments in the code.'
          stop
       end if
       close(13)
    else
       close(13)
       ! ios.eq.0: Effenberger's file is available
       call readEffenberger()
    end if
    call Write_ReadingInput('delwidth_full.dat',1)


    do i=0,nm
       mass=float(i)*dm+massmin
       do j=0,np
          do k=0,nrho
             rho=float(k)*drho+rhoMin !*     rho in units of rho0

             select case (deltaSwitch)
             case (1) ! use Oset self energies

                call deloset(mass,rho*rhoNull,imsig2,imsig3,  imsigq)
                xcoll1=2.*imsig3+gcoll1(i,j,k)*(1.+5./3.*(mass-1.076)*rho)
                xcoll2=2.*imsigq

             case (2) ! Spreading potential

                xcoll1=2.*delspread*rho
                xcoll2=0.

             case (3) ! use Oset self energies

                call deloset(mass,rho*rhoNull,imsig2,imsig3, imsigq)
                xcoll1=2.*(imsig2+imsig3)
                xcoll2=2.*imsigq

             case (4)

                g1pi(i,j,k)=g1pi(i,j,k)*(1.+betad*rho)
                xcoll1=gcoll1(i,j,k)*(1.+alphad*rho)
                xcoll2=0.

             end select

             medwidth(i,j,k) = g1pi(i,j,k) + xcoll1 + xcoll2

          end do
       end do
    end do

    if (debug) then
       write(33,*)
       do i=0,nm
          do j=0,np
             do k=0,nrho
                write(33,'(4E14.5)') massMin+dm*i,pMin+dp*j,rhoMin+drho*rhoNull*k,g1pi(i,j,k)
             end do
          end do
       end do
    end if

!    call PrintFullWidth()

    return
  end subroutine InitDelwidth


  !****************************************************************************
  !****************************************************************************
  subroutine readEffenberger()
    use inputGeneral, only: path_to_input
    use constants, only: rhoNull

    integer :: ios
    integer :: nm2,np2,nrho2
    real    :: dm_effe,dp_effe,drho_effe
    integer :: ii,jj,kk
    integer :: i,j,k
    integer :: i_effe,j_effe,k_effe

    real , Allocatable, dimension (:,:,:) :: gcoll1_effe, gcoll2_effe ,g1pi_effe  ! Fields for collisional width
    real , Allocatable, dimension (:) ::  mass_effe,impuls_effe,dens_effe  ! Fields for collisional width
    real :: mass, impuls,dens

    open(13,file=trim(path_to_input)//'/delwidth_full.dat', status='old',iostat=ios)
    if (ios.ne.0) then
       write(*,*) 'File delwidth_full.dat not available'
       stop
    end if

    read(13,*) nm2,np2,nrho2

    allocate(gcoll1_effe(1:nm2,1:np2,1:nrho2))
    allocate(gcoll2_effe(1:nm2,1:np2,1:nrho2))
    allocate(g1pi_effe(1:nm2,1:np2,1:nrho2))
    allocate(dens_effe(1:nrho2))
    allocate(mass_effe(1:nm2))
    allocate(impuls_effe(1:np2))

    allocate(gcoll1(0:nm,0:np,0:nrho))
    allocate(gcoll2(0:nm,0:np,0:nrho))

    do i=1,nm2
       do j=1,np2
          do k=1,nrho2
             read(13,*) mass_effe(i),impuls_effe(j),dens_effe(k),g1pi_effe(i,j,k),gcoll1_effe(i,j,k),gcoll2_effe(i,j,k)
          end do
       end do
    end do

    ! Make ouput
    open(40,file='delwidth_effe.dat',iostat=ios)
    write(40,*) '# Masses   ', mass_effe
    write(40,*) '# Momenta  ', impuls_effe
    write(40,*) '# Densities', dens_effe
    dm_effe=mass_effe(2)-mass_effe(1)
    drho_effe=dens_effe(2)-dens_effe(1)
    dp_effe=impuls_effe(2)-impuls_effe(1)
    do i=1,nm2
       mass=mass_effe(1)+(i-1)*dm_effe
       do j=1,np2
          impuls=impuls_effe(1)+(j-1)*dp_effe
          do k=1,nrho2
             dens=dens_effe(1)+(k-1)*drho_effe
             write(40,'(6F15.8)') mass,impuls,dens, g1pi_effe(i,j,k),gcoll1_effe(i,j,k),gcoll2_effe(i,j,k)
          end do
       end do
    end do
    close(40)

    ! Assign to real fields
    do i=0,nm
       mass=i*dm+massmin
       do j=0,np
          impuls=j*dp+pMin
          do k=0,nrho
             dens=k*drho*rhoNull+rhoMin
             ! Search closest point in the grid which we read in
             i_effe=ubound(mass_effe,dim=1)
             do ii=lbound(mass_effe,dim=1),ubound(mass_effe,dim=1)
                if (mass_effe(ii).gt.mass) then
                   if (ii.eq.lbound(mass_effe,dim=1)) then
                      i_effe=ii
                   else if (abs(mass_effe(ii)-mass).gt.abs(mass_effe(ii-1)-mass)) then
                      i_effe=ii-1
                   else
                      i_effe=ii
                   end if
                   exit
                end if
             end do
             j_effe=ubound(impuls_effe,dim=1)
             do jj=lbound(impuls_effe,dim=1),ubound(impuls_effe,dim=1)
                if (impuls_effe(jj).gt.impuls) then
                   if (jj.eq.lbound(impuls_effe,dim=1)) then
                      j_effe=jj
                   else if (abs(impuls_effe(jj)-impuls).gt.abs(impuls_effe(jj-1)-impuls)) then
                      j_effe=jj-1
                   else
                      j_effe=jj
                   end if
                   exit
                end if
             end do

             k_effe=ubound(dens_effe,dim=1)
             do kk=lbound(dens_effe,dim=1),ubound(dens_effe,dim=1)
                if (dens_effe(kk).gt.dens) then
                   if (kk.eq.lbound(dens_effe,dim=1)) then
                      k_effe=kk
                   else if (abs(dens_effe(kk)-dens).gt.abs(dens_effe(kk-1)-dens)) then
                      k_effe=kk-1
                   else
                      k_effe=kk
                   end if
                   exit
                end if
             end do
!             write(*,*) mass, mass_effe(i_effe),impuls,impuls_effe(j_effe),dens,dens_effe(k_effe)
             ! Evalutate the new grid at the points of the closest grid points of the readin grid.
             g1pi  (i,j,k)=g1pi_effe  (i_effe,j_effe,k_effe)
             gcoll1(i,j,k)=gcoll1_effe(i_effe,j_effe,k_effe)
             gcoll2(i,j,k)=gcoll2_effe(i_effe,j_effe,k_effe)
          end do
       end do
    end do


  end subroutine readEffenberger

  !****************************************************************************
  !****s* deltaWidth/get_g1P
  ! NAME
  ! subroutine get_g1P()
  ! PURPOSE
  ! Evaluates the partial width of a delta decaying into a pion nucleon
  ! systemfor a whole array and writes result to file as a function of mass,
  ! momentum and density. If it finds the information already stored in a
  ! file, it just reads in this file
  ! NOTES
  ! This subroutine calculates the in-medium width of a resonance in a fermi
  ! gas model.
  !****************************************************************************
  subroutine get_g1Pi()

    use random
    use idTable, only: delta
    use constants, only: pi, rhoNull, mN, mPi
    use lorentzTrafo, only: lorentz
    use baryonWidth, only: fullWidthBaryon
    use inputGeneral, only: path_to_input
    use output

    integer, parameter :: nmonte1=10000     ! Monte Carlo Steps for free width

    real, dimension( 1:3) :: beta
    real, dimension( 0:3) :: pcmn

    real :: mass1, mass2

    real :: mass, rho , gtot, enuk

    real :: eres, pres, pfermi!, efermi
    real :: pcms, cost, sum
    integer :: mom_index, mass_index, rho_index

    integer :: i,ios
    real :: massMin_file,rhoMin_file,pMin_file,dm_file,drho_file,dp_file,nm_file,nrho_file,np_file

    mass1=mN
    mass2=mPi

    ! Check whether we stored the information in an earlier run
    call Write_ReadingInput('delwidth_1pi.dat',0)
    open(13,file=trim(path_to_input)//'/delwidth_1pi.dat', status='old',iostat=ios)
    if (ios.eq.0) then ! File is available
       ! -> Check that it uses the same dp, drho, dm and so forth...
       read(13,*)    massMin_file,rhoMin_file,pMin_file,dm_file,drho_file,dp_file,nm_file,nrho_file,np_file
       if ( (abs(massMin_file-massMin).lt.0.0001).and.&
            &(abs(rhoMin_file-rhoMin).lt.0.0001).and. &
            &(abs(pMin_file-pMin).lt.0.0001).and.&
            &(abs(dm_file-dm).lt.0.0001).and.&
            &(abs(drho_file-drho).lt.0.0001).and.&
            &(abs(dp_file-dp).lt.0.0001).and.&
            &(abs(nm_file-nm).lt.0.0001).and.&
            &(abs(np_file-np).lt.0.0001).and.&
            &(abs(nrho_file-nrho).lt.0.0001) &
            & ) then
          ! Input file is suited
          write(*,*) 'deltaWidth/get_g1P: Reading input file ',trim(path_to_input)//'/delwidth_1pi.dat'
          massloop_input: do mass_index=0,nm
             momentumloop_input: do mom_index=0,np
                densityLoop_input : do rho_Index=0,nrho
                   read(13,iostat=ios,fmt='(E15.6)') g1pi(mass_index,mom_index,rho_index)
                   if (ios.ne.0) then
                      write(*,*) 'Error in reading delwidth_1pi.dat, seems corrupt! Delete file! Then start again!',ios
                      write(*,*) mass_index,mom_index,rho_index
                      stop
                   end if
                end do densityLoop_input
             end do momentumloop_input
          end do massloop_input
          close(13)
          call Write_ReadingInput('delwidth_1pi.dat',1)
          return ! Go back and don't calculate since readin worked out
      else
        write(*,*) 'Old input file:', massMin_file,rhoMin_file,pMin_file,dm_file,drho_file,&
             & dp_file,nm_file,nrho_file,np_file
       end if
    end if
    close(13)
    call Write_ReadingInput('delwidth_1pi.dat',1)

    ! Calculate the field and write the results to delwidth_1pi.dat, such that it's available next time
    write(*,*) 'Calcute and write input file for the next time...'
    open(13,file=trim(path_to_input)//'/delwidth_1pi.dat.NEW', iostat=ios)
    write(13,'(6E15.5,3I8)')   massMin,rhoMin,pMin,dm,drho,dp,nm,nrho,np
    massloop: do mass_index=0,nm

       mass=float(mass_index)*dm+massmin

       ! Total width in vacuum
       gtot= FullWidthBaryon(delta,mass)

       if (mod(Mass_index,10).eq.0) write(*,*) mass_index,'/',nm

       momentumloop: do mom_index=0,np
          pres=float(mom_index)*dp+pMin
          eres=sqrt(mass**2+pres**2)
          densityLoop : do rho_Index=0,nrho
             rho=float(rho_index)*drho*rhoNull+rhoMin

             pfermi=(3./2.*pi**2*rho)**(1./3.)*0.197
             !efermi=sqrt(mass1**2+pfermi**2)

             !#################################################################
             !  spontaneous decay
             !  integrate over decay angle
             !#################################################################

             sum=0.
             pcms=sqrt(max(0.,(mass**2-(mass1+mass2)**2)*   (mass**2-(mass1-mass2)**2)/4./mass**2))

             beta(1)=0.
             beta(2)=0.
             beta(3)=pres/eres

             do i=1,nmonte1
                cost=(rn()-0.5)*2.
                pcmn(1)=pcms*sqrt(max(1.-cost**2,0.))
                pcmn(2)=0.
                pcmn(3)=pcms*cost
                enuk=sqrt(mass1**2+pcms**2)
                pcmn(0)=enuk

                call lorentz(-beta,pcmn)
                ! Check Pauli Blocking :
                if (Dot_Product(pcmn(1:3),pcmn(1:3)).gt.pfermi**2) sum=sum+1.
             end do

             g1pi(mass_index,mom_index,rho_index)=sum/float(nmonte1)*gtot
             write(13,'(E15.6)') g1pi(mass_index,mom_index,rho_index)
          end do densityLoop
       end do momentumloop
    end do massloop
    close(13)
    write(*,*) 'stop!'
    stop
    return
  end subroutine get_g1Pi

!!$ subroutine evaluateGColl()
!!$             do l=1,nmonte2
!!$                debugflag2=.false.
!!$                id1=2
!!$                id2=1
!!$                idi=1
!!$                idj=2
!!$                is1=3
!!$                is2=1
!!$                iz1=1
!!$                iz2=nint(rn(0))
!!$                izt=iz1+iz2
!!$                em1=mass
!!$                em2=mass1
!!$
!!$                !*momentum of incoming nucleon
!!$                pnuk=(rn(0))**(1./3.)*pfermi
!!$                cost=(rn(0)-0.5)*2.
!!$                pnuki(1)=pnuk*sqrt(max(1.-cost**2,0.))
!!$                pnuki(2)=0.
!!$                pnuki(3)=pnuk*cost
!!$                enuk=sqrt(mass1**2+pnuk**2)
!!$                beta(1)=pnuki(1)/(enuk+eres)
!!$                beta(2)=0.
!!$                beta(3)=(pnuki(3)+pres)/(enuk+eres)
!!$                do k2=1,3
!!$                   pcm(k2)=pnuki(k2)
!!$                end do
!!$                ecm=enuk
!!$                call lorentz(beta(1),beta(2),beta(3),pcm(1),pcm(2),   pcm(3),ecm)
!!$                prcm=sqrt(pcm(1)**2+pcm(2)**2+pcm(3)**2)
!!$                srts=sqrt((enuk+eres)**2-pnuki(1)**2-  (pnuki(3)+pres)**2)
!!$                pinitial2=(srts**2-(mass+mass1)**2)*(srts**2-(mass-  mass1)**2)/4./srts**2
!!$                pinitial=sqrt(pinitial2)
!!$                if(abs(pinitial2-pcm(1)**2-pcm(2)**2-pcm(3)**2).gt. 1e-05) then
!!$                   write(*,*)'problems in delwidth',beta,pnuki,pcm, srts,mass,mass1,enuk,eres,pinitial2
!!$                   stop
!!$                end if
!!$                vrel=pinitial*srts/eres/enuk
!!$                indexs=min(nint((srts-sigs0)/delsigs),nsmax)
!!$
!!$                srtfree=srts
!!$                if(debugflag2) write(*,*)'vor N D',beta,pcm,pinitial,
!!$                &                 srts,mass1,mass
!!$
!!$                !*N Delta -> N N
!!$                !*masses of outgoing particles
!!$                id3=1
!!$                id4=1
!!$                em3=rmass
!!$                em4=rmass
!!$                iscatt=1
!!$                !*angle
!!$                cost=2.*(rn(0)-0.5)
!!$                call dimido(dsdo,srts,mass,cost)
!!$                sigma1=2.*dsdo*8./3.*iifac(1,1,iz3,iz4,0)
!!$                !*rotation
!!$                phi=rn(0)*2.*pi
!!$                sint=sqrt(max(1.-cost**2,0.))
!!$                pscattc(1)=sint*cos(phi)
!!$                pscattc(2)=sint*sin(phi)
!!$                pscattc(3)=cost
!!$                costk=pcm(3)/prcm
!!$                if(pcm(2).ne.0.or.pcm(1).ne.0) then
!!$                   phik=atan2(pcm(2),pcm(1))
!!$                else
!!$                   phik=0.
!!$                end if
!!$                call rotate(cost,phik,pscattc,pscatt)
!!$                !*check pauli-blocking
!!$                pcms=sqrt(srts**2/4.-mass1**2)
!!$                do k2=1,3
!!$                   pcmn2(1,k2)=pscatt(k2)*pcms
!!$                   pcmn2(2,k2)=-pcmn2(1,k2)
!!$                end do
!!$                ecm1=sqrt(mass1**2+pcms**2)
!!$                ecm2=sqrt(mass1**2+pcms**2)
!!$
!!$                call lorentz(-beta(1),-beta(2),-beta(3),pcmn2(1,1), pcmn2(1,2),pcmn2(1,3),ecm1)
!!$                call lorentz(-beta(1),-beta(2),-beta(3),pcmn2(2,1), pcmn2(2,2),pcmn2(2,3),ecm2)
!!$
!!$                !*checks
!!$                if(abs(ecm1+ecm2-enuk-eres).gt.1e-04) then
!!$                   write(*,*)'problems energy conservation in  delwidth',ecm1,ecm2,enuk,eres
!!$                   stop
!!$                end if
!!$
!!$                block1=1.
!!$                if(ecm1.lt.efermi.or.ecm2.lt.efermi) block1=0.
!!$
!!$                !*N Delta -> N Delta
!!$                !*masses of outgoing particles
!!$                id3=1
!!$                id4=2
!!$                iscatt=2
!!$                minmass=rmass+pmass
!!$                maxmass=srts-rmass
!!$                em4=rn(0)*(maxmass-minmass)+minmass
!!$                em3=rmass
!!$                call manley(em4,2,0,0,ratio,gamtot,0,404)
!!$                spectral=2./pi*em4**2*gamtot/((em4**2- barprop1(2,1)**2)**2+gamtot**2*em4**2)
!!$                pcms2=(srts**2-(em3+em4)**2)*(srts**2-(em3-em4)**2)   /4./srts**2
!!$                pcms=sqrt(max(pcms2,0.))
!!$                csrt=srts-2.*rmass
!!$                if(pinitial.gt.0) then
!!$                   sigma2=(35.0/(1.+csrt*100.0)+20.0)/pinitial* iifac(1,2,iz3,iz4,0)*pcms*spectral*    (maxmass-minmass)
!!$                else
!!$                   sigma2=0.
!!$                end if
!!$
!!$                !*angle
!!$                cost=(rn(0)-0.5)*2.
!!$                sint=sqrt(max(1.-cost**2,0.))
!!$                phi=2.*pi*rn(0)
!!$                pscatt(3)=cost
!!$                pscatt(2)=sint*sin(phi)
!!$                pscatt(1)=sint*cos(phi)
!!$                !*check pauli-blocking
!!$                do k2=1,3
!!$                   pcmn2(1,k2)=pscatt(k2)*pcms
!!$                   pcmn2(2,k2)=-pcmn2(1,k2)
!!$                end do
!!$                ecm1=sqrt(mass1**2+pcms**2)
!!$                ecm2=sqrt(em4**2+pcms**2)
!!$
!!$                if(debugflag2) write(*,*)'vor lorentz'
!!$                call lorentz(-beta(1),-beta(2),-beta(3),pcmn2(1,1),     pcmn2(1,2),pcmn2(1,3),ecm1)
!!$                call lorentz(-beta(1),-beta(2),-beta(3),pcmn2(2,1),     pcmn2(2,2),pcmn2(2,3),ecm2)
!!$
!!$                !*checks
!!$                if(abs(ecm1+ecm2-enuk-eres).gt.1e-04) then
!!$                   write(*,*)'problems energy conservation in     delwidth 2',ecm1,ecm2,enuk,eres
!!$                   stop
!!$                end if
!!$
!!$                block2=1.
!!$                if(ecm1.lt.efermi) block2=0.
!!$
!!$                sigtot=sigma1+sigma2
!!$
!!$                if(sigtot.gt.sigd) then
!!$                   sigma1=sigma1*sigd/sigtot
!!$                   sigma2=sigma2*sigd/sigtot
!!$                end if
!!$
!!$                gcoll1(i,j,k)=gcoll1(i,j,k)+    eres/mass*vrel*rho*(0.197)**3*    sigma1*block1/0.389       /float(nmonte2)
!!$                gcoll2(i,j,k)=gcoll2(i,j,k)+    eres/mass*vrel*rho*(0.197)**3*    sigma2*block2/0.389       /float(nmonte2)
!!$
!!$             end do
!!$
!!$             medwidth(i,j,k)=g1pi(i,j,k)+gcoll1(i,j,k)+gcoll2(i,j,k)
!!$ end subroutine evaluateGColl()

  !****************************************************************************
  !****s* deltaWidth/deloset
  ! NAME
  ! subroutine deloset(mass,rho,imsig2,imsig3,imsigq)
  ! PURPOSE
  ! Evaluates Delta Width according to Oset and Salcedo
  ! NPA 468 (1987) 631-652
  !****************************************************************************
  subroutine deloset(mass,rho,imsig2,imsig3,imsigq)
    use constants, only: rhoNull, pi, mN, mPi

    real, intent(in) ::  mass,rho
    real, intent(out) :: imsig2,imsig3,imsigq

    real :: ompi,x,kf,T
    integer :: i,j
    real a(3,5),b(5)

    ! See table 2 in the paper:
    data ((a(i,j),j=1,5),i=1,3)    &
    &     / -5.19,   1.06, -13.46,  0.382, -0.038,   &
         &  15.35,  -6.64,  46.17, -1.322,  0.205,     &
         &   2.06,  22.66, -20.34,  1.466,  0.613/


    kf=(3./2.*pi**2*rho)**(1./3.)*0.197

    ! Energy of pion (see Effenberger Phd thesis page 248):
    ompi=(mass**2-mN**2-mpi**2)/(3./5.*kf**2/2./mN+mN)/2.
    ! T/m of pion:
    T=min(0.45,ompi-mPi)  ! Osets calculation only valid up to about here!  
    ! At 450 MeV the collisional Delta width becomes small and is frozen to this
    ! value 
    T=max(0.0,T)   !ensures that at threshold m = mpi + mN T is not negative
                   !because of average kinetic energy of nucleon
    
    x=T/mPi

    do i=1,5
       b(i)=0.
       do j=1,3
          b(i)=b(i)+a(j,i)*x**(-j+3)
       end do
    end do

    if (b(1)*rho.gt.0) then
       imsigq=max(b(1)*(rho/rhoNull)**b(4),0.)/1000.
    else
       imsigq=0.
    end if
    if (b(2)*rho.gt.0) then
       imsig2=max(b(2)*(rho/rhoNull)**b(5),0.)/1000.
    else
       imsig2=0.
    end if
    if (b(3)*rho.gt.0) then
       imsig3=max(b(3)*(rho/rhoNull)**(2.*b(5)),0.)/1000.
    else
       imsig3=0.
    end if

  end subroutine deloset

  !****************************************************************************
  !****is* deltaWidth/deltaOset_ND_ND
  ! NAME
  ! real function deltaOset_ND_ND(part)
  ! PURPOSE
  ! Evaluates ND->ND cross sections based
  ! on Delta Width according to Oset and Salcedo  NPA 468 (1987) 631-652
  ! INPUTS
  ! * type(particle),intent(in),dimension(1:2) :: part ! scattering particles
  ! OUTPUT
  ! Cross section in mb
  !****************************************************************************
  real function deltaOset_ND_ND(part)
    use particleDefinition

    type(particle),intent(in),dimension(1:2) :: part

    deltaOset_ND_ND=cross(part,1)

  end function deltaOset_ND_ND

  !****************************************************************************
  !****is* deltaWidth/deltaOset_ND_NN
  ! NAME
  ! real function deltaOset_ND_NN(part)
  ! PURPOSE
  ! Evaluates ND->NN cross sections based
  ! on Delta Width according to Oset and Salcedo  NPA 468 (1987) 631-652
  ! INPUTS
  ! * type(particle),intent(in),dimension(1:2) :: part ! scattering particles
  ! OUTPUT
  ! Cross section in mb
  !****************************************************************************
  real function deltaOset_ND_NN(part)
    use particleDefinition

    type(particle),intent(in),dimension(1:2) :: part

    deltaOset_ND_NN=cross(part,2)

  end function deltaOset_ND_NN

  !****************************************************************************
  !****is* deltaWidth/cross
  ! NAME
  ! function cross(part,med,switch) result(sigma)
  ! PURPOSE
  ! Evaluates ND->Nd (switch=1) or ND->NN (switch=2) cross sections based
  ! on Delta Width according to Oset and Salcedo  NPA 468 (1987) 631-652
  ! INPUTS
  ! * integer :: switch
  ! * type(particle),intent(in),dimension(1:2) :: part ! scattering particles
  ! OUTPUT
  ! Cross section in mb
  !****************************************************************************
  function cross(part,switch) result(sigma)
    use relativeVelocity
    use dichteDefinition
    use densityModule, only: densityAT
    use particleDefinition
    use constants, only:GeVSquared_times_mb
    use idTable, only: delta

    real :: sigma
    integer :: switch
    type(particle),intent(in),dimension(1:2) :: part ! scattering particles
    real :: mass, imsig2,imsig3,imsigq  ! ,pdel
    real, dimension(1:3) :: pdel_vec
    real :: velo_dens_p, velo_dens_n

    integer :: charge
    real, dimension(1:3) :: pos_delta,pos_nuc
    type(dichte) :: dens_atNuc,dens_atDelta
    real :: rho,rho_p, rho_n

    if (((sum(part%charge).gt.2).or.(sum(part%charge).lt.0)).and.(switch.eq.2)) then
       ! No charge conservation possible in ND->NN: Delta++ p or Delta-n initial state
       sigma=0.
       return
    end if

    sigma=0.
    if (part(1)%ID.eq.delta) then
       mass=part(1)%mass
!        pdel=absMom(part(1))
       pdel_vec=part(1)%momentum(1:3)
       charge=part(1)%charge
       pos_delta=part(1)%position
       pos_nuc  =part(2)%position
    else if (part(2)%ID.eq.delta) then
       mass=part(2)%mass
!        pdel=absMom(part(2))
       pdel_vec=part(2)%momentum(1:3)
       charge=part(2)%charge
       pos_delta=part(2)%position
       pos_nuc  =part(1)%position
    else
       write(*,*) "Wrong ID's in deltaOset_ND_ND. STOP", part%ID
       stop
    end if

    ! We use the density at the point of the delta to evaluate gamma.
    dens_atDelta=densityAt(pos_delta)
    rho=dens_atDelta%proton(0)+dens_atDelta%neutron(0)
    call deloset(mass,rho,imsig2,imsig3,imsigq)

    ! We use the density at the nuc position to evaluate v_rel(average) since we average over the
    ! Fermi sea at the nuc position.
    ! NOTE: For full ensemble both positions become the same, but not for parallel ensemble!
    dens_atNuc=densityAt(pos_nuc)
    rho=dens_atNuc%proton(0)+dens_atNuc%neutron(0)
    rho_p=dens_atNuc%proton(0)
    rho_n=dens_atNuc%neutron(0)

    if (switch.eq.1) then       ! ***********Quasielastic: ND-> ND ******
       velo_dens_p=vrel_times_rho(rho_p,mass,pdel_vec)
       velo_dens_n=vrel_times_rho(rho_n,mass,pdel_vec)
       if ((velo_dens_n+velo_dens_p).gt.1E-15) sigma=2.*imsigq/(velo_Dens_p+velo_Dens_n) !in GeV**(-2)

    else if (switch.eq.2) then  ! ***********Absorptive : ND-> NN********
       ! ATTENTION:
       ! Delta^++ and Delta^- can only be absorbed on neutrons OR protons
       ! Henceforth we need to evaluate sigma
       ! such that gamma=sigma*rho_{proton/neutron}*v.
       if (charge.eq.2) then
          ! Delta^++ can only be absorbed on a neutron
          velo_dens_n=vrel_times_rho(rho_n,mass,pdel_vec)
          if (velo_dens_n.gt.1E-15) sigma=2.*imsig2/velo_Dens_n  !in GeV**(-2)
       else if (charge.eq.-1) then
          ! Delta^- can only be absorbed on a proton
          velo_dens_p=vrel_times_rho(rho_p,mass,pdel_vec)
          if (velo_dens_p.gt.1E-15) sigma=2.*imsig2/velo_Dens_p   !in GeV**(-2)
       else
          velo_dens_p=vrel_times_rho(rho_p,mass,pdel_vec)
          velo_dens_n=vrel_times_rho(rho_n,mass,pdel_vec)
          if ((velo_dens_n+velo_dens_p).gt.1E-15)  sigma=2.*imsig2/(velo_Dens_p+velo_Dens_n) !in GeV**(-2)
       end if
    else
       write(*,*) 'Wrong switch in deltaWidth/cross, switch'
       stop
    end if

    sigma=sigma/GeVSquared_times_mb !convert to mb

    ! Debugging:
    !write(300+switch,*) sqrts(part),sigma
  end function cross


  !****************************************************************************
  !****is* deltaWidth/PrintFullWidth
  ! NAME
  ! subroutine PrintFullWidth
  ! PURPOSE
  ! Print the full width in order to plot it with gnuplot as 2D plot:
  ! * version 1: (M,P) for different rho, file "DeltaFullWidth.gp1.dat"
  ! * version 2: (M,rho) for different P, file "DeltaFullWidth.gp2.dat"
  ! In order to select the third variable, you have to use the parameter
  ! 'index' in gnuplot.
  !****************************************************************************
!   subroutine PrintFullWidth
!     use constants, only : rhoNull
!
!     integer :: iM,iP,iRho
!     real :: M,P,rho
!
!     open(112,file="DeltaFullWidth.ALL.dat", status="unknown")
!     write(112,*) rhoMin,drho,nrho
!     write(112,*) massMin,dM,nM
!     write(112,*) pMin,dp,nP
!     write(112,*) deltaSwitch
!     close(112)
!
!     open(112,file="DeltaFullWidth.gp1.dat", status="unknown")
!     do iRho=0,nrho
!        rho = rhoMin+drho*rhoNull*iRho
!        do iM=0,nM
!           M = massMin+dm*iM
!           do iP=0,nP
!              P = pMin+dp*iP
!              write(112,'(1P,4E14.5)') P,M,rho,medWidth(iM,iP,iRho)
!           end do
!           write(112,*)
!        end do
!        write(112,*)
!        write(112,*)
!     end do
!     close(112)
!
!     open(112,file="DeltaFullWidth.gp2.dat", status="unknown")
!     do iP=0,nP
!        P = pMin+dp*iP
!        do iRho=0,nrho
!           rho = rhoMin+drho*rhoNull*iRho
!           do iM=0,nM
!              M = massMin+dm*iM
!              write(112,'(1P,4E14.5)') P,M,rho,medWidth(iM,iP,iRho)
!           end do
!           write(112,*)
!        end do
!        write(112,*)
!        write(112,*)
!     end do
!     close(112)
!
!   end subroutine PrintFullWidth


  !****************************************************************************
  !****s* deltaWidth/GetMaxQ_Delta
  ! NAME
  ! subroutine GetMaxQ_Delta(mass0,gamma0,rho,BinM,BinMaxQ)
  ! PURPOSE
  ! Calculate the maximal values of the Q weight for bins according BinM
  ! INPUTS
  ! * integer :: ID -- ID of resonance
  ! * real :: mass0 -- pole mass
  ! * real :: gamma0 -- width at pole mass
  ! * real :: rho -- density
  ! * real, dimension(:) :: BinM -- array with boundaries for M binning
  ! OUTPUT
  ! * real, dimension(:) :: BinMaxQ -- the maximal Q values for each bin.
  !
  ! NOTES
  ! * The size of BinMaxQ has to be at least the size of BinM minus 1.
  ! * It first calculates Q at the boundaries, then it iterates over the
  !   tabulated width values in order to take into account, that the Q
  !   value may be larger inbetween the boundaries.
  ! * if the Q value is maximal at the upper bound, we store its value as
  !   -Q.
  !
  !****************************************************************************
  subroutine GetMaxQ_Delta(mass0,gamma0,rho,BinM,BinMaxQ)
    use MassAssInfoDefinition, only: MassAssInfoQ
    use constants, only: rhoNull

    real, intent(in) :: mass0,gamma0,rho
    real, dimension(:),intent(in)  :: BinM
    real, dimension(:),intent(out) :: BinMaxQ

    integer :: iB, iM,iM1,iM2, iMom, iRho
    real :: mass,gamma, Q

    iRho = Max(Min(NINT((rho-rhoMin)/(drho*rhoNull)),nrho),0)
    iMom = 0

    ! 1) calculate Q at the boundaries:

    do iB=1,size(BinM)
       iM = Max(Min(NINT((BinM(iB)  -massMin)/dm),nm),0)
       BinMaxQ(iB) = 0
       do iMom=0,nP
          gamma = medWidth(iM,iMom,iRho)
          Q = MassAssInfoQ(mass0,gamma0,BinM(iB),gamma)
          if (Q.gt.BinMaxQ(iB)) BinMaxQ(iB) = Q
       end do

    end do

    do iB=1,size(BinM)-1
       if (BinMaxQ(iB+1).gt.BinMaxQ(iB)) BinMaxQ(iB) = -BinMaxQ(iB+1) ! sign!!
    end do

    ! 2) calculate Q between the boundaries:

    do iB=1,size(BinM)-1
       iM1 = Max(Min(NINT((BinM(iB)  -massMin)/dm),nm),0)
       iM2 = Max(Min(NINT((BinM(iB+1)-massMin)/dm),nm),0)

       iMom = 0

       do iM=iM1+1,iM2
          mass  = iM*dm+massMin
          do iMom=0,nP
             gamma = medWidth(iM,iMom,iRho)
             Q = MassAssInfoQ(mass0,gamma0,mass,gamma)
             if (Q.gt.abs(BinMaxQ(iB))) BinMaxQ(iB) = Q
             Q = MassAssInfoQ(mass0,gamma0,mass+dm/2,gamma) ! respect bin width
             if (Q.gt.abs(BinMaxQ(iB))) BinMaxQ(iB) = Q
          end do
       end do
    end do

  end subroutine GetMaxQ_Delta

  !****************************************************************************
  !****s* deltaWidth/GetRhoValue
  ! NAME
  ! subroutine  GetRhoValue(rhomin,drho,nrho)
  ! PURPOSE
  ! return the variables concerning the density of the parametrization
  ! INPUTS
  ! none
  ! OUTPUT
  ! * real :: rhomin, drho
  ! * integer :: nrho
  !****************************************************************************
  subroutine GetRhoValue(rhomin_,drho_,nrho_)
    real,intent(out) :: rhomin_,drho_
    integer, intent(out) :: nrho_
    rhomin_ = rhomin
    drho_ = drho
    nrho_ = nrho
  end subroutine GetRhoValue

  !****************************************************************************
  !****s* deltaWidth/GetRhoBin
  ! NAME
  ! subroutine GetRhoBin(rho,iRho)
  ! PURPOSE
  ! return the bin nummer connected to a certain density
  ! INPUTS
  ! * real :: rho -- the density
  ! OUTPUT
  ! integer :: iRho -- the bin nummer
  !****************************************************************************
  subroutine GetRhoBin(rho,iRho)
    use constants, only: rhoNull
    real,intent(in) :: rho
    integer, intent(out) :: iRho
    iRho =Max(Min(NINT((rho-rhoMin)/(drho*rhoNull)),nrho),0)
  end subroutine GetRhoBin

end module deltaWidth
