!******************************************************************************
!****m* /yukawa
! NAME
! module yukawa
! PURPOSE
! This module contains the implementation of the Yukawa potential.
! Only minor changes were made with respect to the old BUU code.
!******************************************************************************
module yukawa

  implicit none
  private


  !****************************************************************************
  !****g* yukawa/yukawaFlag
  ! PURPOSE
  ! Switches Yukawa potential on/off
  ! SOURCE
  logical, save :: yukawaFlag = .false.
  ! NOTES
  ! can be set in jobcard, namelist yukawa
  !****************************************************************************


  !****************************************************************************
  !****g* yukawa/smu
  ! PURPOSE
  ! Yukawa mass in fm**(-1). (range of potential)
  ! SOURCE
  real, save, public :: smu = 2.175
  ! NOTES
  ! can be set in jobcard, namelist yukawa
  !****************************************************************************

  !****************************************************************************
  !****g* yukawa/debug
  ! PURPOSE
  ! Switches debugging output on/off
  ! SOURCE
  logical, parameter :: debug = .false.
  !****************************************************************************


  integer,save :: maxx = 0
  integer,save :: maxz = 0
  real, save   :: ddgrid

  real, dimension(:,:,:), allocatable, save :: yup
  real, dimension(:,:,:), allocatable, save :: rhob_4

  real, save :: vzero                ! strength of potential
  real, dimension(1:2), save :: rmy  ! min and max eigenvalues

  logical, save :: initFlag = .true.


  public :: getYukawaFlag, getYukawaAlpha
  public :: updateYukawa
  public :: yukPot
  public :: tuneYukawaMass
  public :: cleanup


contains


  !****************************************************************************
  !****s* yukawa/tuneYukawaMass
  ! NAME
  ! subroutine tuneYukawaMass(amplitude)
  ! PURPOSE
  ! Tunes the yukawa mass "smu" such that the radii oszillations get
  ! minimalized. Output to "YukawaTuning.dat".
  ! INPUTS
  ! real, intent(in) :: amplitude ! Minimal radius of the nucleus in one run
  !****************************************************************************
  subroutine tuneYukawaMass(amplitude)
    real, intent(in) :: amplitude ! Amplitude of the nucleus oszillation

    real,save :: amplitude_prev
    real,save :: smu_prev
    logical :: firstTime=.true.
    real :: smu_new,steigung

    if (initFlag) call yuint

    if (firstTime) then
       smu_new=smu*1.10
       firsttime=.false.
       if (debug) open(17,file='YukawaTuning.dat')
    else
       if (debug) open(17,file='YukawaTuning.dat',position='append')
       ! Newton-Verfahren
       steigung=(amplitude-amplitude_prev)/(smu-smu_prev)
       if (steigung.lt.0.000001) then
          write(*,*) 'Gradient too small in tuneYukawaMass', steigung, smu, smu_prev, amplitude, amplitude_prev
          stop
       end if
       smu_new=smu-amplitude/steigung
    end if

    if (debug) then
       ! OUTPUT
       write(17,'(3(A,E12.5))') 'smu last run:', smu,'Amplitude last run:', amplitude, '=>smu set to:', smu_new
       close(17)
    end if

    smu_prev=smu
    amplitude_prev=amplitude
    smu=smu_new

  end subroutine tuneYukawaMass


  !****************************************************************************
  !****f* yukawa/getYukawaFlag
  ! NAME
  ! logical function getYukawaFlag()
  ! PURPOSE
  ! returns yukawaFlag
  !****************************************************************************
  logical function getYukawaFlag()

    if (initFlag) call yuint

    getYukawaFlag=yukawaFlag

  end function getYukawaFlag


  !****************************************************************************
  !****f* yukawa/getYukawaAlpha
  ! NAME
  !   real function getYukawaAlpha()
  ! PURPOSE
  ! returns the normalization for the mean field potential due to
  ! a Yukawa term in the  Hamiltonian
  ! NOTES
  ! alpha_Skyrme -> alpha_Skyrme-alpha_Yukawa
  !****************************************************************************
  real function getYukawaAlpha()
    use constants, only: pi, rhoNull

    if (initFlag) call yuint

    if (yukawaFlag) then
       getYukawaAlpha=abs(4.0*pi/smu**3*vzero*rhoNull)
    else
       getYukawaAlpha=0.  !if yukawa is switched off,
                          !then no renormalization necessary
    end if
  end function getYukawaAlpha


  !****************************************************************************
  !****f* yukawa/init
  ! NAME
  ! subroutine init
  ! PURPOSE
  ! initializes input
  !****************************************************************************
  subroutine init
    use output, only: Write_ReadingInput

    !**************************************************************************
    !****n* yukawa/Yukawa
    ! NAME
    ! NAMELIST /yukawa/
    ! PURPOSE
    ! Includes the input switches:
    ! * yukawaFlag
    ! * smu
    !**************************************************************************
    NAMELIST /yukawa/ yukawaFlag,smu

    logical, save :: readinFlag = .true.
    integer :: ios

    if (readinFlag) then
       call Write_ReadingInput('yukawa',0)
       rewind(5)
       read(5,nml=yukawa,iostat=ios)
       call Write_ReadingInput('yukawa',0,ios)
       write(*,*) 'yukawaFlag  = ',yukawaFlag,'.'
       write(*,*) 'yukawa Mass = ',smu,' fm^-1 .'
       call Write_ReadingInput('yukawa',1)
       readinFlag=.false.
    end if

  end subroutine init


  subroutine cleanup
    if (allocated(rhob_4)) deallocate(rhob_4)
    if (allocated(yup)) deallocate(yup)
  end subroutine


  !****************************************************************************
  !****f* yukawa/yuint
  ! NAME
  ! subroutine yuint
  ! PURPOSE
  ! Calculate the initial value of yukawa variables.
  !****************************************************************************
  subroutine yuint
    use constants, only: pi
    use output, only: Write_InitStatus
    use densitymodule, only: get_realGridSpacing, getGridPoints, densityField

    real, parameter    :: alp = 0.33  ! minimum eigenvalue
    real, parameter    :: bet = 5.7   ! maxmum eigenvalue
    integer, parameter :: in  = 10    ! number of iterations to determine Yukawa potential at first initialization

    integer, dimension(1:3) :: gridPoints
    real,    dimension(1:3) :: gridSpacing
    real :: gam, bety, alpy

    initFlag = .false.  ! set initFlag to false after first initilization

    call init

    if (.not. yukawaFlag) return

    call Write_InitStatus("Yukawa",0)

    call cleanup()

    gridPoints  = getGridPoints()
    gridSpacing = get_realGridSpacing()

    ddgrid = gridSpacing(1)
    maxx   = gridPoints(1)
    maxz   = gridPoints(3)

    allocate (rhob_4(-gridpoints(1):gridPoints(1),-gridPoints(2):gridPoints(2),-gridPoints(3):gridPoints(3)))
    allocate (yup(-maxx:maxx,-maxx:maxx,-maxz:maxz))

    if ((abs(gridSpacing(1)-gridSpacing(2)).gt.0.0001).or.(Abs(gridSpacing(1)-gridSpacing(3)).gt.0.0001) &
         &.or.(abs(gridSpacing(2)-gridSpacing(3)).gt.0.0001).or.(gridPoints(1).ne.gridPoints(2))) then
       write(*,*) 'Error in Yukawa'
       write(*,*) 'Need grid which is symmetric in XY-direction, with same grid spacing in all directions!'
       write(*,*) 'grid spacing=', gridSpacing(:)
       write(*,*) 'number of grid points=', gridpoints(:)
       stop
    end if

    alpy=alp+1.0/3.0*smu**2*ddgrid**2
    bety=bet+1.0/3.0*smu**2*ddgrid**2
    gam=(alpy*bety)**0.25*sqrt((alpy+bety)/2.0)
    rmy(1)=gam-sqrt(gam**2-alpy*bety)
    rmy(2)=gam+sqrt(gam**2-alpy*bety)


    if (debug) then
       open(20,File='YukawaOut.dat')

       write(20,'(''c:'',43x,''==== parameters of yukawa '','' ======= ================'')')
       write(20,'(''c:'',43x,''1) input parameters'')')
       write(20,'(''c:'',45x,''   in    ='',i4,'' initial iteratsumn'','' number'')') in
       write(20,'(''c:'',43x,''2) min and max eigen values'')')
       write(20,'(''c:'',45x,''   alpy  ='',f7.3,''   bety  ='',f7.3, &
            &   /''c:'',45x,''   rmy(1)='',e13.5,''  rmy(2)='',e13.5)') alpy,bety,rmy(1),rmy(2)
       write(20,'(''c:'',43x,''============================='',''========================='',15(/))')

       close(20)
    end if

    vzero=-0.378*smu**3/2.175**3

    ! initial value of yup(:,:,:) from zero-range solution:

    rhob_4 = densityField%baryon(0)

    yup(-maxx:maxx,-maxx:maxx,-maxz:maxz) = 4.0*pi/smu**3*vzero * rhob_4(-maxx:maxx,-maxx:maxx,-maxz:maxz)

    call updateYukawa(.false.,2*in)

    call Write_InitStatus("Yukawa",1)

  end subroutine yuint


  !****************************************************************************
  !****f* yukawa/updateYukawa
  ! NAME
  ! subroutine updateYukawa(start, numberIterationsIn)
  ! PURPOSE
  ! Update Yukawa potential, calculating the yukawa potential from (m) to
  ! (m+1) of the iteration.
  ! Evaluates the Yukawa potential according to the testparticle density
  ! (therefore we need the calculation of the dynamic density in updateDensity).
  !
  ! INPUTS
  ! * logical :: start   ---  .true. if the fields should be reinitialized
  ! * integer, optional :: numberIterationsIn   ---   number of iterations
  !
  ! NOTES
  ! we declare this subroutine recursive, since it call yuinit, which again
  ! calls this subroutine
  !****************************************************************************
  recursive subroutine updateYukawa(start, numberIterationsIn)
    use constants, only: pi
    use output, only: DoPR
    use densitymodule, only: densityField !,densitySwitch

    logical, intent(in) :: start
    integer, intent(in), optional :: numberIterationsIn

    integer :: mxx,mxz,maxm,j,ii,jj,i,ix,iy,iz,ni
    real, allocatable :: w(:), g(:), ee(:)
    real :: rml(2), rmr(2)    ! acceleration parameters, 1:even, 2:odd
    real :: yuppre, ci, ai
    integer,parameter :: numberIterations=2

    if (initFlag.or.start) call yuint

    if (.not.yukawaFlag) return

    mxx=maxx-1
    mxz=maxz-1
    maxm=(2*mxz+1)*(2*mxx+1)**2

    allocate(w(0:maxm))
    allocate(g(0:maxm))
    allocate(ee(1:maxm))

    if (DoPr(2)) write(*,*) 'Updating Yukawa'

!    If(densitySwitch.ne.1) then
!       write(*,*) 'Error in Yukawa'
!       Write(*,*) 'Need calculation of dynamic density. densitySwitch must be set to 1!'
!       Stop
!    end if

    if (present(numberIterationsIn)) then
       ni=numberIterationsIn
    else
       ni=numberIterations
    end if


    rhob_4=densityField%baryon(0)

    rml(1:2)=rmy(1:2)+1.0/3.0*smu**2*ddgrid**2
    rmr(1:2)=rmy(1:2)-2.0/3.0*smu**2*ddgrid**2

    do j=1,ni
       ii=1
       jj=j/2*2
       if (jj.eq.j) ii=2

       w(0)=0.0
       g(0)=0.0
       i=0
       do  iy=-mxx,mxx
          do  iz=-mxz,mxz
             do  ix=-mxx,mxx
                i=i+1
                ee(i) = vzero/smu*4.0*pi*rhob_4(ix,iy,iz)*ddgrid**2  &
                      + (rmr(ii)-4.0)*yup(ix,iy,iz)                  &
                      + yup(ix, iy+1, iz  )                          &
                      + yup(ix, iy-1, iz  )                          &
                      + yup(ix, iy  , iz+1)                          &
                      + yup(ix, iy  , iz-1)
                if (ix==-mxx) then
                   ee(i)=ee(i)+yup(ix-1,iy ,iz)
                   ai=0.0
                else
                   ai=-1.0
                end if
                if (ix==mxx) then
                   ee(i)=ee(i)+yup(ix+1,iy ,iz )
                   ci=0.0
                else
                   ci=-1.0
                end if
                g(i)=(ee(i)-ai*g(i-1))/(2.0+rml(ii)-ai*w(i-1))
                w(i)=ci/(2.0+rml(ii)-ai*w(i-1))
             end do
          end do
       end do


       i=i+1
       yuppre=0.0
       do  iy=mxx,-mxx,-1
          do  iz=mxz,-mxz,-1
             do  ix=mxx,-mxx,-1
                i=i-1
                yup(ix,iy,iz)=g(i)-w(i)*yuppre
                yuppre=yup(ix,iy,iz)
             end do
          end do
       end do

       i=i-1
       do  iz=-mxz,mxz
          do  ix=-mxx,mxx
             do  iy=-mxx,mxx
                i=i+1
                ee(i) = vzero/smu*4.0*pi*rhob_4(ix,iy,iz)*ddgrid**2  &
                      + (rmr(ii)-4.0)*yup(ix,iy,iz)                  &
                      + yup(ix+1, iy, iz  )                          &
                      + yup(ix-1, iy, iz  )                          &
                      + yup(ix  , iy, iz+1)                          &
                      + yup(ix  , iy, iz-1)
                if (iy==-mxx) then
                   ee(i)=ee(i)+yup(ix ,iy-1,iz )
                   ai=0.0
                else
                   ai=-1.0
                end if
                if (iy==mxx) then
                   ee(i)=ee(i)+yup(ix ,iy-1,iz )
                   ci=0.0
                else
                   ci=-1.0
                end if
                g(i)=(ee(i)-ai*g(i-1))/(2.0+rml(ii)-ai*w(i-1))
                w(i)=ci/(2.0+rml(ii)-ai*w(i-1))
             end do
          end do
       end do

       i=i+1
       yuppre=0.0
       do  iz=mxz,-mxz,-1
          do  ix=mxx,-mxx,-1
             do  iy=mxx,-mxx,-1
                i=i-1
                yup(ix,iy,iz)=g(i)-w(i)*yuppre
                yuppre=yup(ix,iy,iz)
             end do
          end do
       end do


       i=i-1
       do  ix=-mxx,mxx
          do  iy=-mxx,mxx
             do  iz=-mxz,mxz
                i=i+1
                ee(i) = vzero/smu*4.0*pi*rhob_4(ix,iy,iz)*ddgrid**2  &
                      + (rmr(ii)-4.0)*yup(ix,iy,iz)                  &
                      + yup(ix+1, iy  , iz)                          &
                      + yup(ix-1, iy  , iz)                          &
                      + yup(ix  , iy+1, iz)                          &
                      + yup(ix  , iy-1, iz)
                if (iz==-mxz) then
                   ee(i)=ee(i)+yup(ix ,iy ,iz-1)
                   ai=0.0
                else
                   ai=-1.0
                end if
                if (iz==mxz) then
                   ee(i)=ee(i)+yup(ix ,iy ,iz+1)
                   ci=0.0
                else
                   ci=-1.0
                end if
                g(i)=(ee(i)-ai*g(i-1))/(2.0+rml(ii)-ai*w(i-1))
                w(i)=ci/(2.0+rml(ii)-ai*w(i-1))
             end do
          end do
       end do

       i=i+1
       yuppre=0.0
       do  ix=mxx,-mxx,-1
          do  iy=mxx,-mxx,-1
             do  iz=mxz,-mxz,-1
                i=i-1
                yup(ix,iy,iz)=g(i)-w(i)*yuppre
                yuppre=yup(ix,iy,iz)
             end do
          end do
       end do
    end do
    deallocate(w,g,ee)

  end subroutine updateYukawa


  !****************************************************************************
  !****f* yukawa/yukpot
  ! NAME
  ! real function yukpot(rhap)
  ! PURPOSE
  ! Interface between yukawa routine and potential routine.
  ! Initializes Yukawa at first call.
  ! Returns the yukawa potential dependent on position.
  ! INPUTS
  ! * real, dimension(1:3),intent(in) :: rhap   ---   position
  !****************************************************************************
  real function yukpot(rhap)

    real, dimension(1:3),intent(in) :: rhap

    integer :: ix,iy,iz

    if (initFlag) call yuint

    yukpot = 0.

    if (.not.yukawaFlag) return

    ix=nint(rhap(1)/ddgrid)
    iy=nint(rhap(2)/ddgrid)
    iz=nint(rhap(3)/ddgrid)

    if (abs(ix)<maxx .and. abs(iy)<maxx .and. abs(iz)<maxz) &
      yukpot = yup(ix,iy,iz)

  end function yukpot


end module yukawa
