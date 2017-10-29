!******************************************************************************
!****m* /formfactors_A_piPlus_neutron
! NAME
! module formfactors_A_pPlus_neutron
!
! PURPOSE
! This module formfactors_A reads in the form factors for gamma* N -> N pi out of text files
! and includes modules to return the form factors as functions of theta, s and Q^2.
!
! For the isospin channel gamma p -> n pi^+.
!
! INPUTS
! None
!******************************************************************************
module formfactors_A_piPlus_neutron
  implicit none

  private

  public :: getA_piPlus_neutron, cleanUp

  logical, save :: debug=.false.

  integer , parameter :: maxGridPoints=100000 ! maximal number of gridpoints in the input dat-file
  character(60) :: modName='formfactors_A_pi+_neutron'

  ! Fields to save the grid obtained by MAID (initialized at first call):
  complex, save , dimension(:,:), Allocatable :: A                   ! The amplitude field
  real,    save , dimension(:)  , Allocatable :: W, theta, Q2, s, t  ! Fields to store W,theta,q^2,s and t

contains

  !****************************************************************************
  !****f* formfactors_A_piPlus_neutron/getA_piPlus_neutron
  ! NAME
  ! function getA_piPlus_neutron(thetaIn,sIn,QSquaredIn)
  !
  ! PURPOSE
  ! This function returns the invariant Amplitudes A_1, ... A_6 of the MAID analysis.
  !
  ! INPUTS
  ! * real  ::  thetaIn    ! Theta scattering angle of the pion relative to q in the CM system of the hadronic vertex
  ! * real  ::  QSquaredIn ! Q^2=-q^mu q_mu  of the gamma
  ! * real  ::  sIn        ! Mandelstam s of the hadronic vertex: gamma* N-> pi N
  ! OUTPUT
  ! complex, dimension(1:6) :: getA   ! The complex amplitudes A_1 to A_6 in units of:
  ! * [GeV**-2] for A(1)
  ! * [GeV**-4] for A(2) and A(5)
  ! * [GeV**-3] for A(3), A(4) and A(5)
  !****************************************************************************
  function getA_piPlus_neutron (thetaIn,sIn,QSquaredIn) result(getA)
    use output, only: Write_InitStatus
    use formfactors_A_input, only: ValidateVars

    complex, dimension(1:6) :: getA
    real, intent(in) :: thetaIn, QSquaredIn,sIn

    logical, save :: initFlag=.true.

    real,save :: W_max, theta_max, q2_max                              ! maximal values for W,theta,q^2
    real,    save :: dW,dTheta,dQ2                                     ! steps size for W,theta,q^2
    integer, save :: numSteps_q,numSteps_t,numSteps_w                  ! number of different values for q^2,theta, w

    integer :: i,k,j,jj
    real :: WIn,Qin2
    integer,dimension(1:8) :: n_index
    logical,parameter:: linearInterpolation=.true.
    real :: factor, factorsum
    real, dimension(1:8) :: theta_cube, q2_cube, w_cube


    ! (1) Read amplitudes from file
    if (initFlag) then
       call Write_InitStatus('MAID Amplitudes for gamma p -> n pi+',0)
       call init_datFile_pipn
       initFlag=.false.
       write(*,'(A,3F9.5)') 'dW,dTheta,dQ2 of the MAID-File:',dW,dTheta,dQ2
       write(*,'(A,3I9)')   'Number of steps in W,Theta,Q2 :',numsteps_w,numSteps_t,numSteps_q
       write(*,'(A,F9.5)')  'W_max:',W_max
       call Write_InitStatus('MAID Amplitudes for gamma p -> n pi+',1)
    end if

    Win=sqrt(sin)
    Qin2=qSquaredIn

    ! (2) Check that input values are not out of bounds:
    if (.not. ValidateVars (win, thetaIn, Qin2, w(1)-dw/2, w_Max, theta(1), theta_max, q2(1), q2_max, modName)) then
      getA = 0.
      return
    end if

    ! (3) Return form factors:
    if (linearInterpolation.and.(.not.(win.lt.w(1)))) then
       ! Interpolate in a linear manner between the closest grid points to get the amplitude at (thetaIn,WIn,QsquaredIn)
       n_Index=get_neighbors(thetaIn,WIn,Qin2)
       theta_cube=theta(n_index)
       w_cube=w(n_index)
       q2_cube=q2(n_index)
       do
          getA=0.
          factorsum=0
          do i=1,8
             factor=1.
             factor=factor*( 1.  - abs(theta_cube(i)-thetaIn   )/dTheta     )
             factor=factor*( 1.  - abs(w_cube(i)    -wIn       )/dW         )
             factor=factor*( 1.  - abs(q2_cube(i)   -Qin2      )/dQ2        )
             factorsum=factorsum+factor
             if (debug) write(*,*) i, factor,factorsum, Real(A(1,n_index(i)))
             getA=getA+A(:,n_index(i))*factor
          end do
          if (abs(factorsum-1.).gt.0.001) then
             write(*,*) 'Crucial error in linearInterpolation of getA: stop',factorsum,'in',modname
             if (debug) stop
             debug=.true.
             do jj=1,8
                write(*,*) theta_cube(jj),w_cube(jj),q2_cube(jj)
             end do
             write(*,*) 'Input  :', thetaIn,Win,'(s=',sin,')',qSquaredIn,Qin2
             write(*,*) 'Maximal:',theta_max, w_max,q2_max
             write(*,*) 'Minimal:',theta(1),W(1),q2(1)
             write(*,*) 'Stop!'
          else
             exit
          end if
       end do
    else
       ! Just get closest grid point and take the value at this point
       i=get_i(thetaIn,WIn,Qin2,.false.)
       getA=A(:,i)
    end if


    ! (4) Debugging output
    if (debug) then
       write(*,'(A,F12.5)') 's=',sIn
       write(*,'(A,F12.5)') 'W=',sqrt(sIN)
       write(*,'(A,2F12.5)') 'Q^2=',qSquaredIn,Qin2
       write(*,'(A,F12.5)') 'theta=',theta(i)
       do k=1,6
          write(*,'(A,I6,A,2E19.12)') 'A_{',k,'}(theta)=',A(k,i)
          select case (k)
          case (1)
             write(*,'(A,2E19.12,A)') '=',A(k,i)*0.197**2, 'fm^2'
          case (2,5)
             write(*,'(A,2E19.12,A)') '=',A(k,i)*0.197**4, 'fm^4'
          case default
             write(*,'(A,2E19.12,A)') '=',A(k,i)*0.197**3, 'fm^3'
          end select
       end do
    end if

  contains

    !**************************************************************************
    !****if* getA_piPlus_neutron/get_i
    ! NAME
    ! integer function get_i(tin,Win,Qin,int_nint_switch)
    !
    ! PURPOSE
    ! Evaluates index which corresponds to the grid point which lies closest to the input values:
    ! * For int_nint_switch=.false. it's the overall closest point
    ! * For int_nint_switch=.true.  it's the closest point which is also in all dimensions closer to 0 than the input point.
    !**************************************************************************
    integer function get_i(tin,Win,Qin,int_nint_switch)
      real, intent(in) :: Tin,Win,Qin
      logical, intent(in) :: int_nint_switch
      integer :: t_index, W_index, Q_index

      if (int_nint_switch) then
         t_index=int((Tin-theta(1))/dTheta)+1
         w_index=int((Win-w(1))/dW)
         q_index=int((Qin-q2(1))/dQ2)
      else
         t_index=nint((Tin-theta(1))/dTheta)+1
         w_index=nint((Win-w(1))/dW)
         q_index=nint((Qin-q2(1))/dQ2)
      end if
      t_index=min(max(1,t_index),numsteps_t)
      w_index=min(max(0,w_index),numsteps_w)
      q_index=min(max(0,q_index),numsteps_q)

      ! Loop order is : q,w ,t . Therefore:
      get_i=q_index*numSteps_t*numSteps_w+w_index*numSteps_t+t_index

      if (debug) then
         write(*,*) 'theta=',tIn, theta(get_i),dTheta
         write(*,*) 'q**2 =',qIn, q2(get_i)   ,dQ2
         write(*,*) 'W    =',wIn, w(get_i)    ,dW
      end if

    end function get_i

    !**************************************************************************
    !****if* getA_piPlus_neutron/get_neighbors
    ! NAME
    ! function get_neighbors(tin,Win,Qin)
    !
    ! PURPOSE
    ! Returns the eight grid points which are defining the smallest cube in which the input
    ! point (tin,win,Qin) can be locked into. These are the neighboring grid-points of our
    ! input point.
    ! OUTPUT
    ! integer,dimension(1:8) :: get_neighbors
    !**************************************************************************
    function get_neighbors(tin,Win,Qin)
      real, intent(in) :: Tin,Win,Qin
      integer,dimension(1:8) :: get_neighbors
      real :: tin_new,win_new,qin_new
      integer :: closest
      real , parameter :: eps=0.000000001
      get_neighbors=0

      ! (1) Check first whether the input point is an exact grid-point in some dimension.
      ! If it is exact, then move it a little bit by some negligible epsilon
      ! to make the algorithm work properly.
      ! Take care that we do this by respecting the boudaries. Points at the upper boundary must
      ! be shifted downwards, points at the lower boundary upwards.

      ! Get closest grid point:
      closest=get_i(tin,Win,Qin,.false.)

      ! Compare each dimension:
      if (abs(theta(closest)-tin).lt.eps) then
         tin_new=tin-eps
         if (abs(tin-theta(1)).lt.eps) then
            tin_new=theta(1)+eps
         end if
         if (abs(tin-theta_max).lt.eps) then
            tin_new=theta_max-eps
         end if
      else
         tin_new=tin
      end if
      if (abs(w(closest)-win).lt.eps) then
         win_new=win-eps
         if (abs(win-w(1)).lt.eps) then
            win_new=w(1)+eps
         end if
         if (abs(win-w_max).lt.eps) then
            win_new=w_max-eps
         end if
      else
         win_new=win
      end if
      if (abs(q2(closest)-qin).lt.eps) then
         qin_new=qin-eps
         if (abs(qin-q2(1)).lt.eps) then
            qin_new=q2(1)+eps
         end if
         if (abs(qin-q2_max).lt.0) then
            qin_new=q2_max-eps
         end if
      else
         qin_new=qin
      end if


      ! (2) get neighbors
      do i=0,1
         do j=0,1
            do k=0,1
               get_neighbors(4*i+2*j+k+1)=get_i(tin_new+dtheta*i,Win_new+dW*j,Qin_new+dQ2*k,.true.)
            end do
         end do
      end do
    end function get_neighbors


    !**************************************************************************
    !****if* getA_piPlus_neutron/init_datFile
    ! NAME
    ! subroutine init_datFile
    !
    ! PURPOSE
    ! Reads in the MAID input file, writes into FF_A.dat the control output. It expects an input file with the following format:
    ! * first some comment lines (doesn't matter how many)
    ! * Data lines with the following order : W, ThetaCM , Q^2 ,s ,t, Re(A_1), Im(A_1),...Re(A_6), Im(A_6)
    !
    ! The loop structure to generate this input file should be :
    ! * thetaCM-Loop, W-Loop, Q^2-Loop
    !
    ! OUPTUT
    ! Initializes the fields:
    ! * complex, save , dimension(:,:), Allocatable :: A
    ! * real,    save , dimension(:)  , Allocatable :: W, theta, Q2, s, t
    ! * real,save :: W_max, theta_max, q2_max             ! maximaal values for W,theta,q^2
    ! * real,    save :: dW,dTheta,dQ2                    ! steps size for W,theta,q^2
    ! * integer, save :: numSteps_q,numSteps_t,numSteps_w ! number of different values for q^2,theta, w
    !**************************************************************************
    subroutine init_datFile_pipn
      use formfactors_A_input, only: get_Filename
      use bzip

      integer :: ios,i,k
      real, dimension(1:12,1:maxGridPoints) :: spalte
      real, dimension(1:maxGridPoints) :: WRead, ThetaRead , Q2Read ,sRead ,tRead
      logical,parameter :: convertToGeV=.true.

      type(bzFile) :: f
      character(len=300) :: line
      integer :: ll

      f = bzOpenR(trim(get_Filename('pipn')))

      numsteps_t=1
      numsteps_q=1
      numsteps_w=1
      i=1
      do
         if (i<2) ll = 0
         call bzReadLine(f,line,ll)
         if (ll==0) cycle

         read(line(1:ll),*,IOSTAT=IOS)  WREAD(i), ThetaREAD(i) , Q2READ(i) ,sREAD(i) ,tREAD(i), spalte(1:12,i)

         if (ios.eq.0) then
            ! Start determine grid spacings
            if (i.eq.1) then
               w_max=WREAD(i)
               theta_max= ThetaREAD(i)
               q2_max=  Q2READ(i)
            else if (i.gt.1) then
               if (WREAD(i)-W_max.gt.0.000001) then
                  dW=WREAD(i)-W_max
                  W_max=WREAD(i)
                  numsteps_w=numsteps_w+1
               end if
               if (ThetaREAD(i)-Theta_Max.gt.0.000001) then
                  dTheta=ThetaREAD(i)-Theta_max
                  theta_max=ThetaREAD(i)
                  numsteps_t=numsteps_t+1
               end if
               if (Q2READ(i)-Q2_max.gt.0.000001) then
                  dq2=Q2READ(i)-Q2_max
                  Q2_max=Q2READ(i)
                  numsteps_q=numsteps_q+1
               end if
            end if
            ! End determine grid spacings
            i=i+1
         else
            write(*,'(i3,i3,A)') i, ios, line(1:ll)
         end if
         if (f%EOF) exit
      end do
      i=i-1

      call bzCloseR(f)

      allocate(A(1:6,1:i))

      allocate(W(1:i))
      w=wread(1:i)/1000.
      W_max=w_max/1000.
      dw=dw/1000.

      allocate(theta(1:i))
      theta=thetaRead(1:i)

      allocate(Q2(1:i))
      q2=q2read(1:i)

      allocate(s(1:i))
      s=sread(1:i)

      allocate(t(1:i))
      t=tread(1:i)


      do k=1,i
         A(1,k)=CMPLX(Spalte(1,k),Spalte(2,k))
         A(2,k)=CMPLX(Spalte(3,k),Spalte(4,k))
         A(3,k)=CMPLX(Spalte(5,k),Spalte(6,k))
         A(4,k)=CMPLX(Spalte(7,k),Spalte(8,k))
         A(5,k)=CMPLX(Spalte(9,k),Spalte(10,k))
         A(6,k)=CMPLX(Spalte(11,k),Spalte(12,k))
      end do

      if (debug) then
         ! Control output
         open(101,file='FF_A_pi+neutron.dat')
         do i= lbound(A,dim=2),ubound(A,dim=2)
            write(100+1,'(17E12.4)')     W(i), Theta(i) , Q2(i) ,s(i) ,t(i), (A(k,i),k=1,6)
         end do
         close(101)
      end if
      if (convertToGeV) then
         ! Convert everything from fm to GeV^-1
         ! 197 MeV fm =1 => fm=1/(0.197 GeV)
         A(1,:)=A(1,:)/(0.197**2)
         A(2,:)=A(2,:)/(0.197**4)
         A(3,:)=A(3,:)/(0.197**3)
         A(4,:)=A(4,:)/(0.197**3)
         A(5,:)=A(5,:)/(0.197**4)
         A(6,:)=A(6,:)/(0.197**3)
      end if
    end subroutine init_datFile_pipn
  end function getA_piPlus_neutron


  subroutine cleanUp
    if (allocated(A)) deallocate(A)
    if (allocated(W)) deallocate(W)
    if (allocated(theta)) deallocate(theta)
    if (allocated(Q2)) deallocate(Q2)
    if (allocated(s)) deallocate(s)
    if (allocated(t)) deallocate(t)
  end subroutine


end module formfactors_A_piPlus_neutron
