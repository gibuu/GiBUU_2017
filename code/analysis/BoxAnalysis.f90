!******************************************************************************
!****m* /BoxAnalysis
! NAME
! module BoxAnalysis
! PURPOSE
!******************************************************************************
module BoxAnalysis

  use histf90
  use histMPf90
  use histMC
  use TmunuDefinition, only: tTmunuNmu, fillTmunu, headTmunu

  implicit none

  private

  public :: DoBoxAnalysisTime

  type(histogram), save :: hMassRho
  type(histogramMC), save :: hMCMomPion

  integer, parameter :: nSet = 5
  type(histogramMP), dimension(nSet), save :: hMP_ESet, hMP_pSet
  logical, dimension(nSet), save :: useSet

  type(tTmunuNmu), dimension(:), allocatable, save :: arrTmunuNmu
  type(tTmunuNmu), dimension(2), save :: arrTmunuNmu_hadr

  logical, save :: initFlag=.true.
  logical, parameter :: do_Tmunu_pirho=.false.


  !****************************************************************************
  !****g* BoxAnalysis/do_Tmunu
  ! SOURCE
  logical,save :: do_Tmunu=.false.
  ! PURPOSE
  ! Switch for Tmunu output. default: Only one file for all ensemble!
  ! you may change this with the flag perEnsemble_Tmunu
  !****************************************************************************

  !****************************************************************************
  !****g* BoxAnalysis/perEnsemble_Tmunu
  ! SOURCE
  logical,save :: perEnsemble_Tmunu=.false.
  ! PURPOSE
  ! Switch for Tmunu output. One file per ensemble!
  !
  ! NOTES
  ! this may slow down the execution dramatically, since huge output to the
  ! hard drive is induced.
  ! You may observe this, if e.g the cpu load drops permanently to 30%.
  ! Thus: switch it on, only if you want it!
  !****************************************************************************

  !****************************************************************************
  !****g* BoxAnalysis/do_P
  ! SOURCE
  logical,save :: do_P=.false.
  ! PURPOSE
  ! Switch for dN/p^2 dp output
  !****************************************************************************

  !****************************************************************************
  !****g* BoxAnalysis/do_velrel
  ! SOURCE
  logical,save :: do_velrel=.false.
  ! PURPOSE
  ! Switch for calculating velrel
  !****************************************************************************

  !****************************************************************************
  !****g* BoxAnalysis/do_Cumulants
  ! SOURCE
  logical,save :: do_Cumulants=.false.
  ! PURPOSE
  ! Switch for calculating cumulants
  !****************************************************************************

contains
  !****************************************************************************
  !****s* BoxAnalysis/DoBoxAnalysisTime
  ! NAME
  ! subroutine DoBoxAnalysisTime(realPart,timestep)
  ! PURPOSE
  !****************************************************************************
  subroutine DoBoxAnalysisTime(realPart,timestep)

    use CallStack, only: TRACEBACK
    use densityModule, only: gridsize
    use history, only: history_getParents
    use output, only: Write_InitStatus, intToChar4 !, WriteParticleVector
    use particleDefinition
    use collisionReporter, only: CR_write

    type(particle),dimension(:,:),intent(in), target :: realPart
    integer, intent(in) :: timestep

    integer, save :: nHist
    real, save :: boxVol ! the volume of the box in fm^3

    integer, save :: nEns,nPart
    real, save :: mulfak

    integer :: i,j, iID, iCh, iSet
    real :: mom0,mom, mass
    type(particle), POINTER :: pPart
    integer :: parents(1:3)

    type(tTmunuNmu) :: TmunuNmu

    real, dimension(10) :: countQ, countNch
    real, dimension(10,2) :: countNpi

    if (initFlag) then
       call Write_InitStatus('BoxAnalysis',0)
       initFlag=.false.

       call readInput

       if ((do_Tmunu) .and. (perEnsemble_Tmunu) .and. (nEns > 9999)) then
          call TRACEBACK("BoxAnalysis: Tmunu not prepared for more than 9999 ensembles.")
       end if

       nEns  = size(realPart,dim=1)
       nPart = size(realPart,dim=2)
       boxVol = 8.*gridsize(1)*gridsize(2)*gridsize(3)
       mulfak = 1.0/(nEns*boxVol)

       !----- mass distribution: -----
       call CreateHist(hMassRho, "mass(rho)", 0.0, 2.5, 0.01)
       call CreateHistMC(hMCMomPion, "momentum(pion)", 0.0, 2.5, 0.02, 6)
       hMCMomPion%yDesc(1:6) = (/ "original  ",  &
            "rho       ", "sigma     ", "other dec ", &
            "pi pi     ", "other coll" /)

       !----- multiplicities: -----

       do iSet=1,nSet
          if (.not.useSet(iSet)) cycle
          call CreateHistMP(hMP_ESet(iSet), "dN/pE dE",  0.0, 2.5, 0.02, iSet)
          call CreateHistMP(hMP_pSet(iSet), "dN/p^2 dp", 0.0, 2.5, 0.02, iSet)
          open(123,file="BoxAnalysis_Mult_Set"//achar(48+iSet)//".dat", status="unknown")
          call WriteHistMP_Names(iSet,123)
          close(123)
       end do

       !----- hydro tensors: -----
       if (do_Tmunu) then
          allocate( arrTmunuNmu(nEns) )

          open(123,file="BoxAnalysis_Tmunu.dat", status="unknown")
          write(123,'(A)') headTmunu
          close(123)
          if (do_Tmunu_pirho) then
             open(123,file="BoxAnalysis_Tmunu.pion.dat", status="unknown")
             write(123,'(A)') headTmunu
             close(123)
             open(123,file="BoxAnalysis_Tmunu.rho.dat", status="unknown")
             write(123,'(A)') headTmunu
             close(123)
          end if

          if (perEnsemble_Tmunu) then
             do i=1,nEns
                open(123,file="BoxAnalysis_Tmunu."//intToChar4(i)//".dat", status="unknown")
                write(123,'(A)') headTmunu
                close(123)
             end do
          end if
       end if

       call Write_InitStatus('BoxAnalysis',1)
    end if

    do iSet=1,nSet
       if (.not.useSet(iSet)) cycle
       call ClearHistMP(hMP_pSet(iSet))
       call ClearHistMP(hMP_ESet(iSet))
    end do

    call ClearHistMC(hMCMomPion)
    call ClearHist(hMassRho)

    if (do_Tmunu) then
       arrTmunuNmu = TmunuNmu ! set all values to zero
       arrTmunuNmu_hadr = TmunuNmu ! set all values to zero
    end if

    if (do_cumulants) then
       countQ = 0.0
       countNch = 0.0
       countNpi = 0.0
    end if

    !
    ! accumulate data:
    !

    do i=1,nEns
       do j=1,nPart
          pPart => realPart(i,j)
          if (pPart%Id <  0) exit
          if (pPart%Id <= 0) cycle

          mom = absMom(pPart)
          mom0 = pPart%momentum(0)

          select case (pPart%ID)
          case (101)
             parents = history_getParents(pPart%history)
             if (parents(2) == 0) then
                select case (parents(1))
                case (0)
                   iCh = 1
                case (103)
                   iCh = 2
                case (104)
                   iCh = 3
                case default
                   iCh = 4
                end select
             else
                if (parents(1)==101 .and. parents(2)==101) then
                   iCh = 5
                else
                   iCh = 6
                end if
             end if
             call AddHistMC(hMCMomPion, mom, iCh, 1.0/(mom**2))
          case (103)
             mass = sqrtS(pPart)
             call AddHist(hMassRho, mass, 1.0)
          end select

          do iSet=1,nSet
             if (.not.useSet(iSet)) cycle
             call AddHistMP(hMP_pSet(iSet), pPart, mom, 1.0/(mom**2), 1.0)
             call AddHistMP(hMP_ESet(iSet), pPart, mom0, 1.0/(mom0*mom), 1.0)
          end do

          ! fill Tmunu and Jmu:
          if (do_Tmunu) then
             call fillTmunu(TmunuNmu, pPart)
             if (do_Tmunu_pirho) then
                select case (pPart%ID)
                case (101)
                   call fillTmunu(arrTmunuNmu_hadr(1), pPart)
                case (103)
                   call fillTmunu(arrTmunuNmu_hadr(2), pPart)
                end select
             end if
             if (perEnsemble_Tmunu) call fillTmunu(arrTmunuNmu(i), pPart)
          end if

          if (do_cumulants) call calcCumulants1

       end do
    end do

    !
    ! produce output:
    !

    if (timestep>0) then ! this is a real timestep
       call doOutputTimestep
    else
       call doOutputFinal
    end if



  contains

    !**************************************************************************
    subroutine calcCumulants1

      real, dimension(3) :: scalePos
      real :: mVol
      integer :: ii,iVol

      ! scale the position to the box extends:
      scalePos(1:3) = abs(pPart%position(1:3))/gridsize(1:3)

      ! the minimal box volume with all coordinates included:
      mVol = maxval(scalePos)**3

      ! select minimal volume bin:
      iVol = int(mVol*10)+1

      if ((iVol<1).or.(iVol>10)) then
         write(*,*) 'mVol,iVol: ',mVol,iVol
         write(*,*) 'ooops!!!!!!!!!'
      end if

      if (pPart%charge /= 0) then
         countQ(iVol:10) = countQ(iVol:10) + pPart%charge
         countNch(iVol:10) = countNch(iVol:10) + 1.0
      end if

      if (pPart%ID==101) then
         if (pPart%charge==-1) then
            countNpi(iVol:10,1) = countNpi(iVol:10,1) + 1.0
         else if (pPart%charge==1) then
            countNpi(iVol:10,2) = countNpi(iVol:10,2) + 1.0
         end if
      end if

    end subroutine calcCumulants1

    !**************************************************************************
!!$    subroutine calcCumulants2
!!$
!!$      integer :: ii
!!$      real :: h
!!$
!!$      do ii=1,10
!!$         call AveragerAdd(AveQ(ii), countQ(ii))     ! Q = sum_i q_i
!!$         call AveragerAdd(AveNch(ii), countNch(ii)) ! Nch = sum_i (q_i/=0?1:0)
!!$         call AveragerAdd(AveNpi(ii,1), countNpi(ii,1)) ! = N(pi-)
!!$         call AveragerAdd(AveNpi(ii,2), countNpi(ii,2)) ! = N(pi+)
!!$
!!$         h = countNch(ii)
!!$         if (h > 0.0) then
!!$            h = countQ(ii)/h
!!$         else
!!$            h = 10000.0
!!$         end if
!!$         call AveragerAdd(AveF(ii), h ) ! = Q/Nch
!!$
!!$         h = sum(countNpi(ii,:))
!!$         if (h > 0.0) then
!!$            h = (countNpi(ii,2)-countNpi(ii,1))/h
!!$         else
!!$            h = 10000.0
!!$         end if
!!$         call AveragerAdd(AveFpi(ii), h ) ! = (N(pi+)-N(pi-))/(N(pi+)+N(pi+))
!!$
!!$         h = countNpi(ii,2)
!!$         if (h > 0.0) then
!!$            h = countNpi(ii,1)/h
!!$         else
!!$            h = 10000.0
!!$         end if
!!$         call AveragerAdd(AveR(ii), h ) ! = N(pi+)/N(pi-)
!!$
!!$         write(7000+ii,*) timestep,CountNpi(ii,1)
!!$         write(7100+ii,*) timestep,CountNpi(ii,2)
!!$
!!$         flush(7000+ii)
!!$         flush(7100+ii)
!!$
!!$      end do
!!$
!!$      do ii=1,9
!!$         call AveragerWrite(5000+ii,timestep*1.0,AveNpi(ii,1))
!!$         call AveragerWrite(5100+ii,timestep*1.0,AveNpi(ii,2))
!!$         call AveragerWrite(5200+ii,timestep*1.0,AveFpi(ii))
!!$         call AveragerWrite(5300+ii,timestep*1.0,AveR(ii))
!!$         flush(5000+ii)
!!$         flush(5100+ii)
!!$         flush(5200+ii)
!!$         flush(5300+ii)
!!$      end do
!!$
!!$      if (mod(timestep,10)==0) then
!!$         rewind(6001)
!!$         rewind(6002)
!!$         rewind(6003)
!!$         rewind(6004)
!!$         do ii=1,10
!!$            call AveragerWrite(6001,ii*0.1,AveNpi(ii,1))
!!$            call AveragerWrite(6002,ii*0.1,AveNpi(ii,2))
!!$            call AveragerWrite(6003,ii*0.1,AveFpi(ii))
!!$            call AveragerWrite(6004,ii*0.1,AveR(ii))
!!$         end do
!!$         flush(6001)
!!$         flush(6002)
!!$         flush(6003)
!!$         flush(6004)
!!$      end if
!!$
!!$    end subroutine calcCumulants2

    !**************************************************************************
    subroutine doOutputTimestep
      if (do_P) then
         if (mod(timestep,5)==1) then
            do iSet=1,nSet
               if (.not.useSet(iSet)) cycle

               call WriteHistMP(hMP_pSet(iSet), &
                    file='p_Set'//achar(48+iSet)//'_'//intTochar4(timestep)//'.dat', &
                    add=1e-20, mul=mulfak, iColumn=1)

               call WriteHistMP(hMP_ESet(iSet), &
                    file='E_Set'//achar(48+iSet)//'_'//intTochar4(timestep)//'.dat', &
                    add=1e-20, mul=mulfak, iColumn=1)

            end do

            call WriteHist(hMassRho, file='massRho_'//intTochar4(timestep)//'.dat', add=1e-20, mul=mulfak)
            ! call WriteParticleVector('parts_'//intTochar(timestep),realPart)
            ! call WriteHistMC(hMCMomPion, file='MomPion_'//intTochar4(timestep)//'.dat', add=1e-20, mul=mulfak)
         end if
      end if


      do iSet=1,nSet
         if (.not.useSet(iSet)) cycle

         open(123,file="BoxAnalysis_Mult_Set"//achar(48+iSet)//".dat",&
              status="old",position='append')
         write(123,'(i11,1P,100E12.4,0P)') timestep, &
              sum(hMP_ESet(iSet)%yVal(:,:,2),dim=2)*mulfak, &
              0. ! 0 for historical reasons
         close(123)
      end do


      if (do_Tmunu) then

         open(123,file="BoxAnalysis_Tmunu.dat",status="old",position='append')
         write(123,'(i11,1P,100E14.6,0P)') timestep, &
              & TmunuNmu%Tmunu(:)*mulfak, &
              & TmunuNmu%Nmu(:)*mulfak, &
              & TmunuNmu%Jmu(:)*mulfak
         close(123)


         if (do_Tmunu_pirho) then
            open(123,file="BoxAnalysis_Tmunu.pion.dat",status="old",position='append')
            write(123,'(i11,1P,100E14.6,0P)') timestep, &
                 & arrTmunuNmu_hadr(1)%Tmunu(:)*mulfak, &
                 & arrTmunuNmu_hadr(1)%Nmu(:)*mulfak, &
                 & arrTmunuNmu_hadr(1)%Jmu(:)*mulfak
            close(123)

            open(123,file="BoxAnalysis_Tmunu.rho.dat",status="old",position='append')
            write(123,'(i11,1P,100E14.6,0P)') timestep, &
                 & arrTmunuNmu_hadr(2)%Tmunu(:)*mulfak, &
                 & arrTmunuNmu_hadr(2)%Nmu(:)*mulfak, &
                 & arrTmunuNmu_hadr(2)%Jmu(:)*mulfak
            close(123)
         end if

         if (perEnsemble_Tmunu) then
            ! since we print all information for every ensemble, we must not divide by nEns here
            do i=1,nEns
               open(123,file="BoxAnalysis_Tmunu."//intToChar4(i)//".dat", status="old",position='append')
               write(123,'(i11,1P,100E12.4,0P)') timestep, &
                    & arrTmunuNmu(i)%Tmunu(:)/boxVol, &
                    & arrTmunuNmu(i)%Nmu(:)/boxVol, &
                    & arrTmunuNmu(i)%Jmu(:)/boxVol
               close(123)
            end do
         end if
      end if

      if (do_velrel) then
         ! calculate average of velrel
         call CalcAverageVelRel(realPart,timestep,nEns,nPart)
      end if

!      if (do_cumulants) call calcCumulants2

      call cR_Write(1)
    end subroutine doOutputTimestep

    !**************************************************************************
    subroutine doOutputFinal

      do iSet=1,nSet
         if (.not.useSet(iSet)) cycle

         open(123,file="BoxAnalysis_Final_Mult_Set"//achar(48+iSet)//".dat",&
              status="unknown")
         call WriteHistMP_Names(iSet,123)
         write(123,'(i11,1P,100E12.4,0P)') timestep, &
              sum(hMP_ESet(iSet)%yVal(:,:,2),dim=2)*mulfak, &
              0. ! 0 for historical reasons
         close(123)

         call WriteHistMP(hMP_pSet(iSet), &
              file='p_Set'//achar(48+iSet)//'_final.dat', &
              add=1e-20, mul=mulfak, iColumn=1)

         call WriteHistMP(hMP_ESet(iSet), &
              file='E_Set'//achar(48+iSet)//'_final.dat', &
              add=1e-20, mul=mulfak, iColumn=1)
      end do


    end subroutine doOutputFinal

  end subroutine DoBoxAnalysisTime

  !****************************************************************************
  !****s* BoxAnalysis/CalcAverageVelRel
  ! NAME
  ! subroutine CalcAverageVelRel(realPart,timestep,nEns,nPart)
  ! PURPOSE
  ! calculate the average v_rel of all particles with each other. Due to
  ! speed reasons, it may be a good idea to correlate only particles in the
  ! same ensemble, but this s only a approximation.
  !****************************************************************************
  subroutine CalcAverageVelRel(realPart,timestep,nEns,nPart)
    use particleDefinition
!    use lorentzTrafo, only: eval_sigmaBoost

    type(particle),dimension(:,:),intent(in), target :: realPart
    integer, intent(in) :: timestep
    integer, intent(in) :: nEns,nPart

    integer :: iEns1,iEns2, iPart1,iPart2
    type(particle), POINTER :: pPart1, pPart2
    real :: sum0,sum1,velrel,m1,m2,s
    real :: ptot(0:3)
!    real :: vrel
!    real, dimension(1:3) :: vrel_vector


    ! due to speed reasons, I only correlate particles in the same ensemble

!    write(*,*) 'calculating velrel....'
    sum0 = 0.0
    sum1 = 0.0
    do iEns1=1,nEns
       do iPart1=1,nPart
          pPart1 => realPart(iEns1,iPart1)
          if (pPart1%Id <  0) exit
          if (pPart1%Id <= 0) cycle
          m1 = pPart1%mass**2

          iEns2 = iEns1
          do iPart2=iPart1+1,nPart

             pPart2 => realPart(iEns2,iPart2)
             if (pPart2%Id <  0) exit
             if (pPart2%Id <= 0) cycle
             m2 = pPart2%mass**2

             ptot = pPart1%momentum + pPart2%momentum
             s = ptot(0)**2-sum(ptot(1:3)**2)
             velrel = sqrt( max(0.0,(s-m1-m2)**2/4-m1*m2) )/( pPart1%momentum(0)*pPart2%momentum(0) )

!             vrel_vector=pPart1%velocity-pPart2%velocity
!             vrel = sqrt(Dot_product(vrel_vector,vrel_vector))

!             write(*,*) velrel,vrel,velrel/vrel,eval_sigmaBoost(pPart1%momentum,pPart2%momentum)

!             write(*,*) (velrel-vrel*eval_sigmaBoost(pPart1%momentum,pPart2%momentum))/velrel

             sum0 = sum0 + 1
             sum1 = sum1 + velrel

          end do
       end do
    end do

    write(*,'(A,f8.5,f15.0,i13)') 'Average velrel: ',sum1/sum0, sum0, timestep

  end subroutine CalcAverageVelRel

  !****************************************************************************
  !****s* BoxAnalysis/readInput
  ! NAME
  ! subroutine readInput
  !
  ! PURPOSE
  ! Reads input in jobcard out of namelist "BoxAnalysis"
  !****************************************************************************
  subroutine readInput

    use output, only: Write_ReadingInput

    !**************************************************************************
    !****n* BoxAnalysis/BoxAnalysis
    ! NAME
    ! NAMELIST BoxAnalysis
    ! PURPOSE
    ! Includes the switches:
    ! * do_Tmunu
    ! * perEnsemble_Tmunu
    ! * do_P
    ! * do_velrel
    ! * do_Cumulants
    ! * useSet
    !*************************************************************************
    NAMELIST /BoxAnalysis/ &
         do_Tmunu, do_P, do_velrel, perEnsemble_Tmunu, &
         do_cumulants, useSet
    integer :: ios

    useSet = (/ .false., .false., .false., .true., .true. /)

    call Write_ReadingInput('BoxAnalysis',0)
    rewind(5)
    read(5,nml=BoxAnalysis,IOSTAT=ios)
    call Write_ReadingInput('BoxAnalysis',0,ios)

    write(*,*) '  do Tmunu: ',do_Tmunu,'   perEnsemble: ',perEnsemble_Tmunu
    write(*,*) '  do P:     ',do_P
    write(*,*) '  do velrel:',do_velrel
    write(*,*) '  do cumlants: ',do_cumulants
    write(*,*) '  use MP set : ',useSet

    call Write_ReadingInput('BoxAnalysis',1)
  end subroutine readInput

end module BoxAnalysis
