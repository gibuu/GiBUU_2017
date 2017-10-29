program createDataTable
  use ParamEP
  use inputGeneral
  use particleProperties, only: initParticleProperties
  use particleDefinition
  !use Coll_gammaN_low
  use output
  use constants
  use CollTools

  implicit none

  integer, parameter :: nPart=3, nTypes=3

  real :: W,Q2,eps,XSp,XS,nu,XS2,XS3,Gamma
  integer :: i,j,iq,iiq
  real,          dimension(0:nTypes)       :: XS_s
  real ::  rMC, fT, xB

  type(particle)   :: inPart               ! incoming nucleon
  type(particle),dimension(10) :: outPart  ! outgoing particles

  real :: XPQ(-25:25)
  integer,save :: nMC = 10000

  NAMELIST /datatable/ eps,nMC

  call readinputGeneral
  call initParticleProperties
  call SetSomeDefaults_PY

  call CalcParamEP(1.5,1.0,0.5, XS) ! dummy in order to read input

  inPart%ID    = 1
  inPart%Charge= 1
  inPart%mass  = 0.938
  inPart%momentum(0) = 0.938

  call WriteParticle(6,1,1,inpart)

  call PYGIVE("MSTP(51)=")

  ! step 0: read input

  eps = 0.95
  
  call Write_ReadingInput('datatable',0)
  rewind(5)
  read(5,nml=datatable)
  write(*,*) 'eps=',eps
  write(*,*) 'nMC=',nMC
  call Write_ReadingInput('datatable',1)

  ! step 1: Q2=0

  Q2 = 0.0
  do i=1,52
     nu = 0.12+i*0.03
     W = sqrt(0.938**2+2*0.938*nu)
     call CalcParamEP(W,Q2,0.0, XS)
     write(119,'(4f12.4,1P,e12.4)') W,nu,Q2,0.0, XS 
  end do

  ! F2 parametrizations

!  call GenerateF2
!  stop


  ! calculate DIS cross section

!  call GenerateDIS
!  stop


  ! step 2: Param-Data-Table
  ! step 3: MC-Data-Table

  do j=1,45
!  do j=22,45

     Q2=j*0.1

     do i=110,200,2
     !  do i=190,200,2
        W = real(i)*0.01
        if (W.ge.2.0) W = 1.99

        nu = (W**2+Q2-0.938**2)/(2*0.938)

        !        write(*,*) 

        call CalcParamEP(W,Q2,eps, XSp)
        Gamma = Flux_Bosted(W,Q2,eps)
        write(121,'(4f12.4,1P,99e12.4)') W,nu,Q2,eps, XSp, Gamma

        rMC = 0.0
        XS_s = 0.0
!        call CalcXS_MC(inPart,W,Q2,eps,nMC, XS_s,rMC)
        write(122,'(4f12.4,1P,3e12.4,0P,f12.4,1P,99e12.4)') Q2,nu,W,eps, XSp, fT, XS_s(0), rMC, XS_s(1:3), Gamma
        write(6,'(A,4f12.4,1P,3e12.4,0P,f12.4,1P,99e12.4)') '==>',Q2,nu,W,eps, XSp, fT, XS_s(0), rMC, XS_s(1:3),Gamma


     end do
     write(*,*)
     write(121,*)
     write(122,*)

  end do




contains
  subroutine GenerateF2


    do j=1,45
       Q2=j*0.1 

       do i=174,300,2
          !  do i=190,200,2
          W = real(i)*0.01
          !        if (W.eq.2.0) W = 1.99

          xB = Q2/(W**2-0.938**2+Q2)
          nu = (W**2+Q2-0.938**2)/(2*0.938)

          call PYPDFU(2212,xB,Q2,XPQ)

          rMC = 0.0
          iiq = 1
          do iq=1,6
             rMC = rMC+ (iiq/3.0)**2 * xpq(iq)
             iiq = 3-iiq
          end do
          rMC = rMC * 4 * pi**2 * alphaQED/(Q2*(1-xB)) / GeVSquared_times_mb * 1000 ! now in mub
          rMC = rMC * (1+(4*0.938**2*Q2)/(Q2+W**2-0.938**2)**2) ! F2 ~= sigmaT+sigmaL = sigmaTot

          call CalcParamEP_ALLM(W,Q2, XS)
          call CalcParamEP_ALLM97(W,Q2, XS2)

          call CalcParamEP(W,Q2,0.999, XS3)

          write(6,'(A,4f12.4,1P,10e12.4)') '==>',Q2,nu,W,xB, rMC,XS,XS2,XS3

          write(124,'(4f12.4,1P,10e12.4)') Q2,nu,W,xB, rMC,XS,XS2,XS3
       end do

       write(*,*)
       write(124,*)
    end do

  end subroutine GenerateF2


  subroutine GenerateDIS
    do j=1,45
       !  do j=22,45

       Q2=j*0.1

       do i=110,200,2
          !  do i=190,200,2
          W = real(i)*0.01
          if (W.ge.2.0) W = 1.99

          xB = Q2/(W**2-0.938**2+Q2)
          nu = (W**2+Q2-0.938**2)/(2*0.938)

          call PYPDFU(2212,xB,Q2,XPQ)

          rMC = 0.0
          iiq = 1
          do iq=1,6
             rMC = rMC+ (iiq/3.0)**2 * xpq(iq)
             iiq = 3-iiq
          end do
          rMC = rMC * 4 * pi**2 * alphaQED/(Q2*(1-xB)) / GeVSquared_times_mb * 1000 ! now in mub
          rMC = rMC * Q2**2/(Q2+0.77**2)**2

          !        rMC = rMC*(1.0 - (W**2/(W**2+Q2))**3) ! additional factor

          write(123,'(4f12.4,1P,3e12.4,0P,f12.4,1P,10e12.4)') Q2,nu,W,0.0, rMC, xB
          write(6,'(A,4f12.4,1P,3e12.4,0P,f12.4,1P,10e12.4)') '==>',Q2,nu,W,0.0, rMC, xB

       end do
       write(*,*)
       write(123,*)
    end do

  end subroutine GenerateDIS


!!$  subroutine CalcXS_MC
!!$
!!$    use constants, only : pi
!!$
!!$    use PyVP, only : CalcFlux
!!$
!!$    real :: E_li,E_lf,theta_lf
!!$    integer :: iMC,i
!!$
!!$    type(particle),dimension(nTypes,nPart) :: outPart_s
!!$    logical,       dimension(nTypes)       :: flagOK_s
!!$    !    real,          dimension(nTypes)       :: XS_s
!!$
!!$    real,dimension(2) :: XS2
!!$    real :: x,y,dum
!!$
!!$    call CalcElectronFromPhoton(W,Q2,eps, E_li,E_lf,theta_lf)
!!$    x = Q2/(2*0.938*nu)
!!$    y = nu/E_li
!!$    call CalcFlux(E_li,Q2,y,x, fT,dum)
!!$    fT = fT / (1e3* pi/(E_lf*E_li))
!!$
!!$    !    write(*,'(1P,7e14.3)') W,Q2,eps,E_li,E_lf,theta_lf,fT
!!$
!!$    outPart(:)%ID = 0
!!$    do i=1,nTypes
!!$       outPart_s(i,:) = outPart(1:nPart) ! save presets
!!$    enddo
!!$
!!$
!!$    XS_s = 0
!!$
!!$    ! ==== 1: RESONANCES ====
!!$
!!$
!!$    call GenEvent_Resonance(inPart,outPart_s(1,:),flagOK_s(1),XS_s(1), &
!!$         & E_li,E_lf,theta_lf)
!!$
!!$    !    write(*,*) "1:",flagOK_s(1),XS_s(1)
!!$    !    stop
!!$
!!$    ! ==== 2: 1Pi Background ====
!!$
!!$    XS2 = 0.0
!!$    iMC = 0
!!$    do 
!!$       XS_s(2) = 0.
!!$       call GenEvent_1Pi_back( inPart,outPart_s(2,:),flagOK_s(2),XS_s(2), &
!!$            & E_li,E_lf,theta_lf, W,Q2,eps)
!!$       !       write(*,*) "2:",flagOK_s(2),XS_s(2)
!!$       if (flagOK_s(2)) XS2=XS2+(/ 1.0,XS_s(2)/)
!!$
!!$       iMC = iMC+1
!!$
!!$       if (iMC.eq.10*nMC) exit
!!$       if (XS2(1).eq.nMC) exit
!!$    end do
!!$
!!$    rMC = XS2(1)/iMC
!!$
!!$    if (XS2(1).gt.0.0) XS2(2)=XS2(2)/XS2(1) * 4*pi
!!$
!!$
!!$    XS_s(2) = XS2(2)
!!$
!!$
!!$    ! ==== 3: 2Pi Background ====
!!$    call GenEvent_2Pi_back( inPart,outPart_s(3,:),flagOK_s(3),XS_s(3), &
!!$         & E_li,E_lf,theta_lf, W,Q2,eps, fT)
!!$
!!$    ! ==== RESCALE ALL CROSS SECTIONS ===
!!$
!!$    !... dividing by flux in order to get sigma^{gamma*} (in mub):
!!$    XS_s = XS_s / fT
!!$
!!$    !    write(*,*) ">>",XS_s
!!$    XS =  sum(    XS_s )
!!$
!!$
!!$  end subroutine CalcXS_MC

end program createDataTable
