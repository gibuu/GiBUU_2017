program DoPiN
  use output
  use particleDefinition
  use particleProperties
  use Coll_Fritiof
  use Coll_Pythia
  use CollTools
  !use rhoMassParameter, only : srtFreeRhoMass
  use hadronFormation, only : forceInitFormation
  use random
  use histMPf90
  use version

  implicit none

  COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
  integer N,NPAD,K
  double precision P,V
  SAVE /PYJETS/
  
  COMMON/PYINT1/MINT(400),VINT(400)
  integer MINT
  double precision VINT
  SAVE /PYINT1/
  
  COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
  integer MSTP,MSTI
  double precision PARP,PARI
  SAVE /PYPARS/
  
  logical GetSwitchPythiaHermes

  integer :: i
  real :: r


  integer iEV,nEV, nEv10, nEv100
  !      parameter (NEV=10000000) !
!  parameter (NEV=1000000) ! 9m 12s (w=4), 4m 50s (w=0)
  parameter (NEV=100000) ! 30s
  !    parameter (NEV=1000)

  type(histogramMP),save :: HMP_pT

  integer iH
  double precision hhh

  nEv100 = nEv / 100
  nEv10 = nEv / 10

  call PrintVersion

  call forceInitFormation
  call InitParticleProperties


  !--------------------
  !... Initializations:
  !--------------------

  call SetSomeDefaults_PY

  call PYGIVE('MSTU(12) =0')

!!$      call PYGIVE('MSEL=1')

      call PYGIVE('MSEL=2')
!!$
!!$      call PYGIVE('MSEL=0')
!!$      call PYGIVE('MSUB(11)=0') !  11 f + f' -> f + f' (QCD)     
!!$      call PYGIVE('MSUB(12)=0') !  12 f + fbar -> f' + fbar'     
!!$      call PYGIVE('MSUB(13)=0') !  13 f + fbar -> g + g          
!!$      call PYGIVE('MSUB(28)=0') !  28 f + g -> f + g             
!!$      call PYGIVE('MSUB(53)=0') !  53 g + g -> f + fbar          
!!$      call PYGIVE('MSUB(68)=0') !  68 g + g -> g + g             
!!$      call PYGIVE('MSUB(91)=1') !  91 Elastic scattering         
!!$      call PYGIVE('MSUB(92)=1') !  92 Single diffractive (XB)    
!!$      call PYGIVE('MSUB(93)=1') !  93 Single diffractive (AX)    
!!$      call PYGIVE('MSUB(94)=1') !  94 Double  diffractive        
!!$      call PYGIVE('MSUB(95)=1') !  95 Low-pT scattering      
!!$      call PYGIVE('MSUB(96)=0') !  96 Semihard QCD 2 -> 2




!      call PYGIVE('MSTP(111) = 0') ! master switch fragmentation/decay

!      call PYGIVE('MSTP(91) = 0') ! switch 'primord. kT'
!      call PYGIVE('MSTP(61) = 0') ! switch 'ini.state.radiation'
!      call PYGIVE('MSTP(71) = 0') ! switch 'fin.state.radiation'

!      call PYGIVE('PARP(91)=2.5') ! <kT>: 62->1.05, 200->2.5
!      call PYGIVE('PARP(93)=10') ! upper cut off of kT-distr

!      call PYGIVE('CKIN(3)=1.5d0') ! pThat_min
      call PYGIVE('CKIN(3)=0.0d0') ! pThat_min
!      call PYGIVE('CKIN(4)=1.5d0') ! pThat_max

      call PYGIVE('PARP(2) = 1d0') ! lowest c.m. energy

!      call PYGIVE('MSTJ(21)=0') ! particle decay on/off
 
                               ! forbid the decay of:
      call PYGIVE('MDCY(C111,1)=0') ! pi0
      call PYGIVE('MDCY(C311,1)=0') ! K0
      call PYGIVE('MDCY(C310,1)=0') ! K_S0

      call PYGIVE('PARP(111) = 0d0') ! otherwise problems at low sqrts


!      call PYGIVE('CKIN(3)=1.5d0') ! pThat_min
!      call PYGIVE('MSTP(142) = 1') ! use reweighted events


      VINT(48) = 2d0 
      call PYEVWT(hhh)
      write(*,*) '>>> if Reweighting, then according  pThat^',&
           & 2*log(hhh)/log(2d0)



!!$      CALL PYINIT('CMS','p','n', 200d0)

!!$      CALL PYINIT('FIXT','pi-','p',50.00000D0)
!!$      CALL PYINIT('FIXT','pi-','p',62.94627D0)
!!$      CALL PYINIT('FIXT','pi-','p',79.24466D0)   
!!$      CALL PYINIT('FIXT','pi-','p',99.76311D0)   
!!$      CALL PYINIT('FIXT','pi-','p',125.5943D0)   
!!$      CALL PYINIT('FIXT','pi-','p',158.1139D0)  
!!$      CALL PYINIT('FIXT','pi-','p',199.0536D0)   
!!$      CALL PYINIT('FIXT','pi-','p',250.5936D0)   
!!$      CALL PYINIT('FIXT','pi-','p',315.4787D0)   
!!$      CALL PYINIT('FIXT','pi-','p',397.1642D0)   
!!$      CALL PYINIT('FIXT','pi-','p',500.0000D0)
!!$
!!$      CALL PYINIT('FIXT','pi+','p',50.00000D0)
!!$      CALL PYINIT('FIXT','pi+','p',62.94627D0)
!!$      CALL PYINIT('FIXT','pi+','p',79.24466D0)   
!!$      CALL PYINIT('FIXT','pi+','p',99.76311D0)   
!!$      CALL PYINIT('FIXT','pi+','p',125.5943D0)   
!!$      CALL PYINIT('FIXT','pi+','p',158.1139D0)  
!!$      CALL PYINIT('FIXT','pi+','p',199.0536D0)   
!!$      CALL PYINIT('FIXT','pi+','p',250.5936D0)   
!!$      CALL PYINIT('FIXT','pi+','p',315.4787D0)   
!!$      CALL PYINIT('FIXT','pi+','p',397.1642D0)   
!!$      CALL PYINIT('FIXT','pi+','p',500.0000D0)
!!$
!!$      CALL PYINIT('FIXT','pi+','n',50.00000D0)
!!$      CALL PYINIT('FIXT','pi+','n',62.94627D0)
!!$      CALL PYINIT('FIXT','pi+','n',79.24466D0)   
!!$      CALL PYINIT('FIXT','pi+','n',99.76311D0)   
!!$      CALL PYINIT('FIXT','pi+','n',125.5943D0)   
!!$      CALL PYINIT('FIXT','pi+','n',158.1139D0)  
!!$      CALL PYINIT('FIXT','pi+','n',199.0536D0)   
!!$      CALL PYINIT('FIXT','pi+','n',250.5936D0)   
!!$      CALL PYINIT('FIXT','pi+','n',315.4787D0)   
!!$      CALL PYINIT('FIXT','pi+','n',397.1642D0)   
!!$      CALL PYINIT('FIXT','pi+','n',500.0000D0)

      call PYGIVE('PARP(104)=0.30')            ! Min.energy for XS definition


      call PYGIVE('MSTP(122) =1')

!      CALL PYINIT('FIXT','pi-','p',515.0000D0)
!      srtFreeRhoMass = 31.102d0

!      CALL PYINIT('FIXT','p','p',3.809D0)
      CALL PYINIT('FIXT','p','p',1.809D0)

      write(*,*) 'UseHermes:',GetSwitchPythiaHermes()

      call CreateHistMP(HMP_pT, 'pT', 0.,10.,0.2, 2)

      !--------------
      !... MAIN LOOP:
      !--------------

      do iEV=1,NEV
         CALL WriteStatusMC(iEv,nEv,nEv100,nEv10)
   
         call GetJetsetVecINIT

         CALL PYEVNT

!         call PYLIST(2)
!         call PYGIVE('MSTI(1)=')

!         stop
!         call AddPartons(HMP_pT,PARI(10))
!         stop
      end do

      !-------------------
      !... Final Writeout:
      !-------------------

      write(6,*) 'PARI(1)=',PARI(1),'  <-- XS'
      write(6,*) 'MSTI(5)=',MSTI(5),'  <-- #(generated events)'
      write(6,*) 'PARI(2)=',PARI(2),' =?= ',PARI(1)/MSTI(5)
      write(6,*) 

      call PYSTAT(1)


      rewind(201)
      call WriteHistMP(HMP_pT,201,mul=PARI(2) ,add=1e-20)

      write(*,*) 'Histograms written.'

contains
  
  subroutine WriteStatusMC(iEV, nEV, nEV_100, nEv_10)
    IMPLICIT NONE
    integer iEV, nEV, nEV_100, nEv_10
    
    
    integer n_Perc
    
    if (nEV.lt.100) return  ! dont write status
    if (mod(iEV,nEV_100) .eq. 0) then
       n_Perc = iEV/nEV_100
       
!         write(*,*) 'DoPiN : ',n_Perc,'%'

       open(12,file='DoPiN.run',status='unknown')
       write(12,*) 'DoPiN.run : ',n_Perc,'%'
       close(12)
        
!!$         if (mod(iEV,nEV_10) .eq. 0) then
!!$            call ReportTotal
!!$            open(12,file='CHQLHL.run.rep',status='unknown')
!!$            call DumpHistos(12)
!!$            close(12)
!!$
!!$            call DumpRANDOM
!!$         endif

    endif

  end subroutine WriteStatusMC

  subroutine AddPartons(H,w)
    use histMPf90

    implicit none
    type(histogramMP) :: H
    real              :: w

    COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
    integer N,NPAD,K
    double precision P,V
    SAVE /PYJETS/

    integer i,iH
    double precision beta,fak,pT,yCM,yBin
    
    yBin = 0.75d0
    
    beta = -(P(1,3)+P(2,3))/(P(1,4)+P(2,4))
    fak = log((1+beta)/(1-beta))

!    write(*,*) 'fak=',fak
    
    do i=1,N
       if (K(i,1).gt.10) cycle ! not stable
       if (abs(K(i,2)).lt.100) cycle

       

       yCM = (fak+log((P(i,4)+P(i,3))/(P(i,4)-P(i,3))))/2
       !         write(*,*) i,yCM
       pT = sqrt(P(i,1)**2+P(i,2)**2)

!       write(*,*) i,ycm,pt

       if (abs(yCM).gt.yBin) cycle ! not in y-Bin

       if ((K(i,2).eq.111) .and. (pT > 3.5d0)) then
          write(*,*) 'particle ',i,yCM,pT
          call PYLIST(2)
       endif

       call AddHistMP(H,K(i,2),pT, w/(pT*(2*yBin)*6.28))  ! E ds/d^3p !!!

!       stop
    end do


  end subroutine AddPartons
    


end program DoPiN
