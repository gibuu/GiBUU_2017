program CalcPT

! run this e.g. with: ... < ~/CHI/QC/BUU/NEW/WASTI/CSC2/N_JLAB5/00/14559/HiLep.NuQ2planeXS.000.dat

  implicit none

  integer :: iQ2,iNu,iZh,iPT,i
  real,dimension(5) :: v
  character*(100) :: BUF
  real :: pT,mP,mX, W
  real, dimension(0:2,100) :: S
  real, dimension(2,100) :: Sx
  real, dimension(0:11,200) :: SS
  real :: SS0


  mp = 0.938
!  mX = 0.938
!  mX = 1.232 ! Delta
  mX = 1.076 ! N+Pi
!  mX = 2.0

  S = 0.0
  Sx(1,:) = 99.9
  Sx(2,:) = -99.9

  SS0 = 0.0
  SS = 0.0

!!$  do iZh=1,100
!!$     i = int((iZh*0.01)*10)
!!$     write(*,*) izH,iZH*0.01,i
!!$  end do

  v(1) = 2.300! Q2
  v(2) = 3.83019! nu
  w = mP**2-v(1)+2*mP*v(2)
  w = sqrt(w)
  do iZh=1,120
     pT = CalcPT_1(W, v(2),v(1),iZh*0.01)
     write(*,*) iZh*0.01,pT
  enddo
  stop
     

  do iQ2=1,51
     do iNu=1,27
        read(*,*) v
        if (v(3).lt.1e-5) cycle


        w = mP**2-v(1)+2*mP*v(2)
        if (w.lt.0.0) cycle
        w = sqrt(w)


        do iZh=1,100
           pT = CalcPT_1(W, v(2),v(1),iZh*0.01)

           do ipT=1,200
              if (ipT*0.02 .lt. pT**2) then
                 SS0 = SS0 + v(3)
                 i = int((iZh*0.01)*10)
                 SS(i,iPT) = SS(i,iPT)+v(3)
              endif
           enddo

           

           S(0,iZh) =  S(0,iZh)+v(3)
           S(1,iZh) =  S(1,iZh)+v(3)*pT
           S(2,iZh) =  S(2,iZh)+v(3)*pT**2

           if (pT.lt.Sx(1,iZh)) Sx(1,iZh) = pT
           if (pT.gt.Sx(2,iZh)) Sx(2,iZh) = pT

        enddo


     end do
     read(*,'(A)') BUF

  end do

  do iZh=1,100
     write(*,'(3f12.5)') iZh*0.01,S(1,iZh)/(S(0,iZh)+1e-20),&
          & sqrt(S(2,iZh)/(S(0,iZh)+1e-20)-(S(1,iZh)/(S(0,iZh)+1e-20))**2)

     write(11,'(3f12.5)') iZh*0.01,CalcPT_1(2.1978, 3.5211, 1.7265,iZh*0.01)

     write(12,'(3f12.5)') iZh*0.01,Sx(1,iZh),Sx(2,iZh)

  enddo

  do iZH=0,9
     do ipT=1,200
        write(101+iZH,'(1f12.5,12e12.3)') iPT*0.02, SS(iZH,iPT)/SS0+1e-20
     enddo
  enddo


  do ipT=1,200
     write(13,'(1f12.5,12e12.3)') iPT*0.02, SS(0:9,iPT)/SS0
  enddo


contains
  real function CalcPT_1(W,nu,Q2,zH)
    implicit none
    real, intent(in) :: W,nu,Q2,zH

    real :: x, a,b,c

    a = (w**2-mX**2)/(2*w)
    b=  a * (mP+nu)/(w*nu)
    c=  a * sqrt(nu**2+Q2)/(w*nu)

    x= (zH-b)/c
    
    CalcPT_1 = 0.0
    if (abs(x).gt.1.0) return
    CalcPT_1 = a*sqrt(1-x**2)


  end function CalcPT_1


end program CalcPT
