!******************************************************************************
! program to print out the Bosted
!
! run it with:
! ./testBosted.x < jobBosted
!******************************************************************************
program testBosted
  use ParamEP
  use inputGeneral
  use particleProperties
  use output

  implicit none

  real :: W,Q2,eps,XS,Gamma
  real :: Ein,theta,thetadeg,nu
  real, save :: pi = 3.14159265
  integer :: i,ios    
  
  NAMELIST /epkin/ Ein,thetadeg

  call readinputGeneral
  call initParticleProperties
  
  rewind(5,IOSTAT=ios)
  read(5,nml=epkin,IOSTAT=ios)
  call Write_ReadingInput("epkin",0,ios) 
!  Ein = 5.498
!  thetadeg = 20.45
  
  theta = thetadeg * pi/180. 
  write(119,"('#',1X,'Ein =',f6.3,3X,'theta(deg) =',f6.2,3X,'theta(rad) =',f6.4,/)") Ein, thetadeg, theta
  
  Q2=1
  eps=0.5
  call CalcParamEP(1.5,Q2,eps, XS) ! dummy in order to read input
  
  
  
  write(119,"('#',4X,'W',9X,'Q2',8X,'nu',8X,'eps',8X,'Gamma',8X,'XS')")
   
  do i=110,300
     W = i*0.01
     
     nu = Ein - Ef(W,Ein,theta)
     Q2 = Qsq(W,Ein,theta) 
     eps = (1 + 2*qvec2(W,Ein,theta)/Q2 * tan(theta/2.)**2)**(-1)
     
     call CalcParamEP(W,Q2,eps, XS) 
     
! Now multiply sum of longitudinal and transverse cross section with factor Gamma
! to get dd X-section      
     
     Gamma = (nu - Q2/(2*0.93822))/(2*pi**2*Q2) * (Ein - nu)/Ein &
           & * 1./(1 - eps) * 1./137.
     XS = XS * Gamma
                
! XS is in mubar/(sr * GeV)        
     write(119,'(4f10.4,2e14.4)') W,Q2,nu,eps, Gamma, XS
  end do
  
CONTAINS

REAL FUNCTION Qsq(W,Ein,theta)
Real W, Ein,theta

Qsq = 4 * Ein *Ef(W,Ein,theta) * sin(theta/2)**2

END FUNCTION Qsq

REAL FUNCTION Ef(W,Ein,theta)
Real W,Ein,theta,M

M = 0.93822

Ef = (M**2 + 2*M*Ein - W**2)/(2*M + 4*Ein*sin(theta/2.)**2 )

END FUNCTION Ef

REAL FUNCTION qvec2(W,Ein,theta)
Real W,Ein,theta,me

me = 0.511

qvec2 = Ef(W,Ein,theta)**2 + Ein**2 - 2*me**2 - 2*Sqrt(Ef(W,Ein,theta)**2 - me**2) &
      &  * Sqrt(Ein**2 - me**2) * cos(theta)
      
End FUNCTION qvec2      
  
end program testBosted
