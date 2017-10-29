

program test
  use helicityAmplitudes
  use particleProperties, only: initParticleProperties

!  call writeOutFORMFACTOR()
  call initParticleProperties
  call cross_sections
end program test


subroutine cross_sections
  use helicityAmplitudes
  use idTable
  use particleProperties, only: hadron
  use baryonWidthMedium, only: WidthBaryonMedium, partialWidthBaryonMedium
  use mediumDefinition, only: medium, vacuum
  USE minkowski, only: SP
  use constants, only: mN
  implicit none
  integer, parameter :: targetCharge=0
  real,parameter :: QSquared=0.
  real :: W
  real,parameter :: W_max=4.
  real,parameter :: W_min=0.
  real,parameter :: W_delta=0.01
  real :: A1,A3,S1
  integer :: resonanceID
  real,dimension(nucleon+1:F37_1950) :: sigma_T,sigma_L
  real,dimension(nucleon+1:F37_1950) :: sigma_T_pi,sigma_L_pi
  real :: sigma_T_total,sigma_L_total
  real :: sigma_T_total_pi,sigma_L_total_pi
  real, dimension(0:3) :: q,pRes,pres_pole
  integer :: i
  real :: W_dependence_pi,W_dependence,width,widthOnPole, gammaNPi, gamma_NGamma, k_W, k_R

  
  W_loop : do i=0,nint((W_max-W_min)/W_delta)
     W=W_min+float(i)*W_delta
     q=0
     q(0)=(W**2+QSquared-mN**2)/2./mN
     if(q(0)**2+QSquared.lt.0) stop 'weird error'
     q(3)=sqrt(q(0)**2+QSquared)
     pRes(1:3)=q(1:3)
     pRes(0)=sqrt(W**2+Dot_Product(pRes(1:3),pRes(1:3)))

     write(31,*) W, QSquared, SP(q,q),sqrt(SP(pres,pres))

     resonance_loop: Do resonanceId=nucleon+1,F37_1950
        call get_helicityAmplitudes(targetCharge,resonanceId,QSquared,A1,A3,S1)
        write(21,*) A1, A3,S1
        pRes_pole(1:2)=0
        pRes_pole(3)=sqrt(((W**2-QSquared-mN**2)/2./mN)**2-QSquared)
        pRes_pole(0)=sqrt(hadron(resonanceID)%mass**2+Dot_Product(pRes(1:3),pRes(1:3)))

        width=WidthBaryonMedium(resonanceID,W,pRes,vacuum)
        gammaNPi=partialWidthBaryonMedium(resonanceID,W,.false.,pion,nucleon,pres,vacuum)

        k_W=(W**2-mN**2)/2./W
        k_R=(hadron(resonanceID)%mass**2-mN**2)/2./hadron(resonanceID)%mass
        gamma_NGamma=(k_W/k_R)**getN(resonanceID)*((0.5**2+k_R**2)/(0.5**2+k_W**2))
        gamma_NGamma=1


        W_dependence_pi=W**2*width*gammaNPi/((W**2-hadron(resonanceID)%mass**2)**2+W**2*width**2)*gamma_NGamma
        W_dependence=W**2*width**2/((W**2-hadron(resonanceID)%mass**2)**2+W**2*width**2)*gamma_NGamma
        
        widthOnPole=WidthBaryonMedium(resonanceID,hadron(resonanceId)%mass,pRes,vacuum)
        sigma_T(resonanceID)=2.*mN/(hadron(resonanceId)%mass*widthOnPole)*(A1**2+A3**2)
        sigma_L(resonanceID)=4.*mN/(hadron(resonanceId)%mass*widthOnPole)& 
& *QSQuared/Dot_Product(q(1:3),q(1:3))*S1**2

        sigma_T(resonanceID)=sigma_T(resonanceID)*W_dependence
        sigma_L(resonanceID)=sigma_L(resonanceID)*W_dependence
        sigma_T_pi(resonanceID)=sigma_T(resonanceID)*W_dependence_pi
        sigma_L_pi(resonanceID)=sigma_L(resonanceID)*W_dependence_pi
        
     end do resonance_loop
     ! Convert :  GeV^-2=0.197**2 fm^2=0.197**2 10-30 m^2 = 0.197**2* 10mb=0.197**2 *10^4 muB
     sigma_T=sigma_T*0.197**2*1E4
     sigma_L=sigma_L*0.197**2*1E4
     sigma_T_pi=sigma_T_pi*0.197**2*1E4
     sigma_L_pi=sigma_L_pi*0.197**2*1E4

     sigma_T_total=sum(sigma_T)
     sigma_L_total=sum(sigma_L)
     sigma_T_total_pi=sum(sigma_T_pi)
     sigma_L_total_pi=sum(sigma_L_pi)



     write(12,'(40E15.4)')   W,q(0),sigma_T_total, sigma_L_total,sigma_T_total_pi, sigma_L_total_pi  
     write(11,'(40E15.4)')   W,q(0),sigma_T
     write(*,*)   W,sigma_T_total, sigma_L_total, sigma_L
  end do W_loop

  contains 
    integer function getN(ID)
      use idtable
      implicit none
      integer, intent(in) :: ID
      SELECT CASE(ID)
      case(Delta)
         getN=2
      case(P11_1440)
         getN=1
      case(D13_1520)
         getN=2
      case(S11_1535)
         getN=3
      case(S11_1650) 
         getN=4
      case(F15_1680)
         getN=3
      case(D33_1700)
         getN=4
      case DEFAULT
         getN=1
      end select

 end function getN

  
end subroutine cross_sections
