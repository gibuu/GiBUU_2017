
Program testNucD
  ! test NucD module
  use NucD
  implicit double precision (A-H,O-Z)
  !c
  real :: rd(3),ad(3),rhod(3),RMS(3),RCHRG(3)
  !c
  integer :: Anuc, Znuc
  Anuc=103
  Znuc=45
  !c
  call DFS(Anuc,Znuc,rd,ad,rhod,RMS,RCHRG)
  !c
  dRMS=RMS(2)-RMS(1)
  write(6,6001)RMS(2),RMS(1),dRMS
6001 format(//28x,'Root-Mean-Square Radii and Neutron Skin:'/&
       &12x,' sqrt{<r^2>}_n : ',f10.4,' [fm]'/                    &
       &12x,' sqrt{<r^2>}_p : ',f10.4,' [fm]'/                    &
       &12x,' diff(n - p)   : ',f10.4,' [fm]')
!  write(6,6005)(Rqq(i,Anuc,Znuc),aqq(i,Anuc,Znuc),i=1,2)
!6005 format(//12x,' Rp : ',F10.4,' [fm]',3x,'ap : ',F10.4,' [fm]'/&
!       &12x,' Rn : ',F10.4,' [fm]',3x,'an : ',F10.4,' [fm]')  
  !c

!  call DFS(Anuc,Znuc,rd,ad,rhod,RMS,RCHRG)
!  write(*,*) Anuc,Znuc, rd, ad, rhod


  call DFS(16,8,rd,ad,rhod,RMS,RCHRG)
  write(111,'("A=",I3,"Z=",I3,3("&",2F9.3))') 16,8, rhod(1:2),rd(1:2), ad(1:2)
  call DFS(40,20,rd,ad,rhod,RMS,RCHRG)
  write(111,'("A=",I3,"Z=",I3,3("&",2F9.3))') 40,20, rhod(1:2),rd(1:2), ad(1:2)
  call DFS(56,26,rd,ad,rhod,RMS,RCHRG)
  write(111,'("A=",I3,"Z=",I3,3("&",2F9.3))') 56,26, rhod(1:2),rd(1:2), ad(1:2)
  call DFS(103,45,rd,ad,rhod,RMS,RCHRG)
  write(111,'("A=",I3,"Z=",I3,3("&",2F9.3))') 103,45, rhod(1:2),rd(1:2), ad(1:2)
  call DFS(197,79,rd,ad,rhod,RMS,RCHRG)
  write(111,'("A=",I3,"Z=",I3,3("&",2F9.3))') 197,79, rhod(1:2),rd(1:2), ad(1:2)
  call DFS(207,82,rd,ad,rhod,RMS,RCHRG)
  write(111,'("A=",I3,"Z=",I3,3("&",2F9.3))') 207,82, rhod(1:2),rd(1:2), ad(1:2)


  stop
end Program testNucD
