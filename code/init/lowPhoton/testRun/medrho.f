      subroutine check_med
      implicit none
      include "common"
      integer i,j,k
      real rho,ptot(0:3),eg,betacm(3),pcn(3),en,srts,sigr,mass,
     &     rhopot
      do i=1,10
         rho=(float(i)-0.5)/10.*rho0
         do j=1,200
            eg=0.8+float(j)/100.
            ptot(0)=eg+rmass
            ptot(1)=0.
            ptot(2)=0.
            ptot(3)=eg
            do k=1,3
               betacm(k)=ptot(k)/ptot(0)
               pcn(k)=0.
            end do
            en=rmass
            call lorentz(betacm(1),betacm(2),betacm(3),pcn(1),
     &           pcn(2),pcn(3),en)
            srts=sqrt(ptot(0)**2-ptot(1)**2-ptot(2)**2-ptot(3)**2)
*            write(*,*)'vor medrho',srts,rho,pcn(1),pcn(2),pcn(3),
*     &           betacm
            call medrho(srts,sigr,betacm,pcn(1),pcn(2),pcn(3),rho,
     &           mass,rhopot,103)
*            write(700,*)rho,eg,sigr
         end do
*         write(700,*)
      end do
      return
      end
            
      subroutine medrho(srts,sigr,betacm,pcnx,pcny,pcnz,rho,mass,
     &     rhopot,ires)
      implicit none 
      include "common"
      include "cominput"
      include "comwidth"
      real srts,sigr,matrix,pinitial,pfinal2,rn
      integer i,j,ires,k

      real mass1,minmass,maxmass,psa,dma,mass,            
     &     gamtot,ratio(numchr),spectral,pfinal,pcnx,pcny,pcnz,
     &     xout,pcm,betacm(3)
      real intfac,y,ymax,ymin,integr,pst,pst2,x
      real mres0,gamres0,bwmes
      real ratios(numchs),ratiom(nmesch2+nmesch3)
      integer nm2,nang,nmm
      real dy,dya
      parameter (nmm=200,dy=2.*pi/nmm,nang=20)
      real pscatt(3),pscattc(4),costk,phik,phi,cost,pout(3),erho,
     &     ps(nmm),rhopot
      logical flag
      
      real rho

      rhores=rho

*      write(*,*)'in medrho',srts,rho,pcnx,pcny,pcnz
      pinitial=(srts**2-rmass**2)/2./srts
      pcm=sqrt(pcnx**2+pcny**2+pcnz**2)
      mass1=rmass

*      ires=3+idbmax            
      if(ipotvec.ge.1) then
         call vecspot(rho,0.,rhopot,ires-idbmax)
      else
         rhopot=0.
      end if
*integrate over mass distribution
      if(ires.le.nres+1) then
         minmass=barprop1(1,1)+mesprop1(1,1)
         mres0=barprop1(ires,1)
         gamres0=barprop1(ires,2)
      else if(ires.le.nres+nsres+1) then
         minmass=sbarprop1(1,1)+mesprop1(1,1)         
         mres0=sbarprop1(ires-nres-1,1)
         gamres0=sbarprop1(ires-nres-1,2)
      else if(ires.gt.idbmax) then
         minmass=2.*mesprop1(1,1)-rhopot         
         mres0=mesprop1(ires-idbmax,1)
         gamres0=mesprop1(ires-idbmax,2)
         if(ires.eq.105.and.iwsc2.ge.2) then
            gamres0=gamres0+0.08*rho/rho0
         end if
      else
         write(*,*)'problems in massint2: wrong ires'
         stop
      end if
      maxmass=srts-mass1-rhopot
      if(maxmass.le.minmass) then
         pst=0.
*         ps(1)=0.
*         ps(2)=0.
*         ps(3)=0.
*         maxcontr=0.
*      else if(gamres0.lt.1e-05) then
*         if(srts.gt.mres0+mass1) then
*            ps(1)=sqrt((srts**2-(mres0+mass1)**2)*(srts**2-
*     &           (mres0-mass1)**2)/4./srts**2)
*         else
*            ps(1)=0
*         end if
*         ps(2)=0.
*         ps(3)=0.
      else
         pst=0.
*         write(*,*)maxmass,minmass,mres0,gamres0
         ymax=2.*atan((maxmass-mres0)
     &        /gamres0*2.)
*         write(*,*)ymax
         ymin=2.*atan((minmass-mres0)
     &        /gamres0*2.)
*         write(*,*)ymin
         nm2=max(int((ymax-ymin)/dy),1)
*         write(*,*)nm2
         dya=(ymax-ymin)/float(nm2)
*         write(*,*)dya,i,nm2
         do i=1,nm2
*            write(*,*)i
            y=ymin+(float(i)-0.5)*dya
*            write(*,*)i,y
            mass=.5*tan(y/2.)*gamres0+mres0
            mass=min(max(mass,minmass),maxmass)
            
*            write(*,*)srts,mass,mass1
            pfinal=sqrt((srts**2-(mass+rhopot+mass1)**2)*
     &           (srts**2-(mass+rhopot-
     &           mass1)**2)/4./srts**2)
            if(ires.eq.103) then
               matrix=0.16
            else if(ires.eq.105) then
               matrix=0.08*pfinal**2/(2.*(srts-1.73)**2+pfinal**2)
            else
               write(*,*)'wrong ires in medrho',ires
               stop
            end if
*            write(*,*)'pfinal=',pfinal
            ps(i)=0.
            do j=1,nang
*               write(*,*)mass,j
               call vecmesa(srts,3,mass+rhopot,cost)
               phi=2.*pi*rn(iseed)
*determine rho-meson momentum in lab
               if((pcny.ne.0.).or.(pcnx.ne.0.)) then 
                  phik=atan2(pcny,pcnx)
               else
                  phik=0.
               end if
*               write(*,*)pcm
               if(pcm.gt.0.) then
                  costk=pcnz/pcm
               else
                  costk=1.
               end if
*     cm frame:
               pscattc(1)=cos(phi)*sqrt(max(1.-cost**2,0.))
               pscattc(2)=sin(phi)*sqrt(max(1.-cost**2,0.))
               pscattc(3)=cost
*     rotation to lab:
*               write(*,*)cost,pscattc
               call rotate(costk,phik,pscattc,pscatt)
               xout=(srts**2-(mass+rhopot+rmass)**2)*(srts**2-
     &              (mass+rhopot-rmass)**2)/4./srts**2
               xout=sqrt(max(xout,0.))
*               write(*,*)'xout',xout,pscatt
               do k=1,3
                  pout(k)=pscatt(k)*xout
               end do
               erho=sqrt((mass+rhopot)**2+xout**2)
               call lorentz(-betacm(1),-betacm(2),-betacm(3),
     &              pout(1),pout(2),pout(3),erho)

               pres=sqrt(pout(1)**2+pout(2)**2+pout(3)**2)

               if(ires.le.nres+1) then
                  call manley(mass,ires,0,0,ratio,gamtot,0)
               else if(ires.le.nres+nsres+1) then
                  call swidth(mass,ires-nres-1,gamtot,ratios)
                  ratio(1)=ratios(1)
                  ratio(4)=0
               else 
*                  write(*,*)'vor bwmes',mass
                  if((iwsc2.ge.1.and.ires.eq.103).or.(iwsc2.ge.2.and.
     &                 ires.eq.105)) then
                     gamtot=bwmes((mass+rhopot)**2,
     &                    ires-idbmax,0,0,0,ratiom,1)
                  else 
                     gamtot=bwmes((mass+rhopot)**2
     &                    ,ires-idbmax,0,0,0,ratiom,0)
                  end if
*                  write(*,*)'nach bwmes',gamtot
                  ratio(1)=0
                  ratio(4)=0
               end if
               spectral=2./pi*mass**2*gamtot/((mass**2-
     &              mres0**2)**2+gamtot**2*mass**2)
            
               intfac=gamres0/((mass-mres0)**2
     &              +gamres0**2/4.)


               ps(i)=ps(i)+matrix*pfinal*spectral*dya/
     &              intfac/float(nang)
            end do
            pst=pst+ps(i)
         end do
      end if
      sigr=1000.*pst/pinitial/srts**2
*select mass
      x=rn(iseed)
      flag=.true.
      i=1
      pst2=0.
      do while(flag)
         pst2=pst2+ps(i)
         if(pst2.gt.x*pst) then
            flag=.false.
         else if(i.eq.nm2) then
            write(*,*)'problems in medrho',i,nm2,pst2,pst
            i=1
         else
            i=i+1
         end if
      end do
      y=ymin+(float(i)-0.5)*dya
      mass=.5*tan(y/2.)*gamres0+mres0
      mass=min(max(mass,minmass),maxmass)
      return
      end

