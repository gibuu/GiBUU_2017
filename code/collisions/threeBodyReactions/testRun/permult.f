      subroutine permult
      !subroutine for calculation of transition rates for N N Pi-> N N
      !Date : 14/10/03
      !Author :Oliver Buss
      !Note : This version is based on M.Effenbergers code.

      implicit none
      include "common"
      include "cominput"
      include "commcoit"
      include "comstat"

      !global variables:

      integer, dimension(1:2,0:1) :: istore !Array where Nuksearch saves nucleon id's into
      integer :: jj,init,irun  
      logical :: debug            ! Flag to switch on debug write statements
      logical :: Flag             
      real, dimension(1:3) :: g   ! Vector with decay width for different channels
                                  ! 1 : pp pion->x
                                  ! 2 : np pion->x
                                  ! 3 : nn pion->x
      real :: gtot                ! total decay width
      real :: rn                  ! random number function             
      integer :: i1,i2            ! nucleon id's, used by finalstate 
      real :: rhop,rhon           ! proton&neutron density
      logical,parameter ::   erzeugFlag=.true.

      debug=.false.
      If (debug) write(*,*)'in permult'
      !loop over all parallel runs*
      if(isstatea.eq.1)   return         !Switch for NNPi->NN

      do irun=1,num
         init=(irun-1)*maxper
         do jj=init+1,init+maxper
            flag=.false.
            if(idp(jj,1).ne.idbmax+1) cycle    !Only pions are relevant
            call NukSearch                     !Searches for nucleons in reach of pion
            If (flag) cycle                    !No Nucleons found
            call width                         !Evaluates decay width 
            if (flag) cycle                    !pion charge excludes possible channel
            call auswuerfeln                   !makes monte-carlo decision
            if (flag) cycle                    !no decay
            call finalstate                    !Generates Final State
         end do
      end do

      contains !internal procedures
      
***********************************************************************************
      subroutine NukSearch
            !Stores nearest nucleonic neighbours of pion into field istore
            implicit none

            integer,dimension (0:1) :: inucl
            real:: x,y,z,rho,radius,dist
            integer :: ivec
            integer :: j,i1,i,istart
            real :: derivd
            real, dimension(0:3) :: deriv

            if (debug) Print *, '*****************************'
            
            !Density at pion position
            x=rp(1,jj)
            y=rp(2,jj)
            z=rp(3,jj)
            ivec=4
            call dichte(x,y,z,deriv,derivd,ivec)
            rho=deriv(0)
            rhop=(deriv(0)+derivd)/2.
            rhon=(deriv(0)-derivd)/2.
            if(rho.lt.1e-03) then             !If density too small, then no absorption
               flag=.true.
               return
            end if

            !Radius to look for nucleons
            radius=Min(rmesdec*(rho0/rho)**(1./3.),5.0)
            If (debug) Print *,"Radius=", radius
            !Search for nucleons               
            inucl(0)=0          ! count protons found 
            inucl(1)=0          ! count neutrons found
            do i=0,1
               do j=1,2
                  istore(i,j)=0 !array i.o. to store id of found nuleons
               end do
            end do 
            !randomize starting point for loop over particles
            istart=int(rn(iseed)*masstot)
            do j=1,masstot
               i1=mod(j+istart,masstot)+1+(irun-1)*maxpe+maxbb
               if(id(i1,1).ne.1)   cycle
               if(idp(jj,3).eq.i1) cycle  
               dist=sqrt((x-r(1,i1))**2+(y-r(2,i1))**2+(z-r(3,i1))**2)
               if(dist.le.radius) then
                  if(((id(i1,2).eq.0).and.(inucl(0).le.1)).or.
     &                 ((id(i1,2).eq.1).and.(inucl(1).le.1))) then
                     inucl(id(i1,2))=inucl(id(i1,2))+1
                     istore(inucl(id(i1,2)),id(i1,2))=i1
                  end if
               end if
               if((inucl(1).eq.2).and.(inucl(0).eq.2)) exit !enough nucleons are found
            end do   
            if((inucl(1)+inucl(0)).lt.2) flag=.true. !two nucleons could not be found             
           If (debug) Print *, 'nach Nuksuche',inucl,flag  
           If (debug) then            
              Open(973,File="./oliver/NukSuche.dat")        
              if(inucl(1)+inucl(0).lt.2) then 
                 write(973,*)'inucl <2',inucl(1),inucl(0),radius,x,y,z
              else if((inucl(1)+inucl(0)).eq.2) then
                 write(973,*)'inucl =2',inucl(1),inucl(0),radius,x,y,z
              else if((inucl(1)+inucl(0)).eq.3) then 
                 !three nucleons have been found             
                 write(973,*)'inucl =3',inucl(1),inucl(0),radius,x,y,z
              else if((inucl(1)+inucl(0)).eq.4) then
                 write(973,*)'inucl =4',inucl(1),inucl(0),radius,x,y,z
              end if
            end if
      end subroutine NukSearch


************************************************************************************
      subroutine width
            !stores decay width into field g
            implicit none
            integer :: qpion
            real :: e1,e2,e3
            real :: momPion,gammaPion
            logical,dimension (3) :: kanalFlag
            integer :: counter
            real, dimension(3) :: srtfreek
            integer :: qtot
            real :: rho2,pcmf
            integer :: imat, index
            integer :: k,i
            
            qPion=idp(jj,2)

            Do i=1,3
               g(i)=0
               kanalFlag(i)=.false.
            end do
            gtot=0
            If (debug) Print *, gtot
            !Look for possible reaction channels
            If ((istore(1,1).ne.0).and.(istore(2,1).ne.0)
     &            .and.(qPion.ne.1)) then !pp-Kanal
               kanalflag(1)=.true.
            end if
            If ((istore(1,1).ne.0).and.(istore(1,0).ne.0)) then !pn-Kanal
               kanalflag(2)=.true.
            end if
            If ((istore(1,0).ne.0).and.(istore(2,0).ne.0)
     &            .and.(qpion.ne.-1)) then ! nn-Kanal
               kanalflag(3)=.true.
            end if
            If (debug) Print *, kanalflag

            !Count channels which are open
            counter=0
            Do k=1,3
               if(kanalflag(k)) counter=counter+1
            end do
            if (debug) Print  *, counter
            if (counter.eq.0) then
             If(debug) then
                Print *, 'Keine Nukleonen in subroutine width'
     &                 , istore(1,0), istore(2,0),istore(1,1)
     &                 , istore(2,1),qpion
                Print *, 'bzw. Keine Ladungserhaltung mit gefundenen'
     &                ,'Nukleonen moeglich'
             end if
             flag=.true.
             return
            end if
            If (Debug) Print *, pp(1,jj),pp(2,jj),pp(3,jj)
            momPion=SQRT(pp(1,jj)**2+pp(2,jj)**2+pp(3,jj)**2)
            If(Debug) Print *, mompion,SQRT(mompion**2+0.140**2)
            If (iOsetPionGamma.eq.1.and.momPion.lt.0.140) then !Absorption nach NPA 554 by Oset et al
               call OsetPionImaginaryMain(jj,gammaPion) 
               Do k=1,3  !Grobe Aufteilung in Kanäle
                 if (kanalflag(k)) g(k)=gammapion/float(counter)
               end do
               gtot=gammaPion
               If (debug) Print *, 'osetdecay',gtot,g
            else    !Detailed Balance
                do k=1,3        ! loop over all three possible channels 
                  imat=0
                  if((k.eq.1).and.kanalFlag(1)) then  !pp   
                     i1=istore(1,1)
                     i2=istore(2,1)
                     !rho2=ratp**2*rho**2
                     rho2=rhop**2/2.
                     if(qpion.eq.0) imat=1
                     if(qpion.eq.-1) imat=2
                  else if((k.eq.2).and.kanalFlag(2)) then !pn
                     i1=istore(1,1)
                     i2=istore(1,0)
                     !rho2=ratn*ratp*rho**2
                     rho2=rhop*rhon
                     if(abs(qpion).eq.1) imat=4
                     if(qpion.eq.0) imat=3
                  else if((k.eq.3).and.kanalFlag(3)) then !nn   
                     i1=istore(1,0)
                     i2=istore(2,0)
                     !rho2=ratn**2*rho**2
                     rho2=rhon**2/2.
                     if(qpion.eq.0) imat=1
                     if(qpion.eq.1) imat=2
                  end if
                  if(kanalflag(k)) then 
                     qtot=id(i1,2)+id(i2,2)+qpion
                     if((qtot.lt.3).and.(qtot.gt.-1)) then 
                        e1=sqrt(e(i1)**2+p(1,i1)**2+p(2,i1)**2+
     &                       p(3,i1)**2)
                        e2=sqrt(e(i2)**2+p(1,i2)**2+p(2,i2)**2+
     &                       p(3,i2)**2)
                        If (debug) Print *,'pion', ep(jj),mompion
                        e3=sqrt(ep(jj)**2+momPion**2)
                        If (debug) Print *,'energien', e1,e2,e3
                        srtfreek(k)=sqrt((e3+e1+e2)**2-
     &                       (pp(1,jj)+p(1,i1)+p(1,i2))**2-
     &                       (pp(2,jj)+p(2,i1)+p(2,i2))**2-
     &                       (pp(3,jj)+p(3,i1)+p(3,i2))**2)
                        If(debug) Print *, 'srtfreek',srtfreek(k)
                        if(srtfreek(k).lt.(2.*rmass+pmass)) then
                           Write(*,*)'Problems with srtfree'
     &                                ,srtfreek(k)
                           stop
                        end if
                        index=min(max(nint((srtfreek(k)-
     &                       sigs0)/delsigs),0),nsmax)
                        if (debug)Print *, 'index=',index
                        if((imat.le.0).or.(imat.ge.5)) then
                           write(*,*)'problems in mesdec imat',imat
                           stop
                        end if
                        pcmf=sqrt(srtfreek(k)**2/4.-rmass**2)
                        IF(debug) print *, pcmf
                        g(k)=matrixs(imat,index)*pcmf/
     &                       srtfreek(k)/4.
     &                       /pi*rho2/8./e1/e2/e3*(0.197)**6
                        !factor 2.57672 because matrixs contains still the mb of cross section
                        g(k)= g(k)*2.57672 
                        gtot=g(k)+gtot 
                        If (debug) Print *, 'detailed Balance',gtot,g
                     end if
                  end if
                end do
            end if 
            If (debug) Print *, 'nach width',gtot
      end subroutine width 
*******************************************************************************************
      subroutine auswuerfeln
              !chooses nucleons to make absorption with
              implicit none

              real :: w,rnxx 
              integer :: kanal

              !Monte-Carlo decision if particle is going to be absorbed
              w=exp(-dt*gtot/0.197)
              if (debug) Print *,w
              !dt anstatt dt0, da Zerfallsbreite im Laborsystem berechnet wurde
              if (Debug) Print *, "in auswuerfeln",gtot,g
              rnxx=rn(iseed)

              if(rnxx.le.w) then !No absorption
                 If (Debug) Print *, rnxx,'<',w,'No absorption'
                 flag=.true.
                 return
              end if
              If (Debug) Print *, rnxx,'>',w,'Absorption'
              

              if(rnxx.le.(g(1)/gtot)) then
                 kanal=1
              else if(rnxx.le.(g(1)+g(2))/gtot) then
                 kanal=2
              else
                 kanal=3
              end if
              if(kanal.eq.1) then        !pp   
                     i1=istore(1,1)
                     i2=istore(2,1)
              else if(kanal.eq.2) then   !pn
                     i1=istore(1,1)
                     i2=istore(1,0)
              else if(kanal.eq.3) then   !nn   
                     i1=istore(1,0)
                     i2=istore(2,0)
              end if
              if((i1.eq.0).or.(i2.eq.0)) then
                   write(*,*)'probleme in wuerfeln',g(1),g(2),g(3),gtot
                   stop
              end if
              if(i1.eq.i2) then
                  write(*,*)'probleme in wuerfeln,i1=i2',i1,i2
                  stop
              end if
              
              if (debug) Print *, kanal,g,gtot

      end subroutine auswuerfeln
 
*************************************************************
******************************************************************
      subroutine finalstate
          !determine the kinematics of outgoing nucleons, which stem from nucleons
          !choosen in auswuerfeln*
          implicit none

          logical ::paulflag,flagiter,forflag,leadp
          real,dimension(0:3) :: p1, ptot
          real ::pot1, mass1, enuk1,enuk2,srtfree,epion
          real :: etotal,phase
          integer :: k,j,hole,j3
          integer ::niter, qtot,ntag
          real :: cost, sint,phi,xout,rnxx,eout
          real,dimension (1:3,1:2) :: pout
          real,dimension (1:2) ::massout,potout
          integer,dimension(1:2) :: idout
          integer :: numparf,collid2,collid3
          real :: x,y,z
          integer, dimension(2) :: izout,nucProd
          real :: weightPion !perweight of pion
              do k=1,3
                 p1(k)=pp(k,jj)
              end do
              mass1=ep(jj)
              pot1=upotp(jj)
              p1(0)=sqrt(ep(jj)**2+pp(1,jj)**2+pp(2,jj)**2+
     &        pp(3,jj)**2)
              weightPion=perweight(jj)

              enuk1=sqrt(e(i1)**2+p(1,i1)**2+p(2,i1)**2+
     &                   p(3,i1)**2)
              enuk2=sqrt(e(i2)**2+p(1,i2)**2+p(2,i2)**2+
     &                       p(3,i2)**2)
         
              srtfree=sqrt((p1(0)+enuk1+enuk2)**2-
     &                       (pp(1,jj)+p(1,i1)+p(1,i2))**2-
     &                       (pp(2,jj)+p(2,i1)+p(2,i2))**2-
     &                       (pp(3,jj)+p(3,i1)+p(3,i2))**2)

              do k=1,3
                 betlrf(1,k)=betlrfboo(i1,k)
                 betlrf(2,k)=betlrfboo(i2,k)
                 r1(k)=r(k,i1)
                 r2(k)=r(k,i2)
              end do

              j01=betlrfboo(i1,4)
              j02=betlrfboo(i2,4)
              j01n=betlrfboo(i1,5)
              j02n=betlrfboo(i2,5)

              if(ipotcrosw.eq.0) then
                 epion=sqrt(ep(jj)**2+pp(1,jj)**2+pp(2,jj)**2+
     &                      pp(3,jj)**2)
                 enuk1=sqrt(e(i1)**2+p(1,i1)**2+p(2,i1)**2+
     &                   p(3,i1)**2)
                 enuk2=sqrt(e(i2)**2+p(1,i2)**2+p(2,i2)**2+
     &                   p(3,i2)**2)
              else
                 epion=sqrt((ep(jj)+upotp(jj))**2+pp(1,jj)**2
     &                     +pp(2,jj)**2+pp(3,jj)**2)
                 enuk1=sqrt((e(i1)+upot(i1))**2+p(1,i1)**2
     &                      +p(2,i1)**2+p(3,i1)**2)
                 enuk2=sqrt((e(i2)+upot(i2))**2+p(1,i2)**2
     &                      +p(2,i2)**2+p(3,i2)**2)
              end if
              etotal=enuk1+enuk2+epion
              do k=1,3
                  ptot(k)=p1(k)+p(k,i1)+p(k,i2)
                  betacm(k)=ptot(k)/etotal
              end do
              srts=sqrt(etotal**2-ptot(1)**2-ptot(2)**2-
     &                 ptot(3)**2) 
              If (debug) Print*, "srts vorher=",srts
              em3=rmass
              em4=rmass
              id3=1
              id4=1
              qtot=id(i1,2)+id(i2,2)+idp(jj,2)
              

*     charge assignmnet
              If (Debug) Print *,'total charge', qtot
              If (Debug) Print *,'nucleonic charges', id(i1,2),id(i2,2)
              If (Debug) Print *,'pion charge', idp(jj,2)
              
              if(qtot.eq.1) then
                  rnxx=rn(iseed)
                  if(rnxx.lt.0.5) then
                     iz3=1
                     iz4=0
                  else
                     iz3=0
                     iz4=1
                  end if
               else
                  iz3=nint(qtot/2.)
                  iz4=nint(qtot/2.)
               end if
               If (Debug) Print *,'iz3,iz4',iz3,iz4
               

                  
               flagiter=.false.
               niter=0
               do while(.not.flagiter)
                     niter=niter+1
                     cost=2.*(rn(iseed)-0.5)
                     sint=sqrt(max(1.-cost**2,0.))
                     phi=2.*pi*rn(iseed)
                     pscatt(1)=sint*cos(phi)
                     pscatt(2)=sint*sin(phi)
                     pscatt(3)=cost
                     if(ipotcrosw.ge.1) then
                        call siter(flagiter,xout)
                     else
                        xout=sqrt((srtfree**2-4.*rmass**2)/4.)
                        enuk3=sqrt(rmass**2+xout**2)
                        enuk4=sqrt(rmass**2+xout**2)
                        flagiter=.true.
                     end if
                     if(niter.gt.nitermax) then
                        write(*,*)'problems in kinematik: niter.gt.
     &                       nitermax',srts,srtfree,cost,sint,phi
                        stop
                     end if
                end do

                !lorentz boost to lab frame
                do k=1,3
                     pout(k,1)=pscatt(k)*xout
                     pout(k,2)=-pout(k,1)
                end do
                call lorentz(-betacm(1),-betacm(2),-betacm(3),
     &                 pout(1,1),pout(2,1),pout(3,1),enuk3)

                call lorentz(-betacm(1),-betacm(2),-betacm(3),
     &                 pout(1,2),pout(2,2),pout(3,2),enuk4)
                massout(1)=em3
                massout(2)=em4
                potout(1)=u3
                potout(2)=u4
                idout(1)=1
                idout(2)=1
                izout(1)=iz3
                izout(2)=iz4
                numparf=2



               !check for pauli-blocking
               if(ipauli.ge.1) then
                  paulflag=.true.
                  j=1
                  ntag=0
                  do while(paulflag)                     
                     if(idout(j).eq.1) then
                        eout=sqrt((massout(j)+potout(j))**2+
     &                       pout(1,j)**2+pout(2,j)**2+
     &                       pout(3,j)**2)
                        x=rp(1,jj)
                        y=rp(2,jj)
                        z=rp(3,jj)
                        call pauli(ntag,phase,x,y,z,
     &                       pout(1,j),pout(2,j),pout(3,j),eout) 
                     end if
                     j=j+1
                     if(j.gt.numparf.or.ntag.eq.-1) paulflag=.false.
                  end do
                  !check for pauli blocking
                  if(ntag.eq.-1) then
                     If (Debug) Print *, "PauliBlocking"
                     return
                  end if
               end if

               !collision was allowed by pauli principle

               !settings
               call pionGridAbs(jj,2)
               idp(jj,1)=0
               

               if(erzeugFlag) then
                 do j=1,numparf
                  Do hole=init+1,init+maxper
                       if (idp(hole,1).eq.0) then
                          j3=hole
                          exit
                       end if
                       if (hole.eq.init+maxper) then
                          Print *, "No Hole in particle vector"
                          stop
                       end if
                  end do
                  If (Debug) Print *,jj,'Loch',j3,init+maxper 
                  forflag=.false.
                  leadp=.true.
                  collid2=i1
                  collid3=i2
                  call setidp(j3,pout(1,j),pout(2,j),
     &                          pout(3,j),
     &                          massout(j),potout(j),idout(j),
     &                          izout(j),
     &                          collid2,collid3,leadp,forflag,x,
     &                          y,z)
                  perweight(j3)=weightPion
                  If (debug) then
                      Print *, "Final states:",idp(j3,1),idp(j3,2)
                      Print *, "Viererimpuls", pp(0,j3),pp(1,j3)
     &                                         ,pp(2,j3),pp(3,j3)
                      Print *, "Ort", rp(1,j3),rp(2,j3),rp(3,j3)
                      Print *, "Breiten ", g
                      Print *, "Ursprungsteilchen",id(i1,2),id(i2,2)
                  end if
                  if(debug) then !Check energy conservation
                    if(j.eq.1) then
                       enuk1=sqrt((ep(j3)+upotp(j3))**2+pp(1,j3)**2
     &                     +pp(2,j3)**2+pp(3,j3)**2)
                       do k=1,3
                         ptot(k)=pp(k,j3)
                       end do
                       Print *, 'enuk1,ptot',enuk1,pp(1,j3),pp(2,j3)
     &                                      ,pp(3,j3)
                    else
                     enuk2=sqrt((ep(j3)+upotp(j3))**2+pp(1,j3)**2
     &                     +pp(2,j3)**2+pp(3,j3)**2)
                     do k=1,3
                         ptot(k)=ptot(k)+pp(k,j3)
                     end do
                     Print *, 'enuk2,ptot',enuk2,pp(1,j3),pp(2,j3)
     &                                      ,pp(3,j3)
                    end if
                  end if
                 end do
               end if
               If (debug) then
                 etotal=enuk1+enuk2
                 Print *, etotal, ptot(1),ptot(2),ptot(3)
                 srts=sqrt(etotal**2-ptot(1)**2-ptot(2)**2-
     &                 ptot(3)**2) 

                 Print *, "srts nachher=",srts
                stop
               end if
               
               ldecp(1+idbmax,-1,0)=ldecp(1+idbmax,-1,0)+
     &                              perweight(jj)
       end subroutine finalstate
************************************************************************

       end Subroutine permult
