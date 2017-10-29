      subroutine gamabsper                 
      implicit none
      include "common"
      include "cominput"
      include "commcoit"
      include "comwidth" 
      real rn
      integer nres2,ntheta,k
      parameter (nres2=4,ntheta=50)
      integer numtryt,ifra,numtry,maxnum,numpro(0:13),numcrt
      integer ibeg,i2,ka,k5,i5,k3
      real egamma,pz1,sigabs
      integer irun,i,qnuk,kmin,kmax,j,k2,i3
      real px2,py2,pz2,etn,pdel,pci2,pcf2,qpion,pzt,et
      real gamma,traf,pnbeta,transf,pcnz,pcny,
     &     pcnx,pcm
      real sigmapi,x,dsigdo(4,0:nres2),dsigdot(4),dsdot2(4),dsdot3,
     &     probka(4,-ntheta:ntheta-1),prob(-ntheta:ntheta-1,4,0:nres2)
     &     ,dsigma(-ntheta:ntheta-1)
      real sig2(nres2),sig2pi(0:3),sigx,sigx2,sigtot2(maxens),sigtot,
     &     sigtoto,meffr,potr,ener

      real x5,x6,x3,x4,sum,sum2,cost,phi,sint,phik,costk,
     &     encm,p1l,p2l,p3l
      real enuk,epion2,epion1,p3(3,3),eout2(3),ptotal,
     &     phase,x7
      real pr00,pr10
      real em5,u5
      integer ib,ib1,ib2
      integer id5,iz5
      integer ntag,inu,inuk,isu,isu2,idm,idch
      integer num2(nmes,-1:1),angbins,maxpt,nth,npt,maxtk,ntk
      real dpt,dtk

      parameter (angbins=20,maxpt=150,dpt=0.05,dtk=0.02,maxtk=200)
      real dist2,theta

      real tsigma(nmes,-1:1),asigma(nmes,-1:1,0:angbins-1),
     &     msigma(nmes,-1:1,0:maxpt),ptrans,
     &     ptot2,tkin,tksigma(nmes,-1:1,0:maxtk),
     &     a2sigma(nmes,-1:1,0:angbins-1)
      real dprodt(0:maxtk)
      integer numcre(10,nmes,-1:1),maxncre,inumcre,pertyp2
*two pion
      integer anzp,maxanzp
      parameter (maxanzp=1500)
      integer ipis(maxanzp,6),imass
      real ppis(3,maxanzp),eps(maxanzp),ener1,ener2,invmass
      integer gidps(maxanzp)
      real pervs(maxanzp)

      integer iii
      integer nmass
      real dmass
      parameter (dmass=0.01,nmass=100)
      real t2sigma(-3:2),t2sigdm(-3:2,0:nmass),t2sigdmr(-3:2,0:nmass)


      integer maxanz,ioutput
      parameter (maxanz=40,ioutput=0)
      integer id1s(maxanz),id2s(maxanz),id5s(maxanz),id6s(maxanz),
     &     anz
      real pnukl(3,maxanz),rnukl(3,maxanz),potnuk(maxanz)
      real finmass(2),scapot(2),meff,velon(0:3),velog(3),vrel,
     &     pinitial,vr(3),spot,resmass
      logical flagps,flagiter
      integer niter
      real xout,phist(-ntheta:ntheta-1),
     &     pscattc(3),psratio,etotal,ptot(0:3)
      
      logical debugflag

*fritiof      
      integer maxb,qp,id1,id2
      real pcm2(0:3),beta(3),em1,em2,prcm,sigf
      real xpom,epspom,yreg,etareg
      parameter (xpom=0.071,yreg=0.12,etareg=0.46,epspom=0.075)
*shadowing
      real r22,costr,shadfac,boffif,sigs
      integer ishadow2
*
      logical frflag,vecflag,forflag,leflag,paulflag,iflag,jflag,
     &     leadp,eflag
      real sigvec(3),sigveci(3),sigvect,sigveex,sigvein
      real pout(3,nprom),massout(nprom),potout(nprom)
      integer idout(nprom),izout(nprom),numparf
      integer sumid,j3,i1
      real massmi,eout,massr,masso,poto,pr,massvac

      integer ipert,ipert2
      logical perflag

      real pfinal2,pfinalo,gamtot,bwmes,ratiom(nmesch2+nmesch3),
     &     spectralo,spectraln,pso,psn
*strangeness
      real sigst(3),sigstt
      

      save numtryt, ifra, numtry

*PM+++++++++++++++++++++++++++++++++++++++++++++++++++!
      integer::rCh                                    !
      real,dimension(0:3)::pi1,pi2,po1,po2,po3        !
*PM+++++++++++++++++++++++++++++++++++++++++++++++++++!


!      logical, parameter :: pascalTwoPi=.true.
      logical, parameter :: OnlyTwoPi=.true.
      real, dimension(1:3) ::ort
      integer, dimension(1:3) ::ladung
      logical,parameter :: oliverDebug=.false.

      interface
         subroutine DreierEnergie(p3,ptot,beta,ladung,Ort)
         implicit none
         real, dimension(0:3), intent(in) :: ptot
         real, dimension(1:3,1:3), intent(inout) :: p3 
         real, dimension(1:3), intent(in) :: Ort,beta
         integer, dimension(1:3),intent(in) :: ladung
         end subroutine DreierEnergie
      end interface
      interface
         subroutine EventGenerator(r,pi1,pi2,po1,po2,po3)
         implicit none
         include "cominput"
         real,dimension(0:3)::pi1,pi2,po1,po2,po3
         integer::r
         end subroutine EventGenerator
      end interface

      write(*,*)'*************************************'
      write(*,*)'In GamAbsPer! Elab=',elab
      write(*,*)'*************************************'
      write(*,*)

      Print *, "energieTwoPi =",energieTwoPi

      write(*,*)'vor loop ueber ensembles'
      egamma=elab
      pz1=egamma                                      
                                                                       
      maxnum = 1                                                       
      
      numtry=0

      do i=0,13
        numpro(i)=0
      end do

      write(*,*)'hier'
*      debugflag=.true.
      debugflag=.false.
* --- loop over parallel runs:                                         

      sigabs=0.

      ishadow2=ishadow
      if(egamma.lt.0.8.and.ishadow.eq.1) then
         write(*,*)'egamma lt 1. => ishadow2=0'
         ishadow2=0
      end if
      if(ishadow2.eq.1) then
         call inidichte2(massta,mstapr)
         call inipath
      end if

      do 1000 irun = 1, num
         if(debugflag) write(*,*)irun
         if(oliverDebug) write(*,*)'In Ensemble:', irun
         numcrt = 0
         ipert=0
*     ----- loop over testparticles in one run:                            
         
         ibeg=(irun-1)*maxpe+maxbb                                       
         do 100 i=ibeg+1,ibeg+masstot                               
            if(oliverDebug) write(*,*) 'Teilchen:',i
            pertyp2=0
            if(debugflag) write(*,*)i
            if(ipert.gt.maxper) then
               write(*,*)'perturbative vector too small in gamabs',
     &              ipert,irun,i
               stop
            end if
                                                               
            numtry=numtry+1
            px2=p(1,i)                                                 
            py2=p(2,i)                                                 
            pz2=p(3,i)
            qnuk=id(i,2)
            kmin=-2*qnuk+3
            kmax=kmin+1
            if(ipotcrosw.eq.0) then
               meff=rmass
            else
               meff=rmass+upot(i)
            end if
            etn=sqrt(meff**2+px2**2+py2**2+pz2**2)
            pdel=sqrt((pz2+pz1)**2+px2**2+py2**2)
            srts=sqrt((egamma+etn)**2-pdel**2)
            ptot(0)=egamma+etn  !Total energy in calculation frame
            ptot(1)=px2
            ptot(2)=py2
            ptot(3)=pz2+pz1

            j01=betlrfboo(i,4)
            j02=betlrfboo(i,4)
            j01n=betlrfboo(i,5)
            j02n=betlrfboo(i,5)

*     Lorentz-Trafo into cm-frame
            pzt=pz2+pz1
            et=egamma+etn
            betacm(3)=pzt/et
            betacm(2)=py2/et
            betacm(1)=px2/et
            pcnx=px2
            pcny=py2
            pcnz=pz2
            encm=etn
            call lorentz(betacm(1),betacm(2),betacm(3),
     &           pcnx,pcny,pcnz,encm)
            pcm=sqrt(pcnx**2+pcny**2+pcnz**2)
            
            pci2=(srts**2-meff**2)**2/(4.*srts**2)
            if(abs(sqrt(pci2)-pcm).gt.1e-05) then
               write(*,*)'problems in gamabs',sqrt(pci2),pcm
            end if
            srtfree=sqrt(rmass**2+pci2)+sqrt(pci2)  !=E(Nukleon)+E(gamma)
                                                    !Gesamtimpuls=0 in CM-Frame

            if(srtfree.le.rmass+pmass) goto 100

            if(iphoch.ne.3) then
*     determine cross sections for vector meson production
               call vecmes(srtfree,sigvec,sigveci,1,j01)
*
               if(iwsc2.ge.1) then
                  call medrho(srtfree,sigvec(1),betacm,pcnx,pcny,pcnz,
     &                 j01,massr,potr,103)
                  if(iwsc2.ge.2) then
                     call medrho(srtfree,sigvec(2),betacm,pcnx,pcny,
     &                    pcnz,j01,masso,poto,105)
                  end if
               end if

               sigvect=sigvec(1)+sigvec(2)+sigvec(3)
               sigveex=sigvect
**achtung
*               do k2=1,3
*                  sigveci(k2)=0.
*               end do
**               

*add inclusive vector meson production
               sigvein=0.
               if(srtfree.lt.thresgn) then
                  sigvein=sigveci(1)+sigveci(2)+sigveci(3)
                  sigvect=sigvect+sigvein
               end if
            end if
*determine strangeness cross sections
            call photostrange(srtfree,sigst)
            sigstt=sigst(1)+sigst(2)+sigst(3)
*            write(600,*)srtfree,sigveci(1),sigvec(1),sigvec(2),
*     &           sigveci(2)

            forflag=.false.
            leflag=.false.
            sigtot=0
            sigtoto=0
*     fritiof
            if(srtfree.gt.thresgn.and.ifritzi2.eq.1.and.
     &           iphoch.ne.2) then
               qp=0
*pcm2 is center of mass momentum of rho meson (first particle in fritiof)
               pcm2(1)=-pcnx
               pcm2(2)=-pcny
               pcm2(3)=-pcnz
*id1 for fritiof (initialize photon)
               id1=idbmax+nmes+1
               id2=1
               em1=0.
               em2=e(i)
               frflag=.true.
*determine total cross section and subtract contributions that are
*treated independent (vector mesons, strangeness)
               sigtoto=(xpom*srtfree**(2.*epspom)+
     &              yreg*srtfree**(-2.*etareg))*1000.-sigvect-sigstt
               do while(frflag)
                  frflag=.false.
                  call fritzi(pcm2,srtfree,em1,em2,id1,id2,qp,qnuk,
     &                 sigf,
     &                 prcm,betacm,idout,izout,massout,pout,potout,
     &                 numparf,.false.,.true.)
*     dont allow exclusive vector meson production because this is already
*     taken into account
                  if(numparf.eq.2.and.min(idout(1),idout(2)).eq.1) 
     &                 then
                     sumid=idout(1)+idout(2)-idbmax
                     if(sumid.eq.4.or.sumid.eq.6.or.sumid.eq.8) 
     &                    then
                        frflag=.true.
                     end if                     
                  end if

*     dont allow K Kbar N production via fritiof
                  if(numparf.eq.3.and.
     &                 min(idout(1),idout(2),idout(3)).eq.1) then
                     sumid=idout(1)+idout(2)+idout(3)-2*idbmax
                     if(sumid.eq.22.and.(idout(1).eq.110.or.
     &                    idout(2).eq.110.or.idout(3).eq.110)) then
                        frflag=.true.
                     end if
                  end if

                  pertyp2=13
                  do j=1,numparf
                     if(idout(j).eq.idbmax+3) then
                        pertyp2=4
                     else if(idout(j).eq.idbmax+5) then
                        pertyp2=5
                     else if(idout(j).eq.idbmax+7) then
                        pertyp2=6
                     end if
                  end do

                  forflag=.true.
               end do
            else if(iphoch.ne.2) then
               pinitial=sqrt(pci2)
               pcf2=((srtfree**2-rmass**2-pmass**2)**2-
     &              4.*rmass**2*pmass**2)/(4.*srtfree**2)          
               qpion=sqrt(max(pcf2,0.))
*     Integration ueber Streuwinkel

*     initialization of 1 pion routine
               if(debugflag) write(*,*)'vor inidsdo',qnuk,srts 
               call inidsdo(qnuk,srtfree)
               if(debugflag) write(*,*)'nach inidsdo' 
               
               sigmapi=0.

               if(ipotcrosw.ge.1) then
*     set variables needed in siter
                  id3=1
                  id4=1+idbmax
*avoid charge sum here
                  iz3=qnuk
*
                  do j=1,2
                     do k2=1,3
                        betlrf(j,k2)=betlrfboo(i,k2)
                     end do
                  end do
                  j02=j01
                  do k2=1,3
                     r1(k2)=r(k2,i)
                     r2(k2)=r(k2,i)
                  end do
                  em3=rmass
                  em4=pmass
*     determine masses of produced resonances
                  call detmass(ptot,resmass,spot,1,qnuk)
                  finmass(1)=resmass
                  scapot(1)=spot
                  call detmass(ptot,resmass,spot,2,qnuk)
                  finmass(2)=resmass
                  scapot(2)=spot
               else
                  do k2=1,2
                     finmass(k2)=srtfree
                     scapot(k2)=0.
                  end do
               end if
               if(irrespar.eq.3) then
                  do k2=1,2
                     finmass(k2)=srtfree
                     scapot(k2)=0
                  end do
               end if

               do i3=-ntheta,ntheta-1
                  x=(float(i3)+0.5)/float(ntheta)
*     monte-carlo for phi-integration sufficient
                  if(ipotcrosw.ge.4) then
                     phi=rn(iseed)*2.*pi
                     phist(i3)=phi

*     determine scattering direction in lab frame
                     if((pcny.ne.0.).or.(pcnx.ne.0.)) then 
                        phik=atan2(pcny,pcnx)
                     else
                        phik=0.
                     end if
                     costk=pcnz/pcm
*     cm frame:
                     pscattc(1)=cos(phi)*sqrt(1.-x**2)
                     pscattc(2)=sin(phi)*sqrt(1.-x**2)
                     pscattc(3)=x
*     rotation to lab:
                     call rotate(costk,phik,pscattc,pscatt)

                     call siter(flagiter,xout)
                  else
                     xout=qpion
                     flagiter=.true.
                  end if
                  
                  if(flagiter) then
                     pres=sqrt(ptot(1)**2+ptot(2)**2+ptot(3)**2)
                     rhores=j01
*     write(401,*)srtfree,pres,rhores
                     call dsdo(srtfree,x,dsigdo,dsigdot,finmass,
     &                    pinitial)
*     
                     if(finmass(1).le.rmass+pmass.or.
     &                    irrespar.eq.3) then
                        do k2=kmin,kmax
                           do i2=2,4
                              dsigdo(k2,0)=dsigdo(k2,0)+
     &                             dsigdo(k2,i2)
                              dsigdo(k2,i2)=0.
                           end do
                        end do
                     end if
                     if(finmass(2).le.rmass+pmass.or.
     &                    irrespar.eq.3) then
                        do k2=kmin,kmax
                           dsigdo(k2,0)=dsigdo(k2,0)+dsigdo(k2,1)
                           dsigdo(k2,1)=0.
                        end do
                     end if

                     dsdot3=0.
                     do k2=kmin,kmax
                        dsdot2(k2)=0.
                        do j=0,nres2
                           if(j.ne.3) then
*     exclude N(1535) because of eta production
                              dsdot2(k2)=dsdot2(k2)+dsigdo(k2,j)
                           end if
                        end do
                        dsdot3=dsdot3+dsigdot(k2)
                     end do
                     do k2=kmin,kmax
                        probka(k2,i3)=dsigdot(k2)/dsdot3
*     probability for N(1535)
                        prob(i3,k2,3)=dsigdo(k2,3)/
     &                       dsigdot(k2)
                        do i2=0,nres2
                           if(i2.ne.3) then
                              prob(i3,k2,i2)=dsigdo(k2,i2)/
     &                             dsdot2(k2)*(1.-prob(i3,k2,3))
                           end if
                        end do
                     end do
                     if(ipotcrosw.ge.4) then
*     phase space correction of dsdot3
                        call psfac(flagps,xout,psratio)                   
                        dsdot3=psratio*dsdot3
*     write(303,*)u3,u4,upot(i),upot(i)-u3,psratio
***************
                     end if
                     dsigma(i3)=dsdot3*2.*pi/float(ntheta)*389.
                     sigmapi=sigmapi+dsigma(i3)
                  else
                     dsigma(i3)=0.
                  end if
               end do
               if(debugflag) write(*,*)'vor 2pi'
*     
*     2 Pionen/eta  Anteil
               call gamma2pi(qnuk,srtfree,sig2,sig2pi,finmass,onlyTwoPi)
               if(oliverdebug) write(*,*)'nach 2pi-Querschnitt'
               sigx=sig2(1)+sig2(2)+sig2(3)+sig2(4)
               sigx2=sig2pi(0)
               sigtoto=sigmapi+sigx+sigx2
               leflag=.true.
            end if
            
            if(iphoch.eq.3) then
               sigvect=0
            end if

            sigtot=sigtoto+sigvect+sigstt

            
            if(sigtot.lt.1e-05) goto 100

            vecflag=.false.

            x7=rn(iseed)
            if(debugflag) write(*,*)x7,sigtot,sigvect
*Monte-Carlo starts
            
            ! Falls nur 2Pi events produziert werden sollen,
            ! dann nicht würfeln, sondern direkt in 2Pi
            ! Produktion einsteigen :
            if(OnlyTwoPi) Goto 238


            if(x7.le.sigveex/sigtot) then
*     vector meson production
               vecflag=.true.
               forflag=.false.
               leflag=.false.
               numparf=2
               idout(1)=1
               if(x7.le.sigvec(1)/sigtot) then
*     rho
                  idout(2)=3+idbmax
                  massmi=2*pmass
                  pertyp2=1
               else if(x7.le.(sigvec(1)+sigvec(2))/sigtot) then
*     omega
                  idout(2)=5+idbmax
                  massmi=3.*pmass
                  pertyp2=2
               else
*     phi
                  massmi=3.*pmass
                  idout(2)=7+idbmax
                  pertyp2=3
               end if
               izout(1)=qnuk
               izout(2)=0
               iz3=qnuk
*     determine masses
               massout(1)=rmass

               flagiter=.false.
               niter=0
               do while(.not.flagiter)
                  niter=niter+1
*mass
                  if(ipotvec.ge.1) then
                     call vecspot(j01,0.,spot,idout(2)-idbmax)
                  else
                     spot=0.
                  end if
                  call resmasdi(srtfree,massmi,rmass,-idout(2)+idbmax,
     &                 massout(2),spot)
                  meffr=massout(2)+spot
                  if(iwsc2.ge.1.and.idout(2).eq.3+idbmax) then
                     massout(2)=massr
                     meffr=massr+potr
                  else if(iwsc2.ge.2.and.idout(2).eq.5+idbmax) then
                     massout(2)=masso
                     meffr=masso+poto
                  end if
*
*     scattering angle
                  call vecmesa(srtfree,idout(2),meffr,cost)

                  phi=2.*pi*rn(iseed)

*     determine scattering direction in lab frame
                  if((pcny.ne.0.).or.(pcnx.ne.0.)) then 
                     phik=atan2(pcny,pcnx)
                  else
                     phik=0.
                  end if
                  costk=pcnz/pcm
*     cm frame:
                  pscattc(1)=cos(phi)*sqrt(max(1.-cost**2,0.))
                  pscattc(2)=sin(phi)*sqrt(max(1.-cost**2,0.))
                  pscattc(3)=cost
*     rotation to lab:
                  call rotate(costk,phik,pscattc,pscatt)

                  
                  em3=massout(1)
                  em4=massout(2)

                  if(ipotcrosw.ge.1) then
*     set variables for iteration routines
                     id3=idout(1)
                     id4=idout(2)
                     if(debugflag) write(*,*)'vor siter vec',id3,
     &                    id4,em3,em4,cost
                     call siter(flagiter,xout)
                     if(debugflag) write(*,*)'nach siter',flagiter
                     if(ipotvec.eq.2.and.(idout(2).eq.103.or.
     &                    idout(2).eq.105)) then
*determine first "old" pase space + spectral function factor
*(with potential at p=0)
                        if(.not.flagiter) then
*cross section=0
                           numparf=0
                           goto 100
                        end if
                        pfinal2=(srtfree**2-(rmass+meffr)**2)*
     &                       (srtfree**2-(rmass-
     &                       meffr)**2)/4./srtfree**2
                        if(pfinal2.le.0.) then
                           write(*,*)'problems in gamabsper ipotv',
     &                          srtfree,rmass,meffr,spot,em3,em4
                           stop
                        end if
                        pfinalo=sqrt(pfinal2)
*"old" spectral function
                        gamtot=bwmes(meffr**2
     &                       ,idout(2)-idbmax,0,0,0,ratiom,0)
                        spectralo=2./pi*em4**2*gamtot/((em4**2-
     &                       mesprop1(idout(2)-idbmax,1)**2)**2+
     &                       gamtot**2*em4**2)
*"new" factors
                        gamtot=bwmes((em4+u4)**2
     &                       ,idout(2)-idbmax,0,0,0,ratiom,0)
                        spectraln=2./pi*em4**2*gamtot/((em4**2-
     &                       mesprop1(idout(2)-idbmax,1)**2)**2+
     &                       gamtot**2*em4**2)
                        sigtot=sigtot*xout/pfinalo*spectraln/spectralo
                     end if
                        
*     if(ipotcrosw.ge.2.and.flagiter) then
*     call psfac(flagiter,xout,psratio)
*     end if

                  else
                     xout=(srtfree**2-(em3+em4)**2)*(srtfree**2-
     &                    (em3-em4)**2)/4./srtfree**2
                     if(xout.lt.0) then
                        write(*,*)'fehler in gamabsper xout.lt.0',
     &                       xout,srtfree,em3,em4,idout(1),idout(2),
     &                       id3,id4
                        stop
                     end if
                     xout=sqrt(xout)
                     enuk3=sqrt(em3**2+xout**2)
                     enuk4=sqrt(em4**2+xout**2)
                     flagiter=.true.
                  end if
                  
                  if(niter.gt.nitermax) then
                     write(*,*)'iteration not succesful vector'
                     goto 100
                  end if
               end do
*     lorentz boost
               do k=1,3
                  pout(k,1)=pscatt(k)*xout
                  pout(k,2)=-pout(k,1)
               end do
               if(debugflag) write(*,*)'vor lorentz',betacm,enuk3
               call lorentz(-betacm(1),-betacm(2),-betacm(3),
     &              pout(1,1),pout(2,1),pout(3,1),enuk3)
               
               call lorentz(-betacm(1),-betacm(2),-betacm(3),
     &              pout(1,2),pout(2,2),pout(3,2),enuk4)
               potout(1)=u3
               potout(2)=u4
*
            else if(x7.le.sigvect/sigtot) then
               vecflag=.true.
               forflag=.false.
               leflag=.false.
               numparf=3
               idout(1)=1
               idout(2)=1+idbmax
               if(x7.le.(sigveex+sigveci(1))/sigtot) then
*     rho
                  idout(3)=3+idbmax
                  massmi=2*pmass
                  pertyp2=1
               else if(x7.le.(sigveex+sigveci(1)+sigveci(2))/sigtot) 
     &                 then
*     omega
                  idout(3)=5+idbmax
                  massmi=3.*pmass
                  pertyp2=2
               else
*     phi
                  massmi=3.*pmass
                  idout(3)=7+idbmax
                  pertyp2=3
               end if
               if(rn(iseed).le.0.5) then
                  izout(1)=qnuk
                  izout(2)=0
               else
                  izout(1)=-qnuk+1
                  izout(2)=2*qnuk-1
               end if
               izout(3)=0
*     determine masses
               massout(1)=rmass
               massout(2)=pmass
*mass
               if(ipotvec.ge.1) then
                  call vecspot(j01,0.,spot,idout(3)-idbmax)
               else
                  spot=0.
               end if

*momentum used for the calculation of the width of the produced
*meson (in case of iwsc2.ge.1)
               pr=1.

               call vecmesm(srtfree,massmi-spot,massout(1),massout(2),
     &              idout(3),massout(3),j01,pr,spot)
*determine momenta
               call bopsmom(srtfree,massout(1),massout(2),massout(3)+
     &              spot,p3)
*boost
               do k=1,3
                  if(k.le.2) then
                     eout2(k)=sqrt(massout(k)**2+p3(k,1)**2+
     &                    p3(k,2)**2+p3(k,3)**2)
                  else
                     eout2(k)=sqrt((massout(k)+spot)**2+p3(k,1)**2+
     &                    p3(k,2)**2+p3(k,3)**2)
                  end if

                  call lorentz(-betacm(1),-betacm(2),-betacm(3),
     &                 p3(k,1),p3(k,2),p3(k,3),eout2(k))
               end do

               if(ipotvec.eq.2.and.(idout(3).eq.103.or.
     &              idout(3).eq.105)) then

*determine first "old" phase space + spectral function factor
*(with potential at p=0)
                  call bops3(pso,srtfree,rmass,pmass,massout(3)+spot)
*"old" spectral function
                  gamtot=bwmes((massout(3)+spot)**2,
     &                 idout(3)-idbmax,0,0,0,ratiom,0)
                  spectralo=2./pi*massout(3)**2*gamtot/
     &                 ((massout(3)**2-mesprop1(idout(3)-idbmax,1)**2)
     &                 **2+gamtot**2*massout(3)**2)
*"new" factors
                  pr=sqrt(p3(3,1)**2+p3(3,2)**2+p3(3,3)**2)
                  call vecspot(j01,pr,spot,idout(3)-idbmax)
                  call bops3(psn,srtfree,rmass,pmass,massout(3)+spot)
                  gamtot=bwmes((massout(3)+spot)**2
     &                 ,idout(3)-idbmax,0,0,0,ratiom,0)
                  spectraln=2./pi*massout(3)**2*gamtot/
     &                 ((massout(3)**2-mesprop1(idout(3)-idbmax,1)**2)
     &                 **2+gamtot**2*massout(3)**2)
**achtung
                  sigtot=sigtot*psn/pso*spectraln/spectralo
               end if
               
               do j=1,3
                  do k=1,3
                     pout(k,j)=p3(j,k)
                  end do
               end do
               potout(1)=0.
               potout(2)=0.
               potout(3)=spot
            else if(x7.le.(sigstt+sigvect)/sigtot) then
*strangeness production
               call strangeprod_photo(sigst,qnuk,pout,massout,potout,
     &              idout,izout,numparf)
               if(numparf.eq.0) then
                  write(*,*)'problems in strangeness production'
                  goto 100
               end if
            else if(leflag) then
*     low energy particle production
*     Auswuerfeln was produziert wird
               niter=0
               x5=rn(iseed)
 200           x6=rn(iseed)
               x3=rn(iseed)
               x4=rn(iseed)
               if(x5.le.sigmapi/sigtoto) then
                  sum=0.
                  do i3=-ntheta,ntheta-1
                     sum=sum+dsigma(i3)
                     if(x4.le.sum/sigmapi) then
                        sum2=0.
                        ka=min(int(probka(kmax,i3)/x6),1)+kmin
                        do i2=0,nres2
                           sum2=prob(i3,ka,i2)+sum2
                           if(x3.le.sum2) then
                              if(i2.eq.0) then
*     nucleon + pion wird eventuell produziert
                                 cost=(float(i3)+0.5)/float(ntheta)
                                 if(ipotcrosw.ge.4) then
                                    phi=phist(i3)
                                 else
                                    phi=rn(iseed)*2.*pi
                                 end if
*     determine scattering direction in lab frame
                                 if((pcny.ne.0.).or.(pcnx.ne.0.)) then 
                                    phik=atan2(pcny,pcnx)
                                 else
                                    phik=0.
                                 end if
                                 costk=pcnz/pcm
*     cm frame:
                                 pscattc(1)=cos(phi)*sqrt(1.-cost**2)
                                 pscattc(2)=sin(phi)*sqrt(1.-cost**2)
                                 pscattc(3)=cost
*     rotation to lab:
                                 call rotate(costk,phik,pscattc,
     &                                pscatt)
                                 
*charge assignment before siter (bec. of symmetrie pot.)
                                 if((ka.eq.1).or.(ka.eq.4)) then
                                    iz3=qnuk
                                 else
                                    iz3=abs(qnuk-1)
                                 end if

                                 if(debugflag) write(*,*)'vor ipot'
                                 if(ipotcrosw.ge.1) then
                                    if(debugflag) write(*,*)'vor siter'
                                    call siter(flagiter,xout)
                                    if(debugflag) write(*,*)'nach siter'
                                    if(.not.flagiter) then
                                       if(ipotcrosw.ge.4) then
                                          write(*,*)'problems in gamabs
     &                                         ipotcrosw,flagiter'
                                          stop
                                       end if
                                       niter=niter+1
                                       if(niter.ge.nitermax) then
                                          write(*,*)'problems in gamabs
     &                                         iteration not 
     &                                         succesful'
                                          numparf=0
                                          goto 100
                                       end if
                                       goto 200
                                    end if
                                 else
                                    xout=qpion
                                    enuk3=sqrt(rmass**2+xout**2)
                                 end if
*     Bestimmen des Impulses im lab-frame
                                 p1l=pscatt(1)*xout
                                 p2l=pscatt(2)*xout
                                 p3l=pscatt(3)*xout
                                 call lorentz(-betacm(1),-betacm(2),
     &                                -betacm(3),
     &                                p1l,
     &                                p2l,p3l,enuk3)
*     pion + nucleon wird produziert
                                 if(ipotcrosw.eq.0) then
                                    em3=rmass
                                    u3=0.
                                    id3=1
                                    em4=pmass
                                    u4=0.
                                    id4=1+idbmax
                                 end if


                                 numparf=2
                                 pertyp2=11

                                 idout(1)=id3
                                 idout(2)=id4
                                 izout(1)=iz3
                                 izout(2)=qnuk-iz3
                                 
                                 pout(1,1)=p1l
                                 pout(2,1)=p2l
                                 pout(3,1)=p3l
                                 potout(1)=u3
                                 massout(1)=em3
                                 
                                 pout(1,2)=px2-p1l
                                 pout(2,2)=py2-p2l
                                 pout(3,2)=pzt-p3l
                                 potout(2)=u4
                                 massout(2)=em4
*     check energy conservation                  
                                 etotal=sqrt((em3+u3)**2+p1l**2+
     &                                p2l**2
     &                                +p3l**2)+sqrt((em4+u4)**2+
     &                                (px2-p1l)**2+(py2-p2l)**2+(pzt-
     &                                p3l)**2)
                                 if(abs(etotal-ptot(0)).gt.1e-03) then
                                    write(*,*)'problems gamabs energy
     &                                   conservation',etotal,ptot(0),
     &                                   etotal-ptot(0)
                                 end if
                                 goto 25
                              else
*     Resonanz wird produziert
                                 goto 26
                              end if

                           end if
                        end do
                        write(*,*)'Fehler in gamabs x3,sum,prob'
                        stop
                     end if
                  end do
                  If(.not.onlyTwoPi) then
                    write(*,*) 'Fehler in gamabs x4,dsigma,sigmapi'
                    stop
                  end if
               else if(x5.le.(sigmapi+sigx)/sigtoto) then
*     Resonanz wird produziert
                  if(debugflag) write(*,*)'res prod'
                  sum=0.
                  do i2=1,nres2
                     sum=sum+sig2(i2)
                     if(x6.le.sum/sigx) then
 26                     continue
                        
                        if(debugflag) write(*,*)'vor id'
                        if(i2.eq.1) then
                           id3=2
                        else if(i2.eq.2) then
                           id3=7
                        else if(i2.eq.3) then
                           id3=4
                        else if(i2.eq.4) then
                           id3=16
                        end if
                        
                        iz3=qnuk
                        if(barprop2(id3,1).eq.1) then
                           em3=finmass(1)
                           u3=scapot(1)
                        else
                           em3=finmass(2)
                           u3=scapot(2)
                        end if

                        numparf=1
                        pertyp2=i2+6
                        idout(1)=id3
                        izout(1)=iz3
                        pout(1,1)=px2
                        pout(2,1)=py2
                        pout(3,1)=pzt
                        potout(1)=u3
                        massout(1)=em3
*     energy conservation
                        etotal=sqrt((em3+u3)**2+px2**2+py2**2+pzt**2)
                        if(abs(etotal-ptot(0)).gt.1e-03) then
                           write(*,*)'problems energy cons res',
     &                          etotal,
     &                          ptot(0)
                           write(*,*)em3,u3,finmass,scapot
                           write(*,*)srtfree,ptot(0),px2**2+py2**2+
     &                          pzt**2
                           write(*,*)egamma,px2,py2,pz2
                        end if
                        goto 25
                     end if         
                  end do
                  write(*,*)'Fehler in gamabs sum,x6,sigx'
                  stop
               else 
*     2 pi werden produziert
*     Auswuerfeln der Impulse
 238              continue
                  if(debugflag) write(*,*)'2pi production'
                  !  Auswuerfelnder Ladungen
                  x7=rn(iseed)
                  pr00=sig2pi(3)/sigx2
                  pr10=sig2pi(2)/sigx2
                  if(x7.lt.pr10) then
                     iz3=abs(qnuk-1)
                     iz4=0
                     iz5=qnuk-iz3-iz4
                     If (qnuk.eq.1) then
                        rCh=5
                     else
                        rCh=2
                     end if
                  else if(x7.lt.(pr10+pr00)) then
                     iz3=qnuk
                     iz4=0
                     iz5=0
                     If (qnuk.eq.1) then
                        rCh=3
                     else
                        rCh=6
                     end if
                  else
                     iz3=qnuk
                     iz4=1
                     iz5=-1
                     If (qnuk.eq.1) then
                        rCh=1
                     else
                        rCh=4
                     end if
                  end if
                  !Auswürfeln der Impulse

                  If(PascalTwoPi.and.(srtfree.lt.1.55)) then
                     !  rch=1: gamma p -> pi+ pi- p
                     !  rch=2: gamma p -> pi+ pi0 n
                     !  rch=3: gamma p -> pi0 pi0 p
                     !  rch=4: gamma n -> pi+ pi- n
                     !  rch=5: gamma n -> pi- pi0 p
                     !  rch=6: gamma n -> pi0 pi0 n      
                     pi1(0)=egamma
                     pi1(1)=0.
                     pi1(2)=0.
                     pi1(3)=pi1(0)
                     pi2(1)=p(1,i)
                     pi2(2)=p(2,i)
                     pi2(3)=p(3,i)
                     pi2(0)=sqrt(rmass**2+p(1,i)**2+p(2,i)**2+p(3,i)**2)
                     if(oliverdebug) then
                        write(*,*)'vor eventgenerator'
                        write(*,*)'.srtfree=',srtfree
                     End if
                     ! Check threshold:
                     if(srtfree.le.rmass+2*pmass) goto 100
                     call eventGenerator(rch,pi1,pi2,po1,po2,po3)
                     p3(1,1:3)=po1(1:3)
                     p3(2,1:3)=po2(1:3)
                     p3(3,1:3)=po3(1:3)
                     if(oliverdebug) write(*,*)'nach eventgenerator'
                  else
                     call bopsmom(srtfree,rmass,pmass,pmass,p3)
                     If (energieTwoPi.and.(ipotcrosw.ne.0)) then
                        ladung(1)=iz3
                        ladung(2)=iz4
                        ladung(3)=iz5
                        ort=r(:,i)
                        call DreierEnergie(p3,ptot,-betacm,ladung,Ort)
                                !Print *,"Nach DreierEnergie"
                                !stop
                     else
                        enuk=sqrt(rmass**2+p3(1,1)**2+p3(1,2)**2+p3(1,3)**2)
                        call lorentz(-betacm(1),-betacm(2),-betacm(3)
     &                       ,p3(1,1),p3(1,2),p3(1,3),enuk)
                                ! 2 pion + nucleon wird produziert
                                ! Lorentz-Trafo fuer Pionen
                        epion1=sqrt(pmass**2+p3(2,1)**2+p3(2,2)**2+p3(2,3)
     &                       **2)
                        call lorentz(-betacm(1),-betacm(2),-betacm(3),
     &                       p3(2,1),p3(2,2),
     &                       p3(2,3),epion1)

                        epion2=sqrt(pmass**2+p3(3,1)**2+p3(3,2)**2+p3(3,3)
     &                       **2)
                        call lorentz(-betacm(1),-betacm(2),-betacm(3),
     &                       p3(3,1),p3(3,2),
     &                       p3(3,3),epion2)
                                ! test for momentum conservation
                        if(ipotcrosw.eq.0) then
                           do k2=1,3
                              ptotal=0.
                              do k3=1,3
                                 ptotal=ptotal+p3(k3,k2)
                              end do
                              if(abs(ptotal-ptot(k2)).ge.1e-03) then
                                 write(*,*)'problems momentum conservation
     &                                gamabs',k2,p3,ptot
                              end if
                           end do
                                !energy conservation
                           etotal=enuk+epion1+epion2
                           if(abs(etotal-ptot(0)).gt.1e-03) then
                              write(*,*)'problems energy con in 2pi',
     &                             enuk,epion1,epion2,ptot(0),srtfree
                           end if
                        end if
                     end if
                  end if

                  
                  
                  em3=rmass
                  u3=0.
                  id3=1

                  em4=pmass
                  u4=0.
                  id4=1+idbmax

                  em5=pmass
                  u5=0.
                  id5=1+idbmax

                  numparf=3
                  pertyp2=12
                  idout(1)=id3
                  idout(2)=id4
                  idout(3)=id5
                  izout(1)=iz3
                  izout(2)=iz4
                  izout(3)=iz5
                  
                  massout(1)=em3
                  massout(2)=em4
                  massout(3)=em5
                  potout(1)=u3
                  potout(2)=u4
                  potout(3)=u5
                  
                  do k=1,3
                     do j=1,3
                        pout(k,j)=p3(j,k)
                     end do
                  end do
                  goto 25
               end if
 25            continue
            end if
            if(ishadow2.eq.1) then
               r22=sqrt(r(1,i)**2+r(2,i)**2+r(3,i)**2)
               costr=r(3,i)/r22
*sigs in mb! because of boffif
               sigs=xpom*srtfree**(2.*epspom)+
     &              yreg*srtfree**(-2.*etareg)
               shadfac=boffif(r22,costr,egamma,sigs)
               sigtot=sigtot*shadfac
*     write(*,*)'nach shadow',shadfac,sigtot
            end if
               
            if(numparf.eq.0) then
               write(*,*)'problems in gamabs numparf.eq.0'
               stop
            end if
*     check for pauli
            if(ipauli.ge.1) then
               if(debugflag) write(*,*)'vor pauli'
*               if(vecflag) write(*,*)'vor pauli'
               paulflag=.true.
               j=1
               ntag=0
               do while(paulflag)                     
                  if(idout(j).eq.1) then
                     eout=sqrt((massout(j)+potout(j))**2+
     &                    pout(1,j)**2+pout(2,j)**2+
     &                    pout(3,j)**2)
*                     if(vecflag) write(*,*)eout-massout(j),j01
                     call pauli(ntag,phase,r(1,i),r(2,i),r(3,i),
     &                    pout(1,j),pout(2,j),pout(3,j),eout) 
*                     if(vecflag) write(*,*)ntag
                  end if
                  j=j+1
                  if(j.gt.numparf.or.ntag.eq.-1) paulflag=.false.
               end do
*     check for pauli blocking
               if(ntag.eq.-1) then
                  goto 100
               end if
            end if

*     set ids of produced particles
            if(debugflag) write(*,*)'vor loop over produced'

            do j=1,numparf
               If ((numparf.ne.3).and.OnlyTwoPi) goto 100
               perflag=.false.
               if(idout(j).ne.1) then
                  perflag=.true.
               else
                  ener=sqrt(massout(j)**2+pout(1,j)**2+pout(2,j)**2
     &                 +pout(3,j)**2)
                  if(ener-rmass.gt.minenp) perflag=.true.
               end if
               if(perflag) then

                  ipert=ipert+1
                  ipert2=(irun-1)*maxper+ipert
*     important settings
                  if(debugflag) write(*,*)'vor setid',pout(1,j),
     &                 pout(2,j),pout(3,j)
                  call setidp(ipert2,pout(1,j),pout(2,j),pout(3,j),
     &                 massout(j),potout(j),idout(j),izout(j),i,0,
     &                 .false.,forflag,r(1,i),r(2,i),r(3,i))
                  perweight(ipert2)=sigtot
                  if (onlyTwoPi) perweight(ipert2)=sigx2

                  gidp(ipert2)=i
                  pertyp(ipert2)=pertyp2
                  if(pertyp2.eq.0) then
*                     write(*,*)'problems in gamabsper',pertyp2,
*     &                    numparf,idout(1),idout(2),idout(3)
                  end if

*cheating potential, which sets omegas on mass shell if they travel outside
*                  if(iwsc2.eq.3.and.idout(j).eq.105) then
*                     cheatpar(ipert2)=(ep(ipert2)-mesprop1(5,1))/j01
*                     ep(ipert2)=mesprop1(5,1)
*                     if(cheatpar(ipert2)*rho0+ep(ipert2)+
*     &                    upotp(ipert2).lt.3.*pmass) then
*                        cheatpar(ipert2)=(3.*pmass-ep(ipert2)-
*     &                       upotp(ipert2))/rho0
*                     end if
*                     upotp(ipert2)=upotp(ipert2)+cheatpar(ipert2)*j01
*                  end if
*fake potential for rho's and omega's
                  if(iwsc22.eq.1.and.(idout(j).eq.103.or.idout(j).
     &                 eq.105).and.vecflag) then
                     if(idout(j).eq.103) then
                        massmi=2.*pmass
                     else
                        massmi=3.*pmass
                     end if
                     if(numparf.eq.2) then
                        call resmasdi(srtfree,massmi,rmass,
     &                       -idout(j)+idbmax,massvac,potout(j))
                     else if(numparf.eq.3) then
                        call vecmesm(srtfree,massmi,rmass,
     &                       pmass,idout(j),massvac,0.,0.,potout(j))
                     end if
                     cheatpar(ipert2)=(ep(ipert2)-massvac)/j01
                     ep(ipert2)=massvac
**achtung
                     if(cheatpar(ipert2)*rho0+ep(ipert2)+
     &                    upotp(ipert2).lt.massmi) then
                        cheatpar(ipert2)=(massmi-ep(ipert2)-
     &                       upotp(ipert2))/rho0
                     end if

                     upotp(ipert2)=upotp(ipert2)+cheatpar(ipert2)*j01
                     if(ep(ipert2)+upotp(ipert2).lt.massmi) then
                        write(*,*)'problems in gamabsper',numparf,
     &                       massout(j),potout(j),ep(ipert2),
     &                       upotp(ipert2),massvac
                     end if
                  end if
*make sure that no rhos with effective mass lt 2*pmass are produced by fritiof
                  if(forflag.and.idout(j).eq.103) then
                     if(ipotvec.ge.1) then
                        pr=sqrt(pout(1,j)**2+pout(2,j)**2+
     &                       pout(3,j)**2)
                        call vecspot(j01,pr,spot,idout(j)-idbmax)
                     else
                        spot=0.
                     end if
                     if(massout(j)+spot.le.2.*pmass) then
                        ep(ipert2)=massout(j)-spot
                     end if
                  end if
**achtung
*                  if(idout(j).eq.105.and.massout(j).gt.0.50.and.
*     &                 massout(j).lt.0.55) then
*                     write(650,*)massout(j),sigtot,potout(j),numparf,
*     &                    srtfree,i
*                  end if

                  if(debugflag) write(*,*)'nach setid'
*     statistical settings
                  if(abs(idout(j)).le.idbmax) then
                     idp(ipert2,6)=2
                  else
                     do k5=1,3
                        rpiep(k5,ipert2)=r(k5,i)
                     end do

                     idp(ipert2,4)=1
                     if(numparf.eq.2.and.leflag) then
                        idp(ipert2,5)=-6
                     else if(numparf.eq.3.and.leflag) then
                        idp(ipert2,5)=-7
                     else
                        idp(ipert2,5)=-8
                     end if

                     rpiep(4,ipert2)=j01        
                     rpiep(6,ipert2)=1.0
                  end if
               end if
            end do
*     more statistics
            if(numparf.eq.1) then
               numpro((i2-1)*2+3+id(i,2))=numpro((i2-1)*2+3+
     &              id(i,2))+1
            else if(numparf.eq.2.and.leflag) then
               numpro(izout(2)+1)=numpro(izout(2)+1)+1
            else if(leflag) then
               numpro(11)=numpro(11)+1
            else if(vecflag) then
               numpro(12)=numpro(12)+1
            else
               numpro(13)=numpro(13)+1
            end if

            numcrt=numcrt+1
*            if(numcrt.ge.maxnum) goto 1000
            if(debugflag) write(*,*)'vor 100'
            sigabs=sigabs+sigtot
 100        continue
 1000    continue


      write(198,*)egamma,sigabs/float(numtry)
      
      write(199,*)egamma,masstot,numtry
      write(199,*)numpro(0),numpro(1),numpro(2)
      write(199,*)numpro(3),numpro(4)
      write(199,*)numpro(5),numpro(6)
      write(199,*)numpro(7),numpro(8)
      write(199,*)numpro(9),numpro(10)
      write(199,*)numpro(11),numpro(12),numpro(13)
*      write(122,*)'0.',numpro(0)+numpro(1)+numpro(2)+2*numpro(11),
*     &     '0','0','0'
*      write(117,*)'0.',numpro(3)+numpro(4),numpro(5)+numpro(6),
*     &     numpro(9)+numpro(10),numpro(7)+numpro(8)
      return


***********************************************************************
*                                                                     *
      entry gamausper(isu,isu2)           
*                                                                     *
*  purpose: folds yield of final mesons with propabilities of         *
*           initial meson or resonance production                     *
*                                                                     *
*   input variables: mass   (mass of target)                          *
*                    num    (number of parallel runs)                 *
*                    ...                                              *
*   output: rest of them                                              *
*                                                                     *
*ioutput=0: nur cross sections,ohne fort.300
*       =1: alles
*       =2: nur Mesonen
*       =3: nur Nukleonen
*       =4: nur 2-Nukleonen
***********************************************************************
      write(*,*)'in gamausper',numtry
      
      if(isu2.eq.1) then
        numtryt=0
        ifra=0
        do i=1,nmes
           do k2=-1,1
              num2(i,k2)=0
              tsigma(i,k2)=0.
              do j=0,angbins-1
                 asigma(i,k2,j)=0.
                 a2sigma(i,k2,j)=0.
              end do
              do j=0,maxpt
                 msigma(i,k2,j)=0.
              end do
              do j=0,maxtk
                 tksigma(i,k2,j)=0.
              end do
           end do
        end do
        do i=-3,2
           t2sigma(i)=0.
           do j=0,nmass
              t2sigdm(i,j)=0.
              t2sigdmr(i,j)=0.
           end do
        end do
        do i=1,10
           do j=1,nmes
              do k2=-1,1
                 numcre(i,j,k2)=0
              end do
           end do
        end do
        do i=0,maxtk
           dprodt(i)=0.
        end do
      end if
      if(ioutput.gt.0) then
        if(isu2.eq.1) then
          write(300,*)egamma,masstot,num,isubs,ioutput
        end if
        write(300,*)isu,numtry
      end if
*      write(*,*)'vor loop'
* --- loop over parallel runs:                                         
      do irun=1,num
         anzp=0
         ibeg=(irun-1)*maxper
         do i=ibeg+1,ibeg+maxper
            if(idp(i,1).eq.0)  goto 600
            if(idp(i,1).gt.idbmax) then
*mesons
               idm=idp(i,1)-idbmax
               idch=idp(i,2)
               if(idm.gt.nmes) then
                  write(*,*)'Fehler in gamout id',idp(i,1)
                  stop
               end if
               ptrans=sqrt(pp(1,i)**2+pp(2,i)**2)
*****************************
*     angular distribution
               theta=atan2(ptrans,pp(3,i))
               nth=int(theta/pi*angbins)
*     write(*,*)'ptrans,theta,nth',ptrans,theta,nth
*     impulse distribution
               ptot2=sqrt(pp(1,i)**2+pp(2,i)**2+pp(3,i)**2)
               npt=min(int(ptot2/dpt),maxpt)
               if(npt.eq.maxpt) then
                  write(*,*)'gamabsper,maxpt too small',ptot2,dpt,
     &                 maxpt
               end if
*     kinetic energy distribution
               tkin=sqrt(ep(i)**2+ptot2**2)-ep(i)
               ntk=min(int(tkin/dtk),maxtk)

               num2(idm,idch)=num2(idm,idch)+1

               tsigma(idm,idch)=tsigma(idm,idch)+
     &              perweight(i)/float(numtry*isubs)
*     write(*,*)tsigma(ii)
               asigma(idm,idch,nth)=perweight(i)/(2.*pi*
     &              (cos(float(nth)/
     &              angbins*
     &              pi)-cos(float(nth+1)/angbins*pi)))/float(numtry*
     &              isubs)+
     &              asigma(idm,idch,nth)
               a2sigma(idm,idch,nth)=perweight(i)/float(numtry*isubs)/
     &              (pi/float(angbins))+a2sigma(idm,idch,nth)
*     write(*,*)asigma(ii,nth)
               msigma(idm,idch,npt)=perweight(i)/dpt/
     &              float(numtry*isubs)+
     &              msigma(idm,idch,npt)
*     write(*,*)msigma(ii,npt)
               tksigma(idm,idch,ntk)=perweight(i)/dtk/
     &              float(numtry*isubs)+
     &              tksigma(idm,idch,ntk)
*     Absorption distribution
               if(idp(i,4).gt.maxncre) maxncre=idp(i,4)
               inumcre=min(idp(i,4),10)
               if(inumcre.eq.0) then
                  write(*,*)'problems gamaus inumcre',idp(i,4),
     &                 idp(i,5),idp(i,6),idp(i,1)
                  inumcre=1
               end if
               numcre(inumcre,idm,idch)=numcre(inumcre,idm,idch)+1
*     two pion production
               if(idm.eq.1) then
*     
                  anzp=anzp+1
                  if(anzp.gt.maxanzp) then
                     write(*,*)'problems in gamabs: too many pions'
                     write(*,*)irun,idp(i,2),pp(1,i),pp(2,i),pp(3,i)
                     stop
                  end if 
                  
                  ipis(anzp,2)=idp(i,2)
                  ipis(anzp,3)=idp(i,3)
                  ipis(anzp,6)=idp(i,6)
                  do i5=1,3
                     ppis(i5,anzp)=pp(i5,i)
                  end do
                  gidps(anzp)=gidp(i)
                  pervs(anzp)=perweight(i)
                  eps(anzp)=ep(i)
                  if(anzp.ge.2) then
*     sum up cross sectiona
                     do i5=1,anzp-1
                        if(gidps(i5).eq.gidp(i)) then
*     check
                           if(abs(perweight(i)-pervs(i5)).gt.1e-04)
     &                          then
                              write(*,*)'problems gamaus',gidps(i5),
     &                             gidp(i),perweight(i),pervs(i5)
                              stop
                           end if
                           iii=idp(i,2)+ipis(i5,2)
                           if(iii.eq.0.and.idp(i,2).ne.0) then
*     pi+ pi-
                              iii=-3
                           end if
                           t2sigma(iii)=t2sigma(iii)+perweight(i)/
     &                          float(numtry*isubs)
*invariant mass distribution
                           ener1=sqrt(eps(i5)**2+ppis(1,i5)**2+
     &                          ppis(2,i5)**2+ppis(3,i5)**2)
                           ener2=sqrt(ep(i)**2+pp(1,i)**2+pp(2,i)**2+
     &                          pp(3,i)**2)
                           invmass=sqrt((ener1+ener2)**2-
     &                          (ppis(1,i5)+pp(1,i))**2-
     &                          (ppis(2,i5)+pp(2,i))**2-
     &                          (ppis(3,i5)+pp(3,i))**2)
                           imass=max(min(int((invmass-2.*pmass)/dmass)
     &                          ,nmass),0)
                           t2sigdm(iii,imass)=t2sigdm(iii,imass)+
     &                          perweight(i)/dmass/float(numtry*isubs)
                           if(idp(i,3).eq.ipis(i5,3).and.
     &                          idp(i,6).eq.103.and.
     &                          ipis(i5,6).eq.103) then
                              t2sigdmr(iii,imass)=t2sigdmr(iii,imass)+
     &                          perweight(i)/dmass/float(numtry*isubs)
*                              if(abs(rpiep(5,i)-invmass).gt.1e-03) then
*                                 write(*,*)'problems gamaus'
*                                 write(*,*)'masses',i,i5,invmass,
*     &                             rpiep(5,i),idp(i,3)
*                              end if
                           end if
                        end if
                     end do
                  end if
               end if
            else
*bayrons
               if(idp(i,1).eq.1.and.idp(i,2).eq.1) then
*proton
                  ptot2=sqrt(pp(1,i)**2+pp(2,i)**2+pp(3,i)**2)
*     kinetic energy distribution
                  tkin=sqrt(ep(i)**2+ptot2**2)-ep(i)
                  ntk=min(int(tkin/dtk),maxtk)
                  dprodt(ntk)=dprodt(ntk)
     &                 +perweight(i)/dtk/float(numtry*isubs)
               end if
            end if
* -    -(end of loop over all particles)-
 600     end do
 
      end do
      numtryt=numtryt+numtry

      write(230,*) egamma,isu,numtry
      if(isu2.ne.isubs) return
      write(210,301)egamma,tsigma(1,-1),tsigma(1,0),tsigma(1,1),
     &     tsigma(2,0),tsigma(10,1),tsigma(11,-1)
 301  format(7(e12.4,1x))
      write(200,*)egamma,masstot,numtryt,num*isubs,ifra            
      write(200,*)tsigma(1,-1),tsigma(1,0),tsigma(1,1),tsigma(2,0)
      write(200,*)num2(1,-1),num2(1,0),num2(1,1),num2(2,0)  
      write(201,*)'#',egamma,masstot,num
      write(209,*)'#',egamma,masstot,num
      do i=0,angbins-1
        write(201,302)(float(i)+0.5)*180./angbins,asigma(1,-1,i),
     &        asigma(1,0,i),asigma(1,1,i),asigma(2,0,i),
     &        asigma(10,1,i),asigma(11,-1,i)
        write(209,302)(float(i)+0.5)*180./angbins,a2sigma(1,-1,i),
     &       a2sigma(1,0,i),a2sigma(1,1,i),a2sigma(2,0,i),
     &       a2sigma(10,1,i),a2sigma(11,-1,i)
 302    format(7(e12.4,1x))
      end do
      write(202,*)'#',egamma,masstot,num
      do i=0,maxpt
        write(202,302)(float(i)+0.5)*dpt,msigma(1,-1,i),
     &        msigma(1,0,i),msigma(1,1,i),msigma(2,0,i),
     &        msigma(10,1,i),msigma(11,-1,i)
      end do
      write(204,*)egamma,masstot,numtryt,maxncre
      do i=1,10
        write(204,*)i,numcre(i,1,-1),numcre(i,1,0),numcre(i,1,1),
     &        numcre(i,2,0)
      end do
      write(208,*)'#',egamma,masstot,num
      do i=0,maxtk
         write(208,303) (float(i)+0.5)*dtk,tksigma(1,-1,i),
     &        tksigma(1,0,i),tksigma(1,1,i),
     &        tksigma(2,0,i),tksigma(10,1,i),tksigma(11,-1,i),
     &        dprodt(i)
      end do
 303  format(8(e12.4,1x))
      write(240,*)egamma,t2sigma(-3),t2sigma(0),t2sigma(-3)+
     &     t2sigma(-2)+t2sigma(-1)+t2sigma(0)+t2sigma(1)+t2sigma(2)
      write(241,*)'#',egamma,masstot,num
      do i=0,nmass
         write(241,304)(float(i)+0.5)*dmass+2.*pmass,
     &        (t2sigdm(j,i),j=-3,2)
      end do
      write(241,*)
      do i=0,nmass
         write(241,304)(float(i)+0.5)*dmass+2.*pmass,
     &        (t2sigdmr(j,i),j=-3,2)
      end do
 304  format(7(e12.4,1x))
      return                                                    
      end

