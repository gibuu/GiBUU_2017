c=================================================================
c MomPart.F:   Provide Some Routines To Handle /VarsMomPart/
c              Kai Gallmeister, 19.10.2001...28.2.2005
c=================================================================
c
c provided routines:
c ==================
c
c subroutine MP_WriteVersion(no)
c
c subroutine MP_Set4(i, m, x,y,z,E)
c subroutine MP_Set3(i, m, x,y,z)
c
c double precision function MP_M(iPart)
c double precision function MP_P(iPart,i)
c
c subroutine MP_Write(no, iMin,iMax)
c subroutine MP_Copy(i1,i2)
c subroutine MP_CopyMom(i1,i2,fak)
c subroutine MP_AddMom(i1,i2,fak)
c subroutine MP_V2Vars(i, mT,pT,y,phi)
c
c subroutine MP_CalcROBO(i1,i2, sqrts, theta, phi, beta)
c subroutine MP_ROBO(iMin,iMax,theta,phi,beX,beY,beZ)
c subroutine MP_ROBO_Inv(iMin,iMax,theta,phi,beX,beY,beZ)
c
c double precision function MP_SqrtS(i1,i2)
c double precision function MP_ScalProd4(i1,i2)
c double precision function MP_ScalProd3(i1,i2)
c double precision function MP_CosAngle3(i1,i2)
c
c double precision function MP_Datanh2(x1,x2)
c
c subroutine MP_TransformLL(i1,i2,o1,o2)
c
c-----------------------------------------------------------------
c uses the internal arrays MP(4,nVarsMomPart) and MPm(nVarsMomPart),
c while the last entry is reserved for internal use.
c The user is free to use all the lines 1...(nVarsMomPart-1) for
c his own purposes.
c-----------------------------------------------------------------
c
c 28.02.05: v1.1, TransformLL

c=================================================================
c Write out the current version

      subroutine MP_WriteVersion(no)
      IMPLICIT NONE
      integer no


      integer nVarsMomPart
      parameter (nVarsMomPart = 51)
      COMMON /VarsMomPart/
     $     MP(4,nVarsMomPart),  ! (direction [x,y,z,E],particle)
     $     MPm(nVarsMomPart)    ! mass of particle
      double precision MP, MPm
      SAVE /VarsMomPart/

      write(no,*) '** MomPart v1.1, 28.02.2005 **',
     $     nVarsMomPart-1,' arrays'
      end

c=================================================================
c Set the contents for particle "i", energy explicitly given

      subroutine MP_Set4(i, m, x,y,z,E)
      IMPLICIT NONE
      integer i
      double precision x,y,z,E,m


      integer nVarsMomPart
      parameter (nVarsMomPart = 51)
      COMMON /VarsMomPart/
     $     MP(4,nVarsMomPart),  ! (direction [x,y,z,E],particle)
     $     MPm(nVarsMomPart)    ! mass of particle
      double precision MP, MPm
      SAVE /VarsMomPart/

      if (i.ge.nVarsMomPart) then
         write(*,*) 'MP_Set: iLine=',i,'>',nVarsMomPart-1,'! Stop!'
         stop
      endif

      MP(1,i) = x
      MP(2,i) = y
      MP(3,i) = z
      MP(4,i) = E
      MPm(i) = m
      
      end

c=================================================================
c Set the contents for particle "i", energy calculated

      subroutine MP_Set3(i, m, x,y,z)
      IMPLICIT NONE
      integer i
      double precision x,y,z,m

      call MP_Set4(i,m,x,y,z,sqrt(x**2+y**2+z**2+m**2))
      
      end

c=================================================================
c Get the Mass of particle "iPart"

      double precision function MP_M(iPart)

      integer nVarsMomPart
      parameter (nVarsMomPart = 51)
      COMMON /VarsMomPart/
     $     MP(4,nVarsMomPart),  ! (direction [x,y,z,E],particle)
     $     MPm(nVarsMomPart)    ! mass of particle
      double precision MP, MPm
      SAVE /VarsMomPart/
      MP_M = MPm(iPart)
      return
      end

c=================================================================
c Get the "i"th Momentum Component of particle "iPart"

      double precision function MP_P(iPart,i)

      integer nVarsMomPart
      parameter (nVarsMomPart = 51)
      COMMON /VarsMomPart/
     $     MP(4,nVarsMomPart),  ! (direction [x,y,z,E],particle)
     $     MPm(nVarsMomPart)    ! mass of particle
      double precision MP, MPm
      SAVE /VarsMomPart/
      MP_P = MP(i,iPart)
      return
      end

c=================================================================
c Write the contents for particle "iMin" to "iMax" to file "no"

      subroutine MP_Write(no, iMin,iMax)
      IMPLICIT NONE
      integer iMin,iMax, no


      integer nVarsMomPart
      parameter (nVarsMomPart = 51)
      COMMON /VarsMomPart/
     $     MP(4,nVarsMomPart),  ! (direction [x,y,z,E],particle)
     $     MPm(nVarsMomPart)    ! mass of particle
      double precision MP, MPm
      SAVE /VarsMomPart/

      character*5 T(9)
      integer i,j

      double precision eps,x(9)
      parameter (eps = 1e-15)

      data T /'x','y','z','E','m', 'mT','pT','y','phi'/

      write(no,11) (T(i),i=1,9)
      do i=iMin,iMax
         do j=1,4
            x(j) = MP(j,i)
         enddo
         x(5) = MPm(i)
         call MP_V2Vars(i, x(6),x(7),x(8),x(9))

         do j=1,9
            if (abs(x(j)).lt.eps) x(j) = 0.0
         enddo

         write(no,10) i, (x(j),j=1,9)
      enddo

 10   FORMAT (i2.2,': ',5(1P,e12.4),' -> ',4(1P,e12.4))
 11   FORMAT ('    ',5A12,' -> ',4A12)


      end

c=================================================================
c Copy Particle i1 to Particle i2 in the MP-Array

      subroutine MP_Copy(i1,i2)
      IMPLICIT NONE
      integer i1,i2


      integer nVarsMomPart
      parameter (nVarsMomPart = 51)
      COMMON /VarsMomPart/
     $     MP(4,nVarsMomPart),  ! (direction [x,y,z,E],particle)
     $     MPm(nVarsMomPart)    ! mass of particle
      double precision MP, MPm
      SAVE /VarsMomPart/

      integer j
      do j=1,4
         MP(j,i2) = MP(j,i1)
      enddo
      MPm(i2) = MPm(i1)
      end

c=================================================================
c Copy Particle i1 to Particle i2 in the MP-Array

      subroutine MP_CopyMom(i1,i2,fak)
      IMPLICIT NONE
      integer i1,i2
      double precision fak


      integer nVarsMomPart
      parameter (nVarsMomPart = 51)
      COMMON /VarsMomPart/
     $     MP(4,nVarsMomPart),  ! (direction [x,y,z,E],particle)
     $     MPm(nVarsMomPart)    ! mass of particle
      double precision MP, MPm
      SAVE /VarsMomPart/

      integer j
      do j=1,4
         MP(j,i2) = MP(j,i1)*fak
      enddo
      end

c=================================================================
c Add Momenta of Particle i1 to Particle i2 in the MP-Array
c (multiply them with factor fak before adding) 

      subroutine MP_AddMom(i1,i2,fak)
      IMPLICIT NONE
      integer i1,i2
      double precision fak


      integer nVarsMomPart
      parameter (nVarsMomPart = 51)
      COMMON /VarsMomPart/
     $     MP(4,nVarsMomPart),  ! (direction [x,y,z,E],particle)
     $     MPm(nVarsMomPart)    ! mass of particle
      double precision MP, MPm
      SAVE /VarsMomPart/

      integer j
      do j=1,4
         MP(j,i2) = MP(j,i2) + fak * MP(j,i1)
      enddo
      end

c=================================================================
c Calculate some kinematical variables

      subroutine MP_V2Vars(i, mT,pT,y,phi)
      IMPLICIT NONE
      double precision mT,pT,y,phi
      integer i


      integer nVarsMomPart
      parameter (nVarsMomPart = 51)
      COMMON /VarsMomPart/
     $     MP(4,nVarsMomPart),  ! (direction [x,y,z,E],particle)
     $     MPm(nVarsMomPart)    ! mass of particle
      double precision MP, MPm
      SAVE /VarsMomPart/

      double precision MP_Datanh2

      pT  = sqrt(MP(1,i)**2 + MP(2,i)**2)
      phi = atan2( MP(2,i), MP(1,i) )
      y   = MP_Datanh2( MP(3,i), MP(4,i) )
      mT  = sqrt(max(0d0,MP(4,i)**2 - MP(3,i)**2))

      end

c=================================================================
c Calculate the Params. for a Poincare-Transform to the CMS
c
c if (i1!=i2), then the params are for the rest system of the pair,
c otherwise it iss for the particle i1.
c (also sqrt(s) is calculated).
c
c theta, phi and beta are given in such a way, that calling ROBO_MP
c with these parameters for particles in some rest-frame, the momenta
c of these particles are afterwards in (i1,i2) frame.
c i.e.: if i1,i2 are given in the Lab frame, you can create some
c particles in the CM frame and call ROBO_MP(...,theta,phi,beta) to
c transform these created particles in the original lab frame
c
c everythin clear? "get the power, read the source!" - especially the
c comments and the out-commented commands!

      subroutine MP_CalcROBO(i1,i2, sqrts, theta, phi, beta)
      IMPLICIT NONE
      integer i1,i2
      double precision sqrts, theta, phi, beta(3)


      integer nVarsMomPart
      parameter (nVarsMomPart = 51)
      COMMON /VarsMomPart/
     $     MP(4,nVarsMomPart),  ! (direction [x,y,z,E],particle)
     $     MPm(nVarsMomPart)    ! mass of particle
      double precision MP, MPm
      SAVE /VarsMomPart/

      double precision SHR
      integer i

c calculate beta, sqrts:

      if (i1.ne.i2) then        ! 2 particle system

         SHR = MP(4,i1)+MP(4,i2)
         sqrts = SHR**2
         do i=1,3
            beta(i) = MP(i,i1)+MP(i,i2)
            sqrts = sqrts-beta(i)**2
            beta(i) = beta(i)/SHR
         enddo
         sqrts=sqrt(sqrts)

      else                      ! 1 particle system

         sqrts = MPm(i1)
         do i=1,3
            beta(i) = MP(i,i1)/MP(4,i1)
         enddo      

      endif

c      write(6,*) 'Start:'
c      call MP_Write(6, 1,2)

      CALL MP_Copy(i1,nVarsMomPart) ! save the original vector

      CALL MP_ROBO(i1,i1, 0d0,0d0, -beta(1),-beta(2),-beta(3))
c      call MP_Write(6, i1,i1)

      phi = atan2( MP(2,i1), MP(1,i1))

      CALL MP_ROBO(i1,i1, 0d0,-phi, 0d0,0d0,0d0)
c      call MP_Write(6, i1,i1)
      theta = atan2( MP(1,i1), MP(3,i1) )

c      CALL MP_ROBO(i1,i1, -theta,0d0, 0d0,0d0,0d0)
c      call MP_Write(6, i1,i1)

c      write(6,1000) 'cms(A,B): ',theta,phi,beta(1),beta(2),beta(3)
c      call MP_Write(6, i1,i1)

c      write(6,*) 'Back:'
c      CALL MP_ROBO(i1,i1, theta,phi, beta(1),beta(2),beta(3))
c      call MP_Write(6, i1,i1)

      CALL MP_Copy(nVarsMomPart,i1) ! restore the original vector
c      write(6,*) 'Restored:'
c      call MP_Write(6, i1,i1)

 1000 FORMAT(A,6D11.3)

      end

c=================================================================
c Performs rotations and boosts.
c (Taken from PYTHIA)
c
c Does the RoBo for all particles from "iMin" to "iMax"

      subroutine MP_ROBO(iMin,iMax,theta,phi,beX,beY,beZ)
      IMPLICIT NONE
      integer iMin,iMax
      double precision theta,phi,beX,beY,beZ


      integer nVarsMomPart
      parameter (nVarsMomPart = 51)
      COMMON /VarsMomPart/
     $     MP(4,nVarsMomPart),  ! (direction [x,y,z,E],particle)
     $     MPm(nVarsMomPart)    ! mass of particle
      double precision MP, MPm
      SAVE /VarsMomPart/

      double precision ROT(3,3),PR(3),DP(4)
      double precision DBX,DBY,DBZ, DB, EPS1, DGA,DBP,DGABP
      integer i,j

      parameter (EPS1=1D0-1D-12)

C...Rotate, typically from z axis to direction (theta,phi).

      if (theta**2+phi**2 .gt. 1D-20) then
         ROT(1,1)=cos(theta)*cos(phi)
         ROT(1,2)=-sin(phi)
         ROT(1,3)=sin(theta)*cos(phi)
         ROT(2,1)=cos(theta)*sin(phi)
         ROT(2,2)=cos(phi)
         ROT(2,3)=sin(theta)*sin(phi)
         ROT(3,1)=-sin(theta)
         ROT(3,2)=0D0
         ROT(3,3)=cos(theta)
         do i=iMin,iMax
            do j=1,3
               PR(j)=MP(j,i)
            enddo
            do j=1,3
               MP(j,i)=ROT(j,1)*PR(1)+ROT(j,2)*PR(2)+ROT(j,3)*PR(3)
            enddo
         enddo
      endif
 
C...Boost, typically from rest to momentum/energy=beta.

      DB = beX**2+beY**2+beZ**2
      if (DB .gt. 1D-20) then
         DBX=beX
         DBY=beY
         DBZ=beZ
         DB=SQRT(DB)

         if (DB.gt.EPS1) then
C...  Rescale boost vector if too close to unity.
            write (*,*) '(MP_ROBO:) boost vector too large'
            DBX=DBX*(EPS1/DB)
            DBY=DBY*(EPS1/DB)
            DBZ=DBZ*(EPS1/DB)
            DB=EPS1
         endif

         DGA=1D0/SQRT(1D0-DB**2)
         do i=IMIN,IMAX
            do j=1,4
               DP(j)=MP(j,i)
            enddo
            DBP=DBX*DP(1)+DBY*DP(2)+DBZ*DP(3)
            DGABP=DGA*(DGA*DBP/(1D0+DGA)+DP(4))
            MP(1,i)=DP(1)+DGABP*DBX
            MP(2,i)=DP(2)+DGABP*DBY
            MP(3,i)=DP(3)+DGABP*DBZ
            MP(4,i)=DGA*(DP(4)+DBP)
         enddo
      endif
      
      return
      end
c=================================================================
c Performs inverse boosts and inverse rotation.
c
c Does the inverse RoBo for all particles from "iMin" to "iMax"

      subroutine MP_ROBO_Inv(iMin,iMax,theta,phi,beX,beY,beZ)
      IMPLICIT NONE
      integer iMin,iMax
      double precision theta,phi,beX,beY,beZ

      call MP_ROBO(iMin,iMax,0d0,0d0,-beX,-beY,-beZ)
      call MP_ROBO(iMin,iMax,0d0,-phi,0d0,0d0,0d0)
      call MP_ROBO(iMin,iMax,-theta,0d0,0d0,0d0,0d0)
      end

c=================================================================
c Calculate the Invariant Mass of one/two 4-Vectors

      double precision function MP_SqrtS(i1,i2)
      IMPLICIT NONE
      integer i1,i2
  

      integer nVarsMomPart
      parameter (nVarsMomPart = 51)
      COMMON /VarsMomPart/
     $     MP(4,nVarsMomPart),  ! (direction [x,y,z,E],particle)
     $     MPm(nVarsMomPart)    ! mass of particle
      double precision MP, MPm
      SAVE /VarsMomPart/
 
      if (i1.ne.i2) then
         MP_SqrtS =sqrt( (MP(4,i1)+MP(4,i2))**2  
     $        -((MP(1,i1)+MP(1,i2))**2+(MP(2,i1)+MP(2,i2))**2
     $        +(MP(3,i1)+MP(3,i2))**2) )
      else
         MP_SqrtS =sqrt( MP(4,i1)**2  
     $        -(MP(1,i1)**2+MP(2,i1)**2+MP(3,i1)**2) )
      endif
      return
      end

c=================================================================
c Calculate the scalar product of two 4-Vectors

      double precision function MP_ScalProd4(i1,i2)
      IMPLICIT NONE
      integer i1,i2
  

      integer nVarsMomPart
      parameter (nVarsMomPart = 51)
      COMMON /VarsMomPart/
     $     MP(4,nVarsMomPart),  ! (direction [x,y,z,E],particle)
     $     MPm(nVarsMomPart)    ! mass of particle
      double precision MP, MPm
      SAVE /VarsMomPart/
 
      MP_ScalProd4 = MP(4,i1)*MP(4,i2)
     $     -(MP(1,i1)*MP(1,i2)+MP(2,i1)*MP(2,i2)+MP(3,i1)*MP(3,i2))
      return
      end

c=================================================================
c Calculate the scalar product of two 3-Vectors

      double precision function MP_ScalProd3(i1,i2)
      IMPLICIT NONE
      integer i1,i2
  

      integer nVarsMomPart
      parameter (nVarsMomPart = 51)
      COMMON /VarsMomPart/
     $     MP(4,nVarsMomPart),  ! (direction [x,y,z,E],particle)
     $     MPm(nVarsMomPart)    ! mass of particle
      double precision MP, MPm
      SAVE /VarsMomPart/
 
      MP_ScalProd3 =
     $     MP(1,i1)*MP(1,i2)+MP(2,i1)*MP(2,i2)+MP(3,i1)*MP(3,i2)
      return
      end

c=================================================================
c Calculate the Cosine of the angle between two 3-Vectors

      double precision function MP_CosAngle3(i1,i2)
      IMPLICIT NONE
      integer i1,i2
  

      integer nVarsMomPart
      parameter (nVarsMomPart = 51)
      COMMON /VarsMomPart/
     $     MP(4,nVarsMomPart),  ! (direction [x,y,z,E],particle)
     $     MPm(nVarsMomPart)    ! mass of particle
      double precision MP, MPm
      SAVE /VarsMomPart/
      
      double precision MP_ScalProd3
 
      MP_CosAngle3 =  MP_ScalProd3(i1,i2)
     $     /sqrt(MP_ScalProd3(i1,i1)*MP_ScalProd3(i2,i2))
      return
      end

c=================================================================

      double precision function MP_Datanh2(x1,x2)
      IMPLICIT NONE
      double precision x1,x2

      double precision eps,res, d,s,h

      PARAMETER (eps=1d-10, res=11.8594991d0 ) 
                                ! res = 0.5d0*log((2d0-eps)/eps)

      d = x2-x1
      s = x2+x1
      h = x2*eps

      if (x2.ge.0d0) then

         if (d.le.h) then
            MP_Datanh2 = res      ! here also x1==x2==0 !!!
         else if (s.le.h) then
            MP_Datanh2 = -res
         else
            MP_Datanh2 = 0.5d0*log(s/d)
         endif 

      else

         if (d.ge.h) then
            MP_Datanh2 = res
         else if (s.ge.h) then
            MP_Datanh2 = -res
         else
            MP_Datanh2 = 0.5d0*log(s/d)
         endif 

      endif
      return
      end

c=================================================================
c Convert vectors i1,i2 to lightcone vectors o1,o2

      subroutine MP_TransformLL(i1,i2,o1,o2)
      IMPLICIT NONE
      integer i1,i2, o1,o2


      integer nVarsMomPart
      parameter (nVarsMomPart = 51)
      COMMON /VarsMomPart/
     $     MP(4,nVarsMomPart),  ! (direction [x,y,z,E],particle)
     $     MPm(nVarsMomPart)    ! mass of particle
      double precision MP, MPm
      SAVE /VarsMomPart/

      double precision DHKC, DHKS, DHK1,DHK2, DHM1,DHM2
      integer j

      double precision MP_M     ! prototype
      double precision MP_ScalProd4 ! prototype 


      DHM1 = MP_M(i1)**2
      DHM2 = MP_M(i2)**2
      DHKC = MP_ScalProd4(i1,i2)
      DHKS = sqrt(DHKC**2-DHM1*DHM2)
      DHK1 = 0.5D0*((DHM2+DHKC)/DHKS-1D0)
      DHK2 = 0.5D0*((DHM1+DHKC)/DHKS-1D0)

      do j=1,4
         MPm(o1) = 0d0
         MPm(o2) = 0d0
         
         MP(j,o1) = (1d0+DHK1)*MP(j,i1) - DHK2*MP(j,i2)
         MP(j,o2) = (1d0+DHK2)*MP(j,i2) - DHK1*MP(j,i1)
      enddo


      end


c=================================================================
