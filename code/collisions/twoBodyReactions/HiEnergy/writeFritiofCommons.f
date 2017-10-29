      subroutine writeFritiofCommons(iFile)

      implicit none
      integer iFile


      integer KSZJ,KSZ1,KSZ2

      PARAMETER (KSZJ=4000,KSZ1=20, KSZ2=300)

      COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)
      integer KFR
      real VFR
      SAVE /FRPARA1/

      COMMON/FRINTN0/PLI0(2,4),AOP(KSZ1),IOP(KSZ1),NFR(KSZ1)
      integer IOP,NFR
      real PLI0,AOP
      SAVE /FRINTN0/

      COMMON/FRINTN3/IDN(2,KSZ2),FMN(2,KSZ2),NUC(2,3000)
      integer IDN,NUC
      real FMN
      SAVE /FRINTN3/

      COMMON/FRCODES/IPT(2),PACD(27),NNUC(27),NPROT(27),KCD(27)
     >           ,RO1(27,2),EXMA(9,2)
      integer IPT,NNUC,NPROT,KCD
      real RO1,EXMA
      character PACD*4
      SAVE /FRCODES/
      
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      integer MSTU,MSTJ
      real PARU,PARJ
      SAVE /LUDAT1/


      integer i

      integer iCall
      save iCall
      data iCall /0/

      write(iFile+iCall,1000) 'FRPARA1, FRINTN0'
      write(iFile+iCall,1000)
     $     'KFR(i),VFR(i),IOP(i),NFR(i),AOP(i)'
      do i=1,KSZ1
         write(iFile+iCall,1001) i,KFR(i),VFR(i),IOP(i),NFR(i),AOP(i)
      end do

      write(iFile+iCall,1000)
     $     'LUDAT1:MSTU(200),PARU(200),MSTJ(200),PARJ(200)'
      do i=1,200
         write(iFile+iCall,1002) i,MSTU(i),PARU(i),MSTJ(i),PARJ(i)
      enddo


      iCall = iCall+1


 1000 FORMAT (//,'=== ',A,' ===')
 1001 FORMAT (i5,':',i12,g12.5,2i12,g12.5)
 1002 FORMAT (i5,':',i12,g12.5,i12,g12.5)

      end

