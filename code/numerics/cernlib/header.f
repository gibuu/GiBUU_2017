***************************************************************************************************************************
* This subroutine is part of CERNLIB
* (C) Copyright CERN except where explicitly stated otherwise. 
* Permission to use and/or redistribute this work is granted under the terms of the GNU General Public License. 
***************************************************************************************************************************
*
* $Id: header.F,v 1.1.1.1 1996/04/01 15:01:12 mclareni Exp $
*
* $Log: header.F,v $
* Revision 1.1.1.1  1996/04/01 15:01:12  mclareni
* Mathlib gen
*
*
      SUBROUTINE HEADER(CODE,MODE)
C     This routine prints a page header with a testing routine name
C     message.
      CHARACTER*(*) CODE
*
* $Id: iorc.inc,v 1.1.1.1 1996/04/01 15:01:31 mclareni Exp $
*
* $Log: iorc.inc,v $
* Revision 1.1.1.1  1996/04/01 15:01:31  mclareni
* Mathlib gen
*
*
*
* iorc.inc
*
      COMMON/IOLUNS/LIN,LOUT
      COMMON/GTSTAT/NTEST,NFAIL,IRC
      NTEST=NTEST+1
      IF(MODE.EQ.1) PRINT      1000, NTEST,CODE
                    WRITE(LOUT,1001) CODE
      RETURN
1000  FORMAT(' Test#',I3,' ( ',A,' ):     started')
1001  FORMAT('1',25X,30('*')/26X,'**   Testing Routine ',A,3X,'**'/
     +       26X,30('*'))
      END
