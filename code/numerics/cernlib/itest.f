***************************************************************************************************************************
* This subroutine is part of CERNLIB
* (C) Copyright CERN except where explicitly stated otherwise. 
* Permission to use and/or redistribute this work is granted under the terms of the GNU General Public License. 
***************************************************************************************************************************
*
* $Id: itest.F,v 1.1.1.1 1996/04/01 15:01:12 mclareni Exp $
*
* $Log: itest.F,v $
* Revision 1.1.1.1  1996/04/01 15:01:12  mclareni
* Mathlib gen
*
*
      FUNCTION ITEST(CODE,TEST)
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
      CHARACTER*(*) CODE
      LOGICAL TEST
      IF(TEST) THEN
        PRINT 1000,NTEST,CODE
        IF (LOUT .NE. 6) WRITE(LOUT,1000) NTEST,CODE
        ITEST=0
      ELSE
        PRINT 1001,NTEST,CODE
        IF (LOUT .NE. 6) WRITE(LOUT,1001) NTEST,CODE
        ITEST=1
      ENDIF
      NFAIL=NFAIL+ITEST
1000  FORMAT(' Test#',I3,' ( ',A,' ):     completed successfully')
1001  FORMAT(' Test#',I3,' ( ',A,' ): *** failed ***')
      RETURN
      END
