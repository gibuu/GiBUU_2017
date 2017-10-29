***************************************************************************************************************************
* This subroutine is part of CERNLIB
* (C) Copyright CERN except where explicitly stated otherwise. 
* Permission to use and/or redistribute this work is granted under the terms of the GNU General Public License. 
***************************************************************************************************************************
*
* $Id: pagend.F,v 1.1.1.1 1996/04/01 15:01:12 mclareni Exp $
*
* $Log: pagend.F,v $
* Revision 1.1.1.1  1996/04/01 15:01:12  mclareni
* Mathlib gen
*
*
      SUBROUTINE PAGEND(CODE)
C     This subroutine prints a page end message
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
                    WRITE(LOUT,1001) CODE
      RETURN
1001  FORMAT(/26X,30('*')/26X,'**   End of Test of  ',A,3X,'**'/
     +       26X,30('*'))
      END
