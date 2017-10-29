      subroutine SetSwitchLUSTRF(inFlag)
      IMPLICIT NONE
      logical inFlag

      common /DataSwitchLUSTRF/ DoModif
      logical DoModif
      save /DataSwitchLUSTRF/

      data DoModif /.FALSE./

      DoModif = inFlag

      end



      SUBROUTINE LUSTRF(IP)
      IMPLICIT NONE
      integer IP

      common /DataSwitchLUSTRF/ DoModif
      logical DoModif
      save /DataSwitchLUSTRF/

      if (DoModif) then
         call LUSTRF_Modif(IP)
      else
         call LUSTRF_orig(IP)
      endif
      

      end
