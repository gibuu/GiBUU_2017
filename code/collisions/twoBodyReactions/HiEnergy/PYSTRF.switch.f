      subroutine SetSwitchPYSTRF(inFlag)
      IMPLICIT NONE
      logical inFlag

      common /DataSwitchPYSTRF/ DoModif
      logical DoModif
      save /DataSwitchPYSTRF/

      data DoModif /.FALSE./

      DoModif = inFlag

      end



      SUBROUTINE PYSTRF(IP)
      IMPLICIT NONE
      integer IP

      common /DataSwitchPYSTRF/ DoModif
      logical DoModif
      save /DataSwitchPYSTRF/

      if (DoModif) then
         call PYSTRF_Modif(IP)
      else
         call PYSTRF_orig(IP)
      endif
      

      end
