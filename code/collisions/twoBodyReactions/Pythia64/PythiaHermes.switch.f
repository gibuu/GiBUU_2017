      subroutine SetSwitchPythiaHermes(inFlag)
      IMPLICIT NONE
      logical inFlag

c$$$      common /DataSwitchPythiaHermes/ UseHermes
c$$$      logical UseHermes
c$$$      save /DataSwitchPythiaHermes/
c$$$
c$$$      data UseHermes /.FALSE./
c$$$
c$$$      UseHermes = inFlag

      if (inFlag) then
         write(*,*) 'Error: We do not allow for ',
     $        'PYTHIA6.4 with HERMES tuning.'
         write(*,*) '(HERMES tuned v6.2 and not v6.4)'
         write(*,*) 'STOP!'
         stop
      endif
      end


      logical function GetSwitchPythiaHermes()
      IMPLICIT NONE

c$$$      common /DataSwitchPythiaHermes/ UseHermes
c$$$      logical UseHermes
c$$$      save /DataSwitchPythiaHermes/
c$$$
c$$$      write(*,*) '####### GetSwitchPythiaHermes:',UseHermes
c$$$
c$$$      GetSwitchPythiaHermes = UseHermes

      GetSwitchPythiaHermes = .false.
      return
      end


c$$$      SUBROUTINE PYDIFF
c$$$      IMPLICIT NONE
c$$$
c$$$      common /DataSwitchPythiaHermes/ UseHermes
c$$$      logical UseHermes
c$$$      save /DataSwitchPythiaHermes/
c$$$
c$$$      if (UseHermes) then
c$$$         call PYDIFF_hermes
c$$$      else
c$$$         call PYDIFF_orig
c$$$      endif
c$$$      end
c$$$
c$$$
c$$$      SUBROUTINE PYGAGA(IGAGA,WTGAGA)
c$$$      IMPLICIT NONE
c$$$      integer iGAGA
c$$$      double precision WTGAGA
c$$$
c$$$      common /DataSwitchPythiaHermes/ UseHermes
c$$$      logical UseHermes
c$$$      save /DataSwitchPythiaHermes/
c$$$
c$$$      if (UseHermes) then
c$$$         call PYGAGA_hermes(IGAGA,WTGAGA)
c$$$      else
c$$$         call PYGAGA_modif(IGAGA,WTGAGA) ! modif !
c$$$      endif
c$$$      end
c$$$
c$$$
c$$$      SUBROUTINE PYSIGH(NCHN,SIGS)
c$$$      IMPLICIT NONE
c$$$      integer NCHN
c$$$      double precision SIGS
c$$$
c$$$      common /DataSwitchPythiaHermes/ UseHermes
c$$$      logical UseHermes
c$$$      save /DataSwitchPythiaHermes/
c$$$
c$$$      if (UseHermes) then
c$$$         call PYSIGH_hermes(NCHN,SIGS)
c$$$      else
c$$$         call PYSIGH_modif(NCHN,SIGS)
c$$$      endif
c$$$      end
c$$$
c$$$
c$$$      SUBROUTINE PYXTOT
c$$$      IMPLICIT NONE
c$$$
c$$$      common /DataSwitchPythiaHermes/ UseHermes
c$$$      logical UseHermes
c$$$      save /DataSwitchPythiaHermes/
c$$$
c$$$      if (UseHermes) then
c$$$         call PYXTOT_hermes
c$$$      else
c$$$         call PYXTOT_orig
c$$$      endif
c$$$      end
c$$$
