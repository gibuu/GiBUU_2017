c       Options/extend_source
       Logical Function beta_filter(dataset_name,xm,pt,y,weight_cut)
c
c Written:  18-NOV-1996 by R.Bossingham
c
c Purpose:  Provide a trivial interface to the DLS filter
c
c
       Implicit none
c
       Character*3 dataset_name
       Character*3 dataset_name_old
c
       Logical l_first
       Logical CHECK_POINT
c
       Real weight_cut
       Real xm, pt, y
c
       Include 'osbk_control.inc'
c
       Data dataset_name_old /'   '/
       Data l_first /.true./
c
c Begin execution:
c
       if (l_first .or. dataset_name .ne. dataset_name_old) then
          l_first = .false.
          l_V93 = .false.
          l_95A = .false.
          l_95P = .false.
c
          if (dataset_name.eq.'95A') then
             l_95A = .true.
          elseif (dataset_name.eq.'95P') then
             l_95P = .true.
          else
             write (6,700) dataset_name
 700         format (//,' *** Fatal error ***',/,
     &                  ' Illegal dataset name; use 95A or 95P',/,
     &                  ' Terminating.',//)
             stop
          endif
c
          dataset_name_old = dataset_name
c
          call TABLE_INITIALIZATION ! Initialize acceptance arrays
       endif
c
       if (weight_cut.le.0.) then
          write (6,710) weight_cut
 710      format (//,' *** Warning ***',/,
     &               ' Weight cut =',E14.7,//)
          beta_filter  = .FALSE.
          return
       endif
c
       beta_filter = CHECK_POINT(xm,pt,y,weight_cut)
       return
       end
c
c       Options/extend_source
       Subroutine TABLE_INITIALIZATION
c
c Initialize acceptance arrays.
c
       Implicit none
c
       Integer DLS_ACCEPTANCE
       Integer iacc
       Integer ic
       Integer npair_conf(0:4), npair_acc(0:4)
c
       Real pt,y
       Real wt,wterr
c
c Begin execution:
c
       do ic=0,4
          npair_conf(ic) = 0
          npair_acc(ic) = 0
       enddo
c
c Initialization call, only:
       iacc = DLS_ACCEPTANCE(-1.0,pt,y,wt,wterr) ! for initialization
c
       return
       end
c
c       Options/extend_source
       Logical Function CHECK_POINT(xm,pt,y,weight_cut)
c
       Implicit none
c
       Integer DLS_ACCEPTANCE
       Integer iacc
       Integer ic
c
       Real weight_cut
       Real wt,wterr
       Real xm,pt,y
c
       Include 'dls_filter.inc'
c
c Begin execution:
c
       iacc = DLS_ACCEPTANCE(xm,pt,y,wt,wterr)
c
       CHECK_POINT = .FALSE.
c
       if (iacc.eq.0) return
c
       do ic=1,4
          if (wt_config(ic).lt.weight_cut) then
             CHECK_POINT = .TRUE.
             return
          endif
       enddo
c
       return
       end
