c     version b - howard matis

c       options/extend_source
       integer function dls_acceptance(m,pt,y,wt,wt_err)
c
       implicit none
c
       integer im, ipt, iy
       integer ic
       integer n1, n2, n3
c
       logical l_fine
       logical l_init
c
       real acc_nc(4), dacc_nc(4), acc_tot, dacc_tot
       real dm, dpt, dy
       real m, pt, y
       real rm, rpt, ry
       real wt, wt_err, wt_n0
c
       include 'dls_filter.inc'
c
       save l_init
c
*
*--> explicit mass cut
*
       real mcut_low
       data mcut_low/0.050/
       data l_init /.false./
c
c begin execution; set defaults for returned variables:
c
       dls_acceptance = 0
       wt = -999.
       wt_err = 0.
       nic = 0
       do ic=1,4
          iacc_config(ic) = 0
          wt_config(ic) = -99999.
          wterr_config(ic) = 0.
       enddo
c
c initialization needed?
       if (.not.l_init) then
          call setup_filter
          l_init = .true.
       endif
c
c initialization call only?
       if (m.lt.0) return
       if (m.lt.mcut_low)return
c
c check outer boundaries:
       if ( m.lt.bot_b(1) .or.  m.gt.top_b(1) .or.
     &     pt.lt.bot_b(2) .or. pt.gt.top_b(2) .or.
     &      y.lt.bot_b(3) .or.  y.gt.top_b(3)) return
c
c use fine-grained tables, if we are not too close to their edges;
c it is assumed (but not required) that the lower edges in m and pt are
c the same for the fine and coarse tables.
       l_fine = ( m.ge.bot_a(1)          .and.
     &            m.le.top_a(1)-del_a(1) .and.
     &           pt.ge.bot_a(2)          .and.
     &           pt.le.top_a(2)-del_a(2) .and.
     &            y.ge.bot_a(3)+del_a(3) .and.
     &            y.le.top_a(3)-del_a(3) )
c
c find number of bin widths to first grid points:
       if (l_fine) then
          rm  =  (m-bot_a(1))/del_a(1) - 0.5
          rpt = (pt-bot_a(2))/del_a(2) - 0.5
          ry  =  (y-bot_a(3))/del_a(3) - 0.5
          n1 = nm_a
          n2 = npt_a
          n3 = ny_a
       else
          rm  =  (m-bot_b(1))/del_b(1) - 0.5
          rpt = (pt-bot_b(2))/del_b(2) - 0.5
          ry  =  (y-bot_b(3))/del_b(3) - 0.5
          n1 = nm_b
          n2 = npt_b
          n3 = ny_b
       endif
c
       im  = ifix( rm) + 1
       ipt = ifix(rpt) + 1
       iy  = ifix( ry) + 1
c
       if (l_fine) then
          im  = min( max(1,im),  nbin_a(1)-1)
          ipt = min( max(1,ipt), nbin_a(2)-1)
          iy  = min( max(1,iy),  nbin_a(3)-1)
       else
          im  = min( max(1,im),  nbin_b(1)-1)
          ipt = min( max(1,ipt), nbin_b(2)-1)
          iy  = min( max(1,iy),  nbin_b(3)-1)
       endif
c
       dm  =  rm - float(im-1)
       dpt = rpt - float(ipt-1)
       dy  =  ry - float(iy-1)
c
       do ic=1,4
          acc_nc(ic) = 0.
          dacc_nc(ic) = 0.
c
          if (l_fine) then
             call compute_accb(im,ipt,iy,ic, dm,dpt,dy, wt,wt_err,wt_n0,
     &                         n1,n2,n3, acc_a, err_acc_a, l_acc_a)
          else
             call compute_accb(im,ipt,iy,ic, dm,dpt,dy, wt,wt_err,wt_n0,
     &                         n1,n2,n3, acc_b, err_acc_b, l_acc_b)
          endif
c
c check limits, since interpolation with empty bins and extrapolation
c by a half bin can occur:
          if (wt.gt.weight_max .or. wt_n0.lt.n0_min) then ! apply limits
             wt     = 0.
             wt_err = 0.
             wt_n0  = 0.
          endif
c
c finite acceptance?
          if (wt.gt.0.) then
             dls_acceptance = 1
             nic = nic + 1
             iacc_config(ic) = 1
             acc_nc(ic) = 1./wt
             dacc_nc(ic) = wt_err*acc_nc(ic)**2
             wt_config(ic) = wt
             wterr_config(ic) = wt_err
          endif
       enddo
c
c return total acceptance and the unscaled individual acceptances:
       if (dls_acceptance.eq.1) then
          acc_tot  = acc_nc(1) + acc_nc(2) + acc_nc(3) + acc_nc(4)
          dacc_tot = sqrt(dacc_nc(1)**2 + dacc_nc(2)**2 +
     &                    dacc_nc(3)**2 + dacc_nc(4)**2)
          wt = 1./acc_tot
          wt_err = dacc_tot*wt**2
       endif

       return
       end

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

c       options/extend_source
       subroutine setup_filter
c
c modified: 17-jul-1995 by r.bossingham
c           support the geant acceptance (con) table format directly,
c           as opposed to the semi-digested format produced by cjn.
c           this form of the table is simply a list of the numbers of events
c           accepted, followed by a list of the numbers of events started.
c
c   in all discussions the pair ordering for fields,charges and configurations
c   will be:
c                   (right arm, left arm)
c
c     acc_tables   conf#       topology ( for charges=(+,+) )
c         mm         1         bend right  = r
c         mp         2         bend in     = i
c         pm         3         bend out    = o
c         pp         4         bend left   = l
c
c
c remember    1=right, 2=left !!!!!!!!!!
c
c the full matrix ! huan's abcd added by wkw
c
c
c                    iq
c                     1         2         3         4
c                |  (-,-)   | (+,+)   | (+,-)   | (-,+)   |
c         -------------------------------------------------
c   ib   1 (+,+) |  1(r) ab | 4(l) cd | 3(o) cb | 2(i) ad |
c        2 (+,-) |  2(i) aa | 3(o) cc | 4(l) ca | 1(r) ac |
c        3 (-,+) |  3(o) bb | 2(i) dd | 1(r) db | 4(l) bd |
c        4 (-,-) |  4(l) ba | 1(r) dc | 2(i) da | 3(o) bc |
c
c
c    iq*ib       conf #     bend    acc_tables
c    -----------------------------------------
c     (-,-)       1          r         mm
c     (-,+)       2          i         mp
c     (+,-)       3          o         pm
c     (+,+)       4          l         pp
c
c
c******************************************************************************
c
       implicit none
c
       include 'dls_filter.inc'
       include 'osbk_control.inc'
       include 'osbk_init.inc'          !needed only once per program
c
       character*132 acc_file_list_v93(1:4)
       character*132 acc_file_list_95a_fine(1:4)
       character*132 acc_file_list_95a(1:4)
       character*132 acc_file_list_95p_fine(1:4)
       character*132 acc_file_list_95p(1:4)
c
       integer i, j, k
       integer ic
c
c---- acc-filetypes = old ---- old table
c                   = pbp ---- the point by point
c                   = mm_ --- the configuration 1 table
c                   = mp_ --- the configuration 2 table
c                   = pm_ --- the configuration 3 table
c                   = pp_ --- the configuration 4 table
c                   = nul --- the readin acceptance is used!
c                   = 000 --- all are zeroed
c                   = con --- read in all 5 tables and do correction!
c
c chuck naudet's 1993 tables:
       data acc_file_list_v93 /
     &                  'dls_table_bend_right.acc',  !conf=1
     &                  'dls_table_bend_in.acc',     !conf=2
     &                  'dls_table_bend_out.acc',    !conf=3
     &                  'dls_table_bend_left.acc'/   !conf=4
c
c r.bossingham's 1995 tables with a 1300 mev/c momentum cut (for a-a collisions):
       data acc_file_list_95a_fine /   ! fine-grain tables
     &                  'table_bend_rgt_fine.dat',   !conf=1
     &                  'table_bend_in_fine.dat',    !conf=2
     &                  'table_bend_out_fine.dat',   !conf=3
     &                  'table_bend_lft_fine.dat'/   !conf=4
c
       data acc_file_list_95a /        ! coarse-grain tables
     &                  'table_bend_rgt.dat',        !conf=1
     &                  'table_bend_in.dat',         !conf=2
     &                  'table_bend_out.dat',        !conf=3
     &                  'table_bend_lft.dat'/        !conf=4
c
c r.bossingham's 1995 tables with a 2000 mev/c momentum cut (for p-p, p-d collisions):
       data acc_file_list_95p_fine /   ! fine-grain tables
     &                  'table_bend_rgt_fine.dat',   !conf=1
     &                  'table_bend_in_fine.dat',    !conf=2
     &                  'table_bend_out_fine.dat',   !conf=3
     &                  'table_bend_lft_fine.dat'/   !conf=4
c
       data acc_file_list_95p /        ! coarse-grain tables
     &                  'table_bend_rgt.dat',        !conf=1
     &                  'table_bend_in.dat',         !conf=2
     &                  'table_bend_out.dat',        !conf=3
     &                  'table_bend_lft.dat'/        !conf=4
c
c begin execution; set to read all tables:
       acc_filetype = 'con'
c
c sets min. number of events in bin for acceptance calculation
c in one configuration:
       acc_precision = 0.5
       n0_min = 1./acc_precision**2
c
c sets maximum weight in bin for one configuration:
c      weight_max = 2500.
       weight_max = 1000.
c
       if (acc_filetype.ne.'con' .and. acc_filetype.ne.'con') return
c
       do ic=1,4
c
c choose acceptance tables with the accf data card:
c   accf 1 = l_v93 = cjn format
c   accf 2 = l_95a = dls_sim format, 1300 mev/c cut
c   accf 3 = l_95p = dls_sim format, 2000 mev/c cut
c
          if (l_v93) then
             acc_file = '93_'//acc_file_list_v93(ic)
c
          elseif (l_95a) then
             acc_file = '95a_'//acc_file_list_95a(ic)
c
          elseif (l_95p) then
             acc_file = '95p_'//acc_file_list_95p(ic)
c
          else
             stop 'setup_filter: acceptance table format not specified'
          endif
c
          call init_acc         ! clear arrays, etc.
          call read_acceptance  ! read file
          call smooth1          ! smooth configuration
          call clean1           ! remove ragged edges and isolated points
c
c require that all tables of a given type have the same binning:
          if (ic.eq.1) then
             if (nbin1(1).gt.nm_b .or.
     &           nbin1(2).gt.npt_b .or.
     &           nbin1(3).gt.ny_b) then
                write (6,*) 'setup_filter: coarse input arrays too big'
                stop
             endif
c
             do i=1,3
                nbin_b(i) = nbin1(i)
                bot_b(i)  = bot1(i)
                top_b(i)  = top1(i)
                del_b(i)  = del1(i)
             enddo
          else
             do i=1,3
                if (nbin_b(i).ne.nbin1(i) .or.
     &              bot_b(i).ne.bot1(i) .or.
     &              top_b(i).ne.top1(i)) then
                   stop 'setup_filter error: coarse tables inconsistent'
                endif
             enddo
          endif
c
c fill coarse acceptance tables for four configurations:
          do i=1,nbin_b(1)
             do j=1,nbin_b(2)
                do k=1,nbin_b(3)
                   acc_b(i,j,k,ic)     = acc(i,j,k)
                   err_acc_b(i,j,k,ic) = err_acc(i,j,k)
                   l_acc_b(i,j,k,ic)   = l_acc(i,j,k)
                end do          ! k=1,nbin_b(3)
             end do             ! j=1,nbin_b(2)
          end do                ! i=1,nbin_b(1)
c
       end do                   ! ic=1,4
c
c record # configs at each m-pt-y cell:
       do i=1,nbin_b(1)
          do j=1,nbin_b(2)
             do k=1,nbin_b(3)
                configs_b(i,j,k) = 0.
                do ic=1,4
                   if (l_acc_b(i,j,k,ic))
     &                configs_b(i,j,k)=configs_b(i,j,k)+1.
                enddo           ! ic=1,4
             enddo              ! k=1,nbin_b(3)
          enddo                 ! j=1,nbin_b(2)
       enddo                    ! i=1,nbin_b(1)
c
c load the fine-grained tables for the 95a or 95p options only:
       if (.not.(l_95a.or.l_95p)) then
          do i=1,3
             nbin_a(i) = 0
             bot_a(i)  = 0.
             top_a(i)  = 0.
          enddo
          return
       endif
c
       do ic=1,4
          if (l_95a) then
             acc_file = '95a_'//acc_file_list_95a_fine(ic)
c
          elseif (l_95p) then
             acc_file = '95p_'//acc_file_list_95p_fine(ic)
c
          else
             stop ' setup_filter: internal logic error.'
          endif
c
          call init_acc         ! clear arrays, etc.
          call read_acceptance  ! read file
          call smooth1          ! smooth configuration
          call clean1           ! remove ragged edges and isolated points
c
c require that all tables of a given type have the same limits:
          if (ic.eq.1) then
             if (nbin1(1).gt.nm_a .or.
     &           nbin1(2).gt.npt_a .or.
     &           nbin1(3).gt.ny_a) then
                write (6,*) ' setup_filter: fine input arrays too big.'
                stop
             endif
c
             do i=1,3
                nbin_a(i) = nbin1(i)
                bot_a(i)  = bot1(i)
                top_a(i)  = top1(i)
                del_a(i)  = del1(i)
             enddo
          else
             do i=1,3
                if (nbin_a(i).ne.nbin1(i) .or.
     &              bot_a(i).ne.bot1(i) .or.
     &              top_a(i).ne.top1(i)) then
                   stop ' setup_filter error: fine tables inconsistent'
                endif
             enddo
          endif
c
c fill fine acceptance tables for four configurations:
          do i=1,nbin_a(1)
             do j=1,nbin_a(2)
                do k=1,nbin_a(3)
                   acc_a(i,j,k,ic)     = acc(i,j,k)
                   err_acc_a(i,j,k,ic) = err_acc(i,j,k)
                   l_acc_a(i,j,k,ic)   = l_acc(i,j,k)
                end do          ! k=1,nbin_a(3)
             end do             ! j=1,nbin_a(2)
          end do                ! i=1,nbin_a(1)
c
       end do                   ! ic=1,4
c
c record # configs at each m-pt-y cell:
       do i=1,nbin_a(1)
          do j=1,nbin_a(2)
             do k=1,nbin_a(3)
                configs_a(i,j,k) = 0.
                do ic=1,4
                   if (l_acc_a(i,j,k,ic))
     &                configs_a(i,j,k)=configs_a(i,j,k)+1.
                enddo           ! ic=1,4
             enddo              ! k=1,nbin_a(3)
          enddo                 ! j=1,nbin_a(2)
       enddo                    ! i=1,nbin_a(1)
c
       return
       end

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

c       options/extend_source
       subroutine init_acc
c
       implicit none
c
       include 'dls_filter.inc'
c
       integer i,k,j
c
c initialize acceptance and normalisation arrays
       do k=1,ny
          do j=1,npt
             do i=1,nm
                acc(i,j,k)    = 0.
                estart(i,j,k) = 0.
             enddo
          enddo
       enddo
c
       return
       end

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

c       options/extend_source
       subroutine read_acceptance
c
       use inputGeneral, only: path_To_Input
       implicit none
c
       include 'dls_filter.inc'
       include 'osbk_control.inc'

       character*132 ligne, file_name, junk1*45, junk2*6
c
       integer i,j,k
       integer n_and, n_gen, nmax, ngood_cell
c
       real acc_max, e, xp
c
c begin execution:
c
       file_name = trim(path_To_Input) // '/dls/' // acc_file
       if (acc_filetype.eq.'000') return
c
 200   continue
       open (unit=14, file=file_name, status='old', err=999)
c
       if (acc_filetype.ne.'pbp') then
          write (6,210) file_name
 210      format (/,' reading acceptance file: ',a50)
c
          read (14,'(a132)',end=999) ligne
          write (6,'(a132)') ligne
c
          if (l_v93) then
             read (14,211) (bot1(i),top1(i),nbin1(i), i=1,3)
 211         format (3(f6.3,f6.3,i4))
          else
             read (14,212) (bot1(i),top1(i),nbin1(i), i=1,3)
 212         format (3(f5.3,1x,f5.3,1x,i2,1x))
          endif
c
          do i=1,3
             del1(i) = (top1(i)-bot1(i)) / float(nbin1(i))
          end do
c
          ngood_cell = 0
          nmax = 0
          write (6,214)  bot1, top1, del1, nbin1
 214      format (27x, '      m       pt        y',/,
     &            ' acceptance, lower limits =',3f9.5,/,
     &            ' acceptance, upper limits =',3f9.5,/,
     &            ' acceptance, bin widths   =',3f9.5,/,
     &            ' acceptance, upper indices=',3i9,/)
c
          if (l_v93) then
             read (14,'(a132)') ligne
             write (6,'(a132)') ligne
          endif
c
          do i=1,nbin1(1)
             if (l_v93) read (14,'(a132)') ligne
             do j=1,nbin1(2)
                if (l_v93) then
                   read (14,*) (estart(i,j,k), k=1,nbin1(3))
                else
                   read (14,*) (acc(i,j,k),    k=1,nbin1(3))
                endif
             end do
          end do
c
          read (14,'(a132)') ligne
          write (6,'(a132)') ligne
c
          if (.not.l_v93) then
             read (14,212) (bot1(i),top1(i),nbin1(i), i=1,3)
          endif
c
          do i=1,nbin1(1)
             if (l_v93) read (14,'(a132)') ligne
             do j=1,nbin1(2)
                if (l_v93) then
                   read (14,*) (acc(i,j,k),    k=1,nbin1(3))
                else
                   read (14,*) (estart(i,j,k), k=1,nbin1(3))
                endif
                do k=1,nbin1(3)
                   nmax = nmax + 1
                   if (l_v93) then
                      acc(i,j,k) = acc(i,j,k)/1000.
                   else
                      if (estart(i,j,k).ne.0.) then
                         acc(i,j,k) = acc(i,j,k)/estart(i,j,k)
                      else
                         write (6,*) 'read_acceptance: bin not filled'
                         acc(i,j,k) = 0.
                      endif
                   endif
                   if (acc(i,j,k).ge.acc_max) acc_max=acc(i,j,k)
                   if (estart(i,j,k).gt.0) then
                      if (acc(i,j,k).gt.0) then
                         err_acc(i,j,k) = sqrt(acc(i,j,k)/estart(i,j,k))
                         if (acc(i,j,k).gt.1./weight_max) then
                            ngood_cell = ngood_cell + 1
                         endif
                      else
                         err_acc(i,j,k) = 1./estart(i,j,k)
                      end if
                   else
                      err_acc(i,j,k) = 0.
                   end if
c
c initially, set all bins to be accepted:
                   l_acc(i,j,k) = .true.
                end do
             end do
          end do
c
       else
c
c     pbp by wkw and cjn
          do i=1,1000000
             read (14,401,end=402) junk1,n_gen,junk2,n_and
             if (n_gen.gt.0) then
                e = float(n_and)/float(n_gen)
                w_pbp(i) = 1./e
c     i use poissons statistics. cjn
                if (e.gt.0. .and. e.lt.1)
     1               err_w_pbp(i) = sqrt((1-e)/(float(n_gen)*e**3))
             else
                w_pbp(i) = 0
                err_w_pbp(i) = 0.
             endif
          enddo
 401      format (a45,i7,a6,i4)
 402      continue
       end if
c
       write (6,*) ' the number of cubic acceptance cells = ',nmax
       write (6,*) ' the number of cubic filled cells     = ',ngood_cell
c
       xp = float(ngood_cell)*100./float(nmax)
       write(6,*) '                       percentage =',xp
c
       close(14)
       return
c
c error handling:
 999   continue
       write (*,*)' error accessing acceptance table:'
       write (*,'(3x,a80)') file_name
       stop
       end

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

c       options/extend_source
       subroutine smooth1
c
c smooth the most recently read acceptance file:
c
       implicit none

       include 'dls_filter.inc'
c
       integer ii
       integer j, jj, jjj
       integer jm1, jp1
       integer k, kk, kkk
       integer km1, kp1
       integer l, ll, lll
       integer lm1, lp1
       integer nsum_j(-1:1)
       integer nsum_k(-1:1)
       integer nsum_l(-1:1)
c
       real an
       real asum
       real asum_j(-1:1)
       real asum_k(-1:1)
       real asum_l(-1:1)
       real bsum
c
c begin execution;
c turn off averaging at array limits:
c
       do j=1,nbin1(1)
          if (j.eq.1 .or. j.eq.nbin1(1)) then
             jm1 = j
             jp1 = j
          else
             jm1 = j - 1
             jp1 = j + 1
          endif
          do k=1,nbin1(2)
             if (k.eq.1 .or. k.eq.nbin1(2)) then
                km1 = k
                kp1 = k
             else
                km1 = k - 1
                kp1 = k + 1
             endif
             do l=1,nbin1(3)
                if (l.eq.1 .or. l.eq.nbin1(3)) then
                   lm1 = l
                   lp1 = l
                else
                   lm1 = l - 1
                   lp1 = l + 1
                endif
c
c initial summing variables:
                do ii=-1,1
                   asum_j(ii) = 0.
                   asum_k(ii) = 0.
                   asum_l(ii) = 0.
c
                   nsum_j(ii) = 0
                   nsum_k(ii) = 0
                   nsum_l(ii) = 0
                enddo
c
c sum events and filled bins in selected bin and its neighbors:
                jjj = jm1 - j
                do jj=jm1,jp1
                   kkk = km1 - k
                   do kk=km1,kp1
                      lll = lm1 - l
                      do ll=lm1,lp1
c
c "acceptance" = events/estart -> "an" = number of events = acceptance*estart:
                         an = acc(jj,kk,ll)*estart(jj,kk,ll)
c
                         if (an.gt.0.) then
                            asum_j(jjj) = asum_j(jjj) + an
                            asum_k(kkk) = asum_k(kkk) + an
                            asum_l(lll) = asum_l(lll) + an
                            nsum_j(jjj) = nsum_j(jjj) + 1
                            nsum_k(kkk) = nsum_k(kkk) + 1
                            nsum_l(lll) = nsum_l(lll) + 1
                         endif
                         lll = lll + 1
                      enddo
                      kkk = kkk + 1
                   enddo
                   jjj = jjj + 1
                enddo
c
c to reduce problems at the edge of the acceptance, we only average
c in directions in which there are some filled elements on both sides:
c
                if (asum_j(-1).eq.0. .or. asum_j(+1).eq.0.) then
                   jm1 = j
                   jp1 = j
                endif
c
                if (asum_k(-1).eq.0. .or. asum_k(+1).eq.0.) then
                   km1 = k
                   kp1 = k
                endif
c
                if (asum_l(-1).eq.0. .or. asum_l(+1).eq.0.) then
                   lm1 = l
                   lp1 = l
                endif
c
c try to avoid interpolating across an edge:
c                if (nsum_j(0).lt.3) then
c                   km1 = k
c                   kp1 = k
c                   lm1 = l
c                   lp1 = l
c                endif
c
c                if (nsum_k(0).lt.3) then
c                   jm1 = j
c                   jp1 = j
c                   lm1 = l
c                   lp1 = l
c                endif
c
c                if (nsum_l(0).lt.3) then
c                   jm1 = j
c                   jp1 = j
c                   km1 = k
c                   kp1 = k
c                endif
c
                asum = 0.
                bsum = 0.
                do jj=jm1,jp1
                   do kk=km1,kp1
                      do ll=lm1,lp1
                         asum = asum + acc(jj,kk,ll)*estart(jj,kk,ll)
                         bsum = bsum + estart(jj,kk,ll)
                      enddo
                   enddo
                enddo
c
c store in temporary arrays to avoid disrupting the originals:
                atemp(j,k,l) = asum/bsum
                if (asum.gt.0.) then
                   etemp(j,k,l) = sqrt(asum)/bsum
                else
                   etemp(j,k,l) = 1./bsum
                endif
c
             end do
          end do
       end do
c
c reload the original arrays with the new values:
       do j=1,nbin1(1)
          do k=1,nbin1(2)
             do l=1,nbin1(3)
                acc(j,k,l)     = atemp(j,k,l)
                err_acc(j,k,l) = etemp(j,k,l)
             end do
          end do
       end do
c
       return
       end

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

c       options/extend_source
       subroutine clean1
c
c remove low-accuracy and isolated points from the most recently read
c acceptance table.
c
       implicit none
c
       include 'dls_filter.inc'
c
       integer i
       integer j, jj, jm1, jp1
       integer k, kk, km1, kp1
       integer l, ll, lm1, lp1
       integer msum, nsum
c
       real asum
       real ave_acc
       real rel_err
c
c begin execution;
c turn off m-pt edge and isolated points w/ low accuracy or high weight:
       do j=1,nbin1(1)
          jm1 = max(1,       j-1)
          jp1 = min(nbin1(1),j+1)
          do k=1,nbin1(2)
             km1 = max(1,       k-1)
             kp1 = min(nbin1(2),k+1)
             do l=1,nbin1(3)
                l_acc(j,k,l) = .false.
                if (acc(j,k,l).gt.0.) then
                   msum = 0
                   nsum = 0
                   do jj=jm1,jp1
                      do kk=km1,kp1
                         msum = msum + 1
c
                         if (acc(jj,kk,l).gt.0.) then
                            rel_err = err_acc(jj,kk,l)/acc(jj,kk,l)
                            if (1./rel_err**2.ge.n0_min
     &                          .and. 1./acc(jj,kk,l).le.weight_max)
     &                         nsum = nsum + 1
                         endif
c
                      enddo
                   enddo
c
                   if (float(nsum)/float(msum).gt.0.66) then
                      l_acc(j,k,l) = .true.
                   else
                      rel_err = err_acc(j,k,l)/acc(j,k,l)
                      l_acc(j,k,l) = 1./rel_err**2.ge.n0_min
     &                     .and. 1./acc(j,k,l).le.weight_max
                   endif
                endif
             enddo
          enddo
       enddo
c
c iterate turn off of m-pt edge/isolated points w/ low accuracy or high
c (or infinite) weight:
       do i=1,3
          do j=1,nbin1(1)
             jm1 = max(1,       j-1)
             jp1 = min(nbin1(1),j+1)
             do k=1,nbin1(2)
                km1 = max(1,       k-1)
                kp1 = min(nbin1(2),k+1)
                do l=1,nbin1(3)
                   ltemp(j,k,l) = l_acc(j,k,l)
                   if (l_acc(j,k,l)) then
                      msum = 0
                      nsum = 0
                      do jj=jm1,jp1
                         do kk=km1,kp1
                            msum = msum + 1
                            if (l_acc(jj,kk,l)) nsum=nsum+1
                         enddo
                      enddo
                      if (float(nsum)/float(msum).lt.0.66) then
                         rel_err = err_acc(j,k,l)/acc(j,k,l)
                         ltemp(j,k,l) = 1./rel_err**2.ge.n0_min
     &                        .and. 1./acc(j,k,l).le.weight_max
                      endif
c
c discard isolated points in m-pt:
                      if (msum.lt.6) then
                         if (nsum.lt.2) ltemp(j,k,l)=.false.
                      elseif (msum.lt.9) then
                         if (nsum.lt.3) ltemp(j,k,l)=.false.
                      else
                         if (nsum.lt.5) ltemp(j,k,l)=.false.
                      endif
                   endif
                enddo
             enddo
          enddo
c
c save the temporary array:
          do j=1,nbin1(1)
             do k=1,nbin1(2)
                do l=1,nbin1(3)
                   l_acc(j,k,l) = ltemp(j,k,l)
                end do
             end do
          end do
c
       enddo   ! do i=1,3
c
c eliminate edge points where the weight is rising rapidly:
       do j=1,nbin1(1)
          jm1 = max(1,       j-1)
          jp1 = min(nbin1(1),j+1)
          do k=1,nbin1(2)
             km1 = max(1,       k-1)
             kp1 = min(nbin1(2),k+1)
             do l=1,nbin1(3)
                lm1 = max(1,       l-1)
                lp1 = min(nbin1(3),l+1)
c
                if (l_acc(j,k,l)) then
                   asum = 0.
                   msum = 0
                   nsum = 0
                   do jj=jm1,jp1
                      do kk=km1,kp1
                         do ll=lm1,lp1
                            msum = msum + 1
                            asum = asum + acc(jj,kk,ll)
                            if (l_acc(jj,kk,ll)) then
                               nsum = nsum + 1
                            endif
                         enddo
                      enddo
                   enddo
c
                   if (float(nsum)/float(msum).lt.0.67) then
                      ave_acc = asum/float(msum)
                      ltemp(j,k,l) =  (acc(j,k,l).gt.0.5*ave_acc)
                   endif
c
c discard isolated points:
                   if (msum.lt.18) then
                      ltemp(j,k,l) = (nsum.ge.3) ! corner or edge
                   elseif (msum.lt.27) then
                      ltemp(j,k,l) = (nsum.ge.4) ! face
                   else
                      ltemp(j,k,l) = (nsum.ge.7) ! bulk
                   endif
                endif
             end do
          end do
       end do
c
c record the temporary copy:
       do j=1,nbin1(1)
          do k=1,nbin1(2)
             do l=1,nbin1(3)
                l_acc(j,k,l) = ltemp(j,k,l)
             end do
          end do
       end do
c
c eliminate isolated points (we require 3-7 bins out of 8-27);
c we iterate, in case large changes occur during the first pass:
       do i=1,3
          do j=1,nbin1(1)
             jm1 = max(1,       j-1)
             jp1 = min(nbin1(1),j+1)
             do k=1,nbin1(2)
                km1 = max(1,       k-1)
                kp1 = min(nbin1(2),k+1)
                do l=1,nbin1(3)
                   lm1 = max(1,       l-1)
                   lp1 = min(nbin1(3),l+1)
c
                   if (l_acc(j,k,l)) then
                      msum= 0
                      nsum = 0
                      do jj=jm1,jp1
                         do kk=km1,kp1
                            do ll=lm1,lp1
                               msum = msum + 1
                               if (l_acc(jj,kk,ll)) nsum=nsum+1
                            enddo
                         enddo
                      enddo
c
                      if (msum.lt.18) then
                         ltemp(j,k,l) = (nsum.ge.3)     ! corner or edge
                      elseif (msum.lt.27) then
                         ltemp(j,k,l) = (nsum.ge.4)     ! face
                      else
                         ltemp(j,k,l) = (nsum.ge.7)     ! bulk
                      endif
                   endif
c
                   if (ltemp(j,k,l)) then
                      msum= 0
                      nsum = 0
                      do jj=jm1,jp1
                         do kk=km1,kp1
                            msum = msum + 1
                            if (l_acc(jj,kk,l)) nsum=nsum+1
                         enddo
                      enddo
                      if (msum.lt.6) then
                         ltemp(j,k,l) = (nsum.ge.2) ! 2x2, 2x3
                      elseif (msum.lt.9) then
                         ltemp(j,k,l) = (nsum.ge.3) ! 2x2, 2x3
                      else
                         ltemp(j,k,l) = (nsum.ge.5) ! 3x3
                      endif
                   endif
c
                end do
             end do
          end do
c
c reload the original arrays:
          do j=1,nbin1(1)
             do k=1,nbin1(2)
                do l=1,nbin1(3)
                   l_acc(j,k,l) = ltemp(j,k,l)
                end do
             end do
          end do
c
       enddo
c
c trim the the edges again:
c
c turn off m-pt edge/isolated points w/ low accuracy or high weight:
       do i=1,3
          do j=1,nbin1(1)
             jm1 = max(1,       j-1)
             jp1 = min(nbin1(1),j+1)
             do k=1,nbin1(2)
                km1 = max(1,       k-1)
                kp1 = min(nbin1(2),k+1)
                do l=1,nbin1(3)
                   ltemp(j,k,l) = l_acc(j,k,l)
                   if (l_acc(j,k,l)) then
                      msum = 0
                      nsum = 0
                      do jj=jm1,jp1
                         do kk=km1,kp1
                            msum = msum + 1
                            if (l_acc(jj,kk,l)) nsum=nsum+1
                         enddo
                      enddo
c
c discard low-weight and high-error points at edges in m-pt:
                      if (float(nsum)/float(msum).lt.0.66) then
                         rel_err = err_acc(j,k,l)/acc(j,k,l)
                         ltemp(j,k,l) = 1./rel_err**2.ge.n0_min
     &                        .and. 1./acc(j,k,l).le.weight_max
                      endif
c
c discard isolated points in m-pt:
                      if (msum.lt.6) then
                         if (nsum.lt.2) ltemp(j,k,l)=.false.
                      elseif (msum.lt.9) then
                         if (nsum.lt.3) ltemp(j,k,l)=.false.
                      else
                         if (nsum.lt.5) ltemp(j,k,l)=.false.
                      endif
                   endif
                enddo
             enddo
          enddo
c
c save the temporary array:
          do j=1,nbin1(1)
             do k=1,nbin1(2)
                do l=1,nbin1(3)
                   l_acc(j,k,l) = ltemp(j,k,l)
                end do
             end do
          end do
       enddo
c
c eliminate edge points where the weight is rising rapidly:
       do j=1,nbin1(1)
          jm1 = max(1,       j-1)
          jp1 = min(nbin1(1),j+1)
          do k=1,nbin1(2)
             km1 = max(1,       k-1)
             kp1 = min(nbin1(2),k+1)
             do l=1,nbin1(3)
                lm1 = max(1,       l-1)
                lp1 = min(nbin1(3),l+1)
c
                if (l_acc(j,k,l)) then
                   asum = 0.
                   msum = 0
                   nsum = 0
                   do jj=jm1,jp1
                      do kk=km1,kp1
                         do ll=lm1,lp1
                            msum = msum + 1
                            asum = asum + acc(jj,kk,ll)
                            if (l_acc(jj,kk,ll)) then
                               nsum = nsum + 1
                            endif
                         enddo
                      enddo
                   enddo
c
                   if (float(nsum)/float(msum).lt.0.66) then
                      ave_acc = asum/float(msum)
                      ltemp(j,k,l) =  (acc(j,k,l).gt.0.5*ave_acc)
                   endif
c
c discard isolated points in m-pt:
                   if (msum.lt.6) then
                      if (nsum.lt.2) ltemp(j,k,l)=.false.
                   elseif (msum.lt.9) then
                      if (nsum.lt.3) ltemp(j,k,l)=.false.
                   else
                      if (nsum.lt.5) ltemp(j,k,l)=.false.
                   endif
                endif
             end do
          end do
       end do
c
c save the temporary array:
       do j=1,nbin1(1)
          do k=1,nbin1(2)
             do l=1,nbin1(3)
                l_acc(j,k,l) = ltemp(j,k,l)
             end do
          end do
       end do
c
c make a final pass to eliminate isolated points only:
       do i=1,3
          do j=1,nbin1(1)
             jm1 = max(1,       j-1)
             jp1 = min(nbin1(1),j+1)
             do k=1,nbin1(2)
                km1 = max(1,       k-1)
                kp1 = min(nbin1(2),k+1)
                do l=1,nbin1(3)
                   lm1 = max(1,       l-1)
                   lp1 = min(nbin1(3),l+1)
c
                   if (l_acc(j,k,l)) then
                      asum = 0.
                      msum = 0
                      nsum = 0
                      do jj=jm1,jp1
                         do kk=km1,kp1
                            msum = msum + 1
                            if (l_acc(jj,kk,l)
     &                          .and. acc(jj,kk,l).gt.0.) then
                               rel_err = err_acc(jj,kk,l)/acc(jj,kk,l)
                               if (1./rel_err**2.ge.n0_min
     &                              .and. 1./acc(jj,kk,l).le.weight_max)
     &                            nsum = nsum + 1
                            endif
                         enddo
                      enddo
c
                      if (msum.lt.6) then
                         if (nsum.lt.2) ltemp(j,k,l)=.false.
                      elseif (msum.lt.9) then
                         if (nsum.lt.3) ltemp(j,k,l)=.false.
                      else
                         if (nsum.lt.5) ltemp(j,k,l)=.false.
                      endif
                   endif
                end do
             end do
          end do
c
c save the temporary array:
          do j=1,nbin1(1)
             do k=1,nbin1(2)
                do l=1,nbin1(3)
                   l_acc(j,k,l) = ltemp(j,k,l)
                end do
             end do
          end do
       enddo
c
       return
       end

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

c       options/extend_source
       subroutine compute_accb(im,ipt,iy,ic, dm,dpt,dy, wt,wt_err,wt_n0,
     &                         n1,n2,n3, acc_i, err_acc_i, l_acc_i)
c
       implicit none
c
       integer ic
       integer iim, iipt, iiy
       integer im, ipt, iy
       integer im1, ipt1, iy1
       integer n1, n2, n3
c
       logical l_acc_i(n1,n2,n3,4)
c
       real a
       real a1, a2, a3, a4, b1, b2
       real acc_i(n1,n2,n3,4)
       real dm, dpt, dy
       real err_acc_i(n1,n2,n3,4)
       real ra
       real ra1, ra2, ra3, ra4, rb1, rb2
       real wt, wt_err, wt_n0
c
c begin execution; set defaults:
c
       a      = 0.
       ra     = 0.
       wt     = 0.
       wt_n0  = 0.
       wt_err = 99999999.
c
c check the sanity of array indices:
       if (im .lt.1 .or. im .ge.n1 .or.
     &     ipt.lt.1 .or. ipt.ge.n2 .or.
     &     iy .lt.1 .or. iy .ge.n3 .or.
     &     ic .lt.1 .or. ic .gt.4) return
c
c check that extrapolation does not exceed one half bin beyond the
c selected interpolation cube:
       if (abs( dm-0.5).gt.1. .or.
     &     abs(dpt-0.5).gt.1. .or.
     &     abs( dy-0.5).gt.1.) return
c
       im1  = im  + 1
       ipt1 = ipt + 1
       iy1  = iy  + 1
c
c find the bin in which our point sits; begin by assuming the lowest bins:
       iim  = im
       iipt = ipt
       iiy  = iy
       if (dm .gt.0.5) iim =im1
       if (dpt.gt.0.5) iipt=ipt1
       if (dy .gt.0.5) iiy =iy1
c
c verify that the point is within the defined acceptance:
       if (.not.l_acc_i(iim,iipt,iiy,ic)) return
c
c interpolate in mass:
       a1 = acc_i(im,ipt ,iy, ic)*(1.-dm) + acc_i(im1,ipt, iy, ic)*dm
       a2 = acc_i(im,ipt1,iy, ic)*(1.-dm) + acc_i(im1,ipt1,iy, ic)*dm
       a3 = acc_i(im,ipt ,iy1,ic)*(1.-dm) + acc_i(im1,ipt, iy1,ic)*dm
       a4 = acc_i(im,ipt1,iy1,ic)*(1.-dm) + acc_i(im1,ipt1,iy1,ic)*dm
c
c we have already smoothed using adjacent points, so we treat (as a reasonable
c approximation) the errors as correlated and do not add them in quadrature:
       ra1 = err_acc_i(im,ipt, iy, ic)*(1.-dm)
     &     + err_acc_i(im1,ipt, iy, ic)*dm
       ra2 = err_acc_i(im,ipt1,iy, ic)*(1.-dm)
     &     + err_acc_i(im1,ipt1,iy, ic)*dm
       ra3 = err_acc_i(im,ipt, iy1,ic)*(1.-dm)
     &     + err_acc_i(im1,ipt, iy1,ic)*dm
       ra4 = err_acc_i(im,ipt1,iy1,ic)*(1.-dm)
     &     + err_acc_i(im1,ipt1,iy1,ic)*dm
c
c since a slight extrapolation is allowed, guarantee that values are positive:
        a1 = max(a1,0.)
        a2 = max(a2,0.)
        a3 = max(a3,0.)
        a4 = max(a4,0.)
       ra1 = max(ra1,0.)
       ra2 = max(ra2,0.)
       ra3 = max(ra3,0.)
       ra4 = max(ra4,0.)
c
c interpolate in pt:
        b1 =  a1*(1.-dpt) +  a2*dpt
        b2 =  a3*(1.-dpt) +  a4*dpt
       rb1 = ra1*(1.-dpt) + ra2*dpt
       rb2 = ra3*(1.-dpt) + ra4*dpt
c
c guarantee that values are positive:
        b1 = max(b1,0.)
        b2 = max(b2,0.)
       rb1 = max(rb1,0.)
       rb2 = max(rb2,0.)
c
c interpolate in rapidity:
        a =  b1*(1.-dy) +  b2*dy
       ra = rb1*(1.-dy) + rb2*dy
c
c guarantee that values are positive:
        a = max(a,0.)
       ra = max(ra,0.)
c
c convert to weights:
       if (a.gt.0. .and. ra.gt.0. .and.
     &     a*a.gt.ra*1.e-38 .and. ra.gt.a*1.e-19) then
          wt     = 1./a
          wt_err = (ra/a)/a
          wt_n0  = (a/ra)**2
       else
          wt     = 0.
          wt_err = 0.
          wt_n0  = 0.
       endif
c
       return
       end
