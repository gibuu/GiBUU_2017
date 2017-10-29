c DLS_FILTER.INC
c
c Input file information:
      Character*132 acc_file
      Character*3 acc_filetype
      Common /FILE/ acc_filetype, acc_file
c
c Cuts applied to acceptance:
      Real acc_precision
      Real n0_min
      Real weight_max
      Common /REQURD/ acc_precision, n0_min, weight_max
c
c Configuration-related information returned:
      Integer iacc_config(4)
      Integer nic
      Real wt_config(4), wterr_config(4)
      Common /ICPASS/ iacc_config, nic, wt_config, wterr_config
c
c Input and manipulation arrays:
      Integer nm, npt, ny
      Parameter (nm=55, npt=48, ny=49)
c
      Integer nbin1(3)
      Logical l_acc(nm,npt,ny)
      Real acc(nm,npt,ny)
      Real bot1(3), top1(3), del1(3)
      Real err_acc(nm,npt,ny)
      Real estart(nm,npt,ny)
      Common /ACC_INPUT/ acc, err_acc, l_acc, estart,
     &                   nbin1, bot1, top1, del1
c
      Logical ltemp
      Real atemp, etemp
      Common /TSPACE/ atemp(nm,npt,ny), etemp(nm,npt,ny),
     &                ltemp(nm,npt,ny)
c
c Fine-grained acceptance probability arrays:
      Integer nm_a, npt_a, ny_a
      Parameter (nm_a=18, npt_a=20, ny_a=49)
c
      Integer nbin_a(3)
      Logical l_acc_a(nm_a,npt_a,ny_a,4)
      Real acc_a(nm_a,npt_a,ny_a,4)
      Real bot_a(3), top_a(3), del_a(3)
      Real configs_a(nm_a,npt_a,ny_a)
      Real err_acc_a(nm_a,npt_a,ny_a,4)
      Common /ACC_FINE/ acc_a, err_acc_a, l_acc_a, configs_a,
     &                  nbin_a, bot_a, top_a, del_a
c
c Coarse-grained acceptance probability arrays:
      Integer nm_b, npt_b, ny_b
      Parameter (nm_b=55, npt_b=48, ny_b=32)
c
      Integer nbin_b(3)
      Logical l_acc_b(nm_b,npt_b,ny_b,4)
      Real acc_b(nm_b,npt_b,ny_b,4)
      Real bot_b(3), top_b(3), del_b(3)
      Real configs_b(nm_b,npt_b,ny_b)
      Real err_acc_b(nm_b,npt_b,ny_b,4)
      Common /ACC_COARSE/ acc_b, err_acc_b, l_acc_b, configs_b,
     &                    nbin_b, bot_b, top_b, del_b
c
      Integer n_evt_max
      Parameter (n_evt_max=200000)
      Real w_pbp(n_evt_max), err_w_pbp(n_evt_max)
      Common /WEIGHT_COM/ w_pbp, err_w_pbp
c
c End of DLS_FILTER.INC
