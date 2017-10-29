*
*       OSBK_INIT.INC
*
*       The default information is included here in data statements
*       It should only be included once per program
*
        data dsf_file/'  '/
        data ntup_file/'  '/
        data stat_file/'statistics.dat'/


        data m_min_dist, pt_min_dist, y_min_dist / 0.025, 0., 0.4/
        data m_max_dist, pt_max_dist, y_max_dist / 1.4, 1.2, 2./
        data im_max, ipt_max, iy_max / 55, 48, 32/
        data dmass,dpt,dy /0.025, 0.025, 0.05/


        data ndiv_sam/8/
        data read_scaler/.true./
        data mult_cut/.false./
        data apply_wt/.false./
        data proj_bends/.false./
        data wtcut/1000./

*
*--> mapping of pair type to bend configuration
*
        data ibcf_os_pcf/2,3,1,4,4,1,3,2/   !maps OS pair conf to acc. conf
        data ibcf_ss_pcf/4,1,3,2,2,3,1,4/   !maps SS pair conf to acc. conf


*
*--> Acceptance table format
*
*
        Data l_V93 /.false./
        Data l_95A /.false./
        Data l_95P /.false./
*
* ----------------- End of osbk_init.inc --------------------
