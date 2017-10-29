*
*       OSBK_CONTROL.INC
*
*       file to include the control information for running the package

        character*76 dsf_file,ntup_file,stat_file

        common/files/dsf_file,ntup_file,stat_file

        real*4 m_min_dist,pt_min_dist,y_min_dist
        real*4 m_max_dist,pt_max_dist,y_max_dist
        real*4 dmass,dpt,dy
        integer*4 im_max,ipt_max,iy_max
        integer*4 im1,im2,ipt1,ipt2,iy1,iy2

        common/cut_limits/m_min_dist,pt_min_dist,y_min_dist,
     &                    m_max_dist,pt_max_dist,y_max_dist,
     &                    im_max,ipt_max,iy_max,
     &                    dmass,dpt,dy,
     &                    im1,im2,ipt1,ipt2,iy1,iy2


        integer*4  ndiv_sam,ma_max,ma_min
        logical read_scaler,mult_cut,apply_wt,proj_bends
        real*4  wtcut


        common/options/ndiv_sam,read_scaler,proj_bends,
     &                 mult_cut,ma_max,ma_min,
     &                 apply_wt,wtcut


        common/cnf_map/ibcf_os_pcf,ibcf_ss_pcf
        integer ibcf_os_pcf(8),ibcf_ss_pcf(8)
*
*--> Acceptance table format
*
        Logical l_V93, l_95A, l_95P
        Common /ACC_FILE_TYPE/ l_V93, l_95A, l_95P
*
* ----------------- End of osbk_control.inc --------------------
