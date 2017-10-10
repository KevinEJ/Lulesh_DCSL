CalcKine_kerenl = [   
'NULL'     , 'Real_t_deltatime_h', 'Real_t_volo',  'Real_t_delv',  'Real_t_vdov', 'Real_t_arealg',
'Real_t_v' , 'Real_t_dxx'        , 'Real_t_dyy' ,  'Real_t_dzz' ,  'Real_t_x'   , 'Real_t_volume',
'Real_t_B' , 'Real_t_D'   ]  
CalcVolumeFor_kernel = [
'Real_t_volo'    , 'Real_t_v'    , 'Real_t_x'      , 'Real_t_p'      ,  'Real_t_q'      ,  'Real_t_ss',  
'Real_t_elemMass', 'Real_t_hourg', 'Real_t_fx_elem', 'Real_t_fy_elem',  'Real_t_fz_elem',  'Real_t_dn',
'Real_t_dvd'     , 'Real_t_hg'   , 'Real_t_coeff'  , 'Real_t_sig'    ,  'Real_t_volume ']

ApplyMPAndUV_kernel = [
'Real_t_delv',  'Real_t_v'   ,  'Real_t_p'    ,  'Real_t_q'     ,  'Real_t_ss'   ,  'Real_t_refdens',
'Real_t_emin',  'Real_t_pmin',  'Real_t_e_cut',  'Real_t_p_cut' ,  'Real_t_q_cut',  'Real_t_v_cut'  ,
'Real_t_ql'  ,  'Real_t_qq'  ,  'Real_t_e'    ,  'Real_t_ss4o3' ,  'Real_t_comp' ,  'Real_t_work'   ,
'Real_t_ssc' ,  'Real_t_qtilde'  ] 

all_targets = [
'NULL',
'Real_t_deltatime_h', 
'Real_t_volo',
'Real_t_delv',
'Real_t_vdov',
'Real_t_arealg',
'Real_t_v',
'Real_t_dxx',
'Real_t_dyy',
'Real_t_dzz',
'Real_t_x',
'Real_t_volume',
'Real_t_B',
'Real_t_D', 
# ==============================
'Real_t_p',
'Real_t_q',
'Real_t_ss',
'Real_t_elemMass',
'Real_t_hourg',
'Real_t_fx_elem',
'Real_t_fy_elem',
'Real_t_fz_elem',
'Real_t_dn',
'Real_t_dvd',
'Real_t_hg',
'Real_t_coeff',
'Real_t_sig', 
#=============================
'Real_t_refdens',
'Real_t_emin',
'Real_t_pmin',
'Real_t_e_cut',
'Real_t_p_cut',
'Real_t_q_cut',
'Real_t_v_cut',
'Real_t_ql',
'Real_t_qq',
'Real_t_e',
'Real_t_ss4o3',
'Real_t_comp',
'Real_t_work',
'Real_t_ssc' ,  'Real_t_qtilde'  ] 
