#ifndef TYPES
#define TYPES
  
#define DOUBLE_PRECISION 1

typedef float Real_t_deltatime_h ; //CalcKine_kerenl, 
typedef float Real_t_volo        ; //CalcKine_kerenl, CalcVolumeFor_kernel 
typedef float Real_t_delv        ; //CalcKine_kerenl, ApplyMPAndUV_kernel 
typedef float Real_t_vdov        ; //CalcKine_kerenl, 
typedef float Real_t_arealg      ; //CalcKine_kerenl, 
typedef float Real_t_v           ; //CalcKine_kerenl, CalcVolumeFor_kernel ApplyMPAndUV_kernel 
typedef double Real_t_dxx         ; //CalcKine_kerenl, 
typedef float Real_t_dyy         ; //CalcKine_kerenl, 
typedef float Real_t_dzz         ; //CalcKine_kerenl, 
typedef float Real_t_x           ; //CalcKine_kerenl,  CalcVolumeFor_kernel 
typedef double Real_t_volume      ; //CalcKine_kerenl,  CalcVolumeFor_kernel 
typedef double Real_t_B           ; //CalcKine_kerenl, 
typedef double Real_t_D           ; //CalcKine_kerenl, 
//***********************************
typedef float Real_t_p            ; //CalcVolumeFor_kernel ApplyMPAndUV_kernel 
typedef float Real_t_q            ; //CalcVolumeFor_kernel ApplyMPAndUV_kernel 
typedef float Real_t_ss           ; //CalcVolumeFor_kernel ApplyMPAndUV_kernel 
typedef float Real_t_elemMass     ; //CalcVolumeFor_kernel 
typedef double Real_t_hourg        ; //CalcVolumeFor_kernel 
typedef float Real_t_fx_elem      ; //CalcVolumeFor_kernel 
typedef float Real_t_fy_elem      ; //CalcVolumeFor_kernel 
typedef float Real_t_fz_elem      ; //CalcVolumeFor_kernel 
typedef double Real_t_dn           ; //CalcVolumeFor_kernel 
typedef double Real_t_dvd          ; //CalcVolumeFor_kernel
typedef double Real_t_hg           ; //CalcVolumeFor_kernel 
typedef float Real_t_coeff        ; //CalcVolumeFor_kernel 
typedef float Real_t_sig          ; //CalcVolumeFor_kernel
//***********************************
typedef float Real_t_refdens     ; //ApplyMate_kernel
typedef float Real_t_emin        ; //ApplyMate_kerne;
typedef double Real_t_pmin        ; //ApplyMate_kerne; 
typedef double Real_t_e_cut       ; //ApplyMate_kerne; 
typedef float Real_t_p_cut       ; //ApplyMate_kerne; 
typedef float Real_t_q_cut       ; //ApplyMate_kerne; 
typedef double Real_t_v_cut       ; //ApplyMate_kerne; 
typedef double Real_t_ql          ; //ApplyMate_kerne; 
typedef float Real_t_qq          ; //ApplyMate_kerne; 
typedef float Real_t_e           ; //ApplyMate_kerne; 
typedef float Real_t_ss4o3       ; //ApplyMate_kerne;
typedef float Real_t_comp        ; //ApplyMate_kerne;
typedef double Real_t_work        ; //ApplyMate_kerne;
typedef float Real_t_ssc         ; //ApplyMate_kerne;           
typedef float Real_t_qtilde      ; //ApplyMate_kerne;           

typedef Real_t_volume Real_t_detJ ; 
typedef Real_t_x Real_t_delv_xi   ; 
typedef Real_t_x Real_t_delv_eta  ; 
typedef Real_t_x Real_t_delv_zeta ; 
typedef Real_t_x Real_t_delx_xi   ; 
typedef Real_t_x Real_t_delx_eta  ; 
typedef Real_t_x Real_t_delx_zeta ; 
typedef Real_t_v Real_t_vnew      ;
typedef Real_t_v Real_t_eosvmin      ;
typedef Real_t_v Real_t_eosvmax      ;

//typedef double Real_t_y           ;
//typedef double Real_t_z           ;
//typedef double Real_t_xd          ; 
//typedef double Real_t_yd          ;
//typedef double Real_t_zd          ;
//typedef double Real_t_delv_xi     ; 
//typedef double Real_t_delv_eta    ; 
//typedef double Real_t_delv_zeta   ; 
//typedef double Real_t_delx_xi     ; 
//typedef double Real_t_delx_eta    ; 
//typedef double Real_t_delx_zeta   ; 
//typedef double Real_t_vnew        ; //CalcKine_kerenl, ApplyMPAndUV_kernel 
/*
Vector_d<Real_t_e> e ;            
Vector_d<Real_t_p> p ;            
Vector_d<Real_t_q> q ;            
Vector_d<Real_t_ql> ql ;          
Vector_d<Real_t_qq> qq ;          
Vector_d<Real_t_v> v ;           
Vector_d<Real_t_volo> volo ;        
Vector_d<Real_t_delv> delv ;        
Vector_d<Real_t_vdov> vdov ;        
Vector_d<Real_t_arealg> arealg ;    
Vector_d<Real_t_ss> ss ;           
Vector_d<Real_t_elemNass> elemMass ;
Vector_d<Real_t_vnew>* vnew ;       
Vector_d<Real_t_delv_xi>* delv_xi ; 
Vector_d<Real_t_delv_eta>* delv_eta 
Vector_d<Real_t_delv_zeta>* delv_zet
Vector_d<Real_t_delx_xi>* delx_xi ; 
Vector_d<Real_t_delx_eta>* delx_eta 
Vector_d<Real_t_delx_zeta>* delx_zet
Vector_d<Real_t_dxx>* dxx ;         
Vector_d<Real_t_dyy>* dyy ;
Vector_d<Real_t_dzz>* dzz ;
Vector_d<Real_t_x> x ;            
Vector_d<Real_t_y> y ;
Vector_d<Real_t_z> z ;
Vector_d<Real_t_xd> xd ;           
Vector_d<Real_t_yd> yd ;
Vector_d<Real_t_zd> zd ;
Vector_d<Real_t_xdd> xdd ;          
Vector_d<Real_t_ydd> ydd ;
Vector_d<Real_t_zdd> zdd ;
Vector_d<Real_t_fx> fx ;           
Vector_d<Real_t_fy> fy ;
Vector_d<Real_t_fz> fz ;
Vector_d<Real_t_nodalMass> nodalMass
Vector_h<Real_t_h_nodalMass> h_nodal
Real_t_d_delv_xi   *d_delv_xi ;     
Real_t_d_delv_eta  *d_delv_eta ;
Real_t_d_delv_zeta *d_delv_zeta ;
Real_t_d_x *d_x ;            
Real_t_d_y *d_y ;
Real_t_d_z *d_z ;
Real_t_d_xd *d_xd ;         
Real_t_d_yd *d_yd ;
Real_t_d_zd *d_zd ;
Real_t_d_fx *d_fx ;           
Real_t_d_fy *d_fy ;
Real_t_d_fz *d_fz ;
Real_t_d_delv_xi&   get_delv_xi(Inde
Real_t_d_delv_eta&  get_delv_eta(Ind
Real_t_d_delv_zeta& get_delv_zeta(In
Real_t_d_x& get_x(Index_t idx) { ret
Real_t_d_y& get_y(Index_t idx) { ret
Real_t_d_z& get_z(Index_t idx) { ret
Real_t_d_xd& get_xd(Index_t idx) { r
Real_t_d_yd& get_yd(Index_t idx) { r
Real_t_d_zd& get_zd(Index_t idx) { r
Real_t_d_fx& get_fx(Index_t idx) { r
Real_t_d_fy& get_fy(Index_t idx) { r
Real_t_d_fz& get_fz(Index_t idx) { r
Real_t_h_nodalMass& get_nodalMass(In
Real_t_dtfixed dtfixed ;            
Real_t_deltatimemultlb deltatimemult
Real_t_deltatimemultub deltatimemult
Real_t_stoptime stoptime ;          
Real_t_dtmax dtmax ;                
Real_t_dthydro_h*   dthydro_h;      
Real_t_dtcourant_h* dtcourant_h;    
Real_t_time_h time_h ;              
Real_t_deltatime_h deltatime_h ;    
Real_t_u_cut u_cut ;                
Real_t_hg_coef hgcoef ;             
Real_t_q_stop qstop ;               
Real_t_monoq_max_slope monoq_max_slo
Real_t_monoq_limiter_mult monoq_limi
Real_t_e_cut e_cut ;                
Real_t_p_cut p_cut ;                
Real_t_ss4o3 ss4o3 ;
Real_t_q_cut q_cut ;                
Real_t_v_cut v_cut ;                
Real_t_qlc_monoq qlc_monoq ;        
Real_t_qqc_monoq qqc_monoq ;        
Real_t_qqc qqc ;
Real_t_eosvmax eosvmax ;
Real_t_eosvmin eosvmin ;
Real_t_pmin pmin ;                 
Real_t_emin emin ;                 
Real_t_dvovmax dvovmax ;            
Real_t_refdens refdens ;    */

#endif

#ifndef EJ_Time_Start
#define EJ_Time_Start\
    clock_t t = clock() ;\
    timeval EJstart;\
    gettimeofday(&EJstart, NULL) ;\
    cudaEvent_t EJJstart, EJJstop;\
    cudaEventCreate(&EJJstart);\
    cudaEventCreate(&EJJstop);\
    cudaEventRecord(EJJstart);
#endif
#ifndef EJ_Time_End
#define EJ_Time_End\
    cudaEventRecord(EJJstop);\
    cudaEventSynchronize(EJJstop);\
    float milliseconds = 0;\
    cudaEventElapsedTime(&milliseconds, EJJstart, EJJstop);\
    t = clock() - t ;
#endif

