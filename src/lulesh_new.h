#include "vector.h"
#include "typedef.h"
#define LULESH_SHOW_PROGRESS 0
#define DOUBLE_PRECISION 0
//#define SAMI 


#if USE_MPI
#include <mpi.h>

/*
   define one of these three symbols:

   SEDOV_SYNC_POS_VEL_NONE
   SEDOV_SYNC_POS_VEL_EARLY
   SEDOV_SYNC_POS_VEL_LATE
*/

// TODO: currently we support only early sync!
#define SEDOV_SYNC_POS_VEL_EARLY 1
#endif

enum {
  VolumeError = -1,
  QStopError = -2,
  LFileError = -3
} ;

/* Could also support fixed point and interval arithmetic types */
//typedef float        real4 ;
typedef float        real4 ;
//typedef double       real8 ;
typedef double real8 ;

typedef int    Index_t ; /* array subscript and loop index */
typedef int    Int_t ;   /* integer representation */
#ifdef DOUBLE_PRECISION
typedef real8  Real_t ;  /* floating point representation */
#else
typedef real4  Real_t ;  /* floating point representation */
#endif

class Domain
{

public: 

  void sortRegions(Vector_h<Int_t>& regReps_h, Vector_h<Index_t>& regSorted_h);
  void CreateRegionIndexSets(Int_t nr, Int_t balance);


  Index_t max_streams;
  std::vector<cudaStream_t> streams;

  /* Elem-centered */

  Vector_d<Index_t> matElemlist ; /* material indexset */
  Vector_d<Index_t> nodelist ;    /* elemToNode connectivity */

  Vector_d<Index_t> lxim ;        /* element connectivity through face */
  Vector_d<Index_t> lxip ;
  Vector_d<Index_t> letam ;
  Vector_d<Index_t> letap ;
  Vector_d<Index_t> lzetam ;
  Vector_d<Index_t> lzetap ;

  Vector_d<Int_t> elemBC ;        /* elem face symm/free-surf flag */

  Vector_d<Real_t_e> e ;            /* energy */

  Vector_d<Real_t_p> p ;            /* pressure */

  Vector_d<Real_t_q> q ;            /* q */
  Vector_d<Real_t_ql> ql ;           /* linear term for q */
  Vector_d<Real_t_qq> qq ;           /* quadratic term for q */

  Vector_d<Real_t_v> v ;            /* relative volume */

  Vector_d<Real_t_volo> volo ;         /* reference volume */
  Vector_d<Real_t_delv> delv ;         /* m_vnew - m_v */
  Vector_d<Real_t_vdov> vdov ;         /* volume derivative over volume */

  Vector_d<Real_t_arealg> arealg ;       /* char length of an element */
  
  Vector_d<Real_t_ss> ss ;           /* "sound speed" */

  Vector_d<Real_t_elemNass> elemMass ;     /* mass */

  Vector_d<Real_t_vnew>* vnew ;         /* new relative volume -- temporary */

  Vector_d<Real_t_delv_xi>* delv_xi ;      /* velocity gradient -- temporary */
  Vector_d<Real_t_delv_eta>* delv_eta ;
  Vector_d<Real_t_delv_zeta>* delv_zeta ;

  Vector_d<Real_t_delx_xi>* delx_xi ;      /* coordinate gradient -- temporary */
  Vector_d<Real_t_delx_eta>* delx_eta ;
  Vector_d<Real_t_delx_zeta>* delx_zeta ;

  Vector_d<Real_t_dxx>* dxx ;          /* principal strains -- temporary */
  Vector_d<Real_t_dyy>* dyy ;
  Vector_d<Real_t_dzz>* dzz ;

  /* Node-centered */

  Vector_d<Real_t_x> x ;            /* coordinates */
  Vector_d<Real_t_y> y ;
  Vector_d<Real_t_z> z ;

  Vector_d<Real_t_xd> xd ;           /* velocities */
  Vector_d<Real_t_yd> yd ;
  Vector_d<Real_t_zd> zd ;


  Vector_d<Real_t_xdd> xdd ;          /* accelerations */
  Vector_d<Real_t_ydd> ydd ;
  Vector_d<Real_t_zdd> zdd ;

  Vector_d<Real_t_fx> fx ;           /* forces */
  Vector_d<Real_t_fy> fy ;
  Vector_d<Real_t_fz> fz ;

  Vector_d<Real_t_nodalMass> nodalMass ;    /* mass */
  Vector_h<Real_t_h_nodalMass> h_nodalMass ;    /* mass - host */

  /* device pointers for comms */
  Real_t_d_delv_xi   *d_delv_xi ;      /* velocity gradient -- temporary */
  Real_t_d_delv_eta  *d_delv_eta ;
  Real_t_d_delv_zeta *d_delv_zeta ;

  Real_t_d_x *d_x ;            /* coordinates */
  Real_t_d_y *d_y ;
  Real_t_d_z *d_z ;

  Real_t_d_xd *d_xd ;           /* velocities */
  Real_t_d_yd *d_yd ;
  Real_t_d_zd *d_zd ;

  Real_t_d_fx *d_fx ;           /* forces */
  Real_t_d_fy *d_fy ;
  Real_t_d_fz *d_fz ;

  /* access elements for comms */
  Real_t_d_delv_xi&   get_delv_xi(Index_t idx) { return d_delv_xi[idx] ; }
  Real_t_d_delv_eta&  get_delv_eta(Index_t idx) { return d_delv_eta[idx] ; }
  Real_t_d_delv_zeta& get_delv_zeta(Index_t idx) { return d_delv_zeta[idx] ; }

  Real_t_d_x& get_x(Index_t idx) { return d_x[idx] ; }
  Real_t_d_y& get_y(Index_t idx) { return d_y[idx] ; }
  Real_t_d_z& get_z(Index_t idx) { return d_z[idx] ; }

  Real_t_d_xd& get_xd(Index_t idx) { return d_xd[idx] ; }
  Real_t_d_yd& get_yd(Index_t idx) { return d_yd[idx] ; }
  Real_t_d_zd& get_zd(Index_t idx) { return d_zd[idx] ; }

  Real_t_d_fx& get_fx(Index_t idx) { return d_fx[idx] ; }
  Real_t_d_fy& get_fy(Index_t idx) { return d_fy[idx] ; }
  Real_t_d_fz& get_fz(Index_t idx) { return d_fz[idx] ; }

  // host access
  Real_t_h_nodalMass& get_nodalMass(Index_t idx) { return h_nodalMass[idx] ; }

  /* Boundary nodesets */

  Vector_d<Index_t> symmX ;       /* symmetry plane nodesets */
  Vector_d<Index_t> symmY ;        
  Vector_d<Index_t> symmZ ;
   
  Vector_d<Int_t> nodeElemCount ;
  Vector_d<Int_t> nodeElemStart;
  Vector_d<Index_t> nodeElemCornerList ;

  /* Parameters */

  Real_t_dtfixed dtfixed ;               /* fixed time increment */
  Real_t_deltatimemultlb deltatimemultlb ;
  Real_t_deltatimemultub deltatimemultub ;
  Real_t_stoptime stoptime ;              /* end time for simulation */
  Real_t_dtmax dtmax ;                 /* maximum allowable time increment */
  Int_t cycle ;                  /* iteration count for simulation */

  Real_t_dthydro_h*   dthydro_h;             /* hydro time constraint */ 
  Real_t_dtcourant_h* dtcourant_h;           /* courant time constraint */
  Index_t* bad_q_h;              /* flag to indicate Q error */
  Index_t* bad_vol_h;            /* flag to indicate volume error */

  /* cuda Events to indicate completion of certain kernels */
  cudaEvent_t time_constraint_computed;

  Real_t_time_h time_h ;               /* current time */
  Real_t_deltatime_h deltatime_h ;          /* variable time increment */

  Real_t_u_cut u_cut ;                /* velocity tolerance */
  Real_t_hg_coef hgcoef ;               /* hourglass control */
  Real_t_q_stop qstop ;                /* excessive q indicator */
  Real_t_monoq_max_slope monoq_max_slope ;
  Real_t_monoq_limiter_mult monoq_limiter_mult ;   
  Real_t_e_cut e_cut ;                /* energy tolerance */
  Real_t_p_cut p_cut ;                /* pressure tolerance */
  Real_t_ss4o3 ss4o3 ;
  Real_t_q_cut q_cut ;                /* q tolerance */
  Real_t_v_cut v_cut ;                /* relative volume tolerance */
  Real_t_qlc_monoq qlc_monoq ;            /* linear term coef for q */
  Real_t_qqc_monoq qqc_monoq ;            /* quadratic term coef for q */
  Real_t_qqc qqc ;
  Real_t_eosvmax eosvmax ;
  Real_t_eosvmin eosvmin ;
  Real_t_pmin pmin ;                 /* pressure floor */
  Real_t_emin emin ;                 /* energy floor */
  Real_t_dvovmax dvovmax ;              /* maximum allowable volume change */
  Real_t_refdens refdens ;              /* reference density */

   Index_t m_colLoc ;
   Index_t m_rowLoc ;
   Index_t m_planeLoc ;
   Index_t m_tp ;

   Index_t&  colLoc()             { return m_colLoc ; }
   Index_t&  rowLoc()             { return m_rowLoc ; }
   Index_t&  planeLoc()           { return m_planeLoc ; }
   Index_t&  tp()                 { return m_tp ; }

  Index_t sizeX ;
  Index_t sizeY ;
  Index_t sizeZ ;
  Index_t maxPlaneSize ;
  Index_t maxEdgeSize ;

  Index_t numElem ;
  Index_t padded_numElem ; 

  Index_t numNode;
  Index_t padded_numNode ; 

  Index_t numSymmX ; 
  Index_t numSymmY ; 
  Index_t numSymmZ ; 

  Index_t octantCorner;

   // Region information
   Int_t numReg ; //number of regions (def:11)
   Int_t balance; //Load balance between regions of a domain (def: 1)
   Int_t  cost;  //imbalance cost (def: 1)
   Int_t*   regElemSize ;   // Size of region sets
   Vector_d<Int_t> regCSR;  // records the begining and end of each region
   Vector_d<Int_t> regReps; // records the rep number per region
   Vector_d<Index_t> regNumList;    // Region number per domain element
   Vector_d<Index_t> regElemlist;  // region indexset 
   Vector_d<Index_t> regSorted; // keeps index of sorted regions
   
   //
   // MPI-Related additional data
   //

   Index_t m_numRanks;
   Index_t& numRanks() { return m_numRanks ; }

   void SetupCommBuffers(Int_t edgeNodes);
   void BuildMesh(Int_t nx, Int_t edgeNodes, Int_t edgeElems, Int_t domNodes, Int_t padded_domElems, Vector_h<Real_t> &x_h, Vector_h<Real_t> &y_h, Vector_h<Real_t> &z_h, Vector_h<Int_t> &nodelist_h);

   // Used in setup
   Index_t m_rowMin, m_rowMax;
   Index_t m_colMin, m_colMax;
   Index_t m_planeMin, m_planeMax ;

#if USE_MPI   
   // Communication Work space 
   Real_t *commDataSend ;
   Real_t *commDataRecv ;

   Real_t *d_commDataSend ;
   Real_t *d_commDataRecv ;

   // Maximum number of block neighbors 
   MPI_Request recvRequest[26] ; // 6 faces + 12 edges + 8 corners 
   MPI_Request sendRequest[26] ; // 6 faces + 12 edges + 8 corners 
#endif

};

typedef Real_t& (Domain::* Domain_member )(Index_t) ;

// Assume 128 byte coherence
// Assume Real_t is an "integral power of 2" bytes wide
#define CACHE_COHERENCE_PAD_REAL (128 / sizeof(Real_t))

#define CACHE_ALIGN_REAL(n) \
   (((n) + (CACHE_COHERENCE_PAD_REAL - 1)) & ~(CACHE_COHERENCE_PAD_REAL-1))

// MPI Message Tags
#define MSG_COMM_SBN      1024
#define MSG_SYNC_POS_VEL  2048
#define MSG_MONOQ         3072

#define MAX_FIELDS_PER_MPI_COMM 6

// cpu-comms
void CommRecv(Domain& domain, Int_t msgType, Index_t xferFields,
              Index_t dx, Index_t dy, Index_t dz,
              bool doRecv, bool planeOnly);
void CommSend(Domain& domain, Int_t msgType,
              Index_t xferFields, Domain_member *fieldData,
              Index_t dx, Index_t dy, Index_t dz,
              bool doSend, bool planeOnly);
void CommSBN(Domain& domain, Int_t xferFields, Domain_member *fieldData);
void CommSyncPosVel(Domain& domain);
void CommMonoQ(Domain& domain);

// gpu-comms
void CommSendGpu(Domain& domain, Int_t msgType,
              Index_t xferFields, Domain_member *fieldData,
              Index_t dx, Index_t dy, Index_t dz,
              bool doSend, bool planeOnly, cudaStream_t stream);
void CommSBNGpu(Domain& domain, Int_t xferFields, Domain_member *fieldData, cudaStream_t *streams);
void CommSyncPosVelGpu(Domain& domain, cudaStream_t *streams);
void CommMonoQGpu(Domain& domain, cudaStream_t stream);
