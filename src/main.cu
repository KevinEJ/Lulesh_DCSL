
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <sstream>

#include <util.h>
#include <sm_utils.inl>
#include <cuda.h>
#include <allocator.h>
#include "cuda_profiler_api.h"

#ifdef USE_MPI
#include <mpi.h>
#endif

#include <sys/time.h>
#include <unistd.h>

#include "lulesh.h"



int main(int argc, char *argv[])
{
  if (argc < 3) {
    printUsage(argv);
    exit( LFileError );
  }
  
  if ( strcmp(argv[1],"-u") != 0 && strcmp(argv[1],"-s") != 0 ) 
  {
    printUsage(argv);
    exit( LFileError ) ;
  }
  //EJ 
  int num_iters = -1;
  /*if (argc == 5) {
    num_iters = atoi(argv[4]);
  }*/
  if( argc < 5 || strcmp(argv[3],"-b")!=0){
    printf( "usage ./lulesh -s elems -n blocksize \n" ) ;
    exit( LFileError ) ;
  }
  global_block_size = atoi(argv[4]) ;
  /*if( argc < 7 || strcmp(argv[5],"-i")!=0){
    printf( "usage ./lulesh -s elems -n blocksize -i iters \n" ) ;
    exit( LFileError ) ;
  }*/
  if( argc == 7 && strcmp(argv[5],"-i")==0){
    num_iters = atoi(argv[6]) ;
  }
  //EJ end
  bool structured = ( strcmp(argv[1],"-s") == 0 );

  Int_t numRanks ;
  Int_t myRank ;

#if USE_MPI   
  Domain_member fieldData ;

  MPI_Init(&argc, &argv) ;
  MPI_Comm_size(MPI_COMM_WORLD, &numRanks) ;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank) ;
#else
  numRanks = 1;
  myRank = 0;
#endif

  cuda_init(myRank);

  /* assume cube subdomain geometry for now */
  Index_t nx = atoi(argv[2]);

  Domain *locDom ;

  // Set up the mesh and decompose. Assumes regular cubes for now
  Int_t col, row, plane, side;
  InitMeshDecomp(numRanks, myRank, &col, &row, &plane, &side);

  // TODO: change default nr to 11
  Int_t nr = 11;
  Int_t balance = 1;
  Int_t cost = 1;

  // TODO: modify this constructor to account for new fields
  // TODO: setup communication buffers
  locDom = NewDomain(argv, numRanks, col, row, plane, nx, side, structured, nr, balance, cost); 

#if USE_MPI   
   // copy to the host for mpi transfer
   locDom->h_nodalMass = locDom->nodalMass;

   fieldData = &Domain::get_nodalMass;

   // Initial domain boundary communication 
   CommRecv(*locDom, MSG_COMM_SBN, 1,
            locDom->sizeX + 1, locDom->sizeY + 1, locDom->sizeZ + 1,
            true, false) ;
   CommSend(*locDom, MSG_COMM_SBN, 1, &fieldData,
            locDom->sizeX + 1, locDom->sizeY + 1, locDom->sizeZ + 1,
            true, false) ;
   CommSBN(*locDom, 1, &fieldData) ;

   // copy back to the device
   locDom->nodalMass = locDom->h_nodalMass;

   // End initialization
   MPI_Barrier(MPI_COMM_WORLD);
#endif

  cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);

  /* timestep to solution */
  int its=0;

  if (myRank == 0) {
    if (structured)
      printf("Running until t=%f, Problem size=%dx%dx%d\n",locDom->stoptime,nx,nx,nx);
    else 
      printf("Running until t=%f, Problem size=%d \n",locDom->stoptime,locDom->numElem);
  }

  cudaProfilerStart();

#if USE_MPI   
   double start = MPI_Wtime();
#else
   timeval start;
   gettimeofday(&start, NULL) ;
#endif

  while(locDom->time_h < locDom->stoptime)
  {
    // this has been moved after computation of volume forces to hide launch latencies
    //TimeIncrement(locDom) ;

    LagrangeLeapFrog(locDom) ;

    checkErrors(locDom,its,myRank);

    #if LULESH_SHOW_PROGRESS
     if (myRank == 0) 
	 printf("cycle = %d, time = %e, dt=%e\n", its+1, double(locDom->time_h), double(locDom->deltatime_h) ) ;
    #endif
    its++;
    if (its == num_iters) break;
  }

  // make sure GPU finished its work
  cudaDeviceSynchronize();

// Use reduced max elapsed time
   double elapsed_time;
#if USE_MPI   
   elapsed_time = MPI_Wtime() - start;
#else
   timeval end;
   gettimeofday(&end, NULL) ;
   elapsed_time = (double)(end.tv_sec - start.tv_sec) + ((double)(end.tv_usec - start.tv_usec))/1000000 ;
#endif

   double elapsed_timeG;
#if USE_MPI   
   MPI_Reduce(&elapsed_time, &elapsed_timeG, 1, MPI_DOUBLE,
              MPI_MAX, 0, MPI_COMM_WORLD);
#else
   elapsed_timeG = elapsed_time;
#endif

  cudaProfilerStop();

  if (myRank == 0) 
    VerifyAndWriteFinalOutput(elapsed_timeG, *locDom, its, nx, numRanks);

#ifdef SAMI
  DumpDomain(locDom) ;
#endif
  cudaDeviceReset();

#if USE_MPI
   MPI_Finalize() ;
#endif

  return 0 ;
}
