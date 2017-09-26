#include "lulesh.h"

void
Domain::SetupCommBuffers(Int_t edgeNodes)
{
  // allocate a buffer large enough for nodal ghost data 
  maxEdgeSize = MAX(this->sizeX, MAX(this->sizeY, this->sizeZ))+1 ;
  maxPlaneSize = CACHE_ALIGN_REAL(maxEdgeSize*maxEdgeSize) ;
  maxEdgeSize = CACHE_ALIGN_REAL(maxEdgeSize) ;

  // assume communication to 6 neighbors by default 
  m_rowMin = (m_rowLoc == 0)        ? 0 : 1;
  m_rowMax = (m_rowLoc == m_tp-1)     ? 0 : 1;
  m_colMin = (m_colLoc == 0)        ? 0 : 1;
  m_colMax = (m_colLoc == m_tp-1)     ? 0 : 1;
  m_planeMin = (m_planeLoc == 0)    ? 0 : 1;
  m_planeMax = (m_planeLoc == m_tp-1) ? 0 : 1;

#if USE_MPI   
  // account for face communication 
  Index_t comBufSize =
    (m_rowMin + m_rowMax + m_colMin + m_colMax + m_planeMin + m_planeMax) *
    maxPlaneSize * MAX_FIELDS_PER_MPI_COMM ;

  // account for edge communication 
  comBufSize +=
    ((m_rowMin & m_colMin) + (m_rowMin & m_planeMin) + (m_colMin & m_planeMin) +
     (m_rowMax & m_colMax) + (m_rowMax & m_planeMax) + (m_colMax & m_planeMax) +
     (m_rowMax & m_colMin) + (m_rowMin & m_planeMax) + (m_colMin & m_planeMax) +
     (m_rowMin & m_colMax) + (m_rowMax & m_planeMin) + (m_colMax & m_planeMin)) *
    maxPlaneSize * MAX_FIELDS_PER_MPI_COMM ;

  // account for corner communication 
  // factor of 16 is so each buffer has its own cache line 
  comBufSize += ((m_rowMin & m_colMin & m_planeMin) +
                 (m_rowMin & m_colMin & m_planeMax) +
                 (m_rowMin & m_colMax & m_planeMin) +
                 (m_rowMin & m_colMax & m_planeMax) +
                 (m_rowMax & m_colMin & m_planeMin) +
                 (m_rowMax & m_colMin & m_planeMax) +
                 (m_rowMax & m_colMax & m_planeMin) +
                 (m_rowMax & m_colMax & m_planeMax)) * CACHE_COHERENCE_PAD_REAL ;

  this->commDataSend = new Real_t[comBufSize] ;
  this->commDataRecv = new Real_t[comBufSize] ;

  // pin buffers
  cudaHostRegister(this->commDataSend, comBufSize*sizeof(Real_t), 0);
  cudaHostRegister(this->commDataRecv, comBufSize*sizeof(Real_t), 0);

  // prevent floating point exceptions 
  memset(this->commDataSend, 0, comBufSize*sizeof(Real_t)) ;
  memset(this->commDataRecv, 0, comBufSize*sizeof(Real_t)) ;

  // allocate shadow GPU buffers
  cudaMalloc(&this->d_commDataSend, comBufSize*sizeof(Real_t));
  cudaMalloc(&this->d_commDataRecv, comBufSize*sizeof(Real_t));
  
  // prevent floating point exceptions 
  cudaMemset(this->d_commDataSend, 0, comBufSize*sizeof(Real_t));
  cudaMemset(this->d_commDataRecv, 0, comBufSize*sizeof(Real_t));
#endif
}



void Domain::BuildMesh(Int_t nx, Int_t edgeNodes, Int_t edgeElems, Int_t domNodes, Int_t padded_domElems, Vector_h<Real_t> &x_h, Vector_h<Real_t> &y_h, Vector_h<Real_t> &z_h, Vector_h<Int_t> &nodelist_h)
{
  Index_t meshEdgeElems = m_tp*nx ;

  x_h.resize(domNodes);
  y_h.resize(domNodes);
  z_h.resize(domNodes);

  // initialize nodal coordinates 
  Index_t nidx = 0 ;
  Real_t tz = Real_t(1.125)*Real_t(m_planeLoc*nx)/Real_t(meshEdgeElems) ;
  for (Index_t plane=0; plane<edgeNodes; ++plane) {
    Real_t ty = Real_t(1.125)*Real_t(m_rowLoc*nx)/Real_t(meshEdgeElems) ;
    for (Index_t row=0; row<edgeNodes; ++row) {
      Real_t tx = Real_t(1.125)*Real_t(m_colLoc*nx)/Real_t(meshEdgeElems) ;
      for (Index_t col=0; col<edgeNodes; ++col) {
        x_h[nidx] = tx ;
        y_h[nidx] = ty ;
        z_h[nidx] = tz ;
        ++nidx ;
        // tx += ds ; // may accumulate roundoff... 
        tx = Real_t(1.125)*Real_t(m_colLoc*nx+col+1)/Real_t(meshEdgeElems) ;
      }
      // ty += ds ;  // may accumulate roundoff... 
      ty = Real_t(1.125)*Real_t(m_rowLoc*nx+row+1)/Real_t(meshEdgeElems) ;
    }
    // tz += ds ;  // may accumulate roundoff... 
    tz = Real_t(1.125)*Real_t(m_planeLoc*nx+plane+1)/Real_t(meshEdgeElems) ;
  }

  x = x_h;
  y = y_h;
  z = z_h;

  nodelist_h.resize(padded_domElems*8);

  // embed hexehedral elements in nodal point lattice 
  Index_t zidx = 0 ;
  nidx = 0 ;
  for (Index_t plane=0; plane<edgeElems; ++plane) {
    for (Index_t row=0; row<edgeElems; ++row) {
      for (Index_t col=0; col<edgeElems; ++col) {
        nodelist_h[0*padded_domElems+zidx] = nidx                                       ;
        nodelist_h[1*padded_domElems+zidx] = nidx                                   + 1 ;
        nodelist_h[2*padded_domElems+zidx] = nidx                       + edgeNodes + 1 ;
        nodelist_h[3*padded_domElems+zidx] = nidx                       + edgeNodes     ;
        nodelist_h[4*padded_domElems+zidx] = nidx + edgeNodes*edgeNodes                 ;
        nodelist_h[5*padded_domElems+zidx] = nidx + edgeNodes*edgeNodes             + 1 ;
        nodelist_h[6*padded_domElems+zidx] = nidx + edgeNodes*edgeNodes + edgeNodes + 1 ;
        nodelist_h[7*padded_domElems+zidx] = nidx + edgeNodes*edgeNodes + edgeNodes     ;
        ++zidx ;
        ++nidx ;
      }
      ++nidx ;
    }
    nidx += edgeNodes ;
  }

  nodelist = nodelist_h;
}


/*******************	to support region	*********************/
void Domain::sortRegions(Vector_h<Int_t>& regReps_h, Vector_h<Index_t>& regSorted_h)
{
  Index_t temp;
  Vector_h<Index_t> regIndex;
  regIndex.resize(numReg);
  for(int i = 0; i < numReg; i++)
	regIndex[i] = i;

  for(int i = 0; i < numReg-1; i++)
	for(int j = 0; j < numReg-i-1; j++)
		if(regReps_h[j] < regReps_h[j+1])
		{
                  temp = regReps_h[j];
                  regReps_h[j] = regReps_h[j+1];
                  regReps_h[j+1] = temp;

		  temp = regElemSize[j];
		  regElemSize[j] = regElemSize[j+1];
		  regElemSize[j+1] = temp;

                  temp = regIndex[j];
                  regIndex[j] = regIndex[j+1];
                  regIndex[j+1] = temp;
		}
  for(int i = 0; i < numReg; i++)
        regSorted_h[regIndex[i]] = i;
}


void Domain::CreateRegionIndexSets(Int_t nr, Int_t b)
{
#if USE_MPI   
   Index_t myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank) ;
   srand(myRank);
#else
   srand(0);
   Index_t myRank = 0;
#endif
   numReg = nr;
   balance = b;

   regElemSize = new Int_t[numReg];
   Index_t nextIndex = 0;

   Vector_h<Int_t> regCSR_h(regCSR.size());  // records the begining and end of each region
   Vector_h<Int_t> regReps_h(regReps.size()); // records the rep number per region
   Vector_h<Index_t> regNumList_h(regNumList.size());    // Region number per domain element
   Vector_h<Index_t> regElemlist_h(regElemlist.size());  // region indexset 
   Vector_h<Index_t> regSorted_h(regSorted.size()); // keeps index of sorted regions

   //if we only have one region just fill it
   // Fill out the regNumList with material numbers, which are always
   // the region index plus one 
   if(numReg == 1) {
      while (nextIndex < numElem) {
         regNumList_h[nextIndex] = 1;
         nextIndex++;
      }
      regElemSize[0] = 0;
   }
   //If we have more than one region distribute the elements.
   else {
      Int_t regionNum;
      Int_t regionVar;
      Int_t lastReg = -1;
      Int_t binSize;
      Int_t elements;
      Index_t runto = 0;
      Int_t costDenominator = 0;
      Int_t* regBinEnd = new Int_t[numReg];
      //Determine the relative weights of all the regions.
      for (Index_t i=0 ; i<numReg ; ++i) {
         regElemSize[i] = 0;
         costDenominator += POW((i+1), balance);  //Total cost of all regions
         regBinEnd[i] = costDenominator;  //Chance of hitting a given region is (regBinEnd[i] - regBinEdn[i-1])/costDenominator
      }
      //Until all elements are assigned
      while (nextIndex < numElem) {
         //pick the region
         regionVar = rand() % costDenominator;
         Index_t i = 0;
         while(regionVar >= regBinEnd[i])
            i++;
         //rotate the regions based on MPI rank.  Rotation is Rank % NumRegions
         regionNum = ((i + myRank) % numReg) + 1;
         // make sure we don't pick the same region twice in a row
         while(regionNum == lastReg) {
            regionVar = rand() % costDenominator;
            i = 0;
            while(regionVar >= regBinEnd[i])
               i++;
            regionNum = ((i + myRank) % numReg) + 1;
         }
         //Pick the bin size of the region and determine the number of elements.
         binSize = rand() % 1000;
         if(binSize < 773) {
           elements = rand() % 15 + 1;
         }
         else if(binSize < 937) {
           elements = rand() % 16 + 16;
         }
         else if(binSize < 970) {
           elements = rand() % 32 + 32;
         }
         else if(binSize < 974) {
           elements = rand() % 64 + 64;
         }
         else if(binSize < 978) {
           elements = rand() % 128 + 128;
         }
         else if(binSize < 981) {
           elements = rand() % 256 + 256;
         }
         else
            elements = rand() % 1537 + 512;
         runto = elements + nextIndex;
         //Store the elements.  If we hit the end before we run out of elements then just stop.
         while (nextIndex < runto && nextIndex < numElem) {
            regNumList_h[nextIndex] = regionNum;
            nextIndex++;
         }
         lastReg = regionNum;
      }
   }
   // Convert regNumList to region index sets
   // First, count size of each region 
   for (Index_t i=0 ; i<numElem ; ++i) {
      int r = regNumList_h[i]-1; // region index == regnum-1
      regElemSize[r]++;
   }

   Index_t rep;
   // Second, allocate each region index set
   for (Index_t r=0; r<numReg ; ++r) {
       if(r < numReg/2)
         rep = 1;
       else if(r < (numReg - (numReg+15)/20))
         rep = 1 + cost;
       else
         rep = 10 * (1+ cost);
       regReps_h[r] = rep;
   }

   sortRegions(regReps_h, regSorted_h);

   regCSR_h[0] = 0;
   // Second, allocate each region index set
   for (Index_t i=1 ; i<numReg ; ++i) {
      regCSR_h[i] = regCSR_h[i-1] + regElemSize[i-1];
   }

   // Third, fill index sets
   for (Index_t i=0 ; i<numElem ; ++i) {
      Index_t r = regSorted_h[regNumList_h[i]-1];       // region index == regnum-1
      regElemlist_h[regCSR_h[r]] = i;
      regCSR_h[r]++;
   }

   // Copy to device
   regCSR =  regCSR_h;  // records the begining and end of each region
   regReps =  regReps_h; // records the rep number per region
   regNumList =  regNumList_h;    // Region number per domain element
   regElemlist = regElemlist_h;  // region indexset 
   regSorted = regSorted_h; // keeps index of sorted regions

} // end of create function
