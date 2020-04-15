#include<stdio.h>

#include "read_amr_c.h"

int main()
{
  printf("read_amr_c starting\n");

  // this would be needed for parallel execution only
  wrapamr_init_mpi(); 

  // read header if only a partial set of variables are to be read
  // and the set is based on information stored in the header.
  // Usually one can call wrapamr_read_file directly.
  // the first integer should be larger than the length of the filename
  // the last integer is 1 for verbose, 0 otherwise
  wrapamr_read_header("data/3d__all_3_t00000010_n0000059.idl",200,0);

  printf("read_amr_c read header\n");

  int nDim;
  wrapamr_get_ndim(&nDim);
  printf("nDim= %i\n", nDim);

  int nVar;
  wrapamr_get_nvar(&nVar);
  printf("nVar= %i\n", nVar);
  fflush(stdout);

  // read data from file. 
  // the first integer should be larger than the length of the filename
  // the second integer is 1 to read in header and data, 0 for reading data only.
  // the last integer is 1 for verbose, 0 otherwise
  wrapamr_read_file("data/3d__all_3_t00000010_n0000059.idl",200,0,1);

  printf("read_amr_c read data file\n");

  int l;
  char NameVar[500];
  wrapamr_get_namevar(NameVar, &l);
  printf("length = %i NameVar= %s\n", l, NameVar);

  char NameUnit[500];
  wrapamr_get_nameunit(NameUnit, &l);
  printf("length = %i NameUnit= %s\n", l, NameUnit);

  // Get coordinate limits. Note that these have nDim elements
  double CoordMin_D[nDim];
  double CoordMax_D[nDim];
  wrapamr_get_domain(CoordMin_D, CoordMax_D);
  int i;
  for (i=0; i<nDim; i++){
    printf("i, CoordMin[i], CoordMax[i]= %i %f, %f\n", i, CoordMin_D[i], CoordMax_D[i]);
  }

  // Get interpolated value at some position near the right side
  double x_D[nDim];
  double State_V[nVar];
  int iFound;
  printf("x_D=");
  for (i=0; i<nDim; i++){
    x_D[i] = 0.1* CoordMin_D[i] + 0.9*CoordMax_D[i];
    printf("%f ",x_D[i]);
  }
  printf("\n");

  wrapamr_get_data_serial(x_D, State_V, &iFound);

  printf("State_V=");
  for (i=0; i<nVar; i++){
    printf("%f ",State_V[i]);
  }
  printf("\n");
  printf("iFound= %i\n", iFound);

  // Get weight (index 0) and interpolated value (index 1..nVar)
  // This works in parallel mode too:
  //    have to add up contributions for all processors
  double WeightState_V[nVar+1];

  wrapamr_get_data(x_D, WeightState_V, &iFound);

  printf("WeightState_V=");
  for (i=0; i<nVar+1; i++){
    printf("%f ",WeightState_V[i]);
  }
  printf("\n");
  printf("iFound= %i\n", iFound);

  // Get weight, interpolated value and cell size.
  double CellSize_D[nDim];

  wrapamr_get_data_cell(x_D, WeightState_V, CellSize_D, &iFound);

  // Divide by the weight
  for (i=1; i<nVar+1; i++){
    WeightState_V[i] /= WeightState_V[0];
  }
  
  printf("WeightState_V=");
  for (i=0; i<nVar+1; i++){
    printf("%f ",WeightState_V[i]);
  }
  printf("\n");
  printf("CellSize_D=");
  for (i=0; i<nDim; i++){
    printf("%f ",CellSize_D[i]);
  }
  printf("\n");
  printf("iFound= %i\n", iFound);

  // Clean storage
  wrapamr_clean();

  return 0;
}
