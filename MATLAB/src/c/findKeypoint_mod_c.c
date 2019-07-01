


#include <mex.h>
#include <math.h>
// #include <vector.h>
#include <memory.h>
#include "matrix.h"
#include <matrix.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  if (nrhs != 2){
    mexErrMsgTxt("Two inputs are required");
  }

  if (nlhs != 2){
    mexErrMsgTxt("Two outputs are required");
  }

  // Initialization
  const mwSize *dimsDoG;
  const mwSize *dimsNeigh;
  mwSize *dimsCell;
  const mxArray *NeighCell = prhs[1];
  const mxArray *DoG = prhs[0];
  mxArray *NeighborsCell;
  double *NeighborsElem;
  mwIndex jvert;
  mwIndex jelem;
  mwIndex jDoG;
  mwIndex jneigh;
  mwSize NumNeighbors;
  double *DoGPr;
  double CurrentDoG;
  double NeighDoG;
  mwIndex CurrentNeigh;
  mwIndex jlevel;
  int LessThan;
  int GreatThan;
  unsigned long int count = 0;

  double *LocationIndexOut;
  double *LevelIndexOut;


  dimsDoG = mxGetDimensions(prhs[0]);
  dimsNeigh = mxGetDimensions(prhs[1]);

  unsigned long int *LocationIndex = malloc(sizeof(unsigned long int)*dimsDoG[0]*dimsDoG[1]);
  unsigned long int *LevelIndex = malloc(sizeof(unsigned long int)*dimsDoG[0]*dimsDoG[1]);


  // mexPrintf("DoG Dims [%d, %d]\n", dimsDoG[0], dimsDoG[1]);
  // mexPrintf("Neigh Dims [%d, %d]\n", dimsNeigh[0], dimsNeigh[1]);


  DoGPr = mxGetPr(DoG);

  // loop over the vertex locations
  for (jvert=0; jvert < dimsNeigh[0]; jvert++){
  // for (jvert=0; jvert < 10; jvert++){

    NeighborsCell = mxGetCell(NeighCell, jvert);
    NumNeighbors = mxGetNumberOfElements(NeighborsCell);
    NeighborsElem = mxGetPr(NeighborsCell);


    // loop over the DoG levels
    for (jDoG=1; jDoG < dimsDoG[1]-1; jDoG++){
    // for (jDoG=1; jDoG < 10; jDoG++){

      LessThan = 1;
      GreatThan = 1;
      CurrentDoG = DoGPr[dimsDoG[0]*jDoG + jvert];
      // CurrentNeigh = NeighborsElem[jneigh];
      // mexPrintf("Current DoG is %0.8f\n", CurrentDoG);

      // mexPrintf("Current DoG %0.8f \n", CurrentDoG);
      // mexPrintf("Before Level DoG %0.8f, Next Level DoG %0.8f \n", DoGPr[dimsDoG[0]*(jDoG-1) + jvert], DoGPr[dimsDoG[0]*(jDoG+1) + jvert]);
      //
      // for (jlevel = 0; jlevel < 3; jlevel++){
      //   for (jneigh = 0; jneigh < NumNeighbors; jneigh++){
      //     mexPrintf("Neigh Val %0.8f \n", DoGPr[dimsDoG[0]*(jDoG+jlevel-1) + (int)NeighborsElem[jneigh]-1]);
      //   }
      // }
        //if (NumNeighbors > 7){
            

          // find if less than all neighbors
          if (CurrentDoG < DoGPr[dimsDoG[0]*(jDoG-1) + jvert]){
            if (CurrentDoG < DoGPr[dimsDoG[0]*(jDoG+1) + jvert]){
                for (jneigh = 0; jneigh < NumNeighbors; jneigh++){
                  NeighDoG = DoGPr[dimsDoG[0]*(jDoG+1-1) + (int)NeighborsElem[jneigh]-1];
                  if (CurrentDoG > NeighDoG){
                    LessThan = 0;
                    break;
                  }
                  if (!LessThan){
                    break;
                  }
                }
            } else{ LessThan = 0; }
          } else{ LessThan = 0; }


      // find if greater than all Neighbors
          if (CurrentDoG > DoGPr[dimsDoG[0]*(jDoG-1) + jvert]){
            if (CurrentDoG > DoGPr[dimsDoG[0]*(jDoG+1) + jvert]){
                for (jneigh = 0; jneigh < NumNeighbors; jneigh++){
                  NeighDoG = DoGPr[dimsDoG[0]*(jDoG+1-1) + (int)NeighborsElem[jneigh]-1];
                  if (CurrentDoG < NeighDoG){
                    GreatThan = 0;
                    break;
                  }
                  if (!GreatThan){
                    break;
                  }
                }
            } else{ GreatThan = 0; }
          } else{ GreatThan = 0; }

          if (LessThan || GreatThan){
            LocationIndex[count] = jvert + 1;
            LevelIndex[count] = jDoG + 1;
            count++;
            // mexPrintf("Keypoint at Vertex %d,  Level %d \n", jvert+1, jDoG+1);
          }


      // mexPrintf("DoG value DoG[%d,%d] = %0.8f\n", jvert+1, jDoG+1, CurrentDoG);

      // loop over the elements in each cell
      // for (jelem=0; jelem < NumElem; jelem++){}
       // } //else{ mexPrintf("Num Neighbors %d\n", NumNeighbors); }

    }

      // mexPrintf("cell %d element is %d \n", jvert, NumElem);
      // mexPrintf("cell %d element is %g \n", jvert, NeighborsElem[jelem]);

  }


  // Place keypoint locations and levels for output

  plhs[0] = mxCreateDoubleMatrix(count, 1, mxREAL );
  plhs[1] = mxCreateDoubleMatrix(count, 1, mxREAL );

  LocationIndexOut = mxGetPr(plhs[0]);
  LevelIndexOut = mxGetPr(plhs[1]);

  for (mwIndex c = 0; c < count; c++){
    LocationIndexOut[c] = LocationIndex[c];
    LevelIndexOut[c] = LevelIndex[c];
  }

  // Release memory of malloc variables
  free(LocationIndex);
  free(LevelIndex);


  // printf("Size of DoG is (%d,%d).\n ", dimsDoG[0], dimsDoG[1]);
  // printf("Size of Neigh is (%d,%d).\n ", dimsNeigh[0], dimsNeigh[1]);

 // mexPrintf("count %lu \n", count);
}
