


#include <mex.h>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <Eigen/Eigen>
#include "performBDFDiffusion_hpp.hpp"


// #include <math.h>
// #include <memory.h>
// #include <matrix.h>
// #include <stdlib.h>
// #include <stdio.h>
// #include <string.h>



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  const mxArray *pItL = prhs[0];
  const mxArray *pSignal = prhs[1];
  const mxArray *pNumIter = prhs[2];

  // double *pItL11 = mxGetPr(pItL1);
  double *pSig = mxGetPr(prhs[1]);

  int NumIter = (int)*mxGetPr(pNumIter);
  int nRowsSignal = (int)mxGetM(pSignal);
  int nColsSignal = (int)mxGetN(pSignal);


  // mexPrintf("nRowsSignal %d\n", nRowsSignal);
  // mexPrintf("nColsSignal %d\n", nColsSignal);

  const mwSize *pDimsItL = mxGetDimensions(pItL);
  // mexPrintf("pDimsItL %d\n", pDimsItL[0]);

  // const mwSize *pDimsSignal = mxGetDimensions(pSignal);

  // Eigen::SparseMatrix<float> ItL = *pItL1;
  // Eigen::VectorXf SignalIn = &pSignal;



  // mexPrintf("ItLMat size [%d,%d]\n\n", ItLmat.rows(), ItLmat.cols());

  // mexPrintf("NumIter %d \n", NumIter);
  // mexPrintf("Num Dim ItL [%d, %d] \n", pDimsItL[0],pDimsItL[1]);
  // mexPrintf("Num Dim Signal [%d, %d]\n", pDimsSignal[0],pDimsSignal[1]);

  Eigen::MatrixXd SignalOut(nRowsSignal, NumIter);

  Eigen::Map<Eigen::VectorXd> SignalIn (pSig, nRowsSignal, nColsSignal);

  // mexPrintf("SignalIn %0.4f \n", SignalIn(0));

  if (pDimsItL[0] == 1){
    const mxArray *pItL0 = mxGetCell(pItL, 0);
    MatlabSparse ItL0 = matlab_to_eigen_sparse(pItL0);
    SignalOut = runBDF1(ItL0, SignalIn, NumIter);

  }
  else if (pDimsItL[0] == 2){
    const mxArray *pItL0 = mxGetCell(pItL, 0);
    const mxArray *pItL1 = mxGetCell(pItL, 1);
    MatlabSparse ItL0 = matlab_to_eigen_sparse(pItL0);
    MatlabSparse ItL1 = matlab_to_eigen_sparse(pItL1);

    SignalOut = runBDF2(ItL0, ItL1, SignalIn, NumIter);
  }
  else if (pDimsItL[0] == 3){
    const mxArray *pItL0 = mxGetCell(pItL, 0);
    const mxArray *pItL1 = mxGetCell(pItL, 1);
    const mxArray *pItL2 = mxGetCell(pItL, 2);
    MatlabSparse ItL0 = matlab_to_eigen_sparse(pItL0);
    MatlabSparse ItL1 = matlab_to_eigen_sparse(pItL1);
    MatlabSparse ItL2 = matlab_to_eigen_sparse(pItL2);

    SignalOut = runBDF3(ItL0, ItL1, ItL2, SignalIn, NumIter);

  }
  else if (pDimsItL[0] == 4){
    const mxArray *pItL0 = mxGetCell(pItL, 0);
    const mxArray *pItL1 = mxGetCell(pItL, 1);
    const mxArray *pItL2 = mxGetCell(pItL, 2);
    const mxArray *pItL3 = mxGetCell(pItL, 3);
    MatlabSparse ItL0 = matlab_to_eigen_sparse(pItL0);
    MatlabSparse ItL1 = matlab_to_eigen_sparse(pItL1);
    MatlabSparse ItL2 = matlab_to_eigen_sparse(pItL2);
    MatlabSparse ItL3 = matlab_to_eigen_sparse(pItL3);

    SignalOut = runBDF4(ItL0, ItL1, ItL2, ItL3, SignalIn, NumIter);
  }


  

  plhs[0] = mxCreateDoubleMatrix(nRowsSignal, NumIter, mxREAL);
  Eigen::Map<Eigen::MatrixXd> DiffusedSignal (mxGetPr(plhs[0]), nRowsSignal, NumIter);
  DiffusedSignal = SignalOut;



  // mwSize rows = SignalOut.rows();
  // mwSize cols = SignalOut.cols();
  // plhs[0] = mxCreateDoubleMatrix(rows, cols, mxREAL);
  // Eigen::Map<Eigen::MatrixXd> map(mxGetPr(plhs[0]), rows, cols); // Map the array
  // map = SignalOut;
}
