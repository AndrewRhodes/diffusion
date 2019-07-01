


#include <mex.h>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <Eigen/Eigen>
#include <unsupported/Eigen/MatrixFunctions>
#include "performExpMDiffusion_hpp.hpp"




void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  const mxArray *pLBM = prhs[0];
  const mxArray *pSignal = prhs[1];
  const mxArray *pNumIter = prhs[2];

  // double *pItL11 = mxGetPr(pItL1);
  double *pSig = mxGetPr(prhs[1]);

  int NumIter = (int)*mxGetPr(pNumIter);
  int nRowsSignal = (int)mxGetM(pSignal);
  int nColsSignal = (int)mxGetN(pSignal);


  // mexPrintf("nRowsSignal %d\n", nRowsSignal);
  // mexPrintf("nColsSignal %d\n", nColsSignal);

  const mwSize *pDimsLBM = mxGetDimensions(pLBM);
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


  // const mxArray *pLBM = mxGetCell(pLBM, 0);
  MatlabSparse LBM = matlab_to_eigen_sparse(pLBM);


  SignalOut = runExpM(LBM, SignalIn, NumIter);



  plhs[0] = mxCreateDoubleMatrix(nRowsSignal, NumIter, mxREAL);
  Eigen::Map<Eigen::MatrixXd> DiffusedSignal (mxGetPr(plhs[0]), nRowsSignal, NumIter);
  DiffusedSignal = SignalOut;


}
