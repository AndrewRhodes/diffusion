

#ifndef BDF_HPP_
#define BDF_HPP_

#include <Eigen/Eigen>
#include <Eigen/IterativeLinearSolvers>


Eigen::MatrixXd runBDF1(Eigen::SparseMatrix<double> ItL0, Eigen::VectorXd SignalIn, int numIter){

  int numRows = SignalIn.rows();
  Eigen::MatrixXd SignalOut(numRows, numIter);
  SignalOut.col(0) = SignalIn;

  Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;
  solver.compute(ItL0);

  for (int j = 0; j < numIter - 1; j++ ){
    SignalOut.col(j+1) = solver.solve(SignalOut.col(j));
  }

  return SignalOut;
}




Eigen::MatrixXd runBDF2(Eigen::SparseMatrix<double> ItL0, Eigen::SparseMatrix<double> ItL1, Eigen::VectorXd SignalIn, int numIter){

  int numRows = SignalIn.rows();
  Eigen::MatrixXd SignalOut(numRows, numIter);
  SignalOut.col(0) = SignalIn;

  Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;


  for (int j = 0; j < numIter - 1; j++ ){
    if (j==0){
      solver.compute(ItL0);
      SignalOut.col(j+1) = solver.solve(SignalOut.col(j));
    } else if(j==1){
      solver.compute(ItL1);
      SignalOut.col(j+1) = solver.solve((4/3)*SignalOut.col(j) + (1/3)*SignalOut.col(j-1));
    } else {
      SignalOut.col(j+1) = solver.solve((4/3)*SignalOut.col(j) + (1/3)*SignalOut.col(j-1));
    }

  }

  return SignalOut;

}


Eigen::MatrixXd runBDF3(Eigen::SparseMatrix<double> ItL0, Eigen::SparseMatrix<double> ItL1, Eigen::SparseMatrix<double> ItL2, Eigen::VectorXd SignalIn, int numIter){

  int numRows = SignalIn.rows();
  Eigen::MatrixXd SignalOut(numRows, numIter);
  SignalOut.col(0) = SignalIn;

  Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;


  for (int j = 0; j < numIter - 1; j++ ){
    if (j==0){
      solver.compute(ItL0);
      SignalOut.col(j+1) = solver.solve(SignalOut.col(j));
    } else if (j==1){
      solver.compute(ItL1);
      SignalOut.col(j+1) = solver.solve((4/3)*SignalOut.col(j) + (1/3)*SignalOut.col(j-1));
    } else if (j==2){
      solver.compute(ItL2);
      SignalOut.col(j+1) = solver.solve((18/11)*SignalOut.col(j) - (9/11)*SignalOut.col(j-1) + (2/11)*SignalOut.col(j-2));
    } else {
      SignalOut.col(j+1) = solver.solve((18/11)*SignalOut.col(j) - (9/11)*SignalOut.col(j-1) + (2/11)*SignalOut.col(j-2));
    }

  }

  return SignalOut;


}



Eigen::MatrixXd runBDF4(Eigen::SparseMatrix<double> ItL0, Eigen::SparseMatrix<double> ItL1, Eigen::SparseMatrix<double> ItL2, Eigen::SparseMatrix<double> ItL3, Eigen::VectorXd SignalIn, int numIter){

  int numRows = SignalIn.rows();
  Eigen::MatrixXd SignalOut(numRows, numIter);
  SignalOut.col(0) = SignalIn;

  Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;

  for (int j = 0; j < numIter - 1; j++ ){
    if (j==0){
      solver.compute(ItL0);
      SignalOut.col(j+1) = solver.solve(SignalOut.col(j));
    } else if (j==1){
      solver.compute(ItL1);
      SignalOut.col(j+1) = solver.solve((4/3)*SignalOut.col(j) + (1/3)*SignalOut.col(j-1));
    } else if (j==2){
      solver.compute(ItL2);
      SignalOut.col(j+1) = solver.solve((18/11)*SignalOut.col(j) - (9/11)*SignalOut.col(j-1) + (2/11)*SignalOut.col(j-2));
    } else if (j==3){
      solver.compute(ItL3);
      SignalOut.col(j+1) = solver.solve((48/25)*SignalOut.col(j) - (36/25)*SignalOut.col(j-1) + (16/25)*SignalOut.col(j-2) - (3/25)*SignalOut.col(j-3));
    } else {
      SignalOut.col(j+1) = solver.solve((48/25)*SignalOut.col(j) - (36/25)*SignalOut.col(j-1) + (16/25)*SignalOut.col(j-2) - (3/25)*SignalOut.col(j-3));
    }

  }

  return SignalOut;



}





typedef Eigen::SparseMatrix<double,Eigen::ColMajor,std::make_signed<mwIndex>::type> MatlabSparse;


Eigen::Map<MatlabSparse >
matlab_to_eigen_sparse(const mxArray * mat)
{
    mxAssert(mxGetClassID(mat) == mxDOUBLE_CLASS,
    "Type of the input matrix isn't double");
    mwSize     m = mxGetM (mat);
    mwSize     n = mxGetN (mat);
    mwSize    nz = mxGetNzmax (mat);
    /*Theoretically fails in very very large matrices*/
    mxAssert(nz <= std::numeric_limits< std::make_signed<mwIndex>::type>::max(),
    "Unsupported Data size."
    );
    double  * pr = mxGetPr (mat);
    MatlabSparse::StorageIndex* ir = reinterpret_cast<MatlabSparse::StorageIndex*>(mxGetIr (mat));
    MatlabSparse::StorageIndex* jc = reinterpret_cast<MatlabSparse::StorageIndex*>(mxGetJc (mat));
    Eigen::Map<MatlabSparse> result (m, n, nz, jc, ir, pr);
    return result;
}

mxArray*
eigen_to_matlab_sparse(const Eigen::Ref<const MatlabSparse,Eigen::StandardCompressedFormat>& mat)
{
    mxArray * result = mxCreateSparse (mat.rows(), mat.cols(), mat.nonZeros(), mxREAL);
    const MatlabSparse::StorageIndex* ir = mat.innerIndexPtr();
    const MatlabSparse::StorageIndex* jc = mat.outerIndexPtr();
    const double* pr = mat.valuePtr();

    mwIndex * ir2 = mxGetIr (result);
    mwIndex * jc2 = mxGetJc (result);
    double  * pr2 = mxGetPr (result);

    for (mwIndex i = 0; i < mat.nonZeros(); i++) {
        pr2[i] = pr[i];
        ir2[i] = ir[i];
    }
    for (mwIndex i = 0; i < mat.cols() + 1; i++) {
        jc2[i] = jc[i];
    }
    return result;
}




#endif // BDF_HPP_
