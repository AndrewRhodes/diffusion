


#ifndef EXPM_HPP_
#define EXPM_HPP_

#include <Eigen/Eigen>
#include <unsupported/Eigen/MatrixFunctions>


Eigen::MatrixXd runExpM(Eigen::SparseMatrix<double> LBM, Eigen::VectorXd SignalIn, int numIter){

  int numRows = SignalIn.rows();
  Eigen::MatrixXd SignalOut(numRows, numIter);
  SignalOut.col(0) = SignalIn;

  Eigen::MatrixXd dLBM;
  Eigen::MatrixXd expLBM;

  dLBM = Eigen::MatrixXd(LBM);
  expLBM = dLBM.exp();


//  for (int j = 0; j < numIter - 1; j++ ){
//   SignalOut.col(j+1) =  expLBM * SignalOut.col(j);
//  }

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



#endif // EXPM_HPP_
