#include <mex.h>
#include <math.h>
#include <vector>
#include <memory.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <assert.h>
#include <Eigen/Eigen>
#include "meshGausMatrix.h"

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    /*  check for proper number of arguments */
	/* NOTE: You do not need an else statement when using mexErrMsgTxt
	 *       within an if statement, because it will never get to the else
	 *       statement if mexErrMsgTxt is executed. (mexErrMsgTxt breaks you out of
	 *       the MEX-file)
	*/
	if(nrhs < 1){
		mexErrMsgTxt("At least one inputs required.");
	}
	if(nrhs > 2){
		mexErrMsgTxt("At most two inputs required.");
	}
	// if(nlhs != 4){
	// 	mexErrMsgTxt("Four output required.");
	// }


    /* check to make sure the first input argument is a string */
	if( mxGetClassID(prhs[0]) != mxCHAR_CLASS){
		mexErrMsgTxt("Input x must be a string.");
	}
	char* filename;
	filename = mxArrayToString(prhs[0]);

  unsigned int htype = 0, dtype = 0;
	bool atype = false;
  double hs = 1.6, rho = 3.0;

  if(nrhs == 2){
    const mxArray *opt = prhs[1];
    mxClassID category = mxGetClassID(opt);

    if(category != mxSTRUCT_CLASS){
			mexErrMsgTxt("Opt input must be a structure");
		}

    // mwSize m_elements;
    // m_elements = mxGetNumberOfElements(opt);
    // mexPrintf("m_elements %d\n", m_elements);

    // mwIndex i;
    int n_field, i_field, ht, dt;
    const char *field_name;
    const mxArray *pfield_array;

    n_field = mxGetNumberOfFields(opt);
    // mexPrintf("n_field %d\n", n_field);

    for(mwIndex i = 0; i < n_field; i++){
      field_name = mxGetFieldNameByNumber(opt, i);
      pfield_array = mxGetFieldByNumber(opt, 0, i);

      if(strcmp(field_name, "hs") == 0){
        mxClassID cat = mxGetClassID(pfield_array);
        hs = (double)(*mxGetPr(pfield_array));
      }
      else if(strcmp(field_name, "rho") == 0){
        rho = (double)(*mxGetPr(pfield_array));
      }
      else if(strcmp(field_name, "htype") == 0){
        htype = (int)(*mxGetPr(pfield_array));
      }
      else if(strcmp(field_name, "dtype") == 0){
        dtype = (int)(*mxGetPr(pfield_array));
      }
			else if(strcmp(field_name, "atype") == 0){
				atype = (bool)(*mxGetPr(pfield_array));
			}

    }

  }

  mexPrintf("dtype: %d, htype: %d, hs: %.2f, rho: %.2f\n", dtype, htype, hs, rho);

	double *I, *J, *S, *ph;
	double h;
  vector<double> Sv, Av;
  vector<unsigned int> Iv, Jv;
	unsigned int nv, nelem;

	if(atype){ // Use area weights
		if(dtype == 0){// dtype == 'euclidean'
			generate_mesh_gaus_matrix(filename, htype, hs, rho, h, Iv, Jv, Sv, Av, nv);
		}
		else if(dtype == 1){//dtype == 'geodesic'
			generate_mesh_gaus_geod_matrix(filename, htype, hs, rho, h, Iv, Jv, Sv, Av, nv);
		}
	}
	else{ // Do not use area weights
	  if(dtype == 0){// dtype == 'euclidean'
	    generate_mesh_gaus_matrix(filename, htype, hs, rho, h, Iv, Jv, Sv, nv);
	  }
	  else if(dtype == 1){//dtype == 'geodesic'
	    generate_mesh_gaus_geod_matrix(filename, htype, hs, rho, h, Iv, Jv, Sv, nv);
	  }
	}

	nelem = Iv.size();

	// mexPrintf("nv: %d \n", nv);

	if(atype){
		plhs[0] = mxCreateDoubleMatrix(nelem, 1, mxREAL);
		plhs[1] = mxCreateDoubleMatrix(nelem, 1, mxREAL);
		plhs[2] = mxCreateDoubleMatrix(nelem, 1, mxREAL);
		plhs[3] = mxCreateDoubleMatrix(nv, 1, mxREAL);
		plhs[4] = mxCreateDoubleMatrix(1, 1, mxREAL);

		double *I, *J, *S, *A, *ph;
		I = mxGetPr(plhs[0]);
		J = mxGetPr(plhs[1]);
		S = mxGetPr(plhs[2]);
		A = mxGetPr(plhs[3]);
		ph = mxGetPr(plhs[4]);

		*ph = h;

		for(mwSize i = 0; i < nelem; i++){
			I[i] = Iv[i];
			J[i] = Jv[i];
			S[i] = Sv[i];
		}

		for(mwSize i = 0; i < nv; i++){
			A[i] = Av[i];
		}
	}
	else{
		plhs[0] = mxCreateDoubleMatrix(nelem, 1, mxREAL);
		plhs[1] = mxCreateDoubleMatrix(nelem, 1, mxREAL);
		plhs[2] = mxCreateDoubleMatrix(nelem, 1, mxREAL);
		plhs[3] = mxCreateDoubleMatrix(1, 1, mxREAL);

	//
		I = mxGetPr(plhs[0]);
		J = mxGetPr(plhs[1]);
		S = mxGetPr(plhs[2]);
		ph = mxGetPr(plhs[3]);

		*ph = h;

		for(mwSize i = 0; i < nelem; i++){
			I[i] = Iv[i];
			J[i] = Jv[i];
			S[i] = Sv[i];
		}
	}

	mexPrintf("h: %f\n", h);

  mxFree(filename);
}
