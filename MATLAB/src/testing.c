/* ========================================================================
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *========================================================================*/

#include <mex.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// #include <tgmath.h>

typedef double myfloat;


myfloat FP0 = 0.0;
myfloat FP1 = 1.0;
myfloat FP2 = 2.0;
myfloat FP3 = 3.0;
myfloat FP4 = 4.0;
myfloat FP5 = 5.0;
myfloat FPhalf = 0.5;


#define myerr(s) { mexErrMsgTxt(s); exit(EXIT_FAILURE); }





unsigned long *pVertices, *pFaces;

unsigned long numVertices, numFaces;





/*  the gateway routine.  */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[] )
{
    
//     const mxArray *pFaceCell;
//     const mxArray *pFacesArray;
    
    
    double *pVertices, *pFaces, *pKeypoints, *pFaceCell;
    unsigned long int nkp, nf;
    myfloat cpx, cpy, cpz, dd;
    myfloat x, y, z;
    unsigned long curFace;
    unsigned long numKeypoints;
    unsigned long NumFaces;
    const mwSize *dims;
    mxArray *pFaceCellElem;
    double *pFaceInd;
    
    
    mexPrintf("nrhs %d ", nrhs);
    
    pVertices = mxGetPr(prhs[0]);
    pFaces = mxGetPr(prhs[1]);
    pFaceCell = mxGetPr(prhs[2]);
    pKeypoints = mxGetPr(prhs[3]);
    
    
    numVertices = (unsigned long) mxGetM(prhs[0]);
    numFaces = (unsigned long) mxGetM(prhs[1]);
    numKeypoints = (unsigned long) mxGetM(prhs[3]);
    
    mexPrintf("num vert %d\n", numVertices);
    mexPrintf("num face %d\n", numFaces);
    mexPrintf("num kp %d\n", numKeypoints);
    
    dims = mxGetDimensions(prhs[2]);
    mexPrintf("cell dim %d \n", dims[0]);
    
    
    for (nkp=0; nkp<numKeypoints; nkp++){
        
        mexPrintf("nkp %d\n", nkp);
        
        pFaceCellElem = mxGetCell(prhs[2], nkp);
        
//         mxGetUint64s
        NumFaces = mxGetNumberOfElements(pFaceCellElem);
        mexPrintf("num f %d\n", NumFaces);
        
        pFaceInd = mxGetPr(pFaceCellElem);
        
        for (nf=0;nf<NumFaces;nf++){
            
            
            mexPrintf("ind %g\n", pFaceInd[nf]);
            
            curFace =pFaceInd[nf];
            
            x = pKeypoints[nkp];
            y = pKeypoints[numKeypoints + nkp];
            z = pKeypoints[2*numKeypoints + nkp];
            
            mexPrintf("x %f, %f, %f\n",x,y,z);

            
            FindClosestPointToOneTri(x, y, z, curFace, &cpx, &cpy, &cpz, &dd);
            
            
            mexPrintf("keypoint %d: ", nkp);
            mexPrintf("distance^2 %g :", dd);
            
        }
        
    }
    
    
    
    return;
}