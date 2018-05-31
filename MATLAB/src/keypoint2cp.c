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





double *pVertices, *pFaces;
// 
unsigned long numVertices, numFaces;

// int DEBUG_LEVEL;


/*
 * Functions for max/min floating point.
 * (these were inline in C99 but not sure how to do that in ansi C)
 */
myfloat float_max(myfloat a, myfloat b)
{
  return (a > b) ? a : b;
}

myfloat float_min(myfloat a, myfloat b)
{
  return (a < b) ? a : b;
}




/*
 * Project the point (c1,c2,c3) onto the line segment specified by
 * (p1,p2,p3) and (q1,q2,q3).
 * (code by Steve Ruuth)
 */
void ProjectOnSegment(myfloat *c1, myfloat *c2, myfloat *c3, \
		      myfloat p1, myfloat p2, myfloat p3, \
		      myfloat q1, myfloat q2, myfloat q3)
{
  myfloat lambda, lambda_star;
  myfloat cmp1, cmp2, cmp3;
  myfloat qmp1, qmp2, qmp3;

  cmp1 = *c1-p1;
  cmp2 = *c2-p2;
  cmp3 = *c3-p3;
  qmp1 =  q1-p1;
  qmp2 =  q2-p2;
  qmp3 =  q3-p3;

  lambda = (cmp1*qmp1+cmp2*qmp2+cmp3*qmp3)/(qmp1*qmp1+qmp2*qmp2+qmp3*qmp3);
  lambda_star = float_max(FP0, float_min(lambda, FP1));

  *c1 = p1+lambda_star*qmp1;
  *c2 = p2+lambda_star*qmp2;
  *c3 = p3+lambda_star*qmp3;
  
//   return;
}






/*
 * Closest point and distance from a point (a1,a2,a3) to a triangle
 * indexed by face_index.  Returns the *squared* distance and the
 * closest point in (c1,c2,c3).  Uses global vars `face' and
 * `vertex'.  (Based on code by Steve Ruuth)
 * TODO: this code may not be robust to degenerate triangles (line
 * segments and points).  More testing required.
 */
void FindClosestPointToOneTri(myfloat a1, myfloat a2, myfloat a3, \
				 unsigned long face_index, \
				 myfloat *c1, myfloat *c2, myfloat *c3)
{
  unsigned long index_p, index_q, index_r;
  myfloat a11, a12, a22, b1, b2;
  myfloat i11, i12, i22, factor;
  myfloat q1, q2, q3;
  myfloat r1, r2, r3;
  myfloat lambda, mu;
  myfloat dd;

  
  /* obtain the indices to the three vertices */
  index_p = (long) pVertices[ (long) pFaces[face_index] ];
  index_q = (long) pVertices[ (long) pFaces[numFaces + face_index] ]; 
  index_r = (long) pVertices[ (long) pFaces[2*numFaces + face_index] ]; 
  
//     mexPrintf("Here %lu \n", face_index);
  mexPrintf("Here %ld \n", index_p);
  return;

//   index_p = face[face_index].v1;
//   index_q = face[face_index].v2;
//   index_r = face[face_index].v3;

  /* translate so the p is at the origin */
  a1 -= pVertices[ index_p ];
  a2 -= pVertices[ numVertices + index_p ];
  a3 -= pVertices[ 2*numVertices + index_p ];
//   a1 -= vertex[index_p].x;
//   a2 -= vertex[index_p].y;
//   a3 -= vertex[index_p].z;
  
  q1 = pVertices[index_q] - pVertices[ index_p ];
  q2 = pVertices[numVertices + index_q] - pVertices[ numVertices +index_p ];
  q3 = pVertices[2*numVertices + index_q] - pVertices[ 2*numVertices +index_p ];
//   q1 = vertex[index_q].x-vertex[index_p].x;
//   q2 = vertex[index_q].y-vertex[index_p].y;
//   q3 = vertex[index_q].z-vertex[index_p].z;
  
  r1 = pVertices[index_r] - pVertices[index_p];
  r2 = pVertices[numVertices + index_r] - pVertices[numVertices + index_p];
  r3 = pVertices[2*numVertices + index_r] * pVertices[2*numVertices + index_p];  
//   r1 = vertex[index_r].x-vertex[index_p].x;
//   r2 = vertex[index_r].y-vertex[index_p].y;
//   r3 = vertex[index_r].z-vertex[index_p].z;


  
  /* evaluate the various matrix entries */
  a11 = q1*q1+q2*q2+q3*q3;
  a12 = q1*r1+q2*r2+q3*r3;
  a22 = r1*r1+r2*r2+r3*r3;
  b1  = a1*q1+a2*q2+a3*q3;
  b2  = a1*r1+a2*r2+a3*r3;
  


  /* find the inverse matrix and solve for lambda and mu */
  factor = FP1/(a11*a22-a12*a12);
  i11 = a22*factor;
  i12 = -a12*factor;
  i22 = a11*factor;
  lambda = i11*b1+i12*b2;
  mu     = i12*b1+i22*b2;

  *c1 = lambda * q1 + mu * r1;

  *c2 = lambda * q2 + mu * r2;
  *c3 = lambda * q3 + mu * r3;
  


  if ((lambda<0) && (mu<0) && (lambda+mu<=1)) {
    *c1 = *c2 = *c3 = FP0;
  } else if ((lambda>=0) && (mu<0) && (lambda+mu<=1)) {
    ProjectOnSegment(c1,c2,c3,FP0,FP0,FP0,q1,q2,q3);
  } else if ((lambda>=0) && (mu<0) && (lambda+mu>1)) {
    *c1 = q1;
    *c2 = q2;
    *c3 = q3;
  } else if ((lambda>=0) && (mu>=0) && (lambda+mu>1)) {
    ProjectOnSegment(c1,c2,c3,q1,q2,q3,r1,r2,r3);
  } else if ((lambda<0) && (mu>=0) && (lambda+mu>1)) {
    *c1 = r1;
    *c2 = r2;
    *c3 = r3;
  } else if ((lambda<0) && (mu>=0) && (lambda+mu<=1)) {
    ProjectOnSegment(c1,c2,c3,r1,r2,r3,FP0,FP0,FP0);
  } else if ((lambda>=0) && (mu>=0) && (lambda+mu<=1)) {
    /* TODO: do nothing, what case is this? */
  } else {
//     dbg_printf(0, "Error: non-enumerated case, can this happen?\n");
// #ifdef EXTENDEDPRECISION
//     dbg_printf(0, "lambda mu %Lg %Lg\n", lambda, mu);
// #else
//     dbg_printf(0, "lambda mu %g %g\n", lambda, mu);
//     dbg_printf(0, "factor %g\n", factor);
//     dbg_printf(0, "a11,a22,a12=%g,%g,%g\n", a11, a22, a12);
//     dbg_printf(0, "det?=%g\n", a11*a22 - a12*a12);
// #endif
//     myerr("Error in CP to tri, unanticipated case");
  }

  /* Calculate distance */
  /* Note: dd is dist squared! */
  /* TODO: HORRIBLE HACK 2010-07-28, this was to deal with a ply file
     with degenerate triangles. */
  /*if (isinf(factor)) {*/
  /* TODO, just gets worse and worse, where is isinf on windows? */
  if ((factor > 1e20) || (factor < -1e20)) {
    dd = 10000.0;
    myerr("'factor' is infinite: panic!");
  } else {
    dd  = (a1-*c1)*(a1-*c1)+(a2-*c2)*(a2-*c2)+(a3-*c3)*(a3-*c3);
  }


  /* Shift everything back */
  *c1 += pVertices[index_p];
  *c2 += pVertices[numVertices + index_p];
  *c3 += pVertices[2*numVertices + index_p];
//   *c1+= vertex[index_p].x;
//   *c2+= vertex[index_p].y;
//   *c3+= vertex[index_p].z;
//   return (dd);  /* return square distance */
  return;
}











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
    myfloat *pFaceInd;
    
    
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
        
//         mexPrintf("nkp %d\n", nkp);
        
        pFaceCellElem = mxGetCell(prhs[2], nkp);
        
        NumFaces = mxGetNumberOfElements(pFaceCellElem);
//         mexPrintf("num f %d\n", NumFaces);
        
        pFaceInd = mxGetPr(pFaceCellElem);
        
        for (nf=0;nf<NumFaces;nf++){
            
            
            mexPrintf("ind %g\n", pFaceInd[nf]);
            
            curFace = pFaceInd[nf];

            
            x = pKeypoints[nkp];
            y = pKeypoints[numKeypoints + nkp];
            z = pKeypoints[2*numKeypoints + nkp];
            
//             mexPrintf("x %f, %f, %f\n",x,y,z);

            
            FindClosestPointToOneTri(x, y, z, curFace, &cpx, &cpy, &cpz);
            
            
//             mexPrintf("keypoint : %d \n", nkp);
//             mexPrintf("distance^2:  %g \n", dd);
            
        }
        
    }
    
    
    
    return;
}


