#include <string.h> /* Needed for the ceil() prototype. */
#include <stdio.h>
#include <math.h>
#include "mex.h"

/* If you are using a compiler that equates NaN to be zero, you 
 *  * must compile this example using the flag  -DNAN_EQUALS_ZERO.
 *   * For example:
 *    *
 *     *     mex -DNAN_EQUALS_ZERO fulltosparse.c
 *      *
 *       * This will correctly define the IsNonZero macro for your C
 *        * compiler.
 *         */

#if defined(NAN_EQUALS_ZERO)
#define IsNonZero(d) ((d) != 0.0 || mxIsNaN(d))
#else
#define IsNonZero(d) ((d) != 0.0)
#endif

void mexFunction(
		int nlhs,       mxArray *plhs[],
		int nrhs, const mxArray *prhs[]
		)
{
	if (nrhs != 3) {
		mexErrMsgTxt("Three input argument required.");
	} 
	if (nlhs > 1) {
		mexErrMsgTxt("Too many output arguments.");
	}
	
	int n = mxGetN(prhs[0]);
	int m = mxGetM(prhs[0]); 
	int l = mxGetM(prhs[1]); 
	double* x = mxGetPr(prhs[0]);
	double* c = mxGetPr(prhs[1]);
	double* theta = mxGetPr(prhs[2]);
	double alpha = theta[1]; 
	double gamma = 0.5*theta[0];
	plhs[0] = mxCreateDoubleMatrix(m, l, mxREAL); 
	double* out = mxGetPr(plhs[0]); 
	
	int i, j, k; 

	for (i = 0; i < m; i++) {
		for (j = 0; j < l; j++) {
			double dist = 0.0; 
			for (k = 0; k < n; k++) {
				dist += (x[k*m+i] - c[k*l+j])*(x[k*m+i] - c[k*l+j]); 
			}
			out[j*m+i] = alpha*exp(-dist*gamma); 
		}
	}
}
