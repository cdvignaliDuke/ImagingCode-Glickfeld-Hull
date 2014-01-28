#include <mex.h>

// [out]=chopvec(in,os,n);
void mexFunction(int nlhs,       mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{

    mwSize frame_dims[]={2,2};
    size_t in_len, os_len, elsiz;
    
    unsigned short i,j;
    
    unsigned short *pin, *pout;
    unsigned long *pos, *plen;

    size_t bytes_to_copy, input_offset,output_offset;
    
    /* Check for proper number of input and output arguments */    
    if (nrhs != 3) {
	mexErrMsgTxt("Three input arguments required.");
    }
    
    if (nlhs > 1){
	mexErrMsgTxt("Too many output arguments.");
    }
    
    if (mxGetClassID(prhs[0])!= mxINT16_CLASS) {
        mexErrMsgTxt("first argument should be int16");
    }

    if (mxGetClassID(prhs[1])!= mxUINT32_CLASS) {
        mexErrMsgTxt("second argument should be uint32");
    }

    if (mxGetClassID(prhs[2])!= mxUINT32_CLASS) {
        mexErrMsgTxt("third argument should be uint32");
    }
    
    /* Get the number of elements in the input argument */
    in_len = mxGetNumberOfElements(prhs[0]);
    os_len = mxGetNumberOfElements(prhs[1]);

    pin = (signed short *) mxGetPr(prhs[0]);    

    pos = (unsigned long *) mxGetPr(prhs[1]);    
    plen = (unsigned long *) mxGetPr(prhs[2]);

    /* Declare output array */    
    frame_dims[0] = plen[0]; // number of data points per line
    frame_dims[1] = (mwSize) os_len; // number of lines
    
    plhs[0] = mxCreateNumericArray(2, frame_dims, mxINT16_CLASS, mxREAL);
    pout = (signed short *) mxGetPr(plhs[0]);
        
    elsiz = mxGetElementSize(plhs[0]);
    bytes_to_copy = frame_dims[0] * elsiz;
      
    for(j=0;j<os_len;j++){
        
        input_offset = pos[j];
        output_offset = j*frame_dims[0];
        
        // memcpy about 2x as fast 
        memcpy(pout + output_offset,pin + input_offset,bytes_to_copy);
/*
        for (i=0;i<bytes_to_copy;i++) {
            *(pout + output_offset + i) = *(pin + input_offset + i);
        }
*/
    }    
    
    return;
}

/*

 nd = 32;
 in = uint16([0:nd-1]);
 os = uint32([2 5 8]);
 n = uint32(3);
 [out]=chopvec(in,os,n);
 
 *
 */
