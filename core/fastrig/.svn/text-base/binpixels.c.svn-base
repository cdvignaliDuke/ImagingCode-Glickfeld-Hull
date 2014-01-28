#include <mex.h>

// [frame]=resampleframe(raw,xii,width)

void mexFunction(int nlhs,       mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{

    mwSize *in_dims, new_dims[2]={1,1};
    size_t in_len, os_len, elsiz;
    
    unsigned short i,j;
    
    signed short *pin, *pinCol;

    signed short *pout,*poutCol;
    unsigned short *pxii,*pwidth;
    unsigned short count = 0, last_xii = 0; 
            
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

    if (mxGetClassID(prhs[1])!= mxUINT16_CLASS) {
        mexErrMsgTxt("second argument should be uint16");
    }

    if (mxGetClassID(prhs[2])!= mxUINT16_CLASS) {
        mexErrMsgTxt("third argument should be uint16");
    }
    
    in_dims =mxGetDimensions(prhs[0]);    
    pin = (signed short *) mxGetPr(prhs[0]);          
    pxii = (unsigned short *) mxGetPr(prhs[1]);    

    pwidth = (unsigned short *) mxGetPr(prhs[2]);    
    
    new_dims[0]=in_dims[0];
    new_dims[1]=*pwidth;
      
    plhs[0] = mxCreateNumericArray(2, new_dims, mxINT16_CLASS, mxREAL);    
    pout = (signed short *) mxGetPr(plhs[0]);

    
    count = 0;
    for(j=0;j<in_dims[1];j++){ // loop over input columns
	pinCol = pin+j*new_dims[0];
	poutCol = pout +(pxii[j]-1)*new_dims[0];
        
        // sum next column
        for(i=0;i<in_dims[0];i++){ // loop over rows
            //printf("%i %i\n", input_offset, output_offset);

		*(poutCol) = *(poutCol) + *(pinCol);
		poutCol++;
		pinCol++;
        }
        
        count = count + 1;
        
        // normalize
        if (j + 1 < in_dims[1] & pxii[j+1]>pxii[j] & count > 0) {
if (count>1){
	     poutCol = pout +(pxii[j]-1)*new_dims[0];
             for(i=0;i<in_dims[0];i++){
                  *(poutCol) = *(poutCol) / count;
                  poutCol++;
             }
}            
             count = 0;            
        }            
    }

    j = in_dims[1]-1;
    poutCol = pout +(pxii[j]-1)*new_dims[0];
if (count>1){
    for(i=0;i<in_dims[0];i++){

//        output_offset = (pxii[j]-1)*new_dims[0]+i;
//        *(pout + output_offset) = *(pout + output_offset) / count;
                  *(poutCol) = *(poutCol) / count;
                  poutCol++;
    }
}
    
    return;
}

/*
*/
