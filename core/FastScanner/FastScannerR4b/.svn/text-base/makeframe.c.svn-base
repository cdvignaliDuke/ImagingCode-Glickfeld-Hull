#include <mex.h>

// [frame]=makeframe(raw,linestarts,xii,width,linelen,delay)
//delay=maxdelay-delay

void mexFunction(int nlhs,       mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{

    mwSize frame_dims[]={2,2};
    unsigned short i,j;
    signed short *pin, *pout, *pinCol, *poutCol;
    unsigned short *pxii,*pwidth;
    unsigned long *pos, *plen, *del;
    unsigned short count = 0;
    size_t line_n, ind_n, hf_line_n;

    /* Check for proper number of input and output arguments */    
    if (nrhs != 6) 
    {
	mexErrMsgTxt("Six input arguments required.");
    }
    
    if (nlhs > 1)
    {
	mexErrMsgTxt("Too many output arguments.");
    }
    
    if (mxGetClassID(prhs[0])!= mxINT16_CLASS) 
    {
        mexErrMsgTxt("first argument should be int16");
    }

    if (mxGetClassID(prhs[1])!= mxUINT32_CLASS)  
    {
        mexErrMsgTxt("second argument should be uint16");
    }

    if (mxGetClassID(prhs[2])!= mxUINT16_CLASS) 
    {
        mexErrMsgTxt("third argument should be uint16");
    }

    if (mxGetClassID(prhs[3])!= mxUINT16_CLASS) 
    {
        mexErrMsgTxt("Fourth argument should be uint16");
    }
    if (mxGetClassID(prhs[4])!= mxUINT32_CLASS) 
    {
        mexErrMsgTxt("Fourth argument should be uint16");
    }

hf_line_n=mxGetNumberOfElements(prhs[1]);
    line_n = mxGetNumberOfElements(prhs[1])*2;
    ind_n = mxGetNumberOfElements(prhs[2]);
    pin = (signed short *) mxGetPr(prhs[0]);
    pos = (unsigned long *) mxGetPr(prhs[1]);
    pxii = (unsigned short *) mxGetPr(prhs[2]);
    pwidth = (unsigned short *) mxGetPr(prhs[3]);
plen = (unsigned long *) mxGetPr(prhs[4]);
del = (unsigned long *) mxGetPr(prhs[5]);

    /* Declare output array */    
    frame_dims[0] = (mwSize) line_n; 
    frame_dims[1] = *pwidth; 

    plhs[0] = mxCreateNumericArray(2, frame_dims, mxINT16_CLASS, mxREAL);    
    pout = (signed short *) mxGetPr(plhs[0]);

    count = 0;
    // loop over input indx
    for(j=0;j<ind_n-1;j++)
    { 
        //loop overlines

        for(i=0;i<hf_line_n;i++)
        {
            *(pout+(pxii[j]-1)*line_n+i*2) = *(pout+(pxii[j]-1)*line_n+i*2) + *(pin+pos[i]+*(del)+j);
//            *(pout+(*(pwidth)-(pxii[j]-1)-1)*line_n+i*2+1) = *(pout+(*(pwidth)-(pxii[j]-1)-1)*line_n+i*2+1) + *(pin+pos[i]+*(plen)-1+*(del)+j);
            *(pout+(*(pwidth)-pxii[j])*line_n+i*2+1) = *(pout+(*(pwidth)-pxii[j])*line_n+i*2+1) + *(pin+pos[i]+*(plen)-1+*(del)+j);

        }

        count = count + 1;

//        // normalize
        if (j + 1 < ind_n & pxii[j+1]>pxii[j] & count > 0) 
        {
            if (count>1)
            {
                 for(i=0;i<hf_line_n;i++)
                 {
                      *(pout+(pxii[j]-1)*line_n+i*2) = *(pout+(pxii[j]-1)*line_n+i*2) / count;
            		*(pout+(*(pwidth)-(pxii[j]-1)-1)*line_n+i*2+1) = *(pout+(*(pwidth)-(pxii[j]-1)-1)*line_n+i*2+1) / count;
                 }
            }            
             count = 0;            
        }   
    }

    j = ind_n-1;

    if (count>1)
    {
        for(i=0;i<hf_line_n;i++)
        {
           *(pout+(pxii[j]-1)*line_n+i*2) = *(pout+(pxii[j]-1)*line_n+i*2) / count;
           *(pout+(*(pwidth)-(pxii[j]-1)-1)*line_n+i*2+1) = *(pout+(*(pwidth)-(pxii[j]-1)-1)*line_n+i*2+1) / count;
         }
     }  
    return;
}




