function stack = read_stack_2 (fname, frame_ind, nbinning) 
% read tif files (8-bit) into 3 dimensional matrix.
% pixel value will be between 0 and 1.
% fname: e.g. 'C:\share\OG031210\OG031210_6\OG031210_6_vstim12_t'
% frame_ind: e.g. [0:799] or [0:319,400:799] etc.
% nbinning: binning factor for x-y dimension
% this program requires numextCR & add2d2dCR to run. 
% there are four formats to use this function:
%   read_stack(fname, nframes): read nframes without binning or depth correction
%   read_stack(fname, nframes, nbinning): read nframes with binning without depth correction


if nargin <= 2
    nbinning=1;
end;

    
filename=[fname,numextCR(frame_ind(1),3),'_ch01.tif'];
tmp=double(imread(filename,'tif'))/255;
dim=size(tmp);
stack=zeros(dim(1)/nbinning, dim(2)/nbinning, size(frame_ind,1));
if nbinning > 1
    stack(:,:,1)=add2d2dCR(tmp,nbinning)/(nbinning^2);
else
    stack(:,:,1)=tmp;
end;
n=1;
for fr = frame_ind
        filename=[fname,numextCR(fr,3),'_ch01.tif'];
        tmp=double(imread(filename))/255;
        if nbinning > 1
            tmp=add2d2dCR(tmp,nbinning)/(nbinning^2);
        end;
        stack(:,:,n)=tmp;
        n=n+1;
end;
