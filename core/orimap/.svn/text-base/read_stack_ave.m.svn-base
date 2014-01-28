function [stack, xsize, ysize]= read_stack_ave (fname, nframes_per_run, run_inds, nbinning)
% modified from Clay's og0312110CR4.m
% read tif files into 3 dimensional matrix.
% tif images are averaged across multiple runs.
% pixel value will be between 0 and 1.
% fname: e.g. 'C:\share\OG031210\OG031210_6\OG031210_6_vstim12_t'
% nframes_per_run: number of frames in one run (one repetition). 
% run_inds: indexes of runs to be averaged. zero based. e.g. [0:9], [0:3,5:9], or [0,3,5,8]
% nbinning: binning factor for x-y dimension
% this program requires numextCR & add2d2dCR to run. 
% there are three formats to use this function:
%   read_stack(fname, nframes_per_run): read images in the first run, without averaging or binning
%   read_stack(fname, nframes_per_run, run_inds): average across runs index by run_inds, without binning
%   read_stack(fname, nframes_per_run, run_inds, nbinning): average across runs index by run_inds, with binning

if nargin <= 2
    run_inds=[0];
end;
if nargin <= 3
    nbinning=1;
end;
    
filename=[fname,numextCR(0,3),'_ch00.tif'];
tmp=imread(filename,'tif');
dim=size(tmp);
xsize=dim(1)/nbinning;
ysize=dim(2)/nbinning;

stack=zeros(xsize, ysize, nframes_per_run);

for i = 0:nframes_per_run-1
  for j = run_inds
    fileind=j*nframes_per_run + i
    filename=[fname,numextCR(fileind,3),'_ch00.tif'];
    tmp=double(imread(filename))/255;
    if nbinning > 1 
        tmp=add2d2dCR(tmp,nbinning)/(nbinning^2);
    end
    stack(:,:,i+1)=stack(:,:,i+1)+tmp;
  end
end
stack = stack./length(run_inds) ;
