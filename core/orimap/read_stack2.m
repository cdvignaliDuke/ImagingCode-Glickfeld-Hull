function stack = read_stack2 (fname, channel, nframes, nbinning, depth_correct, slice_thickness)
% read tif files (8-bit) into 3 dimensional matrix.
% pixel value will be between 0 and 1.
% fname: e.g. 'C:\share\OG031210\OG031210_6\OG031210_6_vstim12_t'
% nframes: number of frames to read. Zero to nframes-1 will be read.
% nbinning: binning factor for x-y dimension
% this program requires numextCR & add2d2dCR to run. 
% there are four formats to use this function:
%   read_stack(fname): read only 0th frame
%   read_stack(fname, nframes): read nframes without binning or depth correction
%   read_stack(fname, nframes, nbinning): read nframes with binning without depth correction
%   read_stack(fname, nframes, nbinning, depth_correct (,slice_thickness)): read nframes with binning and depth correction (
% depth_correct: correct attenuation of signal intensity in a deep region
%   assuming attenuation is exp(-depth(with a unit of slice_thichness)/depth_correct)
%   if slice thickness is 1 micron (default value) & set depth_correct as 100, then attenuation is exp(-depth(micron)/100)
%   if slice thickness is not 1 micron, you can specify the slice thickness. if you do not specify, slice thickness = 1 micron.
%   if you take image from deep to shallow, then you can set negative depth_correct value
%   if you use this option, pixel values would be not necessarily between 0 and 1



if nargin <= 1
    nframes = 1;
end;
if nargin <=2
    channel=0;
end;
if nargin <= 3
    nbinning=1;
end;
if nargin <= 4
    depth_correct=0;
end;
if nargin <=5
    slice_thickness =1;
end;

    
filename=[fname,numextCR(0,3),'_ch',numextCR(channel,2),'.tif'];
tmp=double(imread(filename,'tif'))/255;
dim=size(tmp);
stack=zeros(dim(1)/nbinning, dim(2)/nbinning, nframes);
if nbinning > 1
    stack(:,:,1)=add2d2dCR(tmp,nbinning)/(nbinning^2);
else
    stack(:,:,1)=tmp;
end;
if nframes > 1
    for fr = 1:nframes-1
        filename=[fname,numextCR(fr,3),'_ch00.tif'];
        tmp=double(imread(filename))/255;
        if nbinning > 1
            tmp=add2d2dCR(tmp,nbinning)/(nbinning^2);
        end;
        if depth_correct ~= 0
            tmp=tmp.*exp(fr*slice_thickness/depth_correct);
        end;         
        stack(:,:,fr+1)=tmp;            
    end;
end;