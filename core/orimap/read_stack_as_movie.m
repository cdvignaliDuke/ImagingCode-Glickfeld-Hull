function mov = read_stack_as_movie (fname, channel, nframes, color, nbinning, depth_correct, slice_thickness)
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
    channel=0;
end;
if nargin <=2
    nframes=1;
end;
if nargin <=3
    color=[0,0,0];
end;

if nargin <= 4
    nbinning=1;
end;
if nargin <= 5
    depth_correct=0;
end;
if nargin <=6
    slice_thickness =1;
end;
    
filename=[fname,numextCR(0,3),'_ch',numextCR(channel,2),'.tif'];
tmp=double(imread(filename,'tif'))/255;
dim=size(tmp);
if nbinning > 1
    tmp=add2d2dCR(tmp,nbinning)/(nbinning^2);
end

tmp2=zeros(size(tmp,1), size(tmp,2), 3);
tmp2(:,:,1)=tmp.*color(1);
tmp2(:,:,2)=tmp.*color(2);
tmp2(:,:,3)=tmp.*color(3);

mov=im2frame(tmp2);

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
        tmp2(:,:,1)=tmp.*color(1);
        tmp2(:,:,2)=tmp.*color(2);
        tmp2(:,:,3)=tmp.*color(3);

        mov=[mov,im2frame(tmp2)];            
    end;
end;