function [labelimg, bw3] = imFindCells (in, params)
%IMFINDCELLS Find cells from image
%   [labelimg, bw3] = imFindCells (in, params)
%
% in: input 2-D image
%
% bw3:(logical) binary image of all of the cell regions
% labelimg:(double) same as bw3, but regions are indexed
%
% params: structure with following fields
%
%   java_hist_eq: default false: this is now for debugging
%   min_area: minimum cell area in pixels
%   win_size: length in pixels of sides of square neighborhood used for
%       adaptive histogram. 
%       suggestion: (cell diameter) * 3
%   noise_radius: radius in pixels of image 'speckle'
%       suggestion: 1 for 128x128 or 256x256
%                   2 for 512x512
%   junk_radius: radius of bright 'junk' in between cell bodies
%       suggestion: 1 if cell area is any smaller than 100 pixels
%   ratio_th: percentage of cell area in 2-D image.
%       This sets the threshold for candidate cell areas after hist eq
%       which is expressed as a percentage of the histogram
%       suggestion: 0.1 to 0.2 (10% - 20%)
%   min_center: minimum disk area in pixels that a cell center must encompass
%       suggestion: (1/2 to 3/4) of min_area 
%   con_ratio: minimum concavity ratio {smaller ratio means greater
%       concavity will be allowed in objects} 
%       ratio = (radius of bridge between objects)/(radius of smallest object)
%       in general: higher ratios increase breakup 
%       suggestion: 0.5 to 1, lower values if breakup_factor is high
%   breakup_factor: # of iterations of bottom hat subtraction preceding
%       breakup based on concavity 
%       suggestion: 1 to 5 (greater than 10 will be ignored)
%   fast_breakup: (0 or 1) can decrease breakup processing time, but may
%       be slightly less effective. 
%   clear_border: (0,1,2) eliminates all objects at the border of the image
%       1 means do this before morphological analysis
%       2 means do it after
%       suggestion: 0 unless you have bright cells right at edge, then 2
%   show_flag: (0 or 1) displays images at relevant processing steps 
%
%
% modification of find_cellKO 
%   by SY(addition of postprocessing) 
%   by AK (addition of more morphological analysis) 
%
% 03/03/08 VB Cleaned up. Rescales before uint8 cast. Parameters passed
%   as struct params 
% 18/03/08 VB Histogram equalization performed on uin16
%
%$Id: imFindCellsV1.m 204 2008-04-30 19:07:26Z histed $

%%%%%%% default values %%%%%%%%%%

if ~isfield(params, 'do_hist_eq')
    params.do_hist_eq=1;
end

if ~isfield(params, 'java_hist_eq')
    params.java_hist_eq=0;
end

if ~isfield(params,'win_size')
    params.win_size=16;
end

if ~isfield(params, 'reverse_borders')
    params.reverse_borders=0;
end


in = imScale(in,[],[],'uint16');

si = size(in);
savein = in;

if params.show_flag>0
    figure;
    imshow(in);
    title('Input');
end

%% Perform adaptive histogram equalization 
if params.do_hist_eq
    if params.java_hist_eq
        % use java filter_rank
        if params.reverse_borders==0
            in = imAdaptHistEqJava(in, params.win_size);
        else
            in = imAdaptHistEqJavaRevBorders(in, params.win_size);
        end
    else
        % use MATLAB image processing toolbox
        Ntiles_x=ceil(si(1)./params.win_size);
        Ntiles_y=ceil(si(2)./params.win_size);
        if params.reverse_borders==0
            in = adapthisteq(in, ...
                             'clipLimit',0.02,...
                             'numTiles',[Ntiles_x Ntiles_y]); 
        else
            in = adapthisteqRevBorders(in, Ntiles_x, Ntiles_y, 0.02);
        end
    end
end    

in = imScale(in,[],[],'uint8'); % WARNING CODE BELOW ASSUMES UINT8


%h = fspecial('gaussian',[10 10],1);
%in = imfilter(in,h);

saveadapt = in;

if params.show_flag>0
    figure;
    imshow(in);
    title('Hist Equalization');
end

% clear small particles
clear_radius = params.noise_radius -1;
if clear_radius > 0;
    se = strel('disk', clear_radius);
    in = imopen(in, se);
else
    se = strel('disk', 1);
    in = imopen(in, se);
end

if params.show_flag>0
    figure;
    imshow(in);
    title('Removed Noise');
end

% **** effects marginal ******
% fill holes in cells
se = strel('disk', (params.noise_radius));
in = imclose(in, se);

if params.show_flag>0
    figure;
    imshow(in);
    title('Filled holes');
end

% eliminate medium size 'junk' forming bridges between cells

% REM by AK 9/21/04
%se = strel('disk',round(params.junk_radius));
%nhood = getnhood(se);
%idx = find(nhood == 1);
%filtwin = size(nhood);
%thresh = 80;
%in = nlfilter(in, filtwin, @adapterodeAK, idx, thresh);

%if params.show_flag>0
%    figure;
%    imshow(in);
%    colormap(gray);
%    title('temp');
%end

%in = imdilate(in, se);

se = strel('disk', params.junk_radius);
in = imopen(in, se);

if params.show_flag>0
    figure;
    imshow(in);
    colormap(gray); %axis square;
    title('Junk Removed');
end

% threshold at relatively low value to capture all regions (cellular area)
% even if neighboring cells are connected
sorted=sort(reshape(in,si(1)*si(2),1));
th=sorted(round(si(1)*si(2)*(1-params.ratio_th)),1);
regionmask = (in>=th);

if params.show_flag>0
    figure;
    imshow(regionmask);
    title('Region Mask');
end

% strong clear border objects
if params.clear_border==1;
    regionmask = imclearborder(regionmask,4);
    if params.show_flag>0
        figure;
        imshow(regionmask);
        title('Border Constraint - strong');
    end
end

regionmask = bwareaopen(regionmask, params.min_area);

if params.show_flag>0
    figure;
    imshow(regionmask);
    title('Area and Border Constraints');
end

shrink = immultiply(in,regionmask);

%subtract bottom hat to enhance small dark features
disp('Preforming bottom hat subtraction...please wait');
if (params.breakup_factor>10)
    params.breakup_factor = 10;
end
r = sqrt(params.min_area/pi);
se = strel('disk', round(r));
for i=1:round(params.breakup_factor)
    w = warning('off', 'MATLAB:intMathOverflow');  % neg cast to 0 on subtract
    warning('off', 'MATLAB:intConvertOverflow');  % bothat can raise this warn
    bh = imbothat(shrink,se);
    shrink = imsubtract(shrink,bh);
    warning(w);
end
disp('Bottom hat subtraction complete');

if params.show_flag>0
    figure;
    imshow(shrink);
    title('-BottomHat Gray');
end

%auto-threshold
level = graythresh(shrink);
shrinkmask = im2bw(shrink,level);

if params.show_flag>0
    figure;
    imshow(shrinkmask);
    title('-BottomHat Mask');
end


%remove small bright spots that probably represent small pieces of one cell
%shrinkmask = bwareaopen(shrinkmask, round(params.min_center/2));

% break up bright areas based on concavity and preform no merge dilation
bw3 = imBreakCells(regionmask,shrinkmask,params.min_area, params.fast_breakup, params.con_ratio, params.show_flag);

if params.show_flag>0
    figure;
    imshow(bw3);
    title('Breakup');
end

% reapply area constraint
bw3 = bwareaopen(bw3,params.min_area, 4);

if params.show_flag>0
    figure;
    imshow(bw3);
    title('Total Area Constraints');
end

% remove objects that do not contain the minimum center disk
[templabel num] = bwlabel(bw3);
r2 = sqrt(params.min_center/pi);
se = strel('disk',round(r2));
shrinkmask = imerode(bw3,se);

if params.show_flag>0
    figure;
    imshow(shrinkmask);
    title('Center Area Constraint');
end

for i=1:num
    idx = find(templabel==i);
    if (sum(shrinkmask(idx))==0)
        bw3(idx) = 0;
    end
end

%weak clear border objects
if params.clear_border==2;
    bw3 = imclearborder(bw3,4);
    if params.show_flag>0
        figure;
        imshow(bw3);
        title('Border Constraint - weak');
    end
end

%find labels
labelimg=bwlabel(bw3,4);

%open original image in ImageJ for inspection
if params.do_manual == 1
    wd = java.lang.System.getProperty('user.dir');
    java.lang.System.setProperty('user.dir', 'C:\Program Files\ImageJ');
    IJ = ij.ImageJ;
    bw3ij = subArray2IjAK(savein);
    bw3ij.show;
end

%select objects to be eliminated
if params.do_manual == 1
    shadeimg = imShade(saveadapt, bw3);
    H = figure;
    imshow(shadeimg);
    title('Object Elimination');
    [Y X] = getpts(H);
    for i=1:length(X)
        val = labelimg(round(X(i)),round(Y(i)));
        if val~=0
            idx = find(labelimg==val);
            bw3(idx) = 0;
        end
    end
    close(H);
end


%select positions to add cells

if params.do_manual == 1
    se = strel('disk',1);
    tempmask = logical(zeros(si(1),si(2)));
    r2 = sqrt(params.min_area/pi);
    se2 = strel('disk',round(r2));

    shadeimg = imShadeCells(saveadapt, bw3);
    H = figure;
    imshow(shadeimg);
    title('Object Addition');
    [Y X] = getpts(H);
    for i=1:length(X)
        tempmask(:) = 0;
        dilateorg = imdilate(bw3,se);
        tempmask(round(X(i)),round(Y(i))) = 1;
        tempmask = imdilate(tempmask,se2);
        tempmask = tempmask & ~dilateorg;
        bw3 = bw3 | tempmask;
    end
    close(H);
end


if params.do_manual == 1
%close ImageJ
bw3ij.hide;
ij.WindowManager.closeAllWindows;
IJ.quit;
java.lang.System.setProperty('user.dir', wd);
end

if params.show_flag>0
    figure;
    imshow(bw3);
    title('Output');
end

if params.show_flag>0
    shadeimg = imShade(saveadapt, bw3);
    figure;
    imshow(shadeimg);
    title('Shaded');
    shadeimg = imShade(savein, bw3);
    figure;
    imshow(shadeimg);
    title('Shaded');
    
end

%relabel
labelimg=bwlabel(bw3,4);

%%%%%%%%%%%%%%%%

function out = subArray2IjAK(in)

in = uint8(in);
si = size(in);
ip = ij.process.ByteProcessor(si(2), si(1));
Pix = reshape((in'),(si(2)*si(1)),1);
ip.setPixels(Pix);
out = ij.ImagePlus('',ip);
