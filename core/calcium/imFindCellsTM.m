function [labelimg, hp_img] = imFindCellsTM (img, template_diam, r_threshold, cell_diam, finemode, height_th)
%IMFINDCELLSTM Template matching based cell segmentation
% [LABELIMG, HP_IMG] = IMFINDCELLSTM (IMG, TEMPLATE_DIAM, R_THRESHOLD, CELL_DIAM) where
%  IMG
%  TEMPLATE_DIAM is diamter of difference of Gaussians in pixels and
%  R_THRESHOLD cell detection threshold as correlation coefficient (default 0.5)
%  CELL_DIAM is diameter used for dilation. If it is not specified,
%       CELL_DIAM = TEMPLATE_DIAM
%  finemode: if it is 1, imTemplateMatch will be used instead of
%       normxcorr2. It will be slower.
%  height_th: if this threshold is given, only cells with pixel values higher than the threshold survive.
%        Give the threshold as a top percentile (recommendation: 20)
% by Kenichi Ohki  04.28.2008
% 

% Changelog
% - 08/10/28 vb added template_diam and modified cell_diam

if nargin<3
    r_threshold=0.5;
end

if nargin < 4
    cell_diam = template_diam;
end

if nargin < 5
    finemode=0;
end

if nargin<6
    height_th=[];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% high pass filtering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   
hp_filter=fspecial('gaussian',ceil(template_diam*5),template_diam);

hp_img=img./imFilter2(hp_filter,img);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% template matching with DOG functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

display('template matching....');

template=DOG(template_diam/2);

if finemode==1
    mask=makeDisk(template_diam,size(template,1));
    corr_map=imTemplateMatch(hp_img,template,mask);
else
    corr_map=normxcorr2sm(hp_img,template);
end

%figure;imagesq(img);
%figure;imagesq(hp_img);colorbar
%figure;imagesq(corr_map)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% detect cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cells=(corr_map>r_threshold);
cells=bwlabel(cells);
fprintf('%i cells detected.\n',max(cells(:)));

%%%%%%%%%%%%%%%%%
% dilate cells
%%%%%%%%%%%%%%%%%

labelimg=imDilate2(cells,cell_diam);

%%%%%%%%%%%%%%%%%%%%
% height threshold
%%%%%%%%%%%%%%%%%%%%

if isempty(height_th)
    return;
end

th_value=getPercentile(hp_img,100-height_th);

pix_values=getCellBasedIndexFromPixelMap(hp_img, labelimg);

labelimg=deleteLabels(labelimg, find(pix_values<=th_value));

return;


function labelimg=deleteLabels(labelimg, toBeDeleted)
Ncells=max(labelimg(:));
count=1;
for i=1:Ncells
    if isempty(find(toBeDeleted==i))
        labelimg(find(labelimg==i))=count;
        count=count+1;
    else
        labelimg(find(labelimg==i))=0;
    end
end
return;


function ker = DOG(r)
r1=r/1.9227;         % zero-crossing at r
r2=r1*2;

dim=ceil(r2*5);
ker=fspecial('gaussian',dim,r1)-fspecial('gaussian',dim,r2);
return;

function out = makeDisk(radius, dim);
center=floor(dim/2)+1;
out=zeros(dim);
for x=1:dim
    for y=1:dim
        if norm([x,y]-[center,center])<radius
            out(x,y)=1;
        end
    end
end
return;

function out = normxcorr2sm (in, template)
% same as normxcorr2    
% but select valid pixels (as in conv2)

    imsiz = size(in);
    corr_map_full = normxcorr2(template,in);
    lb = fix((size(corr_map_full)-size(in))/2+1);
    ub = imsiz + lb - 1;
    out = corr_map_full(lb(1):ub(1),lb(2):ub(2));
    return;
    

function out=imTemplateMatch(in, template, mask);

% find features similar to a "template" in input image
% template matching is performed in a mask (aperture)
% output is correlation coefficient map
% 
% Kenichi Ohki 04.28.2008

dim=size(in);
dim2=size(template);
mx=floor(dim2(1)/2);
px=dim2(1)-mx-1;
my=floor(dim2(2)/2);
py=dim2(2)-my-1;

if nargin<3
    mask=ones(dim2(1),dim2(2));
end

out=zeros(dim(1),dim(2));

for x=1:dim(1)
    for y=1:dim(2)
        window_x=[max(1,x-mx):min(dim(1),x+px)];
        window_y=[max(1,y-my):min(dim(2),y+py)];
        window_x2=window_x-x+mx+1;
        window_y2=window_y-y+my+1;
        a=in(window_x,window_y);
        b=template(window_x2,window_y2);
        k=find(mask(window_x2,window_y2));
        a=a(k);
        b=b(k);
        r=corrcoef(a,b);
        out(x,y)=r(1,2);
    end
end
return;


        
       
        
            


