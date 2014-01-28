function RegInfo = doCellFind4D(RegInfo)

is_rect=RegInfo.All.is_rect;
cell_dia_pixel=RegInfo.All.cell_dia_um./(RegInfo.All.FOV(1)/RegInfo.All.MATsize(1));
sigma = round(cell_dia_pixel); %halfwidth of gaussian
RegInfo.All.sigma=sigma;

% optargin = size(varargin,2);
% if optargin == 1
%     avg = varargin{1};
% else
%     avg = RegInfo.All.avg_up;
% end    

avg = RegInfo.All.avg_up;
NoDataMask=~(RegInfo.All.mask);
Lx=size(avg,1);

if RegInfo.All.is_rect == 1
    NoDataMask(Lx/2 + [1:Lx/2],:,:) = 1;
end

SKIPLASTSCANLINE=1;
if SKIPLASTSCANLINE == 1
    if is_rect == 1
        NoDataMask(Lx/2,:,:) = 1;
    else
        NoDataMask(Lx,:,:) = 1;
    end
end

stack_avg_norm = contrastnorm_3DMA(avg, sigma, is_rect, NoDataMask);

if ~isfield(RegInfo.All, 'binarymask_up')
    RegInfo.All.binarymask_up=[];
end    
%find cell Mask
Im_contrast = .8;
labelmask_up = imCellEditInteractive_3DMA_090210(stack_avg_norm,RegInfo.All.binarymask_up,[],1,Im_contrast);

binarymask_up=(labelmask_up>0);
%reapply NoDataMask
binarymask_up=binarymask_up.*~NoDataMask;

RegInfo.All.binarymask_up=binarymask_up;
[labelimg_up,nCells] = bwlabeln(binarymask_up,6);

RegInfo.All.labelimg_up=labelimg_up;

labelimg_temp=labelimg_up(:,:,[1:ceil(size(binarymask_up,3)./2)]*2-1);

labelimg=zeros(size(labelimg_temp));

% reindex after removal of some 1-sub-plane cells
counter=1;
for i=1:nCells
    cind=find(labelimg_temp==i);
    if ~isempty(cind)
        [x,y,z]=ind2sub(size(labelimg_temp),cind);
        for j=1:length(x)
            labelimg(x(j),y(j),z(j))=counter;
        end    
        counter=counter+1;
    end
end    

nCells=counter-1;

binarymask=binarymask_up(:,:,[1:ceil(size(binarymask_up,3)./2)]*2-1);

%[labelimg,nCells] = bwlabeln(binarymask,6);

strr=regionprops(labelimg,'Centroid', 'Area');
centroidtmp={strr.Centroid};
Centroid = cell2mat(centroidtmp')';
Area = [strr.Area]; 

RegInfo.Cell.avg_norm = stack_avg_norm;
RegInfo.Cell.binarymask = binarymask;
RegInfo.Cell.labelimg = labelimg;
RegInfo.Cell.nCells = nCells;
RegInfo.Cell.centroid = Centroid;
RegInfo.Cell.area = Area;

alloutfname=RegInfo.All.alloutfname;

RegInfo.Code.doCellFind4D=getCode('doCellFind4D');

save([alloutfname,'RegInfo.mat'],'RegInfo');