function RegInfo=doNPsubPrep4D(RegInfo,cells_per_thread,multicoredir)

cell_dia_um = RegInfo.All.cell_dia_um;
np_radius_um = RegInfo.All.np_radius_um;
np_radius_pixel=np_radius_um./(RegInfo.All.FOV(1)/RegInfo.All.MATsize(1));
cell_dia_pixel=cell_dia_um./(RegInfo.All.FOV(1)/RegInfo.All.MATsize(1));
z_skew=(RegInfo.All.FOV(3)/RegInfo.All.MATsize(3))/(RegInfo.All.FOV(1)/RegInfo.All.MATsize(1));

np_radius = round(np_radius_pixel);
%erode = ceil(0.1*cell_dia_pixel);
space = ceil(0.5*cell_dia_pixel);

RegInfo.All.np_radius=np_radius;
RegInfo.All.z_skew=z_skew;
%RunInfo.erode=erode;
RegInfo.All.space=space; 

binarymask = RegInfo.Cell.binarymask;
centroid = RegInfo.Cell.centroid;
shiftmask=RegInfo.All.mask;
shiftmask=shiftmask(:,:,[1:ceil(size(shiftmask,3)./2)]*2-1);

nCells = size(centroid,2)
nChunks = ceil(nCells/cells_per_thread)

for i=1:nChunks
        start = (i-1)*(cells_per_thread)+1;
        if i == nChunks
            stop = nCells;
        else
            stop = start + (cells_per_thread) - 1;
        end
        cell_ind=start:stop;
        parameterCell{1,i} = {binarymask,cell_ind,shiftmask,centroid,z_skew,np_radius,space};
end

maxMasterEvaluations = 0;
out = startmulticoremaster(@np_mask4D, parameterCell, multicoredir, maxMasterEvaluations);

for i=1:(nChunks-1)
    npmasks(:,:,:,:,i)=out{1,i};
end

si=size(npmasks);
npmasks=reshape(npmasks,si(1),si(2),si(3),si(4)*si(5));
si=size(npmasks);
si2=size(out{1,nChunks});
if length(si2)>3
    fvols=si2(4);
else
    fvols=1;
end    
npmasks(:,:,:,(si(4)+1):(si(4)+fvols))=out{1,nChunks};

RegInfo.Cell.npmasks=npmasks;

RegInfo.Code.doNPsubPrep4D=getCode('doNPsubPrep4D');

save([RegInfo.All.alloutfname,'RegInfo.mat'],'RegInfo');


%   % this is a hardcoded name of a new directory that is created in rootroot...
%   % perhaps it can be passed as a varargin argument later....
%   procroot=[  rootroot,'\dataproc'];
%   procdir = [procroot,'\',RunInfo.expdir];
%   
%   load([RunInfo.procoutfname{1},'_cells.mat']);  %get cell definitions and info
% 
%   strcat(procdir,'\',RunInfo.fname{1},'_cells.mat')
% 
%   nhood = ones(erode+1,erode+1);
%   if exist('old_binarymask')
%       binarymask = imerode(old_binarymask,nhood);
%   else
%       old_binarymask = binarymask;
%       binarymask = imerode(binarymask,nhood);
%   end
%   old_labelimg = labelimg;
%   labelimg=bwlabel(binarymask,4);
%   for i=1:max(max(old_labelimg))
%     templabel = labelimg(:);
%     cind=find(old_labelimg(:)==i);
%     newnum = unique(templabel(cind));
%     h = find(newnum>0);
%     finalnum=newnum(h);
%     if length(finalnum)==1
%         cellremap(i)=finalnum;
%     else
%         cellremap(i)=0;
%     end
%   end  
%   [strr,centroid,diam,eccentricity,extent,majoraxis,minoraxis] = get_region_statsCR(labelimg);
% 
% 
%   npmasks = np_maskAK(old_binarymask,centroid,np_radius, space);
% 
%  save(strcat(procdir,'\',RunInfo.fname{1},'_cells.mat'),'labelimg','binarymask','old_binarymask','npmasks','strr','centroid','diam','eccentricity','extent','majoraxis','minoraxis','cellremap');  
%  