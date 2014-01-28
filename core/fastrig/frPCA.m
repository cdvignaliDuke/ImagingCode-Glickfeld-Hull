function frPCA(newdir,varargin)
%FRPCA
% FRPCA(NEWDIR)
% FRPCA(...,VARARGIN) where varargin is list of property pairs 
% 

nFramesPerFile = 1000;

fprintf('%s: ',newdir);

defaultopts = {'Overwrite',false,'k',6};
options = cell2struct(defaultopts(2:2:end),defaultopts(1:2:end),2);

for iarg = 1:2:length(varargin)
    options.(varargin{iarg}) =varargin{iarg+1}
end

dirs = frGetDirs;

imdir = fullfile(dirs.registered,newdir);

if exist(imdir) ~= 7
    str = sprintf('directory ''%s'' does not exist',imdir);
    error(str);
end

%% recursively process children directories
list = dir(imdir);list(1:2)=[]; % deletes curent and parent dir entry
sel = find([list.isdir]);

if ~any(strmatchl('green',{list(sel).name}))
    for index = sel
        frPCA(fullfile(newdir,list(index).name),varargin{:});
    end
    return;
end
    
%% path names
[animal,name]=fileparts(newdir);
expt = frGetExpt(animal,name);

%%
greenTargetPn = expt.dirs.pcagreenpn;
[dummy,greenTargetFn] = fileparts(expt.filenames.pca);

list = dir(fullfile(greenTargetPn,'*.tif'));
if length(list) && ~options.Overwrite
    fprintf('ignored. Already exists.\n');
	return;
end

if length(list) options.Overwrite
    fprintf('deleted existing files');
    delete(fullfile(greenTargetPn,'*.tif'));    
    delete(fullfile(greenTargetPn,'*.mat'));
end

fprintf('\n');

list = dir(fullfile(expt.dirs.reggreenpn,'*.tif')); 
%nfiles = length(list);
   
%% pca 
stack = single(readtiff(expt.dirs.reggreenpn));
[ny,nx,nt]=size(stack);
stack = reshape(stack,ny*nx,nt);
av = mean(stack,2);
stack = bsxfun(@minus,stack,av);

fprintf(1,'Calculating principal components\n');
[u,s,v]=pca(stack,options.k);
% U = reshape(u,[ny,nx,options.k]);

%% load shifts
load(fullfile(expt.dirs.reggreenpn,expt.filenames.shifts))
shifts =[];
for ind = 1:length(outputs)
    shifts = [shifts;outputs{ind}(:,3:4)];
end

rs =[];
for ic = 1:options.k
    r1 = corrcoef(shifts(:,1),v(:,ic));
    r2 = corrcoef(shifts(:,2),v(:,ic));
    rs(ic,:)=[r1(1,2) r2(1,2)];
end

ic = any(abs(rs')>.15);
s0 = diag(diag(s).*ic(:));


v2 = tcRemoveArtifacts(bsxfun(@minus,v,mean(v)),bsxfun(@minus,shifts,mean(shifts)));
figure;subplot(211);tcOffsetPlot(v);subplot(212);tcOffsetPlot(v2);

% v2 = tcRemoveArtifacts(v(:,2),shifts(:,2));
% figure;subplot(211);tcOffsetPlot(shifts);subplot(212);tcOffsetPlot(v);
% figure;subplot(211);tcOffsetPlot(tcCycleAverage(shifts,2304));subplot(212);tcOffsetPlot(tcCycleAverage(v,2304))

%%
figure;imagesq(U(:,:,1))
%%

fprintf(1,'Removing the following components ');
fprintf(1,'%i ',find(ic));
fprintf(1,'\n');

nfiles = ceil(nt/nFramesPerFile);

for iFile =1:nfiles    
    sel = (iFile-1)*nFramesPerFile+1:min(iFile*nFramesPerFile,nt);
    nFrames = length(sel);
    
    rec = u*s0*v(sel,:)';
    
    out = stack(:,sel)-rec;
    out = bsxfun(@plus,out,av);
    out = uint16(reshape(out,[ny,nx,nFrames]));
      
    fn = sprintf('%s%06d.tif',greenTargetFn,iFile);
    disp(fn);
    writetiff(out,fullfile(greenTargetPn,fn),'uint16');
end

save(fullfile(greenTargetPn,'usv.mat'),'u','s','v');

return;
