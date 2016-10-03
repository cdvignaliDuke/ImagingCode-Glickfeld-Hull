awFSAVdatasets_naive100ms
iexp = 5;
% step 1: align drifting grating dataset
fname = fullfile('Z:\data',expt(iexp).mouse,expt(iexp).folder,expt(iexp).date,expt(iexp).dirtuning,[expt(iexp).dirtuning '_000_000']);
sbxalignmaster(fname);

% step 2: select ROIs in drifting grating dataset using scanbox tool.
sbxsegmentflood

%step 3: select ROIs in drifting grating dataset using dF/F and 
% analysis folder
try
    filedir = fullfile('Z:\analysis\',expt(iexp).mouse,expt(iexp).folder,expt(iexp).date,expt(iexp).dirtuning);
    cd(filedir);
catch
    filedir = fullfile('Z:\analysis\',expt(iexp).mouse,expt(iexp).folder);
    cd(filedir)
    mkdir(expt(iexp).date,expt(iexp).dirtuning)
    filedir = fullfile('Z:\analysis\',expt(iexp).mouse,expt(iexp).folder,expt(iexp).date,expt(iexp).dirtuning);
    cd(filedir);
end

% load MWorks file
mworks = ['data-' 'i' expt(iexp).SubNum '-' expt(iexp).date '-' expt(iexp).dirtuning_time]; 
load (fullfile('Y:\home\andrew\Behavior\Data',mworks));

% load raw data
% load([fname '.mat']);
% data = squeeze(sbxread(fname,0,info.config.frames));
% data_avg = info.aligned.m;

%load scanbox aligned data
load([fname '.mat']);
global info
for jj = 1:info.config.frames-1
    z = double(sbxreadpacked(fname,jj-1,1));
    z2 = circshift(z,info.aligned.T(jj,:));
    data_align(1,:,:,jj) = z2;
end
data = squeeze(data_align);
data_avg = mean(data(:,:,1:100),3);
figure;imagesc(data_avg);colormap gray
save('sbx_aligned_img','data_avg')

%downsample and register image to mean from sbx aligned data
load(fullfile('Z:\analysis\',expt(iexp).mouse,expt(iexp).folder,expt(iexp).date,expt(iexp).dirtuning,'sbx_aligned_img.mat'));
load([fname '.mat']);
cd(fname(1:end-11));
fstr = fname(end-10:end);
data = squeeze(sbxread(fstr,0,info.config.frames));
down = 10;
data_down = stackGroupProject(data,down);
clear data
data_sub = data_down-min(min(min(data_down,[],1),[],2),[],3);
clear data_down
[out data_reg] = stackRegister(data_sub, data_avg);
clear data_sub

% get dF/F
nON = double(input.nScansOn)./down;
nOFF = double(input.nScansOff)./down;
nStim = double(input.gratingDirectionStepN);
nRep = size(data_reg,3)./((nON+nOFF)*nStim);
nTrials = (nStim.*nRep);
if (mod(nRep,1)) >0
    nframes = floor(nRep)*((nON+nOFF)*nStim);
    data_reg = data_reg(:,:,1:nframes);
    nRep = size(data_reg,3)./((nON+nOFF)*nStim);
    nTrials = (nStim.*nRep);
end
nOFF_ind = zeros(1,(nOFF*nStim*nRep));
nOFF_1 = zeros(1,nRep*nStim);
start = 1;
for iStim = 1:(nRep*nStim)
    nOFF_ind(1, start:start+nOFF-1) = 1+((iStim-1)*(nOFF+nON)):nOFF + ((iStim-1)*(nOFF+nON));
    nOFF_1(1,iStim) = 1+((iStim-1)*(nOFF+nON));
    start = start+nOFF;
end

nON_ind = setdiff(1:size(data_reg,3),nOFF_ind);
nON_avg = mean(data_reg(:,:,nON_ind),3);
nOFF_avg = mean(data_reg(:,:,nOFF_ind),3);

dF_data = bsxfun(@minus,data_reg, nOFF_avg);
dFoverF = bsxfun(@rdivide,dF_data,nOFF_avg);
% maxDFoverF = max(dFoverF,[],3);
% maxDFoverF(1:6,:) = 0;
maxDFoverF(:,[1:34 758:end]) = 0;

figure; imagesq(maxDFoverF); colormap(gray)
writetiff(maxDFoverF, 'maxDFoverF_meanAlign');

%% combine sbx segmentation with max df over F segmentation
load('-mat',[fname '.segment']);
sbx_mask_cell = mask;
maxDFoverF_minusSbxMask = maxDFoverF;
maxDFoverF_minusSbxMask(sbx_mask_cell > 0) = min(min(maxDFoverF));

bwout = imCellEditInteractive(maxDFoverF_minusSbxMask);
mask_cell = bwlabel(bwout);

figure;colormap hot; subplot(1,2,1);imagesq(sbx_mask_cell);subplot(1,2,2);imagesq(mask_cell);

cell_ind = mask_cell(mask_cell > 0 );
mask_cell_all = mask_cell;
mask_cell_all(cell_ind) = mask_cell_all(cell_ind)+max(max(sbx_mask_cell));
mask_cell_all = bsxfun(@plus,sbx_mask_cell,mask_cell_all);

figure;colormap hot; subplot(1,3,1);imagesq(sbx_mask_cell);subplot(1,3,2);imagesq(mask_cell);subplot(1,3,3);imagesq(mask_cell_all);

%% get timecourses
data_TC = stackGetTimeCourses(data_reg,mask_cell_all);
figure; tcOffsetPlot(data_TC)

%neuropil subtraction
nCells = size(data_TC,2);
buf = 4;
np = 6;
neuropil = imCellNeuropil(mask_cell_all,buf,np);

npTC = zeros(size(data_TC));
for i = 1:size(data_TC,2)
    tempNPmask = squeeze(neuropil(:,:,i));
    if sum(sum(tempNPmask)) > 0
    npTC(:,i) = stackGetTimeCourses(data_reg,tempNPmask);
    end
end

   %get weights by maximizing skew
ii= 0.01:0.01:1;
x = zeros(length(ii), nCells);
for i = 1:100
    x(i,:) = skewness(data_TC-tcRemoveDC(npTC*ii(i)));
end
[max_skew ind] =  max(x,[],1);
   % skew(buf,:) = max_skew;
np_w = 0.01*ind;
data_TC_subNP = data_TC-bsxfun(@times,tcRemoveDC(npTC),np_w);

%% save data
save(fullfile('Z:\Analysis',expt(iexp).mouse,expt(iexp).folder,expt(iexp).date,expt(iexp).dirtuning,'sbxaligned_data'),'mask_cell_all','mask_cell','neuropil','data_TC','data_TC_subNP','npTC');
save(fullfile('Z:\Analysis',expt(iexp).mouse,expt(iexp).folder,expt(iexp).date,expt(iexp).dirtuning,'sbxaligned_timecourses'),'data_TC_subNP')