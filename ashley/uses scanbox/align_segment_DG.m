awFSAVdatasets
iexp = 4;
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
    mkdir(date,ImgFolder)
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

%downsample and register image to mean from 
down = 10;
nON = double(input.nScansOn)./down;
nOFF = double(input.nScansOff)./down;
nStim = double(input.gratingDirectionStepN);
data_down = stackGroupProject(data,down);
clear data
data_sub = data_down-min(min(min(data_down,[],1),[],2),[],3);
clear data_down
data_reg = data_sub;

% get dF/F
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
maxDFoverF = max(dFoverF,[],3);
figure; imagesq(maxDFoverF); colormap(gray)
writetiff(maxDFoverF, 'maxDFoverF_meanAlign');

bwout = imCellEditInteractive(maxDFoverF);
mask_cell = bwlabel(bwout);
%% combine sbx segmentation with max df over F segmentation
