%1. try not downsampling the data (don' do stackGroupProject)--- tried?
%2. try with a bigger dataset (more frames, like 300000).
%   The way Court did it: take first 6000 frames of the 300000, and take every
%   other frame --> feed these 3000 to PCA. 
%3. try different nPCAs. Does having a bigger number help/hurt? Is there any
%   difference b/t taking the first 100 of 1000 principle component vs. just
%   give it nPCA = 100 in the first place? 
%4. look at your registered movies. How Court did it was to average every
%   500 frames,and choose the one that looks the sharpest as the reference
%   (the ave image is the reference)
% mask the edges after registration?






%% SECTION ONE - assign pathnames and datasets to be analyzed/written. 
clear;
%NEED TO UPDATE THIS SO IT ACCESSES SPREADSHEET INSTEAD OF JUST WRITING IN THE NAMES
sessions = '1023_img1016'; 
%ID = '1016';
image_source_base  = 'Z:\Data\2photon\'; %location of permanently stored image files for retreiving meta data
image_analysis_base    = 'Z:\Analysis\2photon_test\'; %stores the data on crash in the movingDots analysis folder
%image_source = [image_source_base, sessions,'\',ID,'\'];
image_source = [image_source_base, sessions,'\'];
image_analysis_dest = [image_analysis_base, sessions, '\'];
%% 
cd(image_source);
order = '000';
%file = [ID '_000_' order];
file = [sessions '_000_' order];
%% Motion correct
%First, find a stack of frames that do not shift, this is the refrence for motion correct.
% read part of the sbx file, write into tiff, this allows you to look at
% the data in imageJ. Look at the tiff stack, find a stack of about 100
% frames that doesn't shift by eyeballing. 

frame_start = 100;
nframes = 500;
load ([file,'.mat']);
imgread = squeeze(sbxread(file,0,1));
writetiff(imgread,[image_analysis_base,sessions,'\' sessions '_'...
    order '_tiff_', num2str(frame_start), '_', num2str(nframes)]);
%motion register (use stack register, needs images that need to be registered 
% and reference of what these images need to be aligned to. Output: how
% much each frame moved and registed images
stable_int = (180:330);
img_ref = mean(imgread(:,:,stable_int),3);
writetiff(img_ref,[image_analysis_base,sessions,'\' sessions '_'...
    order '_avetiffref_', num2str(stable_int(1)), '_', num2str(stable_int(end))]);
[rgs_out,img_rgs] = stackRegister(imgread, img_ref);
% look at the part of the registered image and save the new array
idx = (1:2:nframes);
writetiff(img_rgs(:,:,idx),[image_analysis_base,sessions,'\' sessions '_'...
    order '_rgstr_tiff_', num2str(frame_start), '_', num2str(nframes)]);
save([image_analysis_base,sessions,'\' sessions '_' order '_img_rgstr.mat'], 'img_rgstr', 'rgstr_out', 'img_ref');


%% PCA: reduce the dimensions that we don't need (the dimensions that do not
%contribute much to the varience of the data)
%first, down sample the raw data by averaging every 3 frames (reduces the
%temporal resolution but Jake said it's helpful for the PCA analysis

%img_ave = stackGroupProject(img_rgstr, 3); %average together every 3 frames
%imgsize = size(img_ave,3);
imgsize = size (img_rgs,3);
nPCA = 300; % the command window tells you "First #PCs contain #% of variance
            % if the nPCA is too good (contains all of the variances, ICA
            % will give you warinings (matrix is cloase to singular or
            % badly saled) and the # of convergence in rounds will be very
            % low like 3.
%run the PCA
[mixedsig_PCA, mixedfilters_PCA, CovEvals_PCA, ~, ~, ~] = CellsortPCA_2P_Jin(img_ave,[1 imgsize], nPCA,[], []);
% if don't down sample, first 300PCs contain 34.3% of the variance. If down
% sample, first 300PCs contain 91.7% of the variance
% visualize 
figure;
for i = 1:100
subplot(10,10,i); imagesc(mixedfilters_PCA(:,:,i));
colormap gray
end
figure; plot(CovEvals_PCA); % looks like can take 2-60th PCA or sth
save([image_analysis_base,sessions, '\', sessions '_' order '_PCA_variables_dsamped_', num2str(nPCA),'.mat'], 'mixedsig_PCA', 'mixedfilters_PCA', 'CovEvals_PCA', 'nPCA');
% save figure?

%% ICA: seperates independent spatial and temporal components
PCuse =       1:size(mixedsig_PCA,1);% try 3:70
mu =          0.5; % weight of temporal info in spatio-teporal ICA
nIC =         299; % cannot be bigger than nPCA. does it hurt if it's too small?
ica_A_guess = [];
termtol =      1e-6;
maxrounds =   400;
npw =         264;
nph =         796;
[ica_sig, mixedfilters_ICA, ica_A, numiter] = CellsortICA_2P(mixedsig_PCA,...
                                                             mixedfilters_PCA,...
                                                             CovEvals_PCA, ...
                                                             PCuse,...
                                                             mu,...
                                                             nIC,...
                                                             [],...
                                                             termtol,...
                                                             maxrounds);
                                               %convergence in 74 rounds
 icasig = permute(mixedfilters_ICA,[2,3,1]);% make it the same dimension as the raw data (pixel by pixel by time)
 ICA_variables.mu = mu; ICA_variables.nIC = nIC;  ICA_variables.termtol = termtol; ICA_variables.naxrounds = maxrounds;
 ICA_variables.npw = npw;  ICA_variables.nph = nph;
 save([image_analysis_base,sessions,'\' sessions '_' order ...
     '_ICA_variables_for_nPCA_', num2str(nPCA), '.mat'], 'ica_sig', 'icasig', 'ICA_variables');
 writetiff(icasig, [image_analysis_base,sessions,'\' sessions '_'...
    order '_icasig_for_nPCA', num2str(nPCA)]);
 figure; imagesc(sum(icasig,3));
 savefig([image_analysis_base,sessions,'\' sessions '_'...
    order '_icasig_sum_', num2str(nPCA)]);

 
%select which ICs to keep based on morphology
%disp(['Beginning IC selection for ', sessions, ' ',order, ' nPCA=', num2str(nPCA)])
%IC_use = IC_manual_check(icasig);
%save([out_dir, 'ICA_variables_for_nPCA_', num2str(nPCA), '.mat'], 'IC_use', '-append');

%% smooth and threshold the mask
%stack filter - acts as a low pass spatial filter, reduces noise. It will apply the Gaussian filter, and filter out low spatial frequency noise such as small blips and bloops in the ICs which do not belong to dendrites.
icasig_filt = stackFilter(icasig);

%set threshold a threshold for which pixels to include in a given dendrite's mask.
nIC = size(icasig_filt, 3);
cluster_threshold = 97; % this is using the top 3 percent of the fluorescence values, so brightest 3% is yes (1), and the rest is no (0)
mask_cell = zeros(size(icasig_filt));
sm_logical = zeros(npw,nph);
%bwimgcell = zeros(size(icasig2));
for ic = 1:nIC
    %convert to a binary mask (0 and 1)
    icasig_filt(:,:,ic) = imclearborder(icasig_filt(:,:,ic));
    sm_logical((icasig_filt(:,:,ic)> mean([max(prctile(icasig_filt(:,:,ic),cluster_threshold,1)) max(prctile(icasig_filt(:,:,ic),cluster_threshold,2))])))=1;
    sm_logical((icasig_filt(:,:,ic)<=mean([max(prctile(icasig_filt(:,:,ic),cluster_threshold,1)) max(prctile(icasig_filt(:,:,ic),cluster_threshold,2))])))=0;
    sm_logical = logical(sm_logical);
    %bwlabel identifies unconnected objects within a single mask. So 0=background 1=Object#1 2=Object#2
    if sum(sum(sm_logical,2),1) <51
        sm_logical = 0;
    end
    mask_cell(:,:,ic) = bwlabel(sm_logical); 
end
save([image_analysis_base,sessions,'\' sessions '_' order, '_nPCA', ...
    num2str(nPCA), '_mu', num2str(mu), '_nIC', num2str(nIC), '_thresh', ...
    num2str(cluster_threshold), '_mask_cell.mat'], 'mask_cell');

%visualize the masks
figure;
figure('rend', 'painters', 'pos', [50 150 (796*1.5) (264*1.5)]); imagesc(sum(mask_cell,3));...
    title([sessions, order, '_', ' nPCA ', num2str(nPCA), ' mu ', num2str(mu), ' nIC ', num2str(nIC)]);
savefig([image_analysis_base,sessions,'\' sessions '_' order, '_nPCA', num2str(nPCA),...
    '_mu', num2str(mu), '_nIC', num2str(nIC), '_thresh', num2str(cluster_threshold), '_mask_cell_sum.fig']);

figure;
for i = 1:100
subplot(10,10,i); imagesc(mask_cell(:,:,i));
colormap gray
end
savefig([image_analysis_base,sessions,'\' sessions '_' order, '_nPCA', num2str(nPCA),...
    '_mu', num2str(mu), '_nIC', num2str(nIC), '_thresh', num2str(cluster_threshold), '_mask_cell_indi.fig']);


%% seprate overlapping masks and combine the ones that have a high correlation
%split individual masks, remove small masks, deal with overlapping
mask_final = processMask(mask_cell);
mask_raw = reshape(mask_final, npw, nph);
figure; imagesc(mask_raw); truesize;
savefig([image_analysis_base,sessions,'\' sessions '_' order, '_nPCA', num2str(nPCA),...
    '_mu', num2str(mu), '_nIC', num2str(nIC), '_thresh', num2str(cluster_threshold), '_mask_process.fig']);

% combine highly correlated ICs to one
threshold = 0.8;
[ ~, mask3D, ~] = finalMask_Jin(img_rgs, mask_final, threshold);
figure; imagesc(sum(mask3D,3)); truesize; % got the same thing as the figure above
savefig([image_analysis_base,sessions,'\' sessions '_' order, '_nPCA', num2str(nPCA),...
    '_mu', num2str(mu), '_nIC', num2str(nIC), '_thresh', num2str(cluster_threshold), '_mask_final.fig']);

% get TCs
nmask = size(mask3D,3);
FrameRate = 30;
tc_avg = getTC(img_rgs, mask3D, nmask);

rgstr_sum = sum(img_rgs,3);
plotTC_Jin(tc_avg, mask3D, rgstr_sum, 1:size(tc_avg,2), FrameRate);
savefig([image_analysis_base,sessions,'\' sessions '_' order, '_nPCA', num2str(nPCA),...
    '_mu', num2str(mu), '_nIC', num2str(nIC), '_thresh', num2str(cluster_threshold), '_TC.fig']);
mask_flat = plotMask_Jin(mask3D,1); %second input is either 0 or 1, if 1, makes a figure. 
savefig([image_analysis_base,sessions,'\' sessions '_' order, '_nPCA', num2str(nPCA),...
    '_mu', num2str(mu), '_nIC', num2str(nIC), '_thresh', num2str(cluster_threshold), '_mask_plotMask.fig']);

data_corr = corrcoef(tc_avg);
figure; fig = imagesc(data_corr);
savefig([image_analysis_base,sessions,'\' sessions '_' order, '_nPCA', num2str(nPCA),...
    '_mu', num2str(mu), '_nIC', num2str(nIC), '_thresh', num2str(cluster_threshold), '_corrcoef.fig']);






