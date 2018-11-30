%1. try not downsampling the data (don' do stackGroupProject)--- tried
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

frame_start = 0;
nframes = 1000;
imgread = squeeze(sbxread(file,frame_start,nframes));
writetiff(imgread,[image_analysis_base,sessions,'\' sessions '_'...
    order '_tiff_', num2str(frame_start), '_', num2str(nframes)]);
%motion register (use stack register, needs images that need to be registered 
% and reference of what these images need to be aligned to. Output: how
% much each frame moved and registed images
stable_int = (180:330);
img_ref = mean(imgread(:,:,stable_int),3);
writetiff(img_ref,[image_analysis_base,sessions,'\' sessions '_'...
    order '_avetiffref_', num2str(stable_int(1)), '_', num2str(stable_int(end))]);
[rgstr_out,img_rgstr] = stackRegister(imgread, img_ref);
% look at the part of the registered image and save the new array
idx = (1:2:nframes);
writetiff(img_rgstr(:,:,idx),[image_analysis_base,sessions,'\' sessions '_'...
    order '_rgstr_tiff_', num2str(frame_start), '_', num2str(nframes)]);
save([image_analysis_base,sessions,'\' sessions '_' order '_img_rgstr.mat'], 'img_rgstr', 'rgstr_out', 'img_ref');


%% PCA: reduce the dimensions that we don't need (the dimensions that do not
%contribute much to the varience of the data)
%first, down sample the raw data by averaging every 3 frames (reduces the
%temporal resolution but Jake said it's helpful for the PCA analysis

%img_ave = stackGroupProject(img_rgstr, 3); %average together every 3 frames
%imgsize = size(img_ave,3);
nPCA = 300; % the command window tells you "First #PCs contain #% of variance
            % if the nPCA is too good (contains all of the variances, ICA
            % will give you warinings (matrix is cloase to singular or
            % badly saled) and the # of convergence in rounds will be very
            % low like 3.
%run the PCA
[mixedsig_PCA, mixedfilters_PCA, CovEvals_PCA, ~, ~, ~] = CellsortPCA_2P_Jin(img_rgstr,[1 imgsize], nPCA,[], []);
% visualize 
figure;
for i = 1:100
subplot(10,10,i); imagesc(mixedfilters_PCA(:,:,i));
colormap gray
end
figure; plot(CovEvals_PCA); % looks like can take 2-60th PCA or sth
save([image_analysis_base,sessions, '\', sessions '_' order '_PCA_variables_', num2str(nPCA),'.mat'], 'mixedsig_PCA', 'mixedfilters_PCA', 'CovEvals_PCA', 'nPCA');
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
                                               %convergence in 101 rounds
 icasig = permute(mixedfilters_ICA,[2,3,1]);
                    ICA_variables.mu = mu; ICA_variables.nIC = nIC;  ICA_variables.termtol = termtol; ICA_variables.naxrounds = maxrounds;
                    ICA_variables.npw = npw;  ICA_variables.nph = nph;
                    save([out_dir, 'ICA_variables_for_nPCA_', num2str(nPCA), '.mat'], 'ica_sig', 'icasig', 'ICA_variables');
%                     writetiff(icasig, [out_dir, 'icasig_for_nPCA_', num2str(nPCA)]);
                    figure; imagesc(sum(icasig,3));
                    savefig([out_dir, 'icasig_sum_', num2str(nPCA)]);                                          


