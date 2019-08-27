%main_CRP_comb_days

clear
file_info_CRP_all;
behav_dir = 'Z:\Data\2P_imaging\behavior\';
crp_dir = 'Z:\Analysis\Cue_reward_pairing_analysis\2P\PCA_ICA_recalibration';
match_dir_base = 'Z:\Analysis\Cue_reward_pairing_analysis\2P\PCA_ICA_recalibration\FoV_matching';
reRunPCA = 0; PCA_iteration =1; f_of_v =2; 

for sub = [6,14]%[3,6,14]%
    for learning_day = 1:2
        %assign condition
        if learning_day ==1
            days_subset = days_1;
        elseif learning_day ==2
            days_subset = days_post;
        end
        
        %collect session/mouse information
        session = days_subset{sub};
        session_date = days_subset{sub}(1:6);
        if session(end-2) == 'g'
            mouse_num = ['9', session(end-1:end)];
            mouse_ID = ['img', session(end-1:end)];
        elseif session(end-2) =='0'
            mouse_num = session(end-2:end);
            mouse_ID = ['img', session(end-2:end)];
        end
        session
        
        %set pathnames
        data_dir = fullfile('Z:\Data\2P_imaging',session, mouse_ID, '\');
        for rID = 1:2
            if exist([data_dir, mouse_ID, '_000_', runID{rID}, '.sbx'], 'file') == 2
                break
            end
        end
        if strcmp(session, '170426_img91') | strcmp(session, '170522_img89')
            rID=2;
        end
        out_dir =  fullfile(crp_dir, [session_date, '_', runID{rID}, '_', mouse_ID], '\');
        data_dir = fullfile('Z:\Data\2P_imaging',session, mouse_ID, '\');
        
        %load img_reg x2
        load([out_dir, '\img_reg.mat']);
        disp(['img_reg data loaded']);
        
        %trim img_reg
        if learning_day==1
            img_reg_d1_1 = img_reg([match_fov_y_d1{sub}],[match_fov_x_d1{sub}],:);
            img_reg_d1_2 = img_reg([match_fov_y_d1_2{sub}],[match_fov_x_d1_2{sub}],:);
        else learning_day==2
            img_reg_pl_1 = img_reg([match_fov_y_pl{sub}],[match_fov_x_pl{sub}],:);
            img_reg_pl_2 = img_reg([match_fov_y_pl_2{sub}],[match_fov_x_pl_2{sub}],:);
        end
        clear img_reg reg_out img_reg_laser_on_ind_conserv img_nframes
    end
    clear out_dir
    
    %assign the output directories
    match_dir = [match_dir_base, '\', mouse_ID, '_D1PL\'];
    if ~exist(match_dir)
        mkdir(match_dir);
    end
    match_dir_fov_1 = [match_dir, 'fov_1\'];
    if ~exist(match_dir_fov_1);
        mkdir(match_dir_fov_1);
    end
    match_dir_fov_2 = [match_dir, 'fov_2\'];
    if ~exist(match_dir_fov_2);
        mkdir(match_dir_fov_2);
    end

    %% Motion registration (skip)
%     img_ref = mean(img_reg_d1(:,:,[1:2000]),3);
%     [reg_out_pl_2_d1, img_reg_pl_2_d1] = stackRegister(img_reg_pl, img_ref);
%     stack_reg_corr.d1_2_PL_b4_reg = corr2( mean(img_reg_pl(:,:,[1:2000]),3), img_ref);
%     stack_reg_corr.d1_2_PL = corr2( mean(img_reg_pl_2_d1(:,:,[1:2000]),3),   img_ref);
%     stack_reg_corr.d1_2_d1 = corr2( mean(img_reg_pl_2_d1(:,:,[1:2000]),3),   mean(img_reg_pl_2_d1(:,:,[2001:4000]),3))
%     save([match_dir, 'img_reg_D1PL.mat'], 'img_reg_pl_2_d1',  'img_reg_d1' ,'img_ref', 'stack_reg_corr', '-v7.3');

save([match_dir_fov_1, 'img_reg_D1PL.mat'], 'img_reg_d1_1', 'img_reg_pl_1');
save([match_dir_fov_2, 'img_reg_D1PL.mat'], 'img_reg_d1_2', 'img_reg_pl_2');
disp('img_reg variables saved')

figure;
subplot(2,1,1);   % imagesc(mean(img_reg_pl_2_d1(:,:,[1:2000]),3));
imagesc(mean(img_reg_pl_1(:,:,[1:2000]),3));
title([mouse_ID, ' day N']);
subplot(2,1,2); %imagesc(img_ref);
imagesc(mean(img_reg_d1_1(:,:,[1:2000]),3));  title('day 1'); suptitle([mouse_ID, ' fov 1']);
savefig([match_dir_fov_1, 'D1PL_fov_1.fig']);%savefig([match_dir, 'd1_stackReg_PL.fig']);
figure;
subplot(2,1,1);   % imagesc(mean(img_reg_pl_2_d1(:,:,[1:2000]),3));
imagesc(mean(img_reg_pl_2(:,:,[1:2000]),3));
title([mouse_ID, ' day N']);
subplot(2,1,2); %imagesc(img_ref);
imagesc(mean(img_reg_d1_2(:,:,[1:2000]),3));  title('day 1'); suptitle([mouse_ID, ' fov 2']);
savefig([match_dir_fov_2, 'D1PL_fov_2.fig']);%savefig([match_dir, 'd1_stackReg_PL.fig']);

%% compute PRINCIPAL COMPONENTS
img_fn = [session_date, '_', mouse_ID];
nPCA_d1 = all_nPCA_d1dn{sub};
nPCA_pl = all_nPCA_d1dn{sub};
if reRunPCA==1
    load([match_dir, 'img_reg_D1PL.mat']);
    if PCA_iteration ==2
        load([match_dir, 'mask_cell.mat']);
    else
        load([match_dir, 'mask_cell', num2str(PCA_iteration-1), '.mat']);
    end
    nPCA_d1 = nPCA_d1-round(nPCA_d1/3);
    nPCA_pl = nPCA_pl-round(nPCA_pl/3);
    sz_d1 = size(img_reg_d1);  sz_pl = size(img_reg_pl_2_d1);
    %get a logical version of the summed 3d masks and convert to linear indexing
    mask_ind_d1 = find(logical(sum(mask_cell_d1,3)));
    mask_ind_pl = find(logical(sum(mask_cell_pl,3)));
    %blackout those pixels in the registered movies
    img_reg_d1 = reshape(img_reg_d1, sz_d1(1)*sz_d1(2), sz_d1(3));
    img_reg_pl_2_d1 = reshape(img_reg_pl_2_d1, sz_pl(1)*sz_pl(2), sz_pl(3));
    img_reg_d1([mask_ind_d1],:) = 0;
    img_reg_pl_2_d1([mask_ind_pl],:) = 0;
    img_reg_d1 = reshape(img_reg_d1, sz_d1(1), sz_d1(2), sz_d1(3));
    img_reg_pl_2_d1 = reshape(img_reg_pl_2_d1, sz_pl(1), sz_pl(2), sz_pl(3));
end

%avg every X frames together 
if f_of_v >1
    img_pca_d1 = stackGroupProject(img_reg_d1_1, 3); %average together every 3 frames
    img_pca_pl = stackGroupProject(img_reg_pl_1, 3); %average together every 3 frames
    img_pca_d1_2 = stackGroupProject(img_reg_d1_2, 3); %average together every 3 frames
    img_pca_pl_2 = stackGroupProject(img_reg_pl_2, 3); %average together every 3 frames
else
    img_pca_d1 = stackGroupProject(img_reg_d1, 3); %average together every 3 frames
    img_pca_pl = stackGroupProject(img_reg_pl_2_d1, 3); %average together every 3 frames
end

%run the PCA
[mixedsig_PCA_d1, mixedfilters_PCA_d1, CovEvals_PCA_d1, ~, ~, ~] = CellsortPCA_2P(img_pca_d1,[1 size(img_pca_d1,3)], nPCA_d1,[], match_dir, img_fn, []);
[mixedsig_PCA_pl, mixedfilters_PCA_pl, CovEvals_PCA_pl, ~, ~, ~] = CellsortPCA_2P(img_pca_pl,[1 size(img_pca_pl,3)], nPCA_pl,[], match_dir, img_fn, []);
nph = length(match_fov_y_d1{sub}); npw = length(match_fov_x_pl{sub});
if reRunPCA==1
    save([match_dir, 'PCA_variables', num2str(PCA_iteration),'.mat'], 'mixedsig_PCA_d1', 'mixedfilters_PCA_d1', 'CovEvals_PCA_d1', 'mixedsig_PCA_pl', 'mixedfilters_PCA_pl', 'CovEvals_PCA_pl', 'nPCA_d1', 'nPCA_pl', 'npw', 'nph');
elseif f_of_v > 1
    save([match_dir_fov_1, 'PCA_variables.mat'], 'mixedsig_PCA_d1', 'mixedfilters_PCA_d1', 'CovEvals_PCA_d1', 'mixedsig_PCA_pl', 'mixedfilters_PCA_pl', 'CovEvals_PCA_pl', 'nPCA_d1', 'nPCA_pl', 'npw', 'nph');
    [mixedsig_PCA_d1, mixedfilters_PCA_d1, CovEvals_PCA_d1, ~, ~, ~] = CellsortPCA_2P(img_pca_d1_2,[1 size(img_pca_d1_2,3)], nPCA_d1,[], match_dir, img_fn, []);
    [mixedsig_PCA_pl, mixedfilters_PCA_pl, CovEvals_PCA_pl, ~, ~, ~] = CellsortPCA_2P(img_pca_pl_2,[1 size(img_pca_pl_2,3)], nPCA_pl,[], match_dir, img_fn, []);
    nph = length(match_fov_y_d1{sub}); npw = length(match_fov_x_pl{sub});
    save([match_dir_fov_2, 'PCA_variables.mat'], 'mixedsig_PCA_d1', 'mixedfilters_PCA_d1', 'CovEvals_PCA_d1', 'mixedsig_PCA_pl', 'mixedfilters_PCA_pl', 'CovEvals_PCA_pl', 'nPCA_d1', 'nPCA_pl', 'npw', 'nph');
else
    save([match_dir, 'PCA_variables.mat'], 'mixedsig_PCA_d1', 'mixedfilters_PCA_d1', 'CovEvals_PCA_d1', 'mixedsig_PCA_pl', 'mixedfilters_PCA_pl', 'CovEvals_PCA_pl', 'nPCA', 'npw', 'nph');
end
continue 

load([match_dir_fov_1, 'PCA_variables.mat'], 'mixedfilters_PCA_d1', 'mixedfilters_PCA_pl');
 writetiff(mixedfilters_PCA_d1, [match_dir_fov_1, 'PCuse_d1']);
 writetiff(mixedfilters_PCA_pl, [match_dir_fov_1, 'PCuse_pl']);
 load([match_dir_fov_2, 'PCA_variables.mat'], 'mixedfilters_PCA_d1', 'mixedfilters_PCA_pl');
 writetiff(mixedfilters_PCA_d1, [match_dir_fov_2, 'PCuse_d1']);
 writetiff(mixedfilters_PCA_pl, [match_dir_fov_2, 'PCuse_pl']);
continue

load([match_dir, 'PCA_variables', num2str(PCA_iteration),'.mat']);
writetiff(mixedfilters_PCA_d1, [match_dir_fov_1, 'PCuse_d1_', num2str(PCA_iteration)]);
writetiff(mixedfilters_PCA_pl, [match_dir_fov_1, 'PCuse_pl_', num2str(PCA_iteration)]);
continue

%% Compute independent components
%Load variables and assign parameters for ICA
%load([match_dir, 'PCA_variables.mat']);
% load([match_dir, 'PCA_variables', num2str(PCA_iteration),'.mat']);
for fov_num = 1:f_of_v
    if fov_num ==1
        load([match_dir_fov_1, 'PCA_variables.mat']);
    elseif fov_num ==2
        load([match_dir_fov_2, 'PCA_variables.mat']);
    end
    if reRunPCA==1
        PCuse_d1 = [1:size(mixedsig_PCA_d1,1)-5];
        PCuse_pl = [1:size(mixedsig_PCA_pl,1)-5];
    else
        PCuse_d1 = all_PCuse_d1{sub};
        PCuse_pl = all_PCuse_pl{sub};
    end
    nIC_d1 = round(length(PCuse_d1)*0.9);
    nIC_pl = round(length(PCuse_pl)*0.9);
    mu = 0.97; % spatial temporal mix
    termtol = 0.00001;
    maxrounds = 1000;
    
    %Day 1 ICA
    [ica_sig_d1, mixedfilters_ICA_d1, ica_A_d1, numiter_d1] = CellsortICA_2P(mixedsig_PCA_d1, mixedfilters_PCA_d1, CovEvals_PCA_d1, PCuse_d1, mu, nIC_d1, [], termtol, maxrounds);
    icasig_d1 = permute(mixedfilters_ICA_d1,[2,3,1]);
    %Post learning ICA
    [ica_sig_pl, mixedfilters_ICA_pl, ica_A_pl, numiter_pl] = CellsortICA_2P(mixedsig_PCA_pl, mixedfilters_PCA_pl, CovEvals_PCA_pl, PCuse_pl, mu, nIC_pl, [], termtol, maxrounds);
    icasig_pl = permute(mixedfilters_ICA_pl,[2,3,1]);
    %saveoutputs
    ICA_variables.mu = mu; ICA_variables.nIC_d1 = nIC_d1; ICA_variables.nIC_pl = nIC_pl; ICA_variables.termtol = termtol; ICA_variables.naxrounds = maxrounds;
    ICA_variables.npw = npw;  ICA_variables.nph = nph;
    %save([match_dir, 'ICA_variables.mat'], 'ica_sig_d1', 'icasig_d1', 'ica_sig_pl', 'icasig_pl', 'ICA_variables');
    %save([match_dir, 'ICA_variables', num2str(PCA_iteration), '.mat'], 'ica_sig_d1', 'icasig_d1', 'ica_sig_pl', 'icasig_pl', 'ICA_variables');
    save([match_dir, 'ICA_variables', num2str(PCA_iteration), '.mat'], 'ica_sig_d1', 'icasig_d1', 'ica_sig_pl', 'icasig_pl', 'ICA_variables');
    if fov_num ==1
         save([match_dir_fov_1, 'ICA_variables.mat'], 'ica_sig_d1', 'icasig_d1', 'ica_sig_pl', 'icasig_pl', 'ICA_variables');
    elseif fov_num ==2
         save([match_dir_fov_2, 'ICA_variables.mat'], 'ica_sig_d1', 'icasig_d1', 'ica_sig_pl', 'icasig_pl', 'ICA_variables');
    end
    
    %plot icasig so you can see how much the ICs are overlapping
    figure;
    subplot(2,1,1); imagesc(sum(icasig_d1,3)); title([mouse_ID, ' ICs day 1']);
    subplot(2,1,2); imagesc(sum(icasig_pl,3)); title([mouse_ID, ' ICs post learning']);
    %savefig([match_dir, 'icasig_sum_d1pl']);
    %savefig([match_dir, 'icasig_sum_d1pl', num2str(PCA_iteration)]);
    if fov_num ==1
         savefig([match_dir_fov_1, 'icasig_sum_d1pl']);
    elseif fov_num ==2
         savefig([match_dir_fov_2, 'icasig_sum_d1pl']);
    end
end
continue
   
%     nph=ICA_variables.nph; npw = ICA_variables.npw; mu = ICA_variables.mu;

%% Process the masks from the PCA/ICA

%select which ICs to keep based on morphology
%load([match_dir, 'ICA_variables.mat']);
%load([match_dir, 'ICA_variables', num2str(PCA_iteration),'.mat']);
for fov_num =1:f_of_v
    if fov_num ==1
        load([match_dir_fov_1, 'ICA_variables.mat']);
    elseif fov_num ==2
        load([match_dir_fov_2, 'ICA_variables.mat']);
    end
    disp(['Beginning IC selection for ', mouse_ID])
    IC_use_d1 = IC_manual_check(icasig_d1);
    IC_use_pl = IC_manual_check(icasig_pl);
    %save([match_dir, 'ICA_variables.mat'], 'IC_use_d1', 'IC_use_pl', '-append');
    %save([match_dir, 'ICA_variables', num2str(PCA_iteration), '.mat'], 'IC_use_d1', 'IC_use_pl', '-append');
    if fov_num ==1
        save([match_dir_fov_1, 'ICA_variables.mat'], 'IC_use_d1', 'IC_use_pl', '-append');
    elseif fov_num ==2
        save([match_dir_fov_2, 'ICA_variables.mat'], 'IC_use_d1', 'IC_use_pl', '-append');
    end
    figure; subplot(2,1,1); imagesc(sum(icasig_d1(:,:,[find(IC_use_d1)]),3)); title([mouse_ID, ' day 1']);
    subplot(2,1,2); imagesc(sum(icasig_pl(:,:,[find(IC_use_pl)]),3)); title([mouse_ID, ' day N']);
end
 continue

for fov_num =1:f_of_v
    if fov_num ==1
        load([match_dir_fov_1, 'ICA_variables.mat']);
    elseif fov_num ==2
        load([match_dir_fov_2, 'ICA_variables.mat']);
    end
%stack filter - acts as a low pass spatial filter. It will filter out low spatial frequency noise such as small blips and bloops in the ICs which do not belong to dendrites.
icasig_d1 = icasig_d1(:,:,[find(IC_use_d1)]);
icasig_pl = icasig_pl(:,:,[find(IC_use_pl)]);
icasig2_d1 = stackFilter(icasig_d1);
icasig2_pl = stackFilter(icasig_pl);
%assign variables and allocate memory
nIC_d1 = size(icasig_d1,3);
nIC_pl = size(icasig_pl,3);
icasig3_d1 = icasig2_d1;
icasig3_pl = icasig2_pl;
                    
%set threshold a threshold for which pixels to include in a given dendrite's mask.
cluster_threshold = 93; % 90:100; %97- img90 97- img91 97- img070
mask_cell_d1 = zeros(size(icasig2_d1));
mask_cell_pl = zeros(size(icasig2_pl));
% sm_logical_d1 = zeros( ICA_variables.nph, ICA_variables.npw);
% sm_logical_pl = zeros( ICA_variables.nph, ICA_variables.npw);
sm_logical_d1 = zeros( size(mask_cell_d1,1), size(mask_cell_d1,2));
sm_logical_pl = zeros( size(mask_cell_pl,1), size(mask_cell_pl,2));
%bwimgcell = zeros(size(icasig2));
for ic = 1:nIC_d1
    %convert to a binary mask
    icasig3_d1(:,:,ic) = imclearborder(icasig2_d1(:,:,ic));
    sm_logical_d1((icasig3_d1(:,:,ic)> mean([max(prctile(icasig3_d1(:,:,ic),cluster_threshold,1)) max(prctile(icasig3_d1(:,:,ic),cluster_threshold,2))])))=1;
    sm_logical_d1((icasig3_d1(:,:,ic)<=mean([max(prctile(icasig3_d1(:,:,ic),cluster_threshold,1)) max(prctile(icasig3_d1(:,:,ic),cluster_threshold,2))])))=0;
    sm_logical_d1 = logical(sm_logical_d1);
    %bwlabel identifies unconnected objects within a single mask. So 0=background 1=Object#1 2=Object#2
%     if sum(sum(sm_logical_d1,2),1) <51
%         sm_logical_d1 = 0;
%     end
    mask_cell_d1(:,:,ic) = bwlabel(sm_logical_d1);
end 
for ic = 1:nIC_pl
    %convert to a binary mask
    icasig3_pl(:,:,ic) = imclearborder(icasig2_pl(:,:,ic));
    sm_logical_pl((icasig3_pl(:,:,ic)> mean([max(prctile(icasig3_pl(:,:,ic),cluster_threshold,1)) max(prctile(icasig3_pl(:,:,ic),cluster_threshold,2))])))=1;
    sm_logical_pl((icasig3_pl(:,:,ic)<=mean([max(prctile(icasig3_pl(:,:,ic),cluster_threshold,1)) max(prctile(icasig3_pl(:,:,ic),cluster_threshold,2))])))=0;
    sm_logical_pl = logical(sm_logical_pl);
    %bwlabel identifies unconnected objects within a single mask. So 0=background 1=Object#1 2=Object#2
%     if sum(sum(sm_logical_pl,2),1) <51
%         sm_logical_pl = 0;
%     end
    mask_cell_pl(:,:,ic) = bwlabel(sm_logical_pl);
end

%
%save([match_dir, 'mask_cell',num2str(PCA_iteration),'.mat'], 'mask_cell_d1', 'mask_cell_pl');
figure; subplot(2,1,1); imagesc(sum(mask_cell_d1,3)); title([mouse_ID, 'day 1']);
subplot(2,1,2); imagesc(sum(mask_cell_pl,3)); title([mouse_ID, 'day N']);
if fov_num ==1
    savefig([match_dir_fov_1, 'mask_cell.fig']);
    save([match_dir_fov_1, 'mask_cell.mat'], 'mask_cell_d1', 'mask_cell_pl');
elseif fov_num ==2
    save([match_dir_fov_2, 'mask_cell.mat'], 'mask_cell_d1', 'mask_cell_pl');
end
%load([match_dir, 'mask_cell.mat']);
figure; subplot(1,3,1); imagesc(logical(sum(mask_cell_d1,3))); title('day 1 mask');
subplot(1,3,2); imagesc(logical(sum(mask_cell_pl,3))); title('day N mask');
subplot(1,3,3); imshowpair( logical(sum(mask_cell_d1,3)), logical(sum(mask_cell_pl,3))  ); title(['green=day1, purple=dayN, gray = overlap']);
suptitle(mouse_ID);
%savefig([match_dir, 'PCAICA_mask_overlap_d1pl_',num2str(PCA_iteration),'.fig']);
if fov_num ==1
    savefig([match_dir_fov_1, 'PCAICA_mask_overlap_d1pl.fig']);
elseif fov_num ==2
    savefig([match_dir_fov_2, 'PCAICA_mask_overlap_d1pl.fig']);
end
end
%%

%                      load([out_dir, 'img_reg']);
%                     load([out_dir, 'ICA_variables_for_nPCA_', num2str(nPCA), '.mat']);
%                     mu = ICA_variables.mu; npw = ICA_variables.npw; nph = ICA_variables.nph; nIC=ICA_variables.nIC;
%                     mask_cell_file = dir(fullfile(out_dir, ['nPCA_' num2str(nPCA) '_mu_' num2str(mu) '_nIC_' '*' '_thresh_' num2str(cluster_threshold) '_mask_cell.mat']));
%                     load([out_dir, mask_cell_file.name]);
       
%                     %split individual masks if their blobs are not touching
%                     mask_cell_out = blob_buster_CRP(mask_cell);
%                     %consolidate highly correlated masks  and %combine masks that are overlapping
%                     threshold = 0.80; % correlation threshold
%                     mask3D_all = combine_corr_masks(img_reg, mask_cell_out, threshold);
%                     %cut masks <200 pixels
%                     too_small_ind = find( squeeze(sum(sum(mask3D_all,2),1))<200 );
%                     mask3D = mask3D_all;
%                     mask3D(:,:,too_small_ind) = [];
%                     %implement buffer
%                     mask_3D_buffed = make_mask_buffer(mask3D);
%                     mask2D = plotMask(mask_3D_buffed, 0,0, 1);
%                     title([session_date, ' ', mouse_ID, ' nPCA =', num2str(nPCA)]);
%                     savefig([out_dir, session, '_nPCA_', num2str(nPCA), '_post_process.fig']);
%                     save([out_dir, session, '_nPCA_', num2str(nPCA), '_post_process_outputs.mat'], 'mask2D', 'mask_3D_buffed');
%                     continue                     
%        
% %                     mask_2D_buffed = plotMask(mask_3D_buffed, 0,0, 0);
% %                     figure; image(mask_2D_buffed);
% %                     title([session, ' final mask with buffer']);
% %                     savefig([out_dir, 'nPCA_', num2str(nPCA), '_mu_', num2str(mu), '_nIC_', num2str(nIC), ...
% %                             '_thresh_', num2str(cluster_threshold), '_final_buffer', '.fig']);
              
end
