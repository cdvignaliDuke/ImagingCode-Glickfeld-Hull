% main: process movie as following:
% 1. registeration -- find the most stable frame out of 100
% 2. PCA -- #PCA depends on how dense dendrite population
% 3. ICA -- #IC is less than #PCA and with 30HZ, the temperal-spatial ratio is set to be 0.97
% 4. Apply gaussian filter to ICA
% 5. Threshold ICA signal to get binary mask
% 6. Process mask (break multiple components on single layer to multiple, combine highly correlated ones and remove small ones)
% 7. Extract timecourse (manually inspect TC and remove bad ones) and create a singal layer mask for display

clear
%file_info_CRP2;
file_info_CRP_all;
usFacs = 100;
behav_dir = 'Z:\Data\2P_imaging\behavior\';
crp_dir = 'Z:\Analysis\Cue_reward_pairing_analysis\2P\PCA_ICA_recalibration';
doRerun = 0; useWarp = 0; doPCA = 0; useGPU = 0; checkImg = 0;
for session_subset = 2
    if session_subset ==1
        stable_int = stable_int_1;
        all_nPCA = all_nPCA_1;
        all_PCuse = all_PCuse_1;
        all_PCuse2 = all_PCuse2_1;
        days_subset = days_1;
        nPCA_to_use = nPCA_to_use_1;
    elseif session_subset ==2
        stable_int = stable_int_post;
        all_nPCA = all_nPCA_post;
        all_PCuse = all_PCuse_post;
        all_PCuse2 = all_PCuse2_post;
        days_subset = days_post;
        nPCA_to_use = nPCA_to_use_post;
    end

for sub = 2:15 %  %[2, 4, 5, 8, 9, 10, 12, 15]     %[1:2, 4:(length(stable_int)-2)]; %[2,5, 9, 22, 23, 31, 3, 6, 10, 25, 27, 32]
    if sub<=2 && session_subset ==2
        continue
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
    out_dir =  fullfile(crp_dir, [session_date, '_', runID{rID}, '_', mouse_ID], '\');
    data_dir = fullfile('Z:\Data\2P_imaging',session, mouse_ID, '\');
    if ~exist(out_dir)
        mkdir(out_dir);
    end
            
            %% Extract pockel cell data
%             
%             %set directories
%             subMat = dir([behav_dir, '*', mouse_num, '-', session_date, '*']); %sets directory as the behavior file
%             
%             %load the behavior data
%             load([behav_dir, subMat.name]);
%             
%             %identify times of "lever press" and trial end
%             cStart = cell2mat(input.cLeverDown);
%             cEnd = cell2mat(input.cTrialEnd);
%             
%             %load imaging matfile
%             config_fn = dir(fullfile(data_dir,['*' runID{rID} '.mat']));
%             load([data_dir, config_fn.name]);
%             
%             if exist([data_dir, mouse_ID, '_000_', runID{rID}, '_realtime.mat'], 'file') ~= 2 %checks for data about pockel cell
%                 %for some reason Ziye overwrote ttl_log with this new version derived from the bx data. It makes no sense to I renamed it to ttl_log2 so it wont overwrite ttl_log
%                 ttl_log2 = zeros(info.config.frames,1);
%                 for jj = 1:length(cStart)
%                     ttl_log2(cStart(jj):cEnd(jj)) = 1;
%                 end
%                 save([data_dir, mouse_ID, '_', runID{rID}, '_realtime.mat'], 'ttl_log2');  %overwrites existing ttl_log data
%             end
            
            %% Motion registration or check for exisint motion reg outputs
%             %             if exist([out_dir, 'Reg_out.mat'],'file') == 2
%             %                 load([out_dir, 'Reg_out.mat']);
%             %                 [~,img_reg]=stackRegister_MA(img,img(:,:,1),usFacs,double(reg_out));
%             %             elseif exist([out_dir, 'img_reg.mat'],'file') == 2
%             %                 load([out_dir, 'img_reg.mat']);
%             %             else
%             %             [img_reg, laser_on_ind_conserv, nt] = motionCorrect(data_dir, runID{rID}, out_dir); %uses stackRegister, automatically saves img_reg and Reg_out
%             config_fn = dir(fullfile(data_dir,['*' runID{rID} '.mat']));
%             laser_power_fn = dir(fullfile(data_dir,['*' runID{rID} '_realtime.mat']));
%             [img, ~, ~] = loadFile( data_dir, runID{rID}, [], []);
%             [img_mat_file, laser_power_vec_ttl] = get_laser_power_data(data_dir, config_fn, laser_power_fn);
%             if isempty(laser_power_vec_ttl)
%                 laser_on_ind_conserv = 1:nt;
%             else
%                 laser_on_ind_conserv = conservative_laser_on(laser_power_vec_ttl);
%             end
%             laser_on_ind = find(laser_power_vec_ttl);
%             img_ref = mean(img(:,:,stable_int{sub}),3);
%             [reg_out, img_reg] = stackRegister(img(:,:,laser_on_ind_conserv(find(laser_on_ind_conserv<=size(img,3)))), img_ref);
%             save([out_dir, 'img_reg.mat'], 'img_reg', 'reg_out', 'img_ref', 'laser_on_ind_conserv', '-v7.3');
%             %end
% %            continue

            if doRerun == 0 
                 %% compute PRINCIPAL COMPONENTS 
                cd(data_dir); 
%                  load([out_dir, 'img_reg']);
%                   [npw, nph, ~] = size(img_reg);
               for this_nPCA = 1:2
%                      nPCA = all_nPCA{sub}(this_nPCA);
%                     %img_pca = img_reg(:,:,1:2:end); % downsample in time by 2 or 5
%                     img_pca = stackGroupProject(img_reg, 3); %average together every 3 frames
%                     nf = size(img_pca,3);
%                     img_fn = [session_date, '_', runID{rID}, '_', mouse_ID];
%                     
%                     %run the PCA
%                     [mixedsig_PCA, mixedfilters_PCA, CovEvals_PCA, ~, ~, ~] = CellsortPCA_2P(img_pca,[1 nf], nPCA,[], out_dir, img_fn, []);
%                     save([out_dir, 'PCA_variables_', num2str(nPCA),'.mat'], 'mixedsig_PCA', 'mixedfilters_PCA', 'CovEvals_PCA', 'nPCA');
%                     continue
%                     
                    % to view PCs
%                      load([out_dir, 'PCA_variables_', num2str(nPCA),'.mat']);
%                    writetiff(mixedfilters_PCA, [out_dir, 'PCuse_for_nPCA_', num2str(nPCA)]);
%                     continue
                    
                    %manually enter the PC #s to be used 
%                     if this_nPCA == 1
%                         PCuse = all_PCuse{sub};
%                     elseif this_nPCA ==2 
%                         PCuse = all_PCuse2{sub};
%                     end
                    
                    %% Compute independent components
%                     %mus = [0.97, 0.7];
%                     %all_nIC = {[PCuse], [PCuse(1):(PCuse(end)-50)]};
%                     %for this_mu = mus
%                       mu = 0.97; % spatial temporal mix
%                     %for this_nIC = [1:2]
%                     nIC = round(length(PCuse)*0.9); %length(all_nIC{this_nIC});%300;  %400- img90 100- img91
%                     termtol = 0.00001;
%                     maxrounds = 1000;
%                     nph = 796; npw = 264;
%                     
%                     %run ICA and save outputs
%                     [ica_sig, mixedfilters_ICA, ica_A, numiter] = CellsortICA_2P(mixedsig_PCA, mixedfilters_PCA, CovEvals_PCA, PCuse, mu, nIC, [], termtol, maxrounds);
%                     icasig = permute(mixedfilters_ICA,[2,3,1]);
%                     ICA_variables.mu = mu; ICA_variables.nIC = nIC;  ICA_variables.termtol = termtol; ICA_variables.naxrounds = maxrounds;
%                     ICA_variables.npw = npw;  ICA_variables.nph = nph;
%                     save([out_dir, 'ICA_variables_for_nPCA_', num2str(nPCA), '.mat'], 'ica_sig', 'icasig', 'ICA_variables');
% %                     writetiff(icasig, [out_dir, 'icasig_for_nPCA_', num2str(nPCA)]);
%                     figure; imagesc(sum(icasig,3));
%                     savefig([out_dir, 'icasig_sum_', num2str(nPCA)]);
%                     continue
%                     load([out_dir, 'ICA_variables_for_nPCA_', num2str(nPCA), '.mat']);
%                       load([out_dir, 'img_reg']);
%                       nph=ICA_variables.nph; npw = ICA_variables.npw; mu = ICA_variables.mu;
                    %% Process the masks from the PCA/ICA
                    
                    %select which ICs to keep based on morphology
%                     disp(['Beginning IC selection for ', session_date, ' ', mouse_ID, ' nPCA=', num2str(nPCA)])
%                     IC_use = IC_manual_check(icasig);
%                     save([out_dir, 'ICA_variables_for_nPCA_', num2str(nPCA), '.mat'], 'IC_use', '-append');
%                      continue
                    
                    %stack filter - acts as a low pass spatial filter. It will filter out low spatial frequency noise such as small blips and bloops in the ICs which do not belong to dendrites.
%                     icasig = icasig(:,:,[find(IC_use)]);
%                     icasig2 = stackFilter(icasig);
%                     %icasig2 = icasig;
%                     %assign variables and allocate memory
%                     nIC = size(icasig,3);
%                     icasig3 = icasig2;
                    
                    %set threshold a threshold for which pixels to include in a given dendrite's mask.
%                     cluster_threshold = 96; % 90:100; %97- img90 97- img91 97- img070
%                     mask_cell = zeros(size(icasig2));
%                     sm_logical = zeros(npw,nph);
%                     %bwimgcell = zeros(size(icasig2));
%                     for ic = 1:nIC
%                         %convert to a binary mask
%                         icasig3(:,:,ic) = imclearborder(icasig2(:,:,ic));
%                         sm_logical((icasig3(:,:,ic)> mean([max(prctile(icasig3(:,:,ic),cluster_threshold,1)) max(prctile(icasig3(:,:,ic),cluster_threshold,2))])))=1;
%                         sm_logical((icasig3(:,:,ic)<=mean([max(prctile(icasig3(:,:,ic),cluster_threshold,1)) max(prctile(icasig3(:,:,ic),cluster_threshold,2))])))=0;
%                         sm_logical = logical(sm_logical);
%                         %bwlabel identifies unconnected objects within a single mask. So 0=background 1=Object#1 2=Object#2
%                         if sum(sum(sm_logical,2),1) <51
%                             sm_logical = 0;
%                         end
%                         mask_cell(:,:,ic) = bwlabel(sm_logical);
%                     end
%                     %close all
                    
%                     save([out_dir, 'nPCA_', num2str(nPCA), '_mu_', num2str(mu), '_nIC_', num2str(nIC), '_thresh_', num2str(cluster_threshold), '_mask_cell.mat'], 'mask_cell')
%                     figure('rend', 'painters', 'pos', [50 150 (796*1.5) (264*1.5)]); imagesc(sum(mask_cell,3)); title([session_date, '_', mouse_ID, ' nPCA ', num2str(nPCA), ' mu ', num2str(mu), ' nIC ', num2str(nIC)]);
%                     savefig([out_dir, session, '_nPCA_', num2str(nPCA), '_mu_', num2str(mu), '_nIC_', num2str(nIC), '_thresh_', num2str(cluster_threshold), '_mask_cell.fig']);
                    
%                     load([out_dir, 'img_reg']);
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
% %                     continue                     
       
%                     mask_2D_buffed = plotMask(mask_3D_buffed, 0,0, 0);
%                     figure; image(mask_2D_buffed);
%                     title([session, ' final mask with buffer']);
%                     savefig([out_dir, 'nPCA_', num2str(nPCA), '_mu_', num2str(mu), '_nIC_', num2str(nIC), ...
%                             '_thresh_', num2str(cluster_threshold), '_final_buffer', '.fig']);
                   end
%                end
                % end
%                  continue
                
                
                %% Calculate TCs
                if sub >3
                nPCA = nPCA_to_use{sub};
                load([out_dir, 'img_reg']);
                load([out_dir, 'ICA_variables_for_nPCA_', num2str(nPCA), '.mat']);
                mu = ICA_variables.mu; npw = ICA_variables.npw; nph = ICA_variables.nph; nIC=ICA_variables.nIC;
                load([out_dir, session, '_nPCA_', num2str(nPCA), '_post_process_outputs.mat'], 'mask2D', 'mask_3D_buffed');
                if min(sum(sum(mask_3D_buffed,1),2)) < 50
                    pixel_sizes = squeeze(sum(sum(mask_3D_buffed,1),2));
                    too_small_ind = find(pixel_sizes <= 50);
                    mask_3D_buffed(:,:,[too_small_ind]) = [];
                    %mask2D(:,:,[too_small_ind]) = [];
                    save([out_dir, session, '_nPCA_', num2str(nPCA), '_post_process_outputs.mat'], 'mask2D', 'mask_3D_buffed');
                end
                nmask = size(mask_3D_buffed,3);
                FrameRate = 30;
                tc_avg = getTC(img_reg, mask_3D_buffed, nmask);
                 save([out_dir, session_date, '_', mouse_ID, '_nPCA_', num2str(nPCA), '_Results.mat'], 'tc_avg'); % 'mask_raw', 'mask_flat', 'mask3D', 'mask_final', 'data_corr'
                end
                 %                 
%                 saveData = 0;
%                 reg_sum = sum(img_reg,3);
%                 %plotTC(tc_avg, mask_3D_buffed, reg_sum, 1:size(tc_avg,2), FrameRate, out_dir, saveData);
%                 %mask_flat = plotMask(mask_3D_buffed, saveData, out_dir, 1);
%                 data_corr = corrcoef(tc_avg);
%                 sz = [npw, nph];
%                 mask_final = mask_3D_buffed;
%                 mask_final_color = mask2D;
%                 
%                 save([out_dir, dates{sub}, '_', mouseID{sub}, '_nPCA_', num2str(nPCA), 'Results.mat'], 'tc_avg', 'mask_final', 'mask_final_color', 'data_corr', 'sz');
%                 
%                 %             out_dir2  = fullfile('Z:','home','jake','Analysis','Cue_reward_pairing_analysis','2P',[dates{sub}, '_', runID{rID}, '_', mouseID{sub}], '\');
%                 %             load([out_dir2, 'Results.mat']);
%                 %             nmask = size(mask3D,3);
%                 %             FrameRate = 30;
%                 %             tc_avg = getTC(img_reg, mask3D, nmask);
%             else
%                 load([out_dir, 'Results.mat']);
%                 nmask = size(mask3D,3);
%                 FrameRate = 30;
%                 tc_avg = getTC(img_reg, mask3D, nmask);
%                 data_corr = corrcoef(tc_avg);
            end
%             
%              continue
%             tc_temp = zeros(nt, nmask);
%             tc_temp(laser_on_ind_conserv,:) = tc_avg;
%             tc_avg = tc_temp;
%             sz = [npw, nph];
%             
            
%             save([out_dir, 'ROI_TCs.mat'],'tc_avg', 'mask_final', 'sz');
%             %clearvars -except sub usFacs docombineMask doRerun useGPU outdir
%             %close all
            

            nPCA = nPCA_to_use{sub};
            load([out_dir, session_date, '_', mouse_ID, '_nPCA_', num2str(nPCA), '_post_process_outputs.mat']);
            load([out_dir, session_date, '_', mouse_ID, '_nPCA_', num2str(nPCA), '_Results.mat']);
            
            subMat = dir([behav_dir, '*', mouse_num, '-', session_date, '*']);
            
            load([behav_dir, subMat.name]);
            dest = out_dir;
            
            getTC_events;
            CuePair_2P_TC_quantification;
%             close all
%             clear
                         %end
        
    %end
end
end

