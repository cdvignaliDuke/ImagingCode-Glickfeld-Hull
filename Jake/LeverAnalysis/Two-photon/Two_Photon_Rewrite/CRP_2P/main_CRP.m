% main: process movie as following:
% 1. registeration -- find the most stable frame out of 100
% 2. PCA -- #PCA depends on how dense dendrite population, 100-200
% 3. ICA -- #IC is less than #PCA and with 30HZ, the temperal-spatial ratio
% is set to be ~1
% 4. Apply gaussian filter to ICA
% 5. Threshold ICA signal to get binary mask
% 6. Process mask (break multiple components on single layer to multiple,
% combine highly correlated ones and remove small ones)
% 7. Extract timecourse (manually inspect TC and remove bad ones) and create a singal layer mask for display

clear
file_info_CRP2;
% file_info;  stable_int = {[1250:1500]}; dates = date;
usFacs = 100;
behav_dir = 'Z:\Data\2P_imaging\behavior\';
crp_dir = 'Z:\Analysis\Cue_reward_pairing_analysis\2P\FoV_mathcing_PCA_ICA';
docombineMask = 0; doRerun = 0; useWarp = 0; doPCA = 0;
useGPU = 0; checkImg = 0;

for sub =[1:2, 4:(length(stable_int)-2)]; %[2,5, 9, 22, 23, 31, 3, 6, 10, 25, 27, 32]
    %     days_pair = days_all;
    %if docombineMask == 1   %combines masks for pre and post learning days
    
    %         rID = 1;
    %         data_dir = fullfile('Z:\Data\2P_imaging',[dates{sub} '_' mouseID{sub}], mouseID{sub}, '\');
    %         out_dir  = fullfile(crp_dir, [dates{sub}, '_', runID{rID}, '_', mouseID{sub}],'\');
    %         %         data_dir = fullfile('Z:\home\jake\Data\2P_imaging',[days_pair{sub}(1:6) '_' days_pair{sub}(end-4:end)], days_pair{sub}(end-4:end), '\');
    %         %         out_dir = fullfile('Z:','home','jake','Analysis','Cue_reward_pairing_analysis','2P',[days_pair{sub}(1:6), '_', runID{rID}, '_', days_pair{sub}(end-4:end)], '\');
    %         %         while ~exist(out_dir)
    %         %             rID = rID + 1;
    %         %             out_dir = fullfile('Z:','home','jake','Analysis','Cue_reward_pairing_analysis','2P',[days_pair{sub}(1:6), '_', runID{rID}, '_', days_pair{sub}(end-4:end)], '\');
    %         %         end
    %         %%%%%%
    %         %         config_fn = dir(fullfile(data_dir,['*' runID{rID} '.mat']));
    %         %         load([data_dir,config_fn.name]);
    %         %         nframes = info.config.frames;
    %         %         img_fn   = dir(fullfile(data_dir,['*' runID{rID} '.sbx']));
    %         %         [~,img_fn,~] = fileparts(img_fn.name);
    %         %         [m,disp] = sbxalign_nonrigid([data_dir, img_fn],1:nframes);
    %
    %         if checkImg == 1 %creates movies of the first and last 1000 frames motion registered to the start of imaging
    %             [img, ~, ~] = loadFile( data_dir, runID{rID}, [], 1000);
    %             [npw, nph, nframes1] = size(img);
    %             target = mean(img(:,:,500:600),3);
    %             [reg_out_f1000, f1000_reg] = stackRegister(img(:,:,1:1000), target);
    %             meanReg = mean(f1000_reg,3);
    %             save([out_dir, 'mean_1000_reg_first.mat'], 'meanReg', 'reg_out_f1000', 'f1000_reg');
    %             clear img
    %             %             fig=figure;
    %             %             imagesc(meanRandreg); truesize
    %             %             saveas(fig, [out_dir 'mean_1000_reg_first.fig']);
    %             %             print([out_dir 'mean_1000_reg_first.pdf'], '-dpdf');
    %
    %             [img, ~, ~] = loadFile( data_dir, runID{rID}, 119001, 1000);
    %             [reg_out_l1000, l1000_reg] = stackRegister(img, target);
    %             meanReg = mean(l1000_reg,3);
    %             save([out_dir, 'mean_1000_reg_last.mat'], 'meanReg', 'reg_out_l1000', 'l1000_reg');
    %             clear img
    %             %             fig=figure;
    %             %             imagesc(meanRandreg); truesize
    %             %             saveas(fig, [out_dir 'mean_1000_reg_last.fig']);
    %             %             print([out_dir 'mean_1000_reg_last.pdf'], '-dpdf');
    %
    %             Regfl = cat(3, f1000_reg, l1000_reg);
    %             writetiff(Regfl, [out_dir, 'checkFistLast1000Img.tiff']);
    %         elseif checkImg == 0
    %             %[img_reg, laser_on_ind_conserv] = motionCorrect(data_dir, runID{rID}, out_dir);
    %             %mkdir(out_dir)
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
    %             img_ref = mean(img(:,:,stable_int{sub}));
    %             [reg_out, img_reg] = stackRegister(img(:,:,laser_on_ind_conserv), img_ref);
    %             %     img_reg = img_reg(:,:,laser_on_ind_conserv);
    %             save([out_dir, 'img_reg.mat'], 'img_reg', '-v7.3');
    %             save([out_dir, 'Reg_out.mat'], 'reg_out','img_ref', 'laser_on_ind_conserv', 'nt');
    %             %             out_dir2 = fullfile('Z:','home','jake','Analysis','Cue_reward_pairing_analysis','2P',[days_pair{sub+1}(1:6), '_', runID{rID}, '_', days_pair{sub+1}(end-4:end)], '\');
    %             %             while ~exist(out_dir2)
    %             %                 rID = rID + 1;
    %             %                 out_dir2 = fullfile('Z:','home','jake','Analysis','Cue_reward_pairing_analysis','2P',[days_pair{sub+1}(1:6), '_', runID{rID}, '_', days_pair{sub+1}(end-4:end)], '\');
    %             %             end
    %             %             data_dir2 = fullfile('Z:\home\jake\Data\2P_imaging',[days_pair{sub+1}(1:6) '_' days_pair{sub+1}(end-4:end)], days_pair{sub+1}(end-4:end), '\');
    %             out_dir2  = fullfile(crp_dir, [dates{sub+1}, '_', runID{rID}, '_', mouseID{sub+1}],'\');
    %             data_dir2 = fullfile('Z:\Data\2P_imaging',[dates{sub+1} '_' mouseID{sub+1}], mouseID{sub+1}, '\');
    %             [img_reg2, laser_on_ind_conserv2] = motionCorrect(data_dir2, runID{rID}, out_dir2);
    %             [npw, nph, nframes2] = size(img_reg2);
    %
    %             threshold = 0.8; % correlation threshold
    %
    %             mask_final = reshape(labels, 1, npw*nph);
    %             [ ~, mask3D, ~] = finalMask(img_reg(:,:,1:4:end), mask_final, threshold, out_dir);
    %             [ ~, mask3D_2, ~] = finalMask(registered(:,:,1:4:end), mask_final, threshold, out_dir2);
    %
    %             %% Plotting TCs
    %             nmask = size(mask3D,3);
    %             FrameRate = 30;
    %             saveData = 1;
    %             if skip_firstrun == 1
    %                 tc_avg = getTC(img_reg2, mask3D, nmask);
    %                 reg_sum = mean(img_reg2(:,:,2000:2100),3);
    %                 plotTC(tc_avg, mask3D, reg_sum, 1:size(tc_avg,2), FrameRate, out_dir_comb, saveData);
    %             else
    %                 reg_sum = sum(img_reg,3);
    %                 nmask = size(mask3D,3);
    %                 tc_avg = getTC(img_reg, mask3D, nmask);
    %                 plotTC(tc_avg, mask3D, reg_sum, 1:size(tc_avg,2), FrameRate, out_dir, saveData);
    %                 data_corr = corrcoef(tc_avg);
    %                 mask_final = processMask_CRP(mask3D);
    %                 %%%%%%%
    %                 reg_sum = sum(registered,3);
    %                 mask3D = mask3D_2;
    %                 nmask = size(mask3D,3);
    %                 tc_avg = getTC(registered, mask3D, nmask);
    %                 plotTC(tc_avg, mask3D, reg_sum, 1:size(tc_avg,2), FrameRate, out_dir2, saveData);
    %                 data_corr = corrcoef(tc_avg);
    %                 mask_final = processMask_CRP(mask3D);
    %             end
    %             %         plotTC(tc_avg, mask3D, reg_sum, 1:size(tc_avg,2), FrameRate, out_dir_comb, saveData);
    %
    %             mask_raw = reshape(mask_final, npw, nph);
    %             sz = [npw, nph];
    %
    %             if skip_firstrun == 0
    %                 save([out_dir_comb, 'Results.mat'], 'tc_avg', 'tc_avg2', 'mask_raw', 'mask_flat', 'mask3D', 'mask_final', 'reg_sum');
    %                 save([out_dir_comb, 'ROI_TCs.mat'],'tc_avg', 'tc_avg2', 'mask_final', 'sz');
    %             else
    %                 save([out_dir, 'Results.mat'], 'tc_avg', 'mask_raw', 'mask3D', 'mask_final', 'reg_sum');
    %                 save([out_dir, 'ROI_TCs.mat'],'tc_avg', 'mask_final', 'sz');
    %             end
    %         end
    %         clearvars -except sub usFacs docombineMask doRerun useGPU doPCA skip_firstrun checkImg
    %         close all
    %
    %else  %if you are not combining masks from multiple days..=====================================================================================================
    
    %set pathnames    
    data_dir = fullfile('Z:\Data\2P_imaging',[dates{sub} '_' mouseID{sub}], mouseID{sub}, '\');
    for rID = 1:2
        if exist([data_dir, mouseID{sub}, '_000_', runID{rID}, '.sbx'], 'file') == 2
            break
        end
    end
    out_dir =  fullfile(crp_dir, [dates{sub}, '_', runID{rID}, '_', mouseID{sub}], '\');
    data_dir = fullfile('Z:\Data\2P_imaging',[dates{sub} '_' mouseID{sub}], mouseID{sub}, '\');
    if ~exist(out_dir)
        mkdir(out_dir);
    end
            %           load([out_dir, 'img_reg.mat']);
            
%             %% Extract pockel cell data
%             
%             %set directories
%             mouse = mouseID{sub};
%             subMat = dir([behav_dir, '*', behavID{sub}, '-', dates{sub}, '*']); %sets directory as the behavior file
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
%             if exist([data_dir, mouseID{sub}, '_000_', runID{rID}, '_realtime.mat'], 'file') ~= 2 %checks for data about pockel cell
%                 %for some reason Ziye overwrote ttl_log with this new version derived from the bx data. It makes no sense to I renamed it to ttl_log2 so it wont overwrite ttl_log
%                 ttl_log2 = zeros(info.config.frames,1);
%                 for jj = 1:length(cStart)
%                     ttl_log2(cStart(jj):cEnd(jj)) = 1;
%                 end
%                 save([data_dir, mouse, '_', runID{rID}, '_realtime.mat'], 'ttl_log2');  %overwrites existing ttl_log data
%             end
%             
%             %% Motion registration or check for exisint motion reg outputs
%             %             if exist([out_dir, 'Reg_out.mat'],'file') == 2
% %                 load([out_dir, 'Reg_out.mat']);
% %                 [~,img_reg]=stackRegister_MA(img,img(:,:,1),usFacs,double(reg_out));
% %             elseif exist([out_dir, 'img_reg.mat'],'file') == 2
% %                 load([out_dir, 'img_reg.mat']);
% %             else
%                 [img_reg, laser_on_ind_conserv, nt] = motionCorrect(data_dir, runID{rID}, out_dir); %uses stackRegister, automatically saves img_reg and Reg_out
%                 config_fn = dir(fullfile(data_dir,['*' runID{rID} '.mat']));
%                 laser_power_fn = dir(fullfile(data_dir,['*' runID{rID} '_realtime.mat']));
%                 [img, ~, ~] = loadFile( data_dir, runID{rID}, [], []);
%                 [img_mat_file, laser_power_vec_ttl] = get_laser_power_data(data_dir, config_fn, laser_power_fn);
%                 if isempty(laser_power_vec_ttl)
%                     laser_on_ind_conserv = 1:nt;
%                 else
%                     laser_on_ind_conserv = conservative_laser_on(laser_power_vec_ttl);
%                 end
%                 laser_on_ind = find(laser_power_vec_ttl);
%                 img_ref = mean(img(:,:,stable_int{sub}),3);
%                 [reg_out, img_reg] = stackRegister(img(:,:,laser_on_ind_conserv), img_ref);
%                 save([out_dir, 'img_reg.mat'], 'img_reg', 'reg_out', 'img_ref', 'laser_on_ind_conserv', '-v7.3');  
%              %end
%             %writetiff(img_reg(:,:,1:500), [out_dir, 'img_reg_1_1500']);
%             
            if doRerun == 0 
                 %% compute PRINCIPAL COMPONENTS 
                cd(data_dir); 
%                 [npw, nph, ~] = size(img_reg);
                %nPCA = 400; %100 for old datasets, 500 for newer
                load([out_dir, 'img_reg']);
                for this_nPCA = 1:2
                    nPCA = all_nPCA{sub}(this_nPCA);
%                     %img_pca = img_reg(:,:,1:2:end); % downsample in time by 2 or 5
%                     img_pca = stackGroupProject(img_reg, 3);
%                     nf = size(img_pca,3);
%                     img_fn = [dates{sub}, '_', runID{rID}, '_', mouseID{sub}];
%                     
%                     %run the PCA
%                     [mixedsig_PCA, mixedfilters_PCA, CovEvals_PCA, ~, ~, ~] = CellsortPCA_2P(img_pca,[1 nf], nPCA,[], out_dir, img_fn, []);
%                     %             [mixedsig_PCA, mixedfilters, CovEvals, ~, ~, ~] = CellsortPCA2(img_reg,[1 nFrames], nPCA, [], out_dir, img_fn, []);
%                     PCuse = [1:nPCA]; %only use 85% of the PCs
%                     save([out_dir, 'PCA_variables_', num2str(nPCA),'.mat'], 'mixedsig_PCA', 'mixedfilters_PCA', 'CovEvals_PCA', 'nPCA');
%                     
%                     % to view PCs
%                     load([out_dir, 'PCA_variables_', num2str(nPCA),'.mat']);
%                     %writetiff(mixedfilters_PCA, [out_dir, 'PCuse_for_nPCA_', num2str(nPCA)]);
%                     
%                     %manually enter the PC #s to be used 
%                     if this_nPCA == 1
%                         PCuse = all_PCuse{sub};
%                     elseif this_nPCA ==2 
%                         PCuse = all_PCuse2{sub};
%                     end
%                     %[PCuse] = CellsortChoosePCs_2P(mixedfilters_PCA);
%                     
%                     %% Compute independent components
%                     %mus = [0.97, 0.7];
%                     %all_nIC = {[PCuse], [PCuse(1):(PCuse(end)-50)]};
%                     %for this_mu = mus
                     mu = 0.97; % spatial temporal mix
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
%                     writetiff(icasig, [out_dir, 'icasig_for_nPCA_', num2str(nPCA)]);
                    %figure; imagesc(sum(icasig,3));
                    %savefig([out_dir, 'icasig_sum']);
                    %continue
                    load([out_dir, 'ICA_variables_for_nPCA_', num2str(nPCA), '.mat']);
%                     load([out_dir, 'img_reg']);
                    nph=ICA_variables.nph; npw = ICA_variables.npw;
                    %% Process the masks from the PCA/ICA
                    
                    %select which ICs to keep based on morphology
%                     disp(['Beginning IC selection for ', dates{sub}, ' ', mouseID{sub}, ' nPCA=', num2str(nPCA)])
%                     IC_use = IC_manual_check(icasig);
%                     save([out_dir, 'ICA_variables_for_nPCA_', num2str(nPCA), '.mat'], 'IC_use', '-append');
%                     continue
                    
                    %stack filter - acts as a low pass spatial filter. It will filter out low spatial frequency noise such as small blips and bloops in the ICs which do not belong to dendrites.
                    icasig = icasig(:,:,[find(IC_use)]);
                    icasig2 = stackFilter(icasig);
                    %icasig2 = icasig;
                    %assign variables and allocate memory
                    nIC = size(icasig,3);
                    icasig3 = icasig2;
                    %                 mask_cell = zeros(size(icasig2));
                    %                 sm_logical = zeros(npw,nph);
                    
                    %set threshold a threshold for which pixels to include in a given dendrite's mask.
                    cluster_threshold = 96 % 90:100; %97- img90 97- img91 97- img070
                        mask_cell = zeros(size(icasig2));
                        sm_logical = zeros(npw,nph);
                        %bwimgcell = zeros(size(icasig2));
                        for ic = 1:nIC
                            %convert to a binary mask
                            icasig3(:,:,ic) = imclearborder(icasig2(:,:,ic));
                            sm_logical((icasig3(:,:,ic)> mean([max(prctile(icasig3(:,:,ic),cluster_threshold,1)) max(prctile(icasig3(:,:,ic),cluster_threshold,2))])))=1;
                            sm_logical((icasig3(:,:,ic)<=mean([max(prctile(icasig3(:,:,ic),cluster_threshold,1)) max(prctile(icasig3(:,:,ic),cluster_threshold,2))])))=0;
                            sm_logical = logical(sm_logical);
                            %bwlabel identifies unconnected objects within a single mask. So 0=background 1=Object#1 2=Object#2
                            if sum(sum(sm_logical,2),1) <51
                                sm_logical = 0;
                            end
                            mask_cell(:,:,ic) = bwlabel(sm_logical);
                            
                            % an older thresholding method? usese a semi manual interface to label cells.
                            %bwimgcell(:,:,ic) = imCellEditInteractive(icasig2(:,:,ic),[]);
                            %mask_cell(:,:,ic) = bwlabel(bwimgcell);
                            %close all
                        end
                        %close all
                        
                        %                     %save thresholded ICs
                        %                     IC_cut = find(squeeze(sum(sum(icasig3,2),1))'==0);
                        %                     icasig3(:,:,[IC_cut]) = [];
                        %                     IC_ind = [1:nIC];
                        %                     IC_ind(IC_cut) = [];
                        %                     figure; imagesc(sm_logical); title(['summed IC mask. Threshold = ', num2str(cluster_threshold)]);
                        %                     save([out_dir, 'automated_mask_thresh_', num2str(cluster_threshold)], 'sm_logical', 'icasig3');
                        %                             mask_final = processMask(mask_cell);
                        %                             mask_final = reshape(mask_final, npw, nph);
                        out_dir_spec = [out_dir, 'manual_PC_IC_checks\'];
                        if ~exist(out_dir_spec)
                            mkdir(out_dir_spec);
                        end
                        save([out_dir_spec, 'nPCA_', num2str(nPCA), '_mu_', num2str(mu), '_nIC_', num2str(nIC), ...
                            '_thresh_', num2str(cluster_threshold), '_mask_cell'], 'mask_cell')
                        %                             figure; imagesc(sum(mask_cell,3));
                        %                             title([dates{sub}, ' ', mouseID{sub}, ' mask post process stage 1']);
                        %                             savefig([out_dir_spec, 'nPCA_', num2str(nPCA), '_mu_', num2str(mu), '_nIC_', num2str(nIC), ...
                        %                                 '_thresh_', num2str(cluster_threshold), '.fig']);
                        
                        
                        %                 %manual mask edit
                        %                 output_mask = manual_mask_editor(bwimgcell, img_reg);
                        %                 save([out_dir, 'pix_mask_data'], 'output_mask', 'bwimgcell');
                        %                 plot_mask_contours(reshape(output_mask, npw,nph,size(output_mask,2)), []);
                        %                 savefig([out_dir, 'pixel mask after manual editing in manual_mask_editor'])
                        %
                        %                 figure; imagesc(sum(mask_cell,3)); title([ ' binary_ica_mask_automated nICs=', num2str(size(mask_cell,3))]);
                        %                 savefig([out_dir, 'binary_ica_mask_automated']);
                        %                 writetiff(mask_cell, [out_dir, 'binary_ica_mask_automated']);
                        %                 continue
                        
                        %splits or merges overlapping pixel masks, plots output
                        mask_cell(:,:, [find(sum(sum(mask_cell,2),1) < 51)] ) = [];
                        mask_final = processMask_CRP(mask_cell);
                        mask_raw = reshape(mask_final, npw, nph);
                        %                             figure; subplot(2,2,1); imagesc(max(img_reg(:,:,1:1500),[],3)); title([dates{sub}, ' ', mouseID{sub} 'motion registered max intensity 1500 frames']);
                        %                             subplot(2,2,2); imagesc(mask_raw); title('post thresholding')%truesize
                        
                        %consolidates highly correlated timecourses
                        threshold = 0.80; % correlation threshold
                        [ ~, mask3D_all, ~] = finalMask(img_reg(:,:,1:10:end), mask_final, threshold, out_dir);
                        [mask2D] = plotMask(mask3D_all,0,0,0);
                        %                             subplot(2,2,3); image(mask2D); title('after consolidating correlated masks');
                        
                        %                             %if two masks have overlapping pixels then combine them
                        %                             mask3D_all = reshape(mask3D_all, npw*nph, size(mask3D_all,3));
                        %                             for mask_num = 1:size(mask3D_all,2)
                        %                                 this_mask_ind = find(mask3D_all(:,mask_num));
                        %                                 if sum(sum( mask3D_all(this_mask_ind,[mask_num+1:end]) ,2),1) > 1 %check for masks with overlapping pixels. Exclude current mask since there will be 100% overlap
                        %                                     overlapping_masks = find(sum(mask3D_all(this_mask_ind,:),1)); %find indeces of the masks which overlap with this_mask
                        %                                     for this_overlapping_masks = overlapping_masks;
                        %                                         mask3D_all( find(mask3D_all(:,this_overlapping_masks)), mask_num )= 1; %add the overlapping mask to the current mask
                        %                                         mask3D_all(:,this_overlapping_masks) = 0;  %erase the overlapping mask from mask3D_all
                        %                                     end
                        %                                 end
                        %                             end
                        %                             mask3D_all = reshape(mask3D_all, npw, nph, size(mask3D_all,2));
                        
                        %Remove masks < 200 pixels from mask3D
                        too_small_ind = find( squeeze(sum(sum(mask3D_all(:,:,:),2),1))<200 );
                        mask3D = mask3D_all;
                        mask3D(:,:,too_small_ind) = [];
                        mask2D = plotMask(mask3D, 0,0, 0);
                        figure; image(mask2D);
                        title([dates{sub}, ' ', mouseID{sub}, ' final mask']);
                        savefig([out_dir_spec, 'nPCA_', num2str(nPCA), '_mu_', num2str(mu), '_nIC_', num2str(nIC), ...
                            '_thresh_', num2str(cluster_threshold), '_final', '.fig']);
                        save( [out_dir_spec, 'nPCA_', num2str(nPCA), '_mu_', num2str(mu), '_nIC_', num2str(nIC), ...
                            '_thresh_', num2str(cluster_threshold), '_mask3D.mat'], 'mask3D');
                        
                        mask_3D_buffed = make_mask_buffer(mask3D);
                        %
                        %                             too_small_ind = find( squeeze(sum(sum(mask_3D_buffed(:,:,:),2),1))<200 );
                        %                             mask_3D_buffed(:,:,too_small_ind) = [];
                        mask_2D_buffed = plotMask(mask_3D_buffed, 0,0, 0);
                        figure; image(mask_2D_buffed);
                        title([dates{sub}, ' ', mouseID{sub}, ' final mask with buffer']);
                        savefig([out_dir_spec, 'nPCA_', num2str(nPCA), '_mu_', num2str(mu), '_nIC_', num2str(nIC), ...
                            '_thresh_', num2str(cluster_threshold), '_final_buffer', '.fig']);
                        %
                   % end
                end
                % end
                continue
                            %% Plotting TCs
                            nmask = size(mask3D,3);
                            FrameRate = 30;
                            tc_avg = getTC(img_reg, mask3D, nmask);
                            
                            saveData = 0;
                            reg_sum = sum(img_reg,3);
                            plotTC(tc_avg, mask3D, reg_sum, 1:size(tc_avg,2), FrameRate, out_dir, saveData);
                            mask_flat = plotMask(mask3D, saveData, out_dir);

                %
                data_corr = corrcoef(tc_avg);
                %                     figure; fig = imagesc(data_corr);
                %                     saveas(fig, [out_dir, 'data_corr.fig']);
                %                     print([out_dir, 'data_corr.eps'],'-depsc')
                
                mask_final = processMask_CRP(mask3D);
                sz = [npw, nph];
                save([out_dir, 'Results.mat'], 'tc_avg', 'mask_raw', 'mask_flat', 'mask3D', 'mask_final', 'data_corr');
                save([out_dir, 'ROI_TCs.mat'],'tc_avg', 'mask_final', 'sz');
                
                %             out_dir2  = fullfile('Z:','home','jake','Analysis','Cue_reward_pairing_analysis','2P',[dates{sub}, '_', runID{rID}, '_', mouseID{sub}], '\');
                %             load([out_dir2, 'Results.mat']);
                %             nmask = size(mask3D,3);
                %             FrameRate = 30;
                %             tc_avg = getTC(img_reg, mask3D, nmask);
            else
                load([out_dir, 'Results.mat']);
                nmask = size(mask3D,3);
                FrameRate = 30;
                tc_avg = getTC(img_reg, mask3D, nmask);
                data_corr = corrcoef(tc_avg);
            end
            
            continue
            tc_temp = zeros(nt, nmask);
            tc_temp(laser_on_ind_conserv,:) = tc_avg;
            tc_avg = tc_temp;
            sz = [npw, nph];
            
            save([out_dir, 'Results.mat'], 'tc_avg', 'mask_raw', 'mask_flat', 'mask3D', 'mask_final', 'data_corr');
            save([out_dir, 'ROI_TCs.mat'],'tc_avg', 'mask_final', 'sz');
            %clearvars -except sub usFacs docombineMask doRerun useGPU outdir
            %close all
            
            mouse = mouseID{sub};
            datesID = dates;
            dates = datesID{sub};
            
            subMat = dir([behav_dir, '*', behavID{sub}, '-', dates(1:6), '*']);
            
            load([behav_dir, subMat.name]);
            dest = out_dir;
            
            getTC_events;
            CuePair_2P_TC_quantification;
%             close all
%             clear
                         %end
        
    %end
end


