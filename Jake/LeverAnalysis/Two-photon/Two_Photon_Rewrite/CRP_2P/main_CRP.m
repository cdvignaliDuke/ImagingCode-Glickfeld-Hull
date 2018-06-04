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
file_info_CRP;
usFacs = 10;
behav_dir = 'Z:\Data\2P_imaging\behavior\';
crp_dir = 'Z:\Analysis\Cue_reward_pairing_analysis\2P';
docombineMask = 0; doRerun = 0; useWarp = 0; doPCA = 0;
useGPU = 0; checkImg = 0;

for sub = 1 %[2,5, 9, 22, 23, 31, 3, 6, 10, 25, 27, 32]
    %     days_pair = days_all;
    if docombineMask == 1   %combines masks for pre and post learning days
        rID = 1;
        data_dir = fullfile('Z:\Data\2P_imaging',[dates{sub} '_' mouseID{sub}], mouseID{sub}, '\');
        out_dir  = fullfile(crp_dir, [dates{sub}, '_', runID{rID}, '_', mouseID{sub}],'\');
        %         data_dir = fullfile('Z:\home\jake\Data\2P_imaging',[days_pair{sub}(1:6) '_' days_pair{sub}(end-4:end)], days_pair{sub}(end-4:end), '\');
        %         out_dir = fullfile('Z:','home','jake','Analysis','Cue_reward_pairing_analysis','2P',[days_pair{sub}(1:6), '_', runID{rID}, '_', days_pair{sub}(end-4:end)], '\');
        %         while ~exist(out_dir)
        %             rID = rID + 1;
        %             out_dir = fullfile('Z:','home','jake','Analysis','Cue_reward_pairing_analysis','2P',[days_pair{sub}(1:6), '_', runID{rID}, '_', days_pair{sub}(end-4:end)], '\');
        %         end
        %%%%%%
        %         config_fn = dir(fullfile(data_dir,['*' runID{rID} '.mat']));
        %         load([data_dir,config_fn.name]);
        %         nframes = info.config.frames;
        %         img_fn   = dir(fullfile(data_dir,['*' runID{rID} '.sbx']));
        %         [~,img_fn,~] = fileparts(img_fn.name);
        %         [m,disp] = sbxalign_nonrigid([data_dir, img_fn],1:nframes);
        
        if checkImg == 1 %creates movies of the first and last 1000 frames motion registered to the start of imaging
            [img, ~, ~] = loadFile( data_dir, runID{rID}, [], 1000);
            [npw, nph, nframes1] = size(img);
            target = mean(img(:,:,500:600),3);
            [reg_out_f1000, f1000_reg] = stackRegister(img(:,:,1:1000), target);
            meanReg = mean(f1000_reg,3);
            save([out_dir, 'mean_1000_reg_first.mat'], 'meanReg', 'reg_out_f1000', 'f1000_reg');
            clear img
            %             fig=figure;
            %             imagesc(meanRandreg); truesize
            %             saveas(fig, [out_dir 'mean_1000_reg_first.fig']);
            %             print([out_dir 'mean_1000_reg_first.pdf'], '-dpdf');
            
            [img, ~, ~] = loadFile( data_dir, runID{rID}, 119001, 1000);
            [reg_out_l1000, l1000_reg] = stackRegister(img, target);
            meanReg = mean(l1000_reg,3);
            save([out_dir, 'mean_1000_reg_last.mat'], 'meanReg', 'reg_out_l1000', 'l1000_reg');
            clear img
            %             fig=figure;
            %             imagesc(meanRandreg); truesize
            %             saveas(fig, [out_dir 'mean_1000_reg_last.fig']);
            %             print([out_dir 'mean_1000_reg_last.pdf'], '-dpdf');
            
            Regfl = cat(3, f1000_reg, l1000_reg);
            writetiff(Regfl, [out_dir, 'checkFistLast1000Img.tiff']);
        elseif checkImg == 0
            [img_reg, laser_on_ind_conserv] = motionCorrect(data_dir, runID{rID}, out_dir);
            %             out_dir2 = fullfile('Z:','home','jake','Analysis','Cue_reward_pairing_analysis','2P',[days_pair{sub+1}(1:6), '_', runID{rID}, '_', days_pair{sub+1}(end-4:end)], '\');
            %             while ~exist(out_dir2)
            %                 rID = rID + 1;
            %                 out_dir2 = fullfile('Z:','home','jake','Analysis','Cue_reward_pairing_analysis','2P',[days_pair{sub+1}(1:6), '_', runID{rID}, '_', days_pair{sub+1}(end-4:end)], '\');
            %             end
            %             data_dir2 = fullfile('Z:\home\jake\Data\2P_imaging',[days_pair{sub+1}(1:6) '_' days_pair{sub+1}(end-4:end)], days_pair{sub+1}(end-4:end), '\');
            out_dir2  = fullfile(crp_dir, [dates{sub+1}, '_', runID{rID}, '_', mouseID{sub+1}],'\');
            data_dir2 = fullfile('Z:\Data\2P_imaging',[dates{sub+1} '_' mouseID{sub+1}], mouseID{sub+1}, '\');
            [img_reg2, laser_on_ind_conserv2] = motionCorrect(data_dir2, runID{rID}, out_dir2);
            [npw, nph, nframes2] = size(img_reg2);
            
            threshold = 0.8; % correlation threshold
            
            mask_final = reshape(labels, 1, npw*nph);
            [ ~, mask3D, ~] = finalMask(img_reg(:,:,1:4:end), mask_final, threshold, out_dir);
            [ ~, mask3D_2, ~] = finalMask(registered(:,:,1:4:end), mask_final, threshold, out_dir2);
            
            %% Plotting TCs
            nmask = size(mask3D,3);
            FrameRate = 30;
            saveData = 1;
            if skip_firstrun == 1
                tc_avg = getTC(img_reg2, mask3D, nmask);
                reg_sum = mean(img_reg2(:,:,2000:2100),3);
                plotTC(tc_avg, mask3D, reg_sum, 1:size(tc_avg,2), FrameRate, out_dir_comb, saveData);
            else
                reg_sum = sum(img_reg,3);
                nmask = size(mask3D,3);
                tc_avg = getTC(img_reg, mask3D, nmask);
                plotTC(tc_avg, mask3D, reg_sum, 1:size(tc_avg,2), FrameRate, out_dir, saveData);
                data_corr = corrcoef(tc_avg);
                mask_final = processMask_CRP(mask3D);
                %%%%%%%
                reg_sum = sum(registered,3);
                mask3D = mask3D_2;
                nmask = size(mask3D,3);
                tc_avg = getTC(registered, mask3D, nmask);
                plotTC(tc_avg, mask3D, reg_sum, 1:size(tc_avg,2), FrameRate, out_dir2, saveData);
                data_corr = corrcoef(tc_avg);
                mask_final = processMask_CRP(mask3D);
            end
            %         plotTC(tc_avg, mask3D, reg_sum, 1:size(tc_avg,2), FrameRate, out_dir_comb, saveData);
            
            mask_raw = reshape(mask_final, npw, nph);
            sz = [npw, nph];
            
            if skip_firstrun == 0
                save([out_dir_comb, 'Results.mat'], 'tc_avg', 'tc_avg2', 'mask_raw', 'mask_flat', 'mask3D', 'mask_final', 'reg_sum');
                save([out_dir_comb, 'ROI_TCs.mat'],'tc_avg', 'tc_avg2', 'mask_final', 'sz');
            else
                save([out_dir, 'Results.mat'], 'tc_avg', 'mask_raw', 'mask3D', 'mask_final', 'reg_sum');
                save([out_dir, 'ROI_TCs.mat'],'tc_avg', 'mask_final', 'sz');
            end
        end
        clearvars -except sub usFacs docombineMask doRerun useGPU doPCA skip_firstrun checkImg
        close all
        
    else  %if you are not combining masks from multiple days..=====================================================================================================
        for rID = 1
            %set pathnames
            usFacs = 100; % upsample factor for stackRegister
            out_dir =  fullfile(crp_dir, [dates{sub}, '_', runID{rID}, '_', mouseID{sub}], '\');
            data_dir = fullfile('Z:\Data\2P_imaging',[dates{sub} '_' mouseID{sub}], mouseID{sub}, '\');
            if ~exist(out_dir)
                mkdir(out_dir);
            end
          
            %% Extract pockel cell data
            if exist([data_dir, mouseID{sub}, '_', runID{rID}, '_', runID{rID}, '_realtime.mat'], 'file') ~= 2 %checks for data about pockel cell
                %set directories
                mouse = mouseID{sub};
                subMat = dir([behav_dir, '*', behavID{sub}, '-', dates{sub}, '*']); %sets directory as the behavior file 
                
                %load the behavior data
                load([behav_dir, subMat.name]);
                
                %identify times of "lever press" and trial end
                cStart = cell2mat(input.cLeverDown);
                cEnd = cell2mat(input.cTrialEnd);
                
                %load imaging matfile
                config_fn = dir(fullfile(data_dir,['*' runID{rID} '.mat']));
                load([data_dir, config_fn.name]);
                
                %for some reason Ziye overwrote ttl_log with this new version derived from the bx data. It makes no sense to I renamed it to ttl_log2 so it wont overwrite ttl_log 
                ttl_log2 = zeros(info.config.frames,1);
                for jj = 1:length(cStart)
                    ttl_log2(cStart(jj):cEnd(jj)) = 1;
                end
                save([data_dir, mouse, '_', runID{rID}, '_realtime.mat'], 'ttl_log2');  %overwrites existing ttl_log data
            end
            
            %Motion registration or check for exisint motion reg outputs
            if exist([out_dir, 'Reg_out.mat'],'file') == 2
                load([out_dir, 'Reg_out.mat']);
                [~,img_reg]=stackRegister_MA(img,img(:,:,1),usFacs,double(reg_out));
            elseif exist([out_dir, 'img_reg.mat'],'file') == 2
                load([out_dir, 'img_reg.mat']);
            else
                [img_reg, laser_on_ind_conserv, nt] = motionCorrect(data_dir, runID{rID}, out_dir); %uses stackRegister, automatically saves img_reg and Reg_out
            end
            
            if doRerun == 0 
                %% compute PRINCIPAL COMPONENTS 
                cd(data_dir); 
                [npw, nph, ~] = size(img_reg);
                nPCA = 100; %100 for old datasets, 500 for newer
                img_pca = img_reg(:,:,1:2:end); % downsample in time by 2 or 5
                nf = size(img_pca,3); 
                img_fn = [date, '_', runID{rID}, '_', mouseID{sub}];
                
                %run the PCA
                [mixedsig, mixedfilters, CovEvals, ~, ~, ~] = CellsortPCA_2P(img_pca,[1 nf], nPCA,[], out_dir, img_fn, []);
                %             [mixedsig, mixedfilters, CovEvals, ~, ~, ~] = CellsortPCA2(img_reg,[1 nFrames], nPCA, [], out_dir, img_fn, []);
                PCuse = 1:nPCA;
                
                % to view PCs
                %[PCuse] = CellsortChoosePCs_2P(mixedfilters);
                
                %% Compute independent components
                mu = 0.98; % spatial temporal mix
                nIC = 50;  %400- img90 100- img91
                termtol = 0.00001;
                maxrounds = 1000;
                
                %run ICA and save outputs
                [ica_sig, mixedfilters, ica_A, numiter] = CellsortICA_2P(mixedsig, mixedfilters, CovEvals, PCuse, mu, nIC, [], termtol, maxrounds);
                icasig = permute(mixedfilters,[2,3,1]);
                save([out_dir, 'ICA.mat'], 'ica_sig', 'icasig');
                
                %% Process the masks from the PCA/ICA
                %assign variables and allocate memory
                nIC = size(icasig,3);
                icasig = stackFilter(icasig);
                mask_cell = zeros(size(icasig));
                sm_logical = zeros(npw,nph);
                
                %set threshold a threshold for which pixels to include in a given dendrite's mask. 
                cluster_threshold = 97; %97- img90 97- img91 97- img070
                for ic = 1:nIC
                    icasig(:,:,ic) = imclearborder(icasig(:,:,ic));
                    sm_logical((icasig(:,:,ic)> mean([max(prctile(icasig(:,:,ic),cluster_threshold,1)) max(prctile(icasig(:,:,ic),cluster_threshold,2))])))=1;
                    sm_logical((icasig(:,:,ic)<=mean([max(prctile(icasig(:,:,ic),cluster_threshold,1)) max(prctile(icasig(:,:,ic),cluster_threshold,2))])))=0;
                    sm_logical = logical(sm_logical);
                    mask_cell(:,:,ic) = bwlabel(sm_logical); %run bwlabel to convert all mask %pixels to 1 and all others to 0
                    
                    % an older thresholding method? usese a semi manual interface to label cells. 
                    %                     bwimgcell = imCellEditInteractive(icasig(:,:,ic),[]);
                    %                     mask_cell(:,:,ic) = bwlabel(bwimgcell);
                    %                     close all
                end
                
                %splits or merges overlapping pixel masks, plots output
                mask_final = processMask_CRP(mask_cell);
                mask_raw = reshape(mask_final, npw, nph);
                figure;imagesc(mask_raw);truesize
                
                %consolidates highly correlated timecourses
                threshold = 0.8; % correlation threshold
                [ ~, mask3D, ~] = finalMask(img_reg(:,:,1:10:end), mask_final, threshold, out_dir);
                
                %% Plotting TCs
                nmask = size(mask3D,3);
                FrameRate = 30;
                tc_avg = getTC(img_reg, mask3D, nmask);
                
                saveData = 1;
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
            %             end
        end
    end
end


