% main: process movie as following:
% 1. registeration -- find the most stable frame out of 100
% 2. PCA -- #PCA depends on how dense dendrite population, 100-200
% 3. ICA -- #IC is less than #PCA and with 30HZ, the temperal-spatial ratio
% is set to be ~1
% 4. Apply gaussian filter to ICA
% 5. Threshold ICA signal to get binary mask
% 6. Process mask (break multiple components on single layer to multiple,
% combine highly correlated ones and remove small ones)
% 7. Extrack timecourse (manually inspect TC and remove bad ones) and create a single layer mask for display

clear
file_info;

usFacs = 100;
behav_dir = 'Z:\Data\2P_imaging\behavior\';
for sub = [1] %size(mouseID,2) 
    for rID = 1
        file_info;
        out_dir  = fullfile('Z:', 'Analysis','Cue_reward_pairing_analysis','2P',[date{sub}, '_', runID{rID}, '_', mouseID{sub}],'\');
        
        [img, skip_run, img_fn] = loadFile(sub, rID);
        if length(size(img))==4;  %two channels were collected...
            img2 = img(2,:,:,:);   %red channel
            img= squeeze(img(1,:,:,:));  %green channel
        end
        
        if skip_run == 1
            continue
        end
        
        if ~exist(out_dir)
            mkdir(out_dir);
        end
        %
        % not cropping images anymore
        %             if exist([out_dir, 'ROI_xy.mat'],'file') == 2
        %                 load([out_dir, 'ROI_xy.mat']);
        %             else
        %                 [ROI_x, ROI_y] = get_2P_ROI(img); % get the ROI -  must be a rectangle
        %                 save([out_dir 'ROI_xy.mat'],  'ROI_x', 'ROI_y');
        %             end
        %             %
        %             img = img(ROI_x,ROI_y,:);
        %             img = img(:,:,22976:end);
        %% Motion registration
        %find the most stable frame among randomly selected frames and motion register to that frame
        [npw, nph, nt] = size(img);
        %if reg_out already exists use that instead
        if exist([out_dir, 'Reg_out.mat'],'file') == 2
            load([out_dir, 'Reg_out.mat']);
            [~,img_reg]=stackRegister_MA(img,img(:,:,1),usFacs,reg_out);
        else
            
            %restrict frame selection to frames with laser on
            [img_mat_file, laser_power_vec_ttl] = get_laser_power_data(sub, rID);
            if isempty(laser_power_vec_ttl)
                laser_power_vec_ttl = ones(1,nt);
            else
                laser_on_ind_conserv = conservative_laser_on(laser_power_vec_ttl);
            end
            laser_on_ind = find(laser_power_vec_ttl);   %frame numbers of frames with laser power on determined by _realtime ttl_log
            frame_nums_for_ref30 = laser_on_ind_conserv(randi([1,length(laser_on_ind_conserv)],1,30));
            frame_nums_for_samp100 = laser_on_ind_conserv(round(linspace(1,length(laser_on_ind_conserv))));
            
            %select 30 random frames from throughout the movie
            ref30 = img(:,:,frame_nums_for_ref30);
            
            %motion register each of the 30 random frames to 100 frames from the movie. Find the one with the lowerst dshift
            samp100 = img(:,:,frame_nums_for_samp100);
            dshift = [];
            for r = 1:size(ref30,3)
                [reg_out,aa] = stackRegister(samp100, ref30(:,:,r));
                dshift = [dshift;mean(((reg_out(:,3).^2)+(reg_out(:,4).^2)))];
            end
            
            %pick the frame which had the lowest dshift and motion register full movie to that frame
            min_f = find(dshift == min(dshift));
            img_ref = ref30(:,:,min_f);
            [reg_out, img_reg] = stackRegister(img, img_ref);
            if ~exist(out_dir)
                mkdir(out_dir);
            end
            save([out_dir, 'Reg_out.mat'], 'reg_out','img_ref', 'laser_on_ind');
            clear img img2
        end
        
        %% Compute principle components 
        nPCA = 500; %100 for old datasets, 500 for newer
        img_pca = img_reg(:,:,laser_on_ind_conserv); %only run PCA on frames with laser power on
        if size(img_pca, 3) > 65000 % downsample in time by 2 or 5
            img_pca = img_pca(:,:,1:2:end);
        end
        nf = size(img_pca,3);
        [mixedsig, mixedfilters, CovEvals, ~, ~, ...
            ~] = CellsortPCA2(img_pca,[1 nf], nPCA,[], out_dir, img_fn, []);
        PCuse = 1:nPCA;
        % to view PCs
        %             [PCuse] = CellsortChoosePCs(mixedfilters);
        
        %% Compute independent components
        %set variable values
        mu = 0.98; % spatial temporal mix
        nIC = 300;  %400- img90 100- img91
        termtol = 0.00001;
        maxrounds = 1000;
        
        %run ICA and save outputs
        [ica_sig, mixedfilters, ica_A, numiter] = CellsortICA(mixedsig, ...
            mixedfilters, CovEvals, PCuse, mu, nIC, [], termtol, maxrounds);
        icasig = permute(mixedfilters,[2,3,1]);
        save([out_dir, 'ICA.mat'], 'ica_sig', 'icasig');
        
        nIC = size(icasig,3);
        icasig = stackFilter(icasig);
        
        %% Create mask and process ICA outputs
        cluster_threshold = 95; %97- img90 97- img91   XXXXXXXXXXXX
        threshold = 0.8; % correlation threshold
        mask_cell = zeros(size(icasig));
        sm_logical = zeros(npw,nph);
        
        for ic = 1:nIC
            icasig(:,:,ic) = imclearborder(icasig(:,:,ic));
            sm_logical((icasig(:,:,ic)>mean([max(prctile(icasig(:,:,ic),cluster_threshold,1)) max(prctile(icasig(:,:,ic),cluster_threshold,2))])))=1;
            sm_logical((icasig(:,:,ic)<=mean([max(prctile(icasig(:,:,ic),cluster_threshold,1)) max(prctile(icasig(:,:,ic),cluster_threshold,2))])))=0;
            sm_logical = logical(sm_logical);
            mask_cell(:,:,ic) = bwlabel(sm_logical);
        end
        
        %                         mask_cell = zeros(size(icasig));
        %                         for ic = 1:nIC
        %                             bwimgcell = imCellEditInteractive(icasig(:,:,ic),[]);
        %                             mask_cell(:,:,ic) = bwlabel(bwimgcell);
        %                             close all
        %                         end
        mask_final = processMask(mask_cell);
        mask_raw = reshape(mask_final, npw, nph);
        figure; imagesc(mask_raw); truesize;
        
        [ ~, mask3D, ~] = finalMask(img_reg(:,:,1:10:end), mask_final, threshold, out_dir);
        
        %% Plotting TCs
        nmask = size(mask3D,3);
        FrameRate = 30;
        tc_avg = getTC(img_reg, mask3D, nmask);
        
        % check bad TC
        %             ICbad = [];
        %             saveData = 0;
        %             kl = floor(nmask/5);
        %             for k = 1:kl
        %                 plotTC(tc_avg, mask3D, reg_sum, (k-1)*5+1:k*5, FrameRate, out_dir,saveData);
        %                 ICbad_input = input('Number of bad IC ', 's');
        %                 ICbad = [ICbad str2num(ICbad_input)];
        %             end
        %             if mod(nmask,5) ~= 0
        %                 plotTC(tc_avg, mask3D, reg_sum, 5*kl+1:nmask, FrameRate, out_dir,saveData);
        %                 ICbad_input = input('Number of bad IC ', 's');
        %                 ICbad = [ICbad str2num(ICbad_input)];
        %             end
        %             close all
        %             tc_avg(:,ICbad) = []; mask3D(:,:,ICbad) = [];
        
        saveData = 1;
        reg_sum = sum(img_reg,3);
        %plot TCs with and without the frames without laser power
        plotTC(tc_avg, mask3D, reg_sum, 1:size(tc_avg,2), FrameRate, out_dir, saveData);
        plotTC(tc_avg(laser_on_ind,:), mask3D, reg_sum, 1:length(laser_on_ind), FrameRate, out_dir, 0);
        mask_flat = plotMask(mask3D, saveData, out_dir);
        
        data_corr = corrcoef(tc_avg);
        figure; fig = imagesc(data_corr);
        saveas(fig, [out_dir, 'data_corr.fig']);
        print([out_dir, 'data_corr.eps'],'-depsc')
        
        mask_final = processMask(mask3D);
        sz = [npw, nph];
        save([out_dir, 'Results.mat'], 'tc_avg', 'mask_raw', 'mask_flat', 'mask3D', 'mask_final', 'data_corr');
        save([out_dir, 'ROI_TCs.mat'],'tc_avg', 'mask_final', 'sz');
        close all
        
        mouse = mouseID{sub};
        dateID = date;
        date = dateID{sub};
        
        subMat = dir([behav_dir '*' mouse(end-2:end) '*' date '*']);
        
        load([behav_dir, subMat.name]);
        dest = out_dir;
        
        getTC_events;
        CuePair_2P_TC_quantification;
        close all
        %             clear
    end
end

% %             smwidth = 3;
% %             thresh = 6;
% %             arealims = 200;
% %             plotting = 1;
% %             weight = 0;
% %
% %             [ica_segments, segmentlabel, segcentroid] = CellsortSegmentation(mixedfilters, smwidth, thresh, arealims, plotting, weight);
% % %
% %             subtractmean = 0; normal = 0;
% %
% %             cell_sig = CellsortApplyFilter(img_reg(:,:,1:4:end), ica_segments, [], movm, subtractmean, normal);
% %
% %             [data_corr, sm, mask_flat] = finalMask2( cell_sig', ica_segments, threshold, out_dir);
% %
% % %             ica_seg = permute(sm, [3,1,2]);
% % %
% % %             cell_sig = CellsortApplyFilter(img_reg, ica_seg, [], movm, subtractmean);
% %
% %             deconvtau = 0;
% %             spike_thresh = 2;
% %             normalization = 0;
% %             dt = 0.1;
% %             [~, spt, spc] = CellsortFindspikes(cell_sig, spike_thresh, dt, deconvtau, normalization);
% %
% %             %% Show results
% %             ICbad = [];
% %             for i = 1:nIC/5
% %                 figure;
% %                 CellsortICAplot('contour', ica_seg, cell_sig, movm, [], dt, 1, 2, [(i-1)*5+1:i*5], spt, spc);
% %                 ICbad_input = input('Number of bad IC ', 's');
% %                 ICbad = [ICbad str2num(ICbad_input)];
% %             end
% %             close all
% %
% %             ica_seg(ICbad2,:,:) = []; cell_sig(ICbad2,:) = [];
% %             [~, spt, spc] = CellsortFindspikes(cell_sig, spike_thresh, dt, deconvtau, normalization);
% %             figure;CellsortICAplot('contour', ica_seg, cell_sig, movm, [], dt, 1, 1, 1:size(ica_seg,1) , spt, spc);
% %
% % %             save([out_dir, 'Results_new.mat'], 'cell_sig', 'ica_seg', 'mask_flat', 'data_corr');
% %
% %
% %
% %             cell_sig = cell_sig';
