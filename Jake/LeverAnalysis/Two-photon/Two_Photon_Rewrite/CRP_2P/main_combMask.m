% main: process movie as following:
% 1. registeration -- find the most stable frame out of 100
% 2. PCA -- #PCA depends on how dense dendrite population, 100-200
% 3. ICA -- #IC is less than #PCA and with 30HZ, the temperal-spatial ratio
% is set to be ~1
% 4. Apply gaussian filter to ICA
% 5. Threshold ICA signal to get binary mask
% 6. Process mask (break multiple components on single layer to multiple,
% combine highly correlated ones and remove small ones)
% 7. Extrack timecourse (manually inspect TC and remove bad ones) and create a singal layer mask for display

% clear
file_info_CRP;
usFacs = 100;
behav_dir = 'Z:\home\andrew\Behavior\Data\';
docombineMask = 0; doRerun = 1;
useGPU = 1;

for sub = [29, 33, 21,  15, 28] % : size(mouseID,2)
    if docombineMask == 1
        rID = 1;
        out_dir = fullfile('Z:','home','jake','Analysis','Cue_reward_pairing_analysis','2P',[days_pair{sub}(1:6), '_', runID{rID}, '_', days_pair_mouse{sub}], '\');
        while ~exist(out_dir)
            rID = rID + 1;
            out_dir = fullfile('Z:','home','jake','Analysis','Cue_reward_pairing_analysis','2P',[days_pair{sub}(1:6), '_', runID{rID}, '_', days_pair_mouse{sub}], '\');
        end
        
        data_dir = fullfile('Z:\home\jake\Data\2P_imaging',[days_pair{sub}(1:6) '_' days_pair_mouse{sub}], days_pair_mouse{sub}, '\');
        
        %%%%%%
%         config_fn = dir(fullfile(data_dir,['*' runID{rID} '.mat']));
%         load([data_dir,config_fn.name]);
%         nframes = info.config.frames;
%         img_fn   = dir(fullfile(data_dir,['*' runID{rID} '.sbx']));
%         [~,img_fn,~] = fileparts(img_fn.name);
%         [m,disp] = sbxalign_nonrigid([data_dir, img_fn],1:nframes);
        
        [img, ~, ~] = loadFile( data_dir, runID{rID});
        
        nframes1 = size(img,3);
        
        if exist([out_dir, 'Reg_out.mat'],'file') == 2
            out1 = load([out_dir, 'Reg_out.mat']);
            img_reg = [];
            if useGPU == 1
                for iF = 1 : floor(nframes1/1000)
                    [~,img_reg_temp]=stackRegister_MA(gpuArray(img(:,:,1 : 1000)),gpuArray(img(:,:,1)),usFacs,...
                        gpuArray(out1.reg_out(1 : 1000,:)));
                    img(:,:,1:1000) = [];
                    img_reg = cat(3, img_reg, double(gather(img_reg_temp)));
                    clear img_reg_temp
                end
                
            else
                [~,img_reg]=stackRegister_MA(img,img(:,:,1),usFacs,double(out1.reg_out));
            end
        end
        
        rID = 1;
        out_dir2 = fullfile('Z:','home','jake','Analysis','Cue_reward_pairing_analysis','2P',[days_pair{sub+1}(1:6), '_', runID{rID}, '_', days_pair_mouse{sub+1}], '\');
        while ~exist(out_dir)
            rID = rID + 1;
            out_dir2 = fullfile('Z:','home','jake','Analysis','Cue_reward_pairing_analysis','2P',[days_pair{sub+1}(1:6), '_', runID{rID}, '_', days_pair_mouse{sub+1}], '\');
        end
        data_dir2 = fullfile('Z:\home\jake\Data\2P_imaging',[days_pair{sub+1}(1:6) '_' days_pair_mouse{sub+1}], days_pair_mouse{sub+1}, '\');
        
        [img2, ~, ~] = loadFile( data_dir2, runID{rID});
        
        
        if exist([out_dir2, 'Reg_out.mat'],'file') == 2
            out2 = load([out_dir2, 'Reg_out.mat']);
            %             [~,img_reg2]=stackRegister_MA(img2,img2(:,:,1),usFacs,double(out2.reg_out));
        end
        
        img_reg2 = [];
        if useGPU == 1
            for iF = 1 : floor(size(img2,3)/1000)
                [~,img_reg_temp]=stackRegister_MA(gpuArray(img2(:,:,1 : 1000)),gpuArray(img2(:,:,1)),usFacs,...
                    gpuArray(out2.reg_out(1 : 1000,:)));
                img2(:,:,1:1000) = [];
                img_reg2 = cat(3, img_reg2, double(gather(img_reg_temp)));
                clear img_reg_temp
                %                     reset(gpuDevice(1))
            end
            
        else
            [~,img_reg2]=stackRegister_MA(img2,img2(:,:,1),usFacs,double(reg_out2));
        end
        
        [D, img_comb_ref] = imregdemons(out2.img_ref,out1.img_ref, [32 16 8 4],'AccumulatedFieldSmoothing',2.5,'PyramidLevels',4);
        
        [img_reg2] = imwarp(img_reg2, D);
        
        
    else
        for rID = 1:2
            file_info_CRP;
            
            usFacs = 100; % upsample factor
            out_dir = fullfile('Z:','home','jake','Analysis','Cue_reward_pairing_analysis','2P',[dates{sub}, '_', runID{rID}, '_', mouseID{sub}], '\', '_rerun\');
            out_dir1 = fullfile('Z:','home','jake','Analysis','Cue_reward_pairing_analysis','2P',[dates{sub}, '_', runID{rID}, '_', mouseID{sub}], '\');
            
            data_dir = fullfile('Z:\home\jake\Data\2P_imaging',[dates{sub} '_' mouseID{sub}], mouseID{sub}, '\');
            
            [img, skip_run, img_fn] = loadFile(data_dir, runID{rID});
            
            if skip_run == 1 || size(img,3) < 10000
                continue
            else
                if ~exist(out_dir)
                    mkdir(out_dir);
                end
                %%
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
                %%
                [npw, nph, nt] = size(img);
                
                if exist([out_dir, 'Reg_out.mat'],'file') == 2
                    load([out_dir, 'Reg_out.mat']);
                    [~,img_reg]=stackRegister_MA(img,img(:,:,1),usFacs,double(reg_out));
                else
                    % find most stable reference image among the 100 randomly
                    % selected images
                    ref30 = img(:,:,randi([1,nt],1,30));
                    
                    sf = 200*100+10;
                    samp100 = img(:,:,11:200:sf);
                    dshift = [];
                    for r = 1:size(ref30,3)
                        
                        [reg_out,aa] = stackRegister(samp100, ref30(:,:,r));
                        dshift = [dshift;mean(((reg_out(:,3).^2)+(reg_out(:,4).^2)))];
                        
                    end
                    
                    min_f = find(dshift == min(dshift));
                    img_ref = ref30(:,:,min_f);
                    %                 img_gpu = gpuArray(img(:,:,1:3000));
                    %                 img_ref_g = gpuArray(img_ref);
                    
                    if useGPU == 1
                        img_reg = []; reg_out = [];
                        for iF = 1 : floor(nt/1000)
                            [reg_out_temp,img_reg_temp]=stackRegister(gpuArray(img(:,:,1:1000)),gpuArray(img_ref));
                            img(:,:,1:1000) = [];
                            img_reg = cat(3, img_reg, double(gather(img_reg_temp)));
                            clear img_reg_temp
                            reg_out = [reg_out; gather(reg_out)];
                            %                     reset(gpuDevice(1))
                        end
                    else
                        [reg_out, img_reg] = stackRegister(img, img_ref);
                    end
                    
                    save([out_dir, 'Reg_out.mat'], 'reg_out','img_ref');
                    clear img
                end
                
                if doRerun == 0
                    nPCA = 500; %100 for old datasets, 500 for newer
                    img_pca = img_reg(:,:,1:2:end); % downsample in time by 2 or 5
                    nf = size(img_pca,3);
                    [mixedsig, mixedfilters, CovEvals, ~, ~, ...
                        ~] = CellsortPCA2(img_pca,[1 nf], nPCA,[], out_dir, img_fn, []);
                    %             [mixedsig, mixedfilters, CovEvals, ~, ~, ~] = CellsortPCA2(img_reg,[1 nFrames], nPCA, [], out_dir, img_fn, []);
                    
                    PCuse = 1:nPCA;
                    %             % to view PCs
                    %             %             [PCuse] = CellsortChoosePCs(mixedfilters);
                    %
                    %             %% Compute independent components
                    mu = 0.98; % spatial temporal mix
                    nIC = 200;  %400- img90 100- img91
                    termtol = 0.00001;
                    maxrounds = 1000;
                    
                    [ica_sig, mixedfilters, ica_A, numiter] = CellsortICA(mixedsig, ...
                        mixedfilters, CovEvals, PCuse, mu, nIC, [], termtol, maxrounds);
                    
                    icasig = permute(mixedfilters,[2,3,1]);
                    
                    save([out_dir, 'ICA.mat'], 'ica_sig', 'icasig');
                    
                    nIC = size(icasig,3);
                    
                    icasig = stackFilter(icasig);
                    
                    mask_cell = zeros(size(icasig));
                    sm_logical = zeros(npw,nph);
                    cluster_threshold = 97; %97- img90 97- img91
                    
                    for ic = 1:nIC
                        icasig(:,:,ic) = imclearborder(icasig(:,:,ic));
                        sm_logical((icasig(:,:,ic)>mean([max(prctile(icasig(:,:,ic),cluster_threshold,1)) max(prctile(icasig(:,:,ic),cluster_threshold,2))])))=1;
                        sm_logical((icasig(:,:,ic)<=mean([max(prctile(icasig(:,:,ic),cluster_threshold,1)) max(prctile(icasig(:,:,ic),cluster_threshold,2))])))=0;
                        sm_logical = logical(sm_logical);
                        mask_cell(:,:,ic) = bwlabel(sm_logical);
                    end
                    %
                    %             %                         mask_cell = zeros(size(icasig));
                    %             %                         for ic = 1:nIC
                    %             %                             bwimgcell = imCellEditInteractive(icasig(:,:,ic),[]);
                    %             %                             mask_cell(:,:,ic) = bwlabel(bwimgcell);
                    %             %                             close all
                    %             %                         end
                    mask_final = processMask(mask_cell);
                    
                    mask_raw = reshape(mask_final, npw, nph);
                    figure;imagesc(mask_raw);truesize
                    
                    threshold = 0.8; % correlation threshold
                    
                    [ ~, mask3D, ~] = finalMask(img_reg(:,:,1:10:end), mask_final, threshold, out_dir);
                    
                    %% Plotting TCs
                    nmask = size(mask3D,3);
                    FrameRate = 30;
                    tc_avg = getTC(img_reg, mask3D, nmask);
                    
                    % check bad TC
                    %             ICbad = [];
                    %             saveData = 0;
                    %             kl = floormask/5);
                    %             for k = 1:kl
                    %                 plotTC(tc_avg, mask3D, reg_sum, (k-1)*5+1:k*5, FrameRate, out_dir,saveData);
                    %                 ICbad_input = input('Number of bad IC ', 's');
                    %                 ICbad = [ICbad str2num(ICbad_input)];
                    %
                    %             end
                    %             if mod(nmask,5) ~= 0
                    %                plotTC(tc_avg, mask3D, reg_sum, 5*kl+1:nmask, FrameRate, out_dir,saveData);
                    %                ICbad_input = input('Number of bad IC ', 's');
                    %                ICbad = [ICbad str2num(ICbad_input)];
                    %             end
                    %             close all
                    tc_avg(:,ICbad) = []; mask3D(:,:,ICbad) = [];
                    saveData = 1;
                    reg_sum = sum(img_reg,3);
                    plotTC(tc_avg, mask3D, reg_sum, 1:size(tc_avg,2), FrameRate, out_dir, saveData);
                    mask_flat = plotMask(mask3D, saveData, out_dir);
                    %
                    data_corr = corrcoef(tc_avg);
                    figure; fig = imagesc(data_corr);
                    saveas(fig, [out_dir, 'data_corr.fig']);
                    print([out_dir, 'data_corr.eps'],'-depsc')
                    
                    mask_final = processMask(mask3D);
                    
                    %             out_dir2  = fullfile('Z:','home','jake','Analysis','Cue_reward_pairing_analysis','2P',[dates{sub}, '_', runID{rID}, '_', mouseID{sub}], '\');
                    %             load([out_dir2, 'Results.mat']);
                    %             nmask = size(mask3D,3);
                    %             FrameRate = 30;
                    %             tc_avg = getTC(img_reg, mask3D, nmask);
                else
                    load([out_dir1, 'Results.mat']);
                    nmask = size(mask3D,3);
                    FrameRate = 30;
                    tc_avg = getTC(img_reg, mask3D, nmask);
                    data_corr = corrcoef(tc_avg);
                end
                
                sz = [npw, nph];
                save([out_dir1, 'Results.mat'], 'tc_avg', 'mask_raw', 'mask_flat', 'mask3D', 'mask_final', 'data_corr');
                save([out_dir1, 'ROI_TCs.mat'],'tc_avg', 'mask_final', 'sz');
                save([out_dir, 'Results.mat'], 'tc_avg', 'mask_raw', 'mask_flat', 'mask3D', 'mask_final', 'data_corr');
                save([out_dir, 'ROI_TCs.mat'],'tc_avg', 'mask_final', 'sz');
                close all
                
                %             mouse = mouseID{sub};
                %             datesID = dates;
                %             dates = datesID{sub};
                %
                %             subMat = dir([behav_dir '*' '9' mouse(end-1:end) '*' dates '*']);
                %
                %             load([behav_dir, subMat.name]);
                %             dest = out_dir;
                %
                %             getTC_events;
                %             CuePair_2P_TC_quantification;
                %             close all
                %             clear
            end
        end
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
