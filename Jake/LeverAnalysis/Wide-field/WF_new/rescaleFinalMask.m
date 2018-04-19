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

clear
file_info_WF;


for sub = 6%size(mouseID,2)
    sub
    for rID = 1
        file_info_WF;
        usFacs = 100;
        out_dir  = fullfile('Z:','home','jake','Analysis','WF Lever Analysis','PCA_ICA_output',date{sub},'\');
        WF_dir = fullfile('Z:','home','jake','Analysis','WF Lever Analysis','Meta-kmeans_output_dir', date{sub}, '\');
        %         [img, skip_run, img_fn] = loadFile(sub, rID);
        
        
        if exist([out_dir, 'ROI_xy.mat'],'file') == 2
            load([out_dir, 'ROI_xy.mat']);
        elseif exist([WF_dir, date{sub}, 'reg_out.mat'], 'file') == 2
            load([WF_dir, date{sub}, 'reg_out.mat']);
        else
            [ROI_x, ROI_y] = get_2P_ROI(img); % get the ROI -  must be a rectangle
            save([out_dir 'ROI_xy.mat'],  'ROI_x', 'ROI_y');
        end

        load([out_dir, 'Results.mat'], 'tc_avg', 'mask_raw', 'mask_flat', 'mask3D', 'mask_final', 'data_corr', 'bkMask');
        load([out_dir, 'ROI_TCs.mat'],'tc_avg', 'mask_final', 'sz');
        
        img = zeros(1002, 1004);
        img(ROI_x, ROI_y) = reshape(mask_final, sz(1), sz(2));
        mask_final = reshape(img, 1, 1002*1004);
        sz = [1002,1004];
        save([out_dir, 'Results.mat'], 'tc_avg', 'mask_raw', 'mask_flat', 'mask3D', 'mask_final', 'data_corr', 'bkMask');
        save([out_dir, 'ROI_TCs.mat'],'tc_avg', 'mask_final', 'sz');
        close all
        
        
        %
        %             getTC_events;
        %             HAD_2P_TC_quantification;
%         close all
        clear
        
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
