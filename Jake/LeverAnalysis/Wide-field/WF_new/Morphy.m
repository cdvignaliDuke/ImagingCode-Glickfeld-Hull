% I=double(imgMasked(:,:,245));
% hy = fspecial('sobel');
% hx = hy';
% Iy = imfilter(I, hy, 'replicate');
% Ix = imfilter(I, hx, 'replicate');
% gradmag = sqrt(Ix.^2 + Iy.^2);
% figure
% imshow(gradmag,[]), title('Gradient magnitude (gradmag)');truesize
%
% sim = gradmag;
% % Perform Sobel Edge detection
% a = edge(gradmag, 'sobel');
% % Find the mask of the object
% mask = gradmag > 1.1*mean(mean(gradmag(a)));
%
% % Shrink the object
% se = strel('square', 5);
% mask_s = imerode(mask, se);
%
% % Remove the inner boundary of the object
% mask(mask_s) = false;
%
% % Slightly enlarge now to ensure outer boundary is removed
% mask = imdilate(mask, se);
%
% % Create new image by removing the boundaries of the
% % edge image
% sim(mask) = 0;
%
% [Gx, Gy] = imgradientxy(I);
% Gy(~sim) = 0;
%
% sim2 = stackFilter(-Gy,1.5);
% % Show the result
% figure; imagesc(sim2);
%
% Ithreshold = 93;
% sm_logical = zeros(size(gradmag));
% sm_logical (sim2 > mean([max(prctile(sim2,Ithreshold,1)) max(prctile(sim2,Ithreshold,2))])) = 1;
% figure;imagesc(sm_logical)
clear
file_info;


for sub = 3:size(mouseID,2)
    for rID = 1
        file_info_WF;
        usFacs = 100;
        PCA_dir  = fullfile('Z:','home','jake','Analysis','WF Lever Analysis','PCA_ICA_output',date{sub},'\');
        WF_dir = fullfile('Z:','home','jake','Analysis','WF Lever Analysis','Meta-kmeans_output_dir', date{sub}, '\');
        out_dir = fullfile('Z:','home','jake','Analysis','WF Lever Analysis','Morphology_output',date{sub},'\');
        %         [img, skip_run, img_fn] = loadFile(sub, rID);
        img_dir = fullfile('Z:','home','jake','Data','WidefieldImaging','GCaMP', date{sub}, '\');
        
        img = readtiff(img_dir);
        
        if ~exist(out_dir)
            mkdir(out_dir);
        end
        
        if exist([out_dir, 'ROI_xy.mat'],'file') == 2
            load([out_dir, 'ROI_xy.mat']);
        elseif exist([WF_dir, date{sub}, 'reg_out.mat'], 'file') == 2
            load([WF_dir, date{sub}, 'reg_out.mat']);
        end
        
        img = img(ROI_x,ROI_y,:);
        
        load([PCA_dir, 'Results.mat'], 'bkMask');
        
        imgMasked = bsxfun(@times, img, cast(bkMask, 'like', img));
        
        figure;imagesc(imgMasked(:,:,1));
        [x,y,nt] = size(imgMasked);
        h = impoly;
        comp_mask = (createMask(h));
        img1 = bsxfun(@times, imgMasked, cast(comp_mask, 'like', imgMasked));
        
        mask1 = getCells(img1, 95);
        
        img1 = bsxfun(@times, imgMasked, cast(~comp_mask, 'like', imgMasked));
        
        mask2 = getCells(img1, 95);
        
        mask_final = bwlabel(mask1 + mask2);
        
        mask_final = reshape(mask_final, 1, x*y);
        
        threshold = 0.85;
        
        [ ~, mask3D, ~] = finalMask(imgMasked(:,:,1:nt), mask_final, threshold, out_dir);
        
        %% Plotting TCs
        nmask = size(mask3D,3);
        FrameRate = 100;
        tc_avg = getTC(imgMasked, mask3D, nmask);
        
        saveData = 1;
        plotTC(tc_avg, mask3D, [], 1:size(tc_avg,2), FrameRate, out_dir, saveData);
        mask_flat = plotMask(mask3D, saveData, out_dir);
        
        data_corr = corrcoef(tc_avg);
        figure; fig = imagesc(data_corr);
        saveas(fig, [out_dir, 'data_corr.fig']);
        print([out_dir, 'data_corr.eps'],'-depsc')
        
        mask_final = processMask(mask3D);
        sz = [x, y];
        save([out_dir, 'Results.mat'], 'tc_avg', 'mask_flat', 'mask3D', 'mask_final', 'data_corr', 'bkMask');
        save([out_dir, 'ROI_TCs.mat'],'tc_avg', 'mask_final', 'sz');
        close all
        
      
        clear
        
    end
end
% bwimgcell = imCellEditInteractive(gradmag,[]);
% bwimgcell = imCellEditInteractive(I,[]);

% overlay segmentation on original image
% reg_mean = zeros(1002,1004);
% reg_mean(ROI_x,ROI_y) = imgMasked(:,:,245);
% c = zeros(1002,1004);
% c(ROI_x, ROI_y) = cellMask;
% figure;imshow(mat2gray(reg_mean))
% hold on;contour(c, 'Color', [1 1 1])



