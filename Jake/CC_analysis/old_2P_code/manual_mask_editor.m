function output_mask = manual_mask_editor(bwimgcell, img_reg);
% script meant to provide an environment where one can pseudo-manually
% remove overlapping ROIs. 

%define colors for plotting 
color_spec = repmat([1,0.4,0.7;   1,1,0;   0,1,0;  1,0,1;   1,0.5,0;   0.2,0.6,1,;  1,1,1], 100, 1);  %1,0,0;   excluding red bc mask of interest is plotted in red

%remove masks with too few pixels 
rm_mask_ind = find( squeeze(sum(sum(bwimgcell,2),1)) < 50);
bwimgcell(:,:,[rm_mask_ind]) = [];

%identify overlapping pixels
[nph, npw, nIC] = size(bwimgcell);
bwimgcell = reshape(bwimgcell, nph*npw, nIC); %convert to 2D matrix
overlap_ind = sum(bwimgcell,2);
overlap_ind(find(overlap_ind < 2)) = 0;
overlap_ind = logical(overlap_ind); %1 for each pixel belonging to more than one mask. 0 for all others

%reshape movie
img_reg = reshape(img_reg, nph*npw, size(img_reg,3));

%sort ROIs by the fraction of their pixels which are overlapping with other ROIs
each_mask_size = squeeze(sum(bwimgcell,1))';
overlap_ind_all = bwimgcell.*repmat(overlap_ind,1,nIC); 
overlap_size = sum(overlap_ind_all,1)';
overlap_ratio = overlap_size./each_mask_size;
[~, ratio_sort_ind] = sort(overlap_ratio, 1, 'descend');
bwimgcell = bwimgcell(:,[ratio_sort_ind]);
overlap_ind_all = overlap_ind_all(:,ratio_sort_ind);
overlap_pix = sum(sum(overlap_ind_all,2),1);

while overlap_pix > 0
    %serially select masks of interest.
    overlap_sum = sum(overlap_ind_all,1);
    this_mask = find(overlap_sum > 0,1,'first');
    
    %find masks which overlap with this one. If none then continue
    overlap_this_mask = repmat(bwimgcell(:,this_mask),1,size(bwimgcell,2))  .*  bwimgcell; 
    overlap_this_mask = find(sum(overlap_this_mask,1));
    overlap_this_mask = overlap_this_mask(overlap_this_mask~=this_mask);
    if isempty(overlap_this_mask)
        continue
    end
    
    %plot mask of interest in red. Plot overlapping ROIs in various colors with
    %numbers. Plot all other ROIs in white.
    figure('rend', 'painters', 'pos', [1900 550 1250 500]); imagesc(zeros(nph,npw)); hold on;
    for mask_num = 1:size(overlap_ind_all,2)  %all masks
        contour( reshape(bwimgcell(:,mask_num), nph,npw), 'Color', [1,1,1]);
    end
    color_spec_ind = 1;
    for mask_num = overlap_this_mask %overlapping masks
        contour( reshape(bwimgcell(:,mask_num), nph,npw), 'Color', color_spec(color_spec_ind,:));
        color_spec_ind = color_spec_ind+1;
    end
    contour(reshape(bwimgcell(:,this_mask), nph,npw), 'Color', [1,0,0]); %plot mask of interest
    title(['contour map: red = mask of interest. Colors = overlapping masks']);
    
    %extract from mask of interst and neighbors
    this_mask_TC = mean( img_reg( find(bwimgcell(:,this_mask)) ,:) ,1);
    neighbor_TCs = NaN(length(overlap_this_mask), size(img_reg,2));
    for mask_num = 1:length(overlap_this_mask)
        neighbor_TCs(mask_num,:) = mean( img_reg( find(bwimgcell(:,overlap_this_mask(mask_num))) ,:) ,1);
    end
    
    %plot TCs
    figure('rend', 'painters', 'pos', [1900 30 1250 500]); 
    plot(this_mask_TC([1:8000]), 'Color', 'r'); hold on;
    offset_val = diff([min(min(neighbor_TCs)), max(max(neighbor_TCs))])/2;
    color_spec_ind = 1;
    for mask_num = 1:size(neighbor_TCs,1)
        plot(neighbor_TCs(mask_num,[1:8000]) - (offset_val*mask_num), 'Color', color_spec(color_spec_ind,:)); 
        color_spec_ind = color_spec_ind+1;
    end
    
    %calculate correltations between that mask and its neighbors
    activity_corrs = NaN(1,size(neighbor_TCs,1));
    for mask_num = 1:size(neighbor_TCs,1)
        activity_corrs(mask_num) = corr(this_mask_TC', neighbor_TCs(mask_num,:)');
    end
    title(['temporal correlations between the mask of interest and each of its neighbors: ', num2str(activity_corrs)]);
    next_mask_val = 0;
    
    %manually select which rois to
    for mask_num = 1:length(overlap_this_mask)
        %check for exit condition
        if next_mask_val == 1;
            break
        end
        next_mask_val = 0;
        
        %find #of overlapping pixels
        this_ovelap_num = bwimgcell(:,overlap_this_mask(mask_num)).*bwimgcell(:,this_mask);
        disp(['Overlapping neighbor #', num2str(mask_num)]); 
        disp(['Color: ', num2str(color_spec(mask_num,:))]); 
        disp(['Number of overlapping pixels=', num2str(sum(this_ovelap_num))]); 
        disp(['Correlation=', num2str(activity_corrs(mask_num))]);
        
        release_val=0;
        while release_val == 0
            out_str = input(['For neighboring mask #', num2str(mask_num), ' select an option: cut, combine, or keep:   '], 's');
            % 1) cut
            if strcmp(out_str, 'cut')
                release_val =1;
                bwimgcell(:,overlap_this_mask(mask_num)) = 0;
            % 2) combine
            elseif strcmp(out_str, 'com')
                release_val = 1;
                bwimgcell(:,this_mask) = logical( bwimgcell(:,this_mask) + bwimgcell(:,overlap_this_mask(mask_num)));
                bwimgcell(:,overlap_this_mask(mask_num)) = 0;
            % 3) keep both but remove overlapping pixels
            elseif strcmp(out_str, 'keep')
                release_val = 1;
                overlap_region =  find( (bwimgcell(:,this_mask) + bwimgcell(:,overlap_this_mask(mask_num))) == 2);
                bwimgcell(overlap_region,this_mask) =0;
                bwimgcell(overlap_region,overlap_this_mask(mask_num)) = 0;
            % 4) ignore 
            elseif strcmp(out_str, 'ign')
                release_val=1;
            % 5) cut mask of interest and DELETE all its pixels from other  masks
            elseif strcmp(out_str, '5')
                release_val=1;
                bwimgcell( [find(bwimgcell(:, this_mask))],: ) = 0;
                next_mask_val = 1;
            % 6) cut mask of interest and KEEP pixels in other masks 
            elseif strcmp(out_str, '6')
                 bwimgcell(:,this_mask) = 0;
                 release_val=1;
                 next_mask_val = 1;
            end
        end
    end
     
    %update overlap_ind 
    overlap_ind = sum(bwimgcell,2);
    overlap_ind(find(overlap_ind < 2)) = 0;
    overlap_ind = logical(overlap_ind); %1 for each pixel belonging to more than one mask. 0 for all others
    %update overlap_ind_all
    each_mask_size = squeeze(sum(bwimgcell,1))';
    overlap_ind_all = bwimgcell.*repmat(overlap_ind,1,nIC);
    overlap_size = sum(overlap_ind_all,1)';
    overlap_ratio = overlap_size./each_mask_size;
    [~, ratio_sort_ind] = sort(overlap_ratio, 1, 'descend');
    bwimgcell = bwimgcell(:,[ratio_sort_ind]);
    overlap_ind_all = overlap_ind_all(:,ratio_sort_ind);
    %update overlap_pix
    overlap_pix = sum(sum(overlap_ind_all,2),1);
    
    close all
end

output_mask = bwimgcell;
end


