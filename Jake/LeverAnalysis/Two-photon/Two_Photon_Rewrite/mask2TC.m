function mask2TC(sz, mask_cell, img_reg)
%consolidates all ROIs within IC into single ROI
thresh = 0.8; %correlation threshold for calling two dendrites one thing
nIC  = size(mask_cell,3);
mask_cell_temp = zeros(sz(1)*sz(2), nIC);
for ic = sel
    if length(unique(reshape(mask_cell(:,:,ic),[1 sz(1)*sz(2)])))>2
        data_tc_temp = stackGetTimeCourses(img_reg,mask_cell(:,:,ic));
        data_corr_temp = corrcoef(data_tc_temp);
        ind_rem = 1:length(unique(reshape(mask_cell(:,:,ic),[1 sz(1)*sz(2)])))-1;
        for i = 1:length(unique(reshape(mask_cell(:,:,ic),[1 sz(1)*sz(2)])))-1
            ind = ind_rem(find(data_corr_temp(min(ind_rem,[],2),ind_rem)>thresh));
            if length(ind)>1
                for ii = ind
                    if i == 1
                        mask_cell_temp(find(mask_cell(:,:,ic)== ii),ic) = 1;
                    else
                        mask_cell_temp(find(mask_cell(:,:,ic)== ii),nIC) = 1;
                    end
                end
            else
                if i == 1
                    mask_cell_temp(find(mask_cell(:,:,ic)== i),ic) = 1;
                else
                    mask_cell_temp(find(mask_cell(:,:,ic)== ind),nIC) = 1;
                end
            end
            ind_rem = ind_rem(~ismember(ind_rem,ind));
            if ~isempty(ind_rem)
                cat(3, mask_cell_temp, zeros(size(mask_cell_temp(:,:,1))));
                nIC = nIC+1;
            else
                break
            end
        end
    else
        mask_cell_temp(find(mask_cell(:,:,ic)),ic) = 1;
    end
end

%get preliminary timecourses for segregating and grouping ROIs
data_tc = zeros(size(img_down,3), nIC);
for ic = sel;
    if sum(mask_cell_temp(:,ic),1)>0
        data_tc(:,ic) = stackGetTimeCourses(img_reg, reshape(mask_cell_temp(:,ic), [sz(1) sz(2)]));
    end
end
data_corr = corrcoef(data_tc);
%figure; imagesc(data_corr)

%consolidate timecourses that are highly correlated
[i j] = find(and(data_corr>thresh,data_corr<1));
n = size(i,1);
if n>1
    for ii = 1:(n/2)
        ind = find(mask_cell_temp(:,i(ii)));
        mask_cell_temp(ind,j(ii)) = 1;
        mask_cell_temp(ind,i(ii)) = 0;
    end
end

%finds overlapping pixels of ROIs and based on correlations decides whether
%to group them or to split them- if splitting, then overlapping pixels are
%eliminated from both ROIs

mask_overlap = zeros(1,sz(1)*sz(2));
mask_all = zeros(1,sz(1)*sz(2));
count = 0;
for ic = 1:nIC
    ind_new = find(mask_cell_temp(:,ic))';
    if length(ind_new)>1
        ind_old = find(mask_all);
        overlap = ismember(ind_old,ind_new);
        ind_both = find(overlap);
        if length(ind_both)>1
            ic_match = unique(mask_all(ind_old(ind_both)));
            for im = 1:length(ic_match)
                if data_corr(ic, ic_match(im))> thresh
                    count = count+1;
                    mask_all(ind_new) = ic_match(im);
                else
                    mask_all(ind_new) = ic;
                    mask_all(ind_old(ind_both)) = 0;
                    mask_overlap(ind_old(ind_both)) = 1;
                end
            end
        else
            mask_all(ind_new) = ic;
        end
    end
end
figure; imagesc(reshape(mask_all,[sz(1) sz(2)]))


% removes ICs smaller than 200 pixels, renumbers ROIs so in continuous ascending order
start = 1;
mask_final = zeros(size(mask_all));
for ic = 1:max(mask_all,[],2)
    ind = find(mask_all==ic);
    if length(ind)<200
        mask_overlap(find(mask_all==ic)) = 1;
        mask_all(ind) = 0;
    end
    ind = find(mask_all==ic);
    if length(ind)>0
        mask_final(ind)=start;
        start= start+1;
    end
end

figure; imagesc(reshape(mask_final,[sz(1) sz(2)]))

% print([dest '_mask_final.eps'], '-depsc');
% print([dest '_mask_final.pdf'], '-dpdf');





%% 5. Extract timecourses
data_tc = stackGetTimeCourses(img, reshape(mask_final,[sz(1) sz(2)]));
data_corr = corrcoef(data_tc);
sz = size(mask_cell);
clear input

%use this graph to eliminate specific TCs based on there plots. Need to
%edit data_tc, mask_final
shift=0;
figure;
for i = 1:size(data_tc,2)
    plot(data_tc(:,i)+shift); hold on
    if ismember(i, [5,10,15,20,25,30,35,40,45,50])
        hline(shift);
    end
    shift = shift+6000;
end