function [ data_corr, sm, sm_dsum] = finalMask(img_reg,mask_final, threshold, out_dir)
nmask = max(max(mask_final));
[npw, nph, ~] = size(img_reg);
%nmask = size(sel,2);
sm = zeros(npw,nph,nmask);

for ii=1:nmask
    sm_placeholder = zeros(npw,nph);
    sm_placeholder(mask_final==ii) = 1;
    sm(:,:,ii)= sm_placeholder;
end

tc_avg = getTC(img_reg, sm, nmask);

data_corr = corrcoef(tc_avg);
%  figure; imagesc(data_corr);

%consolidate timecourses that are highly correlated
[i, j] = find(and(data_corr>threshold,data_corr<1));
n = size(i,1);
sm_temp = reshape(sm,npw*nph,nmask);
remove_2 =[];
if n>1
    comp = [j(1); i(j == j(1))];
    total = unique(i);
    rest_ind = ismember(total, comp);
    rest_comp = total(~rest_ind);
    for k = 2:length(comp)
        sm_temp(find(sm(:,:,comp(k))),comp(1)) = 1;
        sm(:, :, comp(1)) = reshape(sm_temp(:,comp(1)), npw ,nph);
    end
    if length(rest_comp) == 2
        sm_temp(find(sm(:,:,rest_comp(2))),rest_comp(1)) = 1;
        sm(:, :, rest_comp(1)) = reshape(sm_temp(:,rest_comp(1)), npw ,nph);
        remove_2 = rest_comp(2);
    elseif length(rest_comp) > 2
        for k = 1:length(rest_comp)
            if k > length(rest_comp)
                break
            else
                rest_comp1 = rest_comp(k);
                rest_comp2 = j(i == rest_comp1);
                rest_comp(ismember(rest_comp,rest_comp2))=[];
                for kk = 1:length(rest_comp2)
                    sm_temp(find(sm(:,:,rest_comp2(kk))),rest_comp1) = 1;
                    sm(:, :, rest_comp1) = reshape(sm_temp(:,rest_comp1), npw ,nph);
                end
                remove_2 = [remove_2; rest_comp2];
            end
        end
    end
    remove_2 = [comp(2:end);remove_2]; 
    sm(:, :, remove_2) = [];
end
saveData = 0;
sm_dsum = plotMask(sm,saveData);
end