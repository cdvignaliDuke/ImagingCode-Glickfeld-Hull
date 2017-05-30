function [ data_corr, sm, sm_flat] = finalMask2( cell_sig, ica_segments, threshold, out_dir)
% nmask = max(max(mask_final));
% %nmask = size(sel,2);
% sm = zeros(npw,nph,nmask);
%
% for ii=1:nmask
%     sm_placeholder = zeros(npw,nph);
%     sm_placeholder(mask_final==ii) = 1;
%     sm(:,:,ii)= sm_placeholder;
% end
sm = permute(ica_segments,[2,3,1]);
[npw, nph, nmask] = size(sm);

% mask_overlap = zeros(1,npw*nph);
% mask_all = zeros(1,npw*nph);
% sm_temp = reshape(sm,npw*nph,nmask);
ic = 1;
while ic < nmask
    for kc = ic+1:nmask
        sm_temp = sm(:,:,kc);
        sm_temp(sm(:,:,ic)&sm(:,:,kc)) = 0;
        sm(:,:,kc) = sm_temp;
        if length(find(sm(:,:,kc))) < 200
            sm(:,:,kc) = [];
            nmask = nmask - 1;
        end
    end
    
    ic = ic + 1;
end

data_corr = corrcoef(cell_sig);
figure; fig = imagesc(data_corr);
% saveas(fig, [out_dir, 'data_corr.fig']);
% print([out_dir, 'data_corr.eps'],'-depsc')


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
        temp = sm(:,:,comp(k));
        sm_temp(find(temp),comp(1)) = temp(find(temp));
        sm(:, :, comp(1)) = reshape(sm_temp(:,comp(1)), npw ,nph);
    end
    if length(rest_comp) == 2
        temp = sm(:,:,rest_comp(2));
        sm_temp(find(temp),rest_comp(1)) = temp(find(temp));
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
                    temp = sm(:,:,rest_comp2(kk));
                    sm_temp(find(temp),rest_comp1) = temp(find(temp));
                    sm(:, :, rest_comp1) = reshape(sm_temp(:,rest_comp1), npw ,nph);
                end
                remove_2 = [remove_2; rest_comp2];
            end
        end
    end
    remove_2 = [comp(2:end);remove_2];
    sm(:, :, remove_2) = [];
end

sm_flat = sum(sm,3);

figure;
fig = image(sm_flat);
axis image;
% saveas(fig, [out_dir, 'mask_2.fig']);
% print([out_dir, 'mask_2.eps'],'-depsc')
end