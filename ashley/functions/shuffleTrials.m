function [tc1_shuff, tc2_shuff] = shuffleTrials(tc1, tc2);
all_tr = cat(3,tc1,tc2);
nt = size(all_tr,3);

ind1_shuff = randsample(nt,floor(nt/2));
ind2_shuff = setdiff(1:nt,ind1_shuff);

tc1_shuff = mean(all_tr(:,:,ind1_shuff),3);
tc2_shuff = mean(all_tr(:,:,ind2_shuff),3);
end