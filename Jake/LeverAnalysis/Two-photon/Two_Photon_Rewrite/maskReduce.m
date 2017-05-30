
[bd_cell] = bwboundaries(mask_raw,'noholes');
mask_reduce = mask_raw;
ind_size = size(mask_raw);
mask_bd = zeros(ind_size(1), ind_size(2));
for ic = 1:size(bd_cell,1)
    bd_indx = sub2ind(ind_size, bd_cell{ic}(:,1), bd_cell{ic}(:,2));
    mask_reduce(bd_indx) = 0;
    mask_bd(bd_indx) = 20;
end

img_reg2 = img_reg(:,:,1:10:end);
img_new = zeros(ind_size(1),ind_size(2),size(img_reg2,3));

for nf = 1:size(img_reg2,3)
    temp = imfuse(img_reg2(:,:,nf),mask_bd);
    img_new(:,:,nf) = sum(temp,3);
end

writetiff(img_new, 'mask_movie');

mask_final = reshape(mask_raw, 1, ind_size(1)*ind_size(2));

nmask = max(max(mask_final));

sm = zeros(ind_size(1), ind_size(2), nmask);



for ii=1:nmask
    sm_placeholder = zeros(ind_size(1),ind_size(2));
    sm_placeholder(mask_final==ii) = 1;
    sm(:,:,ii)= sm_placeholder;
end

tc_new = getTC(img_reg, sm, nmask);

nCells = size(tc_avg,2);
data_tc = tc_avg;
npTC = zeros(size(data_tc));

buf = 0;
np = 2;
neuropil = squeeze(imCellNeuropil(mask_final,buf,np));
np_mask = zeros(sz(1),sz(2),nCells);
for i = 1:nCells
%     ind_both = and(neuropil(:,i),mask_overlap');
%     neuropil(ind_both,i) = 0;
    np_mask(:,:,i) = reshape(neuropil(:,i),[sz(1) sz(2)]);
    npTC(:,i) = stackGetTimeCourses(img_reg,np_mask(:,:,i));
end

%get weights by maximizing skew
ii= 0.01:0.01:1;
x = zeros(length(ii), nCells);
tc_avg = tsmovavg(data_tc,'s',1,1);
np_avg = tsmovavg(npTC,'s',1,1);
for i = 1:100
    x(i,:) = skewness(tc_avg-tcRemoveDC(np_avg.*ii(i)));
end
[max_skew ind] =  max(x,[],1);
np_w = 0.01*ind;
npSubTC = data_tc-bsxfun(@times,tcRemoveDC(npTC),np_w);
save([dest '_npSubTCs.mat'],'npSubTC',  'neuropil', 'np_w');    