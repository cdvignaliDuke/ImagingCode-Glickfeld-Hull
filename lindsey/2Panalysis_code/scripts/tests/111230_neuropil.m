mouse = 'Y13';
date = '110509';
userun = [1:2];

count_protocol = 1;
blanks = 1;
run = 0;

nCond =25;
dir = 1;
P = 2;
nON = 12;
nOFF = 12;
nPlanes = 1;
begin = 1;
TFSFetc = [1:2];
pre_win = [7 12];
post_win = [13 24];


base = 'G:\users\lindsey\analysisLG\active mice';
outDir = fullfile(base, mouse,date);

fn_stack = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_sorted.tif']);
stack = readtiff(fn_stack);
stack = single(stack);
siz = size(stack);
ntrials = siz(3)/(nOFF+nON);
base_std = zeros(siz(1),siz(2),ntrials);
for itrial = 1:ntrials
    base_std(:,:,itrial)= std(stack(:,:,7+(itrial-1)*(nON+nOFF):12+(itrial-1)*(nON+nOFF)),[],3);
end
fn_mask = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_axon_mask.mat']);
load(fn_mask)

noise = mean(base_std,3);
noise_long = reshape(noise, [siz(1)*siz(2),1]);
fn_dF = fullfile(outDir, 'analysis',[date '_' mouse '_run' num2str(userun) '_stack_dF_all.tif']);
stack_dF = readtiff(fn_dF);
stack_dF_long = reshape(stack_dF, [siz(1)*siz(2) nCond+1]);
percent = [.25 .5 1 2 5 10 20 50 75 100];
r_long = reshape(r,[siz(1)*siz(2),1]);
r_sorted= sort(r_long);
mask = zeros(siz(1)*siz(2), 10);
resp_avg = zeros(26,10);
base_std = zeros(1,10);
snr = zeros(1,10);

for iper = 1:10
    rank = (siz(1)*siz(2))-ceil(((percent(iper)/100)*(siz(1)*siz(2))-1));
    thresh = r_sorted(rank);
    ind = find(r>thresh);
    mask(ind,iper) =1;
    resp_avg(:,iper) =mean(stack_dF_long(ind,:),1)*100;
    base_std(:,iper) = mean(noise_long(ind,:),1);
    snr(:,iper) = max(resp_avg(:,iper),[],1)/base_std(:,iper);
end

figure
for iper = 1:10
    subplot(4,3,iper)
    mask_sq = reshape(mask(:,iper),[siz(1),siz(2)]);
    imagesq(mask_sq);
end


stack_sorted_np = zeros(size(stack_sorted));
avg = mean(mean(mean(stack_sorted,3),2),1);
for iframe = 1:size(stack_sorted,3)
    np_av = mean(mean(stack_sorted(:,:,iframe),2),1);
    stack_sorted_np(:,:,iframe) = stack_sorted(:,:,iframe)-np_av+avg;
end

stack_dF_np = zeros(siz(1), siz(2), nCond+1);
start = 0;
for iCond = 1:Cond;
    nRep = length(Big_Seqposition(iCond).ind);
    rep_dF = zeros(siz(1), siz(2),nRep);
    for iRep = 1:nRep
        rep_base = mean(stack_sorted_np(:,:,start+pre_win(1):start+pre_win(2)),3);
        rep_resp = mean(stack_sorted_np(:,:,start+post_win(1):start+post_win(2)),3);
        rep_dF(:,:,iRep) = (rep_resp-rep_base)./rep_base;
        start = ((nOFF+nON)/nPlanes)+start;
    end
    stack_dF_np(:,:,iCond) = mean(rep_dF,3);
end
fn_out = fullfile(outDir, 'analysis',[date '_' mouse '_run' num2str(userun) '_stack_dF_np.tif']);
writetiff(stack_dF_np, fn_out);

siz = size(r);
axon_mask_pt2 = zeros(size(r));
axon_mask_pt2(find(r>0.2)) = 1;
nhood = ones(2,2);
axon_mask_dil = imdilate(axon_mask_pt2,nhood);
axon_mask_dil(1:5,:) = 1;
axon_mask_dil(:,1:5) = 1;
axon_mask_dil(siz(1)-5:siz(1),:) = 1;
axon_mask_dil(:,siz(2)-5:siz(2)) = 1;

npMA = find(axon_mask_dil ==0);
stack_sorted_npMA = zeros(size(stack_sorted));
for iframe = 1:size(stack_sorted,3)
    stack_reshape = squeeze(reshape(stack_sorted(:,:,iframe),[siz(1)*siz(2) 1]));
    npMA_av = mean(stack_reshape(npMA),1);
    stack_sorted_npMA(:,:,iframe) = stack_sorted(:,:,iframe)-npMA_av+avg;
end

stack_dF_npMA = zeros(siz(1), siz(2), nCond+1);
start = 0;
for iCond = 1:Cond;
    nRep = length(Big_Seqposition(iCond).ind);
    rep_dF = zeros(siz(1), siz(2),nRep);
    for iRep = 1:nRep
        rep_base = mean(stack_sorted_npMA(:,:,start+pre_win(1):start+pre_win(2)),3);
        rep_resp = mean(stack_sorted_npMA(:,:,start+post_win(1):start+post_win(2)),3);
        rep_dF(:,:,iRep) = (rep_resp-rep_base)./rep_base;
        start = ((nOFF+nON)/nPlanes)+start;
    end
    stack_dF_npMA(:,:,iCond) = mean(rep_dF,3);
end


fn_out = fullfile(outDir, 'analysis',[date '_' mouse '_run' num2str(userun) '_stack_dF_npMA.tif']);
writetiff(stack_dF_npMA, fn_out);

fn_dF = fullfile(outDir, 'analysis',[date '_' mouse '_run' num2str(userun) '_stack_dF_all.tif']);
fn_dF_np = fullfile(outDir, 'analysis',[date '_' mouse '_run' num2str(userun) '_stack_dF_np.tif']);
fn_dF_npMA = fullfile(outDir, 'analysis',[date '_' mouse '_run' num2str(userun) '_stack_dF_npMA.tif']);

corr = [0.4:.1:.8];
resp_avg = zeros(26,5);
resp_avg_np = zeros(26,5);
resp_avg_npMA = zeros(26,5);
stack_dF = readtiff(fn_dF);
stack_dF_np = readtiff(fn_dF_np);
stack_dF_npMA = readtiff(fn_dF_npMA);
siz = size(stack_dF);
stack_dF_long = reshape(stack_dF,[siz(1)*siz(2) nCond+1]);
stack_dF_np_long = reshape(stack_dF_np,[siz(1)*siz(2) nCond+1]);
stack_dF_npMA_long = reshape(stack_dF_npMA,[siz(1)*siz(2) nCond+1]);
beg = 1;
figure;
    for icorr =1:5    
        axon_mask_ind = find(r>.2);
        subplot(5,4,beg)
        imagesq(axon_mask(:,:,icorr));
        resp_avg_temp  = zeros(26,1);
        resp_avg_temp_np  = zeros(26,1);
        resp_avg_temp_npMA  = zeros(26,1);
        start = 1;
        if dir == 2
            for iCond = 1:25
                resp_avg_temp(iCond,:) = mean(mean(stack_dF_long(axon_mask_ind,start:start+1),2),1);
                resp_avg_temp_np(iCond,:) = mean(mean(stack_dF_np_long(axon_mask_ind,start:start+1),2),1);
                resp_avg_temp_npMA(iCond,:) = mean(mean(stack_dF_npMA_long(axon_mask_ind,start:start+1),2),1);
                start = start+2;
            end
            resp_avg_temp(26,:) = mean(stack_dF_long(axon_mask_ind,end),1);
            resp_avg_temp_np(26,:) = mean(stack_dF_np_long(axon_mask_ind,end),1);
            resp_avg_temp_npMA(26,:) = mean(stack_dF_npMA_long(axon_mask_ind,end),1);
        elseif dir ==1
            resp_avg_temp = mean(stack_dF_long(axon_mask_ind,:),1)';
            resp_avg_temp_np = mean(stack_dF_np_long(axon_mask_ind,:),1)';
            resp_avg_temp_npMA = mean(stack_dF_npMA_long(axon_mask_ind,:),1)';
        end
        resp_avg(:,icorr) = resp_avg_temp;
        resp_avg_np(:,icorr) = resp_avg_temp_np;
        resp_avg_npMA(:,icorr) = resp_avg_temp_npMA;
        resp_avg_sq = reshape(resp_avg_temp(1:25,:), [5 5])';
        resp_avg_sq_np = reshape(resp_avg_temp_np(1:25,:), [5 5])';
        resp_avg_sq_npMA = reshape(resp_avg_temp_npMA(1:25,:), [5 5])';
        subplot(5,4,beg+1)
        imagesq(resp_avg_sq);colormap(gray); colorbar
        subplot(5,4,beg+2)
        imagesq(resp_avg_sq_np);colormap(gray); colorbar
        subplot(5,4,beg+3)
        imagesq(resp_avg_sq_npMA);colormap(gray); colorbar
        beg = beg+4;
    end
    subplot(5,4,2); title('Raw');  subplot(5,4,3); title('NP sub');  subplot(5,4,2); title('NP MA'); 
    
    resp_avg_np = mean(stack_dF_long(npMA,:),1)';
    resp_avg_np_sq = reshape(resp_avg_np(1:25,:), [5 5])';
    npMA_mask_long = zeros(siz(1)*siz(2), 1);
    npMA_mask_long(npMA) = 1;
    mask_long = zeros(siz(1)*siz(2), 1);
    mask_long(axon_mask_ind) = 1;
    npMA_mask = reshape(npMA_mask_long, [siz(1) siz(2)]);
    resp_avg_all = mean(stack_dF_long,1)';
    resp_avg_all_sq = reshape(resp_avg_all(1:25,:), [5 5])';
    mask = reshape(mask_long, [siz(1) siz(2)]);
    figure; subplot(4,2,1); imagesq(ones(size(r)));  title('whole field'); subplot(4,2,2);  imagesq(resp_avg_all_sq); colorbar;  
    subplot(4,2,3);  imagesq(mask);  title('axon mask r>0.2');  subplot(4,2,4);  imagesq(resp_avg_sq); colorbar;
    subplot(4,2,5); imagesq(axon_mask(:,:,1));  title('axon mask r>0.4'); subplot(4,2,6);  imagesq(resp_avg_sq); colorbar;
    subplot(4,2,7);  imagesq(npMA_mask);  title('NP mask'); subplot(4,2,8);  imagesq(resp_avg_np_sq);colorbar; 
    colormap(gray);
    
    
r_np = zeros(siz(1), siz(2));
r_npMA = zeros(siz(1), siz(2));
b= 5;
for iy = b+1:240-b
    fprintf('.');
    for ix = b+1:256-b
            sub_np = stack_dF_np(iy-1:iy+1,ix-1:ix+1,:);
            sub_npMA = stack_dF_npMA(iy-1:iy+1,ix-1:ix+1,:);
            r_np(iy,ix)=mean(triu2vec(corrcoef(reshape(sub_np,[3*3,size(stack_dF_np,3)])'),1));
            r_npMA(iy,ix)=mean(triu2vec((corrcoef(reshape(sub_npMA,[3*3,size(stack_dF_npMA,3)])')),1));
    end;
end;

axon_mask_np = zeros(siz(1),siz(2), 5);
axon_mask_npMA = zeros(siz(1),siz(2), 5);
beg = 1;
figure;
    for icorr =1:5    
        axon_mask_ind_np = find(r_np>corr(icorr));
        axon_mask_ind_npMA = find(r_npMA>corr(icorr));
        axon_mask_temp_np = zeros(siz(1)*siz(2),1);
        axon_mask_temp_np(axon_mask_ind_np,1) = 1;
        axon_mask_np(:,:,icorr) = reshape(axon_mask_temp_np, [siz(1) siz(2)]);
        axon_mask_temp_npMA = zeros(siz(1)*siz(2),1);
        axon_mask_temp_npMA(axon_mask_ind_npMA,1) = 1;
        axon_mask_npMA(:,:,icorr) = reshape(axon_mask_temp_npMA, [siz(1) siz(2)]);
        
        resp_avg_temp_np  = zeros(26,1);
        resp_avg_temp_npMA  = zeros(26,1);
        start = 1;
        if dir == 2
            for iCond = 1:25
                resp_avg_temp_np(iCond,:) = mean(mean(stack_dF_np_long(axon_mask_ind_np,start:start+1),2),1);
                resp_avg_temp_npMA(iCond,:) = mean(mean(stack_dF_npMA_long(axon_mask_ind_npMA,start:start+1),2),1);
                start = start+2;
            end
            resp_avg_temp_np(26,:) = mean(stack_dF_np_long(axon_mask_ind_np,end),1);
            resp_avg_temp_npMA(26,:) = mean(stack_dF_npMA_long(axon_mask_ind_npMA,end),1);
        elseif dir ==1
            resp_avg_temp_np = mean(stack_dF_np_long(axon_mask_ind_np,:),1)';
            resp_avg_temp_npMA = mean(stack_dF_npMA_long(axon_mask_ind_npMA,:),1)';
        end
        resp_avg_np(:,icorr) = resp_avg_temp_np;
        resp_avg_npMA(:,icorr) = resp_avg_temp_npMA;
        resp_avg_sq_np = reshape(resp_avg_temp_np(1:25,:), [5 5])';
        resp_avg_sq_npMA = reshape(resp_avg_temp_npMA(1:25,:), [5 5])';
        subplot(5,4,beg)
        imagesq(axon_mask_np(:,:,icorr));
        subplot(5,4,beg+1)
        imagesq(resp_avg_sq_np);colormap(gray); colorbar
        subplot(5,4,beg+2)
        imagesq(axon_mask_npMA(:,:,icorr));
        subplot(5,4,beg+3)
        imagesq(resp_avg_sq_npMA);colormap(gray); colorbar
        beg = beg+4;
    end

     subplot(5,4,2); title('NP sub');  subplot(5,4,4); title('NP MA'); 

axon_percent = zeros(1,5);
for icorr = 1:5
    axon_mask_ind = find(r>corr(icorr));
    axon_percent(:,icorr) = size(axon_mask_ind,1)/(siz(1)*siz(2));
end
    
  subplot(5,4,1); title('3.8%');
  subplot(5,4,5); title('1.1%'); 
  subplot(5,4,9); title('0.4%'); 
  subplot(5,4,13); title('0.1%'); 
  subplot(5,4,17); title('0.01%');
  
  subplot(5,4,2); title('Raw');
  subplot(5,4,3); title('NP sub'); 
  subplot(5,4,4); title('NP MA'); 
 
