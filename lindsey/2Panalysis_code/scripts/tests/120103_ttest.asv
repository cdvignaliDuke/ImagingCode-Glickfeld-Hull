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
stack_sorted = readtiff(fn_stack);
siz = size(stack_sorted);
stack_sorted_np = zeros(size(stack_sorted));
avg = mean(mean(mean(stack_sorted,3),2),1);
for iframe = 1:size(stack_sorted,3)
    np_av = mean(mean(stack_sorted(:,:,iframe),2),1);
    stack_sorted_np(:,:,iframe) = stack_sorted(:,:,iframe)-np_av+avg;
end

fn_stack = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_sorted_npsub.tif']);
writetiff(stack_sorted_np, fn_stack);

fn_mask = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_axon_mask.mat']);
load(fn_mask);
fn = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_reps.mat']);
load(fn);

nReps = sum(stim_reps(1,:));

roi_stim = zeros(siz(1), siz(2), nOFF+nON,nReps);
start = 1;
rep = 1;
for iCond = 1:nCond+1; 
    nRep = stim_reps(iCond);
    for iRep = 1:nRep;
        roi_stim(:,:,:, rep) = stack_sorted_np(:,:,start:start-1+nON+nOFF);
        start = start+nON+nOFF;
        rep = rep+1;
    end
end
resp_off = squeeze(mean(roi_stim(:,:,pre_win(1):pre_win(2),:),3));
resp_on = squeeze(mean(roi_stim(:,:,post_win(1):post_win(2),:),3));

fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_resp.mat']);
save(fn_out, 'resp_on', 'resp_off');

alphaB = .05./(60000);
Info_ttest_mat = zeros(siz(1),siz(2),nCond);

b= 5;

f1 = fspecial('average');
siz = size(resp_on);
resp_on_long = reshape(resp_on, [siz(1) siz(2)*siz(3)]);
resp_off_long = reshape(resp_off, [siz(1) siz(2)*siz(3)]);

resp_on_sm = reshape(filter2(f1,resp_on_long),[siz(1) siz(2) siz(3)]);
resp_off_sm = reshape(filter2(f1,resp_off_long),[siz(1) siz(2) siz(3)]);

for iy = b+1:240-b
    fprintf([num2str(iy) ' '])
    for ix = b+1:256-b
        start = 1;
        p_ttestB = zeros(1,1,nCond);
        for iCond = 1:nCond
            nRep = stim_reps(1,iCond);
            [h_ttestB1,p_ttestB1] = ttest(sub_off(iy,ix,start:start-1+nRep),sub_on(iy,ix,start:start-1+nRep),alphaB,'left');
            p_ttestB(1,1,iCond) = p_ttestB1;
            start = start+nRep;
        end
    Info_ttest_mat(iy,ix,:) = p_ttestB;
    end
end
ttest_mask = 1-squeeze(min(Info_ttest_mat,[],3));
ttest_mask_npix = min(Info_ttest_mat,[],3) < alphaB;

siz = size(Info_ttest_mat);
Info_ttest_mat_long = reshape(Info_ttest_mat, [siz(1) siz(2)*siz(3)]);
ttest_mask = reshape(filter2(f1,Info_ttest_mat_long), [siz(1) siz(2) siz(3)]);
figure; imagesq(ttest_mask);

ttest_mask_nCond(1:5,1:end) = 0;
ttest_mask_nCond(1:end, 1:5) = 0;
ttest_mask_nCond(1:end, 251:end) = 0;
ttest_mask_nCond(235:end,1:end) = 0;
figure
for iCond = 1:nCond
    subplot(5,5,iCond)
    imagesq(Info_ttest_mat(:,:,iCond));
end
fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_ttest_mask.mat']);
save(fn_out, 'ttest_mask', 'Info_ttest_mat');

siz = size(ttest_mask);
figure; imagesq(corr_mask);
colorbar
r_sorted = sort(reshape(r, [siz(1)*siz(2), 1]));
ttest_percent = sum(reshape(ttest_mask, [siz(1)*siz(2) 1]),1)/(siz(1)*siz(2));

percent = [0.25 0.5 1 2 5 10 15 20 25 40];     
corr_mask_long = zeros(size(r_sorted));
for iper = 1:10;
    rank = (siz(1)*siz(2))-ceil((5.72/100*(siz(1)*siz(2))-1));
    thresh = r_sorted(rank);
    ind = r>thresh;
    corr_mask_long(ind) = 1;
end
    corr_mask = reshape(corr_mask_long,[siz(1) siz(2)]);
    
    figure;
    for iper = 1:10
        subplot(4,3,iper)
        imagesq(corr_mask(:,:,iper));
    end
    

overlay = 2*corr_mask + ttest_mask;
figure; 
subplot(3,1,1)
imagesq(ttest_mask)
colorbar;
title('ttest')
subplot(3,1,2)
imagesq(corr_mask)
title('correlation')
colorbar;
subplot(3,1,3)
imagesq(overlay);
title('merge')
colorbar;

figure; imagesq(r); colorbar

fn_dF = fullfile(outDir, 'analysis',[date '_' mouse '_run' num2str(userun) '_stack_dF_all.tif']);
stack_dF = readtiff(fn_dF);
stack_dF_long = reshape(stack_dF,[siz(1)*siz(2) nCond+1]);
corr_mask_ind = find(corr_mask ==1);
ttest_mask_ind = find(ttest_mask ==1);
resp_avg_corr = mean(stack_dF_long(corr_mask_ind,:),1)';
resp_avg_ttest = mean(stack_dF_long(ttest_mask_ind,:),1)';

figure; subplot(3,1,1);  imagesq(ttest_mask_npix); title('pixel correction');subplot(3,1,2);  imagesq(ttest_mask_nCond);title('condition correction'); subplot(3,1,3); imagesq(ttest_mask_smooth); title('smoothed condition correction');
figure; subplot(3,1,1); imagesq(ttest_mask_smooth); title('smoothed condition correction'); subplot(3,1,2);  imagesq(corr_mask); title('matched number of pixels from correlation'); subplot(3,1,3); imagesq(overlay); title('overlay');
