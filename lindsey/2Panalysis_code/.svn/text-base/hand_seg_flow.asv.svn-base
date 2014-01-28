%hand segmentation
clear all

%% expt params
date = '110812';
mouse = 'AC39';

userun = [8:10];
count_protocol = 1;

blanks = 1;
run = 0;

nCond = 50;
%%
P = 2;
nON = 12;
nOFF = 12;
nPlanes = 1;
begin = 7;
TFSFetc = [1:3];
pre_win = [1 6];
post_win = [7 18];

base = 'G:\users\lindsey\analysisLG\active mice';
outDir = fullfile(base, mouse,date);

%%
zoloto = '\\zoloto\bigstorlab\Lindsey\fastrig\2011\May';

for iRun = 1:length(userun)
    fn_in = fullfile(zoloto, [date '_' mouse], ['run' num2str(userun(iRun))], ['run' num2str(userun(iRun)) '_green']);
    stack = readtiff(fn_in);
    stack_dec = stackGroupProject(stack,12);
    clear stack

    %fix lag on decimated stack
    [outStack,lagP] = stackFixBidirPhase(stack_dec);
    clear outStack
    clear stack_dec

    lagP_interp = round(interp(lagP,12));

    %fix lag on original stack
    list = dir(fullfile(fn_in,'*.tif')); 
    nfiles = length(list);
    start = 1;
    for ifile = 1:nfiles;
        stack = readtiff(fn_in, ifile);
        z = size(stack,3);
        stack_shiftcorr = imBidirShift(stack, lagP_interp(:,start:start+z-1));
        fn_out = fullfile(outDir, 'shiftcorr', [date '_' mouse '_run' num2str(userun(iRun))],[date '_' mouse '_run' num2str(userun(iRun)) sprintf('_%06d.tif',ifile)]);
        writetiff(stack_shiftcorr, fn_out);
        start = start+z;
    end

    %decimate lag-shifted stack
    fn_corr = fullfile(outDir, 'shiftcorr', [date '_' mouse '_run' num2str(userun(iRun))]);
    stack_shiftcorr = readtiff(fn_corr);
    stack_shiftcorr_dec = stackGroupProject(stack_shiftcorr,12);
    fn_out = fullfile(outDir, [date '_' mouse '_run' num2str(userun(iRun)) '_stack_shiftcorr_dec.tif']);
    writetiff(stack_shiftcorr_dec, fn_out);
    clear stack_shiftcorr

    %register decimated stack
    if iRun ==1
        av = mean(stack_shiftcorr_dec(:,:,stable),3);
    end
    [out reg] = stackRegister(stack_shiftcorr_dec,av,10);
    fn_out = fullfile(outDir, [date '_' mouse '_run' num2str(userun(iRun)) '_dec_reg.tif']);
    writetiff(reg, fn_out);
    fn_out = fullfile(outDir, [date '_' mouse '_run' num2str(userun(iRun)) '_reg_out.mat']);
    save(fn_out, 'out');
    clear reg
    clear stack_shiftcorr_dec
end

%% 
eval(['PARAMS_' date '_' mouse]);
resort_seq_only

stack_sort

%%
fn_stack = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_sorted.tif']);
PrinCompLG(fn_stack);

%% figure principal components
fn_pcs = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_sorted_pca_usv.mat']);
pcs = load(fn_pcs);
nt = size(pcs.v,1);

figure;
sm = stackFilter(pcs.U,1.5);

ax=[];
for pc = 1:16;
    ax(pc)=subplot(4,4,pc);
    imagesc(sm(:,:,pc));
    
    %imstretch(sm(:,:,pc),[.5 .99],1.5);
end;
colormap gray;

%% compute independent components
PCuse = [2:75];
mu = .5;
nIC = 32;
termtol = 1e-6;
maxrounds = 400;
mixedsig = pcs.v';
mixedfilters = pcs.U;
CovEvals = diag(pcs.s).^2;
[ica_sig, ica_filters, ica_A, numiter] = CellsortICA(mixedsig, ...
    mixedfilters, CovEvals, PCuse, mu, nIC,[],termtol,maxrounds);

% f0 = mean(reg,3);
dt = 1/frGetFrameRate;
tt = [0:nt-1]/frGetFrameRate;

% plot independent filters
sel = 1:nIC;
cs = permute(ica_filters,[2,3,1]);
sm = stackFilter(cs,1.5);
figure;
ind = 1;
for ic = sel
    if ind == 17
        figure;
        ind = 1;
    end
    subplot(4,4,ind);
    imstretch(sm(:,:,ic),[.5 .99],1.5);
    ind = ind+1;
    text(.8,.1,num2str(ic),'fontsize',12,'color','w','fontweight','bold','unit','norm');
end;
fn_ica = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_sorted_ica.mat']);
save(fn_ica, 'ica_sig', 'ica_filters', 'nIC', 'PCuse');

%% consolidate ics
r_corr_filter= zeros(nIC);
sm_r = reshape(sm, [240*256 nIC]);
for ic1 = 1:nIC
    for ic2 = 1:nIC
        r_corr_filter(ic1,ic2) = triu2vec(corrcoef(sm_r(:,ic1),sm_r(:,ic2)));
    end
end
figure
imagesc(r_corr_filter)
axis image;colorbar; colormap hot

[i j] = find(r_corr_filter>.3 & r_corr_filter<.9);
ic_remove_list = [];
ica_sig_consol = ica_sig;
for a = 1:length(i)
    ica_sig_temp = sum(ica_sig([i(a) j(a)],:),1);
    if i(a)>j(a)
        ic_remove = i(a);
        ica_sig_consol(i(a),:) = ica_sig_temp;
    else
        ic_remove = j(a);
        ica_sig_consol(j(a),:) = ica_sig_temp;
    end
    ic_remove_list = [ic_remove_list ic_remove];
end
ic_remove_list = unique(ic_remove_list);
%% average ic timecourses
roi_avg = ica_sig_consol';
roi_sort_norm


resp_avg_up = resp_avg(:,1:2:nCond);
resp_avg_down = resp_avg(:,2:2:nCond);
resp_avg_both = zeros(size(resp_avg_up));
start = 0;
for iCond = 1:nCond/2
resp_avg_both(:,iCond) = mean(resp_avg(:,iCond+start:iCond+start+1),2);
start = start+1;
end

figure;
sub = 4;
for ic = 1:16;
    resp_sq = reshape(resp_avg_up(ic,:),[size(uvar,1) size(uvar,1)])';
    subplot(sub, sub ,ic)
    imagesc(resp_sq); title(num2str(ic)); colormap gray; axis image;
end

figure;
sub = 4;
for ic = 1:16;
    resp_sq = reshape(resp_avg_down(ic,:),[size(uvar,1) size(uvar,1)])';
    subplot(sub, sub ,ic)
    imagesc(resp_sq); title(num2str(ic)); colormap gray; axis image;
end

figure;
sub = 4;
for ic = 1:16;
    resp_sq = reshape(resp_avg_both(ic,:),[size(uvar,1) size(uvar,1)])';
    subplot(sub, sub ,ic)
    imagesc(resp_sq); title(num2str(ic)); colormap gray; axis image;
end

resp_sum = sum(resp_avg_both(1:16,:),1);
resp_sq= reshape(resp_sum,[size(uvar,1) size(uvar,1)])';
figure;
imagesc(resp_sq); colormap gray; axis image;title('sum');

nCells = 16;
Nshuf = 0;
SFTF_fit_LG
%% choose ics
%correlation of spatial filters
ttest_USE = find(H_ttest ==1);
for ic= 1:length(ic_remove_list)
    rem = ic_remove_list(ic);
    ind = find(ttest_USE == rem);
    ttest_USE(ind)= [];
end

figure; tcOffsetPlot(ica_sig_consol(ttest_USE,:)');

ic_USE = ttest_USE(1:end-2);

ic_USE_sum = sum(resp_avg(ic_USE,:),1);
ic_USE_sum_sq = reshape(ic_USE_sum(:,1:nCond),[6 6])';
figure; 
sub = ceil(sqrt(length(ic_USE)+1));
for ic = 1:length(ic_USE)
    subplot(sub,sub,ic)
    resp_sq = reshape(resp_avg(ic_USE(ic),1:nCond),[6 6])';
    imagesc(resp_sq); colormap gray;  colorbar; axis image; title(num2str(ic_USE(ic)));
end
subplot(sub,sub,ic+1)
imagesc(ic_USE_sum_sq); colormap gray; colorbar; axis image; title('sum');

%% find average frames for ics
fn_stack = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_sorted.tif']);
stack = readtiff(fn_stack);
[x y z] = size(stack);
ic_avg = zeros(x,y,nIC);
for ic = 1:nIC
    ind = find(ica_sig(ic,:)>0.05);
    if size(ind)>0;
       ic_avg(:,:,ic) = mean(stack(:,:,ind),3); 
    end
end

ic_stack = zeros(x,y,2*nIC);
start = 1;
for ic = 1:nIC
    ic_stack(:,:,start) = sm(:,:,ic)./max(max(sm(:,:,ic),[],2),[],1);
    ic_stack(:,:,start+1) = ic_avg(:,:,ic)./max(max(ic_avg(:,:,ic),[],2),[],1);
    start = start+2;
end
fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_ic_stack.tif']);
writetiff(ic_stack, fn_out);


%% make rois and get TCs
pts =4;
manual_rois

fn = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_sorted.tif']);
roi_TCs

%% make stack without signal correlations
fn = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_sorted.tif']);
stack = readtiff(fn);
Stack_uncorr   
fn = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_uncorr.tif']);
roi_TCs

%% do correllations of rois with pca
fn = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_uncorr.tif']);
PrinCompLG(fn);
fn_pcs = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_uncorr_pca_usv.mat']);
pcs = load(fn_pcs);


nRoi = size(roi_avg,2);
r = zeros(300, nRoi);
im_corr = zeros(240, 256, nRoi);
for iRoi = 1:nRoi;
    for ix = 1:300
       temp = corrcoef(pcs.v(:,ix),roi_dF(:,iRoi));
       r(ix,iRoi)=temp(1,2);
    end;
    im = zeros(240*256, 1);
    u = reshape(pcs.U, [240*256 300]);
    for ix = 1:300
       im = im + r(ix,iRoi)*u(:,ix);
    end
    im_corr_pca(:,:,iRoi) = reshape(im,[240 256]);
end

ind = 1;
figure;
for iRoi = 1:nRoi;
    if ind == 17
        figure;
        ind = 1;
    end
    subplot(4,4,ind)
    imstretch(im_corr_pca(:,:,iRoi));
    axis square;colormap hot;
    ind = ind+1;
end

%% do correlations from data
fn = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_uncorr.tif']);
stack = readtiff(fn);
sz = size(stack);
stack2 = reshape(stack,[sz(1)*sz(2) sz(3)]);
clear stack
nRoi = size(roi_dF,2);
im_corr = zeros(nRoi, sz(1)*sz(2));
for iRoi = 1:nRoi
   im = zeros(1,240*256);
   for ix = 2:240*256-1
        if mod(ix,1000)==0
        fprintf('%i\n',ix);
        end;
    temp = corrcoef(double(mean(stack2(ix-1:ix+1,:),1)),roi_dF(:,iRoi))';
    im(ix)=temp(1,2);
    end;
    im_corr(iRoi,:) = im;
end


im_corr2 = reshape(im_corr,[nRoi sz(1) sz(2)]);
ind = 1;
figure;
for iRoi = 1:nRoi;
    if ind == 17
        figure;
        ind = 1;
    end
    subplot(4,4,ind)
    imstretch(squeeze(im_corr2(iRoi,:,:)));
    axis square;colormap gray;
    ind = ind+1;
end

fn = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_uncorr_r.mat']);
save(fn,'im_corr2');
%% digitize data
fn = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_uncorr.tif']);
stack = readtiff(fn);

sz = size(stack);

roi_dF_dec = tcDecimate(roi_dF, 5);
event_plot = zeros(nRoi, size(roi_dF_dec,1));
events_all = zeros(nRoi, size(roi_dF_dec,1));
for iRoi = 1:nRoi;
    var = std(roi_dF_dec(:,iRoi));
    ind = find(roi_dF_dec(:,iRoi)>(2.5*var));
    events = zeros(1,size(roi_dF_dec,1));
    events(1,ind) = iRoi;
    event_plot(iRoi,:) = events;
    events_all(iRoi,:) = events;
    ind2 = find(event_plot(iRoi,:) == 0);
    event_plot(iRoi,ind2) = NaN;
end

for ic = 1:max(roi_ic)
    ind = find(roi_ic==ic);
    if length(ind) > 0
        figure;
        for iRoi = 1:length(ind);
        roiA = ind(iRoi);
        if length(ind) == 1
            subplot(1,1,iRoi)
        else
            subplot(ceil(length(ind)/2),2,iRoi)
        end
        plot(squeeze(roi_dF_dec(:,roiA)))
        hold on
        plot(event_plot(roiA,:)-roiA,'*r')
        title([num2str(ic) ' roi' num2str(roiA)]);
        end
    end
end


figure;
for iRoi = 1:nRoi
    plot(event_plot(iRoi,:), '*r')
    hold on
end

figure;
for iRoi1 = 1:nRoi
    for iRoi2 = 1:nRoi
        r_corr_spikes(iRoi1, iRoi2) = triu2vec(corrcoef(events_all(iRoi1,:),events_all(iRoi2,:)));
    end
end
imagesq(r_corr_spikes,[0 1]); axis image; colormap hot; colorbar;

figure;
start =1;
for ic = 1:max(roi_ic)
    ind = find(roi_ic==ic);
    if length(ind)>1
        r_corr_thresh = r_corr_spikes(start:start+length(ind)-1,start:start+length(ind)-1);
        thresh = find(r_corr_thresh<0.65);
        r_corr_thresh(thresh) = 0;
        start = start+length(ind);
        subplot(ceil(sqrt(max(roi_ic))),ceil(sqrt(max(roi_ic))),ic)
        imagesq(r_corr_thresh,[0 1]);
        colormap hot;
        colorbar;
        title([num2str(ind(1)) '-' num2str(ind(end))]);
    elseif length(ind) == 1
        r_corr_thresh = r_corr_spikes(start:start+length(ind)-1,start:start+length(ind)-1);
        subplot(ceil(sqrt(max(roi_ic))),ceil(sqrt(max(roi_ic))),ic)
        imagesq(r_corr_thresh, [0 1]);
        colormap hot;
        colorbar;
        title(num2str(ind(1)));
        start = start+1;
    end
end    

%% consolidate rois
for ic = 1:max(roi_ic)
end
%% sort and average responses
fn_rois = fullfile(outDir, [date '_' mouse '_run' num2str(userun) '_roiTCs.mat']);
load(fn_rois);

roi_sort_norm

figure;
start = 1;
for ic = 1:nIC
    if start == 17
        figure;
        start=1;
    end
    resp_avg_sq = reshape(resp_avg(ic,1:nCond),[6 6])';
    subplot(4,4, start)
    imagesc(resp_avg_sq); axis image, colormap gray;
    if H_ttest(ic,1) == 1
        title([num2str(ic) '**']);
    else
        title([num2str(ic)]);
    end
    start = 1+start;
end
    

ind = find(H_ttest==1);    
figure
tcoffsetplot(ica_sig(ind,:)')

col = 'k';
Make_TC_subplots

fn = fullfile(base,mouse,date,'analysis',[date '_' mouse '_run' num2str(userun) 'roi_peak_resp.mat']);
save(fn, 'resp_avg', 'resp_std');

%plot tuning curves
var_plot
edit var_SFTF_plot

%% Use ROIs to find segs with local correlations

%find local correlations from raw dataset
stack_pn = fullfile(outDir,[date '_' mouse '_run' num2str(userun) '_sorted.tif']);
stack = readtiff(stack_pn);
Local_corr

fn_out = fullfile(outDir,[date '_' mouse '_run' num2str(userun) 'manualroi_localcorr_thresh4pt5.mat']);
save(fn_out, 'r');

%remove average correlations and find local correlations from residuals
stack_pn = fullfile(outDir,[date '_' mouse '_run' num2str(userun) '_sorted.tif']);
stack = readtiff(stack_pn);
Stack_uncorr

stack_pn = fullfile(outDir,[date '_' mouse '_run' num2str(userun) '_stack_uncorr.tif']);
stack = readtiff(stack_pn);
Local_corr
fn_out = fullfile(outDir,[date '_' mouse '_run' num2str(userun) 'manualroi_localcorr_uncorr.mat']);
save(fn_out, 'r');

%threshold and make masks from local correlations
siz = size(stack);
roi_mask = zeros(siz(1)*siz(2), size(r,3));
stack_ghost = zeros(siz(1),siz(2));
for iRoi = 1:size(r,3);
    thresh = 0.7*max(max(r(:,:,iRoi),[],1),[],2);
    mask = find(r(:,:,iRoi)>thresh);
    roi_mask(mask, iRoi) = 1;
end

fn_out = fullfile(outDir,[date '_' mouse '_run' num2str(userun) 'localcorr_mask.mat']);
save(fn_out, 'mask');

roi_masks = reshape(roi_mask,[siz(1) siz(2) size(r,3)]);
subs = ceil(sqrt(size(r,3)));
figure;
for iRoi = 1:size(r,3);
    subplot(subs, subs, iRoi);
    imagesq(roi_masks(:,:,iRoi));
    colormap gray
end

%get timecourses from masks
stack_pn = fullfile(outDir,[date '_' mouse '_run' num2str(userun) '_sorted.tif']);
stack = readtiff(stack_pn);

mask_avg = zeros(siz(3),size(r,3));
stack_r = reshape(stack,[siz(1)*siz(2) siz(3)]);
for iRoi = 1:size(r,3);
    mask = find(roi_mask(:,iRoi)==1);
    mask_avg(:,iRoi) = squeeze(mean(stack_r(mask,:),1));
end

fn_out = fullfile(outDir,[date '_' mouse '_run' num2str(userun) 'localcorr_mask_TCs.mat']);
save(fn_out, 'mask_avg');
mask_dF = bsxfun(@minus, mask_avg, mean(mask_avg,1));
figure;
tcOffsetPlot(mask_dF);
%sort and average responses
roi_avg = mask_avg;
% roi_sort_norm
fn = fullfile(base,mouse,date,'analysis',[date '_' mouse '_run' num2str(userun) 'mask_peak_resp.mat']);
save(fn, 'resp_avg', 'resp_std');

%plot tuning curves
var_plot



