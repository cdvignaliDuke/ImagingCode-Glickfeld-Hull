base = 'G:\users\lindsey\analysisLG\active mice';
reg = 'G:\users\lindsey\regLG';
date = '110517';
mouse = 'Y13';
expt = 'run6';

% don't save as 32 bit tiff
%stack = readtiff(fullfile(base, mouse,date,[date '_' mouse '_' expt '_dec_reg.tif']));
stack = readtiff(fullfile(reg, [date '_' mouse], expt, [expt '_green']));
siz = size(stack);
stack_d6 = double(stackGroupProject(stack,6));
clear stack
a = zeros(240,256);
b = 5; % border

%compute local correlation
fprintf('Computation local correlation:');
for iy = b+1:240-b
    fprintf('.');
    for ix = b+1:256-b
            sub = stack_d(iy-1:iy+1,ix-1:ix+1,:);
            a(iy,ix)=mean(triu2vec((corrcoef(reshape(sub,[3*3,siz(3)])')),1));
    end;
end;
figure;
imagesq(a);

fn_out = fullfile(base, mouse,date,[date '_' mouse '_' expt '_localcorr.m']);
save(fn_out,'a');

%find discrete regions of interest
roi_image = zeros(size(sm_2));
nrois = 100;
rois = zeros(nrois,2);
a_copy = sm_2;
for roi = 1:nrois;
  [x y] = find(a_copy == max(max(a_copy(b+1:240-b,b+1:256-b))));
  rois(roi,:) = [x y];
  a_copy(x-5:x+5,y-5:y+5) = 0;
  roi_image(x,y) = 1; 
end
figure;
imagesq(roi_image);

%compute correlation for each region of interest
r = zeros(240,256,nrois);
for iroi = 1:10;
    roi = stack_d(rois(iroi,1)-3:rois(iroi,1)+3,rois(iroi,2)-3:rois(iroi,2)+3,:);
    for iy = b+1:240-b
        fprintf('.');
        for ix = b+1:256-b
            sub = stack_d(iy-3:iy+3,ix-3:ix+3,:);
            r(iy,ix,iroi)= triu2vec(corrcoef(roi,sub));
        end;
    end;
end;

fn_out = fullfile(base, mouse,date,[date '_' mouse '_' expt '_roicorr.m']); 
save(fn_out,'r');

figure;
id = 1;
for iroi = 1:10;
    subplot(3,4,id);
    imagesq(r(:,:,iroi));
    id = id+1;
end;

%remove noise
n = 0;
for iroi = 1:nrois;
    r_smooth(:,:,iroi) = reshape(smooth(r(b+1:240-b,b+1:256-b,iroi)),230,246);
    r_base(:,:,iroi) = r_smooth(:,:,iroi)-min(min(r_smooth(:,:,iroi)));
    r_norm(:,:,iroi) = r_base(:,:,iroi)./max(max(r_base(:,:,iroi)));
    r_single = reshape(r_norm(:,:,iroi),230*246,1);
    r_order = sort(r_single);
    bottom = mean(r_order(1:(230*246*.5)));
    top = mean(r_order(ceil(230*246*.99):end));
    if bottom./top < 0.5;
        n = 1+n;
        r_sig(:,:,n) = r_norm(:,:,iroi);
    end;
end;
r_signal=zeros(240,256,n);
r_signal(b+1:240-b,b+1:256-b,:)=r_sig;

figure;
for iroi = 1:size(r_signal,3);
    subplot(ceil(sqrt(size(r_signal,3))),ceil(sqrt(size(r_signal,3))),iroi);
    imagesq(r_signal(:,:,iroi));
end


%consolidate regions of interest
clear axons
r_copy = r_signal;
for iroi = 1:size(r_signal,3);
    clear r_corr
    for iroi2 = 1:size(r_copy,3);
        r_corr(iroi, iroi2) = triu2vec(corrcoef(r_copy(:,:,1),r_copy(:,:,iroi2)));    
    end;
    same = find(r_corr(iroi,:)>0.8);
    clear r_corr
    axons(:,:,iroi) = mean(r_copy(:,:,same),3);
    figure;
    subplot(5,5,1);
    imagesq(axons(:,:,iroi));
    for icomp = 1:length(same);
        subplot(5,5,icomp+1);
        imagesq(r_copy(:,:,same(icomp)));
    end
    for iroi2 = 1:size(r_copy,3);
        r_corr(iroi, iroi2) = triu2vec(corrcoef(axons(:,:,iroi),r_copy(:,:,iroi2)));
    end
    same2 = find(r_corr(iroi,:)>0.8);
    axons(:,:,iroi) = mean(r_copy(:,:,same2),3);
    figure;
    subplot(5,5,1);
    imagesq(axons(:,:,iroi));
    for icomp = 1:length(same2);
        subplot(5,5,icomp+1);
        imagesq(r_copy(:,:,same2(icomp)));
    end
    for iroi2 = 1:size(r_copy,3);
        r_corr(iroi, iroi2) = triu2vec(corrcoef(axons(:,:,iroi),r_copy(:,:,iroi2)));
    end
    r_copy(:,:,same2)=[];
end;

figure;
for iaxon = 1:size(axons,3);
    subplot(ceil(sqrt(size(axons,3))),ceil(sqrt(size(axons,3))), iaxon);
    imagesq(axons_norm(:,:,iaxon));
end

fn_out = fullfile(base, mouse,date,[date '_' mouse '_' expt '_axons.m']);
save(fn_out,'axons');

%extract mask
axons_thresh = zeros(240*256,size(axons,3));
axons_dig= zeros(240*256,size(axons,3));
for iaxon = 1:size(axons,3);
    clear dig
    axons_norm(:,:,iaxon) = axons(:,:,iaxon)./max(max(axons(:,:,iaxon)));
    dig = find(axons_norm(:,:,iaxon)>0.45);
    axons_thresh(dig(1:end,:),iaxon) = 1;
end;
axons_dig = reshape(axons_thresh, [240 256 size(axons,3)]);
figure;
for iaxon = 1:size(axons_dig,3);
    subplot(ceil(sqrt(size(axons_dig,3))),ceil(sqrt(size(axons_dig,3))), iaxon);
    imagesq(axons_dig(:,:,iaxon));
end;

axons_unique = zeros(240*256,size(axons_dig,3));
for iaxon = 1:size(axons_dig,3);
    clear unique
    axons_each(:,:,iaxon) = axons_dig(:,:,iaxon)-(sum(axons_dig,3)-axons_dig(:,:,iaxon));
    unique = find(axons_each(:,:,iaxon) == 1); 
    axons_unique(unique,iaxon) = iaxon;
end
axons_unique = reshape(axons_unique,240,256,size(axons_dig,3));
for iaxon = 1:size(axons_dig,3);
    subplot(ceil(sqrt(size(axons_dig,3))),ceil(sqrt(size(axons_dig,3))), iaxon);
    imagesq(axons_unique(:,:,iaxon));
end;
axons_all = sum(axons_unique,3);
figure;
imagesq(axons_all(:,:,1));

fn_out = fullfile(base, mouse,date,[date '_' mouse '_' expt '_axonmasks.m']);
save(fn_out,'axons_all');

timeCourses = stackGetTimeCourses(stack,axons_all);
figure;
tcOffsetplot(timeCourses(:,:));

%% hand segment
base = 'G:\users\lindsey\analysisLG\active mice';
reg = 'G:\users\lindsey\regLG';
date = '110517';
mouse = 'Y13';
userun = [6];

%make regions of interest
points = a(:,2:3);
n = length(points);
npix = 9;
roi = zeros(length(points)*npix,2);
for ipoint = 1:n;
    for ix = 1:3;
        for iy =1:3;
            roi(1+iy-1+(3*(ix-1))+((ipoint-1)*9),:) = [points(ipoint,2)-1+(iy-1),points(ipoint,1)-1+(ix-1)];
        end
    end
end

fn_out = fullfile(base, mouse, date, [date '_' mouse '_run' num2str(userun) '_manualrois.mat']);
save(fn_out, 'roi');

%display regions of interest
figure;
ind = zeros(npix*4,n/4);
subs = ceil(sqrt(n/4));
for iroi = 1:n/4
    for ipix = 1:36
        ind(ipix,iroi) = sub2ind([240 256],roi(ipix+((iroi-1)*36),1),roi(ipix+((iroi-1)*36),2));
    end
    fov = zeros(240,256);
    a = ind(:,iroi);
    fov(a)=1;
    subplot(subs,subs,iroi);    
    imagesq(fov);
    colormap('hot');
end

%extract time courses
fn = fullfile(base, mouse, date, [date '_' mouse '_run' num2str(userun) 'stack_sorted.tif']);
stack = readtiff(fn);
for iroi = 1:n/4;
    for ipix = 1:36;
        rois(ipix,iroi,:) = stack(roi(ipix+((iroi-1)*36),1),roi(ipix+((iroi-1)*36),2),:);
    end
end

%
roi_avg = mean(rois,1);
roi_avg_all = mean(roi_avg,3);
roi_dF = bsxfun(@minus, roi_avg, roi_avg_all);
roi_dF = squeeze(roi_dF)';
tt = [0:(1/2.67):(1/2.67)*size(roi_avg,3)];
figure;
tcOffsetPlot(roi_dF);

b = 5;
for iroi = 1:n/4;
    for iy = b+1:240-b
        fprintf('.');
        for ix = b+1:256-b
            sub = stack(iy-3:iy+3,ix-3:ix+3,:);
            sub_avg = mean(mean(sub,2),1);
            r(iy,ix,iroi)= triu2vec(corrcoef(roi_avg(:,iroi),sub_avg));
        end;
    end;
end;

figure;
ind = 1;
for iroi = 1:n/4;
    subplot(subs,subs,ind); 
    imagesq(r(:,:,iroi));
    colormap('hot');
    ind =ind+1;
    colorbar;
end

fn_out = fullfile(base, mouse, date, [date '_' mouse '_run' num2str(userun) 'uncorr_manualroi_localcorrelations.mat']);
save(fn_out, 'r');

%% resort data and remove signal correlations
nON = 12;
nOFF = 12;
nRep = 10;
nCond = 8;
nPlanes = 1;

fn = fullfile(base, mouse, date, [date '_' mouse '_run' num2str(userun) '_dec_reg.tif']);
stack = readtiff(fn);

seqfile = [date '_' mouse '_run' num2str(userun) '_Big_Seqposition.mat'];
load(fullfile(base, mouse, date,'analysis',seqfile));

outDir = fullfile(base,mouse,date,'analysis');
stack_sort

siz = size(stack_sorted);
stack_avg = zeros(siz(1), siz(2), nON+nOFF, nCond+1);

start = 0;
for iCond = 1:nCond+1;
    for iFrame = 1:(nOFF+nON)
    stack_sorted_avg(:,:,iFrame,iCond) = mean(stack_sorted(:,:,iFrame+start:nOFF+nON:start+((nOFF+nON)*stim_reps(iCond))),3);
    end
    start = (nOFF+nON)*stim_reps(iCond)+start;
end
stack_sorted_avg = uint16(stack_sorted_avg);

stack_uncorr = zeros(siz,'uint16');
start = 1;
for iCond = 1:nCond;
    nRep = stim_reps(iCond);
    for it = 1:nRep
        for iframe = 1:nOFF+nON;
            stack_uncorr(:,:,start) = stack_sorted(:,:,start)- stack_sorted_avg(:,:,iframe,iCond);
            start = start+1;
        end
    end
end


blank_avg = zeros(240, 256, 1);
blank_avg = mean(stack_sorted(:,:,start:end),3);
blank_avg = uint16(blank_avg);

nblanks = stim_reps(end);
for it = 1:nblanks*(nOFF+nON);
    stack_uncorr(:,:,start) = stack_sorted(:,:,start)- blank_avg;
    start = start+1;
end

fn_out = fullfile(base, mouse, date, 'analysis', [date '_' mouse '_run' num2str(userun) '_stack_uncorr.tif']);
writetiff(stack_uncorr, fn_out)

%% Tuning curves
fn = fullfile(base,mouse,date,'analysis',[date '_' mouse '_run' num2str(userun) '_peak_resp.mat']);
Var = 3;
for iS = 1:length(Seqposition)
    vars(iS)=Seqposition(iS).TFSFetc(Var);
end
uvars = unique(vars);
nvars = length(uvars);

figure;
subs = ceil(sqrt(size(roi_dF,2)));
for iRoi = 1:size(roi_dF,2);
    subplot(subs, subs, iRoi);
    errorbar(uvars, resp_avg(iRoi,1:nCond),resp_std(iRoi,1:nCond));
    hold on
    errorbar(-1,resp_avg(iRoi,end),resp_std(iRoi,end),'*k')
end

