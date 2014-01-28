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
%make regions of interest
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

%display regions of interest
figure;
ind = zeros(npix*4,n/4);
for iroi = 1:n/4
    for ipix = 1:36
        ind(ipix,iroi) = sub2ind([240 256],roi(ipix+((iroi-1)*36),1),roi(ipix+((iroi-1)*36),2));
    end
    fov = zeros(240,256);
    a = ind(:,iroi);
    fov(a)=1;
    subplot(2,2,iroi);    
    imagesq(fov);
    colormap('hot');
end

%extract time courses
for iroi = 1:n/4;
    for ipix = 1:36;
        rois(ipix,iroi,:) = stack_d12(roi(ipix+((iroi-1)*36),1),roi(ipix+((iroi-1)*36),2),:);
    end
end

%
roi_avg = mean(rois,1);
b = 5;
for iroi = 1:n/4;
    for iy = b+1:240-b
        fprintf('.');
        for ix = b+1:256-b
            sub = stack_d6(iy-3:iy+3,ix-3:ix+3,:);
            sub_avg = mean(mean(sub,2),1);
            r(iy,ix,iroi)= triu2vec(corrcoef(roi_avg(:,iroi,:),sub_avg));
        end;
    end;
end;

figure;
ind = 1;
for iroi = 1:n/4;
    subplot(2,2,ind); 
    imagesq(r(:,:,iroi));
    colormap('hot');
    ind =ind+1;
    colorbar;
end

fn_out = fullfile(base, mouse,date,[date '_' mouse '_' expt 'manualrois.mat']);
save(fn_out,'roi','r');

figure;
tc = tcDecimate(squeeze(roi_avg),2);
base = mean(tc,2);
tc_base = bsxfun(@minus,tc,base);
tcoffsetplot(tc_base');

%% resort data and remove signal correlations
nON = 144;
nOFF = 144;
nRep = 10;
nCond = 8;
nPlanes = 1;
seqfile = [date '_' mouse '_' expt '_dec_Seqposition.mat'];
load(fullfile(base, mouse, date,'analysis',seqfile));
stack_sorted = zeros(siz,'uint16');
for iRep = 1:nRep;
    for iCond = 1:nCond;
        ind = Seqposition(iCond).ind(iRep);
        stack_sorted(:,:,1+((iRep-1)*nCond*((nOFF+nON)/nPlanes))+((iCond-1)*((nOFF+nON)/nPlanes)):((nOFF+nON)/nPlanes)+((iRep-1)*nCond*((nOFF+nON)/nPlanes))+((iCond-1)*((nOFF+nON)/nPlanes))) = stack(:,:,1+(ind-1)*((nOFF+nON)/nPlanes):((nOFF+nON)/nPlanes)+(ind-1)*((nOFF+nON)/nPlanes));
    end
end
stack_sorted_avg = zeros(240,256,nCond*(nOFF+nON),'uint16');
for iframe = 1:nCond*(nOFF+nON)
    stack_sorted_avg(:,:,iframe) = mean(stack_sorted(:,:,iframe:nCond*(nOFF+nON):nCond*(nOFF+nON)*nRep),3);
end

stack_uncorr = zeros(siz,'uint16');
for irep = 1:nRep;
    for iframe = 1:nCond*(nOFF+nON);
        stack_uncorr(:,:,iframe+((irep-1)*(nCond*(nOFF+nON)))) = stack_sorted(:,:,iframe+((irep-1)*(nCond*(nOFF+nON))))- stack_sorted_avg(:,:,iframe);
    end
end

stack_avg = zeros(240, 256, 1);
stack_avg = mean(stack(:,:,23041:end),3);
stack_avg = uint16(stack_avg);

nblanks = length(Seqposition(end).ind);
for iblank = 1:nblanks;
    ind = Seqposition(end).ind(iblank);
    blanks(:,:,1+(iblank-1)*((nOFF+nON)/nPlanes):((nOFF+nON)/nPlanes)+(iblank-1)*((nOFF+nON)/nPlanes)) = stack(:,:,1+(ind-1)*((nOFF+nON)/nPlanes):((nOFF+nON)/nPlanes)+(ind-1)*((nOFF+nON)/nPlanes));
end
blanks_dF = bsxfun(@minus,blanks,stack_avg);
stack_uncorr(:,:,23041:end) = blanks_dF;

fn = fullfile(reg, [date '_' mouse], [expt '_uncorr'], [expt '_uncorr_green']);
writetiffseq(stack_uncorr, fn, 'uint16',216);

fn = fullfile(reg, [date '_' mouse], [expt '_uncorr'], [expt '_uncorr_green']);
stack_uncorr = readtiff(fn);
stack_uncorr_d6 = stackGroupProject(stack_uncorr,6);