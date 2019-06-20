%%mouse and cell IDs for axon imaging
rc = behavConstsAV;
mouse = '1202';
expDate = '190110';
dataFolder = '003';
i = 3;
ii = 1500;

fnout = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate);
load(fullfile(fnout, dataFolder, [mouse '_' expDate '_input.mat']))
load(fullfile(fnout, dataFolder, [mouse '_' expDate '_lbub_fits.mat']))
load(fullfile(fnout, dataFolder, [mouse '_' expDate '_Fit_struct_sub.mat']))
load(fullfile(fnout, dataFolder, [mouse '_' expDate '_mask_cell.mat']))
load(fullfile(fnout, dataFolder, [mouse '_' expDate '_TCs.mat']))
%%
Az = celleqel2mat_padded(input.tGratingAzimuthDeg);
El = celleqel2mat_padded(input.tGratingElevationDeg);
Azs = unique(Az);
Els = unique(El);
nAz = length(Azs);
nEl = length(Els);
dAz = median(diff(Azs));
dEl = median(diff(Els));
Az_vec00 = Azs(1):(dAz/10):Azs(end);
El_vec00 = Els(1):(dEl/10):Els(end);
[AzAz00,ElEl00]=meshgrid(Az_vec00,El_vec00);
[AzAz, ElEl] = meshgrid(Azs,Els);
grid2.AzAz = AzAz;
grid2.ElEl = ElEl;
grid2.AzAz00 = AzAz00;
grid2.ElEl00 = ElEl00;

%% example time courses
nOff = input.nScansOff;
nOn = input.nScansOn;
tframes = size(data_tc,1);
ntrials = length(Az);
nCells = size(data_tc,2);
frameRateHz = 15.5;
nframes = ceil(frameRateHz.*2);
baseWin = ceil(frameRateHz.*0.5);
trial_tc = nan(nframes,nCells,ntrials-1);
stimOn = nOff+1:nOff+nOn:tframes;
for it = 1:ntrials-1
	trial_tc(:,:,it) =data_tc(stimOn(it)-baseWin:stimOn(it)+nframes-baseWin-1,:);
end
trial_f = mean(trial_tc(1:baseWin,:,:),1);
trial_dfof = (trial_tc-trial_f)./trial_f;

tt = (1-baseWin:nframes-baseWin).*(1000/frameRateHz);
figure;
iC = goodfit_ind(i);
Els_flip = fliplr(Els);
start = 1;
for iEl = 1:nEl
    ind_El = find(El(1:ntrials-1) == Els_flip(iEl));
    for iAz = 1:nAz
        ind_Az = find(Az(1:ntrials-1) == Azs(iAz));
        ind = intersect(ind_El,ind_Az);
        subplot(7,7,start)
        shadedErrorBar(tt',mean(trial_dfof(:,iC,ind),3),std(trial_dfof(:,iC,ind),[],3)./sqrt(length(ind)),'r');
        ylim([-0.2 0.6])
        axis square
        if start<49
        axis off
        end
        start = start+1;
    end
end
suptitle(['i1202 190110 - Cell#' num2str(iC)])
print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Manuscripts\sizeTuning\Figures5-6\i1202_190110_Cell#' num2str(iC) '_retTC.pdf'],'-dpdf','-bestfit')

figure;
iC = goodfit_ind(ii);
Els_flip = fliplr(Els);
start = 1;
for iEl = 1:nEl
    ind_El = find(El(1:ntrials-1) == Els_flip(iEl));
    for iAz = 1:nAz
        ind_Az = find(Az(1:ntrials-1) == Azs(iAz));
        ind = intersect(ind_El,ind_Az);
        subplot(7,7,start)
        shadedErrorBar(tt',mean(trial_dfof(:,iC,ind),3),std(trial_dfof(:,iC,ind),[],3)./sqrt(length(ind)),'b');
        ylim([-0.2 0.6])
        axis square
        axis off
        start = start+1;
    end
end
suptitle(['i1202 190110 - Cell#' num2str(iC)])
print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Manuscripts\sizeTuning\Figures5-6\i1202_190110_Cell#' num2str(iC) '_retTC.pdf'],'-dpdf','-bestfit')


    
%% example RFs and fits
c_red = zeros(64,3);
c_red(:,1) = 0:1/63:1;
c_blue = zeros(64,3);
c_blue(:,3) = 0:1/63:1;

figure; 
iC = goodfit_ind(i);
subplot(2,2,1)
x  = Fit_struct_sub(iC).True.s_.data;
imagesq(x);
title(num2str(chop(max(x(:)),3)))
colormap gray
subplot(2,2,2)
y  = Fit_struct_sub(iC).True.s_.k2_plot_oversamp;
imagesq(y);
hold on
a = interp1(grid2.AzAz00(1,:),1:61,lbub_fits(iC,4,4));
e = interp1(flipud(grid2.ElEl00(:,1)),1:61,lbub_fits(iC,5,4));
a1 = interp1(grid2.AzAz00(1,:),1:61,lbub_fits(iC,2,4));
e1 = interp1(flipud(grid2.ElEl00(:,1)),1:61,lbub_fits(iC,3,4));
azero = interp1(grid2.AzAz00(1,:),1:61,0);
ezero = interp1(flipud(grid2.ElEl00(:,1)),1:61,0);
scatter(a, e, 'ok')
ellipse(abs(a1-azero), abs(e1-ezero), 0, a, e,'r');
colormap gray
suptitle(['i1202 190110 - Cell#' num2str(iC)])
print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Manuscripts\sizeTuning\Figures5-6\i1202_190110_Cell#' num2str(iC) '_ret.pdf'],'-dpdf','-bestfit')

figure; 
iC = goodfit_ind(ii);
subplot(2,2,1)
x  = Fit_struct_sub(iC).True.s_.data;
imagesq(x);
title(num2str(chop(max(x(:)),3)))
colormap gray
subplot(2,2,2)
y  = Fit_struct_sub(iC).True.s_.k2_plot_oversamp;
imagesq(y);
hold on
a = interp1(grid2.AzAz00(1,:),1:61,lbub_fits(iC,4,4));
e = interp1(flipud(grid2.ElEl00(:,1)),1:61,lbub_fits(iC,5,4));
a1 = interp1(grid2.AzAz00(1,:),1:61,lbub_fits(iC,2,4));
e1 = interp1(flipud(grid2.ElEl00(:,1)),1:61,lbub_fits(iC,3,4));
azero = interp1(grid2.AzAz00(1,:),1:61,0);
ezero = interp1(flipud(grid2.ElEl00(:,1)),1:61,0);
scatter(a, e, 'ok')
ellipse(abs(a1-azero), abs(e1-ezero), 0, a, e,'b');
colormap gray
suptitle(['i1202 190110 - Cell#' num2str(iC)])
print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Manuscripts\sizeTuning\Figures5-6\i1202_190110_Cell#' num2str(iC) '_ret.pdf'],'-dpdf','-bestfit')

%% ret for field of view
fprintf('Making bouton response map \n')
mask_label = bwlabel(mask_all);
retMap_El = NaN(size(mask_cell));
retMap_Az = retMap_El;
for i=1:length(goodfit_ind)
    ind = find(mask_label == goodfit_ind(i));
    retMap_El(ind) = lbub_fits(goodfit_ind(i),5,4);
    retMap_Az(ind) = lbub_fits(goodfit_ind(i),4,4);
end

imAlpha=ones(size(retMap_El));
imAlpha(isnan(retMap_El))=0; % set all unmasked pixels to alpha=0

figure;clf;
colormap jet
imagesc(retMap_El,'AlphaData',imAlpha)
title('Retinotopy of goodfit cells by El')
h = colorbar;
ylabel(h,'El (deg)','Rotation',270.0,'VerticalAlignment','bottom')
set(gca,'color',0*[1 1 1]);
truesize
print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Manuscripts\sizeTuning\Figures5-6\' mouse '_' expDate '_retMapEl.pdf'], '-dpdf','-bestfit')


figure;clf;
colormap jet
imagesc(retMap_Az,'AlphaData',imAlpha)
title('Retinotopy of goodfit cells by Az')
h = colorbar;
ylabel(h,'Az (deg)','Rotation',270.0,'VerticalAlignment','bottom')
set(gca,'color',0*[1 1 1]); 
truesize

print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Manuscripts\sizeTuning\Figures5-6\' mouse '_' expDate '_retMapAz.pdf'], '-dpdf','-bestfit')
