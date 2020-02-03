mouse = 'i1312';
date = '200120';
ImgFolder = strvcat('003');
ref_date = '200118';
ref_run = strvcat('002');
nrun = size(ImgFolder,1);
frame_rate = 15.5;
run_str = catRunName(ImgFolder, nrun);
ref_str = catRunName(ref_run, size(ref_run,1));
fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P';

%% 
oriTuning_D1 = load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_oriTuningAndFits.mat']));
oriTuning_D2 = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_oriTuningAndFits.mat']));

%%
goodfit_D1 = find(oriTuning_D1.fitReliability<22.5);
goodfit_D2 = find(oriTuning_D2.fitReliability<22.5);
goodfit = intersect(goodfit_D1,goodfit_D2);

[maxResp_D1 prefOri_D1] = max(squeeze(oriTuning_D1.vonMisesFitAllCellsAllBoots(:,1,:)),[],1);
[maxResp_D2 prefOri_D2] = max(squeeze(oriTuning_D2.vonMisesFitAllCellsAllBoots(:,1,:)),[],1);

figure; 
subplot(2,2,1)
scatter(maxResp_D1(goodfit),maxResp_D2(goodfit))
axis square
xlim([0 4])
ylim([0 4])
refline(1,0)
xlabel('Day 1 max dF/F')
ylabel('Day 2 max dF/F')

prefOri_diff = abs(prefOri_D1(goodfit)-prefOri_D2(goodfit));
prefOri_diff(find(prefOri_diff>90)) = 180-prefOri_diff(find(prefOri_diff>90));
subplot(2,2,2)
hist(prefOri_diff)
xlabel('D1 vs D2 Pref ori')
ylabel('Number of cells')

subplot(2,2,3)
scatter(maxResp_D1(goodfit),prefOri_diff)
xlabel('Day 1 max dF/F')
ylabel('D1 vs D2 Pref ori')
axis square

maxResp_diff = abs((maxResp_D1(goodfit)-maxResp_D2(goodfit))./(maxResp_D1(goodfit)+maxResp_D2(goodfit)));
subplot(2,2,4)
scatter(maxResp_diff,prefOri_diff)
xlabel('D1 vs D2 max dF/F')
ylabel('D1 vs D2 Pref ori')
axis square

suptitle([mouse ' ' ref_date ' vs ' date '- n = ' num2str(length(goodfit)) ' cells'])
print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_compAcrossDaysTuning.pdf']),'-dpdf', '-bestfit')
save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_compAcrossDaysTuning.mat']), 'goodfit', 'prefOri_D1', 'maxResp_D1', 'prefOri_D2', 'maxResp_D2', 'prefOri_diff', 'maxResp_diff');

%%
load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_transform.mat']));
load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']));

cell_stats = regionprops(mask_cell);
nCells = size(cell_stats,1);
%select region of image
%ROI = [250:450; 100:300]; %i1313
ROI = [150:350; 300:500]; %i1312
cell_list = intersect(goodfit, unique(mask_cell(ROI(1,:),ROI(2,:))));
figure;
subplot(2,2,1)
imagesc(ref);
for iC = 1:length(cell_list)
    iCell = cell_list(iC);
    if length(find(mask_cell == iCell)) == length(find(mask_cell(ROI(1,:),ROI(2,:)) == iCell))
        text(cell_stats(iCell).Centroid(1), cell_stats(iCell).Centroid(2), num2str(iCell), 'Color', 'white',...
            'Fontsize', 6, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle') 
        hold on
    else
        cell_list(iC) = NaN;
    end
end
xlim([ROI(2,1) ROI(2,end)])
ylim([ROI(1,1) ROI(1,end)])
axis square
axis off
title('Day 1')

subplot(2,2,2)
imagesc(reg2ref);
for iC = 1:length(cell_list)
    iCell = cell_list(iC);
    if ~isnan(iCell)
        text(cell_stats(iCell).Centroid(1), cell_stats(iCell).Centroid(2), num2str(iCell), 'Color', 'white',...
            'Fontsize', 6, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle') 
        hold on
    end
end
xlim([ROI(2,1) ROI(2,end)])
ylim([ROI(1,1) ROI(1,end)])
axis square
axis off
title('Day 2')

subplot(2,2,3)
imagesc(data_dfof_max_ref);
for iC = 1:length(cell_list)
    iCell = cell_list(iC);
    if ~isnan(iCell)
        text(cell_stats(iCell).Centroid(1), cell_stats(iCell).Centroid(2), num2str(iCell), 'Color', 'white',...
            'Fontsize', 6, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle') 
        hold on
    end
end
xlim([ROI(2,1) ROI(2,end)])
ylim([ROI(1,1) ROI(1,end)])
axis square
axis off

subplot(2,2,4)
imagesc(reg2ref_dfof);
for iC = 1:length(cell_list)
    iCell = cell_list(iC);
    if ~isnan(iCell)
        text(cell_stats(iCell).Centroid(1), cell_stats(iCell).Centroid(2), num2str(iCell), 'Color', 'white',...
            'Fontsize', 6, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle') 
        hold on
    end
end
xlim([ROI(2,1) ROI(2,end)])
ylim([ROI(1,1) ROI(1,end)])
axis square
axis off

print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_cellAlignZoom.pdf']),'-dpdf', '-bestfit')

nC = sum(~isnan(cell_list));
figure;
start = 1;
for iC = 1:length(cell_list)
    if ~isnan(cell_list(iC))
        iCell = cell_list(iC);
        subplot(5, 6, start)
        plot(0:180, oriTuning_D1.vonMisesFitAllCellsAllBoots(:,1,iCell))
        hold on
        errorbar(0:22.5:180-22.5, oriTuning_D1.avgResponseEaOri(iCell,:), oriTuning_D1.semResponseEaOri(iCell,:))
        title(['Cell ' num2str(iCell) '- ' num2str(prefOri_D1(iCell)) ' deg'])
        subplot(5, 6, start+1)
        plot(0:180, oriTuning_D2.vonMisesFitAllCellsAllBoots(:,1,iCell))
        hold on
        errorbar(0:22.5:180-22.5, oriTuning_D2.avgResponseEaOri(iCell,:), oriTuning_D2.semResponseEaOri(iCell,:))
        title(['Cell ' num2str(iCell) '- ' num2str(prefOri_D2(iCell)) ' deg'])
        start = start+2;
    end
end

print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_cellAlignZoomTuning.pdf']),'-dpdf', '-bestfit')

    
    