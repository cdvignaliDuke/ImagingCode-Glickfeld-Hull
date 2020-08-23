clc; clear all; close all;
doRedChannel = 0;
ds = 'CrossOriRandPhase_ExptList';
eval(ds)
rc = behavConstsAV;
frame_rate = 30;
nexp = size(expt,2);
dspace_all = [];
dphase_all = [];
%%
for iexp = 1:nexp
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    area = expt(iexp).img_loc{1};
    ImgFolder = expt(iexp).coFolder;
    time = expt(iexp).coTime;
    nrun = length(ImgFolder);
    run_str = catRunName(cell2mat(ImgFolder), nrun);

    LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';

    fprintf([mouse ' ' date '\n'])

    %% load data

    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']))
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_phaseFits.mat']))
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']))

    %% 
    mod_ind = find(p_anova_all<0.05);
    phas_mask = zeros(size(mask_cell));
    nCells = size(amp_hat_all,1);
    pha = rad2deg(pha_hat_all);
    while find(pha<0)
        pha(find(pha<0)) = pha(find(pha<0))+360;
    end
    while find(pha>360)
        pha(find(pha>360)) = pha(find(pha>360))-360;
    end
    for iCell = 1:nCells
        if find(mod_ind == iCell)
            ind = find(mask_cell==iCell);
            phas_mask(ind) = pha(iCell,:);
        end
    end
    figure;
    imagesc(phas_mask)
    colorbar
    truesize
    title([mouse ' ' date])
    print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_phaseMap.pdf']),'-dpdf','-bestfit')

    stats = regionprops(mask_cell);
    centroids = reshape([stats.Centroid],[nCells 2]);
    nMod = length(mod_ind);
    d_space = nan(nMod);
    d_phase = nan(nMod);
    for i = 1:nMod
        iC = mod_ind(i);
        for j = 1:nMod
            jC = mod_ind(j);
            if j>i
                d_space(i,j) = pdist([centroids(iC,:); centroids(jC,:)], 'euclidean');
                d_phase(i,j) = circ_dist(pha_hat_all(iC,:),pha_hat_all(jC,:));
            end
        end
    end
    dphase_temp = reshape(rad2deg(abs(d_phase)),[nMod.*nMod 1]);
    dphase = dphase_temp;
    dphase(isnan(dphase_temp)) = [];
    dspace = reshape(d_space,[nMod.*nMod 1]);
    dspace(isnan(dphase_temp)) = [];

    figure; 
    subplot(2,2,1)
    scatter(dspace,dphase)
    xlabel('Spatial distance (pix)')
    ylabel('Phase distance (deg)')

    [n edges bins] = histcounts(dspace,5);
    phase_avg = zeros(2,length(n));
    space_avg = zeros(2,length(n));
    for i = 1:length(n)
        ind = find(bins==i);
        phase_avg(1,i) = mean(dphase(ind,:),1);
        space_avg(1,i) = mean(dspace(ind,:),1);
        phase_avg(2,i) = std(dphase(ind,:),[],1)./sqrt(n(i));
        space_avg(2,i) = std(dspace(ind,:),[],1)./sqrt(n(i));
    end
    subplot(2,2,2)
    errorbar(space_avg(1,:), phase_avg(1,:), phase_avg(2,:),phase_avg(2,:),space_avg(2,:),space_avg(2,:))
    xlabel('Spatial distance (pix)')
    ylabel('Phase distance (deg)')
    ylim([0 200])
    print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_phaseDist.pdf']),'-dpdf','-bestfit')

    save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_phaseDist.mat']),'dphase','dspace')
    dphase_all = [dphase_all; dphase];
    dspace_all = [dspace_all; dspace];
end
summaryDir = fullfile(LG_base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary');

figure;
subplot(2,2,1)
scatter(dspace_all,dphase_all)
xlabel('Spatial distance (pix)')
ylabel('Phase distance (deg)')

[n edges bins] = histcounts(dspace_all,5);
phase_all_avg = zeros(2,length(n));
space_all_avg = zeros(2,length(n));
for i = 1:length(n)
    ind = find(bins==i);
    phase_all_avg(1,i) = mean(dphase_all(ind,:),1);
    space_all_avg(1,i) = mean(dspace_all(ind,:),1);
    phase_all_avg(2,i) = std(dphase_all(ind,:),[],1)./sqrt(n(i));
    space_all_avg(2,i) = std(dspace_all(ind,:),[],1)./sqrt(n(i));
end
subplot(2,2,2)
errorbar(space_avg(1,:), phase_avg(1,:), phase_avg(2,:),phase_avg(2,:),space_avg(2,:),space_avg(2,:))
xlabel('Spatial distance (pix)')
ylabel('Phase distance (deg)')
ylim([0 200])

print(fullfile(summaryDir, 'phaseDistSummary.pdf'),'-dpdf','-bestfit')

