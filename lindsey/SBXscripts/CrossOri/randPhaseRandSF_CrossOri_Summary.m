clc; clear all; close all;
doRedChannel = 0;
ds = 'CrossOriRandPhase_randSF_ExptList';
eval(ds)
rc = behavConstsAV;
frame_rate = 30;
nexp = size(expt,2);
seed = rng;
nSF =5;

LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
summaryDir = fullfile(LG_base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary');

p_anova_all = [];
plaidSI_all = [];
testPI_all = [];
RsqSF_all = [];
prefSF_all = [];
b_hat_all = [];
amp_hat_all = [];
amp_hat_shuf_all = [];
resp_ind_all = [];
max_sf_ind_all = [];

totCells = 0;

for iexp = 1:nexp
    
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    area = expt(iexp).img_loc{1};
    ImgFolder = expt(iexp).coFolder;
    nrun = length(ImgFolder);
    run_str = catRunName(cell2mat(ImgFolder), nrun);

    fprintf(['Mouse: ' mouse ' Date: ' date '\n'])

    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dataStim.mat']))
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']))
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']))
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_pupil.mat']))
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_phaseFits.mat']))
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_SfFits.mat']))
    if iexp == 1
        resp_ind_sf_all = cell(1,nSF);
    end
    
    p_anova_all = [p_anova_all; p_anova];
    plaidSI_all = [plaidSI_all; plaidSI];
    testPI_all = [testPI_all; testPI];
    RsqSF_all = [RsqSF_all; RsqSF'];
    prefSF_all = [prefSF_all; prefSF'];
    b_hat_all = [b_hat_all; b_hat];
    amp_hat_all = [amp_hat_all; amp_hat];
    amp_hat_shuf_all = [amp_hat_shuf_all; amp_hat_shuf];
    resp_ind_all = [resp_ind_all; resp_ind+totCells];
    max_sf_ind_all = [max_sf_ind_all; max_sf_ind'];
    for isf = 1:nSF
        resp_ind_sf_all{isf} = [resp_ind_sf_all{isf}; resp_ind_sf{isf}+totCells];
    end
    totCells = totCells+size(p_anova,1);
end
    
    figure;
    movegui('center');
    legstr = [];
    for isf = 1:nSF
        subplot(1,2,1)
        cdfplot(plaidSI_all(resp_ind_sf_all{isf},isf))
        xlim([-1 1])
        hold on
        subplot(1,2,2)
        cdfplot(testPI_all(resp_ind_sf_all{isf},isf))
        xlim([0 1])
        hold on
        legstr = [legstr; [num2str(SFs(isf)) '- n = ' num2str(length(resp_ind_sf_all{isf}))]];
    end
    subplot(1,2,1)
    xlabel('Suppression index')
    subplot(1,2,2)
    legend(legstr,'location', 'northwest')
    xlabel('Preference index')
    print(fullfile(summaryDir, 'randPhaseRandSF_SuppAndPrefIndex_cdfs.pdf'),'-dpdf','-fillpage');
    
    figure;
    movegui('center');
    legstr1 = [];
    legstr2 = [];
    for isf = 1:nSF
        subplot(1,2,1)
        cdfplot(testPI_all(resp_ind_sf_all{isf},isf))
        legstr1 = [legstr1; [num2str(SFs(isf)) '- n = ' num2str(length(resp_ind_sf_all{isf}))]];
        xlim([0 1])
        hold on
        subplot(1,2,2)
        cdfplot(testPI_all(intersect(find(max_sf_ind_all==isf),resp_ind_sf_all{isf}),isf))
        xlim([0 1])
        hold on
        legstr2 = [legstr2; [num2str(SFs(isf)) '- n = ' num2str(length(intersect(find(max_sf_ind_all==isf),resp_ind_sf_all{isf})))]];
    end
    subplot(1,2,1)
    xlabel('Preference index')
    legend(legstr1,'location', 'northwest')
    title('All')
    subplot(1,2,2)
    legend(legstr2,'location', 'northwest')
    xlabel('Preference index')
    title('Peak')
    print(fullfile(summaryDir, 'randPhaseRandSF_PrefIndex_cdfs.pdf'),'-dpdf','-fillpage');

    figure; movegui('center');
    RsqSF_all(find(RsqSF_all<0)) = 0;
    subplot(2,1,1); histogram(RsqSF_all(resp_ind_all)); vline(0.8)
    ind = intersect(find(RsqSF_all>0.8),resp_ind_all); xlabel('Rsq')
    subplot(2,1,2); histogram(prefSF_all(ind),nSF)
    set(gca, 'XTick', log2(SFs), 'XTickLabels', SFs)
    vline(mean(prefSF_all(ind),1))
    xlabel('SF (cpd)')
    print(fullfile(summaryDir, 'randPhaseRandSF_SFprefs.pdf'),'-dpdf','-fillpage');

    figure;  movegui('center');
    amp_diff_all = amp_hat_all-amp_hat_shuf_all;
    subplot(2,2,1)
    legstr = [];
    for isf = 1:nSF
        if length(resp_ind_sf_all{isf})>1
            if sum(~isnan(amp_diff_all(resp_ind_sf_all{isf},isf)),1)>1
            cdfplot(amp_diff_all(resp_ind_sf_all{isf},isf))
            hold on
            legstr = [legstr; [num2str(SFs(isf)) '- n = ' num2str(length(resp_ind_sf_all{isf}))]];
            end
        end
    end
    xlabel('Sine Amp: Data-Shuf')
    legend(legstr,'location', 'northwest')
    title('All responsive')
    subplot(2,2,2)
    legstr = [];
    for isf = 1:nSF
        ind = intersect(find(max_sf_ind_all==isf),resp_ind_sf_all{isf});
        if length(ind)>1
            if sum(~isnan(amp_diff_all(ind,isf)),1)>1
            cdfplot(amp_diff_all(ind,isf))
            hold on
            legstr = [legstr; [num2str(SFs(isf)) '- n = ' num2str(length(ind))]];
            end
        end
    end
    xlabel('Sine Amp: Data-Shuf')
    legend(legstr,'location', 'northwest')
    title('Peak SF')
    subplot(2,2,3)
    legstr = [];
    for isf = 1:nSF
        if length(resp_ind_sf_all{isf})>1
            if sum(~isnan(b_hat_all(resp_ind_sf_all{isf},isf)),1)>1
            cdfplot(b_hat_all(resp_ind_sf_all{isf},isf))
            hold on
            legstr = [legstr; [num2str(SFs(isf)) '- n = ' num2str(length(resp_ind_sf_all{isf}))]];
            end
        end
    end
    xlabel('Sine Baseline')
    subplot(2,2,4)
    legstr = [];
    for isf = 1:nSF
        ind = intersect(find(max_sf_ind_all==isf),resp_ind_sf_all{isf});
        if length(ind)>1
            if sum(~isnan(b_hat_all(ind,isf)),1)>1
            cdfplot(b_hat_all(ind,isf))
            hold on
            legstr = [legstr; [num2str(SFs(isf)) '- n = ' num2str(length(ind))]];
            end
        end
    end
    xlabel('Sine Baseline')
    
    print(fullfile(summaryDir, 'randPhaseRandSF_diffAmpByStimSF_cdfs.pdf'),'-dpdf','-bestfit');
    

    figure;  movegui('center');
    for isf = 1:nSF
        subplot(3,2,isf)
        for iCell = 1:totCells
            if find(resp_ind_sf_all{isf} == iCell)
                if p_anova_all(iCell,isf)<0.05
                    scatter(2.^prefSF_all(iCell), amp_diff_all(iCell,isf),'or')
                else
                    scatter(2.^prefSF_all(iCell), amp_diff_all(iCell,isf),'ok')
                end
                hold on
            end
        end
        ylim([-0.2 0.7])
        xlim([0.02 0.32])
        ylabel('Sine Amp: Data-Shuf')
        xlabel('Preferred SF')
        title(['Stim SF: ' num2str(SFs(isf))])
    end
    subplot(3,2,6)
    for iCell = 1:totCells
        if find(resp_ind_all == iCell)
            if p_anova_all(iCell,max_sf_ind_all(iCell))<0.05
                scatter(2.^prefSF_all(iCell), amp_diff_all(iCell,max_sf_ind_all(iCell)),'or')
            else
                scatter(2.^prefSF_all(iCell), amp_diff_all(iCell,max_sf_ind_all(iCell)),'ok')
            end
            hold on
        end
    end
    ylabel('Sine Amp: Data-Shuf')
    xlabel('Preferred SF')
    title('At Peak SF')
    ylim([-0.2 0.7])
    xlim([0.02 0.32])
    print(fullfile(summaryDir, 'randPhaseRandSF_diffAmpByPrefSF_scatters.pdf'),'-dpdf','-fillpage');
    
    figure;  movegui('center');
    for isf = 1:nSF
        subplot(3,2,isf)
        for iCell = 1:totCells
            if find(resp_ind_sf_all{isf} == iCell)
                scatter(2.^prefSF_all(iCell), b_hat_all(iCell,isf),'ok')
                hold on
            end
        end
        ylim([-1 1])
        xlim([0.02 0.32])
        ylabel('Sine Baseline')
        xlabel('Preferred SF')
        title(['Stim SF: ' num2str(SFs(isf))])
    end
    subplot(3,2,6)
    for iCell = 1:totCells
        if find(resp_ind_all == iCell)
            scatter(2.^prefSF_all(iCell), b_hat_all(iCell,max_sf_ind_all(iCell)),'ok')
            hold on
        end
    end
    ylabel('Sine Baseline')
    xlabel('Preferred SF')
    title('At Peak SF')
    ylim([-1 1])
    xlim([0.02 0.32])
    print(fullfile(summaryDir, 'randPhaseRandSF_baseByPrefSF_scatters.pdf'),'-dpdf','-fillpage');
