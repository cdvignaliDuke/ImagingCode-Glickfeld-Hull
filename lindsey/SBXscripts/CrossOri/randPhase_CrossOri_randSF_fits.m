clc; clear all; close all;
doRedChannel = 0;
ds = 'CrossOriRandPhase_randSF_ExptList';
eval(ds)
rc = behavConstsAV;
frame_rate = 30;
nexp = size(expt,2);
seed = rng;
for iexp = 2:nexp

    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    area = expt(iexp).img_loc{1};
    ImgFolder = expt(iexp).coFolder;
    time = expt(iexp).coTime;
    nrun = length(ImgFolder);
    run_str = catRunName(cell2mat(ImgFolder), nrun);

    LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';

    fprintf(['2P imaging sine fitting analysis\nSelected data:\nMouse: ' mouse '\nDate: ' date '\nExperiments:\n'])
    for irun=1:nrun
        fprintf([ImgFolder{irun} ' - ' time{irun} '\n'])
    end

    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dataStim.mat']))
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']))
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']))
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_pupil.mat']))
    
    nCells = size(resp_cell{end,end,end,end},1);
    nTrials = size(stimCon_all,2);
    trial_n = zeros(nMaskCon,nStimCon,nMaskPhas,nSF);
    trialInd = cell(nMaskCon,nStimCon,nMaskPhas,nSF);
    for isf = 1:nSF
        ind_sf = find(SF_all == SFs(isf));
        for im = 1:nMaskCon
            ind_mask = find(maskCon_all == maskCons(im));
            for it = 1:nStimCon
                ind_stim = find(stimCon_all == stimCons(it));
                ind_sm = intersect(ind_mask,ind_stim);
                if it>1 & im>1
                    for ip = 1:nMaskPhas
                        ind_phase = find(maskPhas_all == maskPhas(ip));
                        ind = intersect(ind_phase,intersect(ind_sf,ind_sm));
                        trialInd{im,it,ip,isf} = ind;
                        trial_n(im,it,ip,isf) = length(ind);
                    end
                else
                    trialInd{im,it,1,isf} = intersect(ind_sf,ind_sm);
                    trial_n(im,it,isf) = length(trialInd{im,it,1,isf});
                end
            end
        end
    end

    p_anova = nan(nCells,nSF);
    b_hat = nan(nCells,nSF); 
    amp_hat = nan(nCells,nSF); 
    per_hat = nan(nCells,nSF); 
    pha_hat = nan(nCells,nSF); 
    sse = nan(nCells,nSF); 
    R_square = nan(nCells,nSF);
    yfit = nan(nCells,length(0:1:359),nSF);
    p_anova_shuf = nan(nCells,nSF);
    b_hat_shuf = nan(nCells,nSF); 
    amp_hat_shuf = nan(nCells,nSF); 
    per_hat_shuf = nan(nCells,nSF); 
    pha_hat_shuf = nan(nCells,nSF); 
    sse_shuf = nan(nCells,nSF); 
    R_square_shuf = nan(nCells,nSF);
    yfit_shuf = nan(nCells,length(0:1:359),nSF);
            
    eye_n = nan(nSF,nMaskPhas);
    phase_range = 0:1:359;
    if ~exist('centroid_dist_sf')
        centroid_dist_sf = cell(1,nSF);
        centroid_med_sf = cell(1,nSF);
        ind = find(~isnan(centroid_stim(1,:)));
        centroid_med = findMaxNeighbors(centroid_stim(:,ind),2);
        for isf = 1:nSF
            ind_sf = intersect(find(SF_all == SFs(isf)),find(~isnan(centroid_stim(1,:))));
            centroid_med_sf{isf} = findMaxNeighbors(centroid_stim(:,ind_sf),2);
            centroid_dist_sf{isf} = sqrt((centroid_stim(1,:)-centroid_med(1)).^2 + (centroid_stim(2,:)-centroid_med(2)).^2);
        end
        save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_pupil.mat']), 'rect', 'Area', 'Centroid', 'SNR', 'Val', 'frame_rate' , 'rad_mat_start','centroid_mat_start', 'cStimOn', 'rad_base','rad_stim','centroid_base', 'centroid_stim', 'centroid_dist','centroid_med', 'centroid_dist_sf', 'centroid_med_sf');
    end
    
    %less than 2 deg
    for isf = 1:nSF
        im = 2;
        it = 2;
        [memb ind_test] = ismember(trialInd{1,it,1,isf},find(centroid_dist_sf{isf}<2));
        [memb ind_mask] = ismember(trialInd{im,1,1,isf},find(centroid_dist_sf{isf}<2));
        test_avg = mean(resp_cell{1,it,1,isf}(:,find(ind_test)),2);
        test_avg_rect = test_avg;
        test_avg_rect(find(test_avg<0)) = 0;
        mask_avg = mean(resp_cell{im,1,1,isf}(:,find(ind_mask)),2);
        mask_avg_rect = mask_avg;
        mask_avg_rect(find(mask_avg<0)) = 0;
        resp_avg = nan(nCells,nMaskPhas);
        resp_all = [];
        stim_all = [];
        for ip = 1:nMaskPhas
            [memb ind] = ismember(trialInd{im,it,ip,isf},find(centroid_dist_sf{isf}<2));
            resp_all = [resp_all resp_cell{im,it,ip,isf}(:,find(ind))];
            stim_all = [stim_all ip.*ones(size(resp_cell{im,it,ip,isf}(1,find(ind))))];
            resp_avg(:,ip) = mean(resp_cell{im,it,ip,isf}(:,find(ind)),2);
        end

        resp_all_rect = resp_all;
        resp_all_rect(find(resp_all<0)) = 0;
        SI_all = (resp_all_rect-(test_avg_rect+mask_avg_rect))./(resp_all_rect+(test_avg_rect+mask_avg_rect));
        resp_avg_rect = resp_avg;
        resp_avg_rect(find(resp_avg<0)) = 0;
        SI_avg = (resp_avg_rect-(test_avg_rect+mask_avg_rect))./(resp_avg_rect+(test_avg_rect+mask_avg_rect));
        [eye_n(isf,:) edges bin] = histcounts(stim_all,[1:5]);
        fprintf([num2str(SFs(isf)) ' SF Eye-n: ' num2str(eye_n(isf,:)) '\n'])
        if sum(squeeze(eye_n(isf,:))<4)==0
            figure;
            start = 1;
            n = 1;
            for iCell = 1:nCells
                if start>25
                    suptitle([mouse ' ' date '- Mask ' num2str(im) ' Test ' num2str(it) ' SF ' num2str(SFs(isf))])
                    print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_phaseFits_SI' num2str(n) '_M' num2str(im) 'T' num2str(it) 'SF' num2str(isf) '_LT2degEyeMvmt.pdf']), '-dpdf','-fillpage')
                    figure;
                    start = 1;
                    n = n+1;
                end
                p_anova(iCell,isf) = anova1(resp_all(iCell,:), stim_all,'off');
                if max(SI_avg(iCell,:),[],2)>min(SI_avg(iCell,:),[],2)
                    [b_hat(iCell,isf), amp_hat(iCell,isf), per_hat(iCell,isf),pha_hat(iCell,isf),sse(iCell,isf),R_square(iCell,isf)] = sinefit(deg2rad(maskPhas(stim_all)),SI_all(iCell,:));
                    subplot(5,5,start)
                    scatter(maskPhas(stim_all),SI_all(iCell,:));
                    hold on
                    scatter(maskPhas,SI_avg(iCell,:))
                    yfit(iCell,:,isf) = b_hat(iCell,isf)+amp_hat(iCell,isf).*(sin(2*pi*deg2rad(phase_range)./per_hat(iCell,isf) + 2.*pi/pha_hat(iCell,isf)));
                    plot(phase_range, yfit(iCell,:,isf));
                    title(['Rsq = ' num2str(chop(R_square(iCell,isf),2)) '; p = ' num2str(chop(p_anova(iCell,isf),2))])
                else
                    b_hat(iCell,isf) = max(SI_avg(iCell,:),[],2);
                    amp_hat(iCell,isf) = 0;
                    per_hat(iCell,isf) = NaN;
                    pha_hat(iCell,isf) = NaN;
                    sse(iCell,isf) = 0;
                    R_square(iCell,isf)= 0;
                    yfit(iCell,:,isf) = nan(length(phase_range),1);
                end
                start = start+1;
            end
            suptitle([mouse ' ' date '- Mask ' num2str(im) ' Test ' num2str(it) ' SF ' num2str(SFs(isf))])
            print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_phaseFits_SI' num2str(n) '_M' num2str(im) 'T' num2str(it) 'SF' num2str(isf) '_LT2degEyeMvmt.pdf']), '-dpdf','-fillpage')
            
            rng(seed);
            stim_all_shuf = stim_all(randperm(length(stim_all)));
            figure;
            start = 1;
            n = 1;
            for iCell = 1:nCells
                if start>25
                    suptitle([mouse ' ' date '- Shuffled- SF ' num2str(SFs(isf)) ])
                    print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_phaseFits_SI_SF' num2str(isf) '_shuffled' num2str(n) '_LT2degEyeMvmt.pdf']), '-dpdf','-fillpage')
                    figure;
                    start = 1;
                    n = n+1;
                end
                p_anova_shuf(iCell,isf) = anova1(resp_all(iCell,:), stim_all_shuf,'off');
                if max(SI_avg(iCell,:),[],2)>min(SI_avg(iCell,:),[],2)
                    [b_hat_shuf(iCell,isf), amp_hat_shuf(iCell,isf), per_hat_shuf(iCell,isf),pha_hat_shuf(iCell,isf),sse_shuf(iCell,isf),R_square_shuf(iCell,isf)] = sinefit(deg2rad(maskPhas(stim_all_shuf)),SI_all(iCell,:));
                    subplot(5,5,start)
                    scatter(maskPhas(stim_all_shuf),SI_all(iCell,:));
                    hold on
                    scatter(maskPhas,SI_avg(iCell,:))
                    yfit_shuf(iCell,:,isf) = b_hat_shuf(iCell,isf)+amp_hat_shuf(iCell,isf).*(sin(2*pi*deg2rad(phase_range)./per_hat_shuf(iCell,isf) + 2.*pi/pha_hat_shuf(iCell,isf)));
                    plot(phase_range, yfit_shuf(iCell,:,isf));
                    title(['Rsq = ' num2str(chop(R_square_shuf(iCell,isf),2)) '; p = ' num2str(chop(p_anova_shuf(iCell,isf),2))])
                else
                    b_hat_shuf(iCell,isf) = max(SI_avg(iCell,:),[],2);
                    amp_hat_shuf(iCell,isf) = 0;
                    per_hat_shuf(iCell,isf) = NaN;
                    pha_hat_shuf(iCell,isf) = NaN;
                    sse_shuf(iCell,isf) = 0;
                    R_square_shuf(iCell,isf)= 0;
                    yfit_shuf(iCell,:,isf) = nan(length(phase_range),1);
                end
                start = start+1;
            end
            suptitle([mouse ' ' date '- Shuffled- SF ' num2str(SFs(isf)) ])
            print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_phaseFits_SI_SF' num2str(isf) '_shuffled' num2str(n) '_LT2degEyeMvmt.pdf']), '-dpdf','-fillpage')
        end
    end
    save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_phaseFits.mat']),'seed', 'yfit', 'b_hat', 'amp_hat', 'per_hat', 'pha_hat', 'sse', 'R_square',  'p_anova', 'yfit_shuf', 'b_hat_shuf', 'amp_hat_shuf', 'per_hat_shuf', 'pha_hat_shuf', 'sse_shuf', 'R_square_shuf', 'p_anova_shuf','trial_n', 'trialInd', 'eye_n')



%%

    if nSF >1

        sum(p_anova(iCell,:)<0.05,1)

        plaid_resp = zeros(nCells,nSF);
        mask_resp = zeros(nCells,nSF);
        test_resp = zeros(nCells,nSF);
        plaidSI = zeros(nCells,nSF);
        testPI = zeros(nCells,nSF);

        for isf = 1:nSF
            plaid_resp(:,isf) = mean(resp_cell{end,end,1,isf},2);
            mask_resp(:,isf) =  mean(resp_cell{end,1,1,isf},2);
            test_resp(:,isf) =  mean(resp_cell{1,end,1,isf},2);
            plaid_resp(find(plaid_resp(:,isf)<0),isf) = 0;
            mask_resp(find(mask_resp(:,isf)<0),isf) = 0;
            test_resp(find(test_resp(:,isf)<0),isf) = 0;
            plaidSI(:,isf) = (plaid_resp(:,isf)-(mask_resp(:,isf)+test_resp(:,isf))) ./ (plaid_resp(:,isf) + mask_resp(:,isf) + test_resp(:,isf));
            testPI(:,isf) = abs((test_resp(:,isf)-mask_resp(:,isf)) ./ (mask_resp(:,isf)+test_resp(:,isf)));
        end

        figure;
        for isf = 1:nSF
            subplot(1,2,1)
            cdfplot(plaidSI(resp_ind_sf{isf},isf))
            xlim([-1 1])
            hold on
            subplot(1,2,2)
            cdfplot(testPI(resp_ind_sf{isf},isf))
            xlim([0 1])
            hold on
        end
        subplot(1,2,1)
        xlabel('Suppression index')
        subplot(1,2,2)
        legend(num2str(SFs'),'location', 'northwest')
        xlabel('Preference index')

        data_dfof_stim = zeros(nCells,3,nSF);
        for isf = 1:nSF
            data_dfof_stim(:,1,isf) = mean(resp_cell{end,1,1,isf},2);
            data_dfof_stim(:,2,isf) = mean(resp_cell{1,end,1,isf},2);
            data_dfof_stim(:,3,isf) = mean(resp_cell{end,end,1,isf},2);
        end

        [max_dir_val max_dir_ind] = max(squeeze(mean(data_dfof_stim,3)),[],2);
        max_sf_val = zeros(1,nCells);
        max_sf_ind = zeros(1,nCells);
        fit_out = cell(1,nCells);
        g_fit = cell(1,nCells);
        prefSF = zeros(1,nCells);
        % figure; movegui('center')
        % start = 1;
        for iCell = 1:nCells
        %     if start >49
        %         figure;
        %         start = 1;
        %     end
            [max_sf_val(1,iCell) max_sf_ind(1,iCell)] = max(data_dfof_stim(iCell,max_dir_ind(iCell),:),[],3);
            [fit_out{iCell} g_fit{iCell}] =fit(log2(SFs)',squeeze(data_dfof_stim(iCell,max_dir_ind(iCell),:)),'gauss1','Lower',[0 log2(SFs(1)) 0],'Upper',[Inf log2(SFs(end)) Inf]);
        %     subplot(7,7,start)
        %     plot(log2(SFs),squeeze(data_dfof_stim(iCell,max_dir_ind(iCell),:)),'o')
        %     hold on
        %     plot(fit_out{iCell})
            prefSF(1,iCell)= fit_out{iCell}.b1;
            RsqSF(1,iCell)= g_fit{iCell}.rsquare;
        %     if find(h_all == iCell)
        %         title('Sig')
        %     end
        %     start = start+1;
        %     legend('off')
            if rem(iCell, 10) == 0
                fprintf([num2str(iCell) '\n'])
            end
        end
        save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_SFfits.mat']),'fit_out','g_fit','prefSF','RsqSF','max_sf_ind','max_dir_ind','plaidSI','testPI')

        figure; movegui('center');
        RsqSF(find(RsqSF<0)) = 0;
        subplot(2,1,1); histogram(RsqSF(resp_ind)); vline(0.8)
        ind = intersect(find(RsqSF>0.8),resp_ind); xlabel('Rsq')
        subplot(2,1,2); histogram(prefSF(ind),nSF)
        set(gca, 'XTick', log2(SFs), 'XTickLabels', SFs)
        vline(mean(prefSF(ind),2))
        xlabel('SF (cpd)')

        figure;
        amp_diff = amp_hat-amp_hat_shuf;
        legstr = [];
        for isf = 1:nSF
            if length(resp_ind_sf{isf})>1
                if sum(~isnan(amp_diff(resp_ind_sf{isf},isf)),1)>1
                cdfplot(amp_diff(resp_ind_sf{isf},isf))
                hold on
                legstr = [legstr; [num2str(SFs(isf)) '- n = ' num2str(length(resp_ind_sf{isf}))]];
                end
            end
        end
        xlabel('Sine Amp: Data-Shuf')
        legend(legstr,'location', 'southeast')
        title('')
        print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_diffAmpByStimSF_cdfs.pdf']),'-dpdf','-bestfit');


        figure;
        movegui('center')
        for isf = 1:nSF
            subplot(3,2,isf)
            for iCell = 1:nCells
                if find(resp_ind_sf{isf} == iCell)
                    if p_anova(iCell,isf)<0.05
                        scatter(2.^prefSF(iCell), amp_diff(iCell,isf),'or')
                    else
                        scatter(2.^prefSF(iCell), amp_diff(iCell,isf),'ok')
                    end
                    hold on
                end
            end
            ylim([-0.2 0.5])
            xlim([0.02 0.32])
            ylabel('Sine Amp: Data-Shuf')
            xlabel('Preferred SF')
            title(['Stim SF: ' num2str(SFs(isf))])
        end
        subplot(3,2,6)
        for iCell = 1:nCells
            if find(resp_ind == iCell)
                if p_anova(iCell,max_sf_ind(1,iCell))<0.05
                    scatter(2.^prefSF(iCell), amp_diff(iCell,max_sf_ind(1,iCell)),'or')
                else
                    scatter(2.^prefSF(iCell), amp_diff(iCell,max_sf_ind(1,iCell)),'ok')
                end
                hold on
            end
        end
        ylabel('Sine Amp: Data-Shuf')
        xlabel('Preferred SF')
        title('At Peak SF')
        ylim([-0.2 0.5])
        xlim([0.02 0.32])
        print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_diffAmpByPrefSF_scatter.pdf']),'-dpdf','-fillpage');
    end
    close all
end