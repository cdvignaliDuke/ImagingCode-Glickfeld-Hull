mouse_mat = strvcat('i671','i698','i699');
date_mat{1} = strvcat('170127','170131');
date_mat{2} = strvcat('170923','170930','171001','171002','171003','171004');
date_mat{3} = strvcat('170923','170930','171013');
run_mat{1} = {'002-004','003-006'};
run_mat{2} = {'001-003','002-003', '002', '002-003','001-002','001'};
run_mat{3} = {'001', '001','001-002'};
dir_mat{1} = {'005','007'};
dir_mat{2} = {nan,nan,nan,nan,nan,nan};
dir_mat{3} = {nan,nan,nan};

%%

% fRate = cell(size(mouse_mat,1),1);
% lRate = cell(size(mouse_mat,1),1);
% ntrials = cell(size(mouse_mat,1),1);
% base_trials = cell(size(mouse_mat,1),1);
% targ_trials = cell(size(mouse_mat,1),1);
% base_cells = cell(size(mouse_mat,1),1);
% targ_cells = cell(size(mouse_mat,1),1);
% used_expt = [];
% base_resp_avg_all = [];
% base_resp_off_avg_all = [];
% base_resp_off_cyc_avg_all = [];
% targ_resp_off_avg_all = [];
% targ_resp_off_delta_avg_all = [];
% targ_resp_off_delta_outcome_avg_all = [];
% targ_resp_off_delta_outcome_M_avg_all = [];
t = 1;
for imouse = 1:size(mouse_mat,1)
    mouse = mouse_mat(imouse,:);
    for iexp = 1:size(date_mat{imouse},1)
        date = date_mat{imouse}(iexp,:);
        run = run_mat{imouse}{iexp};
        load(fullfile('Z:\home\lindsey\Analysis\2P',[date '_' mouse],[date '_' mouse '_runs-' run],[date '_' mouse '_runs-' run '_input.mat']));
        load(fullfile('Z:\home\lindsey\Analysis\2P',[date '_' mouse],[date '_' mouse '_runs-' run],[date '_' mouse '_runs-' run '_dfofData.mat']));
        load(fullfile('Z:\home\lindsey\Analysis\2P',[date '_' mouse],[date '_' mouse '_runs-' run],[date '_' mouse '_runs-' run '_stimData.mat']));

        tGratingDir = celleqel2mat_padded(input.tGratingDirectionDeg);
        nTrials = size(tGratingDir,2);
        dirs = unique(tGratingDir);
        ind_easy = find(tGratingDir == dirs(end));
        reactTimes = celleqel2mat_padded(input.reactTimesMs);
        SIx = zeros(1,nTrials);
        MIx = zeros(1,nTrials);
        holdTimes = celleqel2mat_padded(input.holdTimesMs);
        tCyc = celleqel2mat_padded(input.tCyclesOn);
        nCyc = celleqel2mat_padded(input.nCyclesOn);
        
%         ntrials{imouse}(iexp,:) = sum(FIx+MIx+SIx);
%         fRate{imouse}(iexp,:) = sum(FIx)./(sum(FIx+MIx+SIx));
%         lRate{imouse}(iexp,:) = sum(SIx(ind_easy))./(sum(SIx(ind_easy))+sum(MIx(ind_easy)));

        
        base_win = 21:23;
        resp_win = 29:31;
        react_win = [200 550];
        sz = size(data_dfof);
        base_dfof = nan(sz);
        nOn = unique(cell2mat(input.nFramesOn));
        FAIx = zeros(1,nTrials);
        FA_resp = nan(sz(2),sz(4));
        preFA_resp = nan(sz(2),sz(4));
        FAInt = nan(1,nTrials);
        preFAInt = nan(1,nTrials);
        FACyc = nan(1,nTrials);
        for itrial = 1:nTrials
            if holdTimes(itrial)>700
                if tCyc(itrial) == nCyc(itrial)
                    totTime = (sum(tFramesOff(itrial,1:nCyc(itrial)),2)+(nCyc(itrial).*nOn)).*(1000/frameRateHz);
                    if holdTimes(itrial)-totTime>react_win(1) & holdTimes(itrial)-totTime<react_win(2)
                        SIx(itrial)=1;
                    elseif holdTimes(itrial)-totTime>react_win(2)
                        MIx(itrial) =1;
                    else
                        totTime = (sum(tFramesOff(itrial,1:nCyc(itrial)-1),2)+((nCyc(itrial)-1).*nOn)).*(1000/frameRateHz);
                        if holdTimes(itrial)-totTime>react_win(1) & holdTimes(itrial)-totTime<react_win(2)
                            if nCyc(itrial)-3 > 0
                                FAIx(itrial)=1;
                                FACyc(itrial) = nCyc(itrial)-1;
                                FA_resp(:,itrial) = squeeze(mean(data_dfof(resp_win,:,nCyc(itrial)-1,itrial),1)-mean(data_dfof(base_win,:,nCyc(itrial)-1,itrial),1));
                                preFA_resp(:,itrial) = squeeze(mean(data_dfof(resp_win,:,nCyc(itrial)-2,itrial),1)-mean(data_dfof(base_win,:,nCyc(itrial)-2,itrial),1));
                                FAInt(itrial) = tFramesOff(itrial,nCyc(itrial)-2);
                                preFAInt(itrial) = tFramesOff(itrial,nCyc(itrial)-3);
                            end
                        end
                    end
                else
                    totTime = (sum(tFramesOff(itrial,1:tCyc(itrial)),2)+(tCyc(itrial).*nOn)).*(1000/frameRateHz);
                    if holdTimes(itrial)-totTime>react_win(1) & holdTimes-totTime(itrial)<react_win(2)
                        if tCyc(itrial)-2 > 0
                            FAIx(itrial)=1;
                            FACyc(itrial) = tCyc(itrial);
                            FA_resp(:,itrial) = squeeze(mean(data_dfof(resp_win,:,tCyc(itrial),itrial),1)-mean(data_dfof(base_win,:,tCyc(itrial),itrial),1));
                            preFA_resp(:,itrial) = squeeze(mean(data_dfof(resp_win,:,tCyc(itrial)-1,itrial),1)-mean(data_dfof(base_win,:,tCyc(itrial)-1,itrial),1));
                            FAInt(itrial) = tFramesOff(itrial,tCyc(itrial)-1);
                            preFAInt(itrial) = tFramesOff(itrial,tCyc(itrial)-2);
                        end
                    else
                        totTime = (sum(tFramesOff(itrial,1:tCyc(itrial)-1),2)+((tCyc(itrial)-1).*nOn)).*(1000/frameRateHz);
                        if holdTimes(itrial)-totTime>react_win(1) & holdTimes(itrial)-totTime<react_win(2)
                            if tCyc(itrial)-3 > 0
                                FAIx(itrial)=1;
                                FACyc(itrial) = tCyc(itrial)-1;
                                FA_resp(:,itrial) = squeeze(mean(data_dfof(resp_win,:,tCyc(itrial)-1,itrial),1)-mean(data_dfof(base_win,:,tCyc(itrial)-1,itrial),1));
                                preFA_resp(:,itrial) = squeeze(mean(data_dfof(resp_win,:,tCyc(itrial)-2,itrial),1)-mean(data_dfof(base_win,:,tCyc(itrial)-2,itrial),1));
                                FAInt(itrial) = tFramesOff(itrial,tCyc(itrial)-2);
                                preFAInt(itrial) = tFramesOff(itrial,tCyc(itrial)-3);
                            end
                        end
                    end
                end
            end
        end
        
        MIx_bx = strcmp(input.trialOutcomeCell,'ignore');
        MIx(find(MIx_bx))=1;
        SIx(find(MIx_bx))=0;
        
        [h_base, p_base] = ttest(squeeze(mean(data_dfof(resp_win,:,1,:),1)),squeeze(mean(data_dfof(base_win,:,1,:),1)),'dim',2,'alpha', 0.01);
        good_ind_base = find(h_base);
        data_dfof_targ = nan(sz(1),sz(2),sz(4));
        targCyc = double(nCyc)+1;
        FIx = zeros(1,nTrials);
        FIx(setdiff(1:nTrials,[find(MIx) find(SIx)])) = 1;
        targCyc(find(FIx)) = NaN;
        expt(t).targInt = nan(1,nTrials);
        expt(t).lastBaseInt = nan(1,nTrials);
        for itrial = 1:nTrials
            if ~isnan(targCyc(itrial))
                data_dfof_targ(:,:,itrial) = data_dfof(:,:,targCyc(itrial),itrial);
                expt(t).targInt(1,itrial) = tFramesOff(itrial,targCyc(itrial)-1);
                expt(t).lastBaseInt(1,itrial) = tFramesOff(itrial,targCyc(itrial)-2);
            end
        end
        targCyc(find(FAIx)) = FACyc(find(FAIx));
        expt(t).targInt(find(FAIx)) = FAInt(find(FAIx));
        expt(t).lastBaseInt(find(FAIx)) = preFAInt(find(FAIx));
        targ_resp = squeeze(bsxfun(@minus,data_dfof_targ,mean(data_dfof_targ(base_win,:,:),1)));
        deltas = unique(tGratingDir);
        nDelta = length(deltas);
        h_targ = zeros(nCells,nDelta);
        for i = 1:nDelta
            ind = find(tGratingDir==deltas(i));
            [h_targ(:,i), p_targ] = ttest(squeeze(mean(data_dfof_targ(resp_win,:,ind),1)),squeeze(mean(data_dfof_targ(base_win,:,ind),1)),'dim',2,'alpha', 0.05./(nDelta-1));
        end
        good_ind_targ = find(sum(h_targ,2));
        good_ind = unique([good_ind_base; good_ind_targ]);

        expt(t).firstBaseResp = squeeze(mean(data_dfof(resp_win,:,1,:),1)- mean(data_dfof(base_win,:,1,:),1));
        expt(t).targetResp = squeeze(mean(data_dfof_targ(resp_win,:,:),1)-mean(data_dfof_targ(base_win,:,:),1));
        expt(t).targetResp(:,find(FAIx)) = FA_resp(:,find(FAIx));
        expt(t).lastBaseResp = nan(size(expt(t).firstBaseResp));
        expt(t).trialOutcome = cell(1,nTrials);
        for itrial = 1:nTrials
            if SIx(itrial)
                expt(t).lastBaseResp(:,itrial) = squeeze(mean(data_dfof(resp_win,:,targCyc(itrial)-1,itrial),1)-mean(data_dfof(base_win,:,targCyc(itrial)-1,itrial),1));
                expt(t).trialOutcome{itrial} = 's';
            elseif MIx(itrial)
                expt(t).lastBaseResp(:,itrial) = squeeze(mean(data_dfof(resp_win,:,targCyc(itrial)-1,itrial),1)-mean(data_dfof(base_win,:,targCyc(itrial)-1,itrial),1));
                expt(t).trialOutcome{itrial} = 'm';
            elseif FAIx(itrial)
                expt(t).lastBaseResp(:,itrial) = preFA_resp(:,itrial);
                expt(t).trialOutcome{itrial} = 'f';
            else
                expt(t).trialOutcome{itrial} = nan;
            end
        end
        expt(t).targetOrientation = targetDelta;
        expt(t).targetOrientation(find(FIx)) = nan;
        expt(t).targetOrientation(find(FAIx)) = 0;
        expt(t).signifResponsiveCells = zeros(1,sz(2));
        expt(t).signifResponsiveCells(1,good_ind) = 1;
        expt(t).lastBaseCycN= targCyc-1;
        expt(t).mouse = mouse;
        expt(t).date = date;
        expt(t).info = 'cells x trials';
        
        ind = find(expt(t).targetOrientation==90);
        resp_90 = nanmean(expt(t).targetResp(:,ind),2);
        resp_0 = nanmean(expt(t).firstBaseResp,2);
        expt(t).baseTargRatio = (resp_0-resp_90)./(resp_0+resp_90);
        
        if ~isnan(dir_mat{imouse}{iexp})
            dir_run = dir_mat{imouse}{iexp};
            load(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_runs-' dir_run], [date '_' mouse '_runs-' dir_run '_oriFit.mat']))
            expt(t).oriTuning = max_ori(:,1);
            expt(t).oriTuningTheta90 = theta_90;
        else
            expt(t).oriTuning = [];
            expt(t).oriTuningTheta90 = [];
        end
        
        t = 1+t;
    end
end
save(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P\Adaptation', 'randIntBehavData.mat'),'expt')
%% figures 
        
        %all cells
        figure;
        longHoldTimes = find(celleqel2mat_padded(input.holdTimesMs)>500);
        for iCell = 1:nCells 
            subplot(3,2,1)
            plot(tt(1:40),nanmean(mean(base_resp(1:40,iCell,1,longHoldTimes),2),4));
            hold on
            hline(0)
            vline(0)
            xlim([tt(1) tt(40)])
            subplot(3,2,2)
            plot(tt(1:40),nanmean(mean(targ_resp(1:40,iCell,:),2),3));
            hold on
            hline(0)
            vline(0)
            xlim([tt(1) tt(40)])
        end
        subplot(3,2,1)
        title('First Base Resp')
        subplot(3,2,2)
        title('Targ Resp')
        subplot(3,2,3)
        leg_str = [];
        for icyc = 1:maxCyc
            if sum(~isnan(base_resp(1,1,icyc,longHoldTimes)),4) > 20
                plot(tt(1:40),nanmean(mean(base_resp(1:40,:,icyc,longHoldTimes),2),4));
                hold on
                leg_str = [leg_str num2str(icyc)];
            end
        end
        hline(0)
        vline(0)
        xlim([tt(1) tt(40)])
        title('Avg Base')
        legend(leg_str','Location','northwest')
        subplot(3,2,4)
        plot(tt(1:40),nanmean(mean(targ_resp(1:40,:,:),2),3));
        hline(0)
        vline(0)
        xlim([tt(1) tt(40)])
        title('Avg Targ')
        subplot(3,2,5)
        for ioff = 1:noff
            base_resp_off = [];
            for icyc = 2:maxCyc
                ind = find(tFramesOff(:,icyc-1) == offs(ioff));
                base_resp_off = cat(4,base_resp_off,base_resp(:,:,icyc,ind));
            end
            plot(tt(1:40),nanmean(mean(base_resp_off(1:40,:,:,:),2),4))
            hold on
        end
        hline(0)
        vline(0)
        xlim([tt(1) tt(40)])
        legend(num2str(chop(offs.*(1000/frameRateHz),3)),'Location','northwest')
        title('Avg Base ISI')
        subplot(3,2,6)
        for ioff = 1:noff
            ind = find(targInt == offs(ioff));
            plot(tt(1:40),nanmean(mean(targ_resp(1:40,:,ind),2),3));
            hold on
        end
        hline(0)
        vline(0)
        xlim([tt(1) tt(40)])
        title('Avg Targ ISI')
        suptitle([mouse ' ' date '- ' num2str(nTrials{imouse}(iexp,:)) ' trials - FA rate: ' num2str(chop(fRate{imouse}(iexp,:),2)) '; Lapse rate: ' num2str(chop(lRate{imouse}(iexp,:),2))])
        print(fullfile('Z:\home\lindsey\Analysis\2P\Adaptation\Behavior',[date '_' mouse '_runs-' run '_summary.pdf']),'-dpdf', '-fillpage')
    
        %responsive cells
        figure;
        holdTimes = celleqel2mat_padded(input.holdTimesMs);
        longHoldTimes = find(celleqel2mat_padded(input.holdTimesMs)>500);
        for iCell = 1:length(good_ind_base)
            iC = good_ind_base(iCell);
            subplot(3,2,1)
            plot(tt(1:40),nanmean(mean(base_resp(1:40,iC,1,longHoldTimes),2),4));
            hold on
            hline(0)
            vline(0)
            xlim([tt(1) tt(40)])
        end
        for iCell = 1:length(good_ind_targ)
            subplot(3,2,2)
            iC = good_ind_targ(iCell);
            plot(tt(1:40),nanmean(mean(targ_resp(1:40,iC,:),2),3));
            hold on
            hline(0)
            vline(0)
            xlim([tt(1) tt(40)])
        end
        subplot(3,2,1)
        title('First Base Resp')
        subplot(3,2,2)
        title('Targ Resp')
        subplot(3,2,3)
        leg_str = [];
        base_resp_avg = nan(40,length(good_ind_base),10);
        for icyc = 1:maxCyc
            if sum(~isnan(base_resp(1,1,icyc,longHoldTimes)),4) > 20
                plot(tt(1:40),nanmean(mean(base_resp(1:40,good_ind_base,icyc,longHoldTimes),2),4));
                hold on
                leg_str = [leg_str num2str(icyc)];
                base_resp_avg(:,:,icyc) = nanmean(base_resp(1:40,good_ind_base,icyc,longHoldTimes),4);
            end
        end
        hline(0)
        vline(0)
        xlim([tt(1) tt(40)])
        title('Avg Base')
        legend(leg_str','Location','northwest')
        subplot(3,2,4)
        plot(tt(1:40),nanmean(mean(targ_resp(1:40,good_ind_targ,:),2),3));
        hline(0)
        vline(0)
        xlim([tt(1) tt(40)])
        title('Avg Targ')
        subplot(3,2,5)
        base_cells{imouse}(1,iexp) = length(good_ind_base);
        base_resp_off_avg = nan(40, length(good_ind_base), noff);
        for ioff = 1:noff
            base_resp_off = [];
            for icyc = 3:maxCyc
                ind = find(tFramesOff(:,icyc-1) == offs(ioff));
                base_resp_off = cat(4,base_resp_off,base_resp(:,:,icyc,ind));
            end
            if sum(~isnan(base_resp_off(1,1,1,:)))>20
                base_resp_off_avg(:,:,ioff) = nanmean(base_resp_off(1:40,good_ind_base,:,:),4);
                plot(tt(1:40),nanmean(mean(base_resp_off(1:40,good_ind_base,:,:),2),4))
                hold on
                base_trials{imouse}(iexp,ioff) = sum(~isnan(base_resp_off(1,1,1,:)));
            end
        end
        hline(0)
        vline(0)
        xlim([tt(1) tt(40)])
        legend(num2str(chop(offs.*(1000/frameRateHz),3)),'Location','northwest')
        title('Avg Base ISI')
        base_resp_off_cyc_avg = nan(40, length(good_ind_base), noff,10);
        for ioff = 1:noff
            base_resp_off = [];
            for icyc = 2:maxCyc
                ind = find(tFramesOff(:,icyc-1) == offs(ioff));
                if length(ind)>20
                     base_resp_off_cyc_avg(:,:,ioff,icyc) = nanmean(base_resp(1:40,good_ind_base,icyc,ind),4);
                end
            end
        end
        subplot(3,2,6)
        deltas = [23 90];
        nDelta = length(deltas);
        targ_cells{imouse}(1,iexp) = length(good_ind_targ);
        targ_resp_off_avg = nan(40,length(good_ind_targ),noff);
        targ_resp_off_delta_avg = nan(40,length(good_ind_targ),noff,nDelta);
        targ_resp_off_delta_outcome_avg = nan(40,length(good_ind_targ),noff,nDelta,2);
        targ_resp_off_delta_outcome_M_avg = nan(40,length(good_ind_targ_miss),noff,nDelta,2);
        for ioff = 1:noff
            ind = find(targInt == offs(ioff));
            plot(tt(1:40),nanmean(mean(targ_resp(1:40,good_ind_targ,ind),2),3));
            hold on
            %if length(ind) > 10
                targ_resp_off_avg(:,:,ioff) = nanmean(targ_resp(1:40,good_ind_targ,ind),3);
                for idelta = 1:nDelta
                    ind2 = intersect(ind, find(tGratingDir == deltas(idelta)));
                    %if length(ind2) > 5
                        targ_resp_off_delta_avg(:,:,ioff,idelta) = nanmean(targ_resp(1:40,good_ind_targ,ind2),3);
                        ind_SIx = intersect(ind2, find(SIx));
                        ind_MIx = intersect(ind2, find(MIx));
                        %if length(ind_SIx) > 5 
                            targ_resp_off_delta_outcome_avg(:,:,ioff,idelta,1) = nanmean(targ_resp(1:40,good_ind_targ,ind_SIx),3);
                            targ_resp_off_delta_outcome_M_avg(:,:,ioff,idelta,1) = nanmean(targ_resp(1:40,good_ind_targ_miss,ind_SIx),3);
                        %end
                        %if length(ind_MIx) > 5 
                            targ_resp_off_delta_outcome_avg(:,:,ioff,idelta,2) = nanmean(targ_resp(1:40,good_ind_targ,ind_MIx),3);
                            targ_resp_off_delta_outcome_M_avg(:,:,ioff,idelta,2) = nanmean(targ_resp(1:40,good_ind_targ_miss,ind_MIx),3);
                        %end
                    %end
                end
            %end
        end
        for ioff = 1:noff
            ind = find(targInt == offs(ioff));
            for idelta = 1:nDelta
            	ind2 = intersect(ind, find(tGratingDir == deltas(idelta)));
                ind_SIx = intersect(ind2, find(SIx));
                ind_MIx = intersect(ind2, find(MIx));
                targ_trials{imouse}(iexp,ioff,idelta,1) = length(ind_SIx);
                targ_trials{imouse}(iexp,ioff,idelta,2) = length(ind_MIx);
            end
        end
        hline(0)
        vline(0)
        xlim([tt(1) tt(40)])
        title('Avg Targ ISI')
        suptitle([mouse ' ' date '- ' num2str(length(good_ind_base)) ' base resp cells; ' num2str(length(good_ind_targ)) ' targ resp cells'])
        print(fullfile('Z:\home\lindsey\Analysis\2P\Adaptation\Behavior',[date '_' mouse '_runs-' run '_goodCellSummary.pdf']),'-dpdf', '-fillpage')
        
        save(fullfile('Z:\home\lindsey\Analysis\2P',[date '_' mouse],[date '_' mouse '_runs-' run],[date '_' mouse '_runs-' run '_baseSummary.mat']),'base_resp_avg','base_resp_off_avg','base_resp_off_cyc_avg');
        save(fullfile('Z:\home\lindsey\Analysis\2P',[date '_' mouse],[date '_' mouse '_runs-' run],[date '_' mouse '_runs-' run '_targSummary.mat']),'targ_resp_off_avg','targ_resp_off_delta_avg', 'targ_resp_off_delta_outcome_avg');
        if lRate{imouse}(iexp,:) > 0.85 & fRate{imouse}(iexp,:) < 0.4
            used_expt = [used_expt; [mouse '_' date]];
            base_resp_avg_all = cat(2,base_resp_avg_all, base_resp_avg);
            base_resp_off_avg_all = cat(2,base_resp_off_avg_all, base_resp_off_avg);
            base_resp_off_cyc_avg_all = cat(2,base_resp_off_cyc_avg_all, base_resp_off_cyc_avg);
            targ_resp_off_avg_all = cat(2,targ_resp_off_avg_all, targ_resp_off_avg);
            targ_resp_off_delta_avg_all = cat(2,targ_resp_off_delta_avg_all, targ_resp_off_delta_avg);
            targ_resp_off_delta_outcome_avg_all = cat(2,targ_resp_off_delta_outcome_avg_all, targ_resp_off_delta_outcome_avg);
            targ_resp_off_delta_outcome_M_avg_all = cat(2,targ_resp_off_delta_outcome_M_avg_all,targ_resp_off_delta_outcome_M_avg);
        end
    end
end
save(fullfile('Z:\home\lindsey\Analysis\2P\Adaptation\Behavior\behavSummary.mat'),'mouse_mat','date_mat', 'run_mat','fRate','lRate','nTrials','base_trials','targ_trials','used_expt','base_resp_avg_all','base_resp_off_avg_all','base_resp_off_cyc_avg_all','targ_resp_off_avg_all','targ_resp_off_delta_avg_all','targ_resp_off_delta_outcome_avg_all');

figure; 
plot(tt(1:40),squeeze(nanmean(base_resp_avg_all(:,:,1),2)))
hold on
plot(tt(1:40),squeeze(nanmean(nanmean(targ_resp_off_avg_all(:,:,:),2),3)))
title(['Resp win: ' num2str(tt(resp_win))])
print(fullfile('Z:\home\lindsey\Analysis\2P\Adaptation\Behavior\allExptBaseAndTargResp_longwin.pdf'),'-dpdf', '-fillpage')


figure; 
subplot(2,2,1)
plot(tt(1:40),squeeze(nanmean(bsxfun(@minus,base_resp_avg_all,mean(base_resp_avg_all(base_win,:,:),1)),2)))
vline((resp_win-19)*(1000/frameRateHz))
xlabel('Time from stim (ms)')
ylabel('dF/F')
title('Avg by cycle #')
subplot(2,2,2)
plot(tt(1:40),squeeze(nanmean(bsxfun(@minus,base_resp_off_avg_all,mean(base_resp_off_avg_all(base_win,:,:),1)),2)))
vline((resp_win-19)*(1000/frameRateHz))
xlabel('Time from stim (ms)')
ylabel('dF/F')
title('Avg by ISI')
legend(num2str(chop(offs.*(1000/frameRateHz),2)),'Location','southeast')
subplot(2,2,3)
n = squeeze(sum(~isnan(base_resp_avg_all(1,:,:)),2));
errorbar(1:size(base_resp_avg_all,3),squeeze(nanmean(mean(base_resp_avg_all(resp_win,:,:),1),2)),squeeze(nanstd(mean(base_resp_avg_all(resp_win,:,:),1),[],2))./sqrt(n),'-o')
ylim([0 inf])
xlabel('Cycle #')
ylabel('dF/F')
title('Avg by cycle #')
subplot(2,2,4)
n = squeeze(sum(~isnan(base_resp_off_avg_all(1,:,:)),2));
errorbar(offs.*(1000/frameRateHz),squeeze(nanmean(mean(base_resp_off_avg_all(resp_win,:,:),1),2)),squeeze(nanstd(mean(base_resp_off_avg_all(resp_win,:,:),1),[],2))./sqrt(n),'-o')
xlabel('ISI (ms)')
ylabel('dF/F')
title('Avg by ISI')
xlim([0 1000])
ylim([0 inf])
suptitle(['Random interval behavior summary- n = ' num2str(size(base_resp_avg_all,2)) ' cells'])
print(fullfile('Z:\home\lindsey\Analysis\2P\Adaptation\Behavior\allExptCycAndIsiSummary.pdf'),'-dpdf', '-fillpage')

ttAxisTick = find(ismember(floor(tt),-600:100:600));
ttAxisLabel = floor(tt(ttAxisTick));
figure;
subplot(3,2,1)
imagesc(permute(base_resp_avg_all(:,:,1),[2 1]))
colormap(brewermap([],'*RdBu'));
figXAxis([],'time (ms)',[],ttAxisTick,ttAxisLabel)
figYAxis([],'cell #',[])
clim([-0.1 0.1])
title('First Stim Response')
subplot(3,2,3)
imagesc(permute(nanmean(nanmean(targ_resp_off_delta_outcome_avg_all(:,:,:,:,1),3),4),[2 1]))
colormap(brewermap([],'*RdBu'));
figXAxis([],'time (ms)',[],ttAxisTick,ttAxisLabel)
figYAxis([],'cell #',[])
clim([-0.1 0.1])
title('Target Response- all resp- Hits')
subplot(3,2,5)
imagesc(permute(nanmean(nanmean(targ_resp_off_delta_outcome_avg_all(:,:,:,:,2),3),4),[2 1]))
colormap(brewermap([],'*RdBu'));
figXAxis([],'time (ms)',[],ttAxisTick,ttAxisLabel)
figYAxis([],'cell #',[])
clim([-0.1 0.1])
title('Target Response- all resp- Misses')
subplot(3,2,2)
imagesc(permute(base_resp_avg_all(:,:,1),[2 1]))
colormap(brewermap([],'*RdBu'));
figXAxis([],'time (ms)',[],ttAxisTick,ttAxisLabel)
figYAxis([],'cell #',[])
clim([-0.1 0.1])
title('First Stim Response')
subplot(3,2,4)
imagesc(permute(nanmean(nanmean(targ_resp_off_delta_outcome_M_avg_all(:,:,:,:,1),3),4),[2 1]))
colormap(brewermap([],'*RdBu'));
figXAxis([],'time (ms)',[],ttAxisTick,ttAxisLabel)
figYAxis([],'cell #',[])
clim([-0.1 0.1])
title('Target Response- miss resp- Hits')
subplot(3,2,6)
imagesc(permute(nanmean(nanmean(targ_resp_off_delta_outcome_M_avg_all(:,:,:,:,2),3),4),[2 1]))
colormap(brewermap([],'*RdBu'));
figXAxis([],'time (ms)',[],ttAxisTick,ttAxisLabel)
figYAxis([],'cell #',[])
clim([-0.1 0.1])
title('Target Response- miss resp- Misses')
print(fullfile('Z:\home\lindsey\Analysis\2P\Adaptation\Behavior\allExptBaseAndTargResp.pdf'),'-dpdf', '-fillpage')

figure; 
[n,n2] = subplotn(size(base_resp_avg_all,2));
for iCell = 1:size(base_resp_avg_all,2)
    subplot(n,n2,iCell)
    plot(tt(1:40),squeeze(bsxfun(@minus,base_resp_off_avg_all(:,iCell,:),mean(base_resp_off_avg_all(base_win,iCell,:),1))))
    title(num2str(iCell))
end
suptitle('All cells average response by ISI- blue- 250; red- 500; yellow- 750')
print(fullfile('Z:\home\lindsey\Analysis\2P\Adaptation\Behavior\allExptAllCellIsiSummary.pdf'),'-dpdf', '-fillpage')

figure;
for icyc = 2:7
    subplot(2,6,icyc-1)
    plot(tt(1:40),squeeze(nanmean(base_resp_off_cyc_avg_all(:,:,:,icyc),2)))
    title(['Cycle ' num2str(icyc)])
    ylim([-0.05 0.05])
    subplot(2,6,icyc+5)
    if icyc == 2
        ylabel('dF/F')
        xlabel('Time (ms)')
    end
    n = sum(~isnan(base_resp_off_cyc_avg_all(1,:,1,icyc)),2);
    errorbar(offs.*(1000/frameRateHz),squeeze(nanmean(mean(base_resp_off_cyc_avg_all(resp_win,:,:,icyc),1),2)),squeeze(nanstd(mean(base_resp_off_cyc_avg_all(resp_win,:,:,icyc),1),[],2))./sqrt(n),'-o');
    ylim([0 0.04])
    xlim([0 1000])
    if icyc == 2
        ylabel('dF/F')
        xlabel('ISI (ms)')
    end
end
suptitle('All cells average response by cyc by ISI- blue- 250; red- 500; yellow- 750')
print(fullfile('Z:\home\lindsey\Analysis\2P\Adaptation\Behavior\allExptIsiByCycSummary.pdf'),'-dpdf', '-fillpage')

figure; 
subplot(1,3,1)
plot(tt(1:40),squeeze(nanmean(targ_resp_off_avg_all,2)))
title('All targets')
xlabel('Time (ms)')
ylabel('dF/F')
ylim([-.02 0.1])
for idelta = 1:nDelta
    subplot(1,3,1+idelta)
    plot(tt(1:40),squeeze(nanmean(targ_resp_off_delta_avg_all(:,:,:,idelta),2)))
    title([num2str(deltas(idelta)) ' deg; n = ' num2str(unique(sum(~isnan(targ_resp_off_delta_avg_all(1,:,:,idelta)),2))')])
    xlabel('Time (ms)')
    ylabel('dF/F')
    ylim([-.02 0.1])
end

suptitle('All cells average response to target by ISI- blue- 250; red- 500; yellow- 750')
print(fullfile('Z:\home\lindsey\Analysis\2P\Adaptation\Behavior\allExptTargByISIByDeltaSummary.pdf'),'-dpdf', '-fillpage')

figure;
for idelta = 1:nDelta
    for i = 1:2
        subplot(2,2,i+ (nDelta.*(idelta-1)))
        plot(tt(1:40),squeeze(nanmean(targ_resp_off_delta_outcome_avg_all(:,:,:,idelta,i),2)))
        if i == 1
            str = 'Hits';
        else
            str = 'Misses';
        end
        xlabel('Time (ms)')
        ylabel('dF/F')
        ylim([-.02 0.1])
        title([num2str(deltas(idelta)) ' deg ' str])
    end
end
suptitle('All cells average response to target by outcome, ori and ISI- blue- 250; red- 500; yellow- 750')
print(fullfile('Z:\home\lindsey\Analysis\2P\Adaptation\Behavior\allExptTargByISIByDeltaByOutcomeSummary.pdf'),'-dpdf', '-fillpage')

figure;
ttAxisTick = find(ismember(floor(tt),-500:500:500));
ttAxisLabel = floor(tt(ttAxisTick));
for idelta = 1:nDelta
    for i = 1:2
        subplot(2,2,i+ (nDelta.*(idelta-1)))
        imagesc(permute(nanmean(targ_resp_off_delta_outcome_avg_all(:,:,:,idelta,i),3),[2 1]))
        colormap(brewermap([],'*RdBu'));
        figXAxis([],'time (ms)',[],ttAxisTick,ttAxisLabel)
        figYAxis([],'cell #',[])
        clim([-0.1 0.1])
        if i == 1
            str = 'Hits';
        else
            str = 'Misses';
        end
        title([num2str(deltas(idelta)) ' deg ' str])
    end
end
suptitle('Response by target and outcome- average all ISI')
print(fullfile('Z:\home\lindsey\Analysis\2P\Adaptation\Behavior\allExptTargByDeltaAndOutcomeSummary.pdf'),'-dpdf', '-fillpage')

    
