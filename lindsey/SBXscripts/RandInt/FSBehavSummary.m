%%
mouse_mat = {'i613','i614','i668'};
date_mat{1} = {'150511'};
date_mat{2} = {'150623','150626'};
date_mat{3} = {'161031'};
run_mat{1} = {'001-003'};
run_mat{2} = {'002-004','001-003'};
run_mat{3} = {'001-003'};
dir_mat{1} = {'005'};
dir_mat{2} = {'006','005'};
dir_mat{3} = {'006'};

%%

t = 1;
for imouse = 1:size(mouse_mat,2)
    mouse = mouse_mat{imouse};
    for iexp = 1:size(date_mat{imouse},2)
        date = date_mat{imouse}{iexp};
        run = run_mat{imouse}{iexp};
        
        load(fullfile('Z:\home\lindsey\Analysis\2P',[date '_' mouse],[date '_' mouse '_runs-' run],[date '_' mouse '_runs-' run '_input.mat']));
        load(fullfile('Z:\home\lindsey\Analysis\2P',[date '_' mouse],[date '_' mouse '_runs-' run],[date '_' mouse '_runs-' run '_dfofData.mat']));
        load(fullfile('Z:\home\lindsey\Analysis\2P',[date '_' mouse],[date '_' mouse '_runs-' run],[date '_' mouse '_runs-' run '_stimData.mat']));

        b2 = celleqel2mat_padded(input.tBlock2TrialNumber);
        tGratingDir = celleqel2mat_padded(input.tGratingDirectionDeg);
        tGratingDir(find(b2))=[];
        dirs = unique(tGratingDir);
        ind_easy = find(tGratingDir == dirs(end));
        reactTimes = celleqel2mat_padded(input.reactTimesMs);
        reactTimes(find(b2))=[];
        holdTimes = celleqel2mat_padded(input.holdTimesMs);
        holdTimes(find(b2))=[];
        tCyc = celleqel2mat_padded(input.tCyclesOn);
        tCyc(find(b2))= [];
        nCyc = celleqel2mat_padded(input.nCyclesOn);
        nCyc(find(b2))= [];
        
        data_dfof(:,:,:,find(b2)) = [];
        sz = size(data_dfof);
        nTrials = sz(4);
        react_win = [200 550];
        base_dfof = nan(sz);
        if iscell(input.nFramesOn)
            nOn = unique(cell2mat(input.nFramesOn));
            nOff = unique(cell2mat(input.nFramesOff));
        else
            nOn = unique(input.nFramesOn);
            nOff = unique(input.nFramesOff);
        end
        
        SIx = zeros(1,nTrials);
        MIx = zeros(1,nTrials);
        FAIx = zeros(1,nTrials);
        FA_resp = nan(sz(2),sz(4));
        preFA_resp = nan(sz(2),sz(4));
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
                        end
                    else
                        totTime = (sum(tFramesOff(itrial,1:tCyc(itrial)-1),2)+((tCyc(itrial)-1).*nOn)).*(1000/frameRateHz);
                        if holdTimes(itrial)-totTime>react_win(1) & holdTimes(itrial)-totTime<react_win(2)
                            if tCyc(itrial)-3 > 0
                                FAIx(itrial)=1;
                                FACyc(itrial) = tCyc(itrial)-1;
                                FA_resp(:,itrial) = squeeze(mean(data_dfof(resp_win,:,tCyc(itrial)-1,itrial),1)-mean(data_dfof(base_win,:,tCyc(itrial)-1,itrial),1));
                                preFA_resp(:,itrial) = squeeze(mean(data_dfof(resp_win,:,tCyc(itrial)-2,itrial),1)-mean(data_dfof(base_win,:,tCyc(itrial)-2,itrial),1));
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
        for itrial = 1:nTrials
            if ~isnan(targCyc(itrial))
                data_dfof_targ(:,:,itrial) = data_dfof(:,:,targCyc(itrial),itrial);
            end
        end
        targCyc(find(FAIx)) = FACyc(find(FAIx));
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
        expt(t).targetOrientation = tGratingDir;
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
            load(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_runs-' dir_run], ['oriTuningAndFits.mat']))
            [~,max_ori] = max(vonMisesFitAllCells,[],1);
            expt(t).oriTuning = max_ori(1,:)-1;
            expt(t).oriTuningTheta90 = fitReliability;
        else
            expt(t).oriTuning = [];
            expt(t).oriTuningTheta90 = [];
        end
        
        t = 1+t;
    end
end
save(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P\Adaptation', 'FSBehavData.mat'),'expt')