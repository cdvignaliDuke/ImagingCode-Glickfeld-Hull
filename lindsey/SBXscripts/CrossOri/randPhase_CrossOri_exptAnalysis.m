clc; clear all; close all;
doRedChannel = 0;
ds = 'CrossOriRandPhase_ExptList';
iexp = 1; 
rc = behavConstsAV;
eval(ds)

frame_rate = 30;

%%
mouse = expt(iexp).mouse;
date = expt(iexp).date;
area = expt(iexp).img_loc{1};
ImgFolder = expt(iexp).coFolder;
time = expt(iexp).coTime;
nrun = length(ImgFolder);
run_str = catRunName(cell2mat(ImgFolder), nrun);

LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
%LG_base = '\\CRASH.dhe.duke.edu\data\home\lindsey';

fprintf(['2P imaging retinotopy analysis\nSelected data:\nMouse: ' mouse '\nDate: ' date '\nExperiments:\n'])
for irun=1:nrun
    fprintf([ImgFolder{irun} ' - ' time{irun} '\n'])
end

%% Test stim analysis
% load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']))
% load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dataStim.mat']))
% load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))
% if doRedChannel == 0
%     red_cells = [];
% end

prewin_frames = frame_rate;
nFramesOn = unique(celleqel2mat_padded(input.nStimOneFramesOn));
postwin_frames = unique(celleqel2mat_padded(input.nStimOneFramesOn));
tt = (1-prewin_frames:postwin_frames).*(1/frame_rate);
data_resp = nan(prewin_frames+postwin_frames,nCells,nTrials);
data_f = nan(1,nCells,nTrials);

for itrial = 1:nTrials
    if cStimOn(itrial) + postwin_frames < sz(3)
        data_resp(:,:,itrial) = npSub_tc(cStimOn(itrial)-prewin_frames:cStimOn(itrial)+postwin_frames-1,:);
    end
end

data_f = mean(data_resp(1:prewin_frames,:,:),1);
data_dfof_tc = (data_resp-data_f)./data_f;

maskPhas_all = celleqel2mat_padded(input.tMaskOneGratingPhaseDeg);
maskPhas = unique(maskPhas_all);
nMaskPhas = length(maskPhas);


nMaskCon = length(maskCons);
nStimCon = length(stimCons);

resp_cell = cell(nMaskCon,nStimCon,nMaskPhas);
trialsperstim = zeros(nMaskCon,nStimCon,nMaskPhas);
h_resp =zeros(nCells,nMaskCon,nStimCon);
p_resp =zeros(nCells,nMaskCon,nStimCon);
resp_win = prewin_frames+5:prewin_frames+nFramesOn;
data_dfof_con_ph_tc_avg = nan(prewin_frames+postwin_frames, nCells, nMaskCon, nStimCon, nMaskPhas);
for im = 1:nMaskCon
    ind_mask = find(maskCon_all == maskCons(im));
    for it = 1:nStimCon
        ind_stim = find(stimCon_all == stimCons(it));
        if it>1 & im>1
            for ip = 1:nMaskPhas
                ind_phas = find(maskPhas_all == maskPhas(ip));
                ind = intersect(ind_phas, intersect(ind_stim,ind_mask));
                trialsperstim(im,it,ip) = length(ind);
                resp_cell{im,it,ip} = squeeze(mean(data_dfof_tc(resp_win,:,ind),1));
                data_dfof_con_ph_tc_avg(:,:,im,it,ip) = squeeze(nanmean(data_dfof_tc(:,:,ind),3));
            end
        else
            ind = intersect(ind_mask,ind_stim);
            trialsperstim(im,it,1) = length(ind);
            resp_cell{im,it,1} = squeeze(mean(data_dfof_tc(resp_win,:,ind),1));
            data_dfof_con_ph_tc_avg(:,:,im,it,1) = squeeze(nanmean(data_dfof_tc(:,:,ind),3));
            [h_resp(:,im,it), p_resp(:,im,it)] = ttest2(resp_cell{im,it,1},resp_cell{1,1,1},'dim',2,'tail','right','alpha',0.05./(nMaskCon*nStimCon));
        end
    end
end

test_resp = resp_cell{1,nStimCon};
mask_resp = resp_cell{nMaskCon,1};
plaid_resp = resp_cell{nMaskCon,nStimCon};
null_resp = resp_cell{1,1};

[h_resptest p_resptest] =  ttest2(test_resp,null_resp,'dim',2,'tail','right');
[h_respmask p_respmask] = ttest2(mask_resp,null_resp,'dim',2,'tail','right');
[h_respplaid p_respplaid] = ttest2(plaid_resp,null_resp,'dim',2,'tail','right');
resptest_ind = find(h_resptest);
respmask_ind = find(h_respmask);
respplaid_ind = find(h_respplaid);

resptest_ind = setdiff(resptest_ind,red_cells);
respmask_ind = setdiff(respmask_ind,red_cells);
respplaid_ind = setdiff(respplaid_ind,red_cells);

resp_ind = unique([resptest_ind; respmask_ind; respplaid_ind]);
resp_ind_nomask = unique([resptest_ind; respplaid_ind]);

[h_preftest p_preftest] = ttest2(test_resp,mask_resp,'dim',2,'tail','right');
[h_prefmask p_prefmask] = ttest2(test_resp,mask_resp,'dim',2,'tail','left');
[h_prefplaid1 p_prefplaid1] = ttest2(plaid_resp, test_resp ,'dim',2,'tail','right');
[h_prefplaid2 p_prefplaid2] = ttest2(plaid_resp, mask_resp ,'dim',2,'tail','right');
preftest_ind = intersect(resp_ind,find(h_preftest));
prefmask_ind = intersect(resp_ind,find(h_prefmask));
prefplaid_ind = intersect(resp_ind,intersect(find(h_prefplaid1),find(h_prefplaid2)));
prefplaid_nonlinear = find(mean(mask_resp,2)+mean(test_resp,2)<mean(plaid_resp,2));
prefplaidonly_ind =  intersect(prefplaid_ind,prefplaid_nonlinear);
preftestonly_ind = setdiff(preftest_ind,prefplaidonly_ind);
prefmaskonly_ind = setdiff(prefmask_ind,prefplaidonly_ind);

preftestonly_ind = setdiff(preftestonly_ind,red_cells);
prefmaskonly_ind = setdiff(prefmaskonly_ind,red_cells);
prefplaidonly_ind = setdiff(prefplaidonly_ind,red_cells);

save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']), 'data_dfof_con_ph_tc_avg', 'resp_cell', 'data_dfof_tc', 'resp_ind', 'preftest_ind', 'prefmask_ind', 'prefplaid_ind', 'preftestonly_ind', 'prefmaskonly_ind', 'prefplaidonly_ind', 'tt', 'frame_rate');
save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']), 'prewin_frames', 'postwin_frames', 'resp_win', 'stimCon_all', 'maskCon_all', 'maskPhas_all', 'stimCons', 'maskCons', 'maskPhas', 'nStimCon', 'nMaskCon', 'nMaskPhas');

%% plots
tf = input.stimOneGratingTemporalFreqCPS;
for i = 1:nCells
    start = 3;
    for im = 1:nMaskCon
        for it = 1:nStimCon
            if im==1 & it == 2 & find(resptest_ind == i)
                temp  = squeeze(data_dfof_con_ph_tc_avg(:,i,:,:,:));
                max_val = nanmax(temp(:));
                min_val = nanmin(temp(:));
                figure;
                subplot(6,2,1)
                test_alone = data_dfof_con_ph_tc_avg(:,i,im,it,1);
                plot(tt,test_alone);
                hold on
                ylim([min_val max_val])
                xlim([-1.5 6.5])
                tt2 = double((1-(frame_rate./2):frame_rate./tf)).*(1/frame_rate);
                test_alone_cycle = zeros(size(tt2,2),3);
                for ic = 1:3
                    test_alone_cycle(:,ic) = test_alone((prewin_frames-frame_rate./2)+((ic).*(frame_rate./tf)):...
                        (prewin_frames+frame_rate./tf)+((ic).*(frame_rate./tf))-1,:);
                end
                
                subplot(6,2,2)
                plot(tt2, mean(test_alone_cycle,2)-mean(test_alone_cycle(frame_rate./2,:),2))
                hold on
            end
            if im==3 & it == 1 & find(resptest_ind == i)
                subplot(6,2,1)
                mask_alone = data_dfof_con_ph_tc_avg(:,i,im,it,1);
                plot(tt,mask_alone);
                mask_alone_cycle = zeros(size(tt2,2),3);
                for ic = 1:3
                    mask_alone_cycle(:,ic) = mask_alone((prewin_frames-frame_rate./2)+((ic).*(frame_rate./tf)):...
                        (prewin_frames+frame_rate./tf)+((ic).*(frame_rate./tf))-1,:);
                end
                ylim([min_val max_val])
                xlim([-1.5 6.5])
                title(['Test- ' num2str(stimCons(2)) ' ;Mask- '  num2str(maskCons(3))])
                subplot(6,2,2)
                plot(tt2, mean(mask_alone_cycle,2)-mean(mask_alone_cycle(frame_rate./2,:),2))
                ylim([min_val/2 max_val/2])
                xlim([-.5 1])
                title('Cycle average')
            end
            if im>1 & it>1
                for ip = 1:nMaskPhas
                    if im == 3 & it == 2 & find(resptest_ind == i)
                        subplot(6,2,start)
                        plot(tt,test_alone);
                        hold on
                        plot(tt,data_dfof_con_ph_tc_avg(:,i,im,it,ip));
                        temp =data_dfof_con_ph_tc_avg(:,i,im,it,ip);
                        phase_cycle = zeros(size(tt2,2),3);
                        for ic = 1:3
                            phase_cycle(:,ic) = temp((prewin_frames-frame_rate./2)+((ic).*(frame_rate./tf)):...
                                (prewin_frames+frame_rate./tf)+((ic).*(frame_rate./tf))-1,:);
                        end
                        ylim([min_val max_val])
                        xlim([-1.5 6.5])
                        title(['Test- ' num2str(stimCons(it)) ';Mask- '  num2str(maskPhas(ip))])
                        subplot(6,2,start+1)
                        plot(tt2, mean(test_alone_cycle,2)-mean(test_alone_cycle(frame_rate./2,:),2))
                        hold on
                        plot(tt2, mean(phase_cycle,2)-mean(phase_cycle(frame_rate./2,:),2))
                        ylim([min_val max_val/4])
                        xlim([-.5 1])
                        start = start+2;
                        subplot(6,2,12)
                        plot(tt2, mean(phase_cycle,2)-mean(phase_cycle(frame_rate./2,:),2))
                        hold on
                        ylim([min_val max_val/4])
                        xlim([-.5 1])
                        
                    end
                end
            end
        end
    end
    if find(resptest_ind == i)
        suptitle(['Cell ' num2str(i)])
        print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_compTestMaskByPhase_Cell' num2str(i) '.pdf']), '-dpdf', '-bestfit')
    end
end


figure;
start = 1;
prefTest_resp = zeros(nMaskCon,nStimCon,nMaskPhas,2);
for im = 1:nMaskCon
    ind_mask = find(maskCon_all == maskCons(im));
    for it = 1:nStimCon
        ind_stim = find(stimCon_all == stimCons(it));
        ind = intersect(ind_stim,ind_mask);
        subplot(nMaskCon,nStimCon,start)
        if im>1 & it>1
            for ip = 1:nMaskPhas
                ind_phas = find(maskPhas_all == maskPhas(ip));
                ind_p = intersect(ind,ind_phas);
                plot(tt,nanmean(nanmean(data_dfof_tc(:,preftest_ind,ind_p),2),3))
                prefTest_resp(im,it,ip,1) = nanmean(nanmean(resp_cell{im,it,ip}(preftest_ind,:),1),2);
                prefTest_resp(im,it,ip,2) = nanstd(nanmean(resp_cell{im,it,ip}(preftest_ind,:),2),[],1)./sqrt(length(preftest_ind));
                hold on
            end
        else
            plot(tt,nanmean(nanmean(data_dfof_tc(:,preftest_ind,ind),2),3))
        end
        title(['T: ' num2str(stimCons(it)) '; M: ' num2str(maskCons(im))])
        ylim([-0.1 .8])
        xlabel('Time (s)')
        start = start+1;
    end
end
suptitle(['Test Pref (n = ' num2str(length(preftest_ind)) ') resp by phase'])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCbyPhase_preftest.pdf']),'-dpdf','-bestfit');

figure;
start = 1;
prefPlaid_resp = zeros(nMaskCon,nStimCon,nMaskPhas,2);
for im = 1:nMaskCon
    ind_mask = find(maskCon_all == maskCons(im));
    for it = 1:nStimCon
        ind_stim = find(stimCon_all == stimCons(it));
        ind = intersect(ind_stim,ind_mask);
        subplot(nMaskCon,nStimCon,start)
        if im>1 & it>1
            for ip = 1:nMaskPhas
                ind_phas = find(maskPhas_all == maskPhas(ip));
                ind_p = intersect(ind,ind_phas);
                plot(tt,nanmean(nanmean(data_dfof_tc(:,prefplaid_ind,ind_p),2),3))
                prefPlaid_resp(im,it,ip,1) = nanmean(nanmean(resp_cell{im,it,ip}(prefplaid_ind,:),1),2);
                prefPlaid_resp(im,it,ip,2) = nanstd(nanmean(resp_cell{im,it,ip}(prefplaid_ind,:),2),[],1)./sqrt(length(prefplaid_ind));
                hold on
            end
        else
            plot(tt,nanmean(nanmean(data_dfof_tc(:,prefplaid_ind,ind),2),3))
        end
        title(['T: ' num2str(stimCons(it)) '; M: ' num2str(maskCons(im))])
        ylim([-0.1 .8])
        xlabel('Time (s)')
        start = start+1;
    end
end
suptitle(['Plaid Pref (n = ' num2str(length(prefplaid_ind)) ') resp by phase'])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCbyPhase_prefplaid.pdf']),'-dpdf','-bestfit');


%%
prefTestOrPlaid = unique([preftest_ind;prefplaid_ind]);
nc = length(prefTestOrPlaid);
if nc>25
    nc = 25;
end
[n n2] = subplotn(nc);
figure; 
for iC = 1:nc
    iCell = prefTestOrPlaid(iC);
    subplot(n,n2,iC)
    plot(tt, data_dfof_con_ph_tc_avg(:,iCell,1,2,1))
    ylabel('dF/F')
    xlabel('Time (s)')
    title(num2str(iCell))
end
suptitle([num2str(nc) ' test or plaid preferring cells- Test = ' num2str(stimCons(2)) ' Mask = ' num2str(maskCons(1))]) %'- n = ' num2str(length(prefTestOrPlaid))])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_testTCs_preftestorplaid_first25.pdf']),'-dpdf','-fillpage');

figure; 
for iC = 1:nc
    iCell = prefTestOrPlaid(iC);
    subplot(n,n2,iC)
    for ip =1:4;
        plot(tt, data_dfof_con_ph_tc_avg(:,iCell,3,2,ip))
        hold on
    end
    ylabel('dF/F')
    xlabel('Time (s)')
    title(num2str(iCell))
end
suptitle([num2str(nc) ' test or plaid preferring cells- Test = ' num2str(stimCons(2)) ' Mask = ' num2str(maskCons(3))]) %'- n = ' num2str(length(prefTestOrPlaid))])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_plaidTCsByPhase_preftestorplaid_first25.pdf']),'-dpdf','-fillpage');


figure;
[n n2] = subplotn(length(prefTestOrPlaid));
for iC = 1:nc
    iCell = prefTestOrPlaid(iC);
    subplot(n,n2,iC)
    test_avg = mean(resp_cell{1,2,1}(iCell,:),2);
    mask_avg = mean(resp_cell{end,1,1}(iCell,:),2);
    test_sem = mean(resp_cell{1,2,1}(iCell,:),2)./sqrt(size(resp_cell{1,end,1}(iCell,:),2));
    mask_sem = mean(resp_cell{end,1,1}(iCell,:),2)./sqrt(size(resp_cell{end,1,1}(iCell,:),2));
    shadedErrorBar(1:360,repmat(test_avg,[1 360]),repmat(test_sem,[1 360]));
    hold on
    shadedErrorBar(1:360,repmat(mask_avg,[1 360]),repmat(mask_sem,[1 360]));
    resp_all = [];
    stim_all = [];
    for ip = 1:nMaskPhas
        resp_avg(1,ip) = mean(resp_cell{end,2,ip}(iCell,:),2);
        resp_sem(1,ip) = std(resp_cell{end,2,ip}(iCell,:),[],2)./sqrt(size(resp_cell{end,2,ip}(iCell,:),2));
        resp_all = [resp_all resp_cell{end,2,ip}(iCell,:)];
        stim_all = [stim_all ip.*ones(size(resp_cell{end,2,ip}(iCell,:)))];
    end
    errorbar(maskPhas,resp_avg,resp_sem)
    p = anova1(resp_all, stim_all,'off');
    title(num2str(chop(p,2)))
    ylabel('dF/F')
    xlabel('Phase (deg)')
end
suptitle([num2str(nc) ' test or plaid preferring cells- Test = ' num2str(stimCons(2)) ' Mask = ' num2str(maskCons(3))]) %'- n = ' num2str(length(prefTestOrPlaid))])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respByPhase_preftestorplaid_first25.pdf']),'-dpdf','-fillpage');
 
figure;
start = 1;
for it = 2:3
    for im = 2:3
        subplot(2,2,start)
        test_avg = mean(mean(resp_cell{1,it,1}(prefTestOrPlaid,:),2),1);
        mask_avg = mean(mean(resp_cell{im,1,1}(prefTestOrPlaid,:),2),1);
        test_sem = std(mean(resp_cell{1,it,1}(prefTestOrPlaid,:),2),[],1)./sqrt(length(prefTestOrPlaid));
        mask_sem = std(mean(resp_cell{im,1,1}(prefTestOrPlaid,:),2),[],1)./sqrt(length(prefTestOrPlaid));
        shadedErrorBar(1:360,repmat(test_avg,[1 360]),repmat(test_sem,[1 360]));
        hold on
        shadedErrorBar(1:360,repmat(mask_avg,[1 360]),repmat(mask_sem,[1 360]));
        resp_avg_all = zeros(length(prefTestOrPlaid),nMaskPhas);
        for ip = 1:nMaskPhas
            resp_avg(1,ip) = mean(mean(resp_cell{im,it,ip}(prefTestOrPlaid,:),2),1);
            resp_sem(1,ip) = std(mean(resp_cell{im,it,ip}(prefTestOrPlaid,:),2),[],1)./sqrt(length(prefTestOrPlaid));
            resp_avg_all(:,ip) = mean(resp_cell{im,it,ip}(prefTestOrPlaid,:),2);
        end
        errorbar(maskPhas,resp_avg,resp_sem)
        ylabel('dF/F')
        xlabel('Phase (deg)')
        ylim([-0.1 0.4])
        [p tbl] = anova1(resp_avg_all,[],'off');
        title(['Test = ' num2str(stimCons(it)) '; Mask = ' num2str(maskCons(im)) '; p = ' num2str(chop(p,3))])
        start = start+1;
    end
end
suptitle(['All test or plaid preferring cells- n = ' num2str(length(prefTestOrPlaid))])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgRespByPhase_preftestorplaid.pdf']),'-dpdf','-bestfit');

figure;
start = 1;
for it = 2:3
    for im = 2:3
        subplot(2,2,start)
        test_avg = mean(mean(resp_cell{1,it,1}(prefplaid_ind,:),2),1);
        mask_avg = mean(mean(resp_cell{im,1,1}(prefplaid_ind,:),2),1);
        test_sem = std(mean(resp_cell{1,it,1}(prefplaid_ind,:),2),[],1)./sqrt(length(prefplaid_ind));
        mask_sem = std(mean(resp_cell{im,1,1}(prefplaid_ind,:),2),[],1)./sqrt(length(prefplaid_ind));
        shadedErrorBar(1:360,repmat(test_avg,[1 360]),repmat(test_sem,[1 360]));
        hold on
        shadedErrorBar(1:360,repmat(mask_avg,[1 360]),repmat(mask_sem,[1 360]));
        resp_avg_all = zeros(length(prefplaid_ind),nMaskPhas);
        for ip = 1:nMaskPhas
            resp_avg(1,ip) = mean(mean(resp_cell{im,it,ip}(prefplaid_ind,:),2),1);
            resp_sem(1,ip) = std(mean(resp_cell{im,it,ip}(prefplaid_ind,:),2),[],1)./sqrt(length(prefplaid_ind));
            resp_avg_all(:,ip) = mean(resp_cell{im,it,ip}(prefplaid_ind,:),2);
        end
        errorbar(maskPhas,resp_avg,resp_sem)
        ylabel('dF/F')
        xlabel('Phase (deg)')
        ylim([-0.1 0.6])
        [p tbl] = anova1(resp_avg_all,[],'off');
        title(['Test = ' num2str(stimCons(it)) '; Mask = ' num2str(maskCons(im)) '; p = ' num2str(chop(p,3))])
        start = start+1;
    end
end
suptitle(['All plaid preferring cells- n = ' num2str(length(prefplaid_ind))])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgRespByPhase_prefplaid.pdf']),'-dpdf','-bestfit');

figure;
start = 1;
for it = 2:3
    for im = 2:3
        subplot(2,2,start)
        test_avg = mean(mean(resp_cell{1,it,1}(preftest_ind,:),2),1);
        mask_avg = mean(mean(resp_cell{im,1,1}(preftest_ind,:),2),1);
        test_sem = std(mean(resp_cell{1,it,1}(preftest_ind,:),2),[],1)./sqrt(length(preftest_ind));
        mask_sem = std(mean(resp_cell{im,1,1}(preftest_ind,:),2),[],1)./sqrt(length(preftest_ind));
        shadedErrorBar(1:360,repmat(test_avg,[1 360]),repmat(test_sem,[1 360]));
        hold on
        shadedErrorBar(1:360,repmat(mask_avg,[1 360]),repmat(mask_sem,[1 360]));
        resp_avg_all = zeros(length(preftest_ind),nMaskPhas);
        for ip = 1:nMaskPhas
            resp_avg(1,ip) = mean(mean(resp_cell{im,it,ip}(preftest_ind,:),2),1);
            resp_sem(1,ip) = std(mean(resp_cell{im,it,ip}(preftest_ind,:),2),[],1)./sqrt(length(preftest_ind));
            resp_avg_all(:,ip) = mean(resp_cell{im,it,ip}(preftest_ind,:),2);
        end
        errorbar(maskPhas,resp_avg,resp_sem)
        ylabel('dF/F')
        xlabel('Phase (deg)')
        ylim([-0.1 0.6])
        [p tbl] = anova1(resp_avg_all,[],'off');
        title(['Test = ' num2str(stimCons(it)) '; Mask = ' num2str(maskCons(im)) '; p = ' num2str(chop(p,3))])
        start = start+1;
    end
end
suptitle(['All test preferring cells- n = ' num2str(length(preftest_ind))])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgRespByPhase_preftest.pdf']),'-dpdf','-bestfit');

