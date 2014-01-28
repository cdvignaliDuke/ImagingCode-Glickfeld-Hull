fn = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_reps.mat']);
load(fn);
fn = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_roi_avg.mat']);
load(fn);
if run == 1;
    fn = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_stim_reps_run.mat']);
    load(fn);
    fn = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_stim_reps_norun.mat']);
    load(fn);
end

nReps = sum(stim_reps(1,:));

roi_stim = zeros(nON+nOFF, size(roi_avg,2),nReps);
start = 1;
rep = 1;
if blanks == 1
    for iCond = 1:nCond+1; 
        nRep = stim_reps(iCond);
        for iRep = 1:nRep;
            roi_stim(:,:, rep) = roi_avg(start:start-1+nON+nOFF,:);
            start = start+nON+nOFF;
            rep = rep+1;
        end
    end
else
    for iCond = 1:nCond; 
        nRep = stim_reps(iCond);
        for iRep = 1:nRep;
            roi_stim(:,:, rep) = roi_avg(start:start-1+nON+nOFF,:);
            start = start+nON+nOFF;
            rep = rep+1;
        end
    end
end

resp_off = mean(roi_stim(pre_win(1):pre_win(2),:,:),1);
resp_on = squeeze(mean(roi_stim(post_win(1):post_win(2),:,:),1));
if size(roi_stim,2)==1
    resp_on =resp_on';
end

if run == 1;
    ind_fn = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_running_ind.mat']);
    load(ind_fn);
    resp_off_run = mean(roi_stim(pre_win(1):pre_win(2),:,run_ind),1);
    resp_off_norun = mean(roi_stim(pre_win(1):pre_win(2),:,norun_ind),1);
    resp_on_run = squeeze(mean(roi_stim(post_win(1):post_win(2),:,run_ind),1));
    resp_on_norun = squeeze(mean(roi_stim(post_win(1):post_win(2),:,norun_ind),1));
end

TC_norm = bsxfun(@minus, roi_stim, resp_off);
TC_dF = 100*(bsxfun(@rdivide, TC_norm,resp_off));
if size(roi_stim,2)>1
    resp_dF = squeeze(mean(double(TC_dF(post_win(1):post_win(2),:,:)),1));
    resp_norm = squeeze(mean(double(TC_norm(post_win(1):post_win(2),:,:)),1));
else
    resp_dF = mean(double(TC_dF(post_win(1):post_win(2),:,:)),1);
    resp_norm = squeeze(mean(double(TC_norm(post_win(1):post_win(2),:,:)),1));
end


fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_resp_POST' num2str(post_win) '.mat']);
save(fn_out, 'resp_on', 'resp_off', 'TC_norm', 'TC_dF', 'resp_dF', 'resp_norm');
if run ==1
    fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_resp_POST' num2str(post_win) '_run_norun.mat']);
    save(fn_out,'resp_on_run', 'resp_off_run', 'resp_on_norun', 'resp_off_norun');
end



resp_avg = [];
resp_sem = [];
TC_avg = [];
TC_sem = [];
TC_norm_avg = [];
TC_norm_sem = [];

if run == 1;
    nRunning = 3;
    ind_fn = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_running_ind.mat']);
    load(ind_fn);
else
    nRunning = 1;
end
for iRunning = 1:nRunning;
    start =1;
    if blanks == 1;
        for iCond = 1:nCond+1;
            if iRunning == 1;
                nRep = stim_reps(1,iCond);
                ind = 1:nReps;
            elseif iRunning == 2;
                nRep = stim_reps_run(1,iCond);
                ind = run_ind;
            elseif iRunning == 3;
                nRep = stim_reps_norun(1,iCond);
                ind = norun_ind;
            end
            TC_avg(:,:,iCond,iRunning) = mean(double(TC_dF(:,:,ind(start:start-1+nRep))),3);
            TC_sem(:,:,iCond,iRunning) = std(TC_dF(:,:,ind(start:start-1+nRep)),[],3)./sqrt(size(TC_dF(:,:,ind(start:start-1+nRep)),3));
            TC_norm_avg(:,:,iCond,iRunning) = mean(double(TC_norm(:,:,ind(start:start-1+nRep))),3);
            TC_norm_sem(:,:,iCond,iRunning) = std(TC_norm(:,:,ind(start:start-1+nRep)),[],3)./sqrt(size(TC_norm(:,:,ind(start:start-1+nRep)),3));
            if P ==1
                resp_avg(:,iCond,iRunning) = mean(mean(double(TC_dF(post_win(1):post_win(2),:,ind(start:start-1+nRep))),1),3);
                resp_sem(:,iCond,iRunning) = std(mean(TC_dF(post_win(1):post_win(2),:,ind(start:start-1+nRep)),1),[],3)./sqrt(size(TC_dF(post_win(1):post_win(2),:,ind(start:start-1+nRep)),3));
            elseif P==2
                resp_avg(:,iCond,iRunning) = mean(mean(double(TC_norm(post_win(1):post_win(2),:,ind(start:start-1+nRep))),1),3);
                resp_sem(:,iCond,iRunning) = std(mean(TC_norm(post_win(1):post_win(2),:,ind(start:start-1+nRep)),1),[],3)./sqrt(size(TC_norm(post_win(1):post_win(2),:,ind(start:start-1+nRep)),3));
            end
                start= start+nRep;
        end
    else
        for iCond = 1:nCond;
            if iRunning == 1;
                nRep = stim_reps(1,iCond);
                ind = 1:nReps;
            elseif iRunning == 2;
                nRep = stim_reps_run(1,iCond);
                ind = run_ind;
            elseif iRunning == 3;
                nRep = stim_reps_norun(1,iCond);
                ind = norun_ind;
            end
            TC_avg(:,:,iCond,iRunning) = mean(double(TC_dF(:,:,ind(start:start-1+nRep))),3);
            TC_sem(:,:,iCond,iRunning) = std(TC_dF(:,:,ind(start:start-1+nRep)),[],3)./sqrt(size(TC_dF(:,:,ind(start:start-1+nRep)),3));
            TC_norm_avg(:,:,iCond,iRunning) = mean(double(TC_norm(:,:,ind(start:start-1+nRep))),3);
            TC_norm_sem(:,:,iCond,iRunning) = std(TC_norm(:,:,ind(start:start-1+nRep)),[],3)./sqrt(size(TC_norm(:,:,ind(start:start-1+nRep)),3));
            if P ==1
                resp_avg(:,iCond,iRunning) = mean(mean(double(TC_dF(post_win(1):post_win(2),:,ind(start:start-1+nRep))),1),3);
                resp_sem(:,iCond,iRunning) = std(mean(TC_dF(post_win(1):post_win(2),:,ind(start:start-1+nRep)),1),[],3)./sqrt(size(TC_dF(post_win(1):post_win(2),:,ind(start:start-1+nRep)),3));
            elseif P==2
                resp_avg(:,iCond,iRunning) = mean(mean(double(TC_norm(post_win(1):post_win(2),:,ind(start:start-1+nRep))),1),3);
                resp_sem(:,iCond,iRunning) = std(mean(TC_norm(post_win(1):post_win(2),:,ind(start:start-1+nRep)),1),[],3)./sqrt(size(TC_norm(post_win(1):post_win(2),:,ind(start:start-1+nRep)),3));
            end
            start= start+nRep;
        end
    end
end

fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_resp_norm_POST' num2str(post_win) '.mat']);
save(fn_out, 'resp_avg', 'resp_sem', 'TC_avg', 'TC_sem', 'TC_norm_avg', 'TC_norm_sem');
