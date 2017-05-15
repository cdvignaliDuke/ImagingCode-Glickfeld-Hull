clear all
close all
ds = 'movDotsSpeedTun_V1';
slct_expt = [1];
%%
rc = behavConstsAV;
eval(ds)
for iexp = slct_expt
    subNum = expt(iexp).SubNum;
    mouse = expt(iexp).mouse;
    expDate = expt(iexp).date;
    down = expt(iexp).downSampleRate;
    for irun = 1:expt(iexp).nrun
        runFolder = expt(iexp).runs(irun,:);
        runTime = expt(iexp).time_mat(irun,:);
        fName = [runFolder '_000_000'];
        
        % load cell timecourses and mworks file
        fn = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging', expDate,runFolder);
        load(fullfile(fn,'timecourses.mat'))
        input = Load_SBXdataPlusMWorksData(subNum,expDate,runTime,mouse,runFolder,fName);
        
        % img params
        nfr = size(data_tc_subnp,1);
        nc = size(data_tc_subnp,2);
        frRateHz = expt(iexp).frame_rate/down;
        
        % vis stim params
        on = input.nScansOn./down;
        off = input.nScansOff./down;
        
        tSpeed = chop(cell2mat(input.tDotSpeedDPS),3);
        [tSpeed_ind, speeds] = findgroups(tSpeed);
        ntrials = length(tSpeed);
        nstim = length(speeds);
        
        % analysis params
        resp_delay = round(frRateHz*0.2);
        base_win = off-frRateHz+1:off;
        resp_win = off+resp_delay:off+on; 
        
        % align length of experiment with last full-length trial
        if (off+on)*ntrials > nfr 
            ntrials = floor(nfr./(off+on));
            [tSpeed_ind speeds] = findgroups(tSpeed(1:ntrials));
            data = data_tc_subnp(1:((off+on)*ntrials),:);
        elseif (off+on)*ntrials < nfr
            ntrials = length(tSpeed);
            [tSpeed_ind speeds] = findgroups(tSpeed(1:ntrials));
            data = data_tc_subnp(1:((off+on)*ntrials),:);
        else
            data = data_tc_subnp;
        end
        
        % reshape data into each trial
        data_tr = reshape(data,off+on,ntrials,nc);
        
        % dF/F
        f = mean(data_tr(base_win,:,:),1);
        dff = bsxfun(@rdivide, bsxfun(@minus, data_tr, f), f);
        
        %% tuning curves and response reliability
        
        % sort trials by speed, test response reliability
        trial_tc = cell(1,nstim);
        resp_ttest = zeros(nstim,nc);
        for istim = 1:nstim
            ind = tSpeed_ind == istim;
            trial_tc{istim} = dff(:,ind,:);
            resp_ttest(istim,:) = ttest(mean(dff(resp_win,ind,:),1),mean(dff(base_win,ind,:),1),'dim',2,'tail','right');%
        end;
        
        % average response and time-course for each cell, for each speed
        speed_tc = cellfun(@(x) squeeze(mean(x,2)), trial_tc,'unif',0);
        speed_resp = cell2mat(cellfun(@(x) mean(x(resp_win,:),1)-mean(x(base_win,:),1), speed_tc,'unif',0)');
        speed_resp_err = cell2mat(cellfun(@(x) ste(x(resp_win,:),1), speed_tc,'unif',0)');
        
        % rectify tuning curves at non-significant speeds for each cell
        speedTuning_rect = speed_resp;
        speedTuning_rect(resp_ttest == 0) = 0;
        speedTuning_rect = reshape(speedTuning_rect,nstim,nc);
        speed_err_rect = speed_resp_err;
        speed_err_rect(resp_ttest == 0) = 0;
        speed_err_rect = reshape(speed_err_rect,nstim,nc);
        
        % plot example cells tuning 
        exCells = randsample(nc,16);
        
        figure; setFigParams4Print('landscape')
        suptitle({[expt(iexp).dataset_exp{irun} ' speed tuning'];'ns responses in gray'})
        for iplot = 1:16
            subplot(4,4,iplot)
            
            tempresp = speed_resp(:,exCells(iplot));
            temperr = speed_resp_err(:,exCells(iplot));
            tempttest = logical(resp_ttest(:,exCells(iplot)));
            h = errorbar(speeds,tempresp,temperr,'ko');
            hold on
            h = errorbar(speeds(~tempttest),tempresp(~tempttest),temperr(~tempttest),'ko');
            h.Color = [0.7 0.7 0.7];
            figXAxis(h.Parent,'speed',[0 speeds(end)+10],[1,10,50],[1,10,50])
            ylabel('dF/F')
%             h.Parent.TickDir = 'out';
            h.Parent.Box = 'off';
            h.Parent.XScale = 'log';
            title(['cell ' num2str(exCells(iplot))])
        end
        print(fullfile(fn,'exCellsTuning'),'-dpdf','-fillpage')
        
    end
end