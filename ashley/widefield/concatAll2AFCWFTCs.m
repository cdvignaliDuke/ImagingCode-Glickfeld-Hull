function concatAll2AFCWFTCs(mouse)
    rc = behavConstsWF(mouse);
    xd = frm_xls2frm(rc.indexFilename, [], rc.indexTextCols);
    
    dateStr_all = [];
    SIx_all = [];
    MIx_all = [];
    FIx_all = [];
    NIx_all = [];
    tLeftTrial_all = [];
    tGratingContrast_all = [];
    tProbLeft_all = [];
    data_stim_tc_all = [];
    data_dec_tc_all = [];
    for indexRowN = 1:xd.nRows
        % do the one we found
        dateStr = xd.DateStr{indexRowN};
        fprintf('Starting subject %s, date %s\n', ...
        mouse, dateStr);
        load(fullfile(rc.structOutput, [mouse '_' dateStr '_tc.mat']));
        timeStr = num2str(eval(['xd.MatFileRun1' '(indexRowN)']));
        fn_mworks = fullfile(rc.pathStr,['data-' mouse '-' dateStr '-' timeStr '.mat']);
        input = mwLoadData(fn_mworks, [], []);

        tLeftTrial = celleqel2mat_padded(input.tLeftTrial);
        tProbLeft = celleqel2mat_padded(input.tStimProbAvgLeft);
        tGratingContrast = celleqel2mat_padded(input.tGratingContrast);
        cStimOn = celleqel2mat_padded(input.cStimOn);
        cDecision = celleqel2mat_padded(input.cDecision);
        SIx = strcmp(input.trialOutcomeCell, 'success');
        FIx = strcmp(input.trialOutcomeCell, 'incorrect');
        MIx = strcmp(input.trialOutcomeCell, 'ignore');
        NIx = celleqel2mat_padded(input.isNoGo);
        nTrials = length(tLeftTrial);
        sz = size(roiTC);
        data_stim_tc = nan(sz(2),40,nTrials);
        data_dec_tc = nan(sz(2),40,nTrials);
        for itrial = 1:nTrials
            if cStimOn(itrial)+30<sz(1)
                data_stim_tc(:,:,itrial) = roiTC(cStimOn(itrial)-10:cStimOn(itrial)+29,:)';
            end
            if cDecision(itrial)+30<sz(1)
                data_dec_tc(:,:,itrial) = roiTC(cDecision(itrial)-10:cDecision(itrial)+29,:)';
            end
        end
        dateStr_all = strvcat(dateStr_all, dateStr);
        SIx_all = [SIx_all SIx];
        MIx_all = [MIx_all MIx];
        FIx_all = [FIx_all FIx];
        NIx_all = [NIx_all NIx];
        tLeftTrial_all = [tLeftTrial_all tLeftTrial];
        tGratingContrast_all = [tGratingContrast_all tGratingContrast];
        tProbLeft_all = [tProbLeft_all tProbLeft];
        data_stim_tc_all = cat(3,data_stim_tc_all,data_stim_tc);
        data_dec_tc_all = cat(3,data_dec_tc_all,data_dec_tc);
    end
    save(fullfile(rc.structOutput, [mouse '_concatTC.mat']), 'dateStr_all', 'SIx_all', 'MIx_all', 'FIx_all', 'NIx_all', 'tLeftTrial_all', 'tGratingContrast_all', 'tProbLeft_all', 'data_stim_tc_all', 'data_dec_tc_all');
end
        
        
        