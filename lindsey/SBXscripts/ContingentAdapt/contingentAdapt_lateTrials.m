%contingent adaptation experiments
mouse_mat = strvcat('i1303','i1103','i1304','i1303','i1304','i1103','i1303');
date_mat = strvcat('190501','190508','190508','190619','190619','190621','190621');
run_mat = strvcat('002','002','002','003','003','002','002');

LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
nexp = size(mouse_mat,1);

for iexp = 1:nexp
    mouse = mouse_mat(iexp,:);
    date = date_mat(iexp,:);
    run_str = catRunName(run_mat(iexp,:), 1);
    fprintf([mouse ' ' date '\n'])
    
    fn = fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]);
    load(fullfile(fn,[date '_' mouse '_' run_str '_TCs.mat']));
    load(fullfile(fn,[date '_' mouse '_' run_str '_dataStim.mat']));
    
    blockSize = find(diff(ind_noadapt)>1);
    
    ind_noadapt_late = [];
    ind_contadapt_late = [];
    ind_asynadapt_late = [];
    ind_noadapt_early = [];
    ind_contadapt_early = [];
    ind_asynadapt_early = [];
    
    ncyc = floor(length(ind_noadapt)./blockSize);
    for i = 1:ncyc
        ind_noadapt_late = [ind_noadapt_late ind_noadapt(1+ (i-1).*blockSize + blockSize./2:blockSize + (i-1).*blockSize)];
        ind_noadapt_early = [ind_noadapt_early ind_noadapt((i-1).*blockSize + 1:blockSize./2 + (i-1).*blockSize)];
    end
    ncyc = floor(length(ind_asynadapt)./blockSize);
    for i = 1:ncyc
        ind_asynadapt_late = [ind_asynadapt_late ind_asynadapt(1+ (i-1).*blockSize + blockSize./2:blockSize + (i-1).*blockSize)];
        ind_asynadapt_early = [ind_asynadapt_early ind_asynadapt((i-1).*blockSize + 1:blockSize./2 + (i-1).*blockSize)];
    end
    ncyc = floor(length(ind_contadapt)./blockSize);
    for i = 1:ncyc
        ind_contadapt_late = [ind_contadapt_late ind_contadapt(1+ (i-1).*blockSize + blockSize./2:blockSize + (i-1).*blockSize)];
        ind_contadapt_early = [ind_contadapt_early ind_contadapt((i-1).*blockSize + 1:blockSize./2 + (i-1).*blockSize)];
    end
    
    prewin_frames = 15;
    postwin_frames = 45;
    data_test = nan(prewin_frames+postwin_frames,nCells,nTrials);

    for itrial = 1:nTrials
        if cTest(itrial) + postwin_frames < sz(3)
            data_test(:,:,itrial) = npSub_tc(cTest(itrial)-prewin_frames:cTest(itrial)+postwin_frames-1,:);
        end
    end

    data_f = mean(data_test(1:prewin_frames,:,:),1);
    data_dfof_tc = (data_test-data_f)./data_f;

    
    nMask = length(maskCons);
    nTest = length(testCons);
    noadapt_resp_late_cell = cell(nMask,nTest);
    asynadapt_resp_late_cell = cell(nMask,nTest);
    contadapt_resp_late_cell = cell(nMask,nTest);
    noadapt_resp_early_cell = cell(nMask,nTest);
    asynadapt_resp_early_cell = cell(nMask,nTest);
    contadapt_resp_early_cell = cell(nMask,nTest);
    for im = 1:nMask
        ind_mask = find(maskCon == maskCons(im));
        for it = 1:nTest
            ind_test = find(testCon == testCons(it));
            ind = intersect(ind_noadapt_late,intersect(ind_test,ind_mask));
            if length(ind)>1
                noadapt_resp_late_cell{im,it} = squeeze(mean(data_dfof_tc(prewin_frames+6:prewin_frames+21,:,ind),1));
            elseif length(ind)>0
                noadapt_resp_late_cell{im,it} = squeeze(mean(data_dfof_tc(prewin_frames+6:prewin_frames+21,:,ind),1))';
            else
                noadapt_resp_late_cell{im,it} = nan(nCells,1);
            end
            ind = intersect(ind_noadapt_early,intersect(ind_test,ind_mask));
            if length(ind)>1
                noadapt_resp_early_cell{im,it} = squeeze(mean(data_dfof_tc(prewin_frames+6:prewin_frames+21,:,ind),1));
            elseif length(ind)>0
                noadapt_resp_early_cell{im,it} = squeeze(mean(data_dfof_tc(prewin_frames+6:prewin_frames+21,:,ind),1))';
            else
                noadapt_resp_early_cell{im,it} = nan(nCells,1);
            end
            ind = intersect(ind_asynadapt_late,intersect(ind_test,ind_mask));
            if length(ind)>1
                asynadapt_resp_late_cell{im,it} = squeeze(mean(data_dfof_tc(prewin_frames+6:prewin_frames+21,:,ind),1));
            elseif length(ind)>0
                asynadapt_resp_late_cell{im,it} = squeeze(mean(data_dfof_tc(prewin_frames+6:prewin_frames+21,:,ind),1))';
            else
                asynadapt_resp_late_cell{im,it} = nan(nCells,1);
            end
            ind = intersect(ind_asynadapt_early,intersect(ind_test,ind_mask));
            if length(ind)>1
                asynadapt_resp_early_cell{im,it} = squeeze(mean(data_dfof_tc(prewin_frames+6:prewin_frames+21,:,ind),1));
            elseif length(ind)>0
                asynadapt_resp_early_cell{im,it} = squeeze(mean(data_dfof_tc(prewin_frames+6:prewin_frames+21,:,ind),1))';
            else
                asynadapt_resp_early_cell{im,it} = nan(nCells,1);
            end
            ind = intersect(ind_contadapt_late,intersect(ind_test,ind_mask));
            if length(ind)>1
                contadapt_resp_late_cell{im,it} = squeeze(mean(data_dfof_tc(prewin_frames+6:prewin_frames+21,:,ind),1));
            elseif length(ind)>0
                contadapt_resp_late_cell{im,it} = squeeze(mean(data_dfof_tc(prewin_frames+6:prewin_frames+21,:,ind),1))';
            else
                contadapt_resp_late_cell{im,it} = nan(nCells,1);
            end
            ind = intersect(ind_contadapt_early,intersect(ind_test,ind_mask));
            if length(ind)>1
                contadapt_resp_early_cell{im,it} = squeeze(mean(data_dfof_tc(prewin_frames+6:prewin_frames+21,:,ind),1));
            elseif length(ind)>0
                contadapt_resp_early_cell{im,it} = squeeze(mean(data_dfof_tc(prewin_frames+6:prewin_frames+21,:,ind),1))';
            else
                contadapt_resp_early_cell{im,it} = nan(nCells,1);
            end
        end
    end
    save(fullfile(fn,[date '_' mouse '_' run_str '_allCellRespLate.mat']),'noadapt_resp_late_cell','asynadapt_resp_late_cell','contadapt_resp_late_cell','noadapt_resp_early_cell','asynadapt_resp_early_cell','contadapt_resp_early_cell');
end