mouse_mat = strvcat('i745', 'i923','i842','i843','i843');
date_mat = strvcat('171210', '180205','180402','180406','180416');
run_str_mat = {'runs-002','runs-002-003','runs-002','runs-002','runs-002'};
noff_mat = [3 3 3 3 1];

nexp = size(mouse_mat,1);
mouse_str = [];
for iexp = 1:nexp
    mouse_str = [mouse_str mouse_mat(iexp,:)];
end

LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
%%
for iexp = 1:nexp
    mouse = mouse_mat(iexp,:);
    date = date_mat(iexp,:);
    run_str = run_str_mat{iexp};
    load(fullfile([LG_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']))
    load(fullfile([LG_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dfofData.mat']))
    load(fullfile([LG_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))
    
    if input.doRandCon
        baseCon = nan(maxCyc,nTrials);
        for itrial = 1:nTrials
            baseCon(:,itrial) = input.tBaseGratingContrast{itrial}(1:tCyc(itrial));
        end
    end


    targCon = celleqel2mat_padded(input.tGratingContrast);
    cons = unique(targCon);
    ncon = length(cons);

    ind_base = find(baseCon(1,:)==max(cons));
    [h1,p1] = ttest(squeeze(mean(data_dfof(base_win,:,1,ind_base),1))',squeeze(mean(data_dfof(resp_win,:,1,ind_base),1))','tail','left','alpha',0.05);
    [h,p] = ttest(squeeze(mean(data_dfof(base_win,:,1,:),1))',squeeze(mean(data_dfof(resp_win,:,1,:),1))','tail','left','alpha',0.05);
    good_ind_temp = find(h1+h);

    ht1 = zeros(nDelta,nCells);
    pt1 = zeros(nDelta,nCells);


    targCyc = unique(celleqel2mat_padded(input.tCyclesOn))+1;
    if length(targCyc==1)
        if nDelta>1
            for idelta = 1:nDelta
                ind = find(targetDelta == deltas(idelta));
                for iCell = 1:nCells
                    [ht1(idelta,iCell), pt1(idelta,iCell)] = ttest(squeeze(nanmean(data_dfof(base_win,iCell,targCyc,ind),1)),squeeze(nanmean(data_dfof(resp_win,iCell,targCyc,ind),1)),'tail','left','alpha',0.05./(nDelta-1));
                end
            end
        end
        [ht,pt] = ttest(squeeze(mean(data_dfof(base_win,:,targCyc,:),1))',squeeze(mean(data_dfof(resp_win,:,targCyc,:),1))','tail','left','alpha',0.05);
    end
    good_ind_targ = find(sum([ht;ht1],1));

    good_ind_base = good_ind_temp;
    good_ind = unique([good_ind_base good_ind_targ]);
    good_int = intersect(good_ind_base, good_ind_targ);
    good_sub = find(ismember(good_ind_base,good_int));

    save(fullfile(LG_base,'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respSig.mat']), 'p','h','p1','h1','pt','ht','pt1','ht1','good_ind','good_ind_base','good_ind_targ')

    resp = cell(ncon,noff);
    for icon = 1:ncon
        for ioff = 1:noff
            resp{icon,ioff} = [];
            for itrial = 1:nTrials
                ind = intersect(find(baseCon(3:end,itrial)==cons(icon)),find(tFramesOff(itrial,2:end-1)==offs(ioff)));
                if length(ind)>0
                    resp{icon,ioff} = cat(3,resp{icon,ioff},data_dfof(:,:,ind+2,itrial));
                end
            end
        end
    end
    resp_off = cell(1,noff);
    for ioff = 1:noff
        resp_off{ioff} = [];
        for icon = 1:ncon
            resp_off{ioff} = cat(3,resp_off{ioff},resp{icon,ioff});
        end
    end

    resp_con = cell(1,ncon);
    for icon = 1:icon
        resp_con{icon} = [];
        for ioff = 1:noff
            resp_con{icon} = cat(3,resp_con{icon},resp{icon,ioff});
        end
    end

    % amplitude distributions

    resp_con_amp = cell(ncon,nDelta+1);
    cell_con = zeros(ncon,nDelta+1);
    cell_ori = zeros(ncon,nDelta+1);
    for icon = 1:ncon
        resp_con_amp{icon,1} = squeeze(mean(resp_con{1,icon}(resp_win,:,:),1)-mean(resp_con{1,icon}(base_win,:,:),1));
        cell_con(icon,1) = cons(icon);
        cell_ori(icon,1) = 0;
        for itarg = 1:nDelta
            ind = intersect(find(tGratingDir == deltas(itarg)), find(targCon == cons(icon)));
            resp_con_amp{icon,itarg+1} = squeeze(mean(data_dfof(resp_win,:,end,ind),1)-mean(data_dfof(base_win,:,end,ind),1));
            cell_con(icon,itarg+1) = cons(icon);
            cell_ori(icon,itarg+1) = deltas(itarg);
        end
    end

    save(fullfile(LG_base,'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respDist.mat']),'resp_con_amp','cell_con','cell_ori','offs', 'good_ind_base','good_ind_targ')
end