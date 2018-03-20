mouse_mat = strvcat('i671','i698','i699');
date_mat{1} = strvcat('170127','170131');
date_mat{2} = strvcat('171001','171002','171003');
date_mat{3} = strvcat('171013');
run_mat{1} = {'002-004','003-006'};
run_mat{2} = {'002', '002-003','001-002'};
run_mat{3} = {'001-002'};
dir_mat{1} = {'005','007'};
dir_mat{2} = {'003','004','003'};
dir_mat{3} = {'003'};
time_mat{1} = {'1253','1550'};
time_mat{2} = {'1339','1444','1758'};
time_mat{3} = {'1139'};

%%
p = 1;
for imouse = 2:size(mouse_mat,1)
    mouse = mouse_mat(imouse,:);
    for iexp = 1:size(date_mat{imouse},1)
        date = date_mat{imouse}(iexp,:);
        run = run_mat{imouse}{iexp};
        load(fullfile('Z:\home\lindsey\Analysis\2P',[date '_' mouse],[date '_' mouse '_runs-' run],[date '_' mouse '_runs-' run '_reg_shifts.mat']));
        load(fullfile('Z:\home\lindsey\Analysis\2P',[date '_' mouse],[date '_' mouse '_runs-' run],[date '_' mouse '_runs-' run '_mask_cell.mat']));
        dir_run = dir_mat{imouse}{iexp};
        
        if mouse=='i699'
            fName = ['\\CRASH.dhe.duke.edu\data\home\andrew\Behavior\Data\data-i698-' date '-' time_mat{imouse}{iexp} '.mat'];
        else
        	fName = ['\\CRASH.dhe.duke.edu\data\home\andrew\Behavior\Data\data-' mouse '-' date '-' time_mat{imouse}{iexp} '.mat'];
        end
        load(fName);
        if ~exist(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_runs-' dir_run]))
            mkdir(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_runs-' dir_run]))
        end
        save(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_runs-' dir_run], [date '_' mouse '_runs-' dir_run '_input.mat']), 'input')
        
        if mouse == 'i671'
            CD = ['Z:\home\ashley\data\AW71\two-photon imaging\' date '\' dir_mat{imouse}{iexp}];
        else
            CD = ['Z:\home\lindsey\Data\2P_images\' mouse '\' date '_' mouse '\' dir_mat{imouse}{iexp}];
        end
        cd(CD);
        imgMatFile = [dir_mat{imouse}{iexp} '_000_000.mat'];
        load(imgMatFile);
        
        nframes = [input.counterValues{end}(end) info.config.frames];
        data = sbxread(imgMatFile(1,1:11),0,min(nframes));
        if size(data,1)== 2
            data = data(1,:,:,:);
        end
        data = squeeze(data);
        [out, data_reg] = stackRegister(data,data_avg);
        save(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_runs-' dir_run], [date '_' mouse '_runs-' dir_run '_reg_shifts.mat']), 'out', 'data_avg')
        
        load(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_runs-' run], [date '_' mouse '_runs-' run '_mask_cell.mat']))
        save(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_runs-' dir_run], [date '_' mouse '_runs-' dir_run '_mask_cell.mat']),'mask_cell','mask_np')
        clear data
        
        % neuropil subtraction
        down = 5;
        sz = size(data_reg);

        data_tc = stackGetTimeCourses(data_reg, mask_cell);
        data_reg_down  = stackGroupProject(data_reg,down);
        data_tc_down = stackGetTimeCourses(data_reg_down, mask_cell);
        nCells = size(data_tc,2);
        np_tc = zeros(sz(3),nCells);
        np_tc_down = zeros(floor(sz(3)./down), nCells);
        for i = 1:nCells
             np_tc(:,i) = stackGetTimeCourses(data_reg,mask_np(:,:,i));
             np_tc_down(:,i) = stackGetTimeCourses(data_reg_down,mask_np(:,:,i));
             fprintf(['Cell #' num2str(i) '%s/n']) 
        end
        %get weights by maximizing skew
        ii= 0.01:0.01:1;
        x = zeros(length(ii), nCells);
        for i = 1:100
            x(i,:) = skewness(data_tc_down-tcRemoveDC(np_tc_down*ii(i)));
        end
        [max_skew ind] =  max(x,[],1);
        np_w = 0.01*ind;
        npSub_tc = data_tc-bsxfun(@times,tcRemoveDC(np_tc),np_w);
        clear data_reg data_reg_down

        save(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_runs-' dir_run], [date '_' mouse '_runs-' dir_run '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc')
        clear data_tc data_tc_down np_tc np_tc_down mask_np mask_cell
        
        ori = celleqel2mat_padded(input.tGratingDirectionDeg);
        oris = unique(ori);
        nori = length(oris);
        ntrials = length(ori);
        [nFrames nCells] = size(npSub_tc);
        
        if isfield(input,'doLever')
            cStimOn = celleqel2mat_padded(input.cFirstStim);
            nOn = input.nFramesOn;
            if iscell(nOn)
                nOn = unique(celleqel2mat_padded(nOn));
            end
            nOff = nOn*4;
        else
            nOn = input.nScansOn;
            nOff = input.nScansOff;
            cStimOn = nOff+1:nOff+nOn:nFrames;
        end
        
        pre_win = 10;
        data_trial = nan(pre_win+nOn+nOff,nCells,ntrials);
        for itrial = 1:ntrials
            if nFrames > cStimOn(itrial)+nOn+nOff
                data_trial(:,:,itrial) = npSub_tc(cStimOn(itrial)-pre_win+1:cStimOn(itrial)+nOn+nOff,:);
            end
        end
        data_f = mean(data_trial(1:pre_win,:,:),1);
        data_df = bsxfun(@minus,data_trial,data_f);
        data_dfof = bsxfun(@rdivide,data_df,data_f);
        
        save(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_runs-' dir_run], [date '_' mouse '_runs-' dir_run '_dfofData.mat']), 'data_dfof','pre_win')
        save(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_runs-' dir_run], [date '_' mouse '_runs-' dir_run '_stimData.mat']), 'ori','oris','nori','nOn','ntrials')

        subplot(2,3,p)
        plot(nanmean(data_dfof,3))
        vline([pre_win pre_win+6 pre_win+nOn pre_win+nOn+6])
        p = p+1;
    end
end

%%
for imouse = 1:size(mouse_mat,1)
    mouse = mouse_mat(imouse,:);
    for iexp = 1:size(date_mat{imouse},1)
        date = date_mat{imouse}(iexp,:);
        dir_run = dir_mat{imouse}{iexp};
        run = run_mat{imouse}{iexp};
        load(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_runs-' dir_run], [date '_' mouse '_runs-' dir_run '_dfofData.mat']))
        load(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_runs-' dir_run], [date '_' mouse '_runs-' dir_run '_stimData.mat']))
        nCells = size(data_dfof,2);
        [n n2] = subplotn(nCells);
        ind0 = [find(ori==0) find(ori==180)];
        ind90 = [find(ori==90) find(ori==270)];
        figure;
        for iCell = 1:nCells
            subplot(n,n2,iCell)
            plot(nanmean(data_dfof(:,iCell,ind0),3))
            hold on
            plot(nanmean(data_dfof(:,iCell,ind90),3))
        end
        suptitle([mouse ' ' date '- tuning'])
        
        if mouse == 'i671'
            load(fullfile('Z:\home\lindsey\Analysis\2P',[date '_' mouse],[date '_' mouse '_runs-' run],[date '_' mouse '_runs-' run '_trialData.mat']));
            load(fullfile('Z:\home\lindsey\Analysis\2P',[date '_' mouse],[date '_' mouse '_runs-' run],[date '_' mouse '_runs-' run '_stimData.mat']));
            load(fullfile('Z:\home\lindsey\Analysis\2P',[date '_' mouse],[date '_' mouse '_runs-' run],[date '_' mouse '_runs-' run '_input.mat']));
        else
            load(fullfile('Z:\home\lindsey\Analysis\2P',[date '_' mouse],[date '_' mouse '_runs-' run],[date '_' mouse '_runs-' run '_dfofData.mat']));
            load(fullfile('Z:\home\lindsey\Analysis\2P',[date '_' mouse],[date '_' mouse '_runs-' run],[date '_' mouse '_runs-' run '_stimData.mat']));
            load(fullfile('Z:\home\lindsey\Analysis\2P',[date '_' mouse],[date '_' mouse '_runs-' run],[date '_' mouse '_runs-' run '_input.mat']));
        end
        
        targCyc = celleqel2mat_padded(input.nCyclesOn)+1;
        FIx = strcmp(input.trialOutcomeCell,'failure');
        targCyc(FIx) = nan;
        sz = size(data_dfof);
        data_dfof_targ = nan(sz(1),sz(2),sz(4));
        for itrial = 1:sz(4)
            if ~isnan(targCyc(itrial))
                data_dfof_targ(:,:,itrial) = data_dfof(:,:,targCyc(itrial),itrial);
            end
        end
        
        figure;
        ind90 = find(targetDelta ==90);
        for iCell = 1:nCells
            subplot(n,n2,iCell)
            plot(mean(data_dfof(:,iCell,1,:),4)-mean(mean(data_dfof(1:10,iCell,1,:),4),1))
            hold on
            plot(nanmean(data_dfof_targ(:,iCell,ind90),3)-mean(nanmean(data_dfof_targ(1:10,iCell,ind90),3),1))
        end
        suptitle([mouse ' ' date '- behav'])
    end
end
        
%%
for imouse = 1
    mouse = mouse_mat(imouse,:);
    for iexp = 1:size(date_mat{imouse},1)
        date = date_mat{imouse}(iexp,:);
        dir_run = dir_mat{imouse}{iexp};
        load(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_runs-' dir_run], [date '_' mouse '_runs-' dir_run '_dfofData.mat']))
        load(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_runs-' dir_run], [date '_' mouse '_runs-' dir_run '_stimData.mat']))
        
        dir = ori;
        dirs = oris;
        ndir = nori;
        ori(find(dir>179)) = ori(find(dir>179))-180;
        oris = unique(ori);
        nori = length(oris);
        base_win = 1:pre_win;
        resp_win = pre_win+10:pre_win+nOn+10;
        data_dfof_base = squeeze(mean(data_dfof(base_win,:,:),1));
        data_dfof_resp = squeeze(mean(data_dfof(resp_win,:,:),1));
        [h p] = ttest(data_dfof_base,data_dfof_resp,'dim',2,'tail','left');
        nCells = size(data_dfof,2);
        h_ori = zeros(nCells,nori);
        p_ori = zeros(nCells,nori);
        
        for iori = 1:nori
            ind = find(ori==oris(iori));
            [h_ori(:,iori) p_ori(:,iori)] = ttest(data_dfof_base(:,ind),data_dfof_resp(:,ind),'dim',2,'tail','left','alpha',0.05./(nori-1));
        end
        good_ind = unique([find(h); find(sum(h_ori,2))]);
        
        
        nboot = 1000;
        theta_smooth = 0:1:180;
        data_ori = zeros(nCells,nori,nboot+1);
        k_hat = nan(nCells,nboot+1);
        sse = nan(nCells,nboot+1);
        R_square = nan(nCells,nboot+1);
        y_fit = nan(length(theta_smooth),nCells,nboot+1);
        max_ori = nan(nCells,nboot+1);
        for iboot = 1:nboot+1
            fprintf('%d\n', iboot)
            for iori = 1:nori
                ind = find(ori == oris(iori));
                if iboot == 1
                    ind_use = ind;
                else
                    ind_use = ind(randsample(1:length(ind),length(ind),1));
                end
                data_ori(:,iori,iboot) = nanmean(data_dfof_resp(:,ind_use)-data_dfof_base(:,ind_use),2);
            end
            for iCell = 1:length(good_ind)    
                iC = good_ind(iCell);
                [b_hat, k_hat(iC,iboot), R_hat,u_hat,sse(iC,iboot),R_square(iC,iboot)] = miaovonmisesfit_ori(deg2rad(oris),squeeze(data_ori(iC,:,iboot)));
                y_fit(:,iC,iboot) = b_hat+R_hat.*exp(k_hat(iC,iboot).*(cos(2.*(deg2rad(theta_smooth)-u_hat))-1));
                [y_max max_ori(iC,iboot)] = max(y_fit(:,iC,iboot),[],1);
                max_ori(iC,iboot) = max_ori(iC,iboot)-1;
            end
        end
        
        ori_boot_diff = abs(max_ori(:,1)-max_ori(:,2:nboot+1));
        ori_boot_diff(find(ori_boot_diff>90))= 180-ori_boot_diff(find(ori_boot_diff>90));
        ori_boot_sort = sort(ori_boot_diff,2,'ascend');
        theta_90 = ori_boot_sort(:,900);
        
        save(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_runs-' dir_run], [date '_' mouse '_runs-' dir_run '_oriFit.mat']), 'y_fit','max_ori','theta_90','data_ori','nboot','good_ind','oris','base_win','resp_win')
    end
end


        

        

