PPWFdatasets
for iexp = 1:size(expt,2)
    date = expt(iexp).date;
    mouse = ['i' expt(iexp).SubNum];
    disp([date ' ' mouse])
    time_mat = expt(iexp).time_mat;
    nrun = expt(iexp).runs;
    g_ind = expt(iexp).g_ind;
    frame_rate = expt(iexp).frame_rate;
    if expt(iexp).obj > 10
        area_list = expt(iexp).img_loc;
    else
        area_list = strvcat('LM', 'AL', 'V1', 'RL', 'PM');
    end
    nROI = size(area_list,1);
    
    col_mat = strvcat('k', 'b');

    fn_base = fullfile(anal_root, [date '_' mouse]);
    fn_tc = fullfile(fn_base, [date '_' mouse '_' expt_name '.mat']);
    fn_mask = fullfile(fn_base, [date '_' mouse '.mat']);
    
    if exist(fn_tc)
        load(fn_tc)
        load(fn_mask)
    else
        data =[];
        for irun = 1:nrun;
            fn = fullfile(data_root, [date '_' mouse], [date '_' mouse '_' expt_name '_' num2str(irun)], [date '_' mouse '_' expt_name '_' num2str(irun) '_MMStack.ome.tif']);
            data_temp = readtiff(fn);
            data = cat(3,data, data_temp);
            clear data_temp
        end
        sz = size(data);
        tframes = 1:sz(3);
        gframes = g_ind:g_ind:sz(3);
        bframes = setdiff(tframes, gframes);
        bdata = double(data(:,:,bframes));
        gdata = double(data(:,:,gframes));
        bdata= reshape(bdata,[sz(1)*sz(2) length(bframes)]);
        bdata_int = zeros(sz(1)*sz(2), sz(3));
        for i = 1:sz(1)*sz(2)
            bdata_int(i,:) = interp1(bframes,bdata(i,:),tframes,'previous');
        end
        gdata= reshape(gdata,[sz(1)*sz(2) length(gframes)]);
        gdata_int = zeros(sz(1)*sz(2), sz(3));
        for i = 1:sz(1)*sz(2)
            gdata_int(i,:) = interp1(gframes,gdata(i,:),tframes,'previous');
        end
        bdata_int = reshape(bdata_int, [sz(1) sz(2) sz(3)]);
        gdata_int = reshape(gdata_int, [sz(1) sz(2) sz(3)]);

        %load mworks .mat file
        fn_mworks = fullfile('\\CRASH.dhe.duke.edu\data\home\andrew\Behavior\Data', ['data-' mouse '-' date '-' time_mat '.mat']);
        load(fn_mworks)
        stimOff  = cell2mat(input.tStimOffTimeMs);
        ntrials = size(stimOff,2);
        cLeverDown = cell2mat(input.cLeverDown);
        cTargetOn = cell2mat(input.cTargetOn);
        holdTimes = cell2mat(input.holdTimesMs);
        preframes = 15;
        postframes = 45;
        bdata_base = NaN(sz(1), sz(2), preframes+postframes, ntrials);
        gdata_base = NaN(sz(1), sz(2), preframes+postframes, ntrials);
        bdata_test = NaN(sz(1), sz(2), preframes+postframes, ntrials);
        gdata_test = NaN(sz(1), sz(2), preframes+postframes, ntrials);
        for itrial = 1:ntrials
            if holdTimes(itrial) < 10000;
                if ((cLeverDown(itrial)-preframes)>=1) & ((cLeverDown(itrial)+postframes)<=sz(3))
                    bdata_base(:,:,:,itrial) = bdata_int(:,:,cLeverDown(itrial)-preframes:cLeverDown(itrial)+postframes-1);
                    gdata_base(:,:,:,itrial) = gdata_int(:,:,cLeverDown(itrial)-preframes:cLeverDown(itrial)+postframes-1);
                end
                if ((cTargetOn(itrial)-preframes)>=1) & ((cTargetOn(itrial)+postframes)<=sz(3))
                    bdata_test(:,:,:,itrial) = bdata_int(:,:,cTargetOn(itrial)-preframes:cTargetOn(itrial)+postframes-1);
                    gdata_test(:,:,:,itrial) = gdata_int(:,:,cTargetOn(itrial)-preframes:cTargetOn(itrial)+postframes-1);
                end
            end
        end

        bdata_base_F = mean(bdata_base(:,:,1:preframes,:),3);
        bdata_base_dF = bsxfun(@minus, bdata_base, bdata_base_F);
        bdata_base_dFoverF = bsxfun(@rdivide, bdata_base_dF, bdata_base_F);
        gdata_base_F = mean(gdata_base(:,:,1:preframes,:),3);
        gdata_base_dF = bsxfun(@minus, gdata_base, gdata_base_F);
        gdata_base_dFoverF = bsxfun(@rdivide, gdata_base_dF, gdata_base_F);

        bdata_test_dF = bsxfun(@minus, bdata_test, bdata_base_F);
        bdata_test_dFoverF = bsxfun(@rdivide, bdata_test_dF, bdata_base_F);
        gdata_test_dF = bsxfun(@minus, gdata_test, gdata_base_F);
        gdata_test_dFoverF = bsxfun(@rdivide, gdata_test_dF, gdata_base_F);

        clear data bdata gdata bdata_int gdata_int bdata_trials gdata_trials bdata_trials_F bdata_trials_dF gdata_trials_F gdata_trials_dF    

        bdata_base_gsub = bdata_base_dFoverF-gdata_base_dFoverF;
        bdata_test_gsub = bdata_test_dFoverF-gdata_test_dFoverF;
        
        if isdir(fn_base)
            if exist(fullfile(fn_base, [date '_' mouse '.mat']), 'file')
                load(fullfile(fn_base, [date '_' mouse '.mat']));
            else
                roiPolySelect
                save(fullfile(fn_base, [date '_' mouse '.mat']), 'mask_cell', 'area_list', 'bdata_avg')
            end
        else
            mkdir(fn_base)     
            roiPolySelect
            save(fullfile(fn_base, [date '_' mouse '.mat']), 'mask_cell', 'area_list', 'bdata_avg')
        end

        base_tc_sub = reshape(stackGetTimeCourses(reshape(bdata_base_gsub,[sz(1) sz(2) (preframes+postframes)*ntrials]),mask_cell), [(preframes+postframes) ntrials size(area_list,1)]);
        test_tc_sub = reshape(stackGetTimeCourses(reshape(bdata_test_gsub,[sz(1) sz(2) (preframes+postframes)*ntrials]),mask_cell), [(preframes+postframes) ntrials size(area_list,1)]);

        offs = unique(stimOff);
        test_tc_sub_off = cell(1,length(offs));
        base_tc_sub_off = cell(1,length(offs));
        resp_win = 23:26;
        base_win = 17:20;
        
        for ia = 1:size(area_list,1)
            for io = 1:length(offs)
                ind1 = find(stimOff == offs(io));
                test_tc_sub_off{io} = nanmean(test_tc_sub(resp_win,ind1,:),1)-nanmean(test_tc_sub(base_win,ind1,:),1);
                base_tc_sub_off{io} = nanmean(base_tc_sub(resp_win,ind1,:),1)-nanmean(base_tc_sub(base_win,ind1,:),1);
            end
        end
                
        save(fn_tc, 'base_tc_sub', 'test_tc_sub', 'test_tc_sub_off', 'base_tc_sub_off', 'offs', 'stimOff', 'base_win', 'resp_win')
    end
    
    base_all = zeros(1,length(offs));
    for io = 1:length(offs)
        base_all(1,io) = nanmean(base_tc_sub_off{io},2).*1.25;
    end
    base_max = max(base_all,[],2);
    
    for ia = 1:size(area_list,1)
        figure;
        for io = 1:length(offs)
            ind1 = find(stimOff == offs(io));
            subplot(2,3,io)
            shadedErrorBar(1:size(base_tc_sub,1), nanmean(base_tc_sub(:,ind1,ia),2), nanstd(base_tc_sub(:,ind1,ia),[],2)./sqrt(sum(~isnan(squeeze(base_tc_sub(1,ind1,1))),2)),'-k');
            hold on
            shadedErrorBar(1:size(test_tc_sub,1), nanmean(test_tc_sub(:,ind1,ia),2)-nanmean(nanmean(test_tc_sub(base_win,ind1,ia),2),1), nanstd(test_tc_sub(:,ind1,1),[],2)./sqrt(sum(~isnan(squeeze(test_tc_sub(1,ind1,1))),2)),'-r');
            title(num2str(offs(io)))
            ylim([-0.01 base_max])
        end
        suptitle([mouse ' ' date '- Area ' area_list(ia,:) '- Paired pulse TCs'])   
        print(fullfile(fn_base, [date '_' mouse '_' area_list(ia,:) '_' expt_name '_TCs.pdf']), '-dpdf');
    end

    for ia = 1:size(area_list,1)
        figure;
        for io = 1:length(offs)
            norm = bsxfun(@rdivide, test_tc_sub_off{io}, nanmean(base_tc_sub_off{io},2));
            errorbar(offs(io), nanmean(norm(:,:,ia),2), nanstd(norm(:,:,ia),[],2)./sqrt(size(norm,2)), 'ok')
            hold on
        end
        ylim([0 1.5])
        xlabel('Time (sec)')
        ylabel('Normalized amplitude')
        suptitle([mouse ' ' date '- Area ' area_list(ia,:) '- PP recovery'])   
        print(fullfile(fn_base, [date '_' mouse '_' area_list(ia,:) '_' expt_name '_change_resp.pdf']), '-dpdf');
    end
end
 
    
    
   
    