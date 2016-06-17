FSWFdatasets
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
        gratingDirectionDeg = cell2mat(input.gratingDirectionDeg);
        holdTimes = cell2mat(input.holdTimesMs);
        stimOff  = cell2mat(input.tStimOffTimeMs);
        ntrials = size(gratingDirectionDeg,2);
        cStimOn = cell2mat(input.cFirstStim);
        preframes = 15;
        postframes = 120;
        bdata_trials = NaN(sz(1), sz(2), preframes+postframes, ntrials);
        gdata_trials = NaN(sz(1), sz(2), preframes+postframes, ntrials);
        for itrial = 1:ntrials
            if holdTimes(itrial) < 10000;
                if ((cStimOn(itrial)-preframes)>=1) & ((cStimOn(itrial)+postframes)<=sz(3))
                    bdata_trials(:,:,:,itrial) = bdata_int(:,:,cStimOn(itrial)-preframes:cStimOn(itrial)+postframes-1);
                    gdata_trials(:,:,:,itrial) = gdata_int(:,:,cStimOn(itrial)-preframes:cStimOn(itrial)+postframes-1);
                end
            end
        end

        bdata_trials_F = mean(bdata_trials(:,:,1:preframes,:),3);
        bdata_trials_dF = bsxfun(@minus, bdata_trials, bdata_trials_F);
        bdata_trials_dFoverF = bsxfun(@rdivide, bdata_trials_dF, bdata_trials_F);
        gdata_trials_F = mean(gdata_trials(:,:,1:preframes,:),3);
        gdata_trials_dF = bsxfun(@minus, gdata_trials, gdata_trials_F);
        gdata_trials_dFoverF = bsxfun(@rdivide, gdata_trials_dF, gdata_trials_F);

        clear data bdata gdata bdata_int gdata_int bdata_trials gdata_trials bdata_trials_F bdata_trials_dF gdata_trials_F gdata_trials_dF

        bdata_trials_gsub = bdata_trials_dFoverF-gdata_trials_dFoverF;

    %     writetiff(squeeze(mean(bdata_trials_gsub,4)), 'C:\Users\lindsey\Desktop\160308_all.tif')

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
    
        tc_sub = reshape(stackGetTimeCourses(reshape(bdata_trials_gsub,[sz(1) sz(2) (preframes+postframes)*ntrials]),mask_cell), [(preframes+postframes) ntrials size(area_list,1)]);

        dirs=  unique(gratingDirectionDeg);
        offs = unique(stimOff);
        tc_sub_dir = cell(length(offs),length(dirs),6,size(area_list,1));
        tc_sub_base = cell(length(offs),5,size(area_list,1));
        
        intervals = unique(cell2mat(input.cTargetOn)-cell2mat(input.cStimOn));
        for io = 1:length(intervals)
            resp_win{io} = [(23:intervals(io):(intervals(io)*5)+23)' (26:intervals(io):(intervals(io)*5)+26)'];
            base_win{io} = [(17:intervals(io):(intervals(io)*5)+17)' (20:intervals(io):(intervals(io)*5)+20)'];
        end

        for ia = 1:size(area_list,1)  
            for io = 1:length(offs)
                ind1 = find(stimOff == offs(io));
                base_temp = base_win{io};
                resp_temp = resp_win{io};
                for id = 1:length(dirs)
                    ind = intersect(ind1, find(gratingDirectionDeg == dirs(id)));
                    indn(id) = length(ind);
                    for ip = 1:6
                        tc_sub_dir{io,id,ip,ia} = nanmean(tc_sub(resp_temp(ip,1):resp_temp(ip,2),ind,:),1)-nanmean(tc_sub(base_temp(ip,1):base_temp(ip,2),ind,:),1);
                    end
                end
                for ip = 1:5
                    tc_sub_base{io,ip,ia} = nanmean(tc_sub(resp_temp(ip,1):resp_temp(ip,2),ind1,:),1)-nanmean(tc_sub(base_temp(ip,1):base_temp(ip,2),ind1,:),1);
                end
            end
        end

        save(fullfile(fn_base, [date '_' mouse '_' expt_name '.mat']), 'tc_sub', 'tc_sub_dir', 'tc_sub_base', 'dirs', 'offs', 'stimOff', 'gratingDirectionDeg', 'resp_win', 'base_win')
    end

    [n, n2] = subplotn(size(dirs,2)*size(offs,2));
    for ia = 1:size(area_list,1)  
        figure;
        start=0;
        for io = 1:length(offs)
            ind1 = find(stimOff == offs(io));
            base_temp = base_win{io};
            resp_temp = resp_win{io};
            for id = 1:length(dirs)
                ind = intersect(ind1, find(gratingDirectionDeg == dirs(id)));
                indn(id) = length(ind);
                subplot(n,n2,id+start)
                shadedErrorBar(1:size(tc_sub,1), nanmean(tc_sub(:,ind,ia),2), nanstd(tc_sub(:,ind,ia),[],2)./sqrt(length(ind)));
    %             plot(1:size(tc_sub,1), mean(tc_sub(:,ind,1),2), col_mat(id,:));
                vline(base_temp(:,1),'--k')
                vline(base_temp(:,2),'--k')
                vline(resp_temp(:,1))
                vline(resp_temp(:,2))
                title([num2str(dirs(id)) ' deg; ' num2str(offs(io)) ' ms off'])
                hold on
            end
            start = start+length(dirs);
        end
        suptitle([mouse ' ' date '- Area ' area_list(ia,:)])   
        print(fullfile(fn_base, [date '_' mouse '_FS5_' area_list(ia,:) '.pdf']), '-dpdf');
    end
    
    [n, n2] = subplotn(length(offs));
    for ia = 1:size(area_list,1)  
        figure;
        for io = 1:length(offs)
            ind1 = find(stimOff == offs(io));
            base_temp = base_win{io};
            resp_temp = resp_win{io};
            subplot(n,n2,io)
            for id = 1:length(dirs)
                ind = intersect(ind1, find(gratingDirectionDeg == dirs(id)));
                indn(id,io) = length(ind);
                shadedErrorBar(1:size(tc_sub,1), nanmean(tc_sub(:,ind,ia),2), nanstd(tc_sub(:,ind,ia),[],2)./sqrt(length(ind)), col_mat(id,:));
                hold on;
                title([num2str(offs(io)) ' ms off'])
                hold on
            end
        end
        suptitle([mouse ' ' date '- Area ' area_list(ia,:) '; Black: 45 deg; Blue: 90 deg'])   
        print(fullfile(fn_base, [date '_' mouse '_FS5_' area_list(ia,:) '_overlay.pdf']), '-dpdf');
    end

    figure;
    [n,n2] = subplotn(size(area_list,1));
    for ia = 1:size(area_list,1)
        subplot(n,n2,ia)
        for io = 1:length(offs)
            temp = bsxfun(@rdivide, tc_sub_base{io,5,ia}, nanmean(tc_sub_base{io,1,ia},2));
            shadedErrorBar(dirs,repmat(nanmean(temp(:,:,ia),2), [1, length(dirs)]),repmat(nanstd(temp(:,:,ia),[],2)./sqrt(size(temp,2)), [1, length(dirs)]), col_mat(io,:));
            hold on
            for id = 1:length(dirs)
                temp = bsxfun(@rdivide, tc_sub_dir{io,id,6,ia}, nanmean(tc_sub_base{io,1,ia},2));
                errorbar(dirs(id),nanmean(temp(:,:,ia),2),nanstd(temp(:,:,ia),[],2)./sqrt(size(temp,2)), ['-o' col_mat(io,:)])
                hold on
            end
        end
        title(['Area ' area_list(ia,:)])
        ylim([0 1.5])
    end
    suptitle([mouse ' ' date '- Black: 250ms; Blue: 500ms'])
    print(fullfile(fn_base, [date '_' mouse '_FS5_change_resp.pdf']), '-dpdf');
end
    
% d = [];    
% group_o = [];
% group_d = [];
% ia = 2;
%     for io = 1:2
%         for id = 1:2
%             temp = squeeze(tc_sub_dir{io,id,6});
%             group_o = [group_o;  double(offs(io)).*ones(size(temp,1),1)];
%             group_d = [group_d;  double(dirs(id)).*ones(size(temp,1),1)];
%             d = [d; temp(:,ia)];
%         end
%     end
%     [p, tab, stats] = anovan(d, {group_o,group_d}, 'model', 2, 'varnames', {'Interval', 'Direction'});

            
