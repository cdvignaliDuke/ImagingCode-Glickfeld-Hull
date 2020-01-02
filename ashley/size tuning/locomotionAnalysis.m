clear all
area_list = {'LM', 'AL', 'PM'};
n_run_sz_all = [];
n_norun_sz_all = [];
n_run_ret_all = [];
n_norun_ret_all = [];
mouse_all = [];
speed_cutoff = 6; 
for iarea = 1:size(area_list,2);
ds = ['szTuning_axons_' area_list{iarea}];
eval(ds)
nexp = size(expt,2);
wheel_ret = zeros(1,nexp);
wheel_sz = zeros(1,nexp);
rc = behavConstsAV;
n_run_sz = nan(nexp,8);
n_norun_sz = nan(nexp,8);
n_run_ret = nan(nexp,7,2);
n_norun_ret = nan(nexp,7,2);
for iexp = 1:nexp
    mouse = expt(iexp).mouse;
    expDate = expt(iexp).date;
    retRun = expt(iexp).retinotopyFolder;
    szRun = expt(iexp).sizeTuningFolder;
    fnout = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate);
    run_resp_all{iarea} = [];
    norun_resp_all{iarea} = [];
        
    if exist(fullfile(fnout, cell2mat(retRun), [mouse '_' expDate '_input.mat']))
        load(fullfile(fnout, cell2mat(retRun), [mouse '_' expDate '_input.mat']))
        if input.doWheelSpeed == 1
            ntrials = length(input.tGratingDiameterDeg);
            for itrial = 1:ntrials
                if sum(input.wheelSpeedValues{itrial},2)>0
                    wheel_ret(1,iexp) = 1;
                    break
                end
            end
        end
    end

    if wheel_ret(1,iexp)

        mouse = expt(iexp).mouse;
        expDate = expt(iexp).date;
        retRun = expt(iexp).retinotopyFolder;
        szRun = expt(iexp).sizeTuningFolder;
        fnout = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate);
        load(fullfile(fnout, cell2mat(retRun), [mouse '_' expDate '_input.mat']))

        ntrials = length(input.tGratingDirectionDeg);
        nOn = input.nScansOn;
        nOff = input.nScansOff;
        
        wheel_speed = zeros(1,ntrials);
        for itrial = 1:ntrials
            fr_num_start = (itrial.*nOff) + ((itrial-1).*nOn);
            fr_num_end = (itrial.*nOff) + (itrial.*nOn);
            fr_time_start = input.counterTimesUs{itrial}(find(input.counterValues{itrial}==fr_num_start));
            fr_time_end = input.counterTimesUs{itrial}(find(input.counterValues{itrial}==fr_num_end));
            ind = intersect(find(input.wheelSpeedTimesUs{itrial}>=fr_time_start),find(input.wheelSpeedTimesUs{itrial}<=fr_time_end));
            wheel_speed(:,itrial) = mean(input.wheelSpeedValues{itrial}(ind),2);
        end

        figure; 
        subplot(3,1,1)
        plot(wheel_speed)
        ylabel('Wheel speed')
        xlabel('Trial')


        run_ind = find(wheel_speed>=speed_cutoff);
        norun_ind = find(wheel_speed<speed_cutoff);

        az_mat = celleqel2mat_padded(input.tGratingAzimuthDeg);
        azs = unique(az_mat);
        nAz = length(azs);
        for iaz = 1:nAz
            n_run_ret(iexp,iaz,1) = length(find(az_mat(run_ind)==azs(iaz)));
            n_norun_ret(iexp,iaz,1) = length(find(az_mat(norun_ind)==azs(iaz)));
        end

        subplot(3,1,2)
        plot(azs,n_run_ret(iexp,:,1)./length(run_ind),'-o')
        hold on
        plot(azs,n_norun_ret(iexp,:,1)./length(norun_ind),'-o')
        xlabel('Azimuth (deg)')
        ylabel('Fract. Trials')
        ylim([0 1])
        legend(['Run- n = ' num2str(length(run_ind))], ['Stationary- n = ' num2str(length(norun_ind))])

        el_mat = celleqel2mat_padded(input.tGratingElevationDeg);
        els = unique(el_mat);
        nEl = length(els);
        for iel = 1:nEl
            n_run_ret(iexp,iel,2) = length(find(el_mat(run_ind)==els(iel)));
            n_norun_ret(iexp,iel,2) = length(find(el_mat(norun_ind)==els(iel)));
        end

        subplot(3,1,3)
        plot(els,n_run_ret(iexp,:,2)./length(run_ind),'-o')
        hold on
        plot(els,n_norun_ret(iexp,:,2)./length(norun_ind),'-o')
        xlabel('Elevation (deg)')
        ylabel('Fract. Trials')
        ylim([0 1])
        legend(['Run- n = ' num2str(length(run_ind))], ['Stationary- n = ' num2str(length(norun_ind))])

        suptitle([mouse ' ' expDate ' Ret'])
        print(fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\SizeTuning\Locomotion', [mouse '_' expDate '_locomotionRet_2cps.pdf']),'-dpdf','-fillpage')
    end
    
    if exist(fullfile(fnout, cell2mat(szRun), [mouse '_' expDate '_input.mat']))
        load(fullfile(fnout, cell2mat(szRun), [mouse '_' expDate '_input.mat']))
        if input.doWheelSpeed == 1
            ntrials = length(input.tGratingDiameterDeg);
            for itrial = 1:ntrials
                if sum(input.wheelSpeedValues{itrial},2)>0
                    wheel_sz(1,iexp) = 1;
                    break
                end
            end
        end
    end
    
    if wheel_sz(1,iexp)
        ntrials = length(input.tGratingDirectionDeg);
        nOn = input.nScansOn;
        nOff = input.nScansOff;
        wheel_speed = zeros(1,ntrials);
        for itrial = 1:ntrials
            fr_num_start = (itrial.*nOff) + ((itrial-1).*nOn);
            fr_num_end = (itrial.*nOff) + (itrial.*nOn);
            fr_time_start = input.counterTimesUs{itrial}(find(input.counterValues{itrial}==fr_num_start));
            fr_time_end = input.counterTimesUs{itrial}(find(input.counterValues{itrial}==fr_num_end));
            ind = intersect(find(input.wheelSpeedTimesUs{itrial}>=fr_time_start),find(input.wheelSpeedTimesUs{itrial}<=fr_time_end));
            wheel_speed(:,itrial) = mean(input.wheelSpeedValues{itrial}(ind),2);
        end

        figure; 
        subplot(2,1,1)
        plot(wheel_speed)
        ylabel('Wheel speed')
        xlabel('Trial')


        run_ind = find(wheel_speed>=speed_cutoff);
        norun_ind = find(wheel_speed<speed_cutoff);
        
        sz_mat = celleqel2mat_padded(input.tGratingDiameterDeg);
        szs = unique(sz_mat);
        nSize = length(szs);
        for isz = 1:nSize
            n_run_sz(iexp,isz) = length(find(sz_mat(run_ind)==szs(isz)));
            n_norun_sz(iexp,isz) = length(find(sz_mat(norun_ind)==szs(isz)));
        end
        
        if (length(run_ind)./length(wheel_speed))>0.1
            load(fullfile(fnout, cell2mat(szRun), [mouse '_' expDate '_Tuning.mat']))
            nCells = size(sizeTune,2);
            run_resp = nan(nCells,nSize);
            norun_resp = nan(nCells,nSize);
            for isz = 1:nSize
                ind = find(sz_mat == szs(isz));
                run_ind_sz = [];
                for i = 1:length(run_ind)
                    run_ind_sz = [run_ind_sz find(ind == run_ind(i))];
                end
                if length(run_ind_sz)>0
                    for iCell = 1:nCells
                        run_resp(iCell,isz) = mean(sizeTune{isz,iCell}(run_ind_sz,:),1);
                    end
                end
                norun_ind_sz = [];
                for i = 1:length(norun_ind)
                    norun_ind_sz = [norun_ind_sz find(ind == norun_ind(i))];
                end
                if length(norun_ind_sz)>0
                    for iCell = 1:nCells
                        norun_resp(iCell,isz) = mean(sizeTune{isz,iCell}(norun_ind_sz,:),1);
                    end
                end
            end
            load(fullfile(fnout, cell2mat(szRun), [mouse '_' expDate '_lbub_fits.mat']))
            goodfit_ind_size = intersect(goodfit_ind_size,find(cellDists<10));
            figure;
            errorbar(szs, nanmean(run_resp(goodfit_ind_size,:),1), nanstd(run_resp(goodfit_ind_size,:)./sqrt(length(goodfit_ind_size)),[],1), '-o')
            hold on
            errorbar(szs, nanmean(norun_resp(goodfit_ind_size,:),1), nanstd(norun_resp(goodfit_ind_size,:)./sqrt(length(goodfit_ind_size)),[],1), '-o')
            legend('run', 'stationary')
            xlabel('Size (deg)')
            ylabel('dF/F')
            title([mouse ' ' expDate])
            print(fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\SizeTuning\Locomotion', [mouse '_' expDate '_locomotionSzResp_2cps.pdf']),'-dpdf','-fillpage')
            
            
            run_resp_all{iarea} = [run_resp_all{iarea}; run_resp];
            norun_resp_all{iarea} = [norun_resp_all{iarea}; norun_resp];
            
        end
            
            
        subplot(2,1,2)
        plot(szs,n_run_sz(iexp,:)./length(run_ind),'-o')
        hold on
        plot(szs,n_norun_sz(iexp,:)./length(norun_ind),'-o')
        xlabel('Size (deg)')
        ylabel('Fract. Trials')
        ylim([0 1])
        xlim([0 100])
        legend(['Run- n = ' num2str(length(run_ind))], ['Stationary- n = ' num2str(length(norun_ind))])
        suptitle([mouse ' ' expDate ' Size'])
        print(fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\SizeTuning\Locomotion', [mouse '_' expDate '_locomotionSz_2cps.pdf']),'-dpdf','-fillpage')
    end

    close all
    mouse_all = [mouse_all; mouse];
end
save(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\SizeTuning\Locomotion\','trialsPerCond_locomotion_' area_list{iarea} '_2cps.mat'],'n_run_sz', 'n_norun_sz', 'n_run_ret', 'n_norun_ret');
n_run_sz_all = [n_run_sz_all; n_run_sz];
n_norun_sz_all = [n_norun_sz_all; n_norun_sz];
n_run_ret_all = [n_run_ret_all; n_run_ret];
n_norun_ret_all = [n_norun_ret_all; n_norun_ret];
end

save(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\SizeTuning\Locomotion\','trialsPerCond_locomotionSummary_2cps.mat'],'n_run_sz_all', 'n_norun_sz_all', 'n_run_ret_all', 'n_norun_ret_all');

fractrun_sz = sum(n_run_sz_all,2)./(sum(n_run_sz_all,2) + sum(n_norun_sz_all,2));
fractrun_sz_avg = nanmean(fractrun_sz,1);
fractrun_sz_sem = nanstd(fractrun_sz,[],1)./sqrt(sum(~isnan(fractrun_sz),1));

fractrun_szs = n_run_sz_all./(n_run_sz_all + n_norun_sz_all);
figure; 
plot(szs', fractrun_szs', '-c')
hold on
errorbar(szs, nanmean(fractrun_szs,1), nanstd(fractrun_szs,[],1)./sqrt(sum(~isnan(fractrun_szs(:,1)),1)),'-ok')
title(['Fraction Running trials for Size- n = ' num2str(sum(~isnan(fractrun_sz(:,1)),1))])
print(fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\SizeTuning\Locomotion', 'locomotionSzSummary_2cps.pdf'),'-dpdf','-fillpage')

        
fractrun_ret = sum(n_run_ret_all(:,:,1),2)./(sum(n_run_ret_all(:,:,1),2) + sum(n_norun_ret_all(:,:,1),2));
fractrun_ret_avg = nanmean(fractrun_ret,1);
fractrun_ret_sem = nanstd(fractrun_ret,[],1)./sqrt(sum(~isnan(fractrun_ret),1));

fractrun_az = n_run_ret_all(:,:,1)./(n_run_ret_all(:,:,1) + n_norun_ret_all(:,:,1));
fractrun_el = n_run_ret_all(:,:,2)./(n_run_ret_all(:,:,2) + n_norun_ret_all(:,:,2));
figure; errorbar(1:7, nanmean(fractrun_az,1), nanstd(fractrun_az,[],1)./sqrt(sum(~isnan(fractrun_az(:,1)),1)),'-o')
hold on
errorbar(1:7, nanmean(fractrun_el,1), nanstd(fractrun_el,[],1)./sqrt(sum(~isnan(fractrun_el(:,1)),1)),'-o')
title(['Fraction Running trials for Retinotopy- n = ' num2str(sum(~isnan(fractrun_el(:,1)),1))])

mouse_all = str2num(mouse_all);
expt_run_ind = find(~isnan(fractrun_sz));
expt_mouse_ind = unique(mouse_all(expt_run_ind,:));
col_mat = strvcat('k', 'b', 'g', 'm');
figure;
for i = 1:size(expt_mouse_ind,1)
    mouse_temp = expt_mouse_ind(i,:);
    ind = [];
    for ii = 1:size(mouse_all,1)
        if find(mouse_all(ii,:)==mouse_temp) & ~isnan(fractrun_sz(ii))
            ind = [ind ii]; 
        end
    end
    if length(ind)>1
        for ii = 1:length(ind)-1
            for iii = ii+1:length(ind)
                temp = [fractrun_sz(ind(ii)) fractrun_sz(ind(iii))];
                scatter(min(temp,[],2),max(temp,[],2), ['o' col_mat(i)])
                hold on
            end
        end
    end
end
xlim([0 0.3])
ylim([0 0.3])
xlabel('Fraction running trials')
ylabel('Fraction running trials')

