%load data
data = readtiff('Z:\home\miaomiao\other\TA_demo\data\170905_i720_contrast_orientation_1_MMStack_Pos0.ome.tif','single');
load('\\CRASH.dhe.duke.edu\data\home\andrew\Behavior\Data\data-i720-170905-1449.mat')

%extract stimulus variables
rightGratingContrast = celleqel2mat_padded(input.rightGratingContrast);
leftGratingContrast = celleqel2mat_padded(input.leftGratingContrast);
nTrials = length(rightGratingContrast);
cons = unique(rightGratingContrast);
ncon = length(cons);

%extract stimulus timing
cStimOn = celleqel2mat_padded(input.cStimOn);

%extract trial epochs
sz = size(data);
data_mat = zeros(sz(1),sz(2),20,nTrials);
for itrial = 1:nTrials
    data_mat(:,:,:,itrial) = data(:,:,cStimOn(itrial)-10:cStimOn(itrial)+9);
end

%transform into dF/F
data_f = mean(data_mat(:,:,1:9,:),3);
data_df = bsxfun(@minus, data_mat, data_f);
data_dfof = bsxfun(@rdivide, data_df, data_f);

clear data_mat data_f data_df

%measure response to high contrast stim
highRConInd = find(rightGratingContrast == cons(end));
noRConInd = find(rightGratingContrast == cons(1));
highRConResp = mean(data_dfof(:,:,:,highRConInd),4);
noRConResp = mean(data_dfof(:,:,:,noRConInd),4);
highConResp = highRConResp-noRConResp;
writetiff(highConResp,'C:\Users\lindsey\Desktop\highCon.tif');

%extract timecourses
highConAvg = mean(highConResp(:,:,13:17),3);
bwout = imCellEditInteractive(highConAvg);
mask = bwlabel(bwout);
figure; imagesc(mask);
area_list = strvcat('V1','LM','LI','PM');
data_tc = stackGetTimeCourses(data,mask);
clear data_mat data

%extract trial epochs
data_tc_mat = zeros(40,size(data_tc,2),nTrials);
for itrial = 1:nTrials
    data_tc_mat(:,:,itrial) = data_tc(cStimOn(itrial)-10:cStimOn(itrial)+29,:);
end

%transform into dF/F
data_tc_f = mean(data_tc_mat(1:9,:,:),1);
data_tc_df = bsxfun(@minus,data_tc_mat,data_tc_f);
data_tc_dfof = bsxfun(@rdivide,data_tc_df,data_tc_f);

%% contrast only
%find contrast response on right when no stim on left

sz = size(data_tc_mat);
data_tc_con = zeros(sz(1),sz(2),ncon,2);
data_tc_con_avg = zeros(sz(2),ncon,2);
data_tc_con_sem = zeros(sz(2),ncon,2);
n = zeros(2,ncon);
tt = [-9:30]*100;
for iarea = 1:4
    figure;
    a = [0 1];
    for i = 1:2
        indL = find(leftGratingContrast == a(i));
        for icon = 1:ncon
            ind = intersect(indL,find(rightGratingContrast == cons(icon)));
            n(i,icon) = length(ind);
            data_tc_con(:,:,icon,i) = mean(data_tc_dfof(:,:,ind),3);
            subplot(2,2,1+a(i))
            plot(tt,data_tc_con(:,iarea,icon,i))
            hold on
            data_tc_con_avg(:,icon,i) = mean(data_tc_con(13:20,:,icon,i),1);
            data_tc_con_sem(:,icon,i) = std(mean(data_tc_dfof(13:20,:,ind),1),[],3)./sqrt(length(ind));
        end
        ylim([-0.03 0.2])
        title(['Left contrast = ' num2str(a(i))])
        subplot(2,2,3+a(i))
        errorbar(cons, data_tc_con_avg(iarea,:,i),data_tc_con_sem(iarea,:,i))
        xlim([-0.1 1.1])
        ylim([-0.02 0.2])
    end
    suptitle(area_list(iarea,:))
end

%find responsed based on contrast differences
diffGratingContrast = rightGratingContrast-leftGratingContrast;
diffs = unique(diffGratingContrast);
ndiff = length(diffs);

data_tc_diff = zeros(sz(1),sz(2),ndiff);
data_tc_diff_avg = zeros(sz(2),ndiff);
data_tc_diff_sem = zeros(sz(2),ndiff);
for iarea = 1:4
    figure;
    for idiff = 1:ndiff
        ind = find(diffGratingContrast == diffs(idiff));
        data_tc_diff(:,:,idiff) = mean(data_tc_dfof(:,:,ind),3);
        subplot(1,2,1)
        plot(tt,data_tc_diff(:,iarea,idiff))
        hold on
        data_tc_diff_avg(:,idiff) = mean(data_tc_diff(10:14,:,idiff),1);
        data_tc_diff_sem(:,idiff) = std(mean(data_tc_dfof(10:14,:,ind),1),[],3)./sqrt(length(ind));
    end
        ylim([-0.03 0.02])
        subplot(1,2,2)
        errorbar(diffs, data_tc_diff_avg(iarea,:),data_tc_diff_sem(iarea,:))
        ylim([-0.02 0.02])
    suptitle(area_list(iarea,:))
end

%% contrast plus orientation

leftGratingOrientation = celleqel2mat_padded(input.leftGratingDirectionDeg);
oris = unique(leftGratingOrientation);
sz = size(data_tc_mat);
ncon_I = 2;
nori_I = length(oris);
data_tc_ori = zeros(sz(1),sz(2),ncon,ncon_I,nori_I);


data_tc_ori_avg = zeros(sz(2),ncon,ncon_I,nori_I);
data_tc_ori_sem = zeros(sz(2),ncon,ncon_I,nori_I);

n = zeros(ncon, ncon_I,nori_I);
tt = [-9:30]*100;
%for iarea = 1:4
    figure;
    for icon = 1:ncon
        ind = find(rightGratingContrast == cons(icon));
        ii = 1;
        for iIcon = 1:ncon_I
            ind2 = intersect(ind, find(leftGratingContrast == iIcon-1));
            for iori = 1:nori_I
                if iIcon == 1
                    ind3 = ind2;
                    if iori == 2
                        continue
                    end
                else
                    ind3 = intersect(ind2, find(leftGratingOrientation == oris(iori)));
                end
                n(icon,iIcon,iori) = length(ind3);
                data_tc_ori(:,:,icon,iIcon,iori) = mean(data_tc_dfof(:,:,ind3),3);
                data_tc_ori_avg(:,icon,iIcon,iori) = squeeze(mean(mean(data_tc_dfof(13:30,:,ind3),1),3));
                data_tc_ori_sem(:,icon,iIcon,iori) = squeeze(std(mean(data_tc_dfof(13:30,:,ind3),1),[],3)./sqrt(length(ind3)));  
                subplot(3,ncon,icon)
                plot(tt, data_tc_ori(:,1,icon,iIcon,iori));
                xlim([tt(1) tt(end)])
                ylim([-0.02 0.1])
                title(num2str(cons(icon)))
                if icon == 1
                    legend({'noIpsi', 'vertical', 'horizontal'})
                end
                hold on
                subplot(3,ncon,icon+3)
                errorbar(ii, data_tc_ori_avg(1,icon,iIcon,iori),data_tc_ori_sem(1,icon,iIcon,iori),'-o')
                hold on
                xlim([0 4])
                ylim([-0.02 0.1])
                ii = ii +1;
            end
            
        end
    end
    
    str = {'noIpsi', 'vertical', 'horizontal'};
    start = 7;
    ii = 1;
    for iIcon = 1:ncon_I
        for iori = 1:nori_I
            if iIcon == 1 & iori == 2
                continue
            else
                subplot(3,3,start)
                errorbar(cons, data_tc_ori_avg(1,:,iIcon,iori), data_tc_ori_sem(1,:,iIcon,iori),'-o')
                start = start+1;
            end
            ylim([-0.02 0.1])
            xlim([-0.2 1.2])
            xlabel('Contrast')
            title(str{ii})
            ii = ii+1;
        end
    end
    
    figure;
    for iIcon = 1:ncon_I
        for iori = 1:nori_I
            if iIcon == 1 & iori == 2
                continue
            else
                errorbar(cons, data_tc_ori_avg(1,:,iIcon,iori), data_tc_ori_sem(1,:,iIcon,iori),'-o')
                hold on;
            end
            ylim([-0.02 0.1])
            xlim([-0.2 1.2])
            xlabel('Contrast')
        end
    end
    text(0,0.08,['p-noIpsiVsHor: ' num2str(p_noIpsiVsHor)])
    text(0,0.07,['p-noIpsiVsVert: ' num2str(p_noIpsiVsVert)])
    text(0,0.06,['p-HorVsVert: ' num2str(p_HorVsVert)])
    
    iarea = 1;
    ind_con = find(rightGratingContrast == cons(1));
    iIcon = 1;
    ind_noIpsi = intersect(ind_con, find(leftGratingContrast == iIcon-1));
    iIcon = 2;
    ind_Ipsi = intersect(ind_con, find(leftGratingContrast == iIcon-1));
    iori = 1;
    ind_Vert = intersect(ind_Ipsi, find(leftGratingOrientation == oris(1)));
    ind_Hor = intersect(ind_Ipsi, find(leftGratingOrientation == oris(2)));
    
    [h_noIpsiVsVert p_noIpsiVsVert] = ttest2(mean(data_tc_dfof(13:30,1,ind_noIpsi),1),mean(data_tc_dfof(13:30,1,ind_Vert),1));
        [h_noIpsiVsHor p_noIpsiVsHor] = ttest2(mean(data_tc_dfof(13:30,1,ind_noIpsi),1),mean(data_tc_dfof(13:30,1,ind_Hor),1));
        [h_HorVsVert p_HorVsVert] = ttest2(mean(data_tc_dfof(13:30,1,ind_Hor),1),mean(data_tc_dfof(13:30,1,ind_Vert),1));

    
        for iIcon = 1:ncon_I
            ind2 = intersect(ind, find(leftGratingContrast == iIcon-1));
            for iori = 1:nori_I
                if iIcon == 1
                    ind3 = ind2;
                    if iori == 2
                        continue
                    end
                else
                    ind3 = intersect(ind2, find(leftGratingOrientation == oris(iori)));
                end
                n(icon,iIcon,iori) = length(ind3);
                data_tc_ori(:,:,icon,iIcon,iori) = ttest(mean(data_tc_dfof(13:30,iarea,ind),1);
                data_tc_ori_avg(:,icon,iIcon,iori) = squeeze(mean(mean(data_tc_dfof(13:30,iarea,ind3),1),3));
    %end

    
            
                    
            








    