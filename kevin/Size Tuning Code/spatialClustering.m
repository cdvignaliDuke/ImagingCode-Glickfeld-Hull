%% spatial analysis of FOVs


% test code
x =1:300;
y=1:400;

[xx yy] = meshgrid(x,y);
a = zeros(size(xx));

ind = (xx-69).^2 + (yy-250).^2 < 17^2;
imagesc(ind)

xbar = mean(xx(ind(:)));
ybar = mean(yy(ind(:)));
hold on
plot(xbar,ybar,'x')

%% load csv to define exp

clear all; clc;
%expfile = '\\CRASH.dhe.duke.edu\data\home\kevin\Code\Ai9x_experiment_list.txt';
%expfile = 'C:\Users\kevin\Documents\Repositories\ImagingCode-Glickfeld-Hull\kevin\Size Tuning Code\Ai9x_experiment_list.txt';
expfile = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\kevin\Code\Ai9x_experiment_list.txt';
fID = fopen(expfile);
head = textscan(fID,'%s%s%s%s%s',1,'delimiter',',');
head = vertcat(head{:});
temp = textscan(fID,'%s%s%s%s%s','delimiter',',','HeaderLines',1);
temp = horzcat(temp{:});
expdata = cell2table(temp,'VariableNames',head);
nExp = size(expdata,1);
%isvalid = ones(1,nExp);
%expdata = addvars(expdata,isvalid);

fprintf(['Size tuning spatial analysis - by KM, Glickfeld Lab\nLoading ' num2str(nExp) ' experiments\n'])

%%
iExp = 39;
distComb3_all = []; distComb11_all = []; distComb12_all = []; distComb22_all = []; distCombSame_all = [];
areaTag_all =[];areaTag11_all =[];areaTag12_all =[];areaTag22_all =[];areaTagS_all =[];
nM_all=[];nComb_all=[];
for iExp = 1:nExp
    date = expdata.date{iExp};
    mouse = expdata.mouse{iExp};
    ret_str = 'runs-002';
    run_str = expdata.run_str{iExp};
    fprintf('Experiment %d/%d: mouse %s date %s\n',iExp,nExp,mouse,date)
    
    % load masks
    load(fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ret_str], [date '_' mouse '_' ret_str '_mask_cell.mat']))
    fprintf('Loaded: cell masks')
    
    %load goodfit_ind (ret)
    filename = fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ret_str], [date '_' mouse '_' ret_str '_lbub_fits.mat']);
    if ~exist(filename, 'file')
        fprintf([[date '_' mouse '_' ret_str '_lbub_fits.mat'] ' not found! Please remove from list\n'])
    end
    load(filename, 'lbub_fits', 'goodfit_ind')
    nCellsR = size(lbub_fits,1);
    lbub_fits_ret = lbub_fits;
    cellAz = lbub_fits(goodfit_ind,4,4);
    cellEl = lbub_fits(goodfit_ind,5,4);
    sigmax = lbub_fits(goodfit_ind,2,4);
    sigmay = lbub_fits(goodfit_ind,3,4);
    RFsize = sqrt(2*log(2))*geo_mean([sigmax sigmay],2);
    fprintf(', lbub_fits(ret)')
    
    % load goodfit_ind_size, sizeFits
    % sizeFitResults_SP
    filename = fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_sizeFitResults_SP.mat']);
    if ~exist(filename, 'file')
        fprintf([[date '_' mouse '_' run_str '_sizeFitResults_SP.mat'] ' not found! Please remove from list\n'])
    end
    load(filename, 'sizeFits')
    fprintf(', sizeFitResults_SP')
    
    %cellDists
    filename = fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_sizeTuneData.mat']);
    if ~exist(filename, 'file')
        fprintf([[date '_' mouse '_' run_str '_sizeTuneData.mat'] ' not found! Please remove from list\n'])
    end
    load(filename, 'sizeTune', 'cellDists')
    fprintf(', cellDists')
    
    % lbub_fits (size)
    filename = fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_lbub_fits.mat']);
    if ~exist(filename, 'file')
        fprintf([[date '_' mouse '_' run_str '_lbub_fits.mat'] ' not found! Please remove from list\n'])
    end
    load(filename, 'lbub_fits', 'goodfit_ind_size')
    fprintf(', lbub_fits(size) - done\n')
    
    switch expdata.area{iExp}
        case 'V1'
            RFcutoff = 10;
        case {'LM', 'AL'}
            RFcutoff = 15;
        case 'PM'
            RFcutoff = 20;
    end
    
    y=1:size(mask_cell,1);
    x=1:size(mask_cell,2);
    [xx yy] = meshgrid(x,y);
    
    %figure(1);clf;
    nCell1 = max(mask_cell(:));
    nCell2 = length(goodfit_ind);
    retMask = mask_cell*0;
    for i=1:nCell2
        retMask(mask_cell==goodfit_ind(i)) = i;
    end
    nCell3 = length(goodfit_ind_size);
    
    if size(sizeFits,2)==6
        ism2 = [sizeFits(goodfit_ind_size,5).Ftest];
    else
        ism2 = [sizeFits(goodfit_ind_size,4).Ftest];
    end
    
    % shuffle model identities
    ism2 = ism2(randperm(length(ism2)));
    
    sizeMask = mask_cell*0;modelMap = sizeMask;
    for i=1:nCell3
        sizeMask(retMask==goodfit_ind_size(i)) = i;
%         if ism2(i)
%             modelMap(retMask==goodfit_ind_size(i))=2;
%         else
%             modelMap(retMask==goodfit_ind_size(i))=1;
%         end
    end
    %subplot(2,3,[1 2 3])
    %imagesc(modelMap)
    %colormap(gca,[[0 0 0];[1 0 0];[0 0 1]])
    %title(['Exp: ' num2str(iExp) ', area: ' char(expdata.area{iExp})])
    %hold on
    xbar = 0*[1:nCell3]; ybar = xbar;
    for iCell=1:nCell3
        xb = mean(xx(sizeMask(:)==iCell));
        yb = mean(yy(sizeMask(:)==iCell));
        %plot(xb,yb,'xk')
        xbar(iCell) = xb;
        ybar(iCell) = yb;
    end
    nComb = nchoosek(nCell3,2);
    distComb3 = 0*[1:nComb]; distComb11=[]; distComb12 = []; distComb22 = []; distCombSame = [];
    areaTag=[]; areaTag11=[]; areaTag12=[]; areaTag22=[]; areaTagS=[];
    count = 1;
    ar = find(strcmp(expdata.area{iExp},areas));
    for iCell = 1:nCell3
        if cellDists(goodfit_ind_size(iCell))>RFcutoff
            continue
        end
        for jCell = iCell+1:nCell3
            if cellDists(goodfit_ind_size(jCell))>RFcutoff
                continue
            end
            dist = norm([xbar(iCell) - xbar(jCell), ybar(iCell) - ybar(jCell)]);
            distComb3(count) = dist;
            areaTag = [areaTag ar];
            count=count+1;
            if ism2(iCell) & ism2(jCell)
                distComb22 = [distComb22 dist];
                distCombSame = [distCombSame dist];
                areaTag22 = [areaTag22 ar];
                areaTagS = [areaTagS ar];
            elseif ~ism2(iCell) & ~ism2(jCell)
                distComb11 = [distComb11 dist];
                distCombSame = [distCombSame dist];
                areaTag11 = [areaTag11 ar];
                areaTagS = [areaTagS ar];
            else
                distComb12 = [distComb12 dist];
                areaTag12 = [areaTag12 ar];
            end
        end
    end
%     
%     subplot(2,3,4)
%     histogram(distComb3,50)
%     subplot(2,3,5)
%     histogram(distCombSame,50,'FaceColor','k','facealpha',0.1)
%     hold on
%     histogram(distComb11,50,'FaceColor','g','facealpha',0.3)
%     histogram(distComb22,50,'FaceColor','r','facealpha',0.3)
%     subplot(2,3,6)
%     histogram(distCombSame,50,'FaceColor','g')
%     hold on
%     histogram(distComb12,50,'FaceColor','r')
    
    distComb3_all = [distComb3_all distComb3];
    distComb11_all = [distComb11_all distComb11];
    distComb12_all = [distComb12_all distComb12];
    distComb22_all = [distComb22_all distComb22];
    distCombSame_all = [distCombSame_all distCombSame];
    areaTag_all = [areaTag_all areaTag];
    areaTag11_all = [areaTag11_all areaTag11];
    areaTag12_all = [areaTag12_all areaTag12];
    areaTag22_all = [areaTag22_all areaTag22];
    areaTagS_all = [areaTagS_all areaTagS];
    
    nM_all = [nM_all [sum(~ism2); sum(ism2)]];
    nComb_all = [nComb_all nComb];
    fprintf('Exp: %d done\n',iExp)
    %pause
end
%% simulate random points for expected distance
est = [];
for iTr = 1:100000;
    est(iTr) = sqrt((randi(length(x),1)-randi(length(x),1))^2 + (randi(length(y),1)-randi(length(y),1))^2);
end
mean(est)

%% plot all exps
scaleFactor = geo_mean([581/size(retMask,1),1030/size(retMask,2)])./1.7; %um per pixel
figure(2);clf;
% subplot(1,3,1)
% h = cdfplot(distComb3_all);
% hold on
% line([mean(distComb3_all),mean(distComb3_all)],[0 1],'color','b','LineWidth',2)
% text(500, .2,['all: ' num2str(mean(distComb3_all))])
% text(500, .1,['n=' num2str(length(distComb3_all))])
% title('All cell pairs, all exps')
% xlabel('Pairwise dist (microns)')
% grid off;set(gca,'box','off','TickDir','out')
subplot(1,2,1)
h = cdfplot(distCombSame_all*scaleFactor);
hold on
cdfplot(distComb12_all*scaleFactor);
line([mean(distCombSame_all)*scaleFactor,mean(distCombSame_all)*scaleFactor],[0 1],'color','b','LineWidth',2)
line([mean(distComb12_all)*scaleFactor,mean(distComb12_all)*scaleFactor],[0 1],'color','r','LineStyle','--','LineWidth',2)
title('Same vs Diff pairs, all exps')
xlabel('Pairwise dist (microns)')
text(300, .5,['same: ' num2str(mean(distCombSame_all)*scaleFactor)])
text(300, .4,['n=' num2str(length(distCombSame_all))])
text(300, .2,['diff: ' num2str(mean(distComb12_all)*scaleFactor)])
text(300, .1,['n=' num2str(length(distComb12_all))])
legend('same','diff')
grid off;set(gca,'box','off','TickDir','out')
subplot(1,2,2)
h = cdfplot(distComb11_all*scaleFactor);
hold on
cdfplot(distComb22_all*scaleFactor);
line([mean(distComb11_all)*scaleFactor,mean(distComb11_all)*scaleFactor],[0 1],'color','b','LineWidth',2)
line([mean(distComb22_all)*scaleFactor,mean(distComb22_all)*scaleFactor],[0 1],'color','r','LineStyle','--','LineWidth',2)
title('Same only pairs by model, all exps')
xlabel('Pairwise dist (microns)')
legend('m1-m1','m2-m2')
text(300, .5,['m1-m1: ' num2str(mean(distComb11_all)*scaleFactor)])
text(300, .4,['n=' num2str(length(distComb11_all))])
text(300, .2,['m2-m2: ' num2str(mean(distComb22_all)*scaleFactor)])
text(300, .1,['n=' num2str(length(distComb22_all))])
grid off;set(gca,'box','off','TickDir','out')
%%
figure(3);clf;
ism2_all = [sizeFits_all(:,4).Ftest];
for iA=1:4
    switch iA
        case 1
            RFcutoff = 10;
        case {2 3}
            RFcutoff = 15;
        case 4
            RFcutoff = 20;
    end
    ind = intersect(find(areaInd(goodfit_ind_all)==iA),goodfit_ind_size_all);
    ind = intersect(ind,find(cellDists_all<=RFcutoff));
    xMax=20;
    subplot(4,3,3*iA-2);
    cdfplot(cellDists_all(ind));
    xlim([0 xMax])
    title(['RF-stim distance: Area ' char(areas(iA))])
    text(xMax/2,0.5,['Mean: ' num2str(mean(cellDists_all(ind)))])
    text(xMax/2,0.25,['n=' num2str(length(ind))])
    ind1 = intersect(ind,find(~ism2_all));
    subplot(4,3,3*iA-1);
    cdfplot(cellDists_all(ind1));
    xlim([0 xMax])
    title(['M1 only in ' char(areas(iA))])
    text(xMax/2,0.5,['Mean: ' num2str(mean(cellDists_all(ind1)))])
    text(xMax/2,0.25,['n=' num2str(length(ind1))])
    ind2 = intersect(ind,find(ism2_all));
    subplot(4,3,3*iA);
    cdfplot(cellDists_all(ind2));
    xlim([0 xMax])
    title(['M2 only in ' char(areas(iA))])
    text(xMax/2,0.5,['Mean: ' num2str(mean(cellDists_all(ind2)))])
    text(xMax/2,0.25,['n=' num2str(length(ind2))])
end
%%
figure(4);clf;
for iA=1:length(areas)
    ind = find(areaTag_all == iA);
    ind11 = find(areaTag11_all == iA);
    ind12 = find(areaTag12_all == iA);
    ind22 = find(areaTag22_all == iA);
    indS = find(areaTagS_all == iA);
    
    subplot(4,2,2*iA-1)
    h = cdfplot(distCombSame_all(indS)*scaleFactor);
    hold on
    cdfplot(distComb12_all(ind12)*scaleFactor);
    line([mean(distCombSame_all(indS))*scaleFactor,mean(distCombSame_all(indS))*scaleFactor],[0 1],'color','b','LineWidth',2)
    line([mean(distComb12_all(ind12))*scaleFactor,mean(distComb12_all(ind12))*scaleFactor],[0 1],'color','r','LineStyle','--','LineWidth',2)
    title(['Same vs Diff pairs, area: ' char(areas(iA))])
    xlabel('Pairwise dist (microns)')
    text(300, .5,['same: ' num2str(mean(distCombSame_all(indS))*scaleFactor)])
    text(300, .4,['n=' num2str(length(distCombSame_all(indS)))])
    text(300, .2,['diff: ' num2str(mean(distComb12_all(ind12))*scaleFactor)])
    text(300, .1,['n=' num2str(length(distComb12_all(ind12)))])
    legend('same','diff')
    grid off;set(gca,'box','off','TickDir','out')
    subplot(4,2,2*iA)
    h = cdfplot(distComb11_all(ind11)*scaleFactor);
    hold on
    cdfplot(distComb22_all(ind22)*scaleFactor);
    line([mean(distComb11_all(ind11))*scaleFactor,mean(distComb11_all(ind11))*scaleFactor],[0 1],'color','b','LineWidth',2)
    line([mean(distComb22_all(ind22))*scaleFactor,mean(distComb22_all(ind22))*scaleFactor],[0 1],'color','r','LineStyle','--','LineWidth',2)
    title(['Same only pairs by model, area: ' char(areas(iA))])
    xlabel('Pairwise dist (microns)')
    legend('m1-m1','m2-m2')
    text(300, .5,['m1-m1: ' num2str(mean(distComb11_all(ind11))*scaleFactor)])
    text(300, .4,['n=' num2str(length(distComb11_all(ind11)))])
    text(300, .2,['m2-m2: ' num2str(mean(distComb22_all(ind22))*scaleFactor)])
    text(300, .1,['n=' num2str(length(distComb22_all(ind22)))])
    grid off;set(gca,'box','off','TickDir','out')
end

%% number of m1 cells by area
% start with just pm

expIndi = find(cellfun(@(x) strcmp(x,'PM'), expdata.area, 'UniformOutput', 1));
fracM1 = nM_all(1,:)./sum(nM_all,1);

figure(1);clf;
subplot(1,2,1)
histogram(fracM1(expIndi),10);
xlim([0 1])
xlabel('Frac M1 cells')
subplot(1,2,2)
plot(fracM1(expIndi),nM_all(1,expIndi),'k.')
xlim([0 1])
xlabel('Frac M1 cells')
ylabel('Number M1 cells')

%% stats
[h p] = lillietest(distCombSame_all); % all fail
ranksum(distComb11_all,distComb22_all);

% estimate distance in um
geo_mean([(581/512),1030/796])*mean(distComb12_all)

geo_mean([(581/512),1030/796])*std(distComb12_all)

%% example FOVs
exFOV = [19 34 45 49 3 8 20 42 17 26 35 39 4 18 50 51];
figure(5);clf;
for iFOV = 1:length(exFOV)
    iExp = exFOV(iFOV);
    date = expdata.date{iExp};
    mouse = expdata.mouse{iExp};
    ret_str = 'runs-002';
    run_str = expdata.run_str{iExp};
    fprintf('Experiment %d/%d: mouse %s date %s\n',iExp,nExp,mouse,date)
    
    % load masks
    load(fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ret_str], [date '_' mouse '_' ret_str '_mask_cell.mat']))
    fprintf('Loaded: cell masks')
    
    %load goodfit_ind (ret)
    filename = fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ret_str], [date '_' mouse '_' ret_str '_lbub_fits.mat']);
    if ~exist(filename, 'file')
        fprintf([[date '_' mouse '_' ret_str '_lbub_fits.mat'] ' not found! Please remove from list\n'])
    end
    load(filename, 'lbub_fits', 'goodfit_ind')
    nCellsR = size(lbub_fits,1);
    lbub_fits_ret = lbub_fits;
    cellAz = lbub_fits(goodfit_ind,4,4);
    cellEl = lbub_fits(goodfit_ind,5,4);
    sigmax = lbub_fits(goodfit_ind,2,4);
    sigmay = lbub_fits(goodfit_ind,3,4);
    RFsize = sqrt(2*log(2))*geo_mean([sigmax sigmay],2);
    fprintf(', lbub_fits(ret)')
    
    % load goodfit_ind_size, sizeFits
    % sizeFitResults_SP
    filename = fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_sizeFitResults_SP.mat']);
    if ~exist(filename, 'file')
        fprintf([[date '_' mouse '_' run_str '_sizeFitResults_SP.mat'] ' not found! Please remove from list\n'])
    end
    load(filename, 'sizeFits')
    fprintf(', sizeFitResults_SP')
    
    %cellDists
    filename = fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_sizeTuneData.mat']);
    if ~exist(filename, 'file')
        fprintf([[date '_' mouse '_' run_str '_sizeTuneData.mat'] ' not found! Please remove from list\n'])
    end
    load(filename, 'sizeTune', 'cellDists')
    fprintf(', cellDists')
    
    % lbub_fits (size)
    filename = fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_lbub_fits.mat']);
    if ~exist(filename, 'file')
        fprintf([[date '_' mouse '_' run_str '_lbub_fits.mat'] ' not found! Please remove from list\n'])
    end
    load(filename, 'lbub_fits', 'goodfit_ind_size')
    fprintf(', lbub_fits(size) - done\n')
    
    switch expdata.area{iExp}
        case 'V1'
            RFcutoff = 10;
        case {'LM', 'AL'}
            RFcutoff = 15;
        case 'PM'
            RFcutoff = 20;
    end
    
    nCell1 = max(mask_cell(:));
    nCell2 = length(goodfit_ind);
    retMask = mask_cell*0;
    for i=1:nCell2
        retMask(mask_cell==goodfit_ind(i)) = i;
    end
    nCell3 = length(goodfit_ind_size);
    
    if size(sizeFits,2)==6
        ism2 = [sizeFits(goodfit_ind_size,5).Ftest];
    else
        ism2 = [sizeFits(goodfit_ind_size,4).Ftest];
    end
    
    sizeMask = mask_cell*0;modelMap = sizeMask;
    for i=1:nCell3
        sizeMask(retMask==goodfit_ind_size(i)) = i;
        if ism2(i)
            modelMap(retMask==goodfit_ind_size(i))=2;
            if cellDists(goodfit_ind_size(i))>RFcutoff
                modelMap(retMask==goodfit_ind_size(i))=0;
            end
        else
            modelMap(retMask==goodfit_ind_size(i))=1;
            if cellDists(goodfit_ind_size(i))>RFcutoff
                modelMap(retMask==goodfit_ind_size(i))=0;
            end
        end
    end
    subplot(4,4,iFOV)
    imagesc(modelMap)
    colormap(gca,[[0 0 0];[1 0 0];[0 0 1]]) %;[1 0 1];[0 1 0]
    title(['Exp: ' num2str(iExp) ', area: ' char(expdata.area{iExp})])


end