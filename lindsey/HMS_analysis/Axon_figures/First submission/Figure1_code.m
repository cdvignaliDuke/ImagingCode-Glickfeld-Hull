fig_base = '\\zmey\storlab\users\Lindsey\Projects\HVAs\Manuscript\Figures_2012';
fig = 1;

col = strvcat('c', 'k', 'r');
P = 2;
matrix = 'SF5xTF5';
inj = 'V1';
sum_base = 'G:\users\lindsey\analysisLG\experiments';

fn_summary = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'all_fits.mat');
load(fn_summary);
fn_good = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'good_fits.mat');
load(fn_good);

fn_FOV = 'G:\users\lindsey\analysisLG\active mice\M14\111117\111117_M14_run1_FOV.tif';
epi_FOV = double(readtiff(fn_FOV));
epi_FOV_log = log(epi_FOV);
figure;
subplot(1,2,1)
imagesc(epi_FOV);
axis image
subplot(1,2,2)
imagesc(epi_FOV_log);
axis image

mouse = 'M14';
date = '111127';
userun = [1:4];
iArea = 1;
iexp = 12;

base = 'G:\users\lindsey\analysisLG\active mice';
outDir = fullfile(base, mouse,date);

local_max = [];
i = [];
j = [];
fn_local =  fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_local_max.mat']);
load(fn_local);

fn_stim = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_allstim.mat']);
load(fn_stim);
FOV = mean(stim_off,3);

fn_df =  fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_stack_dF_np.tif']);
stack = readtiff(fn_df);

siz = size(stack);

dF_nodir = zeros(siz(1), siz(2), 25);
start = 1;
for iCond = 1:25
    dF_nodir(:,:,iCond) = mean(stack(:,:,start:start+1),3);
    start = start+2;
    subplot(5,5,iCond)
    imagesq(dF_nodir(:,:,iCond));
end
dF_max = max(dF_nodir,[],3);

dF_nodir_gamma = zeros(size(dF_nodir));
for iCond = 1:25
    dF_nodir_gamma(:,:,iCond) = imadjust(dF_nodir(:,:,iCond),[],[],1/2.2);
end

dF_max_gamma = imadjust(dF_max,[],[],1/2.2);
FOV_gamma = imadjust(FOV./max(max(FOV,[],2),[],1),[],[],1/2.2);

ind = [];
n = all_fits(iArea).expt(iexp).n(1);
for iCell = 1:n
    pos = all_fits(iArea).expt(iexp).bouton(iCell).pos;
    if 134<pos(1) & pos(1)<211
        if 99<pos(2) & pos(2)<176
            ind = [ind iCell];
        end
    end
end

figure;
imagesq(dF_max);
colormap(gray);
hold on
n = length(ind);
for iCell = 1:n
    pos = all_fits(iArea).expt(iexp).bouton(ind(iCell)).pos;
    text(pos(2), pos(1), num2str(ind(iCell)),'Color', 'w');
    hold on
end
xlim([100 175])
ylim([135 210])

ex_cells = [284 297 255 321];
ex_col = strvcat('c', 'g', 'r', 'y');

figure;
subplot(2,2,1)
imagesq(FOV);
caxis([0 500])
title('FOV')
colormap(gray)
subplot(2,2,2)
imagesq(dF_nodir_gamma(:,:,11));
caxis([0 .75])
title('dF/F')
colormap(gray)
subplot(4,4,9)
imagesq(dF_nodir_gamma(:,:,1));
title('1Hz; 0.32 cpd')
caxis([0 .75])
subplot(4,4,10)
imagesq(dF_nodir_gamma(:,:,5));
caxis([0 .75])
title('15Hz; 0.32 cpd')
subplot(4,4,13)
imagesq(dF_nodir_gamma(:,:,21));
caxis([0 .75])
title('1Hz; 0.02 cpd')
subplot(4,4,14)
imagesq(dF_nodir_gamma(:,:,25));
caxis([0 .75])
title('15Hz; 0.02 cpd')
subplot(2,2,4)
imagesq(dF_max_gamma);
caxis([0 1])
title('Maximum dF/F')
colormap(gray)
hold on
x = [101 175];
y = [136 210];
x_n = length(x(1):1:x(2));
plot(x(1)*ones(1,x_n), y(1):1:y(2), ':w');
hold on
plot(x(2)*ones(1,x_n), y(1):1:y(2), ':w');
hold on
plot(x(1):1:x(2), y(1)*ones(1,x_n), ':w');
hold on
plot(x(1):1:x(2), y(2)*ones(1,x_n), ':w');
fn_out = fullfile(fig_base, ['Figure' num2str(fig)], [matrix '_' num2str(P) 'P_' inj '_M14_FOV.ps']);
        print(gcf, '-depsc', fn_out);
        
        
figure
subplot(2,2,1)
imagesq(dF_max_gamma);
colormap(gray);
caxis([0 1]);
xlim(x);
ylim(y);
subplot(2,2,2)
imagesq(dF_max_gamma);
colormap(gray);
caxis([0 1]);
hold on
n = length(ind);
xlim(x);
ylim(y);
for iCell = 1:n
    pos = all_fits(iArea).expt(iexp).bouton(ind(iCell)).pos;
    if all_fits(iArea).expt(iexp).bouton(iCell).goodfit == 1;
        scatter(pos(2), pos(1), 100, [0.7 0.7 0.7])
    else
        scatter(pos(2), pos(1), 100, 'k')
    end
    hold on
end
title('Significantly responsive boutons')
for iCell = 1:length(ex_cells)
    pos = all_fits(iArea).expt(iexp).bouton(ex_cells(iCell)).pos;
    scatter(pos(2), pos(1), 100,[ex_col(iCell,:)])
    hold on
end
fn_out = fullfile(fig_base, ['Figure' num2str(fig)], [matrix '_' num2str(P) 'P_' inj '_M14_FOV_boutons.ps']);
        print(gcf, '-depsc2', fn_out);


%timecourses tuning
%resort data with time window starting at pt 7
nON = 12;
nOFF = 12;
nCond = 50;
nPlanes = 1;
P = 2;
begin = 7;
pre_win = [1 6];
post_win = [7 18];

if length(userun) == 1
    seqfile =[date '_' mouse '_run' num2str(userun) '_Seqposition.mat'];
    load(fullfile(outDir,'analysis',seqfile));
    Big_Seqposition = Seqposition;
else
    seqfile = [date '_' mouse '_run' num2str(userun) '_Big_Seqposition.mat'];
    load(fullfile(outDir,'analysis',seqfile));
end

stack = [];
if P ==1
    for iRun = 1:length(userun);
        substack = readtiff(fullfile(base,mouse,date, [date '_' mouse '_run' num2str(userun(iRun)) '.tif']));
        stack = cat(3, stack, substack);
    end
elseif P==2
    if nPlanes ==1
        for iRun = 1:length(userun);
            substack = uint16(readtiff(fullfile(base,mouse,date, [date '_' mouse '_run' num2str(userun(iRun)) '_dec_reg.tif'])));
            stack = cat(3, stack, substack);
        end
    else
        for iRun = 1:length(userun);
            substack = uint16(readtiff(fullfile(base,mouse,date, [date '_' mouse '_run' num2str(userun(iRun)) '_plane' num2str(iPlane) '_reg.tif'])));
            stack = cat(3, stack, substack);
        end
    end
end
clear('substack');

stack_sorted = zeros(size(stack) ,'uint16');

%resort stimuli
start = 1;
for iCond = 1:nCond;
    nRep = length(Big_Seqposition(iCond).ind);
    for iRep = 1:nRep;
        ind = Big_Seqposition(iCond).ind(iRep);
        if begin-1+((nOFF+nON)/nPlanes)+((ind-1)*((nOFF+nON)/nPlanes)) > size(stack_sorted,3);
            stack_sorted(:,:,start:start-1+((nOFF+nON)/nPlanes)-begin+1) = stack(:,:,begin+((ind-1)*((nOFF+nON)/nPlanes)):((nOFF+nON)/nPlanes)+((ind-1)*((nOFF+nON)/nPlanes)));
            stack_sorted(:,:,start+((nOFF+nON)/nPlanes)-begin+1:start+((nOFF+nON)/nPlanes)-1) = stack(:,:,1:begin-1);
        else
        stack_sorted(:,:,start:start-1+((nOFF+nON)/nPlanes)) = stack(:,:,begin+((ind-1)*((nOFF+nON)/nPlanes)):begin-1+((nOFF+nON)/nPlanes)+((ind-1)*((nOFF+nON)/nPlanes)));
        end
        start = start+((nOFF+nON)/nPlanes);
    end
end

%resort blanks
nblanks = length(Big_Seqposition(end).ind);
for iblank = 1:nblanks;
    ind = Big_Seqposition(end).ind(iblank);
        if begin-1+((nOFF+nON)/nPlanes)+((ind-1)*((nOFF+nON)/nPlanes)) > size(stack_sorted,3);
            stack_sorted(:,:,start:start-1+((nOFF+nON)/nPlanes)-begin+1) = stack(:,:,begin+((ind-1)*((nOFF+nON)/nPlanes)):((nOFF+nON)/nPlanes)+((ind-1)*((nOFF+nON)/nPlanes)));
            stack_sorted(:,:,start+((nOFF+nON)/nPlanes)-begin+1:start+((nOFF+nON)/nPlanes)-1) = stack(:,:,1:begin-1);
        else
        stack_sorted(:,:,start:start-1+((nOFF+nON)/nPlanes)) = stack(:,:,begin+((ind-1)*((nOFF+nON)/nPlanes)):begin-1+((nOFF+nON)/nPlanes)+((ind-1)*((nOFF+nON)/nPlanes)));
        end
    start = start+((nOFF+nON)/nPlanes);
end
if start<size(stack_sorted,3);
    stack_sorted(:,:,start:end)=[];
end

fn = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_sorted_new.tif']);
writetiff(stack_sorted,fn);

x = 1:5;
y = 5:-1:1;
[x_grid y_grid] = meshgrid(x,y);
x_grid_long = reshape(x_grid', [25 1]);
y_grid_long = reshape(y_grid', [25 1]);

fn = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_reps.mat']);
load(fn);
for iCell = 2:4
    pos = all_fits(iArea).expt(iexp).bouton(ex_cells(iCell)).pos;
    suby = pos(1)-1:pos(1)+1;
    subx = pos(2)-1:pos(2)+1;
    TC = squeeze(mean(mean(stack_sorted(suby,subx,:),2),1));

    TC_reps = zeros(24, stim_reps(:,1), nCond);
    start = 1;
    for iCond = 1:nCond
        nReps = stim_reps(:,iCond);
        for iRep = 1:nReps
            TC_reps(:,iRep,iCond) = TC(start:start+(nON+nOFF)-1,:);
            start= start+(nON+nOFF);
        end
    end

    TC_dFoverF = zeros(size(TC_reps));
    TC_base = mean(TC_reps(pre_win(1):pre_win(2),:,:),1);
    TC_dF =  bsxfun(@minus, TC_reps, TC_base);
    for iCond = 1:nCond
        nReps = stim_reps(:,iCond);
        for iRep = 1:nReps
            TC_dFoverF(:,iRep,iCond) = TC_dF(:,iRep,iCond)./TC_base(:,iRep,iCond);
        end
    end
    dFoverF_avg_90 = squeeze(mean(mean(TC_dFoverF(post_win(1):post_win(2),:,1:2:end),1),2));
    dFoverF_avg_270 = squeeze(mean(mean(TC_dFoverF(post_win(1):post_win(2),:,2:2:end),1),2));
    alldir_dFoverF = [dFoverF_avg_90; dFoverF_avg_270];
    
    TC_dFoverF_avg = zeros(nON+nOFF, nCond/2);
    TC_dFoverF_sem = zeros(nON+nOFF, nCond/2);
    start = 1;
    for iCond = 1:nCond/2
        TC_dFoverF_tempA =  squeeze(TC_dFoverF(:,:,start));
        TC_dFoverF_tempB =  squeeze(TC_dFoverF(:,:,start+1));
        TC_dFoverF_dirs = [TC_dFoverF_tempA TC_dFoverF_tempB];
        TC_dFoverF_avg(:,iCond)= mean(TC_dFoverF_dirs,2);
        TC_dFoverF_sem(:,iCond)= std(TC_dFoverF_dirs,[],2)./sqrt(size(TC_dFoverF_dirs,1));
        start = start+2;
    end
    
    figure;
    a = 4:1:8;
    site = [a a+12 a+24 a+36 a+48];
    for iCond = 1:nCond/2
        subplot(12,12,site(1, iCond))
        shadedErrorBar(9/24:9/24:9, TC_dFoverF_avg(:,iCond),TC_dFoverF_sem(:,iCond),'k');
        set(gca, 'FontSize', 6);
        ylim([-.1 0.7])
        xlim([0 10])
        box off
        if iCond>1
            axis off
        end
    end

     fit = reshape(all_fits(iArea).expt(iexp).bouton(ex_cells(iCell)).plotfit',1,25);
     for iCond = 1:25
        subplot(3,5,12)
        scatter(x_grid_long(iCond,1),y_grid_long(iCond,1), ((dFoverF_avg_90(iCond,:)./max(alldir_dFoverF,[],1))*40)^2,'.k')
        hold on
        axis square
        axis off
        xlim([0 6])
        ylim([0 6])
        title('90 deg');
        subplot(3,5,13)
        scatter(x_grid_long(iCond,1),y_grid_long(iCond,1), ((dFoverF_avg_270(iCond,:)./max(alldir_dFoverF,[],1))*40)^2,'.k')
        hold on
        axis square
        axis off
        xlim([0 6])
        ylim([0 6])
        title('270 deg');
        subplot(3,5,14)
        scatter(x_grid_long(iCond,1),y_grid_long(iCond,1), ((fit(:,iCond)./max(alldir_dFoverF,[],1))*40)^2,'.k');
        hold on
        axis square
        axis off
        title('Fit'); 
        xlim([0 6])
        ylim([0 6])
     end
     suptitle(['Cell #' num2str(ex_cells(iCell))]);
     fn_out = fullfile(fig_base, ['Figure' num2str(fig)], [matrix '_' num2str(P) 'P_' inj '_M14_tuning_Cell#' num2str(ex_cells(iCell)) '.ps']);
        print(gcf, '-depsc2', fn_out);
end

