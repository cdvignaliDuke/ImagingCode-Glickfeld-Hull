P = 2;
matrix = 'SF5xTF5';
inj = 'V1';
sum_base = 'G:\users\lindsey\analysisLG\experiments';
anal_base = '\\zoloto\bigstorlab\Lindsey\Analysis\120205';
mouse = 'Y13';
date = '110509';
userun = [1:2];

base = 'G:\users\lindsey\analysisLG\active mice';
outDir = fullfile(base, mouse,date);

fn_local =  fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_local_max.mat']);
load(fn_local);

fn_stim = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_allstim.mat']);
load(fn_stim);
FOV = mean(stim_off,3);

fn_df =  fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_stack_dF_np.tif']);
stack = readtiff(fn_df);
dF_max = max(stack,[],3);

figure;
subplot(2,2,1)
imagesq(FOV);
MAX = max(max(FOV,[],1),[],2);
caxis([0 MAX])
title('Field of view (F)')
subplot(4,4,3)
imagesq(stack(:,:,1))
MAX = max(max(dF_max,[],1),[],2);
caxis([0 .2*MAX])
title('1Hz; 0.32 cpd')
subplot(4,4,4)
imagesq(stack(:,:,5));
caxis([0 .2*MAX])
title('15Hz; 0.32 cpd')
subplot(4,4,7)
imagesq(stack(:,:,21));
caxis([0 .2*MAX])
title('1Hz; 0.02 cpd')
subplot(4,4,8)
imagesq(stack(:,:,25));
caxis([0 .2*MAX])
title('15Hz; 0.02 cpd')
subplot(2,2,3)
imagesq(dF_max)
caxis([0 .5*MAX])
title('Maximum dF/F')
colormap(gray)
subplot(2,2,4)
imagesq(dF_max);
caxis([0 .5*MAX])
colormap(gray)
hold on
scatter(j,i,2.5,'b');
title('Significantly responsive boutons')


fn_out = fullfile(anal_base, [matrix '_' num2str(P) 'P_' inj '_Y13_example_FOV.pdf']);
        print(gcf, '-dpdf', fn_out);


fn_summary = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'all_fits.mat');
load(fn_summary);


iArea = 1;
iexp =1;
P = 2;
matrix = 'SF5xTF5';
inj = 'V1';
base = 'G:\users\lindsey\analysisLG\active mice';
sum_base = 'G:\users\lindsey\analysisLG\experiments';
anal_base = '\\zoloto\bigstorlab\Lindsey\Analysis\120205';
mouse = 'Y13';
date = '110509';
userun = [1:2];
outDir = fullfile(base, mouse, date);
iCell = 17;
nON = 12;
nOFF = 12;
nCond = 25;
fn_reps= fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_reps.mat']);
load(fn_reps);
fn_stack = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_sorted.tif']);
stack_sorted = readtiff(fn_stack);

pos = all_fits(iArea).expt(iexp).bouton(iCell).pos;
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
TC_base = mean(TC_reps(7:12,:,:),1);
TC_dF =  bsxfun(@minus, TC_reps, TC_base);
for iCond = 1:nCond
    nReps = stim_reps(:,iCond);
    for iRep = 1:nReps
        TC_dFoverF(:,iRep,iCond) = TC_dF(:,iRep,iCond)./TC_base(:,iRep,iCond);
    end
end

TC_dFoverF_avg = squeeze(mean(TC_dFoverF,2));
TC_dFoverF_sem = squeeze(std(TC_dFoverF,[],2)./sqrt(size(TC_dFoverF,2)));
figure;
site = [1 2 3 4 5 13 14 15 16 17 25 26 27 28 29 37 38 39 40 41 49 50 51 52 53];
for iCond = 1:nCond
    subplot(12,12,site(1, iCond))
    shadedErrorBar(9/24:9/24:9, TC_dFoverF_avg(:,iCond),TC_dFoverF_sem(:,iCond),'k');
    set(gca, 'FontSize', 6);
    ylim([-.5 2])
    xlim([0 10])
    box off
    if iCond>1
        axis off
    end
end

ex_cells = [17 54 143];

start=19;
ax_col = strvcat('r', 'g', 'c');
for iCell = 1:3
    subplot(6,6,start)
    imagesq(all_fits(iArea).expt(iexp).bouton(ex_cells(iCell)).dFoverF);
    set(gca, 'Xcolor',ax_col(iCell,:),'Ycolor',ax_col(iCell,:), 'Ticklength', [0 0], 'LineWidth', 2)
    if iCell == 1;
        title('data')
    end
    subplot(6,6,start+1)
    imagesq(all_fits(iArea).expt(iexp).bouton(ex_cells(iCell)).plotfit);
    set(gca, 'Ticklength', [0 0])
    if iCell == 1;
        title('fit')
    end
    colormap(gray)
    start = start+6;
end
subplot(2,2,2)
imagesq(dF_max);
caxis([0 .5*MAX])
colormap(gray)
hold on
n = all_fits(iArea).expt(iexp).n(1);
for iCell = 1:n
    if all_fits(iArea).expt(iexp).bouton(iCell).goodfit == 1 
        subplot(2,2,2)
        scatter(all_fits(iArea).expt(iexp).bouton(iCell).pos(2),all_fits(iArea).expt(iexp).bouton(iCell).pos(1) ,4,'b');
        hold on
        subplot(2,2,4)
        scatter(all_fits(iArea).expt(iexp).bouton(iCell).TF_fit, all_fits(iArea).expt(iexp).bouton(iCell).SF_fit,4,'b');
        hold on
    else
        subplot(2,2,2)
        scatter(all_fits(iArea).expt(iexp).bouton(iCell).pos(2),all_fits(iArea).expt(iexp).bouton(iCell).pos(1) ,4,[0.3 0.3 0.3]);
    end
end
subplot(2,2,2)
for iCell = 1:3
    scatter(all_fits(iArea).expt(iexp).bouton(ex_cells(iCell)).pos(2), all_fits(iArea).expt(iexp).bouton(ex_cells(iCell)).pos(1),4,ax_col(iCell,:));
end
subplot(2,2,4)
xlim([0 16])
ylim([0 0.35])
xlabel('TF(Hz)')
ylabel('SF(cpd)')
for iCell = 1:3
    scatter(all_fits(iArea).expt(iexp).bouton(ex_cells(iCell)).TF_fit, all_fits(iArea).expt(iexp).bouton(ex_cells(iCell)).SF_fit,4,ax_col(iCell,:));
end

fn_out = fullfile(anal_base, [matrix '_' num2str(P) 'P_' inj '_Y13_example_scatter.pdf']);
        print(gcf, '-dpdf', fn_out);