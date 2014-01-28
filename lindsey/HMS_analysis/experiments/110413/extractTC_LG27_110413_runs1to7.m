%% set path name
setpath;

global USER_PROFILE
USER_PROFILE = 'lindsey'; % 

protdir = dirs.images;
%%
newdir = '110413_LG27\runs1to7';
expt = frGetExpt(newdir);
%% register
frRegister(newdir,'Overwrite',true,'Oversampling',10,...
               'DoRecurse',false,'Engine','subpixel');
%% file properties
stack = single(readtiff(expt.dirs.reggreenpn));
nON = 150;
nOFF = 150;
nCond = 7;
nStim = 4;
epoch = (nON + nOFF).*nStim;
[a b z] = size(stack);

%% get active cells
rep = z/epoch;
down = 100;
dec = zeros(a, b, z/down);
dec = stackGroupProject(stack, down);
all = zeros(a, b, epoch/down);
for time = 1:epoch/down;
    all(:,:,time) = mean(dec(:,:,time:epoch/down:end),3);
end
avg = mean(dec,3);
dFoverF = zeros(a, b, epoch/down);
for time = 1:epoch/down;
    dFoverF(:,:,time) = (all(:,:,time)-avg)./avg;
end
%dF_dec = stackGroupProject(dFoverF, 10);
%max_proj = max(dF_dec, [], 3);
max_proj = max(dFoverF, [], 3);
figure;
imagesc(max_proj);

clear dec
%% get active cells during on period
rep = z/epoch;
stim_matrix = zeros(nON*rep*nStim, 1);
baseline_matrix = zeros((nOFF/3)*rep*nStim,1);

for stim = 1:rep*nStim
    stim_matrix(1+((stim-1)*nON):nON+((stim-1)*nON))= (nOFF+1+((stim-1)*(nON+nOFF)):(nON+nOFF)+((stim-1)*(nON+nOFF)));
    baseline_matrix(1+((stim-1)*(nOFF/3)):(nOFF/3)+((stim-1)*(nOFF/3))) =1+(2*nOFF/3) + ((stim-1)*(nON+nOFF)):nOFF +((stim-1)*(nON+nOFF));
end

stim_all = zeros(a, b, nON*rep*nStim);
stim_all = stack(:,:,stim_matrix);
stim_istim_irep = stackGroupProject(stim_all,nON);

baseline_all = zeros(a, b, (nOFF/3)*rep*nStim);
baseline_all = stack(:,:,baseline_matrix);
baseline_istim_irep = stackGroupProject(baseline_all,nOFF/3);

down = 50;
baseline_avg = mean(baseline_all,3);
stim_down = zeros(a, b, (nON*rep*nStim)/down);
stim_down = stackGroupProject(stim_all,down);

dFoverF_down = zeros(a, b, nON*rep*nStim)/down;
for time = 1:(nON*rep*nStim)/down;
    dFoverF_down(:,:,time) = (stim_down(:,:,time)-baseline_istim_irep(:,:,ceil(time/3)))./baseline_istim_irep(:,:,ceil(time/3));
end
dFoverF_avg = zeros(a, b, nStim*(nON/down));
for stim = 1:nStim*(nON/down);
    dFoverF_avg(:,:,stim) = mean(dFoverF_down(:,:,stim:nStim*(nON/down):end),3);
end

max_proj_on = max(dFoverF_avg,[],3);
figure;
imagesc(max_proj_on);
%% find cell masks
% cell find code from MH
bwimgcell = imCellEditInteractive(max_proj,[]);
mask_cell = bwlabel(bwimgcell);
image(mask_cell);
save(fullfile(expt.dirs.analrootpn,'LG27_runs1to7_masks.mat'),'mask_cell');

%% cell time courses
%reload
mask_cell = load(fullfile(expt.dirs.analrootpn,'LG27_runs1to7_masks.mat'));
mask_cell = mask_cell.mask_cell;
stack = single(readtiff(expt.dirs.reggreenpn));

allcells = max(unique(mask_cell));
stack_dec = stackGroupProject(stack, 10);

%get time courses and remove low frequencies
timeCourses = stackGetTimeCourses(stack_dec,mask_cell);
timeCourses_lowcut = tcLowCut (timeCourses, 1000, 'gaussian', 1);



for stim = 1:rep*nStim
    stim_matrix_down(1+((stim-1)*nON/10):nON/10+((stim-1)*nON/10))= (nOFF/10+1+((stim-1)*(nON+nOFF)/10):(nON+nOFF)/10+((stim-1)*(nON+nOFF)/10));
    baseline_matrix_down(1+((stim-1)*(nOFF/30)):(nOFF/30)+((stim-1)*(nOFF/30))) =1+(2*nOFF/30) + ((stim-1)*(nON+nOFF)/10):nOFF/10 +((stim-1)*(nON+nOFF)/10);
end

%get dF/F for trial
baseline_all = mean(timeCourses_lowcut(baseline_matrix_down,:));
dF_all = bsxfun(@minus,timeCourses_lowcut,baseline_all);
ratio_all = bsxfun(@rdivide,dF_all,baseline_all)*100;

av = tcCycleAverage(ratio_all,epoch/10);
%plot overlay of all trials with average
figure;
for iCell = 1:allcells;
    subplot(ceil(sqrt(allcells)),ceil(sqrt(allcells)),iCell)
    for trial = 1:rep;
        plot(ratio_all(((1+((trial-1)*epoch/10)):epoch/10*trial), iCell), 'c');
        hold on;
    end;
    plot(av(:, iCell),'k');
    axis([0 120 -20 150]);
    xlabel('Time (s)');
    ylabel('dF/F');hold on;
end

%plot overlay of trials by condition
for iCell = 1:allcells;
    figure;
    title(['Cell ' num2str(iCell)])
    for cond = 1:nCond;
        ratio_all_cond = zeros(epoch/10, rep/nCond);
        for trial = 1+((rep/nCond)*(cond-1)):cond*(rep/nCond);
            subplot(3,3,cond);
            plot(ratio_all(((1+((trial-1)*(epoch/10))):epoch*trial/10), iCell), 'c');
            hold on;
            ratio_all_cond(1:(epoch/10),trial-((cond-1)*(rep/nCond))) = ratio_all(((1+((trial-1)*(epoch/10))):epoch*trial/10), iCell);
        end;
        av_ratio_all_cond = mean(ratio_all_cond,2);
        plot(av_ratio_all_cond,'k')
    end;
    subplot(3,3,1);
    title('motor at 10 rpm');
    subplot(3,3,2);
    title('motor at 0 rpm');
    subplot(3,3,3);
    title('motor at 20 rpm'); 
    subplot(3,3,4);
    title('motor at 5 rpm');
    subplot(3,3,5);
    title('freely rotating');
    subplot(3,3,6);
    title('freely rotating; 25 deg');
    subplot(3,3,7);
    title('motor at 10 rpm; 25 deg');
end

%% Quantify amplitude and stdev
mask_cell = load(fullfile(expt.dirs.analrootpn,'LG27_runs1to4_masks.mat'));
mask_cell = mask_cell.mask_cell;
stack = single(readtiff(expt.dirs.reggreenpn));
allcells = max(unique(mask_cell));
timeCourses = stackGetTimeCourses(stack,mask_cell);

%measure dF/F
dFoverF = zeros(rep/nCond,nStim,nCond,allcells);
for iCell = 1:allcells;
    for cond = 1:nCond
        for stim = 1:nStim
            for trial = 1:(rep/nCond);
                 stim_avg = mean(timeCourses(1+nON+((trial-1)*epoch)+((stim-1)*(nOFF+nON))+((cond-1)*epoch*(rep/nCond)):nOFF+nON+((trial-1)*epoch)+((stim-1)*(nOFF+nON))+((cond-1)*epoch*(rep/nCond)),iCell));
                 baseline_avg = mean(timeCourses(1+(2*nOFF/3)+((trial-1)*epoch)+((stim-1)*(nOFF+nON))+((cond-1)*epoch*(rep/nCond)):nOFF+((trial-1)*epoch)+((stim-1)*(nOFF+nON))+((cond-1)*epoch*(rep/nCond)),iCell));
                 dFoverF(trial,stim,cond,iCell)= (stim_avg-baseline_avg)./baseline_avg;
            end
        end
    end
end

dFoverF_avg = mean(dFoverF,1)*100;
dFoverF_stdev = nanstd(dFoverF,[],1)*100;
dFoverF_sem = dFoverF_stdev/sqrt(rep/nCond); 

%plot dF/F
for iCell = 1:allcells;
    figure;
    for iCond = 1:cond;
        subplot(3,3,iCond)
        errorbar([0 45 90 135], dFoverF_avg(:,:,iCond,iCell), dFoverF_sem(:,:,iCond,iCell));
        axis([-15 150 -.5 75]);
        xlabel('Orientation (Deg)');
        ylabel('dF/F');
        hold on;
    end
end

figure;
for iCell = 1:allcells;
    subplot(ceil(sqrt(allcells)),ceil(sqrt(allcells)),iCell)
    errorbar([0 45 90 135], dFoverF_avg(:,:,1,iCell), dFoverF_sem(:,:,1,iCell),'om-','markerfacecolor','m');
    hold on;
    errorbar([0 45 90 135], dFoverF_avg(:,:,2,iCell), dFoverF_sem(:,:,2,iCell),'ok-','markerfacecolor','k');
    hold on;
    errorbar([0 45 90 135], dFoverF_avg(:,:,3,iCell), dFoverF_sem(:,:,3,iCell),'or-','markerfacecolor','r');
    hold on;
    errorbar([0 45 90 135], dFoverF_avg(:,:,4,iCell), dFoverF_sem(:,:,4,iCell),'ob-','markerfacecolor','b');
    hold on;
    errorbar([0 45 90 135], dFoverF_avg(:,:,5,iCell), dFoverF_sem(:,:,5,iCell),'og-','markerfacecolor','g');
    axis([-15 150 -5 75]);
    xlabel('Orientation (Deg)');
    ylabel('dF/F');
end

%% Determine effect of conditions

for iCell = 1:allcells
    peak = max(max(dFoverF_avg(:,:,:,iCell)));
    pos = find(dFoverF_avg(:,:,:,iCell) == peak);
    cond = ceil(pos/nStim);
    stim = ((pos/nStim)-(cond-1))*nStim;
    dFoverF_norm(:,:,:,iCell) = dFoverF_avg(:,stim,:,iCell)/peak;
    dFoverF_var(:,:,:,iCell) = dFoverF_stdev(:,stim,:,iCell);
end

dFoverF_norm_cond = squeeze(dFoverF_norm);
dFoverF_norm_cond_avg = mean(dFoverF_norm_cond,2);
dFoverF_norm_cond_stdev = nanstd(dFoverF_norm_cond,[],2);
dFoverF_norm_cond_sem = dFoverF_norm_cond_stdev/sqrt(allcells);

dFoverF_var_cond = squeeze(dFoverF_var);
dFoverF_var_cond_avg = mean(dFoverF_var_cond,2);
dFoverF_var_cond_stdev = nanstd(dFoverF_var_cond,[],2);
dFoverF_var_cond_sem = dFoverF_var_cond_stdev/sqrt(allcells);

subplot(2,1,1)
errorbar([1:7], dFoverF_norm_cond_avg, dFoverF_norm_cond_sem, 'ok-','markerfacecolor','k');
xlabel('Condition');
ylabel('Normalized dF/F');

subplot(2,1,2)
errorbar([1:7], dFoverF_var_cond_avg, dFoverF_var_cond_sem, 'ok-','markerfacecolor','k');
xlabel('Condition');
ylabel('Standard Deviation');
            
            