%% set path name
setpath;

global USER_PROFILE
USER_PROFILE = 'lindsey'; % 

protdir = dirs.images;
%%
newdir = '110413_LG27\runs5to7';
expt = frGetExpt(newdir);
%%
frRegister(newdir,'Overwrite',true,'Oversampling',10,...
               'DoRecurse',false,'Engine','subpixel');
stack = single(readtiff(expt.dirs.reggreenpn));

%% file properties
nON = 150;
nOFF = 150;
nCond = 4;
epoch = (nON + nOFF).*nCond;
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

%% find cell masks
% cell find code from MH
bwimgcell = imCellEditInteractive(max_proj,[]);
mask_cell = bwlabel(bwimgcell);
image(mask_cell);
save(fullfile(expt.dirs.analrootpn,'LG27_runs1to4_masks.mat'),'mask_cell');

%% cell time courses
mask_cell = load(fullfile(expt.dirs.analrootpn,'LG27_runs1to4_masks.mat'));
allcells = max(unique(mask_cell));
substack = stack(:,:,:);
substack_dec = stackGroupProject(substack, 10);
%get time courses and remove low frequencies
timeCourses = stackGetTimeCourses(substack_dec,mask_cell);
timeCourses_lowcut = tcLowCut (timeCourses, 200, 'gaussian', 1);
av = tcCycleAverage(timeCourses_lowcut,epoch/10);
%get dF/F for trial
baseline_all = mean(timeCourses_lowcut(1.5*nOFF/10):(nOFF+nON)/10:end,:));
dF_all = bsxfun(@minus,timeCourses_lowcut,baseline_all);
ratio_all = bsxfun(@rdivide,dF_all,baseline_all)*100;
%get dF/F for average
baseline = mean(av(1.5*nOFF/10:(nON+nOFF)/10:end,:));
dF = bsxfun(@minus,av,baseline);
ratio = bsxfun(@rdivide,dF,baseline)*100;

%plot overlay of all trials with average
for iCell = 1:allcells;
    figure;
        for trial = 1:rep;
        plot(ratio_all(((1+((trial-1)*(epoch/10))):epoch*trial/10), iCell), 'c');
        hold on;
        end;
    plot(ratio(:, iCell),'k');
    hold on;
end

%plot overlay of trials by condition
for iCell = 1:allcells;
    figure;
    title(['Cell ' num2str(iCell)])
    for cond = 1:nCond;
        ratio_all_cond = zeros(epoch/10, rep/nCond);
        for trial = 1+((rep/nCond)*(cond-1)):cond*(rep/nCond);
            subplot(2,2,cond);
            plot(ratio_all(((1+((trial-1)*(epoch/10))):epoch*trial/10), iCell), 'c');
            hold on;
            ratio_all_cond(1:(epoch/10),trial-((cond-1)*(rep/nCond))) = ratio_all(((1+((trial-1)*(epoch/10))):epoch*trial/10), iCell);
        end;
        av_ratio_all_cond = mean(ratio_all_cond,2);
        plot(av_ratio_all_cond,'k')
    end;
    subplot(2,2,1);
    title('motor at 4V');
    subplot(2,2,2);
    title('motor at 0V');
    subplot(2,2,3);
    title('motor at 8V'); 
    subplot(2,2,4);
    title('motor at 2V')
end

