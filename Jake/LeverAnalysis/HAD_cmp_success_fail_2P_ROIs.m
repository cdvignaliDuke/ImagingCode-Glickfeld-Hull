%expt info
date = '150703';
run = '_000_000';
mouse = 'img24';

%output directory
out_base = 'Z:\home\lindsey\Analysis\2P\Jake';
run_name = [date '_' mouse '_run' run(length(run)-2:end)];
out_path = fullfile(out_base,run_name);
dest =  fullfile(out_path,run_name);

image_dest  = [dest '_ROI.tif'];
img = readtiff(image_dest);

img_down = stackGroupProject(img,10);

%prep for pca
global stack
stack = single(img_down);
defaultopts = {'nComp',300,'BorderWidth',4};
options = cell2struct(defaultopts(2:2:end),defaultopts(1:2:end),2);
[ny,nx,nt]=size(stack);
roi = imborder([ny,nx],options.BorderWidth,0); 
fprintf('Masking edges... ');
stack= bsxfun(@times,stack,single(roi));
fprintf('Done\n');
% compute thin svd using randomized algorithm
pcs = stackFastPCA(1,options.nComp);
% save extracted components 
fprintf('Saving principal components\n');
save([dest '_pca_usv.mat'],'-struct','pcs');

%visualize pca components
nt = size(pcs.v,1);
figure;
sm = stackFilter(pcs.U,1.5);
ax=[];
for pc = 1:25;                   % in order to visualize additional PCs simply alter the range (e.g. 26:50) Then subtract the appropriate amount from pc in the next line
    ax(pc)=subplot(5,5,pc);
    imagesc(sm(:,:,[pc]));
end;
colormap gray;

%compute independent components
PCuse = [2:100];
mu = 0;
nIC = 32;
termtol = 1e-6;
maxrounds = 400;
mixedsig = pcs.v';
mixedfilters = pcs.U;
CovEvals = diag(pcs.s).^2;
[ica_sig, ica_filters, ica_A, numiter] = CellsortICA(mixedsig, ...
    mixedfilters, CovEvals, PCuse, mu, nIC,[],termtol,maxrounds);

dt = 1/frGetFrameRate;
tt = [0:nt-1]/frGetFrameRate;


%% TC amd ROI code
cs = permute(ica_filters,[2,3,1]);
sm = stackFilter(cs,1.5);
figure;
ind = 1;
sel = [1:32];    
for ic = sel
    subplot(8,4,ind);                 %change here too
    imstretch(sm(:,:,ic),[.5 .99],1.5);
    ind = ind+1;
    text(.8,.1,num2str(ic),'fontsize',12,'color','w','fontweight','bold','unit','norm');
end;

%segment from ICs
mask_cell = zeros(size(sm));
for ic = sel
    bwimgcell = imCellEditInteractive(sm(:,:,ic),[]);
    mask_cell(:,:,ic) = bwlabel(bwimgcell);
    close all
end

mask_cell_temp = reshape(mask_cell,[sz(1)*sz(2) nIC]);
for ic = sel
    ind = find(mask_cell_temp(:,ic));
    mask_cell_temp(ind,ic)=1;
end
mask_cell_temp = reshape(mask_cell_temp,[sz(1)*sz(2) nIC]);

data_tc = zeros(size(img_down,3), nIC);
for ic = sel;
    if sum(mask_cell_temp(:,ic),1)>0
        data_tc(:,ic) = stackGetTimeCourses(img_down, reshape(mask_cell_temp(:,ic), [sz(1) sz(2)]));
    end
end
data_corr = corrcoef(data_tc);

sz = size(img);
mask_all = zeros(1,sz(1)*sz(2));
count = 0;
for ic = 1:nIC
    ind_new = find(mask_cell_temp(:,ic))';
    if length(ind_new)>1
        ind_old = find(mask_all);
        overlap = ismember(ind_old,ind_new);
        ind_both = find(overlap);
        if length(ind_both)>1
            ic_match = unique(mask_all(ind_old(ind_both)));
            for im = 1:length(ic_match)
                if data_corr(ic, ic_match(im))> 0.8
                    count = count+1;
                    mask_all(ind_new) = ic_match(im);
                else
                    mask_all(ind_new) = ic;
                    mask_all(ind_old(ind_both)) = 0;
                end
            end
        else
             mask_all(ind_new) = ic;
        end
    end
end
figure; imagesc(reshape(mask_all,[sz(1) sz(2)]))

start = 1;
mask_final = zeros(size(mask_all));
for ic = 1:max(mask_all,[],2)
    ind = find(mask_all==ic);
    if length(ind)>0
        mask_final(ind)=start;
        start= start+1;
    end
end
    
data_tc = stackGetTimeCourses(img, reshape(mask_final,[sz(1) sz(2)]));
save([dest '_ROI_TCs.mat'],'data_tc', 'mask_final');

%load frame and lever info
data_dest = [dest '_parse_behavior.mat'];
load(data_dest)
frame_info_dest = [dest '_frame_times.mat'];
load(frame_info_dest);
b_data.input = input; clear input;

%Obtain a df/f TC from baseline times
data_tc = data_tc'; 
startT = round(b_data.input.counterTimesUs{1}(1)./1000);
tc_dfoverf = zeros(size(data_tc));    %this could be problematic due to the frame skipping issue
first_baseline = find(~isnan(lever.baseline_timesMs(1,:)),1, 'first');    %find the first trial / baseline_timesMs window that is not NaN
F_range = [];
for iT=2:length(lever.baseline_timesMs)-1;    %this could be problematic due to unremoved NaNs
    if ~isnan(lever.baseline_timesMs(1,iT));
        F_range = frame_info.counter(lever.baseline_timesMs(1,iT)):frame_info.counter(lever.baseline_timesMs(2,iT));
    elseif isempty(F_range)
        F_range = frame_info.counter(lever.baseline_timesMs(1,first_baseline)):frame_info.counter(lever.baseline_timesMs(2,first_baseline));
    end
    F_avg= mean(data_tc(:,F_range),2);
    t_range = frame_info.counter(cell2mat(b_data.input.tThisTrialStartTimeMs(iT))-startT):frame_info.counter(cell2mat(b_data.input.tThisTrialStartTimeMs(iT+1))-startT);
    t_df = bsxfun(@minus, double(data_tc(:,t_range)), F_avg);
    t_dfoverf = bsxfun(@rdivide, t_df, F_avg);
    tc_dfoverf(:,t_range) = t_dfoverf;
end 

% ---- do simple movie analysis
func = @mean;
pre_frames = 5;
post_frames = 10;

ts = (-pre_frames:post_frames)*1000/round(double(Sampeling_rate));
tot_frame = pre_frames + post_frames+1;

%successes
use_ev_success = trial_outcome.success_time;
if strcmp(b_data.input.trialOutcomeCell{1}, 'success')
    use_ev_success(1) = [];
elseif strcmp(b_data.input.trialOutcomeCell{end}, 'success')
    use_ev_success(end) = [];
end
%----- uncomment to use only events w/o lever press after release
%     use_ev_success = remove_events_by_lever_state(use_ev_success,  ...
%         lever.state, 10,ceil(post_frames*1000/Sampeling_rate), 0);
%------ uncomment to use event only w/o lever press before release time
%     use_ev_success = remove_events_by_lever_state(use_ev_success,  ...
%         lever.state, -ceil(pre_frames*1000/Sampeling_rate),0, 1);
%
success_movie = trigger_movie_by_event(tc_dfoverf, frame_info, ...
    use_ev_success, pre_frames, post_frames);
avg_success = squeeze(func(success_movie,1));
sem_success = squeeze(std(success_movie,1)./sqrt(size(success_movie,1)));

avg_success_all = mean(avg_success,1);
sem_success_all = std(avg_success,1)./sqrt(size(avg_success,1));

%failures
use_ev_fail = trial_outcome.early_time;
if strcmp(b_data.input.trialOutcomeCell{1}, 'failure')
    use_ev_fail(1) = [];
elseif strcmp(b_data.input.trialOutcomeCell{end}, 'failure')
    use_ev_fail(end) = [];
end
%----- uncomment to use only events w/o lever press after release
%     use_ev_fail = remove_events_by_lever_state(use_ev_fail,  ...
%         lever.state, 10,ceil(post_frames*1000/Sampeling_rate), 0);
%------ uncomment to use event only w/o lever press before release time
%     use_ev_fail = remove_events_by_lever_state(use_ev_fail,  ...
%         lever.state, -ceil(pre_frames*1000/Sampeling_rate),0, 1);
%

% -----trigger movie by early release
fail_movie = trigger_movie_by_event(tc_dfoverf, frame_info, ...
    use_ev_fail, pre_frames, post_frames);
avg_fail = squeeze(func(fail_movie,1));
sem_fail = squeeze(std(fail_movie,1)./sqrt(size(fail_movie,1)));

%average of all ROIs
avg_fail_all = mean(avg_fail,1);
sem_fail_all = std(avg_fail,1)./sqrt(size(avg_fail,1));
tt =-pre_frames:post_frames;
figure; errorbar(tt,avg_success_all, sem_success_all,'k')
hold on;
errorbar(tt,avg_fail_all, sem_fail_all,'r')

%average by ROI
nCells = size(data_tc,1);
z = ceil(sqrt(nCells));
figure;
for ic = 1:nCells
    subplot(z,z,ic)
    errorbar(tt,avg_success(ic,:), sem_success(ic,:),'k')
    hold on;
    errorbar(tt,avg_fail(ic,:), sem_fail(ic,:),'r')
    ylim([-0.03 0.05])
end

save([dest '_resp_by_outcome.mat'],'fail_movie','success_movie');
