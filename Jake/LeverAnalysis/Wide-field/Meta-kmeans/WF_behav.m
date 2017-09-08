%
% i = 14;
% figure;plot(input.lickometerTimesUs{1,i}/1000, ones(size(input.lickometerValues{1,i},2),1),'*')
% hold on;
% plot(input.holdStartsMs{1,i}, ones(size(input.holdStartsMs{1,i},2),1),'r*');
% plot(input.holdStartsMs{1,i} + input.holdTimesMs{1,i}, ones(size(input.holdStartsMs{1,i},2),1),'k*');
%
% min_diff1 = [];
% min_diff2 = [];
% for i = 5:220
%     if ~isempty(input.lickometerTimesUs{1,i})
%         press_lick = input.lickometerTimesUs{1,i}/1000 - input.holdStartsMs{1,i};
%         rel_lick = input.lickometerTimesUs{1,i}/1000 - (input.holdStartsMs{1,i} + input.holdTimesMs{1,i});
%
%         min_diff1 = [min_diff1; press_lick(press_lick == min(abs(press_lick)))];
%         min_diff2 = [min_diff2; rel_lick(rel_lick==min(abs(rel_lick)))];
%
%     end
% end
% %
% sum(min_diff1 < 101)/size(min_diff1,1)
%
% sum(min_diff2 < 101)/size(min_diff2,1)


% subID = '160904_img55'; %'160722_img53';  % '160606_img46';
% load(['\\crash.dhe.duke.edu\data\home\jake\Analysis\WF Lever Analysis\kmeans_output_dir\', subID, '\', subID,'_downsampled_cropped_motionCorr.mat']);
% 
% % load('Z:\home\andrew\Behavior\Data\data-i946-160606-1531.mat');
% % load('Z:\home\andrew\Behavior\Data\data-i953-160722-1700.mat');
% load('Z:\home\andrew\Behavior\Data\data-i955-160904-1736.mat');

sFrame_r = []; fFrame_r = []; sFrame_pr = []; fFrame_pr = [];
lFrame = []; rFrame = []; iFrame = []; fail_trial = []; succ_trial = [];

trialS = max(min(find(cell2mat(input.cLeverDown)~=0)), min(find(cell2mat(input.cLeverUp)~=0))); trialE = length(input.cLeverDown);
randTrial = randi([trialS, trialE], 1, trialE - trialS + 1);

[x,y,nMovie] = size(downsampled_movie);
df_f_norm = normalizeTC(downsampled_movie);
df_f_frame = reshape(df_f_norm, x, y, nMovie);

for i = trialS:trialE
    if strcmp(input.trialOutcomeCell{i}, 'success')
        succ_trial = [succ_trial; i];
    end
    if strcmp(input.trialOutcomeCell{i}, 'failure')
        
        fail_trial = [fail_trial; i];
    end
end

% rand_i = fail_trial(randperm(length(fail_trial)))';
% rand_s = datasample(succ_trial,length(rand_i));
% 
% for i = randTrial
%     if input.cLeverUp{i} >= nMovie
%         continue;
%     else
%     % get lick frames after release > 0.5s and match ITI without licking
%      if ~isempty(input.lickometerTimesUs{1,i})
%          rel_lick = input.lickometerTimesUs{1,i}/1000 - input.holdStartsMs{1,i} - input.holdTimesMs{i};
%          if rel_lick(1) < 0
%              rel_lick = (rel_lick(rel_lick >=0));
%          end
%          if ~isempty(rel_lick) && length(rel_lick) >= 5
%              lFrame = cat(3, lFrame, downsampled_movie(:,:, input.cLeverUp{i} + floor(rel_lick(1)/100) : input.cLeverUp{i} + floor(rel_lick(end)/100)));
%              cWaitForPress = floor((input.tStartTrialWaitForPressTimeMs{1,i} - input.tThisTrialStartTimeMs{1,i})/100);
%              itiLen = length(input.cLeverUp{i} + floor(rel_lick(1)/100) : input.cLeverUp{i} + floor(rel_lick(end)/100));
%                           
%              cPresslick = floor((input.lickometerTimesUs{1,i}/1000 - input.tStartTrialWaitForPressTimeMs{1,i})/100);
%              
%              if  ~isempty(cPresslick) && cPresslick(end) < 0
%                  cPresslick = abs(cPresslick(end));
%                  
%                  if length(cWaitForPress - cPresslick + 2 + input.cTrialStart{i} : cWaitForPress + input.cTrialStart{i}) > itiLen
%                      iFrame = cat(3, iFrame, downsampled_movie(:,:, input.cTrialStart{i} + cWaitForPress - itiLen + 1 : cWaitForPress + input.cTrialStart{i}));
%                  else
%                     iFrame = cat(3, iFrame, downsampled_movie(:,:, cWaitForPress - cPresslick(end) + 2 + input.cTrialStart{i} : cWaitForPress + input.cTrialStart{i}));
%                  end
%              else
%                  iFrame = cat(3, iFrame, downsampled_movie(:,:, input.cTrialStart{i} + cWaitForPress - itiLen + 1 : cWaitForPress + input.cTrialStart{i}));
%              end
%          
%              
%          end
%          
%      end
%     end
% end
% 
% j = 1; endFrame = 1000/100;
% 
% for i = rand_i%randTrial  %2:220
%     
%     % get success and fail after release and match length
%     if input.cLeverUp{rand_s(j)} >= nMovie || input.cLeverUp{i} >= nMovie
%         continue;
%     else
%         sFrame_r = cat(3, sFrame_r, downsampled_movie(:,:,input.cLeverUp{rand_s(j)} : input.cLeverUp{rand_s(j)} + endFrame));
% 
%         fFrame_r = cat(3, fFrame_r, downsampled_movie(:,:,input.cLeverUp{i} : input.cLeverUp{i} + endFrame));
%    
%     % get success and fail between press and release and match length
%     fn = length(input.cLeverDown{i} : input.cLeverUp{i});
%     
%     sFrame_pr = cat(3, sFrame_pr, downsampled_movie(:,:,input.cLeverDown{rand_s(j)} : input.cLeverDown{rand_s(j)} + fn -1));
%     
%     fFrame_pr = cat(3, fFrame_pr, downsampled_movie(:,:,input.cLeverDown{i} : input.cLeverUp{i}));
%     
%     j = j+1;
%     end
% end

sFrame = []; fFrame = [];
startFrame = 100/100;
endFrame = 300/100;
% get all success trial one frame before release and 3 frames after release
for i = succ_trial'
    if input.cLeverUp{i} >= nMovie-1
        continue;
    else
        
        sFrame = cat(3, sFrame, df_f_frame(:,:,input.cLeverUp{i} - startFrame : input.cLeverUp{i} + endFrame));
    end
end

for i = fail_trial'
    if input.cLeverUp{i} >= nMovie
        continue;
    else
        fFrame = cat(3, fFrame, df_f_frame(:,:,input.cLeverUp{i} - startFrame : input.cLeverUp{i} + endFrame));
    end
end

% % outdir = ['\\crash.dhe.duke.edu\data\home\jake\Analysis\WF Lever Analysis\Meta-kmeans_output_dir\',subID, '\movie_ITIvsLick\'];
% outdir = [kmeans_output_dir, '\movie_ITIvsLick\'];
% if ~exist(outdir)
%     mkdir(outdir);
% end
% figName = 'cluster_lick_match_rand.fig';
% fileName = 'lick_clutstering_match_rand.mat';
% runMetaK(lFrame, outdir, figName, fileName)
% 
% outdir = ['\\crash.dhe.duke.edu\data\home\jake\Analysis\WF Lever Analysis\Meta-kmeans_output_dir\',subID, '\movie_ITIvsLick\'];
% if ~exist(outdir)
%     mkdir(outdir);
% end
% figName = 'cluster_iti_match_rand.fig';
% fileName = 'iti_clutstering_match_rand.mat';
% runMetaK(iFrame, outdir, figName, fileName)
% 
% outdir = ['\\crash.dhe.duke.edu\data\home\jake\Analysis\WF Lever Analysis\Meta-kmeans_output_dir\',subID, '\movie_SvsF_afterR\'];
% if ~exist(outdir)
%     mkdir(outdir);
% end
% figName = 'cluster_success_match_rand.fig';
% fileName = 'success_clutstering_match_rand.mat';
% runMetaK(sFrame_r, outdir, figName, fileName)
% 
% outdir = ['\\crash.dhe.duke.edu\data\home\jake\Analysis\WF Lever Analysis\Meta-kmeans_output_dir\',subID, '\movie_SvsF_afterR\'];
% if ~exist(outdir)
%     mkdir(outdir);
% end
% figName = 'cluster_fail_match_rand.fig';
% fileName = 'fail_clutstering_match_rand.mat';
% runMetaK(fFrame_r, outdir, figName, fileName)
%     
% outdir = ['\\crash.dhe.duke.edu\data\home\jake\Analysis\WF Lever Analysis\Meta-kmeans_output_dir\',subID, '\movie_SvsF_P2R\'];
% if ~exist(outdir)
%     mkdir(outdir);
% end
% figName = 'cluster_success_match_rand.fig';
% fileName = 'success_clutstering_match_rand.mat';
% runMetaK(sFrame_pr, outdir, figName, fileName)
% 
% outdir = ['\\crash.dhe.duke.edu\data\home\jake\Analysis\WF Lever Analysis\Meta-kmeans_output_dir\',subID, '\movie_SvsF_P2R\'];
% if ~exist(outdir)
%     mkdir(outdir);
% end
% figName = 'cluster_fail_match_rand.fig';
% fileName = 'fail_clutstering_match_rand.mat';
% runMetaK(fFrame_pr, outdir, figName, fileName)
        
    % get success and fail whole length movie
    %     if strcmp(input.trialOutcomeCell{i}, 'success')
    %         sFrame = cat(3, sFrame, downsampled_movie(:,:,input.cTrialStart{i} : input.cTrialEnd{i}));
    %     end
    %
    %     if strcmp(input.trialOutcomeCell{i}, 'failure')
    %         fFrame = cat(3, fFrame, downsampled_movie(:,:,input.cTrialStart{i} : input.cTrialEnd{i}));
    %     end
    
    % get lick frames before press
%     if ~isempty(input.lickometerTimesUs{1,i})
%         press_lick = input.lickometerTimesUs{1,i}/1000 - input.holdStartsMs{1,i};
%         if  press_lick(1) < 0
%             press_lick = abs(press_lick(press_lick<-400));
%             if ~isempty(press_lick)
%                 lFrame = cat(3, lFrame, downsampled_movie(:,:, input.cLeverDown{i} - floor(press_lick(1)/100) : input.cLeverDown{i} - floor(press_lick(end)/100)));
%                 
%             end
%         end
%     end
    

    % get release to trialend frame
    %rFrame = cat(3, rFrame, downsampled_movie(:,:, input.cLeverUp{i} - 1 : input.cTrialEnd{i}));
    % get release around lever release
%     rFrame = cat(3, rFrame, downsampled_movie(:,:, input.cLeverUp{i} - 4 : input.cLeverUp{i} + 4));
    
    % get ITI without licking
%     cWaitForPress = floor((input.tStartTrialWaitForPressTimeMs{1,i} - input.tThisTrialStartTimeMs{1,i})/100);
%     if ~isempty(input.lickometerTimesUs{1,i})
%         cPresslick = floor((input.lickometerTimesUs{1,i}/1000 - input.tStartTrialWaitForPressTimeMs{1,i})/100);
%         
%         if  ~isempty(press_lick) && press_lick(end) < 0
%             press_lick = abs(press_lick(end));
%             
%             iFrame = cat(3, iFrame, downsampled_movie(:,:, cWaitForPress - cPresslick + 2 + input.cTrialStart{i} : cWaitForPress + input.cTrialStart{i}));
%                 
%         else
%             iFrame = cat(3, iFrame, downsampled_movie(:,:, input.cTrialStart{i} : cWaitForPress + input.cTrialStart{i}));
%         end
%     else
%         iFrame = cat(3, iFrame, downsampled_movie(:,:, input.cTrialStart{i} : cWaitForPress + input.cTrialStart{i}));
%     end



