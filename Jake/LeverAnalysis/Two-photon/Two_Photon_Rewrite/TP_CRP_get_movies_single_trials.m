% PLOT example movies from sbx read
clear all; clear global; clear cd;
img_dir = 'Z:\Data\2P_imaging\';
index_dir = 'Z:\Analysis\Cue_reward_pairing_analysis\2P\ITI licking\';
mouseID = 'img92';
dateID =  '170428';
fname = 'img92_000_000';

%% load index of events and sbx file

cd([img_dir dateID '_' mouseID, '\', mouseID]);
data = squeeze(sbxread(fname,0,10000));
if length(size(data)) ==4
    data = squeeze(data(1,:,:,:));
end
data_avg = mean(data,3);

%% motion correct 
%determine frame nums to use for selecting stable frame
sz = size(data); 
laser_on_ind = [1:sz(3)];
frame_nums_for_ref30 = laser_on_ind(randi([1,size(laser_on_ind,2)],1,30));
frame_nums_for_samp100 = laser_on_ind(round(linspace(1,size(laser_on_ind,2))));

% select 30 random frames from throughout the movie
ref30 = data(:,:,frame_nums_for_ref30);

%motion register each of the 30 random frames to 100 frames from the movie. Find the one with the lowerst dshift
samp100 = data(:,:,frame_nums_for_samp100);
dshift = [];
for r = 1:size(ref30,3)
    [reg_out,aa] = stackRegister(samp100, ref30(:,:,r));
    dshift = [dshift; mean(((reg_out(:,3).^2)+(reg_out(:,4).^2)))];
end

%pick the frame which had the lowest dshift and motion register full movie to that frame
min_f = find(dshift == min(dshift));
if length(min_f) > 1
    min_f = min_f(1,1);
end
img_ref = ref30(:,:,min_f);

%% plot individual trials and save as tiffs

load([index_dir, dateID, '_', mouseID, '\ITI_licking_outputs'], 'bout_start_inx');
if length(bout_start_inx)>10
    bout_num = 1:length(bout_start_inx);  %round(linspace(1, length(bout_start_inx), 8));
    for iii = bout_num;
        frame_num = bout_start_inx(iii);
        data = squeeze(sbxread(fname,[frame_num-30],91));
        [reg_out, img_reg] = stackRegister(data, img_ref);
        if iii == 1
            all_clips = NaN(size(img_reg,1), size(img_reg,2), size(img_reg,3), length(bout_start_inx));
        end
        all_clips(:,:,:,iii) = img_reg;
        %writetiff(img_reg, [index_dir, dateID, '_', mouseID, '\motion_corr_trial_', num2str(iii)]);
    end
end
mean_movie = nanmean(all_clips,4);
writetiff(mean_movie, [index_dir, dateID, '_', mouseID, '\motion_corr_mean_bout']);
    
    
    
    
    



