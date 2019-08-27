%find the number of trials and the duration in minutes of each trial used
%in the cc paper

clear;
%file_info_CRP_all
%WF_CRP_list_of_days;
%days_1 = {'170321_img86', '170408_img88', '170605_img95', '170612_img96', '170628_img98', '170705_img99'}
%days_post = {'170325_img86', '170415_img88', '170614_img95', '170622_img96', '170704_img98', '170710_img99'};

 days_1 = {'170321_img86',   '170408_img88', '170513_img89', '170416_img90',  '170417_img91', '170420_img92', ... 
     '170510_img93', '170524_img94', '170605_img95', '170612_img96', '170628_img98', '170705_img99', '170921_img044', ...
     '171113_img050', '171227_img067', '180104_img070', '180322_img077', '180507_img081', '180425_img084', '180509_img085', ...
     '181213_img087', '181214_img088', '181214_img089', '181213_img091', '190213_img092', '190215_img093', '190310_img094'};  %'170408_img87' first day of imaging but animal quit licking part way through
 
 days_post = {'170325_img86', '170415_img88', '170519_img89', '170426_img90',  '170425_img91', '170428_img92', ... 
     '170518_img93', '170529_img94', '170614_img95', '170622_img96', '170704_img98', '170710_img99' '170926_img044', ...
    '171122_img050', '180104_img067', '180108_img070',  '180403_img077', '180518_img081', '180502_img084', '180517_img085', ...
    '181217_img087', '181219_img088', '181218_img089', '181217_img091', '190217_img092', '190219_img093', '190315_img094'}; %'170326_img86',  are also 10%omission post learning

days_UR = {'170506_img90', '170426_img91', '170429_img92', '170519_img93', '170522_img89', '170530_img94', '181218_img087', ...
    '181219_img089', '181220_img091', '170327_img86', '170417_img88', '170623_img96', '170705_img98', '170711_img99', '190218_img092', '190220_img093', '190316_img094'};
 
data_dir = 'Y:\home\jake\Data\WidefieldImaging\GCaMP\behavior\';
trial_num_1 = [];
trial_num_N = [];
trial_num_UR = [];
sess_dur_1 = [];
sess_dur_N = [];
sess_dur_UR = [];

%% find the duration and # of trials for imaged sessions

for ii = 1:length(days_UR)
        if ~isempty(days_UR{ii})
                b_data = get_bx_data(data_dir, days_UR{ii});
                trial_num_UR = [trial_num_UR, b_data.trialSinceReset];
                sess_dur_UR = [sess_dur_UR, round((b_data.tThisTrialStartTimeMs{end} - b_data.tThisTrialStartTimeMs{1})/1000)/60];
        end
end

for ii = 1:length(days_1)
        if ~isempty(days_1{ii})
                b_data = get_bx_data(data_dir, days_1{ii});
                trial_num_1 = [trial_num_1, b_data.trialSinceReset];
                sess_dur_1 = [sess_dur_1, round((b_data.tThisTrialStartTimeMs{end} - b_data.tThisTrialStartTimeMs{1})/1000)/60];
        end
end

for ii = 1:length(days_post)
        if ~isempty(days_post{ii})
                b_data = get_bx_data(data_dir, days_post{ii});
                trial_num_N = [trial_num_N, b_data.trialSinceReset];
                sess_dur_N = [sess_dur_N, round((b_data.tThisTrialStartTimeMs{end} - b_data.tThisTrialStartTimeMs{1})/1000)/60];
        end
end

disp(['minimum and maximum number of trials: min=', num2str(min([trial_num_1, trial_num_N, trial_num_UR])),' max=', num2str(max([trial_num_1, trial_num_N, trial_num_UR]))]);
[trial_num_mean, trial_num_sem] = get_mean_and_sem([trial_num_1, trial_num_N, trial_num_UR]);
disp(['mean trial number = ', num2str(trial_num_mean), ' +/- ', num2str(trial_num_sem), ' trials per session'])

disp(['shortest and longest session durations: min=', num2str(min([sess_dur_1, sess_dur_N, sess_dur_UR])),' max=', num2str(max([sess_dur_1, sess_dur_N, sess_dur_UR]))]);
[sess_dur_mean, sess_dur_sem] = get_mean_and_sem([sess_dur_1, sess_dur_N, sess_dur_UR]);
disp(['mean session duration = ', num2str(sess_dur_mean), ' +/- ', num2str(sess_dur_sem), ' minutes session'])

%% calculate same numbers for training sessions

cd('Y:\home\jake\Data\WidefieldImaging\GCaMP\behavior\');
trials_per_sess = [];
sess_dur_all=[];
 for mouse_num = 1:length(days_1)
     trials_this_animal = [];
     sess_dur_this_animal=[];
     %determine 3 digit mouse id
     if days_1{mouse_num}(end-2) == 'g'
         mouse_id = ['i9', days_1{mouse_num}(end-1:end)];
     elseif days_1{mouse_num}(end-2) =='0'
          mouse_id = ['i', days_1{mouse_num}(end-2:end)];
     end
     
     %find all bx files fore that mosue
     all_bx_files = dir(['*', mouse_id, '*']);
     for sess_num = 1:length(all_bx_files)
         this_date = str2num(all_bx_files(sess_num).name(11:16));
         if this_date > str2num(days_1{mouse_num}(1:6)) & this_date < str2num(days_post{mouse_num}(1:6))
             b_data = load([data_dir, all_bx_files(sess_num).name]);  b_data = b_data.input;
                trials_this_animal = [trials_this_animal, b_data.trialSinceReset];
                sess_dur_this_animal = [sess_dur_this_animal, round((b_data.tThisTrialStartTimeMs{end} - b_data.tThisTrialStartTimeMs{1})/1000)/60];
         end
     end
     trials_per_sess = [trials_per_sess, trials_this_animal];
     sess_dur_all = [sess_dur_all, sess_dur_this_animal];
 end

[trials_per_sess_mean, trials_per_sess_sem] = get_mean_and_sem(trials_per_sess);
disp('Training sessions only:')
disp(['mean trial number = ', num2str(trials_per_sess_mean), ' +/- ', num2str(trials_per_sess_sem), ' trials per session'])

[sess_dur_all_mean, sess_dur_all_sem] = get_mean_and_sem(sess_dur_all);
disp(['mean session duration = ', num2str(sess_dur_all_mean), ' +/- ', num2str(sess_dur_all_sem), ' trials per session'])


