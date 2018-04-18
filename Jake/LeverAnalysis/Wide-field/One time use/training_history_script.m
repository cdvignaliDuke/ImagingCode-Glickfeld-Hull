% need to determine the # of training days before the first imaging session

min_corr_day = [];
bx_dir = ['Z:\Data\WidefieldImaging\GCaMP\behavior\BehaviorHistory\'];

%Normal(2):
date_2P = {'150703', '150704','150703', '150704', '151016', '151018', '151216', '160202','160203', '160205', '160319', '160320','160324',   '160327',  '160607', ...
    '160609', '160723', '160901', '161102', '161103', '171115', '171122', '180102', '180103', '180105', '180107', '180110', '180109', '180111', '180109'};
mouseID = {'img24','img24', 'img25', 'img25', 'img30', 'img30', 'img32', 'img36', 'img36','img36', 'img38', 'img38', 'img41','img41', 'img46', 'img46', 'img53', ...
    'img55', 'img59', 'img59', 'img043', 'img044', 'img060', 'img060',  'img063', 'img063', 'img064', 'img064', 'img065', 'img065'};

days_2P = {};
for session_num = 1:length(date_2P)
    days_2P{session_num} = strcat(date_2P{session_num}, '_',mouseID{session_num});
end
% %Omission:
% date={'180104',  '180109', '180110'}
% mouseI={'img060', 'img063', 'img065'}
% percent corr = 88.93,  61.37,  75.00

%Normal(WF):
WF_plotting_lists_of_days;
days_WF = days;
mouseIDs_WF = {'img29', 'img30', 'img32', 'img35', 'img36', 'img38', 'img41', 'img46', 'img53', 'img55'}; %'150718_img27', '150719_img27', '150716_img28', '150717_img28', 

%find the training day of the first time each animal achieved 60% correct
all_IDs_unique = unique([mouseIDs_WF, mouseID]);
for this_mouse=1:length(all_IDs_unique)
    %get mouse #
    training_day_num = 0;
    current_mouse = all_IDs_unique{this_mouse};
    if strcmp(current_mouse(4), '0');
        subj_num = ['i', current_mouse(end-2:end)];
    else
        subj_num = ['i9', current_mouse(end-1:end)];
    end
    
    %set directory
    current_pathname = [bx_dir, all_IDs_unique{this_mouse}];
    cd(current_pathname);
    training_history = dir; 
    while strcmp(training_history(1).name, '.') | strcmp(training_history(1).name, '..')
        training_history(1) = [];
    end
    if strcmp(training_history(end).name, 'pdfs')
        training_history(end) = [];
    end
    
    %assert that files are in chronological order
    assert(str2num(training_history(1).name(11:16))<=str2num(training_history(2).name(11:16))); 
    assert(str2num(training_history(2).name(11:16))<=str2num(training_history(3).name(11:16)));
    assert(str2num(training_history(end-1).name(11:16))<=str2num(training_history(end).name(11:16)));
    
    for curr_day = 1:length(training_history)
        %load session data
        bx_data = load([current_pathname, '\', training_history(curr_day).name]);
        bx_data = bx_data.input;
        if ~isfield(bx_data, 'juiceTimesMsCell')
            continue
        end
        training_day_num = training_day_num+1;
        
        %find percent correct
        trial_num = length(bx_data.juiceTimesMsCell);
        num_corr = sum(~cellfun(@isempty, bx_data.juiceTimesMsCell));  % find all the correct trials
        percent_corr = num_corr/trial_num;
        
        %keep searching till %corr>60% and min hold duration met
        if percent_corr < 0.6 | bx_data.randReqHoldMaxMs < 900
            %if both conditions are never met
            if curr_day == length(training_history)
                min_corr_day = [min_corr_day, NaN];
            end
            continue
        end
            
        %log the training day number
        min_corr_day = [min_corr_day, training_day_num];
        break
    end
end

%determine the number of trials on imaging days 
num_trials_per_session = [];
multi_file_sessions{1} = [];
all_days = [days_WF, days_2P];
bx_dir = ['Z:\Data\WidefieldImaging\GCaMP\behavior\BehaviorHistory\all animals\'];
for session_num = 1:length(all_days)
    %check for multiple days
    if all_days{session_num}(end-2) == 'g'
        bfile = dir([bx_dir 'data-i' '*' all_days{session_num}(end-1:end) '-' all_days{session_num}(1:6) '*' ]);
    elseif all_days{session_num}(end-2) =='0'
         bfile = dir([bx_dir 'data-i' '*' all_days{session_num}(end-2:end) '-' all_days{session_num}(1:6) '*' ]);
    end
    
    
    if length(bfile) > 1
        multi_file_sessions = [multi_file_sessions, all_days{session_num}];
        num_trials_per_session = [num_trials_per_session, NaN];
        continue
    end
    
    bx_data = get_bx_data(bx_dir, all_days{session_num});
    num_trials = length(bx_data.trialOutcomeCell);
    num_trials_per_session = [num_trials_per_session, num_trials];
end


% %determine the training day of each imaging session
% for sessions_num=1:length(all_days)
%     current_session = all_days(sessions_num);
%     subj_name = all_days{sessions_num}(8:end);
%     subj_num  = ['i9', all_days{sessions_num}(end-1:end)];
%     session_date = all_days{sessions_num}(1:6);
%     current_pathname = [bx_dir, num2str(subj_name)];
%     cd(current_pathname);
%     training_history = dir;
% end













