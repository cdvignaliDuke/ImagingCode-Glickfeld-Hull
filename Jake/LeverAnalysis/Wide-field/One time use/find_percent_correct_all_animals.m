%scipt for finding %corr and peak %corr and average %corr of any given
%imaging day
clear
bx_dir = ['Z:\Data\WidefieldImaging\GCaMP\behavior\BehaviorHistory\'];
summary_dir = ['Z:\Analysis\WF Lever Analysis\LeverSummaryFolder\'];

%list of all imaging days post img28
WF_plotting_lists_of_days;
days_WF = days;

date_2P = {'150703', '150704','150703', '150704', '151016', '151018', '151216', '160202','160203', '160205', '160319', '160320','160324',   '160327',  '160607', ...
    '160609', '160723', '160901', '161102', '161103', '171115', '171122', '180102', '180103', '180105', '180107', '180110', '180109', '180111', '180109'};
mouseID = {'img24','img24', 'img25', 'img25', 'img30', 'img30', 'img32', 'img36', 'img36','img36', 'img38', 'img38', 'img41','img41', 'img46', 'img46', 'img53', ...
    'img55', 'img59', 'img59', 'img043', 'img044', 'img060', 'img060',  'img063', 'img063', 'img064', 'img064', 'img065', 'img065'};
mouseIDs_WF = {'img29', 'img30', 'img32', 'img35', 'img36', 'img38', 'img41', 'img46', 'img53', 'img55'}; %'150718_img27', '150719_img27', '150716_img28', '150717_img28', 
all_IDs_unique = unique([mouseIDs_WF, mouseID]);

%make a list of all 2P sessions and all imaging sessions
for session_num = 1:length(date_2P)
    this_session = [date_2P{session_num}, '_', mouseID{session_num}];
    if session_num ==1
        sessions_2P{1} = this_session;
    else
        sessions_2P = [sessions_2P, this_session];
    end
end
all_sessions = [sessions_2P, days];

%find the first day of imaging for each animal
first_date_all=[];
for this_mouse=1:length(all_IDs_unique)
    this_mouse_dates=[];
    this_mouse_sessions_inx = find(~cellfun(@isempty, strfind(all_sessions, all_IDs_unique{this_mouse}))); %find the sessions for each mouse
    for session_num = 1:length(this_mouse_sessions_inx) %determine which session happened first
        this_mouse_dates = [this_mouse_dates, str2num(all_sessions{this_mouse_sessions_inx(session_num)}(1:6))];    
    end
    first_date = min(this_mouse_dates);
    first_date_all = [first_date_all, first_date];
    if this_mouse==1
        first_session_all{1} = [num2str(first_date), '_', all_IDs_unique{this_mouse}];
    else
        first_session_all = [first_session_all, [num2str(first_date), '_', all_IDs_unique{this_mouse}] ];
    end
end

%find the peak %correct for each mouse up to and including the first day of imaging 
peak_corr_all = [];
peak_corr_day_all = [];
for session_num = 1:length(first_session_all)
    %get mouse number
    current_session = first_session_all{session_num};
    if strcmp(current_session(end-2), '0');
        subj_num = ['i', current_session(end-2:end)];
    elseif strcmp(current_session(end-2), 'g');
        subj_num = ['i9', current_session(end-1:end)];
    end
    
    %set directory
    current_pathname = [bx_dir, current_session(8:end)];
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
    
    
    training_day_num = 0;
    curr_peak_corr = 0;
    curr_training_day_of_peak = 0;
    for curr_day = 1:length(training_history)
        %if the day in question happens after the first imaging day then stop looking 
        if first_date_all(session_num) <  str2num(training_history(curr_day).name(11:16));
            break
        end
        
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
        
        %if there are a minimum of 50 trials and  a minimum rand hold of 900 then log the %corr
        if trial_num>=50 & bx_data.randReqHoldMaxMs > 900
            if curr_peak_corr < percent_corr
                curr_peak_corr = percent_corr;
                curr_training_day_of_peak = training_day_num;
            end
        end
       
    end
    
    %log the peak %corr and the training day
    peak_corr_all = [peak_corr_all, curr_peak_corr];
    peak_corr_day_all = [peak_corr_day_all, curr_training_day_of_peak]; 
    % 3/28/18  mean %corr= 79.67 +/-0.23%    Day=26.05 +/-2.42 
end










