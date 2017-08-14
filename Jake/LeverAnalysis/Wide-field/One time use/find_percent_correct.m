%scipt for finding %corr and peak %corr and average %corr of any given
%imaging day
if 0>1
clear
summary_dir = ['Z:\Analysis\WF Lever Analysis\LeverSummaryFolder\'];
%look through bx data where the date is < target day
% divide successful trials by the total # of trials 

%list of all imaging days post img28
all_days = {'151009_img30', '151011_img30', '151013_img30', '151015_img29','151019_img30','151021_img29','151022_img29','151026_img30','151210_img32','151211_img32','151212_img32','151222_img32','160129_img36','160129_img35','160131_img35','160131_img36','160205_img35','160207_img35','160207_img36','160208_img35','160208_img36','160209_img36','160228_img36','160314_img38','160315_img38','160319_img41','160320_img41','160512_img47','160515_img44','160515_img47','160515_img48','160516_img47','160606_img46','160704_img54','160722_img53', '160725_img53'}; 
for ii = 1:length(all_days)
    clear Ncorr NtooFast Nfail Nfidget Nlapse correct tooFast fail fidget lapse
   %find # of correct trials 
   if exist([summary_dir, all_days{ii}, '_success.mat'], 'file') == 2
       correct = load([summary_dir, all_days{ii}, '_success.mat']);
       Ncorr = size(correct.success_roi, 1);
   end
   
   %find # of tooFast trials
   if exist([summary_dir, all_days{ii}, '_tooFast.mat'], 'file') == 2
       tooFast = load([summary_dir, all_days{ii}, '_tooFast.mat']);
       if ndims(tooFast.tooFast_roi) ==3
           NtooFast = size(tooFast.tooFast_roi, 1);
       elseif ndims(tooFast.tooFast_roi) ==2
           NtooFast = 1; 
       elseif isempty(tooFast.tooFast_roi)
           NtooFast = 0; 
       end
   else
       NtooFast = 0;
   end
   
   %find # of early trials
   if exist([summary_dir, all_days{ii}, '_fail.mat'], 'file') == 2
       fail = load([summary_dir, all_days{ii}, '_fail.mat']);
       Nfail = size(fail.fail_roi, 1);
   end
   
   %find # of fidget trials
   if exist([summary_dir, all_days{ii}, '_fidget.mat'], 'file') == 2
       fidget = load([summary_dir, all_days{ii}, '_fidget.mat']);
       Nfidget = size(fidget.fidget_roi, 1);
   else 
       Nfidget = 0;
   end
   
   %find # of lapse trials
   if exist([summary_dir, all_days{ii}, '_lapse.mat'], 'file') == 2
       lapse = load([summary_dir, all_days{ii}, '_lapse.mat']);
       if ndims(lapse.lapse_roi) ==3
           Nlapse = size(lapse.lapse_roi, 1);
       elseif ndims(lapse.lapse_roi) ==2
           Nlapse = 1;
       elseif isempty(lapse.lapse_roi)
           Nlapse = 0;
       end
   else
       Nlapse = 0; 
   end
   if exist('Ncorr')
   Ntotal = Ncorr + NtooFast + Nfail  +Nfidget + Nlapse;
   percent_corr = Ncorr/ Ntotal;
   disp([all_days{ii},'--', num2str(round(percent_corr*100)) '--percent correct--'])
   else
       disp([all_days{ii}, ' No data'])
   end
end
end 

%% find the peak %correct and the avg %correct of the previous 3 days 
clear
%list of all imaging days post img28
all_days = {'151009_img30', '151011_img30', '151013_img30', '151015_img29','151019_img30','151021_img29','151022_img29','151026_img30','151210_img32','151211_img32','151212_img32','151222_img32','160129_img36','160129_img35','160131_img35','160131_img36','160205_img35','160207_img35','160207_img36','160208_img35','160208_img36','160209_img36','160228_img36','160314_img38','160315_img38','160319_img41','160320_img41','160512_img47','160515_img44','160515_img47','160515_img48','160516_img47','160606_img46','160704_img54','160722_img53', '160725_img53'}; 
all_days = {'160904_img55'};
bx_dir = ['Z:\Data\WidefieldImaging\GCaMP\behavior\BehaviorHistory\'];
for ii = 1:length(all_days)
    current_session = all_days(ii);
    subj_name = all_days{ii}(8:end);
    subj_num  = ['i9', all_days{ii}(end-1:end)];
    session_date = all_days{ii}(1:6);
    current_pathname = [bx_dir, num2str(subj_name)];
    cd(current_pathname);
    training_history = dir;
    for iii = size(training_history,1):-1:1  %remove data which are just dots. Go in reverse order to avoid errors in shifting indeces
        if size(training_history(iii).name,2) < 5
            training_history(iii) = [];
        end
    end
    %find peak %corr
    %this section assumes that the files will be in numerical order by date
    %with the earliest sessions first 
    assert(str2num(training_history(1).name(11:16))<=str2num(training_history(2).name(11:16))); 
    assert(str2num(training_history(2).name(11:16))<=str2num(training_history(3).name(11:16)));
    assert(str2num(training_history(end-1).name(11:16))<=str2num(training_history(end).name(11:16)));
    %find the target day's bx file
    recent_percent_corr = [];
    peak_percent_corr = 0; 
    for iii = 1:size(training_history,1) %search through all the training days for this animal previous to the current imaging day
        if str2num(training_history(iii).name(11:16))==str2num(session_date);
            target_file_inx = iii;
            break
        end
    end
    for iii = target_file_inx:-1:1
        load([bx_dir, subj_name, '\', training_history(iii).name]);
        bx_data = input;
        curr_percent_corr = length(cell2mat(bx_data.juiceTimesMsCell))/length(bx_data.juiceTimesMsCell);
        if curr_percent_corr > peak_percent_corr & length(bx_data.juiceTimesMsCell)>100; %have to have done at least 100 trials that day.
            peak_percent_corr = curr_percent_corr;   %log the peak percent correct for this animal previous to the current imaging session
        end
        if iii==target_file_inx-1 | iii == target_file_inx-2 | iii == target_file_inx-3 % look at the 3 training days previous to the current imaging day
            recent_percent_corr = [recent_percent_corr, curr_percent_corr];
        end
    end
    disp([all_days(ii), ' peak % correct=', num2str(peak_percent_corr), ' recent % correct=', num2str(mean(recent_percent_corr))])
end




