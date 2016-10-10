clear

DATA_DIR =  'C:\Users\jake\TempData\';
FRAME_TIME_DIR = 'C:\Users\jake\TempData\';
BEHAVE_DIR = 'Z:\Data\WidefieldImaging\GCaMP\behavior\';
ANALYSIS_DIR ='Z:\Analysis\LeverAnalysis\';

days = {'150716_img28', '150717_img28', '151021_img29', '151022_img29', '151009_img30', '151011_img30', '151211_img32', '151212_img32', '160129_img35', '160131_img35', '160129_img36','160131_img36', '160314_img38', '160315_img38', '160319_img41', '160320_img41'}; %'150718_img27', '150719_img27',
%days = {'160606_img46'};

for kk=1:length(days)
    destySucc = strcat(ANALYSIS_DIR, 'LeverSummaryFolder\', days{kk}, '_success');
    destyFail = strcat(ANALYSIS_DIR, 'LeverSummaryFolder\', days{kk}, '_fail');
    destyFidget = strcat(ANALYSIS_DIR, 'LeverSummaryFolder\', days{kk}, '_fidget');
    destyTooFast = strcat(ANALYSIS_DIR, 'LeverSummaryFolder\', days{kk}, '_tooFast');
    load(destySucc);
    load(destyFail);
    load(destyFidget);
    load(destyTooFast);
    if exist(strcat(ANALYSIS_DIR, 'LeverSummaryFolder\', days{kk}, '_lapse.mat'));
        destyLapse = strcat(ANALYSIS_DIR, 'LeverSummaryFolder\', days{kk}, '_lapse');
        load(destyLapse);
        if length(size(lapse_roi))==3
            num_lapse = size(lapse_roi,1);
        else
            num_lapse = 1;
        end
    else 
        num_lapse=0;
    end
    total_trials = [size(success_roi,1) + size(fail_roi,1)+ size(fidget_roi,1)+ size(tooFast_roi,1)+ num_lapse];
    disp(['day/animal: ' num2str(days{kk})])
    %disp(['% correct trials = ' num2str(size(success_roi,1)/total_trials)])
    %disp(['% early trials = ' num2str(size(fail_roi,1)/total_trials)])
    disp(['% fidget trials = ' num2str(size(fidget_roi,1)/total_trials)])
   % disp(['% tooFast correct = ' num2str(size(tooFast_roi,1)/total_trials)])
    %disp(['% lapse trials = ' num2str(num_lapse/total_trials)])
end