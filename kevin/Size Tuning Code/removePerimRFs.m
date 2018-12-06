%% code to re-write ret indices after removing RFs along perimeter
% strategy:
% load each experiment successively
% load:
% input or Els/Azs
% lbub_fits (ret)
% goodfit_ind, goodfit_ind_size
% if lbub_fits(goodfit_ind,4,4) == max/minAzs and (i,5,4) for els
% then mark goodfit_ind2 at zero
% also check if index within goodfit_ind matches value in goodfit_ind_size
% in which case, remove that value from goodfit_ind_size as well
% also need to remove rows for size tuning data arrays

%% load csv to define exp

clear all; clc;
expfile = '\\CRASH.dhe.duke.edu\data\home\kevin\Code\Ai9x_fixret_list.txt';
fID = fopen(expfile);
head = textscan(fID,'%s%s%s%s%s',1,'delimiter',',');
head = vertcat(head{:});
temp = textscan(fID,'%s%s%s%s%s','delimiter',',','HeaderLines',1);
temp = horzcat(temp{:});
expdata = cell2table(temp,'VariableNames',head);
nExp = size(expdata,1);
%isvalid = ones(1,nExp);
%expdata = addvars(expdata,isvalid);

fprintf(['Size-tuning visual-area comparison analysis - by KM, Glickfeld Lab\nLoading ' num2str(nExp) ' experiments\n'])

%% load each experiment and concatenate data

fprintf('\nBegin loading each experiment data...\n')

for i=1:nExp
    fprintf(['Exp: ' num2str(i) '/' num2str(nExp) '...'])
    date = expdata.date{i};
    mouse = expdata.mouse{i};
    fprintf(' mouse:%s, date:%s...',mouse,date)
    
    % lbub_fits
    fprintf('Loading: lbub_fits')
    run_str = 'runs-002'; %expdata.run_str{i};
    filename = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_lbub_fits.mat']);
    if ~exist(filename, 'file')
        fprintf([[date '_' mouse '_' run_str '_lbub_fits.mat'] ' not found! Please remove from list\n'])
    end
    load(filename)
    
    fprintf(', Azs/Els')
    filename = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOVs.mat']);

    if ~exist(filename, 'file')
        fprintf([[date '_' mouse '_' run_str '_FOVs.mat'] ' not found! Please remove from list\n'])
    end
    load(filename, 'Azs', 'Els')
    
    fprintf(['\n  #Good cells: initial = ' num2str(length(goodfit_ind))])
    
    goodfit_ind2 = zeros(size(goodfit_ind));
    for i=1:length(goodfit_ind)
        if sum(round(lbub_fits(goodfit_ind(i),4,4))==[min(Azs) max(Azs)])
            continue
        elseif sum(round(lbub_fits(goodfit_ind(i),5,4))==[min(Els) max(Els)])
            continue
        end
        goodfit_ind2(i) = goodfit_ind(i);
    end
    goodfit_ind2(goodfit_ind2==0) = [];
    goodfit_ind = goodfit_ind2;
    fprintf([', final = ' num2str(length(goodfit_ind)) ' - Saving fixed fits\n'])
    
    fn_out = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_lbub_fits.mat']);
    save(fn_out, 'lbub_fits', 'lbub_diff', 'goodfit_ind', 'resp_ind')
    
end
fprintf(['\nFinished checking all ' num2str(nExp) ' experiments.\n'])
