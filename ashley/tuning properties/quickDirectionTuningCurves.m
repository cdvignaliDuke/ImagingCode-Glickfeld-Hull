SubNum = expt(iexp).SubNum;
mouse = expt(iexp).mouse;
expDate = expt(iexp).date;
title_str = [mouse '-' expDate];
dirFolder = expt(iexp).dirtuning;
dirTime = expt(iexp).dirtuning_time;
down = 10;

fntun = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging', expDate, dirFolder);
fn = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging', expDate);

%% load tuning data
fName = [dirFolder '_000_000'];
[input] = Load_SBXdataPlusMWorksData(SubNum,expDate,dirTime,mouse,dirFolder,fName);  

load(fullfile(fntun,'timecourses.mat'))

%% some expt variables
if iscell(input.nScansOn)
    on = unique(cell2mat(input.nScansOn))./down;
    off = unique(cell2mat(input.nScansOff))./down;
else    
    on = input.nScansOn./down;
    off = input.nScansOff./down;
end
tDirections = cell2mat(input.tGratingDirectionDeg);
% [dir_ind directions] = findgroups(tDirections);
ncells = size(data_tc_subnp,2);
nfr = size(data_tc_subnp,1);
%% number of trials
ntrials = length(tDirections);
if (off+on)*ntrials > nfr 
    ntrials = floor(nfr./(off+on));
    data_tc_subnp = data_tc_subnp(1:((off+on).*ntrials),:);
    tDirections = tDirections(1:ntrials);
elseif (off+on)*ntrials < nfr
    data_tc_subnp = data_tc_subnp(1:((off+on).*ntrials),:);
end
[dir_ind directions] = findgroups(tDirections);
nstim = length(directions);

%% trial dF/F
data_tr = reshape(data_tc_subnp,on+off,ntrials,ncells);

off_ind = floor(off/2):off;

dFF = getDFoverF(data_tr,off_ind);

%% stim response 

resp_ind = off+1:off+on;

resp_dir = zeros(ncells,nstim);
resp_dir_err = zeros(ncells,nstim);
for istim = 1:nstim
    
   ind = find(dir_ind == istim);
   d = mean(dFF(resp_ind,ind,:,1));
   resp_dir(:,istim) = squeeze(mean(d,2));
   resp_dir_err(:,istim) = squeeze(std(d,[],2)./sqrt(length(ind)));
    
end

%% plot rand sample of tuning curves

exCells = randsample(ncells,9);

figure; setFigParams4Print('landscape')
suptitle({title_str;'example cells'})
for iplot = 1:9
    
   subplot(3,3,iplot)
   h = errorbar(directions,resp_dir(exCells(iplot),:),resp_dir_err(exCells(iplot),:),'ko-');
   h.Parent.XTick = directions;
   xlabel('direction')
   ylabel('dF/F')
   xlim([directions(1)-10 directions(end)+10]);
   title(num2str(exCells(iplot)))
    
end
print(['Z:\Analysis\FSAV Summaries\awFSAVdatasets_V1\quick tuning\' title_str '_exCellsTun'],'-dpdf','-fillpage')
