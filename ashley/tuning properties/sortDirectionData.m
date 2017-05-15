function [resp, data_trReshape,nstim] = sortDirectionData(data,down,input)

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
% ncells = size(data,2);
xpix = size(data,2);
ypix = size(data,1);
nfr = size(data,3);
%% number of trials
ntrials = length(tDirections);
if (off+on)*ntrials > nfr 
    ntrials = floor(nfr./(off+on));
    data = data(:,:,((off+on).*ntrials));
    tDirections = tDirections(1:ntrials);
elseif (off+on)*ntrials < nfr
    data = data(:,:,((off+on).*ntrials));
end
[dir_ind directions] = findgroups(tDirections);
nstim = length(directions);

%% analysis window
pre_win = floor(0.75*off):off;
post_win = off+1:off+on;

%% trial reshaped data matrix
data_trReshape = reshape(data,ypix,xpix,on+off,ntrials);

%% trial type sorted data structure
resp = struct;
for istim = 1:nstim
   resp(istim).on = [];
   resp(istim).off = [];
   trial_ind = dir_ind == istim;
   stim_data_temp = data_trReshape(:,:,:,trial_ind);
   resp(istim).off = squeeze(mean(stim_data_temp(:,:,pre_win,:),3));
   resp(istim).on = squeeze(mean(stim_data_temp(:,:,post_win,:),3));
end

end
% off_ind = floor(off/2):off;

% dFF = getDFoverF(data_tr,off_ind);

% %% stim response 
% 
% resp_ind = off+1:off+on;
% 
% resp_dir = zeros(ncells,nstim);
% resp_dir_err = zeros(ncells,nstim);
% for istim = 1:nstim
%     
%    ind = find(dir_ind == istim);
%    d = mean(dFF(resp_ind,ind,:,1));
%    resp_dir(:,istim) = squeeze(mean(d,2));
%    resp_dir_err(:,istim) = squeeze(std(d,[],2)./sqrt(length(ind)));
%     
% end
