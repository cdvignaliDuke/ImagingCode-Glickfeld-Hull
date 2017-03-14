%% some variables from mworks
if iscell(input.nScansOn)
    on = unique(cell2mat(input.nScansOn))./down;
    off = unique(cell2mat(input.nScansOff))./down;
else
    on = input.nScansOn./down;
    off = input.nScansOff./down; 
end
off_starts = 1:off+on:size(data_reg,3);
on_starts = off+1:off+on:size(data_reg,3);

tDir = cell2mat(input.tGratingDirectionDeg);
[tDir_ind dir] = findgroups(tDir);
ntrials = length(tDir);
nstim = length(dir);

%% dF/F by trial

if (off+on)*ntrials > nfr 
    ntrials = floor(nfr./(off+on));
    [tDir_ind dir] = findgroups(tDir(1:ntrials));
    data_reg = data_reg(:,:,1:((off+on).*ntrials));
    data_tr = reshape(data_reg,ypix,xpix,off+on,ntrials);
elseif (off+on)*ntrials < nfr
    ntrials = length(tDir);
    [tDir_ind dir] = findgroups(tDir(1:ntrials));
    data_reg = data_reg(:,:,1:((off+on).*ntrials));
    data_tr = reshape(data_reg,ypix,xpix,off+on,ntrials);
else
    data_tr = reshape(data_reg,ypix,xpix,off+on,ntrials);
end

F = mean(data_tr(:,:,off-6:off,:),3);
dF = bsxfun(@minus,data_tr,F);
dFF = bsxfun(@rdivide,dF,F);

maxDFF = max(squeeze(mean(dFF(:,:,off+1:off+on,:),3)),[],3);

%% ******choose crop parameters*******
figure;colormap gray; imagesc(maxDFF)

% adjust maxDFF by percentage of median
adj_low = 6.*median(maxDFF(:));
maxDFF_adj = maxDFF;
maxDFF_adj(maxDFF_adj > adj_low) = adj_low;
figure;colormap gray; imagesc(maxDFF_adj)