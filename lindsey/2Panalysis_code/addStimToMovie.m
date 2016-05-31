data_avg = mean(data(:,:,1:100),3);
[data_out, data_reg] = stackRegister(data,data_avg);
save('C:\Users\lindsey\Desktop\150623_AW14_002_reg.mat', 'data_reg', 'input','-v7.3');
%% flashing stim
base_fr_list = [];
targ_fr_list = [];
nOn = input.nFramesOn;
nOff = input.nFramesOff;
tCyc = cell2mat(input.tCyclesOn);
cFS = celleqel2mat_padded(input.cFirstStim);
cTarg = celleqel2mat_padded(input.cTargetOn);
tGratingDir = cell2mat(input.tGratingDirectionDeg);

F_frames = [];
for it = 1:length(tGratingDir)
	F_frames = [F_frames cFS(it)-15:cFS(it)];
end
dataF = mean(data_reg(:,:,F_frames),3);

data_df = bsxfun(@minus,double(data_reg),dataF);
data_dfof = bsxfun(@rdivide,data_df,dataF);

for it = 1:length(tCyc)
    ind = cFS(it):nOn+cFS(it)-1;
    base_fr_list = [base_fr_list ind];
    for icyc = 2:tCyc(it)
        ind = cFS(it)+((icyc-1)*(nOn+nOff)):cFS(it)+((icyc-1)*(nOn+nOff))+nOn-1;
        base_fr_list = [base_fr_list ind];
    end
    ind = cTarg(it):nOn+cTarg(it)-1;
    targ_fr_list = [targ_fr_list ind];
end

info.S = sparseint;

for i = 1:size(data_dfof,3)
    dataSquish(:,:,i) = data_dfof(:,:,i)*info.S;
end

data_dfof_squish = data_dfof.*info.S; 
data_dfof_sub = data_dfof(6:255,151:650,:);
mean_lum = makeGabor(50,5,0,0,0);
data_dfof_sub(1:50, 1:50, :) = repmat(mean_lum, [1 1 size(data_dfof_sub,3)]);
gabor = makeGabor(50,5,0,0,1);
data_dfof_sub(1:50, 1:50, base_fr_list) = repmat(gabor, [1 1 length(base_fr_list)]);

targ_fr_list(isnan(targ_fr_list)) = [];
gabor = makeGabor(50,5,tGratingDir(ceil(i/nOn)),0,1);
data_dfof_sub(1:50, 1:50, targ_fr_list) = repmat(gabor, [1 1 length(targ_fr_list)]);

save('C:\Users\lindsey\Desktop\150623_AW14_002_movie.mat', 'data_dfof_sub','-v7.3');
%% direction  
% single trial
tGratingDir = cell2mat(input.tGratingDirectionDeg);
nOn = input.nScansOn;
nOff = input.nScansOff;

F_frames = [];
for it = 1:length(tGratingDir)
    if it*(nOn+nOff) < length(data_reg)
        F_frames = [F_frames (nOff/2)+((nOff+nOn).*(it-1))+1:((nOn+nOff).*(it-1))+nOff];
    end
end
dataF = mean(data_reg(:,:,F_frames),3);

data_df = bsxfun(@minus,double(data_reg),dataF);
data_dfof = bsxfun(@rdivide,data_df,dataF);

fr_list = zeros(length(tGratingDir), nOn);

for it = 1:length(tGratingDir)
    fr_list(it,:) = nOff+((nOff+nOn).*(it-1))+1:(nOn+nOff).*(it);
end

data_dfof_sub = data_dfof(6:255,151:650,:);
mean_lum = makeGabor(50,5,0,0,0);
data_dfof_sub(1:50, 1:50, :) = repmat(mean_lum, [1 1 size(data_dfof_sub,3)]);
for it = 1:length(tGratingDir)
    if it*(nOn+nOff) < length(data_dfof_sub)
        for i = 1:nOn
            phase = rem(i,15)./15;
            gabor = makeGabor(50,5,tGratingDir(it),phase,1);
            ind = fr_list(it,i);
            data_dfof_sub(1:50, 1:50, ind) = gabor;
        end
    end
end

% average
dirs = unique(tGratingDir);
sz = size(data_dfof);
resp_avg = zeros(sz(1), sz(2), nOn+nOff, length(dirs));
for idir = 1:length(dirs)
    ind = find(tGratingDir == dirs(idir));
    resp = zeros(sz(1), sz(2), nOn+nOff, length(ind));
    for i = 1:length(ind)
        it = ind(i);
        if it*(nOn+nOff) < length(data_dfof)
            resp(:,:,:,i) = data_dfof(:,:,(nOff/2)+((nOff+nOn).*(it-1))+1:(nOn+nOff).*(it));
        end
    end
    resp_avg(:,:,:,idir) = squeeze(mean(resp,4));
end

resp_avg = reshape(resp_avg, [sz(1) sz(2) (nOff+nOn).*length(dirs)]);
fr_list = zeros(length(dirs),nOn);
for idir = 1:length(dirs)
    fr_list(idir,:) = (nOff/2)+((nOff+nOn).*(idir-1))+1:((nOn+nOff).*(idir))-nOff/2;
end
resp_avg_sub = resp_avg(6:255,151:650,:);
for idir = 1:length(dirs)
    for i = 1:nOn
        phase = rem(i,15)./15;
        gabor = makeGabor(50,5,dirs(idir),phase);
        ind = fr_list(idir,i);
        resp_avg_sub(1:50, 1:50, ind) = gabor;
    end
end

        