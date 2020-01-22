mouse = 'i1306';
date = '191016';
run = '003';
time = '1325';


%% load data
CD = fullfile('Z:\All_Staff\home\grace\2P_Imaging', [mouse], [date '_' mouse], run); %#ok<NBRAK>
cd(CD);


load([run '_000_000.mat'])

load(fullfile('Z:\All_Staff\Behavior\Data',['data-' mouse '-' date '-' time]));
 nframes = min([info.config.frames input.counterValues{end}(end)],[],2);
 data = sbxread([run '_000_000'],0,nframes);
 
 data = squeeze(data);
 
 %% choose stable image
 
 figure;
 t = floor(14000./6);
 for i = 1:6
     subplot(2,3,i)
     imagesc(mean(data(:,:,1+((i-1).*t):500+1+((i-1).*t)),3))
     title(num2str(1+((i-1).*t)))
 end
 
 i = 2;
 data_avg = mean(data(:,:,1+((i-1).*t):500+1+((i-1).*t)),3);
 
 %% register data 
 
 [out data_reg] = stackRegister(data,data_avg);
 
 fnout= fullfile('Z:\All_Staff\home\grace\Analysis\2P',[date '_' mouse], [date '_' mouse '_runs-003']);
 mkdir(fnout);
 save(fullfile(fnout,[date '_' mouse '_' run '_input.mat']),'input');
 save(fullfile(fnout,[date '_' mouse '_' run '_reg.mat']),'out','data_avg');
 
 %% 
 tGratingDirection = celleqel2mat_padded(input.tGratingDirectionDeg);
 nTrials = size(tGratingDirection,2);
 nOn = input.nScansOn;
 nOff = input.nScansOff;
 nFrames = nOn+nOff;
 sz = size(data_reg);
 
 data_trial = permute(reshape(data_reg,[sz(1) sz(2) nFrames nTrials]),[1 2 4 3]);
 data_f = mean(data_trial(:,:,:,ceil(nOff/2):nOff),4);
 data_dfof= bsxfun(@rdivide,bsxfun(@minus,double(data_trial),data_f),data_f);
 data_dfof_avg = squeeze(mean(data_dfof(:,:,:,nOff+1:nOff+nOn),4));
 
 dirs = unique(tGratingDirection);
 nDir = length(dirs);
 
 data_dfof_dir = zeros(sz(1),sz(2),nDir);
 
 figure;
 [n n2] = subplotn(nDir);
 for idir = 1:nDir
     ind = find(tGratingDirection == dirs(idir));
     data_dfof_dir(:,:,idir) = mean(data_dfof(:,:,ind),3);
     subplot(n,n2,idir)
     imagesc(data_dfof_dir(:,:,idir));
 end
 
 
 
     
     
 
     