%get path names
date = '191120';
ImgFolder = strvcat('003');
run = strvcat('000');
mouse = 'i1306';
nrun = size(ImgFolder,1);
frame_rate = 15.5;
run_str = catRunName(ImgFolder, nrun);

%open TC data
CD = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P\' date '_' mouse '\' date '_' mouse '_' run_str];
cd(CD);
imgMatFile = [date '_' mouse '_' run_str '_TCs.mat'];
load(imgMatFile);
imgMatFile2 = [date '_' mouse '_' run_str '_input.mat'];
load(imgMatFile2);

%create necessary variables
tGratingDir = celleqel2mat_padded(input.tGratingDirectionDeg);
dirs = unique(tGratingDir);
nDir = length(dirs);
nOn = input.nScansOn;
nOff = input.nScansOff;
nFrames = nOn+nOff;
nCells = size(npSub_tc,2);
nTrials = size(tGratingDir,2);

%transform into dF/F for each trial
trial_tc = reshape(npSub_tc,[nFrames nTrials nCells]);
trial_f = mean(trial_tc(nOff/2:nOff,:,:),1);
trial_dfof = (trial_tc-trial_f)./trial_f;

%divide into trial types and plot timecourses
dfof_dir = zeros(nFrames,nCells,nDir);
tuning = zeros(nCells,nDir);
for iCell = 1:nCells
    [n n2] = subplotn(nDir);
    figure;
    for idir = 1:nDir
        subplot(2.*n, n2,idir)
        ind = find(tGratingDir==dirs(idir));
        dfof_dir(:,:,idir) = squeeze(mean(trial_dfof(:,ind,:),2));
        plot(dfof_dir(:,iCell,idir))
        title(num2str(dirs(idir)))
        ylabel('dF/F')
        xlabel('Frames')
    end 
    %plot average response for each trial type
    subplot(2,1,2)
    tuning(iCell,:) = squeeze(mean(dfof_dir(nOff+1:end,iCell,:),1));
    plot(dirs,tuning(iCell,:),'-o')
    xlabel('Direction (deg)')
    ylabel('dF/F')
    suptitle(['Cell #' num2str(iCell)]) 
    if exist(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P\' date '_' mouse '\' date '_' mouse '_' run_str '\directionTuning'])
        print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P\' date '_' mouse '\' date '_' mouse '_' run_str '\directionTuning\directionTuning_Cell' num2str(iCell)],'-dpdf','-fillpage') 

    else mkdir(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P\' date '_' mouse '\' date '_' mouse '_' run_str '\directionTuning'])
        print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P\' date '_' mouse '\' date '_' mouse '_' run_str '\directionTuning\directionTuning_Cell' num2str(iCell)],'-dpdf','-fillpage')
    end
end

    