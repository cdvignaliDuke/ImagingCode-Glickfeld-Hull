date = '180914';
mouse = 'i1103';
ImgFolder = {'001','002','003'};
time = {'1221','1244','1252'};
doReg = 0;
nrun = size(ImgFolder,2);
rc = behavConstsAV;

for irun = 3
    CD = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Data\2P_images\' date '_' mouse '\' ImgFolder{irun}];
    cd(CD);
    imgMatFile = [ImgFolder{irun} '_000_000.mat'];
    load(imgMatFile);

    nframes = info.config.frames;
    data = sbxread([ImgFolder{irun} '_000_000'],0,nframes);
    fprintf(['Loaded ' num2str(nframes) ' frames \r\n'])
    
    fName = ['\\CRASH.dhe.duke.edu\data\home\andrew\Behavior\Data\data-' mouse '-' date '-' time{irun} '.mat'];
    load(fName);
    
    nOn = input.nScansOn;
    nOff = input.nScansOff;
    ntrials = size(input.tGratingDirectionDeg,2);
    
    if size(data,1) == 2
        data = data(1,:,:,:);
    end
    data = squeeze(data);
  
    if doReg
    data_avg = mean(data(:,:,1000:1500),3);
    [out data_reg] = stackRegister(data,data_avg);
    data = data_reg;
    clear data_reg
    end
    
    sz = size(data);
    data = data(:,:,1:(nOn+nOff)*ntrials);
    data_tr = reshape(data,[sz(1), sz(2), nOn+nOff, ntrials]);
    data_f = mean(data_tr(:,:,nOff/2:nOff,:),3);
    data_df = bsxfun(@minus, double(permute(data_tr,[1 2 4 3])), squeeze(data_f)); 
    data_dfof = permute(bsxfun(@rdivide,data_df, squeeze(data_f)),[1 2 4 3]); 
    clear data_f data_df data_tr

    Az = celleqel2mat_padded(input.tGratingAzimuthDeg);
    El = celleqel2mat_padded(input.tGratingElevationDeg);
    Azs = unique(Az);
    Els = unique(El);
    if min(Els,[],2)<0
        Els = fliplr(Els);
    end
    nStim = length(Azs).*length(Els);
    Stims = [];
    data_dfof_avg = zeros(sz(1), sz(2), nOn+nOff, length(Azs).*length(Els));
    start = 1;
    for iEl = 1:length(Els)
        ind1 = find(El == Els(iEl));
        for iAz = 1:length(Azs)
            Stims = [Stims; Els(iEl) Azs(iAz)];
            ind2 = find(Az == Azs(iAz));
            ind = intersect(ind1,ind2);
            data_dfof_avg(:,:,:,start) = mean(data_dfof(:,:,:,ind),4);
            start = start +1;
        end
    end
    clear data_dfof
    myfilter = fspecial('gaussian',[20 20], 0.5);
    data_dfof_avg_all = squeeze(mean(imfilter(data_dfof_avg(:,:,nOff:nOff+nOn,:),myfilter),3));
    
    figure; 
    a = data_dfof_avg_all(:,:,[1 4 7]);
    b = data_dfof_avg_all(:,:,[2 5 8]);
    c = data_dfof_avg_all(:,:,[3 6 9]);
    subplot(3,1,1)
    image(a)
    title([mouse ' ' date ' Run ' num2str(ImgFolder{irun}) ' retinotopy'])
    subplot(3,1,2)
    image(b)
    subplot(3,1,3)
    image(c)
    print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\' date '_' mouse '\' date '_' mouse '_' ImgFolder{irun} '\' date '_' mouse '_' ImgFolder{irun} '_retRGB_all.pdf'], '-dpdf','-fillpage')

    figure;
    a = sum(data_dfof_avg_all(:,:,[1 4 7]),3);
    b = sum(data_dfof_avg_all(:,:,[2 5 8]),3);
    c = sum(data_dfof_avg_all(:,:,[3 6 9]),3);
    d = cat(3,a,b,c);
    image(d)
    truesize
    title([mouse ' ' date ' Run ' num2str(ImgFolder{irun}) ' retinotopy- Sum of El'])
    print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\' date '_' mouse '\' date '_' mouse '_' ImgFolder{irun} '\' date '_' mouse '_' ImgFolder{irun} '_retRGB_sum.pdf'], '-dpdf','-bestfit')

    figure
    a = max(data_dfof_avg_all(:,:,[1 4 7]),[],3);
    b = max(data_dfof_avg_all(:,:,[2 5 8]),[],3);
    c = max(data_dfof_avg_all(:,:,[3 6 9]),[],3);
    d = cat(3,a,b,c);
    image(d)
    truesize
    title([mouse ' ' date ' Run ' num2str(ImgFolder{irun}) ' retinotopy- Max of El'])
    
    print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\' date '_' mouse '\' date '_' mouse '_' ImgFolder{irun} '\' date '_' mouse '_' ImgFolder{irun} '_retRGB_max.pdf'], '-dpdf','-bestfit')

    figure;
    d2 = imfilter(imadjust(d,[0 0 0; .4 .4 .4],[]),myfilter);
    image(d2)
    truesize
    title([mouse ' ' date ' Run ' num2str(ImgFolder{irun}) ' retinotopy- Max of El'])
    print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\' date '_' mouse '\' date '_' mouse '_' ImgFolder{irun} '\' date '_' mouse '_' ImgFolder{irun} '_retRGB_conFilt.pdf'], '-dpdf','-bestfit')
end
