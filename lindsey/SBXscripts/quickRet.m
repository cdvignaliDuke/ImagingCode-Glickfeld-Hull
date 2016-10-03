function quickRet(date,mouse,ImgFolder,time,subnum)

nrun = size(ImgFolder,1);
rc = behavConstsAV;

run_str = ['runs-' ImgFolder(1,:)];
if nrun>1
    run_str = [run_str '-' ImgFolder(nrun,:)];
end

data = [];
clear temp
for irun = 1:nrun
    if strcmp(rc.name,'ashle')
    CD = ['Y:\home\ashley\data\' mouse '\two-photon imaging\' date '\' ImgFolder(irun,:)];
    else
    CD = ['Z:\home\lindsey\Data\2P_images\' date '_' mouse '\' ImgFolder(irun,:)];
    end
    cd(CD);
    imgMatFile = [ImgFolder(irun,:) '_000_000.mat'];
    load(imgMatFile);

    nframes = info.config.frames;
    data_temp = sbxread([ImgFolder(irun,:) '_000_000'],0,nframes);
    
    
    if strcmp(rc.name,'ashle')        
    fName = ['\\CRASH.dhe.duke.edu\data\home\andrew\Behavior\Data\data-i' subnum '-' date '-' time(irun,:) '.mat'];
    else
    fName = ['\\CRASH.dhe.duke.edu\data\home\andrew\Behavior\Data\data-' mouse '-' date '-' time(irun,:) '.mat'];
    end
    load(fName);
    temp(irun) = input;
    
    nOn = temp(irun).nScansOn;
    nOff = temp(irun).nScansOff;
    ntrials = size(temp(irun).tGratingDirectionDeg,2);
    
    data_temp = squeeze(data_temp);
%     if nframes>ntrials*(nOn+nOff)
%         data_temp = data_temp(:,:,1:ntrials*(nOn+nOff));
%     elseif nframes<ntrials*(nOn+nOff)
%         temp(irun) = trialChopper(temp(irun),1:ceil(nframes./(nOn+nOff)));
%     end
        
    data = cat(3,data,data_temp);
end
clear data_temp
input = concatenateDataBlocks(temp);


    nOn = input.nScansOn;
    nOff = input.nScansOff;
    ntrials = size(input.tGratingDirectionDeg,2);
    
   %use if pmt 2 was saved
%     data = squeeze(data(1,:,:,:));
    
    sz = size(data);
    data_tr = reshape(data,[sz(1), sz(2), nOn+nOff, ntrials]);
    data_f = mean(data_tr(:,:,nOff/2:nOff,:),3);
    data_df = bsxfun(@minus, double(data_tr), data_f); 
    data_dfof = bsxfun(@rdivide,data_df, data_f); 
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

    img_avg_resp = zeros(1,nStim);
    figure; 
    for i = 1:nStim; 
        subplot(length(Els),length(Azs),i); 
        imagesc(data_dfof_avg_all(:,:,i)); 
        colormap(gray)
        title(num2str(Stims(i,:)))
        img_avg_resp(i) = mean(mean(mean(data_dfof_avg_all(:,:,i),3),2),1);
        %clim([0 max(data_dfof_avg_all(:))./2])
    end
    figure
    heatmap = imagesc(reshape(img_avg_resp,length(Els),length(Azs)));
    heatmap.Parent.YTick = 1:length(Els);
    heatmap.Parent.YTickLabel = strread(num2str(Els),'%s');    
    heatmap.Parent.XTick = 1:length(Azs);
    heatmap.Parent.XTickLabel = strread(num2str(Azs),'%s');
    xlabel('Azimuth');
    ylabel('Elevation');
    colorbar
    caxis([-0.1 0.1])
    
