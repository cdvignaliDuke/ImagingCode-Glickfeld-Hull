

%bookkeeping
data= readtiff('C:\Users\lindsey\Desktop\ratV1_run5.tif');
data= data(:,:,1:4240);
siz = size(data);
nfon = 20;
nfoff = 60;
nftot = nfon+nfoff;
nrep = siz(3)./(nfon+nfoff);
nstim = 1;

%extract and avg trials
data_rep = zeros(siz(1),siz(2),nftot,nrep,nstim);
for irep = 1:nrep
    for istim = 1:nstim
        data_rep(:,:,:,irep,istim) = data(:,:,1+((nftot).*(irep-1)):(nftot).*irep);
    end
end

data_avg = zeros(siz(1),siz(2),nftot,nstim);
for istim = 1:nstim
    data_avg(:,:,:,istim) = squeeze(mean(data_rep(:,:,:,:,istim),4));
end
%flipping data so on comes before off
data_avg = cat(3,data_avg(:,:,1+nfoff:end,:),data_avg(:,:,1:nfoff,:));

%calc df/f
for istim = 1:nstim
    data_F = squeeze(mean(data_avg(:,:,nftot-9:end,istim),3));
    data_dF(:,:,:,istim) = bsxfun(@minus,squeeze(data_avg(:,:,:,istim)),data_F);
    data_dFoverF(:,:,:,istim) = bsxfun(@rdivide,squeeze(data_dF(:,:,:,istim)), data_F);
    data_dFoverF_down(:,:,:,istim) = stackGroupProject(squeeze(data_dFoverF(:,:,:,istim)),5);
end

%plot response timecourse
for istim = 1:nstim
    figure; 
    for i = 1:size(data_dFoverF_down,3)
        subplot(4,4,i)
        imagesc(data_dFoverF_down(:,:,i,istim));
        colormap(gray)
        axis off
        title(num2str(i*5));
    end
end

%find averageing windows
for istim = 1:nstim
    figure;
    n1 = [6 11 16];
    n2 = [30 35 40];
    start =1;
    for is=1:3
        for ie=1:3
            subplot(3,3,start)
            imagesc(squeeze(mean(data_dFoverF(:,:,n1(is):n2(ie),istim),3)));
            colormap(gray)
            axis off
            title(num2str([n1(is) n2(ie)]))
            start = 1+start;
        end
    end
end