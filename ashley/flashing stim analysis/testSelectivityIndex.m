function [siTest95CI,siShuffTest] = testSelectivityIndex(...
    visData,audData)
    d_av = cat(1,visData,audData);
    nBoot = 1000;
    nc = size(d_av,2);
    n_vis = size(visData,1);
    n_aud = size(d_av,1) - n_vis;
    trueID = cat(2,ones(1,n_vis),zeros(1,n_aud));
    si_boot = nan(nc,nBoot);
    siShuff_boot = nan(nc,nBoot);
    for iboot = 1:nBoot
        sampInd = randsample(n_vis+n_aud,n_vis+n_aud,1);
        randID = trueID(randperm(n_vis+n_aud));
        trueInd_temp = trueID(sampInd);
        randInd_temp = randID(sampInd);
        d_av_temp = d_av(sampInd,:);
        vis_temp = d_av_temp(trueInd_temp==1,:);
        aud_temp = d_av_temp(trueInd_temp==0,:);
        si_boot(:,iboot) = getSelectivityIndex(vis_temp,aud_temp);
        vis_shuff_temp = d_av_temp(randInd_temp==1,:);
        aud_shuff_temp = d_av_temp(randInd_temp==0,:);
        siShuff_boot(:,iboot) = getSelectivityIndex(vis_shuff_temp,aud_shuff_temp);
    end

    si_bootPctInc = sum(si_boot>0,2)./nBoot;
    siTest95CI = si_bootPctInc > 0.95 | si_bootPctInc < 0.05;
    siShuff_bootPctInc = sum(siShuff_boot>0,2)./nBoot;
    siShuffTest = siShuff_bootPctInc > 0.95 | siShuff_bootPctInc < 0.05;
end

