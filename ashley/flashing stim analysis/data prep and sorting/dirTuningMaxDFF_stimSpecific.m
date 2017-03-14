dFF_crop = dFF;
dFF_crop(:,xcrop,:,:) = 0;
dFF_crop(ycrop,:,:,:) = 0;

dFF_dirmax = zeros(ypix,xpix,nstim);
for istim = 1:nstim
   dFF_dirmax(:,:,istim) = max(squeeze(mean(dFF_crop(:,:,off+1:off+on,tDir_ind == istim),3)),[],3); 
end

