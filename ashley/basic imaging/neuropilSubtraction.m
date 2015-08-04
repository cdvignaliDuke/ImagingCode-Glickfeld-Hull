% First, I used your imCellNeuropil to create neuropil masks- I chose a buffer of 4 
% pixels and a neuropil width of 6 pixels.  I played around with this a bit- the buffer 
% size makes very little difference for most cells- you just need a few pixels of space 
% to avoid subtracting too much cell signal. The width mostly affects the noise- 6 pixels allows 
% sufficient averaging to not be too noisy.
% 
% Then I tried to find the optimal weight of the neuropil to subtract.  
% I did this by maximizing the skewness of the timecourse- this is how Vincent did it in his 2011 paper, 
% and the logic is as follows: an ideal calcium signal should be very sparse with the majority of points 
% being around the baseline and a few going positive. Skewness of the distribution is a good measure of 
% this: if you shift some of those baseline points to being more positive (as you would if you had 
% neuropil contamination) you decrease skewness; conversely, if you subtract too much neuropil, then
% values will start to go negative, also decreasing the skewness.  So an ideal weight for the
% neuropil will create maximum skewness. 
% 
% Then I just subtracted the neuropil according to the weight for each cell.


buf = 4;
np = 6;
neuropil = imCellNeuropil(mask_cell,buf,np);
% for i = 1:nCells
%     npTC(:,i) = stackGetTimeCourses(data_reg,neuropil(:,:,i));
% end
npTC = stackGetTimeCourses(data_reg,neuropil);

%get weights by maximizing skew
ii= 0.01:0.01:1;
x = zeros(length(ii), nCells);
for i = 1:100
    x(i,:) = skewness(data_TC-tcRemoveDC(npTC*ii(i)));
end
[max_skew ind] =  max(x,[],1);
skew(buf,:) = max_skew;
np_w = 0.01*ind;
npSubTC = data_TC-bsxfun(@times,tcRemoveDC(npTC),np_w);

fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
cd(fileSave);
save('neuropil.mat','neuropil','npTC','npSubTC');
