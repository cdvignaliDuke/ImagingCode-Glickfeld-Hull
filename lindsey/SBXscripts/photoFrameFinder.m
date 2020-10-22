function [photoLoc photoFrames] = photoFrameFinder(photoData);

    % uses ephys output from scanbox to find frames stim onsets
    
    %Clock
    clock = photoData(1:2:end);
    clockSize = size(clock);
    clockAvg = mean(clock);
    clockStd = std(clock);
    clockAvgStd = clockAvg-4*clockStd;

    %Photodiode
    photo = photoData(2:2:end);

    clockLoc = [];
    for i = 2:clockSize(1)
        if(clock(i) < clockAvgStd && clock(i-1) > clockAvgStd)
            clockLoc = [clockLoc i];
        end
    end

    photo_smooth = smooth(photo,500);
    photo_high = movmax(photo_smooth, 100000);
    photo_smooth = photo_smooth-photo_high;
    photoSize = size(photo_smooth);
    photoMax = max(photo_smooth(clockLoc(1):photoSize(1)),[],1);
    photoMin = min(photo_smooth(clockLoc(1):photoSize(1)),[],1);
    photoAvg50= photoMax-(0.5*(photoMax-photoMin));
    photoAvg25 = photoMax-(0.25*(photoMax-photoMin));
    photoAvg75 = photoMax-(0.75*(photoMax-photoMin));
    photoAvg1 = photoMax-(0.01*(photoMax-photoMin));

    photoLoc = [];
    photoFrames = [];
    n = 0;
    for i = 2:length(clockLoc)
        photoAmp_min = min(photo_smooth(clockLoc(i-1):clockLoc(i)),[],1);
        photoAmp_prev = photo_smooth(clockLoc(i-1));
        
        if n == 0 && photoAmp_min<photoAvg75 && photoAmp_prev>photoAvg75
            ind = find(photo_smooth(clockLoc(i):clockLoc(i+1))<photoAvg75);
            if length(ind)>100
                photoLoc = [photoLoc clockLoc(i)];
                photoFrames = [photoFrames i];
                n=1+n;
            end
        elseif n>0 && max(photo_smooth(photoLoc(n):clockLoc(i-1)),[],1) > photoAvg1 && i-photoFrames(n) > 4
            if photoAmp_min<photoAvg25 && photoAmp_prev>photoAvg25 
                photoLoc = [photoLoc clockLoc(i)];
                photoFrames = [photoFrames i];
                n=1+n;
            end
        end
    end
end