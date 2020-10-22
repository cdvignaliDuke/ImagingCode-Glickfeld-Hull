clear;
filename = 'C:\Duke\Research\Lindsey\Cover\001\001_000_001.ephys';
fileID = fopen(filename, 'r', 'ieee-le');
if fileID == -1, error('Cannot open file: %s', filename); end
format = 'uint32';
Data = fread(fileID, Inf, format);
fclose(fileID);
figure(1); clf; plot(Data(1:2:end)); title('Clock');
figure(2); clf; plot(Data(2:2:end)); title('Photodiode');

%% Data Analysis
%Clock
clock = Data(1:2:end);
clockLoc = [];
clockSize = size(clock);
clockAvg = mean(clock);
clockStd = std(clock);
clockAvgStd = clockAvg-4*clockStd;

%Photodiode
photo = Data(2:2:end);
photoLoc = [];
photoSize = size(photo);

figure(1); clf; plot(clock); title('Clock'); hold on; hline(clockAvg); hline(clockAvg-5*clockStd);

for i = 2:clockSize(1)
    if(clock(i) < clockAvgStd && clock(i-1) > clockAvgStd)
        clockLoc = [clockLoc i];
    end
end

% j = 1;
% begin = 0;
% while j < photoSize(1)
%     if(photo(j) > 3.2075*10^9)
%         begin = j;
%         j = size(photo) + 1;
%     else
%         j = j + 1;
%     end
% end

photoAvg = mean(photo(clockLoc(1):photoSize(1)));
photoStd = std(photo(clockLoc(1):photoSize(1)));
photoAvgStd = photoAvg-photoStd;
figure(2); clf; plot(photo); title('Photodiode'); hline(photoAvg); hline(photoAvgStd);

for i = clockLoc(1):photoSize(1) %begin:size(photo)
    if(photo(i) < photoAvgStd && photo(i-1) > photoAvgStd && (i+400) < photoSize(1) && photo(i+400) < photoAvgStd)
        photoLoc = [photoLoc i];
    end
end

photoFrames = [];
photoLocSize = size(photoLoc);
clockLocSize = size(clockLoc);
for i = 1:photoLocSize(2)
    for j = 1:clockLocSize(2)
        if(photoLoc(i) > clockLoc(j) && photoLoc(i) < clockLoc(j+1))
            photoFrames = [photoFrames j];
        end
    end
end
figure(3); clf; plot(Data(1:2:end)); hold on; plot(Data(2:2:end), 'r-'); plot(photoLoc, photo_smooth(photoLoc), 'om');
