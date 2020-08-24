openNSx
%% check if the channel ID is correct for poly shanks 
sum(double(NS6.MetaTags.ChannelID(1:32,1))-(1:1:32)')
if iscell(NS6.Data)
   NS6.Data = cell2mat(NS6.Data);
end
if sum(double(NS6.MetaTags.ChannelID(1:32,1))-(1:1:32)')==0

data=NS6.Data(1:end-3,:);

% for disabled channels of shank probes 
else
IDs = NS6.MetaTags.ChannelID; 
data= zeros(32,size(NS6.Data,2)); 
for i = 1:(length(IDs)-3)
    data(IDs(i),:) = NS6.Data(i,:);
    
end 

end
% 
% % deal with 3 shanks
% IDs = NS6.MetaTags.ChannelID; 
% data= zeros(32,size(NS6.Data,2)); 
% for i = 1:(length(IDs)-3)
%     data(IDs(i),:) = NS6.Data(i,:);
%     
% end 





data = data(:);

newFilename = ['D:\Miao_temp\spyking-circus\190220\D006\Data\' 'D006.dat'];

% Opening the output file for saving
FIDw = fopen(newFilename, 'w+', 'ieee-le');

% Writing data into file
disp('Writing the converted data into the new .dat file...');
fwrite(FIDw, data, 'int16');
fclose(FIDw);