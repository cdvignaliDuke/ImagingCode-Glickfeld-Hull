
filelist='exp_list_scmr4_2.xls';
s=readXLSfilelist(filelist);

outDirectory = 'E:\users\kenichi\scmr4\pixelOrimaps';
InDirectory = '\\zmey\storlab\data\Soumya\scmr4\rtif'
pstiRootDirectory =  'E:\users\kenichi\scmr4\pixelmaps'
nbinning=1;                 % spatial binning factor
lowpass = 2;
for i=1:length(s)
    if isempty(findstr(s(i).stim_code,'Img_ori'))
        continue;
    end

    
    pstiDirectory=[pstiRootDirectory,'\',s(i).filename];
    s(i).nbinning=nbinning;
    s(i).lowpass=lowpass;1
    load(fullfile(pstiDirectory, 'psti.mat'));
    getPixOriMaps(fullfile(InDirectory, s(i).filename), psti, outDirectory, s(i));
%    psti=getColorMaps(fullfile(InDirectory, s(i).fname), [], outDirectory, s(i));
end
