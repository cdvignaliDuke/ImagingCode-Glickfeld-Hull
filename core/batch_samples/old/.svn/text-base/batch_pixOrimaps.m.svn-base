
filelist='exp_list_scmr6_2.xls';
s=readXLSfilelist(filelist);

outDirectory = 'E:\users\kenichi\scmr6\pixelOrimaps';
InDirectory = '\\zmey\storlab\data\Soumya\scmr6\rtif'
pstiRootDirectory =  '\\zmey\storlab\data\Soumya\scmr6\pixelmaps'
nbinning=1;                 % spatial binning factor
lowpass = 2;
for i=1:length(s)
    if isempty(findstr(s(i).stim_code,'ori'))
        continue;
    end

    
    pstiDirectory=[pstiRootDirectory,'\',s(i).filename];
    s(i).nbinning=nbinning;
    s(i).lowpass=lowpass;
    load(fullfile(pstiDirectory, 'psti.mat'));
    getPixOriMaps(fullfile(InDirectory, s(i).filename), psti, outDirectory, s(i));
%    psti=getColorMaps(fullfile(InDirectory, s(i).fname), [], outDirectory, s(i));
end
