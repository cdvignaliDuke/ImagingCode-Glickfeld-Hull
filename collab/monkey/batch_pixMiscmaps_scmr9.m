
filelist='exp_list_scmr9_2_misc.xls';
s=readXLSfilelist(filelist);

outDirectory = 'E:\users\kenichi\scmr9\pixelMiscmaps';
InDirectory = '\\zmey\storlab\data\Soumya\scmr9\rtif'
pstiRootDirectory =  'E:\users\kenichi\scmr9\pixelmaps'
nbinning=1;                 % spatial binning factor
lowpass = 2;
for i=1:length(s)
    pstiDirectory=[pstiRootDirectory,'\',s(i).filename];
    s(i).nbinning=nbinning;
    s(i).lowpass=lowpass;1
    load(fullfile(pstiDirectory, 'psti.mat'));
    getPixOriMaps(fullfile(InDirectory, s(i).filename), psti, outDirectory, s(i));
%    psti=getColorMaps(fullfile(InDirectory, s(i).fname), [], outDirectory, s(i));
end
