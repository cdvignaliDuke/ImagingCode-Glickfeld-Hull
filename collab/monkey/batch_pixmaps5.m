
filelist='exp_list_scmr5_2.xls';
s=readXLSfilelist(filelist);

outDirectory = 'E:\users\kenichi\scmr5\pixelmaps';
InDirectory = '\\zmey\storlab\data\Soumya\scmr5\rtif'

nbinning=1;                 % spatial binning factor
lowpass = 2;
for i=1:length(s)
    s(i).nbinning=nbinning;
    s(i).lowpass=lowpass;
    psti=getPixMaps(fullfile(InDirectory, s(i).filename), [], outDirectory, s(i));
%    psti=getColorMaps(fullfile(InDirectory, s(i).fname), [], outDirectory, s(i));
end
