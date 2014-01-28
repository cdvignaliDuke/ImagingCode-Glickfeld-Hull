
filelist='exp_list_scmr7_2.xls';
s=readXLSfilelist(filelist);

outDirectory = 'E:\users\kenichi\scmr7\pixelmaps';
InDirectory = '\\zmey\storlab\data\Soumya\scmr7\rtif'

nbinning=1;                 % spatial binning factor
lowpass = 2;
for i=37:38
%for i=41:length(s)
    s(i).nbinning=nbinning;
    s(i).lowpass=lowpass;
    psti=getPixMaps(fullfile(InDirectory, s(i).filename), [], outDirectory, s(i));
%    psti=getColorMaps(fullfile(InDirectory, s(i).fname), [], outDirectory, s(i));
end
