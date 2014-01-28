
filelist='exp_list_scmr8_2.xls';
s=readXLSfilelist(filelist);

outDirectory = 'E:\users\kenichi\scmr8\pixelmaps';
InDirectory = '\\zmey\storlab\data\Soumya\scmr8\rtif'

nbinning=1;                 % spatial binning factor
lowpass = 2;
for i=76:length(s)
    s(i).nbinning=nbinning;
    s(i).lowpass=lowpass;
    psti=getPixMaps(fullfile(InDirectory, s(i).filename), [], outDirectory, s(i));
%    psti=getColorMaps(fullfile(InDirectory, s(i).fname), [], outDirectory, s(i));
end
