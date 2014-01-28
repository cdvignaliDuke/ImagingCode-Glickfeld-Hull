
filelist='exp_list_scmr9_2.xls';
s=readXLSfilelist(filelist);

outDirectory = 'E:\users\kenichi\scmr9\pixelmaps';
InDirectory = '\\zmey\storlab\data\Soumya\scmr9\rtif'

nbinning=1;                 % spatial binning factor
lowpass = 2;
for i=53:length(s)
    if ~isempty(findstr(s(i).stim_code,'hartley'))
        continue;
    end

    s(i).nbinning=nbinning;
    s(i).lowpass=lowpass;
    psti=getPixMaps(fullfile(InDirectory, s(i).filename), [], outDirectory, s(i));
%    psti=getColorMaps(fullfile(InDirectory, s(i).fname), [], outDirectory, s(i));
end
