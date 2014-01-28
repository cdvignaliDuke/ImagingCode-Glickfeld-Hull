
filelist='exp_list_scmr7_2_conecycle.xls';
s=readXLSfilelist(filelist);

outDirectory = 'E:\users\kenichi\conemaps\scmr7\';
InDirectory = 'E:\users\kenichi\scmr7\pixelmaps\'

for i=1:length(s)
    load(fullfile(InDirectory,s(i).filename,'dir_dF_sm.mat'));
    conemap=getConemap(dir_dF_sm);
    imwrite(conemap, fullfile(outDirectory, [s(i).filename, '.bmp']));
end
