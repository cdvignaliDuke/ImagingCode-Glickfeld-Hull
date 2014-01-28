
filelist='exp_list_scmr6_2_conecycle.xls';
s=readXLSfilelist(filelist);

outDirectory = 'E:\users\kenichi\conemaps\scmr6\';
InDirectory = 'E:\users\kenichi\scmr6\pixelmaps\'

for i=1:length(s)
    load(fullfile(InDirectory,s(i).filename,'dir_dF_sm.mat'));
    conemap=getConemap(dir_dF_sm);
    imwrite(conemap, fullfile(outDirectory, [s(i).filename, '.bmp']));
end
