
filelist='exp_list_scmr4_2_conecycle.xls';
s=readXLSfilelist(filelist);

outDirectory = 'E:\users\kenichi\conemaps\scmr4\';
InDirectory = 'E:\users\kenichi\scmr4\pixelmaps\'

for i=1:length(s)
    load(fullfile(InDirectory,s(i).filename,'dir_dF_sm.mat'));
    conemap=getConemap(dir_dF_sm);
    imwrite(conemap, fullfile(outDirectory, [s(i).filename, '.bmp']));
end
