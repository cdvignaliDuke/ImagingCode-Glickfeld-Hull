% batch convert all MPDs

clear all

MPdir='\\zstorlab\v1msq\data\Soumya\scmr2'
target_dir='\\zstorlab2\v1msq\users\Soumya\data\scmr2'
mkdir(target_dir, 'tif');
[filelist, nFiles] = GetMPFileNames(MPdir);
for i=1:nFiles
    filelist{i}
    params=ConvertMP2tif (filelist{i},MPdir,[target_dir,'\tif']);
end

clear all
