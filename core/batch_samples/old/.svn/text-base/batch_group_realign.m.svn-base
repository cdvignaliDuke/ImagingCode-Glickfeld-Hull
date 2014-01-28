
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% batch file to realign images
%   data_dir
%       |-tif   xxxx.tif (original stack)
%       |-rtif  rxxxx.tif (realigned stack)
%
% start a few slaves (startRegSlave.m) before running this batch
%
%       Kenichi Ohki  04.28.08


filelist='exp_list_scmr6_3.xls';
s=readXLSfilelist(filelist);

groups=[s.group];
Ngroups=max(groups);

mcoreDIR='E:\users\kenichi\multcore';
% set this to your multicore directory (avoid to use my folder!)

overwrite=0; % do not overwrite

data_dir='\\zmey\storlab\data\Soumya\scmr6\';

tif_dir=[data_dir,'tif\'];
rtif_dir=[data_dir,'rtif\'];
 
% [filelist, nFiles] = GetTifFileNames(tif_dir);


for i=15:Ngroups
%for i=1:Ngroups
    run_index=find(groups==i);
    if isempty(run_index)
        continue;
    end
    
    % read target image
    n=run_index(floor((length(run_index)+1)/2));
    filename=[s(n).filename,'_avg.tif'];
    target=readtiff(fullfile(rtif_dir,filename));
    
    groupname=['group',num2str(i)];

    for j=1:length(run_index)
        rfname=s(run_index(j)).filename;
        filename=s(run_index(j)).filename;
        filename=filename(2:length(filename)); % remove 'r'
        fname=fullfile(tif_dir, [filename,'.tif']);
        outfname=fullfile(rtif_dir, [rfname,'_',groupname]);
        if exist(outfname)==2 && overwrite~=1
             continue;
        end
        stackRealignTifStack(fname, outfname, mcoreDIR,target)
    end
end
clear all

