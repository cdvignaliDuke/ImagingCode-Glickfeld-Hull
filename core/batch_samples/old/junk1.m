
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
data_dir='\\zmey\storlab\data\Soumya\scmr6\';

tif_dir=[data_dir,'tif\'];
rtif_dir=[data_dir,'rtif\'];
 

for i=1:Ngroups
    run_index=find(groups==i);
    if isempty(run_index)
        continue;
    end
    
    
    groupname=['group',num2str(i)];
    out_dir=[data_dir,groupname,'\rtif\'];

    for j=1:length(run_index)
        rfname=s(j).filename;
        infname=fullfile(out_dir, [rfname,'.tif']);
        outfname=fullfile(out_dir,[rfname,'_',groupname,'.tif']);
        rename(infname, outfname);
    end
end
clear all

