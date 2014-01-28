function exptnames=frGetAllExpts(newdir)
% EXPTNAMES=FRGETALLEXPTS(NEWDIR)

dirs = frGetDirs;

exptdir = fullfile(dirs.images,newdir);

if exist(exptdir) ~= 7
    exptdir = fullfile(dirs.backup,'images',newdir);
end

if exist(exptdir) ~= 7    
    str = sprintf('directory ''%s'' does not exist',exptdir);
    error(str);
end

%% recursively process children directories
list = dir(exptdir);list(1:2)=[]; % deletes current and parent dir entry
sel = find([list.isdir]);

index = 1;
for idir = sel
    exptnames{index}=sprintf('%s\\%s',newdir,list(idir).name);
    index = index + 1;
end

return;