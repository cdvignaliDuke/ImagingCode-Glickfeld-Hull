function imgs = makestims(protname,s,targetdir)
%MAKESTIMS
% MAKESTIMS(PROTNAME,S,TARGETDIR)
% IMGS = MAKESTIMS(PROTNAME,S,TARGETDIR)

if nargin > 2
    if ~exist(targetdir)
        mkdir(targetdir);
    end
end

if isstr(protname) 
    [parthstr,name,ext]=fileparts(protname);
    if isempty(ext)
        protname = [protname '.prot'];
    end
else
    error('PROTNAME should be a valid path');
end

if exist(protname)
    [pathstr,name,ext]=fileparts(protname);    
    copyfile(protname,fullfile(targetdir,[name ext]));
    fprintf(1,'Loaded %s\n',protname);
    protname = readprot(protname);
else
    error(sprintf('protname file name %s does not exist',protname));
end

fprintf(1,'Generating stimuli...\n');

mystimfunc = str2func(protname.StimulusType);
imgs = mystimfunc(protname, s);

if nargin > 2
    writestims(imgs,targetdir);
end

return;

