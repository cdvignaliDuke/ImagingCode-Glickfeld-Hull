function imgs = domakestims(protdir,stimdir,scrn,protname);
%DOMAKESTIMS
% DOMAKESTIMS(PROTDIR,STIMDIR,SCRN);
% DOMAKESTIMS(PROTDIR,STIMDIR,SCRN,PROTNAME);

if nargin < 4
    list = dir(fullfile(protdir,'*.prot'));
    
    fprintf('Found %i prot files\n',length(list));
    
    for ilist = 1:length(list);
        imgs = domakestims(protdir,stimdir,scrn,list(ilist).name);
    end;
    
    return;
end

[pathstr,name,ext]=fileparts(protname);

protfn= fullfile(protdir,protname);
targetpn = fullfile(stimdir,name);
imgs = makestims(protfn,scrn,targetpn);

return;