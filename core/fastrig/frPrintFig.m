function frPrintFig(expt,figs,subdir)
% FRPRINTFIG(EXPT,FIGS,SUBDIR)

if isempty(figs)
    return;
end

if exist(expt.dirs.figsOut)~=7
    mkdir(expt.dirs.figsOut);
end

if nargin < 3
    pn = expt.dirs.figsOut;
else
    pn = fullfile(expt.dirs.figsOut,subdir);
    if exist(pn)~=7
        mkdir(pn);
    end
end

% % maximize docked figures container
% drawnow;
% desktop=com.mathworks.mde.desk.MLDesktop.getInstance;
% container=desktop.getGroupContainer('Figures').getTopLevelAncestor;
% container.setMaximized(1); % or 0 to return to normal
%  
for ifig = 1:length(figs)
    name = get(figs(ifig),'name');

    fn = fullfile(pn,[name '.jpeg']);
    
    print(figs(ifig),'-djpeg90',fn);
    fprintf('wrote %s\n',fn);
    
%     pn = fullfile(expt.dirs.figsOut,[name '.pdf']);
%     disp(pn);
%     print(figs(ifig),'-dpdf',pn);
    
end

% container.setMaximized(0); % or 0 to return to normal

return;