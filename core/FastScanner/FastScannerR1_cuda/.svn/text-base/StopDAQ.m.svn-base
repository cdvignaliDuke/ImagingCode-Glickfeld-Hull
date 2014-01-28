function p = StopDAQ
p =0;

global fs;

stop(fs.DAQ.ai);
%stop(fs.DAQ.aoAcqTrigger);
%putsample(fs.DAQ.aoAcqTrigger,[5]); % receiver triggers on negative edge

stop(fs.DAQ.aoPockels);
if fs.DAQ.BlockEnds==1
    set(fs.handles.sliderPCellMin,'Enable','on');
    set(fs.handles.sliderPCellMax,'Enable','on');
    set(fs.handles.txtPCellMin,'Enable','on');
    set(fs.handles.txtPCellMax,'Enable','on');
end

% for estim
fs.DAQ.aoCurrentTriggerValue=0;

SetPCell;

% if ~fs.DAQ.Running
%     stop(fs.DAQ.aoExtTrigger);
%     putsample(fs.DAQ.aoExtTrigger,[5]); % receiver triggers on negative edge
% end

p = 1;

fs.DAQ.started = 0;

return