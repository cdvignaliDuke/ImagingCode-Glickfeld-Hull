function GenericKeyPress_FastScanner
global fs;
[aaa,bbb]=gcbo;
val = double(get(bbb,'CurrentCharacter'));
if isempty(val)
    return;
end
switch val
    case 111 %o change time delay by -1
        fs.DAQ.delay=fs.DAQ.delay-1;
        set(fs.handles.txtDelay,'String',fs.DAQ.delay);  
    case 112 %p chnage time delay by +1
        fs.DAQ.delay=fs.DAQ.delay+1;
        set(fs.handles.txtDelay,'String',fs.DAQ.delay);  
end
