function GenericKeyPress_FastDisplay
global fd;
[aaa,bbb]=gcbo;
%val = double(get(gcbo,'CurrentCharacter'));
val = double(get(bbb,'CurrentCharacter'));
if isempty(val)
    return;
end
switch val
    case 111 %o change time delay by -0.5
        fd.timing.FastMirrorTriggerDelay=fd.timing.FastMirrorTriggerDelay-0.5;
        set(fd.handles.txtFastMirrorTriggerDelay,'String',fd.timing.FastMirrorTriggerDelay);
        FastDisplay('txtFastMirrorTriggerDelay_Callback',fd.handles.txtFastMirrorTriggerDelay,[],fd.handles)
    case 112 %p chnage time delay by +0.5
        fd.timing.FastMirrorTriggerDelay=fd.timing.FastMirrorTriggerDelay+0.5;
        set(fd.handles.txtFastMirrorTriggerDelay,'String',fd.timing.FastMirrorTriggerDelay);
        FastDisplay('txtFastMirrorTriggerDelay_Callback',fd.handles.txtFastMirrorTriggerDelay,[],fd.handles)
end
