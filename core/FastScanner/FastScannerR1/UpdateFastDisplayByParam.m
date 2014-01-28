function UpdateFastDisplayByParam
global fd;

%fd.img.Ch(1).min
set(fd.handles.sliderMinChannel1,'Value',fd.img.Ch(1).min);
FastDisplay('sliderMinChannel1_Callback',fd.handles.sliderMinChannel1, [], fd.handles);

%fd.img.Ch(1).max
set(fd.handles.sliderMaxChannel1,'Value',fd.img.Ch(1).max);
FastDisplay('sliderMaxChannel1_Callback',fd.handles.sliderMaxChannel1, [], fd.handles);

%fd.img.Ch(2).min
set(fd.handles.sliderMinChannel2,'Value',fd.img.Ch(2).min);
FastDisplay('sliderMinChannel2_Callback',fd.handles.sliderMinChannel2, [], fd.handles);

%fd.img.Ch(2).max
set(fd.handles.sliderMaxChannel2,'Value',fd.img.Ch(2).max);
FastDisplay('sliderMaxChannel2_Callback',fd.handles.sliderMaxChannel2, [], fd.handles);

%fd.img.Ch(1).use
set(fd.handles.chkChannel1,'Value',fd.img.Ch(1).use);
FastDisplay('chkChannel1_Callback',fd.handles.chkChannel1, [], fd.handles);

%fd.img.Ch(2).use
set(fd.handles.chkChannel2,'Value',fd.img.Ch(2).use);
FastDisplay('chkChannel2_Callback',fd.handles.chkChannel2, [], fd.handles);

%fd.img.XSizePix
set(fd.handles.txtXSize,'String',fd.img.XSizePix);

%fd.img.YSizePix
set(fd.handles.txtYSize,'String',fd.img.YSizePix);

%fd.timing.FastMirrorTriggerDelay
set(fd.handles.txtFastMirrorTriggerDelay,'String',fd.timing.FastMirrorTriggerDelay);

%fd.img.AvgOnDisplay
set(fd.handles.txtAvgOnDisplay,'String',fd.img.AvgOnDisplay);



