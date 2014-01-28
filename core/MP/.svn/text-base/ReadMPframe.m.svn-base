function frame = ReadMPframe (mpfile, frame_num, channel, params)
if frame_num < 1 | frame_num > params.FrameCount
    display('frame number is out of range');
    return;
end
frame=invoke(mpfile, 'ReadFrameData', channel,frame_num);
frame=reshape(frame,params.FrameHeight,params.FrameWidth);
frame=frame';
