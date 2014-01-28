function closeShutter
global fs


fs.shutter.shutterIsOpen=0;
putvalue(fs.DAQ.shutterLine, fs.shutter.closedLevel);