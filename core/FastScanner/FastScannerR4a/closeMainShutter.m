function closeMainShutter
global fs

fs.shutter.MainshutterIsOpen=0;
putvalue(fs.DAQ.MainshutterLine, fs.shutter.MainclosedLevel);