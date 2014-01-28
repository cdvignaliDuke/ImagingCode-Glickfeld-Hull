function openMainShutter
global fs


fs.shutter.MainshutterIsOpen=1;
putvalue(fs.DAQ.MainshutterLine, 1-fs.shutter.MainclosedLevel);