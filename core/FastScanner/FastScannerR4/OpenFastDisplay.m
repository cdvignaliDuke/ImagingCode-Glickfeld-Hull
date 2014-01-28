function p=OpenFastDisplay(varargin)
%function called on activation of GUI
p=0;
global fd;

IniParameters('FastDisplay');

%ReadIniForFS('c:\projects\fastscanner\FastDisplay.ini','FastDisplay');
ReadIniForFS('FastDisplay.ini','FastDisplay');

fd.iniDone=1;

UpdateFastDisplayByParam
p=1;