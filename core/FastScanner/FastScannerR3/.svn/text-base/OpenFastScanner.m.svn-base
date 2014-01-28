function p=OpenFastScanner
%function called on activation of GUI

global fs;

IniParameters('FastScanner');

%ReadIniForFS('C:\projects\fastscanner\FastScanner.ini','FastScanner');
ReadIniForFS('FastScanner.ini','FastScanner');

SetupDAQ;

fs.iniDone=1;

fs.pCell.LUT=ReadPockelsVoltageLUT(fs.pCell.LUTfile);

UpdateFastScannerByParam;
SetPCell;