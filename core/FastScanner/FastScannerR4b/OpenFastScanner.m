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

%IniStage;  % mh 090130
if IniStage
    set(fs.handles.txtZPos,'Enable','on');
    set(fs.handles.sliderZPos,'Enable','on');
 %   set(fs.handles.txtStartZ,'Enable','on');
    set(fs.handles.txtStepZ,'Enable','on');
end