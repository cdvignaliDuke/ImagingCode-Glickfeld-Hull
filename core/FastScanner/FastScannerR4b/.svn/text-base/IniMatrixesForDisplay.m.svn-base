function p=IniMatrixesForDisplay
p=0;

global fd;


if fd.img.Ch(1).use
    set(fd.img.Ch(1).ax,'XLim',[1 fd.datasize2]);
    set(fd.img.Ch(1).ax,'YLim',[1 fd.datasize1]);
    set(fd.img.Ch(1).imH, 'EraseMode', 'none', 'CDataMapping','scaled');
end
if fd.img.Ch(2).use
    set(fd.img.Ch(2).ax,'XLim',[1 fd.datasize2]);
    set(fd.img.Ch(2).ax,'YLim',[1 fd.datasize1]);
    set(fd.img.Ch(2).imH, 'EraseMode', 'none', 'CDataMapping','scaled');
end

fd.img.Ch(1).CData =zeros(fd.datasize1,fd.datasize2);
fd.img.Ch(2).CData =zeros(fd.datasize1,fd.datasize2);

fd.img.Ch(1).CDataV2 =zeros(fd.datasize1*fd.datasize2,1);
fd.img.Ch(2).CDataV2 =zeros(fd.datasize1*fd.datasize2,1);

fd.dataShow0=zeros(fd.datasize2,fd.datasize1,3);
fd.dataShow=zeros(fd.datasize1,fd.datasize2,3);

fd.img.Ch(1).CData2D =zeros(fd.datasize1*fd.datasize2,fd.img.AvgOnDisplay);
fd.img.Ch(2).CData2D =zeros(fd.datasize1*fd.datasize2,fd.img.AvgOnDisplay);



