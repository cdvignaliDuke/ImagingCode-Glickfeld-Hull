exp = '160606_img46\';
FolderName = ['\\crash.dhe.duke.edu\data\home\jake\Data\WidefieldImaging\GCaMP\',exp];
cd(FolderName);
D = dir(FolderName);
Num = length(D(not([D.isdir])));

for i = 1:Num
    if i == 1
        FileTif{i}='160606_img46_MMStack.ome.tif';
        
    else
        FileTif{i}=['160606_img46_MMStack_' ,num2str(i-1) ,'.ome.tif'];
    end
    
    InfoImage=imfinfo(FileTif{i});
    mImage=InfoImage(1).Width;
    nImage=InfoImage(1).Height;
    NumberImages(i)=length(InfoImage);
end

FinalImage=zeros(nImage,mImage,sum(NumberImages),'uint16');

for fi = 2:Num
    TifLink = Tiff(FileTif{fi}, 'r');
    for i=1:NumberImages(fi)
        TifLink.setDirectory(i);
        if fi > 1
            FinalImage(:,:,i+NumberImages(fi-1))=TifLink.read();
        else
            FinalImage(:,:,i)=TifLink.read();
        end
    end
    TifLink.close();
end