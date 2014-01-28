function [filelist, nFiles] = GetFileList(dirs, type)

switch (type)
    case ('data')
        [filelist, nFiles] = GetTifFileNames(dirs.tif);
    case ('rdata')
        [filelist, nFiles] = GetTifFileNames(dirs.rtif,[],[],'avg');
    case ('labelimg')
        [filelist, nFiles] = GetMatFileNames(dirs.label,[],'labelimg');  
        for i=1:length(filelist)
            filelist{i}=GetBaseFname(filelist{i},'','_labelimg');
        end
    otherwise
        filelist=[];
        nFiles=0;
end
        
        
        
