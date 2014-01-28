function out=checkDAD_directory(targetdir,startind, endind,nch,usech,threshold)
% out=checkDAD_directory(targetdir,startind, endind,nch,usech,threshold)

[basedir,exptdir]=fileparts(targetdir);

%threshold=0.01;
out = [];
format compact
lastpoints=[];
for ind=startind:endind
    %filename = fullfile(targetdir, [exptdir num2str(ind) '.dad']);
    filename = makeDadName(fullfile(targetdir,exptdir), ind);
    filename
    fid= fopen(filename,'r');

    if fid<0
        display(['Unable to open ' filename]);
        return
    end
    
    [out, count]=fread(fid,inf,'uint16');
    fclose(fid);
    out=reshape(out,nch,count/nch);

    chkcann=[lastpoints out(usech,:)];

    out=[];
    lastpoints=chkcann(end-500:end);
    
    ma=max(chkcann(1:3000));
    mi=min(chkcann(1:3000));
    
    chkcann=(chkcann-mi)/(ma-mi);
    
    d=diff(chkcann);

    pind=find(abs(d)> threshold);
    if ~isempty(pind)
        pind
        figure
        plot(chkcann(pind(1)-500:pind(1)+500));
    end
    chkcann=[];    
    d=[];
end

return;
