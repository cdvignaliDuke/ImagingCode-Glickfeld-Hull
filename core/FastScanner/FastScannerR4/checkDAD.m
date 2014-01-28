function out=checkDAD(filename,nch,usech)



fid= fopen(filename,'r');
fseek(fid, 0, 'bof');
c=0;
%while 1
c=c+1;
    [a, N]=fread(fid,[nch,1000000],'uint16');
%    [a, N]=fread(fid,[1000000],'uint16');

%     figure
%     plot(a');
    if N==0
%        break
    end
    a=a';
out=a(:,usech);

d=diff(out);

% if any(abs(d)> 30)
%     x=1
%     c
% end




%end
    fclose(fid);
    
