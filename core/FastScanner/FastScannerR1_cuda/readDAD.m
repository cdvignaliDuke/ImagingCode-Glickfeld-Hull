function out=readDAD(filename,pstart,pend,nch,usech)



    fid= fopen(filename,'r');
fseek(fid, (pstart-1)*2*nch, 'bof');
    a=(fread(fid,[nch,pend-pstart+1],'int16'))';
out=a(:,usech);
    fclose(fid);
figure
plot(diff(out));
figure;
plot(out)