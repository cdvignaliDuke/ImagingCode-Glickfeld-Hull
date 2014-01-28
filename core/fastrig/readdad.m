function out = readdad(inputfn)
% out = readdad(inputfn)

fid= fopen(inputfn,'r');
if fid<0
    display(['Unable to open ' inputfn]);
end

% for some reason reading native uint16 from disk very slow
tic;
[out, count]=fread(fid,inf,'uint16=>int16');
t = toc;

fclose(fid);

fprintf(1,'Read %2.1f MB in %2.1f s = %2.1f MB/s\n',count*2/1024^2,t,2*count/1024^2/t);

out = reshape(out,4,count/4);

% out = out(chan,:);

return;
