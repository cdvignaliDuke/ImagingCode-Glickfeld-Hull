function [prot]= readprot(filename)
% [prot]= readprot(filename)

fid = fopen(filename,'r');

if fid==-1 
    error('file does not exist');
end

tline = 0;

while true
    line = fgetl(fid); 
    tline = tline + 1;
    div = find(line=='=');
    if ~length(div)
        break;
    end   
    array = textscan(line(1:div-1),'%s');
    name = array{1}{1};
    
    % this allows of white spaces inside strings
    value = line(div+1:end);
    value = fliplr(deblank(fliplr(deblank(value))));% deblank on both ends
    
    if isempty(str2num(value));
        prot.(name) = value;
    else
        prot.(name) = str2num(value);
    end
end

array=textscan(line,'%s');
parnames=array{1};
npars=length(parnames);

index = 1;
pars = zeros(1,npars);

while true
    line = fgetl(fid);tline = tline +1;
    if line ==-1 
        break;
    end

    % ignore empty lines, need to be made more robust
    if isempty(line) 
        continue;
    end
    array=textscan(line,'%f');
    if length(array{1})==npars
        pars(index,:)=array{1}';
    else
        msg = sprintf('READPROT: Parse error at line %i in file \n%s\n',tline,filename);
        error(msg);
    end
    index = index + 1;
end

fclose(fid);

[pline , npars] = size(pars);

for iline = 1:pline
    for ipar = 1:npars
        prot.pars(iline).(parnames{ipar}) = pars(iline,ipar);
    end
end

prot.nstim = max([prot.pars.n]);

return;