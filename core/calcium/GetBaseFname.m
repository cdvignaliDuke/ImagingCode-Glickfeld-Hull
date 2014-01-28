
function fname = GetBaseFname (infname, prefix, suffix)

% prefix: default='' (do not change). If you want to add 'r' '+r'. 
%   If you want to delete 'r', '-r'
% suffix = '', or '_tcourse' etc., which will be removed from the file
% names
%   5.13.2008. Kenichi Ohki

if nargin < 2
    prefix ='';
end

if nargin < 3
    suffix ='';
end

[path,fname,ext]=fileparts(infname);

if ~isempty(prefix)
    if prefix(1)=='+'
        prefix(1)=[];
        fname=[prefix,fname];
    else if prefix(1)=='-'
            prefix(1)=[];
            if ~strcmpi(prefix, fname(1:length(prefix)))
                display('prefix does not match');
                fname=[];
                return;
            end
            fname(1:length(prefix))=[];
        else 
            display('first character should be + or -');
            fname=[];
            return;
        end
    end
end

if ~isempty(suffix)
    l=findstr(fname,suffix);
    if isempty(l)
        display('cannot find suffix in input filename');
        fname=[];
        return;
    end
    fname=fname(1:l-1);
end

end


