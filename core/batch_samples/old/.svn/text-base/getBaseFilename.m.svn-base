function outfname = getBaseFilename(infname, insuffix)

% get base of file name
% e.g. suffix = '', or '_tcourse' etc.
%
%   5.13.2008. Kenichi Ohki

[path,fname,ext]=fileparts(infname);

if ~isempty(insuffix)
    l=findstr(fname,insuffix);
    if isempty(l)
        display('cannot find suffix in input filename');
        outfname=fname;
        return;
    end
    outfname=fname(1:l-1);
end
