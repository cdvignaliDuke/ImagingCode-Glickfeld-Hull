function matchc(arg1)

if nargin < 1
    cl = clim;
    set(gca,'clim',max(abs([cl(:)]))*[-1 1]);
else
    cl = clim(arg1);
    set(arg1,'clim',max(abs([cl{:}]))*[-1 1]);
end

return
