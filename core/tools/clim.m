function cl = clim(arg1,arg2)
%CLIM
%CLIM(LIMS)
%CLIM(AX,LIMS)

if nargin ==0
    cl = get(gca,'clim');
elseif nargin ==1 
    if all(ishandle(arg1))
        cl = get(arg1,'clim');
    else 
        set(gca,'clim',arg1);
    end;
elseif nargin == 2
    set(arg1,'clim',arg2);
end
    
return;
   
