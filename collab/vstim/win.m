function mask = win(pars,s,wintype)
% mask = win(pars,s,wintype)

if nargin < 3
    wintype = 'rect';
end

switch wintype
    
    case 'rect' 
        
        if pars.width > 0
            maskX = abs(s.X - pars.xc)<=pars.width/2.0;
        else
            maskX = ones(s.Resolution);
        end
        
        if pars.height > 0
            maskY = abs(s.Y - pars.yc)<=pars.height/2.0;        
        else
            maskY = ones(s.Resolution);
        end
        
        mask = maskX&maskY;
        
    case  'gaus'
    
        xx2 = 0.5 * (s.X - pars.xc).^2 / (pars.width / 2.0).^2;
        yy2 = 0.5 * (s.Y - pars.yc).^2 / (pars.height / 2.0).^2;                
        mask = exp(-xx2 - yy2);  % peak = 1.0
        
    otherwise
        error('Unknown Window Type');
end

return;