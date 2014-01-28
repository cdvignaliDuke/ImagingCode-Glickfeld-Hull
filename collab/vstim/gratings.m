function imgs = gratings(prot, s)
%GRATINGS 
% IMGS = GRATINGS(PROT , S)
%

% change log
%
% 09/07/23 vb fixed white blanks bug. watch out. phase of square wave
% gratings has changed.
% 

tic;

if isfield(prot,'DecimationRatio');
    s.DecimationRatio = prot.DecimationRatio;
end

if ~isfield(prot,'WindowType');
    prot.WindowType = 'rect';
end

if ~isfield(prot,'GratingType');
    prot.GratingType = 'Sinusoidal';
end
    
if strcmp(prot.GratingType,'SquareWave')
   squareWaveFlag = true;
else
   squareWaveFlag = false;
end

dt = s.DecimationRatio/s.RefreshRate;

total = 0;

% each protocol consists of multiple stimuli
for istim = 1:prot.nstim    
    sel = find([prot.pars.n]==istim);
    stimlen = sum(floor([prot.pars(sel).dur]/dt));
    imgs{istim} = zeros([s.Resolution stimlen],'uint8');   
    
    index = 1;
    % each stimulus consists of multiple contiguous segments    
    for iseg= 1:length(sel);
        pars = prot.pars(sel(iseg));
        disp(pars);

        nframes = round(pars.dur/dt);
        tt = [0:nframes-1]*dt;
               
        amps = pars.c*cos(2*pi*pars.flick*tt);
        phs = 2*pi*pars.drift*tt;
        
        % compute angular frequency
        % add pi/2 keeps constant phase across directions
        xrads = 2*pi*pars.sf*cos(pars.ori/180*pi)*(s.X-pars.xc);
        yrads = 2*pi*pars.sf*sin(pars.ori/180*pi)*(s.Y-pars.yc);     
        rads = xrads + yrads + pi/2 + pars.ph/180*pi; 

        mask = win(pars,s,prot.WindowType);
        
        
        for iframe = 1:nframes
            angles = rads+phs(iframe);
            % interval [-1,1] as double;
            stim = amps(iframe)*sin(angles).*mask;
            
            if squareWaveFlag
                 stim(find(stim>0)) = 1;
                 stim(find(stim==0)) = 0;
                 stim(find(stim<0)) = -1;
                 
%                stim = round(stim);
            end
            
            imgs{istim}(:,:,index)= round((stim+1)*127.5);            
            index = index + 1;
        end % iframe
%        fprintf(1,'.');
    end % iseg
    total = total + stimlen;
    
end % istim

t = toc;
s = whos('imgs');

fprintf(1,'Generated %i grating bitmaps (%2.1f MB/s)\n',total,s.bytes/1024/1024/t);

return;

%% test
s=screen;
x0 = 0; y0=0;c = 1.0; ori= 0; sf = 0.05;ph = 0;tf = 4;
pars = drift([x0 y0 c ori sf ph],s,tf,1,1);
imgs = gratings(pars, s);
playstim(imgs,5,60)
