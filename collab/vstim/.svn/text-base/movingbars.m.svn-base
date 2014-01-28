function imgs = movingbars(prot, s)
% imgs = gratings(prot , s)

% dur * 

tic;

if isfield(prot,'DecimationRatio');
    s.DecimationRatio = prot.DecimationRatio;
end

if ~isfield(prot,'WindowType');
    prot.WindowType = 'rect';
end

if ~isfield(prot,'BarType');
    prot.GratingType = 'Gaussian';
end
    
if strcmp(prot.BarType,'SquareWave')
   squareWaveFlag = true;
else
   squareWaveFlag = false;
end

frameDur = s.DecimationRatio/s.RefreshRate;

total = 0;

xrange = [min(s.xdeg) max(s.xdeg)];
yrange = [min(s.ydeg) max(s.ydeg)];

xyrange =  sqrt(xrange.^2+yrange.^2);xyrange(1)=-xyrange(1);

% each protocol consists of multiple stimuli
for istim = 1:prot.nstim    
    sel = find([prot.pars.n]==istim);
    stimlen = sum(floor([prot.pars(sel).dur]/frameDur));
    imgs{istim} = zeros([s.Resolution stimlen],'uint8');   
    
    index = 1;
    % each stimulus consists of multiple contiguous segments    
    for iseg= 1:length(sel);
        pars = prot.pars(sel(iseg));
        disp(pars);

        nframes = round(pars.dur/frameDur);
        
        % bar positions along axis of motion
        pos = linspace(- nframes*frameDur * pars.drift/2,...
                         nframes*frameDur * pars.drift/2,nframes);

        % project all points in visual space onto axis of motion
        projection = cos(pars.ori/180*pi)*(s.X-pars.xc) + ...
                     sin(pars.ori/180*pi)*(s.Y-pars.yc);     

        mask = win(pars,s,prot.WindowType);
        
        for iframe = 1:nframes            
           % interval [-1,1]
            stim = pars.c*exp(-(projection+pos(iframe)).^2/2/pars.sigma^2);
            if squareWaveFlag
             stim = round(stim);
            end
            
            % interval [0,255]
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

