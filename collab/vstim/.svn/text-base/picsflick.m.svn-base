function imgs = picsflick(prot, s)

tic;

if isfield(prot,'DecimationRatio');
    s.DecimationRatio = prot.DecimationRatio;
end

if ~isfield(prot,'WindowType');
    prot.WindowType = 'rect';
end

dt = s.DecimationRatio/s.RefreshRate;

total = 0;

for istim = 1:prot.nstim    
    sel = find([prot.pars.n]==istim);
    stimlen = sum(floor([prot.pars(sel).dur]/dt));
    imgs{istim} = zeros([s.Resolution stimlen],'uint8');   
    
    index = 1;
    for iseg= 1:length(sel);
        pars = prot.pars(sel(iseg));
        
        nframes = floor(pars.dur/dt);
        tt = [0:nframes-1]*dt;
               
        amps = pars.c*cos(2*pi*pars.flick*tt);    
        flash = round((sin(2*pi*pars.flash*tt)+1)/2);
        mask = win(pars,s,prot.WindowType);
        
        fn = sprintf('pic%03d.bmp',pars.pic);
        pic = (double(imread(fullfile(prot.SourceDir,fn)))-128)/128.0;
        
        if ndims(pic)>2
            pic(:,:,2:3)=[];
        end
        
        for iframe = 1:nframes
            stim = flash(iframe)*amps(iframe)*pic.*mask;
            imgs{istim}(:,:,index)= round((stim+1)*127.5);            
            index = index + 1;
        end % iframe
        
    end % iseg
    
    total = total + stimlen;
    
end % istim

t = toc;
s = whos('imgs');

fprintf(1,'Generated %i grating bitmaps (%2.1f MB/s)\n',total,s.bytes/1024/1024/t);

return;


