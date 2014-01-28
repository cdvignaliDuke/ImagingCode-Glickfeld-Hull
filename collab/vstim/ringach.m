function imgs3 = ringach(prot, s)
%RINGACH
%IMGS3 = RINGACH(PROT, S)

% if ~isfield(prot.pars,'tf_0')
%     prot.pars.tf_0 = 0;
% end

%%
X = s.X(1:prot.DownsamplingFactor:end,1:prot.DownsamplingFactor:end);
Y = s.Y(1:prot.DownsamplingFactor:end,1:prot.DownsamplingFactor:end);
res = s.Resolution / prot.DownsamplingFactor;

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

refreshRate = s.RefreshRate/s.DecimationRatio;

dt = s.DecimationRatio/s.RefreshRate;
imgs3 = cell(1,prot.nstim);

fprintf(1,'RINGACH\n');

for istim = 1:prot.nstim    
    pars = prot.pars(istim);
    
    % create stimulus matrix
    sfs = linspace(pars.sf1,pars.sf2,pars.nsfs);
    oris = [0:pars.noris-1]*360/pars.noris;
    phs = [0:pars.nphs-1]*360/pars.nphs;
    
    nt = round(pars.non/pars.drift/dt);
    tt = [0:nt-1]*dt;
    wt = 2*pi*pars.drift*tt;
    
    imgs = cell(pars.nsfs,pars.noris,pars.nphs);
    
    for isf = 1:pars.nsfs
        sf = sfs(isf);
        for iori = 1:pars.noris
            ori = oris(iori);
            xrads = 2*pi*sf*cos(ori/180*pi)*X;
            yrads = 2*pi*sf*sin(ori/180*pi)*Y;     

            % add pi/2 keeps constant phase across directions
            rads = xrads + yrads + pi/2 ; 
                
            for iph = 1:pars.nphs
                ph = phs(iph);
                
                imgs{isf,iori,iph}= zeros([res nt],'uint8');   
                
                for it = 1:nt
                    angles = rads+ph/180*pi+wt(it);
                    
                    % interval [-1,1] as double;
                    stim = pars.c*sin(angles);

                    if squareWaveFlag
                         stim(find(stim>0)) = 1;
                         stim(find(stim==0)) = 0;
                         stim(find(stim<0)) = -1;
                    end

                    imgs{isf,iori,iph}(:,:,it)= round((stim+1)*127.5);            
                end % it
            end % iph
        end % iori
    end %isf
    
    siz = size(imgs);
    imgs = imgs(:);
    
    nb = round(pars.noff/pars.drift/dt);
    imgs{end+1}=ones([res nb],'uint8')*128;
    
    blank_ind = length(imgs);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%% begin code affecting state of random generator 
    
    % create stimulus sequence 
    rand('twister',pars.seed);    
    %seq = randperm(prod(size(imgs)));
    
    % quasi random sequences cover more uniformily 
    % than pseudo random ones, particularly important 
    % for short sequences 
    
    x = lhsdesign(pars.seqlen,3); % interval [0,1]
    x = ceil(bsxfun(@times,x,siz));
    
    %[length(find(diff(x(:,1))==0))     length(find(diff(x(:,2))==0))]
    
    % randomly permute successive stimuli with same sf or ori
    sel = find(diff(sfs(x(:,1))) == 0 ...
             | diff(oris(x(:,2)))== 0 ...
             | abs(diff(oris(x(:,2))))== 180);
    iter = 1;
    while length(sel)>1 & iter < 100
        p = randperm(length(sel));
        x(sel(p),:)= x(sel,:);
        %[length(find(diff(x(:,1))==0))     length(find(diff(x(:,2))==0))]
        sel = find(diff(sfs(x(:,1))) == 0 ...
                 | diff(oris(x(:,2)))== 0 ...
                 | abs(diff(oris(x(:,2))))== 180);
        iter = iter + 1;
    end
    
    if iter >= 100; warning('max scrambling iterations reached');end;
    
    % %%% end code affecting state random generator 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    seq = sub2ind(siz,x(:,1),x(:,2),x(:,3));
    
    % interleave blank stimuli
    if pars.inter > 0
        pad_len = ceil(pars.seqlen/pars.inter)*pars.inter;
        seq(end+1:pad_len) = blank_ind;
        seq2 = cat(1,reshape(seq,[pars.inter pad_len/pars.inter]),ones(1,pad_len/pars.inter)*blank_ind);
        seq2 = seq2(:);
    else
        seq2 = seq(:);
    end
    
    imgs2 = cell(1,1,length(seq2(:)));
    imgs2(:) = imgs(seq2(:));
    imgs3{istim} = cell2mat(imgs2);
    siz2 = size(imgs3{istim});
    
    if pars.tf_0 > 0
        b = fir1(21,pars.tf_0/(refreshRate/2));
        c = reshape(imgs3{istim},[siz2(1)*siz2(2),siz2(3)]);
        c = imfilter(c,b);
        c = reshape(c,[siz2(1),siz2(2),siz2(3)]);
        imgs3{istim} = c;
    end

    % last frame a blank per definition
    imgs3{istim}(:,:,end)=128;
    
    fprintf(' Stim %i, %i sfs %i oris %i phs = %i stims. %i stim -> %i pictures\n',...
              istim,pars.nsfs,pars.noris,pars.nphs,...
              prod([pars.nsfs,pars.noris,pars.nphs]),...
              pars.seqlen,siz2(3));
    
end % istim

%%
return;

%%
prot = readprot('I:\users\vincent\stimulation\prots\stim090723\ringach.prot');

