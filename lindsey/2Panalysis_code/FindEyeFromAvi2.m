function eyefit=FindEyeFromAvi2( aviName,ROI,subsampl, size_range_1, size_range_2,dilation_range, motion_range, contrast_tsh,smooth,squeeze_coeff,saveFrames,frameindA)
% find position of reflrction and pupl using fast radial transform for avi data

%aviName - name of avi file
% ROI - 4 element array specifing ROI ( [x1 y1 x2 y2])
%size_range_1 - range of reflection radii to test
%size_range_2 - range of pupl radii to test
%dilation_range - reasonble upper limit on pupl dilation in one frame time
%motion_range - reasonble upper limit on reflection and pupl motion in one frame time 
%contrast_tsh - treshold for image gradient (small gradient is ignored)
%smooth - smoothing parameter for postprocessing
%subsampl - image subsample 
%squeeze_coeff - compensation for image elongation, for horiztontaly elongated images squeeze_coeff<1
%saveFrames - if ==1 then save images with positions marked
%frameindA - frames to analize, (all if empty)




% eyefit=FindEyeFromAvi2('E:\Lindsey\110426\long\110426_LONG.avi',[221 1 460 240],2, 12:18, 35:60,1,5, 1,3,0.5,1,1:500);

ind_chunk=500; %for memory limited computers

if nargin<3
    subsampl=1;
end
if nargin<4 || isempty(size_range_1)
    size_range_1=5:20;
end

if nargin<5 || isempty(size_range_2)
    size_range_2=10:50;
end

if nargin<6 || isempty(dilation_range)
    dilation_range=1;
end

if nargin<7 || isempty(motion_range)
    motion_range=10;
end

if nargin<8 || isempty(contrast_tsh)
    contrast_tsh=0;
end

if nargin<9 || isempty(smooth)
    smooth=1;
end
if nargin<10
    squeeze_coeff=1;
end

if nargin<11
    saveFrames=0;
end

warning off all;
fileinfo = aviinfo(aviName);
warning on all;

if nargin<12 || isempty(frameindA)
    frameindA=1:fileinfo. NumFrames;
end


    % smoothing filter
    sp_filter=fspecial('gaussian', round(smooth*4), smooth);
    %make filter round - to make sure that filter does not create
    % any vertical/horizonatl artifacts
    f_sz=size(sp_filter);
    [meshX, meshY]=meshgrid(-(f_sz(1)-1)/2:(f_sz(1)-1)/2, -(f_sz(2)-1)/2:(f_sz(2)-1)/2);
    sp_filter((meshX.^2+meshY.^2)>(f_sz(1)/2).^2)=0;
    sp_filter=sp_filter/sum(sp_filter(:)); %renormalize filter (sum of all filter pixels


frameindA(frameindA>fileinfo. NumFrames)=[];

n_chunks=ceil(numel(frameindA)/ind_chunk);
eyefit=[];
    inip=[];
    inir=[];
    inim=[];
    inima=[];
eyefit_ind_shift=0;
for n_ch=0:n_chunks-1
    display(['Reading and processing chunk ' num2str(n_ch+1)]);
    
    tic
    frameind=frameindA(n_ch*ind_chunk+1:min((n_ch+1)*ind_chunk,numel(frameindA)));
    frameind(frameind>fileinfo. NumFrames)=[];
    % read avi file
    warning off all;
    movD=aviread(aviName,frameind);
    warning on all;
    sz=length(movD);
    dsz=size(movD(1).cdata);
    
    
    
    
    if nargin<2 || isempty(ROI)
        stack=zeros(dsz(1),dsz(2),sz);
        for fr=1:length(frameind)
            stack(:,:,fr)=double(movD(fr).cdata);
        end
    else
        stack=zeros(ROI(4)-ROI(2)+1,ROI(3)-ROI(1)+1,sz);
        for fr=1:length(frameind)
            stack(:,:,fr)=double(movD(fr).cdata(ROI(2):ROI(4),ROI(1):ROI(3)));
        end
    end
    

    
    for fr=1:length(frameind)
        stack(:,:,fr)=conv2(stack(:,:,fr),sp_filter,'same');
    end
    %filter end
    
    toc
    tic
    fit=FindEyeFromStack( stack,subsampl, size_range_1, size_range_2,dilation_range, motion_range, contrast_tsh,squeeze_coeff,inip,inir,inim,inima);
    
    inip=fit(end).pos;
    inir=fit(end).ra;
    inim=fit(end).minmax;
    inima=fit(end).m;

    
    eyefit=[eyefit fit];
    tottime=toc;
    
    display(['Time=' num2str(tottime) ' , avg time=' num2str(tottime/length(frameind)) ' sec/frame']);
    
    
    subsz=size(stack);
    %Save frames
    display(['Saving chunk ' num2str(n_ch+1)]);
    tic
    if saveFrames
        sl=strfind(aviName,'\');
        dirname=aviName(1:sl(end)-1);
        warning off all;
        mkdir(dirname,'fit');
        mkdir(dirname,'debug');
        warning on all;
        
        for fr=1:length(frameind)
            fname=[dirname '\fit\fit_' num2str(frameind(fr)) '.tiff'];
            
            if nargin<2 || isempty(ROI)
                subframe=movD(fr).cdata;
            else
                subframe=movD(fr).cdata(ROI(2):ROI(4),ROI(1):ROI(3));
            end
            posR=round(eyefit(fr+eyefit_ind_shift).pos);
            ra=eyefit(fr+eyefit_ind_shift).ra;
            fitgood=eyefit(fr+eyefit_ind_shift).fitgood;
            if fitgood(1)==1
                subframe(posR(2), max(posR(1)-ra(1),1):min(posR(1)+ra(1),subsz(2)))=0;
                subframe(max(posR(2)-round(ra(1)*squeeze_coeff),1):min(posR(2)+round(ra(1)*squeeze_coeff),subsz(1)), posR(1))=0;
            else
                subframe(max(posR(2)-2,1):min(posR(2)+2,subsz(1)), max(posR(1)-2,1):min(posR(1)+2,subsz(2)))=0;
            end
            
           if fitgood(2)==1
                subframe(posR(4), round(max(posR(3)-ra(2),1)):round(min(posR(3)+ra(2),subsz(2))))=255;
                subframe(max(posR(4)-round(ra(2)*squeeze_coeff),1):min(posR(4)+round(ra(2)*squeeze_coeff),subsz(1)), posR(3))=255;
            else
                subframe(max(posR(4)-2,1):min(posR(4)+2,subsz(1)), max(posR(3)-2,1):min(posR(3)+2,subsz(2)))=255;
            end
            imwrite(subframe,fname,'tif','Compression','none');
        end
%         for fr=1:length(frameind)
%             fname=[dirname '\debug\fit_' num2str(frameind(fr)) '.tiff'];
%             imwrite(eyefit(fr).debugImg,fname,'tif','Compression','none');
%         end
        
        
    end
    eyefit_ind_shift=eyefit_ind_shift+length(frameind);
    toc
    clear stack;
    clear movD;
end