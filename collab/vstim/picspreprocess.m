function picspreprocess(sourcedir,targetdir,s,options)
% picspreprocess(sourcedir,targetdir,s, options)
% options = [pad,nfiles,tile,inv] where 
%  pad = 0 means output file names are zero padded
%  nfiles number of files to process
%  tile input picture to to get more output pictures
%  inv = 1 also saves inverse version of pictures

if nargin < 4
    options = [0 inf 0 0];
end

if options(1) % no padding, to be read directly by presentation
%    fmt = 'pic%d.bmp';
    fmt = 'stim%i frame0.bmp';
    ipic = 0;
else % zero-padded
    fmt = 'pic%03d.bmp';
    ipic = 1;
end

siz = [240 320];

list = dir(fullfile(sourcedir,'*.tif'));

fs = 1/mean(s.PixelSize) % samples  / deg
wn = [.5]/(fs/2);
b = fir1(24,wn);
h = ftrans2(b);
    
if ~exist(targetdir)
    mkdir(targetdir);
else
    delete(fullfile(targetdir,'*.bmp'));
end

for index = 1:min(length(list),options(2));
    fn = list(index).name;   
    fprintf(1,'Reading %s\n',fn);
    rgb = imread(fullfile(sourcedir,fn)); % imagesc(rgb)
    [ny,nx,dummy]=size(rgb);
    
    if options(3)
        tilesiz = floor([ny nx]./siz);
    else
        tilesiz = [1 1];
    end

    img = rgb2gray(rgb); % imagesc(img) ; colormap gray;
	imgf = imfilter(im2double(img),h); % imagesc(img3) ; colormap gray;
    
    for irow = 1:tilesiz(1)
        for icol = 1:tilesiz(2)
            iy = (irow-1)*siz(1)+1:irow*siz(1);
            ix = (icol-1)*siz(2)+1:icol*siz(2);   
            
            img1 = img(iy,ix);

            % discard frames with saturated pixels
            if length(find(img1(:)>253))/prod(siz) > 0.01 % threshold 1%
                fprintf(1,'Dropped this frame because saturated\n');
                continue;
            end
                
            img2 = imadjust(imgf(iy,ix),stretchlim(imgf(iy,ix)),[]);
            img3 = histeq(img2); % equalize stimulus mean and contrast
        	img3 = imfilter(img3,h); % imagesc(img3) ; colormap gray;
            
            % write final picture
            img4 = repmat(im2uint8(img3),[1 1 3]);   
            fn = sprintf(fmt,ipic);   
            imwrite(img4,fullfile(targetdir,fn),'bmp'); 
            fprintf(1,'Writing %s\n',fn);
            ipic = ipic + 1;
            imshow(img4);drawnow;
            
            % write inverse picture
            if options(4)
                img4 = repmat(im2uint8(-img3+1),[1 1 3]);   
                fn = sprintf(fmt,ipic);   
                imwrite(img4,fullfile(targetdir,fn),'bmp'); 
                fprintf(1,'Writing %s\n',fn);
                ipic = ipic + 1;
                imshow(img4);drawnow;
            end
        end
    end
end

% bright flash
img3 = 1*ones(siz);
img4 = repmat(im2uint8(img3),[1 1 3]);            
fn = sprintf(fmt,ipic);   
fprintf(1,'Writing %s\n',fn);
imwrite(img4,fullfile(targetdir,fn),'bmp'); 
ipic = ipic + 1;

% dark flash
img3 = 0*ones(siz);
img4 = repmat(im2uint8(img3),[1 1 3]);            
fn = sprintf(fmt,ipic);   
fprintf(1,'Writing %s\n',fn);
imwrite(img4,fullfile(targetdir,fn),'bmp'); 
ipic = ipic + 1;

% last picture is blank
img3 = 0.5*ones(siz);
img4 = repmat(im2uint8(img3),[1 1 3]);            
fn = sprintf(fmt,ipic);   
fprintf(1,'Writing %s\n',fn);
imwrite(img4,fullfile(targetdir,fn),'bmp'); 
 
return