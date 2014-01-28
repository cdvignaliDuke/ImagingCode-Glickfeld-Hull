function [imgs2,freqs,oris] = hartly(width,s)
%HARTLY Hartly transform stimuli
%[IMGS2,FREQS,ORIS] = HARTLY(WIDTH,S)
% where S is SCREEN structure
% see SCREEN

%% hartly select based on sf
siz = max(s.Resolution);
Fs =mean(1./s.PixelSize);
[wx,wy]=freqspace(siz,'meshgrid');

%% calc sf and oris
w = wx+j*wy;
freqs = abs(w)*Fs;
oris = round(unwrap(angle(w))/pi*180);

os = round((siz-width)/2);
sel= [1:width]+os;
oris = oris(sel,sel);
freqs = freqs(sel,sel);

%% calc
ax=[];

imgs = {};
for iy = 1:width
    for ix = 1:width
        ind=(iy-1)*width+ix;
        im = zeros(siz,siz);
        im(iy+os,ix+os)=1;
        X = fft2(double(fftshift(im)));
        Y = im2double((real(X)-imag(X))/sqrt(2)+1)/2.0;
        imgs{ind}=Y(1:s.Resolution(1),1:s.Resolution(2));
    end
end

%% sort per sf and ori
w = w(sel,sel)';
[vals,ord]=sort(w(:));
imgs = imgs(ord);

%% add negative images
imgs2 ={};

imgs2{1}= im2uint8(1.0*ones(s.Resolution));
imgs2{2}= im2uint8(0.0*ones(s.Resolution));

ind2 = 3;
for ind = 2:length(imgs);
    imgs2{ind2}=im2uint8(imgs{ind});
    ind2 = ind2+1;
    if ind > 1
        imgs2{ind2}=im2uint8(-imgs{ind}+1);
        ind2 = ind2+1;
    end
end

imgs2{end+1}= im2uint8(0.5*ones(s.Resolution));

%% display
for ind = 1:length(imgs2)
    ax(ind) = subplot(width*2,width+1,ind);
    imshow(imgs2{ind});colormap gray;box off;axis off;
end

set(ax,'plotboxaspectratio',[fliplr(s.Resolution) 1]);

return;