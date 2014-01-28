function eyefit=FindEyeFromStack( stack,subsampl, size_range_1, size_range_2,dilation_range, motion_range, contrast_tsh,squeeze_coeff,inip,inir,inim,inima)
% runs FastRadialTransformOneFrame on stack

if nargin<2
    subsampl=1;
end
if nargin<3 || isempty(size_range_1)
    size_range_1=5:20;
end

if nargin<4 || isempty(size_range_2)
    size_range_2=10:50;
end

if nargin<5 || isempty(dilation_range)
    dilation_range=1;
end

if nargin<6 || isempty(motion_range)
    motion_range=10;
end

if nargin<7 || isempty(contrast_tsh)
    contrast_tsh=0;
end

if nargin<8
    squeeze_coeff=1;
end


sz=size(stack);
%tic;

if nargin<9 || isempty(inim) || isempty(inip) || isempty(inir) || isempty(inima)
    mamax=-inf;
    mimin=inf;
    [pos ra ma mi debugImg]=FastRadialTransformOneFrame(stack(:,:,1), size_range_1, size_range_2, [],[],contrast_tsh,[],[],subsampl,8,squeeze_coeff);
    startfr=2;
    
%     ma=ma*subsampl.^2/length(size_range_1);
%     mi=mi*subsampl.^2/length(size_range_2);
    
    ma=ma/length(size_range_1);
    mi=mi/length(size_range_2);

    
    mamax=max(mamax,ma);
    mimin=min(mimin,mi);
    
    
    eyefit(1).pos=pos;
    eyefit(1).ra=ra;
    eyefit(1).m=[ma mi];
    eyefit(1).minmax=[mamax mimin];
    eyefit(1).fitgood=[1 1];
    eyefit(1).debugImg=debugImg;
    
    
else
    mamax=inim(1);
    mimin=inim(2);
    ma=inima(1);
    mi=inima(2);
    pos=inip;
    ra=inir;
    startfr=1;
end



%toc

% do for the rest of frames
for fr=startfr:sz(3)
%    fr
%    tic;
    
    % do fastradialtransform only for small region around previous position and
    % radius
    ma0=ma;
    mi0=mi;
    ra0=ra;
    pos0=pos;
    mamax=mamax*0.999;
    mimin=mimin*0.999;
    ra(1)=max(ra(1),size_range_1(1));
    ra(1)=min(ra(1),size_range_1(end));
    ran1=[ra(1)-1:ra(1)+1];
    ran2=[ra(2)-dilation_range:ra(2)+dilation_range];
    [pos ra ma mi debugImg]=FastRadialTransformOneFrame(stack(:,:,fr), ran1,ran2, pos,[motion_range motion_range],contrast_tsh,[],[],subsampl,8,squeeze_coeff);
    
%     ma=ma*subsampl.^2/length(ran1);
%     mi=mi*subsampl.^2/length(ran2);
    ma=ma/length(ran1);
    mi=mi/length(ran2);

    mamax=max(mamax,ma);
    mimin=min(mimin,mi);
%     [length(ran1) length(ran2)]
%     [ma ma0 mi mi0 mamax mimin]
    % pos
    % ra
    fitgood(1:2)=1;
   if mi>mimin*0.0001 && mod(fr,2)==0 %very bad, assuming it is blinking (every second frame) 
%    if ma<mamax*0.015 && mi>mimin*0.015 %very bad, assuming it is blinking
        fitgood(1:2)=0;
        ra=ra0;
        pos=pos0;
    else
        if ma<ma0*0.3 || mi>mi0*0.4 || ma<mamax*0.15 || mi>mimin*0.15 || isnan(ma) || isnan(mi)% if fit not good then repeat with larger range of corresponding parameters
%        if ma<ma0*0.3 || mi>mi0*0.4 || ma<mamax*0.10 || mi>mimin*0.10 || isnan(ma) || isnan(mi)
            %        [ma ma0 mi mi0]
            if ma<ma0*0.3 || ma<mamax*0.15 || isnan(ma)
                ran1=max(ra0(1)-5,size_range_1(1)):min(ra0(1)+5,size_range_1(end));
                motion_range_1=round(ra0(1)*1.0);
            else
                %ran1=ra0(1);
%                display('Good reflection, 1');
                motion_range_1=motion_range;
            end
            
            if mi>mi0*0.4 || mi>mimin*0.15 || isnan(mi)
                ran2=ra0(2)-dilation_range*5:ra0(2)+dilation_range*5;
                motion_range_2=round(ra0(2)*1.0);
            else
                %ran2=ra0(2);
%                display('Good pupl, 1');
                motion_range_2=motion_range;
            end

%            [pos ra ma mi debugImg]=FastRadialTransformOneFrame(stack(:,:,fr), ran1,ran2, pos0,[motion_range_1 motion_range_2],contrast_tsh,[],[],1,8,squeeze_coeff);
            [pos ra ma mi debugImg]=FastRadialTransformOneFrame(stack(:,:,fr), ran1,ran2, pos0,[motion_range_1 motion_range_2],contrast_tsh,[],[],subsampl,8,squeeze_coeff);
            
            %         ma=ma*subsampl.^2/length(ran1);
            %         mi=mi*subsampl.^2/length(ran2);
            
            ma=ma/length(ran1);
            mi=mi/length(ran2);
            
            mamax=max(mamax,ma);
            mimin=min(mimin,mi);
            %         pos
            %         ra
%             [length(ran1) length(ran2)]
%             [ma ma0 mi mi0 mamax mimin]
        end
        
        if ma<ma0*0.3 || mi>mi0*0.4  || ma<mamax*0.15 || mi>mimin*0.15 || isnan(ma) || isnan(mi)% if fit not good then repeat with even larger range of corresponding parameters
        % version - do it only for pupl
%        if mi>mi0*0.4  || mi>mimin*0.15 || isnan(mi)% if fit not good then repeat with even larger range of corresponding parameters
        
            %        [ma ma0 mi mi0]
            if ma<ma0*0.3 || ma<mamax*0.15 || isnan(ma)
                ran1=size_range_1;
                motion_range_1=inf;
            else
%                display('Good reflection, 2');
%               ran1=ra0(1);
           end
            
            if mi>mi0*0.4 || mi>mimin*0.15 || isnan(mi)
                ran2=size_range_2;
                motion_range_2=inf;
            else
%                display('Good pupl, 2');
%                ran2=ra0(2);
            end
%            [pos ra ma mi debugImg]=FastRadialTransformOneFrame(stack(:,:,fr), ran1,ran2, [],[],contrast_tsh,[],[],subsampl,8,squeeze_coeff);
%            [pos ra ma mi debugImg]=FastRadialTransformOneFrame(stack(:,:,fr), ran1,ran2, pos0,[motion_range_1 inf],contrast_tsh,[],[],subsampl,16,squeeze_coeff);
            [pos ra ma mi debugImg]=FastRadialTransformOneFrame(stack(:,:,fr), ran1,ran2, pos0,[motion_range_1 motion_range_2],contrast_tsh,[],[],subsampl,8,squeeze_coeff);
%             ma=ma*subsampl.^2/length(ran1);
%             mi=mi*subsampl.^2/length(ran2);
                    ma=ma/length(ran1);
                    mi=mi/length(ran2);
            
            mamax=max(mamax,ma);
            mimin=min(mimin,mi);
            %         pos
            %         ra
%             [length(ran1) length(ran2)]
%             [ma ma0 mi mi0 mamax mimin]
            if ma<ma0*0.3  || ma<mamax*0.15
                fitgood(1)=0;
            else
                fitgood(1)=1;
            end
            
            if mi>mi0*0.4  || mi>mimin*0.15
                fitgood(2)=0;
            else
                fitgood(2)=1;
            end
        end
    end
    
    eyefit(fr).pos=pos;
    eyefit(fr).ra=ra;
    eyefit(fr).m=[ma mi];
    eyefit(fr).minmax=[mamax mimin];
    eyefit(fr).fitgood=fitgood;
    eyefit(fr).debugImg=debugImg;
%    toc
end


