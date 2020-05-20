%% MASSIVE SUBPLOT with corr
figure;
start = 1;
for iCell = 1:15
    width = 24; height = 24;
    xCenter = round(cell_stats(iCell).Centroid(2));
    yCenter = round(cell_stats(iCell).Centroid(1));
    xLeft = (xCenter - width/2);
    yBottom = (yCenter - height/2);
    if xLeft > 12 && xLeft < 488 && yBottom > 12 && yBottom < 772
    subplot(15,27,start)
    a = pix_3hz1(xLeft:(xLeft+width),yBottom:(height+yBottom));
    imagesc(a)
    pos = get(gca, 'Position');
    pos(1) = 0.025;
    pos(3) = 0.02;
    set(gca, 'Position', pos)
    axis off
    axis square
    end
    xCenter2 = round(cell_stats2(iCell).Centroid(2));
    yCenter2 = round(cell_stats2(iCell).Centroid(1));
    xLeft2 = (xCenter2 - width/2);
    yBottom2 = (yCenter2 - height/2);
    if xLeft2 > 12 && xLeft2 < 488 && yBottom2 > 12 && yBottom2 < 772
    subplot(15,27,start+1)
    x = pix_fgta2(xLeft2:(xLeft2+width),yBottom2:(height+yBottom2));
    [reg shift] = shift_opt(x,a,5);
    imagesc(reg)
    pos = get(gca, 'Position');
    pos(1) = 0.05;
    pos(3) = 0.02;
    set(gca, 'Position', pos)
    axis off
    axis square
    r = triu2vec(corrcoef(a(:),reg(:)));
    til_str = num2str(chop(r,2));
    title(til_str,'FontSize',6);
    end
%     xCenter3 = round(cell_stats3(iCell).Centroid(2));
%     yCenter3 = round(cell_stats3(iCell).Centroid(1));
%     xLeft3 = (xCenter3 - width/2);
%     yBottom3 = (yCenter3 - height/2);
    if xLeft2 > 12 && xLeft2 < 488 && yBottom2 > 12 && yBottom2 < 772
    subplot(15,27,start+2)
y = pix_fgta2(xLeft2:(xLeft2+width),yBottom2:(height+yBottom2));
    imagesc(y)
    pos = get(gca, 'Position');
    pos(1) = 0.075;
    pos(3) = 0.02;
    set(gca, 'Position', pos)
    axis square
    axis off
    r = triu2vec(corrcoef(a(:),y(:)));
    til_str = num2str(chop(r,2));
    title(til_str,'FontSize',6);
    end
    
    xCenter = round(cell_stats(iCell+15).Centroid(2));
    yCenter = round(cell_stats(iCell+15).Centroid(1));
    xLeft = (xCenter - width/2);
    yBottom = (yCenter - height/2);
    if xLeft > 12 && xLeft < 488 && yBottom > 12 && yBottom < 772
    subplot(15,27,start+3)
    a = pix_3hz1(xLeft:(xLeft+width),yBottom:(height+yBottom)); 
    imagesc(a)
    pos = get(gca, 'Position');
    pos(1) = 0.125;
    pos(3) = 0.02;
    set(gca, 'Position', pos)
    axis off
    axis square
    end
    xCenter2 = round(cell_stats2(iCell+15).Centroid(2));
    yCenter2 = round(cell_stats2(iCell+15).Centroid(1));
    xLeft2 = (xCenter2 - width/2);
    yBottom2 = (yCenter2 - height/2);
    if xLeft2 > 12 && xLeft2 < 488 && yBottom2 > 12 && yBottom2 < 772
    subplot(15,27,start+4)
    x = pix_fgta2(xLeft2:(xLeft2+width),yBottom2:(height+yBottom2));
    [reg shift] = shift_opt(x,a,5);
    imagesc(reg)
    pos = get(gca, 'Position');
    pos(1) = 0.15;
    pos(3) = 0.02;
    set(gca, 'Position', pos)
    axis off
    axis square
    r = triu2vec(corrcoef(a(:),reg(:)));  
    til_str = num2str(chop(r,2));
    title(til_str,'FontSize',6);
    end
%     xCenter3 = round(cell_stats3(iCell+15).Centroid(2));
%     yCenter3 = round(cell_stats3(iCell+15).Centroid(1));
%     xLeft3 = (xCenter3 - width/2);
%     yBottom3 = (yCenter3 - height/2);
    if xLeft2 > 12 && xLeft2 < 488 && yBottom2 > 12 && yBottom2 < 772
    subplot(15,27,start+5)
    y = pix_fgta2(xLeft2:(xLeft2+width),yBottom2:(height+yBottom2));
    imagesc(y)
    pos = get(gca, 'Position');
    pos(1) = 0.175;
    pos(3) = 0.02;
    set(gca, 'Position', pos)
    axis square
    axis off
    r = triu2vec(corrcoef(a(:),y(:)));
    til_str = num2str(chop(r,2));
    title(til_str,'FontSize',6);
    end
    
    xCenter = round(cell_stats(iCell+30).Centroid(2));
    yCenter = round(cell_stats(iCell+30).Centroid(1));
    xLeft = (xCenter - width/2);
    yBottom = (yCenter - height/2);
    if xLeft > 12 && xLeft < 488 && yBottom > 12 && yBottom < 772
    subplot(15,27,start+6)
    a = pix_3hz1(xLeft:(xLeft+width),yBottom:(height+yBottom)); 
    imagesc(a)
    pos = get(gca, 'Position');
    pos(1) = 0.225;
    pos(3) = 0.02;
    set(gca, 'Position', pos)
    axis off
    axis square
    end
    xCenter2 = round(cell_stats2(iCell+30).Centroid(2));
    yCenter2 = round(cell_stats2(iCell+30).Centroid(1));
    xLeft2 = (xCenter2 - width/2);
    yBottom2 = (yCenter2 - height/2);
    if xLeft2 > 12 && xLeft2 < 488 && yBottom2 > 12 && yBottom2 < 772
    subplot(15,27,start+7)
    x = pix_fgta2(xLeft2:(xLeft2+width),yBottom2:(height+yBottom2));
    [reg shift] = shift_opt(x,a,5);
    imagesc(reg)
    pos = get(gca, 'Position');
    pos(1) = 0.25;
    pos(3) = 0.02;
    set(gca, 'Position', pos)
    axis off
    axis square
    r = triu2vec(corrcoef(a(:),reg(:)));  
    til_str = num2str(chop(r,2));
    title(til_str,'FontSize',6);
    end
%     xCenter3 = round(cell_stats3(iCell+30).Centroid(2));
%     yCenter3 = round(cell_stats3(iCell+30).Centroid(1));
%     xLeft3 = (xCenter3 - width/2);
%     yBottom3 = (yCenter3 - height/2);
    if xLeft2 > 12 && xLeft2 < 488 && yBottom2 > 12 && yBottom2 < 772
    subplot(15,27,start+8)
    y = pix_fgta2(xLeft2:(xLeft2+width),yBottom2:(height+yBottom2));
    imagesc(y)
    pos = get(gca, 'Position');
    pos(1) = 0.275;
    pos(3) = 0.02;
    set(gca, 'Position', pos)
    axis square
    axis off
    r = triu2vec(corrcoef(a(:),y(:)));
    til_str = num2str(chop(r,2));
    title(til_str,'FontSize',6);
    end
    
    xCenter = round(cell_stats(iCell+45).Centroid(2));
    yCenter = round(cell_stats(iCell+45).Centroid(1));
    xLeft = (xCenter - width/2);
    yBottom = (yCenter - height/2);
    if xLeft > 12 && xLeft < 488 && yBottom > 12 && yBottom < 772
    subplot(15,27,start+9)
    a = pix_3hz1(xLeft:(xLeft+width),yBottom:(height+yBottom)); 
    imagesc(a)
    pos = get(gca, 'Position');
    pos(1) = 0.325;
    pos(3) = 0.02;
    set(gca, 'Position', pos)
    axis off
    axis square
    end
    xCenter2 = round(cell_stats2(iCell+45).Centroid(2));
    yCenter2 = round(cell_stats2(iCell+45).Centroid(1));
    xLeft2 = (xCenter2 - width/2);
    yBottom2 = (yCenter2 - height/2);
    if xLeft2 > 12 && xLeft2 < 488 && yBottom2 > 12 && yBottom2 < 772
    subplot(15,27,start+10)
    x = pix_fgta2(xLeft2:(xLeft2+width),yBottom2:(height+yBottom2));
    [reg shift] = shift_opt(x,a,5);
    imagesc(reg)
    pos = get(gca, 'Position');
    pos(1) = 0.35;
    pos(3) = 0.02;
    set(gca, 'Position', pos)
    axis off
    axis square
    r = triu2vec(corrcoef(a(:),reg(:)));  
    til_str = num2str(chop(r,2));
    title(til_str,'FontSize',6);
    end
%     xCenter3 = round(cell_stats3(iCell+45).Centroid(2));
%     yCenter3 = round(cell_stats3(iCell+45).Centroid(1));
%     xLeft3 = (xCenter3 - width/2);
%     yBottom3 = (yCenter3 - height/2);
    if xLeft2 > 12 && xLeft2 < 488 && yBottom2 > 12 && yBottom2 < 772
    subplot(15,27,start+11)
    y = pix_fgta2(xLeft2:(xLeft2+width),yBottom2:(height+yBottom2));
    imagesc(y)
    pos = get(gca, 'Position');
    pos(1) = 0.375;
    pos(3) = 0.02;
    set(gca, 'Position', pos)
    axis square
    axis off
    r = triu2vec(corrcoef(a(:),y(:)));
    til_str = num2str(chop(r,2));
    title(til_str,'FontSize',6);
    end
    
    xCenter = round(cell_stats(iCell+60).Centroid(2));
    yCenter = round(cell_stats(iCell+60).Centroid(1));
    xLeft = (xCenter - width/2);
    yBottom = (yCenter - height/2);
    if xLeft > 12 && xLeft < 488 && yBottom > 12 && yBottom < 772
    subplot(15,27,start+12)
    a = pix_3hz1(xLeft:(xLeft+width),yBottom:(height+yBottom)); 
    imagesc(a)
    pos = get(gca, 'Position');
    pos(1) = 0.425;
    pos(3) = 0.02;
    set(gca, 'Position', pos)
    axis off
    axis square
    end
    xCenter2 = round(cell_stats2(iCell+60).Centroid(2));
    yCenter2 = round(cell_stats2(iCell+60).Centroid(1));
    xLeft2 = (xCenter2 - width/2);
    yBottom2 = (yCenter2 - height/2);
    if xLeft2 > 12 && xLeft2 < 488 && yBottom2 > 12 && yBottom2 < 772
    subplot(15,27,start+13)
    x = pix_fgta2(xLeft2:(xLeft2+width),yBottom2:(height+yBottom2));
    [reg shift] = shift_opt(x,a,5);
    imagesc(reg)
    pos = get(gca, 'Position');
    pos(1) = 0.45;
    pos(3) = 0.02;
    set(gca, 'Position', pos)
    axis off
    axis square
    r = triu2vec(corrcoef(a(:),reg(:)));  
    til_str = num2str(chop(r,2));
    title(til_str,'FontSize',6);
    end
%     xCenter3 = round(cell_stats3(iCell+60).Centroid(2));
%     yCenter3 = round(cell_stats3(iCell+60).Centroid(1));
%     xLeft3 = (xCenter3 - width/2);
%     yBottom3 = (yCenter3 - height/2);
    if xLeft2 > 12 && xLeft2 < 488 && yBottom2 > 12 && yBottom2 < 772
    subplot(15,27,start+14)
    y = pix_fgta2(xLeft2:(xLeft2+width),yBottom2:(height+yBottom2));
    imagesc(y)
    pos = get(gca, 'Position');
    pos(1) = 0.475;
    pos(3) = 0.02;
    set(gca, 'Position', pos)
    axis square
    axis off
    r = triu2vec(corrcoef(a(:),y(:)));
    til_str = num2str(chop(r,2));
    title(til_str,'FontSize',6);
    end
       
    xCenter = round(cell_stats(iCell+75).Centroid(2));
    yCenter = round(cell_stats(iCell+75).Centroid(1));
    xLeft = (xCenter - width/2);
    yBottom = (yCenter - height/2);
    if xLeft > 12 && xLeft < 488 && yBottom > 12 && yBottom < 772
    subplot(15,27,start+15)
    a = pix_3hz1(xLeft:(xLeft+width),yBottom:(height+yBottom)); 
    imagesc(a)
    pos = get(gca, 'Position');
    pos(1) = 0.525;
    pos(3) = 0.02;
    set(gca, 'Position', pos)
    axis off
    axis square
    end
    xCenter2 = round(cell_stats2(iCell+75).Centroid(2));
    yCenter2 = round(cell_stats2(iCell+75).Centroid(1));
    xLeft2 = (xCenter2 - width/2);
    yBottom2 = (yCenter2 - height/2);
    if xLeft2 > 12 && xLeft2 < 488 && yBottom2 > 12 && yBottom2 < 772
    subplot(15,27,start+16)
    x = pix_fgta2(xLeft2:(xLeft2+width),yBottom2:(height+yBottom2));
    [reg shift] = shift_opt(x,a,5);
    imagesc(reg)
    pos = get(gca, 'Position');
    pos(1) = 0.55;
    pos(3) = 0.02;
    set(gca, 'Position', pos)
    axis off
    axis square
    r = triu2vec(corrcoef(a(:),reg(:)));  
    til_str = num2str(chop(r,2));
    title(til_str,'FontSize',6);
    end
%     xCenter3 = round(cell_stats3(iCell+75).Centroid(2));
%     yCenter3 = round(cell_stats3(iCell+75).Centroid(1));
%     xLeft3 = (xCenter3 - width/2);
%     yBottom3 = (yCenter3 - height/2);
    if xLeft2 > 12 && xLeft2 < 488 && yBottom2 > 12 && yBottom2 < 772
    subplot(15,27,start+17)
    y = pix_fgta2(xLeft2:(xLeft2+width),yBottom2:(height+yBottom2));
    imagesc(y)
    pos = get(gca, 'Position');
    pos(1) = 0.575;
    pos(3) = 0.02;
    set(gca, 'Position', pos)
    axis square
    axis off
    r = triu2vec(corrcoef(a(:),y(:)));
    til_str = num2str(chop(r,2));
    title(til_str,'FontSize',6);
    end
       
    xCenter = round(cell_stats(iCell+90).Centroid(2));
    yCenter = round(cell_stats(iCell+90).Centroid(1));
    xLeft = (xCenter - width/2);
    yBottom = (yCenter - height/2);
    if xLeft > 12 && xLeft < 488 && yBottom > 12 && yBottom < 772
    subplot(15,27,start+18)
    a = pix_3hz1(xLeft:(xLeft+width),yBottom:(height+yBottom)); 
    imagesc(a)
    pos = get(gca, 'Position');
    pos(1) = 0.625;
    pos(3) = 0.02;
    set(gca, 'Position', pos)
    axis off
    axis square
    end
    xCenter2 = round(cell_stats2(iCell+90).Centroid(2));
    yCenter2 = round(cell_stats2(iCell+90).Centroid(1));
    xLeft2 = (xCenter2 - width/2);
    yBottom2 = (yCenter2 - height/2);
    if xLeft2 > 12 && xLeft2 < 488 && yBottom2 > 12 && yBottom2 < 772
     subplot(15,27,start+19)
    x = pix_fgta2(xLeft2:(xLeft2+width),yBottom2:(height+yBottom2));
    [reg shift] = shift_opt(x,a,5);
    imagesc(reg)
    pos = get(gca, 'Position');
    pos(1) = 0.65;
    pos(3) = 0.02;
    set(gca, 'Position', pos)
    axis off
    axis square
    r = triu2vec(corrcoef(a(:),reg(:)));  
    til_str = num2str(chop(r,2));
    title(til_str,'FontSize',6);
    end
%     xCenter3 = round(cell_stats3(iCell+90).Centroid(2));
%     yCenter3 = round(cell_stats3(iCell+90).Centroid(1));
%     xLeft3 = (xCenter3 - width/2);
%     yBottom3 = (yCenter3 - height/2);
    if xLeft2 > 12 && xLeft2 < 488 && yBottom2 > 12 && yBottom2 < 772
    subplot(15,27,start+20)
    y = pix_fgta2(xLeft2:(xLeft2+width),yBottom2:(height+yBottom2));
    imagesc(y)
    pos = get(gca, 'Position');
    pos(1) = 0.675;
    pos(3) = 0.02;
    set(gca, 'Position', pos)
    axis square
    axis off
    r = triu2vec(corrcoef(a(:),y(:)));
    til_str = num2str(chop(r,2));
    title(til_str,'FontSize',6);
    end
    
    xCenter = round(cell_stats(iCell+105).Centroid(2));
    yCenter = round(cell_stats(iCell+105).Centroid(1));
    xLeft = (xCenter - width/2);
    yBottom = (yCenter - height/2);
    if xLeft > 12 && xLeft < 488 && yBottom > 12 && yBottom < 772
     subplot(15,27,start+21)
    a = pix_3hz1(xLeft:(xLeft+width),yBottom:(height+yBottom)); 
    imagesc(a)
    pos = get(gca, 'Position');
    pos(1) = 0.725;
    pos(3) = 0.02;
    set(gca, 'Position', pos)
    axis off
    axis square
    end
    xCenter2 = round(cell_stats2(iCell+105).Centroid(2));
    yCenter2 = round(cell_stats2(iCell+105).Centroid(1));
    xLeft2 = (xCenter2 - width/2);
    yBottom2 = (yCenter2 - height/2);
    if xLeft2 > 12 && xLeft2 < 488 && yBottom2 > 12 && yBottom2 < 772
    subplot(15,27,start+22)
    x = pix_fgta2(xLeft2:(xLeft2+width),yBottom2:(height+yBottom2));
    [reg shift] = shift_opt(x,a,5);
    imagesc(reg)
    pos = get(gca, 'Position');
    pos(1) = 0.75;
    pos(3) = 0.02;
    set(gca, 'Position', pos)
    axis off
    axis square
    r = triu2vec(corrcoef(a(:),reg(:)));  
    til_str = num2str(chop(r,2));
    title(til_str,'FontSize',6);
    end
%     xCenter3 = round(cell_stats3(iCell+105).Centroid(2));
%     yCenter3 = round(cell_stats3(iCell+105).Centroid(1));
%     xLeft3 = (xCenter3 - width/2);
%     yBottom3 = (yCenter3 - height/2);
    if xLeft2 > 12 && xLeft2 < 488 && yBottom2 > 12 && yBottom2 < 772
    subplot(15,27,start+23)
    y = pix_fgta2(xLeft2:(xLeft2+width),yBottom2:(height+yBottom2));
    imagesc(y)
    pos = get(gca, 'Position');
    pos(1) = 0.775;
    pos(3) = 0.02;
    set(gca, 'Position', pos)
    axis square
    axis off
    r = triu2vec(corrcoef(a(:),y(:)));
    til_str = num2str(chop(r,2));
    title(til_str,'FontSize',6);
    end
    
    xCenter = round(cell_stats(iCell+120).Centroid(2));
    yCenter = round(cell_stats(iCell+120).Centroid(1));
    xLeft = (xCenter - width/2);
    yBottom = (yCenter - height/2);
    if xLeft > 12 && xLeft < 488 && yBottom > 12 && yBottom < 772
    subplot(15,27,start+24)   
    a = pix_3hz1(xLeft:(xLeft+width),yBottom:(height+yBottom)); 
    imagesc(a)
    pos = get(gca, 'Position');
    pos(1) = 0.825;
    pos(3) = 0.02;
    set(gca, 'Position', pos)
    axis off
    axis square
    end
    xCenter2 = round(cell_stats2(iCell+120).Centroid(2));
    yCenter2 = round(cell_stats2(iCell+120).Centroid(1));
    xLeft2 = (xCenter2 - width/2);
    yBottom2 = (yCenter2 - height/2);
    if xLeft2 > 12 && xLeft2 < 488 && yBottom2 > 12 && yBottom2 < 772
    subplot(15,27,start+25)
    x = pix_fgta2(xLeft2:(xLeft2+width),yBottom2:(height+yBottom2));
    [reg shift] = shift_opt(x,a,5);
    imagesc(reg)
    pos = get(gca, 'Position');
    pos(1) = 0.85;
    pos(3) = 0.02;
    set(gca, 'Position', pos)
    axis off
    axis square
    r = triu2vec(corrcoef(a(:),reg(:)));  
    til_str = num2str(chop(r,2));
    title(til_str,'FontSize',6);
    end
    if xLeft2 > 12 && xLeft2 < 488 && yBottom2 > 12 && yBottom2 < 772
%     xCenter3 = round(cell_stats3(iCell+120).Centroid(2));
%     yCenter3 = round(cell_stats3(iCell+120).Centroid(1));
%     xLeft3 = (xCenter3 - width/2);
%     yBottom3 = (yCenter3 - height/2);
    subplot(15,27,start+26)
    y = pix_fgta2(xLeft2:(xLeft2+width),yBottom2:(height+yBottom2));
    imagesc(y)
    pos = get(gca, 'Position');
    pos(1) = 0.875;
    pos(3) = 0.02;
    set(gca, 'Position', pos)
    axis square
    axis off
    r = triu2vec(corrcoef(a(:),y(:)));
    til_str = num2str(chop(r,2));
    title(til_str,'FontSize',6);
    end
    
    start = start+27;
end
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], ['AcrossAllDays'], ['CellMaps'], [ref_date '_' mouse '_' ref_str '_map_stackReg.pdf']),'-dpdf', '-bestfit')
