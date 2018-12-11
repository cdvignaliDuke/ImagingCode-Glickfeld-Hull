%test to fix model selection problem

% bad examples: 50_4, 54_4, 4_1, 20_4(last one not thresholding)
% good examples: 4_4(m1), 5_4(m2)
combos = [4 1;
    4 4;
    5 4;
    50 4;
    54 4];
combos = [1:nCells; 4*ones(1,63)]';
combos = [20 4];

for i=1:size(combos,1)
    iCell = combos(i,1);
    iCon = combos(i,2);
    
    szs0 = sizeFits(iCell,iCon).szs0;
    data = sizeFits(iCell,iCon).data;
    c1 = sizeFits(iCell,iCon).fit1.c1;
    c2 = sizeFits(iCell,iCon).fit2.c2;
    cex = c2(1:3) + [0 c2(5) 0];
    cin = c2(4:6) + [0 0 c2(3)];
    
    figure(1);clf;
    subplot(1,3,1)
    plot(szs0,data,'.b')
    hold on
    plot(szRng, m1(c1,szRng),'-k')
    hold off
    title(['Cell ' num2str(iCell) ' con ' num2str(iCon) ' m1' sizeFits(iCell,iCon).Fstr1])
    subplot(1,3,2)
    plot(szs0,data,'.b')
    hold on
    plot(szRng, m2(c2,szRng),'-k')
    hold off
    title(['Cell ' num2str(iCell) ' con ' num2str(iCon) ' m2' sizeFits(iCell,iCon).Fstr2])
    subplot(1,3,3)
    plot(szs0,data,'.b')
    hold on
    plot(szRng, m1(cex,szRng),'--g')
    plot(szRng, m1(cin,szRng),'--r')
    text(7.5,sizeFits(iCell,iCon).maxResp2/2,['Ai:' num2str(cin(1)) ', ki:' num2str(cin(2))])
    text(7.5,sizeFits(iCell,iCon).maxResp2,['Max:' num2str(sizeFits(iCell,iCon).maxResp2) ', %in:' num2str(100*cin(1)/sizeFits(iCell,iCon).maxResp2)])
    hold off
    title('m2 components')
    pause
end


%% examine ki only
for i=1:numel(sizeFits)
    c2 = sizeFits(i).fit2.c2;
    ki(i) = c2(4);
end
kisub = ki([sizeFits(:).Ftest]); % only Ftest=1

figure(2);clf;
subplot(2,2,1)
histogram(ki,100)
title('ki')
subplot(2,2,2)
histogram(log(ki),100)
title('log(ki)')
subplot(2,2,3)
histogram(kisub,100)
title('ki, m2 only')
subplot(2,2,4)
histogram(log(kisub),100)
title('log(ki), m2 only')

%% explore ki in the model
figure(3);clf;
hold on
for i = 1:5
    plot(szRng,m1([1, 10^(-i), 30],szRng))
end
hold off
legend('1','2','3','4','5')