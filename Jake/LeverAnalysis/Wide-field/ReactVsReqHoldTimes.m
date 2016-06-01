clear
BEHAVE_DIR = 'Z:\Data\WidefieldImaging\GCaMP\behavior\';
days = {'150518_img24', '150519_img24', '150518_img25', '150517_img25', '150716_img27', '150718_img27', '150716_img28', '150717_img28', '151021_img29', '151022_img29', '151009_img30', '151011_img30', '151211_img32', '151212_img32', '160129_img35', '160131_img35', '160129_img36','160131_img36', '160314_img38', '160315_img38', '160319_img41', '160320_img41'}; %'150718_img27', '150719_img27',
%days = {'150716_img28'};
figure
hold on; 
for kk=1:length(days)
    day_name  =  days{kk};
    bfile = dir([BEHAVE_DIR 'data-*i9' days{kk}(end-1:end) '-' days{kk}(1:6) '*' ]);
    behave_dest = [BEHAVE_DIR bfile.name];
    b_data = load(behave_dest);
    if b_data.input.randReqHoldMaxMs < 4500;
        continue
    end
    days{kk}
    reactTimesMs = cell2mat(b_data.input.reactTimesMs);
    reqHoldTimeMs = round(cell2mat(b_data.input.tTotalReqHoldTimeMs));
    scatter(reqHoldTimeMs, reactTimesMs);
end
refline(0,0);

