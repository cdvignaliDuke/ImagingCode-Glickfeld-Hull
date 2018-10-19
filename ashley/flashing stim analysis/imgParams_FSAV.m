%% colors
cueColor = {[0 0 0];[.5 .5 1]};
AVColor = {[0 0 0];[0.5 0.5 0.5]};
hiLoColor = {[0.5 0.5 0.5];[0 0 0]};

%% example cells
if strcmp(ds,'FSAV_attentionV1')
    exampleCell_1 = 418; %738 % first-stim responsive
    exampleCell_2 = 1223;%1269;%1269 % late responsive
    exampleCell_3 = 543;%386; % late suppressed
    
    attnExCell_1 = 366; %first-stim responsive and modulated by attention(+V)
elseif strcmp(ds,'FSAV_V1_100ms_naive')    
    attnExCell_1 = 506;%603; %first-stim responsive 
end