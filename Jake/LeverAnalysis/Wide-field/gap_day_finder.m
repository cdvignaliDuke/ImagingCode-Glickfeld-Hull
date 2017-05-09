function gap_inx = gap_day_finder(days);
%find the gaps where training days were not consecutive. For a list of days this function convets them to numbers and returns the values between the
%indeces of non-consecutive days. e.g. 2.5 if days 2 and 3 were not back to back. 

%establish order of days in a year. 
order_of_days = cat(2, [170101:170131],[170201:170228],[170301:170331], [170401:170430], [170501:170531], [170601:170630], ...
    [170701:170731], [170801:170831], [170901:170930], [171001,171031], [171101:171130], [171201:171231]);
gap_inx = [];

for ii = 1:[length(days)-1]
    %convert days info to numbers. Get index of analyzed day
    this_day = str2double(days{ii}(1:6));
    next_day = str2double(days{ii+1}(1:6));
    this_day_inx = find(order_of_days==this_day);
    
    %ensure that the dates listed are included in the list of days
    assert(ismember(this_day, order_of_days));
    assert(ismember(next_day, order_of_days));
    
    %determine if current day and next day are consecutive
    if next_day == order_of_days(this_day_inx+1);
        continue
    else
        gap_inx = [gap_inx, ii+0.5];
    end
end