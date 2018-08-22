function caldt = getDateFromNumberForm(dateNum)
% dateNum can be a scalar or a vector

n = length(dateNum);
caldt = [];
for i = 1:n
    dateChar = char(num2str(dateNum(i)));

    yr = num2str(str2num(dateChar(1:2))+2000);
    m = dateChar(3:4);
    d = dateChar(5:6);
    
    caldt = cat(2,caldt,...
        datetime(datestr(datenum([yr '-' m '-' d],'yyyy-mm-dd'),29)));
end


end