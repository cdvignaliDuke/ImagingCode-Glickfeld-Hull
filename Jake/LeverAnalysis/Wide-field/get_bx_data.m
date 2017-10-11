function b_data = get_bx_data(bx_source, day);
%if statement to distinguish between 900s and 000s and find bx file
if day(end-2) == 'g'
    bfile = dir([bx_source 'data-i' '*' day(end-1:end) '-' day(1:6) '*' ]);
elseif day(end-2) =='0'
    bfile = dir([bx_source 'data-i' '*' day(end-2:end) '-' day(1:6) '*' ]);
end

%load behavior file
behave_dest = [bx_source bfile.name];
assert(length(bfile)) = 1;
b_data = load(behave_dest);
b_data = b_data.input;
end