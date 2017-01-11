function b_data = get_bx_data(bx_source, day);
%find and load behavior file
bfile = dir([bx_source 'data-*i9' day(end-1:end) '-' day(1:6) '*' ]);
behave_dest = [bx_source bfile.name];
assert(length(bfile)) = 1;
b_data = load(behave_dest);
b_data = b_data.input;
end