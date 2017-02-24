function isDivisResult = divisible(a,b)
%returns 1 if a is evenly divisible by b, 0 if not
if rem(a,b)==0 
isDivisResult = 1;
else
isDivisResult = 0;
end
end