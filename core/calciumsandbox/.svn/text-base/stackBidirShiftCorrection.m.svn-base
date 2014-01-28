function [outstack,shift]=stackBidirShiftCorrection(stack);
sz=size(stack);
shift=zeros(sz(3),1);
outstack=zeros(sz(1),sz(2),sz(3));
for i=1:sz(3)
    shift(i)=-CalcDelay(stack(:,:,i))*2;
    outstack(:,:,i)=imBidirShift(stack(:,:,i),round(shift(i)));
    if mod(i,100)==0
        fprintf('%d ',i);
    end
end
