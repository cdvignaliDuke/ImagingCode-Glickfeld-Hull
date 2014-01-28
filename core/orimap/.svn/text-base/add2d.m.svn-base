function Aout = add2d(A, Bx)

if Bx == 1
    Aout = A;
else
    x=size(A,1);
    y=size(A,2);
    Aout=zeros(x/Bx,y);
    for i=1:x/Bx
        for j=1:Bx
            Aout(i,:) = Aout(i,:) + A((i-1)*Bx+j,:);
        end
    end
end
