function [order,groups] = corrsort(r)
%CORRSORT
% [ORDER,GROUPS] = CORRSORT(R)

siz = size(r);
n = siz(1);

dist = -r;
for i = 1:n
dist(i,i) = inf;
end
for k = 1:n-1
[min1,ind1] = min(dist);
[min2,ind2] = min(min1);
r = ind1(ind2);
c = ind2;
% Segment to order samples here
if k == 1
 groups = zeros(round(n/2),n);
 groups(1,1:2) = [c r]; gi = 1;
else
 % does r belong to an existing group?
 [zr1,zr2] = find(groups==r);
 % does c belong to an existing group?
 [zc1,zc2] = find(groups==c);
 % If neither c nor r belong to a group they form their own
 if isempty(zr1)   %r doesn't belong to a group
   if isempty(zc1) %c doesn't belong to a group
     gi = gi+1;
     groups(gi,1:2) = [c r];
   else   % r doesn't belong but c does, add r to group c
     sgc = size(find(groups(zc1(1),:)));   %how big is group c
     % Figure out what side to add to
     cgc = groups(zc1(1),1:sgc(2));
     [mindg,inddg] = min([dist(cgc(1),r) dist(cgc(sgc(2)),r)]);
     if inddg == 2
       groups(zc1(1),sgc(2)+1) = r;
     else
       groups(zc1(1),1:sgc(2)+1) = [r groups(zc1(1),1:sgc(2))];
     end
   end
 else   %r does belong to a group
   if isempty(zc1) %c doesn't belong to a group, add c to group r
     sgr = size(find(groups(zr1(1),:)));   %how big is group r
     % Figure out what side to add to
     cgr = groups(zr1(1),1:sgr(2));
     [mindg,inddg] = min([dist(cgr(1),c) dist(cgr(sgr(2)),c)]);
     if inddg == 2
       groups(zr1(1),sgr(2)+1) = c;
     else
       groups(zr1(1),1:sgr(2)+1) = [c groups(zr1(1),1:sgr(2))];
     end
   else  %both c and r belong to groups, add group c to group r
     sgr = size(find(groups(zr1(1),:)));  %size of group r
     sgc = size(find(groups(zc1(1),:)));  %size of group c
     % Figure out what side to add to
     cgc = groups(zc1(1),1:sgc(2));  % current group c
     cgr = groups(zr1(1),1:sgr(2));  % current group r
     [mindg,inddg] = min([dist(cgc(1),cgr(1)) dist(cgc(1),cgr(sgr(2))) ...
            dist(cgc(sgc(2)),cgr(1)) dist(cgc(sgc(2)),cgr(sgr(2)))]);
     if inddg == 1
       % flop group c and add to the left of r
       groups(zr1(1),1:sgr(2)+sgc(2)) = [cgc(sgc(2):-1:1) cgr];
     elseif inddg == 2
       % add group c to the right of group r
       groups(zr1(1),sgr(2)+1:sgr(2)+sgc(2)) = cgc;
     elseif inddg == 3
       % add group c to the left of group r
       groups(zr1(1),1:sgr(2)+sgc(2)) = [cgc cgr];
     else
       % flop group c and add to the right of group r
       groups(zr1(1),1:sgr(2)+sgc(2)) = [cgr cgc(sgc(2):-1:1)];
     end
     groups(zc1,:) = zeros(1,n);
   end
 end
end
dist(r,c) = inf;
dist(c,r) = inf;
z1 = find(dist(r,:)==inf);
z2 = find(dist(c,:)==inf);
z1n = z1(find(z1~=r));
z1n = z1n(find(z1n~=c));
z2n = z2(find(z2~=c));
z2n = z2n(find(z2n~=r));
z = [z1 z2];
sz = size(z);
for j = 1:max(sz);
 for k = 1:max(sz);
   dist(z(j),z(k)) = inf;
 end
end
end

order = groups(find(groups(:,1)),:);

return;