function [m_out dim_shift]=fliprotND(m_in,curr_order, new_order)

% function rotates, flips N dimentional matrix from current dimention/direction order
% to new dimention/direction order

% for 3 D matrixes:
% up Y  :  1
% down Y : -1
% rigth : -2
% left : 2
% up Z : 3
% down Z : -3

% for example to reorientd 3 D stack with orientation:
% posterior is rigth, lateral is down Y,pia is up Z.
% to
% pia is up Y,posterior is up Z,lateral is left.
% one need to specify curr_order=[-2 -1 3]  new_order=[3 2 1] (using
% arbitrary order [post lat pia])

dim_shift=zeros(size(curr_order));
flip_key=zeros(size(curr_order));
for dim=1:length(curr_order)
    a=find(abs(curr_order)==abs(new_order(dim)));
    dim_shift(dim)=a;
%    flip_key(dim)=curr_order(a)/new_order(dim);
    flip_key(dim)=sign(curr_order(dim)/new_order(dim));
end
% dim_shift
% flip_key
    m_out = permute(m_in,dim_shift);

if length(flip_key)>0 && flip_key(1)==-1
    m_out(:,:,:,:,:)=m_out(end:-1:1,:,:,:,:);
end

if length(flip_key)>1 && flip_key(2)==-1
    m_out(:,:,:,:,:)=m_out(:,end:-1:1,:,:,:);
end

if length(flip_key)>2 && flip_key(3)==-1
    m_out(:,:,:,:,:)=m_out(:,:,end:-1:1,:,:);
end

if length(flip_key)>3 && flip_key(4)==-1
    m_out(:,:,:,:,:)=m_out(:,:,:,end:-1:1,:);
end

if length(flip_key)>4 && flip_key(5)==-1
    m_out(:,:,:,:,:)=m_out(:,:,:,:,end:-1:1);
end




