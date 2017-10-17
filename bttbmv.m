function [z,nflops]=bttbmv(T,x)
% block Toeplitz with Toeplitz blocks matvec
% T is (2*m-1) x (2*n-1), x is (n*m) x nrhs
%
% BTTB is of dimension N=m*n each block is of size m*m
%
% T(:,j) is the vector (t_-(m-1),...,t_-1,t_0,t_1,...,t_m-1)
% generating T_j-n
%
% Xin Ye, April 2016

[m,n]=size(T);
m=(m+1)/2;
n=(n+1)/2;

nrhs=size(x,2);

% extend to a BCCB matrix of size (2*m-1)*(2*n-1)
C=[T(m:2*m-1,:);T(1:m-1,:)];

C=fliplr(C);
C=[C(:,n:2*n-1),C(:,1:n-1)];

v=zeros((2*m-1)*(2*n-1),nrhs);
for i=1:n
    v((i-1)*(2*m-1)+1:(i-1)*(2*m-1)+m,:)=x((i-1)*m+1:i*m,:);
end

[v,nflops]=bccbmv(C,v);

z=zeros(size(x));
for i=1:n
    z((i-1)*m+1:i*m,:)=v((i-1)*(2*m-1)+1:(i-1)*(2*m-1)+m,:);
end