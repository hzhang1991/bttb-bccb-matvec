function [z,nflops]=bccbmv(C,x)
% block circulant matrix with circulant blocks matvec
% C is m x n, x is (m*n) x nrhs
%
% BCCB is of size N=m*n, each blocks is of size m
% C(:,j) contains the first column of (1,j) block in the BCCB matrix
%
% Xin Ye, april 2016
[m,n]=size(C);
nrhs=size(x,2);
nflops=0;

d=[C(:,1),fliplr(C(:,2:n))];
d=fft2(d);
nflops=nflops+5*m*n*(log(m)+log(n))/log(2);

z=zeros(size(x));

for i=1:nrhs
    y=reshape(x(:,i),m,n);
    y=fft2(y); % sqrt(m*n)*F*x
    nflops=nflops+5*m*n*(log(m)+log(n))/log(2);
    
    y=d.*y; % sqrt(m*n)*D*F*x
    nflops=nflops+numel(y);
    
    y=ifft2(y); % F^H*D*F*x
    nflops=nflops+5*m*n*(log(m)+log(n))/log(2);
    
    z(:,i)=reshape(y,m*n,1);
end