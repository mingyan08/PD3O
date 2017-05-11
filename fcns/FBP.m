function [u0,F0]=FBP(P,N,option)

[nr,nt]=size(P);
M=N;
if nargin<3
    option='Ram-Lak';
end
angle=[0:nt-1]*180./nt;
u0=flipud(iradon(P,angle,option));
nsha=(size(u0,1)-M)/2;
nshb=(size(u0,2)-N)/2;
u0=u0(nsha+1:nsha+N,nshb+1:nshb+N);

u0=u0([end,1:end-1],[end,1:end-1]); 

if nargout>1
    F0=fftshift(fft2(fftshift(u0)));
end