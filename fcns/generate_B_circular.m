function [B1,B2]=generate_B_circular(M,N)

Dx=sparse(eye(N));
Dx_1=-ones(N-1,1);
Dx_1=sparse(diag(Dx_1,-1));
Dx(1,N)=-1;
Dx=Dx+Dx_1;
IM=sparse(eye(M));
B1=kron(Dx,IM);%x

Dy=sparse(eye(M));
Dy_1=-ones(M-1,1);
Dy_1=sparse(diag(Dy_1,-1));
Dy(1,M)=-1;
Dy=Dy+Dy_1;
IN=sparse(eye(N));
B2=kron(IN,Dy);
%B=[B1;B2];
