close all
clear

%% 4.30 Consider the two degree of freedom system in Figure P4.30.
%  Write the equations of motion in vector form and compute the following
m1=1; %kg
m2=4; %kg
k1=240; %N/m
k2=300; %N/m
M=[m1, 0;0, m2]; 
K=[k1+k2,-k2;-k2,k2];
Ktil=M^(-1/2)*K*M^(-1/2);
[P,d]=eig(Ktil);
Lamd=transpose(P)*Ktil*P;
%a natural frequencies
wn1=sqrt(Lamd(1,1));
wn2=sqrt(Lamd(2,2));
%b mode shapes
D=inv(M)*K;
[v,d]=eig(D);
u1=v(:,1);
u2=v(:,2);
%c eigenvalues
lamda1=d(1);
lamda2=d(4);
%d eigenvectors
[P2,d2]=eigs(Ktil);
v1=P2(:,1);
v2=P2(:,2);
%e show mode shapes are not orthogonal
orthif0_1=dot(lamda1,lamda2);
%f show eigenvectors are orthogonal
orthif0_2=dot(v1,v2);
%g show modeshapes and eigenvectors are related by M^(-1/2)
d_confirm=M^(-1/2)*P2 %this isn't right
%h write the equations of motion in modal coordinates



%% 4.59

%% 4.62

%% 4.68

%% 4.71


save
