close all
clear

%% P1
syms t s
A=[0,1;-1,2];
B=[0;1];
C=[2, 1];
D=[0];

x0=[1;-2];

solution=initial(

sol=inv(s*eye(2)-A);
U=laplace(exp(-2*t));
X=sol*x0+sol*B*U;
x=ilaplace(X);
t1=0:.01:10;
x_soln=subs(x,t,t1);
plot(t1,x_soln(1,:))


%% P2 
m1=1;
m2=2;
m3=3;
b1=.4;
b2=.8;
b3=1.2;
b4=1.6;
k1=10;
k2=20;
k3=30;
k4=40;

u2=[3.*ones(1,1000);2.*ones(1,1000);1.*ones(1,1000)];

A=[0 1 0 0 0 0; -(k1+k2)/m1, -(b1+b2)/m1,k2/m1,b2/m1,0 0;0 0 0 1 0 0;k2/m2,b2/m2,-(k2+k3)/m2,-(b2+b3)/m2, k3/m2, b3/m2; 0 0 0 0 0 1; 0 0, k3/m3, b3/m3, -(k3+k4)/m3, -(b3+b4)/m3];
B=[0 0 0; 1/m1 0 0; 0 0 0; 0 1/m2 0; 0 0 0; 0 0 1/m3];
C=[0 1 0 0 0 0; 0 0 0 1 0 0; 0 0 0 0 0 1];
D=[0];

t2=linspace(0,50,1000);
x_0=[0;0;0;0;0;0];

mysys=ss(A,B,C,D);
[ysim,tsim,xsim]=lsim(mysys,u2,t2,x_0);
plot(t2,ysim(:,1),t2,ysim(:,2),t2,ysim(:,3));
legend('y1','y2','y3');



save