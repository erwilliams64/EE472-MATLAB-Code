%% Exam 1 Williams EE 472

clear
close all

%% P1
%a
A1=[0 1;-13,-4];
B1=[-1 0;0 2];
C1=[1 0];
D1=[0];
sys1=ss(A1,B1,C1,D1);

%b
syms t
x1_0=[1 0.5]';
t1=0:.001:1;

u1= [cos(.02*pi()*t);exp(-.002*t)];
u1_sym = double(subs(u1,t1));

figure
[y1_out, t1_out, x1_out] = lsim(sys1,u1_sym,t1,x1_0);
plot(t1,y1_out(:,1));
title('Linear w/ IC Pos.')

%c
[t1_out,y1_out2]= ode45(@odefcn2,t1,x1_0);

figure
plot(t1_out,y1_out2(:,1),t1,y1_out(:,1))
title('Non-Linear w/ IC Pos.')


%% P2
load p2system.mat

%a

% There are 5 states, 2 inputs, and 3 outputs by inspection. These were
% deterimed by the size of the A, B, and C matricies respectively

%b
t2=0:.1:1;
u2=[exp(-2*t2);exp(-2*t2)];
x2_0=[0 0 0 0 0]';
u2_t=[u2(1).*ones(1,length(t2));u2(2).*ones(1,length(t2))];

syms s t;
X2=inv(s*eye(5,5)-p2system.A)*x2_0+inv(s*eye(5,5)-p2system.A)*p2system.B*laplace(t/t*u2_t);
x2=ilaplace(X2);
Y_s=p2system.C*X2+p2system.D*laplace(t/t*u2_t);
y_t2=ilaplace(Y_s);

figure
plot(t2,double(subs(y_t2(1),t,t2)),t2,double(subs(y_t2(2),t,t2)),t2,double(subs(y_t2(3),t,t2)));

%c
controllableif5=rank(ctrb(p2system.A,p2system.B))
observableif5=rank(obsv(p2system.A,p2system.C))
% yes system is controllable and observable

%% P3
A3=[0 1;-8,-6];
x3_0=[1 -1]';
u3=[0;0];
t3=0:.1:6;

sys3=ss(A3,[0,0]',[0,0],[0]); %The values of B,C,D do not matter

figure
[y3,t3,x3]=initial(sys3,x3_0,t3);
plot(t3,x3(:,1),t3,x3(:,2));

save

