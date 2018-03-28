close all
clear

%P1
%a

x_0=[0 0 0 0];
t= 0:.001:1;

[t,x]= ode45(@odefcn,t,x_0);

figure()
subplot(2,2,1)
plot(t,x(:,1));
title('Non-Linear Zero cond. Pos.')

subplot(2,2,2)
plot(t,x(:,2));
title('Non-Linear Zero cond. Vel.')

subplot(2,2,3)
plot(t,x(:,3));
title('Non-Linear Zero cond. Ang.')

subplot(2,2,4)
plot(t,x(:,4));
title('Non-Linear Zero cond. Torque')

%b
m=.9048;
r=.03;
g=9.81;
J=.0676;
J_b=.000326;


u=.002.*cos(2*pi.*t);

A=[0 1 0 0; 0 0 -m/(J_b/(r^2)+m)*g 0; 0 0 0 1; -m*g/(m*x_0(1)^2+J+J_b) 0 0 0];
B=[0;0;0; 1/(m*x_0(1)^2+J+J_b)];
C=[1 0 0 0];
D=[0];

my_sys =ss(A,B,C,D);
[y2 t2 x2] = lsim(my_sys, u, t, x_0);

figure()
subplot(2,2,1)
plot(t2,x2(:,1));
title('Linear Zero cond. Pos.')

subplot(2,2,2)
plot(t2,x2(:,2));
title('Linear Zero cond. Vel.')

subplot(2,2,3)
plot(t2,x2(:,3));
title('Linear Zero cond. Ang.')

subplot(2,2,4)
plot(t2,x2(:,4));
title('Linear Zero cond. Torque')

%c
x_0=[.25 0 .01 0];
t= 0:.001:1;

[t,x3]= ode45(@odefcn,t,x_0);
[y3 t3 x4] = lsim(my_sys, u, t, x_0);

figure()
subplot(3,1,1)
plot(t,x3(:,1));
title('Non-Linear w/ IC Pos.')

subplot(3,1,2)
plot(t3,x4(:,1));
title('Linear w/ IC Pos.')

subplot(3,1,3)
plot(t3,x4(:,1)-x3(:,1))
title('Difference w/ IC Pos.')

%P2
%a
load('something.mat');
%b
syms Ra La s Kt Ke J b t
G = (1/(Ra+La*s))*Kt*(1/(J*s^2+b*s));
H = Ke*s;
my_tf = G/(1+G*H); %#ok<NASGU>

my_tf=subs(G,[Ra, La, Kt, Ke, J, b], [4.4,.004,.47,.495,3.58*10^-6,.01]);
[my_tf_num,my_tf_den]=numden(my_tf);
tfnum=sym2poly(my_tf_num);
tfden=sym2poly(my_tf_den);
[A2, B2, C2, D2] = tf2ss(tfnum,tfden);
motor_ss=ss(A2, B2, C2, D2);

[pos2,t2,vel2]=lsim(motor_ss, vin, tsim);

figure()
subplot(2,1,1);
plot(t2,pos*2*pi);
title('measured position');

subplot(2,1,2);
plot(vel);
title('measured velocity');

figure()
subplot(2,1,1);
plot(t2,pos2);
title('modeled position');

subplot(2,1,2);
plot(diff(pos2*2*pi)/diff(t2));
title('modeled velocity');

%The main differance between the two velocity plots is that the measured
%system cannot physically act like the beautiful square wave velocity
%function predicted. There is also noise, though it is not that bad.

save
