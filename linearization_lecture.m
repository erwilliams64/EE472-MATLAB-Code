%% Simulate a DC Motor Model in Matlab
%  There is no control here, the model is open loop, we solve the system
%  using ilaplace and compare to the solution found using lsim, they are
%  equivalent
clear
close all
%% Method 1, Inverse Laplace
syms Ra La s Kt Ke J b t
G = (1/(Ra+La*s))*Kt*(1/(J*s^2+b*s));
H = Ke*s;
my_tf = G/(1+G*H);
% Now give actual values
my_tf = subs(my_tf,[Ra La Kt Ke J b],[1.4 .006 4.375 .00867 1 .5]); 
% To apply a step response, multiply by 1/s (the step Laplace Transform)
my_tf_step = my_tf/s;
output_time = ilaplace(my_tf_step,s,t);
t_num = [0:.1:500];
theta_out = double(subs(output_time,t,t_num));
figure
subplot(2,2,1)
plot(t_num,theta_out);
title('Inverse Laplace, Symbolic');
%% Method 2, Declare a transfer function and use step
[my_tf_num, my_tf_den] = numden(my_tf);
my_tf_num = sym2poly(my_tf_num); 
my_tf_den = sym2poly(my_tf_den); 
my_tf2 = tf(my_tf_num,my_tf_den);
subplot(2,2,2)
step(my_tf2,500)
title('Transfer Function');
% Method 3, Declare a state space representation and use lsim
[A2 B2 C2 D2] = tf2ss(my_tf_num,my_tf_den);
motor_ss = ss(A2,B2,C2,D2);
mystep = ones(1,length(t_num));
[ys2 ts2 xs2] = lsim(motor_ss,mystep,t_num,[0;0;0]);
subplot(2,2,3)
plot(ts2,ys2)
title('State Space');
% Method 4, Type in the state space representation we found by hand and use lsim
motor_ss2 = ss([0 1 0; 0 0 1; 0 -116.6667 -233.8333],[0;0;1],[729.1667 0 0],0);
[ys3 ts3 xs3] = lsim(motor_ss2,mystep,t_num,[0;0;0]);
subplot(2,2,4)
plot(ts3,ys3)
title('State Space Direct');