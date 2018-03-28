%% EE/ME 472 Mini Lab 05, Full State Feedback Control with Observer on hardware! NOW with an Integrator!
% This lab takes place with two parts.  We will use our motor model we
% found from Mini Lab 01 and design a full state feedback position
% controller. We will compare this to what is achievable using
% straight proportional or PD control.  We will also implement our
% controller on the arduino Minseg motor with a full state observer.
clear
close all
%% Here is the motor model from Mini Lab 01
Ra = 4.4;
La = .004;
Kt = .47;
Ke = .495;
J = 3.58e-6; 
b = .01;
% Here is the transfer function for voltage (input) to position (radians)
motor_tf = tf(Kt,[J*La Ra*J+La*b Ra*b+Kt*Ke 0]);
[A, B, C, D] = tf2ss(Kt,[J*La Ra*J+La*b Ra*b+Kt*Ke 0]);
motor_ss = ss(A,B,C,D);
Ts = 0.003;
motorDiscrete = c2d(motor_ss,Ts);
% First Test the observability and controllability of the motor
ctr_rank = rank(ctrb(motor_ss));
obs_rank = rank(obsv(motor_ss));
% These both check out as rank 3, so we are okay
% Lets check the stability as well
open_loop_eigs = eigs(motor_ss.A);
open_loop_eigs_discrete = eigs(motorDiscrete.A);
% We do have a pole at zero, but the other two eigenvalues all have
% negative real parts

%% Problem 1-2 dlqr, Testing on the model
% We will use dlqr to choose the gains
Q= motor_ss.C'*2*motor_ss.C; % this gives us a nice large value for Q
R=.4;
K= dlqr(A,B,Q,R);

% First Load the desired reference (provided)
load reference
% Now build a linear model, you don't have to include the observer
ss_feedback = ss(motorDiscrete.A-motorDiscrete.B*K,motorDiscrete.B,motorDiscrete.C,motorDiscrete.D,Ts);
step(ss_feedback)
% Check the eigenvalues/poles of the new state feedback system
newEigs = eigs(ss_feedback.A);

% Lets scale the reference as per the notes
Nbar=-inv(C*inv(motor_ss.A-motor_ss.B*K)*motor_ss.B);

% The following contains data from the proportional controller with a
% feedback gain of 12, it is nice for comparison
load RealMotorData003_PF12.mat % Measured Position for the proportional feedback controller

% Observer, we still need to use place and implement an observer, this was
% all part of MiniLab04, here is my observer code.
poles_desired = [-20.5 -1946.64804469274+3940.77346066813*1i -1946.64804469274-3940.77346066813*1i];
obs_poles = poles_desired.*2; % Say 10 times as fast
pDiscObs=exp(obs_poles.*Ts);  %poles in z-domain
L = place(motorDiscrete.A',motorDiscrete.C',pDiscObs)';
A_big = [motorDiscrete.A -motorDiscrete.B*K; L*motorDiscrete.C motorDiscrete.A-motorDiscrete.B*K-L*motorDiscrete.C];    %[A-B*K zeros(size(A)); L*C-B*K A-L*C];
B_big = [motorDiscrete.B*Nbar;motorDiscrete.B*Nbar];
C_big = [motorDiscrete.C zeros(size(motorDiscrete.C))];

%% Problem 3 Implement the dlqr gains on the MINSEG
% Load the results from the observer and full state feedback (with LQR gains) implemented on
% the MinSeg
load FSB_LQR_Minseg.mat
%% Plotting (Problems 1-4)
% Here is the model simulated with the observer, we want to test it with
% the gains we found using dlqr
BigSys = ss(A_big,B_big,C_big,0,Ts);
[y_testobs, t_testobs, x_testobs] = lsim(BigSys,ref,tsim,[0 0 0 0 0 0]);
x_cell = num2cell(x_testobs(:,1:3),2);
mulfunc = @(x) K*x';
x_scaled = cellfun(mulfunc,x_cell);
model_input = ref.*Nbar-x_scaled;
figure('Name','Full State Feedback Model versus Reality');
subplot(1,2,1)
plot(t_testobs,y_testobs,tsim,ref,tsim,pos);
legend('LQR-G Model','Reference','LQR-G Measured');
subplot(1,2,2);
plot(tsim,vin,tsim,model_input)
legend('Actual Input','Model Input');
xlabel('Time');
ylabel('Volts');

%% Problem The steady state result did not work too well. Lets try an LQI
% The model prediction shows a good steady state agreement, but it doesn't
% work so well in reality, lets try an LQI which is more robust to modeling
% errors

%% Continuous Model
% I have had some difficulty getting this to work on the
% discrete model, the dcgain command is useful to find the dc gain of a
% system, this is the model I use to try and set gains
 % Q = 
 % R = 
 % [KC, Sc, ec] = lqi(motor_ss,Q,R);
 
 Kint = -KC(4);
 K2 = KC(1:3);
 NewA = [motor_ss.A-motor_ss.B*K2 motor_ss.B*Kint; -motor_ss.C zeros(size(motorDiscrete.C,1),size(motorDiscrete.B,2))];
 NewB = [0;0;0;1];
 NewC = [motor_ss.C 0];
 motorcont_LQI = ss(NewA,NewB,NewC,0);
 motorcont_LQI_D = c2d(motorcont_LQI,Ts);
 [y_testobs, t_testobs, ~] = lsim(motorcont_LQI_D,ref,tsim,[0 0 0 0]);
 tsim_ideal=tsim;
 figure('Name','LQRI Model Prediction');
 plot(t_testobs,y_testobs,tsim,ref);
 legend('LQI Model','Reference');
%% Discrete Model LQI
% This is what I am implementing on the MinSeg!
% It produces the same poles as the continuous model, you can check using
% exp(ec.*Ts), should be the same as ed
 [Kbig Sd, ed] = lqi(motorDiscrete,Q,R);
 K=-Kbig;
 Kint = -Kbig(4);
 K2 = Kbig(1:3);

%% Plotting for problem 5 
% The proof is in the pudding
x_cell = num2cell(x_testobs(:,1:4),2);
mulfunc = @(x) -K*x';
x_scaled = cellfun(mulfunc,x_cell);
model_input = ref-x_scaled;
load FSB_LQRI_Minseg.mat
figure('Name','LQRI versus Reality');
subplot(1,2,1)
plot(t_testobs,y_testobs,tsim,ref,tsim,pos);
hold on
load RealMotorData003_PF12.mat % Fun to see if we are better
plot(tsim,pos);
legend('LQRI Model','Reference','LQRI Measured','Prop 12');
subplot(1,2,2);
load FSB_LQRI_Minseg.mat
plot(tsim,vin,tsim_ideal,model_input)
legend('Actual Input','Model Input');
xlabel('Time');
ylabel('Volts');