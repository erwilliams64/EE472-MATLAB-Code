%% EE/ME 472 Mini Lab 04, Full State Feedback Control with Observer on hardware!
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
%% First Test the observability and controllability of the motor
controllableif3=rank(ctrb(A,B)); %System is controllable
observableif3=rank(obsv(A,C)); %system is observable
%% Lets check the stability as well
stability=eigs(A); %all eigenvalues are negative so stable
%% Lets calculate the gains for FSB
% First choose the desired poles, lets be nice and conservative
poles_desired = [-20.5 -1946.64804469274+3940.77346066813*1i -1946.64804469274-3940.77346066813*1i]; % Continuous, PF with a gain of 12
pDisc=exp(poles_desired.*Ts); %poles in z-domain

K=place(A,B,poles_desired);
%% Problem 1a, Now find the gain
% place command, use the discrete model
Kdisk=place(motorDiscrete.A,motorDiscrete.B,pDisc); % very similar K values

%% Problem 1b Lets check our newfangled controller using our model, then we can compare it to reality
% First Load the desired reference (provided)
load reference
% Build the linear model with the FSB included
t = 0:Ts:5;
u = sin(t);
motor_ssfb=ss(A-B*K,B,C,D);
[y1,t1,x1]=lsim(motor_ssfb,u,t);
% Check the eigenvalues/poles of the new state feedback system
eigs(motor_ssfb.A); %all nice and negative

[y_test, t_test, x_test] = lsim(motor_ssfb,ref,tsim,[0 0 0]);
figure('Name','Problem 1')
subplot(1,2,1)
plot(t_test,y_test,tsim,ref)
legend('Model Output','Reference Output');
title('Model with no G Matrix');
xlabel('Time');
ylabel('Angle (radians)');
% Lets scale the reference as per the notes
G=-inv(C*inv(A-B*Kdisk)*B);
% Redefine the system with a scaled input
scaled_ss=ss(motorDiscrete.A-motorDiscrete.B*Kdisk,G*motorDiscrete.B,motorDiscrete.C,motorDiscrete.D, Ts);
% Simulate and plot results
[y_test, t_test, x_test] = lsim(scaled_ss,ref,tsim,[0 0 0]);
subplot(1,2,2)
plot(t_test,y_test,tsim,ref)
legend('Model Output','Reference Output');
title('Model with G Matrix');
xlabel('Time');
ylabel('Angle (radians)');

% Now load the results of a test run using the proportional controller in
% simulink
load motordata.mat
% Also load the proportional feedback model for comparison
load MotorModelPFB.mat % Loads IOTransfer_r2y, proportional feedback model with gain of 12

% Now simulate the model and compare it to the actual performance
t2=linspace(0,45,14977);
[y_pf,t_pf,x2]=lsim(IOTransfer_r2y,squeeze(vin),t2);
figure('Name','Proportional Feedback Model versus Reality');
plot(t_test,y_test,'.',t_pf,y_pf,tsim,squeeze(pos),tsim,ref)
legend('State Feedback Model','Proportional Feedback Model','Measured PF','Reference');

%% We have one more problem before we can implement this on the actual Minseg, if we are using only wheel encoders, we don't know all the states, in fact we only know the output.
% We could add sensors, there are rate gyros and accelerometers on the
% minseg, but instead lets take the opportunity to see if we can use an
% observer in tandem with our full state feedback controller, we will
% design the observer in Matlab, test it on the model, then try and
% implement the controller and observer on the Minseg.
obs_poles = poles_desired.*2; % Say 10 times as fast
pDiscObs=exp(obs_poles.*Ts);  %poles in z-domain
% use place on the discrete model to place the discrete poles
L=place(motorDiscrete.A',motorDiscrete.C',pDiscObs)';
% Create a linear model of the system using FSB and the observer!
Aobs=[motorDiscrete.A, -motorDiscrete.B*K;
    L*motorDiscrete.C, motorDiscrete.A-motorDiscrete.B*K-L*motorDiscrete.C];
obs_ss=ss(Aobs,[motorDiscrete.B;motorDiscrete.B],[motorDiscrete.C zeros(size(motorDiscrete.C))],motorDiscrete.D, Ts);

%% Load the results from the observer and full state feedback implemented in
% simulink on the Minseg!
load motordata_withobserver.mat;



%% Simulate the whole thing with lsim and compare to the measured
% performance
t3=linspace(0,45,15001);
[y_testobs,t_testobs,x_testobs]=lsim(obs_ss,G*squeeze(ref)',t3);
x_est=x_testobs(:,4:6)';
model_input=G*ref'-Kdisk*x_est;
figure('Name','Full State Feedback Model versus Reality');
subplot(1,2,1)
plot(t_testobs,y_testobs,tsim,ref,tsim,squeeze(pos));
legend('FSB Model','Reference','FSB Measured');
subplot(1,2,2);

plot(tsim,vin,tsim,model_input)
legend('Actual Input','Model Input');
xlabel('Time');
ylabel('Volts');
 
save;
