%% EE 472 Controlling Via Direct Pole Placement
clear 
close all

% The following model is of an inverted pendulum robot, linearized around
% the upright position

A = [0 1 0 0; 62.0193 -44.5897 0 -2123.32; 0 0 0 1; 6.09908 -10.1911 0 -485.289];
B = [0; -90.0275; 0; -20.5759];
C = eye(size(A));
D = zeros(size(A,1),1);

% There are four states
% x = [angle, angle_vel, x, x_vel]
% Where the angle is angle of the pendulum robot from vertical and x is the
% travelling distance of the wheel
% This is a SIMO system, we have one input, the voltage to the motor
% The output vector is just the states

% The goal of the system is to maintain the angle of the robot.  To this
% end, determine a set of pole locations that will allow the robot to stand
% up.  Then determine a feedback gain matrix that will achieve those poles


% Step one, controllability and observability.  We don't even know if we
% can place the poles where we want

controllableif4=rank(ctrb(A,B)); %System is controllable
observableif4=rank(obsv(A,C)); %system is observable

% Okay, now look at stability, what are the eigenvalues for the system?

stability=eigs(A); %[-529.934816422993;5.71642965208417;-5.66031322909107;0]
                   % One pole is positive so not stable
                   
% We have two eigenvalues with real parts that are non-negative, so the
% sytem is not stable.
% Verify that the system is not stable by simulating the system with no
% input, and an initial condition of [-.02; 0; 0; 0], plot the angle with
% respect to time, remember +- pi/2 means its on the ground, and the system
% will probably be unrecoverable for much less than that, remember this is
% a linearized model, so our dynamics are only accurate pretty close to
% upright
x_0=[-.02; 0; 0; 0];
t=0:.1:5;
u=zeros([length(t),1]);

ss1=ss(A,B,C,D);
[y1,t1,x1]=lsim(ss1,u,t,x_0);
figure('Name','simulated system with no controller, no input, and starting with an initial condition')
plot(t,x1);
% We haven't discussed choosing pole locations (although we did quite a bit
% in 471, but in order for the controller to compensate it must act faster
% than the system dynamics, so lets move the poles to all negative real
% parts, and perhaps 10 times as fast
poles_minseg = [-1065, -3.6+0.5*1i, -0.4, -3.6-0.5*1i]; % optimized poles
% poles_minseg=[-5000;-50;-51;-1]; %chosen poles

% Now find the gains
K=place(A,B,poles_minseg);
% Now lets simulate an impulse input and calculate the controller effort
% and the predicted response
ss2=ss(A-B*K,B,C,D);
[y2,t2,x2]=impulse(ss2,1);

figure('name','minseg angle and controller effort')
subplot(2,1,1)
plot(t2,y2(:,1))% shows alpha
subplot(2,1,2)
plot(t2,K(1)*x2(:,1))%controller effort?

%% This part takes forever to run!
% opt=gramOptions('TimeIntervals', [t(1) t(end)]);
% my_grammie=gram(ss2,'c',opt);
% ts=t;
% 
% syms t
% u_2sec = -ss2.B'*expm(ss2.A'*(ts(end)-t))*inv(my_grammie)*(expm(ss2.A*t(end))*x_0-[0;0;0;0]);
% t_num = linspace(0,ts(end),10);
% u_sym = double(subs(u_2sec,t_num));
% subplot(2,1,2)
% plot(ts,u_sym)%controller effort

%% 
% Check the eigenvalues/poles of the new state feedback system
eigs2=eigs(A-B*K);
% They check out
% In order to test this controller, lets simulate the system for an impulse
% input to the angle, and see what the controller has to do to return it
x_3=[-.02; 0; 0; 0];
[y3,t3,x3]=impulse(ss(A-B*K,B+x_3,C,D),1);
figure('name','input angle with impulse')
plot(t3,y3)

% That might not be totally fair, lets try no disturbance, just an initial
% condition of -.02 radians
[y4,t4,x4]=initial(ss(A-B*K,B,C,D),x_3);
figure('name','input angle with no impulse')
plot(t4,y4) %stabalizing at around 8 seconds

% Even with a pretty small disturbance we need to apply 500 volts to our
% little lego motor

% Experiment with the gains, try just a stabilizing controller, check the
% inputs required, increase the gains until you get a reasonable voltage,
% around 10V or so (USB power supply)

%% Stabilizing controller

% Look at the time it takes this controller to stabilize, it is probably
% much too slow, it is not much faster than the falling dynamic itself


%% Pole Values experimentally found to be ok


% Will the robot balance
save