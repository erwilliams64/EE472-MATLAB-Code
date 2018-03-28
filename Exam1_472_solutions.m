%% EE 472
% Exam 1 Solutions
clear 
close all
%% Problem 1
A = [0 1 0; -1 0 -8.6; 0 0 -2/.2];
B = [0 0; 3 -9; 1 0];
C = [ 1 0 0; 0 0 1];
eigs(A);
B = [0 0; 3 -9; 1 0];
C = [1 0 0; 0 0 1];
p1ss = ss(A,B,C,0);
t_num = 0:.01:10;
u = [.07*cos(pi*t_num); 3*exp(-t_num).*cos(pi*t_num)]';
[y, t_out, x] = lsim(p1ss,u,t_num,[0;.3;4]);
figure
plot(t_out,y(:,1),t_out,y(:,2))
title('Problem 1 Outputs versus Time');
xlabel('Time');
ylabel('Outputs');
% find rank of each matrix to see control vs. obs

%% Problem 2
p2system4 = rss(5,3,2);
% save('p2system','p2system');
% load p2system
% % Is this system controllable, observable?
% % How many states, how many inputs, how many outputs 
x2_0 = [.1 .1 .1 .1 .1]';
t_num = 0:.1:3;
syms t
syms s
inputs = laplace(cos(t));
mexp = inv((s*eye(size(p2system4.A))-p2system4.A));
y_laplace = p2system4.C*mexp*x2_0+(p2system4.C*mexp*p2system4.B)*[inputs; inputs];
y2_time = ilaplace(y_laplace);
y2_num = double(subs(y2_time,t_num));
figure
plot(t_num,y2_num(1,:),t_num,y2_num(2,:),t_num,y2_num(3,:));
title('Problem 2: All Outputs');
legend('Output 1','Output 2','Output 3');
%% Problem Three
% Give a minimal realization of the following system, how many states are
% required to represent it, is it observable and controllable

p3tm = [tf([1 0],[1 -2 1]) tf(1,[1 -1]);tf(-6,[1 2 -3]) tf(1,[1 3])] ;
p3ss = ss(p3tm);
% save('p3ss','p3ss');
% load p3ss
p3obsv = rank(obsv(p3ss));
p3ctr = rank(ctrb(p3ss));
[p3minreal, U] = minreal(p3ss); %U=inv(T)
t3 = linspace(0,1,1000)';
u3 = [cos(4*pi*t3) sin(4*pi*t3)];
[y3, t_out, x3] = lsim(p3minreal,u3,t3,[0;0;0]);
figure
subplot(1,2,1)
plot(t_out,y3(:,1),t3,x3(:,2))
title('Problem 3 MinReal');
xlabel('Time');
legend('First Output','Second State');
subplot(1,2,2)
[y3, t_out, x3] = lsim(p3ss,u3,t3,[0;0;0;0;0;0]);
plot(t_out,y3(:,1),t_out,x3(:,2))
title('Problem 3 Original');
xlabel('Time');
legend('First Output','Second State');

% Of the original six states how many are controllable and observable, how
% many are controllable but not observable, how many are observable but not
% controllable, and how many are neither controllable, nor observable
% For this lets find the Kalman Decomposition

% Now I can use the following to get the Kalman Decomposition
% Remember the inverse of an orthogonal matrix is the same as its transpose
Ak = U*p3ss.A*U';
Bk = U*p3ss.B;
Ck = p3ss.C*U';
Dk = 0;
Ak(abs(Ak)<1e-10)=0;
Bk(abs(Bk)<1e-10)=0;
Ck(abs(Ck)<1e-10)=0;
% Ramp up i, until you get a drop in either rank
i=6;
numctr = rank(ctrb(Ak(1:i,1:i),Bk(1:i,:)))
numobs = rank(obsv(Ak(1:i,1:i),Ck(:,1:i)))
% 4th state is controllable but not observable
% 5th state is controllable but not observable
% 6th state is observable but not controllable
% No states are neither controllable nor observable;

%% Problem 4
% Create a state space system for the following nonlinear function.
% Linearize around x = [1 0];




