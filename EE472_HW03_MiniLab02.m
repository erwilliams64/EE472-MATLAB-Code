
%% Mini Lab 2
% State Space Solutions in Matlab
% Simulate the two mass translational system in Matlab, with the following
% input, time steps, and initial conditions:
% t_num = 0:.1:50;
% x_0 = [.1 0 .1 0]'; 
% u1 = 20.*[1 ones(1,length(t_num)-1)];
% u2 = 10.*[1 ones(1,length(t_num)-1)];
% Solve this system in three ways:
%   1. Using the Matlab function integral, with everything in the time
%   domain 
%   2. Using linear algebra in the Laplace domain
%   3. Using the Matlab function lsim
% For all time domain solutions use the Matlab function expm
% What are the eigenvalues of the A matrix?
% Convert the solution to modal and companion form, check the eigenvalues
% of both new A matrices and compare to the eigenvalues of the original A
% matrix

clear
close all
%% Matlab MIMO Implementation
% Two mass translation system
m1 = 40;
b1 = 20;
k1 = 400;
m2 = 20;
b2 = 10;
k2 = 200;
A_2mts = [0 1 0 0; -(k1+k2)/m1 -(b1+b2)/m1 k2/m1 b2/m1; 0 0 0 1; k2/m2 b2/m2 -k2/m2 -b2/m2];
B_2mts = [0 0; 1/m1 0; 0 0; 0 1/m2];
C_2mts = [1 0 0 0; 0 0 1 0];
D_2mts = [0 0; 0 0];
ss_2mts = ss(A_2mts,B_2mts,C_2mts,D_2mts);

t_num = 0:.1:50;
x_0 = [.1 0 .1 0]'; 
u1 = 20;
u2 = 10;
u_t=[u1.*ones(1,length(t_num));u2.*ones(1,length(t_num))];


%% LSIM
[y_2mts t_2mts x_2mts]=lsim(ss_2mts, u_t,t_num);
figure
plot(t_2mts,y_2mts(:,1),t_2mts,y_2mts(:,2));
legend('m1','m2')

%% Integral
% Force ourselves to do this once in the time domain using integral, then probably never again
% t_cur=zeros(1,length(t_num));
% f=@(tao,u_t) C_2mts*exp(A_2mts*(t_cur-tao))*B_2mts*u_t;
% 
% 
% for i=1:length(t_num);
%     t_cur=t_num(i,1);
%     y_t(1,i)=C_2mts*exp(A_2mts*(t_cur-t_num(1)))*x_0+integral(f,t_num(1),t_cur)+D_2mts*u_t;
%     i=i+1;
% end

%% Laplace
% X(s) = (sI-A)^(-1)x_0+(sI-A)^(-1)*B*U(s)

syms s t;
X_s=inv(s*eye(4,4)-A_2mts)*x_0+inv(s*eye(4,4)-A_2mts)*B_2mts*laplace(t/t*u_t);
x_t=ilaplace(X_s);
Y_s=C_2mts*X_s+D_2mts*laplace(t/t*u_t);
y_t2=ilaplace(Y_s);
plot(t_num,subs(y_t2(1),t,t_num),t_num,subs(y_t2(2),t,t_num));




%% Convert this system to both Controller Canonical and Modal Form then compare eigenvalues

modal_sys = canon(ss_2mts,'modal');
CC_sys = canon(ss_2mts,'companion');
eigs_modal = eigs(modal_sys.A);
eigs_cc = eigs(CC_sys.A);
eigs_original = eigs(ss_2mts.A);

save;
