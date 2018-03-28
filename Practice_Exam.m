%% Exam 1 Practice

clear
close all

%% P1
A1=[0 1 0; -1 0 -8.6; 0 0 -10];
B1=[0 0;3 -9; 1 0];
C1=[1 0 0;0 0 1];
D1=[0];
x1_0=[0 0.3 4]';

syms t
u1= [.07.*cos(pi()*t);3.*exp(-t).*cos(pi()*t)];
sys1=ss(A1,B1,C1,D1);
t1=0:.01:10;
u1_sym = double(subs(u1,t1));

%plotting outputs
[y1_out, t1_out, x1_out] = lsim(sys1,u1_sym,t1,x1_0);
plot(t1,y1_out(:,1),t1,y1_out(:,2))


%% P2
p2system = rss(5,3,2);
% % Is this system controllable, observable?

% % How many states, how many inputs, how many outputs 
x2_0 = [.1 .1 .1 .1 .1]';
t_num = 0:.1:3;


X_s=inv(s*eye(4,4)-A_2mts)*x_0+inv(s*eye(4,4)-A_2mts)*B_2mts*laplace(t/t*u_t);
x_t=ilaplace(X_s);
Y_s=C_2mts*X_s+D_2mts*laplace(t/t*u_t);
y_t2=ilaplace(Y_s);
plot(t_num,subs(y_t2(1),t,t_num),t_num,subs(y_t2(2),t,t_num));
%% P3

%% P4

save;
