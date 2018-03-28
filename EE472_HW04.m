clear 
close all

%% P1

A1=[8,-5,10;0,-1,1;-8,5,-9];
B1=[-1,0,1]';
C1=[1,-2,4];
D1=[0];
T=[5,0,-3;0,2,-2;-4,2,2];
sys1=ss(A1,B1,C1,D1);

x_0=[1,2,0]';
tsys1=ss2ss(sys1,T);
xt_0=inv(T)*x_0;

confirmA1=T*A1*inv(T); % Why is this different?

figure 
subplot(2,1,1);
step(sys1);
subplot(2,1,2);
step(tsys1);

% figure
% [y1_out, t1_out, x1_out] = lsim(tsys1,[1,0,0].*ones(3,1001)',0:.01:10);
% Above should be plot of first state through time for both systems


%% P2
load p2data(1);


% Is this system controllable?
cmat_rank = rank(ctrb(ss_2mts.A,ss_2mts.B));

opt=gramOptions('TimeIntervals', [t(1) t(end)]);
my_grammie=gram(ss_2mts,'c',opt);
ts=t;

syms t
u_2sec = -ss_2mts.B'*expm(ss_2mts.A'*(ts(end)-t))*inv(my_grammie)*(expm(ss_2mts.A*t(end))*xdes(:,1)-xdes(:,2));
t_num = linspace(0,ts(end),1000);
u_sym = double(subs(u_2sec,t_num));

[y_out, t_out, x_out] = lsim(ss_2mts,u_sym,t_num,xdes(:,1));

figure
plot(t_num,x_out(:,1),t_num,x_out(:,2),t_num,u_sym);
title('Controller Effort');
xlabel('Time')
legend('State 1','State 2','Input');


%% P3
A3=[0,1,0,0;0,0,1,0;0,0,0,1;-250,-255,-91,-15];
B3=[0,0,0,1]';
C3=[25,8,1,0];
p3ss=ss(A3,B3,C3,0);

sz=size(A3);
P=zeros(sz(1));
Q=zeros(sz(1));

%Making controlability matrix
for i=1:sz(1);
    P(:,i)=A3^i*B3;
end
Prank=rank(P);

%Making observability matrix
for i=1:sz(1);
    Q(:,i)=C3*A3^i;
end
Qrank=rank(Q);
mrp3ss=minreal(p3ss);

figure 
subplot(2,1,1);
step(p3ss);
subplot(2,1,2);
step(mrp3ss);

save