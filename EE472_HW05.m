close all
clear

%% P1
A=[0,1;1,0];
B=[-1,1]';
C=[0,1];

%a
Aeigs = eigs(A);
myJordan = jordan(A);
syms t
stab_check = expm(myJordan*t); %there is a matrix exponential that does not tend towards zero so this system is not asymptotically stable

%b
H=C*expm(A*t)*B;
H_s=laplace(H);
Poles=poles(H_s); %1 pole at -1 shows system is BIBO stable

%c
[minsys,u]=minreal(ss(A,B,C,0));
[A2,B2,C2,D2]=ssdata(minsys);

myJordan2 = jordan(A2);
syms t
stab_check2 = expm(myJordan2*t);%well this shows a matrix exponential that tends towards zero! yay

H2=C2*expm(A2*t)*B2;
H2_s=laplace(H2);
Poles2=poles(H2_s); %1 pole at -1 shows system is BIBO stable

%% P2
A3=[8 -5 10;0 -1 1; -8 5 -9];
B3=[0 0; 0 1; 1 0];
C3=[1 0 0];
sys3=ss(A3,B3,C3,0);
%a is system unstable?
A3eigs = eigs(A3); %result is -2, 2i, -2i meaning system is undamped
%bounded unstable

%b
t3=linspace(0,10,1000);
u3=[1;2];
u3_t=u3.*ones(2,size(t3,2));
[y3,t3,x3]=lsim(sys3,u3_t,t3);

%c
figure
plot3(x3(:,1),x3(:,2),x3(:,3))
%peridocity suggests bounded 
%%
save