function [ out ] = state_trans_integral(tau_in,tf)
% Time domain solution to be integrated
% it was tricky to work the input in using a one line function, so I am 
% declaring this one as its own function
syms tau
U = [20;10];
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
% We have to calculate the state transition matrix, lets cheat and let
% Matlab do it, this requires the symbolic toolbox
% There are some ways to do this using all the eigenvalues or a spectral
% decomposition, usually we don't care because we can just solve for it in
% the Laplace domain
state_trans = expm(A_2mts*(tf-tau));
out = double(C_2mts*subs(state_trans,tau,tau_in)*B_2mts*U);

end

