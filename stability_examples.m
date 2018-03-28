% Stability Examples
% EE 472
clear
close all
% Jordan Normal Form Examples
% This matrix has four eigenvalues,2,2,1,1 
A = [ 12,  32,   66,  116;-25, -76, -164, -294; 21,  66,  143,  256; -6, -19,  -41,  -73];
Aeigs = eigs(A);
myJordan = jordan(A);
syms t
stab_check = expm(myJordan*t);
% You can look at the matrix exponential here and see that the first term
% corresponding the eigenvalue 1, is e^t, this is not bounded, nor does it
% trend to zero as t goes to infinity, this system is not internally stable
% We knew that anyway, just from looking at the eigenvalues, but seeing it
% in the Jordan from tells us why we can just look at the eigenvalues

% We can do lyapunov stability in Matlab as well
P_stability = lyap(A',eye(size(A)));
% Note we do A' since the matlab version solves Ax+XA'=-Q
% P_stability must be positive definate 
eigs(P_stability)
% It is not (negative eigenvalues), so our A matrix is not stable