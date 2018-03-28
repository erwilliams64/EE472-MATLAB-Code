function [dxdt]= odefcn(t,x)
m=.9048;
r=.03;
g=9.81;
J=.0676;
J_b=.000326;

u=.002.*cos(2*pi.*t);

dxdt(1)=x(2);
dxdt(2)=(m/(J_b/(r^2)+m))*(x(1)*x(4)^2-(g*sin(x(3))));
dxdt(3)=x(4);
dxdt(4)=((-2*m*x(1)*x(2)*x(4))-(m*g*x(1)*cos(x(3)))+u)/((m*(x(1)^2))+J+J_b);

dxdt=dxdt';
end

