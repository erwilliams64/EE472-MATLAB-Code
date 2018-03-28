function [dxdt]= odefcn2(t,x1)
u= [cos(.02*pi()*t);exp(-.002*t)];

dxdt(1)=x1(2)-u(1);
dxdt(2)=-2*x1(2)-2*x1(1)*x1(2)-4*x1(1)^3+2*u(2);

dxdt=dxdt';
end

