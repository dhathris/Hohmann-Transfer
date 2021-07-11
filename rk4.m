%% Fixed-step Runge-Kutta 4th order integrator
function [t,x] = rk4(h,first,last,x0,MU)
% Fixed-step Runge-Kutta 4th order integrator
x=zeros((last-first)/h,size(x0,2));
t=zeros((last-first)/h,1);
x(1,:)=transpose(x0);
t(1,1)=first;
index=1;
step=first+h;
while step <= last
    xdot=two_body_eom(step,transpose(x(index,:)),MU);
    k1=h*xdot;
    xdot=two_body_eom(step+(h/2),transpose(x(index,:))+(k1/2),MU);
    k2=h*xdot;
    xdot=two_body_eom(step+(h/2),transpose(x(index,:))+(k2/2),MU);
    k3=h*xdot;
    xdot=two_body_eom(step+h,transpose(x(index,:))+k3,MU);
    k4=h*xdot;
    index=index+1;
    x(index,:)=transpose(transpose(x(index-1,:))+(1/6)*(k1+2*k2+2*k3+k4));
    t(index)=step;
    step=step+h;
end
end