function xdot = two_body_eom(t,x,MU)
xdot=zeros(6,1);
r=norm([x(1) x(2) x(3)]);
xdot(1)=x(4);
xdot(2)=x(5);
xdot(3)=x(6);
xdot(4)=-MU/r^3*x(1);
xdot(5)=-MU/r^3*x(2);
xdot(6)=-MU/r^3*x(3);
end

