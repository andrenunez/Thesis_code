r1 = 0.5;
r2 = 0.1;
d1 = 0;
d2 = 0;
b = 3;
tmax = 10^5;
t_cond = [0,10^5];
init_cond = [1,1,0.5];

[t,y] = ode45(@(t,y) f(t,y,r1,r2,d1,d2,b),t_cond,init_cond);
prop = (y(:,1)+y(:,3))./(1+ y(:,1) + y(:,2));
plot(t,prop)


function dydt = f(t,y,r1,r2,d1,d2,b)
    dydt = zeros(3,1);
    dydt(1) = 0.5*(y(2)*(-2*b*r2*y(2) + b*(r2-1)*(y(3)-2) + 2*d1 - 2*d2)     + b*y(1)*((r1-1)*(y(3)-1) - 2*r1*y(2)));
    dydt(2) = 0.5*(y(1)*(b*(-r1*(2*y(1) + y(3) + 1) + y(3) + 1) + 2*d1-2*d2) + b*y(2)*((1-r2)*y(3) - 2*r2*y(1)));
    dydt(3) = 0.5*b*(r1*y(1)*(y(3)-1) + r2*y(2)*y(3));
end