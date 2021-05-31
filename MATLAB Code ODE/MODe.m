
r1 = 0.5;
r2 = 1;
b = 2;
m = 1e-6;

y0 = [1,1,1,1];

[t,y] = ode45(@(t,y) f(t,y,r1,r2,b,m),[0 10^4],y0);

prop  = (y(:,1)+y(:,2))./(y(:,1)+y(:,2)+y(:,3)+y(:,4));
prop2 = (y(:,3)+y(:,4))./(y(:,1)+y(:,2)+y(:,3)+y(:,4));
plot(t,prop)
hold on;
plot(t,prop2)
legend("Trait 1, r_{1} = 0.5","Trait 2, r_{2} = 1",'Location','best')
xlabel("Time (t)")
ylabel("Proportion of total population")

function yprime = f(t,y,r1,r2,b,m)
    dydt = zeros(4,1);
    denom = 1/(y(1)+y(3));
    dydt(1) = (1-m)*(b*denom)*(y(2)*y(1)*r1 +     0.5*y(2)*y(3)*r1 +     0.5*y(1)*y(4)*r2) + m*(b*denom)*(y(3)*y(4)*r2 +     0.5*y(2)*y(3)*r1 +   0.5*y(1)*y(4)*r2);
    dydt(2) = (1-m)*(b*denom)*(y(2)*y(1)*(1-r1) + 0.5*y(2)*y(3)*(1-r1) + 0.5*y(1)*y(4)*(1-r2)) + m*(b*denom)*(y(3)*y(4)*(1-r2) + 0.5*y(2)*y(3)*(1-r1) + 0.5*y(1)*y(4)*(1-r2));
    dydt(3) = (1-m)*(b*denom)*(y(3)*y(4)*r2 +     0.5*y(2)*y(3)*r1 +   0.5*y(1)*y(4)*r2) + m*(b*denom)*(y(2)*y(1)*r1 +     0.5*y(2)*y(3)*r1 +     0.5*y(1)*y(4)*r2);
    dydt(4) = (1-m)*(b*denom)*(y(3)*y(4)*(1-r2) + 0.5*y(2)*y(3)*(1-r1) + 0.5*y(1)*y(4)*(1-r2)) + m*(b*denom)*(y(2)*y(1)*(1-r1) + 0.5*y(2)*y(3)*(1-r1) + 0.5*y(1)*y(4)*(1-r2));
    
    yprime = zeros(4,1);
    
    yprime(1) = dydt(1)/(y(1)+y(3)) - y(1)/((y(1)+y(3))^2)*(dydt(1)+dydt(3));
    yprime(2) = dydt(2)/(y(1)+y(3)) - y(2)/((y(1)+y(3))^2)*(dydt(1)+dydt(3));
    yprime(3) = dydt(3)/(y(1)+y(3)) - y(3)/((y(1)+y(3))^2)*(dydt(1)+dydt(3));
    yprime(4) = dydt(4)/(y(1)+y(3)) - y(4)/((y(1)+y(3))^2)*(dydt(1)+dydt(3));
end