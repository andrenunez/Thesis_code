r1 = 0.5;
r2 = 0.9;
d1 = 1;
d2 = 2;
b = 3;

y0 = [1,1,1,1];

[t,y] = ode45(@(t,y) f(t,y,r1,r2,d1,d2,b),[0 10^5],y0);

prop = (y(:,1)+y(:,2))./(y(:,1)+y(:,2)+y(:,3)+y(:,4));
prop2 =(y(:,3)+y(:,4))./(y(:,1)+y(:,2)+y(:,3)+y(:,4));

plot(t,prop);
% hold on;
% plot(t,prop2);
% legend("Trait 1, r_{1} = 0.5,","Trait 2, r_{2} = 0.1",'Location','bestoutside')
% title("d_{m} = 2, d_{f} = 1")

function yprime = f(t,y,r1,r2,d1,d2,b)
    dydt = zeros(4,1);
    denom = 1/(y(1)+y(3));
    dydt(1) = -d1*y(1) + (b*denom)*(y(2)*y(1)*r1 +     0.5*y(2)*y(3)*r1 +     0.5*y(1)*y(4)*r2);
    dydt(2) = -d2*y(2) + (b*denom)*(y(2)*y(1)*(1-r1) + 0.5*y(2)*y(3)*(1-r1) + 0.5*y(1)*y(4)*(1-r2));
    dydt(3) = -d1*y(3) + (b*denom)*(y(3)*y(4)*r2 +     0.5*y(2)*y(3)*r1 +   0.5*y(1)*y(4)*r2);
    dydt(4) = -d2*y(4) + (b*denom)*(y(3)*y(4)*(1-r2) + 0.5*y(2)*y(3)*(1-r1) + 0.5*y(1)*y(4)*(1-r2));
    
    yprime = zeros(4,1);
    
    yprime(1) = dydt(1)/(y(1)+y(3)) - y(1)/((y(1)+y(3))^2)*(dydt(1)+dydt(3));
    yprime(2) = dydt(2)/(y(1)+y(3)) - y(2)/((y(1)+y(3))^2)*(dydt(1)+dydt(3));
    yprime(3) = dydt(3)/(y(1)+y(3)) - y(3)/((y(1)+y(3))^2)*(dydt(1)+dydt(3));
    yprime(4) = dydt(4)/(y(1)+y(3)) - y(4)/((y(1)+y(3))^2)*(dydt(1)+dydt(3));
end