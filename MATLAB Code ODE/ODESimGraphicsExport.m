y0 = [1,1,0.5];
r1 = 0.5;
r2 = 1;
b = 1e-3;
mu = 0.2;

[t,y] = ode45(@(t,y) f(t,y,r1,r2,b,mu),[0 10^7],y0);

propf1 = y(:,1)./((y(:,1) + y(:,2)));
d = propf1;
propf2 = y(:,2)./((y(:,1) + y(:,2)));
prop1 = (propf1+y(:,3))./(y(:,1)+y(:,2)+y(:,3) + 1-y(:,3));
prop2 = (propf2+(1-y(:,3)))./(y(:,1)+y(:,2)+y(:,3) + 1-y(:,3));

h = subplot(2,1,1);

%figure(1)
plot(t,prop1,'linewidth',1.5)
hold on;
plot(t,prop2,'linewidth',1.5)
%plot(t,(y(:,3)))
legend("Trait 1, r_{1} = 0.5","Trait 2, r_{2} = 1",'Location','best')
title("Subplot 1: trait 1 dominates")
xlabel("Time (t)")
ylabel("Proportion of total population")
%legend("$\frac{f_{1}}{f_1+f_2}$","$\frac{m_{1}}{m_1+m_2}$",'Interpreter','latex')


% subplot(3,1,2)
% y0 = [1,1,0.5];
% r1 = 0.9;
% r2 = 0.5;
% b = 1e-3;
% mu = 0.2;
% 
% [t,y] = ode45(@(t,y) f(t,y,r1,r2,b,mu),[0 10^7],y0);
% propf1 = y(:,1)./((y(:,1) + y(:,2)));
% propf2 = y(:,2)./((y(:,1) + y(:,2)));
% prop1 = (propf1+y(:,3))./(y(:,1)+y(:,2)+y(:,3) + 1-y(:,3));
% prop2 = (propf2+(1-y(:,3)))./(y(:,1)+y(:,2)+y(:,3) + 1-y(:,3));
% 
% figure(2)
% plot(t,prop1,'linewidth',1.5)
% hold on;
% plot(t,prop2,'linewidth',1.5)
% plot(t,(y(:,3)))
% legend("Trait 1, r_{1} = 0.9","Trait 2, r_{2} = 0.5",'Location','best')
% title("Subplot 2: trait 2 dominates")
% xlabel("Time (t)")
% ylabel("Proportion of total population")

subplot(2,1,2)
y0 = [0.5,0.5,0.5];
r1 = 0.6;
r2 = 0.4;
b = 1e-3;
mu = 0.2;


[t,y] = ode45(@(t,y) f(t,y,r1,r2,b,mu),[0 10^7],y0);
propf1 = y(:,1)./((y(:,1) + y(:,2)));
propf2 = y(:,2)./((y(:,1) + y(:,2)));
prop1 = (propf1+y(:,3))./(y(:,1)+y(:,2)+y(:,3) + 1-y(:,3));
prop2 = (propf2+(1-y(:,3)))./(y(:,1)+y(:,2)+y(:,3) + 1-y(:,3));

%figure(3)
plot(t,prop1,'linewidth',1.5)
hold on;
plot(t,prop2,'linewidth',1.5)
%plot(t,(y(:,3)))
legend("Trait 1, r_{1} = 0.6","Trait 2, r_{2} = 0.4",'Location','best')
title("Subplot 2: coexistence")
xlabel("Time (t)")
ylabel("Proportion of total population")

%print(gcf, '-dpdf', 'test.pdf')
myAxes=findobj(h,'Type','Axes');
exportgraphics(myAxes,'2TraitSim.pdf');


function dydt = f(t,y,r1,r2,b,mu)
    dydt = zeros(3,1);
    dydt(1) = (b/2)*((1-r1)*y(1)*y(3)     + (1-r2)*y(2)*y(3)-2*r1*y(1)^2       + (1-r1)*y(1) - 2*r2*y(2)*y(1) );
    dydt(2) = (b/2)*((1-r2)*y(2)*(1-y(3)) + (1-r1)*y(1)*(1-y(3)) - 2*r2*y(2)^2 + (1-r2)*y(2) - 2*y(1)*y(2)*r1 );
    dydt(3) = (b/2)*(-r1*y(1)*y(3)-r2*y(2)*y(3) + r1*y(1));
end