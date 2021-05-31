%Cases: %

% Simulation 1

                    % m1,f1 ,m2, f2
% Initial Conditions [5, 5 ,15, 15]
% r1 = 0.5, r2 = 0,0.1,...,1
% d_males = 1, d_females = 2
% time steps = 10^5

%Result, trait 1 dominates regardless unless, trait 2 = 0.5, then we
%converge to the initial condition%

% Simulation 2

                    % m1,f1 ,m2, f2
% Initial Conditions [15, 15 ,5, 5]
% r1 = 0.5, r2 = 0,0.1,...,1
% d_males = 1, d_females = 2
% time steps = 10^5

%Result, trait 1 dominates regardless unless, trait 2 = 0.5, then we
%converge to the initial condition%

% Simulation 3
                    % m1,f1 ,m2, f2
% Initial Conditions [5, 5 ,5, 5]
% r1 = 0.5, r2 = 0,0.1,...,1
% d_males = 0.5, d_females = 2
% time steps = 10^5

%Result, trait 1 dominates regardless unless, trait 2 = 0.5, then we
%converge to the initial condition%

% Simulation 4
                    % m1,f1 ,m2, f2
% Initial Conditions [5, 5 ,5, 5]
% r1 = 0.5, r2 = 0,0.1,...,1
% d_males = 2, d_females = 0.5
% time steps = 10^5

%Result, trait 1 dominates regardless unless, trait 2 = 0.5, then we
%converge to the initial condition%

r1 = 0.5;
d1 = 20;
d2 = 1;
b = 2;
tmax = 10^6;


start = 0;
stop = 1;
step = 0.1;

figure(1)
for i = start:step:stop
    r2 = i;
    y0 = [10^2,10^2,10^2,10^2];
    [t,y] = ode45(@(t,y) f(t,y,r1,r2,d1,d2,b),[0 tmax],y0);
    prop = (y(:,1)+y(:,2))./(y(:,1)+y(:,2)+y(:,3)+y(:,4));
    prop2 =(y(:,3)+y(:,4))./(y(:,1)+y(:,2)+y(:,3)+y(:,4));
    plot(t,prop)
    hold on;
end

k = start:step:stop;
for z=1:length(k)
  leg{z}=sprintf('Trait 2 %.2f',k(z));
end
legend(leg,'location','bestoutside')

function yprime = f(t,y,r1,r2,d1,d2,b)
    dydt = zeros(4,1);
    denom = 1/(y(1)+y(3));
    T = (y(1) + y(2) + y(3) + y(4));
    mu = 0.1;
    
    dydt(1) = -d1*y(1) + (b*denom)*(y(2)*y(1)*r1     + 0.5*y(2)*y(3)*r1     + 0.5*y(1)*y(4)*r2) ;
    dydt(2) = -d2*y(2) + (b*denom)*(y(2)*y(1)*(1-r1) + 0.5*y(2)*y(3)*(1-r1) + 0.5*y(1)*y(4)*(1-r2)) ;
    dydt(3) = -d1*y(3) + (b*denom)*(y(3)*y(4)*r2     + 0.5*y(2)*y(3)*r1     + 0.5*y(1)*y(4)*r2) ;
    dydt(4) = -d2*y(4) + (b*denom)*(y(3)*y(4)*(1-r2) + 0.5*y(2)*y(3)*(1-r1) + 0.5*y(1)*y(4)*(1-r2)) ;
    
    yprime = zeros(4,1);
    
    yprime(1) = dydt(1)/(y(1)+y(3)) - y(1)/((y(1)+y(3))^2)*(dydt(1)+dydt(3));
    yprime(2) = dydt(2)/(y(1)+y(3)) - y(2)/((y(1)+y(3))^2)*(dydt(1)+dydt(3));
    yprime(3) = dydt(3)/(y(1)+y(3)) - y(3)/((y(1)+y(3))^2)*(dydt(1)+dydt(3));
    yprime(4) = dydt(4)/(y(1)+y(3)) - y(4)/((y(1)+y(3))^2)*(dydt(1)+dydt(3));
end