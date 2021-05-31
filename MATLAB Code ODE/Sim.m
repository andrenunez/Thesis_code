function [y,t,val] = Sim

%Number of Traits
n = 50;

%Initial conditions males
m0 = repmat(1/n,1,n);

%Initial onditions females
f0 = repmat(1,1,n);
y0 = [m0,f0];

%Birth rate%
b = 1e-2;

%Initialise all the traits%
r = 1:n;
for l = 1:n
    r(l) = l/n;
end

% returns all pairs where the average of the sum 
% (i+j)/2 produces a valid trait%

v = [1:n,1:n];
C = nchoosek(v,2);
combs = unique(C,'rows');
[t,y] = compute(b,r,n,combs,y0);


%Get total population of females%
tot = ones(length(y(:,1)),1);
parfor m = n+1:2*n
   tot = tot + y(:,m); 
end

%Calculate trait with r = 0.5%
mid = (n)/2;
val = r(mid);

%Plot the proportion of the population that the trait value with r = 0.5
%makes up%
plot(t,(y(:,mid) + y(:,n+mid))./(1+tot))

%Compute ODE%
function [t,y] = compute(b,r,n,combs,y0)
    %options = odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4 1e-5]);
    [t,y] = ode23(@(t,y) fun(t,y,b,r,n,combs),[0 10^4   ],y0);
end

%ODE Function%
function yprime = fun(t,y,b,r,n,combs)
    
    %Initialise 2n dimensional array%
    dydt = zeros(length(y),1);
    yprime = zeros(length(y),1);
    
    %Get the sum of all males%
    denom = 0;
    parfor i = 1:n
       denom = denom + y(i);
    end
    
    %Add all valid pairs of mating combinations to the relevant y(i)
    
    for i = 1:n
        
        %Loop through all the possible mating combinations%
        for z = 1:length(combs)
           temp = combs(z,:);
           j = temp(1);
           k = temp(2);
           
           %Check if the mating combination produces an individual with
           %trait i%
           if (j+k)/2 == i
            dydt(i) = dydt(i) + b*y(j)*y(k+n)*r(k)*(1/denom);
            dydt(i+n) = dydt(i+n) + b*y(j)*y(k+n)*(1-r(k))*(1/denom);
           end
              
           if floor((j+k)/2) == i && floor((j+k)/2) ~= ceil((j+k)/2)
             dydt(i) = dydt(i) + b*0.5*y(j)*y(k+n)*r(k)*(1/denom);
             dydt(i+n) = dydt(i+n) + b*0.5*y(j)*y(k+n)*(1-r(k))*(1/denom);
           end
              
           if ceil((j+k)/2) == i && floor((j+k)/2) ~= ceil((j+k)/2)
             dydt(i) = dydt(i) + b*0.5*y(j)*y(k+n)*r(k)*(1/denom);
             dydt(i+n) = dydt(i+n) + b*0.5*y(j)*y(k+n)*(1-r(k))*(1/denom);
           end
        end
    end
    
    %Get the sum of all the differential equations relating to males%
    % i.e. z'%
    total = 0;
    parfor i = 1:n
       total = total+dydt(i); 
    end
    
    %Perform the transformation to reduce to a 2n-1 dimensional ODE%
    parfor i = 1:2*n
        yprime(i) = (dydt(i)/denom) - (y(i)/denom^2)*total;
    end
    
end

end
