t = 500;
times = 1:t;

y1 = zeros(1,t);
y2 = zeros(1,t);
y3 = zeros(1,t);
y4 = zeros(1,t);

y1(1) = 10;
y2(1) = 10;
y3(1) = 10;
y4(1) = 10;

b = 2;
r2 = 1;
r1 = 0.5;
m = 0.1;


for n = 2:t
    denom = (1/(y1(n-1)+y3(n-1)));
    y1(n) = y1(n-1) + (1-m)*(b*denom)*(y2(n-1)*y1(n-1)*r1 +     0.5*y2(n-1)*y3(n-1)*r1 +     0.5*y1(n-1)*y4(n-1)*r2)     + m*(b*denom)*(y3(n-1)*y4(n-1)*r2 +     0.5*y2(n-1)*y3(n-1)*r1 +   0.5*y1(n-1)*y4(n-1)*r2);
    y2(n) = y2(n-1) + (1-m)*(b*denom)*(y2(n-1)*y1(n-1)*(1-r1) + 0.5*y2(n-1)*y3(n-1)*(1-r1) + 0.5*y1(n-1)*y4(n-1)*(1-r2)) + m*(b*denom)*(y3(n-1)*y4(n-1)*(1-r2) + 0.5*y2(n-1)*y3(n-1)*(1-r1) + 0.5*y1(n-1)*y4(n-1)*(1-r2));
    y3(n) = y3(n-1) + (1-m)*(b*denom)*(y3(n-1)*y4(n-1)*r2 +     0.5*y2(n-1)*y3(n-1)*r1 +   0.5*y1(n-1)*y4(n-1)*r2)       + m*(b*denom)*(y2(n-1)*y1(n-1)*r1 +     0.5*y2(n-1)*y3(n-1)*r1 +     0.5*y1(n-1)*y4(n-1)*r2);
    y4(n) = y4(n-1) + (1-m)*(b*denom)*(y3(n-1)*y4(n-1)*(1-r2) + 0.5*y2(n-1)*y3(n-1)*(1-r1) + 0.5*y1(n-1)*y4(n-1)*(1-r2)) + m*(b*denom)*(y2(n-1)*y1(n-1)*(1-r1) + 0.5*y2(n-1)*y3(n-1)*(1-r1) + 0.5*y1(n-1)*y4(n-1)*(1-r2));
end

prop = (y1+y2)./(y1+y2+y3+y4);
prop2 = (y3+y4)./(y1+y2+y3+y4);

plot(times,prop)
hold on;
plot(times,prop2)