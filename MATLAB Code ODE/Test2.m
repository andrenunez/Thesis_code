n = 2;
v = [1:n,1:n];
C = nchoosek(v,2);
answer = unique(C,'rows');