%aafsput(r,q,s,sigma,t,N)
function aafsput(r,q,s,sigma,t,N)
tic
dt=t/N;
u=exp(sigma*sqrt(dt));
d=1/u;
disc = exp(-r*dt);
p = (exp((r-q)*dt)-d)/(u-d);

S = repmat(s,2^N,1);
A = repmat(s,2^N,1);
for i=1:1:N
    UP=repmat(u,2^(N-i),1);
    DN=repmat(d,2^(N-i),1);
    combine = [UP;DN];
    combine = repmat(combine,2^(i-1),1);
    S = S.*combine;
    A = A+S;
   
end
A = A./(N+1);
V = max(A-S,0);

for i=N:-1:1
    up_ind=1:2:2^i;
    dn_ind=2:2:2^i;
    V = disc*(p*V(up_ind)+(1-p)*V(dn_ind));
end    
aafput=V
toc
return
