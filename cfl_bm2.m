%  cfl_bm2(r,q,s,smin,sigma,t,N)
%N is number of time periods
function cfl_bm2(r,q,s,smin,sigma,t,N)
x=min(s,smin)/s;
dt=t/N;
u=exp(sigma*sqrt(dt));
d=1/u;
k=floor(x/log(u));
disc = exp(-r*dt);
p = (exp((r-q)*dt)-d)/(u-d);

if k>=1
    k+N
n=max(k+n,0);
zeroes=sum(n<=0);
n=n(1:length(n)-zeroes+1);
W=ones(1,length(n))- d'.^n; %Value of option at terminal
for i=length(W)-1:-1:0
    temp = disc*(p*u*W(i(1)+1)+(1-p)*d*W(i(1)+2));
    j=1:1:i;
    W=disc*(p*u*W(j)+(1-p)*d*W(j+2));
    W(i(1)+1)=temp;
end
W(1)*s
return