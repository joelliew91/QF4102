%  cfl_bm(r,q,s,sigma,t,N)
%N is number of time periods
function cfl_bm(r,q,s,sigma,t,N)
dt=t/N;
u=exp(sigma*sqrt(dt));
d=1/u;
disc = exp(-r*dt);
p = (exp((r-q)*dt)-d)/(u-d);

n=N:-1:0;
    W=ones(1,N+1)- d'.^n; %Value of option at terminal
for i=N-1:-1:0
    k=1:1:i;
    temp=disc*(p*u*W(i(1)+1)+(1-p)*d*W(i(1)+2));
    W=disc*(p*u*W(k)+(1-p)*d*W(k+2));
    W(i(1)+1)=temp;
end
W(1)*s
return