%  cfl_bm_am(r,q,s,smin,sigma,t,N)
%N is number of time periods
function cfl_bm_am(r,q,s,smin,sigma,t,N)
x=min(s,smin)/s;
dt=t/N;
u=exp(sigma*sqrt(dt));
d=1/u;
k=abs(ceil(log(x)/log(u)));
disc = exp(-r*dt);
p = (exp((r-q)*dt)-d)/(u-d);

j=(N+k+1):-1:max(0,k-N);%Set indices for Terminal value
setx=d.^j;
W=ones(1,length(j))- setx; %Value of option at terminal

for i=N-1:-1:0
    payoff_length=length(W);
    m=k+i+2:-1:0;
    m=m(1:payoff_length);
    V=ones(1,length(m))-d.^m;
    W = max(V,W);
    if k-i<=0 %Check for terminal lower bound
        temp = disc*(p*u*W(length(W)-1)+(1-p)*d*W(length(W))); %Calculate lower bound if it exists
    end
    j=1:1:length(W)-2;
    W=disc*(p*u*W(j)+(1-p)*d*W(j+2));
    if k-i<=0
        W(length(W)+1)=temp; %Input the lower bound value
    end
end
z=length(setx);
cfl_bm_am=(x-setx(z-k-1))/(setx(z-k)-setx(z-k-1))*(s*W(2))+(setx(z-k)-x)/(setx(z-k)-setx(z-k-1))*(s*W(1))
return
