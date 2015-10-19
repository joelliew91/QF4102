function OptVal=btm_LBFScallrunmin(Smin,S0,q,r,sigma,T,N)  % btm_LBFScallrunmin(4.75,5,0.01,0.02,0.35,0.5,10000)
dt=T/N; dx=sigma*sqrt(dt);          
u=exp(dx); d=1/u;
df=exp(-r*dt);                  % discount factor
p=(exp((r-q)*dt)-d)/(u-d);      % risk neutral probability
X0=min(S0,Smin)/S0;
k=abs(ceil(log(X0)/log(u)));
% initialization
for i=0:N
    A{i+1}=1:i+2+k;               % introduce cell array of N+1 cells to store values
end
i=max(k-N,0):1:N+1+k;
X=u.^(-i);
W=max(1-X,0);
A{N+1}=W;                       % assign W to N+1 cell 
% Backward iterations through time
for n=N-1:-1:0
    if k-n<1
        A{n+1}(1)=df*(p*u*A{n+2}(2)+(1-p)*d*A{n+2}(1));
    end    
    for j=max(k-n,1):1:k+1+n
        A{n+1}(j+1)=df*(p*u*A{n+2}(j+2)+(1-p)*d*A{n+2}(j));
    end
end
OptVal=(X0-X(k+2))/(X(k+1)-X(k+2))*(S0*A{1}(k+1))+(X(k+1)-X0)/(X(k+1)-X(k+2))*(S0*A{1}(k+2));
