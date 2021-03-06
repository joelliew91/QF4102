% Sample BTM program for European vanilla call options
% call syntax: OptVal=btm_EurCall(S0,X,r,T,sigma,q,N)
function OptVal=btm_EurCall(S0,X,r,T,sigma,q,N)
% set up lattice parameters
dt=T/N; dx=sigma*sqrt(dt);
u=exp(dx); d=1/u;
df=exp(-r*dt);     % discount factor 
p=(exp((r-q)*dt)-d)/(u-d);  % risk-neutral probability
% initialization
j = 1:1:N+1;  % range of index for price states
V=max(S0*u.^(2*j-N-2)-X,0);
% backward recursive through time
for n=N-1:-1:0   
   j = 1:1:n+1; 
   V=df*(p*V(j+1)+(1-p)*V(j));
end
OptVal=V(1);


