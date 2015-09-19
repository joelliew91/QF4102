% call syntax: OptVal=btm_EurCall(S0,X,r,T,sigma,q,N)
function OptVal=btm_Amcall(S0,X,r,T,sigma,q,N)
% set up lattice parameters
dt=T/N; dx=sigma*sqrt(dt);
u=exp(dx); d=1/u;
df=exp(-r*dt);     % discount factor 
p=(exp((r-q)*dt)-d)/(u-d);  % risk-neutral probability
% initialization
j = 1:1:N+1;  % range of index for price states
V=max(S0*u.^(2*j-N-2)-X,0);
V = ones(length(p),1)*V;
f = size(V);
p = repmat(p',1,f(2));
% backward recursive through time
for n=N:-1:1   
    j = 1:1:n;
    k = 1:1:n+1;
    V1 = max(S0*u.^(2*k-n-2)-X,0);
    f = size(V);
    V1 = ones(f(1),1)*V1;
    V = max(V,V1);
    V=df*(V(:,j+1).*p(:,j+1)+(1-p(:,j)).*V(:,j));
end
OptVal = V;
return

