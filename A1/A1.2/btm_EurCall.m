% Sample BTM program for European vanilla call options
% call syntax: OptVal=btm_EurCall(S0, X, r, T, sigma, q, N)
function OptVal=btm_EurCall(S0,X,r,T,sigma,q,N,H)
% set up lattice parameters
dt=T/N; dx=sigma*sqrt(dt);
u=exp(dx); d=1/u;
df=exp(-r*dt); % discount factor
p=(exp((r-q)*dt)-d)/(u-d); % risk?neutral probability
% initialization
j = 1:1:N+1;
V = max(S0*u.^(2*j-N-2) ,0);% to ensure that if the Share price at tiem t is above the barrier we calculate the pay off
for i=1:1:N+1
    if V(i)> H
        V(i)=V(i)-X;
    else
        V(i)=0;
    end
end
% backward recursive through time
for n = N-1:-1:0; % for each time 
j = 1:1:n+1;
V=df*(p*V(j+1)+ (1-p)*V(j));
for jj=1:1:n+1;% for each share price
    if S0*u.^(2*jj-n-2)<H
        V(jj)=0;
    end
end
end
%
OptVal=V(1);
return;