% Sample BTM program for European vanilla call options
% call syntax: OptVal=btm_EurCall(S0, X, r, T, sigma, q, N)
function OptVal=btm_amCall(S0,X,r,T,sigma,q,N,H)
% set up lattice parameters
dt=T/N; dx=sigma*sqrt(dt);
u=exp(dx); d=1/u;
df=exp(-r*dt); % discount factor
p=(exp((r-q)*dt)-d)/(u-d); % risk?neutral probability
% initialization
j = 1:1:N+1;
V = S0*u.^(2*j-N-2);
for i=1:1:N+1% to check at the final Share price above the barrier , if it is find the payoff else it is zero.
    if V(i)> H
        V(i)=V(i)-X;
    else
        V(i)=0;
    end
end
% backward recursive through time
for n = N-1:-1:0;
j = 1:1:n+1;
V=df*(p*V(j+1)+ (1-p)*V(j));


for jj=1:1:n+1;% at time n, we want to check whether the share price is above or below the barrier. if it is above, we use a vector VV to store the value of excercising it at that time. We use VVV to store the value that we have gotten before to ensure that is no mismatch in matrix size later on.
    if S0*u.^(2*jj-n-2)>H
        
        
        VV(jj)=S0*u.^(2*jj-n-2)-X;
        VVV(jj)=V(jj);
    else
        VV(jj)=0;% since the share price is below the barrier that the option is out and the value becomes zero
        VVV(jj)=0;

    end
end
V=max(VVV,VV);% we want the max of the either the value of holding the option at the point of time or excersing the option at that time
end

%
OptVal=V(1);
return;