%fsg_aafsput(r,q,s,sigma,t,N,rho)
function fsg_aafsput(r,q,s,sigma,t,N,rho)
dt=t/N;
u=exp(sigma*sqrt(dt));
d=1/u;
disc = exp(-r*dt);
p = (exp((r-q)*dt)-d)/(u-d);
dy = rho*sigma*sqrt(dt);


ind = N:-2:-N;
S = s*u.^ind;
k=N/rho;
ind = k:-1:-k;

S_hist = repmat(s,2^N,1);
A_hist = repmat(s,2^N,1);
for i=1:1:N
    UP=repmat(u,2^(N-i),1);
    DN=repmat(d,2^(N-i),1);
    combine = [UP;DN];
    combine = repmat(combine,2^(i-1),1);
    S_hist = [S_hist S_hist(:,i).*combine];
    A_hist = [A_hist (A_hist(:,i)*i+S_hist(:,i+1))/(i+1)];
end
k_floor = floor(log(A_hist/s)/dy);
A_floor = exp(k_floor*dy)*s;
A_ceiling = exp(dy*(k_floor+1))*s;
V_ceiling = max(0,A_ceiling(:,N+1)-s)
V_floor = max(0,A_floor(:,N+1)-s)
for i=N:-1:0
    V = 
end
return
