%fsg_aafsput(r,q,s,sigma,t,N,rho)
function fsg_aafsput1(r,q,s,sigma,t,N,rho)
dt=t/N;
u=exp(sigma*sqrt(dt));
d=1/u;
disc = exp(-r*dt);
p = (exp((r-q)*dt)-d)/(u-d);
dy = rho*sigma*sqrt(dt);

ind = N:-2:-N;
S = s*u.^ind;
k=N/rho;


S_hist = repmat(s,2^N,1);
A_floor_hist = repmat(s,2^N,1);
A_ceil_hist = repmat(s,2^N,1);
k_floor = zeros(2^N,1);
k_ceil = zeros(2^N,1);
A_hist = repmat(s,2^N,1);
for i=1:1:N
    UP=repmat(u,2^(N-i),1);
    DN=repmat(d,2^(N-i),1);
    combine = [UP;DN];
    combine = repmat(combine,2^(i-1),1);
    S_hist = [S_hist S_hist(:,i).*combine];
    A_hist = [A_hist (A_hist(:,i)*i+S_hist(:,i+1))/(i+1)];
    AF =  (A_floor_hist(:,i)*i+S_hist(:,i+1))/(i+1);
    k_floor = [k_floor floor(log(AF/s)/dy)];
    A_floor = exp(k_floor(:,i+1)*dy)*s;
    A_floor_hist = [A_floor_hist A_floor];
    AC = (A_ceil_hist(:,i)*i+S_hist(:,i+1))/(i+1);
    k_ceil = [k_ceil ceil(log(AC/s)/dy)];
    A_ceil = exp(k_ceil(:,i+1)*dy)*s;
    A_ceil_hist = [A_ceil_hist A_ceil];
end

V_up=max(A_ceil_hist(:,N+1)-S_hist(:,N+1),0);
V_dn=max(A_floor_hist(:,N+1)-S_hist(:,N+1),0);
A_floor_hist(:,N+1)
fac = A_hist(:,N+1) - A_floor_hist(:,N+1);
fac = fac./(A_ceil_hist(:,N+1)-A_floor_hist(:,N+1));
V =V_up.*fac+V_dn.*(1-fac)

for i=N:-1:1
    up_ind=1:2:2^i;
    dn_ind=2:2:2^i;
    fac_ind = 1:2^(N-i+1):2^N;
    V=disc*(p*V(up_ind)+(1-p)*V(dn_ind))
    fac = A_hist(fac_ind,i+1) - A_floor_hist(fac_ind,i+1);
    fac = fac./(A_ceil_hist(fac_ind,i+1)-A_floor_hist(fac_ind,i+1)) 
end

return
