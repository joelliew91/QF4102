%fsg_aafsput(r,q,s,sigma,t,N,rho)
function fsg_aafsput(r,q,s,sigma,t,N,rho)
dt=t/N;
u=exp(sigma*sqrt(dt));
d=1/u;
disc = exp(-r*dt);
p = (exp((r-q)*dt)-d)/(u-d);

ind = N:-2:-N;
S = s*u.^ind;
k=N/rho;
ind = k:-1:-k;

%A = s*u.^ind


return
