%fsg_agput_am_linear(r,q,s,x,running_ave,sigma,t,N,rho)
function fsg_agput_am_linear(r,q,s,x,running_ave,sigma,t,N,rho)
dt=t/N;
u=exp(sigma*sqrt(dt));
d=1/u;
disc = exp(-r*dt);
p = (exp((r-q)*dt)-d)/(u-d);
dy = rho*sigma*sqrt(dt);
m=1/rho;

for i=0:N
final = (i+1)*m:-1:-(1+i)*m;
final_s = i:-2:-i;
S{i+1} = s*u.^final_s;
A{i+1} = s*exp(final*dy);
A{i+1}
end

return