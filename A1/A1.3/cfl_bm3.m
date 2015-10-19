%  cfl_bm2(r,q,s,smin,sigma,t,N)
%  N is number of time periods
%  if the call option is issued at t=0, put smin=s
function output=cfl_bm3(r,q,s,smin,sigma,t,N)
x=min(s,smin)/s;
dt=t/N;
u=exp(sigma*sqrt(dt));
d=1/u;
k_u=ceil(log(x)/log(d));
k_d=k_u-1;
disc = exp(-r*dt);
p = (exp((r-q)*dt)-d)/(u-d);

for i=0:1:N
j=k_u+i:-1:max(0,k_d-i); %Set indices for Terminal value,with an uppper and lower limit incase the smin != s
x_grid{i+1} = d.^j;
x_grid{i+1} = x_grid{i+1}'; %Set the value of x for each dt
end
W{N+1}=max(1-x_grid{N+1},0); %Value of Option at t=T

for i=N:-1:1
    if k_d-i+1<=0 %Check for lower bound i.e. the 1's at the bottom of the graph
        l = length(x_grid{i+1}); % To get the last index of the node
        temp = disc*(p*u*W{i+1}(l-1)+(1-p)*d*W{i+1}(l)); %Calculate boundary option value
    end
    j=1:1:length(x_grid{i+1})-2;
    W{i}=disc*(p*u*W{i+1}(j)+(1-p)*d*W{i+1}(j+2)); %interior option values
    if k_d-i+1<=0
        W{i}(length(W{i})+1)=temp; %Input the lower bound value
    end
end
den = x_grid{1}(2)-x_grid{1}(1); %denominator of the linear interpolation factor
fac = (x - x_grid{1}(1))/den; %factor of the linear interpolation method
t=fac*W{1}(2)+(1-fac)*W{1}(1);
output=s*t;
end
