%  cfl_bm_am(r,q,s,smin,sigma,t,N)
%N is number of time periods
function output=cfl_bm_am(r,q,s,smin,sigma,t,N)
x=min(s,smin)/s;
dt=t/N;
u=exp(sigma*sqrt(dt));
d=1/u;
k_u=ceil(log(x)/log(d));
k_d = k_u-1;
disc = exp(-r*dt);
p = (exp((r-q)*dt)-d)/(u-d);

for i=0:1:N
j=(i+k_u):-1:max(0,k_d-i);%Set indices for each dt
x_grid{i+1} = d.^j;%set the x grid values
W{i+1}=ones(1,length(j))- x_grid{i+1}; %Value of option for exercising it immediately
end
for i=N:-1:1
    if k_d-i+1<=0 %Check for lower bound i.e. the 1's at the bottom of the graph
        l = length(x_grid{i+1});% To get the last index of the node
        temp = disc*(p*u*W{i+1}(l-1)+(1-p)*d*W{i+1}(l)); %Calculate boundary option value
    end
    j=1:1:length(W{i+1})-2;
    temp_W=disc*(p*u*W{i+1}(j)+(1-p)*d*W{i+1}(j+2)); %interior option values
    if k_d-i+1<=0
        temp_W(length(temp_W)+1)=temp; %Input the lower bound value if it exists
    end
    W{i} = max(W{i},temp_W); %Account for the option of exercising early ie compare between 
                             %current option value and payoff for
                             %exercising now
end
den = x_grid{1}(2)-x_grid{1}(1); %denominator of the linear interpolation factor
fac = (x - x_grid{1}(1))/den; %factor of the linear interpolation method
t=fac*W{1}(2)+(1-fac)*W{1}(1);
output=s*t;
end
