%  cfl_bm2(r,q,s,smin,sigma,t,N)
%N is number of time periods
function cfl_bm3(r,q,s,smin,sigma,t,N)
x=min(s,smin)/s;
dt=t/N;
u=exp(sigma*sqrt(dt));
d=1/u;
k_u=ceil(log(x)/log(d));
k_d=k_u-1;
disc = exp(-r*dt);
p = (exp((r-q)*dt)-d)/(u-d);

for i=0:1:N
j=k_u+i:-1:max(0,k_d-i); %Set indices for Terminal value
setx=d.^j;
x_grid{i+1} = setx;
x_grid{i+1} = x_grid{i+1}';
end
W{N+1}=max(1-x_grid{N+1},0);
%Value of option at terminal
for i=N:-1:1
    if k_d-i<0 %Check for terminal lower bound
        l = length(x_grid{i+1});
        temp = disc*(p*u*W{i+1}(l-1)+(1-p)*d*W{i+1}(l)); %Calculate lower bound if it exists
    end
    j=1:1:length(x_grid{i+1})-2;
    W{i}=disc*(p*u*W{i+1}(j)+(1-p)*d*W{i+1}(j+2));
    if k_d-i<0
        W{i}(length(W{i})+1)=temp; %Input the lower bound value
    end
end
W{1};
x_grid{1};
den = x_grid{1}(1)-x_grid{1}(2);
fac = (x - x_grid{1}(2))/den;
t=fac*W{1}(1)+(1-fac)*W{1}(2);
cfl_bm3=s*t
return
