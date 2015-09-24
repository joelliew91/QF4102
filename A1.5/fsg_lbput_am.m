%fsg_lbput_am(r,q,s,x,running_min,sigma,t,N)
function fsg_lbput_am(r,q,s,x,running_min,sigma,t,N)
tic
dt=t/N;
u=exp(sigma*sqrt(dt));
d=1/u;
disc = exp(-r*dt);
p = (exp((r-q)*dt)-d)/(u-d);
dy = sigma*sqrt(dt);
initial = ceil(log(running_min/s)/dy);
for i=0:N
    final = i:-1:-i;
    final_s = i:-2:-i;    
    S{i+1} = (s*u.^final_s)';
    A{i+1} = (min(s*u.^final,running_min))';
    A{i+1} = A{i+1}(sum(A{i+1}==running_min):length(A{i+1}));
    length(A{i+1});
end

for i=1:N+1
    W{i} = max(x-A{N+1},0);
end

for i=N:-1:1
    for j=1:1:i
        A_u=min(A{i+1},S{i}(j)*u);
        A_d=min(A{i+1},S{i}(j)*d);
        ind_A_u = log(A_u/s)/dy;
        ind_A_d = log(A_d/s)/dy;
        if ind_A_d(1)==initial
            ind_A_d = [1;1-round(ind_A_d(2:length(ind_A_d)))+initial];
        else
            ind_A_d = 1-round(ind_A_d(1:length(ind_A_d)))+initial;
        end
        if ind_A_u(1)==initial
            ind_A_u = [1;1-round(ind_A_u(2:length(ind_A_u)))+initial];
        else
            ind_A_u = 1 - round(ind_A_u(1:length(ind_A_u)))+initial;
        end
        ind_A_u;
        length(W{j});
        
        if length(A{i+1})==1
            ind_A_u = 1;
            ind_A_d = 1;
        end
        V_up = max(W{j}(ind_A_u),x-A{i+1}(ind_A_u));
        V_dn = max(W{j+1}(ind_A_d),x-A{i+1}(ind_A_d));
        F{j} = disc*(p*V_up+(1-p)*V_dn);
    end
    W=F;

end
W{1}(length(W{1}))
toc
return 