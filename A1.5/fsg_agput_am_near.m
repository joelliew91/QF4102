%fsg_agput_am_near(r,q,s,x,running_ave,running_ave_time,sigma,t,N,rho)
function fsg_agput_am_near(r,q,s,x,running_ave,running_ave_time,sigma,t,N,rho)
tic
dt=t/N;
u=exp(sigma*sqrt(dt));
d=1/u;
disc = exp(-r*dt);
p = (exp((r-q)*dt)-d)/(u-d);
dy = rho*sigma*sqrt(dt);
m=1/rho;
initial = floor(log(running_ave/s)/dy);

for i=0:N
final = max(initial,0)+i*m:-1:-i*m+min(initial,0);
final_s = i:-2:-i;
S{i+1} = s*u.^final_s;
A{i+1} = s*exp(final*dy);
end
for i=1:N+1
    W{i} = max(x-A{N+1},0)';
end
for i=N:-1:1
    for j=1:1:i
        A_u = exp((log(S{i}(j)*u)+(running_ave_time+i)*log(A{i}))/(running_ave_time+i+1));
        A_d = exp((log(S{i}(j)*d)+(running_ave_time+i)*log(A{i}))/(running_ave_time+i+1));
        id_A_u = log(A_u/s)/dy;
        id_A_d = log(A_d/s)/dy;
        id_A_u_floor = floor(log(A_u/s)/dy);
        id_A_d_floor = floor(log(A_d/s)/dy);
        id_A_d_floor(id_A_d-id_A_d_floor >0.5)=id_A_d_floor(id_A_d-id_A_d_floor >0.5)+1;
        id_A_u_floor(id_A_u-id_A_u_floor >0.5)=id_A_u_floor(id_A_u-id_A_u_floor >0.5)+1;
        
        V_up = max(W{j}(m*i-id_A_u_floor+max(initial,0)),(x-s*exp(id_A_u_floor*dy))');
        V_dn= max(W{j+1}(m*i-id_A_d_floor+max(initial,0)),(x-s*exp(id_A_d_floor*dy))');
        
        F{j}=disc*(p*V_up+(1-p)*V_dn);

    end
    W = F;
end
res = abs(running_ave - A{1});

V = max(x-running_ave,W{1}(find(res==min(res))))
toc
return
