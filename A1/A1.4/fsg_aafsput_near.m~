%fsg_aafsput(r,q,s,sigma,t,N,rho)
function fsg_aafsput_near(r,q,s,sigma,t,N,rho)
tic
dt=t/N;
u=exp(sigma*sqrt(dt));
d=1/u;
disc = exp(-r*dt);
p = (exp((r-q)*dt)-d)/(u-d);
dy = rho*sigma*sqrt(dt);
m=1/rho;

for i=0:N
final = i*m:-1:-i*m;
final_s = i:-2:-i;
S{i+1} = s*u.^final_s;
A{i+1} = s*exp(final*dy);
end
for i=1:N+1
    W{i} = max(A{N+1}-S{N+1}(i),0)';
end
for i=N:-1:1
    for j=1:1:i
        A_u = (S{i}(j)*u+i*A{i})/(i+1);
        A_d = (S{i}(j)*d+i*A{i})/(i+1);
        id_A_u = log(A_u/s)/dy;
        id_A_d = log(A_d/s)/dy;
        id_A_u_floor = floor(log(A_u/s)/dy);
        id_A_d_floor = floor(log(A_d/s)/dy);
        
        id_A_d_floor(id_A_d-id_A_d_floor >0.5)=id_A_d_floor(id_A_d-id_A_d_floor >0.5)+1;
        id_A_u_floor(id_A_u-id_A_u_floor >0.5)=id_A_u_floor(id_A_u-id_A_u_floor >0.5)+1;
        id_A_u_floor;
        
        V_up = W{j}(m*i-id_A_u_floor+1);
        V_dn= W{j+1}(m*i-id_A_d_floor+1);
        
        F{j}=disc*(p*V_up+(1-p)*V_dn);

    end
    W = F;
end
W{1}
toc
return
