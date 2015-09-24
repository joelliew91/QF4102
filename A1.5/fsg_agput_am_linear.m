%fsg_agput_am_linear(r,q,s,x,running_ave,running_ave_time,sigma,t,N,rho)
function fsg_agput_am_linear(r,q,s,x,running_ave,running_ave_time,sigma,t,N,rho)
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
final = [initial+i*m:-1:i*m+1 i*m:-1:-i*m -i*m-1:-1:-i*m+initial];
final_s = i:-2:-i;
S{i+1} = s*u.^final_s;
A{i+1} = s*exp(final*dy);
end
for i=1:N+1
    W{i} = max(x-A{N+1},0)';
end

for i=N:-1:1
    for j=1:1:i
        A_u = exp((log(S{i}(j)*u)+(running_ave_time+i)*log(A{i}))/(running_ave_time+i+1))';
        A_d = exp((log(S{i}(j)*d)+(running_ave_time+i)*log(A{i}))/(running_ave_time+i+1))';
        id_A_u_floor = floor(log(A_u/s)/dy);
        id_A_d_floor = floor(log(A_d/s)/dy);
        
        x_up_h = A{i+1}(m*i-id_A_u_floor+max(initial,0))';
        x_up_l = A{i+1}(m*i-id_A_u_floor+max(initial,0)+1)';
        V_up_l = max(W{j}(m*i-id_A_u_floor+1+max(initial,0)),-x_up_l+x);
        V_up_h = max(W{j}(m*i-id_A_u_floor+max(initial,0)),-x_up_h+x);     

        fac = x_up_h - x_up_l;
        fac = (A_u-x_up_l)./fac;
        V_up = fac.*V_up_h+(1-fac).*V_up_l;


        x_dn_h = A{i+1}(m*i-id_A_d_floor+max(initial,0))';
        x_dn_l = A{i+1}(m*i-id_A_d_floor+1+max(initial,0))';
        V_dn_l = max(W{j+1}(m*i-id_A_d_floor+1+max(initial,0)),-x_dn_l+x);
        V_dn_h = max(W{j+1}(m*i-id_A_d_floor+max(initial,0)),-x_dn_h+x);   
        
        fac = x_dn_h - x_dn_l;
        fac = (A_d-x_dn_l)./fac;
        V_dn= fac.*V_dn_h+(1-fac).*V_dn_l;
        
        F{j}=disc*(p*V_up+(1-p)*V_dn);
        
    end
    W = F;
end
if initial<=0
fl = length(A{1});
V_h = max(W{1}(fl-1),-A{1}(fl-1)+x);
V_l = max(W{1}(fl),-A{1}(fl)+x);
fac = A{1}(fl-1)-A{1}(fl);
fac = (running_ave-A{1}(fl))/fac;
V = fac*V_h+(1-fac)*V_l;
elseif initial>0
V_h = max(W{1}(1),-A{1}(1)+x);
V_l = max(W{1}(2),-A{1}(2)+x);
fac = A{1}(1)-A{1}(2);
fac = (running_ave-A{1}(2))/fac;
V = fac*V_h+(1-fac)*V_l;
end

V
toc
return