%fsg_aafsput_quad(r,q,s,sigma,t,N,rho)
function fsg_agput_am_quad(r,q,s,x,running_ave,running_ave_time,sigma,t,N,rho)
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
    final = [1+initial+i*m:-1:i*m+1 i*m:-1:-i*m -i*m-1:-1:-i*m+initial-1];
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
        

        x_up_ll = A{i+1}(m*i-id_A_u_floor+2+max(initial,0))';
        x_up_h = A{i+1}(m*i-id_A_u_floor+max(initial,0))';
        x_up_l = A{i+1}(m*i-id_A_u_floor+1+max(initial,0))';        
        V_up_ll = max(W{j}(m*i-id_A_u_floor+2+max(initial,0)),x-x_up_ll);
        V_up_l = max(W{j}(m*i-id_A_u_floor+1+max(initial,0)),x-x_up_l);
        V_up_h = max(W{j}(m*i-id_A_u_floor+max(initial,0)),x-x_up_h);
        part1 = V_up_ll.*((x_up_l-A_u).*(x_up_h-A_u))./((x_up_l - x_up_ll).*(x_up_h-x_up_ll));
        part2 = V_up_l.*(x_up_ll-A_u).*(x_up_h-A_u)./((x_up_ll-x_up_l).*(x_up_h-x_up_l));
        part3 = V_up_h.*(x_up_ll-A_u).*(x_up_l-A_u)./((x_up_ll-x_up_h).*(x_up_l-x_up_h));
        V_up = part1+part2+part3;


        x_dn_ll = A{i+1}(m*i-id_A_d_floor+2+max(initial,0))';
        x_dn_l = A{i+1}(m*i-id_A_d_floor+1+max(initial,0))';
        x_dn_h = A{i+1}(m*i-id_A_d_floor+max(initial,0))';       
        V_dn_ll = max(W{j+1}(m*i-id_A_d_floor+2+max(initial,0)),x-x_dn_ll);
        V_dn_l = max(W{j+1}(m*i-id_A_d_floor+1+max(initial,0)),x-x_dn_l);
        V_dn_h = max(W{j+1}(m*i-id_A_d_floor+max(initial,0)),x-x_dn_h);
        part1 = V_dn_ll.*(x_dn_l-A_d).*(x_dn_h-A_d)./((x_dn_l-x_dn_ll).*(x_dn_h-x_dn_ll));
        part2 = V_dn_l.*(x_dn_ll-A_d).*(x_dn_h-A_d)./((x_dn_ll-x_dn_l).*(x_dn_h-x_dn_l));
        part3 = V_dn_h.*(x_dn_ll-A_d).*(x_dn_l-A_d)./((x_dn_ll-x_dn_h).*(x_dn_l-x_dn_h));
        V_dn= part1+part2+part3;
        
        F{j}=disc*(p*V_up+(1-p)*V_dn);
        
    end
    W = F;
end

if initial<=0
fl = length(A{1})-1;
x_l=A{1}(fl);
x_ll=A{1}(fl+1);
x_h=A{1}(fl-1);
V_h = max(W{1}(fl-1),x-A{1}(fl-1));
V_ll= max(W{1}(fl+1),x-A{1}(fl+1));
V_l = max(W{1}(fl),x-A{1}(fl));
part1 = V_ll.*(x_l-running_ave).*(x_h-running_ave)./((x_l-x_ll).*(x_h-x_ll));
part2 = V_l.*(x_ll-running_ave).*(x_h-running_ave)./((x_ll-x_l).*(x_h-x_l));
part3 = V_h.*(x_ll-running_ave).*(x_l-running_ave)./((x_ll-x_h).*(x_l-x_h));
V = part1+part2+part3;
elseif initial>0
fl = 2;
x_l=A{1}(fl);
x_ll=A{1}(fl+1);
x_h=A{1}(fl-1);
V_h = max(W{1}(fl-1),x-A{1}(fl-1));
V_ll= max(W{1}(fl+1),x-A{1}(fl+1));
V_l = max(W{1}(fl),x-A{1}(fl));
part1 = V_ll.*(x_l-running_ave).*(x_h-running_ave)./((x_l-x_ll).*(x_h-x_ll));
part2 = V_l.*(x_ll-running_ave).*(x_h-running_ave)./((x_ll-x_l).*(x_h-x_l));
part3 = V_h.*(x_ll-running_ave).*(x_l-running_ave)./((x_ll-x_h).*(x_l-x_h));
V = part1+part2+part3;
end
V
toc
return
