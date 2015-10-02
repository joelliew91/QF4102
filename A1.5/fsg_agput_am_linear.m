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
final = [initial+i*m:-1:i*m+1 i*m:-1:-i*m -i*m-1:-1:-i*m+initial];% To ensure that the grid includes the 
                                                                  % floor of the index of running_ave
final_s = i:-2:-i;
S{i+1} = s*u.^final_s; % Set up the historical Stock Price Grid
A{i+1} = s*exp(final*dy);% Set up the FSG for A
end
for i=1:N+1
    W{i} = max(x-A{N+1},0)';                % Calculate the terminal Payoff value 
                                            % for each leaf at the terminal branch
end

for i=N:-1:1
    for j=1:1:i
        A_u = exp((log(S{i}(j)*u)+(running_ave_time+i)*log(A{i}))/(running_ave_time+i+1))'; % Value of A at t=i+1 if S goes up 
        A_d = exp((log(S{i}(j)*d)+(running_ave_time+i)*log(A{i}))/(running_ave_time+i+1))'; % Value of A at t=i+1 if S goes dn
        id_A_u_floor = floor(log(A_u/s)/dy);% Locating the index of A_u wrt to the W grid
        id_A_d_floor = floor(log(A_d/s)/dy);% Locating the index of A_d wrt to the W grid
        
        x_up_h = A{i+1}(m*i-id_A_u_floor+max(initial,0))'; % Corresponding upper bound value of x
        x_up_l = A{i+1}(m*i-id_A_u_floor+max(initial,0)+1)';% Corresponding lower bound value of x
        V_up_l = max(W{j}(m*i-id_A_u_floor+1+max(initial,0)),-x_up_l+x);% Corresponding lower bound value of Option for A_u
                                                                        % with the option to exercise early                                                  
        V_up_h = max(W{j}(m*i-id_A_u_floor+max(initial,0)),-x_up_h+x);  % Corresponding upper bound value of Option for A_u   
                                                                        % with the option to exercise early
        fac = x_up_h - x_up_l; % Linear Interpolation
        fac = (A_u-x_up_l)./fac;
        V_up = fac.*V_up_h+(1-fac).*V_up_l; % Value of V for S UP branch after interpolation


        x_dn_h = A{i+1}(m*i-id_A_d_floor+max(initial,0))'; % Corresponding upper bound value of x
        x_dn_l = A{i+1}(m*i-id_A_d_floor+1+max(initial,0))';% Corresponding lower bound value of x
        V_dn_l = max(W{j+1}(m*i-id_A_d_floor+1+max(initial,0)),-x_dn_l+x);% Corresponding lower bound value of Option for A_d
                                                                          % with the option to exercise early
        V_dn_h = max(W{j+1}(m*i-id_A_d_floor+max(initial,0)),-x_dn_h+x);  % Corresponding upper bound value of Option for A_d
                                                                          % with the option to exercise early 
        
        fac = x_dn_h - x_dn_l; % Linear interpolation
        fac = (A_d-x_dn_l)./fac;
        V_dn= fac.*V_dn_h+(1-fac).*V_dn_l; % Value of V for S DN branch after interpolation
        
        F{j}=disc*(p*V_up+(1-p)*V_dn); % Value of V 
        
    end
    W = F;
end
% Since the option is initiate when t != 0, need to
% interpolate once more
if initial<=0 % Check if running_ave is smaller than the current price 
fl = length(A{1});                      
V_h = max(W{1}(fl-1),-A{1}(fl-1)+x); % If smaller, take the last 2 values in the W{1} and A{1} grid    
V_l = max(W{1}(fl),-A{1}(fl)+x);     % for linear interpolation
fac = A{1}(fl-1)-A{1}(fl);
fac = (running_ave-A{1}(fl))/fac;
V = fac*V_h+(1-fac)*V_l;
elseif initial>0
V_h = max(W{1}(1),-A{1}(1)+x); % If bigger, take the first 2 values in the W{1} and A{1} grid
V_l = max(W{1}(2),-A{1}(2)+x); % for linear interpolation
fac = A{1}(1)-A{1}(2);
fac = (running_ave-A{1}(2))/fac;
V = fac*V_h+(1-fac)*V_l;
end

V
toc
return