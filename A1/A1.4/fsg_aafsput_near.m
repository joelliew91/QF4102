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
S{i+1} = s*u.^final_s;      % Set up the historical Stock Price Grid
A{i+1} = s*exp(final*dy);   % Set up the FSG for A
end
for i=1:N+1
    W{i} = max(A{N+1}-S{N+1}(i),0)';        % Calculate the terminal Payoff value 
                                            % for each leaf at the terminal branch
end
for i=N:-1:1
    for j=1:1:i
        A_u = (S{i}(j)*u+i*A{i})/(i+1);     % Value of A at t=i+1 if S goes up
        A_d = (S{i}(j)*d+i*A{i})/(i+1);     % Value of A at t=i+1 if S goes dn
        id_A_u = log(A_u/s)/dy;             % Get the index of the A_u
        id_A_d = log(A_d/s)/dy;             % Get the index of the A_d
        id_A_u_floor = floor(log(A_u/s)/dy);% Locating the index of A_u wrt to the W grid
        id_A_d_floor = floor(log(A_d/s)/dy);% Locating the index of A_d wrt to the W grid
        
        id_A_d_floor(id_A_d-id_A_d_floor >=0.5)=id_A_d_floor(id_A_d-id_A_d_floor >=0.5)+1; %Checking for nearest point
        id_A_u_floor(id_A_u-id_A_u_floor >=0.5)=id_A_u_floor(id_A_u-id_A_u_floor >=0.5)+1; %Checking for nearest point
        
        V_up = W{j}(m*i-id_A_u_floor+1);    % Value of V for S UP branch after interpolation
        V_dn= W{j+1}(m*i-id_A_d_floor+1);   % Value of V for S DN branch after interpolation  
        
        F{j}=disc*(p*V_up+(1-p)*V_dn);      % Value of V

    end
    W = F;
end
W{1}
toc
return
