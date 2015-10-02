%fsg_aafsput_quad(r,q,s,sigma,t,N,rho)
function fsg_aafsput_quad(r,q,s,sigma,t,N,rho)
tic
dt=t/N;
u=exp(sigma*sqrt(dt));
d=1/u;
disc = exp(-r*dt);
p = (exp((r-q)*dt)-d)/(u-d);
dy = rho*sigma*sqrt(dt);
m=1/rho;

for i=0:N
    final = i*m:-1:-i*m-1;
    final_s = i:-2:-i;
    S{i+1} = s*u.^final_s;  % Set up the historical Stock Price Grid
    A{i+1} = s*exp(final*dy); % Set up the FSG for A
end

for i=1:N+1
    W{i} = max(A{N+1}-S{N+1}(i),0)';        % Calculate the terminal Payoff value 
                                            % for each leaf at the terminal branch
end

for i=N:-1:1
    for j=1:1:i
        A_u = ((S{i}(j)*u+i*A{i})/(i+1))';  % Value of A at t=i+1 if S goes up
        A_d = ((S{i}(j)*d+i*A{i})/(i+1))';  % Value of A at t=i+1 if S goes dn
        id_A_u_floor = floor(log(A_u/s)/dy);% Locating the index of A_u wrt to the W grid
        id_A_d_floor = floor(log(A_d/s)/dy);% Locating the index of A_d wrt to the W grid
        
        V_up_ll = W{j}(m*i-id_A_u_floor+2); % Corresponding lower-lower bound value of Option for A_u
        V_up_l = W{j}(m*i-id_A_u_floor+1);  % Corresponding lower bound value of Option for A_u
        V_up_h = W{j}(m*i-id_A_u_floor);    % Corresponding upper bound value of Option for A_u 
        x_up_ll = A{i+1}(m*i-id_A_u_floor+2)';% Corresponding lower-lower bound value of x
        x_up_h = A{i+1}(m*i-id_A_u_floor)'; % Corresponding upper bound value of x
        x_up_l = A{i+1}(m*i-id_A_u_floor+1)';% Corresponding lower bound value of x
        part1 = V_up_ll.*((x_up_l-A_u).*(x_up_h-A_u))./((x_up_l - x_up_ll).*(x_up_h-x_up_ll));%Quad Interpolation
        part2 = V_up_l.*(x_up_ll-A_u).*(x_up_h-A_u)./((x_up_ll-x_up_l).*(x_up_h-x_up_l));
        part3 = V_up_h.*(x_up_ll-A_u).*(x_up_l-A_u)./((x_up_ll-x_up_h).*(x_up_l-x_up_h));
        V_up = part1+part2+part3;           % Value of V for S UP branch after interpolation

        V_dn_ll = W{j+1}(m*i-id_A_d_floor+2);% Corresponding lower-lower bound value of Option for A_d
        V_dn_l = W{j+1}(m*i-id_A_d_floor+1); % Corresponding lower bound value of Option for A_d
        V_dn_h = W{j+1}(m*i-id_A_d_floor);   % Corresponding upper bound value of Option for A_d 
        x_dn_ll = A{i+1}(m*i-id_A_d_floor+2)';% Corresponding lower-lower bound value of x
        x_dn_l = A{i+1}(m*i-id_A_d_floor+1)';% Corresponding lower bound value of x
        x_dn_h = A{i+1}(m*i-id_A_d_floor)';% Corresponding upper bound value of x
        part1 = V_dn_ll.*(x_dn_l-A_d).*(x_dn_h-A_d)./((x_dn_l-x_dn_ll).*(x_dn_h-x_dn_ll));%Quad Interpolation
        part2 = V_dn_l.*(x_dn_ll-A_d).*(x_dn_h-A_d)./((x_dn_ll-x_dn_l).*(x_dn_h-x_dn_l));
        part3 = V_dn_h.*(x_dn_ll-A_d).*(x_dn_l-A_d)./((x_dn_ll-x_dn_h).*(x_dn_l-x_dn_h));
        V_dn= part1+part2+part3;            % Value of V for S DN branch after interpolation
        
        F{j}=disc*(p*V_up+(1-p)*V_dn);      % Value of V 
        
    end
    W = F;
end
W{1}(1)
toc
return
