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
    S{i+1} = (s*u.^final_s)';               % Set up the historical Stock Price Grid
    A{i+1} = (min(s*u.^final,running_min))';% Set up the FSG for A
    A{i+1} = A{i+1}(sum(A{i+1}==running_min):length(A{i+1}));% Remove all duplicates of running_min 
end

for i=1:N+1
    W{i} = max(x-A{N+1},0);                 % Calculate the terminal Payoff value 
                                            % for each leaf at the terminal branch
end

for i=N:-1:1
    for j=1:1:i
        A_u=min(A{i+1},S{i}(j)*u);          % Value of A at t=i+1 if S goes up
        A_d=min(A{i+1},S{i}(j)*d);          % Value of A at t=i+1 if S goes dn
        ind_A_u = log(A_u/s)/dy;            % Get the index of the A_u
        ind_A_d = log(A_d/s)/dy;            % Get the index of the A_d
        if ind_A_d(1)==initial
            ind_A_d = [1;1-round(ind_A_d(2:length(ind_A_d)))+initial]; % If first index of A_d equal to initial,
                                                                       % can ignore the first element and 
                                                                       % change the remaining elements to point to 
                                                                       % the correct corresponding W grid
        else
            ind_A_d = 1-round(ind_A_d(1:length(ind_A_d)))+initial; % Recalibrate if the first index of A_d 
                                                                   % if not equal to initial so as to point 
                                                                   % to the correct corresponding W grid
                                                                   
        end
        if ind_A_u(1)==initial
            ind_A_u = [1;1-round(ind_A_u(2:length(ind_A_u)))+initial]; % If first index of A_u equal to initial,
                                                                       % can ignore the first element and 
                                                                       % change the remaining elements to point to 
                                                                       % the correct corresponding W grid
        else
            ind_A_u = 1 - round(ind_A_u(1:length(ind_A_u)))+initial;   % Recalibrate if the first index of A_u 
                                                                       % if not equal to initial so as to point 
                                                                       % to the correct corresponding W grid
        end

        if length(A{i+1})==1 % To account for cases if there is only 1 element and that the index do not point 
            ind_A_u = 1;     % outside 
            ind_A_d = 1;
        end
        V_up = max(W{j}(ind_A_u),x-A{i+1}(ind_A_u));    % Value of V for S UP branch
        V_dn = max(W{j+1}(ind_A_d),x-A{i+1}(ind_A_d));  % Value of V for S DN branch
        F{j} = disc*(p*V_up+(1-p)*V_dn);                % Value of V
    end
    W=F;

end
W{1}(length(W{1}))
toc
return 