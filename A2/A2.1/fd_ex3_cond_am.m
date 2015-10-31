%v = fd_ex3(r,q,s,x,sig,t)
function v = fd_ex3_cond_am(r,q,s,x,sig,t)
    smax = 3*x;
    h = 0.01;
    I = round(smax/h);
    dt = 1/((sig*I)^2+r*I); % get dt such that it fulfils the monotonicity condition
    N = round(t/dt);
    V_grid = zeros(I+1,N+1);
    % Initiate Boundary values
    V_grid(I+1,:) = h*I*exp(q*(t-(0:dt:t)))-x*exp(-r*(t-(0:dt:t)));
    V_grid(1,:) =0;
    %
    V_grid(:,N+1) = max((0:I)*h - x,0); % Initiate Extreme values
    i = (1:I-1)';
    isq = i.^2;
    coeff_p1=(0.5*sig^2*isq+(r-q)*i)*dt/(1+r*dt);
    coeff_0=(1-sig^2*isq*dt-(r-q)*i*dt)/(1+r*dt);
    coeff_m1=(0.5*sig^2*isq)*dt/(1+r*dt);
    
    ishift=1;
    j = (2:I)';
    for n=N:-1:1 % backward time recursive
        V_grid(j,n+1) = max(V_grid(j,n+1),h*(j-1)-x);
        V_grid(i+ishift,n)=coeff_m1.*V_grid(i-1+ishift,n+1)+coeff_0 .*V_grid(i+ishift,n+1)+coeff_p1.*V_grid(i+1+ishift,n+1);
    end;
    v=V_grid(round(s/h)+ishift,1);

end