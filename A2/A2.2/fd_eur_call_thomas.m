% fd_eur_call_thomas(r,q,s,x,sig,t,dt,dx)
function v = fd_eur_call_thomas(r,q,s,x,sig,t,dt,dx)
    I = round(12/dx);
    N = round(t/dt);
    VGrid = zeros(I+1,N+1);
    VGrid(I+1,:) = 0;
    VGrid(1,:) = exp(6) - x*exp(-r*(t-(0:dt:t)));
    
    i = (1:I-1)';
    i_sq = i.^2;
    
    grid_x = (6:-dx:-6)';
    VGrid(:,N+1) = max(exp(grid_x)-x,0);
    a = repmat(-0.5*dt*(sig^2/dx^2+(r-q-sig^2/2)/dx),I-1,1);
    b = repmat(dt*(r+1/dt+sig^2/dx^2),I-1,1);
    c = repmat(0.5*dt*((r-q-sig^2/2)/dx-sig^2/dx^2),I-1,1);
    
    ishift = 1;
    

    for j=N:-1:1
        Rhs = VGrid(i+ishift,j+1);
        Rhs(1) = Rhs(1) - a(1)*VGrid(0+ishift,j);
        Rhs(I-1) = Rhs(I-1) - c(I-1)*VGrid(I+ishift,j);
        n = length(c);
        alpha = zeros(n-1,1);
        beta= zeros(n-1,1);
        alpha(1) = c(1)/b(1);
        beta(1) = Rhs(1)/b(1);
        for k=2:n-1
            alpha(k) = c(k)/(b(k)-alpha(k-1)*a(k));
        end
        for k=2:n-1
            beta(k) = (Rhs(k)-beta(k-1)*a(k))/(b(k)-alpha(k-1)*a(k));
        end
        VGrid(n+ishift,j) = beta(n-1);
        for k=n-1:-1:1
            VGrid(k+ishift,j) = beta(k)-alpha(k)*VGrid(k+ishift+1,j);
        end

    end
    X = log(s)+6;
    ind = round(X/dx);
    v = VGrid(I+1-ind,1);
end
