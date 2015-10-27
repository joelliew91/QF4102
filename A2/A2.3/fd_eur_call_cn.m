% fd_eur_call_cn(r,q,s,x,sig,t,dt,dx)
function v = fd_eur_call_cn(r,q,s,x,sig,t,dt,dx)
    I = round(12/dx);
    N = round(t/dt);
    VGrid = zeros(I+1,N+1);
    VGrid(I+1,:) = 0;
    VGrid(1,:) = exp(6) - x*exp(-r*(t-(0:dt:t)));
    
    i = (1:I-1)';
    
    grid_x = (6:-dx:-6)';
    VGrid(:,N+1) = max(exp(grid_x)-x,0);

    a = repmat(-0.25*(sig^2/dx^2+(r-q-sig^2/2)/dx),I-1,1);
    b = repmat(0.5*(sig^2/dx^2+r),I-1,1);
    c = repmat(0.25*(-sig^2/dx^2+(r-q-sig^2/2)/dx),I-1,1);
    
    A = spdiags([c,b,a],-1:1,I-1,I-1)';
    
    ishift = 1;
    
    for j=N:-1:1
        Rhs = zeros(I-1,1);
        Rhs(1) = Rhs(1) - a(1) * VGrid(1,j);
        Rhs(I-1) = Rhs(I-1) - c(I-1) * VGrid(I+1,j);
        A_left = eye(I-1)/dt + A;
        A_right = (eye(I-1)/dt - A) * VGrid(2:I,j+1) + Rhs;
        VGrid(ishift+i,j) = linsolve(A_left,A_right);
    end
    
    X = log(s)+6;
    ind = round(X/dx);
    v = VGrid(I+1-ind,1);
end