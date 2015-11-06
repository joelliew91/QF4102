% fd_am_put_sor1(r,q,s,x,sig,t,dt,dx,w)
function v = fd_am_put_sor1(r,q,s,x,sig,t,dt,dx,w)
    I = round(12/dx);
    N = round(t/dt);
    VGrid = zeros(I+1,N+1);
    VGrid(I+1,:) = x*exp(-r*(t-(0:dt:t)));
    VGrid(1,:) = 0;
    
    a = dt*(-sig^2/(2*dx^2)-(r-q-sig^2/2)/dx);
    b = dt*(r+1/dt+(sig/dx)^2+(r-q-sig^2/2)/dx);
    c = -dt*sig^2/(2*dx^2);
    
    rep_a = repmat(dt*(-sig^2/(2*dx^2)-(r-q-sig^2/2)/dx),I-1,1);
    rep_b = repmat(dt*(r+1/dt+(sig/dx)^2+(r-q-sig^2/2)/dx),I-1,1);
    rep_c = repmat(-dt*sig^2/(2*dx^2),I-1,1);
    
    A = spdiags([rep_c,rep_b,rep_a],-1:1,I-1,I-1)';
    
    ishift = 1;
    z = (1:I-1)';
    grid_x = (6:-dx:-6)';
    VGrid(:,N+1) = max(x-exp(grid_x),0);

    for j=N:-1:1
        Rhs = VGrid(z+ishift,j+1);
        Rhs(1) = Rhs(1) - a * VGrid(1,j);
        Rhs(I-1) = Rhs(I-1) - c * VGrid(I+1,j);
        prev_u = max(0,x-exp((I:-1:2)'*dx));
        u_gs = zeros(I-1,1);
        flag = 1;
        while(flag)
            u_gs(1) = 1/b * (-c*prev_u(2)+Rhs(1));
            u(1) = max((1-w)*prev_u(1)+w*u_gs(1),max(x*exp(r)-exp(dx*I-q),0)); 
            for i=2:I-2
                u_gs(i) = 1/b * (-a*u(i-1)-c*prev_u(i+1)+Rhs(i));
                u(i) = max((1-w)*prev_u(i)+w*u_gs(i),max(0,x*exp(r)-exp(dx*(I-i+1)-q)));
            end
            
            u_gs(I-1) = 1/b*(-a*u_gs(I-2)-Rhs(I-1));
            u(I-1) = max((1-w)*prev_u(I-1)+w*u_gs(I-1),max(x*exp(r)-exp(dx*2-q),0));
            
            if(abs(u'-prev_u)<0.000001)
                flag = 0;
            else
                prev_u = u';
                u_gs = zeros(I-1,1);
            end
        end
        VGrid(z+ishift,j) = u';
    end
    
    X = log(s)+6;
    ind = round(X/dx);
    v = VGrid(I+1-ind,1);
    
end