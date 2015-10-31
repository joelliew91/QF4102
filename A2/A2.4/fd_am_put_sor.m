% fd_am_put_sor(r,q,s,x,sig,t,dt,dx,w)
function v = fd_am_put_sor(r,q,s,x,sig,t,dt,dx,w)
    I = round(12/dx);
    N = round(t/dt);
    VGrid = zeros(I+1,N+1);
    VGrid(I+1,:) = x*exp(-r*(t-(0:dt:t)));
    VGrid(1,:) = 0;
    
    a = dt*(-sig^2/(2*dx^2)-(r-q-sig^2/2)/dx); % for the upper value
    b = dt*(r+1/dt+(sig/dx)^2+(r-q-sig^2/2)/dx); % for the centre value
    c = -dt*sig^2/(2*dx^2); % for the lower value
    
    ishift = 1;
    z = (1:I-1)';
    VGrid(:,N+1) = max(x-exp((I:-1:0)'*dx-6),0);
    Pay_off = VGrid(z+ishift,N+1);
    for j=N:-1:1
        Rhs = VGrid(z+ishift,j+1);
        Rhs(1) = Rhs(1) - a * VGrid(1,j);
        Rhs(I-1) = Rhs(I-1) - c * VGrid(I+1,j);
        prev_u = max(0,x-exp((I-1:-1:1)'*dx-6));
        flag = 1;
        while(flag)
            u_gs(1) = 1/b * (-c*prev_u(2)+Rhs(1));
            u(1) = max((1-w)*prev_u(1)+w*u_gs(1),Pay_off(1)); 
            for i=2:I-2
                u_gs(i) = 1/b * (-a*u(i-1)-c*prev_u(i+1)+Rhs(i));
                u(i) = max((1-w)*prev_u(i)+w*u_gs(i),Pay_off(i));
            end
            
            u_gs(I-1) = 1/b*(-a*u(I-2)+Rhs(I-1));
            u(I-1) = max((1-w)*prev_u(I-1)+w*u_gs(I-1),Pay_off(I-1));
            
            if(abs(u'-prev_u)<0.001)
                flag = 0;
            else
                prev_u = u';
            end
        end
        VGrid(z+ishift,j) = u';
    end
    
    X = log(s)+6;
    ind = round(X/dx);
    v = VGrid(I+1-ind,1);

    
end
