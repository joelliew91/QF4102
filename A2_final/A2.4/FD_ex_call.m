% Finite difference - explicit scheme for call vanilla call options
% calling syntax:
% v=FD_ex_call(S0,X,r,T,sig,dt,ds,q)
function v=FD_ex_call(S0,X,r,T,sig,dt,ds,q)

Smax=3*X;
I=round(Smax/ds);

N=floor(T/dt);
NN=ceil(sig^2*(I^2)*T+r*I*T);
if NN>N
N=NN;
dt=T/N;
N=round(T/dt);

end
Vgrid=zeros(I+1,N+1);% finite difference grid


%boundary conditions
Vgrid(1,:)=0; % S=0
Vgrid(I+1,:)=((I+1)*ds-X)*exp(-r*(T-(0:dt:T)));

%Terminal condtion
Vgrid(:,N+1)=max((0:I)*ds-X,0);
i=(1:I-1)';
isq=i.^2;
%explicit Scheme 3
coeff_p1=(0.5*sig^2*isq+(r-q)*i)*dt/(1+r*dt);
coeff_0=(1-sig^2*isq*dt-(r-q)*i*dt)/(1+r*dt);
coeff_m1=(0.5*sig^2*isq)*dt/(1+r*dt);

%backward time recursive
ishift=1;

for n=N:-1:1
    Vgrid(i+ishift,n)=max(coeff_m1.*Vgrid(i-1+ishift,n+1)+coeff_0.*Vgrid(i+ishift,n+1)+coeff_p1.*Vgrid(i+1+ishift,n+1),(i)*ds-X);
    %american option with early excercise.
    
    %Vgrid(i+ishift,n)=coeff_m1.*Vgrid(i-1+ishift,n+1)+coeff_0.*Vgrid(i+ishift,n+1)+coeff_p1.*Vgrid(i+1+ishift,n+1);
    %vanilla call
end;
v=Vgrid(round(S0/ds)+ishift,1);

return;


