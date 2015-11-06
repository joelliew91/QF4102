%Finite Difference - Fully implicit Scheme for vanilla call options with
%crank-Nicolson scheme
%Calling syntex:
%v=FD_im_crank(S0,X,r,T,sig,ds,dt,q)
function v=FD_im_crank(S0,X,r,T,sig,ds,dt,q)

N=T/dt;
I=12/ds;
Vgrid=zeros(I+1,N+1); %finite difference grid

%boundary conditions

Vgrid(1,:)=0; % since exp(-6) is really small
Vgrid(I+1,:)=(exp(6)-X)*exp(-r*(T-(0:dt:T)));

%terminal conditions

Vgrid(:,N+1)=max(exp(-6:ds:6)-X,0);

%coeffcients for Vn
var=sig^2;
alpha=r-q-var/2;
coeff_m1=-dt*(var/(4*ds^2) -alpha/(4*ds));
coeff_0=1+dt*var/(2*ds^2) +dt*r/2;
coeff_p1=-dt*(var/(4*ds^2) +alpha/(4*ds));

%coeffcients for Vn+1
ccoeff_m1=dt*(var/(4*ds^2) -alpha/(4*ds));
ccoeff_0=1-dt*var/(2*ds^2) -dt*r/2;
ccoeff_p1=dt*(var/(4*ds^2) +alpha/(4*ds));

%setting the matrix for AVn
a=ones(I-1,1)*coeff_m1;
b=ones(I-1,1)*coeff_0;
c=ones(I-1,1)*coeff_p1;
A=spdiags([c,b,a],-1:1,I-1,I-1)'+zeros(I-1,I-1);

%setting the matrix for AA(Vn+1)
aa=ones(I-1,1)*ccoeff_m1;
bb=ones(I-1,1)*ccoeff_0;
cc=ones(I-1,1)*ccoeff_p1;
AA=spdiags([cc,bb,aa],-1:1,I-1,I-1)'+zeros(I-1,I-1)';

ii=(1:I-1);
%back ward time recursion
for m=N:-1:1
rr=Vgrid(ii+1,m+1);%set up the right hand side
rr=AA*rr;
rr(1)=rr(1)-coeff_m1*Vgrid(1,m)+ccoeff_m1*Vgrid(1,m+1);

rr(I-1)=rr(I-1)-coeff_p1*Vgrid(I+1,m)+ccoeff_p1*Vgrid(I+1,m+1);
new=linsolve(A,rr);
Vgrid(ii+1,m)=new;
end;

find=floor((log(S0)+6)/ds)+1;
if find - (log(S0)+6)/ds+1>0.5
    find=find+1;
end;
v=Vgrid(find,1);


end
