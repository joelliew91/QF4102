%Finite Difference - Fully implicit Scheme for vanilla call options
%Calling syntex:
%v=FD_newSOR_putvector(S0,X,r,T,sig,ds,dt,q,w)
function v=FD_newSOR_putvector(S0,X,r,T,sig,ds,dt,q,w)
N=T/dt;
I=12/ds;
Vgrid=zeros(I+1,N+1); %finite difference grid
error=10^(-10);
II=0;
%boundary conditions

Vgrid(1,:)=(X-exp(-6))*exp(-r*(T-(0:dt:T)));
Vgrid(I+1,:)=0;%since X-exp(6) is most likely negative number

%terminal conditions

Vgrid(:,N+1)=max(X-exp(-6:ds:6),0)';

%coeffcients
var=sig^2;
alpha=r-q-var/2;
dss=ds^2;
coeff_m1=-dt*(var/(2*dss));
coeff_0=1+dt*var/(dss) +dt*r+dt*alpha/ds;
coeff_p1=-dt*(var/(2*dss) +alpha/(ds));

%coeff_m1=-dt*(var/(2*ds^2))+dt*alpha/ds;
%coeff_0=1+dt*var/(ds^2) +dt*r-dt*alpha/ds;
%coeff_p1=-dt*(var/(2*ds^2));

% setting up matrix to solve the inital guess
a=ones(I-1,1)*coeff_m1;
b=ones(I-1,1)*coeff_0;
c=ones(I-1,1)*coeff_p1;
A=spdiags([c,b,a],-1:1,I-1,I-1)'+zeros(I-1,I-1);

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   




ii=(1:I-1);

%back ward time recursion
for m=N:-1:1
rr=Vgrid(ii+1,m+1);%set up the right hand side
rr(1)=rr(1)-coeff_m1*Vgrid(1,m);
rr(I-1)=rr(I-1)-coeff_p1*Vgrid(I+1,m);
new0=zeros(I-1,1);
new1=max(X-exp((-6+ds):ds:(6-ds))',0);
new1=linsolve(A,rr);
while (new1-new0)'*(new1-new0)>error
%while sum(new1-new0)>error
II=II+1;
    new0=new1;
    new1(1)=max((1-w)*new0(1)+(w/coeff_0)*(-coeff_p1*new0(2)+rr(1)),(X-exp(-6+ds)));
   
    for j=2:I-2
  
        new1(j)=max((1-w)*new0(j)+(w/coeff_0)*(-coeff_m1*new1(j-1)-coeff_p1*new0(j+1)+rr(j)),max((X-exp(-6+(j)*ds)),0));
    end;
    new1(I-1)=max((1-w)*new0(I-1)+(w/coeff_0)*(-coeff_m1*new1(I-2)+rr(I-1)),(X-exp(6-ds)));

   
end;


Vgrid(ii+1,m)=new1;
end;
II
find=floor((log(S0)+6)/ds);
find1=(log(S0)+6)/ds;
%if find - (log(S0)+6)/ds+1>0.5
%    find=find+1;
%end;
%v=(1-find1+find)*Vgrid(find+1,1)+(find1-find)*Vgrid(find+2,1);
v=Vgrid(:,1);


end
