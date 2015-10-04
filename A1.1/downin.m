% downin(S,X,r,T,sigma,q,H)
function d=downin(S,X,r,T,sigma,q,H)
d1=(log(S/X)+(r-q+sigma^2/2)*T)/sigma/sqrt(T);
d2=d1-sigma*sqrt(T);
c=exp(-q*T)*S.*normcdf(d1)-exp(-r*T)*X*normcdf(d2);% to find the value of a call option
lamda=(r-q-sigma^2)/(sigma^2);
x1=log(S/H)/(sigma*sqrt(T))+lamda*sigma*sqrt(T);
y1=log(H./S)/(sigma*sqrt(T))+lamda*sigma*sqrt(T);
Cdo=S.*normcdf(x1)*exp(-q*T)-X*exp(-r*T)*normcdf(x1-sigma*sqrt(T))-S*exp(-q*T).*((H./S).^(2*lamda)).*normcdf(y1)+X*exp(-r*T)*((H./S).^(2*lamda-2)).*normcdf(y1-sigma*sqrt(T));% the value of a down and out option
Cdi=c-Cdo;% value of a down and in option
d=Cdi;
return;

