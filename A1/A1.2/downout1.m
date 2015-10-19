% downout1(S,X,r,T,sigma,q,H)
function c=downout1(S,X,r,T,sigma,q,H)
S0=S;
lambda=(r-q+sigma^2/2)/(sigma^2); %to get value of lambda ,x1, y1
x1=log(S0/H)/sigma/sqrt(T)+lambda*sigma*sqrt(T);
y1=log(H/S0)/sigma/sqrt(T)+lambda*sigma*sqrt(T);
c=S0*normcdf(x1)*exp(-q*T)-X*exp(-r*T)*normcdf(x1-sigma*sqrt(T))-S0*exp(-q*T)*((H/S0)^(2*lambda))*normcdf(y1)+X*exp(-r*T)*((H/S0)^(2*lambda-2))*normcdf(y1-sigma*sqrt(T)); % to true value of a downout call with the parameter provided. 

return;
