% BS_call(S,X,r,T,sigma,q)
function c =BS_call(S,X,r,T,sigma,q)
d1=(log(S/X)+(r-q+(sigma^2)/2)*T)/(sigma*sqrt(T));
d2=d1-sigma*sqrt(T);
c=exp(-q*T)*S.*normcdf(d1)-exp(-r*T)*X*normcdf(d2);
return
