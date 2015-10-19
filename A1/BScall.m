%BScall(r,t,s,strike,q,vol)
function x = BScall(r,t,s,strike,q,vol)
moneyness = s/strike;
d1 = (log(moneyness)+(r-q+(vol^2)/2)*t)./(vol*sqrt(t));
d2 = d1- vol*sqrt(t);
x =s.*exp(-q*t).*normcdf(d1)-strike*exp(-r*t)*normcdf(d2);
return