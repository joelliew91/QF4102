% bs_call(S,X,r,t,sigma,q)
function y=bs_call(S,X,r,t,sigma,q)
d1=(log(S/X)+(r-q+sigma*sigma/2)*t)/sigma/sqrt(t);
d2=d1-sigma*sqrt(t);
y=X*exp(-r*t)*normcdf(-d2)-S*exp(-q*t)*normcdf(-d1);
end