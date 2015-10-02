%cfl(r,q,s,sigma,t)
%if option is initiated at t=0, set smin=s
%else input the min of s
function cfl(r,q,s,smin,sigma,t)
a1=(log(s/smin)+t*(r-q+0.5*sigma^2))/(sigma*sqrt(t)); %a1,a2,a3,y1 follows the notation from John Hull's book
a2=a1-sigma*sqrt(t);                                  %(8th edition) taken from page 582 
a3=(log(s/smin)+(-r+q+0.5*sigma^2)*t)/(sigma*sqrt(t));
y1=-(log(s/smin)*2*(r-q-0.5*sigma^2))/sigma^2;
part1=s*exp(-q*t)*normcdf(a1);
part2=s*exp(-q*t)*normcdf(-a1)*(sigma^2)/(2*(r-q));
part3=smin*exp(-r*t)*(normcdf(a2)-exp(y1)*normcdf(-a3)*(sigma^2)/(2*(r-q)));
cfl=part1-part2-part3
return