%Calldi(s,x,t,r,q,vol,h)
function cdi =Calldi(s,x,t,r,q,vol,h)
lam = (r-q+vol.^2)/(2*vol.^2);
d1 = log(s/h)/(vol*sqrt(t))+lam*vol*sqrt(t);
d2 = d1 - vol*sqrt(t);
y1 = log(h./s)/(vol*sqrt(t))+lam*vol*sqrt(t);
y2 = y1 - vol*sqrt(t);
C = s.*normpdf(d1)*exp(-q*t)-x*exp(-r*t).*normpdf(d2);
Cdo = C - s.*normpdf(y1).*(exp(-q*t)*(h./s).^(2*lam))+x.*exp(-r*t)*(h./s).^(2*lam-2).*normpdf(y2);
Cdi = C - Cdo;
cdi = Cdi;
return