% Sample Monte-Carlo simulation program for European vanilla call options
%
% call syntax: MC_MA_EurCall(X,no_samples)
%
function v5=MC_MA_EurCall(X,no_samples)

%Setting up all the All the parameters need 

%correlations
p12=0.86;
p13=.22;
p23=-.25;
%starting share price
S0=.95;
S1=1.25;
S2=.46;
%time to maturity and risk free rate
r=0.03;
T=0.5;
%volality
sig=.36;
sig1=.25;
sig2=.47;
%dividends yield.
q=.02;
q1=.03;
q2=.01;


mu=r-q-sig^2/2;

%finding the closed from exact values for the three vanilla call option prices.
d1=(log(S0/X)+(r-q+(sig^2)/2)*T)/(sig*sqrt(T));
d2=d1-sig*sqrt(T);
cc=exp(-q*T)*S0.*normcdf(d1)-exp(-r*T)*X*normcdf(d2);

d1=(log(S1/X)+(r-q1+(sig1^2)/2)*T)/(sig1*sqrt(T));
d2=d1-sig1*sqrt(T);
cc1=exp(-q1*T)*S1.*normcdf(d1)-exp(-r*T)*X*normcdf(d2);

d1=(log(S2/X)+(r-q2+(sig2^2)/2)*T)/(sig2*sqrt(T));
d2=d1-sig2*sqrt(T);
cc2=exp(-q2*T)*S2.*normcdf(d1)-exp(-r*T)*X*normcdf(d2);

for i=1:20 % for 20 runs
%generation of random numbers
epsv=randn(no_samples,1);
epsv1=randn(no_samples,1);
epsv2=randn(no_samples,1);
%to generate correalted set of numbers
eepsv1=p12*epsv+sqrt(1-p12^2)*epsv1;
eepsv2=p13*epsv+(p23-p13*p12)/sqrt(1-p12^2)*epsv1+sqrt((1+2*p23*p12*p13-p12^2-p13^2-p23^2)/(1-p12^2))*epsv2;


mu=r-q-sig^2/2;
mu1=r-q1-sig1^2/2;
mu2=r-q2-sig2^2/2;

% terminal prices in a vector
ST=S0*exp(mu*T+epsv*sig*sqrt(T));

ST1=S1*exp(mu1*T+eepsv1*sig1*sqrt(T));

ST2=S2*exp(mu2*T+eepsv2*sig2*sqrt(T));

STmax=max(ST,max(ST1,ST2));
%optval to store the result for each run.
optval(i)=exp(-r*T)*mean(max(STmax-X,0));

optval2(i)=(exp(-r*T)/3)*(mean(max(ST1-X,0))+mean(max(ST2-X,0))+mean(max(ST-X,0)));
end;

v(1)=mean(optval);
v(2)=std(optval);
A=cov(optval,optval2);
pp=A(1,2);
Beta=pp/(std(optval2)^2);
v1(1)=mean(optval2);
v1(2)=std(optval2);

new1=optval-Beta*(optval2-(cc+cc1+cc2)/3);

v5(1)=mean(new1);
v5(2)=std(new1);




end


