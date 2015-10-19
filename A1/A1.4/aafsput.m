%aafsput(r,q,s,sigma,t,N)
function aafsput(r,q,s,sigma,t,N)
tic
dt=t/N;
u=exp(sigma*sqrt(dt));
d=1/u;
disc = exp(-r*dt);
p = (exp((r-q)*dt)-d)/(u-d);

S = repmat(s,2^N,1); %Create the size of the values of Stock at t=T
A = repmat(s,2^N,1); %Create the size of the values of A at t = T
for i=1:1:N
    UP=repmat(u,2^(N-i),1);
    DN=repmat(d,2^(N-i),1);
    combine = [UP;DN];
    combine = repmat(combine,2^(i-1),1); % This is to create a suitable Up and 
                                         % down factor to do an
                                         % element-wise multiplication to
                                         % get S(t+dt). For eg from S(0) to
                                         % S(dt), half of the final
                                         % elements would have went Up and
                                         % the other half went down. From
                                         % S(dt) to S(2dt), half of the
                                         % half elements (that went up) will
                                         % continue UP, while the other
                                         % half of half go Down.
    S = S.*combine;
    A = A+S;
   
end
A = A./(N+1); % Average the A among the time points
V = max(A-S,0); % Terminal Payoff value

for i=N:-1:1
    up_ind=1:2:2^i;
    dn_ind=2:2:2^i;
    V = disc*(p*V(up_ind)+(1-p)*V(dn_ind)); % Discount the Expected Value back to t=0
end    
aafput=V
toc
return
