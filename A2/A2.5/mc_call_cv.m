% mc_call_cv(n,runs,x)
function result=mc_call_cv(n,runs,x)
    r = 0.03;t = 0.5;
    s1 = 0.95;sig1 = 0.36;q1 = 0.02;
    s2 = 1.25;sig2 = 0.25;q2 = 0.03;
    s3 = 0.46;sig3 = 0.47;q3 = 0.01;
    rho = [1,0.86,0.22;0.86,1,-0.25;0.22,-0.25,1];
    L = chol(rho)';
    d11 = (log(s1/x)+(r-q1+sig1^2/2*t))/(sig1*sqrt(t));
    d12 = d11 - sig1*sqrt(t);
    d21 = (log(s2/x)+(r-q2+sig2^2/2*t))/(sig2*sqrt(t));
    d22 = d21 - sig2*sqrt(t);
    d31 = (log(s3/x)+(r-q3+sig3^2/2*t))/(sig3*sqrt(t));
    d32 = d31 - sig3*sqrt(t);
    bs_call1 = s1*normcdf(d11) - normcdf(d12)*x*exp(-r*t);
    bs_call2 = s2*normcdf(d21) - normcdf(d22)*x*exp(-r*t);
    bs_call3 = s3*normcdf(d31) - normcdf(d32)*x*exp(-r*t);
    e_f = (bs_call1+bs_call2+bs_call3)/3;
    for i=1:runs
        v = randn(3,n);
        e = (L*v)';
        e1 = e(:,1);
        e2 = e(:,2);
        e3 = e(:,3);
        p1 = s1*exp((r-q1-sig1^2/2)*t+sig1*e1*sqrt(t));
        p2 = s2*exp((r-q2-sig2^2/2)*t+sig2*e2*sqrt(t));
        p3 = s3*exp((r-q3-sig3^2/2)*t+sig3*e3*sqrt(t));
        S_A = max(max(p1,p2),p3);
        A_mean = mean(exp(-r*t)*max(S_A-x,0));
        diff_A = S_A - A_mean;
        S_B = (max(p1-x,0)+max(p2-x,0)+ max(p3-x,0))/3;
        B_mean = mean(exp(-r*t)*S_B);
        diff_B = S_B - B_mean;
        Beta = sum(diff_B.*diff_A)/sum(diff_B.^2);
        val(1) = A_mean- Beta*(B_mean - e_f); 
    end
    result(1) = mean(val);
    result(2) = std(val);
end