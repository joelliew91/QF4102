% mc_call(n,runs,x)
function result=mc_call(n,runs,x)
    r = 0.03;t = 0.5;
    s1 = 0.95;sig1 = 0.36;q1 = 0.02;
    s2 = 1.25;sig2 = 0.25;q2 = 0.03;
    s3 = 0.46;sig3 = 0.47;q3 = 0.01;
    rho = [1,0.86,0.22;0.86,1,-0.25;0.22,-0.25,1];
    L = chol(rho)';
    for i=1:runs
        v = randn(3,n);
        e = (L*v)';
        e1 = e(:,1);
        e2 = e(:,2);
        e3 = e(:,3);
        p1 = s1*exp((r-q1-sig1^2/2)*t+sig1*e1*sqrt(t));
        p2 = s2*exp((r-q2-sig2^2/2)*t+sig2*e2*sqrt(t));
        p3 = s3*exp((r-q3-sig3^2/2)*t+sig3*e3*sqrt(t));
        S = max(max(p1,p2),p3);
        val(i) = mean(exp(-r*t)*max(S-x,0));
    end
    result(1) = mean(val);
    result(2) = std(val);
end