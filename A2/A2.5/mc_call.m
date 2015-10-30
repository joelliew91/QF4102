% mc_call(n,runs,x)
function result=mc_call(n,runs,x)
    r = 0.03;t = 0.5;
    s1 = 0.95;sig1 = 0.36;q1 = 0.02;
    s2 = 1.25;sig2 = 0.25;q2 = 0.03;
    s3 = 0.46;sig3 = 0.47;q3 = 0.01;
    for i=1:runs
        v = randn(n,3);
        e1 = v(:,1);
        e2 = 0.86*e1+sqrt(1-0.86^2)*v(:,2);
        e3 = 0.22*e1+(-0.25-0.22*0.86)/sqrt(1-0.86^2)*e2+sqrt((1+2*(-0.25)*0.86*0.22-0.86^2-0.22^2-0.25^2)/(1-0.86^2))*v(:,3);
        p1 = s1*exp((r-q1-sig1^2/2)*t+sig1*e1*sqrt(t));
        p2 = s2*exp((r-q2-sig2^2/2)*t+sig2*e2*sqrt(t));
        p3 = s3*exp((r-q3-sig3^2/2)*t+sig3*e3*sqrt(t));
        S = max(max(p1,p2),p3);
        val(i) = mean(max(S-x,0));
    end
    result(1) = mean(val);
    result(2) = std(val);
end