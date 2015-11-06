function v=compare(x,x2)
n=length(x);

for i=2:1:n-1
    vv(i)=x(i)-x2(2*i);
end
v=max(abs(vv));
return;
