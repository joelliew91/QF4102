function v=compare(x,x2)
nn=length(x);

for i=2:nn-1
    vv(i)=x(i)-x2(2*i);
end
v=max(abs(vv));
return;